module mem_globaer
  
  use mem_precision, only: &
       small_pc, &  !INTENT(IN)
       one          !INTENT(IN)

  use mem_aerad, only: &
       nelem,    &  !INTENT(IN)
       ngroup,   &  !INTENT(IN)
       ngas,     &  !INTENT(IN)
       nsolute      !INTENT(IN)

  implicit none
  !  @(#) globaer.h  McKie  Oct-1995
  !  This is the global include file for the Ames Aerosol model.
  !  This file is intended to be included in all major model
  !  source code modules that need access to global variables.
  !  All constants and variables in this file are assumed to be
  !  available for use in all major source code modules.
  
  !  Note:  Some support source code modules that do not need any
  !  model variables might not include this file.
  
  !  This file contains the following:
  
  !   Symbolic constant declarations & definitions.
  !   Common block declarations.
  
  !  Note:  When the brief definitions given here are not
  !  entirely clear, pointers to the source code routines
  !  with more complete definitions are given in brackets: {}
  
  
  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  
  !  Declare and define symbolic constants
  !   (Implicit typing in effect unless explicitly specified)
  !  Define text string name of this model
  
  character (LEN=*),parameter :: prognam= 'CARMA'
  
  
  !  Define version tag string for this version of the model
  
  character (LEN=*),parameter :: progtag = '2.01' 
  
  
  !  Define flag to indicate if debugging mode is active
  
  logical,parameter :: debug = .false.
  
  
  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  
  
  !  Start of user-defined symbolic constants
  
  
  !  Define particle number concentration [ # / cm^3 ]
  !  used to decide whether to bypass microphysical processes:
  !  set it to SMALL_PC to never bypass the calculations.
  
  double precision,parameter :: few_pc = small_pc * 1.e5 
  
  
  !  Define core fraction (for core mass and second moment) used
  !  when particle number concentrations are limited to SMALL_PC
  
  double precision,parameter :: fix_coref = one * 0.01 
  
  
  !  End of user-defined symbolic constants
  
  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  
  
  !  The remaining symbolic constants will need no attention from most
  !  users (unless extending the model capabilities or simulating other than
  !  a terrestrial atmosphere)
  
  
  !  Define # vertical grid box boundaries, including top & bottom
  
  integer ::  nzradp1  != nz_rad + 1 
  
  
  !  Define # components of variables with 3 spatial dimensions
  !   (used for collapsing first few dimensions in some calcs)
  
  integer :: nxy       != nx * ny 
  integer :: nxyz      != nxy * nz 
  integer :: nxyzp1    != nxy * nzp1 
  integer :: nxyzrad   != nxy * nz_rad 
  integer :: nxyzradp1 != nxy * nzradp1 
  
  
  !  Define # components of pc() in first 4 and 5 dimensions
  !   (used for collapsing dimensions of <pc> in some calcs)
  
  integer :: npc4      != nxyz * nbin 
  integer :: npc5      != npc4 * nelem 
  
  
  !  Define integer safety marker value placed at end of common blocks
  
  integer,parameter :: isafety = 12345
  
  
  !  Define character safety marker value placed at end of common blocks
  
  character (LEN=5),parameter :: csafety = '12345'
  
  
  !  Define values of flag used for specification of
  !  horizontal transport algorithm
  
  integer,parameter :: i_ppm = 0 
  integer,parameter :: i_galerkin = 1 
  
  
  !  Define values of flag used for vertical transport
  !  boundary conditions
  
  integer,parameter :: i_fixed_conc = 0 
  integer,parameter :: i_flux_spec = 1 
  
  
  !  Define values of flag used for particle element
  !  composition specification
  
  integer,parameter :: i_bcinorg = 0 
  !      INTEGER,PARAMETER :: I_ORG = 0 
  !      INTEGER,PARAMETER :: I_WATER = 1 
  !      INTEGER,PARAMETER :: I_ICE = 2 
  !      INTEGER,PARAMETER :: I_MIXEDWAT = 3 
  
  
  !  Define values of flag used for particle element
  !  type specification
  
  integer,parameter :: i_involatile = 0 
  integer,parameter :: i_volatile = 1 
  integer,parameter :: i_coremass = 2 
  integer,parameter :: i_volcore = 3 
  integer,parameter :: i_core2mom = 4 
  
  !KML Define the mass percentil of the BC core
  real,parameter :: pcore = 10. 
  
  !  Define values of flag used for nucleation process
  !  specification
  
  integer,parameter :: i_dropact = 0 
  integer,parameter :: i_aerfreeze = 1 
  integer,parameter :: i_dropfreeze = 2 
  integer,parameter :: i_mixedfreeze = 3 
  integer,parameter :: i_mixedmelt = 4 
  integer,parameter :: i_icemelt = 5 
  
  
  !  Define values of flag used specify direction in
  !  horizontal transport calculations
  
  integer,parameter :: idirx = 0 
  integer,parameter :: idiry = 1 
  
  
  !  Define values of symbols used to specify horizontal & vertical grid type.
  !   Grid selection is made by defining each of the variables
  !   <igridv> and <igridh> to one of the grid types known to the model.
  
  !   Possible values for igridv:
  !       I_CART cartesian
  !       I_SIG sigma
  
  !    Possible values for igridh:
  !       I_CART   cartesian
  !       I_LL     longitude_latitude
  !       I_LC     lambert_conformal
  !       I_PS     polar_stereographic
  !       I_ME     mercator
  
  integer,parameter :: i_cart = 1 
  integer,parameter :: i_sig = 2 
  integer,parameter :: i_ll = 3 
  integer,parameter :: i_lc = 4 
  integer,parameter :: i_ps = 5 
  integer,parameter :: i_me = 6 
  
  
  !  Define values of flag used to specify calculation of solar zenith angle
  
  integer,parameter :: i_fixed = 0 
  integer,parameter :: i_diurnal = 1 
  
  
  !  Define triple-point temperature (K)
  
  real,parameter :: t0 = 273.16D+0 
  
  
  !  Define circle constant [ unitless ]
  
  real,parameter :: pi = 3.14159265358979D+0 
  
  
  !  Define degrees to radian & radian to degrees factors [unitless]
  
  real,parameter :: deg2rad = pi / 180. 
  real,parameter :: rad2deg = 180. / pi 
  
  
  !  Define acceleration of gravity near Earth surface [ cm/s^2 ]
  
  real,parameter :: grav = 980.6D+0 
  
  
  !  Define planet equatorial radius [ cm ]
  
  real,parameter :: rearth = 6.37D+8 
  
  
  !  Define avogadro's number [ # particles / mole ]
  
  real,parameter :: avg = 6.02252D+23 
  
  
  !  Define Boltzmann's constant [ erg / deg_K ]
  
  real,parameter :: bk = 1.38054D-16 
  
  
  !  Define Loschmidt's number [ mole / cm^3, @ STP ]
  
  real,parameter :: alos = 2.68719D+19 
  
  
  !  Define molecular weight of dry air [ g / mole ]
  
  real,parameter :: wtmol_air = 28.966D+0 
  
  
  !  Define reference pressure, e.g. for potential temp calcs [ dyne / cm^2 ]
  
  real,parameter :: pref = 1000.d+3 
  
  
  !  Define conversion factor for mb to cgm [ dyne / cm^2 ] units
  
  real,parameter :: rmb2cgs = 1000.d+0 
  
  
  !  Define universal gas constant [ erg / deg_K / mole ]
  
  real,parameter :: rgas = 8.31430D+07 
  
  
  !  Define gas constant for dry air [ erg / deg_K / mole ]
  
  real,parameter :: r_air = rgas / wtmol_air 
  
  
  !  Define number of seconds per the planet's day [ s / d ]
  
  real,parameter :: scday = 86400.d+0 
  
  
  !  Define specific heat at constant pres of dry air [ cm^2 / s^2 / deg_K ]
  
  real,parameter :: cp = 1.004D+7 
  
  
  !  Define mass density of liquid water [ g / cm^3 ]
  
  real,parameter :: rho_w = 1.d+0 
  
  
  !  Define mass density of water ice [ g / cm^3 ]
  
  real,parameter :: rho_i = 0.93D+0 
  
  
  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  
  !  Declare global common blocks for model startup control
  !   (These variables are defined fresh for each new run, both
  !    for cold starts and restarts, and do not get dumped into
  !    the output restart file at the end of a run.  They are
  !    all defined at the beginning of the init routine.  They
  !    control the type of initialization that will be performed,
  !    as well as control things that can change from run to run
  !    within a single simulation.)
  
  !   ibtime    Beginning timestep index for this run
  !   ietime    Ending timestep index for this run
  !   endtime   Total simulation time for this run
  !   nprint    # of timesteps between print reports (used when > 0)
  !   nhist     # of timesteps between history output (used when > 0)
  !   nrest     # of timesteps between restart output (used when > 0)
  !   pprint    time period between print reports  (used when nprint < 0)
  !   phist     time period between history outputs (used when nhist < 0)
  !   prest     time period between restart outputs (used when nrest < 0)
  !   khist     Counter for # of history timepoints output this run
  !   kstep
  !   prtofil   Name of output print file
  !   resifil   Name of input restart file
  !   resofil   Name of output restart file
  !   hisofil   Name of output history file
  !   stepofil  Name of time-step diagnostics file
  !   radofil   Name of radiation submodel print output file
  !   do_print  .t. if print output during timestepping is desired
  !   do_hist   .t. if history output during timestepping is desired
  !   do_rest   .t. if restart output during timestepping is desired
  
  character (LEN=50) :: prtofil
  character (LEN=50) :: resifil
  character (LEN=50) :: resofil
  character (LEN=50) :: hisofil
  character (LEN=50) :: stepofil
  character (LEN=50) :: radofil
  logical :: do_print
  logical :: do_hist
  logical :: do_rest
  real :: btime
  integer :: ietime
  real :: endtime
  real :: pprint
  real :: phist
  real :: prest
  integer :: nprint
  integer :: nhist
  integer :: nrest
  integer :: khist
  integer :: kstep
  
  
  
  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  
  !  Declare global common blocks for model grid
  
  !   dom_llx    Domain limit, lower left x coordinate     {initatm}
  !   dom_lly    Domain limit, lower left y coordiante     {initatm}
  !   dom_urx    Domain limit, upper right x coordinate     {initatm}
  !   dom_ury    Domain limit, upper right y coordinate     {initatm}
  !   rlon0      center longitude for LC, PS, LL projections    {initatm}
  !   rlat0      true latitude for ME projection      {initatm}
  !   rlat1      #1 true latitude for LC projection      {initatm}
  !   rlat2      #2 true latitude for LC projection      {initatm}
  !   hemisph    +1.=southern, -1.=northern hemisphere for PS, LC proj {initatm}
  !   igridv     flag to specify desired vertical grid coord system    {initatm}
  !   igridh     flag to specify desired horizontal grid coord system  {initatm}
  !   zl      Altitude at top of layer       {initatm}
  !   zlold      Altitude at top of layer at start of time step
  !   zc      Altitude at layer mid-point      {initatm}
  !   zcold      Altitude at layer mid-point at start of time step
  !   xc      Horizontal position at center of box     {initatm}
  !   yc      Horizontal position at center of box     {initatm}
  !   xl      Horizontal position at lower edge of box     {initatm}
  !   yl      Horizontal position at lower edge of box     {initatm}
  !   xu      Horizontal position at upper edge of box     {initatm}
  !   yu      Horizontal position at upper edge of box     {initatm}
  !   dx      Horizontal grid spacing       {initatm}
  !   dy      Horizontal grid spacing       {initatm}
  !   dz      Thickness of vertical layers      {initatm}
  !   xmet       Horizontal ds/dx (ds is metric distance)     {initatm}
  !   ymet       Horizontal ds/dy (ds is metric distance)     {initatm}
  !   zmet       Vertical ds/dz (ds is metric distance)     {initatm}
  !   zmetl      Vertical ds/dz at beginning of time-step
  !   rlon,rlat  Longitude, latitude [deg] at each horiz grid point    {initatm}
  !   gridname   Text description of horiz & vert grid coord system    {initatm}
  !   iaer1      Safety marker for common block aer1
  
  character (LEN=80) :: gridname
  character (LEN=5) :: caer1s
  real :: dom_llx
  real :: dom_lly
  real :: dom_urx
  real :: dom_ury
  real :: rlon0
  real :: rlat0
  real :: rlat1
  real :: rlat2
  real :: hemisph
  !!real, allocatable :: zl(:,:,:)
  !!real, allocatable :: zc(:,:,:)
  !!real, allocatable :: zlold(:)
  !!real, allocatable :: zcold(:)
  !!real, allocatable :: xc(:,:,:)
  real, allocatable :: xl(:,:,:)
  !!real, allocatable :: xu(:,:,:)
  !!real, allocatable :: yc(:,:,:)
  !!real, allocatable :: yl(:,:,:)
  !!real, allocatable :: yu(:,:,:)
  real, allocatable :: dx(:,:,:)
  real, allocatable :: dy(:,:,:)
  !!real, allocatable :: dz(:,:,:)
  !!real, allocatable :: xmet(:,:,:)
  !!real, allocatable :: ymet(:,:,:)
  !!real, allocatable :: zmet(:,:,:)
  !!real, allocatable :: zmetl(:,:,:)
  real, allocatable :: rlon(:,:)
  real, allocatable :: rlat(:,:)

  integer :: igridv
  integer :: igridh
  integer :: iaer1    
  
  ! Declare alias names for grid stuff with first 2, 3 dimensions treated linearly
  
  !!real, allocatable :: zl2(:,:)
  !!real, allocatable :: zl3(:)
  !!real, allocatable :: zc2(:,:)
  !!real, allocatable :: zc3(:)
  !!real, allocatable :: dz2(:,:)
  !!real, allocatable :: dz3(:)
  !!real, allocatable :: xc2(:,:)
  !!real, allocatable :: xc3(:)
  !!real, allocatable :: yc2(:,:)
  !!real, allocatable :: yc3(:)
  !!real, allocatable :: xl2(:,:)
  !!real, allocatable :: xl3(:)
  !!real, allocatable :: yl2(:,:)
  !!real, allocatable :: yl3(:)
  !!real, allocatable :: xu2(:,:)
  !!real, allocatable :: xu3(:)
  !!real, allocatable :: yu2(:,:)
  !!real, allocatable :: yu3(:)
  real, allocatable :: dx2(:,:)
  real, allocatable :: dx3(:)
  !!real, allocatable :: dy2(:,:)
  !!real, allocatable :: dy3(:)
  !!real, allocatable :: xmet2(:,:)
  !!real, allocatable :: xmet3(:)
  !!real, allocatable :: ymet2(:,:)
  !!real, allocatable :: ymet3(:)
  !!real, allocatable :: zmet2(:,:)
  !!real, allocatable :: zmet3(:)
  !!real, allocatable :: zmetl2(:,:)
  !!real, allocatable :: zmetl3(:)
  
  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  
  !  Declare global common blocks for model option & control variables
  
  !   time        Simulation time at end of current timestep [s]
  !   dtime       Timestep size [s]
  !   dtmin       Minimum time-step
  !   dtmax       Maximum time-step
  !   dpctol      Maximum change in particle concentrations that will be tolerated
  !   dgstol      Maximum change in gas concentrations that will be tolerated
  !   conmax      Minumum relative concentration to consider in varstep  {prestep}
  !   itime       Timestep index at end of current timestep
  !   igelem      Groups to which elements belong     {setupbins}
  !   itype       Particle type specification array     {setupbins}
  !   icomp       Particle compound specification array    {setupbins}
  !   nelemg      Number of elements in group
  !   ncore       Number of core elements (itype = 2) in group
  !   icorelem    Core elements (itype = 2) in group
  !   ienconc     Particle number conc. element for group    {setupbins}
  !   ishape      Describes particle shape for group
  !   icoag       Coagulation mapping array      {setupcoag}
  !   icoagelem   Coagulation element mapping array     {setupcoag}
  !   ifall       Fall velocity options      {setupvfall}
  !   icoagop     Coagulation kernel options     {setupckern}
  !   icollec     Gravitational collection options        {setupckern}
  !   itbnd       Top boundary condition        {vertical}
  !   ibbnd       Bottom boundary condition        {vertical}
  !   itbnd_pc    Top boundary condition flag for particles      {init}
  !   ibbnd_pc    Bottom boundary condition flag for particles     {init}
  !   itbnd_gc    Top boundary condition flag for gas      {init}
  !   ibbnd_gc    Bottom boundary condition flag for gas      {init}
  !   itbnd_ptc   Top boundary condition flag for potential temp.     {init}
  !   ibbnd_ptc   Bottom boundary condition flag for potential temp.    {init}
  !   ihoradv     Specification of horizontal advection algorithm     {init}
  !   do_coag     If .true. then do coagulation       {init}
  !   do_grow     If .true. then do condensational growth and evap.     {init}
  !   do_thermo   If .true. then do solve thermodynamics equation     {init}
  !   do_vtran    If .true. then do vertical transport      {init}
  !   do_ew       If .true. then do east-west transport      {init}
  !   do_ns       If .true. then do north-south transport      {init}
  !   do_ccoef    If .true. then calculate coefficients for htran     {htranglk}
  !   do_varstep  If .true then use variable time-step      {init}
  !   do_step     If .false. then step was too big, so don't advance    {varstep}
  !   do_error    If .true. then do error trapping for debugging     {init}
  !   do_netcdf   If .true. then output history in netcdf file format   {init}
  !   do_parcel   If .true. then do parcel simulation      {init}
  !   ncdf_file   Netcdf file handle integer for internal netcdf use {outhis_ncdf}
  !   sec_mom     If .true. then core second moment (itype = 3) used   {setupgrow}
  !   igrowgas    Gas that condenses into a particle element     {setupgrow}
  !   inucgas     Gas that nucleates a particle group      {setupnuc}
  !   if_nuc      Nucleation conditional array       {setupaer}
  !   inucproc    Nucleation conditional array       {setupaer}
  !   nnuc2elem   Number of elements that nucleate to element     {setupnuc}
  !   inuc2elem   Nucleation transfers particles into element inuc2elem {setupnuc}
  !   ievp2elem   Total evap. transfers particles into group ievp2elem  {setupnuc}
  !   ievp2bin    Total evap. transfers particles into bin ievp2bin     {setupnuc}
  !   inuc2bin    Nucleation transfers particles into bin inuc2bin      {setupnuc}
  !   isolelem    Index of solute for each particle element      {setupnuc}
  !   ix       Current index for spatial grid, general east-west direction
  !   iy       Current index for spatial grid, general north-south direction
  !   iz       Current index for spatial grid, vertical direction
  !   ixy       Current index for spatial grid, linearized 2-D ix,iy (horizontal)
  !   ixyz        Current index for spatial grid, linearized 3-D ix,iy,iz
  !   ntsubsteps  Number of time substeps for fast microphysics processes
  !   maxsubsteps Maximum number of time substeps allowed
  !   minsubsteps Maximum number of time substeps allowed
  !   iaer2       Safety marker for common block aer2

  !   simtitle   Model simulation title string
  !   elemname   Names of particle elements
  !   groupname  Names of particle elements
  !   gasname    Names of gas species
  !   solname    Names of solutes
  !   caer2s     Safety marker for common block aer2s
  
  logical :: do_coag
  logical :: do_grow
  logical :: do_thermo
  logical :: do_ew
  logical :: do_ns
  logical :: do_vtran
  logical :: do_varstep
  logical :: do_step
  logical :: do_ccoef
  logical :: do_error
  !!logical, allocatable :: if_sec_mom(:)
  !!logical, allocatable :: if_nuc(:,:)
  logical :: do_netcdf
  logical :: do_parcel
  character (LEN=80) :: simtitle

  real,parameter,dimension(nelem)                ::  &
                     rhocore=(/1.85/)
  real,parameter,dimension(nelem)                ::  & 
                     rhoshell=(/1.31/)
  real,parameter,dimension(nelem)                ::  &
                     rhoelem=(/100.*rhocore(1)*rhoshell(1)/ &
		     ((100-pcore)*rhocore(1)+pcore*rhoshell(1))/)
  !  This array specIFies the type of each element:
  !     I_INVOLATILE is CN (involatile) number concentration [#/cm^3]
  !     I_VOLATILE   is water droplet (volatile) number concentration [#/cm^3]
  !     I_COREMASS   is core mass concentration [g/cm^3re]
  !     I_VOLCORE    is core mass concentration [g/cm^3] of a volatile core
  !     I_CORE2MOM   is core second moment (mass^2) [g^2/cm^3]
  integer,parameter,dimension(nelem)             :: & 
                    itype=(/I_INVOLATILE/)
  !  This array specifies the composition of each element:
  integer,parameter,dimension(nelem)             :: &
                    icomp=(/I_BCinORG/)
  !  Name for each element
  character (LEN=50),parameter,dimension(nelem)  ::  &
                    elemname=(/'Black-carbon in organic matter'/)
  !  Name for each group
  character (LEN=50),parameter,dimension(ngroup) :: &
                    groupname=(/'Biomass burning particles'/)
  !  Number of elements in each group (elements in a group can include
  !  particle number concentration, mass concentrations of cores, and second
  !  moments of core mass distributions).
  integer,parameter,dimension(ngroup)            :: &
                    nelemg=(/1/)
  !  Minimum radius for each group (used below to calculate <rmassmin>) [cm]
  real,parameter,dimension(ngroup)               :: &
                    rmin=(/1.e-6/)
  !  Ratio of particle mass between successive bins (one for each group)
  real,parameter,dimension(ngroup)               :: &
                    rmrat=(/1.902/)
  !  The values of <ishape> and <eshape> determine particle geometry
  !  (one for each group):
  !
  !    <ishape> = 1: spherical
  !    <ishape> = 2: hexagonal prisms or plates
  !    <ishape> = 3: circular disks, cylinders, or spheroids
  integer,parameter,dimension(ngroup)            :: & 
                    ishape=(/1/)
  !    <eshape> = particle length/diameter
  real,parameter,dimension(ngroup)               :: & 
                    eshape=(/1./)
  !  IF <is_grp_ice> = .true. THEN the particle group is an ice crystal,
  !  ELSE the particle group is liquid (or DOes not grow).  This array
  !  is used to select the appropriate ventilation factors in setupgkern.f.
  logical,parameter,dimension(ngroup)            :: &
                    is_grp_ice=(/.false./)
  !  IF <is_grp_mixed> = .true. THEN the particle group is a mixed ice/liquid
  !  hydrometeor.  This array is used to select processes (such as core
  !  melting) that only occur in mixed particles.
  logical,parameter,dimension(ngroup)            :: &
                    is_grp_mixed=(/.false./)
		    
  character (LEN=20) :: gasname(ngas)
  character (LEN=20) :: solname(nsolute)
  character (LEN=5) :: caer2s
  real :: time
  real :: dtime
  real :: dtmin
  real :: dtmax
  real :: dpctol
  real :: dgstol
  real :: conmax
  real :: time_nuc
  real :: period_nuc
  integer :: maxsubsteps
  integer :: minsubsteps
  real :: dtime_save
  integer :: itime
  integer :: ix
  integer :: iy
  integer :: iz
  integer :: ixy
  integer :: ixyz
  integer :: ntsubsteps
  integer, allocatable :: igelem(:)
  integer, allocatable :: ncore(:)
  integer, allocatable :: ienconc(:)
  !!integer, allocatable :: icoag(:,:)
  !!integer, allocatable :: icoagelem(:,:)
  integer :: ifall
  integer :: icoagop
  integer :: icollec
  integer :: itbnd
  integer :: ibbnd
  integer :: itbnd_pc
  integer :: ibbnd_pc
  integer :: itbnd_gc
  integer :: ibbnd_gc
  integer :: itbnd_ptc
  integer :: ibbnd_ptc
  integer :: ihoradv
  integer :: ncdf_file
  !!integer, allocatable :: imomelem(:)
  !!integer, allocatable :: inucproc(:,:)
  !!integer, allocatable :: igrowgas(:)
  !!integer, allocatable :: inucgas(:)
  !!integer, allocatable :: nnuc2elem(:)
  !!integer, allocatable :: inuc2elem(:,:)
  !!integer, allocatable :: ievp2elem(:)
  !!integer, allocatable :: inuc2bin(:,:,:)
  !!integer, allocatable :: isolelem(:)
  !!integer, allocatable :: ievp2bin(:,:,:)
  !!integer, allocatable :: icorelem(:,:)
  !!integer, allocatable :: nnucelem(:)
  !!integer, allocatable :: nnucbin(:,:,:)
  !!integer, allocatable :: inucelem(:,:)
  !!integer, allocatable :: inucbin(:,:,:,:)
  integer :: iaer2  
  
  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  
  !  Declare global common blocks for particle grid structure

  !   rmin      Radius of particle in first bin [g]
  !   rmassmin  Mass of particle in first bin [g]
  !   rmrat     Ratio of masses of particles in consecutive bins  {setupaer}
  !   r     Radius bins [cm]
  !   rcore     Core Radius bins [cm]
  !   rmass     Mass bins [g]
  !   rmasscore Mass bin core [g]
  !   vol     Particle volume [cm^3]
  !   dr     Width of bins in radius space [cm]
  !   dm     Width of bins in mass space [g]
  !   dv     Width of bins in volume space [cm^3]
  !   rmassup   Upper bin boundary mass [g]
  !   rmasscoreup Upper bin boundary core mass [g]
  !   rup     Upper bin boundary radius [cm]
  !   rcoreup   Upper bin boundary core radius [cm]
  !   rlow      Lower bin boundary radius [cm]
  !   diffmass  Difference between <rmass> values
  !   rhop      Mass density of particle groups [g/cm^3]
  !   rhoelem   Mass density of particle elements [g/cm^3]
  !   rhocore      Mass density of core particle elements [g/cm^3]
  !   rhoshell Mass density of shell particle elements [g/cm^3]
  !   eshape    Ratio of particle length / diameter
  !   iaer3     Safety marker for common block aer3
  
  real, allocatable :: rmassmin(:)
  real, allocatable :: r(:,:)
  real, allocatable :: rmass(:,:)
  real, allocatable :: rmasscore(:,:)
  real, allocatable :: rcore(:,:)
  real, allocatable :: vol(:,:)
  real, allocatable :: dr(:,:)
  real, allocatable :: dm(:,:)
  real, allocatable :: dv(:,:)
  real, allocatable :: rmassup(:,:)
  real, allocatable :: rup(:,:)
  real, allocatable :: rmasscoreup(:,:)
  real, allocatable :: rcoreup(:,:)
  real, allocatable :: rlow(:,:)
  real, allocatable :: diffmass(:,:,:,:)
  !!real, allocatable :: rhop(:,:,:,:,:)
  !!real, allocatable :: rhopcore(:,:,:,:,:)
  !!real, allocatable :: rhshell(:,:,:,:,:)
  integer :: iaer3
  
  !KLF$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$&$$$$$$
  
  ! N   Number of primary-particles in one agglomerate, considering the volume
  !     equivalent     {setupvf}
  ! rf  Outer radius   {setupvf}
  ! ra  Projected area radius  {setupvf}
  ! rmc Continuum-regime mobility radius {setupvf}
  ! rmt Transition-regime mobility radius {setupvf}
  ! rkna Adjusted sphere Knudsen number {setupvf}
  ! dfrac Fractal dimension   {setupvf}
  integer, allocatable :: n(:,:)
  real, allocatable :: rf(:,:)
  !!real, allocatable :: ra(:,:)
  !!real, allocatable :: rmc(:,:)
  !!real, allocatable :: rmt(:,:)
  !!real, allocatable :: rkna (:,:,:)
  !!real, allocatable :: dfrac(:,:)
  integer :: iaer3n  
  
  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  
  !  Declare global common blocks for primary model state variables
  
  !kml
  !   r0        Initial count median diameter of the lognormal particle size 
  !             distribution [cm]
  !   rsig      Initial geometric standard deviation of the particle diameter
  !   totm      Innitial Total mass particle concentration [g/cm3]
  !   pc       Particle concentration     {initaer}
  !   gc       Gas concentration [g/x_units/y_units/z_units]  {initgas}
  !   ptc       Potential temperature concentration [g K/x_units/y_units/z_units]
  !   iaer4       Safety marker for common block aer4
  
  real, allocatable :: r0(:,:,:)
  real, allocatable :: rsig(:,:,:)
  real, allocatable :: totm(:,:,:)
  !REAL, allocatable :: pc(:,:,:,:,:)
  !REAL, allocatable :: gc(:,:,:,:)
  !!real, allocatable :: ptc(:,:,:)
  integer :: iaer4
  
  !  Declare alias names for pc() and rhop() with first 2, 3, 4, 5 dimensions
  !  treated linearly
  
  !REAL, allocatable :: pc2(:,:,:,:)
  !REAL, allocatable :: pc3(:,:,:)
  !REAL, allocatable :: pc4(:,:)
  !REAL, allocatable :: pc5(:)
  !!real, allocatable :: rhop2(:,:,:,:)
  real, allocatable :: rhop3(:,:,:)
  real, allocatable :: rhopcore3(:,:,:)
    
  !EQUIVALENCE( pc2, pc )
  !EQUIVALENCE( pc3, pc )
  !EQUIVALENCE( pc4, pc )
  !EQUIVALENCE( pc5, pc )

  !EQUIVALENCE( rhop2, rhop )
  !EQUIVALENCE( rhop3, rhop )
  !EQUIVALENCE( rhopcore3, rhopcore )

  !  Declare alias name for gc() and ptc() with first 2, 3, 4, 5 dimensions treated linearly
  
  !REAL, allocatable :: gc2(:,:,:)
  !REAL, allocatable :: gc3(:,:)
  !!real, allocatable :: ptc2(:,:)
  !!real, allocatable :: ptc3(:)
  
  !EQUIVALENCE( gc2, gc )
  !EQUIVALENCE( gc3, gc )

  !EQUIVALENCE( ptc2, ptc )
  !EQUIVALENCE( ptc3, ptc )  

  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  
  !  Declare global common blocks for secondary model variables
  
  !   pcl       Particle concentration at beginning of time-step
  !   pcmax       Maximum concentration for each element   {prestep}
  !   pconmax     Maximum particle concentration for each grid point
  !   gcl       Gas concentration at beginning of time-step
  !   ptcl        Potential temperature concentration at beginning of time-step
  !   d_pc        Change in particle concentration due to transport
  !   d_gc        Change in gas concentration due to transport
  !   d_ptc       Change in potential temperature concentration due to transport
  !   cvert       Temporary storage for vertical transport
  !   cvert_tbnd  Temporary storage for vertical transport
  !   cvert_bbnd  Temporary storage for vertical transport
  !   divcor      Correction term for vertical divergence
  !   chor        Temporary storage for horizontal transport
  !   dhor        Temporary storage of grid spacing for horizontal transport
  !   coaglg      Total particle loss rate due to coagulation for group
  !   coagpe      Particle production due to coagulation
  !   rnuclg      Total particle loss rate due to nucleation for group
  !   rnucpe      Particle production due to nucleation
  !   growlg      Total particle loss rate due to growth for group
  !   growle      Partial particle loss rate due to growth for element
  !   growpe      Particle production due to growth
  !   evaplg      Total particle loss rate due to evaporation for group
  !   evapls      Partial particle loss rate due to evaporation for element
  !   evappe      Particle production due to evaporation
  !   gasprod     Gas production term
  !   rlheat      Latent heating rate [deg_K/s]
  !   vertdifu    Upward vertical flux at grid boundary
  !   vertdifd    Downward vertical flux at grid boundary
  !   ftopgas     Downward gas flux across top boundary of model
  !   fbotgas     Upward gas flux across bottom boundary of model
  !   ftoppart    Downward particle flux across top boundary of model
  !   fbotpart    Upward flux particle across bottom boundary of model
  !   ftop        Downward flux across top boundary of model
  !   fbot        Upward flux across bottom boundary of model
  !   pc_topbnd   Particle concentration assumed just above the top boundary
  !   pc_botbnd   Particle concentration assumed just below the bottom boundary
  !   gc_topbnd   Gas concentration assumed just above the top boundary
  !   gc_botbnd   Gas concentration assumed just below the bottom boundary
  !   ptc_topbnd  Thermodynamic variable value assumed just above the top boundary
  !   ptc_botbnd  Thermodynamic variable value assumed just below the bottom 
  !               boundary
  !   cmf         Core mass fraction in a droplet
  !   totevap     .true. if droplets are totally evaporating to CN
  !   inucmin     Index of smallest particle nucleated; used for tot. evap.
  !   inucstep    Index of smallest particle nucleated during time step
  !   iaer5       Safety marker for common block aer5
  
  !!real, allocatable :: pcl(:,:,:)
  !!real, allocatable :: gcl(:,:)
  !!real, allocatable :: ptcl(:)
  !!real, allocatable :: d_pc(:,:,:)
  !!real, allocatable :: d_gc(:,:)
  !!real, allocatable :: d_ptc(:)
  !!real, allocatable :: pcmax(:)
  !!real, allocatable :: cvert(:)
  !!real, allocatable :: divcor(:)
  !!real, allocatable :: chor(:)
  real              :: cvert_tbnd 
  real              :: cvert_bbnd 
  !!real, allocatable :: dhor(:)
  !!real, allocatable :: pconmax(:,:)
  !!real, allocatable :: coaglg(:,:,:)
  !!real, allocatable :: coagpe(:,:,:)
  !!real, allocatable :: rnuclg(:,:,:)
  !!real, allocatable :: rnucpe(:,:)
  !!real, allocatable :: growlg(:,:)
  !!real, allocatable :: growpe(:,:)
  !!real, allocatable :: evaplg(:,:)
  !!real, allocatable :: evappe(:,:)
  !!real, allocatable :: gasprod(:)
  real              :: rlheat
  !!real, allocatable :: vertdifd(:)
  !!real, allocatable :: vertdifu(:)
  !!real, allocatable :: ftopgas(:,:)
  !!real, allocatable :: fbotgas(:,:)
  !!real, allocatable :: ftoppart(:,:,:)
  !!real, allocatable :: fbotpart(:,:,:)
  real              :: ftop
  real              :: fbot
  !!real, allocatable :: pc_topbnd(:,:,:)
  !!real, allocatable :: pc_botbnd(:,:,:)
  !!real, allocatable :: gc_topbnd(:,:)
  !!real, allocatable :: gc_botbnd(:,:)
  !!real, allocatable :: ptc_topbnd(:)
  !!real, allocatable :: ptc_botbnd(:)
  !!real, allocatable :: cmf(:,:)
  !!logical, allocatable :: totevap(:,:)
  !!integer, allocatable :: inucmin(:)
  !!integer, allocatable :: inucstep(:)
  integer              :: iaer5
  
  !  Declare alias name for <evappe> to collapse into one dimension
  
  !!real, allocatable :: evappe5(:,:)
  
  !EQUIVALENCE( evappe5, evappe )
  
  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  
  !  Declare global common blocks for coagulation kernels and bin
  !  pair mapping
  
  !   ck0  Constant coagulation kernel       {setupaer}
  !   grav_e_coll0  Constant value for collection effic.  {setupaer}
  !   ckernel Coagulation kernels [/cm^3/s]       {setupckern}
  !   pkernel Coagulation production variables      {setupcoag}
  !   volx   Coagulation subdivision variable      {setupcoag}
  !   ilow   Bin pairs for coagulation production  {setupcoag}
  !   jlow   Bin pairs for coagulation production  {setupcoag}
  !   iup  Bin pairs for coagulation production  {setupcoag}
  !   jup  Bin pairs for coagulation production  {setupcoag}
  !   npairl Bin pair indices        {setupcoag}
  !   npairu Bin pair indices        {setupcoag}
  !   iaer6  Safety marker for common block aer6
  
  real :: ck0
  real :: grav_e_coll0

  !!real, allocatable :: cbr(:,:,:,:,:)
  !!real, allocatable :: ccd(:,:,:,:,:)
  !!real, allocatable :: cgr(:,:,:,:,:)
  !!real, allocatable :: tim(:,:,:,:,:)
  !!real, allocatable :: tsc(:,:,:,:,:)
  !!real, allocatable :: ckernel(:,:,:,:,:)
  !!real, allocatable :: pkernel(:,:,:,:,:,:,:)
  !!real, allocatable :: volx(:,:,:,:,:)
  !!integer, allocatable :: ilow(:,:,:)
  !!integer, allocatable :: jlow(:,:,:)
  !!integer, allocatable :: iup(:,:,:)
  !!integer, allocatable :: jup(:,:,:)
  !!integer, allocatable :: npairl(:,:)
  !!integer, allocatable :: npairu(:,:)

  integer :: iaer6
  
  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$&$$$$$$$$$$$$$
  
  !  Declare global common blocks for coagulation group pair mapping
  
  !   iglow      Group pairs for coagulation production  {setupcoag}
  !   iglow      Group pairs for coagulation production  {setupcoag}
  !   jglow      Group pairs for coagulation production  {setupcoag}
  !   igup       Group pairs for coagulation production  {setupcoag}
  !   jgup       Group pairs for coagulation production  {setupcoag}
  !   iaer7      Safety marker for common block aer7
  
  !!integer, allocatable :: iglow(:,:,:)
  !!integer, allocatable :: jglow(:,:,:)
  !!integer, allocatable :: igup(:,:,:)
  !!integer, allocatable :: jgup(:,:,:)

  integer :: iaer7
  
  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  
  !  Declare global common blocks for particle fall velocities, transport
  !  rates, and coagulation kernels
  
  !   rkn     knudsen number       {setupvfall}
  !   bpm     Corrects for non-sphericity and non-continuum effects {setupvfall}
  !   bpma      Corrects for non-sphericity and non-continuum effects for an
  !      adjusted sphere {setupvfall}
  !   vf     Fall velocities at layer mid-pt     {setupvfall}
  !   vtrans    Net vertical transport rate at layer boundary   {vertical}
  !   re     Reynolds' number based on <vfall>     {setupvfall}
  !   vf_const  Constant vertical fall velocity when ifall=0   {setupaer}
  !   vertadvu  Upward mass transport rate due to advection    {vertran}
  !   vertadvd  Downward mass transport rate due to advection   {vertran}
  !   htrans    Net horizontal transport rate at layer boundary   {horizont}
  !   hdiff     Horizontal diffusion coefficient at layer boundary    {horizont}
  !   ca     Coefficient for Galerkin horizontal transport   {glkcoef}}
  !   cb     Coefficient for Galerkin horizontal transport   {glkcoef}}
  !   cc     Coefficient for Galerkin horizontal transport   {glkcoef}}
  !   cd     Coefficient for Galerkin horizontal transport   {glkcoef}}
  !   ce     Coefficient for Galerkin horizontal transport   {glkcoef}}
  !   cg     Coefficient for Galerkin horizontal transport   {glkcoef}}
  !   rmfp      mean free path of air molecules [cm]    {setupvf)
  !   iaer8     Safety marker for common block aer8
  
  !!real, allocatable :: rkn(:,:,:)
  !!real, allocatable :: bpm(:,:,:)
  !!real, allocatable :: bpma(:,:,:)
  !!real, allocatable :: vf(:,:,:)
  !!real, allocatable :: vtrans(:)
  !!real, allocatable :: re(:,:,:)

  real :: vf_const

  !!real, allocatable :: vertadvu(:)
  !!real, allocatable :: vertadvd(:)
  !!real, allocatable :: htrans(:)
  !!real, allocatable :: hdiff(:)
  !!real, allocatable :: ca(:,:,:)
  !!real, allocatable :: cb(:,:,:)
  !!real, allocatable :: cd(:,:,:)
  !!real, allocatable :: ce(:,:,:)
  !!real, allocatable :: cf(:,:,:)
  !!real, allocatable :: cg(:,:,:)
  !!real, allocatable :: rmfp(:)

  integer :: iaer8
  
  !  Declare alias names for <rhostar> with first 2, 3 dimensions treated linearly
  
  !REAL, allocatable :: rhostar2(:,:)
  !REAL , allocatable:: rhostar3(:)
  
  !EQUIVALENCE( rhostar2, rhostar )
  !EQUIVALENCE( rhostar3, rhostar )
  
  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  
  !  Declare global common blocks for atmospheric structure.
  !  Note: air density, winds, and diffusion coefficients are in scaled units
  
  !   zbot      height of the bottom of the model [cm]      {initatm}
  !   p_surf    surface pressure [dyne/cm^2]       {initatm}
  !   p_top     Atmospheric pressure at top of model domain [dyne/cm^2] {initatm}
  !   p     Atmospheric pressure at layer mid-pt [dyne/cm^2]     {initatm}
  !   rhoa      Air density at layer mid-pt [g/x_units/y_units/z_units] {initatm}
  !   t     Air temperature at layer mid-pt [deg_K]      {initatm}
  !   t_surf    Air temperature at surface [deg_K]        {initatm}
  !   rmu     Air viscosity at layer mid-pt [g/cm/s]      {initatm}
  !   thcond    Thermal conductivity of dry air [erg/cm/sec/deg_K]      {initatm}
  !   w     Vertical wind speed at layer boundary [z_units/s]     {initatm}
  !   u     East-west wind speed at layer center [x_units/s]     {initatm}
  !   v     North-south wind speed at layer center [y_units/s]      {initatm}
  !   dkz     Vert diffusion coef at layer boundary [z_units^2/s]     {initatm}
  !   dkx     Horiz x diffusion coeff at layer boundary [x_units^2/s] {initatm}
  !   dky     Horiz y diffusion coeff at layer boundary [y_units^2/s] {initatm}
  !   told      Temperature at beginning of time-step
  !   pold      Pressure at beginning of time-step
  !   rhoaold   Air density at beginning of time-step
  !   iaer9     Safety marker for common block aer9
  
  real :: zbot

  !REAL, allocatable :: p(:,:,:)
  !REAL, allocatable :: rhoa(:,:,:)
  !REAL, allocatable :: p_surf(:,:)
  !REAL, allocatable :: p_top(:,:)
  !REAL, allocatable :: t(:,:,:)
  !!real, allocatable :: told(:)
  !!real, allocatable :: pold(:)
  !!real, allocatable :: rhoaold(:)
  !!real, allocatable :: rmu(:)
  !!real, allocatable :: thcond(:)
  real, allocatable :: w(:,:,:)
  !REAL, allocatable :: zmetold(:)
  real, allocatable :: u(:,:,:)
  real, allocatable :: v(:,:,:)
  real, allocatable :: t_surf(:,:)
  !!real, allocatable :: dkz(:,:,:)
  !!real, allocatable :: dkx(:,:,:)
  !!real, allocatable :: dky(:,:,:)

  integer :: iaer9  
  
  !  Declare alias names for atm stuff with first 2, 3 dimensions treated linearly
  
  !!real, allocatable :: p2(:,:)
  !!real, allocatable :: p3(:)
  !!real, allocatable :: t2(:,:)
  !!real, allocatable :: t3(:)
  !!real, allocatable :: u2(:,:)
  !!real, allocatable :: u3(:)
  !!real, allocatable :: v2(:,:)
  !!real, allocatable :: v3(:)
  !!real, allocatable :: w2(:,:)
  !!real, allocatable :: w3(:)
  !!real, allocatable :: dkx2(:,:)
  !!real, allocatable :: dkx3(:)
  !!real, allocatable :: dky2(:,:)
  !!real, allocatable :: dky3(:)
  !!real, allocatable :: dkz2(:,:)
  !!real, allocatable :: dkz3(:)
  !!real, allocatable :: rhoa2(:,:)
  !!real, allocatable :: rhoa3(:)
  !REAL, allocatable :: p_surf2(:)
  !!real, allocatable :: p_top2(:)
  
  !EQUIVALENCE( zl2, zl )
  !EQUIVALENCE( zl3, zl )
  !EQUIVALENCE( zc2, zc )
  !EQUIVALENCE( zc3, zc )
  !EQUIVALENCE( dz2, dz )
  !EQUIVALENCE( dz3, dz )
  !EQUIVALENCE( xc2, xc )
  !EQUIVALENCE( xc3, xc )
  !EQUIVALENCE( xl2, xl )
  !EQUIVALENCE( xl3, xl )
  !EQUIVALENCE( xu2, xu )
  !EQUIVALENCE( xu3, xu )
  !EQUIVALENCE( yc2, yc )
  !EQUIVALENCE( yc3, yc )
  !EQUIVALENCE( yl2, yl )
  !EQUIVALENCE( yl3, yl )
  !EQUIVALENCE( yu2, yu )
  !EQUIVALENCE( yu3, yu )
  !EQUIVALENCE( dx2, dx )
  !EQUIVALENCE( dx3, dx )
  !EQUIVALENCE( dy2, dy )
  !EQUIVALENCE( dy3, dy )
  
  !EQUIVALENCE( p2, p )
  !EQUIVALENCE( p3, p )
  !EQUIVALENCE( rhoa2, rhoa )
  !EQUIVALENCE( rhoa3, rhoa )
  !EQUIVALENCE( t2, t )
  !EQUIVALENCE( t3, t )

!  EQUIVALENCE( u2, u )
!  EQUIVALENCE( u3, u )
!  EQUIVALENCE( v2, v )
!  EQUIVALENCE( v3, v )
!  EQUIVALENCE( w2, w )
!  EQUIVALENCE( w3, w )
!  EQUIVALENCE( dkx2, dkx )
!  EQUIVALENCE( dkx3, dkx )
!  EQUIVALENCE( dky2, dky )
!  EQUIVALENCE( dky3, dky )
!  EQUIVALENCE( dkz2, dkz )
!  EQUIVALENCE( dkz3, dkz )
!  EQUIVALENCE( xmet2, xmet )
!  EQUIVALENCE( xmet3, xmet )
!  EQUIVALENCE( ymet2, ymet )
!  EQUIVALENCE( ymet3, ymet )
!  EQUIVALENCE( zmet2, zmet )
!  EQUIVALENCE( zmet3, zmet )
!  EQUIVALENCE( zmetl2, zmetl )
!  EQUIVALENCE( zmetl3, zmetl )

  !EQUIVALENCE( zmetold2, zmetold )
  !EQUIVALENCE( zmetold3, zmetold )
  !EQUIVALENCE( p_surf2, p_surf )
  !EQUIVALENCE( p_top2, p_top )
  
  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  
  !  Declare global common blocks for condensational growth parameters
  
  !   gwtmol    Molecular weight for gases [g/mol]    {setupgrow}
  !   diffus    Diffusivity of gas in air [cm^2/s]    {setupgrow}
  !   rlhe      Latent heat of evaporation for gas [cm^2/s^2] {setupgrow}
  !   rlhe      Latent heat of ice melting for gas [cm^2/s^2] {setupgrow}
  !   pvapl     Saturation vapor pressure over water [dyne/cm^2] {vaporp}
  !   pvapi     Saturation vapor pressure over ice [dyne/cm^2] {vaporp}
  !   surfctwa  Surface tension of water-air interface  {setupgkern}
  !   surfctiw  Surface tension of water-ice interface  {setupgkern}
  !   surfctia  Surface tension of ice-air interface  {setupgkern}
  !   akelvin   Exponential arg. in curvature term for growth {setupgkern}
  !   akelvini  Curvature term for ice    {setupgkern}
  !   ft     Ventilation factor      {setupgkern}
  !   gro     Growth kernel [UNITS?]    {setupgkern}
  !   gro1      Growth kernel conduction term [UNITS?]  {setupgkern}
  !   gro2      Growth kernel radiation term [UNITS?]  {setupgkern}
  !   gvrat     For converting particle growth to gas loss [UNITS?] {setupgkern}
  !   supsatl   Supersaturation of vapor w.r.t. liquid water [dimless]
  !   supsati   Supersaturation of vapor w.r.t. ice [dimless]
  !   supsatlold Supersaturation (liquid) before time-step    {prestep}
  !   supsatiold Supersaturation (ice) before time-step    {prestep}
  !   scrit     Critical supersaturation for nucleation [dimless] {setupnuc}
  !   sol_ions  Number of ions solute dissociates into  {setupnuc}
  !   solwtmol  Molecular weight of solute     {setupnuc}
  !   rhosol    Mass density of solute    {setupnuc}
  !   rlh_nuc   Nucleation latent heat    {setupaer}
  !   iaer10   Safety marker for common block aer10
  
  !!real, allocatable :: gwtmol(:)
  !!real, allocatable :: diffus(:,:)
  !!real, allocatable :: rlhe(:,:)
  !!real, allocatable :: rlhm(:,:)
  !!real, allocatable :: pvapl(:,:,:,:)
  !!real, allocatable :: pvapi(:,:,:,:)
  !!real, allocatable :: surfctwa(:)
  !!real, allocatable :: surfctiw(:)
  !!real, allocatable :: surfctia(:)
  !!real, allocatable :: akelvin(:,:)
  !!real, allocatable :: akelvini(:,:)
  !!real, allocatable :: ft(:,:,:)
  !!real, allocatable :: gro(:,:,:)
  !!real, allocatable :: gro1(:,:,:)
  !!real, allocatable :: gro2(:,:)
  !!real, allocatable :: gvrat(:,:,:)
  !!real, allocatable :: supsatl(:,:,:,:)
  !!real, allocatable :: supsati(:,:,:,:)
  !!real, allocatable :: supsatlold(:,:)
  !!real, allocatable :: supsatiold(:,:)
  !!real, allocatable :: scrit(:,:,:)
  !!real, allocatable :: sol_ions(:)
  !!real, allocatable :: solwtmol(:)
  !!real, allocatable :: rhosol(:)
  !!real, allocatable :: rlh_nuc(:,:)

  integer :: iaer10
  
  !  Declare alias names for pvapl(), pvapi(), supsatl(), and supsati() with
  !  first 2, 3 dimensions treated linearly
  
  !!real, allocatable :: pvapl2(:,:,:)
  !!real, allocatable :: pvapl3(:,:)
  !!real, allocatable :: pvapi2(:,:,:)
  !!real, allocatable :: pvapi3(:,:)
  !!real, allocatable :: supsatl2(:,:,:)
  !!real, allocatable :: supsatl3(:,:)
  !!real, allocatable :: supsati2(:,:,:)
  !!real, allocatable :: supsati3(:,:)
  
  !EQUIVALENCE( pvapl2, pvapl )
  !EQUIVALENCE( pvapl3, pvapl )
  !EQUIVALENCE( pvapi2, pvapi )
  !EQUIVALENCE( pvapi3, pvapi )
  !EQUIVALENCE( supsatl2, supsatl )
  !EQUIVALENCE( supsatl3, supsatl )
  !EQUIVALENCE( supsati2, supsati )
  !EQUIVALENCE( supsati3, supsati )
  
  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  
  !  Declare global common blocks for nucleation parameters
  
  !   adelf    Coefficient for phase change activation energy  {freznuc}
  !   bdelf    Coefficient for phase change activation energy  {freznuc}
  !   prenuc   Pre-exponential factore for frezzing nucleation {freznuc}
  !   rmiv     Contact angle for ice/nucleus interface    {freznuc}
  !   iaer11   Safety marker for common block rad1
  
  real :: adelf
  real :: bdelf
  real :: prenuc
  real :: rmiv
  integer :: iaer11
  
  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  
  !  Declare global common blocks for radiative transfer parameters
  
  !   do_rad      If .true. then do radiative transfer      {init}
  !   do_solar    If .true. then do solar calculations      {init}
  !   do_ir       If .true. then do infrared calculations      {init}
  !   nrad        # of timesteps between radiation calcs (used when > 0)
  !   prad        Time period between radiation calcs (used when nrad < 0)
  !   isolar_zen  =I_FIXED: fixed, =I_DIURNAL: calculated every time step
  !   u0       cos( solar_zenith_angle )
  !   u0_fixed    Fixed value of cos( solar_zenith_angle )
  !   rad_start   solar time corresponding to <time> = 0 [s]
  !   zsin,zcos   sin and cos terms for solar zenith angle precalculation
  !   wave        Bin-center wavelengths [micron]
  !   radheat     Radiative heating rate [deg_K/s]
  !   qrad        Particle heating rate [deg_K/s]
  !   alb_tomi    Spectrally-integrated albedo at top-of-model
  !   alb_toai    Spectrally-integrated albedo at top-of-atmosphere
  !   alb_toa     Spectrally-resolved albedo at top-of-atmosphere
  !   opd       Spectrally-resolved optical depth
  !   fsl_up      Solar upwelling flux [W m^-2]
  !   fsl_dn      Solar downwelling flux [W m^-2]
  !   fir_up      Infrared upwelling flux [W m^-2]
  !   fir_dn      Infrared downwelling flux [W m^-2]
  !   irad1       Safety marker for common block rad1
  
  real              :: u0
  real              :: u0_fixed
  real              :: rad_start

  real, allocatable :: z_sin(:,:)
  real, allocatable :: z_cos(:,:)
  real, allocatable :: wave(:)
  real, allocatable :: qrad(:,:,:,:,:)
  !!real, allocatable :: radheat(:,:,:)
  real, allocatable :: alb_tomi(:,:)
  real, allocatable :: alb_toai(:,:)
  real, allocatable :: alb_toa(:,:,:)
  real, allocatable :: opd(:,:,:)
  !!real, allocatable :: fsl_up(:,:,:)
  !!real, allocatable :: fsl_dn(:,:,:)
  !!real, allocatable :: fir_up(:,:,:)
  !!real, allocatable :: fir_dn(:,:,:)

  logical :: do_rad
  logical :: do_solar
  logical :: do_ir
  integer :: nrad
  real :: prad
  integer :: isolar_zen
  integer :: irad1
    
  !  Declare alias names for radiative transfer parameters
  !  with first 2, 3 dimensions treated linearly
  
  !!real, allocatable :: radheat2(:,:)
  !!real, allocatable :: radheat3(:)
  !!real, allocatable :: qrad2(:,:,:,:)
  !!real, allocatable :: qrad3(:,:,:)
  !!real, allocatable :: alb_tomi2(:)
  !!real, allocatable :: alb_toai2(:)
  !!real, allocatable :: alb_toa2(:,:)
  !!real, allocatable :: opd2(:,:)
  !!real, allocatable :: fsl_up2(:,:)
  !!real, allocatable :: fsl_up3(:)
  !!real, allocatable :: fsl_dn2(:,:)
  !!real, allocatable :: fsl_dn3(:)
  !!real, allocatable :: fir_up2(:,:)
  !!real, allocatable :: fir_up3(:)
  !!real, allocatable :: fir_dn2(:,:)
  !!real, allocatable :: fir_dn3(:)
  
  !EQUIVALENCE( radheat2, radheat )
  !EQUIVALENCE( radheat3, radheat )
  !EQUIVALENCE( qrad2, qrad )
  !EQUIVALENCE( qrad3, qrad )
  !EQUIVALENCE( alb_tomi2, alb_tomi )
  !EQUIVALENCE( alb_toai2, alb_toai )
  !EQUIVALENCE( alb_toa2, alb_toa )
  !EQUIVALENCE( opd2, opd )
  !EQUIVALENCE( fsl_up2, fsl_up )
  !EQUIVALENCE( fsl_up3, fsl_up )
  !EQUIVALENCE( fsl_dn2, fsl_dn )
  !EQUIVALENCE( fsl_dn3, fsl_dn )
  !EQUIVALENCE( fir_up2, fir_up )
  !EQUIVALENCE( fir_up3, fir_up )
  !EQUIVALENCE( fir_dn2, fir_dn )
  !EQUIVALENCE( fir_dn3, fir_dn )

  contains

    subroutine initial_definitions_globaer()

      use mem_aerad, only: &
           nz_rad,   &  !INTENT(IN)
           nx,       &  !INTENT(IN)
           ny,       &  !INTENT(IN)
           nz,       &  !INTENT(IN)
           nzp1,     &  !INTENT(IN)
           nbin,     &  !INTENT(IN)
           nelem,    &  !INTENT(IN)
           nwave,    &  !INTENT(IN)
           nsol,     &  !INTENT(IN)
           ngroup,   &  !INTENT(IN)
           ngas,     &  !INTENT(IN)
           nxorny,   &  !INTENT(IN)
           nxornyp1, &  !INTENT(IN)
           nsolute      !INTENT(IN)

      implicit none

      ! Variables definitions
      !  Define # vertical grid box boundaries, including top & bottom
      nzradp1   = nz_rad + 1
      !  Define # components of variables with 3 spatial dimensions
      !   (used for collapsing first few dimensions in some calcs)
      nxy       = nx * ny
      nxyz      = nxy * nz
      nxyzp1    = nxy * nzp1
      nxyzrad   = nxy * nz_rad
      nxyzradp1 = nxy * nzradp1
      !  Define # components of pc() in first 4 and 5 dimensions
      !   (used for collapsing dimensions of <pc> in some calcs)
      npc4      = nxyz * nbin
      npc5      = npc4 * nelem

      ! Arrays allocation
      !!allocate(zl(nx,ny,nzp1))
      !!allocate(zc(nx,ny,nz))
      !!allocate(zlold(nxyzp1))
      !!allocate(zcold(nxyz))
      !!allocate(xc(nx,ny,nz))
      allocate(xl(nx,ny,nz))
      !!allocate(xu(nx,ny,nz))
      !!allocate(yc(nx,ny,nz))
      !!allocate(yl(nx,ny,nz))
      !!allocate(yu(nx,ny,nz))
      allocate(dx(nx,ny,nz))
      allocate(dy(nx,ny,nz))
      !!allocate(dz(nx,ny,nz))
      !!allocate(xmet(nx,ny,nz))
      !!allocate(ymet(nx,ny,nz))
      !!allocate(zmet(nx,ny,nz))
      !!allocate(zmetl(nx,ny,nz))
      allocate(rlon(nx,ny))
      allocate(rlat(nx,ny))
      !!allocate(zl2(nxy,nzp1))
      !!allocate(zl3(nxyzp1))
      !!allocate(zc2(nxy,nz))
      !!allocate(zc3(nxyz))
      !!allocate(dz2(nxy,nz))
      !!allocate(dz3(nxyz))
      !!allocate(xc2(nxy,nz))
      !!allocate(xc3(nxyz))
      !!allocate(yc2(nxy,nz))
      !!allocate(yc3(nxyz))
      !!allocate(xl2(nxy,nz))
      !!allocate(xl3(nxyz))
      !!allocate(yl2(nxy,nz))
      !!allocate(yl3(nxyz))
      !!allocate(xu2(nxy,nz))
      !!allocate(xu3(nxyz))
      !!allocate(yu2(nxy,nz))
      !!allocate(yu3(nxyz))
      allocate(dx2(nxy,nz))
      allocate(dx3(nxyz))
      !!allocate(dy2(nxy,nz))
      !!allocate(dy3(nxyz))
      !!allocate(xmet2(nxy,nz))
      !!allocate(xmet3(nxyz))
      !!allocate(ymet2(nxy,nz))
      !!allocate(ymet3(nxyz))
      !!allocate(zmet2(nxy,nz))
      !!allocate(zmet3(nxyz))
      !!allocate(zmetl2(nxy,nz))
      !!allocate(zmetl3(nxyz))
      !!allocate(if_sec_mom(ngroup))
      !!allocate(if_nuc(nelem,nelem))
      allocate(igelem(nelem))
      allocate(ncore(ngroup))
      allocate(ienconc(ngroup))
      !!allocate(icoag(ngroup,ngroup))
      !!allocate(icoagelem(nelem,ngroup))
      !!allocate(imomelem(ngroup))
      !!allocate(inucproc(nelem,nelem))
      !!allocate(igrowgas(nelem))
      !!allocate(inucgas(ngroup))
      !!allocate(nnuc2elem(nelem))
      !!allocate(inuc2elem(5,nelem))
      !!allocate(ievp2elem(nelem))
      !!allocate(inuc2bin(nbin,ngroup,ngroup))
      !!allocate(isolelem(nelem))
      !!allocate(ievp2bin(nbin,ngroup,ngroup))
      !!allocate(icorelem(nelem,nelem))
      !!allocate(nnucelem(nelem))
      !!allocate(nnucbin(ngroup,nbin,ngroup))
      !!allocate(inucelem(nelem,nelem*ngroup))
      !!allocate(inucbin(nbin*ngroup,ngroup,nbin,ngroup))
      allocate(rmassmin(ngroup))
      allocate(r(nbin,ngroup))
      allocate(rmass(nbin,ngroup))
      allocate(rmasscore(nbin,ngroup))
      allocate(rcore(nbin,ngroup))
      allocate(vol(nbin,ngroup))
      allocate(dr(nbin,ngroup))
      allocate(dm(nbin,ngroup))
      allocate(dv(nbin,ngroup))
      allocate(rmassup(nbin,ngroup))
      allocate(rup(nbin,ngroup))
      allocate(rmasscoreup(nbin,ngroup))
      allocate(rcoreup(nbin,ngroup))
      allocate(rlow(nbin,ngroup))
      allocate(diffmass(nbin,ngroup,nbin,ngroup))
      !!allocate(rhop(nx,ny,nz,nbin,ngroup))
      !!allocate(rhopcore(nx,ny,nz,nbin,ngroup))
      !!allocate(rhshell(nx,ny,nz,nbin,ngroup))
      allocate(n(nbin,ngroup))
      allocate(rf(nbin,ngroup))
      !!allocate(ra(nbin,ngroup))
      !!allocate(rmc(nbin,ngroup))
      !!allocate(rmt(nbin,ngroup))
      !!allocate(rkna (nz,nbin,ngroup))
      !!allocate(dfrac(ngroup,nelem))
      allocate(r0(nx,ny,nz))
      allocate(rsig(nx,ny,nz))
      allocate(totm(nx,ny,nz))
      !ALLOCATE(pc(nx,ny,nz,nbin,nelem))
      !ALLOCATE(gc(nx,ny,nz,ngas))
      !!allocate(ptc(nx,ny,nz))
      !ALLOCATE(pc2(nxy,nz,nbin,nelem))
      !ALLOCATE(pc3(nxyz,nbin,nelem))
      !ALLOCATE(pc4(npc4,nelem))
      !ALLOCATE(pc5(npc5))
      !!allocate(rhop2(nxy,nz,nbin,ngroup))
      allocate(rhop3(nxyz,nbin,ngroup))
      allocate(rhopcore3(nxyz,nbin,ngroup))
      !ALLOCATABLE(gc2(nxy,nz,ngas))
      !ALLOCATABLE(gc3(nxyz,ngas))
      !!allocate(ptc2(nxy,nz))
      !!allocate(ptc3(nxyz))
      !!allocate(pcl(nxyz,nbin,nelem))
      !!allocate(gcl(nxyz,ngas))
      !!allocate(ptcl(nxyz))
      !!allocate(d_pc(nxyz,nbin,nelem))
      !!allocate(d_gc(nxyz,ngas))
      !!allocate(d_ptc(nxyz))
      !!allocate(pcmax(nelem))
      !!allocate(cvert(nz))
      !!allocate(divcor(nz))
      !!allocate(chor(nxorny))
      !!allocate(dhor(nxorny))
      !!allocate(pconmax(nxyz,ngroup))
      !!allocate(coaglg(nxyz,nbin,ngroup))
      !!allocate(coagpe(nxyz,nbin,nelem))
      !!allocate(rnuclg(nbin,ngroup,ngroup))
      !!allocate(rnucpe(nbin,nelem))
      !!allocate(growlg(nbin,ngroup))
      !!allocate(growpe(nbin,nelem))
      !!allocate(evaplg(nbin,ngroup))
      !!allocate(evappe(nbin,nelem))
      !!allocate(gasprod(ngas))
      !!allocate(vertdifd(nzp1))
      !!allocate(vertdifu(nzp1))
      !!allocate(ftopgas(nxy,ngas))
      !!allocate(fbotgas(nxy,ngas))
      !!allocate(ftoppart(nxy,nbin,nelem))
      !!allocate(fbotpart(nxy,nbin,nelem))
      !!allocate(pc_topbnd(nxy,nbin,nelem))
      !!allocate(pc_botbnd(nxy,nbin,nelem))
      !!allocate(gc_topbnd(nxy,ngas))
      !!allocate(gc_botbnd(nxy,ngas))
      !!allocate(ptc_topbnd(nxy))
      !!allocate(ptc_botbnd(nxy))
      !!allocate(cmf(nbin,ngroup))
      !!allocate(totevap(nbin,ngroup))
      !!allocate(inucmin(ngroup))
      !!allocate(inucstep(ngroup))
      !!allocate(evappe5(npc5,ngas))
      !!allocate(cbr(nz,nbin,nbin,ngroup,ngroup))
      !!allocate(ccd(nz,nbin,nbin,ngroup,ngroup))
      !!allocate(cgr(nz,nbin,nbin,ngroup,ngroup))
      !!allocate(tim(nz,nbin,nbin,ngroup,ngroup))
      !!allocate(tsc(nz,nbin,nbin,ngroup,ngroup))
      !!allocate(ckernel(nz,nbin,nbin,ngroup,ngroup))
      !!allocate(pkernel(nz,nbin,nbin,ngroup,ngroup,ngroup,6))
      !!allocate(volx(ngroup,ngroup,ngroup,nbin,nbin))
      !!allocate(ilow(ngroup,nbin,nbin*nbin))
      !!allocate(jlow(ngroup,nbin,nbin*nbin))
      !!allocate(iup(ngroup,nbin,nbin*nbin))
      !!allocate(jup(ngroup,nbin,nbin*nbin))
      !!allocate(npairl(ngroup,nbin))
      !!allocate(npairu(ngroup,nbin))
      !!allocate(iglow(ngroup,nbin,nbin*nbin))
      !!allocate(jglow(ngroup,nbin,nbin*nbin))
      !!allocate(igup(ngroup,nbin,nbin*nbin))
      !!allocate(jgup(ngroup,nbin,nbin*nbin))
      !!allocate(rkn(nz,nbin,ngroup))
      !!allocate(bpm(nz,nbin,ngroup))
      !!allocate(bpma(nz,nbin,ngroup))
      !!allocate(vf(nzp1,nbin,ngroup))
      !!allocate(vtrans(nzp1))
      !!allocate(re(nz,nbin,ngroup))
      !!allocate(vertadvu(nzp1))
      !!allocate(vertadvd(nzp1))
      !!allocate(htrans(nxornyp1))
      !!allocate(hdiff(nxorny))
      !!allocate(ca(2,nx,ny))
      !!allocate(cb(2,nx,ny))
      !!allocate(cd(2,nx,ny))
      !!allocate(ce(2,nx,ny))
      !!allocate(cf(2,nx,ny))
      !!allocate(cg(2,nx,ny))
      !!allocate(rmfp(nz))
      !ALLOCATE(rhostar2(nxy,nz))
      !ALLOCATE(rhostar3(nxyz))
      !ALLOCATE(p(nx,ny,nz))
      !ALLOCATE(rhoa(nx,ny,nz))
      !ALLOCATE(p_surf(nx,ny))
      !ALLOCATE(p_top(nx,ny))
      !ALLOCATE(t(nx,ny,nz))
      !!allocate(told(nxyz))
      !!allocate(pold(nxyz))
      !!allocate(rhoaold(nxyz))
      !!allocate(rmu(nz))
      !!allocate(thcond(nz))
      allocate(w(nx,ny,nzp1))
      !ALLOCATE(zmetold(nxyz))
      allocate(u(nx,ny,nz))
      allocate(v(nx,ny,nz))
      allocate(t_surf(nx,ny))
      !!allocate(dkz(nx,ny,nzp1))
      !!allocate(dkx(nx,ny,nzp1))
      !!allocate(dky(nx,ny,nzp1))
      !!allocate(p2(nxy,nz))
      !!allocate(p3(nxyz))
      !!allocate(t2(nxy,nz))
      !!allocate(t3(nxyz))
      !!allocate(u2(nxy,nzp1))
      !!allocate(u3(nxyzp1))
      !!allocate(v2(nxy,nzp1))
      !!allocate(v3(nxyzp1))
      !!allocate(w2(nxy,nzp1))
      !!allocate(w3(nxyzp1))
      !!allocate(dkx2(nxy,nzp1))
      !!allocate(dkx3(nxyzp1))
      !!allocate(dky2(nxy,nzp1))
      !!allocate(dky3(nxyzp1))
      !!allocate(dkz2(nxy,nzp1))
      !!allocate(dkz3(nxyzp1))
      !!allocate(rhoa2(nxy,nz))
      !!allocate(rhoa3(nxyz))
      !ALLOCATE(p_surf2(nxy))
      !!allocate(p_top2(nxy))
      !!allocate(gwtmol(ngas))
      !!allocate(diffus(nz,ngas))
      !!allocate(rlhe(nz,ngas))
      !!allocate(rlhm(nz,ngas))
      !!allocate(pvapl(nx,ny,nz,ngas))
      !!allocate(pvapi(nx,ny,nz,ngas))
      !!allocate(surfctwa(nz))
      !!allocate(surfctiw(nz))
      !!allocate(surfctia(nz))
      !!allocate(akelvin(nz,ngas))
      !!allocate(akelvini(nz,ngas))
      !!allocate(ft(nz,nbin,ngroup))
      !!allocate(gro(nz,nbin,ngroup))
      !!allocate(gro1(nz,nbin,ngroup))
      !!allocate(gro2(nz,ngroup))
      !!allocate(gvrat(nbin,nelem,ngas))
      !!allocate(supsatl(nx,ny,nz,ngas))
      !!allocate(supsati(nx,ny,nz,ngas))
      !!allocate(supsatlold(nxyz,ngas))
      !!allocate(supsatiold(nxyz,ngas))
      !!allocate(scrit(nz,nbin,ngroup))
      !!allocate(sol_ions(nsolute))
      !!allocate(solwtmol(nsolute))
      !!allocate(rhosol(nsolute))
      !!allocate(rlh_nuc(nelem,nelem))
      !!allocate(pvapl2(nxy,nz,ngas))
      !!allocate(pvapl3(nxyz,ngas))
      !!allocate(pvapi2(nxy,nz,ngas))
      !!allocate(pvapi3(nxyz,ngas))
      !!allocate(supsatl2(nxy,nz,ngas))
      !!allocate(supsatl3(nxyz,ngas))
      !!allocate(supsati2(nxy,nz,ngas))
      !!allocate(supsati3(nxyz,ngas))
      !!allocate(radheat2(nxy,nz))
      !!allocate(radheat3(nxyz))
      !!allocate(qrad2(nxy,nz,nbin,ngroup))
      !!allocate(qrad3(nxyz,nbin,ngroup))
      !!allocate(alb_tomi2(nxy))
      !!allocate(alb_toai2(nxy))
      !!allocate(alb_toa2(nxy,nsol))
      !!allocate(opd2(nxy,nwave))
      !!allocate(fsl_up2(nxy,nzradp1))
      !!allocate(fsl_up3(nxyzradp1))
      !!allocate(fsl_dn2(nxy,nzradp1))
      !!allocate(fsl_dn3(nxyzradp1))
      !!allocate(fir_up2(nxy,nzradp1))
      !!allocate(fir_up3(nxyzradp1))
      !!allocate(fir_dn2(nxy,nzradp1))
      !!allocate(fir_dn3(nxyzradp1))
      allocate(z_sin(nx,ny))
      allocate(z_cos(nx,ny))
      allocate(wave(nwave))
      allocate(qrad(nx,ny,nz,nbin,ngroup))
      !!allocate(radheat(nx,ny,nz))
      allocate(alb_tomi(nx,ny))
      allocate(alb_toai(nx,ny))
      allocate(alb_toa(nx,ny,nsol))
      allocate(opd(nx,ny,nwave))
      !!allocate(fsl_up(nx,ny,nzradp1))
      !!allocate(fsl_dn(nx,ny,nzradp1))
      !!allocate(fir_up(nx,ny,nzradp1))
      !!allocate(fir_dn(nx,ny,nzradp1))

    end subroutine initial_definitions_globaer

    ! *************************************************************************

    subroutine final_definitions_globaer()
      implicit none

      ! Arrays allocation
      !!deallocate(zl)
      !!deallocate(zc)
      !!deallocate(zlold)
      !!deallocate(zcold)
      !!deallocate(xc)
      deallocate(xl)
      !!deallocate(xu)
      !!deallocate(yc)
      !!deallocate(yl)
      !!deallocate(yu)
      deallocate(dx)
      deallocate(dy)
      !!deallocate(dz)
      !!deallocate(xmet)
      !!deallocate(ymet)
      !!deallocate(zmet)
      !!deallocate(zmetl)
      deallocate(rlon)
      deallocate(rlat)
      !!deallocate(zl2)
      !!deallocate(zl3)
      !!deallocate(zc2)
      !!deallocate(zc3)
      !!deallocate(dz2)
      !!deallocate(dz3)
      !!deallocate(xc2)
      !!deallocate(xc3)
      !!deallocate(yc2)
      !!deallocate(yc3)
      !!deallocate(xl2)
      !!deallocate(xl3)
      !!deallocate(yl2)
      !!deallocate(yl3)
      !!deallocate(xu2)
      !!deallocate(xu3)
      !!deallocate(yu2)
      !!deallocate(yu3)
      deallocate(dx2)
      deallocate(dx3)
      !!deallocate(dy2)
      !!deallocate(dy3)
      !!deallocate(xmet2)
      !!deallocate(xmet3)
      !!deallocate(ymet2)
      !!deallocate(ymet3)
      !!deallocate(zmet2)
      !!deallocate(zmet3)
      !!deallocate(zmetl2)
      !!deallocate(zmetl3)
      !!deallocate(if_sec_mom)
      !!deallocate(if_nuc)
      deallocate(igelem)
      deallocate(ncore)
      deallocate(ienconc)
      !!deallocate(icoag)
      !!deallocate(icoagelem)
      !!deallocate(imomelem)
      !!deallocate(inucproc)
      !!deallocate(igrowgas)
      !!deallocate(inucgas)
      !!deallocate(nnuc2elem)
      !!deallocate(inuc2elem)
      !!deallocate(ievp2elem)
      !!deallocate(inuc2bin)
      !!deallocate(isolelem)
      !!deallocate(ievp2bin)
      !!deallocate(icorelem)
      !!deallocate(nnucelem)
      !!deallocate(nnucbin)
      !!deallocate(inucelem)
      !!deallocate(inucbin)
      deallocate(rmassmin)
      deallocate(r)
      deallocate(rmass)
      deallocate(rmasscore)
      deallocate(rcore)
      deallocate(vol)
      deallocate(dr)
      deallocate(dm)
      deallocate(dv)
      deallocate(rmassup)
      deallocate(rup)
      deallocate(rmasscoreup)
      deallocate(rcoreup)
      deallocate(rlow)
      deallocate(diffmass)
      !!deallocate(rhop)
      !!deallocate(rhopcore)
      !!deallocate(rhshell)
      deallocate(n)
      deallocate(rf)
      !!deallocate(ra)
      !!deallocate(rmc)
      !!deallocate(rmt)
      !!deallocate(rkna)
      !!deallocate(dfrac)
      deallocate(r0)
      deallocate(rsig)
      deallocate(totm)
      !DEALLOCATE(pc)
      !DEALLOCATE(gc)
      !!deallocate(ptc)
      !DEALLOCATE(pc2)
      !DEALLOCATE(pc3)
      !DEALLOCATE(pc4)
      !DEALLOCATE(pc5)
      !!deallocate(rhop2)
      deallocate(rhop3)
      deallocate(rhopcore3)
      !ALLOCATABLE(gc2)
      !ALLOCATABLE(gc3)
      !!deallocate(ptc2)
      !!deallocate(ptc3)
      !!deallocate(pcl)
      !!deallocate(gcl)
      !!deallocate(ptcl)
      !!deallocate(d_pc)
      !!deallocate(d_gc)
      !!deallocate(d_ptc)
      !!deallocate(pcmax)
      !!deallocate(cvert)
      !!deallocate(divcor)
      !!deallocate(chor)
      !!deallocate(dhor)
      !!deallocate(pconmax)
      !!deallocate(coaglg)
      !!deallocate(coagpe)
      !!deallocate(rnuclg)
      !!deallocate(rnucpe)
      !!deallocate(growlg)
      !!deallocate(growpe)
      !!deallocate(evaplg)
      !!deallocate(evappe)
      !!deallocate(gasprod)
      !!deallocate(vertdifd)
      !!deallocate(vertdifu)
      !!deallocate(ftopgas)
      !!deallocate(fbotgas)
      !!deallocate(ftoppart)
      !!deallocate(fbotpart)
      !!deallocate(pc_topbnd)
      !!deallocate(pc_botbnd)
      !!deallocate(gc_topbnd)
      !!deallocate(gc_botbnd)
      !!deallocate(ptc_topbnd)
      !!deallocate(ptc_botbnd)
      !!deallocate(cmf)
      !!deallocate(totevap)
      !!deallocate(inucmin)
      !!deallocate(inucstep)
      !!deallocate(evappe5)
      !!deallocate(cbr)
      !!deallocate(ccd)
      !!deallocate(cgr)
      !!deallocate(tim)
      !!deallocate(tsc)
      !!deallocate(ckernel)
      !!deallocate(pkernel)
      !!deallocate(volx)
      !!deallocate(ilow)
      !!deallocate(jlow)
      !!deallocate(iup)
      !!deallocate(jup)
      !!deallocate(npairl)
      !!deallocate(npairu)
      !!deallocate(iglow)
      !!deallocate(jglow)
      !!deallocate(igup)
      !!deallocate(jgup)
      !!deallocate(rkn)
      !!deallocate(bpm)
      !!deallocate(bpma)
      !!deallocate(vf)
      !!deallocate(vtrans)
      !!deallocate(re)
      !!deallocate(vertadvu)
      !!deallocate(vertadvd)
      !!deallocate(htrans)
      !!deallocate(hdiff)
      !!deallocate(ca)
      !!deallocate(cb)
      !!deallocate(cd)
      !!deallocate(ce)
      !!deallocate(cf)
      !!deallocate(cg)
      !!deallocate(rmfp)
      !DEALLOCATE(rhostar2)
      !DEALLOCATE(rhostar3)
      !DEALLOCATE(p)
      !DEALLOCATE(rhoa)
      !DEALLOCATE(p_surf)
      !DEALLOCATE(p_top)
      !DEALLOCATE(t)
      !!deallocate(told)
      !!deallocate(pold)
      !!deallocate(rhoaold)
      !!deallocate(rmu)
      !!deallocate(thcond)
      deallocate(w)
      !DEALLOCATE(zmetold)
      deallocate(u)
      deallocate(v)
      deallocate(t_surf)
      !!deallocate(dkz)
      !!deallocate(dkx)
      !!deallocate(dky)
      !!deallocate(p2)
      !!deallocate(p3)
      !!deallocate(t2)
      !!deallocate(t3)
      !!deallocate(u2)
      !!deallocate(u3)
      !!deallocate(v2)
      !!deallocate(v3)
      !!deallocate(w2)
      !!deallocate(w3)
      !!deallocate(dkx2)
      !!deallocate(dkx3)
      !!deallocate(dky2)
      !!deallocate(dky3)
      !!deallocate(dkz2)
      !!deallocate(dkz3)
      !!deallocate(rhoa2)
      !!deallocate(rhoa3)
      !DEALLOCATE(p_surf2)
      !!deallocate(p_top2)
      !!deallocate(gwtmol)
      !!deallocate(diffus)
      !!deallocate(rlhe)
      !!deallocate(rlhm)
      !!deallocate(pvapl)
      !!deallocate(pvapi)
      !!deallocate(surfctwa)
      !!deallocate(surfctiw)
      !!deallocate(surfctia)
      !!deallocate(akelvin)
      !!deallocate(akelvini)
      !!deallocate(ft)
      !!deallocate(gro)
      !!deallocate(gro1)
      !!deallocate(gro2)
      !!deallocate(gvrat)
      !!deallocate(supsatl)
      !!deallocate(supsati)
      !!deallocate(supsatlold)
      !!deallocate(supsatiold)
      !!deallocate(scrit)
      !!deallocate(sol_ions)
      !!deallocate(solwtmol)
      !!deallocate(rhosol)
      !!deallocate(rlh_nuc)
      !!deallocate(pvapl2)
      !!deallocate(pvapl3)
      !!deallocate(pvapi2)
      !!deallocate(pvapi3)
      !!deallocate(supsatl2)
      !!deallocate(supsatl3)
      !!deallocate(supsati2)
      !!deallocate(supsati3)
      !!deallocate(radheat2)
      !!deallocate(radheat3)
      !!deallocate(qrad2)
      !!deallocate(qrad3)
      !!deallocate(alb_tomi2)
      !!deallocate(alb_toai2)
      !!deallocate(alb_toa2)
      !!deallocate(opd2)
      !!deallocate(fsl_up2)
      !!deallocate(fsl_up3)
      !!deallocate(fsl_dn2)
      !!deallocate(fsl_dn3)
      !!deallocate(fir_up2)
      !!deallocate(fir_up3)
      !!deallocate(fir_dn2)
      !!deallocate(fir_dn3)
      deallocate(z_sin)
      deallocate(z_cos)
      deallocate(wave)
      deallocate(qrad)
      !!deallocate(radheat)
      deallocate(alb_tomi)
      deallocate(alb_toai)
      deallocate(alb_toa)
      deallocate(opd)
      !!deallocate(fsl_up)
      !!deallocate(fsl_dn)
      !!deallocate(fir_up)
      !!deallocate(fir_dn)

    end subroutine final_definitions_globaer

end module mem_globaer
