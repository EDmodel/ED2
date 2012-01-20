module mem_globaer

  use mem_aerad, only: &
       nelem,    &  !INTENT(IN)
       ngroup,   &  !INTENT(IN)
       ngas,     &  !INTENT(IN)
       nsolute      !INTENT(IN)
  use grid_dims, only : str_len
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
  
  
  
  !  Define core fraction (for core mass and second moment) used
  !  when particle number concentrations are limited to SMALL_PC
  
  real, parameter :: fix_coref = 0.01 
  
  !  Define small particle number concentration [ # / x_units / y_units / z_units ]
  real, parameter :: small_pc = 1.e-20
  
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
  
  character (LEN=str_len) :: gridname
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
  character (LEN=str_len) :: simtitle

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

  integer :: iaer6
  
  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$&$$$$$$$$$$$$$
  
  !  Declare global common blocks for coagulation group pair mapping
  
  !   iaer7      Safety marker for common block aer7

  integer :: iaer7
  
  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  
  !  Declare global common blocks for particle fall velocities, transport
  !  rates, and coagulation kernels
  !   vf_const  Constant vertical fall velocity when ifall=0   {setupaer}
  !   iaer8     Safety marker for common block aer8

  real :: vf_const

  integer :: iaer8
  
  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  
  !  Declare global common blocks for atmospheric structure.
  !  Note: air density, winds, and diffusion coefficients are in scaled units
  
  !   zbot      height of the bottom of the model [cm]      {initatm}
  !   t_surf    Air temperature at surface [deg_K]        {initatm}
  !   w     Vertical wind speed at layer boundary [z_units/s]     {initatm}
  !   u     East-west wind speed at layer center [x_units/s]     {initatm}
  !   v     North-south wind speed at layer center [y_units/s]      {initatm}
  !   iaer9     Safety marker for common block aer9
  
  real :: zbot
  real, allocatable :: w(:,:,:)
  real, allocatable :: u(:,:,:)
  real, allocatable :: v(:,:,:)
  real, allocatable :: t_surf(:,:)

  integer :: iaer9  
  integer :: iaer10
 
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
  real, allocatable :: alb_tomi(:,:)
  real, allocatable :: alb_toai(:,:)
  real, allocatable :: alb_toa(:,:,:)
  real, allocatable :: opd(:,:,:)

  logical :: do_rad
  logical :: do_solar
  logical :: do_ir
  integer :: nrad
  real :: prad
  integer :: irad1
    
  !  Declare alias names for radiative transfer parameters
  !  with first 2, 3 dimensions treated linearly
  

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
      allocate(xl(nx,ny,nz))
      allocate(dx(nx,ny,nz))
      allocate(dy(nx,ny,nz))
      allocate(rlon(nx,ny))
      allocate(rlat(nx,ny))
      allocate(dx2(nxy,nz))
      allocate(dx3(nxyz))
      allocate(igelem(nelem))
      allocate(ncore(ngroup))
      allocate(ienconc(ngroup))
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
      allocate(n(nbin,ngroup))
      allocate(rf(nbin,ngroup))
      allocate(r0(nx,ny,nz))
      allocate(rsig(nx,ny,nz))
      allocate(totm(nx,ny,nz))
      allocate(rhop3(nxyz,nbin,ngroup))
      allocate(rhopcore3(nxyz,nbin,ngroup))
      allocate(w(nx,ny,nzp1))
      allocate(u(nx,ny,nz))
      allocate(v(nx,ny,nz))
      allocate(t_surf(nx,ny))
      allocate(z_sin(nx,ny))
      allocate(z_cos(nx,ny))
      allocate(wave(nwave))
      allocate(qrad(nx,ny,nz,nbin,ngroup))
      allocate(alb_tomi(nx,ny))
      allocate(alb_toai(nx,ny))
      allocate(alb_toa(nx,ny,nsol))
      allocate(opd(nx,ny,nwave))

    end subroutine initial_definitions_globaer

    ! *************************************************************************

    subroutine final_definitions_globaer()
      implicit none

      ! Arrays allocation
      deallocate(xl)
      deallocate(dx)
      deallocate(dy)
      deallocate(rlon)
      deallocate(rlat)
      deallocate(dx2)
      deallocate(dx3)
      deallocate(igelem)
      deallocate(ncore)
      deallocate(ienconc)
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
      deallocate(n)
      deallocate(rf)
      deallocate(r0)
      deallocate(rsig)
      deallocate(totm)
      deallocate(rhop3)
      deallocate(rhopcore3)
      deallocate(w)
      deallocate(u)
      deallocate(v)
      deallocate(t_surf)
      deallocate(z_sin)
      deallocate(z_cos)
      deallocate(wave)
      deallocate(qrad)
      deallocate(alb_tomi)
      deallocate(alb_toai)
      deallocate(alb_toa)
      deallocate(opd)

    end subroutine final_definitions_globaer

end module mem_globaer
