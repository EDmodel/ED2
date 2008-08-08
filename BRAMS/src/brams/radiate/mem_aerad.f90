module mem_aerad
  !  @(#) aerad.h  McKie  Oct-1995
  !  This is the include file for the interface between the aerosol
  !  microphysics and radiative transfer components of CARMA.
 
  !  Global symbolic constants are defined and common blocks are
  !  declared.
 
  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
 
  use mem_grid_dim_defs, only: maxz ! intent(in)

  !  Start of user-defined symbolic constants
 
  !  Define # grid pts in x, y, z directions
  integer,parameter :: nx =  1 
  integer,parameter :: ny =  1 
  !	 parameter( NZ = 1 )
  !	 parameter( NZ = NZPMAX )  >>>>>>>>>> no futuro mude isto
!  INTEGER,PARAMETER :: nz = 32 !versao Q2002,Q2003 e Q2004
!  INTEGER,PARAMETER :: nz = 39  !nz=nzp-1
!  INTEGER,PARAMETER :: nz = 42 ! Versao operacional atual
  integer :: nz != maxz-1
  !  Define maximum of NX or NY
  integer,parameter :: nxorny = nx
  !  Define # x, y direction grid box boundaries
  integer,parameter :: nxp1 = nx + 1
  integer,parameter :: nyp1 = ny + 1
  !  Define maximum of NXP1 or NYP1
  integer,parameter :: nxornyp1 = nxp1
  !  Define # particle radius bins
  integer,parameter :: nbin = 30
  !  Define # particle elements
  integer,parameter :: nelem = 1
  !  Define # particle groups
  integer,parameter :: ngroup = 1
  !  Define # solutes
  integer,parameter :: nsolute = 1
  !  Define # gases
  integer,parameter :: ngas = 1
  !  Define # solar wavelength bins
  !	parameter( NSOL = 26 )
  integer,parameter :: nsol = 32
  !  Define # infrared wavelength bins
  integer,parameter :: nir = 18
  !  Define total # wavelength bins
  integer,parameter :: nwave = nsol + nir
  !  Define # layers in rad xfer model domain underlying aerosol model domain
  integer,parameter :: nz_below = 0
  !  End of user-defined symbolic constants
 
  !  The remaining symbolic constants will need no attention from most
  !  users
 
  !  Define # layers in radiation model
!  INTEGER,PARAMETER :: nz_rad = nz + nz_below 
  integer :: nz_rad
  !  Define # vertical grid boundaries
!  INTEGER,PARAMETER :: nzp1 = nz + 1 
  integer :: nzp1
  !  Define logical unit number for output print file
  integer,parameter :: lunoprt = 10 
  !  Define logical unit number for time step info output
  integer,parameter :: lunostep = 11 
  !  Define logical unit number for input restart file
  integer,parameter :: lunires = 12 
  !  Define logical unit number for output restart file
  integer,parameter :: lunores = 13 
  !  Define logical unit number for output history file
  integer,parameter :: lunohis = 14 
  !  Define logical unit number for input and output of Mie coefficients
  integer,parameter :: lunmie = 15 
  !  Define logical unit number for print output from radiation submodel
  integer,parameter :: lunorad = 16 
!!!!!!!!!/kml
  integer,parameter :: iprocopio = 1 
!!!!!!!!!/kml
 
  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
 
  !  Declare common blocks for input to radiative transfer code
 
  !   is_grp_ice     =.true. means group is ice crystals
  !   r_aerad	     radius mid-pts from aerosol grids [cm]
  !   rcore_aerad    core radius mid-pts from aerosol grids [cm]
  !   rup_aerad      upper radii from aerosol grids [cm]
  !   rcoreup_aerad  upper core radius  from aerosol grids [cm]
  !   p_aerad	     pressure [dyne/cm^2]
  !   t_aerad	     temperature [K]
  !   pc_aerad       particle concentration [#/cm^3]
  !   qv_aerad       water vapor mixing ratio [g/g]
  !   tabove_aerad   blackbody temperature for downwelling IR flux into model [K]
  !   ptop_aerad     pressure at top of aerosol model domain [dyne/cm^2]
  !   pbot_aerad     pressure at bottom of aerosol model domain [dyne/cm^2]
  !   u0_aerad       cosine of solar zenith angle [dimensionless]
  !   sfc_alb_aerad  surface albedo [dimensionless]
  !   emisir_aerad   surface IR emissivity [dimensionless]
  !   tsfc_aerad     surface temperature [K]
  !   h2ocol_aerad   water vapor column above model domain [g/cm^2]
  !   isl_aerad      =1 means do solar calculations
  !   ir_aerad       =1 means do infrared calculations
  !   do_below       =.true. means include radiative layers below model domain
  !   ir_above_aerad =1 means include downwelling flux into top of model domain
  
  real, allocatable :: r_aerad(:,:)
  real, allocatable :: rcore_aerad(:,:)
  real, allocatable :: rup_aerad(:,:)
  real, allocatable :: rcoreup_aerad(:,:)
  !REAL, allocatable :: p_aerad(:)
  !REAL, allocatable :: pc_aerad(:,:,:)  
  !REAL, allocatable :: t_aerad(:)
  real, allocatable :: c_aerad(:,:,:)
  !REAL, allocatable :: qv_aerad(:)
  real :: tabove_aerad
  real :: ptop_aerad
  real :: pbot_aerad
  real :: u0_aerad
  real :: sfc_alb_aerad
  real :: emisir_aerad
  real :: tsfc_aerad
  real :: h2ocol_aerad
  logical, allocatable :: is_grp_ice_aerad(:)
  logical :: do_below
  real :: sl_aerad
  integer :: ir_aerad
  integer :: ir_above_aerad
  !INTEGER :: isl_aerad
  integer :: iaerad1
 
 
  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
 
  !  Declare common blocks for output from radiative transfer code
 
  !   heati_aerad    infrared heating rates [K/s]
  !   heats_aerad    solar heating rates [K/s]
  !   qrad_aerad     particle radiative heating rates [K/s]
  !   alb_tomi_aerad spectrally-integrated albedo at top-of-model
  !   alb_toai_aerad spectrally-integrated albedo at top-of-atmosphere
  !   alb_toa_aerad  spectrally-resolved albedo at top-of-atmosphere
  !   opd_aerad      spectrally-resolved optical depth
  !   fsl_up_aerad   solar upwelling flux [W m^-2]
  !   fsl_dn_aerad   solar downwelling flux [W m^-2]
  !   fir_up_aerad   infrared upwelling flux [W m^-2]
  !   fir_dn_aerad   infrared downwelling flux [W m^-2]
 
  real, allocatable :: wave_aerad(:)
  !REAL, allocatable :: heati_aerad(:)
  !REAL, allocatable :: heats_aerad(:)
  real, allocatable :: qrad_aerad(:,:,:)
  real :: alb_tomi_aerad
  real :: alb_toai_aerad
  real, allocatable :: alb_toa_aerad(:)
  real, allocatable :: opd_aerad(:)
  real, allocatable :: fsl_up_aerad(:)
  real, allocatable :: fsl_dn_aerad(:)
  real, allocatable :: fir_up_aerad(:)
  real, allocatable :: fir_dn_aerad(:)

  integer :: iaerad2

contains

  subroutine initial_definitions_aerad()

    implicit none

    ! Defining variables
    nz = maxz-1
    nz_rad = nz + nz_below
    nzp1 = nz + 1

    ! Allocating arrays

    allocate(r_aerad(nbin,ngroup))
    allocate(rcore_aerad(nbin,ngroup))
    allocate(rup_aerad(nbin,ngroup))
    allocate(rcoreup_aerad(nbin,ngroup))
    !ALLOCATE(p_aerad(nz_rad))
    !ALLOCATE(pc_aerad(nz_rad,nbin,ngroup))
    !ALLOCATE(t_aerad(nz_rad))
    allocate(c_aerad(nz_rad,nbin,ngroup))
    !allocate(qv_aerad(nz_rad))
    allocate(is_grp_ice_aerad(ngroup))
    allocate(wave_aerad(nwave+1))
    !allocate(heati_aerad(nz_rad))
    !allocate(heats_aerad(nz_rad))
    allocate(alb_toa_aerad(nsol))
    allocate(opd_aerad(nwave))
    allocate(qrad_aerad(nbin,nz_rad,ngroup))
    allocate(fsl_up_aerad(nz_rad+1))
    allocate(fsl_dn_aerad(nz_rad+1))
    allocate(fir_up_aerad(nz_rad+1))
    allocate(fir_dn_aerad(nz_rad+1))

  end subroutine initial_definitions_aerad

  ! **************************************************************************

  subroutine final_definitions_aerad

    implicit none

    ! Deallocating arrays
    deallocate(r_aerad)
    deallocate(rcore_aerad)
    deallocate(rup_aerad)
    deallocate(rcoreup_aerad)
    !deALLOCATE(p_aerad)
    !deALLOCATE(pc_aerad)
    !deALLOCATE(t_aerad)
    deallocate(c_aerad)
    !deallocate(qv_aerad)
    deallocate(is_grp_ice_aerad)
    deallocate(wave_aerad)
    !deallocate(heati_aerad)
    !deallocate(heats_aerad)
    deallocate(alb_toa_aerad)
    deallocate(opd_aerad)
    deallocate(qrad_aerad)
    deallocate(fsl_up_aerad)
    deallocate(fsl_dn_aerad)
    deallocate(fir_up_aerad)
    deallocate(fir_dn_aerad)

  end subroutine final_definitions_aerad

end module mem_aerad
