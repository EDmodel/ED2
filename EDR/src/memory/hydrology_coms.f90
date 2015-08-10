Module hydrology_coms


  ! See "initialize_hydro_coms" for default values

  implicit none

  integer :: useRUNOFF
  integer :: useTOPMODEL

  integer :: HydroOutputPeriod !! multiples of dtlsm

  !! tuning parameter for TOPMODEL moisture redistribution 
  !! characteristic length (tau)
  real :: MoistRateTuning

  !! threshold for considering a soil layer saturated, units: m3/m3
  real :: MoistSatThresh

  ! For estimating water table depth, constant multiplied by 1/dt 
  ! to set maximum rate of change in watertable estimate
  real :: Moist_dWT

  ! threshold proportion of surface water that must be liquid for surface runoff to be computed
  real :: FracLiqRunoff

  ! maximum allowed overland runoff velocity (m/s)
  real :: runoff_vmax

  ! reference value for grass surface roughness
  real :: GrassLAImax

  real :: inverse_runoff_time

  !! defined elsewhere: infiltration_method, dewmax

end Module hydrology_coms
