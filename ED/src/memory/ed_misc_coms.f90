Module ed_misc_coms

   use ed_max_dims, only: str_len,maxpvars,str_len_short,maxgrds

   implicit none

   type simtime
      integer :: year
      integer :: month
      integer :: date
      integer :: hour
      integer :: min
      integer :: sec
      real    :: time
      integer :: ifirst
   end type simtime
   type(simtime) :: current_time
   type(simtime) :: end_time

   character(len=str_len) :: expnme
   character(len=str_len) :: runtype

   integer :: itimea
   integer :: iyeara
   integer :: imontha
   integer :: idatea

   integer :: itimez
   integer :: iyearz
   integer :: imonthz
   integer :: idatez

   integer :: itimeh
   integer :: iyearh
   integer :: imonthh
   integer :: idateh

   real :: dtlsm
   real :: radfrq

   integer :: ifoutput
   integer :: idoutput
   integer :: imoutput
   integer :: iqoutput
   integer :: iyoutput
   integer :: itoutput
   integer :: isoutput
   integer :: iclobber

   integer :: unitfast
   integer :: unitstate
   real ::  frqstate
   real ::  frqfast
   real ::  frqsum
   real ::  outstate
   real ::  outfast

   integer :: ndcycle

   integer :: nrec_fast
   integer :: nrec_state
   integer :: irec_fast
   integer :: irec_state
   type(simtime) :: out_time_fast
   type(simtime) :: out_time_state

   character(len=str_len), dimension(maxgrds) :: sfilin
   character(len=str_len) ::ffilout 
   character(len=str_len) ::sfilout
   integer :: ied_init_mode
   
   character(len=str_len) :: thsums_database

   !---------------------------------------------------------------------------------------!
   !    Maximum distance to the current polygon that we still consider the file grid point !
   ! to be representative of the polygon for thermal sums.                                 !
   !---------------------------------------------------------------------------------------!
   real    :: max_thsums_dist


   !---------------------------------------------------------------------------------------!
   !      Alternative method for mixing 1 grid and POI's.  Only use the grid if their is   !
   ! NOT an POI  within a user specified resolution.  Remember, this assumes there is only !
   ! 1 gridded file, and it is the first file when ied_init_mode is set to 99  (Developer  !
   ! use only).                                                                            !
   !---------------------------------------------------------------------------------------!
   real    :: max_poi99_dist
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !      This variable is used for the history start initialisation.  This sets the       !
   ! maximum acceptable distance between the expected polygon and the polygon found in the !
   ! history file.  Units: m.                                                              !
   !---------------------------------------------------------------------------------------!
   real    :: max_poihist_dist
   !---------------------------------------------------------------------------------------!

   integer :: integration_scheme


   ! Control parameters for printing. Read in the namelist
   integer :: iprintpolys
   integer :: npvars
   character(len=str_len), dimension(maxpvars) :: printvars
   character(len=str_len), dimension(maxpvars) :: pfmtstr
   integer :: ipmin
   integer :: ipmax
   
   
   ! XML configuration file
   character(len=str_len) :: iedcnfgf

   ! XML event file
   character(len=str_len) :: event_file

   integer :: burnin          !! number of years to ignore demography when starting a run

   integer :: outputMonth     !! month to output annual files

   integer :: restart_target_year    !! year to read when parsing pss/css with multiple years

   integer :: use_target_year        !! flag specifying whether to search for a target year in pss/css

   ! flags to turn on/off intersite variability in edaphic variables
   ! easier than redefining site files
   integer :: vary_elev 
   integer :: vary_rad
   integer :: vary_hyd  

   ! soil biogeochem initial conditions (over-rides patch files)
   ! useful for data assimilation & sensitivity analysis
   real    :: init_fsc 
   real    :: init_stsc 
   real    :: init_ssc 
   real    :: init_stsl 
   real    :: init_fsn 
   real    :: init_msn 

   ! Logical Switches for various memory structures

   logical :: fast_diagnostics       !! If ifoutput,idoutput,and imoutput are zero, then
                                     !! there is no need to integrate fast flux diagnostics

   ! Namelist option to attach metadata to HDF5 output files 0=no, 1=yes

   integer :: attach_metadata
   
   !---------------------------------------------------------------------------------------!
   !     Age and Size classes.                                                             !
   !---------------------------------------------------------------------------------------!
   real    :: maxdbh ! Maximum DBH to be divided in classes 
   real    :: maxage ! Maximum age to be divided in classes
                     ! In both cases, if the value exceeds the maximum, they will all
                     !    go to the last class.

   real    :: ddbhi  ! Inverse of DBH class bin size
   real    :: dagei  ! Inverse of age class bin size.
   !---------------------------------------------------------------------------------------!


   !----- Minimum site area that we will allocate (unused if ied_init_mode is 3 or 4). ----!
   real    :: min_site_area
   !---------------------------------------------------------------------------------------!


   !----- Namelist option for allometry scheme. -------------------------------------------!
   integer :: iallom ! 0 -- Original ED-2.1 allometry
                     ! 1 -- DBH -> AGB Tree allometry based on Baker et al. (2004)
                     !      keep original ED-2.1 Bl/Bd ratio
                     ! 2 -- DBH -> AGB Tree allometry based on Baker et al. (2004)
                     !      keep original ED-2.1 Bl
                     ! 3 -- Same as 2, root profile as in Kenzo et al. (2008)
                     ! 4 -- Same as 2, root profile defined in a simple equation that
                     !      puts roots at 0.5 m when the height is 0.15m, and 5.0 m when
                     !      the height is 35.0m.
   !---------------------------------------------------------------------------------------!

end module ed_misc_coms
