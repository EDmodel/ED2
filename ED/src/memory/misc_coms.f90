Module misc_coms


  use max_dims, only: str_len,maxpvars,str_len_short

  implicit none

  type simtime
     integer year
     integer month
     integer date
     real time
     integer ifirst
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

  real :: dtlsm
  real :: radfrq

  integer :: ifoutput
  integer :: idoutput
  integer :: imoutput
  integer :: iyoutput
  integer :: isoutput
  integer :: iclobber
  real ::  frqstate
  real ::  frqfast
  real ::  frqsum
  real ::  outstate
  real ::  outfast
  integer :: nrec_fast
  integer :: nrec_state
  integer :: irec_fast
  integer :: irec_state
  type(simtime) :: out_time_fast
  type(simtime) :: out_time_state

  character(len=str_len) :: sfilin
  character(len=str_len) ::ffilout 
  character(len=str_len) ::sfilout
  integer :: ied_init_mode
  
  character(len=str_len) :: ed_inputs_dir
  integer :: integration_scheme


  ! Control parameters for printing. Read in the namelist
  integer :: iprintpolys
  integer :: npvars
  character(len=str_len_short), dimension(maxpvars) :: printvars
  character(len=str_len_short), dimension(maxpvars) :: pfmtstr
  integer :: ipmin
  integer :: ipmax
  
  
  ! XML configuration file
  character(len=str_len) :: iedcnfgf

end Module misc_coms
