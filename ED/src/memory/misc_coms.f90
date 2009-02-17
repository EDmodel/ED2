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

  ! XML event file
  character(len=str_len) :: event_file

  ! integrator error statistics
  integer(kind=8),dimension(1:50,1:2):: integ_err
  logical,parameter :: record_err = .true.

contains
  function err_label(i) result(lab)
    integer :: i
    character(10) :: lab
    select case(i)
    case(1)  
       lab = "can_temp"
    case(2)  
       lab = "can_shv"
    case(3)  
       lab = "can_co2"
    case(4)  
       lab = "soil_h20_1"
    case(5)  
       lab = "soil_h20_2"
    case(6)  
       lab = "soil_h20_3"
    case(7)  
       lab = "soil_h20_4"
    case(8)  
       lab = "soil_h20_5"
    case(9)  
       lab = "soil_h20_6"
    case(10) 
       lab = "soil_h20_7"
    case(11) 
       lab = "soil_h20_8"
    case(12) 
       lab = "soil_h20_9"
    case(13) 
       lab = "soil_h2010"
    case(14) 
       lab = "soil_h2011"
    case(15) 
       lab = "soil_h2012"
    case(16) 
       lab = "soil_en_01"
    case(17) 
       lab = "soil_en_02"
    case(18) 
       lab = "soil_en_03"
    case(19) 
       lab = "soil_en_04"
    case(20) 
       lab = "soil_en_05"
    case(21) 
       lab = "soil_en_06"
    case(22) 
       lab = "soil_en_07"
    case(23) 
       lab = "soil_en_08"
    case(24) 
       lab = "soil_en_09"
    case(25) 
       lab = "soil_en_10"
    case(26) 
       lab = "soil_en_11"
    case(27) 
       lab = "soil_en_12"
    case(28) 
       lab = "sfcW_en_01"
    case(29) 
       lab = "sfcW_en_02"
    case(30) 
       lab = "sfcW_en_03"
    case(31) 
       lab = "sfcW_en_04"
    case(32) 
       lab = "sfcW_en_05"
    case(33) 
       lab = "sfcW_mas_1"
    case(34) 
       lab = "sfcW_mas_2"
    case(35) 
       lab = "sfcW_mas_3"
    case(36) 
       lab = "sfcW_mas_4"
    case(37) 
       lab = "sfcW_mas_5"
    case(38) 
       lab = "virt_heat"
    case(39) 
       lab = "virt_water"
    case(40) 
       lab = "veg_water"
    case(41) 
       lab = "veg_energy"
    case(42) 
       lab = "soilTlfliq"
    case(43) 
       lab = "soilfraclq"
    case(44) 
       lab = "sfcH2OTemp"
    case(45) 
       lab = "soilTemp"
    case(46) 
       lab = "vegTemp"
    case default 
       lab = "unspec"
    end select
    
  end function err_label


end Module misc_coms
