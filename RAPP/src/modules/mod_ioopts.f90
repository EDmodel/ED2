!------------------------------------------------------------------------------------------!
!  Module ioopts: contains the list of input/output basic settings.                        !
!------------------------------------------------------------------------------------------!
module mod_ioopts
   use mod_maxdims , only : maxstr   & ! intent(in)
                          , maxfiles & ! intent(in)
                          , maxgrds  & ! intent(in)
                          , maxtimes ! ! intent(in)
   use mod_time    , only : time_stt
   implicit none
   
   !---------------------------------------------------------------------------------------!
   !    Namelist variables                                                                 !
   !---------------------------------------------------------------------------------------!
   !----- Paths, prefixes, and types of inputs and outputs --------------------------------!
   character(len=maxstr) :: intype  ! Type of input file you want
   character(len=maxstr) :: inpath  ! Path of the input dataset
   character(len=maxstr) :: outpref ! Prefix of output file

   !----- Time information ----------------------------------------------------------------!
   integer               :: iyeara  ! First year to process
   integer               :: iyearz  ! Last year to process
   real                  :: inpfrq  ! Input variable frequency
   real                  :: radfrq  ! Output frequency for radiation fluxes

   !----- Height information. -------------------------------------------------------------!
   real                  :: ref_hgt  ! Reference height.

   !----- Output edge information. --------------------------------------------------------!
   real                  :: lonw     ! Westernmost longitude
   real                  :: lone     ! Easternmost longitude
   real                  :: lats     ! Southernmost latitude
   real                  :: latn     ! Northernmost latitude
   integer               :: edgeoff  ! # of edge points not to be included in the output.
   !---------------------------------------------------------------------------------------!

   !----- Ratio between time variables ... ------------------------------------------------!
   integer               :: radratio    ! Inpfrq / Radfrq
   integer               :: ninputs_day ! 1 day  / Inpfrq
   integer               :: nsteps      ! Inpfrq / Dtinc
   integer               :: nrads       ! Radfrq / Dtinc
   real                  :: nstepsi     ! Dtinc  / Inpfrq
   real                  :: nradsi      ! Dtinc  / Radfrq
   !---------------------------------------------------------------------------------------!

   !----- Time map ------------------------------------------------------------------------!
   type(time_stt)        , dimension(maxtimes)                  :: outtimes
   integer                                                      :: nouttimes
   integer               , dimension(maxfiles,maxtimes,maxgrds) :: time_address
   !---------------------------------------------------------------------------------------!

   !----- Parameters ----------------------------------------------------------------------!
   integer               , parameter           :: nvars = 11 ! # of variables
   integer               , parameter           :: missflg_int   = 999999
   real                  , parameter           :: missflg_real  = -9.9692100E+36 !1.E+34
   character(len=maxstr) , parameter           :: missflg_char  = '___'
   real(kind=8)          , parameter           :: missflg_dble  = -9.96920996838687e+36 -9.96921E36 !1.D+34 
   integer               , parameter           :: rapp_absent   = -999
   integer               , parameter           :: rapp_badtype  = -888
   integer               , parameter           :: rapp_partly   = -777
   integer               , parameter           :: rapp_libvar   = -666
   integer               , parameter           :: rapp_notfound = -555
   !---------------------------------------------------------------------------------------!

end module mod_ioopts
