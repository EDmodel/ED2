!==========================================================================================!
!==========================================================================================!
!   module micro_coms.f90 (former micphys.h).  This module contains the main microphysics  !
! parameters and dimensions for Ramspost.                                                  !
!------------------------------------------------------------------------------------------!
module micro_coms
   use rpost_dims, only : str_len
   !---------------------------------------------------------------------------------------!
   !     The product [(nthz-1)  * dthz ] must equal 25.0.                                  !
   !     The product [(nrhhz-1) * drhhz] must equal 0.18.                                  !
   !     The product [(ntc-1)   * dtc  ] must equal 20.0.                                  !
   !     The product [(ndnc-1)  * ddnc ] must equal 20.e-6.                                !
   !---------------------------------------------------------------------------------------!

   integer, parameter :: nthz   =   26
   integer, parameter :: nrhhz  =   10
   integer, parameter :: ngam   = 5000
   integer, parameter :: ninc   =  201
   integer, parameter :: ndns   =   15
   integer, parameter :: ntc    =   21
   integer, parameter :: ndnc   =   11
   integer, parameter :: nd1cc  =   30
   integer, parameter :: nd1cr  =   15
   integer, parameter :: nr2cr  =   10
   integer, parameter :: nd2cr  =   30
   integer, parameter :: nr2rr  =   20
   integer, parameter :: nd2rr  =   20
   integer, parameter :: nccn   =    6
   integer, parameter :: nak    =   10
   integer, parameter :: ncc    =    7
   integer, parameter :: nsup   =   11
   integer, parameter :: ntemp  =   16
   integer, parameter :: ncat   =    7
   integer, parameter :: nhcat  =   15
   integer, parameter :: npairc =   93
   integer, parameter :: npairr =  131
   integer, parameter :: nembc  =   20
   real   , parameter :: dtc    =    1.
   real   , parameter :: ddnc   =    2.e-6
   real   , parameter :: dthz   =    1.
   real   , parameter :: drhhz  =    0.02

   !---------------------------------------------------------------------------------------!
   !     Variables to be filled by rcio.f90.                                               !
   !---------------------------------------------------------------------------------------!
   integer                :: iccnflg
   integer                :: ifnflg
   integer                :: icloud
   integer                :: irain
   integer                :: ipris
   integer                :: isnow
   integer                :: iaggr
   integer                :: igraup
   integer                :: ihail
   integer                :: mkcoltab

   real                   :: cparm
   real                   :: rparm
   real                   :: pparm
   real                   :: sparm
   real                   :: aparm
   real                   :: gparm
   real                   :: hparm
   real                   :: rictmin
   real                   :: rictmax
   real                   :: dps
   real                   :: dps2
   real                   :: d1min
   real                   :: r2min
   real                   :: d2min
   real                   :: d1max
   real                   :: r2max
   real                   :: d2max
   real                   :: d1ecc
   real                   :: d1ecr
   real                   :: r2ecr
   real                   :: r2err
   real                   :: colf
   real                   :: pi4dt
   real                   :: sedtime0
   real                   :: sedtime1

   real, dimension(nhcat) :: shape
   real, dimension(nhcat) :: cfmas
   real, dimension(nhcat) :: pwmas
   real, dimension(nhcat) :: cfvt
   real, dimension(nhcat) :: pwvt
   real, dimension(nhcat) :: dpsmi
   real, dimension(nhcat) :: cfden
   real, dimension(nhcat) :: pwden
   real, dimension(nhcat) :: cfemb0
   real, dimension(nhcat) :: cfen0
   real, dimension(nhcat) :: pwemb0
   real, dimension(nhcat) :: pwen0
   real, dimension(nhcat) :: vtfac
   real, dimension(nhcat) :: frefac1
   real, dimension(nhcat) :: frefac2
   real, dimension(nhcat) :: cfmasft
   real, dimension(nhcat) :: dnfac
   real, dimension(nhcat) :: sipfac
   real, dimension(nhcat) :: pwmasi
   real, dimension(nhcat) :: ch1
   real, dimension(nhcat) :: ch3
   real, dimension(nhcat) :: cdp1
   real, dimension(nhcat) :: pwvtmasi

   character(len=str_len) :: coltabfn
end module micro_coms
!==========================================================================================!
!==========================================================================================!







