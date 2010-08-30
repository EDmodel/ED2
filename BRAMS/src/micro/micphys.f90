!====================================== Change Log ========================================!
! 5.0.0                                                                                    !
!                                                                                          !
!==========================================================================================!
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved                       !
!  Regional Atmospheric Modeling System - RAMS                                             !
!==========================================================================================!
!==========================================================================================!
module micphys

   use grid_dims, only: nzpmax, maxgrds

   !---------------------------------------------------------------------------------------!
   !     The product [(nthz-1)  * dthz ] must equal 25.0.                                  !
   !     The product [(nrhhz-1) * drhhz] must equal 0.18.                                  !
   !     The product [(ntc-1)   * dtc  ] must equal 20.0.                                  !
   !     The product [(ndnc-1)  * ddnc ] must equal 20.e-6.                                !
   !---------------------------------------------------------------------------------------!
   integer, parameter ::    &
                    nthz     =   26  & ! # of temp values spanning haze nucleation table
                   ,nrhhz    =   10  & ! # of R.H. values spanning haze nucleation table
                   ,ngam     = 8000  & ! # of values in incomplete gamma function table
                   ,ninc     =  201  & ! 
                   ,ndns     =   15  & ! 
                   ,ntc      =   21  & ! 
                   ,ndnc     =   11  & ! 
                   ,nd1cc    =   30  & ! Dimension of cloud-cloud 
                   ,nd1cr    =   15  & ! 
                   ,nr2cr    =   10  & ! 
                   ,nd2cr    =   30  & ! 
                   ,nr2rr    =   20  & ! 
                   ,nd2rr    =   20  & ! 
                   ,ncat     =    7  & ! # of hydrometeor categories
                   ,nhcat    =   15  & ! # of hydrometeor categories including ice habits
                   ,npairc   =   93  & ! # of pairs of species in number collection table
                   ,npairr   =  131  & ! # of pairs of species in mass collection table
                   ,nembc    =   20  & !
                   ,nembfall =   20  & !
                   ,maxkfall =    4  ! !


   real           , parameter ::      & !
                    dtc      = 1.     & !
                   ,ddnc     = 2.e-6  & !
                   ,dthz     = 1.     & !
                   ,drhhz    =  .02   ! !

   real(kind=8)   , parameter ::           & !
                    dtc8     = dble(dtc)   & !
                   ,ddnc8    = dble(ddnc)  & !
                   ,dthz8    = dble(dthz)  & !
                   ,drhhz8   = dble(drhhz) ! !


   !---------------------------------------------------------------------------------------!
   !    Namelist-related variables                                                         !
   !---------------------------------------------------------------------------------------!
   integer :: icloud   ! Method to determine cloud
   integer :: irain    ! Method to determine rain
   integer :: ipris    ! Method to determine pristine ice.
   integer :: isnow    ! Method to determine snow.
   integer :: iaggr    ! Method to determine aggregates.
   integer :: igraup   ! Method to determine graupel.
   integer :: ihail    ! Method to determine hail.
   integer :: mkcoltab ! Use collection table or generate a new one. 
   integer :: lpw      ! Lowest point in W.
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !    Nucleation levels.                                                                 !
   !---------------------------------------------------------------------------------------!
   integer :: k1cnuc,k2cnuc,k1pnuc,k2pnuc
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !    Flags for hydrometeor attitudes.                                                   !
   !---------------------------------------------------------------------------------------!
   integer, dimension(ncat)        :: jnmb              ! Integer flag
   logical, dimension(ncat)        :: availcat,progncat ! Logical flag.

   integer, dimension(10)          :: k1,k2,k3          ! Lower and upper levels
   integer, dimension(nhcat,nhcat) :: ipairc,ipairr
   integer, dimension(31,100,2)    :: jhabtab
   integer, dimension(nzpmax,ncat) :: jhcat,ict1,ict2

   !---------------------------------------------------------------------------------------!
   !   Variables that will be explained at some point...                                   !
   !---------------------------------------------------------------------------------------!
   real                                :: cparm,rparm,pparm,sparm,aparm,gparm,hparm
   real                                :: rictmin,rictmax,dps,dps2,d1min,r2min,d2min,d1max
   real                                :: r2max,d2max,d1ecc,d1ecr,r2ecr,r2err,colf,pi4dt
   real                                :: sedtime0,sedtime1

   real, dimension(ncat)               :: emb0,emb1,gnu,parm,emb0log,emb1log,dict,rxmin

   real, dimension(nhcat)              :: shapefac,cfmas,pwmas,cfvt,pwvt,dpsmi,cfden,pwden
   real, dimension(nhcat)              :: cfemb0,cfen0,pwemb0,pwen0,vtfac,frefac1,frefac2
   real, dimension(nhcat)              :: cfmasft,dnfac,sipfac,cfmasi,pwmasi,ch1,ch3,cdp1
   real, dimension(nhcat)              :: pwvtmasi,emb2,cxmin

   real, dimension(nzpmax)             :: tair,tairc,tairstr,til,rvstr,press,exner
   real, dimension(nzpmax)             :: rhoa,rhoi,rtot,rvap,rliq,rice,qhydm
   real, dimension(nzpmax)             :: rvlsair,rvisair,thrmcon
   real, dimension(nzpmax)             :: vapdif,dynvisc,rdynvsci,denfac,dn0i,colfacr
   real, dimension(nzpmax)             :: colfacr2,colfacc,colfacc2,sumuy,sumuz,sumvr
   real, dimension(nzpmax)             :: scrmic1,scrmic2,scrmic3,cccnx,cifnx
   real, dimension(nzpmax)             :: dsed_thil,totcond,thil,pottemp,vertvelo,rloss
   real, dimension(nzpmax)             :: enloss,rfall,cfall,qrfall,theiv
   real, dimension(nzpmax,ncat)        :: rx,cx,qr,qx,tx,emb,vterm,vap,ttest,wct1
   real, dimension(nzpmax,ncat)        :: wct2,sb,sd,se,sf,sg,sh,sm,ss,su,sw,sy,sz

   real, dimension(nzpmax,2)           :: tref,rvsref,rvsrefp
   real, dimension(nzpmax,9)           :: sa
   real, dimension(nzpmax,10)          :: eff

   real, dimension(nzpmax,ncat,ncat)   :: rxfer,qrxfer,enxfer
   real, dimension(nhcat)              :: dispemb0,dispemb0i,dispemb1,ch2

   real, dimension(nembc,nembc,npairc) :: coltabc
   real, dimension(nembc,nembc,npairr) :: coltabr

   real, dimension(nrhhz,nthz)         :: frachz
   real, dimension(ndnc,ntc,maxgrds)   :: fracc
   real, dimension(4)                  :: gamm,gamn1
   real, dimension(ngam,3)             :: gam
   real, dimension(ngam,2)             :: gaminc
   real, dimension(ngam)               :: gamsip13,gamsip24
   real, dimension(ninc)               :: rmlttab
   real, dimension(ninc,nhcat)         :: enmlttab
   real, dimension(ninc,ndns)          :: shedtab
   real, dimension(2)                  :: sc,sl,sq
   real, dimension(7)                  :: sj,pcprx,accpx

   real, dimension(nd1cc)              :: r1tabcc,c1tabcc,c2tabcc
   real, dimension(nd1cr,nr2cr,nd2cr)  :: r1tabcr,c1tabcr
   real, dimension(nr2rr,nd2rr)        :: c2tabrr

   character(len=256)                  :: coltabfn


end module micphys
!==========================================================================================!
!==========================================================================================!
