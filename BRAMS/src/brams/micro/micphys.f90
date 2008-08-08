!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################


Module micphys

  use grid_dims

  !--------------------------------------------------------------------------
  !     The product [(nthz-1)  * dthz ] must equal 25.0.
  !     The product [(nrhhz-1) * drhhz] must equal 0.18.
  !     The product [(ntc-1)   * dtc  ] must equal 20.0.
  !     The product [(ndnc-1)  * ddnc ] must equal 20.e-6.

  integer, parameter :: nthz=26,nrhhz=10,ngam=5000,ninc=201  &
       ,ndns=15,ntc=21,ndnc=11  &
       ,nd1cc=30,nd1cr=15,nr2cr=10,nd2cr=30,nr2rr=20  &
       ,nd2rr=20  &
       ,ncat=7,nhcat=15,npairc=93,npairr=131,nembc=20
  real, parameter    :: dtc=1.,ddnc=2.e-6 ,dthz=1.,drhhz=.02
  !--------------------------------------------------------------------------
  integer :: level,icloud,irain,ipris,isnow,iaggr,igraup,ihail  &
       ,mkcoltab
  integer, dimension(ncat)        :: jnmb
  integer, dimension(nhcat,nhcat) :: ipairc,ipairr
  integer, dimension(31,100,2)    :: jhabtab
  integer, dimension(nzpmax,ncat) :: jhcat,ict1,ict2

  real :: cparm,rparm,pparm,sparm,aparm,gparm,hparm,rictmin,rictmax  &
       ,dps,dps2,d1min,r2min,d2min,d1max,r2max,d2max  &
       ,d1ecc,d1ecr,r2ecr,r2err,colf,pi4dt,sedtime0,sedtime1

  real, dimension(ncat)  :: emb0,emb1,gnu,parm,emb0log,emb1log,dict,rxmin

  ! Changing name of variabel SHAPE to VAR_SHAPE
  ! To avoid confusing with SHAPE FUNCTION (intrinsic)
  ! ALF
  real, dimension(nhcat) :: var_shape,cfmas,pwmas,cfvt,pwvt,dpsmi,cfden,pwden &
       ,cfemb0,cfen0,pwemb0,pwen0,vtfac,frefac1,frefac2  &
       ,cfmasft,dnfac,sipfac,pwmasi,ch1,ch3,cdp1,pwvtmasi,emb2

  real, dimension(nzpmax) :: tair,tairc,tairstrc,til,rvstr,press,pitot  &
       ,rliq,rice,qhydm,rvlsair,rvisair,rvs0,thrmcon  &
       ,vapdif,dynvisc,rdynvsci,denfac,dn0i,colfacr  &
       ,colfacr2,colfacc,colfacc2,sumuy,sumuz,sumvr  &
       ,scrmic1,scrmic2,scrmic3,cccnx,cifnx
  real, dimension(nzpmax,ncat) :: rx,cx,qr,qx,tx,emb,vterm,vap,ttest,wct1  &
       ,wct2,sb,sd,se,sf,sg,sh,sm,ss,su,sw,sy,sz
  real, dimension(nzpmax,2)  :: tref,rvsref,rvsrefp
  real, dimension(nzpmax,9)  :: sa
  real, dimension(nzpmax,10) :: eff

  real, dimension(nzpmax,ncat,ncat) :: rxfer,qrxfer,enxfer
  real, dimension(nhcat,maxgrds)    :: dispemb0,dispemb1,ch2

  real, dimension(nembc,nembc,npairc) :: coltabc
  real, dimension(nembc,nembc,npairr) :: coltabr

  real, dimension(nrhhz,nthz)       :: frachz
  real, dimension(ndnc,ntc,maxgrds) :: fracc
  real, dimension(4)                :: gamm,gamn1
  real, dimension(ngam,3)           :: gam
  real, dimension(ngam,2)           :: gaminc
  real, dimension(ngam)             :: gamsip13,gamsip24
  real, dimension(ninc)             :: rmlttab
  real, dimension(ninc,nhcat)       :: enmlttab
  real, dimension(ninc,ndns)        :: shedtab
  real, dimension(2)                :: sc,sk,sl
  real, dimension(7)                :: sj,pcprx,accpx

  real, dimension(nd1cc)             :: r1tabcc,c1tabcc,c2tabcc
  real, dimension(nd1cr,nr2cr,nd2cr) :: r1tabcr,c1tabcr
  real, dimension(nr2rr,nd2rr)       :: c2tabrr

  character(len=256) :: coltabfn
  ! Modif. by ALF

  !---------------------------------------------------------------------------


end Module micphys
