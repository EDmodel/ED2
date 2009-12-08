!f90
!############################# Change Log ##################################
! 4.3.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!  Mission Research Corporation / *ASTeR Division
!###########################################################################

!--------------------------------------------------------------------------
!     The product [(nthz-1)  * dthz ] must equal 25.0.
!     The product [(nrhhz-1) * drhhz] must equal 0.18.
!     The product [(ntc-1)   * dtc  ] must equal 20.0.
!     The product [(ndnc-1)  * ddnc ] must equal 20.e-6.

integer, parameter :: nthz=26,nrhhz=10,ngam=5000,ninc=201  &
                     ,ndns=15,ntc=21,ndnc=11  &
                     ,nd1cc=30,nd1cr=15,nr2cr=10,nd2cr=30,nr2rr=20  &
                     ,nd2rr=20,nccn=6,nak=10,ncc=7,nsup=11,ntemp=16  &
                     ,ncat=7,nhcat=15,npairc=93,npairr=131,nembc=20
real, parameter    :: dtc=1.,ddnc=2.e-6 ,dthz=1.,drhhz=.02
!--------------------------------------------------------------------------
integer :: iccnflg,ifnflg,icloud,irain,ipris,isnow,iaggr,igraup,ihail  &
          ,mkcoltab
!integer, dimension(ncat)        :: jnmb
!integer, dimension(nhcat,nhcat) :: ipairc,ipairr
!integer, dimension(31,100,2)    :: jhabtab
!integer, dimension(nzpmax,ncat) :: jhcat,ict1,ict2

real :: cparm,rparm,pparm,sparm,aparm,gparm,hparm,rictmin,rictmax  &
       ,dps,dps2,d1min,r2min,d2min,d1max,r2max,d2max  &
       ,d1ecc,d1ecr,r2ecr,r2err,colf,pi4dt,sedtime0,sedtime1

real, dimension(nhcat) :: shape,cfmas,pwmas,cfvt,pwvt,dpsmi,cfden,pwden  &
                         ,cfemb0,cfen0,pwemb0,pwen0,vtfac,frefac1,frefac2  &
                         ,cfmasft,dnfac,sipfac,pwmasi,ch1,ch3,cdp1,pwvtmasi

character(len=80) :: coltabfn

common /micphys/cfmas,pwmas








