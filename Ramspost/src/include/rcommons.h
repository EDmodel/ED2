!############################# Change Log ##################################
! 4.3.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!  Mission Research Corporation / *ASTeR Division
!###########################################################################

!-------------------------------------------------------------------------------                                                                     *
!     Common block include file for the   R A M S   model
!------------------------------------------------------------------------------- 
!
include 'rconfig.h'
include 'micphys.h'

!-------------------------------------------------------------------------------
!        Parameters for some array dimensions

integer, parameter :: nxyzpm=nzpmax*nxpmax*nypmax  &
                     ,maxdimp=maxdim+1,nstyp=12,nvtyp=30  &
                     ,nkeep=90,nke=nkeep,nintgm=12,maxsndg=200  &
                     ,maxvarf=200,maxsstf=200

!-------------------------------------------------------------------------------
!        COMMON block includes

!-------------------------------------------------------------------------------
character(len=64) :: expnme
common /allch/ expnme
!-------------------------------------------------------------------------------
integer                     :: itopo,initial,impl,iadvl,iadvf,lonrad,ngrids  &
                              ,lsflg,ibnd,jbnd,icorflg,iexev,imassflx,ilwrtyp,iswrtyp,iref  &
                              ,jref,ihtran,nfpt,nsndg,ideltat,nacoust,iflag  &
                              ,ntopsmth,izflat,iyear1,imonth1,idate1,ihour1  &
                              ,itime1,isfcl,ihorgrad
integer, dimension(maxgrds) :: idiffk
common /all/ itopo,initial,impl,iadvl,iadvf,lonrad,ngrids  &
            ,lsflg,ibnd,jbnd,icorflg,ilwrtyp,iswrtyp,iref,jref  &
            ,ihtran,nfpt,nsndg,ideltat,nacoust,iflag,ntopsmth  &
            ,izflat,iyear1,imonth1,idate1,ihour1,itime1  &
            ,idiffk,isfcl,ihorgrad
!-------------------------------------------------------------------------------
integer                     :: naddsc,nestz1,nestz2,nzg,nzs,n_pft,iversion,npatch  &
                              ,nvegpat,nclouds
integer, dimension(maxgrds) :: nnqparm,nndtrat,nstratx,nstraty  &
                              ,ngbegun,nnacoust,nxtnest,nnsttop,nnstbot  &
                              ,nnxp,nnyp,nnzp,ninest,njnest,nknest
integer, dimension(nzpmax)  :: nstratz1,nstratz2
integer, dimension(nvtyp)   :: kroot
common /all/ nnqparm,nndtrat,nstratx,nstraty,nstratz1,nstratz2  &
            ,ngbegun,nnacoust,nxtnest,nnsttop,nnstbot  &
            ,nnxp,nnyp,nnzp,ninest,njnest,nknest  &
            ,naddsc,nestz1,nestz2,nzg,nzs,n_pft,iversion,npatch  &
            ,nvegpat,nclouds,kroot
!-------------------------------------------------------------------------------
integer, parameter                     :: maxsched=200,maxschent=5
integer                                :: nsubs
integer, dimension(maxsched,maxschent) :: isched
common/schedule/nsubs,isched
!-------------------------------------------------------------------------------
real                     :: brunt,wcldbs,drtcon,rmin,radfrq,distim,seatmp  &
                           ,ubmin,eps,albedo,dthcon,rmax  &
                           ,cphas,dtlong,topref,sspct,polelat,polelon
real, dimension(maxclouds):: confrq
real, dimension(maxgrds) :: platn,plonn,centlat,centlon  &
                           ,zkhkm,xkhkm,cflxy,cflz,csz,csx,akmin
integer                  :: nhemgrd2
common /all/ brunt,wcldbs,drtcon,rmin,radfrq,distim,seatmp  &
            ,confrq,cflxy,cflz,csz,csx,rmax,akmin  &
            ,ubmin,eps,albedo,xkhkm,zkhkm,dthcon  &
            ,centlat,centlon,cphas,dtlong,topref,sspct  &
            ,polelat,polelon,platn,plonn,nhemgrd2
!-------------------------------------------------------------------------------
integer                     :: nhemt,nhemu,nhemv  
integer, dimension(4,maxhp) :: ihem1tt,jhem1tt,ihem1uu,jhem1uu  &
                              ,ihem1uv,jhem1uv,ihem1vu,jhem1vu  &
                              ,ihem1vv,jhem1vv
integer, dimension(maxhp)   :: ihem2tt,jhem2tt,ihem2uu,jhem2uu  &
                              ,ihem2uv,jhem2uv,ihem2vu,jhem2vu  &
                              ,ihem2vv,jhem2vv
real, dimension(maxhp)      :: hlatt,hlatu,hlatv,hlont,hlonu,hlonv
real, dimension(4,maxhp)    :: whem1tt,whem1uu,whem1uv,whem1vu,whem1vv
common /hemisphere/ nhemt,nhemu,nhemv  &
                   ,ihem1tt,jhem1tt,whem1tt  &
                   ,ihem1uu,jhem1uu,whem1uu  &
                   ,ihem1uv,jhem1uv,whem1uv  &
                   ,ihem1vu,jhem1vu,whem1vu  &
                   ,ihem1vv,jhem1vv,whem1vv  &
                   ,ihem2tt,jhem2tt,ihem2uu,jhem2uu  &
                   ,ihem2uv,jhem2uv,ihem2vu,jhem2vu,ihem2vv,jhem2vv  &
                   ,hlatt,hlatu,hlatv,hlont,hlonu,hlonv
!-------------------------------------------------------------------------------
real, dimension(nzpmax,maxgrds) :: u01dn,v01dn,pi01dn,th01dn,dn01dn,rt01dn
common /reference1d/ u01dn,v01dn,pi01dn,th01dn,dn01dn,rt01dn
!-------------------------------------------------------------------------------
integer              :: ipsflg,itsflg,irtsflg,iusflg
real, dimension(maxsndg) :: us,vs,ts,thds,ps,hs,rts
common /init_sounding/ ipsflg,itsflg,irtsflg,iusflg,us,vs,ts,thds,ps,hs,rts
!-------------------------------------------------------------------------------
real, dimension(nstyp)        :: slden,slcpd,slbs,slcond  &
                                ,slcons,slmsts,slpots,ssand,sclay  &
                                ,sorgan,sporo,soilcp,slfc,emisg
real, dimension(nvtyp)        :: albedv,emisv,vglai,vgdlai,vgfrac,vgdfrac  &
                                ,vegzo,vgdisp 
real                          :: cmin,corg,cwat,cair,cka,ckw
real, dimension(nzgmax)       :: slz
real, dimension(nzgmax,nvtyp) :: root
common /soil_veg/ slden,slcpd,slbs,slcond,slcons,slz,slmsts  &
                 ,slpots,ssand,sclay,sorgan,sporo,soilcp,slfc,emisg  &
                 ,albedv,emisv,vglai,vgdlai,vgfrac,vgdfrac,vegzo,vgdisp  &
                 ,root,cmin,corg,cwat,cair,cka,ckw
!-------------------------------------------------------------------------------
real                     :: time,ztop,dzrat,dzmax
real, dimension(maxgrds) :: deltazn,deltaxn,deltayn,dtlongn,dimove,djmove  &
                           ,gridu,gridv
real, dimension(nzpmax)  :: zz
common /time_grid / time,deltazn,deltaxn,deltayn,ztop,zz,dzrat,dzmax  &
                   ,dtlongn,dimove,djmove,gridu,gridv
!-------------------------------------------------------------------------------
character(len=8), dimension(50) :: plfmt,pltit,iplfld
common /prtchr/plfmt,pltit,iplfld
!-------------------------------------------------------------------------------
integer                :: nplt
integer, dimension(50) :: ixsctn,iplvect,isbval,iaa,iab,joa,job,naavg,noavg
real                   :: frqprt
real, dimension(50)    :: plconlo,plconhi,plconin
common /prtcom/ nplt,ixsctn,iplvect,isbval,iaa,iab,joa,job,plconin  &
               ,naavg,noavg,plconlo,plconhi,frqprt
!-------------------------------------------------------------------------------
character (len=10) :: runtype
character(len=1)   :: timeunit
character (len=32) :: vtabcust
character(len=80)  :: hfilin,afilin,hfilout,afilout,sfcfiles
character(len=20)  :: xlite,ylite,zlite
common /filchr/ hfilin,hfilout,afilin,afilout  &
               ,runtype,timeunit,sfcfiles,vtabcust,xlite,ylite,zlite
!-------------------------------------------------------------------------------
integer :: ioutput,iinput,iopunt,kwrite,ihistdel,iclobber
real    :: frqhis,frqanl,timstr,avgtim,frqlite,frqmean,frqboth  
common /files/ frqhis,frqanl,ioutput,iinput,timstr,iopunt,kwrite  &
              ,ihistdel,avgtim,frqlite,frqmean,frqboth,iclobber
!-------------------------------------------------------------------------------

integer, dimension(maxgrds) :: itoptflg,isstflg,ivegtflg,isoilflg  &
                              ,nofilflg,itopsflg,iz0flg
real                        :: z0fact
real, dimension(maxgrds)    :: z0max,toptenh,toptwvl
common /topocom/ itoptflg,isstflg,ivegtflg,isoilflg,nofilflg  &
                ,toptenh,toptwvl,itopsflg,iz0flg,z0max,z0fact
!-------------------------------------------------------------------------------
character(len=80), dimension(maxgrds) :: itoptfn,isstfn,ivegtfn,isoilfn
common /topccom/ itoptfn,isstfn,ivegtfn,isoilfn
!-------------------------------------------------------------------------------
integer :: ngridsh
common /hisgrd/ ngridsh
!-------------------------------------------------------------------------------
integer :: level,nqparm,itopbc
common /option/ level,nqparm,itopbc
!-------------------------------------------------------------------------------
integer, dimension(maxgrds) :: nnx,nnx1,nnx2,nny,nny1,nny2,nnz,nnz1  &
                              ,nnxyzp,nnxysp,nnxyp
common /grpnts/ nnx,nnx1,nnx2,nny,nny1,nny2,nnz,nnz1  &
               ,nnxyzp,nnxysp,nnxyp
!-------------------------------------------------------------------------------
real, dimension(nzpmax,maxgrds) :: htn,ht2n,ht4n,hwn,hw2n,hw4n
common /sigmaz/ htn,ht2n,ht4n,hwn,hw2n,hw4n
!-------------------------------------------------------------------------------
real, dimension(nzpmax,maxgrds) :: dztn,dzmn,dzt2n,dzm2n,ztn,zmn
real, dimension(nxpmax,maxgrds) :: xtn,xmn
real, dimension(nypmax,maxgrds) :: ytn,ymn
common /spacing/ dztn,dzmn,dzt2n,dzm2n,xtn,xmn,ytn,ymn,ztn,zmn
!-------------------------------------------------------------------------------
real, dimension(nxpmax,6,6,2,maxgrds) :: advwtx
real, dimension(nypmax,6,6,2,maxgrds) :: advwty
real, dimension(nzpmax,6,6,2,maxgrds) :: advwtz
common /advctn/ advwtx,advwty,advwtz
!-------------------------------------------------------------------------------
integer                 :: nslcon,nvgcon
real                    :: zrough,pctlcon
real, dimension(nzgmax) :: stgoff,slmstr
common /soil_veg2/ stgoff,slmstr,nslcon,zrough,pctlcon,nvgcon
!-------------------------------------------------------------------------------
integer :: nxp,nx,nx1,nx2,nyp,ny,ny1,ny2,nzp,nzpp,nz,nz1  &
          ,nxyzp,nxyp,nxysp,nscl,nsttop,nstbot,ndtrat,jdim
common /dompts/ nxp,nx,nx1,nx2,nyp,ny,ny1,ny2,nzp,nzpp,nz,nz1  &
               ,nxyzp,nxyp,nxysp,nscl,nsttop,nstbot,ndtrat,jdim
!-------------------------------------------------------------------------------
real                    :: deltax,deltay,deltaz
real, dimension(nzpmax) :: ht,ht2,ht4,hw,hw2,hw4,zt,zm,dzt,dzm,dzt2,dzm2
real, dimension(nxpmax) :: xt,xm
real, dimension(nypmax) :: yt,ym
common /grid/ deltax,deltay,deltaz,ht,ht2,ht4,hw,hw2,hw4  &
             ,xt,xm,yt,ym,zt,zm,dzt,dzm,dzt2,dzm2
!-------------------------------------------------------------------------------
real, dimension(nzpmax,nypmax,4) :: cphx
real, dimension(nzpmax,nxpmax,4) :: cphy
common /bndcon/cphx,cphy
!-------------------------------------------------------------------------------

character(len=80)                     :: varfil1,varfil2,varfpfx
character(len=80), dimension(maxvarf) :: varfil
common /varchr/ varfil1,varfil2,varfil,varfpfx
!-------------------------------------------------------------------------------
real                    :: vtime1,vtime2,vwait1,vwaittot
real,dimension(maxvarf) :: vtime
integer                 :: nvarf
common /var2/ vtime1,vtime2,vtime,nvarf,vwait1,vwaittot
!-------------------------------------------------------------------------------
character(len=80)                             :: sstfpfx
character(len=80), dimension(maxgrds)         :: sstfil1,sstfil2
character(len=80), dimension(maxsstf,maxgrds) :: vsstfil,sstfil
common /sstchr/ vsstfil,sstfil,sstfil1,sstfil2,sstfpfx
!-------------------------------------------------------------------------------
integer                          :: iupdsst
integer, dimension(maxgrds)      :: nvsstf,nsstf,isstf,isstrecy
integer, dimension(maxsstf)      :: iyearvs,imonthvs,idatevs,ihourvs
real ,dimension(maxgrds)         :: ssttime1,ssttime2
real ,dimension(maxsstf,maxgrds) :: ssttime
common /sst2/ iyearvs,imonthvs,idatevs,ihourvs  &
             ,ssttime,ssttime1,ssttime2  &
             ,nvsstf,nsstf,isstf,iupdsst,isstrecy
!-------------------------------------------------------------------------------
real, dimension(maxdimp)    :: vctr1 ,vctr2 ,vctr11,vctr12
integer, dimension(maxdimp) :: ivctr 
common /vctemp/ vctr1 ,vctr2 ,vctr11 ,vctr12
!-------------------------------------------------------------------------------
integer :: isstp,istp,initfld
real    :: timmax,dts,dtlt,dtlv
common /cntrlr/ timmax,isstp,istp,dts,dtlt,dtlv,initfld
!-------------------------------------------------------------------------------
real, dimension(nzpmax) :: pi,p0,temp,prt,rc,thet,rvls
real                    :: pfactr,tnudlat,tnudcent,tnudtop,znudtop
integer                 :: nudlat
common /stuff/ pi,p0,temp,prt,rc,thet,rvls,pfactr  &
              ,nudlat,tnudlat,tnudcent,tnudtop,znudtop

!-------------------------------------------------------------------------------
!      I/O table commons and information

!-------------------------------------------------------------------------------
integer :: ngrid,ngridc,ngrido,iscr1,iscr2
common /ioinfo/ ngrid,ngridc,ngrido,iscr1,iscr2
!-------------------------------------------------------------------------------
integer                    :: memsize,iounit,maxpro,memscr,memind  &
                             ,iogrid,maxpts,maxnzp,maxnxp,maxnyp,i2dvar
integer,dimension(maxgrds) :: memgrd
common/ioparm/memsize,iounit,maxpro,memscr,memind  &
             ,iogrid,memgrd,maxpts,maxnzp,maxnxp,maxnyp,i2dvar

!-------------------------------------------------------------------------------
integer                     :: marker3  &
                              ,iup     ,iuc     ,ivp     ,ivc     ,iwp  &
                              ,iwc     ,ipp     ,ipc     ,ithp    ,irtp  &
                              ,ircp    ,irrp    ,irpp    ,irsp    ,irap  &
                              ,irgp    ,irhp    ,iccp    ,icrp    ,icpp  &
                              ,icsp    ,icap    ,icgp    ,ichp    ,icccnp  &
                              ,icifnp  ,ihkm    ,iq2     ,iq6     ,iq7  &
                              ,irv     ,itheta  ,itkep   ,itklp   ,ithsrc  &
                              ,irtsrc  ,ifthrd  ,ivarup  ,ivarvp  ,ivartp  &
                              ,ivarrp  ,ivaruf  ,ivarvf  ,ivartf  ,ivarrf  &
                              ,ivarwts ,ipi0    ,idn0    ,ith0    ,ivkm  &
                              ,ivkh    ,idn0u   ,idn0v &
                              ! add the averaged 3d variables starting at 
                              ! index 54 --> 86  
                              ,iupm     ,ivpm   ,iwpm    ,ippm    ,ircpm  &
                              ,irrpm    ,irppm  ,irspm   ,irapm   ,irgpm  &
                              ,irhpm    ,iccpm  ,icrpm   ,icppm   ,icspm  &
                              ,icapm    ,icgpm  ,ichpm   ,icccnpm ,icifnpm  &
                              ,ihkmm    ,iq2m   ,iq6m    ,iq7m    ,irvm  &
                              ,ithetam  ,itkepm ,itklpm  ,ithsrcm ,irtsrcm  &
                              ,ifthrdm  ,ivkmm  ,ivkhm
integer, dimension(maxsclr) :: isclp

! would normally add isclpm(maxsclr) here, but don't have an index
! to assign in vtable, so lets not allow averaging of added scalars

common/index3/ marker3  &
              ,iup     ,iuc     ,ivp     ,ivc     ,iwp  &
              ,iwc     ,ipp     ,ipc     ,ithp    ,irtp  &
              ,ircp    ,irrp    ,irpp    ,irsp    ,irap  &
              ,irgp    ,irhp    ,iccp    ,icrp    ,icpp  &
              ,icsp    ,icap    ,icgp    ,ichp    ,icccnp  &
              ,icifnp  ,ihkm    ,iq2     ,iq6     ,iq7  &
              ,irv     ,itheta  ,itkep   ,itklp   ,ithsrc  &
              ,irtsrc  ,ifthrd  ,ivarup  ,ivarvp  ,ivartp  &
              ,ivarrp  ,ivaruf  ,ivarvf  ,ivartf  ,ivarrf  &
              ,ivarwts ,ipi0    ,idn0    ,ith0    ,ivkm  &
              ,ivkh    ,idn0u   ,idn0v &
              ,iupm     ,ivpm   ,iwpm    ,ippm    ,ircpm  &
              ,irrpm    ,irppm  ,irspm   ,irapm   ,irgpm  &
              ,irhpm    ,iccpm  ,icrpm   ,icppm   ,icspm  &
              ,icapm    ,icgpm  ,ichpm   ,icccnpm ,icifnpm  &
              ,ihkmm    ,iq2m   ,iq6m    ,iq7m    ,irvm  &
              ,ithetam  ,itkepm ,itklpm  ,ithsrcm ,irtsrcm  &
              ,ifthrdm  ,ivkmm  ,ivkhm  &
              ,isclp

!-------------------------------------------------------------------------------
integer                     :: marker3m  &
                              ,iut     ,ivt     ,iwt     ,ipt     ,itht  &
                              ,irtt    ,irct    ,irrt    ,irpt    ,irst  &
                              ,irat    ,irgt    ,irht    ,icct    ,icrt  &
                              ,icpt    ,icst    ,icat    ,icgt    ,icht  &
                              ,icccnt  ,icifnt  ,idum1t  ,itket   ,itklt  &
                              ,ivt3da  ,ivt3db  ,ivt3dc  ,ivt3dd  ,ivt3de  &
                              ,ivt3df  ,ivt3dg  ,ivt3dh  ,ivt3di  ,ivt3dj  &
                              ,ivt3dk  ,ivt3dl  ,ivt3dm  ,ivt3dn  ,ivt3do  &
                              ,ivt3dp
integer, dimension(maxsclr) :: isclt
common/indmx3/ marker3m  &
              ,iut     ,ivt     ,iwt     ,ipt     ,itht  &
              ,irtt    ,irct    ,irrt    ,irpt    ,irst  &
              ,irat    ,irgt    ,irht    ,icct    ,icrt  &
              ,icpt    ,icst    ,icat    ,icgt    ,icht  &
              ,icccnt  ,icifnt  ,idum1t  ,itket   ,itklt  &
              ,ivt3da  ,ivt3db  ,ivt3dc  ,ivt3dd  ,ivt3de  &
              ,ivt3df  ,ivt3dg  ,ivt3dh  ,ivt3di  ,ivt3dj  &
              ,ivt3dk  ,ivt3dl  ,ivt3dm  ,ivt3dn  ,ivt3do  &
              ,ivt3dp  ,isclt
             
!-------------------------------------------------------------------------------
integer :: marker2  &
          ,itopt    ,itopu    ,itopv    ,itopm    ,irtgt  &
          ,irtgu    ,irtgv    ,irtgm    ,if13t    ,if13u  &
          ,if13v    ,if13m    ,if23t    ,if23u    ,if23v  &
          ,if23m    ,idxu     ,idxv     ,idxt     ,idxm  &
          ,idyu     ,idyv     ,idyt     ,idym     ,ifmapu  &
          ,ifmapv   ,ifmapt   ,ifmapm   ,ifmapui  ,ifmapvi  &
          ,ifmapti  ,ifmapmi  ,iglat    ,iglon    ,iuw  &
          ,ivw      ,iwfz     ,itfz     ,iqfz     ,iaccpr  &
          ,iaccpp   ,iaccps   ,iaccpa   ,iaccpg   ,iaccph  &
          ,ipcprr   ,ipcprp   ,ipcprs   ,ipcpra   ,ipcprg  &
          ,ipcprh   ,ipcpg    ,iqpcpg   ,idpcpg   ,iaconpr  &
          ,iconprr  ,irshort  ,irlong   ,irlongup ,ialbedt  &
          ,ivarp    ,ivarp2   ,ifcoru   ,ifcorv   ,ivt2da  &
          ,ivt2db   ,ivt2dc   ,ivt2dd   ,ivt2de   ,ivt2df  &
          ,icosz    ,itopzo  &
          ! add the averaged 2d variables starting at
          ! index 73 --> 103  
          ,itoptm   ,iglatm   ,iglonm   ,iuwm    ,ivwm  &
          ,iwfzm    ,itfzm    ,iqfzm    ,iaccprm ,iaccppm  &
          ,iaccpsm  ,iaccpam  ,iaccpgm  ,iaccphm ,ipcprrm  &
          ,ipcprpm  ,ipcprsm  ,ipcpram  ,ipcprgm ,ipcprhm  &
          ,ipcpgm   ,iqpcpgm  ,idpcpgm  ,iaconprm,iconprrm  &
          ,irshortm ,irlongm  ,irlongupm,ialbedtm,icoszm  &
          ,topzom
common/index2/ marker2  &
              ,itopt    ,itopu    ,itopv    ,itopm    ,irtgt  &
              ,irtgu    ,irtgv    ,irtgm    ,if13t    ,if13u  &
              ,if13v    ,if13m    ,if23t    ,if23u    ,if23v  &
              ,if23m    ,idxu     ,idxv     ,idxt     ,idxm  &
              ,idyu     ,idyv     ,idyt     ,idym     ,ifmapu  &
              ,ifmapv   ,ifmapt   ,ifmapm   ,ifmapui  ,ifmapvi  &
              ,ifmapti  ,ifmapmi  ,iglat    ,iglon    ,iuw  &
              ,ivw      ,iwfz     ,itfz     ,iqfz     ,iaccpr  &
              ,iaccpp   ,iaccps   ,iaccpa   ,iaccpg   ,iaccph  &
              ,ipcprr   ,ipcprp   ,ipcprs   ,ipcpra   ,ipcprg  &
              ,ipcprh   ,ipcpg    ,iqpcpg   ,idpcpg   ,iaconpr  &
              ,iconprr  ,irshort  ,irlong   ,irlongup ,ialbedt  &
              ,ivarp    ,ivarp2   ,ifcoru   ,ifcorv   ,ivt2da  &
              ,ivt2db   ,ivt2dc   ,ivt2dd   ,ivt2de   ,ivt2df  &
              ,icosz    ,itopzo  &
              ,itoptm   ,iglatm   ,iglonm   ,iuwm    ,ivwm  &
              ,iwfzm    ,itfzm    ,iqfzm    ,iaccprm ,iaccppm  &
              ,iaccpsm  ,iaccpam  ,iaccpgm  ,iaccphm ,ipcprrm  &
              ,ipcprpm  ,ipcprsm  ,ipcpram  ,ipcprgm ,ipcprhm  &
              ,ipcpgm   ,iqpcpgm  ,idpcpgm  ,iaconprm,iconprrm  &
              ,irshortm ,irlongm  ,irlongupm,ialbedtm,icoszm  &
              ,topzom

!-------------------------------------------------------------------------------
integer :: marker4s  &
          ,itgp     ,iwgp     ,ischar   ,igsf  &
          ! add the averaged 3d soil variables starting at
          ! index 5 --> 8  
          ,itgpm    ,iwgpm    ,ischarm  ,igsfm
common/index4s/ marker4s  &
               ,itgp    ,iwgp     ,ischar   ,igsf  &
               ,itgpm   ,iwgpm    ,ischarm  ,igsfm

!-------------------------------------------------------------------------------
integer :: iscc,iscp,isct
common/indtrc/ iscc,iscp,isct

save






