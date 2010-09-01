!==========================================================================================!
!==========================================================================================!
!    Ramspost main memory module.  This replaces rcommons.h, so we can avoid common blocks !
! that are very easy to mess up configurations.                                            !
!------------------------------------------------------------------------------------------!
module rpost_coms
   use rpost_dims

   !---------------------------------------------------------------------------------------!
   !     A long list of variables...                                                       !
   !---------------------------------------------------------------------------------------!
   character(len=64)                      :: expnme

   integer                                :: itopo
   integer                                :: initial
   integer                                :: impl
   integer                                :: iadvl
   integer                                :: iadvf
   integer                                :: lonrad
   integer                                :: ngrids
   integer                                :: lsflg
   integer                                :: ibnd
   integer                                :: jbnd
   integer                                :: icorflg
   integer                                :: iexev
   integer                                :: imassflx
   integer                                :: ilwrtyp
   integer                                :: iswrtyp
   integer                                :: iref
   integer                                :: jref
   integer                                :: ihtran
   integer                                :: nfpt
   integer                                :: nsndg
   integer                                :: ideltat
   integer                                :: nacoust
   integer                                :: iflag
   integer                                :: ntopsmth
   integer                                :: izflat
   integer                                :: iyear1
   integer                                :: imonth1
   integer                                :: idate1
   integer                                :: ihour1
   integer                                :: itime1
   integer                                :: isfcl
   integer                                :: istar
   integer                                :: ihorgrad
   integer, dimension(maxgrds)            :: idiffk
   integer                                :: naddsc
   integer                                :: nestz1
   integer                                :: nestz2
   integer                                :: nzg
   integer                                :: nzs
   integer                                :: n_pft
   integer                                :: iversion
   integer                                :: npatch
   integer                                :: nvegpat
   integer                                :: nclouds
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   integer, dimension(maxgrds)            :: nnqparm
   integer, dimension(maxgrds)            :: nndtrat
   integer, dimension(maxgrds)            :: nstratx
   integer, dimension(maxgrds)            :: nstraty
   integer, dimension(maxgrds)            :: ngbegun
   integer, dimension(maxgrds)            :: nnacoust
   integer, dimension(maxgrds)            :: nxtnest
   integer, dimension(maxgrds)            :: nnsttop
   integer, dimension(maxgrds)            :: nnstbot
   integer, dimension(maxgrds)            :: nnxp
   integer, dimension(maxgrds)            :: nnyp
   integer, dimension(maxgrds)            :: nnzp
   integer, dimension(maxgrds)            :: ninest
   integer, dimension(maxgrds)            :: njnest
   integer, dimension(maxgrds)            :: nknest
   integer, dimension(nzpmax)             :: nstratz1
   integer, dimension(nzpmax)             :: nstratz2
   integer, dimension(nvtyp)              :: kroot
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   integer, parameter                     :: maxsched=200,maxschent=5
   integer                                :: nsubs
   integer, dimension(maxsched,maxschent) :: isched
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   real                                   :: brunt
   real                                   :: wcldbs
   real                                   :: drtcon
   real                                   :: rmin
   real                                   :: radfrq
   real                                   :: distim
   real                                   :: seatmp
   real                                   :: eps
   real                                   :: albedo
   real                                   :: dthcon
   real                                   :: rmax
   real                                   :: cphas
   real                                   :: dtlong
   real                                   :: topref
   real                                   :: sspct
   real                                   :: polelat
   real                                   :: polelon
   real, dimension(maxclouds)             :: confrq
   real, dimension(maxgrds)               :: platn
   real, dimension(maxgrds)               :: plonn
   real, dimension(maxgrds)               :: centlat
   real, dimension(maxgrds)               :: centlon
   real, dimension(maxgrds)               :: zkhkm
   real, dimension(maxgrds)               :: xkhkm
   real, dimension(maxgrds)               :: cflxy
   real, dimension(maxgrds)               :: cflz
   real, dimension(maxgrds)               :: csz
   real, dimension(maxgrds)               :: csx
   real, dimension(maxgrds)               :: akmin
   integer                                :: nhemgrd2
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   integer                     :: nhemt,nhemu,nhemv  
   integer, dimension(4,maxhp) :: ihem1tt,jhem1tt,ihem1uu,jhem1uu  &
                                 ,ihem1uv,jhem1uv,ihem1vu,jhem1vu  &
                                 ,ihem1vv,jhem1vv
   integer, dimension(maxhp)   :: ihem2tt,jhem2tt,ihem2uu,jhem2uu  &
                                 ,ihem2uv,jhem2uv,ihem2vu,jhem2vu  &
                                 ,ihem2vv,jhem2vv
   real, dimension(maxhp)      :: hlatt,hlatu,hlatv,hlont,hlonu,hlonv
   real, dimension(4,maxhp)    :: whem1tt,whem1uu,whem1uv,whem1vu,whem1vv
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   real, dimension(nzpmax,maxgrds) :: u01dn,v01dn,pi01dn,th01dn,dn01dn,rt01dn
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   integer              :: ipsflg,itsflg,irtsflg,iusflg
   real, dimension(maxsndg) :: us,vs,ts,thds,ps,hs,rts
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   real, dimension(nstyp)        :: slden,slcpd,slbs,slcond  &
                                   ,slcons,slmsts,slpots,ssand,sclay  &
                                   ,sorgan,sporo,soilcp,slfc,emisg
   real, dimension(nvtyp)        :: albedv,emisv,vglai,vgdlai,vgfrac,vgdfrac  &
                                   ,vegzo,vgdisp 
   real                          :: cmin,corg,cwat,cair,cka,ckw
   real, dimension(nzgmax)       :: slz
   real, dimension(nzgmax,nvtyp) :: root
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   real                     :: time,ztop,dzrat,dzmax
   real, dimension(maxgrds) :: deltazn,deltaxn,deltayn,dtlongn,dimove,djmove  &
                              ,gridu,gridv
   real, dimension(nzpmax)  :: zz
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   character(len=8), dimension(50) :: plfmt,pltit,iplfld
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   integer                :: nplt
   integer, dimension(50) :: ixsctn,iplvect,isbval,iaa,iab,joa,job,naavg,noavg
   real                   :: frqprt
   real   , dimension(50) :: plconlo,plconhi,plconin
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   character (len=10)     :: runtype
   character(len=1)       :: timeunit
   character (len=32)     :: vtabcust
   character(len=80)      :: hfilin,afilin,hfilout,afilout,sfcfiles
   character(len=20)      :: xlite,ylite,zlite
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   integer :: ioutput,iinput,iopunt,kwrite,ihistdel,iclobber
   real    :: frqhis,frqanl,timstr,avgtim,frqlite,frqmean,frqboth  
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   integer, dimension(maxgrds) :: itoptflg,isstflg,ivegtflg,isoilflg  &
                                 ,nofilflg,itopsflg,iz0flg
   real                        :: z0fact
   real, dimension(maxgrds)    :: z0max,toptenh,toptwvl
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   character(len=80), dimension(maxgrds) :: itoptfn,isstfn,ivegtfn,isoilfn
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   integer :: ngridsh
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   integer :: level,nqparm,itopbc
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   integer, dimension(maxgrds) :: nnx,nnx1,nnx2,nny,nny1,nny2,nnz,nnz1  &
                                 ,nnxyzp,nnxysp,nnxyp
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   real, dimension(nzpmax,maxgrds) :: htn,ht2n,ht4n,hwn,hw2n,hw4n
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   real, dimension(nzpmax,maxgrds) :: dztn,dzmn,dzt2n,dzm2n,ztn,zmn
   real, dimension(nxpmax,maxgrds) :: xtn,xmn
   real, dimension(nypmax,maxgrds) :: ytn,ymn
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   real, dimension(nxpmax,6,6,2,maxgrds) :: advwtx
   real, dimension(nypmax,6,6,2,maxgrds) :: advwty
   real, dimension(nzpmax,6,6,2,maxgrds) :: advwtz
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   integer                 :: nslcon,nvgcon
   real                    :: zrough,pctlcon
   real, dimension(nzgmax) :: stgoff,slmstr
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   integer :: nxp,nx1,nx2,nyp,ny1,ny2,nzp,nzpp,nz1  &
             ,nxyzp,nxyp,nxysp,nscl,nsttop,nstbot,ndtrat,jdim
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   real                    :: deltax,deltay,deltaz
   real, dimension(nzpmax) :: ht,ht2,ht4,hw,hw2,hw4,zt,zm,dzt,dzm,dzt2,dzm2
   real, dimension(nxpmax) :: xt,xm
   real, dimension(nypmax) :: yt,ym
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   real, dimension(nzpmax,nypmax,4) :: cphx
   real, dimension(nzpmax,nxpmax,4) :: cphy
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   character(len=80)                     :: varfil1,varfil2,varfpfx
   character(len=80), dimension(maxvarf) :: varfil
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   real                    :: vtime1,vtime2,vwait1,vwaittot
   real,dimension(maxvarf) :: vtime
   integer                 :: nvarf
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   character(len=80)                             :: sstfpfx
   character(len=80), dimension(maxgrds)         :: sstfil1,sstfil2
   character(len=80), dimension(maxsstf,maxgrds) :: vsstfil,sstfil
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   integer                          :: iupdsst
   integer, dimension(maxgrds)      :: nvsstf,nsstf,isstf,isstrecy
   integer, dimension(maxsstf)      :: iyearvs,imonthvs,idatevs,ihourvs
   real ,dimension(maxgrds)         :: ssttime1,ssttime2
   real ,dimension(maxsstf,maxgrds) :: ssttime
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   real   , dimension(maxdimp) :: vctr1 ,vctr2 ,vctr11,vctr12
   integer, dimension(maxdimp) :: ivctr 
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   integer :: isstp,istp,initfld
   real    :: timmax,dts,dtlt,dtlv
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   real, dimension(nzpmax) :: pi,p0,temp,prt,rc,thet,rvls
   real                    :: pfactr,tnudlat,tnudcent,tnudtop,znudtop
   integer                 :: nudlat
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   integer :: ngrid,ngridc,ngrido,iscr1,iscr2
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   integer                    :: memsize,iounit,maxpro,memscr,memind  &
                                ,iogrid,maxpts,maxnzp,maxnxp,maxnyp,i2dvar
   integer,dimension(maxgrds) :: memgrd
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
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
                                ,iupm     ,ivpm   ,iwpm    ,ippm    ,ircpm  &
                                ,irrpm    ,irppm  ,irspm   ,irapm   ,irgpm  &
                                ,irhpm    ,iccpm  ,icrpm   ,icppm   ,icspm  &
                                ,icapm    ,icgpm  ,ichpm   ,icccnpm ,icifnpm  &
                                ,ihkmm    ,iq2m   ,iq6m    ,iq7m    ,irvm  &
                                ,ithetam  ,itkepm ,itklpm  ,ithsrcm ,irtsrcm  &
                                ,ifthrdm  ,ivkmm  ,ivkhm
   integer, dimension(maxsclr) :: isclp
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
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
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
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
             ,itoptm   ,iglatm   ,iglonm   ,iuwm    ,ivwm  &
             ,iwfzm    ,itfzm    ,iqfzm    ,iaccprm ,iaccppm  &
             ,iaccpsm  ,iaccpam  ,iaccpgm  ,iaccphm ,ipcprrm  &
             ,ipcprpm  ,ipcprsm  ,ipcpram  ,ipcprgm ,ipcprhm  &
             ,ipcpgm   ,iqpcpgm  ,idpcpgm  ,iaconprm,iconprrm  &
             ,irshortm ,irlongm  ,irlongupm,ialbedtm,icoszm  &
             ,topzom
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   integer :: marker4s  &
            ,itgp     ,iwgp     ,ischar   ,igsf  &
            ,itgpm    ,iwgpm    ,ischarm  ,igsfm
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   integer :: iscc,iscp,isct
   !---------------------------------------------------------------------------------------!
end module rpost_coms
!==========================================================================================!
!==========================================================================================!



