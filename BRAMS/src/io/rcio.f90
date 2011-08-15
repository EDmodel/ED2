!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

subroutine commio (cfile,io,iun)

   use mem_all
   use therm_lib, only: level
   use mem_mass, only: imassflx, iexev
   use leaf_coms, only : ubmin             & ! intent(inout)
                       , ugbmin            & ! intent(inout)
                       , ustmin            & ! intent(inout)
                       , leaf_gamm => gamm & ! intent(inout)
                       , leaf_gamh => gamh & ! intent(inout)
                       , tprandtl          & ! intent(inout)
                       , ribmax            & ! intent(inout)
                       , leaf_maxwhc       & ! intent(inout)
                       , min_patch_area    & ! intent(inout)
                       , kroot             & ! intent(inout)
                       , nvtyp             & ! intent(inout)
                       , nvtyp_teb         & ! intent(inout)
                       , nstyp             & ! intent(inout)
                       , slden             & ! intent(inout)
                       , slcpd             & ! intent(inout)
                       , slbs              & ! intent(inout)
                       , slcond            & ! intent(inout)
                       , slcons            & ! intent(inout)
                       , slmsts            & ! intent(inout)
                       , slpots            & ! intent(inout)
                       , ssand             & ! intent(inout)
                       , sclay             & ! intent(inout)
                       , sorgan            & ! intent(inout)
                       , sporo             & ! intent(inout)
                       , soilwp            & ! intent(inout)
                       , soilcp            & ! intent(inout)
                       , slfc              & ! intent(inout)
                       , emisg             & ! intent(inout)
                       , emisv             & ! intent(inout)
                       , root              & ! intent(inout)
                       , cmin              & ! intent(inout)
                       , corg              & ! intent(inout)
                       , cwat              & ! intent(inout)
                       , cair              & ! intent(inout)
                       , cka               & ! intent(inout)
                       , ckw               ! ! intent(inout)
   use grell_coms, only : iupmethod,radius,depth_min,cap_maxs,cld2prec        &
                         ,zkbmax,zcutdown,z_detr,max_heat,closure_type        &
                         ,maxens_lsf,maxens_eff,maxens_cap
   use turb_coms, only : nna, nnb, nnc
   implicit none
   integer :: iun
   character(len=*) :: io,cfile

   !     This routine reads or writes the history and analysis file common blocks.

   integer, external :: cio_i,cio_f,cio_f8,cio_c,cio_i_sca,cio_f_sca,cio_f8_sca,cio_c_sca
   character(len=2) :: cng
   integer :: irw,ie,ng,max_nnzp

   IF(IO.EQ.'READ') irw=1
   IF(IO.EQ.'WRITE') irw=2
   !      print*,'in commio',cfile,' ',io,' ',iun

   ie=cio_i_sca(iun,irw,'iversion',iversion,1)
   ie=cio_c_sca(iun,irw,'expnme',expnme,1)
   ie=cio_i_sca(iun,irw,'ioutput',ioutput,1)
   ie=cio_i_sca(iun,irw,'if_adap',if_adap,1)
   ie=cio_i_sca(iun,irw,'ngrids',ngrids,1)
   ie=cio_i_sca(iun,irw,'nzg',nzg,1)
   ie=cio_i_sca(iun,irw,'nzs',nzs,1)
   ie=cio_i_sca(iun,irw,'nclouds',nclouds,1)
   ie=cio_i_sca(iun,irw,'naddsc',naddsc,1)
   ie=cio_f8_sca(iun,irw,'time',time,1)
   ie=cio_f_sca(iun,irw,'ztop',ztop,1)
   ie=cio_f_sca(iun,irw,'polelat',polelat,1)
   ie=cio_f_sca(iun,irw,'polelon',polelon,1)
   ie=cio_f_sca(iun,irw,'dzrat',dzrat,1)
   ie=cio_f_sca(iun,irw,'dzmax',dzmax,1)

   ie=cio_i(iun,irw,'nnxp',nnxp,ngrids)
   ie=cio_i(iun,irw,'nnyp',nnyp,ngrids)
   ie=cio_i(iun,irw,'nnzp',nnzp,ngrids)
   max_nnzp=maxval(nnzp(1:ngrids))

   ie=cio_i(iun,irw,'nndtrat',nndtrat,ngrids)
   ie=cio_i(iun,irw,'nstratx',nstratx,ngrids)
   ie=cio_i(iun,irw,'nstraty',nstraty,ngrids)
   ie=cio_i(iun,irw,'ngbegun',ngbegun,ngrids)
   ie=cio_i(iun,irw,'nnacoust',nnacoust,ngrids)
   ie=cio_i(iun,irw,'nxtnest',nxtnest,ngrids)
   ie=cio_i(iun,irw,'nnsttop',nnsttop,ngrids)
   ie=cio_i(iun,irw,'nnstbot',nnstbot,ngrids)
   ie=cio_i(iun,irw,'ninest',ninest,ngrids)
   ie=cio_i(iun,irw,'njnest',njnest,ngrids)
   ie=cio_i(iun,irw,'nknest',nknest,ngrids)
   ie=cio_i(iun,irw,'idiffk',idiffk,ngrids)
   ie=cio_f(iun,irw,'gridu',gridu,ngrids)
   ie=cio_f(iun,irw,'gridv',gridv,ngrids)
   ie=cio_f(iun,irw,'nna',nna,2)
   ie=cio_f(iun,irw,'nnb',nnb,2)
   ie=cio_f(iun,irw,'nnc',nnc,5)
   ie=cio_f(iun,irw,'akmin',akmin,ngrids)
   ie=cio_f(iun,irw,'akmax',akmax,ngrids)
   ie=cio_f(iun,irw,'hgtmin',hgtmin,ngrids)
   ie=cio_f(iun,irw,'hgtmax',hgtmax,ngrids)
   ie=cio_f(iun,irw,'csz',csz,ngrids)
   ie=cio_f(iun,irw,'csx',csx,ngrids)
   ie=cio_f(iun,irw,'xkhkm',xkhkm,ngrids)
   ie=cio_f(iun,irw,'zkhkm',zkhkm,ngrids)
   ie=cio_f(iun,irw,'centlat',centlat,ngrids)
   ie=cio_f(iun,irw,'centlon',centlon,ngrids)
   ie=cio_f(iun,irw,'dimove',dimove,ngrids)
   ie=cio_f(iun,irw,'djmove',djmove,ngrids)
   ie=cio_f(iun,irw,'deltazn',deltazn,ngrids)
   ie=cio_f(iun,irw,'deltaxn',deltaxn,ngrids)
   ie=cio_f(iun,irw,'deltayn',deltayn,ngrids)
   ie=cio_f(iun,irw,'platn',platn,ngrids)
   ie=cio_f(iun,irw,'plonn',plonn,ngrids)
   ie=cio_f(iun,irw,'zz',zz,nnzp(1))

   ie=cio_i_sca(iun,irw,'nestz1',nestz1,1)
   ie=cio_i_sca(iun,irw,'nestz2',nestz2,1)
   ie=cio_i(iun,irw,'nstratz1',nstratz1,nnzp(1))
   ie=cio_i(iun,irw,'nstratz2',nstratz2,nnzp(1))

   do ng=1,ngrids
      write(cng,1) ng
   1       format(i2.2)
      ie=cio_f(iun,irw,'xmn'//cng,xmn(1,ng),nnxp(ng))
      ie=cio_f(iun,irw,'xtn'//cng,xtn(1,ng),nnxp(ng))
      ie=cio_f(iun,irw,'ymn'//cng,ymn(1,ng),nnyp(ng))
      ie=cio_f(iun,irw,'ytn'//cng,ytn(1,ng),nnyp(ng))
      ie=cio_f(iun,irw,'zmn'//cng,zmn(1,ng),nnzp(ng))
      ie=cio_f(iun,irw,'ztn'//cng,ztn(1,ng),nnzp(ng))
      ie=cio_f(iun,irw,'dzmn'//cng,dzmn(1,ng),nnzp(ng))
      ie=cio_f(iun,irw,'dztn'//cng,dztn(1,ng),nnzp(ng))
      ie=cio_f(iun,irw,'u01dn'//cng,u01dn(1,ng),nnzp(ng))
      ie=cio_f(iun,irw,'v01dn'//cng,v01dn(1,ng),nnzp(ng))
      ie=cio_f(iun,irw,'pi01dn'//cng,pi01dn(1,ng),nnzp(ng))
      ie=cio_f(iun,irw,'th01dn'//cng,th01dn(1,ng),nnzp(ng))
      ie=cio_f(iun,irw,'dn01dn'//cng,dn01dn(1,ng),nnzp(ng))
      ie=cio_f(iun,irw,'rt01dn'//cng,rt01dn(1,ng),nnzp(ng))
      ie=cio_f(iun,irw,'co201dn'//cng,co201dn(1,ng),nnzp(ng))
   enddo

   ie=cio_i(iun,irw,'kroot',kroot,nvtyp+nvtyp_teb)

   ie=cio_i    (iun,irw,'nnqparm' ,nnqparm,ngrids)
   ie=cio_i    (iun,irw,'ndeepest',ndeepest,ngrids)
   ie=cio_i    (iun,irw,'nshallowest',nshallowest,ngrids)
   ie=cio_f_sca(iun,irw,'wcldbs',wcldbs,1)
   ie=cio_f_sca(iun,irw,'confrq',confrq,1)
   ie=cio_i_sca(iun,irw,'iupmethod',iupmethod,1)
   ie=cio_f    (iun,irw,'radius',radius,nclouds)
   ie=cio_f    (iun,irw,'depth_min',depth_min,nclouds)
   ie=cio_f_sca(iun,irw,'cap_maxs',cap_maxs,1)
   ie=cio_f_sca(iun,irw,'cld2prec',cld2prec,1)
   ie=cio_f_sca(iun,irw,'zkbmax',zkbmax,1)
   ie=cio_f_sca(iun,irw,'zcutdown',zcutdown,1)
   ie=cio_f_sca(iun,irw,'z_detr',z_detr,1)
   ie=cio_f_sca(iun,irw,'max_heat',max_heat,1)
   ie=cio_c_sca(iun,irw,'closure_type',closure_type,1)
   ie=cio_i_sca(iun,irw,'maxens_lsf',maxens_lsf,1)
   ie=cio_i_sca(iun,irw,'maxens_eff',maxens_eff,1)
   ie=cio_i_sca(iun,irw,'maxens_cap',maxens_cap,1)

   ie=cio_i_sca(iun,irw,'itopo',itopo,1)
   ie=cio_i_sca(iun,irw,'initial',initial,1)
   ie=cio_i_sca(iun,irw,'impl',impl,1)
   ie=cio_i_sca(iun,irw,'iinput',iinput,1)
   ie=cio_i_sca(iun,irw,'jdim',jdim,1)
   ie=cio_i_sca(iun,irw,'iadvl',iadvl,1)
   ie=cio_i_sca(iun,irw,'iadvf',iadvf,1)
   ie=cio_i_sca(iun,irw,'lonrad',lonrad,1)
   ie=cio_i_sca(iun,irw,'lsflg',lsflg,1)
   ie=cio_i_sca(iun,irw,'ibnd',ibnd,1)
   ie=cio_i_sca(iun,irw,'jbnd',jbnd,1)
   ie=cio_i_sca(iun,irw,'icorflg',icorflg,1)
   ![MLO
   ie=cio_i_sca(iun,irw,'iexev',iexev,1)
   ie=cio_i_sca(iun,irw,'imassflx',imassflx,1)
   !MLO]
   ie=cio_i_sca(iun,irw,'ilwrtyp',ilwrtyp,1)
   ie=cio_i_sca(iun,irw,'iswrtyp',iswrtyp,1)
   ie=cio_i_sca(iun,irw,'icumfdbk',icumfdbk,1)
   ie=cio_i_sca(iun,irw,'iref',iref,1)
   ie=cio_i_sca(iun,irw,'jref',jref,1)
   ie=cio_i_sca(iun,irw,'ihtran',ihtran,1)
   ie=cio_i_sca(iun,irw,'nfpt',nfpt,1)
   ie=cio_i_sca(iun,irw,'nsndg',nsndg,1)
   ie=cio_i_sca(iun,irw,'ideltat',ideltat,1)
   ie=cio_i_sca(iun,irw,'nacoust',nacoust,1)
   ie=cio_i_sca(iun,irw,'iflag',iflag,1)
   ie=cio_i_sca(iun,irw,'ntopsmth',ntopsmth,1)
   ie=cio_i_sca(iun,irw,'izflat',izflat,1)
   ie=cio_i_sca(iun,irw,'iyear1',iyeara,1)
   ie=cio_i_sca(iun,irw,'imonth1',imontha,1)
   ie=cio_i_sca(iun,irw,'idate1',idatea,1)
   ie=cio_i_sca(iun,irw,'itime1',itimea,1)
   ie=cio_i_sca(iun,irw,'iyear2',iyearz,1)
   ie=cio_i_sca(iun,irw,'imonth2',imonthz,1)
   ie=cio_i_sca(iun,irw,'idate2',idatez,1)
   ie=cio_i_sca(iun,irw,'itime2',itimez,1)
   ie=cio_i_sca(iun,irw,'isfcl',isfcl,1)
   ie=cio_i_sca(iun,irw,'istar',istar,1)
   ie=cio_i_sca(iun,irw,'igrndvap',igrndvap,1)
   ie=cio_f_sca(iun,irw,'ubmin',ubmin,1)
   ie=cio_f_sca(iun,irw,'ugbmin',ugbmin,1)
   ie=cio_f_sca(iun,irw,'ustmin',ustmin,1)
   ie=cio_f_sca(iun,irw,'gamm',leaf_gamm,1)
   ie=cio_f_sca(iun,irw,'gamh',leaf_gamh,1)
   ie=cio_f_sca(iun,irw,'tprandtl',tprandtl,1)
   ie=cio_f_sca(iun,irw,'ribmax',ribmax,1)
   ie=cio_f_sca(iun,irw,'leaf_maxwhc',leaf_maxwhc,1)
   ie=cio_f_sca(iun,irw,'min_patch_area',min_patch_area,1)
   ie=cio_i_sca(iun,irw,'ico2',ico2,1)
   ie=cio_f    (iun,irw,'co2con',co2con,max_nnzp)
   ie=cio_i_sca(iun,irw,'npatch',npatch,1)
   ie=cio_i_sca(iun,irw,'nvegpat',nvegpat,1)
   ie=cio_i_sca(iun,irw,'level',level,1)
   ie=cio_i_sca(iun,irw,'irain',irain,1)
   ie=cio_i_sca(iun,irw,'ipris',ipris,1)
   ie=cio_i_sca(iun,irw,'isnow',isnow,1)
   ie=cio_i_sca(iun,irw,'iaggr',iaggr,1)
   ie=cio_i_sca(iun,irw,'igraup',igraup,1)
   ie=cio_i_sca(iun,irw,'icloud',icloud,1)
   ie=cio_i_sca(iun,irw,'ihail',ihail,1)

   ie=cio_i_sca(iun,irw,'ibruvais',ibruvais,1)
   ie=cio_i_sca(iun,irw,'ibotflx',ibotflx,1)
   ie=cio_f_sca(iun,irw,'brunt',brunt,1)
   ie=cio_f_sca(iun,irw,'drtcon',drtcon,1)
   ie=cio_f_sca(iun,irw,'rmin',rmin,1)
   ie=cio_f_sca(iun,irw,'radfrq',radfrq,1)
   ie=cio_f_sca(iun,irw,'distim',distim,1)
   ie=cio_f_sca(iun,irw,'seatmp',seatmp,1)
   ie=cio_f_sca(iun,irw,'rmax',rmax,1)
   ie=cio_f_sca(iun,irw,'eps',eps,1)
   ie=cio_f_sca(iun,irw,'albedo',albedo,1)
   ie=cio_f_sca(iun,irw,'dthcon',dthcon,1)
   ie=cio_f_sca(iun,irw,'cphas',cphas,1)
   ie=cio_f_sca(iun,irw,'topref',topref,1)
   ie=cio_f_sca(iun,irw,'sspct',sspct,1)
   ie=cio_f_sca(iun,irw,'rparm',rparm,1)
   ie=cio_f_sca(iun,irw,'pparm',pparm,1)
   ie=cio_f_sca(iun,irw,'sparm',sparm,1)
   ie=cio_f_sca(iun,irw,'aparm',aparm,1)
   ie=cio_f_sca(iun,irw,'gparm',gparm,1)
   ie=cio_f_sca(iun,irw,'cparm',cparm,1)
   ie=cio_f_sca(iun,irw,'hparm',hparm,1)

   ie=cio_f(iun,irw,'cfmas',cfmas,nhcat)
   ie=cio_f(iun,irw,'pwmas',pwmas,nhcat)

   ie=cio_f(iun,irw,'us',us,maxsndg)
   ie=cio_f(iun,irw,'vs',vs,maxsndg)
   ie=cio_f(iun,irw,'ts',ts,maxsndg)
   ie=cio_f(iun,irw,'thds',thds,maxsndg)
   ie=cio_f(iun,irw,'ps',ps,maxsndg)
   ie=cio_f(iun,irw,'hs',hs,maxsndg)
   ie=cio_f(iun,irw,'co2s',co2s,maxsndg)

   ie=cio_f(iun,irw,'slden',slden,nstyp)
   ie=cio_f(iun,irw,'slcpd',slcpd,nstyp)
   ie=cio_f(iun,irw,'slbs',slbs,nstyp)
   ie=cio_f(iun,irw,'slcond',slcond,nstyp)
   ie=cio_f(iun,irw,'slcons',slcons,nstyp)
   ie=cio_f(iun,irw,'slmsts',slmsts,nstyp)
   ie=cio_f(iun,irw,'slpots',slpots,nstyp)
   ie=cio_f(iun,irw,'ssand',ssand,nstyp)
   ie=cio_f(iun,irw,'sclay',sclay,nstyp)
   ie=cio_f(iun,irw,'sorgan',sorgan,nstyp)
   ie=cio_f(iun,irw,'sporo',sporo,nstyp)
   ie=cio_f(iun,irw,'soilwp',soilwp,nstyp)
   ie=cio_f(iun,irw,'soilcp',soilcp,nstyp)
   ie=cio_f(iun,irw,'slfc',slfc,nstyp)
   ie=cio_f(iun,irw,'emisg',emisg,nstyp)

   ie=cio_f(iun,irw,'emisv',emisv,nvtyp+nvtyp_teb)

   ie=cio_f(iun,irw,'root',root,nzgmax*(nvtyp+nvtyp_teb))
   ie=cio_f(iun,irw,'slz',slz,nzg)

   ie=cio_f_sca(iun,irw,'cmin',cmin,1)
   ie=cio_f_sca(iun,irw,'corg',corg,1)
   ie=cio_f_sca(iun,irw,'cwat',cwat,1)
   ie=cio_f_sca(iun,irw,'cair',cair,1)
   ie=cio_f_sca(iun,irw,'cka',cka,1)
   ie=cio_f_sca(iun,irw,'ckw',ckw,1)

   return
end subroutine commio 

!---------------------------------------------------------

subroutine cio_pos_file(iun,cstr,ierr)
implicit none
integer :: iun,ierr
character(len=*) :: cstr
character(len=256) :: line,csearch

integer :: nl,nc,iend,ilen

!      print*,'cio_pos:',iun,cstr

iend=0
1    continue
do nl=1,1000000
   read(iun,10,end=100) line
10      format(a)
   ilen=len(cstr)
   csearch='__'//cstr(1:ilen)
   nc=index(line,csearch(1:ilen+2) )
!         print*,'cio_pos:',nl,nc,line
   if(nc.eq.1) then
      ierr=0
!            print*,'---- Name found on header file:',cstr
      return
   endif
enddo

100  continue
if(iend.eq.1) then
   ierr=1
   print*,'---- Name NOT found on header file:',cstr
   rewind(iun)
   return
endif
rewind(iun)
iend=1
goto 1

end

!---------------------------------------------------------

integer function cio_i(iun,irw,cstr,ia,n)
implicit none
integer :: iun,irw,n
integer ia(*)
character(len=*) :: cstr
character(len=256) :: string
integer :: nn,i

if (irw.eq.1) then
   call cio_pos_file (iun,cstr,cio_i)
   if(cio_i.eq.1) return
   read(iun,*) nn
   read(iun,*) (ia(i),i=1,nn)
elseif(irw.eq.2) then
   write(iun,20) cstr
20      format('__',a)
   write(iun,*) n
   write(iun,11) (ia(i),i=1,n)
11      format(i6)
   cio_i=0
endif

return
end function cio_i

!---------------------------------------------------------

integer function cio_f(iun,irw,cstr,ia,n)
implicit none
integer :: iun,irw,n
real ia(*)
character(len=*) :: cstr
character(len=256) :: string
integer :: nn,i

if (irw.eq.1) then
   call cio_pos_file (iun,cstr,cio_f)
   if(cio_f.eq.1) return
   read(iun,*) nn
   read(iun,*) (ia(i),i=1,nn)
elseif(irw.eq.2) then
   write(iun,20) cstr
20      format('__',a)
   write(iun,*) n
   write(iun,11) (ia(i),i=1,n)
11      format(e16.8)
   cio_f=0
endif

return
end function cio_f

!---------------------------------------------------------

integer function cio_f8(iun,irw,cstr,ia,n)
implicit none
integer :: iun,irw,n
real(kind=8) :: ia(*)
character(len=*) :: cstr
character(len=256) :: string
integer :: nn,i

if (irw.eq.1) then
   call cio_pos_file (iun,cstr,cio_f8)
   if(cio_f8.eq.1) return
   read(iun,*) nn
   read(iun,*) (ia(i),i=1,nn)
elseif(irw.eq.2) then
   write(iun,20) cstr
20      format('__',a)
   write(iun,*) n
   write(iun,11) (ia(i),i=1,n)
11      format(e24.16)
   cio_f8=0
endif

return
end function cio_f8

!---------------------------------------------------------

integer function cio_c(iun,irw,cstr,ia,n)
  implicit none
  integer :: iun,irw,n
  character(len=*) :: ia(*)
  character(len=*) :: cstr
  character(len=256) :: string
  integer :: nn,i

  if (irw.eq.1) then
     call cio_pos_file (iun,cstr,cio_c)
     if(cio_c.eq.1) return
     read(iun,*) nn
     read(iun,10) (ia(i),i=1,nn)
  elseif(irw.eq.2) then
     write(iun,20) cstr
  20      format('__',a)
     write(iun,*) n
     write(iun,10) (ia(i),i=1,n)
  10      format(a)
     cio_c=0
  endif

  return
end function cio_c




integer function cio_i_sca(iun,irw,cstr,ia,n)
implicit none
integer :: iun,irw,n
integer ia
character(len=*) :: cstr
character(len=256) :: string
integer :: nn,i

if (n /= 1) then
  write (*,*) 'You should not have scalar call for this variable: '//cstr
  cio_i_sca=1
  stop 'cio_i_sca'
end if

if (irw.eq.1) then
   call cio_pos_file (iun,cstr,cio_i_sca)
   if(cio_i_sca.eq.1) return
   read(iun,*) nn
   if (nn /= 1) then
     write (*,*) 'The variable '//trim(cstr)//' should be scalar but I found a vector instead!!!'
     cio_i_sca=1
     stop 'cio_i_sca'
   end if
   read(iun,*) ia
elseif(irw.eq.2) then
   write(iun,20) cstr
20      format('__',a)
   write(iun,*) n
   write(iun,11) ia
11      format(i6)
   cio_i_sca=0
endif

return
end function cio_i_sca

!---------------------------------------------------------

integer function cio_f_sca(iun,irw,cstr,ia,n)
implicit none
integer :: iun,irw,n
real ia
character(len=*) :: cstr
character(len=256) :: string
integer :: nn,i

if (n /= 1) then
  write (*,*) 'You should not have scalar call for this variable: '//cstr
  cio_f_sca=1
  stop 'cio_f_sca'
end if

if (irw.eq.1) then
   call cio_pos_file (iun,cstr,cio_f_sca)
   if(cio_f_sca.eq.1) return
   read(iun,*) nn
   if (nn /= 1) then
     write (*,*) 'The variable '//trim(cstr)//' should be scalar but I found a vector instead!!!'
     cio_f_sca=1
     stop 'cio_f_sca'
   end if
   read(iun,*) ia
elseif(irw.eq.2) then
   write(iun,20) cstr
20      format('__',a)
   write(iun,*) n
   write(iun,11) ia
11      format(e16.8)
   cio_f_sca=0
endif

return
end function cio_f_sca

!---------------------------------------------------------

integer function cio_f8_sca(iun,irw,cstr,ia,n)
implicit none
integer :: iun,irw,n
real(kind=8) :: ia
character(len=*) :: cstr
character(len=256) :: string
integer :: nn,i

if (n /= 1) then
  write (*,*) 'You should not have scalar call for this variable: '//trim(cstr)
  cio_f8_sca=1
  stop 'cio_f8_sca'
end if

if (irw.eq.1) then
   call cio_pos_file (iun,cstr,cio_f8_sca)
   if(cio_f8_sca.eq.1) return
   read(iun,*) nn
   if (nn /= 1) then
     write (*,*) 'The variable '//trim(cstr)//' should be scalar but I found a vector instead!!!'
     cio_f8_sca=1
     stop 'cio_f8_sca'
   end if
   read(iun,*) ia
elseif(irw.eq.2) then
   write(iun,20) cstr
20      format('__',a)
   write(iun,*) n
   write(iun,11) ia
11      format(e24.16)
   cio_f8_sca=0
endif

return
end function cio_f8_sca

!---------------------------------------------------------

integer function cio_c_sca(iun,irw,cstr,ia,n)
implicit none
integer :: iun,irw,n
character(len=*) :: ia
character(len=*) :: cstr
character(len=256) :: string
integer :: nn,i

if (n /= 1) then
  write (*,*) 'You should not have scalar call for this variable: '//trim(cstr)
  cio_c_sca=1
  stop 'cio_c_sca'
end if

if (irw.eq.1) then
   call cio_pos_file (iun,cstr,cio_c_sca)
   if(cio_c_sca.eq.1) return
   read(iun,*) nn
   if (nn /= 1) then
     write (*,*) 'The variable '//trim(cstr)//' should be scalar but I found a vector instead!!!'
     cio_c_sca=1
     stop 'cio_c_sca'
   end if
   read(iun,10) ia
elseif(irw.eq.2) then
   write(iun,20) cstr
20      format('__',a)
   write(iun,*) n
   write(iun,10) ia
10      format(a)
   cio_c_sca=0
endif

return
end function cio_c_sca


