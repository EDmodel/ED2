!*************************************************************
subroutine RAMS_anal_init(nfile,fnames,file_prefix,          &
      dep_zlev,iep_nx,iep_ny,iep_nz,iep_ng,iep_np,iep_nc,    &
      iep_stdate,iep_step,iep_ngrids)

use an_header
include 'rcommons.h'

parameter (maxfiles=2000)

!SRF
dimension dep_zlev(nzpmax,maxgrds),iep_nx(maxgrds),iep_ny(maxgrds),  &
      iep_nz(maxgrds),iep_stdate(6),iep_step(6)
!SRF

character*240 fnames(maxfiles),fnamesp(maxfiles)  &
     ,file_prefix*(*),fpref*256
common/rams_data/ ftimes(maxfiles),nfgpnts(4,maxgrds,maxfiles)  &
     ,nfgrids(maxfiles),ifdates(maxfiles),iftimes(maxfiles)  &
     ,flevels(nzpmax,maxgrds,maxfiles),startutc,httop  &
     ,fdelx(maxgrds,maxfiles),fdely(maxgrds,maxfiles)
common /mem/ memsize4

dimension fdata(*),idata(*)
character*(*) cdata(*)
dimension mondays(13)
data mondays/31,28,31,30,31,30,31,31,30,31,30,31,31/



maxmem=0
fpref=file_prefix
nc=lastchar(fpref)+1
nfile=-1
fpref=file_prefix
fpref(nc:)='*-head.txt'
print*,'RAMS_filelist searching for: ',fpref(1:nc+10)

 call RAMS_filelist(fnames,fpref,nfile)

! construct arrays of various stuff

do nfn=1,nfile
   open(10,file=fnames(nfn),form='formatted')

   read(10,*) nvbtab

print*,fnames(nfn)(1:lastchar(fnames(nfn))),nvbtab

   allocate (anal_table(nvbtab))
   do nv=1,nvbtab
      read(10,*) anal_table(nv)%string   &
                ,anal_table(nv)%npointer  &
                ,anal_table(nv)%idim_type  &
                ,anal_table(nv)%ngrid  &
                ,anal_table(nv)%nvalues


!      print*, anal_table(nv)%string   &
!                ,anal_table(nv)%npointer  &
!                ,anal_table(nv)%idim_type  &
!                ,anal_table(nv)%ngrid  &
!                ,anal_table(nv)%nvalues


   enddo
   call commio('ANAL','READ',10)
   close(10)
   ftimes(nfn)=time
   istrhrs=nint(float(itime1)/100.+0.0001)
   strtim=istrhrs+float(itime1-istrhrs*100)/60.
   startutc=strtim
!   print*,'X',nfn,ftimes(nfn)
   deallocate (anal_table)
enddo

call RAMS_fltsort(nfile,ftimes,fnames)

do nfn=1,nfile
   open(10,file=fnames(nfn),form='formatted')
   read(10,*) nvbtab
   allocate (anal_table(nvbtab))
   do nv=1,nvbtab
      read(10,*) anal_table(nv)%string   &
                ,anal_table(nv)%npointer  &
                ,anal_table(nv)%idim_type  &
                ,anal_table(nv)%ngrid  &
                ,anal_table(nv)%nvalues
   enddo
   call commio('ANAL','READ',10)
   close(10)
   ftimes(nfn)=time
 if(nfn.eq.1) then
   iep_ngrids=ngrids
   do n=1,ngrids
    maxmem=max(maxmem,nnxp(n)*nnyp(n)*max(nnzp(n),npatch*nzg,nnzp(n)*nclouds))
    iep_nx(n)=nnxp(n)
    iep_ny(n)=nnyp(n)
    iep_nz(n)=nnzp(n)
    iep_ng   =nzg
    iep_np   =npatch
    iep_nc   =nclouds
    do nn=1,nnzp(n)
      dep_zlev(nn,n)=ztn(nn,n)
    enddo
   enddo
endif
   itme=nint(strtim*3600+time)
   iadddays=itme/86400
   izhours =(itme-iadddays*86400)/3600
   izmin   =(itme-iadddays*86400-izhours*3600)/60
   izsec   =itme-iadddays*86400-izhours*3600-izmin*60

   iftimes(nfn)= izhours
   iftimes(nfn)= iftimes(nfn) *100 + izmin
   iftimes(nfn)= iftimes(nfn) *100 + izsec
   if(iyear1.gt.50.and.iyear1.lt.100) iyear1=iyear1+1900
   if(iyear1.lt.50) iyear1=iyear1+2000
   mondays(2)=28+(1-min(1,mod(iyear1,4)))

   iiyear=iyear1
   iidate=idate1+iadddays
   iimon=imonth1
   101 if(iidate.gt.mondays(iimon)) then
   iidate=iidate-mondays(iimon)
   iimon=iimon+1
   goto 101
   endif
      102 if(iimon.gt.12) then
      iiyear=iiyear+1
      iimon=iimon-12
      go to 102
   endif
   
   ifdates(nfn)=iiyear       *100  + iimon
   ifdates(nfn)=ifdates(nfn) *100  + iidate

   nfgrids(nfn)=ngrids

   do ng=1,ngrids
      nfgpnts(1,ng,nfn)=nnxp(ng)
      nfgpnts(2,ng,nfn)=nnyp(ng)
      nfgpnts(3,ng,nfn)=nnzp(ng)
      nfgpnts(4,ng,nfn)=nzg
      fdelx(ng,nfn)=DELTAXN(NG)
      fdely(ng,nfn)=DELTAYN(NG)
         !print*,ng,nnxp(ng),nnyp(ng),nnzp(ng),nzg

      do k=1,nnzp(ng)
         flevels(k,ng,nfn)=ztn(k,ng)
     enddo
   enddo
   httop=zmn(nnzp(1)-1,1)

   close(10)
!!!!!!!!!!
!SRF : nao desaloque no final, sera usado em outra rotina
   if(nfn.lt.nfile) deallocate (anal_table)
!!!!!!!!!!
enddo

memsize4=maxmem

return
!srf

entry RAMS_get_time_init(nfl,iyear,imonth,idate,ihour,imin)
 iyear =iyear1
 imonth=imonth1
 idate =idate1
 ihour =int(float(iftimes(nfl))/10000.)
 imin  =int(float(iftimes(nfl)-10000*ihour)/100.)
return

entry RAMS_get_time_step(iistep,hunit,nfiles)
 if(nfiles.eq.1) then
  iistep = 1
  hunit  = 3
  return
 endif
 iistep=iftimes(2)-iftimes(1)
 
 if(iistep .ge. 6000) then 
  iistep=int(float(iistep)/10000.)
  hunit=3 ! horas
  return
 endif
 if(iistep.ge.60) then
  iistep=int(float(iistep)/100.)
  hunit=2 !min
  return
 endif
 if(iistep.lt.60) then
  hunit=1 !seg
  return
 endif

return

entry RAMS_get_idata(nopt,nfl,ngr,idata,nval)

if(nopt.eq.0)then
   nval=4
   do k=1,nval
      idata(k)=nfgpnts(k,ngr,nfl)
   enddo
elseif(nopt.eq.1)then
   nval=1
   idata(1)=nfgrids(nfl)
elseif(nopt.eq.2)then
   nval=1
   idata(1)=ifdates(nfl)
elseif(nopt.eq.3)then
   nval=1
   idata(1)=iftimes(nfl)
endif
return

entry RAMS_get_fdata(nopt,nfl,ngr,fdata,nval)

if(nopt.eq.0)then
   nval=nfgpnts(3,ngr,nfl)
   do k=1,nval
      fdata(k)=flevels(k,ngr,nfl)
   enddo
elseif(nopt.eq.1)then
   nval=1
   fdata(1)=ftimes(nfl)
elseif(nopt.eq.2)then
   nval=1
   fdata(1)=startutc
   print*,'RAMS_get_fdata: startutc',startutc
elseif(nopt.eq.3)then
   nval=1
   fdata(1)=httop
elseif(nopt.eq.4)then
   nval=1
   fdata(1)=fdelx(ngr,nfl)
elseif(nopt.eq.5)then
   nval=1
   fdata(1)=fdely(ngr,nfl)
endif

return

entry RAMS_get_cdata(nopt,nfl,cdata,nval)

if(nopt.eq.0)then
   cdata(1)=fnames(nfl)
   nval=1
elseif(nopt.eq.1)then
   cdata(1)=fnamesp(nfl)
   nval=1
endif

return

end
!***************************************************************************

integer function RAMS_getvar(stringg,itype,ngrd,a,b,flnm)

use an_header

implicit none
!include 'interface.h'

real :: a(*),b(*)
integer :: itype,ngrd,il,lastchar,ill
character*(*) flnm,cgrid*1,flng*240,errmsg*120,stringg,string*20
logical there
integer :: ierr_getvar,ifound,ni,npts,iword
common /getvar/ierr_getvar,ifound

!print*,'----------------------------'
!print*,stringg,string
!print*,'----------------------------'
!!!!!!!!!!!!
!para leitura de analises de medias
il=lastchar(flnm)
!print*,il,flnm
!print*,'---- ',flnm(il-18:il-18)


ill=lastchar(stringg)

if(flnm(il-18:il-18).eq.'M') then
string=stringg(1:ill)//'M'
!print*,ill,string(1:ill+1),' ',string,' ',char(0)

else

string=stringg(1:ill)
endif

do ni=1,nvbtab


   if(string.eq.anal_table(ni)%string.and.ngrd.eq.anal_table(ni)%ngrid) then
      write(cgrid,'(i1)') ngrd
      flng=flnm//'-g'//cgrid//'.vfm'//char(0)
      inquire(file=flng,exist=there)
      if(.not.there) then
         errmsg='File not found - '//flng
         call error_mess(errmsg)
         return
      endif
      npts=anal_table(ni)%nvalues
      itype=anal_table(ni)%idim_type
      iword=anal_table(ni)%npointer
            
            print*,'------------------------------------------'
            print*,'Get Var: ',anal_table(ni)%string,' itype=', itype
            print*,'String size =',npts,' Pointer=',iword
            print*,'------------------------------------------'
      call RAMS_c_open(flng,'r'//char(0))
      call vfirecr(10,a,npts,'LIN',b,iword)
      call RAMS_c_close()


      RAMS_getvar=0
      ifound=ifound+1
      return

   endif
enddo

errmsg='Variable not available in this run - '//string
call error_mess(errmsg)
RAMS_getvar=1
ierr_getvar=1

return
end
!************
subroutine lixo(a,n)
dimension a(n)
do i=1,n
print*,i,a(i)
enddo
!stop
return
end
!************
!
!***************************************************************************
!############################# Change Log ##################################
! 2.3.0.3
!
! 000830 CJT RAMS_varlib ##
!            Removed interface.h ##
! 000829 CJT RAMS_varlib ##
!            Changed from passing arguments to using an_header module ##
! 000828 MJB RAMS_varlib ##
!            Replaced c dynamic allocations to f90. ##
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!  Mission Research Corporation / *ASTeR Division
!###########################################################################

subroutine RAMS_varlib(cvar,n1,n2,n3,n4,n5,n6,ngrd,flnm  &
                 ,cdname,cdunits,ivar_type,a,b,a2,a6)
use rconstants
implicit none
include 'rcommons.h'
integer :: n1,n2,n3,ngrd,n4,n5,n6

integer :: memsiz4

common /mem/memsiz4

character*(*) cvar,flnm,cdname,cdunits

real :: a(*),b(*),a2(*),a6(*)

real, allocatable :: c(:),d(:),e(:),f(:),h(:)
integer :: lv,lv2,idim_type,irecind,irecsize,irecsizep,ind,ispec
integer :: memsave4,ierr,kp

integer, external :: RAMS_getvar, lastchar, irfree, iralloc
integer :: ierr_getvar,ifound,ivar_type
real  f1(n1,n2,n3),f2(n1,n2,n3)
common /getvar/ierr_getvar,ifound
data memsave4/0/

real, allocatable :: pv1(:,:,:),pv2(:,:,:),pv3(:,:,:),pv4(:,:,:)
real, allocatable :: pv5(:,:,:)

integer, parameter :: nwave=50
write(unit=*,fmt=*) memsiz4,memsave4

if (memsiz4 > memsave4) then
   if (allocated(c)) deallocate (c)
   if (allocated(d)) deallocate (d)
   if (allocated(e)) deallocate (e)
   if (allocated(f)) deallocate (f)
   if (allocated(h)) deallocate (h)
   allocate(c(memsiz4))
   allocate(d(memsiz4))
   allocate(e(memsiz4))
   allocate(f(memsiz4))
   allocate(h(memsiz4))
   memsave4 = memsiz4
endif

lv=lastchar(cvar)
lv2=min(lv,index(cvar,':')-1)  ! for HYPACT fields

ivar_type=0
ierr_getvar=0
ierr=0
ifound=0

! 3D VELOCITY AND VORTICITY VARIABLES

if(cvar(1:lv).eq.'u') then
   ivar_type=3
   ierr= RAMS_getvar('UP',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname='u'
   cdunits='m/s'

elseif(cvar(1:lv).eq.'v') then
   ivar_type=3
   ierr= RAMS_getvar('VP',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname='v'
   cdunits='m/s'

elseif(cvar(1:lv).eq.'zitheta') then
    ivar_type=2
    ierr= RAMS_getvar('THETA',idim_type,ngrd,c,b,flnm)
    ierr_getvar = ierr_getvar + ierr
    ierr= RAMS_getvar('RCP',idim_type,ngrd,e,b,flnm)
    ierr_getvar = ierr_getvar + ierr
    call get_ZItheta(n1,n2,n3,a,c,e,ngrd)
    cdname='Height PBL'
    cdunits='m -sigmaz'

elseif(cvar(1:lv).eq.'ue') then
   ivar_type=3
   ierr= RAMS_getvar('UP',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr= RAMS_getvar('VP',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_rotate(n1,n2,n3,a,c,ngrd)
   cdname='ue'
   cdunits='m/s'

elseif(cvar(1:lv).eq.'ve') then
   ivar_type=3
   ierr= RAMS_getvar('VP',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr= RAMS_getvar('UP',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_rotate(n1,n2,n3,c,a,ngrd)
   cdname='ve'
   cdunits='m/s'

elseif(cvar(1:lv).eq.'ue_avg') then
   ivar_type=3
   ierr=RAMS_getvar('UP',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr=RAMS_getvar('VP',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_rotate(n1,n2,n3,a,c,ngrd)
   call RAMS_comp_avgu(n1,n2,n3,a)
   cdname='ue_avg'
   cdunits='m/s'

elseif(cvar(1:lv).eq.'ve_avg') then
   ivar_type=3
   ierr=RAMS_getvar('VP',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr=RAMS_getvar('UP',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_rotate(n1,n2,n3,c,a,ngrd)
   call RAMS_comp_avgv(n1,n2,n3,a)
   cdname='ve_avg'
   cdunits='m/s'

elseif(cvar(1:lv).eq.'w') then
   ivar_type=3
   ierr= RAMS_getvar('WP',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname='w'
   cdunits='m/s'

elseif(cvar(1:lv).eq.'wcms') then
   ivar_type=3
   ierr= RAMS_getvar('WP',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_wcms(n1,n2,n3,a)
   cdname='w'
   cdunits='cm/s'

elseif(cvar(1:lv).eq.'w_avg') then
   ivar_type=3
   ierr=RAMS_getvar('WP',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_avgw(n1,n2,n3,a)
   cdname='w_avg'
   cdunits='m/s'

elseif(cvar(1:lv).eq.'speed') then
   ivar_type=3
   ierr= RAMS_getvar('UP',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr= RAMS_getvar('VP',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_speed(n1,n2,n3,a,c)
   cdname='speed'
   cdunits='m/s'

elseif(cvar(1:lv).eq.'speed_mph') then
   ivar_type=3
   ierr= RAMS_getvar('UP',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr= RAMS_getvar('VP',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_speed(n1,n2,n3,a,c)
   call RAMS_comp_mults(n1,n2,n3,a,2.237)
   cdname='speed'
   cdunits='mph'

elseif(cvar(1:lv).eq.'tempc2m') then
   ivar_type=2
   ierr=RAMS_getvar('UP',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr=RAMS_getvar('VP',idim_type,ngrd,d,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_speed(n1,n2,n3,c,d)
   ierr=RAMS_getvar('THETA',idim_type,ngrd,d,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr=RAMS_getvar('PI',idim_type,ngrd,f,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr=RAMS_getvar('TOPT',idim_type,ngrd,e,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   allocate (pv1(nnxp(ngrd),nnyp(ngrd),npatch))
   allocate (pv2(nnxp(ngrd),nnyp(ngrd),npatch))
   allocate (pv3(nnxp(ngrd),nnyp(ngrd),npatch))
   allocate (pv4(nnxp(ngrd),nnyp(ngrd),npatch))
   allocate (pv5(nnxp(ngrd),nnyp(ngrd),npatch))
   ierr=RAMS_getvar('USTAR'      ,idim_type,ngrd,pv1,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr=RAMS_getvar('VEG_ROUGH'  ,idim_type,ngrd,pv2,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr=RAMS_getvar('CAN_THETA'   ,idim_type,ngrd,pv3,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr=RAMS_getvar('CAN_PRSS'   ,idim_type,ngrd,pv4,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_theta2temp(n1,n2,npatch,pv3,pv4)

   ierr=RAMS_getvar('PATCH_AREA' ,idim_type,ngrd,pv4,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr=RAMS_getvar('TSTAR'      ,idim_type,ngrd,pv5,b,flnm)
   ierr_getvar = ierr_getvar + ierr
 
    call RAMS_reduced_temp(nnxp(ngrd),nnyp(ngrd),nnzp(ngrd),npatch       &
                          ,a,c,pv1,pv5  &
			  ,2.,ztn(2,ngrd),pv2,pv4,pv3,d,f,e  &
                          ,zmn(nnzp(1)-1,1))
   deallocate (pv1,pv2,pv3,pv4,pv5)
   call RAMS_comp_tempK(n1,n2,1,a,f) !converte de temp. potencial para
                                     !temperatura termodinamica    
   call RAMS_comp_tempC(n1,n2,1,a)
   cdname='temp - 2m AGL;'
   cdunits='C'

!

elseif(cvar(1:lv).eq.'speed10m') then

   ivar_type=2
   ierr=RAMS_getvar('UP',idim_type,ngrd  &
         ,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr=RAMS_getvar('VP',idim_type,ngrd  &
         ,d,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_speed(n1,n2,n3,c,d)
   ierr=RAMS_getvar('THETA',idim_type,ngrd  &
         ,d,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr=RAMS_getvar('PI',idim_type,ngrd  &
         ,f,b,flnm)
   ierr_getvar = ierr_getvar + ierr
!        Get topo
   ierr= RAMS_getvar('TOPT',idim_type,ngrd  &
         ,e,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   allocate (pv1(nnxp(ngrd),nnyp(ngrd),npatch) )
   allocate (pv2(nnxp(ngrd),nnyp(ngrd),npatch) )
   allocate (pv3(nnxp(ngrd),nnyp(ngrd),npatch) )
   allocate (pv4(nnxp(ngrd),nnyp(ngrd),npatch) )

!           Get ustar
   ierr = RAMS_getvar('USTAR',idim_type,ngrd,pv1,b,flnm)        
   ierr_getvar = ierr_getvar + ierr
!           Get net roughness
   ierr = RAMS_getvar('PATCH_ROUGH',idim_type,ngrd   &
        ,pv2,b,flnm)        
   ierr_getvar = ierr_getvar + ierr
!           Get patch canopy temperature
   ierr = RAMS_getvar('CAN_THETA',idim_type,ngrd   &
        ,pv3,b,flnm)       
   ierr_getvar = ierr_getvar + ierr
   ierr = RAMS_getvar('CAN_PRSS',idim_type,ngrd   &
        ,pv4,b,flnm)       
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_theta2temp(n1,n2,npatch,pv3,pv4)
!           Get % coverage
   ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
        ,pv4,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_reduced_wind(nnxp(ngrd),nnyp(ngrd)  &
      ,nnzp(ngrd),npatch  &
      ,a,c,pv1,10.,ztn(2,ngrd)  &
      ,pv2,pv4,pv3,d,f,e  &
      ,zmn(nnzp(1)-1,1))

   deallocate (pv1,pv2,pv3,pv4)
   cdname='speed - 10m AGL'
   cdunits='m/s'

elseif(cvar(1:lv).eq.'direction') then
   ivar_type=3
   ierr= RAMS_getvar('UP',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr= RAMS_getvar('VP',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_dir(n1,n2,n3,a,c,ngrd)
   cdname='direction'
   cdunits='deg'

elseif(cvar(1:lv).eq.'pvap') then
   ivar_type=3
   ierr= RAMS_getvar('PI',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_press(n1,n2,n3,c)                   ! c is pressure [hPa]
   ierr= RAMS_getvar('RV',idim_type,ngrd,a,b,flnm)    ! a is mixing ratio
   call RAMS_comp_noneg(n1,n2,n3,a)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_pvap(n1,n2,n3,c,a)                  ! a is vapour pressure [hPa]
   cdname='water vapour pressure'
   cdunits='hPa'

elseif(cvar(1:lv).eq.'alpha') then
   ivar_type=3
   ierr= RAMS_getvar('PI',idim_type,ngrd,c,b,flnm)    ! c is Exner function
   ierr_getvar = ierr_getvar + ierr
   ierr= RAMS_getvar('THETA',idim_type,ngrd,a,b,flnm) ! a is potential temperature [K]
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_tempK(n1,n2,n3,a,c)                 ! a is temperature [K]
   ierr= RAMS_getvar('RV',idim_type,ngrd,d,b,flnm)    ! d is mixing ratio [kg/kg]
   call RAMS_comp_noneg(n1,n2,n3,d)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_press(n1,n2,n3,c)                   ! c is pressure [hPa]
   call RAMS_comp_mults(n1,n2,n3,c,100.)              ! c is pressure [Pa]
   call RAMS_comp_pvap(n1,n2,n3,c,d)                  ! d is vapour pressure [hPa]
   call RAMS_comp_spvol(n1,n2,n3,a,d,c)               ! a is specific volume [m3/kg]
   cdname='specific volume'
   cdunits='m3/kg'

elseif(cvar(1:lv).eq.'solenoidx') then
   ivar_type=3
   ierr= RAMS_getvar('PI',idim_type,ngrd,c,b,flnm)    ! c is Exner function
   ierr_getvar = ierr_getvar + ierr
   ierr= RAMS_getvar('THETA',idim_type,ngrd,e,b,flnm) ! e is potential temperature [K]
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_tempK(n1,n2,n3,e,c)                 ! e is temperature [K]
   ierr= RAMS_getvar('RV',idim_type,ngrd,d,b,flnm)    ! d is mixing ratio [kg/kg]
   call RAMS_comp_noneg(n1,n2,n3,d)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_press(n1,n2,n3,c)                   ! c is pressure [hPa]
   call RAMS_comp_mults(n1,n2,n3,c,100.)              ! c is pressure [Pa]
   call RAMS_comp_pvap(n1,n2,n3,c,d)                  ! d is vapour pressure [hPa]
   call RAMS_comp_spvol(n1,n2,n3,e,d,c)               ! e is specific volume [m3/kg]
   ierr= RAMS_getvar('TOPT',idim_type,ngrd,h,b,flnm)  ! h is topography [m]
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_solenoidx(n1,n2,n3,e,c,a,h,ngrd)
   cdname='x-solenoid term'
   cdunits='rad/s2'

elseif(cvar(1:lv).eq.'solenoidy') then
   ivar_type=3
   ierr= RAMS_getvar('PI',idim_type,ngrd,c,b,flnm)    ! c is Exner function
   ierr_getvar = ierr_getvar + ierr
   ierr= RAMS_getvar('THETA',idim_type,ngrd,e,b,flnm) ! e is potential temperature [K]
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_tempK(n1,n2,n3,e,c)                 ! e is temperature [K]
   ierr= RAMS_getvar('RV',idim_type,ngrd,d,b,flnm)    ! d is mixing ratio [kg/kg]
   call RAMS_comp_noneg(n1,n2,n3,d)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_press(n1,n2,n3,c)                   ! c is pressure [hPa]
   call RAMS_comp_mults(n1,n2,n3,c,100.)              ! c is pressure [Pa]
   call RAMS_comp_pvap(n1,n2,n3,c,d)                  ! d is vapour pressure [hPa]
   call RAMS_comp_spvol(n1,n2,n3,e,d,c)               ! e is specific volume [m3/kg]
   ierr= RAMS_getvar('TOPT',idim_type,ngrd,h,b,flnm)  ! h is topography [m]
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_solenoidy(n1,n2,n3,e,c,a,h,ngrd)
   cdname='y-solenoid term'
   cdunits='rad/s2'

elseif(cvar(1:lv).eq.'relvortx') then
   ivar_type=3
   ierr= RAMS_getvar('VP',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr= RAMS_getvar('WP',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr= RAMS_getvar('TOPT',idim_type,ngrd,d,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_relvortx(n1,n2,n3,a,c,b,d,ngrd)
   cdname='x-vorticity'
   cdunits='rad/s'

elseif(cvar(1:lv).eq.'relvorty') then
   ivar_type=3
   ierr= RAMS_getvar('UP',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr= RAMS_getvar('WP',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr= RAMS_getvar('TOPT',idim_type,ngrd,d,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_relvorty(n1,n2,n3,a,c,b,d,ngrd)
   cdname='y-vorticity'
   cdunits='rad/s'

elseif(cvar(1:lv).eq.'relvortz') then
   ivar_type=3 
   ierr= RAMS_getvar('UP',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr= RAMS_getvar('VP',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr= RAMS_getvar('TOPT',idim_type,ngrd,d,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_relvortz(n1,n2,n3,a,c,b,d,ngrd)
   cdname='relative z-vorticity'
   cdunits='rad/s'

elseif(cvar(1:lv).eq.'absvortz') then
   ivar_type=3
   ierr= RAMS_getvar('UP',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr= RAMS_getvar('VP',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr= RAMS_getvar('TOPT',idim_type,ngrd,d,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_totvortz(n1,n2,n3,a,c,b,d,ngrd)
   cdname='absolute z-vorticity'
   cdunits='rad/s'

elseif(cvar(1:lv).eq.'potvortz') then
   ivar_type=3
   ierr= RAMS_getvar('UP',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr= RAMS_getvar('VP',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr= RAMS_getvar('TOPT',idim_type,ngrd,d,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_totvortz(n1,n2,n3,a,c,b,d,ngrd)
   call RAMS_comp_dn0(n1,n2,n3,e,b,c,d,ngrd)

   ierr= RAMS_getvar('THETA',idim_type,ngrd,b,e,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_potvortz(n1,n2,n3,a,b,c,e,d,ngrd)
   cdname='potential z-vorticity'
   cdunits='rad/s'

elseif(cvar(1:lv).eq.'horiz_div') then
   ivar_type=3
   ierr= RAMS_getvar('WP',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_horizdiv(n1,n2,n3,a)
   cdname='horizontal divergence'
   cdunits='/s'

! 3D THERMODYNAMIC PROPERTIES OF AIR

elseif(cvar(1:lv).eq.'pi') then
   ivar_type=3
   ierr= RAMS_getvar('PI',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname='Exner function'
   cdunits='J/(kg K)'

elseif(cvar(1:lv).eq.'press') then
   ivar_type=3
   ierr= RAMS_getvar('PI',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_press(n1,n2,n3,a)
   cdname='pressure'
   cdunits='mb'

elseif(cvar(1:lv).eq.'theta') then
   ivar_type=3
   ierr= RAMS_getvar('THETA',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname='potential temp'
   cdunits='K'

elseif(cvar(1:lv).eq.'thil') then
   ivar_type=3
   ierr= RAMS_getvar('THP',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname='ice-liquid potential temp'
   cdunits='K'

elseif(cvar(1:lv).eq.'co2') then
   ivar_type=3
   ierr= RAMS_getvar('CO2P',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname='carbon dioxide mixing ratio'
   cdunits='umol/mol'

elseif(cvar(1:lv).eq.'dn0') then
   ivar_type=3
   ierr= RAMS_getvar('TOPT',idim_type,ngrd,e,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_dn0(n1,n2,n3,c,b,a,e,ngrd)
   cdname='ref density'

   cdunits='kg/m3'

elseif(cvar(1:lv).eq.'pi0') then
   ivar_type=3
   ierr= RAMS_getvar('TOPT',idim_type,ngrd,e,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_dn0(n1,n2,n3,a,b,c,e,ngrd)
   cdname='ref Exner func'
   cdunits='J/(kg K)'

elseif(cvar(1:lv).eq.'th0') then
   ivar_type=3
   ierr= RAMS_getvar('TOPT',idim_type,ngrd,e,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_dn0(n1,n2,n3,b,a,c,e,ngrd)
   cdname='reference virtual potential temp'
   cdunits='K'

elseif(cvar(1:lv).eq.'pert_pressure') then
   ivar_type=3
   ierr= RAMS_getvar('TOPT',idim_type,ngrd,e,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_dn0(n1,n2,n3,c,a,b,e,ngrd)
   ierr= RAMS_getvar('PI',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   if (ierr.eq.0) call RAMS_comp_ppress(n1,n2,n3,a,c)
   cdname='pert pressure'
   cdunits='mb'

elseif(cvar(1:lv).eq.'tempk') then
   ivar_type=3
   ierr= RAMS_getvar('THETA',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr= RAMS_getvar('PI',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_tempK(n1,n2,n3,a,c)
   cdname='temperature'
   cdunits='K'

elseif(cvar(1:lv).eq.'tempc') then
   ivar_type=3
   ierr= RAMS_getvar('THETA',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr= RAMS_getvar('PI',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_tempK(n1,n2,n3,a,c)
   call RAMS_comp_tempC(n1,n2,n3,a)
   cdname='temperature'
   cdunits='C'

elseif(cvar(1:lv).eq.'tempf') then
   ivar_type=3
   ierr= RAMS_getvar('THETA',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr= RAMS_getvar('PI',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_tempK(n1,n2,n3,a,c)
   call RAMS_comp_tempF(n1,n2,n3,a)
   cdname='temperature'
   cdunits='F'

!<Demerval

elseif(cvar(1:lv).eq.'theta_e') then
   ivar_type=3
   ierr= RAMS_getvar('RV',idim_type,ngrd,e,b,flnm)
   ierr_getvar = ierr_getvar + ierr
!!  e=rvap
!
   ierr= RAMS_getvar('THETA',idim_type,ngrd,f1,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr= RAMS_getvar('PI',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_tempK(n1,n2,n3,f1,c)
!
!! f1=tempk
!
   ivar_type=3
   ierr= RAMS_getvar('RV',idim_type,ngrd,f2,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr= RAMS_getvar('PI',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr= RAMS_getvar('THETA',idim_type,ngrd,d,b,flnm)
   ierr_getvar = ierr_getvar + ierr
!
   call RAMS_comp_rh(n1,n2,n3,f2,c,d)
!
!
!! f2=umidade relativa %
!
   ierr= RAMS_getvar('THETA',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
!! a=theta
!
   call RAMS_comp_thete(n1,n2,n3,a,e,f1,f2)
   cdname='Equiv pot temp;'
   cdunits='K'


elseif(cvar(1:lv).eq.'thetae_iv' .or. cvar(1:lv).eq.'theiv') then
   ivar_type=3
   call RAMS_comp_tempK(n1,n2,n3,f,c)
   ierr= RAMS_getvar('THP',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr= RAMS_getvar('PI',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_press(n1,n2,n3,c)
   ierr= RAMS_getvar('THETA',idim_type,ngrd,f,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr= RAMS_getvar('RV',idim_type,ngrd,d,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr= RAMS_getvar('RTP',idim_type,ngrd,e,b,flnm)
   ierr_getvar = ierr_getvar + ierr

   call RAMS_comp_thetaeiv(n1,n2,n3,a,f,c,d,e)
   cdname='Ice-vapour equivt pot temp'
   cdunits='K'

elseif(cvar(1:lv).eq.'theta_v') then
   ivar_type=3
   ierr= RAMS_getvar('THETA',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr= RAMS_getvar('RV',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr= RAMS_getvar('RTP',idim_type,ngrd,d,b,flnm)
   ierr_getvar = ierr_getvar + ierr

   call RAMS_comp_thetv(n1,n2,n3,a,c,d)
   cdname='virtual pot temp'
   cdunits='K'

! 3D MOISTURE MASS MIXING RATIOS AND HUMIDITY

elseif(cvar(1:lv).eq.'rv') then
   ivar_type=3
   ierr= RAMS_getvar('RV',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   if(ierr.eq.0) then
      call RAMS_comp_mults(n1,n2,n3,a,1.e3)
      call RAMS_comp_noneg(n1,n2,n3,a)
   endif
   cdname='vapor mix ratio'
   cdunits='g/kg'

elseif(cvar(1:lv).eq.'cloud') then
   ivar_type=3
   ierr= RAMS_getvar('RCP',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   if(ierr.eq.0) then
      call RAMS_comp_mults(n1,n2,n3,a,1.e3)
      call RAMS_comp_noneg(n1,n2,n3,a)
   endif
   cdname='cloud mix ratio'
   cdunits='g/kg'

elseif(cvar(1:lv).eq.'rain') then
   ivar_type=3
   ierr= RAMS_getvar('RRP',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   if(ierr.eq.0) then
      call RAMS_comp_mults(n1,n2,n3,a,1.e3)
      call RAMS_comp_noneg(n1,n2,n3,a)
   endif
   cdname='rain mix ratio'
   cdunits='g/kg'

elseif(cvar(1:lv).eq.'pristine') then
   ivar_type=3
   ierr= RAMS_getvar('RPP',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   if(ierr.eq.0) then
      call RAMS_comp_mults(n1,n2,n3,a,1.e3)
      call RAMS_comp_noneg(n1,n2,n3,a)
   endif
   cdname='pristine mix ratio'
   cdunits='g/kg'

elseif(cvar(1:lv).eq.'snow') then
   ivar_type=3
   ierr= RAMS_getvar('RSP',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   if(ierr.eq.0) then
      call RAMS_comp_mults(n1,n2,n3,a,1.e3)
      call RAMS_comp_noneg(n1,n2,n3,a)
   endif
   cdname='snow mix ratio'
   cdunits='g/kg'

elseif(cvar(1:lv).eq.'aggregates' .or. cvar(1:lv).eq.'agg') then
   ivar_type=3
   ierr= RAMS_getvar('RAP',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   if(ierr.eq.0) then
      call RAMS_comp_mults(n1,n2,n3,a,1.e3)
      call RAMS_comp_noneg(n1,n2,n3,a)
   endif
   cdname='aggregate mix ratio'
   cdunits='g/kg'

elseif(cvar(1:lv).eq.'graupel') then
   ivar_type=3
   ierr= RAMS_getvar('RGP',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   if(ierr.eq.0) then
      call RAMS_comp_mults(n1,n2,n3,a,1.e3)
      call RAMS_comp_noneg(n1,n2,n3,a)
   endif
   cdname='graupel mix ratio'
   cdunits='g/kg'

elseif(cvar(1:lv).eq.'hail') then
   ivar_type=3
   ierr= RAMS_getvar('RHP',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   if(ierr.eq.0) then
      call RAMS_comp_mults(n1,n2,n3,a,1.e3)
      call RAMS_comp_noneg(n1,n2,n3,a)
   endif
   cdname='hail mix ratio'
   cdunits='g/kg'

elseif(cvar(1:lv).eq.'liquid') then
   ivar_type=3
   call RAMS_comp_zero(n1,n2,n3,a)
   ierr= RAMS_getvar('RCP',idim_type,ngrd,c,b,flnm)
   if(ierr.eq.0) call RAMS_comp_accum(n1,n2,n3,a,c)
   ierr= RAMS_getvar('RRP',idim_type,ngrd,c,b,flnm)
   if(ierr.eq.0) call RAMS_comp_accum(n1,n2,n3,a,c)

   ierr= RAMS_getvar('RGP',idim_type,ngrd,c,b,flnm)
   if(ierr.eq.0) then
      ierr= RAMS_getvar('Q6',idim_type,ngrd,d,b,flnm)
      if(ierr.eq.0) then
         call RAMS_comp_fracliq(n1,n2,n3,d)
         call RAMS_comp_mult(n1,n2,n3,c,d)
      endif
      call RAMS_comp_accum(n1,n2,n3,a,c)
   endif

   ierr= RAMS_getvar('RHP',idim_type,ngrd,c,b,flnm)
   if(ierr.eq.0) then
      ierr= RAMS_getvar('Q7',idim_type,ngrd,d,b,flnm)
      if(ierr.eq.0) then
         call RAMS_comp_fracliq(n1,n2,n3,d)
         call RAMS_comp_mult(n1,n2,n3,c,d)
      endif
      call RAMS_comp_accum(n1,n2,n3,a,c)
    endif

   call RAMS_comp_mults(n1,n2,n3,a,1.e3)
   call RAMS_comp_noneg(n1,n2,n3,a)
   cdname='liquid mix ratio'
   cdunits='g/kg'

elseif(cvar(1:lv).eq.'ice') then
   ivar_type=3
   call RAMS_comp_zero(n1,n2,n3,a)
   ierr= RAMS_getvar('RPP',idim_type,ngrd,c,b,flnm)
   if(ierr.eq.0) call RAMS_comp_accum(n1,n2,n3,a,c)
   ierr= RAMS_getvar('RSP',idim_type,ngrd,c,b,flnm)
   if(ierr.eq.0) call RAMS_comp_accum(n1,n2,n3,a,c)
   ierr= RAMS_getvar('RAP',idim_type,ngrd,c,b,flnm)
   if(ierr.eq.0) call RAMS_comp_accum(n1,n2,n3,a,c)

   ierr= RAMS_getvar('RGP',idim_type,ngrd,c,b,flnm)
   if(ierr.eq.0) then
      ierr= RAMS_getvar('Q6',idim_type,ngrd,d,b,flnm)
      if(ierr.eq.0) then
         call RAMS_comp_fracice(n1,n2,n3,d)
         call RAMS_comp_mult(n1,n2,n3,c,d)
      endif
      call RAMS_comp_accum(n1,n2,n3,a,c)
   endif

   ierr= RAMS_getvar('RHP',idim_type,ngrd,c,b,flnm)
   if(ierr.eq.0) then
      ierr= RAMS_getvar('Q7',idim_type,ngrd,d,b,flnm)
      if(ierr.eq.0) then
         call RAMS_comp_fracice(n1,n2,n3,d)
         call RAMS_comp_mult(n1,n2,n3,c,d)
      endif
      call RAMS_comp_accum(n1,n2,n3,a,c)
   endif

   call RAMS_comp_mults(n1,n2,n3,a,1.e3)
   call RAMS_comp_noneg(n1,n2,n3,a)
   cdname='ice mix ratio'
   cdunits='g/kg'

elseif(cvar(1:lv).eq.'total_cond') then
   ivar_type=3
   call RAMS_comp_zero(n1,n2,n3,a)
   ierr= RAMS_getvar('RCP',idim_type,ngrd,c,b,flnm)
   if(ierr.eq.0) call RAMS_comp_accum(n1,n2,n3,a,c)
   ierr= RAMS_getvar('RRP',idim_type,ngrd,c,b,flnm)
   if(ierr.eq.0) call RAMS_comp_accum(n1,n2,n3,a,c)
   ierr= RAMS_getvar('RPP',idim_type,ngrd,c,b,flnm)
   if(ierr.eq.0) call RAMS_comp_accum(n1,n2,n3,a,c)
   ierr= RAMS_getvar('RSP',idim_type,ngrd,c,b,flnm)
   if(ierr.eq.0) call RAMS_comp_accum(n1,n2,n3,a,c)
   ierr= RAMS_getvar('RAP',idim_type,ngrd,c,b,flnm)
   if(ierr.eq.0) call RAMS_comp_accum(n1,n2,n3,a,c)
   ierr= RAMS_getvar('RGP',idim_type,ngrd,c,b,flnm)
   if(ierr.eq.0) call RAMS_comp_accum(n1,n2,n3,a,c)
   ierr= RAMS_getvar('RHP',idim_type,ngrd,c,b,flnm)
   if(ierr.eq.0) call RAMS_comp_accum(n1,n2,n3,a,c)

   call RAMS_comp_mults(n1,n2,n3,a,1.e3)
   call RAMS_comp_noneg(n1,n2,n3,a)
   cdname='cloud mix ratio'
   cdunits='g/kg'

elseif(cvar(1:lv).eq.'rtotal') then
   ivar_type=3
   ierr= RAMS_getvar('RV',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr= RAMS_getvar('RCP',idim_type,ngrd,c,b,flnm)
   if(ierr.eq.0) call RAMS_comp_accum(n1,n2,n3,a,c)
   ierr= RAMS_getvar('RRP',idim_type,ngrd,c,b,flnm)
   if(ierr.eq.0) call RAMS_comp_accum(n1,n2,n3,a,c)
   ierr= RAMS_getvar('RPP',idim_type,ngrd,c,b,flnm)
   if(ierr.eq.0) call RAMS_comp_accum(n1,n2,n3,a,c)
   ierr= RAMS_getvar('RSP',idim_type,ngrd,c,b,flnm)
   if(ierr.eq.0) call RAMS_comp_accum(n1,n2,n3,a,c)
   ierr= RAMS_getvar('RAP',idim_type,ngrd,c,b,flnm)
   if(ierr.eq.0) call RAMS_comp_accum(n1,n2,n3,a,c)
   ierr= RAMS_getvar('RGP',idim_type,ngrd,c,b,flnm)
   if(ierr.eq.0) call RAMS_comp_accum(n1,n2,n3,a,c)
   ierr= RAMS_getvar('RHP',idim_type,ngrd,c,b,flnm)
   if(ierr.eq.0) call RAMS_comp_accum(n1,n2,n3,a,c)

   call RAMS_comp_mults(n1,n2,n3,a,1.e3)
   call RAMS_comp_noneg(n1,n2,n3,a)
   cdname='total mix ratio'
   cdunits='g/kg'

elseif(cvar(1:lv).eq.'rtotal_orig') then
   ivar_type=3
   ierr= RAMS_getvar('RTP',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_mults(n1,n2,n3,a,1.e3)
   call RAMS_comp_noneg(n1,n2,n3,a)
   cdname='orig rtotal'
   cdunits='g/kg'

elseif(cvar(1:lv).eq.'dewptk') then
   ivar_type=3
   ivar_type=3
   ierr= RAMS_getvar('RV',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr= RAMS_getvar('PI',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr= RAMS_getvar('THETA',idim_type,ngrd,d,b,flnm)
   ierr_getvar = ierr_getvar + ierr

   call RAMS_comp_dewK(n1,n2,n3,a,c,d)
   call RAMS_comp_tempK(n1,n2,n3,a,c)
   cdname='dewpoint temp'
   cdunits='K'

elseif(cvar(1:lv).eq.'dewptf') then
   ivar_type=3
   ierr= RAMS_getvar('RV',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr= RAMS_getvar('PI',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr= RAMS_getvar('THETA',idim_type,ngrd,d,b,flnm)
   ierr_getvar = ierr_getvar + ierr

   call RAMS_comp_dewK(n1,n2,n3,a,c,d)
   call RAMS_comp_tempF(n1,n2,n3,a)
   cdname='dewpoint temp'
   cdunits='F'

elseif(cvar(1:lv).eq.'dewptc') then
   ivar_type=3
   ierr= RAMS_getvar('RV',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr= RAMS_getvar('PI',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr= RAMS_getvar('THETA',idim_type,ngrd,d,b,flnm)
   ierr_getvar = ierr_getvar + ierr

   call RAMS_comp_dewK(n1,n2,n3,a,c,d)
   call RAMS_comp_tempC(n1,n2,n3,a)
   cdname='dewpoint temp'
   cdunits='C'

elseif(cvar(1:lv).eq.'rh') then
   ivar_type=3
   ierr= RAMS_getvar('RV',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr= RAMS_getvar('PI',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr= RAMS_getvar('THETA',idim_type,ngrd,d,b,flnm)
   ierr_getvar = ierr_getvar + ierr

   call RAMS_comp_rh(n1,n2,n3,a,c,d)
   call RAMS_comp_noneg(n1,n2,n3,a)
   cdname='relative humidity'
   cdunits='pct'

elseif(cvar(1:lv).eq.'clear_frac') then
   ivar_type=2
   ierr= RAMS_getvar('RV',idim_type,ngrd,b,a,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr= RAMS_getvar('PI',idim_type,ngrd,c,a,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr= RAMS_getvar('THETA',idim_type,ngrd,d,a,flnm)
   ierr_getvar = ierr_getvar + ierr

   call RAMS_comp_rh(n1,n2,n3,b,c,d)
   call RAMS_comp_noneg(n1,n2,n3,b)
   
   call cldfraction(n1,n2,n3,a,c,b)
   
   cdname='clear sky'
   cdunits='frac'

! 3D HYDROMETEOR, CCN, CN, Dep N, AND NONHYGROSCOPIC AEROSOL NUMBER CONCEN

elseif(cvar(1:lv).eq.'cloud_concen_mg') then
   ivar_type=3
! variable 18 is iccp
   ierr= RAMS_getvar('CCP',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_mults(n1,n2,n3,a,1.e-6)
   call RAMS_comp_noneg(n1,n2,n3,a)
   cdname='cloud concen'
   cdunits='#/mg'

elseif(cvar(1:lv).eq.'rain_concen_kg') then
   ivar_type=3
   ierr= RAMS_getvar('CRP',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_noneg(n1,n2,n3,a)
   cdname='rain concen'
   cdunits='#/kg'

elseif(cvar(1:lv).eq.'pris_concen_kg') then
   ivar_type=3
   ierr= RAMS_getvar('CPP',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_noneg(n1,n2,n3,a)
   cdname='pristine concen'
   cdunits='#/kg'

elseif(cvar(1:lv).eq.'snow_concen_kg') then
   ivar_type=3
   ierr= RAMS_getvar('CSP',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_noneg(n1,n2,n3,a)
   cdname='snow concen'
   cdunits='#/kg'

elseif(cvar(1:lv).eq.'agg_concen_kg') then
   ivar_type=3
   ierr= RAMS_getvar('CAP',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_noneg(n1,n2,n3,a)
   cdname='aggregate concen'
   cdunits='#/kg'

elseif(cvar(1:lv).eq.'graup_concen_kg') then
   ivar_type=3
   ierr= RAMS_getvar('CGP',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_noneg(n1,n2,n3,a)
   cdname='graupel concen'
   cdunits='#/kg'

elseif(cvar(1:lv).eq.'hail_concen_kg') then
   ivar_type=3
   ierr= RAMS_getvar('CHP',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_noneg(n1,n2,n3,a)
   cdname='hail concen'
   cdunits='#/kg'

elseif(cvar(1:lv).eq.'cloud_concen_cm3') then
   ivar_type=3
   ierr= RAMS_getvar('CCP',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr= RAMS_getvar('TOPT',idim_type,ngrd,e,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_dn0(n1,n2,n3,b,c,d,e,ngrd)
   call RAMS_comp_mult(n1,n2,n3,a,d)
   call RAMS_comp_mults(n1,n2,n3,a,1.e-6)
   call RAMS_comp_noneg(n1,n2,n3,a)
   cdname='cloud concen'
   cdunits='#/cm3'

elseif(cvar(1:lv).eq.'rain_concen_m3') then
   ivar_type=3
   ierr= RAMS_getvar('CRP',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr= RAMS_getvar('TOPT',idim_type,ngrd,e,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_dn0(n1,n2,n3,b,c,d,e,ngrd)
   call RAMS_comp_mult(n1,n2,n3,a,d)
   call RAMS_comp_noneg(n1,n2,n3,a)
   cdname='rain concen'
   cdunits='#/m3'

elseif(cvar(1:lv).eq.'pris_concen_m3') then
   ivar_type=3
   ierr= RAMS_getvar('CPP',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr= RAMS_getvar('TOPT',idim_type,ngrd,e,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_dn0(n1,n2,n3,b,c,d,e,ngrd)
   call RAMS_comp_mult(n1,n2,n3,a,d)
   call RAMS_comp_noneg(n1,n2,n3,a)
   cdname='pristine concen'
   cdunits='#/m3'

elseif(cvar(1:lv).eq.'snow_concen_m3') then
   ivar_type=3
   ierr= RAMS_getvar('CSP',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr= RAMS_getvar('TOPT',idim_type,ngrd,e,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_dn0(n1,n2,n3,b,c,d,e,ngrd)
   call RAMS_comp_mult(n1,n2,n3,a,d)
   call RAMS_comp_noneg(n1,n2,n3,a)
   cdname='snow concen'
   cdunits='#/m3'

elseif(cvar(1:lv).eq.'agg_concen_m3') then
   ivar_type=3
   ierr= RAMS_getvar('CAP',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr= RAMS_getvar('TOPT',idim_type,ngrd,e,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_dn0(n1,n2,n3,b,c,d,e,ngrd)
   call RAMS_comp_mult(n1,n2,n3,a,d)
   call RAMS_comp_noneg(n1,n2,n3,a)
   cdname='aggregates concen'
   cdunits='#/m3'

elseif(cvar(1:lv).eq.'graup_concen_m3') then
   ivar_type=3
   ierr= RAMS_getvar('CGP',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr= RAMS_getvar('TOPT',idim_type,ngrd,e,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_dn0(n1,n2,n3,b,c,d,e,ngrd)
   call RAMS_comp_mult(n1,n2,n3,a,d)
   call RAMS_comp_noneg(n1,n2,n3,a)
   cdname='graupel concen'
   cdunits='#/m3'

elseif(cvar(1:lv).eq.'hail_concen_m3') then
   ivar_type=3
   ierr= RAMS_getvar('CHP',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr= RAMS_getvar('TOPT',idim_type,ngrd,e,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_dn0(n1,n2,n3,b,c,d,e,ngrd)
   call RAMS_comp_mult(n1,n2,n3,a,d)
   call RAMS_comp_noneg(n1,n2,n3,a)
   cdname='hail concen'
   cdunits='#/m3'

elseif(cvar(1:lv).eq.'ccn_concen') then
   ivar_type=3
   ierr= RAMS_getvar('CCCNP',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_noneg(n1,n2,n3,a)
   call RAMS_comp_mults(n1,n2,n3,a,1.e-6)
   cdname='ccn1 concen'
   cdunits='#/mg'

elseif(cvar(1:lv).eq.'ifn_conc') then
   ivar_type=3
   ierr= RAMS_getvar('CIFNP',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_noneg(n1,n2,n3,a)
   cdname='CN mix ratio'
   cdunits='#/kg'

! 3D HYDROMETEOR DIAMETERS

elseif(cvar(1:lv).eq.'cloud_diam') then
   ivar_type=3
   ierr= RAMS_getvar('RCP',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr= RAMS_getvar('CCP',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_hydrodiam(n1,n2,n3,a,c,cfmas(1),pwmas(1))
   call RAMS_comp_mults(n1,n2,n3,a,1.e6)
   call RAMS_comp_noneg(n1,n2,n3,a)
   cdname='cloud diam'
   cdunits='microns'

elseif(cvar(1:lv).eq.'rain_diam') then
   ivar_type=3
   ierr= RAMS_getvar('RRP',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr= RAMS_getvar('CRP',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_hydrodiam(n1,n2,n3,a,c,cfmas(2),pwmas(2))
   call RAMS_comp_mults(n1,n2,n3,a,1.e3)
   call RAMS_comp_noneg(n1,n2,n3,a)
   cdname='rain diam'
   cdunits='mm'

elseif(cvar(1:lv).eq.'pris_diam') then
   ivar_type=3
   ierr= RAMS_getvar('RPP',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr= RAMS_getvar('CPP',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
! more general case: write habit to anal file for cfmas & pwmas index
   call RAMS_comp_hydrodiam(n1,n2,n3,a,c,cfmas(3),pwmas(3))
   call RAMS_comp_mults(n1,n2,n3,a,1.e6)
   call RAMS_comp_noneg(n1,n2,n3,a)
   cdname='pristine diam'
   cdunits='microns'

elseif(cvar(1:lv).eq.'snow_diam') then
   ivar_type=3
   ierr= RAMS_getvar('RSP',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr= RAMS_getvar('CSP',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
! more general case: write habit to anal file for cfmas & pwmas index
   call RAMS_comp_hydrodiam(n1,n2,n3,a,c,cfmas(4),pwmas(4))
   call RAMS_comp_mults(n1,n2,n3,a,1.e3)
   call RAMS_comp_noneg(n1,n2,n3,a)
   cdname='snow diam'
   cdunits='mm'

elseif(cvar(1:lv).eq.'agg_diam') then
   ivar_type=3
   ierr= RAMS_getvar('RAP',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr= RAMS_getvar('CAP',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_hydrodiam(n1,n2,n3,a,c,cfmas(5),pwmas(5))
   call RAMS_comp_mults(n1,n2,n3,a,1.e3)
   call RAMS_comp_noneg(n1,n2,n3,a)
   cdname='aggregates diam'
   cdunits='mm'

elseif(cvar(1:lv).eq.'graup_diam') then
   ivar_type=3
   ierr= RAMS_getvar('RGP',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr= RAMS_getvar('CGP',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_hydrodiam(n1,n2,n3,a,c,cfmas(6),pwmas(6))
   call RAMS_comp_mults(n1,n2,n3,a,1.e3)
   call RAMS_comp_noneg(n1,n2,n3,a)
   cdname='graupel diam'
   cdunits='mm'

elseif(cvar(1:lv).eq.'hail_diam') then
   ivar_type=3
   ierr= RAMS_getvar('RHP',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr= RAMS_getvar('CHP',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_hydrodiam(n1,n2,n3,a,c,cfmas(7),pwmas(7))
   call RAMS_comp_mults(n1,n2,n3,a,1.e3)
   call RAMS_comp_noneg(n1,n2,n3,a)
   cdname='hail diam'
   cdunits='mm'

! 3D HYDROMETEOR TEMPERATURE, THERMAL ENERGY, LIQUID WATER FRACTION

elseif(cvar(1:lv).eq.'q2') then
   ivar_type=3
   ierr= RAMS_getvar('Q2',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname='q2'
   cdunits='J/kg'

elseif(cvar(1:lv).eq.'q6') then
   ivar_type=3
   ierr= RAMS_getvar('Q6',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname='q6'
   cdunits='J/kg'

elseif(cvar(1:lv).eq.'q7') then
   ivar_type=3
   ierr= RAMS_getvar('Q7',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname='q7'
   cdunits='J/kg'

elseif(cvar(1:lv).eq.'rain_temp') then
   ivar_type=3
   ierr= RAMS_getvar('Q2',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_raintemp(n1,n2,n3,a)
   cdname='rain temperature'
   cdunits='K'

elseif(cvar(1:lv).eq.'graup_temp') then
   ivar_type=3
   ierr= RAMS_getvar('Q6',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_qtcpcp(n1,n2,n3,a)
   cdname='graupel temperature'
   cdunits='C'

elseif(cvar(1:lv).eq.'hail_temp') then
   ivar_type=3
   ierr= RAMS_getvar('Q7',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_qtcpcp(n1,n2,n3,a)
   cdname='hail temperature'
   cdunits='C'

elseif(cvar(1:lv).eq.'rain_air_tempdif') then
   ivar_type=3
   ierr= RAMS_getvar('Q2',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_raintemp(n1,n2,n3,a)
   ierr= RAMS_getvar('THETA',idim_type,ngrd,d,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr= RAMS_getvar('PI',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_tempK(n1,n2,n3,d,c)
   call RAMS_comp_tempC(n1,n2,n3,d)
   call RAMS_comp_subt(n1,n2,n3,a,d)
   cdname='rain-air temp'
   cdunits='K'

elseif(cvar(1:lv).eq.'graup_air_tempdf') then
   ivar_type=3
   ierr= RAMS_getvar('Q6',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_qtcpcp(n1,n2,n3,a)
   ierr= RAMS_getvar('THETA',idim_type,ngrd,d,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr= RAMS_getvar('PI',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_tempK(n1,n2,n3,d,c)
   call RAMS_comp_tempC(n1,n2,n3,d)
   call RAMS_comp_subt(n1,n2,n3,a,d)
   cdname='graupel-air temp'
   cdunits='K'

elseif(cvar(1:lv).eq.'hail_air_tempdif') then
   ivar_type=3
   ierr= RAMS_getvar('Q7',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_qtcpcp(n1,n2,n3,a)
   ierr= RAMS_getvar('THETA',idim_type,ngrd,d,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr= RAMS_getvar('PI',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_tempK(n1,n2,n3,d,c)
   call RAMS_comp_tempC(n1,n2,n3,d)
   call RAMS_comp_subt(n1,n2,n3,a,d)
   cdname='hail-air temp'
   cdunits='K'

elseif(cvar(1:lv).eq.'graup_fracliq') then
   ivar_type=3
   ierr= RAMS_getvar('Q6',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_fracliq(n1,n2,n3,a)
   cdname='graupel liq frac'
   cdunits=' '

elseif(cvar(1:lv).eq.'hail_fracliq') then
   ivar_type=3
   ierr= RAMS_getvar('Q7',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_fracliq(n1,n2,n3,a)
   cdname='hail liq frac'
   cdunits=' '

! 3D MISCELLANEOUS FIELDS

!Demerval[
elseif(cvar(1:lv).eq.'pw') then
   ivar_type=2
   ierr= RAMS_getvar('RV',idim_type,ngrd,e,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr= RAMS_getvar('TOPT',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
         
   call RAMS_comp_noneg(n1,n2,n3,e)
   call RAMS_comp_pw(n1,n2,n3,a,e,c,ngrd)     
   cdname='Precipitable water'
   cdunits='cm' 

!Demerval]


elseif(cvar(1:lv).eq.'geo') then
   ivar_type=3
   ierr= RAMS_getvar('TOPT',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_z(n1,n2,n3,a,c,ngrd)
   cdname='geopotential height'
   cdunits='m'

elseif(cvar(1:lv).eq.'tke') then
   ivar_type=3
   ierr= RAMS_getvar('TKEP',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_noneg(n1,n2,n3,a)
   cdname='turb kinetic energy'
   cdunits='m2/s2'
!-------------
!trace gases

!     ; etc. for scalars

elseif(cvar(1:lv).eq.'TKUO') then
   ivar_type=3
   ierr= RAMS_getvar('DUM3',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_transf_ppb_day(n1,n2,n3,a)
   cdname='CO tend conc due conv trans'
   cdunits='ppb/day'


elseif(cvar(1:lv).eq.'thsrc') then
   ivar_type=6
   ierr= RAMS_getvar('THSRC',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_mults(n1,n2,n3*n6,a,86400.)
   call get_cumulus(n1,n2,n3,n6,a,a6)
   cdname='deep conv heat rate'
   cdunits='K/day'

elseif(cvar(1:lv).eq.'rtsrc') then
   ivar_type=6
   ierr= RAMS_getvar('RTSRC',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_mults(n1,n2,n3*n6,a,86400.)
   call RAMS_comp_mults(n1,n2,n3*n6,a,1000.)
   call get_cumulus(n1,n2,n3,n6,a,a6)
   cdname='deep conv moist rate'
   cdunits='g/kg/day'

elseif(cvar(1:lv).eq.'curidp') then
   ivar_type=3
   ierr= RAMS_getvar('d3005',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_mults(n1,n2,n3,a,86400.)
   call RAMS_comp_mults(n1,n2,n3,a,1000.)
   cdname='conv liquid/ice rate'
   cdunits='g/kg/day'

elseif(cvar(1:lv).eq.'fthrd') then
   ivar_type=3
   ierr= RAMS_getvar('FTHRD',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_mults(n1,n2,n3,a,86400.)
   cdname='rad heat rate'
   cdunits='K/day'

elseif(cvar(1:lv).eq.'fthrd_lw') then
   ivar_type=3
   ierr= RAMS_getvar('FTHRD_LW',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_mults(n1,n2,n3,a,86400.)
   cdname='rad heat rate (longwave)'
   cdunits='K/day'

elseif(cvar(1:lv).eq.'khh') then
   ivar_type=3
   ierr= RAMS_getvar('HKH',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname='horiz diffusion coeff'
   cdunits='m2/s'

elseif(cvar(1:lv).eq.'khv') then
   ivar_type=3
   ierr= RAMS_getvar('VKH',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname='vert diffusion coeff'
   cdunits='m2/s'

! 2D SURFACE PRECIPITATION

!      elseif(cvar(1:lv).eq.'accpc') then
!         ivar_type=2
!         ierr= RAMS_getvar('ACCPC',idim_type,ngrd,a,b,flnm)
!         cdname='accum fog precip'
!         cdunits='kg/m2'

elseif(cvar(1:lv).eq.'accpr' .or. cvar(1:lv).eq.'liqpcp') then
   ivar_type=2
   ierr= RAMS_getvar('ACCPR',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   if (cvar(1:lv) == 'accpr') then
      cdname='accum rain'
   else
      cdname='purely liquid precip'
   end if
   cdunits='kg/m²'

elseif(cvar(1:lv).eq.'accpp') then
   ivar_type=2
   ierr= RAMS_getvar('ACCPP',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname='accum pristine'
   cdunits='kg/m2'

elseif(cvar(1:lv).eq.'accps') then
   ivar_type=2
   ierr= RAMS_getvar('ACCPS',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname='accum snow'
   cdunits='kg/m2'

elseif(cvar(1:lv).eq.'accpa') then
   ivar_type=2
   ierr= RAMS_getvar('ACCPA',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname='accum aggregates'
   cdunits='kg/m2'

elseif(cvar(1:lv).eq.'accpg') then
   ivar_type=2
   ierr= RAMS_getvar('ACCPG',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname='accum graupel'
   cdunits='kg/m2'

elseif(cvar(1:lv).eq.'accph') then
   ivar_type=2
   ierr= RAMS_getvar('ACCPH',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname='accum hail'
   cdunits='kg/m2'

elseif(cvar(1:lv).eq.'totpcp' .or. cvar(1:lv).eq.'totpcp_in' .or.  &
       cvar(1:lv).eq.'precip' .or. cvar(1:lv).eq.'precip_in') then
   ivar_type=2
   call RAMS_comp_zero(n1,n2,1,a)
   ierr= RAMS_getvar('ACCPR',idim_type,ngrd,c,b,flnm)
   if(ierr.eq.0) call RAMS_comp_accum(n1,n2,1,a,c)
   ierr= RAMS_getvar('ACCPP',idim_type,ngrd,c,b,flnm)
   if(ierr.eq.0) call RAMS_comp_accum(n1,n2,1,a,c)
   ierr= RAMS_getvar('ACCPS',idim_type,ngrd,c,b,flnm)
   if(ierr.eq.0) call RAMS_comp_accum(n1,n2,1,a,c)
   ierr= RAMS_getvar('ACCPA',idim_type,ngrd,c,b,flnm)
   if(ierr.eq.0) call RAMS_comp_accum(n1,n2,1,a,c)
   ierr= RAMS_getvar('ACCPG',idim_type,ngrd,c,b,flnm)
   if(ierr.eq.0) call RAMS_comp_accum(n1,n2,1,a,c)
   ierr= RAMS_getvar('ACCPH',idim_type,ngrd,c,b,flnm)
   if(ierr.eq.0) call RAMS_comp_accum(n1,n2,1,a,c)

   if (cvar(1:lv).eq.'precip'.or.cvar(1:lv).eq.'precip_in') then
      ierr= RAMS_getvar('ACONPR',idim_type,ngrd,c,b,flnm)
      if(ierr.eq.0) call RAMS_comp_accum(n1,n2,1,a,c)
      cdname='total accum precip'
   else
      cdname='total resolved precip'
   endif

   if(cvar(1:lv).eq.'totpcp'.or.cvar(1:lv).eq.'precip') then
      cdunits='kg/m²'
   else
      call RAMS_comp_mults(n1,n2,n3,a,.03937)
      cdunits='in liq'
   endif
   call RAMS_comp_noneg(n1,n2,1,a)

elseif(cvar(1:lv).eq.'icepcp') then
   ivar_type=2
   call RAMS_comp_zero(n1,n2,1,a)
   ierr= RAMS_getvar('ACCPP',idim_type,ngrd,c,b,flnm)
   if(ierr.eq.0) call RAMS_comp_accum(n1,n2,1,a,c)
   ierr= RAMS_getvar('ACCPS',idim_type,ngrd,c,b,flnm)
   if(ierr.eq.0) call RAMS_comp_accum(n1,n2,1,a,c)
   ierr= RAMS_getvar('ACCPA',idim_type,ngrd,c,b,flnm)
   if(ierr.eq.0) call RAMS_comp_accum(n1,n2,1,a,c)

   cdname='purely ice precip'
   cdunits='kg/m²'
   call RAMS_comp_noneg(n1,n2,1,a)

elseif(cvar(1:lv).eq.'mixpcp') then
   ivar_type=2
   call RAMS_comp_zero(n1,n2,1,a)
   ierr= RAMS_getvar('ACCPG',idim_type,ngrd,c,b,flnm)
   if(ierr.eq.0) call RAMS_comp_accum(n1,n2,1,a,c)
   ierr= RAMS_getvar('ACCPH',idim_type,ngrd,c,b,flnm)
   if(ierr.eq.0) call RAMS_comp_accum(n1,n2,1,a,c)

   cdname='mixed (ice/liq) precip'
   cdunits='kg/m²'
   call RAMS_comp_noneg(n1,n2,1,a)

elseif(cvar(1:lv).eq.'pcprr') then
   ivar_type=2
   ierr= RAMS_getvar('PCPRR',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_mults(n1,n2,1,a,3600.)
   cdname='rain precip rate'
   cdunits='mm/hr liq equiv'

elseif(cvar(1:lv).eq.'pcprp') then
   ivar_type=2
   ierr= RAMS_getvar('PCPRP',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_mults(n1,n2,1,a,3600.)
   cdname='pristine precip rate'
   cdunits='mm/hr liq equiv'

elseif(cvar(1:lv).eq.'psprs') then
   ivar_type=2
   ierr= RAMS_getvar('PCPRS',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_mults(n1,n2,1,a,3600.)
   cdname='snow precip rate'
   cdunits='mm/hr liq equiv'

elseif(cvar(1:lv).eq.'pcpra') then
   ivar_type=2
   ierr= RAMS_getvar('PCPRA',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_mults(n1,n2,1,a,3600.)
   cdname='aggregates precip rate'
   cdunits='mm/hr liq equiv'

elseif(cvar(1:lv).eq.'pcprg') then
   ivar_type=2
   ierr= RAMS_getvar('PCPRG',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname='graupel precip rate'
   cdunits='mm/hr liq equiv'

elseif(cvar(1:lv).eq.'pcprh') then
   ivar_type=2
   ierr= RAMS_getvar('PCPRH',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_mults(n1,n2,1,a,3600.)
   cdname='hail precip rate'
   cdunits='mm/hr liq equiv'

elseif(cvar(1:lv).eq.'pcpg') then
   ivar_type=2
   ierr= RAMS_getvar('PCPG',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname='pcpg'
   cdunits='kg/m2'

elseif(cvar(1:lv).eq.'qpcpg') then
   ivar_type=2
   ierr= RAMS_getvar('QPCPG',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname='qpcpg'
   cdunits='J/m2'

elseif(cvar(1:lv).eq.'dpcpg') then
   ivar_type=2
   ierr= RAMS_getvar('DPCPG',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname='dpdpg'
   cdunits='m'

elseif(cvar(1:lv).eq.'pcprate'.or.cvar(1:lv).eq.'pcprate_in') then
   ivar_type=2
   call RAMS_comp_zero(n1,n2,1,a)
   ierr= RAMS_getvar('PCPRR',idim_type,ngrd,c,b,flnm)
   if(ierr.eq.0) call RAMS_comp_accum(n1,n2,1,a,c)
   ierr= RAMS_getvar('PCPRP',idim_type,ngrd,c,b,flnm)
   if(ierr.eq.0) call RAMS_comp_accum(n1,n2,1,a,c)
   ierr= RAMS_getvar('PCPRS',idim_type,ngrd,c,b,flnm)
   if(ierr.eq.0) call RAMS_comp_accum(n1,n2,1,a,c)
   ierr= RAMS_getvar('PCPRA',idim_type,ngrd,c,b,flnm)
   if(ierr.eq.0) call RAMS_comp_accum(n1,n2,1,a,c)
   ierr= RAMS_getvar('PCPRG',idim_type,ngrd,c,b,flnm)
   if(ierr.eq.0) call RAMS_comp_accum(n1,n2,1,a,c)
   ierr= RAMS_getvar('PCPRH',idim_type,ngrd,c,b,flnm)
   if(ierr.eq.0) call RAMS_comp_accum(n1,n2,1,a,c)
   call RAMS_comp_noneg(n1,n2,1,a)
   cdname='resolved precip rate'

   if(cvar(1:lv).eq.'pcprate'.or.cvar(1:lv).eq.'precipr') then
      call RAMS_comp_mults(n1,n2,1,a,3600.)
      cdunits='mm/hr'
   elseif(cvar(1:lv).eq.'pcprate_in'.or.cvar(1:lv).eq.'precipr_in') then
      call RAMS_comp_mults(n1,n2,1,a,141.732)
      cdunits='in/hr'
   endif

elseif(cvar(1:lv).eq.'conpcp' .or. cvar(1:lv).eq.'conprr') then
   ivar_type=9
   ierr= RAMS_getvar('CONPRR',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_mults(n1,n2,1,a,3600.)
   call RAMS_comp_noneg(n1,n2,1,a)
   cdname='convective pcp rate'
   cdunits='mm/hr'

elseif(cvar(1:lv).eq.'acccon') then
   ivar_type=2
   ierr= RAMS_getvar('ACONPR',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_noneg(n1,n2,1,a)
   cdname='accum convective pcp'
   cdunits='mm'




elseif(cvar(1:lv).eq.'cape') then
   ivar_type=2
!- rel hum (e)
   ierr= RAMS_getvar('RV',idim_type,ngrd,e,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr= RAMS_getvar('PI',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr= RAMS_getvar('THETA',idim_type,ngrd,d,b,flnm)
   ierr_getvar = ierr_getvar + ierr

   call RAMS_comp_rh(n1,n2,n3,e,c,d)
   call RAMS_comp_noneg(n1,n2,n3,e)
!- tempk (d)
   ierr= RAMS_getvar('THETA',idim_type,ngrd,d,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr= RAMS_getvar('PI',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_tempK(n1,n2,n3,d,c)
!- press (c)
   ierr= RAMS_getvar('PI',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_press(n1,n2,n3,c)
!- cape
   call cape_cine(n1,n2,n3,c,d,e,a,'cape',-9.99e33)

   cdname='cape'
   cdunits='J/kg'

elseif(cvar(1:lv).eq.'cine') then
   ivar_type=2
!- rel hum (e)
   ierr= RAMS_getvar('RV',idim_type,ngrd,e,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr= RAMS_getvar('PI',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr= RAMS_getvar('THETA',idim_type,ngrd,d,b,flnm)
   ierr_getvar = ierr_getvar + ierr

   call RAMS_comp_rh(n1,n2,n3,e,c,d)
   call RAMS_comp_noneg(n1,n2,n3,e)
!- tempk (d)
   ierr= RAMS_getvar('THETA',idim_type,ngrd,d,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr= RAMS_getvar('PI',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_tempK(n1,n2,n3,d,c)
!- press (c)
   ierr= RAMS_getvar('PI',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_press(n1,n2,n3,c)
!- cape
   call cape_cine(n1,n2,n3,c,d,e,a,'cine',-9.99e33)

   cdname='cine'
   cdunits='J/kg'

! Tropopause values
![ML - Altura da tropopausa
elseif(cvar(1:lv).eq.'ztropop') then
   ivar_type=2
   ierr=RAMS_getvar('THETA',idim_type,ngrd,e,b,flnm)  ! e= temperatura potencial
   ierr_getvar = ierr_getvar + ierr
   call acha_ztropop(n1,n2,n3,a,e,ngrd)               ! a= altura da tropopausa
   cdname='Tropopause height'
   cdunits='m - sigmaz'
!ML]   

![ML - Temperatura da tropopausa
elseif(cvar(1:lv).eq.'ttropop') then
   ivar_type=2
   ierr=RAMS_getvar('THETA',idim_type,ngrd,c,b,flnm)  ! c= temperatura potencial
   ierr_getvar = ierr_getvar + ierr
   ierr=RAMS_getvar('THETA',idim_type,ngrd,e,b,flnm)  ! e= temperatura potencial
   ierr_getvar = ierr_getvar + ierr
   ierr=RAMS_getvar('PI',idim_type,ngrd,d,b,flnm)     ! d= funçao de Exner
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_tempK(n1,n2,n3,e,d)                 ! e= temperatura em Kelvin
   call acha_ttropop(n1,n2,n3,a,c,e,ngrd)             ! a= temperatura da tropopausa
   cdname='Tropopause temperature'
   cdunits='K'
!ML]   

![ML - Pressao na tropopausa
elseif(cvar(1:lv).eq.'ptropop') then
   ivar_type=2
   ierr=RAMS_getvar('THETA',idim_type,ngrd,c,b,flnm)  ! c= temperatura potencial
   ierr_getvar = ierr_getvar + ierr
   ierr=RAMS_getvar('PI',idim_type,ngrd,e,b,flnm)     ! e= funçao de Exner
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_press(n1,n2,n3,e)                   ! e= pressao
   call acha_ptropop(n1,n2,n3,a,c,e,ngrd)             ! a= temperatura da tropopausa
   cdname='Tropopause pressure'
   cdunits='hPa'
!ML]   


![Marcos Parâmetros da convecção para uso em STILT

elseif(cvar(1:lv).eq.'afxu') then
   ivar_type=3
   ierr= RAMS_getvar('AFXUB',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname='Mean area zonal flux'
   cdunits='kg/m²/s'

elseif(cvar(1:lv).eq.'afxv') then
   ivar_type=3
   ierr= RAMS_getvar('AFXVB',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname='Mean area meridional flux'
   cdunits='kg/m²/s'

elseif(cvar(1:lv).eq.'afxw') then
   ivar_type=3
   ierr= RAMS_getvar('AFXWB',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname='Mean area vertical flux'
   cdunits='kg/m²/s'

elseif(cvar(1:lv).eq.'cfxup_deep') then
   ivar_type=3
   ierr= RAMS_getvar('CFXUP1',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
!   call RAMS_comp_noneg(n1,n2,n3,a)
   cdname='Conv. upward flux - deep'
   cdunits='kg/m²/s'

elseif(cvar(1:lv).eq.'cfxdn_deep') then
   ivar_type=3
   ierr= RAMS_getvar('CFXDN1',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
!   call RAMS_comp_nopos(n1,n2,n3,a)
   cdname='Conv. downward flux - deep'
   cdunits='kg/m²/s'

elseif(cvar(1:lv).eq.'cfxup_shal') then
   ivar_type=3
   ierr= RAMS_getvar('CFXUP2',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
!   call RAMS_comp_noneg(n1,n2,n3,a)
   cdname='Conv. upward flux - shallow'
   cdunits='kg/m²/s'

elseif(cvar(1:lv).eq.'efxup_deep') then
   ivar_type=3
   ierr= RAMS_getvar('EFXUP1',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
!   call RAMS_comp_noneg(n1,n2,n3,a)
   cdname='Updraft entrainment flux - deep'
   cdunits='kg/m²/s'

elseif(cvar(1:lv).eq.'efxdn_deep') then
   ivar_type=3
   ierr= RAMS_getvar('EFXDN1',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
!   call RAMS_comp_noneg(n1,n2,n3,a)
   cdname='Downdraft entrainment flux - deep'
   cdunits='kg/m²/s'

elseif(cvar(1:lv).eq.'efxup_shal') then
   ivar_type=3
   ierr= RAMS_getvar('EFXUP2',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
!   call RAMS_comp_noneg(n1,n2,n3,a)
   cdname='Updraft entrainment flux - shallow'
   cdunits='kg/m²/s'

elseif(cvar(1:lv).eq.'dfxup_deep') then
   ivar_type=3
   ierr= RAMS_getvar('DFXUP1',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
!   call RAMS_comp_noneg(n1,n2,n3,a)
   cdname='Updraft detrainment flux - deep'
   cdunits='kg/m²/s'

elseif(cvar(1:lv).eq.'dfxdn_deep') then
   ivar_type=3
   ierr= RAMS_getvar('DFXDN1',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
!   call RAMS_comp_noneg(n1,n2,n3,a)
   cdname='Downdraft detrainment flux - deep'
   cdunits='kg/m²/s'

elseif(cvar(1:lv).eq.'dfxup_shal') then
   ivar_type=3
   ierr= RAMS_getvar('DFXUP2',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
!   call RAMS_comp_noneg(n1,n2,n3,a)
   cdname='Updraft detrainment flux - shallow'
   cdunits='kg/m²/s'

elseif(cvar(1:lv).eq.'cfxup') then
   ivar_type=6
   ierr= RAMS_getvar('CFXUP',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_noneg(n1,n2,n3*n6,a)
   call get_cumulus(n1,n2,n3,n6,a,a6)
   cdname='Convective upward flux'
   cdunits='kg/m²/s'

elseif(cvar(1:lv).eq.'cfxdn') then
   ivar_type=6
   ierr= RAMS_getvar('CFXDN',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_nopos(n1,n2,n3*n6,a)
   call get_cumulus(n1,n2,n3,n6,a,a6)
   cdname='Convective downward flux'
   cdunits='kg/m²/s'

elseif(cvar(1:lv).eq.'dfxup') then
   ivar_type=6
   ierr= RAMS_getvar('DFXUP',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_noneg(n1,n2,n3*n6,a)
   call get_cumulus(n1,n2,n3,n6,a,a6)
   cdname='Detrainment upward flux'
   cdunits='kg/m²/s'

elseif(cvar(1:lv).eq.'dfxdn') then
   ivar_type=6
   ierr= RAMS_getvar('DFXDN',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_noneg(n1,n2,n3*n6,a)
   call get_cumulus(n1,n2,n3,n6,a,a6)
   cdname='Detrainment upward flux'
   cdunits='kg/m²/s'

elseif(cvar(1:lv).eq.'efxup') then
   ivar_type=6
   ierr= RAMS_getvar('EFXUP',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_noneg(n1,n2,n3*n6,a)
   call get_cumulus(n1,n2,n3,n6,a,a6)
   cdname='Entrainment upward flux'
   cdunits='kg/m²/s'

elseif(cvar(1:lv).eq.'efxdn') then
   ivar_type=6
   ierr= RAMS_getvar('EFXDN',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_noneg(n1,n2,n3*n6,a)
   call get_cumulus(n1,n2,n3,n6,a,a6)
   cdname='Entrainment upward flux'
   cdunits='kg/m²/s'

! Extra turbulence parameters

elseif(cvar(1:lv).eq.'tkem') then
   ivar_type=3
   ierr= RAMS_getvar('TKEPB',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
!   call RAMS_comp_noneg(n1,n2,n3,a)
   cdname='Mean turbulent kinetic energy'
   cdunits='m²/s²'

elseif(cvar(1:lv).eq.'ltscale') then
   ivar_type=3
   ierr= RAMS_getvar('TL',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
!   call RAMS_comp_noneg(n1,n2,n3,a)
   cdname='Mean Lagrangean time scale'
   cdunits='s'

elseif(cvar(1:lv).eq.'sigw') then
   ivar_type=3
   ierr= RAMS_getvar('SIGW',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
!   call RAMS_comp_noneg(n1,n2,n3,a)
   cdname='Mean vertical velocity standard deviation'
   cdunits='m/s'

elseif(cvar(1:lv).eq.'sm') then
   ivar_type=3
   ierr= RAMS_getvar('SSM',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
!   call RAMS_comp_noneg(n1,n2,n3,a)
   cdname='Eddy diffusity coefficient for momentum'
   cdunits='n/d'

elseif(cvar(1:lv).eq.'sh') then
   ivar_type=3
   ierr= RAMS_getvar('SSH',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
!   call RAMS_comp_noneg(n1,n2,n3,a)
   cdname='Eddy diffusity coefficient for heat'
   cdunits='n/d'

elseif(cvar(1:lv).eq.'lturb') then
   ivar_type=3
   ierr= RAMS_getvar('LTURB',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
!   call RAMS_comp_noneg(n1,n2,n3,a)
   cdname='Turbulence scale'
   cdunits='m'

elseif(cvar(1:lv).eq.'rich') then
   ivar_type=3
   ierr= RAMS_getvar('RICH',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
!   call RAMS_comp_noneg(n1,n2,n3,a)
   cdname='Richardson number'
   cdunits='n/d'

elseif(cvar(1:lv).eq.'rf') then
   ivar_type=3
   ierr= RAMS_getvar('RIFLUX',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
!   call RAMS_comp_noneg(n1,n2,n3,a)
   cdname='Flux Richardson number'
   cdunits='n/d'

elseif(cvar(1:lv).eq.'gm') then
   ivar_type=3
   ierr= RAMS_getvar('GGM',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
!   call RAMS_comp_noneg(n1,n2,n3,a)
   cdname='Dimensionless square of mean shear'
   cdunits='n/d'

elseif(cvar(1:lv).eq.'gh') then
   ivar_type=3
   ierr= RAMS_getvar('GGH',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
!   call RAMS_comp_noneg(n1,n2,n3,a)
   cdname='Dimensionless of neg Brunt-Vaisala frequency'
   cdunits='n/d'

elseif(cvar(1:lv).eq.'snowdepth') then
    ivar_type=2
    ierr= RAMS_getvar('SNOW_DEPTH',idim_type,ngrd,a,b,flnm)
    ierr_getvar = ierr_getvar + ierr
    cdname='Depth of the snow layer'
    cdunits='m'

elseif(cvar(1:lv).eq.'snowdepth') then
    ivar_type=2
    ierr= RAMS_getvar('SNOW_DEPTH',idim_type,ngrd,a,b,flnm)
    ierr_getvar = ierr_getvar + ierr
    cdname='Depth of the snow layer'
    cdunits='m'

elseif(cvar(1:lv).eq.'snowdepth') then
    ivar_type=2
    ierr= RAMS_getvar('SNOW_DEPTH',idim_type,ngrd,a,b,flnm)
    ierr_getvar = ierr_getvar + ierr
    cdname='Depth of the snow layer'
    cdunits='m'

elseif(cvar(1:lv).eq.'pblhgt') then
    ivar_type=2
    ierr= RAMS_getvar('PBLHGT',idim_type,ngrd,a,b,flnm)
    ierr_getvar = ierr_getvar + ierr
    cdname='PBL Depth'
    cdunits='m'

elseif(cvar(1:lv).eq.'lmo') then
    ivar_type=2
    ierr= RAMS_getvar('LMO',idim_type,ngrd,a,b,flnm)
    ierr_getvar = ierr_getvar + ierr
    cdname='Obukhov lenght scale'
    cdunits='m'

elseif(cvar(1:lv).eq.'akscal') then
    ivar_type=2
    ierr= RAMS_getvar('AKSCAL',idim_type,ngrd,a,b,flnm)
    ierr_getvar = ierr_getvar + ierr
    cdname='Scaling factor'
    cdunits=''

elseif(cvar(1:lv).eq.'m2') then
    ivar_type=3
    ierr= RAMS_getvar('MSQU',idim_type,ngrd,a,b,flnm)
    ierr_getvar = ierr_getvar + ierr
    cdname='Square of momentum frequency'
    cdunits='Hz²'

elseif(cvar(1:lv).eq.'n2') then
    ivar_type=3
    ierr= RAMS_getvar('BRUNT',idim_type,ngrd,a,b,flnm)
    ierr_getvar = ierr_getvar + ierr
    cdname='Square of Brunt-Vaisala frequency'
    cdunits='Hz²'

elseif(cvar(1:lv).eq.'tke2') then
    ivar_type=3
    ierr= RAMS_getvar('TKE2',idim_type,ngrd,a,b,flnm)
    ierr_getvar = ierr_getvar + ierr
    cdname='Level 2 TKE'
    cdunits='m²/s²'

!Marcos]


! Vertically-integrated atmospheric moisture

elseif(cvar(1:lv).eq.'vertint_rt' .or. cvar(1:lv).eq.'vertint_cond') then
   ivar_type=2

   ierr= RAMS_getvar('TOPT',idim_type,ngrd,e,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_dn0(n1,n2,n3,c,b,d,e,ngrd)

   if (cvar(1:lv).eq.'vertint_rt') then
      ierr= RAMS_getvar('RV',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='vertint total water'
   else
      call RAMS_comp_zero(n1,n2,n3,a)
      cdname='vertint condensate'
   endif

   ierr= RAMS_getvar('RCP',idim_type,ngrd,c,b,flnm)
   if(ierr.eq.0) call RAMS_comp_accum(n1,n2,n3,a,c)
   ierr= RAMS_getvar('RRP',idim_type,ngrd,c,b,flnm)
   if(ierr.eq.0) call RAMS_comp_accum(n1,n2,n3,a,c)
   ierr= RAMS_getvar('RPP',idim_type,ngrd,c,b,flnm)
   if(ierr.eq.0) call RAMS_comp_accum(n1,n2,n3,a,c)
   ierr= RAMS_getvar('RSP',idim_type,ngrd,c,b,flnm)
   if(ierr.eq.0) call RAMS_comp_accum(n1,n2,n3,a,c)
   ierr= RAMS_getvar('RAP',idim_type,ngrd,c,b,flnm)
   if(ierr.eq.0) call RAMS_comp_accum(n1,n2,n3,a,c)
   ierr= RAMS_getvar('RGP',idim_type,ngrd,c,b,flnm)
   if(ierr.eq.0) call RAMS_comp_accum(n1,n2,n3,a,c)
   ierr= RAMS_getvar('RHP',idim_type,ngrd,c,b,flnm)
   if(ierr.eq.0) call RAMS_comp_accum(n1,n2,n3,a,c)

   call RAMS_comp_mult(n1,n2,n3,a,d)
   call RAMS_comp_vertint(n1,n2,n3,a,e,ngrd)

   cdunits='mm'


! 2D SURFACE HEAT, MOISTURE, MOMENTUM AND RADIATIVE FLUXES

elseif(cvar(1:lv).eq.'SFLUX_T') then
   ivar_type=2
   ierr= RAMS_getvar('SFLUX_T',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname='SFLUX_T'
   cdunits='m'

elseif(cvar(1:lv).eq.'SFLUX_R') then
   ivar_type=2
   ierr= RAMS_getvar('SFLUX_R',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname='SFLUX_R'
   cdunits='m'

elseif(cvar(1:lv).eq.'uw') then
   ivar_type=2
   ierr= RAMS_getvar('UW',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname='uw'
   cdunits='m'

elseif(cvar(1:lv).eq.'vw') then
   ivar_type=2
   ierr= RAMS_getvar('VW',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname='vw'
   cdunits='m'

elseif(cvar(1:lv).eq.'SFLUX_W') then
   ivar_type=2
   ierr= RAMS_getvar('SFLUX_W',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname='SFLUX_W'
   cdunits='m'

elseif(cvar(1:lv).eq.'SFLUX_C') then
   ivar_type=2
   ierr= RAMS_getvar('SFLUX_C',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname='SFLUX_C'
   cdunits='m'

elseif(cvar(1:lv).eq.'hflxca') then
   ivar_type=2
   ierr= RAMS_getvar('SFLUX_T',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr= RAMS_getvar('TOPT',idim_type,ngrd,e,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_dn0(n1,n2,n3,b,c,d,e,ngrd)
   call RAMS_comp_mult(n1,n2,1,a,d)
   call RAMS_comp_mults(n1,n2,1,a,1004.)
   cdname='sfc sens heat flx canopy to atmosphere'
   cdunits='W/m2'

elseif(cvar(1:lv).eq.'co2flux') then
   ivar_type=2
   ierr= RAMS_getvar('SFLUX_C',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr= RAMS_getvar('TOPT',idim_type,ngrd,e,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_dn0(n1,n2,n3,b,c,d,e,ngrd)
   call RAMS_comp_mult(n1,n2,1,a,d)
   call RAMS_comp_mults(n1,n2,1,a,mmdryi)
   cdname='CO2 flux'
   cdunits='umol/m2/s'

elseif(cvar(1:lv).eq.'lfsensi' .or. cvar(1:lv).eq.'lfsensi_ps') then

   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
!!
   ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
        ,a(irecind),b,flnm)

   irecind = irecind + irecsize
!!
   ierr = RAMS_getvar('LF_SENSI',idim_type,ngrd   &
        ,a(irecind),b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_mults(n1,n2,1,a(irecind),1004.)

   if(cvar(1:lv).eq.'lfsensi') then
      ivar_type = 7
   else
      ivar_type = 2
      call RAMS_comp_patchsum(nnxp(ngrd),nnyp(ngrd),1,npatch  &
         ,a(irecind),a(1),b)
   endif

   cdname=' Sensible heat'
   cdunits='W/m²'


elseif(cvar(1:lv).eq.'lfevap' .or. cvar(1:lv).eq.'lfevap_ps') then

   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   if (cvar(1:lv).eq.'lfevap_ps') then
      ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
           ,a(irecind),b,flnm)

      irecind = irecind + irecsize
   end if
   ierr = RAMS_getvar('LF_EVAP',idim_type,ngrd   &
        ,a(irecind),b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_mults(n1,n2,1,a(irecind),2.5e6)

   if(cvar(1:lv).eq.'lfevap') then
      ivar_type = 7
   else
      ivar_type = 2
      call RAMS_comp_patchsum(nnxp(ngrd),nnyp(ngrd),1,npatch  &
         ,a(irecind),a(1),b)
   endif

   cdname=' Evaporation'
   cdunits='W/m²'

elseif(cvar(1:lv).eq.'lftransp' .or. cvar(1:lv).eq.'lftransp_ps') then

   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   if(cvar(1:lv).eq.'lftransp_ps') then
      ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
           ,a(irecind),b,flnm)

      irecind = irecind + irecsize
   end if
   ierr = RAMS_getvar('LF_TRANSP',idim_type,ngrd   &
        ,a(irecind),b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_mults(n1,n2,1,a(irecind),2.5e6)

   if(cvar(1:lv).eq.'lftransp') then
      ivar_type = 7
   else
      ivar_type = 2
      call RAMS_comp_patchsum(nnxp(ngrd),nnyp(ngrd),1,npatch  &
         ,a(irecind),a(1),b)
   endif

   cdname=' Transpiration'
   cdunits='W/m²'

elseif(cvar(1:lv).eq.'qwflxca') then
   ivar_type=2
   ierr= RAMS_getvar('SFLUX_R',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr= RAMS_getvar('TOPT',idim_type,ngrd,e,b,flnm)
   call RAMS_comp_dn0(n1,n2,n3,b,c,d,e,ngrd)
   call RAMS_comp_mult(n1,n2,1,a,d)
   call RAMS_comp_mults(n1,n2,1,a,2.5e6)
   cdname='water flux from canopy to atmosphere'
   cdunits='W/m2'

elseif(cvar(1:lv).eq.'etrans') then
   ivar_type=2
   ierr= RAMS_getvar('SFLUX_R',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr= RAMS_getvar('TOPT',idim_type,ngrd,e,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_dn0(n1,n2,n3,b,c,d,e,ngrd)
   call RAMS_comp_mult(n1,n2,1,a,d)
!                 Divide by water density to get depth and 
!                   convert units from m/s to mm/hour (3600./1000.)
   call RAMS_comp_mults(n1,n2,1,a,3.6)
   cdname='evapo-transpiration'
   cdunits='mm/hour'

elseif(cvar(1:lv).eq.'etrans_in') then
   ivar_type=2
   ierr= RAMS_getvar('SFLUX_R',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr= RAMS_getvar('TOPT',idim_type,ngrd,e,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_dn0(n1,n2,n3,b,c,d,e,ngrd)
   call RAMS_comp_mult(n1,n2,1,a,d)
!                 Divide by water density to get depth and 
!                   convert units from m/s to in/hour (39.37 * 3600./1000.)
   call RAMS_comp_mults(n1,n2,n3,a,141.732)
   cdname='evapo-transpiration'
   cdunits='in/hour'

elseif(cvar(1:lv).eq.'umom_flx') then
   ivar_type=2
   ierr= RAMS_getvar('UW',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr= RAMS_getvar('TOPT',idim_type,ngrd,e,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_dn0(n1,n2,n3,b,c,d,e,ngrd)
   call RAMS_comp_mult(n1,n2,1,a,d)
   cdname='sfc u-momentum flx'
   cdunits='Pa'

elseif(cvar(1:lv).eq.'vmom_flx') then
   ivar_type=2
   ierr= RAMS_getvar('VW',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr= RAMS_getvar('TOPT',idim_type,ngrd,e,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_dn0(n1,n2,n3,b,c,d,e,ngrd)
   call RAMS_comp_mult(n1,n2,1,a,d)
   cdname='sfc v-momentum flx'
   cdunits='Pa'

elseif(cvar(1:lv).eq.'wmom_flx') then
   ivar_type=2
   ierr= RAMS_getvar('SFLUX_W',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr= RAMS_getvar('TOPT',idim_type,ngrd,e,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_dn0(n1,n2,n3,b,c,d,e,ngrd)
   call RAMS_comp_mult(n1,n2,1,a,d)
   cdname='sfc w-momentum flx'
   cdunits='Pa'

elseif(cvar(1:lv).eq.'bowen') then
   ivar_type=2
   ierr= RAMS_getvar('SFLUX_T',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr= RAMS_getvar('SFLUX_R',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_bowen(n1,n2,1,a,c)
   cdname='bowen ratio'
   cdunits=' '

elseif(cvar(1:lv).eq.'cosz') then
   ivar_type=2
   ierr= RAMS_getvar('COSZ',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname='co-sine of zenith angle'
   cdunits=' '

elseif(cvar(1:lv).eq.'rshort') then
   ivar_type=2
   ierr= RAMS_getvar('RSHORT',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname='rshort'
   cdunits='W/m2'

elseif(cvar(1:lv).eq.'rlong') then
   ivar_type=2
   ierr= RAMS_getvar('RLONG',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname='rlong'
   cdunits='W/m2'

elseif(cvar(1:lv).eq.'rlongup') then
   ivar_type=2
   ierr= RAMS_getvar('RLONGUP',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname='rlongup'
   cdunits='W/m2'

elseif(cvar(1:lv).eq.'albedt') then
   ivar_type=2
   ierr= RAMS_getvar('ALBEDT',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname='albedt'
   cdunits=' '

!2D misc
elseif(cvar(1:lv).eq.'qsc1') then
   ivar_type=2
   ierr= RAMS_getvar('DUM1',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call D3toD2(n1,n2,n3,2,a,c)
   cdname='qsc1'
   cdunits='????'



! 2D TOPOGRAPHY AND GEOGRAPHIC VALUES

elseif(cvar(1:lv).eq.'topoa') then
   ivar_type=2
   ierr= RAMS_getvar('TOPTA',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname='topo'
   cdunits='m'

elseif(cvar(1:lv).eq.'topo') then
   ivar_type=2
   ierr= RAMS_getvar('TOPT',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname='topo'
   cdunits='m'

elseif(cvar(1:lv).eq.'lat') then
   ivar_type=2
   ierr= RAMS_getvar('GLAT',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname='latitude'
   cdunits='deg'

elseif(cvar(1:lv).eq.'lon') then
   ivar_type=2
   ierr= RAMS_getvar('GLON',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname='longitude'
   cdunits='deg'

elseif(cvar(1:lv).eq.'mynum') then
   ivar_type=2
   ierr= RAMS_getvar('MYNUM',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname='node ID'
   cdunits=' '

! 2D MISCELLANEOUS FIELDS

elseif(cvar(1:lv).eq.'slp_rams') then
   ivar_type=2
   ierr= RAMS_getvar('TOPT',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_z(n1,n2,n3,c,a,ngrd)

   ierr= RAMS_getvar('PI',idim_type,ngrd,d,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr= RAMS_getvar('THETA',idim_type,ngrd,e,b,flnm)
   ierr_getvar = ierr_getvar + ierr

   call RAMS_comp_slpress(n1,n2,n3,e,d,c,a)
   cdname='sea level pressure'
   cdunits='mb'

elseif(cvar(1:lv).eq.'sea_press') then
   ivar_type=2
   ierr= RAMS_getvar('TOPT',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr= RAMS_getvar('PI',idim_type,ngrd,d,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr= RAMS_getvar('THETA',idim_type,ngrd,e,b,flnm)
   ierr_getvar = ierr_getvar + ierr

! get terrain in c(1+ioffc)
!         ierr= RAMS_getvar(1,ind,nind,2,ngrd,c(1+ioffc),b,flnm)
! get Exner function in d(1+ioffd)
!         ierr= RAMS_getvar(7,ind,nind,3,ngrd,d(1+ioffd),b,flnm)
! get theta in 1+ioffe
!         ierr= RAMS_getvar(32,ind,nind,3,ngrd,e(1+ioffe),b,flnm)

         call RAMS_comp_slpmm5(n1,n2,n3,e,d,c,a)
         cdname='sea level pressure;'
         cdunits='mb;'



elseif(cvar(1:lv).eq.'sfc_div') then
   ivar_type=2
   ierr= RAMS_getvar('WP',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_sfcdiv(n1,n2,n3,a,ngrd)
   cdname='surface divergence'
   cdunits='1/s'

! Special use of sst: acquired for patch #1 even where no water exists

elseif(cvar(1:lv).eq.'sst') then
   ivar_type=2
   ierr = RAMS_getvar('SOIL_ENERGY',idim_type,ngrd   &
        ,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   kp = nzg
   call rams_fill_sst(n1,n2,nzg*npatch,kp,a,c)

!   call RAMS_comp_tempC(n1,n2,1,a)
   cdname='water temperature'
   cdunits='C'


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! LEAF2 variables section

! If want a horiz plot, specify a string like 'tgpatch'; it will
!   return i,j,ip array.
! Specify a new ivar_type, not corresponding to anal file var type.  With
!   horiz plot, get back into iplt.  If have this var type, don't slice.
! Need replacement for rams3to2d because windowing is done in there.
! Replacement would window but not slice.
! Then, if want xz (vert cross section) have name like tgpatch_vert.
! This would return entire 4d array from hvlib.f.
! Then we have to slice and window with yet another replacement to rams3to2d.

! nkk is the record number, where n is the LEAF field number (1, 2, 3, or 4)
! and kk is the k level.
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


elseif(cvar(1:lv).eq.'pfarea') then

   ivar_type = 7
   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
        ,a(irecind),b,flnm)
   ierr_getvar = ierr_getvar + ierr

   irecind = irecind + irecsize
   ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
        ,a(irecind),b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname='patch fractional area'
   cdunits=''

elseif(cvar(1:lv).eq.'soil_z0_p' .or. cvar(1:lv).eq.'soil_z0_ps') then

   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   if(cvar(1:lv).eq.'soil_z0_ps') then
      ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
           ,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr
      irecind = irecind + irecsize
   end if
   ierr = RAMS_getvar('SOIL_Z0',idim_type,ngrd   &
        ,a(irecind),b,flnm)
   ierr_getvar = ierr_getvar + ierr

   if(cvar(1:lv).eq.'soil_z0_p') then
      ivar_type = 7
   else
      ivar_type = 2
      call RAMS_comp_patchsum_l(nnxp(ngrd),nnyp(ngrd),1,npatch  &
         ,a(irecind),a(1),b)
   endif

   cdname='soil roughness'
   cdunits='m'

elseif(cvar(1:lv).eq.'vtype' .or. cvar(1:lv).eq.'vtype_bp') then

   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch

   if(cvar(1:lv).eq.'vtype_bp') then
      ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
          ,a(irecind),b,flnm)

      irecind = irecind + irecsize
   end if

   ierr = RAMS_getvar('LEAF_CLASS',idim_type,ngrd   &
        ,a(irecind),b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_vegclass(irecsize,1,1,a(irecind))
   if (cvar(1:lv).eq.'vtype') then
      ivar_type = 7
   else
      ivar_type = 2
      call RAMS_comp_bigpatch(nnxp(ngrd),nnyp(ngrd),1,npatch  &
         ,a(irecind),a(1),b)
   endif
   cdname='vegetation class'
   cdunits='#'

elseif(cvar(1:lv).eq.'ndvi' .or. cvar(1:lv).eq.'ndvi_bp') then

   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch

   if(cvar(1:lv).eq.'ndvi_bp') then
      ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
          ,a(irecind),b,flnm)

      irecind = irecind + irecsize
   end if

   ierr = RAMS_getvar('VEG_NDVIC',idim_type,ngrd   &
        ,a(irecind),b,flnm)
   ierr_getvar = ierr_getvar + ierr

   if (cvar(1:lv).eq.'ndvi') then
      ivar_type = 7
   else
      ivar_type = 2
      call RAMS_comp_bigpatch(nnxp(ngrd),nnyp(ngrd),1,npatch  &
         ,a(irecind),a(1),b)
   endif

   cdname='ndvi'
   cdunits='#'


elseif(cvar(1:lv).eq.'qveg_class_p' .or. cvar(1:lv).eq.'qveg_class_bp') then

   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   if(cvar(1:lv).eq.'qveg_class_bp') then
      ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
           ,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr
   end if

   irecind = irecind + irecsize
   ierr = RAMS_getvar('DATQ_CLASS',idim_type,ngrd   &
        ,a(irecind),b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_vegclass(irecsize,1,1,a(irecind))

   if (cvar(1:lv).eq.'qveg_class_p') then
      ivar_type = 7
   else
      ivar_type = 2
      call RAMS_comp_bigpatch(nnxp(ngrd),nnyp(ngrd),1,npatch  &
         ,a(irecind),a(1),b)
   endif

   cdname='q vegetation class'
   cdunits='#'

elseif(cvar(1:lv).eq.'vegfrac' .or. cvar(1:lv).eq.'veg_fracarea_ps')then

   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch

   if (cvar(1:lv).eq.'veg_fracarea_ps')then
      ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
           ,a(irecind),b,flnm)

      irecind = irecind + irecsize
   end if

   ierr = RAMS_getvar('VEG_FRACAREA',idim_type,ngrd   &
        ,a(irecind),b,flnm)
   ierr_getvar = ierr_getvar + ierr

   if(cvar(1:lv).eq.'vegfrac') then
      ivar_type = 7
   else
      ivar_type = 2
      call RAMS_comp_patchsum_l(nnxp(ngrd),nnyp(ngrd),1,npatch  &
         ,a(irecind),a(1),b)
   endif

   cdname='vegetation frac area'
   cdunits=''

elseif(cvar(1:lv).eq.'land') then

   ivar_type = 2
   ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_1minus(nnxp(ngrd),nnyp(ngrd),1,a)
   cdname='land frac area'
   cdunits=''

elseif(cvar(1:lv).eq.'lai' .or. cvar(1:lv).eq.'veglai' .or. cvar(1:lv).eq.'lai_ps') then

   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   if (cvar(1:lv).eq.'lai_ps') then
      ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
           ,a(irecind),b,flnm)

      irecind = irecind + irecsize
   end if

   ierr = RAMS_getvar('VEG_LAI',idim_type,ngrd   &
        ,a(irecind),b,flnm)
   ierr_getvar = ierr_getvar + ierr

   if(cvar(1:lv).eq.'lai' .or. cvar(1:lv).eq.'veglai') then
      ivar_type = 7
   else
      ivar_type = 2
      call RAMS_comp_patchsum_l(nnxp(ngrd),nnyp(ngrd),1,npatch  &
         ,a(irecind),a(1),b)
   endif

   cdname='green leaf area index'
   cdunits=''


elseif(cvar(1:lv).eq.'tai' .or. cvar(1:lv).eq.'tai_ps') then

   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch

   if (cvar(1:lv).eq.'tai_ps') then
      ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
           ,a(irecind),b,flnm)

      irecind = irecind + irecsize
   end if

   ierr = RAMS_getvar('VEG_TAI',idim_type,ngrd   &
        ,a(irecind),b,flnm)
   ierr_getvar = ierr_getvar + ierr

   if(cvar(1:lv).eq.'tai') then
      ivar_type = 7
   else
      ivar_type = 2
      call RAMS_comp_patchsum_l(nnxp(ngrd),nnyp(ngrd),1,npatch  &
         ,a(irecind),a(1),b)
   endif

   cdname=' total leaf area index'
   cdunits=''

elseif(cvar(1:lv).eq.'z0' .or. cvar(1:lv).eq.'z0_ps') then

   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch

   if(cvar(1:lv).eq.'z0_ps') then
      ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
           ,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr
   end if

   irecind = irecind + irecsize
   ierr = RAMS_getvar('PATCH_ROUGH',idim_type,ngrd   &
        ,a(irecind),b,flnm)
   ierr_getvar = ierr_getvar + ierr

   if(cvar(1:lv).eq.'z0_p') then
      ivar_type = 7
   else
      ivar_type = 2
      call RAMS_comp_patchsum(nnxp(ngrd),nnyp(ngrd),1,npatch  &
         ,a(irecind),a(1),b)
   endif

   cdname='roughness'
   cdunits='m'


elseif(cvar(1:lv).eq.'net_z0_p' .or. cvar(1:lv).eq.'net_z0_ps') then

   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   if(cvar(1:lv).eq.'net_z0_ps') then
      ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
           ,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr
   end if

   irecind = irecind + irecsize
   ierr = RAMS_getvar('NET_Z0',idim_type,ngrd   &
        ,a(irecind),b,flnm)
   ierr_getvar = ierr_getvar + ierr

   if(cvar(1:lv).eq.'net_z0_p') then
      ivar_type = 7
   else
      ivar_type = 2
      call RAMS_comp_patchsum(nnxp(ngrd),nnyp(ngrd),1,npatch  &
         ,a(irecind),a(1),b)
   endif

   cdname='net roughness'
   cdunits='m'

elseif(cvar(1:lv).eq.'vegz0' .or. cvar(1:lv).eq.'vegz0_ps') then

   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch

   if(cvar(1:lv).eq.'vegz0_ps') then
      ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
           ,a(irecind),b,flnm)

      irecind = irecind + irecsize
   end if


   ierr = RAMS_getvar('VEG_ROUGH',idim_type,ngrd   &
        ,a(irecind),b,flnm)
   ierr_getvar = ierr_getvar + ierr

   if(cvar(1:lv).eq.'vegz0') then
      ivar_type = 7
   else
      ivar_type = 2
      call RAMS_comp_patchsum_l(nnxp(ngrd),nnyp(ngrd),1,npatch  &
         ,a(irecind),a(1),b)
   endif

   cdname='vegetation roughness'
   cdunits='m'

elseif(cvar(1:lv).eq.'vegdisp' .or. cvar(1:lv).eq.'veg_disp_ps') then

   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch

   if(cvar(1:lv).eq.'veg_disp_ps') then
      ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
           ,a(irecind),b,flnm)

      irecind = irecind + irecsize
   end if

   ierr = RAMS_getvar('VEG_DISP',idim_type,ngrd   &
        ,a(irecind),b,flnm)
   ierr_getvar = ierr_getvar + ierr

   if(cvar(1:lv).eq.'vegdisp') then
      ivar_type = 7
   else
      ivar_type = 2
      call RAMS_comp_patchsum_l(nnxp(ngrd),nnyp(ngrd),1,npatch  &
         ,a(irecind),a(1),b)
   endif

   cdname='vegetation displacement height'
   cdunits='m'

elseif(cvar(1:lv).eq.'patch_wetind') then

   ivar_type = 7
   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
        ,a(irecind),b,flnm)
   ierr_getvar = ierr_getvar + ierr

   irecind = irecind + irecsize
   ierr = RAMS_getvar('WET_INDEX',idim_type,ngrd   &
        ,a(irecind),b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname='patch wetness index'
   cdunits=''

elseif(cvar(1:lv).eq.'snowlevels') then

   ivar_type = 7
   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
        ,a(irecind),b,flnm)
   ierr_getvar = ierr_getvar + ierr

   irecind = irecind + irecsize
   ierr = RAMS_getvar('KSNOW',idim_type,ngrd   &
        ,a(irecind),b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname='number of snow levels'
   cdunits='#'

elseif(cvar(1:lv).eq.'grnd_mixrat_p' .or. cvar(1:lv).eq.'grnd_mixrat_ps') then

   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   if (cvar(1:lv).eq.'grnd_mixrat_ps') then
      ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
           ,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr
   end if

   irecind = irecind + irecsize
   ierr = RAMS_getvar('SFC_RS',idim_type,ngrd   &
        ,a(irecind),b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_mults(n1,n2,npatch,a(irecind),1.e3)

   if(cvar(1:lv).eq.'grnd_mixrat_p') then
      ivar_type = 7
   else
      ivar_type = 2
      call RAMS_comp_patchsum_l(nnxp(ngrd),nnyp(ngrd),1,npatch  &
         ,a(irecind),a(1),b)
   endif

   cdname='ground mixing ratio'
   cdunits='g/kg'

elseif(cvar(1:lv).eq.'soil_mixrat_p' .or. cvar(1:lv).eq.'soil_mixrat_ps') then

   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   if (cvar(1:lv).eq.'soil_mixrat_ps') then
      ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
           ,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr
   end if

   irecind = irecind + irecsize
   ierr = RAMS_getvar('SOIL_RS',idim_type,ngrd   &
        ,a(irecind),b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_mults(n1,n2,npatch,a(irecind),1.e3)

   if(cvar(1:lv).eq.'soil_mixrat_p') then
      ivar_type = 7
   else
      ivar_type = 2
      call RAMS_comp_patchsum_l(nnxp(ngrd),nnyp(ngrd),1,npatch  &
         ,a(irecind),a(1),b)
   endif

   cdname='soil mixing ratio'
   cdunits='g/kg'

elseif(cvar(1:lv).eq.'lwater_p' .or. cvar(1:lv).eq.'lwater_ps') then

   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   if (cvar(1:lv).eq.'lwater_ps') then
      ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
           ,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr
   irecind = irecind + irecsize
   end if

   ierr = RAMS_getvar('VEG_WATER',idim_type,ngrd   &
        ,a(irecind),b,flnm)
   ierr_getvar = ierr_getvar + ierr

   if(cvar(1:lv).eq.'lwater_p') then
      ivar_type = 7
   else
      ivar_type = 2
      call RAMS_comp_patchsum_l(nnxp(ngrd),nnyp(ngrd),1,npatch  &
         ,a(irecind),a(1),b)
   endif

    cdname='leaf water'
   cdunits='kg/m2'

elseif(cvar(1:lv).eq.'rib' .or. cvar(1:lv) .eq. 'rib_ps') then

   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   if (cvar(1:lv).eq.'rib_ps') then
      ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
           ,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr
      irecind = irecind + irecsize
   end if

   ierr = RAMS_getvar('VEG_ROUGH',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr

   ierr= RAMS_getvar('UP',idim_type,ngrd,d,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr= RAMS_getvar('VP',idim_type,ngrd,e,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_speed(n1,n2,n3,d,e)
   
   ierr = RAMS_getvar('THETA',idim_type,ngrd,e,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr = RAMS_getvar('RTP',idim_type,ngrd,f,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr = RAMS_getvar('RV',idim_type,ngrd,h,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_thetv(n1,n2,n3,e,f,h)
   
   ierr = RAMS_getvar('CAN_THETA',idim_type,ngrd,f,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr = RAMS_getvar('CAN_RVAP',idim_type,ngrd,h,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_thetv(n1,n2,npatch,f,h,h)

   ierr = RAMS_getvar('TOPT',idim_type,ngrd,h,b,flnm)
   ierr_getvar = ierr_getvar + ierr

   call RAMS_comp_richardson(n1,n2,n3,npatch,a(irecind),c,d,e,f,h,ngrd)

   if(cvar(1:lv).eq.'rib') then
      ivar_type = 7
   else
      ivar_type = 2
      call RAMS_comp_patchsum(nnxp(ngrd),nnyp(ngrd),1,npatch  &
         ,a(irecind),a(1),b)
   endif

   cdname='bulk richardson number'
   cdunits='n/d'



elseif(cvar(1:lv).eq.'rvcan' .or. cvar(1:lv).eq.'rvcan_ps') then

   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   if (cvar(1:lv).eq.'rvcan_ps') then
      ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
           ,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr
      irecind = irecind + irecsize
   end if

   ierr = RAMS_getvar('CAN_RVAP',idim_type,ngrd   &
        ,a(irecind),b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_mults(n1,n2,npatch,a(irecind),1.e3)

   if(cvar(1:lv).eq.'rvcan') then
      ivar_type = 7
   else
      ivar_type = 2
      call RAMS_comp_patchsum(nnxp(ngrd),nnyp(ngrd),1,npatch  &
         ,a(irecind),a(1),b)
   endif

   cdname='canopy mixing ratio'
   cdunits='g/kg'

elseif(cvar(1:lv).eq.'co2can' .or. cvar(1:lv).eq.'co2can_ps') then

   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   if (cvar(1:lv).eq.'co2can_ps') then
      ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
           ,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr
      irecind = irecind + irecsize
   end if

   ierr = RAMS_getvar('CAN_CO2',idim_type,ngrd   &
        ,a(irecind),b,flnm)
   ierr_getvar = ierr_getvar + ierr

   if(cvar(1:lv).eq.'co2can') then
      ivar_type = 7
   else
      ivar_type = 2
      call RAMS_comp_patchsum(nnxp(ngrd),nnyp(ngrd),1,npatch  &
         ,a(irecind),a(1),b)
   endif

   cdname='CO2 mixing ratio'
   cdunits='umol/mol'

elseif(cvar(1:lv).eq.'gpp_p' .or. cvar(1:lv).eq.'gpp') then

   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   if (cvar(1:lv).eq.'gpp') then
      ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
           ,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr
      irecind = irecind + irecsize
   end if

   ierr = RAMS_getvar('GPP',idim_type,ngrd,a(irecind),b,flnm)
   ierr_getvar = ierr_getvar + ierr

   if(cvar(1:lv).eq.'gpp_p') then
      ivar_type = 7
   else
      ivar_type = 2
      call RAMS_comp_patchsum(nnxp(ngrd),nnyp(ngrd),1,npatch  &
         ,a(irecind),a(1),b)
   endif

   cdname='Gross Primary Production'
   cdunits='umol/m2/s'

elseif(cvar(1:lv).eq.'plresp_p' .or. cvar(1:lv).eq.'plresp') then

   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   if (cvar(1:lv).eq.'plresp') then
      ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr
      irecind = irecind + irecsize
   end if

   ierr = RAMS_getvar('PLRESP',idim_type,ngrd,a(irecind),b,flnm)
   ierr_getvar = ierr_getvar + ierr

   if(cvar(1:lv).eq.'plresp_p') then
      ivar_type = 7
   else
      ivar_type = 2
      call RAMS_comp_patchsum(nnxp(ngrd),nnyp(ngrd),1,npatch  &
         ,a(irecind),a(1),b)
   endif

   cdname='Plant respiration'
   cdunits='umol/m2/s'

elseif(cvar(1:lv).eq.'resphet_p' .or. cvar(1:lv).eq.'resphet') then

   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   if (cvar(1:lv).eq.'resphet') then
      ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr
      irecind = irecind + irecsize
   end if

   ierr = RAMS_getvar('RESPHET',idim_type,ngrd,a(irecind),b,flnm)
   ierr_getvar = ierr_getvar + ierr

   if(cvar(1:lv).eq.'resphet_p') then
      ivar_type = 7
   else
      ivar_type = 2
      call RAMS_comp_patchsum(nnxp(ngrd),nnyp(ngrd),1,npatch  &
         ,a(irecind),a(1),b)
   endif

   cdname='Heterotrophic respiration'
   cdunits='umol/m2/s'

elseif(cvar(1:lv).eq.'h_p' .or. cvar(1:lv).eq.'h') then

   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   if (cvar(1:lv).eq.'h') then
      ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr
      irecind = irecind + irecsize
   end if

   ierr = RAMS_getvar('SENSIBLE',idim_type,ngrd,a(irecind),b,flnm)
   ierr_getvar = ierr_getvar + ierr

   if(cvar(1:lv).eq.'h_p') then
      ivar_type = 7
   else
      ivar_type = 2
      call RAMS_comp_patchsum(nnxp(ngrd),nnyp(ngrd),1,npatch  &
         ,a(irecind),a(1),b)
   endif

   cdname='Sensible heat flux'
   cdunits='W/m2'

elseif(cvar(1:lv).eq.'evap_p' .or. cvar(1:lv).eq.'evap') then

   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   if (cvar(1:lv).eq.'evap') then
      ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr
      irecind = irecind + irecsize
   end if

   ierr = RAMS_getvar('EVAP',idim_type,ngrd,a(irecind),b,flnm)
   ierr_getvar = ierr_getvar + ierr

   if(cvar(1:lv).eq.'evap_p') then
      ivar_type = 7
   else
      ivar_type = 2
      call RAMS_comp_patchsum(nnxp(ngrd),nnyp(ngrd),1,npatch  &
         ,a(irecind),a(1),b)
   endif

   cdname='Latent heat flux due to evaporation'
   cdunits='W/m2'

elseif(cvar(1:lv).eq.'transp_p' .or. cvar(1:lv).eq.'transp') then

   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   if (cvar(1:lv).eq.'transp') then
      ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr
      irecind = irecind + irecsize
   end if

   ierr = RAMS_getvar('TRANSP',idim_type,ngrd,a(irecind),b,flnm)
   ierr_getvar = ierr_getvar + ierr

   if(cvar(1:lv).eq.'transp_p') then
      ivar_type = 7
   else
      ivar_type = 2
      call RAMS_comp_patchsum(nnxp(ngrd),nnyp(ngrd),1,npatch  &
         ,a(irecind),a(1),b)
   endif

   cdname='Latent heat flux due to transpiration'
   cdunits='W/m2'

elseif(cvar(1:lv).eq.'le_p' .or. cvar(1:lv).eq.'le') then

   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   if (cvar(1:lv).eq.'le') then
      ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr
      irecind = irecind + irecsize
   end if

   ierr = RAMS_getvar('TRANSP',idim_type,ngrd,a(irecind),b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr = RAMS_getvar('EVAP',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_accum(n1,n2,npatch,a(irecind),c)

   if(cvar(1:lv).eq.'le_p') then
      ivar_type = 7
   else
      ivar_type = 2
      call RAMS_comp_patchsum(nnxp(ngrd),nnyp(ngrd),1,npatch  &
         ,a(irecind),a(1),b)
   endif

   cdname='Latent heat flux'
   cdunits='W/m2'

elseif(cvar(1:lv).eq.'tveg' .or. cvar(1:lv).eq.'tveg_ps') then

   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   if(cvar(1:lv).eq.'tveg_ps') then
      ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
           ,a(irecind),b,flnm)

      irecind = irecind + irecsize
   end if
   
   ierr = RAMS_getvar('CAN_THETA',idim_type,ngrd,e,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr = RAMS_getvar('CAN_PRSS',idim_type,ngrd,d,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_theta2temp(n1,n2,npatch,e,d)
   call RAMS_comp_tempC(n1,n2,npatch,e)
   
   ierr = RAMS_getvar('VEG_ENERGY',idim_type,ngrd,a(irecind),b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr = RAMS_getvar('VEG_WATER',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr = RAMS_getvar('VEG_HCAP',idim_type,ngrd,d,b,flnm)
   ierr_getvar = ierr_getvar + ierr

   call RAMS_comp_tvegc(n1,n2,npatch,a(irecind),c,d,e)

   !----- Filling first patch with SST. -------------------------------!
   ierr = RAMS_getvar('SOIL_ENERGY',idim_type,ngrd,e,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   kp = nzg
   call rams_fill_sst(n1,n2,nzg*npatch,kp,a(irecind),e)

   if(cvar(1:lv).eq.'tveg') then
      ivar_type = 7
   else
      ivar_type = 2
      call RAMS_comp_patchsum(nnxp(ngrd),nnyp(ngrd),1,npatch  &
         ,a(irecind),a(1),b)
   endif

   cdname='vegetation temperature'
   cdunits='C'

elseif(cvar(1:lv).eq.'tcan' .or. cvar(1:lv).eq.'tcan_ps') then

   irecind = 1
   if(cvar(1:lv).eq.'tcan_ps') then
      irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
      ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
           ,a(irecind),b,flnm)
   end if
   irecind = irecind + irecsize
   ierr = RAMS_getvar('CAN_THETA',idim_type,ngrd,a(irecind),b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr = RAMS_getvar('CAN_PRSS',idim_type,ngrd,c,b,flnm)
   call RAMS_comp_theta2temp(n1,n2,npatch,a(irecind),c)
   call RAMS_comp_tempC(n1,n2,npatch,a(irecind))

   if(cvar(1:lv).eq.'tcan') then
      ivar_type = 7
   else
      ivar_type = 2
      call RAMS_comp_patchsum(nnxp(ngrd),nnyp(ngrd),1,npatch  &
         ,a(irecind),a(1),b)
   endif

   cdname='canopy temperature'
   cdunits='C'


elseif(cvar(1:lv).eq.'pcan' .or. cvar(1:lv).eq.'pcan_ps') then

   irecind = 1
   if(cvar(1:lv).eq.'pcan_ps') then
      irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
      ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
           ,a(irecind),b,flnm)
   end if
   irecind = irecind + irecsize
   ierr = RAMS_getvar('CAN_PRSS',idim_type,ngrd   &
        ,a(irecind),b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_mults(n1,n2,npatch,a(irecind),1.e-2)

   if(cvar(1:lv).eq.'pcan') then
      ivar_type = 7
   else
      ivar_type = 2
      call RAMS_comp_patchsum(nnxp(ngrd),nnyp(ngrd),1,npatch  &
         ,a(irecind),a(1),b)
   endif

   cdname='canopy press'
   cdunits='hPa'

elseif(cvar(1:lv).eq.'thcan' .or. cvar(1:lv).eq.'thcan_ps') then

   irecind = 1
   if(cvar(1:lv).eq.'thcan_ps') then
      irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
      ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
           ,a(irecind),b,flnm)
   end if
   irecind = irecind + irecsize
   ierr = RAMS_getvar('CAN_THETA',idim_type,ngrd   &
        ,a(irecind),b,flnm)
   ierr_getvar = ierr_getvar + ierr

   if(cvar(1:lv).eq.'thcan') then
      ivar_type = 7
   else
      ivar_type = 2
      call RAMS_comp_patchsum(nnxp(ngrd),nnyp(ngrd),1,npatch  &
         ,a(irecind),a(1),b)
   endif

   cdname='canopy potential temperature'
   cdunits='K'

! sib - stuffs
!itb...src_co2
elseif(cvar(1:lv).eq.'src_co2') then
   ivar_type=2
   ierr= RAMS_getvar('SRC_CO2',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   if(ierr.eq.0) then
      call RAMS_comp_mults(n1,n2,n3,a,1.e6)
      call RAMS_comp_noneg(n1,n2,n3,a)
   endif
   cdname='CO2 flux'
   cdunits='umol/m**2/sec'

  elseif(cvar(1:lv).eq.'CO2_SIB') then
   ivar_type=3
   ierr= RAMS_getvar('SCLP001',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname='CO2 Concentration'
   cdunits='ppm'



elseif(cvar(1:lv).eq.'pco2ap' ) then
   ivar_type=2
   ierr= RAMS_getvar('pco2ap',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
    cdname='CAS CO2'
   cdunits='Pa'

elseif(cvar(1:lv).eq.'pco2m' ) then
   ivar_type=2
   ierr= RAMS_getvar('pco2m',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
    cdname='REF LEVEL CO2'
   cdunits='Pa'


!itb...
elseif(cvar(1:lv).eq.'rst') then
   ivar_type=2
   ierr= RAMS_getvar('rst',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname='stomatal resistance'
   cdunits='sec/meter'

!...itb - NEW DIAGNOSTICS...

elseif(cvar(1:lv).eq.'fss') then
   ivar_type=2
   ierr= RAMS_getvar('fss',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname='sensible heat flux'
   cdunits='W/m^2'

elseif(cvar(1:lv).eq.'fws') then
   ivar_type=2
   ierr= RAMS_getvar('fws',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname='latent heat flux'
   cdunits='kg H2O/m^2/sec'

elseif(cvar(1:lv).eq.'assimn') then
   ivar_type=2
   ierr= RAMS_getvar('assimn',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname='canopy net assimilation'
   cdunits='mol/m^2/sec'

elseif(cvar(1:lv).eq.'respg') then
   ivar_type=2
   ierr= RAMS_getvar('respg',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname='ground respiration'
   cdunits='mol/m^2/sec'

elseif(cvar(1:lv).eq.'rstfac1') then
   ivar_type=2
   ierr= RAMS_getvar('rstfac1',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname='stress factor 1-leaf to CAS humidity'
   cdunits='(-)'

elseif(cvar(1:lv).eq.'rstfac2') then
   ivar_type=2
   ierr= RAMS_getvar('rstfac2',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname='stress factor 2-soil moisture'
   cdunits='(-)'

elseif(cvar(1:lv).eq.'rstfac3') then
   ivar_type=2
   ierr= RAMS_getvar('rstfac3',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname='stress factor 3-temperature'
   cdunits='(-)'

elseif(cvar(1:lv).eq.'rstfac4') then
   ivar_type=2
   ierr= RAMS_getvar('rstfac4',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname='stress factor 4-combination of factors 1-3'
   cdunits='(-)'

elseif(cvar(1:lv).eq.'ect') then
   ivar_type=2
   ierr= RAMS_getvar('ect',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname='canopy transpiration'
   cdunits='W/m^2'

elseif(cvar(1:lv).eq.'eci') then
   ivar_type=2
   ierr= RAMS_getvar('eci',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname='canopy interception evaporation'
   cdunits='W/m^2'

elseif(cvar(1:lv).eq.'egi') then
   ivar_type=2
   ierr= RAMS_getvar('egi',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname='ground interception evaporation'
   cdunits='W/m^2'

elseif(cvar(1:lv).eq.'egs') then
   ivar_type=2
   ierr= RAMS_getvar('egs',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname='top soil layer evaporation'
   cdunits='W/m^2'

elseif(cvar(1:lv).eq.'hc') then
   ivar_type=2
   ierr= RAMS_getvar('hc',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname='canopy sensible heat flux'
   cdunits='W/m^2'

elseif(cvar(1:lv).eq.'hg') then
   ivar_type=2
   ierr= RAMS_getvar('hg',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname='ground sensible heat flux'
   cdunits='W/m^2'


elseif(cvar(1:lv).eq.'capac1' ) then
   ivar_type=2
   ierr= RAMS_getvar('capac1',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname='VEGETATION INTERCEPTION STORE'
   cdunits='kg/m^2'


elseif(cvar(1:lv).eq.'capac2'.or. cvar(1:lv).eq.'capac2_ps' ) then

   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   if (cvar(1:lv).eq.'capac2_ps' ) then
      ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
           ,a(irecind),b,flnm)

      irecind = irecind + irecsize
   end if

   ierr = RAMS_getvar('capac2',idim_type,ngrd   &
        ,a(irecind),b,flnm)
   ierr_getvar = ierr_getvar + ierr

   if(cvar(1:lv).eq.'capac2') then
      ivar_type = 7
   else
      ivar_type = 2
      call RAMS_comp_patchsum_l(nnxp(ngrd),nnyp(ngrd),1,npatch  &
         ,a(irecind),a(1),b)
   endif

   call RAMS_comp_noneg(n1,n2,n3,a)

   cdname='GROUND INTERCEPTION STORE'
   cdunits='kg/m^2'



elseif(cvar(1:lv).eq.'ustar' .or. cvar(1:lv).eq.'ustar_ps') then

   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch

   if (cvar(1:lv).eq.'ustar_ps') then
      ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
           ,a(irecind),b,flnm)
      irecind = irecind + irecsize
   end if

   ierr = RAMS_getvar('USTAR',idim_type,ngrd,a(irecind),b,flnm)
   ierr_getvar = ierr_getvar + ierr

   if(cvar(1:lv).eq.'ustar') then
      ivar_type = 7
   else
      ivar_type = 2
      call RAMS_comp_patchsum(nnxp(ngrd),nnyp(ngrd),1,npatch  &
         ,a(irecind),a(1),b)
   endif

   cdname='ustar'
   cdunits='m/s'

elseif(cvar(1:lv).eq.'tstar' .or. cvar(1:lv).eq.'tstar_ps') then

   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch

   if (cvar(1:lv).eq.'tstar_ps') then
      ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
           ,a(irecind),b,flnm)
      irecind = irecind + irecsize
   end if

   ierr = RAMS_getvar('TSTAR',idim_type,ngrd,a(irecind),b,flnm)
   ierr_getvar = ierr_getvar + ierr

   if(cvar(1:lv).eq.'tstar') then
      ivar_type = 7
   else
      ivar_type = 2
      call RAMS_comp_patchsum(nnxp(ngrd),nnyp(ngrd),1,npatch  &
         ,a(irecind),a(1),b)
   endif

   cdname='tstar'
   cdunits='K'

elseif(cvar(1:lv).eq.'rstar' .or. cvar(1:lv).eq.'rstar_ps') then

   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch

   if(cvar(1:lv).eq.'rstar_ps') then
      ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   & 
           ,a(irecind),b,flnm)                           
      irecind = irecind + irecsize
   end if

   ierr = RAMS_getvar('RSTAR',idim_type,ngrd   &
        ,a(irecind),b,flnm)
   ierr_getvar = ierr_getvar + ierr

   if(cvar(1:lv).eq.'rstar') then
      ivar_type = 7
   else
      ivar_type = 2
      call RAMS_comp_patchsum(nnxp(ngrd),nnyp(ngrd),1,npatch  &
         ,a(irecind),a(1),b)
   endif

   cdname='rstar'
   cdunits='kg/kg'

elseif(cvar(1:lv).eq.'cstar' .or. cvar(1:lv).eq.'cstar_ps') then

   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch

   if(cvar(1:lv).eq.'cstar_ps') then
      ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   & 
           ,a(irecind),b,flnm)                           
      irecind = irecind + irecsize
   end if

   ierr = RAMS_getvar('CSTAR',idim_type,ngrd   &
        ,a(irecind),b,flnm)
   ierr_getvar = ierr_getvar + ierr

   if(cvar(1:lv).eq.'cstar') then
      ivar_type = 7
   else
      ivar_type = 2
      call RAMS_comp_patchsum(nnxp(ngrd),nnyp(ngrd),1,npatch  &
         ,a(irecind),a(1),b)
   endif

   cdname='cstar'
   cdunits='umol/mol'

elseif(cvar(1:lv).eq.'snow_depth_p' .or. cvar(1:lv).eq.'snow_depth_ps') then

   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch

   if (cvar(1:lv).eq.'snow_depth_ps') then
      ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
           ,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr
   end if

   irecind = irecind + irecsize
   ierr = RAMS_getvar('SNOW_DEPTH',idim_type,ngrd   &
        ,a(irecind),b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_sum_snowlayers(nnxp(ngrd)*nnyp(ngrd),nzs,npatch,a(irecind))

   if(cvar(1:lv).eq.'snow_depth_p') then
      ivar_type = 7
   else
      ivar_type = 2
      call RAMS_comp_patchsum_l(nnxp(ngrd),nnyp(ngrd),1,npatch  &
         ,a(irecind),a(1),b)
   endif

   cdname='snow depth'
   cdunits='m'

elseif(cvar(1:lv).eq.'snowcover_p' .or. cvar(1:lv).eq.'snowcover_ps') then

   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   if (cvar(1:lv).eq.'snowcover_ps') then
      ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
           ,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr
   end if

   irecind = irecind + irecsize
   ierr = RAMS_getvar('SNOW_MOIST',idim_type,ngrd   &
        ,a(irecind),b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_sum_snowlayers(nnxp(ngrd)*nnyp(ngrd),nzs,npatch,a(irecind))

   if(cvar(1:lv).eq.'snow_depth_p') then
      ivar_type = 7
   else
      ivar_type = 2
      call RAMS_comp_patchsum_l(nnxp(ngrd),nnyp(ngrd),1,npatch  &
         ,a(irecind),a(1),b)
   endif

   cdname='snowcover'
   cdunits='kg/m2'

elseif(cvar(1:lv).eq.'sltex' .or. cvar(1:lv).eq.'sltex_bp') then

   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   if (cvar(1:lv).eq.'sltex_bp') then
      ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
           ,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr
   end if
   irecind = irecind + irecsize
   ierr = RAMS_getvar('SOIL_TEXT',idim_type,ngrd   &
        ,a(irecind),b,flnm)
   ierr_getvar = ierr_getvar + ierr

   if (cvar(1:lv).eq.'sltex') then
      ivar_type = 8
   else
      ivar_type = 10
      call RAMS_comp_bigpatch(nnxp(ngrd),nnyp(ngrd),nzg,npatch  &
         ,a(irecind),a(1),b)
   endif

   cdname='soil textural class'
   cdunits='#'

elseif(cvar(1:lv).eq.'soilq' .or. cvar(1:lv).eq.'soilq_ps') then

   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   
   if(cvar(1:lv).eq.'soilq_ps') then
      ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
           ,a(irecind),b,flnm)
      irecind = irecind + irecsize
   end if
   ierr = RAMS_getvar('SOIL_ENERGY',idim_type,ngrd   &
        ,a(irecind),b,flnm)
   ierr_getvar = ierr_getvar + ierr

   call get_leaf_soil(n1,n2,n3,n4,n5,a(irecind),a2)

   if (cvar(1:lv).eq.'soilq') then
      ivar_type = 8
   else
      ivar_type = 10
      call RAMS_comp_patchsum_l(nnxp(ngrd),nnyp(ngrd),nzg,npatch  &
         ,a(irecind),a(1),b)
   endif

   cdname='soil q'
   cdunits='J/m3'

elseif(cvar(1:lv).eq.'smoist' .or. cvar(1:lv).eq.'smoist_ps') then

   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   if(cvar(1:lv).eq.'smoist_ps') then
      ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
           ,a(irecind),b,flnm)
      irecind = irecind + irecsize
   end if
   ierr = RAMS_getvar('SOIL_WATER',idim_type,ngrd   &
        ,a(irecind),b,flnm)
   ierr_getvar = ierr_getvar + ierr

   call get_leaf_soil(n1,n2,n3,n4,n5,a(irecind),a2)

   if (cvar(1:lv).eq.'smoist') then
      ivar_type = 8
   else
      ivar_type = 10
      call RAMS_comp_patchsum_l(nnxp(ngrd),nnyp(ngrd),nzg,npatch  &
         ,a(irecind),a(1),b)
   endif

   cdname='soil moisture'
   cdunits='m3/m3'


elseif(cvar(1:lv).eq.'tsoil' .or. cvar(1:lv).eq.'tsoil_ps') then

   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch

   if(cvar(1:lv).eq.'tsoil' .or. cvar(1:lv).eq.'tsoil_ps') then
      ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
           ,a(irecind),b,flnm)
      irecind = irecind + irecsize
   end if

   ierr = RAMS_getvar('SOIL_ENERGY',idim_type,ngrd   &
        ,a(irecind),b,flnm)
   ierr_getvar = ierr_getvar + ierr

   ierr = RAMS_getvar('SOIL_WATER',idim_type,ngrd   &
        ,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr = RAMS_getvar('SOIL_TEXT',idim_type,ngrd   &
        ,d,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_copysst(n1,n2,nzg,a(irecind))

   irecsizep = nnxp(ngrd) * nnyp(ngrd) * nzg
   call RAMS_comp_qwtc(n1,n2,nzg*(npatch-1),a(irecind+irecsizep)  &
      ,c(1+irecsizep),d(1+irecsizep))

   call get_leaf_soil(n1,n2,n3,n4,n5,a(irecind),a2)

   if (cvar(1:lv).eq.'tsoil') then
      ivar_type = 8
   else
      ivar_type = 10
      call RAMS_comp_patchsum(nnxp(ngrd),nnyp(ngrd),nzg,npatch  &
         ,a(irecind),a(1),b)
   endif

   cdname='soil/sea temp'
   cdunits='C'

elseif(cvar(1:lv).eq.'smfrac' .or. cvar(1:lv).eq.'smfrac_ps') then

   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   if (cvar(1:lv).eq.'smfrac_ps') then
      ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
           ,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr
   end if
   irecind = irecind + irecsize
   ierr = RAMS_getvar('SOIL_WATER',idim_type,ngrd   &
        ,a(irecind),b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr = RAMS_getvar('SOIL_TEXT',idim_type,ngrd   &
        ,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   
   
   call rams_comp_slmstf(n1,n2,nzg*npatch,a(irecind),c)

   call get_leaf_soil(n1,n2,n3,n4,n5,a(irecind),a2)

   if (cvar(1:lv).eq.'smfrac') then
      ivar_type = 8
   else
      ivar_type = 10
      call RAMS_comp_patchsum_l(nnxp(ngrd),nnyp(ngrd),nzg,npatch  &
         ,a(irecind),a(1),b)
   endif

   cdname='soil moisture frac'
   cdunits='m3/m3'

elseif(cvar(1:lv).eq.'sfcw_temp') then
   ivar_type=2
   ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr = RAMS_getvar('SFCWATER_ENERGY',idim_type,ngrd,d,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr = RAMS_getvar('SFCWATER_MASS',idim_type,ngrd,e,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr = RAMS_getvar('SFCWATER_NLEV',idim_type,ngrd,f,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_sfcwmeantemp(n1,n2,nzs,npatch,c,d,e,f,a)
   ierr_getvar = ierr_getvar + ierr
   cdname='Pond/snow mean temperature'
   cdunits='C'

elseif(cvar(1:lv).eq.'sfcw_mass') then
   ivar_type=2
   ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr = RAMS_getvar('SFCWATER_MASS',idim_type,ngrd,e,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr = RAMS_getvar('SFCWATER_NLEV',idim_type,ngrd,f,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_sfcwinteg(n1,n2,nzs,npatch,c,e,f,a)
   call RAMS_comp_noneg(n1,n2,1,a)
   ierr_getvar = ierr_getvar + ierr
   cdname='Pond/snow mass'
   cdunits='kg/m2'

elseif(cvar(1:lv).eq.'sfcw_depth') then
   ivar_type=2
   ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr = RAMS_getvar('SFCWATER_DEPTH',idim_type,ngrd,d,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr = RAMS_getvar('SFCWATER_NLEV',idim_type,ngrd,f,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_sfcwinteg(n1,n2,nzs,npatch,c,d,f,a)
   call RAMS_comp_noneg(n1,n2,1,a)
   ierr_getvar = ierr_getvar + ierr
   cdname='Pond/snow depth'
   cdunits='m'

elseif(cvar(1:lv).eq.'leaf2_moisture') then

! These values should somehow be scaled across soil, snow, vegetation, and canopy air
! using calls to rams_comp_snownorm, rams_comp_vegnorm, and rams_comp_cannorm,
! which are not yet completed.

   ivar_type = 10

   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
        ,a(irecind),b,flnm)
   ierr_getvar = ierr_getvar + ierr

   irecind = irecind + irecsize
   irecsize = nnxp(ngrd) * nnyp(ngrd) * nzg * npatch
   ierr = RAMS_getvar('SOIL_WATER',idim_type,ngrd   &
        ,a(irecind),b,flnm)
   ierr_getvar = ierr_getvar + ierr
   ierr = RAMS_getvar('SOIL_TEXT',idim_type,ngrd   &
        ,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call rams_comp_slmstf(n1,n2,nzg*npatch,a(irecind),c)

   irecind = irecind + irecsize
   irecsize = nnxp(ngrd) * nnyp(ngrd) * nzs * npatch
   ierr = RAMS_getvar('SNOW_MOIST',idim_type,ngrd   &
        ,a(irecind),b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call rams_comp_snownorm(irecsize,1,1,a(irecind),c)

   irecind = irecind + irecsize
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   ierr = RAMS_getvar('VEG_WATER',idim_type,ngrd   &
        ,a(irecind),b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call rams_comp_vegnorm(irecsize,1,1,a(irecind),c)

   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   irecind = irecind + irecsize
   ierr = RAMS_getvar('CAN_RV',idim_type,ngrd   &
        ,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call rams_comp_cannorm(irecsize,1,1,a(irecind),c)

   ! also get pcpg, and vapor for k=2

   cdname='leaf2 moisture frac'
   cdunits='m3/m3'


! CATT

elseif(cvar(1:lv).eq.'CO') then
   ivar_type=3
   ierr= RAMS_getvar('SCLP001',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_noneg(n1,n2,n3,a)
   call RAMS_comp_mults(n1,n2,n3,a,1.e-6)  ! converte de mg/kg para kg/kg
  call RAMS_transf_ppb(n1,n2,n3,a)
   cdname='CO Concentration'
   cdunits='ppb'

   elseif(cvar(1:lv).eq.'src1') then
   ivar_type=3
   ierr= RAMS_getvar('scrsc001',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_mults(n1,n2,n3,a,1.e-6)  ! converte de mg/day para kg/day
   cdname='emission 1'
   cdunits='kg/m2/day'

   elseif(cvar(1:lv).eq.'src2') then
   ivar_type=3
   ierr= RAMS_getvar('scrsc002',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_mults(n1,n2,n3,a,1.e-6)  ! converte de mg/day para kg/day
   cdname='emission 2'
   cdunits='kg/m2/day'

   elseif(cvar(1:lv).eq.'src3') then
   ivar_type=3
   ierr= RAMS_getvar('scrsc003',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_mults(n1,n2,n3,a,1.e-6)  ! converte de mg/day para kg/day
   cdname='emission 3'
   cdunits='kg/m2/day'

   elseif(cvar(1:lv).eq.'src4') then
   ivar_type=3
   ierr= RAMS_getvar('scrsc004',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_mults(n1,n2,n3,a,1.e-6)  ! converte de mg/day para kg/day
   cdname='emission 4'
   cdunits='kg/m2/day'

   elseif(cvar(1:lv).eq.'src5') then
   ivar_type=3
   ierr= RAMS_getvar('scrsc005',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_mults(n1,n2,n3,a,1.e-6)  ! converte de mg/day para kg/day
   cdname='emission 5'
   cdunits='kg/m2/day'

   elseif(cvar(1:lv).eq.'src6') then
   ivar_type=3
   ierr= RAMS_getvar('scrsc006',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_mults(n1,n2,n3,a,1.e-6)  ! converte de mg/day para kg/day
   cdname='emission 6'
   cdunits='kg/m2/day'

    elseif(cvar(1:lv).eq.'src7') then
   ivar_type=3
   ierr= RAMS_getvar('scrsc007',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_mults(n1,n2,n3,a,1.e-6)  ! converte de mg/day para kg/day
   cdname='emission 7'
   cdunits='kg/m2/day'

   elseif(cvar(1:lv).eq.'src8') then
   ivar_type=3
   ierr= RAMS_getvar('scrsc008',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_mults(n1,n2,n3,a,1.e-6)  ! converte de mg/day para kg/day
   cdname='emission 5'
   cdunits='kg/m2/day'
  
elseif(cvar(1:lv).eq.'COstc') then
   ivar_type=3
   ierr= RAMS_getvar('SCLP002',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_noneg(n1,n2,n3,a)
   call RAMS_comp_mults(n1,n2,n3,a,1.e-6)  ! converte de mg/kg para kg/kg
  call RAMS_transf_ppb(n1,n2,n3,a)
   cdname='CO Conc. without conv. transp'
   cdunits='ppb'

elseif(cvar(1:lv).eq.'COANT') then
   ivar_type=3
   ierr= RAMS_getvar('SCLP004',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_noneg(n1,n2,n3,a)
   call RAMS_comp_mults(n1,n2,n3,a,1.e-6)  ! converte de mg/kg para kg/kg
   call RAMS_transf_ppb(n1,n2,n3,a)
   cdname='CO Concentration ANTRO'
   cdunits='ppb'

elseif(cvar(1:lv).eq.'COTOT') then
   ivar_type=3
   ierr= RAMS_getvar('SCLP005',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_noneg(n1,n2,n3,a)
   call RAMS_comp_mults(n1,n2,n3,a,1.e-6)  ! converte de mg/kg para kg/kg
  call RAMS_transf_ppb(n1,n2,n3,a)
   cdname='CO Conc ANTRO+BB'
   cdunits='ppb'

elseif(cvar(1:lv).eq.'PM25') then
   ivar_type=3
   ierr= RAMS_getvar('SCLP003',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_noneg(n1,n2,n3,a)
   call RAMS_comp_mults(n1,n2,n3,a,1.e-6)  ! converte de mg/kg para kg/kg
!air density
   ierr= RAMS_getvar('TOPT',idim_type,ngrd,e,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_dn0(n1,n2,n3,b,c,d,e,ngrd)

   call RAMS_transf_ugm3(n1,n2,n3,a,d)
   call RAMS_comp_noneg(n1,n2,n3,a)
   cdname='PM25 Concentration'
   cdunits='ug/m3'


elseif(cvar(1:lv).eq.'PMINT') then
   ivar_type=2
   ierr= RAMS_getvar('SCLP003',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_noneg(n1,n2,n3,a)
   call RAMS_comp_mults(n1,n2,n3,a,1.e-6)  ! converte de mg/kg para kg/kg
!air density
   ierr= RAMS_getvar('TOPT',idim_type,ngrd,e,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_dn0(n1,n2,n3,b,c,d,e,ngrd)
   call RAMS_comp_mult(n1,n2,n3,a,d) !Unit: kg[pm25]/m3
   call RAMS_comp_vertint(n1,n2,n3,a,e,ngrd) ! Unit: kg[pm25]/m2
   call RAMS_comp_mults(n1,n2,n3,a,1.e+9)  ! converte de kg/m2 para ug/m2

   cdname='PM25 vert int'
   cdunits='ug/m2'

! ------------------ AOT ------------------ 
! WAVE / 0.256, 0.280, 0.296, 0.319, 0.335, 0.365, 0.420, 0.482,
!	 0.598, 0.690, 0.762, 0.719, 0.813, 0.862, 0.926, 1.005,
!	 1.111, 1.333, 1.562, 1.770, 2.051, 2.210, 2.584, 3.284,
!	 3.809, 4.292,
!	 4.546, 4.878, 5.128, 5.405, 5.714, 6.061, 6.452, 6.897,
!	 7.407, 8.333, 9.009, 10.309,12.500,13.889,16.667,
!	 20.000, 26.316, 35.714, 62.50  		       
elseif(cvar(1:lv).eq.'aot256') then
   ivar_type=2
   ierr= RAMS_getvar('AOT',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call D3toD2(n1,n2,nwave,1,a,c)
   cdname='AOT 256nm'
   cdunits=' '

elseif(cvar(1:lv).eq.'aot296') then
   ivar_type=2
   ierr= RAMS_getvar('AOT',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call D3toD2(n1,n2,nwave,3,a,c)
   cdname='AOT 296nm'
   cdunits=' '

elseif(cvar(1:lv).eq.'aot335') then
   ivar_type=2
   ierr= RAMS_getvar('AOT',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call D3toD2(n1,n2,nwave,5,a,c)
   cdname='AOT 335nm'
   cdunits=' '

elseif(cvar(1:lv).eq.'aot420') then
   ivar_type=2
   ierr= RAMS_getvar('AOT',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call D3toD2(n1,n2,nwave,7,a,c)
   cdname='AOT 420nm'
   cdunits=' '

elseif(cvar(1:lv).eq.'aot482') then
   ivar_type=2
   ierr= RAMS_getvar('AOT',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call D3toD2(n1,n2,nwave,8,a,c)
   cdname='AOT 482nm'
   cdunits=' '


elseif(cvar(1:lv).eq.'aot598') then
   ivar_type=2
   ierr= RAMS_getvar('AOT',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call D3toD2(n1,n2,nwave,9,a,c)
   cdname='AOT 598nm'
   cdunits=' '

elseif(cvar(1:lv).eq.'aot690') then
   ivar_type=2
   ierr= RAMS_getvar('AOT',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call D3toD2(n1,n2,nwave,10,a,c)
   cdname='AOT 690nm'
   cdunits=' '

elseif(cvar(1:lv).eq.'aot500') then
   ivar_type=2
   ierr= RAMS_getvar('AOT',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call D3toD2(n1,n2,nwave,11,a,c)
   cdname='AOT 500nm'
   cdunits=' '

elseif(cvar(1:lv).eq.'aot550') then
   ivar_type=2
   ierr= RAMS_getvar('AOT',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call D3toD2(n1,n2,nwave,12,a,c)
   cdname='AOT 550nm'
   cdunits=' '


elseif(cvar(1:lv).eq.'secog') then
   ivar_type=2
   ierr= RAMS_getvar('DUM1',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call D3toD2(n1,n2,n3,2,a,c)
   call RAMS_comp_mults(n1,n2,n3,a,1.e-6)  ! converte de mg/day para kg/day
   cdname='GOES-8 ABBA CO emission'
   cdunits='kg/m2/day'


elseif(cvar(1:lv).eq.'secod') then
   ivar_type=2
   ierr= RAMS_getvar('DUM1',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call D3toD2(n1,n2,n3,11,a,c)
   call RAMS_comp_mults(n1,n2,n3,a,1.e-6)  ! converte de mg/day para kg/day
   cdname='Duncan CO emission'
   cdunits='kg/m2/day'

elseif(cvar(1:lv).eq.'secoant') then
   ivar_type=2
   ierr= RAMS_getvar('DUM1',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call D3toD2(n1,n2,n3,11,a,c)
   call RAMS_comp_mults(n1,n2,n3,a,1.e-6)  ! converte de mg/day para kg/day
   cdname='Antropogenic CO emission'
   cdunits='kg/m2/day'

elseif(cvar(1:lv).eq.'secoe') then
   ivar_type=2
   ierr= RAMS_getvar('DUM1',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call D3toD2(n1,n2,n3,14,a,c)
   call RAMS_comp_mults(n1,n2,n3,a,1.e-6)  ! converte de mg/day para kg/day
   cdname='EDGAR CO emission'
   cdunits='kg/m2/day'


elseif(cvar(1:lv).eq.'scco') then
   ivar_type=2
   ierr= RAMS_getvar('QSC1',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname='Massa de CO emitida'
   cdunits='kg/(m2 day)'

elseif(cvar(1:lv).eq.'scpm25') then
   ivar_type=2
   ierr= RAMS_getvar('QSC2',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname='Massa de PM25 emitida'
   cdunits='kg/(m2 day)'

elseif(cvar(1:lv).eq.'sccofe') then
   ivar_type=2
   ierr= RAMS_getvar('QSC3',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname='Massa de CO FWB - EDGAR emitida'
   cdunits='kg/(m2 day)'

elseif(cvar(1:lv).eq.'sccoae') then
   ivar_type=2
   ierr= RAMS_getvar('QSC4',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname='Massa de CO AWB - EDGAR emitida'
   cdunits='kg/(m2 day)'

elseif(cvar(1:lv).eq.'sccobbe') then
   ivar_type=2
   ierr= RAMS_getvar('QSC5',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname='Massa de CO BB - EDGAR emitida'
   cdunits='kg/(m2 day)'

elseif(cvar(1:lv).eq.'sccod') then
   ivar_type=2
   ierr= RAMS_getvar('QSC9',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname='Massa de CO Duncan emitida'
   cdunits='kg/(m2 day)'

elseif(cvar(1:lv).eq.'sccol') then
   ivar_type=2
   ierr= RAMS_getvar('QSC3',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname='Massa de CO emitida -logan'
   cdunits='kg/(m2 day)'

elseif(cvar(1:lv).eq.'sccoant') then
   ivar_type=2
   ierr= RAMS_getvar('QSC9',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname='Massa de CO emitida -ANTROPO'
   cdunits='kg/(m2 day)'

elseif(cvar(1:lv).eq.'pwv') then
   ivar_type=2
   ierr= RAMS_getvar('RV',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
!air density
   ierr= RAMS_getvar('TOPT',idim_type,ngrd,e,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_dn0(n1,n2,n3,b,c,d,e,ngrd) ! d=dens_ar
   call RAMS_comp_mult(n1,n2,n3,a,d)         ! aqui a=rv*dens_ar
   call RAMS_comp_vertint(n1,n2,n3,a,e,ngrd) ! agua em kg/m^2
   call RAMS_comp_mults(n1,n2,n3,a,0.1) !converte para cm = 1 kg/m^2 * 100 cm/m / (1000 kg/m^3 dens_agua)
   cdname='precipitable water vapor'
   cdunits='cm'




elseif(cvar(1:lv).eq.'TKUO') then
   ivar_type=3
   ierr= RAMS_getvar('DUM3',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_transf_ppb_day(n1,n2,n3,a)
   cdname='CO tend conc due conv trans'
   cdunits='ppb/day'

elseif(cvar(1:lv).eq.'TKUOSH') then
   ivar_type=3
   ierr= RAMS_getvar('DUM8',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_transf_ppb_day(n1,n2,n3,a)
   cdname='CO tend conc due Shallow conv trans'
   cdunits='ppb/day'



! ------------------------ Stilt-RAMS coupling------------
elseif(cvar(1:lv).eq.'afxu') then
   ivar_type=3
   ierr= RAMS_getvar('AFXU',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname='adv u flux'
   cdunits='kg/m^2s'

elseif(cvar(1:lv).eq.'afxub') then
   ivar_type=3
   ierr_getvar = ierr_getvar + ierr
   ierr= RAMS_getvar('AFXUB',idim_type,ngrd,a,b,flnm)
   cdname='averaged adv u flux'
   cdunits='kg/m^2s'

elseif(cvar(1:lv).eq.'afxv') then
   ivar_type=3
   ierr= RAMS_getvar('AFXV',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname='adv v flux'
   cdunits='kg/m^2s'

elseif(cvar(1:lv).eq.'afxvb') then
   ivar_type=3
   ierr= RAMS_getvar('AFXVB',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname='averaged adv v flux'
   cdunits='kg/m^2s'

elseif(cvar(1:lv).eq.'afxw') then
   ivar_type=3
   ierr= RAMS_getvar('AFXW',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname='adv w flux'
   cdunits='kg/m^2s'

elseif(cvar(1:lv).eq.'afxwb') then
   ivar_type=3
   ierr= RAMS_getvar('AFXWB',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname='averaged adv W flux'
   cdunits='kg/m^2s'


elseif(cvar(1:lv).eq.'sigw') then
   ivar_type=3
   ierr= RAMS_getvar('SIGW',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname=' sigma W'
   cdunits=''

elseif(cvar(1:lv).eq.'sigwb') then
   ivar_type=3
   ierr= RAMS_getvar('SIGWB',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname='averaged sigma W'
   cdunits='m/s'

elseif(cvar(1:lv).eq.'tlb') then
   ivar_type=3
   ierr= RAMS_getvar('TLB',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname='averaged Lagr timescale'
   cdunits='s'

elseif(cvar(1:lv).eq.'tl') then
   ivar_type=3
   ierr= RAMS_getvar('TL',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname='Lagr timescale'
   cdunits='s'

elseif(cvar(1:lv).eq.'tkeb') then
   ivar_type=3
   ierr= RAMS_getvar('TKEPB',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_noneg(n1,n2,n3,a)
   cdname='average turb kinetic energy'
   cdunits='m2/s2'



!------------Grell cumulus scheme --------------------------

elseif(cvar(1:lv).eq.'wdm1') then
   ivar_type=2
   ierr= RAMS_getvar('wetdep001',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_mults(n1,n2,n3,a,1.e-6)  ! converte de mg/kg para kg/kg
   cdname='Wet deposition mass tracer 1'
   cdunits='kg/m2'


elseif(cvar(1:lv).eq.'wdm3') then
   ivar_type=2
   ierr= RAMS_getvar('wetdep003',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_mults(n1,n2,n3,a,1.e-6)  ! converte de mg/kg para kg/kg
   cdname='Wet deposition mass tracer 3'
   cdunits='kg/m2'


elseif(cvar(1:lv).eq.'cuprliq') then
   ivar_type=6
   ierr= RAMS_getvar('CUPRLIQ',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_noneg(n1,n2,n3*n6,a)
   call RAMS_comp_mults(n1,n2,n3*n6,a,1000.)
   call get_cumulus(n1,n2,n3,n6,a,a6)
   cdname='Conv. water mixing ratio'
   cdunits='g/kg'

elseif(cvar(1:lv).eq.'cuprice') then
   ivar_type=6
   ierr= RAMS_getvar('CUPRICE',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call RAMS_comp_noneg(n1,n2,n3*n6,a)
   call RAMS_comp_mults(n1,n2,n3*n6,a,1000.)
   call get_cumulus(n1,n2,n3,n6,a,a6)
   cdname='Conv. water mixing ratio'
   cdunits='g/kg'

elseif(cvar(1:lv).eq.'areadn') then
   ivar_type=9
   ierr= RAMS_getvar('AREADN',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname='Downdraft relative area'
   cdunits=''

elseif(cvar(1:lv).eq.'areaup') then
   ivar_type=9
   ierr= RAMS_getvar('AREAUP',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname='Updraft relative area'
   cdunits=''


elseif(cvar(1:lv).eq.'ierr') then
   ivar_type=9
   ierr= RAMS_getvar('XIERR',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname='Grell''s error flag:'
   cdunits=' '

elseif(cvar(1:lv).eq.'upmf') then
   ivar_type=9
   ierr= RAMS_getvar('UPMF',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname='updraft mass flux'
   cdunits='kg/m²/s'

elseif(cvar(1:lv).eq.'dnmf') then
   ivar_type=9
   ierr= RAMS_getvar('DNMF',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname='downdraft mass flux'
   cdunits='kg/m²/s'

elseif(cvar(1:lv).eq.'zktop') then
   ivar_type=9
   ierr= RAMS_getvar('ZKTOP',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname='Cloud top'
   cdunits='m'

elseif(cvar(1:lv).eq.'zkbcon') then
   ivar_type=9
   ierr= RAMS_getvar('ZKBCON',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname='Cloud bottom'
   cdunits='m'

elseif(cvar(1:lv).eq.'zjmin') then
   ivar_type=9
   ierr= RAMS_getvar('ZJMIN',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname='Downdraft originating height'
   cdunits='m'

elseif(cvar(1:lv).eq.'zk22') then
   ivar_type=9
   ierr= RAMS_getvar('ZK22',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname='Updraft originating height'
   cdunits='m'

elseif(cvar(1:lv).eq.'zkdet') then
   ivar_type=9
   ierr= RAMS_getvar('ZKDT',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname='Top of downdraft detrainment'
   cdunits='m'

elseif(cvar(1:lv).eq.'edt') then
   ivar_type=9
   ierr= RAMS_getvar('EDT',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname='Downdraft/updraft ratio'
   cdunits=' '

elseif(cvar(1:lv).eq.'aadn') then
   ivar_type=9
   ierr= RAMS_getvar('AADN',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname='Dndraft cloud work function'
   cdunits='J/kg'

elseif(cvar(1:lv).eq.'aaup') then
   ivar_type=9
   ierr= RAMS_getvar('AAUP',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname='Updraft cloud work function'
   cdunits='J/kg'


elseif(cvar(1:lv).eq.'xave') then
   ivar_type=2
   ierr= RAMS_getvar('DUM5',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call D3toD2(n1,n2,n3,2,a,c)
   cdname='X_AVE'
   cdunits=' '

elseif(cvar(1:lv).eq.'xavec1') then
   ivar_type=2
   ierr= RAMS_getvar('DUM5',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call D3toD2(n1,n2,n3,16,a,c)
   cdname='X_AVE Capmax'
   cdunits=' '

elseif(cvar(1:lv).eq.'xavec2') then
   ivar_type=2
   ierr= RAMS_getvar('DUM5',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call D3toD2(n1,n2,n3,17,a,c)
   cdname='X_AVE Capmax'
   cdunits=' '

elseif(cvar(1:lv).eq.'xavec3') then
   ivar_type=2
   ierr= RAMS_getvar('DUM5',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call D3toD2(n1,n2,n3,18,a,c)
   cdname='X_AVE Capmax'
   cdunits=' '

!****************

elseif(cvar(1:lv).eq.'xff0') then
   ivar_type=2
   ierr= RAMS_getvar('d2003',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname='XFF0 for deep'
   cdunits=' '
elseif(cvar(1:lv).eq.'xff0sh') then
   ivar_type=2
   ierr= RAMS_getvar('d2002',idim_type,ngrd,a,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   cdname='XFF0 for shallow'
   cdunits=' '

elseif(cvar(1:lv).eq.'prgr1') then
   ivar_type=2
   ierr= RAMS_getvar('d3004',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call D3toD2(n1,n2,n3,16,a,c)
   cdname=' precip closure 1 large cap'
   cdunits=' mm/h'

elseif(cvar(1:lv).eq.'prgr2') then
   ivar_type=2
   ierr= RAMS_getvar('d3004',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call D3toD2(n1,n2,n3,17,a,c)
   cdname=' precip closure 1 medium cap'
   cdunits=' mm/h'

elseif(cvar(1:lv).eq.'prgr3') then
   ivar_type=2
   ierr= RAMS_getvar('d3004',idim_type,ngrd,c,b,flnm)
   call D3toD2(n1,n2,n3,18,a,c)
   cdname=' precip closure 1 low cap'
   cdunits=' mm/h'

elseif(cvar(1:lv).eq.'prw1') then
   ivar_type=2
   ierr= RAMS_getvar('d3004',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call D3toD2(n1,n2,n3,19,a,c)
   cdname=' precip closure 2 large cap'
   cdunits=' mm/h'

elseif(cvar(1:lv).eq.'prw2') then
   ivar_type=2
   ierr= RAMS_getvar('d3004',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call D3toD2(n1,n2,n3,20,a,c)
   cdname=' precip closure 2 medium cap'
   cdunits=' mm/h'

elseif(cvar(1:lv).eq.'prw3') then
   ivar_type=2
   ierr= RAMS_getvar('d3004',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call D3toD2(n1,n2,n3,21,a,c)
   cdname=' precip closure 2 low cap'
   cdunits=' mm/h'

elseif(cvar(1:lv).eq.'prmc1') then
   ivar_type=2
   ierr= RAMS_getvar('d3004',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call D3toD2(n1,n2,n3,22,a,c)
   cdname=' precip closure 3 large cap'
   cdunits=' mm/h'

elseif(cvar(1:lv).eq.'prmc2') then
   ivar_type=2
   ierr= RAMS_getvar('d3004',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call D3toD2(n1,n2,n3,23,a,c)
   cdname=' precip closure 3 medium cap'
   cdunits=' mm/h'

elseif(cvar(1:lv).eq.'prmc3') then
   ivar_type=2
   ierr= RAMS_getvar('d3004',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call D3toD2(n1,n2,n3,24,a,c)
   cdname=' precip closure 3 low cap'
   cdunits=' mm/h'


elseif(cvar(1:lv).eq.'prst1') then
   ivar_type=2
   ierr= RAMS_getvar('d3004',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call D3toD2(n1,n2,n3,25,a,c)
   cdname=' precip closure 4 large cap'
   cdunits=' mm/h'

elseif(cvar(1:lv).eq.'prst2') then
   ivar_type=2
   ierr= RAMS_getvar('d3004',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call D3toD2(n1,n2,n3,26,a,c)
   cdname=' precip closure 4 medium cap'
   cdunits=' mm/h'

elseif(cvar(1:lv).eq.'prst3') then
   ivar_type=2
   ierr= RAMS_getvar('d3004',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call D3toD2(n1,n2,n3,27,a,c)
   cdname=' precip closure 4 low cap'
   cdunits=' mm/h'


elseif(cvar(1:lv).eq.'pras1') then
   ivar_type=2
   ierr= RAMS_getvar('d3004',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call D3toD2(n1,n2,n3,28,a,c)
   cdname=' precip closure 5 large cap'
   cdunits=' mm/h'

elseif(cvar(1:lv).eq.'pras2') then
   ivar_type=2
   ierr= RAMS_getvar('d3004',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call D3toD2(n1,n2,n3,29,a,c)
   cdname=' precip closure 5 medium cap'
   cdunits=' mm/h'

elseif(cvar(1:lv).eq.'pras3') then
   ivar_type=2
   ierr= RAMS_getvar('d3004',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call D3toD2(n1,n2,n3,30,a,c)
   cdname=' precip closure 5 low cap'
   cdunits=' mm/h'

   
elseif(cvar(1:lv).eq.'xstd') then
   ivar_type=2
   ierr= RAMS_getvar('DUM5',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call D3toD2(n1,n2,n3,3,a,c)
   cdname='X_STD'
   cdunits=' '

elseif(cvar(1:lv).eq.'xske') then
   ivar_type=2
   ierr= RAMS_getvar('DUM5',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call D3toD2(n1,n2,n3,4,a,c)
   cdname='x_ske'
   cdunits=' '

elseif(cvar(1:lv).eq.'xcur') then
   ivar_type=2
   ierr= RAMS_getvar('DUM5',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call D3toD2(n1,n2,n3,5,a,c)
   cdname='x_cur'
   cdunits=' '



elseif(cvar(1:lv).eq.'xmbgr') then
   ivar_type=2
   ierr= RAMS_getvar('DUM5',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call D3toD2(n1,n2,n3,6,a,c)
   cdname='xmbgr'
   cdunits=' '


elseif(cvar(1:lv).eq.'xmbw') then
   ivar_type=2
   ierr= RAMS_getvar('DUM5',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call D3toD2(n1,n2,n3,7,a,c)
   cdname='xmbw'
   cdunits=' '

elseif(cvar(1:lv).eq.'xmbmc') then
   ivar_type=2
   ierr= RAMS_getvar('DUM5',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call D3toD2(n1,n2,n3,8,a,c)
   cdname='xmbmc'
   cdunits=' '

elseif(cvar(1:lv).eq.'xmbst') then
   ivar_type=2
   ierr= RAMS_getvar('DUM5',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call D3toD2(n1,n2,n3,9,a,c)
   cdname='xmbst'
   cdunits=' '

elseif(cvar(1:lv).eq.'xmbas') then
   ivar_type=2
   ierr= RAMS_getvar('DUM5',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call D3toD2(n1,n2,n3,10,a,c)
   cdname='xmbas'
   cdunits=' '




elseif(cvar(1:lv).eq.'prgr') then
   ivar_type=2
   ierr= RAMS_getvar('DUM5',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call D3toD2(n1,n2,n3,11,a,c)
   cdname='prgr'
   cdunits=' '


elseif(cvar(1:lv).eq.'prw') then
   ivar_type=2
   ierr= RAMS_getvar('DUM5',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call D3toD2(n1,n2,n3,12,a,c)
   cdname='prw'
   cdunits=' '

elseif(cvar(1:lv).eq.'prmc') then
   ivar_type=2
   ierr= RAMS_getvar('DUM5',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call D3toD2(n1,n2,n3,13,a,c)
   cdname='prmc'
   cdunits=' '

elseif(cvar(1:lv).eq.'prst') then
   ivar_type=2
   ierr= RAMS_getvar('DUM5',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call D3toD2(n1,n2,n3,14,a,c)
   cdname='prst'
   cdunits=' '

elseif(cvar(1:lv).eq.'pras') then
   ivar_type=2
   ierr= RAMS_getvar('DUM5',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call D3toD2(n1,n2,n3,15,a,c)
   cdname='pras'
   cdunits=' '

elseif(cvar(1:lv).eq.'um') then
   ivar_type=2
   ierr= RAMS_getvar('DUM5',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call D3toD2(n1,n2,n3,24,a,c)
   cdname='u mean'
   cdunits='m/s '

elseif(cvar(1:lv).eq.'vm') then
   ivar_type=2
   ierr= RAMS_getvar('DUM5',idim_type,ngrd,c,b,flnm)
   ierr_getvar = ierr_getvar + ierr
   call D3toD2(n1,n2,n3,25,a,c)
   cdname='v mean'
   cdunits='m/s '



else

   print*,'Variable name not found in hvlib.f - ',cvar(1:lv)
   ivar_type=0

endif

if(ierr_getvar > 0) then
  print *, 'Not all variables needed for computing '//cvar(1:lv)//' are available...' 
  ivar_type=0
end if
return
end

!*******************************************************************************
!############################# Change Log ##################################
! 2.3.0.1
!
! 000830 CJT rams_comp ##
!            Corrected reference to "cpr" and "r" to "cpor" and "rgas" ##
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!  Mission Research Corporation / *ASTeR Division
!###########################################################################

subroutine RAMS_comp(n1,n2,n3,n4,n5,n6)

use somevars
use rconstants
use therm_lib, only : qtk, qwtk, dewfrostpoint, rslif, virtt, thetaeiv

dimension a(n1,n2,n3),b(n1,n2,n3),c(n1,n2,n3),d(n1,n2,n3),e(n1,n2,n3),o(n1,n2,n3),topt(n1,n2)
dimension a2(n1,n2,n4,n5),a6(n1,n2,n3,n6)
real f1(n1,n2,n3),f2(n1,n2,n3)
dimension theta(n1,n2,n3),pp(n1,n2,n3),slp(n1,n2),z(n1,n2,n3)
include 'rcommons.h'
dimension slmsts0(12)
data slmsts0/0.395, 0.410, 0.435, 0.485, 0.451, 0.420  &
            ,0.477, 0.476, 0.426, 0.492, 0.482, 0.863/
real, parameter, dimension(12) :: myslcpd=(/1.465e6, 1.407e6, 1.344e6, 1.273e6, 1.214e6, 1.177e6, &
                                            1.319e6, 1.227e6, 1.177e6, 1.151e6, 1.088e6,  .874e6 /)
real :: temptemp,fracliq
integer :: nsoil
!SRF
!  PMAR, PMCO, PMC02 = pesos moleculares do ar, CO, CO2
data PMAR/28.96/
data PMCO/28./
data PMCO2/44./

entry RAMS_transf_ppb_day(n1,n2,n3,a)
!  PMAR e PMCO = pesos moleculares do ar e do CO
!  TRANSFORMACAO DE kg[CO]/kg[AR] para ppb (PARTE POR BILHAO)
     do k=1,n3
 	do j=1,n2
 	   do i=1,n1

     a(i,j,k)=a(i,j,k)*(PMAR/PMCO)*1.E+9*1.e-6*86400.

 	   enddo
 	enddo
     enddo
return

entry RAMS_transf_ppb(n1,n2,n3,a)
!  TRANSFORMACAO DE kg[CO]/kg[AR] para ppb (PARTE POR BILHAO)
     do k=1,n3
 	do j=1,n2
 	   do i=1,n1

     a(i,j,k)=a(i,j,k)*(PMAR/PMCO)*1.E+9

 	   enddo
 	enddo
     enddo
return

entry RAMS_transf_ppm(n1,n2,n3,a)
!  TRANSFORMACAO DE kg[CO2]/kg[AR] para ppm (PARTE POR MILHAO)
     do k=1,n3
 	do j=1,n2
 	   do i=1,n1
 	      a(i,j,k)=a(i,j,k)*(PMAR/PMCO2)*1.E+6
 	   enddo
 	enddo
     enddo
return

entry RAMS_transf_ugm3(n1,n2,n3,a,b)
!   Transformacao de conc  de kg/kg para ug/m3
    do k=1,n3
      do j=1,n2
    	do i=1,n1
        !print*,a(i,j,k),b(i,j,k)
	a(i,j,k)=a(i,j,k)*b(i,j,k)*1.e+9
    
    	enddo
      enddo
    enddo
return

entry get_ZI(n1,n2,n3,a,c,ngrd)
tkemin =1.E-2
      
      do i=1,n1
         do j=1,n2
          a(i,j,1) = (ztn(2,ngrd)+ ztn(1,ngrd))
         enddo
      enddo

      do i=1,n1
        do j=1,n2
         do k = 2,n3-1
!
!         print*,k,c(i,j,k),ztn(k,ngrd)
         if(c(i,j,k).lt.tkemin) then
            kzi = k
            a(i,j,1)=0.5*(ztn(kzi,ngrd)+ ztn(kzi-1,ngrd))
	    go to 500
	 endif
	 enddo
 500     continue
!         print*,i,j,a(i,j,1)	  
        enddo
      enddo
return

entry RAMS_comp_nebulosidade(n1,n2,n3,a,b,c,e,ngrd)
!
! b=cloud
! c=dn0
! e=topo
       do j=1,n2
           do i=1,n1
         rnebu=0.
         rodzint=0
         do k=2,n3-1
          dz=(ztn(k,ngrd)-ztn(k-1,ngrd))*(1.-e(i,j,1)/zmn(nnzp(1)-1,1))
          rnebu   = rnebu   + b(i,j,k)*c(i,j,k)*dz	  
          rodzint = rodzint +	       c(i,j,k)*dz	  
         enddo
              a(i,j,1) = 1000.*rnebu/(rodzint+1.e-5)
	 enddo
       enddo
return

!entry RAMS_get_leaf(n1,n2,n5,a)
!       do j=1,n2
!        do i=1,n1
!         do ip=1,n5
!            a(i,j,ip)=a(i,j,ip)
!         enddo
!        enddo
!       enddo
!return

entry RAMS_comp_vegclass2(n1,n2,n3,a)
      do i=1,n1
      a(i,1,1) = a(i,1,1) + .5
      enddo
return

!entry RAMS_comp_vegclass(n1,n2,n5,a)
!       do j=1,n2
!        do i=1,n1
!         do ip=1,n5
!!            print*,i,j,ip,' veg=',a(i,j,ip)
!            a(i,j,ip)=float(int(a(i,j,ip)+0.5))
!!            print*,i,j,ip,' veg=',a(i,j,ip)
!         enddo
!        enddo
!       enddo
!return

entry D3toD2(n1,n2,n3,klevel,a,c)
  do j=1,n2
      do i=1,n1
	a(i,j,1)=c(i,j,klevel)
!	print*,i,j,klevel,c(i,j,klevel)
      enddo
  enddo
return


entry up_to_tp(n1,n2,n3,a,b)
         do k=1,n3
            do j=1,n2
               a(1,j,k)=b(1,j,k)
               do i=2,n1
               a(i,j,k)=0.5* (b(i,j,k) + b(i-1,j,k))
               enddo
            enddo
         enddo
return

entry vp_to_tp(n1,n2,n3,a,b)
         do k=1,n3
            do i=1,n1
               a(i,1,k)=b(i,1,k)
               do j=2,n2
               a(i,j,k)=0.5* (b(i,j,k) + b(i,j-1,k))
               enddo
            enddo
         enddo
return

entry wp_to_tp(n1,n2,n3,a,b)
         do j=1,n2
           do i=1,n1
              a(i,j,1)=b(i,j,1)
              do k=2,n3
               a(i,j,k)=0.5* (b(i,j,k) + b(i,j,k-1))
              enddo
            enddo
         enddo
return

entry set_undef(n1,n2,n3,a,c)
  do j=1,n2
      do i=1,n1
	a(i,j,1)=c(i,j,1)
!        if(a(i,j,1) .lt. 0.0001) a(i,j,1) = -9.99e33
!	print*,i,j,klevel,c(i,j,klevel)
      enddo
  enddo
return
!SRF
!***********************************************************

entry RAMS_comp_vals(n1,n2,n3,a)
   print*,'==================== values =========================='
   do k=1,n3
      do j=1,n2
         do i=1,n1
            if(a(i,j,k).ne.0.) write(*,'(3i3,e14.6)') i,j,k  &
                                                     ,a(i,j,k)
         enddo
      enddo
   enddo
   print*,'======================================================'
return

entry RAMS_comp_tot(n1,n2,n3,a)
   tot=0.
   do k=1,n3
      do j=1,n2
         do i=1,n1
            tot=tot+a(i,j,k)
         enddo
      enddo
   enddo
   write(*,'(a,e12.6)') '-> total- ',tot
return

entry RAMS_comp_maxval(n1,n2,n3,a)
   zmax=-1.0e30
   do k=1,n3
      do j=1,n2
         do i=1,n1
            if(a(i,j,k).gt.zmax) then
               zmax=a(i,j,k)
               maxx=i
               maxy=j
               maxz=k
            endif
         enddo
      enddo
   enddo
   write(*,'(a,e12.6,a,3i3)') '-> max- ',zmax,' at i,j,k-',maxx  &
                              ,maxy,maxz
return

entry RAMS_comp_minval(n1,n2,n3,a)
   zmin=1.0e30
   do k=1,n3
      do j=1,n2
         do i=1,n1
            if(a(i,j,k).lt.zmin) then
               zmin=a(i,j,k)
               minx=i
               miny=j
               minz=k
            endif
         enddo
      enddo
   enddo
   write(*,'(a,e12.6,a,3i3)') '-> min- ',zmin,' at i,j,k-',minx  &
                              ,miny,minz
return

entry RAMS_comp_zero(n1,n2,n3,a)
   do k=1,n3
      do j=1,n2
         do i=1,n1
            a(i,j,k)=0.
         enddo
      enddo
   enddo
return

entry RAMS_comp_1minus(n1,n2,n3,a)
   do k=1,n3
      do j=1,n2
         do i=1,n1
            a(i,j,k)=1.-a(i,j,k)
         enddo
      enddo
   enddo
return

entry RAMS_comp_mults(n1,n2,n3,a,s)
   do k=1,n3
      do j=1,n2
         do i=1,n1
            a(i,j,k)=a(i,j,k) * s
         enddo
      enddo
   enddo
return


entry RAMS_comp_accum(n1,n2,n3,a,b)
   do k=1,n3
      do j=1,n2
         do i=1,n1
            a(i,j,k)=a(i,j,k)+b(i,j,k)
         enddo
      enddo
   enddo
return

entry RAMS_comp_noneg(n1,n2,n3,a)
   do k=1,n3
      do j=1,n2
         do i=1,n1
            a(i,j,k)=max(a(i,j,k),0.)
         enddo
      enddo
   enddo
return

entry RAMS_comp_nopos(n1,n2,n3,a)
   do k=1,n3
      do j=1,n2
         do i=1,n1
            a(i,j,k)=min(a(i,j,k),0.)
         enddo
      enddo
   enddo
return

entry RAMS_comp_subt(n1,n2,n3,a,b)
   do k=1,n3
      do j=1,n2
         do i=1,n1
            a(i,j,k)=a(i,j,k)-b(i,j,k)
         enddo
      enddo
   enddo
return

entry RAMS_comp_mult(n1,n2,n3,a,b)
   do k=1,n3
      do j=1,n2
         do i=1,n1
            a(i,j,k)=a(i,j,k)*b(i,j,k)
         enddo
      enddo
   enddo
return

!Demerval [
entry RAMS_comp_pw(n1,n2,n3,a,e,c,ngrd)

  do j=1,n2
    do i=1,n1
       a(i,j,1)=0.
    enddo
  enddo  
  do j=1,n2
     do i=1,n1
       do k=2,n3
            a(i,j,1)=a(i,j,1) + e(i,j,k)*c(i,j,k)*    &
             (ztn(k,ngrd)-ztn(k-1,ngrd))*             &
             (1.-c(i,j,k)/zmn(nnzp(1)-1,1))*0.0001
       enddo
     enddo
  enddo          

return
!Demerval ]

entry RAMS_comp_z(n1,n2,n3,a,c,ngrd)
   do k=1,n3
      do j=1,n2
         do i=1,n1
            a(i,j,k)=c(i,j,1)  &
                 +ztn(k,ngrd)*(1.-c(i,j,1)/zmn(nnzp(1)-1,1))
         enddo
      enddo
   enddo
return

entry RAMS_comp_rotate(n1,n2,n3,a,b,ngrd)
   do k=1,n3
      do j=1,n2
         do i=1,n1
            call xy_ll(qlat,qlon,platn(ngrd),plonn(ngrd)  &
               ,xtn(i,ngrd),ytn(j,ngrd))
            u=a(i,j,k)
            v=b(i,j,k)
            call uvtoueve(u,v,a(i,j,k),b(i,j,k)  &
                         ,qlat,qlon,platn(ngrd),plonn(ngrd))
         enddo
      enddo
   enddo
return

entry RAMS_comp_tempK(n1,n2,n3,a,b)
   do k=1,n3
      do j=1,n2
         do i=1,n1
            a(i,j,k)=a(i,j,k)*b(i,j,k)/cp
         enddo
      enddo
   enddo
return

entry RAMS_comp_press(n1,n2,n3,a)
   do k=1,n3
      do j=1,n2
         do i=1,n1
            a(i,j,k)=(a(i,j,k)/cp)**cpor*p00*.01
         enddo
      enddo
   enddo
return

entry RAMS_comp_tempC(n1,n2,n3,a)
   do k=1,n3
      do j=1,n2
         do i=1,n1
            a(i,j,k)=a(i,j,k)-t00
         enddo
      enddo
   enddo
return

entry RAMS_comp_theta2temp(n1,n2,n3,a,b)
   do k=1,n3
      do j=1,n2
         do i=1,n1
            a(i,j,k)=a(i,j,k) * (p00i * b(i,j,k)) ** rocp
         enddo
      enddo
   enddo
return

entry RAMS_comp_tempF(n1,n2,n3,a)
   do k=1,n3
      do j=1,n2
         do i=1,n1
            a(i,j,k)=(a(i,j,k)-273.16)*1.8+32.
         enddo
      enddo
   enddo
return

entry RAMS_comp_wcms(n1,n2,n3,a)
   do k=1,n3
      do j=1,n2
         do i=1,n1
            a(i,j,k)=a(i,j,k)*100.
         enddo
      enddo
   enddo
return

entry RAMS_comp_avgw(n1,n2,n3,a)
   do k=n3,2,-1
      do j=1,n2
         do i=1,n1
            a(i,j,k)=0.5*(a(i,j,k)+a(i,j,k-1))
         enddo
      enddo
   enddo
return

entry RAMS_comp_avgu(n1,n2,n3,a)
   do k=1,n3
      do j=1,n2
         do i=n1,2,-1
            a(i,j,k)=0.5*(a(i,j,k)+a(i-1,j,k))
         enddo
      enddo
   enddo
return

entry RAMS_comp_avgv(n1,n2,n3,a)
   do k=1,n3
      do j=n2,2,-1
         do i=1,n1
            a(i,j,k)=0.5*(a(i,j,k)+a(i,j-1,k))
         enddo
      enddo
   enddo
return

entry RAMS_comp_sfcdiv(n1,n2,n3,a,ngrd)
   do j=1,n2
      do i=1,n1
         a(i,j,1)=-(a(i,j,2)-a(i,j,1))*dztn(2,ngrd)
      enddo
   enddo
return

entry RAMS_comp_rt(n1,n2,n3,a)
   do k=1,n3
      do j=1,n2
         do i=1,n1
            a(i,j,k)=max(0.,a(i,j,k))*1000.
         enddo
      enddo
   enddo
return

entry RAMS_comp_hr_pcprate(n1,n2,n3,a,b,c)
   do k=1,n3
      do j=1,n2
         do i=1,n1
            a(i,j,k)=(b(i,j,k)-c(i,j,k))
         enddo
      enddo
   enddo
return

entry RAMS_comp_speed(n1,n2,n3,a,b)
   do k=1,n3   
      do j=1,n2
         do i=1,n1
            a(i,j,k)=sqrt(a(i,j,k)**2+b(i,j,k)**2)
         enddo
      enddo
   enddo
return

entry RAMS_comp_dir(n1,n2,n3,a,b,ngrd)
   do k=1,n3
      do j=1,n2
         do i=1,n1
            call xy_ll(qlat,qlon,platn(ngrd),plonn(ngrd)  &
               ,xtn(i,ngrd),ytn(j,ngrd))
            u=a(i,j,k)
            v=b(i,j,k)
            call uvtoueve(u,v,a(i,j,k),b(i,j,k)  &
                         ,qlat,qlon,platn(ngrd),plonn(ngrd))
            call winddf(a(i,j,k),ff,a(i,j,k),b(i,j,k))
         enddo
      enddo
   enddo
return

entry RAMS_comp_dewK(n1,n2,n3,a,b,c)
   do k=1,n3
      do j=1,n2
         do i=1,n1
            xpress=(b(i,j,k)/cp)**cpor*p00
            xtemp=c(i,j,k)*b(i,j,k)/cp
            xwatsat=rslif(xpress,xtemp)
            a(i,j,k)=dewfrostpoint(xpress,min(a(i,j,k),xwatsat) )
         enddo
      enddo
   enddo
return

!<Demerval
entry RAMS_comp_thete(n1,n2,n3,a,e,f1,f2)
   do k=1,n3
      do j=1,n2
         do i=1,n1
         TL=55.+1./( 1./(f1(i,j,k)-55.) -  &
            log(f2(i,j,k)/100.)/2840.)
            a(i,j,k)=a(i,j,k)*exp((alvl*e(i,j,k))/(cp*TL))
         enddo
      enddo
   enddo
return

!entry RAMS_comp_thete(n1,n2,n3,a,b,c)
!   do k=1,n3
!      do j=1,n2
!         do i=1,n1
!            xpress=(b(i,j,k)/cp)**cpor*p00
!            xtemp=c(i,j,k)*b(i,j,k)/cp
!            xwatsat=rslif(xpress,xtemp)
!            a(i,j,k)=c(i,j,k)*exp( alvl*xwatsat  &
!                 /(cp*dewfrostpoint(xpress,min(a(i,j,k),xwatsat) )) )
!         enddo
!      enddo
!   enddo
!return

!Demerval>

entry RAMS_comp_thetv(n1,n2,n3,a,b,c)
   do k=1,n3
      do j=1,n2
         do i=1,n1
            a(i,j,k)=virtt(a(i,j,k),b(i,j,k),c(i,j,k))
         enddo
      enddo
   enddo
return

entry RAMS_comp_thetaeiv(n1,n2,n3,a,b,c,d,e)
   do k=1,n3
      do j=1,n2
         do i=1,n1
            a(i,j,k)=thetaeiv(a(i,j,k),c(i,j,k),b(i,j,k),d(i,j,k),e(i,j,k),12,.true.)
         enddo
      enddo
   enddo
return

entry RAMS_comp_bowen(n1,n2,n3,a,b)
   do k=1,n3
      do j=1,n2
         do i=1,n1
            a(i,j,k)=a(i,j,k)/max(1.e-12,b(i,j,k))*1004./2.5e6
         enddo
      enddo
   enddo
return

entry RAMS_comp_rh(n1,n2,n3,a,b,c)
   do k=1,n3
      do j=1,n2
         do i=1,n1
            xtemp=c(i,j,k)*b(i,j,k)/cp
            xpress=(b(i,j,k)/cp)**cpor*p00
            a(i,j,k)=100.*min(1.  &
                 ,max(0.,a(i,j,k)/rslif(xpress,xtemp)))
         enddo
      enddo
   enddo
return

entry RAMS_comp_watsat(n1,n2,n3,a,b,c)
   do k=1,n3
      do j=1,n2
         do i=1,n1
            c(i,j,k)=c(i,j,k)*a(i,j,k)/cp
            b(i,j,k)=(b(i,j,k)/cp)**cpor*p00
            a(i,j,k)=rslif(b(i,j,k),c(i,j,k))
         enddo
      enddo
   enddo
return

entry RAMS_comp_pvap(n1,n2,n3,a,b) 
   ! a is pressure in Pa or hPa
   ! b is mixing ratio in kg/kg
   do k=1,n3
      do j=1,n2
         do i=1,n1
            b(i,j,k)=b(i,j,k)*a(i,j,k)/(ep+b(i,j,k))
         enddo
      enddo
   enddo
return
   
entry RAMS_comp_spvol(n1,n2,n3,a,b,c)
!a is temperature [K]
!b is vapour pressure in [Pa]
!c is pressure [Pa]
   do k=1,n3
      do j=1,n2
         do i=1,n1
            a(i,j,k)=rgas*a(i,j,k)/(c(i,j,k)-(1-ep)*b(i,j,k))
         enddo
      enddo
   enddo
return


entry RAMS_comp_vegclass(n1,n2,n3,a)
   do i=1,n1
      a(i,1,1) = float(nint(a(i,1,1)))
!      print*,i,int(a(i,1,1)) 
   enddo
   return

entry RAMS_comp_horizdiv(n1,n2,n3,a)
   do k=2,n3
      do j=1,n2
         do i=1,n1
            a(i,j,k)=-(a(i,j,k)-a(i,j,k-1))*dztn(k,ngrd)
         enddo
      enddo
   enddo
return

entry RAMS_comp_vertint(n1,n2,n3,a,topt,ngrd)
   ztop = zmn(nnzp(1)-1,1)
   do j = 1,n2
      do i = 1,n1
         rtgt = 1. - topt(i,j) / ztop
         a(i,j,1) = 0.
         do k = 2,n3-1
            a(i,j,1) = a(i,j,1) + a(i,j,k) * (zmn(k,ngrd)-zmn(k-1,ngrd)) * rtgt
         enddo
      enddo
   enddo
return

entry RAMS_comp_ppress(n1,n2,n3,a,c)
   do k=1,n3
      do j=1,n2
         do i=1,n1
            a(i,j,k) = 1000. * (a(i,j,k)/cp) ** cpor  &
                     - 1000. * (c(i,j,k)/cp) ** cpor
         enddo
      enddo
   enddo
return

entry RAMS_comp_raintemp(n1,n2,n3,a)
   do k=1,n3
      do j=1,n2
         do i=1,n1
            if (a(i,j,k) > 0.) then
               a(i,j,k) = tsupercool + a(i,j,k) * cliqi - t00
            else
               a(i,j,k) = -9.99e33
            end if
         enddo
      enddo
   enddo
return

entry RAMS_comp_qwtc(n1,n2,n3,a,b,c)

!MLO: Input: a = Energy  in J/m3 so we don't need to multiply...
!            b = Water      m3/m3
!            c = Soil class 
!     Output: a = Temperature in Celsius
!             b = Fraction in liquid phase

   do k=1,n3
      do j=1,n2
         do i=1,n1

            nsoil = nint(c(i,j,k))
            dryhcap = myslcpd(nsoil)
            b(i,j,k) = b(i,j,k) * 1.e3
            temptemp = 6.0
            fracliq = 4.0
            call qwtk(a(i,j,k),b(i,j,k),dryhcap,temptemp,fracliq)
            !write (unit=*,fmt='(5(a,1x,es12.5,1x))') &
            !  'QW=',a(i,j,k),'W=',b(i,j,k),'SCP=',dryhcap,'TEMPK=',temptemp,'FLIQ=',fracliq
            a(i,j,k) = temptemp-t00
            b(i,j,k) = fracliq

         end do
      end do
   end do
return

entry RAMS_comp_copysst(n1,n2,n3,a)
   do k=1,n3
      do j=1,n2
         do i=1,n1
            call qtk(a(i,j,n3),temptemp,fracliq)
            a(i,j,k) = temptemp - t00
         enddo
      enddo
   enddo
return

entry RAMS_comp_qtcpcp(n1,n2,n3,a)
   do k=1,n3
      do j=1,n2
         do i=1,n1
            if (a(i,j,k) > 0.) then
               call qtk(a(i,j,k),temptemp,fracliq)
               a(i,j,k) = temptemp - t00
            else
               a(i,j,k) = -9.99e33
            end if
         end do
      end do
   end do
return

entry RAMS_comp_fracliq(n1,n2,n3,a)
   do k=1,n3
      do j=1,n2
         do i=1,n1
            call qtk(a(i,j,k),temptemp,fracliq)
            a(i,j,k) = fracliq
         end do
      end do
   enddo
return

entry RAMS_comp_fracice(n1,n2,n3,a)
   do k=1,n3
      do j=1,n2
         do i=1,n1
            call qtk(a(i,j,k),temptemp,fracliq)
            a(i,j,k) = 1.0 - fracliq
         enddo
      enddo
   enddo
return

entry RAMS_comp_hydrodiam(n1,n2,n3,a,c,ccfmas,ppwmas)
   rpwmas = 1. / ppwmas
   do k=1,n3
      do j=1,n2
         do i=1,n1
            if(a(i,j,k) .gt. 1.e-10 .and. c(i,j,k).gt.1.e-10)then
               a(i,j,k) = (a(i,j,k) / (c(i,j,k) * ccfmas))**rpwmas
            else
               a(i,j,k) = 0.
            endif
         enddo
      enddo
   enddo
return

entry RAMS_comp_slmstf(n1,n2,n3,a,c)
   do i=1,n1
     do j=1,n2
        do k=1,n3
           nsoil = nint(c(i,j,k))
           if (nsoil /= 0) then
              a(i,j,k) = a(i,j,k) / max(1.e-6,slmsts0(nsoil))
           else
              a(i,j,k) = 1.
           end if
        end do
     end do
   enddo
return

entry rams_sum_snowlayers(n1,n2,n3,a)
   do ip=1,n3
      do k=2,n2
         do ij=1,n1
            a(ij,1,ip) = a(ij,1,ip) + a(ij,k,ip)
         enddo
      enddo
   enddo
return

entry rams_fill_sst(n1,n2,n3,kp,a,c)
   do j=1,n2
      do i = 1,n1
         call qtk(c(i,j,kp),temptemp,fracliq)
         a(i,j,1) = temptemp-t00
      enddo
   enddo
return

entry rams_comp_pcpgnorm(n1,n2,n3,a,c)
return

entry rams_comp_vapnorm(n1,n2,n3,a,c)
return

entry rams_comp_snownorm(n1,n2,n3,a,c)
return

entry rams_comp_vegnorm(n1,n2,n3,a,c)
return

entry rams_comp_cannorm(n1,n2,n3,a,c)
return


entry acha_ztropop(n1,n2,n3,a,c,ngrd)
  a(1:n1,1:n2,1)=0.
  estratosfera=0.018 ! 15 K/Km
  do i=1,n1
    do j=1,n2
      malhakz: do k=3,n3
        dtheta  = (c(i,j,k)-c(i,j,k-1))/(ztn(k,ngrd)-ztn(k-1,ngrd))
        if (dtheta.gt.estratosfera) exit malhakz
      end do malhakz
      a(i,j,1)=0.5*(ztn(k,ngrd)+ztn(k-1,ngrd))
    end do
  end do
return

entry acha_ttropop(n1,n2,n3,a,c,e,ngrd)
  a(1:n1,1:n2,1)=0.
  estratosfera=0.018 ! 15 K/Km
  do i=1,n1
    do j=1,n2
      malhakt: do k=3,n3
        dtheta  = (c(i,j,k)-c(i,j,k-1))/(ztn(k,ngrd)-ztn(k-1,ngrd))
        if (dtheta.gt.estratosfera) exit malhakt
      end do malhakt
      a(i,j,1)=0.5*(e(i,j,k)+e(i,j,k-1))
    end do
  end do
return

entry acha_ptropop(n1,n2,n3,a,c,e,ngrd)
  a(1:n1,1:n2,1)=0.
  estratosfera=0.018 ! 15 K/Km
  do i=1,n1
    do j=1,n2
      malhakp: do k=3,n3
        dtheta  = (c(i,j,k)-c(i,j,k-1))/(ztn(k,ngrd)-ztn(k-1,ngrd))
        if (dtheta.gt.estratosfera) exit malhakp
      end do malhakp
      a(i,j,1)=exp(0.5*(log(e(i,j,k))+log(e(i,j,k-1))))
    end do
  end do
return


entry get_ZItheta(n1,n2,n3,a,c,e,ngrd)
      do i=1,n1
        do j=1,n2
          a(i,j,1) = 0.
        enddo
      enddo

      do i=1,n1
        do j=1,n2

          do k=3,n3-5            
            dtheta=c(i,j,k)-c(i,j,k-1)
            x=0.6
            if(dtheta.gt.x.or.abs(e(i,j,k)).gt.1.e-5) go to 1881
          enddo
 1881     continue
          kzi = k
          if(abs(e(i,j,k)).gt.1.e-5) then

            a(i,j,1) =0.5*(ztn(kzi,ngrd)+ ztn(kzi-1,ngrd))

          else
        
            a(i,j,1) =0.5*(ztn(kzi,ngrd)+ ztn(kzi-1,ngrd))
          endif

          if(a(i,j,1).lt.0..or.kzi.le.2) a(i,j,1)=0.
       enddo
     enddo

   return



entry RAMS_comp_pbl (n1,n2,n3,a,c,ngrd)

tkethrsh=0.001   ! tke threshold for PBL height in m2/s2
do j=1,n2
   do i=1,n1
      do k=1,n3
        a(i,j,k)=0.
      enddo
   enddo
enddo

do j=1,n2
   do i=1,n1
      pblht=0.
      do k=2,n3
         pblht=ztn(k,ngrd)*(1.-c(i,j,1)/zmn(nnzp(1)-1,1))
!DSM         if(a(i,j,k).le.tkethrsh) goto 10
      enddo
      10 continue
      do k=1,n3
        a(i,j,k)=pblht
      enddo
   enddo
enddo

return

entry RAMS_comp_etrans(n1,n2,n3,a,b,a2d)
   do j=1,n2
      do i=1,n1
         temp1=a(i,j,1)*b(i,j,1)/cp
         press1=(b(i,j,1)/cp)**cpor*p00
         dens=press1/(rgas*temp1)
         if(i.eq.5.and.j.eq.5) then
            print*,'============++++++'
            print*,temp1,press1,dens,a2d(i,j)
         endif
         a(i,j,1)=a2d(i,j)*dens*1.e-3*39.37*3600.
         do k=2,n3
            a(i,j,k)=a(i,j,1)
         enddo
      enddo
   enddo
return

entry RAMS_comp_slpress(n1,n2,n3,theta,pp,z,slp)
!
!     This subroutine calculates the pressure at level zlev. it
!       is hardwired here to calculate mean sea level pressure,
!       but can be easily changed to calculate pressure at any level
!       by just changing zlev.
!     a standard atmosphere lapse rate of 6.5 c/km is used to
!       interpolate temperature down from level 2 in the model.
!

!c      rlap=-.0065  ! standard temp lapse rate
   rlap=.0025     ! approx standard theta lapse rate
   zlev=0.

   do j=1,n2
      do i=1,n1
         levloop: do k=2,n3
            if(z(i,j,k).ge.zlev) then
               ktop=k
               kbot=k-1
               exit levloop
            endif
         end do levloop

         ddz=zlev-z(i,j,kbot)
         if(zlev.lt.z(i,j,kbot))then
            thbar=(theta(i,j,kbot)-.5*ddz*rlap)
         else
            thbar=.5*(theta(i,j,kbot)+theta(i,j,ktop))
         endif
         slp(i,j)=pp(i,j,kbot)-ddz*sl_g/thbar
         slp(i,j)=(slp(i,j) * cpi)**cpor*p00
      end do
   end do
return

entry RAMS_comp_ctprof(n1,n2,n3,a,b,ngrd)
   do i=1,n1
      do j=1,n2
         kmax=0
         do k=1,n3
            if(a(i,j,k).ge.0.0001.and.b(i,j,k).ge.0.99)kmax=k
         enddo
         if(kmax.gt.2)then
            a(i,j,1)=ztn(kmax,ngrd)
         else
            a(i,j,1)=0.0
         endif
      enddo
   enddo
return

end


subroutine RAMS_comp_multap(n1,n2,n3,n4,a,b)
dimension a(n1,n2,n4),b(n1,n2,n3)
   do k=1,n4
      do j=1,n2
         do i=1,n1
            a(i,j,k)=a(i,j,k)*b(i,j,1)
         enddo
      enddo
   enddo
return
end
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine get_leaf_soil(n1,n2,n3,n4,n5,a,a2)
   implicit none
   integer, intent(in) :: n1,n2,n3,n4,n5
   real, dimension(n1,n2,n4,n5), intent(out) :: a2
   real, dimension(n1,n2,n4*n5), intent(in) :: a
   integer :: kip, k,i,j,ip
   kip=0
   do ip=1,n5
      do k=1,n4
         kip=kip+1
         do i=1,n1
            do j=1,n2
               a2(i,j,k,ip)=a(i,j,kip)
            end do
         end do
      end do
   end do
   return
end subroutine get_leaf_soil
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine get_cumulus(n1,n2,n3,n6,a,a6)
   implicit none
   integer, intent(in) :: n1,n2,n3,n6
   real, dimension(n1,n2,n3,n6), intent(out) :: a6
   real, dimension(n1,n2,n3*n6), intent(in) :: a
   integer :: kip, k,i,j,ip

   kip=0
   do ip=1,n6
      do k=1,n3
         kip=kip+1
         do i=1,n1
            do j=1,n2
               a6(i,j,k,ip)=a(i,j,kip)
            end do
         end do
      end do
   end do
   return
end subroutine get_cumulus
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine RAMS_comp_richardson(n1,n2,n3,np,rib,z0,speed,thetav_atm,thetav_can,topt,ngrd)

   use somevars
   use rconstants
   use therm_lib, only : qtk, qwtk, dewfrostpoint, rslif, virtt, thetaeiv
   implicit none
   integer                     , intent(in)    :: n1,n2,n3,np
   real   , dimension(n1,n2,np), intent(inout) :: rib
   real   , dimension(n1,n2,np), intent(in)    :: z0,thetav_can
   real   , dimension(n1,n2,n3), intent(in)    :: speed,thetav_atm
   real   , dimension(n1,n2)   , intent(in)    :: topt
   integer                     , intent(in)    :: ngrd
   integer                                     :: i,j,p
   real                                        :: zedtop,zagl,spd

   zedtop = myzmn(mynnzp(1)-1,1)
   do j=1,n2
      do i=1,n1
         zagl=myztn(2,ngrd)*(1.-topt(i,j)/zedtop)
         do p=1,np
            spd = max(speed(i,j,2),0.65)
            rib(i,j,p) = grav * (zagl-z0(i,j,p)) * (thetav_atm(i,j,2)-thetav_can(i,j,p))     &
                       / (0.5 * (thetav_can(i,j,2)+thetav_can(i,j,p)) * spd**2)
         end do
      end do
   end do
   return
end subroutine RAMS_comp_richardson
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine RAMS_comp_dn0(n1,n2,n3,a,b,c,topt,ngrd)

   use somevars
   use rconstants
   use therm_lib, only : qtk, qwtk, dewfrostpoint, rslif, virtt, thetaeiv
   implicit none
   integer                     , intent(in)    :: n1,n2,n3
   real   , dimension(n1,n2,n3), intent(inout) :: a,b,c
   real   , dimension(n1,n2)   , intent(in)    :: topt
   integer                     , intent(in)    :: ngrd
   integer                                     :: i,j,k
   real   , dimension(n3)                      :: scratch2,scratch11,scratch12
   real                                        :: zedtop,c1,c2,c3

   zedtop = myzmn(mynnzp(1)-1,1)
   do j=1,n2
      do i=1,n1
         do k=1,n3
           scratch2(k)=myztn(k,ngrd)*(1.-topt(i,j)/zedtop)+topt(i,j)
         enddo
         call htint(n3,mypi01dn(1,ngrd),myztn(1,ngrd),n3,scratch11,scratch2)
         call htint(n3,myth01dn(1,ngrd),myztn(1,ngrd),n3,scratch12,scratch2)       
         do k=1,n3
            b(i,j,k)=scratch12(k)
         enddo
         a(i,j,n3) = scratch11(n3)

         c1=grav*2.*(1.-topt(i,j)/zedtop)
         c2=(1-cpor)
         c3=cp**c2
         do k=n3-1,1,-1
            a(i,j,k)=a(i,j,k+1) +c1/((b(i,j,k)+b(i,j,k+1))*mydzmn(k,ngrd))
         enddo
         do k=1,n3
            c(i,j,k)=(c3*p00)/(rdry*b(i,j,k)*a(i,j,k)**c2)
         enddo

      enddo
   enddo
   return
end subroutine RAMS_comp_dn0
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine RAMS_comp_relvortx(n1,n2,n3,a,b,c,topt,ngrd)

   use somevars
   use rconstants
   use therm_lib, only : qtk, qwtk, dewfrostpoint, rslif, virtt, thetaeiv
   implicit none
   integer                     , intent(in)    :: n1,n2,n3
   real   , dimension(n1,n2,n3), intent(inout) :: a,b,c
   real   , dimension(n1,n2)   , intent(in)    :: topt
   integer                     , intent(in)    :: ngrd
   integer                                     :: i,j,k,j1,j2,k1,k2
   real                                        :: factor
   real   , dimension(n1+n2+n3)                :: dum1,dum2

   factor = myztn(1,ngrd) / myztn(2,ngrd)
   do j=1,n2
      do i=1,n1
         a(i,j,1) = a(i,j,2) * factor
      enddo
   enddo

   call gradr(n1,n2,n3,2,n1-1,1,n2-1,b,c,'ydir','wpnt',topt  &
      ,myxmn(:,ngrd),myxtn(:,ngrd),myymn(:,ngrd),myytn(:,ngrd)  &
      ,myzmn(:,ngrd),myztn(:,ngrd),mydeltayn(ngrd)  &
      ,mydzmn(:,ngrd),mydztn(:,ngrd),dum1,dum2,myzmn(mynnzp(1)-1,1)  &
      ,myjdim,myihtran)
   call gradr(n1,n2,n3,2,n1-1,1,n2-1,a,b,'zdir','vpnt',topt  &
      ,myxmn(:,ngrd),myxtn(:,ngrd),myymn(:,ngrd),myytn(:,ngrd)  &
      ,myzmn(:,ngrd),myztn(:,ngrd),mydeltayn(ngrd)  &
      ,mydzmn(:,ngrd),mydztn(:,ngrd),dum1,dum2,myzmn(mynnzp(1)-1,1)  &
      ,myjdim,myihtran)

   do k=1,n3
      do j=2,n2-1
         do i=1,n1
            b(i,j,k) = c(i,j,k) - b(i,j,k)
         enddo
      enddo
   enddo
   do j = 1,n2
      do i =1,n1
         do k =1,n3
            a(i,j,k) = 0.
         end do
      end do
   end do
   do j = 2,n2-1
      j1 = max(j-1,2)
      j2 = min(j,n2-2)
      do i = 2,n1-1
         do k = 1,n3
            k1 = max(k-1,1)
            k2 = min(k,n3-1)
            a(i,j,k) =0.25 * (b(i,j1,k1) + b(i,j1,k2)  &
                            + b(i,j2,k1) + b(i,j2,k2))
         enddo
      enddo
   enddo
   
  return
end subroutine RAMS_comp_relvortx
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine RAMS_comp_relvorty(n1,n2,n3,a,b,c,topt,ngrd)

   use somevars
   use rconstants
   use therm_lib, only : qtk, qwtk, dewfrostpoint, rslif, virtt, thetaeiv
   implicit none
   integer                     , intent(in)    :: n1,n2,n3
   real   , dimension(n1,n2,n3), intent(inout) :: a,b,c
   real   , dimension(n1,n2)   , intent(in)    :: topt
   integer                     , intent(in)    :: ngrd
   integer                                     :: i,j,k,i1,i2,k1,k2
   real                                        :: factor
   real   , dimension(n1+n2+n3)                :: dum1,dum2

   factor = myztn(2,ngrd) / myztn(3,ngrd)
   do j=1,n2
      do i=1,n1
         a(i,j,1) = a(i,j,2) * factor
      enddo
   enddo

   call gradr(n1,n2,n3,1,n1-1,2,n2-1,b,c,'xdir','wpnt',topt  &
      ,myxmn(:,ngrd),myxtn(:,ngrd),myymn(:,ngrd),myytn(:,ngrd)  &
      ,myzmn(:,ngrd),myztn(:,ngrd),mydeltayn(ngrd)  &
      ,mydzmn(:,ngrd),mydztn(:,ngrd),dum1,dum2,myzmn(mynnzp(1)-1,1)  &
      ,myjdim,myihtran)
   call gradr(n1,n2,n3,1,n1-1,2,n2-1,a,b,'zdir','upnt',topt  &
      ,myxmn(:,ngrd),myxtn(:,ngrd),myymn(:,ngrd),myytn(:,ngrd)  &
      ,myzmn(:,ngrd),myztn(:,ngrd),mydeltayn(ngrd)  &
      ,mydzmn(:,ngrd),mydztn(:,ngrd),dum1,dum2,myzmn(mynnzp(1)-1,1)  &
      ,myjdim,myihtran)

   do k=1,n3
      do j=1,n2
         do i=1,n1
            b(i,j,k) = b(i,j,k) - c(i,j,k)
         enddo
      enddo
   enddo

   do j = 1,n2
      do i =1,n1
         do k =1,n3
            a(i,j,k) = 0.
         end do
      end do
   end do

   do j = 2,n2-1
      do i = 2,n1-1
         i1 = max(i-1,2)
         i2 = min(i,n1-2)
         do k = 1,n3
            k1 = max(k-1,1)
            k2 = min(k,n3-1)
            a(i,j,k) = 0.25 * (b(i1,j,k1) + b(i1,j,k2)  &
                             + b(i2,j,k1) + b(i2,j,k2))
         enddo
      enddo
   enddo


   return
end subroutine RAMS_comp_relvorty
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine RAMS_comp_relvortz(n1,n2,n3,a,b,c,topt,ngrd)

   use somevars
   use rconstants
   use therm_lib, only : qtk, qwtk, dewfrostpoint, rslif, virtt, thetaeiv
   implicit none
   integer                     , intent(in)    :: n1,n2,n3
   real   , dimension(n1,n2,n3), intent(inout) :: a,b,c
   real   , dimension(n1,n2)   , intent(in)    :: topt
   integer                     , intent(in)    :: ngrd
   integer                                     :: i,j,k,i1,i2,j1,j2
   real                                        :: factor
   real   , dimension(n1+n2+n3)                :: dum1,dum2

   factor = myztn(1,ngrd) / myztn(2,ngrd)
   do j=1,n2
      do i=1,n1
         a(i,j,1) = a(i,j,2) * factor
         b(i,j,1) = b(i,j,2) * factor
      enddo
   enddo

   call gradr(n1,n2,n3,1,n1-1,1,n2-1,b,c,'xdir','vpnt',topt  &
      ,myxmn(:,ngrd),myxtn(:,ngrd),myymn(:,ngrd),myytn(:,ngrd)  &
      ,myzmn(:,ngrd),myztn(:,ngrd),mydeltayn(ngrd)  &
      ,mydzmn(:,ngrd),mydztn(:,ngrd),dum1,dum2,myzmn(mynnzp(1)-1,1)  &
      ,myjdim,myihtran)

   call gradr(n1,n2,n3,1,n1-1,1,n2-1,a,b,'ydir','upnt',topt  &
      ,myxmn(:,ngrd),myxtn(:,ngrd),myymn(:,ngrd),myytn(:,ngrd)  &
      ,myzmn(:,ngrd),myztn(:,ngrd),mydeltayn(ngrd)  &
      ,mydzmn(:,ngrd),mydztn(:,ngrd),dum1,dum2,myzmn(mynnzp(1)-1,1)  &
      ,myjdim,myihtran)

   do k=1,n3
      do j=1,n2
         do i=1,n1
            b(i,j,k) = c(i,j,k) - b(i,j,k)
         enddo
      enddo
   enddo

   do j = 1,n2
      j1 = max(j-1,1)
      j2 = min(j,n2-1)
      do i = 1,n1
         i1 = max(i-1,1)
         i2 = min(i,n1-1)
         do k = 1,n3
            a(i,j,k) = 0.25 * (b(i1,j1,k) + b(i1,j2,k)  &
                             + b(i2,j1,k) + b(i2,j2,k))
         enddo
      enddo
   enddo
   return
end subroutine RAMS_comp_relvortz
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine RAMS_comp_totvortz(n1,n2,n3,a,b,c,topt,ngrd)

   use somevars
   use rconstants
   use therm_lib, only : qtk, qwtk, dewfrostpoint, rslif, virtt, thetaeiv
   implicit none
   integer                     , intent(in)    :: n1,n2,n3
   real   , dimension(n1,n2,n3), intent(inout) :: a,b,c
   real   , dimension(n1,n2)   , intent(in)    :: topt
   integer                     , intent(in)    :: ngrd
   integer                                     :: i,j,k,i1,i2,j1,j2
   real                                        :: factor,omega2,fcor,xlon,xlat
   real   , dimension(n1+n2+n3)                :: dum1,dum2

   factor = myztn(1,ngrd) / myztn(2,ngrd)
   do j=1,n2
      do i=1,n1
         a(i,j,1) = a(i,j,2) * factor
         b(i,j,1) = b(i,j,2) * factor
      enddo
   enddo

   call gradr(n1,n2,n3,1,n1-1,1,n2-1,b,c,'xdir','vpnt',topt  &
      ,myxmn(:,ngrd),myxtn(:,ngrd),myymn(:,ngrd),myytn(:,ngrd)  &
      ,myzmn(:,ngrd),myztn(:,ngrd),mydeltayn(ngrd)  &
      ,mydzmn(:,ngrd),mydztn(:,ngrd),dum1,dum2,myzmn(mynnzp(1)-1,1)  &
      ,myjdim,myihtran)
   call gradr(n1,n2,n3,1,n1-1,1,n2-1,a,b,'ydir','upnt',topt  &
      ,myxmn(:,ngrd),myxtn(:,ngrd),myymn(:,ngrd),myytn(:,ngrd)  &
      ,myzmn(:,ngrd),myztn(:,ngrd),mydeltayn(ngrd)  &
      ,mydzmn(:,ngrd),mydztn(:,ngrd),dum1,dum2,myzmn(mynnzp(1)-1,1)  &
      ,myjdim,myihtran)

   do k=1,n3
      do j=1,n2
         do i=1,n1
            b(i,j,k) = c(i,j,k) - b(i,j,k)
         enddo
      enddo
   enddo

   do j = 1,n2
      j1 = max(j-1,1)
      j2 = min(j,n2-1)
      do i = 1,n1
         i1 = max(i-1,1)
         i2 = min(i,n1-1)
         do k = 1,n3
            a(i,j,k) = 0.25 * (b(i1,j1,k) + b(i1,j2,k)  &
                             + b(i2,j1,k) + b(i2,j2,k))
         enddo
      enddo
   enddo

   omega2 = 2. * omega
   do j = 1,n2
      do i = 1,n1
         call xy_ll(xlat,xlon,myplatn(ngrd),myplonn(ngrd),myxtn(i,ngrd),myytn(j,ngrd))
         fcor = omega2 * sin(xlat * pio180)
         do k = 1,n3
            a(i,j,k) = a(i,j,k) + fcor
         enddo
      enddo
   enddo
   return
end subroutine RAMS_comp_totvortz
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine RAMS_comp_potvortz(n1,n2,n3,a,b,c,e,topt,ngrd)
   use somevars
   use rconstants
   use therm_lib, only : qtk, qwtk, dewfrostpoint, rslif, virtt, thetaeiv
   implicit none
   integer                     , intent(in)    :: n1,n2,n3
   real   , dimension(n1,n2,n3), intent(inout) :: a,b,c,e
   real   , dimension(n1,n2)   , intent(in)    :: topt
   integer                     , intent(in)    :: ngrd
   integer                                     :: i,j,k
   real   , dimension(n1+n2+n3)                :: dum1,dum2

   call gradr(n1,n2,n3,1,n1-1,1,n2-1,b,e,'zdir','tpnt',topt  &
      ,myxmn(:,ngrd),myxtn(:,ngrd),myymn(:,ngrd),myytn(:,ngrd)  &
      ,myzmn(:,ngrd),myztn(:,ngrd),mydeltayn(ngrd)  &
      ,mydzmn(:,ngrd),mydztn(:,ngrd),dum1,dum2,myzmn(mynnzp(1)-1,1)  &
      ,myjdim,myihtran)

   do k=1,n3
      do j=1,n2
         do i=1,n1
            a(i,j,k) = a(i,j,k) * e(i,j,k) / (grav * c(i,j,k))
         enddo
      enddo
   enddo
   return
end subroutine RAMS_comp_potvortz
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine RAMS_comp_solenoidx(n1,n2,n3,a,b,e,topt,ngrd)
   use somevars
   use rconstants
   use therm_lib, only : qtk, qwtk, dewfrostpoint, rslif, virtt, thetaeiv
   implicit none
   integer                     , intent(in)    :: n1,n2,n3
   real   , dimension(n1,n2,n3), intent(inout) :: a,b,e
   real   , dimension(n1,n2)   , intent(in)    :: topt
   integer                     , intent(in)    :: ngrd
   real   , dimension(n1,n2,n3)                :: d,o
   integer                                     :: i,j,k
   real   , dimension(n1+n2+n3)                :: dum1,dum2
   do k=1,n3
      do j=1,n2
         do i=1,n1
            e(i,j,k) =0.
         enddo
      enddo
   enddo
   ! a - specific volume
   ! b - pressure
   ! o - derivative
   ! d - derivative
   ! e - sum - solenoid term
   call gradr(n1,n2,n3,2,n1-1,2,n2-1,a,o,'zdir','tpnt',topt  &
      ,myxmn(:,ngrd),myxtn(:,ngrd),myymn(:,ngrd),myytn(:,ngrd)  &
      ,myzmn(:,ngrd),myztn(:,ngrd),mydeltayn(ngrd)  &
      ,mydzmn(:,ngrd),mydztn(:,ngrd),dum1,dum2,myzmn(mynnzp(1)-1,1)  &
      ,myjdim,myihtran)
   ! o is da/dz
   call gradr(n1,n2,n3,2,n1-1,2,n2-1,b,d,'ydir','tpnt',topt  &
      ,myxmn(:,ngrd),myxtn(:,ngrd),myymn(:,ngrd),myytn(:,ngrd)  &
      ,myzmn(:,ngrd),myztn(:,ngrd),mydeltayn(ngrd)  &
      ,mydzmn(:,ngrd),mydztn(:,ngrd),dum1,dum2,myzmn(mynnzp(1)-1,1)  &
      ,myjdim,myihtran)
   !d is dp/dy

   do k=1,n3
      do j=2,n2-1
         do i=2,n1-1
            e(i,j,k) =o(i,j,k)*d(i,j,k)
         enddo
      enddo
   enddo


   call gradr(n1,n2,n3,2,n1-1,2,n2-1,b,o,'zdir','tpnt',topt  &
      ,myxmn(:,ngrd),myxtn(:,ngrd),myymn(:,ngrd),myytn(:,ngrd)  &
      ,myzmn(:,ngrd),myztn(:,ngrd),mydeltayn(ngrd)  &
      ,mydzmn(:,ngrd),mydztn(:,ngrd),dum1,dum2,myzmn(mynnzp(1)-1,1)  &
      ,myjdim,myihtran)
   ! d is dp/dz

   call gradr(n1,n2,n3,2,n1-1,2,n2-1,a,d,'ydir','tpnt',topt  &
      ,myxmn(:,ngrd),myxtn(:,ngrd),myymn(:,ngrd),myytn(:,ngrd)  &
      ,myzmn(:,ngrd),myztn(:,ngrd),mydeltayn(ngrd)  &
      ,mydzmn(:,ngrd),mydztn(:,ngrd),dum1,dum2,myzmn(mynnzp(1)-1,1)  &
      ,myjdim,myihtran)
   !c is da/dy


   do k=1,n3
      do j=2,n2-1
         do i=2,n1-1
            e(i,j,k) = e(i,j,k)-o(i,j,k)*d(i,j,k)
         enddo
      enddo
   enddo
  return
end subroutine RAMS_comp_solenoidx
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine RAMS_comp_solenoidy(n1,n2,n3,a,b,e,topt,ngrd)
   use somevars
   use rconstants
   use therm_lib, only : qtk, qwtk, dewfrostpoint, rslif, virtt, thetaeiv
   implicit none
   integer                     , intent(in)    :: n1,n2,n3
   real   , dimension(n1,n2,n3), intent(inout) :: a,b,e
   real   , dimension(n1,n2)   , intent(in)    :: topt
   integer                     , intent(in)    :: ngrd
   real   , dimension(n1,n2,n3)                :: d,o
   integer                                     :: i,j,k
   real   , dimension(n1+n2+n3)                :: dum1,dum2
   do k=1,n3
      do j=1,n2
         do i=1,n1
            e(i,j,k) =0.
         enddo
      enddo
   enddo
   ! a - specific volume
   ! b - pressure
   ! o - derivative
   ! d - derivative
   ! e - sum - solenoid term
   call gradr(n1,n2,n3,2,n1-1,2,n2-1,a,o,'xdir','tpnt',topt  &
      ,myxmn(:,ngrd),myxtn(:,ngrd),myymn(:,ngrd),myytn(:,ngrd)  &
      ,myzmn(:,ngrd),myztn(:,ngrd),mydeltayn(ngrd)  &
      ,mydzmn(:,ngrd),mydztn(:,ngrd),dum1,dum2,myzmn(mynnzp(1)-1,1)  &
      ,myjdim,myihtran)
   !c is da/dx
   call gradr(n1,n2,n3,2,n1-1,2,n2-1,b,d,'zdir','tpnt',topt  &
      ,myxmn(:,ngrd),myxtn(:,ngrd),myymn(:,ngrd),myytn(:,ngrd)  &
      ,myzmn(:,ngrd),myztn(:,ngrd),mydeltayn(ngrd)  &
      ,mydzmn(:,ngrd),mydztn(:,ngrd),dum1,dum2,myzmn(mynnzp(1)-1,1)  &
      ,myjdim,myihtran)
   !p is dp/dz

   do k=1,n3
      do j=2,n2-1
         do i=2,n1-1
            e(i,j,k) =o(i,j,k)*d(i,j,k)
         enddo
      enddo
   enddo

   call gradr(n1,n2,n3,2,n1-1,2,n2-1,b,o,'xdir','tpnt',topt  &
      ,myxmn(:,ngrd),myxtn(:,ngrd),myymn(:,ngrd),myytn(:,ngrd)  &
      ,myzmn(:,ngrd),myztn(:,ngrd),mydeltayn(ngrd)  &
      ,mydzmn(:,ngrd),mydztn(:,ngrd),dum1,dum2,myzmn(mynnzp(1)-1,1)  &
      ,myjdim,myihtran)
   !c is dp/dx

   call gradr(n1,n2,n3,2,n1-1,2,n2-1,a,d,'zdir','tpnt',topt  &
      ,myxmn(:,ngrd),myxtn(:,ngrd),myymn(:,ngrd),myytn(:,ngrd)  &
      ,myzmn(:,ngrd),myztn(:,ngrd),mydeltayn(ngrd)  &
      ,mydzmn(:,ngrd),mydztn(:,ngrd),dum1,dum2,myzmn(mynnzp(1)-1,1)  &
      ,myjdim,myihtran)
   ! d is da/dz


   do k=1,n3
      do j=2,n2-1
         do i=2,n1-1
            e(i,j,k) = e(i,j,k)-o(i,j,k)*d(i,j,k)
         enddo
      enddo
   enddo
  return
end subroutine RAMS_comp_solenoidy
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine RAMS_comp_sfcwmeantemp(n1,n2,ns,np,a,b,c,d,e)
   use rconstants
   use therm_lib, only: qtk
   implicit none 
   integer :: n1,n2,ns,np,nlev,i,j,ip,k
   real, dimension(n1,n2,np)    :: a,d,e
   real, dimension(n1,n2,ns,np) :: b,c
   real :: temptemp,fracliq,snowarea,xmasstot
   real, parameter :: undef=-9.99e33
   !a  area
   !b energy
   !c  mass
   !d  nlev
   !e integrated value
   do j=1,n2
     do i=1,n1
        e(i,j,1) = 0.
        if (a(i,j,1) <= 0.99) then
           do ip=2,np
              xmasstot = 0.
              nlev=nint(d(i,j,ip))
              e(i,j,ip) = 0.
              if (nlev > 0 ) then
                 do k=1,nlev
                    if (c(i,j,k,ip) > 1.e-6) then
                       xmasstot  = xmasstot + c(i,j,k,ip)
                       call qtk(b(i,j,k,ip),temptemp,fracliq)
                       e(i,j,ip) = e(i,j,ip) + temptemp*c(i,j,k,ip)
                    end if
                 end do
                 if (xmasstot > 1.e-6) e(i,j,ip) = e(i,j,ip) / xmasstot - t00
              else
                 e(i,j,ip) = 0.
              end if
           end do

           snowarea = 0.
           do ip=2,np
              if( nint(d(i,j,ip)) > 0) then
                 snowarea = snowarea + a(i,j,ip)
                 e(i,j,1) = e(i,j,1) + e(i,j,ip)*a(i,j,ip)
              end if
           end do
           if (snowarea > 0.) then
              e(i,j,1) = e(i,j,1) / snowarea
           else
              e(i,j,1) = undef
           end if
        else
           e(i,j,1) = undef
        end if
     end do
   end do
   return
end subroutine RAMS_comp_sfcwmeantemp
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine RAMS_comp_sfcwinteg(n1,n2,ns,np,a,c,d,e)
   implicit none 
   integer :: n1,n2,ns,np,nlev,i,j,ip,k
   real, dimension(n1,n2,np)    :: a,d,e
   real, dimension(n1,n2,ns,np) :: c
   real, parameter :: undef=-9.99e33
   !a  area                                                   
   !c mass/depth                                             
   !d  nlev                                                   
   !e  integrated value                                       
      do j=1,n2                                               
        do i=1,n1                                             
           if (a(i,j,1) > 0.99) then                          
              e(i,j,1) = -9.99e33                                
           else
              e(i,j,1) = 0.                                  
              do ip=2,np                                     
                 nlev=nint(d(i,j,ip))
                 e(i,j,ip) = 0.                               
                 do k=1,nlev                                  
                    e(i,j,ip) = e(i,j,ip) + c(i,j,k,ip)      
                 end do                                       
                 e(i,j,1) = e(i,j,1) + e(i,j,ip)*a(i,j,ip)    
              end do                                          
              e(i,j,1) = e(i,j,1) / (1.-a(i,j,1))            
           end if                                             
        end do                                                
      end do                                                  
   return        
end subroutine RAMS_comp_sfcwinteg
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This routine is for quantities defined for all patches.                               !
!------------------------------------------------------------------------------------------!
subroutine RAMS_comp_patchsum(n1,n2,n3,n4,a,f,psum)
   use somevars, only : mynbig
   implicit none
   integer                         , intent(in)    :: n1,n2,n3,n4
   real   , dimension(n1,n2,n3,n4) , intent(in)    :: a
   real   , dimension(n1,n2,mynbig), intent(inout) :: f
   real   , dimension(n1,n2,n3)    , intent(out)   :: psum
   integer                                         :: i,j,k,ip

   do k = 1,n3
      do j = 1,n2
         do i = 1,n1
            psum(i,j,k) = 0.
            do ip = 1,n4
               psum(i,j,k) = psum(i,j,k) + f(i,j,ip) * a(i,j,k,ip)
            end do
         end do
      end do
   end do

   !---------------------------------------------------------------------------------------!
   !      Copy psum into f, which was passed in as a(1).                                   !
   !---------------------------------------------------------------------------------------!
   do k = 1,n3
      do j = 1,n2
         do i = 1,n1
            f(i,j,k) = psum(i,j,k)
         end do
      end do
   end do
   return
end subroutine RAMS_comp_patchsum
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!      This routine is for quantities that are not defined for water  patches.             !
!------------------------------------------------------------------------------------------!
subroutine RAMS_comp_patchsum_l(n1,n2,n3,n4,a,f,psum)
   use somevars, only : mynbig
   implicit none
   integer, intent(in) :: n1,n2,n3,n4
   real   , dimension(n1,n2,n3,n4) , intent(in)    :: a
   real   , dimension(n1,n2,mynbig), intent(inout) :: f
   real   , dimension(n1,n2,n3)    , intent(out)   :: psum
   integer                                         :: i,j,k,ip


   do k = 1,n3
      do j = 1,n2
         do i = 1,n1
            if (f(i,j,1) < .991) then
               psum(i,j,k) = 0.
               do ip = 2,n4
                  psum(i,j,k) = psum(i,j,k) + f(i,j,ip) * a(i,j,k,ip) / (1. - f(i,j,1))
               end do
            else
               psum(i,j,k) = a(i,j,k,2)
            end if
         end do
      end do
   end do

   !---------------------------------------------------------------------------------------!
   !      Copy psum into f, which was passed in as a(1).                                   !
   !---------------------------------------------------------------------------------------!
   do k = 1,n3
      do j = 1,n2
         do i = 1,n1
            f(i,j,k) = psum(i,j,k)
         end do
      end do
   end do
   return
end subroutine RAMS_comp_patchsum_l
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!      Extract value from largest patch.                                                   !
!------------------------------------------------------------------------------------------!
subroutine RAMS_comp_bigpatch(n1,n2,n3,n4,a,f,bpat)
   use somevars, only : mynbig
   implicit none
   integer, intent(in) :: n1,n2,n3,n4
   real   , dimension(n1,n2,n3,n4) , intent(in)    :: a
   real   , dimension(n1,n2,mynbig), intent(inout) :: f
   real   , dimension(n1,n2,n3)    , intent(out)   :: bpat
   integer                                         :: i,j,k,ip
   do k = 1,n3
      do j = 1,n2
         do i = 1,n1
            if (f(i,j,2) >= f(i,j,1)) then
               bpat(i,j,k) = a(i,j,k,2)
            else
               bpat(i,j,k) = a(i,j,k,1)
            end if
         end do
      end do
   end do

   !---------------------------------------------------------------------------------------!
   !      Copy bpat into f, which was passed in as a(1).                                   !
   !---------------------------------------------------------------------------------------!

   do k = 1,n3
      do j = 1,n2
         do i = 1,n1
            f(i,j,k) = bpat(i,j,k)
         end do
      end do
   end do
   return
end subroutine RAMS_comp_bigpatch
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine RAMS_comp_tvegc(n1,n2,n3,a,b,c,e)
   use rconstants, only : t00
   use therm_lib , only : qwtk
   implicit none
   integer, intent(in) :: n1,n2,n3
   real, dimension(n1,n2,n3), intent(in)    :: c,e
   real, dimension(n1,n2,n3), intent(inout) :: a,b
   integer :: i,j,k
   real :: temptemp,fracliq
!MLO: Input:  a = Veg. Energy        in     J/m2
!             b = Veg. Water         in    kg/m2
!             c = Veg. Heat capacity in   J/m2/K
!             e = Canopy Theta       in  Celsius
!     Output: a = Temperature in Celsius
!             b = Fraction in liquid phase

   do k=1,n3
      do j=1,n2
         do i=1,n1
            !------------------------------------------------------------------------ 
            !   Compute tveg only if there is enough heat capacity, otherwise assign
            ! canopy temperature.
            !------------------------------------------------------------------------ 
            if (c(i,j,k) > 10.) then
               call qwtk(a(i,j,k),b(i,j,k),c(i,j,k),temptemp,fracliq)
               a(i,j,k) = temptemp-t00
               b(i,j,k) = fracliq
            else
               a(i,j,k) = e(i,j,k)
               b(i,j,k) = 0.5
            end if
         end do
      end do
   end do
   return
end subroutine RAMS_comp_tvegc

subroutine RAMS_comp_5050(n1,n2,n3,a,d)
real a(n1,n2),d(n1,n2,n3)

do j = 1,n2
   do i = 1,n1
      a(i,j) = .5 * (a(i,j) + d(i,j,2))
   enddo
enddo

return
end subroutine RAMS_comp_5050
!------------------subroutine to calculate cloud fraction
subroutine cldfraction(n1,n2,n3,frac,pi,rh)
use rconstants
implicit none

integer :: i,j,k,kmax,n1,n2,n3
real :: frac(n1,n2),pi(n1,n2,n3),rh(n1,n2,n3)
real, allocatable::rhc(:), cs(:)
real :: kappai,c_1,c_2,c_junk,pop2,csmax

c_1     = 2.
c_junk  = 3.
c_2     = c_junk**0.5
kappai = (1./.286)


allocate (rhc(n3),cs(n3) )
      print*,'+++++++:',n1,n2,n3

do j=1,n2
   do i=1,n1
      frac(i,j) = 0.
      csmax = 0.
      kmax  = 0
      do k = 1, n3
         rhc(k)= 0.
         cs(k)= 0.
      enddo

      do k = 1, n3
         pop2 = (pi(i,j,k)/pi(i,j,2))**kappai

         rhc(k) = 100. - (100.*c_1*pop2)*  &
               (1.-pop2)*(1.+c_2*(pop2-0.5))

         if(rh(i,j,k) .ge. rhc(k))then
            if(rhc(k).eq.100.)rhc(k)=rhc(k)+0.0000001
            cs(k) = ( (rh(i,j,k)-rhc(k))/(100.-rhc(k)) ) **2. 
         else
            cs(k) = 0.
         endif
         if(cs(k).gt.csmax)then
            csmax=cs(k)
            kmax = k
         endif
         frac(i,j) = frac(i,j) + cs(k)*(1./float(k))
      if(i==20.and.j==20) print*,'+++++++:',k,pi(i,j,k),rh(i,j,k),frac(i,j)
      enddo
      
      csmax=max(csmax,0.)

!      frac(i,j) = 1.-min(1.,max(0.,csmax))
      frac(i,j) = 1.-min(1.,max(0.,frac(i,j)))  ! actually returns 
                                                ! clear sky fraction

   enddo
enddo

deallocate (rhc,cs)

return
end
!***************************************************************************

subroutine RAMS_reduced_temp (n1,n2,n3,n4,tempnew,speed,ustar  &
                             ,tstar,znew,zold,zrough,patfrac  &
                             ,cantemp,theta,pi,topo,ztop)
        
use rconstants
implicit none

integer :: n1,n2,n3,n4,i,j,np
real :: tempnew(n1,n2),speed(n1,n2,n3),ustar(n1,n2,n4),znew,zold  &
       ,zrough(n1,n2,n4),patfrac(n1,n2,n4),cantemp(n1,n2,n4)  &
       ,theta(n1,n2,n3),pi(n1,n2,n3),topo(n1,n2),ztop,tstar(n1,n2,n4)

real :: richno,rtgt,zagl,rtemp,rtempw,z0,a2,spd,cantheta,sfcpi,fh
!srf: consistent with Louis 1981
fh = 1.

do j=1,n2
   do i=1,n1
      
      rtgt=1.-topo(i,j)/ztop
      zagl=zold*rtgt
      sfcpi=.5*(pi(i,j,1)+pi(i,j,2))

      
      rtempw=0.
      
      do np=1,n4
      
         z0=max(zrough(i,j,np),0.001)
         if(np==1) z0=.001
         spd=max(speed(i,j,2),.25)
         cantheta=cantemp(i,j,np)*cp/sfcpi

!-srf
!         richno=g*zagl*(theta(i,j,2)-cantheta)  &
!                     /(theta(i,j,2)*spd**2)
         richno=grav*zagl*(theta(i,j,2)-cantheta)  &
                     /(.5*(theta(i,j,2)+cantheta)*spd**2)
         a2 = (vonk / log(znew / z0)) ** 2

        !print*,richno,spd,theta(i,j,2),zagl


         if(richno.gt.0.) then
            rtemp=cantheta                            &
!srf         +(ustar(i,j,np)*tstar(i,j,np)*0.74)/(a2*spd)  &
             +(ustar(i,j,np)*tstar(i,j,np)*fh  )/(a2*spd)  &
                    *(1.+15.*richno*sqrt(1+5*richno))  
            
	    
            !print*,rtemp, cantheta,theta(i,j,2)
	    rtemp=min(max(rtemp, cantheta),theta(i,j,2))

        else
            rtemp=cantheta                              &
!srf         +((ustar(i,j,np)*tstar(i,j,np)*0.74)/(a2*spd))  &
             +((ustar(i,j,np)*tstar(i,j,np)*fh  )/(a2*spd))  &
                              / (1.- 15.*richno/(1.+75.*a2   &
                              * sqrt(-znew*richno/z0)))

            
            !print*,rtemp, cantheta,theta(i,j,2)
	    rtemp=max(min(rtemp, cantheta),theta(i,j,2))


         endif
         
         !if((i==50.and.j==25)) then
         !   print*,'====tempf2m:',i,j
         !   print*,np,patfrac(i,j,np),cantheta
         !   print*,np,ustar(i,j,np),zrough(i,j,np),tstar(i,j,np)
         !  print*,np,theta(i,j,2),speed(i,j,2),rtemp
         !endif
         
         rtempw=rtempw+rtemp*patfrac(i,j,np)

  
      enddo
      
     tempnew(i,j)=rtempw ! temperatura potencial

      
   enddo
enddo

return
end
!-----------------------------------------------------------------------

subroutine RAMS_reduced_wind(n1,n2,n3,n4,velnew,speed,ustar &
         ,znew,zold,zrough,patfrac,cantemp,theta,pi,topo,ztop)
use rconstants
implicit none
integer :: n1,n2,n3,n4,i,j,np
real :: velnew(n1,n2),speed(n1,n2,n3),ustar(n1,n2,n4),znew,zold  &
          ,zrough(n1,n2,n4),patfrac(n1,n2,n4),cantemp(n1,n2,n4)  &
          ,theta(n1,n2,n3),pi(n1,n2,n3),topo(n1,n2),ztop

real:: richno,rtgt,zagl,rwind,rwindw,z0,a2,spd,cantheta,sfcpi



do j=1,n2
   do i=1,n1
      
      rtgt=1.-topo(i,j)/ztop
      zagl=zold*rtgt
      sfcpi=.5*(pi(i,j,1)+pi(i,j,2))
      
      rwindw=0.
      
      do np=1,n4
      
         z0=zrough(i,j,np)
         if(np==1) z0=.001
         spd=max(speed(i,j,2),.25)
         cantheta=cantemp(i,j,np)*cp/sfcpi

         richno=grav*zagl*(theta(i,j,2)-cantheta)  &
                      /(theta(i,j,2)*spd**2)
         a2 = (vonk / log(znew / z0)) ** 2

         if(richno.gt.0.) then
            rwind=sqrt(ustar(i,j,np)**2/a2   &
                     *(1.+10.*richno/sqrt(1+5*richno)) )
         else
            rwind=sqrt( ustar(i,j,np)**2/a2  &
                / (1.- 10.*richno/(1.+75.*a2  &
                              * sqrt(-znew*richno/z0))))
         endif
         
         rwind=max(min(rwind,speed(i,j,2)),0.)
         
         !if(i==50.and.j==25) then
         !   print*,'====speed10m'
         !   print*,np,patfrac(i,j,np),cantemp(i,j,np)
         !   print*,np,ustar(i,j,np),zrough(i,j,np)
         !   print*,np,theta(i,j,2),speed(i,j,2),rwind
         !endif
         
         rwindw=rwindw+rwind*patfrac(i,j,np)
      
      enddo
      
      velnew(i,j)=rwindw
      

   enddo
enddo

return
end

!***************************************************************************
!-------------------------------------------------------------------------
!*rmc Will Cheng's code for calculating slp with mm5's GRAPH method
! ------- added for calculating SLP from MM5 algorithm ------

      subroutine RAMS_comp_slpmm5(n1,n2,n3,theta,pp,z,slp)

!    The subroutine calculates SLP from an algorithm taken from
!    GRAPH, a post-processing packing of MM5 V3.3
!
!    Input: theta - potential temperature (K)         3D
!           pp    - Exner function        (J/kg K)    3D
!           z     - terrain               (m)         2D
!
!    Ouput: SLP   - sea-level pressure    (hPa)       2D

! ------ define dimension of arrays ----------

      dimension theta(n1,n2,n3), pp(n1,n2,n3), z(n1,n2),slp(n1,n2)

! ------ input variables to GRAPH subroutine ----------

      dimension sfp(n1,n2), ts(n1,n2), t_mm5(n1,n2,n3-1), p_mm5(n1,n2,n3-1)

! -----------------------------------------------------

      cp = 1004
      rgas = 287
      cpor = cp / rgas
      p00 = 1.e5

      do j = 1,n2
        do i = 1,n1
!! calculate surface pressure
        sfp(i,j) = (0.5*(pp(i,j,1)+pp(i,j,2))/cp)**cpor*p00*.01
!! calculate surface temp
        ts(i,j) = (0.5/cp)*(theta(i,j,1)*pp(i,j,1)+&
                            theta(i,j,2)*pp(j,j,2))
        enddo
      enddo

      do k = 2,n3
        kk = n3-k+1
        do j = 1,n2
!! flip array upside down for input to GRAPH subroutine
          do i = 1,n1
            t_mm5(i,j,kk) = theta(i,j,k)*pp(i,j,k)/cp
            p_mm5(i,j,kk) = (pp(i,j,k)/cp)**cpor*p00*.01
          enddo
        enddo
      enddo

      call SEAPRS_0(t_mm5,p_mm5,z,sfp,ts,n1,n2,n3-1,slp)

      return
      end


!------------------------------------------------------------------------
      SUBROUTINE SEAPRS_0(T,PP,TER,SFP,TS,IMX,JMX,KX,SLP)
!
!     SECTION  DIAGNOSTIC
!     PURPOSE  COMPUTES SEA LEVEL PRESSURE FROM THE RULE
!              T1/T2=(P1/P2)**(GAMMA*R/G).
!
!     *** LEVELS GO FROM TOP-DOWN ***
!
!     INPUT       T        TEMPERATURE (Kelvin)                3D
!                 TER      TERRAIN     (m)                     2D
!                 SFP      SURFACE PRESSURE (hPa)              2D
!                 IMX      DOT POINT DIMENSION N-S
!                 JMX      DOT POINT DIMENSION E-W
!                 KX       NUMBER OF VERTICAL LEVELS
!
!     OUTPUT      SLP      SEA LEVEL PRESSURE (hPa)            2D
!
      DIMENSION T(IMX,JMX,KX), PP(IMX,JMX,KX),&
                PS(IMX,JMX)  ,SFP(IMX,JMX) , &
                TER(IMX,JMX)
      DIMENSION PL(IMX,JMX),T0(IMX,JMX),TS(IMX,JMX),&
                XKLEV(IMX,JMX)
      DIMENSION SLP(IMX,JMX)
      PARAMETER (R=287.04,G=9.8,GAMMA=6.5E-3)
      PARAMETER (TC=273.16+17.5) ! T CRITICAL IN PSFC/PSLV
      PARAMETER (PCONST=100.)
!
      LOGICAL L1,L2,L3,L4
!
!
!
!
!     ... SEA LEVEL PRESSURE
!
      XTERM=GAMMA*R/G
!
!     ... COMPUTE PRESSURE AT PCONST MB ABOVE SURFACE (PL)
!
      KUPTO=KX/2
99    CONTINUE
      DO 100 J=1,JMX
      DO 100 I=1,IMX
         PL(I,J)=SFP(I,J)-PCONST
         XKLEV(I,J)=0.
100   CONTINUE
!
!     ... FIND 2 LEVELS ON SIGMA SURFACES SURROUNDING PL AT EACH I,J
!
      DO 150 J=1,JMX
      DO 150 I=1,IMX
         DO 125 K=KX-1,KUPTO,-1
            XK=FLOAT(K)
            XKHOLD=XKLEV(I,J)
!srf            XKLEV(I,J)=CVMGT(XK,XKHOLD,   &
            XKLEV(I,J)=merge(XK,XKHOLD,   &
              (((PP(I,J,K)).LT.PL(I,J)) .AND.  &
               ((PP(I,J,K+1)).GE.PL(I,J))))
125      CONTINUE
         IF(XKLEV(I,J).LT.1.) THEN
            PRINT *,'ERROR FINDING PRESSURE LEVEL ',PCONST,' MB ',&
                   'ABOVE THE SURFACE'
            PRINT *,'LAST K LEVEL =',KUPTO
            IF(KUPTO.NE.1) THEN
               PRINT *,'TRYING AGAIN WITH KUPTO=1'
               KUPTO=1
               GOTO 99
            ELSE
               PRINT *,'I,J=',I,J
               PRINT *,'PL=',PL(I,J)
               PRINT *,'PSFC=',SFP(I,J)
               STOP
            END IF
         END IF
150   CONTINUE
!
!     ... GET TEMPERATURE AT PL (TL), EXTRAPOLATE T AT SURFACE (TS)
!         AND T AT SEA LEVEL (T0) WITH 6.5 K/KM LAPSE RATE
!
      DO 200 J=1,JMX
      DO 200 I=1,IMX
         KLO=NINT(XKLEV(I,J))+1
         KHI=NINT(XKLEV(I,J))
         PLO=PP(I,J,KLO)
         PHI=PP(I,J,KHI)
         TLO=T(I,J,KLO)
         THI=T(I,J,KHI)
         TL=THI-(THI-TLO)*ALOG(PL(I,J)/PHI)/ALOG(PLO/PHI)
         TS(I,J)=TL*(SFP(I,J)/PL(I,J))**XTERM
         TBAR=(TS(I,J)+TL)*0.5
         HL=TER(I,J)-R/G*ALOG(PL(I,J)/SFP(I,J))*TBAR
         T0(I,J)=TL+GAMMA*HL
200   CONTINUE
!
!     ... CORRECT SEA LEVEL TEMPERATURE IF TOO HOT
!
      DO 400 J=1,JMX
      DO 400 I=1,IMX
         L1=T0(I,J).LT.TC
         L2=TS(I,J).LE.TC
         L3=.NOT.L1
         T0HOLD=T0(I,J)
!srf         T0(I,J)=CVMGT(T0HOLD,&
!srf           CVMGT(TC,TC-0.005*(TS(I,J)-TC)**2,L2.AND.L3),L1.AND.L2)
         T0(I,J)=merge(T0HOLD,&
           merge(TC,TC-0.005*(TS(I,J)-TC)**2,L2.AND.L3),L1.AND.L2)
400   CONTINUE
!
!     ... COMPUTE SEA LEVEL PRESSURE
!
      DO 600 J=1,JMX
      DO 600 I=1,IMX
         SLP(I,J)=SFP(I,J)*EXP(2.*G*TER(I,J)/(R*(TS(I,J)+T0(I,J))))
600   CONTINUE
      RETURN
      END

!-------------------------------------------------------------------------

