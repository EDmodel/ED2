!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

subroutine geodat(n2,n3,datr,hfn,ofn,vt2da,vt2db,ngr,vnam)

use mem_grid
use io_params
use rconstants, only: spcon
use grid_dims, only : str_len

implicit none
integer :: n2,n3,ngr
real ::  vt2da(*),vt2db(*),datr(n2,n3)
character(len=str_len) :: hfn,ofn,title
character(len=3) :: vnam

integer :: lb,iblksizo,no,isbego,iwbego,iodim,mof,niq,njq,np
real :: offlat,offlon,deltallo,deltaxq,deltayq,deltaxp,deltayp,erad

real,allocatable:: dato(:)

LB=len_trim(HFN)
if(LB.le.0) then
   print*,'==================================================='
   print*,'|  Problem in GEODAT, Input data prefix incorrect !'
   print*,'|  Grid :',ngrid
   print*,'|  File prefix:',HFN
   print*,'==================================================='
   stop 'GEODAT-file'
endif

if((vnam(1:2).eq.'TO'.or.vnam(1:2).eq.'ZO').and.  &
   (ITOPSFLG(NGR).eq.1.and.TOPTENH(NGR).gt.1.)) then
   print*,'==================================================='
   print*,'|  Problem in GEODAT, Silhouette weight too high !'
   print*,'|  Grid :',NGR
   print*,'|  Weight (range TOPTENH=0-1):',TOPTENH(NGR)
   print*,'==================================================='
   stop 'GEODAT'
endif

!     configure grid specs for raw data to rams grid (R) transfer

!     raw data grid (O)
if(vnam.eq.'TOD'.or.vnam.eq.'ZOD') then  !using dted data
   iblksizo=5 ! 5 degree squares
   no=600     ! 30" data over a 5 degree square
   isbego=-90
   iwbego=0
   offlat=0.
   offlon=0.
   DELTALLO=FLOAT(IBLKSIZO)/FLOAT(NO)
else
   TITLE=HFN(1:LB)//'HEADER'
   LB=len_trim(TITLE)
   call rams_f_open(29,title(1:lb),'FORMATTED','OLD','READ',0)
   READ(29,*)IBLKSIZO,NO,ISBEGO,IWBEGO,offlat,offlon
   CLOSE(29)
   DELTALLO=FLOAT(IBLKSIZO)/FLOAT(NO-1)
endif

iodim=max(100000,4*no*no)
MOF=IODIM/(NO*NO)

   allocate(dato(iodim+mof+mof))

!     temp grid (Q) - smoothing only applied to topo
if(vnam(1:2).eq.'TO') then
      DELTAXQ=0.5*TOPTWVL(NGR)*DELTAXN(NGR)
   DELTAYQ=0.5*TOPTWVL(NGR)*DELTAYN(NGR)
else
      DELTAXQ=DELTAXN(NGR)
   DELTAYQ=DELTAYN(NGR)
endif
NIQ=INT(FLOAT(NNXP(NGR)-1)*DELTAXN(NGR)/DELTAXQ)+4
NJQ=INT(FLOAT(NNYP(NGR)-1)*DELTAYN(NGR)/DELTAYQ)+4

!     interpollated raw data grid (P)
NP=MIN(10,MAX(1,INT(DELTAXQ/(DELTALLO*spcon))))
DELTAXP=DELTAXQ/FLOAT(NP)
DELTAYP=DELTAYQ/FLOAT(NP)

CALL SFCOPQR(NO,MOF,NP,NIQ,NJQ,N2,N3,XTN(1,NGR),YTN(1,NGR)  &
     ,platn(ngr),plonn(ngr)  &
     ,ERAD,DELTALLO,DELTAXP,DELTAYP,DELTAXQ,DELTAYQ,IBLKSIZO  &
     ,ISBEGO,IWBEGO,DATO(1),VT2DA,VT2DB,DATR  &
     ,OFN,offlat,offlon,VNAM,NGR,itopsflg(ngr),iz0flg(ngr))

deallocate(dato)

RETURN
END

!     ******************************************************************

subroutine sfcopqr(no,mof,np,niq,njq,n2,n3,xt,yt,platn,plonn  &
     ,erad,deltallo,deltaxp,deltayp,deltaxq,deltayq,iblksizo  &
     ,isbego,iwbego,dato,datp,datq,datr  &
     ,ofn,offlat,offlon,vnam,ngr,itopsflg,iz0flg)

use teb_spm_start, only: TEB_SPM
use grid_dims, only : str_len

implicit none
integer :: no,mof,np,niq,njq,n2,n3,iblksizo,isbego,iwbego,ngr  &
          ,itopsflg,iz0flg
real :: dato(no,no,mof),datp(np,np),datq(niq,njq),datr(n2,n3)  &
         ,xt(n2),yt(n3)
real :: erad,deltallo,deltaxp,deltayp,deltaxq,deltayq,offlat,offlon  
character(len=str_len) :: ofn,title3
character(len=3) :: title1,vnam
character(len=4) :: title2
logical l1,l2
integer,parameter :: maxmiss=1000
character(len=str_len) :: fnmiss(maxmiss)
real, allocatable :: sdq(:,:),shaq(:,:),sdr(:,:),datre(:,:)
real, allocatable :: iso(:),iwo(:)

real :: platn,plonn,xcentr,ycentr,glatp,glonp,rio_full,rjo_full  &
       ,xq,yq,xp,yp,wio1,wio2,wjo1,wjo2,sha,rha,rh2,sh,rh &
       ,xq1,yq1,xr,yr,rval,diff,difflcl
integer :: nmiss,nono,nofr,iof,iq,jq,ip,jp,iwoc,isoc,io1,io2,jo1,jo2 &
          ,lb,nn,isocpt,isocpo,iwocpo,iwocph,iwocpt,io_full,jo_full &
          ,iofr,jofr,ir,jr,is,js,i,j



allocate (sdq(niq,njq),shaq(niq,njq),sdr(n2,n3),datre(n2,n3))
allocate (iso(mof),iwo(mof))

nmiss=0

nono=no*no
XCENTR=0.5*(XT(1)+XT(N2))
YCENTR=0.5*(YT(1)+YT(N3))
NOFR=0
DO IOF=1,MOF
   ISO(IOF)=0
   IWO(IOF)=0
ENDDO
DO JQ=1,NJQ
   DO IQ=1,NIQ
      XQ=(FLOAT(IQ)-0.5*FLOAT(NIQ+1))*DELTAXQ+XCENTR
      YQ=(FLOAT(JQ)-0.5*FLOAT(NJQ+1))*DELTAYQ+YCENTR
      DO JP=1,NP
         DO IP=1,NP
            XP=XQ+(FLOAT(IP)-0.5*FLOAT(NP+1))*DELTAXP
            YP=YQ+(FLOAT(JP)-0.5*FLOAT(NP+1))*DELTAYP

            call xy_ll(GLATP,GLONP,platn,plonn,xp,yp)

            glatp = max(-89.9999,min(89.9999,glatp - offlat))
            glonp = glonp - offlon

            if (glonp >=  180.) glonp = glonp - 360.
            if (glonp <= -180.) glonp = glonp + 360.

            rio_full = (glonp - float(iwbego)) / deltallo
            rjo_full = (glatp - float(isbego)) / deltallo

            io_full = int(rio_full)
            jo_full = int(rjo_full)

            iwoc = (io_full / (no-1)) * iblksizo + iwbego
            isoc = (jo_full / (no-1)) * iblksizo + isbego

            wio2 = rio_full - float(io_full)
            wjo2 = rjo_full - float(jo_full)
           
            wio1 = 1. - wio2
            wjo1 = 1. - wjo2

            io1 = mod(io_full,no-1) + 1
            jo1 = mod(jo_full,no-1) + 1

            io2 = io1 + 1
            jo2 = jo1 + 1
            
            DO IOFR=1,NOFR
               JOFR=IOFR
               IF(ISO(IOFR).EQ.ISOC.AND.IWO(IOFR).EQ.IWOC)GO TO 10
            ENDDO

!                  not using dted data
            if(vnam.ne.'TOD'.and.vnam.ne.'ZOD') then
               ISOCPT=ABS(ISOC)/10
               ISOCPO=ABS(ISOC)-ISOCPT*10
               IWOCPH=ABS(IWOC)/100
               IWOCPT=(ABS(IWOC)-IWOCPH*100)/10
               IWOCPO=ABS(IWOC)-IWOCPH*100-IWOCPT*10
               IF(ISOC.GE.0) THEN
                  WRITE(TITLE1,'(2I1,A1)')ISOCPT,ISOCPO,'N'
               ELSE
                  WRITE(TITLE1,'(2I1,A1)')ISOCPT,ISOCPO,'S'
               ENDIF
               IF(IWOC.GE.0) THEN
                  WRITE(TITLE2,'(3I1,A1)')IWOCPH,IWOCPT,IWOCPO,'E'
               ELSE
                  WRITE(TITLE2,'(3I1,A1)')IWOCPH,IWOCPT,IWOCPO,'W'
               ENDIF
               LB=len_trim(OFN)
               TITLE3=OFN(1:LB)//TITLE1//TITLE2
               LB=len_trim(TITLE3)
               INQUIRE(FILE=TITLE3(1:LB),EXIST=L1,OPENED=L2)

               IF(.NOT.L1)THEN
                  do nn=1,nmiss
                     if(TITLE3(1:LB).eq.fnmiss(nn)) goto 302
                  enddo
                  nmiss=nmiss+1
                  fnmiss(nmiss)=TITLE3(1:LB)
302                    continue
                  DATP(IP,JP)=0.
                  GOTO 20
               ENDIF
            ENDIF

            IF(NOFR.GE.MOF) THEN
               DO IOF=1,MOF
                  ISO(IOF)=0
                  IWO(IOF)=0
               ENDDO
               NOFR=0
            ENDIF
            NOFR=NOFR+1
            JOFR=NOFR

!                 using dted data
            if(vnam.eq.'TOD'.or.vnam.eq.'ZOD') then
               call dted(no,ofn,isoc,iwoc,dato(1,1,nofr))
            else
               call rams_f_open  &
                    (29,TITLE3(1:LB),'FORMATTED','OLD','READ',0)
               CALL VFIREC(29,DATO(1,1,NOFR),NONO,'LIN')
               CLOSE(29)
            endif

            ISO(NOFR)=ISOC
            IWO(NOFR)=IWOC

10          CONTINUE

            datp(ip,jp)=wio1*(wjo1*dato(io1,jo1,jofr)   &
                             +wjo2*dato(io1,jo2,jofr))  &
                       +wio2*(wjo1*dato(io2,jo1,jofr)   &
                             +wjo2*dato(io2,jo2,jofr))
                                                          
20          CONTINUE
         ENDDO
      ENDDO

!           std dev for envelope orog and topo based zo schemes
      SHA=0.
      RHA=0.
      RH2=0.
      DO JP=1,NP
         SH=0.
         RH=0.
         thisloop: DO IP=1,NP
            !------------------------------------------------------------------------------!
            !   No, this doesn't make any sense but if I don't put this "cycle" it gives   !
            ! floating point exception when it attempts to compute the max and it crashes. !
            !------------------------------------------------------------------------------!
            if (datp(ip,jp) < 1.e-16 ) cycle thisloop
            !------------------------------------------------------------------------------!

            SH=MAX(SH,DATP(IP,JP))
            RH=RH+DATP(IP,JP)
            RH2=RH2+DATP(IP,JP)**2
         END DO thisloop

         SHA=SHA+SH/(2.*FLOAT(NP))
         RHA=RHA+RH
      ENDDO
      DATQ(IQ,JQ)=RHA/FLOAT(NP*NP)
      SDQ(IQ,JQ)=SQRT(max(0.,RH2/NP**2-DATQ(IQ,JQ)**2))
      DO IP=1,NP
         SH=0.
         DO JP=1,NP
            SH=MAX(SH,DATP(IP,JP))
         ENDDO
         SHA=SHA+SH/(2.*FLOAT(NP))
      ENDDO
      SHAQ(IQ,JQ)=SHA
      
   ENDDO
!         print*,'finished sfcopqr row jq = ',jq
ENDDO

!     envelope and zo schemes

if((vnam.eq.'TOP'.or.vnam.eq.'TOD').and.  &
   (ITOPSFLG.eq.2.or.ITOPSFLG.eq.3).and.  &
   NP*NP.lt.8) print*,'Warning - '  &
   ,'trying to calc a std dev for: ',NP*NP,' points'
if((vnam.eq.'ZOT'.or.vnam.eq.'ZOD').and.  &
   IZ0FLG.eq.1.and.NP*NP.lt.8) print*,'Warning - '  &
   ,'trying to calc a std dev for: ',NP*NP,' points'

if(vnam.eq.'TOP'.or.vnam.eq.'TOD')  &
   CALL TOPOQ(NIQ,NJQ,DELTAXQ,DELTAYQ,DATQ,SDQ,SHAQ,DATRE  &
           ,NGR,N2,N3)
if(vnam.eq.'ZOT'.or.vnam.eq.'ZOD')  &
   CALL ZOQ(NIQ,NJQ,DATQ,SDQ,NGR)

XQ1=(1.-0.5*FLOAT(NIQ+1))*DELTAXQ+XCENTR
YQ1=(1.-0.5*FLOAT(NJQ+1))*DELTAYQ+YCENTR
DO JR=1,N3
   DO IR=1,N2
      XR=(XT(IR)-XQ1)/DELTAXQ+1.
      YR=(YT(JR)-YQ1)/DELTAYQ+1.
      CALL GDTOST(DATQ,NIQ,NJQ,XR,YR,RVAL)

!           envelope orog and zo schemes

      if(vnam.eq.'ZOT'.or.vnam.eq.'ZOD') then
         DATR(IR,JR)=MAX(0.,RVAL)
!               print*,'z0r',IR,JR,DATR(IR,JR)
      else
         if (TEB_SPM==1) then
            if(vnam.ne.'FUS')then
               DATR(IR,JR)=MAX(0.,RVAL)
            else
               datr(ir,jr)=rval
            endif
         else
            DATR(IR,JR)=MAX(0.,RVAL)
         endif
      endif
      
      CALL GDTOST(SDQ,NIQ,NJQ,XR,YR,RVAL)
      SDR(IR,JR)=MAX(0.,RVAL)
   ENDDO
ENDDO

if(nmiss.gt.0) then
   print*,'-----------------------------------------------------'
   print*,'Input physiographical data file processing:'
   print*,'-----------------------------------------------------'
   print*,'  Input data blocks not found '  &
        ,' (data assumed to be zero):'
   do nn=1,nmiss
      print*,trim(fnmiss(nn))
   enddo
   print*,'-----------------------------------------------------'
endif

!     check to find the largest change in topo height

if(vnam.eq.'TOP'.or.vnam.eq.'TOD') then
   diff=0.
   difflcl=0.
   is=-999
   js=-999
   do j=2,n3-1
      do i=2,n2-1
!               print*,'1=',difflcl,max(difflcl,abs(datr(i,j)-datr(i-1,j)))
         difflcl=max(difflcl,abs(datr(i,j)-datr(i-1,j)))
!               print*,'2=',difflcl,max(difflcl,abs(datr(i,j)-datr(i+1,j)))
         difflcl=max(difflcl,abs(datr(i,j)-datr(i+1,j)))
!               print*,'3=',difflcl,max(difflcl,abs(datr(i,j)-datr(i,j-1)))
         difflcl=max(difflcl,abs(datr(i,j)-datr(i,j-1)))
!               print*,'4=',difflcl,max(difflcl,abs(datr(i,j)-datr(i,j+1)))
         difflcl=max(difflcl,abs(datr(i,j)-datr(i,j+1)))
         if(abs(diff-difflcl).gt.1.) then
            is=i
            js=j
         endif
         diff=max(diff,difflcl)
      enddo
   enddo
   write(6,100) ' Max d(topo) on grid @i,j=',ngr,is,js,diff
100     format(a,3i4,f8.1)
endif

deallocate(SDQ,SHAQ,SDR,DATRE)
deallocate(ISO,IWO)

RETURN
END

!**********************************************************************

subroutine topoq(niq,njq,deltaxq,deltayq,datq,sdq,shaq,datre  &
                ,ngr,n2,n3)
                
use io_params
                
implicit none
integer :: niq,njq,ngr,n2,n3
real :: deltaxq,deltayq
real :: datq(niq,njq),sdq(niq,njq),shaq(niq,njq),datre(n2,n3)

integer :: iq,jq,jmin,imin,ire,jre,imax,jmax
real :: rad,count,total,remax,remin,average

!     orographic schemes

if(ITOPSFLG(ngr).lt.0) then                         ! No orography
   do jq=1,njq
      do iq=1,niq
         datq(iq,jq)=0.
!            print*,'None',iq,jq,datq(iq,jq)
      enddo
   enddo
   print *,'No orography'

elseif(ITOPSFLG(ngr).lt.0) then                      ! Average
   print *,'No orography enhancement applied'

   elseif(ITOPSFLG(ngr).eq.1) then                   ! Silhouette
      do jq=1,njq
         do iq=1,niq
            datq(iq,jq)=SHAQ(IQ,JQ)*toptenh(ngr)  &
                       +DATQ(IQ,JQ)*(1.-toptenh(ngr))
!                  print*,'Silhouette',iq,jq,datq(iq,jq)
         enddo
      enddo
      print *,'Silhouette Orography applied with'
      print *,'weighting = ',toptenh(ngr)

   elseif(ITOPSFLG(ngr).eq.2) then                   ! Envelope
      do jq=1,njq
         do iq=1,niq
            datq(iq,jq)=datq(iq,jq)+toptenh(ngr)*sdq(iq,jq)
!                  print*,'EO',iq,jq,datq(iq,jq)
         enddo
      enddo
      print *,'Envelope Orography applied with'
      print *,'enhancement = ',toptenh(ngr),' x std dev'

   else if(ITOPSFLG(ngr).ge.3) then                  ! Reflected Envelope

!        the radius we want to search for the current pts relative
!        height should correspond well to half the filtering wavelength
!        used on the topo (toptwvl)

         Rad=toptwvl(ngr)/2
      do jq=1,njq
         do iq=1,niq
            datre(iq,jq)=datq(iq,jq)
         enddo
      enddo
      do jq=1,njq
         do iq=1,niq
            count=0.
            total=0.
            remax=datre(iq,jq)
            remin=datre(iq,jq)
            jmin=jq-nint(Rad)
            imin=iq-nint(Rad)
            jmax=jq+nint(Rad)
            imax=iq+nint(Rad)
            do jre=max(1,jmin),min(njq,jmax)
               do ire=max(1,imin),min(niq,imax)
                  if((float((iq-ire)))**2  &
                    +(float((jq-jre)))**2.le.Rad**2) then
                   count=count+1.
                  total=total+datre(ire,jre)
                  remax=max(remax,datre(ire,jre))
                  remin=min(remin,datre(ire,jre))
               endif
            enddo
            enddo
         average=total/count
         if(remax.ne.remin)  &
                  datq(iq,jq)=datre(iq,jq)+(datre(iq,jq)-average)/  &
            ((remax-remin)/2)*toptenh(ngr)*sdq(iq,jq)
!               print*,'REO',iq,jq,datre(iq,jq),sdq(iq,jq),datq(iq,jq)
!               print*,'avg,n',average,count,remax,remin
      enddo
   enddo
   print *,'Reflected Envelope Orography applied with'
   print *,'enhancement = ',toptenh(ngr),' x std dev'
   print *,'and search radius (grid points) = ',Rad
endif

RETURN
END

!**********************************************************************

SUBROUTINE ZOQ(NIQ,NJQ,DATQ,SDQ,NGR)

use io_params
use mem_leaf

implicit none
integer :: NIQ,NJQ,NGR
real :: DATQ(NIQ,NJQ),SDQ(NIQ,NJQ)

integer :: iq,jq

!     topo base roughness length.


do jq=1,njq
   do iq=1,niq
      if(ITOPSFLG(ngr).lt.0) then  ! No orography
         datq(iq,jq)=zrough
      elseif(iz0flg(ngr).eq.1) then
         datq(iq,jq)=min(z0fact*sdq(iq,jq),z0max(NGR))
      else
         datq(iq,jq)=zrough
      endif
!            print*,'z0',iq,jq,datq(iq,jq)
   enddo
enddo
if(ITOPSFLG(ngr).lt.0) then  ! No orography
   print *,'No orography'
else
   print *,'Subgrid terrain roughness applied with'
   print *,'factor  = ',z0fact,' x std dev'
   print *,'maximum = ',z0max(NGR)
endif

RETURN
END

!**********************************************************************

subroutine dted(no,pathname,lat,lon,dato)
use grid_dims, only : str_len
implicit none
integer :: no,lat,lon
real :: dato(no,no)
character(len=str_len) :: pathname

!     Let's try and bypass all the bookeeping and just read the file
!     Note that the latitude bands are 5 degrees less than those
!     specified in the original code from Sarma since we are not
!     considering the max latitude but rather the start latitude.

character(len=str_len) :: fname
integer :: ifact,notfnd

ifact = 6
if(lat .ge. 0)then
   if(lat .le. 75) ifact = 4
   if(lat .le. 70) ifact = 3
   if(lat .le. 65) ifact = 2
   if(lat .le. 45) ifact = 1
else
   ifact = 1
   if(lat .le. -45) ifact = 2
   if(lat .le. -65) ifact = 3
   if(lat .le. -70) ifact = 4
   if(lat .le. -75) ifact = 6
end if
call dtedint(no,ifact,lon,lat,notfnd,pathname,dato)

return
end

! ----------------------------------------------------------------------

subroutine dtedint(no,iwres,lon,lat,notfnd,pathname,dato)
use grid_dims, only : str_len
implicit none
integer :: no,iwres,lon,lat,notfnd
real :: dato(no,no)
character*(*) pathname

!     This routine (call by readdtdt) determines the filename
!     of the dted terrain file, based on the latitude and longitude,
!     and calls the c routine that reads the file.
!     uses shell script to read a data set that has been compressed
!     Most recent development Paul Boris.


!     new input data holding arrays:

real readin2(360000)
character(len=str_len) ::  newname1
character(len=5) :: fmtstr
character*3  degns,degew
character*12 dtedfile
character*4 subdir
character(len=str_len) :: extdted,callname
real :: dtedwk(600,600)
integer :: no_blanks,i,num_chrs,ifile,lc,isave,ndtx,ndty,ik,jk,ibtre,j,id
real :: rvaln,wt
real, external :: readdted1

include 'interface.h'

no_blanks = 1
do i = 1, len(pathname)
   if(pathname(i:i) .ne. ' ') no_blanks = i
enddo

num_chrs = no_blanks + 17
write(fmtstr,'(a2,i2,a1)') '(a',num_chrs,')'

write(degew,'(i3)') abs(lon)
if( degew(1:1) .eq. ' ') degew(1:1) = '0'
if( degew(2:2) .eq. ' ') degew(2:2) = '0'
write(degns,'(i3)') abs(lat)
if( degns(1:1) .eq. ' ') degns(1:1) = '0'
if( degns(2:2) .eq. ' ') degns(2:2) = '0'
if(lon .lt. 0) then
   dtedfile(5:5) = 'w'
else
   dtedfile(5:5) = 'e'
endif
dtedfile(6:8) = degew(1:3)
if(lat .lt. 0) then
   dtedfile(1:1) = 's'
else
   dtedfile(1:1) = 'n'
endif
dtedfile(2:4) = degns(1:3)

!     Change from filenames .mvk.gz <---> .mgz
!      dtedfile(9:12) = '.mgz'
dtedfile(9:12) = '.mvk'

!     open the dted data file and read the terrain data:

!      print*, 'OPENING DTED DATA FILE: ',dtedfile
subdir='/'//dtedfile(1:1)//dtedfile(5:5)//'/'
newname1 = pathname(1:no_blanks)//subdir

ifile=42
open(ifile,file='namefils.out')
write(ifile,fmt='(a12)') dtedfile
close (ifile)

!     Call all the decompress stuff from this Fortran code.
lc=len_trim(newname1)

!     Change from filenames .mvk.gz <---> .mgz
!      extdted=newname1(1:lc)//dtedfile
extdted=newname1(1:lc)//dtedfile//'.gz'

lc=len_trim(extdted)
write(newname1,1002)extdted(1:lc)
1002 format('cp ',a,' tmp.gz')
print*,newname1(1:index(newname1,'tmp.gz')+5)
call system(newname1)
call system('chmod ugo+rw tmp.gz')
call system('gzip -d tmp.gz')
write(newname1,1003)dtedfile
1003 format('mv tmp ',a)
call system(newname1)

call azero(no*no,dato)

isave = iwres
rvaln = READDTED1(iwres,readin2)
if(iwres .eq. -100000) then
   notfnd = 1
   return
endif
iwres = isave

ndtx = 600/iwres
ndty = 600

do jk=1,ndty
   do ik=1,ndtx
     ibtre = (jk-1)*ndtx + ik
     dtedwk(ik,jk) = readin2(ibtre)
   enddo
enddo

write(newname1,1004)dtedfile
1004 format('rm -f ',a)
call system(newname1)

do j=1,no
   do i=1,no-iwres+1
     id=(i+iwres-1)/iwres
     wt=mod(real(i+iwres-1),real(iwres))/real(iwres)
     dato(i,j)=(1.-wt)*dtedwk(id,j)+wt*dtedwk(id+1,j)
   enddo
enddo

return
end
