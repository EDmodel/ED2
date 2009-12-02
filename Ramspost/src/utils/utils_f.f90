!############################# Change Log ##################################
! 1.0.0.7
!
! 010428 MJB RAMS_getvar ##
!            Fixed to allow use of averaged analysis files. ##
! 010223 MJB rams_f_open ##
!            Only close file if opened. ##
! 001106 MJB RAMS_getvar ##
!            Removed char(0) from filename printing. ##
! 001002 MJB makefnam parsefnam ##
!            Replaced index calls with f90 intrinsic len_trim. ##
! 000921 MJB RAMS_getvar ##
!            Printing filename and variable here instead of in the c
!            routine. ##
! 000829 CJT RAMS_getvar ##
!            Changed from passing arguments to using an_header module. ##
! 000828 MJB RAMS_getvar ##
!            Added " include 'interface.h' " for NT. ##
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!  Mission Research Corporation / *ASTeR Division
!###########################################################################

subroutine aminmax(n1,n2,n3,a)
dimension a(n1,n2,n3)

do k=1,n1
   xmin=1e20
   xmax=-1e20
   do j=1,n3
      do i=1,n2
         xmin=min(xmin,a(k,i,j))
         xmax=max(xmax,a(k,i,j))
      enddo
   enddo
   print*,'Max,min:',k,xmin,xmax
enddo
return
end

!***************************************************************************

subroutine RAMS_mm(indata,ni1,omin,omax)
real indata(ni1),omin,omax

omax=indata(1)
omin=indata(1)
   do i=2,ni1
      omax=max(indata(i),omax)
      omin=min(indata(i),omin)
   enddo

return
end

!***************************************************************************

subroutine RAMS2_mm(m1m,m2m,indata,ni1,vmin,vmax)
real indata(m1m,m2m),vmin,vmax

vmax=-1.E30
vmin=1.E31
   do ini=1,ni1
      if(indata(1,ini).lt.1.e20) then
         vmax=max(indata(1,ini),vmax)
         vmin=min(indata(1,ini),vmin)
      endif
   enddo

return
end

!***************************************************************************

function walltime(wstart)

call system_clock(count=ii,count_rate=ir)
walltime=float(ii)/float(ir) - wstart
return
end

FUNCTION cputime(w1)
  call timing(2,cc)
  cputime=cc
  fsecs=72559200.
  w1=walltime(fsecs)
RETURN
END

!***************************************************************************

subroutine rearrange(nzp,nxp,nyp,a,b)

dimension a(nzp,nxp,nyp),b(nxp,nyp,nzp)

do i=1,nxp
   do j=1,nyp
      do k=1,nzp
         b(i,j,k)=a(k,i,j)
      enddo
   enddo
enddo
return
end

!***************************************************************************

subroutine unarrange(nzp,nxp,nyp,a,b)

dimension a(nxp,nyp,nzp),b(nzp,nxp,nyp)

do i=1,nxp
   do j=1,nyp
      do k=1,nzp
         b(k,i,j)=a(i,j,k)
      enddo
   enddo
enddo
return
end

!***************************************************************************

integer function RAMS_getvar (string,itype,ngrd,a,b,flnm)

use an_header

implicit none
include 'interface.h'
real :: a(*),b(*)
integer :: itype,ngrd
character*(*) flnm,cgrid*1,flng*80,errmsg*120,string
logical there
integer :: ierr_getvar,ifound,ni,npts,iword
common /getvar/ierr_getvar,ifound
integer :: lastchar
write(*,*) 'panis et circenses'
write(*,*) 'nvbtab=',nvbtab
do ni=1,nvbtab

   if((string.eq.anal_table(ni)%string.or.  &
      string//'M'.eq.anal_table(ni)%string).and.  &
      ngrd.eq.anal_table(ni)%ngrid) then
      
      if(string//'M'.eq.anal_table(ni)%string) string=string//'M'
   
      write(cgrid,'(i1)') ngrd
      flng=flnm//'-g'//cgrid//'.vfm'
      print*,' C_open - ',flng(1:len_trim(flng)),'   ',string
      flng=flng(1:len_trim(flng))//char(0)

      inquire(file=flng,exist=there)
      if(.not.there) then
         errmsg='File not found - '//flng
         call error_mess(errmsg)
         return
      endif

      npts=anal_table(ni)%nvalues
      itype=anal_table(ni)%idim_type
      iword=anal_table(ni)%npointer

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

!***************************************************************************

subroutine makefnam (fname,prefix,tinc,iyr,imn,idy,itm,type,post,fmt)

! creates standard timestamped filename

implicit none

integer iyr, imn, idy, itm
integer oyr, omn, ody, otm
integer ib1,ib2
real tinc
character*(*) fname,prefix,post
character dstring*40,fmt*3,type*1

!print*,iyr,imn,idy,itm,tinc
call date_add_to(iyr,imn,idy,itm,tinc,'s',oyr,omn,ody,otm)
!print*,oyr,omn,ody,otm

write(dstring,100) '-',type,'-',oyr,'-',omn,'-',ody,'-',otm
100 format(3a1,i4.4,a1,i2.2,a1,i2.2,a1,i6.6)

ib1=len_trim(prefix)
fname=prefix(1:ib1)//dstring(1:20)
if (len_trim(post).gt.0) then
   ib1=len_trim(fname)
   ib2=len_trim(post)
   fname=fname(1:ib1)//'-'//post(1:ib2)
endif
ib1=len_trim(fname)
fname=fname(1:ib1)//'.'//fmt(1:3)

return
end

!***************************************************************************

subroutine parsefnam (fname,prefix,iyr,imn,idy,itm,type,post,fmt)

! breaks standard timestamped filename into component parts

! MJW - 5/4/00 have to be careful since prefix may include '-'
!    So lets assume that the filename will always have the form
!    prfx_incl_dashs-typ-yyyy-mm-dd-hhmmss[-*].fmt, & look backwards
!    We might have at most 5 dashes within and following datestring
!    if there is a 'post' character string

implicit none
integer iyr, imn, idy, itm
integer ib1,ib2,n,lch,ndash
parameter (ndash=5)
integer idash(ndash)
character*(*) fname,prefix
character dstring*40,post*10,fmt*3,type*1

ib1=index(fname,'.',.TRUE.)
fmt=fname(ib1+1:ib1+3)

lch=len_trim(fname)
do n=1,ndash
   idash(n)=index(fname(1:lch),'-',.TRUE.)
   lch=idash(n)-1
enddo

! Check to see if a post exists by checking the position of -'s.
post=''
if(idash(4)==idash(5)+5 .and. idash(3)==idash(4)+3 .and.  &
   idash(2)==idash(3)+3 .and. idash(1)==idash(2)+7 )  &
   post=fname(idash(1)+1:ib1-1)

if(len_trim(post)>0)idash(1:4)=idash(2:5)

prefix=fname(1:idash(4)-3)
read(fname(idash(4)-1:idash(1)+6),100) type,iyr,imn,idy,itm
100 format(a1,1x,i4,1x,i2,1x,i2,1x,i6)

return
end

!***************************************************************************

subroutine rams_f_open (iunit,filenm,formt,stat,act,iclob)

! replaces old jclopen and jclget
! files are overwritten unless iclob (ICLOBBER) set to 1

implicit none

integer iunit,iclob
character*(*) filenm,formt,stat,act
logical exans,opnd
integer, external :: lastchar

!close(iunit)
inquire(FILE=filenm,EXIST=exans,OPENED=opnd)
!print*,iclob,':',act,':',exans,':',opnd
if(opnd) close(iunit)
if(exans.and.iclob.eq.0.and.  &
     (act(1:5).eq.'WRITE'.or.act(1:5).eq.'write')) then
   print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   print*,'!!!   trying to open file name :'
   print*,'!!!       ',filenm
   print*,'!!!   but it already exists. run is ended.'
   print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   stop 'rams_f_open'
endif

!print*,'filenm,formt,stat= ',filenm(1:lastchar(filenm)),' ',formt,' ',stat
open(iunit,STATUS=stat,FILE=filenm(1:lastchar(filenm)),FORM=formt)

print*,'F_open - ',filenm(1:lastchar(filenm))

return
end

!***************************************************************************

SUBROUTINE JCL
CHARACTER*(*) FILENM,FORMT,FILENL,FILENP
CHARACTER CFNAME*16,TEXTSTR*40
LOGICAL EXANS

!***************************************************************************

ENTRY JCLOPEN(IUNIT,FILENM,FORMT,IPRNT)

! This routine opens a new file with the file name of FILENM
! to write on and assigns it unit number IUNIT with format
! of FORMT (either FORMATTED or UNFORMATTED), first checking
! if file already exists.

INQUIRE(FILE=FILENM,EXIST=EXANS)
IF(EXANS) THEN
   PRINT*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   PRINT*,'!!!   Trying to OPEN file name :'
   PRINT*,'!!!       ',FILENM
   PRINT*,'!!!   But it already exists. Run is ended.'
   PRINT*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   STOP 'JCLOPEN'
ENDIF

OPEN(IUNIT,STATUS='NEW',FILE=FILENM,FORM=FORMT)
RETURN

!***************************************************************************

ENTRY JCLGET(IUNIT,FILENM,FORMT,IPRNT)

! This routine access an existing file with the file name of FILENM
! and assigns it unit number IUNIT.

IF(IPRNT.EQ.1) THEN
PRINT*,' Opening input unit ',IUNIT,' file name ',FILENM
PRINT*,'         format  ',FORMT
ENDIF
OPEN(IUNIT,STATUS='OLD',FILE=FILENM,FORM=FORMT)

RETURN

!***************************************************************************

ENTRY JCLDISP(IUNIT,FILENL,FILENP,IOUT,IPRNT)

! This routine disposes an attached file with the file name of FILENM.
! Disposition can be dependent on the IOUT flag.

CLOSE(IUNIT)

RETURN
END

!***************************************************************************

FUNCTION RAMRAN(idum)

! random number generator with [0,1] uniform distribution 
! by Knuth subtractive method

INTEGER idum
INTEGER MBIG,MSEED,MZ
!REAL MBIG,MSEED,MZ
REAL ran3,FAC
PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1./MBIG)
!PARAMETER (MBIG=4000000.,MSEED=1618033.,MZ=0.,FAC=1./MBIG)
INTEGER i,iff,ii,inext,inextp,k
INTEGER mj,mk,ma(55)
!REAL mj,mk,ma(55)
SAVE iff,inext,inextp,ma
DATA iff /0/
if(idum.lt.0.or.iff.eq.0)then
  iff=1
  mj=MSEED-iabs(idum)
  mj=mod(mj,MBIG)
  ma(55)=mj
  mk=1
  do 11 i=1,54
    ii=mod(21*i,55)
    ma(ii)=mk
    mk=mj-mk
    if(mk.lt.MZ)mk=mk+MBIG
    mj=ma(ii)
  11 continue
  do 13 k=1,4
    do 12 i=1,55
      ma(i)=ma(i)-ma(1+mod(i+30,55))
      if(ma(i).lt.MZ)ma(i)=ma(i)+MBIG
    12 continue
  13 continue
  inext=0
  inextp=31
  idum=1
endif
inext=inext+1
if(inext.eq.56)inext=1
inextp=inextp+1
if(inextp.eq.56)inextp=1
mj=ma(inext)-ma(inextp)
if(mj.lt.MZ)mj=mj+MBIG
ma(inext)=mj
ramran=mj*FAC

return
END







