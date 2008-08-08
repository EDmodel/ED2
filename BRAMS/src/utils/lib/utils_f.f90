!############################# Change Log ##################################
! 2.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

subroutine aminmax(n1,n2,n3,a)
  implicit none
  integer :: n1,n2,n3
  real :: a(n1,n2,n3)
  integer :: k,i,j
  real :: xmax,xmin

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
end subroutine aminmax

!***************************************************************************

subroutine RAMS_mm(indata,ni1,omin,omax)
  implicit none
  integer :: ni1
  real :: indata(ni1),omin,omax
  integer :: i

  omax=indata(1)
  omin=indata(1)
  do i=2,ni1
     omax=max(indata(i),omax)
     omin=min(indata(i),omin)
  enddo

  return
end subroutine RAMS_mm

!***************************************************************************

subroutine RAMS2_mm(m1m,m2m,indata,ni1,vmin,vmax)
  implicit none
  integer :: m1m,m2m,ni1
  real :: indata(m1m,m2m),vmin,vmax
  integer :: ini

  vmax=-1.E30
  vmin=1.E31
  do ini=1,ni1
     if(indata(1,ini).lt.1.e20) then
        vmax=max(indata(1,ini),vmax)
        vmin=min(indata(1,ini),vmin)
     endif
  enddo

  return
end subroutine RAMS2_mm

!***************************************************************************

real function walltime(wstart)
  implicit none
  real :: wstart
  integer :: ii,ir

  call system_clock(count=ii,count_rate=ir)
  walltime=float(ii)/float(ir) - wstart
  return
end function walltime

real function cputime(w1)
  implicit none
  real :: w1
  real :: cc,fsecs
  real, external :: walltime

  call timing(2,cc)
  cputime=cc
  fsecs=72559200.
  w1=walltime(fsecs)
  return
end function cputime

!***************************************************************************

subroutine rearrange(nzp,nxp,nyp,a,b)
  implicit none
  integer :: nzp,nxp,nyp
  real :: a(nzp,nxp,nyp),b(nxp,nyp,nzp)
  integer :: k,i,j

  do i=1,nxp
     do j=1,nyp
        do k=1,nzp
           b(i,j,k)=a(k,i,j)
        enddo
     enddo
  enddo
  return
end subroutine rearrange

!***************************************************************************

subroutine unarrange(nzp,nxp,nyp,a,b)
  implicit none
  integer :: nzp,nxp,nyp
  real :: a(nxp,nyp,nzp),b(nzp,nxp,nyp)
  integer :: k,i,j

  do i=1,nxp
     do j=1,nyp
        do k=1,nzp
           b(k,i,j)=a(i,j,k)
        enddo
     enddo
  enddo
  return
end subroutine unarrange

!***************************************************************************

subroutine makefnam (fname,prefix,tinc,iyr,imn,idy,itm,type,post,fmt)

  ! creates standard timestamped filename

  implicit none

  integer :: iyr, imn, idy, itm
  character(len=*) ::  fname,prefix,post
  character(len=*) :: fmt
  character(len=1) :: type

  real(kind=8) :: tinc
  integer :: oyr, omn, ody, otm
  integer :: ib1,ib2
  character(len=40) :: dstring

  !   print*,iyr,imn,idy,itm,tinc
  if(tinc == 0.) then
     oyr=iyr ; omn=imn ; ody=idy ; otm=itm
  else
     call date_add_to(iyr,imn,idy,itm,tinc,'s',oyr,omn,ody,otm)
     !   print*,oyr,omn,ody,otm
  endif

  write(dstring,100) '-',type,'-',oyr,'-',omn,'-',ody,'-',otm
100 format(3a1,i4.4,a1,i2.2,a1,i2.2,a1,i6.6)

  ib1=len_trim(prefix)
  fname=prefix(1:ib1)//dstring(1:20)
  if (post(1:1) /= '$') then
     ib1=len_trim(fname)
     ib2=len_trim(post)
     fname=fname(1:ib1)//'-'//post(1:ib2)
  endif
  ib1=len_trim(fname)
  fname=fname(1:ib1)//'.'//fmt(1:3)

  return
end subroutine makefnam

!***************************************************************************

subroutine parsefnam (fname,prefix,iyr,imn,idy,itm,type,post,fmt)

  ! breaks standard timestamped filename into component parts

  ! MJW - 5/4/00 have to be careful since prefix may include '-'
  !    So lets assume that the filename will always have the form
  !    prfx_incl_dashs-typ-yyyy-mm-dd-hhmmss[-*].fmt, & look backwards
  !    We might have at most 5 dashes within and following datestring
  !    if there is a 'post' character string

  implicit none
  integer :: iyr, imn, idy, itm
  character(len=*) :: fname,prefix
  character :: post*10,fmt*3,type*1

  integer, parameter :: ndash=5
  integer :: idash(ndash)
  character(len=40) :: dstring
  integer :: ib1,ib2,n,lch  !,ndash  !Declare two times

  ib1=index(fname,'.',.true.)
  fmt=fname(ib1+1:ib1+3)

  lch=len_trim(fname)
  do n=1,ndash
     idash(n)=index(fname(1:lch),'-',.true.)
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
end subroutine parsefnam

!***************************************************************************

subroutine rams_f_open(iunit, filenm, formt, stat, act, iclob)

  ! replaces old jclopen and jclget
  ! files are overwritten unless iclob (ICLOBBER) set to 1

  implicit none

  integer :: iunit, iclob
  character(len=*) :: filenm, formt, stat, act
  logical :: exans,opnd

  inquire(FILE=filenm,EXIST=exans)

  if(exans.and.iclob.eq.0.and.  &
       (act(1:4).eq.'WRIT'.or.act(1:4).eq.'writ')) then
     print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     print*,'!!!   trying to open file name :'
     print*,'!!!       ',filenm
     print*,'!!!   but it already exists. run is ended.'
     print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     stop 'rams_f_open - exists'
  endif

  open(iunit,STATUS=stat,FILE=filenm(1:len_trim(filenm)),FORM=formt)
!!$  print*,'F_open - ',filenm(1:len_trim(filenm))

  return
end subroutine rams_f_open

!***************************************************************************

real function RAMRAN(idum)
  implicit none

  ! random number generator with [0,1] uniform distribution 
  ! by Knuth subtractive method

  integer idum
  integer MBIG,MSEED,MZ
  !REAL MBIG,MSEED,MZ
  real ran3,FAC
  parameter (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1./MBIG)
  !PARAMETER (MBIG=4000000.,MSEED=1618033.,MZ=0.,FAC=1./MBIG)
  integer i,iff,ii,inext,inextp,k
  integer mj,mk,ma(55)
  !REAL mj,mk,ma(55)
  save iff,inext,inextp,ma
  data iff /0/
  if(idum.lt.0.or.iff.eq.0)then
     iff=1
     mj=MSEED-iabs(idum)
     mj=mod(mj,MBIG)
     ma(55)=mj
     mk=1
     do i=1,54 !11 i=1,54
        ii=mod(21*i,55)
        ma(ii)=mk
        mk=mj-mk
        if(mk.lt.MZ)mk=mk+MBIG
        mj=ma(ii)
!11   CONTINUE
     enddo
     do k=1,4 !13 k=1,4
        do i=1,55 !12 i=1,55
           ma(i)=ma(i)-ma(1+mod(i+30,55))
           if (ma(i).lt.MZ) ma(i)=ma(i)+MBIG
!12         CONTINUE
        enddo
!13         CONTINUE
     enddo
     inext=0
     inextp=31
     idum=1
  endif
  inext=inext+1
  if (inext.eq.56) inext=1
  inextp=inextp+1
  if (inextp.eq.56) inextp=1
  mj=ma(inext)-ma(inextp)
  if (mj.lt.MZ) mj=mj+MBIG
  ma(inext)=mj
  ramran=mj*FAC

  return
end function RAMRAN








