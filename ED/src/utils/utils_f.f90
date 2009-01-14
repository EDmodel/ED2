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
  if(tinc == 0.d0) then
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
  integer :: ib1,n,lch

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
  logical :: exans

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
