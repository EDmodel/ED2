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
