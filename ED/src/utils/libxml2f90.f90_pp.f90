
!===========================================================================
! Copyright (C) 2005, Alexander Poddey <alexander.poddey@gmx.net>
!===========================================================================
! for bugreports please use subject libxml2f90:bug
! for feature requests please use subject libxml2f90:feature
!===========================================================================
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!===========================================================================!
!===========================================================================
! This General Public License does *NOT* permit incorporating this 
! program/code into ->proprietary<- programs.
! If you want to do so, please contact the author!
!===========================================================================
!===========================================================================

!VERSION 2.00

!====================================
!=========== TO DO ==================
!
!



module libxml2f90_module
  integer(4)                  ::  nfilmax=424242  !
  integer(4)                  ::  nfilstart=4242
  integer(4)                  ::  nfilmaxused=4242
  character(256), allocatable ::  stringa(:)
  character(1),allocatable    ::  readstring(:)
  character(1),allocatable    ::  tempstringa(:)
  integer(4)                  ::  filelines
  integer(4)                  ::  lbact
  character(32)               ::  default_llid='CNTL'
  integer(4)                  ::  xmlformat=3
  integer(4)                  ::  arraystep=2000
  integer(4)                  ::  indstep=2
  integer(4),allocatable      ::  lineposa(:)
  logical(4)                  ::  ttransform_paw=.false.
  logical(4)                  ::  tpaw=.false.
  logical(4)                  ::  twrite_paw=.false.
  logical(4)                  ::  trmquotes=.false.
  logical(4)                  ::  trmcomma=.false.

  !0=xml without bl. rem., 1=xml with blank removal, 
  !2=flexible without blank removal, 3=flex with bl. remov.
end module libxml2f90_module

module libxml2f90_strings_module
  PUBLIC
  INTERFACE OPERATOR (-)
     MODULE PROCEDURE lowercase
  END INTERFACE
  INTERFACE OPERATOR (+) 
     MODULE PROCEDURE uppercase
  END INTERFACE

  INTERFACE OPERATOR (.itos.)
     MODULE PROCEDURE i_2_s
  END INTERFACE
  INTERFACE OPERATOR (.rtos.)
     MODULE PROCEDURE r_2_s
  END INTERFACE
contains
  function lowercase(old) result(new)
    implicit none
    character(*), intent(in):: old
    character(len(old))     :: new
    integer(4)              :: i,isvar
    !.......................................
    new=old 
    do i=1,len(trim(old))
       isvar=iachar(old(i:i))
       if(isvar.ge.65.and.isvar.le.90) new(i:i)=achar(isvar+32)
    enddo
  end function lowercase

  function uppercase(old) result(new)
    implicit none
    character(*), intent(in):: old
    character(len(old))     :: new
    integer(4)              :: i,isvar
    !............................................ 
    new=old 
    do i=1,len(trim(old))
       isvar=iachar(old(i:i))
       if(isvar.ge.97.and.isvar.le.122) new(i:i)=achar(isvar-32)
    enddo
  end function uppercase

  function i_2_s(i) result(string)
    implicit none
    integer(4),intent(in)        :: i
    character(256)               :: string
    integer(4)                   :: ibase
    real(8)                      :: rbase=10_8
    integer(4)                   :: ibasemax=10
    integer(4)                   :: ifl
    real(8)                      :: rbj
    real(8)                      :: iact
    integer(4)                   :: iwrite
    logical(4)                   :: tzero
    
    string=''
    iwrite=1
    iact=abs(real(i,kind=8))
    tzero=.true.
    
    !checks
    if(dble(i).gt.rbase**(ibasemax-1)) then
       print*, 'ERROR STOP IN INTEGER TO STRING'
       print*, 'YOUR INTEGER NUMBER IS LARGER THAN IMPLEMENTED'
    end if
    
    if(i.lt.0) then
       string(1:1)='-'
       iwrite=2
    end if
    
    do ibase=1,ibasemax 
       rbj=rbase**(ibasemax-ibase)
       ifl=floor(iact/rbj)
       iact=iact-real(ifl,kind=8)*rbj
       
       select case (ifl)
       case (0)
          if(tzero) then !we do not keep leading zeros
             if(rbj.eq.1.d0) then !we keep the 0 for 10**0
                string(iwrite:iwrite)='0'
                iwrite=iwrite+1
             end if
          else  !we have a value .neq. 0
             string(iwrite:iwrite)='0'
             iwrite=iwrite+1
          end if
       case (1)
          string(iwrite:iwrite)='1'
          iwrite=iwrite+1
          tzero=.false.
       case (2)
          string(iwrite:iwrite)='2'
          iwrite=iwrite+1
          tzero=.false.
       case (3)
          string(iwrite:iwrite)='3'
          iwrite=iwrite+1
          tzero=.false.
       case (4)
          string(iwrite:iwrite)='4'
          iwrite=iwrite+1
          tzero=.false.
       case (5)
          string(iwrite:iwrite)='5'
          iwrite=iwrite+1
          tzero=.false.
       case (6)
          string(iwrite:iwrite)='6'
          iwrite=iwrite+1
          tzero=.false.
       case (7)
          string(iwrite:iwrite)='7'
          iwrite=iwrite+1
          tzero=.false.
       case (8)
          string(iwrite:iwrite)='8'
          iwrite=iwrite+1
          tzero=.false.
       case (9)
          string(iwrite:iwrite)='9'
          iwrite=iwrite+1
          tzero=.false.
       end select
    end do
  end function i_2_s

  function r_2_s(i) result(string)
    implicit none
    real(8),intent(in)           :: i
    character(256)               :: string
    integer(4)                   :: ibase
    real(8)                      :: rbase=10.d0
    integer(4)                   :: ibasemax=10
    integer(4)                   :: digits=13
    integer(4)                   :: j,ifl,izero,k
    real(8)                      :: rbj
    real(8)                      :: iact
    integer(4)                   :: iwrite
    logical(4)                   :: tzero,tfzero
    
    string=''
    iwrite=1
    iact=abs(i)
    tzero=.true.
    
    !checks
    if(i.gt.rbase**(ibasemax-1)) then
       print*, 'ERROR STOP IN INTEGER TO STRING'
       print*, 'YOUR INTEGER NUMBER IS LARGER THAN IMPLEMENTED'
    end if
    
    if(i.lt.0.d0) then
       string(1:1)='-'
       iwrite=2
    end if
    
    do ibase=1,ibasemax 
       rbj=rbase**(ibasemax-ibase)
       ifl=floor(iact/rbj)
       iact=iact-real(ifl,kind=8)*rbj
       
       select case (ifl)
       case (0)
          if(tzero) then !we do not keep leading zeros
             if(rbj.eq.1.d0) then !we keep the 0 for 10**0
                string(iwrite:iwrite)='0'
                iwrite=iwrite+1
             end if
          else  !we have a value .neq. 0
             string(iwrite:iwrite)='0'
             iwrite=iwrite+1
          end if
       case (1)
          string(iwrite:iwrite)='1'
          iwrite=iwrite+1
          tzero=.false.
       case (2)
          string(iwrite:iwrite)='2'
          iwrite=iwrite+1
          tzero=.false.
       case (3)
          string(iwrite:iwrite)='3'
          iwrite=iwrite+1
          tzero=.false.
       case (4)
          string(iwrite:iwrite)='4'
          iwrite=iwrite+1
          tzero=.false.
       case (5)
          string(iwrite:iwrite)='5'
          iwrite=iwrite+1
          tzero=.false.
       case (6)
          string(iwrite:iwrite)='6'
          iwrite=iwrite+1
          tzero=.false.
       case (7)
          string(iwrite:iwrite)='7'
          iwrite=iwrite+1
          tzero=.false.
       case (8)
          string(iwrite:iwrite)='8'
          iwrite=iwrite+1
          tzero=.false.
       case (9)
          string(iwrite:iwrite)='9'
          iwrite=iwrite+1
          tzero=.false.
       end select
    end do

    string(iwrite:iwrite)='.'
    iwrite=iwrite+1
    
    tzero=.false.
    tfzero=.false.
    izero=0
    do j=1,digits
       rbj=10.d0**(-j)
       ifl=floor(iact/rbj)
       iact=real(iact,kind=8)-real(ifl,kind=8)*real(rbj,kind=8)
       
       select case (ifl)
       case (0)
          if(tzero) then !we have a subsequent zero
             tfzero=.true.
             izero=izero+1
          else  !we have the first zero
             string(iwrite:iwrite)='0'
             iwrite=iwrite+1
             tzero=.true.
          end if
       case (1)
          if(tfzero) then
             !write the foregoing zeros
             do k=1,izero
                string(iwrite:iwrite)='0'
                iwrite=iwrite+1
             end do
             izero=0
             tfzero=.false.
             tzero=.false.
          end if
          string(iwrite:iwrite)='1'
          iwrite=iwrite+1
          tzero=.false.
       case (2)
          if(tfzero) then
             !write the foregoing zeros
             do k=1,izero
                string(iwrite:iwrite)='0'
                iwrite=iwrite+1
             end do
             izero=0
             tfzero=.false.
             tzero=.false.
          end if
          string(iwrite:iwrite)='2'
          iwrite=iwrite+1
          tzero=.false.
       case (3)
          if(tfzero) then
             !write the foregoing zeros
             do k=1,izero
                string(iwrite:iwrite)='0'
                iwrite=iwrite+1
             end do
             izero=0
             tfzero=.false.
             tzero=.false.
          end if
          string(iwrite:iwrite)='3'
          iwrite=iwrite+1
          tzero=.false.
       case (4)
          if(tfzero) then
             !write the foregoing zeros
             do k=1,izero
                string(iwrite:iwrite)='0'
                iwrite=iwrite+1
             end do
             izero=0
             tfzero=.false.
             tzero=.false.
          end if
          string(iwrite:iwrite)='4'
          iwrite=iwrite+1
          tzero=.false.
       case (5)
          if(tfzero) then
             !write the foregoing zeros
             do k=1,izero
                string(iwrite:iwrite)='0'
                iwrite=iwrite+1
             end do
             izero=0
             tfzero=.false.
             tzero=.false.
          end if
          string(iwrite:iwrite)='5'
          iwrite=iwrite+1
          tzero=.false.
       case (6)
          if(tfzero) then
             !write the foregoing zeros
             do k=1,izero
                string(iwrite:iwrite)='0'
                iwrite=iwrite+1
             end do
             izero=0
             tfzero=.false.
             tzero=.false.
          end if
          string(iwrite:iwrite)='6'
          iwrite=iwrite+1
          tzero=.false.
       case (7)
          if(tfzero) then
             !write the foregoing zeros
             do k=1,izero
                string(iwrite:iwrite)='0'
                iwrite=iwrite+1
             end do
             izero=0
             tfzero=.false.
             tzero=.false.
          end if
          string(iwrite:iwrite)='7'
          iwrite=iwrite+1
          tzero=.false.
       case (8)
          if(tfzero) then
             !write the foregoing zeros
             do k=1,izero
                string(iwrite:iwrite)='0'
                iwrite=iwrite+1
             end do
             izero=0
             tfzero=.false.
             tzero=.false.
          end if
          string(iwrite:iwrite)='8'
          iwrite=iwrite+1
          tzero=.false.
       case (9)
          if(tfzero) then
             !write the foregoing zeros
             do k=1,izero
                string(iwrite:iwrite)='0'
                iwrite=iwrite+1
             end do
             izero=0
             tfzero=.false.
             tzero=.false.
          end if
          string(iwrite:iwrite)='9'
          iwrite=iwrite+1
          tzero=.false.
       end select
    end do
    

  end function r_2_s

  
end module libxml2f90_strings_module


module ll_module
  type ll_type
     character(32)            :: LL_ID
     character(32)            :: TAG
     character(32)            :: ID
     character(1),dimension(:),pointer :: VALUE
!     real(8),dimension(:),pointer     :: row
     type(ll_type),pointer    :: NEXT_LL
     type(ll_type),pointer    :: FIRST_TAG
     type(ll_type),pointer    :: NEXT_TAG
     type(ll_type),pointer    :: UP_TAG
     type(ll_type),pointer    :: DOWN_TAG
     type(ll_type),pointer    :: FIRST_ID
     type(ll_type),pointer    :: NEXT_ID   
     type(ll_type),pointer    :: NEXT_PID   
  end type ll_type

  type(ll_type),pointer       :: LL_ROOT
  type(ll_type),pointer       :: THIS,THISTEMP,THISTMP1

  logical(4)                  :: INITIALIZED=.false.
  logical(4)                  :: TCASESENSITIVE=.true.


  interface operator(.xmleq.)
     module procedure xmlequals
  end interface

contains
  function xmlequals(string1,string2) result(teq)
    use libxml2f90_strings_module
    character(*),intent(in)    :: string1
    character(*),intent(in)    :: string2
    logical(4)                 :: teq
    
    if(tcasesensitive) then
       teq=(trim(adjustl(string1)).eq.trim(adjustl(string2)))
    else
       !make the strings uppercase for comparison
       teq=(+trim(adjustl(string1)).eq.+trim(adjustl(string2)))
    end if
  end function xmlequals


end module ll_module


module libxml2f90_interface_module
  interface libxml2f90__addid
     module procedure libxml2f90_ll_addidr8
     module procedure libxml2f90_ll_addidr8a
     module procedure libxml2f90_ll_addidc8
  end interface

  
  interface libxml2f90__getid
     subroutine libxml2f90__ll_getr8(id,size_,val)
       implicit none
       character(*),intent(in)         :: id
       integer(4),intent(in)           :: size_
       real(8),intent(out)             :: val(size_)
     end subroutine libxml2f90__ll_getr8
     subroutine libxml2f90__ll_getr8_(id,val)
       implicit none
       character(*),intent(in)         :: id
       real(8),intent(out)             :: val
     end subroutine libxml2f90__ll_getr8_
     subroutine libxml2f90__ll_getc8(id,size_,val)
       implicit none
       character(*),intent(in)         :: id
       integer(4),intent(in)           :: size_
       complex(8),intent(out)          :: val(size_)
     end subroutine libxml2f90__ll_getc8
     subroutine libxml2f90__ll_getc8_(id,val)
       implicit none
       character(*),intent(in)         :: id
       complex(8),intent(out)          :: val
     end subroutine libxml2f90__ll_getc8_
     subroutine libxml2f90__ll_geti4(id,size_,val)
       implicit none
       character(*),intent(in)         :: id
       integer(4),intent(in)           :: size_
       integer(4),intent(out)          :: val(size_)
     end subroutine libxml2f90__ll_geti4
     subroutine libxml2f90__ll_geti4_(id,val)
       implicit none
       character(*),intent(in)         :: id
       integer(4),intent(out)          :: val
     end subroutine libxml2f90__ll_geti4_
     subroutine libxml2f90__ll_getl4(id,size_,val)
       implicit none
       character(*),intent(in)         :: id
       integer(4),intent(in)           :: size_
       logical(4),intent(out)          :: val(size_)
     end subroutine libxml2f90__ll_getl4
     subroutine libxml2f90__ll_getl4_(id,val)
       implicit none
       character(*),intent(in)         :: id
       logical(4),intent(out)          :: val
     end subroutine libxml2f90__ll_getl4_
     subroutine libxml2f90__ll_getstring(id,size_,val)
       implicit none
       character(*),intent(in)         :: id
       integer(4),intent(in)           :: size_
       character(*),intent(out)        :: val(size_)
     end subroutine libxml2f90__ll_getstring
     subroutine libxml2f90__ll_getstring_(id,val)
       implicit none
       character(*),intent(in)         :: id
       character(*),intent(out)        :: val
     end subroutine libxml2f90__ll_getstring_
  end interface

  interface libxml2f90__getpid
     subroutine libxml2f90__ll_getpr8(id,size_,val)
       implicit none
       character(*),intent(in)         :: id
       integer(4),intent(in)           :: size_
       real(8),intent(out)             :: val(size_)
     end subroutine libxml2f90__ll_getpr8
     subroutine libxml2f90__ll_getpr8_(id,val)
       implicit none
       character(*),intent(in)         :: id
       real(8),intent(out)             :: val
     end subroutine libxml2f90__ll_getpr8_
     subroutine libxml2f90__ll_getpc8(id,size_,val)
       implicit none
       character(*),intent(in)         :: id
       integer(4),intent(in)           :: size_
       complex(8),intent(out)          :: val(size_)
     end subroutine libxml2f90__ll_getpc8
     subroutine libxml2f90__ll_getpc8_(id,val)
       implicit none
       character(*),intent(in)         :: id
       complex(8),intent(out)          :: val
     end subroutine libxml2f90__ll_getpc8_
     subroutine libxml2f90__ll_getpi4(id,size_,val)
       implicit none
       character(*),intent(in)         :: id
       integer(4),intent(in)           :: size_
       integer(4),intent(out)          :: val(size_)
     end subroutine libxml2f90__ll_getpi4
     subroutine libxml2f90__ll_getpi4_(id,val)
       implicit none
       character(*),intent(in)         :: id
       integer(4),intent(out)          :: val
     end subroutine libxml2f90__ll_getpi4_
     subroutine libxml2f90__ll_getpl4(id,size_,val)
       implicit none
       character(*),intent(in)         :: id
       integer(4),intent(in)           :: size_
       logical(4),intent(out)          :: val(size_)
     end subroutine libxml2f90__ll_getpl4
     subroutine libxml2f90__ll_getpl4_(id,val)
       implicit none
       character(*),intent(in)         :: id
       logical(4),intent(out)          :: val
     end subroutine libxml2f90__ll_getpl4_
     subroutine libxml2f90__ll_getpstring(id,size_,val)
       implicit none
       character(*),intent(in)         :: id
       integer(4),intent(in)           :: size_
       character(*),intent(out)        :: val(size_)
     end subroutine libxml2f90__ll_getpstring
     subroutine libxml2f90__ll_getpstring_(id,val)
       implicit none
       character(*),intent(in)         :: id
       character(*),intent(out)        :: val
     end subroutine libxml2f90__ll_getpstring_
  end interface

contains
    
  subroutine libxml2f90_ll_addidr8(id,value)
    use libxml2f90_strings_module
    implicit none
    character(*),intent(in)        :: id
    real(8),intent(in)             :: value
    integer(4)                     :: size_ 
    character(32)                  :: ch
    character, dimension(32)       :: chvec
    integer                        :: s
    !............................................
    
    !convert the value to a string
    ch=r_2_s(value)
    
    size_=len(trim(adjustl(ch)))
    do s=1,size_
      chvec(s)=ch(s:s)
    end do
    call libxml2f90_ll_addid(id,size_,chvec)
    
    return
  end subroutine libxml2f90_ll_addidr8
  
  
  subroutine libxml2f90_ll_addidr8a(id,n,value)
    use libxml2f90_strings_module
    implicit none
    character(*),intent(in)        :: id
    integer(4),intent(in)          :: n
    real(8),intent(in)             :: value(n)
    character(32)                  :: ch(n)
    integer(4)                     :: i,l,ipos
    integer(4)                     :: size_
    character(1),allocatable       :: ch1(:)
    !............................................
    ipos=0
    l=0
    !convert
    do i=1,n
       ch(i)=r_2_s(value(i))
       l=l+len(trim(adjustl(ch(i))))
    end do
    
    allocate(ch1(l+n))!n: a blank between the numbers
    
    do i=1,n
       do l=1,len(trim(adjustl(ch(i))))
          ipos=ipos+1
          ch1(ipos)=ch(i)(l:l)
       end do
       ipos=ipos+1 !a blank between the numbers
       ch1(ipos)=' '
    end do
    
    call libxml2f90_ll_addid(id,l+n,ch1)
    deallocate(ch1)
    return
  end subroutine libxml2f90_ll_addidr8a


  subroutine libxml2f90_ll_addidc8(id,value)
    use libxml2f90_strings_module
    implicit none
    character(*),intent(in)        :: id
    complex(8),intent(in)          :: value
    integer(4)                     :: size_ 
    character(32)                  :: ch
    !............................................
    
!!__    !convert the value to a string
!!__    ch=r_2_s(value)
!!__    
!!__    size_=len(trim(adjustl(ch)))
!!__    call libxml2f90_ll_addid(id,size_,ch)
!!__    
!!__    return
  end subroutine libxml2f90_ll_addidc8
 


end module libxml2f90_interface_module



subroutine how_to_read_a_file()
  implicit none
  integer(4)                      :: i
!  logical(4)                      :: tread

  call libxml2f90__readin_file('config.cntl','CNTL')

  call libxml2f90__ll_report('CNTL',6,.true.)

  call libxml2f90__ll_selectlist('CNTL')
  call libxml2f90__ll_exist('DOWN','tag1',i)
  call libxml2f90__ll_selecttag('DOWN','tag1',2)

end subroutine how_to_read_a_file


subroutine test()
!this is for testing purpose
implicit none
call how_to_read_a_file()
end subroutine test



subroutine libxml2f90__get_fileunit(file,nfil)
  !=======================================================
  !===== RETURNS THE FILEUNIT FOR FILE IF OPENED,
  !===== OTHERWISE 0
  !=======================================================
  use libxml2f90_module
  IMPLICIT NONE
  integer(4),intent(out)           :: nfil
  character(*),intent(in)          :: file
  character(256)                   :: afile
  logical(4)                       :: topened
  integer(4)                       :: ifil
  integer(4)                       :: rb,lb

  nfil=0
  do ifil=nfilstart,nfilmax
     inquire(UNIT=ifil,OPENED=topened,NAME=afile)

     !some compilers return the whole path!
     !get the right boundary
     rb=scan(afile,' ')
     if(rb.eq.0) then
        rb=len(afile)
     end if

     lb=rb-len(file)

     if(topened.and.file.eq.afile(lb:rb)) then
        nfil=ifil
        exit
     end if
  end do
  return
end subroutine libxml2f90__get_fileunit


subroutine libxml2f90__getunit(nfil)
  !=======================================================
  !===== RETURNS THE NEXT FREE FILEUNIT
  !=======================================================
  use libxml2f90_module
  IMPLICIT NONE
  integer(4),intent(out)           :: nfil
  logical(4)                       :: topened
  
  do nfil=nfilstart,nfilmax
     inquire(UNIT=nfil,OPENED=topened)
     if(.not.topened)exit
  end do
  if(nfil.eq.nfilmax) then
     print*,"ERROR STOP IN LIBXML2F90 - GETUNIT: REACHED MAXIMUM AVAILABLE # OF FILEUNITS"
     stop  
  end if
  
  if(nfil.gt.nfilmaxused)nfilmaxused=nfil

  return
end subroutine libxml2f90__getunit

     subroutine libxml2f90__openfile(id,file,nfil)
       use libxml2f90_module
       IMPLICIT NONE
       character(*),intent(in)           :: id! 
       character(*),intent(in)           :: file 
       
       integer(4),intent(out)              :: nfil ! file id
       logical(4)                          :: EX !does the file exist
       !..........................................................

       if(id.eq.'READ') then
          !test if inputfile exists
          INQUIRE (FILE=trim(adjustl(file)), EXIST = EX)
          IF ( EX ) THEN
             !assign the file unit
             call libxml2f90__getunit(nfil)

          !open inputfile
             OPEN (nfil, FILE=trim(adjustl(file)), STATUS="OLD", &
                  &ACCESS="SEQUENTIAL", ACTION="READ", POSITION="REWIND")
             !
             !PRINT *,trim(adjustl(file))," OPENED"
             !
             
          else
             print*,"ERROR STOP: ",trim(adjustl(file))," NOT FOUND"
             stop
          end IF
       else if(id.eq.'WRITE') then
          !assign the file unit
          call libxml2f90__getunit(nfil)
          
          OPEN (nfil, FILE=trim(adjustl(file)), STATUS="UNKNOWN", &
               &ACCESS="SEQUENTIAL", ACTION="WRITE", POSITION="REWIND")

       else if(id.eq.'APPEND') then
          !assign the file unit
          call libxml2f90__getunit(nfil)
          
          OPEN (nfil, FILE=trim(adjustl(file)), STATUS="UNKNOWN", &
               &ACCESS="SEQUENTIAL", ACTION="WRITE", POSITION="APPEND")         

       else
          print*, 'ID not recognized in libxml2f90__openfile: ',trim(adjustl(id))
          stop
       end if
       return
     end subroutine libxml2f90__openfile
        
     subroutine libxml2f90__closeall()
       use libxml2f90_module
       IMPLICIT NONE
       integer(4)       :: i
       logical(4)       :: topened
       !.................................

       do i=nfilstart,nfilmaxused
          inquire(UNIT=i,OPENED=topened)
          if(.not.topened)close(i)
       end do
     end subroutine libxml2f90__closeall

     subroutine libxml2f90__closefile(file)
       use libxml2f90_module
       IMPLICIT NONE
       character(*),intent(in)       :: file
       integer(4)                    :: nfil
       !.................................

       call libxml2f90__get_fileunit(trim(adjustl(file)),nfil)
       if(nfil.eq.0) then
          print*,'ERROR STOP IN libxml2f90__closefile: the file ',trim(adjustl(file))
          print*,'CAN NOT BE CLOSED'
          stop
       end if
       close(nfil)

       return
     end subroutine libxml2f90__closefile


     subroutine libxml2f90__flush(nfil)
       use libxml2f90_module
       IMPLICIT NONE
       integer(4),intent(in)       :: nfil
       character(256)              :: file,access,form,action,pad,direct
       character(256)              :: status,blank,position,delim
       logical(4)                  :: opened
       integer(4)                  :: recl
       !....................................
       inquire(UNIT=nfil,NAME=file,OPENED=opened,ACCESS=access,FORM=form,&
            &RECL=recl,BLANK=blank,POSITION=position,ACTION=action,&
            &DELIM=delim,PAD=pad)

       if(.not.opened) then
          print*,'CAN NOT FLUSH A NON OPENED FILE'
       else
          close(nfil)
          if(trim(adjustl(position)).eq.'UNDEFINED')position='APPEND'
          open(UNIT=nfil,FILE=trim(adjustl(file)),STATUS='OLD',&
               &ACCESS=trim(adjustl(access)),FORM=trim(adjustl(form)),&
               &RECL=recl,BLANK=trim(adjustl(blank)),&
               &POSITION=trim(adjustl(position)),ACTION=trim(adjustl(action)),&
               &DELIM=trim(adjustl(delim)),PAD=trim(adjustl(pad)))
!!__          print*,'flushed ',trim(adjustl(file))
!!__          print*,opened
!!__          print*,trim(adjustl(access))
!!__          print*,trim(adjustl(form))
!!__          print*,recl
!!__          print*,trim(adjustl(blank))
!!__          print*,trim(adjustl(position))
!!__          print*,trim(adjustl(action))
!!__          print*,trim(adjustl(delim))
!!__          print*,trim(adjustl(pad))
       end if
     end subroutine libxml2f90__flush



subroutine libxml2f90__set_default_ll_id(llid)
use libxml2f90_module
implicit none
character(*),intent(in)             :: llid
!..............................
default_llid=llid
end subroutine libxml2f90__set_default_ll_id




!=========================================================
subroutine libxml2f90__readin_file(file,ll_id)
  !=======================================================
  !===== READS THE FILE INTO 
  !===== READINSTRING,a character(1) array
  !===== INDEPENDENT FROM THE ACTUAL LINELENGTH!
  !===== CUTS TABS AND UNNECESSARY BLANKS ETC.
  !===== SETS THE LL_ID AND CALLS THE PARSER
  !=======================================================
  implicit none
  character(*),intent(in)       :: file
  character(*),intent(in)       :: ll_id
  integer(4)                    :: nfil  !config-file unit
  !.................................
  
  !open file
  call libxml2f90__openfile('READ',file,nfil)
  !read file
  call libxml2f90_readin_file(nfil,ll_id)
  close(nfil)
end subroutine libxml2f90__readin_file



!=========================================================
subroutine libxml2f90__readin_nfil(nfil,ll_id)
  implicit none
  integer(4),intent(in)        :: nfil
  character(*),intent(in)       :: ll_id
  !.................................
  
  call libxml2f90_readin_file(nfil,ll_id)

end subroutine libxml2f90__readin_nfil



!=========================================================
subroutine libxml2f90_readin_file(nfil,ll_id)
  !=======================================================
  !===== READS THE FILE INTO 
  !===== READINSTRING,a character(1) array
  !===== INDEPENDENT FROM THE ACTUAL LINELENGTH!
  !===== CUTS TABS AND UNNECESSARY BLANKS ETC.
  !===== SETS THE LL_ID AND CALLS THE PARSER
  !=======================================================
  use libxml2f90_module
  use libxml2f90_strings_module
  implicit none
  integer,intent(in)            :: nfil
  character(*),intent(in)       :: ll_id
  integer(4)                    :: eos,is,nwrite,i
  character(256)                :: string !read in 256 chunks
  logical(4)                    :: tblank,tnewline
  integer(4),allocatable        :: tempia(:)
  !.................................
  
  is = 0

  !===========read file into stringarry
  filelines=0
  if(allocated(readstring)) deallocate(readstring)
  allocate(readstring(arraystep))
  deallocate(readstring)
  allocate(readstring(arraystep))
  if(allocated(lineposa)) deallocate(lineposa)
  allocate(lineposa(arraystep))
  lineposa(:)=0

  nwrite=1
  tblank=.false.

  do !over the lines
     tnewline=.true.
     do !the chunks per line
        READ (nfil, FMT="(A)" , ADVANCE="NO" ,SIZE=is,IOSTAT= eos) string !read a chunk
        !=======================================
        !process that chunk (is=0 @ end of file!)
        !=======================================

        if(eos.eq.-1) then
           !the end of file leads to double reading of the last line
           !-->if eos=-1 we have end of file (-2 end of line)
           exit !the readin do-loop
        end if

        
        !test if the line starts with '<#' -> we use the line as comment
        ! we also skip lines containing <?xml and <!DOCTYPE
        if(tnewline) then
           tnewline=.false.
           !find first nonblank 
           i=verify(string,' ')
           if(i.gt.0) then
              if(string(i:i).eq.'<'.and.(string(i+1:i+1).eq.'#'&
                   &.or.+string(i+1:i+4).eq.'?XML'&
                   &.or.+string(i+1:i+8).eq.'!DOCTYPE')) then
                 do
                    if(eos.eq.-1.or.eos.eq.-2) exit !end of line (-2) or file (-1) reached 
                    READ (nfil, FMT="(A)" , ADVANCE="NO" ,SIZE=is,IOSTAT= eos) string !read a chunk
                 end do
                 tnewline=.true. !we skip the whole line --> the next read uses a new line
                 !remember the lineposition (for the comment lines also!)
                 !----------------
                 filelines=filelines+1
                 !keep the information about how many character(1) 
                 !entries come from which line. this allows us to 
                 !inform the user in which line syntax errors occur
                 
                 !check the size of the array
                 if(size(lineposa).lt.filelines) then
                    allocate(tempia(filelines-1))
                    tempia(:)=lineposa(:)
                    deallocate(lineposa)
                    allocate(lineposa(filelines+arraystep))
                    lineposa(1:filelines-1)=tempia(1:filelines-1)
                    deallocate(tempia)
                 end if
                 lineposa(filelines)=nwrite
                 !--------------

                 cycle
              end if
           end if
        end if


        !ensure that readstring is large enough
        if(nwrite+is.gt.ubound(readstring,1)) then
            allocate(tempstringa(ubound(readstring,1)))
           tempstringa(:)=readstring(:)
           deallocate(readstring)
           allocate(readstring(ubound(tempstringa,1)+arraystep))
           readstring(1:ubound(tempstringa,1))=tempstringa(1:ubound(tempstringa,1))
           deallocate(tempstringa)
        end if

        !now transfer the chunk=> readstring and cut unneeded characters
        do i=1,is !if is=0 we skip
           if(iachar(string(i:i)).lt.33.or.&
                &(trmquotes.and.iachar(string(i:i)).eq.34).or.&
                &(trmquotes.and.iachar(string(i:i)).eq.39).or.&
                &(trmcomma.and.iachar(string(i:i)).eq.44)) then 
              !special characters including blank or "/' if should be removed or 
              if(xmlformat.eq.1.or.xmlformat.eq.3) then !blank removal
                 if(.not.tblank) then
                    !write a blank instead if we do not have a blank at the last position 
                    readstring(nwrite)=' '
                    nwrite=nwrite+1
                    tblank=.true.
                 end if !else we drop the character
              else !no blank removal,just replace special characters
                 !write a blank instead
                 readstring(nwrite)=' '
                 nwrite=nwrite+1
              end if
           else if(string(i:i).eq.'=') then != corresponds to a blank on the left side
              !we do not need blanks left of a = if we use blankremoval:
              if(readstring(nwrite-1).eq.' '&
                   &.and.(xmlformat.eq.1.or.xmlformat.eq.3)) then 
                 !we can remove a blank to the left of a =
                 readstring(nwrite-1)=string(i:i)
                 !nwrite has to stay the same!
              else
                 !no blanks on the left side or no blank removal
                 !-> just copy the =
                 readstring(nwrite)=string(i:i)
                 tblank=.true.
                 nwrite=nwrite+1
              end if
           else if(string(i:i).eq.'>') then !> corresponds to a blank on the left side
              readstring(nwrite)=string(i:i)
              tblank=.true.
              nwrite=nwrite+1
           else if(nwrite.gt.1) then
              if(string(i:i).eq.'<'.and.readstring(nwrite-1).eq.' '&
                   &.and.(xmlformat.eq.1.or.xmlformat.eq.3)) then 
                 !we can remove a blank to the left of a <
                 readstring(nwrite-1)=string(i:i)
                 !nwrite has to stay the same!
              else    !no blank, no special character, no '='
                 readstring(nwrite)=string(i:i)
                 tblank=.false.
                 nwrite=nwrite+1
              end if
           else    !no blank, no special character, no '='
              readstring(nwrite)=string(i:i)
              tblank=.false.
              nwrite=nwrite+1
           end if
        end do
        if(eos.eq.-1) exit !end of file (-1) reached
        if(eos.eq.-2) then !end of line (-2) or file (-1) reached
           !add a blank for the linebreak if there is no blank at the end
           if(readstring(nwrite-1).ne.' ') then
              readstring(nwrite)=' '
              tblank=.true.
              nwrite=nwrite+1
           end if
           exit
        end if
     end do
     filelines=filelines+1
     !keep the information about how many character(1) 
     !entries come from which line. this allows us to 
     !inform the user in which line syntax errors occur
     
     !check the size of the array
     if(size(lineposa).lt.filelines) then
        allocate(tempia(filelines-1))
        tempia(:)=lineposa(:)
        deallocate(lineposa)
        allocate(lineposa(filelines+arraystep))
        lineposa(1:filelines-1)=tempia(1:filelines-1)
        deallocate(tempia)
     end if

     lineposa(filelines)=nwrite
     
     if(eos.eq.-1) exit 
  end do

  !cut the unneeded rest
  allocate(tempstringa(nwrite-2))!-2 because we have a blank at the end
!  tempstringa(:)=readstring(:)
  tempstringa(1:(nwrite-2))=readstring(1:(nwrite-2))
  deallocate(readstring)
  allocate(readstring(ubound(tempstringa,1)))
  readstring(:)=tempstringa(:)
  deallocate(tempstringa)
  

  if(ttransform_paw) call libxml2f90_transform_paw()
  call libxml2f90__set_default_ll_id(ll_id)
  call libxml2f90_parse_file()

end subroutine libxml2f90_readin_file


!=========================================================
subroutine libxml2f90__findinchara(dim,chara,start,toright,nsigns,signs,pos,sign)
!finds the character of chara (array) to the right or left, returns the first occurence
!search start 1 sign right of start!
  implicit none
  integer(4),intent(in)         :: dim
  character(1),intent(in)       :: chara(dim)
  integer(4),intent(in)         :: start !start search here
  logical(4),intent(in)         :: toright !.true.: search next, .false. search to left
  integer(4),intent(in)         :: nsigns !the number of signs to search for
  character(1),intent(in)       :: signs(nsigns) !the signs in arrayformat
  integer(4),intent(out)        :: pos !the position we found one of the signs
  integer(4),intent(out)        :: sign !the found sign. =0 if we didn't find any
  !note: start=0 to start from the first char to right
  !      start=dim+1 to start from the last char to left

  integer(4)                    :: i,j
  !..........................

  pos=0
  sign=0
  if(toright) then
     do i=start+1, dim
        do j=1,nsigns
           if(chara(i).eq.signs(j)) then
              !we found a sign
              sign=j
              pos=i
              exit
           end if
        end do
        if(sign.ne.0) exit !we already found a sign
     end do
 
  else
     !to left
     do i=1, start-1
        do j=1,nsigns
           if(chara(start-i).eq.signs(j)) then
              !we found a sign
              sign=j
              pos=start-i
              exit
           end if
        end do
        if(sign.ne.0) exit !we already found a sign
     end do
 
     
  end if

  return
end subroutine libxml2f90__findinchara



!=================================================
subroutine libxml2f90__setformat(xmlformat_)
  use libxml2f90_module
  implicit none
  integer(4),intent(in)        :: xmlformat_
  !.........................................
  if(xmlformat.le.3) then
     xmlformat=xmlformat_
  else
     print*,'ERROR STOP IN LINKLIST: SETFORMAT: FORMAT NOT SUPPORTED: ',xmlformat_
  end if
  return
end subroutine libxml2f90__setformat



subroutine libxml2f90_parse_file()
  !=======================================================
  !===== PARSES THE readstring INTO A LINKLIST REFERENCED
  !===== BY DEFAULT_LLID
  !===== depends on xmlformat 
  !=======================================================
  use libxml2f90_module
  use ll_module
  implicit none
  integer(4)                      :: e !# of elements in readstring
  integer(4)                      :: p !the actual position
  integer(4)                      :: pcl!the last >
  integer(4)                      :: fp !the position we found a sign
  integer(4)                      :: isign !the of the sign found
  integer(4)                      :: llb,lb,lbv,rb,rrb,plb !some bounds
  logical(4)                      :: tflex=.false.
  character(32)                   :: dummy32
  character(32)                   :: sh32   !the tag for short hand notation
  logical(4)                      :: tsh    !do we have a short hand notation?
  character(1),allocatable        :: val(:)
  integer(4)                      :: j,line
  character,dimension(2)          :: ctmp_le,ctmp_ge,ctmp_gt_sp
  character,dimension(1)          :: ctmp_gt,ctmp_eq,ctmp_sp
 !..................................
  tsh=.false.
  tflex=.false.


  !init the linklist (don't forget to set the default_llid first)
  call libxml2f90_ll_add_list(default_llid)
  pcl=0
  e=ubound(readstring,1) !the # of elements
  p=0 !we start the search from here to right

  ctmp_le = (/'<','='/)
  ctmp_ge = (/'>','='/)
  ctmp_gt_sp = (/'>',' '/)
  ctmp_gt = (/'>'/)
  ctmp_eq = (/'='/)
  ctmp_sp = (/' '/)
  do
     call libxml2f90__findinchara(e,readstring,p,.true.,2,ctmp_le,fp,isign)
     if(isign.eq.0) exit !we reached the end


!print*,tflex,'tflex',':',readstring(p+1:fp),':'
     if(isign.eq.1) then
        !<
        !Did we have intermediate text without = ? keep it!
         !deleted the +1 for both pcl at 10.10.2005 because the diference was length-1!!!
         if(.not.tflex.and.(((fp-1)-(pcl)).gt.1.or.(((fp-1)-(pcl)).eq.1&
             &.and.(readstring(pcl+1).ne.' ')))) then 
           if(allocated(val)) deallocate(val)
           allocate(val((fp)-(pcl+1)))
           do j=1,(fp)-(pcl+1)          !rb-1-(lb+1)+1
            val(j)=readstring(pcl+j)
           end do
           lbact=fp !for the errorinfo
           call libxml2f90_ll_addpureid(size(val),val)
        end if
        
        !is it a comment?
        if(readstring(fp+1).eq.'!'.and.readstring(fp+2).eq.'-'.and.readstring(fp+3).eq.'-') then
           lb=fp
           !remember the actual position for line information on error
           lbact=lb+2

           !we have a comment <!--
           !search for the closing -->
           call find_closing_comment(e,readstring,lb,lbact,p,pcl,tflex) !enters recursion if nec.
           cycle
        else !it is no comment: process as usual
           lb=fp
           lbact=lb+2
           !find the closing >
           call libxml2f90__findinchara(e,readstring,lb,.true.,1,ctmp_gt,fp,isign)
           if(isign.eq.0) then
              call libxml2f90_getline(lbact,line)
              print*,'ERROR STOP IN PARSING: CANNOT FIND A CLOSING > (line ',line,')!'
              stop
           else if(isign.eq.1) then !found the '>'
              
              if(readstring(fp-1).eq.'/') then !we have a short hand notation <tag />
                 tsh=.true.
                 rb=fp-1
              else
                 rb=fp
              end if
              
              pcl=fp
              tflex=.false.
              p=rb  !for the next loop;remember, we do one step at the beginning of the search
              if(readstring(lb+1).eq.'/') then 
                 !we have a closing tag
                 !transfer it into a character(32) variable
                 dummy32=''
                 do j=1,rb-1-(lb+2)+1
                    dummy32(j:j)=readstring(lb+1+j)
                 end do
                 lbact=lb !for the line information
                 call libxml2f90_ll_closetag(trim(adjustl(dummy32)))
              else  
                 !we have an opening tag that can contain pid's!
                 
                 ! do we have pid's inside?
                 call libxml2f90__findinchara(rb,readstring(1:rb),lb,.true.,1,ctmp_eq,fp,isign)
                 if(isign.eq.0) then !no pid's inside the opening tag
                    dummy32=''
                    do j=1,rb-1-(lb+1)+1
                       dummy32(j:j)=readstring(lb+j)
                    end do

                    if(tsh) sh32=dummy32 !keep the tag for sh notation
                    call libxml2f90_ll_opentag(trim(adjustl(dummy32)))
                 else if(isign.eq.1) then !pid's inside 
                    !=== seperate the tag from the pid's ===
                    call libxml2f90__findinchara(e,readstring,lb,.true.,1,ctmp_sp,fp,isign)
                    if(isign.eq.0) then
                       call libxml2f90_getline(lb+2,line)
                       print*,'ERROR STOP IN LINKLIST: CANNOT FIND RIGHTBOUND FOR TAG&
                            & CONTAINING PIDs! THE WHOLE THING IS: ',readstring(lb:fp)
                       print*,'(line ',line,')!'
                       stop
                    else if(isign.eq.1) then
                       dummy32=''
                       do j=1,fp-1-(lb+1)+1
                          dummy32(j:j)=readstring(lb+j)
                       end do
                       if(tsh) sh32=dummy32 !keep the tag for sh notation
                       call libxml2f90_ll_opentag(trim(adjustl(dummy32)))
                       lb=fp-1 !think about this!!! 
                    end if
                    
                    !=== process the pid's ====
                    !enter a subloop to process the pid's
                    plb=lb
                    do
                       call libxml2f90__findinchara(rb,readstring(1:rb),plb,.true.,1,ctmp_sp,fp,isign)
                       if(isign.eq.0) exit !no pid's left,exit the loop
                       !we have found another pid 
                       !=============================
                       !=
                       lb=fp
!++++++++++++++
!                       !do we have a ' ' left of the = ? that's not allowed (parsing is not possible)
!                       if(readstring(lb-1).eq.' ') then
!                          print*,"ERROR STOP IN PARSING: IT'S NOT ALLOWED TO HAVE A ' ' NEXT (LEFT SIDE) TO AN ="
!                          stop
!                       end if
!++++++++++++++
                       !find next sign to the right for the valuestring
                       call libxml2f90__findinchara(e,readstring,lb,.true.,2,ctmp_ge,fp,isign)
                       if(isign.eq.0) then
                          print*,'ERROR STOP IN PARSING: CANNOT FIND A RIGHT BOUND FOR A ='
                          call libxml2f90_getline(lb+2,line)
                          print*,'(line ',line,')!'
                          stop
                       else if(isign.eq.1) then !found a >
                          if(readstring(fp-1).eq.'/') then !we have a short hand notation <tag pid=val/>
                             rrb=fp-1
                             plb=rrb-1
                          else
                             rrb=fp
                             plb=rrb-1
                          end if
                       else if(isign.eq.2) then
                          != -> go to the left to find the correct boundary
                          rrb=fp
                          plb=rrb-1!check this!
                          call libxml2f90__findinchara(e,readstring,rrb,.false.,1,ctmp_sp,fp,isign)
                          if(isign.eq.0) then
                             print*,'ERROR STOP IN PARSING: CANNOT FIND DELIMITER BETWEEN TWO ='
                             call libxml2f90_getline(lb+2,line)
                             print*,'(line ',line,')!'
                             stop
                          else if(isign.eq.1) then
                             !found the delemiter
                             plb=rrb-1 !this sould be exactly 1 before the next=
                             rrb=fp
                          end if
                       end if
                       !now we have the value between lb+1 and rrb-1
                       if(lb.gt.rrb) then
                          print*,'ERROR STOP IN PARSING: THE RIGHT DELIMITER FOR A VALUE IS LEFT OF THE LEFT DELIMITER!!!'
                          call libxml2f90_getline(lb+2,line)
                          print*,'(line ',line,')!'
                          stop
                       end if
                       !keep the lb for the value in lbv because lb can move to the left if 
                       !we have a ' ' left of the =
                       lbv=lb

                       !+++ let's find the id +++
!++++++++++++++++
                       !cut the ' ' left of the '='               
                       do
                          if(readstring(lb-1).eq.' ') then !lb-1 because the search starts one position to the left
                             lb=lb-1
                          else
                             exit
                          end if
                       end do
!++++++++++++++++

                       call libxml2f90__findinchara(e,readstring,lb,.false.,2,ctmp_gt_sp,fp,isign)
                       llb=fp
                       if(isign.eq.0) then
                          print*,'ERROR STOP IN PARSING: CANNOT FIND A LEFT BOUND FOR A ='
                          call libxml2f90_getline(lb+2,line)
                          print*,'(line ',line,')!'
                          stop
                       else if(isign.eq.1) then
                          print*,'ERROR STOP IN LINKLIST: > IS NOT ALLOWED AS LEFT DELIMITER&
                               & FOR A PID!'
                          call libxml2f90_getline(fp,line)
                          print*,'(line ',line,')!'

                       else if(isign.eq.2) then
                          llb=fp
                       end if
                       !we now have the id between llb+1 and lb-1
                       if(llb.gt.lb) then
                          print*,'ERROR STOP IN PARSING: THE RIGHT DELIMITER FOR AN ID IS LEFT OF THE LEFT DELIMITER!!!'
                          call libxml2f90_getline(lb+2,line)
                          print*,'(line ',line,')!'
                          stop
                       end if
                       
                       !write the value/id to the linklist
                       !transfer it into a character(32)/ character(256) variable
                       dummy32=''
                       do j=1,lb-1-(llb+1)+1          !rb-1-(lb+1)+1
                          dummy32(j:j)=readstring(llb+j)
                       end do
                       
                       if(allocated(val)) deallocate(val)
                       allocate(val(rrb-1-(lbv+1)+1))
                       do j=1,rrb-1-(lbv+1)+1          !rb-1-(lb+1)+1
                          val(j)=readstring(lbv+j)
                       end do
                       lbact=lb !for the errorinformation
                       call libxml2f90_ll_addpid(trim(adjustl(dummy32)),size(val),val)
                    end do!from process pid's
                 end if !from pid's inside
                 !close the shorthand notation tag <tag pid=val /> if necessary
                 if(tsh) then 
                    call libxml2f90_ll_closetag(trim(adjustl(sh32)))
                    tsh=.false.
                 end if
              end if !from we have an opening tag
           end if ! we found a > as right delim. of a <
        end if !from it is no comment
     else if(isign.eq.2) then

       != between to tags (no pid)
        !cycle if we read pure Xml (no id=value between tags)
        if(xmlformat.le.1) then !0 and 1 are pure xml (without and with blank removal)
           p=fp
           cycle
        end if
        tflex=.true.
        lb=fp
!++++++++++++
!        !do we have a ' ' left of the = ? that's not allowed (parsing is not possible)
!        if(readstring(lb-1).eq.' ') then
!           print*,"ERROR STOP IN PARSING: IT'S NOT ALLOWED TO HAVE A ' ' NEXT (LEFT SIDE) TO AN ="
!           stop
!        end if
!++++++++++++

        !find next sign to the right for the valuestring
        call libxml2f90__findinchara(e,readstring,lb,.true.,2,ctmp_le,fp,isign)
        if(isign.eq.0) then
           print*,'ERROR STOP IN PARSING: CANNOT FIND A RIGHT BOUND FOR A ='
           call libxml2f90_getline(lb+2,line)
           print*,'(line ',line,')!'
           stop
        else if(isign.eq.1) then
           rrb=fp
           p=rrb-1
        else if(isign.eq.2) then
           != -> go to the left to find the correct boundary
           rrb=fp
           p=rrb-1!check this!
           call libxml2f90__findinchara(e,readstring,rrb,.false.,1,ctmp_sp,fp,isign)
           if(isign.eq.0) then
              print*,'ERROR STOP IN PARSING: CANNOT FIND DELIMITER BETWEEN TWO ='
              call libxml2f90_getline(lb+2,line)
              print*,'(line ',line,')!'
              stop
           else if(isign.eq.1) then
              !found the delemiter
              p=rrb-1 !this sould be exactly 1 before the next=
              rrb=fp
           end if
        end if
        !now we have the value between lb+1 and rrb-1
        if(lb.gt.rrb) then
           print*,'ERROR STOP IN PARSING: THE RIGHT DELIMITER FOR A VALUE IS LEFT OF THE LEFT DELIMITER!!!'
              call libxml2f90_getline(lb+2,line)
              print*,'(line ',line,')!'
           stop
        end if
        !keep the lb for the value in lbv because lb can move to the left if 
        !we have a ' ' left of the =
        lbv=lb

        !+++ let's find the id +++
!++++++++++++++++
        !cut the ' ' left of the '='               
        do
           if(readstring(lb-1).eq.' ') then !lb-1 because the search starts one position to the left
              lb=lb-1
           else
              exit
           end if
        end do
!++++++++++++++++

        !let's find the id
        call libxml2f90__findinchara(e,readstring,lb,.false.,2,ctmp_gt_sp,fp,isign)
        llb=fp
        if(isign.eq.0) then
           print*,'ERROR STOP IN PARSING: CANNOT FIND A LEFT BOUND FOR A ='
           call libxml2f90_getline(lb+2,line)
           print*,'(line ',line,')!'
           stop
        else if(isign.eq.1) then
           !llb=fp already set above
        else if(isign.eq.2) then
           llb=fp
        end if
        !we now have the id between llb+1 and lb-1
        if(llb.gt.lb) then
           print*,'ERROR STOP IN PARSING: THE RIGHT DELIMITER FOR AN ID IS LEFT OF THE LEFT DELIMITER!!!'
           call libxml2f90_getline(lb+2,line)
           print*,'(line ',line,')!'
           stop
        end if
        
        !write the value/id to the linklist
        !transfer it into a character(32)/ character(256) variable
        dummy32=''
        do j=1,lb-1-(llb+1)+1          !rb-1-(lb+1)+1
           dummy32(j:j)=readstring(llb+j)
        end do

        if(allocated(val)) deallocate(val)
        allocate(val(rrb-1-(lbv+1)+1))
        do j=1,rrb-1-(lbv+1)+1          !rb-1-(lb+1)+1
           val(j)=readstring(lbv+j)
        end do
        call libxml2f90_ll_addid(trim(adjustl(dummy32)),size(val),val)
     end if !end from = sign
  end do
  
  if(allocated(val)) deallocate(val)
  return

end subroutine libxml2f90_parse_file


recursive subroutine find_closing_comment(e,readstring,lb,lbact,p,pcl,tflex)
  implicit none
  integer(4),intent(in)                  :: e 
  character(1),intent(in)                :: readstring(e)
  integer(4),intent(inout)               :: lb
  integer(4),intent(inout)               :: lbact
  integer(4),intent(inout)               :: p
  integer(4),intent(inout)               :: pcl
  logical(4),intent(inout)               :: tflex 
  integer(4)                             :: isign
  integer(4)                             :: lbact_used
  integer(4)                             :: fp
  integer(4)                             :: line
  character,dimension(2)                 :: ctmp_lt_gt

  ctmp_lt_gt = (/'>','<'/)

  do
     call libxml2f90__findinchara(e,readstring,lb,.true.,2,ctmp_lt_gt,fp,isign)
     lb=fp
     if(isign.eq.0) then
        call libxml2f90_getline(lbact,line)
        print*,'ERROR STOP IN PARSING: CANNOT FIND A CLOSING --> (line',line,')!'
        stop
     else if(isign.eq.1.and.readstring(fp-1).eq.'-'.and.readstring(fp-2).eq.'-') then    
        !found the closing --> with no intermediate <!--
        p=fp
        pcl=fp
        tflex=.false.
        exit ! this do loop
     else if(isign.eq.2.and.readstring(fp+1).eq.'!'.and.readstring(fp+2).eq.'-'&
          &.and.readstring(fp+3).eq.'-') then
        lb=fp
        lbact_used=lb+2
        !we have an intermediate <!--
        call find_closing_comment(e,readstring,lb,lbact_used,p,pcl,tflex) !recursion
     end if
  end do
  
  return
end subroutine find_closing_comment

!**********************************************************************
!**********************************************************************
!***********                 LINKLIST
!**********************************************************************
!**********************************************************************

subroutine libxml2f90__ll_selectlist(ll_id)
!===================================================================
!===== SELECTS A LINKLIST AND JUMPS ON THE FIRST TAGLAYER
!===== ERROR STOP IF IT IS EMPTY
!===================================================================
use ll_module
implicit none
character(*),intent(in)                :: ll_id
logical(4)                             :: tjump
!..............................................
if(.NOT.INITIALIZED) then
   stop 'ERROR STOP IN LINKLIST SELECT: NO LIST EXISTS!'
else
   this=>ll_root
   do
      if(trim(THIS%LL_ID).eq.trim(LL_ID)) then
         call libxml2f90__ll_down(tjump)
         if(.not.tjump)then
            print*,"ERROR STOP IN LINKLIST: THE LINKLIST IS EMPTY&
                 & AND CAN THEREFORE NOT BE SELECTED ",ll_id
            stop
         end if
         exit
      else
         if(associated(THIS%NEXT_LL)) then
         THIS=>THIS%NEXT_LL
         else
            print*,'ERROR STOP IN LINKLIST: THERE IS NO LIST CALLED: ',LL_ID
            stop
         end if
      end if
   end do
end if
!print*,'SELECTED LIST: ',THIS%LL_ID

end subroutine libxml2f90__ll_selectlist



!================================================================
subroutine libxml2f90__ll_exist(layer,tag,count)
!select the actual layer first!
use ll_module
implicit none
character(*),intent(in)         :: layer
character(*),intent(in)         :: tag
integer(4),intent(out)          :: count !how often is the tag present in the !next! layer
!........................................
count=0
thistemp=>this
if(layer.eq.'UP') then
   if(.not.associated(this%first_id%up_tag)) return 
       !note: this%uptag is not associated for 
       !the first layer beneath the linklistlayer
       !why do we not test up_tag%next_tag? read the comment
       !for down_tag
   this=>this%first_id%up_tag%first_tag%next_tag
else if(layer.eq.'DOWN') then
   if(.not.associated(this%first_id%down_tag)) return !just test
   !just test %down_tag here because %down_tag%next_tag 
   !exists if %down_tag exists, but test %down_tag%next_tag will 
   !give an addresserror if down_tag does not exist!
   this=>this%first_id%down_tag%next_tag
else if(layer.eq.'ACT') then
   if(.not.associated(this%first_id).and..not.associated(this%first_tag)) return !just test
   if(associated(this%first_id)) then
      if(.not.associated(this%first_id%first_tag)) return !just test
      if(.not.associated(this%first_id%first_tag%next_tag)) return !just test
      this=>this%first_id%first_tag%next_tag
   else !first_tag is associated
      if(.not.associated(this%first_tag%next_tag)) return !just test
      this=>this%first_tag%next_tag
   end if
else
   print*,'ERROR STOP IN LINKLIST: LL_EXIST ONLY ACCEPTS <UP>, <DOWN> OR <ACT>&
        &AND NOT: ',layer
   stop
end if

do
   if(this%tag.xmleq.tag) count=count+1
   if(associated(this%next_tag)) then
      this=>this%next_tag
   else
      exit
   end if
end do

this=>thistemp
return
end subroutine libxml2f90__ll_exist



!================================================================
subroutine libxml2f90__ll_selecttag(layer,tag,count)
use ll_module
implicit none
character(*),intent(in)         :: layer
character(*),intent(in)         :: tag
integer(4),intent(in)           :: count 
!logical(4),intent(out)          :: tselect=.false. !did we manage to select the tag
integer(4)                      :: count_
!........................................
count_=0
if(count.lt.1) then
   print*, 'ERROR IN LINKLIST SELECTTAG: COUNT MUST BE > 0!'
   stop
end if
if(layer.eq.'DOWN') then
   !go one layer down
   if(.not.associated(this%first_id%down_tag)) then !read the comment @ exist code
      print*, 'ERROR STOP IN LINKLIST SELECTTAG: NO LAYER BELOW THE ACTUAL!'
      stop
   else
      this=>this%first_id%down_tag%next_tag 
      do
         if(this%tag.xmleq.tag) then
            count_=count_+1
            if(count_.eq.count) exit
         end if
         if(associated(this%next_tag)) then
            this=>this%next_tag
         else
            print*,'ERROR INLINKLIST SELECTTAG: COULD NOT FIND THE ',count,'. tag: ',tag
            stop
         end if
      end do
   end if
else if(layer.eq.'ACT') then !read the comment @ exist code
   !the actual layer
   if(.not.associated(this%first_id%first_tag)) then
      print*, 'ERROR STOP IN LINKLIST SELECTTAG: NO TAG IN THE ACTUAL LAYER!'
      stop
   else
      this=>this%first_id%first_tag%next_tag 
      do
         if(this%tag.xmleq.tag) then
            count_=count_+1
            if(count_.eq.count) exit
         end if
         if(associated(this%next_tag)) then
            this=>this%next_tag
         else
            print*,'ERROR INLINKLIST SELECTTAG: COULD NOT FIND THE ',count,'. tag: ',tag,' IN THE ACTUAL LAYER'
            stop
         end if
      end do
   end if
else if(layer.eq.'UP') then
   !the layer above the actual
   if(.not.associated(this%first_id%up_tag).or.&
        &.not.associated(this%first_id%up_tag%first_tag)) then 
      !read the comment @ exist code
      print*, 'ERROR STOP IN LINKLIST SELECTTAG: NO TAG IN THE LAYER ABOVE THE ACTUAL!' 
      !this should not happen ;-)
      stop
   else
      this=>this%first_id%up_tag%first_tag%next_tag 
      do
         if(this%tag.xmleq.tag) then
            count_=count_+1
            if(count_.eq.count) exit
         end if
         if(associated(this%next_tag)) then
            this=>this%next_tag
         else
            print*,'ERROR INLINKLIST SELECTTAG: COULD NOT FIND THE ',count,'. tag: ',tag,' IN THE LAYER ABOVE &
                 &THE ACTUAL!'
            stop
         end if
      end do
   end if
else
   print*,'ERROR STOP IN LINKLIST: LL_SELECTTAG ONLY ACCEPTS <UP>, <DOWN> OR <ACT>&
        &AND NOT: ',layer
   stop
end if
return
end subroutine libxml2f90__ll_selecttag



!=========================================
!========= these call the generic interface
!=========================================

!============================================
subroutine libxml2f90__existid(id,texist)
  implicit none
  character(*),intent(in)        :: id
  logical(4),intent(out)         :: texist
  !............................
  texist=.false.
  call libxml2f90_existid(1,id,texist)!1=id 0=pid
  return
end subroutine libxml2f90__existid

!============================================
subroutine libxml2f90__existpid(id,texist)
  implicit none
  character(*),intent(in)        :: id
  logical(4),intent(out)         :: texist
  !............................
  texist=.false.
  call libxml2f90_existid(0,id,texist)!1=id 0=pid
  return
end subroutine libxml2f90__existpid


!================================================================
subroutine libxml2f90__ll_getr8(id,size_,val)
implicit none
character(*),intent(in)         :: id
integer(4),intent(in)           :: size_
real(8),intent(out)             :: val(size_)
!........................................
call libxml2f90_ll_getr8(1,id,size_,val)
return
end subroutine libxml2f90__ll_getr8

!================================================================
subroutine libxml2f90__ll_getr8_(id,val)
implicit none
character(*),intent(in)         :: id
real(8),intent(out)             :: val
real(8)                         :: vala(1)
!........................................
call libxml2f90_ll_getr8(1,id,1,vala)
val=vala(1)
return
end subroutine libxml2f90__ll_getr8_


!================================================================
subroutine libxml2f90__ll_getc8(id,size_,val)
implicit none
character(*),intent(in)         :: id
integer(4),intent(in)           :: size_
complex(8),intent(out)          :: val(size_)
!........................................
call libxml2f90_ll_getc8(1,id,size_,val)
return
end subroutine libxml2f90__ll_getc8

!================================================================
subroutine libxml2f90__ll_getc8_(id,val)
implicit none
character(*),intent(in)         :: id
complex(8),intent(out)          :: val
complex(8)                      :: vala(1)
!........................................
call libxml2f90_ll_getc8(1,id,1,vala)
val=vala(1)
return
end subroutine libxml2f90__ll_getc8_


!================================================================
subroutine libxml2f90__ll_geti4(id,size_,val)
implicit none
character(*),intent(in)         :: id
integer(4),intent(in)           :: size_
integer(4),intent(out)          :: val(size_)
!........................................
call libxml2f90_ll_geti4(1,id,size_,val)
return
end subroutine libxml2f90__ll_geti4

!================================================================
subroutine libxml2f90__ll_geti4_(id,val)
implicit none
character(*),intent(in)         :: id
integer(4),intent(out)          :: val
integer(4)                      :: vala(1)
!........................................
call libxml2f90_ll_geti4(1,id,1,vala)
val=vala(1)
return
end subroutine libxml2f90__ll_geti4_



!================================================================
subroutine libxml2f90__ll_getl4(id,size_,val)
implicit none
character(*),intent(in)         :: id
integer(4),intent(in)           :: size_
logical(4),intent(out)          :: val(size_)
!........................................
call libxml2f90_ll_getl4(1,id,size_,val)
return
end subroutine libxml2f90__ll_getl4

!================================================================
subroutine libxml2f90__ll_getl4_(id,val)
implicit none
character(*),intent(in)         :: id
logical(4),intent(out)          :: val
logical(4)                      :: vala(1)
!........................................
call libxml2f90_ll_getl4(1,id,1,vala)
val=vala(1)
return
end subroutine libxml2f90__ll_getl4_


!================================================================
subroutine libxml2f90__ll_getch(id,size_,val)
implicit none
character(*),intent(in)         :: id
integer(4),intent(in)           :: size_
character(1),intent(out)        :: val(size_)
!........................................
call libxml2f90_ll_getch(1,id,size_,val)
return
end subroutine libxml2f90__ll_getch

!================================================================
!    MLO.  Added this function for the case that val is a scalar in the call.  If it is
!          not, and we use libxml2f90__ll_getch, then the interface check will fail
subroutine libxml2f90__ll_getch_scal(id,size_,val)
implicit none
character(len=*),intent(in)        :: id
integer(4),intent(in)              :: size_
character(len=size_),intent(out)   :: val
character(len=1), dimension(size_) :: val_loc

!........................................
call libxml2f90_ll_getch(1,id,size_,val_loc)
write (val(1:size_),'(a)') val_loc(1:size_)
return
end subroutine libxml2f90__ll_getch_scal

!================================================================
subroutine libxml2f90__ll_getstring(id,size_,val)
implicit none
character(*),intent(in)         :: id
integer(4),intent(in)           :: size_
character(*),intent(out)        :: val(size_)
!........................................
call libxml2f90_ll_getstring(1,id,size_,val)
return
end subroutine libxml2f90__ll_getstring

!================================================================
subroutine libxml2f90__ll_getstring_(id,val)
implicit none
character(*),intent(in)         :: id
character(*),intent(out)        :: val
character(512)                  :: vala(1)
!........................................
call libxml2f90_ll_getstring(1,id,1,vala)
val=trim(adjustl(vala(1)))
return
end subroutine libxml2f90__ll_getstring_

!================================================================
subroutine libxml2f90__ll_getsize(id,size_)
implicit none
character(*),intent(in)         :: id
integer(4),intent(out)          :: size_
!........................................
call libxml2f90_ll_getsize(1,id,size_)
return
end subroutine libxml2f90__ll_getsize


!================================================================
subroutine libxml2f90__ll_getpr8(id,size_,val)
implicit none
character(*),intent(in)         :: id
integer(4),intent(in)           :: size_
real(8),intent(out)             :: val(size_)
!........................................
call libxml2f90_ll_getr8(0,id,size_,val)
return
end subroutine libxml2f90__ll_getpr8

!================================================================
subroutine libxml2f90__ll_getpr8_(id,val)
implicit none
character(*),intent(in)         :: id
real(8),intent(out)             :: val
real(8)                         :: vala(1)
!........................................
call libxml2f90_ll_getr8(0,id,1,vala)
val=vala(1)
return
end subroutine libxml2f90__ll_getpr8_


!================================================================
subroutine libxml2f90__ll_getpc8(id,size_,val)
implicit none
character(*),intent(in)         :: id
integer(4),intent(in)           :: size_
complex(8),intent(out)          :: val(size_)
!........................................
call libxml2f90_ll_getc8(0,id,size_,val)
return
end subroutine libxml2f90__ll_getpc8

!================================================================
subroutine libxml2f90__ll_getpc8_(id,val)
implicit none
character(*),intent(in)         :: id
complex(8),intent(out)          :: val
complex(8)                      :: vala(1)
!........................................
call libxml2f90_ll_getc8(0,id,1,vala)
val=vala(1)
return
end subroutine libxml2f90__ll_getpc8_


!================================================================
subroutine libxml2f90__ll_getpi4(id,size_,val)
implicit none
character(*),intent(in)         :: id
integer(4),intent(in)           :: size_
integer(4),intent(out)          :: val(size_)
!........................................
call libxml2f90_ll_geti4(0,id,size_,val)
return
end subroutine libxml2f90__ll_getpi4

!================================================================
subroutine libxml2f90__ll_getpi4_(id,val)
implicit none
character(*),intent(in)         :: id
integer(4),intent(out)          :: val
integer(4)                      :: vala(1)
!........................................
call libxml2f90_ll_geti4(0,id,1,vala)
val=vala(1)
return
end subroutine libxml2f90__ll_getpi4_


!================================================================
subroutine libxml2f90__ll_getpl4(id,size_,val)
implicit none
character(*),intent(in)         :: id
integer(4),intent(in)           :: size_
logical(4),intent(out)          :: val(size_)
!........................................
call libxml2f90_ll_getl4(0,id,size_,val)
return
end subroutine libxml2f90__ll_getpl4

!================================================================
subroutine libxml2f90__ll_getpl4_(id,val)
implicit none
character(*),intent(in)         :: id
logical(4),intent(out)          :: val
logical(4)                      :: vala(1)
!........................................
call libxml2f90_ll_getl4(0,id,1,vala)
val=vala(1)
return
end subroutine libxml2f90__ll_getpl4_


!================================================================
subroutine libxml2f90__ll_getpch(id,size_,val)
implicit none
character(*),intent(in)         :: id
integer(4),intent(in)           :: size_
character(1),intent(out)        :: val(size_)
!........................................
call libxml2f90_ll_getch(0,id,size_,val)
return
end subroutine libxml2f90__ll_getpch

!================================================================
subroutine libxml2f90__ll_getpstring(id,size_,val)
implicit none
character(*),intent(in)         :: id
integer(4),intent(in)           :: size_
character(*),intent(out)        :: val(size_)
!........................................
call libxml2f90_ll_getstring(0,id,size_,val)
return
end subroutine libxml2f90__ll_getpstring

!================================================================
subroutine libxml2f90__ll_getpstring_(id,val)
implicit none
character(*),intent(in)         :: id
character(*),intent(out)        :: val
character(512)                  :: vala(1)
!........................................
call libxml2f90_ll_getstring(0,id,1,vala)
val=trim(adjustl(vala(1)))
return
end subroutine libxml2f90__ll_getpstring_

!================================================================
subroutine libxml2f90__ll_getpsize(id,size_)
implicit none
character(*),intent(in)         :: id
integer(4),intent(out)          :: size_
!........................................
call libxml2f90_ll_getsize(0,id,size_)
return
end subroutine libxml2f90__ll_getpsize

!========= these call the generic interface
!=========================================


!============================================
!================= general readinterface
!============================================

subroutine libxml2f90_existid(wid,id,texist)
  use ll_module
  implicit none
  integer(4),intent(in)          :: wid
  character(*),intent(in)        :: id
  logical(4),intent(out)         :: texist
  !............................
  texist=.false.

  THIS=>THIS%FIRST_ID
  if(wid.eq.0) then !read pids
     do
        if(associated(THIS%NEXT_PID)) then
           THIS=>THIS%NEXT_PID
           if(THIS%ID.xmleq.ID) then
              texist=.true.
              exit
           end if
        else
           exit
        end if
     end do
  else if(wid.eq.1) then !read ids
     do
        if(associated(THIS%NEXT_ID)) then
           THIS=>THIS%NEXT_ID
           if(THIS%ID.xmleq.ID) then
              texist=.true.
              exit
           end if
        else
           exit
        end if
     end do
  else
     stop 'libxml2f90_ll_get interface read method not allowed'
  end if
  return
end subroutine libxml2f90_existid



!================================================================
!=================== THE GET SAFE ROUTINES ======================
!================================================================

!================================================================
subroutine libxml2f90_getsafer8(id,csize,chara,size_,value)
  implicit none
  character(*),intent(in)        :: id
  integer(4),intent(in)          :: csize    
  character(1),intent(in)        :: chara(csize)
  integer(4),intent(in)          :: size_
  real(8),intent(out)            :: value(size_) 
  integer(4)                     :: i,isign,blocks,left,right
  character(256)                 :: string
  integer(4)                     :: stringsize=256
  character,dimension(1)         :: ctmp_sp 
  ctmp_sp = (/' '/)

  !=======================================
  !        check for unallowed signs
  !=======================================
  ! for real these are everything except
  ! space, e,d,.,+,-
  do i=1,csize
     !check for unallowed signs
     if(scan(chara(i),'1234567890.dD eE+-').ne.1) then
        print*,'ERROR STOP: UNALLOWED SIGN WHILE TRYING TO READ&
             &NUMERICAL INPUT: ',id,chara(i)
     end if
  end do


  !=======================================
  !           split in blocks
  !=======================================
  left=1
  right=1
  blocks=0
 
  do
     !find the first non blank
     do
        if(chara(left).eq.' '.and.left.lt.csize) then
           left=left+1
        else
           exit
        end if
     end do
     
     call libxml2f90__findinchara(csize,chara,left,.true.,1,ctmp_sp,right,isign)
     if(isign.eq.0) right=csize !no blank found - we use the right bound
     if(right.lt.left) then
        print*,'ERROR STOP IN LINKELIST: GETSAFE: RIGHTBOUND IS SMALLER THAN LEFTBOUND'
        stop
     end if
     
     if(left.eq.right.and.scan(chara(left),' .dD eE+-').eq.1) then
        if(chara(left).eq.' ') then
           !we skip the blank
        else
           print*,'ERROR STOP IN LINKLIST: GETSAFE NUMERICAL INPUT CONSISTS SOLELY OF A ',chara(left)
           stop
        end if
     else
        !we read
        if(blocks+1.gt.size_) then
           print*,'ERROR STOP IN LINKLIST: GETSAFE: THE NUMERICAL INPUT IS LARGER THAN THE &
                &PROVIDED ARRAY! THE ID WAS: ',trim(adjustl(id))
           stop
        else
           !if(allocated(val))deallocate(val)
           !allocate(val(right-left+1))
           if((right-left+1).gt.stringsize) then
              print*,'ERROR STOP IN LINKLIST: GETSAFE: THE INTERNAL FIXED STRING SIZE IS &
                   &TOO SMALL! U CAN CHANGE THE VALUE IN libxml2f90_getsafe*. ID WAS: ',trim(adjustl(id))
              stop
           end if
           call libxml2f90_tostring((right-left+1),chara(left:right),string)
           blocks=blocks+1
           string=trim(adjustl(string))
           read(string,*)value(blocks)
        end if
     end if
     left=right
     if(right.eq.csize) exit !the loop
  end do
  if(blocks.lt.size_) then
     print*,'ERROR STOP IN LINKLIST: GETSAFE: # OF BLOCKS READ IS SMALLER THAN&
          & USER WANTS TO READ'
     stop
  end if
  return
end subroutine libxml2f90_getsafer8

!================================================================
subroutine libxml2f90_getsafec8(id,csize,chara,size_,value)
  implicit none
  character(*),intent(in)        :: id
  integer(4),intent(in)          :: csize    
  character(1),intent(in)        :: chara(csize)
  integer(4),intent(in)          :: size_
  complex(8),intent(out)         :: value(size_) 
  integer(4)                     :: i,isign,blocks,left,right
  character(256)                 :: string
  integer(4)                     :: stringsize=256
  character, dimension(1)        :: ctmp_sp 
  ctmp_sp = (/' '/)
  
  !=======================================
  !        check for unallowed signs
  !=======================================
  ! for real these are everything except
  ! space, e,d,.,+,-
  do i=1,csize
     !check for unallowed signs
     if(scan(chara(i),'1234567890.dD eE+-').ne.1) then
        print*,'ERROR STOP: UNALLOWED SIGN WHILE TRYING TO READ&
             &NUMERICAL INPUT: ',id,chara(i)
     end if
  end do


  !=======================================
  !           split in blocks
  !=======================================
  left=1
  right=1
  blocks=0
 
  do
     !find the first non blank
     do
        if(chara(left).eq.' '.and.left.lt.csize) then
           left=left+1
        else
           exit
        end if
     end do
     
     call libxml2f90__findinchara(csize,chara,left,.true.,1,ctmp_sp,right,isign)
     if(isign.eq.0) right=csize !no blank found - we use the right bound
     if(right.lt.left) then
        print*,'ERROR STOP IN LINKELIST: GETSAFE: RIGHTBOUND IS SMALLER THAN LEFTBOUND'
        stop
     end if
     
     if(left.eq.right.and.scan(chara(left),' .dD eE+-').eq.1) then
        if(chara(left).eq.' ') then
           !we skip the blank
        else
           print*,'ERROR STOP IN LINKLIST: GETSAFE NUMERICAL INPUT CONSISTS SOLELY OF A ',chara(left)
           stop
        end if
     else
        !we read
        if(blocks+1.gt.size_) then
           print*,'ERROR STOP IN LINKLIST: GETSAFE: THE NUMERICAL INPUT IS LARGER THAN THE &
                &PROVIDED ARRAY! THE ID WAS: ',trim(adjustl(id))
           stop
        else
           !if(allocated(val))deallocate(val)
           !allocate(val(right-left+1))
           if((right-left+1).gt.stringsize) then
              print*,'ERROR STOP IN LINKLIST: GETSAFE: THE INTERNAL FIXED STRING SIZE IS &
                   &TOO SMALL! U CAN CHANGE THE VALUE IN libxml2f90_getsafe*. ID WAS: ',trim(adjustl(id))
              stop
           end if
           call libxml2f90_tostring((right-left+1),chara(left:right),string)
           blocks=blocks+1
           string=trim(adjustl(string))
           read(string,*)value(blocks)
        end if
     end if
     left=right
     if(right.eq.csize) exit !the loop
  end do
  if(blocks.lt.size_) then
     print*,'ERROR STOP IN LINKLIST: GETSAFE: # OF BLOCKS READ IS SMALLER THAN&
          & USER WANTS TO READ'
     stop
  end if
  return
end subroutine libxml2f90_getsafec8

!================================================================
subroutine libxml2f90_getsafei4(id,csize,chara,size_,value)
  implicit none
  character(*),intent(in)        :: id
  integer(4),intent(in)          :: csize    
  character(1),intent(in)        :: chara(csize)
  integer(4),intent(in)          :: size_
  integer(4),intent(out)         :: value(size_) 
  integer(4)                     :: i,isign,blocks,left,right
  character(256)                 :: string
  integer(4)                     :: stringsize=256
  character,dimension(1)         :: ctmp_sp 
  ctmp_sp = (/' '/)
  
  !=======================================
  !        check for unallowed signs
  !=======================================
  do i=1,csize
     !check for unallowed signs
     if(scan(chara(i),'1234567890 +-').ne.1) then
        print*,'ERROR STOP: UNALLOWED SIGN WHILE TRYING TO READ&
             &NUMERICAL INPUT: ',id,chara(i)
     end if
  end do


  !=======================================
  !           split in blocks
  !=======================================
  left=1
  right=1
  blocks=0
 
  do
     !find the first non blank
     do
        if(chara(left).eq.' '.and.left.lt.csize) then
           left=left+1
        else
           exit
        end if
     end do
     
     call libxml2f90__findinchara(csize,chara,left,.true.,1,ctmp_sp,right,isign)
     if(isign.eq.0) right=csize !no blank found - we use the right bound
     if(right.lt.left) then
        print*,'ERROR STOP IN LINKELIST: GETSAFE: RIGHTBOUND IS SMALLER THAN LEFTBOUND'
        stop
     end if
     
     if(left.eq.right.and.scan(chara(left),' +-').eq.1) then
        if(chara(left).eq.' ') then
           !we skip the blank
        else
           print*,'ERROR STOP IN LINKLIST: GETSAFE NUMERICAL INPUT CONSISTS SOLELY OF A ',chara(left)
           stop
        end if
     else
        !we read
        if(blocks+1.gt.size_) then
           print*,'ERROR STOP IN LINKLIST: GETSAFE: THE NUMERICAL INPUT IS LARGER THAN THE &
                &PROVIDED ARRAY! THE ID WAS: ',trim(adjustl(id))
           stop
        else
           !if(allocated(val))deallocate(val)
           !allocate(val(right-left+1))
           if((right-left+1).gt.stringsize) then
              print*,'ERROR STOP IN LINKLIST: GETSAFE: THE INTERNAL FIXED STRING SIZE IS &
                   &TOO SMALL! U CAN CHANGE THE VALUE IN libxml2f90_getsafe*. ID WAS: ',trim(adjustl(id))
              stop
           end if
           call libxml2f90_tostring((right-left+1),chara(left:right),string)
           blocks=blocks+1
           string=trim(adjustl(string))
           read(string,*)value(blocks)
        end if
     end if
     left=right
     if(right.eq.csize) exit !the loop
  end do
  if(blocks.lt.size_) then
     print*,'ERROR STOP IN LINKLIST: GETSAFE: # OF BLOCKS READ IS SMALLER THAN&
          & USER WANTS TO READ'
     stop
  end if
  return
end subroutine libxml2f90_getsafei4

!================================================================
subroutine libxml2f90_getsafel4(id,csize,chara,size_,value)
  implicit none
  character(*),intent(in)        :: id
  integer(4),intent(in)          :: csize    
  character(1),intent(in)        :: chara(csize)
  integer(4),intent(in)          :: size_
  logical(4),intent(out)         :: value(size_) 
  integer(4)                     :: i,isign,blocks,left,right
  character(256)                 :: string
  integer(4)                     :: stringsize=256
  character, dimension(1)        :: ctmp_sp 
  ctmp_sp = (/' '/)
  
  !=======================================
  !        check for unallowed signs
  !=======================================
  do i=1,csize
     !check for unallowed signs
     if(scan(chara(i),' .truefalsTRUEFALS').ne.1) then
        print*,'ERROR STOP: UNALLOWED SIGN WHILE TRYING TO READ&
             &NUMERICAL INPUT: ',id,' ',chara(i)
     end if
  end do


  !=======================================
  !           split in blocks
  !=======================================
  left=1
  right=1
  blocks=0
 
  do
     !find the first non blank
     do
        if(chara(left).eq.' '.and.left.lt.csize) then
           left=left+1
        else
           exit
        end if
     end do
     
     call libxml2f90__findinchara(csize,chara,left,.true.,1,ctmp_sp,right,isign)
     if(isign.eq.0) right=csize !no blank found - we use the right bound
     if(right.lt.left) then
        print*,'ERROR STOP IN LINKELIST: GETSAFE: RIGHTBOUND IS SMALLER THAN LEFTBOUND'
        stop
     end if
     
     if(left.eq.right.and.scan(chara(left),' .RUEALSrueals').eq.1) then
        if(chara(left).eq.' ') then
           !we skip the blank
        else
           print*,'ERROR STOP IN LINKLIST: GETSAFE NUMERICAL INPUT CONSISTS SOLELY OF A ',chara(left)
           stop
        end if
     else
        !we read
        if(blocks+1.gt.size_) then
           print*,'ERROR STOP IN LINKLIST: GETSAFE: THE NUMERICAL INPUT IS LARGER THAN THE &
                &PROVIDED ARRAY! THE ID WAS: ',trim(adjustl(id))
           stop
        else
           !if(allocated(val))deallocate(val)
           !allocate(val(right-left+1))
           if((right-left+1).gt.stringsize) then
              print*,'ERROR STOP IN LINKLIST: GETSAFE: THE INTERNAL FIXED STRING SIZE IS &
                   &TOO SMALL! U CAN CHANGE THE VALUE IN libxml2f90_getsafe*. ID WAS: ',trim(adjustl(id))
              stop
           end if
           call libxml2f90_tostring((right-left+1),chara(left:right),string)
           blocks=blocks+1
           string=trim(adjustl(string))
           read(string,*)value(blocks)
        end if
     end if
     left=right
     if(right.eq.csize) exit !the loop
  end do
  if(blocks.lt.size_) then
     print*,'ERROR STOP IN LINKLIST: GETSAFE: # OF BLOCKS READ IS SMALLER THAN&
          & USER WANTS TO READ'
     stop
  end if
  return
end subroutine libxml2f90_getsafel4


!================================================================
!=============== END OF THE GET SAFE ROUTINES ===================
!================================================================



!================================================================
subroutine libxml2f90_ll_getr8(wid,id,size_,val)
use ll_module
implicit none
integer(4),intent(in)           :: wid! 0=pid,1=id
character(*),intent(in)         :: id
integer(4),intent(in)           :: size_
real(8),intent(out)             :: val(size_)
character,dimension(:),allocatable  :: lvalue
!........................................

THIS=>THIS%FIRST_ID
if(wid.eq.0) then !read pids
   do
      if(associated(THIS%NEXT_PID)) then
         THIS=>THIS%NEXT_PID
         if(THIS%ID.xmleq.ID) then
            !the safe way:
            call libxml2f90_getsafer8(trim(id),size(THIS%VALUE),THIS%VALUE,size_,val)
            exit
         end if
      else
         print*,'ERROR STOP IN LINKLIST: CAN NOT GET R8 VALUE FOR PID: ',trim(adjustl(id))
         stop
      end if
   end do
else if(wid.eq.1) then !read ids
   do
      if(associated(THIS%NEXT_ID)) then
         THIS=>THIS%NEXT_ID
         if(THIS%ID.xmleq.ID) then
            !the safe way:
            allocate(lvalue(size(THIS%VALUE)))
            lvalue = THIS%VALUE
            call libxml2f90_getsafer8(trim(id),size(THIS%VALUE),LVALUE,size_,val)
            deallocate(lvalue)
            exit
         end if
      else
         print*,'ERROR STOP IN LINKLIST: CAN NOT GET R8 VALUE FOR ID: ',trim(adjustl(id))
         stop
      end if
   end do
else
   stop 'libxml2f90_ll_get interface read method not allowed'
end if
return
end subroutine libxml2f90_ll_getr8


!================================================================
subroutine libxml2f90_ll_getc8(wid,id,size_,val)
use ll_module
implicit none
integer(4),intent(in)           :: wid! 0=pid,1=id
character(*),intent(in)         :: id
integer(4),intent(in)           :: size_
complex(8),intent(out)          :: val(size_)
!character(30000)                :: string='' !solve more elegant!!!
!........................................

THIS=>THIS%FIRST_ID
if(wid.eq.0) then !read pids
   do
      if(associated(THIS%NEXT_PID)) then
         THIS=>THIS%NEXT_PID
         if(THIS%ID.xmleq.ID) then
            !the safe way:
            call libxml2f90_getsafec8(trim(id),size(THIS%VALUE),THIS%VALUE,size_,val)

            !alex: this has to be solved more elegant!
 !           if(len(string).lt.size(THIS%VALUE)) stop'ERROR VALUE IS LARGER THAN IMPLEMENTED'
 !           call libxml2f90_tostring(size(THIS%VALUE),THIS%VALUE,string)
 !           read(string,*)val
            exit
         end if
      else
         print*,'ERROR STOP IN LINKLIST: CAN NOT GET C8 VALUE FOR ID: ',trim(adjustl(id))
         stop
      end if
   end do
else if(wid.eq.1) then !read pids
   do
      if(associated(THIS%NEXT_ID)) then
         THIS=>THIS%NEXT_ID
         if(THIS%ID.xmleq.ID) then
            call libxml2f90_getsafec8(trim(id),size(THIS%VALUE),THIS%VALUE,size_,val)
           !alex: this has to be solved more elegant!
!            if(len(string).lt.size(THIS%VALUE)) stop'ERROR VALUE IS LARGER THAN IMPLEMENTED'
!            call libxml2f90_tostring(size(THIS%VALUE),THIS%VALUE,string)
!            read(string,*)val
            exit
         end if
      else
         print*,'ERROR STOP IN LINKLIST: CAN NOT GET C8 VALUE FOR ID: ',trim(adjustl(id))
         stop
      end if
   end do
else
   stop 'libxml2f90_ll_get interface read method not allowed'
end if

return
end subroutine libxml2f90_ll_getc8


!================================================================
subroutine libxml2f90_ll_geti4(wid,id,size_,val)
use ll_module
implicit none
integer(4),intent(in)           :: wid! 0=pid,1=id
character(*),intent(in)         :: id
integer(4),intent(in)           :: size_
integer(4),intent(out)          :: val(size_)
character,dimension(:),allocatable  :: lvalue
!character(30000)                :: string='' !solve more elegant!!!
!........................................

THIS=>THIS%FIRST_ID
if(wid.eq.0) then !read pids
   do
      if(associated(THIS%NEXT_PID)) then
         THIS=>THIS%NEXT_PID
         if(THIS%ID.xmleq.ID) then
            call libxml2f90_getsafei4(trim(id),size(THIS%VALUE),THIS%VALUE,size_,val)
            !alex: this has to be solved more elegant!
!            if(len(string).lt.size(THIS%VALUE)) stop'ERROR VALUE IS LARGER THAN IMPLEMENTED'
!            call libxml2f90_tostring(size(THIS%VALUE),THIS%VALUE,string)
!            read(string,*)val
            exit
         end if
      else
         print*,'ERROR STOP IN LINKLIST: CAN NOT GET I4 VALUE FOR ID: ',trim(adjustl(id))
         stop
      end if
   end do
else if(wid.eq.1) then !read pids
   do
      if(associated(THIS%NEXT_ID)) then
         THIS=>THIS%NEXT_ID
         if(THIS%ID.xmleq.ID) then
            allocate(lvalue(size(THIS%VALUE)))
            lvalue = THIS%VALUE
            call libxml2f90_getsafei4(trim(id),size(THIS%VALUE),LVALUE,size_,val)
            deallocate(lvalue)
            !alex: this has to be solved more elegant!
!            if(len(string).lt.size(THIS%VALUE)) stop'ERROR VALUE IS LARGER THAN IMPLEMENTED'
!            call libxml2f90_tostring(size(THIS%VALUE),THIS%VALUE,string)
!            read(string,*)val
            exit
         end if
      else
         print*,'ERROR STOP IN LINKLIST: CAN NOT GET I4 VALUE FOR ID: ',trim(adjustl(id))
         stop
      end if
   end do
else
   stop 'libxml2f90_ll_get interface read method not allowed'
end if

return
end subroutine libxml2f90_ll_geti4



!================================================================
subroutine libxml2f90_ll_getl4(wid,id,size_,val)
use ll_module
implicit none
integer(4),intent(in)           :: wid! 0=pid,1=id
character(*),intent(in)         :: id
integer(4),intent(in)           :: size_
logical(4),intent(out)          :: val(size_)
!character(30000)                :: string='' !solve more elegant!!!
!........................................

THIS=>THIS%FIRST_ID
if(wid.eq.0) then !read pids
   do
      if(associated(THIS%NEXT_PID)) then
         THIS=>THIS%NEXT_PID
         if(THIS%ID.xmleq.ID) then
            call libxml2f90_getsafel4(trim(id),size(THIS%VALUE),THIS%VALUE,size_,val)
            !alex: this has to be solved more elegant!
 !           if(len(string).lt.size(THIS%VALUE)) stop'ERROR VALUE IS LARGER THAN IMPLEMENTED'
 !           call libxml2f90_tostring(size(THIS%VALUE),THIS%VALUE,string)
 !           read(string,*)val
            exit
         end if
      else
         print*,'ERROR STOP IN LINKLIST: CAN NOT GET L4 VALUE FOR ID: ',trim(adjustl(id))
         stop
      end if
   end do
else if(wid.eq.1) then !read pids
   do
      if(associated(THIS%NEXT_ID)) then
         THIS=>THIS%NEXT_ID
         if(THIS%ID.xmleq.ID) then
            call libxml2f90_getsafel4(trim(id),size(THIS%VALUE),THIS%VALUE,size_,val)
            !alex: this has to be solved more elegant!
            !if(len(string).lt.size(THIS%VALUE)) stop'ERROR VALUE IS LARGER THAN IMPLEMENTED'
            !call libxml2f90_tostring(size(THIS%VALUE),THIS%VALUE,string)
            !read(string,*)val
            exit
         end if
      else
         print*,'ERROR STOP IN LINKLIST: CAN NOT GET L4 VALUE FOR ID: ',trim(adjustl(id))
         stop
      end if
   end do
else
   stop 'libxml2f90_ll_get interface read method not allowed'
end if

return
end subroutine libxml2f90_ll_getl4


!================================================================
subroutine libxml2f90_ll_getch(wid,id,size_,val)
use ll_module
implicit none
integer(4),intent(in)           :: wid! 0=pid,1=id
character(*),intent(in)         :: id
integer(4),intent(in)           :: size_
character(1),intent(out)        :: val(size_)
!........................................

THIS=>THIS%FIRST_ID
if(wid.eq.0) then !read pids
   do
      if(associated(THIS%NEXT_PID)) then
         THIS=>THIS%NEXT_PID
         if(THIS%ID.xmleq.ID) then
            if(size(this%value).ne.size_) then
               print*,'ERROR STOP IN LINKLIST: SIZE OF CH(1) ARRAY IN GETCH DOES NOT FIT: '&
                    &,SIZE_,size(this%value)
            end if
            val(:)=THIS%VALUE(:)
            exit
         end if
      else
         print*,'ERROR STOP IN LINKLIST: CAN NOT GET CH VALUE FOR ID: ',trim(adjustl(id))
         stop
      end if
   end do
else if(wid.eq.1) then !read pids
   do
      if(associated(THIS%NEXT_ID)) then
         THIS=>THIS%NEXT_ID
         if(THIS%ID.xmleq.ID) then
            if(size(this%value).ne.size_) then
               print*,'ERROR STOP IN LINKLIST: SIZE OF CH(1) ARRAY IN GETCH DOES NOT FIT: '&
                    &,SIZE_,size(this%value)
            end if
            val(:)=THIS%VALUE(:)
            exit
         end if
      else
         print*,'ERROR STOP IN LINKLIST: CAN NOT GET CH VALUE FOR ID: ',trim(adjustl(id))
         stop
      end if
   end do
else
   stop 'libxml2f90_ll_get interface read method not allowed'
end if

return
end subroutine libxml2f90_ll_getch

!================================================================
subroutine libxml2f90_ll_getstring(wid,id,size_,val)
use ll_module
implicit none
integer(4),intent(in)           :: wid! 0=pid,1=id
character(*),intent(in)         :: id
integer(4),intent(in)           :: size_
character(*),intent(out)        :: val(size_)
!........................................

THIS=>THIS%FIRST_ID
if(wid.eq.0) then !read pids
   do
      if(associated(THIS%NEXT_PID)) then
         THIS=>THIS%NEXT_PID
         if(THIS%ID.xmleq.ID) then
            if(size(this%value).gt.(size_*len(val(1)))) then
               print*,'ERROR STOP IN LINKLIST: SIZE OF STRING ARRAY IN GETSTRING &
                    &IS TOO SMALL: ',SIZE_*len(val(1)),size(this%value)
               stop
            end if
            call libxml2f90_tostringa(size(this%value),this%value,size_,val)
            exit
         end if
      else
         print*,'ERROR STOP IN LINKLIST: CAN NOT GET CH VALUE FOR ID: ',trim(adjustl(id))
         stop
      end if
   end do
else if(wid.eq.1) then !read pids
   do
      if(associated(THIS%NEXT_ID)) then
         THIS=>THIS%NEXT_ID
         if(THIS%ID.xmleq.ID) then
            if(size(this%value).gt.(size_*len(val(1)))) then
               print*,'ERROR STOP IN LINKLIST: SIZE OF STRING ARRAY IN GETSTRING &
                    &IS TOO SMALL: ',SIZE_*len(val(1)),size(this%value)
               stop
            end if
            call libxml2f90_tostringa(size(this%value),this%value,size_,val)
            exit
         end if
      else
         print*,'ERROR STOP IN LINKLIST: CAN NOT GET CH VALUE FOR ID: ',trim(adjustl(id))
         stop
      end if
   end do
else
   stop 'libxml2f90_ll_get interface read method not allowed'
end if

return
end subroutine libxml2f90_ll_getstring

!================================================================
subroutine libxml2f90_ll_getsize(wid,id,size_)
use ll_module
implicit none
integer(4),intent(in)           :: wid! 0=pid,1=id
character(*),intent(in)         :: id
integer(4),intent(out)          :: size_
!........................................

THIS=>THIS%FIRST_ID
if(wid.eq.0) then !read pids
   do
      if(associated(THIS%NEXT_PID)) then
         THIS=>THIS%NEXT_PID
         if(THIS%ID.xmleq.ID) then
            size_=size(this%value)
            exit
         end if
      else
         print*,'ERROR STOP IN LINKLIST: CAN NOT GET SIZE OF VALUE FOR ID: ',trim(adjustl(id))
         stop
      end if
   end do
else if(wid.eq.1) then !read pids
   do
      if(associated(THIS%NEXT_ID)) then
         THIS=>THIS%NEXT_ID
         if(THIS%ID.xmleq.ID) then
            size_=size(this%value)
            exit
         end if
      else
         print*,'ERROR STOP IN LINKLIST: CAN NOT GET SIZE OF VALUE FOR ID: ',trim(adjustl(id))
         stop
      end if
   end do
else
   stop 'libxml2f90_ll_get interface read method not allowed'
end if

return
end subroutine libxml2f90_ll_getsize

!=================end general readinterface
!============================================

!================================================================
subroutine libxml2f90_tostring(size_,val,string)
  use ll_module
  implicit none
  integer(4),intent(in)              :: size_
  character(1),intent(in)            :: val(size_)
  character(*),intent(out)           :: string
  integer(4)                         :: i
  !.............................................
  string=''
  do i=1,size_
     string(i:i)=val(i)
  end do
  return
end subroutine libxml2f90_tostring

!================================================================
subroutine libxml2f90_tostringa(size_,val,sizea,stringa)
  use ll_module
  implicit none
  integer(4),intent(in)              :: size_
  character(1),intent(in)            :: val(size_)
  integer(4),intent(in)              :: sizea
  character(*),intent(out)           :: stringa(sizea)
  integer(4)                         :: i,j,k
  logical(4)                         :: texit
  !.............................................
  texit=.false.
  stringa(:)=''
  do i=1,sizea                     !loop over the array layers
     do j=1,len(stringa(1))        !loop over the layer itself
        k=(i-1)*len(stringa(1))+j  !get the actual character(1) value
        if(k.gt.size_) then
           !(we do not need this warning)
           !print*,'WARNING IN LINKLIST: TOSTRINGA: &
           !     &YOUR STRINGA IS LARGER THEN CH(1) ARRAY'
           texit=.true.
           exit
        end if
        stringa(i)(j:j)=val(k)
     end do
     if(texit) exit
  end do
  if(k.lt.size_) then
     print*,'ERROR STOP IN LINKLIST: TOSTRINGA: &
          &YOUR STRINGA IS SMALLER THEN CH(1) ARRAY'
     stop
  end if
  return
end subroutine libxml2f90_tostringa



!!__!================================================================
!!__subroutine libxml2f90$ll_getvalue(id,value,tget)
!!__use ll_module
!!__implicit none
!!__character(*),intent(in)         :: id
!!__character(*),intent(out)        :: value
!!__logical(4),intent(out)          :: tget !did we manage to find that id
!!__!........................................
!!__tget=.false.
!!__THIS=>THIS%FIRST_ID
!!__do
!!__   if(associated(THIS%NEXT_ID)) then
!!__      THIS=>THIS%NEXT_ID
!!__!print*,'debug-',trim(this%id)
!!__!print*,'debug--',trim(this%value)
!!__      if(trim(THIS%ID).eq.trim(ID)) then
!!__         value=this%value
!!__         tget=.true.
!!__         exit
!!__      end if
!!__   else
!!__      exit
!!__   end if
!!__end do
!!__
!!__return
!!__end subroutine libxml2f90$ll_getvalue



!================================================================
subroutine libxml2f90_ll_add_list(ll_id)
  use ll_module
  implicit none
  character(*),intent(in)       :: ll_id !the id of the first linklist
  !.........................
  if(.NOT.INITIALIZED) then
     call libxml2f90_ll_initlist(ll_id)
     return
  end if

  THIS=>LL_ROOT
  do
     if(associated(THIS%NEXT_LL)) then
        THIS=>THIS%NEXT_LL
        if(trim(THIS%LL_ID).eq.trim(LL_ID)) then
           print*,'ERROR STOP IN LINKLIST: YOU CAN NOT USE AN ID TWICE: ',LL_ID
           stop
        end if
     else
        !add the new ll
        allocate(THIS%NEXT_LL)
        THIS=>THIS%NEXT_LL
        THIS%LL_ID=ll_ID
        THIS%FIRST_ID=>THIS

        NULLIFY(THIS%NEXT_LL)   
        NULLIFY(THIS%FIRST_TAG) 
        NULLIFY(THIS%NEXT_TAG)  
        NULLIFY(THIS%UP_TAG)    
        NULLIFY(THIS%DOWN_TAG)  
        NULLIFY(THIS%NEXT_ID)
        NULLIFY(THIS%NEXT_PID)
        exit
     end if
  end do
end subroutine libxml2f90_ll_add_list


!================================================================
subroutine libxml2f90_ll_initlist(ll_id)
  use ll_module
  implicit none
  character(*),intent(in)       :: ll_id !the id of the first linklist
  !...........................

  allocate(LL_ROOT)
  THIS=>LL_ROOT
  THIS%LL_ID=ll_ID
  THIS%FIRST_ID=>THIS

  NULLIFY(THIS%NEXT_LL)   
  NULLIFY(THIS%FIRST_TAG) 
  NULLIFY(THIS%NEXT_TAG)  
  NULLIFY(THIS%UP_TAG)    
  NULLIFY(THIS%DOWN_TAG)  
  NULLIFY(THIS%NEXT_ID)
  NULLIFY(THIS%NEXT_PID)
    INITIALIZED=.true.  
end subroutine libxml2f90_ll_initlist


!================================================================
subroutine libxml2f90_ll_inittag()
  use ll_module
  implicit none
  !...........................
  !SELECT THE FIRST_ID BECAUSE THIS ONE'S CONNECTED UP/DOWN
  !THIS=>THIS%FIRST_ID this is done in the calling routine!!!

  allocate(THIS%DOWN_TAG)
  THIS%DOWN_TAG%UP_TAG=>THIS

  THIS=>THIS%DOWN_TAG
  THIS%FIRST_TAG=>THIS
  THIS%FIRST_ID=>THIS
  
  NULLIFY(THIS%NEXT_LL)   
  NULLIFY(THIS%NEXT_TAG)  
  NULLIFY(THIS%DOWN_TAG)  
  NULLIFY(THIS%NEXT_ID)  
  NULLIFY(THIS%NEXT_PID)

  THIS%LL_ID=THIS%UP_TAG%LL_ID
  THIS%TAG='' 

  !JUMP ONE UP, BECAUSE THE OPEN_TAG JUMPS ONE DOWN!!!
  THIS=>THIS%UP_TAG
end subroutine libxml2f90_ll_inittag




!================================================================
subroutine libxml2f90_ll_opentag(tag)
  use ll_module
  implicit none
  character(*),intent(in)       :: tag !
  !...........................
  !check the size of the incoming tag
  if(len(trim(adjustl(tag))).gt.32) then
     print*,'THE TAG CAN NOT BE LARGER THAN 32'
  end if

  !SELECT THE FIRST_ID BECAUSE THIS ONE'S CONNECTED UP/DOWN
  THIS=>THIS%FIRST_ID

  if(.not.associated(THIS%DOWN_TAG)) then
     call libxml2f90_ll_inittag()
!?     return
  end if

  THIS=>THIS%DOWN_TAG
  do
     if(associated(THIS%NEXT_TAG)) then
        THIS=>THIS%NEXT_TAG
        !in XML we can use the same tag more than once!!!
        !we have to be careful while browsing the linklist and descending to lower layers!
        !if(trim(THIS%TAG).eq.trim(TAG)) then
           !print*,'ERROR STOP IN LINKLIST: YOU CAN NOT USE THE SAME TAG TWICE IN THE SAME LAYER: ',TAG
           !stop
        !end if
     else
        !add the new TAG
        allocate(THIS%NEXT_TAG)
        THIS%NEXT_TAG%FIRST_TAG=>THIS%FIRST_TAG
        THIS%NEXT_TAG%LL_ID=THIS%ll_ID
        THIS%NEXT_TAG%UP_TAG=>THIS%UP_TAG
        THIS=>THIS%NEXT_TAG
        THIS%FIRST_ID=>THIS
        THIS%TAG=TAG


        NULLIFY(THIS%NEXT_LL)   
        NULLIFY(THIS%NEXT_TAG)  
        NULLIFY(THIS%DOWN_TAG)  
        NULLIFY(THIS%NEXT_ID)
        NULLIFY(THIS%NEXT_PID)
        exit
     end if
  end do
end subroutine libxml2f90_ll_opentag


subroutine libxml2f90_ll_closetag(tag)
  use ll_module
  implicit none
  character(*),intent(in)       :: tag !
  integer(4)                    :: line !corresponds to the fileline
  !...........................
  !JUMP TO THE FIRST ID BECAUSE IT'S CONNECTED UPWARDS

  THIS=>THIS%FIRST_ID
  if(.not.(trim(THIS%TAG).xmleq.trim(tag))) then
     print*, 'ERROR STOP IN LINKLIST: YOU CAN NOT CLOSE THIS TAG: ',TAG

     !lineinformation
     call libxml2f90_error_getline(line)
     if(line.gt.0) then
        !print line information only if it makes sense
        print *, '(line ',line,')!'
     end if

     stop
  end if
  
  THIS=>THIS%UP_TAG
end subroutine libxml2f90_ll_closetag


!====================================================================  
subroutine libxml2f90_ll_addid(id,size_,value)
  use ll_module
  implicit none
  character(*),intent(in)             :: id !
  integer(4),intent(in)               :: size_
  character(1),intent(in)             :: value(size_)
  integer(4)                          :: line
  !...........................
  if(len(trim(adjustl(id))).gt.32) then
     print*,'THE ID CAN NOT BE LARGER THAN 32'
  end if

  THIS=>THIS%FIRST_ID
  do
     if(associated(THIS%NEXT_ID)) then
        THIS=>THIS%NEXT_ID
        if(trim(THIS%ID).xmleq.trim(ID)) then
           print*,'ERROR STOP IN LINKLIST: YOU CAN NOT USE THE SAME ID TWICE IN THE SAME TAG: ',ID
           !lineinformation
           call libxml2f90_error_getline(line)
           if(line.gt.0) then
              !print line information only if it makes sense
              print *, '(line ',line,')!'
           end if

           stop
        end if
     else
        !add the new ID
        allocate(THIS%NEXT_ID)
        THIS%NEXT_ID%FIRST_ID=>THIS%FIRST_ID
        THIS%NEXT_ID%LL_ID=THIS%ll_ID

        THIS%NEXT_ID%TAG=THIS%TAG
        THIS=>THIS%NEXT_ID
        THIS%ID=ID
        allocate(THIS%VALUE(size_))
        THIS%VALUE(:)=VALUE(:)
  
        NULLIFY(THIS%NEXT_LL)   
        NULLIFY(THIS%FIRST_TAG)
        NULLIFY(THIS%NEXT_TAG)  
        NULLIFY(THIS%UP_TAG)  
        NULLIFY(THIS%DOWN_TAG)  
        NULLIFY(THIS%NEXT_ID)
        NULLIFY(THIS%NEXT_PID)
        exit
     end if
  end do
end subroutine libxml2f90_ll_addid

!====================================================================  
subroutine libxml2f90__ll_edit_id(id,new_id,value)
  implicit none
  character(*),intent(in)             :: id !
  character(*),intent(in)             :: new_id !
  character(*),intent(in)             :: value
  character(1), allocatable, dimension(:) :: valvec
  integer :: l,lenval
  !............................................
  lenval=len(value)
  allocate(valvec(lenval))
  do l=1,lenval
    valvec(l)=value(l:l)
  end do
  !................................
  call libxml2f90_ll_edit_pid(id,new_id,lenval,valvec)
  deallocate(valvec)
  return
end subroutine libxml2f90__ll_edit_id


!====================================================================  
subroutine libxml2f90_ll_edit_id(id,new_id,size_,value)
  !===============================================
  !edit the id or the vale; keep the value for size_=0
  !===============================================
  use ll_module
  implicit none
  character(*),intent(in)             :: id !
  character(*),intent(in)             :: new_id !
  integer(4),intent(in)               :: size_
  character(1),intent(in)             :: value(size_)
  logical(4)                          :: tfound=.false.
  !...........................

  if(id.ne.new_id) then !we change the id -> 
     ! check if the new one's already used
       THIS=>THIS%FIRST_ID
       do
          if(associated(THIS%NEXT_ID)) then
             THIS=>THIS%NEXT_ID
             if(trim(THIS%ID).xmleq.trim(NEW_ID)) then
                print*,'ERROR STOP IN LINKLIST: LL_EDIT: YOU CAN NOT USE&
                     & THE SAME ID TWICE IN THE SAME TAG: ',NEW_ID
                stop
             end if
          else
             exit
          end if
       end do
    end if

  THIS=>THIS%FIRST_ID
  do
     if(associated(THIS%NEXT_ID)) then
        THIS=>THIS%NEXT_ID
        if(trim(THIS%ID).xmleq.trim(ID)) then
           THIS%ID=NEW_ID
           if(size_.ne.0) then
              if(associated(this%value)) nullify(this%value)
              allocate(this%value(size_))
              this%value(:)=''
              this%value=value
              tfound=.true.
              exit
           end if
        end if
     else
        exit
     end if
  end do
  if(.not.tfound) then
     print*,'ERROR STOP IN LINKLIST: EDIT_ID: ID NOT FOUND ',ID
     stop
  end if
end subroutine libxml2f90_ll_edit_id




!====================================================================  
subroutine libxml2f90__ll_edit_pid(id,new_id,value)
  implicit none
  character(*),intent(in)             :: id !
  character(*),intent(in)             :: new_id !
  character(*),intent(in)             :: value
  character(1), allocatable, dimension(:) :: valvec
  integer :: l,lenval
  !............................................
  lenval=len(value)
  allocate(valvec(lenval))
  do l=1,lenval
    valvec(l)=value(l:l)
  end do
  !................................
  call libxml2f90_ll_edit_pid(id,new_id,lenval,valvec)
  deallocate(valvec)
  return


end subroutine libxml2f90__ll_edit_pid


!====================================================================  
subroutine libxml2f90_ll_edit_pid(id,new_id,size_,value)
  !===============================================
  !edit the id or the vale; keep the value for size_=0
  !===============================================
  use ll_module
  implicit none
  character(*),intent(in)             :: id !
  character(*),intent(in)             :: new_id !
  integer(4),intent(in)               :: size_
  character(1),intent(in)             :: value(size_)
  logical(4)                          :: tfound=.false.
  !...........................

  if(id.ne.new_id) then !we change the id -> 
     ! check if the new one's already used
       THIS=>THIS%FIRST_ID
       do
          if(associated(THIS%NEXT_PID)) then
             THIS=>THIS%NEXT_PID
             if(trim(THIS%ID).xmleq.trim(NEW_ID)) then
                print*,'ERROR STOP IN LINKLIST: LL_EDIT: YOU CAN NOT USE&
                     & THE SAME ID TWICE IN THE SAME TAG: ',NEW_ID
                stop
             end if
          else
             exit
          end if
       end do
    end if

  THIS=>THIS%FIRST_ID
  do
     if(associated(THIS%NEXT_PID)) then
        THIS=>THIS%NEXT_PID
        if(trim(THIS%ID).xmleq.trim(ID)) then
           THIS%ID=NEW_ID
           if(size_.ne.0) then
              if(associated(this%value)) nullify(this%value)
              allocate(this%value(size_))
              this%value(:)=''
              this%value=value
              tfound=.true.
              exit
           end if
        else
           exit
        end if
     end if
  end do
  if(.not.tfound) then
    print*,'ERROR STOP IN LINKLIST: EDIT_PID: ID NOT FOUND ',ID
     stop
  end if

end subroutine libxml2f90_ll_edit_pid


!====================================================================  
subroutine libxml2f90_ll_addpureid(size_,value)
  use ll_module
  implicit none
  integer(4),intent(in)               :: size_
  character(1),intent(in)             :: value(size_)
  character(32)                       :: id
  integer(4)                          :: line
  !...........................
  THIS=>THIS%FIRST_ID
  ID=THIS%TAG
  do
     if(associated(THIS%NEXT_ID)) then
        THIS=>THIS%NEXT_ID
        if(trim(THIS%ID).xmleq.trim(ID)) then
           print*,'ERROR STOP IN LINKLIST: PURE XML: YOU CAN NOT USE THE SAME ID TWICE IN THE SAME TAG: ',ID
           !lineinformation
           call libxml2f90_error_getline(line)
           if(line.gt.0) then
              !print line information only if it makes sense
              print *, '(line ',line,')!'
           end if
           
           stop
        end if
     else
        !add the new ID
        allocate(THIS%NEXT_ID)
        THIS%NEXT_ID%FIRST_ID=>THIS%FIRST_ID
        THIS%NEXT_ID%LL_ID=THIS%ll_ID

        THIS%NEXT_ID%TAG=THIS%TAG
        THIS=>THIS%NEXT_ID
        THIS%ID=ID
        allocate(THIS%VALUE(size_))
        THIS%VALUE(:)=VALUE(:)
  
        NULLIFY(THIS%NEXT_LL)   
        NULLIFY(THIS%FIRST_TAG)
        NULLIFY(THIS%NEXT_TAG)  
        NULLIFY(THIS%UP_TAG)  
        NULLIFY(THIS%DOWN_TAG)  
        NULLIFY(THIS%NEXT_ID)
        NULLIFY(THIS%NEXT_PID)
        exit
     end if
  end do
end subroutine libxml2f90_ll_addpureid




!====================================================================  
subroutine libxml2f90_ll_addpid(id,size_,value)
!==================================================
!===== add a 'property' id/value pair to the list.
!===== property means that this pair belongs inside 
!===== a tag 
!==================================================
  use ll_module
  implicit none
  character(*),intent(in)             :: id !
  integer(4),intent(in)               :: size_
  character(1),intent(in)             :: value(size_)
  integer(4)                          :: line
  !...........................
  if(len(trim(adjustl(id))).gt.32) then
     print*,'THE PID CAN NOT BE LARGER THAN 32'
  end if

  THIS=>THIS%FIRST_ID
  do
     if(associated(THIS%NEXT_PID)) then
        THIS=>THIS%NEXT_PID
        if(trim(THIS%ID).xmleq.trim(ID)) then
           print*,'ERROR STOP IN LINKLIST: YOU CAN NOT USE THE SAME PID TWICE IN THE SAME TAG: ',ID
           !lineinformation
           call libxml2f90_error_getline(line)
           if(line.gt.0) then
              !print line information only if it makes sense
              print *, '(line ',line,')!'
           end if

           stop
        end if
     else
        !add the new ID
        allocate(THIS%NEXT_PID)
        THIS%NEXT_PID%FIRST_ID=>THIS%FIRST_ID
        THIS%NEXT_PID%LL_ID=THIS%ll_ID

        THIS=>THIS%NEXT_PID
        THIS%ID=ID
        allocate(THIS%VALUE(size_))
        THIS%VALUE(:)=VALUE(:)
  
        NULLIFY(THIS%NEXT_LL)   
        NULLIFY(THIS%FIRST_TAG)
        NULLIFY(THIS%NEXT_TAG)  
        NULLIFY(THIS%UP_TAG)  
        NULLIFY(THIS%DOWN_TAG)  
        NULLIFY(THIS%NEXT_ID)
        NULLIFY(THIS%NEXT_PID)
        exit
     end if
  end do
end subroutine libxml2f90_ll_addpid



!================================================================
!========   public ll routines


!================================================================
subroutine libxml2f90__ll_add_list(ll_id)
  implicit none
  character(*),intent(in)       :: ll_id !
  !.........................
  if(len(trim(adjustl(ll_id))).gt.32) then
     print*,'THE LISTs NAME CAN NOT BE LARGER THAN 32'
  end if

  call libxml2f90_ll_add_list(ll_id)
end subroutine libxml2f90__ll_add_list

!================================================================
subroutine libxml2f90__ll_opentag(tag)
  implicit none
  character(*),intent(in)           :: tag
  !........................................
  call libxml2f90_ll_opentag(tag)

end subroutine libxml2f90__ll_opentag

!================================================================
subroutine libxml2f90__ll_addid(id,value)
  implicit none
  character(*),intent(in)        :: id
  character(*),intent(in)         :: value
  character(1), allocatable, dimension(:) :: valvec
  integer :: l,lenval
  lenval=len(value)
  allocate(valvec(lenval))
  do l=1,lenval
    valvec(l)=value(l:l)
  end do
  !............................................
  call libxml2f90_ll_addid(id,lenval,valvec)
  deallocate(valvec)
  return
end subroutine libxml2f90__ll_addid

!================================================================
subroutine libxml2f90__ll_addpid(id,value)
  implicit none
  character(*),intent(in)        :: id
  character(*),intent(in)         :: value
  character(1), allocatable, dimension(:) :: valvec
  integer :: l,lenval
  lenval=len(value)
  allocate(valvec(lenval))
  do l=1,lenval
    valvec(l)=value(l:l)
  end do
  !............................................
  call libxml2f90_ll_addpid(id,lenval,valvec)
  deallocate(valvec)
  return
end subroutine libxml2f90__ll_addpid


!================================================================
subroutine libxml2f90__ll_closetag(tag)
  implicit none
  character(*),intent(in)           :: tag
  !........................................
  call libxml2f90_ll_closetag(tag)

end subroutine libxml2f90__ll_closetag


!====================================================================  
subroutine libxml2f90__ll_report(ll_id,nfil,theader)
  use ll_module
  use libxml2f90_module, only: twrite_paw 
  implicit none
  character(*),intent(in)       :: ll_id !the id of the linklist
  integer(4),intent(in)         :: nfil !the file unit we report to
  logical(4),intent(in)         :: theader
  integer(4)                    :: nind  !indentation

  !...........................
  ! 
  if(.NOT.INITIALIZED) then
     print*,'ERROR STOP IN LINKLIST REPORT: THERE IS NO LIST AT ALL!'
     stop
  end if

  !select the right list
  THIS=>LL_ROOT
  do
     if(trim(THIS%LL_ID).eq.trim(LL_ID)) then
        !we found the right list
        exit
     end if
     if(associated(THIS%NEXT_LL)) then
        THIS=>THIS%NEXT_LL
     else
        print*,'ERROR STOP IN LINKLIST REPORT: THERE IS NO LIST CALLED ',LL_ID
        stop
     end if
  end do
  if(theader) then
     write(nfil,*)'<!--'
     write(nfil,*)'=================================================='
     write(nfil,*)'====           LINKLIST REPORT               ====='
     write(nfil,*)'=================================================='
     write(nfil,*)''
     write(nfil,*)'LL_ID: ','    '//'<<<<<<<<>'//trim(LL_ID)//'<>>>>>>>>'
     write(nfil,*)'-->'
     write(nfil,*)''
  end if
  !this already points to the root of the list
  THIS=>THIS%FIRST_ID
  
  !the first tag on the first level
  if(.not.associated(this%down_tag)) then
     !print*,'no tag below'
     return
  else
     thistmp1=>this%down_tag%next_tag
     this=>this%down_tag%next_tag !the first tag is empty by default
  end if
!  nind=-5
  nind=-1
  call libxml2f90_ll_report_rec_wrap(LL_ID,nfil,nind)

  !write a !EOB in the last line for PAW output
  if(twrite_paw) then
     write(nfil,*)'!EOB'
  end if

  this=>thistmp1
end subroutine libxml2f90__ll_report


!====================================================================


subroutine libxml2f90_ll_report_rec_wrap(LL_ID,nfil,nind)
  implicit none
  character(*),intent(in)       :: ll_id !the id of the linklist
  integer(4),intent(in)         :: nfil !the file unit we report to
  integer(4),intent(inout)      :: nind  !indentation

  call libxml2f90_ll_report_rec(LL_ID,nfil,nind)
  return
end subroutine libxml2f90_ll_report_rec_wrap


!====================================================================


subroutine libxml2f90_ll_report_rec(LL_ID,nfil,nind)
  use ll_module
  use libxml2f90_module, only: twrite_paw,indstep 
  implicit none
  character(*),intent(in)       :: ll_id !the id of the linklist
  integer(4),intent(in)         :: nfil !the file unit we report to
  integer(4),intent(inout)      :: nind  !indentation
  character(1)                  :: indch(512)
  integer(4)                    :: pos,sign
  type(ll_type),pointer         :: THISTEMPLOCAL
  character,dimension(46)       :: ctmp_alpha
  ctmp_alpha = (/'a','b','c','g','h','i','j','k','l',&
              &'m','n','o','p','q','r','s','u','v','w','x','y','z','A','B','C','D','G','H','I','J','K','L','M','N',&
              &'O','P','Q','R','S','U','V','W','X','Y','Z','_'/)

  !...........................
  indch(:)=' '

  !increase the indentation
!  nind=nind+5

  !print the opentag including pid's
  !  write(nfil)''
  if(.not.associated(THIS%FIRST_ID%NEXT_PID)) then
     !no pids: write the tag as usual
     if(twrite_paw) then
        write(nfil,*)indch(1:nind),'!',trim(this%tag)
     else
        write(nfil,*)indch(1:nind),'<',trim(this%tag),'>'
     end if
  else !we have pids!
     if(twrite_paw) then
        write(nfil,*)indch(1:nind),'!',trim(this%tag)
     else
        write(nfil,*)indch(1:nind),'<',trim(this%tag)
     end if
     nind=nind+indstep
     !write the pids
     THIS=>THIS%FIRST_ID%NEXT_PID !the first is empty by default
     do
        !write the actual id/value
        if(twrite_paw) then
           !we need ' around the strings
           call libxml2f90__findinchara(size(THIS%VALUE),THIS%VALUE,0,.true.,53,ctmp_alpha,pos,sign)
           if(pos.gt.0.and.(THIS%VALUE(1).ne."'".and.THIS%VALUE(size(THIS%VALUE)).ne."'")) then
              !if we have a string without ''
              write(nfil,*)indch(1:nind),trim(adjustl(THIS%ID)),'=',"'",THIS%VALUE,"'"
           else
              write(nfil,*)indch(1:nind),trim(adjustl(THIS%ID)),'=',THIS%VALUE
           end if
        else
           write(nfil,*)indch(1:nind),trim(adjustl(THIS%ID)),'=',THIS%VALUE
        end if
        !jump to the next id if possible
        if(associated(THIS%NEXT_PID)) then
           THIS=>THIS%NEXT_PID
        else
           exit
        end if
     end do
     !write the >
     if(.not.twrite_paw) then
        write(nfil,*)indch(1:nind),'>'
     end if
     nind=nind-indstep
  end if
  !  write(nfil)''
  nind=nind+indstep+1

  !print the id's
  if(associated(THIS%FIRST_ID%NEXT_ID)) then
     THIS=>THIS%FIRST_ID%NEXT_ID !the first is empty by default

     do
        !write the actual id/value
        if(trim(adjustl(THIS%ID)).eq.trim(adjustl(THIS%TAG))) then
           write(nfil,*)indch(1:nind),THIS%VALUE
        else !write id=value
           if(twrite_paw) then
              !we need ' around the strings
              call libxml2f90__findinchara(size(THIS%VALUE),THIS%VALUE,0,.true.,53,ctmp_alpha,pos,sign)
              if(pos.gt.0.and.(THIS%VALUE(1).ne."'".and.THIS%VALUE(size(THIS%VALUE)).ne."'")) then    
                 write(nfil,*)indch(1:nind),trim(adjustl(THIS%ID)),'=',"'",THIS%VALUE,"'"
              else
                 write(nfil,*)indch(1:nind),trim(adjustl(THIS%ID)),'=',THIS%VALUE
              end if
           else
              !as usual
              write(nfil,*)indch(1:nind),trim(adjustl(THIS%ID)),'=',THIS%VALUE
           end if
        end if
        !jump to the next id if possible
        if(associated(THIS%NEXT_ID)) then
           THIS=>THIS%NEXT_ID
        else
           exit
        end if
     end do
!  nind=nind-2
  end if

  if(.not.associated(this%first_id%down_tag)) then
     !we have no tag below
  else
     thistemplocal=>this
     this=>this%first_id%down_tag%next_tag !the first tag is empty by default
     call libxml2f90_ll_report_rec_wrap(LL_ID,nfil,nind)
     this=>thistemplocal
  end if


  !print closing tag
  nind=nind-(indstep+1)
 ! write(nfil,*)''
  this=>this%first_id
  if(twrite_paw) then
     write(nfil,*)indch(1:nind),'!END'
  else
     write(nfil,*)indch(1:nind),'</',trim(this%tag),'>'
  end if
!  write(nfil,*)''
  !nind=nind-5

  !the next tag on this layer
   if(.not.associated(this%next_tag)) then
     !we have no next tag 
  else
     this=>this%next_tag%first_id
     call libxml2f90_ll_report_rec_wrap(LL_ID,nfil,nind)
  end if

end subroutine libxml2f90_ll_report_rec



!====================================================================  
subroutine libxml2f90__ll_up(tjump)
  use ll_module
  implicit none
  logical(4),intent(out)      :: tjump
  !...........................
  tjump=.false.
  if(.not.associated(this%first_id%up_tag%next_tag))return
  this=>this%first_id%up_tag%next_tag
  tjump=.true.
  return
end subroutine libxml2f90__ll_up


!====================================================================  
subroutine libxml2f90__ll_down(tjump)
  use ll_module
  implicit none
  logical(4),intent(out)      :: tjump
  !...........................
  tjump=.false.
  if(.not.associated(this%first_id%down_tag%next_tag))return
  this=>this%first_id%down_tag%next_tag
  tjump=.true.
  return
end subroutine libxml2f90__ll_down













!**********************************************************************
!**********************************************************************
!***********                  END LINKLIST
!**********************************************************************
!**********************************************************************




!================================================================
!actually not used, but eventually useful
subroutine libxml2f90_parse_find_char(charstart,linestart,charend,lineend,char2find,ichar,line)
  !RETURNS THE ICHAR AND LINE IN STRINGA WHERE CHAR2FIND CAN BE FOUND
  !ICHAR=0 MEANS THAT CHAR2FIND WAS NOT FOUND
  use libxml2f90_module
  implicit none
  integer(4),intent(in)                ::  charstart,linestart,charend,lineend
  character(1),intent(in)              ::  char2find !
  integer(4),intent(out)               ::  ichar,line !the position we found char2find
  !......................................

  ichar=0
  line=0
  do line=linestart, lineend
     if(line.eq.linestart) then
        ichar=scan(stringa(line)(charstart:charend),char2find)
        if(ichar.ne.0) then
           !we found 
           ichar=charstart-1+ichar
           exit
        end if
     else
        ichar=scan(stringa(line)(1:charend),char2find)
        if(ichar.ne.0) then
           !we found 
           exit
        end if
     end if
  end do
  return
end subroutine libxml2f90_parse_find_char



!Note for PID's: we use THIS%FIRST_ID as root and THIS%NEXT_PID to navigate.
!we use id and value as wit ususal id's!




subroutine libxml2f90_transform_paw()
  use libxml2f90_module
  use libxml2f90_strings_module
  implicit none
  character(32),allocatable   :: idarray(:)
  integer(4)                  :: iid
  integer(4)                  :: iread,iwrite,i,c,j
  character(1),allocatable    :: newstringa(:)
  character(3)                :: tmpch
  !..............................................

  allocate(idarray(arraystep))!2000 subsequent id's should be enough ;-)
  idarray(:)=''
  iid=0
  allocate(newstringa(arraystep))
  newstringa(:)='o'
  iwrite=1
  iread=1
  do 
     !cases:
     !exlamation mark
     if(readstring(iread).eq.'!') then
        do i=1,3
           tmpch(i:i)=readstring(iread+i)
        end do
        tmpch=uppercase(tmpch)
        if(tmpch.eq.'EOB'.and.tpaw) then
           !-------------------------------
           !---skip this tag (closes the file in paw files)
           !-------------------------------
           iread=iread+3
 
        else if(tmpch.eq.'END') then
           !-------------------------------
           !---write a closing tag
           !-------------------------------
           
           !check wheter we have an opening tag left
           if(iid.lt.1) then
              print*,'ERROR STOP IN LBXML2F90: TRANSFORM_PAW:'
              print*,'NO TAG LEFT FOR CLOSING !END !'
              stop
           end if

           !---check writestringa's size
           if(size(newstringa).lt.(iwrite+len(trim(adjustl(idarray(iid))))+2)) then
              allocate(tempstringa(size(newstringa)))
              tempstringa(:)=newstringa(:)
              deallocate(newstringa)
              allocate(newstringa(size(tempstringa)+arraystep))
              newstringa(:)=''
              newstringa(1:size(tempstringa))=tempstringa(:)
              deallocate(tempstringa)
           end if

           newstringa(iwrite)='<'
           newstringa(iwrite+1)='/'
           iwrite=iwrite+2
!           allocate(tempstringa())
           do i=1,len(trim(adjustl(idarray(iid))))
              newstringa(iwrite)=idarray(iid)(i:i)
              iwrite=iwrite+1
           end do
           newstringa(iwrite)='>'
           iwrite=iwrite+1
           iid=iid-1
           iread=iread+4
       else
           !-------------------------------
           !---find the opening tag's id
           !-------------------------------
           c=0 !the correction factor
           do i=iread+1,size(readstring)
             if(readstring(i).eq.' '.or.(i.eq.size(readstring)).or.(readstring(i).eq.'!')) then
                 !we have the tag @ iread+1:i-1
                 if(i.eq.size(readstring)) c=-1
!                 if(readstring(i).eq.'!') c=1
                 !we have the tag @ iread+1:i
  !check the c and '!'            

                 !---check writestringa's size
                 if(size(newstringa).lt.(iwrite+i-iread+2)) then
                    allocate(tempstringa(size(newstringa)))
                    tempstringa(:)=newstringa(:)
                    deallocate(newstringa)
                    allocate(newstringa(size(tempstringa)+arraystep))
                    newstringa(:)=''
                    newstringa(1:size(tempstringa))=tempstringa(:)
                    deallocate(tempstringa)
                 end if

                 !write the opening tag
                 newstringa(iwrite)='<'
                 iwrite=iwrite+1
                 newstringa(iwrite:iwrite+i-2-(iread+1)+1-c)=readstring(iread+1:i-2-c)
                 newstringa(iwrite+i-2-(iread+1)+2-c)='>'
                 iwrite=iwrite+i-2-(iread+1)+2-c+1

                 !keep the id for the closing tag
                 idarray(iid+1)=''
                 do j=1,i-2-c-(iread+1)+2
                    idarray(iid+1)(j:j)=readstring(iread+j)
                 end do
                 iid=iid+1
                 iread=i+1-c
                 exit
              end if
              c=0
           end do
           iread=iread-1
      !  if(readstring())
        end if !from '!'
     else if(readstring(iread).eq."'".and.tpaw) then 
        !---------------------------
        !remove "'"
        !---------------------------
        iread=iread+1
     else if(((readstring(iread).eq.">").or.(readstring(iread).eq."<")).and.tpaw) then 
        !---------------------------
        !remove 
        !---------------------------
        iread=iread+1
     else 
        !---------------------------
        !just a normal character
        !---------------------------

        !---check writestringa's size
        if(size(newstringa).lt.(iwrite+1)) then
           allocate(tempstringa(size(newstringa)))
           tempstringa(:)=newstringa(:)
           deallocate(newstringa)
           allocate(newstringa(size(tempstringa)+arraystep))
           newstringa(:)=''
           newstringa(1:size(tempstringa))=tempstringa(:)
           deallocate(tempstringa)
        end if
        
        newstringa(iwrite)=readstring(iread)
        iread=iread+1
        iwrite=iwrite+1
     end if
     !check for exit
     if(iread.ge.size(readstring)) exit
  end do


  !get rid of the unneeded rest
  allocate(tempstringa(iwrite-1))!-1 check this!
  tempstringa(1:iwrite-1)=newstringa(1:iwrite-1)
  deallocate(newstringa)
  deallocate(readstring)
  allocate(readstring(ubound(tempstringa,1)))
  readstring(:)=tempstringa(:)
  deallocate(tempstringa)
  deallocate(idarray)
  return

end subroutine libxml2f90_transform_paw


subroutine libxml2f90__settransform_exm(ttransform_paw_)
  use libxml2f90_module
  implicit none
  logical(4),intent(in)          :: ttransform_paw_
  !.................................................
  ttransform_paw=ttransform_paw_
  return
end subroutine libxml2f90__settransform_exm

subroutine libxml2f90__setwrite_exm(twrite_paw_)
  use libxml2f90_module
  implicit none
  logical(4),intent(in)          :: twrite_paw_
  !.................................................
  twrite_paw=twrite_paw_
  return
end subroutine libxml2f90__setwrite_exm

 subroutine libxml2f90__set_paw(tpaw_)
  use libxml2f90_module
  implicit none
  logical(4),intent(in)          :: tpaw_
  !.................................................
  tpaw=tpaw_
  !paw files are exm structured
  if(tpaw_) call libxml2f90__settransform_exm(.true.)
  return
end subroutine libxml2f90__set_paw

subroutine libxml2f90__set_rmcomma(trmcomma_)
  use libxml2f90_module
  implicit none
  logical(4),intent(in)          :: trmcomma_
  !.................................................
  trmcomma=trmcomma_
  return
end subroutine libxml2f90__set_rmcomma

subroutine libxml2f90__set_rmquotes(trmquotes_)
  use libxml2f90_module
  implicit none
  logical(4),intent(in)          :: trmquotes_
  !.................................................
  trmquotes=trmquotes_
  return
end subroutine libxml2f90__set_rmquotes


subroutine libxml2f90__set_casesensitive(tcasesensitive_)
  use ll_module
  implicit none
  logical(4),intent(in)          :: tcasesensitive_
  !.................................................
  tcasesensitive=tcasesensitive_
  
  !comment:
  !tcasesensitive is kept in ll_module and not in libxml2f90_module 
  !because it is solely used in ll routines
  return
end subroutine libxml2f90__set_casesensitive

subroutine libxml2f90_getline(ipos,line)
  !returns the line ipos belongs to
  !if not found 0
  use libxml2f90_module
  implicit none
  integer(4),intent(in)       :: ipos
  integer(4),intent(out)      :: line
  integer(4)                  :: iline

  line=0
  do iline =1,filelines
     if(lineposa(iline).ge.ipos) then
        line=iline
        exit
     end if
  end do
  return
end subroutine libxml2f90_getline


subroutine libxml2f90_error_getline(line)
  !provides the getline routine for the linklist
  !solved in this way, because it allows to easily comment
  !the call in the linkedlist routines for stand alone usage
  use libxml2f90_module
  implicit none
  integer(4),intent(out)       :: line
  call libxml2f90_getline(lbact,line)
  return
end subroutine libxml2f90_error_getline
