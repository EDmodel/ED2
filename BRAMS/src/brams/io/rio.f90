!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################


subroutine history_start(name_name)

  ! This routine initializes the model from the history file

  use grid_dims
  use var_tables
  use io_params
  use mem_grid
  use ref_sounding

!srf- bug at history runtype with carma radiation
  USE mem_aerad, ONLY: nwave

  implicit none

  character (len=*) :: name_name

  integer :: ngrids1,ioutput1  &
       ,nnxp1(maxgrds),nnyp1(maxgrds),nnzp1(maxgrds),nzg1,nzs1,npatch1

  integer :: iyr,imn,idy,itm,ie,maxarr,ngr,nc
  character (len=80) :: hnameinh,prefix
  character (len=2) :: cng
  integer, external :: cio_i,cio_f,cio_i_sca,cio_f_sca,cio_f8_sca
  integer,save :: iunhd=11


  ! Open the input history header file and read some of the info.

  nc=len_trim(hfilin)
  hnameinh=hfilin(1:nc-9)//'.vfm'

  call rams_f_open(iunhd,hfilin,'FORMATTED','OLD','READ',0)

  ie=cio_i_sca(iunhd,1,'ngrids',ngrids1,1)
  ngridsh=ngrids1
  ie=cio_i(iunhd,1,'nnxp',nnxp1,ngrids1)
  ie=cio_i(iunhd,1,'nnyp',nnyp1,ngrids1)
  ie=cio_i(iunhd,1,'nnzp',nnzp1,ngrids1)
  ie=cio_i_sca(iunhd,1,'npatch',npatch1,1)
  ie=cio_i_sca(iunhd,1,'nzg',nzg1,1)
  ie=cio_i_sca(iunhd,1,'nzs',nzs1,1)
  ie=cio_i_sca(iunhd,1,'ioutput',ioutput1,1)
  ie=cio_f8_sca(iunhd,1,'time',time,1)

  ! Get the 1-d reference state

  do ngr=1,ngridsh
     write(cng,1) ngr
1    format(i2.2)
     ie=cio_f(iunhd,1,'u01dn'//cng,u01dn(1,ngr),nnzp(ngr))
     ie=cio_f(iunhd,1,'v01dn'//cng,v01dn(1,ngr),nnzp(ngr))
     ie=cio_f(iunhd,1,'pi01dn'//cng,pi01dn(1,ngr),nnzp(ngr))
     ie=cio_f(iunhd,1,'th01dn'//cng,th01dn(1,ngr),nnzp(ngr))
     ie=cio_f(iunhd,1,'dn01dn'//cng,dn01dn(1,ngr),nnzp(ngr))
     ie=cio_f(iunhd,1,'rt01dn'//cng,rt01dn(1,ngr),nnzp(ngr))
  enddo

  ! Put these into regular arrays (for moving grids)
  ie=cio_i(iunhd,1,'ninest',ninest,ngrids1)
  ie=cio_i(iunhd,1,'njnest',njnest,ngrids1)

  ! Check time on file for time requested

  if(dabs(time-timstr).gt..1*dtlong)then
     print*,' !!! History start error                     !!!'
     print*,' !!! Requested time does not match file time !!!'
     print*,' !!! Requested time, file time               !!!'
     print*,' !!! TIMSTR,time,dtlong=',timstr,time,dtlong
     print*,' !!! TIMSTR(m,h,d)=',time/60.,time/3600.,time/86400.
     stop 'HISTORY_START time error'
  endif

  ! Find maximum size of any array on history file. Allocate scratch array of
  ! this size.

  maxarr=0
  do ngr=1,ngridsh
     maxarr=max(maxarr,nnxp1(ngr)*nnyp1(ngr)*nnzp1(ngr)  &
          ,nnxp1(ngr)*nnyp1(ngr)*nzg1*npatch1 &
          ,nnxp1(ngr)*nnyp1(ngr)*nzs1*npatch1 &
	  ,nnxp1(ngr)*nnyp1(ngr)*nwave)
  enddo

  ! read stuff here

  call hist_read(maxarr,hnameinh,iunhd,ioutput1)

  print*,'back from read'
  close(iunhd)


  return
end subroutine history_start

!******************************************************************************

subroutine hist_read(maxarr,hnamein,iunhd,iout1)

  use an_header
  use var_tables

  implicit none

  integer :: maxarr,iunhd,iout1

  include 'interface.h'

  character (len=*) :: hnamein
  integer :: ngr,npts,nptsh,nv,nvh,i
  character type*1,post*10,fmt*3
  real, allocatable :: scr(:)
  integer :: inhunt=10

  type (head_table), allocatable,save :: hr_table(:)

  allocate (scr(maxarr))


  !  Read variable header info

  rewind(iunhd)

  read(iunhd,*) nvbtab
  allocate (hr_table(nvbtab))
  do nv=1,nvbtab
     read(iunhd,*)  hr_table(nv)%string   &
          ,hr_table(nv)%npointer  &
          ,hr_table(nv)%idim_type  &
          ,hr_table(nv)%ngrid  &
          ,hr_table(nv)%nvalues
  enddo


  ! Open history data file

  call rams_f_open(inhunt,hnamein,'UNFORMATTED','OLD','READ',0)

  ! Loop through all variables
  do nvh=1,nvbtab
     ! Read a variable
     nptsh=hr_table(nvh)%nvalues
     read(inhunt)(scr(i),i=1,nptsh)

     !  See if this variable is active in the current run
     ngr=hr_table(nvh)%ngrid
     if(ngr > nvgrids) cycle

     do nv = 1,num_var(ngr)
        npts=vtab_r(nv,ngr)%npts
        if(hr_table(nvh)%string == vtab_r(nv,ngr)%name) then
           if(nptsh /= npts) then
              print*,'Grid point number mismatch on history field:',  &
                   vtab_r(nv,ngr)%name,npts,nptsh
              stop 'History read number points error'
           endif

	   print 33,'History filling grid: ',ngr,nv,vtab_r(nv,ngr)%name,npts
33         format(a25,2i5,3x,a18,i10)

           call atob(npts,scr(1),vtab_r(nv,ngr)%var_p)
           exit
        endif
     enddo

  enddo

  ! Close the input history file

  close(inhunt)

  deallocate(scr,hr_table)

  return
end subroutine hist_read

!******************************************************************************

subroutine hiswrt(restart)

  use an_header
  use var_tables
  use mem_scratch
  use mem_grid
  use io_params

  implicit none

  character*(*) restart

  include 'interface.h'

  ! This routine writes the chosen variables on the history file.

  !character*80 hnamel,hnamelh
  character(len=256) :: hnamel,hnamelh !Changed by ALF

  integer :: nv,nwordh,ngr,nvcnt

  integer, save :: iohunt=10, ncall=0,ncall_head=0,nvtoth=0
  character(len=128), save :: hnameold,hnameoldh

  type (head_table), allocatable,save :: hw_table(:)

  character(len=10) :: c0, c1
  logical :: hereitis

  
  if (ioutput == 0) return

  if (ncall_head == 0) then
     !  Find total number of fields to be written
     do ngr=1,ngrids
        do nv = 1,num_var(ngr)
           if (vtab_r(nv,ngr)%ihist == 1) nvtoth=nvtoth+1
        enddo
     enddo
     allocate (hw_table(nvtoth))
     ncall_head=1
  endif

  ! open a new output file.
  if(restart.eq.'yes') THEN
     call makefnam(hnamel,hfilout,time,iyeara,imontha,idatea,itimea*100  &
          ,'R','$','vfm')
     call makefnam(hnamelh,hfilout,time,iyeara,imontha,idatea,itimea*100  &
          ,'R','head','txt')
  else
     call makefnam(hnamel,hfilout,time,iyeara,imontha,idatea,itimea*100  &
          ,'H','$','vfm')
     call makefnam(hnamelh,hfilout,time,iyeara,imontha,idatea,itimea*100  &
          ,'H','head','txt')
  endif


  call rams_f_open(iohunt,hnamel,'UNFORMATTED','REPLACE','WRITE',iclobber)

  ! Loop through each nest

  nwordh=0
  nvcnt=0
  do ngr=1,ngrids

     !  Loop through the main variable table and write hist flagged variables

     do nv = 1,num_var(ngr)
        if ( vtab_r(nv,ngr)%ihist == 1) then
!!$           print*,'his write:',vtab_r(nv,ngr)%name,vtab_r(nv,ngr)%npts
           call writebin(iohunt,vtab_r(nv,ngr)%var_p,vtab_r(nv,ngr)%npts)
           nwordh=nwordh+vtab_r(nv,ngr)%npts
           nvcnt=nvcnt+1
           hw_table(nvcnt)%string=vtab_r(nv,ngr)%name
           hw_table(nvcnt)%npointer=0
           hw_table(nvcnt)%idim_type=vtab_r(nv,ngr)%idim_type
           hw_table(nvcnt)%ngrid=ngr
           hw_table(nvcnt)%nvalues=vtab_r(nv,ngr)%npts
        endif
     enddo

  enddo

  write(c0,"(f10.1)") time
  write(c1,"(i10)") nwordh
  write(*,"(/,a)") " === History write at Sim time "//trim(adjustl(c0))//" ==="
  write(*,"(a,/)") " === wrote "//trim(adjustl(c1))//" words to file "//&
       &trim(adjustl(hnamel))//" ==="
!!$  print 12,time,nwordh,hnamel
!!$12 format(/,1X,80('*'),/,'  History  write  ','  Time = ',F9.1  &
!!$       ,'      Total words written - ',I9,/,'      File name - ',A60,/,1X,80('*'))

  ! Close history file

  close(iohunt)

  ! Write the COMMON info out to the header file.

  call rams_f_open(iohunt,hnamelh,'FORMATTED','REPLACE','WRITE',iclobber)
  write(iohunt,110) nvcnt
  do nv=1,nvcnt
     write(iohunt,120) hw_table(nv)%string   &
          ,hw_table(nv)%npointer  &
          ,hw_table(nv)%idim_type  &
          ,hw_table(nv)%ngrid  &
          ,hw_table(nv)%nvalues
  enddo

110 format(i6)
120 format(a16,1x,i12,i3,i3,1x,i9)
  call commio('HIST','WRITE',iohunt)
  close(iohunt)

  ! DO NOT remove the old history file if doing a restart or if IFLAG is set

  if(ihistdel == 1) then
     if(ncall == 0) then
        hnameold = hnamel
        hnameoldh = hnamelh
     endif
     if(ncall == 1 .and. iflag == 0) then
        inquire (file=trim(hnameold),exist=hereitis)
        if (hereitis) then
           open (unit=76,file=trim(hnameold))
           close (unit=76,status='delete')
        end if
        inquire (file=trim(hnameoldh),exist=hereitis)
        if (hereitis) then
           open (unit=76,file=trim(hnameoldh))
           close (unit=76,status='delete')
        end if
        hnameold = hnamel
        hnameoldh = hnamelh
     endif
     ncall = 1
  endif

  return
end subroutine hiswrt

!******************************************************************************

subroutine rams_aprep_p (n1,a,b,c)
  implicit none
  integer :: n1
  real :: a(*),b(*),c(*)

  integer :: i

  do i=1,n1
     c(i)=a(i)+b(i)
  enddo

  return
end subroutine rams_aprep_p

!******************************************************************************

subroutine rams_aprep_hkh(n1,hkm,vkh,dn0,scr1,idiffk,xkhkm)
  implicit none
  integer :: n1,idiffk
  real :: xkhkm
  real, dimension(*) :: hkm,vkh,dn0,scr1
  integer :: ind

  if (idiffk <= 3 .or. idiffk == 7) then
     do ind = 1,n1
        scr1(ind) = hkm(ind) * xkhkm / dn0(ind)
     enddo
  elseif (idiffk >= 4 .and. idiffk /= 7) then
     do ind = 1,n1
        scr1(ind) = vkh(ind) / dn0(ind)
     enddo
  endif

  return
end subroutine rams_aprep_hkh

!******************************************************************************

subroutine rams_aprep_vkh(n1,vkh,dn0,vt3dd)
  implicit none
  integer :: n1
  real :: vkh(*),dn0(*),vt3dd(*)
  integer :: ind

  do ind = 1,n1
     vt3dd(ind) = vkh(ind) / dn0(ind)
  enddo

  return
end subroutine rams_aprep_vkh

!******************************************************************************

subroutine anlwrt(restart,vtype)

  use an_header
  use var_tables
  use mem_scratch
  use mem_basic
  use mem_turb
  use mem_grid
  use io_params

  use mem_cuparm, only: nclouds

  ![MLO - For  CARMA
   use mem_aerad, only: nwave
  !MLO]

  implicit none

  include 'interface.h'

  ! This routine writes the chosen variables on the analysis file.

  character*(*) restart,vtype

  character(len=128) :: anamel,anamelh
  character(len=2)  :: cgrid
  character(len=25) :: subaname
  character(len=16) :: varn
  character(len=1)  :: vnam
  character(len=10) :: c0
  logical exans
  integer, save :: ioaunt=10,ncall_head=0,nvtota=0,nvtotl=0  &
       ,nvtot
  integer :: ngr,nv,nvcnt,lenl,npointer,n3d,indwrt,n2d
  real(kind=8) :: timeold

  type (head_table), allocatable,save :: aw_table(:)

  !To avoid rewriting analysis when it is a history restart. It will write if it 
  !   doesn't find all the analysis (header and grid) there, though
  logical :: found
  
  found = .false.

  if (ioutput .eq. 0) return

  if (ncall_head == 0) then
     !  Find total number of fields to be written
     do ngr=1,ngrids
        do nv = 1,num_var(ngr)
           if ( vtab_r(nv,ngr)%ianal == 1) nvtota=nvtota+1
           if ( vtab_r(nv,ngr)%ilite == 1) nvtotl=nvtotl+1
        enddo
     enddo
     nvtot=max(nvtota,nvtotl)
     allocate (aw_table(nvtot))
     ncall_head=1
  endif


  timeold=time
  if(vtype == 'MEAN'.or.vtype == 'BOTH') time=min(time,time-avgtim/2.)


  ! Construct header file name

  if(vtype == 'INST') vnam='A'
  if(vtype == 'LITE') vnam='L'
  if(vtype == 'MEAN') vnam='M'
  if(vtype == 'BOTH') vnam='B'
  call makefnam(anamelh,afilout,time,iyeara,imontha,idatea,  &
       itimea*100,vnam,'head','txt')


! Here is a point that is called just at history start. It checks whether the analysis
! are already there. If they are, it will return without writing anything
  if (restart == 'yes') then
    inquire (file=trim(anamelh),exist=found)
    if (found) then
      gridloop: do ngr=1,ngrids
          write(cgrid,'(a1,i1)') 'g',ngr
          call makefnam(anamel,afilout,time,iyeara,imontha,idatea,  &
                        itimea*100,vnam,cgrid,'vfm')
          inquire(file=trim(anamel),exist=found)
          if (.not.found) exit gridloop
       end do gridloop
    end if
    !So if any of the files were missing, I continue, otherwise, I leave the subroutine
    if (found) return
  end if

  ! Loop through each nest

  nvcnt=0

  do ngr=1,ngrids

     write(cgrid,'(a1,i1)') 'g',ngr
     call makefnam(anamel,afilout,time,iyeara,imontha,idatea,  &
          itimea*100,vnam,cgrid,'vfm')

     lenl = len_trim(anamel)

     inquire(file=anamel,exist=exans)
     if(exans.and.iclobber.eq.0) then
        print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        print*,'!!!   trying to open file name :'
        print*,'!!!       ',anamel
        print*,'!!!   but it already exists. run is ended.'
        print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        stop 'anlwrt'
     endif

     call rams_c_open(anamel(1:lenl)//char(0),'w'//char(0))
     npointer=0

     !  Loop through the main variable table and write those variables
     !     with the correct flag set

     do nv = 1,num_var(ngr)
        !print*,'rio:',vtab_r(nv,ngr)%name,vtab_r(nv,ngr)%ilite

! Writing instantaneous analysis
       if ((vtype == 'INST' .and. vtab_r(nv,ngr)%ianal == 1) .or. &
           (vtype == 'LITE' .and. vtab_r(nv,ngr)%ilite == 1)) then

           varn= vtab_r(nv,ngr)%name
!--- First check whether the variable is a special case -----------------------------!
           if(varn == 'PP') then
              ! Output total Exner function
              call RAMS_aprep_p (nnxyzp(ngr),vtab_r(nv,ngr)%var_p &
                   ,basic_g(ngr)%pi0(1,1,1),scratch%scr1(1) )
              varn='PI'
              !  Rearrange 3-d variables to (i,j,k)
              call rearrange(nnzp(ngr),nnxp(ngr),nnyp(ngr)  &
                   ,scratch%scr1(1),scratch%scr2(1))
           elseif(varn == 'HKM') then
              ! Convert to HKM to HKH (note that VKH is HKH for Deardorff)
              call RAMS_aprep_hkh (nnxyzp(ngr),vtab_r(nv,ngr)%var_p  &
                   ,turb_g(ngr)%vkh(1,1,1),basic_g(ngr)%dn0(1,1,1)  &
                   ,scratch%scr1(1),idiffk(ngr),xkhkm(ngr))
              varn='HKH'
              !  Rearrange 3-d variables to (i,j,k)
              call rearrange(nnzp(ngr),nnxp(ngr),nnyp(ngr)  &
                   ,scratch%scr1(1),scratch%scr2(1))
           elseif(varn == 'VKH') then
              ! Un-density weight VKH
              call RAMS_aprep_vkh (nnxyzp(ngr),vtab_r(nv,ngr)%var_p  &
                   ,basic_g(ngr)%dn0(1,1,1),scratch%scr1(1))
              !  Rearrange 3-d variables to (i,j,k)
              call rearrange(nnzp(ngr),nnxp(ngr),nnyp(ngr)  &
                   ,scratch%scr1(1),scratch%scr2(1))
!---- Now the regular variables ----------------------------------------------------------!
           elseif(vtab_r(nv,ngr)%idim_type == 3) then
              !  Rearrange 3-d variables to (i,j,k)
              call rearrange(nnzp(ngr),nnxp(ngr),nnyp(ngr)  &
                   ,vtab_r(nv,ngr)%var_p,scratch%scr2(1))
           elseif(vtab_r(nv,ngr)%idim_type == 4) then
              !  Rearrange 4-d leaf%soil variables to (i,j,k,ip)
              call rearrange_p(nnxp(ngr),nnyp(ngr),nzg,npatch  &
                   ,vtab_r(nv,ngr)%var_p,scratch%scr2(1))
           elseif(vtab_r(nv,ngr)%idim_type == 5) then
              !  Rearrange 4-d leaf%sfcwater variables to (i,j,k,ip)
              call rearrange_p(nnxp(ngr),nnyp(ngr),nzs,npatch  &
                   ,vtab_r(nv,ngr)%var_p,scratch%scr2(1))
           !  Rearrange 4-d cuparm variables to (i,j,k,icld)
           elseif (vtab_r(nv,ngr)%idim_type == 8) then
              call rearrange_p(nnxp(ngr),nnyp(ngr),nnzp(ngr),nclouds  &
                   ,vtab_r(nv,ngr)%var_p,scratch%scr2(1))
!    For types 2, and 6-9 we don't need to change the order, but I need to copy to         !
!  scr2, so I will use a dum function to do that                                           !
           ! Copy 2-d (i,j,1)
           elseif (vtab_r(nv,ngr)%idim_type == 2) then
              call rearrange_dum(nnxp(ngr),nnyp(ngr),1       &
                   ,vtab_r(nv,ngr)%var_p,scratch%scr2(1))
           ! Copy 3-d (i,j,ip)
           elseif (vtab_r(nv,ngr)%idim_type == 6) then
              call rearrange_dum(nnxp(ngr),nnyp(ngr),npatch  &
                   ,vtab_r(nv,ngr)%var_p,scratch%scr2(1))
           ! Copy 3-d (i,j,kwave)
           elseif (vtab_r(nv,ngr)%idim_type == 7) then
              call rearrange_dum(nnxp(ngr),nnyp(ngr),nwave  &
                   ,vtab_r(nv,ngr)%var_p,scratch%scr2(1))
           ! Copy 3-d (i,j,ip)
           elseif (vtab_r(nv,ngr)%idim_type == 9) then
              call rearrange_dum(nnxp(ngr),nnyp(ngr),nclouds  &
                   ,vtab_r(nv,ngr)%var_p,scratch%scr2(1))
           end if
           nvcnt=nvcnt+1
           aw_table(nvcnt)%string=varn
           aw_table(nvcnt)%npointer=npointer
           aw_table(nvcnt)%idim_type=vtab_r(nv,ngr)%idim_type
           aw_table(nvcnt)%ngrid=ngr
           aw_table(nvcnt)%nvalues=vtab_r(nv,ngr)%npts
           call vforecr(ioaunt,scratch%scr2(1),vtab_r(nv,ngr)%npts  &
                ,18,scratch%scr1(1),scratch%scr1(1),'LIN',npointer)
       elseif((vtype == 'MEAN' .and. vtab_r(nv,ngr)%ianal == 1) .or. &
              (vtype == 'BOTH' .and. vtab_r(nv,ngr)%ilite == 1)) then

           varn= vtab_r(nv,ngr)%name

!--- First check whether the variable is a special case -----------------------------!
           if(varn == 'PP') then
              ! Output total Exner function
              call RAMS_aprep_p (nnxyzp(ngr),vtab_r(nv,ngr)%var_m &
                   ,basic_g(ngr)%pi0(1,1,1),scratch%scr1(1) )
              varn='PI'
              !  Rearrange 3-d variables to (i,j,k)
              call rearrange(nnzp(ngr),nnxp(ngr),nnyp(ngr)  &
                   ,scratch%scr1(1),scratch%scr2(1))
           elseif(varn == 'HKM') then
              ! Convert to HKM to HKH (note that VKH is HKH for Deardorff)
              call RAMS_aprep_hkh (nnxyzp(ngr),vtab_r(nv,ngr)%var_m  &
                   ,turb_g(ngr)%vkh(1,1,1),basic_g(ngr)%dn0(1,1,1)  &
                   ,scratch%scr1(1),idiffk(ngr),xkhkm(ngr))
              varn='HKH'
              !  Rearrange 3-d variables to (i,j,k)
              call rearrange(nnzp(ngr),nnxp(ngr),nnyp(ngr)  &
                   ,scratch%scr1(1),scratch%scr2(1))
           elseif(varn == 'VKH') then
              ! Un-density weight VKH
              call RAMS_aprep_vkh (nnxyzp(ngr),vtab_r(nv,ngr)%var_m  &
                   ,basic_g(ngr)%dn0(1,1,1),scratch%scr1(1))
              !  Rearrange 3-d variables to (i,j,k)
              call rearrange(nnzp(ngr),nnxp(ngr),nnyp(ngr)  &
                   ,scratch%scr1(1),scratch%scr2(1))
!---- Now the regular variables -----------------------------------------------------------!
           elseif(vtab_r(nv,ngr)%idim_type == 3) then
              !  Rearrange 3-d variables to (i,j,k)
              call rearrange(nnzp(ngr),nnxp(ngr),nnyp(ngr)  &
                   ,vtab_r(nv,ngr)%var_m,scratch%scr2(1))
           elseif(vtab_r(nv,ngr)%idim_type == 4) then
              !  Rearrange 4-d leaf%soil variables to (i,j,k,ip)
              call rearrange_p(nnxp(ngr),nnyp(ngr),nzg,npatch  &
                   ,vtab_r(nv,ngr)%var_m,scratch%scr2(1))
           elseif(vtab_r(nv,ngr)%idim_type == 5) then
              !  Rearrange 4-d leaf%sfcwater variables to (i,j,k,ip)
              call rearrange_p(nnxp(ngr),nnyp(ngr),nzs,npatch  &
                   ,vtab_r(nv,ngr)%var_m,scratch%scr2(1))
              !  Rearrange 4-d cuparm% variables to (i,j,k,ip)
           elseif(vtab_r(nv,ngr)%idim_type == 8) then
              call rearrange_p(nnxp(ngr),nnyp(ngr),nnzp(ngr),nclouds  &
                   ,vtab_r(nv,ngr)%var_m,scratch%scr2(1))
!    For types 2, 6, 7, and 9 we don't need to change the order, but I need to copy to     !
!  scr2, so I will use a dum function to do that                                           !
           ! Copy 2-d (i,j,1)
           elseif (vtab_r(nv,ngr)%idim_type == 2) then
              call rearrange_dum(nnxp(ngr),nnyp(ngr),1       &
                   ,vtab_r(nv,ngr)%var_m,scratch%scr2(1))
           ! Copy 3-d (i,j,ip)
           elseif (vtab_r(nv,ngr)%idim_type == 6) then
              call rearrange_dum(nnxp(ngr),nnyp(ngr),npatch  &
                   ,vtab_r(nv,ngr)%var_m,scratch%scr2(1))
           ! Copy 3-d (i,j,kwave)
           elseif (vtab_r(nv,ngr)%idim_type == 7) then
              call rearrange_dum(nnxp(ngr),nnyp(ngr),nwave  &
                   ,vtab_r(nv,ngr)%var_m,scratch%scr2(1))
           ! Copy 3-d (i,j,kwave)
           elseif (vtab_r(nv,ngr)%idim_type == 9) then
              call rearrange_dum(nnxp(ngr),nnyp(ngr),nclouds  &
                   ,vtab_r(nv,ngr)%var_m,scratch%scr2(1))
           end if
           nvcnt=nvcnt+1
           aw_table(nvcnt)%string=varn
           aw_table(nvcnt)%npointer=npointer
           aw_table(nvcnt)%idim_type=vtab_r(nv,ngr)%idim_type
           aw_table(nvcnt)%ngrid=ngr
           aw_table(nvcnt)%nvalues=vtab_r(nv,ngr)%npts
           call vforecr(ioaunt,scratch%scr2(1),vtab_r(nv,ngr)%npts  &
                ,18,scratch%scr1(1),scratch%scr1(1),'LIN',npointer)
       end if



     enddo

     call rams_c_close()
     close(ioaunt)

  enddo

  ! Write the header information out to the file.

  call rams_f_open(ioaunt,anamelh,  &
       'FORMATTED','REPLACE','WRITE',iclobber)

  write(ioaunt,110) nvcnt
  do nv=1,nvcnt
     write(ioaunt,120) aw_table(nv)%string   &
          ,aw_table(nv)%npointer  &
          ,aw_table(nv)%idim_type  &
          ,aw_table(nv)%ngrid  &
          ,aw_table(nv)%nvalues
  enddo

110 format(i6)
120 format(a16,1x,i12,i3,i3,1x,i9)

  call commio('ANAL','WRITE',ioaunt)
  close(ioaunt)
  select case (trim(vtype))
  case('LITE')
     subaname='  Analysis lite write'
  case('MEAN')
     subaname='  Averaged analysis write    '
  case('BOTH')
     subaname='  Averaged analysis lite write    '
  case default
     subaname='  Analysis write         '
  end select

  write(c0,"(f10.1)") time
  write(*,"(/,a)") " === "//trim(adjustl(subaname))//" at Sim time "//trim(adjustl(c0))//" ==="
  write(*,"(a,/)") " === wrote file "//&
       &trim(adjustl(anamel))//" ==="

!!$  print 12,subaname,time,anamel
!!$12 format(/,1X,79('*'),/,  &
!!$       A35,'  Time = ',F9.0  &
!!$       ,/,'      Header file name - ',A60  &
!!$       ,/,1X,79('*'))

  ! Reset the time back to the original
  if(vtype == 'MEAN'.or.vtype == 'BOTH')  time=timeold

  return
end subroutine anlwrt


!******************************************************************************

subroutine rearrange_dum(n2,n3,n4,a,b)
  implicit none
  integer, intent(in)                       :: n2,n3,n4
  real,    intent(in) , dimension(n2,n3,n4) :: a
  real,    intent(out), dimension(n2,n3,n4) :: b
  integer                                   :: i,j,k
  do i=1,n2
    do j=1,n3
      do k=1,n4
        b(i,j,k)=a(i,j,k)
      end do
    end do
  end do
  return
end subroutine rearrange_dum

!******************************************************************************

subroutine rearrange_dum4(n2,n3,n4,n5,a,b)
  implicit none
  integer, intent(in)                          :: n2,n3,n4,n5
  real,    intent(in) , dimension(n2,n3,n4,n5) :: a
  real,    intent(out), dimension(n2,n3,n4,n5) :: b
  integer                                      :: i,j,k,l
  do i=1,n2
    do j=1,n3
      do k=1,n4
        do l=1,n5
          b(i,j,k,l)=a(i,j,k,l)
        end do
      end do
    end do
  end do
  return
end subroutine rearrange_dum4

!******************************************************************************

subroutine rearrange_p(n2,n3,n4,n5,a,b)
  implicit none

  integer :: n2,n3,n4,n5
  real :: a(n4,n2,n3,n5),b(n2,n3,n4,n5)

  integer :: i,j,k,ip

  do ip = 1,n5
     do k = 1,n4
        do j = 1,n3
           do i = 1,n2
              b(i,j,k,ip) = a(k,i,j,ip)
           enddo
        enddo
     enddo
  enddo
  return
end subroutine rearrange_p

!******************************************************************************

subroutine writebin(iun,var,npts)
  implicit none
  real :: var(*)
  integer :: iun,npts,i

  write(iun) (var(i),i=1,npts)

  return
end subroutine writebin
