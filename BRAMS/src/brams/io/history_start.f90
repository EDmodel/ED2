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
  integer, external :: cio_i,cio_f,cio_f8,cio_i_sca,cio_f_sca,cio_f8_sca
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
