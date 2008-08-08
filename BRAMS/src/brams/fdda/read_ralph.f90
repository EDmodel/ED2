!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

subroutine rr_upa_ver (ifile)

! Reads Ralph sfc file version and header

use obs_input

implicit none

integer :: ifile

integer :: imarker,nvar,nh,nt,nw

header(ifile)%head_string(1:max_head_vars)=''
header(ifile)%therm_string(1:max_upa_vars)=''
header(ifile)%wind_string(1:max_upa_vars)=''

read(header(ifile)%iun,*) imarker
rewind header(ifile)%iun

if(imarker.eq.999999) then
   read(header(ifile)%iun,*) imarker,header(ifile)%iver
else
   header(ifile)%iver=1
endif

if(header(ifile)%iver==3) then

   read(header(ifile)%iun,*) header(ifile)%nhead
   do nh=1,header(ifile)%nhead
      read(header(ifile)%iun,*) header(ifile)%head_string(nh)
   enddo
   read(header(ifile)%iun,*) header(ifile)%nvtherm
   do nt=1,header(ifile)%nvtherm
      read(header(ifile)%iun,*) header(ifile)%therm_string(nt)
   enddo
   read(header(ifile)%iun,*) header(ifile)%nvwind
   do nw=1,header(ifile)%nvwind
      read(header(ifile)%iun,*) header(ifile)%wind_string(nw)
   enddo
endif

return
end

!***************************************************************************

subroutine rr_upa_obs (ifile,qcheck,ierr)

! Reads one upper air obs from Ralph upper air file
!   need to detect if there is an additional header

use obs_input

implicit none

integer :: ifile,ierr
character(len=*) :: qcheck

integer :: ntok,iqfl(5)
character(len=80) :: var_string
character(len=16) :: cflags,tokens(100)
character(len=256) :: line
integer :: jd,n,i,j,iq,iqf,kl

ierr=0

if(header(ifile)%iver==1) then

   read(header(ifile)%iun,'(a)',end=20,err=20) line
   if(len_trim(line)==0) goto 10
   call parse(line,tokens,ntok)
   read(tokens(1),*) rupa_obs%jdate
   read(tokens(2),*) rupa_obs%jtime
   read(tokens(3),'(a)') rupa_obs%id
   read(tokens(4),*) rupa_obs%lp
   read(tokens(5),*) rupa_obs%lz
   read(tokens(6),*) rupa_obs%lat
   read(tokens(7),*) rupa_obs%lon
   read(tokens(8),*) rupa_obs%elev

   if(rupa_obs%lp.gt.max_up_levs.or. rupa_obs%lz.gt.max_up_levs) then
      print*,'Rawindsonde read error!'
      print*,'  Number of input levels greater than max_up_levs'
      print*,'  upa_lp,upa_lz,max_up_levs:',rupa_obs%lp  &
               ,rupa_obs%lz, max_up_levs
      stop 'rawin-max_up_levs'
   endif

   do n=1,rupa_obs%lp
      read(header(ifile)%iun,80,end=20,err=20)rupa_obs%p(n),rupa_obs%z(n)  &
                   ,rupa_obs%t(n),rupa_obs%r(n)  &
                 ,((rupa_obs%iqflagsp(n,i,j),j=1,3),i=1,4)
   enddo
   80 format(2f12.3,f10.2,f10.4,2x,4(':',3i1))

   do n=1,rupa_obs%lz
      read(header(ifile)%iun,85,end=20,err=20)rupa_obs%zz(n),rupa_obs%fz(n)  &
                   ,rupa_obs%dz(n)  &
                 ,((rupa_obs%iqflagsz(n,i,j),j=1,3),i=1,3)
   enddo
   85 format(f12.3,2f10.2,2x,3(':',3i1))

elseif(header(ifile)%iver==2.or.header(ifile)%iver==3) then

   read(header(ifile)%iun,'(a)',end=20,err=20) line
   if(len_trim(line)==0) goto 10
   call parse(line,tokens,ntok)
   read(tokens(1),*) rupa_obs%jyear
   read(tokens(2),*) rupa_obs%jmonth
   read(tokens(3),*) rupa_obs%jdate
   read(tokens(4),*) rupa_obs%jtime
   read(tokens(5),'(a)') rupa_obs%id
   read(tokens(6),*) rupa_obs%lp
   read(tokens(7),*) rupa_obs%lz
   read(tokens(8),*) rupa_obs%lat
   read(tokens(9),*) rupa_obs%lon
   read(tokens(10),*) rupa_obs%elev

   IF(rupa_obs%lp > max_up_levs.or.  &
      rupa_obs%lz > max_up_levs) then
      print*,'Rawindsonde read error!'
      print*,'  Number of input levels greater than max_up_levs'
      print*,'  upa_lp,upa_lz,max_up_levs:',rupa_obs%lp,rupa_obs%lz,max_up_levs
      stop 'rawin-max_up_levs'
   endif

!   rupa_obs%jdate=rupa_obs%jyear*10000+rupa_obs%jmonth*100+rupa_obs%jdate

   do n=1,rupa_obs%lp
      read(header(ifile)%iun,*,end=20,err=20)   &
                 rupa_obs%p(n),iqfl(1),rupa_obs%z(n),iqfl(2)  &
                ,rupa_obs%t(n),iqfl(3),rupa_obs%r(n),iqfl(4)
      do i=1,4
         rupa_obs%iqflagsp(n,i,1)=iqfl(i)/100
         rupa_obs%iqflagsp(n,i,2)=mod(iqfl(i)/10,10)
         rupa_obs%iqflagsp(n,i,3)=mod(iqfl(i),10)
      enddo
   enddo

   do n=1,rupa_obs%lz
       read(header(ifile)%iun,*,end=20,err=20)  &
                 rupa_obs%zz(n),iqfl(1),rupa_obs%fz(n),iqfl(2)  &
                ,rupa_obs%dz(n),iqfl(3)
       do i=1,3
          rupa_obs%iqflagsz(n,i,1)=iqfl(i)/100
          rupa_obs%iqflagsz(n,i,2)=mod(iqfl(i)/10,10)
          rupa_obs%iqflagsz(n,i,3)=mod(iqfl(i),10)
       enddo
   enddo

endif

! If desired, check QC flags and set appropriate values to missing.

if(qcheck=='yes') then

   do kl=1,rupa_obs%lp
      iqf=1
      do iq=1,3
         if((rupa_obs%iqflagsp(kl,1,iq).ne.5 .and.  &
             rupa_obs%iqflagsp(kl,1,iq).ne.0 .and.  &
             rupa_obs%iqflagsp(kl,1,iq).ne.8) .or.  &
             rupa_obs%p(kl).lt.-998.) iqf=0
      enddo
      if(iqf.eq.0) rupa_obs%p(kl)=-999.
      
      iqf=1
      do iq=1,3
         if((rupa_obs%iqflagsp(kl,2,iq).ne.5 .and.  &
             rupa_obs%iqflagsp(kl,2,iq).ne.0 .and.  &
             rupa_obs%iqflagsp(kl,2,iq).ne.8) .or.  &
             rupa_obs%z(kl).lt.-998.) iqf=0
      enddo
      if(iqf.eq.0) rupa_obs%z(kl)=-999.

      iqf=1
      do iq=1,3
         if((rupa_obs%iqflagsp(kl,3,iq).ne.5 .and.  &
             rupa_obs%iqflagsp(kl,3,iq).ne.0 .and.  &
             rupa_obs%iqflagsp(kl,3,iq).ne.8) .or.  &
             rupa_obs%t(kl).lt.-998.) iqf=0
      enddo
      if(iqf.eq.0) rupa_obs%t(kl)=-999.

      iqf=1
      do iq=1,3
         if((rupa_obs%iqflagsp(kl,4,iq).ne.5 .and.  &
             rupa_obs%iqflagsp(kl,4,iq).ne.0 .and.  &
             rupa_obs%iqflagsp(kl,4,iq).ne.8) .or.  &
             rupa_obs%r(kl).lt.-998.) iqf=0
      enddo
      if(iqf.eq.0) rupa_obs%r(kl)=-999.

   enddo

   do kl=1,rupa_obs%lz
      iqf=1
      do iq=1,3
         if((rupa_obs%iqflagsz(kl,1,iq).ne.5 .and.  &
             rupa_obs%iqflagsz(kl,1,iq).ne.0 .and.  &
             rupa_obs%iqflagsz(kl,1,iq).ne.8) .or.  &
             rupa_obs%zz(kl).lt.-998.) iqf=0
      enddo
      if(iqf.eq.0) rupa_obs%zz(kl)=-999.

      iqf=1
      do iq=1,3
         if((rupa_obs%iqflagsz(kl,2,iq).ne.5 .and.  &
             rupa_obs%iqflagsz(kl,2,iq).ne.0 .and.  &
             rupa_obs%iqflagsz(kl,2,iq).ne.8) .or.  &
             rupa_obs%fz(kl).lt.-998.) iqf=0
      enddo
      if(iqf.eq.0) rupa_obs%fz(kl)=-999.

      iqf=1
      do iq=1,3
         if((rupa_obs%iqflagsz(kl,3,iq).ne.5 .and.  &
             rupa_obs%iqflagsz(kl,3,iq).ne.0 .and.  &
             rupa_obs%iqflagsz(kl,3,iq).ne.8) .or.  &
             rupa_obs%dz(kl).lt.-998.) iqf=0
      enddo
      if(iqf.eq.0) rupa_obs%dz(kl)=-999.
      
   enddo
   
endif

10 continue

return

20 continue
ierr=1

return
end

!***************************************************************************

subroutine rr_sfc_ver (ifile)

! Reads Ralph sfc file version and header

use obs_input

implicit none

integer :: ifile

integer :: imarker,nvar,nh,ns

header(ifile)%head_string(1:max_head_vars)=''
header(ifile)%sfc_string(1:max_sfc_vars)=''

read(header(ifile)%iun,*) imarker
rewind header(ifile)%iun

if(imarker.eq.999999) then
   read(header(ifile)%iun,*) imarker,header(ifile)%iver
else
   header(ifile)%iver=1
endif

if(header(ifile)%iver==3) then
   read(header(ifile)%iun,*) header(ifile)%nhead
   do nh=1,header(ifile)%nhead
      read(header(ifile)%iun,*) header(ifile)%head_string(nh)
   enddo
endif
read(header(ifile)%iun,*) header(ifile)%nvsfc
do ns=1,header(ifile)%nvsfc
   read(header(ifile)%iun,*) header(ifile)%sfc_string(ns)
enddo
if(header(ifile)%iver.eq.1) read(header(ifile)%iun,*)

return
end

!***************************************************************************

subroutine rr_sfc_obs(ifile,qcheck,ierr)

! Reads one surface obs from Ralph sfc file
!   need to detect if there is an additional header

use obs_input

implicit none

integer :: ifile,ierr
character(len=*) :: qcheck

integer :: ntok,iqfl(5)
character(len=80) :: var_string
character(len=16) :: cflags,tokens(100)
character(len=256) :: line
integer :: jd,nvar,nv,ic,k,iq,iqf,nh,nt

ierr=0
nvar=5
tokens(1:100)=''

if(header(ifile)%iver==1) then

   read(header(ifile)%iun,'(a)',end=20,err=20) line
   if(len_trim(line)==0) goto 10
   call parse(line,tokens,ntok)
   read(tokens(1),*,err=1,end=1) jd
   read(tokens(2),*,err=1,end=1) rsfc_obs%jtime
   read(tokens(3),'(a)',err=1,end=1) rsfc_obs%id
   read(tokens(4),*,err=1,end=1) rsfc_obs%lat
   read(tokens(5),*,err=1,end=1) rsfc_obs%lon
   read(tokens(6),*,err=1,end=1) rsfc_obs%elev
   rsfc_obs%hgt=0.
   rsfc_obs%ihgtflg=1
   read(tokens(7),*,err=1,end=1) rsfc_obs%ff
   read(tokens(8),*,err=1,end=1) rsfc_obs%dd
   read(tokens(9),*,err=1,end=1) rsfc_obs%t
   read(tokens(10),*,err=1,end=1) rsfc_obs%td
   read(tokens(11),*,err=1,end=1) rsfc_obs%p
   read(tokens(12),'(a)',err=1,end=1) cflags
   1 continue

   ic=1
   do nv=1,nvar
      read(cflags(ic:),'(1x,3i1)') (rsfc_obs%iqflags(nv,k),k=1,3)
      ic=ic+4
   enddo
   
   rsfc_obs%jyear=jd/10000
   rsfc_obs%jmonth=mod(jd,10000)/100
   rsfc_obs%jdate=mod(jd,100)

elseif(header(ifile)%iver==2) then

   read(header(ifile)%iun,'(a)',end=20,err=20) line
   if(len_trim(line)==0) goto 10
   call parse(line,tokens,ntok)
   read(tokens(1),*,err=2,end=2) rsfc_obs%jyear
   read(tokens(2),*,err=2,end=2) rsfc_obs%jmonth
   read(tokens(3),*,err=2,end=2) rsfc_obs%jdate
   read(tokens(4),*,err=2,end=2) rsfc_obs%jtime
   read(tokens(5),'(a)',err=2,end=2) rsfc_obs%id
   read(tokens(6),*,err=2,end=2) rsfc_obs%lat
   read(tokens(7),*,err=2,end=2) rsfc_obs%lon
   read(tokens(8),*,err=2,end=2) rsfc_obs%elev
   rsfc_obs%hgt=0.
   rsfc_obs%ihgtflg=1
   read(tokens(9),*,err=2,end=2) rsfc_obs%ff
   read(tokens(10),*,err=2,end=2) iqfl(1)
   read(tokens(11),*,err=2,end=2) rsfc_obs%dd
   read(tokens(12),*,err=2,end=2) iqfl(2)
   read(tokens(13),*,err=2,end=2) rsfc_obs%t
   read(tokens(14),*,err=2,end=2) iqfl(3)
   read(tokens(15),*,err=2,end=2) rsfc_obs%td
   read(tokens(16),*,err=2,end=2) iqfl(4)
   read(tokens(17),*,err=2,end=2) rsfc_obs%p
   read(tokens(18),*,err=2,end=2) iqfl(5)
   2 continue
   
   do nv=1,nvar
      rsfc_obs%iqflags(nv,1)=iqfl(nv)/100
      rsfc_obs%iqflags(nv,2)=mod(iqfl(nv)/10,10)
      rsfc_obs%iqflags(nv,3)=mod(iqfl(nv),10)
   enddo

elseif(header(ifile)%iver==3) then

   read(header(ifile)%iun,'(a)',end=20,err=20) line
   if(len_trim(line)==0) goto 10
   call parse(line,tokens,ntok)
   
   rsfc_obs%hgt=0.
   rsfc_obs%ihgtflg=1
   nt=1
   do nh=1,header(ifile)%nhead
      if(header(ifile)%head_string(nh)=='YEAR') then
         read(tokens(nt),*,err=3) rsfc_obs%jyear
         nt=nt+1
      elseif(header(ifile)%head_string(nh)=='MONTH') then
         read(tokens(nt),*,err=3) rsfc_obs%jmonth
         nt=nt+1
      elseif(header(ifile)%head_string(nh)=='DAY') then
         read(tokens(nt),*,err=3) rsfc_obs%jdate
         nt=nt+1
      elseif(header(ifile)%head_string(nh)=='HOUR') then
         read(tokens(nt),*,err=3) rsfc_obs%jtime
         nt=nt+1
      elseif(header(ifile)%head_string(nh)=='STATION_ID') then
         read(tokens(nt),'(a)',err=3) rsfc_obs%id
         nt=nt+1
      elseif(header(ifile)%head_string(nh)=='LATITUDE') then
         read(tokens(nt),*,err=3) rsfc_obs%lat
         nt=nt+1
      elseif(header(ifile)%head_string(nh)=='LONGITUDE') then
         read(tokens(nt),*,err=3) rsfc_obs%lon
         nt=nt+1
      elseif(header(ifile)%head_string(nh)=='ELEVATION') then
         read(tokens(nt),*,err=3) rsfc_obs%elev
         nt=nt+1
      elseif(header(ifile)%head_string(nh)=='HEIGHT') then
         read(tokens(nt),*,err=3) rsfc_obs%hgt
         nt=nt+1
      elseif(header(ifile)%head_string(nh)=='HEIGHT_FLAG') then
         read(tokens(nt),*,err=3) rsfc_obs%ihgtflg
         nt=nt+1
      else
         print*,'unknown header record'
      endif
   enddo
   
   read(tokens(nt),*,err=3,end=3) rsfc_obs%ff
   read(tokens(nt+1),*,err=3,end=3) iqfl(1)
   read(tokens(nt+2),*,err=3,end=3) rsfc_obs%dd
   read(tokens(nt+3),*,err=3,end=3) iqfl(2)
   read(tokens(nt+4),*,err=3,end=3) rsfc_obs%t
   read(tokens(nt+5),*,err=3,end=3) iqfl(3)
   read(tokens(nt+6),*,err=3,end=3) rsfc_obs%td
   read(tokens(nt+7),*,err=3,end=3) iqfl(4)
   read(tokens(nt+8),*,err=3,end=3) rsfc_obs%p
   read(tokens(nt+9),*,err=3,end=3) iqfl(5)
   3 continue
   
   do nv=1,nvar
      rsfc_obs%iqflags(nv,1)=iqfl(nv)/100
      rsfc_obs%iqflags(nv,2)=mod(iqfl(nv)/10,10)
      rsfc_obs%iqflags(nv,3)=mod(iqfl(nv),10)
   enddo

endif

! If desired, check QC flags and set appropriate values to missing.

if(qcheck=='yes') then
   iqf=1
   do iq=1,3
      if((rsfc_obs%iqflags(1,iq).ne.5 .and.  &
          rsfc_obs%iqflags(1,iq).ne.0 .and.  &
          rsfc_obs%iqflags(1,iq).ne.8) .or. &
          rsfc_obs%ff.lt.-998.) iqf=0
   enddo
   if(iqf.eq.0) rsfc_obs%ff=-999.

   iqf=1
   do iq=1,3
      if((rsfc_obs%iqflags(2,iq).ne.5 .and.  &
          rsfc_obs%iqflags(2,iq).ne.0 .and.  &
          rsfc_obs%iqflags(2,iq).ne.8) .or.  &
          rsfc_obs%dd.lt.-998.) iqf=0
   enddo
   if(iqf.eq.0) rsfc_obs%dd=-999.

   iqf=1
   do iq=1,3
      if((rsfc_obs%iqflags(3,iq).ne.5 .and.  &
          rsfc_obs%iqflags(3,iq).ne.0 .and.  &
          rsfc_obs%iqflags(3,iq).ne.8) .or.  &
          rsfc_obs%t.lt.-998.) iqf=0
   enddo
   if(iqf.eq.0) rsfc_obs%t=-999.

   iqf=1
   do iq=1,3
      if((rsfc_obs%iqflags(4,iq).ne.5 .and.  &
          rsfc_obs%iqflags(4,iq).ne.0 .and.  &
          rsfc_obs%iqflags(4,iq).ne.8) .or.  &
          rsfc_obs%td.lt.-998.) iqf=0
   enddo
   if(iqf.eq.0) rsfc_obs%td=-999.

   iqf=1
   do iq=1,3
      if((rsfc_obs%iqflags(5,iq).ne.5 .and.  &
          rsfc_obs%iqflags(5,iq).ne.0 .and.  &
          rsfc_obs%iqflags(5,iq).ne.8) .or.  &
          rsfc_obs%p.lt.-998.) iqf=0
   enddo
   if(iqf.eq.0) rsfc_obs%p=-999.

endif

10 continue

return

20 continue
ierr=1

return
end


!***************************************************************************

subroutine sfc_data_convert (varn,cvars,nvars)

use obs_input

implicit none

integer :: nvars
real :: varn(nvars)
character(len=*) :: cvars(nvars)

character(len=16) :: cvar
real :: vv
integer :: ll,nv
real,external :: rs

! Convert units and type of sfc input data

do nv=1,nvars
   cvar=cvars(nv)
   ll=len_trim(cvar)
   varn(nv)=-999.
   if(cvar(1:ll)=='ue') then
      ! earth-relative u in m/s
      if(rsfc_obs%dd>-998..and.rsfc_obs%ff>-998.)  &
         call winduv(rsfc_obs%dd,rsfc_obs%ff,varn(nv),vv)
   elseif(cvar(1:ll)=='ve') then  
      ! earth-relative u in m/s
      if(rsfc_obs%dd>-998..and.rsfc_obs%ff>-998.)  &
         call winduv(rsfc_obs%dd,rsfc_obs%ff,vv,varn(nv))
   elseif(cvar(1:ll)=='speed') then  
      ! wind speed in m/s
      if(rsfc_obs%ff>-998.) varn(nv)=rsfc_obs%ff
   elseif(cvar(1:ll)=='direction') then  
      ! wind direction
      if(rsfc_obs%dd>-998.) varn(nv)=rsfc_obs%dd
   elseif(cvar(1:ll)=='tempc') then  
      ! temperature in C
      if(rsfc_obs%t>-998.) varn(nv)=rsfc_obs%t
   elseif(cvar(1:ll)=='tempf') then  
      ! temperature in F
      if(rsfc_obs%t>-998.)  varn(nv)=rsfc_obs%t*1.8+32.
   elseif(cvar(1:ll)=='dewptc') then  
      ! dewpoint in C
      if(rsfc_obs%t>-998.) varn(nv)=rsfc_obs%td
   elseif(cvar(1:ll)=='dewptf') then  
      ! dewpoint in F
      if(rsfc_obs%t>-998.) varn(nv)=rsfc_obs%td*1.8+32.
   elseif(cvar(1:ll)=='press') then  
      ! pressure in mb
      if(rsfc_obs%p>-998.) varn(nv)=rsfc_obs%p*.01
   elseif(cvar(1:ll)=='relhum') then  
      ! rh in percent
      if(rsfc_obs%t>-998. .and. rsfc_obs%td>-998. .and.  &
         rsfc_obs%p>-998.) &
      varn(nv)=100.*min(1.,max(0.  &
               ,rs(rsfc_obs%p,rsfc_obs%td+273.16)  &
               /rs(rsfc_obs%p,rsfc_obs%t+273.16)))
   else
      print*,'UNKNOWN CONVERT VARIABLE in sfc_data_convert !!!!',cvar
      stop 'sfc_data_convert'
   endif
enddo


return
end

!***************************************************************************

subroutine upa_get_profile (varn,nlevels,cvar,ctype)

use obs_input
use rconstants

implicit none

integer :: nlevels
real :: varn(*)
character(len=*) :: cvar,ctype

real :: vv
integer :: ll,nv,k,nlev
real,external :: rs,td

! Convert units and type of sfc input data

   if(ctype(1:1)=='z') nlevels=rupa_obs%lz
   if(ctype(1:1)=='p') nlevels=rupa_obs%lp
   ll=len_trim(cvar)
   do k=1,nlevels
      varn(k)=-999.

      if(cvar(1:ll)=='ue') then
         ! earth-relative u in m/s
         if(rupa_obs%dz(k) > -998..and.rupa_obs%fz(k)>-998.)  &
            call winduv(rupa_obs%dz(k),rupa_obs%fz(k),varn(k),vv)
      elseif(cvar(1:ll)=='ve') then  
         ! earth-relative v in m/s
         if(rupa_obs%dz(k) > -998..and.rupa_obs%fz(k)>-998.)  &
            call winduv(rupa_obs%dz(k),rupa_obs%fz(k),vv,varn(k))
      elseif(cvar(1:ll)=='zz') then  
         ! height of wind levels in m
         if(rupa_obs%zz(k) > -998.) varn(k)=rupa_obs%zz(k)
      elseif(cvar(1:ll)=='speed') then  
         ! wind speed in m/s
         if(rupa_obs%fz(k) > -998.) varn(k)=rupa_obs%fz(k)
      elseif(cvar(1:ll)=='direction') then  
         ! wind direction
         if(rupa_obs%dz(k) > -998.) varn(k)=rupa_obs%dz(k)
      elseif(cvar(1:ll)=='tempc') then  
         ! temperature in C
         if(rupa_obs%t(k)>-998.) varn(k)=rupa_obs%t(k)
      elseif(cvar(1:ll)=='tempf') then  
         ! temperature in F
         if(rupa_obs%t(k)>-998.) varn(k)=rupa_obs%t(k)*1.8+32.
      elseif(cvar(1:ll)=='theta') then  
         ! theta in K
         if(rupa_obs%p(k)>-998..and.rupa_obs%t(k)>-998.) &
            varn(k)=(rupa_obs%t(k)+273.16)*(p00/rupa_obs%p(k))**rocp
      elseif(cvar(1:ll)=='p_pas') then  
         ! pressure in pascals
         if(rupa_obs%p(k)>-998.) varn(k)=rupa_obs%p(k)
      elseif(cvar(1:ll)=='pi') then  
         ! Exner function
         if(rupa_obs%p(k)>-998.) varn(k)=cp*(rupa_obs%p(k)*p00i)**rocp
      elseif(cvar(1:ll)=='dewptc') then  
         ! dewpoint in C
         if(rupa_obs%r(k)>-998..and.rupa_obs%t(k)>-998..and.  &
            rupa_obs%p(k)>-998.) &
            vv=min(1.,max(0.,rupa_obs%r(k)))  &
                 *rs(rupa_obs%p(k),rupa_obs%t(k)+273.16)
            varn(k)=td(rupa_obs%p(k),vv )-273.16
      elseif(cvar(1:ll)=='dewptf') then  
         ! dewpoint in F
         if(rupa_obs%r(k)>-998..and.rupa_obs%t(k)>-998..and.  &
                rupa_obs%p(k)>-998.) &
            vv=min(1.,max(0.,rupa_obs%r(k)))  &
                 *rs(rupa_obs%p(k),rupa_obs%t(k)+273.16)
            varn(k)=(td(rupa_obs%p(k),vv )-273.16)*1.8+32.
      elseif(cvar(1:ll)=='geo') then  
         ! geopotential in m
         if(rupa_obs%z(k)>-998.) varn(k)=rupa_obs%z(k)
      elseif(cvar(1:ll)=='mixrat') then  
         ! vapor in kg/kg
         if(rupa_obs%r(k)>-998..and.rupa_obs%t(k)>-998..and.  &
            rupa_obs%p(k)>-998.) &
            varn(k)=min(1.,max(0.,rupa_obs%r(k)))  &
                 *rs(rupa_obs%p(k),rupa_obs%t(k)+273.16)
       elseif(cvar(1:ll)=='relhum') then  
         ! rh in percent
         if(rupa_obs%r(k)>-998.) varn(k)=100.  &
             *min(1.,max(0.,rupa_obs%r(k)))
      elseif(cvar(1:ll)=='relfrac') then  
         ! rh in fraction
         if(rupa_obs%r(k)>-998.) varn(k)=  &
                                 min(1.,max(0.,rupa_obs%r(k)))
      else
         print*,'UNKNOWN CONVERT VARIABLE in upa_get_profile !!!!',cvar
         stop 'upa_get_profile'
      endif
   enddo

return
end

!***************************************************************************

subroutine ralph_vars (ivok,nvar,cvar,nvr,namev,unitv,cfmt)

implicit none

integer :: ivok,nvar,nvr
character(len=*) :: cvar,namev,unitv,cfmt

if(cvar.eq.'ue') then
   nvr=nvar
   namev='U_WIND_COMPONET'
   unitv='m/s'
   cfmt='2x,f9.2,2x,i4.3'
   
elseif(cvar.eq.'ve') then
   nvr=nvar
   namev='V_WIND_COMPONET'
   unitv='m/s'
   cfmt='2x,f9.2,2x,i4.3'
   
elseif(cvar.eq.'speed') then
   nvr=nvar
   namev='WINDSPEED'
   unitv='m/s'
   cfmt='2x,f9.2,2x,i4.3'
   
elseif(cvar.eq.'direction_e') then
   nvr=nvar
   namev='WIND_DIRECTION'
   unitv='deg'
   cfmt='2x,f7.0,2x,i4.3'
   
elseif(cvar.eq.'tempc') then
   nvr=nvar
   namev='TEMPERATURE'
   unitv='C'
   cfmt='2x,f7.1,2x,i4.3'
   
elseif(cvar.eq.'tempf') then
   nvr=nvar
   namev='TEMPERATURE'
   unitv='F'
   cfmt='2x,f7.1,2x,i4.3'
   
elseif(cvar.eq.'dewptc') then
   nvr=nvar
   namev='DEWPOINT'
   unitv='C'
   cfmt='2x,f7.1,2x,i4.3'
   
elseif(cvar.eq.'dewptf') then
   nvr=nvar
   namev='DEWPOINT'
   unitv='F'
   cfmt='2x,f7.1,2x,i4.3'
   
elseif(cvar.eq.'dewptk') then
   nvr=nvar
   namev='DEWPOINT'
   unitv='K'
   cfmt='2x,f7.1,2x,i4.3'
   
elseif(cvar.eq.'relhum') then
   nvr=nvar
   namev='RELATIVE HUMIDITY'
   unitv='pct'
   cfmt='2x,f7.0,2x,i4.3'
   
elseif(cvar.eq.'press') then
   nvr=nvar
   namev='STN_PRES'
   unitv='Pa'
   cfmt='2x,f10.1,2x,i4.3'
   
elseif(cvar.eq.'sea_press') then
   nvr=nvar
   namev='SLP'
   unitv='Pa'
   cfmt='2x,f10.1,2x,i4.3'
   
!elseif(cvar.eq.'24hr_precip') then
!   nvr=nvar
!   namev='6-HR_PCP'
!   unitv='mm'
!   cfmt='2x,f7.1,2x,i4.3'

!elseif(cvar.eq.'6hr_precip') then
!   nvr=nvar
!   namev='24-HR_PCP'
!   unitv='mm'
!   cfmt='2x,f7.1,2x,i4.3'

!elseif(cvar.eq.'snow_depth_ps') then
!   nvr=nvar
!   namev='SNOW_DEPTH'
!   unitv='m'
!   cfmt='2x,f7.1,2x,i4.3'

elseif(cvar.eq.'cloud_frac') then
   nvr=nvar
   namev='CLOUD_COVER'
   unitv='fraction'
   cfmt='2x,f6.3,2x,i4.3'
   
else
   print*,'remove variable from namelist:',cvar
   stop 'ralph_vars: variable not in Ralph file variable list'
endif

ivok=1

return
end
