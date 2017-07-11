!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

subroutine isenio (inout,iun,n1,n2)

use isan_coms

implicit none

integer :: iun,n1,n2
character(len=*) :: inout

integer :: npts,nlt,nx3,ny3,ninn,l

if(inout.eq.'IN') THEN

   read(iun,920) iyy,imm,idd,ihh,nx3,ny3,ninn,(levth(l),l=1,ninn)
   920 format(7i4,(13i6))
   if(nx3.ne.n1.or.ny3.ne.n2.or.ninn.ne.nisn) then
      print*,'Isentropic stage grid dimensions do not match'
      print*,'   configuration file on read !'
      print*,' File dimens - ',nx3,ny3,ninn
      print*,' Run  dimens - ',n1,n2,nisn
      stop 'IO3-2'
   endif

   npts=n1*n2
   do nlt=1,nisn
      call vfirec(iun,pi_u(1,1,nlt),npts,'LIN')
      call vmissr(pi_u(1,1,nlt),npts,1e30,-998.)
      call vfirec(iun,pi_v(1,1,nlt),npts,'LIN')
      call vmissr(pi_v(1,1,nlt),npts,1e30,-998.)
      call vfirec(iun,pi_s(1,1,nlt),npts,'LIN')
      call vmissr(pi_p(1,1,nlt),npts,1e30,-.5)
      call vfirec(iun,pi_p(1,1,nlt),npts,'LIN')
      call vmissr(pi_s(1,1,nlt),npts,1e30,-.5)
      call vfirec(iun,pi_r(1,1,nlt),npts,'LIN')
      call vmissr(pi_r(1,1,nlt),npts,1e30,-.5)
   enddo

   call vfirec(iun,rs_u,npts,'LIN')
   call vmissr(rs_u,npts,1e30,-998.)
   call vfirec(iun,rs_v,npts,'LIN')
   call vmissr(rs_v,npts,1e30,-998.)
   call vfirec(iun,rs_p,npts,'LIN')
   call vmissr(rs_p,npts,1e30,-.5)
   call vfirec(iun,rs_t,npts,'LIN')
   call vmissr(rs_t,npts,1e30,-.5)
   call vfirec(iun,rs_r,npts,'LIN')
   call vmissr(rs_r,npts,1e30,-.5)
   call vfirec(iun,rs_s,npts,'LIN')
   call vmissr(rs_s,npts,1e30,-.5)
   call vfirec(iun,rs_top,npts,'LIN')
   call vmissr(rs_top,npts,1e30,-.5)
   call vfirec(iun,rs_qual,npts,'LIN')
   call vmissr(rs_qual,npts,1e30,-.5)
   
   call vfirec(iun,rs_slp,npts,'LIN')
   call vmissr(rs_slp,npts,1e30,-.5)
   call vfirec(iun,rs_sfp,npts,'LIN')
   call vmissr(rs_sfp,npts,1e30,-.5)
   call vfirec(iun,rs_sft,npts,'LIN')
   call vmissr(rs_sft,npts,1e30,-.5)
   call vfirec(iun,rs_snow,npts,'LIN')
   call vmissr(rs_snow,npts,1e30,-.5)
   call vfirec(iun,rs_sst,npts,'LIN')
   call vmissr(rs_sst,npts,1e30,-.5)

   print 201,' *****  Isentropic file input *****************'  &
        ,iyear,imonth,idate,ihour,n1,n2,nisn  &
        ,(levth(l),l=1,nisn)
   201 format(//,a,//  &
        ,' *',7X,' Date (year,month,day,hour)  - ',4I5,/  &
        ,' *',7X,' Number of X,Y points        - ',2I5,/  &
        ,' *',7X,' Number of isentropic levels - ',I5,/  &
        ,' *',7X,' Isentropic levels (K)       - '/,(32X,8I5))
   print '(a)',' **********************************************'

endif

if(inout.eq.'out') then

   write(iun,920) iyear,imonth,idate,ihour,n1,n2,nisn  &
        ,(levth(l),l=1,nisn)

   npts=n1*n2
   do nlt=1,nisn
      call vmissw(pi_u(1,1,nlt),npts,pi_scra,1E30,-999.)
      call vforec(iun,pi_scra,npts,18,pi_scrb,'LIN')
      call vmissw(pi_v(1,1,nlt),npts,pi_scra,1E30,-999.)
      call vforec(iun,pi_scra,npts,18,pi_scrb,'LIN')
      call vmissw(pi_p(1,1,nlt),npts,pi_scra,1E30,-1.)
      call vforec(iun,pi_scra,npts,18,pi_scrb,'LIN')
      call vmissw(pi_s(1,1,nlt),npts,pi_scra,1E30,-1.)
      call vforec(iun,pi_scra,npts,18,pi_scrb,'LIN')
      call vmissw(pi_r(1,1,nlt),npts,pi_scra,1E30,-1.)
      call vforec(iun,pi_scra,npts,18,pi_scrb,'LIN')
   ENDDO

   call vmissw(rs_u,npts,pi_scra,1E30,-999.)
   call vforec(iun,pi_scra,npts,18,pi_scrb,'LIN')
   call vmissw(rs_v,npts,pi_scra,1E30,-999.)
   call vforec(iun,pi_scra,npts,18,pi_scrb,'LIN')
   call vmissw(rs_p,npts,pi_scra,1E30,-1.)
   call vforec(iun,pi_scra,npts,18,pi_scrb,'LIN')
   call vmissw(rs_t,npts,pi_scra,1E30,-1.)
   call vforec(iun,pi_scra,npts,18,pi_scrb,'LIN')
   call vmissw(rs_r,npts,pi_scra,1E30,-1.)
   call vforec(iun,pi_scra,npts,18,pi_scrb,'LIN')
   call vmissw(rs_s,npts,pi_scra,1E30,-1.)
   call vforec(iun,pi_scra,npts,18,pi_scrb,'LIN')
   call vmissw(rs_top,npts,pi_scra,1E30,-1.)
   call vforec(iun,pi_scra,npts,18,pi_scrb,'LIN')
   call vmissw(rs_qual,npts,pi_scra,1E30,-1.)
   call vforec(iun,pi_scra,npts,18,pi_scrb,'LIN')
   
   call vmissw(rs_slp,npts,pi_scra,1E30,-1.)
   call vforec(iun,pi_scra,npts,18,pi_scrb,'LIN')
   call vmissw(rs_sfp,npts,pi_scra,1E30,-1.)
   call vforec(iun,pi_scra,npts,18,pi_scrb,'LIN')
   call vmissw(rs_sft,npts,pi_scra,1E30,-1.)
   call vforec(iun,pi_scra,npts,18,pi_scrb,'LIN')
   call vmissw(rs_snow,npts,pi_scra,1E30,-1.)
   call vforec(iun,pi_scra,npts,18,pi_scrb,'LIN')
   call vmissw(rs_sst,npts,pi_scra,1E30,-1.)
   call vforec(iun,pi_scra,npts,18,pi_scrb,'LIN')

   print 201,' *****  Isentropic file written *************'  &
        ,iyear,imonth,idate,ihour,n1,n2,nisn  &
        ,(levth(l),l=1,nisn)

   print 303,igridfl,gobsep,gobrad
   303 format(/,  &
         ' Grid flag (IGRIDFL)               -',I4,/  &
        ,' Grid-obs separation in degrees    -',F5.2,/  &
        ,' Grid-obs radius influence degrees -',F5.2)

endif

return
end

!***************************************************************************

subroutine sigzio (inout,iun,n1,n2)

use isan_coms

implicit none

integer :: iun,n1,n2
character(len=*) :: inout

integer :: npts,nlt,l,ninn,nx3,ny3

if(inout.eq.'IN') then
   read(iun,920) iyy,imm,idd,ihh,nx3,ny3,ninn  &
        ,(sigz(l),l=1,ninn)
   920 format(7i4,(9f8.2))
   if(nx3.ne.n1.or.ny3.ne.n2.or.ninn.ne.nsigz)then
      print*,'Sigma-z grid dimensions do not match'
      print*,'   input data on read !'
      print*,' File  dimensions - ',nx3,ny3,ninn
      print*,' Input dimensions - ',n1,n2,nsigz
      stop 'iO3-2'
   endif

   npts=n1*n2
   do nlt=1,nsigz
      call vfirec(iun,ps_u(1,1,nlt),npts,'LIN')
      call vmissr(ps_u(1,1,nlt),npts,1e30,-998.)
      call vfirec(iun,ps_v(1,1,nlt),npts,'LIN')
      call vmissr(ps_v(1,1,nlt),npts,1e30,-998.)
      call vfirec(iun,ps_p(1,1,nlt),npts,'LIN')
      call vmissr(ps_p(1,1,nlt),npts,1e30,-.5)
      call vfirec(iun,ps_t(1,1,nlt),npts,'LIN')
      call vmissr(ps_t(1,1,nlt),npts,1e30,-.5)
      call vfirec(iun,ps_r(1,1,nlt),npts,'LIN')
      call vmissr(ps_r(1,1,nlt),npts,1e30,-.5)
   enddo

   print 201,' *****  Sigma-z file input *****************'  &
        ,iyear,imonth,idate,ihour,n1,n2,nsigz  &
        ,(sigz(l),l=1,nsigz)
   201 format(//,a,//  &
        ,' *',7X,' Date (year,month,day,hour)  - ',4I5,/  &
        ,' *',7X,' Number of X,Y points        - ',2I5,/  &
        ,' *',7X,' Number of sigma-z levels    - ',I5,/  &
        ,' *',7X,' Sigma-z levels (m)          - '/,(32X,7F8.1))
   print '(a)',' **********************************************'

endif

if(inout.eq.'OUT') then
   write(iun,920) iyear,imonth,idate,ihour,n1,n2,nsigz  &
        ,(sigz(l),l=1,nsigz)

   npts=n1*n2
   do nlt=1,nsigz
      call vmissw(ps_u(1,1,nlt),npts,ps_scra,1E30,-999.)
      call vforec(iun,ps_scra,npts,18,ps_scrb,'LIN')
      call vmissw(ps_v(1,1,nlt),npts,ps_scra,1E30,-999.)
      call vforec(iun,ps_scra,npts,18,ps_scrb,'LIN')
      call vmissw(ps_p(1,1,nlt),npts,ps_scra,1E30,-1.)
      call vforec(iun,ps_scra,npts,18,ps_scrb,'LIN')
      call vmissw(ps_t(1,1,nlt),npts,ps_scra,1E30,-1.)
      call vforec(iun,ps_scra,npts,18,ps_scrb,'LIN')
      call vmissw(ps_r(1,1,nlt),npts,ps_scra,1E30,-1.)
      call vforec(iun,ps_scra,npts,18,ps_scrb,'LIN')
   enddo

   print 201,' *****  Sigma-z file written *************'  &
        ,iyear,imonth,idate,ihour,n1,n2,nsigz   &
        ,(sigz(l),l=1,nsigz)

endif

return
end

!***************************************************************************

subroutine vmissw (af,n,as,fm,fx)
implicit none

integer :: n
real :: af(*),as(*),fm,fx

integer :: i

do i=1,n
   as(i)=af(i)
   if(af(i).ge.fm) as(i)=fx
enddo

return
end

!***************************************************************************

subroutine vmissr (af,n,fm,fx)
implicit none

integer :: n
real :: af(*),fm,fx

integer :: i

do i=1,n
   if(af(i).le.fx) af(i)=fm
enddo

return
end

