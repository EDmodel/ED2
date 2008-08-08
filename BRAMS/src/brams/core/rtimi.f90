!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

subroutine tend0()

use mem_grid
use mem_tend
use var_tables
use node_mod

implicit none

integer :: n,mxyzp

!     This routine simply sets all tendency arrays to zero.

!     First u,v tendencies

mxyzp = mxp * myp * mzp
call azero(mxyzp,tend%ut(1))
call azero(mxyzp,tend%vt(1))
call azero(mxyzp,tend%wt(1))
call azero(mxyzp,tend%pt(1))

!     Now sclrr tendencies

do n = 1,num_scalar(ngrid)
   call azero(mxyzp,scalar_tab(n,ngrid)%var_t)
enddo

return
end

!**************************************************************************

subroutine hadvance(iac)

use mem_grid
use mem_tend
use mem_basic
use mem_scratch
use node_mod

implicit none

integer :: iac

integer :: mxyzp

!     It is here that the Asselin filter is applied.  For the velocities
!     and pressure, this must be done in two stages, the first when
!     IAC=1 and the second when IAC=2.

mxyzp = mxp * myp * mzp
eps = .2

!     For both IAC=1 and IAC=2, call PREDICT for U, V, W, and P.

call predict(mxyzp,basic_g(ngrid)%uc(1,1,1)   &
   ,basic_g(ngrid)%up(1,1,1),tend%ut(1),scratch%vt3da(1),iac,dtlv)

if (icorflg .eq. 1 .or. jdim .eq. 1) then
   call predict(mxyzp,basic_g(ngrid)%vc(1,1,1)  &
      ,basic_g(ngrid)%vp(1,1,1),tend%vt(1),scratch%vt3da(1),iac,dtlv)
endif

call predict(mxyzp,basic_g(ngrid)%wc(1,1,1),basic_g(ngrid)%wp(1,1,1)  &
   ,tend%wt(1),scratch%vt3da(1),iac,dtlv)
call predict(mxyzp,basic_g(ngrid)%pc(1,1,1),basic_g(ngrid)%pp(1,1,1)  &
   ,tend%pt(1),scratch%vt3da(1),iac,dtlv)

return
end


!**************************************************************************

subroutine predict(npts,ac,ap,fa,af,iac,dtlp)

use mem_grid
use node_mod

implicit none

integer :: npts,iac,m
real :: epsu,dtlp
real, dimension(*) :: ac,ap,fa,af

!     For IAC=3, this routine moves the arrays AC and AP forward by
!     1 time level by adding in the prescribed tendency. It also
!     applies the Asselin filter given by:

!              {AC} = AC + EPS * (AP - 2 * AC + AF)

!     where AP,AC,AF are the past, current and future time levels of A.
!     All IAC=1 does is to perform the {AC} calculation without the AF
!     term present.  IAC=2 completes the calculation of {AC} by adding
!     the AF term only, and advances AC by filling it with input AP
!     values which were already updated in ACOUSTC.
!
epsu = eps
if (ngbegun(ngrid) .eq. 0) epsu = 0.5

if (iac .eq. 1) then
   do m = 1,npts
      ac(m) = ac(m) + epsu * (ap(m) - 2. * ac(m))
   enddo
   return
elseif (iac .eq. 2) then
   do m = 1,npts
      af(m) = ap(m)
      ap(m) = ac(m) + epsu * af(m)
   enddo
elseif (iac .eq. 3) then
   do m = 1,npts
      af(m) = ap(m) + dtlp * fa(m)
   enddo
   if (ngrid .eq. 1 .and. ipara .eq. 0) call cyclic_set(nzp,nxp,nyp,af,'T')
   do m = 1,npts
      ap(m) = ac(m) + epsu * (ap(m) - 2. * ac(m) + af(m))
   enddo
endif

do m = 1,npts
  ac(m) = af(m)
enddo

return
end

!**************************************************************************

subroutine predtr()

use mem_grid
use var_tables
use node_mod

implicit none

integer :: mxyzp,n

!   -  Step thermodynamic variables from  t  to  t+1.
!   -  Set top, lateral and bottom boundary conditions on some variables
!        if needed.
!   -  Call adjustment to assure all positive definite quantities
!        remain positive.
!   -  Rediagnose some thermodynamic quantities for use on the small
!        timestep.

!     Update the scalars and apply lateral, top, and bottom boundary
!     conditions.

mxyzp = mxp * myp * mzp

do n = 1,num_scalar(ngrid)
   call update(mxyzp,scalar_tab(n,ngrid)%var_p  &
                    ,scalar_tab(n,ngrid)%var_t, dtlt)
enddo

return
end






