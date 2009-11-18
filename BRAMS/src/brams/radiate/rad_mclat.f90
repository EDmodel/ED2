!===============================================================================
! OLAM version 2.12  

! Copyright (C) 2002-2006; All Rights Reserved; 
! Duke University, Durham, North Carolina, USA 

! Portions of this software are copied or derived from the RAMS software
! package.  The following copyright notice pertains to RAMS and its derivatives,
! including OLAM:  

   !----------------------------------------------------------------------------
   ! Copyright (C) 1991-2006  ; All Rights Reserved ; Colorado State University; 
   ! Colorado State University Research Foundation ; ATMET, LLC 

   ! This software is free software; you can redistribute it and/or modify it 
   ! under the terms of the GNU General Public License as published by the Free
   ! Software Foundation; either version 2 of the License, or (at your option)
   ! any later version. 

   ! This software is distributed in the hope that it will be useful, but
   ! WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
   ! or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
   ! for more details.
 
   ! You should have received a copy of the GNU General Public License along
   ! with this program; if not, write to the Free Software Foundation, Inc.,
   ! 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA 
   ! (http://www.gnu.org/licenses/gpl.html) 
   !----------------------------------------------------------------------------

! It is requested that any scientific publications based on application of OLAM
! include the following acknowledgment:  "OLAM was developed at the 
! Edmund T. Pratt Jr. School of Engineering, Duke University."

! For additional information, including published references, please contact
! the software authors, Robert L. Walko (robert.walko@duke.edu)
! or Roni Avissar (avissar@duke.edu).
!===============================================================================
!MLO - Adapted to ED-BRAMS
subroutine rad_mclat(m1,nrad,koff,glat,rtgt)

  use mem_mclat,   only: sslat, mclat, ypp_mclat, mcol
  use mem_radiate, only: nadd_rad, zmrad
  use rconstants,  only: gordry
  use mem_grid,    only: zm
  use harr_coms,   only: dl,pl,rl,co2l,tl,o3l,zml,ztl,dzl

  implicit none
  integer, intent(in) :: m1
  integer, intent(in) :: nrad     ! # of vertical radiation levels
  integer, intent(in) :: koff     ! offset between model and radiation levels
  real, intent(in)    :: glat
  real, intent(in)    :: rtgt

  integer :: lv, lf, k, kadd
  real    :: deltaz, tavg

  ! Before this subroutine is called, McClatchey soundings have been interpolated
  ! in time between summer and winter values, and spline coefficients for latitudinal
  ! interpolation have been pre-computed.
  !
  ! The following call to spline2 completes the latitudinal spline interpolation 
  ! of a complete vertical column for a specific latitude.
  ! The result is in mcol(lv,lf).

  do lv = 1,33    ! Loop over number of vertical levels
     do lf = 1,6  ! Loop over number of data types

        if (lf == 2) cycle                                    ! no need for pressure
        if (nadd_rad == 0 .and. lf /= 5 .and. lf /= 1) cycle  ! only o3 needed if nadd_rad=0

        call spline2(13,sslat,mclat(1:13,lv,lf),ypp_mclat(1:13,lv,lf),glat,mcol(lv,lf))
     end do
  end do

  ! Model values of dl, pl, tl, rl, zml, and ztl were filled in harr_raddriv
  ! from k = 1 to k = mza - 1 - koff

  if (nadd_rad > 0) then

     if (zmrad < zml(m1-koff-1)) then
        print*, 'Error - top of radiation grid is below the model grid'
        stop    'in rad_mclat'
     endif

     ! Compute heights of added levels for this column.

     deltaz = (zmrad - zml(m1-koff-1)) / real(nadd_rad)

     do k = m1-koff,nrad
        zml(k) = zml(k-1) + deltaz
        ztl(k) = .5 * (zml(k) + zml(k-1))
     enddo

  endif

  ! Compute dzl values.

  do k = 2,nrad
     dzl(k) = zml(k) - zml(k-1)
  enddo
  dzl(1) = zml(2)

  ! Interpolate O3 from Mclatchy sounding to all levels in radiation column,

  call hintrp_cc(33, mcol(1:33,5), mcol(1:33,1), nrad, o3l(1:nrad), ztl(1:nrad))

  if (nadd_rad > 0) then

     ! Interpolate other variables (temperature, density, vapor mixing ratio)
     ! to added levels.

     kadd = m1 - koff
     call hintrp_cc(33, mcol(1:33,3), mcol(1:33,1), nadd_rad, tl(kadd:nrad), ztl(kadd:nrad))
     call hintrp_cc(33, mcol(1:33,4), mcol(1:33,1), nadd_rad, rl(kadd:nrad), ztl(kadd:nrad))
     call hintrp_cc(33, mcol(1:33,6), mcol(1:33,1), nadd_rad, dl(kadd:nrad), ztl(kadd:nrad))

     ! Compute pressure of added levels by hydrostatic integration.
     
     do k = kadd, nrad
        tavg = 0.5 * (tl(k)+tl(k-1))
        pl(k) = pl(k-1) * exp( -gordry * (ztl(k) - ztl(k-1))*rtgt / tavg )
     enddo
     
     !----- Copy the CO2 at the top level to the remainder levels. ------------------------!
     do k = kadd,nrad
        co2l(k) = co2l(kadd-1) * dl(kadd-1) / dl(k)
     end do
  end if

  return
end subroutine rad_mclat

!----------------------------------------------------------------------------------------

subroutine spline1(n,xdat,ydat,yppdat)
  !MLO Imported from OLAM 2.5.1 
  implicit none
  integer , intent(in)                  :: n
  real    , intent(in)   , dimension(n) :: xdat,ydat
  real    , intent(inout), dimension(n) :: yppdat
  !Local variables
  integer :: i,j
  real    :: fx,fy
  real, allocatable, dimension(:) :: scr
  
  allocate(scr(n))
  yppdat(1) = 0.
  yppdat(n) = 0.
  scr(1)    = 0.
  
  do i = 2,n-1
     fx = (xdat(i)-xdat(i-1)) / (xdat(i+1)-xdat(i-1))
     fy = fx * yppdat(i-1) + 2.
     yppdat(i) = (fx - 1.) / fy
     scr(i) = (6.* ((ydat(i+1)-ydat(i)) / (xdat(i+1)-xdat(i))   &
                   -(ydat(i)-ydat(i-1)) / (xdat(i)-xdat(i-1)))  &
                   /(xdat(i+1)-xdat(i-1)) - fx * scr(i-1)) / fy
  end do

  do j=n-1,2,-1
    yppdat(j) = yppdat(j) * yppdat(j+1) + scr(j)
  end do
  
  deallocate(scr)
  
  return
end subroutine spline1


!----------------------------------------------------------------------------------------
subroutine spline2(n,xdat,ydat,yppdat,x,y)
  !MLO Imported from OLAM 2.5.1 
  implicit none
  integer, intent(in)                 :: n
  real   , intent(in)  , dimension(n) :: xdat,ydat,yppdat
  real   , intent(in)                 :: x
  real   , intent(out)                :: y
  
  integer :: j,jhi,jlo
  real    :: a,b,h
  integer :: i

  jlo = 1
  jhi = n
  i=0
  do while (jhi-jlo > 1 .and. i <= n)
     i=i+1
     j = (jhi+jlo) / 2
     if (xdat(j) > x) then
       jhi = j
     else
       jlo = j
     end if
  end do
  if (i > n) then
     write(*,*) '!!!! Infinit loop at spline2 (rad_mclat.f90)'
     write(*,*) 'n=',n,' x=',x,' y=',y
     write(*,*)  'xdat=',xdat
     stop '!!!!! Infinit loop spline2 (rad_mclat.f90)'
  end if
  h = xdat(jhi) - xdat(jlo)
  if (h == 0) stop '»»» Bad xdat input in spline2'
  a = (xdat(jhi)-x) / h
  b = (x-xdat(jlo)) / h
  y = a * ydat(jlo) + b * ydat(jhi) &
    + ((a**3-a) * yppdat(jlo) + (b**3-b) * yppdat(jhi)) * h**2 / 6.

  return
end subroutine spline2

!----------------------------------------------------------------------------------------

subroutine hintrp_cc(na,vctra,eleva,nb,vctrb,elevb)
  !MLO - Imported from OLAM 2.5.1
  ! Interploates vctra values at eleva heights linearly by height vctrb
  ! values at elevb heights. For any elevb heights outside the range of
  ! eleva heights, vctrb values are set equal to the first or last vctra value. 
  implicit none
  
  integer, intent(in)                 :: na,nb
  real,    intent(in) , dimension(na) :: vctra, eleva
  real,    intent(in) , dimension(nb) :: elevb
  real,    intent(out), dimension(nb) :: vctrb

  !Local variables
  integer :: ka,kb
  real :: grada
  
  ka=1
  do kb=1,nb
     do while (ka < na-1 .and. elevb(kb) > eleva(ka+1))
        ka=ka+1
     end do
     grada = (vctra(ka+1) - vctra(ka)) / (eleva(ka+1) - eleva(ka))

     if (elevb(kb) < eleva(1)) then
        vctrb(kb) = vctra(1)
     elseif (elevb(kb) > eleva(na)) then
        vctrb(kb) = vctra(na)
     else
        vctrb(kb) = vctra(ka) + grada * (elevb(kb) - eleva(ka))
     end if
  end do
  return
end subroutine hintrp_cc

