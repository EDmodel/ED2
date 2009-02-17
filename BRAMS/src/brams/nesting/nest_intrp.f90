!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

subroutine fmrefs1d(ngbegin,ngend)

use mem_grid
use ref_sounding
use mem_scratch
use rconstants

implicit none

integer :: ngbegin,ngend,ifm,icm,k
real :: c1,c2

!     Interpolate the fine mesh 1-d reference state variables.

c1 = rgas / (cp - rgas)
c2 = cp * (rgas / p00) ** c1
do ifm = ngbegin,ngend
   icm = nxtnest(ifm)
   if (icm .ge. 1) then
      do k = 1,nnzp(icm)
         vctr1(k) = th01dn(k,icm) * dn01dn(k,icm)
         vctr2(k) = u01dn(k,icm) * dn01dn(k,icm)
         vctr3(k) = v01dn(k,icm) * dn01dn(k,icm)
         vctr4(k) = rt01dn(k,icm) * dn01dn(k,icm)
      enddo

     call eintp(dn01dn(1,icm),dn01dn(1,ifm),maxnzp,1,1,nnzp(ifm)  &
        ,1,1,ifm,1,'t',0,0)
     call eintp(vctr1,vctr5,maxnzp,1,1,nnzp(ifm),1,1,ifm,1,'t',0,0)
     call eintp(vctr2,vctr6,maxnzp,1,1,nnzp(ifm),1,1,ifm,1,'t',0,0)
     call eintp(vctr3,vctr7,maxnzp,1,1,nnzp(ifm),1,1,ifm,1,'t',0,0)
     call eintp(vctr4,vctr8,maxnzp,1,1,nnzp(ifm),1,1,ifm,1,'t',0,0)

      do k = 1,nnzp(ifm)
         th01dn(k,ifm) = vctr5(k) / dn01dn(k,ifm)
         u01dn(k,ifm) = vctr6(k) / dn01dn(k,ifm)
         v01dn(k,ifm) = vctr7(k) / dn01dn(k,ifm)
         rt01dn(k,ifm) = vctr8(k) / dn01dn(k,ifm)
         pi01dn(k,ifm) = c2 * (dn01dn(k,ifm) * th01dn(k,ifm)) ** c1
      enddo
   endif
enddo
return
end

!     *****************************************************************

subroutine fmrefs3d(ifm,mynum)

use mem_basic
use mem_scratch
use mem_grid
use mem_nestb
use rconstants

implicit none

integer :: ifm,icm,l,i,j,k,mynum

real :: c1,c2

!     Interpolate the fine mesh 3-D reference state variables.

icm = nxtnest(ifm)
if (icm .eq. 0) return

!   Don't need B components here, so just pass in bux, etc.

call fmint3(nnzp(icm),nnxp(icm),nnyp(icm),nnzp(ifm),nnxp(ifm),nnyp(ifm)                    &
           ,maxnzp,maxnxp,maxnyp,ifm,icm,nnstbot(ifm),nnsttop(ifm),jdim,1,0,0,'t'          &
           ,basic_g(icm)%dn0,basic_g(ifm)%dn0,basic_g(icm)%dn0,basic_g(ifm)%dn0            &
           ,scratch%scr1,scratch%scr2,grid_g(ifm)%topt,scratch%vt2da                       &
           ,nbounds(ifm)%bux,nbounds(ifm)%buy,nbounds(ifm)%buz,mynum)

call fmint3(nnzp(icm),nnxp(icm),nnyp(icm),nnzp(ifm),nnxp(ifm),nnyp(ifm)                    &
           ,maxnzp,maxnxp,maxnyp,ifm,icm,nnstbot(ifm),nnsttop(ifm),jdim,1,1,0,'t'          &
           ,basic_g(icm)%th0,basic_g(ifm)%th0,basic_g(icm)%dn0,basic_g(ifm)%dn0            &
           ,scratch%scr1,scratch%scr2,grid_g(ifm)%topt,scratch%vt2da                       &
           ,nbounds(ifm)%bux,nbounds(ifm)%buy,nbounds(ifm)%buz,mynum)

c1 = rgas / (cp - rgas)
c2 = cp * (rgas / p00) ** c1
do j = 1,nnyp(ifm)
   do i = 1,nnxp(ifm)
      do k = 1,nnzp(ifm)
         basic_g(ifm)%pi0(k,i,j) = c2  &
            * (basic_g(ifm)%dn0(k,i,j) * basic_g(ifm)%th0(k,i,j)) ** c1
      enddo
   enddo
enddo

call fillscr(1,maxnxp,maxnyp,1,nnxp(icm),nnyp(icm),1,1,scratch%scr1,grid_g(icm)%topt)
call eintp(scratch%scr1,scratch%scr2,1,maxnxp,maxnyp,1,nnxp(ifm),nnyp(ifm),ifm,2,'t',0,0)
call fillvar(1,maxnxp,maxnyp,1,nnxp(ifm),nnyp(ifm),1,1,scratch%scr2,scratch%scr1)

call rtgintrp(nnzp(ifm),nnxp(ifm),nnyp(ifm),basic_g(ifm)%th0,scratch%scr1                  &
             ,grid_g(ifm)%topt,ifm,'t')
call rtgintrp(nnzp(ifm),nnxp(ifm),nnyp(ifm),basic_g(ifm)%pi0,scratch%scr1                  &
             ,grid_g(ifm)%topt,ifm,'t')


! Define dn0u and dn0v

call fill_dn0uv(nnzp(ifm),nnxp(ifm),nnyp(ifm),basic_g(ifm)%dn0,basic_g(ifm)%dn0u           &
               ,basic_g(ifm)%dn0v )

return
end

!*****************************************************************************

subroutine fmdn0(ifm)

use mem_scratch
use mem_basic
use mem_grid

implicit none

integer :: ifm,icm

!     Special vertical interpolation of DN0 must be done after all other
!     3-D reference state and prognostic variables are interpolated.

icm = nxtnest(ifm)
if (icm .eq. 0) return
call fillscr(1,maxnxp,maxnyp,1,nnxp(icm),nnyp(icm),1,1,scratch%scr1,grid_g(icm)%topt)
call eintp(scratch%scr1,scratch%scr2,1,maxnxp,maxnyp,1,nnxp(ifm),nnyp(ifm),ifm,2,'t',0,0)
call fillvar(1,maxnxp,maxnyp,1,nnxp(ifm),nnyp(ifm),1,1,scratch%scr2,scratch%scr1)

call rtgintrp(nnzp(ifm),nnxp(ifm),nnyp(ifm),basic_g(ifm)%dn0,scratch%scr1                  &
             ,grid_g(ifm)%topt,ifm,'t')

! Define dn0u and dn0v

call fill_dn0uv(nnzp(ifm),nnxp(ifm),nnyp(ifm),basic_g(ifm)%dn0,basic_g(ifm)%dn0u           &
               ,basic_g(ifm)%dn0v)

return
end


!*****************************************************************************

subroutine fmint3(m1c,m2c,m3c,n1f,n2f,n3f  &
     ,maxnzp,maxnxp,maxnyp,ifm,icm  &
     ,nstbot,nsttop,jdim,initflg,idwt,irtgflg,vnam  &
     ,varc,varf,dn0xc,dn0xf,scr1,scr2,toptf,vt2da,bx,by,bz,mynum)

  ! fmint3: 
  !  one way nesting from coarser grid to finner grid:
  !  (1) during initialization, interpolate coarser field varc
  !      into finer field varc;
  !  (2) ow, store boundary conditions for finner field on arrays
  !      bx, by and bz, by interpolating coarser field varc

  use grid_dims, only: &
       maxgrds,        &
       nxpmax,         &
       nypmax,         &
       nzpmax

  use mem_grid, only: &
       nxtnest,       &
       ek1, ek2, ek3, &
       ek4, ek5, ek6, &
       ek7,           &
       ei1, ei2, ei3, &
       ei4, ei5, ei6, &
       ei7,           &
       ej1, ej2, ej3, &
       ej4, ej5, ej6, &
       ej7,           &
       ipm,           &
       jpm,           &
       kpm,           &
       nnzp

  implicit none
  integer, intent(in) :: m1c
  integer, intent(in) :: m2c
  integer, intent(in) :: m3c
  integer, intent(in) :: n1f
  integer, intent(in) :: n2f
  integer, intent(in) :: n3f
  integer, intent(in) :: maxnzp
  integer, intent(in) :: maxnxp
  integer, intent(in) :: maxnyp
  integer, intent(in) :: ifm
  integer, intent(in) :: icm
  integer, intent(in) :: nstbot
  integer, intent(in) :: nsttop
  integer, intent(in) :: jdim
  integer, intent(in) :: initflg
  integer, intent(in) :: idwt
  integer, intent(in) :: irtgflg
  integer, intent(in) :: mynum

  real, dimension(m1c,m2c,m3c) :: varc
  real, dimension(n1f,n2f,n3f) :: varf
  real, dimension(*) :: dn0xc,dn0xf,scr1,scr2,toptf,vt2da,bx,by,bz
  character(len=*) :: vnam

  integer :: ivartyp2
  character(len=10) :: c0, c1, c2, c3, c4, c5, c6, c7, c8
  character(len=*), parameter :: h="**(fmint3)**"

  if (icm .eq. 0) return

  ! copy coarser field varc to scr1
  ! scr1, supposed to be dimensioned (maxnzp, maxnxp, maxnyp) is a
  ! scratch area large enough to contain any of finner grid, coarser
  ! grid.

  call fillscr(maxnzp,maxnxp,maxnyp,m1c,m2c,m3c,1,m1c,scr1,varc)

  ! if density weigted, multiply scr1 by coarser density

  if (idwt .eq. 1) then
     call dnswt2(maxnzp,maxnxp,maxnyp,m1c,m2c,m3c  &
          ,scr1,dn0xc,vnam,1)
  end if

  !bob 16 may 2000:  Call parallel routine PAR_BINTP in place of sequential
  !bob    routine EINTP for nesting interpolation if simulation is 3D and
  !bob    for timestepping, not initialization.  This is to help eliminate
  !bob    machine epsilon differences between parallel and sequential runs.

  ! if 3D and not initialization then
  !   fill finner field boundary arrays (bx, by, bz) with interpolated scr1
  ! else
  !   build finner full field scr2 by interpolating scr1
  !   if density weigted and not initialization, then 
  !       divide scr2 by finner density
  !   end if
  !   if initialization then
  !       copy scr2 to varf
  !       if necessary, account for terrain differences on varf
  !   else
  !       fill finner field boundary arrays (bx, by, bz) by copying scr2
  !   end if
  ! end if

  if (jdim .eq. 1 .and. initflg .eq. 0) then
     if (vnam .eq. 'u') ivartyp2 = 1
     if (vnam .eq. 'v') ivartyp2 = 2
     if (vnam .eq. 'w') ivartyp2 = 3
     if (vnam .eq. 't' .and. idwt .eq. 0) ivartyp2 = 4
     if (vnam .eq. 't' .and. idwt .eq. 1) ivartyp2 = 5
     call par_bintp(scr1,scr2,dn0xf,maxnzp,maxnxp,maxnyp,n1f,n1f,n2f,n3f, &
          ifm,ivartyp2,0,0,15,bx,by,bz,mynum)

  else

     call eintp(scr1,scr2,maxnzp,maxnxp,maxnyp,n1f,n2f,n3f,ifm,3,vnam,0,0)

     if (idwt .eq. 1 .and. initflg .ge. 0) then
        call dnswt2(maxnzp,maxnxp,maxnyp,n1f,n2f,n3f,scr2,dn0xf,vnam,2)
     end if

     if (initflg .eq. 1) then
        call fillvar(maxnzp,maxnxp,maxnyp,n1f,n2f,n3f,1,n1f,scr2,varf)
        if(irtgflg.eq.1) then
           call rtgintrp(n1f,n2f,n3f,varf,vt2da,toptf,ifm,vnam)
        end if
     else
        call nstb(maxnzp,maxnxp,maxnyp,n1f,n2f,n3f,scr2,bx,by,bz  &
             ,vnam)
     endif

  endif
end subroutine fmint3


!******************************************************************************

subroutine rtgintrp(n1,n2,n3,fld,vt2da,topt,ngr,vpnt)

  use mem_grid
  use mem_scratch

  implicit none

  integer :: n1,n2,n3,ngr,i,j,k
  real, dimension(n1,n2,n3) :: fld
  real, dimension(n2,n3) :: vt2da,topt
  character(len=*) :: vpnt

  !     Do special vertical interpolation in case terrain on this grid
  !     (topt) is different from what would be interpolated from the
  !     coarser grid (vt2da).

  do j = 1,n3
     do i = 1,n2
        if (vpnt .eq. 'W') then
           do k = 1,n1
              vctr1(k) = zmn(k,ngr) * (1. - vt2da(i,j) / ztop) + vt2da(i,j)
              vctr2(k) = zmn(k,ngr) * (1. - topt(i,j) / ztop) + topt(i,j)
              vctr3(k) = fld(k,i,j)
           enddo
        else
           do k = 1,n1
              vctr1(k) = ztn(k,ngr) * (1. - vt2da(i,j) / ztop) + vt2da(i,j)
              vctr2(k) = ztn(k,ngr) * (1. - topt(i,j) / ztop) + topt(i,j)
              vctr3(k) = fld(k,i,j)
           enddo
        endif
        call htint(n1,vctr3,vctr1,n1,fld(1,i,j),vctr2)
     enddo
  enddo
  return
end subroutine rtgintrp

!*****************************************************************************

subroutine nstb(mx1,mx2,mx3,n1,n2,n3,vt3da,bx,by,bz,vnam)

  use mem_grid

  implicit none
  integer :: mx1,mx2,mx3,n1,n2,n3,nxfm,nyfm,nzfm,i,j,k
  real, dimension(mx1,mx2,mx3) :: vt3da
  real, dimension(n1,n3,2) :: bx
  real, dimension(n1,n2,2) :: by
  real, dimension(n2,n3,2) :: bz
  character(len=*) ::  vnam

  nxfm = n2
  nyfm = n3
  nzfm = n1
  if (vnam .eq. 'u') nxfm = n2 - 1
  if (vnam .eq. 'v') nyfm = (n3 - 2) * jdim + 1
  if (vnam .eq. 'w') nzfm = n1 - 1

  do j = 1,n3
     do k = 1,n1
        bx(k,j,1) = vt3da(k,1,j)
        bx(k,j,2) = vt3da(k,nxfm,j)
     enddo
  enddo

  if (jdim .eq. 1) then
     do i = 1,n2
        do k = 1,n1
           by(k,i,1) = vt3da(k,i,1)
           by(k,i,2) = vt3da(k,i,nyfm)
        enddo
     enddo
  endif

  if (nstbot .eq. 0) then
     do j = 1,n3
        do i = 1,n2
           bz(i,j,1) = vt3da(1,i,j)
        enddo
     enddo
  endif

  if (nsttop .eq. 0) then
     do j = 1,n3
        do i = 1,n2
           bz(i,j,2) = vt3da(nzfm,i,j)
        enddo
     enddo
  endif
  return
end subroutine nstb

!*****************************************************************************


subroutine nstbtnd(m1,m2,m3,ia,iz,ja,jz,ibcon  &
     ,scp,sct,bx,by,bz,vnam,tymeinv,nstbot,nsttop,jdim)

  implicit none
  integer :: m1,m2,m3,ia,iz,ja,jz,ibcon,nstbot,nsttop,jdim,i,j,k  &
       ,nzfm,nxfm,nyfm,incia,inciz,incja,incjz

  real :: tymeinv
  real, dimension(m1,m2,m3) :: scp,sct
  real, dimension(m1,m3,2) :: bx
  real, dimension(m1,m2,2) :: by
  real, dimension(m2,m3,2) :: bz
  character(len=*) :: vnam

  nzfm = m1
  nxfm = iz + 1
  nyfm = jz + 1
  if (vnam .eq. 'w') nzfm =nzfm - 1
  if (vnam .eq. 'u') nxfm =nxfm - 1
  if (vnam .eq. 'v') nyfm =nyfm - 1
  if (jdim .eq. 0)  nyfm = 1

  incia = 0
  inciz = 0
  if (iand(ibcon,1) .ne. 0) incia = 1
  if (iand(ibcon,2) .ne. 0 .and. vnam .ne. 'u') inciz = 1
  incja = 0
  incjz = 0
  if (iand(ibcon,4) .ne. 0) incja = jdim
  if (iand(ibcon,8) .ne. 0 .and. vnam .ne. 'v') incjz = jdim

  if (nstbot .eq. 0) then
     do j = ja,jz
        do i = ia,iz
           sct(1,i,j) = (bz(i,j,1) - scp(1,i,j)) * tymeinv
        enddo
     enddo
  endif
  if (nsttop .eq. 0) then
     do j = ja,jz
        do i = ia,iz
           sct(nzfm,i,j) = (bz(i,j,2) - scp(nzfm,i,j)) * tymeinv
        enddo
     enddo
  endif

  if (iand(ibcon,1) .ne. 0) then
     do j = ja-incja,jz+incjz
        do k = 1,m1
           sct(k,1,j) = (bx(k,j,1) - scp(k,1,j)) * tymeinv
        enddo
     enddo
  endif

  if (iand(ibcon,2) .ne. 0) then
     do j = ja-incja,jz+incjz
        do k = 1,m1
           sct(k,nxfm,j) = (bx(k,j,2) - scp(k,nxfm,j)) * tymeinv
        enddo
     enddo
  endif

  if (jdim .eq. 1) then
     if (iand(ibcon,4) .ne. 0) then
        do i = ia-incia,iz+inciz
           do k = 1,m1
              sct(k,i,1) = (by(k,i,1) - scp(k,i,1)) * tymeinv
           enddo
        enddo
     endif
     if (iand(ibcon,8) .ne. 0) then
        do i = ia-incia,iz+inciz
           do k = 1,m1
              sct(k,i,nyfm) = (by(k,i,2) - scp(k,i,nyfm)) * tymeinv
           enddo
        enddo
     endif
  endif
  return
end subroutine nstbtnd

!*****************************************************************************

subroutine eintp(ac,as,n1m,n2m,n3m,n1f,n2f,n3f,ifm,ndim,ispnt,i0,j0)

  ! eintp:
  !  interpolate one field from coarser into finner grid.
  !  finner grid is identified by ifm; coarser grid by nxtnest(ifm).
  !  ac is input coarser grid field, destroyed during computation.
  !  as is output interpolated finner grid field.
  !
  !  finner grid dimensions are n1f, n2f, n3f.
  !  both ac and as are dimensioned (n1m, n2m, n3m); these dimensions should 
  !  accomodate a coarser grid field as well as a finner grid field.
  !
  !  ispnt indicates the kind of field to interpolate;
  !  if ispnt is one of "u", "v", or "w", special four levels interpolation 
  !  is performed, if required by the interpolation direction.
  !
  !  ndim defines interpolation direction: (ic, jc, kc indicate coarser
  !    grid limits that correspond to finner grid limits)
  !  ndim==1:  interpolate on z:
  !            as(1:n1f, 1:n2f, 1:n3f) from ac(kc1:kc2, 1:n2f, 1:n3f)
  !  ndim==2:  interpolate on x and y:
  !            as(kc1:kc2, 1:n2f, 1:n3f) from ac(kc1:kc2, ic1:ic2, jc1:jc2)
  !  ndim==3:  interpolate on all directions:
  !            as(1:n1f, 1:n2f, 1:n3f) from ac(kc1:kc2, ic1:ic2, jc1:jc2)
  !
  ! CHECK THIS INFO!!! Code very hard to understand !!!

  use mem_grid, only: &
       jdim,          &
       nxtnest,       &
       nnzp,          &
       ipm,           &
       ei1, ei2, ei3, &
       ei4, ei5, ei6, &
       ei7,           &
       jpm,           &
       ej1, ej2, ej3, &
       ej4, ej5, ej6, &
       ej7,           &
       kpm,           &
       ek1, ek2, ek3, &
       ek4, ek5, ek6, &
       ek7

  implicit none

  integer,          intent(in   ) :: n1m
  integer,          intent(in   ) :: n2m
  integer,          intent(in   ) :: n3m
  integer,          intent(in   ) :: n1f
  integer,          intent(in   ) :: n2f
  integer,          intent(in   ) :: n3f
  integer,          intent(in   ) :: ifm
  integer,          intent(in   ) :: ndim
  character(len=*), intent(in   ) :: ispnt
  integer,          intent(in   ) :: i0
  integer,          intent(in   ) :: j0
  real,             intent(inout) :: ac(n1m,n2m,n3m)
  real,             intent(inout) :: as(n1m,n2m,n3m)

  integer :: icm
  integer :: k1
  integer :: k2
  integer :: jc1
  integer :: jc
  integer :: kc
  integer :: if
  integer :: jf
  integer :: kf

  ! return if ifm is the coarsest grid of all

  icm = nxtnest(ifm)
  if (icm == 0) return

  ! interpolation vertical limits 

  if(ndim /= 2) then
     k1 = max(1,kpm(2,ifm)-2)
     k2 = min(nnzp(icm),kpm(n1f-1,ifm)+2)
  else
     k1 = 1
     k2 = 1
  end if

  if (ndim /= 1) then

     ! interpolate coarser points of array ac at x direction,
     ! storing results at array as

     if (ispnt .eq. 'u') then
        do jc1 = jpm(1,ifm)-jdim,jpm(n3f,ifm)+jdim
           jc = jc1 - j0
           do if = 1,n2f-1
              do kc = k1,k2
                 as(kc,if,jc) = &
                      ei4(if,ifm) * ac(kc,ipm(if,ifm)-i0-2,jc)  + &
                      ei5(if,ifm) * ac(kc,ipm(if,ifm)-i0-1,jc)  + &
                      ei6(if,ifm) * ac(kc,ipm(if,ifm)-i0,  jc)  + &
                      ei7(if,ifm) * ac(kc,ipm(if,ifm)-i0+1,jc)
              end do
           end do
        end do
     else
        do jc1 = jpm(1,ifm)-jdim,jpm(n3f,ifm)+jdim
           jc = jc1 - j0
           do if = 1,n2f
              do kc = k1,k2
                 as(kc,if,jc) = &
                      ei1(if,ifm) * ac(kc,ipm(if,ifm)-i0-1,jc)  + &
                      ei2(if,ifm) * ac(kc,ipm(if,ifm)-i0  ,jc)  + &
                      ei3(if,ifm) * ac(kc,ipm(if,ifm)-i0+1,jc)
              end do
           end do
        end do
     end if

     if (ispnt .eq. 'v' .and. jdim .eq. 1) then
        do jf = 1,n3f-1
           do if = 1,n2f
              do kc = k1,k2
                 ac(kc,if,jf) = &
                      ej4(jf,ifm) * as(kc,if,jpm(jf,ifm)-j0-2)  + &
                      ej5(jf,ifm) * as(kc,if,jpm(jf,ifm)-j0-1)  + &
                      ej6(jf,ifm) * as(kc,if,jpm(jf,ifm)-j0  )  + &
                      ej7(jf,ifm) * as(kc,if,jpm(jf,ifm)-j0+1)
              end do
           end do
        end do
     else if (jdim .eq. 1) then
        if (ndim .eq. 2) then
           do jf = 1,n3f
              do if = 1,n2f
                 do kc = k1,k2
                    ac(kc,if,jf) = &
                         ej1(jf,ifm) * as(kc,if,jpm(jf,ifm)-j0-1)  + &
                         ej2(jf,ifm) * as(kc,if,jpm(jf,ifm)-j0  )  + &
                         ej3(jf,ifm) * as(kc,if,jpm(jf,ifm)-j0+1)
                 end do
              end do
           end do

           do jf = 1,n3f
              do if = 1,n2f
                 do kf = k1,k2
                    as(kf,if,jf) = ac(kf,if,jf)
                 end do
              end do
           end do
           return
        else
           do jf = 1,n3f
              do if = 1,n2f
                 do kc = k1,k2
                    ac(kc,if,jf) = &
                         ej1(jf,ifm) * as(kc,if,jpm(jf,ifm)-j0-1)  + &
                         ej2(jf,ifm) * as(kc,if,jpm(jf,ifm)-j0  )  + &
                         ej3(jf,ifm) * as(kc,if,jpm(jf,ifm)-j0+1)
                 end do
              end do
           end do
        end if
     else
        if (ndim .eq. 2) then
           !            do if=1,n2f
           !               do kc=k1,k2
           !                  ac(kc,if,1)=as(kc,if,1)
           !               end do
           !            end do
           return
        else
           do if = 1,n2f
              do kc = k1,k2
                 ac(kc,if,1) = as(kc,if,1)
              end do
           end do
        end if
     end if
  end if


  ! at this point, array ac is supposed to be interpolated
  ! at x and y;
  ! store, at array as, the vertical interpolation of ac.


  if (ispnt .eq. 'w') then
     do jf = 1,n3f
        do if = 1,n2f
           do kf = 1,n1f-1
              kc = kpm(kf+1,ifm)
              as(kf,if,jf) = &
                   ek4(kf,ifm) * ac(max(1          ,kc-2),if,jf)  + &
                   ek5(kf,ifm) * ac(                kc-1 ,if,jf)  + &
                   ek6(kf,ifm) * ac(                kc   ,if,jf)  + &
                   ek7(kf,ifm) * ac(min(nnzp(icm)-1,kc+1),if,jf)
           end do
        end do
     end do
  else
     do jf = 1,n3f
        do if = 1,n2f
           do kf = 1,n1f
              as(kf,if,jf) = &
                   ek1(kf,ifm) * ac(kpm(kf,ifm)-1,if,jf)  + &
                   ek2(kf,ifm) * ac(kpm(kf,ifm)  ,if,jf)  + &
                   ek3(kf,ifm) * ac(kpm(kf,ifm)+1,if,jf)
           end do
        end do
     end do
  end if
end subroutine eintp

!     *****************************************************************

subroutine cofnest

  use mem_grid
  use mem_scratch

  implicit none

  integer :: nf,nc,nrat,if,jf,kf,kc
  real :: alpha,et,ev

  do nf = 2,ngrids
     nc = nxtnest(nf)
     if (nc .eq. 0) go to 50
     nrat = nstratx(nf)
     alpha = ((1. / float(nrat)) ** 2 - 1.) / 24.
     do if = 1,nnxp(nf)
        et = -.5 + float(2 * mod(if+nrat-2,nrat) + 1) / (2.0 * float(nrat))
        ev = -.5 + float(mod(if+nrat-2,nrat) + 1) / float(nrat)

        ei1(if,nf) = et * (et - 1.) / 2. + alpha
        ei2(if,nf) = (1. - et * et) - 2. * alpha
        ei3(if,nf) = et * (et + 1.) / 2. + alpha
        ei4(if,nf) = (ev * ev - 0.25) * (1.5 - ev) / 6.0
        ei5(if,nf) = (0.5 - ev) * (2.25 - ev * ev) * 0.5
        ei6(if,nf) = (0.5 + ev) * (2.25 - ev * ev) * 0.5
        ei7(if,nf) = (ev * ev - 0.25) * (1.5 + ev) / 6.0
     enddo

     if (jdim .eq. 1) then
        nrat = nstraty(nf)
        alpha = ((1. / float(nrat)) ** 2 - 1.) / 24.
        do jf = 1,nnyp(nf)
           et = -.5 + float(2 * mod(jf+nrat-2,nrat) + 1) / (2.0 * float(nrat))
           ev = -.5 + float(mod(jf+nrat-2,nrat) + 1) / float(nrat)

           ej1(jf,nf) = et * (et - 1.) / 2. + alpha
           ej2(jf,nf) = (1. - et * et) - 2. * alpha
           ej3(jf,nf) = et * (et + 1.) / 2. + alpha
           ej4(jf,nf) = (ev * ev - 0.25) * (1.5 - ev) / 6.0
           ej5(jf,nf) = (0.5 - ev) * (2.25 - ev * ev) * 0.5
           ej6(jf,nf) = (0.5 + ev) * (2.25 - ev * ev) * 0.5
           ej7(jf,nf) = (ev * ev - 0.25) * (1.5 + ev) / 6.0
        enddo
     endif

     do kc = 1,nnzp(nc)
        vctr1(kc) = 0.
        vctr2(kc) = 0.
        vctr3(kc) = 0.
     enddo

     do kf = 1,nnzp(nf)
        kc = kpm(kf,nf)
        ek1(kf,nf) = (ztn(kf,nf) - ztn(kc,nc)) * (ztn(kf,nf) - ztn(kc+1,nc))  &
             / ((ztn(kc-1,nc) - ztn(kc,nc)) * (ztn(kc-1,nc) - ztn(kc+1,nc)))
        ek2(kf,nf) = (ztn(kf,nf) - ztn(kc-1,nc)) * (ztn(kf,nf) - ztn(kc+1,nc  &
             )) / ((ztn(kc,nc) - ztn(kc-1,nc)) * (ztn(kc,nc) - ztn(kc+1,nc)))
        ek3(kf,nf) = (ztn(kf,nf) - ztn(kc-1,nc)) * (ztn(kf,nf) - ztn(kc,nc))  &
             / ((ztn(kc+1,nc) - ztn(kc-1,nc)) * (ztn(kc+1,nc) - ztn(kc,nc)))
        if (kf .ge. 2 .and. kf .le. nnzp(nf) - 1) then
           vctr1(kc) = vctr1(kc) + (zmn(kf,nf) - zmn(kf-1,nf)) * ek1(kf,nf)
           vctr2(kc) = vctr2(kc) + (zmn(kf,nf) - zmn(kf-1,nf)) * ek2(kf,nf)
           vctr3(kc) = vctr3(kc) + (zmn(kf,nf) - zmn(kf-1,nf)) * ek3(kf,nf)
        endif
     enddo
     !c         if(kpm(2,nf).gt.2)then
     !c            vctr1(1)=vctr1(1+nrz(kpm(2,nf),nc)
     !c            vctr2(1)=vctr2(1)
     !c            vctr3(1)=vctr3(1)
     !c         endif
     !c         if(kpm(nnzp(nf)-1,nf).lt.nnzp(nc)-1)then
     !c            vctr1(nnzp(nf))=vctr1(nnzp(nf)-1)
     !c            vctr2(nnzp(nf))=vctr2(nnzp(nf)-1)
     !c            vctr3(nnzp(nf))=vctr3(nnzp(nf)-1)
     !c         endif
     do kf = 2,nnzp(nf) - 1
        kc = kpm(kf,nf)
        ek1(kf,nf) = ek1(kf,nf) - dztn(kc,nc) * vctr1(kc)
        ek2(kf,nf) = ek2(kf,nf) - dztn(kc,nc) * vctr2(kc) + 1.
        ek3(kf,nf) = ek3(kf,nf) - dztn(kc,nc) * vctr3(kc)
     enddo
     !cccccccccccccccc
     !cc         ek1(1,nf)=ek1(2,nf)
     !cc         ek2(1,nf)=ek2(2,nf)
     !cc         ek3(1,nf)=ek3(2,nf)
     !cccccccccccccccccc

     do kf = 1,nnzp(nf)-1
        kc = kpm(kf+1,nf)
        if (kf .eq. 1 .or. kf .eq. nnzp(nf)-1 .or. kc .ne. kpm(kf,nf)) then
           ek4(kf,nf) = 0.
           ek5(kf,nf) = 1.
           ek6(kf,nf) = 0.
           ek7(kf,nf) = 0.
        else
           ek4(kf,nf) = ek4(kf-1,nf) + (zmn(kf,nf) - zmn(kf-1,nf))  &
                * (-ek1(kf,nf)) * dztn(kc-1,nc)
           ek5(kf,nf) = ek5(kf-1,nf) + (zmn(kf,nf) - zmn(kf-1,nf))  &
                * (ek1(kf,nf) * dztn(kc-1,nc) - ek2(kf,nf) * dztn(kc,nc))
           ek6(kf,nf) = ek6(kf-1,nf) + (zmn(kf,nf) - zmn(kf-1,nf))  &
                * (ek2(kf,nf) * dztn(kc,nc) - ek3(kf,nf) * dztn(kc+1,nc))
           ek7(kf,nf) = ek7(kf-1,nf) + (zmn(kf,nf) - zmn(kf-1,nf))  &
                * ek3(kf,nf) * dztn(kc+1,nc)
        endif
     enddo
     do kf = 1,nnzp(nf)
        fbcf(kf,nf,1) = dztn(kpm(kf,nf),nc) / (float(nstraty(nf)) * dztn(kf,nf))
        fbcf(kf,nf,2) = dztn(kpm(kf,nf),nc) / (float(nstratx(nf)) * dztn(kf,nf))
        fbcf(kf,nf,3) = 1. / float(nstratx(nf) * nstraty(nf))
        fbcf(kf,nf,4) = dztn(kpm(kf,nf),nc)  &
             / (float(nstratx(nf) * nstraty(nf)) * dztn(kf,nf))
     enddo
50   continue
  enddo
  return
end subroutine cofnest


!     *****************************************************************

subroutine fmint4(var1,var2,dn0xc,dn0xf,vt2da,ifm,icm,vpnt,idwt)

  use mem_scratch
  use mem_grid

  implicit none
  integer :: ifm,icm,idwt

  real, dimension(*) :: var1,var2,vt2da,dn0xc,dn0xf
  character(len=*) :: vpnt

  if (icm .eq. 0) return

  call fillscr(maxnzp,maxnxp,maxnyp,nnzp(icm),nnxp(icm),nnyp(icm)  &
              ,1,nnzp(icm),scratch%scr1,var1)


  if (idwt .eq. 1) then
     call dnswt2(maxnzp,maxnxp,maxnyp,nnzp(icm),nnxp(icm),nnyp(icm)  &
                ,scratch%scr1,dn0xc,vpnt,1)
  endif

  call eintp(scratch%scr1,scratch%scr2,maxnzp,maxnxp,maxnyp  &
       ,nnzp(ifm),nnxp(ifm),nnyp(ifm),ifm,3,vpnt,0,0)

  call fillvar(maxnzp,maxnxp,maxnyp,nnzp(ifm),nnxp(ifm),nnyp(ifm)  &
       ,1,nnzp(ifm),scratch%scr2,var2)

  if (idwt .eq. 1) call dnswt2(nnzp(ifm),nnxp(ifm),nnyp(ifm),nnzp(ifm)  &
       ,nnxp(ifm),nnyp(ifm),var2,dn0xf,vpnt,2)

  call rtgintrp(nnzp(ifm),nnxp(ifm),nnyp(ifm),var2,vt2da  &
       ,grid_g(ifm)%topt,ifm,vpnt)

  return
end subroutine fmint4


!     *****************************************************************

subroutine fmint2d(icm,ifm,vpnt,var1,var2)

  use mem_scratch
  use mem_grid

  implicit none

  integer :: ifm,icm
  real, dimension(*) :: var1,var2
  character(len=*) :: vpnt

  if (icm == 0) return

  call fillscr(1,maxnxp,maxnyp,1,nnxp(icm),nnyp(icm),1,1,scratch%scr1,var1)

  call eintp(scratch%scr1,scratch%scr2,1,maxnxp,maxnyp,1,nnxp(ifm),nnyp(ifm)               &
            ,ifm,2,vpnt,0,0)

  call fillvar(1,maxnxp,maxnyp,1,nnxp(ifm),nnyp(ifm),1,1,scratch%scr2,var2)

  return
end subroutine fmint2d





