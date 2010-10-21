!
! Copyright (C) 1991-2004  ; All Rights Reserved ; Colorado State University
! Colorado State University Research Foundation ; ATMET, LLC
! 
! This file is free software; you can redistribute it and/or modify it under the
! terms of the GNU General Public License as published by the Free Software 
! Foundation; either version 2 of the License, or (at your option) any later version.
! 
! This software is distributed in the hope that it will be useful, but WITHOUT ANY 
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
! PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along with this 
! program; if not, write to the Free Software Foundation, Inc., 
! 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
!======================================================================================
!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################

subroutine par_bintp(ac,as,dn0f,n1m,n2m,n3m,m1f  &
     ,m1,m2,m3,ifm,ivnam,i0,j0,ibcon,bx,by,bz,mynum)
  ! par_bintp: 
  !  interpolates coarser field points into finner field boundaries.
  !
  !  That is, interpolate array ac from the coarser grid into the
  !  finner grid and extract boundaries, storing at
  !  finner grid boundary arrays bx, by, bz. 
  !
  !  In fact, the full finner grid field is not built. The algorithm
  !  first interpolate every horizontal plane from the coarser to the
  !  finner grid, keeping the vertical coarser coordinate. Finally, it
  !  interpolates vertically from coarser to finner only at boundaries.
  !
  !  Array ac is destroyed during computation. Array as is a scratch area.
  !
  !  Input argument ibcon is a binary code for which boundary to build.
  !  If bit *i* is set, boundary is built, otherwise it is not build.
  !  Bit 0 is low x boundary:
  !  Bit 1 is high x boundary:
  !  Bit 2 is low y boundary:
  !  Bit 3 is high y boundary:
  !  nstbot controls if low z boundary is to be build (0 means yes)
  !  nsttop controls if high z boundary is to be build (0 means yes)
  !**Old comments (maybe garbage by now):
  ! When calling par_bintp, pass in the following values for its local
  !   index limits ia,iz,ja,jz, which we are calling ia0,iz0,ja0,jz0
  !   in the call:
  !      ia0=iz
  !      if(ivnam.eq.1)ia0=iz-1
  !      iz0=1
  !      ja0=jz
  !      if(ivnam.eq.2)ja0=jz-1
  !      jz0=1
  !      if(iand(ibcon,2).ne.0)iz0=ia0
  !      if(iand(ibcon,1).ne.0)ia0=1
  !      if(iand(ibcon,8).ne.0)jz0=ja0
  !      if(iand(ibcon,4).ne.0)ja0=1
  !**End of old comments

  use mem_grid, only : &
         maxgrds &                   ! maximum number of grids
        ,nxpmax  &                   ! maximum x dimension
        ,nypmax  &                   ! maximum y dimension
        ,nzpmax  &                   ! maximum z dimension
        ,nxtnest & !(maxgrds)        ! next coarser grid number (0 if grid is not nested)
        ,nsttop  &                   ! top z boundary control (0 means build)
        ,nstbot  &                   ! high z boundary control (0 means build)
        ,nnzp    & !(maxgrds)        ! grid points z direction
        ,ipm     & !(nxpmax,maxgrds) ! next coarser grid cell index (icoarser) that contains this finer grid cell
        ,ei1     & !(nxpmax,maxgrds) ! interp. weight for icoarser-1 on 3 points interp
        ,ei2     & !(nxpmax,maxgrds) ! interp. weight for icoarser   on 3 points interp
        ,ei3     & !(nxpmax,maxgrds) ! interp. weight for icoarser+1 on 3 points interp
        ,ei4     & !(nxpmax,maxgrds) ! interp. weight for icoarser-2 on 4 points interp
        ,ei5     & !(nxpmax,maxgrds) ! interp. weight for icoarser-1 on 4 points interp
        ,ei6     & !(nxpmax,maxgrds) ! interp. weight for icoarser   on 4 points interp
        ,ei7     & !(nxpmax,maxgrds) ! interp. weight for icoarser+1 on 4 points interp
        ,jpm     & !(nypmax,maxgrds) ! next coarser grid cell index (jcoarser) that contains this finer grid cell
        ,ej1     & !(nypmax,maxgrds) ! interp. weight for jcoarser-1 on 3 points interp
        ,ej2     & !(nypmax,maxgrds) ! interp. weight for jcoarser   on 3 points interp
        ,ej3     & !(nypmax,maxgrds) ! interp. weight for jcoarser+1 on 3 points interp
        ,ej4     & !(nypmax,maxgrds) ! interp. weight for jcoarser-2 on 4 points interp
        ,ej5     & !(nypmax,maxgrds) ! interp. weight for jcoarser-1 on 4 points interp
        ,ej6     & !(nypmax,maxgrds) ! interp. weight for jcoarser   on 4 points interp
        ,ej7     & !(nypmax,maxgrds) ! interp. weight for jcoarser+1 on 4 points interp
        ,kpm     & !(nzpmax,maxgrds) ! next coarser grid cell index that contains this finer grid cell
        ,ek1     & !(nzpmax,maxgrds) ! interp. weight for kcoarser-1 on 3 points interp
        ,ek2     & !(nzpmax,maxgrds) ! interp. weight for kcoarser   on 3 points interp
        ,ek3     & !(nzpmax,maxgrds) ! interp. weight for kcoarser+1 on 3 points interp
        ,ek4     & !(nzpmax,maxgrds) ! interp. weight for kcoarser-2 on 4 points interp
        ,ek5     & !(nzpmax,maxgrds) ! interp. weight for kcoarser-1 on 4 points interp
        ,ek6     & !(nzpmax,maxgrds) ! interp. weight for kcoarser   on 4 points interp
        ,ek7       !(nzpmax,maxgrds) ! interp. weight for kcoarser+1 on 4 points interp
     
  implicit none
     
  integer, intent(in   ) :: n1m                 ! first dimension of ac and as; should accomodate both grids
  integer, intent(in   ) :: n2m                 ! second dimension of ac and as; should accomodate both grids
  integer, intent(in   ) :: n3m                 ! third dimension of ac and as; should accomodate both grids
  integer, intent(in   ) :: m1f                 ! last vertical level to receive interpolation on finner grid (?)
  integer, intent(in   ) :: m1                  ! first dimension of finner grid
  integer, intent(in   ) :: m2                  ! second dimension of finner grid
  integer, intent(in   ) :: m3                  ! third dimension of finner grid
  integer, intent(in   ) :: ifm                 ! finner grid id
  integer, intent(in   ) :: ivnam               ! 1="u", 2="v", 3="w", 4="t" and idt=0, 5="t" and idt=1
  integer, intent(in   ) :: i0
  integer, intent(in   ) :: j0
  integer, intent(in   ) :: ibcon               ! code for boundary building selection (see above)
  real,    intent(inout) :: ac(n1m,n2m,n3m)     ! input array to be interpolated, destroyed during computation
  real,    intent(inout) :: as(n1m,n2m,n3m)     ! scratch array
  real,    intent(in   ) :: dn0f(m1,m2,m3)      ! finner grid density
  real,    intent(out  ) :: bx(m1,m3,2)         ! finner grid x boundaries
  real,    intent(out  ) :: by(m1,m2,2)         ! finner grid y boundaries
  real,    intent(out  ) :: bz(m2,m3,2)         ! finner grid z boundaries
  integer, intent(in   ) :: mynum               ! ID number for this node 
  ! Internal
  integer :: ia,iz,ja,jz,nc,kc,ic,jc,kf,if,jf,k1,k2,im,jm

  !_____________________________________________________________________
  !
  !    Set various indices. DO NOT CHANGE THE ORDER OF THESE STATEMENTS!!

  ia=m2
  if(ivnam == 1)ia=m2-1
  iz=1

  ja=m3
  if(ivnam == 2)ja=m3-1
  jz=1

  if(ibcon /= 0) then
     if(ibcon /= 1)iz=ia
     if(ibcon /= 2)ia=1
     if(ibcon /= 4)jz=ja
     if(ibcon /= 8)ja=1
  endif
  !_____________________________________________________________________

  ! i0 and j0 now passed in as fm node offsets
  ! ia, iz, ja, and jz refer to locations where the b arrays need to
  !   be filled; they are not the usual limits of prognostic points.
  ! The nstbot=0 and nsttop=0 options do not work with the current ipaths
  !   array definitions.

  nc=nxtnest(ifm)
  k1=max(1,kpm(2,ifm)-2)
  k2=min(nnzp(nc),kpm(m1f-1,ifm)+2)

  ! expand the x direction of ac, storing 
  ! as (k1:k2, ia:iz, jac:jzc) <= ac (k1:k2, iac:izc, jac:jzc)
  ! where iac:izc, jac:jzc define the coarser grid region that
  ! impacts the finner grid 


  if(ivnam == 1)then
     ! for "u" fields
     do jc=jpm(ja+j0,ifm)-2,jpm(jz+j0,ifm)+1
        do if=ia,iz
           im=if+i0
           ic=ipm(im,ifm)
           do kc=k1,k2
              as(kc,if,jc)=ei4(im,ifm)*ac(kc,ic-2,jc)  &
                          +ei5(im,ifm)*ac(kc,ic-1,jc)  &
                          +ei6(im,ifm)*ac(kc,ic  ,jc)  &
                          +ei7(im,ifm)*ac(kc,ic+1,jc)
           enddo
        enddo
     enddo
  else

     ! remaining fields

     do jc=jpm(ja+j0,ifm)-2,jpm(jz+j0,ifm)+1
        do if=ia,iz
           im=if+i0
           ic=ipm(im,ifm)
           do kc=k1,k2
              as(kc,if,jc)=ei1(im,ifm)*ac(kc,ic-1,jc)  &
                          +ei2(im,ifm)*ac(kc,ic  ,jc)  &
                          +ei3(im,ifm)*ac(kc,ic+1,jc)
           enddo
        enddo
     enddo
  endif

  ! expand the y direction of as, storing 
  ! ac (k1:k2, ia:iz, ja:jz) <= as (k1:k2, ia:iz, jac:jzc)
  ! where jac:jzc defines the coarser grid region that
  ! impacts the finner grid 

  if(ivnam.eq.2)then

     ! for "v" fields

     do jf=ja,jz
        jm=jf+j0
        jc=jpm(jm,ifm)
        do if=ia,iz
           do kc=k1,k2
              ac(kc,if,jf)=ej4(jm,ifm)*as(kc,if,jc-2)  &
                          +ej5(jm,ifm)*as(kc,if,jc-1)  &
                          +ej6(jm,ifm)*as(kc,if,jc  )  &
                          +ej7(jm,ifm)*as(kc,if,jc+1)
           enddo
        enddo
     enddo
  else

     ! remaining fields

     do jf=ja,jz
        jm=jf+j0
        jc=jpm(jm,ifm)
        do if=ia,iz
           do kc=k1,k2
              ac(kc,if,jf)=ej1(jm,ifm)*as(kc,if,jc-1)  &
                          +ej2(jm,ifm)*as(kc,if,jc  )  &
                          +ej3(jm,ifm)*as(kc,if,jc+1)
           enddo
        enddo
     enddo
  endif

  ! at this point, ac is interpolated at x and y. 
  ! it remains to interpolate at z and store boundaries
  ! at bx, by, bz arrays.

  if(ivnam == 3)then

     ! for "w" fields

     if(iand(ibcon,1).ne.0) then
        if=ia
        do jf=ja,jz
           do kf=1,m1f-1
              kc=kpm(kf+1,ifm)
              bx(kf,jf,1)=  &
                  (ek4(kf,ifm)*ac(max(1         ,kc-2),if,jf)  &
                 + ek5(kf,ifm)*ac(               kc-1 ,if,jf)  &
                 + ek6(kf,ifm)*ac(               kc   ,if,jf)  &
                 + ek7(kf,ifm)*ac(min(nnzp(nc)-1,kc+1),if,jf))  &
                 / (.5 * (dn0f(kf,if,jf) + dn0f(kf+1,if,jf)))
           enddo
        enddo
     endif

     if(iand(ibcon,2).ne.0) then
        if=iz
        do jf=ja,jz
           do kf=1,m1f-1
              kc=kpm(kf+1,ifm)
              bx(kf,jf,2)=  &
                  (ek4(kf,ifm)*ac(max(1         ,kc-2),if,jf)  &
                 + ek5(kf,ifm)*ac(               kc-1 ,if,jf)  &
                 + ek6(kf,ifm)*ac(               kc   ,if,jf)  &
                 + ek7(kf,ifm)*ac(min(nnzp(nc)-1,kc+1),if,jf))  &
                 / (.5 * (dn0f(kf,if,jf) + dn0f(kf+1,if,jf)))
           enddo
        enddo
     endif

     if(iand(ibcon,4).ne.0) then
        jf=ja
        do if=ia,iz
           do kf=1,m1f-1
              kc=kpm(kf+1,ifm)
              by(kf,if,1)=  &
                  (ek4(kf,ifm)*ac(max(1         ,kc-2),if,jf)  &
                 + ek5(kf,ifm)*ac(               kc-1 ,if,jf)  &
                 + ek6(kf,ifm)*ac(               kc   ,if,jf)  &
                 + ek7(kf,ifm)*ac(min(nnzp(nc)-1,kc+1),if,jf))  &
                 / (.5 * (dn0f(kf,if,jf) + dn0f(kf+1,if,jf)))
           enddo
        enddo
     endif

     if(iand(ibcon,8).ne.0) then
        jf=jz
        do if=ia,iz
           do kf=1,m1f-1
              kc=kpm(kf+1,ifm)
              by(kf,if,2)=  &
                  (ek4(kf,ifm)*ac(max(1         ,kc-2),if,jf)  &
                 + ek5(kf,ifm)*ac(               kc-1 ,if,jf)  &
                 + ek6(kf,ifm)*ac(               kc   ,if,jf)  &
                 + ek7(kf,ifm)*ac(min(nnzp(nc)-1,kc+1),if,jf))  &
                 / (.5 * (dn0f(kf,if,jf) + dn0f(kf+1,if,jf)))
           enddo
        enddo
     endif

     if(nstbot.eq.0)then
        kf=1
        do jf=ja,jz
           do if=ia,iz
              kc=kpm(kf+1,ifm)
              bz(if,jf,1)=  &
                  (ek4(kf,ifm)*ac(max(1         ,kc-2),if,jf)  &
                 + ek5(kf,ifm)*ac(               kc-1 ,if,jf)  &
                 + ek6(kf,ifm)*ac(               kc   ,if,jf)  &
                 + ek7(kf,ifm)*ac(min(nnzp(nc)-1,kc+1),if,jf))  &
                 / (.5 * (dn0f(kf,if,jf) + dn0f(kf+1,if,jf)))
           enddo
        enddo
     endif

     if(nsttop.eq.0)then
        kf=m1f-1
        do jf=ja,jz
           do if=ia,iz
              kc=kpm(kf+1,ifm)
              bz(if,jf,2)=  &
                  (ek4(kf,ifm)*ac(max(1         ,kc-2),if,jf)  &
                 + ek5(kf,ifm)*ac(               kc-1 ,if,jf)  &
                 + ek6(kf,ifm)*ac(               kc   ,if,jf)  &
                 + ek7(kf,ifm)*ac(min(nnzp(nc)-1,kc+1),if,jf))  &
                 / (.5 * (dn0f(kf,if,jf) + dn0f(kf+1,if,jf)))
           enddo
        enddo
     endif

  elseif(ivnam.eq.4)then

     ! for "t" fields without density weigth

     if(iand(ibcon,1).ne.0) then
        if=ia
        do jf=ja,jz
           do kf=1,m1f
              bx(kf,jf,1)=  &
                  (ek1(kf,ifm)*ac(kpm(kf,ifm)-1,if,jf)  &
                 + ek2(kf,ifm)*ac(kpm(kf,ifm)  ,if,jf)  &
                 + ek3(kf,ifm)*ac(kpm(kf,ifm)+1,if,jf))
           enddo
        enddo
     endif

     if(iand(ibcon,2).ne.0) then
        if=iz
        do jf=ja,jz
           do kf=1,m1f
              bx(kf,jf,2)=  &
                  (ek1(kf,ifm)*ac(kpm(kf,ifm)-1,if,jf)  &
                 + ek2(kf,ifm)*ac(kpm(kf,ifm)  ,if,jf)  &
                 + ek3(kf,ifm)*ac(kpm(kf,ifm)+1,if,jf))
           enddo
        enddo
     endif

     if(iand(ibcon,4).ne.0) then
        jf=ja
        do if=ia,iz
           do kf=1,m1f
              by(kf,if,1)=  &
                  (ek1(kf,ifm)*ac(kpm(kf,ifm)-1,if,jf)  &
                 + ek2(kf,ifm)*ac(kpm(kf,ifm)  ,if,jf)  &
                 + ek3(kf,ifm)*ac(kpm(kf,ifm)+1,if,jf))
           enddo
        enddo
     endif

     if(iand(ibcon,8).ne.0) then
        jf=jz
        do if=ia,iz
           do kf=1,m1f
              by(kf,if,2)=  &
                  (ek1(kf,ifm)*ac(kpm(kf,ifm)-1,if,jf)  &
                 + ek2(kf,ifm)*ac(kpm(kf,ifm)  ,if,jf)  &
                 + ek3(kf,ifm)*ac(kpm(kf,ifm)+1,if,jf))
           enddo
        enddo
     endif

     if(nstbot.eq.0)then
        kf=1
        do jf=ja,jz
           do if=ia,iz
              bz(if,jf,1)=  &
                  (ek1(kf,ifm)*ac(kpm(kf,ifm)-1,if,jf)  &
                 + ek2(kf,ifm)*ac(kpm(kf,ifm)  ,if,jf)  &
                 + ek3(kf,ifm)*ac(kpm(kf,ifm)+1,if,jf))
           enddo
        enddo
     endif

     if(nsttop.eq.0)then
        kf=m1f
        do jf=ja,jz
           do if=ia,iz
              bz(if,jf,2)=  &
                  (ek1(kf,ifm)*ac(kpm(kf,ifm)-1,if,jf)  &
                 + ek2(kf,ifm)*ac(kpm(kf,ifm)  ,if,jf)  &
                 + ek3(kf,ifm)*ac(kpm(kf,ifm)+1,if,jf))
           enddo
        enddo
     endif

  else

     ! for "u" fields, "v" fields, or density weigthed "t" fields

     if(iand(ibcon,1).ne.0) then
        if=ia
        do jf=ja,jz
           do kf=1,m1f
              bx(kf,jf,1)=  &
                  (ek1(kf,ifm)*ac(kpm(kf,ifm)-1,if,jf)  &
                 + ek2(kf,ifm)*ac(kpm(kf,ifm)  ,if,jf)  &
                 + ek3(kf,ifm)*ac(kpm(kf,ifm)+1,if,jf))  &
                 / dn0f(kf,if,jf)
           enddo
        enddo
     endif

     if(iand(ibcon,2).ne.0) then
        if=iz
        do jf=ja,jz
           do kf=1,m1f
              bx(kf,jf,2)=  &
                  (ek1(kf,ifm)*ac(kpm(kf,ifm)-1,if,jf)  &
                 + ek2(kf,ifm)*ac(kpm(kf,ifm)  ,if,jf)  &
                 + ek3(kf,ifm)*ac(kpm(kf,ifm)+1,if,jf))  &
                 / dn0f(kf,if,jf)
           enddo
        enddo
     endif

     if(iand(ibcon,4).ne.0) then
        jf=ja
        do if=ia,iz
           do kf=1,m1f
              by(kf,if,1)=  &
                  (ek1(kf,ifm)*ac(kpm(kf,ifm)-1,if,jf)  &
                 + ek2(kf,ifm)*ac(kpm(kf,ifm)  ,if,jf)  &
                 + ek3(kf,ifm)*ac(kpm(kf,ifm)+1,if,jf))  &
                 / dn0f(kf,if,jf)
           enddo
        enddo
     endif

     if(iand(ibcon,8).ne.0) then
        jf=jz
        do if=ia,iz
           do kf=1,m1f
              by(kf,if,2)=  &
                  (ek1(kf,ifm)*ac(kpm(kf,ifm)-1,if,jf)  &
                 + ek2(kf,ifm)*ac(kpm(kf,ifm)  ,if,jf)  &
                 + ek3(kf,ifm)*ac(kpm(kf,ifm)+1,if,jf))  &
                 / dn0f(kf,if,jf)
           enddo
        enddo
     endif

     if(nstbot.eq.0)then
        kf=1
        do jf=ja,jz
           do if=ia,iz
              bz(if,jf,1)=  &
                  (ek1(kf,ifm)*ac(kpm(kf,ifm)-1,if,jf)  &
                 + ek2(kf,ifm)*ac(kpm(kf,ifm)  ,if,jf)  &
                 + ek3(kf,ifm)*ac(kpm(kf,ifm)+1,if,jf))  &
                 / dn0f(kf,if,jf)
           enddo
        enddo
     endif

     if(nsttop.eq.0)then
        kf=m1f
        do jf=ja,jz
           do if=ia,iz
              bz(if,jf,2)=  &
                  (ek1(kf,ifm)*ac(kpm(kf,ifm)-1,if,jf)  &
                 + ek2(kf,ifm)*ac(kpm(kf,ifm)  ,if,jf)  &
                 + ek3(kf,ifm)*ac(kpm(kf,ifm)+1,if,jf))  &
                 / dn0f(kf,if,jf)
           enddo
        enddo
     endif

  endif
  return
end subroutine par_bintp



