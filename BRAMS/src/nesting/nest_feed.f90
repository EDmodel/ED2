!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################


subroutine fdback(ac,af,dc,df,nzc,nxc,nyc,nzf,nxf,nyf,nf,vnam,sumflg)

  use mem_grid, only: &
       nxtnest,       & ! intent(in)
       jdim,          & ! intent(in)
       nstratx,       & ! intent(in)
       nstraty,       & ! intent(in)
       nrz,           & ! intent(in)
       kpm,           & ! intent(in)
       jpm,           & ! intent(in)
       ipm,           & ! intent(in)
       fbcf             ! intent(in)

  implicit none

  ! Arguments:
  integer, intent(in)          :: nzc, nxc, nyc, nzf, nxf, nyf, nf
  character(len=*), intent(in) :: vnam
  real, intent(inout)          :: ac(nzc,nxc,nyc)      !value at Coarse grid
  real, intent(in)             :: dc(nzc,nxc,nyc)      !density at Coarse grid
  real, intent(inout)          :: sumflg(nzc,nxc,nyc)  !Coarse grid
  real, intent(in)             :: af(nzf,nxf,nyf)      !value at Fine grid
  real, intent(in)             :: df(nzf,nxf,nyf)      !density at Fine grid

  ! Local variables:
  integer :: ibeg, jbeg, kbeg, iend, jend, kend, i, j, k, nc, iinc, jinc, &
       kv, ifbcf, ic, jc, kc, if, jf, kf


  nc = nxtnest(nf)
  call azero(nzc*nxc*nyc,sumflg)
  ibeg = 2
  jbeg = 1 + jdim
  kbeg = 2
  iend = nxf - 1
  jend = nyf - jdim
  kend = nzf - 1
  iinc = 1
  jinc = 1
  kv = 0

  if (vnam .eq. 'u') then
     ibeg = 1 + nstratx(nf)
     iend = nxf - 1 - nstratx(nf)
     iinc = nstratx(nf)
     ifbcf = 1
  elseif (vnam .eq. 'v') then
     jbeg = 1 + nstraty(nf) * jdim
     jend = nyf - (1 + nstraty(nf)) * jdim
     jinc = nstraty(nf)
     ifbcf = 2
  elseif (vnam .eq. 'w' .or. vnam .eq. 'terr') then
     if (vnam .eq. 'w') then
        kbeg = 1 + nrz(kpm(2,nf),nf)
        kend = nzf - 1 - nrz(kpm(nzf-1,nf),nf)
        kv = 1
     else
        kbeg = 1
        kend = 1
     endif
     ifbcf = 3
  else
     ifbcf = 4
  endif

!print*,'fdback:',vnam,':',ibeg,iend,ipm(ibeg,nf),ipm(iend,nf)

  kf = kbeg
1 continue
  kc = kpm(kf,nf)
  if (vnam .eq. 'terr') kc = 1

  do jf = jbeg,jend,jinc
     jc = jpm(jf,nf)
     if (vnam .eq. 'p' .or. vnam .eq. 'terr') then
        do if = ibeg,iend,iinc
           ic = ipm(if,nf)
           ac(kc,ic,jc) = ac(kc,ic,jc) * sumflg(kc,ic,jc)  &
                + af(kf,if,jf) * fbcf(kf,nf,ifbcf)
           sumflg(kc,ic,jc) = 1.
        enddo
     elseif (vnam .eq. 'w') then
        do if = ibeg,iend,iinc
           ic = ipm(if,nf)
           ac(kc,ic,jc) = ac(kc,ic,jc) * sumflg(kc,ic,jc)  &
                + af(kf,if,jf) * fbcf(kf,nf,ifbcf)  &
                * (df(kf,if,jf) + df(kf+kv,if,jf))

           sumflg(kc,ic,jc) = 1.
        enddo
     else
        do if = ibeg,iend,iinc
           ic = ipm(if,nf)
           ac(kc,ic,jc) = ac(kc,ic,jc) * sumflg(kc,ic,jc)  &
                + af(kf,if,jf) * fbcf(kf,nf,ifbcf) * df(kf,if,jf)
           ! if(vnam=='u'.and.ic==8.and.jc>=14.and.jc<=18.and.kc==2) then
           !    print*,vnam,':',ac(kc,ic,jc),sumflg(kc,ic,jc),af(kf,if,jf) &
           !          ,fbcf(kf,nf,ifbcf),df(kf,if,jf),kf,if,jf
           ! endif
           sumflg(kc,ic,jc) = 1.
        enddo
     endif
  enddo
  kf = kf + 1
  if (vnam .eq. 'w') kf = kf + nrz(kpm(kf,nf),nf) - 1
  if (kf .le. kend) go to 1

  !      if(vnam.ne.'w'.and.vnam.ne.'terr')then
  !         if(nstbot.eq.1)call botset(nzp,nxp,nyp,ac,vnam)
  !         if(nsttop.eq.1)call topset(nzp,nxp,nyp,ac,ac,vnam)
  !      endif

  do kc = kpm(kbeg,nf),kpm(kend,nf)
     do jc = jpm(jbeg,nf),jpm(jend,nf)

        if (vnam .eq. 'w') then

           do ic = ipm(ibeg,nf),ipm(iend,nf)
              ac(kc,ic,jc) = ac(kc,ic,jc) / (dc(kc,ic,jc) + dc(kc+kv,ic,jc))
           enddo

        elseif (vnam .ne. 'p' .and. vnam .ne. 'terr') then

           do ic = ipm(ibeg,nf),ipm(iend,nf)
              ac(kc,ic,jc) = ac(kc,ic,jc) / dc(kc,ic,jc)
              ! if(vnam=='u'.and.ic==8.and.jc>=14.and.jc<=18.and.kc==2) then
              !  print*,vnam,':',ac(kc,ic-1,jc),ac(kc,ic,jc),dc(kc,ic,jc),ic,jc
              ! endif
           enddo

        endif
     enddo
  enddo

  return
end subroutine fdback

!******************************************************************************

subroutine fdbackp(ivarn,af,buf,mtp,df,dfu,dfv,m1,m2,m3,ifm,icm  &
     ,ia,iz,ja,jz,i0,j0,ibcon,nestratx,nestraty,mynum,i1,i2)

  use mem_grid, only: &
       ipm,           & ! intent(in)
       jpm,           & ! intent(in)
       nrz,           & ! intent(in)
       kpm,           & ! intent(in)
       fbcf             ! intent(in)

  implicit none

  ! Arguments:
  integer, intent(in)  :: ivarn, m1, m2, m3, ifm, icm, ia, iz, ja, jz, &
       i0, j0, ibcon, nestratx, nestraty, mynum
  integer, intent(out) :: mtp
  real, intent(in)     :: af(m1,m2,m3)
  real, intent(in)     :: df(m1,m2,m3)
  real, intent(in)     :: dfu(m1,m2,m3)
  real, intent(in)     :: dfv(m1,m2,m3)
  real, intent(inout)  :: buf(*)
  integer, intent(out) :: i1, i2

  ! Local variables:
  integer :: ibeg, jbeg, kbeg, iend, jend, kend, iinc, jinc, &
       kv, ifbcf, if, jf, kf, indcf

  !     ivarn = variable types 1- u
  !                            2- v
  !                            3- w
  !                            4- p
  !                            5- scalar

  ! Local variables ia, iz, ja, and jz in this subroutine refer not to the
  ! limits of prognostic points but to the limits of points on this fm node
  ! to be averaged to a given cm node.


  ibeg = ia
  jbeg = ja
  kbeg = 2
  iend = iz
  jend = jz
  kend = m1 - 1
  iinc = 1
  jinc = 1
  kv   = 0

  if (ivarn==1) then
     do ibeg = ia+i0,ia+i0+nestratx
        if (mod(ibeg-1,nestratx)==0) exit !go to 11
     enddo
     !11   continue
     i1 = ipm(ibeg,ifm)
     ibeg = ibeg - i0
     iinc = nestratx
     ifbcf = 1
  elseif (ivarn==2) then
     do jbeg = ja+j0,ja+j0+nestraty
        if (mod(jbeg-1,nestraty)==0) exit !go to 10
     enddo
     !10   continue
     i1 = jpm(jbeg,ifm)
     jbeg = jbeg - j0
     jinc = nestraty
     ifbcf = 2
  elseif (ivarn==3) then
     kbeg = 1 + nrz(kpm(2,ifm),ifm)
     kend = m1 - 2
     kv = 1
     ifbcf = 3
  else
     ifbcf = 4
  endif
  
  kf = kbeg
  indcf = 1

  do 

     !1 continue
     if (ivarn==4) then
        do jf = jbeg,jend,jinc
           do if = ibeg,iend,iinc
              buf(indcf) = af(kf,if,jf) * fbcf(kf,ifm,ifbcf)
              indcf = indcf + 1
           enddo
        enddo
     elseif (ivarn==3) then
        do jf = jbeg,jend,jinc
           do if = ibeg,iend,iinc
              buf(indcf) = af(kf,if,jf) * fbcf(kf,ifm,ifbcf)  &
                   * (df(kf,if,jf) + df(kf+kv,if,jf))
              indcf = indcf + 1
           enddo
        enddo
     elseif (ivarn==2) then
        do jf = jbeg,jend,jinc
           i2 = jpm(jf+j0,ifm)
           do if = ibeg,iend,iinc
              buf(indcf) = af(kf,if,jf) * fbcf(kf,ifm,ifbcf)  &
                   * dfv(kf,if,jf)
              indcf = indcf + 1
           enddo
        enddo
     elseif (ivarn==1) then
        do jf = jbeg,jend,jinc
           do if = ibeg,iend,iinc
              i2 = ipm(if+i0,ifm)
              buf(indcf) = af(kf,if,jf) * fbcf(kf,ifm,ifbcf)  &
                   * dfu(kf,if,jf)
              indcf = indcf + 1
           enddo
        enddo
     else
        do jf = jbeg,jend,jinc
           do if = ibeg,iend,iinc
              buf(indcf) = af(kf,if,jf) * fbcf(kf,ifm,ifbcf)  &
                   * df(kf,if,jf)
              indcf = indcf + 1
           enddo
        enddo
     endif
     kf = kf + 1
     if (ivarn==3) kf = kf + nrz(kpm(kf,ifm),ifm) - 1
     
     !if (kf<=kend) go to 1
     if (kf>kend) exit

  enddo

  mtp = indcf - 1

  return
end subroutine fdbackp

