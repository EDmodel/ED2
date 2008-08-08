!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

! MLO - This part used to be in micro, but Harrington radiation needs to have 
!       the tables filled before the microphysics scheme is called by the 1st 
!       time
subroutine micro_1st()
  implicit none
  logical, save :: first=.true.
  
  if (first) then
    first=.false.  
    call micinit()
    call make_autotab()
    call haznuc()
    call tabmelt()
    call tabhab()
  end if
  return
end subroutine micro_1st

!     *****************************************************************

subroutine micro()

  use mem_basic, only : &
       basic_g            ! INTENT(INOUT)

  use mem_micro, only:  &
       micro_g            ! INTENT(INOUT)

  use mem_grid, only:   &
       ngrids,          & ! INTENT(IN)
       ngrid,           & ! INTENT(IN)
       zm,              & ! INTENT(IN)
       dzt,             & ! INTENT(IN)
       dtlt,            & ! INTENT(IN)
       jdim,            & ! INTENT(IN)
       maxnzp,          & ! INTENT(IN)
       time,            & ! INTENT(IN)
       dtlongn,         & ! INTENT(IN)
       zt,              & ! INTENT(IN)
       itimea,          & ! INTENT(IN)
       if_adap,         & ! INTENT(IN)
       grid_g             ! INTENT(IN)

  use mem_radiate, only:&
       radiate_g          ! INTENT(INOUT)

  use node_mod, only :  &
       mmzp,            & ! INTENT(IN)
       mzp,             & ! INTENT(IN)
       mxp,             & ! INTENT(IN)
       myp,             & ! INTENT(IN)
       ja,              & ! INTENT(IN)
       jz,              & ! INTENT(IN)
       ia,              & ! INTENT(IN)
       iz,              & ! INTENT(IN)
       mynum              ! INTENT(IN)

  use micphys, only:    &
       rx, &
       maxgrds,         & ! INTENT(IN)
       level,           & ! INTENT(IN)
       nhcat,           & ! INTENT(IN)
       ch3,             & ! INTENT(OUT)
       pwvt,            & ! INTENT(IN)
       pwmasi,          & ! INTENT(IN)
       cdp1,            & ! INTENT(OUT)
       pwvtmasi,        & ! INTENT(OUT)
       ch2,             & ! INTENT(OUT)
       dispemb1,        & ! INTENT(IN)
       dispemb0,        & ! INTENT(IN)
       icloud,          & ! INTENT(IN)
       irain,           & ! INTENT(IN)
       ipris,           & ! INTENT(IN)
       isnow,           & ! INTENT(IN)
       iaggr,           & ! INTENT(IN)
       igraup,          & ! INTENT(IN)
       ihail              ! INTENT(IN)
       
  use mem_leaf, only : isfcl ! INTENT(IN)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Commented for now, it will be uncommented soon
!  use mem_ed,   only : dtlsm ! INTENT(IN)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  implicit none

  ! Local Variables:

  integer :: nembfall,maxkfall,ngr,lhcat,i,j,k
  integer, dimension(10)  :: k1,k2,k3
  integer, save :: ncall = 0
  integer, save, dimension(15)  :: ncall2g = 0
  real :: dtlti ! rshort, cosz, rlongup, rlong, albedt ! Not used
  type pcp_tab_type
     real, pointer, dimension(:,:,:,:) :: pcpfillc,pcpfillr
     real, pointer, dimension(:,:,:)   :: sfcpcp
  end type pcp_tab_type
  type (pcp_tab_type), save :: pcp_tab(maxgrds)


  if (level .ne. 3) return
  nembfall = 20
  maxkfall = 4

  if(ncall == 0) then
     ncall = 1

     do ngr = 1,ngrids
        allocate (pcp_tab(ngr)%pcpfillc(mmzp(ngr),maxkfall,nembfall,nhcat))
        allocate (pcp_tab(ngr)%pcpfillr(mmzp(ngr),maxkfall,nembfall,nhcat))
        allocate (pcp_tab(ngr)%sfcpcp(maxkfall,nembfall,nhcat))
     enddo

!     call micro_1st()

     do lhcat = 1,nhcat
        ch3(lhcat) = pwvt(lhcat) * pwmasi(lhcat)
        cdp1(lhcat) = pwmasi(lhcat) * (1.5 + .5 * pwvt(lhcat))
        pwvtmasi(lhcat) = pwvt(lhcat) * pwmasi(lhcat)
     enddo
  endif
  if (ncall2g(ngrid) /= 5) then
     ncall2g(ngrid) = 5

     call mksedim_tab(mzp,mxp,myp,ngrid,nembfall,maxkfall,zm,dzt  &
          ,pcp_tab(ngrid)%pcpfillc(1,1,1,1),pcp_tab(ngrid)%pcpfillr(1,1,1,1)  &
	  ,pcp_tab(ngrid)%sfcpcp(1,1,1))

     do lhcat = 1,nhcat
        ch2(lhcat,ngrid) = float(nembfall-1) &
             / log10(dispemb1(lhcat,ngrid) / dispemb0(lhcat,ngrid))
     enddo

     call homfrzcl(dtlt,ngrid)
  endif

  call each_call(mzp,dtlt)
  dtlti = 1. / dtlt
  ngr = ngrid

  do j = ja,jz
     do i = ia,iz

        call range_check(mzp,k1,k2,k3,i,j,grid_g(ngr)%flpw(i,j),micro_g(ngr))

        call mcphys(mzp,k1,k2,k3,i,j,ngrid,jdim,maxnzp             &
             ,nembfall,maxkfall,mynum,dtlt,dtlti,time,zm,dzt                 &
             ,zt,itimea,radiate_g(ngr)                &
             ,basic_g(ngr)%thp     (1,i,j)   ,basic_g(ngr)%theta   (1,i,j)   &
             ,basic_g(ngr)%pp      (1,i,j)   ,basic_g(ngr)%rtp     (1,i,j)   &
             ,basic_g(ngr)%rv      (1,i,j)   ,basic_g(ngr)%wp      (1,i,j)   &
             ,basic_g(ngr)%dn0     (1,i,j)   ,basic_g(ngr)%pi0     (1,i,j)   &
             ,grid_g(ngr)%rtgt     (i,j)     ,grid_g(ngr)%flpw      (i,j)     &
             ,micro_g(ngr)%pcpg    (i,j)     ,micro_g(ngr)%qpcpg   (i,j)     &
             ,micro_g(ngr)%dpcpg   (i,j)     ,pcp_tab(ngr)%pcpfillc(1,1,1,1) &
             ,pcp_tab(ngr)%pcpfillr(1,1,1,1) ,pcp_tab(ngr)%sfcpcp  (1,1,1)   &
             ,grid_g(ngr)%glat     (i,j)     ,grid_g(ngr)%topt     (i,j)     &
             ,if_adap     )

        call copyback(mzp,k1,k2,k3,grid_g(ngr)%flpw(i,j),i,j,micro_g(ngr))

     enddo
  enddo

  ! Computing the average pcpg/qpcpg/dpcpg for ED. This subroutine is in src/mixed/met_update.f90
  if (isfcl == 5) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! IT WILL BE UNCOMMENTED SOON
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     call average_micro_rates_for_ed(mxp,myp,ia,iz,ja,jz,dtlongn(ngr),dtlsm                        &
!                                    ,micro_g(ngr)%pcpg       (1,1) , micro_g(ngr)%qpcpg      (1,1) &
!                                    ,micro_g(ngr)%dpcpg      (1,1) , micro_g(ngr)%mean_pcpg  (1,1) &
!                                    ,micro_g(ngr)%mean_qpcpg (1,1) , micro_g(ngr)%mean_dpcpg (1,1) )
  end if

  return
end subroutine micro

!     *****************************************************************

subroutine mcphys(m1,k1,k2,k3,i,j,ngr,jdim,maxnzp  &
     ,nembfall,maxkfall,mynum,dtlt,dtlti,time,zm,dzt,zt  &
     ,itime1,radiate  &
     ,thp,theta,pp,rtp,rv,wp,dn0,pi0  &
     ,rtgt,flpw,pcpg,qpcpg,dpcpg  &
     ,pcpfillc,pcpfillr,sfcpcp  &
     ,glat,topt,if_adap)

  use mem_radiate, only: &
       radiate_vars,     & ! INTENT(IN)
       iswrtyp,          & ! INTENT(IN)
       ilwrtyp,          & ! INTENT(IN)
       radfrq,           & ! INTENT(IN)
       nadd_rad

  use micphys, only:     &
       nhcat,            & ! INTENT(IN)
       jnmb,             & ! INTENT(IN)
       ncat,             & ! INTENT(IN)
       scrmic1,          & ! INTENT(INOUT)
       scrmic2,          & ! INTENT(INOUT)
       tairc,            & ! INTENT(OUT)
       ch1,              & ! INTENT(OUT)
       cfvt,             & ! INTENT(IN)
       scrmic3,rx          ! INTENT(OUT)
       

  implicit none

  ! Arguments:
  integer, intent(in)                   :: m1, i, j, ngr, jdim, maxnzp, &
       nembfall, maxkfall, mynum
  integer, dimension(10), intent(inout) :: k1, k2, k3
  real, intent(in)                      :: dtlt, dtlti
  real(kind=8) :: time

  real, dimension(m1), intent(in)       :: zm
  real, dimension(m1), intent(in)       :: dzt ! Not used
  integer, intent(in)                   :: itime1 ! Not used
  type (radiate_vars), intent(inout)    :: radiate
  real, dimension(m1), intent(inout)    :: thp
  real, dimension(m1), intent(inout)    :: theta
  real, dimension(m1), intent(in)       :: pp
  real, dimension(m1), intent(inout)    :: rtp
  real, dimension(m1), intent(out)      :: rv
  real, dimension(m1), intent(in)       :: wp
  real, dimension(m1), intent(in)       :: dn0
  real, dimension(m1), intent(in)       :: pi0  
  real, intent(in)                      :: rtgt
  real, intent(in)                      :: flpw
  real, intent(out)                     :: pcpg, qpcpg, dpcpg
  integer, intent(in)                   :: if_adap

  ! Variables needed for Harrington radiation scheme
  real, dimension(m1,maxkfall,nembfall,nhcat), intent(in) :: pcpfillc, pcpfillr
  real, dimension(maxkfall,nembfall,nhcat), intent(in)    :: sfcpcp
  real, intent(in)                                        :: glat, topt
  real, dimension(m1), intent(in)                         :: zt

  ! Local Variables:
  integer :: k,jflag,lcat,icv,icx,mc1,mc2,mc3,mc4  &
       ,mcat,k1cnuc,k2cnuc,k1pnuc,k2pnuc,lhcat,ilpw
  integer, dimension(7)   :: mcats,mivap,mix02
  real, dimension(7)      :: dpcp0
  integer, dimension(9,4) :: mcat1
  integer, dimension(7,2) :: mcat2
  integer, dimension(4)   :: mcat33

  data mcats /0,3,0,0,6,7,10/
  data mcat1 /3,3,3,4,4,4,5,5,6  &
       ,5,6,7,5,6,7,6,7,7  &
       ,5,6,7,5,6,7,6,7,7  &
       ,4,7,8,5,7,8,7,8,8/
  data mcat2 /0,0,0,6,6,7,7  &
       ,0,0,0,2,2,9,9/
  data mcat33 /0,0,4,5/
  data mivap /1,3,4,5,2,6,7/
  data mix02 /3,1,4,5,6,7,2/
  data dpcp0 /.001,.001,.010,.010,.010,.003,.001/

  save

  call thrmstr(m1,k1,k2,flpw,pp(1),thp(1),theta(1),pi0(1),rtp(1),rv(1),i,j)

  call each_column(m1,k1,k2,i,j,flpw,rv(1),dn0(1))

  ! Diagnose hydrometeor mean mass emb, and if necessary, number concentration.

  jflag = 1

  do lcat = 1,7
     if (jnmb(lcat) .ge. 1) then
        call enemb(m1,k1(lcat),k2(lcat),lcat,jflag,dn0(1),i,j)
     endif
  enddo

! MLO - There used to have a call for radiation with Harrington. But let's leave the
!       radiation driver to drive the radiation and the microphysics to drive the 
!       microphysics...

  do lcat = 1,7
     if (jnmb(lcat) .ge. 1) then
        call diffprep(m1,lcat,k1(lcat),k2(lcat),rv(1),dn0(1),i,j,mynum)
     endif
  enddo

  call vapdiff(m1,k1(10),k2(10),rv(1),i,j,mynum)

  do icv = 1,7
     lcat = mivap(icv)

     if (jnmb(lcat) .ge. 1) then
        call vapflux(m1,lcat,i,j,mynum,k1(lcat),k2(lcat),dn0(1),rv(1))
     endif
  enddo

  if (jnmb(4) .ge. 1) then
     call psxfer (m1,min(k1(3),k1(4)),max(k2(3),k2(4)),dn0(1),i,j)
  endif

  jflag = 2
  do lcat = 1,7
     if (jnmb(lcat) .ge. 1) then
        call enemb(m1,k1(lcat),k2(lcat),lcat,jflag,dn0(1),i,j)
        call getict(k1(lcat),k2(lcat),lcat,i,j,mynum)
     endif
  enddo

  call newtemp(m1,k1(10),k2(10),rv(1),theta(1),i,j)

  if (jnmb(2) .ge. 1) then
     call auto_accret(m1,k1(1),k2(1),dn0(1),dtlt,i,j)
  endif

  call effxy(m1,k1,k2,i,j)

  ! Self collection of rain, aggregates, graupel, hail:  number change only

  do lcat = 2,7

     if (lcat /= 3 .and. lcat /= 4) then
       mc1 = mcats(lcat)
       if (jnmb(lcat) >= 5) then
          call cols (m1,lcat,mc1,k1(lcat),k2(lcat),i,j)
       endif
     end if
  enddo

  ! Self collection of pristine ice, snow

  do lcat = 3,4
     mc1 = mcat33(lcat)
     if (jnmb(lcat) .ge. 1 .and. jnmb(5) .ge. 1) then
        call col3344 (m1,lcat,5,mc1,k1(lcat),k2(lcat),i,j)
     endif
  enddo

  ! Collection between pristine ice and snow

  if (jnmb(5) .ge. 1) then
     call col3443 (m1,3,4,5,max(k1(3),k1(4)),min(k2(3),k2(4)),i,j)
  endif

  ! Ice-ice collisions

  do icx = 1,9
     mc1 = mcat1(icx,1)
     mc2 = mcat1(icx,2)
     mc3 = mcat1(icx,3)
     mc4 = mcat1(icx,4)

     if (jnmb(mc1) .ge. 1 .and. jnmb(mc3) .ge. 1) then
        call col1 (m1,mc1,mc2,mc3,mc4,max(k1(mc1),k1(mc2))  &
             ,min(k2(mc1),k2(mc2)),i,j)
     endif
  enddo

  ! Ice-cloud collisions

  do lcat = 4,7
     mc1 = mcat2(lcat,1)
     mc2 = mcat2(lcat,2)

     if (jnmb(lcat) .ge. 1 .and. jnmb(mc1) .ge. 1) then
        call col2 (m1,1,lcat,mc1,mc2,max(k1(1),k1(lcat)),min(k2(1),k2(lcat))  &
             ,dn0(1),dtlt,i,j)
     endif
  enddo

  ! Ice-rain collisions

  do lcat = 3,7
     if (jnmb(lcat) .ge. 1 .and. jnmb(7) .ge. 1) then
        call col3 (m1,2,lcat,7,max(k1(2),k1(lcat)),min(k2(2),k2(lcat)),i,j)
     endif
  enddo

  call colxfers(m1,k1,k2,i,j,scrmic1,scrmic2)

  do mcat = 1,7
     lcat = mix02(mcat)
     if (jnmb(lcat) .ge. 1) call x02(m1,k1,k2,lcat,dn0(1),i,j)
  enddo

  if (jnmb(1) .ge. 1) then
     call cldnuc(m1,k1cnuc,k2cnuc,flpw,rv(1),wp(1),i,j)
  endif

  k1(1) = min(k1(1),k1cnuc)
  k2(1) = max(k2(1),k2cnuc)
  k3(1) = max(k2(1),k3(1))

  if (jnmb(1) .ge. 1) then
     call c03(m1,k1(1),k2(1),1,dn0(1),i,j)
  endif

  if (jnmb(3) .ge. 1) then
     call icenuc(m1,k1(1),k2(1),k1pnuc,k2pnuc,flpw,ngr,rv(1),dn0(1),dtlt,i,j)
  endif

  k1(3) = min(k1(3),k1pnuc)
  k2(3) = max(k2(3),k2pnuc)
  k3(3) = max(k2(3),k3(3))

  do lcat = 3,1,-2
     if (jnmb(lcat) .ge. 1) then
        call pc03(m1,k1(lcat),k2(lcat),lcat,dn0(1),i,j)
     endif
  enddo

  !  Zero out precip arrays.

  pcpg = 0.
  qpcpg = 0.
  dpcpg = 0.

  ! tairc is used here to accumulate changes to thp from sedim

  ilpw = nint(flpw)

  do k = ilpw,m1
     tairc(k) = 0.
  enddo

  do lhcat = 2,nhcat
     ch1(lhcat) = dtlt * cfvt(lhcat) / rtgt
  enddo

  do lcat = 2,7
     if (jnmb(lcat) .ge. 1) then
        ! New subroutine from RAMS 6.0
        call sedim (m1,lcat,ngr,nembfall,maxkfall  &
             ,k1(lcat),k2(lcat),flpw,i,j  &
             ,rtp(1),thp(1),theta(1),dn0(1),dpcp0(lcat)  &
             ,pcpg,qpcpg,dpcpg,dtlti,scrmic1,scrmic2,scrmic3  &
             ,pcpfillc,pcpfillr,sfcpcp,dzt,if_adap)
     endif
  enddo

  do k = ilpw,m1
     thp(k) = thp(k) + tairc(k)
  enddo

  return
end subroutine mcphys

!******************************************************************************

subroutine copyback(m1,k1,k2,k3,flpw,i,j,micro)

  use mem_micro, only: &
       micro_vars        ! INTENT(IN) ! Only a type structure

  use micphys, only:   &
       jnmb,           & ! INTENT(IN)
       rx,             & ! INTENT(IN)
       cx,             & ! INTENT(IN)
       qx,             & ! INTENT(IN)
       accpx,          & ! INTENT(IN)
       pcprx             ! INTENT(IN)

  implicit none

  ! Arguments:
  integer, intent(in)              :: m1
  type (micro_vars), intent(inout) :: micro
  integer, dimension(10), intent(in) :: k1, k2, k3
  integer, intent(in) :: i,j
  integer :: ilpw
  real, intent(in) :: flpw


  ilpw=nint(flpw)
  
  if (jnmb(1) >= 1) then
     call ae1kmic(ilpw,k3(1),micro%rcp(1,i,j),rx(1,1))
     if (jnmb(1) >= 5) call ae1kmic(ilpw,k3(1),micro%ccp(1,i,j),cx(1,1))
  endif

  if (jnmb(2) >= 1) then
     call ae1kmic(ilpw,k2(10),micro%rrp(1,i,j),rx(1,2))
     call ae1kmic(ilpw,k2(10),micro%q2(1,i,j),qx(1,2))
     micro%accpr(i,j) = micro%accpr(i,j) + accpx(2)
     micro%pcprr(i,j) = pcprx(2)
     if (jnmb(2) >= 5) call ae1kmic(ilpw,k3(1),micro%crp(1,i,j),cx(1,2))
  endif

  if (jnmb(3) >= 1) then
     call ae1kmic(ilpw,k3(3),micro%rpp(1,i,j),rx(1,3))
     micro%accpp(i,j) = micro%accpp(i,j) + accpx(3)
     micro%pcprp(i,j) = pcprx(3)
     if (jnmb(3) >= 5) call ae1kmic(ilpw,k3(3),micro%cpp(1,i,j),cx(1,3))
  endif

  if (jnmb(4) >= 1) then
     call ae1kmic(ilpw,k2(10),micro%rsp(1,i,j),rx(1,4))
     micro%accps(i,j) = micro%accps(i,j) + accpx(4)
     micro%pcprs(i,j) = pcprx(4)
     if (jnmb(4) >= 5) call ae1kmic(ilpw,k2(10),micro%csp(1,i,j),cx(1,4))
  endif

  if (jnmb(5) >= 1) then
     call ae1kmic(ilpw,k2(10),micro%rap(1,i,j),rx(1,5))
     micro%accpa(i,j) = micro%accpa(i,j) + accpx(5)
     micro%pcpra(i,j) = pcprx(5)
     if (jnmb(5) >= 5) call ae1kmic(ilpw,k2(10),micro%cap(1,i,j),cx(1,5))
  endif

  if (jnmb(6) >= 1) then
     call ae1kmic(ilpw,k2(10),micro%rgp(1,i,j),rx(1,6))
     call ae1kmic(ilpw,k2(10),micro%q6(1,i,j),qx(1,6))
     micro%accpg(i,j) = micro%accpg(i,j) + accpx(6)
     micro%pcprg(i,j) = pcprx(6)
     if (jnmb(6) >= 5) call ae1kmic(ilpw,k2(10),micro%cgp(1,i,j),cx(1,6))
  endif

  if (jnmb(7) >= 1) then
     call ae1kmic(ilpw,k2(10),micro%rhp(1,i,j),rx(1,7))
     call ae1kmic(ilpw,k2(10),micro%q7(1,i,j),qx(1,7))
     micro%accph(i,j) = micro%accph(i,j) + accpx(7)
     micro%pcprh(i,j) = pcprx(7)
     if (jnmb(7) >= 5) call ae1kmic(ilpw,k2(10),micro%chp(1,i,j),cx(1,7))
  endif

  return
end subroutine copyback


