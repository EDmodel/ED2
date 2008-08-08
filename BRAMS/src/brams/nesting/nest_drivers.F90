!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################


subroutine nstbdriv()

  use mem_tend, only: &
       tend

  use var_tables, only: &
       num_scalar,      &
       scalar_tab

  use mem_basic, only: &
       basic_g

  use mem_nestb, only: &
       nbounds

  use node_mod, only: &
       mxp,           &
       myp,           &
       mzp,           &
       ia,            &
       iz,            &
       ja,            &
       jz,            &
       ibcon

  use mem_grid, only: &
       dtlv,          &
       dtlt,          &
       nndtrat,       &
       isstp,         &
       ngrid,         &
       nstbot,        &
       nsttop,        &
       jdim

  implicit none


  integer :: n
  real :: tymeinvv,tymeinvs
#if USE_INTERF
  interface 
    subroutine nstbtnd(m1,m2,m3,ia,iz,ja,jz,ibcon,scp,sct,bx,by,&
                       bz,vnam,tymeinv,nstbot,nsttop,jdim)
      integer :: m1,m2,m3,ia,iz,ja,jz,ibcon
      real, dimension(m1,m2,m3) :: scp,sct
      real, dimension(m1,m3,2) :: bx,by,bz
      character(*) :: vnam
      real :: tymeinv
      integer :: nstbot,nsttop,jdim
    end subroutine nstbtnd
  end interface 
#endif
  tymeinvv = 1.0 / (dtlv * 0.5 * float(nndtrat(ngrid)+2-isstp))
  tymeinvs = 1.0 / (dtlt * float(nndtrat(ngrid)+1-isstp))

  call nstbtnd(mzp,mxp,myp,ia,iz,ja,jz,ibcon   &
       ,basic_g(ngrid)%up,tend%ut              &
       ,nbounds(ngrid)%bux,nbounds(ngrid)%buy  &
       ,nbounds(ngrid)%buz                     &
       ,'u',tymeinvv,nstbot,nsttop,jdim)

  call nstbtnd(mzp,mxp,myp,ia,iz,ja,jz,ibcon   &
       ,basic_g(ngrid)%vp,tend%vt              &
       ,nbounds(ngrid)%bvx,nbounds(ngrid)%bvy  &
       ,nbounds(ngrid)%bvz                     &
       ,'v',tymeinvv,nstbot,nsttop,jdim)

  call nstbtnd(mzp,mxp,myp,ia,iz,ja,jz,ibcon   &
       ,basic_g(ngrid)%wp,tend%wt              &
       ,nbounds(ngrid)%bwx,nbounds(ngrid)%bwy  &
       ,nbounds(ngrid)%bwz                     &
       ,'w',tymeinvv,nstbot,nsttop,jdim)

  call nstbtnd(mzp,mxp,myp,ia,iz,ja,jz,ibcon   &
       ,basic_g(ngrid)%pp,tend%pt              &
       ,nbounds(ngrid)%bpx,nbounds(ngrid)%bpy  &
       ,nbounds(ngrid)%bpz                     &
       ,'p',tymeinvv,nstbot,nsttop,jdim)

  do n = 1,num_scalar(ngrid)
     call nstbtnd(mzp,mxp,myp,ia,iz,ja,jz,ibcon                     &
          ,scalar_tab(n,ngrid)%var_p,scalar_tab(n,ngrid)%var_t      &
          ,nbounds(ngrid)%bsx(:,:,:,n),nbounds(ngrid)%bsy(:,:,:,n)  &
          ,nbounds(ngrid)%bsz(:,:,:,n)                              &
          ,'t',tymeinvs,nstbot,nsttop,jdim)
  enddo
  return
end subroutine nstbdriv

!******************************************************************************
!******************************************************************************

!from initlz:
!      call prgintrp(nnzp(icm),nnxp(icm),nnyp(icm),nnzp(icm),nnxp(icm)
!     +   ,nnyp(icm),0,0,ifm,1,b)
!from model:
!      call prgintrp(nnzp(icm),nnxp(icm),nnyp(icm),mmzp(icm),mmxp(icm)
!     +   ,mmyp(icm),mi0(icm),mj0(icm),ifm,0,b)
!from parallel model:
!      call prgintrp(nnzp(icm),nnxp(icm),nnyp(icm),mmzp(icm),mmxp(icm)
!     +   ,mmyp(icm),mi0(icm),mj0(icm),ifm,-1,b)

subroutine prgintrp(n1c,n2c,n3c,m1c,m2c,m3c,i0c,j0c,ifm,initflg,mynum)

  use var_tables, only: &
       num_scalar,      &
       scalar_tab

  use mem_scratch, only: &
       scratch

  use mem_basic, only: &
       basic_g

  use mem_grid, only: &
       nnxp,          &
       nnyp,          &
       nnzp,          &
       maxnxp,        &
       maxnyp,        &
       maxnzp,        &
       grid_g,        &
       nxtnest,       &
       nnstbot,       &
       nnsttop,       &
       icorflg,       &
       jdim

  use mem_nestb, only: &
       nbounds

  implicit none

  integer, intent(in) :: n1c       ! unused!!!
  integer, intent(in) :: n2c
  integer, intent(in) :: n3c
  integer, intent(in) :: m1c
  integer, intent(in) :: m2c
  integer, intent(in) :: m3c
  integer, intent(in) :: i0c       ! unused!!!
  integer, intent(in) :: j0c       ! unused!!!
  integer, intent(in) :: initflg
  integer, intent(in) :: ifm
  integer, intent(in) :: mynum     ! unused!!!

  integer :: icm
  integer :: nf
  integer :: nc
  character(len=10) :: c0, c1
  character(len=*), parameter :: h="**(prgintrp)**"

  !     Temporarily fill VT2DA with interpolated topography from coarser grid

  icm = nxtnest(ifm)
  if (icm .eq. 0) return
  if (initflg .eq. 1) then
     call fillscr(1,maxnxp,maxnyp,1,n2c,n3c,1,1  &
          ,scratch%scr1(1),grid_g(icm)%topt(1,1))
     call eintp(scratch%scr1(1),scratch%scr2(1)  &
          ,1,maxnxp,maxnyp,1,nnxp(ifm),nnyp(ifm),ifm,2,'t',0,0)
     call fillvar(1,maxnxp,maxnyp,1,nnxp(ifm),nnyp(ifm),1,1  &
          ,scratch%scr2(1),scratch%vt2da(1))
  endif

  !     Interpolate atmospheric variables

  call fmint3(m1c,m2c,m3c,nnzp(ifm),nnxp(ifm),nnyp(ifm),maxnzp,maxnxp,maxnyp  &
       ,ifm,icm,nnstbot(ifm),nnsttop(ifm),jdim,initflg,1,1,'u'  &
       ,basic_g(icm)%uc(1,1,1),basic_g(ifm)%uc(1,1,1)  &
       ,basic_g(icm)%dn0u(1,1,1),basic_g(ifm)%dn0u(1,1,1)  &
       ,scratch%scr1(1),scratch%scr2(1)  &
       ,grid_g(ifm)%topt(1,1),scratch%vt2da(1)  &
       ,nbounds(ifm)%bux(1,1,1),nbounds(ifm)%buy(1,1,1)  &
       ,nbounds(ifm)%buz(1,1,1),mynum)

  if (jdim .eq. 1 .or. icorflg .eq. 1)  &
       call fmint3(m1c,m2c,m3c,nnzp(ifm),nnxp(ifm),nnyp(ifm),maxnzp,maxnxp,maxnyp  &
       ,ifm,icm,nnstbot(ifm),nnsttop(ifm),jdim,initflg,1,1,'v'  &
       ,basic_g(icm)%vc(1,1,1),basic_g(ifm)%vc(1,1,1)  &
       ,basic_g(icm)%dn0v(1,1,1),basic_g(ifm)%dn0v(1,1,1)  &
       ,scratch%scr1(1),scratch%scr2(1)  &
       ,grid_g(ifm)%topt(1,1),scratch%vt2da(1)  &
       ,nbounds(ifm)%bvx(1,1,1),nbounds(ifm)%bvy(1,1,1)  &
       ,nbounds(ifm)%bvz(1,1,1),mynum)

  call fmint3(m1c,m2c,m3c,nnzp(ifm),nnxp(ifm),nnyp(ifm),maxnzp,maxnxp,maxnyp  &
       ,ifm,icm,nnstbot(ifm),nnsttop(ifm),jdim,initflg,1,1,'w'  &
       ,basic_g(icm)%wc(1,1,1),basic_g(ifm)%wc(1,1,1)  &
       ,basic_g(icm)%dn0(1,1,1),basic_g(ifm)%dn0(1,1,1)  &
       ,scratch%scr1(1),scratch%scr2(1)  &
       ,grid_g(ifm)%topt(1,1),scratch%vt2da(1)  &
       ,nbounds(ifm)%bwx(1,1,1),nbounds(ifm)%bwy(1,1,1)  &
       ,nbounds(ifm)%bwz(1,1,1),mynum)

  call fmint3(m1c,m2c,m3c,nnzp(ifm),nnxp(ifm),nnyp(ifm),maxnzp,maxnxp,maxnyp  &
       ,ifm,icm,nnstbot(ifm),nnsttop(ifm),jdim,initflg,0,1,'t'  &
       ,basic_g(icm)%pc(1,1,1),basic_g(ifm)%pc(1,1,1)  &
       ,basic_g(icm)%dn0(1,1,1),basic_g(ifm)%dn0(1,1,1)  &
       ,scratch%scr1(1),scratch%scr2(1)  &
       ,grid_g(ifm)%topt(1,1),scratch%vt2da(1)  &
       ,nbounds(ifm)%bpx(1,1,1),nbounds(ifm)%bpy(1,1,1)  &
       ,nbounds(ifm)%bpz(1,1,1),mynum)

  do nf = 1,num_scalar(ifm)
     do nc = 1,num_scalar(icm)
        if (scalar_tab(nf,ifm)%name == scalar_tab(nc,icm)%name) then
           call fmint3(m1c,m2c,m3c,nnzp(ifm),nnxp(ifm),nnyp(ifm)  &
                ,maxnzp,maxnxp,maxnyp,ifm,icm  &
                ,nnstbot(ifm),nnsttop(ifm),jdim,initflg,1,1,'t'  &
                ,scalar_tab(nc,icm)%var_p &
                ,scalar_tab(nf,ifm)%var_p  &
                ,basic_g(icm)%dn0(1,1,1),basic_g(ifm)%dn0(1,1,1)  &
                ,scratch%scr1(1),scratch%scr2(1)  &
                ,grid_g(ifm)%topt(1,1),scratch%vt2da(1)  &
                ,nbounds(ifm)%bsx(1,1,1,nf),nbounds(ifm)%bsy(1,1,1,nf)  &
                ,nbounds(ifm)%bsz(1,1,1,nf),mynum)
        endif
     enddo
  enddo
end subroutine prgintrp

!******************************************************************************

subroutine nstfeed(ifm,icm)

  use var_tables, only: &
       num_scalar,      &
       scalar_tab

  use mem_basic, only: &
       basic_g

  use mem_scratch, only: &
       scratch

  use node_mod, only: &
       mia,           &
       miz,           &
       mja,           &
       mjz,           &
       ibcon

  use mem_grid, only: &
       nnzp,          &
       nnxp,          &
       nnyp,          &
       nnstbot,       &
       nnsttop

  implicit none

  integer :: ifm,icm

  integer :: nf,nc
  character(len=10) :: c0, c1
  character(len=*), parameter :: h="**(nstfeed)**"

  !     feed back the finer mesh to the coarser mesh.

  if (icm .eq. 0) return

  call fdback(basic_g(icm)%uc   (1,1,1)  ,basic_g(ifm)%uc   (1,1,1)  &
       ,basic_g(icm)%dn0u (1,1,1)  ,basic_g(ifm)%dn0u (1,1,1)  &
       ,nnzp(icm),nnxp(icm),nnyp(icm),nnzp(ifm),nnxp(ifm),nnyp(ifm)    &
       ,ifm,'u',scratch%scr1(1))
  if (nnstbot(icm) .eq. 1) then
     call botset(nnzp(icm),nnxp(icm),nnyp(icm)  &
          ,mia(icm),miz(icm),mja(icm),mjz(icm),ibcon  &
          ,basic_g(icm)%uc(1,1,1),'U')
  end if
  if (nnsttop(icm) .eq. 1) then
     call topset(nnzp(icm),nnxp(icm),nnyp(icm)  &
          ,mia(icm),miz(icm),mja(icm),mjz(icm),ibcon  &
          ,basic_g(icm)%uc(1,1,1),basic_g(icm)%uc(1,1,1),'U')
  end if

  call fdback(basic_g(icm)%vc   (1,1,1)  ,basic_g(ifm)%vc   (1,1,1)  &
       ,basic_g(icm)%dn0v (1,1,1)  ,basic_g(ifm)%dn0v (1,1,1)  &
       ,nnzp(icm),nnxp(icm),nnyp(icm),nnzp(ifm),nnxp(ifm),nnyp(ifm)    &
       ,ifm,'v',scratch%scr1(1))
  if (nnstbot(icm) .eq. 1) then
     call botset(nnzp(icm),nnxp(icm),nnyp(icm)  &
          ,mia(icm),miz(icm),mja(icm),mjz(icm),ibcon  &
          ,basic_g(icm)%vc(1,1,1),'V')
  end if
  if (nnsttop(icm) .eq. 1) then
     call topset(nnzp(icm),nnxp(icm),nnyp(icm)  &
          ,mia(icm),miz(icm),mja(icm),mjz(icm),ibcon  &
          ,basic_g(icm)%vc(1,1,1),basic_g(icm)%vc(1,1,1),'V')
  end if

  call fdback(basic_g(icm)%wc  (1,1,1)  ,basic_g(ifm)%wc  (1,1,1)  &
       ,basic_g(icm)%dn0 (1,1,1)  ,basic_g(ifm)%dn0 (1,1,1)  &
       ,nnzp(icm),nnxp(icm),nnyp(icm),nnzp(ifm),nnxp(ifm),nnyp(ifm)  &
       ,ifm,'w',scratch%scr1(1))

  call fdback(basic_g(icm)%pc  (1,1,1)  ,basic_g(ifm)%pc   (1,1,1)  &
       ,basic_g(icm)%dn0 (1,1,1)  ,basic_g(ifm)%dn0  (1,1,1)  &
       ,nnzp(icm),nnxp(icm),nnyp(icm),nnzp(ifm),nnxp(ifm),nnyp(ifm)   &
       ,ifm,'p',scratch%scr1(1))
  if (nnstbot(icm) .eq. 1) then
     call botset(nnzp(icm),nnxp(icm),nnyp(icm)  &
          ,mia(icm),miz(icm),mja(icm),mjz(icm),ibcon  &
          ,basic_g(icm)%pc(1,1,1),'P')
  end if
  if (nnsttop(icm) .eq. 1) then
     call topset(nnzp(icm),nnxp(icm),nnyp(icm)  &
          ,mia(icm),miz(icm),mja(icm),mjz(icm),ibcon  &
          ,basic_g(icm)%pc(1,1,1),basic_g(icm)%pc(1,1,1),'P')
  end if

  do nf = 1,num_scalar(ifm)
     do nc = 1,num_scalar(icm)
        if (scalar_tab(nf,ifm)%name == scalar_tab(nc,icm)%name) then
           call fdback(scalar_tab(nc,icm)%var_p,scalar_tab(nf,ifm)%var_p  &
                ,basic_g(icm)%dn0(1,1,1),basic_g(ifm)%dn0(1,1,1)  &
                ,nnzp(icm),nnxp(icm),nnyp(icm),nnzp(ifm),nnxp(ifm)  &
                ,nnyp(ifm),ifm,'t',scratch%scr1(1))
           if (nnstbot(icm) .eq. 1) then
              call botset(nnzp(icm),nnxp(icm),nnyp(icm)  &
                   ,mia(icm),miz(icm),mja(icm),mjz(icm),ibcon,scalar_tab(nc,icm)%var_p,'T')
           end if
           if (nnsttop(icm) .eq. 1) then
              call topset(nnzp(icm),nnxp(icm),nnyp(icm)                  &
                   ,mia(icm),miz(icm),mja(icm),mjz(icm),ibcon            &
                   ,scalar_tab(nc,icm)%var_p,scalar_tab(nc,icm)%var_p,'T')
           end if
        endif
     enddo
  enddo
end subroutine nstfeed

