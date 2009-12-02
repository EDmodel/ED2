!     *****************************************************************
subroutine gradr(m1,m2,m3,ia,iz,ja,jz,vc3da,vc3db,dir,gpnt,topt &
  ,xm,xt,ym,yt,zm,zt,deltay,dzm,dzt,vctr1,vctr2,ztop,jdim,ihtran)

  implicit none
  integer         , intent(in)                          :: m1,m2,m3,ia,iz,ja,jz
  real            , intent(in)    , dimension(m1,m2,m3) :: vc3da
  real            , intent(inout) , dimension(m1,m2,m3) :: vc3db
  character(len=*), intent(in)                          :: dir,gpnt
  character(len=6)                                      :: optyp
  real                            , dimension(m1,m2)    :: topt
  real                            , dimension(m1)       :: xm,xt
  real                            , dimension(m2)       :: ym,yt
  real                            , dimension(m3)       :: zm,zt,dzm,dzt
  real                            , dimension(*)        :: vctr1,vctr2
  real                                                  :: deltay,ztop
  integer                                               :: jdim,ihtran
     
   optyp = 'gradnt'
   call rams_grad(m1,m2,m3,ia,iz,ja,jz,vc3da,vc3db,dir,gpnt,optyp,topt &
        ,xm,xt,ym,yt,zm,zt,deltay,dzm,dzt,vctr1,vctr2,ztop,jdim,ihtran)
   return
end subroutine gradr



subroutine divcartr(m1,m2,m3,ia,iz,ja,jz,vc3da,vc3db,dir,gpnt,topt &
        ,xm,xt,ym,yt,zm,zt,deltay,dzm,dzt,vctr1,vctr2,ztop,jdim)

  implicit none
  integer, intent(in) :: m1,m2,m3,ia,iz,ja,jz
  real, intent(in)    :: vc3da(m1,m2,m3)
  real, intent(inout) :: vc3db(m1,m2,m3)
  character(len=*), intent(in) :: dir,gpnt
  character(len=6) :: optyp
  real :: topt(m1,m2),xm(*),xt(*),ym(*),yt(*),zm(*),zt(*),dzm(*),dzt(*), &
          vctr1(*),vctr2(*)
  real :: deltay,ztop
  integer :: jdim,ihtran

     
   optyp = 'divcrt'
   call rams_grad(m1,m2,m3,ia,iz,ja,jz,vc3da,vc3db,dir,gpnt,optyp,topt &
        ,xm,xt,ym,yt,zm,zt,deltay,dzm,dzt,vctr1,vctr2,ztop,jdim,1)
   return
end subroutine

subroutine divstarr(m1,m2,m3,ia,iz,ja,jz,vc3da,vc3db,dir,gpnt,topt &
        ,xm,xt,ym,yt,zm,zt,deltay,dzm,dzt,vctr1,vctr2,ztop,jdim)

  implicit none
  integer, intent(in) :: m1,m2,m3,ia,iz,ja,jz
  real, intent(in)    :: vc3da(m1,m2,m3)
  real, intent(inout) :: vc3db(m1,m2,m3)
  character(len=*), intent(in) :: dir,gpnt
  character(len=6) :: optyp
  real :: topt(m1,m2),xm(*),xt(*),ym(*),yt(*),zm(*),zt(*),dzm(*),dzt(*), &
          vctr1(*),vctr2(*)
  real :: deltay,ztop
  integer :: jdim,ihtran

     
   optyp = 'divsrt'
   call rams_grad(m1,m2,m3,ia,iz,ja,jz,vc3da,vc3db,dir,gpnt,optyp,topt &
        ,xm,xt,ym,yt,zm,zt,deltay,dzm,dzt,vctr1,vctr2,ztop,jdim,1)
   return

end subroutine divstarr

subroutine rams_grad(m1,m2,m3,ia,iz,ja,jz,vc3da,vc3db,dir,gpnt,optyp,topt &
        ,xm,xt,ym,yt,zm,zt,deltay,dzm,dzt,vctr1,vctr2,ztop,jdim,ihtran)
  implicit none
  integer, intent(in) :: m1    &
                       , m2    &
                       , m3    &
                       , ia    &
                       , iz    &
                       , ja    &
                       , jz

  real, intent(in)    :: vc3da(m1,m2,m3)

  real, intent(inout) :: vc3db(m1,m2,m3)

  character(len=*), intent(in) :: dir,gpnt

  character(len=6) :: optyp
  real, dimension(:), allocatable :: e
  real :: topt(m1,m2),xm(*),xt(*),ym(*),yt(*),zm(*),zt(*),dzm(*),dzt(*), &
          vctr1(*),vctr2(*)
  real, dimension(m1) :: ht,hw

  integer::isizexy,isizeall,iglat,iglon,ifmapu,ifmapv,ifmapt,ifmapm,ifmapui,ifmapvi&
         ,ifmapti,ifmapmi,idxu,idxv,idxt,idxm,idyu,idyv                          &
         ,idyt,idym,itopu,itopv,itopm,irtgt,irtgu,irtgv,irtgm,if13u              &
         ,if13v,if13t,if13m,if23u,if23v,if23t,if23m,ioff,ierr                    &
         ,ihtran,jdim,iaddr,jaa,jzz!
  integer, external :: irfree
  real :: polelat,polelon,deltay,ztop


      jaa = ja
      jzz = jz
      if(jdim.eq.0) then
         jaa = 1
         jzz = 1
      endif
      isizexy = m1 * m2
      isizeall = isizexy * 33
      ioff=0
      allocate(e(isizeall))
      iglat   = ioff    + 1
      iglon   = iglat   + isizexy
      ifmapu  = iglon   + isizexy
      ifmapv  = ifmapu  + isizexy
      ifmapt  = ifmapv  + isizexy
      ifmapm  = ifmapt  + isizexy
      ifmapui = ifmapm  + isizexy
      ifmapvi = ifmapui + isizexy
      ifmapti = ifmapvi + isizexy
      ifmapmi = ifmapti + isizexy
      idxu    = ifmapmi + isizexy
      idxv    = idxu    + isizexy
      idxt    = idxv    + isizexy
      idxm    = idxt    + isizexy
      idyu    = idxm    + isizexy
      idyv    = idyu    + isizexy
      idyt    = idyv    + isizexy
      idym    = idyt    + isizexy
      itopu   = idym    + isizexy
      itopv   = itopu   + isizexy
      itopm   = itopv   + isizexy
      irtgt   = itopm   + isizexy
      irtgu   = irtgt   + isizexy
      irtgv   = irtgu   + isizexy
      irtgm   = irtgv   + isizexy
      if13u   = irtgm   + isizexy
      if13v   = if13u   + isizexy
      if13t   = if13v   + isizexy
      if13m   = if13t   + isizexy
      if23u   = if13m   + isizexy
      if23v   = if23u   + isizexy
      if23t   = if23v   + isizexy
      if23m   = if23t   + isizexy
      
      if(ihtran == 1)then
         call polarst(m1,m2,e(iglat),e(iglon),e(ifmapu),e(ifmapv) &
           ,e(ifmapt),e(ifmapm),e(ifmapui),e(ifmapvi) &
           ,e(ifmapti),e(ifmapmi),xm,xt,ym,yt,polelat,polelon,ihtran)
      else
         call ae0(m1*m2,e(ifmapu),1.)
         call ae0(m1*m2,e(ifmapv),1.)
         call ae0(m1*m2,e(ifmapt),1.)
         call ae0(m1*m2,e(ifmapm),1.)
         call ae0(m1*m2,e(ifmapui),1.)
         call ae0(m1*m2,e(ifmapvi),1.)
         call ae0(m1*m2,e(ifmapti),1.)
         call ae0(m1*m2,e(ifmapmi),1.)
      endif

      call grdspc(m1,m2,e(idxu),e(idxv),e(idxt),e(idxm) &
        ,e(idyu),e(idyv),e(idyt),e(idym) &
        ,e(ifmapu),e(ifmapv),e(ifmapt),e(ifmapm) &
        ,xm,xt,ym,yt,jdim,deltay)

      call transfm(m1,m2,m3,topt,e(itopu),e(itopv),e(itopm) &
        ,e(irtgt),e(irtgu),e(irtgv),e(irtgm),e(if13u),e(if13v) &
        ,e(if13t),e(if13m),e(if23u),e(if23v),e(if23t) &
        ,e(if23m),e(idxu),e(idxv),e(idxt),e(idxm),e(idyu),e(idyv)&
        ,e(idyt),e(idym),xm,xt,ym,yt,zm,zt,jdim,ztop,ht,hw)



      if(dir.eq.'xdir')then
         if(gpnt.eq.'upnt')then
            call gradxu(m1,m2,m3,ia,iz,jaa,jzz  &
                 ,optyp,vc3da,vc3db,vctr1,e(irtgu)  &
                 ,e(irtgt),e(idxt),dzt  &
                 ,e(ifmapui),e(ifmapt)  &
                 ,e(if13t)  &
                 ,hw,vctr2,'t',jdim)
         elseif(gpnt.eq.'vpnt')then
            call gradxt(m1,m2,m3,ia,iz,jaa,jzz  &
                 ,optyp,vc3da,vc3db,vctr1,e(irtgv)  &
                 ,e(irtgm),e(idxm),dzt  &
                 ,e(ifmapvi),e(ifmapm)  &
                 ,e(if13m)  &
                 ,hw,vctr2,'t',jdim)
         elseif(gpnt.eq.'wpnt')then
            call gradxt(m1,m2,m3,ia,iz,jaa,jzz  &
                 ,optyp,vc3da,vc3db,vctr1,e(irtgt)  &
                 ,e(irtgu),e(idxu),dzm  &
                 ,e(ifmapti),e(ifmapu)  &
                 ,e(if13u)  &
                 ,ht,vctr2,'w',jdim)
         elseif(gpnt.eq.'tpnt')then
            call gradxt(m1,m2,m3,ia,iz,jaa,jzz  &
                 ,optyp,vc3da,vc3db,vctr1,e(irtgt)  &
                 ,e(irtgu),e(idxu),dzt  &
                 ,e(ifmapti),e(ifmapu)  &
                 ,e(if13u)  &
                 ,hw,vctr2,'t',jdim)
         elseif(gpnt.eq.'npnt')then
            call gradxt(m1,m2,m3,ia,iz,jaa,jzz  &
                 ,optyp,vc3da,vc3db,vctr1,e(irtgv)  &
                 ,e(irtgm),e(idxm),dzm  &
                 ,e(ifmapvi),e(ifmapm)  &
                 ,e(if13m)  &
                 ,ht,vctr2,'w',jdim)
         elseif(gpnt.eq.'opnt')then
            call gradxu(m1,m2,m3,ia,iz,jaa,jzz  &
                 ,optyp,vc3da,vc3db,vctr1,e(irtgu)  &
                 ,e(irtgt),e(idxt),dzm  &
                 ,e(ifmapui),e(ifmapt)  &
                 ,e(if13t)  &
                 ,ht,vctr2,'w',jdim)
         elseif(gpnt.eq.'ppnt')then
            call gradxu(m1,m2,m3,ia,iz,jaa,jzz  &
                 ,optyp,vc3da,vc3db,vctr1,e(irtgm)  &
                 ,e(irtgv),e(idxv),dzt  &
                 ,e(ifmapmi),e(ifmapv)  &
                 ,e(if13v)  &
                 ,hw,vctr2,'t',jdim)
         elseif(gpnt.eq.'mpnt')then
            call gradxu(m1,m2,m3,ia,iz,jaa,jzz  &
                 ,optyp,vc3da,vc3db,vctr1,e(irtgm)  &
                 ,e(irtgv),e(idxv),dzm  &
                 ,e(ifmapmi),e(ifmapv)  &
                 ,e(if13v)  &
                 ,ht,vctr2,'w',jdim)
         endif
      elseif(dir.eq.'ydir')then
         if(gpnt.eq.'upnt')then
            call gradyt(m1,m2,m3,ia,iz,jaa,jzz  &
                 ,optyp,vc3da,vc3db,vctr1,e(irtgu)  &
                 ,e(irtgm),e(idym),dzt  &
                 ,e(ifmapui),e(ifmapm)  &
                 ,e(if23m)  &
               ,hw,vctr2,'t',jdim)
         elseif(gpnt.eq.'vpnt')then
            call gradyv(m1,m2,m3,ia,iz,jaa,jzz  &
                 ,optyp,vc3da,vc3db,vctr1,e(irtgv)  &
                 ,e(irtgt),e(idyt),dzt  &
                 ,e(ifmapvi),e(ifmapt)  &
                 ,e(if23t)  &
                 ,hw,vctr2,'t',jdim)
         elseif(gpnt.eq.'wpnt')then
            call gradyt(m1,m2,m3,ia,iz,jaa,jzz  &
                 ,optyp,vc3da,vc3db,vctr1,e(irtgt)  &
                 ,e(irtgv),e(idyv),dzm  &
                 ,e(ifmapti),e(ifmapv)  &
                 ,e(if23v)  &
                 ,ht,vctr2,'w',jdim)
         elseif(gpnt.eq.'tpnt')then
            call gradyt(m1,m2,m3,ia,iz,jaa,jzz  &
                 ,optyp,vc3da,vc3db,vctr1,e(irtgt)  &
                 ,e(irtgv),e(idyv),dzt  &
                 ,e(ifmapti),e(ifmapv)  &
                 ,e(if23v)  &
                 ,hw,vctr2,'t',jdim)
         elseif(gpnt.eq.'npnt')then
            call gradyv(m1,m2,m3,ia,iz,jaa,jzz  &
                 ,optyp,vc3da,vc3db,vctr1,e(irtgv)  &
                 ,e(irtgt),e(idyt),dzm  &
                 ,e(ifmapvi),e(ifmapt)  &
                 ,e(if23t)  &
                 ,ht,vctr2,'w',jdim)
         elseif(gpnt.eq.'opnt')then
            call gradyt(m1,m2,m3,ia,iz,jaa,jzz  &
                 ,optyp,vc3da,vc3db,vctr1,e(irtgu)  &
                 ,e(irtgm),e(idym),dzm  &
                 ,e(ifmapui),e(ifmapm)  &
                 ,e(if23m)  &
                 ,ht,vctr2,'w',jdim)
         elseif(gpnt.eq.'ppnt')then
            call gradyv(m1,m2,m3,ia,iz,jaa,jzz  &
                 ,optyp,vc3da,vc3db,vctr1,e(irtgm)  &
                 ,e(irtgu),e(idyu),dzt  &
                 ,e(ifmapmi),e(ifmapu)  &
                 ,e(if23u)  &
                 ,hw,vctr2,'t',jdim)
         elseif(gpnt.eq.'mpnt')then
            call gradyv(m1,m2,m3,ia,iz,jaa,jzz  &
                 ,optyp,vc3da,vc3db,vctr1,e(irtgm)  &
                 ,e(irtgu),e(idyu),dzm  &
                 ,e(ifmapmi),e(ifmapu)  &
                 ,e(if23u)  &
                 ,ht,vctr2,'w',jdim)
         endif
      elseif(dir.eq.'zdir')then
         if(gpnt.eq.'upnt')then
            call gradzt(m1,m2,m3,ia,iz,jaa,jzz,vc3da,vc3db  &
                 ,e(irtgu),dzm)
         elseif(gpnt.eq.'vpnt')then
            call gradzt(m1,m2,m3,ia,iz,jaa,jzz,vc3da,vc3db  &
             ,e(irtgv),dzm)
         elseif(gpnt.eq.'wpnt')then
            call gradzw(m1,m2,m3,ia,iz,jaa,jzz,vc3da,vc3db  &
                 ,e(irtgt),dzt)
         elseif(gpnt.eq.'tpnt')then
            call gradzt(m1,m2,m3,ia,iz,jaa,jzz,vc3da,vc3db  &
                 ,e(irtgt),dzm)
         elseif(gpnt.eq.'npnt')then
            call gradzw(m1,m2,m3,ia,iz,jaa,jzz,vc3da,vc3db  &
                 ,e(irtgv),dzt)
         elseif(gpnt.eq.'opnt')then
            call gradzw(m1,m2,m3,ia,iz,jaa,jzz,vc3da,vc3db  &
                 ,e(irtgu),dzt)
         elseif(gpnt.eq.'ppnt')then
            call gradzt(m1,m2,m3,ia,iz,jaa,jzz,vc3da,vc3db  &
                 ,e(irtgm),dzm)
         elseif(gpnt.eq.'mpnt')then
            call gradzw(m1,m2,m3,ia,iz,jaa,jzz,vc3da,vc3db  &
                 ,e(irtgm),dzt)
         endif
      endif

      deallocate(e)

      return
      end

!     ******************************************************************
!
!     this is a general subroutine which computes any component of the
!     gradient or divergence of vc3da and stores it in vc3db.

subroutine gradxu(m1,m2,m3,ia,iz,ja,jz  &
     ,optyp,vc3da,vc3db,vc1da,rtge,rtgc  &
     ,dx,dz,fmapi,fmap,fq,hq,hq4,lev,jd)

  implicit none

  integer, intent(in) :: m1  &
                        ,m2  &
                        ,m3  &
                        ,ia  &
                        ,iz  &
                        ,ja  &
                        ,jz  &
                        ,jd

  real, intent(in) :: vc3da(m1,m2,m3)   &
                    , rtge(m1,m2)       &
                    , rtgc(m1,m2)       &
                    , dx(m1,m2)         &
                    , fmap(m1,m2)       &
                    , fmapi(m1,m2)      &
                    , dz(m3)            &
                    , fq(m1,m2)         &
                    , hq(*)

  real, intent(inout)  :: vc3db(m1,m2,m3)   &
                        , vc1da(*)          &
                        , hq4(*)

  character(len=*), intent(in) :: optyp,lev
  
  integer :: i,j,k
  
  if(optyp.eq.'gradnt')then
     do j=ja,jz
        do i=ia,iz
           do k=1,m3
              vc3db(i,j,k)=(vc3da(i,j,k)*rtge(i,j)  &
                   -vc3da(i-1,j,k)*rtge(i-1,j))  &
                   *dx(i,j)/rtgc(i,j)
           enddo
        enddo
     enddo
  else
     do j=ja,jz
        do i=ia,iz
           do k=1,m3
              vc3db(i,j,k)=(vc3da(i,j,k)*rtge(i,j)  &
                   *fmapi(i,j)  &
                   -vc3da(i-1,j,k)*rtge(i-1,j)  &
                   *fmapi(i-1,j))  &
                   *dx(i,j)/rtgc(i,j)*fmap(i,j)
           enddo
        enddo
     enddo
  endif

  if(optyp.ne.'divstr')then
     if(lev.eq.'w')then
        do k=1,m3
           hq4(k)=0.25*hq(k)
        enddo
     else
        do k=2,m3
           hq4(k)=0.25*hq(k-1)
        enddo
     endif

   ! **(jp)** quebra aninhamento, permitindo vetorizacao
   ! **(jp)** dos ultimos aninhamentos em i,j. o primeiro
   ! **(jp)** continua vetorizado em k. 
!!$   do j=ja,jz
!!$      do i=ia,iz
!!$         do k=2,m1
!!$            vc1da(k)=hq4(k)*(vc3da(k,i,j)+vc3da(k-1,i,j)  &
!!$                 +vc3da(k,i-1,j)+vc3da(k-1,i-1,j))
!!$         enddo
!!$         do k=2,m1-1
!!$            vc3db(k,i,j)=vc3db(k,i,j)  &
!!$                 +fq(i,j)*dz(k)*(vc1da(k+1)-vc1da(k))
!!$         enddo
!!$         vc3db(1,i,j)=vc3db(2,i,j)
!!$         if(lev.eq.'w')vc3db(m1-1,i,j)=vc3db(m1-2,i,j)
!!$         if(lev.eq.'t')vc3db(m1,i,j)=vc3db(m1-1,i,j)
!!$      enddo
!!$   enddo
     do j=ja,jz
        do i=ia,iz
           do k=2,m3
              vc1da(k)=hq4(k)*(vc3da(i,j,k)+vc3da(i,j,k-1)  &
                   +vc3da(i-1,j,k)+vc3da(i-1,j,k-1))
           enddo
           do k=2,m3-1
              vc3db(i,j,k)=vc3db(i,j,k)  &
                   +fq(i,j)*dz(k)*(vc1da(k+1)-vc1da(k))
           enddo
        enddo
     enddo
     do j=ja,jz
        do i=ia,iz
           vc3db(i,j,1)=vc3db(i,j,2)
        enddo
     enddo
     if(lev.eq.'w')then
        do j=ja,jz
           do i=ia,iz
              vc3db(i,j,m3-1)=vc3db(i,j,m3-2)
           end do
        end do
     else if(lev.eq.'t')then
        do j=ja,jz
           do i=ia,iz
              vc3db(i,j,m3)=vc3db(i,j,m3-1)
           end do
        end do
     end if
     ! **(jp)** fim modificacao
  endif
  
  return
end subroutine gradxu

subroutine gradxt(m1,m2,m3,ia,iz,ja,jz  &
     ,optyp,vc3da,vc3db,vc1da,rtge,rtgc  &
     ,dx,dz,fmapi,fmap,fq,hq,hq4,lev,jd)

  implicit none

  integer, intent(in) :: m1   &
                       , m2   &
                       , m3   &
                       , ia   &
                       , iz   & 
                       , ja   &
                       , jz   &
                       , jd

  real, intent(in)    :: vc3da(m1,m2,m3)  &
                       , rtge(m1,m2)      &
                       , rtgc(m1,m2)      &
                       , dx(m1,m2)        &
                       , fmap(m1,m2)      &
                       , fmapi(m1,m2)     &
                       , dz(*)            &
                       , fq(m1,m2)        &
                       , hq(*)

  real, intent(inout) :: vc1da(*)         &
                       , vc3db(m1,m2,m3)  &
                       , hq4(*)           

  character(len=*), intent(in) :: optyp,lev
  

  integer :: i,j,k
  
  if(optyp.eq.'gradnt')then
     do j=ja,jz
        do i=ia,iz
           do k=1,m3
              vc3db(i,j,k)=(vc3da(i+1,j,k)*rtge(i+1,j)  &
                   -vc3da(i,j,k)*rtge(i,j))  &
                   *dx(i,j)/rtgc(i,j)
           enddo
        enddo
     enddo
  else
     do j=ja,jz
        do i=ia,iz
           do k=1,m3
              vc3db(i,j,k)=(vc3da(i+1,j,k)*rtge(i+1,j)  &
                   *fmapi(i+1,j)  &
                   -vc3da(i,j,k)*rtge(i,j)  &
                   *fmapi(i,j))  &
                   *dx(i,j)/rtgc(i,j)*fmap(i,j)
           enddo
        enddo
     enddo
  endif

  if(optyp.ne.'divstr')then
     if(lev.eq.'w')then
        do k=1,m3
           hq4(k)=0.25*hq(k)
        enddo
     else
        do k=2,m3
           hq4(k)=0.25*hq(k-1)
        enddo
     endif

   ! **(jp)** quebra aninhamento, permitindo vetorizacao
   ! **(jp)** dos ultimos aninhamentos em i,j. o primeiro
   ! **(jp)** continua vetorizado em k. 
!!$   do j=ja,jz
!!$      do i=ia,iz
!!$         do k=2,m1
!!$            vc1da(k)=hq4(k)*(vc3da(k,i,j)+vc3da(k-1,i,j)  &
!!$                 +vc3da(k,i+1,j)+vc3da(k-1,i+1,j))
!!$         enddo
!!$         do k=2,m1-1
!!$            vc3db(k,i,j)=vc3db(k,i,j)  &
!!$                 +fq(i,j)*dz(k)*(vc1da(k+1)-vc1da(k))
!!$         enddo
!!$         vc3db(1,i,j)=vc3db(2,i,j)
!!$         if(lev.eq.'w')vc3db(m1-1,i,j)=vc3db(m1-2,i,j)
!!$         if(lev.eq.'t')vc3db(m1,i,j)=vc3db(m1-1,i,j)
!!$      enddo
!!$   enddo
     do j=ja,jz
        do i=ia,iz
           do k=2,m3
              vc1da(k)=hq4(k)*(vc3da(i,j,k)+vc3da(i,j,k-1)  &
                   +vc3da(i+1,j,k)+vc3da(i+1,j,k-1))
           enddo
           do k=2,m3-1
              vc3db(i,j,k)=vc3db(i,j,k)  &
                   +fq(i,j)*dz(k)*(vc1da(k+1)-vc1da(k))
           enddo
        enddo
     enddo
     do j=ja,jz
        do i=ia,iz
           vc3db(i,j,1)=vc3db(i,j,2)
        enddo
     enddo
     if (lev.eq.'W') then
        do j=ja,jz
           do i=ia,iz
              vc3db(i,j,m3-1)=vc3db(i,j,m3-2)
           enddo
        enddo
     else if(lev.eq.'T') then
        do j=ja,jz
           do i=ia,iz
              vc3db(i,j,m3)=vc3db(i,j,m3-1)
           enddo
        enddo
     end if
     ! **(jp)** fim modificacao
  endif
  
  return
end subroutine gradxt

!

subroutine gradyv(m1,m2,m3,ia,iz,ja,jz  &
     ,optyp,vc3da,vc3db,vc1da,rtge,rtgc  &
     ,dy,dz,fmapi,fmap,fq,hq,hq4,lev,jd)

  implicit none

  integer, intent(in) :: m1   &
                       , m2   &
                       , m3   &
                       , ia   &
                       , iz   &
                       , ja   &
                       , jz   &
                       , jd

  real, intent(in)  :: vc3da(m1,m2,m3)   &
                     , rtge(m1,m2)       &
                     , rtgc(m1,m2)       &
                     , dy(m1,m2)         &
                     , fmap(m1,m2)       &
                     , fmapi(m1,m2)      &
                     , dz(*)             &
                     , fq(m1,m2)         &
                     , hq(*)

  real, intent(inout) :: vc3db(m1,m2,m3)   &
                       , hq4(*)            &
                       , vc1da(*)


  character(len=*) :: optyp,lev

  integer :: i,j,k

  if(optyp.eq.'gradnt')then
     do j=ja,jz
        do i=ia,iz
           do k=1,m3
              vc3db(i,j,k)=(vc3da(i,j,k)*rtge(i,j)  &
                   -vc3da(i,j-jd,k)*rtge(i,j-jd))  &
                   *dy(i,j)/rtgc(i,j)
           enddo
        enddo
     enddo
  else
     do j=ja,jz
        do i=ia,iz
           do k=1,m3
              vc3db(i,j,k)=(vc3da(i,j,k)*rtge(i,j)  &
                   *fmapi(i,j)  &
                   -vc3da(i,j-jd,k)*rtge(i,j-jd)  &
                   *fmapi(i,j-jd))  &
                   *dy(i,j)/rtgc(i,j)*fmap(i,j)
           enddo
        enddo
     enddo
  endif

  if(optyp.ne.'divstr')then
     if(lev.eq.'w')then
        do k=1,m3
           hq4(k)=0.25*hq(k)
        enddo
     else
        do k=2,m3
           hq4(k)=0.25*hq(k-1)
        enddo
     endif

   ! **(jp)** quebra aninhamento, permitindo vetorizacao
   ! **(jp)** dos ultimos aninhamentos em i,j. o primeiro
   ! **(jp)** continua vetorizado em k. 
!!$   do j=ja,jz
!!$      do i=ia,iz
!!$         do k=2,m1
!!$            vc1da(k)=hq4(k)*(vc3da(k,i,j)+vc3da(k-1,i,j)  &
!!$                 +vc3da(k,i,j-jd)+vc3da(k-1,i,j-jd))
!!$         enddo
!!$         do k=2,m1-1
!!$            vc3db(k,i,j)=vc3db(k,i,j)  &
!!$                 +fq(i,j)*dz(k)*(vc1da(k+1)-vc1da(k))
!!$         enddo
!!$         vc3db(1,i,j)=vc3db(2,i,j)
!!$         if(lev.eq.'w')vc3db(m1-1,i,j)=vc3db(m1-2,i,j)
!!$         if(lev.eq.'t')vc3db(m1,i,j)=vc3db(m1-1,i,j)
!!$      enddo
!!$   enddo
     do j=ja,jz
        do i=ia,iz
           do k=2,m3
              vc1da(k)=hq4(k)*(vc3da(i,j,k)+vc3da(i,j,k-1)  &
                   +vc3da(i,j-jd,k)+vc3da(i,j-jd,k-1))
           enddo
           do k=2,m3-1
              vc3db(i,j,k)=vc3db(i,j,k)  &
                   +fq(i,j)*dz(k)*(vc1da(k+1)-vc1da(k))
           enddo
        enddo
     enddo
     do j=ja,jz
        do i=ia,iz
           vc3db(i,j,1)=vc3db(i,j,2)
        enddo
     enddo
     if (lev.eq.'w') then
        do j=ja,jz
           do i=ia,iz
              vc3db(i,j,m3-1)=vc3db(i,j,m3-2)
           enddo
        enddo
     else if (lev.eq.'t') then
        do j=ja,jz
           do i=ia,iz
              vc3db(i,j,m3)=vc3db(i,j,m3-1)
           enddo
        enddo
     end if
     ! **(jp)** fim modificacao
  endif

  return
end subroutine gradyv

subroutine gradyt(m1,m2,m3,ia,iz,ja,jz  &
     ,optyp,vc3da,vc3db,vc1da,rtge,rtgc  &
     ,dy,dz,fmapi,fmap,fq,hq,hq4,lev,jd)

 implicit none

  integer, intent(in) :: m1    &
                       , m2    &
                       , m3    &
                       , ia    &
                       , iz    &
                       , ja    &
                       , jz    &
                       , jd

  real, intent(in) :: vc3da(m1,m2,m3)   &
                    , rtge(m1,m2)       &
                    , rtgc(m1,m2)       &
                    , dy(m1,m2)         &
                    , fmap(m1,m2)       &
                    , fmapi(m1,m2)      &
                    , dz(*)             &
                    , fq(m1,m2)         &
                    , hq(*)

  real, intent(inout) :: vc3db(m1,m2,m3)   &
                       , hq4(*)            &
                       , vc1da(*)          

  character(len=*), intent(in) :: optyp, lev
  
  integer :: i,j,k

  if(optyp.eq.'gradnt')then
     do j=ja,jz
        do i=ia,iz
           do k=1,m3
              !write(*,'(a)') '------------------------------------------'
              !write(*,'(3(a,1x,i5,1x))') ' i=',i,'j=',j,'k=',k
              !write(*,'(a,1x,es12.5)') ' + vc3db_ijk=',vc3db(i,j,k)
              !write(*,'(a,1x,es12.5)') ' + vc3da_ij+k=',vc3da(i,j+jd,k)
              !write(*,'(a,1x,es12.5)') ' + rtge_ij+k=',rtge(i,j+jd)
              !write(*,'(a,1x,es12.5)') ' + vc3da_ijk=',vc3da(i,j,k)
              !write(*,'(a,1x,es12.5)') ' + rtge_ijk=',rtge(i,j)
              !write(*,'(a,1x,es12.5)') ' + dy=',dy(i,j)
              !write(*,'(a,1x,es12.5)') ' + rtgc=',rtgc(i,j)
              !write(*,'(a)') '------------------------------------------'
              !write(*,'(a)') ' '
              
              vc3db(i,j,k)=(vc3da(i,j+jd,k)*rtge(i,j+jd)  &
                   -vc3da(i,j,k)*rtge(i,j))*dy(i,j)/rtgc(i,j)
           enddo
        enddo
     enddo
  else
     do j=ja,jz
        do i=ia,iz
           do k=1,m3
              vc3db(i,j,k)=(vc3da(i,j+jd,k)*rtge(i,j+jd)  &
                   *fmapi(i,j+jd)  &
                   -vc3da(i,j,k)*rtge(i,j)  &
                   *fmapi(i,j))  &
                   *dy(i,j)/rtgc(i,j)*fmap(i,j)
           enddo
        enddo
     enddo
  endif

  if(optyp.ne.'divstr')then
     if(lev.eq.'w')then
        do k=1,m3
           hq4(k)=0.25*hq(k)
        enddo
     else
        do k=2,m3
           hq4(k)=0.25*hq(k-1)
        enddo
     endif

   ! **(jp)** quebra aninhamento, permitindo vetorizacao
   ! **(jp)** dos ultimos aninhamentos em i,j. o primeiro
   ! **(jp)** continua vetorizado em k. 
!!$   do j=ja,jz
!!$      do i=ia,iz
!!$         do k=2,m1
!!$            vc1da(k)=hq4(k)*(vc3da(k,i,j)+vc3da(k-1,i,j)  &
!!$                 +vc3da(k,i,j+jd)+vc3da(k-1,i,j+jd))
!!$         enddo
!!$         do k=2,m1-1
!!$            vc3db(k,i,j)=vc3db(k,i,j)  &
!!$                 +fq(i,j)*dz(k)*(vc1da(k+1)-vc1da(k))
!!$         enddo
!!$         vc3db(1,i,j)=vc3db(2,i,j)
!!$         if(lev.eq.'w')vc3db(m1-1,i,j)=vc3db(m1-2,i,j)
!!$         if(lev.eq.'t')vc3db(m1,i,j)=vc3db(m1-1,i,j)
!!$      enddo
!!$   enddo
     do j=ja,jz
        do i=ia,iz
           do k=2,m3
              vc1da(k)=hq4(k)*(vc3da(i,j,k)+vc3da(i,j,k-1)  &
                   +vc3da(i,j+jd,k)+vc3da(i,j+jd,k-1))
           enddo
           do k=2,m3-1
              vc3db(i,j,k)=vc3db(i,j,k)  &
                   +fq(i,j)*dz(k)*(vc1da(k+1)-vc1da(k))
           enddo
        enddo
     enddo
     do j=ja,jz
        do i=ia,iz
           vc3db(i,j,1)=vc3db(i,j,2)
        enddo
     enddo
     if (lev.eq.'w') then
        do j=ja,jz
           do i=ia,iz
              vc3db(i,j,m3-1)=vc3db(i,j,m3-2)
           enddo
        enddo
     else if (lev.eq.'t') then
        do j=ja,jz
           do i=ia,iz
              vc3db(i,j,m3)=vc3db(i,j,m3-1)
           enddo
        enddo
     end if
     ! **(jp)** fim modificacao
  endif
  
  return
end subroutine gradyt

subroutine gradzw(m1,m2,m3,ia,iz,ja,jz,vc3da,vc3db,rtgc,dz)

  implicit none

  integer, intent(in) :: m1  &
                       , m2  &
                       , m3  &
                       , ia  &
                       , iz  &
                       , ja  &
                       , jz

  real, intent(in)    :: vc3da(m1,m2,m3)  &
                       , rtgc(m1,m2)      &
                       , dz(*)

  real, intent(inout) :: vc3db(m1,m2,m3) 
  
  integer :: i,j,k

  do j=ja,jz
     do i=ia,iz
        do k=2,m3
           vc3db(i,j,k)=(vc3da(i,j,k)-vc3da(i,j,k-1))*dz(k)  &
                /rtgc(i,j)
        enddo
     enddo
  enddo
  return
end subroutine gradzw

subroutine gradzt(m1,m2,m3,ia,iz,ja,jz,vc3da,vc3db,rtgc,dz)

  implicit none

  integer :: m1,m2,m3,ia,iz,ja,jz

  real, intent(in)    :: vc3da(m1,m2,m3) &
                       , rtgc(m1,m2)     &
                       , dz(*)

  real, intent(inout) :: vc3db(m1,m2,m3)

  integer :: i,j,k

  do j=ja,jz
     do i=ia,iz
        do k=1,m3-1
           vc3db(i,j,k)=(vc3da(i,j,k+1)-vc3da(i,j,k))*dz(k)  &
                /rtgc(i,j)
        enddo
     enddo
  enddo
  return
end subroutine gradzt



!     ******************************************************************
!
      subroutine xytops(x,y,pla,plo,erad)
!
!     this convert x,y-polar stereographic coordinates to
!     lat/lon values.
!     longitude:   0 - 360  ; positive to the east
!     latitude : -90 -  90  ; positive for northern hemisphere
!     it is assumed that the x-axis point towards the east why the
!     longitude is rotated relative to the 'standard pol.ste.' location
!            

      pi180=3.14159/180.
      vdist = erad*2.0
!
!     calculate distance from (0,0) and define (0,0) as 90,0 (90,-90 in
!     the rotated system)
!
      dist=(x**2+y**2)**0.5

      if(dist.eq.0) then
         pla= 90.0
         plo=-90.0
      else
!
!     calculate the latitude by means of atan
!
         pla=atan(dist/vdist)/pi180
         pla=90.0-2.0*pla
!
!     calculate the longitude taking the directions into account
!
         if(x.eq.0.0) then
            if(y.gt.0.0) then
               plo= 90.0
            else
               plo=-90.0
            end if
         else
            if(x.gt.0.0) then
               plo=atan(y/x)/pi180
            else
               plo=atan(y/x)/pi180+180.0
            end if
         end if
      end if
!
!     rotate the longitude
!
      plo=amod(plo+450.0,360.0)
      return
      end
!-------------------------------------------------------------------
!      subroutine ll_xy(xlat,xlon,plat,plon,x,y)
!      call getops(pla,plo,xlat,xlon,plat,plon)
!      call pstoxy(x,y,pla,plo,6376000.)
!      return
!      end
!
!
!      subroutine xy_ll(xlat,xlon,plat,plon,x,y)
!      call xytops(x,y,pla,plo,6376000.)
!      call pstoge(pla,plo,xlat,xlon,plat,plon)
!      return
!      end
!
!      function hg_xy(nih,xh,hg)
!      dimension xh(*)
!      ix=int(hg)
!      hg_xy=xh(ix)+(xh(ix+1)-xh(ix))*(hg-ix)
!      return
!      end
!
!      function hifromx(nih,xh,x)
!      dimension xh(*)
!
!!      print*,'hifromx-',x,nih,xh(1),xh(nih)
!      if(x.lt.xh(1).or.x.gt.xh(nih))then
!         print*, 'x, y, or z value exceeds hgrid limits'
!         stop 'xtohi'
!      endif
!
!      ilow=1
!      ihigh=nih
!      do m=2,11
!         imid=(ilow+ihigh)/2
!         msw=int(max(0.,min(1.5,1.e20*(x-xh(imid)))))
!         ilow=ilow*(1-msw)+imid*msw
!         ihigh=ihigh*msw+imid*(1-msw)
!      enddo
!      hifromx=float(ilow)+(x-xh(ilow))/(xh(ihigh)-xh(ilow))
!      end


!
!      function hg_z(nih,njh,nkh,xh,yh,zh,hgx,hgy,hgz,topth,nh1,nh2,ztop
!     +     ,topo)
!      dimension topth(nh1,nh2),xh(*),yh(*),zh(*)
!
!      ihp=int(hgx)
!      jhp=int(hgy)
!      qh20=hgx-float(ihp)
!      qh02=hgy-float(jhp)
!      topo=(1.0-qh02)*((1.0-qh20)*topth(ihp  ,jhp  )
!     +     +qh20 *topth(ihp+1,jhp  ))
!     +     +qh02 *((1.0-qh20)*topth(ihp  ,jhp+1)
!     +     +qh20 *topth(ihp+1,jhp+1))
!      rtg=1.-topo/ztop
!
!
!      iz=int(hgz)
!      hg_z=topo+(zh(iz)+(zh(iz+1)-zh(iz))*(hgz-iz))*rtg
!      print*,'hg_z-',hgx,hgy,hgz,topo,hg_z,topth(ihp,jhp)
!      return
!      end
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine uvtouv(up,vp,ue,ve,x,y,dir,rlat,wlon1,erad)
!
!     converts u,v in polar stereographic space to u,v in normal space
!     for dir='stog' else the convertion will be from normal to polar
!     stereographic space.
!     the conversion is made by transforming to dd,ff and rotate dd
!     after which the new dd,ff is converted back to u,v
!     note: dd is not the normal dd but the angle between the x-axis
!           and the wind vector
!
!     tsp 03 july 89
!
      character*4 dir
      pi180=3.14159/180.0
!
!     determine the location in pla,plo- and gla,glo coordinates
!
      call xytops(x,y,pla,plo,erad)
      call pstoge(pla,plo,gla,glo,rlat,wlon1)
!
!     determine the rotaion angle alpha which depends on the sign of rlat
!
      hsign=1.0
      if(rlat.lt.0.0) hsign=-1.0
      arga1=cos(pla*pi180)*cos(hsign*gla*pi180)
      if(arga1.ne.0.0.and.((x**2+y**2).gt.1.0)) then
         arga2=(sin(hsign*rlat*pi180)- &
               sin(pla*pi180)*sin(hsign*gla*pi180))/arga1
      else
         arga2=1.0
      end if
      if(arga2.gt.1.0) then
         arga2=1.0
      else if(arga2.lt.-1.0) then
         arga2=-1.0
      end if
      alpha=acos(arga2)
!
!     determine rotation angle and sin/cos values between pol.ste. and
!     the geo. system
!
      rotang=(plo+270.0)*pi180
      rsin = sin(rotang)
      rcos = cos(rotang)
!
      if(dir.eq.'stog') then
!
!     convert up,vp to us,vs (ie from x,y directions in pol.ste. to
!     east, north in the pla, plo pol.ste. rotated system
!
        us   =-rsin*up + rcos*vp
        vs   =-rcos*up - rsin*vp
!
!     convert us,vs to dd,ff with dd <-180, 180>
!         
         ff=sqrt(us**2+vs**2)
         if(us.eq.0.0.and.vs.ge.0.0) then
            dd=pi180*90.0
         else if(us.eq.0.0.and.vs.lt.0.0) then
            dd=pi180*270.0
         else if(us.gt.0.0) then
            dd=atan(vs/us)
         else
            dd=atan(vs/us)
            dd=dd+pi180*180.0
         end if
!
!     the actual rotation depends on rlat and plo
!
         if(rlat.ge.0.0) then
            if(plo.lt.180.0) then
               dd=dd+alpha
            else
               dd=dd-alpha
            end if
         else
            if(plo.lt.180.0) then
               dd=dd-alpha+pi180*180.0
            else
               dd=dd+alpha+pi180*180.0
            end if
         end if
!
!     conversion from dd,ff to ue,ve
!
         ue=ff*cos(dd)
         ve=ff*sin(dd)
!
!     this corresponds to the conversion from normal to polar stereographic
!         
      else
!         
!     convert ue,ve to dd,ff with dd <-180, 180>
!         
         ff=sqrt(ue**2+ve**2)
         if(ue.eq.0.0.and.ve.ge.0.0) then
            dd=pi180*90.0
         else if(ue.eq.0.0.and.ve.lt.0.0) then
            dd=pi180*270.0
         else if(ue.gt.0.0) then
            dd=atan(ve/ue)
         else
            dd=atan(ve/ue)
            dd=dd+pi180*180.0
         end if
!
!     the actual rotation depends on rlat and plo
!
         if(rlat.ge.0.0) then
            if(plo.lt.180.0) then
               dd=dd-alpha
            else
               dd=dd+alpha
            end if
         else
            if(plo.lt.180.0) then
               dd=dd+alpha-pi180*180.0
            else
               dd=dd-alpha-pi180*180.0
            end if
         end if
!
!     conversion from dd,ff to us,vs
!
         us=ff*cos(dd)
         vs=ff*sin(dd)
!
!     convert us,vs to up,vp (ie from east, north directions in pol.ste. to
!     x,y in the pol.ste. rotated system
!
        up   =-rsin*us - rcos*vs
        vp   = rcos*us - rsin*vs
      end if
!
      return
      end
!
!<~>~<~<~>~<~<~>~<~<~>~<~<~>~<~<~>~<~<~>~<~<~>~<~<~>~<~<~>~<~<~>~<~<~>~<~<~>~<~<~>~<~>!
!
subroutine transfm(n2,n3,nzp,topt,topu,topv,topm,rtgt,rtgu,rtgv,rtgm  &
     ,f13u,f13v,f13t,f13m,f23u,f23v,f23t,f23m  &
     ,dxu,dxv,dxt,dxm,dyu,dyv,dyt,dym,xm,xt,ym,yt,zm,zt,jdim,ztop,ht,hw)


  implicit none

  integer :: n2,n3,nzp,jdim
  real, dimension(n2,n3) :: topt,topu,topv,topm,rtgt,rtgu,rtgv,rtgm  &
       ,f13u,f13v,f13t,f13m,f23u,f23v,f23t,f23m  &
       ,dxu,dxv,dxt,dxm,dyu,dyv,dyt,dym
  real, dimension(*) :: xm,xt,ym,yt,zm,zt,ht,hw
  real :: ztop

  integer :: iztflag=0,i,j,k

  !     this routine computes the coordinate transformation constants
  !     based on the topographical values of topt.

  do k = 1,nzp
     ht(k)  = zt(k) / ztop - 1.
     hw(k)  = zm(k) / ztop - 1.
  enddo

  do j = 1,n3
     do i = 1,n2
        rtgt(i,j) = 1. - topt(i,j) / ztop
        rtgu(i,j) = 1. - topu(i,j) / ztop
        rtgv(i,j) = 1. - topv(i,j) / ztop
        rtgm(i,j) = 1. - topm(i,j) / ztop
     enddo
  enddo
  do j = 1,n3
     do i = 2,n2
        f13t(i,j) = (topu(i,j) - topu(i-1,j)) * dxt(i,j) / rtgt(i,j)
        f13v(i,j) = (topm(i,j) - topm(i-1,j)) * dxv(i,j) / rtgv(i,j)
     enddo
     do i = 1,n2-1
        f13u(i,j) = (topt(i+1,j) - topt(i,j)) * dxu(i,j) / rtgu(i,j)
        f13m(i,j) = (topv(i+1,j) - topv(i,j)) * dxm(i,j) / rtgm(i,j)
     enddo
     f13t(1,j)  = f13u(1,j)
     f13v(1,j)  = f13m(1,j)
     f13u(n2,j) = f13t(n2,j)
     f13m(n2,j) = f13v(n2,j)
  enddo

  do i = 1,n2
     do j = 2,n3
        f23t(i,j) = (topv(i,j) - topv(i,j-jdim)) * dyt(i,j) / rtgt(i,j)
        f23u(i,j) = (topm(i,j) - topm(i,j-jdim)) * dyu(i,j) / rtgu(i,j)
     enddo
     do j = 1,n3-1
        f23v(i,j) = (topt(i,j+jdim) - topt(i,j)) * dyv(i,j) / rtgv(i,j)
        f23m(i,j) = (topu(i,j+jdim) - topu(i,j)) * dym(i,j) / rtgm(i,j)
     enddo
     if (jdim == 1) then
        f23t(i,1)  = f23v(i,1)
        f23u(i,1)  = f23m(i,1)
        f23v(i,n3) = f23t(i,n3)
        f23m(i,n3) = f23u(i,n3)
     endif
  enddo
  return
end subroutine transfm
!
!<~>~<~<~>~<~<~>~<~<~>~<~<~>~<~<~>~<~<~>~<~<~>~<~<~>~<~<~>~<~<~>~<~<~>~<~<~>~<~<~>~<~>!
!
subroutine polarst(n2,n3,glat,glon,fmapu,fmapv,fmapt,fmapm  &
     ,fmapui,fmapvi,fmapti,fmapmi,xm,xt,ym,yt,platn,plonn,ihtran)

  implicit none
  real, parameter :: erad=6367000.
  integer :: n2,n3
  real, dimension(n2,n3) :: glat,glon,fmapu,fmapv,fmapt,fmapm  &
       ,fmapui,fmapvi,fmapti,fmapmi
  real, dimension(*) :: xm,xt
  real, dimension(*) :: ym,yt
  integer :: ihtran
  real :: platn,plonn

  integer :: i,j
  real :: c1,xm2,xt2,ym2,yt2
  !  Calculates map factors and inverse map factors at u,v,t,m-points and
  !  geographical lat/lon at t-points for a given polar stereographic grid

  c1 = (2. * erad) ** 2
  do j = 1,n3
     do i = 1,n2
        xm2 = xm(i) * xm(i)
        xt2 = xt(i) * xt(i)
        ym2 = ym(j) * ym(j)
        yt2 = yt(j) * yt(j)

        fmapt(i,j) = 1. + (xt2 + yt2) / c1
        fmapu(i,j) = 1. + (xm2 + yt2) / c1
        fmapv(i,j) = 1. + (xt2 + ym2) / c1
        fmapm(i,j) = 1. + (xm2 + ym2) / c1

        fmapui(i,j) = 1.0 / fmapu(i,j)
        fmapvi(i,j) = 1.0 / fmapv(i,j)
        fmapti(i,j) = 1.0 / fmapt(i,j)
        fmapmi(i,j) = 1.0 / fmapm(i,j)

        call xy_ll(glat(i,j),glon(i,j),platn,plonn,xt(i),yt(j))

        !       write(6,344)i,j,fmapt(i,j),fmapm(i,j),glat(i,j),glon(i,j)
        ! 344   format('polst:i,j,fmt,fmm,glt,gln',2i4,4e12.3)

     enddo
  enddo

  if (ihtran == 0) then
     call ae0(n2*n3,fmapu,1.)
     call ae0(n2*n3,fmapv,1.)
     call ae0(n2*n3,fmapt,1.)
     call ae0(n2*n3,fmapm,1.)
     call ae0(n2*n3,fmapui,1.)
     call ae0(n2*n3,fmapvi,1.)
     call ae0(n2*n3,fmapti,1.)
     call ae0(n2*n3,fmapmi,1.)
  endif
end subroutine polarst
!
!<~>~<~<~>~<~<~>~<~<~>~<~<~>~<~<~>~<~<~>~<~<~>~<~<~>~<~<~>~<~<~>~<~<~>~<~<~>~<~<~>~<~>!
!
!subroutine xy_ll (qlat,qlon,polelat,polelon,x,y)
! implicit none
! real, intent(out) :: qlat
! real, intent(out) :: qlon
! real, intent(in ) :: polelat
! real, intent(in ) :: polelon
! real, intent(in ) :: x
! real, intent(in ) :: y

! real, parameter :: erad=6.367e6
! real, parameter :: erad2=1.2734e7
! real, parameter :: pi180=3.14159265/180.

! real :: sinplat
! real :: cosplat
! real :: sinplon
! real :: cosplon
! real :: x3p
! real :: y3p
! real :: z3p
! real :: z3q
! real :: x3q
! real :: y3q
! real :: xq
! real :: yq
! real :: zq
! real :: t
! real :: d
! real :: alpha
! real :: r3q


! ! Evaluate sine and cosine of latitude and longitude of pole point p.

! sinplat = sin(polelat * pi180)
! cosplat = cos(polelat * pi180)
! sinplon = sin(polelon * pi180)
! cosplon = cos(polelon * pi180)

! ! Compute (x3,y3,z3) coordinates of the pole point where the origin is the
! ! center of the earth, the z axis is the north pole, the x axis is the
! ! equator and prime meridian, and the y axis is the equator and 90 E.

! x3p = erad * cosplat * cosplon
! y3p = erad * cosplat * sinplon
! z3p = erad * sinplat

! ! Compute distance d from given point R on the polar stereographic plane
! ! to the pole point P:

! d = sqrt (x ** 2 + y ** 2)

! ! Compute angle QCP where C is the center of the Earth.  This is twice
! ! angle QAP where A is the antipodal point.  Angle QAP is the same as
! ! angle RAP:

! alpha = 2. * atan2(d,erad2)

! ! Compute zq, the height of Q relative to the polar stereographic plane:

! zq = erad * (cos(alpha) - 1.)

! ! Compute the parameter t which is the the distance ratio AQ:AR

! t = (erad2 + zq) / erad2

! ! Compute xq and yq, the x and y coordinates of Q in polar stereographic space:

! xq = t * x
! yq = t * y

! ! Transform location of Q from (x,y,z) coordinates to (x3,y3,z3):

! x3q = x3p - xq * sinplon - yq * cosplon * sinplat  &
!      + zq * cosplat * cosplon
! y3q = y3p + xq * cosplon - yq * sinplon * sinplat  &
!      + zq * cosplat * sinplon
! z3q = z3p + yq * cosplat + zq * sinplat

! ! Compute the latitude and longitude of Q:

! qlon = atan2(y3q,x3q) / pi180
! r3q = sqrt(x3q ** 2 + y3q ** 2)
! qlat = atan2(z3q,r3q) / pi180

!end subroutine xy_ll
!
!<~>~<~<~>~<~<~>~<~<~>~<~<~>~<~<~>~<~<~>~<~<~>~<~<~>~<~<~>~<~<~>~<~<~>~<~<~>~<~<~>~<~>!
!
subroutine grdspc(n2,n3,dxu,dxv,dxt,dxm,dyu,dyv,dyt,dym  &
        ,fmapu,fmapv,fmapt,fmapm,xm,xt,ym,yt,jdim,deltay)
  implicit none

  integer :: n2,n3
  real, dimension(n2,n3) :: dxu,dxv,dxt,dxm,dyu,dyv,dyt,dym  &
       ,fmapu,fmapv,fmapt,fmapm
  integer :: jdim
  real :: deltay
  real, dimension(*) :: xm,xt,yt,ym

  integer :: i,j

  do j = 1,n3
     do i = 1,n2-1
        dxu(i,j) = fmapu(i,j) / (xt(i+1)-xt(i))
        dxm(i,j) = fmapm(i,j) / (xt(i+1)-xt(i))
     enddo
     dxu(n2,j)=dxu(n2-1,j)*fmapu(n2,j)/fmapu(n2-1,j)
     dxm(n2,j)=dxm(n2-1,j)*fmapm(n2,j)/fmapm(n2-1,j)
     do i = 2,n2
        dxv(i,j)=fmapv(i,j)/(xm(i)-xm(i-1))
        dxt(i,j)=fmapt(i,j)/(xm(i)-xm(i-1))
     enddo
     dxv(1,j)=dxv(2,j)*fmapv(1,j)/fmapv(2,j)
     dxt(1,j)=dxt(2,j)*fmapt(1,j)/fmapt(2,j)
  enddo

  if (jdim == 1) then
     do i = 1,n2
        do j = 1,n3-1
           dyv(i,j)=fmapv(i,j)/(yt(j+1)-yt(j))
           dym(i,j)=fmapm(i,j)/(yt(j+1)-yt(j))
        enddo
        dyv(i,n3)=dyv(i,n3-1)*fmapv(i,n3)/fmapv(i,n3-1)
        dym(i,n3)=dym(i,n3-1)*fmapm(i,n3)/fmapm(i,n3-1)
        do j = 2,n3
           dyu(i,j)=fmapu(i,j)/(ym(j)-ym(j-1))
           dyt(i,j)=fmapt(i,j)/(ym(j)-ym(j-1))
        enddo
        dyu(i,1)=dyu(i,2)*fmapu(i,1)/fmapu(i,2)
        dyt(i,1)=dyt(i,2)*fmapt(i,1)/fmapt(i,2)
     enddo
  else
     do i=1,n2
        do j=1,n3
           dyu(i,j)=1./deltay
           dyv(i,j)=1./deltay
           dyt(i,j)=1./deltay
           dym(i,j)=1./deltay
        enddo
     enddo
  endif
end subroutine grdspc


!     ******************************************************************
!
      SUBROUTINE PSTOXY(X,Y,PLA,PLO,ERAD)
!
!     This program convert polar stereographic coordinates to x,y ditto
!     longitude:   0 - 360  ; positive to the east
!     latitude : -90 -  90  ; positive for northern hemisphere
!     it is assumed that the x-axis point towards the east and
!     corresponds to longitude = 0
!
!     TSP 20/06-89
!
!     constants and functions
!            
      FAC(PLA) = ERAD*2.0/(1.0+SIN(PLA*PI180))*COS(PLA*PI180)
      XC(PLA,PLO) = FAC(PLA)*COS(PLO*PI180)
      YC(PLA,PLO) = FAC(PLA)*SIN(PLO*PI180)      
      PI180=3.14159/180.0
!
!     Calculate the coordinates
!
      X = XC(PLA,PLO)
      Y = YC(PLA,PLO)
!
      RETURN
      END
!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE PSTOGE(PLA,PLO,GLAT,GLON,RLAT,WLON1)
!
!     Convert polar stereographic coordinates to geographical lat/lon
!     ditto with the pol.ste. pole at rlat,wlon1 (these names are
!     used to be compatible to the ramsin-parameters)
!     longitude:   0 ; 360 positive east (on input)
!               -180 ; 180 positive east (on output)
!     latitude : -90 ;  90 posive on northern hemisphere
!     It is assumed that the polar stereographic coordinates have been
!     rotated to the standard format with 0 degrees longitude along wlon1
!
!     TSP 21 JUNE 89
!
!     set flag for n/s hemisphere
!
      PI180=3.14159/180.
!      
      IF(RLAT.GE.0.0) THEN
         HSIGN= 1.0
      ELSE
         HSIGN=-1.0
      END IF
!
!     test for a n/s pole case
!
      IF(RLAT.EQ.90.0) THEN
	 GLAT=PLA
         GLON=MOD(PLO+WLON1,360.0)
         GO TO 2000
      END IF
      IF(RLAT.EQ.-90.0) THEN
         GLAT=-PLA
         GLON=MOD(PLO+WLON1,360.0)
         GO TO 2000
      END IF
!
!     test for longitude on 'greenwich or date line'
!
      IF(PLO.EQ.0) THEN
         GLAT=RLAT-90.0+PLA
         IF(GLAT.LT.-90.0) THEN
            GLAT=-180.0-GLAT
            GLON=MOD(WLON1+180.0,360.0)
         ELSE
            GLON=WLON1
         END IF
         GO TO 2000
      END IF      
      IF(PLO.EQ.180.0) THEN
         GLAT=RLAT+90.0-PLA
         IF(GLAT.GT.90.0) THEN
            GLAT=180.0-GLAT
            GLON=MOD(WLON1+180.0,360.0)
         ELSE
            GLON=WLON1
         END IF
         GO TO 2000         
      END IF
!
!     Determine longitude distance relative to wlon1 so it belongs to
!     the absolute interval 0 - 180
!
      ARGU1=PLO
      IF(PLO.GT.180.0) ARGU1 = PLO-360.0
!
!     Get the latitude, the help circle BB and the longitude by first
!     calculating the argument and legalize it - then take the inverse fct.
!
      IF(HSIGN.GT.0.0) THEN
         ARG2A = SIN(PLA*PI180)*SIN(HSIGN*RLAT*PI180)+ &
             COS(PLA*PI180)*COS(RLAT*PI180)*COS((180.0-ARGU1)*PI180)
      ELSE
        ARG2A = SIN(PLA*PI180)*SIN(HSIGN*RLAT*PI180)+ &
             COS(PLA*PI180)*COS(RLAT*PI180)*COS(ARGU1*PI180)
      END IF
      ARG2A = MIN(ARG2A, 1.0)
      ARG2A = MAX(ARG2A,-1.0)
      GLAT  = HSIGN*ASIN(ARG2A)
!
      IF(HSIGN.GT.0.0) THEN
         ARG2A = COS(RLAT*PI180)*SIN(PLA*PI180)+ &
             SIN(RLAT*PI180)*COS(PLA*PI180)*COS(ARGU1*PI180)
      ELSE
        ARG2A = COS(RLAT*PI180)*SIN(PLA*PI180)+ &
            SIN(-RLAT*PI180)*COS(PLA*PI180)*COS((180.0-ARGU1)*PI180)
      END IF
      ARG2A = MIN(ARG2A, 1.0)
      ARG2A = MAX(ARG2A,-1.0)      
      BB    = ACOS(ARG2A)
!
      ARG2A = COS(GLAT)*COS(BB)/(1.0-SIN(GLAT)**2)
      ARG2A = MIN(ARG2A, 1.0)
      ARG2A = MAX(ARG2A,-1.0)      
      GLON  = ACOS(ARG2A)
!     
!     convert the radians to degrees 
!
        GLAT = GLAT/PI180
        GLON = GLON/PI180
!
!       the solution is symmetric so the direction must be if'ed
!
        IF(ARGU1.LT.0.0) THEN
           GLON = 360.0-GLON
        END IF
        GLON=AMOD(GLON+WLON1,360.0)
!
 2000 CONTINUE
!
!     the resultant longitude must be in the interval from -180, 180
!      
      IF(GLON.GT.180.0) GLON=GLON-360.0
      RETURN
      END
!
!     ******************************************************************
!
!     ****************************************************************
!
      SUBROUTINE GETOPS(PLA,PLO,GLAT,GLON,RLAT,WLON1)
!
!     Convert geographical lat/lon coordinates to polar stereographic
!     ditto with the pol.ste. pole at RLAT,WLON1 (these names are
!     used to be compatible to the RAMSIN-parameters)
!     longitude:-180 ; 180 positive east (on input)
!              :   0 ; 360 positive east (on output)
!     latitude : -90 ;  90 posive on northern hemisphere
!     The result is rotated 270 degrees relative to 'standard pol.ste.'
!     WLON1 is defined in the same way as the input
!     approach so as to get the x-axis to point towards the east, and the
!     y-axis towards the north along 0 degrees (at NP south along 180)
!
!     TSP 20/06-89
      double precision pi180,c1,c2,c3,c4,c5,c6,c7,arg2a,bb,pla1,alpha &
        ,plo1,pla90,argu2
!
!     constants
!
      c1=1.
      PI180 = dasin(c1)/90.
!
!     Set flag for N/S hemisphere and convert longitude to <0 ; 360> interval
!
      IF(RLAT.GE.0.0) THEN
         HSIGN= 1.0
      ELSE
         HSIGN=-1.0
      END IF
      GLOR=GLON
      IF(GLOR.LT.0.0) GLOR=360.0+GLOR
      RWLON1=WLON1
      IF(RWLON1.LT.0.0) RWLON1=360.0+WLON1
!
!     Test for a N/S pole case
!
      IF(RLAT.EQ.90.0) THEN
         PLA=GLAT
         PLO=AMOD(GLOR+270.0,360.0)
         GO TO 2000
      END IF
      IF(RLAT.EQ.-90.0) THEN
         PLA=-GLAT
         PLO=AMOD(GLOR+270.0,360.0)
         GO TO 2000
      END IF
!
!     Test for longitude on 'Greenwich or date line'
!
      IF(GLOR.EQ.RWLON1) THEN
         IF(GLAT.GT.RLAT) THEN
            PLA=90.0-GLAT+RLAT
            PLO=90.0
         ELSE
            PLA=90.0-RLAT+GLAT
            PLO=270.0
         END IF
         GO TO 2000
      END IF      
      IF(AMOD(GLOR+180.0,360.0).EQ.RWLON1) THEN
         PLA=RLAT-90.0+GLAT
         IF(PLA.LT.-90.0) THEN
            PLA=-180.0-PLA
            PLO=270.0
         ELSE
            PLO= 90.0
         END IF
         GO TO 2000         
      END IF
!
!     Determine longitude distance relative to RWLON1 so it belongs to
!     the absolute interval 0 - 180
!
      ARGU1 = GLOR-RWLON1
      IF(ARGU1.GT. 180.0) ARGU1 = ARGU1-360.0
      IF(ARGU1.LT.-180.0) ARGU1 = ARGU1+360.0
!
!     1. Get the help circle BB and angle ALPHA (legalize arguments)
!
      c2=glat*pi180
      c3=argu1*pi180
      ARG2A = dCOS(c2)*dCOS(c3)
      ARG2A = dMAX1(ARG2A,-c1)
      ARG2A = dMIN1(ARG2A, c1)         
      BB    = dACOS(ARG2A)
!
      c4=hsign*glat*pi180
      ARG2A = dSIN(c4)/dSIN(BB)
      ARG2A = dMAX1(ARG2A,-c1)
      ARG2A = dMIN1(ARG2A, c1)
      ALPHA = dASIN(ARG2A)
!
!     2. Get PLA and PLO (still legalizing arguments)
!
      c5=rlat*pi180
      c6=hsign*rlat*pi180
      ARG2A = dCOS(c5)*dCOS(BB)+ &
             dSIN(c6)*dSIN(c4)
      ARG2A = dMAX1(ARG2A,-c1)
      ARG2A = dMIN1(ARG2A, c1)         
      PLA1   = dASIN(ARG2A)
!
      ARG2A = dSIN(BB)*dCOS(ALPHA)/dCOS(PLA1)
      ARG2A = dMAX1(ARG2A,-c1)
      ARG2A = dMIN1(ARG2A, c1)
      PLO1   = dASIN(ARG2A)
!
!    Test for passage of the 90 degree longitude (duallity in PLO)
!         Get PLA for which PLO=90 when GLAT is the latitude
!
      ARG2A = dSIN(c4)/dSIN(c6)
      ARG2A = dMAX1(ARG2A,-c1)
      ARG2A = dMIN1(ARG2A, c1)         
      PLA90 = dASIN(ARG2A)
!
!         Get help arc BB and angle ALPHA
!
      ARG2A = dCOS(c5)*dSIN(PLA90)
      ARG2A = dMAX1(ARG2A,-c1)
      ARG2A = dMIN1(ARG2A, c1)
      BB    = dACOS(ARG2A)

      ARG2A = dSIN(c4)/dSIN(BB)
      ARG2A = dMAX1(ARG2A,-c1)
      ARG2A = dMIN1(ARG2A, c1)        
      ALPHA = dASIN(ARG2A)
!
!         Get GLOLIM - it is nesc. to test for the existence of solution
!
      ARGU2  = dCOS(c2)*dCOS(BB)/ &
                 (1.-dSIN(c4)*dSIN(BB)*dSIN(ALPHA))
      IF(dABS(ARGU2).GT.c1) THEN
      GLOLIM = 999.0
      ELSE
        GLOLIM = dACOS(ARGU2)/PI180
      END IF
!
!     Modify (if nesc.) the PLO solution
!
      IF((ABS(ARGU1).GT.GLOLIM.AND.GLAT.LE.RLAT).OR. &
        GLAT.GT.RLAT) THEN
            PLO1 = PI180*180.0 - PLO1
      END IF
!
!     The solution is symmetric so the direction must be if'ed
!
      IF(ARGU1.LT.0.0) THEN
         PLO1 = -PLO1
      END IF
!
!     Convert the radians to degrees
!
      PLA = PLA1/PI180        
      PLO = PLO1/PI180
!
!     To obtain a rotated value (ie so x-axis in pol.ste. points east)
!     add 270 to longitude
!
      PLO=AMOD(PLO+270.0,360.0)
!
 2000 CONTINUE      
      RETURN
      END                                  

