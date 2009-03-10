!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

subroutine grad(m1,m2,m3,ia,iz,ja,jz  &
     ,vc3da,vc3db,dir,gpnt)
  implicit none

  integer, INTENT(IN) :: m1   &
                       , m2   &
                       , m3   &
                       , ia   &
                       , iz   &
                       , ja   &
                       , jz

  real, INTENT(IN)    :: vc3da(m1,m2,m3)
                  
  real, INTENT(INOUT) :: vc3db(m1,m2,m3)

  character(len=*), INTENT(IN) :: dir,gpnt

  character(len=6) :: optyp

  optyp='GRADNT'

  call rams_grad(m1,m2,m3,ia,iz,ja,jz,VC3DA,VC3DB,DIR,GPNT,optyp)

  return
end subroutine grad

!------------------------------------------------------------------------------------

subroutine divcart(m1,m2,m3,ia,iz,ja,jz,vc3da,vc3db,dir,gpnt)
  implicit none

  integer, INTENT(IN) :: m1   &
                       , m2   &
                       , m3   &
                       , ia   &
                       , iz   &
                       , ja   &
                       , jz

  real, INTENT(IN)    :: vc3da(m1,m2,m3)

  real, INTENT(INOUT) :: vc3db(m1,m2,m3)

  character(len=*), INTENT(IN) :: dir,gpnt

  character(len=6) :: optyp

  optyp='DIVCRT'

  call rams_grad(m1,m2,m3,ia,iz,ja,jz,VC3DA,VC3DB,DIR,GPNT,optyp)

  return
end subroutine divcart

!------------------------------------------------------------------------------------

subroutine divstar(m1,m2,m3,ia,iz,ja,jz,vc3da,vc3db,dir,gpnt)
  implicit none

  integer, INTENT(IN) :: m1    &
                       , m2    &
                       , m3    &
                       , ia    &
                       , iz    &
                       , ja    &
                       , jz

  real, INTENT(IN)    :: vc3da(m1,m2,m3)

  real, INTENT(INOUT) :: vc3db(m1,m2,m3)

  character(len=*), INTENT(IN) :: dir,gpnt

  character(len=6) :: optyp

  optyp='DIVSTR'

  call rams_grad(m1,m2,m3,ia,iz,ja,jz,VC3DA,VC3DB,DIR,GPNT,optyp)

  return
end subroutine divstar

!------------------------------------------------------------------------------------

subroutine rams_grad(m1,m2,m3,ia,iz,ja,jz,vc3da,vc3db,dir,gpnt,optyp)

  use mem_grid, only : jdim   &  !INTENT(IN)
                     , dzt    &  !INTENT(IN)
                     , hw     &  !INTENT(IN)
                     , ht     &  !INTENT(IN)
                     , dzm    &  !INTENT(IN)
                     , ngrid  &  !INTENT(IN)
                     , grid_g    !INTENT(IN)


  use mem_scratch, only : vctr1    &   ! INTENT(INOUT)
                        , vctr2        ! INTENT(INOUT)

  implicit none

  integer, INTENT(IN) :: m1    &
                       , m2    &
                       , m3    &
                       , ia    &
                       , iz    &
                       , ja    &
                       , jz

  real, INTENT(IN)    :: vc3da(m1,m2,m3)

  real, INTENT(INOUT) :: vc3db(m1,m2,m3)

  character(len=*), INTENT(IN) :: dir,gpnt

  character(len=6), INTENT(IN) :: optyp

  integer :: jaa,jzz

  jaa=ja
  jzz=jz
  if(jdim.eq.0) then
     jaa=1
     jzz=1
  endif

  IF(DIR.EQ.'XDIR')THEN
     IF(GPNT.EQ.'UPNT')THEN
        CALL GRADXU(m1,m2,m3,ia,iz,jaa,jzz                &
             ,OPTYP,VC3DA,VC3DB,VCTR1,GRID_G(NGRID)%RTGU  &
             ,GRID_G(NGRID)%RTGT,GRID_G(NGRID)%DXT,DZT    &
             ,GRID_G(NGRID)%FMAPUI,GRID_G(NGRID)%FMAPT    &
             ,GRID_G(NGRID)%F13T,HW,VCTR2,'T',JDIM)
     ELSEIF(GPNT.EQ.'VPNT')THEN
        CALL GRADXT(m1,m2,m3,ia,iz,jaa,jzz  &
             ,OPTYP,VC3DA,VC3DB,VCTR1,GRID_G(NGRID)%RTGV  &
             ,GRID_G(NGRID)%RTGM,GRID_G(NGRID)%DXM,DZT  &
             ,GRID_G(NGRID)%FMAPVI,GRID_G(NGRID)%FMAPM  &
             ,GRID_G(NGRID)%F13M  &
             ,HW,VCTR2,'T',JDIM)
     ELSEIF(GPNT.EQ.'WPNT')THEN
        CALL GRADXT(m1,m2,m3,ia,iz,jaa,jzz  &
             ,OPTYP,VC3DA,VC3DB,VCTR1,GRID_G(NGRID)%RTGT  &
             ,GRID_G(NGRID)%RTGU,GRID_G(NGRID)%DXU,DZM  &
             ,GRID_G(NGRID)%FMAPTI,GRID_G(NGRID)%FMAPU  &
             ,GRID_G(NGRID)%F13U  &
             ,HT,VCTR2,'W',JDIM)
     ELSEIF(GPNT.EQ.'TPNT')THEN
        CALL GRADXT(m1,m2,m3,ia,iz,jaa,jzz  &
             ,OPTYP,VC3DA,VC3DB,VCTR1,GRID_G(NGRID)%RTGT  &
             ,GRID_G(NGRID)%RTGU,GRID_G(NGRID)%DXU,DZT  &
             ,GRID_G(NGRID)%FMAPTI,GRID_G(NGRID)%FMAPU  &
             ,GRID_G(NGRID)%F13U  &
             ,HW,VCTR2,'T',JDIM)
     ELSEIF(GPNT.EQ.'NPNT')THEN
        CALL GRADXT(m1,m2,m3,ia,iz,jaa,jzz  &
             ,OPTYP,VC3DA,VC3DB,VCTR1,GRID_G(NGRID)%RTGV  &
             ,GRID_G(NGRID)%RTGM,GRID_G(NGRID)%DXM,DZM  &
             ,GRID_G(NGRID)%FMAPVI,GRID_G(NGRID)%FMAPM  &
             ,GRID_G(NGRID)%F13M  &
             ,HT,VCTR2,'W',JDIM)
     ELSEIF(GPNT.EQ.'OPNT')THEN
        CALL GRADXU(m1,m2,m3,ia,iz,jaa,jzz  &
             ,OPTYP,VC3DA,VC3DB,VCTR1,GRID_G(NGRID)%RTGU  &
             ,GRID_G(NGRID)%RTGT,GRID_G(NGRID)%DXT,DZM  &
             ,GRID_G(NGRID)%FMAPUI,GRID_G(NGRID)%FMAPT  &
             ,GRID_G(NGRID)%F13T  &
             ,HT,VCTR2,'W',JDIM)
     ELSEIF(GPNT.EQ.'PPNT')THEN
        CALL GRADXU(m1,m2,m3,ia,iz,jaa,jzz  &
             ,OPTYP,VC3DA,VC3DB,VCTR1,GRID_G(NGRID)%RTGM  &
             ,GRID_G(NGRID)%RTGV,GRID_G(NGRID)%DXV,DZT  &
             ,GRID_G(NGRID)%FMAPMI,GRID_G(NGRID)%FMAPV  &
             ,GRID_G(NGRID)%F13V  &
             ,HW,VCTR2,'T',JDIM)
     ELSEIF(GPNT.EQ.'MPNT')THEN
        CALL GRADXU(m1,m2,m3,ia,iz,jaa,jzz  &
             ,OPTYP,VC3DA,VC3DB,VCTR1,GRID_G(NGRID)%RTGM  &
             ,GRID_G(NGRID)%RTGV,GRID_G(NGRID)%DXV,DZM  &
             ,GRID_G(NGRID)%FMAPMI,GRID_G(NGRID)%FMAPV  &
             ,GRID_G(NGRID)%F13V  &
             ,HT,VCTR2,'W',JDIM)
     ENDIF
  ELSEIF(DIR.EQ.'YDIR')THEN
     IF(GPNT.EQ.'UPNT')THEN
        CALL GRADYT(m1,m2,m3,ia,iz,jaa,jzz  &
             ,OPTYP,VC3DA,VC3DB,VCTR1,GRID_G(NGRID)%RTGU  &
             ,GRID_G(NGRID)%RTGM,GRID_G(NGRID)%DYM,DZT  &
             ,GRID_G(NGRID)%FMAPUI,GRID_G(NGRID)%FMAPM  &
             ,GRID_G(NGRID)%F23M  &
           ,HW,VCTR2,'T',JDIM)
     ELSEIF(GPNT.EQ.'VPNT')THEN
        CALL GRADYV(m1,m2,m3,ia,iz,jaa,jzz  &
             ,OPTYP,VC3DA,VC3DB,VCTR1,GRID_G(NGRID)%RTGV  &
             ,GRID_G(NGRID)%RTGT,GRID_G(NGRID)%DYT,DZT  &
             ,GRID_G(NGRID)%FMAPVI,GRID_G(NGRID)%FMAPT  &
             ,GRID_G(NGRID)%F23T  &
             ,HW,VCTR2,'T',JDIM)
     ELSEIF(GPNT.EQ.'WPNT')THEN
        CALL GRADYT(m1,m2,m3,ia,iz,jaa,jzz  &
             ,OPTYP,VC3DA,VC3DB,VCTR1,GRID_G(NGRID)%RTGT  &
             ,GRID_G(NGRID)%RTGV,GRID_G(NGRID)%DYV,DZM  &
             ,GRID_G(NGRID)%FMAPTI,GRID_G(NGRID)%FMAPV  &
             ,GRID_G(NGRID)%F23V  &
             ,HT,VCTR2,'W',JDIM)
     ELSEIF(GPNT.EQ.'TPNT')THEN
        CALL GRADYT(m1,m2,m3,ia,iz,jaa,jzz  &
             ,OPTYP,VC3DA,VC3DB,VCTR1,GRID_G(NGRID)%RTGT  &
             ,GRID_G(NGRID)%RTGV,GRID_G(NGRID)%DYV,DZT  &
             ,GRID_G(NGRID)%FMAPTI,GRID_G(NGRID)%FMAPV  &
             ,GRID_G(NGRID)%F23V  &
             ,HW,VCTR2,'T',JDIM)
     ELSEIF(GPNT.EQ.'NPNT')THEN
        CALL GRADYV(m1,m2,m3,ia,iz,jaa,jzz  &
             ,OPTYP,VC3DA,VC3DB,VCTR1,GRID_G(NGRID)%RTGV  &
             ,GRID_G(NGRID)%RTGT,GRID_G(NGRID)%DYT,DZM  &
             ,GRID_G(NGRID)%FMAPVI,GRID_G(NGRID)%FMAPT  &
             ,GRID_G(NGRID)%F23T  &
             ,HT,VCTR2,'W',JDIM)
     ELSEIF(GPNT.EQ.'OPNT')THEN
        CALL GRADYT(m1,m2,m3,ia,iz,jaa,jzz  &
             ,OPTYP,VC3DA,VC3DB,VCTR1,GRID_G(NGRID)%RTGU  &
             ,GRID_G(NGRID)%RTGM,GRID_G(NGRID)%DYM,DZM  &
             ,GRID_G(NGRID)%FMAPUI,GRID_G(NGRID)%FMAPM  &
             ,GRID_G(NGRID)%F23M  &
             ,HT,VCTR2,'W',JDIM)
     ELSEIF(GPNT.EQ.'PPNT')THEN
        CALL GRADYV(m1,m2,m3,ia,iz,jaa,jzz  &
             ,OPTYP,VC3DA,VC3DB,VCTR1,GRID_G(NGRID)%RTGM  &
             ,GRID_G(NGRID)%RTGU,GRID_G(NGRID)%DYU,DZT  &
             ,GRID_G(NGRID)%FMAPMI,GRID_G(NGRID)%FMAPU  &
             ,GRID_G(NGRID)%F23U  &
             ,HW,VCTR2,'T',JDIM)
     ELSEIF(GPNT.EQ.'MPNT')THEN
        CALL GRADYV(m1,m2,m3,ia,iz,jaa,jzz  &
             ,OPTYP,VC3DA,VC3DB,VCTR1,GRID_G(NGRID)%RTGM  &
             ,GRID_G(NGRID)%RTGU,GRID_G(NGRID)%DYU,DZM  &
             ,GRID_G(NGRID)%FMAPMI,GRID_G(NGRID)%FMAPU  &
             ,GRID_G(NGRID)%F23U  &
             ,HT,VCTR2,'W',JDIM)
     ENDIF
  ELSEIF(DIR.EQ.'ZDIR')THEN
     IF(GPNT.EQ.'UPNT')THEN
        CALL GRADZT(m1,m2,m3,ia,iz,jaa,jzz,VC3DA,VC3DB  &
             ,GRID_G(NGRID)%RTGU,DZM)
     ELSEIF(GPNT.EQ.'VPNT')THEN
        CALL GRADZT(m1,m2,m3,ia,iz,jaa,jzz,VC3DA,VC3DB  &
         ,GRID_G(NGRID)%RTGV,DZM)
     ELSEIF(GPNT.EQ.'WPNT')THEN
        CALL GRADZW(m1,m2,m3,ia,iz,jaa,jzz,VC3DA,VC3DB  &
             ,GRID_G(NGRID)%RTGT,DZT)
     ELSEIF(GPNT.EQ.'TPNT')THEN
        CALL GRADZT(m1,m2,m3,ia,iz,jaa,jzz,VC3DA,VC3DB  &
             ,GRID_G(NGRID)%RTGT,DZM)
     ELSEIF(GPNT.EQ.'NPNT')THEN
        CALL GRADZW(m1,m2,m3,ia,iz,jaa,jzz,VC3DA,VC3DB  &
             ,GRID_G(NGRID)%RTGV,DZT)
     ELSEIF(GPNT.EQ.'OPNT')THEN
        CALL GRADZW(m1,m2,m3,ia,iz,jaa,jzz,VC3DA,VC3DB  &
             ,GRID_G(NGRID)%RTGU,DZT)
     ELSEIF(GPNT.EQ.'PPNT')THEN
        CALL GRADZT(m1,m2,m3,ia,iz,jaa,jzz,VC3DA,VC3DB  &
             ,GRID_G(NGRID)%RTGM,DZM)
     ELSEIF(GPNT.EQ.'MPNT')THEN
        CALL GRADZW(m1,m2,m3,ia,iz,jaa,jzz,VC3DA,VC3DB  &
             ,GRID_G(NGRID)%RTGM,DZT)
     ENDIF
  ENDIF

  RETURN
END subroutine rams_grad

!     ******************************************************************
!
!     This is a general subroutine which computes any component of the
!     gradient or divergence of VC3DA and stores it in VC3DB.

subroutine gradxu(m1,m2,m3,ia,iz,ja,jz  &
     ,optyp,vc3da,vc3db,vc1da,rtge,rtgc  &
     ,dx,dz,fmapi,fmap,fq,hq,hq4,lev,jd)

  implicit none

  integer, INTENT(IN) :: m1  &
                        ,m2  &
                        ,m3  &
                        ,ia  &
                        ,iz  &
                        ,ja  &
                        ,jz  &
                        ,jd

  real, INTENT(IN) :: vc3da(m1,m2,m3)   &
                    , rtge(m2,m3)       &
                    , rtgc(m2,m3)       &
                    , dx(m2,m3)         &
                    , fmap(m2,m3)       &
                    , fmapi(m2,m3)      &
                    , dz(m1)            &
                    , fq(m2,m3)         &
                    , hq(*)

  real, INTENT(INOUT)  :: vc3db(m1,m2,m3)   &
                        , vc1da(*)          &
                        , hq4(*)

  character(len=*), INTENT(IN) :: optyp,lev
  
  integer :: i,j,k
  
  IF(OPTYP.EQ.'GRADNT')THEN
     DO J=ja,jz
        DO I=ia,iz
           DO K=1,m1
              VC3DB(K,I,J)=(VC3DA(K,I,J)*RTGE(I,J)  &
                   -VC3DA(K,I-1,J)*RTGE(I-1,J))  &
                   *DX(I,J)/RTGC(I,J)
           ENDDO
        ENDDO
     ENDDO
  ELSE
     DO J=ja,jz
        DO I=ia,iz
           DO K=1,m1
              VC3DB(K,I,J)=(VC3DA(K,I,J)*RTGE(I,J)  &
                   *FMAPI(I,J)  &
                   -VC3DA(K,I-1,J)*RTGE(I-1,J)  &
                   *FMAPI(I-1,J))  &
                   *DX(I,J)/RTGC(I,J)*FMAP(I,J)
           ENDDO
        ENDDO
     ENDDO
  ENDIF

  IF(OPTYP.NE.'DIVSTR')THEN
     IF(LEV.EQ.'W')THEN
        DO K=1,m1
           HQ4(K)=0.25*HQ(K)
        ENDDO
     ELSE
        DO K=2,m1
           HQ4(K)=0.25*HQ(K-1)
        ENDDO
     ENDIF

   ! **(JP)** quebra aninhamento, permitindo vetorizacao
   ! **(JP)** dos ultimos aninhamentos em i,j. O primeiro
   ! **(JP)** continua vetorizado em k. 
!!$   DO J=ja,jz
!!$      DO I=ia,iz
!!$         DO K=2,m1
!!$            VC1DA(K)=HQ4(K)*(VC3DA(K,I,J)+VC3DA(K-1,I,J)  &
!!$                 +VC3DA(K,I-1,J)+VC3DA(K-1,I-1,J))
!!$         ENDDO
!!$         DO K=2,m1-1
!!$            VC3DB(K,I,J)=VC3DB(K,I,J)  &
!!$                 +FQ(I,J)*DZ(K)*(VC1DA(K+1)-VC1DA(K))
!!$         ENDDO
!!$         VC3DB(1,I,J)=VC3DB(2,I,J)
!!$         IF(LEV.EQ.'W')VC3DB(m1-1,I,J)=VC3DB(m1-2,I,J)
!!$         IF(LEV.EQ.'T')VC3DB(m1,I,J)=VC3DB(m1-1,I,J)
!!$      ENDDO
!!$   ENDDO
     do j=ja,jz
        do i=ia,iz
           do k=2,m1
              vc1da(k)=hq4(k)*(vc3da(k,i,j)+vc3da(k-1,i,j)  &
                   +vc3da(k,i-1,j)+vc3da(k-1,i-1,j))
           enddo
           if (optyp /= 'GRADNT') vc1da(2) = 0. 
           do k=2,m1-1
              vc3db(k,i,j)=vc3db(k,i,j)  &
                   +fq(i,j)*dz(k)*(vc1da(k+1)-vc1da(k))
           enddo
        enddo
     enddo
     do j=ja,jz
        do i=ia,iz
           vc3db(1,i,j)=vc3db(2,i,j)
        enddo
     enddo
     if(lev.eq.'W')then
        do j=ja,jz
           do i=ia,iz
              vc3db(m1-1,i,j)=vc3db(m1-2,i,j)
           end do
        end do
     else if(lev.eq.'T')then
        do j=ja,jz
           do i=ia,iz
              vc3db(m1,i,j)=vc3db(m1-1,i,j)
           end do
        end do
     end if
     ! **(JP)** fim modificacao
  ENDIF
  
  RETURN
END subroutine gradxu

SUBROUTINE GRADXT(m1,m2,m3,ia,iz,ja,jz  &
     ,OPTYP,VC3DA,VC3DB,VC1DA,RTGE,RTGC  &
     ,DX,DZ,FMAPI,FMAP,FQ,HQ,HQ4,LEV,JD)

  implicit none

  integer, INTENT(IN) :: m1   &
                       , m2   &
                       , m3   &
                       , ia   &
                       , iz   & 
                       , ja   &
                       , jz   &
                       , jd

  real, INTENT(IN)    :: VC3DA(m1,m2,m3)  &
                       , RTGE(m2,m3)      &
                       , RTGC(m2,m3)      &
                       , DX(m2,m3)        &
                       , FMAP(m2,m3)      &
                       , FMAPI(m2,m3)     &
                       , DZ(*)            &
                       , FQ(m2,m3)        &
                       , HQ(*)

  real, INTENT(INOUT) :: VC1DA(*)         &
                       , VC3DB(m1,m2,m3)  &
                       , HQ4(*)           

  character(len=*), INTENT(IN) :: OPTYP,LEV
  

  integer :: i,j,k
  
  IF(OPTYP.EQ.'GRADNT')THEN
     DO J=ja,jz
        DO I=ia,iz
           DO K=1,m1
              VC3DB(K,I,J)=(VC3DA(K,I+1,J)*RTGE(I+1,J)  &
                   -VC3DA(K,I,J)*RTGE(I,J))  &
                   *DX(I,J)/RTGC(I,J)
           ENDDO
        ENDDO
     ENDDO
  ELSE
     DO J=ja,jz
        DO I=ia,iz
           DO K=1,m1
              VC3DB(K,I,J)=(VC3DA(K,I+1,J)*RTGE(I+1,J)  &
                   *FMAPI(I+1,J)  &
                   -VC3DA(K,I,J)*RTGE(I,J)  &
                   *FMAPI(I,J))  &
                   *DX(I,J)/RTGC(I,J)*FMAP(I,J)
           ENDDO
        ENDDO
     ENDDO
  ENDIF

  IF(OPTYP.NE.'DIVSTR')THEN
     IF(LEV.EQ.'W')THEN
        DO K=1,m1
           HQ4(K)=0.25*HQ(K)
        ENDDO
     ELSE
        DO K=2,m1
           HQ4(K)=0.25*HQ(K-1)
        ENDDO
     ENDIF

   ! **(JP)** quebra aninhamento, permitindo vetorizacao
   ! **(JP)** dos ultimos aninhamentos em i,j. O primeiro
   ! **(JP)** continua vetorizado em k. 
!!$   DO J=ja,jz
!!$      DO I=ia,iz
!!$         DO K=2,m1
!!$            VC1DA(K)=HQ4(K)*(VC3DA(K,I,J)+VC3DA(K-1,I,J)  &
!!$                 +VC3DA(K,I+1,J)+VC3DA(K-1,I+1,J))
!!$         ENDDO
!!$         DO K=2,m1-1
!!$            VC3DB(K,I,J)=VC3DB(K,I,J)  &
!!$                 +FQ(I,J)*DZ(K)*(VC1DA(K+1)-VC1DA(K))
!!$         ENDDO
!!$         VC3DB(1,I,J)=VC3DB(2,I,J)
!!$         IF(LEV.EQ.'W')VC3DB(m1-1,I,J)=VC3DB(m1-2,I,J)
!!$         IF(LEV.EQ.'T')VC3DB(m1,I,J)=VC3DB(m1-1,I,J)
!!$      ENDDO
!!$   ENDDO
     do j=ja,jz
        do i=ia,iz
           do k=2,m1
              vc1da(k)=hq4(k)*(vc3da(k,i,j)+vc3da(k-1,i,j)  &
                   +vc3da(k,i+1,j)+vc3da(k-1,i+1,j))
           enddo
           if (optyp /= 'GRADNT') vc1da(2) = 0. 
           do k=2,m1-1
              vc3db(k,i,j)=vc3db(k,i,j)  &
                   +fq(i,j)*dz(k)*(vc1da(k+1)-vc1da(k))
           enddo
        enddo
     enddo
     do j=ja,jz
        do i=ia,iz
           vc3db(1,i,j)=vc3db(2,i,j)
        enddo
     enddo
     if (lev.eq.'W') then
        do j=ja,jz
           do i=ia,iz
              vc3db(m1-1,i,j)=vc3db(m1-2,i,j)
           enddo
        enddo
     else if(lev.eq.'T') then
        do j=ja,jz
           do i=ia,iz
              vc3db(m1,i,j)=vc3db(m1-1,i,j)
           enddo
        enddo
     end if
     ! **(JP)** fim modificacao
  ENDIF
  
  RETURN
END SUBROUTINE GRADXT

!

SUBROUTINE GRADYV(m1,m2,m3,ia,iz,ja,jz  &
     ,OPTYP,VC3DA,VC3DB,VC1DA,RTGE,RTGC  &
     ,DY,DZ,FMAPI,FMAP,FQ,HQ,HQ4,LEV,JD)

  implicit none

  integer, INTENT(IN) :: m1   &
                       , m2   &
                       , m3   &
                       , ia   &
                       , iz   &
                       , ja   &
                       , jz   &
                       , jd

  real, INTENT(IN)  :: VC3DA(m1,m2,m3)   &
                     , RTGE(m2,m3)       &
                     , RTGC(m2,m3)       &
                     , DY(m2,m3)         &
                     , FMAP(m2,m3)       &
                     , FMAPI(m2,m3)      &
                     , DZ(*)             &
                     , FQ(m2,m3)         &
                     , HQ(*)

  real, INTENT(INOUT) :: VC3DB(m1,m2,m3)   &
                       , HQ4(*)            &
                       , VC1DA(*)


  character(len=*) :: OPTYP,LEV

  integer :: i,j,k

  IF(OPTYP.EQ.'GRADNT')THEN
     DO J=ja,jz
        DO I=ia,iz
           DO K=1,m1
              VC3DB(K,I,J)=(VC3DA(K,I,J)*RTGE(I,J)  &
                   -VC3DA(K,I,J-jd)*RTGE(I,J-jd))  &
                   *DY(I,J)/RTGC(I,J)
           ENDDO
        ENDDO
     ENDDO
  ELSE
     DO J=ja,jz
        DO I=ia,iz
           DO K=1,m1
              VC3DB(K,I,J)=(VC3DA(K,I,J)*RTGE(I,J)  &
                   *FMAPI(I,J)  &
                   -VC3DA(K,I,J-jd)*RTGE(I,J-jd)  &
                   *FMAPI(I,J-jd))  &
                   *DY(I,J)/RTGC(I,J)*FMAP(I,J)
           ENDDO
        ENDDO
     ENDDO
  ENDIF

  IF(OPTYP.NE.'DIVSTR')THEN
     IF(LEV.EQ.'W')THEN
        DO K=1,m1
           HQ4(K)=0.25*HQ(K)
        ENDDO
     ELSE
        DO K=2,m1
           HQ4(K)=0.25*HQ(K-1)
        ENDDO
     ENDIF

   ! **(JP)** quebra aninhamento, permitindo vetorizacao
   ! **(JP)** dos ultimos aninhamentos em i,j. O primeiro
   ! **(JP)** continua vetorizado em k. 
!!$   DO J=ja,jz
!!$      DO I=ia,iz
!!$         DO K=2,m1
!!$            VC1DA(K)=HQ4(K)*(VC3DA(K,I,J)+VC3DA(K-1,I,J)  &
!!$                 +VC3DA(K,I,J-jd)+VC3DA(K-1,I,J-jd))
!!$         ENDDO
!!$         DO K=2,m1-1
!!$            VC3DB(K,I,J)=VC3DB(K,I,J)  &
!!$                 +FQ(I,J)*DZ(K)*(VC1DA(K+1)-VC1DA(K))
!!$         ENDDO
!!$         VC3DB(1,I,J)=VC3DB(2,I,J)
!!$         IF(LEV.EQ.'W')VC3DB(m1-1,I,J)=VC3DB(m1-2,I,J)
!!$         IF(LEV.EQ.'T')VC3DB(m1,I,J)=VC3DB(m1-1,I,J)
!!$      ENDDO
!!$   ENDDO
     do j=ja,jz
        do i=ia,iz
           do k=2,m1
              vc1da(k)=hq4(k)*(vc3da(k,i,j)+vc3da(k-1,i,j)  &
                   +vc3da(k,i,j-jd)+vc3da(k-1,i,j-jd))
           enddo
           if (optyp /= 'GRADNT') vc1da(2) = 0. 
           do k=2,m1-1
              vc3db(k,i,j)=vc3db(k,i,j)  &
                   +fq(i,j)*dz(k)*(vc1da(k+1)-vc1da(k))
           enddo
        enddo
     enddo
     do j=ja,jz
        do i=ia,iz
           vc3db(1,i,j)=vc3db(2,i,j)
        enddo
     enddo
     if (lev.eq.'W') then
        do j=ja,jz
           do i=ia,iz
              vc3db(m1-1,i,j)=vc3db(m1-2,i,j)
           enddo
        enddo
     else if (lev.eq.'T') then
        do j=ja,jz
           do i=ia,iz
              vc3db(m1,i,j)=vc3db(m1-1,i,j)
           enddo
        enddo
     end if
     ! **(JP)** fim modificacao
  ENDIF

  RETURN
END SUBROUTINE GRADYV

SUBROUTINE GRADYT(m1,m2,m3,ia,iz,ja,jz  &
     ,OPTYP,VC3DA,VC3DB,VC1DA,RTGE,RTGC  &
     ,DY,DZ,FMAPI,FMAP,FQ,HQ,HQ4,LEV,JD)

 implicit none

  integer, INTENT(IN) :: m1    &
                       , m2    &
                       , m3    &
                       , ia    &
                       , iz    &
                       , ja    &
                       , jz    &
                       , jd

  real, INTENT(IN) :: VC3DA(m1,m2,m3)   &
                    , RTGE(m2,m3)       &
                    , RTGC(m2,m3)       &
                    , DY(m2,m3)         &
                    , FMAP(m2,m3)       &
                    , FMAPI(m2,m3)      &
                    , DZ(*)             &
                    , FQ(m2,m3)         &
                    , HQ(*)

  real, INTENT(INOUT) :: VC3DB(m1,m2,m3)   &
                       , HQ4(*)            &
                       , VC1DA(*)          

  character(len=*), INTENT(IN) :: OPTYP, LEV
  
  integer :: i,j,k

  IF(OPTYP.EQ.'GRADNT')THEN
     DO J=ja,jz
        DO I=ia,iz
           DO K=1,m1
              VC3DB(K,I,J)=(VC3DA(K,I,J+jd)*RTGE(I,J+jd)  &
                   -VC3DA(K,I,J)*RTGE(I,J))  &
                   *DY(I,J)/RTGC(I,J)
           ENDDO
        ENDDO
     ENDDO
  ELSE
     DO J=ja,jz
        DO I=ia,iz
           DO K=1,m1
              VC3DB(K,I,J)=(VC3DA(K,I,J+jd)*RTGE(I,J+jd)  &
                   *FMAPI(I,J+jd)  &
                   -VC3DA(K,I,J)*RTGE(I,J)  &
                   *FMAPI(I,J))  &
                   *DY(I,J)/RTGC(I,J)*FMAP(I,J)
           ENDDO
        ENDDO
     ENDDO
  ENDIF

  IF(OPTYP.NE.'DIVSTR')THEN
     IF(LEV.EQ.'W')THEN
        DO K=1,m1
           HQ4(K)=0.25*HQ(K)
        ENDDO
     ELSE
        DO K=2,m1
           HQ4(K)=0.25*HQ(K-1)
        ENDDO
     ENDIF

   ! **(JP)** quebra aninhamento, permitindo vetorizacao
   ! **(JP)** dos ultimos aninhamentos em i,j. O primeiro
   ! **(JP)** continua vetorizado em k. 
!!$   DO J=ja,jz
!!$      DO I=ia,iz
!!$         DO K=2,m1
!!$            VC1DA(K)=HQ4(K)*(VC3DA(K,I,J)+VC3DA(K-1,I,J)  &
!!$                 +VC3DA(K,I,J+jd)+VC3DA(K-1,I,J+jd))
!!$         ENDDO
!!$         DO K=2,m1-1
!!$            VC3DB(K,I,J)=VC3DB(K,I,J)  &
!!$                 +FQ(I,J)*DZ(K)*(VC1DA(K+1)-VC1DA(K))
!!$         ENDDO
!!$         VC3DB(1,I,J)=VC3DB(2,I,J)
!!$         IF(LEV.EQ.'W')VC3DB(m1-1,I,J)=VC3DB(m1-2,I,J)
!!$         IF(LEV.EQ.'T')VC3DB(m1,I,J)=VC3DB(m1-1,I,J)
!!$      ENDDO
!!$   ENDDO
     do j=ja,jz
        do i=ia,iz
           do k=2,m1
              vc1da(k)=hq4(k)*(vc3da(k,i,j)+vc3da(k-1,i,j)  &
                   +vc3da(k,i,j+jd)+vc3da(k-1,i,j+jd))
           enddo
           if (optyp /= 'GRADNT') vc1da(2) = 0. 
           do k=2,m1-1
              vc3db(k,i,j)=vc3db(k,i,j)  &
                   +fq(i,j)*dz(k)*(vc1da(k+1)-vc1da(k))
           enddo
        enddo
     enddo
     do j=ja,jz
        do i=ia,iz
           vc3db(1,i,j)=vc3db(2,i,j)
        enddo
     enddo
     if (lev.eq.'W') then
        do j=ja,jz
           do i=ia,iz
              vc3db(m1-1,i,j)=vc3db(m1-2,i,j)
           enddo
        enddo
     else if (lev.eq.'T') then
        do j=ja,jz
           do i=ia,iz
              vc3db(m1,i,j)=vc3db(m1-1,i,j)
           enddo
        enddo
     end if
     ! **(JP)** fim modificacao
  ENDIF
  
  RETURN
END SUBROUTINE GRADYT

SUBROUTINE GRADZW(m1,m2,m3,ia,iz,ja,jz,VC3DA,VC3DB,RTGC,DZ)

  implicit none

  integer, INTENT(IN) :: m1  &
                       , m2  &
                       , m3  &
                       , ia  &
                       , iz  &
                       , ja  &
                       , jz

  real, INTENT(IN)    :: VC3DA(m1,m2,m3)  &
                       , RTGC(m2,m3)      &
                       , DZ(*)

  real, INTENT(INOUT) :: VC3DB(m1,m2,m3) 
  
  integer :: i,j,k

  DO J=ja,jz
     DO I=ia,iz
        DO K=2,m1
           VC3DB(K,I,J)=(VC3DA(K,I,J)-VC3DA(K-1,I,J))*DZ(K)  &
                /RTGC(I,J)
        ENDDO
     ENDDO
  ENDDO
  RETURN
END SUBROUTINE GRADZW

SUBROUTINE GRADZT(m1,m2,m3,ia,iz,ja,jz,VC3DA,VC3DB,RTGC,DZ)

  implicit none

  integer :: m1,m2,m3,ia,iz,ja,jz

  real, INTENT(IN)    :: VC3DA(m1,m2,m3) &
                       , RTGC(m2,m3)     &
                       , DZ(*)

  real, INTENT(INOUT) :: VC3DB(m1,m2,m3)

  integer :: i,j,k

  DO J=ja,jz
     DO I=ia,iz
        DO K=1,m1-1
           VC3DB(K,I,J)=(VC3DA(K+1,I,J)-VC3DA(K,I,J))*DZ(K)  &
                /RTGC(I,J)
        ENDDO
     ENDDO
  ENDDO
  RETURN
END SUBROUTINE GRADZT
