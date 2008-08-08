!##############################################################
!
! Programmed by  Enio Pereira de Souza
! Adapted    by  Alvaro Luiz Fazenda    (for V.5.04)
!
!##############################################################
subroutine souza_cupar_driver()

  ! USE Modules for 5.0
  use mem_basic
  use mem_micro
  use mem_grid
  use mem_turb
  use mem_tend
  use node_mod, only : MXP,   &   ! INTENT(IN)
       MYP,                   &   ! INTENT(IN)
       MZP,                   &   ! INTENT(IN)
       IA,                    &   ! INTENT(IN)
       IZ,                    &   ! INTENT(IN)
       JA,                    &   ! INTENT(IN)
       JZ,                    &   ! INTENT(IN)
       I0,                    &   ! INTENT(IN)  ! Rever função
       J0                         ! INTENT(IN)  ! Rever função

  use mem_cuparm, only : confrq,cuparm_g,nclouds,cptime       ! INTENT(IN)

  implicit none

  ! INCLUDE 'rcommons.h' ! Not necessary in V.5.04
  !      INCLUDE 'rnode.h'  ! Not necessary in V.5.x

  integer :: I, J
  integer :: icld
!  TFZMIN = turb_g(ngrid)%sflux_t(1,1)
!  TFZMAX = turb_g(ngrid)%sflux_t(1,1)
  icld = nclouds ! Just to make it similar to other methods

  !

  call AZERO(mxp*myp*mzp,cuparm_g(ngrid)%thsrc(1,1,1,icld))
  call AZERO(mxp*myp*mzp,cuparm_g(ngrid)%rtsrc(1,1,1,icld))
  call AZERO(mxp*myp,cuparm_g(ngrid)%upmf(1,1,icld))

  call SHCUPAR(mzp,mxp,myp,ia,iz,ja,jz,i0,j0,                   &
       basic_g(ngrid)%wp(1,1,1), basic_g(ngrid)%theta(1,1,1),   &
       basic_g(ngrid)%pp(1,1,1), basic_g(ngrid)%pi0(1,1,1),     &
       basic_g(ngrid)%dn0(1,1,1), basic_g(ngrid)%rv(1,1,1),     &
       cuparm_g(ngrid)%thsrc(1,1,1,icld),                       &
       cuparm_g(ngrid)%rtsrc(1,1,1,icld),                       &
       cuparm_g(ngrid)%upmf(1,1,icld), grid_g(ngrid)%rtgt(1,1), &
       turb_g(ngrid)%sflux_t(1,1),                              & 
       turb_g(ngrid)%sflux_r(1,1),                              &
       turb_g(ngrid)%vkh(1,1,1),                                &
       micro_g(ngrid)%rcp(1,1,1))

       ! turb_g()%vkh     == vkkh
       ! turb_g()%sflux_r == QFZ
       ! turb_g()%sflux_t == TFZ
       ! micro_g()%rcp    == rcp

  call ACCUM(mxp*myp*mzp, tend%tht(1), cuparm_g(ngrid)%thsrc(1,1,1,icld))
  call ACCUM(mxp*myp*mzp, tend%tht(1), cuparm_g(ngrid)%rtsrc(1,1,1,icld))

  return
end subroutine souza_cupar_driver

!##############################################################
subroutine SHCUPAR(m1,m2,m3,ia,iz,ja,jz,i0,j0,                 &
     WP,THETA,PP,PI0,DN0,RV,                                   &
     THSRCSH,RTSRCSH,SHMF,RTGT,TFZ,QFZ,KHV,RCLOUD)

  use conv_coms, only : WCON,   &   ! INTENT(OUT)
       THTCON,                  &   ! INTENT(OUT)
       PICON,                   &   ! INTENT(OUT)
       DNCON,                   &   ! INTENT(OUT)
       ZCON,                    &   ! INTENT(OUT)
       ZZCON,                   &   ! INTENT(OUT)
       ICPRTFL,                 &   ! INTENT(OUT)   ! Maybe local variable
       ICPLTFL,                 &   ! INTENT(OUT)   ! Maybe local variable
       IGO,                     &   ! INTENT(OUT)   ! Maybe local variable
       KLCL                         ! INTENT(IN)
  use shcu_vars_const, only : ENTF,   &   ! INTENT(OUT)
       ALHF,                    &   ! INTENT(OUT)
       QVCON,                   &   ! INTENT(OUT)
       AKVD,                    &   ! INTENT(OUT)
       CL_CON,                  &   ! INTENT(OUT)
       DTDT,                    &   ! INTENT(IN)
       WC,                      &   ! INTENT(IN)
       DRDT                         ! INTENT(IN)

  use mem_grid, only : time,    &   ! INTENT(IN)
       ZT,                      &   !INTENT(IN)
       ZM                           !INTENT(IN/OUT)

  implicit none

  integer, intent(IN) :: m1, m2, m3, ia, iz, ja, jz, i0, j0
  real, intent(IN)    :: WP(m1,m2,m3), THETA(m1,m2,m3),        &
       PP(m1,m2,m3), PI0(m1,m2,m3), DN0(m1,m2,m3),             &
       RV(m1,m2,m3)
  real, intent(OUT)   :: THSRCSH(m1,m2,m3), RTSRCSH(m1,m2,m3), &
       SHMF(m2,m3)
  real, intent(IN)    :: RTGT(m2,m3), TFZ(m2,m3), QFZ(m2,m3),  &
       KHV(m1,m2,m3), RCLOUD(m1,m2,m3)

  integer :: ICPCNT = 0
  integer :: IPRTFRQ, I, J, K

!  TFZMIN = TFZ(1,1)
!  TFZMAX = TFZ(1,1)


  ICPRTFL=0
  IPRTFRQ=8
  ICPLTFL=0
  ICPCNT=ICPCNT+1
  if(mod(ICPCNT-IPRTFRQ+1,IPRTFRQ).eq.0) then
     ICPRTFL=1
  endif

  do J=ja,jz
     do I=ia,iz
        !	        
        !
        ! Definition of the enthalpy flux ENTF (K*m/s) and
        ! latent energy flux (m/s)*(g/g)
        !
        ENTF=min(.50,TFZ(I,J))
        ALHF=min(.00028,QFZ(I,J))
        !
        do K=1,m1
           WCON(K)=WP(K,I,J)
           THTCON(K)=THETA(K,I,J)
           PICON(K)=(PP(K,I,J)+PI0(K,I,J))
           DNCON(K)=DN0(K,I,J) 
           !EPS
           !	At this point, we change rv by qv as required by the
           !	basic equations.
           !EPS            
           QVCON(K) =RV(K,I,J)/(1+RV(K,I,J))            
           ZCON(K)  =ZT(K) *RTGT(I,J)
           ZZCON(K) =ZM(K) *RTGT(I,J)
           AKVD(K)  =KHV(K,I,J)
           CL_CON(K)=RCLOUD(K,I,J)
           THSRCSH(K,I,J)=0.
           RTSRCSH(K,I,J)=0.
        enddo
        !
        if(ENTF.le.0.0) then
           !	  WRITE(11,*)'SH= ',ENTF  
           GO TO 1000
        endif
        !        	   	  
        IGO=1
        !
        call SHCU_ENV(m1-1) 
        !
        if(IGO.ne.0) call CL_TOP
        ! 
        if(IGO.ne.0) call W_SHALLOW(I,J,TIME)
        !
        if(IGO.ne.0) call SH_RATES
        ! 
        if(IGO.ne.0) then
           !
           call SH2MOD(m1)

           do K=2,m1-1
              THSRCSH(K,I,J)=DTDT(K)
              RTSRCSH(K,I,J)=DRDT(K) 
           enddo
           SHMF(i,j) = DNCON(KLCL)*WC(KLCL)
           !
        endif
        !
1000    continue
        !
     enddo
  enddo
  !
  return
end subroutine SHCUPAR
!
!     ******************************************************************
!
subroutine SHCU_ENV(NZ)  

  use conv_coms, only : NKP,    &   ! INTENT(IN)  ! Parameter
       ZC,                      &   ! INTENT(OUT)
       ZE,                      &   ! INTENT(OUT)
       PKE,                     &   ! INTENT(OUT)
       THE,                     &   ! INTENT(IN)
       THVE,                    &   ! INTENT(OUT)
       TE,                      &   ! INTENT(OUT)
       PE,                      &   ! INTENT(OUT)
       RHOE,                    &   ! INTENT(OUT)
       DZLOW,                   &   ! INTENT(OUT) ! Talvez var.local
       DZHIGH,                  &   ! INTENT(OUT) ! Talvez var.local
       ZMID,                    &   ! INTENT(OUT) ! Talvez var.local
       CDZMIN,                  &   ! INTENT(OUT) ! Rever função.
       ZCON,                    &   ! INTENT(IN)
       KMT,                     &   ! INTENT(OUT)
       WCON,                    &   ! INTENT(IN/OUT)
       ZZCON,                   &   ! INTENT(IN)
       WPE,                     &   ! INTENT(IN/OUT)
       THTCON,                  &   ! INTENT(IN)
       PICON,                   &   ! INTENT(IN)
       KCON,                    &   ! INTENT(OUT)  ! Talvez var.local
       TLCL,                    &   ! INTENT(IN/OUT)
       PLCL,                    &   ! INTENT(IN/OUT)
       DZLCL,                   &   ! INTENT(IN/OUT)
       KLCL,                    &   ! INTENT(OUT)
       IGO                          ! INTENT(OUT)
  use shcu_vars_const, only : QVE,    &   ! INTENT(IN/OUT)
       DSE,                     &   ! INTENT(OUT)
       UHE,                     &   ! INTENT(OUT)
       EVAPS,                   &   ! INTENT(OUT) ! Maybe local var.?
       QVSE,                    &   ! INTENT(OUT) ! Maybe local var.?
       UHES,                    &   ! INTENT(OUT)
       RHE,                     &   ! INTENT(OUT) ! **
       GAMMA,                   &   ! INTENT(OUT)
       UHC,                     &   ! INTENT(OUT)
       DELZ,                    &   ! INTENT(OUT)
       DLDZBY2,                 &   ! INTENT(OUT) ! Maybe local var.?
       DSC,                     &   ! INTENT(OUT)
       DSC0,                    &   ! INTENT(OUT)
       QVC,                     &   ! INTENT(OUT)
       WLC,                     &   ! INTENT(OUT)
       QVCON,                   &   ! INTENT(IN/OUT)
       AKVD,                    &   ! INTENT(IN/OUT)
       CL_CON,                  &   ! INTENT(IN/OUT)
       CL_PE,                   &   ! INTENT(IN/OUT)
       G,                       &   ! INTENT(IN)  ! Parameter
       CP,                      &   ! INTENT(IN)  ! Parameter
       CPR,                     &   ! INTENT(IN)  ! Parameter
       P00,                     &   ! INTENT(IN)  ! Parameter
       R,                       &   ! INTENT(IN)  ! Parameter
       KZI,                     &   ! INTENT(OUT) ! Maybe local var.?
       ALVL,                    &   ! INTENT(IN)  ! Parameter
       AKVDE                        ! INTENT(IN/OUT)

  implicit none

  integer, intent(IN) :: NZ

  real :: DLAMB(NKP)  ! , HZ(NKP)  ! not used
  !       Basic constants
  real :: CONST1, CONST2, ES00, EPSLON, UMMEPS, TA0, CONST3,  &
       C0, DLAMB0, ZREF, ZNZ, TLLL, PLLL, ZLLL, RLLL, DZLLL, DZDD, DTHV
  integer :: NKMID, K


!  R=287.
!  CP=1004.
!  RCP=.286
!  CPR=3.4965
!  ALVL=2.5E6
!  ALIV=2.834E6
!  P00=1E5
!  G=9.8
!  AKLV=2340.6
!  AKIV=2825.7
!  ------------All variables above have been defined in USE shcu_vars_const
  DZLOW=200.
  DZHIGH=500.
  !srf      ZMID=3000.
  ZMID=4000.
!  GAMMAD=.00976
  CONST1=17223003.15
  CONST2=29.65
  ES00=611.2
  EPSLON=.622
  UMMEPS=.378
  TA0=273.15
  CONST3=17.67
  !
  CDZMIN=3000.
  !EPS
  C0=0.0
  !_modified in 01/22/99      DLAMB=0.003
  DLAMB0=1.E-06
  ZREF=700.
  !
  !
  !           INTERPOLATE MODEL SOUNDING (ENVIRONMENT) TO HIGHER
  !             RESOLUTION GRID
  !
  NKMID=ZMID/DZLOW+1
  ZC(1)=0.
  do K=2,NKMID
     ZC(K)=ZC(K-1)+DZLOW
  enddo
  do K=NKMID+1,NKP
     ZC(K)=ZC(K-1)+DZHIGH
  enddo
  ZE(1)=0.
  do K=2,NKP
     ZE(K)=(ZC(K)+ZC(K-1))*.5
  enddo
  !                   FIND MODEL TOP ON CONVECTIVE GRID
  ZNZ=ZCON(NZ)
  do K=NKP,1,-1
     if(ZE(K).lt.ZNZ)GO TO 13
  enddo
  stop ' ENVIR STOP 12'
13 continue
  KMT=K
  !                   DO ACTUAL INTERPOLATION
  !
  call HTINT(NZ,WCON,ZZCON,KMT,WPE,ZE)
  call HTINT(NZ,THTCON,ZCON,KMT,THE,ZE)
  call HTINT(NZ,QVCON,ZCON,KMT,QVE,ZE)
  call HTINT(NZ,AKVD,ZCON,KMT,AKVDE,ZE)
  call HTINT(NZ,CL_CON,ZCON,KMT,CL_PE,ZE)
  !
  do K=1,KMT
     QVE(K)=max(QVE(K),1E-8)
  enddo
  !
  !         COMPUTE THETA V, THETA E, AND GET PRESSURE PROFILE
  !
  PKE(1)=PICON(1)
  do K=1,KMT
     THVE(K)=THE(K)*(1.+.61*QVE(K))
  enddo
  !
  do K=2,KMT
     PKE(K)=PKE(K-1)-G*2.*(ZE(K)-ZE(K-1))/(THVE(K)+THVE(K-1))
  enddo
  do K=1,KMT
     TE(K)=THE(K)*PKE(K)/CP
     PE(K)=(PKE(K)/CP)**CPR*P00
     RHOE(K)=PE(K)/(R*TE(K)*(1.+.61*QVE(K)))
  enddo
  !
  !EPS    FIND THE MAIN SOURCE LEVEL OF THE UPDRAFT. WE WILL ASSUME
  !EPS	THAT PARCELS START FROM THE FIRST LAYER DEFINED BY K=2
  !
  KCON=2
  !
  !         FIND THE LCL OF A LAYER AVERAGE AROUND THE SOURCE LEVEL
  !
  !77 CONTINUE !Never used
  !      TLLL=(TE(KCON)+TE(KCON+1)+TE(KCON-1))/3.
  TLLL=(TE(KCON)+TE(KCON-1))/2.      
  PLLL=PE(KCON)
  !      RLLL=(QVE(KCON)+QVE(KCON+1)+QVE(KCON-1))/3.
  RLLL=(QVE(KCON)+QVE(KCON-1))/2.      
  ZLLL=ZE(KCON)
  !
  call LCL(TLLL,PLLL,RLLL,TLCL,PLCL,DZLCL)
  !
  !         FIND THE CLOSEST LEVEL ON THE CONVECTIVE GRID TO THE LCL
  !
  DZLLL=1E20
  do K=1,KMT
     DZDD=abs(ZE(K)-(ZLLL+DZLCL))
     if(DZDD.lt.DZLLL)then
        DZLLL=DZDD
        KLCL=K
     endif
  enddo
  !
  !     Determination of Zi, the PBL top as the level where
  !     the turbulent diffusion coefficient (AKVDE) vanishes. The loop
  !     starts at k=3 because sometimes AKVDE is zero at k=1 or 2.
  !
  !      AKVMIN=1.
  !      DO K=3,KMT
  !        IF(AKVDE(K).LT.AKVMIN.or.abs(cl_pe(k)).gt.1.e-5) THEN
  !          KZI=K-1
  !          GO TO 78
  !        ENDIF
  !      ENDDO
  !
  !      Determination of Zi, the PBL top as the level where
  !      the virtual potential temperature THVE increases by
  !      more than DTHV as compared to the previous level
  !
  DTHV=0.5
  do K=3,KMT
     if((THVE(K)-(THVE(K-1))).gt.DTHV) then
  !eps      KZI=K-1
        KZI=K
        exit
     endif
  enddo
  !
!78 CONTINUE !Not used, commented line 509
  !
  !     Following Wilde et. al. (1985) shallow cumulus onset occurs
  !     when the top of the entrainment zone overlaps the base of 
  !     the LCL zone. 
  !
  !_Comentado em 10/03/99 para teste  de sensibilidade
  if(KZI.lt.KLCL) then
     IGO=0
     return
  endif
  !
  !EPS
  !	determination of environment variables
  !
  do K=1,KMT
     DSE(K)=CP*TE(K)+G*ZE(K)
     UHE(K)=DSE(K)+ALVL*QVE(K)
     EVAPS(K)=ES00*exp(CONST3*(TE(K)-TA0)/(TE(K)-CONST2))
     QVSE(K)=EPSLON*EVAPS(K)/(PE(K)-UMMEPS*EVAPS(K))
     UHES(K)=DSE(K)+ALVL*QVSE(K)
     RHE(K)=QVE(K)/QVSE(K)
     !
     GAMMA(K)=CONST1*PE(K)*QVSE(K)**2/EVAPS(K)
     GAMMA(K)=GAMMA(K)/((TE(K)-CONST2)*(TE(K)-CONST2))
  enddo
  !
  !     calculating the cloud moist static energy profile. We assume
  !     that entrainment occurs only above the LCL.
  !      
  UHC(1)=UHE(1)
  UHC(2)=UHE(2)
  DELZ(2)=ZE(2)-ZE(1)
  !
  !     Option 1: the parcel is entrained since its source level
  !
  !C      DO K=3,KMT
  !C          DLAMB(K)=EXP(LOG(DLAMB0)+2.3*(ZE(K)-ZE(1))/ZREF)            
  !C          DLDZBY2(K)=MIN(1.,DLAMB(K))
  !C          DLDZBY2(K)=DLDZBY2(K)*DELZ(K)/2
  !     
  !C        UHC(K)=(UHC(K-1)-DLDZBY2(K)*(UHC(K-1)-UHE(K)-UHE(K-1)))/  &
       !C       (1+DLDZBY2(K))
  !C      ENDDO
  !CCC
  !
  !     Option 2: the parcel rises without mixing till its lcl and 
  !     starts being entrained from the lcl up to the cloud top.
  !
  if(KLCL.ge.3) then
     do K=3,KLCL
        !      IF(KLCL.GE.2) THEN
        !         DO K=2,KLCL         
        UHC(K)=UHE(2)
        DELZ(K)=ZE(K)-ZE(K-1)
     enddo
     !
     do K=KLCL+1,KMT
        DELZ(K)=ZE(K)-ZE(K-1)
        DLAMB(K)=exp(log(DLAMB0)+2.3*(ZE(K)-ZE(1))/ZREF)            
        DLDZBY2(K)=min(1.,DLAMB(K))
        DLDZBY2(K)=DLDZBY2(K)*DELZ(K)/2
        !     
        UHC(K)=(UHC(K-1)-DLDZBY2(K)*(UHC(K-1)-UHE(K)-UHE(K-1)))   &
             /(1+DLDZBY2(K))
     enddo
  else
     IGO=0
     !        WRITE(11,*)'C',KLCL
     return
  endif
  !
  !	calculating the in-cloud variables
  !
  do K=1,KMT
     DSC(K) =DSE(K)+(UHC(K)-UHES(K))/(1+GAMMA(K))
     DSC0(K)=DSE(K)+(UHE(2)-UHES(K))/(1+GAMMA(K))
     QVC(K) =QVSE(K)+GAMMA(K)*(UHC(K)-UHES(K))/(ALVL*(1+GAMMA(K)))
     WLC(K)=0.0
  enddo
  !
  do K=KLCL+1,KMT    
     WLC(K)=WLC(K-1)-(QVC(K)-QVC(K-1))-DLAMB(K)*              &
          (QVC(K)-QVE(K))*DELZ(K)+(C0-DLAMB(K))*WLC(K-1)*DELZ(K)
     ! FOR DEBUG:
     !* 250 Floating-point data overflow PROG=shcu_env ELN=599(4003dc44c)
     !* 253 Invalid operation PROG=shcu_env ELN=599(4003dc44c)
     !WLC(K)=WLC(K-1)
     !WLC(K)=WLC(K) - (QVC(K)-QVC(K-1))
     !TEMPOR1 = DLAMB(K)*(QVC(K)-QVE(K))
     !TEMPOR1 = TEMPOR1*DELZ(K)
     !WLC(K)=WLC(K) - TEMPOR1
     !WLC(K)=WLC(K) + (C0-DLAMB(K))*WLC(K-1)*DELZ(K)

     !IF (WLC(K)<0) EXIT

     WLC(K)=min(QVC(KLCL),WLC(K))
     WLC(K)=max(.1E-12,WLC(K))

  enddo
  !
  return
end subroutine SHCU_ENV
!
!     ******************************************************************
!
subroutine SH2MOD(m1)

  use conv_coms, only : KMT,  &   ! INTENT(IN)
       QVCT1,                 &   ! INTENT(OUT)
       RHOE,                  &   ! INTENT(IN)
       PKE,                   &   ! INTENT(IN)
       QVCT2,                 &   ! INTENT(OUT)
       QVCT3,                 &   ! INTENT(OUT)
       ZC,                    &   ! INTENT(IN)
       QVCT4,                 &   ! INTENT(OUT)
       ZZCON,                 &   ! INTENT(IN)
       DNCON,                 &   ! INTENT(IN)
       PICON                      ! INTENT(IN)
  use shcu_vars_const, only : DTDT, &   ! INTENT(IN/OUT)
       ALVL,                        &   ! INTENT(IN) ! Parameter
       DRDT                             ! INTENT(IN/OUT)

  use mem_scratch, only : VCTR5,  & ! INTENT(IN/OUT)
       VCTR6                        ! INTENT(IN/OUT)

  implicit none

  integer, intent(IN) :: m1

  ! Local variables
  integer :: K, TFTC, TFRC, TFTM, TFRM, FTRES, FRRES

  ! External Function to Sum a array
  real, external :: ssum

  !        Compute integrated heating and moistening tendencies
  !
  do K=2,KMT
     QVCT1(K)=RHOE(K)*DTDT(K)*PKE(K)
     QVCT2(K)=RHOE(K)*ALVL*DRDT(K)
     QVCT3(K)=(ZC(K)-ZC(K-1))*QVCT1(K)
     QVCT4(K)=(ZC(K)-ZC(K-1))*QVCT2(K)
     !        print*,'0',QVCT1(K),QVCT2(k),k,DRDT(K),zc(k)
  enddo
  TFTC=SSUM(KMT-1,QVCT3(2),1)
  TFRC=SSUM(KMT-1,QVCT4(2),1)
  !
  !         Transfer tendencies to model grid
  !
  !new--------
  call vertmap2(qvct1,zc,kmt,vctr5,zzcon,m1-1)
  call vertmap2(qvct2,zc,kmt,vctr6,zzcon,m1-1)

  do K=2,m1-1
     VCTR5(K)=VCTR5(K)*(ZZCON(K)-ZZCON(K-1))
     VCTR6(K)=VCTR6(K)*(ZZCON(K)-ZZCON(K-1))
  enddo
  !new--------

  !old--------
  !      DO K=1,m1
  !        VCTR5(K)=0.
  !        VCTR6(K)=0.
  !      ENDDO
  !
  !      DZLFT=0.
  !      L=2
  !      DO K=2,m1-1
  !        IF(DZLFT.NE.0.) THEN
  !          VCTR5(K)=VCTR5(K)+QVCT1(L)*DZLFT
  !          VCTR6(K)=VCTR6(K)+QVCT2(L)*DZLFT
  !          L=L+1
  !        ENDIF
  !   60   CONTINUE
  !        IF(ZC(L).LE.ZZCON(K)) THEN
  !          VCTR5(K)=VCTR5(K)+QVCT1(L)*(ZC(L)-ZC(L-1))
  !          VCTR6(K)=VCTR6(K)+QVCT2(L)*(ZC(L)-ZC(L-1))
  !          L=L+1
  !          DZLFT=0.
  !          GO TO 60
  !        ELSE
  !          VCTR5(K)=VCTR5(K)+QVCT1(L)*(ZZCON(K)-ZC(L-1))
  !          VCTR6(K)=VCTR6(K)+QVCT2(L)*(ZZCON(K)-ZC(L-1))
  !          DZLFT=ZC(L)-ZZCON(K)
  !        ENDIF
  !      ENDDO
  !
  !old--------
  !
  !         Make sure the transfer from the convective grid to the model
  !           grid happened correctly.
  !
  TFTM=SSUM(m1-2,VCTR5(2),1)
  TFRM=SSUM(m1-2,VCTR6(2),1)
  !
  FTRES=TFTM-TFTC
  FRRES=TFRM-TFRC
  if(abs(FTRES).gt..01*abs(TFTC)) then
     print*,' Energy error in grid tranfser in convective param.'
     print*,' TFTM,TFTC ',TFTM,TFTC
     print*,'ERRO SHCU'
  endif
  !
  !         Change energy tendencies to temperature and mixing ratio
  !           tendencies.
  !
  do K=2,m1-1
     DTDT(K)=VCTR5(K)/((ZZCON(K)-ZZCON(K-1))*DNCON(K)*PICON(K))
     DRDT(K)=VCTR6(K)/((ZZCON(K)-ZZCON(K-1))*DNCON(K)*ALVL)
     !        print*,'1',DTDT(k),DRDT(k),zzcon(k)
  enddo

  !
  return
end subroutine SH2MOD

!
!     ******************************************************************
!      
subroutine CL_TOP

  use conv_coms, only : NKP,  &   ! INTENT(IN) ! Parameter
       KMT,                   &   ! INTENT(IN)
       TE,                    &   ! INTENT(IN)
       ZE,                    &   ! INTENT(IN)
       KLCL,                  &   ! INTENT(IN)
       IGO                        ! INTENT(OUT)
  use shcu_vars_const, only : DSCV, &   ! INTENT(OUT)
       DSC,                         &   ! INTENT(IN)
       CP,                          &   ! INTENT(IN) ! Parameter
       QVC,                         &   ! INTENT(IN)
       WLC,                         &   ! INTENT(IN)
       DSEV,                        &   ! INTENT(OUT)
       DSE,                         &   ! INTENT(IN)
       QVE,                         &   ! INTENT(IN)
       DSC0V,                       &   ! INTENT(OUT) ! Maybe local var.?
       DSC0,                        &   ! INTENT(IN)
       DSC0VM,                      &   ! INTENT(OUT) ! Maybe local var.?
       ENTF,                        &   ! INTENT(IN)
       G,                           &   ! INTENT(IN)  ! Parameter
       CAPE,                        &   ! INTENT(OUT)
       DELZ,                        &   ! INTENT(IN)
       KTOP                             ! INTENT(OUT)

  implicit none

  !
  !EPS      
  !
  !

  ! Local variables
  real :: DSCVM(NKP), DSEVM(NKP), TEM(NKP), EMP(NKP)
  real :: GAMMAD, BUOY1, BUOY2
  integer :: K, L

  !       Basic constants

!  R=287.
!  CP=1004.5
!  RCP=.286
!  CPR=3.4965
!  ALVL=2.50E6
!  ALIV=2.837E6
!  P00=1E5
!  G=9.806
!  AKLV=2340.6
!  AKIV=2825.7
!  All variables above have been defined in USE shcu_vars_const

  GAMMAD=.00976
  !EPS
  !     Determination of the cloud top based on integrated cloud buoyancy. 
  !     First, determination of virtual static energy profiles 
  !     Sv=S+0.608Cp<T>q-Cp<T>l
  !
  do K=1,KMT
     DSCV(K)=DSC(K)+CP*TE(K)*(.608*QVC(K)-WLC(K))
     DSEV(K)=DSE(K)+.608*CP*TE(K)*QVE(K)
     DSC0V(K)=DSC0(K)+.608*CP*TE(K)*QVE(K)
     !       DSC0V(K)=DSC0(K)+CP*TE(K)*(.608*QVC(K)-WLC(K))
  enddo
  !
  do K=2,KMT
     DSCVM(K)=(DSCV(K)+DSCV(K-1))/2
     DSEVM(K)=(DSEV(K)+DSEV(K-1))/2
     DSC0VM(K)=(DSC0V(K)+DSC0V(K-1))/2
     TEM(K)=(TE(K)+TE(K-1))/2
  enddo
  !
  !     determination of the integrated bouyancy between the surface and 
  !      the LCL (see Albrecht et al. 1986, Eq. A5)
  !     ENTF is the enthalpy flux at surface in (K*m/s)
  !
  if(ENTF.gt.0.0) then
     BUOY1=ENTF*(1+.608*QVE(1))
     BUOY1=G*ZE(KLCL)*BUOY1/TE(1)
     BUOY1=1.3333*BUOY1**.6667
  else
     BUOY1=.0
     IGO=0
     !        WRITE(11,*)'D'
     return
  endif
  !
  !     checking if the parcel is able to sustain positive buoyancy one
  !     level above the LCL
  !
  CAPE=.0
  BUOY2=.0
  L=KLCL+1
  BUOY2=BUOY2+GAMMAD*(DSCVM(L)-DSEVM(L))*DELZ(L)/TEM(L)
  !                           
  if((BUOY1+BUOY2).le.0.0) then
     IGO=0 
     !       WRITE(11,*)'BUOY'      
     return
  endif
  !           
  !     Integration continues till the level of zero buoyancy is found.
  !     The cloud top is assumed to be one level below.
  !
  !      
  do K=KLCL+2,KMT                 
     BUOY2=BUOY2+GAMMAD*(DSCVM(K)-DSEVM(K))*DELZ(K)/TEM(K)
     if((BUOY1+BUOY2).le.0.0) then
        KTOP=K-1
        GO TO 88
     endif
  enddo
88 continue      
  !
  do K=KLCL,KTOP
     EMP(K)=DSC0VM(K)-DSEVM(K)
     EMP(K)=max(0.,EMP(K))
     CAPE=CAPE+GAMMAD*EMP(K)*DELZ(K)/TEM(K)
  enddo
  !
  return
end subroutine CL_TOP
!
!     ******************************************************************
! 
subroutine W_SHALLOW(IP,JP,TIME)

  use conv_coms, only : TE,   &   ! INTENT(IN)
       KLCL,                  &   ! INTENT(IN)
       IGO,                   &   ! INTENT(OUT)
       RHOE,                  &   ! INTENT(IN)
       ZE                         ! INTENT(IN)
  use shcu_vars_const, only : KTOP, &   ! INTENT(IN)
       EFIC,                        &   ! INTENT(OUT) ! Maybe local var.?
       CAPE,                        &   ! INTENT(IN)
       CP,                          &   ! INTENT(IN)
       ENTF,                        &   ! INTENT(IN)  ! Parameter
       ALVL,                        &   ! INTENT(IN)  ! Parameter
       ALHF,                        &   ! INTENT(IN)
       DCAPE,                       &   ! INTENT(OUT) ! Maybe local var.?
       TCAPE,                       &   ! INTENT(OUT) ! Maybe local var.?
       WC                               ! INTENT(OUT)

  implicit none
 
  integer, intent(IN) :: IP, JP
  real(kind=8), intent(IN)    :: TIME

  !      
  !     The vertical velocity at cloud base is calculated according
  !     to the heat engine framework as defined by Renno'  and
  !     Ingersoll, 1996 Eq(42)
  !

  ! Local variables
  real :: CMI, TCOLD, THOT, FIN, SIGWSHB
  integer :: ICOUNT, K

!  CP=1004.   ! defined in USE shcu_vars_const
  CMI=16.
!  EMIS=.9         ! not used
!  STEFAN=5.67E-8  ! not used
  !
  TCOLD=TE(KLCL)
  !*      TCOLD=TE(2)
  ICOUNT=1
  do K=KLCL+1,KTOP
     !*      DO K=3,KTOP
     TCOLD=TCOLD+TE(K)
     ICOUNT=ICOUNT+1
  enddo
  !
  !     TCOLD is the temperature of the environment, averaged
  !     between the first level and the cloud top.
  !      
  TCOLD=TCOLD/float(ICOUNT)
  THOT=TE(2)
  !
  EFIC=(THOT-TCOLD)/THOT
  !
  if(EFIC.le.0.0.or.CAPE.lt.20.0)then
     IGO=0 
     !        WRITE(11,*)'Eta= ',EFIC,'CAP= ',CAPE       
     return
  endif
  !
  !     the effective vertical velocity at cloud base is calculated 
  !     according to the heat engine framework as deffined by Renno'
  !     and Ingersoll, 1996 Eq.(34)
  !
  FIN=RHOE(2)*(CP*ENTF+ALVL*ALHF)
  !
  if(FIN.le.50.0) then
     IGO=0
     !        WRITE(11,*)'F'
     return
  endif
  !      
  !     TCAPE=2*CAPE
  !
  !     Calculating BCAPE
  !
  !      WSTAR=ENTF*(1+0.608*QVE(1))
  !      WSTAR=(G*ZE(KLCL)*WSTAR/TE(1))**.333
  !      BCAPE=CMI*WSTAR**2
  !
  !     Calculating DCAPE. Here we simply assume DCAPE=0.5CAPE
  !
  !eps DCAPE=1.0*CAPE
  DCAPE=0.
  !
  !      TCAPE=CAPE+DCAPE+BCAPE
  TCAPE=CAPE+DCAPE     
  !
  if(TCAPE.lt.50.0) then
     IGO=0
     return
  endif
  !
  SIGWSHB=EFIC*FIN/(RHOE(KLCL)*TCAPE)
  !      SIGWSHB=SIGWSHB+WPE(KLCL)
  !
  if(SIGWSHB.le.0.0) then
     IGO=0
     !	WRITE(11,*)'G'
     return
  endif
  !
  !
  !     calculating the effetive vertical velocity inside the cloud
  !     interpolating from SIGWSHB in the cloud base to zero
  !     on the cloud top.
  ! 
  !   
  do K=KLCL,KTOP
     WC(K)=SIGWSHB*((ZE(KTOP)-ZE(K))/(ZE(KTOP)-ZE(KLCL)))
  enddo
  !
  return
end subroutine W_SHALLOW
!
!     ******************************************************************
!
subroutine SH_RATES

  use conv_coms, only : NKP,  &   ! INTENT(IN) ! Parameter
       KMT,                   &   ! INTENT(IN)
       KLCL,                  &   ! INTENT(IN)
       ZE,                    &   ! INTENT(IN)
       PKE                        ! INTENT(IN)
  use shcu_vars_const, only : DRDT, &   ! INTENT(OUT)
       DTDT,                        &   ! INTENT(OUT)
       KTOP,                        &   ! INTENT(IN)
       WC,                          &   ! INTENT(IN/OUT)
       DSC,                         &   ! INTENT(IN)
       ALVL,                        &   ! INTENT(IN)  ! Parameter
       WLC,                         &   ! INTENT(IN)
       DSE,                         &   ! INTENT(IN)
       QVC,                         &   ! INTENT(IN)
       QVE,                         &   ! INTENT(IN)
       DQDT,                        &   ! INTENT(OUT) ! Maybe local var.?
       ALHF,                        &   ! INTENT(IN)
       DELZ                             ! INTENT(IN)

  implicit none

! Local variables                                                
  real :: WSSC(NKP), WQSC(NKP), DSDT(NKP)
  real :: WRSUP, WRBASE
  integer :: K, ISUBCL

  do K=1,KMT
     WSSC(K)=0.
     WQSC(K)=0.
     DSDT(K)=0.
     DRDT(K)=0.
     DTDT(K)=0.
  enddo
  !
  !     calculating the transports w's' and w'r'
  !      
  do K=KLCL+1,KTOP
     WSSC(K)=WC(K)*(DSC(K)-ALVL*WLC(K)-DSE(K))
     WQSC(K)=WC(K)*(QVC(K)+WLC(K)-QVE(K))
  enddo
  !     calculating heating and moistening rates due to
  !     shallow non-precipitating cumulus
  !
  do K=KLCL+1,KTOP-1
     DQDT(K)=-(WQSC(K+1)-WQSC(K-1))/(ZE(K+1)-ZE(K-1))
     DRDT(K)=DQDT(K)/(1-QVE(K))**2
     DSDT(K)=-(WSSC(K+1)-WSSC(K-1))/(ZE(K+1)-ZE(K-1))
     DTDT(K)=DSDT(K)/PKE(K)
  enddo
  !
  !     ISUBCL is a flag corresponding to removal (=1) or not (=0)
  !     of moisture from the boudary layer by shallow clouds
  !
  ISUBCL=0
  !                                                           
  if(ISUBCL.eq.0)return
  !
  !     for calculating the removal of moisture from the PBL, we 
  !     consider that the latent heat flux is linearly interpolated between
  !     the surface and the cloud base
  !                            
  WRSUP=ALHF
  WRBASE=WC(KLCL)*(QVC(KLCL)-QVE(KLCL))
  !
  do K=1,KLCL
     WC(K)  =(WRBASE*ZE(K)+(WRSUP*(ZE(KLCL)-ZE(K))))/(QVE(K)*ZE(KLCL))
     WQSC(K)= WC(K)*(QVC(K)-QVE(K))
  enddo
  !
  do K=2,KLCL
     DRDT(K)=-(WQSC(K)-WQSC(K-1))/DELZ(K)
  enddo
  !

!  PRINT *, "Saindo de SH_RATES"

  return
end subroutine SH_RATES
!                           
!*********************************************************************
