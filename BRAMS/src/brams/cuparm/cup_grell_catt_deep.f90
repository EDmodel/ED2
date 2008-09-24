!---------------------------GRELL CUMULUS SCHEME---------------------------
subroutine CUPARTH_CATT( &
     mynum, & 	!01
     mgmxp, &	!02
     mgmyp, &	!03
     mgmzp, &	!04 
     m1   , &	!05
     m2   , &	!06
     m3   , &	!07
     ia   , &	!08
     iz   , &	!09
     ja   , &	!10
     jz   , &	!11
     i0   , &	!12
     j0   , &	!13
     maxiens, &	!14
     iens   , &	!15
     ngrid  , &	!16
     ngrids_cp, & !17
     icbase ,   & !17¼
     depth_min, & !17½
     cap_maxs,  & !17¾
     dtime  , &	!18
     time   , &	!19
     !
     ua     , &	!20
     va     , &	!21
     wa     , &	!22
     theta, &	!23
     pp     , &	!24
     pi0    , &	!25
     dn0    , &	!26
     rv   , &	!27
     kpbl   , & !27½
     tke    , &	!28
     tkmin , &  !29
     rcp    , & !30
     topt   , &	!31
     rtgt   ,  &!32
     !
     tht    , &	!33
     rtt    , &	!34
     !
     pt     , &	!35
     outtem , &	!36
     outrt  , &	!37
     precip , &	!38
     sgrell1_3d, & !39
     sgrell2_3d, & !40
     sgrell3_3d, & !41
     sgrell1_2d, & !42
     ierr4d , &	!43
     jmin4d , &	!44
     kdet4d , &	!45
     k224d  , &	!46
     kbcon4d, &	!47
     ktop4d , &	!48
     kpbl4d , & !49
     kstabi4d, &!50
     kstabm4d, &!51
     xmb4d, &	!52
     edt4d, &	!53
     !
     zcup5d, & !54
     pcup5d, & !55
     enup5d, & !56
     endn5d, & !57
     deup5d, & !58
     dedn5d, & !59
     zup5d, &  !60
     zdn5d, &  !61
     !Lufla	p_lw5d,     &   !57
     prup5d, & !62
     clwup5d, & !63
     tup5d, &   !64	
     !insercao de variaveis para salvar nas analises para uso em
     !cup_direction2
     upmf  , & !65
     dnmf  , & !66
     xierr, &  !67
     xktop, &  !68
     xkbcon, & !69
     xk22, &   !70 Lufla
     xjmin, &  !71
     xkdt , &  !72
     xiact_p, &!73
     xiact_c, &!74
     confrq,  &!75
     frqanl, & !76
     deltaxy2,&!77
     patch_area,&!78
     npat,      &!79
     level)      !80

  !srf   mgmxp, mgmyp, mgmzp sao usadas alocar memoria para as
  !      variaveis da parametrizacao do Grell.

  ! USE Modules for Grell Parameterization
  use mem_grell_param, only : maxens,  & !INTENT(IN)
       maxens2,                        & !INTENT(IN)
       maxens3,                        & !INTENT(IN)
       ensdim,                         & !INTENT(IN)
       icoic                             !INTENT(IN)

  use mem_scratch2_grell, only : &
       massflx,                  &       !INTENT(OUT)
       iact_old_gr,              &       !INTENT(OUT)
       aa0,                      &       !INTENT(OUT)
       xland,                    &       !INTENT(OUT)
       kdt,                      &       !INTENT(OUT)
       iact_gr,                  &       !INTENT(OUT)
       kdet,                     &       !INTENT(OUT)
       pret,                     &       !INTENT(OUT)
       mconv,                    &       !INTENT(OUT)
       umean,                    &       !INTENT(OUT)
       vmean,                    &       !INTENT(OUT)
       pmean,                    &       !INTENT(OUT)
       TER11,                    &       !INTENT(OUT)
       PSUR,                     &       !INTENT(OUT)
       PO,                       &       !INTENT(OUT)
       US_Grell,                 &       !INTENT(OUT)
       VS_Grell,                 &       !INTENT(OUT)
       OMEG,                     &       !INTENT(OUT)
       T,                        &       !INTENT(OUT)
       Q,                        &       !INTENT(OUT)
       TKEG,                     &       !INTENT(OUT)
       RCPG,                     &       !INTENT(OUT)
       TN,                       &       !INTENT(OUT)
       QO,                       &       !INTENT(OUT)
       P,                        &       !INTENT(OUT)
       OUTT,                     &       !INTENT(OUT)
       OUTQ,                     &       !INTENT(OUT)
       OUTQC,                    &       !INTENT(OUT)
       DIRECTION,                &       !INTENT(OUT)
       PBLIDX,                   &       !INTENT(OUT)
       massfln                           !INTENT(?)

![ MLO - Changed to rconstants, to remove Phys_const.o
  use rconstants, only: cp, p00, tcrit, g, cpor,pi1,onerad !INTENT(IN)
!MLO]
  implicit none
  integer, intent(in) :: mgmxp, mgmyp, mgmzp, ngrid, ngrids_cp
  integer, intent(in) :: maxiens ! maxens,maxens2,maxens3,ensdim
  integer, intent(in) :: iens     
  integer, intent(in) :: icbase    ! Method to determine the cloud base; 
  real,    intent(in) :: depth_min ! Minimum depth that the cloud should have [m]
  real,    intent(in) :: cap_maxs  ! Maximum depth of capping inversion [mb]

  integer, intent(in) :: m1, m2, m3, ia, iz, ja, jz, i0, j0, mynum, npat, level

  integer :: j1, j2 !Local

  real(kind=8), intent(in) :: time

  integer, dimension(m2,m3) :: kpbl
  real, dimension(m1,m2,m3), intent(in) :: ua, va, wa, pp, pi0, dn0,   &
       tht, rtt, pt, tke, rcp, theta,rv

  real, dimension(m1,m2,m3), intent(inout) :: outtem, outrt

  real, dimension(m1,m2,m3), intent(inout) :: sgrell1_3d, sgrell2_3d, sgrell3_3d

  real, dimension(m2,m3), intent(inout)    :: sgrell1_2d

  real, dimension(m2,m3,npat), intent(in)  :: patch_area

  real, intent(in) :: tkmin, confrq, frqanl, deltaxy2

  real :: dti, tscl_KF !Local

  real, dimension(m2,m3), intent(in)    :: topt
  real, dimension(m2,m3), intent(inout) :: PRECIP
  real, dimension(m2,m3), intent(in)    :: rtgt !Not used

  real, dimension(m2,m3) :: upmf, dnmf, xierr, xktop, xkbcon, xjmin,  &
                            xkdt, xiact_p, xiact_c,xk22

  !-------salva parametros da CUP para uso no transporte convectivo:
  integer, dimension(m2,m3,maxiens) :: ierr4d, jmin4d,  &
       kdet4d, k224d, kbcon4d, ktop4d, kpbl4d, kstabi4d, kstabm4d

  real,dimension(m2,m3,maxiens) :: xmb4d, edt4d

  real,dimension(m1,m2,m3,maxiens) :: enup5d,  &
                                                         endn5d,  &
							 deup5d,  &
							 dedn5d,  &
							 zup5d,   &
							 zdn5d,   & 
            						 prup5d,  & 
            						 clwup5d, & 
            						 tup5d,   & 
            						 zcup5d,  & 
            						 pcup5d
  !variaveis locais:

  integer :: istart, iend, i, j, k, mix, mjx, mkx, kr, m !kk
  real    :: vspeed, dp, dtime, dq, cpdTdt, exner


  !----------------------------------------------------------------------
  ISTART = ia
  IEND   = iz
  j1     = ja
  j2     = jz
  MKX    = m1 - 1    !MKX nao deve ser igual a m1
  MIX    = m2            
  MJX    = m3            

!- for ensemble average
  dti = confrq/frqanl
  
  do j=1,m3  ! loop em todo dominio para passar informacoes da fronteira
     do i=1,m2 ! dos nodes
        massflx(i,j)     = dnmf(i,j)
        !Lufla iact_gr(i,j)     = INT(xiact_c(i,j))
        iact_old_gr(i,j) = int(xiact_p(i,j))
        !        print*,massflx(i,j),iact_gr(i,j),iact_old_gr(i,j)
     enddo
  enddo

  !srf- out/2004
  ! A. Betts suggestion for the KF timescale < DELTAX  / 50 m/s
    tscl_KF =  sqrt(deltaxy2) / 50.  !units: sec
                                     ! use 25 m/s for shallow
    

  !- loop externo : j
  do J=j1,j2

     do I = ISTART,IEND
        aa0(i)           =0.
        xland(i,j)       = patch_area(i,j,1) ! land < 1 /water=1 flag
        if(xland(i,j) < 0.95) xland(i,j) = 0.
        iact_gr(i,j)     = 0
        kdt(i,j)         = 0
        precip(i,j)      = 0.
     enddo


     !--- Prepare input, erase output

     do I = ISTART,IEND        
        kdet(i)  =2
        pret(i)  =0.
        mconv(i) =0.
        umean(i) =0.
        vmean(i) =0.
        pmean(i) =0.
     enddo

     !------- Transfere valores do RAMS para o eschema
     do K=1,MKX
        kr = K + 1          ! nivel K da grade do Grell corresponde ao
        ! nivel K + 1 do RAMS
        do I = ISTART,IEND

           TER11(I)= topt(i,j)
           ! Pressure in mbar
           PSUR(I) = .5*( ((pp(1,i,j)+pi0(1,i,j))/cp)**cpor*p00 +  &
                          ((pp(2,i,j)+pi0(2,i,j))/cp)**cpor*p00 )*1.e-2
           PO(I,K) = ((pp(kr,i,j)+pi0(kr,i,j))/cp)**cpor*p00*1.e-2

           US_Grell(I,K) = .5*( ua(kr,i,j) + ua(kr,i-1,j) )
           VS_Grell(I,K) = .5*( va(kr,i,j) + va(kr,i,j-1) )
           OMEG(I,K)   = -g*dn0(kr,i,j)*.5*( wa(kr,i,j)+wa(kr-1,i,j) )

           T(I,K)  = theta(kr,i,j)*(pp(kr,i,j)+pi0(kr,i,j))/cp
           Q(I,K)  = rv(kr,i,j)

           !- variables for PBL top height
           pblidx(i) = kpbl(i,j)
	   TKEG(I,K) = TKE(kr,i,j)
	   RCPG(I,K) = RCP(kr,i,j)
           !        Calcula tendencia projetada na temperatura em funcao 
           !        das tendencias de theta e PI : cp*T=Pi*Theta
           exner= pp(kr,i,j)+pi0(kr,i,j)

           !        cpdTdt= exner*tht(kr,i,j) + theta(kr,i,j)*pt(kr,i,j)
           cpdTdt  = exner*tht(kr,i,j)
           ! assumindo PT(KR,I,J) << exner*THT(KR,I,J)/theta

           !        Temperatura projetada se a conveccao nao ocorrer
           TN(I,K) = T(I,K) + ( cpdTdt/cp )*dtime

           !        Umidade projetada se a conveccao nao ocorrer
           QO(I,K) = Q(I,K) +   rtt(kr,i,j)*dtime

           !-------

           P(I,K)  = PO(I,K)

           if((PSUR(I)-P(I,K)).gt.150.and.P(I,K).gt.300.)then   
              DP       = -.5*(P(I,K+1)-P(I,K-1))        
              UMEAN(I) = UMEAN(I)+US_Grell(I,K)*DP 
              VMEAN(I) = VMEAN(I)+VS_Grell(I,K)*DP
              PMEAN(I) = PMEAN(I)+DP
           endif

           if(TN(I,K).lt.200.)    TN(I,K) = T(I,K)
           if(QO(I,K).lt.1.E-08)  QO(I,K) = 1.E-08

           OUTT(I,K)  = 0. !- Tendencia no campo de temperatura 
	                   !  associada aos cumulus           
           OUTQ(I,K)  = 0. !- Tendencia na razao de mist. de vapor d'agua
	                   !  assoc. aos cumulus
           OUTQC(I,K) = 0. !- Tendencia na razao de mistura de agua de nuvem e/ou gelo
                           !  associada aos cumulus
        enddo
     enddo

     do I = ISTART,IEND
        UMEAN(I)=UMEAN(I)/PMEAN(I)
        VMEAN(I)=VMEAN(I)/PMEAN(I)
        VSPEED=sqrt(UMEAN(I)*UMEAN(I)+VMEAN(I)*VMEAN(I))
        DIRECTION(I)=(atan2(UMEAN(I),VMEAN(I))+pi1)*onerad
        if(DIRECTION(I).gt.360.)DIRECTION(I)=DIRECTION(I)-360.
        if(VSPEED.lt.5.)DIRECTION(I)=9999.
        !sgrell1_3d(24,i,j)=UMEAN(I)
        !sgrell1_3d(25,i,j)=VMEAN(I)
     enddo

     do K=2,MKX-1
        do I = ISTART,IEND
           dq=.5*(q(i,k+1)-q(i,k-1))
           !- convergencia de umidade da coluna (omega em pa/s)
           mconv(i)=mconv(i)+omeg(i,k)*dq/g 
        enddo
     enddo
     do I = ISTART,IEND
        if(mconv(i) .lt. 0.)  mconv(i) = 0.
     enddo

     !---  CUMULUS PARAMETERIZATION
     !srf- aqui se deve colocal o loop no ensemble dependente do tipo de cumulus
     ! iens =1

     call CUP_enss_catt(ngrid,mynum,m1,m2,m3,i0,j0,                     &
          mgmxp,mgmyp,mgmzp,maxiens,maxens,maxens2,maxens3,             &
	  ensdim,icoic, icbase, depth_min, cap_maxs,j,iens,istart,iend,mix,mjx,mkx,                  &
	  massfln,massflx,iact_gr,iact_old_gr,xland,ter11,              &
	  aa0, t, q, tn, qo, po, pret, p, outt, outq, outqc,            &
          dtime, psur, us_grell, vs_grell,kdet,                         &
	  tcrit,time,mconv, omeg, direction, pblidx, tkeg,rcpg,tkmin,   &
            ierr4d(1:m2,j,iens),   jmin4d(1:m2,j,iens),           &
            kdet4d(1:m2,j,iens),    k224d(1:m2,j,iens),           &
           kbcon4d(1:m2,j,iens),   ktop4d(1:m2,j,iens),           &
            kpbl4d(1:m2,j,iens),                               &
          kstabi4d(1:m2,j,iens), kstabm4d(1:m2,j,iens),           &
             xmb4d(1:m2,j,iens),    edt4d(1:m2,j,iens),           &
          zcup5d(1:m1,1:m2,j,iens), pcup5d(1:m1,1:m2,j,iens),           &
          enup5d(1:m1,1:m2,j,iens), endn5d(1:m1,1:m2,j,iens),           &
          deup5d(1:m1,1:m2,j,iens), dedn5d(1:m1,1:m2,j,iens),           &
           zup5d(1:m1,1:m2,j,iens),  zdn5d(1:m1,1:m2,j,iens),           &
          prup5d(1:m1,1:m2,j,iens),clwup5d(1:m1,1:m2,j,iens),           &
           tup5d(1:m1,1:m2,j,iens),                                     &
          upmf,dnmf,xierr,xktop,xkbcon,xk22,xjmin,xkdt,xiact_p,xiact_c, &
          sgrell1_3d,sgrell2_3d,sgrell1_2d,dti,tscl_KF)

     !--- Output 
     !srf out/2004: coupling with RAMS microphysics 
!    if(level == 1) then
     if(level > 0) then
      do K=1,MKX-1 
        kr = K + 1
         do I = ISTART,IEND
           ! Converte tend da temperatura (OUTT) em tend de theta (OUTTEM)
           ! cp*T=Pi*Theta => cp dT/dt = Theta*dPi/dt + Pi*dTheta/dt,
           ! assumindo dPi/dt (=pt(kr,i,j)) << (exner/theta)*dTheta/dt:
           ! Exner's function = pp(kr,i,j)+pi0(kr,i,j)
           exner          = pp(kr,i,j) + pi0(kr,i,j)
           ! tendencia do Theta   devida a conv profunda
           outtem(kr,i,j) = CP/exner   * OUTT(I,K) 
           ! tendencia do Rtotal  devida a conv profunda
           outrt(kr,i,j)  = OUTQ(I,K)  + OUTQC(I,K) 
         enddo
       enddo
!    elseif(level > 1) then 
     elseif(level > 1111111) then
      do K=1,MKX-1 
        kr = K + 1
         do I = ISTART,IEND
           ! Converte tend da temperatura (OUTT) em tend de theta (OUTTEM)
           ! cp*T=Pi*Theta => cp dT/dt = Theta*dPi/dt + Pi*dTheta/dt,
           ! assumindo dPi/dt (=pt(kr,i,j)) << (exner/theta)*dTheta/dt:
           ! Exner's function = pp(kr,i,j)+pi0(kr,i,j)
           exner          = pp(kr,i,j) + pi0(kr,i,j)
           ! tendencia do Theta  devida a conv profunda
           outtem    (kr,i,j) = CP/exner * OUTT(I,K) 
           ! tendencia da vapor d'agua devida a conv profunda
           outrt     (kr,i,j) = OUTQ(I,K)  
           ! tendencia da agua condensada devida a conv profunda
	   sgrell3_3d(kr,i,j)= OUTQC(I,K) 
         enddo
       enddo
      endif

     do I = ISTART,IEND
        PRECIP(I,J)=PRET(I)

     enddo



     do I = ISTART,IEND
        if(precip(i,j).le.0.)then
           iact_gr(i,j) =0
           precip(i,j)  =0.
           do k=1,mkx
              kr = k + 1
              outtem(kr,i,j)    = 0.! comente p/ testes com conservacao (c0=0.)              
              outrt(kr,i,j)     = 0.! comente p/ testes com conservacao (c0=0.)              
	      sgrell3_3d(kr,i,j)= 0.! comente p/ testes com conservacao (c0=0.) 
           enddo
           do k=1,ensdim
              massfln(i,j,k)=0.
           enddo
        else
           iact_gr(i,j)=1
        endif
     enddo

     !--- Salva nas analises

     do I = ISTART,IEND
       !massflx(i,j) = dnmf(i,j)
        xiact_c(i,j) = float(IACT_GR(I,J))
       !xiact_p(i,j) = float(iact_old_gr(I,J))
     enddo

  enddo     ! loop externo - j -


!--- shallow convection

end subroutine CUPARTH_CATT

!--------------------------------------------------------------------
!
!

subroutine CUP_enss_catt(ngrid, mynum, m1, m2, m3, i0, j0,                      &
     	       mgmxp, mgmyp, mgmzp, maxiens, maxens, maxens2, maxens3, ensdim,  &
     	       icoic, icbase,depth_min,cap_maxs, j, iens, ISTART, IEND, mix, mjx, mkx, massfln, massflx,	&
     	       iact_gr, iact_old_gr, xland, Z1,                          &
     	       AAEQ, T, Q, TN, QO, PO, PRE, P, OUTT, OUTQ, OUTQC, DTIME, &
     	       PSUR, US, VS, KDET, TCRIT, time, mconv, omeg, direction,  &
	        pblidx,tkeg,rcpg,tkmin,  & 
     	       ierr4d, jmin4d, kdet4d, k224d, kbcon4d,ktop4d, kpbl4d,    &
     	       kstabi4d, kstabm4d, xmb4d, edt4d,                         &
	       zcup5d,pcup5d, enup5d, endn5d, deup5d, dedn5d,  &
     	       zup5d, zdn5d, prup5d,clwup5d,tup5d,             & 
     	       upmf, dnmf, xierr, xktop, xkbcon, xk22, xjmin,xkdt, xiact_p, xiact_c, &
               sgrell1_3d,sgrell2_3d,sgrell1_2d,dti,tscl_KF)

!- USE Modules for Grell Cumulus Parameterization
  use mem_scratch3_grell
![MLO - switching to rconstants
  use rconstants, only: rgas,cp,rm,p00,g,cpor,pkdcut,day_sec,alvl
!MLO]
  implicit none
  integer maxiens,maxens,maxens2,maxens3,ensdim
  integer mix,mjx,mkx,mgmxp, mgmyp, mgmzp
  integer nall,nens,iens,iedt    !,ktau
  integer nens3
  integer icoic
  real(kind=8) :: time
  integer,intent(in) :: icbase
  real, intent(in) :: depth_min
  real, intent(in) :: cap_maxs

  !--- Input variables -----------------------------
  !
  ! basic environmental input includes moisture convergence (mconv)
  ! omega (omeg), windspeed (us,vs), and a flag (aaeq) to turn off
  ! convection for this call only and at that particular gridpoint
  !
  integer pblidx(mgmxp)
  real mconv(mgmxp), Z1(mgmxp), direction(mgmxp), AAEQ(mgmxp),  &
       pre(mgmxp), PSUR(mgmxp)
  real T(mgmxp,mgmzp), Q(mgmxp,mgmzp), TN(mgmxp,mgmzp), QO(mgmxp,mgmzp),   &
       P(mgmxp,mgmzp), PO(mgmxp,mgmzp), US(mgmxp,mgmzp), VS(mgmxp,mgmzp),  &
       omeg(mgmxp,mgmzp), tkeg(mgmxp,mgmzp),rcpg(mgmxp,mgmzp)

  real massfln(mgmxp,mgmyp,ensdim)
  real massflx(mgmxp,mgmyp), xland(mgmxp,mgmyp)
  integer iact_gr(mgmxp,mgmyp), iact_old_gr(mgmxp,mgmyp), kdet(mgmxp)
  integer ngrid, mynum, i0, j0, m1, m2, m3

  !-------CU_P parameters for the convective transport:
  integer, dimension(m2) :: ierr4d, jmin4d, kdet4d, k224d, kbcon4d,  &
       ktop4d, kpbl4d, kstabi4d, kstabm4d

  real, dimension(m2) :: xmb4d, edt4d

  real, dimension(m1,m2) ::            &
                    zcup5d, pcup5d,          &
                    enup5d, endn5d, deup5d,  &
                    dedn5d,  zup5d,  zdn5d,  & 
                    prup5d,clwup5d,   tup5d 

  !------Variables saved in RAMS Analisys
  ! use (m1,m2,m3) para dimensionar os vetores que sao
  ! escritos nas analises do RAMS
  real, dimension(m2,m3) :: upmf, dnmf, xierr, xktop, xkbcon, xjmin,    &
       xkdt, xiact_p, xiact_c, xk22 
  real, dimension(m1,m2,m3) :: sgrell1_3d,sgrell2_3d
  real, dimension(m2,m3)    :: sgrell1_2d
  real tkmin,dti,tscl_KF
  !
  !--- Work variables - Allocatable in this point --
  !
  integer fquasi, fstab, fmconv, iresult
  integer ki, m       
  integer I, J, K, ISTART, IEND 

  real mbdt

  !--- Output variables ----------------------------
  ! outt   = output temp tendency (per s)
  ! outq   = output q tendency (per s)
  ! outqc  = output qc tendency (per s)
  ! pre    = output precip

  real OUTQC(mgmxp,mgmzp), OUTT(mgmxp,mgmzp), OUTQ(mgmxp,mgmzp)

  real dz, tcrit, dtime  
  real dellaqsum, dellaqcsum, dp

  !--- New entrainment/detrainment related stuff --------------------

  real mentr_rate, mentrd_rate, entr_rate, radius,              &  !entrd_rate,  &
       massfld, zcutdown, edtmax, edtmin, zkbmax,    &
       z_detr, zktop, dh

  real :: fdiur, tke_start
  integer :: kstart,kstart_way 
  integer,parameter :: i_cup_dir_flag = 0
  
  !---   the following are your basic environmental
  !	 variables. They carry a "_cup" if they are
  !	 on model cloud levels (staggered). They carry
  !	 an "o"-ending (z becomes zo), if they are the forced
  !	 variables. They are preceded by x (z becomes xz)
  !	 to indicate modification by some typ of cloud
  !
  ! z  	        = heights of model levels
  ! q  	        = environmental mixing ratio
  ! qes         = environmental saturation mixing ratio
  ! t           = environmental temp
  ! p           = environmental pressure
  ! he          = environmental moist static energy
  ! hes         = environmental saturation moist static energy
  ! z_cup       = heights of model cloud levels
  ! q_cup       = environmental q on model cloud levels
  ! qes_cup     = saturation q on model cloud levels
  ! t_cup       = temperature (Kelvin) on model cloud levels
  ! p_cup       = environmental pressure
  ! he_cup      = moist static energy on model cloud levels
  ! hes_cup     = saturation moist static energy on model cloud levels
  ! gamma_cup   = gamma on model cloud levels
  ! hcd         = moist static energy in downdraft
  ! zd          = normalized downdraft mass flux
  ! dby         = buoancy term
  ! entr        = entrainment rate
  ! zd          = downdraft normalized mass flux
  ! entr        = entrainment rate
  ! hcd         = h in model cloud
  ! bu          = buoancy term
  ! zd          = normalized downdraft mass flux
  ! gamma_cup   = gamma on model cloud levels
  ! mentr_rate  = entrainment rate
  ! qcd         = cloud q (including liquid water) after entrainment
  ! qrch        = saturation q in cloud
  ! pwd         = evaporate at that level
  ! pwev        = total normalized integrated evaoprate (I2)
  ! entr        = entrainment rate
  ! z1          = terrain elevation
  ! entr        = downdraft entrainment rate
  ! jmin        = downdraft originating level
  ! kdet        = level above ground where downdraft start detraining
  ! psur        = surface pressure
  ! z1          = terrain elevation
  ! pr_ens      = precipitation ensemble
  ! xf          = mass flux ensembles
  ! massfln     = downdraft mass flux ensembles used in next timestep
  ! omeg        = omega from large scale model
  ! mconv       = moisture convergence from large scale model
  ! zd          = downdraft normalized mass flux
  ! zu          = updraft normalized mass flux
  ! dir         = "storm motion"
  ! mbdt        = arbitrary numerical parameter
  ! dtime       = dt over which forcing is applied
  ! iact_gr_old = flag to tell where convection was active
  ! kbcon       = LFC of parcel from k22
  ! k22         = updraft originating level
  ! icoic       = flag if only want one closure (usually set to zero!)
  ! dby         = buoancy term
  ! ktop        = cloud top (output)
  ! xmb         = total base mass flux
  ! hc          = cloud moist static energy
  ! hkb         = moist static energy at originating level
  ! mentr_rate  = entrainment rate
  ! cd          = detrainment function for updraft
  ! cdd         = detrainment function for downdraft
  ! dellat      = change of temperature per unit mass flux of cloud ensemble
  ! dellaq      = change of q per unit mass flux of cloud ensemble
  ! dellaqc     = change of qc per unit mass flux of cloud ensemble
  ! aa0         = cloud work function for downdraft
  ! edt         = epsilon
  ! aa0         = cloud work function without forcing effects
  ! aa1         = cloud work function with forcing effects
  ! xaa0        = cloud work function with cloud effects (ensemble dependent)
  ! edt         = epsilon

  !----------------------------------------------------------------
  !
  !     if(ktau.gt.3.and.ktau.lt.7) ...
  !
  !
  !--- specify entrainmentrate and detrainmentrate
  !
  !     radius=14000.-float(iens)*2000.
  radius=12000.
  !      radius=5000.
  fquasi=1
  fstab=0  
  fmconv=0
  !
  !--- gross entrainment rate (these may be changed later on in the
  !--- program, depending what your detrainment is!!)
  !
  entr_rate=.2/radius
  !
  !--- entrainment of mass
  !
  !      mentrd_rate=0.
  mentrd_rate=entr_rate
  mentr_rate =entr_rate
  !
  !--- initial detrainmentrates
  !
  do k=1,mkx
     do i=istart,iend
        cd(i,k)  = 0.1*entr_rate
        !if(iens.gt.4)then
        !   cd(i,k)=(0.3+float(iens-4)*.15)*entr_rate
        !endif
        cdd(i,k) = 0.
     enddo
  enddo
  !
  !--- max/min allowed value for epsilon (ratio downdraft base mass flux/updraft)
  !    base mass flux
  !
  edtmax=.95
  edtmin=.2
  !
  cap_max_increment(:)=20.

  !--- initialize cap_max 
  do i=istart,iend          
     cap_max(i)=cap_maxs    
  end do                    

  do I=ISTART,IEND
     aa0(i)=0.
     aa1(i)=0.
     aad(i)=0.
     kstabm(i)=mkx-2
     if (aaeq(i).lt.0.) then
        ierr(i)=20
     else
        IERR(i)=0
        XIERR(i,j)=0.
        cwf(i,j)=0.
        pwf(i,j)=0.
        pwdf(i,j)=0.
        eddt(i,j)=0.
        ! xktop(i,j)=0.
        ! xkbas(i,j)=0.
        xmass(i,j)=0.
        predb(i,j)=0.
     endif
       ierr2(i)=ierr(i)
       ierr3(i)=ierr(i)
  enddo


  !--- first check for upstream convection
  if(i_cup_dir_flag == 1 ) then

   do i=istart,iend
     if (ierr(i).eq.0) then
        iresult=0
        massfld=0.
        call cup_direction2(i,j,direction,iact_old_gr,mix,mjx,  &
             mgmxp,mgmyp,massflx,iresult,ensdim,0,0,maxens3,massfld)

  !--- increase cap_max if is there upstream convection
        if (iresult.eq.1) then
           cap_max(i)=cap_max(i)+20.
        endif
     endif
   enddo

  endif

  !--- max height(m) above ground where updraft air can originate

  zkbmax=4000.

  !--- height(m) above which no downdrafts are allowed to originate

  zcutdown=3000.

  !--- depth(m) over which downdraft detrains all its mass

  z_detr=1250.

  !--- MBDT parameter 
  !srf - for capmax ensemble mbdt_ens is constant (=4e-3*dtime)
  !      independent of nens
  do nens=1,maxens
     mbdt_ens(nens)=(float(2)-3.)*dtime*1.e-3+dtime*5.E-03
    !mbdt_ens(nens)=(float(nens)-1.5)*dtime*2.e-3+dtime*5.E-03
  enddo
  do nens=1,maxens2
    !edt_ens(nens)=.7-float(nens)*.1
     edt_ens(nens)=.95-float(nens)*.01
  enddo
  !     if(j.eq.jpr)then
  !       print *,'radius ensemble ',iens,radius
  !       print *,mbdt_ens
  !       print *,edt_ens
  !     endif
  !
  !--- environmental conditions, FIRST HEIGHTS
  !
  do i=istart,iend
     if(ierr(i).ne.20)then
        do k=1,maxens*maxens2*maxens3
           xf_ens(  i,j,(iens-1)*maxens*maxens2*maxens3+k)= 0.
           pr_ens(  i,j,(iens-1)*maxens*maxens2*maxens3+k)= 0.
           outt_ens(i,j,(iens-1)*maxens*maxens2*maxens3+k)= 0.
        enddo
     endif
  enddo

  !--- calculate moist static energy, heights, qes
  !
  call cup_env(j, z, qes, he, hes, t, q, p, z1, mix, mgmxp,  &
       mkx, mgmzp, istart, iend, psur, ierr, tcrit, 0)
  call cup_env(j, zo, qeso, heo, heso, tn, qo, po, z1, mix,mgmxp, &
       mkx, mgmzp, istart, iend, psur, ierr, tcrit, 0)
  !
  !--- environmental values on cloud levels
  !
  call cup_env_clev(j, t, qes, q, he, hes, z, p, qes_cup,  &
       q_cup, he_cup, hes_cup, z_cup, p_cup, gamma_cup, t_cup, psur, &
       mix, mgmxp, mkx, mgmzp, istart, iend, ierr, z1)
  call cup_env_clev(j, tn, qeso, qo, heo, heso, zo, po,    &
       qeso_cup, qo_cup, heo_cup, heso_cup, zo_cup, po_cup,          &
       gammao_cup, tn_cup, psur, mix, mgmxp, mkx, mgmzp, istart,     &
       iend, ierr, z1)

  do i=istart,iend
     if (ierr(i).eq.0) then

        do k=1,mkx
           if (zo_cup(i,k).gt.zkbmax+z1(i)) then
              kbmax(i)=k
              !GO TO 25
              exit
           endif
        enddo
        !25      CONTINUE

        !--- level where detrainment for downdraft starts

        do k=1,mkx
           if (zo_cup(i,k) .gt. z_detr+z1(i)) then
              kdet(i)=k
              !GO TO 26
              exit
           endif
        enddo
        !26      CONTINUE

     endif
  enddo

  !--- DETERMINE LEVEL WITH HIGHEST MOIST STATIC ENERGY CONTENT - K22

  !   here, ,should just take top of PBL for shallow clouds, also maybe for
  !   deep clouds in tropics: kstart=level(pbltop)
  !srf-fev2003
  kstart_way = 0

!srf- TESTAR NO FUTURO - 
  if (icbase == 2 .and. pblidx(i) == 0) then
     !	  New way to define K22
     !----   Determine PBL top using TKE (TKEG) and liquid water mixing ratio (RCPG)
     call get_zi(mix,mgmxp,mkx,mgmzp,istart,iend,j,ierr,kzi,&
          TKEG,RCPG,zo,z1,tkmin)

     do i=istart,iend
        !srf-14-fev-2003
        !A segunda forma produz uma cobertura de shallow mais representativa
        !Em caso de alta resolucao vertical tente a versao original (k22=kzi)
        !	     IF(ierr(i).eq.0) k22(i) = kzi(i)
        if(ierr(i) == 0) k22(i) = max(2, kzi(i) - 1)
        !	    IF(ierr(i).eq.0) PRINT*,k22(i)
     enddo
  elseif (icbase == 2 .and. pblidx(i) /= 0) then
     k22(i) = pblidx(i)
  else !(icbase == 1)
     !	 Old way to define k22
     kzi(1:mgmxp) = 1
     kstart = 3
     call maximi(heo_cup,mix,mgmxp,mkx,mgmzp,kstart,kbmax,k22,istart,iend,ierr)
  endif
  do I=ISTART,IEND
     if (ierr(I).eq.0.) then
        if (K22(I).ge.KBMAX(i)) ierr(i)=2
     endif
  enddo


  !--- DETERMINE THE LEVEL OF CONVECTIVE CLOUD BASE  - KBCON
  !--- srf: out2004
!srf- TESTAR NO FUTURO - 
  !--- new version that includes the PBL mean for hkbo
  !--- not working yet!
    call cup_kbcon_catt(cap_max_increment,1,k22,kbcon,heo_cup,heso_cup, &
         hkbo,kzi,mix,mgmxp,mkx,mgmzp,istart,iend,ierr,kbmax,po_cup,cap_max, &
	 j)
  !- srf: out2004
   do i=istart,iend
    if(ierr(i).eq.0) hkb(i)=hkbo(i)
   enddo

  !srf 30-jan-2002. This routine blow up the model
  !CALL cup_kbcon_cin(1, k22,kbcon, heo_cup, heso_cup, z, tn_cup,  &
  !     qeso_cup, mix, mgmxp, mkx, mgmzp, istart, iend, ierr,      &
  !     kbmax, po_cup, cap_max)

  !--- Increase detrainment in stable layers

  call MINIMI(HEso_cup, mix, mgmxp, mkx, mgmzp, Kbcon, kstabm,  &
       kstabi, ISTART, IEND, ierr)


  do I=ISTART,IEND
     if (ierr(I).eq.0.) then
        if (kstabm(i)-1.gt.kstabi(i)) then
           do k=kstabi(i),kstabm(i)-1
	      cd(i,k)=cd(i,k-1)+1.5*entr_rate
              if (iens.gt.4) then
                 cd(i,k) = cd(i,k-1)+float(iens-4)*entr_rate/  &
                           float(kstabm(i)-kstabi(i))
              else
                 cd(i,k)=cd(i,k)
              endif
              if (cd(i,k).gt.10.0*entr_rate) cd(i,k)=10.0*entr_rate
           enddo
        endif
     endif
  enddo

  !--- Calculate incloud moist static energy

  call cup_up_he(k22, hkb, z_cup, cd, mentr_rate, he_cup, hc,   &
       mix, mgmxp, mkx, mgmzp, kbcon, ierr, istart, iend, dby,  &
       he, hes_cup)
  call cup_up_he(k22, hkbo, zo_cup, cd, mentr_rate, heo_cup, hco,  &
       mix, mgmxp, mkx, mgmzp, kbcon, ierr, istart, iend, dbyo,    &
       heo, heso_cup)

  !srf
  call cup_ktop(1, dbyo, kbcon, ktop, mix, mgmxp, mkx, mgmzp, istart, &
       iend, ierr)

  do I=ISTART,IEND
     kzdown(i)=0
     if (ierr(i).eq.0) then
        zktop = (zo_cup(i,ktop(i))-z1(i))*.6
        zktop = min(zktop+z1(i),zcutdown+z1(i))
        do k=1,mkx
           if (zo_cup(i,k).gt.zktop) then
              kzdown(i)=k
              !GO TO 37
              exit
           endif
        enddo
     endif
   enddo


  !--- Downdraft originating level - JMIN
  call MINIMI(HEso_cup, mix, mgmxp, mkx, mgmzp, K22, kzdown,  &
       JMIN, ISTART, IEND, ierr)

  do I=ISTART,IEND
     !  DO 100 I=ISTART,IEND
     if (ierr(I).eq.0.) then

        !--- Check whether it would have buoyancy, if there where
        !--- no entrainment/detrainment

101     continue
        if (jmin(i)-1 .lt. kdet(i)  ) kdet(i) = jmin(i)-1
        if (jmin(i)   .ge. ktop(i)-1) jmin(i) = ktop(i)-2
        ki=jmin(i)
        hcdo(i,ki) = heso_cup(i,ki)
        DZ         = Zo_cup(i,Ki+1)-Zo_cup(i,Ki)
        dh         = dz*(HCDo(i,Ki)-heso_cup(i,ki))
        dh         = 0.

        do k=ki-1,1,-1
           hcdo(i,k) = heso_cup(i,jmin(i))
           DZ        = Zo_cup(i,K+1)-Zo_cup(i,K)
           dh        = dh+dz*(HCDo(i,K)-heso_cup(i,k))
           if (dh .gt. 0.) then
              jmin(i)=jmin(i)-1
              if (jmin(i) .gt. 3     ) then
                 GO TO 101
              else if (jmin(i) .le. 3) then
                 ierr(i)=9
                 GO TO 100
              endif
           endif
        enddo

        if (JMIN(I).le.3) then
           ierr(i)=4
        endif

     endif
100  continue
  enddo


  !--- Must have at least depth_min m between cloud convective
  !    base and cloud top.

  do i=istart,iend
     if (ierr(i).eq.0.) then
        if (-zo_cup(i,kbcon(i))+zo_cup(i,ktop(i)) .lt. depth_min) ierr(i)=6
     endif
  enddo

  !--- Normalized updraft mass flux profile

  call cup_up_nms(zu, z_cup, mentr_rate, cd, kbcon, ktop,  &
       mix, mgmxp, mkx, mgmzp, istart, iend, ierr, k22)
  call cup_up_nms(zuo, zo_cup, mentr_rate, cd, kbcon, ktop,  &
       mix, mgmxp, mkx, mgmzp, istart, iend, ierr, k22)

  !--- Normalized downdraft mass flux profile,also work on
  !    bottom detrainment
  !--- in this routine

  call cup_dd_nms(zd, z_cup, cdd, mentrd_rate, jmin, ierr,    &
       mix, mgmxp, mkx, mgmzp, istart, iend, 0, kdet, z1)
  call cup_dd_nms(zdo, zo_cup, cdd, mentrd_rate, jmin, ierr,  &
       mix, mgmxp, mkx, mgmzp, istart, iend, 1, kdet, z1)

  !--- Downdraft moist static energy

  call cup_dd_he(hes_cup, zd, hcd, z_cup, cdd, mentrd_rate,  &
       jmin, ierr, mix, mgmxp, mkx, mgmzp, istart, iend,he,  &
       kdet, dbyd, he_cup)
  call cup_dd_he(heso_cup, zdo, hcdo, zo_cup, cdd,           &
       mentrd_rate, jmin, ierr, mix, mgmxp, mkx, mgmzp,      &
       istart, iend,heo, kdet, dbydo, he_cup)

  !--- Calculate moisture properties of downdraft

  call cup_dd_moisture(j, zd, hcd, hes_cup, qcd, qes_cup,  &
       pwd, q_cup, z_cup, cdd, mentrd_rate, jmin, ierr,    &
       gamma_cup, pwev, mix, mgmxp, mkx, mgmzp, istart,    &
       iend, bu, qrcd, q, he, hc, t_cup, 2)
  call cup_dd_moisture(j, zdo, hcdo, heso_cup, qcdo,       &
       qeso_cup, pwdo, qo_cup, zo_cup, cdd, mentrd_rate,   &
       jmin, ierr, gammao_cup, pwevo, mix, mgmxp, mkx,     &
       mgmzp, istart, iend, bu, qrcdo, qo, heo, hco, tn_cup, 1)

  !--- Calculate moisture properties of updraft

  call cup_up_moisture(ierr, z_cup, qc, qrc, pw, pwav,     &
       kbcon, ktop, mix, mgmxp, mkx, mgmzp, istart, iend,  &
       cd, dby, mentr_rate, q, GAMMA_cup, zu, qes_cup,     &
       k22, q_cup)
  call cup_up_moisture(ierr, zo_cup, qco, qrco, pwo, pwavo,  &
       kbcon, ktop, mix, mgmxp, mkx, mgmzp, istart, iend,    &
       cd, dbyo,                                             &
       mentr_rate, qo, GAMMAo_cup, zuo, qeso_cup, k22, qo_cup)

  !--- Calculate workfunctions for updrafts

  call cup_up_aa0(aa0, z, zu, dby, GAMMA_CUP, t_cup, kbcon,  &
       ktop, mix, mgmxp, mkx, mgmzp, istart, iend, ierr)
  call cup_up_aa0(aa1, zo, zuo, dbyo, GAMMAo_CUP, tn_cup,    &
       kbcon, ktop, mix, mgmxp, mkx, mgmzp, istart, iend, ierr)
  do i=istart,iend
     if (ierr(i).eq.0) then
        if (aa1(i).eq.0.) ierr(i)=17
     endif
  enddo

  !--- Determine downdraft strength in terms of windshear

  call cup_dd_edt(ierr, us, vs, zo, ktop, kbcon, edt, po,   &
       pwavo, pwevo, mix, mgmxp, mkx, mgmzp, istart, iend,  &
       edtmax, edtmin, maxens2, edtc, vshear, sdp, vws)

!- srf - big loop starts here! ---------------------------------------------

  do iedt=1,maxens2 !orig: DO 250 iedt=1,maxens2


!-----print
!     if(j==21) then
!     print*,'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
!     print*,' I E D T =', iedt
!     endif
!-----print

     do i=istart,iend
        if (ierr(i).eq.0) then
           edt(i)  = edtc(i,iedt)
           edto(i) = edtc(i,iedt)
           edtx(i) = edtc(i,iedt)
        endif
     enddo
     do k=1,mkx
        do i=istart,iend
            dellat_ens(i,k,iedt) = 0.
            dellaq_ens(i,k,iedt) = 0.
           dellaqc_ens(i,k,iedt) = 0.
               pwo_ens(i,k,iedt) = 0.
        enddo
     enddo

     !--- Change per unit mass that a model cloud would modify the environment

     !--- 1. in bottom layer

     call cup_dellabot_catt(heo_cup, ierr, zo_cup, po_cup,                &
                            hcdo, edto, zdo, cdd, heo, mix, mgmxp, mkx, mgmzp,  &
                            istart, iend, dellah, 1, j, mentrd_rate, zo)
     call cup_dellabot_catt(qo_cup, ierr, zo_cup,                 &
                            po_cup, qrcdo, edto, zdo, cdd, qo, mix, mgmxp,  &
                            mkx, mgmzp, istart, iend, dellaq, 2, j,         &
                            mentrd_rate, zo)

     !--- 2. everywhere else

     call cup_dellas_catt(ierr, zo_cup, po_cup, hcdo, edto, zdo,             &
                          cdd, heo, mix, mgmxp, mkx, mgmzp, istart, iend,    &
                          dellah, 1, j, mentrd_rate, zuo, cd, hco, ktop,     &
                          k22, kbcon, mentr_rate, jmin, heo_cup, kdet, k22,  &
                          'deep')

     !-- Take out cloud liquid water for detrainment

     do k=1,mkx
        do i=istart,iend
           scr1(i,k)   = 0.
           dellaqc(i,k)= 0.
           
	   if (ierr(i).eq.0) then
           
	      scr1(i,k)=qco(i,k)-qrco(i,k)

              if (k.eq.ktop(i)-0) dellaqc(i,k) = .01*zuo(i,ktop(i))*  &
                                  qrco(i,ktop(i))*g/(po_cup(i,k)-po_cup(i,k+1))

              if (k.lt.ktop(i).and.k.gt.kbcon(i)) then
                 
		 dz = zo_cup(i,k+1)-zo_cup(i,k)
                 dellaqc(i,k) = .01*g*cd(i,k)*dz*zuo(i,k)*            &
                                .5*(  qrco(i,k) +   qrco(i,k+1) )/    &
                                   (po_cup(i,k) - po_cup(i,k+1) )
              endif
           endif
        enddo
     enddo

     call cup_dellas_catt(ierr, zo_cup, po_cup, qrcdo, edto, zdo,  &
          cdd, qo, mix, mgmxp, mkx, mgmzp, istart, iend,      &
          dellaq, 2, j, mentrd_rate, zuo, cd, scr1, ktop,     &
          k22, kbcon, mentr_rate, jmin, qo_cup, kdet, k22,    &
          'deep')

     !--- Using dellas, calculate changed environmental profiles

        mbdt=mbdt_ens(1)
        do i=istart,iend
           do k=1,maxens
              xaa0_ens(i,k)=0.
          enddo
        enddo

        do k=1,mkx-1
           do i=istart,iend
              dellat(i,k)=0.
              if (ierr(i).eq.0) then
                 XHE     (I,K) = DELLAH(I,K)*MBDT + HEO(I,K)
                 XQ      (I,K) = DELLAQ(I,K)*MBDT +  QO(I,K)
                 DELLAT  (I,K) = (1./cp)*(DELLAH(I,K)-alvl*DELLAQ(I,K))
                 XT_Grell(I,K) = DELLAT(I,K)*MBDT +  TN(I,K) !(XT_Grell == old XT)
                 if (XQ(I,K).le.0.) XQ(I,K)=1.E-08
              endif
           enddo
        enddo
        !
        do i=istart,iend
           if (ierr(i).eq.0) then
              XHE     (I,mkx)  = HEO(I,mkx)
              XQ      (I,mkx)  = QO (I,mkx)
              XT_Grell(I,mkx)  = TN (I,mkx)
              if (XQ(I,mkx).le.0.) XQ(I,mkx) = 1.E-08
           endif
        enddo

        !--- Calculate moist static energy, heights, qes

        call cup_env(j, xz, xqes, xhe, xhes, xt_grell,  &
             xq, po, z1, mix, mgmxp, mkx, mgmzp, istart,    &
             iend, psur, ierr, tcrit, 2)

        !--- Environmental values on cloud levels

        call cup_env_clev(j, xt_grell, xqes, xq, xhe,   &
             xhes, xz, po, xqes_cup, xq_cup, xhe_cup,       &
             xhes_cup, xz_cup, po_cup, gamma_cup, xt_cup,   &
             psur, mix, mgmxp, mkx, mgmzp, istart, iend,    &
             ierr, z1)

        !**************************** Static Control

        !--- Moist static energy inside cloud

        do i=istart,iend
           if (ierr(i).eq.0) then
              xhkb(i)=xhe(i,k22(i))
           endif
        enddo
        call cup_up_he(k22, xhkb, xz_cup, cd, mentr_rate,  &
             xhe_cup, xhc, mix, mgmxp, mkx, mgmzp, kbcon,  &
             ierr, istart, iend, xdby, xhe, xhes_cup)

        !--- Normalized mass flux profile

        call cup_up_nms(xzu, xz_cup, mentr_rate, cd, kbcon,  &
             ktop, mix, mgmxp, mkx, mgmzp, istart, iend,     &
             ierr, k22)
        call cup_dd_nms(xzd, xz_cup, cdd, mentrd_rate, jmin, &
             ierr, mix, mgmxp, mkx, mgmzp, istart, iend, 1,  &
             kdet, z1)

        !--- Moisture downdraft

        call cup_dd_he(xhes_cup, xzd, xhcd, xz_cup, cdd,     &
             mentrd_rate, jmin, ierr, mix, mgmxp, mkx,       &
             mgmzp, istart, iend, xhe, kdet, dbyd, xhe_cup)
        call cup_dd_moisture(j, xzd, xhcd, xhes_cup, xqcd,   &
             xqes_cup, xpwd, xq_cup, xz_cup, cdd,            &
             mentrd_rate, jmin, ierr, gamma_cup, xpwev,mix,  &
             mgmxp, mkx, mgmzp, istart, iend, bu, xqrcd, xq, &
             xhe, xhc, xt_cup,3)

        !--- Moisture updraft

        call cup_up_moisture(ierr, xz_cup, xqc, xqrc, xpw,   &
             xpwav, kbcon, ktop, mix, mgmxp, mkx, mgmzp,     &
             istart, iend, cd, xdby, mentr_rate, xq,         &
             GAMMA_cup, xzu, xqes_cup, k22, xq_cup)
        !
        !--- Workfunctions for updraft
        !
        call cup_up_aa0(xaa0, xz, xzu, xdby, GAMMA_CUP,   &
             xt_cup, kbcon, ktop, mix, mgmxp, mkx, mgmzp, &
             istart, iend,ierr)

        !--- Workfunctions for downdraft
        !call cup_dd_aa0(edtx,ierr,xaa0,jmin,gamma_cup,xt_cup, &
        !     xhcd,xhes_cup,xz,mix,mgmxp,mkx,mgmzp,istart,iend,xzd)

!------------------------- loop at  cap_max ensemble ----------
      do nens=1,maxens

!----print
!      if(j==21) then 
!              print*,'##################################################'
!              print*,'===CAP MAX ENS ==============',nens
!	      endif
!------print

        do i=istart,iend 
           if (ierr(i).eq.0) then

              xaa0_ens(i,nens) = xaa0(i)
              nall=  (iens-1)*maxens3*maxens*maxens2  &
     		   + (iedt-1)*maxens3*maxens	    &
      		   + (nens-1)*maxens3

              do k=1,mkx
                 if (k.le.ktop(i)) then
                    do nens3=1,maxens3

!- pr_ens is the total precipitation on the surface associated to each ensemble member

                       if (nens3.eq.7) then
                   !--- b=0
                          pr_ens(i,j,nall+nens3) = pr_ens(i,j,nall+nens3) +  &
                                                   pwo(i,k)+edto(i)*pwdo(i,k)
                   !--- b=beta
                       else if (nens3.eq.8) then
                          pr_ens(i,j,nall+nens3) = pr_ens(i,j,nall+nens3) +  &
                                                   pwo(i,k)
                   !--- b=beta/2
                       else if (nens3.eq.9) then
                          pr_ens(i,j,nall+nens3) = pr_ens(i,j,nall+nens3) +  &
                                                   pwo(i,k)+.5*edto(i)*pwdo(i,k)
                       else
                          pr_ens(i,j,nall+nens3) = pr_ens(i,j,nall+nens3) +  &
                                                   pwo(i,k)+edto(i)*pwdo(i,k)
                       endif





                    enddo
                 endif
              enddo

              if(pr_ens(i,j,nall+7) .lt. 0.) then
                 ierr(i) = 18
                 do nens3=1,maxens3
                    pr_ens(i,j,nall+nens3)= 0.
                 enddo
	      endif

	      do nens3=1,maxens3
	         if(pr_ens(i,j,nall+nens3) .lt. 0.) then
                    pr_ens(i,j,nall+nens3)= 0.
	         endif
	      enddo
              do nens3=1,maxens3
                 outt_ens(i,j,nall+nens3) = dellat(i,1)
              enddo

           endif
        enddo

     enddo
!---------------------end of loop at  cap_max ensemble ----------

     !--- Large scale forcing

     !   here, ,should just take top of PBL for shallow clouds, also maybe for
     !   deep clouds in tropics: kstart=level(pbltop)
     !kstart = 3
     !call MAXIMI(HE_CUP, mix, mgmxp, mkx, mgmzp,kstart,   &
     !     KBMAX, K22x, ISTART, IEND, ierr)
     
     do I=ISTART,IEND
        if (ierr(i).eq.0) then
           !if (K22x(I).ge.KBMAX(i)) ierr(i)=998
	   k22x(i)=k22(i) 
        endif
        ierr2(i)=ierr(i) 
        ierr3(i)=ierr(i) 
     enddo

     !--- Determine the level of convective cloud base - KBCON

     call cup_kbcon_catt(cap_max_increment,2,k22x,kbconx,heo_cup,heso_cup, &
          hkbo,kzi,mix,mgmxp,mkx,mgmzp,istart,iend,ierr2,kbmax,po_cup,cap_max,&
	  j)

     call cup_kbcon_catt(cap_max_increment,3,k22x,kbconx,heo_cup,heso_cup, &
          hkbo,kzi,mix,mgmxp,mkx,mgmzp,istart,iend,ierr3,kbmax,po_cup,cap_max,&
	  j)

     if(ensdim.eq.1)then 

        call cup_forcing_ens_1_catt(aa0,aa1,xaa0_ens,mbdt_ens,dtime,    &
             xmb,ierr,mix,mgmxp,mjx,mgmyp,mkx,mgmzp,istart,iend,xf_ens, &
             j,fquasi,fstab,'deeps',xland,maxens,iens,iedt,maxens2,     &
             maxens3,mconv,omeg,zdo,kbcon,zuo,pr_ens,edto,aad,kbcon,    &
             massflx,iact_old_gr,direction,ensdim, &
             massfln,xff_ens3, xk,ierr2,ierr3, i_cup_dir_flag)

     elseif(maxens3.eq.16) then

        call cup_forcing_ens_16_catt(aa0,aa1,xaa0_ens,mbdt_ens,dtime,     &
             xmb,ierr,mix,mgmxp,mjx,mgmyp,mkx,mgmzp,istart,iend,xf_ens,j, &
             fquasi,fstab,'deeps',xland,maxens,iens,iedt,maxens2,         &
             maxens3,mconv,omeg,zdo,kbcon,zuo,pr_ens,edto,aad,kbcon,       &
             massflx,iact_old_gr,direction,ensdim,                        &
             massfln,massfld,iresult,xff_ens3, xk,p_cup,ktop,icoic,ierr2, &
             ierr3,ngrid, &
	     sgrell1_2d,m2,m3, tscl_KF, i_cup_dir_flag)       

     endif

     do k=1,mkx
        do i=istart,iend
           if (ierr(i).eq.0) then
              dellat_ens (i,k,iedt)  =  dellat(i,k)
              dellaq_ens (i,k,iedt)  =  dellaq(i,k)
              dellaqc_ens(i,k,iedt) =  dellaqc(i,k)
              pwo_ens    (i,k,iedt)  = pwo(i,k) + edt(i)*pwdo(i,k)

!srf-print-------------------------
!      if(i==47 .and. j==21) then
!       if(k==1)print*,'==========================================================='
!       if(k==1)print*,'=======         pr_ens:   IEDT             ================'
!	if(k==1)   print*,'IEDT K=',iedt,k,edt(i)
!	   print*,k,pwo(i,k),pwdo(i,k), pwo_ens    (i,k,iedt)
!       if(k==mkx)print*,'========================================================='
!      endif
!srf-print-------------------------




           else 
              dellat_ens (i,k,iedt) = 0.
              dellaq_ens (i,k,iedt) = 0.
              dellaqc_ens(i,k,iedt) = 0.
              pwo_ens    (i,k,iedt) = 0.
           endif
        enddo
     enddo
     !250  CONTINUE
  enddo

  !--- FEEDBACK

  call cup_output_ens_catt(xf_ens,ierr,dellat_ens,dellaq_ens,dellaqc_ens,   &
       outt,outq,outqc,pre,pwo_ens,xmb,ktop,mix,mgmxp,mjx,mgmyp,mkx,mgmzp,  &
       istart,iend,j,'deep',maxens2,maxens,iens,pr_ens,outt_ens,    &
       maxens3,ensdim,icoic,massfln,xfac1,xfac_for_dn, &
       sgrell1_3d,sgrell2_3d,m1,m2,m3,ngrid,dti)


  do i=istart,iend
     pre(i) = max(pre(i),0.)
  enddo
 

  !
  !---------------------done----------------------------------------------------------------

  !---Salva parametros para uso no transporte convectivo:
  do i=istart,iend
     ierr4d(i)   = ierr(i)
     jmin4d(i)   = jmin(i)
     kdet4d(i)   = kdet(i)
     k224d(i)    = k22(i)
     kbcon4d(i)  = kbcon(i)
     ktop4d(i)   = ktop(i)
     kstabi4d(i) = kstabi(i)
     kstabm4d(i) = kstabm(i)
     !kpbl4d(i)   = kpbl(i)
     kpbl4d(i)   = k22(i) ! por enquanto kpbl==k22
     xmb4d(i)    = xmb(i)

     do k=1,mkx
        if (iens.eq.1) then
           zcup5d(k,i) = zo_cup(i,k)
           pcup5d(k,i) = po_cup(i,k)
        endif
        enup5d(k,i) = mentr_rate
        endn5d(k,i) = mentrd_rate
        deup5d(k,i) =   cd(i,k)
        dedn5d(k,i) =  cdd(i,k)
         zup5d(k,i)  = zuo(i,k)
         zdn5d(k,i)  = zdo(i,k)
        ! p-lw5d is the ratio between precip 
        ! rate/liq water content after rainout which
        ! unit is:  kg[ar]/(m^2s)	   
        ! It's save for use at wet-deposition scheme
        ! (in-cloud)
	prup5d (k,i) = xmb(i)*pwo(i,k) !only for upfradt
	clwup5d(k,i) =  qrco(i,k)      !only for upfradt
	tup5d  (k,i) = t_cup(i,k)  !>>> em verdade deveria ser a temperatura da parcela
                                   !>>> de ar no updraft e _NAO_ a temperatura ambiente
     enddo
  enddo

  !-------Salva parametros nas analises do RAMS

  do I=ISTART,IEND
     xierr(i,j)=float(ierr(i))
     if (ierr(i).eq.0) then
        upmf(i,j)   = xmb(i)

        ! downdraft mass flux averaged
        dnmf(i,j)  = 0.
        do k=1,ensdim
           dnmf(i,j) = dnmf(i,j) + massfln(i,j,k)
        enddo
        dnmf(i,j) = dnmf(i,j)/float(ensdim)
        !recalcula o parametro edt4d para o transporte convectivo 
        !modifique posteriormente para uso direto do fluxo de massa do downdraft
        !dnmf(i,j) = edt4d(i)*xmb(i)    ! downdraft mass flux averaged
        edt4d(i)  = dnmf(i,j) / ( upmf(i,j) + 1.e-16 )
        !     
        xktop(i,j ) = float(ktop(i) )
        xkbcon(i,j) = float(kbcon(i))
        xkdt(i,j  ) = float(kdet(i) )
        xjmin(i,j ) = float(jmin(i) )
	xk22(i,j  ) = float(k22(i)  )    
     elseif (ierr(i).ne.0.and.ierr(i).ne.20) then
        upmf(i,j)   = 0.
        dnmf(i,j)   = 0.
        xktop(i,j)  = 0.
        xkbcon(i,j) = 0.
        xkdt(i,j)   = 0.
        xjmin(i,j)  = 0.
        xk22(i,j)   = 0. 
     endif
  enddo

end subroutine CUP_enss_catt
!--------------------------------------------------------------------------------

subroutine cup_output_ens_catt(xf_ens,ierr,dellat,dellaq,dellaqc,           &
           outt,outq,outqc,pre,pwo_ens,xmb,ktop,mix,mgmxp,mjx,mgmyp,        &
           mkx,mgmzp,istart,iend,j,name,maxens2,maxens,iens,pr_ens, &
           outt_ens,maxens3,ensdim,icoic,massfln,xfac1,xfac_for_dn,         &
           sgrell1_3d,sgrell2_3d,m1,m2,m3,ngrid, dti)
 

  use cup_output_vars, only: &
       xmb_ave,xmb_std,xmb_ske,xmb_cur, &
        pr_ave, pr_std, pr_ske, pr_cur, &
	 x_ave,  x_std,  x_ske,  x_cur, &
	 x_ave_cap,                     &
	 x_ave_cap1,x_ave_cap2,x_ave_cap3,x_ave_cap4,x_ave_cap5, &
        cup_output_vars_alloc,          &
        alloc_cup_output_vars 
  use rconstants, only: day_sec

  implicit none
  ! xf = ensemble mass fluxes
  ! pr_ens = surface precipitation ensembles
  ! dellat = change of temperature per unit mass flux of cloud ensemble
  ! dellaq = change of q per unit mass flux of cloud ensemble
  ! dellaqc = change of qc per unit mass flux of cloud ensemble
  ! outt = output temp tendency (per s)
  ! outq   = output q tendency (per s)
  ! outqc  = output qc tendency (per s)
  ! pre    = output precip
  ! xmb    = total base mass flux
  ! xfac1  = correction factor
  ! pwo_ens = pw -epsilon*pd (ensemble dependent)
  ! ierr error value, maybe modified in this routine
  !
  character (LEN=*) name  !CHARACTER *(*) name

  integer mix,mjx,mkx,istart,iend,mgmxp,mgmyp,mgmzp,    &
          ensdim,i,k,j,maxens2,n,maxens,        &
          icoic,m1,m2,m3,ngrid,ncount,iens,maxens3  

  integer ktop(mgmxp),ierr(mgmxp)

  real outtes,ddtes,dtt,dtq,dtqc,dtpw

  real :: xf_ens(mgmxp,mgmyp,ensdim), pr_ens(mgmxp,mgmyp,ensdim) & 
       ,outt_ens(mgmxp,mgmyp,ensdim),massfln(mgmxp,mgmyp,ensdim)

  real outt(mgmxp,mgmzp),outq(mgmxp,mgmzp),outqc(mgmxp,mgmzp),       &
       dellat(mgmxp,mgmzp,maxens2), dellaq(mgmxp,mgmzp,maxens2),     &
       pwo_ens(mgmxp,mgmzp,maxens2),dellaqc(mgmxp,mgmzp,maxens2),     &
       pre(mgmxp),xmb(mgmxp),xfac1(mgmxp),xfac_for_dn(mgmxp) 

  real, dimension(m1,m2,m3) :: sgrell1_3d,sgrell2_3d
  real dti
  real tunning
  data tunning /0./

  ! Statistics properties
  integer i_use_stat_prop
  data  i_use_stat_prop/1/
 
  do K=1,MKX
     do I=ISTART,IEND
        outt(i,k) = 0.
        outq(i,k) = 0.
       outqc(i,k) = 0.
     enddo
  enddo

  do I=ISTART,IEND
     pre(i)  =0.
     xmb(i)  =0.
     xfac1(i)=1.
     xfac_for_dn(i) = 1.
  enddo
  !
  !--- Calculate ensemble average mass fluxes
  !
  !Lufla begin -------------
! if(i_use_stat_prop == 1 .and. icoic == 0 ) then
  if(i_use_stat_prop == 1                  ) then

     if(.not. cup_output_vars_alloc)  &
          call alloc_cup_output_vars(mgmxp,maxens,maxens3)


     !1st time: Updraft mass fluxes statistics properties
     !          sgrell1_3d and sgrell2_3d arrays will store them for output
     call massflx_stats_catt( &
          xf_ens , &	      !01
          mix    , &	    !02
          mgmxp  , &	    !03
          mjx    , &	    !04
          mgmyp  , &	    !05
          ensdim , &	    !06
          maxens , &	    !07
          maxens2, & 	     !08
          maxens3, &	     !09
          xmb_ave, &	     !10
          !                  
          xmb_std, &	     !11
          xmb_cur, &	     !12
          xmb_ske, &	     !13
          x_ave  , &	   !14
          x_std  , &	   !15
          x_ske  , &	   !16
          x_cur  , &	   !17
          x_ave_cap, &	     !18
          j,       &	     !19
          istart,  &	     !20
          !
          iend,    &	     !21
          ierr,    &	     !22
          m1,      & 	     !23
          m2,      & 	     !24
          m3,      & 	     !25
          sgrell1_3d, &	     !26
	  sgrell2_3d, &   
	  1,	      &
	  x_ave_cap1, &
	  x_ave_cap2, &
	  x_ave_cap3, &
	  x_ave_cap4, &
	  x_ave_cap5,  &
	  dti)	    

     !2nd time: precipitation statistics properties
     !sgrell1_3d array will store them for output
     call massflx_stats_catt( &           
          pr_ens, &          !01
          mix, &	     !02
          mgmxp, &	     !03
          mjx, &	     !04
          mgmyp, &	     !05
          ensdim, &	     !06
          maxens, &	     !07
          maxens2, & 	     !08
          maxens3, &	     !09
          pr_ave, &	     !10
          !
          pr_std, &	     !11
          pr_cur, &	     !12
          pr_ske, &	     !13
          x_ave, &	     !14
          x_std, &	     !15
          x_ske, &	     !16
          x_cur, &	     !17
          x_ave_cap, &	     !18
          j, &		     !19
          istart, &	     !20
          !
          iend, &	     !21
          ierr, &	     !22
          m1, & 	     !23
          m2, & 	     !24
          m3, & 	     !25
          sgrell1_3d, &	     !26
	  sgrell2_3d, &   
	  2,	      &
	  x_ave_cap1, &
	  x_ave_cap2, &
	  x_ave_cap3, &
	  x_ave_cap4, &
	  x_ave_cap5, &
	  dti)	    

     do i=istart,iend

        if(ierr(i).eq.0)then

           ncount = 0
           do n=(iens-1)*maxens2*maxens*maxens3+1, &
                    iens*maxens2*maxens*maxens3
                pr_ens(i,j,n) =   pr_ens(i,j,n)*xf_ens(i,j,n)
              outt_ens(i,j,n) = outt_ens(i,j,n)*xf_ens(i,j,n)

              xmb(i) = xmb(i) + xf_ens(i,j,n)
              ncount = ncount + 1
           enddo
           xfac_for_dn(i) = xmb(i)/float(ncount)

           !tunning process
           xmb(i) = xmb_ave(i) - tunning * xmb_std(i)
           if(xmb(i).lt.0.) xmb(i)= .1 * xmb_ave(i)

           xfac1(i) = xmb(i)
           !srf - inconsistency at downdraft mass flux due tunning process
           !      on xmb - correction factor:
           xfac_for_dn(i) = xmb(i)/(xfac_for_dn(i) + 1.e-16)

           if(xmb_ave(i) .lt. 1.e-10) ierr(i) = 11

        endif
     enddo

  else

     !----  Simple average       
     do I=ISTART,IEND
        ncount=0
        xmb(i)=0.

        if(ierr(i).eq.0)then

           do n=(iens-1)*maxens2*maxens*maxens3+1,&
	         iens   *maxens2*maxens*maxens3
                pr_ens(i,j,n) =   pr_ens(i,j,n)*xf_ens(i,j,n)
              outt_ens(i,j,n) = outt_ens(i,j,n)*xf_ens(i,j,n)
              !
              !srf- esta restricao gera incompatibilidade com o calculo de 
              !     do fluxo de massa do  downdraft calculado fora
              !IF(xf_ens(i,j,n).GT.0.)THEN !Lufla
              if(xf_ens(i,j,n).ge.0.)then
                 xmb(i) = xmb(i) + xf_ens(i,j,n)
                 ncount = ncount + 1 
             endif
           enddo

           if(ncount.gt.0)then
              xmb(i)=xmb(i)/float(ncount)
           else
              xmb(i)=0.
              ierr(i)=13
           endif
        endif
        xfac1(i)=xmb(i) 
     enddo
  end if
  !
  !-- Now do feedback
  !
  ddtes=250.
  !     if(name.eq.'shal')ddtes=500.
!-srf from WRF
! ddtes=500.
!     if(name.eq.'shal')ddtes=200.
!
  do k=1,mkx  
     do i=istart,iend
        dtt  =0.
        dtq  =0.
        dtqc =0.
        dtpw =0.

        if(ierr(i).eq.0.and.k.le.ktop(i))then
           do n=1,maxens2
              dtt  = dtt  +  dellat(i,k,n)
              dtq  = dtq  +  dellaq(i,k,n)
              dtqc = dtqc + dellaqc(i,k,n)
              dtpw = dtpw + pwo_ens(i,k,n)

           enddo

           outtes = dtt*XMB(I)*day_sec/float(maxens2)

           if (outtes .gt. 2.*ddtes .and. k.gt.2) then
              XMB(I) = 2.*ddtes/outtes * xmb(i)
              outtes = 1.*ddtes
           endif

           if (outtes .lt. -ddtes)                then
              XMB(I) = -ddtes/outtes * xmb(i)
              outtes = -ddtes
           endif

           if (outtes .gt. .5*ddtes .and. k.le.2) then
              XMB(I) =    ddtes/outtes * xmb(i)
              outtes = .5*ddtes
           endif

            outt(i,k) = outt(i,k) +xmb(i)*dtt /float(maxens2)
            outq(i,k) = outq(i,k) +xmb(i)*dtq /float(maxens2)
           outqc(i,k) =outqc(i,k) +xmb(i)*dtqc/float(maxens2)
           pre(i)      = pre(i)   +xmb(i)*dtpw/float(maxens2) ! unit : kg[liq water]/(m^2 s)
        endif
     enddo
  enddo


  do I=ISTART,IEND
     if(ierr(i).eq.0)then
!srf 20fev2003
!          xfac1(i) = xmb(i)/xfac1(i)  = sometimes getting NAN
        xfac1(i)=xmb(i)/(xfac1(i)+1.e-16)
        do k=(iens-1)*maxens2*maxens*maxens3+1,iens*maxens2*maxens*maxens3
            massfln(i,j,k)=  massfln(i,j,k)*xfac1(i)*xfac_for_dn(i)
             pr_ens(i,j,k)=   pr_ens(i,j,k)*xfac1(i)
           outt_ens(i,j,k)= outt_ens(i,j,k)*xfac1(i)
        enddo
     endif
  enddo

end subroutine cup_output_ens_catt
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

subroutine massflx_stats_catt( &
     xf_ens, &       !01
     mix, &          !02
     mgmxp, &        !03
     mjx, &          !04
     mgmyp, &        !05
     ensdim, &       !06
     maxens, &       !07
     maxens2, &      !08
     maxens3, &      !09
    !
     xt_ave, &       !10
     xt_std, &       !11
     xt_cur, &	  !12
     xt_ske, &	  !13
     x_ave,  &	  !14
     x_std,  &	  !15
     x_ske,  &	  !16
     x_cur,  &	  !17
     x_ave_cap, &	  !18
     j, &		  !19
     istart, &	  !20
     !
     iend, &         !21
     ierr, &	  !22
     m1, & 	  !23
     m2, & 	  !24
     m3, & 	  !25
     sgrell1_3d, &	  !26
     sgrell2_3d, &   
     itest,      &  !27
     x_ave_cap1, &
     x_ave_cap2, &
     x_ave_cap3, &
     x_ave_cap4, &
     x_ave_cap5, &
     dti)

  implicit none

  integer :: mix,mgmxp,mjx,mgmyp,ensdim,maxens   !,mkx
  integer :: maxens2,maxens3,i,j,k,istart,iend,num,ielem,ncount
  integer :: kk,m1,m2,m3,num2,num3,iedt,itest
  integer :: ierr(mgmxp)
  real :: xf_ens(mgmxp,mgmyp,ensdim)
  real :: xt_ave(mgmxp),xt_std(mgmxp),xt_ske(mgmxp),xt_cur(mgmxp)
  real :: x_ave(mgmxp,maxens3),x_std(mgmxp,maxens3)
  real :: x_ske(mgmxp,maxens3),x_cur(mgmxp,maxens3)
  real :: x_ave_cap(mgmxp,maxens)

  real :: small_number,rxxx


!08-03-2004
!   extra arrays for outputting precip from cap ensembles in dependence of closures
!
  real :: x_ave_cap1(mgmxp,maxens), x_ave_cap2(mgmxp,maxens) &
         ,x_ave_cap3(mgmxp,maxens), x_ave_cap4(mgmxp,maxens) &
         ,x_ave_cap5(mgmxp,maxens)
  real dti,pw_cap

!- for output at analysis only
  real, dimension(m1,m2,m3) :: sgrell1_3d,sgrell2_3d

  num =maxens *maxens2 ! old way: ensdim/maxens3
  num2=maxens2*maxens3 ! old way: ensdim/maxens
  num3=maxens2         ! old way: num2/maxens3
  small_number=1.e-6
 
 
  do k=1,maxens
     do i=istart,iend
        x_ave_cap (i,k)=0.
        x_ave_cap1(i,k)=0.
        x_ave_cap2(i,k)=0.
        x_ave_cap3(i,k)=0.
        x_ave_cap4(i,k)=0.
        x_ave_cap5(i,k)=0.
     enddo
  enddo

  do k=1,maxens3
     do i=istart,iend
        x_ave(i,k)=0.
        x_std(i,k)=0.
        x_ske(i,k)=0.
        x_cur(i,k)=0.
     enddo
  enddo

  do i=istart,iend
     xt_ave(i)=0.
     xt_std(i)=0.
     xt_ske(i)=0.
     xt_cur(i)=0.
  enddo



!srf-                The mean for each cap_max (1,2 and 3) :
!     average in maxen2 and maxens3 ensemble elements for each
!     maxens (cap_max) ensemble element 
  do iedt=1,maxens2
     do   k=1,maxens
        do kk=1,maxens3
           do i=istart,iend
              if(ierr(i).eq.0) &
                 x_ave_cap(i,k) = x_ave_cap(i,k) + &
		                  xf_ens(i,j,maxens3*(k-1)+(iedt-1)*maxens*maxens3+kk)
           enddo
        enddo
     enddo
  enddo

!srf-  The mean for each cap_max mean on the closure group (grell, AS, KF, etc)
!                                     and the precip efficiency
!      average in maxen2 and maxens3 ensemble elements for each
!      maxens (cap_max) ensemble element 
! x_ave_cap1(k) = fluxo de massa media sobre os tres elementos
!                 de fechamento grell e sobre os tres elementos 
!                 da eficiencia de precipitacao para cada 
!                 cap_max (k=1,2,3)
  do iedt=1,maxens2
    do   k=1,maxens
      do i=istart,iend
       if(ierr(i).eq.0) then
          x_ave_cap1(i,k) = x_ave_cap1(i,k) +                        & 
     	   (xf_ens(i,j,maxens3*(k-1)+(iedt-1)*maxens*maxens3+1)+     & 
     	    xf_ens(i,j,maxens3*(k-1)+(iedt-1)*maxens*maxens3+2)+     & 
     	    xf_ens(i,j,maxens3*(k-1)+(iedt-1)*maxens*maxens3+3))*.33333

          x_ave_cap2(i,k) = x_ave_cap2(i,k) +                        &
     	    (xf_ens(i,j,maxens3*(k-1)+(iedt-1)*maxens*maxens3+4)+    &
     	     xf_ens(i,j,maxens3*(k-1)+(iedt-1)*maxens*maxens3+5)+    &
     	     xf_ens(i,j,maxens3*(k-1)+(iedt-1)*maxens*maxens3+6))*.33333

          x_ave_cap3(i,k) = x_ave_cap3(i,k) +                       &
     	    (xf_ens(i,j,maxens3*(k-1)+(iedt-1)*maxens*maxens3+7)+   &
     	     xf_ens(i,j,maxens3*(k-1)+(iedt-1)*maxens*maxens3+8)+   &
     	     xf_ens(i,j,maxens3*(k-1)+(iedt-1)*maxens*maxens3+9))*.33333

          x_ave_cap4(i,k) = x_ave_cap4(i,k) +                       &
     	    (xf_ens(i,j,maxens3*(k-1)+(iedt-1)*maxens*maxens3+10)+  &
     	     xf_ens(i,j,maxens3*(k-1)+(iedt-1)*maxens*maxens3+11)+  &
     	     xf_ens(i,j,maxens3*(k-1)+(iedt-1)*maxens*maxens3+12))*.33333

          x_ave_cap5(i,k) = x_ave_cap5(i,k) +                        &
     	   (xf_ens(i,j,maxens3*(k-1)+(iedt-1)*maxens*maxens3+13)+    &
!    	     xf_ens(i,j,maxens3*(k-1)+(iedt-1)*maxens*maxens3+14)+
     	    xf_ens(i,j,maxens3*(k-1)+(iedt-1)*maxens*maxens3+15)+    &
!    	     xf_ens(i,j,maxens3*(k-1)+(iedt-1)*maxens*maxens3+16))*.25
     	    xf_ens(i,j,maxens3*(k-1)+(iedt-1)*maxens*maxens3+16))*.33333

       endif
    enddo
   enddo
  enddo

  do k=1,maxens
     do i=istart,iend
        if(ierr(i).eq.0) then
	   x_ave_cap (i,k)=x_ave_cap (i,k)/float(num2)
           x_ave_cap1(i,k)=x_ave_cap1(i,k)/float(num3)
           x_ave_cap2(i,k)=x_ave_cap2(i,k)/float(num3)
           x_ave_cap3(i,k)=x_ave_cap3(i,k)/float(num3)
           x_ave_cap4(i,k)=x_ave_cap4(i,k)/float(num3)
           x_ave_cap5(i,k)=x_ave_cap5(i,k)/float(num3)

	 endif
     enddo
  enddo


!srf- mass flux average for each closure:
!     average in maxens and maxens2 ensemble elements for each
!     maxens3 (closure) ensemble element
!     
  do kk=1,num ! num = maxens * maxens2
     do k=1,maxens3
        do i=istart,iend
           if(ierr(i).eq.0)then
              x_ave(i,k)=x_ave(i,k) + xf_ens(i,j,maxens3*(kk-1)+k)
           endif
        enddo
     enddo
  enddo

  do k=1,maxens3
     do i=istart,iend
        if(ierr(i).eq.0) x_ave(i,k) = x_ave(i,k)/float(num)
     enddo
  enddo

!srf- total average in maxens, maxen2 and maxens3 ensemble elements
!new way
  do k=1,maxens3
     do i=istart,iend
        if(ierr(i).eq.0) xt_ave(i) = xt_ave(i) + x_ave(i,k)/float(maxens3)
     enddo
  enddo


  !--- now do std, skewness,curtosis
  !srf- stats properties in maxens and maxens2 ensemble elements for each
  !     maxens3 (closure) ensemble element (for each closure)
 do kk=1,num
     do k=1,maxens3
        do i=istart,iend
           if(ierr(i).eq.0.and.x_ave(i,k).gt.0.)then
!-----------------------    bug- underflow
              rxxx       = max( small_number, xf_ens(i,j,maxens3*(kk-1)+k)-x_ave(i,k))
             !x_std(i,k) = x_std(i,k)+(xf_ens(i,j,maxens3*(kk-1)+k)-x_ave(i,k))**2
              x_std(i,k) = x_std(i,k)+(rxxx)**2
	      
             !x_ske(i,k) = x_ske(i,k)+(xf_ens(i,j,maxens3*(kk-1)+k)-x_ave(i,k))**3
              x_ske(i,k) = x_ske(i,k)+(rxxx)**3
	      
	     !x_cur(i,k) = x_cur(i,k)+(xf_ens(i,j,maxens3*(kk-1)+k)-x_ave(i,k))**4
	      x_cur(i,k) = x_cur(i,k)+(rxxx)**4
!-orig
!             x_std(i,k) = x_std(i,k)+(xf_ens(i,j,maxens3*(kk-1)+k)-x_ave(i,k))**2
!             x_ske(i,k) = x_ske(i,k)+(xf_ens(i,j,maxens3*(kk-1)+k)-x_ave(i,k))**3
!	      x_cur(i,k) = x_cur(i,k)+(xf_ens(i,j,maxens3*(kk-1)+k)-x_ave(i,k))**4
!-----------------------
              !
              xt_std(i)  = xt_std(i) +(xf_ens(i,j,maxens3*(kk-1)+k)-xt_ave(i) )**2

           endif
        enddo
     enddo
  enddo



  do k=1,maxens3
     do i=istart,iend
        if(ierr(i).eq.0.and.xt_ave(i).gt.0.)then
           !        xt_std(i) = xt_std(i)+(x_ave(i,k)-xt_ave(i))**2
           xt_ske(i) = xt_ske(i)+(x_ave(i,k)-xt_ave(i))**3
	   
!bug- underflow
            if(x_ave(i,k)-xt_ave(i) > small_number) &
            xt_cur(i) = xt_cur(i)+(x_ave(i,k)-xt_ave(i))**4
!           xt_cur(i) = xt_cur(i)+(x_ave(i,k)-xt_ave(i))**4
        endif
     enddo
  enddo
  do k=1,maxens3
     do i=istart,iend
!srf-21/02/2005
!       if(ierr(i).eq.0.and.x_std(i,k).gt.0.)then
        if(ierr(i).eq.0.and.x_std(i,k).gt.small_number)then
           x_std(i,k) = x_std(i,k)/float(num)
           x_std(i,k) = sqrt(x_std(i,k))
           if(xt_std(i) .gt. 0. ) then
              x_ske(i,k) = x_ske(i,k)/float(num)/x_std(i,k)**3
              x_cur(i,k) = x_cur(i,k)/float(num)/x_std(i,k)**4
	      
           endif
        endif
        !       print*,'                               '
        !       print*,'Some statistics at gridpoint i,j, ierr',i,j,ierr(i)
        !       print*,'statistics for closure number ',k
        !       print*,'Average= ',x_ave(i,k),'  Std= ',x_std(i,k)
        !       print*,'Skewness= ',x_ske(i,k),' Curtosis= ',x_cur(i,k)
        !       print*,'                               '
     enddo
  enddo

  do i=istart,iend
     if(ierr(i).eq.0.and.xt_ave(i).gt.0.)then
        !           xt_std(i)=xt_std(i)/float(maxens3)
        xt_std(i)=xt_std(i)/float(num*maxens3)
        xt_std(i)=sqrt(xt_std(i))
!srf-21/02/2005
!       if(xt_std(i) .gt. 0. ) then
        if(xt_std(i) .gt. small_number ) then
           xt_ske(i)=xt_ske(i)/float(maxens3)/xt_std(i)**3
           xt_cur(i)=xt_cur(i)/float(maxens3)/xt_std(i)**4
        endif
        !       print*,' 			      '
        !       print*,'Total ensemble independent statistics at i j =',i,j
        !       print*,'Average= ',xt_ave(i),'  Std= ',xt_std(i)
        !       print*,'Skewness= ',xt_ske(i),' Curtosis= ',xt_cur(i)
        !       print*,' 			      '
     endif
     !srf -- store at analysis files (for output)
     !-------------------For Deep Cumulus -----------------------
     if(itest.eq.1)then
        if(maxens3.ne.16 .or. m1 .lt. 30) then
           print*,'**** maxens3 for deep different of 16 *****'
           print*,'****               OR                 *****'
           print*,'****       NZP   less than 30         *****'
           stop
        endif
        !this the total mean mass flux (average over all ensemble)
        sgrell1_3d(2,i,j) = xt_ave(i)
        sgrell1_3d(3,i,j) = xt_std(i)
        sgrell1_3d(4,i,j) = xt_ske(i)
        sgrell1_3d(5,i,j) = xt_cur(i)

        !this the mass flux for the closure type 1  - Grell
        sgrell1_3d(6,i,j) = .333*(x_ave(i,1)+x_ave(i,2)+x_ave(i,3))
        !this the mass flux for the closure type 4  - low level omega
        sgrell1_3d(7,i,j) = .333*(x_ave(i,4)+x_ave(i,5)+x_ave(i,6))
        !this the mass flux for the closure type 7  - moisture convergence
        sgrell1_3d(8,i,j) = .333*(x_ave(i,7)+x_ave(i,8)+x_ave(i,9))
        !this the mass flux for the closure type 10  - Stability closure (FC-KF)
        sgrell1_3d(9,i,j) = .333*(x_ave(i,10)+x_ave(i,11)+x_ave(i,12))
        !this the mass flux for the closure type 13  - Arakawa-Schubert
        sgrell1_3d(10,i,j) = .25*(x_ave(i,13)+x_ave(i,14)+x_ave(i,15)+x_ave(i,16))

        !this the average mass flux for each cap_max ensemble element
        do k=1,maxens
           sgrell1_3d(16+(k-1),i,j) = x_ave_cap(i,k)
        enddo

         sgrell2_3d(1,i,j)  =x_ave_cap1(i,1)
         sgrell2_3d(2,i,j)  =x_ave_cap1(i,2)
         sgrell2_3d(3,i,j)  =x_ave_cap1(i,3)
         sgrell2_3d(4,i,j)  =x_ave_cap2(i,1)
         sgrell2_3d(5,i,j)  =x_ave_cap2(i,2)
         sgrell2_3d(6,i,j)  =x_ave_cap2(i,3)
         sgrell2_3d(7,i,j)  =x_ave_cap3(i,1)
         sgrell2_3d(8,i,j)  =x_ave_cap3(i,2)
         sgrell2_3d(9,i,j)  =x_ave_cap3(i,3)
         sgrell2_3d(10,i,j) =x_ave_cap4(i,1)
         sgrell2_3d(11,i,j) =x_ave_cap4(i,2)
         sgrell2_3d(12,i,j) =x_ave_cap4(i,3)
         sgrell2_3d(13,i,j) =x_ave_cap5(i,1)
         sgrell2_3d(14,i,j) =x_ave_cap5(i,2)
         sgrell2_3d(15,i,j) =x_ave_cap5(i,3)
        !
     elseif(itest.eq.2)then
!     sgrell1_3d(11,i,j) = .333*(x_ave(i,1)+x_ave(i,2)+x_ave(i,3))
!     sgrell1_3d(12,i,j) = .333*(x_ave(i,4)+x_ave(i,5)+x_ave(i,6))
!     sgrell1_3d(13,i,j) = .333*(x_ave(i,7)+x_ave(i,8)+x_ave(i,9))
!     sgrell1_3d(14,i,j) = .333*(x_ave(i,10)+x_ave(i,11)+x_ave(i,12))
!     sgrell1_3d(15,i,j) = .25*(x_ave(i,13)+x_ave(i,14)+x_ave(i,15)+x_ave(i,16))

     sgrell1_3d(11,i,j) = .333*(x_ave(i,1)+x_ave(i,2)+x_ave(i,3)) &
     			      *sgrell1_3d(6,i,j)*3600.
     sgrell1_3d(12,i,j) = .333*(x_ave(i,4)+x_ave(i,5)+x_ave(i,6)) &
     			      *sgrell1_3d(7,i,j)*3600.
     sgrell1_3d(13,i,j) = .333*(x_ave(i,7)+x_ave(i,8)+x_ave(i,9)) &
     			      *sgrell1_3d(8,i,j)*3600.
     sgrell1_3d(14,i,j) = .333*(x_ave(i,10)+x_ave(i,11)+x_ave(i,12)) &
     			      *sgrell1_3d(9,i,j)*3600.
     sgrell1_3d(15,i,j) = .250*(x_ave(i,13)+x_ave(i,14)+x_ave(i,15)+x_ave(i,16)) &
     			      *sgrell1_3d(10,i,j)*3600. 

!    convert from kg/m2/s to kg/m2/hour
     sgrell2_3d(1,i,j)  =x_ave_cap1(i,1)*sgrell2_3d(1 ,i,j)*3600.
     sgrell2_3d(2,i,j)  =x_ave_cap1(i,2)*sgrell2_3d(2 ,i,j)*3600.
     sgrell2_3d(3,i,j)  =x_ave_cap1(i,3)*sgrell2_3d(3 ,i,j)*3600.
     sgrell2_3d(4,i,j)  =x_ave_cap2(i,1)*sgrell2_3d(4 ,i,j)*3600.
     sgrell2_3d(5,i,j)  =x_ave_cap2(i,2)*sgrell2_3d(5 ,i,j)*3600.
     sgrell2_3d(6,i,j)  =x_ave_cap2(i,3)*sgrell2_3d(6 ,i,j)*3600.
     sgrell2_3d(7,i,j)  =x_ave_cap3(i,1)*sgrell2_3d(7 ,i,j)*3600.
     sgrell2_3d(8,i,j)  =x_ave_cap3(i,2)*sgrell2_3d(8 ,i,j)*3600.
     sgrell2_3d(9,i,j)  =x_ave_cap3(i,3)*sgrell2_3d(9 ,i,j)*3600.
     sgrell2_3d(10,i,j) =x_ave_cap4(i,1)*sgrell2_3d(10,i,j)*3600.
     sgrell2_3d(11,i,j) =x_ave_cap4(i,2)*sgrell2_3d(11,i,j)*3600.
     sgrell2_3d(12,i,j) =x_ave_cap4(i,3)*sgrell2_3d(12,i,j)*3600.
     sgrell2_3d(13,i,j) =x_ave_cap5(i,1)*sgrell2_3d(13,i,j)*3600.
     sgrell2_3d(14,i,j) =x_ave_cap5(i,2)*sgrell2_3d(14,i,j)*3600.
     sgrell2_3d(15,i,j) =x_ave_cap5(i,3)*sgrell2_3d(15,i,j)*3600.

     sgrell2_3d(16,i,j)  =sgrell2_3d(16,i,j)+ sgrell2_3d(1 ,i,j)*dti 
     sgrell2_3d(17,i,j)  =sgrell2_3d(17,i,j)+ sgrell2_3d(2 ,i,j)*dti 
     sgrell2_3d(18,i,j)  =sgrell2_3d(18,i,j)+ sgrell2_3d(3 ,i,j)*dti 
     sgrell2_3d(19,i,j)  =sgrell2_3d(19,i,j)+ sgrell2_3d(4 ,i,j)*dti 
     sgrell2_3d(20,i,j)  =sgrell2_3d(20,i,j)+ sgrell2_3d(5 ,i,j)*dti 
     sgrell2_3d(21,i,j)  =sgrell2_3d(21,i,j)+ sgrell2_3d(6 ,i,j)*dti 
     sgrell2_3d(22,i,j)  =sgrell2_3d(22,i,j)+ sgrell2_3d(7 ,i,j)*dti 
     sgrell2_3d(23,i,j)  =sgrell2_3d(23,i,j)+ sgrell2_3d(8 ,i,j)*dti 
     sgrell2_3d(24,i,j)  =sgrell2_3d(24,i,j)+ sgrell2_3d(9 ,i,j)*dti 
     sgrell2_3d(25,i,j)  =sgrell2_3d(25,i,j)+ sgrell2_3d(10,i,j)*dti 
     sgrell2_3d(26,i,j)  =sgrell2_3d(26,i,j)+ sgrell2_3d(11,i,j)*dti 
     sgrell2_3d(27,i,j)  =sgrell2_3d(27,i,j)+ sgrell2_3d(12,i,j)*dti 
     sgrell2_3d(28,i,j)  =sgrell2_3d(28,i,j)+ sgrell2_3d(13,i,j)*dti 
     sgrell2_3d(29,i,j)  =sgrell2_3d(29,i,j)+ sgrell2_3d(14,i,j)*dti 
     sgrell2_3d(30,i,j)  =sgrell2_3d(30,i,j)+ sgrell2_3d(15,i,j)*dti 
        !---------------For Shallow Cumulus ---------------------------------
    elseif(itest.eq.3)then
        if(maxens3.ne.10) then
           print*,'**** maxens3 for shallow different of 10 *****'
           stop
        endif
        !this the mass flux for the closure type 1  - Grell

        sgrell1_3d(21,i,j) = .333*(x_ave(i,1)+x_ave(i,2)+x_ave(i,3))
        !this the mass flux for the closure type 4 -  Arakawa-Schubert
        sgrell1_3d(22,i,j) = .025*(x_ave(i,4)+x_ave(i,5)+x_ave(i,6)+x_ave(i,7) )
        !this the mass flux for the closure type 7 - Stability closure (FC-KF)
        sgrell1_3d(23,i,j) = .333*(x_ave(i,8)+x_ave(i,9)+x_ave(i,10))
     endif
  enddo

end subroutine massflx_stats_catt

!--------------------------------------------------------------------------
subroutine cup_forcing_ens_1_catt(aa0,aa1,xaa0,mbdt,dtime,xmb,ierr, &
     mix,mgmxp,mjx,mgmyp,mkx,mgmzp,istart,iend,xf,j,fquasi, &
     fstab,name,xland,maxens,iens,iedt,maxens2,maxens3, &
     mconv,omeg,zd,k22,zu,pr_ens,edt,aad,kbcon,massflx, &
     iact_old_gr,dir,ensdim,massfln,xff_ens3, xk,ierr2,ierr3,i_cup_dir_flag)
  implicit none
  character *(*) name

  integer i,istart,iend,mix,mgmxp,mjx,mgmyp,mkx,mgmzp,j,maxens,maxens3
  integer ensdim,ierr(mgmxp),iens,nall,iedt,maxens2,i_cup_dir_flag
  integer ierr2(mgmxp),ierr3(mgmxp)
  integer k22(mgmxp),kbcon(mgmxp),iact_old_gr(mgmxp,mgmyp)
  real aa0(mgmxp),aa1(mgmxp),xaa0(mgmxp,maxens),xmb(mgmxp)
  real mbdt(maxens),dtime,edt(mgmxp),aad(mgmxp),dir(mgmxp)  !,dxxf
  real xf(mgmxp,mgmyp,ensdim),xland(mgmxp,mgmyp)
  real pr_ens(mgmxp,mgmyp,ensdim)
  real mconv(mgmxp),omeg(mgmxp,mgmzp),zd(mgmxp,mgmzp),zu(mgmxp,mgmzp)
  real xff_ens3(maxens3),xk(maxens),xff0      !,xff1,xff2,xff3,xff
  real massflx(mgmxp,mgmyp)
  real massfln(mgmxp,mgmyp,ensdim)
  real massfld
  integer fquasi,fstab,nens,ne,n,nens3,iresult,iresultd,iresulte
  nens=0

  !--- LARGE SCALE FORCING
  !
  !  DO 100 I=ISTART,IEND
  do i=istart,iend
     xmb(i)=0.
     if(name.eq.'deeps'.and.ierr(i).gt.995)then
        aa0(i)=0.
        ierr(i)=0
     endif
     if(ierr(i).eq.0)then
        !--- treatment different for this closure
        if(name.eq.'deeps')then
           xff0= (AA1(I)-AA0(I))/DTIME
           xff_ens3(1)=(AA1(I)-AA0(I))/dtime
           !--- more original Arakawa-Schubert (climatologic value of aa0)
           !
           !     xff_ens3(2)=max(0.,(AA1(I)-1500.)/dtime)
           !     xff_ens3(2)=xff_ens3(1)
           !
           !--- omeg is in bar/s, mconv done with omeg in Pa/s
           !     more like Brown (1979), or Frank-Cohen (199?)
           !     xff_ens3(3)=-omeg(i,k22(i))/9.81
           !     xff_ens3(4)=-omeg(i,kbcon(i))/9.81
           !--- more like Krishnamurti et al.
           !     xff_ens3(5)=mconv(i)
           !     xff_ens3(6)=mconv(i)
           !--- more like Fritsch Chappel or Kain Fritsch (plus triggers)
           !     xff_ens3(7)=AA1(I)/(60.*30.)
           !     xff_ens3(8)=AA1(I)/(60.*40.)
           do nens=1,maxens
              XK(nens)=(XAA0(I,1)-AA1(I))/MBDT(2)
              !print*,nens,XK(nens),XAA0(I,nens)
              if(xk(nens).le.0.and.xk(nens).gt.-1.e-9) xk(nens)=-1.e-9
              if(xk(nens).gt.0.and.xk(nens).lt.1.e-9) xk(nens)=1.e-9
           enddo
           !--- add up all ensembles
           !DO 350 ne=1,maxens
           do ne=1,maxens
              !--- for every xk, we have maxens3 xffs
              !--- iens is from outermost ensemble (most expensive!
              !--- iedt (maxens2 belongs to it)
              !--- is from second, next outermost, not so expensive
              !--- so, for every outermost loop, we have maxens*maxens2*3
              !--- ensembles!!! nall would be 0, if everything is on first
              !--- loop index, then ne would start counting, then iedt, then iens....
              iresultd=0
              iresulte=0
              iresulte=1
              nall= (iens-1)*maxens3*maxens*maxens2+ &
	           (iedt-1)*maxens3*maxens+ &
                   (  ne-1)*maxens3
              !--- check for upwind convectionc
              if(i_cup_dir_flag == 1 ) then
               iresult=0
               massfld=0.
               call cup_direction2(i,j,dir,iact_old_gr,mix,mjx, &
                    mgmxp,mgmyp,massflx,iresult,ensdim,1,nall, &
                    maxens3,massfld)

               if(XK(ne).lt.0.and.xff0.gt.0.) iresultd=1
               iresulte=max(iresult,iresultd)
	      endif
              
	      if(iresulte.eq.1)then
                 !--- special treatment for stability closures
                 if(xff0.gt.0.)then
                    xf(i,j,nall+1)=max(0.,-xff_ens3(1)/xk(ne))+massfld
                 else
                    xf(i,j,nall+1)=massfld
                 endif
                 !--- store new for next time step
                 do nens3=1,maxens3
                    massfln(i,j,nall+nens3)=edt(i)*xf(i,j,nall+nens3)
                    massfln(i,j,nall+nens3)=max(0.,massfln(i,j,nall+nens3))

                 enddo
              endif
           end do
           !350       continue
           !go to 100
           cycle
        end if
     elseif(ierr(i).ne.20.and.ierr(i).ne.0)then
        do n=1,ensdim
           xf(i,j,n)=0.
           massfln(i,j,n)=0.
        enddo
     end if
  end do
  !100   continue
end subroutine cup_forcing_ens_1_catt

!--------------------------------------------------------------------------

subroutine cup_dellas_catt(ierr, z_cup, p_cup, hcd, edt, zd, cdd, he, mix,   &
     mgmxp, mkx, mgmzp, istart, iend, della, itest, j, mentrd_rate, zu, &
     cd, hc, ktop, k22, kbcon, mentr_rate, jmin, he_cup, kdet, kpbl,    &
     name)
  use rconstants, only: g
  implicit none
  character (LEN=*) name  !CHARACTER *(*) name
  integer mix, mgmxp, mkx, mgmzp, i, k, istart, iend,  &
       itest, j
  real z_cup(mgmxp,mgmzp), p_cup(mgmxp,mgmzp), hcd(mgmxp,mgmzp),  &
       zd(mgmxp,mgmzp), cdd(mgmxp,mgmzp), he(mgmxp,mgmzp),        &
       della(mgmxp,mgmzp), hc(mgmxp,mgmzp), cd(mgmxp,mgmzp),      &
       zu(mgmxp,mgmzp), he_cup(mgmxp,mgmzp)
  real edt(mgmxp)
  integer kbcon(mgmxp), ktop(mgmxp), k22(mgmxp), jmin(mgmxp)
  integer ierr(mgmxp), kdet(mgmxp), kpbl(mgmxp)
  real entdo, dp, dz, mentrd_rate,                     &  !detdo1, detdo2, 
       mentr_rate, subin, detdo, entup, detup,         &
       subdown, entdoj, entupk, detupk, totmas
  real xsum,xsumt,hesum

  do k=2,mkx
     do I=ISTART,IEND
        della(i,k) = 0.
     enddo
  enddo
  !  xsum=0.
  !  hesum=0.
  !  xsumt=0.

  do K=2,MKX-1
     !DO 100 K=2,MKX-1
     do I=ISTART,IEND
        !DO 100 I=ISTART,IEND
        !IF (ierr(i).NE.0) GO TO 100
        if (ierr(i).ne.0) cycle
        !IF (K.GT.KTOP(I)) GO TO 100
        if (K.gt.KTOP(I)) cycle
        !
        !--- Specify detrainment of downdraft, has to be consistent
        !--- with zd calculations in soundd.
        !
        dz    = Z_cup(I,K+1)-Z_cup(I,K)
        detdo = edt(i)*CDD(i,K)   *dz*zd(i,k+1)
        entdo = edt(i)*mentrd_rate*dz*zd(i,k+1)
        subin = zu(i,k+1)-zd(i,k+1)*edt(i)
        entup = 0.
        detup = 0.
        if (k.ge.kbcon(i).and.k.lt.ktop(i)) then
           entup = mentr_rate*dz*zu(i,k)
           detup = CD(i,K+1) *dz*zu(i,k)
        endif
        subdown = ( zu(i,k)-zd(i,k)*edt(i) )
        entdoj  = 0.
        entupk  = 0.
        detupk  = 0.

        if (k.eq.jmin(i)) then
           entdoj  = zd(i,k)*edt(i)
        endif

        if (k.eq.k22(i)-1) then
!       if (k.eq.kpbl(i) ) then
           entupk  = zu(i,kpbl(i))
        endif

        if (k.gt.kdet(i)) then
           detdo   = 0.
        endif

        if (k.eq.ktop(i)-0) then
           detupk  = zu(i,ktop(i))
           subin   = 0.
        endif
        if (k.lt.kbcon(i)) then
           detup   = 0.
        endif
        !
        !--- Changed due to subsidence and entrainment
        !
        totmas =subin-subdown+detup-entup-entdo + detdo-entupk-entdoj+detupk
        if (abs(totmas).gt.1.e-6) then
           print *, '**totmass deep********', i, j, k, totmas, name
           !print *, kpbl(i), k22(i), kbcon(i), ktop(i)
           !          print *,'updr stuff = ',subin,
           !    1      subdown,detup,entup,entupk,detupk
           !          print *,'dddr stuff = ',entdo,
           !    1      detdo,entdoj
           stop
        endif

        !srf         dp =  100.*( p_cup(i,k-1)-p_cup(i,k) )
        dp =  100.*( p_cup(i,k)-p_cup(i,k+1) )
!        della(i,k)=(subin  *he_cup(i,k+1) - subdown*he_cup(i,k  ) +         &
!             	    detup*.5*( HC(i,K+1)+ HC(i,K)) +			    &
!             	    detdo*.5*(HCD(i,K+1)+HCD(i,K)) -			    &
!             	    entup*he(i,k) - entdo*he(i,k) -			    &
!             	    entupk*he_cup(i,k22(i)) -  entdoj*he_cup(i,jmin(i)) +   &
!             	    detupk*hc(i,ktop(i))                                    &
!		    )*g/dp
         della(i,k)=(                                 & 
     		    subin  *he_cup(i,k+1)	      & 
     		  - subdown*he_cup(i,k  )	      & 
     		  + detup*.5*( HC(i,K+1)+ HC(i,K))    & 
     		  + detdo*.5*(HCD(i,K+1)+HCD(i,K))    & 
     		  - entup*he(i,k)		      & 
     		  - entdo*he(i,k)		      & 
     		  - entupk*he_cup(i,k22(i))	      & 
     		  - entdoj*he_cup(i,jmin(i))	      & 
     		  + detupk*hc(i,ktop(i))	      & 
     		   )*g/dp

     enddo
  enddo
end subroutine cup_dellas_catt

!--------------------------------------------------------------------------

subroutine cup_dellabot_catt(he_cup, ierr, z_cup, p_cup, hcd, edt,   &
     zd, cdd, he, mix, mgmxp, mkx, mgmzp, istart, iend, della, itest, j,  &
     mentrd_rate, z)
  use rconstants, only: g
  implicit none
  integer mix, mgmxp, mkx, mgmzp, i, istart, iend, itest, j
  real z_cup(mgmxp,mgmzp), p_cup(mgmxp,mgmzp), hcd(mgmxp,mgmzp),   &
       zd(mgmxp,mgmzp), cdd(mgmxp,mgmzp), he(mgmxp,mgmzp),         &
       della(mgmxp,mgmzp), he_cup(mgmxp,mgmzp), z(mgmxp,mgmzp), edt(mgmxp)
  integer ierr(mgmxp), m

!-local variables in this routine
  real detdo1, detdo2, entdo, dp, dz, mentrd_rate, subin, detdo

  do i=istart,iend
     !DO 100 i=istart,iend
     della(i,1)=0.
     !IF (ierr(i).NE.0) GO TO 100
     if (ierr(i).ne.0) cycle
     dz        =       z_cup(i,2)-z_cup(i,1)
     dp        = 100.*(p_cup(i,1)-p_cup(i,2))
     detdo1    = edt(i)*zd(i,2)*cdd(i,1)*dz
     detdo2    = edt(i)*zd(i,1)
     entdo     = edt(i)*zd(i,2)*mentrd_rate*dz
     subin     =-edt(I)*zd(i,2)
     detdo     = detdo1+detdo2-entdo+subin

!Checking mass conservation
     if(abs(detdo) > 1.e-6) then
       write (unit=*,fmt='(a)')                 '---------- Subroutine cup_dellabot ----------'
       write(unit=*, fmt='(4(a,1x,i3,1x))')     '         K= ',1,'   I=',i,'   J=',j
       write(unit=*, fmt='(2(a,1x,es10.3,1x))') '  subin=  ',    subin,'entdo=  ',entdo
       write(unit=*, fmt='(2(a,1x,es10.3,1x))') '  detdo1= ',   detdo1,'detdo2= ',detdo2
       write(unit=*, fmt='(1(a,1x,es10.3,1x))') '  totmas= ',   detdo
       write(unit=*, fmt='(a)')                 '---------------------------------------------'
       write(unit=*, fmt='(a)')                 'The model will stop since it is not conserving mass...' 
       stop
     end if

     DELLA(I,1)= (detdo1*.5*(hcd(i,1)+hcd(i,2)) +  &
          detdo2*hcd(i,1) + subin*he_cup(i,2) -    &
          entdo*he(i,1))*g/dp

     !100  CONTINUE
  enddo
end subroutine cup_dellabot_catt

!------------------------------------------------------------

subroutine cup_forcing_ens_16_catt(aa0,aa1,xaa0,mbdt,dtime,xmb,ierr,   &
     mix,mgmxp,mjx,mgmyp,mkx,mgmzp,istart,iend,xf,j,fquasi,       &
     fstab,name,xland,maxens,iens,iedt,maxens2,maxens3,   &
     mconv,omeg,zd,k22,zu,pr_ens,edt,aad,kbcon,massflx,		  &
     iact_old_gr,dir,ensdim,massfln,massfld,iresult,xff_ens3,xk,  &
     p_cup,ktop,icoic,ierr2,ierr3,ngrid, &
     sgrell1_2d,m2,m3,tscl_KF,i_cup_dir_flag)  

  !
  ! ierr error value, maybe modified in this routine
  ! pr_ens = precipitation ensemble
  ! xf     = mass flux ensembles
  ! massfln = downdraft mass flux ensembles used in next timestep
  ! omeg = omega from large scale model
  ! mconv = moisture convergence from large scale model
  ! zd      = downdraft normalized mass flux
  ! zu      = updraft normalized mass flux
  ! aa0     = cloud work function without forcing effects
  ! aa1     = cloud work function with forcing effects
  ! xaa0    = cloud work function with cloud effects (ensemble dependent)
  ! edt     = epsilon
  ! dir     = "storm motion"
  ! mbdt    = arbitrary numerical parameter
  ! dtime   = dt over which forcing is applied
  ! iact_gr_old = flag to tell where convection was active
  ! kbcon       = LFC of parcel from k22
  ! k22         = updraft originating level
  ! icoic       = flag if only want one closure (usually set to zero!)
  ! name        = deep or shallow convection flag
  !

  use rconstants, only: g
  implicit none
  character (LEN=*) name

  integer k,i,istart,iend,mix,mgmxp,mjx,mgmyp,mkx,mgmzp,j,        &
          maxens,maxens3,ngrid

  !------ ensemble 3 dimension = 16
  integer mkxcrt,kclim
  parameter (mkxcrt=15)
  real pcrit(mkxcrt),acrit(mkxcrt),acritt(mkxcrt),aclim1,         &
       aclim2,aclim3,aclim4
  data pcrit/850.,800.,750.,700.,650.,600.,550.,500.,450.,400.,   &
       350.,300.,250.,200.,150./
  data acrit/.0633,.0445,.0553,.0664,.075,.1082,.1521,.2216,      &
       .3151,.3677,.41,.5255,.7663,1.1686,1.6851/
  !  GDAS derived acrit
  data acritt/.203,.515,.521,.566,.625,.665,.659,.688,            &
       .743,.813,.886,.947,1.138,1.377,1.896/
  integer ktop(mgmxp)
  real p_cup(mgmxp,mgmzp)

  real dec_fudge
  data  dec_fudge /0.3/
  !------

  integer fquasi,fstab,nens,ne,n,nens3,iresult,iresultd,          &
          iresulte,icoic
  integer ensdim,iens,nall,iedt,maxens2,m2,m3
  integer k22(mgmxp),kbcon(mgmxp),ierr(mgmxp),ierr2(mgmxp),ierr3(mgmxp) !Lufla
  integer iact_old_gr(mgmxp,mgmyp)

  real aa0(mgmxp),aa1(mgmxp),xaa0(mgmxp,maxens),xmb(mgmxp)
  real mbdt(maxens),dtime,edt(mgmxp),aad(mgmxp),dir(mgmxp)
  real     xf(mgmxp,mgmyp,ensdim),xland(mgmxp,mgmyp)
  real pr_ens(mgmxp,mgmyp,ensdim)
  real mconv(mgmxp),omeg(mgmxp,mgmzp),zd(mgmxp,mgmzp),            &
       zu(mgmxp,mgmzp)
  real xff_ens3(maxens3),xk(maxens),xff0   !,xff1,xff2,xff3,xff
  real massflx(mgmxp,mgmyp)
  real massfln(mgmxp,mgmyp,ensdim)
  real xomg,massfld,a1

  real, dimension(m2,m3)    :: sgrell1_2d
  real tscl_KF
  integer i_cup_dir_flag  
  nens=0

  !--- LARGE SCALE FORCING
  !
  !DO 100 I=ISTART,IEND
  do I=ISTART,IEND
     xmb(i)=0.
     if(name.eq.'deeps'.and.ierr(i).gt.995)then
        aa0(i) =0.
        ierr(i)=0
     endif
     if(ierr(i).eq.0)then
        !Added for ensemble 3 with dimension = 16
        kclim=0
        do k=mkxcrt,1,-1
           if(p_cup(i,ktop(i)).lt.pcrit(k))then
              kclim=k
              GO TO 9
           endif
        enddo
        if(p_cup(i,ktop(i)).gt.pcrit(1))kclim=1
9       continue
        k= max(kclim-1,1)
!srf: out2004
!        aclim1= acrit(kclim)*1.e3
!        aclim2= acrit(k)*1.e3
!        aclim3= acritt(kclim)*1.e3
!        aclim4= acritt(k)*1.e3
        aclim1= (acrit(kclim) -dec_fudge*acrit(kclim) )*1.e3
        aclim2= (acrit(k)     -dec_fudge*acrit(k)     )*1.e3
        aclim3= (acritt(kclim)-dec_fudge*acritt(kclim))*1.e3
        aclim4= (acritt(k)    -dec_fudge* acritt(k)   )*1.e3
!srf: --fim----

        !
        !--- Treatment different for this closure
        !
        if(name.eq.'deeps')then
        !
        !---- Grell's closure
        !
           xff0       =  (AA1(I)-AA0(I))/dtime
           xff_ens3(1)=    xff0
           xff_ens3(2)= .9*xff_ens3(1)
           xff_ens3(3)=1.1*xff_ens3(1)
	   	   

           !
           !     More like Brown (1979), or Frank-Cohen (199?)
           !
           !---  omeg is in bar/s, mconv done with omeg in Pa/s
           xff_ens3(4)=     -omeg(i,k22(i))/g       
           xff_ens3(5)=     -omeg(i,kbcon(i))/g     
           xff_ens3(6)=     -omeg(i,1)/g            

           do k=2,kbcon(i)-1
              xomg = -omeg(i,k)/g           
              if(xomg .gt. xff_ens3(6)) xff_ens3(6)=xomg
           enddo
           !
           !--- More like Krishnamurti et al.
           !
           xff_ens3(7)=    mconv(i)
           xff_ens3(8)= .9*mconv(i)
           xff_ens3(9)=1.1*mconv(i)
           !
           !--- More like Fritsch Chappel or Kain Fritsch (plus triggers)
!srf: A. Betts suggestion 
           xff_ens3(10)=AA1(I)/(tscl_KF    )
           xff_ens3(11)=AA1(I)/(tscl_KF*0.9)
           xff_ens3(12)=AA1(I)/(tscl_KF*0.8)
           !   
           !--- More original Arakawa-Schubert (climatologic value of aa0)
           !
           xff_ens3(13)=max(0.,(AA1(I)-aclim1)/dtime)
           !srf- dec/2002 - later xff_ens3(14) will be equal to xff_ens3(13)
           xff_ens3(14)=max(0.,(AA1(I)-aclim2)/dtime)
           xff_ens3(15)=max(0.,(AA1(I)-aclim3)/dtime)
           xff_ens3(16)=max(0.,(AA1(I)-aclim4)/dtime)


           do nens=1,maxens
!srf          XK(nens)=(XAA0(I,nens)-AA1(I))/MBDT(nens)
              XK(nens)=(XAA0(I,nens)-AA1(I))/MBDT(1   )
              if(xk(nens).le.0.and.xk(nens).gt.-1.e-9) xk(nens)=-1.e-9
              if(xk(nens).gt.0.and.xk(nens).lt.+1.e-9) xk(nens)=+1.e-9
           enddo
           !
           !--- Add up all ensembles
           !
           do ne=1,maxens
              nall= (iens-1)*maxens3*maxens*maxens2  &
                   +(iedt-1)*maxens3*maxens          &
                   +(ne  -1)*maxens3

              !          observe the mass flux calculation:
              !--------------------------------------------------------------!
              !           ne   |     ierr     | mass flux		     !
              !           1    |     ierr =0  |  mf1 = xff_ens3 / xk (ne)    !
              !           1    |     ierr >0  |  mf1 =  0		     !
              !           2    |     ierr2=0  |  mf2 = mf1		     !
              !           2    |     ierr2>0  |  mf2 =  0		     !
              !           3    |     ierr3=0  |  mf3 = mf1		     !
              !           3    |     ierr3>0  |  mf3 =  0		     !
              !							             !
              ! xk(ne) is the same for any 'ne'.			     !
              !--------------------------------------------------------------!
              ! if ierr2 > 0 (convection was not permited for that cap_max)
              ! then equal to zero the mass flux for the second member of the ensemble (maxens)
              if(ne.eq.2 .and. ierr2(i).gt.0)then
                 do nens3=1,maxens3
                         xf(i,j,nall+nens3)=0.
                    massfln(i,j,nall+nens3)=0.
                 enddo
                 cycle
              endif
              ! if ierr3 > 0 (convection was not permited for that cap_max)
              ! then equal to zero the mass flux for the third member of the ensemble (maxens)
              if(ne.eq.3 .and. ierr3(i).gt.0)then
                 do nens3=1,maxens3
                         xf(i,j,nall+nens3)=0.
                    massfln(i,j,nall+nens3)=0.
                 enddo
                 cycle
              endif
              !
              !--- for every xk, we have maxens3 xffs
              !--- iens is from outermost ensemble (most expensive!
              !
              !--- iedt (maxens2 belongs to it)
              !--- is from second, next outermost, not so expensive
              !
              !--- so, for every outermost loop, we have maxens*maxens2*3
              !--- ensembles!!! nall would be 0, if everything is on first
              !--- loop index, then ne would start counting, then iedt, then iens....
              !
              iresultd=0
              iresulte=0
              !
              !--- check for upwind convection
              if(i_cup_dir_flag == 1 ) then
	      !
                  iresult=0  
                  massfld=0. 
                  call cup_direction2(i,j,dir,iact_old_gr,mix,mjx,  &
                       mgmxp,mgmyp,massflx,iresult,ensdim,1,nall,   &
                       maxens3,massfld)
                  print*,'cup_dir: i j massfld= ',i,j,massfld
                  if(XK(ne).lt.0.and.xff0.gt.0.)iresultd=1
                  iresulte=max(iresult,iresultd)
	          iresulte=1             	       
	       else
	          massfld=0.
	          iresulte=1             	       
	       endif
              
	       if(iresulte.eq.1)then
                 !
                 !--- Special treatment for stability closures
                 !

                 if(xff0.gt.0.)then
                    xf(i,j,nall+1) =max(0., -xff_ens3( 1) /xk(ne))+massfld
                    xf(i,j,nall+2) =max(0., -xff_ens3( 2) /xk(ne))+massfld
                    xf(i,j,nall+3) =max(0., -xff_ens3( 3) /xk(ne))+massfld
                    xf(i,j,nall+13)=max(0., -xff_ens3(13) /xk(ne))+massfld
                    xf(i,j,nall+14)=max(0., -xff_ens3(14) /xk(ne))+massfld
                    xf(i,j,nall+15)=max(0., -xff_ens3(15) /xk(ne))+massfld
                    xf(i,j,nall+16)=max(0., -xff_ens3(16) /xk(ne))+massfld
                 else
                    xf(i,j,nall+1) =massfld
                    xf(i,j,nall+2) =massfld
                    xf(i,j,nall+3) =massfld
                    xf(i,j,nall+13)=massfld
                    xf(i,j,nall+14)=massfld
                    xf(i,j,nall+15)=massfld
                    xf(i,j,nall+16)=massfld
                 endif
                 !
                 !--- if iresult.eq.1, following independent of xff0
                 !
                 xf(i,j,nall+4)=max(0.,xff_ens3(4)+massfld)
                 xf(i,j,nall+5)=max(0.,xff_ens3(5)+massfld)
                 xf(i,j,nall+6)=max(0.,xff_ens3(6)+massfld)

                a1 = max(1.e-9,pr_ens(i,j,nall+7))
                xf(i,j,nall+7)=max(0.,xff_ens3(7)/a1)
     
                a1 = max(1.e-9,pr_ens(i,j,nall+8))
                xf(i,j,nall+8)=max(0.,xff_ens3(8)/a1)
     
                a1 = max(1.e-9,pr_ens(i,j,nall+9))
                xf(i,j,nall+9)=max(0.,xff_ens3(9)/a1)
		 
                 if(XK(ne).lt.0.)then
                    xf(i,j,nall+10)=max(0.,-xff_ens3(10)/xk(ne))+massfld
                    xf(i,j,nall+11)=max(0.,-xff_ens3(11)/xk(ne))+massfld
                    xf(i,j,nall+12)=max(0.,-xff_ens3(12)/xk(ne))+massfld
                 else
                    xf(i,j,nall+10)=massfld
                    xf(i,j,nall+11)=massfld
                    xf(i,j,nall+12)=massfld
                 endif
                 !==============
                 !05-12-2002
                 !srf - forcing 14 is too bad, use the same for 13:
                 !A&S (14) = A&S (13)
                 xf(i,j,nall+14)=xf(i,j,nall+13)
                 !
                 !----- 1d closure ensemble -------------
                 if(icoic.ge.1)then
                    do nens3=1,maxens3                      
                       xf(i,j,nall+nens3)=xf(i,j,nall+icoic)
                    enddo
                 endif
                 !
                 !--- store new for next time step
                 !
                 do nens3=1,maxens3
                    massfln(i,j,nall+nens3)=edt(i)*xf(i,j,nall+nens3)
                    massfln(i,j,nall+nens3)=max(0.,massfln(i,j,nall+nens3))
                 enddo
              endif
           enddo
           cycle !go to 100
        endif
     elseif(ierr(i).ne.20.and.ierr(i).ne.0)then
        do n=1,ensdim
           xf(i,j,n)=0.
           massfln(i,j,n)=0.
        enddo
     endif
     !100  CONTINUE
  enddo
  return
end subroutine cup_forcing_ens_16_catt
!--------------------------------------------------------------------
