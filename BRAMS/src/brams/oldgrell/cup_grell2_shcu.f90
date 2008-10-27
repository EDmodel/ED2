!srf - fev-2003
!---------------------------GRELL SHALLOW CUMULUS SCHEME---------------------------
subroutine CUPARTH_shal(mynum,mgmxp,mgmyp,mgmzp,m1,m2,m3,ia,iz,ja,jz,i0,j0      &
                       ,maxiens,iens,iupmethod,depth_min,cap_maxs,radius,zkbmax &
                       ,zcutdown,z_detr,DTIME,time,UA,VA,WA,THETA               &
                       ,PP,PI0,DN0,RV,kpbl,TKE,TKMIN,RCP,topo,RTGT,THT,RTT,PT   &
                       ,OUTTEM,OUTRT,ierr4d,jmin4d,kdet4d,k224d,kbcon4d,ktop4d  &
                       ,kpbl4d,kstabi4d,kstabm4d,xmb4d,edt4d,zcup5d,pcup5d      &
                       ,enup5d,endn5d,deup5d,dedn5d,zup5d,zdn5d,prup5d,clwup5d  &
                       ,tup5d,upmf,xierr,xktop,xkbcon,xk22,xierr_dp             &
                       )
                       
                         
  use mem_grell_param, only : maxens=>maxens_sh,  & !INTENT(IN)
       maxens2=>maxens2_sh,                        & !INTENT(IN)
       maxens3=>maxens3_sh,                        & !INTENT(IN)
       ensdim=>ensdim_sh,                          & !INTENT(IN)
       icoic=>icoic_sh                             !INTENT(IN)

  use rconstants, only: rgas,cp,rm,p00,t00,g,cpor

  use mem_scratch2_grell_sh

  implicit none
  integer mgmxp,mgmyp,mgmzp
  integer maxiens,iens,iupmethod
  real,    intent(in) :: depth_min
  real,    intent(in) :: cap_maxs
  real,    intent(in) :: radius
  real,    intent(in) :: zkbmax
  real,    intent(in) :: zcutdown
  real,    intent(in) :: z_detr
  !
  !------------------------------------ RAMS vectors:
  integer m1,m2,m3,ia,iz,ja,jz,i0,j0,j1,j2,mynum  ! ,ibcon
  real(kind=8) :: time
  real         :: tkmin
  real, dimension(m1,m2,m3) :: ua,va,wa
  real, dimension(m1,m2,m3) :: theta,pp,pi0,dn0,rv,tht,rtt, &
       pt,OUTTEM,OUTRT,tke,rcp!,SGRELL
  real, dimension(m2,m3) :: topo,rtgt,upmf,xierr,xktop,xkbcon,xk22, &
       xierr_dp
  integer, dimension(m2,m3) :: kpbl

  !-------Salva parametros da CUP para uso no transporte convectivo:
  integer, dimension(m2,m3,maxiens) :: ierr4d,jmin4d, &
                                       kdet4d,k224d,kbcon4d,ktop4d,kpbl4d, &
                                       kstabi4d,kstabm4d
  real,dimension(m2,m3,maxiens)     :: xmb4d,edt4d
  real,dimension(m1,m2,m3,maxiens)  :: enup5d,endn5d, &
       deup5d,dedn5d,zup5d,zdn5d,prup5d,clwup5d,tup5d,zcup5d,pcup5d

  !------------------------------------- variaveis locais:
  integer :: istart,iend,i,j,k,mix,mjx,mkx,kr
  real    :: dtime,cpdTdt,exner


  !Ensemble dimension
  !srf- 04-fev-2003 - mudanca no calculo DO 'ensdim' 
  integer, save :: iens_tmp
  data iens_tmp /1/
!  ensdim=       1*maxens*maxens2*maxens3
!  ensdim= maxiens*maxens*maxens2*maxens3
  !
  !srf- ICOIC is used for choice a specIFic closure for shallow cumulus
  ! icoic = 0 -> ensemble (all closures)
  ! icoic = 1 -> Grell
  ! icoic = 4 -> Arakawa-Schubert
  ! icoic = 8 -> like Fritsch Chappel or Kain Fritsch
  ! icoic = 8


  ISTART = ia
  IEND  = iz
  j1     = ja
  j2     = jz
  MKX    = m1 - 1    !MKX nao deve ser igual a m1
  MIX    = m2            
  MJX    = m3            


  ! Loop externo : j
  do J=j1,j2
     do I = ISTART,IEND
        AA0(I)=0.
        XLAND(I,J) = 0. ! land/water flag - not in use
     enddo
     !------- Transfere valores DO RAMS para o eschema
     do K=1,MKX
        kr = K + 1	  ! nivel K da grade DO Grell corresponde ao nivel K + 1 DO RAMS
        do I = ISTART,IEND
           z1(I)= topo(i,j)
           PSUR(I) = .5*( ((pp(1,i,j)+pi0(1,i,j))/cp)**cpor*p00 + &
                ((pp(2,i,j)+pi0(2,i,j))/cp)**cpor*p00 )*1.e-2 ! Pressure in mbar
           PO(I,K) = ((pp(kr,i,j)+pi0(kr,i,j))/cp)**cpor*p00*1.e-2      ! Pressure in mbar
           US(I,K) = .5*( ua(kr,i,j) + ua(kr,i-1,j) )
           VS(I,K) = .5*( va(kr,i,j) + va(kr,i,j-1) )
           OMEG(I,K)   = -g*dn0(kr,i,j)*.5*( wa(kr,i,j)+wa(kr-1,i,j) )
           T(I,K)  = theta(kr,i,j)*(pp(kr,i,j)+pi0(kr,i,j))/cp
           Q(I,K)  = rv(kr,i,j)
           !variables for PBL top height
           PBLIDX(I) = KPBL(i,j)
           TKEG(I,K) = TKE(kr,i,j)
           RCPG(I,K) = RCP(kr,i,j)	 
           !Calcula tendencia projetada na temperatura em funcao 
           !das tendencias de theta e PI : cp*T=Pi*Theta
           exner= pp(kr,i,j)+pi0(kr,i,j)	 
           !cpdTdt= exner*tht(kr,i,j) + theta(kr,i,j)*pt(kr,i,j)
           cpdTdt  = exner*tht(kr,i,j)    !assuminDO PT(KR,I,J) << exner*THT(KR,I,J)/theta
           !Temperatura projetada se a conveccao nao ocorrer
           TN(I,K) = T(I,K) + ( cpdTdt/cp )*dtime         
           !Umidade projetada se a conveccao nao ocorrer
           QO(I,K) = Q(I,K) +   rtt(kr,i,j)*dtime
           !Atribuicoes DO eschema
           P(I,K)  = PO(I,K)
           !srf PSUR(I) = 0.5*(PO(I,1)+PO(I,2))
           if(TN(I,K).lt.200.)    TN(I,K) = T(I,K)
           if(QO(I,K).lt.1.E-08)  QO(I,K) = 1.E-08
           OUTT(I,K) = 0.   !Tendencia no campo de temperatura associada aos cumulus
           OUTQ(I,K) = 0.   !Tendencia na razao de mist. de vapor d'agua assoc. aos cumulus
        enddo
     enddo
     !---  CUMULUS PARAMETERIZATION
           !iens = 2              !--- shallow convection
!     write (*,*) 'Vou entrar na sub-rotina cup_enss_shal...'
     call CUP_enss_shal(mynum,m1,m2,m3,i0,                                   & ! 06
                        j0,mgmxp,mgmyp,mgmzp,maxiens,maxens,                 & ! 12
                        maxens2,maxens3,ensdim,icoic,iupmethod,depth_min,       & ! 18
                        cap_maxs,radius,zkbmax,zcutdown,z_detr,j,iens_tmp, & !
                        ISTART,IEND,mix,                 & ! 24
                        mjx,mkx,xland,z1,AA0,T,                              & ! 30
                        Q,TN,QO,PO,P,OUTT,                                   & ! 36
                        OUTQ,DTIME,PSUR,T00,time,OMEG,                     & ! 42
                        PBLIDX,TKEG,RCPG,tkmin,                              & ! 46
                        ierr4d(1,j,iens),       jmin4d(1,j,iens),      & ! 48
                        kdet4d(1,j,iens),       k224d(1,j,iens),       & ! 50
                        kbcon4d(1,j,iens),      ktop4d(1,j,iens),      & ! 52
                        kpbl4d(1,j,iens),       kstabi4d(1,j,iens),    & ! 54
                        kstabm4d(1,j,iens),     xmb4d(1,j,iens),       & ! 56
                        edt4d(1,j,iens),        zcup5d(1,1,j,iens), & ! 58
                        pcup5d(1,1,j,iens),  enup5d(1,1,j,iens), & ! 60
                        endn5d(1,1,j,iens),  deup5d(1,1,j,iens), & ! 62
                        dedn5d(1,1,j,iens),  zup5d(1,1,j,iens),  & ! 64
                        zdn5d(1,1,j,iens),   prup5d(1,1,j,iens), & ! 66
                        clwup5d(1,1,j,iens), tup5d(1,1,j,iens),  & ! 68
                        upmf,xierr,xktop,xkbcon,xk22,xierr_dp)!,sgrell         ! 74

     !------------------------- Output   ----------------------------------
     do K=1,MKX-1 
        kr = K + 1
        do I = ISTART,IEND
           !      Converte tendencia da temperatura (OUTT) em tendencia de theta (OUTTEM)
           !      cp*T=Pi*Theta => cp dT/dt = Theta*dPi/dt + Pi*dTheta/dt,
           !      assuminDO dPi/dt (=pt(kr,i,j)) << (exner/theta)*dTheta/dt:
           !      Exner's function = pp(kr,i,j)+pi0(kr,i,j)
           exner          =   pp(kr,i,j) + pi0(kr,i,j)
           !
           outtem(kr,i,j) =   CP/exner   *  OUTT(I,K) ! tendencia DO Theta  devida aos cumulus
           outrt(kr,i,j) =   OUTQ(I,K)               ! tendencia DO Rtotal devida aos cumulus
           !
        enddo
     enddo
  enddo     ! loop externo - j -
end subroutine CUPARTH_shal
!===========================================================================================!
!===========================================================================================!






!===========================================================================================!
!===========================================================================================!
subroutine CUP_enss_shal(mynum,m1,m2,m3,i0,                                & !
                         j0,mgmxp,mgmyp,mgmzp,maxiens,maxens,              & !
                         maxens2,maxens3,ensdim,icoic,iupmethod,depth_min, & !
                         cap_maxs,radius,zkbmax,zcutdown,z_detr,           & !
                         j,iens,ISTART,IEND,mix,                           & !
                         mjx,mkx,xland,Z1,AAEQ,T,                          & !
                         Q,TN,QO,PO,P,OUTT,                                & !
                         OUTQ,DTIME,PSUR,TCRIT,time,omeg,                  & !
                         PBLIDX,TKEG,RCPG,tkmin,                           & !
                         ierr4d,jmin4d,kdet4d,k224d,kbcon4d,ktop4d,        & !
                         kpbl4d,kstabi4d,kstabm4d,xmb4d,edt4d,zcup5d,      & !
                         pcup5d,enup5d,endn5d,deup5d,dedn5d,zup5d,         & !
                         zdn5d,prup5d,clwup5d,tup5d,                       & !
                         upmf,xierr,xktop,xkbcon,xk22,xierr_dp)            ! !
     

  use mem_scratch3_grell_sh
  use rconstants, only: g,cpi,alvl

  implicit none
  integer mynum,i0,j0,m1,m2,m3
  integer maxiens,maxens,maxens2,maxens3,ensdim
  integer mix,mjx,mkx,mgmxp, mgmyp, mgmzp
  integer nens,iens,iedt,icoic   !,nall,nens3,izero,ktau
  real(kind=8) :: time
  real :: tkmin
  !--- Input variables -----------------------------
  real    mconv(mgmxp),Z1(mgmxp),AAEQ(mgmxp),PSUR(mgmxp)  !,pre(mgmxp),direction(mgmxp)
  real    T(mgmxp,mgmzp),Q(mgmxp,mgmzp),TN(mgmxp,mgmzp),QO(mgmxp,mgmzp), &
       P(mgmxp,mgmzp),PO(mgmxp,mgmzp),omeg(mgmxp,mgmzp),TKEG(mgmxp,mgmzp), &
       RCPG(mgmxp,mgmzp)

  integer, dimension(mgmxp) :: PBLIDX
!  real, dimension(m1,m2,m3) :: sgrell
  real xland(mgmxp,mgmyp)

  !-------Salva parametros da CUP para uso no transporte convectivo:
  !real  edt_average
  integer,dimension(m2) :: ierr4d,jmin4d,kdet4d,k224d,kbcon4d,ktop4d
  integer,dimension(m2) :: kstabi4d,kstabm4d,kpbl4d
  real,dimension(m2) :: xmb4d,edt4d
  real,dimension(m1,m2) :: zcup5d,pcup5d,enup5d,endn5d,deup5d,dedn5d,&
       zup5d,zdn5d,prup5d,clwup5d,tup5d
  !------Variables saved in RAMS Analisys
  ! use (m1,m2,m3) para dimensionar os vetores que sao escritos nas analises DO RAMS
  real, dimension(m2,m3) :: upmf,xierr,xktop,xkbcon,xk22,xierr_dp
  integer, intent(in) :: iupmethod
  real, intent(in) :: depth_min
  real, intent(in) :: cap_maxs
  real, intent(in) :: radius
  real, intent(in) :: zkbmax
  real, intent(in) :: zcutdown
  real, intent(in) :: z_detr

  integer fquasi  !,fstab,iresult,fmconv
  integer kstart  !,ki,ip1,jp1,m
  integer I,J,K,ISTART,IEND    !,KK
  !--- Output variables ----------------------------
  real OUTT(mgmxp,mgmzp),OUTQ(mgmxp,mgmzp)
  real tcrit,dtime  !,outtes,outteq,xl,pbcdIF,dz
  real cap_inc,dellaqsum,dp,mbdt

  !--- New entrainment/detrainment related stuff --------------------
  real mentr_rate,entr_rate  !,dh,zktop


  cap_inc=0.
  !--- gross entrainment rate (these may be changed later on in the
  !--- program, depending what your detrainment is!!)
  entr_rate=.2/radius
  !--- entrainment of mass
  mentr_rate =entr_rate
  !--- initial detrainmentrates
  ! strong lateral mixing, entrtainment=dedtrainment,
  ! constant mass flux with height
  do k=1,mkx
     do i=istart,iend
        cd(i,k)=entr_rate
     enddo
  enddo

  do I=ISTART,IEND
     cap_max(i)=cap_maxs
     aa0(i)=0.
     aa1(i)=0.

     ! MLO - May 2008 - Commenting to test the "switch"
     !srf -fev 2003
     !don't permite shallow IF there is deep convection
     !if(int(xierr_dp(i,j)) == 0 ) then 
     !   ierr(i)  =20
     !else
        ierr(i)  =0
        xierr(i,j)=0.
     !endif
  enddo

  do nens=1,maxens
     mbdt_ens(nens)=(float(nens)-3.)*dtime*1.e-3+dtime*5.E-03
  enddo

  !--- environmental conditions, FIRST HEIGHTS
  do i=istart,iend
     if(ierr(i) /= 20)then
        do k=1,maxens*maxens2*maxens3
           xf_ens(  i,j,(iens-1)*maxens*maxens2*maxens3+k)= 0.
        enddo
     endif
  enddo

  !--- calculate moist static energy, heights, qes
  call cup_env(j,z,qes,he,hes,t,q,p,z1,mix,mgmxp, &
       mkx,mgmzp,istart,iend,psur,ierr,tcrit,0)
  call cup_env(j,zo,qeso,heo,heso,tn,qo,po,z1,mix,mgmxp, &
       mkx,mgmzp,istart,iend,psur,ierr,tcrit,0)

  !--- environmental values on cloud levels
  call cup_env_clev(j,t,qes,q,he,hes,z,p,qes_cup, &
       q_cup,he_cup,hes_cup,z_cup,p_cup,gamma_cup,t_cup,psur, &
       mix,mgmxp,mkx,mgmzp,istart,iend,ierr,z1)
  call cup_env_clev(j,tn,qeso,qo,heo,heso,zo,po, &
       qeso_cup,qo_cup,heo_cup,heso_cup,zo_cup,po_cup, &
       gammao_cup,tn_cup,psur,mix,mgmxp,mkx,mgmzp,istart, &
       iend,ierr,z1)

  do i=istart,iend
     if(ierr(i) == 0)then
        do k=1,mkx
           if(zo_cup(i,k).gt.zkbmax+z1(i))then
              kbmax(i)=k
              ! PRINT*,'kbmax', zkbmax,z1(i),zo_cup(i,k),kbmax(i)
              exit
           endif
        enddo
     endif
  enddo

  !--- DETERMINE LEVEL WITH HIGHEST MOIST STATIC ENERGY CONTENT - K22
  !
  !Gr-dec2002
  !   here, ,should just take top of PBL for shallow clouds, also maybe for
  !   deep clouds in tropics: k22=level(pbltop)
  !
  !srf-fev2003

  if (iupmethod == 2 .and. pblidx(i) == 0) then
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
  elseif (iupmethod == 2 .and. pblidx(i) /= 0) then
     k22(i) = pblidx(i)
  else !(iupmethod == 1)
     !	 Old way to define k22
     kstart=2
     call maximi(heo_cup,mix,mgmxp,mkx,mgmzp,kstart,kbmax,k22,istart,iend,ierr)
  endif

  do I=ISTART,IEND
     if (ierr(I).eq.0) then
        if(K22(I) >= KBMAX(i))ierr(i)=2
     endif
  enddo

  !--- DETERMINE THE LEVEL OF CONVECTIVE CLOUD BASE  - KBCON
  call cup_kbcon(1, k22, kbcon, heo_cup, heso_cup, mix, mgmxp,  &
       mkx, mgmzp, istart, iend, ierr, kbmax, po_cup, cap_max)

  !--- Calculate incloud moist static energy
  call cup_up_he(k22,hkb,z_cup,cd,mentr_rate,he_cup,hc,mix,mgmxp,mkx, &
       mgmzp,kbcon,ierr,istart,iend,dby,he,hes_cup)
  call cup_up_he(k22,hkbo,zo_cup,cd,mentr_rate,heo_cup,hco, &
       mix,mgmxp,mkx,mgmzp,kbcon,ierr,istart,iend,dbyo, &
       heo,heso_cup)


  call cup_ktop(1,dbyo,kbcon,ktop,mix,mgmxp,mkx,mgmzp,istart,iend,ierr)

  !--- Normalized updraft mass flux profile
  call cup_up_nms(zu,z_cup,mentr_rate,cd,kbcon,ktop,mix,mgmxp, &
       mkx,mgmzp,istart,iend,ierr,k22)
  call cup_up_nms(zuo,zo_cup,mentr_rate,cd,kbcon,ktop,mix,mgmxp, &
       mkx,mgmzp,istart,iend,ierr,k22)

  !--- Calculate moisture properties of updraft
  call cup_up_moisture(ierr,z_cup, qc, qrc, pw, pwav,kbcon,ktop,mix,&
       mgmxp,mkx,mgmzp,istart,iend,cd,dby,mentr_rate, &
       q,GAMMA_cup, zu, qes_cup, k22,q_cup)
  call cup_up_moisture(ierr,zo_cup,qco,qrco,pwo,pwavo,kbcon,ktop,mix, &
       mgmxp,mkx,mgmzp,istart,iend,cd,dbyo,mentr_rate, &
       qo,GAMMAo_cup,zuo,qeso_cup,k22,qo_cup)

  !--- Calculate workfunctions for updrafts
  call cup_up_aa0(aa0,z, zu, dby, GAMMA_CUP, t_cup,kbcon,ktop,mix, &
       mgmxp,mkx,mgmzp,istart,iend,ierr)
  call cup_up_aa0(aa1,zo,zuo,dbyo,GAMMAo_CUP,tn_cup,kbcon,ktop,mix,&
       mgmxp,mkx,mgmzp,istart,iend,ierr)

  do i=istart,iend
     if(ierr(i) == 0 .and. aa1(i).eq.0.) then
       ierr(i)=17
     endif
  enddo

  iedt=1
  do k=1,mkx
     do i=istart,iend
        dellat_ens(i,k,iedt)=0.
        dellaq_ens(i,k,iedt)=0.
     enddo
  enddo


  !--- changes per unit mass from shallow convection
  call cup_dellas_shallow(ierr,zo_cup,po_cup,heo,mix,mgmxp,mkx, &
       mgmzp,istart,iend,dellah,1,j,zuo,cd,hco, &
       ktop,k22,kbcon,mentr_rate,heo_cup, &
       'shallow')
  call cup_dellas_shallow(ierr,zo_cup,po_cup,qo ,mix,mgmxp,mkx,mgmzp,&
       istart,iend,dellaq,2,j,zuo, cd,qco,ktop,k22, &
       kbcon,mentr_rate,qo_cup,'shallow')

  !--- Using dellas, calculate changed environmental profiles
  do nens=1,maxens
     mbdt=mbdt_ens(nens)
     do i=istart,iend
        xaa0_ens(i,nens)=0.
     enddo
     do k=1,mkx-1
        do i=istart,iend
           dellat(i,k)=0.
           if(ierr(i).eq.0)then
              XHE(I,K) =      DELLAH(I,K)*MBDT + HEO(I,K)
              XQ(I,K) =	      DELLAQ(I,K)*MBDT +  QO(I,K)
              DELLAT(I,K) =  cpi*(DELLAH(I,K)-alvl*DELLAQ(I,K))
              XT(I,K) =	      DELLAT(I,K)*MBDT +  TN(I,K)
              if(XQ(I,K).le.0.)XQ(I,K)=1.E-08

           endif
        enddo
     enddo

     do i=istart,iend
        if(ierr(i).eq.0)then
           XHE(I,mkx)  = HEO(I,mkx)
           XQ(I,mkx)  =  QO(I,mkx)
           XT(I,mkx)  =  TN(I,mkx)
           if(XQ(I,mkx).le.0.) XQ(I,mkx)  =  1.E-08
        endif
     enddo

     !--- Calculate moist static energy, heights, qes
     call cup_env(j,xz,xqes,xhe,xhes,xt,xq,po,z1,mix,mgmxp, &
          mkx,mgmzp,istart,iend,psur,ierr,tcrit,2)

     !--- Environmental values on cloud levels

     call cup_env_clev(j,xt,xqes,xq,xhe,xhes,xz,po,xqes_cup, &
          xq_cup,xhe_cup,xhes_cup,xz_cup,po_cup,gamma_cup,xt_cup, &
          psur,mix,mgmxp,mkx,mgmzp,istart,iend,ierr,z1)
     !**************************** Static Control ************
     !--- Moist static energy inside cloud
     do i=istart,iend
        if(ierr(i) == 0)then
           xhkb(i)=xhe(i,k22(i))
        endif
     enddo

     call cup_up_he(k22,xhkb,xz_cup,cd,mentr_rate,xhe_cup,xhc,&
          mix,mgmxp,mkx,mgmzp,kbcon,ierr,istart,iend,xdby,xhe,&
          xhes_cup)

     !--- Normalized mass flux profile
     call cup_up_nms(xzu,xz_cup,mentr_rate,cd,kbcon,ktop, mix,mgmxp,&
          mkx,mgmzp,istart,iend,ierr,k22)

     !--- Moisture updraft
     call cup_up_moisture(ierr,xz_cup,xqc,xqrc,xpw,xpwav,&
          kbcon,ktop,mix,mgmxp,mkx,mgmzp,istart,iend,cd,xdby,&
          mentr_rate,xq,GAMMA_cup,xzu,xqes_cup,k22,xq_cup)

     !--- Workfunctions for updraft
     call cup_up_aa0(xaa0,xz,xzu,xdby,GAMMA_CUP,xt_cup,kbcon, &
          ktop,mix,mgmxp,mkx,mgmzp,istart,iend,ierr)

     !srf-feb-2003
     do i=istart,iEND
        if(ierr(i).eq.0) xaa0_ens(i,nens)=xaa0(i)
     enddo

  enddo

  !--------- LARGE SCALE FORCING  -----------------------------------
  call cup_forcing_ens_shal(aa0,aa1,xaa0_ens,mbdt_ens,dtime,xmb, &
       ierr,mix,mgmxp,mjx,mgmyp,mkx,mgmzp,istart,iend,xf_ens,j, &
       'shallow',xland,maxens,iens,maxens2,maxens3, &
       ensdim,p_cup,ktop,icoic,iedt)

  do k=1,mkx
     do i=istart,iend
        if(ierr(i) == 0)then
           dellat_ens(i,k,iedt) =  dellat(i,k)
           dellaq_ens(i,k,iedt) =  dellaq(i,k)
        else 
           dellat_ens(i,k,iedt) = 0.
           dellaq_ens(i,k,iedt) = 0.
        endif
     enddo
  enddo
  !250  continue

  !------------------------------ FEEDBACK ------------------------------------
  call cup_output_ens_shal(xf_ens,ierr,dellat_ens,dellaq_ens, &
       outt,outq,xmb,ktop,mix,mgmxp,mjx,mgmyp,mkx,mgmzp,istart, &
       iend,j,'shallow',maxens2,maxens,iens, &
       maxens3,ensdim,xfac1,m1,m2,m3)!,sgrell

  !---------------------------done shallow cumulus scheme -------------

  !-------Salva parametros nas analises DO RAMS
  do i=istart,iend
     xierr(i,j)=float(ierr(i))
     if(ierr(i) == 0)then
        ! Shallow updraft mass flux
        upmf(i,j)= xmb(i)
        xktop(i,j)= float(ktop(i))
        xkbcon(i,j)= float(kbcon(i))
        xk22(i,j)= float(k22(i))
     elseif(ierr(i) /= 20)then
        upmf(i,j) = 0.
        xktop(i,j) = 0.
        xkbcon(i,j) = 0.
        xk22(i,j) = 0.
     endif
  enddo

  do i=istart,iend
     ierr4d(i) =  ierr(i)
     jmin4d(i) =  1 ! no DOwndrafts
     kdet4d(i) =  1 ! no DOwndrafts
     k224d(i) =   k22(i)
     kbcon4d(i) = kbcon(i)
     ktop4d(i) =  ktop(i)
     kstabi4d(i) =kstabi(i)
     kstabm4d(i) =kstabm(i)
     !	 kpbl4d(i) =  kpbl(i)
     kpbl4d(i) =  k22(i) 
     xmb4d(i) =  xmb(i)
     edt4d(i)=0.  ! no DOwndrafts
     upmf(i,j) = 0.
     if(ierr(i) == 0) upmf(i,j) = xmb(i)
     !updraft mass flux averaged
     do k=1,mkx
        !srf- neste ponto iens=iens_tmp = 1
        if(iens == 1) then
           zcup5d(k,i) = zo_cup(i,k)
           pcup5d(k,i) = po_cup(i,k)       
        endif
        enup5d(k,i) = mentr_rate
        endn5d(k,i) = 0.          ! no DOwndrafts
        deup5d(k,i) =  cd(i,k)
        dedn5d(k,i) = 0.          ! no DOwndrafts
        zup5d(k,i) = zuo(i,k)
        zdn5d(k,i) = 0.          ! no DOwndrafts 
        !?? It's save for use at wet-deposition scheme (in-cloud)
        prup5d(k,i) = 0.          ! no precip for shallow
        clwup5d(k,i) = qrco(i,k)   ! cloud liq  water  - only for upfradt
        tup5d(k,i) = t_cup(i,k)   !>>> em verdade deveria ser a temperatura da parcela

     enddo
  enddo
end subroutine CUP_enss_shal

!--------------------------------------------------------------------

subroutine cup_dellas_shallow(ierr,z_cup,p_cup,he,mix,mgmxp,mkx,mgmzp, &
     istart,iend,della,itest,j,zu,cd,hc,ktop,k22,kbcon,mentr_rate, &
     he_cup,name)
  use rconstants, only: g
  implicit none
  character *(*) name
  integer mix,mgmxp,mkx,mgmzp,i,k,istart,iend,itest,j
  real z_cup(mgmxp,mgmzp), p_cup(mgmxp,mgmzp),he(mgmxp,mgmzp)
  real della(mgmxp,mgmzp), hc(mgmxp,mgmzp), cd(mgmxp,mgmzp)
  real zu(mgmxp,mgmzp),he_cup(mgmxp,mgmzp)
  integer kbcon(mgmxp),ktop(mgmxp),k22(mgmxp),ierr(mgmxp)
  real dp,dz,mentr_rate,subin,entup,detup,subDOwn,entupk,detupk,totmas
  !real xsum,xsumt

  do K=2,MKX
     do I=ISTART,IEND
        della(i,k)=0.
     enddo
  enddo
  !       xsum=0.
  !       xsumt=0.
  !
  do k=2,mkx-1
     do i=istart,iend
        if (ierr(i).ne.0) cycle
        if (k.gt.ktop(i)) cycle
        !
        !--- SpecIFy detrainment of DOwndraft, has to be consistent
        !--- with zd calculations in soundd.
        !
        dz    = Z_cup(i,k+1)-Z_cup(i,k)
        subin = zu(i,k+1)
        entup = 0.
        detup = 0.
        if(k.ge.kbcon(i).and.k.lt.ktop(i))then
           entup = mentr_rate*dz*zu(i,k)
           detup = cd(i,k+1) *dz*zu(i,k)
        endif
        subDOwn = zu(i,k)
        entupk  = 0.
        detupk  = 0.

        if(k == k22(i)-1)then
           !srf-fev2003  entupk  = zu(i,kpbl(i))
           entupk  = zu(i,k22(i))
        endif

        if(k == ktop(i))then
           detupk  = zu(i,ktop(i))
           subin   = 0.
        endif

        if(k < kbcon(i))then
           detup   = 0.
        endif
        !
        !--- Changed due to subsidence and entrainment
        !
        totmas=subin-subDOwn+detup-entup-entupk+detupk
        if(abs(totmas).gt.1.e-6)then
           print *,'**TOTMAS*******',i,j,k,totmas,name
           print *,k22(i),kbcon(i),ktop(i)
           print *,'updr stuff = ',subin,subDOwn,detup,entup,entupk,detupk
           stop
        endif

        dp =  100.*( p_cup(i,k)-p_cup(i,k+1) )
        della(i,k)=(subin*he_cup(i,k+1)-subDOwn*he_cup(i,k)+ &
             detup*.5*( HC(i,K+1)+ HC(i,K))-entup*he(i,k)-&
             entupk*he_cup(i,k22(i))+detupk*hc(i,ktop(i)))*g/dp
     enddo
  enddo
end subroutine cup_dellas_shallow

!------------------------------------------------------------

subroutine cup_forcing_ens_shal(aa0,aa1,xaa0,mbdt,dtime,xmb,ierr, &
     mix,mgmxp,mjx,mgmyp,mkx,mgmzp,istart,iend,xf,j, &
     name,xland,maxens,iens,maxens2,maxens3, &
     ensdim,p_cup,ktop,icoic,iedt)

  implicit none
  character *(*) name

  integer k,i,istart,iend,mix,mgmxp,mjx,mgmyp,mkx,mgmzp,j
  integer maxens,maxens3,ensdim,iens,nall,maxens2

  !------ ensemble 3 dimension = 10
  integer kclim
  !!srf- new interpolated PARAMETER:
  integer,parameter :: mkxcrt=25
  real pcrit(mkxcrt),acrit(mkxcrt),acritt(mkxcrt),aclim1
  real aclim2,aclim3,aclim4
  data pcrit/ 850., 837.5, 825., 812.5, 800., 787.5, &
       775., 762.5, 750., 737.5, 725., 712.5, &
       700., 687.5, 675., 662.5, 650., 637.5, &
       625., 612.5, 600., 550. , 500., 450., &
       400./
  data acrit/ 6.323E-02,  5.795E-02, 5.390E-02,  5.236E-02, &
       4.450E-02,  4.965E-02, 5.000E-02,  4.983E-02, &
       5.530E-02,  5.289E-02, 6.080E-02,  5.883E-02, &
       6.640E-02,  6.766E-02, 7.070E-02,  7.937E-02, &
       7.500E-02,  9.396E-02, 0.108,      0.111, &
       0.130,      0.152,     0.221,      0.315, &
       0.368/
  data acritt/0.203,  0.299, 0.359,  0.403, 0.515,  0.478, &
       0.518,  0.530, 0.521,  0.565, 0.543,  0.588, &
       0.566,  0.602, 0.596,  0.611, 0.625,  0.619, &
       0.645,  0.627, 0.665,  0.659, 0.688,  0.743, &
       0.813/

  integer ktop(mgmxp),iedt
  real p_cup(mgmxp,mgmzp)
  integer ierr(mgmxp)  
  real aa0(mgmxp),aa1(mgmxp),xaa0(mgmxp,maxens),xmb(mgmxp)
  real mbdt(maxens),dtime
  real     xf(mgmxp,mgmyp,ensdim),xland(mgmxp,mgmyp)
  real xff_ens3(maxens3),xk(maxens),xff0  !,xff1,xff2,xff3,xff
  integer nens,ne,n,iresultd,iresulte,icoic  !,nens3,iresult

  nens=0

  !--- LARGE SCALE FORCING
  !
  do I=ISTART,IEND
     xmb(i)=0.
     if(ierr(i).eq.0)then
        kclim=0
        do k=mkxcrt,1,-1
           if(p_cup(i,ktop(i)).lt.pcrit(k))then
              kclim=k
              GO TO 9
           endif
        enddo
        if(p_cup(i,ktop(i)).gt.pcrit(1))kclim=1
9       continue
        kclim=min(max(kclim,1),mkxcrt)
        k= max(kclim-1,1)
        aclim1= acrit(kclim)*1.e3
        aclim2= acrit(k)*1.e3
        aclim3= acritt(kclim)*1.e3
        aclim4= acritt(k)*1.e3
        !
        !---- Grell's closure
        !
        xff0	   =  (AA1(I)-AA0(I))/dtime
        xff_ens3(1)=  (AA1(I)-AA0(I))/dtime
        xff_ens3(2)= .9*xff_ens3(1)
        xff_ens3(3)=1.1*xff_ens3(1)

        !--- More original Arakawa-Schubert 
        xff_ens3(4)=max(0.,(AA1(I)-aclim1)/dtime)
        xff_ens3(5)=max(0.,(AA1(I)-aclim2)/dtime)
        xff_ens3(6)=max(0.,(AA1(I)-aclim3)/dtime)
        xff_ens3(7)=max(0.,(AA1(I)-aclim4)/dtime)

        !--- More like Fritsch Chappel or Kain Fritsch (plus triggers)
        xff_ens3(8) =AA1(I)/(60.*20.)
        xff_ens3(9) =AA1(I)/(60.*30.)
        xff_ens3(10)=AA1(I)/(60.*40.)

        do nens=1,maxens
           XK(nens)=(XAA0(I,nens)-AA1(I))/MBDT(nens)
           if(xk(nens).le.0.and.xk(nens).gt.-1.e-9) xk(nens)=-1.e-9
           if(xk(nens).gt.0.and.xk(nens).lt.+1.e-9) xk(nens)=+1.e-9
        enddo

        !--- Add up all ensembles
        do ne=1,maxens
           !
           !--- for every xk, we have maxens3 xffs
           !--- iens is from outermost ensemble (most expensive!
           !
           !--- iedt (maxens2 belongs to it)
           !--- is from second, next outermost, not so expensive
           !
           !--- so, for every outermost loop, we have maxens*maxens2*3
           !--- ensembles!!! nall would be 0, IF everything is on first
           !--- loop index, THEN ne would start counting, THEN iedt, THEN iens....
           !
           iresultd=0
           iresulte=0

           nall= (iens-1)*maxens3*maxens*maxens2+(iedt-1)*maxens3*maxens+ &
                (ne  -1)*maxens3
           !
           !--- Special treatment for stability closures
           !


           if(xff0.gt.0.and.xk(ne).lt.0.)then
              xf(i,j,nall+1) =max(0.,-xff_ens3(1) /xk(ne))
              xf(i,j,nall+2) =max(0.,-xff_ens3(2) /xk(ne))
              xf(i,j,nall+3) =max(0.,-xff_ens3(3) /xk(ne))
              !			PRINT*,'XFSH',xf(i,j,nall+1),-xff_ens3(1),xk(ne)
           endif
           !240 Subscript error  array=xf size=243810 subscript=243944 eln=1054       
           if(XK(ne).lt.0.)then
              xf(i,j,nall+4) =max(0.,-xff_ens3(4) /xk(ne))
              xf(i,j,nall+5) =max(0.,-xff_ens3(5) /xk(ne))
              xf(i,j,nall+6) =max(0.,-xff_ens3(6) /xk(ne))
              xf(i,j,nall+7) =max(0.,-xff_ens3(7) /xk(ne))
              xf(i,j,nall+8) =max(0.,-xff_ens3(8) /xk(ne))
              xf(i,j,nall+9) =max(0.,-xff_ens3(9) /xk(ne))
              xf(i,j,nall+10)=max(0.,-xff_ens3(10)/xk(ne))
           endif
           if(icoic.ge.1)then
              xf(i,j,nall+1) =xf(i,j,nall+icoic)
              xf(i,j,nall+2) =xf(i,j,nall+icoic)
              xf(i,j,nall+3) =xf(i,j,nall+icoic)
              xf(i,j,nall+4) =xf(i,j,nall+icoic)
              xf(i,j,nall+5) =xf(i,j,nall+icoic)
              xf(i,j,nall+6) =xf(i,j,nall+icoic)
              xf(i,j,nall+7) =xf(i,j,nall+icoic)
              xf(i,j,nall+8) =xf(i,j,nall+icoic)
              xf(i,j,nall+9) =xf(i,j,nall+icoic)
              xf(i,j,nall+10)=xf(i,j,nall+icoic)
           endif
        end do
        cycle
     elseif(ierr(i).ne.20.and.ierr(i).ne.0)then
        do n=1,ensdim
           xf(i,j,n)=0.
        enddo
     endif
  end do

end subroutine cup_forcing_ens_shal

!--------------------------------------------------------------------

subroutine cup_output_ens_shal(xf,ierr,dellat,dellaq,outt,outq,xmb, &
     ktop,mix,mgmxp,mjx,mgmyp,mkx,mgmzp, &
     istart,iend,j,name,maxens2,maxens, &
     iens,maxens3,ensdim,xfac1, &
     m1,m2,m3)!,sgrell
  use rconstants, only : day_sec
  implicit none
  character *(*) name
  integer mix,mjx,mkx,istart,iend,mgmxp,mgmyp,mgmzp,ensdim,i,k,j,n
  integer m1,m2,m3,maxens,maxens2
  real	    xf(mgmxp,mgmyp,ensdim)
  real outt(mgmxp,mgmzp),outq(mgmxp,mgmzp),dellat(mgmxp,mgmzp,maxens2)
  real dellaq(mgmxp,mgmzp,maxens2),xmb(mgmxp),xfac1(mgmxp)
  integer ktop(mgmxp),ierr(mgmxp),ncount,iens,maxens3  !,nens3
  real outtes,ddtes,dtt,dtq

  ! Statistics properties
  integer i_use_stat_prop
  data  i_use_stat_prop/0/

  do K=1,MKX
     do I=ISTART,IEND
        outt(i,k)=0.
        outq(i,k)=0.
     enddo
  enddo
  do I=ISTART,IEND
     xmb(i)=0.
     xfac1(i)=1.
  enddo

     !----  Simple average 
     do I=ISTART,IEND
        ncount=0
        xmb(i)=0.
        if(ierr(i).eq.0)then
           do n=(iens-1)*maxens2*maxens*maxens3+1, &
                iens*maxens2*maxens*maxens3
              if(xf(i,j,n).gt.0.)then
                 xmb(i)=xmb(i)+xf(i,j,n)
                 ncount=ncount+1
              endif
           enddo
           if(ncount.gt.0)then
              xmb(i)=xmb(i)/float(ncount)
           else
              xmb(i)=0.
              ierr(i)=13
           endif
        endif
     enddo

  !-- now DO feedback
  ddtes=250.
  if(name.eq.'shallow')ddtes=500.

  do K=1,MKX
     do I=ISTART,IEND
        dtt=0.
        dtq=0.
        if(ierr(i).eq.0.and.k.le.ktop(i))then
           do n=1,maxens2
              dtt=dtt+dellat(i,k,n)
              dtq=dtq+dellaq(i,k,n)
           enddo
           outtes=dtt*XMB(I)*day_sec/float(maxens2)        
           if(outtes .gt. 2.*ddtes .and. k.gt.2) then
              xmb(i)= 2.*ddtes/outtes * xmb(i)
              outtes=    ddtes
           endif
           if(outtes .lt. -ddtes) then
              XMB(I)= -ddtes/outtes * xmb(i)
              outtes= -ddtes
           endif
           if(outtes .gt. .5*ddtes.and.k.le.2) then
              XMB(I)=	ddtes/outtes * xmb(i)
              outtes=.5*ddtes
           endif
           OUTT(I,K)= XMB(I) *dtt / float(maxens2)
           OUTQ(I,K)= XMB(I) *dtq / float(maxens2)

        endif
     enddo
  enddo
  do i=istart,iend
     if(ierr(i).eq.0)then
        !srf 20fev2003
        xfac1(i) = xmb(i)/(xfac1(i)+1.e-8)
     endif
  enddo
end subroutine cup_output_ens_shal
!---------------------------------------------------------------------------------------------
