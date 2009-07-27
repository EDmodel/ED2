!------------------------------------------------------------------------------------------!
!   Deep convection parameterization, following Grell scheme                               !
!                                                                                          !
!------------------------------------------------------------------------------------------!
subroutine cuparth(mynum,mgmxp,mgmyp,mgmzp,m1,m2,m3,ia,iz,ja,jz,i0,j0,maxiens,iens         &
                  ,iupmethod,depth_min,cap_maxs,radius,zkbmax,zcutdown,z_detr,dtime,time   &
                  ,ua,va,wa,theta,pp,pi0,dn0,rv,kpbl,tke,tkmin,rcp,ht,rtgt,tht,rtt         &
                  ,pt,outtem,outrt,precip,ierr4d,jmin4d,kdet4d,k224d,kbcon4d,ktop4d,kpbl4d &
                  ,kstabi4d,kstabm4d,xmb4d,edt4d,zcup5d,pcup5d,enup5d,endn5d,deup5d,dedn5d &
                  ,zup5d,zdn5d,prup5d,clwup5d,tup5d,upmf,dnmf,xierr,xktop,xkbcon,xk22      &
                  ,xjmin,xkdt,xiact_p,xiact_c)

  ! USE Modules for Grell Parameterization
  use mem_grell_param, only : maxens,  & !INTENT(IN)
       maxens2,                        & !INTENT(IN)
       maxens3,                        & !INTENT(IN)
       ensdim,                         & !INTENT(IN)
       icoic                             !INTENT(IN)

  use mem_scratch2_grell
  use rconstants, only: rgas, cp, rm, p00, t00, g, cpor, pi1,onerad

  !srf   mgmxp, mgmyp, mgmzp sao usadas alocar memoria para as
  !      variaveis da parametrizacao do Grell.

  implicit none
  integer, intent(in) :: mgmxp     ! Number of points in x-direction;
  integer, intent(in) :: mgmyp     ! Number of points in y-direction;
  integer, intent(in) :: mgmzp     ! Number of points in z-direction;
  integer, intent(in) :: iupmethod    ! Method to determine the cloud base; 
  real,    intent(in) :: depth_min ! Minimum depth that the cloud should have [m]
  real,    intent(in) :: cap_maxs  ! Maximum depth of capping inversion [mb]
  real,    intent(in) :: radius
  real,    intent(in) :: zkbmax
  real,    intent(in) :: zcutdown
  real,    intent(in) :: z_detr
  integer :: maxiens               ! Number of different clouds;

  integer :: iens                  ! Flag to denote deep (1) or shallow (2) convection

  !  this is the lowest level, then comes ensemble 1!
  !       ensdim=maxiens*maxens*maxens2*maxens3 !Ensemble dimension
  !INTEGER icoic
  !
  !srf- ICOIC is used for choice a specific closure
  ! icoic = 0 -> ensemble (all closures)  [en]
  ! icoic = 1 -> Grell                    [gr]
  ! icoic = 4 -> low level omega          [lo]
  ! icoic = 7 -> moisture convergence     [mc]
  ! icoic =10 -> like Fritsch Chappel or Kain Fritsch [sc]
  ! icoic =13 -> Arakawa-Schubert         [as]
  !PARAMETER (icoic=1)

  ! Defined in mem_grell_param

!!!!!!!!!!!!


  !------------------------------------ RAMS vectors:
  integer m1,m2,m3,ia,iz,ja,jz,i0,j0,ibcon,j1,j2,mynum
  real(kind=8) :: time
  real         :: dtime
  real         :: tkmin
  real, dimension(m1,m2,m3) :: ua, va, wa, pp, pi0, dn0,   &
       tht, rtt, pt, outtem, outrt,tke,rcp,theta,rv
  real, dimension(m2,m3)    :: ht,rtgt,precip

  real, dimension(m2,m3)    :: upmf, dnmf, xierr, xktop, xkbcon, xjmin,  &
       xkdt, xiact_p, xiact_c, xk22
  integer, dimension(m2,m3) :: kpbl
  !-------Salva parametros da CUP para uso no transporte convectivo:
  integer, dimension(m2,m3,maxiens)    :: ierr4d, jmin4d,  &
       kdet4d, k224d, kbcon4d, ktop4d, kpbl4d, kstabi4d, kstabm4d

  real, dimension(m2,m3,maxiens)       :: xmb4d, edt4d

  real,dimension(m1,m2,m3,maxiens) :: enup5d,  &
       endn5d, deup5d, dedn5d, zup5d, zdn5d
  real,dimension(m1,m2,m3,maxiens) :: pcup5d
  real,dimension(m1,m2,m3,maxiens) :: prup5d  !Lufla
  real,dimension(m1,m2,m3,maxiens) :: clwup5d !Lufla
  real,dimension(m1,m2,m3,maxiens) :: tup5d   !Lufla
  real,dimension(m1,m2,m3,maxiens) :: zcup5d  !Lufla
  !------------------------------------- variaveis locais:

  integer :: kk, istart, iend, i, j, k, mix, mjx, mkx, kr, m
  real    :: vspeed, dp, dq, cpdTdt, exner

  ! Init. arrays P and PO with zeros - ALF

  do k=1,mgmzp
     do i=1,mgmxp
        p(i,k)  = 0.
        po(i,k) = 0.
     enddo
  enddo

  !----------------------------------------------------------------------
  istart = ia
  iend   = iz
  j1     = ja
  j2     = jz
  mkx    = m1 - 1    !mkx nao deve ser igual a m1
  mix    = m2
  mjx    = m3


  do j=1,m3  ! loop em todo dominio para passar informacoes da fronteira
     do i=1,m2 ! dos nodes
        massflx(i,j)     = dnmf(i,j)
        iact_gr(i,j)     = int(xiact_c(i,j))
        iact_old_gr(i,j) = int(xiact_p(i,j))
     enddo
  enddo

  ! Loop externo : j

  do j=j1,j2

     do i = istart,iend
        aa0(i)=0.
        xland(i,j) = 0. ! land/water flag - not in use
     enddo

     do i = istart,iend
        iact_gr(i,j)     = 0    !verificar se isto esta' correto
        iact_old_gr(i,j) = 0    !verificar se isto esta' correto
        kdt(i,j)         = 0
        precip(i,j)      = 0.
        massflx(i,j)     = 0.   !verificar se isto esta' correto
     enddo

     !--- Prepare input, erase output

     do i = istart,iend
        kdet(i)  =2
        pret(i)  =0.
        mconv(i) =0.
        umean(i) =0.
        vmean(i) =0.
        pmean(i) =0.
     enddo

     !------- Transfere valores do RAMS para o eschema
     do k=1,mkx
        kr = k + 1          ! nivel K da grade do Grell corresponde ao
        ! nivel K + 1 do RAMS
        do i = istart,iend

           ter11(i)= ht(i,j)
           psur(i) = .5*( ((pp(1,i,j)+pi0(1,i,j))/cp)**cpor*p00 +  &
                ((pp(2,i,j)+pi0(2,i,j))/cp)**cpor*p00 )*1.e-2
           ! Pressure in mbar

           po(i,k) = ((pp(kr,i,j)+pi0(kr,i,j))/cp)**cpor*p00*1.e-2
           ! Pressure in mbar
           us_grell(i,k) = .5*( ua(kr,i,j) + ua(kr,i-1,j) )
           vs_grell(i,k) = .5*( va(kr,i,j) + va(kr,i,j-1) )
           omeg(i,k)   = -g*dn0(kr,i,j)*.5*( wa(kr,i,j)+wa(kr-1,i,j) )

           t(i,k)  = theta(kr,i,j)*(pp(kr,i,j)+pi0(kr,i,j))/cp
           q(i,k)  = rv(kr,i,j)
           !variables for PBL top height
           pblidx(i) = kpbl(i,j)
           tkeg(i,k) = tke(kr,i,j)
           rcpg(i,k) = rcp(kr,i,j)

           !        Calcula tendencia projetada na temperatura em funcao
           !        das tendencias de theta e PI : cp*T=Pi*Theta
           exner= pp(kr,i,j)+pi0(kr,i,j)

           !        cpdTdt= exner*tht(kr,i,j) + theta(kr,i,j)*pt(kr,i,j)
           cpdTdt  = exner*tht(kr,i,j)
           ! assumindo PT(KR,I,J) << exner*THT(KR,I,J)/theta

           !        Temperatura projetada se a conveccao nao ocorrer
           tn(i,k) = t(i,k) + ( cpdtdt/cp )*dtime

           !        Umidade projetada se a conveccao nao ocorrer
           qo(i,k) = q(i,k) +   rtt(kr,i,j)*dtime


           !------- Atribuicoes do esquema

           p(i,k)  = po(i,k)
           !srf	 PSUR(I) = 0.5*(PO(I,1)+PO(I,2))
           if((psur(i)-p(i,k)).gt.150.and.p(i,k).gt.300.)then
              dp       = -.5*(p(i,k+1)-p(i,k-1))
              umean(i) = umean(i)+us_grell(i,k)*dp
              vmean(i) = vmean(i)+vs_grell(i,k)*dp
              pmean(i) = pmean(i)+dp
           endif

           if(tn(i,k).lt.200.)    tn(i,k) = t(i,k)
           if(qo(i,k).lt.1.e-08)  qo(i,k) = 1.e-08

           outt(i,k) = 0.
           !Tendencia no campo de temperatura associada aos cumulus
           outq(i,k) = 0.
           !Tendencia na razao de mist. de vapor d'agua assoc. aos cumulus
           outqc(i,k) = 0.
           !Tendencia na razao de mistura de agua de nuvem e/ou gelo
           ! associada aos cumulus
        enddo
     enddo

     do i = istart,iend
        umean(i)=umean(i)/pmean(i)
        vmean(i)=vmean(i)/pmean(i)
        vspeed=sqrt(umean(i)*umean(i)+vmean(i)*vmean(i))
        direction(i)=(atan2(umean(i),vmean(i))+pi1)*onerad
        if(direction(i) > 360.)direction(i)=direction(i)-360.
        if(vspeed < 5.)direction(i)=9999.
     enddo

     do k=2,mkx-1
        do i = istart,iend
           dq=.5*(q(i,k+1)-q(i,k-1))
           ! mconv(i)=mconv(i)+1.e5*omeg(i,k)*dq/g
           ! Convergencia de umidade da coluna
           mconv(i)=mconv(i)+omeg(i,k)*dq/g
           ! Convergencia de umidade da coluna (omega em Pa/s)
        enddo
     enddo
     do i = istart,iend
       if(mconv(i).lt. 0.)  mconv(i) = 0.
     enddo

     !---  CUMULUS PARAMETERIZATION
     !srf- aqui se deve colocal o loop no ensemble dependente do tipo
     !     de cumulus
     !iens =1

     call cup_enss(mynum, m1, m2, m3, i0, j0,              &
          mgmxp, mgmyp, mgmzp, maxiens, maxens, maxens2,          &
          maxens3, ensdim, icoic, iupmethod, depth_min, cap_maxs, &
          radius,zkbmax,zcutdown,z_detr,  j, iens, istart, iend,  &
          mix, mjx, mkx, massfln, massflx, iact_gr, iact_old_gr,  &
          xland, ter11, aa0, t, q, tn, qo, po, pret, p, outt,     &
          outq, outqc, dtime, psur, us_grell, vs_grell, kdet, t00, time,    &
          mconv, omeg, pblidx,tkeg,rcpg,tkmin, direction,        &
          !Insercao de variaveis para salvar  propriedades dos cumulus
          ierr4d(1,j,iens), jmin4d(1,j,iens),        &
          kdet4d(1,j,iens), k224d(1,j,iens),         &
          kbcon4d(1,j,iens), ktop4d(1,j,iens),       &
          kpbl4d(1,j,iens),                                &
          kstabi4d(1,j,iens), kstabm4d(1,j,iens),    &
          xmb4d(1,j,iens), edt4d(1,j,iens),          &
          zcup5d(1,1,j,iens),     &
	  pcup5d(1,1,j,iens),     &
          enup5d(1,1,j,iens), endn5d(1,1,j,iens),    &
          deup5d(1,1,j,iens), dedn5d(1,1,j,iens),    &
          zup5d(1,1,j,iens), zdn5d(1,1,j,iens),      &
          prup5d(1,1,j,iens),     &
	  clwup5d(1,1,j,iens),    &
          tup5d(1,1,j,iens),      &
          !Insercao de variaveis para salvar  nas analises e  Cup_DIRECTION2
          upmf, dnmf, xierr, xktop, xkbcon, xk22, xjmin, xkdt,         &
          xiact_p, xiact_c)


     !--- Output

     do k=1,mkx-1
        kr = k + 1
        do i = istart,iend
           !      Converte tendencia da temperatura (OUTT) em
           !      tendencia de theta (OUTTEM)
           !      cp*T=Pi*Theta => cp dT/dt = Theta*dPi/dt + Pi*dTheta/dt,
           !      assumindo dPi/dt (=pt(kr,i,j)) << (exner/theta)*dTheta/dt:
           !      Exner's function = pp(kr,i,j)+pi0(kr,i,j)
           exner          = pp(kr,i,j) + pi0(kr,i,j)
           outtem(kr,i,j) = cp/exner   * outt(i,k)
           ! tendencia do Theta  devida aos cumulus
           outrt(kr,i,j)  = outq(i,k)  + outqc(i,k)
           ! tendencia do Rtotal devida aos cumulus
        enddo
     enddo


     do i = istart,iend
        precip(i,j)=pret(i)
     enddo

     do i = istart,iend
        if(precip(i,j).le.0.)then
           iact_gr(i,j) =0
           precip(i,j)  =0.
           do k=1,mkx
              kr = k + 1
              outtem(kr,i,j) = 0.
              ! comente para testes com conservacao (c0=0.)
              outrt(kr,i,j) =0.
              ! comente para testes com conservacao (c0=0.)
           enddo
           do k=1,ensdim
              massfln(i,j,k)=0.
           enddo
        else
           iact_gr(i,j)=1
        endif
     enddo

     !--- Salva nas analises

     do i = istart,iend
        massflx(i,j) = dnmf(i,j)
        xiact_c(i,j) = float(iact_gr(i,j))
        xiact_p(i,j) = float(iact_old_gr(i,j))
     enddo

  enddo     ! loop externo - j -

  !--- shallow convection

  return
end subroutine cuparth

!--------------------------------------------------------------------
!
!

subroutine cup_enss(mynum, m1, m2, m3, i0, j0,                 &
     mgmxp, mgmyp, mgmzp, maxiens, maxens, maxens2, maxens3, ensdim,  &
     icoic, iupmethod,depth_min,cap_maxs,radius,zkbmax,zcutdown,z_detr &
     , j, iens, istart, iend, mix, mjx, mkx, massfln, massflx,   &
     iact_gr, iact_old_gr, xland, z1, aaeq, t, q, tn, qo, po, pre, p, &
     outt, outq, outqc, dtime, psur, us, vs, kdet, tcrit, time,       &
     mconv, omeg, pblidx,tkeg,rcpg,tkmin, &
     direction, ierr4d, jmin4d, kdet4d, k224d, kbcon4d,  &
     ktop4d, kpbl4d, kstabi4d, kstabm4d, xmb4d, edt4d,                &
     zcup5d,        &
     pcup5d,        &
     enup5d, endn5d, deup5d, dedn5d, zup5d, zdn5d,            &
     prup5d,clwup5d,tup5d, & !Lufla    &
     upmf, dnmf, xierr, xktop, xkbcon, xk22, xjmin, xkdt, xiact_p, xiact_c)

  ! USE Modules for Grell Parameterization
  use mem_scratch3_grell
  use rconstants, only: g, day_sec,alvl,cp
  
  implicit none
  integer maxiens,maxens,maxens2,maxens3,ensdim
  integer mix,mjx,mkx,mgmxp, mgmyp, mgmzp
  integer nall,nens,iens,iedt,ktau
  integer nens3
  integer izero
  integer icoic
  integer,intent(in) :: iupmethod
  real, intent(in) :: depth_min
  real, intent(in) :: cap_maxs
  real, intent(in) :: radius
  real, intent(in) :: zkbmax
  real, intent(in) :: zcutdown
  real, intent(in) :: z_detr





  integer, save :: ialloc
  data ialloc/0/
  real(kind=8) :: time

  !--- Input variables -----------------------------
  real tkmin
  real mconv(mgmxp), z1(mgmxp), direction(mgmxp), aaeq(mgmxp),  &
       pre(mgmxp), psur(mgmxp)
  real t(mgmxp,mgmzp), q(mgmxp,mgmzp), tn(mgmxp,mgmzp), qo(mgmxp,mgmzp),   &
       p(mgmxp,mgmzp), po(mgmxp,mgmzp), us(mgmxp,mgmzp), vs(mgmxp,mgmzp),  &
       omeg(mgmxp,mgmzp),tkeg(mgmxp,mgmzp),rcpg(mgmxp,mgmzp)
  integer pblidx(mgmxp)
  real massfln(mgmxp,mgmyp,ensdim)
  real massflx(mgmxp,mgmyp), xland(mgmxp,mgmyp)
  integer iact_gr(mgmxp,mgmyp), iact_old_gr(mgmxp,mgmyp), kdet(mgmxp)
  integer mynum, i0, j0, m1, m2, m3

  !-------Salva parametros da CUP para uso no transporte convectivo:
  real  edt_average

  integer, dimension(m2) :: ierr4d, jmin4d, kdet4d, k224d, kbcon4d,  &
                            ktop4d, kpbl4d, kstabi4d, kstabm4d

  real, dimension(m2) :: xmb4d, edt4d

  real, dimension(m1,m2) ::  zcup5d,pcup5d
  real, dimension(m1,m2) ::  prup5d,clwup5d,tup5d
  real, dimension(m1,m2) :: enup5d, endn5d, deup5d
  real, dimension(m1,m2) :: dedn5d, zup5d, zdn5d



  !------Variables saved in RAMS Analisys
  ! use (m1,m2,m3) para dimensionar os vetores que sao
  ! escritos nas analises do RAMS
  real, dimension(m2,m3) :: upmf, dnmf, xierr, xktop, xkbcon, xjmin,    &
       xkdt, xiact_p, xiact_c,xk22

  !
  !--- Work variables - Allocatable in this point --
  !
  integer fquasi, fstab, fmconv, iresult
  integer ki, ip1, jp1, m
  integer i, j, k, istart, iend, kk

  real mbdt
  !srf - added on 03-mar-2002

  !--- Output variables ----------------------------

  real outqc(mgmxp,mgmzp), outt(mgmxp,mgmzp), outq(mgmxp,mgmzp)

  !--- Variables on cloud levels ------------------------------------

  real dz, tcrit, xl, pbcdif, outtes, dtime, outteq
  !srf_tmp
  real dellaqsum, dellaqcsum, dp
  !srf_tmp

  !--- New entrainment/detrainment related stuff --------------------

  real mentr_rate, mentrd_rate, entr_rate, entrd_rate,  &
       massfld, edtmax, edtmin, zktop, dh

  !----------------------End of memory allocation ----

  !
  !     if(ktau.gt.3.and.ktau.lt.7) ...
  !
  !
  !--- specify entrainment rate and detrainment rate
  !

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
        !cd(i,k) = (0.5-float(iens)*.1)*entr_rate
        !cd(i,k) = 0.
        cdd(i,k) = 0.
     enddo
  enddo
  !
  !--- max/min allowed value for epsilon
  !    (ratio downdraft base mass flux/updraft
  !    base mass flux
  !
  edtmax=.95
  edtmin=.2

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
        xmass(i,j)=0.
        predb(i,j)=0.
     endif
  enddo

  !--- first check for upstream convection

  do i=istart,iend
     if (ierr(i) == 0) then
        iresult=0
        massfld=0.
        !srf  CALL cup_direction2(i,j,direction,iact_old_gr,mix,mjx,  &
        !srf  mgmxp,mgmyp,massflx,iresult,ensdim,0,0,maxens3,massfld)

        cap_max(i)=cap_maxs
        if (iresult == 1) then
           cap_max(i)=cap_maxs+20.
        endif
     endif
  enddo

  do nens=1,maxens
     mbdt_ens(nens)=(float(nens)-3.)*dtime*1.e-3+dtime*5.E-03
     !        mbdt_ens(nens)=(float(nens)-1.5)*dtime*2.e-3+dtime*5.E-03
  enddo
  do nens=1,maxens2
     !        edt_ens(nens)=.7-float(nens)*.1
     edt_ens(nens)=.95-float(nens)*.01
  enddo
  !--- environmental conditions, FIRST HEIGHTS
  !

  do i=istart,iend
     if(ierr(i) /= 20)then
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
  call cup_env(j, zo, qeso, heo, heso, tn, qo, po, z1, mix,  &
       mgmxp, mkx, mgmzp, istart, iend, psur, ierr, tcrit, 0)
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

        kbmaxloop: do k=1,mkx
           if (zo_cup(i,k).gt.zkbmax+z1(i)) then
              kbmax(i)=k
              exit kbmaxloop
           endif
        end do kbmaxloop

        !--- level where detrainment for downdraft starts

        kdetloop: do k=1,mkx
           if (zo_cup(i,k).gt.z_detr+z1(i)) then
              kdet(i)=k
              exit kdetloop
           endif
        end do kdetloop

     endif
  enddo

  !--- DETERMINE LEVEL WITH HIGHEST MOIST STATIC ENERGY CONTENT - K22
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
     call maximi(heo_cup,mix,mgmxp,mkx,mgmzp,3,kbmax,k22,istart,iend,ierr)
  endif
  do I=ISTART,IEND
     if (ierr(I).eq.0.) then
        if (K22(I).ge.KBMAX(i)) ierr(i)=2
     endif
  enddo

  !--- DETERMINE THE LEVEL OF CONVECTIVE CLOUD BASE  - KBCON

  !Grell sugested to test cup_kbcon -18dec2001

  call cup_kbcon(1, k22, kbcon, heo_cup, heso_cup, mix, mgmxp,  &
       mkx, mgmzp, istart, iend, ierr, kbmax, po_cup, cap_max)


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

  !--- Determine cloud top - KTOP
  !
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
              exit
           endif
        enddo
     endif
  enddo

  !--- Downdraft originating level - JMIN

  call MINIMI(HEso_cup, mix, mgmxp, mkx, mgmzp, K22, kzdown,  &
       JMIN, ISTART, IEND, ierr)

  do I=ISTART,IEND
     if (ierr(I).eq.0.) then

        !--- Check whether it would have buoyancy, if there where
        !--- no entrainment/detrainment

101     continue
        if (jmin(i)-1.lt.KDET(I)) kdet(i)=jmin(i)-1
        if (jmin(i).ge.Ktop(I)-1) jmin(i)=ktop(i)-2
        ki=jmin(i)


        !     hcdo(i,ki)=heo_cup(i,ki)
        if (jmin(i) <= 3) then
           ierr(i) = 9
           go to 100
        end if
        hcdo(i,ki) = heso_cup(i,ki)
        DZ         = Zo_cup(i,Ki+1)-Zo_cup(i,Ki)
        dh         = dz*(HCDo(i,Ki)-heso_cup(i,ki))
        dh         = 0.

        do k=ki-1,1,-1
           !        hcdo(i,k)=heo_cup(i,jmin(i))
           hcdo(i,k) = heso_cup(i,jmin(i))
           DZ        = Zo_cup(i,K+1)-Zo_cup(i,K)
           dh        = dh+dz*(HCDo(i,K)-heso_cup(i,k))
           if (dh.gt.0.) then
              jmin(i)=jmin(i)-1
              if (jmin(i).gt.3) then
                 GO TO 101
              else if (jmin(i).le.3) then
                 ierr(i)=9
                 !CYCLE  !GO TO 100
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
     if (ierr(I).eq.0.) then
        if (-zo_cup(I,KBCON(I))+zo_cup(I,KTOP(I)).lt.  &
             depth_min) then
           ierr(i)=6
        endif
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
        if (aa1(i).eq.0.)  ierr(i)=17
     endif
  enddo

  !--- Determine downdraft strength in terms of windshear

  ! Init. array EDTC with Zeros - ALF

  do k=1,maxens2
     do i=1,mgmxp
        EDTC(i,k) = 0.
     enddo
  enddo

  call cup_dd_edt(ierr, us, vs, zo, ktop, kbcon, edt, po,   &
       pwavo, pwevo, mix, mgmxp, mkx, mgmzp, istart, iend,  &
       edtmax, edtmin, maxens2, edtc, vshear, sdp, vws)

  !srf - Big loop starts here!

  do iedt=1,maxens2
     !DO 250 iedt=1,maxens2
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
           dellaqc_ens(i,k,iedt)= 0.
           pwo_ens(i,k,iedt)    = 0.
        enddo
     enddo

     do I=ISTART,IEND
        aad(i)=0.
     enddo

     call cup_dd_aa0(edto, ierr, aad, jmin, gammao_cup, tn_cup, &
          hcdo, heso_cup, zo, mix, mgmxp, mkx, mgmzp, istart,   &
          iend, zdo)

     !--- Change per unit mass that a model cloud would
     !    modify the environment

     !--- 1. in bottom layer
     call cup_dellabot(heo_cup, ierr, zo_cup, po_cup,   &
          hcdo, edto, zdo, cdd, heo, mix, mgmxp, mkx, mgmzp,  &
          istart, iend, dellah, 1, j, mentrd_rate, zo)
     call cup_dellabot(qo_cup, ierr, zo_cup,    &
          po_cup, qrcdo, edto, zdo, cdd, qo, mix, mgmxp,  &
          mkx, mgmzp, istart, iend, dellaq, 2, j,         &
          mentrd_rate, zo)

     !--- 2. everywhere else

     call cup_dellas(ierr, zo_cup, po_cup, hcdo, edto, zdo,  &
          cdd, heo, mix, mgmxp, mkx, mgmzp, istart, iend,    &
          dellah, 1, j, mentrd_rate, zuo, cd, hco, ktop,     &
          k22, kbcon, mentr_rate, jmin, heo_cup, kdet, k22,  &
          'deep')

     !-- Take out cloud liquid water for detrainment

     do k=1,mkx
        do i=istart,iend
           scr1(i,k)=0.
           dellaqc(i,k)=0.
           if (ierr(i).eq.0) then
              scr1(i,k)=qco(i,k)-qrco(i,k)

              if (k.eq.ktop(i)-0)                    &
                   dellaqc(i,k)=.01*zuo(i,ktop(i))*  &
                   qrco(i,ktop(i))*g/(po_cup(i,k)-po_cup(i,k+1))

              if (k.lt.ktop(i).and.k.gt.kbcon(i)) then
                 dz = zo_cup(i,k+1)-zo_cup(i,k)
                 dellaqc(i,k) = .01*g*cd(i,k)*dz*zuo(i,k)*  &
                      .5*(qrco(i,k)+qrco(i,k+1))/              &
                      (po_cup(i,k  )-po_cup(i,k+1))
              endif
           endif
        enddo
     enddo

     call cup_dellas(ierr, zo_cup, po_cup, qrcdo, edto, zdo,  &
          cdd, qo, mix, mgmxp, mkx, mgmzp, istart, iend,      &
          dellaq, 2, j, mentrd_rate, zuo, cd, scr1, ktop,     &
          k22, kbcon, mentr_rate, jmin, qo_cup, kdet, k22,    &
          'deep')


     !--- Using dellas, calculate changed environmental profiles

     do nens=1,maxens
        !DO 200 nens=1,maxens
        mbdt=mbdt_ens(nens)
        do i=istart,iend
           xaa0_ens(i,nens)=0.
        enddo
        do k=1,mkx-1
           do i=istart,iend
              dellat(i,k)=0.
              if (ierr(i).eq.0) then
                 XHE(I,K)   = DELLAH(I,K)*MBDT + HEO(I,K)
                 XQ(I,K)    = DELLAQ(I,K)*MBDT +  QO(I,K)
                 DELLAT(I,K)= (1./cp)*(DELLAH(I,K)-alvl*DELLAQ(I,K))
                 XT_Grell(I,K)    = DELLAT(I,K)*MBDT +  TN(I,K)
                 if (XQ(I,K).le.0.) XQ(I,K)=1.E-08
              endif
           enddo
        enddo
        !
        do i=istart,iend
           if (ierr(i).eq.0) then
              XHE(I,mkx) = HEO(I,mkx)
              XQ(I,mkx)  = QO(I,mkx)
              XT_Grell(I,mkx)  = TN(I,mkx)
              if (XQ(I,mkx).le.0.) XQ(I,mkx) = 1.E-08
           endif
        enddo

        !--- Calculate moist static energy, heights, qes

        call cup_env(j,xz, xqes, xhe, xhes, xt_grell,  &
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
        !
        !
        do i=istart,iend
           if (ierr(i).eq.0) then
              xaa0_ens(i,nens) = xaa0(i)
              nall = (iens-1)*maxens3*maxens*maxens2 +  &
                   (iedt-1)*maxens*maxens3 + (nens-1)*maxens3

              do k=1,mkx
                 if (k.le.ktop(i)) then
                    do nens3=1,maxens3
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
                          pr_ens(i,j,nall+nens3) =  pr_ens(i,j,nall+nens3) +  &
                               pwo(i,k)+edto(i)*pwdo(i,k)
                       endif
                    enddo
                 endif
              enddo

              do nens3=1,maxens3
                 outt_ens(i,j,nall+nens3) = dellat(i,1)
              enddo
           endif
        enddo
     enddo

     !--- LARGE SCALE FORCING

     !------- CHECK wether aa0 should have been zero
     if (iupmethod == 2) then
        do i=istart,iend
          k22x(i) = k22(i)
        end do
     elseif (iupmethod == 1) then
        call MAXIMI(HE_CUP, mix, mgmxp, mkx, mgmzp, 3,   &
             KBMAX, K22x, ISTART, IEND, ierr)
     end if

     do I=ISTART,IEND
        if (ierr(i).eq.0) then
           if (K22x(I).ge.KBMAX(i)) ierr(i)=998
        endif
     enddo

     !--- DETERMINE THE LEVEL OF CONVECTIVE CLOUD BASE - KBCON

!srf - 31-jan-2003 - ensemble of closures

     call cup_forcing_ens_16(aa0,aa1,xaa0_ens,mbdt_ens,dtime,          &
          xmb,ierr,mix,mgmxp,mjx,mgmyp,mkx,mgmzp,istart,iend,xf_ens,j, &
          fquasi,fstab,'deeps',xland,maxens,iens,iedt,maxens2, &
          maxens3,mconv,omeg,zdo,kbcon,zuo,pr_ens,edto,aad,kbcon,       &
          massflx,iact_old_gr,direction,ensdim,                        &
          massfln,massfld,iresult,xff_ens3, xk,p_cup,ktop,icoic)

     do k=1,mkx
        do i=istart,iend
           if (ierr(i).eq.0) then
              dellat_ens(i,k,iedt)  =  dellat(i,k)
              dellaq_ens(i,k,iedt)  =  dellaq(i,k)
              dellaqc_ens(i,k,iedt) = dellaqc(i,k)
              pwo_ens(i,k,iedt)  = pwo(i,k)+edt(i)*pwdo(i,k)
           else
              dellat_ens(i,k,iedt)  = 0.
              dellaq_ens(i,k,iedt)  = 0.
              dellaqc_ens(i,k,iedt) = 0.
              pwo_ens(i,k,iedt)     = 0.
           endif
        enddo
     enddo
  enddo

  !--- FEEDBACK


  call cup_output_ens(xf_ens,ierr,dellat_ens,dellaq_ens,dellaqc_ens,   &
       outt,outq,outqc,pre,pwo_ens,xmb,ktop,mix,mgmxp,mjx,mgmyp,mkx,mgmzp,  &
       istart,iend,j,'deep',maxens2,maxens,iens,pr_ens,outt_ens,    &
       maxens3,ensdim,massfln,xfac1)

  do I=ISTART,IEND
     PRE(I) = max(PRE(I),0.)
  enddo
  !
  !---------------------done------------------------
  !
  !---Salva parametros da CUP para uso no transporte
  !   convectivo:
  do i=istart,iend
     ierr4d(i)   = ierr(i)
     jmin4d(i)   = jmin(i)
     kdet4d(i)   = kdet(i)
     k224d(i)    = k22(i)
     kbcon4d(i)  = kbcon(i)
     ktop4d(i)   = ktop(i)
     kstabi4d(i) = kstabi(i)
     kstabm4d(i) = kstabm(i)
     kpbl4d(i)   = k22(i)    ! por enquanto kpbl==k22
     xmb4d(i)    = xmb(i)
     !srf - Media no ensemble 2 do parametro edt
     edt_average = 0.
     do iedt=1,maxens2
        edt_average = edt_average + edtc(i,iedt)
     enddo

     edt4d(i) = edt_average/float(maxens2)
     do k=1,mkx
        if (iens.eq.1) then
           zcup5d(k,i) = zo_cup(i,k)
           pcup5d(k,i) = po_cup(i,k)
        endif
        enup5d(k,i) = mentr_rate
        endn5d(k,i) = mentrd_rate
        deup5d(k,i) = cd(i,k)
        dedn5d(k,i) = cdd(i,k)
        zup5d(k,i)  = zuo(i,k)
        zdn5d(k,i)  = zdo(i,k)
	prup5d(k,i) = xmb(i)*pwo(i,k) !only for upfradt
	clwup5d(k,i) = qrco(i,k)      !only for upfradt
	tup5d(k,i) = t_cup(i,k)       !>>> em verdade deveria ser a temperatura da parcela
        !>>> de ar no updraft e _NAO_ a temperatura ambiente

     enddo
  enddo

  !-------Salva parametros nas analises do RAMS

  do I=ISTART,IEND
     xierr(i,j)=float(ierr(i))
     if (ierr(i) == 0) then
        upmf(i,j)   = xmb(i)
        ! downdraft mass flux averaged
	dnmf(i,j)  = 0.
        do k=1,ensdim
           dnmf(i,j) = dnmf(i,j) + massfln(i,j,k)
        enddo
        dnmf(i,j) = dnmf(i,j)/float(ensdim)
        !recalcula o parametro edt4d para o transporte convectivo
        !modifique posteriormente para uso direto do fluxo de massa do downdraft
        edt4d(i)  = dnmf(i,j) / ( upmf(i,j) + 1.e-16 )
        xktop(i,j)  = float(ktop(i))
        xkbcon(i,j) = float(kbcon(i))
        xkdt(i,j)   = float(kdet(i))
        xjmin(i,j)  = float(jmin(i))
     elseif (ierr(i) /= 20) then
        upmf(i,j)   = 0.
        dnmf(i,j)   = 0.
        xktop(i,j)  = 0.
        xkbcon(i,j) = 0.
        xkdt(i,j)   = 0.
        xjmin(i,j)  = 0.
     endif
  enddo

  !srf----------------------------------------

  return
end subroutine CUP_enss


!---------------------------------------------
subroutine cup_dellas(ierr, z_cup, p_cup, hcd, edt, zd, cdd, he, mix,   &
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
  real detdo1, detdo2, entdo, dp, dz, mentrd_rate,  &
       mentr_rate, subin, detdo, entup, detup,         &
       subdown, entdoj, entupk, detupk, totmas
  real xsum, xsumt
  i = istart
  do K=2,MKX
     do I=ISTART,IEND
        della(i,k) = 0.
     enddo
  enddo
  xsum=0.
  !      hesum=0.
  xsumt=0.

  do K=2,MKX-1
  !DO 100 K=2,MKX-1
     do I=ISTART,IEND
        if (ierr(i) /= 0) cycle
        if (K > KTOP(I)) cycle
        !
        !--- Specify detrainment of downdraft,
        !    has to be consistent
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

        if (k == jmin(i))  entdoj  = zd(i,k)*edt(i)

        if (k == k22(i)-1) entupk  = zu(i,kpbl(i))

        if (k > kdet(i)) detdo   = 0.

        if (k.eq.ktop(i)-0) then
           detupk  = zu(i,ktop(i))
           subin   = 0.
        endif
        if (k < kbcon(i)) detup   = 0.

        !
        !--- Changed due to subsidence and entrainment
        !
        totmas =subin-subdown+detup-entup-entdo + detdo-entupk-entdoj+detupk
        if (abs(totmas) > 1.e-6) then
           write (unit=*,fmt='(a)')                 '----------- Subroutine Cup_dellas -----------'
           write(unit=*, fmt='(3(a,1x,i3,1x))')     '  K= ',k,'           I=',i,'           J=',j
           write(unit=*, fmt='(2(a,1x,es10.3,1x))') '  subin=  ',    subin,'subdown=',subdown
           write(unit=*, fmt='(2(a,1x,es10.3,1x))') '  detup=  ',    detup,'entup=  ',entup
           write(unit=*, fmt='(2(a,1x,es10.3,1x))') '  entdo=  ',    entdo,'detdo=  ',detdo
           write(unit=*, fmt='(2(a,1x,es10.3,1x))') '  entupk= ',   entupk,'detupk= ',detupk
           write(unit=*, fmt='(2(a,1x,es10.3,1x))') '  zu(k)=  ',  zu(i,k),'zd(k)=  ',zd(i,k)
           write(unit=*, fmt='(2(a,1x,es10.3,1x))') '  zu(kr)= ',zu(i,k+1),'zd(kr)= ',zd(i,k+1)
           write(unit=*, fmt='(2(a,1x,es10.3,1x))') '  entdoj= ',   entdoj,'edt=    ',edt(i)
           write(unit=*, fmt='(1(a,1x,es10.3,1x))') '  totmas= ',   totmas
           write(unit=*, fmt='(a)')                 '---------------------------------------------'
           write(unit=*, fmt='(a)')                 'The model will stop since it is not conserving mass...' 
           stop
        endif


        !srf         dp =  100.*( p_cup(i,k-1)-p_cup(i,k) )
        dp =  100.*( p_cup(i,k)-p_cup(i,k+1) )
        della(i,k)=(subin  *he_cup(i,k+1) - subdown*he_cup(i,k  ) +  &
             detup*.5*( HC(i,K+1)+ HC(i,K)) +                        &
             detdo*.5*(HCD(i,K+1)+HCD(i,K)) -                        &
             entup*he(i,k) - entdo*he(i,k) -                         &
             entupk*he_cup(i,k22(i)) -  entdoj*he_cup(i,jmin(i)) +   &
             detupk*hc(i,ktop(i)))*g/dp
     enddo
  enddo

  return
end subroutine cup_dellas


!-----------------------------------------
subroutine cup_dellabot(he_cup, ierr, z_cup, p_cup, hcd, edt,   &
     zd, cdd, he, mix, mgmxp, mkx, mgmzp, istart, iend, della, itest, j,  &
     mentrd_rate, z)
  use rconstants, only: g
  implicit none
  integer mix, mgmxp, mkx, mgmzp, i, istart, iend, itest, j
  real z_cup(mgmxp,mgmzp), p_cup(mgmxp,mgmzp), hcd(mgmxp,mgmzp),   &
       zd(mgmxp,mgmzp), cdd(mgmxp,mgmzp), he(mgmxp,mgmzp),         &
       della(mgmxp,mgmzp), he_cup(mgmxp,mgmzp), z(mgmxp,mgmzp), edt(mgmxp)
  integer ierr(mgmxp), m
  real detdo1, detdo2, entdo, dp, dz, mentrd_rate, subin, detdo

  do i=istart,iend
     della(i,1)=0.
     if (ierr(i) /= 0) cycle
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
     !
     !srf-----print-------

!100  CONTINUE
  enddo
  return
end subroutine cup_dellabot


!------------------------------------------------------------
subroutine cup_forcing_ens_16(aa0,aa1,xaa0,mbdt,dtime,xmb,ierr,   &
     mix,mgmxp,mjx,mgmyp,mkx,mgmzp,istart,iend,xf,j,fquasi,       &
     fstab,name,xland,maxens,iens,iedt,maxens2,maxens3,   &
     mconv,omeg,zd,k22,zu,pr_ens,edt,aad,kbcon,massflx,		  &
     iact_old_gr,dir,ensdim,massfln,massfld,iresult,xff_ens3,xk,  &
     p_cup,ktop,icoic)
  use rconstants, only: g
  implicit none
  character (LEN=*) name  !CHARACTER *(*) name

  integer k,i,istart,iend,mix,mgmxp,mjx,mgmyp,mkx,mgmzp,j,        &
       maxens,maxens3
  integer ensdim,iens,nall,iedt,maxens2

  !------ ensemble 3 dimension = 16
  integer :: kclim
  integer, parameter :: mkxcrt=15
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
  !------

  integer k22(mgmxp),kbcon(mgmxp),ierr(mgmxp)
  integer iact_old_gr(mgmxp,mgmyp)

  real aa0(mgmxp),aa1(mgmxp),xaa0(mgmxp,maxens),xmb(mgmxp)
  real mbdt(maxens),dtime,dxxf,edt(mgmxp),aad(mgmxp),dir(mgmxp)
  real     xf(mgmxp,mgmyp,ensdim),xland(mgmxp,mgmyp)
  real pr_ens(mgmxp,mgmyp,ensdim)
  real mconv(mgmxp),omeg(mgmxp,mgmzp),zd(mgmxp,mgmzp),            &
       zu(mgmxp,mgmzp)
  real xff_ens3(maxens3),xk(maxens),xff,xff1,xff2,xff3,xff0
  real massflx(mgmxp,mgmyp)
  real massfln(mgmxp,mgmyp,ensdim)
  real xomg,massfld
  integer fquasi,fstab,nens,ne,n,nens3,iresult,iresultd,          &
       iresulte,icoic
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
              go to 9
           endif
        enddo
        if(p_cup(i,ktop(i)).gt.pcrit(1))kclim=1
9       continue
        kclim = min(mkxcrt,max(kclim,1))
        k= max(kclim-1,1)
        aclim1= acrit(kclim)*1.e3
        aclim2= acrit(k)*1.e3
        aclim3= acritt(kclim)*1.e3
        aclim4= acritt(k)*1.e3
        !
        !--- Treatment different for this closure
        !
        if(name.eq.'deeps')then
           !
           xff0       =  (AA1(I)-AA0(I))/dtime
           xff_ens3(1)=  (AA1(I)-AA0(I))/dtime
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
              xomg=     -omeg(i,k)/g   
              if(xomg.gt.xff_ens3(6)) xff_ens3(6)= 0. !xomg
           enddo
           !
           !--- More like Krishnamurti et al.
           !
           xff_ens3(7)=      mconv(i)
           xff_ens3(8)=   .9*mconv(i)
           xff_ens3(9)=  1.1*mconv(i)
           !
           !--- More like Fritsch Chappel or Kain Fritsch (plus triggers)
           !srf - changed at dec/2002 - greater timescale instab. removal
           !xff_ens3(10)=AA1(I)/(60.*20.)
           !xff_ens3(11)=AA1(I)/(60.*30.)
           !xff_ens3(12)=AA1(I)/(60.*40.)
           xff_ens3(10)=AA1(I)/(60.*50.)
           xff_ens3(11)=AA1(I)/(60.*60.)
           xff_ens3(12)=AA1(I)/(60.*70.)

           !
           !--- More original Arakawa-Schubert (climatologic value of aa0)
           !
           xff_ens3(13)=max(0.,(AA1(I)-aclim1)/dtime)
           xff_ens3(14)=max(0.,(AA1(I)-aclim2)/dtime)
           xff_ens3(15)=max(0.,(AA1(I)-aclim3)/dtime)
           xff_ens3(16)=max(0.,(AA1(I)-aclim4)/dtime)

           do nens=1,maxens
              XK(nens)=(XAA0(I,nens)-AA1(I))/MBDT(nens)
              if(xk(nens).le.0.and.xk(nens).gt.-1.e-9) xk(nens)=-1.e-9
              if(xk(nens).gt.0.and.xk(nens).lt.+1.e-9) xk(nens)=+1.e-9
           enddo
           !
           !--- Add up all ensembles
           !
           do ne=1,maxens
              !
              !--- for every xk, we have maxens3 xffs
              !--- iens is from outermost ensemble (most expensive!
              !
              !--- iedt (maxens2 belongs to it)
              !--- is from second, next outermost, not so expensive
              !
              !--- so, for every outermost loop, we have maxens*maxens2*3
              !--- ensembles!!! nall would be 0, if everything is on first
              !--- loop index, then ne would start counting, then iedt,
              !--- then iens....
              !
              iresultd=0
              iresulte=0
              nall=(iens-1)*maxens3*maxens*maxens2  &
                   +(iedt-1)*maxens*maxens3         &
                   +(ne-1)*maxens3
              !
              !--- check for upwind convection
              !
              !iresult=0
              !massfld=0.
              !call cup_direction2(i,j,dir,iact_old_gr,mix,mjx,  &
              !     mgmxp,mgmyp,massflx,iresult,ensdim,1,nall,   &
              !     maxens3,massfld)
              if(XK(ne).lt.0.and.xff0.gt.0.)iresultd=1
              iresulte=max(iresult,iresultd)
              iresulte=1
              if(iresulte.eq.1)then
                 !
                 !--- Special treatment for stability closures
                 !

                 if(xff0.gt.0.)then
                    xf(i,j,nall+1) =max(0., -xff_ens3(1)/xk(ne))+massfld
                    xf(i,j,nall+2) =max(0., -xff_ens3(2)/xk(ne))+massfld
                    xf(i,j,nall+3) =max(0., -xff_ens3(3)/xk(ne))+massfld
                    xf(i,j,nall+13)=max(0.,-xff_ens3(13)/xk(ne))+massfld
                    xf(i,j,nall+14)=max(0.,-xff_ens3(14)/xk(ne))+massfld
                    xf(i,j,nall+15)=max(0.,-xff_ens3(15)/xk(ne))+massfld
                    xf(i,j,nall+16)=max(0.,-xff_ens3(16)/xk(ne))+massfld
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
                 xf(i,j,nall+7)=max(0.,xff_ens3(7)/pr_ens(i,j,nall+7))
                 xf(i,j,nall+8)=max(0.,xff_ens3(8)/pr_ens(i,j,nall+8))
                 xf(i,j,nall+9)=max(0.,xff_ens3(9)/pr_ens(i,j,nall+9))
                 if(XK(ne).lt.0.)then
                    xf(i,j,nall+10)=max(0.,-xff_ens3(10)/xk(ne))+massfld
                    xf(i,j,nall+11)=max(0.,-xff_ens3(11)/xk(ne))+massfld
                    xf(i,j,nall+12)=max(0.,-xff_ens3(12)/xk(ne))+massfld
                 else
                    xf(i,j,nall+10)=massfld
                    xf(i,j,nall+11)=massfld
                    xf(i,j,nall+12)=massfld
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
                    xf(i,j,nall+11)=xf(i,j,nall+icoic)
                    xf(i,j,nall+12)=xf(i,j,nall+icoic)
                    xf(i,j,nall+13)=xf(i,j,nall+icoic)
                    xf(i,j,nall+14)=xf(i,j,nall+icoic)
                    xf(i,j,nall+15)=xf(i,j,nall+icoic)
                    xf(i,j,nall+16)=xf(i,j,nall+icoic)
                 endif
                 !==============
                 !05-12-2002
                 !srf - forcing 14 is too bad, use the same for 13:
                 !A&S (14) = A&S (13)
                 xf(i,j,nall+14)=xf(i,j,nall+13)
                 !==============
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
end subroutine cup_forcing_ens_16
!--------------------------------------------------------------------

!--------------------------------------------------------------------
subroutine cup_output_ens(xf,ierr,dellat,dellaq,dellaqc,    &
     outtem,outq,outqc,pre,pw,xmb,ktop,mix,mgmxp,mjx,mgmyp,   &
     mkx,mgmzp,istart,iend,j,name,nx,nx2,iens,pr_ens, &
     outt_ens,maxens3,ensdim,massfln,xfac1)
  use rconstants, only: day_sec
  implicit none
  character (LEN=*) name  !CHARACTER *(*) name

  integer mix,mjx,mkx,istart,iend,mgmxp,mgmyp,mgmzp, &
       ensdim,i,k,j,nx,n,nx2,m

  real       xf(mgmxp,mgmyp,ensdim), pr_ens(mgmxp,mgmyp,ensdim), &
       outt_ens(mgmxp,mgmyp,ensdim),massfln(mgmxp,mgmyp,ensdim)

  real outtem(mgmxp,mgmzp),outq(mgmxp,mgmzp),outqc(mgmxp,mgmzp), &
       dellat(mgmxp,mgmzp,nx), dellaq(mgmxp,mgmzp,nx),           &
       pw(mgmxp,mgmzp,nx),dellaqc(mgmxp,mgmzp,nx),               &
       pre(mgmxp),xmb(mgmxp),xfac1(mgmxp)

  integer ktop(mgmxp),ierr(mgmxp),ncount,iens,maxens3,nens3
  real outtes,ddtes,dtt,dtq,dtqc,dtpw

  do K=1,MKX
     do I=ISTART,IEND
        outtem(i,k) = 0.
        outq(i,k) = 0.
        outqc(i,k) = 0.
     enddo
  enddo

  do I=ISTART,IEND
     pre(i)  =0.
     xmb(i)  =0.
     xfac1(i)=1.
  enddo
  !
  !--- Calculate ensemble average mass fluxes
  !
  do I=ISTART,IEND
     ncount=0
     xmb(i)=0.

     if(ierr(i).eq.0)then

        do n=(iens-1)*nx*nx2*maxens3+1,iens*nx*nx2*maxens3
           pr_ens(i,j,n) =   pr_ens(i,j,n)*xf(i,j,n)
           outt_ens(i,j,n) = outt_ens(i,j,n)*xf(i,j,n)

           if(xf(i,j,n).gt.0.)then
              xmb(i) = xmb(i) + xf(i,j,n)
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

  enddo
  !
  !-- Now do feedback
  !
  ddtes=250.
  !     if(name.eq.'shal')ddtes=500.
  do I=ISTART,IEND
     xfac1(i)=xmb(i)
  enddo
  do K=1,MKX
     do I=ISTART,IEND
        dtt  =0.
        dtq  =0.
        dtqc =0.
        dtpw =0.

        if(ierr(i).eq.0.and.k.le.ktop(i))then
           do n=1,nx
              dtt  = dtt  +  dellat(i,k,n)
              dtq  = dtq  +  dellaq(i,k,n)
              dtqc = dtqc + dellaqc(i,k,n)
              dtpw = dtpw +      pw(i,k,n)

              !srf-----print-------
              !	   if(k.eq.1) then
              !            write(6,'(a1,78a1)') ' ',('*',m=1,78)
              !	    print*,'nx j i k PREC dtpw  pw '
              !	   endif
              !          write(6,'(4i4,3e12.4)') nx,j,i,k,pre(i),dtpw,pw(i,k,n)
              !     	   if(k.eq.ktop(i)) then
              !            write(6,'(a1,78a1)') ' ',('*',m=1,78)
              !           endif
              !srf-----print-------
           enddo
           outtes = dtt*XMB(I)*day_sec/float(nx)

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

           OUTTEM(I,K) = OUTTEM(I,K) +XMB(I)*dtt /float(nx)
           OUTQ(I,K) =   OUTQ(I,K) +XMB(I)*dtq /float(nx)
           OUTQC(I,K) =  OUTQC(I,K) +XMB(I)*dtqc/float(nx)
           PRE(I)      = PRE(I)      +XMB(I)*dtpw/float(nx)
           ! unit : kg[liq water]/(m^2 s)

           !
           !
           !	   if(k.eq.ktop(i)) then
           !            write(6,'(a1,78a1)') ' ',('-',m=1,78)
           !           endif
           !           if(k.eq.mkx) write(6,'(a1,78a1)') ' ',('-',m=1,78)
           !          endif
           !srf-----print-------

        endif
     enddo
  enddo

  do I=ISTART,IEND
     if(ierr(i).eq.0)then
        xfac1(i)=xmb(i)/xfac1(i)
        do k=1,ensdim
           massfln(i,j,k)=massfln(i,j,k)*xfac1(i)
        enddo
     endif
  enddo

  return
end subroutine cup_output_ens
!------------------------------------------------------------
