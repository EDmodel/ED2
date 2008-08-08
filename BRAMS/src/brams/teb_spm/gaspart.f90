!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

subroutine le_fontes(ng,n1,n2,n3,np,ia,iz,ja,jz)

  use mem_grid
  use mem_leaf
  use mem_basic
  use mem_gaspart

  use rconstants
  implicit none


  integer :: ng,n1,n2,n3,np,ig,kl,i,j,ia,iz,ja,jz

  integer, parameter :: ngases=6
  integer, dimension(ngases) :: kgas
  integer :: len1
  character(len=10) :: GAS(ngases),tracer
  data ( gas(ig),ig=1,ngases) /'NO', 'NO2', 'PM25', 'CO', 'SO2', 'VOCS'/
  data (kgas(ig),ig=1,ngases) /  1,    4,     7,     10,   13,     16/



  !------------ Allocation table of qsc = dum1------------------------
  ! Distribute sources in order to avoid highly localized sources:
  ! (somente definidas em 1 ponto de grade)
  ! kgas(1) = 1,2,3    => NO
  ! kgas(2) = 4,5,6    => NO2
  ! kgas(3) = 7,8,9    => PM25
  ! kgas(4) = 10,11,12 => CO
  ! kgas(5) = 13,14,15 => SO2 
  ! kgas(6) = 16,17,18 => VOCS
  !----------------------------------------------------------------------


  !sources only in the inner grid

  if(ng.ne.ngrids)return

  !Allocate memory in DUM1 vector for sources:

  kl= 3*ngases    ! first kl 3D vector gaspart_g(ng)%gasr(1,1,1) levels are left for 
  ! emission sources calculation


  gaspart_g(ng)%gasr(:,:,:)=0.

  do ig=1,ngases

     !Call dum1_zero(n1,n2,n3,ia,iz,ja,jz,kl,scalar_g(ig,ng)%dum1(1,1,1))  

     tracer=gas(ig)

     len1 = len_trim(tracer)+ 1


     call read_sources_teb(ng,n1,n2,n3,np,ia,iz,ja,jz,gaspart_g(ng)%gasr(1,1,1),gaspart_g(ng)%fusog(1,1),&
          tracer(1:len1),kgas(ig),leaf_g(ng)%G_URBAN(1,1,1),grid_g(ng)%dxt(1,1), &
          grid_g(ng)%dyt(1,1),time)
     call reorganize_sources_teb(n1,n2,n3,ia,iz,ja,jz,gaspart_g(ng)%gasr(1,1,1),kgas(ig))


     call convert_to_misture_ratio_teb(ng,n1,n2,n3,zt,ia,iz,ja,jz,kgas(ig),gaspart_g(ng)%gasr(1,1,1),    &
          basic_g(ng)%dn0(1,1,1),grid_g(ng)%rtgt(1,1), &
          grid_g(ng)%dxt(1,1),grid_g(ng)%dyt(1,1),dzt)


  enddo !end of gases' looping


  return
end subroutine le_fontes
!--------------------------------------------------
subroutine read_sources_teb(ng,n1,n2,n3,np,ia,iz,ja,jz,qsc,fuso,gas,kgas,schar,dxt,dyt,time)

  use mem_emiss, only:EINDNO,EINDNO2,EINDPM,EINDCO,EINDSO2,EINDVOC,  &
       EVEINO,EVEINO2,EVEIPM,EVEICO,EVEISO2,EVEIVOC,  &
       EFSAT,EFSUN,WEEKDAYIN     

  use teb_vars_const, only : RUSHH1,RUSHH2,DAYLIGHT

  implicit none

  integer :: ia,iz,ja,jz,i,j,n1,n2,n3,np,ng,kgas,idays
  real, dimension(n1,n2,n3) :: qsc
  real, dimension(n2,n3,np) ::schar
  real, dimension(n2,n3) :: dxt,dyt,fuso
  character*(*) gas
  character(len=3)cday
  real ::pfat,pfat2,emiss,emiind,area,pft,r_q,tign,ax1,ax2,bx1,bx2
  real ::timeq1,timeq2
  real(kind=8) :: time

  !*****************************************************************************
  !defining the emission rate for each gas/particle
  !*****************************************************************************

  !The following values are provided by CETESB's report for 2001
  !these emission rates are relative to an area of 8 million m2, wich 
  !represent the total area of RMSP. However, the urbanized part of this
  !area is about 1.5 and it will be considered here.

  !    veicular emissions rate in ug/min/m2 (for an area of 8 million of m2)

  !             CO   = 390.6
  !             NOx  =  89.3
  !             NO   = NOx*0.9
  !             NO2  = NOx*0.1
  !             VOC's=  91.0
  !             SO2  =   5.3

  ! In order to consider the diurnal cycle of veicular emission, the emission rates will
  ! have units of kg/m2/day

  !    veicular emissions rate in kg/day/m2 (for an area of 1.5 thousand of m2)

  !             CO   = 3.0173515E-03
  !             NOx  = 6.8566209E-04
  !             NO   = NOx*0.9 =  6.1709585E-04
  !             NO2  = NOx*0.1 =  6.8566209E-05
  !             VOC's= 6.9954334E-04
  !             SO2  = 4.0730592E-05
  !             PM15 = 6.2648396E-06
  !    industrial emissions rate in kg/s/m2 (for an area of 1.5 thousand of m2)
  !
  !  emiind
  !             CO   = 8.1599860E-10
  !             NOx  = 6.8566209E-04
  !             NO   = 2.6636227E-10
  !             NO2  = 2.9595805E-11
  !             VOC's= 2.5367833E-10
  !             SO2  = 3.6149164E-10
  !             PM25 = 4.3421278E-10

  if(gas=='CO')then

     emiss=EVEICO 
     emiind=EINDCO

  endif
  if(gas=='NO')then

     emiss=EVEINO
     emiind=EINDNO


  endif
  if(gas=='NO2')then
     emiss=EVEINO2
     emiind=EINDNO2

  endif
  if(gas=='VOCS')then

     emiss=EVEIVOC
     emiind=EINDVOC

  endif
  if(gas=='SO2')then

     emiss=EVEISO2
     emiind=EINDSO2

  endif
  if(gas=='PM25')then

     emiss=EVEIPM
     emiind=EINDPM

  endif
  emiss=emiss

  pft=emiss/273234.9 !

  !the value of 273234.9 was obtained in order to have the integral for
  !one day = emiss
  !it is used to distribute emissions following the diurnal cycle 

  tign = 0. * 3600.                            !UTC time of ignition 

  !print*,'valor de time=',time

  idays = int((time/3600.)/24.+.00001)  !number of days of simulation
  tign = tign + real(idays)*24.*3600.

  pfat=1.
  call EMFACTOR(WEEKDAYIN,idays,cday)
  if(cday=='SAT')pfat=EFSAT
  if(cday=='SUN')pfat=EFSUN

  !******************************************************************
  !gaussian distribution                                            *
  !                                     timeqi                      *
  !                                      ____                       *
  !                               _     |    |      _               *
  !              1               |      (x-mi)**2    |              *
  !f(x)= ----------------- * EXP | -  -------------- |  i=1,2       *
  !      sigmai* sqrt(2*PI)      |    2*(sigmai**2)  |              *
  !     |__________________|      -  |____________| -               *
  !             axi                    8.5 or 10.6                  *
  !******************************************************************

  ax1= 4.5
  ax2= 9.2

  do i = 1,n2
     do j= 1,n3


        if(nint(SCHAR(I,J,2))/=0)then

           !bx1=10.81
           !bx2=19.0
           bx1=RUSHH1-fuso(i,j)+DAYLIGHT
           bx2=RUSHH2-fuso(i,j)+DAYLIGHT

           timeq1= ( time /3600. - tign/3600.) - bx1
           timeq2= ( time /3600. - tign/3600.) - bx2

           r_q=pft*( ax1*exp(-(timeq1)**2/8.5)   +      &
                ax2*exp(-(timeq2)**2/10.6) )*pfat 

           area=1./(dxt(i,j)*dyt(i,j))
           ! Identificando pontos de grade com a classificacao urbano 1 
           if(nint(SCHAR(I,J,2))==1)then

              pfat2=1.

           endif
           ! Identificando pontos de grade com a classificacao urbano 2 
           if(nint(SCHAR(I,J,2))==2)then

              pfat2=0.33333333

           endif
           ! Identificando pontos de grade com a classificacao urbano 3 
           if(nint(SCHAR(I,J,2))==3)then

              pfat2=0.2

           endif

           qsc(kgas,i,j)= (r_q*area*pfat*pfat2)+(emiind*area)

        endif

     enddo
  enddo


  return
end subroutine read_sources_teb
!-------------------------------------------------------
!
subroutine reorganize_sources_teb(n1,n2,n3,ia,iz,ja,jz,qsc,kgas)

  implicit none
  real, dimension(n1,n2,n3) ::qsc
  integer :: n1,n2,n3,i,j,ii,jj,kgas,ia,iz,ja,jz
  real :: f


  !fator de distribuicao 20% para cada um do 9 primeiros vizinhos
  !( incluindo o proprio site i,j)

  f=0.2
  do j=3,n3-2
     do i=3,n2-2                            

        qsc(kgas+1,i,j) = qsc(kgas+1,i,j) + (1.-f)* qsc(kgas,i,j)  ! ponto grade i,j,k(=2)

        !distribuicao nos 9 sites em torno de i,j    !  j+1  .   .   .
        do jj = j-1,j+1                          !   j   .   .   .  K=2     
           do ii = i-1,i+1                        !  j-1  .   .   .
              !      i-1  i  i+1

              qsc(kgas+1,ii,jj) = qsc(kgas+1,ii,jj) & 
                   + (1./9.) * f * qsc(kgas,i,j)  ! ponto grade i,j,k= 2 e 3

           enddo
        enddo

     enddo
  enddo


  return
end subroutine reorganize_sources_teb
!-------------------------------------------------------------
!
subroutine convert_to_misture_ratio_teb(ng,n1,n2,n3,zt,ia,iz,ja,jz,kgas,qsc, &
     dn0,rtgt,dxt,dyt,dzt)

  implicit none
  integer :: ng,n1,n2,n3,kgas,i,j,k,ia,iz,ja,jz   
  real, dimension(n2,n3) :: dxt,dyt,rtgt
  real, dimension(n1,n2,n3):: dn0,qsc
  real, dimension(n1) ::  dzt,zt
  real :: fcu,vol	
  !

  !Fator de conversao de unidades
  fcu=1.        !=> kg [gas/part] /kg [ar]
  !!fcu =1.e+12   !=> ng [gas/part] /kg [ar]

  !fcu =1.e+6      !=> mg [gas/part] /kg [ar]  

  do j= 1,n3
     do i = 1,n2
        do k=2,2

           vol = 1./(dxt(i,j)*dyt(i,j)*dzt(k)*rtgt(i,j))

           qsc(kgas+1,i,j) = fcu* qsc(kgas+1,i,j)/(vol*dn0(k,i,j))

        enddo
     enddo
  enddo

  return
end subroutine convert_to_misture_ratio_teb

!------------------------------------------------------------
subroutine sources_teb(n1,n2,n3,ia,iz,ja,jz,ig,deltat)

  use mem_grid, only : ngrids

  !USE mem_scalar
  use mem_tend
  use mem_gaspart
  implicit none
  integer :: ia,iz,ja,jz,ig,n1,n2,n3,igas
  real    :: deltat

  if (ig==ngrids)then

     call EMISSAO(n1,n2,n3,ia,iz,ja,jz,gaspart_g(ig),deltat,ig)

  endif

  return
end subroutine sources_teb

!------------------------------------------------------------

subroutine EMISSAO(n1,n2,n3,ia,iz,ja,jz,gaspart,deltat,ng)

  use mem_gaspart
  implicit none

  type(gaspart_vars) gaspart

  integer :: ia,iz,ja,jz,ng,n1,n2,n3,i,j,k

  real :: deltat

  call tendgas(n1,n2,n3,ia,iz,ja,jz,gaspart%pnot  (1),gaspart%gasr(1,1,1),2)
  call tendgas(n1,n2,n3,ia,iz,ja,jz,gaspart%pno2t (1),gaspart%gasr(1,1,1),5)
  call tendgas(n1,n2,n3,ia,iz,ja,jz,gaspart%ppm25t(1),gaspart%gasr(1,1,1),8)
  call tendgas(n1,n2,n3,ia,iz,ja,jz,gaspart%pcot  (1),gaspart%gasr(1,1,1),11)
  call tendgas(n1,n2,n3,ia,iz,ja,jz,gaspart%pso2t (1),gaspart%gasr(1,1,1),14)
  call tendgas(n1,n2,n3,ia,iz,ja,jz,gaspart%pvoct (1),gaspart%gasr(1,1,1),17)

  return
end subroutine EMISSAO
!---------------------------------------------------------------
subroutine TENDGAS(n1,n2,n3,ia,iz,ja,jz,tende,qsc,k)

  implicit none

  integer         ::n1,n2,n3,ia,iz,ja,jz,i,j,k,ii,jj

  real,dimension(n1,n2,n3) :: tende,qsc
  real val


  do i=ia,iz
     do j=ja,jz


        tende(2,i,j)=tende(2,i,j)+qsc(k,i,j)

     enddo
  enddo

  return
end subroutine TENDGAS

!---------------------------------------------------------------
subroutine dum1_zero(n1,n2,n3,ia,iz,ja,jz,kl,d)

  implicit none
  integer :: n1,n2,n3,ia,iz,ja,jz,kl,i,j,k
  real, dimension(n1,n2,n3) :: d

  do j=ja,jz
     do i=ia,iz
        do k=1,kl 
           d(k,i,j) = 0.
        enddo
     enddo
  enddo
  return
end subroutine dum1_zero
!----------------------------------------------------------------
subroutine init_conc1(ictrl,ng,n1,n2,n3,np, &
     G_URBAN,no,           &
     no2,pm25,             &
     co,voc,               &
     so2,so4,aer,zt)

  implicit none

  integer                   :: ictrl,ng,n1,n2,n3,np,i,j,ii,jj,k

  real, dimension(n1,n2,n3) :: no,no2,pm25,co,voc,so2,so4,aer,    &
       nox,no2x,pm25x,cox,vocx,so2x
  real                      :: no0,no20,pm250,co0,voc0,so20

  real, dimension (n2,n3,np):: G_URBAN

  real, dimension (n1)      :: zt

  real, parameter           :: PMAR=28.97

  real pfat,f,expo

  no(:,:,:) =0.			       !NO
  no2(:,:,:) =0.			       !NO2
  pm25(:,:,:) =0.			       !PM25
  co(:,:,:) =0.			       !CO
  so2(:,:,:) =0.			       !SO2
  so4(:,:,:) =0.			       !SO4
  voc(:,:,:) =0.			       !VOCS
  aer(:,:,:) =0.			       !VOCS
  if(ictrl.eq.1) then
     !--------
     !Flag which defines initial condition (IF RUNTYPE = INITIAL)
     do j=1,n3
        do i=1,n2
           pfat=0.0 

           !identify surface type by using G_URBAN 
           if(nint(G_URBAN(i,j,2))/=0.)then

              if(nint(G_URBAN(i,j,2))==1)then
                 pfat=1.
              endif

              if(nint(G_URBAN(i,j,2))==2)then
                 pfat=0.3     !using only 30% of values for urban regions
              endif
              if(nint(G_URBAN(i,j,2))==3)then
                 pfat=0.2     !using only 20% of values for urban regions
              endif
           else
              pfat=0.1     
           endif

           !defining surface backgroud concentrations (typical values for 00 Z in urban regions)
           !can be better difined if using real values (eg. read a concentration file)

           no0= 248.04 * pfat  !NO  (estacao congonhas dia 31/07/99 ug/m3)
           no20=  62.30 * pfat  !NO2 (estacao congonhas dia 31/07/99 ug/m3)
           pm250= 17.28 * pfat  !PM25(estacao santana   dia 31/07/99 ug/m3)
           co0=   1.20 * pfat  !CO  (estacao congonhas dia 31/07/99 ppm)
           so20=  17.50 * pfat  !SO2 (estacao congonhas dia 31/07/99 ug/m3)
           voc0=   0.31 * pfat  !VOCS (Leila)


           do k=2,n1-1
              expo=exp(-zt(k)/zt(n1-1))

              nox  (k,i,j)= (no0*expo)  *(1.e-9)/1.275	   !NO
              no2x (k,i,j)= (no20*expo) *(1.e-9)/1.275	   !NO2
              pm25x(k,i,j)= (pm250*expo)*(1.e-9)/1.275	   !PM25
              cox  (k,i,j)= (co0*expo)  *(28.00/PMAR)*1.e-6 !CO
              so2x (k,i,j)= (so20*expo) *(1.e-9)/1.275	   !SO2
              vocx (k,i,j)= (voc0*expo) *(42.08/PMAR)*1.e-6  !VOCS
           enddo  !end of looping in k(vertical levels)

        enddo  !end of looping in i (longitude)
     enddo   !end of looping in j (latitute)

  endif  !end of if in ictrl

  f=0.2

  !re-distributing sources

  do j=3,n3-2
     do i=3,n2-2                            

        do k=2,n1-1


           no  (k,i,j)  = no  (k,i,j)  + (1.-f)* nox  (k,i,j)
           no2 (k,i,j)  = no2 (k,i,j)  + (1.-f)* no2x (k,i,j)
           pm25(k,i,j)  = pm25(k,i,j)  + (1.-f)* pm25x(k,i,j)
           co  (k,i,j)  = co  (k,i,j)  + (1.-f)* cox  (k,i,j)
           so2 (k,i,j)  = so2 (k,i,j)  + (1.-f)* so2x (k,i,j)
           voc (k,i,j)  = voc (k,i,j)  + (1.-f)* vocx (k,i,j)

           !distribution into the 9 sites around i,j    !  j+1  .   .   .
           do jj = j-1,j+1                          !   j   .   .   .       
              do ii = i-1,i+1                        !  j-1  .   .   .
                 !      i-1  i  i+1

                 no  (k,ii,jj)  =  no  (k,ii,jj)  +(1./9.)*f*nox  (k,i,j)
                 no2 (k,ii,jj)  =  no2 (k,ii,jj)  +(1./9.)*f*no2x (k,i,j)
                 pm25(k,ii,jj)  =  pm25(k,ii,jj)  +(1./9.)*f*pm25x(k,i,j)
                 co  (k,ii,jj)  =  co  (k,ii,jj)  +(1./9.)*f*cox  (k,i,j)
                 so2 (k,ii,jj)  =  so2 (k,ii,jj)  +(1./9.)*f*so2x (k,i,j)
                 voc (k,ii,jj)  =  voc (k,ii,jj)  +(1./9.)*f*vocx (k,i,j)
                 so4 (k,ii,jj)  =  so2 (k,ii,jj)*0.2  
                 aer (k,ii,jj)  =  so2 (k,ii,jj)*0.0 

              enddo
           enddo

        enddo !vertical
     enddo  !longitude
  enddo   !latitude


  return
end subroutine init_conc1
!----------------------------------------------------------------

subroutine init_conc2(ictrl,ng,n1,n2,n3,np,G_URBAN,s7p,s8p,s9p,s10p,s11p,s12p,s13p,zt)

  integer :: ictrl,ng,n1,n2,n3,np,i,j,ii,jj,k
  real, dimension(n1,n2,n3) ::  s7p, s8p, s9p, s10p, s11p, s12p,s13p,      &
       s7p2,s8p2,s9p2,s10p2,s11p2,s12p2,s13p2
  real                      :: s7p0,s8p0,s9p0,s10p0,s11p0,s12p0,s13p0

  real, dimension (n2,n3,np):: G_URBAN

  real, dimension (n1)      :: zt

  real, parameter           :: PMAR=28.97

  real pfat,f,expo


  if(ictrl.eq.1) then
     !--------
     !Flag which defines initial condition (IF RUNTYPE = INITIAL)


     s7p(:,:,:) =0.			       !O3
     s8p(:,:,:) =0.			       !RHCO
     s9p(:,:,:) = 0.                           !HO2
     s10p(:,:,:) = 0.                          !O3P
     s11p(:,:,:) = 0.                          !O1D
     s12p(:,:,:) = 0.                          !HO
     s13p(:,:,:) = 0.                          !RO2


     do j=1,n3
        do i=1,n2


           !identify surface type by using G_URBAN 

           if(nint(G_URBAN(i,j,2))/=0.)then

              if(nint(G_URBAN(i,j,2))==1)then
                 pfat=1.
              endif

              if(nint(G_URBAN(i,j,2))==2)then
                 pfat=0.3     !using only 30% of values for urban regions
              endif
              if(nint(G_URBAN(i,j,2))==3)then
                 pfat=0.2     !using only 20% of values for urban regions
              endif
           else
              pfat=0.1     
           endif

           !defining surface backgroud concentrations (typical values for 00 Z in urban regions)
           !can be better difined if using real values (eg. read a concentration file)

           s7p0=   8.00 * pfat   !O3  (estacao santana   dia 31/07/99 ug/m3)
           s8p0=   0.0187* pfat  !RHCO (valor aproximado segundo Leila,2003, comunicacao pessoal)
           s9p0=   (4.*1e-7)*pfat  !HO2  ppm
           s10p0=  (4.8*1e-10)*pfat !o3p ppm
           s11p0=  (4.8*1e-11)*pfat !o1d ppm
           s12p0=  (8.*1e-7)*pfat  !HO ppm
           s13p0=	(4.*1e-7)*pfat  !RO2 ppm



           do k=2,n1-1
              expo=exp(-zt(k)/zt(n1-1))

              s7p2(k,i,j)= (s7p0*expo)*(1.e-9)/1.275	   !O3
              s8p2(k,i,j)= (s8p0*expo)*(44.05/PMAR)*1.e-6   !RHCO
              s9p2(k,i,j)= (s9p0*expo)*(33.0/PMAR)*1.e-6    !HO2
              s10p2(k,i,j)=(s10p0*expo)*(16.0/PMAR)*1.e-6   !O3P
              s11p2(k,i,j)=(s11p0*expo)*( 16.0/PMAR)*1.e-6  !O1D
              s12p2(k,i,j)=(s12p0*expo)*( 17.0/PMAR)*1.e-6  !HO
              s13p2(k,i,j)=(s13p0*expo)*( 47.0/PMAR)*1.e-6  !RO2


           enddo


        enddo
     enddo

  endif

  f=0.2


  do j=3,n3-2
     do i=3,n2-2                            

        do k=2,n1-1

           s7p(k,i,j) = s7p(k,i,j) + (1.-f)* s7p2(k,i,j)
           s8p(k,i,j) = s8p(k,i,j) + (1.-f)* s8p2(k,i,j)
           s9p(k,i,j) = s9p(k,i,j) + (1.-f)* s9p2(k,i,j)
           s10p(k,i,j) =s10p(k,i,j) + (1.-f)* s10p2(k,i,j)
           s11p(k,i,j) =s11p(k,i,j) + (1.-f)* s11p2(k,i,j)
           s12p(k,i,j) =s12p(k,i,j) + (1.-f)* s12p2(k,i,j)
           s13p(k,i,j) =s13p(k,i,j) + (1.-f)* s13p2(k,i,j)

           !distribution into the 9 sites around i,j    !  j+1  .   .   .
           do jj = j-1,j+1                          !   j   .   .   .  K=2     
              do ii = i-1,i+1                        !  j-1  .   .   .
                 !      i-1  i  i+1

                 s7p(k,ii,jj) =  s7p(k,ii,jj)+(1./9.)*f*s7p2(k,i,j)
                 s8p(k,ii,jj) =  s8p(k,ii,jj)+(1./9.)*f*s8p2(k,i,j)
                 s9p(k,ii,jj) =  s9p(k,ii,jj)+(1./9.)*f*s9p2(k,i,j)
                 s10p(k,ii,jj) = s10p(k,ii,jj)+(1./9.)*f*s10p2(k,i,j)
                 s11p(k,ii,jj) = s11p(k,ii,jj)+(1./9.)*f*s11p2(k,i,j)
                 s12p(k,ii,jj) = s12p(k,ii,jj)+(1./9.)*f*s12p2(k,i,j)
                 s13p(k,ii,jj) = s13p(k,ii,jj)+(1./9.)*f*s13p2(k,i,j)

              enddo
           enddo

        enddo

     enddo
  enddo

  return
end subroutine init_conc2
!----------------------------------------------------------------
!----------------------------------------------------------------
subroutine EMFACTOR(dayin,idays,dayout)
  implicit none
  integer ::i,j,idays,k,id,ndays,irest
  character(len=3) :: cday(7),dayin,dayout
  data ( cday(id),id=1,7) /'SUN', 'MON', 'TUE', 'WED', 'THU', 'FRI', 'SAT'/
  !                                               0      1      2      3   
  !			  4      5      6      7      8      9      10
  !			  11     12     13     14     15     16     17
  !			  18     19     20     21     22     23     24


  ndays=(idays/7)*7  !   (21)
  irest=idays-ndays  !  (3)

  do i=1,7
     if(dayin==cday(i))j=i
  enddo
  k=j+irest

  if(k>7)k=k-7

  dayout=cday(k)

  return
end subroutine EMFACTOR

!###########################################################################

subroutine init_conc_prev(name_name)

  ! This routine initializes gas variables from a previous (day -1) history file

  use grid_dims
  use var_tables
  use io_params
  use mem_grid
  use ref_sounding
  use mem_emiss, only:chemdata_in

  implicit none

  character (len=*) :: name_name

  integer :: ngrids1,ioutput1  &
       ,nnxp1(maxgrds),nnyp1(maxgrds),nnzp1(maxgrds),nzg1,nzs1,npatch1

  integer :: iyr,imn,idy,itm,ie,maxarr,ngr,nc
  character (len=80) :: hnameinh,prefix
  character (len=2) :: cng
  integer, external :: cio_i,cio_f,cio_i_sca
  integer,save :: iunhd=11


  ! Open the input history header file and read some of the info.

  write(*,*)'chemdata_in=',chemdata_in
  !pause

  nc=len_trim(chemdata_in)
  hnameinh=chemdata_in(1:nc-9)//'.vfm'

  call rams_f_open(iunhd,chemdata_in,'FORMATTED','OLD','READ',0)

  ie=cio_i_sca(iunhd,1,'ngrids',ngrids1,1)
  ngridsh=ngrids1
  ie=cio_i(iunhd,1,'nnxp',nnxp1,ngrids1)
  ie=cio_i(iunhd,1,'nnyp',nnyp1,ngrids1)
  ie=cio_i(iunhd,1,'nnzp',nnzp1,ngrids1)
  ie=cio_i_sca(iunhd,1,'npatch',npatch1,1)
  ie=cio_i_sca(iunhd,1,'nzg',nzg1,1)
  ie=cio_i_sca(iunhd,1,'nzs',nzs1,1)
  ie=cio_i_sca(iunhd,1,'ioutput',ioutput1,1)


  ! Find maximum size of any array on history file. Allocate scratch array of
  ! this size.

  maxarr=0
  do ngr=1,ngridsh
     maxarr=max(maxarr,nnxp1(ngr)*nnyp1(ngr)*nnzp1(ngr)  &
          ,nnxp1(ngr)*nnyp1(ngr)*nzg1*npatch1 &
          ,nnxp1(ngr)*nnyp1(ngr)*nzs1*npatch1)
  enddo

  ! read stuff here

  call hist_pol_read(maxarr,hnameinh,iunhd,ioutput1)

  print*,'back from read'
  close(iunhd)


  return
end subroutine init_conc_prev


!*******************************************************************************

subroutine hist_pol_read(maxarr,hnamein,iunhd,iout1)

  use an_header
  use var_tables

  implicit none

  integer :: maxarr,iunhd,iout1

  include 'interface.h'

  character (len=*) :: hnamein
  integer :: ngr,npts,nptsh,nv,nvh,i
  character type*1,post*10,fmt*3
  real, allocatable :: scr(:)
  integer :: inhunt=10

  type (head_table), allocatable,save :: hr_table(:)

  allocate (scr(maxarr))


  !  Read variable header info


  rewind(iunhd)

  read(iunhd,*) nvbtab
  allocate (hr_table(nvbtab))
  do nv=1,nvbtab
     read(iunhd,*)  hr_table(nv)%string   &
          ,hr_table(nv)%npointer  &
          ,hr_table(nv)%idim_type  &
          ,hr_table(nv)%ngrid  &
          ,hr_table(nv)%nvalues

  enddo


  ! Open history data file

  call rams_f_open(inhunt,hnamein,'UNFORMATTED','OLD','READ',0)

  ! Loop through all variables


  do nvh=1,nvbtab
     ! Read a variable
     nptsh=hr_table(nvh)%nvalues

     read(inhunt)(scr(i),i=1,nptsh)

     !  See if this variable is active in the current run
     ngr=hr_table(nvh)%ngrid
     if(ngr > nvgrids) cycle

     do nv = 1,num_var(ngr)
        npts=vtab_r(nv,ngr)%npts
        if(hr_table(nvh)%string == vtab_r(nv,ngr)%name) then
           if(nptsh /= npts) then
              print*,'Grid point number mismatch on history field:',  &
                   vtab_r(nv,ngr)%name,npts,nptsh
              stop 'History read number points error'
           endif

           if(vtab_r(nv,ngr)%name.eq.'PNO'.or.vtab_r(nv,ngr)%name.eq.'PNO2'.or. &
                vtab_r(nv,ngr)%name.eq.'PPM25'.or.vtab_r(nv,ngr)%name.eq.'PCO'.or. &
                vtab_r(nv,ngr)%name.eq.'PVOC'.or.vtab_r(nv,ngr)%name.eq.'PSO2'.or. &
                vtab_r(nv,ngr)%name.eq.'PSO4'.or.vtab_r(nv,ngr)%name.eq.'PAER'.or. &
                vtab_r(nv,ngr)%name.eq.'PVOC'.or.vtab_r(nv,ngr)%name.eq.'PSO2'.or. &
                vtab_r(nv,ngr)%name.eq.'PO3'.or.vtab_r(nv,ngr)%name.eq.'PRHCO'.or. &
                vtab_r(nv,ngr)%name.eq.'PHO2'.or.vtab_r(nv,ngr)%name.eq.'PO3P'.or. &
                vtab_r(nv,ngr)%name.eq.'PO1D'.or.vtab_r(nv,ngr)%name.eq.'PHO'.or. &
                vtab_r(nv,ngr)%name.eq.'PROO'                                 )then

              print 33,'Polutants History filling grid: ',ngr,nv,vtab_r(nv,ngr)%name,npts
33            format(a25,2i5,3x,a18,i10)
              !	 pause
              call atob(npts,scr(1),vtab_r(nv,ngr)%var_p)
              exit
           endif

        endif
     enddo

  enddo

  ! Close the input history file

  close(inhunt)

  deallocate(scr,hr_table)

  return
end subroutine hist_pol_read


