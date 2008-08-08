!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

subroutine ozone(mzp,mxp,myp,ia,iz,ja,jz,ng,deltat)

  use mem_grid
  use mem_basic
  use mem_gaspart
  use mem_radiate
  use rconstants
  use mem_tend
  use var_tables
  implicit none


  real :: deltat
  integer :: mzp,mxp,myp,ia,iz,ja,jz,ng

  !  IF (MOD(time + .001,chemfrq) .LT. dtlt .OR. time .LT. 0.001) THEN

  !     PRINT 90,time,time/3600.+(itime1/100+MOD(itime1,100)/60.)
  !90   FORMAT('   chemical reactions tendencies updated    time =',f10.1,  &
  !          '  UTC TIME (HRS) =',F6.1)

  call chemistry (mzp, mxp, myp, ia, iz, ja, jz,                         &
       gaspart_g(ng)%pno(1,1,1)   ,gaspart_g(ng)%pno2(1,1,1),      &
       gaspart_g(ng)%pco(1,1,1)   ,gaspart_g(ng)%pvoc(1,1,1),      &
       gaspart_g(ng)%po3(1,1,1)   ,gaspart_g(ng)%pso2(1,1,1),      &
       gaspart_g(ng)%pso4(1,1,1)  ,                              &
       gaspart_g(ng)%prhco(1,1,1) ,gaspart_g(ng)%pho2(1,1,1),      &
       gaspart_g(ng)%po3p(1,1,1)  ,gaspart_g(ng)%po1d(1,1,1),      &
       gaspart_g(ng)%pho(1,1,1)   ,gaspart_g(ng)%proo(1,1,1),      &
       basic_g(ng)%theta(1,1,1) ,basic_g(ng)%dn0(1,1,1),       &
       basic_g(ng)%pi0(1,1,1)   ,basic_g(ng)%pp(1,1,1),       &
       basic_g(ng)%rv(1,1,1),                                  &
       radiate_g(ng)%rshort(1,1),radiate_g(ng)%cosz(1,1),      &
       grid_g(ng)%rtgt(1,1),grid_g(ng)%topma (1,1),deltat,cpi,cpor,p00,zt, &
       gaspart_g(ng)%pnot(1)   ,gaspart_g(ng)%pno2t(1),           &
       gaspart_g(ng)%pcot(1)   ,gaspart_g(ng)%pvoct(1),           &
       gaspart_g(ng)%po3t(1)   ,gaspart_g(ng)%pso2t(1),           &
       gaspart_g(ng)%pso4t(1)  ,                                &
       gaspart_g(ng)%prhcot(1) ,gaspart_g(ng)%pho2t(1),           &
       gaspart_g(ng)%po3pt(1)  ,gaspart_g(ng)%po1dt(1),           &
       gaspart_g(ng)%phot(1)   ,gaspart_g(ng)%proot(1))   


  !endif

  return
end subroutine ozone


!##############################################################################
subroutine chemistry(  m1,  m2,   m3, ia, iz, ja, jz,       &
     noi ,        no2i ,          coi ,   & 
     vocsi ,         o3i ,  so2i, so4i, rchoi ,   &
     ho2i ,        O3Pi ,         O1Di ,   &
     HOi ,        RO2i ,                  & 
     theta,    dn0    ,    pi0  ,   pp   ,   &
     rv  , rshort    ,   cosz  ,   rtgt  ,   &
     topt,   dtlt    ,      cpi,   cpor  ,   &
     p00,     zt,                                   &
     not  ,        no2t  ,          cot , & 
     vocst  ,         o3t  , so2t, so4t,rchot  ,&
     ho2t  ,        O3Pt  ,         O1Dt  ,&
     HOt  ,        RO2t)

  use ozone_const

  implicit none

  integer :: m1, m2, m3, ia, iz, ja, jz, i, j, ik,jk

  !variables used to keep time n-1 values

  real  pno,pno2,po3p,po1d,po3,ph2o,pco,pho,pro2,prcho,pvocs,pho2,pso2,pso4,hom

  !chemichal elements and local altitude

  real, dimension(m1,m2,m3) :: alt,   &
       not,no2t,so2t,so4t,cot,vocst,o3t,rchot,ho2t,O3Pt,O1Dt,HOt,RO2t,          &
       noi,no2i,so2i,so4i,coi,vocsi,o3i,rchoi,ho2i,O3Pi,O1Di,HOi,RO2i,          &
       pjno2,pjo3

  real, dimension(m1) :: no, no2, co, vocs, o3, rcho, ho2, O3P, O1D, HO, RO2,h2o,so2,so4

  ! 3-D model's variables used
  real, dimension(m1,m2,m3) :: tempk, dn0,  pi0,  rv,theta,pp,ppi

  ! 2-D model's variables used
  real, dimension(m2,m3) :: rshort,cosz,cosz2,rtgt,topt

  ! 1-D model's variables used
  real, dimension (m1)   :: zt

  ! model's "constants" used

  real                   :: cpi,cpor,p00,avo,temp

  !velocities coefficients

  real, dimension (15)   :: k

  !photolysis rates

  real                   :: j2,j6,kso2

  !production and lost terms used in the implicit scheme

  real prod,plosc

  !time-split determining variables

  real dtll_factor,dtlt,dtll,dtlh,dtlh_factor,dtlti

  integer niter_ozo,inter,niter_ho2,inter2,zera

  !avogadro constant

  AVO=6.0221367E+23

  !space loop

  do j=ja,jz
     do i=ia,iz

        !storing previous values for all gases

        do ik=2,m1

           !calculating absolute temperature using Exner function and Theta

           tempk(ik,i,j)=theta(ik,i,j)*pi0(ik,i,j)*cpi

           !calculating pressure

           ppi(ik,i,j)=((pp(ik,i,j)+pi0(ik,i,j))*cpi)**cpor*p00

           !calculationg the height of each grid point above sea level

           alt(ik,i,j)=zt(ik)*rtgt(i,j)   + topt(i,j)	 

        enddo
     enddo
  enddo

  !calculating time-split for chemical reactions (minimum 1 s)

  niter_ozo = max(1,nint(dtlt/5.e-1))
  dtll_factor = 1. / float(niter_ozo)
  dtll = dtlt * dtll_factor

  do j=ja,jz

     do i=ia,iz


        do ik=1,m1
           no  (ik)=max(0.,noi  (ik,i,j))
           no2 (ik)=max(0.,no2i  (ik,i,j))
           so2 (ik)=max(0.,so2i  (ik,i,j))
           so4 (ik)=max(0.,so4i  (ik,i,j))
           co  (ik)=max(0.,coi  (ik,i,j))
           vocs(ik)=max(0.,vocsi(ik,i,j))
           o3  (ik)=max(0.,o3i  (ik,i,j))
           rcho(ik)=max(0.,rchoi(ik,i,j))
           ho2 (ik)=max(0.,ho2i (ik,i,j))
           o3p (ik)=max(0.,o3pi (ik,i,j))
           o1d (ik)=max(0.,o1di (ik,i,j))
           ho  (ik)=max(0.,hoi  (ik,i,j))
           ro2 (ik)=max(0.,ro2i (ik,i,j))

           h2o(ik) =rv(ik,i,j)

        enddo


        call conv_ppm_rm(m1,o3  (1),48.0 ,28.97,1.e6)
        call conv_ppm_rm(m1,co  (1),28.0 ,28.97,1.e6)
        call conv_ppm_rm(m1,no  (1),30.0 ,28.97,1.e6)
        call conv_ppm_rm(m1,no2 (1),46.0 ,28.97,1.e6)
        call conv_ppm_rm(m1,vocs(1),42.08,28.97,1.e6)
        call conv_ppm_rm(m1,rcho(1),44.0 ,28.97,1.e6)
        call conv_ppm_rm(m1,ho2 (1),33.0 ,28.97,1.e6)
        call conv_ppm_rm(m1,h2o (1),18.0 ,28.97,1.e6)
        call conv_ppm_rm(m1,O3P (1),16.0 ,28.97,1.e6)
        call conv_ppm_rm(m1,O1D (1),16.0 ,28.97,1.e6)
        call conv_ppm_rm(m1,HO  (1),17.0 ,28.97,1.e6)
        call conv_ppm_rm(m1,RO2 (1),47.0 ,28.97,1.e6)


        j2=max(0.001,42.92e-02*cosz(i,j))/60.
        j6=max(0.0000001,2.04E-03*cosz(i,j))/60.

        do ik=2,m1-1

           !calculating k's coefficients

           do jk=1,15
              k(jk)=(A(jk)*((tempk(ik,i,j)/Tref)**B(jk))*exp(-Ea(jk)/(Rcal*tempk(ik,i,j))))/60.
           enddo
           !initiate time-split for reactions in order to avoid numerical instalibity

           do inter=1,niter_ozo

              !keeping values from time n-1
              pno   =no  (ik)
              pno2  =no2 (ik)
              po3p  =o3p (ik)
              po1d  =o1d (ik)
              po3   =o3  (ik) 
              ph2o  =h2o (ik)
              pco   =co  (ik)
              pho   =ho  (ik)
              pro2  =ro2 (ik)
              prcho =rcho(ik)
              pvocs =vocs(ik)
              pho2  =ho2 (ik)  

              !updating concentrations



              !implicit integration of O3P

              o3p(ik) = (po3p+(j2*pno2*dtll))/(1.+(M*k(3)*(o2)*dtll))

              !implicit integration of O1D

              O1D(ik) = (pO1D+(j6* po3*dtll))/(1.+(k(7)*ph2o*dtll))

              !implicit integration of HO

              plosc=(k(8)*pco)+(k(11)*pvocs)+(k(13)*prcho)

              prod=(2.*k(7)*ph2o*pO1D)+(k(9)*pno*pho2)+(k(10)*pho2*po3)

              ho(ik)       = (pho+(prod*dtll))/(1.+(plosc*dtll))


              !implicit integration of ro2

              plosc=(k(12)*pno)+(k(15)* pho2)
              prod=(k(11)*pvocs*pho)

              ro2(ik)   = (pro2+(prod*dtll))/(1.+(plosc*dtll))


              !implicit integration of o3


              plosc=  (k(4)  *  pno)+(k(5)  * pno2)+j6+(k(10) * pho2)
              prod=k(3)  * o2   *   M    *  po3p

              o3(ik)    =   (po3 + dtll*prod)/(1.+dtll*plosc)

              !implicit integration of no

              plosc=(2.*k(1)*pno*o2) +(k(4)*po3)+(k(9)*pho2)+(k(12) * pro2)
              prod=j2*pno2

              no(ik)    =   (pno + dtll*prod)/(1.+dtll*plosc)

              !implicit integration of no2

              plosc=j2+(k(5)  * po3)
              prod=(2.*k(1)*pno*pno*o2)+(k(4)*pno*po3)+(k(9)*pno*pho2)+(k(12) * pro2 * pno)

              no2(ik)    =   (pno2 + dtll*prod)/(1.+dtll*plosc)

              !implicit integration of vocs

              plosc=  k(11)  * pho 
              prod=0.
              vocs(ik)   =   (pvocs + dtll*prod)/(1.+dtll*plosc) 

              !implicit integration of co

              plosc=k(8)  *  pho
              prod=0.
              co(ik)   =     (pco + dtll*prod)/(1.+dtll*plosc) 

              !implicit integration of rcho

              plosc=(k(13) * pho)
              prod=k(12) * pro2 * pno
              rcho(ik)    =   (prcho + dtll*prod)/(1.+dtll*plosc)

              !implicit integration of ho2


              plosc=(k(9)*pno)+(k(10)*po3)+(2.*k(14)*pho2)+(k(15)*pro2)
              prod=(k(8)*pho*pco)+(k(12)*pro2*pno)

              ho2(ik)    =   (pho2 + dtll*prod)/(1.+dtll*plosc)


           enddo !time
        enddo !levels

        !converting units from ppm to kg/kg

        call conv_ppm_rm(m1,o3  (1),28.97,48.0 ,1.e-6)
        call conv_ppm_rm(m1,co  (1),28.97,28.0 ,1.e-6)
        call conv_ppm_rm(m1,no  (1),28.97,30.0 ,1.e-6)
        call conv_ppm_rm(m1,no2 (1),28.97,46.0 ,1.e-6)
        call conv_ppm_rm(m1,vocs(1),28.97,42.08,1.e-6)
        call conv_ppm_rm(m1,rcho(1),28.97,44.0 ,1.e-6)
        call conv_ppm_rm(m1,ho2 (1),28.97,33.0 ,1.e-6)
        call conv_ppm_rm(m1,h2o (1),28.97,18.0 ,1.e-6)
        call conv_ppm_rm(m1,O3P (1),28.97,16.0 ,1.e-6)
        call conv_ppm_rm(m1,O1D (1),28.97,16.0 ,1.e-6)
        call conv_ppm_rm(m1,HO  (1),28.97,17.0 ,1.e-6)
        call conv_ppm_rm(m1,RO2 (1),28.97,47.0 ,1.e-6)

        dtlti=1./dtlt

        ! integration of so4

        do ik=1,m1

           !defining velocity constant for SO2

           kso2=1e-18

           !converting variables to units of s-1

           pso2=(so2(ik)*dn0(ik,i,j)*avo)/(1e-3*64.)
           pso4=(so4(ik)*dn0(ik,i,j)*avo)/(1E-3*96.)

           hom=(ho(ik)*dn0(ik,i,j)*avo)/(1E-3*17.)

           temp=pso4+(kso2*hom*pso2)*dtlt 

           !converting back to kg/kg

           so4(ik)=(temp*(1e-3)*96.)/(dn0(ik,i,j)*avo) 

        enddo

        do ik=1,m1

           so4t (ik,i,j)=so4t (ik,i,j)+((so4 (ik)-so4i (ik,i,j))*dtlti)
           o3t  (ik,i,j)=o3t  (ik,i,j)+((o3  (ik)-o3i  (ik,i,j))*dtlti)
           cot  (ik,i,j)=cot  (ik,i,j)+((co  (ik)-coi  (ik,i,j))*dtlti)
           not  (ik,i,j)=not  (ik,i,j)+((no  (ik)-noi  (ik,i,j))*dtlti)
           no2t (ik,i,j)=no2t (ik,i,j)+((no2 (ik)-no2i (ik,i,j))*dtlti)
           vocst(ik,i,j)=vocst(ik,i,j)+((vocs(ik)-vocsi(ik,i,j))*dtlti)
           rchot(ik,i,j)=rchot(ik,i,j)+((rcho(ik)-rchoi(ik,i,j))*dtlti)
           ho2t (ik,i,j)=ho2t (ik,i,j)+((ho2 (ik)-ho2i (ik,i,j))*dtlti)
           o3pt (ik,i,j)=o3pt (ik,i,j)+((o3p (ik)-o3pi (ik,i,j))*dtlti)
           o1dt (ik,i,j)=o1dt (ik,i,j)+((o1d (ik)-o1di (ik,i,j))*dtlti)
           hot  (ik,i,j)=hot  (ik,i,j)+((ho  (ik)-hoi  (ik,i,j))*dtlti)
           ro2t (ik,i,j)=ro2t (ik,i,j)+((ro2 (ik)-ro2i (ik,i,j))*dtlti)

        enddo

     enddo
  enddo



  return
end subroutine chemistry

subroutine conv_ppm_rm(n1,b,pden,pnum,expo)

  implicit none

  integer n1,k

  real,dimension(n1) :: b

  real pden,pnum,expo,fator


  fator=(pden/pnum)*expo

  do k=1,n1-1
     b(k)=b(k)*fator
  enddo

  return
end subroutine conv_ppm_rm


subroutine posidef(n1,n2,n3,ia,iz,ja,jz,ng,ich,iso)

  use mem_gaspart

  implicit none

  integer :: n1,n2,n3,ia,iz,ja,jz,ng,ich,iso

  call zeragas(n1,n2,n3,ia,iz,ja,jz,gaspart_g(ng),ich,iso)

  return
end subroutine posidef

subroutine zeragas(n1,n2,n3,ia,iz,ja,jz,gaspart,ich,iso)

  use mem_gaspart

  implicit none

  type (gaspart_vars) gaspart

  integer :: n1,n2,n3,ia,iz,ja,jz,i,j,k,ich,iso

  do i=ia,iz

     do j=ja,jz

        do k=1,n1

           if(iso==1)then

              gaspart%pno  (k,i,j)=max(0.,gaspart%pno  (k,i,j))
              gaspart%pno2 (k,i,j)=max(0.,gaspart%pno2 (k,i,j))
              gaspart%pso2 (k,i,j)=max(0.,gaspart%pso2 (k,i,j))
              gaspart%pco  (k,i,j)=max(0.,gaspart%pco  (k,i,j))
              gaspart%ppm25(k,i,j)=max(0.,gaspart%ppm25(k,i,j))
              gaspart%pvoc (k,i,j)=max(0.,gaspart%pvoc (k,i,j))

              if(ich==1)then

                 gaspart%po3  (k,i,j)=max(0.,gaspart%po3  (k,i,j))
                 gaspart%prhco(k,i,j)=max(0.,gaspart%prhco(k,i,j))
                 gaspart%pho2 (k,i,j)=max(0.,gaspart%pho2 (k,i,j))
                 gaspart%po3p (k,i,j)=max(0.,gaspart%po3p (k,i,j))
                 gaspart%po1d (k,i,j)=max(0.,gaspart%po1d (k,i,j))
                 gaspart%pho  (k,i,j)=max(0.,gaspart%pho  (k,i,j))
                 gaspart%proo (k,i,j)=max(0.,gaspart%proo (k,i,j))

              endif

           endif
        enddo
     enddo
  enddo

  return
end subroutine zeragas
