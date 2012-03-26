!--------------------------------------------------------------------
subroutine cup_dd_he(hes_cup,zd,hcd,z_cup,cdd,entr,jmin,ierr,mix,mgmxp,mkx     &
                    ,mgmzp,istart,iend,he,kdet,dby,he_cup)
  implicit none
  integer                         :: mix,mgmxp,mkx,mgmzp,istart,iend,i,k,ki
  integer, dimension(mgmxp)       :: jmin,ierr,kdet
  real,    dimension(mgmxp,mgmzp) :: zd,hcd,z_cup,cdd,he,dby,hes_cup,he_cup
  real                            :: entr,dz
  do k=2,mkx
     do i=istart,iend
        dby(i,k)=0.
        if(ierr(I).eq.0)then
           hcd(i,k)=hes_cup(i,k)
        endif
     enddo
  enddo

  do i=istart,iend
     if(ierr(I).eq.0)then
        k=jmin(i)
        hcd(i,k)=hes_cup(i,k)
        dby(i,k)=hcd(i,jmin(i))-hes_cup(i,k)

        do ki=jmin(i)-1,1,-1
           DZ=Z_cup(i,Ki+1)-Z_cup(i,Ki)
           HCD(i,Ki)=(HCD(i,Ki+1)*(1.-.5*CDD(i,Ki)*DZ) + entr*DZ*HE(i,Ki))/  &
                (1.+entr*DZ-.5*CDD(i,Ki)*DZ)

           dby(i,ki)=HCD(i,Ki)-hes_cup(i,ki)
        enddo
     endif
  enddo
  return
end subroutine cup_dd_he


!--------------------------------------------------------------------
subroutine cup_dd_moisture(j,zd,hcd,hes_cup,qcd,qes_cup,pwd,q_cup,z_cup,cdd    &
                          ,entr,jmin,ierr,gamma_cup,pwev,mix,mgmxp,mkx,mgmzp   &
                          ,istart,iend,bu,qrcd,q,he,hc,t_cup,iloop)
  use rconstants, only : alvl3
  implicit none
  integer                         :: mix,mgmxp,mkx,mgmzp,istart,iend,i,k,ki,j  &
                                    ,iloop
  integer, dimension(mgmxp)       :: jmin,ierr
  real                            :: entr, dz, dqeva, dh
  real,    dimension(mgmxp)       :: pwev, bu
  real,    dimension(mgmxp,mgmzp) :: zd,qcd,pwd,qrcd,hc,t_cup,hes_cup,hcd      &
                                    ,qes_cup,q_cup,z_cup,cdd,gamma_cup,q,he
  do i=istart,iend
     bu(i)=0.
     pwev(i)=0.
  enddo
  do k=1,mkx
     do i=istart,iend
        qcd(i,k)=0.
        qrcd(i,k)=0.
        pwd(i,k)=0.
     enddo
  enddo

  do i=istart,iend
     if(ierr(I).eq.0)then
        k=jmin(i)
        DZ=Z_cup(i,K+1)-Z_cup(i,K)
        qcd(i,k)=q_cup(i,k)
        qrcd(i,k)=qes_cup(i,k)
        pwd(i,jmin(i))=min(0.,qcd(i,k)-qrcd(i,k))
        pwev(i)=pwev(i)+pwd(i,jmin(i))
        qcd(i,k)=qes_cup(i,k)

        DH=HCD(I,k)-HES_cup(I,K)
        bu(i)=dz*dh

        do ki=jmin(i)-1,1,-1
           DZ=Z_cup(i,Ki+1)-Z_cup(i,Ki)
           QCD(i,Ki)=(qCD(i,Ki+1)*(1.-.5*CDD(i,Ki)*DZ) + entr*DZ*q(i,Ki))/  &
                (1.+entr*DZ-.5*CDD(i,Ki)*DZ)

           !--- To be negatively buoyant, hcd should be smaller than hes!

           DH=HCD(I,ki)-HES_cup(I,Ki)
           bu(i)=bu(i)+dz*dh
           QRCD(I,Ki)=qes_cup(i,ki)+(1./alvl3)*(GAMMA_cup(i,ki)/  &
                (1.+GAMMA_cup(i,ki)))*DH
           dqeva=qcd(i,ki)-qrcd(i,ki)
           if(dqeva.gt.0.) dqeva=0.
           pwd(i,ki)=zd(i,ki)*dqeva
           qcd(i,ki)=qrcd(i,ki)
           pwev(i)=pwev(i)+pwd(i,ki)
        enddo
        if(BU(I).ge.0.and.iloop.eq.1)then
           ierr(i)=7
        endif
     endif
  enddo

  return
end subroutine cup_dd_moisture


!--------------------------------------------------------------------
subroutine cup_dd_nms(zd,z_cup,cdd,entr,jmin,ierr,mix,mgmxp,mkx,mgmzp,istart   &
                     ,iend,itest,kdet,z1)
  implicit none
  integer                         :: mix,mgmxp,mkx,mgmzp,istart,iend,i,k,ki,itest
  integer, dimension(mgmxp)       :: jmin,ierr,kdet
  real                            :: entr,dz,a,perc
  real,    dimension(mgmxp)       :: z1
  real,    dimension(mgmxp,mgmzp) :: zd,z_cup,cdd
  
  perc=.03 ! Perc is the percentage of mass left when hitting the ground

  do k=1,mkx
     do i=istart,iend
        zd(i,k)=0.
        if (itest.eq.0) cdd(i,k)=0.
     enddo
  enddo
  a=1.-perc

  do i=istart,iend
     if(ierr(I).eq.0)then
        zd(i,jmin(i))=1.
        !--- Integrate downward, specify detrainment(cdd)!
        do ki=jmin(i)-1,1,-1
           DZ=Z_cup(i,Ki+1)-Z_cup(i,Ki)
           if(ki.le.kdet(i).and.itest.eq.0)then
              cdd(i,ki)=entr+(1.- (a*(z_cup(i,ki)-z1(i))+  &
                   perc*(z_cup(i,kdet(i))-z1(i)) )/        &
                   (a*(z_cup(i,ki+1)-z1(i))+               &
                   perc*(z_cup(i,kdet(i))-z1(i))))/dz
           endif
           zd(i,ki)=zd(i,ki+1)*(1.+(entr-cdd(i,ki))*dz)
        enddo

     endif
  enddo
  return
end subroutine cup_dd_nms


!--------------------------------------------------------------------
subroutine cup_dd_aa0(edt,ierr,aa0,jmin,gamma_cup,t_cup,hcd,hes_cup,z,mix      &
                     ,mgmxp,mkx,mgmzp,istart,iend,zd)
  use rconstants, only : grav, cpdry
  implicit none
  integer                      :: i,k,kk,mix,mgmxp,mkx,mgmzp,istart,iend
  integer, dimension(mgmxp)    :: jmin,ierr
  real                         :: dz
  real, dimension(mgmxp)       :: aa0,edt
  real, dimension(mgmxp,mgmzp) :: gamma_cup,t_cup,qes_cup,z,hes_cup,zd,hcd
  do K=1,MKX-1
     do I=ISTART,IEND
        if (ierr(I).eq.0.and.k.lt.jmin(i)) then
           KK=JMIN(I)-K

           DZ=(Z(I,KK)-Z(I,KK+1))
           AA0(I)=AA0(I)+zd(i,kk)*EDT(I)*DZ*(grav/(cpdry*T_cup(I,KK)))* &
                ((hcd(i,kk)-hes_cup(i,kk))/(1.+GAMMA_cup(i,kk)))
        endif
     enddo
  enddo
  return
end subroutine cup_dd_aa0


!--------------------------------------------------------------------
subroutine cup_dd_edt(ierr,us,vs,z,ktop,kbcon,edt,p,pwav,pwev,mix,mgmxp,mkx    &
                     ,mgmzp,istart,iend,edtmax,edtmin,maxens2,edtc,vshear,sdp  &
                     ,vws)
  implicit none
  integer                           :: i,kk,k,mix,mgmxp,mkx,mgmzp,istart,iend  &
                                      ,maxens2
  integer, dimension(mgmxp)         :: ktop,kbcon,ierr
  real                              :: pef,edtmax,edtmin,einc,pefb,prezk,zkbc
  real,    dimension(mgmxp)         :: pwav,pwev,edt,vshear,sdp,vws
  real,    dimension(mgmxp,mgmzp)   :: us,vs,z,p
  real,    dimension(mgmxp,maxens2) :: edtc

  !--- DETERMINE DOWNDRAFT STRENGTH IN TERMS OF WINDSHEAR

  ! */ calculate an average wind shear over the depth of the cloud

  do i=ISTART,IEND
     edt(i)=0.
     vws(i)=0.
     sdp(i)=0.
     vshear(i)=0.
  enddo
  do kk = 1,mkx-1
     do I=ISTART,IEND
        if(ierr(i).ne.0) cycle
        if (kk .le. min0(ktop(i),mkx-1) .and. kk .ge. kbcon(i)) then
           vws(i) = vws(i)+(abs((us(i,kk+1)-us(i,kk))/(z(i,kk+1)-z(i,kk))) +  &
                abs((vs(i,kk+1)-vs(i,kk))/(z(i,kk+1)-z(i,kk)))) * &
                (p(i,kk) - p(i,kk+1))
           sdp(i) = sdp(i) + p(i,kk) - p(i,kk+1)
        endif
        if (kk .eq. mkx-1)  vshear(i) = 1.e3 * vws(i) / sdp(i)
     enddo
  enddo
  do I=ISTART,IEND
     if (ierr(i).eq.0) then
        pef=(1.591-.639*VSHEAR(I)+.0953*(VSHEAR(I)**2) - .00496*(VSHEAR(I)**3))
        if (pef.gt.edtmax) pef=edtmax
        if (pef.lt.edtmin) pef=edtmin

        !--- cloud base precip efficiency

        zkbc=z(i,kbcon(i))*3.281e-3
        prezk=.02
        if (zkbc.gt.3.) then
           prezk = .96729352+zkbc*(-.70034167+zkbc*(.162179896+zkbc*  &
                (-1.2569798E-2+zkbc*(4.2772E-4-zkbc*5.44E-6))))
        endif
        if (zkbc.gt.25) then
           prezk=2.4
        endif
        pefb = 1./(1.+prezk)
        if (pefb.gt.edtmax) pefb=edtmax
        if (pefb.lt.edtmin) pefb=edtmin
        EDT(I) = 1.-.5*(pefb+pef)
        einc = edt(i)/float(maxens2+1)
        do k=1,maxens2
           edtc(i,k)=edt(i)+float(k-maxens2/2-1)*einc
           edtc(i,k)=edt(i)-float(k)*einc
        enddo
     endif
  enddo
  do I=ISTART,IEND
     if (ierr(i).eq.0) then
        do k=1,maxens2
           EDTC(I,K)=-EDTC(I,K)*PWAV(I)/PWEV(I)
           if (EDTC(I,K).gt.edtmax) EDTC(I,K)=edtmax
           if (EDTC(I,K).lt.edtmin) EDTC(I,K)=edtmin
        enddo
     endif
  enddo

  return
end subroutine cup_dd_edt
