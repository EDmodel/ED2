!----------------------------------------------------------------------
subroutine cup_up_he(k22,hkb,z_cup,cd,entr,he_cup,hc,mix,mgmxp,mkx,mgmzp,kbcon &
                    ,ierr,istart,iend,dby,he,hes_cup)

  implicit none
  ! hc = cloud moist static energy
  ! hkb = moist static energy at originating level
  ! he = moist static energy on model levels
  ! he_cup = moist static energy on model cloud levels
  ! hes_cup = saturation moist static energy on model cloud levels
  ! dby = buoancy term
  ! cd= detrainment function 
  ! z_cup = heights of model cloud levels 
  ! entr = entrainment rate
  integer                         :: i,j,k,mix,mgmxp,mkx,mgmzp,istart,iend
  integer, dimension(mgmxp)       :: kbcon,ierr,k22
  real                            :: entr, dz
  real,    dimension(mgmxp)       :: hkb
  real,    dimension(mgmxp,mgmzp) :: he_cup,hc,z_cup,cd,dby,he,hes_cup

  !--- Moist static energy inside cloud

  do i=istart,iend
     if (ierr(i).eq.0.) then
        hkb(i)=he_cup(i,k22(i))
        do k=1,k22(i)
           hc(i,k)=he_cup(i,k)
           DBY(I,K)=0.
        enddo
        do k=k22(i),kbcon(i)-1
           hc(i,k)=hkb(i)
           DBY(I,K)=0.
        enddo
        k=kbcon(i)
        hc(i,k)=hkb(i)
        DBY(I,Kbcon(i))=Hkb(I)-HES_cup(I,K)
     endif
  enddo
  do k=2,mkx-1
     do i=istart,iend
        if (k.gt.kbcon(i).and.ierr(i).eq.0.) then
           DZ=Z_cup(i,K)-Z_cup(i,K-1)
           HC(i,K)=(HC(i,K-1)*(1.-.5*CD(i,K)*DZ)+entr*      &
                DZ*HE(i,K-1))/(1.+entr*DZ-.5*cd(i,k)*dz)
           DBY(I,K)=HC(I,K)-HES_cup(I,K)
        endif
     enddo

  enddo
  return
end subroutine cup_up_he


!----------------------------------------------------------------------
subroutine cup_up_moisture(ierr,z_cup,qc,qrc,pw,pwav,kbcon,ktop,mix,mgmxp,mkx  &
                          ,mgmzp,istart,iend,cd,dby,mentr_rate,q,GAMMA_cup,zu  &
                          ,qes_cup,k22,qe_cup)
  use rconstants, only : alvl
  implicit none
  ! cd= detrainment function 
  ! q = environmental q on model levels
  ! qe_cup = environmental q on model cloud levels
  ! qes_cup = saturation q on model cloud levels
  ! dby = buoancy term
  ! cd= detrainment function 
  ! zu = normalized updraft mass flux
  ! gamma_cup = gamma on model cloud levels
  ! mentr_rate = entrainment rate
  ! qc = cloud q (including liquid water) after entrainment
  ! qrch = saturation q in cloud
  ! qrc = liquid water content in cloud after rainout
  ! pw = condensate that will fall out at that level
  ! pwav = totan normalized integrated condensate (I1)
  ! c0 = conversion rate (cloud to rain)
  integer                         :: istart,iend,mix,mgmxp,mkx,mgmzp,i,k,iall
  integer, dimension(mgmxp)       :: kbcon,ktop,ierr,k22
  real                            :: radius,dz,dh,qrch,c0,mentr_rate
  real,    dimension(mgmxp)       :: pwav
  real,    dimension(mgmxp,mgmzp) :: q,zu,GAMMA_cup,qe_cup,dby,cd,z_cup        &
                                    ,qes_cup,qc,qrc,pw
 
  iall=0
  c0=.002
 
  !--- No precip for small clouds

  if(mentr_rate.gt.0.)then
     radius=.2/mentr_rate
     if(radius.lt.900.)c0=0.
  endif
  do i=istart,iend
     pwav(i)=0.
  enddo
  do k=1,mkx
     do i=istart,iend
        pw(i,k) =0.
        qc(i,k) =qe_cup(i,k)
        qrc(i,k)=0.
     enddo
  enddo
  do i=istart,iend
     if(ierr(i).eq.0.)then
        do k=k22(i),kbcon(i)-1
           qc(i,k)=qe_cup(i,k22(i))
        enddo
     endif
  enddo

  do K=2,MKX-1
     do I=ISTART,IEND
        if (ierr(i).ne.0)  cycle
        if (K.lt.KBCON(I)) cycle
        if (K.gt.KTOP(I))  cycle
        DZ=Z_cup(i,K)-Z_cup(i,K-1)

        !--- 1. Steady state plume equation, for what could
        !---    be in cloud without condensation
        QC(i,K)=(QC(i,K-1)*(1.-.5*CD(i,K)*DZ)+mentr_rate*        &
             DZ*Q(i,K-1))/(1.+mentr_rate*DZ-.5*cd(i,k)*dz)

        !--- Saturation  in cloud, this is what is allowed to be in it
        QRCH=QES_cup(I,K)+(1./alvl)*(GAMMA_cup(i,k)/(1.+GAMMA_cup(i,k)))*DBY(I,K)

        !--- Liquid water content in cloud after rainout
        QRC(I,K)=(QC(I,K)-QRCH)/(1.+C0*DZ)
        if (qrc(i,k).lt.0.) then
           qrc(i,k)=0.
        endif
        !
        !--- 3.Condensation
        !
        PW(i,k)=c0*dz*QRC(I,K)*zu(i,k)  !unit: kg[liq water]/kg[air]
                                        !unit of c0 is m^-1
        if (iall.eq.1) then
           qrc(i,k)=0.
           pw(i,k)=(QC(I,K)-QRCH)*zu(i,k)
           if (pw(i,k).lt.0.) pw(i,k)=0.
        endif
        !
        !--- Set next level
        !
        QC(I,K)=QRC(I,K)+qrch
        !
        !--- Integrated normalized ondensate
        !
        PWAV(I)=PWAV(I)+PW(I,K)
     enddo
  enddo
  return
end subroutine cup_up_moisture


!----------------------------------------------------------------------
subroutine cup_up_nms(zu,z_cup,entr,cd,kbcon,ktop,mix,mgmxp,mkx,mgmzp,istart   &
                     ,iend,ierr,k22)

  implicit none
  integer                         :: i,k,mix,mgmxp,mkx,mgmzp,istart,iend
  integer, dimension(mgmxp)       :: kbcon,ktop,k22,ierr
  real,    dimension(mgmxp,mgmzp) :: zu,z_cup,cd
  real entr, dz
  do k=1,mkx
     do i=istart,iend
        zu(i,k)=0.
     enddo
  enddo
  do i=istart,iend
     if (ierr(I).eq.0) then
        do k=k22(i),kbcon(i)
           zu(i,k)=1.
        enddo
        do K=KBcon(i)+1,KTOP(i)
           DZ=Z_cup(i,K)-Z_cup(i,K-1)
           ZU(i,K)=ZU(i,K-1)*(1.+(entr-cd(i,k))*DZ)
        enddo
     endif
  enddo
  return
end subroutine cup_up_nms


!----------------------------------------------------------------------
subroutine cup_up_aa0(aa0,z,zu,dby,GAMMA_CUP,t_cup,kbcon,ktop,mix,mgmxp,mkx    &
                     ,mgmzp,istart,iend,ierr)
  use rconstants, only : g, cp
  implicit none
  integer                         :: i,k,mix,mgmxp,mkx,mgmzp,istart,iend
  integer, dimension(mgmxp)       :: kbcon,ktop,ierr
  real                            :: dz, da
  real,    dimension(mgmxp)       :: aa0 
  real,    dimension(mgmxp,mgmzp) :: z,zu,gamma_cup,t_cup,dby
  do I=ISTART,IEND
     aa0(i)=0.
  enddo
  do K=2,MKX-1
     do I=ISTART,IEND
        if (ierr(i).ne.0)  cycle
        if (K.le.KBCON(I)) cycle
        if (K.gt.KTOP(I))  cycle
        DZ = Z(I,K)-Z(I,K-1)
        da = zu(i,k)*DZ*(g/(cp*((T_cup(I,K)))))*DBY(I,K-1)/  &
             (1.+GAMMA_CUP(I,K))
        if (K.eq.KTOP(I).and.da.le.0.) cycle
        AA0(I)=AA0(I)+da
        if (aa0(i).lt.0.) aa0(i)=0.
     enddo
  enddo

  return
end subroutine cup_up_aa0
