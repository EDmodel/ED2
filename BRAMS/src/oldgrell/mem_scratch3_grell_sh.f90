! Module necessary to Grell Cumulus param.
! Scratch variables - module 3

module mem_scratch3_grell_sh

  !TYPE scratch3_grell_vars
  ! adapted in july-15-2002 for 5.x version
  ! 3d dependence (mgmxp,mgmzp,maxens2)
  real, allocatable, dimension(:,:,:) :: dellat_ens,  &
       dellaq_ens,                                    &
       dellaqc_ens,                                   &
       pwo_ens
  ! 3d dependence (mgmxp,mgmyp,ensdim)
  real, allocatable, dimension(:,:,:) :: xf,  &
       xf_ens,                                &
       pr_ens,                                &
       outt_ens

  ! 2d dependence (mgmxp,mgmzp)
  real, allocatable, dimension(:,:) :: HE,      &
       HES,                                     &
       QES,                                     &
       Z,                                       &
       TV,                                      &
       DBY,                                     &
       QC,                                      &
       QRCD,                                    &
       PWD,                                     &
       PW
  real, allocatable, dimension(:,:) :: HEO,     &
       HESO,                                    &
       QESO,                                    &
       ZO,                                      &
       TVO,                                     &
       DBYO,                                    &
       QCO,                                     &
       QRCDO,                                   &
       PWDO,                                    &
       PWO
  real, allocatable, dimension(:,:) :: XHE,     &
       XHES,                                    &
       XQES,                                    &
       XZ,                                      &
       XTV,                                     &
       XT_Grell,                                & ! Substitui XT original
       XQ,                                      &
       XDBY,                                    &
       XQC,                                     &
       XQRCD,                                   &
       XPWD,                                    &
       XPW
  real, allocatable, dimension(:,:) :: hcd,     &
       hcdo,                                    &
       xhcd
  real, allocatable, dimension(:,:) :: qcd,     &
       qcdo,                                    &
       xqcd
  real, allocatable, dimension(:,:) :: dbyd,    &
       dbydo
  real, allocatable, dimension(:,:) :: hc,      &
       hco,                                     &
       xhc,                                     &
       qrc,                                     &
       qrco,                                    &
       xqrc,                                    &
       zu,                                      &
       zuo,                                     &
       xzu,                                     &
       zd,                                      &
       zdo,                                     &
       xzd
  real, allocatable, dimension(:,:) :: DELLAH,  &
       DELLAQ,                                  &
       DELLAT,                                  &
       DELLAQC
  real, allocatable, dimension(:,:) :: qes_cup,  &
       q_cup,                                    &
       he_cup,                                   &
       hes_cup,                                  &
       z_cup,                                    &
       p_cup,                                    &
       gamma_cup,                                &
       t_cup
  real, allocatable, dimension(:,:) :: qeso_cup,  &
       qo_cup,                                    &
       heo_cup,                                   &
       heso_cup,                                  &
       zo_cup,                                    &
       po_cup,                                    &
       gammao_cup,                                &
       tn_cup
  real, allocatable, dimension(:,:) :: xqes_cup,  &
       xq_cup,                                    &
       xhe_cup,                                   &
       xhes_cup,                                  &
       xz_cup,                                    &
       xt_cup, &
       xt
  real, allocatable, dimension(:,:) :: cd,  &
       cdd,                                 &
       scr1
  ! 2d dependence (mgmxp,maxens)
  real, allocatable, dimension(:,:) :: xaa0_ens
  ! 2d dependence (mgmxp,maxens2)
  real, allocatable, dimension(:,:) :: edtc
  ! 2d dependence (mgmxp,mgmyp)
  real, allocatable, dimension(:,:) :: cwf,  &
       pwf,                                  &
       pwdf,                                 &
       eddt,                                 &
       predb,                                &
       xmass

  ! 1d dependence (mgmxp)
  integer, allocatable, dimension(:) :: kzdown,  &
       KBMAX,                                    &
       IERR,                                     &
       K22,                                      &
       KBCON,                                    &
       KB,                                       &
       JMIN,                                     &
       KTOP,                                     &
       kstabi,                                   &
       kstabm,                                   &
       K22x,                                     &
       KBCONx,                                   &
       KBx,                                      &
       KTOPx, &
       kzi
  real, allocatable, dimension(:) :: EDT,  &
       EDTO,                               &
       EDTX,                               &
       AA1,                                &
       AA0,                                &
       XAA0,                               &
       HKB,                                & 
       HKBO,                               &
       aad,                                &
       XHKB,                               &
       QKB,                                &
       QKBO,                               &
       XMB,                                &
       !for CATT
       tkemax
  real, allocatable, dimension(:) :: XPWAV,  &
       XPWEV,                                &
       PWAV,                                 &
       PWEV,                                 &
       PWAVO,                                &
       PWEVO,                                &
       BU,                                   &
       BUO
  ! 1d dependence (maxens3)
  real, allocatable, dimension(:) :: xff_ens3
  ! 1d dependence (maxens)
  real, allocatable, dimension(:) :: xk
  ! 1d dependence (mgmxp)
  real, allocatable, dimension(:) :: xfac1
  real, allocatable, dimension(:) :: cap_max,cap_max_increment
  ! 1d dependence (maxens)
  real, allocatable, dimension(:) :: mbdt_ens
  ! 1d dependence (maxens2)
  real, allocatable, dimension(:) :: edt_ens
  ! 1d dependence (mgmxp)
  !srf variaveis para a rotine "cup_dd_edt"
  real, allocatable, dimension(:) :: vshear,   &
       sdp,                                    &
       vws

  !END TYPE scratch3_grell_vars

contains

  subroutine alloc_scratch3_grell_sh  !(scratch3_grell)

    use mem_grell_param, only : mgmxp,  & ! INTENT(IN)
         mgmyp,                         & ! INTENT(IN)
         mgmzp,                         & ! INTENT(IN)
         maxens,                        & ! INTENT(IN)
         maxens2,                       & ! INTENT(IN)
         maxens3,                       & ! INTENT(IN)
         ensdim                           ! INTENT(IN)

    implicit none
    !TYPE (scratch3_grell_vars) :: scratch3_grell

    integer :: i

    !write(6,'(a1,78a1)') ' ',('_',i=1,78)
    !print*,'In CUPENSS - Alocando memoria'
    !print*,'------SCRATCH3_GRELL---------'
    !tmp      print*,'mgmxp mgmyp mgmzp ensdim ialloc'
    !tmp      print*,mgmxp,mgmyp,mgmzp,ensdim,ialloc
    !tmp      print*,'istart iend mix mjx mkx'
    !tmp      print*,istart,iend,mix,mjx,mkx
    !write(6,'(a1,78a1)') ' ',('_',i=1,78)


    ! 3d dependence (mgmxp,mgmzp,maxens2)
    allocate (dellat_ens (mgmxp,mgmzp,maxens2))
    allocate (dellaq_ens (mgmxp,mgmzp,maxens2))
    allocate (dellaqc_ens(mgmxp,mgmzp,maxens2))
    allocate (pwo_ens    (mgmxp,mgmzp,maxens2))

    ! 3d dependence (mgmxp,mgmyp,ensdim)
    allocate (xf      (mgmxp,mgmyp,ensdim))
    allocate (xf_ens  (mgmxp,mgmyp,ensdim))
    allocate (pr_ens  (mgmxp,mgmyp,ensdim))
    allocate (outt_ens(mgmxp,mgmyp,ensdim))

    ! 2d dependence   (mgmxp,mgmzp)
    allocate (HE      (mgmxp,mgmzp))
    allocate (HES     (mgmxp,mgmzp))
    allocate (QES     (mgmxp,mgmzp))
    allocate (Z       (mgmxp,mgmzp))
    allocate (TV      (mgmxp,mgmzp))
    allocate (DBY     (mgmxp,mgmzp))
    allocate (QC      (mgmxp,mgmzp))
    allocate (QRCD    (mgmxp,mgmzp))
    allocate (PWD     (mgmxp,mgmzp))
    allocate (PW      (mgmxp,mgmzp))
    allocate (HEO     (mgmxp,mgmzp))
    allocate (HESO    (mgmxp,mgmzp))
    allocate (QESO    (mgmxp,mgmzp))
    allocate (ZO      (mgmxp,mgmzp))
    allocate (TVO     (mgmxp,mgmzp))
    allocate (DBYO    (mgmxp,mgmzp))
    allocate (QCO     (mgmxp,mgmzp))
    allocate (QRCDO   (mgmxp,mgmzp))
    allocate (PWDO    (mgmxp,mgmzp))
    allocate (PWO     (mgmxp,mgmzp))
    allocate (XHE     (mgmxp,mgmzp))
    allocate (XHES    (mgmxp,mgmzp))
    allocate (XQES    (mgmxp,mgmzp))
    allocate (XZ      (mgmxp,mgmzp))
    allocate (XTV     (mgmxp,mgmzp))
    allocate (XT_Grell(mgmxp,mgmzp))
    allocate (XQ      (mgmxp,mgmzp))
    allocate (XDBY    (mgmxp,mgmzp))
    allocate (XQC     (mgmxp,mgmzp))
    allocate (XQRCD   (mgmxp,mgmzp))
    allocate (XPWD    (mgmxp,mgmzp))
    allocate (XPW     (mgmxp,mgmzp))
    allocate (hcd     (mgmxp,mgmzp))
    allocate (hcdo    (mgmxp,mgmzp))
    allocate (xhcd    (mgmxp,mgmzp))
    allocate (qcd     (mgmxp,mgmzp))
    allocate (qcdo    (mgmxp,mgmzp))
    allocate (xqcd    (mgmxp,mgmzp))
    allocate (dbyd    (mgmxp,mgmzp))
    allocate (dbydo   (mgmxp,mgmzp))
    allocate (hc      (mgmxp,mgmzp))
    allocate (hco     (mgmxp,mgmzp))
    allocate (xhc     (mgmxp,mgmzp))
    allocate (qrc     (mgmxp,mgmzp))
    allocate (qrco    (mgmxp,mgmzp))
    allocate (xqrc    (mgmxp,mgmzp))
    allocate (zu      (mgmxp,mgmzp))
    allocate (zuo     (mgmxp,mgmzp))
    allocate (xzu     (mgmxp,mgmzp))
    allocate (zd      (mgmxp,mgmzp))
    allocate (zdo     (mgmxp,mgmzp))
    allocate (xzd     (mgmxp,mgmzp))
    allocate (DELLAH  (mgmxp,mgmzp))
    allocate (DELLAQ  (mgmxp,mgmzp))
    allocate (DELLAT  (mgmxp,mgmzp))
    allocate (DELLAQC (mgmxp,mgmzp))
    allocate (qes_cup (mgmxp,mgmzp))
    allocate (q_cup   (mgmxp,mgmzp))
    allocate (he_cup  (mgmxp,mgmzp))
    allocate (hes_cup (mgmxp,mgmzp))
    allocate (z_cup   (mgmxp,mgmzp))
    allocate (p_cup   (mgmxp,mgmzp))
    allocate (gamma_cup(mgmxp,mgmzp))
    allocate (t_cup   (mgmxp,mgmzp))
    allocate (qeso_cup(mgmxp,mgmzp))
    allocate (qo_cup  (mgmxp,mgmzp))
    allocate (heo_cup (mgmxp,mgmzp))
    allocate (heso_cup(mgmxp,mgmzp))
    allocate (zo_cup  (mgmxp,mgmzp))
    allocate (po_cup  (mgmxp,mgmzp))
    allocate (gammao_cup(mgmxp,mgmzp))
    allocate (tn_cup  (mgmxp,mgmzp))
    allocate (xqes_cup(mgmxp,mgmzp))
    allocate (xq_cup  (mgmxp,mgmzp))
    allocate (xhe_cup (mgmxp,mgmzp))
    allocate (xhes_cup(mgmxp,mgmzp))
    allocate (xz_cup  (mgmxp,mgmzp))
    allocate (xt_cup  (mgmxp,mgmzp))
    allocate (xt      (mgmxp,mgmzp))
    allocate (cd      (mgmxp,mgmzp))
    allocate (cdd     (mgmxp,mgmzp))
    allocate (scr1    (mgmxp,mgmzp))

    ! 2d dependence   (mgmxp,maxens)
    allocate (xaa0_ens(mgmxp,maxens))

    ! 2d dependence (mgmxp,maxens2)
    allocate (edtc  (mgmxp,maxens2))

    ! 2d dependence(mgmxp,mgmyp)
    allocate (cwf  (mgmxp,mgmyp))
    allocate (pwf  (mgmxp,mgmyp))
    allocate (pwdf (mgmxp,mgmyp))
    allocate (eddt (mgmxp,mgmyp))
    allocate (predb(mgmxp,mgmyp))
    allocate (xmass(mgmxp,mgmyp))

    ! 1d dependence  (mgmxp)
    allocate (kzdown (mgmxp))
    allocate (KBMAX  (mgmxp))
    allocate (IERR   (mgmxp))
    allocate (K22    (mgmxp))
    allocate (KBCON  (mgmxp))
    allocate (KB     (mgmxp))
    allocate (JMIN   (mgmxp))
    allocate (KTOP   (mgmxp))
    allocate (kstabi (mgmxp))
    allocate (kstabm (mgmxp))
    allocate (K22x   (mgmxp))
    allocate (KBCONx (mgmxp))
    allocate (KBx    (mgmxp))
    allocate (KTOPx  (mgmxp))
    allocate (kzi    (mgmxp))	 
    allocate (EDT    (mgmxp))
    allocate (EDTO   (mgmxp))
    allocate (EDTX   (mgmxp))
    allocate (AA1    (mgmxp))
    allocate (AA0    (mgmxp))
    allocate (XAA0   (mgmxp))
    allocate (HKB    (mgmxp))
    allocate (HKBO   (mgmxp))
    allocate (aad    (mgmxp))
    allocate (XHKB   (mgmxp))
    allocate (QKB    (mgmxp))
    allocate (QKBO   (mgmxp))
    allocate (XMB    (mgmxp))
    allocate (tkemax (mgmxp))
    allocate (XPWAV  (mgmxp))
    allocate (XPWEV  (mgmxp))
    allocate (PWAV   (mgmxp))
    allocate (PWEV   (mgmxp))
    allocate (PWAVO  (mgmxp))
    allocate (PWEVO  (mgmxp))
    allocate (BU     (mgmxp))
    allocate (BUO    (mgmxp))
    allocate (xfac1  (mgmxp))
    allocate (cap_max(mgmxp),cap_max_increment(mgmxp))
    allocate (vshear (mgmxp))
    allocate (sdp    (mgmxp))
    allocate (vws    (mgmxp))


    ! 1d dependence (maxens)
    allocate (xk      (maxens))
    allocate (mbdt_ens(maxens))

    ! 1d dependence (maxens2)
    allocate (edt_ens(maxens2))

    ! 1d dependence (maxens3)
    allocate (xff_ens3(maxens3))

    return
  end subroutine alloc_scratch3_grell_sh

  subroutine dealloc_scratch3_grell_sh !(scratch3_grell)

    implicit none
    !TYPE (scratch3_grell_vars) :: scratch3_grell

    ! 3d dependence (mgmxp,mgmzp,maxens2)
    deallocate (dellat_ens)
    deallocate (dellaq_ens)
    deallocate (dellaqc_ens)
    deallocate (pwo_ens)
    ! 3d dependence (mgmxp,mgmyp,ensdim)
    deallocate (xf)
    deallocate (xf_ens)
    deallocate (pr_ens)
    deallocate (outt_ens)
    deallocate (xf)

    ! 2d dependence (mgmxp,mgmzp)
    deallocate (HE)
    deallocate (HES)
    deallocate (QES)
    deallocate (Z)
    deallocate (TV)
    deallocate (DBY)
    deallocate (QC)
    deallocate (QRCD)
    deallocate (PWD)
    deallocate (PW)
    deallocate (HEO)
    deallocate (HESO)
    deallocate (QESO)
    deallocate (ZO)
    deallocate (TVO)
    deallocate (DBYO)
    deallocate (QCO)
    deallocate (QRCDO)
    deallocate (PWDO)
    deallocate (PWO)
    deallocate (XHE)
    deallocate (XHES)
    deallocate (XQES)
    deallocate (XZ)
    deallocate (XTV)
    deallocate (XT_Grell)
    deallocate (XQ)
    deallocate (XDBY)
    deallocate (XQC)
    deallocate (XQRCD)
    deallocate (XPWD)
    deallocate (XPW)
    deallocate (hcd)
    deallocate (hcdo)
    deallocate (xhcd)
    deallocate (qcd)
    deallocate (qcdo)
    deallocate (xqcd)
    deallocate (dbyd)
    deallocate (dbydo)
    deallocate (hc)
    deallocate (hco)
    deallocate (xhc)
    deallocate (qrc)
    deallocate (qrco)
    deallocate (xqrc)
    deallocate (zu)
    deallocate (zuo)
    deallocate (xzu)
    deallocate (zd)
    deallocate (zdo)
    deallocate (xzd)
    deallocate (DELLAH)
    deallocate (DELLAQ)
    deallocate (DELLAT)
    deallocate (DELLAQC)
    deallocate (qes_cup)
    deallocate (q_cup)
    deallocate (he_cup)
    deallocate (hes_cup)
    deallocate (z_cup)
    deallocate (p_cup)
    deallocate (gamma_cup)
    deallocate (t_cup)
    deallocate (qeso_cup)
    deallocate (qo_cup)
    deallocate (heo_cup)
    deallocate (heso_cup)
    deallocate (zo_cup)
    deallocate (po_cup)
    deallocate (gammao_cup)
    deallocate (tn_cup)
    deallocate (xqes_cup)
    deallocate (xq_cup)
    deallocate (xhe_cup)
    deallocate (xhes_cup)
    deallocate (xz_cup)
    deallocate (xt_cup)
    deallocate (xt)
    deallocate (cd)
    deallocate (cdd)
    deallocate (scr1)
    ! 2d dependence (mgmxp,maxens)
    deallocate (xaa0_ens)
    ! 2d dependence (mgmxp,maxens2)
    deallocate (edtc)
    ! 2d dependence (mgmxp,mgmyp)
    deallocate (cwf)
    deallocate (pwf)
    deallocate (pwdf)
    deallocate (eddt)
    deallocate (predb)
    deallocate (xmass)

    ! 1d dependence (mgmxp)
    deallocate (kzdown)
    deallocate (KBMAX)
    deallocate (IERR)
    deallocate (K22)
    deallocate (KBCON)
    deallocate (KB)
    deallocate (JMIN)
    deallocate (KTOP)
    deallocate (kstabi)
    deallocate (kstabm)
    deallocate (K22x)
    deallocate (KBCONx)
    deallocate (KBx)
    deallocate (KTOPx)
    deallocate (kzi)
    deallocate (EDT)
    deallocate (EDTO)
    deallocate (EDTX)
    deallocate (AA1)
    deallocate (AA0)
    deallocate (XAA0)
    deallocate (HKB)
    deallocate (HKBO)
    deallocate (aad)
    deallocate (XHKB)
    deallocate (QKB)
    deallocate (QKBO)
    deallocate (XMB)
    deallocate (tkemax)
    deallocate (XPWAV)
    deallocate (XPWEV)
    deallocate (PWAV)
    deallocate (PWEV)
    deallocate (PWAVO)
    deallocate (PWEVO)
    deallocate (BU)
    deallocate (BUO)
    ! 1d dependence (maxens3)
    deallocate (xff_ens3)
    ! 1d dependence (maxens)
    deallocate (xk)
    ! 1d dependence 
    deallocate (xfac1)
    deallocate (cap_max, cap_max_increment)
    ! 1d dependence (maxens)
    deallocate (mbdt_ens)
    ! 1d dependence (maxens2)
    deallocate (edt_ens)
    ! 1d dependence 
    !srf variaveis para a rotine "cup_dd_edt"
    deallocate (vshear)
    deallocate (sdp)
    deallocate (vws)

    return
  end subroutine dealloc_scratch3_grell_sh

  subroutine zero_scratch3_grell_sh()

    dellat_ens=0.
    dellaq_ens=0.
    dellaqc_ens=0.
    pwo_ens=0.
    xf=0.
    xf_ens=0.
    pr_ens=0.
    outt_ens=0.
    xf=0.
    HE=0.
    HES=0.
    QES=0.
    Z=0.
    TV=0.
    DBY=0.
    QC=0.
    QRCD=0.
    PWD=0.
    PW=0.
    HEO=0.
    HESO=0.
    QESO=0.
    ZO=0.
    TVO=0.
    DBYO=0.
    QCO=0.
    QRCDO=0.
    PWDO=0.
    PWO=0.
    XHE=0.
    XHES=0.
    XQES=0.
    XZ=0.
    XTV=0.
    XT_Grell=0.
    XQ=0.
    XDBY=0.
    XQC=0.
    XQRCD=0.
    XPWD=0.
    XPW=0.
    hcd=0.
    hcdo=0.
    xhcd=0.
    qcd=0.
    qcdo=0.
    xqcd=0.
    dbyd=0.
    dbydo=0.
    hc=0.
    hco=0.
    xhc=0.
    qrc=0.
    qrco=0.
    xqrc=0.
    zu=0.
    zuo=0.
    xzu=0.
    zd=0.
    zdo=0.
    xzd=0.
    DELLAH=0.
    DELLAQ=0.
    DELLAT=0.
    DELLAQC=0.
    qes_cup=0.
    q_cup=0.
    he_cup=0.
    hes_cup=0.
    z_cup=0.
    p_cup=0.
    gamma_cup=0.
    t_cup=0.
    qeso_cup=0.
    qo_cup=0.
    heo_cup=0.
    heso_cup=0.
    zo_cup=0.
    po_cup=0.
    gammao_cup=0.
    tn_cup=0.
    xqes_cup=0.
    xq_cup=0.
    xhe_cup=0.
    xhes_cup=0.
    xz_cup=0.
    xt_cup=0.
    xt=0.
    cd=0.
    cdd=0.
    scr1=0.
    xaa0_ens=0.
    edtc=0.
    cwf=0.
    pwf=0.
    pwdf=0.
    eddt=0.
    predb=0.
    xmass=0.

    kzdown=0
    KBMAX=0
    IERR=0
    K22=0
    KBCON=0
    KB=0
    JMIN=0
    KTOP=0
    kstabi=0
    kstabm=0
    K22x=0
    KBCONx=0
    KBx=0
    KTOPx=0
    kzi=0

    EDT=0.
    EDTO=0.
    EDTX=0.
    AA1=0.
    AA0=0.
    XAA0=0.
    HKB=0.
    HKBO=0.
    aad=0.
    XHKB=0.
    QKB=0.
    QKBO=0.
    XMB=0.
    tkemax=0.
    XPWAV=0.
    XPWEV=0.
    PWAV=0.
    PWEV=0.
    PWAVO=0.
    PWEVO=0.
    BU=0.
    BUO=0.
    xff_ens3=0.
    xk=0.
    xfac1=0.
    cap_max=0.
    cap_max_increment=0.

    mbdt_ens=0.
    edt_ens=0.
    vshear=0.
    sdp=0.
    vws=0.

  end subroutine zero_scratch3_grell_sh

end module mem_scratch3_grell_sh
