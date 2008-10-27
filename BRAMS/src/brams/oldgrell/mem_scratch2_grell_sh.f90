! Module necessary to Grell Cumulus param.
! Scratch variables - module 2

module mem_scratch2_grell_sh

  !TYPE scratch2_grell_vars

  !                   adapted in july-15-2002 for 5.x version
  !
  ! 3d dependence (mgmxp,mgmyp,ensdim)
  real, allocatable, dimension(:,:,:) :: massfln

  ! 2d dependence (mgmxp,mgmzp)
  real, allocatable, dimension(:,:) ::     &
       T,                                  &
       Q,                                  &
       P,                                  &
       PO,                                 &
       TN,                                 &
       QO,                                 &
       OUTT,                               &
       OUTQ,                               &
       outqc,                              &
       US,                                 & ! Substitui US original
       VS,                                 & ! Substitui VS original
       omeg
  ! 2d dependence (mgmxp, mgmyp)
  integer, allocatable, dimension(:,:) :: KDT, &
       iact_gr,                                &
       iact_old_gr
  real, allocatable, dimension(:,:) :: xland,  &
       massflx,                                &
       tkeg,                                   &
       rcpg


  ! 1d dependence (mgmxp)
  real, allocatable, dimension(:) :: mconv, &
       umean,                               &
       vmean,                               &
       pmean,                               &
       direction
  real, allocatable, dimension(:) :: AA0, &
       PRET,                              &
       PSUR,                              &
       TER11,                             &
       z1
  integer, allocatable, dimension(:) :: KDET, &
                                        PBLIDX

  !END TYPE scratch2_grell_vars

contains

  subroutine alloc_scratch2_grell_sh !(scratch2_grell)

    use mem_grell_param, only : mgmxp,  & ! INTENT(IN)
         mgmyp,                         & ! INTENT(IN)
         mgmzp,                         & ! INTENT(IN)
         ensdim                           ! INTENT(IN)
    use node_mod, only : mynum, &   ! INTENT(IN)
         MXP,                   &   ! INTENT(IN)
         MYP,                   &   ! INTENT(IN)
         MZP,                   &   ! INTENT(IN)
         IA,                    &   ! INTENT(IN)
         IZ,                    &   ! INTENT(IN)
         JA,                    &   ! INTENT(IN)
         JZ,                    &   ! INTENT(IN)
         I0,                    &   ! INTENT(IN)
         J0                         ! INTENT(IN)

    implicit none
    !TYPE (scratch2_grell_vars) :: scratch2_grell

    integer :: i

    !write(6,'(a1,78a1)') ' ',('_',i=1,78)
    !print*,' '
    !print*,'---- In cuparth - alocando memoria--'
    !print*,'------------SCRATCH2_sh-------------'
    !print*,'---------Node configuration---------'
    !print*,'MYNUM_I0____J0____=',mynum,i0,j0
    !print*,'mgmzp_mgmxp_mgmyp_=',mgmzp,mgmxp,mgmyp
    !print*,'m1____m2___m3_____=', mzp, mxp, myp !,m1,m2,m3
    !print*,'ia__iz____ja__jz__=',ia, iz, ja, jz
    !print*,'ensdim_____ialloc_=',ensdim !,ialloc
    !print*,' '
    !write(6,'(a1,78a1)') ' ',('_',i=1,78)

    allocate (massfln(mgmxp, mgmyp, ensdim))

    allocate (mconv    (mgmxp))
    allocate (umean    (mgmxp))
    allocate (vmean    (mgmxp))
    allocate (pmean    (mgmxp))
    allocate (direction(mgmxp))
    allocate (AA0      (mgmxp))
    allocate (KDET     (mgmxp))
    allocate (z1       (mgmxp))
    allocate (PRET     (mgmxp))
    allocate (PSUR     (mgmxp))
    allocate (TER11    (mgmxp))
    allocate (PBLIDX   (mgmxp))

    allocate (T          (mgmxp, mgmzp))
    allocate (Q          (mgmxp, mgmzp))
    allocate (P          (mgmxp, mgmzp))
    allocate (PO         (mgmxp, mgmzp))
    allocate (TN         (mgmxp, mgmzp))
    allocate (QO         (mgmxp, mgmzp))
    allocate (OUTT       (mgmxp, mgmzp))
    allocate (OUTQ       (mgmxp, mgmzp))
    allocate (outqc      (mgmxp, mgmzp))
    allocate (US         (mgmxp, mgmzp))
    allocate (VS         (mgmxp, mgmzp))
    allocate (omeg	 (mgmxp, mgmzp))
    allocate (KDT        (mgmxp, mgmyp))
    allocate (tkeg       (mgmxp, mgmzp))
    allocate (rcpg       (mgmxp, mgmzp))
    allocate (xland      (mgmxp, mgmyp))
    allocate (massflx    (mgmxp, mgmyp))
    allocate (iact_gr    (mgmxp, mgmyp))
    allocate (iact_old_gr(mgmxp, mgmyp))

    return
  end subroutine alloc_scratch2_grell_sh

  subroutine dealloc_scratch2_grell_sh !(scratch2_grell)

    implicit none
    !TYPE (scratch2_grell_vars) :: scratch2_grell

    deallocate (massfln)

    deallocate (mconv)
    deallocate (umean)
    deallocate (vmean)
    deallocate (pmean)
    deallocate (direction)

    deallocate (T)
    deallocate (Q)
    deallocate (P)
    deallocate (PO)
    deallocate (TN)
    deallocate (QO)
    deallocate (OUTT)
    deallocate (OUTQ)
    deallocate (outqc)
    deallocate (PRET)
    deallocate (PSUR)
    deallocate (TER11)
    deallocate (PBLIDX)
    deallocate (US)
    deallocate (VS)
    deallocate (omeg)
    deallocate (AA0)

    deallocate (KDET)
    deallocate (KDT)

    deallocate (xland)
    deallocate (massflx)
    deallocate (iact_gr)
    deallocate (iact_old_gr)

    return
  end subroutine dealloc_scratch2_grell_sh

  subroutine Zero_scratch2_grell_sh()
    massfln=0.
    PBLIDX=0
    T=0.
    Q=0.
    P=0.
    PO=0.
    TN=0.
    QO=0.
    OUTT=0.
    OUTQ=0.
    outqc=0.
    US=0.
    VS=0.
    omeg=0.
    xland=0.
    massflx=0.
    tkeg=0.
    rcpg=0.
    mconv=0.
    umean=0.
    vmean=0.
    pmean=0.
    direction=0.
    AA0=0.
    PRET=0.
    PSUR=0.
    TER11=0.
    z1=0.

    KDET=0
    KDT=0
    iact_gr=0
    iact_old_gr=0

  end subroutine Zero_scratch2_grell_sh

end module mem_scratch2_grell_sh
