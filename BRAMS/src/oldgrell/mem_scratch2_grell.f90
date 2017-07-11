! Module necessary to Grell Cumulus param.
! Scratch variables - module 2

MODULE mem_scratch2_grell

  !TYPE scratch2_grell_vars

     !                   adapted in july-15-2002 for 5.x version
     !
     ! 3d dependence (mgmxp,mgmyp,ensdim)
     REAL, ALLOCATABLE, DIMENSION(:,:,:) :: massfln

     ! 2d dependence (mgmxp,mgmzp)
     REAL, ALLOCATABLE, DIMENSION(:,:) :: T,      &
          Q,                                  &
          P,                                  &
          PO,                                 &
          TN,                                 &
          QO,                                 &
          OUTT,                               &
          OUTQ,                               &
          outqc,                              &
          US_Grell,                           & ! Substitui US original
          VS_Grell,                           & ! Substitui VS original
          omeg
     ! 2d dependence (mgmxp, mgmyp)
     INTEGER, ALLOCATABLE, DIMENSION(:,:) :: KDT, &
          iact_gr,                            &
          iact_old_gr
     REAL, ALLOCATABLE, DIMENSION(:,:) :: xland,  &
          tkeg,                                   &
          rcpg,                                   &
          massflx

     ! 1d dependence (mgmxp)
     REAL, ALLOCATABLE, DIMENSION(:) :: mconv,    &
          umean,                              &
          vmean,                              &
          pmean,                              &
          direction
     REAL, ALLOCATABLE, DIMENSION(:) ::       &
           AA0,                               &
          PRET,                               &
          PSUR,                               &
          TER11
     INTEGER, ALLOCATABLE, DIMENSION(:) :: KDET, &
                                           PBLIDX

  !END TYPE scratch2_grell_vars

CONTAINS

  SUBROUTINE alloc_scratch2_grell !(scratch2_grell)

    USE mem_grell_param, ONLY : mgmxp,  & ! INTENT(IN)
         mgmyp,                         & ! INTENT(IN)
         mgmzp,                         & ! INTENT(IN)
         ensdim                           ! INTENT(IN)
    USE node_mod, ONLY : mynum, &   ! INTENT(IN)
         MXP,                   &   ! INTENT(IN)
         MYP,                   &   ! INTENT(IN)
         MZP,                   &   ! INTENT(IN)
         IA,                    &   ! INTENT(IN)
         IZ,                    &   ! INTENT(IN)
         JA,                    &   ! INTENT(IN)
         JZ,                    &   ! INTENT(IN)
         I0,                    &   ! INTENT(IN)
         J0                         ! INTENT(IN)

    IMPLICIT NONE
    !TYPE (scratch2_grell_vars) :: scratch2_grell

    INTEGER :: i

    ALLOCATE (massfln(mgmxp, mgmyp, ensdim))

    ALLOCATE (mconv    (mgmxp))
    ALLOCATE (umean    (mgmxp))
    ALLOCATE (vmean    (mgmxp))
    ALLOCATE (pmean    (mgmxp))
    ALLOCATE (direction(mgmxp))
    ALLOCATE (PRET     (mgmxp))
    ALLOCATE (PSUR     (mgmxp))
    ALLOCATE (TER11    (mgmxp))
    ALLOCATE (AA0      (mgmxp))
    ALLOCATE (KDET     (mgmxp))
    ALLOCATE (PBLIDX   (mgmxp))

    ALLOCATE (T          (mgmxp, mgmzp))
    ALLOCATE (Q          (mgmxp, mgmzp))
    ALLOCATE (P          (mgmxp, mgmzp))
    ALLOCATE (PO         (mgmxp, mgmzp))
    ALLOCATE (TN         (mgmxp, mgmzp))
    ALLOCATE (QO         (mgmxp, mgmzp))
    ALLOCATE (OUTT       (mgmxp, mgmzp))
    ALLOCATE (OUTQ       (mgmxp, mgmzp))
    ALLOCATE (outqc      (mgmxp, mgmzp))
    ALLOCATE (US_Grell   (mgmxp, mgmzp))
    ALLOCATE (VS_Grell   (mgmxp, mgmzp))
    ALLOCATE (omeg       (mgmxp, mgmzp))
    ALLOCATE (KDT        (mgmxp, mgmyp))
    ALLOCATE (xland      (mgmxp, mgmyp))
    ALLOCATE (massflx    (mgmxp, mgmyp))
    ALLOCATE (iact_gr    (mgmxp, mgmyp))
    ALLOCATE (iact_old_gr(mgmxp, mgmyp))

    ALLOCATE (tkeg       (mgmxp, mgmzp))
    ALLOCATE (rcpg       (mgmxp, mgmzp))

    RETURN
  END SUBROUTINE alloc_scratch2_grell

  SUBROUTINE dealloc_scratch2_grell !(scratch2_grell)

    IMPLICIT NONE
    !TYPE (scratch2_grell_vars) :: scratch2_grell

    DEALLOCATE (massfln)

    DEALLOCATE (mconv)
    DEALLOCATE (umean)
    DEALLOCATE (vmean)
    DEALLOCATE (pmean)
    DEALLOCATE (direction)

    DEALLOCATE (T)
    DEALLOCATE (Q)
    DEALLOCATE (P)
    DEALLOCATE (PO)
    DEALLOCATE (TN)
    DEALLOCATE (QO)
    DEALLOCATE (OUTT)
    DEALLOCATE (OUTQ)
    DEALLOCATE (outqc)
    DEALLOCATE (PRET)
    DEALLOCATE (PSUR)
    DEALLOCATE (TER11)
    DEALLOCATE (US_Grell)
    DEALLOCATE (VS_Grell)
    DEALLOCATE (omeg)
    DEALLOCATE (AA0)

    DEALLOCATE (PBLIDX)
    DEALLOCATE (KDET)
    DEALLOCATE (KDT)

    DEALLOCATE (xland)
    DEALLOCATE (massflx)
    DEALLOCATE (iact_gr)
    DEALLOCATE (iact_old_gr)

    DEALLOCATE (tkeg)
    DEALLOCATE (rcpg)

    RETURN
  END SUBROUTINE dealloc_scratch2_grell

  SUBROUTINE zero_scratch2_grell()
       massfln=0.
       T=0.
       Q=0.
       P=0.
       PO=0.
       TN=0.
       QO=0.
       OUTT=0.
       OUTQ=0.
       outqc=0.
       US_Grell=0.
       VS_Grell=0.
       omeg=0.
       xland=0.
       massflx=0.
       mconv=0.
       umean=0.
       vmean=0.
       pmean=0.
       direction=0.
       AA0=0.
       PRET=0.
       PSUR=0.
       TER11=0.
       KDT=0
       iact_gr=0
       iact_old_gr=0
       KDET=0
       PBLIDX=0
       tkeg=0.
       rcpg=0.
  END SUBROUTINE zero_scratch2_grell

END MODULE mem_scratch2_grell
