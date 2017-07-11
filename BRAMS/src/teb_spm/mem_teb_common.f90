MODULE mem_teb_common

  TYPE teb_common

     ! Variables to be dimensioned by (nxp,nyp)
     
     REAL, POINTER, DIMENSION(:,:) ::              &

     EMIS_TOWN,ALB_TOWN,TS_TOWN

  END TYPE teb_common

  TYPE (teb_common), ALLOCATABLE, target :: tebc_g(:), tebcm_g(:)

CONTAINS
  SUBROUTINE alloc_tebc(tebc,n1,n2,n3,ng)

    IMPLICIT NONE
    TYPE (teb_common) :: tebc
    INTEGER, INTENT(in) :: n1,n2,n3,ng

       ALLOCATE (tebc%EMIS_TOWN(n2,n3),tebc%ALB_TOWN(n2,n3), &
                 tebc%TS_TOWN(n2,n3))

    RETURN
  END SUBROUTINE alloc_tebc


  SUBROUTINE nullify_tebc(tebc)

    IMPLICIT NONE
    TYPE (teb_common) :: tebc

    IF (ASSOCIATED(tebc%EMIS_TOWN))  NULLIFY (tebc%EMIS_TOWN)
    IF (ASSOCIATED(tebc%ALB_TOWN))   NULLIFY (tebc%ALB_TOWN)
    IF (ASSOCIATED(tebc%TS_TOWN))    NULLIFY (tebc%TS_TOWN)

    RETURN
  END SUBROUTINE nullify_tebc

  SUBROUTINE dealloc_tebc(tebc)

    IMPLICIT NONE

    TYPE (teb_common) :: tebc

    IF (ASSOCIATED(tebc%EMIS_TOWN))  DEALLOCATE (tebc%EMIS_TOWN)
    IF (ASSOCIATED(tebc%ALB_TOWN))   DEALLOCATE (tebc%ALB_TOWN)
    IF (ASSOCIATED(tebc%TS_TOWN))    DEALLOCATE (tebc%TS_TOWN)

    RETURN
  END SUBROUTINE dealloc_tebc


  SUBROUTINE filltab_tebc(tebc,tebcm,imean,n1,n2,n3,ng)

    USE var_tables

    IMPLICIT NONE
    TYPE (teb_common) :: tebc,tebcm
    INTEGER, INTENT(in) :: imean,n1,n2,n3,ng
    INTEGER :: npts
    REAL, POINTER :: var,varm

    ! Fill pointers to arrays into variable tables


    npts=n2*n3

       IF (ASSOCIATED(tebc%EMIS_TOWN))  &
         CALL vtables2 (tebc%EMIS_TOWN,tebcm%EMIS_TOWN&
         ,ng, npts, imean,  &
         'EMIS_TOWN :2:hist:anal:lite:mpti:mpt3:mpt1')
       IF (ASSOCIATED(tebc%ALB_TOWN))  &
         CALL vtables2 (tebc%ALB_TOWN,tebcm%ALB_TOWN&
         ,ng, npts, imean,  &
         'ALB_TOWN :2:hist:anal:lite:mpti:mpt3:mpt1')
       IF (ASSOCIATED(tebc%TS_TOWN))  &
         CALL vtables2 (tebc%TS_TOWN,tebcm%TS_TOWN&
         ,ng, npts, imean,  &
         'TS_TOWN :2:hist:anal:lite:mpti:mpt3:mpt1')
    RETURN
  END SUBROUTINE filltab_tebc

END MODULE mem_teb_common
