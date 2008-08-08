MODULE mem_sib_co2

  TYPE sib_co2_vars

     ! Variables to be dimensioned by (nxp, nyp)
     !REAL, POINTER, DIMENSION(:, :, :) :: SRC_CO2 ! dimens.(nxp,nyp,N_CO2)
     REAL, POINTER, DIMENSION(:, :) :: SRC_CO2

  END TYPE sib_co2_vars

  TYPE (sib_co2_vars), ALLOCATABLE :: sib_g(:, :), sibm_g(:, :)

CONTAINS

  SUBROUTINE alloc_sib_co2(sib, n2, n3)

    IMPLICIT NONE
    TYPE (sib_co2_vars) :: sib
    INTEGER, INTENT(in) :: n2, n3

    ! Allocate arrays based on options (if necessary)

    ALLOCATE (sib%SRC_CO2(n2, n3))  !(sib%SRC_CO2(n2,n3,N_CO2))

    RETURN
  END SUBROUTINE alloc_sib_co2

  SUBROUTINE nullify_sib_co2(sib)

    IMPLICIT NONE
    TYPE (sib_co2_vars) :: sib

    IF (ASSOCIATED(sib%SRC_CO2))  NULLIFY (sib%SRC_CO2)

    RETURN
  END SUBROUTINE nullify_sib_co2

  SUBROUTINE dealloc_sib_co2(sib)

    IMPLICIT NONE
    TYPE (sib_co2_vars) :: sib

    IF (ASSOCIATED(sib%SRC_CO2))  DEALLOCATE (sib%SRC_CO2)

    RETURN
  END SUBROUTINE dealloc_sib_co2

  SUBROUTINE filltab_sib_co2(sib, sibm, imean, n2, n3, ng)

    USE var_tables

    IMPLICIT NONE
    TYPE (sib_co2_vars) :: sib, sibm
    INTEGER, INTENT(in) :: imean, n2, n3, ng
    INTEGER :: npts
    ! REAL, POINTER :: var,varm

    ! Fill pointers to arrays into variable tables

    npts=n2*n3   !*N_CO2

    IF (ASSOCIATED(sib%SRC_CO2))  &
         CALL vtables2 (sib%SRC_CO2(1,1),sibm%SRC_CO2(1,1) & !(sib%SRC_CO2(1,1,1),sibm%SRC_CO2(1,1,1) &
         ,ng, npts, imean,  &
         !'SRC_CO2 :3:hist:anal:mpti:mpt3')
         'SRC_CO2 :2:hist:anal:mpti:mpt3')

    RETURN
  END SUBROUTINE filltab_sib_co2

  SUBROUTINE zero_sib_co2(sib, n2, n3)

    IMPLICIT NONE
    TYPE (sib_co2_vars) :: sib
    INTEGER, INTENT(in) :: n2, n3

    INTEGER :: i, j

    DO j=1,n3
       DO i=1,n2
          sib%SRC_CO2(i, j) = 0.
       ENDDO
    ENDDO

  END SUBROUTINE zero_sib_co2

END MODULE mem_sib_co2
