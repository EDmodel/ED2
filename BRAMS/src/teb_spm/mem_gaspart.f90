MODULE mem_gaspart

  TYPE gaspart_vars

     ! Variables to be dimensioned by (nzp,nxp,nyp)
     REAL, POINTER, DIMENSION(:,:,:) :: &
          PCO,PNO,PNO2,PPM25,PVOC,PSO2,PROO,PSO4,PAER, &
	  PO3,PRHCO,PHO2,PO3P,PO1D,PHO,GASR,PEOXID
     ! Variables to be dimensioned by (nzp,nxp)
     REAL, POINTER, DIMENSION(:,:) :: fusog
 
     REAL, POINTER, DIMENSION(:) ::  &
          PCOT,PNOT,PNO2T,PPM25T,PVOCT,PSO2T,PSO4T,PAERT, &
	  PO3T,PRHCOT,PHO2T,PO3PT,PO1DT,PHOT,PROOT
     
  END TYPE gaspart_vars

  TYPE (gaspart_vars), ALLOCATABLE, target :: gaspart_g(:), gaspartm_g(:)

CONTAINS
  SUBROUTINE alloc_gaspart(gaspart,n1,n2,n3,ng)
    USE mem_emiss, ONLY: ichemi
    IMPLICIT NONE
    TYPE (gaspart_vars) :: gaspart
    INTEGER, INTENT(in) :: n1,n2,n3,ng


       ALLOCATE (gaspart%pco(n1,n2,n3))
       ALLOCATE (gaspart%pno(n1,n2,n3))
       ALLOCATE (gaspart%pno2(n1,n2,n3))
       ALLOCATE (gaspart%ppm25(n1,n2,n3))
       ALLOCATE (gaspart%pso2(n1,n2,n3))
       ALLOCATE (gaspart%pvoc(n1,n2,n3))
       ALLOCATE (gaspart%gasr(n1,n2,n3))
       ALLOCATE (gaspart%pso4(n1,n2,n3))
       ALLOCATE (gaspart%paer(n1,n2,n3))
       ALLOCATE (gaspart%PEOXID(n1,n2,n3))
       ALLOCATE (gaspart%fusog(n2,n3))

if(ichemi==1)then

       ALLOCATE (gaspart%po3   (n1,n2,n3))
       ALLOCATE (gaspart%prhco (n1,n2,n3))
       ALLOCATE (gaspart%pho2  (n1,n2,n3))
       ALLOCATE (gaspart%po3p  (n1,n2,n3))
       ALLOCATE (gaspart%po1d  (n1,n2,n3))
       ALLOCATE (gaspart%pho   (n1,n2,n3))
       ALLOCATE (gaspart%proo  (n1,n2,n3))
endif

    RETURN
  END SUBROUTINE alloc_gaspart


  SUBROUTINE nullify_gaspart(gaspart)
    USE mem_emiss, ONLY: ichemi

    IMPLICIT NONE
    TYPE (gaspart_vars) :: gaspart


    IF (ASSOCIATED(gaspart%pco))    NULLIFY (gaspart%pco)
    IF (ASSOCIATED(gaspart%pno))    NULLIFY (gaspart%pno)
    IF (ASSOCIATED(gaspart%pno2))   NULLIFY (gaspart%pno2)
    IF (ASSOCIATED(gaspart%ppm25))  NULLIFY (gaspart%ppm25)
    IF (ASSOCIATED(gaspart%pvoc))   NULLIFY (gaspart%pvoc)
    IF (ASSOCIATED(gaspart%pso2))   NULLIFY (gaspart%pso2)
    IF (ASSOCIATED(gaspart%pso4))   NULLIFY (gaspart%pso4)
    IF (ASSOCIATED(gaspart%paer))   NULLIFY (gaspart%paer)
    IF (ASSOCIATED(gaspart%PEOXID))   NULLIFY (gaspart%PEOXID)
    IF (ASSOCIATED(gaspart%gasr))   NULLIFY (gaspart%gasr)
    IF (ASSOCIATED(gaspart%fusog))   NULLIFY (gaspart%fusog)
if(ichemi==1)then

    IF (ASSOCIATED(gaspart%po3  ))  NULLIFY (gaspart%po3  )
    IF (ASSOCIATED(gaspart%prhco))  NULLIFY (gaspart%prhco)
    IF (ASSOCIATED(gaspart%pho2 ))  NULLIFY (gaspart%pho2 )
    IF (ASSOCIATED(gaspart%po3p ))  NULLIFY (gaspart%po3p )
    IF (ASSOCIATED(gaspart%po1d ))  NULLIFY (gaspart%po1d )
    IF (ASSOCIATED(gaspart%pho  ))  NULLIFY (gaspart%pho  )
    IF (ASSOCIATED(gaspart%proo ))  NULLIFY (gaspart%proo )
endif
    RETURN
  END SUBROUTINE nullify_gaspart

  SUBROUTINE dealloc_gaspart(gaspart)
    USE mem_emiss, ONLY: ichemi

    IMPLICIT NONE
    TYPE (gaspart_vars) :: gaspart


    IF (ASSOCIATED(gaspart%pco))    DEALLOCATE (gaspart%pco)
    IF (ASSOCIATED(gaspart%pno))    DEALLOCATE (gaspart%pno)
    IF (ASSOCIATED(gaspart%pno2))   DEALLOCATE (gaspart%pno2)
    IF (ASSOCIATED(gaspart%ppm25))  DEALLOCATE (gaspart%ppm25)
    IF (ASSOCIATED(gaspart%pvoc))   DEALLOCATE (gaspart%pvoc)
    IF (ASSOCIATED(gaspart%pso2))   DEALLOCATE (gaspart%pso2)
    IF (ASSOCIATED(gaspart%pso4))   DEALLOCATE (gaspart%pso4)
    IF (ASSOCIATED(gaspart%paer))   DEALLOCATE (gaspart%paer)
    IF (ASSOCIATED(gaspart%PEOXID))   DEALLOCATE (gaspart%PEOXID)
    IF (ASSOCIATED(gaspart%gasr))   DEALLOCATE (gaspart%gasr)
    IF (ASSOCIATED(gaspart%fusog))   DEALLOCATE (gaspart%fusog)
if(ichemi==1)then
    IF (ASSOCIATED(gaspart%po3  ))  DEALLOCATE (gaspart%po3  )
    IF (ASSOCIATED(gaspart%prhco))  DEALLOCATE (gaspart%prhco)
    IF (ASSOCIATED(gaspart%pho2 ))  DEALLOCATE (gaspart%pho2 )
    IF (ASSOCIATED(gaspart%po3p ))  DEALLOCATE (gaspart%po3p )
    IF (ASSOCIATED(gaspart%po1d ))  DEALLOCATE (gaspart%po1d )
    IF (ASSOCIATED(gaspart%pho  ))  DEALLOCATE (gaspart%pho  )
    IF (ASSOCIATED(gaspart%proo ))  DEALLOCATE (gaspart%proo )
endif
    RETURN
  END SUBROUTINE dealloc_gaspart


  SUBROUTINE filltab_gaspart(gaspart,gaspartm,imean,n1,n2,n3,ng)

    USE var_tables
    USE mem_emiss, ONLY: ichemi
    IMPLICIT NONE
    TYPE (gaspart_vars) :: gaspart,gaspartm
    INTEGER, INTENT(in) :: imean,n1,n2,n3,ng
    INTEGER :: npts
    REAL, POINTER :: var,varm

    ! Fill pointers to arrays into variable tables

    npts=n2*n3

    IF (ASSOCIATED(gaspart%fusog))  &
         CALL vtables2 (gaspart%FUSOG,gaspartm%FUSOG&
         ,ng, npts, imean,  &
         'FUSOG:2:hist:anal:mpti:mpt3:mpt1')

    npts=n1*n2*n3

    IF (ASSOCIATED(gaspart%pco))  &
         CALL vtables2 (gaspart%PCO,gaspartm%PCO&
         ,ng, npts, imean,  &
         'PCO:3:hist:anal:mpti:mpt3:mpt1')
    IF (ASSOCIATED(gaspart%pno))  &
         CALL vtables2 (gaspart%PNO,gaspartm%PNO&
         ,ng, npts, imean,  &
         'PNO:3:hist:anal:mpti:mpt3:mpt1')
    IF (ASSOCIATED(gaspart%pno2))  &
         CALL vtables2 (gaspart%PNO2,gaspartm%PNO2&
         ,ng, npts, imean,  &
         'PNO2:3:hist:anal:mpti:mpt3:mpt1')
    IF (ASSOCIATED(gaspart%ppm25))  &
         CALL vtables2 (gaspart%PPM25,gaspartm%PPM25&
         ,ng, npts, imean,  &
         'PPM25:3:hist:anal:mpti:mpt3:mpt1')
    IF (ASSOCIATED(gaspart%pvoc))  &
         CALL vtables2 (gaspart%PVOC,gaspartm%PVOC&
         ,ng, npts, imean,  &
         'PVOC:3:hist:anal:mpti:mpt3:mpt1')
    IF (ASSOCIATED(gaspart%pso2))  &
         CALL vtables2 (gaspart%PSO2,gaspartm%PSO2&
         ,ng, npts, imean,  &
         'PSO2:3:hist:anal:mpti:mpt3:mpt1')
    IF (ASSOCIATED(gaspart%pso4))  &
         CALL vtables2 (gaspart%PSO4,gaspartm%PSO4&
         ,ng, npts, imean,  &
         'PSO4:3:hist:anal:mpti:mpt3:mpt1')
    IF (ASSOCIATED(gaspart%paer))  &
         CALL vtables2 (gaspart%PAER,gaspartm%PAER&
         ,ng, npts, imean,  &
         'PAER:3:hist:anal:mpti:mpt3:mpt1')
    IF (ASSOCIATED(gaspart%PEOXID))  &
         CALL vtables2 (gaspart%PEOXID,gaspartm%PEOXID&
         ,ng, npts, imean,  &
         'PEOXID:3:hist:anal:mpti:mpt3:mpt1')
    IF (ASSOCIATED(gaspart%gasr))  &
         CALL vtables2 (gaspart%GASR,gaspartm%GASR&
         ,ng, npts, imean,  &
         'GASR:3:mpti:mpt3:mpt1')

if(ichemi==1)then

    IF (ASSOCIATED(gaspart%po3))  &
         CALL vtables2 (gaspart%PO3,gaspartm%PO3&
         ,ng, npts, imean,  &
         'PO3:3:hist:anal:mpti:mpt3:mpt1')
    IF (ASSOCIATED(gaspart%prhco))  &
         CALL vtables2 (gaspart%PRHCO,gaspartm%PRHCO&
         ,ng, npts, imean,  &
         'PRHCO:3:hist:anal:mpti:mpt3:mpt1')
    IF (ASSOCIATED(gaspart%pho2))  &
         CALL vtables2 (gaspart%PHO2,gaspartm%PHO2&
         ,ng, npts, imean,  &
         'PHO2:3:hist:anal:mpti:mpt3:mpt1')
    IF (ASSOCIATED(gaspart%po3p))  &
         CALL vtables2 (gaspart%PO3P,gaspartm%PO3P&
         ,ng, npts, imean,  &
         'PO3P:3:hist:anal:mpti:mpt3:mpt1')
    IF (ASSOCIATED(gaspart%po1d))  &
         CALL vtables2 (gaspart%PO1D,gaspartm%PO1D&
         ,ng, npts, imean,  &
         'PO1D:3:hist:anal:mpti:mpt3:mpt1')
    IF (ASSOCIATED(gaspart%pho))  &
         CALL vtables2 (gaspart%PHO,gaspartm%PHO&
         ,ng, npts, imean,  &
         'PHO:3:hist:anal:mpti:mpt3:mpt1')
    IF (ASSOCIATED(gaspart%proo))  &
         CALL vtables2 (gaspart%PROO,gaspartm%PROO&
         ,ng, npts, imean,  &
         'PROO:3:hist:anal:mpti:mpt3:mpt1')
endif
    RETURN
  END SUBROUTINE filltab_gaspart

END MODULE mem_gaspart
