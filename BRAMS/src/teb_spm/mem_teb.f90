! Module for urban canopy parameterization (land classes 19 and 21).

MODULE mem_teb

  TYPE teb_vars

     ! Variables to be dimensioned by (3,nxp,nyp)
     REAL, POINTER, DIMENSION(:,:,:) :: &
          T_ROOF,T_ROAD,T_WALL

     ! Variables to be dimensioned by (nxp,nyp)
     
     REAL, POINTER, DIMENSION(:,:) ::              &
     T_CANYON,R_CANYON,TS_ROOF,TS_ROAD,TS_WALL,    &
     TI_ROAD,WS_ROOF,WS_ROAD,TI_BLD,               &
     LE_TRAFFIC,H_TRAFFIC,LE_INDUSTRY,H_INDUSTRY,  &
     T2M_TOWN,R2M_TOWN,fuso

  END TYPE teb_vars

  TYPE (teb_vars), ALLOCATABLE, target :: teb_g(:), tebm_g(:)

CONTAINS

  SUBROUTINE alloc_teb(teb,n1,n2,n3,ng)

    IMPLICIT NONE
    TYPE (teb_vars) :: teb
    INTEGER, INTENT(in) :: n1,n2,n3,ng


       ALLOCATE (teb%T_ROOF(n1,n2,n3),teb%T_ROAD(n1,n2,n3), &
                 teb%T_WALL(n1,n2,n3) )

       ALLOCATE (teb%T_CANYON(n2,n3),teb%R_CANYON(n2,n3), &
                 teb%TS_ROOF(n2,n3),teb%TS_ROAD(n2,n3),teb%TS_WALL(n2,n3), &
                 teb%TI_ROAD(n2,n3),teb%WS_ROOF(n2,n3),teb%WS_ROAD(n2,n3), &
                 teb%TI_BLD(n2,n3),teb%LE_TRAFFIC(n2,n3),teb%H_TRAFFIC(n2,n3), &
		 teb%LE_INDUSTRY(n2,n3),teb%H_INDUSTRY(n2,n3),              &
		 teb%T2M_TOWN(n2,n3),teb%R2M_TOWN(n2,n3),teb%fuso(n2,n3))


    RETURN
  END SUBROUTINE alloc_teb


  SUBROUTINE nullify_teb(teb)

    IMPLICIT NONE
    TYPE (teb_vars) :: teb


    IF (ASSOCIATED(teb%T_ROOF))  NULLIFY (teb%T_ROOF)
    IF (ASSOCIATED(teb%T_ROAD))  NULLIFY (teb%T_ROAD)
    IF (ASSOCIATED(teb%T_WALL))  NULLIFY (teb%T_WALL)

    IF (ASSOCIATED(teb%T_CANYON))  NULLIFY  (teb%T_CANYON)
    IF (ASSOCIATED(teb%R_CANYON))  NULLIFY  (teb%R_CANYON)
    IF (ASSOCIATED(teb%TS_ROOF))   NULLIFY  (teb%TS_ROOF)
    IF (ASSOCIATED(teb%TS_ROAD))   NULLIFY  (teb%TS_ROAD)
    IF (ASSOCIATED(teb%TS_WALL))   NULLIFY  (teb%TS_WALL)
    IF (ASSOCIATED(teb%TI_ROAD))   NULLIFY  (teb%TI_ROAD)
    IF (ASSOCIATED(teb%WS_ROOF))   NULLIFY  (teb%WS_ROOF)
    IF (ASSOCIATED(teb%WS_ROAD))   NULLIFY  (teb%WS_ROAD)
    IF (ASSOCIATED(teb%TI_BLD))    NULLIFY  (teb%TI_BLD)
    IF (ASSOCIATED(teb%LE_TRAFFIC)) NULLIFY (teb%LE_TRAFFIC)
    IF (ASSOCIATED(teb%H_TRAFFIC))  NULLIFY (teb%H_TRAFFIC)
    IF (ASSOCIATED(teb%LE_INDUSTRY)) NULLIFY (teb%LE_INDUSTRY)
    IF (ASSOCIATED(teb%H_INDUSTRY))  NULLIFY (teb%H_INDUSTRY)
    IF (ASSOCIATED(teb%T2M_TOWN))    NULLIFY (teb%T2M_TOWN)
    IF (ASSOCIATED(teb%R2M_TOWN))    NULLIFY (teb%R2M_TOWN)
    IF (ASSOCIATED(teb%fuso))    NULLIFY (teb%fuso)

    RETURN
  END SUBROUTINE nullify_teb

  SUBROUTINE dealloc_teb(teb)

    IMPLICIT NONE
    TYPE (teb_vars) :: teb

    IF (ASSOCIATED(teb%T_ROOF))  DEALLOCATE (teb%T_ROOF)
    IF (ASSOCIATED(teb%T_ROAD))  DEALLOCATE (teb%T_ROAD)
    IF (ASSOCIATED(teb%T_WALL))  DEALLOCATE (teb%T_WALL)

    IF (ASSOCIATED(teb%T_CANYON))  DEALLOCATE  (teb%T_CANYON)
    IF (ASSOCIATED(teb%R_CANYON))  DEALLOCATE  (teb%R_CANYON)
    IF (ASSOCIATED(teb%TS_ROOF))   DEALLOCATE  (teb%TS_ROOF)
    IF (ASSOCIATED(teb%TS_ROAD))   DEALLOCATE  (teb%TS_ROAD)
    IF (ASSOCIATED(teb%TS_WALL))   DEALLOCATE  (teb%TS_WALL)
    IF (ASSOCIATED(teb%TI_ROAD))   DEALLOCATE  (teb%TI_ROAD)
    IF (ASSOCIATED(teb%WS_ROOF))   DEALLOCATE  (teb%WS_ROOF)
    IF (ASSOCIATED(teb%WS_ROAD))   DEALLOCATE  (teb%WS_ROAD)
    IF (ASSOCIATED(teb%TI_BLD))    DEALLOCATE  (teb%TI_BLD)
    IF (ASSOCIATED(teb%LE_TRAFFIC)) DEALLOCATE (teb%LE_TRAFFIC)
    IF (ASSOCIATED(teb%H_TRAFFIC))  DEALLOCATE (teb%H_TRAFFIC)
    IF (ASSOCIATED(teb%LE_INDUSTRY)) DEALLOCATE (teb%LE_INDUSTRY)
    IF (ASSOCIATED(teb%H_INDUSTRY))  DEALLOCATE (teb%H_INDUSTRY)

    IF (ASSOCIATED(teb%T2M_TOWN))    DEALLOCATE (teb%T2M_TOWN)
    IF (ASSOCIATED(teb%R2M_TOWN))    DEALLOCATE (teb%R2M_TOWN)
    IF (ASSOCIATED(teb%fuso))    DEALLOCATE (teb%fuso)

    RETURN
  END SUBROUTINE dealloc_teb


  SUBROUTINE filltab_teb(teb,tebm,imean,n1,n2,n3,ng)

    USE var_tables

    IMPLICIT NONE
    TYPE (teb_vars) :: teb,tebm
    INTEGER, INTENT(in) :: imean,n1,n2,n3,ng
    INTEGER :: npts
    REAL, POINTER :: var,varm

    ! Fill pointers to arrays into variable tables

    npts=n1*n2*n3

    IF (ASSOCIATED(teb%T_ROOF))  &
         CALL vtables2 (teb%T_ROOF,tebm%T_ROOF&
         ,ng, npts, imean,  &
         'T_ROOF :3:hist:anal:lite:mpti:mpt3:mpt1')
    IF (ASSOCIATED(teb%T_ROAD))  &
         CALL vtables2 (teb%T_ROAD,tebm%T_ROAD&
         ,ng, npts, imean,  &
         'T_ROAD :3:hist:anal:lite:mpti:mpt3:mpt1')
    IF (ASSOCIATED(teb%T_WALL))  &
         CALL vtables2 (teb%T_WALL,tebm%T_WALL&
         ,ng, npts, imean,  &
         'T_WALL :3:hist:anal:lite:mpti:mpt3:mpt1')

    npts=n2*n3

       IF (ASSOCIATED(teb%T_CANYON))  &
         CALL vtables2 (teb%T_CANYON,tebm%T_CANYON&
         ,ng, npts, imean,  &
         'T_CANYON :2:hist:anal:lite:mpti:mpt3:mpt1')

       IF (ASSOCIATED(teb%R_CANYON))  &
         CALL vtables2 (teb%R_CANYON,tebm%R_CANYON&
         ,ng, npts, imean,  &
         'R_CANYON :2:hist:anal:lite:mpti:mpt3:mpt1')

       IF (ASSOCIATED(teb%TS_ROOF))  &
         CALL vtables2 (teb%TS_ROOF,tebm%TS_ROOF&
         ,ng, npts, imean,  &
         'TS_ROOF :2:hist:anal:lite:mpti:mpt3:mpt1')

       IF (ASSOCIATED(teb%TS_ROAD))  &
         CALL vtables2 (teb%TS_ROAD,tebm%TS_ROAD&
         ,ng, npts, imean,  &
         'TS_ROAD :2:hist:anal:lite:mpti:mpt3:mpt1')

       IF (ASSOCIATED(teb%TS_WALL))  &
         CALL vtables2 (teb%TS_WALL,tebm%TS_WALL&
         ,ng, npts, imean,  &
         'TS_WALL :2:hist:anal:lite:mpti:mpt3:mpt1')

       IF (ASSOCIATED(teb%TI_ROAD))  &
         CALL vtables2 (teb%TI_ROAD,tebm%TI_ROAD&
         ,ng, npts, imean,  &
         'TI_ROAD :2:hist:anal:lite:mpti:mpt3:mpt1')

       IF (ASSOCIATED(teb%WS_ROOF))  &
         CALL vtables2 (teb%WS_ROOF,tebm%WS_ROOF&
         ,ng, npts, imean,  &
         'WS_ROOF :2:hist:anal:lite:mpti:mpt3:mpt1')

       IF (ASSOCIATED(teb%WS_ROAD))  &
         CALL vtables2 (teb%WS_ROAD,tebm%WS_ROAD&
         ,ng, npts, imean,  &
         'WS_ROAD :2:hist:anal:lite:mpti:mpt3:mpt1')

       IF (ASSOCIATED(teb%TI_BLD))  &
         CALL vtables2 (teb%TI_BLD,tebm%TI_BLD&
         ,ng, npts, imean,  &
         'TI_BLD :2:hist:anal:lite:mpti:mpt3:mpt1')
       IF (ASSOCIATED(teb%LE_TRAFFIC))  &
         CALL vtables2 (teb%LE_TRAFFIC,tebm%LE_TRAFFIC&
         ,ng, npts, imean,  &
         'LE_TRAFFIC :2:hist:anal:lite:mpti:mpt3:mpt1')

       IF (ASSOCIATED(teb%H_TRAFFIC))  &
         CALL vtables2 (teb%H_TRAFFIC,tebm%H_TRAFFIC&
         ,ng, npts, imean,  &
         'H_TRAFFIC :2:hist:anal:lite:mpti:mpt3:mpt1')

       IF (ASSOCIATED(teb%LE_INDUSTRY))  &
         CALL vtables2 (teb%LE_INDUSTRY,tebm%LE_INDUSTRY&
         ,ng, npts, imean,  &
         'LE_INDUSTRY :2:hist:anal:lite:mpti:mpt3:mpt1')

       IF (ASSOCIATED(teb%H_INDUSTRY))  &
         CALL vtables2 (teb%H_INDUSTRY,tebm%H_INDUSTRY&
         ,ng, npts, imean,  &
         'H_INDUSTRY :2:hist:anal:lite:mpti:mpt3:mpt1')

       IF (ASSOCIATED(teb%T2M_TOWN))  &
         CALL vtables2 (teb%T2M_TOWN,tebm%T2M_TOWN&
         ,ng, npts, imean,  &
         'T2M_TOWN :2:hist:anal:lite:mpti:mpt3:mpt1')
	 
       IF (ASSOCIATED(teb%R2M_TOWN))  &
         CALL vtables2 (teb%R2M_TOWN,tebm%R2M_TOWN&
         ,ng, npts, imean,  &
         'R2M_TOWN :2:hist:anal:lite:mpti:mpt3:mpt1')
	 
       IF (ASSOCIATED(teb%fuso))  &
         CALL vtables2 (teb%fuso,tebm%fuso&
         ,ng, npts, imean,  &
         'FUSO :2:hist:anal:lite:mpti:mpt3:mpt1')
	 
    RETURN
  END SUBROUTINE filltab_teb

END MODULE mem_teb
