!!$SUBROUTINE htint (nzz1, vctra, eleva, nzz2, vctrb, elevb)
!!$  IMPLICIT NONE
!!$  INTEGER, INTENT(IN ) :: nzz1
!!$  INTEGER, INTENT(IN ) :: nzz2
!!$  REAL,    INTENT(IN ) :: vctra(nzz1)
!!$  REAL,    INTENT(OUT) :: vctrb(nzz2)
!!$  REAL,    INTENT(IN ) :: eleva(nzz1)
!!$  REAL,    INTENT(IN ) :: elevb(nzz2)
!!$
!!$  INTEGER :: ind1(nzz2)
!!$  INTEGER :: ind2(nzz2)
!!$  REAL    :: weight(nzz2)
!!$
!!$  ! interpolation weight and index
!!$
!!$  CALL htint_index(nzz1, eleva, nzz2, elevb, ind1, ind2, weight)
!!$  !CALL htint_index(nzz1, vctra, eleva, nzz2, vctrb, elevb, ind1, ind2, weight)
!!$
!!$  ! interpolate
!!$
!!$  CALL htint_inter(nzz1, vctra, nzz2, vctrb, ind1, ind2, weight)
!!$END SUBROUTINE htint








SUBROUTINE htint_inter (nzz1, vctra, nzz2, vctrb, ind1, ind2, weight)
  IMPLICIT NONE
  INTEGER, INTENT(IN ) :: nzz1
  INTEGER, INTENT(IN ) :: nzz2
  REAL,    INTENT(IN ) :: vctra(nzz1)
  REAL,    INTENT(OUT) :: vctrb(nzz2)
  INTEGER, INTENT(IN ) :: ind1(nzz2)
  INTEGER, INTENT(IN ) :: ind2(nzz2)
  REAL,    INTENT(IN ) :: weight(nzz2)

  INTEGER :: k

!!$  REAL    :: temp

  DO k=1,nzz2

!!$     temp = weight(k)
!!$     vctrb(k)  = vctra(ind1(k))*(1 - temp) + vctra(ind2(k))*temp

     vctrb(k)  = vctra(ind1(k))*(1 - weight(k)) + vctra(ind2(k))*weight(k)

!!$     vctrb(k)  = vctra(ind1(k))+(vctra(ind2(k))-vctra(ind1(k)))*weight(k)

  END DO
END SUBROUTINE htint_inter





SUBROUTINE htint_index (nzz1, eleva, nzz2, elevb, ind1, ind2, weight)
!!$SUBROUTINE htint_index (nzz1, vctra, eleva, nzz2, vctrb, elevb, ind1, ind2, weight)
  IMPLICIT NONE
  INTEGER, INTENT(IN ) :: nzz1
  INTEGER, INTENT(IN ) :: nzz2

!!$  REAL,    INTENT(IN ) :: vctra(nzz1)
!!$  REAL,    INTENT(OUT) :: vctrb(nzz2)

  REAL,    INTENT(IN ) :: eleva(nzz1)
  REAL,    INTENT(IN ) :: elevb(nzz2)
  INTEGER, INTENT(OUT) :: ind1(nzz2)
  INTEGER, INTENT(OUT) :: ind2(nzz2)
  REAL,    INTENT(OUT) :: weight(nzz2)


!!$  ! FOR DEBUG
!!$  REAL :: vctrc(nzz2), dif
!!$  INTEGER :: kcheck

  INTEGER :: l
  INTEGER :: k
  INTEGER :: kk

!!$  kcheck = 0

  l=1

  DO k=1,nzz2
     DO   !l = l, nzz1-1
        IF (  elevb(k) <  eleva(1) .OR. &
             (elevb(k) >= eleva(l) .AND. elevb(k) <= eleva(l+1))) THEN
           weight(k) = (elevb(k)-eleva(l))/(eleva(l+1)-eleva(l))
           ind1(k)   = l
           ind2(k)   = l+1

!!$           vctrb(k)  = vctra(ind1(k))+(vctra(ind2(k))-vctra(ind1(k)))*weight(k)
!!$
!!$           kcheck = kcheck + 1

           EXIT
        ELSE IF ( elevb(k) >  eleva(nzz1))  THEN
           weight(k) = (elevb(k)-eleva(nzz1))/(eleva(nzz1-1)-eleva(nzz1))
           ind1(k)   = nzz1
           ind2(k)   = nzz1-1

!!$           vctrb(k)  = vctra(ind1(k))+(vctra(ind2(k))-vctra(ind1(k)))*weight(k)
!!$
!!$           kcheck = kcheck + 1

           EXIT
        END IF

        l=l+1
        IF(l == nzz1) THEN
           PRINT *,'htint:nzz1',nzz1
           DO kk=1,l
              PRINT*,'kk,eleva(kk),elevb(kk)',eleva(kk),elevb(kk)
           END DO
           STOP 'htint'
        END IF
     END DO

!!$     vctrb(k)  = vctra(ind1(k))+(vctra(ind2(k))-vctra(ind1(k)))*weight(k)

  END DO

!!$  ! DEBUGGING
!!$
!!$  ! Verifica se todos os pontos foram calculados
!!$  IF (kcheck /= nzz2) PRINT *, "htint=>ERRO>nzz2,kcheck=", nzz2, kcheck
!!$
!!$  ! Verifica diferencas no calculo
!!$  DO k=1,nzz2
!!$     IF ((ind1(k) > nzz1).OR.(ind2(k) > nzz1)) THEN
!!$          PRINT *, "htint=>ERRO>k,nzz1,ind1,ind2", k,nzz1,ind1,ind2
!!$          STOP
!!$     END IF
!!$     vctrc(k)  = vctra(ind1(k))+(vctra(ind2(k))-vctra(ind1(k)))*weight(k)
!!$     dif = vctrb(k) - vctrc(k)
!!$     IF (abs(dif)>1.e-4)  PRINT *, "htint=>k,dif, dif_rel=", &
!!$          k,dif,abs(dif/vctra(k))
!!$  ENDDO

END SUBROUTINE htint_index




SUBROUTINE htint (nzz1, vctra, eleva, nzz2, vctrb, elevb)
  IMPLICIT NONE
  INTEGER, INTENT(IN ) :: nzz1
  INTEGER, INTENT(IN ) :: nzz2
  REAL,    INTENT(IN ) :: vctra(nzz1)
  REAL,    INTENT(OUT) :: vctrb(nzz2)
  REAL,    INTENT(IN ) :: eleva(nzz1)
  REAL,    INTENT(IN ) :: elevb(nzz2)

  INTEGER :: l
  INTEGER :: k
  INTEGER :: kk
  REAL    :: wt

  l=1

  DO k=1,nzz2
     DO
        IF ( (elevb(k) <  eleva(1)) .OR. &
             ((elevb(k) >= eleva(l)) .AND. (elevb(k) <= eleva(l+1))) ) THEN
           wt       = (elevb(k)-eleva(l))/(eleva(l+1)-eleva(l))
           vctrb(k) = vctra(l)+(vctra(l+1)-vctra(l))*wt
           EXIT
        ELSE IF ( elevb(k) >  eleva(nzz1))  THEN
           wt       = (elevb(k)-eleva(nzz1))/(eleva(nzz1-1)-eleva(nzz1))
           vctrb(k) = vctra(nzz1)+(vctra(nzz1-1)-vctra(nzz1))*wt
           EXIT
        END IF

        l=l+1
        IF(l == nzz1) THEN
           PRINT *,'htint:nzz1',nzz1
           DO kk=1,l
              PRINT*,'kk,eleva(kk),elevb(kk)',eleva(kk),elevb(kk)
           END DO
           STOP 'htint'
        END IF
     END DO
  END DO
END SUBROUTINE htint
