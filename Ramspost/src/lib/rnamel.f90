!############################# Change Log ##################################
! 1.0.0.2
!
! 001012 MJB TOKNZE ##
!            Declared passed args with length * instead of 1 to pass array
!            bounds compiler option for REVU. ##
! 001002 MJB CH2INT CH2REAL FILLPG ##
!            Replaced index calls with f90 intrinsic len_trim. ##
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!  Mission Research Corporation / *ASTeR Division
!###########################################################################

!------------------------------------------------------------------
!
!      Namelist reading routines
!
!------------------------------------------------------------------

SUBROUTINE VARSETF (VARN,VAR,IS1,MAXSUB,FNUM,FMIN,FMAX)

CHARACTER*(*) VARN

! Assign real value read in (FNUM) to the corresponding real variable
! in the NAMELIST and do a bounds check

IF(IS1.LE.MAXSUB) THEN
  VAR=FNUM
ELSE
  PRINT 9,VARN(1:lastchar(VARN)),IS1,MAXSUB,FNUM
  9 FORMAT(' -- ERROR --   Input variable - ',A,' - attempting'  &
        ,' to read extra values, values ignored.',/,  &
         ' Subscript, max dimension, value ',2I6,F20.6)
ENDIF

IF(FNUM.LT.FMIN.OR.FNUM.GT.FMAX) THEN
  PRINT 10,VARN(1:lastchar(VARN)),FNUM,FMIN,FMAX
  10 FORMAT(' -- ERROR --   Input variable - ',A,' - set to ',F18.5  &
      ,/,'                 allowable range ',F15.5,' to ',F15.5)
ENDIF

RETURN
END

!***************************************************************************

SUBROUTINE VARSETF2 (VARN,VAR,IS1,MAX1,IS2,MAX2,FNUM,FMIN,FMAX)

CHARACTER*(*) VARN

! Assign real value read in (FNUM) to the corresponding real variable
! in the NAMELIST and do a bounds check

IF(IS1.LE.MAX1.AND.IS2.LE.MAX2) THEN
  VAR = FNUM
ELSE
  PRINT 9,VARN(1:lastchar(VARN)),IS1,MAX1,IS2,MAX2,FNUM
  9 FORMAT(' -- ERROR --   Input variable - ',A,' - attempting'  &
        ,' to read extra values, values ignored.',/,  &
         ' Subscript',I6,', max dimension',I6,/,  &
         ' Subscript',I6,', max dimension',I6,', value ',F20.6)
ENDIF

IF(FNUM.LT.FMIN.OR.FNUM.GT.FMAX) THEN
  PRINT 10,VARN(1:lastchar(VARN)),FNUM,FMIN,FMAX
  10 FORMAT(' -- ERROR --   Input variable - ',A,' - set to ',F18.5  &
      ,/,'                 allowable range ',F15.5,' to ',F15.5)
ENDIF

RETURN
END

!***************************************************************************

SUBROUTINE VARSETI (VARN,IVAR,IS1,MAXSUB,INUM,IMIN,IMAX)

CHARACTER*(*) VARN

! Assign integer value read in (INUM) to the corresponding integer var
! (IVAR) in the NAMELIST and do a bounds check

IF(IS1.LE.MAXSUB) THEN
  IVAR=INUM
ELSE
  PRINT 9,VARN(1:lastchar(VARN)),IS1,MAXSUB,INUM
  9 FORMAT(' -- ERROR --   Input variable - ',A,' - attempting'  &
        ,' to read extra values, values ignored.',/,  &
         ' Subscript, max dimension, value ',2I6,I20)
ENDIF

IF(INUM.LT.IMIN.OR.INUM.GT.IMAX) THEN
  PRINT 10,VARN(1:lastchar(VARN)),INUM,IMIN,IMAX
  10 FORMAT(' -- ERROR --   Input variable - ',A,' - set to ',I10  &
      ,/,'                 allowable range ',I10,' to ',I10)
ENDIF

RETURN
END

!***************************************************************************

SUBROUTINE VARSETI2 (VARN,IVAR,IS1,MAX1,IS2,MAX2,INUM,IMIN,IMAX)

CHARACTER*(*) VARN

! Assign integer value read in (INUM) to the corresponding integer var
! (IVAR) in the NAMELIST and do a bounds check

IF(IS1.LE.MAX2.AND.IS2.LE.MAX2) THEN
  IVAR = INUM
ELSE
  PRINT 9,VARN(1:lastchar(VARN)),IS1,MAX1,IS2,MAX2,INUM
  9 FORMAT(' -- ERROR --   Input variable - ',A,' - attempting'  &
        ,' to read extra values, values ignored.',/,  &
         ' Subscript',I6,', max dimension',I6,/,  &
         ' Subscript',I6,', max dimension',I6,', value ',I10)
ENDIF

IF(INUM.LT.IMIN.OR.INUM.GT.IMAX) THEN
  PRINT 10,VARN(1:lastchar(VARN)),INUM,IMIN,IMAX
  10 FORMAT(' -- ERROR --   Input variable - ',A,' - set to ',I10  &
      ,/,'                 allowable range ',I10,' to ',I10)
ENDIF

RETURN
END

!***************************************************************************

SUBROUTINE VARSETC (VARN,VAR,IS1,MAXSUB,CH,IMIN,IMAX)

CHARACTER*(*) VARN,CH,VAR

! Assign character value read in (CH) to the corresponding character
! variable (VAR) in the NAMELIST and do a bounds check

IF(IS1.LE.MAXSUB) THEN
  VAR=CH
ELSE
  PRINT 9,VARN(1:lastchar(VARN)),IS1,MAXSUB,CH(1:lastchar(CH))
  9 FORMAT(' -- ERROR --   Input variable - ',A,' - attempting'  &
        ,' to read extra values, values ignored.',/,  &
         ' Subscript, max dimension, value ',2I6,A)
ENDIF

LCH=LEN(CH)

IF(LCH.LT.IMIN.OR.LCH.GT.IMAX) THEN
  PRINT 10,VARN(1:lastchar(VARN)),CH(1:lastchar(CH)),IMIN,IMAX
  10 FORMAT(' -- ERROR --   Input variable - ',A,' - set to ',A  &
      ,/,'                 allowable length ',I10,' to ',I10)
ENDIF

RETURN
END

!***************************************************************************

SUBROUTINE VARSETC2 (VARN,VAR,IS1,MAX1,IS2,MAX2,CH,IMIN,IMAX)

CHARACTER*(*) VARN,CH,VAR

! Assign character value read in (CH) to the corresponding character
! variable (VAR) in the NAMELIST and do a bounds check

IF(IS1.LE.MAX1.AND.IS2.LE.MAX2) THEN
  VAR = CH
ELSE
  PRINT 9,VARN(1:lastchar(VARN)),IS1,MAX1,IS2,MAX2  &
       ,CH(1:lastchar(CH))
  9 FORMAT(' -- ERROR --   Input variable - ',A,' - attempting'  &
        ,' to read extra values, values ignored.',/,  &
         ' Subscript',I6,', max dimension',I6,/,  &
         ' Subscript',I6,', max dimension',I6,', value ',A)
ENDIF

LCH=LEN(CH)

IF(LCH.LT.IMIN.OR.LCH.GT.IMAX) THEN
  PRINT 10,VARN(1:lastchar(VARN)),CH(1:lastchar(CH)),IMIN,IMAX
  10 FORMAT(' -- ERROR --   Input variable - ',A,' - set to ',A  &
      ,/,'                 allowable length ',I10,' to ',I10)
ENDIF

RETURN
END

!***************************************************************************

SUBROUTINE VARCHK (VARN,GROUP,NAMELST,INAME,NNAME,INRFLG)

DIMENSION INAME(NNAME)
CHARACTER*(*) VARN,GROUP,NAMELST(NNAME)
COMMON/NREAD/NRFLAG,NFATAL

! This routine checks that all member variables of the NAMELIST
! specified by GROUP have been assigned values.
!
! GROUP   - name of the NAMELIST being checked
! INNAME  - storage vector for counting number of times each member
!           of the NAMELIST has been assigned a value
! NAMELST - character vector containing names of members of NAMELIST
! NNAME   - number of variables in NAMELIST
! VARN    - name of NAMELIST variable being checked

INRFLG=0

IF(VARN.NE.'$END')THEN

  ! Update INAME value corresponding to the NAMELIST variable which has
  ! been assigned a value in the call to NVFILL.

  DO NN=1,NNAME
    IF(VARN.EQ.NAMELST(NN)) THEN
      INAME(NN)=1
      RETURN
    ENDIF
  enddo

  INRFLG=1
  PRINT 20,GROUP(1:lastchar(GROUP)),VARN(1:lastchar(VARN))
  20 FORMAT(' Extra variable in namelist -- ',A,'-- variable name -',A)

ELSE

  ! End of NAMELIST input has been reached - check that all variables ha
  ! been given values

  INRFLG=1
  IF(NFATAL.EQ.0) RETURN
  INOSET=0
  DO NN=1,NNAME
    IF(INAME(NN).EQ.0) THEN
      IF(INOSET.EQ.0) THEN
        PRINT 30,GROUP(1:lastchar(GROUP))
      ENDIF
      INOSET=1
      PRINT 31,NAMELST(NN)(1:lastchar(NAMELST(NN)))
    ENDIF
  enddo

  30 FORMAT(' The following variables have not been set in the -- ',A  &
           ,' -- namelist:')
  31 FORMAT('     variable name: ',A)

ENDIF

RETURN
END

!***************************************************************************

SUBROUTINE CH2INT (STR,INT)

CHARACTER*(*) STR
CHARACTER*8 FORM

! Read integer value INT from character string STR

NC=len_trim(STR)
WRITE(FORM,90)NC
90 FORMAT('(I',I2,')')
READ(STR,FORM)INT
  
RETURN
END

!***************************************************************************

SUBROUTINE CH2REAL (STR,FNUM)

CHARACTER*(*) STR
CHARACTER*8 FORM

! Read real value FNUM from character string STR

NC=len_trim(STR)
WRITE(FORM,90)NC
90 FORMAT('(F',I2,'.0)')
READ(STR,FORM)FNUM
  
RETURN
END

!***************************************************************************

SUBROUTINE CH2CH (STR,CHVAL,NCW)

CHARACTER*(*) STR,CHVAL
CHARACTER*8 FORM

! Remove trailing blanks from character string STR and store remaining
! NCW characters in character string CHVAL

NCSTR=LEN(STR)
DO 10 NC=NCSTR,1,-1
  IF(STR(NC:NC).EQ.' ') GOTO 10
  CHVAL=STR(2:NC)
  NCW=NC-2
  RETURN
10 CONTINUE

RETURN
END

!***************************************************************************

FUNCTION LETTER (STR)

CHARACTER*(*) STR

! First character alpha check - test to see if the first character of
! the string STR is alphabetic: LETTER = 0 if 'no', = 1 if 'yes'.

LETTER=0
IF((STR(1:1).GE.'A'.AND.STR(1:1).LE.'Z').or.  &
   (STR(1:1).GE.'a'.AND.STR(1:1).LE.'z')) LETTER=1
   
RETURN
END

!***************************************************************************

FUNCTION NUMBER (STR)

CHARACTER*(*) STR

! First character number check - test to see if the first character of
! the string STR is numeric:  NUMBER = 0 if 'no', = 1 if 'yes' (includ
! a decimal point or minus sign).

NUMBER=0
IF(STR(1:1).GE.'0'.AND.STR(1:1).LE.'9') NUMBER=1
IF(STR(1:1).EQ.'.'.OR.STR(1:1).EQ.'-') NUMBER=1

RETURN
END

!***************************************************************************

FUNCTION LETINT (STR)

CHARACTER*(*) STR

! First character integer variable check - test to see if the first
! character of STR is an I, J, K, L, M, or N, or i, j, k, l ,m or n:
! LETINT = 0 if 'no', = 1 if 'yes'

LETINT=0
IF((STR(1:1).GE.'I'.AND.STR(1:1).LE.'N').or.  &
   (STR(1:1).GE.'i'.AND.STR(1:1).LE.'n')) LETINT=1
   
RETURN
END

!***************************************************************************

FUNCTION LETQUO (STR)

CHARACTER*(*) STR

! First character quote check - test to see if the first character
! of STR is a quote:  LETQUO = 0 if 'no', = 1 if 'yes'.

LETQUO=0
IF(STR(1:1).EQ.'''') LETQUO=1

RETURN
END

!***************************************************************************

SUBROUTINE FINDGR (IUNIT,GROUP,MAXREC)

CHARACTER*(*) GROUP
CHARACTER*80 LINE
COMMON/NREAD/NRFLAG,NFATAL

! This routine checks to see if any of the first MAXREC lines on input
! unit IUNIT contains the character string GROUP

DO NR=1,MAXREC
  READ(IUNIT,'(A80)',END=100) LINE
  IND=INDEX(LINE,GROUP)
  IF(IND.NE.0) RETURN
enddo
100 CONTINUE
IF(NFATAL.EQ.1) THEN
PRINT *,' Namelist read error -- group not found -- ',GROUP
ENDIF

RETURN
END

!***************************************************************************

SUBROUTINE STRIP (LIN1,LIN2,NC2)

CHARACTER*(*) LIN1,LIN2

! This routine strips blank characters from character string LIN1
! (as well as comments beginning with an '!') and stores the stripped-
! down remainder in character string LIN1.  NC2 is the number of
! characters in LIN2.

NL=LEN(LIN1)

NC2=0
IQUOTE=0
DO NC=1,NL
  IF(IQUOTE.EQ.0) THEN
    IF(LIN1(NC:NC).NE.' ') THEN
      IF(LIN1(NC:NC).EQ.'!') RETURN
      NC2=NC2+1
      LIN2(NC2:NC2)=LIN1(NC:NC)
    ENDIF
  ELSE
    NC2=NC2+1
    LIN2(NC2:NC2)=LIN1(NC:NC)
  ENDIF
  IF(LIN1(NC:NC).EQ.'''') THEN
    IF(IQUOTE.EQ.0) THEN
      IQUOTE=1
    ELSE
      IQUOTE=0
    ENDIF
  ENDIF
enddo

RETURN
END

!***************************************************************************

SUBROUTINE TOKNZE (STR,NCH,TOKENS,NTOK)

PARAMETER (NSEP=4)
CHARACTER*(*) STR
CHARACTER*(*) TOKENS(*)
CHARACTER*1 TOKSEP(NSEP)

! This routine "parses" character string STR into different pieces
! or tokens by looking for one of four possible token separators (TOKS
! STR contains NCH characters.  The number of tokens identified is NTO
! the character string tokens are stored in TOKENS.

DATA TOKSEP/'=',',','(',')'/

NTOK=0
NPT=1
DO 10 NC=1,NCH
  DO 5 NS=1,NSEP
    IF(STR(NC:NC).EQ.TOKSEP(NS))THEN
      IF(NC-NPT.GE.1)THEN
        NTOK=NTOK+1
        TOKENS(NTOK)=STR(NPT:NC-1)
      ENDIF
      NTOK=NTOK+1
      TOKENS(NTOK)=STR(NC:NC)
      NPT=NC+1
      GOTO 10
    ENDIF
  5 CONTINUE
10 CONTINUE

RETURN
END

!***************************************************************************

!------------------------------------------------------------------
!
!      NAMELIST output routines
!
!------------------------------------------------------------------

SUBROUTINE FILLPG (VARNAME,IROW,ICOL,IVALUE,RVALUE,CVALUE,VARTYPE,NVALUE)

! his routine builds an output page for a listing of NAMELIST
! ariable values by inserting NAMELIST variable names and values
! nto a CHARACTER array called PAGE.  One call to this routine
! s required for each variable to be output.
!
!
! VARNAME - name of variable to be inserted on page (character stri
! IROW    - row number of page on which variable name and value are
!              to be inserted
! ICOL    - column of page in which name and value are to be insert
! IVALUE  - vector of values if variable is an integer
! RVALUE  - vector of values if variable is a real
! CVALUE  - vector of values if variable is a character string
! VARTYPE - type of variable to be written (choice of 'I', 'R', or
!              'C')
! NVALUE  - number of values in value vector

CHARACTER*(*)  CVALUE(1)
CHARACTER*132  PAGE(80)
CHARACTER*64   CDUP, FORM
CHARACTER*8    VARNAME
CHARACTER*4    NUMD
CHARACTER*1    VARTYPE, APSTRPH
INTEGER        IVALUE(1), IDUP, NVALUE, ROW, COL, COLP7, COLP9
INTEGER        ICOL, IROW, LEFTCOL, NUMDUP, OLDNUM
REAL           RVALUE(1), RDUP

COMMON /PAGEOUT/ PAGE

DATA  APSTRPH/''''/

ROW = IROW
COL = ICOL
COLP7 = COL + 7
COLP9 = COL + 9
PAGE(ROW)(COL:COLP7) = VARNAME
PAGE(ROW)(COLP9:COLP9) = '='
LEFTCOL = COL + 11
NUMDUP = 1

IF(VARTYPE .EQ. 'I') THEN

   ! Integer variable

   IDUP = IVALUE(1)
   IF(NVALUE .EQ. 1) THEN

      ! Single variable value

      WRITE(FORM,600) IDUP
      LENGTH = 8
      CALL WRITEVL(ROW, COL, LEFTCOL, FORM, LENGTH)
   ELSE

      ! More than one variable value

      DO NV = 2, NVALUE, 1
         OLDNUM = NUMDUP
         IF(IDUP .EQ. IVALUE(NV)) NUMDUP = NUMDUP + 1
         IF(NUMDUP .EQ. OLDNUM) THEN
            IF(NUMDUP .EQ. 1) THEN
               WRITE(FORM,600) IDUP
               LENGTH = 10
            ELSE
               WRITE(FORM,650) NUMDUP, IDUP
               LENGTH = 13
            ENDIF
            CALL WRITEVL(ROW, COL, LEFTCOL, FORM, LENGTH)
            IDUP = IVALUE(NV)
            NUMDUP = 1
         ENDIF
      enddo
      WRITE(FORM,650) NUMDUP, IDUP
      LENGTH = 11
      CALL WRITEVL(ROW, COL, LEFTCOL, FORM, LENGTH)
   ENDIF

ELSEIF(VARTYPE .EQ. 'R') THEN

   ! Real variable

   RDUP = RVALUE(1)
   IF(NVALUE .EQ. 1) THEN

      ! Single variable value

      WRITE(FORM,610) RDUP
      LENGTH = 10
      CALL WRITEVL(ROW, COL, LEFTCOL, FORM, LENGTH)
   ELSE

      ! More than one variable value

      DO NV = 2, NVALUE, 1
         OLDNUM = NUMDUP
         IF(RDUP .EQ. RVALUE(NV)) NUMDUP = NUMDUP + 1
         IF(NUMDUP .EQ. OLDNUM) THEN
            IF(NUMDUP .EQ. 1) THEN
               WRITE(FORM,610) RDUP
               LENGTH = 10
            ELSE
               WRITE(FORM,660) NUMDUP, RDUP
               LENGTH = 13
            ENDIF
            CALL WRITEVL(ROW, COL, LEFTCOL, FORM, LENGTH)
            RDUP = RVALUE(NV)
            NUMDUP = 1
         ENDIF
      enddo
      WRITE(FORM,660) NUMDUP, RDUP
      LENGTH = 13
      CALL WRITEVL(ROW, COL, LEFTCOL, FORM, LENGTH)
   ENDIF

ELSE

   ! Character variable

   CDUP = CVALUE(1)
   IF(NVALUE .EQ. 1) THEN

      ! Single variable value

      LENGTH = len_trim(CDUP)
      !LENGTH = INDEX(CDUP, '        ')
      FORM = APSTRPH//CDUP(1:LENGTH)//APSTRPH
      LENGTH = LENGTH + 2
      CALL WRITEVL(ROW, COL, LEFTCOL, FORM, LENGTH)
      
   ELSE

      ! More than one variable value

      DO NV = 2, NVALUE, 1
         OLDNUM = NUMDUP
         IF(CDUP .EQ. CVALUE(NV)) NUMDUP = NUMDUP + 1
         IF(NUMDUP .EQ. OLDNUM) THEN
            LENGTH = len_trim(CDUP)
            !LENGTH = INDEX(CDUP, '        ')
            IF(CDUP(1:1) .EQ. ' ') LENGTH = 8
            IF(NUMDUP .EQ. 1) THEN
               FORM = APSTRPH//CDUP(1:LENGTH)//APSTRPH
               LENGTH = LENGTH + 2
            ELSE
               WRITE(NUMD,670) NUMDUP
               FORM = NUMD//CDUP(1:LENGTH)//APSTRPH
               LENGTH = LENGTH + 5
            ENDIF
            CALL WRITEVL(ROW, COL, LEFTCOL, FORM, LENGTH)
            CDUP = CVALUE(NV)
            NUMDUP = 1
         ENDIF
      enddo
      LENGTH = len_trim(CDUP)
      !LENGTH = INDEX(CDUP, '        ')
      IF(CDUP(1:1) .EQ. ' ') LENGTH = 8
      WRITE(NUMD,670) NUMDUP
      FORM = NUMD//CDUP(1:LENGTH)//APSTRPH
      LENGTH = LENGTH + 5
      CALL WRITEVL(ROW, COL, LEFTCOL, FORM, LENGTH)
   ENDIF

ENDIF

RETURN

600 FORMAT(I8)
610 FORMAT(G10.3)
!  650 FORMAT(I2,1H*,I8)
!  660 FORMAT(I2,1H*,G10.3)
!  670 FORMAT(I2,2H*')
650 FORMAT(I2,'*',I8)
660 FORMAT(I2,'*',G10.3)
670 FORMAT(I2,'*')

END

!***************************************************************************

SUBROUTINE WRITEVL (ROW, COL, LEFTCOL, FORM, LENGTH)

! This routine inserts one NAMELIST variable value into the
! character array PAGE.  ROW and COL refer the the location in
! PAGE of the first value of the variable being output.  LEFTCOL
! gives the first column for any one of the values of the variable
! being inserted in PAGE.

CHARACTER*132 PAGE(80)
CHARACTER*64  FORM
CHARACTER*1   COMMA
INTEGER       COL, LEFTCOL, LENGTH, RGHTCOL, ROW

COMMON /PAGEOUT/ PAGE

DATA  COMMA/','/

RGHTCOL = LEFTCOL + LENGTH - 1
IF(ROW .LT. 1  .OR.  ROW .GT. 80) STOP 'WRITEVL'
PAGE(ROW)(LEFTCOL:RGHTCOL) = FORM(1:LENGTH)
LEFTCOL = RGHTCOL + 2
PAGE(ROW)(LEFTCOL:LEFTCOL) = COMMA
LEFTCOL = LEFTCOL + 5
IF((LEFTCOL+LENGTH) .GT. 90) THEN
   LEFTCOL = COL + 11
   ROW = ROW + 1
ENDIF

RETURN
END

!***************************************************************************

SUBROUTINE BLANKPG

! This routine sets the output page character variable to all blanks

CHARACTER*132 PAGE(80)
CHARACTER*61 BLANKLN

COMMON /PAGEOUT/ PAGE

DATA BLANKLN/' '/

DO I=1,80
   PAGE(I)=BLANKLN
ENDDO

RETURN
END

!***************************************************************************

SUBROUTINE WRITEPG (GROUP,NUMROW)

! This routine writes out the character array containing the list of
! variable names and values for a specified NAMELIST given by  GROUP

CHARACTER*132 PAGE(80)
CHARACTER*(*) GROUP
CHARACTER*9 FMT
INTEGER NR,NUMROW

COMMON  /PAGEOUT/ PAGE

WRITE(6,600) GROUP
DO NR = 1, NUMROW
   DO NCH=132,1,-1
     IF(PAGE(NR)(NCH:NCH).NE.' ') GOTO 110
   ENDDO
   NCH=1
   110 CONTINUE
   WRITE(FMT,630) NCH
   630 FORMAT('(1X,A',I3,')')
   WRITE(6,FMT) PAGE(NR)(1:NCH)
ENDDO
WRITE(6,620)

600 FORMAT(///// ' NAMELIST  ', A7 / 1X, 120('-') /)
620 FORMAT(1X, 120('-') )

RETURN
END
