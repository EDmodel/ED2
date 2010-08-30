!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

subroutine stable(ityp,key,xmu,tau,beta,omeg)
  implicit none
  integer :: ityp,key
  real :: xmu,tau,beta,omeg

  ! PROVIDES TWO-STREAM MODEL PARAMETERS FOR A CLOUD LAYER IN THE
  ! UV-VIS OR NEAR-IR SPECTRUM ( I.E. BACKSCATTER FRACTION AND SINGLE
  ! SCATTER ALBEDO )
  !
  ! FOLLOWING
  !
  ! STEPHENS,G.L.,1978.
  ! RADIATION PROFILES IN EXTENDED WATER CLOUDS. II.PARAMETERIZATION S
  ! JOUR. ATMOS. SCI.,35,2123-2132.
  !
  ! BI-DIMENSIONAL LINEAR OR QUADRATIC INTERPOLATION SCHEME ( LAGRANGE
  !
  ! INPUT PARAMETERS
  !
  ! ITYP = 1 FOR     CONSERVATIVE MEDIUM (  UV-VIS SPECTRUM ) - OMEG E
  !      = 2 FOR NON-CONSERVATIVE MEDIUM ( NEAR-IR SPECTRUM ) - OMEG N
  ! KEY = 1 FOR LINEAR INTERPOLATION
  !     = 2 FOR QUADRATIC INTERPOLATION
  ! XMU = COSINE OF SOLAR ZENITH ANGLE
  ! TAU = OPTICAL DEPTH
  !
  ! OUTPUT PARAMETERS
  !
  ! BETA = BACKSCATTER FRACTION
  ! OMEG = SINGLE SCATTER ALBEDO

  real :: LU0,LU1,LU2,LT0,LT1,LT2
  real :: UTAB(10),TTAB(12),TAB(10,12,3)

  integer :: L ! Necessary to Intel Compiler - Alvaro L.Fazenda

  data UTAB/1.0,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0.03/
  data TTAB/1.0,2.0,5.0,10.0,16.0,25.0,40.0,60.0,80.0,100.0,200.0,  &
       500.0/
  data (TAB(L, 1,1),L=1,10)/  &
       .0421,.0557,.0657,.0769,.0932,.1111,.1295,.1407,.1196,.1053/
  data (TAB(L, 2,1),L=1,10)/  &
       .0472,.0615,.0708,.0803,.0924,.1017,.1077,.1034,.0794,.0624/
  data (TAB(L, 3,1),L=1,10)/  &
       .0582,.0692,.0744,.0782,.0815,.0812,.0776,.0680,.0483,.0325/
  data (TAB(L, 4,1),L=1,10)/  &
       .0681,.0726,.0737,.0733,.0723,.0685,.0626,.0527,.0359,.0230/
  data (TAB(L, 5,1),L=1,10)/  &
       .0734,.0738,.0728,.0707,.0680,.0631,.0564,.0465,.0310,.0193/
  data (TAB(L, 6,1),L=1,10)/  &
       .0768,.0744,.0723,.0691,.0653,.0598,.0526,.0427,.0281,.0171/
  data (TAB(L, 7,1),L=1,10)/  &
       .0791,.0749,.0719,.0680,.0636,.0575,.0501,.0402,.0261,.0155/
  data (TAB(L, 8,1),L=1,10)/  &
       .0805,.0752,.0717,.0674,.0627,.0563,.0488,.0389,.0251,.0147/
  data (TAB(L, 9,1),L=1,10)/  &
       .0812,.0754,.0717,.0672,.0622,.0558,.0481,.0382,.0246,.0144/
  data (TAB(L,10,1),L=1,10)/  &
       .0820,.0757,.0717,.0670,.0619,.0553,.0475,.0376,.0241,.0141/
  data (TAB(L,11,1),L=1,10)/  &
       .0831,.0763,.0721,.0672,.0619,.0552,.0473,.0374,.0241,.0141/
  data (TAB(L,12,1),L=1,10)/  &
       .0874,.0800,.0755,.0703,.0647,.0576,.0494,.0392,.0262,.0163/
  data (TAB(L, 1,2),L=1,10)/  &
       .0477,.0627,.0734,.0855,.1022,.1200,.1379,.1465,.1207,.1119/
  data (TAB(L, 2,2),L=1,10)/  &
       .0537,.0690,.0788,.0886,.1003,.1090,.1133,.1065,.0794,.0592/
  data (TAB(L, 3,2),L=1,10)/  &
       .0660,.0769,.0817,.0850,.0871,.0864,.0801,.0688,.0474,.0322/
  data (TAB(L, 4,2),L=1,10)/  &
       .0759,.0793,.0795,.0781,.0757,.0705,.0629,.0516,.0339,.0202/
  data (TAB(L, 5,2),L=1,10)/  &
       .0801,.0787,.0766,.0732,.0689,.0626,.0543,.0434,.0277,.0157/
  data (TAB(L, 6,2),L=1,10)/  &
       .0807,.0759,.0724,.0678,.0625,.0555,.0471,.0368,.0229,.0125/
  data (TAB(L, 7,2),L=1,10)/  &
       .0770,.0700,.0656,.0603,.0545,.0476,.0396,.0302,.0184,.0096/
  data (TAB(L, 8,2),L=1,10)/  &
       .0699,.0621,.0575,.0522,.0466,.0401,.0329,.0248,.0148,.0075/
  data (TAB(L, 9,2),L=1,10)/  &
       .0634,.0556,.0510,.0460,.0408,.0348,.0283,.0211,.0125,.0062/
  data (TAB(L,10,2),L=1,10)/  &
       .0534,.0461,.0420,.0376,.0330,.0279,.0225,.0166,.0097,.0048/
  data (TAB(L,11,2),L=1,10)/  &
       .0415,.0353,.0319,.0283,.0246,.0206,.0165,.0120,.0068,.0032/
  data (TAB(L,12,2),L=1,10)/  &
       .0251,.0208,.0186,.0163,.0140,.0115,.0090,.0064,.0032,.0011/
  data (TAB(L, 1,3),L=1,10)/  &
       .0225,.0222,.0218,.0208,.0199,.0155,.0109,.0059,.0017,.0002/
  data (TAB(L, 2,3),L=1,10)/  &
       .0213,.0200,.0179,.0176,.0156,.0118,.0078,.0038,.0010,.0001/
  data (TAB(L, 3,3),L=1,10)/  &
       .0195,.0166,.0146,.0125,.0096,.0069,.0043,.0021,.0005,.0000/
  data (TAB(L, 4,3),L=1,10)/  &
       .0173,.0138,.0114,.0093,.0070,.0049,.0026,.0013,.0003,.0000/
  data (TAB(L, 5,3),L=1,10)/  &
       .0156,.0111,.0090,.0073,.0052,.0035,.0019,.0009,.0002,.0000/
  data (TAB(L, 6,3),L=1,10)/  &
       .0115,.0088,.0069,.0052,.0038,.0026,.0014,.0007,.00014,.000/
  data (TAB(L, 7,3),L=1,10)/  &
       .0104,.0070,.0053,.0040,.0028,.0020,.0011,.0005,.0001,.0000/
  data (TAB(L, 8,3),L=1,10)/  &
       .0083,.0055,.00425,.0032,.0023,.00145,.0008,.0003,.0001,.00/
  data (TAB(L, 9,3),L=1,10)/  &
       .0069,.0050,.0038,.0028,.0020,.0013,.0007,.00034,.0000,.000/
  data (TAB(L,10,3),L=1,10)/  &
       .0060,.0043,.0035,.0022,.0018,.0011,.0006,.0003,.0000,.0000/
  data (TAB(L,11,3),L=1,10)/  &
       .0044,.0031,.0025,.0016,.0011,.00072,.0004,.00019,.0000,.00/
  data (TAB(L,12,3),L=1,10)/  &
       .0026,.0018,.0014,.0010,.00072,.00048,.00029,.00015,.00,.00/
  !INTEGER :: l,ii,i0,i1,i2,m,jj,j0,j1,j2,k1,k2,n
  integer :: ii,i0,i1,i2,m,jj,j0,j1,j2,k1,k2,n
  real :: u,t,u0,u1,t0,t1,p0,p1,p,u2,t2,p2

  ! BOUNDS CHECK

  U=XMU
  T=TAU
  if(U.gt.UTAB(1)) U=UTAB(1)
  if(U.lt.UTAB(10)) U=UTAB(10)
  if(T.lt.TTAB(1)) T=TTAB(1)
  if(T.gt.TTAB(12)) T=TTAB(12)

  ! SELECT THE XMU INDICIES

  do L=2,10
     II=L
     if(U.ge.UTAB(L)) goto 2
  enddo
2 continue
  if(KEY.ne.1) goto 3
  I0=II-1
  I1=II
  goto 5
3 continue
  if(II.eq.10) goto 4
  I0=II-1
  I1=II
  I2=II+1
  goto 5
4 continue
  I0=II-2
  I1=II-1
  I2=II

  ! SELECT THE TAU INDICIES

5 continue

  do M=2,12
     JJ=M
     if(T.le.TTAB(M)) goto 7
  enddo
7 continue
  if(KEY.ne.1) goto 8
  J0=JJ-1
  J1=JJ
  goto 10
8 continue
  if(JJ.eq.12) goto 9
  J0=JJ-1
  J1=JJ
  J2=JJ+1
  goto 10
9 continue
  J0=JJ-2
  J1=JJ-1
  J2=JJ

  ! SET PARAMETER INDICIES

10 continue
  if(ITYP.eq.1 )K1=1
  if(ITYP.eq.1) K2=1
  if(ITYP.eq.2) K1=2
  if(ITYP.eq.2) K2=2
  if(ITYP.eq.3) K1=3
  if(ITYP.eq.3) K2=3

  ! BRANCH ON SCHEME TYPE

  if(KEY.ne.1) goto 20

  ! LINEAR SCHEME

  ! DESIGNATE INDEX VALUES

  U0=UTAB(I0)
  U1=UTAB(I1)
  T0=TTAB(J0)
  T1=TTAB(J1)

  ! LAGRANGE POLYNOMIALS

  LU0=(U-U1)/(U0-U1)
  LU1=(U-U0)/(U1-U0)
  LT0=(T-T1)/(T0-T1)
  LT1=(T-T0)/(T1-T0)

  ! LOOP OVER THE SCATTERING PARAMETER INDEX

  do N=K1,K2

     ! INTERPOLATING POLYNOMIALS FOR THE FIRST DIMENSION

     P0=TAB(I0,J0,N)*LU0+TAB(I1,J0,N)*LU1
     P1=TAB(I0,J1,N)*LU0+TAB(I1,J1,N)*LU1

     ! INTERPOLATING POLYNOMIAL FOR SECOND DIMENSION

     P=P0*LT0+P1*LT1

     ! ASSIGN THE INDIVIDUAL PARAMETERS

     if(N.eq.1)BETA=P
     if(N.eq.1)OMEG=0
     if(N.eq.2)BETA=P
     if(N.eq.3)OMEG=1-P
  enddo

  return

  ! QUADRATIC SCHEME

  ! DESIGNATE INDEX VALUES

20 continue

  U0=UTAB(I0)
  U1=UTAB(I1)
  U2=UTAB(I2)
  T0=TTAB(J0)
  T1=TTAB(J1)
  T2=TTAB(J2)

  ! LAGRANGE POLYNOMIALS

  LU0=(U-U1)*(U-U2)/((U0-U1)*(U0-U2))
  LU1=(U-U0)*(U-U2)/((U1-U0)*(U1-U2))
  LU2=(U-U0)*(U-U1)/((U2-U0)*(U2-U1))
  LT0=(T-T1)*(T-T2)/((T0-T1)*(T0-T2))
  LT1=(T-T0)*(T-T2)/((T1-T0)*(T1-T2))
  LT2=(T-T0)*(T-T1)/((T2-T0)*(T2-T1))

  ! LOOP OVER THE SCATTERING PARAMETER INDEX

  do N=K1,K2

     ! INTERPOLATING POLYNOMIALS FOR THE FIRST DIMENSION

     P0=TAB(I0,J0,N)*LU0+TAB(I1,J0,N)*LU1+TAB(I2,J0,N)*LU2
     P1=TAB(I0,J1,N)*LU0+TAB(I1,J1,N)*LU1+TAB(I2,J1,N)*LU2
     P2=TAB(I0,J2,N)*LU0+TAB(I1,J2,N)*LU1+TAB(I2,J2,N)*LU2

     ! INTERPOLATING POLYNOMIAL FOR SECOND DIMENSION

     P=P0*LT0+P1*LT1+P2*LT2

     ! ASSIGN THE INDIVIDUAL PARAMETERS

     if(N.eq.1) BETA=P
     if(N.eq.1) OMEG=0
     if(N.eq.2) BETA=P
     if(N.eq.3) OMEG=1-P

  enddo

  return
end subroutine stable

![MLO - Double precision version of stable
subroutine dpstable(ityp,key,xmu,tau,beta,omeg)
  implicit none
  integer :: ityp,key
  real(kind=8) :: xmu,tau,beta,omeg

  ! PROVIDES TWO-STREAM MODEL PARAMETERS FOR A CLOUD LAYER IN THE
  ! UV-VIS OR NEAR-IR SPECTRUM ( I.E. BACKSCATTER FRACTION AND SINGLE
  ! SCATTER ALBEDO )
  !
  ! FOLLOWING
  !
  ! STEPHENS,G.L.,1978.
  ! RADIATION PROFILES IN EXTENDED WATER CLOUDS. II.PARAMETERIZATION S
  ! JOUR. ATMOS. SCI.,35,2123-2132.
  !
  ! BI-DIMENSIONAL LINEAR OR QUADRATIC INTERPOLATION SCHEME ( LAGRANGE
  !
  ! INPUT PARAMETERS
  !
  ! ITYP = 1 FOR     CONSERVATIVE MEDIUM (  UV-VIS SPECTRUM ) - OMEG E
  !      = 2 FOR NON-CONSERVATIVE MEDIUM ( NEAR-IR SPECTRUM ) - OMEG N
  ! KEY = 1 FOR LINEAR INTERPOLATION
  !     = 2 FOR QUADRATIC INTERPOLATION
  ! XMU = COSINE OF SOLAR ZENITH ANGLE
  ! TAU = OPTICAL DEPTH
  !
  ! OUTPUT PARAMETERS
  !
  ! BETA = BACKSCATTER FRACTION
  ! OMEG = SINGLE SCATTER ALBEDO

  real(kind=8) :: LU0,LU1,LU2,LT0,LT1,LT2
  real(kind=8) :: UTAB(10),TTAB(12),TAB(10,12,3)

  integer :: L ! Necessary to Intel Compiler - Alvaro L.Fazenda

  data UTAB/1.0,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0.03/
  data TTAB/1.0,2.0,5.0,10.0,16.0,25.0,40.0,60.0,80.0,100.0,200.0,  &
       500.0/
  data (TAB(L, 1,1),L=1,10)/  &
       .0421,.0557,.0657,.0769,.0932,.1111,.1295,.1407,.1196,.1053/
  data (TAB(L, 2,1),L=1,10)/  &
       .0472,.0615,.0708,.0803,.0924,.1017,.1077,.1034,.0794,.0624/
  data (TAB(L, 3,1),L=1,10)/  &
       .0582,.0692,.0744,.0782,.0815,.0812,.0776,.0680,.0483,.0325/
  data (TAB(L, 4,1),L=1,10)/  &
       .0681,.0726,.0737,.0733,.0723,.0685,.0626,.0527,.0359,.0230/
  data (TAB(L, 5,1),L=1,10)/  &
       .0734,.0738,.0728,.0707,.0680,.0631,.0564,.0465,.0310,.0193/
  data (TAB(L, 6,1),L=1,10)/  &
       .0768,.0744,.0723,.0691,.0653,.0598,.0526,.0427,.0281,.0171/
  data (TAB(L, 7,1),L=1,10)/  &
       .0791,.0749,.0719,.0680,.0636,.0575,.0501,.0402,.0261,.0155/
  data (TAB(L, 8,1),L=1,10)/  &
       .0805,.0752,.0717,.0674,.0627,.0563,.0488,.0389,.0251,.0147/
  data (TAB(L, 9,1),L=1,10)/  &
       .0812,.0754,.0717,.0672,.0622,.0558,.0481,.0382,.0246,.0144/
  data (TAB(L,10,1),L=1,10)/  &
       .0820,.0757,.0717,.0670,.0619,.0553,.0475,.0376,.0241,.0141/
  data (TAB(L,11,1),L=1,10)/  &
       .0831,.0763,.0721,.0672,.0619,.0552,.0473,.0374,.0241,.0141/
  data (TAB(L,12,1),L=1,10)/  &
       .0874,.0800,.0755,.0703,.0647,.0576,.0494,.0392,.0262,.0163/
  data (TAB(L, 1,2),L=1,10)/  &
       .0477,.0627,.0734,.0855,.1022,.1200,.1379,.1465,.1207,.1119/
  data (TAB(L, 2,2),L=1,10)/  &
       .0537,.0690,.0788,.0886,.1003,.1090,.1133,.1065,.0794,.0592/
  data (TAB(L, 3,2),L=1,10)/  &
       .0660,.0769,.0817,.0850,.0871,.0864,.0801,.0688,.0474,.0322/
  data (TAB(L, 4,2),L=1,10)/  &
       .0759,.0793,.0795,.0781,.0757,.0705,.0629,.0516,.0339,.0202/
  data (TAB(L, 5,2),L=1,10)/  &
       .0801,.0787,.0766,.0732,.0689,.0626,.0543,.0434,.0277,.0157/
  data (TAB(L, 6,2),L=1,10)/  &
       .0807,.0759,.0724,.0678,.0625,.0555,.0471,.0368,.0229,.0125/
  data (TAB(L, 7,2),L=1,10)/  &
       .0770,.0700,.0656,.0603,.0545,.0476,.0396,.0302,.0184,.0096/
  data (TAB(L, 8,2),L=1,10)/  &
       .0699,.0621,.0575,.0522,.0466,.0401,.0329,.0248,.0148,.0075/
  data (TAB(L, 9,2),L=1,10)/  &
       .0634,.0556,.0510,.0460,.0408,.0348,.0283,.0211,.0125,.0062/
  data (TAB(L,10,2),L=1,10)/  &
       .0534,.0461,.0420,.0376,.0330,.0279,.0225,.0166,.0097,.0048/
  data (TAB(L,11,2),L=1,10)/  &
       .0415,.0353,.0319,.0283,.0246,.0206,.0165,.0120,.0068,.0032/
  data (TAB(L,12,2),L=1,10)/  &
       .0251,.0208,.0186,.0163,.0140,.0115,.0090,.0064,.0032,.0011/
  data (TAB(L, 1,3),L=1,10)/  &
       .0225,.0222,.0218,.0208,.0199,.0155,.0109,.0059,.0017,.0002/
  data (TAB(L, 2,3),L=1,10)/  &
       .0213,.0200,.0179,.0176,.0156,.0118,.0078,.0038,.0010,.0001/
  data (TAB(L, 3,3),L=1,10)/  &
       .0195,.0166,.0146,.0125,.0096,.0069,.0043,.0021,.0005,.0000/
  data (TAB(L, 4,3),L=1,10)/  &
       .0173,.0138,.0114,.0093,.0070,.0049,.0026,.0013,.0003,.0000/
  data (TAB(L, 5,3),L=1,10)/  &
       .0156,.0111,.0090,.0073,.0052,.0035,.0019,.0009,.0002,.0000/
  data (TAB(L, 6,3),L=1,10)/  &
       .0115,.0088,.0069,.0052,.0038,.0026,.0014,.0007,.00014,.000/
  data (TAB(L, 7,3),L=1,10)/  &
       .0104,.0070,.0053,.0040,.0028,.0020,.0011,.0005,.0001,.0000/
  data (TAB(L, 8,3),L=1,10)/  &
       .0083,.0055,.00425,.0032,.0023,.00145,.0008,.0003,.0001,.00/
  data (TAB(L, 9,3),L=1,10)/  &
       .0069,.0050,.0038,.0028,.0020,.0013,.0007,.00034,.0000,.000/
  data (TAB(L,10,3),L=1,10)/  &
       .0060,.0043,.0035,.0022,.0018,.0011,.0006,.0003,.0000,.0000/
  data (TAB(L,11,3),L=1,10)/  &
       .0044,.0031,.0025,.0016,.0011,.00072,.0004,.00019,.0000,.00/
  data (TAB(L,12,3),L=1,10)/  &
       .0026,.0018,.0014,.0010,.00072,.00048,.00029,.00015,.00,.00/
  !INTEGER :: l,ii,i0,i1,i2,m,jj,j0,j1,j2,k1,k2,n
  integer :: ii,i0,i1,i2,m,jj,j0,j1,j2,k1,k2,n
  real(kind=8) :: u,t,u0,u1,t0,t1,p0,p1,p,u2,t2,p2

  ! BOUNDS CHECK

  U=XMU
  T=TAU
  if(U.gt.UTAB(1)) U=UTAB(1)
  if(U.lt.UTAB(10)) U=UTAB(10)
  if(T.lt.TTAB(1)) T=TTAB(1)
  if(T.gt.TTAB(12)) T=TTAB(12)

  ! SELECT THE XMU INDICIES

  do L=2,10
     II=L
     if(U.ge.UTAB(L)) goto 2
  enddo
2 continue
  if(KEY.ne.1) goto 3
  I0=II-1
  I1=II
  goto 5
3 continue
  if(II.eq.10) goto 4
  I0=II-1
  I1=II
  I2=II+1
  goto 5
4 continue
  I0=II-2
  I1=II-1
  I2=II

  ! SELECT THE TAU INDICIES

5 continue

  do M=2,12
     JJ=M
     if(T.le.TTAB(M)) goto 7
  enddo
7 continue
  if(KEY.ne.1) goto 8
  J0=JJ-1
  J1=JJ
  goto 10
8 continue
  if(JJ.eq.12) goto 9
  J0=JJ-1
  J1=JJ
  J2=JJ+1
  goto 10
9 continue
  J0=JJ-2
  J1=JJ-1
  J2=JJ

  ! SET PARAMETER INDICIES

10 continue
  if(ITYP.eq.1 )K1=1
  if(ITYP.eq.1) K2=1
  if(ITYP.eq.2) K1=2
  if(ITYP.eq.2) K2=2
  if(ITYP.eq.3) K1=3
  if(ITYP.eq.3) K2=3

  ! BRANCH ON SCHEME TYPE

  if(KEY.ne.1) goto 20

  ! LINEAR SCHEME

  ! DESIGNATE INDEX VALUES

  U0=UTAB(I0)
  U1=UTAB(I1)
  T0=TTAB(J0)
  T1=TTAB(J1)

  ! LAGRANGE POLYNOMIALS

  LU0=(U-U1)/(U0-U1)
  LU1=(U-U0)/(U1-U0)
  LT0=(T-T1)/(T0-T1)
  LT1=(T-T0)/(T1-T0)

  ! LOOP OVER THE SCATTERING PARAMETER INDEX

  do N=K1,K2

     ! INTERPOLATING POLYNOMIALS FOR THE FIRST DIMENSION

     P0=TAB(I0,J0,N)*LU0+TAB(I1,J0,N)*LU1
     P1=TAB(I0,J1,N)*LU0+TAB(I1,J1,N)*LU1

     ! INTERPOLATING POLYNOMIAL FOR SECOND DIMENSION

     P=P0*LT0+P1*LT1

     ! ASSIGN THE INDIVIDUAL PARAMETERS

     if(N.eq.1)BETA=P
     if(N.eq.1)OMEG=0
     if(N.eq.2)BETA=P
     if(N.eq.3)OMEG=1-P
  enddo

  return

  ! QUADRATIC SCHEME

  ! DESIGNATE INDEX VALUES

20 continue

  U0=UTAB(I0)
  U1=UTAB(I1)
  U2=UTAB(I2)
  T0=TTAB(J0)
  T1=TTAB(J1)
  T2=TTAB(J2)

  ! LAGRANGE POLYNOMIALS

  LU0=(U-U1)*(U-U2)/((U0-U1)*(U0-U2))
  LU1=(U-U0)*(U-U2)/((U1-U0)*(U1-U2))
  LU2=(U-U0)*(U-U1)/((U2-U0)*(U2-U1))
  LT0=(T-T1)*(T-T2)/((T0-T1)*(T0-T2))
  LT1=(T-T0)*(T-T2)/((T1-T0)*(T1-T2))
  LT2=(T-T0)*(T-T1)/((T2-T0)*(T2-T1))

  ! LOOP OVER THE SCATTERING PARAMETER INDEX

  do N=K1,K2

     ! INTERPOLATING POLYNOMIALS FOR THE FIRST DIMENSION

     P0=TAB(I0,J0,N)*LU0+TAB(I1,J0,N)*LU1+TAB(I2,J0,N)*LU2
     P1=TAB(I0,J1,N)*LU0+TAB(I1,J1,N)*LU1+TAB(I2,J1,N)*LU2
     P2=TAB(I0,J2,N)*LU0+TAB(I1,J2,N)*LU1+TAB(I2,J2,N)*LU2

     ! INTERPOLATING POLYNOMIAL FOR SECOND DIMENSION

     P=P0*LT0+P1*LT1+P2*LT2

     ! ASSIGN THE INDIVIDUAL PARAMETERS

     if(N.eq.1) BETA=P
     if(N.eq.1) OMEG=0
     if(N.eq.2) BETA=P
     if(N.eq.3) OMEG=1-P

  enddo

  return
end subroutine dpstable
