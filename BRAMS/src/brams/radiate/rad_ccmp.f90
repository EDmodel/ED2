!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

!«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»!
![MLO - Changed many things here to double precision to avoid FPE. Changed dimensions to attempt binary reproducibility.
subroutine shradc(nzp,rvr,rtr,dn0r,dzzr,prd,albedo,solar,cosz,fthr,rshort)
  use rconstants, only: cp
  implicit none                  
!----- List of arguments ------------------------------------------------------------------------------------------------!
  integer     , intent(in)                   :: nzp
  real        , intent(in)                   :: albedo,cosz,solar
  real        , intent(out)                  :: rshort
  real        , intent(in)  , dimension(nzp) :: rvr,rtr,dn0r,prd, dzzr
  real        , intent(out) , dimension(nzp) :: fthr

!----- List of parameters -----------------------------------------------------------------------------------------------!
  real(kind=8), parameter                    :: zero=dble(0.),one=dble(1.)
  real(kind=8), parameter                    :: cpcgs=dble(10000.)*dble(cp)
  integer     , parameter                    :: nzmax=200 , iiv1=  1 ,  iiv2=  2 , iiv3= 3
  integer     , parameter                    ::   iv1=  1 ,  iv2=  2 ,   iv3=  3 ,  iv4= 4  ,  iv5=  5 &
                                               ,  iv6=  6 ,  iv7=  7 ,   iv8=  8 ,  iv9= 9  , iv10= 10 &
                                               , iv11= 11 , iv12= 12 ,  iv13= 13 , iv14= 14 , iv15= 15 &
                                               , iv16= 16 , iv17= 17 ,  iv18= 18 , iv19= 19 , iv20= 20 &
                                               , iv21= 21 , iv22= 22 ,  iv23= 23 , iv24= 24 , iv25= 25 &
                                               , iv26= 26 , iv27= 27 ,  iv28= 28 , iaoz= 29 , iv30= 30 &
                                               , iv31= 31 , iv32= 32 ,  irea= 33 , itr1= 34 , itr2= 35 &
                                               , iv36= 36 , iab2= 37 , irefs= 38 , ireb= 39 , iv40= 40
  real(kind=8), parameter, dimension(8)      :: dsfct1=(/0.00004, 0.002, 0.035, 0.377,  1.95,  9.40,  44.6, 190.0/)
  real(kind=8), parameter, dimension(8)      :: dsfct2=(/  0.647,0.0698,0.1443,0.0584,0.0335,0.0225,0.0158,0.0087/)

!----- Local variables ---------------------------------------------------------------------------------------------!
  real(kind=8)                               :: dalbedo,dcosz,dsolar
  real(kind=8),         dimension(nzp)       :: drvr,drtr,ddn0r,dprd,ddzzr
  real(kind=8),         dimension(nzp+1,40)  :: dsc
  real(kind=8)                               :: o3(nzmax)
  real(kind=8)                               :: dalbray   ! effective albedo including rayleigh scatter

  integer                                    :: nzpp,nz,k,k1,nbnd,icldfl,iii
  real(kind=8)                               :: dradc1,dpfct,drabar,drabarbar,dtrsmt
!----- Functions ---------------------------------------------------------------------------------------------------!
  real(kind=8), external                     :: dssum,dcvmgp,dcvmgm
!-------------------------------------------------------------------------------------------------------------------!


  nzpp=nzp+1
  nz=nzp-1

![MLO - Passing all information to double precision
  dalbedo=dble(albedo)
  dcosz=dble(cosz)
  dsolar=dble(solar)
  do k=1,nzp
    drvr(k) = dble(rvr(k) )
    drtr(k) = dble(rtr(k) )
    ddn0r(k)= dble(dn0r(k))
    dprd(k) = dble(prd(k) )
    ddzzr(k)= dble(dzzr(k))
    do iii=1,40
      dsc(k,iii) = 0.
    end do
  end do


  dsc(1,iv1)=0.
  do k=2,nzp
     dsc(k,iv1)=dsc(k-1,iv1)+1./ddzzr(k)
  enddo

  dsc(1,iv11)=0.
  do k=2,nzpp
     dsc(k,iv11)=.4+.4*dexp(dble(-4.))/(1.+dexp((dsc(nzp+2-k,iv1)  &
          *1.e-5-20.)/5.))
     o3(nzp+3-k)=dsc(k,iv11)-dsc(k-1,iv11)
  enddo

  !     Compute ozone absorptance
  !     REF.....LACIS,HANSEN,1974,J.A.S. P118

  dradc1=35./sqrt(1224.*dcosz*dcosz+1.)
  dsc(1,iv2)=0.
  dsc(1,iv3)=0.
  dsc(1,iv4)=0.
  dsc(1,iv5)=0.
  do k=2,nzpp
     dsc(k,iv2)=dradc1*o3(nzp-k+3)
  enddo
  do k=2,nzpp
     dsc(k,iv3)=dssum(k,dsc(1,iv2),1)
     dsc(k,iv4)=.02118*dsc(k,iv3)/(1.+.042*dsc(k,iv3)  &
          +.000323*dsc(k,iv3)*dsc(k,iv3))
     dsc(k,iv5)=1.082*dsc(k,iv3)/((1.+138.6*dsc(k,iv3))**.805)  &
          +.0658*dsc(k,iv3)/(1.+(103.6*dsc(k,iv3))**3)
  enddo
  do k=2,nzpp
     dsc(nzp-k+3,iaoz)=dsc(k,iv4)-dsc(k-1,iv4)+dsc(k,iv5)-dsc(k-1,iv5)
  enddo
  !
  !     Precomputation of reflectance,transmittance,absorptance in cloudy
  !     REF.....STEPHENS,1978,J.A.S.P2123
  !
  !     --- Cloud fractional coverage
  !  Old way all or nothing cloud
  do k=2,nzp
     dsc(k,iv36)=dcvmgp(one,zero,drtr(k)-drvr(k)-1.e-5)
  enddo
  dsc(nzpp,iv36)=dsc(nzp,iv36)
  icldfl=nint(min(1.,dssum(nz,dsc(2,iv36),1)))

  if(icldfl /= 0) then
    !
    !     .75UM is a line of demarcation
    !
    !==================================================
    !  If passing liquid water mixing ratio from cloud1d change this
    do k=2,nzp
       dsc(k,iv1)=1.e4*(drtr(k)-drvr(k))*ddn0r(k)/ddzzr(k)
    enddo
    !==================================================
    dsc(nzpp,iv1)=0.

    !     iv2: tn1
    !     iv3: tn2
    do k=2,nzpp

       !  if W < 10 g/m^2 then use top two equations so tau linearly goes to zero
       !  old way
       !       SC(K,IV2)=.2633*SC(K,IV1)
       !       SC(K,IV3)=.3492*SC(K,IV1)
       !  new way from tripoli
       dsc(k,iv2)=.1833*dsc(k,iv1)
       dsc(k,iv3)=.2234*dsc(k,iv1)
       dsc(k,iv4)=max(10.,dsc(k,iv1))
       dsc(k,iv5)=10.**(.2633+1.7095*dlog(dlog10(dsc(k,iv4))))
       dsc(k,iv6)=10.**(.3492+1.6518*dlog(dlog10(dsc(k,iv4))))
    enddo
    do k=2,nzpp
       dsc(k,iv2)=dcvmgp(dsc(k,iv5),dsc(k,iv2),dsc(k,iv1)-10.)
       dsc(k,iv3)=dcvmgp(dsc(k,iv6),dsc(k,iv3),dsc(k,iv1)-10.)
    enddo
    !
    do k=2,nzpp
       call dpstable(1,1,dcosz,dsc(k,iv2),dsc(k,iv4),dsc(k,iv6))
       call dpstable(2,1,dcosz,dsc(k,iv3),dsc(k,iv5),dsc(k,iv6))
       call dpstable(3,1,dcosz,dsc(k,iv3),dsc(k,iv6),dsc(k,iv6))
       dsc(k,iv6)=min(dsc(k,iv6),.99999)
    enddo
    !
    !     landa.lt..75um
    !
    do k=2,nzpp
       dsc(k,iv30)=dsc(k,iv4)*dsc(k,iv2)
       dsc(k,irea)=dsc(k,iv30)/(dcosz+dsc(k,iv30))
       dsc(k,irea)=min(1.,max(0.,dsc(k,irea)))
       dsc(k,itr1)=1.-dsc(k,irea)
       dsc(k,iv11)=(1.-dsc(k,iv6)+2.*dsc(k,iv5)*dsc(k,iv6))  &
            /(1.-dsc(k,iv6))
       dsc(k,iv12)=sqrt(dsc(k,iv11))
       dsc(k,iv13)=(1.-dsc(k,iv6))*(1.-dsc(k,iv6)+2.*dsc(k,iv5)  &
            *dsc(k,iv6))
       dsc(k,iv13)=sqrt(dsc(k,iv13))*dsc(k,iv3)/dcosz
       dsc(k,iv30)=dexp(dsc(k,iv13))
       dsc(k,iv31)=1./dsc(k,iv30)
       dsc(k,iv14)=(dsc(k,iv12)+1.)*(dsc(k,iv12)+1.)*dsc(k,iv30)  &
            -(dsc(k,iv12)-1.)*(dsc(k,iv12)-1.)*dsc(k,iv31)
       !
       !       landa.gt..75um
       !
       dsc(k,ireb)=(dsc(k,iv11)-1.)*(dsc(k,iv30)-dsc(k,iv31))  &
            /dsc(k,iv14)
       dsc(k,itr2)=4.*dsc(k,iv12)/dsc(k,iv14)
       dsc(k,iab2)=1.-dsc(k,ireb)-dsc(k,itr2)
    enddo
  end if
  !              limit quantities and multiply by appropriate factors
  !                    for the flux computations
  do k=2,nzpp
     dsc(k,iv1)=dsc(k,ireb)-1.
     dsc(k,ireb)=dcvmgp(one,dsc(k,ireb),dsc(k,iv1))
     dsc(k,itr2)=dcvmgp(zero,dsc(k,itr2),dsc(k,iv1))
     dsc(k,iab2)=dcvmgp(zero,dsc(k,iab2),dsc(k,iv1))
     dsc(k,iv30)=.483*dsc(k,iv36)
     dsc(k,iv31)=.517*dsc(k,iv36)
     dsc(k,ireb)=dcvmgm(zero,dsc(k,ireb),dsc(k,ireb))*dsc(k,iv30)
     dsc(k,itr2)=dcvmgm(one,dsc(k,itr2),dsc(k,ireb))*dsc(k,iv30)
     dsc(k,iab2)=dcvmgm(zero,dsc(k,iab2),dsc(k,ireb))*dsc(k,iv30)
     dsc(k,itr1)=dsc(k,itr1)*dsc(k,iv31)
     dsc(k,irea)=dsc(k,irea)*dsc(k,iv31)
  enddo
  !
  !     Compute short wave fluxes
  !
  !-----------------------------------------------------
  !  old way for rayleigh scatter
  !     .07 reflectance due to molecular scattering
  !
  !          Rayleigh scattering
  !             VCT1=.219/(1.+.816*COSZ)*.517/3039.E3

  !     VCT1=3.7257E-8/(1.+.816*COSZ)
  !     DO K=1,NZP
  !       SC(K,IREFS)=PRD(K)*VCT1*(1.-SC(K,IV36))
  !     ENDDO
  !     SC(NZPP,IREFS)=SC(NZP,IREFS)
  !
  !-----------------------------------------------------
  !  tripoli new way include rayleigh scatter in an effective albedo

  !     Lacis and Hansen Least square fit
  !     for Rayleigh reflectance gives:

  !  k1 is the pressure at the first grid level above the surface
  !  his model was eta type coordinate where 1 could be below ground
  k1 = 2             ! this is the first level above the surface
  dpfct=dprd(k1)*1.e-6
  drabar=dpfct*0.219/(1+0.816*dcosz)
  drabarbar=dpfct*0.144
  dalbray=drabar+(1.-drabar)*(1.-drabarbar)*  &
          dalbedo/(1.-drabarbar*dalbedo)

  !     The above represents the effective albedo of the lower atmosphere
  !     due to Rayleigh scattering

  !-----------------------------------------------------
  !     Compute water vapor path
  !
  !  Vapor water path for clear atmosphere?
  do k=2,nzp
     dsc(k,iv40)=drvr(k)*(dprd(k)/1.01325e6)**.86*dradc1  &
          *ddn0r(k)/ddzzr(k)
  enddo
  dsc(nzpp,iv40)=dsc(nzp,iv40)
  !
  !     REF.....STEPHENS,1977
  !
  do k=1,nzpp
     dsc(k,iv3)=0.
     dsc(k,iv4)=0.
     !  old way with rayleigh reflection added in
     !       sc(k,iv12)=sc(k,irea)+sc(k,ireb)+sc(k,irefs)
     !  new way with out rayleigh reflection (albedo is changed)
     dsc(k,iv12)=dsc(k,irea)+dsc(k,ireb) !!!! MLO - this can have FPE Issues, that's why I switched to double precision
     dsc(k,iv30)=1.-dsc(k,iv12)-dsc(k,iab2)-dsc(k,iaoz)
     dsc(k,itr1)=dsc(k,itr1)+dsc(k,itr2)
  enddo

  !  Tak switched to 8 bands
  !                 Loop through 3 "pseudo-bands"
  !     DO NBND=1,3

  !                 Loop through 8 "pseudo-bands"
  do nbnd=1,8

     !
     !     compute ab,ref,trp,trn
     !
     dtrsmt=1.
     do k=2,nzpp

        !  new tripoli way with trans of clear * cloud
        !  because cloud tr does not include water vapor
        dsc(k,iv1)=(dsc(k,itr1)+(1.-dsc(k,iv36)))  &
             *dexp(-dsfct1(nbnd)*dsc(k,iv40))

        !  old chen way with just liquid water effects
        !         sc(k,iv1)=sc(k,itr1)+(1.-sc(k,iv36))
        !    +       *exp(-sfct1(nbnd)*sc(k,iv40))

        dsc(k,iv16)=(1.-dsc(k,iv36))*(1.-dsc(k,iv1))
     enddo

     do k=nzpp,2,-1
        dsc(k,iv13)=dsc(k,iv30)-dsc(k,iv16)*dtrsmt
        dtrsmt=dtrsmt*max(0.,dsc(k,iv13))
     enddo
     do k=2,nzpp
        dsc(k,iv14)=dsc(k,iv30)-dsc(k,iv16)*dtrsmt
        dtrsmt=dtrsmt*max(0.,dsc(k,iv14))
     enddo
     !
     !     REF.....STEPHENS,1979,J.A.S. P1542
     do k=1,nzp
        dsc(k,iv21)=dsc(nzpp-k+1,iv12)
        dsc(k,iv22)=dsc(nzpp-k+1,iv13)
        dsc(k,iv23)=dsc(nzpp-k+1,iv14)
     enddo
     dsc(1,iv26)=0.
     dsc(1,iv28)=dsolar*dcosz
     do k=1,nzp
        dsc(k,iv18)=1./(1.-dsc(k,iv26)*dsc(k,iv21))
        dsc(k+1,iv26)=dsc(k,iv21)+dsc(k,iv22)*dsc(k,iv26)  &
             *dsc(k,iv23)*dsc(k,iv18)
        dsc(k+1,iv27)=dsc(k,iv21)*dsc(k,iv28)*dsc(k,iv18)
        dsc(k+1,iv28)=dsc(k,iv22)*dsc(k,iv28)*dsc(k,iv18)
     enddo

     !  old way with old rayleigh scatter
     !       sc(nzpp,iv24)=sc(nzpp,iv28)/(1.-sc(nzpp,iv26)*albedo)
     !       sc(nzpp,iv25)=albedo*sc(nzpp,iv24)

     !  new way with effective albedo from tripoli
     dsc(nzpp,iv24)=dsc(nzpp,iv28)/(1.-dsc(nzpp,iv26)*dalbray)
     dsc(nzpp,iv25)=dalbray*dsc(nzpp,iv24)
     do k=nzp,1,-1
        dsc(k,iv25)=dsc(k,iv23)*dsc(k+1,iv25)  &
             /(1.-dsc(k,iv26)*dsc(k,iv21))+dsc(k+1,iv27)
     enddo
     do k=2,nzpp
        dsc(k,iv24)=dsc(k,iv26)*dsc(k,iv25)+dsc(k,iv28)
     enddo
     !       iv3=flxu   iv4=flxd
     do k=2,nzpp
        dsc(nzpp+1-k,iv3)=dsc(nzpp+1-k,iv3)+dsfct2(nbnd)*dsc(k,iv25)
        dsc(nzpp+1-k,iv4)=dsc(nzpp+1-k,iv4)+dsfct2(nbnd)*dsc(k,iv24)
     enddo
  enddo
  !
  !     Compute shortwave radiative tendency
  !
  do k=1,nzp
     dsc(k,iv11)=dsc(k,iv3)-dsc(k,iv4)
  enddo
  do k=2,nzp
     fthr(k)=sngl(-(dsc(k,iv11)-dsc(k-1,iv11))*ddzzr(k)/(cpcgs*ddn0r(k)))
  enddo
  rshort=-sngl(dsc(1,iv11))

  return
end subroutine shradc
!«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»!






!«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»!
![MLO - Changed the dimensions of i/o arrays to attempt to achieve binary reproducibility
subroutine lwradc(nzp,rvr,rtr,co2r,dn0r,temprd,prd,dzzr,fthr,rlong)
  !  +--------------------------------------------------------------------
  !  !  Longwave radiation parameterization based on Rodgers and
  !  !  Stephens and discussed in Chen and Cotton (1983).  First written
  !  !  by Chen, later modified by Bjorn.  All variables are in cgs so
  !  !  as to confuse people. The implicit none statement forces all
  !  !  variables to be declared!
  !  !
  !  !  Upward and downward fluxes will be calculated at w points while
  !  !  the heating rates will be at thermo points. The program takes as
  !  !  input the first five arrays passed in:
  !  ! RVR ......... vapor mixing ratio [kg/kg]
  !  ! RTR ......... total water mixing ratio [kg/kg]
  !  ! CO2R ........ CO2 mixing ratio [kg/kg]
  !  ! DN0R ........ density in (cgs units)
  !  ! TEMPRD ...... temperature in K (TEMPRD(1) is surface temperature)
  !  ! PRD ......... pressure at thermo points in (cgs units)
  !  ! DZZR ........ inverse distance between w points (cgs units)
  !  !
  !  !  and uses the following  2 entries as output or scratch space.
  !  !
  !  ! FTHRL ....... longwave heating rate at w point (K/s)
  !  ! RLONG ....... net flux into ground (or lowest level)
  !  !
  !  ! These 15 variables used to be scratch variables passed as arguments, now they became internal. 
  !  ! BB1 ......... source function for water vapor
  !  ! BB2 ......... source function for CO2
  !  ! VPR ......... water vapor path
  !  ! DMR ......... water vapor path * vapor pressure for dimer correction
  !  ! CO2 ......... CO2 path
  !  ! CLD ......... liquid water path
  !  ! EM1 & 2 ..... scratch arrays for emissivities
  !  ! SCR1 2 & 3 .. general scratch arrays
  !  ! FU1,FU2 ..... upwelling fluxes (1-vapor) (2-CO2)
  !  ! FD1,FD2 ..... downwelling fluxes (1-vapor) (2-CO2)
  !  +--------------------------------------------------------------------
  ! The subroutine uses CGS but constants are based on the global ones.
  use rconstants, only : g, cp, stefan, ep , volmoll, mmco2i
  implicit none
!----- List of arguments --------------------------------------------------------------------------------!
  integer, intent(in)                     :: nzp
  real,    intent(in)  , dimension(nzp)   :: rvr,rtr,co2r,dn0r,prd,dzzr,temprd
  real,    intent(out) , dimension(nzp)   :: fthr
  real,    intent(out)                    :: rlong

!----- List of constants --------------------------------------------------------------------------------!
  real(kind=8),    parameter                      :: gcgs=dble(g)*100.
  real(kind=8),    parameter                      :: stefancgs=1000.*dble(stefan)
  real(kind=8),    parameter                      :: prefcgs=dble(1.01325e6)
  real(kind=8),    parameter                      :: cpcgs=dble(cp)*10000.
  real(kind=8),    parameter   , dimension(4)     :: ad=(/   8.857,    -332.8,    14607.,  -261900.           /)
  real(kind=8),    parameter   , dimension(4)     :: au=(/   9.329,    -446.4,      824.,   259700.           /)
  real(kind=8),    parameter   , dimension(5)     :: bd=(/   .6558,    .12175, 1.4976e-2, 1.4981e-3,   .49e-4 /)
  real(kind=8),    parameter   , dimension(5)     :: bu=(/   .5983,    .15068, 3.4041e-2, 6.5535e-3,  4.887e-4/)
  real(kind=8),    parameter   , dimension(5)     :: ed=(/   .2837,    -.1231,    -.1057,    -.0199,  -1.16e-3/)
  real(kind=8),    parameter   , dimension(5)     :: eu=(/  .21699, -9.185e-2, -7.971e-2, -1.502e-2, -8.754e-4/)
  real(kind=8),    parameter                      :: fhl = 0., bndi = 1./200.
  real(kind=8),    parameter                      :: c11=160.87,  c21=-326.5,  c31=-158.22
  real(kind=8),    parameter                      :: c02=74.103,  c12=19.632,  c22=  0.821, c32=-0.11834
  real(kind=8),    parameter                      ::  b1= 7.345,   b2=142.47

!----- Local variables ----------------------------------------------------------------------------------!
  real                                    :: path_fact, pres_wght,sigmat_fact,trans,x,sx2,sx1
  integer                                 :: k,kk,kl,nz,nzpp,lcldbs,lcldtp
  logical                                 :: found
  !- The following arguments are only scratch variables in reality.
  real(kind=8),                  dimension(nzp+1) :: vpr,scr1,scr2,scr3,dmr,co2,cld
  real(kind=8),                  dimension(nzp+1) :: uf1,uf2,df1,df2,em1,em2,bb1,bb2
  real(kind=8),                  dimension(nzp)   :: drvr,drtr,dco2r,ddn0r,dprd,ddzzr,dtemprd
  !
  !     Water vapor band. The vibration rotation and continuum effects
  !     of the water vapor are considered.
  !
  !
  ! also downwelling flux at top of model is given by FHL
  !
  !
  nzpp=nzp+1
  nz=nzp-1
  ![MLO - Transferring information from single precision to double for internal use
  do k=1,nzp
    drvr(k)    = dble(rvr(k)   )
    drtr(k)    = dble(rtr(k)   )
    dco2r(k)   = dble(co2r(k)  ) * dble (volmoll) * dble(mmco2i) ! CO2 is now in cm3_CO2/g_air
    ddn0r(k)   = dble(dn0r(k)  )
    dprd(k)    = dble(prd(k)   )
    ddzzr(k)   = dble(dzzr(k)  )
    dtemprd(k) = dble(temprd(k))
  end do
  !
  ! calculation of optical path lengths
  !
  do k=2,nz
     pres_wght=(dprd(k)/prefcgs)**.86
     path_fact=pres_wght*ddn0r(k)/ddzzr(k)
     vpr(k)= drvr(k)*path_fact
     cld(k)=(drtr(k)-drvr(k))*path_fact
     dmr(k)=vpr(k)*dprd(k)*drvr(k)/(dble(ep)*prefcgs)
     co2(k)=dco2r(k)*path_fact
  enddo
  vpr(nzp)=vpr(nz)
  vpr(nzpp)=vpr(nz)
  cld(nzp)=0.
  cld(nzpp)=0.
  dmr(nzp)=dmr(nz)
  dmr(nzpp)=dmr(nz)
  co2(nzp)=co2(nz)
  co2(nzpp)=co2(nz) !rco2*dble(prd(nzp))/gcgs
  !
  ! computation of black body source functions with weightings given
  ! by sigma t factor of 0.87*(1/567-1/767)*TEMP
  !
  do k=1,nzp
     sigmat_fact=.0004001021*dtemprd(k)
     bb1(k)=(1.-sigmat_fact)*stefancgs*dtemprd(k)**4
     bb2(k)=bb1(k)*sigmat_fact/(1.-sigmat_fact)
  enddo
  bb1(nzpp)=bb1(nzp)
  bb2(nzpp)=bb2(nzp)
  uf1(1)=bb1(1)
  uf2(1)=bb2(1)
  !
  ! here the level of the lowest and highest cloud levels are computed
  ! so as to control the region for which mixed emmissivities need be
  ! computed
  !
  lcldbs=nzpp+1
  cbloop: do k=1,nzp
    if (drtr(k) > drvr(k)) then
      lcldbs=k
      exit cbloop
    end if
  end do cbloop

  lcldtp=0
  ctloop: do k=nzp,1,-1
     if (drtr(k) > drvr(k)) then
       lcldtp=k
       exit ctloop
     end if
  end do ctloop
  !
  ! ----------------------------- Upward computations
  !
  em1(1)=0.
  em2(1)=0.
  do k=2,nzp
     !
     ! Sum the vapor, dimer corrected vapor and co2 path lengths for use
     ! in the emissivity polynomial fits
     !
     scr1(1)=0.
     scr2(1)=0.
     scr3(1)=0.
     do kk=2,k
        scr1(kk)=scr1(kk-1)+vpr(k-kk+2)
        scr2(kk)=scr2(kk-1)+dmr(k-kk+2)
        scr3(kk)=scr3(kk-1)+co2(k-kk+2)
     enddo
     !
     ! Find level of path length of 1E-3 gm/cm^2 and compute upward
     ! emissivities for water vapor.  store in array indexed IUE1
     !
     kl=k+1
     do kk=k,2,-1
       if (scr1(kk) > 1.e-3) kl=kk
     end do

     do kk=2,kl-1
        x=dsqrt(scr1(kk))
        em1(kk)=x*(au(1)+x*(au(2)+x*(au(3)+x*au(4))))
     end do
     do kk=kl,k
        x=dlog(scr1(kk))
        em1(kk)=bu(1)+x*(bu(2)+x*(bu(3)+x*(bu(4)+x*bu(5))))
     end do
     !
     ! Correct vapor emissivities for dimer path length
     !
     found=.false.
     do kk=k,2,-1
        if(scr2(kk) > 1.e-3) then
           kl=kk
           found=.true.
        end if
     end do
     if (found) then
        do kk=kl,k
           x=dlog(min(dble(1.),scr2(kk)))
           em1(kk)=em1(kk)+eu(1)+x*(eu(2)+x*(eu(3)+x*(eu(4)+x*eu(5))))
           em1(kk)=min(dble(1.),em1(kk))
        end do
     end if
     !
     ! Compute upward emissivities for CO2, storing in IUE2, again finding
     ! the level of the critical path length
     !
     found=.false.
     kl=k+1
     do kk=k,2,-1
        if(scr3(kk) > 1.e-2) kl=kk
     end do

     do kk=2,kl-1
        x=dsqrt(scr3(kk))
        em2(kk)=1.-(x*(c11+x*(c21+x*c31)))*bndi
     enddo
     do kk=kl,k
        x=dlog(scr3(kk))
        em2(kk)=1.-(c02+x*(c12+x*(c22+x*c32)))*bndi
     enddo
     !
     ! Calculate the CO2-H2O overlap emissivity using the transmittance of
     ! water vapor given by the exponential form
     !
     do kk=2,k
        trans=dexp(-b1*scr1(kk)/sqrt(1.+b2*scr1(kk)))
        em2(kk)=1.-em2(kk)*trans
     enddo
     !
     ! if at a level greater than that of the lowest cloud compute upward
     ! emissivity for clouds and mixed emissivities as per Goody.
     !
     if(k >= lcldbs)then
        scr1(1)=0.
        do kk=2,k
           scr1(kk)=scr1(kk-1)+cld(k-kk+2)
           scr2(kk)=1.-exp(-.13e4*scr1(kk))
        enddo
        do kk=2,k
           em1(kk)=1.-(1.-em1(kk))*(1.-scr2(kk))
           em2(kk)=1.-(1.-em2(kk))*(1.-scr2(kk))
        enddo
     endif
     !
     ! compute terms in (RTE), yielding net upward fluxes
     !
     sx1=0.
     sx2=0.
     do kk=2,k
        sx1 = sx1 + bb1(k-kk+2)*(em1(kk)-em1(kk-1))
        sx2 = sx2 + bb2(k-kk+2)*(em2(kk)-em2(kk-1))
     enddo
     !
     uf1(k)=uf1(1)*(1.-em1(k))+sx1
     uf2(k)=uf2(1)*(1.-em2(k))+sx2
  enddo
  !
  ! ----------------------------- Downward computations
  !
  em1(1)=0.
  em2(1)=0.
  do k=1,nzp
     !
     ! Sum the vapor, dimer corrected vapor and co2 path lengths for use
     ! in the emissivity polynomial fits
     !
     scr1(1)=0.
     scr2(1)=0.
     scr3(1)=0.
     do kk=2,k+1
        scr1(kk)=scr1(kk-1)+vpr(kk-k+nzp)
        scr2(kk)=scr2(kk-1)+dmr(kk-k+nzp)
        scr3(kk)=scr3(kk-1)+co2(kk-k+nzp)
     enddo
     !
     ! Find level of path length of 1E-3 gm/cm^2 and compute upward
     ! emissivities for water vapor.  store in array indexed IDE1
     !
     kl=k+2
     do kk=k+1,2,-1
        if(scr1(kk) > 1.e-3) kl=kk
     end do


     do kk=2,kl-1
        x=dsqrt(scr1(kk))
        em1(kk)=x*(ad(1)+x*(ad(2)+x*(ad(3)+x*ad(4))))
     enddo
     do kk=kl,k+1
        x=dlog(scr1(kk))
        em1(kk)=bd(1)+x*(bd(2)+x*(bd(3)+x*(bd(4)+x*bd(5))))
     enddo
     !
     ! Correct vapor emissivities for dimer path length
     !
     found=.false.
     do kk=k+1,2,-1
        if(scr2(kk) > 1e-3) then
          kl=kk
          found=.true.
        end if
     end do
     if (found) then
        do kk=kl,k+1
           x=dlog(min(dble(1.),scr2(kk)))
           em1(kk)=em1(kk)+ed(1)+x*(ed(2)+x*(ed(3)+x*(ed(4)+x*ed(5))))
           em1(kk)=min(dble(1.),em1(kk))
        enddo
     end if
     !
     ! Compute upward emissivities for CO2, storing in IUE2, again finding
     ! the level of the critical path length
     !
     kl=kk+2
     do kk=k+1,2,-1
        if(scr3(kk) > 1.e-2) kl=kk
     end do

     do kk=2,kl-1
        x=dsqrt(scr3(kk))
        em2(kk)=1.-(x*(c11+x*(c21+x*c31)))*bndi
     enddo
     do kk=kl,k+1
        x=dlog(scr3(kk))
        em2(kk)=1.-(c02+x*(c12+x*(c22+x*c32)))*bndi
     enddo
     !
     ! Calculate the CO2-H2O overlap emissivity using the transmittance of
     ! water vapor given by the exponential form
     !
     do kk=2,k+1
        trans=dexp(-b1*scr1(kk)/sqrt(1.+b2*scr1(kk)))
        em2(kk)=1.-em2(kk)*trans
     enddo
     !
     ! if at a level less than that of the highest cloud compute downward
     ! emissivity for clouds and mixed emissivity as per Goody.
     !
     if(nzp+2-k <= lcldtp)then
        scr1(1)=0.
        do kk=2,k+1
           scr1(kk)=scr1(kk-1)+cld(kk-k+nzp)
           scr2(kk)=1.-dexp(-.158e4*scr1(kk))
        enddo
        do kk=2,k+1
           em1(kk)=1.-(1.-em1(kk))*(1.-scr2(kk))
           em2(kk)=1.-(1.-em2(kk))*(1.-scr2(kk))
        enddo
     endif
     !
     ! compute terms in radiative transfer equation
     !
     sx1=0.
     sx2=0.
     do kk=2,k+1
        sx1 = sx1 + bb1(kk-k+nzp)*(em1(kk)-em1(kk-1))
        sx2 = sx2 + bb2(kk-k+nzp)*(em2(kk)-em2(kk-1))
     enddo
     !
     ! Compute downward fluxes for water vapor (iuf1) and co2 (iuf2)
     !
     df1(nzpp-k)=fhl*(1.-em1(k+1))+sx1
     df2(nzpp-k)=fhl*(1.-em2(k+1))+sx2
  enddo
  !
  ! ----------------------------- Net flux & heating rates
  !
  scr1(1)=uf1(1)-df1(1)+uf2(1)-df2(1)
  do k=2,nzp
     scr1(k)=uf1(k)-df1(k)+uf2(k)-df2(k)
     fthr(k)=sngl(-(scr1(k)-scr1(k-1))*ddzzr(k)/(cpcgs*ddn0r(k)))
  enddo
  rlong=sngl(df1(1)+df2(1))
  !
  return
end subroutine lwradc
!«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»!






!«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»!
subroutine shradp(nzp,rvr,dn0r,dzr,sc,pird,cosz,albedo  &
     ,solar,fthr,rshort)
  !+----------------------------------------------------------------------!
  !     Shortwave radiation parameterization described in Mahrer and      !
  !     Pielke(1977).                                                     !
  !                                                                       !
  !       Arguments:                                                      !
  !       ----------                                                      !
  !                                                                       !
  !       Input:  NZP    - number of vertical levels                      !
  !               RVR    - water vapor at each level                      !
  !               DN0R   - air density at each level                      !
  !               DZR    - inverse of delta z = 1./(Z(K)-Z(K-1))          !
  !                        where Z is at the vapor levels                 !
  !               PIRD   - Exner function (p/p00)**(R/Cp) at each level   !
  !               SC     - scratch array at least 2*NZP long              !
  !               COSZ   - cosine of the zenith angle                     !
  !               ALBEDO - albedo of the ground surface                   !
  !               SOLAR  - the solar constant                             !
  !                                                                       !
  !       Output: FTHR   - radiation tendency on potential temperture.    !
  !               RSHORT - downward shortwave flux on a flat surface at   !
  !                        the ground                                     !
  !                                                                       !
  !+----------------------------------------------------------------------!
  use rconstants, only : cp
  implicit none
  !----- Arguments: ------------------------------------------------------!
  integer, intent(in)                     :: nzp
  real   , intent(in)                     :: cosz,albedo,solar
  real   , intent(in)  , dimension(nzp)   :: rvr,dn0r,pird,dzr
  real   , intent(out)                    :: rshort
  real   , intent(out) , dimension(nzp)   :: fthr
  real   , intent(out) , dimension(nzp,2) :: sc
  !----- List of constants -----------------------------------------------!
  integer, parameter                      :: iv1=1,iv2=2
  real, parameter                         :: cpcgs=10000.*cp

  integer                                 :: nz,k
  real                                    :: raysct,rdcon1,vabs
  !----- List of functions -----------------------------------------------!
  real, external                          :: ssum
  !-----------------------------------------------------------------------!
  nz=nzp-1

  !     Rayleigh scattering (numerator in SQRT should be
  !        SQRT((.000949*P+.051)/COSZ), but is ignored. (P in mb)

  raysct=1.021-.0824*sqrt(1./cosz)
  !     Vapor path length

  do k=1,nz
     sc(k,iv1)=(rvr(k)*dn0r(k)+rvr(k+1)*dn0r(k+1))*.5/dzr(k)
     sc(k,iv1)=max(sc(k,iv1),1e-10)
  enddo
  do k=1,nz
     sc(k,iv2)=ssum(nzp-k,sc(k,iv1),1)+1e-10
  enddo
  sc(nzp,iv2)=0.

  !     Shortwave heating by vapor absorbtion

  rdcon1=.0231*solar/cpcgs
  do k=2,nz
     fthr(k)=(rdcon1*(sc(k,iv2)/cosz)**(-.7)*cosz*rvr(k))/pird(k)
  enddo
  vabs=.077*(sc(1,iv2)/cosz)**.3

  !     Shortwave on a flat surface

  rshort=max(solar*cosz*(1.-albedo)*(raysct-vabs),0.)

  return
end subroutine shradp
!«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»!






!«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»!
subroutine lwradp(nzp,temprd,rvr,co2r,dn0r,dzzr,pird,sc,fthr,rlong)
  !--------------------------------------------------------------------------------------------------------------------!
  !     Longwave radiation parameterization described in Mahrer and                                                    !
  !     Pielke (1977).  Does not include any cloud effects.                                                            !
  !                                                                                                                    !
  !       Arguments:                                                                                                   !
  !       ----------                                                                                                   !
  !                                                                                                                    !
  !       Input:  NZP    - number of vertical levels                                                                   !
  !               TEMPRD - temperature in Kelvin at each level                                                         !
  !               RVR    - water vapor at each level                                                                   !
  !               CO2R   - CO2 mixing ratio at each level                                                              !
  !               DN0R   - air density at each level                                                                   !
  !               DZZR   - inverse of delta z = 1./(ZZ(K)-ZZ(K-1))                                                     !
  !                        where ZZ is staggered midway between T levels                                               !
  !               PIRD   - Exner function (p/p00)**(R/Cp) at each level                                                !
  !               SC     - scratch array at least 20*NZP long                                                          !
  !                                                                                                                    !
  !       Output: FTHR   - radiation tendency on potential temperture.                                                 !
  !               RLONG  - downward longwave flux at the ground                                                        !
  !                                                                                                                    !
  !--------------------------------------------------------------------------------------------------------------------!
  use rconstants, only : g,cp,stefan,p00,cpor
  implicit none
  integer , intent(in)            :: nzp
  real    , intent(in)   , dimension(nzp)    :: rvr,co2r,dn0r,temprd,dzzr, pird
  real    , intent(out)                      :: rlong
  real    , intent(out)  , dimension(nzp)    :: fthr
  real    , intent(out)  , dimension(nzp,18) :: sc
  !----- List of constants --------------------------------------------------------------------------------------------!
  real    , parameter                         :: p00cgs    =10.    * p00
  real    , parameter                         :: gcgs      =100.   * g
  real    , parameter                         :: stefancgs =1000.  * stefan
  real    , parameter                         :: cpcgs     =10000. * cp
  integer , parameter                         ::   iv1=  1,  iv2=  2,  iv3=  3,  iv4=  4,  iv5=  5,  iv6=  6 &
                                                ,  iv7=  7,  iv8=  8,  iv9=  9, iv10= 10, iv11= 11, iv12= 12 &
                                                , iv13= 13, iv14= 14, iv15= 15, iv16= 16, iv17= 17, iv18= 18
  !----- Local variables ----------------------------------------------------------------------------------------------!
  integer                                     :: nz,nz1,k
  !--------------------------------------------------------------------------------------------------------------------!
  nz=nzp-1
  nz1=nz-1

  !                   Compute upward and downward vapor path
  do k=2,nz
     sc(k,iv1)=rvr(k)*dn0r(k)/dzzr(k)
     sc(k,iv1)=max(sc(k,iv1),1e-10)
  enddo
  sc(nzp,iv1)=sc(nz,iv1)
  !
  sc(1,iv2)=0.
  do k=2,nz
     sc(k,iv2)=sc(k-1,iv2)+sc(k,iv1)
  enddo
  sc(nz,iv3)=sc(nzp,iv1)
  do k=nz1,1,-1
     sc(k,iv3)=sc(k+1,iv3)+sc(k+1,iv1)
  enddo

  !                          water vapor emissivity calculation
  do k=1,nz
     sc(k,iv4)=log10(sc(k,iv2)+1e-30)
     sc(k,iv5)=log10(sc(k,iv3)+1e-30)
  enddo

  do k=1,nz
     if(sc(k,iv4) <= -4.0                        )  sc(k,iv6)=.1129*log10(1.+12.63*sc(k,iv2))
     if(sc(k,iv4) <= -3.0 .and. sc(k,iv4) >  -4.0)  sc(k,iv6)=.104*sc(k,iv4)+.440
     if(sc(k,iv4) <= -1.5 .and. sc(k,iv4) >  -3.0)  sc(k,iv6)=.121*sc(k,iv4)+.491
     if(sc(k,iv4) <= -1.0 .and. sc(k,iv4) >  -1.5)  sc(k,iv6)=.146*sc(k,iv4)+.527
     if(sc(k,iv4) <=  0.0 .and. sc(k,iv4) >  -1.0)  sc(k,iv6)=.161*sc(k,iv4)+.542
     if(                        sc(k,iv4) >   0.0)  sc(k,iv6)=.136*sc(k,iv4)+.542

     if(sc(k,iv5) <= -4.                         )  sc(k,iv7)=.1129*log10(1.+12.63*sc(k,iv3))
     if(sc(k,iv5) <= -3.0 .and. sc(k,iv5) >  -4.0)  sc(k,iv7)=.104*sc(k,iv5)+.440
     if(sc(k,iv5) <= -1.5 .and. sc(k,iv5) >  -3.0)  sc(k,iv7)=.121*sc(k,iv5)+.491
     if(sc(k,iv5) <= -1.0 .and. sc(k,iv5) >  -1.5)  sc(k,iv7)=.146*sc(k,iv5)+.527
     if(sc(k,iv5) <=  0.0 .and. sc(k,iv5) >  -1.0)  sc(k,iv7)=.161*sc(k,iv5)+.542
     if(                        sc(k,iv5) >   0.0)  sc(k,iv7)=.136*sc(k,iv5)+.542
  enddo

  !                           CO2 path lengths and emissivities
  do k=2,nz
     sc(k,iv11)=co2r(k) * gcgs * dn0r(k) / dzzr(k)
  enddo
  sc(nzp,iv11)=co2r(nzp) * pird(nzp)**(cpor)*p00cgs

  sc(1,iv12)=0.
  do k=2,nz
     sc(k,iv12)=sc(k-1,iv12)+sc(k,iv11)
  enddo
  sc(nz,iv13)=sc(nzp,iv11)
  do k=nz1,1,-1
     sc(k,iv13)=sc(k+1,iv13)+sc(k+1,iv11)
  enddo

  do k=1,nz
     sc(k,iv8)=.185*(1.-exp(-.3919*sc(k,iv12)**.4))
     sc(k,iv9)=.185*(1.-exp(-.3919*sc(k,iv13)**.4))
  enddo

  !                        Add CO2 and H2O emissivities, find SIG(T**4)
  do k=1,nzp
     sc(k,iv14)=sc(k,iv8)+sc(k,iv6)
     sc(k,iv15)=sc(k,iv9)+sc(k,iv7)
     sc(k,iv16)=stefancgs*temprd(k)**4
  enddo
  sc(nzp,iv16)=sc(nzp,iv16)*sc(nz,iv15)

  !                       Calculate upward and downward divergences
  do k=2,nz
     sc(k,iv17)=(sc(k,iv16)-sc(1,iv16))*(sc(k,iv14)-sc(k-1,iv14))
     sc(k,iv18)=(sc(nzp,iv16)-sc(k,iv16))*(sc(k,iv15)-sc(k-1,iv15))
  enddo
  sc(1,iv18)=sc(nzp,iv16)*(1.-sc(1,iv15))+sc(2,iv16)*sc(2,iv15)

  do k=2,nz
     fthr(k)=-(sc(k,iv17)+sc(k,iv18))*dzzr(k)/(cpcgs*dn0r(k)*pird(k))
  enddo
  rlong=sc(1,iv18)

  !------------------------------------------------------------------
  !      print 6667,(k,temprd(k),sc(k,iv6),sc(k,iv7),sc(k,iv8),sc(k,iv9)
  !     +  ,sc(k,iv13),fthr(k)*24.*3600.,k=nzp,1,-1)
  ! 6667 format(' longwave-t,emissur,emissdr,emissuc,emissdc,pathd,d/day',/,(i3,7e10.3))
  !      print 6668,(k,temprd(k),sc(k,iv2),sc(k,iv3),sc(k,iv12),sc(k,iv13)
  !     +  ,k=nzp,1,-1)
  ! 6668 format(' longwave-t,up vap,dn vap,up co2, dn co2',/,(i3,5e10.3))
  !      print*,'  longwave down ',-sc(1,iv18)*1e-3

  return
end subroutine lwradp
!«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»!
