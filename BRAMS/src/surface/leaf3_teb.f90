!############################# Change Log ##################################
! 5.0.2
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003, 2006 - All Rights Reserved
!  Brazilian Regional Atmospheric Modeling System - BRAMS
!###########################################################################

!Subroutine to link leaf3 to TEB for landclass = urban type 1 and 
!urban type 2
!Adapted by Edmilson Freitas 
!DATE : Jul 7th 2006
!Last modification:

! TEB
subroutine LEAF3_TEB_INTERFACE(ISTP,ZTSTEPFRC,ZTSTEP,COSZ,ZZREF,  &
                               ZRAT, SOLAR,                       &
                               ZPA,ZTA,ZUA,ZVA,RV,                &
                               ZPCP_IN,fuso,                      &
			       ZT_CANYON,ZR_CANYON,               &
			       ZTS_ROOF,ZTS_ROAD,ZTS_WALL,        &
			       ZTI_ROAD,ZTI_BLD,                  &
			       ZWS_ROOF,ZWS_ROAD,                 &
 			       ZT_ROOF,ZT_ROAD,ZT_WALL,           &
                               ZH_TOWN,ZLE_TOWN,ZEMIS_TOWN,       &
                               ZSFU_TOWN,ZSFV_TOWN,ZTS_TOWN,      &
			       ZALB_TOWN,IG_URBAN,                &
			       ZH_TRAFFIC,ZH_INDUSTRY,            &
			       ZLE_TRAFFIC,ZLE_INDUSTRY,          &
			       T2M,R2M,                           &
			       time,itime1,dpdz,dens)

  use teb_vars_const
  use mem_emiss, only:EFSAT,EFSUN,WEEKDAYIN     
  use therm_lib, only: rslif
  implicit none
  
  ! Declarations of variables
  !INPUT / OUTPUT VARIABLES
  integer, parameter     :: INTEB = 3    ! TEB Vertical dimension
  integer                :: IG_URBAN      ! Flag for urban(1) and suburban(2)
  integer                :: itime1
  real(kind=8)                   :: time
  real                   :: pfat,dpdz
  !
  real                    :: ZEMIS_TOWN! town equivalent emissivity
  real                    :: ZT_CANYON! canyon air temperature
  real                    :: ZR_CANYON! canyon air vapor mixing ratio
  real                    :: ZTS_ROOF ! roof surface temperature
  real                    :: ZTS_ROAD ! road surface temperature
  real                    :: ZTS_WALL ! wall surface temperature
  real                    :: ZTI_ROAD ! road deep temperature
  real,dimension(INTEB)   :: ZT_ROOF  ! roof layers temperatures
  real,dimension(INTEB)   :: ZT_ROAD  ! road layers temperatures
  real,dimension(INTEB)   :: ZT_WALL  ! wall layers temperatures
  real                    :: ZWS_ROOF ! roof water reservoir
  real                    :: ZWS_ROAD ! road water reservoir
  real                    :: T2M      ! Extrapolated 2 m temperature
  real                    :: R2M      ! Extrapolated 2 m specif humidity
  real                    :: fuso     ! local time correction
  
  real                    :: ZTI_BLD      ! INTERNAL BLDING TEMP (K)
  
  
  real	:: ZTSTEP 		! Time step of integration (sec.)
  real	:: COSZ			! Cos(Zenith Angle)
  
  real    :: ZZREF
  !                                    ZZREF = reference height of the first
  !                                            atmospheric level
  
  real	:: ZRG,ZRAT,ZUA,ZVA,ZPS,ZQA,ZTA,ZRR,RV
  real 	:: SOLAR,ZPRECIP
  real    :: ZTSTEPFRC		!time step of the input forcing
  
  !                                    ZRG = incoming solar radiation
  !                                    ZRAT = atmospheric infrared radiation
  !                                    ZRR = rain rate
  !                                    ZTA = atmospheric temperature at level za
  !                                    ZUA, ZVA = wind speeds in the x and y
  !                                               directions at level za
  !                                    ZPS = surface pressure
  !                                    ZQA = atmospheric specific humidity
  !                                          at level za
  !                                    RV  = atmospheric mixing ratio
  !                                          at level za
  !                                    ZPRECIP = precipitation rate
  !
  !       0.9.D   TEB Diagnostics:
  !	        ----------------
  real        :: ZRN_ROOF     ! net radiation over roof
  real        :: ZH_ROOF      ! sensible heat flux over roof
  real        :: ZLE_ROOF     ! latent heat flux over roof
  real        :: ZGFLUX_ROOF  ! flux through the roof
  real        :: ZRUNOFF_ROOF ! runoff over the ground
  real        :: ZRN_ROAD     ! net radiation over road
  real        :: ZH_ROAD      ! sensible heat flux over road
  real        :: ZLE_ROAD     ! latent heat flux over road
  real        :: ZGFLUX_ROAD  ! flux through the road
  real        :: ZRUNOFF_ROAD ! runoff over the ground
  real        :: ZRN_WALL     ! net radiation over wall
  real        :: ZH_WALL      ! sensible heat flux over wall
  real        :: ZLE_WALL     ! latent heat flux over wall
  real        :: ZGFLUX_WALL  ! flux through the wall
  !
  
  !
  !only OUTPUT VARIABLES
  real	::	ZTS_TOWN	!TOWN SFC TEMP
  real	::	ZALB_TOWN	!TOWN EQV. ALBEDO
  real	::	ZU_CANYON	!CANYON HOR. WIND
  real	::	ZRN_TOWN	!TOWN SCALE NET RAD
  real	::	ZH_TOWN		!TOWN AVE. SENS. HEAT FLUX
  real	::	ZLE_TOWN	!TOWN AVE. LAT. HEAT FLUX
  real	::	ZSFU_TOWN	!TOWN SCALE EDDY U-MOM FLUX
  real	::	ZSFV_TOWN	!TOWN SCALE EDDY V-MOM FLUX
  real	::	ZCH_TOWN	!TOWN AVERAGED HEAT TRANSFER
  real	::	ZGFLUX_TOWN	!TOWN SCALE GROUND HEAT STORAGE
  real	::	ZRUNOFF_TOWN	!TOWN SCALE RUNOFF
  
  !  MODEL FORCING VARIABLES
  real	::      ZPCP_IN         ! Precip. Rate (kg/m^2/s)
  
  real	:: 	ZSVF_ROAD	!ROAD SKY VIEW FACTOR
  real	:: 	ZSVF_WALL	!WALL SKY VIEW FACTOR
  real	::      ZCAN_HW_RATIO   !CANYON HEIGHT TO WIDTH RATIO
  
  real    :: ZDIR_SW_RAD  ! incoming direct solar radiation
  !                       ! on an horizontal surface
  real    :: ZSCA_SW_RAD  ! scattered incoming solar rad.
  
  integer I,JDT,JINDT
  integer :: ISTP
  
  real ::       ZZ0_TOWN	    &
       ,	ZBLD		    &
       ,	ZBLD_HEIGHT	    &
       ,	ZBLD_HL_RATIO	    &
       ,	ZALB_ROOF	    &
       ,	ZEMIS_ROOF	    &
       ,	ZALB_ROAD	    &
       ,	ZEMIS_ROAD	    &
       ,	ZALB_WALL	    &
       ,	ZEMIS_WALL	    &
       ,	ZH_TRAFFIC 	    &
       ,	ZH_INDUSTRY 	    &
       ,	ZLE_TRAFFIC 	    &
       ,	ZLE_INDUSTRY 	    &
       ,       ax1,ax2,bx1,bx2     &
       ,       timeq1,timeq2,dens

  real(kind=8) :: tmp_teb,tign
  integer      :: idays
  character(len=3)cday
  
  real :: ZDIRCOSZW, ZTANZEN
  !                  ZDIRCOSZW = Cosinus of the angle between the 
  !                  normal to the surface and the vertical
  !                  ZTANZEN   = tangent of solar zenith angle
  !
  real   :: ZRHOA,  ZEXNA, ZEXNS, ZVMOD,  ZPA, ZTVI, ZAZENIT,ZEXN2,ZP2
  !                  ZRHOA   = air density
  !                  ZEXNA   = Exner function
  !                  ZEXNS   = Exner function
  !                  ZVMOD   = modulus of the wind parallel to the orography
  !                  ZPA     = atmospheric level pressure
  !                  ZTVI    = virtual temperature 
  !                  ZAZENIT = solar zenith angle
  !
  
  ZRG = SOLAR
  
  ZDIRCOSZW = 1.0
  
  ZQA=(1./((1000./rv)+1.))*1000.
  
  
  ZZ0_TOWN        =  Z0_TOWN(IG_URBAN)    
  ZBLD            =  BLD (IG_URBAN) 
  ZBLD_HEIGHT     =  BLD_HEIGHT (IG_URBAN)    
  ZBLD_HL_RATIO   =  BLD_HL_RATIO(IG_URBAN)
  ZALB_ROOF       =  AROOF(IG_URBAN)
  ZEMIS_ROOF      =  EROOF(IG_URBAN)
  ZALB_ROAD       =  AROAD(IG_URBAN)
  ZEMIS_ROAD      =  EROAD(IG_URBAN)
  ZALB_WALL       =  AWALL(IG_URBAN)
  ZEMIS_WALL      =  EWALL(IG_URBAN)
  
  ZH_TRAFFIC      =  HTRAF(IG_URBAN)
  ZH_INDUSTRY     =  HINDU(IG_URBAN)	
  ZLE_TRAFFIC     =  PLETRAF(IG_URBAN)	 
  ZLE_INDUSTRY    =  PLEINDU(IG_URBAN)
  
  !
  !*      1.     Calculate the canyon shape
  !              --------------------------
  !
  ZCAN_HW_RATIO = ZBLD_HL_RATIO*ZBLD /  &
       (1.0 - ZBLD)
  !write(*,*)'ZCAN_HW_RATIO=',ZCAN_HW_RATIO
  !write(*,*)'ZBLD_HL_RATIO=',ZBLD_HL_RATIO
  !pause				
  !!*      2.     Calculate sky view factors
  !              --------------------------
  !
  ZSVF_ROAD  = sqrt(ZCAN_HW_RATIO*ZCAN_HW_RATIO + 1.0) -  &
       ZCAN_HW_RATIO
  ZSVF_WALL = 0.5*(1.0 - ZSVF_ROAD)/ZCAN_HW_RATIO
  
  !*      3.     Calculate town albedo and emissivity
  !              ------------------------------------
  !
  ZALB_TOWN  = (1.0 - ZBLD)*ZALB_ROAD + ZBLD*ZALB_ROOF
  
  ZEMIS_TOWN = (1.0 - ZBLD)*ZEMIS_ROAD*       ZSVF_ROAD   &
       + (1.0 - ZBLD)*ZEMIS_WALL*(1.0 - ZSVF_ROAD)  &
       +        ZBLD *ZEMIS_ROOF
  
  ZRR      = ZPCP_IN
  
  if (ZRG < 0.0)ZRG = 0.0
  
  ! Initialize forcing values:
  
  ZAZENIT = acos(COSZ)
  ZTANZEN = tan(ZAZENIT)
  if(ZRG > 0.0 .and. ZTANZEN < 0.0)then
     ZTANZEN = sqrt(9999.0)
  endif
  
  ! Perform some conversions and final calculations using forcing variables:
  !                
  ZVMOD   = sqrt( ZUA*ZUA + ZVA*ZVA )
  !
  ! Make sure wind magnitude is above a minimum threshold (original =1.0, 
  ! for test purpose = 0.5) :
  !
  ZVMOD   = max(0.5,ZVMOD)
  ZTVI    = ZTA * ( 1.+((XRV/XRD)-1.)*ZQA )
  ZRHOA   = ZPA / XRD / ZTVI
  !            ZRHOA   = dens
  ZEXNA   = (ZPA/100000.)**(XRD/XCPD)
  
  ZPS     = DPDZ*(ZBLD_HEIGHT-ZZREF)+ZPA
  ZEXNS   = (ZPS/100000.)**(XRD/XCPD)
  
  ! Calculate forcing needed by TEB from existing forcing
  ! for now assume all incoming radiation is direct
  !
  ZDIR_SW_RAD = ZRG
  ZSCA_SW_RAD = 0.0
  !
  
  tmp_teb=time+(itime1/100+mod(itime1,100)/60.)*3600
  idays = int((tmp_teb/3600.)/24.)  !number of days of simulation
  tign = dble(idays)*24.*3600.
  
  pfat=1.
  
  call EMFACTOR(WEEKDAYIN,idays,cday)
  if(cday=='SAT')pfat=EFSAT
  if(cday=='SUN')pfat=EFSUN
  
  !Fonte urbana (castanho, 1999 - tese de mestrado)
  
  bx1=RUSHH1-fuso+DAYLIGHT
  bx2=RUSHH2-fuso+DAYLIGHT
  timeq1= real( tmp_teb - tign )/3600. - bx1
  timeq2= real( tmp_teb - tign )/3600. - bx2
    
  !fonte de calor sensivel veicular
  ax2=ZH_TRAFFIC-5.
  ax1=0.63*ax2
  ZH_TRAFFIC=( (ax1*exp(-(timeq1)**2/8.5)   +      &
       ax2*exp(-(timeq2)**2/10.6))*pfat  +  5.)
  
  !fonte de calor latente veicular
  ax2=ZLE_TRAFFIC-5.
  ax1=0.63*ax2
  ZLE_TRAFFIC=( (ax1*exp(-(timeq1)**2/8.5)   +      &
       ax2*exp(-(timeq2)**2/10.6))*pfat  +  5.)
  
  !                              -------------------------------
  ZRN_ROOF=0.
  ZH_ROOF=0.
  ZLE_ROOF=0.
  ZGFLUX_ROOF=0.
  ZRUNOFF_ROOF=0.
  ZRN_ROAD=0.
  ZH_ROAD=0.
  ZLE_ROAD=0.
  ZGFLUX_ROAD=0.
  ZRUNOFF_ROAD=0.
  ZRN_WALL=0.
  ZH_WALL=0
  ZLE_WALL=0.
  ZGFLUX_WALL=0.
  ZRN_TOWN=0.
  ZH_TOWN=0.
  ZLE_TOWN=0.
  ZGFLUX_TOWN=0.
  ZRUNOFF_TOWN=0.
  ZSFU_TOWN=0.
  ZSFV_TOWN=0.
  ZCH_TOWN=0.
  
  !*                        14.B   CALL THE MAIN SUBROUTINE OF TEB
    
  call URBAN(ZTS_TOWN, ZEMIS_TOWN, ZALB_TOWN,                   &
       ZT_CANYON, ZR_CANYON, ZU_CANYON,                         &
       ZTS_ROOF,ZTS_ROAD,ZTS_WALL,ZTI_ROAD,ZTI_BLD,             &
       ZT_ROOF, ZT_ROAD, ZT_WALL, ZWS_ROOF,ZWS_ROAD,            &
       ZPS, ZPA, ZEXNS, ZEXNA, ZTA, ZQA, ZRHOA,                 &
       ZRAT, ZDIR_SW_RAD, ZSCA_SW_RAD, ZTANZEN,                 &
       ZRR, ZZREF, ZDIRCOSZW, ZUA, ZVA, ZVMOD,                  &
       ZH_TRAFFIC, ZLE_TRAFFIC, ZH_INDUSTRY, ZLE_INDUSTRY,      &
       ZTSTEP,                                                  &
       ZZ0_TOWN,                                                &
       ZBLD,ZBLD_HEIGHT,2.*ZBLD_HL_RATIO*ZBLD,ZCAN_HW_RATIO,    &
       ZALB_ROOF, ZEMIS_ROOF,                                   &
       HC_ROOF(1:3),TC_ROOF(1:3),D_ROOF(1:3),                   &
       ZALB_ROAD, ZEMIS_ROAD, ZSVF_ROAD,                        &
       HC_ROAD(1:3),TC_ROAD(1:3),D_ROAD(1:3),                   &
       ZALB_WALL, ZEMIS_WALL, ZSVF_WALL,                        &
       HC_WALL(1:3),TC_WALL(1:3),D_WALL(1:3),                   &
       ZRN_ROOF, ZH_ROOF, ZLE_ROOF, ZGFLUX_ROOF,                &
       ZRUNOFF_ROOF,                                            &
       ZRN_ROAD, ZH_ROAD, ZLE_ROAD, ZGFLUX_ROAD,                &
       ZRUNOFF_ROAD,                                            &
       ZRN_WALL, ZH_WALL, ZLE_WALL, ZGFLUX_WALL,                &
       ZRN_TOWN, ZH_TOWN, ZLE_TOWN, ZGFLUX_TOWN,                &
       ZRUNOFF_TOWN, ZSFU_TOWN, ZSFV_TOWN, ZCH_TOWN)

  !Extrapolating temperature and specific humidity to 2 m height
  
  ZP2    = DPDZ*(2.-ZBLD_HEIGHT)+ZPS
  
  ZEXN2   = (ZP2/100000.)**(XRD/XCPD)
  
  T2M = ZT_WALL(1) * ZEXN2 / ZEXNS
  
  R2M = ZR_CANYON * (rslif(ZP2,T2M) / rslif(ZPS,ZT_WALL(1)))
  
  return
end subroutine LEAF3_TEB_INTERFACE

! TEB
!     ############################################
subroutine INI_TG_PROFILE(PTS,PTI,PTC,PD,PT)
  !   ############################################
  !
  !!****  *INI_TG_PROFILE*  
  !!
  !!    PURPOSE
  !!    -------
  !
  !     Computes the equilibrium of a temperature profile through a
  !     conductive material, given the two extreme surface temperatures.
  !         
  !!      
  !!    AUTHOR
  !!    ------
  !!
  !!	V. Masson           * Meteo-France *
  !!
  !!    MODIFICATIONS
  !!    -------------
  !!      Original 02/11/98 
  !!               17/05/2002 Edmilson Freitas (IAG-USP). Elimination of the
  !!                          the declarations that are not needed. Change from
  !!                          Module to normal subroutine.
  !----------------------------------------------------------------------------
  !
  !*       0.     DECLARATIONS
  !               ------------
  !
  !
  implicit none
  !
  !*      0.1    declarations of arguments
  !
  !INPUT VARIABLES
  !
  real  :: PTS      ! surface temperature
  real  :: PTI      ! internal temperature
  real, dimension(3)  :: PTC      ! thermal conductivity for roof layers
  real, dimension(3)  :: PD       ! depth of roof layers
  !
  !OUTPUT VARIABLES
  !
  real, dimension(3)  :: PT       ! layers temperatures
  !
  !*      0.2    declarations of local variables
  !
  real, dimension(size(PT)) :: ZA ! lower diag.
  real, dimension(size(PT)) :: ZB ! main  diag.
  real, dimension(size(PT)) :: ZC ! upper diag.
  real, dimension(size(PT)) :: ZY ! r.h.s.
  real, dimension(size(PT)) :: ZX ! solution
  !
  real, dimension(0:size(PT)) :: ZMTC_O_D
  ! mean thermal conductivity over distance between 2 layers
  !
  integer :: ILAYER           ! number of roof,road or wall layers
  integer :: JLAYER           ! loop counter
  !----------------------------------------------------------------------------
  !
  !*      1.     Layer thermal properties
  !              ------------------------
  !
  ILAYER = size(PT)
  ZA(:) = 0.
  ZB(:) = 0.
  ZC(:) = 0.
  ZX(:) = 0.
  ZY(:) = 0.
  ZMTC_O_D(:) = 0.
  !
  ZMTC_O_D(0) = 2. * PTC(1) / PD (1)
  !
  do JLAYER=1,ILAYER-1
     ZMTC_O_D(JLAYER) = 2./(   PD(JLAYER  )/PTC(JLAYER  ) &
          + PD(JLAYER+1)/PTC(JLAYER+1) )
  end do
  !
  ZMTC_O_D(ILAYER) = 2. * PTC(ILAYER) &
       / PD (ILAYER)
  !
  !-----------------------------------------------------------------------------
  !
  !*      2.     Surface layer coefficients
  !              --++++++------------------
  !
  !
  ZA(1) =   0.
  
  ZB(1) =   ZMTC_O_D(0)                          &
       + ZMTC_O_D(1)
  
  ZC(1) = - ZMTC_O_D(1)
  !
  ZY(1) =   ZMTC_O_D(0) * PTS
  !
  !
  !----------------------------------------------------------------------------
  !
  !*      3.     Other layers coefficients
  !              -------------------------
  !
  do JLAYER=2,ILAYER-1
     ZA(JLAYER) = - ZMTC_O_D(JLAYER-1)
     
     ZB(JLAYER) =   ZMTC_O_D(JLAYER-1)                                &
          + ZMTC_O_D(JLAYER  )
     
     ZC(JLAYER) = - ZMTC_O_D(JLAYER  )
     !
     ZY(JLAYER) =   0.
  end do
  !
  !-----------------------------------------------------------------------------
  !
  !*      4.     Inside layer coefficients
  !              -------------------------
  !
  ZA(ILAYER) = - ZMTC_O_D(ILAYER-1)
  
  ZB(ILAYER) =   ZMTC_O_D(ILAYER-1)                          &
       + ZMTC_O_D(ILAYER  )
  
  
  ZC(ILAYER) =   0.
  !
  ZY(ILAYER) = ZMTC_O_D(ILAYER) * PTI
  !
  !
  !-----------------------------------------------------------------------------
  !
  !*      5.     Tri-diagonal system resolution
  !              ------------------------------
  !
  call TRID(ZX,ZA,ZB,ZC,ZY,ILAYER)
  !
  PT(:) = ZX(:)
  !
  !-----------------------------------------------------------------------------
  !
  return
end subroutine INI_TG_PROFILE

!-------------------------------------------------------------------------------
! TEB
subroutine TEB_INIT(ng,n1,n2,n3,np,vegt,theta,rv,pi,pp,TROOF,TROAD,TWALL, &
     TIBLD,TIROAD,TCANYON,RCANYON,TSROOF,TSROAD,TSWALL,                   &
     HT,LET,HIN,LEIN,WSROOF,WSROAD,EMISTOWN,ALBTOWN,                      &
     TSTOWN,G_URBAN)
  use rconstants, only: t00
  use teb_vars_const, only: TMINBLD,D_ROAD,TC_ROAD        &
       ,D_WALL,TC_WALL,D_ROOF,TC_ROOF &
       ,NURBTYPE,ILEAFCOD 

  implicit none
  integer :: n1,n2,n3,ng,np,i,j,k,ilf,inp
  
  real, dimension(n1,n2,n3) :: theta,rv,pi,TROOF,TROAD,TWALL,pp
  real, dimension(n2,n3,np) :: vegt,G_URBAN
  real, dimension(n2,n3)    :: TIBLD,TIROAD,TCANYON,RCANYON,TSROOF,TSROAD, &
       TSWALL, HT,LET,HIN,LEIN,WSROOF,WSROAD
  real, dimension(n2,n3) :: EMISTOWN,ALBTOWN,TSTOWN
  
  real :: cpi,hcpi,pis,pl2
  
  HT=0.
  LET=0.
  HIN=0.
  LEIN=0.
  WSROOF=0.
  WSROAD=0.
  TIBLD=0.
  TIROAD=0.
  TCANYON=0.
  RCANYON=0.
  TSROOF=0.
  TSROAD=0.
  TSWALL=0.
  TROOF=0.
  TROAD=0.
  TWALL=0.
  G_URBAN=0.
  
  cpi = 1./1004.
  hcpi = .5 * cpi
  
  do i=1,n2
     
     do j=1,n3
        
        pis = (pp(1,i,j) + pi(1,i,j) + pp(2,i,j) + pi(2,i,j)) * hcpi
        
        pl2 = (pp(2,i,j) + pi(2,i,j)) * cpi
        
        do inp=2,np !patch looping
           
           do ilf=1,NURBTYPE
              
              if(nint(vegt(i,j,inp)) == ILEAFCOD(ILF)) &
                   G_URBAN(i,j,inp)=float(ilf)
              
           enddo
           
           if(nint(G_URBAN(i,j,inp))/=0.)then
              
              TIBLD(i,j)=t00+TMINBLD !internal temperature defined in RAMSIN
              !     TIROAD(i,j)=  (theta(2,i,j)+theta(1,i,j))*0.5*pis  !surface
              TIROAD(i,j)=  281.16  !surface
              TSROOF(i,j) = (theta(2,i,j))*pl2            !model's first level
              TSROAD(i,j) = (theta(2,i,j)+theta(1,i,j))*0.5*pis  !surface
              TSWALL(i,j) = (TSROOF(i,j)+TSROAD(i,j))*0.5  !average surface and first level
              TCANYON(i,j)= TSWALL(i,j)        !average surface and first level
              TSTOWN(i,j)=  TSWALL(i,j)        !average surface and first level
              EMISTOWN(i,j)=0.8878604
              ALBTOWN(i,j)=0.129
              RCANYON(i,j)=(1./((1000./rv(2,i,j))+1.))*1000.
              HT(i,j)=5.
              LET(i,j)=5.
              HIN(i,j)=10.
              LEIN(i,j)=40.
              
              call INI_TG_PROFILE(TSROOF(i,j), TIBLD(i,j), &
                   TC_ROOF(1:3),D_ROOF(1:3),TROOF(2:4,i,j))
              call INI_TG_PROFILE(TSROAD(i,j), TIROAD(i,j), &
                   TC_ROAD(1:3),D_ROAD(1:3),TROAD(2:4,i,j))
              call INI_TG_PROFILE(TSWALL(i,j), TIBLD(i,j), &
                   TC_WALL(1:3),D_WALL(1:3),TWALL(2:4,i,j))
           endif
           
        enddo
        
     enddo
  enddo
  
  return
end subroutine TEB_INIT


! TEB
subroutine TEBC_INIT(ng,n2,n3,np,vegt,G_URBAN,EMIS,ALB,TS)
  
  implicit none
  integer :: n2,n3,ng,np,i,j,k
  
  real, dimension(n2,n3,np) :: vegt,G_URBAN
  real, dimension(n2,n3)    :: EMIS,ALB,TS  
  
  do i=1,n2
     do j=1,n3
        
        G_URBAN(i,j,:)=0.
        EMIS  (i,j)=0.
        ALB   (i,j)=0.
        TS    (i,j)=0.
        
     enddo
  enddo
  
  return
end subroutine TEBC_INIT
