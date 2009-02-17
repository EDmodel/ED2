MODULE rad_carma
  
  USE mem_carma
  
  INTEGER,ALLOCATABLE :: indexi(:)
  INTEGER,ALLOCATABLE :: indexj(:)  
  
  CONTAINS

!kmlnew
  SUBROUTINE radcarma(m1,m2,m3,ia,iz,ja,jz,solfac  &
     ,theta_,pi0_,pp_,rv_,RAIN_,LWL_,IWL_,dn0_,rtp_,fthrd_,rtgt_,f13t_,f23t_ &
     ,glat_,glon_,rshort_,rlong_,albedt_,cosz_,rlongup_ &
     ,mynum,fmapt_,pm_,aot_,xland_)
!kmlnew
    ! CATT
    use catt_start, only: CATT ! INTENT(IN)

    USE mem_grid,   ONLY:  centlon,	& !INTENT()
  			   dzm, 	  & !INTENT()
  			   dzt, 	  & !INTENT()
  			   idatea,	  & !INTENT()
  			   imontha,	  & !INTENT()
  			   itimea,	  & !INTENT()
  			   itopo,	  & !INTENT()
  			   iyeara,	  & !INTENT()
  			   ngrid,	  & !INTENT()
  			   nzp, 	  & !INTENT()
  			   plonn,	  & !INTENT()
  			   time 	    !INTENT()
  
    USE grid_dims, ONLY: nzpmax	    !INTENT()
    USE mem_radiate, ONLY: lonrad	    !INTENT()
  
    USE rconstants,  ONLY: cp,  	  & !INTENT()
  			   cpor,	  & !INTENT()
  			   p00, 	  & !INTENT()
  			   pio180,	  & !INTENT()
  			   stefan	    !INTENT()

    USE mem_globrad, ONLY: rad_data_not_read,raddatfn !,read_rad_data ! not used

    USE mem_aerad, ONLY: ngas,nwave,iprocopio
!kmlnew 
    USE mem_leaf          , only: leaf_g
!kmlnew  
   !LFR->para salvar o estado da memoria
    !LFR  USE rams_rad_state, ONLY: WRITE_radiate_stat
    !LFR->  
  
      IMPLICIT NONE
  
    INTEGER,INTENT(IN) :: m1,m2,m3,ia,iz,ja,jz,mynum
    REAL,INTENT(IN)    :: solfac
  	 
    REAL,INTENT(IN) :: rtgt_(m2,m3)
    REAL,INTENT(IN) :: f13t_(m2,m3)
    REAL,INTENT(IN) :: f23t_(m2,m3)
    REAL,INTENT(IN) :: glat_(m2,m3)
    REAL,INTENT(IN) :: glon_(m2,m3)
    REAL,INTENT(IN) :: cosz_(m2,m3)
    REAL,INTENT(IN) :: albedt_(m2,m3)
    REAL,INTENT(IN) :: fmapt_(m2,m3)
  
    REAL,INTENT(IN) :: pm_(m1,m2,m3)   ! particulate material (kg[pm]/kg[air])
    REAL,INTENT(IN) :: theta_(m1,m2,m3)
    REAL,INTENT(IN) :: pi0_(m1,m2,m3)
    REAL,INTENT(IN) :: pp_(m1,m2,m3)
    REAL,INTENT(IN) :: rv_(m1,m2,m3)
!kmlnew
    REAL,INTENT(IN) :: RAIN_(m2,m3)
    REAL,INTENT(IN) :: LWL_(m1,m2,m3)
    REAL,INTENT(IN) :: IWL_(m1,m2,m3)
    REAL,INTENT(IN) :: xland_(m2,m3)
!kmlnew    
    REAL,INTENT(IN) :: dn0_(m1,m2,m3)
    REAL,INTENT(IN) :: rtp_(m1,m2,m3)
  
    REAL,INTENT(INOUT) :: fthrd_(m1,m2,m3)  
    REAL,INTENT(INOUT) :: rshort_(m2,m3)
    REAL,INTENT(INOUT) :: rlong_(m2,m3)
    REAL,INTENT(INOUT) :: rlongup_(m2,m3)
    REAL,INTENT(INOUT) :: aot_(m2,m3,nwave)
    
    
    REAL :: rtgt((iz-ia+1)*(jz-ja+1))
    REAL :: f13t((iz-ia+1)*(jz-ja+1))
    REAL :: f23t((iz-ia+1)*(jz-ja+1))
    REAL :: glat((iz-ia+1)*(jz-ja+1))
    REAL :: glon((iz-ia+1)*(jz-ja+1))
    REAL :: cosz((iz-ia+1)*(jz-ja+1))
    REAL :: albedt((iz-ia+1)*(jz-ja+1))
    REAL :: fmapt((iz-ia+1)*(jz-ja+1))
  
    REAL :: pm((iz-ia+1)*(jz-ja+1),m1)! particulate material (kg[pm]/kg[air])
    REAL :: theta((iz-ia+1)*(jz-ja+1),m1)
    REAL :: pi0((iz-ia+1)*(jz-ja+1),m1)
    REAL :: pp((iz-ia+1)*(jz-ja+1),m1)
    REAL :: rv((iz-ia+1)*(jz-ja+1),m1)
!kmlnew
    REAL :: RAIN((iz-ia+1)*(jz-ja+1))
    REAL :: LWL((iz-ia+1)*(jz-ja+1),m1)
    REAL :: IWL((iz-ia+1)*(jz-ja+1),m1)
    REAL :: xland((iz-ia+1)*(jz-ja+1))
!kmlnew    
    REAL :: dn0((iz-ia+1)*(jz-ja+1),m1)
    REAL :: rtp((iz-ia+1)*(jz-ja+1),m1)
  
    REAL :: fthrd((iz-ia+1)*(jz-ja+1),m1)   
    REAL :: rshort((iz-ia+1)*(jz-ja+1))     
    REAL :: rlong((iz-ia+1)*(jz-ja+1))      
    REAL :: rlongup((iz-ia+1)*(jz-ja+1))    
    REAL :: aotl((iz-ia+1)*(jz-ja+1),nwave)
      
    REAL :: prd((iz-ia+1)*(jz-ja+1),nzpmax)
    REAL :: temprd((iz-ia+1)*(jz-ja+1),nzpmax+1)
    REAL :: dn0r((iz-ia+1)*(jz-ja+1),nzpmax)
    REAL :: dztr((iz-ia+1)*(jz-ja+1),nzpmax)
    REAL :: pmr((iz-ia+1)*(jz-ja+1),nzpmax)
    REAL :: rvr((iz-ia+1)*(jz-ja+1),nzpmax)
!kmlnew
    REAL :: RAINr((iz-ia+1)*(jz-ja+1))
    REAL :: LWLr((iz-ia+1)*(jz-ja+1),nzpmax)
    REAL :: IWLr((iz-ia+1)*(jz-ja+1),nzpmax)
    REAL :: xlandr((iz-ia+1)*(jz-ja+1))
!kmlnew
    REAL :: fthrl((iz-ia+1)*(jz-ja+1),nzpmax)
    REAL :: fthrs((iz-ia+1)*(jz-ja+1),nzpmax)
    
    REAL,PARAMETER :: fcui=1.e-6     !de mg [gas/part] /kg [ar] para kg/kg
    
    
    REAL :: pird,dzsdx,dzsdy,dlon,a1,a2,dayhr,gglon,dztri
    REAL :: dayhrr,hrangl,sinz,sazmut,slazim,slangl,cosi
    INTEGER :: igas,kk,ik,iend,ij,i,j,k,nzz

    INTEGER :: ncall = 0 

    ! DEBUG-ALF
    real :: sum_test
    !


    iend=(iz-ia+1)*(jz-ja+1) !Size of vector
    
    CALL AllocIndex(ia,ja,iz,jz,1) !1 to alloc auxiliar
    
    CALL C_2d_1d(rtgt_(ia:iz,ja:jz),rtgt,ia,iz,ja,jz,iend)
    CALL C_2d_1d(f13t_(ia:iz,ja:jz),f13t,ia,iz,ja,jz,iend)
    CALL C_2d_1d(f23t_(ia:iz,ja:jz),f23t,ia,iz,ja,jz,iend)
    CALL C_2d_1d(glat_(ia:iz,ja:jz),glat,ia,iz,ja,jz,iend)
    CALL C_2d_1d(glon_(ia:iz,ja:jz),glon,ia,iz,ja,jz,iend)
    CALL C_2d_1d(cosz_(ia:iz,ja:jz),cosz,ia,iz,ja,jz,iend)
    CALL C_2d_1d(albedt_(ia:iz,ja:jz),albedt,ia,iz,ja,jz,iend)
    CALL C_2d_1d(fmapt_(ia:iz,ja:jz),fmapt,ia,iz,ja,jz,iend)

    if (CATT==1) then
       CALL C_3d_2d(pm_(:,ia:iz,ja:jz),pm,m1,ia,iz,ja,jz,iend)
    endif

    CALL C_3d_2d(theta_(:,ia:iz,ja:jz),theta,m1,ia,iz,ja,jz,iend)
    CALL C_3d_2d(pi0_(:,ia:iz,ja:jz),pi0,m1,ia,iz,ja,jz,iend)
    CALL C_3d_2d(pp_(:,ia:iz,ja:jz),pp,m1,ia,iz,ja,jz,iend)
    CALL C_3d_2d(rv_(:,ia:iz,ja:jz),rv,m1,ia,iz,ja,jz,iend)
!kmlnew
    CALL C_3d_2d(LWL_(:,ia:iz,ja:jz),LWL,m1,ia,iz,ja,jz,iend)
    CALL C_3d_2d(IWL_(:,ia:iz,ja:jz),IWL,m1,ia,iz,ja,jz,iend)
    CALL C_2d_1d( RAIN_(ia:iz,ja:jz),RAIN ,ia,iz,ja,jz,iend)
    CALL C_2d_1d(xland_(ia:iz,ja:jz),xland,ia,iz,ja,jz,iend)
!kmlnew

    CALL C_3d_2d(dn0_(:,ia:iz,ja:jz),dn0,m1,ia,iz,ja,jz,iend)
    CALL C_3d_2d(rtp_(:,ia:iz,ja:jz),rtp,m1,ia,iz,ja,jz,iend)
    
    CALL C_2d_1d(rshort_(ia:iz,ja:jz),rshort,ia,iz,ja,jz,iend)       
    CALL C_2d_1d(rlong_(ia:iz,ja:jz),rlong,ia,iz,ja,jz,iend)         
    CALL C_2d_1d(rlongup_(ia:iz,ja:jz),rlongup,ia,iz,ja,jz,iend)     
    CALL C_3d_2d(fthrd_(:,ia:iz,ja:jz),fthrd,m1,ia,iz,ja,jz,iend)    
    CALL Ci_3d_2d(aot_(ia:iz,ja:jz,:),aotl,nwave,ia,iz,ja,jz,iend)
  
      
    ! liga radiacao de onda longa
    ir_aerad = 1
    
!srf ---- otimizacao (comentar o "call end_carma" abaixo)
    !IF (nCALL == 0) THEN
    !  nCALL = 1
    !  CALL init_carma(ia,iz,ja,jz,m1,m2,m3)
    !ENDIF
!    CALL init_carma(ia,iz,ja,jz,m1,m2,m3)
  
    IF(rad_data_not_read) THEN
       rad_data_not_read=.false.
       CALL init_carma(ia,iz,ja,jz,m1,m2,m3)
       CALL setupbins
    END IF
    
    isl_aerad = 0
    nzz = m1 - 1
    DO k = 1,m1
      DO ij=1,iend
  	pird = (pp(ij,k) + pi0(ij,k)) / cp
  	temprd(ij,k) = theta(ij,k) * pird ! air temperature (K)
  	rvr(ij,k) = max(0.,rv(ij,k))
        
	LWL(ij,k) = max(0.,LWL(ij,k))
	IWL(ij,k) = max(0.,IWL(ij,k))
	
  	! Convert the next 7 variables to cgs for now.
  	prd(ij,k) = pird ** cpor * p00 * 10. ! pressure
  	dn0r(ij,k) = dn0(ij,k) * 1.e-3        ! air density
  	dztr(ij,k) = dzt(k) / rtgt(ij) * 1.e-2

        if (CATT==1) then
           pmr(ij,k) = pm(ij,k)*fcui
        else
           pmr(ij,k) = 0.
        endif
	
      END DO
    END DO
    
    DO ij=1,iend
      temprd(ij,1) = (rlongup(ij) / stefan) ** 0.25     
      temprd(ij,nzp+1) = temprd(ij,nzp)
      !  Initialize atmospheric structure. 
      p_surf(ij) = 0.5*(prd(ij,1) + prd(ij,2))    
      p_top(ij)  = prd(ij,m1)		  
      t_surf(ij) = temprd(ij,1)   
    END DO
    DO k=1,m1-1
      DO ij=1,iend
  	! K level in CARMA grid corresponds to K+1 level in BRAMS grid
  	! Transfer values from RAMS grid to CARMA grid
  	p(ij,k)    =	prd(ij,k+1)
  	t(ij,k)    = temprd(ij,k+1)
  	rhoa(ij,k) =   dn0r(ij,k+1)
!kmlnew
	LWLr(ij,k) = LWL(ij,k+1)*dn0r(ij,k+1) * 1.e+3  ![kg/m3]
	IWLr(ij,k) = IWL(ij,k+1)*dn0r(ij,k+1) * 1.e+3  ![kg/m3]
!kmlnew	
      END DO
    END DO
    
 
    DO ik = 1,NZZ
      !  Reverse the vertical index when in cartesian coordinates
       kk = NZZ + 1 - ik
      DO ij=1,iend
  	 t_aerad(ij,kk) = t(ij,ik)
  	 p_aerad(ij,kk) = p(ij,ik)
!kmlnew
	LWL_aerad(ij,kk) = LWLr(ij,ik)
	IWL_aerad(ij,kk) = IWLr(ij,ik)
	
!        if(LWL_aerad(ij,ik).ne.0.) print*,'LWL [g/cm3]=', LWL_aerad(ij,ik), ik
!	if(IWL_aerad(ij,ik).ne.0.) print*,'IWL [g/cm3]=', IWL_aerad(ij,ik), ik
!kmlnew	
	 
      END DO
    END DO

!kmlnew 

    DO ik = 1,NZZ
      DO ij=1,iend
        dztri=1./(dztr(ij,ik) * 1.e+2)
	LWP_aerad(ij,ik) = LWL_aerad(ij,ik) * dztri   ![kg/m2]
	IWP_aerad(ij,ik) = IWL_aerad(ij,ik) * dztri   ![kg/m2]
!	if(IWL_aerad(ij,ik).gt.0.) print*,'ik,IWL=',ik,IWL_aerad(ij,ik),IWP_aerad(ij,ik)
      END DO
    END DO
!kmlnew 
      DO ij=1,iend
       xland_aerad(ij)=xland(ij)
       RAIN_aerad(ij)=RAIN(ij)
      END DO
   
    
    DO ij=1,iend
      tabove_aerad(ij)  = t(ij,nzz)
    END DO
    
    !  Initialize gas concentrations.
    DO igas = 1,ngas
       DO k = 1,m1-1
  	  ! K level in CARMA grid corresponds to K+1 level in BRAMS grid
  	  DO ij=1,iend
  	    ! water vapor concentration
  	    IF( igas .eq. 1 ) gc(ij,k,igas) = rvr(ij,k+1) * dn0r(ij,k+1) 
  	  END DO
       END DO
    END DO
  
    DO ij=1,iend
      ! The shortwave parameterizations are only valid if the cosine
      !    of the zenith angle is greater than .03 .
      IF (cosz(ij) .gt. .03) isl_aerad(ij) = 1
    END DO
      
    CALL initaer(m1,pmr,dn0r,ia,iz,ja,jz,nzpmax)
    
    !  Initialize radiation
![ED2
    CALL initrad(imontha,idatea,iyeara,itimea,time,m1,ia,ja,iz,jz)
!ED2]
    CALL prerad(m1,dztr,fmapt,ia,iz,ja,jz,nzpmax,m2,m3)
    
    CALL radtran(albedt,cosz,m1,m2,m3,ia,iz,ja,jz,aotl(:,11))  !kml2
  
    CALL radtran_to_rams(nzp,m2,m3,fthrl,rlong,fthrs,rshort,aotl,ia,iz,ja,jz,mynum)   
  
    ! 
    ! Modify the DOwnward surface shortwave flux by considering
    !	 the slope of the topography.
  
    DO ij=1,iend
      IF (itopo .eq. 1) THEN
  	dzsdx = f13t(ij) * rtgt(ij)
  	dzsdy = f23t(ij) * rtgt(ij)
  
  	! The y- and x-directions must be true north and east for
  	! this correction. the following rotates the model y/x
  	! to the true north/east.   
  
  	! The following rotation seems to be incorrect,so CALL this instead:
  	! SUBROUTINE uvtoueve(u,v,ue,ve,qlat,qlon,platn(ngrid),plonn(ngrid))
  
  	dlon = (plonn(ngrid) - glon(ij)) * pio180
  	a1 = dzsdx*cos(dlon) + dzsdy * sin(dlon)
  	a2 = -dzsdx*sin(dlon) + dzsdy * cos(dlon)
  	dzsdx = a1
  	dzsdy = a2
  
  	dayhr = real(time / 3600.) + float(itimea/100)  &
  	    + float(mod(itimea,100)) / 60.
  	gglon = glon(ij)
  	IF (lonrad .eq. 0) gglon = centlon(1)
  	dayhrr = mod(dayhr+gglon/15.+24.,24.)
  	hrangl = 15. * (dayhrr - 12.) * pio180
        !srf - evitando SQRT (<0)
        !sinz = sqrt(1. - cosz(ij) ** 2)
  	sinz = sqrt(max(0., (1. - cosz(ij) ** 2)))
  
  	! ALF - Evitando divisao por zero
  	sinz = max(0.000001, sinz)
  	! ALF
  
  	sazmut = asin(max(-1.,min(1.,cdec*sin(hrangl)/sinz)))
  	IF (abs(dzsdx) .lt. 1e-20) dzsdx = 1.e-20
  	IF (abs(dzsdy) .lt. 1e-20) dzsdy = 1.e-20
  	slazim = 1.571 - atan2(dzsdy,dzsdx)
  	slangl = atan(sqrt(dzsdx*dzsdx+dzsdy*dzsdy))
  	cosi = cos(slangl) * cosz(ij) + sin(slangl) * sinz  &
  	     * cos(sazmut-slazim)
  	rshort(ij) = rshort(ij) * cosi / cosz(ij)        
     END IF
    END DO
    
    !print*,'------ radiative heating rates ---- ---'
    DO k = 2,m1-1
      DO ij=1,iend
  	 fthrd(ij,k) = fthrl(ij,k) + fthrs(ij,k)        
!	 if(fthrd(ij,k)*86400. < -10.) &
!	 print*,'IJ K',ij,k,fthrl(ij,k)*86400. , fthrs(ij,k)*86400.
      END DO
    END DO
  
    ! Convert the downward flux at the ground to SI.
  
    !	     rshort(i,j) = rshort(i,j) * 1.e-3 / (1. - albedt(i,j))
    !	     rlong(i,j) = rlong(i,j) * 1.e-3
    DO ij=1,iend
      rshort(ij) = rshort(ij) / (1. - albedt(ij)) 
      rlong(ij) = rlong(ij)                       
      fthrd(ij,1) = fthrd(ij,2)                  
      !print*,'IJ SW LW=',ij,rshort(ij),rlong(ij)
      !call flush(6)
      !if (rlong(ij).lt.10.) print*,'Antes da conversao!!!',ij,rlong(ij)

    END DO
  
    
!srf ---- comentar a linha abaixo    
    !!! CALL end_carma()
  
    CALL C_1d_2d(rshort,rshort_(ia:iz,ja:jz),ia,iz,ja,jz,iend)    
    CALL C_1d_2d(rlong,rlong_(ia:iz,ja:jz),ia,iz,ja,jz,iend)      
    CALL C_1d_2d(rlongup,rlongup_(ia:iz,ja:jz),ia,iz,ja,jz,iend)  
    CALL C_2d_3d(fthrd,fthrd_(:,ia:iz,ja:jz),m1,ia,iz,ja,jz,iend) 
    CALL Ci_2d_3d(aotl,aot_(ia:iz,ja:jz,:),nwave,ia,iz,ja,jz,iend)
  
    CALL AllocIndex(ia,ja,iz,jz,0) !0 to dealloc auxiliar
    
   
  
  END SUBROUTINE radcarma
  
  SUBROUTINE setupbins
    !  This routine evaluates the derived mapping arrays and sets up
    !  the particle size bins.
  
    USE mem_aerad, ONLY: lunoprt, nbin
    USE mem_globaer, ONLY: ngroup,nelem,itype,i_involatile, &
  			   i_volatile,ienconc,igelem,ncore, &
  			   nelemg,i_coremass,i_volcore, &
  			   i_core2mom,ixyz,nxyz,rhop3,rhoelem,rhopcore3, &
  			   rhocore,pi,rmassmin,rmin,rmrat,one,rmass, &
  			   rmasscore,pcore,rmassup,rmasscoreup,dm,vol, &
  			   r,rcore,rup,rcoreup,dr,rlow,diffmass
    
    IMPLICIT NONE
    
    !Local
  
    INTEGER :: igrp
    INTEGER :: ielem
    INTEGER :: j
    INTEGER :: ie
    INTEGER :: ig
    INTEGER :: ibin
    REAL    :: cpi
    REAL    :: vrfact
  
    !  Determine which elements are particle number concentrations
    !  <ienconc(igroup)> is the element corresponding to particle number 
    !  concentration in group <igroup>
    !
    igrp = 0
    DO ielem = 1, NELEM
       IF( itype(ielem) .eq. I_INVOLATILE .or. &
  	    itype(ielem) .eq. I_VOLATILE )THEN
  
  	  igrp = igrp + 1
  	  ienconc(igrp) = ielem
       END IF
    END DO
    !
    !  Determine which group each element belongs to
    !  i.e., <igelem(ielem)> is the group to which element <ielem> belongs
    !
    igrp = 0
    DO ielem = 1, NELEM
       IF( itype(ielem) .eq. I_INVOLATILE .or.       &
  	    itype(ielem) .eq. I_VOLATILE )THEN
  	  igrp = igrp + 1
       END IF
       igelem(ielem) = igrp
    END DO
    !
    !  Particle mass densities (NXYZ*NBIN for each group) -- the user might want
    !  to modIFy this (this code segment DOes not appear in setupaer SUBROUTINE
    !  because <igelem> is not defined until this SUBROUTINE).
    !
    DO ie = 1,NELEM
       ig = igelem(ie)
       DO ibin = 1,NBIN
  	  DO ixyz = 1,NXYZ
  	     rhop3(ixyz,ibin,ig) = rhoelem(ie)
  	     rhopcore3(ixyz,ibin,ig) = rhocore(ie)
  	  END DO
       END DO
    END DO
    !
    !
    !  Set up the particle bins.
    !  For each particle group, the mass of a particle in
    !  bin j is <rmrat> times that in bin j-1
    !
    !	 rmass(NBIN,NGROUP)	=  bin center mass [g]
    !	 r(NBIN,NGROUP) 	=  bin mean (volume-weighted) radius [cm]
    !	 vol(NBIN,NGROUP)	=  bin center volume [cm^3]
    !	 dr(NBIN,NGROUP)	=  bin width in radius space [cm]
    !	 dv(NBIN,NGROUP)	=  bin width in volume space [cm^3]
    !	 dm(NBIN,NGROUP)	=  bin width in mass space [g]
    !
    cpi = 4./3.*PI
  
    DO igrp = 1, NGROUP
  
       rmassmin(igrp) = cpi*rhop3(1,1,igrp)*rmin(igrp)**3
  
       vrfact = ( (3./2./PI/(rmrat(igrp)+1.))**(ONE/3.) )*    &
  	    ( rmrat(igrp)**(ONE/3.) - 1. )
  
       DO j = 1, NBIN
  	  !PRINT *,'LFRDBG->igrp,j,rmassmin(igrp),rmrat(igrp):', &
  	  !	    igrp,j,rmassmin(igrp),rmrat(igrp)
  	  !CALL FLUSH(6)
  	  rmass(j,igrp)   = rmassmin(igrp) * rmrat(igrp)**(j-1)
  	  rmasscore(j,igrp) = pcore/100. * rmass(j,igrp)
  
  	  rmassup(j,igrp) = 2.*rmrat(igrp)/(rmrat(igrp)+1.)*rmass(j,igrp)
  	  rmasscoreup(j,igrp) = pcore/100. * rmassup(j,igrp)
  
  	  dm(j,igrp)	  = 2.*(rmrat(igrp)-1.)/(rmrat(igrp)+1.)*rmass(j,igrp)
  	  vol(j,igrp) = rmass(j,igrp) / rhop3(1,1,igrp)
  
  	  r(j,igrp)	  = ( rmass(j,igrp)/rhop3(1,1,igrp)/cpi )**(ONE/3.)
  	  !PRINT *,'LFRDBG->cpi,ONE,rmasscore(j,igrp),rhopcore3(1,1,igrp):', &
  	  !	    cpi,one,rmasscore(j,igrp),rhopcore3(1,1,igrp)
  	  !CALL FLUSH(6)
  	  rcore(j,igrp)   = ( rmasscore(j,igrp)/rhopcore3(1,1,igrp)/cpi )**(ONE/3.)
  
  	  rup(j,igrp)	  = ( rmassup(j,igrp)/rhop3(1,1,igrp)/cpi )**(ONE/3.)
  	  rcoreup(j,igrp) = ( rmasscoreup(j,igrp)/rhopcore3(1,1,igrp)/cpi )**(ONE/3.)
  	  !PRINT *,'LFRDBG->vrfact,rmass(j,igrp),rhop3(1,1,igrp):', &
  	  !		    vrfact,rmass(j,igrp),rhop3(1,1,igrp)
  
  	  dr(j,igrp)  = vrfact*(rmass(j,igrp)/rhop3(1,1,igrp))**(ONE/3.)
  	  rlow(j,igrp) = rup(j,igrp) - dr(j,igrp)
       END DO
    END DO
    !PRINT *,'LFRDBG->End of setupbins';CALL FLUSH(6)
  
  END SUBROUTINE setupbins
  
  SUBROUTINE initaer(m1,pmr,dn0r,ia,iz,ja,jz,nzpmax)

    use mem_aerad, only: &
         nbin

    USE mem_globaer, ONLY: nelem      , &
  			   igelem     , &
  			   ienconc    , &
  			   small_pc   , &
  			   itype      , &
  			   i_coremass , &
  			   rmass      , &
  			   fix_coref  , &
  			   i_core2mom , &
  			   rhop3      , &
  			   pi	      , &
  			   dr	      , &
  			   r
    
    IMPLICIT NONE
  
    INTEGER,INTENT(IN)  	      :: m1,ia,iz,ja,jz,nzpmax
    REAL   ,INTENT(IN), DIMENSION((iz-ia+1)*(jz-ja+1),nzpmax) :: pmr
    REAL   ,INTENT(IN), DIMENSION((iz-ia+1)*(jz-ja+1),nzpmax) :: dn0r
  
    !Local
    REAL,DIMENSION((iz-ia+1)*(jz-ja+1),m1) :: totm
    REAL,DIMENSION((iz-ia+1)*(jz-ja+1),m1) :: r0
    REAL,DIMENSION((iz-ia+1)*(jz-ja+1),m1) :: rsig
    INTEGER :: ie,ix,iy
    INTEGER :: ielem
    INTEGER :: ig
    INTEGER :: ip,iend,ij
    INTEGER :: j
    INTEGER :: k
    INTEGER :: kr
    INTEGER :: nzz
    REAL    :: arg1
    REAL    :: arg2
    REAL    :: sum((iz-ia+1)*(jz-ja+1))
    REAL    :: totn((iz-ia+1)*(jz-ja+1))
  
    iend=(iz-ia+1)*(jz-ja+1)
      
    nzz = m1 - 1
    !
    !transfere valores da grade do rams para carma
    DO k = 1,nzz
      kr = K + 1     ! nivel K da grade do carma orresponde ao nivel K + 1 DO RAMS
      DO ij=1,iend
  	!    totm = total mass particle concentration (g/cm3) 
  	totm(ij,k)  = pmr(ij,kr) *  dn0r(ij,kr) 
      END DO
    END DO
    !  Initialize particle number densities 
    !  Core mass is assumed to be 100% of particle mass
    DO ielem = 1,nelem
       ig = igelem(ielem)
       ip = ienconc(ig)
       DO j = 1,nbin
  	  DO k = 1,nzz
  	     DO ij = 1,iend
  	       IF( ielem .eq. ip )THEN
  		  !  Particle number concentration [#/cm^3]
  		  pc(ij,k,j,ielem) = SMALL_PC
  	       ELSE IF( itype(ielem) .eq. I_COREMASS )THEN
  		  !  Core mass concentration [g/cm^3]
  		  pc(ij,k,j,ielem) = pc(ij,k,j,ip)*rmass(j,ig) * &
  					FIX_COREF
  	       ELSE IF( itype(ielem) .eq. I_CORE2MOM )THEN
  		  !  Second moment of core mass distribution [ (g/cm^3)^2 ]
  		  pc(ij,k,j,ielem) = pc(ij,k,j,ip) *	      &
  		       (rmass(j,ig)*FIX_COREF)**2
  	       END IF
  	    END DO
  	  END DO
       END DO
    END DO
    !
    !  Initial particle distribution: log-normal size distribution 
    !  for first particle group (which has only one particle element)
    !  in a single column
    ig = 1
    ie = ienconc(ig)
  
    DO k = 1,nzz
       !
       !  Log-normal parameters:
       !  
       !    r0   = number mode radius
       !    rsig = geometric standard deviation
       !    totm = total mass particle concentration (g/cm3) (proveniente DO rams) 
       !
       DO ij=1,iend
  	 r0(ij,k)  = 1.95e-5
  	 rsig(ij,k) = 1.62
  	 totn(ij) = (6. * totm(ij,k)/(rhop3(1,1,1)*PI*r0(ij,k)**3))* &
  	    exp((-9./2)*log(rsig(ij,k))**2)
       END DO
       !  Adjust prefactor to yield particle number concentration <ntot>
       !
     sum = 0.
     DO j = 1,nbin
        DO ij=1,iend
  	  arg1 = dr(j,ig) / ( sqrt(2.*PI) * r(j,ig) * log(rsig(ij,k)) ) 
  	  arg2 = -log( r(j,ig) / r0(ij,k) )**2 / &
    		       ( 2.*log(rsig(ij,k))**2 )

  	  sum(ij)  = sum(ij) + arg1 * exp( arg2 )
        END DO
     END DO
     DO ij=1,iend
       totn(ij) = totn(ij) / sum(ij)
     END DO
     DO j = 1,nbin
        DO ij=1,iend
  	  arg1 = totn(ij) * dr(j,ig) / ( sqrt(2.*PI) * r(j,ig) * &
  		       log(rsig(ij,k)) ) 

  	  arg2 = -log( r(j,ig) / r0(ij,k) )**2 / &
  		       ( 2.*log(rsig(ij,k))**2 )
  	  pc(ij,k,j,ie) = max( arg1 * exp( arg2 ), REAL(SMALL_PC) )
        END DO
     END DO
    END DO
    
  END SUBROUTINE initaer
  
  SUBROUTINE initrad(imonth1,idate1,iyear1,itime1,time_rams,m1,ia,ja,iz,jz)
  
    USE mem_aerad, ONLY: is_grp_ice_aerad,r_aerad, &
  			 rup_aerad,rcore_aerad,rcoreup_aerad, &
  			 ptop_aerad,pbot_aerad,u0_aerad, &
  			 sfc_alb_aerad,emisir_aerad,tsfc_aerad, &
  			 tabove_aerad,wave_aerad,iprocopio&
                         ,nx,ny,nbin,nwave,nsol
  
    USE mem_globaer, ONLY: time,do_solar,do_ir,isolar_zen,i_diurnal, &
  			   rad_start,scday,pi,ix,iy,rlat,u0,ngroup, &
  			   is_grp_ice,ienconc,r,rup,rcore,rcoreup, &
  			   t_surf,wave,z_sin, &
  			   z_cos
    USE mem_globrad, ONLY: imie
    
    USE mem_carma, ONLY: declin
    
    IMPLICIT NONE
    
    INTEGER,INTENT(IN) :: m1,ia,ja,iz,jz
    INTEGER,INTENT(IN) :: imonth1
    INTEGER,INTENT(IN) :: idate1
    INTEGER,INTENT(IN) :: iyear1
    INTEGER,INTENT(IN) :: itime1
    REAL(kind=8)   ,INTENT(IN) :: time_rams
  
    !Local  
    INTEGER :: iday
    INTEGER :: iwave
    INTEGER :: julday
    REAL    :: saz
    REAL    :: wavetemp
   
    !
    !  Define flag to control the calculation of the solar zenith angle:
    !	 <isolar_zen> = I_FIXED: USE fixed value <u0fixed>
    !		      = I_DIURNAL: calculation based on time, day, lat, and lon
    !
    isolar_zen = I_DIURNAL
  
    IF( isolar_zen .eq. I_DIURNAL )THEN
       !
       !
       !  Define values needed for calculation of solar zenith angle:
       !    <iday> is day of year
       !    <rad_start> is solar time corresponding to <time> = 0, in seconds
       !     = 0 means <time> = 0 corresponds to midnight,
       !     = 6 * 3600 means <time> = 0 corresponds to 6 AM
       !    Note: all times are local standard time.
       !
  
       iday =  julday(imonth1,idate1,iyear1)
       iday =  iday + nint(time_rams/86400.)
       rad_start = (float(itime1/100) + float(mod(itime1,100)) / 60.)*scday
       !
       !
       !  Precalculate terms in solar zenith angle computation:
       !    (adapted from original Toon model)
       !    <saz> is solar azimuth angle [rad]
       !    <declin> is solar declination [rad]
       !    <z_sin> is sin term in precalculation
       !    <z_cos> is cos term in precalculation
       !
  
       saz = 2. * PI / 365. * iday 
  
       declin = 0.006918 - 0.399912*cos(saz)	+0.070257*sin(saz)    &
  	    - 0.006758*cos(2.*saz) +0.000907*sin(2.*saz) &
  	    - 0.002697*cos(3.*saz) +0.001480*sin(3.*saz)
  
       !DO ix = 1,NX
       !   DO iy = 1,NY
       !      rlat(ix,iy) = glat
       !      z_sin(ix,iy) = sin(declin) * sin( rlat(ix,iy) * PI/180. )
       !      z_cos(ix,iy) = cos(declin) * cos( rlat(ix,iy) * PI/180. )
       !   END DO
       !END DO
  
    END IF
    !
  
    !
    !  Initialize the radiative transfer model
    CALL setuprad(m1,ia,ja,iz,jz)
    
     IF(imie == 0) THEN
      CALL calcproperties
     END IF

    !  Get radiative wavelengths
    !
    DO iwave = 1,NWAVE
       !
       !
       !  Solar wavelengths in radiative transfer model are bin centers,
       !  infrared are bin edges
       !
       IF( iwave .le. NSOL )THEN
  	  wave(iwave) = wave_aerad(iwave)
       ELSE
  	  wave(iwave) = 0.5*( wave_aerad(iwave) + wave_aerad(iwave+1) )
       END IF
  
    END DO
    !
    !  Switch bins 11 and 12
    !KLF Corrigidos os comp. de onda (11 e 12) para (17 e 18)!!!
    !	   wavetemp = wave(11)
    !	   wave(11) = wave(12)
    !	   wave(12) = wavetemp
  
    wavetemp = wave(17)
    wave(17) = wave(18)
    wave(18) = wavetemp
    !
    !  Return to CALLer with radiation model initialized
    !
  END SUBROUTINE initrad
  
  SUBROUTINE setuprad(m1,ia,ja,iz,jz)
    !	  *********************************************************
    !	  *  Purpose		:  Defines all constants, and	  *
    !	  *			   calculates pressure averaged   *
    !	  *			   absorption coefficients.	  *
    !	  * *******************************************************
    !
  
  
    USE mem_aerad, ONLY: wave_aerad,u0_aerad, &
  			 sfc_alb_aerad,emisir_aerad,ptop_aerad, &
  			 pbot_aerad,tsfc_aerad
  
    USE mem_globrad, ONLY: nlayer,g, &
  			   ntotal,nprob, &
  			   wave, &
  			   nvert,p,t,o2mol,am, &
  			   co2mol,o3c,o3mol,avg,wol,gol, &
  			   tauray,nsolp,ltemp,nsol,alos,aco2,xaco2, &
  			   ao2,xao2,ao3,xao3,xah2o,psh2o, &
  			   psco2,pso2,pso3,pj,o3mixp, &
  			   akh2o,ako3,akco2,nirp,imie, &
  			   corereal,coreimag
     
    IMPLICIT NONE
   
    INTEGER,INTENT(IN) :: m1,ia,ja,iz,jz
    
    REAL,DIMENSION((iz-ia+1)*(jz-ja+1),m1) :: pbar
    REAL,DIMENSION((iz-ia+1)*(jz-ja+1),m1) :: o3mix
    INTEGER :: i,i1,j1
    INTEGER :: ii((iz-ia+1)*(jz-ja+1),nlayer)
    INTEGER :: ik((iz-ia+1)*(jz-ja+1),nlayer)
    INTEGER :: ij
    INTEGER :: j
    INTEGER :: k
    INTEGER :: l
    REAL    :: co2mix
    REAL    :: dp((iz-ia+1)*(jz-ja+1),nlayer)
    REAL    :: o2mix
    REAL,DIMENSION((iz-ia+1)*(jz-ja+1))  :: o3mix2
    REAL    :: pm
    REAL    :: ps((iz-ia+1)*(jz-ja+1),nlayer)
    REAL    :: wvo
    REAL    :: x
    INTEGER :: iend
    
    iend=(iz-ia+1)*(jz-ja+1)
  
    !pbar  - layer average pressure (bars)
    !(note - the top layer is from ptop to 0, so average = ptop/2)
    !press - pressure at edge of layer (dyne/cm^2)
    !dpg   - mass of layer (g / cm**2)
    DO ij=1,iend
      pbar(ij,1)  = p_top(ij)/2.0E6
      press(ij,1) = p_top(ij)
    END DO
    DO  k  = 2,nvert
      DO ij=1,iend 
  	pbar(ij,k)  = p_aerad(ij,k-1)/1.0E6
  	press(ij,k) = (p_aerad(ij,k-1) + p_aerad(ij,k)) * 0.5
  	dpg(ij,k-1) = (press(ij,k)-press(ij,k-1)) / g
      END DO
    END DO
    DO ij=1,iend
      pbar(ij,nlayer)  = p_aerad(ij,nvert)/1.0E6
      press(ij,nlayer) = p_surf(ij)
      dpg(ij,nvert)  = (press(ij,nlayer)-press(ij,nvert)) / g
      !skin temperature
      tt(ij,nlayer) = t_surf(ij)
      !amount of water vapor above model domain (gm / cm**2)
      !From 1976 U.S. Standard Atmosphere, mid-latitude sounding:
      !
      !For z_top = 5 km
      !RDH2O(1) = .13
      !
      !For z_top = 2.6 km
      !RDH2O(1) = .64
      rdh2o(ij,1)   = h2ocol_aerad
      !interpolate temperature from layer center (t) to layer edge (tt)
      tt(ij,1) = t_aerad(ij,1)
    END DO
    DO  k = 2, nvert
      DO ij=1,iend 
  	tt(ij,k) = t_aerad(ij,k-1) * (press(ij,k)/p_aerad(ij,k-1)) ** &
  		    (LOG(t_aerad(ij,k)/t_aerad(ij,k-1))/&
  		     LOG(p_aerad(ij,k)/p_aerad(ij,k-1)))

      END DO
    END DO
    !
    !
    !	  DEFINE MASS MIXING RATIOS. O3MIX TAKEN FROM U.S. STANDARD ATMOS-
    !	  PHERE, MID-LATITUDE SOUNDING
    !
    o2mix	   =   0.22*o2mol/am
    co2mix	   =   3.5E-4*co2mol/am
    !
    !	  OZONE COLUMN ABUNDANCE O3C (#/CM**2) ABOVE PTOP WAS CALCULATED
    !	  FROM THE U.S. STANDARD ATMOSPHERE MID-LATITUDE PROFILE.
    !
    !	  This is for z_top = 5 km
    o3c 	   =   9.02E18
    !
    !	  This is for z_top = 3 km
    !	  O3C		 =   9.2E18
    !
    !	  This is for z_top = 1 km
    !	   O3C  	  =   9.3E18
    !
    !	  CONVERT O3C TO MASS MIXING RATIO O3MIX2.
    !
    DO ij=1,iend
       o3mix2(ij) = o3c*o3mol*g/(p_top(ij)*avg)
    END DO
    !  !
    DO  l	    =	nsolp+1,ntotal
      ltemp(l-nsolp)  =   nprob(l) - nsol
    END DO
    !
    x		       =   alos/avg
    !
    !	  CONVERT SOLAR ABSORPTION COEFFICIENTS TO CM**2/GM.
    !
    DO  l	    =	1,nsolp
      !srf	   ACO2(L)	   =   ACO2(L)/(X*CO2MOL)
      !srf	   AO2(L)	 =   AO2(L)/(X*O2MOL)
      !srf	   AO3(L)	 =   AO3(L)/(X*O3MOL)
      aco2(l)	      =   xaco2(l)/(x*co2mol)
      ao2(l)	      =   xao2(l)/(x*o2mol)
      !  if(l.eq.14) print*,'XAO2=',Xao2(l),X,O2MOL
      ao3(l)	      =   xao3(l)/(x*o3mol)
    END DO
    !
    !	  CALCULATE ABSORPTION COEFFICIENTS
    !
    pah2o=0.0
    paco2=0.0
    pao2=0.0
    pao3=0.0
   
    DO  j =   1,nlayer
      DO  l=  1,nsolp
  	DO ij=1,iend
  	  pah2o(ij,l,j)   =  xah2o(l)*pbar(ij,j)**psh2o(l)
  	  paco2(ij,l,j)   =   aco2(l)*pbar(ij,j)**psco2(l)
  	  pao2(ij,l,j)    =   ao2(l)*pbar(ij,j)**pso2(l)
  	  pao3(ij,l,j)    =   ao3(l)*pbar(ij,j)**pso3(l)
  	END DO
      END DO
    END DO
    !
    
    DO ij=1,iend
      DO  j	    =	1,nlayer
  	DO  i	  =   1,6
  	  ii(ij,j)	     =   i
  	  IF(pbar(ij,j) > pj(i)) EXIT
  	END DO
      END DO
    END DO
    
    ps = 0.0
    DO  j	    =	1,nlayer
      DO ij=1,iend
  	IF( ii(ij,j) == 1 ) ii(ij,j) = 2
  	dp(ij,j)	   =   LOG(pj(ii(ij,j)-1)/pj(ii(ij,j)))
  	IF( pbar(ij,j) > pj(6) )THEN
  	  ik(ij,j) = ii(ij,j) - 1
  	ELSE
  	  ik(ij,j) = ii(ij,j)
  	END IF
  	ps(ij,j) = pbar(ij,j)/pj(ik(ij,j))
  	IF (j /= 1) o3mix(ij,j) =o3mixp(ik(ij,j))*ps(ij,j)**(LOG(o3mixp(ii(ij,j)-1)/o3mixp(ii(ij,j)))/dp(ij,j))
      END DO
    END DO
    
    DO  j	=   1,nlayer
      DO  l	  =   1,31
  	DO ij=1,iend
  	  pah2o(ij,nsolp+l,j) = akh2o(l,ik(ij,j))*ps(ij,j)**(LOG (akh2o(l,ii(ij,j)-1)/akh2o(l,ii(ij,j)))/dp(ij,j))
  	END DO
      END DO
      DO  l	  =   32,35
  	DO ij=1,iend
  	  pah2o(ij,nsolp+l,j) = akh2o(32,ik(ij,j))*ps(ij,j)**(LOG (akh2o(32,ii(ij,j)-1)/akh2o(32,ii(ij,j)))/dp(ij,j))
  	  pao3(ij,nsolp+l,j)  = ako3(l-31,ik(ij,j))*ps(ij,j)**(LOG  &
  			     (ako3(l-31,ii(ij,j)-1)/ako3(l-31,ii(ij,j)))/dp(ij,j))
  	END DO
      END DO
      DO ij=1,iend
  	pah2o(ij,nsolp+36,j)  = akh2o(33,ik(ij,j))*ps(ij,j)**(LOG (akh2o(33,ii(ij,j)-1)/akh2o(33,ii(ij,j)))/dp(ij,j))
      END DO
      DO  l	  =   37,40
  	DO ij=1,iend
  	  paco2(ij,nsolp+l,j) = akco2(1,ik(ij,j))*ps(ij,j)**(LOG (akco2(1,ii(ij,j)-1)/akco2(1,ii(ij,j)))/dp(ij,j))
  	  pah2o(ij,nsolp+l,j) = akh2o(l-3,ik(ij,j))*ps(ij,j)**(LOG  &
  		(akh2o(l-3,ii(ij,j)-1)/akh2o(l-3,ii(ij,j)))/dp(ij,j))
  	END DO
      END DO
      DO  l	  =   41,44
  	DO ij=1,iend
  	  paco2(ij,nsolp+l,j)  =  akco2(2,ik(ij,j))*ps(ij,j)**(LOG (akco2(2,ii(ij,j)-1)/akco2(2,ii(ij,j)))/dp(ij,j))
  	  pah2o(ij,nsolp+l,j)  =  akh2o(l-7,ik(ij,j))*ps(ij,j)**(LOG  &
  		(akh2o(l-7,ii(ij,j)-1)/akh2o(l-7,ii(ij,j)))/dp(ij,j))
  	END DO
      END DO
      DO  l	  =   45,48
  	DO ij=1,iend
  	  paco2(ij,nsolp+l,j) = akco2(3,ik(ij,j))*ps(ij,j)**(LOG (akco2(3,ii(ij,j)-1)/akco2(3,ii(ij,j)))/dp(ij,j))
  	  pah2o(ij,nsolp+l,j) = akh2o(l-11,ik(ij,j))*ps(ij,j)**(LOG  &
  	      (akh2o(l-11,ii(ij,j)-1)/akh2o(l-11,ii(ij,j)))/dp(ij,j))
      END DO
    END DO
    DO  l	=   49,51
      DO ij=1,iend
        paco2(ij,nsolp+l,j) = akco2(4,ik(ij,j))*ps(ij,j)**(LOG (akco2(4,ii(ij,j)-1)/akco2(4,ii(ij,j)))/dp(ij,j))
        pah2o(ij,nsolp+l,j) = akh2o(l-11,ik(ij,j))*ps(ij,j)**(LOG  &
  	      (akh2o(l-11,ii(ij,j)-1)/akh2o(l-11,ii(ij,j)))/dp(ij,j))
      END DO
    END DO
    DO  l	=   52,54
      DO ij=1,iend
        paco2(ij,nsolp+l,j) = akco2(5,ik(ij,j))*ps(ij,j)**(LOG (akco2(5,ii(ij,j)-1)/akco2(5,ii(ij,j)))/dp(ij,j))
        pah2o(ij,nsolp+l,j) = akh2o(l-14,ik(ij,j))*ps(ij,j)**(LOG  &
  	      (akh2o(l-14,ii(ij,j)-1)/akh2o(l-14,ii(ij,j)))/dp(ij,j))
      END DO
    END DO
    DO  l	=   55,57
      DO ij=1,iend
        paco2(ij,nsolp+l,j) = akco2(6,ik(ij,j))*ps(ij,j)**(LOG (akco2(6,ii(ij,j)-1)/akco2(6,ii(ij,j)))/dp(ij,j))
        pah2o(ij,nsolp+l,j) = akh2o(l-17,ik(ij,j))*ps(ij,j)**(LOG  &
  	      (akh2o(l-17,ii(ij,j)-1)/akh2o(l-17,ii(ij,j)))/dp(ij,j))
      END DO
    END DO
    DO  l	=   58,nirp
      DO ij=1,iend
  	pah2o(ij,nsolp+l,j) = akh2o(l-17,ik(ij,j))*ps(ij,j)**(LOG  &
  		(akh2o(l-17,ii(ij,j)-1)/akh2o(l-17,ii(ij,j)))/dp(ij,j))
  	END DO
      END DO
    END DO
  
    DO ij=1,iend
      ! store o3mix2 in o3mix(1)
      o3mix(ij,1) = o3mix2(ij)
    END DO
    
    !	  here we find taugas. it is tauco2+tauo2+tauo3.
    DO  l=1,ntotal
      DO ij=1,iend
  	pm=p_top(ij)/g
  	taugas(ij,l,1) = pm*(o2mix*pao2(ij,l,1)+co2mix* &
  			     paco2(ij,l,1)+o3mix(ij,1)*pao3(ij,l,1)) 
      END DO
    END DO
    DO    j =	2,nlayer
      DO  l    =   1,ntotal
  	 DO ij=1,iend
  	   pm=dpg(ij,j-1)
  	   taugas(ij,l,j) = pm*(o2mix*pao2(ij,l,j)+co2mix* &
  			      paco2(ij,l,j)+o3mix(ij,j)*pao3(ij,l,j)) 
  	 END DO
      END DO
    END DO
  
    !
    !	  wave must be in microns
    !	  calculate rayleigh optical depth parameters.
    !
    DO  l      =   1,ntotal
      wvo	=    wave(nprob(l))
      tauray(l) =   (8.46E-9/wvo**4) * ( 1.+0.0113/wvo**2+0.00013/wvo**4 )
    END DO
  
    !	  we do not include rayleigh scattering in infrared
    DO  j = 1,nvert
      DO  l= 1,ntotal
  	DO ij=1,iend
  	  IF( l <= nsolp ) THEN
  	    paray(ij,l,j+1) = tauray(l)*dpg(ij,j)*g
  	  ELSE
  	    paray(ij,l,j+1) = 0.
  	  END IF
  	END DO
      END DO
      DO  l   =   1,ntotal
  	DO ij=1,iend
  	  IF( l <= nsolp ) THEN
  	    paray(ij,l,1) = tauray(l)*p_top(ij)
  	  ELSE
  	    paray(ij,l,1) = 0.
  	  END IF
  	END DO
      END DO
    END DO
   
  END SUBROUTINE setuprad
  
  SUBROUTINE calcproperties
    ! **********************************************************************
    !
    !		 CALCULATE THE AEROSOL EXTINCTION CROSS SECTIONS
    !
    ! **********************************************************************
    !
    !	  Get <is_grp_ice> and radius grid from interface common block
    !	  and calculate cross-sectional area for each bin.
    !
    USE mem_aerad, ONLY: is_grp_ice_aerad,r_aerad,rcore_aerad,rup_aerad, &
         rcoreup_aerad,lunmie,lunoprt,ir_above_aerad, tabove_aerad
  
    USE mem_globrad, ONLY: ngroup,nrad,core_rad,xsecta,pi, &
         coreup_rad,i_write,i_read,nwave,rmin,rdqext, &
         qscat,qbrqs,nsol,wave,corereal,coreimag,nsolp, &
         weight,ntotal,sol,solfx,nprob,iblackbody_above, &
         t_above,ncount,nlow,plank,sbk
  
    USE mem_globaer, ONLY: r,rcore,rup,rcoreup,is_grp_ice
  
    IMPLICIT NONE
    
    INTEGER :: i
    INTEGER :: ibeyond_spectrum
    INTEGER :: ig
    INTEGER :: i_mie
    INTEGER :: irefr
    INTEGER :: j
    INTEGER :: jj
    INTEGER :: k
    INTEGER :: l
    INTEGER :: mgroup
    INTEGER :: mrad
    INTEGER :: mwave
    INTEGER :: n_thetd
    LOGICAL :: all_ok
    REAL    :: awave
    REAL    :: corerad
    REAL    :: ctbrqs
    REAL    :: ddr
    REAL    :: ddrc
    REAL    :: qextd
    REAL    :: qscatd
    REAL    :: r_real
    REAL    :: rr
    REAL    :: sum
    REAL    :: sum1
    REAL    :: sum2
    REAL    :: t1
    REAL    :: thetd(1)
    REAL    :: tmag
    REAL    :: v
    REAL    :: wvno
  
    CHARACTER(LEN=*),PARAMETER :: &
         lab355='(//,"setuprad: error in weights ",/," ' &
         //'sum of weights for solar =",1pe15.5,/,"' &
         //' sum of weights for ir = ",1pe15.5,/," ' &
         //'total sum =  ",1pe15.5)'
    
    DO ig = 1, ngroup
       DO I = 1, nrad
          xsecta(i,ig) = pi * r(i,ig)**2.
       END DO
    END DO
    !
    !	  Set <i_mie> = I_READ to WRITE the mie coefficients to a data file,
    !		      = I_WRITE to READ them
    !
    i_mie = I_READ
    !i_mie = i_WRITE
    !srf - nao precisa escrever o arquivo mie.data
    !	   i_mie = 555
    
    ! ------------------------------------------------------------------
    ! "mie.data" Parameters definitions
    mwave=50
    mrad=30
    mgroup=1
    rmin=(/0.1000000E-05/)
    !(1:30,1:1,1:50)
    rdqext(1:30,1:1,1:6)=reshape( (/                                 &
         0.1238996E-05, 0.9044850E-03, 0.1058328E-02, 0.1310231E-02, &
         0.1623058E-02, 0.2011839E-02, 0.2492962E-02, 0.3089402E-02, &
         0.3829632E-02, 0.4747202E-02, 0.5886665E-02, 0.7303554E-02, &
         0.9068833E-02, 0.1127584E-01, 0.1405085E-01, 0.1757244E-01, &
         0.2210986E-01, 0.2809986E-01, 0.3631130E-01, 0.4820201E-01, &
         0.6671599E-01, 0.9805208E-01, 0.1554057E+00, 0.2653104E+00, &
         0.4692753E+00, 0.7888026E+00, 0.1204756E+01, 0.1752262E+01, &
         0.2199486E+01, 0.2739768E+01, 0.2356571E-05, 0.8107370E-03, &
         0.9485792E-03, 0.1175354E-02, 0.1457698E-02, 0.1806319E-02, &
         0.2237227E-02, 0.2773230E-02, 0.3436975E-02, 0.4259836E-02, &
         0.5281698E-02, 0.6551373E-02, 0.8131125E-02, 0.1010238E-01, &
         0.1257303E-01, 0.1569186E-01, 0.1967621E-01, 0.2486489E-01, &
         0.3183057E-01, 0.4162014E-01, 0.5628673E-01, 0.8007501E-01, &
         0.1220185E+00, 0.2010520E+00, 0.3521606E+00, 0.6149174E+00, &
         0.9806738E+00, 0.1471468E+01, 0.1975219E+01, 0.2499236E+01, &
         0.4482199E-05, 0.6878105E-03, 0.8052588E-03, 0.9958311E-03, &
         0.1234672E-02, 0.1529180E-02, 0.1895103E-02, 0.2348306E-02, &
         0.2910469E-02, 0.3606677E-02, 0.4470860E-02, 0.5543770E-02, &
         0.6876979E-02, 0.8536757E-02, 0.1060989E-01, 0.1321133E-01, &
         0.1650270E-01, 0.2072231E-01, 0.2624893E-01, 0.3373386E-01, &
         0.4438854E-01, 0.6062073E-01, 0.8745370E-01, 0.1355955E+00, &
         0.2271680E+00, 0.4006258E+00, 0.6897491E+00, 0.1075740E+01, &
         0.1600417E+01, 0.2066453E+01, 0.8525142E-05, 0.5920068E-03, &
         0.6939178E-03, 0.8606558E-03, 0.1066172E-02, 0.1321362E-02, &
         0.1637412E-02, 0.2028739E-02, 0.2514813E-02, 0.3116922E-02, &
         0.3863455E-02, 0.4789023E-02, 0.5938904E-02, 0.7368284E-02, &
         0.9149896E-02, 0.1137746E-01, 0.1417908E-01, 0.1773632E-01, &
         0.2232307E-01, 0.2838602E-01, 0.3671290E-01, 0.4880296E-01, &
         0.6768916E-01, 0.9976597E-01, 0.1586369E+00, 0.2715197E+00, &
         0.4801895E+00, 0.8039033E+00, 0.1225015E+01, 0.1773117E+01, &
         0.1621482E-04, 0.5126877E-03, 0.5990996E-03, 0.7428692E-03, &
         0.9220786E-03, 0.1141519E-02, 0.1414444E-02, 0.1753329E-02, &
         0.2172156E-02, 0.2691440E-02, 0.3335611E-02, 0.4134388E-02, &
         0.5125767E-02, 0.6357079E-02, 0.7889097E-02, 0.9800074E-02, &
         0.1219303E-01, 0.1521039E-01, 0.1905701E-01, 0.2405036E-01, &
         0.3071964E-01, 0.4002312E-01, 0.5382273E-01, 0.7594469E-01, &
         0.1145147E+00, 0.1866735E+00, 0.3250482E+00, 0.5710782E+00, &
         0.9245852E+00, 0.1392398E+01, 0.3084059E-04, 0.4278156E-03, &
         0.4990920E-03, 0.6195636E-03, 0.7684358E-03, 0.9515302E-03, &
         0.1178236E-02, 0.1461325E-02, 0.1809853E-02, 0.2242170E-02, &
         0.2778809E-02, 0.3443943E-02, 0.4269208E-02, 0.5292946E-02, &
         0.6565067E-02, 0.8148359E-02, 0.1012359E-01, 0.1259963E-01, &
         0.1572569E-01, 0.1971983E-01, 0.2492236E-01, 0.3190922E-01, &
         0.4173374E-01, 0.5646301E-01, 0.8037232E-01, 0.1225615E+00, &
         0.2020946E+00, 0.3541155E+00, 0.6180223E+00, 0.9846209E+00/),&
         (/30,1,6/))

    rdqext(1:30,1:1,7:12)=reshape( (/&
         0.5865879E-04, 0.3386404E-03, 0.3943375E-03, 0.4902230E-03, &
         0.6069456E-03, 0.7545453E-03, 0.9325470E-03, 0.1156643E-02, &
         0.1432099E-02, 0.1774964E-02, 0.2199463E-02, 0.2725592E-02, &
         0.3378012E-02, 0.4187341E-02, 0.5191328E-02, 0.6438720E-02, &
         0.7990645E-02, 0.9926975E-02, 0.1235256E-01, 0.1541248E-01, &
         0.1931667E-01, 0.2439146E-01, 0.3118407E-01, 0.4068905E-01, &
         0.5484687E-01, 0.7765548E-01, 0.1176138E+00, 0.1926048E+00, &
         0.3362654E+00, 0.5893933E+00, 0.1115690E-03, 0.2516592E-03, &
         0.2934801E-03, 0.3651093E-03, 0.4550150E-03, 0.5632643E-03, &
         0.6970990E-03, 0.8643277E-03, 0.1070622E-02, 0.1324547E-02, &
         0.1642295E-02, 0.2034783E-02, 0.2521489E-02, 0.3124935E-02, &
         0.3872860E-02, 0.4801260E-02, 0.5953945E-02, 0.7387332E-02, &
         0.9173343E-02, 0.1140706E-01, 0.1421631E-01, 0.1778384E-01, &
         0.2238511E-01, 0.2846927E-01, 0.3683005E-01, 0.4897850E-01, &
         0.6797409E-01, 0.1002687E+00, 0.1595861E+00, 0.2733425E+00, &
         0.2122042E-03, 0.1571589E-03, 0.1840258E-03, 0.2313218E-03, &
         0.2901404E-03, 0.3512961E-03, 0.4390195E-03, 0.5448884E-03, &
         0.6724349E-03, 0.8375491E-03, 0.1036835E-02, 0.1285574E-02, &
         0.1592136E-02, 0.1972729E-02, 0.2444628E-02, 0.3029664E-02, &
         0.3754555E-02, 0.4654815E-02, 0.5771693E-02, 0.7160614E-02, &
         0.8890757E-02, 0.1105248E-01, 0.1376909E-01, 0.1721303E-01, &
         0.2164283E-01, 0.2747480E-01, 0.3543673E-01, 0.4689975E-01, &
         0.6461873E-01, 0.9437821E-01, 0.4036125E-03, 0.4093335E-01, &
         0.5522372E-01, 0.7828710E-01, 0.1187611E+00, 0.1948035E+00, &
         0.3404121E+00, 0.5961001E+00, 0.9566938E+00, 0.1437752E+01, &
         0.1951394E+01, 0.2464083E+01, 0.2949664E+01, 0.3406605E+01, &
         0.3629449E+01, 0.3356592E+01, 0.2565616E+01, 0.1934376E+01, &
         0.2454620E+01, 0.2495687E+01, 0.2086732E+01, 0.2298010E+01, &
         0.2203053E+01, 0.2162956E+01, 0.2135190E+01, 0.2162642E+01, &
         0.2130731E+01, 0.2120723E+01, 0.2111154E+01, 0.2073417E+01, &
         0.3014119E-01, 0.3649874E-01, 0.4848227E-01, 0.6716944E-01, &
         0.9884995E-01, 0.1569090E+00, 0.2682001E+00, 0.4743640E+00, &
         0.7958630E+00, 0.1214206E+01, 0.1762091E+01, 0.2209994E+01, &
         0.2747663E+01, 0.3190050E+01, 0.3513948E+01, 0.3566287E+01, &
         0.2919923E+01, 0.2062144E+01, 0.2168095E+01, 0.2654343E+01, &
         0.2093068E+01, 0.2346398E+01, 0.2183643E+01, 0.2237450E+01, &
         0.2201363E+01, 0.2169828E+01, 0.2122523E+01, 0.2131855E+01, &
         0.2098403E+01, 0.2083777E+01, 0.2825850E-01, 0.3407567E-01, &
         0.4489003E-01, 0.6141419E-01, 0.8881936E-01, 0.1381295E+00, &
         0.2320493E+00, 0.4095595E+00, 0.7030742E+00, 0.1092780E+01, &
         0.1622257E+01, 0.2083007E+01, 0.2637337E+01, 0.3099045E+01, &
         0.3466475E+01, 0.3576084E+01, 0.3138925E+01, 0.2289058E+01, &
         0.2024652E+01, 0.2619606E+01, 0.2271519E+01, 0.2340897E+01, &
         0.2215942E+01, 0.2211942E+01, 0.2123458E+01, 0.2147658E+01, &
         0.2145206E+01, 0.2089916E+01, 0.2105316E+01, 0.2094704E+01/),&
         (/30,1,6/))

    rdqext(1:30,1:1,13:18)=reshape( (/&
         0.2595173E-01, 0.3114211E-01, 0.4062882E-01, 0.5475407E-01, &
         0.7750003E-01, 0.1173315E+00, 0.1920643E+00, 0.3352450E+00, &
         0.5877376E+00, 0.9459963E+00, 0.1422653E+01, 0.1940509E+01, &
         0.2447865E+01, 0.2935348E+01, 0.3393993E+01, 0.3634520E+01, &
         0.3377689E+01, 0.2585614E+01, 0.1928602E+01, 0.2453908E+01, &
         0.2537595E+01, 0.2120584E+01, 0.2355100E+01, 0.2236805E+01, &
         0.2206528E+01, 0.2137048E+01, 0.2110582E+01, 0.2129888E+01, &
         0.2112728E+01, 0.2111520E+01, 0.2456848E-01, 0.2940105E-01, &
         0.3814531E-01, 0.5096052E-01, 0.7121034E-01, 0.1060126E+00, &
         0.1704660E+00, 0.2941915E+00, 0.5193896E+00, 0.8569373E+00, &
         0.1297529E+01, 0.1841395E+01, 0.2305595E+01, 0.2818348E+01, &
         0.3263292E+01, 0.3597307E+01, 0.3494461E+01, 0.2745187E+01, &
         0.2013707E+01, 0.2287942E+01, 0.2655969E+01, 0.2131976E+01, &
         0.2384275E+01, 0.2205176E+01, 0.2206340E+01, 0.2215608E+01, &
         0.2134920E+01, 0.2115479E+01, 0.2095288E+01, 0.2098952E+01, &
         0.2235298E-01, 0.2663950E-01, 0.3427472E-01, 0.4518270E-01, &
         0.6187829E-01, 0.8962012E-01, 0.1396183E+00, 0.2349173E+00, &
         0.4147888E+00, 0.7108101E+00, 0.1102706E+01, 0.1634750E+01, &
         0.2092761E+01, 0.2647524E+01, 0.3107261E+01, 0.3467156E+01, &
         0.3579113E+01, 0.3125170E+01, 0.2247379E+01, 0.2026859E+01, &
         0.2657282E+01, 0.2265156E+01, 0.2332558E+01, 0.2186288E+01, &
         0.2227778E+01, 0.2146957E+01, 0.2169617E+01, 0.2140218E+01, &
         0.2092722E+01, 0.2088475E+01, 0.1921333E-01, 0.2277928E-01, &
         0.2899945E-01, 0.3757729E-01, 0.5010217E-01, 0.6980440E-01, &
         0.1035102E+00, 0.1657181E+00, 0.2851043E+00, 0.5038010E+00, &
         0.8360577E+00, 0.1268737E+01, 0.1815424E+01, 0.2272147E+01, &
         0.2793396E+01, 0.3235432E+01, 0.3570348E+01, 0.3524326E+01, &
         0.2809205E+01, 0.2151345E+01, 0.2281744E+01, 0.2666198E+01, &
         0.2126539E+01, 0.2388358E+01, 0.2219515E+01, 0.2215051E+01, &
         0.2223578E+01, 0.2160036E+01, 0.2126196E+01, 0.2106740E+01, &
         0.1828648E-01, 0.2165078E-01, 0.2748548E-01, 0.3545161E-01, &
         0.4692184E-01, 0.6465408E-01, 0.9444001E-01, 0.1486191E+00, &
         0.2522533E+00, 0.4460807E+00, 0.7561581E+00, 0.1161547E+01, &
         0.1705030E+01, 0.2152755E+01, 0.2702824E+01, 0.3151649E+01, &
         0.3478999E+01, 0.3582580E+01, 0.3043263E+01, 0.2252439E+01, &
         0.2120131E+01, 0.2719469E+01, 0.2260474E+01, 0.2383116E+01, &
         0.2185176E+01, 0.2205244E+01, 0.2167593E+01, 0.2190734E+01, &
         0.2141221E+01, 0.2117435E+01, 0.1705720E-01, 0.2016090E-01, &
         0.2550504E-01, 0.3270831E-01, 0.4289185E-01, 0.5826783E-01, &
         0.8343055E-01, 0.1281675E+00, 0.2128689E+00, 0.3742239E+00, &
         0.6495301E+00, 0.1024606E+01, 0.1532292E+01, 0.2017566E+01, &
         0.2558653E+01, 0.3033171E+01, 0.3455450E+01, 0.3578496E+01, &
         0.3254357E+01, 0.2367896E+01, 0.1937363E+01, 0.2580973E+01, &
         0.2382489E+01, 0.2224710E+01, 0.2240751E+01, 0.2242849E+01, &
         0.2187788E+01, 0.2117765E+01, 0.2155853E+01, 0.2127444E+01/),&
         (/30,1,6/))

    rdqext(1:30,1:1,19:24)=reshape( (/&
         0.1661168E-01, 0.1962307E-01, 0.2479483E-01, 0.3173472E-01, &
         0.4148188E-01, 0.5607228E-01, 0.7971362E-01, 0.1213590E+00, &
         0.1997862E+00, 0.3497850E+00, 0.6111342E+00, 0.9758610E+00, &
         0.1464721E+01, 0.1970491E+01, 0.2492324E+01, 0.2974746E+01, &
         0.3425845E+01, 0.3616955E+01, 0.3324818E+01, 0.2534858E+01, &
         0.2057498E+01, 0.2525669E+01, 0.2476210E+01, 0.2140753E+01, &
         0.2303478E+01, 0.2236078E+01, 0.2172018E+01, 0.2146353E+01, &
         0.2130648E+01, 0.2142246E+01, 0.1598668E-01, 0.1886983E-01, &
         0.2380455E-01, 0.3038575E-01, 0.3954578E-01, 0.5309152E-01, &
         0.7472860E-01, 0.1123202E+00, 0.1824803E+00, 0.3170920E+00, &
         0.5579348E+00, 0.9075783E+00, 0.1368394E+01, 0.1899847E+01, &
         0.2387506E+01, 0.2883463E+01, 0.3339896E+01, 0.3636674E+01, &
         0.3422941E+01, 0.2665677E+01, 0.2101645E+01, 0.2395268E+01, &
         0.2593362E+01, 0.2087997E+01, 0.2348550E+01, 0.2197394E+01, &
         0.2205405E+01, 0.2186869E+01, 0.2118561E+01, 0.2120739E+01, &
         0.1447862E-01, 0.1705971E-01, 0.2144386E-01, 0.2720915E-01, &
         0.3506626E-01, 0.4635068E-01, 0.6373934E-01, 0.9284600E-01, &
         0.1456350E+00, 0.2465072E+00, 0.4357705E+00, 0.7413899E+00, &
         0.1142247E+01, 0.1682719E+01, 0.2132625E+01, 0.2685511E+01, &
         0.3137660E+01, 0.3472792E+01, 0.3583071E+01, 0.3077484E+01, &
         0.2215237E+01, 0.2114642E+01, 0.2624185E+01, 0.2220664E+01, &
         0.2292433E+01, 0.2218276E+01, 0.2214263E+01, 0.2160881E+01, &
         0.2179611E+01, 0.2127588E+01, 0.1328048E-01, 0.1562869E-01, &
         0.1959484E-01, 0.2475764E-01, 0.3168387E-01, 0.4140854E-01, &
         0.5595865E-01, 0.7952227E-01, 0.1210101E+00, 0.1991165E+00, &
         0.3485271E+00, 0.6091268E+00, 0.9733056E+00, 0.1461134E+01, &
         0.1967971E+01, 0.2488624E+01, 0.2971455E+01, 0.3423538E+01, &
         0.3618834E+01, 0.3328263E+01, 0.2539237E+01, 0.2039369E+01, &
         0.2503954E+01, 0.2453909E+01, 0.2115271E+01, 0.2290610E+01, &
         0.2230195E+01, 0.2159497E+01, 0.2120082E+01, 0.2116015E+01, &
         0.1200182E-01, 0.1410724E-01, 0.1764445E-01, 0.2220361E-01, &
         0.2822565E-01, 0.3648769E-01, 0.4846574E-01, 0.6714261E-01, &
         0.9880274E-01, 0.1568201E+00, 0.2680290E+00, 0.4740632E+00, &
         0.7954467E+00, 0.1213648E+01, 0.1761516E+01, 0.2209371E+01, &
         0.2747197E+01, 0.3189625E+01, 0.3513435E+01, 0.3566643E+01, &
         0.2920861E+01, 0.2063354E+01, 0.2168825E+01, 0.2651909E+01, &
         0.2093758E+01, 0.2346431E+01, 0.2182174E+01, 0.2229954E+01, &
         0.2199746E+01, 0.2167265E+01, 0.1181864E-01, 0.1388969E-01, &
         0.1736677E-01, 0.2184242E-01, 0.2774175E-01, 0.3580971E-01, &
         0.4745414E-01, 0.6550952E-01, 0.9593535E-01, 0.1514246E+00, &
         0.2576532E+00, 0.4557129E+00, 0.7698092E+00, 0.1179520E+01, &
         0.1725135E+01, 0.2171919E+01, 0.2718426E+01, 0.3164545E+01, &
         0.3487779E+01, 0.3580781E+01, 0.2997886E+01, 0.2219576E+01, &
         0.2137670E+01, 0.2651932E+01, 0.2223345E+01, 0.2425306E+01, &
         0.2227556E+01, 0.2204142E+01, 0.2194215E+01, 0.2185364E+01/),&
         (/30,1,6/))

    rdqext(1:30,1:1,25:30)=reshape( (/&
         0.1146874E-01, 0.1347448E-01, 0.1683753E-01, 0.2115573E-01, &
         0.2682521E-01, 0.3453235E-01, 0.4556208E-01, 0.6248134E-01, &
         0.9066294E-01, 0.1415598E+00, 0.2386575E+00, 0.4215866E+00, &
         0.7207971E+00, 0.1115563E+01, 0.1650671E+01, 0.2105534E+01, &
         0.2660306E+01, 0.3117507E+01, 0.3468289E+01, 0.3581810E+01, &
         0.3107697E+01, 0.2204360E+01, 0.2083869E+01, 0.2692604E+01, &
         0.2246163E+01, 0.2310285E+01, 0.2205017E+01, 0.2212293E+01, &
         0.2154623E+01, 0.2157726E+01, 0.1036541E-01, 0.1216810E-01, &
         0.1517874E-01, 0.1901641E-01, 0.2399688E-01, 0.3064706E-01, &
         0.3991923E-01, 0.5366334E-01, 0.7567919E-01, 0.1140351E+00, &
         0.1857564E+00, 0.3233098E+00, 0.5682176E+00, 0.9208928E+00, &
         0.1387183E+01, 0.1914258E+01, 0.2408734E+01, 0.2901393E+01, &
         0.3359867E+01, 0.3639189E+01, 0.3413174E+01, 0.2648792E+01, &
         0.2071824E+01, 0.2370969E+01, 0.2552859E+01, 0.2186358E+01, &
         0.2357531E+01, 0.2219193E+01, 0.2180621E+01, 0.2176074E+01, &
         0.1099698E-01, 0.1291546E-01, 0.1612651E-01, 0.2023638E-01, &
         0.2560486E-01, 0.3284562E-01, 0.4309161E-01, 0.5858051E-01, &
         0.8396278E-01, 0.1291468E+00, 0.2147529E+00, 0.3777217E+00, &
         0.6549321E+00, 0.1031459E+01, 0.1541618E+01, 0.2024099E+01, &
         0.2567300E+01, 0.3040632E+01, 0.3457704E+01, 0.3574554E+01, &
         0.3238655E+01, 0.2353804E+01, 0.1952011E+01, 0.2580399E+01, &
         0.2390371E+01, 0.2200474E+01, 0.2244411E+01, 0.2220684E+01, &
         0.2163651E+01, 0.2148918E+01, 0.9705304E-02, 0.1138785E-01, &
         0.1419229E-01, 0.1775317E-01, 0.2234502E-01, 0.2841550E-01, &
         0.3675437E-01, 0.4886510E-01, 0.6778988E-01, 0.9994365E-01, &
         0.1589722E+00, 0.2721637E+00, 0.4813170E+00, 0.8054546E+00, &
         0.1227107E+01, 0.1775224E+01, 0.2224483E+01, 0.2758411E+01, &
         0.3200041E+01, 0.3526455E+01, 0.3557312E+01, 0.2899145E+01, &
         0.2049803E+01, 0.2180376E+01, 0.2742225E+01, 0.2126863E+01, &
         0.2353465E+01, 0.2187709E+01, 0.2224508E+01, 0.2222815E+01, &
         0.9146382E-02, 0.1072809E-01, 0.1336029E-01, 0.1669216E-01, &
         0.2096744E-01, 0.2657457E-01, 0.3418472E-01, 0.4505027E-01, &
         0.6166824E-01, 0.8925748E-01, 0.1389438E+00, 0.2336179E+00, &
         0.4124213E+00, 0.7073137E+00, 0.1098216E+01, 0.1629121E+01, &
         0.2088338E+01, 0.2642952E+01, 0.3103581E+01, 0.3466837E+01, &
         0.3577819E+01, 0.3131474E+01, 0.2266839E+01, 0.2022272E+01, &
         0.2655272E+01, 0.2286961E+01, 0.2383289E+01, 0.2214466E+01, &
         0.2205394E+01, 0.2140238E+01, 0.8506889E-02, 0.9974504E-02, &
         0.1241214E-01, 0.1548789E-01, 0.1941376E-01, 0.2451917E-01, &
         0.3135827E-01, 0.4093947E-01, 0.5523325E-01, 0.7830309E-01, &
         0.1187902E+00, 0.1948592E+00, 0.3405170E+00, 0.5962692E+00, &
         0.9569100E+00, 0.1438057E+01, 0.1951612E+01, 0.2464408E+01, &
         0.2949951E+01, 0.3406850E+01, 0.3629331E+01, 0.3356170E+01, &
         0.2565260E+01, 0.1935007E+01, 0.2454862E+01, 0.2496692E+01, &
         0.2086216E+01, 0.2314208E+01, 0.2219513E+01, 0.2162504E+01/),&
         (/30,1,6/))

    rdqext(1:30,1:1,31:36)=reshape( (/&
         0.7831769E-02, 0.9179345E-02, 0.1141457E-01, 0.1422596E-01, &
         0.1779609E-01, 0.2240106E-01, 0.2849078E-01, 0.3686021E-01, &
         0.4902377E-01, 0.6804754E-01, 0.1003985E+00, 0.1598312E+00, &
         0.2738129E+00, 0.4842010E+00, 0.8094147E+00, 0.1232453E+01, &
         0.1780571E+01, 0.2230534E+01, 0.2762867E+01, 0.3204298E+01, &
         0.3531933E+01, 0.3553329E+01, 0.2890255E+01, 0.2053158E+01, &
         0.2194848E+01, 0.2679376E+01, 0.2134710E+01, 0.2353978E+01, &
         0.2197762E+01, 0.2211876E+01, 0.7078757E-02, 0.8293798E-02, &
         0.1030588E-01, 0.1282875E-01, 0.1601650E-01, 0.2009444E-01, &
         0.2541709E-01, 0.3258747E-01, 0.4271624E-01, 0.5799330E-01, &
         0.8296370E-01, 0.1273095E+00, 0.2112185E+00, 0.3711553E+00, &
         0.6447718E+00, 0.1018571E+01, 0.1524037E+01, 0.2011806E+01, &
         0.2550897E+01, 0.3026439E+01, 0.3453130E+01, 0.3582561E+01, &
         0.3267390E+01, 0.2385726E+01, 0.1939468E+01, 0.2592190E+01, &
         0.2432119E+01, 0.2266550E+01, 0.2245551E+01, 0.2216399E+01, &
         0.5893333E-02, 0.6901805E-02, 0.8567906E-02, 0.1064833E-01, &
         0.1325984E-01, 0.1656435E-01, 0.2080209E-01, 0.2635475E-01, &
         0.3388031E-01, 0.4460318E-01, 0.6096004E-01, 0.8803720E-01, &
         0.1366774E+00, 0.2292520E+00, 0.4044449E+00, 0.6954624E+00, &
         0.1083038E+01, 0.1609830E+01, 0.2073517E+01, 0.2627049E+01, &
         0.3090688E+01, 0.3465830E+01, 0.3572824E+01, 0.3151445E+01, &
         0.2316983E+01, 0.2038591E+01, 0.2610828E+01, 0.2285346E+01, &
         0.2321426E+01, 0.2202826E+01, 0.5025937E-02, 0.5884357E-02, &
         0.7300625E-02, 0.9065186E-02, 0.1127127E-01, 0.1404497E-01, &
         0.1756500E-01, 0.2210026E-01, 0.2808706E-01, 0.3629318E-01, &
         0.4817502E-01, 0.6667232E-01, 0.9797529E-01, 0.1552611E+00, &
         0.2650325E+00, 0.4687851E+00, 0.7881204E+00, 0.1203845E+01, &
         0.1751305E+01, 0.2198479E+01, 0.2739005E+01, 0.3182265E+01, &
         0.3504886E+01, 0.3572367E+01, 0.2939144E+01, 0.2093931E+01, &
         0.2184590E+01, 0.2653599E+01, 0.2120970E+01, 0.2362548E+01, &
         0.4433422E-02, 0.5189727E-02, 0.6437071E-02, 0.7988810E-02, &
         0.9924201E-02, 0.1234918E-01, 0.1540817E-01, 0.1931125E-01, &
         0.2438432E-01, 0.3117430E-01, 0.4067514E-01, 0.5482540E-01, &
         0.7761944E-01, 0.1175483E+00, 0.1924794E+00, 0.3360287E+00, &
         0.5890094E+00, 0.9476252E+00, 0.1424953E+01, 0.1942187E+01, &
         0.2450353E+01, 0.2937539E+01, 0.3395996E+01, 0.3633851E+01, &
         0.3374526E+01, 0.2582175E+01, 0.1926354E+01, 0.2447965E+01, &
         0.2536666E+01, 0.2114620E+01, 0.3825215E-02, 0.4476689E-02, &
         0.5550712E-02, 0.6885700E-02, 0.8548012E-02, 0.1062361E-01, &
         0.1322856E-01, 0.1652461E-01, 0.2075062E-01, 0.2628650E-01, &
         0.3378580E-01, 0.4446473E-01, 0.6074107E-01, 0.8766055E-01, &
         0.1359789E+00, 0.2279065E+00, 0.4019802E+00, 0.6917781E+00, &
         0.1078330E+01, 0.1603768E+01, 0.2068956E+01, 0.2621971E+01, &
         0.3086536E+01, 0.3465486E+01, 0.3571346E+01, 0.3157148E+01, &
         0.2325003E+01, 0.2047926E+01, 0.2618621E+01, 0.2334253E+01/),&
         (/30,1,6/))

    rdqext(1:30,1:1,37:42)=reshape( (/&
         0.3549440E-02, 0.4153976E-02, 0.5149985E-02, 0.6387233E-02, &
         0.7926833E-02, 0.9847006E-02, 0.1225196E-01, 0.1528491E-01, &
         0.1915287E-01, 0.2417619E-01, 0.3089087E-01, 0.4026837E-01, &
         0.5419936E-01, 0.7657288E-01, 0.1156511E+00, 0.1888472E+00, &
         0.3291641E+00, 0.5778280E+00, 0.9332777E+00, 0.1404682E+01, &
         0.1927353E+01, 0.2428208E+01, 0.2918166E+01, 0.3377411E+01, &
         0.3638303E+01, 0.3399381E+01, 0.2617698E+01, 0.1981718E+01, &
         0.2419561E+01, 0.2593332E+01, 0.3034851E-02, 0.3551537E-02, &
         0.4402093E-02, 0.5458310E-02, 0.6770817E-02, 0.8404701E-02, &
         0.1044440E-01, 0.1300307E-01, 0.1623797E-01, 0.2038019E-01, &
         0.2579521E-01, 0.3310764E-01, 0.4347322E-01, 0.5917915E-01, &
         0.8498399E-01, 0.1310288E+00, 0.2183748E+00, 0.3844299E+00, &
         0.6652285E+00, 0.1044529E+01, 0.1559250E+01, 0.2036551E+01, &
         0.2583322E+01, 0.3054318E+01, 0.3461036E+01, 0.3569443E+01, &
         0.3209359E+01, 0.2341575E+01, 0.2014180E+01, 0.2588365E+01, &
         0.2387623E-02, 0.2794191E-02, 0.3462401E-02, 0.4291734E-02, &
         0.5320973E-02, 0.6599757E-02, 0.8191557E-02, 0.1017778E-01, &
         0.1266782E-01, 0.1581217E-01, 0.1983106E-01, 0.2506921E-01, &
         0.3211028E-01, 0.4202449E-01, 0.5691476E-01, 0.8113549E-01, &
         0.1239570E+00, 0.2047747E+00, 0.3591337E+00, 0.6259582E+00, &
         0.9947014E+00, 0.1491055E+01, 0.1988876E+01, 0.2518942E+01, &
         0.2998372E+01, 0.3440423E+01, 0.3601781E+01, 0.3304095E+01, &
         0.2485376E+01, 0.2048215E+01, 0.2057811E-02, 0.2408417E-02, &
         0.2984623E-02, 0.3698977E-02, 0.4585142E-02, 0.5685347E-02, &
         0.7053476E-02, 0.8756823E-02, 0.1088493E-01, 0.1355791E-01, &
         0.1694377E-01, 0.2129343E-01, 0.2700860E-01, 0.3478724E-01, &
         0.4593816E-01, 0.6308043E-01, 0.9170132E-01, 0.1434965E+00, &
         0.2423881E+00, 0.4283417E+00, 0.7306451E+00, 0.1128293E+01, &
         0.1666130E+01, 0.2118356E+01, 0.2672537E+01, 0.3127278E+01, &
         0.3469966E+01, 0.3582943E+01, 0.3093040E+01, 0.2194497E+01, &
         0.1826921E-02, 0.2137052E-02, 0.2648383E-02, 0.3281895E-02, &
         0.4067792E-02, 0.5043134E-02, 0.6254602E-02, 0.7761572E-02, &
         0.9640561E-02, 0.1199300E-01, 0.1495711E-01, 0.1873215E-01, &
         0.2362413E-01, 0.3014096E-01, 0.3919665E-01, 0.5255830E-01, &
         0.7384456E-01, 0.1107295E+00, 0.1794446E+00, 0.3113199E+00, &
         0.5483193E+00, 0.8950609E+00, 0.1350772E+01, 0.1885957E+01, &
         0.2367349E+01, 0.2866825E+01, 0.3320561E+01, 0.3631001E+01, &
         0.3433996E+01, 0.2677806E+01, 0.1663156E-02, 0.1945878E-02, &
         0.2411773E-02, 0.2989146E-02, 0.3704706E-02, 0.4592277E-02, &
         0.5694410E-02, 0.7064310E-02, 0.8770572E-02, 0.1090218E-01, &
         0.1357968E-01, 0.1697154E-01, 0.2132941E-01, 0.2705653E-01, &
         0.3485396E-01, 0.4603667E-01, 0.6323758E-01, 0.9197413E-01, &
         0.1440057E+00, 0.2433691E+00, 0.4301138E+00, 0.7332162E+00, &
         0.1131626E+01, 0.1670127E+01, 0.2121743E+01, 0.2675674E+01, &
         0.3129785E+01, 0.3470570E+01, 0.3583046E+01, 0.3089532E+01/),&
         (/30,1,6/))

    rdqext(1:30,1:1,43:48)=reshape( (/&
         0.1567061E-02, 0.1833278E-02, 0.2271720E-02, 0.2815047E-02, &
         0.3488794E-02, 0.4324771E-02, 0.5361873E-02, 0.6650588E-02, &
         0.8254934E-02, 0.1025699E-01, 0.1276743E-01, 0.1593864E-01, &
         0.1999407E-01, 0.2528441E-01, 0.3240535E-01, 0.4245194E-01, &
         0.5758059E-01, 0.8226308E-01, 0.1260231E+00, 0.2087452E+00, &
         0.3665486E+00, 0.6375951E+00, 0.1009468E+01, 0.1511518E+01, &
         0.2003096E+01, 0.2538950E+01, 0.3016001E+01, 0.3448988E+01, &
         0.3589492E+01, 0.3284389E+01, 0.1488439E-02, 0.1741797E-02, &
         0.2157624E-02, 0.2673772E-02, 0.3314134E-02, 0.4107578E-02, &
         0.5092440E-02, 0.6315934E-02, 0.7837981E-02, 0.9735893E-02, &
         0.1211273E-01, 0.1510853E-01, 0.1892637E-01, 0.2387865E-01, &
         0.3048642E-01, 0.3968957E-01, 0.5331150E-01, 0.7509397E-01, &
         0.1129788E+00, 0.1837380E+00, 0.3194806E+00, 0.5618942E+00, &
         0.9127132E+00, 0.1375636E+01, 0.1905448E+01, 0.2395723E+01, &
         0.2890356E+01, 0.3347702E+01, 0.3638081E+01, 0.3419334E+01, &
         0.1409372E-02, 0.1649454E-02, 0.2044782E-02, 0.2532639E-02, &
         0.3138920E-02, 0.3890525E-02, 0.4823158E-02, 0.5981331E-02, &
         0.7421338E-02, 0.9215627E-02, 0.1146004E-01, 0.1428342E-01, &
         0.1786957E-01, 0.2249676E-01, 0.2861935E-01, 0.3704118E-01, &
         0.4929540E-01, 0.6848914E-01, 0.1011789E+00, 0.1613057E+00, &
         0.2766435E+00, 0.4891382E+00, 0.8161697E+00, 0.1241602E+01, &
         0.1789592E+01, 0.2240942E+01, 0.2770504E+01, 0.3211761E+01, &
         0.3541587E+01, 0.3546342E+01, 0.1330705E-02, 0.1558003E-02, &
         0.1930830E-02, 0.2391852E-02, 0.2964355E-02, 0.3673612E-02, &
         0.4553958E-02, 0.5646588E-02, 0.7004765E-02, 0.8696110E-02, &
         0.1080916E-01, 0.1346239E-01, 0.1682217E-01, 0.2113581E-01, &
         0.2679852E-01, 0.3449527E-01, 0.4550745E-01, 0.6239447E-01, &
         0.9051254E-01, 0.1412796E+00, 0.2381176E+00, 0.4206070E+00, &
         0.7193625E+00, 0.1113713E+01, 0.1648399E+01, 0.2103685E+01, &
         0.2658494E+01, 0.3116058E+01, 0.3468100E+01, 0.3581523E+01, &
         0.1251933E-02, 0.1466046E-02, 0.1816172E-02, 0.2251175E-02, &
         0.2788802E-02, 0.3456600E-02, 0.4284468E-02, 0.5312049E-02, &
         0.6588972E-02, 0.8177907E-02, 0.1016081E-01, 0.1264632E-01, &
         0.1578495E-01, 0.1979593E-01, 0.2502291E-01, 0.3204691E-01, &
         0.4193268E-01, 0.5677203E-01, 0.8089415E-01, 0.1235155E+00, &
         0.2039266E+00, 0.3575470E+00, 0.6234539E+00, 0.9915216E+00, &
         0.1486626E+01, 0.1985795E+01, 0.2514534E+01, 0.2994468E+01, &
         0.3438259E+01, 0.3604452E+01, 0.1174132E-02, 0.1373839E-02, &
         0.1703073E-02, 0.2110218E-02, 0.2614515E-02, 0.3239858E-02, &
         0.4015790E-02, 0.4978259E-02, 0.6173875E-02, 0.7660876E-02, &
         0.9514965E-02, 0.1183530E-01, 0.1475762E-01, 0.1847636E-01, &
         0.2328914E-01, 0.2968762E-01, 0.3855173E-01, 0.5157678E-01, &
         0.7222367E-01, 0.1078229E+00, 0.1739076E+00, 0.3007662E+00, &
         0.5305643E+00, 0.8717529E+00, 0.1318129E+01, 0.1859136E+01, &
         0.2329573E+01, 0.2836743E+01, 0.3284825E+01, 0.3613000E+01/),&
         (/30,1,6/))

    rdqext(1:30,1:1,49:50)=reshape( (/&
         0.1095844E-02, 0.1281458E-02, 0.1588486E-02, 0.1968793E-02, &
         0.2439237E-02, 0.3023065E-02, 0.3747140E-02, 0.4644954E-02, &
         0.5759797E-02, 0.7145516E-02, 0.8871601E-02, 0.1102877E-01, &
         0.1373922E-01, 0.1717497E-01, 0.2159321E-01, 0.2740864E-01, &
         0.3534436E-01, 0.4676276E-01, 0.6439901E-01, 0.9399497E-01, &
         0.1477852E+00, 0.2506478E+00, 0.4432062E+00, 0.7520572E+00, &
         0.1156173E+01, 0.1698893E+01, 0.2147104E+01, 0.2698069E+01, &
         0.3147782E+01, 0.3476962E+01, 0.9974022E-03, 0.1165836E-02, &
         0.1443403E-02, 0.1788278E-02, 0.2216780E-02, 0.2746567E-02, &
         0.3404155E-02, 0.4219844E-02, 0.5231818E-02, 0.6489299E-02, &
         0.8053846E-02, 0.1000574E-01, 0.1245149E-01, 0.1553765E-01, &
         0.1947781E-01, 0.2460348E-01, 0.3147328E-01, 0.4110512E-01, &
         0.5548911E-01, 0.7873266E-01, 0.1195716E+00, 0.1963573E+00, &
         0.3433384E+00, 0.6008126E+00, 0.9627100E+00, 0.1446232E+01, &
         0.1957438E+01, 0.2473066E+01, 0.2957629E+01, 0.3413126E+01/),&
         (/30,1,2/))

    qscat(1:30,1:1,1:6)=reshape( (/&
         0.1535112E-05, 0.8998886E-09, 0.1729664E-08, 0.4070829E-08, &
         0.9589681E-08, 0.2261375E-07, 0.5328723E-07, 0.1255802E-06, &
         0.2960315E-06, 0.6976022E-06, 0.1644110E-05, 0.3875383E-05, &
         0.9135672E-05, 0.2154004E-04, 0.5080069E-04, 0.1198578E-03, &
         0.2829565E-03, 0.6685623E-03, 0.1581542E-02, 0.3746952E-02, &
         0.8891535E-02, 0.2111296E-01, 0.4993847E-01, 0.1159231E+00, &
         0.2536170E+00, 0.4850497E+00, 0.7988749E+00, 0.1241824E+01, &
         0.1650353E+01, 0.2180239E+01, 0.2919784E-05, 0.5829020E-09, &
         0.1120681E-08, 0.2640323E-08, 0.6234001E-08, 0.1468913E-07, &
         0.3459422E-07, 0.8157038E-07, 0.1922175E-06, 0.4529455E-06, &
         0.1067609E-05, 0.2516517E-05, 0.5931983E-05, 0.1398515E-04, &
         0.3297836E-04, 0.7779063E-04, 0.1835863E-03, 0.4335718E-03, &
         0.1024991E-02, 0.2426435E-02, 0.5753419E-02, 0.1366061E-01, &
         0.3240611E-01, 0.7616633E-01, 0.1729649E+00, 0.3576158E+00, &
         0.6282287E+00, 0.1009327E+01, 0.1438690E+01, 0.1939538E+01, &
         0.5553430E-05, 0.3008072E-09, 0.5785822E-09, 0.1359781E-08, &
         0.3207448E-08, 0.7555678E-08, 0.1780665E-07, 0.4196275E-07, &
         0.9891541E-07, 0.2330545E-06, 0.5492852E-06, 0.1294664E-05, &
         0.3051624E-05, 0.7193325E-05, 0.1696003E-04, 0.3999544E-04, &
         0.9435236E-04, 0.2227023E-03, 0.5260573E-03, 0.1243972E-02, &
         0.2945860E-02, 0.6987582E-02, 0.1659290E-01, 0.3932323E-01, &
         0.9199918E-01, 0.2059242E+00, 0.4121182E+00, 0.7000768E+00, &
         0.1114519E+01, 0.1523778E+01, 0.1056262E-04, 0.1667335E-09, &
         0.3215872E-09, 0.7584442E-09, 0.1788107E-08, 0.4214065E-08, &
         0.9933288E-08, 0.2340783E-07, 0.5519188E-07, 0.1300859E-06, &
         0.3066088E-06, 0.7225132E-06, 0.1702958E-05, 0.4013857E-05, &
         0.9462706E-05, 0.2231118E-04, 0.5262020E-04, 0.1241540E-03, &
         0.2931054E-03, 0.6925743E-03, 0.1638427E-02, 0.3881993E-02, &
         0.9212516E-02, 0.2187444E-01, 0.5172065E-01, 0.1198937E+00, &
         0.2612906E+00, 0.4962334E+00, 0.8145505E+00, 0.1259660E+01, &
         0.2009011E-04, 0.9312614E-10, 0.1788817E-09, 0.4216057E-09, &
         0.9970077E-09, 0.2345042E-08, 0.5529286E-08, 0.1303760E-07, &
         0.3071025E-07, 0.7236788E-07, 0.1705414E-06, 0.4019309E-06, &
         0.9472516E-06, 0.2232579E-05, 0.5262621E-05, 0.1240705E-04, &
         0.2925507E-04, 0.6900538E-04, 0.1628382E-03, 0.3845306E-03, &
         0.9089075E-03, 0.2151185E-02, 0.5099573E-02, 0.1210664E-01, &
         0.2873163E-01, 0.6768199E-01, 0.1548432E+00, 0.3259844E+00, &
         0.5861368E+00, 0.9459786E+00, 0.3821139E-04, 0.4500313E-10, &
         0.8613462E-10, 0.2036637E-09, 0.4809451E-09, 0.1131263E-08, &
         0.2664907E-08, 0.6288710E-08, 0.1480967E-07, 0.3488467E-07, &
         0.8223534E-07, 0.1937976E-06, 0.4568121E-06, 0.1076551E-05, &
         0.2537450E-05, 0.5981515E-05, 0.1410117E-04, 0.3325193E-04, &
         0.7843683E-04, 0.1851134E-03, 0.4371810E-03, 0.1033536E-02, &
         0.2446699E-02, 0.5801552E-02, 0.1377500E-01, 0.3267637E-01, &
         0.7678847E-01, 0.1742809E+00, 0.3598650E+00, 0.6311989E+00/),&
         (/30,1,6/))

    qscat(1:30,1:1,7:12)=reshape( (/&
         0.7267806E-04, 0.1765787E-10, 0.3378414E-10, 0.7984936E-10, &
         0.1884185E-09, 0.4455544E-09, 0.1046103E-08, 0.2469826E-08, &
         0.5812942E-08, 0.1370527E-07, 0.3230256E-07, 0.7611979E-07, &
         0.1794039E-06, 0.4228477E-06, 0.9965112E-06, 0.2348792E-05, &
         0.5536366E-05, 0.1305259E-04, 0.3077852E-04, 0.7260025E-04, &
         0.1713272E-03, 0.4045929E-03, 0.9563949E-03, 0.2263770E-02, &
         0.5366995E-02, 0.1274226E-01, 0.3023528E-01, 0.7115984E-01, &
         0.1623103E+00, 0.3391695E+00, 0.1382336E-03, 0.5461978E-11, &
         0.1045634E-10, 0.2480865E-10, 0.5877250E-10, 0.1386328E-09, &
         0.3260322E-09, 0.7689554E-09, 0.1809900E-08, 0.4257048E-08, &
         0.1004131E-07, 0.2366112E-07, 0.5576113E-07, 0.1314203E-06, &
         0.3096752E-06, 0.7299032E-06, 0.1720274E-05, 0.4054981E-05, &
         0.9559084E-05, 0.2253919E-04, 0.5315715E-04, 0.1254218E-03, &
         0.2961037E-03, 0.6996637E-03, 0.1655233E-02, 0.3921884E-02, &
         0.9307324E-02, 0.2209932E-01, 0.5224660E-01, 0.1210630E+00, &
         0.2629204E-03, 0.8523221E-12, 0.1655710E-11, 0.3974937E-11, &
         0.9475196E-11, 0.2160493E-10, 0.5173064E-10, 0.1217817E-09, &
         0.2857532E-09, 0.6785524E-09, 0.1597406E-08, 0.3767277E-08, &
         0.8871726E-08, 0.2090630E-07, 0.4927249E-07, 0.1161264E-06, &
         0.2736318E-06, 0.6449878E-06, 0.1520047E-05, 0.3583024E-05, &
         0.8446612E-05, 0.1991419E-04, 0.4696518E-04, 0.1108043E-03, &
         0.2615690E-03, 0.6179706E-03, 0.1461678E-02, 0.3462458E-02, &
         0.8215324E-02, 0.1950821E-01, 0.5000746E-03, 0.2305856E-02, &
         0.5466960E-02, 0.1297981E-01, 0.3079700E-01, 0.7245710E-01, &
         0.1650817E+00, 0.3440087E+00, 0.6102081E+00, 0.9822114E+00, &
         0.1416896E+01, 0.1905245E+01, 0.2396360E+01, 0.2865409E+01, &
         0.3078600E+01, 0.2821023E+01, 0.2079308E+01, 0.1451826E+01, &
         0.2001947E+01, 0.2059541E+01, 0.1667592E+01, 0.1893693E+01, &
         0.1813272E+01, 0.1784861E+01, 0.1766699E+01, 0.1802252E+01, &
         0.1777354E+01, 0.1773361E+01, 0.1769136E+01, 0.1735761E+01, &
         0.8336438E-03, 0.1607949E-02, 0.3809640E-02, 0.9040533E-02, &
         0.2146645E-01, 0.5076601E-01, 0.1177685E+00, 0.2571919E+00, &
         0.4902766E+00, 0.8061812E+00, 0.1250220E+01, 0.1660402E+01, &
         0.2188359E+01, 0.2646668E+01, 0.2970078E+01, 0.3017053E+01, &
         0.2416580E+01, 0.1578704E+01, 0.1703342E+01, 0.2210760E+01, &
         0.1668151E+01, 0.1938131E+01, 0.1788607E+01, 0.1854164E+01, &
         0.1828766E+01, 0.1805977E+01, 0.1766851E+01, 0.1782201E+01, &
         0.1754907E+01, 0.1744322E+01, 0.6669299E-03, 0.1286015E-02, &
         0.3045615E-02, 0.7224638E-02, 0.1715594E-01, 0.4064883E-01, &
         0.9501187E-01, 0.2120690E+00, 0.4218818E+00, 0.7130325E+00, &
         0.1132582E+01, 0.1539426E+01, 0.2076132E+01, 0.2551370E+01, &
         0.2925028E+01, 0.3023993E+01, 0.2621071E+01, 0.1807222E+01, &
         0.1550507E+01, 0.2171909E+01, 0.1841270E+01, 0.1927457E+01, &
         0.1817936E+01, 0.1826131E+01, 0.1748519E+01, 0.1781892E+01, &
         0.1787518E+01, 0.1738967E+01, 0.1760810E+01, 0.1754529E+01/),&
         (/30,1,6/))

    qscat(1:30,1:1,13:18)=reshape( (/&
         0.4939059E-03, 0.9520435E-03, 0.2253458E-02, 0.5342496E-02, &
         0.1268401E-01, 0.3009751E-01, 0.7084160E-01, 0.1616292E+00, &
         0.3379757E+00, 0.6021811E+00, 0.9701180E+00, 0.1407054E+01, &
         0.1889471E+01, 0.2381627E+01, 0.2852786E+01, 0.3083985E+01, &
         0.2840933E+01, 0.2098689E+01, 0.1445787E+01, 0.2000638E+01, &
         0.2100897E+01, 0.1701578E+01, 0.1950433E+01, 0.1846692E+01, &
         0.1828247E+01, 0.1768301E+01, 0.1749795E+01, 0.1776199E+01, &
         0.1766105E+01, 0.1769362E+01, 0.4058536E-03, 0.7821483E-03, &
         0.1850704E-02, 0.4385977E-02, 0.1041047E-01, 0.2471486E-01, &
         0.5835165E-01, 0.1345477E+00, 0.2890299E+00, 0.5356299E+00, &
         0.8710567E+00, 0.1318777E+01, 0.1752060E+01, 0.2261234E+01, &
         0.2721528E+01, 0.3050559E+01, 0.2949300E+01, 0.2250151E+01, &
         0.1530140E+01, 0.1829173E+01, 0.2215410E+01, 0.1710420E+01, &
         0.1978602E+01, 0.1811628E+01, 0.1824998E+01, 0.1844743E+01, &
         0.1773051E+01, 0.1760529E+01, 0.1747295E+01, 0.1756304E+01, &
         0.2877247E-03, 0.5543064E-03, 0.1310881E-02, 0.3104613E-02, &
         0.7364837E-02, 0.1748890E-01, 0.4143228E-01, 0.9678915E-01, &
         0.2156751E+00, 0.4275572E+00, 0.7205915E+00, 0.1142951E+01, &
         0.1548665E+01, 0.2086365E+01, 0.2559992E+01, 0.2925613E+01, &
         0.3027108E+01, 0.2608634E+01, 0.1765504E+01, 0.1553489E+01, &
         0.2209895E+01, 0.1835341E+01, 0.1919490E+01, 0.1788675E+01, &
         0.1842055E+01, 0.1772597E+01, 0.1804116E+01, 0.1782525E+01, &
         0.1742764E+01, 0.1744599E+01, 0.1639157E-03, 0.3156405E-03, &
         0.7458854E-03, 0.1764762E-02, 0.4181915E-02, 0.9925411E-02, &
         0.2356508E-01, 0.5567071E-01, 0.1286458E+00, 0.2779670E+00, &
         0.5200978E+00, 0.8485439E+00, 0.1296149E+01, 0.1719943E+01, &
         0.2235500E+01, 0.2693254E+01, 0.3024626E+01, 0.2977588E+01, &
         0.2311464E+01, 0.1667730E+01, 0.1821186E+01, 0.2224497E+01, &
         0.1704216E+01, 0.1982319E+01, 0.1825384E+01, 0.1832962E+01, &
         0.1851820E+01, 0.1797334E+01, 0.1770782E+01, 0.1758385E+01, &
         0.1360369E-03, 0.2619239E-03, 0.6188128E-03, 0.1463670E-02, &
         0.3467186E-02, 0.8226552E-02, 0.1953486E-01, 0.4623934E-01, &
         0.1076391E+00, 0.2373873E+00, 0.4609284E+00, 0.7655962E+00, &
         0.1201772E+01, 0.1605729E+01, 0.2142393E+01, 0.2606626E+01, &
         0.2936520E+01, 0.3031624E+01, 0.2534039E+01, 0.1769609E+01, &
         0.1651198E+01, 0.2273975E+01, 0.1833097E+01, 0.1972575E+01, &
         0.1789315E+01, 0.1821071E+01, 0.1794554E+01, 0.1826451E+01, &
         0.1784653E+01, 0.1767444E+01, 0.1044451E-03, 0.2010619E-03, &
         0.4748902E-03, 0.1122807E-02, 0.2658433E-02, 0.6304602E-02, &
         0.1497037E-01, 0.3549854E-01, 0.8326766E-01, 0.1878861E+00, &
         0.3827521E+00, 0.6613459E+00, 0.1058644E+01, 0.1477897E+01, &
         0.1997881E+01, 0.2482627E+01, 0.2914294E+01, 0.3026468E+01, &
         0.2727181E+01, 0.1884788E+01, 0.1458210E+01, 0.2130828E+01, &
         0.1948996E+01, 0.1807852E+01, 0.1839294E+01, 0.1855008E+01, &
         0.1810482E+01, 0.1750667E+01, 0.1797086E+01, 0.1775596E+01/),&
         (/30,1,6/))

    qscat(1:30,1:1,19:24)=reshape( (/&
         0.9441071E-04, 0.1817370E-03, 0.4291996E-03, 0.1014639E-02, &
         0.2401892E-02, 0.5695102E-02, 0.1352202E-01, 0.3207864E-01, &
         0.7541227E-01, 0.1713672E+00, 0.3548766E+00, 0.6246089E+00, &
         0.1003889E+01, 0.1434349E+01, 0.1932783E+01, 0.2422200E+01, &
         0.2884667E+01, 0.3065633E+01, 0.2791466E+01, 0.2049569E+01, &
         0.1575594E+01, 0.2073899E+01, 0.2041062E+01, 0.1722158E+01, &
         0.1899726E+01, 0.1846775E+01, 0.1793869E+01, 0.1778048E+01, &
         0.1770816E+01, 0.1789440E+01, 0.8151573E-04, 0.1569007E-03, &
         0.3704934E-03, 0.8756850E-03, 0.2072421E-02, 0.4912498E-02, &
         0.1166200E-01, 0.2767925E-01, 0.6524293E-01, 0.1495746E+00, &
         0.3165503E+00, 0.5734129E+00, 0.9269114E+00, 0.1370421E+01, &
         0.1830986E+01, 0.2328271E+01, 0.2798606E+01, 0.3087561E+01, &
         0.2882272E+01, 0.2176012E+01, 0.1618124E+01, 0.1939663E+01, &
         0.2155108E+01, 0.1668577E+01, 0.1943499E+01, 0.1805721E+01, &
         0.1826187E+01, 0.1817263E+01, 0.1757304E+01, 0.1767201E+01, &
         0.5565022E-04, 0.1070932E-03, 0.2528019E-03, 0.5972348E-03, &
         0.1412555E-02, 0.3345877E-02, 0.7938247E-02, 0.1885045E-01, &
         0.4463263E-01, 0.1040232E+00, 0.2302090E+00, 0.4500423E+00, &
         0.7507968E+00, 0.1183003E+01, 0.1586547E+01, 0.2124766E+01, &
         0.2591938E+01, 0.2930691E+01, 0.3031670E+01, 0.2565965E+01, &
         0.1732758E+01, 0.1644230E+01, 0.2178035E+01, 0.1792361E+01, &
         0.1880221E+01, 0.1821903E+01, 0.1829716E+01, 0.1787137E+01, &
         0.1815011E+01, 0.1770356E+01, 0.3980661E-04, 0.7659413E-04, &
         0.1807602E-03, 0.4268896E-03, 0.1009173E-02, 0.2388930E-02, &
         0.5664303E-02, 0.1344884E-01, 0.3190569E-01, 0.7501386E-01, &
         0.1705220E+00, 0.3534238E+00, 0.6226878E+00, 0.1001000E+01, &
         0.1432039E+01, 0.1929170E+01, 0.2418807E+01, 0.2882358E+01, &
         0.3067569E+01, 0.2794610E+01, 0.2053816E+01, 0.1557362E+01, &
         0.2052078E+01, 0.2018643E+01, 0.1696542E+01, 0.1886784E+01, &
         0.1840827E+01, 0.1781331E+01, 0.1751741E+01, 0.1756093E+01, &
         0.2681827E-04, 0.5159516E-04, 0.1217334E-03, 0.2873892E-03, &
         0.6790494E-03, 0.1606384E-02, 0.3805922E-02, 0.9031690E-02, &
         0.2144548E-01, 0.5071694E-01, 0.1176591E+00, 0.2569804E+00, &
         0.4899683E+00, 0.8057491E+00, 0.1249726E+01, 0.1659806E+01, &
         0.2187880E+01, 0.2646229E+01, 0.2969582E+01, 0.3017386E+01, &
         0.2417457E+01, 0.1579919E+01, 0.1704027E+01, 0.2208308E+01, &
         0.1668819E+01, 0.1938163E+01, 0.1787129E+01, 0.1846662E+01, &
         0.1827163E+01, 0.1803410E+01, 0.2525160E-04, 0.4857985E-04, &
         0.1146153E-03, 0.2705707E-03, 0.6392652E-03, 0.1512128E-02, &
         0.3582197E-02, 0.8499914E-02, 0.2018366E-01, 0.4776109E-01, &
         0.1110539E+00, 0.2441137E+00, 0.4710058E+00, 0.7794141E+00, &
         0.1218766E+01, 0.1624015E+01, 0.2158341E+01, 0.2620128E+01, &
         0.2944883E+01, 0.3030327E+01, 0.2490729E+01, 0.1736488E+01, &
         0.1670145E+01, 0.2207124E+01, 0.1796909E+01, 0.2016083E+01, &
         0.1831616E+01, 0.1820178E+01, 0.1820695E+01, 0.1821183E+01/),&
         (/30,1,6/))

    qscat(1:30,1:1,25:30)=reshape( (/&
         0.2244690E-04, 0.4318224E-04, 0.1018737E-03, 0.2404698E-03, &
         0.5680734E-03, 0.1343483E-02, 0.3181970E-02, 0.7548691E-02, &
         0.1792551E-01, 0.4245912E-01, 0.9911473E-01, 0.2203728E+00, &
         0.4348921E+00, 0.7303964E+00, 0.1156200E+01, 0.1560783E+01, &
         0.2099243E+01, 0.2570755E+01, 0.2926596E+01, 0.3029954E+01, &
         0.2592838E+01, 0.1722351E+01, 0.1611468E+01, 0.2245607E+01, &
         0.1816836E+01, 0.1897303E+01, 0.1807774E+01, 0.1827002E+01, &
         0.1780776E+01, 0.1792670E+01, 0.1508714E-04, 0.2902153E-04, &
         0.6845398E-04, 0.1615359E-03, 0.3814493E-03, 0.9016190E-03, &
         0.2133902E-02, 0.5058516E-02, 0.1200905E-01, 0.2850071E-01, &
         0.6714708E-01, 0.1536901E+00, 0.3239290E+00, 0.5833729E+00, &
         0.9418291E+00, 0.1383321E+01, 0.1851519E+01, 0.2346708E+01, &
         0.2818617E+01, 0.3089535E+01, 0.2873804E+01, 0.2160190E+01, &
         0.1588495E+01, 0.1916183E+01, 0.2115090E+01, 0.1767522E+01, &
         0.1952622E+01, 0.1827966E+01, 0.1801739E+01, 0.1806674E+01, &
         0.1903640E-04, 0.3661983E-04, 0.8638486E-04, 0.2038823E-03, &
         0.4815583E-03, 0.1138599E-02, 0.2695894E-02, 0.6393607E-02, &
         0.1518182E-01, 0.3599739E-01, 0.8440971E-01, 0.1902647E+00, &
         0.3866871E+00, 0.6665245E+00, 0.1066252E+01, 0.1483992E+01, &
         0.2006419E+01, 0.2490376E+01, 0.2916542E+01, 0.3022475E+01, &
         0.2712411E+01, 0.1870931E+01, 0.1473295E+01, 0.2130482E+01, &
         0.1957074E+01, 0.1783501E+01, 0.1843272E+01, 0.1833044E+01, &
         0.1786526E+01, 0.1781775E+01, 0.1164130E-04, 0.2239145E-04, &
         0.5281003E-04, 0.1246017E-03, 0.2941640E-03, 0.6950784E-03, &
         0.1644360E-02, 0.3896082E-02, 0.9245978E-02, 0.2195383E-01, &
         0.5190634E-01, 0.1203065E+00, 0.2620847E+00, 0.4973831E+00, &
         0.8161713E+00, 0.1261464E+01, 0.1674268E+01, 0.2199425E+01, &
         0.2656993E+01, 0.2982174E+01, 0.3008611E+01, 0.2397170E+01, &
         0.1566271E+01, 0.1716662E+01, 0.2299051E+01, 0.1702490E+01, &
         0.1945370E+01, 0.1792890E+01, 0.1841396E+01, 0.1849848E+01, &
         0.9210702E-05, 0.1771512E-04, 0.4177722E-04, 0.9855754E-04, &
         0.2326356E-03, 0.5495443E-03, 0.1299600E-02, 0.3077847E-02, &
         0.7301237E-02, 0.1733785E-01, 0.4107690E-01, 0.9598324E-01, &
         0.2140418E+00, 0.4249914E+00, 0.7171717E+00, 0.1138275E+01, &
         0.1544475E+01, 0.2081769E+01, 0.2556129E+01, 0.2925340E+01, &
         0.3025771E+01, 0.2614344E+01, 0.1784989E+01, 0.1548554E+01, &
         0.2207918E+01, 0.1856935E+01, 0.1970010E+01, 0.1816692E+01, &
         0.1819728E+01, 0.1765595E+01, 0.6915328E-05, 0.1330029E-04, &
         0.3136207E-04, 0.7397675E-04, 0.1745795E-03, 0.4122800E-03, &
         0.9745910E-03, 0.2306919E-02, 0.5469492E-02, 0.1298585E-01, &
         0.3081126E-01, 0.7249000E-01, 0.1651520E+00, 0.3441308E+00, &
         0.6103705E+00, 0.9824560E+00, 0.1417094E+01, 0.1905561E+01, &
         0.2396655E+01, 0.2865654E+01, 0.3078475E+01, 0.2820626E+01, &
         0.2078964E+01, 0.1452464E+01, 0.2002201E+01, 0.2060558E+01, &
         0.1667076E+01, 0.1909250E+01, 0.1829093E+01, 0.1784407E+01/),&
         (/30,1,6/))

    qscat(1:30,1:1,31:36)=reshape( (/&
         0.4983547E-05, 0.9583965E-05, 0.2259792E-04, 0.5329637E-04, &
         0.1257496E-03, 0.2968784E-03, 0.7014993E-03, 0.1659575E-02, &
         0.3932194E-02, 0.9331833E-02, 0.2215747E-01, 0.5238259E-01, &
         0.1213650E+00, 0.2641168E+00, 0.5003189E+00, 0.8203176E+00, &
         0.1266057E+01, 0.1680060E+01, 0.2204016E+01, 0.2661379E+01, &
         0.2987473E+01, 0.3004859E+01, 0.2388839E+01, 0.1569596E+01, &
         0.1731562E+01, 0.2236392E+01, 0.1710583E+01, 0.1946052E+01, &
         0.1803016E+01, 0.1828850E+01, 0.3336525E-05, 0.6416181E-05, &
         0.1512720E-04, 0.3567183E-04, 0.8414762E-04, 0.1985995E-03, &
         0.4690658E-03, 0.1109019E-02, 0.2625725E-02, 0.6226893E-02, &
         0.1478571E-01, 0.3506286E-01, 0.8226940E-01, 0.1858021E+00, &
         0.3792889E+00, 0.6567879E+00, 0.1051920E+01, 0.1472531E+01, &
         0.1990234E+01, 0.2475642E+01, 0.2911975E+01, 0.3030586E+01, &
         0.2739404E+01, 0.1902393E+01, 0.1459937E+01, 0.2141851E+01, &
         0.1998446E+01, 0.1849737E+01, 0.1843773E+01, 0.1828411E+01, &
         0.1609733E-05, 0.3095464E-05, 0.7296939E-05, 0.1720301E-04, &
         0.4056917E-04, 0.9570643E-04, 0.2259008E-03, 0.5336171E-03, &
         0.1261881E-02, 0.2988344E-02, 0.7088546E-02, 0.1683272E-01, &
         0.3988795E-01, 0.9328347E-01, 0.2085486E+00, 0.4163024E+00, &
         0.7056225E+00, 0.1122294E+01, 0.1530449E+01, 0.2065823E+01, &
         0.2542608E+01, 0.2924463E+01, 0.3020677E+01, 0.2632277E+01, &
         0.1835126E+01, 0.1563670E+01, 0.2162810E+01, 0.1854607E+01, &
         0.1907334E+01, 0.1804355E+01, 0.8536749E-06, 0.1641531E-05, &
         0.3869212E-05, 0.9121159E-05, 0.2150566E-04, 0.5071913E-04, &
         0.1196660E-03, 0.2825053E-03, 0.6674958E-03, 0.1579006E-02, &
         0.3740940E-02, 0.8877230E-02, 0.2107904E-01, 0.4985903E-01, &
         0.1157458E+00, 0.2532728E+00, 0.4845448E+00, 0.7981708E+00, &
         0.1241009E+01, 0.1649390E+01, 0.2179454E+01, 0.2638595E+01, &
         0.2961328E+01, 0.3022736E+01, 0.2434684E+01, 0.1610582E+01, &
         0.1719000E+01, 0.2209665E+01, 0.1695634E+01, 0.1954261E+01, &
         0.5176838E-06, 0.9954088E-06, 0.2346362E-05, 0.5530849E-05, &
         0.1303859E-04, 0.3074600E-04, 0.7252304E-04, 0.1711466E-03, &
         0.4041660E-03, 0.9553821E-03, 0.2261380E-02, 0.5361313E-02, &
         0.1272874E-01, 0.3020329E-01, 0.7108603E-01, 0.1621523E+00, &
         0.3388926E+00, 0.6034029E+00, 0.9719585E+00, 0.1408524E+01, &
         0.1891890E+01, 0.2383881E+01, 0.2854789E+01, 0.3083266E+01, &
         0.2837948E+01, 0.2095346E+01, 0.1443578E+01, 0.1994793E+01, &
         0.2100053E+01, 0.1695583E+01, 0.2871854E-06, 0.5520937E-06, &
         0.1301202E-05, 0.3067068E-05, 0.7230215E-05, 0.1704622E-04, &
         0.4019847E-04, 0.9483170E-04, 0.2238343E-03, 0.5287331E-03, &
         0.1250309E-02, 0.2960899E-02, 0.7023325E-02, 0.1667778E-01, &
         0.3952316E-01, 0.9245398E-01, 0.2068544E+00, 0.4136037E+00, &
         0.7020445E+00, 0.1117285E+01, 0.1526141E+01, 0.2060742E+01, &
         0.2538260E+01, 0.2924151E+01, 0.3019181E+01, 0.2637346E+01, &
         0.1843114E+01, 0.1572628E+01, 0.2170437E+01, 0.1903235E+01/),&
         (/30,1,6/))

    qscat(1:30,1:1,37:42)=reshape( (/&
         0.2130258E-06, 0.4095630E-06, 0.9652463E-06, 0.2275018E-05, &
         0.5362792E-05, 0.1264280E-04, 0.2981126E-04, 0.7031697E-04, &
         0.1659372E-03, 0.3918542E-03, 0.9262424E-03, 0.2192283E-02, &
         0.5197192E-02, 0.1233868E-01, 0.2928065E-01, 0.6895282E-01, &
         0.1575780E+00, 0.3308383E+00, 0.5926470E+00, 0.9557645E+00, &
         0.1395105E+01, 0.1870389E+01, 0.2363955E+01, 0.2836186E+01, &
         0.3088190E+01, 0.2861268E+01, 0.2129971E+01, 0.1498628E+01, &
         0.1965536E+01, 0.2156023E+01, 0.1139642E-06, 0.2191180E-06, &
         0.5163641E-06, 0.1217040E-05, 0.2868713E-05, 0.6762461E-05, &
         0.1594324E-04, 0.3759716E-04, 0.8869241E-04, 0.2093319E-03, &
         0.4944418E-03, 0.1169102E-02, 0.2768247E-02, 0.6565538E-02, &
         0.1559031E-01, 0.3696075E-01, 0.8661259E-01, 0.1948363E+00, &
         0.3941952E+00, 0.6764117E+00, 0.1080671E+01, 0.1495647E+01, &
         0.2022275E+01, 0.2504614E+01, 0.2919853E+01, 0.3017293E+01, &
         0.2684913E+01, 0.1859088E+01, 0.1536330E+01, 0.2138908E+01, &
         0.4368971E-07, 0.8401571E-07, 0.1979443E-06, 0.4665041E-06, &
         0.1099423E-05, 0.2591267E-05, 0.6108387E-05, 0.1440077E-04, &
         0.3395877E-04, 0.8010484E-04, 0.1890504E-03, 0.4464910E-03, &
         0.1055577E-02, 0.2498972E-02, 0.5925733E-02, 0.1407011E-01, &
         0.3337342E-01, 0.7839159E-01, 0.1776647E+00, 0.3656191E+00, &
         0.6387883E+00, 0.1025151E+01, 0.1451272E+01, 0.1958829E+01, &
         0.2446588E+01, 0.2899261E+01, 0.3050108E+01, 0.2773031E+01, &
         0.2001010E+01, 0.1567242E+01, 0.2413002E-07, 0.4641096E-07, &
         0.1093743E-06, 0.2577491E-06, 0.6074086E-06, 0.1431569E-05, &
         0.3374563E-05, 0.7954735E-05, 0.1875470E-04, 0.4422979E-04, &
         0.1043463E-03, 0.2463109E-03, 0.5818865E-03, 0.1376202E-02, &
         0.3259612E-02, 0.7733199E-02, 0.1836365E-01, 0.4348901E-01, &
         0.1014429E+00, 0.2250517E+00, 0.4421329E+00, 0.7401199E+00, &
         0.1169107E+01, 0.1572968E+01, 0.2111608E+01, 0.2581027E+01, &
         0.2928094E+01, 0.3031281E+01, 0.2579799E+01, 0.1712282E+01, &
         0.1497814E-07, 0.2878358E-07, 0.6784342E-07, 0.1598607E-06, &
         0.3767407E-06, 0.8879427E-06, 0.2092868E-05, 0.4933270E-05, &
         0.1163001E-04, 0.2742337E-04, 0.6468265E-04, 0.1526333E-03, &
         0.3604091E-03, 0.8518165E-03, 0.2015840E-02, 0.4778121E-02, &
         0.1134260E-01, 0.2692300E-01, 0.6348782E-01, 0.1457673E+00, &
         0.3096642E+00, 0.5640594E+00, 0.9129618E+00, 0.1358052E+01, &
         0.1811523E+01, 0.2311157E+01, 0.2779208E+01, 0.3082438E+01, &
         0.2892156E+01, 0.2187006E+01, 0.1030113E-07, 0.1980520E-07, &
         0.4668839E-07, 0.1100485E-06, 0.2593487E-06, 0.6112119E-06, &
         0.1440549E-05, 0.3395485E-05, 0.8004393E-05, 0.1887181E-04, &
         0.4450613E-04, 0.1049992E-03, 0.2478533E-03, 0.5855319E-03, &
         0.1384842E-02, 0.3280106E-02, 0.7781909E-02, 0.1847929E-01, &
         0.4376074E-01, 0.1020566E+00, 0.2262809E+00, 0.4440247E+00, &
         0.7426682E+00, 0.1172451E+01, 0.1576190E+01, 0.2114785E+01, &
         0.2583662E+01, 0.2928657E+01, 0.3031442E+01, 0.2576708E+01/),&
         (/30,1,6/))

    qscat(1:30,1:1,43:48)=reshape( (/&
         0.8111209E-08, 0.1559276E-07, 0.3674471E-07, 0.8658950E-07, &
         0.2040660E-06, 0.4809785E-06, 0.1133488E-05, 0.2671544E-05, &
         0.6297743E-05, 0.1484709E-04, 0.3501190E-04, 0.8259039E-04, &
         0.1949215E-03, 0.4603695E-03, 0.1088433E-02, 0.2576899E-02, &
         0.6110877E-02, 0.1451004E-01, 0.3441217E-01, 0.8077721E-01, &
         0.1826787E+00, 0.3740697E+00, 0.6499188E+00, 0.1041741E+01, &
         0.1464437E+01, 0.1978474E+01, 0.2464825E+01, 0.2907833E+01, &
         0.3037615E+01, 0.2755203E+01, 0.6604543E-08, 0.1270030E-07, &
         0.2991808E-07, 0.7050917E-07, 0.1661954E-06, 0.3916392E-06, &
         0.9230437E-06, 0.2175630E-05, 0.5128434E-05, 0.1209001E-04, &
         0.2850855E-04, 0.6724245E-04, 0.1586770E-03, 0.3746897E-03, &
         0.8856207E-03, 0.2095975E-02, 0.4968437E-02, 0.1179495E-01, &
         0.2799396E-01, 0.6597276E-01, 0.1511538E+00, 0.3193898E+00, &
         0.5772527E+00, 0.9326557E+00, 0.1375426E+01, 0.1838931E+01, &
         0.2335359E+01, 0.2806430E+01, 0.3088753E+01, 0.2879159E+01, &
         0.5314572E-08, 0.1022444E-07, 0.2410669E-07, 0.5677523E-07, &
         0.1338119E-06, 0.3153663E-06, 0.7432756E-06, 0.1751925E-05, &
         0.4129471E-05, 0.9734469E-05, 0.2295285E-04, 0.5413433E-04, &
         0.1277282E-03, 0.3015520E-03, 0.7125568E-03, 0.1685775E-02, &
         0.3994397E-02, 0.9479688E-02, 0.2250814E-01, 0.5320236E-01, &
         0.1231845E+00, 0.2675992E+00, 0.5053291E+00, 0.8274198E+00, &
         0.1273818E+01, 0.1690028E+01, 0.2211889E+01, 0.2669053E+01, &
         0.2996819E+01, 0.2998283E+01, 0.4225172E-08, 0.8133153E-08, &
         0.1916917E-07, 0.4515045E-07, 0.1064213E-06, 0.2507710E-06, &
         0.5910236E-06, 0.1392893E-05, 0.3283100E-05, 0.7739130E-05, &
         0.1824681E-04, 0.4303140E-04, 0.1015181E-03, 0.2396298E-03, &
         0.5660809E-03, 0.1338762E-02, 0.3170775E-02, 0.7522085E-02, &
         0.1786233E-01, 0.4231057E-01, 0.9877851E-01, 0.2196952E+00, &
         0.4338378E+00, 0.7289849E+00, 0.1154307E+01, 0.1559031E+01, &
         0.2097416E+01, 0.2569233E+01, 0.2926430E+01, 0.3029643E+01, &
         0.3311696E-08, 0.6376909E-08, 0.1502053E-07, 0.3541976E-07, &
         0.8341913E-07, 0.1966386E-06, 0.4633977E-06, 0.1092155E-05, &
         0.2574293E-05, 0.6068137E-05, 0.1430610E-04, 0.3373469E-04, &
         0.7957646E-04, 0.1878019E-03, 0.4435409E-03, 0.1048594E-02, &
         0.2482403E-02, 0.5886369E-02, 0.1397656E-01, 0.3315251E-01, &
         0.7788367E-01, 0.1765939E+00, 0.3638026E+00, 0.6363935E+00, &
         0.1021569E+01, 0.1448428E+01, 0.1954509E+01, 0.2442554E+01, &
         0.2897095E+01, 0.3052830E+01, 0.2558936E-08, 0.4921040E-08, &
         0.1160372E-07, 0.2734706E-07, 0.6442724E-07, 0.1518185E-06, &
         0.3578024E-06, 0.8431982E-06, 0.1987328E-05, 0.4684338E-05, &
         0.1104337E-04, 0.2603939E-04, 0.6141687E-04, 0.1449215E-03, &
         0.3421815E-03, 0.8086877E-03, 0.1913604E-02, 0.4535340E-02, &
         0.1076550E-01, 0.2555612E-01, 0.6031040E-01, 0.1388400E+00, &
         0.2969850E+00, 0.5466675E+00, 0.8872278E+00, 0.1334345E+01, &
         0.1775121E+01, 0.2280190E+01, 0.2743268E+01, 0.3065531E+01/),&
         (/30,1,6/))

    qscat(1:30,1:1,49:50)=reshape( (/&
         0.1941416E-08, 0.3729086E-08, 0.8796222E-08, 0.2073586E-07, &
         0.4884795E-07, 0.1151409E-06, 0.2713934E-06, 0.6395728E-06, &
         0.1507417E-05, 0.3552973E-05, 0.8375334E-05, 0.1974695E-04, &
         0.4657068E-04, 0.1098731E-03, 0.2593660E-03, 0.6127629E-03, &
         0.1449341E-02, 0.3433179E-02, 0.8145735E-02, 0.1934303E-01, &
         0.4578911E-01, 0.1066269E+00, 0.2353837E+00, 0.4579038E+00, &
         0.7614718E+00, 0.1196600E+01, 0.1600340E+01, 0.2137544E+01, &
         0.2602569E+01, 0.2934595E+01, 0.1327624E-08, 0.2548038E-08, &
         0.5997629E-08, 0.1413148E-07, 0.3332353E-07, 0.7850970E-07, &
         0.1850432E-06, 0.4361592E-06, 0.1027941E-05, 0.2423005E-05, &
         0.5711535E-05, 0.1346519E-04, 0.3175162E-04, 0.7489492E-04, &
         0.1767484E-03, 0.4174066E-03, 0.9867256E-03, 0.2335701E-02, &
         0.5537855E-02, 0.1314833E-01, 0.3119540E-01, 0.7337637E-01, &
         0.1670410E+00, 0.3474123E+00, 0.6147258E+00, 0.9890165E+00, &
         0.1422405E+01, 0.1913994E+01, 0.2404561E+01, 0.2871936E+01/),&
         (/30,1,2/))

    qbrqs(1:30,1:1,1:6)=reshape( (/&
         0.1901999E-05, 0.4716377E-12, 0.7370087E-12, 0.3116265E-12, &
         -0.2044577E-11, 0.1364757E-11, 0.2750097E-11, 0.1206867E-10, &
         0.7804748E-10, 0.2244802E-09, 0.7640885E-09, 0.2786042E-08, &
         0.1003014E-07, 0.3626222E-07, 0.1312045E-06, 0.4747086E-06, &
         0.1717792E-05, 0.6216306E-05, 0.2250494E-04, 0.8146375E-04, &
         0.2948020E-03, 0.1065178E-02, 0.3828139E-02, 0.1352974E-01, &
         0.4564070E-01, 0.1391708E+00, 0.3677734E+00, 0.7288581E+00, &
         0.1024653E+01, 0.1386942E+01, 0.3617603E-05, 0.2991549E-12, &
         -0.5700891E-12, 0.1586302E-11, 0.9922366E-12, 0.1795592E-11, &
         -0.3836203E-11, 0.1123241E-10, 0.2974527E-10, 0.8357062E-10, &
         0.3698330E-09, 0.1441129E-08, 0.5219182E-08, 0.1896117E-07, &
         0.6871815E-07, 0.2483693E-06, 0.8989706E-06, 0.3253124E-05, &
         0.1177489E-04, 0.4262225E-04, 0.1542823E-03, 0.5580288E-03, &
         0.2012467E-02, 0.7188866E-02, 0.2496183E-01, 0.8068376E-01, &
         0.2299847E+00, 0.5466217E+00, 0.8696541E+00, 0.1232296E+01, &
         0.6880682E-05, 0.1916953E-12, 0.4707819E-12, 0.3414841E-12, &
         0.4034255E-12, 0.1418572E-12, 0.3161054E-12, 0.1475615E-11, &
         0.1669863E-10, 0.3005890E-10, 0.1336998E-09, 0.5276825E-09, &
         0.1935742E-08, 0.6909711E-08, 0.2535454E-07, 0.9163686E-07, &
         0.3316923E-06, 0.1200273E-05, 0.4344654E-05, 0.1572453E-04, &
         0.5692060E-04, 0.2060147E-03, 0.7448617E-03, 0.2682767E-02, &
         0.9544524E-02, 0.3276870E-01, 0.1033715E+00, 0.2852662E+00, &
         0.6315490E+00, 0.9307482E+00, 0.1308705E-04, 0.3488645E-12, &
         -0.2643488E-12, 0.3843259E-12, 0.5917927E-12, 0.1269802E-11, &
         -0.1910060E-11, 0.4475100E-11, 0.2189990E-11, 0.1775183E-10, &
         0.7611258E-10, 0.2165611E-09, 0.8112211E-09, 0.2849104E-08, &
         0.1058281E-07, 0.3819969E-07, 0.1382846E-06, 0.5006514E-06, &
         0.1810602E-05, 0.6554102E-05, 0.2372228E-04, 0.8587552E-04, &
         0.3107595E-03, 0.1122709E-02, 0.4033369E-02, 0.1423876E-01, &
         0.4789057E-01, 0.1452604E+00, 0.3812047E+00, 0.7420156E+00, &
         0.2489158E-04, 0.4701444E-13, 0.1535593E-12, 0.3381911E-12, &
         0.5717167E-12, 0.4095736E-12, 0.2617854E-12, 0.3672750E-11, &
         0.3647954E-11, 0.8980333E-11, 0.2740179E-10, 0.1001760E-09, &
         0.3395771E-09, 0.1195295E-08, 0.4350772E-08, 0.1593059E-07, &
         0.5727330E-07, 0.2076147E-06, 0.7506055E-06, 0.2718510E-05, &
         0.9838780E-05, 0.3561608E-04, 0.1289225E-03, 0.4663866E-03, &
         0.1683064E-02, 0.6024596E-02, 0.2104279E-01, 0.6894487E-01, &
         0.2004059E+00, 0.4936804E+00, 0.4734378E-04, 0.1351300E-13, &
         -0.1268058E-12, 0.2985778E-13, 0.3237699E-12, 0.6953172E-13, &
         -0.5454742E-12, 0.2173147E-11, 0.1272297E-11, 0.2537878E-11, &
         0.6615754E-11, 0.2416349E-10, 0.1269510E-09, 0.4061239E-09, &
         0.1475770E-08, 0.5395762E-08, 0.1915448E-07, 0.6942866E-07, &
         0.2513579E-06, 0.9103257E-06, 0.3293702E-05, 0.1192200E-04, &
         0.4315420E-04, 0.1562069E-03, 0.5649814E-03, 0.2037445E-02, &
         0.7276999E-02, 0.2525676E-01, 0.8155691E-01, 0.2321569E+00/),&
         (/30,1,6/))

    qbrqs(1:30,1:1,7:12)=reshape( (/&
         0.9004787E-04, 0.2465319E-14, 0.6956497E-13, 0.7135020E-13, &
         -0.7006778E-13, 0.4401354E-12, 0.4894050E-12, 0.6742820E-12, &
         -0.1286480E-11, 0.2726573E-12, 0.1833192E-11, 0.4639781E-11, &
         0.2582518E-10, 0.1105154E-09, 0.3573153E-09, 0.1309629E-08, &
         0.4652057E-08, 0.1709559E-07, 0.6191839E-07, 0.2242348E-06, &
         0.8105023E-06, 0.2933317E-05, 0.1061740E-04, 0.3843134E-04, &
         0.1391129E-03, 0.5032207E-03, 0.1815523E-02, 0.6493290E-02, &
         0.2262557E-01, 0.7371690E-01, 0.1712710E-03, 0.1129701E-13, &
         -0.3926947E-13, 0.2638894E-13, 0.5981981E-13, 0.1928346E-12, &
         0.1923591E-12, 0.7062969E-12, 0.8035807E-12, 0.1282500E-11, &
         0.7179135E-12, 0.1216088E-11, 0.4708282E-11, 0.2094201E-10, &
         0.5307086E-10, 0.2274950E-09, 0.8008884E-09, 0.2959651E-08, &
         0.1066964E-07, 0.3889165E-07, 0.1402413E-06, 0.5080287E-06, &
         0.1838745E-05, 0.6654002E-05, 0.2408813E-04, 0.8719498E-04, &
         0.3155253E-03, 0.1139886E-02, 0.4094611E-02, 0.1445007E-01, &
         0.3257575E-03, 0.8448067E-14, 0.9466076E-14, 0.6598221E-14, &
         0.5819855E-13, 0.1483693E-12, 0.2723606E-13, 0.1134216E-12, &
         -0.7536919E-12, 0.2585718E-12, 0.3713785E-13, 0.1066277E-11, &
         0.1599707E-12, 0.3999822E-12, 0.3919878E-11, 0.1730125E-10, &
         0.4113725E-10, 0.1999410E-09, 0.6582198E-09, 0.2469145E-08, &
         0.8993994E-08, 0.3211079E-07, 0.1164509E-06, 0.4218663E-06, &
         0.1527693E-05, 0.5526857E-05, 0.2000606E-04, 0.7242112E-04, &
         0.2620910E-03, 0.9472192E-03, 0.3353513E-01, 0.3950379E-04, &
         0.1429898E-03, 0.5172197E-03, 0.1865839E-02, 0.6671186E-02, &
         0.2322445E-01, 0.7551139E-01, 0.2170407E+00, 0.5241100E+00, &
         0.8542780E+00, 0.1208818E+01, 0.1538107E+01, 0.1891108E+01, &
         0.2099686E+01, 0.1850941E+01, 0.1153786E+01, 0.7584199E+00, &
         0.1408991E+01, 0.1537233E+01, 0.1201549E+01, 0.1408215E+01, &
         0.1379510E+01, 0.1359040E+01, 0.1355852E+01, 0.1382048E+01, &
         0.1385486E+01, 0.1393122E+01, 0.1399668E+01, 0.1368301E+01, &
         0.8438055E-05, 0.2306761E-04, 0.8350273E-04, 0.3021766E-03, &
         0.1091762E-02, 0.3922987E-02, 0.1385760E-01, 0.4668255E-01, &
         0.1419962E+00, 0.3740283E+00, 0.7350657E+00, 0.1032138E+01, &
         0.1392090E+01, 0.1722258E+01, 0.1978830E+01, 0.2038919E+01, &
         0.1466292E+01, 0.8217813E+00, 0.1068222E+01, 0.1678929E+01, &
         0.1179439E+01, 0.1465021E+01, 0.1322734E+01, 0.1388140E+01, &
         0.1398959E+01, 0.1402826E+01, 0.1374126E+01, 0.1401267E+01, &
         0.1384118E+01, 0.1380084E+01, 0.6045238E-05, 0.1652407E-04, &
         0.5981760E-04, 0.2164966E-03, 0.7826978E-03, 0.2818324E-02, &
         0.1001874E-01, 0.3432159E-01, 0.1077876E+00, 0.2957576E+00, &
         0.6457256E+00, 0.9422171E+00, 0.1321832E+01, 0.1654439E+01, &
         0.1947359E+01, 0.2037713E+01, 0.1652677E+01, 0.9366399E+00, &
         0.8927257E+00, 0.1617906E+01, 0.1334269E+01, 0.1372657E+01, &
         0.1323403E+01, 0.1372922E+01, 0.1323401E+01, 0.1366046E+01, &
         0.1377430E+01, 0.1352949E+01, 0.1383479E+01, 0.1394892E+01/),&
         (/30,1,6/))

    qbrqs(1:30,1:1,13:18)=reshape( (/&
         0.3857907E-05, 0.1054521E-04, 0.3817201E-04, 0.1381718E-03, &
         0.4998088E-03, 0.1803249E-02, 0.6449907E-02, 0.2247932E-01, &
         0.7327781E-01, 0.2114086E+00, 0.5139958E+00, 0.8473598E+00, &
         0.1197895E+01, 0.1527089E+01, 0.1880320E+01, 0.2101947E+01, &
         0.1867903E+01, 0.1170279E+01, 0.7712408E+00, 0.1393570E+01, &
         0.1562685E+01, 0.1217631E+01, 0.1419034E+01, 0.1369893E+01, &
         0.1385373E+01, 0.1366546E+01, 0.1356363E+01, 0.1383033E+01, &
         0.1387707E+01, 0.1393639E+01, 0.2875744E-05, 0.7860280E-05, &
         0.2845365E-04, 0.1030003E-03, 0.3726821E-03, 0.1345828E-02, &
         0.4827840E-02, 0.1696944E-01, 0.5645249E-01, 0.1680575E+00, &
         0.4297750E+00, 0.7848734E+00, 0.1099823E+01, 0.1440201E+01, &
         0.1776021E+01, 0.2050056E+01, 0.1986646E+01, 0.1322586E+01, &
         0.7947940E+00, 0.1190645E+01, 0.1659948E+01, 0.1198948E+01, &
         0.1479658E+01, 0.1347349E+01, 0.1393301E+01, 0.1416820E+01, &
         0.1365937E+01, 0.1370853E+01, 0.1368950E+01, 0.1383935E+01, &
         0.1718626E-05, 0.4697151E-05, 0.1700411E-04, 0.6155527E-04, &
         0.2227757E-03, 0.8053585E-03, 0.2899469E-02, 0.1030230E-01, &
         0.3524727E-01, 0.1104061E+00, 0.3019377E+00, 0.6538075E+00, &
         0.9490179E+00, 0.1328286E+01, 0.1660758E+01, 0.1947923E+01, &
         0.2039193E+01, 0.1638955E+01, 0.9184783E+00, 0.8948992E+00, &
         0.1630617E+01, 0.1313450E+01, 0.1395140E+01, 0.1328123E+01, &
         0.1377631E+01, 0.1357848E+01, 0.1388879E+01, 0.1391962E+01, &
         0.1361340E+01, 0.1369491E+01, 0.7404256E-06, 0.2023248E-05, &
         0.7321936E-05, 0.2650530E-04, 0.9594424E-04, 0.3471650E-03, &
         0.1253919E-02, 0.4500857E-02, 0.1584823E-01, 0.5295629E-01, &
         0.1588167E+00, 0.4104221E+00, 0.7685978E+00, 0.1076261E+01, &
         0.1422771E+01, 0.1755290E+01, 0.2024872E+01, 0.2010120E+01, &
         0.1377721E+01, 0.8087896E+00, 0.1148468E+01, 0.1692759E+01, &
         0.1168727E+01, 0.1450756E+01, 0.1345646E+01, 0.1390995E+01, &
         0.1421059E+01, 0.1388422E+01, 0.1380926E+01, 0.1383566E+01, &
         0.5598029E-06, 0.1530578E-05, 0.5538715E-05, 0.2004728E-04, &
         0.7256944E-04, 0.2626222E-03, 0.9491406E-03, 0.3413769E-02, &
         0.1209360E-01, 0.4104747E-01, 0.1265943E+00, 0.3394618E+00, &
         0.6988942E+00, 0.9913729E+00, 0.1363254E+01, 0.1694228E+01, &
         0.1954880E+01, 0.2044281E+01, 0.1563253E+01, 0.8779636E+00, &
         0.9851181E+00, 0.1671286E+01, 0.1245029E+01, 0.1393614E+01, &
         0.1312447E+01, 0.1383697E+01, 0.1368982E+01, 0.1410644E+01, &
         0.1386400E+01, 0.1377081E+01, 0.3768961E-06, 0.1029766E-05, &
         0.3728270E-05, 0.1349196E-04, 0.4884070E-04, 0.1767822E-03, &
         0.6393095E-03, 0.2304285E-02, 0.8216626E-02, 0.2838739E-01, &
         0.9074473E-01, 0.2547902E+00, 0.5869321E+00, 0.8975460E+00, &
         0.1271339E+01, 0.1603128E+01, 0.1935510E+01, 0.2057565E+01, &
         0.1761917E+01, 0.1026894E+01, 0.8074358E+00, 0.1538698E+01, &
         0.1437847E+01, 0.1310425E+01, 0.1362680E+01, 0.1381976E+01, &
         0.1356858E+01, 0.1333898E+01, 0.1381830E+01, 0.1388422E+01/),&
         (/30,1,6/))

    qbrqs(1:30,1:1,19:24)=reshape( (/&
         0.3236136E-06, 0.8855458E-06, 0.3204303E-05, 0.1159740E-04, &
         0.4198187E-04, 0.1519609E-03, 0.5496361E-03, 0.1982331E-02, &
         0.7082549E-02, 0.2460569E-01, 0.7962730E-01, 0.2273520E+00, &
         0.5421261E+00, 0.8665863E+00, 0.1227701E+01, 0.1557529E+01, &
         0.1907946E+01, 0.2091759E+01, 0.1820815E+01, 0.1128173E+01, &
         0.7792270E+00, 0.1441805E+01, 0.1500771E+01, 0.1177203E+01, &
         0.1397181E+01, 0.1378007E+01, 0.1360600E+01, 0.1342718E+01, &
         0.1381574E+01, 0.1394904E+01, 0.2598837E-06, 0.7105884E-06, &
         0.2571407E-05, 0.9306488E-05, 0.3368797E-04, 0.1219437E-03, &
         0.4411651E-03, 0.1592339E-02, 0.5703115E-02, 0.1995317E-01, &
         0.6563381E-01, 0.1919264E+00, 0.4775247E+00, 0.8215518E+00, &
         0.1156758E+01, 0.1487721E+01, 0.1835645E+01, 0.2094382E+01, &
         0.1909802E+01, 0.1222000E+01, 0.7729526E+00, 0.1325953E+01, &
         0.1618831E+01, 0.1129606E+01, 0.1441841E+01, 0.1368787E+01, &
         0.1390773E+01, 0.1396055E+01, 0.1353166E+01, 0.1377240E+01, &
         0.1468962E-06, 0.4008892E-06, 0.1451859E-05, 0.5252375E-05, &
         0.1901069E-04, 0.6881658E-04, 0.2490552E-03, 0.9001884E-03, &
         0.3238772E-02, 0.1148521E-01, 0.3908658E-01, 0.1211622E+00, &
         0.3270030E+00, 0.6846489E+00, 0.9770864E+00, 0.1352290E+01, &
         0.1683818E+01, 0.1951200E+01, 0.2043353E+01, 0.1586950E+01, &
         0.8647165E+00, 0.9542888E+00, 0.1668071E+01, 0.1251961E+01, &
         0.1387784E+01, 0.1315097E+01, 0.1377233E+01, 0.1364871E+01, &
         0.1410496E+01, 0.1380941E+01, 0.8878270E-07, 0.2427971E-06, &
         0.8784166E-06, 0.3178303E-05, 0.1150401E-04, 0.4164483E-04, &
         0.1507380E-03, 0.5452277E-03, 0.1966480E-02, 0.7026603E-02, &
         0.2441811E-01, 0.7907017E-01, 0.2259613E+00, 0.5397338E+00, &
         0.8649546E+00, 0.1225237E+01, 0.1554972E+01, 0.1905898E+01, &
         0.2093079E+01, 0.1824446E+01, 0.1131443E+01, 0.7718529E+00, &
         0.1433445E+01, 0.1505972E+01, 0.1179354E+01, 0.1397423E+01, &
         0.1375723E+01, 0.1354594E+01, 0.1339999E+01, 0.1368543E+01, &
         0.4909048E-07, 0.1342343E-06, 0.4857623E-06, 0.1758437E-05, &
         0.6363671E-05, 0.2303496E-04, 0.8338189E-04, 0.3017350E-03, &
         0.1090177E-02, 0.3917343E-02, 0.1383809E-01, 0.4662061E-01, &
         0.1418286E+00, 0.3736582E+00, 0.7347008E+00, 0.1031695E+01, &
         0.1391785E+01, 0.1721950E+01, 0.1978449E+01, 0.2039113E+01, &
         0.1467247E+01, 0.8228753E+00, 0.1067643E+01, 0.1679000E+01, &
         0.1180979E+01, 0.1465477E+01, 0.1323149E+01, 0.1388206E+01, &
         0.1398749E+01, 0.1403252E+01, 0.4491110E-07, 0.1226495E-06, &
         0.4438543E-06, 0.1606513E-05, 0.5813939E-05, 0.2104629E-04, &
         0.7618610E-04, 0.2757049E-03, 0.9963178E-03, 0.3582313E-02, &
         0.1267850E-01, 0.4292409E-01, 0.1317570E+00, 0.3511752E+00, &
         0.7116793E+00, 0.1005011E+01, 0.1373192E+01, 0.1703715E+01, &
         0.1960462E+01, 0.2044763E+01, 0.1533077E+01, 0.8851303E+00, &
         0.1014989E+01, 0.1669717E+01, 0.1245077E+01, 0.1428272E+01, &
         0.1312307E+01, 0.1375859E+01, 0.1382066E+01, 0.1413976E+01/),&
         (/30,1,6/))

    qbrqs(1:30,1:1,25:30)=reshape( (/&
         0.3771462E-07, 0.1028282E-06, 0.3719264E-06, 0.1346401E-05, &
         0.4873332E-05, 0.1763927E-04, 0.6385347E-04, 0.2310964E-03, &
         0.8353932E-03, 0.3006974E-02, 0.1067756E-01, 0.3646917E-01, &
         0.1138473E+00, 0.3100131E+00, 0.6640764E+00, 0.9579683E+00, &
         0.1336368E+01, 0.1668591E+01, 0.1948658E+01, 0.2041109E+01, &
         0.1620370E+01, 0.8906084E+00, 0.9080491E+00, 0.1650102E+01, &
         0.1287251E+01, 0.1404940E+01, 0.1312431E+01, 0.1369794E+01, &
         0.1347291E+01, 0.1383594E+01, 0.2066275E-07, 0.5672286E-07, &
         0.2053315E-06, 0.7420299E-06, 0.2685598E-05, 0.9722087E-05, &
         0.3519086E-04, 0.1273800E-03, 0.4608098E-03, 0.1663016E-02, &
         0.5953574E-02, 0.2080238E-01, 0.6821626E-01, 0.1985453E+00, &
         0.4901721E+00, 0.8306453E+00, 0.1171311E+01, 0.1501209E+01, &
         0.1851867E+01, 0.2099972E+01, 0.1896659E+01, 0.1208096E+01, &
         0.8045296E+00, 0.1344672E+01, 0.1607270E+01, 0.1197700E+01, &
         0.1446757E+01, 0.1368162E+01, 0.1369020E+01, 0.1384346E+01, &
         0.2943546E-07, 0.8037221E-07, 0.2906546E-06, 0.1051397E-05, &
         0.3806524E-05, 0.1377661E-04, 0.4987255E-04, 0.1805111E-03, &
         0.6527653E-03, 0.2352573E-02, 0.8386328E-02, 0.2895008E-01, &
         0.9238102E-01, 0.2587797E+00, 0.5930740E+00, 0.9019188E+00, &
         0.1276954E+01, 0.1608966E+01, 0.1937796E+01, 0.2052794E+01, &
         0.1750075E+01, 0.1011107E+01, 0.8058437E+00, 0.1549847E+01, &
         0.1425072E+01, 0.1307907E+01, 0.1361917E+01, 0.1381650E+01, &
         0.1353116E+01, 0.1336342E+01, 0.1401500E-07, 0.3835394E-07, &
         0.1390823E-06, 0.5033391E-06, 0.1820643E-05, 0.6589878E-05, &
         0.2385175E-04, 0.8634376E-04, 0.3124391E-03, 0.1128766E-02, &
         0.4054963E-02, 0.1431323E-01, 0.4812629E-01, 0.1458958E+00, &
         0.3825957E+00, 0.7433397E+00, 0.1042454E+01, 0.1399161E+01, &
         0.1729491E+01, 0.1988325E+01, 0.2033609E+01, 0.1445695E+01, &
         0.8003992E+00, 0.1082630E+01, 0.1693802E+01, 0.1152685E+01, &
         0.1457584E+01, 0.1329207E+01, 0.1388109E+01, 0.1413331E+01, &
         0.9938845E-08, 0.2700967E-07, 0.9789189E-07, 0.3541029E-06, &
         0.1281403E-05, 0.4637087E-05, 0.1678571E-04, 0.6076441E-04, &
         0.2199208E-03, 0.7950518E-03, 0.2862566E-02, 0.1017336E-01, &
         0.3482668E-01, 0.1092176E+00, 0.2991365E+00, 0.6501684E+00, &
         0.9459308E+00, 0.1325391E+01, 0.1657931E+01, 0.1947675E+01, &
         0.2038492E+01, 0.1645274E+01, 0.9275196E+00, 0.8941895E+00, &
         0.1624768E+01, 0.1322001E+01, 0.1387358E+01, 0.1331443E+01, &
         0.1374707E+01, 0.1341997E+01, 0.6412746E-08, 0.1769780E-07, &
         0.6366894E-07, 0.2303611E-06, 0.8336778E-06, 0.3016920E-05, &
         0.1092051E-04, 0.3952926E-04, 0.1430871E-03, 0.5175775E-03, &
         0.1867124E-02, 0.6675717E-02, 0.2323973E-01, 0.7555697E-01, &
         0.2171555E+00, 0.5243142E+00, 0.8544180E+00, 0.1209036E+01, &
         0.1538329E+01, 0.1891320E+01, 0.2099623E+01, 0.1850582E+01, &
         0.1153475E+01, 0.7582853E+00, 0.1409339E+01, 0.1536727E+01, &
         0.1201235E+01, 0.1408159E+01, 0.1379772E+01, 0.1358475E+01/),&
         (/30,1,6/))

    qbrqs(1:30,1:1,31:36)=reshape( (/&
         0.3923926E-08, 0.1067798E-07, 0.3897862E-07, 0.1409177E-06, &
         0.5099691E-06, 0.1845879E-05, 0.6680979E-05, 0.2418180E-04, &
         0.8753642E-04, 0.3167605E-03, 0.1144346E-02, 0.4110507E-02, &
         0.1450484E-01, 0.4873199E-01, 0.1475266E+00, 0.3861555E+00, &
         0.7467078E+00, 0.1046755E+01, 0.1402116E+01, 0.1732574E+01, &
         0.1992608E+01, 0.2031040E+01, 0.1437798E+01, 0.7941419E+00, &
         0.1091892E+01, 0.1695350E+01, 0.1146589E+01, 0.1452048E+01, &
         0.1327341E+01, 0.1386919E+01, 0.2159660E-08, 0.5853871E-08, &
         0.2139184E-07, 0.7725549E-07, 0.2794547E-06, 0.1011294E-05, &
         0.3659937E-05, 0.1324533E-04, 0.4794678E-04, 0.1735502E-03, &
         0.6276303E-03, 0.2262394E-02, 0.8069280E-02, 0.2789819E-01, &
         0.8931855E-01, 0.2513025E+00, 0.5814862E+00, 0.8937050E+00, &
         0.1266289E+01, 0.1597861E+01, 0.1933207E+01, 0.2062025E+01, &
         0.1771295E+01, 0.1042512E+01, 0.8124586E+00, 0.1524010E+01, &
         0.1449887E+01, 0.1308147E+01, 0.1367361E+01, 0.1380090E+01, &
         0.7260568E-09, 0.2020068E-08, 0.7262182E-08, 0.2583743E-07, &
         0.9367815E-07, 0.3389539E-06, 0.1226493E-05, 0.4437584E-05, &
         0.1606419E-04, 0.5814762E-04, 0.2104592E-03, 0.7609050E-03, &
         0.2740245E-02, 0.9745683E-02, 0.3342814E-01, 0.1052505E+00, &
         0.2897407E+00, 0.6376665E+00, 0.9356289E+00, 0.1315299E+01, &
         0.1647984E+01, 0.1946704E+01, 0.2036762E+01, 0.1665433E+01, &
         0.9483317E+00, 0.8870602E+00, 0.1611704E+01, 0.1349355E+01, &
         0.1338266E+01, 0.1332024E+01, 0.2729792E-09, 0.7788261E-09, &
         0.2781656E-08, 0.1002817E-07, 0.3617201E-07, 0.1306958E-06, &
         0.4733878E-06, 0.1714013E-05, 0.6202667E-05, 0.2245028E-04, &
         0.8127029E-04, 0.2940956E-03, 0.1062638E-02, 0.3819079E-02, &
         0.1349840E-01, 0.4554095E-01, 0.1388998E+00, 0.3671713E+00, &
         0.7282541E+00, 0.1023936E+01, 0.1386446E+01, 0.1716618E+01, &
         0.1972223E+01, 0.2042006E+01, 0.1485218E+01, 0.8435971E+00, &
         0.1052239E+01, 0.1682687E+01, 0.1208819E+01, 0.1468988E+01, &
         0.1138483E-09, 0.3352924E-09, 0.1323467E-08, 0.4765051E-08, &
         0.1696672E-07, 0.6176972E-07, 0.2235700E-06, 0.8091876E-06, &
         0.2928622E-05, 0.1059961E-04, 0.3837255E-04, 0.1388960E-03, &
         0.5024281E-03, 0.1812669E-02, 0.6483227E-02, 0.2259160E-01, &
         0.7361497E-01, 0.2122605E+00, 0.5155379E+00, 0.8483830E+00, &
         0.1199574E+01, 0.1528771E+01, 0.1882021E+01, 0.2101714E+01, &
         0.1865457E+01, 0.1167623E+01, 0.7686450E+00, 0.1395737E+01, &
         0.1559300E+01, 0.1215185E+01, 0.6871950E-10, 0.1435970E-09, &
         0.5198381E-09, 0.1935983E-08, 0.7118665E-08, 0.2559486E-07, &
         0.9233558E-07, 0.3342500E-06, 0.1209471E-05, 0.4377391E-05, &
         0.1584374E-04, 0.5735440E-04, 0.2075856E-03, 0.7505263E-03, &
         0.2703073E-02, 0.9615588E-02, 0.3300178E-01, 0.1040363E+00, &
         0.2868509E+00, 0.6337283E+00, 0.9324754E+00, 0.1312067E+01, &
         0.1644767E+01, 0.1946319E+01, 0.2036625E+01, 0.1671548E+01, &
         0.9529492E+00, 0.8828400E+00, 0.1606683E+01, 0.1356629E+01/),&
         (/30,1,6/))

    qbrqs(1:30,1:1,37:42)=reshape( (/&
         0.4279296E-10, 0.1019830E-09, 0.3495865E-09, 0.1240186E-08, &
         0.4548052E-08, 0.1637553E-07, 0.5894112E-07, 0.2132708E-06, &
         0.7721861E-06, 0.2795869E-05, 0.1012025E-04, 0.3663430E-04, &
         0.1326112E-03, 0.4797276E-03, 0.1731045E-02, 0.6194473E-02, &
         0.2161726E-01, 0.7068184E-01, 0.2048290E+00, 0.5019357E+00, &
         0.8389409E+00, 0.1184580E+01, 0.1513944E+01, 0.1866371E+01, &
         0.2102310E+01, 0.1884416E+01, 0.1191447E+01, 0.7944575E+00, &
         0.1358286E+01, 0.1585731E+01, 0.1320506E-10, 0.3754272E-10, &
         0.1172757E-09, 0.4610952E-09, 0.1747573E-08, 0.6396914E-08, &
         0.2307197E-07, 0.8355009E-07, 0.3025998E-06, 0.1094051E-05, &
         0.3960026E-05, 0.1433091E-04, 0.5187821E-04, 0.1877786E-03, &
         0.6790176E-03, 0.2446735E-02, 0.8717073E-02, 0.3004434E-01, &
         0.9555029E-01, 0.2664733E+00, 0.6046580E+00, 0.9103130E+00, &
         0.1287310E+01, 0.1619669E+01, 0.1941297E+01, 0.2045045E+01, &
         0.1725672E+01, 0.9880461E+00, 0.8178003E+00, 0.1567536E+01, &
         0.5899305E-11, 0.2042627E-10, 0.3764592E-10, 0.1280842E-09, &
         0.4243975E-09, 0.1496634E-08, 0.5502207E-08, 0.1978670E-07, &
         0.7175582E-07, 0.2595964E-06, 0.9389952E-06, 0.3399184E-05, &
         0.1230400E-04, 0.4453738E-04, 0.1612071E-03, 0.5830508E-03, &
         0.2102342E-02, 0.7505785E-02, 0.2602144E-01, 0.8381479E-01, &
         0.2377558E+00, 0.5596477E+00, 0.8785649E+00, 0.1245335E+01, &
         0.1575926E+01, 0.1921143E+01, 0.2079966E+01, 0.1798919E+01, &
         0.1101398E+01, 0.8241425E+00, 0.1203991E-11, 0.5900688E-11, &
         0.1824277E-10, 0.5389456E-10, 0.1671786E-09, 0.5994772E-09, &
         0.2299403E-08, 0.8175137E-08, 0.2941320E-07, 0.1066256E-06, &
         0.3856113E-06, 0.1395808E-05, 0.5051475E-05, 0.1828509E-04, &
         0.6619136E-04, 0.2395459E-03, 0.8658932E-03, 0.3116101E-02, &
         0.1105803E-01, 0.3770426E-01, 0.1173082E+00, 0.3180817E+00, &
         0.6740167E+00, 0.9669973E+00, 0.1344094E+01, 0.1676006E+01, &
         0.1949607E+01, 0.2042465E+01, 0.1603115E+01, 0.8707204E+00, &
         0.3155373E-11, 0.1247106E-11, 0.9425561E-11, 0.2065244E-10, &
         0.7427725E-10, 0.2928719E-09, 0.1081926E-08, 0.3947230E-08, &
         0.1429128E-07, 0.5198186E-07, 0.1883016E-06, 0.6817618E-06, &
         0.2467739E-05, 0.8929537E-05, 0.3232451E-04, 0.1170097E-03, &
         0.4233310E-03, 0.1528165E-02, 0.5475521E-02, 0.1917978E-01, &
         0.6327049E-01, 0.1858341E+00, 0.4656517E+00, 0.8128141E+00, &
         0.1142859E+01, 0.1475338E+01, 0.1820214E+01, 0.2086410E+01, &
         0.1925889E+01, 0.1240782E+01, 0.1531637E-11, 0.2570500E-11, &
         0.1066159E-11, 0.1496788E-10, 0.5043575E-10, 0.1727316E-09, &
         0.6206500E-09, 0.2259410E-08, 0.8238917E-08, 0.2965713E-07, &
         0.1075832E-06, 0.3893423E-06, 0.1409170E-05, 0.5098709E-05, &
         0.1845820E-04, 0.6681245E-04, 0.2417935E-03, 0.8740014E-03, &
         0.3145104E-02, 0.1115909E-01, 0.3803171E-01, 0.1182230E+00, &
         0.3202052E+00, 0.6765814E+00, 0.9693897E+00, 0.1346075E+01, &
         0.1677899E+01, 0.1949971E+01, 0.2042721E+01, 0.1599026E+01/),&
         (/30,1,6/))

    qbrqs(1:30,1:1,43:48)=reshape( (/&
         0.9917659E-12, 0.1381379E-11, 0.3523106E-11, 0.9577097E-11, &
         0.3459986E-10, 0.1411091E-09, 0.4419279E-09, 0.1553517E-08, &
         0.5765216E-08, 0.2065135E-07, 0.7509527E-07, 0.2717656E-06, &
         0.9833307E-06, 0.3558670E-05, 0.1288058E-04, 0.4662503E-04, &
         0.1687613E-03, 0.6103279E-03, 0.2200291E-02, 0.7850773E-02, &
         0.2717153E-01, 0.8719374E-01, 0.2460894E+00, 0.5732117E+00, &
         0.8879268E+00, 0.1258482E+01, 0.1589697E+01, 0.1929189E+01, &
         0.2068972E+01, 0.1783306E+01, 0.3368454E-12, 0.1753294E-11, &
         0.1704585E-12, 0.3925133E-11, 0.2875746E-10, 0.8277545E-10, &
         0.3136528E-09, 0.1162799E-08, 0.4243529E-08, 0.1521952E-07, &
         0.5533235E-07, 0.1997696E-06, 0.7227579E-06, 0.2614614E-05, &
         0.9465122E-05, 0.3426161E-04, 0.1240174E-03, 0.4486567E-03, &
         0.1619298E-02, 0.5798678E-02, 0.2027742E-01, 0.6662144E-01, &
         0.1944622E+00, 0.4824011E+00, 0.8250818E+00, 0.1162402E+01, &
         0.1492890E+01, 0.1841950E+01, 0.2096908E+01, 0.1904410E+01, &
         -0.1283701E-11, 0.2700828E-12, 0.5694661E-11, 0.1144088E-11, &
         0.1337860E-10, 0.5726562E-10, 0.2234361E-09, 0.8505050E-09, &
         0.3057868E-08, 0.1090547E-07, 0.3982935E-07, 0.1442852E-06, &
         0.5222154E-06, 0.1889433E-05, 0.6839028E-05, 0.2475330E-04, &
         0.8960722E-04, 0.3242505E-03, 0.1171335E-02, 0.4206709E-02, &
         0.1483641E-01, 0.4977831E-01, 0.1503363E+00, 0.3922572E+00, &
         0.7523814E+00, 0.1054152E+01, 0.1407214E+01, 0.1737990E+01, &
         0.2000353E+01, 0.2026293E+01, 0.1193975E-11, 0.1670941E-11, &
         0.5259805E-11, 0.4074115E-11, 0.2072773E-10, 0.5170789E-10, &
         0.1838727E-09, 0.6078937E-09, 0.2173212E-08, 0.7777988E-08, &
         0.2828418E-07, 0.1023848E-06, 0.3703800E-06, 0.1340259E-05, &
         0.4847633E-05, 0.1754610E-04, 0.6351911E-04, 0.2298873E-03, &
         0.8310285E-03, 0.2991342E-02, 0.1062301E-01, 0.3629179E-01, &
         0.1133487E+00, 0.3088466E+00, 0.6626130E+00, 0.9566727E+00, &
         0.1335223E+01, 0.1667487E+01, 0.1948545E+01, 0.2040859E+01, &
         -0.1518614E-11, 0.1025309E-11, 0.2377287E-12, 0.8970014E-11, &
         0.3335786E-11, 0.3682587E-10, 0.1141844E-09, 0.4113278E-09, &
         0.1523402E-08, 0.5454212E-08, 0.1970126E-07, 0.7100462E-07, &
         0.2571013E-06, 0.9296294E-06, 0.3366042E-05, 0.1218411E-04, &
         0.4409765E-04, 0.1596165E-03, 0.5773018E-03, 0.2081709E-02, &
         0.7433035E-02, 0.2577847E-01, 0.8309824E-01, 0.2359817E+00, &
         0.5567063E+00, 0.8765479E+00, 0.1242425E+01, 0.1572880E+01, &
         0.1919152E+01, 0.2082179E+01, 0.3449489E-12, 0.1323594E-12, &
         0.2524434E-11, 0.7129730E-11, 0.1099827E-10, 0.2720793E-10, &
         0.9374397E-10, 0.2881573E-09, 0.1014721E-08, 0.3635337E-08, &
         0.1330892E-07, 0.4817363E-07, 0.1743618E-06, 0.6308808E-06, &
         0.2282745E-05, 0.8262761E-05, 0.2990859E-04, 0.1082676E-03, &
         0.3917242E-03, 0.1414399E-02, 0.5071576E-02, 0.1780289E-01, &
         0.5903469E-01, 0.1748266E+00, 0.4436502E+00, 0.7959911E+00, &
         0.1116611E+01, 0.1453352E+01, 0.1792378E+01, 0.2065975E+01/),&
         (/30,1,6/))

    qbrqs(1:30,1:1,49:50)=reshape( (/&
         -0.9806707E-13, 0.1529088E-11, 0.2389664E-12, 0.2416812E-11, &
         -0.3707256E-12, 0.1301249E-10, 0.6405868E-10, 0.1982637E-09, &
         0.7039382E-09, 0.2456551E-08, 0.8758652E-08, 0.3174901E-07, &
         0.1151133E-06, 0.4168501E-06, 0.1507989E-05, 0.5457498E-05, &
         0.1975455E-04, 0.7150998E-04, 0.2587977E-03, 0.9353363E-03, &
         0.3364424E-02, 0.1192217E-01, 0.4049582E-01, 0.1250701E+00, &
         0.3359799E+00, 0.6949816E+00, 0.9873564E+00, 0.1360237E+01, &
         0.1691364E+01, 0.1953644E+01, 0.1310020E-11, 0.7826037E-12, &
         -0.8811858E-12, 0.2869265E-11, 0.2279223E-11, 0.5720607E-12, &
         0.1676336E-10, 0.9887384E-10, 0.3575822E-09, 0.1381762E-08, &
         0.4984234E-08, 0.1798714E-07, 0.6499247E-07, 0.2345360E-06, &
         0.8492439E-06, 0.3073006E-05, 0.1112344E-04, 0.4026809E-04, &
         0.1457573E-03, 0.5272286E-03, 0.1901805E-02, 0.6798243E-02, &
         0.2365158E-01, 0.7678756E-01, 0.2202472E+00, 0.5297817E+00, &
         0.8581598E+00, 0.1214843E+01, 0.1544259E+01, 0.1896759E+01/),&
         (/30,1,2/))

    ! ------------------------------------------------------------------

   
    IF ( i_mie .eq. i_READ ) THEN
       all_ok = (mwave==nwave) .and. (mrad==nrad) .and. (mgroup==ngroup)
       DO ig = 1, ngroup
          !...dbg:
          all_ok = all_ok .and. ( r(1,ig) .eq. rmin(ig) )
       END DO
       !
       IF ( .not. all_ok )THEN
          WRITE(lunoprt,*) ' setuprad: mie.data grid(s) bad: '
          WRITE(lunoprt,*) ' in mie.data, mwave, mrad, mgroup = ', &
               mwave, mrad, mgroup
          WRITE(lunoprt,*) ' and rmin = ', rmin
          stop 1
       END IF
    ELSE
       ! calculate extinction and scattering coefficients
       !
       DO ig = 1,ngroup
          !
          !	select ice/liquid index of refractive index array
          !
          IF( is_grp_ice(ig) )THEN
             irefr = 2
          ELSE
             irefr = 1
          END IF
          !
          !	<thetd> is angle between incident and scattered radiation
          !	<j_thetd> is number of <thetd> values to consider
          !
          thetd = 0.0
          n_thetd = 1
          DO  l=1,nwave
             !
             !kml for biomass burning particles:
             r_real = 1.495
             tmag = 1.0e-6
  
             !	  real = 1.52
             !	  tmag = 0.015
             
             !
             !	 calculate the center of the wavelength interval of an ir interval
             !
             IF( l .le. nsol ) THEN
                awave = wave(l)
             ELSE
                awave = 0.5*(wave(l)+wave(l+1))
             END IF
             wvno     =	2.*pi/(awave*1.0e-4)
             !
             DO i=1,nrad
                IF(i .eq. 1) THEN
                   ddr     = 0.2*(rup(1,ig)-r(1,ig))
                   rr      = r(1,ig)
                   corerad = rcore(1,ig)
                   ddrc    = 0.2*(rcoreup(i,ig)- rcore(1,ig))
                ELSE
                   ddr     = 0.2*(rup(i,ig)-rup(i-1,ig))
                   rr      = rup(i-1,ig)
                   corerad = rcoreup(i-1,ig)
                   ddrc    = 0.2*(rcoreup(i,ig)-rcoreup(i-1,ig))
                END IF
                !
                rdqext(i,ig,l) = 0.0
                qscat(i,ig,l)  = 0.0
                qbrqs(i,ig,l)  = 0.0
                !
                DO j=1,6
                   !
                   !
                   ! limit x=2*pi/wave to no larger 1000 to avoid anguish in the mie code.
                   !
                   IF( wvno*rr .gt. 1000. ) rr = 1000./wvno
                   !  print*,'-------------------------------------------------'
                   !     IF(j.eq.1) rr = 0.1 * 1.e-4
                   !    corerad=0.5*rr
                   !  print*, rr,corerad,rr/corerad
                   !  print*, rr,real,tmag,thetd,n_thetd,qextd,qscatd,ctbrqs,
                   ! 1	       corerad,corereal,coreimag,wvno
                   !  print*,'-------------------------------------------------'
                   !      stop
                   CALL miess(rr,r_real,tmag,thetd,n_thetd,qextd,qscatd,ctbrqs,&
                        corerad,corereal,coreimag,wvno)
                   !     IF(l.eq.7 .or. l.eq.8)
                   ! &      print*,'qex qsc=',qextd,qscatd
                   !
                   rdqext(i,ig,l) = rdqext(i,ig,l)+qextd/6.
                   qscat(i,ig,l)  = qscat(i,ig,l)+qscatd/6.
                   qbrqs(i,ig,l)  = qbrqs(i,ig,l)+ctbrqs/6.
                   rr             = rr+ddr
                   corerad        = corerad + ddrc
                   !
                END DO
             END DO
          END DO
          !
          !	 stop
       END DO	  ! ig=1,ngroup
       !
    END IF
  
    !
    !srf - nao precisa escrever o arquivo mie.data
    !	 i_mie=555
    
    !	   DO ig = 1,ngroup
    !	    DO i = 1,nrad
    !	     DO l = 1,nwave
    !	       print*, 'rdqext(i,ig,l), qscat(i,ig,l), qbrqs(i,ig,l)'
    !	       print*, rdqext(i,ig,l), qscat(i,ig,l), qbrqs(i,ig,l)
    !	     END DO
    !	    END DO
    !	   END DO
  
    !IF ( i_mie .eq. i_WRITE ) THEN
    !  PRINT *,'LFR->006.1.4: Writing mie';CALL flush(6)
    !
    !   ! WRITE extinction and scattering coefficients to data file
    !   !
    !   OPEN(lunmie,file='./mie.data',form='formatted')
    !   !
    !   PRINT *,'LFR->006.1.4.1: Writing nwave,nrad,ngroup',nwave,nrad,ngroup;CALL flush(6)
    !   WRITE(lunmie,*) nwave,nrad,ngroup
    !   PRINT *,'LFR->006.1.4.1: Writing r';CALL flush(6)
    !   DO ig = 1,ngroup
    !	 DO i = 1,nrad
    !	  WRITE(lunmie,*) r(i,ig)
    !	 END DO
    !    END DO
    !
    !    PRINT *,'LFR->006.1.4.2: Writing rdqext,qscat,qbrqs';CALL flush(6)
    !    DO ig = 1,ngroup
    !	 DO i = 1,nrad
    !	  DO l = 1,nwave
    !	    WRITE(lunmie,*) rdqext(i,ig,l), qscat(i,ig,l), qbrqs(i,ig,l)
    !	  END DO
    !	 END DO
    !    END DO
    !
    !    CLOSE(lunmie)
    !    PRINT *,'LFR->006.1.4.3: Mie close';CALL flush(6)
    !
    !  END IF
    !	 stop
    !
    !	WRITE some values to print file
    !
    !	 DO ig = 1, ngroup
    !	   WRITE (lunoprt,500) ig
    !	   DO i = 1, nrad
    !	     sizparm6 = 2.*pi*r(i,ig)/(wave(6)*1.e-4)
    !	     sizparm24 = 2.*pi*r(i,ig)/(wave(24)*1.e-4)
    !	     WRITE(lunoprt,505) i,r(i,ig),rdqext(i,ig,6),sizparm6,
    !	1	rdqext(i,ig,24),sizparm24
    !	   END DO
    !	 END DO
    !
    !500  FORMAT(/," setuprad: igroup = ",i4,//,"   i	 r(cm)   rdqext(6)	x6","	rdqext(24)	x24",/)
    !505  FORMAT(i4,5(1pe11.2))
    !
    ! *********************************************************************
    !
    !				 check sum of weights
    !
    ! **********************************************************************
    !
    sum  = 0.0
    sum1 = 0.0
    sum2 = 0.0
    DO l = 1,nsolp
       sum  = sum+weight(l)
    END DO
    DO l = nsolp+1,ntotal
       sum1 = sum1+weight(l)
    END DO
    sum2	=   sum+sum1
  
    !
    IF ( abs(nwave-sum2) .gt. 1.e-3 ) WRITE(lunoprt,FMT=lab355) sum,sum1,sum2
    !
    DO l = 1,nsolp
       sol(l) = solfx(nprob(l)) * weight(l)
    END DO
    !	 print*, 'wave(l),nprob(l),weight(l),solfx(nprob(l)),sol(l)'
    !	 DO 361 l   =	1,ntotal
    !
    !	   print*,wave(nprob(l)),nprob(l),weight(nprob(l)),solfx(nprob(l)),sol(l)
    !
    ! 361  continue
  
    !
    ! *********************************************************************
    !
    !	compute planck function table. wave is in units of microns.
    !
    ! **********************************************************************
    !
    !	set <iblackbody_above> = 1 to include a source of radiation
    !	at the top of the radiative transfer model DOmain
    
    iblackbody_above = ir_above_aerad
    t_above = tabove_aerad
    !
    !	set <ibeyond_spectrum> = 1 to include blackbody radiation at
    !	wavelengths longer than wave(nwave+1) in plank(nwave+1-nsol) and
    !	at wavelengths shorter than wave(nsol+1) in plank(1)
    !
    ibeyond_spectrum = 1
    !
    IF( ibeyond_spectrum .eq. 1 )THEN 
       DO j = 1,ncount
          plank(nwave+1-nsol,j) = (0.01*float(nlow+j))**4
       END DO
       DO i = nsol+2,nwave
          DO j = 1,ncount
             k = i-nsol
             v = 1.438e4 / wave(i)
             CALL plnk(v,(0.01*float(nlow+j)),plank(k,j))
          END DO
       END DO
    ELSE 
       DO i = nsol+1,nwave+1
          DO j = 1,ncount
             k = i-nsol
             v = 1.438e4 / wave(i)
             CALL plnk(v,(0.01*float(nlow+j)),plank(k,j))
          END DO
       END DO
    END IF
    !
    DO j = 1,ncount
  
       IF( ibeyond_spectrum .eq. 1 )THEN
  
          plank(1,j) = plank(2,j)*sbk/pi
          DO l = nsol+2,nwave
             k = l-nsol
             plank(k,j) = (plank(k+1,j)-plank(k,j))*sbk/pi
          END DO
  
       ELSE
  
          DO l = nsol+1,nwave
             k = l-nsol
             plank(k,j) = (plank(k+1,j)-plank(k,j))*sbk/pi
          END DO
        
       END IF
  
    END DO
    !
  
  END SUBROUTINE calcproperties
  
!KML2!!!!!!!!1
  SUBROUTINE nocalcproperties
    ! **********************************************************************
    !
    !		 CALCULATE THE AEROSOL EXTINCTION CROSS SECTIONS
    !
    ! **********************************************************************
    !
    !	  Get <is_grp_ice> and radius grid from interface common block
    !	  and calculate cross-sectional area for each bin.
    !
    USE mem_aerad, ONLY: ir_above_aerad,tabove_aerad,lunoprt


  
    USE mem_globrad, ONLY: ngroup,nrad,xsecta,pi, &
  			   i_write,i_read,nwave,rmin, &
  			   nsol,wave,nsolp, &
  			   weight,ntotal,sol,solfx,nprob,iblackbody_above, &
  			   t_above,ncount,nlow,plank,sbk
  
    USE mem_globaer, ONLY: r
  
    IMPLICIT NONE
    
    INTEGER :: i
    INTEGER :: ibeyond_spectrum
    INTEGER :: ig
    INTEGER :: j
    INTEGER :: jj
    INTEGER :: k
    INTEGER :: l
    REAL    :: sum
    REAL    :: sum1
    REAL    :: sum2
    REAL    :: v
  
    CHARACTER(LEN=*),PARAMETER :: &
    lab355='(//,"setuprad: error in weights ",/," ' &
  	   //'sum of weights for solar =",1pe15.5,/,"' &
  	   //' sum of weights for ir = ",1pe15.5,/," ' &
  	   //'total sum =  ",1pe15.5)'
    
    DO ig = 1, ngroup
      DO I = 1, nrad
  	xsecta(i,ig) = pi * r(i,ig)**2.
      END DO
    END DO
   

    !
    ! *********************************************************************
    !
    !				 check sum of weights
    !
    ! **********************************************************************
    !
    sum 	=   0.0
    sum1	=   0.0
    sum2	=   0.0
    DO l	=   1,nsolp
      sum	     =   sum+weight(l)
    END DO
    DO l	=   nsolp+1,ntotal
      sum1	 =   sum1+weight(l)
    END DO
    sum2	=   sum+sum1
  
    !
    IF ( abs(nwave-sum2) .gt. 1.e-3 ) WRITE(lunoprt,FMT=lab355) sum,sum1,sum2
    !
    DO  l   =	1,nsolp
      sol(l)  =   solfx(nprob(l)) * weight(l)
    END DO


    ! *********************************************************************
    !
    !	compute planck function table. wave is in units of microns.
    !
    ! **********************************************************************
    !
    !	set <iblackbody_above> = 1 to include a source of radiation
    !	at the top of the radiative transfer model DOmain
    
    iblackbody_above = ir_above_aerad
    t_above = tabove_aerad
    !
    !	set <ibeyond_spectrum> = 1 to include blackbody radiation at
    !	wavelengths longer than wave(nwave+1) in plank(nwave+1-nsol) and
    !	at wavelengths shorter than wave(nsol+1) in plank(1)
    !
    ibeyond_spectrum = 1
    !
    IF( ibeyond_spectrum .eq. 1 )THEN 
      DO j  =	1,ncount
  	plank(nwave+1-nsol,j) = (0.01*float(nlow+j))**4
      END DO
      DO i =   nsol+2,nwave
        DO j  =	1,ncount
  	  k =	i-nsol
  	  v =	1.438e4  /  wave(i)
  	  CALL plnk(v,(0.01*float(nlow+j)),plank(k,j))
  	END DO
      END DO
    ELSE 
      DO i =   nsol+1,nwave+1
        DO j  =	1,ncount
  	  k =	i-nsol
  	  v =	1.438e4  /  wave(i)
  	  CALL plnk(v,(0.01*float(nlow+j)),plank(k,j))
  	END DO
      END DO
    END IF
    !
    DO j   =   1,ncount
  
      IF( ibeyond_spectrum .eq. 1 )THEN
  
  	plank(1,j) = plank(2,j)*sbk/pi
  	DO l  =   nsol+2,nwave
  	  k  =   l-nsol
  	  plank(k,j) = (plank(k+1,j)-plank(k,j))*sbk/pi
  	END DO
  
      ELSE
  
  	DO l  =   nsol+1,nwave
  	  k  =   l-nsol
  	  plank(k,j) = (plank(k+1,j)-plank(k,j))*sbk/pi
  	END DO
  
      END IF
  
    END DO
    !
  
  END SUBROUTINE nocalcproperties


!KML2!!!!!!!!

  SUBROUTINE prerad(m1,dztr,fmapt,ia,iz,ja,jz,nzpmax,m2,m3)
    !
    !  Note that vertical index in radiative transfer model DOmain is reversed
    !  for cartesian coordinates.
    !
    !  Indices <ix> and <iy> are passed through global COMMON block.
    !
    !  INCLUDE global constants and variables
    !
    
    !USE mem_aerad, ONLY: qv_aerad,pc_aerad
    USE mem_aerad, ONLY: nbin
    USE mem_globaer,ONLY: ngroup,ienconc,nelem
    
    IMPLICIT NONE
  
    INTEGER,INTENT(IN)  	      :: m1,m2,m3,ia,iz,ja,jz,nzpmax
    REAL   ,INTENT(IN), DIMENSION((iz-ia+1)*(jz-ja+1),nzpmax) :: dztr
    REAL   ,INTENT(IN), DIMENSION((iz-ia+1)*(jz-ja+1)) :: fmapt
  
    INTEGER :: ibin
    INTEGER :: i,j,ij,iend
    INTEGER :: iep
    INTEGER :: igas
    INTEGER :: igroup
    INTEGER :: k,kk
    INTEGER :: nzz
    REAL    :: xymet
  
    iend=(iz-ia+1)*(jz-ja+1)
  	   
    !  Load profiles of temperature [K], water vapor [g/g], and
    !  aerosol particle concentrations [#/cm^2]
    !
    igas = 1
    !srf
    nzz = m1 - 1
  
    DO k = 1,NZZ
       !  Reverse the vertical index when in cartesian coordinates
       kk = nzz + 1 - k
       ! For radiation code: qv-aerad have g[H20]/g[ar] unit
       DO ij=1,iend
  	 qv_aerad(ij,kk) = gc(ij,k,igas) / rhoa(ij,k)
       END DO
       DO igroup = 1,ngroup
  	  iep = ienconc(igroup)
  	  DO ibin = 1,nbin
  	    DO ij=1,iend 
  	      xymet = fmapt(ij)*fmapt(ij)
  	      pc_aerad(ij,kk,ibin,igroup) = pc(ij,k,ibin,iep) *  &
  		  (1./dztr(ij,k)) / xymet
  	    END DO
  	  END DO
  	END DO
    END DO
    !
  
  END SUBROUTINE prerad
!kmlnew  
  SUBROUTINE radtran(albedt,cosz,m1,m2,m3,ia,iz,ja,jz,aot11)
!kmlnew     
    USE mem_aerad, ONLY: u0_aerad,qrad_aerad, &
  			 alb_toai_aerad,alb_tomi_aerad,alb_toa_aerad, &
  			 fsl_up_aerad,fsl_dn_aerad,fir_up_aerad,fir_dn_aerad, &
                         nir
      
    USE mem_globrad, ONLY: isl,nvert,nlayer, &
  			   ngroup,nrad,u0,nsolp,albedo_sfc, &
  			   emis,ntotal,emisir,ir,irs, &
  			   fdegday,g,scday,qrad,pi,epsilon,xsecta,rdqext, &
  			   nprob,qscat,nwave,tslu,tsld,fupbs,fdownbs, &
  			   fnetbs,nsol,fslu,fsld,alb_toa,alb_tomi,alb_toai, &
  			   solfx, &
  			   tiru,fupbi,fdownbi,fnetbi,firu,xirup
  
    IMPLICIT NONE
  
    INTEGER,INTENT(IN) :: m1,m2,m3,ia,iz,ja,jz
    REAL,INTENT(IN),DIMENSION((iz-ia+1)*(jz-ja+1)) :: albedt,cosz
    REAL  :: aot11((iz-ia+1)*(jz-ja+1))
    INTEGER :: i,i1,j1
    INTEGER :: ig
    INTEGER :: j
    INTEGER :: l,k
    REAL    :: term1((iz-ia+1)*(jz-ja+1))
    INTEGER :: count=0
    INTEGER :: ij,iend
    REAL :: heati((iz-ia+1)*(jz-ja+1),nlayer)
    REAL :: heats((iz-ia+1)*(jz-ja+1),nlayer)
    REAL :: heat((iz-ia+1)*(jz-ja+1),nlayer)
    
    iend=(iz-ia+1)*(jz-ja+1)
  
    heats   =  0.0
    heati   =  0.0
    
    !
    !	  interpolate temperature from layer center (t) to layer edge (tt)
    DO ij=1,iend
      tt(ij,1) = t_aerad(ij,1)
    END DO
    DO  j = 2, nvert
      DO ij=1,iend
  	  tt(ij,j) = t_aerad(ij,j-1) * (press(ij,j)/p_aerad(ij,j-1)) ** &
  			(LOG(t_aerad(ij,j)/t_aerad(ij,j-1))/ &
  		    LOG(p_aerad(ij,j)/p_aerad(ij,j-1)))
      END DO
    END DO
  
    !	  water vapor (g / cm**2)
      DO  j = 2, nlayer
      DO ij=1,iend
  	  rdh2o(ij,j)	= qv_aerad(ij,j-1) * dpg(ij,j-1)
      END DO
    END DO
  
    !	  aerosol concentrations (# / cm**2)
    DO ig = 1, ngroup
      DO  j = 2, nvert
  	DO  i = 1, nrad
  	  DO ij=1,iend
  	    caer(ij,j,i,ig)  = pc_aerad(ij,j-1,i,ig)  !!!!!!
  	  END DO
  	END DO
      END DO
    END DO  
  
    !surface reflectivity and emissivity
    DO  l =  1,nsolp
      DO ij=1,iend
  	rsfx(ij,l) =  albedt(ij)
  	emis(l) =  0.0
      END DO
    END DO
    DO  l =  nsolp+1,ntotal
      DO ij=1,iend
  	emis(l) =  emisir_aerad
  	rsfx(ij,l) = 1.0 - emis(l)
      END DO
    END DO
    
    !set wavelength limits lla and lls based on values of isl and ir
    lla=  ntotal
    lls=  1
  
    DO ij=1,iend
      IF(isl_aerad(ij)  == 0) THEN
  	lls(ij)   =  nsolp+1
      END IF
    END DO
    !
    IF(ir_aerad   == 0) THEN
      DO ij=1,iend
  	lla(ij)  =  nsolp
      END DO
    END IF
     
      !DO ij=1,iend
      !     print*,'AOT11 na radtran=',ij,aot11
      !END DO
	   
    !calculate the optical properties
    CALL oppr(ia,iz,ja,jz,m1,aot11)
    !
    !	  if infrared calculations are required then calculate
    !	  the plank function
    !
    IF(ir_aerad /= 0) THEN
      CALL oppr1(ia,iz,ja,jz,m1)
    END IF
    !
    !	  if no infrared scattering then set index to number of
    !	  solar intervals
    !
    IF(irs == 0) THEN
      lla  =  nsolp
    END IF
    !
    !	  if either solar or infrared scattering calculations are required
    !	  call the two stream code and find the solution
    !
    
  
    CALL twostr(m1,ia,iz,ja,jz)
    
    !DO i1=ia,iz
    !  DO j1=ja,jz
    !	 IF(isl_aerad(i1,j1) /= 0 .OR. irs .NE. 0 ) THEN
    !	   CALL add(i1,j1,ia,iz,ja,jz,cosz,m2,m3)
    !	 END IF
    !  END DO
    !END DO
    CALL add(m1,ia,iz,ja,jz,cosz,m2,m3)
  
    !
    !	  if infrared calculations are required then call newflux1 for
    !	  a more accurate solution
    !
    IF(ir_aerad /= 0) THEN
      CALL newflux1(m1,ia,iz,ja,jz)
    END IF
    
    !	  calculate infrafred and solar heating rates (deg/day),
    DO  j      =  1,nvert
      DO ij=1,iend
  	IF(isl_aerad(ij) /= 0) THEN
  	  term1(ij)	 =  fdegday/(dpg(ij,j)*g)
  	END IF
      END DO
      DO  l =  1,nsolp
  	DO ij=1,iend
  	  IF(isl_aerad(ij) /= 0) THEN
  	      heats(ij,j)   =  heats(ij,j)+(fnet(ij,l,j+1)-fnet(ij,l,j))*term1(ij)
  	  END IF
  	END DO
      END DO
    END DO
    !  
    DO  j      =  1,nvert
      DO ij=1,iend
  	IF(ir_aerad /= 0) THEN
  	  term1(ij)	 =  fdegday/(dpg(ij,j)*g)
  	END IF
      END DO
      DO  l =  nsolp+1,ntotal
  	DO ij=1,iend
  	  IF(ir_aerad /= 0) THEN
  	    heati(ij,j)  =  heati(ij,j)+(directu(ij,l,j+1)-direc(ij,l,j+1)  &
  		       -(directu(ij,l,j)-direc(ij,l,j)) )*term1(ij)
	  
!srf
!            if(heati(ij,j) < -10) then
!	      print*,ij,j,heati(ij,j)
!        if((ii == 41 .and. jj==36) .or.(ii == 42 .and. jj==36) ) then
!	      if(l==nsolp+1)print*,'radtran ii jj l j',ii,jj,l,j
!	      print*,j,directu(ij,l,j+1),direc(ij,l,j+1)  &
 ! 		       ,directu(ij,l,j),direc(ij,l,j),term1(ij)
!	     endif
!srf		       
  
	  
	  
  	  END IF
  	END DO
      END DO
    END DO
    !
    DO j      =  1,nvert
      DO ij=1,iend
  	!     Load heating rates [deg_K/s] into interface common block
  	heat(ij,j)	   =  heats(ij,j)+heati(ij,j)
  	heats_aerad(ij,j) =  heats(ij,j)/scday
  	heati_aerad(ij,j) =  heati(ij,j)/scday
      END DO
    END DO
    
    !DO  j	=  1,nvert
    ! heats(j)   =  0.0
    !  term1	  =  fdegday/(dpg(j)*g)
    !
    !  IF(isl /= 0) THEN
    !	 DO  l    =  1,nsolp
    !	  heats(j)   =  heats(j)+(fnet(l,j+1)-fnet(l,j))*term1
    !		    print*,l,heats(j)
    !	 END DO
    !  END IF
    !
    !  IF(ir /= 0) THEN
    !	 DO  l    =  nsolp+1,ntotal
    !	  heati(j)  =  heati(j)+( directu(l,j+1)-direc(l,j+1)  &
    !	      -(directu(l,j)-direc(l,j)) )*term1
    !	 END DO
    !  END IF
    !  heat(j)       =  heats(j)+heati(j)
    !
    !	  Load heating rates [deg_K/s] into interface common block
    !
    !  heats_aerad(j) =  heats(j)/scday
    !  heati_aerad(j) =  heati(j)/scday
    !
    !END DO
    !
    !	  Here we Calculate (4 * pi * mean_intensity) for the IR.
    !
    !DO j = 1, nvert
    !  DO l = nsolp+1, ntotal
  !	DO ij=1,iend
  !	  IF (ir_aerad /= 0) THEN
  !	    tmi(ij,l,j) = tmiu(ij,l,j)+tmid(ij,l,j)
  !	  END IF
  !	END DO
  !    END DO
  !  END DO
    !
    !	  Here we compute the heating rates for droplets
    !	  (C11 converts W/m^2 to erg/cm^2)
    !
    !qrad = 0.
    !c11 = 1000.
    !IF (ir /= 0) THEN
    !  DO ig = 1, ngroup
    !	 DO i = 1, nrad
    !	   DO j = 1, nvert
    !	    DO l = nsolp+1, ntotal
    !	      x = tmi(l,j)-4.0*pi*ptemp(l,j)
    !	      IF( ABS(x/tmi(l,j)) < epsilon ) x = 0.
    !	      qrad(i,j,ig) = qrad(i,j,ig) + x*c11*xsecta(i,ig) *  &
    !		  (rdqext(i,ig,nprob(l))-qscat(i,ig,nprob(l)))
    !	    END DO
    !	   END DO  ! j=1,nvert
    !	 END DO   ! i=1,nrad
    !  END DO	    ! ig=1,ngroup
    !END IF
    !
    !IF (isl /= 0) THEN
    !  DO ig = 1, ngroup
    !	 DO i = 1, nrad
    !	   DO j = 1, nvert
    !	    DO l = 1, nsolp
    !	      qrad(i,j,ig) = qrad(i,j,ig) + tmi(l,j)*c11*xsecta(i,ig) *  &
    !		  (rdqext(i,ig,nprob(l))-qscat(i,ig,nprob(l)))
    !	    END DO
    !	   END DO  ! j=1,nvert
    !	 END DO   ! i=1,nrad
    !  END DO	    ! ig=1,ngroup
    !END IF
    !
    !	  Load layer averages of droplet heating rates into interface common block
    !
    !DO ig = 1, ngroup
    !  DO i = 1, nrad
    !	 DO j = 1, nvert
    !	   IF (j == nvert) THEN
    !	    qrad_aerad(i,j,ig) = qrad(i,j,ig)
    !	  ELSE IF (j > 1) THEN
    !	    qrad_aerad(i,j-1,ig) = 0.5 * ( qrad(i,j,ig) + qrad(i,j-1,ig) )
    !	  END IF
    !	 END DO    ! j=1,nvert
    !  END DO	  ! i=1,nrad
    !END DO	  ! ig=1,ngroup
    !
    !
    !	Calculate some diagnostic quantities (formerly done in radout.f) and
    !	load them into the interface common block.  None of these presently
    !	influence any microphysical processes -- hence, the following code only
    !	needs to be executed before the aerosol model writes its output.
    !	Not all of the calculated quantities are presently being
    !	loaded into the interface common block.
    !
    !	Load optical depths into interface common block
    !
    !DO i = 1, nwave
    !  opdaerad(i) = uopd(i,nlayer)
    !	print*,'AOT',uopd(i,nlayer),opd(i,nlayer)
    !END DO
   
    !
    !	  <tsLu> and <tsLd> are total upwelling and downwelling solar
    !	  fluxes at top-of-atmosphere
    !
    !tslu = 0.
    !tsld = 0.
    !
    !	  <fupbs>, <fdownbs>, and <fnetbs> are total upwelling, downwelling, and net
    !	  solar fluxes at grid boundaries
    !
    !DO  j = 1, nlayer
    !  fupbs(j) = 0.
    !  fdownbs(j) = 0.
    !  fnetbs(j) = 0.
    !END DO
    !
    !	  <fsLu> and <fsLd> are upwelling, downwelling, and net
    !	  solar fluxes at top-of-atmosphere (spectrally-resolved)
    !
    !	  <alb_toa> is albedo at top-of-atmosphere (spectrally-resolved)
    !
    !DO  i = 1, nsol
    !  fslu(i) = 0.
    !  fsld(i) = 0.
    !  alb_toa(i) = 0.
    !END DO
    !
    !	   <alb_tomi> and <alb_toai> are total solar albedos at top-of-model
    !	   and top-of-atmosphere
    !
    !alb_tomi = 0.
    !alb_toai = 0.
    !
    !	  CALCULATE SOLAR ABSORBED BY GROUND, SOLNET, AND UPWARD AND DOWNWARD
    !	  LONGWAVE FLUXES AT SURFACE, XIRUP AND XIRDOWN (WATTS/M**2)
    !
    solnet  = 0.0
   
    DO  l   =  1,nsolp
      DO ij=1,iend
  	IF (isl_aerad(ij) /= 0) THEN
  	  solnet(ij) = solnet(ij) - fnet(ij,l,nlayer)
    !fp = ck1(l,1)*el2(l,1) - ck2(l,1)*em2(l,1) + cp(l,1)
    !fslu( nprob(l) ) = fslu( nprob(l) ) + fp
    !DO  j = 1, nlayer
    !	  fp = ck1(l,j)*el1(l,j) + ck2(l,j)*em1(l,j) + cpb(l,j)
    !	  fupbs(j) = fupbs(j) + fp
    !	  fnetbs(j) = fnetbs(j) + fnet(l,j)
    !	  IF (l == nsolp) fdownbs(j) = fupbs(j) - fnetbs(j)
    !	 END DO
  	END IF
      END DO
    END DO
      !DO  i = 1, nsol
  	!fsld(i) = u0*solfx(i)
  	!alb_toa(i) = fslu(i)/fsld(i)
  	!tslu = tslu + fslu(i)
  	!tsld = tsld + fsld(i)
      !END DO
    !
      !alb_tomi = fupbs(1)/fdownbs(1)
      !alb_toai = tslu/tsld
    !
    !	   Load albedos into interface common block
    !
      !alb_toai_aerad = alb_toai
      !alb_tomi_aerad = alb_tomi
      !DO i = 1, nsol
      !  alb_toa_aerad(i) = alb_toa(i)
      !END DO
    !
    !	   Load fluxes into interface common block
    !
      !DO j = 1, nlayer
      !  fsl_up_aerad(j) = fupbs(j)
      !  fsl_dn_aerad(j) = fdownbs(j)
      !END DO
   
   ! END DO
    !
    !	  <tiru> is total upwelling infrared flux at top-of-atmosphere;
    !	  <fupbi>, <fdownbi>, and <fnetbi> are total upwelling, downwelling, and net
    !	  infrared fluxes at grid boundaries
    !
    !tiru = 0.
    !DO  j = 1, nlayer
    !  fupbi(j)   =  0.0
    !  fdownbi(j)   =  0.0
    !  fnetbi(j)  =  0.0
    !END DO
    !
    !	  <firu> is upwelling infrared flux at top-of-atmosphere (spectrally-resolved)
    !
    !DO  i = 1, nir
    !  firu(i) = 0.
    !END DO
   
    xirdown = 0.0
    !xirup   = 0.0
   
    DO  l  =  nsolp+1,ntotal
      DO ij=1,iend
  	IF (ir_aerad /= 0) THEN
  	    xirdown(ij) = xirdown(ij) + direc(ij,l,nlayer)

    !	 xirup   = xirup  + directu(l,nlayer)
    !	 firu( nprob(l)-nsol ) = firu( nprob(l)-nsol ) + directu(l,1)
    !	 DO  j = 1, nlayer
    !	  fupbi(j) = fupbi(j) + directu(l,j)
    !	  fdownbi(j) = fdownbi(j) + direc  (l,j)
    !	  fnetbi(j) = fnetbi(j) + directu(l,j) - direc(l,j)
    !	 END DO
  	END IF
      END DO
    END DO
 
   
    !  DO  i = 1, nir
    !	 tiru = tiru + firu(i)
    !  END DO
    !
    !	   Load fluxes into interface common block
    !
    !  DO j = 1, nlayer
    !	 fir_up_aerad(j) = fupbi(j)
    !	 fir_dn_aerad(j) = fdownbi(j)
    !  END DO
  
  END SUBROUTINE radtran
  
  SUBROUTINE oppr(ia,iz,ja,jz,m1,aot11)
    !
    !	  **************************************************************
    !	  *  Purpose		 :  CaLculates optical properties      *
    !	  *			    such as single scattering albedo,  *
    !	  *			    asymmetry parameter, etc.	       *
    !	  *			    This routine is case dependent and *
    !	  *			    wiLL have to be repLaced by the    *
    !	  *			    user.			       *
    !	  *  Subroutines Called  :  None			       *
    !	  *  Input		 :  PAH2O, RDH2O, CO2, O3, ETC         *
    !	  *  Output		 :  TAUL, W0, G0, OPD,         *
    !	  * ************************************************************
    !
    !INCLUDE 'globrad.h'
    !
    !	  W0(NWAVE,NLAYER) : SINGLE SCATTERING ALBEDO *** delta scaled ***
    !	  G0(NWAVE,NLAYER) : ASYMMETRY PARAMETER *** delta scaled ***
    !	  OPD(NWAVE,NLAYER): cumulative OPTICAL DEPTH *** delta scaled ***
    !	  SFL(NWAVE)	   : SOLAR FLUX
    !	 uW0(NWAVE,NLAYER)  : unscaled SINGLE SCATTERING ALBEDO
    !	 uG0(NWAVE,NLAYER)  : unscaled ASYMMETRY PARAMETER
    !	 uTAUL(NWAVE,NLAYER): unscaled OPTICAL DEPTH of layer
    !
    !	  ASSUME THAT P IS SAME ON ALL SIGMA LEVELS. IF PSURFACE
    !	  VARIES A LOT, THEN WE MAY NEED TO CALCULATE TAUO3,
    !	  TAUCO2, TAUO2 FOR EACH POINT.
    !
    !	  NOTE : THE TOP LAYER IS A DUMMY. IT CONTAINS A DEFINED GAS
    !		 AMOUNT. DIFFERENT MODELS WILL REQUIRE DIFFERENT
    !		 TREATMENT OF THIS.
    !	  CALCULATE TOTAL OPTICAL DEPTH INCLUDING GASES. THEN
    !	  GIVEN THE AEROSOL OPTICAL DEPTHS AND CLOUD OPTICAL DEPTHS,
    !	  CALCULATE FINAL OPTICAL PROPERTIES. WE USE A DELTA
    !	  TWO STREAM APPROACH TO FIND W0, SINGLE SCATTERING ALBEDO,
    !	  G0, ASYMMMETRY PARAMETER, TAUL, LAYER OPTICAL DEPTH,
    !	  OPD, CUMULATIVE OPTICAL DEPTH TO BASE OF LAYER.
    !
    USE mem_globrad, ONLY: nlayer,nwave,ngroup,nrad,xsecta, &
  			   rdqext,qscat,qbrqs,ntotal,nprob, &
  			   epsilon,g,ptop,p,q,nsolp,contnm, &
  			   uw0,ug0,ir,ngauss,gangle,ta,tb,  &
			   wa,wb,ga,gb,tia,tib,wia,wib,gia, &
			   gib,alpha,gama,caseE,caseW,  &          !kml2
          		   caseG,wave,imie               !kml2
    
    USE mem_aerad, ONLY: iprocopio
    use rconstants, only: t00
    IMPLICIT NONE
    
    INTEGER,INTENT(IN) :: ia,iz,ja,jz,m1
    REAL,DIMENSION((iz-ia+1)*(jz-ja+1),ntotal,nlayer) :: taua,taus,g01,wol
    REAL,DIMENSION((iz-ia+1)*(jz-ja+1),ntotal,nlayer) :: gol
    INTEGER :: i,i1,j1,kk
    INTEGER :: ig
    INTEGER :: iradgas
    INTEGER :: j
    INTEGER :: l
    REAL    :: cco((iz-ia+1)*(jz-ja+1))
    REAL    :: den 
    REAL    :: denom
    REAL    :: fo
    REAL    :: pcorr
    REAL    :: qcorr
    REAL    :: ttas
    REAL    :: tauh2o((iz-ia+1)*(jz-ja+1),ntotal,nlayer)
    REAL    :: utaul((iz-ia+1)*(jz-ja+1),ntotal,nlayer)
  !  REAL    :: uw0((iz-ia+1)*(jz-ja+1),ntotal,nlayer)
    REAL    :: wot((iz-ia+1)*(jz-ja+1),ntotal)
    REAL    :: got((iz-ia+1)*(jz-ja+1),ntotal)
    INTEGER,DIMENSION((iz-ia+1)*(jz-ja+1)):: idaot
    REAL :: aot11((iz-ia+1)*(jz-ja+1))         !kml2
    INTEGER :: ij,iend,jjj,in
    
!kmlnew
    REAL,DIMENSION((iz-ia+1)*(jz-ja+1),ntotal,nlayer) :: taucld,wcld,gcld 
    REAL,DIMENSION((iz-ia+1)*(jz-ja+1),ntotal,nlayer) :: wolc,woice,worain,gl,gice,grain
    REAL,DIMENSION((iz-ia+1)*(jz-ja+1),ntotal,nlayer) :: DENC
    REAL,DIMENSION((iz-ia+1)*(jz-ja+1),ntotal,nlayer) :: taucldlw,taucldice,taurain
    REAL,DIMENSION((iz-ia+1)*(jz-ja+1),nlayer) :: CORR,REFFI
    
    REAL,DIMENSION((iz-ia+1)*(jz-ja+1),nwave) :: rdqextnew,wonew,gonew
    real X_teste
  
     
    CORR=0.0
    DENC=0.0
    REFFI=0.0
    taucldlw=0.0
    taucldice=0.0
    taurain=0.0
    wolc=0.0
    woice=0.0
    worain=0.0
    gl=0.0
    gice=0.0
    grain=0.0
    taucld=0.0
    wcld=0.0
    gcld=0.0
    rdqextnew=0.0
    wonew=0.0
    gonew=0.0

!kmlnew    
       
    iend=(iz-ia+1)*(jz-ja+1)



      IF (iprocopio == 1 .and. imie == 1) THEN
!            aot11=0.1          !TMP KML2
 
   	    DO ij=1,iend
		idaot(ij) = MAX(MIN(INT(10*((ANINT(10.*aot11(ij))/10.)+0.1)/2.),9),1)

  	        DO  l = 1,nwave
	         rdqextnew(ij,l) = caseE(idaot(ij),l) 
		 wonew(ij,l)	 = caseW(idaot(ij),l) 
                 gonew(ij,l)	 = caseG(idaot(ij),l)
                 !if(l.eq.11) print*,'ext,wo,go=',rdqextnew(ij,l),wonew(ij,l),gonew(ij,l),ij
                END DO
            END DO

!TMP       
           taua=0.0
           DO j=1,nlayer
             DO ig = 1,ngroup
         	DO  i = 1,nrad
         	  DO  l = 1,ntotal
         	    DO ij=1,iend

         		taua(ij,l,j)=taua(ij,l,j)+rdqextnew(ij,nprob(l))*xsecta(i,ig)* &
         				   caer(ij,j,i,ig)
			tauaer(ij,l,j) = MAX(taua(ij,nprob(l),j),REAL(epsilon))		   
         		wol(ij,l,j)    = wonew(ij,nprob(l))
         		gol(ij,l,j)    = gonew(ij,nprob(l))
!               		if(tauaer(ij,l,j).gt.REAL(epsilon)) then
!			  print*,'ij,j,nprob=',ij,l,nprob(l)
!			  print*,'wo,go,ext,taua=',wol(ij,l,j),gol(ij,l,j),rdqextnew(ij,nprob(l))
!			  print*,'caer, taua=',taua(ij,l,j),caer(ij,j,i,ig)
!			endif
         	    END DO
         	  END DO
         	END DO
             END DO

           END DO
      ELSE

     	  taua=0.0
     	  taus=0.0
     	  g01=0.0
     	  DO j=1,nlayer
     	    DO ig = 1,ngroup
     	       DO  i = 1,nrad
     		 DO  l = 1,nwave
     		   DO ij=1,iend
     		     taua(ij,l,j)=taua(ij,l,j)+rdqext(i,ig,l)*xsecta(i,ig)* &
     				  caer(ij,j,i,ig)
     		     taus(ij,l,j)=taus(ij,l,j)+qscat(i,ig,l)*xsecta(i,ig)* &
     				  caer(ij,j,i,ig)
     		     g01(ij,l,j) =g01(ij,l,j) +qbrqs(i,ig,l)*xsecta(i,ig)* &
     				  caer(ij,j,i,ig)
     		   END DO
     		 END DO
     	       END DO
     	    END DO


     	    DO l= 1,ntotal
     	       DO ij=1,iend
     		 tauaer(ij,l,j) = MAX(taua(ij,nprob(l),j),REAL(epsilon))
!     		 tauaer(ij,l,j) =     taua(ij,nprob(l),j)
    		 wol(ij,l,j)	 = taus(ij,nprob(l),j)/tauaer(ij,l,j)
     		 ttas=1.0
     		 IF( wol(ij,l,j) /= 0. ) ttas = taus(ij,nprob(l),j)
     		 gol(ij,l,j)	= g01(ij,nprob(l),j)/ttas
     	       END DO
     	    END DO
     	  END DO
	  imie = 1
      END IF

!kmlnew 

!     imie = 1 

    DO j=1,nlayer
       DO ij=1,iend
       
       IF (xland_aerad(ij).ge..009) THEN
         REFFI(ij,j) =  7.0 * 1.e+3 * LWL_aerad(ij,j) + 5.5 
       ELSE 
         REFFI(ij,j) =  9.5 * 1.e+3 * LWL_aerad(ij,j) + 4.0 
       END IF
       
       CORR(ij,j) = 1.047 - 0.913e-4 * (tt(ij,j)-t00) + 0.203e-3 * &
                   (tt(ij,j)-t00) **2 - 0.106e-4 * (tt(ij,j)-t00) **3
      
       CORR(ij,j) = MAX(CORR(ij,j),REAL((epsilon)))
              
	 
       END DO  
    END DO

   DO j=1,nlayer
      DO l= 1,ntotal
       DO ij=1,iend
        
	IF( j .eq. 1) taurain(ij,l,j) = 0.00018 * RAIN_aerad(ij) * 2000.0
	
	taucldlw(ij,l,j)= 1.e+3 * LWP_aerad(ij,j) *(ta(l)/REFFI(ij,j)+tb(l)/REFFI(ij,j)**2)

	worain(ij,l,j) = 1. - 0.45
	wolc(ij,l,j) = (1. - wa(l)) + wb(l) * REFFI(ij,j)
	gl(ij,l,j)   = ga(l) + gb(l) * REFFI(ij,j)
	grain(ij,l,j) = 0.95

	DENC(ij,l,j)      = 1. / (tia(l) + tib(l) * 1.0e+3 * IWL_aerad(ij,j))
	taucldice(ij,l,j) = CORR(ij,j) * 1.0e+3 * IWP_aerad(ij,j) * DENC(ij,l,j)
	
!srf-evitando bug - 0**0	
        if(wib(l) < 1.e-5) then
	   X_teste = 0.
	else
	   X_teste = (1.0e+3 * IWL_aerad(ij,j))**wib(l)
	endif

!	woice(ij,l,j) =( 1.0 -  wia(l) * (1.0e+3 * IWL_aerad(ij,j))**wib(l)) * & 
	woice(ij,l,j) =( 1.0 -  wia(l) *  X_teste                          ) * & 
	               ( 1.0 - gama(l) * (CORR(ij,j) - 1)/ CORR(ij,j) )
	
!srf-evitando bug - 0**0	
        if(gib(l) < 1.e-5) then
!	   X_teste = 0.
	   X_teste = 1.
	   
	else
	   X_teste = (1.0e+3 * IWL_aerad(ij,j))**gib(l)
	endif
!	gice(ij,l,j) = ( gia(l) * (1.0e+3 * IWL_aerad(ij,j))**gib(l) ) * &
	gice(ij,l,j) = ( gia(l) * X_teste                            ) * &
	               ( 1.0 + alpha(l) * (CORR(ij,j) - 1)/ CORR(ij,j) )
		       
	
	IF (l>=91 .AND. l<=113) THEN
	
	taucldlw(ij,l,j)= 1.0e+3 * LWP_aerad(ij,j) * ta(l) * exp(tb(l)* REFFI(ij,j))
	
	END IF
	
	IF (l>=114 .AND. l<=154) THEN
	
	taucldlw(ij,l,j)=  1.0e+3 * LWP_aerad(ij,j) * ( ta(l) + tb(l)* REFFI(ij,j))
	
	gl(ij,l,j) = 1. - ga(l) * exp( gb(l) * REFFI(ij,j))
	
	END IF
	
	taucld(ij,l,j)= taucldlw(ij,l,j) + taucldice(ij,l,j) + taurain(ij,l,j)

!        if(LWL_aerad(ij,j).gt.0.) print*, 'Liquid', 1.0e+3 * LWL_aerad(ij,j),taucldlw(ij,l,j)
!	if(IWL_aerad(ij,j).gt.0.) print*, 'Ice', 1.0e+3 * IWL_aerad(ij,j),taucldice(ij,l,j)

        IF ( taucld(ij,l,j).gt.epsilon) THEN
	
	  wcld(ij,l,j) =  (wolc(ij,l,j) *  taucldlw(ij,l,j)  + &
	                  woice(ij,l,j) * taucldice(ij,l,j) + &
			  worain(ij,l,j) * taurain(ij,l,j)) &
	              / taucld(ij,l,j)
	  gcld(ij,l,j) = ( wolc(ij,l,j) *  taucldlw(ij,l,j)*   gl(ij,l,j) + &
	                woice(ij,l,j) * taucldice(ij,l,j)* gice(ij,l,j) + &
			worain(ij,l,j) * taurain(ij,l,j) * grain(ij,l,j)) &
		      / (wcld(ij,l,j) * taucld(ij,l,j))	     
        ELSE 
	  wcld(ij,l,j) = 1.0
	  gcld(ij,l,j) = 0.0
        ENDIF
       END DO                           
      END DO
    END DO
!kmlnew 
    
    iradgas = 1 !iradgas = 0: no gas in radiative xfer
    DO  j = 1,nlayer
      kk = MAX( 1, j-1 )
      !
      !   Bergstrom water vapor continuum fix:
      !
      !    <qcorr> is layer average water vapor mixing ratio
      !    <pcorr> is layer average pressure
      !
      !   For layer 0, calculate mixing ratio [g/g] from vapor column [g/cm^2]
      !   and average pressure [dyne/cm^2]
      !
      IF( j == 1 )THEN
        DO ij=1,iend
  	  qcorr = rdh2o(ij,1) * g / p_top(ij)
  	  pcorr = p_aerad(ij,1) / 2.
  	  cco(ij) = EXP(1800./t_aerad(ij,kk))*(qcorr*pcorr/2.87 + pcorr/4610.)
        END DO
      ELSE
        DO ij=1,iend
  	  qcorr = qv_aerad(ij,kk)
  	  pcorr = p_aerad(ij,kk)
  	  cco(ij) = EXP(1800./t_aerad(ij,kk))*(qcorr*pcorr/2.87 + pcorr/4610.)
        END DO
      END IF
  
      DO  l   = 1,ntotal
  	DO ij=1,iend
  	  IF (l>=lls(ij) .AND. l<=lla(ij)) THEN
  	    tauh2o(ij,l,j) = rdh2o(ij,j)*pah2o(ij,l,j)
  	    !	  Bergstrom water vapor continuum fix (next two statements)
  	    IF( l > nsolp+30 .AND. l <= nsolp+36 ) THEN
  	      !kml	   if( L .GT. NSOLP+36 .AND. L .LE. NSOLP+42 ) then
  		  tauh2o(ij,l,j) = rdh2o(ij,j)*pah2o(ij,l,j)*cco(ij)
  	    ELSE
  		  tauh2o(ij,l,j) = tauh2o(ij,l,j)
  	    END IF
  	    IF (l > nsolp+36) tauh2o(ij,l,j) = tauh2o(ij,l,j) + &
  			    cco(ij)*rdh2o(ij,j)*contnm(l-nsolp)
     
  	    taul(ij,l,j)   = tauh2o(ij,l,j)+taugas(ij,l,j)+ &
  			      paray(ij,l,j)+tauaer(ij,l,j)+taucld(ij,l,j)
  
  	    IF (iradgas == 0) taul(ij,l,j) = tauaer(ij,l,j)
  	    IF( taul(ij,l,j) < epsilon ) taul(ij,l,j) = epsilon
  	  END IF
  	END DO
      END DO
  
      DO  l   = 1,ntotal
  	DO ij=1,iend

  	  IF (l>=lls(ij) .AND. l<=lla(ij)) THEN
  	    utaul(ij,l,j)  = taul(ij,l,j)
  	    wot(ij,l)	   = (paray(ij,l,j)+tauaer(ij,l,j)*wol(ij,l,j)+  &
  			 taucld(ij,l,j)*wcld(ij,l,j))/taul(ij,l,j)
  	    IF (iradgas == 0) wot(ij,l) = wol(ij,l,j)
  	    wot(ij,l)	      = MIN(1.-REAL(epsilon),wot(ij,l))
  !	    uw0(ij,l,j)    = wot(ij)
  	    denom     = (paray(ij,l,j)+taucld(ij,l,j)*wcld(ij,l,j)+ &
  			 tauaer(ij,l,j)*wol(ij,l,j))
  	    !IF( denom <= epsilon ) denom = epsilon
  	    IF( denom > epsilon ) THEN
  	      got(ij,l) = ( wcld(ij,l,j)*gcld(ij,l,j)*taucld(ij,l,j) + &
  		      gol(ij,l,j)* wol(ij,l,j)*tauaer(ij,l,j) ) / denom
	      got(ij,l)=max(REAL(epsilon), got(ij,l))
  	    ELSE
  	      got(ij,l) = 0.
  	    END IF
  	    IF (iradgas == 0) got(ij,l) = gol(ij,l,j)
  	  END IF
  	END DO
      END DO

      DO  l   = 1,ntotal
  	DO ij=1,iend
  	  IF (l>=lls(ij) .AND. l<=lla(ij)) THEN
  	    !ug0(l,j)	 = got(ij)
  	    fo        = got(ij,l)**2
  	    den       = 1.-wot(ij,l)*fo
  	    taul(ij,l,j)   = taul(ij,l,j) * den
  	    w0(ij,l,j)      = (1.-fo)*wot(ij,l)/den
  	    g0(ij,l,j)      = got(ij,l)/(1.+got(ij,l))
  	    opd(ij,l,j)    = 0.0
  	    opd(ij,l,j)    = opd(ij,l,kk)+taul(ij,l,j)
  	    uopd(ij,l,j)   = 0.0
  	    uopd(ij,l,j)   = uopd(ij,l,kk)+utaul(ij,l,j)
  	  END IF
  	END DO
      END DO
    END DO
  
    IF(ir_aerad == 1) THEN
      DO   j =   1,nlayer
  	DO  i =   1,ngauss
  	  DO  l =   1,ntotal
  	    DO ij=1,iend
  	      IF(l>=lls(ij) .AND. l<=lla(ij)) &
  		   y3(ij,l,i,j) = EXP(-taul(ij,l,j)/gangle(i))
  	    END DO
  	  END DO
  	END DO
      END DO
    END IF
  
  END SUBROUTINE oppr 
  
  SUBROUTINE oppr1(ia,iz,ja,jz,m1)
    !
    !	  **********************************************************
    !	  *  Purpose		 :  Calculate Planck Function and  *
    !	  *			    and its derivative at ground   *
    !	  *			    and at all altitudes.	   *
    !	  *  Subroutines Called  :  None			   *
    !	  *  Input		 :  TGRND, NLOW, WEIGHT 	   *
    !	  *  Output		 :  PTEMP, PTEMPG, SLOPE	   *
    !	  * ********************************************************
    !
    USE mem_globrad, ONLY: ntotal,tgrnd,nlow,nirp,plank,ltemp,nsolp, &
  			   weight,iblackbody_above,t_above, &
  			   nlayer,ncount
    
    IMPLICIT NONE
    
    INTEGER,INTENT(IN) :: ia,iz,ja,jz,m1
    INTEGER :: it1((iz-ia+1)*(jz-ja+1),nlayer)
    INTEGER :: itg((iz-ia+1)*(jz-ja+1))
    INTEGER :: itp((iz-ia+1)*(jz-ja+1))
    INTEGER :: j
    INTEGER :: kindex
    INTEGER :: l,i1,j1
    REAL :: pltemp1((iz-ia+1)*(jz-ja+1),ntotal)
    REAL :: ptemp2((iz-ia+1)*(jz-ja+1),ntotal,nlayer)
    INTEGER :: ij,iend,i
    
    iend=(iz-ia+1)*(jz-ja+1)
      !
    !	  **************************************
    !	  * CALCULATE PTEMP AND SLOPE	       *
    !	  **************************************
    !
    !	  CALCULATE THE WAVELENGTH DEPENDENT PLANK FUNCTION AT THE GROUND.
  
    DO ij=1,iend
      itg(ij)= ANINT(100.*t_surf(ij)) - nlow
!srf
!      if( itg(ij) < 0. .OR. itg(ij) > ncount) print*,'1-ITG=',itg(ij)
!srf      
    END DO
    DO i=1,nirp
      DO ij=1,iend
        pltemp1(ij,i)=plank(ltemp(i),itg(ij))
      END DO
    END DO
    DO  l =   nsolp+1,ntotal
      DO ij=1,iend
  	 ptempg(ij,l)=   pltemp1(ij,l-nsolp)*weight(l)
      END DO
    END DO
    !
    IF( iblackbody_above /= 0 )THEN
    !	    CALCULATE THE WAVELENGTH DEPENDENT PLANK FUNCTION AT THE TOP
    !	    OF THE MODEL.
      DO ij=1,iend
  	  itp(ij) = ANINT(100.*tabove_aerad(ij)) - nlow
!srf
!      if( itp(ij) < 0. .OR. itp(ij) > ncount) print*,'2-ITP=',itp(ij)
!srf      
      END DO
      DO i=1,nirp
        DO ij=1,iend  	 
	   pltemp1(ij,i)=plank(ltemp(i),itp(ij))
	   !CALL gather(nirp,pltemp1(ij,:),plank(1,itp(ij)),ltemp)
        END DO
      END DO
      DO  l =	nsolp+1,ntotal
  	DO ij=1,iend
  	    ptempt(ij,l)	=   pltemp1(ij,l-nsolp)*weight(l)
  	END DO
      END DO
  
    END IF
    !
    DO  j	     =   1,nlayer
      DO ij=1,iend
  	it1(ij,j) = ANINT(100.*tt(ij,j)) - nlow
!srf
!      if( it1(ij,j) < 0. .OR. it1(ij,j) > ncount) print*,'3-IT1=',it1(ij,j)
!srf      
      END DO
    END DO
    DO  j	     =   1,nlayer
      DO i=1,nirp
        DO ij=1,iend
           ptemp2(ij,i,j)=plank(ltemp(i),it1(ij,j))
        END DO
      END DO
    END DO
    
   ! kindex makes the top layer isothermal. using kindex, find
   ! plank function at bottom of each layer.
   ! note: if you force slope=0, then you have isothermal
   ! layers with tt(j) corresponding to average temperature
   ! of layer and tt(nlayer) should be set to tgrnd.
    DO  j	     =   1,nlayer
      kindex	      = MAX( 1, j-1 )
      DO  l	   = nsolp+1,ntotal
  	DO ij=1,iend
  	  ptemp(ij,l,j)   = ptemp2(ij,l-nsolp,j)*weight(l)
  	  slope(ij,l,j)   = (ptemp(ij,l,j)-ptemp(ij,l,kindex))/ &
  				 taul(ij,l,j)
  	  IF( taul(ij,l,j) <= 1.0E-6 ) slope(ij,l,j) = 0.
  	END DO
      END DO
    END DO
  
    !
  END SUBROUTINE oppr1
  
  SUBROUTINE twostr(m1,ia,iz,ja,jz)
    !
    !	 ******************************************************************
    !	 *  Purpose		:  Defines matrix properties and sets up  *
    !	 *			   matrix coefficients that do not depend *
    !	 *			   on zenith angle or temperature.	  *
    !	 *  Subroutines Called  :  None 				  *
    !	 *  Input		:  W0, G0				  *
    !	 *  Output		:  B1, B2, GAMI, ACON, EL1, AF, ETC	  *
    !	 * ****************************************************************
    !
  
    USE mem_globrad, ONLY: nsolp,sq3,tpi,nlayer,jn,jdble,irs,ntotal
  
    IMPLICIT NONE
  
    INTEGER,INTENT(IN) :: m1,ia,iz,ja,jz
    INTEGER	   :: j
    INTEGER	   :: jd
    INTEGER	   :: l
    REAL,PARAMETER :: two = 2.d0
    INTEGER :: ij,iend
    
    iend=(iz-ia+1)*(jz-ja+1)
      
    DO  l    =  1,ntotal !lls(i1,j1),lla(i1,j1)
      DO ij=1,iend  
  	IF(isl_aerad(ij) /= 0 .OR. irs .NE. 0 ) THEN
  	  IF( l>=lls(ij) .AND. l<=lla(ij)) THEN
  	    IF(l <= nsolp ) THEN
  	      u1i(ij,l) = sq3
  	    ELSE
  	      u1i(ij,l) = two
  	    END IF
  	    !u1s(l)  =  tpi/u1i(l)
  	  END IF
  	END IF
      END DO
    END DO
    !
    !	   here we define layer properties following general scheme
    !	   of meador and weavor. then we set up layer properties
    !	   needed for matrix.
    !
    DO  j =  1,nlayer
      DO  l=  1,ntotal
  	DO  ij=  1,iend
  	  IF(isl_aerad(ij) /= 0 .OR. irs .NE. 0 ) THEN
  	    IF( l>=lls(ij) .AND. l<=lla(ij)) THEN
  	      !these are for two stream and hemispheric means
  	      b1(ij,l,j)   =  0.5*u1i(ij,l)*(2.-w0(ij,l,j)*(1. + g0(ij,l,j)))
  	      b2(ij,l,j)   =  0.5*u1i(ij,l)*w0(ij,l,j)*(1. - g0(ij,l,j))
  	      ak(ij,l,j)   =  SQRT(ABS(b1(ij,l,j)**2 - b2(ij,l,j)**2))
  	      gami(ij,l,j)  =  b2(ij,l,j)/(b1(ij,l,j) + ak(ij,l,j))
  	      ee1(ij,l,j)   =  EXP(-ak(ij,l,j)*taul(ij,l,j))
  	      el1(ij,l,j)   =  1.0 + gami(ij,l,j) *ee1(ij,l,j)
  	      em1(ij,l,j)   =  1.0 - gami(ij,l,j) * ee1(ij,l,j)
  	      el2(ij,l,j)   =  gami(ij,l,j) + ee1(ij,l,j)
  	      em2(ij,l,j)   =  gami(ij,l,j) - ee1(ij,l,j)
  	    END IF
  	  END IF
  	END DO
      END DO
    END DO
    !
    !	  we seek to solve ax(l-1)+bx(l)+ex(l+1) = d.
    !	  l=2n for even l, l=n+1 for odd l. the mean intensity (tmi/4pi)
    !	  and the net flux (fnet) are related to x's as noted in add.
    !	  first we set up the coefficients that are independent of solar
    !	  angle or temparature: a(i),b(i),e(i). d(i) is defined in add.
    !
    j=  0
    DO  jd=  2,jn,2
      j=  j + 1
      DO  l=  1,ntotal
  	DO  ij=  1,iend
  	  IF(isl_aerad(ij) /= 0 .OR. irs .NE. 0 ) THEN
  	    IF( l>=lls(ij) .AND. l<=lla(ij)) THEN
  	      !here are the even matrix elements
  	      af(ij,l,jd)   =  em1(ij,l,j+1)*el1(ij,l,j)- &
  				  em2(ij,l,j+1)*el2(ij,l,j)
  	      bf(ij,l,jd)   =  em1(ij,l,j+1)* em1(ij,l,j)- &
  				  em2(ij,l,j+1)*em2(ij,l,j)
  	      ef(ij,l,jd)   =  el1(ij,l,j+1)*em2(ij,l,j+1) - &
  				  el2(ij,l,j+1)*em1(ij,l,j+1)
  	      !here are the odd matrix elements except for the top.
  	      af(ij,l,jd+1) =  em1(ij,l,j)*el2(ij,l,j)- &
  				  el1(ij,l,j)*em2(ij,l,j)
  	      bf(ij,l,jd+1) =  el1(ij,l,j+1)*el1(ij,l,j) - &
  				  el2(ij,l,j+1)*el2(ij,l,j)
  	      ef(ij,l,jd+1) =  el2(ij,l,j)*em2(ij,l,j+1)- &
  				el1(ij,l,j)*em1(ij,l,j+1)
  	    END IF
  	  END IF
  	END DO
      END DO
    END DO
    !
    !	  HERE ARE THE TOP AND BOTTOM BOUNDARY CONDITIONS AS WELL AS THE
    !	  BEGINNING OF THE TRIDIAGONAL SOLUTION DEFINITIONS. I ASSUME
    !	  NO DIFFUSE RADIATION IS INCIDENT AT UPPER BOUNDARY.
    !
    DO  l=  1,ntotal
      DO  ij=  1,iend
  	IF(isl_aerad(ij) /= 0 .OR. irs .NE. 0 ) THEN
  	  IF( l>=lls(ij) .AND. l<=lla(ij)) THEN
  	    af(ij,l,1)    = 0.0
  	    bf(ij,l,1)    = el1(ij,l,1)
  	    ef(ij,l,1)    = -em1(ij,l,1)
  	    af(ij,l,jdble) = el1(ij,l,nlayer)-rsfx(ij,l)*el2(ij,l,nlayer)
  	    bf(ij,l,jdble) = em1(ij,l,nlayer)-rsfx(ij,l)*em2(ij,l,nlayer)
  	    ef(ij,l,jdble) = 0.0
  	  END IF
  	END IF
      END DO
    END DO
  
  END SUBROUTINE twostr
  
  SUBROUTINE add(m1,ia,iz,ja,jz,cosz,m2,m3)
  
    USE mem_globrad, ONLY: isl,u0,nlayer,nsolp,sq3,sol,epsilon, &
  			   irs,ntotal,u1s,emis,pi,jn,tpi, &
  			   jdble,ndbl
    
    IMPLICIT NONE 
   
    !	  THIS SUBROUTINE FORMS THE MATRIX FOR THE MULTIPLE LAYERS AND
    !	  USES A TRIDIAGONAL ROUTINE TO FIND RADIATION IN THE ENTIRE
    !	  ATMOSPHERE.
   
    !	  ******************************
    !	  *   CALCULATIONS FOR SOLAR   *
    !	  ******************************
    INTEGER,INTENT(IN) :: ia,iz,ja,jz,m1,m2,m3
    REAL,INTENT(IN),DIMENSION((iz-ia+1)*(jz-ja+1)) :: cosz  
    INTEGER :: j,kk
    INTEGER :: jd
    INTEGER :: kindex
    INTEGER :: l
    REAL    :: b4
    REAL    :: c1
    REAL    :: c2
    REAL    :: cm1
    REAL    :: cp1
    REAL    :: du0
    REAL    :: x
    REAL    :: x2
    REAL    :: x3
    REAL    :: x4
    REAL,DIMENSION((iz-ia+1)*(jz-ja+1),nsolp,nlayer) :: direct,el3,ee3,cm
    REAL,DIMENSION((iz-ia+1)*(jz-ja+1),nsolp) :: sfcs
    REAL,DIMENSION((iz-ia+1)*(jz-ja+1),ntotal,ndbl) :: df,as,ds,xk
    INTEGER :: ij,iend,i1,j1
    
    iend=(iz-ia+1)*(jz-ja+1)
  
    DO  j	 =  1,nlayer
      kk = MAX( 1, j-1 )
      DO  l    =  1,nsolp
  	DO ij=1,iend
  	  du0=1./cosz(ij)
  	  IF(isl_aerad(ij) /= 0)  THEN
  	    b3(ij,l,j)     =  0.5*(1.-sq3*g0(ij,l,j)*cosz(ij))
  	    b4         =  1. - b3(ij,l,j)
  	    x2         =  taul(ij,l,j)*du0
  	    ee3(ij,l,j)   =  EXP(-x2)
  	    x3         =  opd(ij,l,j)*du0
  	    el3(ij,l,j)   =  EXP(-x3)*sol(l)
  	    direct(ij,l,j) =  cosz(ij)*el3(ij,l,j)
  	    c1         =  b1(ij,l,j) - du0
  	    IF( ABS(c1) < epsilon ) c1 = SIGN(REAL(epsilon),c1)
  	    c2         =  ak(ij,l,j)*ak(ij,l,j) - du0*du0
  	    IF( ABS(c2) <= epsilon ) c2 = epsilon
  	    cp1        =  w0(ij,l,j)*(b3(ij,l,j)*c1+b4*b2(ij,l,j))/c2
  	    cpb(ij,l,j)    =  cp1 * el3(ij,l,j)
  	    IF( j /= 1 ) THEN
  	      x4 = el3(ij,l,kk)
  	    ELSE
  	      x4 = sol(l)
  	    END IF
  	    cp(ij,l,j)     =  cp1 * x4
  	    cm1        =  ( cp1*b2(ij,l,j) + w0(ij,l,j)*b4 )/c1
  	    cmb(ij,l,j)    =  cm1 * el3(ij,l,j)
  	    cm(ij,l,j)    =  cm1 * x4
  	  END IF
  	END DO
      END DO
    END DO
    !	     CALCULATE SFCS, THE SOURCE AT THE BOTTOM.
    DO  l=  1,nsolp
      DO ij=1,iend
  	IF(isl_aerad(ij) /= 0)  THEN
  	  sfcs(ij,l)=  direct(ij,l,nlayer) * rsfx(ij,l)
  	END IF
      END DO
    END DO
   
    !	  ******************************
    !	  * CALCULATIONS FOR INFRARED. *
    !	  ******************************
    DO  j= 1,nlayer
      DO  l = nsolp+1,ntotal
  	DO ij=1,iend
  	  IF(irs /= 0)  THEN
  	      kindex = MAX(1,j-1)
  	      b3(ij,l,j)     = 1.0/(b1(ij,l,j)+b2(ij,l,j))
  	      cp(ij,l,j)     = (ptemp(ij,l,kindex)+slope(ij,l,j)* &
  				   b3(ij,l,j))*(tpi/u1i(ij,l))
  	      cpb(ij,l,j)    = cp(ij,l,j) + slope(ij,l,j)* &
  				  taul(ij,l,j)*(tpi/u1i(ij,l))
  	      cm(ij,l,j)     = (ptemp(ij,l,kindex)-slope(ij,l,j)* &
  				   b3(ij,l,j))*(tpi/u1i(ij,l))
  	      cmb(ij,l,j)    = cm(ij,l,j) + slope(ij,l,j)* &
  				  taul(ij,l,j)*(tpi/u1i(ij,l))
  	      el3(ij,l,j)    = 0.0
  	      direct(ij,l,j) = 0.0
  	      ee3(ij,l,j)    = 0.0
  	  END IF
  	END DO
      END DO
    END DO
    
    DO  l= nsolp+1,ntotal
      DO ij=1,iend
  	IF(irs /= 0)  THEN
  	  sfcs(ij,l)= emis(l)*ptempg(ij,l)*pi
  	END IF
      END DO
    END DO
   
    j=  0
    DO  jd=  2,jn,2
     j=  j + 1
     DO  l=1,ntotal
       DO ij=1,iend
  	 IF(isl_aerad(ij) /= 0 .OR. irs .NE. 0 ) THEN
  	   IF(l>=lls(ij) .AND. l<=lla(ij)) THEN
  	 !	    HERE ARE THE EVEN MATRIX ELEMENTS
  	   df(ij,l,jd) = (cp(ij,l,j+1) - cpb(ij,l,j))*em1(ij,l,j+1) -  &
  		(cm(ij,l,j+1) - cmb(ij,l,j))*em2(ij,l,j+1)
  	 !	    HERE ARE THE ODD MATRIX ELEMENTS EXCEPT FOR THE TOP.
  	   df(ij,l,jd+1) =  el2(ij,l,j) * (cp(ij,l,j+1)-cpb(ij,l,j)) +  &
  		el1(ij,l,j) * (cmb(ij,l,j) - cm(ij,l,j+1))
  	    END IF
  	  END IF
  	END DO
      END DO
    END DO
   
    !	  HERE ARE THE TOP AND BOTTOM BOUNDARY CONDITIONS AS WELL AS THE
    !	  BEGINNING OF THE TRIDIAGONAL SOLUTION DEFINITIONS. I ASSUME NO
    !	  DIFFUSE RADIATION IS INCIDENT AT THE TOP.
    DO  l=1,ntotal
      DO ij=1,iend
  	IF(isl_aerad(ij) /= 0 .OR. irs .NE. 0 ) THEN
  	  IF(l>=lls(ij).AND. l<=lla(ij)) THEN
  	    df(ij,l,1)   = -cm(ij,l,1)
  	    df(ij,l,jdble) = sfcs(ij,l)+rsfx(ij,l)*cmb(ij,l,nlayer)- &
  			  cpb(ij,l,nlayer)
  	    ds(ij,l,jdble) = df(ij,l,jdble)/bf(ij,l,jdble)
  	    as(ij,l,jdble) = af(ij,l,jdble)/bf(ij,l,jdble)
  	  END IF
  	END IF
      END DO
    END DO
   
    !	  ********************************************
    !	  *	WE SOLVE THE TRIDIAGONAL EQUATIONS   *
    !	  ********************************************
   
    DO  j = 2, jdble
      DO  l=1,ntotal
  	DO ij=1,iend
  	  IF(isl_aerad(ij) /= 0 .OR. irs .NE. 0 ) THEN
  	    IF(l>=lls(ij) .AND. l<=lla(ij)) THEN
  	      x  = 1./(bf(ij,l,jdble+1-j) - ef(ij,l,jdble+1-j)* &
  			  as(ij,l,jdble+2-j))
  	      as(ij,l,jdble+1-j) = af(ij,l,jdble+1-j)*x
  	      ds(ij,l,jdble+1-j) = (df(ij,l,jdble+1-j) - &
  			  ef(ij,l,jdble+1-j) *ds(ij,l,jdble+2-j))*x
  	    END IF
  	  END IF
  	END DO
      END DO
    END DO
   
    DO  l=1,ntotal
      DO ij=1,iend
  	IF(isl_aerad(ij) /= 0 .OR. irs .NE. 0 ) THEN
  	  IF(l>=lls(ij) .AND. l<=lla(ij)) THEN
  	    xk(ij,l,1)    = ds(ij,l,1)
  	  END IF 
  	END IF
      END DO
    END DO
  
    DO  j	= 2, jdble
      DO  l=1,ntotal
  	DO ij=1,iend
  	  IF(isl_aerad(ij) /= 0 .OR. irs .NE. 0 ) THEN
  	    IF(l>=lls(ij) .AND. l<=lla(ij)) THEN
  	      xk(ij,l,j) = ds(ij,l,j) - as(ij,l,j)*xk(ij,l,j-1)
  	    END IF 
  	  END IF
  	END DO
      END DO
    END DO
   
    !  ***************************************************************
    !	  CALCULATE LAYER COEFFICIENTS, NET FLUX AND MEAN INTENSITY
    !  ***************************************************************
      
     DO j = 1,nlayer
       DO  l=1,ntotal
  	 DO ij=1,iend
  	   IF(isl_aerad(ij) /= 0 .OR. irs .NE. 0 ) THEN
  	     IF(l>=lls(ij) .AND. l<=lla(ij)) THEN
  	       ck1(ij,l,j)   = xk(ij,l,2*j-1)
  	       ck2(ij,l,j)   = xk(ij,l,2*j)
	       
  	       fnet(ij,l,j)  = ck1(ij,l,j)  *( el1(ij,l,j) &
  				 -el2(ij,l,j)) + ck2(ij,l,j) * &
  				    ( em1(ij,l,j)-em2(ij,l,j) ) + &
  				    cpb(ij,l,j) - cmb(ij,l,j) - direct(ij,l,j)
 	       tmi(ij,l,j)   =  el3(ij,l,j) + u1i(ij,l) *(ck1(ij,l,j)  *  &
  			  ( el1(ij,l,j) + el2(ij,l,j))   + ck2(ij,l,j) * &
  			  ( em1(ij,l,j)+em2(ij,l,j) ) +  cpb(ij,l,j) + &
  			  cmb(ij,l,j) )
  	     END IF 
  	   END IF
  	 END DO
       END DO
     END DO
   
    END SUBROUTINE add
    
    SUBROUTINE newflux1(m1,ia,iz,ja,jz)
      !
      !     **************************************************************
      !     *  Purpose  	   :  Calculate upward and downward	 *
      !     *			      intensities and fluxes using Gauss *
      !     *			      Quadrature angles and weights.	 *
      !     *  Subroutines Called  :  None				 *
      !     *  Input		   :  PTEMP, SLOPE, Y3, B3, EE1, EE2	 *
      !     *  Output		   :  DINTENT, UINTENT, DIREC, DIRECTU   *
      !     * ************************************************************
      !
      !INCLUDE 'globrad.h'
      USE mem_globrad, ONLY: ntotal,ngauss,nlayer,nsolp,tpi, &
  			     irs,gangle, &
  			     iblackbody_above, &
  			     gratio,gweight,emis
      
      IMPLICIT NONE
      
      INTEGER,INTENT(IN) :: m1,ia,iz,ja,jz
      INTEGER :: i
      INTEGER :: j
      INTEGER :: kindex
      INTEGER :: l
      INTEGER :: m
      REAL    :: ckm
      REAL    :: ckp
      REAL    :: x4
      REAL    :: ya
      REAL    :: yb
  	
      !
      REAL,DIMENSION((iz-ia+1)*(jz-ja+1),ntotal,ngauss,nlayer) :: y1,y2,y4,y8
      REAL,DIMENSION((iz-ia+1)*(jz-ja+1),ntotal,ngauss,nlayer) :: dintent,uintent
      REAL,DIMENSION((iz-ia+1)*(jz-ja+1),ntotal,nlayer) :: a1,a2,a3,a4,a7
      REAL,DIMENSION((iz-ia+1)*(jz-ja+1),ntotal,nlayer) :: y5
      INTEGER :: ij,iend,i1,j1
      
      iend=(iz-ia+1)*(jz-ja+1)
      !
    
      DO  j	  =  1,nlayer
  	kindex     = MAX( 1, j-1 )
  	DO   l  =  nsolp+1,ntotal
  	  DO ij=1,iend
  	    !HERE WE DO NO SCATTERING COEFFICIENTS
  	    a3(ij,l,j) =  ptemp(ij,l,kindex)*tpi
  	    a4(ij,l,j) =  tpi*slope(ij,l,j)
  	    a7(ij,l,j) =  a3(ij,l,j)
  	    y5(ij,l,j) =  a4(ij,l,j)*taul(ij,l,j)
  	  END DO
      END DO
      !HERE WE DO SCATTERING
      DO  l    =  nsolp+1,ntotal
  	DO ij=1,iend
  	  IF(irs /= 0) THEN
  	    x4         =  slope(ij,l,j)*(tpi*b3(ij,l,j)-(tpi/u1i(ij,l)))
  	    a1(ij,l,j) = u1i(ij,l) - ak(ij,l,j)
  	    a2(ij,l,j) = gami(ij,l,j)*(ak(ij,l,j)+u1i(ij,l))
  	    a3(ij,l,j) = a3(ij,l,j)+x4
  	    a7(ij,l,j) = a7(ij,l,j)-x4
  	  END IF
  	END DO
      END DO
    END DO
    !
    !	  CALCULATIONS FOR ALL GAUSS POINTS. HERE WE DO NO SCATTERING COEFFI
    !
    DO j=  1,nlayer
      DO i=  1,ngauss
  	DO  l =  nsolp+1,ntotal
  	  DO ij=1,iend
  	    y1(ij,l,i,j)  =  0.0
  	    y2(ij,l,i,j)  =  0.0
  	    y4(ij,l,i,j)  =  a7(ij,l,j) - a4(ij,l,j)*gangle(i)
  	    y8(ij,l,i,j)  =  a3(ij,l,j) + a4(ij,l,j)*gangle(i)
  	  END DO
  	END DO
  	!HERE WE DO SCATTERING
  	DO  l =  nsolp+1,ntotal
  	  DO ij=1,iend
  	    IF(irs /= 0) THEN
  	      ya=  a1(ij,l,j)*(y3(ij,l,i,j)-ee1(ij,l,j))/ &
  		      (ak(ij,l,j)*gangle(i)-1.)
  	      yb=  a2(ij,l,j)*(1.- ee1(ij,l,j)*y3(ij,l,i,j))/ &
  			(ak(ij,l,j)*gangle(i)+1.)
  	      ckp= ck1(ij,l,j)+ck2(ij,l,j)
  	      ckm= ck1(ij,l,j) -ck2(ij,l,j)
  	      y1(ij,l,i,j) =  ckp*yb+ckm*ya
  	      y2(ij,l,i,j) =  ckp*ya+ ckm*yb
  	    END IF
  	  END DO
  	END DO
      END DO
    END DO
    !
    DO  j	  =  1,nlayer
      DO   l	 =  nsolp+1,ntotal
  	DO ij=1,iend
  	  !tmid(ij,l,j) = 0.0
  	  !tmiu(ij,l,j) = 0.0
  	  direc(ij,l,j)     =  0.0
  	  directu(ij,l,j)   =  0.0
  	END DO
      END DO
    END DO
    !
    !	  DIREC IS DOWNWARD FLUX. DIRECTU IS UPWARD FLUX.
    !	  CALCULATE DINTENT THE DOWNWARD INTENSITY AND DIREC THE DOWNWARD FL
    !
    DO  i = 1,ngauss
      DO  l = nsolp+1,ntotal
  	DO ij=1,iend
  	  IF( iblackbody_above == 1 )THEN
  	    dintent(ij,l,i,1) = ptempt(ij,l)*y3(ij,l,i,1)*tpi + &
  		      y1(ij,l,i,1)+ (1.-y3(ij,l,i,1))*y4(ij,l,i,1)
  	  ELSE
  	    dintent(ij,l,i,1) = (1.-y3(ij,l,i,1))*y4(ij,l,i,1) + &
  		       y1(ij,l,i,1)
  	  END IF
  	  !tmid(ij,l,1) = tmid(ij,l,1)+dintent(ij,l,i,1)*gratio(i)
  	  direc(ij,l,1)= direc(ij,l,1)+dintent(ij,l,i,1)* gweight(i)
  	END DO
      END DO
    END DO
    !
    !	   DINTENT IS DOWNWARD INTENSITY * TPI. DIREC IS THE DOWNWARD FLUX.
    !
    DO j= 2,nlayer
      DO i = 1,ngauss
  	DO l = nsolp+1,ntotal
  	  DO ij=1,iend
  	    dintent(ij,l,i,j)  = dintent(ij,l,i,j-1)*y3(ij,l,i,j) + &
  			   y1(ij,l,i,j)+y5(ij,l,j)+  &
  			     (1.-y3(ij,l,i,j))*y4(ij,l,i,j)
  	    !tmid(ij,l,j)= tmid(ij,l,j)  +dintent(ij,l,i,j)*gratio(i)
  	    direc(ij,l,j)= direc(ij,l,j)+dintent(ij,l,i,j)*gweight(i)
  	  END DO
  	END DO
      END DO
    END DO
    !
    !	  UINTENT IS THE UPWARD INTENSITY * TPI. DIRECTU IS THE UPWARD FLUX.
    !	  ASSUME THAT THE REFLECTIVITY IS LAMBERT.
    !
    DO i =  1,ngauss
      DO l =  nsolp+1,ntotal
  	DO ij=1,iend
  	  uintent(ij,l,i,nlayer)  =  ptempg(ij,l)*emis(l) *tpi+2.* &
  			      rsfx(ij,l)*direc(ij,l,nlayer)
  	  !tmiu(ij,l,nlayer)=  tmiu(ij,l,nlayer)+ &
  	!		     uintent(ij,l,i,nlayer)*gratio(i)
  	  directu(ij,l,nlayer)    =  directu(ij,l,nlayer)+ &
  			      uintent(ij,l,i,nlayer)*gweight(i)
  	END DO
      END DO
    END DO
    !
    DO m= 2,nlayer
      j = nlayer-m+1
      DO i = 1,ngauss
  	DO l = nsolp+1,ntotal
  	  DO ij=1,iend
  	    uintent(ij,l,i,j)= (uintent(ij,l,i,j+1)-y5(ij,l,j+1)) * &
  			  y3(ij,l,i,j+1)+y2(ij,l,i,j+1)+ &
  			(1.-y3(ij,l,i,j+1))*y8(ij,l,i,j+1)
  	    !tmiu(ij,l,j) = tmiu(ij,l,j)+uintent(ij,l,i,j)*gratio(i)
  	    directu(ij,l,j) = directu(ij,l,j) + gweight(i)* &
  			       uintent(ij,l,i,j)
  	  END DO
  	END DO
      END DO
    END DO
    
  END SUBROUTINE newflux1
  
  
  SUBROUTINE plnk(e,t1,d)
    
    !	  ******************************************************
    !	  *  Purpose		 :  Calculate Planck Function  *
    !	  *  Subroutines Called  :  None		       *
    !	  *  Input		 :  WAVE, NCOUNT	       *
    !	  *  Output		 :  PLANK		       *
    !	  * ****************************************************
   
    !  THIS SUBROUTINE COMPUTES THE INTEGRAL OF THE PLANCK FUNCTION BETWEEN
    !  ZERO AND THE SPECIFIED VALUE OF LAMBDA.  THUS (USING XL AS LAMBDA)
    !  WE WANT TO INTEGRATE
    !  R = INTEGRAL(XL=0 TO XL=XLSPEC) ( C1*XL**-5* / (EXP(C2/XL*T)-1) )*DXL
    !  SUBSTITUTING U=C2/(XL*T), THE INTEGRAL BECOMES
    !  R = A CONSTANT TIMES INTEGRAL (USPEC TO INFINITY) OF
    !		 ( U**3 / (EXP(U) - 1) )*DU
    !  THE APPROXIMATIONS SHOWN HERE ARE ON PAGE 998 OF ABRAMOWITZ AND SEGUN
    !  UNDER THE HEADING OF DEBYE FUNCTIONS.  C2 IS THE PRODUCT OF PLANCK'S
    !  CONSTANT AND THE SPEED OF LIGHT DIVIDED BY BOLTZMANN'S CONSTANT.
    !  C2 = 14390 WHEN LAMBDA IS IN MICRONS.
    !  THE FACTOR 0.15399 IS THE RECIPROCAL OF SIX TIMES
    !  THE SUM OF (1/N**2) FOR ALL N FROM ONE TO INFINITY.  IT IS CHOSEN TO
    !  NORMALIZE THE INTEGRAL TO A MAXIMUM VALUE OF UNITY.
    !  RADIATION IN REAL UNITS IS OBTAINED BY MULTIPLYING THE INTEGRAL BY
    !  THE STEFAN-BOLTZMANN CONSTANT TIMES T**4.
    IMPLICIT NONE
   
    REAL				     :: e
    REAL, INTENT(IN)			     :: t1
    REAL, INTENT(OUT)			     :: d
    REAL :: am(5)
    REAL :: v1,a
    INTEGER :: m
   
    d		 =   0.0
    v1  	 =   e/t1
   
    IF (v1 <= 1.) THEN
      d 	=  1.0 - 0.15399*v1**3 *  &
  	  (1./3.-v1/8. + v1**2/60. - v1**4/5040. +  &
  	  v1**6/272160. - v1**8/13305600	 )
    END IF
   
    IF ( v1 > 1. .AND. v1 <= 50.) THEN
      DO  m   =  1,5
  	a	=  FLOAT(m)*v1
  	am(m)	=  0.15399 * EXP(-a)/m**4 * (((a+3.)*a+6.)*a+6.)
      END DO
   
      d 	 =  am(1)+am(2)+am(3)+am(4)+am(5)
    END IF
   
    d		  =  d*t1**4
   
  END SUBROUTINE plnk
  
  
  
  SUBROUTINE miess( ro, rfr, rfi, thetd, jx, qext, qscat,  &
  	  ctbrqs, r, re2, tmag2, wvno  )
    use rconstants, only: pio180
    !
    ! **********************************************************************
    !	 THIS SUBROUTINE COMPUTES MIE SCATTERING BY A STRATIFIED SPHERE,
    !	 I.E. A PARTICLE CONSISTING OF A SPHERICAL CORE SURROUNDED BY A
    !	 SPHERICAL SHELL.  THE BASIC CODE USED WAS THAT DESCRIBED IN THE
    !	 REPORT: " SUBROUTINES FOR COMPUTING THE PARAMETERS OF THE
    !	 ELECTROMAGNETIC RADIATION SCATTERED BY A SPHERE " J.V. DAVE,
    !	 I B M SCIENTIFIC CENTER, PALO ALTO , CALIFORNIA.
    !	 REPORT NO. 320 - 3236 .. MAY 1968 .
    !
    !	 THE MODIFICATIONS FOR STRATIFIED SPHERES ARE DESCRIBED IN
    !	     TOON AND ACKERMAN, APPL. OPTICS, IN PRESS, 1981
    !
    !	 THE PARAMETERS IN THE CALLING STATEMENT ARE DEFINED AS FOLLOWS :
    !	   RO IS THE OUTER (SHELL) RADIUS;
    !	   R  IS THE CORE RADIUS;
    !	   RFR, RFI  ARE THE REAL AND IMAGINARY PARTS OF THE SHELL INDEX
    !	       OF REFRACTION IN THE FORM (RFR - I* RFI);
    !	   RE2, TMAG2  ARE THE INDEX PARTS FOR THE CORE;
    !	       ( WE ASSUME SPACE HAS UNIT INDEX. )
    !	   THETD(J): ANGLE IN DEGREES BETWEEN THE DIRECTIONS OF THE INCIDENT
    !	       AND THE SCATTERED RADIATION.  THETD(J) IS< OR= 90.0
    !	       IF THETD(J) SHOULD HAPPEN TO BE GREATER THAN 90.0, ENTER WITH
    !	       SUPPLEMENTARY VALUE, SEE COMMENTS BELOW ON ELTRMX;
    !	   JX: TOTAL NUMBER OF THETD FOR WHICH THE COMPUTATIONS ARE
    !	       REQUIRED.  JX SHOULD NOT EXCEED IT UNLESS THE DIMENSIONS
    !	       STATEMENTS ARE APPROPRIATEDLY MODIFIED;
    !
    !	   THE DEFINITIONS FOR THE FOLLOWING SYMBOLS CAN BE FOUND IN"LIGHT
    !	       SCATTERING BY SMALL PARTICLES,H.C.VAN DE HULST, JOHN WILEY '
    !	       SONS, INC., NEW YORK, 1957" .
    !	   QEXT: EFFICIENCY FACTOR FOR EXTINCTION,VAN DE HULST,P.14 ' 127.
    !	   QSCAT: EFFICIENCY FACTOR FOR SCATTERING,V.D. HULST,P.14 ' 127.
    !	   CTBRQS: AVERAGE(COSINE THETA) * QSCAT,VAN DE HULST,P.128
    !	   ELTRMX(I,J,K): ELEMENTS OF THE TRANSFORMATION MATRIX F,V.D.HULST
    !	       ,P.34,45 ' 125. I=1: ELEMENT M SUB 2..I=2: ELEMENT M SUB 1..
    !	       I = 3: ELEMENT S SUB 21.. I = 4: ELEMENT D SUB 21..
    !	   ELTRMX(I,J,1) REPRESENTS THE ITH ELEMENT OF THE MATRIX FOR
    !	       THE ANGLE THETD(J).. ELTRMX(I,J,2) REPRESENTS THE ITH ELEMENT
    !	       OF THE MATRIX FOR THE ANGLE 180.0 - THETD(J) ..
    !	   QBS IS THE BACK SCATTER CROSS SECTION.
    !
    !	   IT: IS THE DIMENSION OF THETD, ELTRMX, CSTHT, PI, TAU, SI2THT,
    !	       IT MUST CORRESPOND EXACTLY TO THE SECOND DIMENSION OF ELTRMX.
    !	   IACAP IS THE DIMENSION OF ACAP
    !	       IN THE ORIGINAL PROGRAM THE DIMENSION OF ACAP WAS 7000.
    !	       FOR CONSERVING SPACE THIS SHOULD BE NOT MUCH HIGHER THAN
    !	       THE VALUE, N=1.1*(NREAL**2 + NIMAG**2)**.5 * X + 1
    !	   WVNO: 2*PI / WAVELENGTH
    !
    !	 ALSO THE SUBROUTINE COMPUTES THE CAPITAL A FUNCTION BY MAKING USE O
    !	 DOWNWARD RECURRENCE RELATIONSHIP.
    !
    !	   TA(1): REAL PART OF WFN(1).  TA(2): IMAGINARY PART OF WFN(1).
    !	   TA(3): REAL PART OF WFN(2).  TA(4): IMAGINARY PART OF WFN(2).
    !	   TB(1): REAL PART OF FNA.	TB(2): IMAGINARY PART OF FNA.
    !	   TC(1): REAL PART OF FNB.	TC(2): IMAGINARY PART OF FNB.
    !	   TD(1): REAL PART OF FNAP.	TD(2): IMAGINARY PART OF FNAP.
    !	   TE(1): REAL PART OF FNBP.	TE(2): IMAGINARY PART OF FNBP.
    !	   FNAP, FNBP  ARE THE PRECEDING VALUES OF FNA, FNB RESPECTIVELY.
    ! **********************************************************************
    !
    !
    !	Include implicit declarations
    !
    !INCLUDE 'precision.h'
    !
    !
    !	Define dimensions of local arrays and arrays passed as arguments
    !
    IMPLICIT NONE
    
    INTEGER, PARAMETER :: iacap = 200000
    INTEGER, PARAMETER :: it = 1
    REAL, INTENT(IN)			     :: ro
    REAL, INTENT(IN)			     :: rfr
    REAL, INTENT(IN)			     :: rfi
    INTEGER, INTENT(IN) 		     :: jx
    REAL, INTENT(IN OUT)		     :: thetd(jx)
    REAL, INTENT(OUT)			     :: qext
    REAL, INTENT(OUT)			     :: qscat
    REAL, INTENT(OUT)			     :: ctbrqs
    REAL, INTENT(IN)			     :: r
    REAL, INTENT(IN)			     :: re2
    REAL, INTENT(IN)			     :: tmag2
    REAL, INTENT(IN)			     :: wvno
    !
    !
    !	Declare arguments passed as arrays
    !
   
    !
    !
    !	Declare local variables
    !
    DOUBLE PRECISION, PARAMETER :: epsilon_mie = 1.d-14
   
    DOUBLE COMPLEX :: fnap,   fnbp,   acap(iacap),  &
  	fna,	fnb,	rf,	  rrf, rrfx,   wm1,    fn1,	 fn2,  &
  	tc1,	tc2,	wfn(2),   z(4), k1,	k2,	k3,	  w(3,iacap),  &
  	rc,	u(8),	dh1, dh2,    dh4,    p24h24,   p24h21,  &
  	pstore, hstore, dummy,    dumsq
   
    DOUBLE PRECISION :: t(5), ta(4), tb(2), tc(2), td(2), te(2),  &
  	pi(3,it), tau(3,it), cstht(it), si2tht(it), eltrmx(4,it,2),  &
  	x, x1, x4, y1, y4, rx, sinx1, sinx4, cosx1, cosx4,  &
  	ey1, e2y1, ey4, ey1my4, ey1py4, aa, bb, cc, dd, denom,  &
  	realp, amagp, qbsr, qbsi, rmm
   
    !INTEGER, external :: imag
    INTEGER :: iflag
    INTEGER :: nmx1
    INTEGER :: nmx2
    INTEGER :: n
    INTEGER :: nn
    INTEGER :: m
    INTEGER :: j
    INTEGER :: k
    INTEGER :: i
  
    EQUIVALENCE (fna,tb(1)),(fnb,tc(1)),(fnap,td(1)),(fnbp,te(1))
    !
    !
    !  Some compilers (e.g. absoft) don't the support imag(z) generic intrinsic function.
    !  In that situation, uncomment the following 2 lines to define stmt function to
    !  redefine imag() function to the alternative function.  I.e. change "alt_function"
    !  in following stmt function to the appropriate function that will return the
    !  imaginary part of a double complex for the compiler being used.
    !  For compilers that support imag(z), simply leave following 2 statments commented.
    !  The "c--alt_imag" comment prefix is designed to be recognized by the
    !  automated editing script create_dmiess, so please do not modify this prefix
    !  in the dmiess.f.template file.
    !  -bm  Nov-1999
    !
    !DOUBLE COMPLEX :: z_dum_arg
   
    !imag(z_dum_arg) = DIMAG(z_dum_arg)
    !
   
    !	    print*, RO, RFR, RFI, THETD, JX, QEXT, QSCAT, CTBRQS,
    !	  1		     R, RE2, TMAG2, WVNO
   
    !
    !	IF THE CORE IS SMALL SCATTERING IS COMPUTED FOR THE SHELL ONLY
    !
    iflag = 1
    IF ( r/ro < 1.d-6 )   iflag = 2
    IF ( jx <= it )   GO TO 20
    WRITE( *,7 )
    WRITE( *,6 )
    STOP 30
    20 rf =  CMPLX( rfr,  -rfi )
    rc =  CMPLX( re2, -tmag2 )
    x  =  ro * wvno
    k1 =  rc * wvno
    k2 =  rf * wvno
    k3 =  CMPLX( wvno, 0.0 )
    z(1) =  k2 * ro
    z(2) =  k3 * ro
    z(3) =  k1 * r
    z(4) =  k2 * r
    x1   =  REAL( z(1) )
    x4   =  REAL( z(4) )
    y1   =  DIMAG( z(1) )
    y4   =  DIMAG( z(4) )
    !	   print*,'Z(1)','Z(4)','x1','x4','y1','y4'
    !	   print*,Z(1),Z(4),x1,x4,y1,y4
   
    rrf  =  1.0 / rf
    rx   =  1.0 / x
    rrfx =  rrf * rx
    t(1) =  ( x**2 ) * ( rfr**2 + rfi**2 )
    t(1) =  SQRT( t(1) )
    nmx1 =  1.10 * t(1)
    !
    IF ( nmx1 <= iacap-1 )   GO TO 21
    WRITE(*,8)
    STOP 32
    21 nmx2 = t(1)
    IF ( nmx1 >  150 )   GO TO 22
    nmx1 = 150
    nmx2 = 135
    !
    22 acap( nmx1+1 )  =  ( 0.0,0.0 )
    IF ( iflag == 2 )	GO TO 26
    DO n = 1,3
      w( n,nmx1+1 )  =  ( 0.0,0.0 )
    END DO
    26 CONTINUE
    DO n = 1,nmx1
      nn = nmx1 - n + 1
      acap(nn) = (nn+1) * rrfx - 1.0 / ( (nn+1) * rrfx + acap(nn+1) )
      IF ( iflag == 2 )   GO TO 23
      DO m = 1,3
  	w( m,nn ) = (nn+1) / z(m+1)  - 1.0 / (  (nn+1) / z(m+1)  +  w( m,nn+1 )  )
      END DO
      23 CONTINUE
    END DO
    !
    DO    j = 1,jx
      IF ( thetd(j) < 0.0 )  thetd(j) =  ABS( thetd(j) )
      IF ( thetd(j) > 0.0 )  GO TO 24
      cstht(j)  = 1.0
      si2tht(j) = 0.0
      CYCLE
      24 IF ( thetd(j) >= 90.0 )  GO TO 25
      t(1)	=  pio180 * thetd(j)
      cstht(j)  =  COS( t(1) )
      si2tht(j) =  1.0 - cstht(j)**2
      CYCLE
      25 IF ( thetd(j) > 90.0 )  GO TO 28
      cstht(j)  =  0.0
      si2tht(j) =  1.0
      CYCLE
      28 WRITE( *,5 )  thetd(j)
      WRITE( *,6 )
      STOP 34
    END DO
    !
    DO   j = 1,jx
      pi(1,j)  =  0.0
      pi(2,j)  =  1.0
      tau(1,j) =  0.0
      tau(2,j) =  cstht(j)
    END DO
    !
    ! INITIALIZATION OF HOMOGENEOUS SPHERE
    !
    t(1)   =  COS(x)
    t(2)   =  SIN(x)
    wm1    =  CMPLX( t(1),-t(2) )
    wfn(1) =  CMPLX( t(2), t(1) )
    ta(1)  =  t(2)
    ta(2)  =  t(1)
    wfn(2) =  rx * wfn(1) - wm1
    ta(3)  =  REAL(wfn(2))
    ta(4)  =  DIMAG(wfn(2))
    !	   print*,'WFN(2)','TA(3)','TA(4)'
    !	   print*,WFN(2),TA(3),TA(4)
   
    !
    IF ( iflag == 2 )	GO TO 560
    n = 1
    !
    ! INITIALIZATION PROCEDURE FOR STRATIFIED SPHERE BEGINS HERE
    !
    sinx1   =  SIN( x1 )
    sinx4   =  SIN( x4 )
    cosx1   =  COS( x1 )
    cosx4   =  COS( x4 )
    ey1     =  EXP( y1 )
    e2y1    =  ey1 * ey1
    ey4     =  EXP( y4 )
    ey1my4  =  EXP( y1 - y4 )
    ey1py4  =  ey1 * ey4
    ey1my4  =  EXP( y1 - y4 )
    aa  =  sinx4 * ( ey1py4 + ey1my4 )
    bb  =  cosx4 * ( ey1py4 - ey1my4 )
    cc  =  sinx1 * ( e2y1 + 1.0 )
    dd  =  cosx1 * ( e2y1 - 1.0 )
    denom   =  1.0  +  e2y1 * ( 4.0 * sinx1 * sinx1 - 2.0 + e2y1 )
    realp   =  ( aa * cc  +  bb * dd ) / denom
    amagp   =  ( bb * cc  -  aa * dd ) / denom
    dummy   =  CMPLX( realp, amagp )
    aa  =  sinx4 * sinx4 - 0.5
    bb  =  cosx4 * sinx4
    p24h24  =  0.5 + CMPLX( aa,bb ) * ey4 * ey4
    aa  =  sinx1 * sinx4  -  cosx1 * cosx4
    bb  =  sinx1 * cosx4  +  cosx1 * sinx4
    cc  =  sinx1 * sinx4  +  cosx1 * cosx4
    dd  = -sinx1 * cosx4  +  cosx1 * sinx4
    p24h21  =  0.5 * CMPLX( aa,bb ) * ey1 * ey4  + 0.5 * CMPLX( cc,dd ) * ey1my4
    dh4  =  z(4) / ( 1.0 + ( 0.0,1.0 ) * z(4) )  -  1.0 / z(4)
    dh1  =  z(1) / ( 1.0 + ( 0.0,1.0 ) * z(1) )  -  1.0 / z(1)
    dh2  =  z(2) / ( 1.0 + ( 0.0,1.0 ) * z(2) )  -  1.0 / z(2)
    pstore  =  ( dh4 + n / z(4) )  *  ( w(3,n) + n / z(4) )
    p24h24  =  p24h24 / pstore
    hstore  =  ( dh1 + n / z(1) )  *  ( w(3,n) + n / z(4) )
    p24h21  =  p24h21 / hstore
    pstore  =  ( acap(n) + n / z(1) )  /  ( w(3,n) + n / z(4) )
    dummy   =  dummy * pstore
    dumsq   =  dummy * dummy
    !
    ! NOTE:  THE DEFINITIONS OF U(I) IN THIS PROGRAM ARE NOT THE SAME AS
    !	     THE USUBI DEFINED IN THE ARTICLE BY TOON AND ACKERMAN.  THE
    !	     CORRESPONDING TERMS ARE:
    !	       USUB1 = U(1)			  USUB2 = U(5)
    !	       USUB3 = U(7)			  USUB4 = DUMSQ
    !	       USUB5 = U(2)			  USUB6 = U(3)
    !	       USUB7 = U(6)			  USUB8 = U(4)
    !	       RATIO OF SPHERICAL BESSEL FTN TO SPHERICAL HENKAL FTN = U(8)
    !
    u(1) =  k3 * acap(n)  -  k2 * w(1,n)
    u(2) =  k3 * acap(n)  -  k2 * dh2
    u(3) =  k2 * acap(n)  -  k3 * w(1,n)
    u(4) =  k2 * acap(n)  -  k3 * dh2
    u(5) =  k1 *  w(3,n)  -  k2 * w(2,n)
    u(6) =  k2 *  w(3,n)  -  k1 * w(2,n)
    u(7) =  ( 0.0,-1.0 )  *  ( dummy * p24h21 - p24h24 )
    u(8) =  ta(3) / wfn(2)
    !
    fna  =  u(8) * ( u(1)*u(5)*u(7)  +  k1*u(1)  -  dumsq*k3*u(5) ) /  &
  	( u(2)*u(5)*u(7)  +  k1*u(2)  -  dumsq*k3*u(5) )
    fnb  =  u(8) * ( u(3)*u(6)*u(7)  +  k2*u(3)  -  dumsq*k2*u(6) ) /  &
  	( u(4)*u(6)*u(7)  +  k2*u(4)  -  dumsq*k2*u(6) )
    GO TO 561
    560 tc1  =  acap(1) * rrf  +  rx
    tc2  =  acap(1) * rf   +  rx
    fna  =  ( tc1 * ta(3)  -  ta(1) ) / ( tc1 * wfn(2)  -  wfn(1) )
    fnb  =  ( tc2 * ta(3)  -  ta(1) ) / ( tc2 * wfn(2)  -  wfn(1) )
    !
    561 CONTINUE
    fnap = fna
    fnbp = fnb
    t(1) = 1.50
    !
    !	 FROM HERE TO THE STATMENT NUMBER 90, ELTRMX(I,J,K) HAS
    !	 FOLLOWING MEANING:
    !	 ELTRMX(1,J,K): REAL PART OF THE FIRST COMPLEX AMPLITUDE.
    !	 ELTRMX(2,J,K): IMAGINARY PART OF THE FIRST COMPLEX AMPLITUDE.
    !	 ELTRMX(3,J,K): REAL PART OF THE SECOND COMPLEX AMPLITUDE.
    !	 ELTRMX(4,J,K): IMAGINARY PART OF THE SECOND COMPLEX AMPLITUDE.
    !	 K = 1 : FOR THETD(J) AND K = 2 : FOR 180.0 - THETD(J)
    !	 DEFINITION OF THE COMPLEX AMPLITUDE: VAN DE HULST,P.125.
    !
    tb(1) = t(1) * tb(1)
    tb(2) = t(1) * tb(2)
    tc(1) = t(1) * tc(1)
    tc(2) = t(1) * tc(2)
    DO  j = 1,jx
      eltrmx(1,j,1) = tb(1) * pi(2,j) + tc(1) * tau(2,j)
      eltrmx(2,j,1) = tb(2) * pi(2,j) + tc(2) * tau(2,j)
      eltrmx(3,j,1) = tc(1) * pi(2,j) + tb(1) * tau(2,j)
      eltrmx(4,j,1) = tc(2) * pi(2,j) + tb(2) * tau(2,j)
      eltrmx(1,j,2) = tb(1) * pi(2,j) - tc(1) * tau(2,j)
      eltrmx(2,j,2) = tb(2) * pi(2,j) - tc(2) * tau(2,j)
      eltrmx(3,j,2) = tc(1) * pi(2,j) - tb(1) * tau(2,j)
      eltrmx(4,j,2) = tc(2) * pi(2,j) - tb(2) * tau(2,j)
    END DO
    !
    qext   = 2.0 * ( tb(1) + tc(1))
    qscat  = ( tb(1)**2 + tb(2)**2 + tc(1)**2 + tc(2)**2 ) / 0.75
    ctbrqs = 0.0
    qbsr   = -2.0*(tc(1) - tb(1))
    qbsi   = -2.0*(tc(2) - tb(2))
    rmm    = -1.0
    n = 2
    65 t(1) = 2*n - 1
    t(2) =   n - 1
    t(3) = 2*n + 1
    DO   j = 1,jx
      pi(3,j)  = ( t(1) * pi(2,j) * cstht(j) - n * pi(1,j) ) / t(2)
      tau(3,j) = cstht(j) * ( pi(3,j) - pi(1,j) )  -  &
  	  t(1) * si2tht(j) * pi(2,j)  +  tau(1,j)
    END DO
    !
    ! HERE SET UP HOMOGENEOUS SPHERE
    !
    wm1    =  wfn(1)
    wfn(1) =  wfn(2)
    ta(1)  =  REAL(wfn(1))
    ta(2)  =  DIMAG(wfn(1))
    ta(4)  =  DIMAG(wfn(2))
    wfn(2) =  t(1) * rx * wfn(1)  -  wm1
    ta(3)  =  REAL(wfn(2))
   
    !	   print*,'WFN(1)','TA(1)','TA(2)'
    !	   print*,WFN(1),TA(1),TA(2)
    !
    IF ( iflag == 2 )	GO TO 1000
    !
    ! HERE SET UP STRATIFIED SPHERE
    !
    dh2  =  - n / z(2)  +  1.0 / ( n / z(2) - dh2 )
    dh4  =  - n / z(4)  +  1.0 / ( n / z(4) - dh4 )
    dh1  =  - n / z(1)  +  1.0 / ( n / z(1) - dh1 )
    pstore  =  ( dh4 + n / z(4) )  *  ( w(3,n) + n / z(4) )
    p24h24  =  p24h24 / pstore
    hstore  =  ( dh1 + n / z(1) )  *  ( w(3,n) + n / z(4) )
    p24h21  =  p24h21 / hstore
    pstore  =  ( acap(n) + n / z(1) )  /  ( w(3,n) + n / z(4) )
    dummy   =  dummy * pstore
    dumsq   =  dummy * dummy
    !
    u(1) =  k3 * acap(n)  -  k2 * w(1,n)
    u(2) =  k3 * acap(n)  -  k2 * dh2
    u(3) =  k2 * acap(n)  -  k3 * w(1,n)
    u(4) =  k2 * acap(n)  -  k3 * dh2
    u(5) =  k1 *  w(3,n)  -  k2 * w(2,n)
    u(6) =  k2 *  w(3,n)  -  k1 * w(2,n)
    u(7) =  ( 0.0,-1.0 )  *  ( dummy * p24h21 - p24h24 )
    u(8) =  ta(3) / wfn(2)
    !
    fna  =  u(8) * ( u(1)*u(5)*u(7)  +  k1*u(1)  -  dumsq*k3*u(5) ) /  &
  	( u(2)*u(5)*u(7)  +  k1*u(2)  -  dumsq*k3*u(5) )
    fnb  =  u(8) * ( u(3)*u(6)*u(7)  +  k2*u(3)  -  dumsq*k2*u(6) ) /  &
  	( u(4)*u(6)*u(7)  +  k2*u(4)  -  dumsq*k2*u(6) )
    !
    1000 CONTINUE
    tc1  =  acap(n) * rrf  +  n * rx
    tc2  =  acap(n) * rf   +  n * rx
    fn1  =  ( tc1 * ta(3)  -  ta(1) ) /  ( tc1 * wfn(2) - wfn(1) )
    fn2  =  ( tc2 * ta(3)  -  ta(1) ) /  ( tc2 * wfn(2) - wfn(1) )
    m	 =  wvno * r
    IF ( n < m )   GO TO 1002
    IF ( iflag == 2 )	GO TO 1001
    IF (ABS((fn1-fna)/fn1) < epsilon_mie .AND.  &
  	ABS(  ( fn2-fnb ) / fn2  ) < epsilon_mie  ) iflag = 2
    IF ( iflag == 1 )	GO TO 1002
    1001 fna  =  fn1
    fnb  =  fn2
    !
    1002 CONTINUE
    t(5)  =  n
    t(4)  =  t(1) / ( t(5) * t(2) )
    t(2)  =  (  t(2) * ( t(5) + 1.0 )  ) / t(5)
    !
    ctbrqs  =  ctbrqs  +  t(2) * ( td(1) * tb(1)  +  td(2) * tb(2)  &
  	+	    te(1) * tc(1)  +  te(2) * tc(2) )  &
  	+  t(4) * ( td(1) * te(1)  +  td(2) * te(2) )
    qext    =	qext  +  t(3) * ( tb(1) + tc(1) )
    !	  $	   T(3), TB(1), TC(1), QEXT
    t(4)    =  tb(1)**2 + tb(2)**2 + tc(1)**2 + tc(2)**2
    qscat   =  qscat  +  t(3) * t(4)
    rmm     =  -rmm
    qbsr    =  qbsr + t(3)*rmm*(tc(1) - tb(1))
    qbsi    =  qbsi + t(3)*rmm*(tc(2) - tb(2))
    !
    t(2)    =  n * (n+1)
    t(1)    =  t(3) / t(2)
    k = (n/2)*2
    DO  j = 1,jx
      eltrmx(1,j,1) = eltrmx(1,j,1)+t(1)*(tb(1)*pi(3,j)+tc(1)*tau(3,j))
      eltrmx(2,j,1) = eltrmx(2,j,1)+t(1)*(tb(2)*pi(3,j)+tc(2)*tau(3,j))
      eltrmx(3,j,1) = eltrmx(3,j,1)+t(1)*(tc(1)*pi(3,j)+tb(1)*tau(3,j))
      eltrmx(4,j,1) = eltrmx(4,j,1)+t(1)*(tc(2)*pi(3,j)+tb(2)*tau(3,j))
      IF ( k == n )  THEN
  	eltrmx(1,j,2) =eltrmx(1,j,2)+t(1)*(-tb(1)*pi(3,j)+tc(1)*tau(3,j))
  	eltrmx(2,j,2) =eltrmx(2,j,2)+t(1)*(-tb(2)*pi(3,j)+tc(2)*tau(3,j))
  	eltrmx(3,j,2) =eltrmx(3,j,2)+t(1)*(-tc(1)*pi(3,j)+tb(1)*tau(3,j))
  	eltrmx(4,j,2) =eltrmx(4,j,2)+t(1)*(-tc(2)*pi(3,j)+tb(2)*tau(3,j))
      ELSE
  	eltrmx(1,j,2) = eltrmx(1,j,2)+t(1)*(tb(1)*pi(3,j)-tc(1)*tau(3,j))
  	eltrmx(2,j,2) = eltrmx(2,j,2)+t(1)*(tb(2)*pi(3,j)-tc(2)*tau(3,j))
  	eltrmx(3,j,2) = eltrmx(3,j,2)+t(1)*(tc(1)*pi(3,j)-tb(1)*tau(3,j))
  	eltrmx(4,j,2) = eltrmx(4,j,2)+t(1)*(tc(2)*pi(3,j)-tb(2)*tau(3,j))
      END IF
    END DO
    !
    IF ( t(4) < epsilon_mie )	GO TO 100
    n = n + 1
    DO  j = 1,jx
      pi(1,j)	=   pi(2,j)
      pi(2,j)	=   pi(3,j)
      tau(1,j)  =  tau(2,j)
      tau(2,j)  =  tau(3,j)
    END DO
    fnap  =  fna
    fnbp  =  fnb
    IF ( n <= nmx2 )   GO TO 65
    !	      print*,N,NMX2
    WRITE( *,8 )
    STOP 36
    100 CONTINUE
    DO j = 1,jx
      DO k = 1,2
  	DO    i= 1,4
  	  t(i)  =  eltrmx(i,j,k)
  	END DO
  	eltrmx(2,j,k)  =      t(1)**2  +  t(2)**2
  	eltrmx(1,j,k)  =      t(3)**2  +  t(4)**2
  	eltrmx(3,j,k)  =  t(1) * t(3)  +  t(2) * t(4)
  	eltrmx(4,j,k)  =  t(2) * t(3)  -  t(4) * t(1)
      END DO
    END DO
    t(1)    =	 2.0 * rx**2
    qext    =	qext * t(1)
    qscat   =  qscat * t(1)
    ctbrqs  =  2.0 * ctbrqs * t(1)
   
   
   
    !	   RO IS THE OUTER (SHELL) RADIUS;
    !	   R  IS THE CORE RADIUS;
    !	   RFR, RFI  ARE THE REAL AND IMAGINARY PARTS OF THE SHELL INDEX
    !	       OF REFRACTION IN THE FORM (RFR - I* RFI);
    !	   RE2, TMAG2  ARE THE INDEX PARTS FOR THE CORE;
   
    !
    ! QBS IS THE BACK SCATTER CROSS SECTION
    !
    !	   PIG   = ACOS(-1.0)
    !	   RXP4  = RX*RX/(4.0*PIG)
    !	   QBS   = RXP4*(QBSR**2 + QBSI**2)
    !
    5  FORMAT( 10X,' THE VALUE OF THE SCATTERING ANGLE IS GREATER THAN 90.0 DEGREES. IT IS ', e15.4 )
    6  FORMAT( // 10X, 'PLEASE READ COMMENTS.' // )
    7  FORMAT( // 10X, 'THE VALUE OF THE ARGUMENT JX IS GREATER THAN IT'//)
    8  FORMAT( // 10X, 'THE UPPER LIMIT FOR ACAP IS NOT ENOUGH. SUGGEST GET DETAILED OUTPUT AND MODIFY SUBROUTINE' // )
    !
  END SUBROUTINE miess
  
  
  SUBROUTINE  radtran_to_rams(m1,m2,m3,fthrl,rlong,fthrs,rshort,aotr,ia,iz,ja,jz,mynum)
  
    USE mem_grid   , ONLY: nzpmax	  !INTENT(IN)
    USE mem_globrad, ONLY: nwave,ntotal,nprob
    
    IMPLICIT NONE
  
    INTEGER,INTENT(IN)  	       :: m1,m2,m3,ia,iz,ja,jz,mynum
    REAL,INTENT(OUT)  ,DIMENSION((iz-ia+1)*(jz-ja+1))           :: rshort
    REAL,INTENT(OUT)  ,DIMENSION((iz-ia+1)*(jz-ja+1))           :: rlong
    REAL,INTENT(OUT)  ,DIMENSION((iz-ia+1)*(jz-ja+1),nwave)     :: aotr
    REAL,INTENT(INOUT),DIMENSION((iz-ia+1)*(jz-ja+1),nzpmax)	:: fthrl
    REAL,INTENT(INOUT),DIMENSION((iz-ia+1)*(jz-ja+1),nzpmax)	:: fthrs
    INTEGER :: ij,iend
  
    !Local
  
    REAL, DIMENSION((iz-ia+1)*(jz-ja+1),nzpmax)         :: duml
    REAL, DIMENSION((iz-ia+1)*(jz-ja+1),nzpmax)         :: dums
    REAL, DIMENSION((iz-ia+1)*(jz-ja+1),nwave)          :: dumaot
    REAL, DIMENSION((iz-ia+1)*(jz-ja+1),ntotal,nzpmax)  :: dum2aot
    
  
    INTEGER :: k,j1,i1
    INTEGER :: k1
    INTEGER :: kr
    INTEGER :: l
    INTEGER :: nzz
  
    iend=(iz-ia+1)*(jz-ja+1)
    !nzz = Vertical level number
    nzz = m1 - 1
    aotr=0.0
       
!srf- otimizado
!    DO k = 1,nzz
!       !  Reverse the vertical index when in cartesian coordinates
!      !
!       k1  = nzz + 1 - k
!       DO ij=1,iend
!	duml(ij,k1) = heati_aerad(ij,k)
!	dums(ij,k1) = heats_aerad(ij,k)
!       END DO
!     END DO
!  
!    ! Transfer values from CARMA grid to BRAMS grid
!     DO k=1,m1-1
!       kr = k + 1     ! K level in CARMA grid corresponds to K+1 level in BRAMS grid
!       DO ij=1,iend
!	 fthrl(ij,kr) = duml(ij,k)
!	 fthrs(ij,kr) = dums(ij,k)
!       END DO
!     END DO
   
    ! reverse the vertical and transfer values from CARMA grid to BRAMS grid
    DO k=2,m1
       kr = nzz+2- k 
       DO ij=1,iend
  	 fthrl(ij,k) = heati_aerad(ij,kr)
  	 fthrs(ij,k) = heats_aerad(ij,kr)
!	print*,k,fthrl(ij,k),fthrs(ij,k)
       END DO
    END DO
!srf


     DO ij=1,iend
  	rshort(ij) = solnet(ij)  ! total absorvido pela superficie
        rlong(ij)  = xirdown(ij) ! downward longwave  na superficie
        !if (rlong(ij).lt.10.) print*,'RADTRAN_TO_RAMS!!!',ij,rlong(ij)
    END DO
    
    dumaot=0.0
    !  DATA WAVE / 0.256, 0.280, 0.296, 0.319, 0.335, 0.365, 0.420,
    ! 2   0.690, 0.762, 0.719, 0.813, 0.862, 0.926, 1.005, 1.111,
    ! 3   1.333, 1.562, 1.770, 2.051, 2.210, 2.584, 3.284, 3.809,
    ! 4   4.292, 4.546, 4.878, 5.128, 5.405, 5.714, 6.061, 6.452,
    ! 5   6.897, 7.407, 8.333, 9.009, 10.309, 12.500, 13.889,	 
    ! 6   16.667,20.000, 26.316, 35.714, 62.50  		/

    DO l=1,ntotal
       DO k=1,m1				  
  	 DO ij=1,iend
           dum2aot(ij,nprob(l),k)=tauaer(ij,nprob(l),k)
  	 END DO 				  
       END DO					  
    END DO

    DO l=1,nwave
       DO k=1,m1				  
  	 DO ij=1,iend
  	   dumaot(ij,l) = dumaot(ij,l) + dum2aot(ij,l,k)   
!  	   dumaot(ij,l) = dumaot(ij,l) + tauaer(ij,l,k)   
  	 END DO 				  
       END DO					  
    END DO
    
    DO l=1,nwave
      DO ij=1,iend
  	aotr(ij,l)= dumaot(ij,l)
      END DO
    END DO
  
  END SUBROUTINE radtran_to_rams

!kmlnew  
  SUBROUTINE radcomp_carma(m1,m2,m3,ia,iz,ja,jz,solfac  &
       ,theta,pi0,pp,rv,RAIN,LWL,IWL,dn0,rtp,fthrd  &
       ,rtgt,f13t,f23t,glat,glon,rshort,rlong,albedt,cosz,rlongup  &
       ,mynum,fmapt,pm,patch_area,npat)
!kmlnew      
       USE mem_carma, ONLY: carma
       USE mem_grid , ONLY: ngrid

       ! For specific optimization depending the type of machine
!       use machine_arq, only: machine ! INTENT(IN)

       INTEGER,INTENT(IN) :: m1,m2,m3,ia,iz,ja,jz,mynum,npat
  
       REAL,INTENT(IN)    :: solfac
       REAL,INTENT(IN)    :: theta(m1,m2,m3)
       REAL,INTENT(IN)    :: pi0(m1,m2,m3)
       REAL,INTENT(IN)    :: pp(m1,m2,m3)
       REAL,INTENT(IN)    :: rv(m1,m2,m3)
!kmlnew
       REAL,INTENT(IN)    :: LWL(m1,m2,m3)
       REAL,INTENT(IN)    :: IWL(m1,m2,m3)
       REAL,INTENT(IN)    :: RAIN(m2,m3)
       REAL,INTENT(IN)    :: patch_area(m2,m3,npat)
!kmlnew      
       REAL,INTENT(IN)    :: dn0(m1,m2,m3)
       REAL,INTENT(IN)    :: rtp(m1,m2,m3)
       REAL,INTENT(IN)    :: rtgt(m2,m3)
       REAL,INTENT(IN)    :: f13t(m2,m3)
       REAL,INTENT(IN)    :: f23t(m2,m3)
       REAL,INTENT(IN)    :: glat(m2,m3)
       REAL,INTENT(IN)    :: glon(m2,m3)
       REAL,INTENT(IN)    :: cosz(m2,m3) 
       REAL,INTENT(IN)    :: albedt(m2,m3)
       REAL,INTENT(IN)    :: fmapt(m2,m3)
       REAL,INTENT(IN)    :: pm(m1,m2,m3) 
  
       REAL,INTENT(INOUT) :: rshort(m2,m3)
       REAL,INTENT(INOUT) :: rlong(m2,m3)
!kml       REAL,INTENT(INOUT) :: fthrd(m2,m3)
       REAL,INTENT(INOUT) :: fthrd(m1,m2,m3)
       
       REAL,INTENT(INOUT) :: rlongup(m2,m3)
       !kmlnew
       REAL :: xland(m2,m3)
       
       INTEGER :: ia1,iz1,ja1,jz1,ii,jj,iend,ipat
       
       iend=(iz-ia+1)*(jz-ja+1)

       ! For specific optimization
!       if (machine==1) then
          !-sx6 
!          p_isize=16;p_jsize=16
!       elseif(machine==0) then
          !-cluster ! Generic IA32
          p_isize=2;p_jsize=2
!       endif
       
!srf      DO ipat= 1,2
       DO ii=ia,iz,p_isize
  	 DO jj=ja,jz,p_jsize
  	   ia1=min(ii,iz);iz1=min(ii+p_isize-1,iz)
  	   ja1=min(jj,jz);jz1=min(jj+p_jsize-1,jz)
  	 END DO
       END DO
!srf      END DO           

       DO ii=ia,iz
  	 DO jj=ja,jz
	   xland(ii,jj) = patch_area(ii,jj,2)
  	 END DO
       END DO

      
       DO ii=ia,iz,p_isize
  	 DO jj=ja,jz,p_jsize
  	   ia1=min(ii,iz);iz1=min(ii+p_isize-1,iz)
  	   ja1=min(jj,jz);jz1=min(jj+p_jsize-1,jz)

  	   CALL radcarma(m1,m2,m3,ia1,iz1,ja1,jz1,solfac  &
  			    ,theta,pi0,pp,rv,RAIN,LWL,IWL,dn0,rtp,fthrd  &
  			    ,rtgt,f13t,f23t,glat,glon,rshort &
  			    ,rlong,albedt,cosz,rlongup,mynum  &
  			    ,fmapt,pm,carma(ngrid)%aot,xland)
			    
!xxxxxxxxxxxxxx
!	    if (rlong(ii,jj).lt.10.0) print*,'apos radcarma',ii,jj,rlong(ii,jj)
!
!		 do k=2,m1-1
!		   if ( (fthrd(k,ii,jj))*86400 .lt. -20.0)then
!		    print*,'apos radcarma i j k mynum',ii,jj,k,mynum
!		    print*,'apos radcarma rlong fthrd',rlong(ii,jj),fthrd(k,ii,jj)*86400
!		   endif
!  	     end do
!xxxxxxxxxxxxxx
	    		    
  	 END DO
       END DO
       
       
       
  END SUBROUTINE radcomp_carma
  
    SUBROUTINE AllocIndex(ia,ja,iz,jz,IsAlloc)
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: ia,iz,ja,jz,IsAlloc
      INTEGER :: ij,i1,j1
      
      IF(IsAlloc==1) THEN
  	ALLOCATE(indexi((iz-ia+1)*(jz-ja+1)))
  	ALLOCATE(indexj((iz-ia+1)*(jz-ja+1)))
      
  	ij=0
  	DO i1=ia,iz
  	  DO j1=ja,jz
  	    ij=ij+1
  	    indexi(ij)=i1
  	    indexj(ij)=j1
  	  END DO
  	END DO
      ELSE
  	DEALLOCATE(indexi)
  	DEALLOCATE(indexj)
      END IF
    
    END SUBROUTINE AllocIndex
    
    SUBROUTINE C_2d_1d(A2d,A1d,ia,iz,ja,jz,iend)
  
       INTEGER,INTENT(IN) :: ia,iz,ja,jz,iend
       REAL,INTENT(IN) :: A2d(ia:iz,ja:jz)
       REAL,INTENT(OUT) :: A1d(iend)
       INTEGER :: ij
       
       DO ij=1,iend
  	 A1d(ij)=A2d(indexi(ij),indexj(ij))
       END DO
       
    END SUBROUTINE C_2d_1d
    
    SUBROUTINE C_1d_2d(A1d,A2d,ia,iz,ja,jz,iend)
  
       INTEGER,INTENT(IN) :: ia,iz,ja,jz,iend
       REAL,INTENT(IN) :: A1d(iend)
       REAL,INTENT(OUT) :: A2d(ia:iz,ja:jz)
       INTEGER :: ij
       
       DO ij=1,iend
  	A2d(indexi(ij),indexj(ij))=A1d(ij)
       END DO
       
    END SUBROUTINE C_1d_2d
     
    SUBROUTINE C_3d_2d(A3d,A2d,m,ia,iz,ja,jz,iend)
  
       INTEGER,INTENT(IN) :: ia,iz,ja,jz,iend,m
       REAL,INTENT(IN) :: A3d(m,ia:iz,ja:jz)
       REAL,INTENT(OUT) :: A2d(iend,m)
       INTEGER :: ij,l
       
       DO l=1,m
  	 DO ij=1,iend
  	   A2d(ij,l)=A3d(l,indexi(ij),indexj(ij))
  	 END DO
       END DO
       
    END SUBROUTINE C_3d_2d
  
    SUBROUTINE C_2d_3d(A2d,A3d,m,ia,iz,ja,jz,iend)
  
       INTEGER,INTENT(IN) :: ia,iz,ja,jz,iend,m
       REAL,INTENT(OUT)   :: A3d(m,ia:iz,ja:jz)
       REAL,INTENT(IN)    :: A2d(iend,m)
       INTEGER :: ij,l
       
       DO l=1,m
  	 DO ij=1,iend
  	   A3d(l,indexi(ij),indexj(ij))=A2d(ij,l)
  	 END DO
       END DO
       
    END SUBROUTINE C_2d_3d
    
    SUBROUTINE Ci_3d_2d(A3d,A2d,m,ia,iz,ja,jz,iend)
  
       INTEGER,INTENT(IN) :: ia,iz,ja,jz,iend,m
       REAL,INTENT(IN) :: A3d(ia:iz,ja:jz,m)
       REAL,INTENT(OUT) :: A2d(iend,m)
       INTEGER :: ij,l
       
       DO l=1,m
  	 DO ij=1,iend
  	   A2d(ij,l)=A3d(indexi(ij),indexj(ij),l)
  	 END DO
       END DO
       
    END SUBROUTINE Ci_3d_2d

    SUBROUTINE Ci_2d_3d(A2d,A3d,m,ia,iz,ja,jz,iend)
  
       INTEGER,INTENT(IN) :: ia,iz,ja,jz,iend,m
       REAL,INTENT(OUT)   :: A3d(ia:iz,ja:jz,m)
       REAL,INTENT(IN)    :: A2d(iend,m)
       INTEGER :: ij,l
       
       DO l=1,m
  	 DO ij=1,iend
  	   A3d(indexi(ij),indexj(ij),l)=A2d(ij,l)
  	 END DO
       END DO
       
    END SUBROUTINE Ci_2d_3d

  
END MODULE rad_carma
