!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

!-----------------------------------------------------------------------
!     ******  These routines provide line printer and
!             plotted output of model fields and parameters.
!-----------------------------------------------------------------------

SUBROUTINE PRTOUT()
use mem_leaf
use mem_grid
use io_params

real, allocatable:: scrslb(:)
!
!     Print or plot model fields
!
CHARACTER*8 XLAB,YLAB,JNAM,B,IFMT,F,T
CHARACTER*3 TIT
CHARACTER*1 HORIZ,VERT
CHARACTER*68 TT
DIMENSION ABCISA(MAXDIM),ORDIN(MAXDIM),XPT(MAXDIM),YPT(MAXDIM)
EQUIVALENCE (ABCISA,VCTR1),(ORDIN,VCTR2),(XPT,VCTR3),(YPT,VCTR4)

npsize = max(maxnxp,maxnyp,maxnzp)
npsize = npsize ** 2
allocate(scrslb(npsize))

IADDP=0
DO IPLT=1,NPLT

   call sfcprt(nxp,nyp,nzg,nzs,npatch,leaf_g(ngrid),iplfld(iplt),lprt)
   if (lprt == 1) go to 1000

  NAAVG(IPLT)=1
  NOAVG(IPLT)=1
  IAA(IPLT)=0
  IAB(IPLT)=0
  JOA(IPLT)=0
  JOB(IPLT)=0
  call azero(npsize,scrslb)
  I=1
  J=1
  K=0
  AA=OPTLIB(IPLFLD(IPLT),K,I,J,IPLGRD,F,T)
  IF(PLFMT(IPLT)(1:1).EQ.' '.or.PLFMT(IPLT)(1:1).EQ.char(0))  &
      PLFMT(IPLT)=F
  IF(PLTIT(IPLT)(1:1).EQ.' '.or.PLTIT(IPLT)(1:1).EQ.char(0))  &
      PLTIT(IPLT)=T
  FACT=1./FLOAT(NOAVG(IPLT)*NAAVG(IPLT))
  IIAA = 1-IAA(IPLT)
  JJOA = 1-JOA(IPLT)
  JJST=(JJOA-1)/NOAVG(IPLT)+1
  IIST=(IIAA-1)/NAAVG(IPLT)+1
  NAA1=IIAA+NAAVG(IPLT)/2
  NOA1=JJOA+NOAVG(IPLT)/2
!
  IF(IXSCTN(IPLT).EQ.1)THEN
    J=ISBVAL(IPLT)
    IIAB = NXP+IAB(IPLT)
    JJOB = NZP+JOB(IPLT)
    IF(IPLGRD.EQ.5)JJOB= NZG+JOB(IPLT)
    DO I=IIAA,IIAB
      DO K=JJOA,JJOB
        JJJ=(K-1)/NOAVG(IPLT)+1-(JJST-1)
        III=(I-1)/NAAVG(IPLT)+1-(IIST-1)
        DELI=MIN(I+NAAVG(IPLT)-1,NXP)-I+1
        DELJ=MIN(K+NOAVG(IPLT)-1,NZP)-K+1
        FACT=1./(DELI*DELJ)
        indp = (jjj-1) * nxp + iii
        SCRSLB(indp)=SCRSLB(indp)  &
                +OPTLIB(IPLFLD(IPLT),K,I,J,IPLGRD,F,T)*FACT
      ENDDO
    ENDDO
    XLAB=' x(km)'
!         XLAB=' W long'
    YLAB=' z(m)'
    HORIZ='X'
    VERT='Z'
    JNAM='J='
    ihdim = nxp
    ivdim = nzp
  ELSEIF(IXSCTN(IPLT).EQ.2)THEN
    I=ISBVAL(IPLT)
    IIAB = NYP+IAB(IPLT)
    JJOB = NZP+JOB(IPLT)
    IF(IPLGRD.EQ.5)JJOB= NZG+JOB(IPLT)
    DO J=IIAA,IIAB
      DO K=JJOA,JJOB
        JJJ=(K-1)/NOAVG(IPLT)+1-(JJST-1)
        III=(J-1)/NAAVG(IPLT)+1-(IIST-1)
        DELI=MIN(J+NAAVG(IPLT)-1,NYP)-J+1
        DELJ=MIN(K+NOAVG(IPLT)-1,NZP)-K+1
        FACT=1./(DELI*DELJ)
        indp = (jjj-1) * nyp + iii
        SCRSLB(indp)=SCRSLB(indp)  &
                +OPTLIB(IPLFLD(IPLT),K,I,J,IPLGRD,F,T)*FACT
      ENDDO
    ENDDO
    XLAB=' y(km)'
!         XLAB=' N lat '
    YLAB=' z(m)'
    HORIZ='Y'
    VERT='Z'
    JNAM='I='
    ihdim = nyp
    ivdim = nzp
  ELSEIF(IXSCTN(IPLT).EQ.3)THEN
    K=ISBVAL(IPLT)
    IIAB = NXP+IAB(IPLT)
    JJOB = NYP+JOB(IPLT)
    DO I=IIAA,IIAB
      DO J=JJOA,JJOB
        JJJ=(J-1)/NOAVG(IPLT)+1-(JJST-1)
        III=(I-1)/NAAVG(IPLT)+1-(IIST-1)
        DELI=MIN(I+NAAVG(IPLT)-1,NXP)-I+1
        DELJ=MIN(J+NOAVG(IPLT)-1,NYP)-J+1
        FACT=1./(DELI*DELJ)
        indp = (jjj-1) * nxp + iii
        SCRSLB(indp)=SCRSLB(indp)  &
                +OPTLIB(IPLFLD(IPLT),K,I,J,IPLGRD,F,T)*FACT
      ENDDO
    ENDDO
    XLAB=' x(km)'
!         XLAB=' W long'
    YLAB=' y(km)'
!         YLAB='  N lat '
    HORIZ='X'
    VERT='Y'
    JNAM='K='
    ihdim = nxp
    ivdim = nyp
  ENDIF
!
  LL=(IIAB-IIAA)/NAAVG(IPLT)+1
  MM=(JJOB-JJOA)/NOAVG(IPLT)+1
!
  DO K=1,LL
    IF(IXSCTN(IPLT).EQ.1.OR.IXSCTN(IPLT).EQ.3)THEN
        ABCISA(K)=XT(MIN((K-1)*NAAVG(IPLT)+NAA1,NXP))*1E-3
        IF(IPLGRD.EQ.1)  &
           ABCISA(K)=XM(MIN((K-1)*NAAVG(IPLT)+NAA1,NXP))*1E-3
    ELSE
        ABCISA(K)=YT(MIN((K-1)*NAAVG(IPLT)+NAA1,NYP))*1E-3
        IF(IPLGRD.EQ.2)  &
           ABCISA(K)=YM(MIN((K-1)*NAAVG(IPLT)+NAA1,NYP))*1E-3
    ENDIF
  ENDDO

  DO K=1,MM
    IF(IXSCTN(IPLT).LE.2)THEN
        ORDIN(K)=ZT(MIN((K-1)*NOAVG(IPLT)+NOA1,NZP))
      IF(IPLGRD.EQ.3)  &
        ORDIN(K)=ZM(MIN((K-1)*NOAVG(IPLT)+NOA1,NZP))
      IF(IPLGRD.EQ.5)  &
        ORDIN(K)=SLZ(MIN((K-1)*NOAVG(IPLT)+NOA1,NZG))
    ELSE
        ORDIN(K)=YT(MIN((K-1)*NOAVG(IPLT)+NOA1,NYP))*1E-3
      IF(IPLGRD.EQ.2)  &
        ORDIN(K)=YM(MIN((K-1)*NOAVG(IPLT)+NOA1,NYP))*1E-3
    ENDIF
  ENDDO

  HR=real(TIME/3600)

  WRITE(TT,41) NGRID,PLTIT(IPLT),TIME,HR,JNAM,ISBVAL(IPLT)
41      FORMAT('Grid:',I3,'    Field: ',A8,'   Time: ',F8.1,'s /'  &
  ,F6.2,'h   Slab: ',A2,I3)

 CALL PRT2D(HORIZ,VERT,SCRSLB,ABCISA,ORDIN,ihdim,ivdim  &
              ,1,LL,1,MM,PLFMT(IPLT),TT,XPT,YPT,XLAB,YLAB)

1000 CONTINUE
enddo

deallocate(scrslb)

RETURN
END

!     ******************************************************************

subroutine sfcprt(n2,n3,mzg,mzs,npat,leaf,vnam,lprt)

use mem_leaf
use leaf_coms
use rconstants, only: wdns
use therm_lib, only: qtk, qwtk
type (leaf_vars) :: leaf

dimension tempkk(20),fracliqq(20),area(20)
character*(*) vnam
character*26 vnam2
character*132 line

lprt = 1
maxcols = 15
npages = (n2 + maxcols - 1) / maxcols
leftcols = n2 - (npages - 1) * maxcols

k1 = 1
k2 = 1

    if (vnam == 'soil_water'     ) then
        vnam2 = '[soil_water (m3/m3) ]    '
        k2 = nzg
elseif (vnam == 'soil_energy'    ) then
        vnam2 = '[soil_energy (J/cm3)]    '
        k2 = nzg
elseif (vnam == 'soil_temp'      ) then
        vnam2 = '[soil_temp (K)]          '
        k2 = nzg
elseif (vnam == 'soil_text'      ) then
        vnam2 = '[soil text]              '
        k2 = nzg
elseif (vnam == 'sfcwater_mass'  ) then
        vnam2 = '[sfcwater_mass (kg/m2)]  '
        k2 = nzs
elseif (vnam == 'sfcwater_energy') then
        vnam2 = '[sfcwater_energy (J/g)]  '
        k2 = nzs
elseif (vnam == 'sfcwater_temp'  ) then
        vnam2 = '[sfcwater_temp (K)]      '
        k2 = nzs
elseif (vnam == 'sfcwater_depth' ) then
        vnam2 = '[sfcwater_depth (m)]     '
        k2 = nzs
elseif (vnam == 'ustar'          ) then
        vnam2 = '[ustar (m/s)]            '
elseif (vnam == 'tstar'          ) then
        vnam2 = '[tstar (K)]              '
elseif (vnam == 'rstar'          ) then
        vnam2 = '[rstar (g/kg)]           '
elseif (vnam == 'cstar'          ) then
        vnam2 = '[cstar (umol/mol)]       '
elseif (vnam == 'zeta'           ) then
        vnam2 = '[zeta]                   '
elseif (vnam == 'ribulk'         ) then
        vnam2 = '[ribulk]                 '
elseif (vnam == 'veg_fracarea'   ) then
        vnam2 = '[veg_fracarea]           '
elseif (vnam == 'veg_agb'        ) then
        vnam2 = '[veg_agb]                '
elseif (vnam == 'veg_lai'        ) then
        vnam2 = '[veg_lai]                '
elseif (vnam == 'veg_tai'        ) then
        vnam2 = '[veg_tai]                '
elseif (vnam == 'veg_rough'      ) then
        vnam2 = '[veg rough (m)]          '
elseif (vnam == 'veg_albedo'     ) then
        vnam2 = '[veg_albedo (m)]         '
elseif (vnam == 'veg_height'     ) then
        vnam2 = '[veg_height (m)]         '
elseif (vnam == 'patch_area'     ) then
        vnam2 = '[patch_area]             '
elseif (vnam == 'patch_rough'    ) then
        vnam2 = '[patch rough (m)]        '
elseif (vnam == 'patch_wetind'   ) then
        vnam2 = '[patch wetind]           '
elseif (vnam == 'leaf_class'     ) then
        vnam2 = '[leaf_class]             '
elseif (vnam == 'soil_rough'     ) then
        vnam2 = '[soil rough (m)]         '
elseif (vnam == 'sfcwater_nlev'  ) then
        vnam2 = '[sfcwater_nlev]          '
elseif (vnam == 'stom_resist'    ) then
        vnam2 = '[stom_resist (s/m)]      '
elseif (vnam == 'ground_rsat'    ) then
        vnam2 = '[ground_rsat (g/kg)]     '
elseif (vnam == 'ground_rvap'    ) then
        vnam2 = '[ground_rvap (g/kg)]     '
elseif (vnam == 'ground_temp'    ) then
        vnam2 = '[ground_temp    (K)]     '
elseif (vnam == 'ground_fliq'    ) then
        vnam2 = '[ground_fliq    (%)]     '
elseif (vnam == 'veg_water'      ) then
        vnam2 = '[veg_water (kg/m2)]      '
elseif (vnam == 'veg_hcap'       ) then
        vnam2 = '[veg hcap (J/m2/K)]      '
elseif (vnam == 'veg_energy'     ) then
        vnam2 = '[veg energy (J/m2)]      '
elseif (vnam == 'can_prss'       ) then
        vnam2 = '[can_prss   (Pa)]        '
elseif (vnam == 'can_theta'      ) then
        vnam2 = '[can_theta   (K)]        '
elseif (vnam == 'can_rvap'       ) then
        vnam2 = '[can_rvap (g/kg)]        '
elseif (vnam == 'can_co2'       ) then
        vnam2 = '[can_co2  (umol/mol)]    '
elseif (vnam == 'veg_ndvip'      ) then
        vnam2 = '[veg_ndvip]              '
elseif (vnam == 'veg_ndvic'      ) then
        vnam2 = '[veg_ndvic]              '
elseif (vnam == 'veg_ndvif'      ) then
        vnam2 = '[veg_ndvif]              '
else
   lprt = 0
   return
endif

do k = k1,k2

   do ipage = 1,npages
      nc = 15
      if (ipage == npages) nc = leftcols
      i1 = (ipage - 1) * 15 + 1
      i2 = i1 + nc - 1

      write(6,20)
      write(6,21) vnam2,k,ngrid,time,ipage,npages
      write(6,22) (i,i=i1,i2)

      do j = n3,1,-1
         do ipat = 1,npat
         
            if     (vnam == 'soil_water'     ) then
               call plin(nc,j,ipat,3,1.               &
                  ,leaf%soil_water(k,i1:i2,j,ipat)       &
                  ,leaf%patch_area(i1:i2,j,ipat))
            elseif (vnam == 'soil_energy'    ) then
               call plin(nc,j,ipat,2,1.e-6            &
                  ,leaf%soil_energy(k,i1:i2,j,ipat)      &
                  ,leaf%patch_area(i1:i2,j,ipat))
            elseif (vnam == 'soil_temp'      ) then

               do i = i1,i2
                  if (ipat == 1 .and. k == mzg) then
                     call qtk(leaf%soil_energy(k,i,j,ipat)  &
                        ,tempkk(i+1-i1),fracliqq(i+1-i1))
                  elseif (ipat == 1) then
                     tempk(i+1-i1) = leaf%soil_energy(k,i,j,ipat)
                  else
                     nsoil = nint(leaf%soil_text(k,i,j,ipat))
                     call qwtk(leaf%soil_energy(k,i,j,ipat)       &
                              ,leaf%soil_water (k,i,j,ipat)*wdns  &
                              ,slcpd(nsoil),tempkk(i+1-i1),fracliqq(i+1-i1))
                  endif
               enddo

               call plin(nc,j,ipat,2,1.               &
                  ,tempk(1:nc)                        &
                  ,leaf%patch_area(i1:i2,j,ipat))
            elseif (vnam == 'soil_text'      ) then
               call plin(nc,j,ipat,3,1.                  &
                  ,leaf%soil_text(k,i1:i2,j,ipat)        &
                  ,leaf%patch_area(i1:i2,j,ipat))
            elseif (vnam == 'sfcwater_mass'  ) then
               call plin(nc,j,ipat,3,1.                  &
                  ,leaf%sfcwater_mass(k,i1:i2,j,ipat)    &
                  ,leaf%patch_area(i1:i2,j,ipat))
            elseif (vnam == 'sfcwater_energy') then
               call plin(nc,j,ipat,2,1.e-3               &
                  ,leaf%sfcwater_energy(k,i1:i2,j,ipat)  &
                  ,leaf%patch_area(i1:i2,j,ipat))
            elseif (vnam == 'sfcwater_temp'  ) then
               do i = i1,i2
                  nsoil = nint(leaf%soil_text(k,i,j,ipat))
                  call qwtk(leaf%soil_energy(k,i,j,ipat)       &
                           ,leaf%soil_water (k,i,j,ipat)*wdns  &
                           ,slcpd(nsoil),tempkk(i+1-i1),fracliqq(i+1-i1))
               enddo
               call plin(nc,j,ipat,2,1.               &
                  ,tempk(1:nc)                        &
                  ,leaf%patch_area(i1:i2,j,ipat))
            elseif (vnam == 'sfcwater_depth' ) then
               call plin(nc,j,ipat,3,1.               &
                  ,leaf%sfcwater_depth(k,i1:i2,j,ipat)   &
                  ,leaf%patch_area(i1:i2,j,ipat))
            elseif (vnam == 'ustar'          ) then
               call plin(nc,j,ipat,3,1.               &
                  ,leaf%ustar(i1:i2,j,ipat)              &
                  ,leaf%patch_area(i1:i2,j,ipat))
            elseif (vnam == 'tstar'          ) then
               call plin(nc,j,ipat,3,1.               &
                  ,leaf%tstar(i1:i2,j,ipat)              &
                  ,leaf%patch_area(i1:i2,j,ipat))
            elseif (vnam == 'rstar'          ) then
               call plin(nc,j,ipat,3,1.e3             &
                  ,leaf%rstar(i1:i2,j,ipat)              &
                  ,leaf%patch_area(i1:i2,j,ipat))
            elseif (vnam == 'cstar'          ) then
               call plin(nc,j,ipat,3,1.                  &
                  ,leaf%cstar(i1:i2,j,ipat)              &
                  ,leaf%patch_area(i1:i2,j,ipat))
            elseif (vnam == 'zeta'           ) then
               call plin(nc,j,ipat,3,1.                  &
                  ,leaf%zeta(i1:i2,j,ipat)               &
                  ,leaf%patch_area(i1:i2,j,ipat))
            elseif (vnam == 'ribulk'         ) then
               call plin(nc,j,ipat,3,1.                  &
                  ,leaf%ribulk(i1:i2,j,ipat)             &
                  ,leaf%patch_area(i1:i2,j,ipat))
            elseif (vnam == 'veg_fracarea'   ) then
               call plin(nc,j,ipat,3,1.               &
                  ,leaf%veg_fracarea(i1:i2,j,ipat)       &
                  ,leaf%patch_area(i1:i2,j,ipat))
            elseif (vnam == 'veg_agb'        ) then
               call plin(nc,j,ipat,3,1.               &
                  ,leaf%veg_agb(i1:i2,j,ipat)            &
                  ,leaf%patch_area(i1:i2,j,ipat))
            elseif (vnam == 'veg_lai'        ) then
               call plin(nc,j,ipat,3,1.               &
                  ,leaf%veg_lai(i1:i2,j,ipat)            &
                  ,leaf%patch_area(i1:i2,j,ipat))
            elseif (vnam == 'veg_tai'        ) then
               call plin(nc,j,ipat,3,1.               &
                  ,leaf%veg_tai(i1:i2,j,ipat)            &
                  ,leaf%patch_area(i1:i2,j,ipat))
            elseif (vnam == 'veg_rough'      ) then
               call plin(nc,j,ipat,3,1.               &
                  ,leaf%veg_rough(i1:i2,j,ipat)          & 
                  ,leaf%patch_area(i1:i2,j,ipat))
            elseif (vnam == 'veg_albedo'     ) then
               call plin(nc,j,ipat,3,1.               &
                  ,leaf%veg_albedo(i1:i2,j,ipat)         &
                  ,leaf%patch_area(i1:i2,j,ipat))
            elseif (vnam == 'veg_height'     ) then
               call plin(nc,j,ipat,3,1.               &
                  ,leaf%veg_height(i1:i2,j,ipat)         &
                  ,leaf%patch_area(i1:i2,j,ipat))
            elseif (vnam == 'patch_area'     ) then
               write(6,243)j,ipat,(leaf%patch_area(i,j,ipat),i=i1,i2)
            elseif (vnam == 'patch_rough'    ) then
               call plin(nc,j,ipat,3,1.               &
                  ,leaf%patch_rough(i1:i2,j,ipat)        &
                  ,leaf%patch_area(i1:i2,j,ipat))
            elseif (vnam == 'patch_wetind'   ) then
               call plin(nc,j,ipat,3,1.               &
                  ,leaf%patch_wetind(i1:i2,j,ipat)       &
                  ,leaf%patch_area(i1:i2,j,ipat))
            elseif (vnam == 'leaf_class'    ) then
               call plin(nc,j,ipat,3,1.               &
                  ,leaf%leaf_class(i1:i2,j,ipat)         &
                  ,leaf%patch_area(i1:i2,j,ipat))
            elseif (vnam == 'soil_rough'     ) then
               call plin(nc,j,ipat,3,1.               &
                  ,leaf%soil_rough(i1:i2,j,ipat)         &
                  ,leaf%patch_area(i1:i2,j,ipat))
            elseif (vnam == 'sfcwater_nlev'  ) then
               call plin(nc,j,ipat,3,1.               &
                  ,leaf%sfcwater_nlev(i1:i2,j,ipat)      &
                  ,leaf%patch_area(i1:i2,j,ipat))
            elseif (vnam == 'stom_resist'    ) then
               call plin(nc,j,ipat,3,1.               &
                  ,leaf%stom_resist(i1:i2,j,ipat)        &
                  ,leaf%patch_area(i1:i2,j,ipat))
            elseif (vnam == 'ground_rsat'    ) then
               call plin(nc,j,ipat,3,1.e3             &
                  ,leaf%ground_rsat(i1:i2,j,ipat)        &
                  ,leaf%patch_area(i1:i2,j,ipat))
            elseif (vnam == 'ground_rvap'    ) then
               call plin(nc,j,ipat,3,1.e3             &
                  ,leaf%ground_rvap(i1:i2,j,ipat)        &
                  ,leaf%patch_area(i1:i2,j,ipat))
            elseif (vnam == 'ground_temp'    ) then
               call plin(nc,j,ipat,3,1.               &
                  ,leaf%ground_temp(i1:i2,j,ipat)        &
                  ,leaf%patch_area(i1:i2,j,ipat))
            elseif (vnam == 'ground_fliq'    ) then
               call plin(nc,j,ipat,3,1.e2             &
                  ,leaf%ground_fliq(i1:i2,j,ipat)        &
                  ,leaf%patch_area(i1:i2,j,ipat))
            elseif (vnam == 'veg_water'      ) then
               call plin(nc,j,ipat,3,1.               &
                  ,leaf%veg_water(i1:i2,j,ipat)          &
                  ,leaf%patch_area(i1:i2,j,ipat))
            elseif (vnam == 'veg_hcap'       ) then
               call plin(nc,j,ipat,2,1.               &
                  ,leaf%veg_hcap(i1:i2,j,ipat)           &
                  ,leaf%patch_area(i1:i2,j,ipat))
            elseif (vnam == 'veg_energy'     ) then
               call plin(nc,j,ipat,2,1.               &
                  ,leaf%veg_energy(i1:i2,j,ipat)         &
                  ,leaf%patch_area(i1:i2,j,ipat))
            elseif (vnam == 'can_prss'       ) then
               call plin(nc,j,ipat,3,1.e2             &
                  ,leaf%can_prss(i1:i2,j,ipat)           &
                  ,leaf%patch_area(i1:i2,j,ipat))
            elseif (vnam == 'can_theta'      ) then
               call plin(nc,j,ipat,3,1.e0             &
                  ,leaf%can_theta(i1:i2,j,ipat)          &
                  ,leaf%patch_area(i1:i2,j,ipat))
            elseif (vnam == 'can_co2'       ) then
               call plin(nc,j,ipat,3,1.                  &
                  ,leaf%can_co2(i1:i2,j,ipat)            &
                  ,leaf%patch_area(i1:i2,j,ipat))
            elseif (vnam == 'veg_ndvip'      ) then
               call plin(nc,j,ipat,3,1.               &
                  ,leaf%veg_ndvip(i1:i2,j,ipat)          &
                  ,leaf%patch_area(i1:i2,j,ipat))
            elseif (vnam == 'veg_ndvic'      ) then
               call plin(nc,j,ipat,3,1.               &
                  ,leaf%veg_ndvic(i1:i2,j,ipat)          &
                  ,leaf%patch_area(i1:i2,j,ipat))
            elseif (vnam == 'veg_ndvif'      ) then
               call plin(nc,j,ipat,3,1.               &
                  ,leaf%veg_ndvif(i1:i2,j,ipat)          &
                  ,leaf%patch_area(i1:i2,j,ipat))
            endif

         enddo
         write(6,19)
      enddo
   enddo
   write(6,19)
enddo

19   format(' ')
20   format(6('_______________________'))
21   format(a26,'  k =',i3,'  [grid ',i2,']  [time ',f8.0  &
         ,']  [pg',i3,'/',i2,']')
22   format('_________i=__',20('_____',i3))
242  format('j=',i3,'  pat=',i2,1x,20f8.2)
243  format('j=',i3,'  pat=',i2,1x,20f8.3)
return
end

!     ****************************************************************

subroutine plin(nc,j,ipat,ifmt,factor,field,area)
character*120 line
real, dimension(nc) :: field,area
real :: factor

! Initialize line to all blanks

do i = 1,15
   j1 = 8 * (i-1) + 1
   j2 = j1 + 7
   line(j1:j2) = '        '
enddo

! Fill line with values

if (ifmt .eq. 1) then
   WRITE(line,'(20f8.1)') (field(i)*factor,i=1,nc)
elseif (ifmt .eq. 2) then
   WRITE(line,'(20f8.2)') (field(i)*factor,i=1,nc)
elseif (ifmt .eq. 3) then
   WRITE(line,'(20f8.3)') (field(i)*factor,i=1,nc)
endif

! Blank out values in patches of zero area

do i = 1,nc
   if (area(i) .lt. .009) then
      j1 = 8 * (i-1) + 1
      j2 = j1 + 7
      line(j1:j2) = '        '
   endif
enddo

write(6,341)j,ipat,line
341  format('j=',i3,'  pat=',i2,1x,a120)

return
end

!     ****************************************************************

FUNCTION OPTLIB(VARN,K,I,J,IPLGRD,FMT,TILO)

use mem_all
use var_tables
use rconstants
use therm_lib, only: rehul,rehui,rehuil,virtt
CHARACTER*(*) VARN
CHARACTER*8 FMT,TILO
!
!     'UC'   - Current U wind component (M/S)
!     'VC'   - Current V wind component (M/S)
!     'WC'   - Current W wind component (CM/S)
!     'PC'   - Current pressure (MB)
!     'UP'   - Past U wind component (M/S)
!     'VP'   - Past V wind component (M/S)
!     'WP'   - Past W wind component (CM/S)
!     'PP'   - Past pressure (MB)

!     'THP' - Theta-il (K)
!     'THETA'- THETA(K)
!     'TH0' - Reference state theta_v
!     'DN0' - Reference state density
!     'PI0' - Reference state Exner function
!     'THV'  - Theta-v     (K)
!     'THVP' - Perturbation Theta_v (K)
!     'TEMP' - temperature (K)
!     'TV'   - Virtual temperature (K)
!     'TVP'  - Perturbation virtual temperature (K)

!     'RT'     - Total water mixing ratio (G/KG)
!     'RV'     - Vapor mixing ratio (G/KG)
!     'RC'     - Cloud water mixing ratio (G/KG)
!     'RR'     - Rain mixing ratio (G/KG)
!     'RP'     - Pristine ice mixing ratio (G/KG)
!     'RS'     - Snow mixing ratio (G/KG)
!     'RA'     - Aggregates mixing ratio (G/KG)
!     'RG'     - Graupel mixing ratio (G/KG)
!     'RH'     - Hail mixing ratio (G/KG)
!     'Q2'     - Rain internal energy (cal/G)
!     'Q6'     - Graupel internal energy (cal/G)
!     'Q7'     - Hail internal energy (cal/G)
!     'RL'     - Liquid water mixing ratio (G/KG)
!     'RI'     - Ice mixing ratio (G/KG)
!     'RCOND'  - Total condensate mixing ratio (G/KG)

!     'CC'     - Cloud droplet number concentration (#/KG)
!     'CR'     - Rain number concentration (#/KG)
!     'CP'     - Pristine ice number concentration (#/KG)
!     'CS'     - Snow number concentration (#/KG)
!     'CA'     - Aggregates number concentration (#/KG)
!     'CG'     - Graupel number concentration (#/KG)
!     'CH'     - Hail number concentration (#/KG)

!     'RTP'    - Perturbation total water mixing ratio(G/KG)
!     'RELHUM' - relative humidity (%)
!     'PCPT'   - Total accumulated precipitation (KG/M^2)
!     'TKE'    - Turbulent Kinetic Energy (J/KG)
!     'RSHORT' - Surface downward shortwave radiation (W/M^2)
!     'RLONG'  - Surface downward longwave radiation (W/M^2)
!     'GLAT'   - Latitude (degrees)
!     'GLON'   - Longitude (degrees)
!     'CONPR'  - Convective precipitation rate (KG/M^2/S)
!     'CONP'   - Accumulated convective precipitation (KG/M^2)
!     'CONH'   - Convective precipitation heating rate (K/day)
!     'CONM'   - Convective moistening rate (KG/KG/day)

!     'MICRO'  - GASPRC
!     'SPEED'  - Wind speed (M/S)
!     'BVF'    - Brunt-Vaisala Frequency (100*1/S*S)
!     'RICH'   - Richardson number
!     'VKM'    - Vertical K for momentum (M*M/S)
!     'VKH'    - Vertical K for heat (M*M/S)
!     'HKM'    - Horizontal K for momentum (M*M/S)
!     'HKH'    - Horizontal K for heat (M*M/S)
!     'FTHRD'  - Radiative flux convergence (K/S)
!     'ZI'     - Boundary layer height above ground (M)
!     'SFLUX_U'- Surface momentum flux for U component  (Pa)
!     'SFLUX_V'- Surface momentum flux for V component  (Pa)
!     'SFLUX_W'- Surface momentum flux for W component  (Pa)
!     'SFLUX_T'- Surface heat flux (K*kg/(M*M*S))
!     'SFLUX_R'- Surface moisture flux (KG/(M*M*S))
!
!         Current U COMPONENT
IF(VARN.EQ.'UC') THEN
  IF(K.EQ.0)THEN
    IPLGRD=1
    FMT='0PF7.1'
    TILO='UC(M/S)'
  ELSE
    OPTLIB=basic_g(ngrid)%UC(k,i,j)
  ENDIF
!
!         Current V COMPONENT
ELSEIF(VARN.EQ.'VC') THEN
  IF(K.EQ.0)THEN
    IPLGRD=2
    FMT='0PF7.1'
    TILO='VC(M/S)'
  ELSE
    OPTLIB=basic_g(ngrid)%VC(k,i,j)
  ENDIF
!
!         Current W COMPONENT
ELSEIF(VARN.EQ.'WC') THEN
  IF(K.EQ.0)THEN
    IPLGRD=3
    FMT=' 2PF7.1'
    TILO='WC(CM/S)'
  ELSE
    OPTLIB=basic_g(ngrid)%WC(k,i,j)
  ENDIF
!
!         Current PERTURBATION PRESSURE
ELSEIF(VARN.EQ.'PC') THEN
  IF(K.EQ.0)THEN
    IPLGRD=4
    FMT='-2PF7.1'
    TILO='PC (MB)'
  ELSE
    VAL1=basic_g(ngrid)%PC(k,i,j)
    VAL2=basic_g(ngrid)%DN0(k,i,j)
    VAL3=basic_g(ngrid)%TH0(k,i,j)
!cc          OPTLIB=VAL1*VAL2*VAL3
    OPTLIB=VAL1
  ENDIF
!
!         Past U COMPONENT
ELSEIF(VARN.EQ.'UP') THEN
  IF(K.EQ.0)THEN
    IPLGRD=1
    FMT='0PF7.1'
    TILO='UP(M/S)'
  ELSE
    OPTLIB=basic_g(ngrid)%UP(k,i,j)
  ENDIF
!
!         Past V COMPONENT
ELSEIF(VARN.EQ.'VP') THEN
  IF(K.EQ.0)THEN
    IPLGRD=2
    FMT='0PF7.1'
    TILO='VP(M/S)'
  ELSE
    OPTLIB=basic_g(ngrid)%VP(k,i,j)
  ENDIF
!
!         Past W COMPONENT
ELSEIF(VARN.EQ.'WP') THEN
  IF(K.EQ.0)THEN
    IPLGRD=3
    FMT=' 2PF7.1'
    TILO='WP(CM/S)'
  ELSE
    OPTLIB=basic_g(ngrid)%WP(k,i,j)
  ENDIF
!
!         Past PERTURBATION PRESSURE
ELSEIF(VARN.EQ.'PP') THEN
  IF(K.EQ.0)THEN
    IPLGRD=4
    FMT='-2PF7.1'
    TILO='PP (MB)'
  ELSE
    VAL1=basic_g(ngrid)%PP(k,i,j)
    VAL2=basic_g(ngrid)%DN0(k,i,j)
    VAL3=basic_g(ngrid)%TH0(k,i,j)
    OPTLIB=VAL1*VAL2*VAL3
  ENDIF
!
!           Theta-il
ELSEIF(VARN.EQ.'THP') THEN
  IF(K.EQ.0)THEN
    IPLGRD=4
    FMT=' 0PF7.1'
    TILO='THP'
  ELSE
    OPTLIB=basic_g(ngrid)%thp(k,i,j)
  ENDIF
!
!           THETA
ELSEIF(VARN.EQ.'THETA') THEN
  IF(K.EQ.0)THEN
    IPLGRD=4
    FMT=' 0PF7.1'
    TILO='THETA'
  ELSE
    OPTLIB=basic_g(ngrid)%THETA(k,i,j)
  ENDIF
!
!           Reference state theta_v
ELSEIF(VARN.EQ.'TH0') THEN
  IF(K.EQ.0)THEN
    IPLGRD=4
    FMT=' 0PF7.2'
    TILO='TH0'
  ELSE
    OPTLIB=basic_g(ngrid)%TH0(k,i,j)
  ENDIF
!
!           Reference state density
ELSEIF(VARN.EQ.'DN0') THEN
  IF(K.EQ.0)THEN
    IPLGRD=4
    FMT=' 0PF7.4'
    TILO='DN0'
  ELSE
    OPTLIB=basic_g(ngrid)%DN0(k,i,j)
  ENDIF
!
!           PI0
ELSEIF(VARN.EQ.'PI0') THEN
  IF(K.EQ.0)THEN
    IPLGRD=4
    FMT=' 0PF7.1'
    TILO='PI0'
  ELSE
    OPTLIB=basic_g(ngrid)%PI0(k,i,j)
  ENDIF
!
!           Virtual Potential Temperature
ELSEIF(VARN.EQ.'THV') THEN
  IF(K.EQ.0)THEN
    IPLGRD=4
    FMT=' 0PF7.1'
    TILO='THETA-V'
  ELSE
    VALTHET = basic_g(ngrid)%THETA(k,i,j)
    VALRTP  = basic_g(ngrid)%RTP(k,i,j)
    VALRV   = basic_g(ngrid)%RV(k,i,j)
    OPTLIB  = virtt(VALTHET,VALRV,VALRTP)
  ENDIF
!
!           PERTURBATION THETA V
ELSEIF(VARN.EQ.'THVP') THEN
  IF(K.EQ.0)THEN
    IPLGRD=4
    FMT='0PF7.1'
    TILO='THV''(K)'
  ELSE
    VALTHET = basic_g(ngrid)%THETA(k,i,j)
    VALTH0  = basic_g(ngrid)%TH0(k,i,j)
    VALRTP  = basic_g(ngrid)%RTP(k,i,j)
    VALRV   = basic_g(ngrid)%RV(k,i,j)
    IF(LEVEL.EQ.0) THEN
       VALRTP = 0.
       VALRV  = 0.
    ENDIF
    OPTLIB = virtt(VALTHET,VALRV,VALRTP)-VALTH0
  ENDIF
!
!            TEMPERATURE
ELSEIF(VARN.EQ.'TEMP') THEN
  IF(K.EQ.0)THEN
    IPLGRD=4
    FMT='0PF7.2'
    TILO='TEMP (K)'
  ELSE
    VALTHET = basic_g(ngrid)%THETA(k,i,j)
    VALPP   = basic_g(ngrid)%PP(k,i,j)
    VALPI0  = basic_g(ngrid)%PI0(k,i,j)
    OPTLIB  = VALTHET*(VALPP+VALPI0)/CP
  ENDIF
!
!            TEMPERATURE
ELSEIF(VARN.EQ.'TEMPC') THEN
  IF(K.EQ.0)THEN
    IPLGRD=4
    FMT='0PF7.2'
    TILO='TEMP (C)'
  ELSE
    VALTHET = basic_g(ngrid)%THETA(k,i,j)
    VALPP   = basic_g(ngrid)%PP(k,i,j)
    VALPI0  = basic_g(ngrid)%PI0(k,i,j)
    OPTLIB  = VALTHET*(VALPP+VALPI0)/CP - t00
  ENDIF
!
!           VIRTUAL TEMPERATURE
ELSEIF(VARN.EQ.'TV') THEN
  IF(K.EQ.0)THEN
    IPLGRD=4
    FMT='0PF7.2'
    TILO='TV (K)'
  ELSE
    VALPI0 = basic_g(ngrid)%PI0(k,i,j)
    VALTHET= basic_g(ngrid)%THETA(k,i,j)
    VALPP  = basic_g(ngrid)%PP(k,i,j)
    VALRTP = basic_g(ngrid)%RTP(k,i,j)
    VALRV = 0
    IF (LEVEL.GE.1)  VALRV = basic_g(ngrid)%RV(k,i,j)
    OPTLIB = virtt(VALTHET,VALRV,VALRTP)  &
                 *(VALPP + VALPI0)/CP
  ENDIF
!
!           VIRTUAL TEMPERATURE PERTURBATION
ELSEIF(VARN.EQ.'TVP') THEN
  IF(K.EQ.0)THEN
    IPLGRD=4
    FMT='0PF7.2'
    TILO='TV'' (K)'
  ELSE
    VALPI0 = basic_g(ngrid)%PI0(k,i,j)
    VALTHV0= basic_g(ngrid)%TH0(k,i,j)
    VALTHET= basic_g(ngrid)%THETA(k,i,j)
    VALPP  = basic_g(ngrid)%PP(k,i,j)
    VALRTP = basic_g(ngrid)%RTP(k,i,j)
    VALRV = 0
    IF (LEVEL.GE.1)  VALRV = basic_g(ngrid)%RV(k,i,j)
    OPTLIB = (virtt(VALTHET,VALRV,VALRTP) -VALTHV0)  &
                 *(VALPP + VALPI0)/CP
  ENDIF
!
!           TOTAL WATER MIXING RATIO
ELSEIF(VARN.EQ.'RT') THEN
  IF(K.EQ.0)THEN
    IPLGRD=4
    FMT='3PF7.3'
    TILO='RT(G/KG)'
  ELSE
    OPTLIB=basic_g(ngrid)%rtp(k,i,j)
  ENDIF
!
!           WATER  VAPOR MIXING RATIO
ELSEIF(VARN.EQ.'RV') THEN
  IF(K.EQ.0)THEN
    IPLGRD=4
    FMT='3PF7.3'
    TILO='RV(G/KG)'
  ELSE
    OPTLIB=basic_g(ngrid)%RV(k,i,j)
  ENDIF
!
!           CLOUD WATER MIXING RATIO
ELSEIF(VARN.EQ.'RC') THEN
  IF(K.EQ.0)THEN
    IPLGRD=4
    FMT='3PF7.3'
    TILO='RC(G/KG)'
  ELSE
    OPTLIB=micro_g(ngrid)%RCP(k,i,j)
  ENDIF
!
!           RAIN MIXING RATIO
ELSEIF(VARN.EQ.'RR') THEN
  IF(K.EQ.0)THEN
    IPLGRD=4
    FMT='3PF7.3'
    TILO='RR(G/KG)'
  ELSE
    OPTLIB=micro_g(ngrid)%RRP(k,i,j)
  ENDIF
!
!           PRISTINE CRYSTAL MIXING RATIO
ELSEIF(VARN.EQ.'RP') THEN
  IF(K.EQ.0)THEN
    IPLGRD=4
    FMT='3PF7.3'
    TILO='RP(G/KG)'
  ELSE
    OPTLIB=micro_g(ngrid)%RPP(k,i,j)
  ENDIF
!
!           SNOW MIXING RATIO
ELSEIF(VARN.EQ.'RS') THEN
  IF(K.EQ.0)THEN
    IPLGRD=4
    FMT='3PF7.3'
    TILO='RS(G/KG)'
  ELSE
    OPTLIB=micro_g(ngrid)%RSP(k,i,j)
  ENDIF
!
!           AGGREGATE MIXING RATIO
ELSEIF(VARN.EQ.'RA') THEN
  IF(K.EQ.0)THEN
    IPLGRD=4
    FMT='3PF7.3'
    TILO='RA(G/KG)'
  ELSE
    OPTLIB=micro_g(ngrid)%RAP(k,i,j)
  ENDIF
!
!           GRAUPEL MIXING RATIO
ELSEIF(VARN.EQ.'RG') THEN
  IF(K.EQ.0)THEN
    IPLGRD=4
    FMT='3PF7.3'
    TILO='RG(G/KG)'
  ELSE
    OPTLIB=micro_g(ngrid)%RGP(k,i,j)
  ENDIF
!
!           HAIL MIXING RATIO
ELSEIF(VARN.EQ.'RH') THEN
  IF(K.EQ.0)THEN
    IPLGRD=4
    FMT='3PF7.3'
    TILO='RH(G/KG)'
  ELSE
    OPTLIB=micro_g(ngrid)%RHP(k,i,j)
  ENDIF
!
!           RAIN Q2
ELSEIF(VARN.EQ.'Q2') THEN
  IF(K.EQ.0)THEN
    IPLGRD=4
    FMT='F7.1'
    TILO='Q2(cal/g)'
  ELSE
    OPTLIB=micro_g(ngrid)%Q2(k,i,j)
  ENDIF
!
!           GRAUPEL Q6
ELSEIF(VARN.EQ.'Q6') THEN
  IF(K.EQ.0)THEN
    IPLGRD=4
    FMT='F7.1'
    TILO='Q6(cal/g)'
  ELSE
    OPTLIB=micro_g(ngrid)%Q6(k,i,j)
  ENDIF
!
!           HAIL Q7
ELSEIF(VARN.EQ.'Q7') THEN
  IF(K.EQ.0)THEN
    IPLGRD=4
    FMT='F7.1'
    TILO='Q7(cal/g)'
  ELSE
    OPTLIB=micro_g(ngrid)%Q7(k,i,j)
  ENDIF
!
!           LIQUID WATER MIXING RATIO
ELSEIF(VARN.EQ.'RL') THEN
  IF(K.EQ.0)THEN
    IPLGRD=4
    FMT='3PF7.3'
    TILO='RL(G/KG)'
  ELSE
    TCN=0.
    IF(LEVEL.GE.3) THEN
      IF(IPris.GT.0)  &
            TCN=TCN+micro_g(ngrid)%RPP(k,i,j)
      IF(ISnow.GT.0)  &
           TCN=TCN+micro_g(ngrid)%RSP(k,i,j)
      IF(IAggr.GT.0)  &
           TCN=TCN+micro_g(ngrid)%RAP(k,i,j)
      IF(IGraup.GT.0)  &
           TCN=TCN+micro_g(ngrid)%RGP(k,i,j)
      IF(Ihail.GT.0)  &
           TCN=TCN+micro_g(ngrid)%RHP(k,i,j)
    ENDIF
    OPTLIB=0.
    IF(LEVEL.GE.2)THEN
      VAL1=basic_g(ngrid)%RTP(k,i,j)
      VAL2=basic_g(ngrid)%RV(k,i,j)
      OPTLIB=VAL1-VAL2-TCN
    ENDIF
  ENDIF
!
!           TOTAL ICE MIXING RATIO
ELSEIF(VARN.EQ.'RI') THEN
  IF(K.EQ.0)THEN
    IPLGRD=4
    FMT='3PF7.3'
    TILO='RI(G/KG)'
  ELSE
    OPTLIB=0.
    IF(IPris.GT.0)  &
         OPTLIB=OPTLIB+micro_g(ngrid)%RPP(k,i,j)
    IF(ISnow.GT.0)  &
         OPTLIB=OPTLIB+micro_g(ngrid)%RSP(k,i,j)
    IF(IAggr.GT.0)  &
         OPTLIB=OPTLIB+micro_g(ngrid)%RAP(k,i,j)
    IF(IGraup.GT.0)  &
         OPTLIB=OPTLIB+micro_g(ngrid)%RGP(k,i,j)
    IF(Ihail.GT.0)  &
         OPTLIB=OPTLIB+micro_g(ngrid)%RHP(k,i,j)
  ENDIF
!
!           TOTAL CONDENSATE MIXING RATIO
ELSEIF(VARN.EQ.'RCOND') THEN
  IF(K.EQ.0)THEN
    IPLGRD=4
    FMT='3PF7.3'
    TILO='RD(G/KG)'
  ELSE
    OPTLIB=basic_g(ngrid)%RTP(k,i,j)  &
          -basic_g(ngrid)%RV(k,i,j)
  ENDIF
!
!
!           Cloud droplet CONCENTRATION
ELSEIF(VARN.EQ.'CC') THEN
  IF(K.EQ.0)THEN
    IPLGRD=4
    FMT='0PF7.2'
    TILO='NRAIN'
  ELSE
    OPTLIB=0.
    IF(Icloud.EQ.5)THEN
      OPTLIB=micro_g(ngrid)%CCP(k,i,j)
    ENDIF
  ENDIF
!
!           RAIN CONCENTRATION
ELSEIF(VARN.EQ.'CR') THEN
  IF(K.EQ.0)THEN
    IPLGRD=4
    FMT='0PF7.2'
    TILO='NRAIN'
  ELSE
    OPTLIB=0.
    IF(Irain.EQ.5)THEN
      OPTLIB=micro_g(ngrid)%CRP(k,i,j)
    ENDIF
  ENDIF
!
!           ICE CRYSTAL CONCENTRATION
ELSEIF(VARN.EQ.'CP') THEN
  IF(K.EQ.0)THEN
    IPLGRD=4
    FMT='0PF7.2'
    TILO='NPRIS'
  ELSE
    OPTLIB=0.
    IF(IPris.EQ.5)THEN
      OPTLIB=log10(max(1.,micro_g(ngrid)%CPP(k,i,j)))
    ENDIF
  ENDIF
!
!           snow CONCENTRATION
ELSEIF(VARN.EQ.'CS') THEN
  IF(K.EQ.0)THEN
    IPLGRD=4
    FMT='0PF7.2'
    TILO='NSNOW'
  ELSE
    OPTLIB=0.
    IF(Isnow.EQ.5)THEN
      OPTLIB=micro_g(ngrid)%CSP(k,i,j)
    ENDIF
  ENDIF
!
!           aggregate CONCENTRATION
ELSEIF(VARN.EQ.'CA') THEN
  IF(K.EQ.0)THEN
    IPLGRD=4
    FMT='0PF7.2'
    TILO='NAGGR'
  ELSE
    OPTLIB=0.
    IF(Iaggr.EQ.5)THEN
      OPTLIB=micro_g(ngrid)%CAP(k,i,j)
    ENDIF
  ENDIF
!
!           graupel CONCENTRATION
ELSEIF(VARN.EQ.'CG') THEN

  IF(K.EQ.0)THEN
    IPLGRD=4
    FMT='0PF7.2'
    TILO='NGRAUP'
  ELSE

    OPTLIB=0.
    IF(Igraup.EQ.5)THEN
      OPTLIB=micro_g(ngrid)%CGP(k,i,j)

    ENDIF
  ENDIF
!
!           hail CONCENTRATION
ELSEIF(VARN.EQ.'CH') THEN
  IF(K.EQ.0)THEN
    IPLGRD=4
    FMT='0PF7.2'
    TILO='NHAIL'
  ELSE
    OPTLIB=0.
    IF(Ihail.EQ.5)THEN
      OPTLIB=micro_g(ngrid)%CHP(k,i,j)
    ENDIF
  ENDIF
!
!           PERTURBATION TOTAL WATER
ELSEIF(VARN.EQ.'RTP') THEN
  IF(K.EQ.0)THEN
    IPLGRD=4
    FMT='3PF7.3'
    TILO='RT'' G/KG'
  ELSE
    VAL1=basic_g(ngrid)%RTP(k,i,j)
    VAL2=grid_g(ngrid)%TOPT(i,j)
    DO KK=1,NZP
      VCTR1(KK)=RT01DN(KK,NGRID)
    ENDDO
    CALL TRNCL2(VCTR1,ZT,VAL2,ZTOP,VCTR2,VCTR3,VCTR4,NZ)
    OPTLIB=VAL1-VCTR1(K)
  ENDIF
!
!           SUPERSATURATION/WATER
ELSEIF(VARN.EQ.'SUPSATW') THEN
  IF(K.EQ.0)THEN
    IPLGRD = 4
      FMT='2PF7.3'
      TILO='SupSatw %'
  ELSE
    VALTHET= basic_g(ngrid)%THETA(k,i,j)
    VALPP  = basic_g(ngrid)%PP(k,i,j)
    VALPI0 = basic_g(ngrid)%PI0(k,i,j)
    VALTEMP= VALTHET*(VALPP+VALPI0)/CP
    VALPRS = ((VALPI0 + VALPP)/CP)**CPOR * P00
    VALRV  = basic_g(ngrid)%RV(k,i,j)
    OPTLIB = 100.*REHUL(VALPRS,VALTEMP,VALRV)
    IF(OPTLIB > 100.)THEN
       OPTLIB=OPTLIB-100.
    ELSE
       OPTLIB=0.0
    ENDIF
  ENDIF
!
!           SUPERSATURATION/ICE
ELSEIF(VARN.EQ.'SUPSATI') THEN
  IF(K.EQ.0)THEN
    IPLGRD = 4
      FMT='2PF7.3'
      TILO='SupSati %'
  ELSE
    VALTHET= basic_g(ngrid)%THETA(k,i,j)
    VALPP  = basic_g(ngrid)%PP(k,i,j)
    VALPI0 = basic_g(ngrid)%PI0(k,i,j)
    VALTEMP= VALTHET*(VALPP+VALPI0)/CP
    VALPRS = ((VALPI0 + VALPP)/CP)**CPOR * P00
    VALRV  = basic_g(ngrid)%RV(k,i,j)
    OPTLIB = 100.*REHUI(VALPRS,VALTEMP,VALRV)
    IF(OPTLIB > 100.)THEN
       OPTLIB=OPTLIB-100.
    ELSE
       OPTLIB=0.0
    ENDIF
  ENDIF
!
!           RELATIVE HUMIDITY
ELSEIF(VARN.EQ.'RELHUM') THEN
  IF(K.EQ.0)THEN
    IPLGRD = 4
    FMT    = 'F7.1'
    TILO    = 'Rel Humid'
  ELSE
    VALTHET= basic_g(ngrid)%THETA(k,i,j)
    VALPP  = basic_g(ngrid)%PP(k,i,j)
    VALPI0 = basic_g(ngrid)%PI0(k,i,j)
    VALTEMP= VALTHET*(VALPP+VALPI0)/CP
    VALPRS = ((VALPI0 + VALPP)/CP)**CPOR * P00
    VALRV  = basic_g(ngrid)%RV(k,i,j)
    OPTLIB = 100. * rehuil(VALPRS,VALTEMP,VALRV)
  ENDIF
!
!           TOTAL NUMBER CONCENTRATION OF CCN
!
ELSEIF(VARN.EQ.'DNCCN') THEN
  IF(K.EQ.0)THEN
    IPLGRD=4
    FMT='-6PF7.1'
    TILO='DCCN(#/cc)'
  ELSE
    OPTLIB=0.
  IF(LEVEL.GE.4)THEN
!           DO L=NSH,NB+NSH-1
     DO L=2+2*ICLOUD,1+2*ICLOUD+LN2
       OPTLIB=OPTLIB+micro_g(ngrid)%cccnp(k,i,j)
     ENDDO
     OPTLIB=OPTLIB*basic_g(ngrid)%DN0(k,i,j)
     IF(OPTLIB.LT.1E+04)OPTLIB=0.
  ENDIF
  ENDIF
!
!
!           TOTAL NUMBER CONCENTRATION OF DROPS
!
ELSEIF(VARN.EQ.'DNTOT') THEN
  IF(K.EQ.0)THEN
    IPLGRD=4
    FMT='-6PF7.1'
    TILO='DN(#/cc)'
  ELSE
    OPTLIB=0.
  IF(LEVEL.EQ.4)THEN
     DO L=2,ICLOUD+1
       OPTLIB=OPTLIB+scalar_tab(l,ngrid)%var_p
     ENDDO
     OPTLIB=OPTLIB*basic_g(ngrid)%DN0(k,i,j)
     IF(OPTLIB.LT.1E+04)OPTLIB=0.
  ENDIF
  IF(LEVEL.GT.4)THEN
     DO L=NSH,ICLOUD+NSH-1
       OPTLIB=OPTLIB+scalar_tab(l,ngrid)%var_p
     ENDDO
     OPTLIB=OPTLIB*basic_g(ngrid)%DN0(k,i,j)
     IF(OPTLIB.LT.1E+03)OPTLIB=0.
  ENDIF
  ENDIF
!
!
!           TOTAL MASS CONCENTRATION OF DROPS
!
ELSEIF(VARN.EQ.'DMTOT') THEN
  IF(K.EQ.0)THEN
    IPLGRD=4
    FMT='3PF7.3'
    TILO='DM(G/KG)'
  ELSE
    OPTLIB=0.
  IF(LEVEL.EQ.4)THEN
     DO L=ICLOUD+2,2*ICLOUD+1
        OPTLIB=OPTLIB+scalar_tab(l,ngrid)%var_p
    ENDDO
    ENDIF
  IF(LEVEL.GT.4)THEN
     DO L=ICLOUD+NSH,2*ICLOUD+NSH-1
        OPTLIB=OPTLIB+scalar_tab(l,ngrid)%var_p
    ENDDO
    ENDIF
  ENDIF
!
!
!           TOTAL NUMBER CONCENTRATION OF ACTIVATED ICE NUCLEI
!
ELSEIF(VARN.EQ.'INACT') THEN
  IF(K.EQ.0)THEN
    IPLGRD=4
    FMT='-3PF7.1'
    TILO='AIN(#/l)'
  ELSE
    OPTLIB=0.
  IF(LEVEL.GT.4)THEN
       OPTLIB=OPTLIB+scalar_tab(nsin,ngrid)%var_p
     OPTLIB=OPTLIB*basic_g(ngrid)%DN0(k,i,j)
     IF(OPTLIB.LT.1E+01)OPTLIB=0.
  ENDIF
  ENDIF
!
!
!           TOTAL NUMBER CONCENTRATION OF PRISTINE ICE
!
ELSEIF(VARN.EQ.'PNTOT') THEN
  IF(K.EQ.0)THEN
    IPLGRD=4
    FMT='-3PF7.1'
    TILO='PN(#/l)'
  ELSE
    OPTLIB=0.
  IF(LEVEL.GT.4)THEN
     DO L=NSH+2*ICLOUD,NSH+2*ICLOUD+IPRIS-1
       OPTLIB=OPTLIB+scalar_tab(l,ngrid)%var_p
     ENDDO
     OPTLIB=OPTLIB*basic_g(ngrid)%DN0(k,i,j)
     IF(OPTLIB.LT.1E+01)OPTLIB=0.
  ENDIF
  ENDIF
!
!
!           TOTAL MASS CONCENTRATION OF PRISTINE ICE
!
ELSEIF(VARN.EQ.'PMTOT') THEN
  IF(K.EQ.0)THEN
    IPLGRD=4
    FMT='3PF7.3'
    TILO='PM(G/KG)'
  ELSE
    OPTLIB=0.
  IF(LEVEL.GT.4)THEN
     DO L=NSH+2*ICLOUD+IPRIS,NSH+2*ICLOUD+2*IPRIS-1
        OPTLIB=OPTLIB+scalar_tab(l,ngrid)%var_p
    ENDDO
    ENDIF
  ENDIF
!
!           TOTAL NUMBER CONCENTRATION OF AGGREGATES
!
ELSEIF(VARN.EQ.'ANTOT') THEN
  IF(K.EQ.0)THEN
    IPLGRD=4
    FMT='-3PF7.1'
    TILO='AN(#/l)'
  ELSE
    OPTLIB=0.
  IF(LEVEL.GT.4)THEN
     DO L=NSH+2*ICLOUD+2*IPRIS,NSH+2*ICLOUD+2*IPRIS+IAGGR-1
       OPTLIB=OPTLIB+scalar_tab(l,ngrid)%var_p
     ENDDO
     OPTLIB=OPTLIB*basic_g(ngrid)%DN0(k,i,j)
     IF(OPTLIB.LT.1E+01)OPTLIB=0.
  ENDIF
  ENDIF
!
!
!           TOTAL MASS CONCENTRATION OF AGGREGATES
!
ELSEIF(VARN.EQ.'AMTOT') THEN
  IF(K.EQ.0)THEN
    IPLGRD=4
    FMT='3PF7.3'
    TILO='AM(G/KG)'
  ELSE
    OPTLIB=0.
  IF(LEVEL.GT.4)THEN
     DO L=NSH+2*ICLOUD+2*IPRIS+IAGGR,NSH+2*ICLOUD+2*IPRIS  &
                      +2*IAGGR-1
        OPTLIB=OPTLIB+scalar_tab(l,ngrid)%var_p
    ENDDO
    ENDIF
  ENDIF
!
!           TOTAL NUMBER CONCENTRATION OF GRAUPEL
!
ELSEIF(VARN.EQ.'GNTOT') THEN
  IF(K.EQ.0)THEN
    IPLGRD=4
    FMT='-3PF7.1'
    TILO='GN(#/l)'
  ELSE
    OPTLIB=0.
  IF(LEVEL.GT.4)THEN
     DO L=NSH+2*ICLOUD+2*IPRIS+2*IAGGR,NSH+2*ICLOUD+2*IPRIS  &
                      +2*IAGGR+IGRAUP-1
       OPTLIB=OPTLIB+scalar_tab(l,ngrid)%var_p
     ENDDO
     OPTLIB=OPTLIB*basic_g(ngrid)%DN0(k,i,j)
     IF(OPTLIB.LT.1E+01)OPTLIB=0.
  ENDIF
  ENDIF
!
!
!           TOTAL MASS CONCENTRATION OF GRAUPEL
!
ELSEIF(VARN.EQ.'GMTOT') THEN
  IF(K.EQ.0)THEN
    IPLGRD=4
    FMT='3PF7.3'
    TILO='GM(G/KG)'
  ELSE
    OPTLIB=0.
  IF(LEVEL.GT.4)THEN
     DO L=NSH+2*ICLOUD+2*IPRIS+2*IAGGR+IGRAUP,NSH+2*ICLOUD  &
                      +2*IPRIS+2*IAGGR+2*IGRAUP-1
        OPTLIB=OPTLIB+scalar_tab(l,ngrid)%var_p
    ENDDO
    ENDIF
  ENDIF
!
!           TOTAL NUMBER CONCENTRATION OF ICE
!
ELSEIF(VARN.EQ.'INTOT') THEN
  IF(K.EQ.0)THEN
    IPLGRD=4
    FMT='-3PF7.1'
    TILO='IN(#/l)'
  ELSE
    OPTLIB=0.
  IF(LEVEL.GT.4)THEN
     DO L=NSH+2*ICLOUD,NSH+2*ICLOUD+IPRIS-1
       OPTLIB=OPTLIB+scalar_tab(l,ngrid)%var_p
     ENDDO
     DO L=NSH+2*ICLOUD+2*IPRIS,NSH+2*ICLOUD+2*IPRIS+IAGGR-1
       OPTLIB=OPTLIB+scalar_tab(l,ngrid)%var_p
     ENDDO
     DO L=NSH+2*ICLOUD+2*IPRIS+2*IAGGR,NSH+2*ICLOUD+2*IPRIS  &
                      +2*IAGGR+IGRAUP-1
       OPTLIB=OPTLIB+scalar_tab(l,ngrid)%var_p
     ENDDO
     OPTLIB=OPTLIB*basic_g(ngrid)%DN0(k,i,j)
     IF(OPTLIB.LT.1E+01)OPTLIB=0.
  ENDIF
  ENDIF
!
!
!           TOTAL MASS CONCENTRATION OF ICE
!
ELSEIF(VARN.EQ.'IMTOT') THEN
  IF(K.EQ.0)THEN
    IPLGRD=4
    FMT='3PF7.3'
    TILO='IM(G/KG)'
  ELSE
    OPTLIB=0.
  IF(LEVEL.GT.4)THEN
     DO L=NSH+2*ICLOUD+IPRIS,NSH+2*ICLOUD+2*IPRIS-1
        OPTLIB=OPTLIB+scalar_tab(l,ngrid)%var_p
    ENDDO
     DO L=NSH+2*ICLOUD+2*IPRIS+IAGGR,NSH+2*ICLOUD+2*IPRIS  &
                      +2*IAGGR-1
        OPTLIB=OPTLIB+scalar_tab(l,ngrid)%var_p
    ENDDO
     DO L=NSH+2*ICLOUD+2*IPRIS+2*IAGGR+IGRAUP,NSH+2*ICLOUD  &
                      +2*IPRIS+2*IAGGR+2*IGRAUP-1
        OPTLIB=OPTLIB+scalar_tab(l,ngrid)%var_p
    ENDDO
    ENDIF
  ENDIF
!
!           TOTAL ACCUMULATED PRECIPITATION
ELSEIF(VARN.EQ.'PCPT') THEN
  IF(K.EQ.0)THEN
    IPLGRD=4
    FMT='0PF7.1'
    TILO='TOTPRE'
  ELSE
    OPTLIB=0.
    IF(IRAIN.GE.1)  &
         OPTLIB=OPTLIB+micro_g(ngrid)%ACCPR(i,j)
    IF(IPRIS.GE.1)  &
         OPTLIB=OPTLIB+micro_g(ngrid)%ACCPP(i,j)
    IF(ISNOW.GE.1)  &
         OPTLIB=OPTLIB+micro_g(ngrid)%ACCPS(i,j)
    IF(IAGGR.GE.1)  &
         OPTLIB=OPTLIB+micro_g(ngrid)%ACCPA(i,j)
    IF(IGRAUP.GE.1)  &
         OPTLIB=OPTLIB+micro_g(ngrid)%ACCPG(i,j)
    IF(IHAIL.GE.1)  &
         OPTLIB=OPTLIB+micro_g(ngrid)%ACCPH(i,j)
  ENDIF
!
!           TKE
ELSEIF(VARN.EQ.'TKE') THEN
  IF(K.EQ.0)THEN
    IPLGRD=4
    FMT=' 0PF7.2'
    TILO='TKE'
  ELSE
      OPTLIB=turb_g(ngrid)%TKEP(k,i,j)
  ENDIF
!_STC.....................................................................
!_STC    Extraction for added scalar sclp(1) = tkeps = tke dissipation
!_STC    E-l and E-eps closures
!_STC    (S. Trini Castelli)
!_STC.....................................................................

!           TKEPS (TKE dissipation)
ELSEIF(VARN.EQ.'EPS') THEN
  IF(K.EQ.0)THEN
    IPLGRD=4
    FMT=' 0PF10.7'
    TILO='TKEPS '
  ELSE
    OPTLIB=turb_g(ngrid)%epsp(k,i,j)
  ENDIF
!_STC.....................................................................
!
!
!           Surface downward shortwave radiation
ELSEIF(VARN.EQ.'RSHORT') THEN
  IF(K.EQ.0)THEN
    IPLGRD=4
    FMT='0PF7.2'
    TILO='SW RAD'
  ELSE
    OPTLIB=radiate_g(ngrid)%RSHORT(i,j)
  ENDIF
!
!           SST
ELSEIF(VARN.EQ.'SST') THEN
  IF(K.EQ.0)THEN
    IPLGRD=4
    FMT='0PF7.2'
    TILO='SST (K)'
  ELSE
!    OPTLIB=A(ITGP+(NZG-1)*nxp*nyp)
  ENDIF
!
!           Surface downward longwave radiation
ELSEIF(VARN.EQ.'RLONG') THEN
  IF(K.EQ.0)THEN
    IPLGRD=4
    FMT='0PF7.2'
    TILO='LW RAD'
  ELSE
    OPTLIB=radiate_g(ngrid)%RLONG(i,j)
  ENDIF
!
!           Latitude
ELSEIF(VARN.EQ.'GLAT') THEN
  IF(K.EQ.0)THEN
    IPLGRD=4
    FMT='0PF7.2'
    TILO='Lat'
  ELSE
    OPTLIB=grid_g(ngrid)%GLAT(i,j)
  ENDIF
!
!           Longitude
ELSEIF(VARN.EQ.'GLON') THEN
  IF(K.EQ.0)THEN
    IPLGRD=4
    FMT='0PF7.2'
    TILO='Lon'
  ELSE
    OPTLIB=grid_g(ngrid)%GLON(i,j)
  ENDIF
!
!           Topography
ELSEIF(VARN.EQ.'TOPT') THEN
  IF(K.EQ.0)THEN
    IPLGRD=4
    FMT='0PF7.1'
    TILO='Topo'
  ELSE
    OPTLIB=grid_g(ngrid)%TOPT(i,j)
  ENDIF
!
!           Topma
ELSEIF(VARN.EQ.'TOPMA') THEN
  IF(K.EQ.0)THEN
    IPLGRD=4
    FMT='0PF7.1'
    TILO='Topma'
  ELSE
    OPTLIB=grid_g(ngrid)%TOPMA(i,j)
  ENDIF
!
!           Seatf
ELSEIF(VARN.EQ.'SEATF') THEN
  IF(K.EQ.0)THEN
    IPLGRD=4
    FMT='0PF7.2'
    TILO='Seatf'
  ELSE
    OPTLIB=leaf_g(ngrid)%SEATF(i,j)
  ENDIF
!
!           Seatp
ELSEIF(VARN.EQ.'SEATP') THEN
  IF(K.EQ.0)THEN
    IPLGRD=4
    FMT='0PF7.2'
    TILO='Seatp'
  ELSE
    OPTLIB=leaf_g(ngrid)%SEATP(i,j)
  ENDIF
!
!           CONVECTIVE PRECIP RATE
ELSEIF(VARN.EQ.'CONPR') THEN
  IF(K.EQ.0)THEN
    IPLGRD=4
    FMT='0PF7.2'
    TILO='CON RATE'
  ELSE
    VAL1=0.
    DO ICLD=1,NCLOUDS
       VAL1=VAL1+cuparm_g(ngrid)%CONPRR(i,j,icld)
    END DO
    OPTLIB=VAL1*hr_sec
  ENDIF
!
!           ACCUMULATED CONVECTIVE PRECIP
ELSEIF(VARN.EQ.'CONP') THEN
  IF(K.EQ.0)THEN
    IPLGRD=4
    FMT='0PF7.2'
    TILO='CON PCP'
  ELSE
    OPTLIB=cuparm_g(ngrid)%ACONPR(i,j)
  ENDIF
!
!           CONVECTIVE HEATING RATE
ELSEIF(VARN.EQ.'CONH') THEN
  IF(K.EQ.0)THEN
    IPLGRD=4
    FMT='0PF7.1'
    TILO='CON HEAT'
  ELSE
    VAL1=0.
    DO ICLD=1,NCLOUDS
       VAL1=VAL1+cuparm_g(ngrid)%THSRC(k,i,j,icld)
    END DO
    OPTLIB=VAL1*day_sec
  ENDIF
!
!           CONVECTIVE MOISTENING RATE
ELSEIF(VARN.EQ.'CONM') THEN
  IF(K.EQ.0)THEN
    IPLGRD=4
    FMT='0PF7.1'
    TILO='CON MOIS'
  ELSE
    VAL1=0.
    DO ICLD=1,NCLOUDS
       VAL1=VAL1+cuparm_g(ngrid)%RTSRC(k,i,j,icld)
    END DO
    OPTLIB=VAL1*day_sec
  ENDIF
!
!
!     MICROPHYSICS DEPICTION
!
ELSEIF(VARN.EQ.'MICRO') THEN
  IF(K.EQ.0)THEN
    IPLGRD=4
    FMT='0PF7.0'
    TILO='GASPRC'
  ELSE
    OPTLIB=0.
    IF(LEVEL.GE.2)THEN
      PRTC=0.
      PRTR=0.
      PRTP=0.
      PRTS=0.
      PRTA=0.
      PRTG=0.
      IF(LEVEL.GE.3) THEN
        IF(IRAIN.GE.1)  PRTR=micro_g(ngrid)%RRP(k,i,j)
        IF(IPRIS.GE.1)  PRTP=micro_g(ngrid)%RPP(k,i,j)
        IF(ISNOW.GE.1)  PRTS=micro_g(ngrid)%RSP(k,i,j)
        IF(IAGGR.GE.1)  PRTA=micro_g(ngrid)%RAP(k,i,j)
        IF(IGRAUP.GE.1) PRTG=micro_g(ngrid)%RGP(k,i,j)
        IF(IHAIL.GE.1)  PRTH=micro_g(ngrid)%RHP(k,i,j)
      ENDIF
      PRALL=basic_g(ngrid)%RTP(k,i,j)
      PRVAP=basic_g(ngrid)%RV(k,i,j)
      PRTC=MAX(0.,PRALL-PRVAP-PRTR-PRTP-PRTS-PRTA-PRTG-PRTH)

      PRTC=MIN(9,INT((PRTC+.9999E-3)*1.E3))
      PRTR=MIN(9,INT((PRTR+.9999E-3)*1.E3))
      PRTP=MIN(9,INT((PRTP+.9999E-3)*1.E3))
      PRTS=MIN(9,INT((PRTS+.9999E-3)*1.E3))
      PRTA=MIN(9,INT((PRTA+.9999E-3)*1.E3))
      PRTG=MIN(9,INT((PRTG+.9999E-3)*1.E3))
      PRTH=MIN(9,INT((PRTH+.9999E-3)*1.E3))
      IPRT=PRTC*1.  &
          +PRTR*10.  &
          +PRTP*100.  &
          +PRTS*1000.  &
          +PRTA*10000.  &
          +PRTG*100000.
      OPTLIB=IPRT
    ENDIF
  ENDIF
!
!           WIND SPEED
ELSEIF(VARN.EQ.'SPEED') THEN
   IF (K.EQ.0) THEN
      IPLGRD = 4
      FMT = '0PF7.1'
      TILO = 'SPEED'
   ELSE
      U1 = basic_g(ngrid)%UP(k,i,j)
      U2 = basic_g(ngrid)%UP(k,i-1,j)
      V1 = basic_g(ngrid)%VP(k,i,j)
      V2 = V1
      IF (JDIM.EQ.1) V2 = basic_g(ngrid)%VP(k,i,j-1)
      OPTLIB = SQRT( ((U1+U2)*0.5)**2 + ((V1+V2)*0.5)**2)
   ENDIF
!
!           BRUNT-VAISALA FREQUENCY
!      ELSEIF(VARN.EQ.'BVF') THEN
!           RICHARDSON NUMBER
!      ELSEIF(VARN.EQ.'RICH') THEN
!           VERTICAL EDDY VISCOSITY FOR MOMENTUM
!      ELSEIF(VARN.EQ.'VKM') THEN
!           VERTICAL EDDY VISCOSITY FOR HEAT
!      ELSEIF(VARN.EQ.'VKH') THEN
!           HORIZONTAL EDDY VISCOSITY FOR MOMENTUM
!      ELSEIF(VARN.EQ.'HKM') THEN
!           HORIZONTAL EDDY VISCOSITY FOR HEAT
!      ELSEIF(VARN.EQ.'HKH') THEN
!
!            RADIATIVE FLUX CONVERGENCE
ELSEIF(VARN.EQ.'FTHRD') THEN
   IF(K.EQ.0)THEN
      IPLGRD=4
      FMT='0PF6.3'
      TILO='FTHRD  '
   ELSE
      OPTLIB = radiate_g(ngrid)%FTHRD(k,i,j)
   ENDIF
!
!           PBL HEIGHT  ZI
!      ELSEIF(VARN.EQ.'ZI') THEN
!
!           SURFACE VERTICAL U-MOMENTUM FLUX  SFLUX_U
ELSEIF(VARN.EQ.'SFLUX_U') THEN
   IF(K.EQ.0)THEN
      IPLGRD=4
      FMT='0PF7.3'
      TILO='SFC UM FLUX '
   ELSE
      OPTLIB = turb_g(ngrid)%sflux_u(i,j)
   ENDIF
!
!           SURFACE VERTICAL V-MOMENTUM FLUX  SFLUX_V
ELSEIF(VARN.EQ.'SFLUX_V') THEN
   IF(K.EQ.0)THEN
      IPLGRD=4
      FMT='0PF7.3'
      TILO='SFC VM FLUX '
   ELSE
      OPTLIB = turb_g(ngrid)%sflux_v(i,j)
   ENDIF
!
!           SURFACE VERTICAL W-MOMENTUM FLUX  SFLUX_W
ELSEIF(VARN.EQ.'SFLUX_W') THEN
   IF(K.EQ.0)THEN
      IPLGRD=4
      FMT='0PF7.3'
      TILO='SFC WM FLUX '
   ELSE
      OPTLIB = turb_g(ngrid)%sflux_w(i,j)
   ENDIF
!
!           SURFACE VERTICAL THETA*MASS FLUX  SFLUX_T
ELSEIF(VARN.EQ.'SFLUX_T') THEN
   IF(K.EQ.0)THEN
      IPLGRD=4
      FMT='0PF7.3'
      TILO='SFC THM FLUX '
   ELSE
      OPTLIB = turb_g(ngrid)%sflux_t(i,j)
   ENDIF
!
!           SURFACE VERTICAL VAPOR SFLUX_R
   ELSEIF(VARN.EQ.'SFLUX_R') THEN
   IF(K.EQ.0)THEN
      IPLGRD=4
      FMT='0PF7.3'
      TILO='SFC VAP FLUX '
   ELSE
      OPTLIB = turb_g(ngrid)%sflux_r(i,j)
   ENDIF
!
ENDIF
!
RETURN
END
!
!  *********************************************************************
!
SUBROUTINE PRTOPT(M)

use mem_grid
use mem_scratch
use ref_sounding
use rconstants
use mem_turb
use ref_sounding
use therm_lib, only: rehuil,tv2temp

IF(INITIAL.NE.2)THEN
  do k=1,nsndg
     vctr1(k) = 100. *rehuil(ps(k),ts(k),rts(k))
  end do
  WRITE(M,41)
41      FORMAT(/,'------------------------------SOUNDING INPUT-------'  &
   ,'---------------------------',//,7X,'PS',9X,'HS',7X,'TS',6X  &
   ,'THDS',6X,'US',7X,'VS',7X,'RTS',5X,'REL HUM',/,6X  &
   ,'(Pa)',7X,'(m)',6X,'(K)',6X,'(K)',6X,'(m/s)',4X,'(m/s)'  &
   ,3X,'(kg/kg)',5X,'(%)',/)
  WRITE(M,42)(PS(K),HS(K),TS(K),THDS(K),US(K),VS(K),RTS(K)  &
             ,VCTR1(K),K=1,NSndg)
42      FORMAT(1X,F11.1,F10.1,2F9.2,2F9.2,F10.5,F9.1)
ENDIF
!
DO K=1,NNZP(1)
  VCTR1(K)=P00*(PI01DN(K,1)/CP)**CPOR
  VCTR2(K)=tv2temp(TH01DN(K,1),RT01DN(K,1))
ENDDO
WRITE(M,310)IREF,JREF,TOPREF,(ZTN(K,1),U01DN(K,1),V01DN(K,1)  &
  ,DN01DN(K,1),PI01DN(K,1),VCTR1(K),TH01DN(K,1),VCTR2(K)  &
  ,RT01DN(K,1),CO201DN(K,1),K=1,NNZP(1))
310   FORMAT(/,'--------REFERENCE STATE at I,J=(',I4,',',I4  &
      ,')   SFC ELEV (M)= ',F6.1,'-------------'  &
 ,//,4X,'Z',6X,'U01D',4X,'V01D',4X,'DN01D',4X  &
 ,'PI01D',4X,'PRESS',4X,'TH01D',4X,'THD',6X,'RT01D',4X,'CO201D'  &
 ,/,3X,'(m)',5X,'(m/s)',3X,'(m/s)',2X,'(kg/m3)',2X  &
 ,'(J/kgK)',4X,'(Pa)',5X,'(K)',5X,'(K)',5X,'(kg/kg)',3X,'(ppm)'  &
 ,//,(1X,F7.1,2F8.2,F8.3,F10.2,F10.1,2F8.2,F10.5,F9.3))
!
!STC..................................................................
!  Print out the values chosen for the empirical constants in
!  E-l and E-eps closures (they are the same for all the grids)
!  (S. Trini Castelli)
!STC..................................................................
!
WRITE(6,*) ' '
if (IDIFFK(1).eq.5) then
WRITE(6,*) 'Empirical constants and other parameters for E-l closure'
WRITE(6,105) ' ',C_MI,C_EPS,ALF_THT,ALF_TKE,ALF_EPS
if(IOPZL.eq.1) WRITE(6,*)  &
 'IOPZL=1 Constant asymptotic mixing length',AL0_CONST,' from Ying'
if(IOPZL.eq.2) WRITE(6,*)  &
 'IOPZL=2 Asymptotic mixing length from Mellor-Yamada'
if(IOPZL.eq.3) WRITE(6,*)  &
 'IOPZL=3 Asymptotic mixing length from Zilitinkevich'
endif
if (IDIFFK(1).eq.6) then
WRITE(6,*) 'Empirical constants and other parameters for E-eps closure'
WRITE(6,106) ' ',C_MI,C1_EPS,C2_EPS,ALF_THT,ALF_TKE,ALF_EPS
endif

105  FORMAT(A1,'C_MI=',f4.2,'  C_EPS=',f4.2,'  ALF_THT=',f4.2  &
      ,'  ALF_TKE=',f4.2,'  ALF_EPS=',f4.2,999(A1,/,I13,3I16))
106  FORMAT(A1,'C_MI=',f4.2,'  C1_EPS=',f4.2,'  C2_EPS=',f4.2,  &
     '  ALF_THT=',f4.2,'  ALF_TKE=',f4.2,'  ALF_EPS=', &
      f4.2,999(A1,/,I13,3I16))
!STC.................................................................                        

RETURN
END
!
!     ******************************************************************
!
SUBROUTINE UWCOMP()

use mem_basic
use mem_grid

call newgrid(1)
CALL UWC(NZP,NXP,NYP,basic_g(ngrid)%UP,basic_g(ngrid)%WP,grid_g(ngrid)%DXT)
!
RETURN
END
!
!     ******************************************************************
!
SUBROUTINE UWC(N1,N2,N3,UP,WP,DXT)

use mem_grid
use mem_scratch
use ref_sounding
use rconstants

DIMENSION UP(N1,N2,N3),WP(N1,N2,N3),DXT(N2,N3)
DIMENSION UMN(NZPMAX,NXPMAX),WMN(NZPMAX,NXPMAX),V1(NXPMAX)
!
DO J=1,NYP
   DO I=2,NX
!
      DO K=1,NZ
         UMN(K,I)=.5*(UP(K,I,J)+UP(K,I-1,J))
         WMN(K,I)=.5*(WP(K,I,J)+WP(K-1,I,J))
      ENDDO
      WMN(1,I)=WP(1,I,J)
!
   ENDDO
ENDDO
!
IR=NX
IL=NX
NTPT=IR-IL+1
J=1
DO K=1,NZ
   VCTR1(K)=SSUM(NTPT,UMN(K,IL),NZPMAX)/(NTPT)
   VCTR2(K)=SSUM(NTPT,WMN(K,IL),NZPMAX)/(NTPT)
   VCTR3(K)=TH01DN(K,NGRID)
   VCTR3(K)=1./DN01DN(K,NGRID)
ENDDO
!
PRINT 89
89   FORMAT('    K',7X,'Z',5X,'PCT',6X,'UWFLX',8X,'TVAL',  &
     4X,'UBARP',4X,'WBARP')
DO K=1,NZ
   DO I=IL,IR
      V1(I)=(UMN(K,I)-VCTR1(K))*(WMN(K,I)-VCTR2(K))  &
           /(VCTR3(K)*DXT(I,J))
   ENDDO
   UWFLX=-SSUM(NTPT,V1(IL),1)
   TVAL=pio4*20.*.01957*(10.**2)/VCTR3(1)
   PRINT 90,K,ZT(K),UWFLX/TVAL*100.,UWFLX,TVAL,VCTR1(K)-20.  &
        ,VCTR2(K)
90      FORMAT(I5,F8.0,F8.1,2E12.4,2F9.5)
ENDDO
!
RETURN
END


!     **************************************************************
!
SUBROUTINE PRT2D(HORIZ,VERT,A,X,Y,IX,IY,I1,I2,J1,J2  &
                ,IFMT,TITLE,XX,YY,XLABL,YLABL)
CHARACTER*(*) IFMT,TITLE,XLABL,YLABL,HORIZ,VERT
CHARACTER*8 FMTX,FMTY
CHARACTER*133 B
CHARACTER*16 XLAB,YLAB,FMTB,FMTC
CHARACTER*80 FMTT,FMTT2
DIMENSION A(IX,IY),X(IX),Y(IY),XX(IX),YY(IY)
!
! for 15 columns      NCOL=105
NCOL=84
IFDW=7
MMAX=NCOL/IFDW
!                 See if the field is horizontally homogeneous (along
!                   the abscissa)
NDIFF=0
DO I=I1,I2
  XX(I-I1+1)=X(I)
  DO J=J1,J2
    YY(J-J1+1)=Y(J)
    A(I-I1+1,J-J1+1)=A(I,J)
    IF(A(I,J).NE.A(1,J-J1+1)) NDIFF=1
  ENDDO
ENDDO
XXMAX=MAX(ABS(X(I1)),ABS(X(I2)))
IF(XXMAX.LT.10.) FMTX='F7.3'
IF(XXMAX.GE.10.AND.XXMAX.LT.100.) FMTX='F7.2'
IF(XXMAX.GE.100.AND.XXMAX.LT.1000.) FMTX='F7.1'
IF(XXMAX.GE.1000.) FMTX='F7.0'
YYMAX=MAX(ABS(Y(J1)),ABS(Y(J2)))
IF(YYMAX.LT.10.) FMTY='F7.3'
IF(YYMAX.GE.10.AND.YYMAX.LT.100.) FMTY='F7.2'
IF(YYMAX.GE.100.AND.YYMAX.LT.1000.) FMTY='F7.1'
IF(YYMAX.GE.1000.) FMTY='F7.0'

!                  Set print window and number of pages accordingly
IF(NDIFF.EQ.0) THEN
  II1=1
  II2=1
  JJ1=1
  JJ2=J2-J1+1
  NPAGES=1
  NPPGE=MIN(MMAX,II2)
  IA=II1
  IB=II1
ELSE
  II1=1
  II2=I2-I1+1
  JJ1=1
  JJ2=J2-J1+1
  NPPGE=MIN(MMAX,II2)
  NPAGES=(II2-II1)/MMAX+1
  IA=II1
  IB=IA+NPPGE-1
ENDIF
!                      Set up formats for the coordinate line
IF(II1.EQ.II2)THEN
  FMTT='(/,1X,''CONSTANT ALONG ABSCISSA'',/)'
ELSE
  WRITE(FMTT,11) NPPGE,FMTX(1:4)
11      FORMAT('(/,1X,''PAGE : '',I3,//,2X,A1,''/'',A1,4X,',I3,A4,')')
  WRITE(FMTT2,12) NPPGE,FMTX(1:4)
12      FORMAT('(2X,A1,''/'',A1,4X,',I3,A4,')')
ENDIF
!
XLAB(1:8)=XLABL
YLAB(1:8)=YLABL
!                  Print top banner line
WRITE(6,220)
220 FORMAT(5(/),1X,'+',130('-'))
WRITE(6,21) TITLE,YLAB(1:8),XLAB(1:8)
21 FORMAT(' ! ',A68,'   Ordinate:',A8,'  Abscissa:',A8)
WRITE(6,221)
221 FORMAT(1X,'+',130('-'))
!
!                  Loop through number of pages
!
DO 1000 IPAGE=1,NPAGES
!
  IF(II1.EQ.II2)THEN
    WRITE(6,FMTT)
  ELSE
    WRITE(6,FMTT) IPAGE,VERT,HORIZ,(XX(II)+.00001,II=IA,IB)
  ENDIF
!
  B=' '
  WRITE(6,132)
  DO JJ=JJ1,JJ2
    JROW=JJ2-JJ+JJ1
    ICHST=1
    FMTC='(1X,'//FMTY//',A1)'
    WRITE(B(ICHST:ICHST+8),FMTC) YY(JROW)+.00001,'!'
    ICHST=ICHST+9
    DO II=IA,IB
      IF(A(II,JROW).EQ.0.)THEN
        B(ICHST:ICHST+IFDW-1)=' '
      ELSE
        FMTC='('//IFMT//')'
        WRITE(B(ICHST:ICHST+IFDW-1),FMTC) A(II,JROW)
      ENDIF
      ICHST=ICHST+IFDW
    ENDDO
    FMTC='(A1,'//FMTY//')'
    WRITE(B(ICHST:ICHST+7),FMTC) '!',YY(JROW)+.00001
    ICHEND=ICHST+7
!
    WRITE(FMTB,233) ICHEND
233     FORMAT('(A',I3,')')
    WRITE(6,FMTB) B(1:ICHEND)
!
  ENDDO
!
  WRITE(6,132)
132   FORMAT(8X,'+',105('-'),'+')
!
  IF(II1.NE.II2)THEN
    WRITE(6,FMTT2) VERT,HORIZ,(XX(II)+.00001,II=IA,IB)
  ENDIF
!
  IB=MIN(IB+NPPGE,II2)
  IA=IA+NPPGE
1000  CONTINUE
!
RETURN
END



