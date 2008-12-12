!========================================================================
! Soil Moisture Estimate for NWP Models
! Coded and implemented by Rodrigo Gevaerd and Saulo Freitas
! Ref.: Gevaerd, R. e S. R. Freitas, Estimativa operacional da umidade 
! do solo para iniciacao de modelos de previsao numerica da atmosfera. 
! Parte I: Descricao da metodologia e validacao. Rev. Bras. Meteo.,
! volume especial do LBA, 2007.
!========================================================================

subroutine soil_moisture_init(n1,n2,n3,mzg,mzs,npat,ifm   &
     ,theta,pi0,pp  &
     ,soil_water     ,soil_energy      ,soil_text          &
     ,sfcwater_mass  ,sfcwater_energy  ,sfcwater_depth     &
     ,glat,glon  &
     ,flpw   &
     )

  use mem_grid, only : RUNTYPE,  &          ! INTENT(IN)
       iyeara, imontha, idatea, itimea      ! INTENT(IN)

  use mem_soil_moisture, only : SOIL_MOIST,   &  ! INTENT(IN)
       SOIL_MOIST_FAIL,                       &  ! INTENT(IN)
       usdata_in,                             &  ! INTENT(IN)
       usmodel_in                                ! INTENT(IN)

  use io_params, only : timstr

  use rconstants, only : cpi,alli1000,cliq1000,t00       ! INTENT(IN)

  use leaf_coms, only : soilcp,  & ! INTENT(IN)
       slmsts,                   & ! INTENT(IN)
       slcpd                       ! INTENT(IN)

  use mem_leaf, only : stgoff,   & ! INTENT(IN)
       slmstr,                   & ! INTENT(IN)
       slz                         ! INTENT(IN)

  implicit none
  !---------
  integer :: n1,n2,n3,mzg,mzs,npat,ifm,i,j,k,ipat,nveg,nsoil

  real :: c1,airtemp, pis

  real(kind=8) :: dif_time
  real :: seconds

  real, dimension(n1,n2,n3) :: theta,pi0,pp
  real, dimension(n2,n3)    :: glat,glon
  real, dimension(n2,n3) :: flpw

  real, dimension(mzg,n2,n3,npat) :: soil_water,soil_energy,soil_text
  real, dimension(mzs,n2,n3,npat) :: sfcwater_mass,sfcwater_energy  &
       ,sfcwater_depth

  !---------
  integer :: qi1,qi2,qj1,qj2,ncount,  &
       ii,jj,jc,ic,i1,j1,i2,j2,kk,ifname,k2,ipref,ipref_start,icihourmin

  integer :: n4us
  integer :: nlat, nlon
  real, allocatable :: slz_us(:)  
  

  real ::  latni,latnf,lonni,lonnf,ilatn,ilonn     &
       ,ilats,ilons,latn,lonn,lats,lons,dlatr,dlonr
  logical there,theref
  character (len=256) :: usdata, usmodel
  character (len=20)  :: pref
  character (len=2)  :: cidate,cimon
  character (len=1)  :: cgrid
  character (len=4)  :: ciyear
  character (len=4)  :: cihourmin
  real, allocatable :: api_us(:,:,:),prlat(:,:),prlon(:,:),usdum(:)

  integer :: int_dif_time, idate2, imonth2, iyear2, hourmin, da
  logical :: sair

  logical :: general
  namelist /gradeumso/ latni, latnf, lonni, lonnf, ilatn, ilonn, nlat, nlon

  iyear2  = iyeara
  imonth2 = imontha
  idate2  = idatea
  
  ! Determinacao do tipo de produto de umidade
  do i=256,1,-1
   if(usdata_in(I:I) == '/') then
    ipref_start=i+1
    exit
   endif
  
!   print*,'i=',i
!   print*,usdata(i+1:ifname)
!   stop
!  endif
  enddo


  ! DEFINICAO DA ESPESSURA DAS CAMADAS
  ipref = len_trim(usdata_in)
  pref = usdata_in(ipref_start:ipref)
  !print*,'soil data=', pref

  if (pref == 'SM_v2.') then
    n4us=6                    ! modelo V2 com 6 camadas
    allocate(slz_us(0:n4us))
    slz_us = (/-3.0, -2.0, -1.0, -0.5, -0.25, -0.1, 0. /)    
  elseif ( pref == 'GL_SM.GPCP.' .or.  pref == 'GL_SM.GPNR.') then
    n4us=8                    ! modelo GLSM V2 com 8 camadas
    allocate(slz_us(0:n4us))
    slz_us = (/-4.5, -2.5, -1.75, -1.0, -0.5, -0.25, -0.13, -0.05, 0./)    
  else
     n4us=4                   ! modelo original com 4 camadas
     allocate(slz_us(n4us))
     slz_us = (/ -2.4, -0.4, -0.1, 0. /)  
  endif

  ! COMPOSICAO DO NOME DOS ARQUIVOS DE ENTRADA E SAIDA

  if ((runtype(1:7) == 'history').and.                 &
       ((soil_moist == 'h').or.(soil_moist == 'h').or.  &
        (soil_moist == 'a').or.(soil_moist == 'a'))) then

     dif_time = timstr

     !if (timeunit == 'h') then
        !int_dif_time = floor(dif_time/24.)
     !elseif (timeunit == 'm') then
        !int_dif_time = floor(dif_time/1440.)
     !elseif (timeunit == 's') then
        int_dif_time = floor(dif_time/86400.)
     !endif

     call alt_dia(idatea, imontha, iyeara, int_dif_time,   &
          idate2, imonth2, iyear2)

  else

     int_dif_time = 0

  endif

  if ((soil_moist_fail == 'l').or.(soil_moist_fail == 'l')) then
     ! Looking for until 5 days old file
     da = 5
  else
     da = 1
  endif

  sair = .false.

  do i=1,da

     write(cidate,'(I2.2)') idate2
     write(cimon,'(I2.2)') imonth2
     write(ciyear,'(I4)') iyear2

     write(cgrid,'(I1)') ifm

     ! Calculating the hour of simulation
     if ((itimea >= 0000).and.(itimea < 1200)) then
        hourmin = 0000
        if(pref == 'GL_SM.GPCP.'  .or.  pref == 'GL_SM.GPNR.')        hourmin = 00
     else
        hourmin = 1200
        if(pref == 'GL_SM.GPCP.'  .or.  pref == 'GL_SM.GPNR.')        hourmin = 12
     endif
    
     if(pref == 'GL_SM.GPCP.'  .or.  pref == 'GL_SM.GPNR.')       then
       write(cihourmin,'(I2.2)') hourmin
       icihourmin=2
     else
       write(cihourmin,'(I4.4)') hourmin
       icihourmin=4
     endif
     !ipref=INDEX(usdata_in,' ')
     !pref=usdata_in(ipref-2:ipref-1)
     !ipref = len_trim(usdata_in)
    
!     if (usdata_in(ipref_start:ipref) == 'us'.OR.usdata_in(ipref_start:ipref) == 'SM') then
!        pref = usdata_in(ipref_start:ipref)
!     endif
!     if (usdata_in(ipref_start:ipref) == 'SM_v2.' .or. usdata_in(ipref_start:ipref) == 'GL_SM.GPCP.') then
!        pref = usdata_in(ipref_start:ipref)
!     endif
    
     if (pref == 'us') then
       cihourmin=''
       icihourmin=0
     endif

 
     usdata=trim(usdata_in)//ciyear//cimon//cidate//cihourmin(1:icihourmin)//'.vfm'

     ifname=len_trim(usdata)
     inquire(file=usdata(1:ifname),exist=theref)
    
     if (.not.theref) usdata=trim(usdata_in)//ciyear//cimon//cidate//cihourmin(1:icihourmin)//'.gra'

     usmodel=trim(usmodel_in)//ciyear//cimon//cidate//cihourmin(1:icihourmin)//'_g'//cgrid//'.mod'
  
     print *, 'Looking up Soil Moisture Date Files: '
     print *, '  usdata : ', usdata(1:len_trim(usdata))
     print *, '  usmodel: ', usmodel(1:len_trim(usmodel))

     ifname=len_trim(usmodel)
     inquire(file=usmodel(1:ifname),exist=there)
     if(.not.there) then

        ifname=len_trim(usdata)
        inquire(file=usdata(1:ifname),exist=there)
        if (there) then
           sair = .true.
        else
           print *, '  Files not Found!'
        endif
     else
        sair = .true.
     endif

     if (sair) exit

     call alt_dia(idatea, imontha, iyeara, (INT_DIF_TIME - I),   &
          idate2, imonth2, iyear2)

  enddo

  ifname=len_trim(usmodel)
  inquire(file=usmodel(1:ifname),exist=there)
  if(.not.there) then

     ifname=len_trim(usdata)
     inquire(file=usdata(1:ifname),exist=there)
     if (.not.there) then
        print *, '  usdata : ', usdata(1:ifname)
        print *, '  usmodel: ', usmodel(1:len_trim(usmodel))
         if ((SOIL_MOIST_FAIL == 's').or.(SOIL_MOIST_FAIL == 'S')) then
           print *, '* Heterogenous Soil Moisture Init. ERROR!'
           print *, '  Program will Stop!'
           stop
        else
           print *, '  Homogeneous Soil Moisture Initialization.'
           return
        endif
     endif

     print*,'----------------------------------------------'
     print*,'  Homogenous Soil Moisture initialization in'
     print*,'     point outside the South America'
     print*,'----------------------------------------------'

     c1 = .5 * cpi

     do j = 1,n3
        do i = 1,n2

           k2=nint(flpw(i,j))
           pis = c1 * (pi0(k2-1,i,j) + pi0(k2,i,j)   &
                +  pp(k2-1,i,j) + pp(k2,i,j))
           airtemp = theta(k2,i,j) * pis

           do ipat= 2,npat
              do k = 1,mzg

                 nsoil = nint(soil_text(k,i,j,ipat))
                 !print*, 'nsoil',nsoil,slmsts(nsoil)
                 soil_water(k,i,j,ipat) = max(soilcp(nsoil),  &
                      slmstr(k)*slmsts(nsoil))
                 !              print*,soil_water(k,i,j,ipat)

                 soil_energy(k,i,j,ipat) = (airtemp + stgoff(k))   &
                      * (slcpd(nsoil) + soil_water(k,i,j,ipat)  * cliq1000)  &
                      + soil_water(k,i,j,ipat)  * alli1000
              enddo
           enddo
        enddo
     enddo


     !--------------------------------------------------------- HETEROGENEOUS
     print*,'----------------------------------------------'
     print*,'       Soil Moisture Initialization'
     print*,'----------------------------------------------'
      
     ! DADOS DA GRADE DE PRECIPITACAO

     inquire (file=TRIM(usdata_in)//'_ENT', exist=general)
     if (general) then
       open (unit=93,file=TRIM(usdata_in)//'_ENT',status='old')
       read (unit=93,nml=gradeumso)
       close (unit=93, status='keep')      
     elseif (pref == 'us') then    
       latni=-45            
       latnf=12.9477        
       lonni=-82            
       lonnf=-30.055        
       ilatn=0.0359477      
       ilonn=0.0382513      
       nlat=1613
       nlon=1359
     elseif (pref == 'SM') then
       latni=-50.125    
       latnf=40.125      
       lonni=-120.125    
       lonnf=60.125      
       ilatn=0.250      
       ilonn=0.250      
       nlat=362
       nlon=722
     elseif (pref == 'SM_v2.') then  
       latni=-59.875              !SM v2  TRMM global
       latnf=59.875        
       lonni=-179.875        
       lonnf=179.875        
       ilatn=0.250        
       ilonn=0.250        
       nlat=480
       nlon=1440
     elseif (pref == 'GL_SM.GPCP.') then  
       latni=-89.5            !  GPCP global
       latnf=89.5          
       lonni=-179.5          
       lonnf=179.5          
       ilatn=1.            
       ilonn=1.            
       nlat=180
       nlon=360
      
     elseif (pref == 'GL_SM.GPNR.') then
       latni=-89.875             !  TRMM/NAVY + GPCP global
       latnf=89.875        
       lonni=-179.875        
       lonnf=179.875        
       ilatn=0.250        
       ilonn=0.250        
       nlat=720
       nlon=1440
    
    
     else
       print *, 'Unexpected prefix (',pref,') for soil moisture'
       stop 'Program will STOP!'
     endif

     allocate(prlat(nlon,nlat),prlon(nlon,nlat))
     call api_prlatlon(nlon,nlat,prlat,prlon,ilatn,ilonn,latni,lonni)

     allocate(api_us(nlon,nlat,n4us),usdum(n4us))

     print*,'--------------------------------------',' Grid=',ifm
     print*,'Opening soil moisture data= ',usdata
    
     if (.not.theref) then ! arquivo .gra
    
        if(pref == 'us'.OR.pref == 'SM') then   ! arquivo .gra acesso direto
           open(2,status='OLD',form='unformatted',access='direct', &
              recl=4*nlat*nlon*n4us,file=usdata)
           read(UNIT=2,REC=1) api_us        ! water content
           close(2)
          
        call swap32(api_us,nlat*nlon*n4us) ! Verify before call swap32 - Demerval
        if ((api_us(nlat/2,nlon/2,1).lt.0).or.(api_us(nlat/2,nlon/2,1).gt.1))   &
           call swap32(api_us,nlat*nlon*n4us)
        print*,'--------------------------------------'      
                          
        else     ! arquivo .gra acesso sequencial
           open(2,status='OLD',form='unformatted',file=usdata)  
           do k=1,n4us
                 read(2) ((api_us(i,j,k),i=1,nlon),j=1,nlat) ! wetness
           enddo
           close(2)
        endif
        
     else  ! arquivo .vfm
        
        OPEN(UNIT=2,FILE=usdata,FORM='formatted',STATUS='old')
        CALL vfirec(2,api_us,nlat*nlon*n4us,'LIN')
        close(2)
        
          
     endif
     ! geva 20.10.04
    
    
    
    
     ! loop no dominio do modelo
     do j = 1,n3
        do i = 1,n2

           ! evite pontos fora do dominio da grade de umidade de solo
           if(glat(i,j) .lt. latni .or. glat(i,j) .gt. latnf .or. &
                glon(i,j) .lt. lonni .or. glon(i,j) .gt. lonnf) cycle !GOTO 1111
           !print*,'ij=',i,j
           call interpolacao (glon(i,j),glat(i,j),nlon,nlat,prlat,prlon,  &
                i1,i2,ic,j1,j2,jc)

           if(ic.ge.0 .and. jc .ge. 0) then
              !print*,ic,jc,i,j,ifm
              dlonr=0.5*(glon(n2,j)-glon(1,j))/float(n2-1)
              dlatr=0.5*(glat(i,n3)-glat(i,1))/float(n3-1)
              qi1=int(dlonr/ilonn+0.5)
              qi2=int(dlonr/ilonn+0.5)
              qj1=int(dlatr/ilatn+0.5)
              qj2=int(dlatr/ilatn+0.5)
        

              do k=1,n4us
                 ncount = 0
                 usdum(k)=0.

                 do jj =max(1,jc-qj1),min(nlat,jc+qj2)
                    do ii = max(1,ic-qi1),min(nlon,ic+qi2)
                    
                                
                       if (api_us(ii,jj,k).gt.1.E-5) then
                          do ipat=2,npat
                             ncount = ncount + 1                          
                             if (pref == 'us'.OR.pref == 'SM') then
                                usdum(k) = usdum(k) + api_us(ii,jj,k)   ! umidade lida em m3/m3 (vol.) - v1 (us e SM)
                             endif
                             if (pref == 'SM_v2.'.or. pref == 'GL_SM.GPCP.' .or. pref == 'GL_SM.GPNR.') then
                                nsoil = nint(soil_text(k,i,j,ipat))    
                                usdum(k) = usdum(k) + api_us(ii,jj,k)*slmsts(nsoil)  ! umidade lida em % (armazenamento) - v2 (SM_v2.)
                            
                             endif
                          
                          enddo
                       endif
                    enddo
                 enddo
                 usdum(k) = usdum(k) / (float(ncount) + 1.E-10)

              enddo

              do k = mzg,1,-1
                 do kk = n4us,1,-1
                    if (slz(k).ge.slz_us(kk)) then
                       do ipat=2,npat
                          nsoil = nint(soil_text(k,i,j,ipat))
                          soil_water(k,i,j,ipat) = usdum(kk+1)
                          
                          if(usdum(kk+1) .lt. 1.e-5)   &
                               soil_water(k,i,j,ipat) = slmstr(k)*slmsts(nsoil)!oceano

                          soil_water(k,i,j,ipat) = max(soilcp(nsoil), &
                                                    min(soil_water(k,i,j,ipat), slmsts(nsoil)))

                       enddo
                       goto 222
                    elseif (slz(k).lt.slz_us(1)) then
                       do ipat=2,npat
                          nsoil = nint(soil_text(k,i,j,ipat))
                          soil_water(k,i,j,ipat) = usdum(1)
                          
                          if(usdum(1) .lt. 1.e-5)   &
                               soil_water(k,i,j,ipat) = slmstr(k)*slmsts(nsoil)!oceano

                          soil_water(k,i,j,ipat) = max(soilcp(nsoil), &
                                                    min(soil_water(k,i,j,ipat), slmsts(nsoil)))
                       enddo
                       goto 222
                    endif
                 enddo
222              continue
              enddo
           endif
        enddo
     enddo

     deallocate(api_us,usdum,prlat,prlon)
        
     open(2,status='NEW',form='unformatted',access='direct', &
          recl=4*n2*n3*mzg*npat,file=usmodel)
     write(UNIT=2,REC=1) soil_water
     close(2)

  else

     print*,'--------------------------------------', ' Grid=',ifm
     print*,'Opening soil moisture file= ',usmodel
     open(2,status='OLD',form='unformatted',access='direct', &
          recl=4*n2*n3*mzg*npat,file=usmodel)
     read(UNIT=2,REC=1) soil_water
     close(2)
     print*,'--------------------------------------'
  endif


  !----- recalculate soil_energy
  c1 = .5 * cpi

  do j = 1,n3
     do i = 1,n2

        k2=nint(flpw(i,j))
        pis = c1 * (pi0(k2-1,i,j) + pi0(k2,i,j)   &
             +  pp(k2-1,i,j) + pp(k2,i,j))
        airtemp = theta(k2,i,j) * pis

        do ipat= 2,npat
           do k = 1,mzg
              nsoil = nint(soil_text(k,i,j,ipat))

              soil_energy(k,i,j,ipat) = (airtemp + stgoff(k))   &
                   * (slcpd(nsoil) + soil_water(k,i,j,ipat)  * cliq1000)  &
                   + soil_water(k,i,j,ipat)  * alli1000
           enddo
        enddo
     enddo
  enddo

  return
end subroutine soil_moisture_init

!============================================================================
subroutine swap32(a,n)
!============================================================================

  implicit none

  integer,      intent(in)                  :: n
  real(kind=4), intent(inout), dimension(n) :: a

  !
  !      REVERSE ORDER OF BYTES IN INTEGER*4 WORD, or REAL*4
  !
  integer (kind=4)  ::   ijklmn
  !
  character (len=1) :: jtemp(4)
  character (len=1) :: ktemp
  real(kind=4)    :: r4mold=6. !MLO - Real mold for bit transfer, the number itself has no meaning.
  integer(kind=4) :: i4mold=6  !MLO - Integer mold for bit transfer, the number itself has no meaning.
  !
  ! Local variables
  integer :: i, itemp

  equivalence (jtemp(1),itemp)
  !
  save
  !
  ![MLO - Alternative way to save bit representation and preserve interface between subroutine and call
  do i = 1,n
     ijklmn   = transfer(a(i),i4mold)
     itemp    = ijklmn
     ktemp    = jtemp(4)
     jtemp(4) = jtemp(1)
     jtemp(1) = ktemp
     ktemp    = jtemp(3)
     jtemp(3) = jtemp(2)
     jtemp(2) = ktemp
     ijklmn   = itemp
     a(i)     = transfer(ijklmn,r4mold)
  enddo

  return
end subroutine swap32

!
! prlatlon
!----------------------------------------------------------------
! SUB-ROTINA QUE ESTABELECE LATITUDES E LONGITUDES DOS PONTOS DE  
! GRADE DO CAMPO DE PRECIPITACAO
subroutine api_prlatlon(nlon,nlat,prlat,prlon,ilatn,ilonn,latni,lonni)
  implicit none
  real prlat(nlon,nlat),prlon(nlon,nlat)
  real ilatn,ilonn,latni,lonni
  integer nlon,nlat,i,j
  do j=1,nlat
     do i=1,nlon
        prlon(i,j)=lonni+(i-1)*ilonn
        prlat(i,j)=latni+(j-1)*ilatn
     enddo
  enddo
  return
end subroutine api_prlatlon

!----------------------------------------------------------------
! interpolacao
!----------------------------------------------------------------
! SUB-ROTINA QUE REALIZA INTERPOLACAO ENTRE GRADES (RAMS E UMIDADE DO SOLO)  
subroutine interpolacao (glon,glat,nlon,nlat,prlat,prlon,i1,i2,ic,j1,j2,jc)

  implicit none

  integer :: i1, i2, ic, j1, j2, jc

  real    :: glat, glon

  integer :: nlon, nlat

  real    :: prlat(nlon,nlat), prlon(nlon,nlat)

  ! Local Variables
  real    :: diffx1, diffx2, diffy1, diffy2
  integer :: i, j

  do i=1,nlon
     if (glon.le.prlon(i,1)) exit !GOTO 333
  enddo
!333 CONTINUE
  i2=i
  i1=i-1

  do j=1,nlat
     if (glat.le.prlat(1,j)) exit !GOTO 555
  enddo
!555 CONTINUE
  j2=j
  j1=j-1
              
  diffx1 =    glon - prlon(i1,j1)
  diffx2 = -( glon - prlon(i1,j2) )
  diffy1 =    glat - prlat(i1,j1)
  diffy2=  -( glat - prlat(i2,j1) )

  jc=j1
  ic=i1
  if(diffx1.gt.diffx2) ic=i2
  if(diffy1.gt.diffy2) jc=j2

  if(i1 .lt. 1 .or. i1 .gt. nlon .or. j1 .lt. 1 .or. j1 .gt. nlat) then
     ic=-9999
     jc=-9999
  endif

  return
end subroutine interpolacao
!----------------------------------------------------------------
