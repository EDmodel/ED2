!==========================================================================================!
!==========================================================================================!
!      Soil Moisture Estimate for NWP Models                                               !
!      Coded and implemented by Rodrigo Gevaerd and Saulo Freitas                          !
!                                                                                          !
! Ref.: Gevaerd, R. e S. R. Freitas, Estimativa operacional da umidade do solo para        !
!          iniciacao de modelos de previsao numerica da atmosfera.  Parte I: Descricao da  !
!          metodologia e validacao. Rev. Bras. Meteo., volume especial do LBA, 2007.       !
!------------------------------------------------------------------------------------------!
subroutine soil_moisture_init(n1,n2,n3,mzg,npat,ifm,can_theta,can_prss,glat,glon           &
                             ,soil_water,soil_energy,soil_text)
   use mem_grid          , only : runtype         & ! intent(in)
                                , iyeara          & ! intent(in)
                                , imontha         & ! intent(in)
                                , idatea          & ! intent(in)
                                , itimea          ! ! intent(in)
   use mem_soil_moisture , only : soil_moist      & ! intent(in)
                                , soil_moist_fail & ! intent(in)
                                , usdata_in       & ! intent(in)
                                , usmodel_in      & ! intent(in)
                                , oslmsts         & ! intent(in)
                                , osoilcp         ! ! intent(in)
   use io_params         , only : timstr          ! ! intent(in)
   use rconstants        , only : cpi             & ! intent(in)
                                , rocp            & ! intent(in)
                                , p00             & ! intent(in)
                                , p00i            & ! intent(in)
                                , tsupercool      & ! intent(in)
                                , cicevlme        & ! intent(in)
                                , cliqvlme        & ! intent(in)
                                , t00             & ! intent(in)
                                , t3ple           & ! intent(in)
                                , day_sec         ! ! intent(in)
   use leaf_coms         , only : soilcp          & ! intent(in)
                                , slmsts          & ! intent(in)
                                , slcpd           ! ! intent(in)
   use mem_leaf          , only : stgoff          & ! intent(in)
                                , slmstr          & ! intent(in)
                                , slz             ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                           , intent(in)    :: n1
   integer                           , intent(in)    :: n2
   integer                           , intent(in)    :: n3
   integer                           , intent(in)    :: mzg
   integer                           , intent(in)    :: npat
   integer                           , intent(in)    :: ifm
   real   , dimension(    n2,n3,npat), intent(in)    :: can_theta
   real   , dimension(    n2,n3,npat), intent(in)    :: can_prss
   real   , dimension(n2,n3)         , intent(in)    :: glat,glon
   real   , dimension(mzg,n2,n3,npat), intent(inout) :: soil_water
   real   , dimension(mzg,n2,n3,npat), intent(inout) :: soil_energy
   real   , dimension(mzg,n2,n3,npat), intent(inout) :: soil_text
   !----- Local variables. ----------------------------------------------------------------!
   character (len=256)                               :: usdata, usmodel
   character (len=20)                                :: pref
   character (len=2)                                 :: cidate,cimon
   character (len=1)                                 :: cgrid
   character (len=4)                                 :: ciyear
   character (len=4)                                 :: cihourmin
   integer                                           :: i,j,k,ipat,nveg,nsoil
   integer                                           :: qi1,qi2,qj1,qj2,ncount
   integer                                           :: ii,jj,jc,ic,i1,j1,i2,j2,kk
   integer                                           :: ipref,ipref_start
   integer                                           :: icihourmin
   integer                                           :: n4us
   integer                                           :: nlat, nlon
   integer                                           :: int_dif_time,da
   integer                                           :: idate2,imonth2,iyear2,hourmin
   logical                                           :: there,theref,sair,general
   real(kind=8)                                      :: dif_time
   real, dimension(    :)            , allocatable   :: slz_us,usdum
   real, dimension(:,:,:)            , allocatable   :: api_us
   real, dimension(  :,:)            , allocatable   :: prlat,prlon
   real                                              :: can_temp,tsoil
   real                                              :: latni,latnf,lonni,lonnf
   real                                              :: ilatn,ilonn,ilats,ilons
   real                                              :: latn,lonn,lats,lons
   real                                              :: dlatr,dlonr
   real                                              :: slmrel
   real                                              :: swat_new
   integer                                           :: size_usmodel
   integer                                           :: size_expected
   !----- External functions. -------------------------------------------------------------!
   integer                           , external      :: filesize4
   !----- Namelist in case soil moisture is not a default one. ----------------------------!
   namelist /gradeumso/ latni, latnf, lonni, lonnf, ilatn, ilonn, nlat, nlon
   !---------------------------------------------------------------------------------------!

   iyear2  = iyeara
   imonth2 = imontha
   idate2  = idatea
   
   !----- Determining which kind of soil moisture we are using. ---------------------------!
   ipref_start = index(usdata_in,'/',back=.true.) + 1


   !----- Defining the layer thickness based on the dataset. ------------------------------!
   ipref = len_trim(usdata_in)
   pref  = usdata_in(ipref_start:ipref)

   select case (trim(pref))
   case ('SM_v2.')  !----- V2 model, with six layers. -------------------------------------!
     n4us   = 6
     allocate(slz_us(0:n4us))
     slz_us = (/-3.0, -2.0, -1.0, -0.5, -0.25, -0.1, 0. /)
   case ('GL_SM.GPCP.','GL_SM.GPNR.') !----- GLSM model, with eight layers. ---------------!
     n4us   = 8
     allocate(slz_us(0:n4us))
     slz_us = (/-4.5, -2.5, -1.75, -1.0, -0.5, -0.25, -0.13, -0.05, 0./)
   case default !----- Original model, with four layers. ----------------------------------!
      n4us   = 4
      allocate(slz_us(n4us))
      slz_us = (/ -2.4, -0.4, -0.1, 0. /)
   end select


   !----- Making the input/output file name. ----------------------------------------------!
   if ((runtype(1:7) == 'history') .and.                                                   &
       ((soil_moist == 'h') .or. (soil_moist == 'h') .or.                                  &
        (soil_moist == 'a') .or. (soil_moist == 'a'))     ) then

      dif_time = timstr
      int_dif_time = floor(dif_time/day_sec)
      call alt_dia(idatea,imontha,iyeara,int_dif_time,idate2,imonth2,iyear2)
   else
      int_dif_time = 0
   end if

   !----- Check whether it should look for files up to 5 days older than the initial date. !
   if ((soil_moist_fail == 'l')) then
      da = 5
   else
      da = 1
   end if

   sair = .false.

   filefinder: do i=1,da

      write(cidate,fmt='(i2.2)') idate2
      write(cimon ,fmt='(i2.2)') imonth2
      write(ciyear,fmt='(i4.4)') iyear2
      write(cgrid ,fmt='(i1)'  ) ifm

      !----- Finding the hour of simulation. ----------------------------------------------!
      if ((itimea >= 0000).and.(itimea < 1200)) then
         hourmin = 0000
         if(pref == 'GL_SM.GPCP.'  .or.  pref == 'GL_SM.GPNR.')        hourmin = 00
      else
         hourmin = 1200
         if(pref == 'GL_SM.GPCP.'  .or.  pref == 'GL_SM.GPNR.')        hourmin = 12
      end if
     
      if(pref == 'GL_SM.GPCP.'  .or.  pref == 'GL_SM.GPNR.')       then
        write(cihourmin,fmt='(i2.2)') hourmin
        icihourmin=2
      else
        write(cihourmin,fmt='(i4.4)') hourmin
        icihourmin=4
      end if
     
      if (pref == 'us') then
        cihourmin  = ''
        icihourmin = 0
      endif

      !----- Bind the information and create the file name. -------------------------------!
      usdata=trim(usdata_in)//ciyear//cimon//cidate//cihourmin(1:icihourmin)//'.vfm'

      inquire(file=trim(usdata),exist=theref)
      if (.not.theref) then
         usdata=trim(usdata_in)//ciyear//cimon//cidate//cihourmin(1:icihourmin)//'.gra'
      end if

      usmodel = trim(usmodel_in)//ciyear//cimon//cidate//cihourmin(1:icihourmin)
      usmodel = trim(usmodel)//'_g'//cgrid//'.mod'
   
      write(unit=*,fmt='(a)')       'Using the following soil Moisture files: '
      write(unit=*,fmt='(a,1x,a)')  '  - USDATA :',trim(usdata)
      write(unit=*,fmt='(a,1x,a)')  '  - USMODEL:',trim(usmodel)

      inquire(file=trim(usmodel),exist=there)
      if (.not.there) then
         inquire(file=trim(usdata),exist=there)
         if (there) then
            sair = .true.
         else
            write(unit=*,fmt=*)  'Files were not found!!!'
         endif
      else
         sair = .true.
      end if
      if (sair) exit filefinder
      call alt_dia(idatea, imontha, iyeara,(int_dif_time-i),idate2, imonth2, iyear2)
   end do filefinder

   inquire (file=trim(usmodel),exist=there)

   size_usmodel  = filesize4(trim(usmodel))
   size_expected = 4*n2*n3*mzg*npat

   if (size_usmodel /= size_expected) then

      inquire(file=trim(usdata),exist=there)
      if (.not.there) then
         write(unit=*,fmt='(a,1x,a)')  '  - USDATA :',trim(usdata)
         write(unit=*,fmt='(a,1x,a)')  '  - USMODEL:',trim(usmodel)
         if ((SOIL_MOIST_FAIL == 's').or.(SOIL_MOIST_FAIL == 'S')) then
            call abort_run('Failed initialising heterogenous Soil Moisture!!!'             &
                          ,'soil_moisture_init','soil_moisture_init.f90')
         else
            write(unit=*,fmt='(a)') ' Failed initialising heterogeneous soil moisture...'
            write(unit=*,fmt='(a)') ' Going for a homogeneous initial state.'
            return
         end if
      end if

      write(unit=*,fmt='(a)') '|---------------------------------------------|'
      write(unit=*,fmt='(a)') '|  Homogenous Soil Moisture initialisation on |'
      write(unit=*,fmt='(a)') '|     points outside the input domain         |'
      write(unit=*,fmt='(a)') '|---------------------------------------------|'

      where (slmstr > 1.0)
         slmstr = 1.0
      end where

      hoploop: do ipat= 2,npat
         hojloop: do j = 1,n3
            hoiloop: do i = 1,n2
               !----- Finding canopy temperature for this patch. --------------------------!
               can_temp = can_theta(i,j,ipat) * (p00i * can_prss(i,j,ipat)) ** rocp

               hokloop: do k = 1,mzg
                  nsoil = nint(soil_text(k,i,j,ipat))
                  soil_water(k,i,j,ipat) = soilcp(nsoil)                                   &
                                         + slmstr(k) * (slmsts(nsoil) - soilcp(nsoil))
                  tsoil = can_temp + stgoff(k)
                  if (tsoil >= t3ple) then
                     soil_energy(k,i,j,ipat) = slcpd(nsoil) * tsoil                        &
                                             + soil_water(k,i,j,ipat) * cliqvlme           &
                                             * (tsoil-tsupercool)
                  else
                     soil_energy(k,i,j,ipat) = tsoil  * ( slcpd(nsoil)                     &
                                                        + soil_water(k,i,j,ipat)*cicevlme)
                  end if
               end do hokloop
            end do hoiloop
         end do hojloop
      end do hoploop


      write(unit=*,fmt='(a)') '|------------------------------------------------|'
      write(unit=*,fmt='(a)') '|  Heterogeneous Soil Moisture initialisation on |'
      write(unit=*,fmt='(a)') '|     points within the input domain             |'
      write(unit=*,fmt='(a)') '|------------------------------------------------|'
       
      !----- Defining the domain boundaries. ----------------------------------------------!
      inquire (file=TRIM(usdata_in)//'_ENT', exist=general)
      if (general) then
        open (unit=93,file=TRIM(usdata_in)//'_ENT',status='old')
        read (unit=93,nml=gradeumso)
        close (unit=93, status='keep')      
      else
         select case (trim(pref))
         case ('us')
           latni = -45.0
           latnf = 12.9477
           lonni = -82.0
           lonnf = -30.055
           ilatn = 0.0359477
           ilonn = 0.0382513
           nlat  = 1613
           nlon  = 1359

         case ('SM')
           latni = -50.125
           latnf = 40.125
           lonni = -120.125
           lonnf = 60.125
           ilatn = 0.250
           ilonn = 0.250
           nlat  = 362
           nlon  = 722

         case ('SM_v2.')
           latni = -59.875
           latnf = 59.875
           lonni = -179.875
           lonnf = 179.875
           ilatn = 0.250
           ilonn = 0.250
           nlat  = 480
           nlon  = 1440

         case ('GL_SM.GPCP.')
           latni = -89.5
           latnf = 89.5
           lonni = -179.5
           lonnf = 179.5
           ilatn = 1.
           ilonn = 1.
           nlat  = 180
           nlon  = 360

         case ('GL_SM.GPNR.')
           latni = -89.875
           latnf = 89.875
           lonni = -179.875
           lonnf = 179.875
           ilatn = 0.250
           ilonn = 0.250
           nlat  = 720
           nlon  = 1440

         case default
           call abort_run ('Unexpected soil moisture input prefix ('//trim(pref)//')!!!'   &
                          ,'soil_moisture_init','soil_moisture_init.f90')
         end select
      end if

      allocate(prlat(nlon,nlat),prlon(nlon,nlat))
      call api_prlatlon(nlon,nlat,prlat,prlon,ilatn,ilonn,latni,lonni)

      allocate(api_us(nlon,nlat,n4us),usdum(n4us))


      write(unit=*,fmt='(a)') '------------------------------------------------'
      write(unit=*,fmt='(a,1x,i2)') ' + Grid:',ifm
      write(unit=*,fmt='(a,1x,a)')  ' + File:',trim(usdata)
      write(unit=*,fmt='(a)') '------------------------------------------------'

      if (.not.theref) then 
         !----- GrADS format (.gra) file. -------------------------------------------------!

         if(pref == 'us'.OR.pref == 'SM') then
            open (unit=19,status='OLD',form='unformatted',access='direct'                  &
                 ,recl=4*nlat*nlon*n4us,file=usdata)
            read (unit=19,rec=1) api_us        ! water content
            close(unit=19,status='keep')

            if ((api_us(nlat/2,nlon/2,1) < 0) .or. (api_us(nlat/2,nlon/2,1) > 1.))         &
               call swap32(api_us,nlat*nlon*n4us)

         else     
            !----- Sequential access file. ------------------------------------------------!
            open(unit=19,status='OLD',form='unformatted',file=usdata)
            do k=1,n4us
               read(unit=19) ((api_us(i,j,k),i=1,nlon),j=1,nlat) ! wetness
            end do
            close(unit=19,status='keep')
         end if
         
      else  
         !---- VFM file. ------------------------------------------------------------------!
         open (unit=19,file=usdata,form='formatted',status='old')
         CALL vfirec(19,api_us,nlat*nlon*n4us,'LIN')
         close(unit=19,status='keep')
      end if

      !---- Domain loop. ------------------------------------------------------------------!
      readjloop: do j = 1,n3
         readiloop: do i = 1,n2

            !----- Check whether the point falls within the domain.  If not, skip it. -----!
            if(glat(i,j) < latni .or. glat(i,j) > latnf .or.                               &
               glon(i,j) < lonni .or. glon(i,j) > lonnf) cycle readiloop

            call nearest_neighbour_sm(glon(i,j),glat(i,j),nlon,nlat,prlat,prlon            &
                                     ,i1,i2,ic,j1,j2,jc)

            if(ic >= 0 .and. jc >= 0) then
               dlonr = 0.5 * (glon(n2,j)-glon(1,j)) / real(n2-1)
               dlatr = 0.5 * (glat(i,n3)-glat(i,1)) / real(n3-1)
               qi1 = int(dlonr / ilonn + 0.5)
               qi2 = int(dlonr / ilonn + 0.5)
               qj1 = int(dlatr / ilatn + 0.5)
               qj2 = int(dlatr / ilatn + 0.5)

               do k=1,n4us
                  ncount = 0
                  usdum(k)=0.

                  do jj =max(1,jc-qj1),min(nlat,jc+qj2)
                     do ii = max(1,ic-qi1),min(nlon,ic+qi2)

                        if (api_us(ii,jj,k) > 1.e-5) then
                           do ipat=2,npat
                              ncount = ncount + 1
                              !------------------------------------------------------------!
                              !    Deciding in which units the soil moisture data are.     !
                              !------------------------------------------------------------!
                              select case (trim(pref))
                              case ('us','SM')
                                 !---------------------------------------------------------!
                                 !    Soil moisture in [m3/m3].  Since the current version !
                                 ! of BRAMS-4.0.6 (and ED2-BRAMS) utilises a different     !
                                 ! number, we first normalise the soil moisture then scale !
                                 ! with the current parameter.                             !
                                 !---------------------------------------------------------!
                                 nsoil    = nint(soil_text(k,i,j,ipat))
                                 slmrel   = (api_us(ii,jj,k) - osoilcp(nsoil))             &
                                          / (oslmsts(nsoil)  - osoilcp(nsoil))
                                 swat_new = soilcp(nsoil)                                  &
                                          + slmrel * (slmsts(nsoil)  - soilcp(nsoil))
                                 usdum(k) = usdum(k) + swat_new

                              case ('SM_v2.','GL_SM.GPCP.','GL_SM.GPNR.')
                                 !---------------------------------------------------------!
                                 !     Soil moisture in fraction of saturation.  Because   !
                                 ! the relative position of dry air soil has changed, we   !
                                 ! must scale the range to fall within the [soilcp;slmsts] !
                                 ! interval.                                               !
                                 !---------------------------------------------------------!
                                 nsoil    = nint(soil_text(k,i,j,ipat))
                                 slmrel   = ( api_us(ii,jj,k) * oslmsts(nsoil)             &
                                            - osoilcp(nsoil))                              &
                                          / (oslmsts(nsoil)  - osoilcp(nsoil))
                                 swat_new = soilcp(nsoil)                                  &
                                          + slmrel * (slmsts(nsoil)  - soilcp(nsoil))
                                 usdum(k) = usdum(k) + swat_new
                              end select
                           end do
                        end if
                     end do
                  end do
                  usdum(k) = usdum(k) / (real(ncount) + 1.E-10)
               end do

               kloop: do k = mzg,1,-1
                  kkloop: do kk = n4us,1,-1
                     if (slz(k) >= slz_us(kk)) then
                        do ipat=2,npat
                           nsoil = nint(soil_text(k,i,j,ipat))
                           !----- Only reasonable soil moisture values are accepted. ------!
                           if (usdum(kk+1) < soilcp(nsoil)) then
                              soil_water(k,i,j,ipat) = soilcp(nsoil)
                           elseif (usdum(kk+1) > slmsts(nsoil)) then
                              soil_water(k,i,j,ipat) = slmsts(nsoil)
                           else
                              soil_water(k,i,j,ipat) = usdum(kk+1)
                           end if
                        end do
                        cycle kloop
                     else
                        do ipat=2,npat
                           nsoil = nint(soil_text(k,i,j,ipat))
                           if (usdum(1) < soilcp(nsoil)) then
                              soil_water(k,i,j,ipat) = soilcp(nsoil)
                           elseif (usdum(1) > slmsts(nsoil)) then
                              soil_water(k,i,j,ipat) = slmsts(nsoil)
                           else
                              soil_water(k,i,j,ipat) = usdum(1)
                           end if
                        end do
                        cycle kloop
                     end if
                  end do kkloop
               end do kloop
            end if
         end do readiloop
      end do readjloop

      deallocate(api_us,usdum,prlat,prlon)
      
      !----- Write the soil moisture into the output. -------------------------------------!
      open (unit=19,file=usmodel,status='replace',form='unformatted',access='direct'       &
           ,recl=4*n2*n3*mzg*npat)
      write(unit=19,rec=1) soil_water
      close(unit=19,status='keep')

   else

      write(unit=*,fmt='(a)') '------------------------------------------------'
      write(unit=*,fmt='(a,1x,i2)') ' + Grid:',ifm
      write(unit=*,fmt='(a,1x,a)')  ' + File:',trim(usmodel)
      write(unit=*,fmt='(a)') '------------------------------------------------'

      open (unit=19,file=usmodel,status='old',form='unformatted',access='direct'           &
           ,recl=4*n2*n3*mzg*npat)
      read (unit=19,rec=1) soil_water
      close(unit=19,status='keep')
   end if


   heploop: do ipat= 2,npat
      hejloop: do j = 1,n3
        heiloop: do i = 1,n2

            !----- Find the canopy temperature for this patch. ----------------------------!
            can_temp = can_theta(i,j,ipat) * (p00i * can_prss(i,j,ipat)) ** rocp

            hekloop: do k = 1,mzg
               nsoil = nint(soil_text(k,i,j,ipat))
               !----- Make sure that the soil moisture is bounded... ----------------------!
               soil_water(k,i,j,npat) = max(soilcp(nsoil)                                  &
                                           ,min(soil_water(k,i,j,ipat), slmsts(nsoil)))
               tsoil = can_temp + stgoff(k)

               if (tsoil >= t3ple) then
                  soil_energy(k,i,j,ipat) = slcpd(nsoil) * tsoil                           &
                                          + soil_water(k,i,j,ipat) * cliqvlme              &
                                          * (tsoil-tsupercool)
               else
                  soil_energy(k,i,j,ipat) = tsoil  * ( slcpd(nsoil)                        &
                                                     + soil_water(k,i,j,ipat) * cicevlme)
               end if
            end do hekloop
         end do heiloop
      end do hejloop
   end do heploop

   return
end subroutine soil_moisture_init
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This sub-routine swaps the order of the bytes for real numbers of size 4.            !
!------------------------------------------------------------------------------------------!
subroutine swap32(a,n)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                       , intent(in)    :: n
   real(kind=4)    , dimension(n), intent(inout) :: a
   !----- Local variables. ----------------------------------------------------------------!
   integer (kind=4)                              :: ijklmn
   character(len=1), dimension(4)                :: jtemp
   character(len=1)                              :: ktemp
   real(kind=4)                                  :: r4mould
   integer(kind=4)                               :: i4mould
   integer                                       :: i
   integer                                       :: itemp
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     This will ensure that the bit values in the character and integer variables       !
   ! match.                                                                                !
   !---------------------------------------------------------------------------------------!
   equivalence (jtemp(1),itemp)
   save


   !---------------------------------------------------------------------------------------!
   !     These values themselves are meaningless, they are just moulds that are needed for !
   ! the bit transfer.                                                                     !
   !---------------------------------------------------------------------------------------!
   r4mould = 6.
   i4mould = 6
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !    This method saves the bit representation whilst preserving the interface between   !
   ! this subroutine and its calls throughout the code.                                    !
   !---------------------------------------------------------------------------------------!
   do i = 1,n
      ijklmn   = transfer(a(i),i4mould)
      itemp    = ijklmn
      ktemp    = jtemp(4)
      jtemp(4) = jtemp(1)
      jtemp(1) = ktemp
      ktemp    = jtemp(3)
      jtemp(3) = jtemp(2)
      jtemp(2) = ktemp
      ijklmn   = itemp
      a(i)     = transfer(ijklmn,r4mould)
   end do
   !---------------------------------------------------------------------------------------!

   return
end subroutine swap32
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!      This sub-routine generates the longitudes and latitudes of the soil moisture input  !
! data.                                                                                    !
!------------------------------------------------------------------------------------------!
subroutine api_prlatlon(nlon,nlat,prlat,prlon,dlatn,dlonn,latni,lonni)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                   , intent(in)  :: nlon
   integer                   , intent(in)  :: nlat
   real, dimension(nlon,nlat), intent(out) :: prlon
   real, dimension(nlon,nlat), intent(out) :: prlat
   real                      , intent(in)  :: dlonn
   real                      , intent(in)  :: dlatn
   real                      , intent(in)  :: lonni
   real                      , intent(in)  :: latni
   !----- Local variables. ----------------------------------------------------------------!
   integer                                 :: i
   integer                                 :: j
   !---------------------------------------------------------------------------------------!

   do j=1,nlat
      do i=1,nlon
         prlon(i,j) = lonni + real(i-1) * dlonn
         prlat(i,j) = latni + real(j-1) * dlatn
      end do
   end do
   return
end subroutine api_prlatlon
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!      This sub-routine finds the input data point that is the closest to a given BRAMS    !
! grid point.  This utilises the minimum great circle distance.                            !
!------------------------------------------------------------------------------------------!
subroutine nearest_neighbour_sm(glon,glat,nlon,nlat,prlat,prlon,i1,i2,ic,j1,j2,jc)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                   , intent(in)  :: nlon
   integer                   , intent(in)  :: nlat
   real, dimension(nlon,nlat), intent(in)  :: prlon
   real, dimension(nlon,nlat), intent(in)  :: prlat
   real                      , intent(in)  :: glon
   real                      , intent(in)  :: glat
   integer                   , intent(out) :: i1
   integer                   , intent(out) :: i2
   integer                   , intent(out) :: ic
   integer                   , intent(out) :: j1
   integer                   , intent(out) :: j2
   integer                   , intent(out) :: jc
   !----- Local variables. ----------------------------------------------------------------!
   integer                                 :: i
   integer                                 :: j
   real                                    :: closest
   real                                    :: distnow
   !----- External functions. -------------------------------------------------------------!
   real                      , external    :: dist_gc
   !---------------------------------------------------------------------------------------!

   !----- Find the dataset grid point index that corresponds to the BRAMS longitude. ------!
   lonloop: do i=1,nlon
      if (glon <= prlon(i,1)) exit lonloop
   end do lonloop
   !----- Found the point, 
   i2 = i
   i1 = i-1

   !----- Find the dataset grid point index that corresponds to the BRAMS latitude. -------!
   latloop: do j=1,nlat
      if (glat <= prlat(1,j)) exit latloop
   end do latloop
   j2=j
   j1=j-1


   !---------------------------------------------------------------------------------------!
   !     Flag the grid centre to an absurd number in case the grid point is outside the    !
   ! dataset domain.  If this turns out the be the case, quit without checking the         !
   ! distances.                                                                            !
   !---------------------------------------------------------------------------------------!
   if(i1 < 1 .or. i1 > nlon .or. j1 < 1 .or. j1 > nlat) then
      ic=-9999
      jc=-9999
      return
   end if


   !----- Start the distance of the closest point so it will be switched by the first. ----!
   closest = huge(1.)


   !----- Try the four combinations, and pick up the closest one. -------------------------!
   do j = j1,j2
      do i = i1,i2
         distnow = dist_gc(glon,prlon(i,j),glat,prlat(i,j))

         if (distnow < closest) then
           ic      = i
           jc      = j
           distnow = closest
         end if
      end do
   end do

   return
end subroutine nearest_neighbour_sm
!==========================================================================================!
!==========================================================================================!
