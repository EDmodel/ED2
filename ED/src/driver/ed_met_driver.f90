!==========================================================================================!
!==========================================================================================!
!     This subroutine is the main procedure that will take care of reading the meteoro-    !
! logic data to be used during the run.                                                    !
!------------------------------------------------------------------------------------------!
subroutine read_met_driver_head()
   use ed_max_dims       , only : max_met_vars     ! ! intent(in)
   use met_driver_coms, only : nformats         & ! intent(in)
                             , met_names        & ! intent(out)
                             , met_nlon         & ! intent(out)
                             , met_nlat         & ! intent(out)
                             , met_dx           & ! intent(out)
                             , met_dy           & ! intent(out)
                             , met_xmin         & ! intent(out)
                             , met_ymin         & ! intent(out)
                             , met_nv           & ! intent(out)
                             , met_vars         & ! intent(out)
                             , met_frq          & ! intent(out)
                             , met_interp       & ! intent(out)
                             , ed_met_driver_db & ! intent(out)
                             , no_ll            ! ! intent(out)
   implicit none  
   !----- Local variables -----------------------------------------------------------------!
   logical :: l1
   logical :: yes_lat     ! Logical for determining whether latitude grids are present
   logical :: yes_lon     ! Logical for determining whether longitude grids are present
   integer :: iformat
   integer :: n
   !---------------------------------------------------------------------------------------!


   !----- First thing, let's check whether the meterological data metafile exists. --------!
   inquire(file=trim(ed_met_driver_db),exist=l1)
   if (.not. l1) then
      write (unit=*,fmt='(a)') 'File '//trim(ed_met_driver_db)//' not found!'
      write (unit=*,fmt='(a)') 'Specify ED_MET_DRIVER_DB properly in ED namelist.'
      call fatal_error('Ed_met_driver_db not found!','read_met_driver_head'                &
                      &,'ed_met_driver.f90')
   end if
   !---------------------------------------------------------------------------------------!


   !----- Loading the meterorological data metafile information. --------------------------!
   open(unit=12,file=trim(ed_met_driver_db),form='formatted',status='old')
   read(unit=12,fmt=*)  ! skip header

   !------ Read the number of different file formats. -------------------------------------!
   read(unit=12,fmt=*) nformats
   !------ Allocate the header information for each format --------------------------------!
   allocate(met_names (nformats)              )
   allocate(met_nlon  (nformats)              )
   allocate(met_nlat  (nformats)              )
   allocate(met_dx    (nformats)              )
   allocate(met_dy    (nformats)              )
   allocate(met_xmin  (nformats)              )
   allocate(met_ymin  (nformats)              )
   allocate(met_nv    (nformats)              )
   allocate(met_vars  (nformats, max_met_vars))
   allocate(met_frq   (nformats, max_met_vars))
   allocate(met_interp(nformats, max_met_vars))

   !----- Just to initialize, if lon/lat are both found, it will become .false. -----------!
   no_ll = .true.    
   !----- Read the information for each format. -------------------------------------------!
   do iformat = 1,nformats
      read(unit=12,fmt='(a)')  met_names(iformat)
      read(unit=12,fmt=*)      met_nlon(iformat), met_nlat(iformat), met_dx(iformat)       &
                             , met_dy(iformat)  , met_xmin(iformat), met_ymin(iformat)
      read(unit=12,fmt=*)      met_nv(iformat)
      read(unit=12,fmt=*)      (met_vars(iformat,n)  ,n=1,met_nv(iformat))
      read(unit=12,fmt=*)      (met_frq(iformat,n)   ,n=1,met_nv(iformat))
      read(unit=12,fmt=*)      (met_interp(iformat,n),n=1,met_nv(iformat))
      
      !----- Just making sure that the variable list is case insensitive. -----------------!
      call tolower(met_vars(iformat,1:met_nv(iformat)),met_nv(iformat))
      
      !----- First check - see if lat/lon data are there. ---------------------------------!
      yes_lon = any(met_vars(iformat,1:met_nv(iformat)) == 'lon')
      yes_lat = any(met_vars(iformat,1:met_nv(iformat)) == 'lat')
      

      !----- Check to see if you have both, none or one of each. --------------------------!
      if (yes_lat .and. yes_lon) then
         no_ll = .false.
      elseif (yes_lat .neqv. yes_lon) then
         call fatal_error('You are missing a lat or a lon variable in the met nl'          &
                         ,'read_met_driver_head','ed_met_driver.f90')
      end if
   end do

   close (unit=12,status='keep')

   return
end subroutine read_met_driver_head
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine will initialise the arrays that will temporarily receive the meteoro- !
! logical dataset. It will nullify the pointers, allocate them, then assign a very strange !
! number so in case the variable is not properly initialised the user will notice some-    !
! thing went wrong right away.                                                             !
!------------------------------------------------------------------------------------------!
subroutine init_met_drivers

   use ed_max_dims        , only : max_met_vars
   use met_driver_coms , only : nformats          & ! intent(in)
                              , met_names         & ! intent(out)
                              , met_nlon          & ! intent(out)
                              , met_nlat          & ! intent(out)
                              , met_dx            & ! intent(out)
                              , met_dy            & ! intent(out)
                              , met_xmin          & ! intent(out)
                              , met_ymin          & ! intent(out)
                              , met_nv            & ! intent(out)
                              , met_vars          & ! intent(out)
                              , met_frq           & ! intent(out)
                              , met_interp        & ! intent(out)
                              , ed_met_driver_db  & ! intent(out)
                              , no_ll             & ! intent(out)
                              , have_co2          ! ! intent(out)
   use ed_state_vars   , only : edgrid_g          & ! structure
                              , edtype            & ! structure
                              , polygontype       ! ! structure
   use grid_coms       , only : ngrids            ! ! intent(in)
   use consts_coms     , only : day_sec           ! ! intent(in)
   implicit none
   !----- Local variables -----------------------------------------------------------------!
   type(polygontype), pointer :: cpoly
   type(edtype)     , pointer :: cgrid
   integer                    :: ipy,igr
   integer                    :: iformat
   integer                    :: iv
   integer                    :: mem_size
   real                       :: westedge,eastedge,southedge,northedge
   !---------------------------------------------------------------------------------------!


   !----- Set the lapse rates. ------------------------------------------------------------!
   do igr=1,ngrids
      call setLapseParms(edgrid_g(igr))
   end do

   !----- Read the information for each format. -------------------------------------------!
   have_co2=.false.
   formloop: do iformat = 1,nformats
      !----- Finding the met driver boundaries. -------------------------------------------!
      westedge  = met_xmin(iformat) - 0.5 * met_dx(iformat)
      eastedge  = met_xmin(iformat) + (real(met_nlon(iformat))-0.5) * met_dx(iformat)
      southedge = met_ymin(iformat) - 0.5 * met_dy(iformat)
      northedge = met_ymin(iformat) + (real(met_nlat(iformat))-0.5) * met_dy(iformat)
      gridloop: do igr = 1,ngrids
         cgrid => edgrid_g(igr)

         polyloop: do ipy = 1,cgrid%npolygons

            cpoly => cgrid%polygon(ipy)

            !----- Make sure site falls within file domain. -------------------------------!
            if(no_ll .and. (cgrid%lon(ipy) < westedge .or. cgrid%lat(ipy) < southedge .or. &
                            cgrid%lon(ipy) > eastedge .or. cgrid%lat(ipy) > northedge))    &
            then 
               
               write(unit=*,fmt='(a)') '=================================================='
               write(unit=*,fmt='(a)') ' Polygon lies outside the met driver domain!!!'
               write(unit=*,fmt='(a)') '=================================================='
               write(unit=*,fmt='(a)')           ' + Polygon: '
               write(unit=*,fmt='(a,1x,es12.5)') '   - Longitude : ',cgrid%lon(ipy)
               write(unit=*,fmt='(a,1x,es12.5)') '   - Latitude  : ',cgrid%lat(ipy)
               write(unit=*,fmt='(a,1x,a)')      ' + File: : ',trim(met_names(iformat))
               write(unit=*,fmt='(a)')           ' + Meteorological data boundaries: '
               write(unit=*,fmt='(a,1x,es12.5)') '   - West      : ',westedge
               write(unit=*,fmt='(a,1x,es12.5)') '   - East      : ',eastedge
               write(unit=*,fmt='(a,1x,es12.5)') '   - South     : ',southedge
               write(unit=*,fmt='(a,1x,es12.5)') '   - North     : ',northedge
               call fatal_error('Polygon outside the domain of the meteorological drivers' &
                               ,'init_met_drivers','ed_met_driver.f90')
            end if
            
            !----- Loop over variables. ---------------------------------------------------!
            varloop: do iv = 1,met_nv(iformat)
               
               !---------------------------------------------------------------------------!
               !  Either this is a constant or variable in time.                           !
               !---------------------------------------------------------------------------!
               select case(met_interp(iformat,iv))
               case(0,3)    !----- Variable. ----------------------------------------------!
                  mem_size = nint(day_sec / met_frq(iformat,iv)) * 31
               case(1)     !----- Interpolated.   Read two months in. ---------------------!
                  mem_size = 2 * nint(day_sec / met_frq(iformat,iv)) * 31
               case (2,4)  !----- Constant in time. ---------------------------------------!
                  mem_size = 1
               case default
                  call fatal_error('Invalid met_interp! It should be between 0 and 4.'     &
                                  ,'init_met_drivers','ed_met_driver.f90')
               end select

               !---------------------------------------------------------------------------!
               !     Allocate variables accordingly, and initialise them with a crazy      !
               ! number to make the model crash in case it attempts to use them without    !
               ! assigning an appropriate number.                                          !
               !---------------------------------------------------------------------------!
               select case (trim(met_vars(iformat,iv)))
               case ('nbdsf')   !----- Near IR beam downward shortwave flux. -- [   W/m²] -!
                  nullify(cgrid%metinput(ipy)%nbdsf)
                  allocate(cgrid%metinput(ipy)%nbdsf(mem_size))
                  cgrid%metinput(ipy)%nbdsf = huge(1.)
               case ('nddsf')   !----- Near IR diffuse downward shortwave flux. [   W/m²] -!
                  nullify(cgrid%metinput(ipy)%nddsf)
                  allocate(cgrid%metinput(ipy)%nddsf(mem_size))
                  cgrid%metinput(ipy)%nddsf = huge(1.)
               case ('vbdsf')   !----- Visible beam downward shortwave flux. -- [   W/m²] -!
                  nullify(cgrid%metinput(ipy)%vbdsf)
                  allocate(cgrid%metinput(ipy)%vbdsf(mem_size))
                  cgrid%metinput(ipy)%vbdsf = huge(1.)
               case ('vddsf')   !----- Visible diffuse downward shortwave flux. [   W/m²] -!
                  nullify(cgrid%metinput(ipy)%vddsf)
                  allocate(cgrid%metinput(ipy)%vddsf(mem_size))
                  cgrid%metinput(ipy)%vddsf = huge(1.)
               case ('prate')   !----- Precipitation rate. -------------------- [kg/m²/s] -!
                  nullify(cgrid%metinput(ipy)%prate)
                  allocate(cgrid%metinput(ipy)%prate(mem_size))
                  cgrid%metinput(ipy)%prate = huge(1.)
               case ('dlwrf')   !----- Downwelling longwave radiation. -------- [   W/m²] -!
                  nullify(cgrid%metinput(ipy)%dlwrf)
                  allocate(cgrid%metinput(ipy)%dlwrf(mem_size))
                  cgrid%metinput(ipy)%dlwrf = huge(1.)
               case ('pres')   !----- Air pressure. --------------------------- [   N/m²] -!
                  nullify(cgrid%metinput(ipy)%pres)
                  allocate(cgrid%metinput(ipy)%pres(mem_size))
                  cgrid%metinput(ipy)%pres = huge(1.)
               case ('hgt')     !----- Air pressure. -------------------------- [      m] -!
                  nullify(cgrid%metinput(ipy)%hgt)
                  allocate(cgrid%metinput(ipy)%hgt(mem_size))
                  cgrid%metinput(ipy)%hgt = huge(1.)
               case ('ugrd')    !----- Zonal wind. ---------------------------- [    m/s] -!
                  nullify(cgrid%metinput(ipy)%ugrd)
                  allocate(cgrid%metinput(ipy)%ugrd(mem_size))
                  cgrid%metinput(ipy)%ugrd = huge(1.)
               case ('vgrd')    !----- Meridional wind. ----------------------- [    m/s] -!
                  nullify(cgrid%metinput(ipy)%vgrd)
                  allocate(cgrid%metinput(ipy)%vgrd(mem_size))
                  cgrid%metinput(ipy)%vgrd = huge(1.)
               case ('sh')      !----- Specific humidity. --------------------- [kg/kg_a] -!
                  nullify(cgrid%metinput(ipy)%sh)
                  allocate(cgrid%metinput(ipy)%sh(mem_size))
                  cgrid%metinput(ipy)%sh = huge(1.)
               case ('tmp')     !----- Air temperature. ----------------------- [      K] -!
                  nullify(cgrid%metinput(ipy)%tmp)
                  allocate(cgrid%metinput(ipy)%tmp(mem_size))
                  cgrid%metinput(ipy)%tmp = huge(1.)
               case ('co2')     !----- CO2 mixing ratio. ---------------------- [    ppm] -!
                  have_co2=.true.
                  nullify(cgrid%metinput(ipy)%co2)
                  allocate(cgrid%metinput(ipy)%co2(mem_size))
                  cgrid%metinput(ipy)%co2 = huge(1.)
               case ('lat','lon') !---- Latitude and longitude: skip them. ----------------!
               case default
                  call fatal_error('Invalid met variable'//trim(met_vars(iformat,iv))//'!' &
                                  ,'init_met_drivers'                                &
                                  ,'ed_met_driver.f90')
               end select

            end do varloop
         end do polyloop
      end do gridloop
   end do formloop

   return
end subroutine init_met_drivers
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine will control the HDF5 file reading.                                  !
!------------------------------------------------------------------------------------------!
subroutine read_met_drivers_init

   use ed_state_vars  , only : edgrid_g       & ! structure
                             , edtype         & ! structure
                             , polygontype    ! ! structure
   use hdf5_utils     , only : shdf5_open_f   & ! subroutine
                             , shdf5_close_f  ! ! subroutine
   use met_driver_coms, only : nformats       & ! intent(in)
                             , ishuffle       & ! intent(in)
                             , met_names      & ! intent(in)
                             , met_nv         & ! intent(in)
                             , met_interp     & ! intent(in)
                             , met_frq        & ! intent(in)
                             , metcyc1        & ! intent(in)
                             , metcycf        & ! intent(in)
                             , metyears       ! ! intent(inout)
   use mem_sites      , only : grid_type      ! ! intent(in)
   use ed_misc_coms      , only : current_time   & ! intent(in)
                             , iyeara         & ! intent(in)
                             , iyearz         ! ! intent(in)
   use grid_coms      , only : ngrids         ! ! intent(in)
   use consts_coms    , only : day_sec        ! ! intent(in)
   implicit none
   !----- Local variables -----------------------------------------------------------------!
   type(edtype)      , pointer :: cgrid
   character(len=256)          :: infile
   integer                     :: igr
   integer                     :: year_use
   integer                     :: iformat
   integer                     :: iv
   integer                     :: offset
   integer                     :: m2
   integer                     :: y2
   integer                     :: year_use_2
   integer                     :: nyears
   integer                     :: ncyc
   integer                     :: nfullcyc
   integer                     :: i1stfull
   integer                     :: iyear
   integer, dimension(8)       :: seedtime
   real                        :: runif
   logical                     :: exans
   !----- Local constants -----------------------------------------------------------------!
   character(len=3), dimension(12), parameter :: mname = (/ 'JAN', 'FEB', 'MAR', 'APR'     &
                                                          , 'MAY', 'JUN', 'JUL', 'AUG'     &
                                                          , 'SEP', 'OCT', 'NOV', 'DEC'   /) 
   !----- External functions --------------------------------------------------------------!
   logical         , external  :: isleap
   !----- Variables to be saved. ----------------------------------------------------------!
   logical         , save      :: first_time=.true.
   !---------------------------------------------------------------------------------------!
 


   !---------------------------------------------------------------------------------------!
   !    If this is the first time this subroutine is called, we build the meteorological   !
   ! driver sequence of years.  Often is the case that the user has a meteorological data- !
   ! set that is shorter than the time span he or she is intending to run. In this case we !
   ! need to fill the other years with some data.                                          !
   !---------------------------------------------------------------------------------------!
   if (first_time) then
      first_time = .false.

      !------ Defining the number of years we have meteorological information available. --!
      ncyc   = metcycf - metcyc1 + 1
      nyears = iyearz  - iyeara  + 1

      !------ Allocating the sequence of years. -------------------------------------------!
      allocate(metyears(nyears))      
      
      !----- Initialising the seed for random number generation. --------------------------!
      call random_seed()

      !------------------------------------------------------------------------------------!
      !     Now we must decide how we are going to select the meteorological driver years, !
      ! in particular, the years outside the meteorological driver range.                  !
      !------------------------------------------------------------------------------------!
      select case (ishuffle)
      case (0)
         !---------------------------------------------------------------------------------!
         !     Make a loop over the years, repeating the loop as many times as we need.    !
         ! First we determine how many full cycles fit between iyeara and metcyc1.  Note   !
         ! that this can be negative in case metcyc1 is less than iyeara, but this will    !
         ! work too.
         !---------------------------------------------------------------------------------!
         nfullcyc = floor(real(metcyc1-iyeara)/real(ncyc))
         i1stfull = metcyc1 - nfullcyc * ncyc
         
         year_use   = metcycf - (i1stfull - iyeara)

         do iyear=1,nyears
            if (year_use == metcycf) then
               year_use = metcyc1
            else
               year_use = year_use + 1
            end if

            metyears(iyear) = year_use
         end do

      !------------------------------------------------------------------------------------!
      !    For years outside the met driver range, we pick up a random year.  Because we   !
      ! use the actual sequence of random years given by the code random generator, this   !
      ! sequence will be always the same for a given met driver and first year.            !                !
      !------------------------------------------------------------------------------------!
      case (1)
         do year_use=iyeara,iyearz
            iyear=year_use-iyeara+1
            if (year_use < metcyc1 .or. year_use > metcycf) then
               call random_number(runif)
               metyears(iyear) = metcyc1 + mod(floor(real(ncyc)*runif),ncyc)
            else
               metyears(iyear) = year_use
            end if
         end do

      !------------------------------------------------------------------------------------!
      !    For years outside the met driver range, we pick up a random year. We skip some  !
      ! random numbers using the system clock, so each time we run the model we will have  !
      ! a different sequence of years (similar to the random method used in STILT).        !
      !------------------------------------------------------------------------------------!
      case (2)
         do year_use=iyeara,iyearz
            iyear=year_use-iyeara+1
            if (year_use < metcyc1 .or. year_use > metcycf) then
               call date_and_time(values=seedtime)
               do iv=0,seedtime(8)
                  call random_number(runif)
               end do
               metyears(iyear) = metcyc1 + mod(floor(real(ncyc)*runif),ncyc)
            else
               metyears(iyear) = year_use
            end if
         end do
      end select
   end if
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     We now retrieve the met driver year based on the stored sequence.                 !
   !---------------------------------------------------------------------------------------!
   iyear    = current_time%year-iyeara+1
   year_use = metyears(iyear)


   gridloop: do igr = 1,ngrids

      cgrid => edgrid_g(igr)

      !----- Loop over the different file formats -----------------------------------------!
      formloop: do iformat = 1, nformats
         
         !----- Create the file name and check whether it exists. -------------------------!
         write(infile,fmt='(a,i4.4,a,a)')   trim(met_names(iformat)), year_use             &
                                           ,mname(current_time%month),'.h5'
         inquire(file=trim(infile),exist=exans)
         if (exans) then
            call shdf5_open_f(trim(infile),'R')
         else
            call fatal_error('Cannot open met driver input file '//trim(infile)//'!'       &
                            ,'read_met_drivers_init','ed_met_driver.f90')
         end if
         
         !---------------------------------------------------------------------------------!
         !     The following subroutine determines grid indices of each polygon's match to !
         ! the met data.                                                                   !
         !---------------------------------------------------------------------------------!
         call getll(cgrid,iformat)
         
         !----- Loop over variables. and read the data. -----------------------------------!
         do iv = 1, met_nv(iformat)
            offset = 0
            call read_ol_file(iformat, iv, year_use, mname(current_time%month)          &
                                , current_time%year, offset, cgrid)
         end do

         !-----  Close the HDF5 file. -----------------------------------------------------!
         call shdf5_close_f()
         
         !---------------------------------------------------------------------------------!
         !---------------------------------------------------------------------------------!
         !      For all interpolated variables, we also need the next time.                !
         !---------------------------------------------------------------------------------!
         !------ Find next month and year -------------------------------------------------!
         m2 = current_time%month + 1
         
         !---------------------------------------------------------------------------------!
         !     If this takes us into the next year, take the next year in sequence and     !
         ! reset month to January.                                                         !
         !---------------------------------------------------------------------------------!
         if (m2 == 13) then
            m2 = 1
            y2 = current_time%year + 1
         else 
            !----- Otherwise, use the same year. ------------------------------------------!
            y2 = current_time%year
         end if
         iyear = y2 - iyeara + 1
         year_use_2 = metyears(iyear)
         
         !----- Now, open the file once. --------------------------------------------------!
         write(infile,fmt='(a,i4.4,a,a)')  trim(met_names(iformat)), year_use_2            &
                                          ,mname(m2),'.h5'
         inquire(file=trim(infile),exist=exans)
         if (exans) then
            call shdf5_open_f(trim(infile),'R')
         else
            call fatal_error('Cannot open met driver input file '//trim(infile)//'!' &
                            ,'read_met_drivers_init','ed_met_driver.f90')
         end if
         
         !----- Loop over variables. ------------------------------------------------------!
         varloop: do iv = 1, met_nv(iformat)
            if (met_interp(iformat,iv) == 1) then
               offset = nint(day_sec / met_frq(iformat,iv)) * 31

               !----- Read the file. ------------------------------------------------------!
               call read_ol_file(iformat,iv,year_use_2,mname(m2),y2,offset,cgrid)
            end if
         end do varloop
         
         !----- Close the HDF5 file. ------------------------------------------------------!
         call shdf5_close_f()

      end do formloop
   end do gridloop


   return
end subroutine read_met_drivers_init
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine read_met_drivers

   use ed_state_vars  , only : edgrid_g      & ! structure
                             , edtype        & ! structure
                             , polygontype   ! ! structure
   use mem_sites      , only : grid_type     ! ! structure
   use met_driver_coms, only : nformats      & ! intent(in)
                             , met_names     & ! intent(in)
                             , met_nv        & ! intent(in)
                             , met_interp    & ! intent(in)
                             , met_frq       & ! intent(in)
                             , met_vars      & ! intent(in)
                             , metcyc1       & ! intent(in)
                             , metcycf       ! ! intent(in)
   use ed_misc_coms      , only : current_time  ! ! intent(in)
   use hdf5_utils     , only : shdf5_open_f  & ! subroutine
                             , shdf5_close_f ! ! subroutine
   use ed_state_vars  , only : edtype        & ! structure
                             , edgrid_g      ! ! structure
   use grid_coms      , only : ngrids        ! ! intent(in)
   use consts_coms    , only : day_sec       ! ! intent(in)

   implicit none
   !----- Local variables -----------------------------------------------------------------!
   type(edtype)      , pointer :: cgrid
   character(len=256)          :: infile
   integer                     :: igr
   integer                     :: year_use
   integer                     :: ncyc
   integer                     :: iformat
   integer                     :: iv
   integer                     :: offset
   integer                     :: m2
   integer                     :: y2
   integer                     :: year_use_2
   logical                     :: exans
   !----- Local constants -----------------------------------------------------------------!
   character(len=3), dimension(12), parameter :: mname = (/ 'JAN', 'FEB', 'MAR', 'APR'     &
                                                          , 'MAY', 'JUN', 'JUL', 'AUG'     &
                                                          , 'SEP', 'OCT', 'NOV', 'DEC'   /) 
   !---------------------------------------------------------------------------------------!




   !----- If we need to recycle over years, find the appropriate year to apply. -----------!
   year_use = current_time%year
   ncyc = metcycf - metcyc1 + 1

   !----- If we are after the last year... ------------------------------------------------!
   do while(year_use > metcycf)
      year_use = year_use - ncyc
   end do

   !----- If we are before the first year... ----------------------------------------------!
   do while(year_use < metcyc1)
      year_use = year_use + ncyc
   end do

   gridloop: do igr=1,ngrids

      cgrid => edgrid_g(igr)

      !----- Loop over the different file formats -----------------------------------------!
      formloop: do iformat = 1, nformats
         
         !----- Create the file name and check whether it exists. -------------------------!
         write(infile,'(a,i4.4,a,a)')trim(met_names(iformat)), year_use,   &
              mname(current_time%month),'.h5'

         inquire(file=trim(infile),exist=exans)
         if(exans)then
            call shdf5_open_f(trim(infile),'R')
         else
            call fatal_error('Cannot open met driver input file '//trim(infile)//'!'       &
                            ,'read_met_drivers','ed_met_driver.f90')
         end if
         
         !---------------------------------------------------------------------------------!
         !     The following subroutine determines grid indices of each polygon's match to !
         ! the met data.                                                                   !
         !---------------------------------------------------------------------------------!
         call getll(cgrid,iformat)
         
         !----- Loop over variables. and read the data. -----------------------------------!
         do iv = 1, met_nv(iformat)
            
            !----- Check whether this is an interpolation variable. -----------------------!
            if (met_interp(iformat,iv) /= 1) then
               !----- If not, things are simple.  Just read in the month. -----------------!
               offset = 0
               call read_ol_file(iformat,iv,year_use,mname(current_time%month)          &
                                   ,current_time%year,offset,cgrid)
            else
               !----- Here, just transfer future to current month.  -----------------------!
               call transfer_ol_month(trim(met_vars(iformat,iv)),met_frq(iformat,iv)    &
                                        ,cgrid)
            end if
         end do

         !----- Close the HDF5 file. ------------------------------------------------------!
         call shdf5_close_f()

         !---------------------------------------------------------------------------------!
         !---------------------------------------------------------------------------------!
         !      For all interpolated variables, we also need the next time.                !
         !---------------------------------------------------------------------------------!
         !------ Find next month and year -------------------------------------------------!
         m2 = current_time%month + 1
         y2 = current_time%year
         year_use_2 = year_use
         
         
         !---------------------------------------------------------------------------------!
         !     If this takes us into the next year, increment year and reset month to      !
         ! January.                                                                        !
         !---------------------------------------------------------------------------------!
         if(m2 == 13)then
            m2 = 1
            y2 = current_time%year + 1
            year_use_2 = y2
            
            !----- If we are now after the last year... -----------------------------------!
            do while(year_use_2 > metcycf)
               year_use_2 = year_use_2 - ncyc
            end do
            
            !----- If we are now before the first year... ---------------------------------!
            do while(year_use_2 < metcyc1)
               year_use_2 = year_use_2 + ncyc
            end do
         end if
         
         !----- Now, open the file once. --------------------------------------------------!
         write(infile,fmt='(a,i4.4,a,a)')  trim(met_names(iformat)), year_use_2            &
                                          ,mname(m2),'.h5'
         inquire(file=trim(infile),exist=exans) 
         if (exans) then
            call shdf5_open_f(trim(infile),'R')
         else
            call fatal_error ('Cannot open met driver input file '//trim(infile)//'!'      &
                             ,'read_met_drivers','ed_met_driver.f90')
         end if
      
         !----- Loop over variables. ------------------------------------------------------!
         do iv = 1, met_nv(iformat)
            if(met_interp(iformat,iv) == 1)then
               offset = nint(day_sec / met_frq(iformat,iv)) * 31
               !----- Read the file. ------------------------------------------------------!
               call read_ol_file(iformat,iv,year_use_2,mname(m2),y2,offset,cgrid)
            end if
         end do
         
         !----- Close the HDF5 file. ------------------------------------------------------!
         call shdf5_close_f()
      end do formloop
   end do gridloop


   return
end subroutine read_met_drivers
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine update_met_drivers(cgrid)
  
   use ed_state_vars  , only : edtype        & ! structure
                             , polygontype   ! ! structure
   use met_driver_coms, only : met_frq       & ! intent(in)
                             , nformats      & ! intent(in)
                             , met_nv        & ! intent(in)
                             , met_interp    & ! intent(in)
                             , met_vars      & ! intent(in)
                             , have_co2      & ! intent(in)
                             , initial_co2   ! ! intent(in)
                             , atm_tmp_intercept &
                             , atm_tmp_slope &
                             , prec_intercept &
                             , prec_slope &
                             , humid_scenario
   use ed_misc_coms   , only : current_time  ! ! intent(in)
   use canopy_air_coms, only : ubmin         ! ! intent(in)
   use consts_coms    , only : rdry          & ! intent(in)
                             , cice          & ! intent(in)
                             , cliq          & ! intent(in)
                             , alli          & ! intent(in)
                             , rocp          & ! intent(in)
                             , p00i          & ! intent(in)
                             , cp            & ! intent(in)
                             , day_sec       & ! intent(in)
                             , t00           & ! intent(in)
                             , t3ple         & ! intent(in)
                             , wdnsi         & ! intent(in)
                             , tsupercool    ! ! intent(in)
   use therm_lib      , only : rslif         & ! function
                             , ptqz2enthalpy ! ! function

   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(edtype)     , target  :: cgrid
   !----- Local variables -----------------------------------------------------------------!
   type(polygontype), pointer :: cpoly
   integer                    :: points_per_day
   integer                    :: nday
   integer                    :: np
   integer                    :: ndays_elapsed
   integer                    :: nseconds_elapsed
   integer                    :: mlo
   integer                    :: mhi
   integer                    :: iformat
   integer                    :: iv
   integer                    :: ipy,isi
   real                       :: t1, T0
   real                       :: t2
   real                       :: rvaux, rvsat
   !----- External functions --------------------------------------------------------------!
   logical         , external :: isleap
   !---------------------------------------------------------------------------------------!


   !----- Initialize vels -----------------------------------------------------------------!
   do ipy=1,cgrid%npolygons
      cgrid%met(ipy)%vels = 0.0
   end do
   


   formloop: do iformat = 1, nformats
      
      !----- Loop over variables ----------------------------------------------------------!
      varloop: do iv = 1, met_nv(iformat)
         
         !-----  In these cases there is no time dynamic, simply set mlo to 1 -------------!
         if (met_interp(iformat,iv) /= 1) then
            
            if (met_interp(iformat,iv) == 2 .or. met_interp(iformat,iv) == 4 ) then
               
               mlo = 1
               
            !------------------------------------------------------------------------------!
            !      In these cases there is no interpolation, but there are time varying    !
            ! data.                                                                        !
            !------------------------------------------------------------------------------!
            elseif (met_interp(iformat,iv) == 0 .or. met_interp(iformat,iv) == 3) then
               
               !---------------------------------------------------------------------------!
               !     If this is not an interpolation variable, just find the time point.   !
               !---------------------------------------------------------------------------!
               ndays_elapsed    = current_time%date - 1
               nseconds_elapsed = nint(current_time%time) + ndays_elapsed * nint(day_sec)
               mlo              = int(float(nseconds_elapsed) / met_frq(iformat,iv)) + 1
               
            end if
            
            !----- Find which variable it is, and then fill the sites. --------------------!
            select case (trim(met_vars(iformat,iv)))
            case('nbdsf')   !----- Near IR beam downward shortwave flux. ------ [   W/m²] -!
               do ipy = 1,cgrid%npolygons
                  cgrid%met(ipy)%nir_beam = cgrid%metinput(ipy)%nbdsf(mlo)
               end do

            case('nddsf')   !----- Near IR diffuse downward shortwave flux. --- [   W/m²] -!
               do ipy = 1,cgrid%npolygons
                  cgrid%met(ipy)%nir_diffuse = cgrid%metinput(ipy)%nddsf(mlo)
               end do

            case('vbdsf')   !----- Visible beam downward shortwave flux. ------ [   W/m²] -!
               do ipy = 1,cgrid%npolygons
                  cgrid%met(ipy)%par_beam = cgrid%metinput(ipy)%vbdsf(mlo)
               end do

            case('vddsf')   !----- Visible diffuse downward shortwave flux. --- [   W/m²] -!
               do ipy = 1,cgrid%npolygons
                  cgrid%met(ipy)%par_diffuse = cgrid%metinput(ipy)%vddsf(mlo)
               end do

            case('prate')   !----- Precipitation rate. ------------------------ [kg/m²/s] -!
               do ipy = 1,cgrid%npolygons
                  cgrid%met(ipy)%pcpg = cgrid%metinput(ipy)%prate(mlo)
               end do

            case('dlwrf')   !----- Downwelling longwave radiation. ------------ [   W/m²] -!
               do ipy = 1,cgrid%npolygons
                  cgrid%met(ipy)%rlong = cgrid%metinput(ipy)%dlwrf(mlo)
               end do

            case('pres')    !----- Air pressure. ------------------------------ [   N/m²] -!
               do ipy = 1,cgrid%npolygons
                  cgrid%met(ipy)%prss = cgrid%metinput(ipy)%pres(mlo)
               end do

            case('hgt')     !----- Air pressure. ------------------------------ [      m] -!
               do ipy = 1,cgrid%npolygons
                  cgrid%met(ipy)%geoht = cgrid%metinput(ipy)%hgt(mlo)
               end do

            case('ugrd')    !----- Zonal wind. -------------------------------- [    m/s] -!
               !---------------------------------------------------------------------------!
               !     Here we add the square of the zonal wind, so in the end vels will     !
               ! have (twice) the kinetic energy.  This way if we interpolate in the       !
               ! vertical for the multi-site case, we interpolate energy.                  !
               !---------------------------------------------------------------------------!
               do ipy = 1,cgrid%npolygons
                  cgrid%met(ipy)%vels = cgrid%met(ipy)%vels                                &
                                      + cgrid%metinput(ipy)%ugrd(mlo)**2
               end do

            case('vgrd')    !----- Meridional wind. --------------------------- [    m/s] -!
               !---------------------------------------------------------------------------!
               !     Here we add the square of the meridional wind, so in the end vels     !
               ! will have (twice) the kinetic energy.  This way if we interpolate in the  !
               ! vertical for the multi-site case, we interpolate energy.                  !
               !---------------------------------------------------------------------------!
               do ipy = 1,cgrid%npolygons
                  cgrid%met(ipy)%vels = cgrid%met(ipy)%vels                                &
                                      + cgrid%metinput(ipy)%vgrd(mlo)**2
               end do

            case('sh')      !----- Specific humidity. ------------------------- [kg/kg_a] -!
               !---------------------------------------------------------------------------!
               !     Here we will just get the specific humidity. But soon we will check   !
               ! whether this value is consistent with our thermodynamics (i.e., it will   !
               ! not be supersaturated).  If not we will impose the maximum.               !
               !---------------------------------------------------------------------------!
               do ipy = 1,cgrid%npolygons
                  cgrid%met(ipy)%atm_shv = cgrid%metinput(ipy)%sh(mlo)
               end do

            case('tmp')     !----- Air temperature. --------------------------- [      K] -!
               do ipy = 1,cgrid%npolygons
                  cgrid%met(ipy)%atm_tmp = cgrid%metinput(ipy)%tmp(mlo)
               end do

            case('co2')     !----- CO2 mixing ratio. -------------------------- [    ppm] -!
               do ipy = 1,cgrid%npolygons
                  cgrid%met(ipy)%atm_co2 = cgrid%metinput(ipy)%co2(mlo)
               end do
            end select
            
         !----- In this case, we need to interpolate. -------------------------------------!
         else !if (met_interp(iformat,iv) /= 1) then
            
            ndays_elapsed    = current_time%date - 1
            nseconds_elapsed = nint(current_time%time) + ndays_elapsed * nint(day_sec)
            
            !----- First, get the number of points per day and per month. -----------------!
            points_per_day = nint(day_sec/met_frq(iformat,iv))
            select case (current_time%month)
            case (1,3,5,7,8,10,12)
              nday = 31
            case (4,6,9,11)
              nday = 30

            case (2)
               if (isleap(current_time%year)) then
                  nday = 29
               else
                  nday = 28
               end if
            end select
            np = nday * points_per_day
            
            
            !----- If this is not an interpolation variable, just find the time point. ----!
            ndays_elapsed    = current_time%date - 1
            nseconds_elapsed = nint(current_time%time) + ndays_elapsed * nint(day_sec)
            mlo              = int(float(nseconds_elapsed) / met_frq(iformat,iv)) + 1
            
            
            !----- Get indices ------------------------------------------------------------!
            mlo = int(float(nseconds_elapsed)/met_frq(iformat,iv)) + 1
            mhi = mlo + 1
            
            if(mhi > np)then
               mhi = 1 + nint(day_sec / met_frq(iformat,iv)) * 31
            endif
            
            !----- Get interpolation factors ----------------------------------------------!
            t1 = mod(float(nseconds_elapsed), met_frq(iformat,iv)) /   &
                 met_frq(iformat,iv)
            t2 = 1.0 - t1

            !----- Find the variable and fill the sites. ----------------------------------!
            select case (trim(met_vars(iformat,iv)))
            case('prate')   !----- Precipitation rate. ------------------------ [kg/m²/s] -!
               do ipy = 1,cgrid%npolygons
                  cgrid%met(ipy)%pcpg = cgrid%metinput(ipy)%prate(mhi) * t1                &
                                      + cgrid%metinput(ipy)%prate(mlo) * t2
               end do

            case('dlwrf')   !----- Downwelling longwave radiation. ------------ [   W/m²] -!
               do ipy = 1,cgrid%npolygons
                  cgrid%met(ipy)%rlong = cgrid%metinput(ipy)%dlwrf(mhi) * t1               &
                                       + cgrid%metinput(ipy)%dlwrf(mlo) * t2
               end do

            case('pres')   !----- Air pressure. ------------------------------- [   N/m²] -!
               do ipy = 1,cgrid%npolygons
                  cgrid%met(ipy)%prss = cgrid%metinput(ipy)%pres(mhi) * t1                 &
                                      + cgrid%metinput(ipy)%pres(mlo) * t2
               end do

            case('hgt')     !----- Air pressure. ------------------------------ [      m] -!
               do ipy = 1,cgrid%npolygons
                  cgrid%met(ipy)%geoht = cgrid%metinput(ipy)%hgt(mhi) * t1                 &
                                       + cgrid%metinput(ipy)%hgt(mlo) * t2
               end do

            case('ugrd')    !----- Zonal wind. -------------------------------- [    m/s] -!
               !---------------------------------------------------------------------------!
               !     Here we add the square of the zonal wind, so in the end vels will     !
               ! have (twice) the kinetic energy.  This way if we interpolate in the       !
               ! vertical for the multi-site case, we interpolate energy.                  !
               !---------------------------------------------------------------------------!
               do ipy = 1,cgrid%npolygons
                  cgrid%met(ipy)%vels =  cgrid%met(ipy)%vels                               &
                                      +  ( cgrid%metinput(ipy)%ugrd(mhi) * t1              &
                                         + cgrid%metinput(ipy)%ugrd(mlo) * t2 )**2
               end do

            case('vgrd')    !----- Meridional wind. --------------------------- [    m/s] -!
               !---------------------------------------------------------------------------!
               !     Here we add the square of the meridional wind, so in the end vels     !
               ! will have (twice) the kinetic energy.  This way if we interpolate in the  !
               ! vertical for the multi-site case, we interpolate energy.                  !
               !---------------------------------------------------------------------------!
               do ipy = 1,cgrid%npolygons
                  cgrid%met(ipy)%vels = cgrid%met(ipy)%vels                                &
                                      + ( cgrid%metinput(ipy)%vgrd(mhi) * t1               &
                                        + cgrid%metinput(ipy)%vgrd(mlo) * t2 )**2
               enddo
            case('sh')      !----- Specific humidity. ------------------------- [kg/kg_a] -!
               do ipy = 1,cgrid%npolygons
                  cgrid%met(ipy)%atm_shv = cgrid%metinput(ipy)%sh(mhi) * t1                &
                                         + cgrid%metinput(ipy)%sh(mlo) * t2
               end do

            case('tmp')     !----- Air temperature. --------------------------- [      K] -!
               do ipy = 1,cgrid%npolygons
                  cgrid%met(ipy)%atm_tmp = cgrid%metinput(ipy)%tmp(mhi) * t1               &
                                         + cgrid%metinput(ipy)%tmp(mlo) * t2
               end do

            case('nbdsf')   !----- Near IR beam downward shortwave flux. ------ [   W/m²] -!
               do ipy = 1,cgrid%npolygons
                  cgrid%met(ipy)%nir_beam = cgrid%metinput(ipy)%nbdsf(mhi) * t1            &
                                          + cgrid%metinput(ipy)%nbdsf(mlo) * t2
               end do

            case('nddsf')   !----- Near IR diffuse downward shortwave flux. --- [   W/m²] -!
               do ipy = 1,cgrid%npolygons
                  cgrid%met(ipy)%nir_diffuse = cgrid%metinput(ipy)%nddsf(mhi) * t1         &
                                             + cgrid%metinput(ipy)%nddsf(mlo) * t2
               end do

            case('vbdsf')   !----- Visible beam downward shortwave flux. ------ [   W/m²] -!
               do ipy = 1,cgrid%npolygons
                  cgrid%met(ipy)%par_beam = cgrid%metinput(ipy)%vbdsf(mhi) * t1            &
                                          + cgrid%metinput(ipy)%vbdsf(mlo) * t2
               end do

            case('vddsf')   !----- Visible diffuse downward shortwave flux. --- [   W/m²] -!
               do ipy = 1,cgrid%npolygons
                  cgrid%met(ipy)%par_diffuse = cgrid%metinput(ipy)%vddsf(mhi) * t1         &
                                             + cgrid%metinput(ipy)%vddsf(mlo) * t2
               enddo

            case('co2')     !----- CO2 mixing ratio. -------------------------- [    ppm] -!
               do ipy = 1,cgrid%npolygons
                  cgrid%met(ipy)%atm_co2 = cgrid%metinput(ipy)%co2(mhi) * t1               &
                                         + cgrid%metinput(ipy)%co2(mlo) * t2
               end do

            end select
         end if
      end do varloop
   end do formloop
   
   
   !---------------------------------------------------------------------------------------!
   !     Change from velocity squared to velocity, and compute qpcpg and dpcpg.            !
   !---------------------------------------------------------------------------------------!
   polyloop: do ipy = 1,cgrid%npolygons
         
      !----- CO2 --------------------------------------------------------------------------!
      if (.not.have_co2) cgrid%met(ipy)%atm_co2 = initial_co2

     !! ADJUST MET FOR SIMPLE CLIMATE SCENARIOS
     T0 = cgrid%met(ipy)%atm_tmp
     cgrid%met(ipy)%atm_tmp = atm_tmp_intercept + &
          cgrid%met(ipy)%atm_tmp * atm_tmp_slope
     cgrid%met(ipy)%pcpg = prec_intercept + &
          cgrid%met(ipy)%pcpg * prec_slope
     if(humid_scenario == 1) then 
        !! change SHV to keep RH const
        cgrid%met(ipy)%atm_shv = cgrid%met(ipy)%atm_shv * &
             exp(5415*(1.0/T0 - 1.0/cgrid%met(ipy)%atm_tmp))
     endif


      !------------------------------------------------------------------------------------!
      !     Checking whether the specific humidity provided will not cause supersatur-     !
      ! ation.  We cannot do it using specific humidity because for that we would need     !
      ! density, which needs the actual value to be computed consistently. So here we      !
      ! use mixing ratio, which depends on pressure, then convert back to specific         !
      ! humidity.                                                                          !
      !------------------------------------------------------------------------------------!
      rvaux = cgrid%met(ipy)%atm_shv / (1. - cgrid%met(ipy)%atm_shv)
      rvsat = rslif(cgrid%met(ipy)%prss,cgrid%met(ipy)%atm_tmp)
      rvaux = min(rvaux,rvsat)

      !------ Converting back to specific humidity. ---------------------------------------!
      cgrid%met(ipy)%atm_shv = rvaux / (1. + rvaux)

      !------------------------------------------------------------------------------------!
      !    We now find the Exner function, potential temperature, and the enthalpy.        !
      !------------------------------------------------------------------------------------!
      cgrid%met(ipy)%exner        = cp * (p00i * cgrid%met(ipy)%prss)**rocp
      cgrid%met(ipy)%atm_theta    = cp * cgrid%met(ipy)%atm_tmp / cgrid%met(ipy)%exner

      if (cgrid%met(ipy)%atm_tmp < 180.   .or. cgrid%met(ipy)%atm_tmp > 400.   .or.        &
          cgrid%met(ipy)%atm_shv < 1.e-8  .or. cgrid%met(ipy)%atm_shv > 0.04   .or.        &
          cgrid%met(ipy)%prss    < 40000. .or. cgrid%met(ipy)%prss    > 110000.) then
          write (unit=*,fmt='(a)') '======== Weird polygon properties... ========'
          write (unit=*,fmt='(a,f7.2)') 'POLY_PRSS [ hPa] = ',cgrid%met(ipy)%prss    * 0.01
          write (unit=*,fmt='(a,f7.2)') 'POLY_TEMP [degC] = ',cgrid%met(ipy)%atm_tmp - t00
          write (unit=*,fmt='(a,f7.2)') 'POLY_SHV  [g/kg] = ',cgrid%met(ipy)%atm_shv * 1.e3
          call fatal_error('Non-sense polygon met values!!!'                               &
                          ,'update_met_drivers','ed_met_driver.f90')
      end if

      cgrid%met(ipy)%atm_enthalpy = ptqz2enthalpy(cgrid%met(ipy)%prss                      &
                                                 ,cgrid%met(ipy)%atm_tmp                   &
                                                 ,cgrid%met(ipy)%atm_shv                   &
                                                 ,cgrid%met(ipy)%geoht)

      !------ Apply met to sites, and adjust met variables for topography. ----------------!
      call calc_met_lapse(cgrid,ipy)

      

      cpoly => cgrid%polygon(ipy)
      siteloop: do isi = 1,cpoly%nsites
         
         !----- Vels.  At this point vels is 2*Kinetic Energy, take the square root. ------!
         cpoly%met(isi)%vels = sqrt(max(0.0,cpoly%met(isi)%vels))
         cpoly%met(isi)%vels_stab = max(ubmin,cpoly%met(isi)%vels)
         cpoly%met(isi)%vels_unstab = max(ubmin,cpoly%met(isi)%vels)
         
         !----- CO2.  In case we used the namelist, use that value. -----------------------!
         if (.not.have_co2) cpoly%met(isi)%atm_co2 = initial_co2
         
         !---------------------------------------------------------------------------------!
         !     We now find some derived properties.  In case several sites exist, the      !
         ! lapse rate was applied to pressure, temperature, and mixing ratio.  Then we     !
         ! calculate the Exner function, potential temperature and enthalpy so it will     !
         ! respect the ideal gas law and first law of thermodynamics.                      !
         !---------------------------------------------------------------------------------!
         cpoly%met(isi)%exner        = cp * (p00i * cpoly%met(isi)%prss) **rocp
         cpoly%met(isi)%atm_theta    = cp * cpoly%met(isi)%atm_tmp / cpoly%met(isi)%exner


         if (cpoly%met(isi)%atm_tmp < 180.   .or. cpoly%met(isi)%atm_tmp > 400.   .or.     &
             cpoly%met(isi)%atm_shv < 1.e-8  .or. cpoly%met(isi)%atm_shv > 0.04   .or.     &
             cpoly%met(isi)%prss    < 40000. .or. cpoly%met(isi)%prss    > 110000.) then
             write (unit=*,fmt='(a)') '======== Weird site properties... ========'
             write (unit=*,fmt='(a,f7.2)')                                                 &
                                  'SITE_PRSS [ hPa] = ',cpoly%met(isi)%prss    * 0.01
             write (unit=*,fmt='(a,f7.2)')                                                 &
                                  'SITE_TEMP [degC] = ',cpoly%met(isi)%atm_tmp - t00
             write (unit=*,fmt='(a,f7.2)')                                                 &
                                  'SITE_SHV  [g/kg] = ',cpoly%met(isi)%atm_shv * 1.e3
             call fatal_error('Non-sense site met values!!!'                               &
                             ,'update_met_drivers','ed_met_driver.f90')
         end if



         cpoly%met(isi)%atm_enthalpy = ptqz2enthalpy(cpoly%met(isi)%prss                   &
                                                   ,cpoly%met(isi)%atm_tmp                 &
                                                   ,cpoly%met(isi)%atm_shv                 &
                                                   ,cpoly%met(isi)%geoht)

         !----- Solar radiation -----------------------------------------------------------!
         cpoly%met(isi)%rshort_diffuse = cpoly%met(isi)%par_diffuse                        &
                                       + cpoly%met(isi)%nir_diffuse
         cpoly%met(isi)%rshort         = cpoly%met(isi)%rshort_diffuse                     &
                                       + cpoly%met(isi)%par_beam + cpoly%met(isi)%nir_beam
         
         !---------------------------------------------------------------------------------!
         !     Checking whether the specific humidity provided will not cause supersatur-  !
         ! ation.  This must be done at the site level too, in case various sites exist in !
         ! a polygon.  We cannot do it using specific humidity because for that we would   !
         ! need density, which needs the actual value to be computed consistently. So here !
         ! we use mixing ratio, which depends on pressure, then convert back to specific   !
         ! humidity.                                                                       !
         !---------------------------------------------------------------------------------!
         rvaux = cpoly%met(isi)%atm_shv / (1. - cpoly%met(isi)%atm_shv)
         rvsat = rslif(cpoly%met(isi)%prss,cpoly%met(isi)%atm_tmp)
         rvaux = min(rvaux,rvsat)

         !------ Converting back to specific humidity. ------------------------------------!
         cpoly%met(isi)%atm_shv = rvaux / (1. + rvaux)


         !---------------------------------------------------------------------------------!
         !     Precipitation internal energy and "depth".  The decision on whether it's    !
         ! rain or snow is rather simple: above the triple point is rain, below it's snow. !
         !---------------------------------------------------------------------------------!
         if (cpoly%met(isi)%atm_tmp > t3ple) then
            cpoly%met(isi)%qpcpg = cliq * (cpoly%met(isi)%atm_tmp - tsupercool)            &
                                        * max(0.0, cpoly%met(isi)%pcpg)
            cpoly%met(isi)%dpcpg = max(0.0, cpoly%met(isi)%pcpg * wdnsi)
         else
            cpoly%met(isi)%qpcpg = cice * cpoly%met(isi)%atm_tmp * cpoly%met(isi)%pcpg
            
            !------------------------------------------------------------------------------!
            !  fcn derived from CLM3.0 documentation which is based on                     !
            !  Anderson 1975 NWS Technical Doc # 19                                        !
            !  which I have yet to find   <mcd>                                            !
            !------------------------------------------------------------------------------!
            cpoly%met(isi)%dpcpg = cpoly%met(isi)%pcpg &
               / (50.0+1.5*max(cpoly%met(isi)%atm_tmp-258.15,0.))**1.5
            !-------------------------------------------------------------------------------!
         end if
      end do siteloop
           
   end do polyloop

   return
end subroutine update_met_drivers
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine read_ol_file(iformat, iv, year_use, mname, year, offset, cgrid)

  use ed_state_vars,only:edgrid_g,edtype,polygontype
  use met_driver_coms, only: met_frq, met_nlon, met_nlat, met_vars,   &
       met_xmin, met_dx, met_ymin, met_dy,met_interp,no_ll

  use hdf5_utils, only: shdf5_irec_f,shdf5_info_f
  use mem_sites, only: grid_type
  use ed_state_vars,only:edtype
  use consts_coms, only: day_sec
  
  implicit none

  type(edtype),target :: cgrid
  integer, intent(in) :: iformat
  integer, intent(in) :: iv
  integer, intent(in) :: year_use
  character(len=3), intent(in) :: mname
  integer, intent(in) :: year
  integer, intent(in) :: offset
  integer :: ipy
  integer :: points_per_day
  integer :: nday,nday_dset
  integer :: np,ip,np_dset
  real, allocatable, dimension(:,:,:) :: metvar
  real, allocatable, dimension(:,:,:) :: metvar_dset
  real, allocatable, dimension(:) :: metscalar
  integer :: ndims
  integer, dimension(3) :: idims
  integer :: ilon
  integer :: ilat
  logical, external :: isleap

  ! This subroutine reads in the meteorogolical data
  ! subroutine shdf5_irec_f uses HDF5 routines to extract the data from the archive
  ! The rest of this code manages the allocation and storage of that data.
  ! The extracted data may be from a year, other than the current simulation year.
  ! This happens when we use a range of years to loop data over, typically during
  ! very long simulations. In such cases, the number of days in February will change
  ! during leap years, and the number of days in February during the model year may not match
  ! that of the simulation year

  ! Geoposition variables do not apply
  if (met_vars(iformat,iv).eq.'lat' .or. met_vars(iformat,iv).eq.'lon') return
  
  !  Determine how many times to read

  !  In cases 2 and 4 we are dealing with constant time, only 1 point
  
  if (met_interp(iformat,iv).eq.2 .or. met_interp(iformat,iv).eq.4) then
     np = 1
  else

     ! Allocate for the number of points in the current year.  If the model year
     ! has more data in February, we will not use that last (29th) day.  If the 
     ! February in model time has more days, then we have to recycle the last day
     ! of the data archive if it is not a leap year also.

     points_per_day = nint(day_sec/met_frq(iformat,iv))

     ! Determine the number of points necessary for the model year
     select case (trim(mname))
     case ('APR','JUN','SEP','NOV')
        nday = 30
        nday_dset = 30
     case ('JAN','MAR','MAY','JUL','AUG','OCT','DEC')
        nday = 31
        nday_dset = 31
     case ('FEB')
        if (isleap(year)) then
           nday = 29
        else 
           nday = 28
        end if
        if (isleap(year_use)) then
           nday_dset=29
        else
           nday_dset=28
        end if
     end select 

     np = nday * points_per_day
     np_dset = nday_dset * points_per_day

  endif

  ! Allocate the buffer space
  
  select case (met_interp(iformat,iv))
  !  Case 3, reading in 1 value representing the whole grid,
  !  but over all available times
  case (3) 
     ndims = 1
     idims(1) = np_dset
     
     allocate(metvar(np,met_nlon(iformat),met_nlat(iformat)))
     allocate(metscalar(np_dset))
     call shdf5_irec_f(ndims, idims, trim(met_vars(iformat,iv)),  &
          rvara = metscalar)
     
     if (np < np_dset) then     ! No interpolation is needed, but we
                                ! need two buffers
        do ip=1,np
           metvar(ip,:,:) = metscalar(ip)
        enddo
        
     elseif(np .eq. np_dset) then      ! Datasets are the same, no double allocation
        
        do ip=1,np
           metvar(ip,:,:) = metscalar(ip)
        enddo
        
     else              ! The dataset buffer is smaller than the models
                       ! buffer, so we recycle the last day
        do ip=1,np_dset
           metvar(ip,:,:) = metscalar(ip)
        enddo

        do ip=np_dset+1,np
           ! Then reuse the 28th for the 29th day
           metvar(ip,:,:) = metscalar(ip-points_per_day+1)
        enddo
     endif

     deallocate(metscalar)

  !  Case 4, dont read, but use the value stored in met_freq
  !  as the scalar value
  case (4)
     
     allocate(metvar(1,met_nlon(iformat),met_nlat(iformat)))
     metvar(1,:,:) = met_frq(iformat,iv)
     
  !  Case 2, 2D grid is read in, ie. np = 1
  case (2)

     allocate(metvar(1,met_nlon(iformat),met_nlat(iformat)))

     call shdf5_info_f(trim(met_vars(iformat,iv)),ndims,idims)
     
     ! Check to see if the dimensions are consistent
     if (ndims.ne.2                          &
          .or. idims(1).ne.met_nlon(iformat) &
          .or. idims(2).ne.met_nlat(iformat)) then
        print*,"DATASET: ",trim(met_vars(iformat,iv))
        print*,"DOES NOT HAVE DIMENSIONS THAT MATCH THE"
        print*,"SPECIFIED INPUT, OR LAT/LON GRID"
        print*,ndims,np,met_nlon(iformat),met_nlat(iformat),idims
        call fatal_error('Mismatch between dataset and specified input' &
                        ,'read_ol_file','ed_met_driver.f90')
     endif

     call shdf5_irec_f(ndims, idims, trim(met_vars(iformat,iv)),  &
          rvara = metvar)

  ! Case default (aka 1 or 0), read time-and-space variable array.
  case default

     call shdf5_info_f(trim(met_vars(iformat,iv)),ndims,idims)
     
     ! Check to see if the dimensions are consistent
     if (ndims.ne.3                          &
          .or. idims(1).ne.np_dset           &
          .or. idims(2).ne.met_nlon(iformat) &
          .or. idims(3).ne.met_nlat(iformat)) then
        print*,"DATASET: ",trim(met_vars(iformat,iv))
        print*,"DOES NOT HAVE DIMENSIONS THAT MATCH THE"
        print*,"SPECIFIED INPUT, OR LAT/LON GRID"
        print*,ndims
        print*,np_dset,met_nlon(iformat),met_nlat(iformat)
        print*,idims
        call fatal_error('Mismatch between dataset and specified input' &
                        ,'read_ol_file','ed_met_driver.f90')
     endif

     if (np < np_dset) then     ! No interpolation is needed, but we
                                ! need two buffers
        
        allocate(metvar(np,met_nlon(iformat),met_nlat(iformat)))
        allocate(metvar_dset(np_dset,met_nlon(iformat),met_nlat(iformat)))
        
        call shdf5_irec_f(ndims, idims, trim(met_vars(iformat,iv)),  &
             rvara = metvar_dset)
        
        metvar(1:np,:,:) = metvar_dset(1:np,:,:)
        
        deallocate(metvar_dset)
        
        
     elseif(np .eq. np_dset) then      ! Datasets are the same, no double allocation
        
        allocate(metvar(np,met_nlon(iformat),met_nlat(iformat)))
        
        call shdf5_irec_f(ndims, idims, trim(met_vars(iformat,iv)),  &
             rvara = metvar)
        
     else                       ! The dataset buffer is smaller than the models
                                ! buffer, so we recycle the last day
        
        allocate(metvar(np,met_nlon(iformat),met_nlat(iformat)))
        allocate(metvar_dset(np_dset,met_nlon(iformat),met_nlat(iformat)))
        
        call shdf5_irec_f(ndims, idims, trim(met_vars(iformat,iv)),  &
             rvara = metvar_dset)
        
        metvar(1:np_dset,:,:) = metvar_dset(1:np_dset,:,:)
        
        ! Then reuse the 28th for the 29th day
        metvar(np_dset+1:np,:,:) = metvar_dset(np_dset-points_per_day+1:np_dset,:,:)
        
        deallocate(metvar_dset)

     endif


  end select


  ! The dataset for that variable has been read into an array
  ! Now, the polygons are cycled through, and they are associated
  ! with an index in the grid.  If this is a lat-lon grid, the procedure
  ! is a simple interpolation between the domain endpoints, if this is 
  ! a polar stereo-graphic projection, we make use of the associated
  ! indices saved in memory, ie. (cgrid%ilon(ipy))

  do ipy = 1,cgrid%npolygons
     
     ! Get the indices.  Remember, latitude is flipped.
     if(no_ll) then
        ilon = min( max(1, 1 + nint( (cgrid%lon(ipy) - met_xmin(iformat)) /   &
             met_dx(iformat) ) ), met_nlon(iformat))
        ilat = met_nlat(iformat) -  &
             min( max(1, 1 + nint( (cgrid%lat(ipy) - met_ymin(iformat)) /   &
             met_dy(iformat) ) ), met_nlat(iformat)) + 1
     else
        ilon = cgrid%ilon(ipy)
        ilat = cgrid%ilat(ipy)
     endif
     
     ! Get the time series.
     select case (trim(met_vars(iformat,iv)))
     case('nbdsf')
        cgrid%metinput(ipy)%nbdsf((offset+1):(offset+np)) =   &
             metvar(1:np,ilon,ilat)
     case('nddsf')
        cgrid%metinput(ipy)%nddsf((offset+1):(offset+np)) =   &
             metvar(1:np,ilon,ilat)
     case('vbdsf')
        cgrid%metinput(ipy)%vbdsf((offset+1):(offset+np)) =   &
             metvar(1:np,ilon,ilat)
     case('vddsf')
        cgrid%metinput(ipy)%vddsf((offset+1):(offset+np)) =   &
             metvar(1:np,ilon,ilat)
     case('prate')
        cgrid%metinput(ipy)%prate((offset+1):(offset+np)) =   &
             metvar(1:np,ilon,ilat)
     case('dlwrf')
        cgrid%metinput(ipy)%dlwrf((offset+1):(offset+np)) =   &
             metvar(1:np,ilon,ilat)
     case('pres')
        cgrid%metinput(ipy)%pres((offset+1):(offset+np)) =    &
             metvar(1:np,ilon,ilat)
     case('hgt')
        cgrid%metinput(ipy)%hgt((offset+1):(offset+np))  =    &
             metvar(1:np,ilon,ilat)
     case('ugrd')
        cgrid%metinput(ipy)%ugrd((offset+1):(offset+np)) =    &
             metvar(1:np,ilon,ilat)
     case('vgrd')
        cgrid%metinput(ipy)%vgrd((offset+1):(offset+np)) =    &
             metvar(1:np,ilon,ilat)
     case('sh')
        cgrid%metinput(ipy)%sh((offset+1):(offset+np)) =      &
             metvar(1:np,ilon,ilat)
     case('tmp')
        cgrid%metinput(ipy)%tmp((offset+1):(offset+np)) =     &
             metvar(1:np,ilon,ilat)
     case('co2')
        cgrid%metinput(ipy)%co2((offset+1):(offset+np)) =     &
             metvar(1:np,ilon,ilat)
     end select
     
  enddo
  
  ! Deallocate the buffer
  deallocate(metvar)
  
  return
end subroutine read_ol_file
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine transfer_ol_month(vname, frq, cgrid)
  
  use ed_state_vars,only : edtype
  use consts_coms, only: day_sec

  implicit none

  type(edtype),target :: cgrid
  character*(*) :: vname
  real, intent(in) :: frq
  integer :: mem_size
  integer :: ipy

  mem_size = nint(day_sec / frq) * 31

  do ipy = 1,cgrid%npolygons

     ! Get the time series.
     select case(trim(vname))
     case('nbdsf')
        cgrid%metinput(ipy)%nbdsf(1:mem_size) =   &
             cgrid%metinput(ipy)%nbdsf((mem_size+1):(2*mem_size))
     case('nddsf')
        cgrid%metinput(ipy)%nddsf(1:mem_size) =   &
             cgrid%metinput(ipy)%nddsf((mem_size+1):(2*mem_size))
     case('vbdsf')
        cgrid%metinput(ipy)%vbdsf(1:mem_size) =   &
             cgrid%metinput(ipy)%vbdsf((mem_size+1):(2*mem_size))
     case('vddsf')
        cgrid%metinput(ipy)%vddsf(1:mem_size) =   &
             cgrid%metinput(ipy)%vddsf((mem_size+1):(2*mem_size))
     case('prate')
        cgrid%metinput(ipy)%prate(1:mem_size) =   &
             cgrid%metinput(ipy)%prate((mem_size+1):(2*mem_size))
     case('dlwrf')
        cgrid%metinput(ipy)%dlwrf(1:mem_size) =   &
             cgrid%metinput(ipy)%dlwrf((mem_size+1):(2*mem_size))
     case('pres')
        cgrid%metinput(ipy)%pres(1:mem_size) =   &
             cgrid%metinput(ipy)%pres((mem_size+1):(2*mem_size))
     case('hgt')
        cgrid%metinput(ipy)%hgt(1:mem_size) =   &
             cgrid%metinput(ipy)%hgt((mem_size+1):(2*mem_size))
     case('ugrd')
        cgrid%metinput(ipy)%ugrd(1:mem_size) =   &
             cgrid%metinput(ipy)%ugrd((mem_size+1):(2*mem_size))
     case('vgrd')
        cgrid%metinput(ipy)%vgrd(1:mem_size) =   &
             cgrid%metinput(ipy)%vgrd((mem_size+1):(2*mem_size))
     case('sh')
        cgrid%metinput(ipy)%sh(1:mem_size) =   &
             cgrid%metinput(ipy)%sh((mem_size+1):(2*mem_size))
     case('tmp')
        cgrid%metinput(ipy)%tmp(1:mem_size) =   &
             cgrid%metinput(ipy)%tmp((mem_size+1):(2*mem_size))
     case('co2')
        cgrid%metinput(ipy)%co2(1:mem_size) =   &
             cgrid%metinput(ipy)%co2((mem_size+1):(2*mem_size))
     end select
     

  enddo
  
  return
end subroutine transfer_ol_month
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine match_poly_grid(cgrid,nlon,nlat,lon,lat)

   use ed_state_vars , only : edtype

   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(edtype)              , target        :: cgrid
   real, dimension(nlon,nlat), intent(inout) :: lon,lat
   integer                   , intent(in)    :: nlon,nlat
   !----- Local variables -----------------------------------------------------------------!
   integer                                   :: ilat,ilon
   integer                                   :: ipy
   real                                      :: min_dist
   real                                      :: this_dist
   !----- External function ---------------------------------------------------------------!
   real                       , external     :: dist_gc
   !---------------------------------------------------------------------------------------!


   do ipy = 1,cgrid%npolygons
      min_dist = huge (1.)

      do ilon = 1,nlon
         do ilat = 1,nlat

            if (lon(ilon,ilat) > 180.0)                                                    &
               lon(ilon,ilat) = lon(ilon,ilat) - 360.0

            this_dist = dist_gc(cgrid%lon(ipy), lon(ilon,ilat)                      &
                               ,cgrid%lat(ipy), lat(ilon,ilat) )
            
            if(this_dist < min_dist) then               
               cgrid%ilon(ipy) = ilon
               cgrid%ilat(ipy) = ilat
               min_dist        = this_dist
            endif
         end do
      end do
   end do

   return
end subroutine match_poly_grid
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine getll(cgrid,iformat)

   use met_driver_coms, only: met_nv, &
        met_vars, &
        met_nlon, &
        met_nlat, &
        lat2d,    &
        lon2d,    &
        met_xmin, &
        met_ymin, &
        met_dx,   &
        met_dy,   &
        no_ll
   
   use hdf5_utils,only : shdf5_info_f,shdf5_irec_f
   use ed_state_vars,only:edtype

   implicit none
   
   integer :: iformat
   integer :: ndims
   integer :: d
   integer, dimension(3) :: idims
   type(edtype),target :: cgrid

   ! First check to see if there is lat/lon data in this dataset
   ! if the data exists, load it
   
   if(.not.no_ll) then
         
      !  Get the dimensioning information on latitude
      call shdf5_info_f('lat',ndims,idims)
      
      if(ndims /= 2) then
         write(unit=*,fmt='(a)') 'Number of dimensions of latitude is wrong...'
         write(unit=*,fmt='(a,1x,i5)') 'NDIMS=',ndims
         do d=1,ndims
            write(unit=*,fmt='(a,1x,i5)') '---> ',d,': DIM=',idims(d)
         end do
         call fatal_error ('Not set up to have time varying latitude...' &
                          ,'getll','ed_met_driver.f90')
      endif
      
      !  Transfer the dimensions into the met_nlon array
      met_nlon(iformat) = idims(1)
      met_nlat(iformat) = idims(2)
      
      !  Allocate the latitude array
      allocate(lat2d(idims(1),idims(2)))
      
      !  Read in the latitude array
      call shdf5_irec_f(ndims, idims, 'lat',  &
           rvara = lat2d )
            
      
      !  Get the dimensioning information on longitude
      call shdf5_info_f('lon',ndims,idims)
      
      if(ndims /= 2) then
         write(unit=*,fmt='(a)') 'Number of dimensions of longitude is wrong...'
         write(unit=*,fmt='(a,1x,i5)') 'NDIMS=',ndims
         do d=1,ndims
            write(unit=*,fmt='(a,1x,i5)') '---> ',d,': DIM=',idims(d)
         end do
         call fatal_error ('Not set up to have time varying longitude...' &
                          ,'getll','ed_met_driver.f90')
      endif
      
      !  Allocate the latitude array
      allocate(lon2d(idims(1),idims(2)))
      ndims = 2
      
      !  Read in the latitude array
      call shdf5_irec_f(ndims, idims, 'lon',  &
           rvara = lon2d )
      
      !  Determine the indices of the grid that each polygon sees
      !  returns poly%ilon and poly%ilat
      
      call match_poly_grid(cgrid,met_nlon(iformat),met_nlat(iformat),lon2d,lat2d)
      
      ! Deallocate the lat-lon arrays
      deallocate(lat2d,lon2d)

   endif
   
   return
   
end subroutine getll
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine calc_met_lapse(cgrid,ipy)
  
  use ed_state_vars,only    : edtype,polygontype
  use canopy_radiation_coms,only : rlong_min
  use met_driver_coms, only : lapse_scheme


  implicit none
  integer, intent(in) :: ipy
  integer :: isi
  type(edtype),target       :: cgrid
  type(polygontype),pointer :: cpoly
  
  real            :: ebar   !! mean elevation
  real            :: delE   !! deviation from mean elevation
  real            :: aterr  !! terrestrial area
  real, parameter :: offset=tiny(1.)/epsilon(1.) !! Tiny offset to avoid FPE
  real :: hillshade

  cpoly => cgrid%polygon(ipy)

  if (lapse_scheme == 0) then  !!!!! bypass
     do isi = 1,cpoly%nsites
        hillshade = cpoly%slope(isi)/180.
        cpoly%met(isi)%geoht       = cgrid%met(ipy)%geoht
        cpoly%met(isi)%atm_tmp     = cgrid%met(ipy)%atm_tmp
        cpoly%met(isi)%atm_shv     = cgrid%met(ipy)%atm_shv
        cpoly%met(isi)%prss        = cgrid%met(ipy)%prss
        cpoly%met(isi)%pcpg        = cgrid%met(ipy)%pcpg
        cpoly%met(isi)%par_diffuse = cgrid%met(ipy)%par_diffuse*(1.0-hillshade)
        cpoly%met(isi)%atm_co2     = cgrid%met(ipy)%atm_co2
        cpoly%met(isi)%rlong       = cgrid%met(ipy)%rlong
        cpoly%met(isi)%par_beam    = cgrid%met(ipy)%par_beam
        cpoly%met(isi)%nir_diffuse = cgrid%met(ipy)%nir_diffuse*(1.0-hillshade)
        cpoly%met(isi)%nir_beam    = cgrid%met(ipy)%nir_beam
        cpoly%met(isi)%vels        = cgrid%met(ipy)%vels

!        if ( cpoly%met(isi)%rlong < rlong_min) then
!           print*,"Problems with RLONG A"  
!           print*,cpoly%met(isi)%rlong,cgrid%met(ipy)%rlong,rlong_min
!           call fatal_error('Problems with RLONG A','calc_met_lapse','ed_met_driver.f90')
!        end if
!        if ( cpoly%met(isi)%atm_tmp < 150.0) then
!           print*,cpoly%met(isi)%atm_tmp,cgrid%met(ipy)%atm_tmp
!           call fatal_error('Problems with ATM TEMP A','calc_met_lapse','ed_met_driver.f90')
!        endif
        if ( cpoly%met(isi)%atm_shv < 0.01e-2) then
            cpoly%met(isi)%atm_shv = 0.01e-2
!           print*,cpoly%met(isi)%atm_shv,cgrid%met(ipy)%atm_shv
!           call fatal_error('Problems with ATM MOISTURE A','calc_met_lapse','ed_met_driver.f90')
        endif

     enddo
   
  else

     !pass over sites once to calc preliminary stats
     ebar = 0.0
     aterr = 0.0
     do isi=1,cpoly%nsites
        ebar = ebar + cpoly%area(isi)*cpoly%elevation(isi)
        aterr = aterr + cpoly%area(isi)
     enddo
     ebar = ebar/aterr

     !!second pass, calc lapse rate adjustment
     do isi = 1,cpoly%nsites
        
        delE = cpoly%elevation(isi) - ebar
        
        !! perform linear adjustments
        cpoly%met(isi)%geoht   = cgrid%met(ipy)%geoht   + cgrid%lapse(ipy)%geoht*delE
        cpoly%met(isi)%atm_tmp = cgrid%met(ipy)%atm_tmp + cgrid%lapse(ipy)%atm_tmp*delE
        cpoly%met(isi)%atm_shv = cgrid%met(ipy)%atm_shv + cgrid%lapse(ipy)%atm_shv*delE
        cpoly%met(isi)%prss    = cgrid%met(ipy)%prss    + cgrid%lapse(ipy)%prss*delE
        cpoly%met(isi)%atm_co2 = cgrid%met(ipy)%atm_co2 + cgrid%lapse(ipy)%atm_co2*delE
        cpoly%met(isi)%rlong   = cgrid%met(ipy)%rlong   + cgrid%lapse(ipy)%rlong*delE
        cpoly%met(isi)%par_diffuse = (cgrid%met(ipy)%par_diffuse + &
             cgrid%lapse(ipy)%par_diffuse*delE)*(1.0-hillshade)
        cpoly%met(isi)%par_beam    = cgrid%met(ipy)%par_beam    + cgrid%lapse(ipy)%par_beam*delE
        cpoly%met(isi)%nir_diffuse = (cgrid%met(ipy)%nir_diffuse + &
             cgrid%lapse(ipy)%nir_diffuse*delE)*(1.0-hillshade)
        cpoly%met(isi)%nir_beam    = cgrid%met(ipy)%nir_beam    + cgrid%lapse(ipy)%nir_beam*delE
        cpoly%met(isi)%vels    = cgrid%met(ipy)%vels    + cgrid%lapse(ipy)%vels*delE
        !! PROPORTIONAL ADJUSTMENTS
        cpoly%met(isi)%pcpg    = cgrid%met(ipy)%pcpg*cpoly%pptweight(isi)

        !! note: at this point VELS is vel^2.  Thus this lapse preserves mean wind ENERGY
        !! not wind SPEED
!        if ( cpoly%met(isi)%rlong < 200.0) then
!           print*,"Problems with RLONG A"  
!           print*,cpoly%met(isi)%rlong,cgrid%met(ipy)%rlong
!           call fatal_error('Problems with RLONG A','calc_met_lapse','ed_met_driver.f90')
!        end if
!        if ( cpoly%met(isi)%atm_tmp < 150.0) then
!           print*,cpoly%met(isi)%atm_tmp,cgrid%met(ipy)%atm_tmp
!           call fatal_error('Problems with ATM TEMP A','calc_met_lapse','ed_met_driver.f90')
!        endif
        if ( cpoly%met(isi)%atm_shv < 0.01e-2) then
            cpoly%met(isi)%atm_shv = 0.01e-2
!           print*,cpoly%met(isi)%atm_shv,cgrid%met(ipy)%atm_shv
!           call fatal_error('Problems with ATM MOISTURE A','calc_met_lapse','ed_met_driver.f90')
        endif

        call MetDiagnostics(cpoly,ipy,isi)
        
     enddo
  endif

  return
end subroutine calc_met_lapse
!==========================================================================================!
!==========================================================================================!


subroutine MetDiagnostics(cpoly,ipy,isi)
  use ed_state_vars,only    : edtype,polygontype
  use canopy_radiation_coms,only : rlong_min
  
  implicit none
  integer, intent(in) :: ipy
  integer, intent(in) :: isi
  type(polygontype),target :: cpoly

  logical, parameter :: fixRad = .true.


  if(cpoly%met(isi)%geoht .le. 0.)  then
     print*,cpoly%met(isi)%geoht
     call fatal_error('Problems with GEOHT','MetDiagnostics','ed_met_driver.f90')
  endif
  if(cpoly%met(isi)%atm_tmp .le. 200. .or. cpoly%met(isi)%atm_tmp .ge. 350.)  then
     print*,cpoly%met(isi)%atm_tmp
     call fatal_error('Problems with atm_tmp','MetDiagnostics','ed_met_driver.f90')
  endif
  if(cpoly%met(isi)%atm_shv .le. 0. .or. cpoly%met(isi)%atm_shv .ge. 1.)  then
     print*,cpoly%met(isi)%atm_shv
     call fatal_error('Problems with atm_shv','MetDiagnostics','ed_met_driver.f90')
  endif
  if(cpoly%met(isi)%prss .le. 0.)  then
     print*,cpoly%met(isi)%prss
     call fatal_error('Problems with prss','MetDiagnostics','ed_met_driver.f90')
  endif
  if(cpoly%met(isi)%pcpg .lt. 0. .or. cpoly%met(isi)%pcpg .gt. 0.1)  then
     print*,cpoly%met(isi)%pcpg
     call fatal_error('Problems with precipitation','MetDiagnostics','ed_met_driver.f90')
  endif
  if(cpoly%met(isi)%atm_co2 .le. 100. .or. cpoly%met(isi)%atm_co2 .gt. 2000.)  then
     print*,cpoly%met(isi)%atm_co2
     call fatal_error('Problems with ATM CO2','MetDiagnostics','ed_met_driver.f90')
  endif
  if(cpoly%met(isi)%rlong .le. rlong_min)  then
     print*,cpoly%met(isi)%rlong
     call fatal_error('Problems with LONGWAVE RADIATION','MetDiagnostics','ed_met_driver.f90')
  endif
  if(cpoly%met(isi)%par_diffuse .lt. 0. .or. cpoly%met(isi)%par_diffuse .gt. 1400.)  then
     if(cpoly%met(isi)%par_diffuse .lt. 0. .and. fixRad) then
        cpoly%met(isi)%par_diffuse = 0.0
     else
        print*,cpoly%met(isi)%par_diffuse
        call fatal_error('Problems with DIFFUSE PAR','MetDiagnostics','ed_met_driver.f90')
     endif
  endif
  if(cpoly%met(isi)%par_beam .lt. 0. .or. cpoly%met(isi)%par_beam .gt. 1400.)  then
     if(cpoly%met(isi)%par_beam .lt. 0. .and. fixRad) then
        cpoly%met(isi)%par_beam = 0.0
     else
        print*,cpoly%met(isi)%par_beam
        call fatal_error('Problems with DIRECT BEAM PAR','MetDiagnostics','ed_met_driver.f90')
     end if
  endif
  if(cpoly%met(isi)%nir_diffuse .lt. 0. .or. cpoly%met(isi)%nir_diffuse .gt. 1400.)  then
     if(cpoly%met(isi)%nir_diffuse .lt. 0. .and. fixRad) then
        cpoly%met(isi)%nir_diffuse = 0.0
     else
        print*,cpoly%met(isi)%nir_diffuse
        call fatal_error('Problems with DIFUSE NIR','MetDiagnostics','ed_met_driver.f90')
     endif
  end if

  if(cpoly%met(isi)%nir_beam .lt. 0. .or. cpoly%met(isi)%nir_beam .gt. 1400.)  then
     if(cpoly%met(isi)%nir_beam .lt. 0. .and. fixRad) then
        cpoly%met(isi)%nir_beam = 0.0
     else
        print*,cpoly%met(isi)%nir_beam
        call fatal_error('Problems with DIRECT BEAM NIR','MetDiagnostics','ed_met_driver.f90')
     end if
  endif
  if(cpoly%met(isi)%vels .lt. 0.)  then
     print*,cpoly%met(isi)%vels
     call fatal_error('Problems with WIND VELOCITY','MetDiagnostics','ed_met_driver.f90')
  endif

  return
end subroutine MetDiagnostics




!==========================================================================================!
!==========================================================================================!
subroutine setLapseParms(cgrid)
  
  use ed_state_vars,only:edtype,polygontype
  use met_driver_coms, only: lapse

  implicit none

  integer :: ipy,isi
  type(edtype), target :: cgrid
  type(polygontype),pointer :: cpoly
  real :: ebar, aterr

  !! right now, simply transfer lapse rates from ed_data
  !! in future, could set parms based on spatial maps of parms
  
  do ipy = 1,cgrid%npolygons
     
     cgrid%lapse(ipy)%geoht   = lapse%geoht
     cgrid%lapse(ipy)%vels    = lapse%vels
     cgrid%lapse(ipy)%atm_tmp = lapse%atm_tmp
     cgrid%lapse(ipy)%atm_shv = lapse%atm_shv
     cgrid%lapse(ipy)%prss    = lapse%prss
     cgrid%lapse(ipy)%pcpg    = lapse%pcpg
     cgrid%lapse(ipy)%atm_co2 = lapse%atm_co2
     cgrid%lapse(ipy)%rlong   = lapse%rlong
     cgrid%lapse(ipy)%nir_beam    = lapse%nir_beam
     cgrid%lapse(ipy)%nir_diffuse = lapse%nir_diffuse
     cgrid%lapse(ipy)%par_beam    = lapse%par_beam
     cgrid%lapse(ipy)%par_diffuse = lapse%par_diffuse
     cgrid%lapse(ipy)%pptnorm = lapse%pptnorm

     !!! PRECIP WEIGHTS !!!
     cpoly => cgrid%polygon(ipy)
     ebar = 0.0
     aterr = 0.0
     do isi = 1,cpoly%nsites
        ebar = ebar + cpoly%area(isi)*cpoly%elevation(isi)
        aterr = aterr + cpoly%area(isi)
     enddo
     ebar = ebar/aterr
     do isi = 1,cpoly%nsites
        cpoly%pptweight(isi) = (cgrid%lapse(ipy)%pptnorm + &
             cgrid%lapse(ipy)%pcpg*(cpoly%elevation(isi) - ebar))&
             /cgrid%lapse(ipy)%pptnorm
     enddo

  enddo

  return
end subroutine setLapseParms
!==========================================================================================!
!==========================================================================================!
