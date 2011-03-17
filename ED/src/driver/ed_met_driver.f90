!==========================================================================================!
!==========================================================================================!
!     This subroutine is the main procedure that will take care of reading the meteoro-    !
! logic data to be used during the run.                                                    !
!------------------------------------------------------------------------------------------!
subroutine read_met_driver_head()
   use ed_max_dims    , only : max_met_vars     ! ! intent(in)
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
   allocate(no_ll     (nformats)              )

   !----- Just to initialize, if lon/lat are both found, it will become .false. -----------!
   no_ll(:) = .true.    
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
         no_ll (iformat) = .false.
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

   use ed_max_dims     , only : max_met_vars
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
            if(no_ll(iformat) .and.                                                        &
               ( cgrid%lon(ipy) < westedge .or. cgrid%lat(ipy) < southedge .or.            &
                 cgrid%lon(ipy) > eastedge .or. cgrid%lat(ipy) > northedge)     ) then
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
               case(1,5)    !----- Interpolated.   Read two months in. ---------------------!
                  mem_size = 2 * nint(day_sec / met_frq(iformat,iv)) * 31
               case (2,4)  !----- Constant in time. ---------------------------------------!
                  mem_size = 1
               case default
                  call fatal_error('Invalid met_interp! It should be between 0 and 5.'     &
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
                                  ,'init_met_drivers','ed_met_driver.f90')
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

   use ed_max_dims       , only : str_len        ! ! intent(in)
   use ed_state_vars     , only : edgrid_g       & ! structure
                                , edtype         & ! structure
                                , polygontype    ! ! structure
   use hdf5_utils        , only : shdf5_open_f   & ! subroutine
                                , shdf5_close_f  ! ! subroutine
   use met_driver_coms   , only : nformats       & ! intent(in)
                                , ishuffle       & ! intent(in)
                                , met_names      & ! intent(in)
                                , met_nv         & ! intent(in)
                                , met_interp     & ! intent(in)
                                , met_frq        & ! intent(in)
                                , metcyc1        & ! intent(in)
                                , metcycf        & ! intent(in)
                                , metyears       ! ! intent(inout)
   use mem_polygons      , only : grid_type      ! ! intent(in)
   use ed_misc_coms      , only : current_time   & ! intent(in)
                                , iyeara         & ! intent(in)
                                , iyearz         & ! intent(in)
                                , imonthz        ! ! intent(in)
   use grid_coms         , only : ngrids         ! ! intent(in)
   use consts_coms       , only : day_sec        ! ! intent(in)
   implicit none
   !----- Local variables -----------------------------------------------------------------!
   type(edtype)          , pointer :: cgrid
   character(len=str_len)          :: infile
   integer                         :: igr
   integer                         :: year_use
   integer                         :: iformat
   integer                         :: iv
   integer                         :: offset
   integer                         :: m2
   integer                         :: y2
   integer                         :: year_use_2
   integer                         :: nyears
   integer                         :: ncyc
   integer                         :: nfullcyc
   integer                         :: i1stfull
   integer                         :: iyear
   integer, dimension(8)           :: seedtime
   real                            :: runif
   logical                         :: exans
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

      !------------------------------------------------------------------------------------!
      !     Define the number of years we have meteorological information available.  In   !
      ! case the simulation is supposed to end in December, we must add an extra year, be- !
      ! cause the meteorological forcing may need to access the following month.           !
      !------------------------------------------------------------------------------------!
      ncyc   = metcycf - metcyc1 + 1
      if (imonthz == 12) then
         nyears = iyearz  - iyeara  + 2
      else
         nyears = iyearz  - iyeara  + 1
      end if

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
            print*,"METYEARS",iyear,year_use,metcyc1,metcycf
            metyears(iyear) = year_use
         end do

      !------------------------------------------------------------------------------------!
      !    For years outside the met driver range, we pick up a random year.  Because we   !
      ! use the actual sequence of random years given by the code random generator, this   !
      ! sequence will be always the same for a given met driver and first year.            !
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
            print*,iyear,year_use
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
            call read_ol_file(infile,iformat, iv, year_use, mname(current_time%month)      &
                             ,current_time%year, offset, cgrid)
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
            call fatal_error('Cannot open met driver input file '//trim(infile)//'!'       &
                            ,'read_met_drivers_init','ed_met_driver.f90')
         end if
         
         !----- Loop over variables. ------------------------------------------------------!
         varloop: do iv = 1, met_nv(iformat)

            select case (met_interp(iformat,iv))
            case (1,5)
               offset = nint(day_sec / met_frq(iformat,iv)) * 31

               !----- Read the file. ------------------------------------------------------!
               call read_ol_file(infile,iformat,iv,year_use_2,mname(m2),y2,offset,cgrid)
            end select

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
   use ed_max_dims    , only : str_len       ! ! intent(in)
   use ed_state_vars  , only : edgrid_g      & ! structure
                             , edtype        & ! structure
                             , polygontype   ! ! structure
   use mem_polygons   , only : grid_type     ! ! structure
   use met_driver_coms, only : nformats      & ! intent(in)
                             , met_names     & ! intent(in)
                             , met_nv        & ! intent(in)
                             , met_interp    & ! intent(in)
                             , met_frq       & ! intent(in)
                             , met_vars      & ! intent(in)
                             , metcyc1       & ! intent(in)
                             , metcycf       ! ! intent(in)
   use ed_misc_coms   , only : current_time  ! ! intent(in)
   use hdf5_utils     , only : shdf5_open_f  & ! subroutine
                             , shdf5_close_f ! ! subroutine
   use ed_state_vars  , only : edtype        & ! structure
                             , edgrid_g      ! ! structure
   use grid_coms      , only : ngrids        ! ! intent(in)
   use consts_coms    , only : day_sec       ! ! intent(in)

   implicit none
   !----- Local variables -----------------------------------------------------------------!
   type(edtype)          , pointer :: cgrid
   character(len=str_len)          :: infile
   integer                         :: igr
   integer                         :: year_use
   integer                         :: ncyc
   integer                         :: iformat
   integer                         :: iv
   integer                         :: offset
   integer                         :: m2
   integer                         :: y2
   integer                         :: year_use_2
   logical                         :: exans
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
            select case (met_interp(iformat,iv))
            case (0,2,3,4)
               !----- If not, things are simple.  Just read in the month. -----------------!
               offset = 0
               call read_ol_file(infile,iformat,iv,year_use,mname(current_time%month)      &
                                ,current_time%year,offset,cgrid)
            case (1,5)
               !----- Here, just transfer future to current month.  -----------------------!
               call transfer_ol_month(trim(met_vars(iformat,iv)),met_frq(iformat,iv)       &
                                     ,cgrid)
            end select
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
            select case (met_interp(iformat,iv))
            case (1,5)
               offset = nint(day_sec / met_frq(iformat,iv)) * 31
               !----- Read the file. ------------------------------------------------------!
               call read_ol_file(infile,iformat,iv,year_use_2,mname(m2),y2,offset,cgrid)
            end select
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
  
   use ed_state_vars  , only : edtype            & ! structure
                             , polygontype       ! ! structure
   use met_driver_coms, only : met_frq           & ! intent(in)
                             , nformats          & ! intent(in)
                             , met_nv            & ! intent(in)
                             , met_interp        & ! intent(in)
                             , met_vars          & ! intent(in)
                             , have_co2          & ! intent(in)
                             , initial_co2       & ! intent(in)
                             , atm_tmp_intercept & ! intent(in)
                             , atm_tmp_slope     & ! intent(in)
                             , prec_intercept    & ! intent(in)
                             , prec_slope        & ! intent(in)
                             , humid_scenario    & ! intent(in)
                             , atm_rhv_min       & ! intent(in)
                             , atm_rhv_max       ! ! intent(in)
   use ed_misc_coms   , only : current_time      ! ! intent(in)
   use canopy_air_coms, only : ubmin             ! ! intent(in)
   use consts_coms    , only : rdry              & ! intent(in)
                             , cice              & ! intent(in)
                             , cliq              & ! intent(in)
                             , alli              & ! intent(in)
                             , rocp              & ! intent(in)
                             , p00               & ! intent(in)
                             , p00i              & ! intent(in)
                             , cp                & ! intent(in)
                             , cpi               & ! intent(in)
                             , day_sec           & ! intent(in)
                             , t00               & ! intent(in)
                             , t3ple             & ! intent(in)
                             , wdnsi             & ! intent(in)
                             , toodry            & ! intent(in)
                             , tsupercool        ! ! intent(in)
   use therm_lib      , only : rslif             & ! function
                             , ptrh2rvapil       & ! function
                             , thetaeiv          & ! function
                             , rehuil            & ! function
                             , qtk               ! ! function

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
   integer                    :: mnext
   integer                    :: mprev
   integer                    :: iformat
   integer                    :: iv
   integer                    :: ipy,isi
   real                       :: tprev
   real                       :: tnext
   real                       :: wnext
   real                       :: wprev

   real                       :: rvaux, rvsat, min_shv 

   real                       :: temp0
   real                       :: theta_prev
   real                       :: theta_next
   real                       :: relhum
   real                       :: snden !! snow density (kg/m3)
   real                       :: fice  !! precipication ice fraction
   real                       :: fliq,tsnow
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
         
         select case (met_interp(iformat,iv))
         case (0,2,3,4)
            !------------------------------------------------------------------------------!
            !      If we fall in this part of the if block, this means that the variable   !
            ! doesn't require time interpolation.                                          !
            !------------------------------------------------------------------------------!
            select case(met_interp(iformat,iv))
            case (2,4)
               !---------------------------------------------------------------------------!
               !     Time independent variables, set mprev to 1.                           !
               !---------------------------------------------------------------------------!
               mprev = 1

            case (0,3)
               
               !---------------------------------------------------------------------------!
               !     Time dependent variables, but with no time interpolation.             !
               !---------------------------------------------------------------------------!
               ndays_elapsed    = current_time%date - 1
               nseconds_elapsed = nint(current_time%time) + ndays_elapsed * nint(day_sec)
               mprev            = int(float(nseconds_elapsed) / met_frq(iformat,iv)) + 1
               
            end select
            
            !----- Find which variable it is, and then fill the polygons. -----------------!
            select case (trim(met_vars(iformat,iv)))
            case('nbdsf')   !----- Near IR beam downward shortwave flux. ------ [   W/m²] -!
               do ipy = 1,cgrid%npolygons
                  cgrid%met(ipy)%nir_beam = cgrid%metinput(ipy)%nbdsf(mprev)
               end do

            case('nddsf')   !----- Near IR diffuse downward shortwave flux. --- [   W/m²] -!
               do ipy = 1,cgrid%npolygons
                  cgrid%met(ipy)%nir_diffuse = cgrid%metinput(ipy)%nddsf(mprev)
               end do

            case('vbdsf')   !----- Visible beam downward shortwave flux. ------ [   W/m²] -!
               do ipy = 1,cgrid%npolygons
                  cgrid%met(ipy)%par_beam = cgrid%metinput(ipy)%vbdsf(mprev)
               end do

            case('vddsf')   !----- Visible diffuse downward shortwave flux. --- [   W/m²] -!
               do ipy = 1,cgrid%npolygons
                  cgrid%met(ipy)%par_diffuse = cgrid%metinput(ipy)%vddsf(mprev)
               end do

            case('prate')   !----- Precipitation rate. ------------------------ [kg/m²/s] -!
               do ipy = 1,cgrid%npolygons
                  cgrid%met(ipy)%pcpg = cgrid%metinput(ipy)%prate(mprev)
               end do

            case('dlwrf')   !----- Downwelling longwave radiation. ------------ [   W/m²] -!
               do ipy = 1,cgrid%npolygons
                  cgrid%met(ipy)%rlong = cgrid%metinput(ipy)%dlwrf(mprev)
               end do

            case('pres')    
               !---------------------------------------------------------------------------!
               !     Air pressure [   N/m²] and Exner function [J/kg/K].                   !
               !---------------------------------------------------------------------------!
               do ipy = 1,cgrid%npolygons
                  cgrid%met(ipy)%prss  = cgrid%metinput(ipy)%pres(mprev)
                  cgrid%met(ipy)%exner = cp * (p00i * cgrid%met(ipy)%prss)**rocp
               end do

            case('hgt')     !----- Air pressure. ------------------------------ [      m] -!
               do ipy = 1,cgrid%npolygons
                  cgrid%met(ipy)%geoht = cgrid%metinput(ipy)%hgt(mprev)
               end do

            case('ugrd')    !----- Zonal wind. -------------------------------- [    m/s] -!
               !---------------------------------------------------------------------------!
               !     Here we add the square of the zonal wind, so in the end vels will     !
               ! have (twice) the kinetic energy.  This way if we interpolate in the       !
               ! vertical for the multi-site case, we interpolate energy.                  !
               !---------------------------------------------------------------------------!
               do ipy = 1,cgrid%npolygons
                  cgrid%met(ipy)%vels = cgrid%met(ipy)%vels                                &
                                      + cgrid%metinput(ipy)%ugrd(mprev)**2
               end do

            case('vgrd')    !----- Meridional wind. --------------------------- [    m/s] -!
               !---------------------------------------------------------------------------!
               !     Here we add the square of the meridional wind, so in the end vels     !
               ! will have (twice) the kinetic energy.  This way if we interpolate in the  !
               ! vertical for the multi-site case, we interpolate energy.                  !
               !---------------------------------------------------------------------------!
               do ipy = 1,cgrid%npolygons
                  cgrid%met(ipy)%vels = cgrid%met(ipy)%vels                                &
                                      + cgrid%metinput(ipy)%vgrd(mprev)**2
               end do

            case('sh')      !----- Specific humidity. ------------------------- [kg/kg_a] -!
               !---------------------------------------------------------------------------!
               !     Here we will just get the specific humidity. But soon we will check   !
               ! whether this value is consistent with our thermodynamics (i.e., it will   !
               ! not be supersaturated).  If not we will impose the maximum.               !
               !---------------------------------------------------------------------------!
               do ipy = 1,cgrid%npolygons
                  cgrid%met(ipy)%atm_shv = cgrid%metinput(ipy)%sh(mprev)
               end do

            case('tmp')    
               !---------------------------------------------------------------------------!
               !     The flag is given at the air temperature, but we use the flag for     !
               ! potential temperature                                           [      K] !
               !---------------------------------------------------------------------------!
               do ipy = 1,cgrid%npolygons
                  cgrid%met(ipy)%atm_theta = cgrid%metinput(ipy)%tmp(mprev)                &
                                           * (p00 / cgrid%metinput(ipy)%pres(mprev))       &
                                           ** rocp
               end do

            case('co2')     !----- CO2 mixing ratio. -------------------------- [    ppm] -!
               do ipy = 1,cgrid%npolygons
                  cgrid%met(ipy)%atm_co2 = cgrid%metinput(ipy)%co2(mprev)
               end do
            end select
            
         !----- In this case, we need to interpolate. -------------------------------------!
         case (1,5)
            !------------------------------------------------------------------------------!
            !      If we fall in this part of the if block, this means that the variable   !
            ! must be interpolated over time.                                              !
            !---------------------------------------------------------------------------!

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
            
            
            !----- Find the number of seconds since the beginning of this month. ----------!
            ndays_elapsed    = current_time%date - 1
            nseconds_elapsed = nint(current_time%time) + ndays_elapsed * nint(day_sec)
            
            
            !------------------------------------------------------------------------------!
            ! mprev - Index of the element of the previous metinput data.                  !
            ! mnext - Index of the element of the next metinput data.                      !
            ! tprev = Time in seconds relative to the previous input time.                 !
            !------------------------------------------------------------------------------!
            mprev = int(float(nseconds_elapsed)/met_frq(iformat,iv)) + 1
            tprev = floor(float(nseconds_elapsed)/met_frq(iformat,iv))                     &
                  * met_frq(iformat,iv)
            mnext = mprev + 1
            tnext = tprev + met_frq(iformat,iv)

            !------------------------------------------------------------------------------!
            !     For the particular case where we at the last day of the month, after the !
            ! last met driver data of that month, the next time is the first day of the    !
            ! next month.  This will always be in the next position after day 31, even     !
            ! when the current month is shorter.                                           !
            !------------------------------------------------------------------------------!
            if(mnext > np)then
               mnext = 1 + nint(day_sec / met_frq(iformat,iv)) * 31
            endif
            
            !------------------------------------------------------------------------------!
            !      Get interpolation factors.                                              !
            !------------------------------------------------------------------------------!
            wnext = max(0.0,min(1.0,(real(nseconds_elapsed)-tprev)/met_frq(iformat,iv)))
            wprev = 1.0 - wnext

            !----- Find the variable and fill the sites. ----------------------------------!
            select case (trim(met_vars(iformat,iv)))
            case('prate')   !----- Precipitation rate. ------------------------ [kg/m²/s] -!
               do ipy = 1,cgrid%npolygons
                  cgrid%met(ipy)%pcpg = cgrid%metinput(ipy)%prate(mnext) * wnext           &
                                      + cgrid%metinput(ipy)%prate(mprev) * wprev
               end do

            case('dlwrf')   !----- Downwelling longwave radiation. ------------ [   W/m²] -!
               do ipy = 1,cgrid%npolygons
                  cgrid%met(ipy)%rlong = cgrid%metinput(ipy)%dlwrf(mnext) * wnext          &
                                       + cgrid%metinput(ipy)%dlwrf(mprev) * wprev
               end do

            case('pres')
               !---------------------------------------------------------------------------!
               !     Air pressure [   N/m²] and Exner function [J/kg/K].                   !
               !---------------------------------------------------------------------------!
               do ipy = 1,cgrid%npolygons
                  cgrid%met(ipy)%prss  = cgrid%metinput(ipy)%pres(mnext) * wnext           &
                                       + cgrid%metinput(ipy)%pres(mprev) * wprev
                  cgrid%met(ipy)%exner = cp * (p00i * cgrid%met(ipy)%prss)**rocp
               end do

            case('hgt')     !----- Air pressure. ------------------------------ [      m] -!
               do ipy = 1,cgrid%npolygons
                  cgrid%met(ipy)%geoht = cgrid%metinput(ipy)%hgt(mnext) * wnext            &
                                       + cgrid%metinput(ipy)%hgt(mprev) * wprev
               end do

            case('ugrd')    !----- Zonal wind. -------------------------------- [    m/s] -!
               !---------------------------------------------------------------------------!
               !     Here we add the square of the zonal wind, so in the end vels will     !
               ! have (twice) the kinetic energy.  This way if we interpolate in the       !
               ! vertical for the multi-site case, we interpolate energy.                  !
               !---------------------------------------------------------------------------!
               do ipy = 1,cgrid%npolygons
                  cgrid%met(ipy)%vels =  cgrid%met(ipy)%vels                               &
                                      +  ( cgrid%metinput(ipy)%ugrd(mnext) * wnext         &
                                         + cgrid%metinput(ipy)%ugrd(mprev) * wprev )**2
               end do

            case('vgrd')    !----- Meridional wind. --------------------------- [    m/s] -!
               !---------------------------------------------------------------------------!
               !     Here we add the square of the meridional wind, so in the end vels     !
               ! will have (twice) the kinetic energy.  This way if we interpolate in the  !
               ! vertical for the multi-site case, we interpolate energy.                  !
               !---------------------------------------------------------------------------!
               do ipy = 1,cgrid%npolygons
                  cgrid%met(ipy)%vels = cgrid%met(ipy)%vels                                &
                                      + ( cgrid%metinput(ipy)%vgrd(mnext) * wnext          &
                                        + cgrid%metinput(ipy)%vgrd(mprev) * wprev )**2
               enddo
            case('sh')      !----- Specific humidity. ------------------------- [kg/kg_a] -!
               do ipy = 1,cgrid%npolygons
                  cgrid%met(ipy)%atm_shv = cgrid%metinput(ipy)%sh(mnext) * wnext           &
                                         + cgrid%metinput(ipy)%sh(mprev) * wprev
               end do

            case('tmp')
               

               !---------------------------------------------------------------------------!
               !     The flag is given at the air temperature, but we use the flag for     !
               ! potential temperature                                           [      K] !
               !---------------------------------------------------------------------------!
               do ipy = 1,cgrid%npolygons

                  theta_next = cgrid%metinput(ipy)%tmp(mnext)                              &
                             * (p00 / cgrid%metinput(ipy)%pres(mnext))** rocp
                  theta_prev = cgrid%metinput(ipy)%tmp(mprev)                              &
                             * (p00 / cgrid%metinput(ipy)%pres(mprev))** rocp

                  !----- Interpolate potential temperature. -------------------------------!
                  cgrid%met(ipy)%atm_theta = theta_next * wnext + theta_prev * wprev
               end do

            case('nbdsf')   !----- Near IR beam downward shortwave flux. ------ [   W/m²] -!
               do ipy = 1,cgrid%npolygons
                  cgrid%met(ipy)%nir_beam = cgrid%metinput(ipy)%nbdsf(mnext) * wnext       &
                                          + cgrid%metinput(ipy)%nbdsf(mprev) * wprev
               end do

            case('nddsf')   !----- Near IR diffuse downward shortwave flux. --- [   W/m²] -!
               do ipy = 1,cgrid%npolygons
                  cgrid%met(ipy)%nir_diffuse = cgrid%metinput(ipy)%nddsf(mnext) * wnext    &
                                             + cgrid%metinput(ipy)%nddsf(mprev) * wprev
               end do

            case('vbdsf')   !----- Visible beam downward shortwave flux. ------ [   W/m²] -!
               do ipy = 1,cgrid%npolygons
                  cgrid%met(ipy)%par_beam = cgrid%metinput(ipy)%vbdsf(mnext) * wnext       &
                                          + cgrid%metinput(ipy)%vbdsf(mprev) * wprev
               end do

            case('vddsf')   !----- Visible diffuse downward shortwave flux. --- [   W/m²] -!
               do ipy = 1,cgrid%npolygons
                  cgrid%met(ipy)%par_diffuse = cgrid%metinput(ipy)%vddsf(mnext) * wnext    &
                                             + cgrid%metinput(ipy)%vddsf(mprev) * wprev
               enddo

            case('co2')     !----- CO2 mixing ratio. -------------------------- [    ppm] -!
               do ipy = 1,cgrid%npolygons
                  cgrid%met(ipy)%atm_co2 = cgrid%metinput(ipy)%co2(mnext) * wnext          &
                                         + cgrid%metinput(ipy)%co2(mprev) * wprev
               end do

            end select
         end select
      end do varloop
   end do formloop
   
   
   !---------------------------------------------------------------------------------------!
   !     Change from velocity squared to velocity, and compute qpcpg and dpcpg.            !
   !---------------------------------------------------------------------------------------!
   polyloop: do ipy = 1,cgrid%npolygons
         
      !----- CO2 --------------------------------------------------------------------------!
      if (.not.have_co2) cgrid%met(ipy)%atm_co2 = initial_co2

      !----- Adjust meteorological variables for simple climate scenarios. ----------------!
      cgrid%met(ipy)%atm_tmp      = cpi * cgrid%met(ipy)%atm_theta * cgrid%met(ipy)%exner
      temp0                       = cgrid%met(ipy)%atm_tmp
      
      if (atm_tmp_intercept /= 0.0 .or. atm_tmp_slope /= 1.0) then
         cgrid%met(ipy)%atm_tmp = atm_tmp_intercept                                        &
                                + cgrid%met(ipy)%atm_tmp * atm_tmp_slope
         !----- We must update potential temperature too. ---------------------------------!
         cgrid%met(ipy)%atm_theta = cp * cgrid%met(ipy)%atm_tmp / cgrid%met(ipy)%exner
      end if
      cgrid%met(ipy)%pcpg = max(0.0,prec_intercept + cgrid%met(ipy)%pcpg * prec_slope)
      !------------------------------------------------------------------------------------!


      if (humid_scenario == 1) then 
         !---------------------------------------------------------------------------------!
         !     Update atm_shv so the relative humidity remains the same.  We use the       !
         ! functions from therm_lib.f90 so it is consistent with the rest of the code.     !
         !---------------------------------------------------------------------------------!
         !----- 1. Temporarily convert specific humidity to mixing ratio. -----------------!
         rvaux  = cgrid%met(ipy)%atm_shv / (1. - cgrid%met(ipy)%atm_shv)
         !----- 2. Find relative humidity. ------------------------------------------------!
         relhum = rehuil(cgrid%met(ipy)%prss,temp0,rvaux)
         !----- 3. Find the equivalent relative humidity with the new temperature. --------!
         rvaux  = ptrh2rvapil(relhum,cgrid%met(ipy)%prss,cgrid%met(ipy)%atm_tmp)
         !----- 4. Convert the mixing ratio back to specific humidity. --------------------!
         cgrid%met(ipy)%atm_shv = rvaux / (1. + rvaux)
         !---------------------------------------------------------------------------------!
      end if


      !------------------------------------------------------------------------------------!
      !     Check the relative humidity associated with the current pressure, temperature, !
      ! and specific humidity.  Impose lower and upper bounds as prescribed by the vari-   !
      ! ables atm_rhv_min and atm_rhv_max (from met_driver_coms.f90, and defined at the    !
      ! init_met_params subroutine in ed_params.f90).                                      !
      !------------------------------------------------------------------------------------!
      rvaux  = cgrid%met(ipy)%atm_shv / (1. - cgrid%met(ipy)%atm_shv)
      relhum = rehuil(cgrid%met(ipy)%prss,cgrid%met(ipy)%atm_tmp,rvaux)
      !------------------------------------------------------------------------------------!
      !      Check whether the relative humidity is off-bounds.  If it is, then we re-     !
      ! calculate the mixing ratio and convert to specific humidity.                       !
      !------------------------------------------------------------------------------------!
      if (relhum < atm_rhv_min) then
         relhum = atm_rhv_min
         rvaux  = ptrh2rvapil(relhum,cgrid%met(ipy)%prss,cgrid%met(ipy)%atm_tmp)
         cgrid%met(ipy)%atm_shv = rvaux / (1. + rvaux)
      elseif (relhum > atm_rhv_max) then
         relhum = atm_rhv_max
         rvaux  = ptrh2rvapil(relhum,cgrid%met(ipy)%prss,cgrid%met(ipy)%atm_tmp)
         cgrid%met(ipy)%atm_shv = rvaux / (1. + rvaux)
      end if

      !------------------------------------------------------------------------------------!
      !    We now find the equivalent potential temperature.                               !
      !------------------------------------------------------------------------------------!
      cgrid%met(ipy)%atm_theiv = thetaeiv(cgrid%met(ipy)%atm_theta,cgrid%met(ipy)%prss     &
                                         ,cgrid%met(ipy)%atm_tmp,rvaux,rvaux,1)

      !------ Apply met to sites, and adjust met variables for topography. ----------------!
      call calc_met_lapse(cgrid,ipy)

      !----- Vels.  At this point vels is 2*Kinetic Energy, take the square root. ---------!
      cgrid%met(ipy)%vels = sqrt(max(0.0,cgrid%met(ipy)%vels))

      cpoly => cgrid%polygon(ipy)
      siteloop: do isi = 1,cpoly%nsites
         
         !----- Vels. The site level is also still in kinetic energy form. ----------------!
         cpoly%met(isi)%vels        = sqrt(max(0.0,cpoly%met(isi)%vels))
         cpoly%met(isi)%vels_stab   = max(ubmin,cpoly%met(isi)%vels)
         cpoly%met(isi)%vels_unstab = max(ubmin,cpoly%met(isi)%vels)
         
         !----- CO2.  In case we used the namelist, use that value. -----------------------!
         if (.not.have_co2) cpoly%met(isi)%atm_co2 = initial_co2

         
         !---------------------------------------------------------------------------------!
         !     We now find some derived properties.  In case several sites exist, the      !
         ! lapse rate was applied to pressure, temperature, and mixing ratio.  Then we     !
         ! calculate the Exner function, potential temperature and equivalent potential    !
         ! temperature, so it will respect the ideal gas law and first law of thermo-      !
         ! dynamics.                                                                       !
         !---------------------------------------------------------------------------------!
         cpoly%met(isi)%exner        = cp * (p00i * cpoly%met(isi)%prss) **rocp
         cpoly%met(isi)%atm_theta    = cp * cpoly%met(isi)%atm_tmp / cpoly%met(isi)%exner

         !---------------------------------------------------------------------------------!
         !     Check the relative humidity associated with the current pressure, temper-   !
         ! ature, and specific humidity.  Impose lower and upper bounds as prescribed by   !
         ! the variables atm_rhv_min and atm_rhv_max (from met_driver_coms.f90, and        !
         ! defined at the init_met_params subroutine in ed_params.f90).                    !
         !---------------------------------------------------------------------------------!
         rvaux  = cpoly%met(isi)%atm_shv / (1. - cpoly%met(isi)%atm_shv)
         relhum = rehuil(cpoly%met(isi)%prss,cpoly%met(isi)%atm_tmp,rvaux)
         !---------------------------------------------------------------------------------!
         !      Check whether the relative humidity is off-bounds.  If it is, then we re-  !
         ! calculate the mixing ratio and convert to specific humidity.                    !
         !---------------------------------------------------------------------------------!
         if (relhum < atm_rhv_min) then
            relhum = atm_rhv_min
            rvaux  = ptrh2rvapil(relhum,cpoly%met(isi)%prss,cpoly%met(isi)%atm_tmp)
            cpoly%met(isi)%atm_shv = rvaux / (1. + rvaux)
         elseif (relhum > atm_rhv_max) then
            relhum = atm_rhv_max
            rvaux  = ptrh2rvapil(relhum,cpoly%met(isi)%prss,cpoly%met(isi)%atm_tmp)
            cpoly%met(isi)%atm_shv = rvaux / (1. + rvaux)
         end if

         !----- Find the atmospheric equivalent potential temperature. --------------------!
         cpoly%met(isi)%atm_theiv = thetaeiv(cpoly%met(isi)%atm_theta,cpoly%met(isi)%prss  &
                                            ,cpoly%met(isi)%atm_tmp,rvaux,rvaux,2)

         !----- Solar radiation -----------------------------------------------------------!
         cpoly%met(isi)%rshort_diffuse = cpoly%met(isi)%par_diffuse                        &
                                       + cpoly%met(isi)%nir_diffuse
         cpoly%met(isi)%rshort         = cpoly%met(isi)%rshort_diffuse                     &
                                       + cpoly%met(isi)%par_beam + cpoly%met(isi)%nir_beam



         !---------------------------------------------------------------------------------!
         !  Precipitation internal energy and "depth".                                     !
         !  fcns derived from                                                              !
         !  Jin et al 1999 Hydrol Process. 13:2467-2482 Table 2                            !
         !  [[modified 11/16/09 by MCD]]                                                   !
         !---------------------------------------------------------------------------------!
         if (cpoly%met(isi)%atm_tmp > (t3ple + 2.5)) then
            !----- Rain only. -------------------------------------------------------------!
            fice = 0.0
            cpoly%met(isi)%dpcpg = max(0.0, cpoly%met(isi)%pcpg) *wdnsi

         elseif (cpoly%met(isi)%atm_tmp <= (t3ple + 2.5) .and.                             &
                 cpoly%met(isi)%atm_tmp  > (t3ple + 2.0) ) then
            !------------------------------------------------------------------------------!
            !     60% snow, 40% rain. (N.B. May not be appropriate for sub-tropical        !
            ! regions where the level of the melting layer is higher...).                  !
            !------------------------------------------------------------------------------!
            fice  = 0.6
            snden = 189.0
            cpoly%met(isi)%dpcpg = max(0.0, cpoly%met(isi)%pcpg)                           &
                                 * ((1.0-fice) * wdnsi + fice / snden)

         elseif (cpoly%met(isi)%atm_tmp <= (t3ple + 2.0) .and.                             &
                 cpoly%met(isi)%atm_tmp > t3ple          ) then
            !------------------------------------------------------------------------------!
            !     Increasing the fraction of snow. (N.B. May not be appropriate for        !
            ! sub-tropical regions where the level of the melting layer is higher...).     !
            !------------------------------------------------------------------------------!
            fice  = min(1.0, 1.+(54.62 - 0.2*cpoly%met(isi)%atm_tmp))
            snden = (50.0+1.7*(cpoly%met(isi)%atm_tmp-258.15)**1.5 )
            cpoly%met(isi)%dpcpg = max(0.0, cpoly%met(isi)%pcpg)                           &
                                 * ((1.0-fice) * wdnsi + fice / snden)

         elseif (cpoly%met(isi)%atm_tmp <= t3ple         .and.                             &
                 cpoly%met(isi)%atm_tmp > (t3ple - 15.0) ) then
            !----- Below freezing point, snow only. ---------------------------------------!
            fice  = 1.0
            snden = (50.0+1.7*(cpoly%met(isi)%atm_tmp-258.15)**1.5 )
            cpoly%met(isi)%dpcpg = max(0.0, cpoly%met(isi)%pcpg) / snden

         else ! if (copy%met(isi)%atm_tmp < (t3ple - 15.0)) then
            !----- Below freezing point, snow only. ---------------------------------------!
            fice  = 1.0
            snden = 50.
            cpoly%met(isi)%dpcpg = max(0.0, cpoly%met(isi)%pcpg) / snden

         end if
         !---------------------------------------------------------------------------------!

         !---------------------------------------------------------------------------------!
         !     Set internal energy.  This will be the precipitation times the specific     !
         ! internal energy of water (above or at triple point) multiplied by the liquid    !
         ! fraction plus the specific internal energy of ice (below or at the triple       !
         ! point) multiplied by the ice fraction.                                          !
         !---------------------------------------------------------------------------------!
         cpoly%met(isi)%qpcpg = max(0.0, cpoly%met(isi)%pcpg)                              &
                              * ( (1.0-fice) * cliq * ( max(t3ple,cpoly%met(isi)%atm_tmp)  &
                                                      - tsupercool)                        &
                                + fice * cice * min(cpoly%met(isi)%atm_tmp,t3ple))
         !---------------------------------------------------------------------------------!


      end do siteloop
           
   end do polyloop

   return
end subroutine update_met_drivers
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!      This subroutine reads in the meteorological data.  Subroutine shdf5_irec_f uses     !
! HDF5 routines to extract the data from the archive.  The rest of this code manages the   !
! allocation and storage of that data.  The extracted data may be from a year, other than  !
! the current simulation year.  This happens when we use a range of years to loop data     !
! over, typically during very long simulations.  In such cases, the number of days in      !
! February will change during leap years, and the number of days in February during the    !
! model year may not match that of the simulation year, therefore this must be taken into  !
! account.                                                                                 !
!==========================================================================================!
subroutine read_ol_file(infile,iformat, iv, year_use, mname, year, offset, cgrid)
   use ed_max_dims    , only : str_len       ! ! intent(in)
   use ed_state_vars  , only : edgrid_g      & ! structure
                             , edtype        & ! structure
                             , polygontype   ! ! structure
   use met_driver_coms, only : met_frq       & ! intent(in)
                             , met_nlon      & ! intent(in)
                             , met_nlat      & ! intent(in)
                             , met_vars      & ! intent(in)
                             , met_xmin      & ! intent(in)
                             , met_dx        & ! intent(in)
                             , met_ymin      & ! intent(in)
                             , met_dy        & ! intent(in)
                             , met_interp    & ! intent(in)
                             , no_ll         ! ! intent(in)
   use hdf5_utils     , only : shdf5_irec_f  & ! subroutine
                             , shdf5_info_f  ! ! subroutine
   use mem_polygons   , only : grid_type     ! ! intent(in)
   use consts_coms    , only : day_sec       ! ! intent(in)
  
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   character(len=str_len), intent(in)  :: infile
   type(edtype)          , target      :: cgrid
   integer               , intent(in)  :: iformat
   integer               , intent(in)  :: iv
   integer               , intent(in)  :: year_use
   character(len=3)      , intent(in)  :: mname
   integer               , intent(in)  :: year
   integer               , intent(in)  :: offset
   !----- Local variables. ----------------------------------------------------------------!
   integer                             :: ipy
   integer                             :: points_per_day
   integer                             :: nday
   integer                             :: nday_dset
   integer                             :: np
   integer                             :: ip
   integer                             :: np_dset
   real, dimension(:,:,:), allocatable :: metvar
   real, dimension(:,:,:), allocatable :: metvar_dset
   real, dimension(:)    , allocatable :: metscalar
   integer                             :: ndims
   integer, dimension(3)               :: idims
   integer                             :: ilon
   integer                             :: ilat
   integer                             :: ioa
   integer                             :: ioz
   integer                             :: d
   !----- External functions. -------------------------------------------------------------!
   logical, external                   :: isleap
   !---------------------------------------------------------------------------------------!


   !----- Geoposition variables do not apply. ---------------------------------------------!
   if (met_vars(iformat,iv).eq.'lat' .or. met_vars(iformat,iv).eq.'lon') return
  

   !---------------------------------------------------------------------------------------!
   !  Determine how many times to read.                                                    !
   !---------------------------------------------------------------------------------------!
  
   select case (met_interp(iformat,iv))
   case (2,4)
      !------------------------------------------------------------------------------------!
      !     In cases 2 and 4 we are dealing with constant time, only 1 point.              !
      !------------------------------------------------------------------------------------!
      np = 1
   case default

      !------------------------------------------------------------------------------------!
      !      Allocate for the number of points in the current year. In case this month is  !
      ! February, we must make sure that we account for the possibility of having a mis-   !
      ! match in the number of days in the model and in the dataset.                       !
      ! 1.  If the model year is leap but the dataset is not, we repeat February 28 twice. !
      ! 2.  If the model year is not leap but the dataset is, we drop February 29.         !
      !------------------------------------------------------------------------------------!

      points_per_day = nint(day_sec/met_frq(iformat,iv))

      !----- Determine the number of points necessary for the model year. -----------------!
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

      !!NEVER ASSUME ANYTHING ABOUT THE DATASET

         call shdf5_info_f(met_vars(iformat,iv),ndims,idims)
         nday_dset=idims(1)/points_per_day
 !        write (unit=*,fmt='(2(a,1x))') ' Reading Dataset : ', trim(met_vars(iformat,iv))
 !        write (unit=*,fmt='(2(a,i5))') ' for year        : ', year
 !        write (unit=*,fmt='(2(a,i5))') ' for dataset year: ', year_use
 !        write(unit=*,fmt='(a,1x,i5)')  ' NUMBER DAYS in FEB set to  = ', nday_dset
         

      end select 

      np      = nday      * points_per_day
      np_dset = nday_dset * points_per_day

   end select

   !----- Allocate the buffer space. ------------------------------------------------------!
  
   select case (met_interp(iformat,iv))
   case (3,5) 
      !------------------------------------------------------------------------------------!
      !     Cases 3 and 5, reading in one value representing the entire domain, though     !
      ! they change over time.                                                             !
      !------------------------------------------------------------------------------------!
      ndims = 1
      idims(1) = np_dset
      
      allocate(metvar(np,met_nlon(iformat),met_nlat(iformat)))
      allocate(metscalar(np_dset))
      call shdf5_irec_f(ndims,idims,trim(met_vars(iformat,iv)),rvara = metscalar)
      
      if (np < np_dset) then     
         !----- No interpolation is needed, but we need two buffers. ----------------------!
         do ip=1,np
            metvar(ip,:,:) = metscalar(ip)
         end do

      elseif(np == np_dset) then
         !----- Datasets are the same, no double allocation. ------------------------------!
         do ip=1,np
            metvar(ip,:,:) = metscalar(ip)
         end do

      else
         !---------------------------------------------------------------------------------!
         !     The dataset buffer is smaller than the models buffer, so we recycle the     !
         ! last day.                                                                       !
         !---------------------------------------------------------------------------------!
         do ip=1,np_dset
            metvar(ip,:,:) = metscalar(ip)
         end do

         do ip=np_dset+1,np
            !----- Then reuse the 28th for the 29th day. ----------------------------------!
            metvar(ip,:,:) = metscalar(ip-points_per_day+1)
         end do
      end if

      deallocate(metscalar)

   case (4)
      !------------------------------------------------------------------------------------!
      !      Case 4, we don't read, instead we use the value stored in met_freq as the     !
      ! scalar value.                                                                      !
      !------------------------------------------------------------------------------------!
      allocate (metvar(1,met_nlon(iformat),met_nlat(iformat)))
      metvar(1,:,:) = met_frq(iformat,iv)
      
   case (2)
      !----- Case 2, 2D grid is read in, ie. np = 1. --------------------------------------!

      allocate(metvar(1,met_nlon(iformat),met_nlat(iformat)))

      call shdf5_info_f(trim(met_vars(iformat,iv)),ndims,idims)
      
      !----- Check to see if the dimensions are consistent. -------------------------------!
      if (ndims /= 2 .or. idims(1) /= met_nlon(iformat)                                    &
                     .or. idims(2) /= met_nlat(iformat)) then
         write (unit=*,fmt='(2(a,1x))') ' In file: ',trim(infile)
         write (unit=*,fmt='(2(a,1x))') ' Dataset: ', trim(met_vars(iformat,iv))
         write (unit=*,fmt='(a)')       ' doesn''t have dimensions that match the'
         write (unit=*,fmt='(a)')       ' specified input, or latitude/longitude grid...'
         write (unit=*,fmt='(a,1x,i5)') ' MET_NLON = ',met_nlon(iformat)
         write (unit=*,fmt='(a,1x,i5)') ' MET_NLAT = ',met_nlat(iformat)
         write (unit=*,fmt='(a,1x,i5)') ' NDIMS    = ',ndims
         do d=1,ndims
            write(unit=*,fmt='(a,i5,a,1x,i5)') 'IDIMS (',d,') =', idims(d)
         end do
         call fatal_error('Mismatch between dataset and specified input'                   &
                         ,'read_ol_file','ed_met_driver.f90')
      end if

      call shdf5_irec_f(ndims,idims,trim(met_vars(iformat,iv)),rvara = metvar)

   case default
      !----- Case 1 or 0, we must read time-and-space variable array. ---------------------!
      call shdf5_info_f(trim(met_vars(iformat,iv)),ndims,idims)
      
      ! Check to see if the dimensions are consistent
      if (ndims /= 3 .or. idims(1).ne.np_dset .or. idims(2).ne.met_nlon(iformat)           &
                                              .or. idims(3).ne.met_nlat(iformat)) then
         print*,"DATASET: ",trim(met_vars(iformat,iv))
         print*,"DOES NOT HAVE DIMENSIONS THAT MATCH THE"
         print*,"SPECIFIED INPUT, OR LAT/LON GRID"
         print*,ndims
         print*,np_dset,met_nlon(iformat),met_nlat(iformat)
         print*,idims
         write (unit=*,fmt='(2(a,1x))') ' In file: ',trim(infile)
         write (unit=*,fmt='(2(a,1x))') ' Dataset: ', trim(met_vars(iformat,iv))
         write (unit=*,fmt='(a)')       ' doesn''t have dimensions that match the'
         write (unit=*,fmt='(a)')       ' specified input, or latitude/longitude grid...'
         write (unit=*,fmt='(a,1x,i5)') ' NP_DSET  = ', np_dset
         write (unit=*,fmt='(a,1x,i5)') ' MET_NLON = ', met_nlon(iformat)
         write (unit=*,fmt='(a,1x,i5)') ' MET_NLAT = ', met_nlat(iformat)
         write (unit=*,fmt='(a,1x,i5)') ' NDIMS    = ', ndims
         do d=1,ndims
            write(unit=*,fmt='(a,i5,a,1x,i5)') 'IDIMS (',d,') =', idims(d)
         end do
         call fatal_error('Mismatch between dataset and specified input' &
                         ,'read_ol_file','ed_met_driver.f90')
      end if

      if (np < np_dset) then
         !----- No interpolation is needed, but we need two buffers. ----------------------!
         allocate(metvar(np,met_nlon(iformat),met_nlat(iformat)))
         allocate(metvar_dset(np_dset,met_nlon(iformat),met_nlat(iformat)))
         
         call shdf5_irec_f(ndims, idims, trim(met_vars(iformat,iv)), rvara = metvar_dset)
         metvar(1:np,:,:) = metvar_dset(1:np,:,:)
         deallocate(metvar_dset)
         
         
      elseif (np == np_dset) then
         !----- Datasets are the same, no double allocation needed. -----------------------!
         allocate(metvar(np,met_nlon(iformat),met_nlat(iformat)))
         call shdf5_irec_f(ndims, idims, trim(met_vars(iformat,iv)), rvara = metvar)

      else
         !---------------------------------------------------------------------------------!
         !     The dataset buffer is smaller than the models buffer, so we recycle the     !
         ! last day.                                                                       !
         !---------------------------------------------------------------------------------!
         allocate(metvar(np,met_nlon(iformat),met_nlat(iformat)))
         allocate(metvar_dset(np_dset,met_nlon(iformat),met_nlat(iformat)))

         call shdf5_irec_f(ndims,idims,trim(met_vars(iformat,iv)),rvara=metvar_dset)
         metvar(1:np_dset,:,:) = metvar_dset(1:np_dset,:,:)
         
         !----- Then reuse the 28th for the 29th day. -------------------------------------!
         metvar(np_dset+1:np,:,:) = metvar_dset(np_dset-points_per_day+1:np_dset,:,:)
         deallocate(metvar_dset)
      end if
   end select

   !----- Define some aliases for the indices. --------------------------------------------!
   ioa = offset + 1
   ioz = offset + np

   !---------------------------------------------------------------------------------------!
   !      The dataset for that variable has been read into an array.  Now, the polygons    !
   ! are cycled through, and they are associated with an index in the grid.  If this is a  !
   ! lat-lon grid, the procedure is a simple interpolation between the domain endpoints,   !
   ! if this is a polar stereo-graphic projection, we make use of the associated indices   !
   ! saved in memory, ie. (cgrid%ilon(ipy)).                                               !
   !---------------------------------------------------------------------------------------!
   do ipy = 1,cgrid%npolygons
      
      !----- Get the indices.  Remember, latitude is flipped. -----------------------------!
      if (no_ll(iformat)) then
         ilon = min(max(1                                                                  &
                       ,1 + nint((cgrid%lon(ipy) - met_xmin(iformat)) / met_dx(iformat)))  &
                   ,met_nlon(iformat))
         ilat = met_nlat(iformat)                                                          &
              - min(max(1                                                                  &
                       ,1 + nint((cgrid%lat(ipy) - met_ymin(iformat)) / met_dy(iformat)))  &
                   , met_nlat(iformat)) + 1
      else
         ilon = cgrid%ilon(ipy)
         ilat = cgrid%ilat(ipy)
      end if
      
      !----- Get the time series. ---------------------------------------------------------!
      select case (trim(met_vars(iformat,iv)))
      case('nbdsf')
         cgrid%metinput(ipy)%nbdsf(ioa:ioz) = metvar(1:np,ilon,ilat)
      case('nddsf')
         cgrid%metinput(ipy)%nddsf(ioa:ioz) = metvar(1:np,ilon,ilat)
      case('vbdsf')
         cgrid%metinput(ipy)%vbdsf(ioa:ioz) = metvar(1:np,ilon,ilat)
      case('vddsf')
         cgrid%metinput(ipy)%vddsf(ioa:ioz) = metvar(1:np,ilon,ilat)
      case('prate')
         cgrid%metinput(ipy)%prate(ioa:ioz) = metvar(1:np,ilon,ilat)
      case('dlwrf')
         cgrid%metinput(ipy)%dlwrf(ioa:ioz) = metvar(1:np,ilon,ilat)
      case('pres')
         cgrid%metinput(ipy)%pres(ioa:ioz)  = metvar(1:np,ilon,ilat)
      case('hgt')
         cgrid%metinput(ipy)%hgt(ioa:ioz)   = metvar(1:np,ilon,ilat)
      case('ugrd')
         cgrid%metinput(ipy)%ugrd(ioa:ioz)  = metvar(1:np,ilon,ilat)
      case('vgrd')
         cgrid%metinput(ipy)%vgrd(ioa:ioz)  = metvar(1:np,ilon,ilat)
      case('sh')
         cgrid%metinput(ipy)%sh(ioa:ioz)    = metvar(1:np,ilon,ilat)
      case('tmp')
         cgrid%metinput(ipy)%tmp(ioa:ioz)   = metvar(1:np,ilon,ilat)
      case('co2')
         cgrid%metinput(ipy)%co2(ioa:ioz)   = metvar(1:np,ilon,ilat)
      end select
      
   end do
  
   !----- Deallocate the buffer. ----------------------------------------------------------!
   deallocate(metvar)
  
   return
end subroutine read_ol_file
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine transfers the dataset stored in the following month to the current   !
! month.                                                                                   !
!------------------------------------------------------------------------------------------!
subroutine transfer_ol_month(vname, frq, cgrid)
  
   use ed_state_vars, only : edtype   ! ! variable type
   use consts_coms  , only : day_sec  ! ! intent(in)

   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   type(edtype)    , target     :: cgrid     ! Current grid
   character(len=*), intent(in) :: vname     ! Variable name
   real            , intent(in) :: frq       ! Frequency of update
   !----- Local variables. ----------------------------------------------------------------!
   integer                      :: mem_size  ! Memory size
   integer                      :: ipy       ! Polygon counter
   integer                      :: ica       ! First element of current   month
   integer                      :: icz       ! Last  element of current   month
   integer                      :: ifa       ! First element of following month
   integer                      :: ifz       ! Last  element of following month
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !     The memory size is the size of each month.  The number of days is always 31       !
   ! because we allocate the maximum possible number of days.                              !
   !---------------------------------------------------------------------------------------!
   mem_size = nint(day_sec / frq) * 31

   !---------------------------------------------------------------------------------------!
   !     Compute the indices.                                                              !
   !---------------------------------------------------------------------------------------!
   ica = 1
   icz = mem_size
   ifa = icz + 1
   ifz = 2 * mem_size

   !---------------------------------------------------------------------------------------!
   !      Transfer the time series from following month to current month.                  !
   !---------------------------------------------------------------------------------------!
   select case (trim(vname))
   case ('nbdsf')
      do ipy=1,cgrid%npolygons
         cgrid%metinput(ipy)%nbdsf(ica:icz) = cgrid%metinput(ipy)%nbdsf(ifa:ifz)
      end do

   case ('nddsf')
      do ipy=1,cgrid%npolygons
         cgrid%metinput(ipy)%nddsf(ica:icz) = cgrid%metinput(ipy)%nddsf(ifa:ifz)
      end do

   case ('vbdsf')
      do ipy=1,cgrid%npolygons
         cgrid%metinput(ipy)%vbdsf(ica:icz) = cgrid%metinput(ipy)%vbdsf(ifa:ifz)
      end do

   case ('vddsf')
      do ipy=1,cgrid%npolygons
         cgrid%metinput(ipy)%vddsf(ica:icz) = cgrid%metinput(ipy)%vddsf(ifa:ifz)
      end do

   case ('prate')
      do ipy=1,cgrid%npolygons
         cgrid%metinput(ipy)%prate(ica:icz) = cgrid%metinput(ipy)%prate(ifa:ifz)
      end do

   case ('dlwrf')
      do ipy=1,cgrid%npolygons
         cgrid%metinput(ipy)%dlwrf(ica:icz) = cgrid%metinput(ipy)%dlwrf(ifa:ifz)
      end do

   case ('pres')
      do ipy=1,cgrid%npolygons
         cgrid%metinput(ipy)%pres(ica:icz)  = cgrid%metinput(ipy)%pres(ifa:ifz)
      end do

   case ('hgt')
      do ipy=1,cgrid%npolygons
         cgrid%metinput(ipy)%hgt(ica:icz)   = cgrid%metinput(ipy)%hgt(ifa:ifz)
      end do

   case ('ugrd')
      do ipy=1,cgrid%npolygons
         cgrid%metinput(ipy)%ugrd(ica:icz)  = cgrid%metinput(ipy)%ugrd(ifa:ifz)
      end do

   case ('vgrd')
      do ipy=1,cgrid%npolygons
         cgrid%metinput(ipy)%vgrd(ica:icz)  = cgrid%metinput(ipy)%vgrd(ifa:ifz)
      end do

   case ('sh')
      do ipy=1,cgrid%npolygons
         cgrid%metinput(ipy)%sh(ica:icz)    = cgrid%metinput(ipy)%sh(ifa:ifz)
      end do

   case ('tmp')
      do ipy=1,cgrid%npolygons
         cgrid%metinput(ipy)%tmp(ica:icz)   = cgrid%metinput(ipy)%tmp(ifa:ifz)
      end do

   case ('co2')
      do ipy=1,cgrid%npolygons
         cgrid%metinput(ipy)%co2(ica:icz)   = cgrid%metinput(ipy)%co2(ifa:ifz)
      end do

   end select

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


   polyloop: do ipy = 1,cgrid%npolygons
      min_dist = huge (1.)

      lonloop: do ilon = 1,nlon
         latloop: do ilat = 1,nlat

            if (lon(ilon,ilat) > 180.0) lon(ilon,ilat) = lon(ilon,ilat) - 360.0

            this_dist = dist_gc(cgrid%lon(ipy), lon(ilon,ilat)                             &
                               ,cgrid%lat(ipy), lat(ilon,ilat) )

            if(this_dist < min_dist) then
               cgrid%ilon(ipy) = ilon
               cgrid%ilat(ipy) = ilat
               min_dist        = this_dist
            end if

         end do latloop
      end do lonloop

   end do polyloop

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
   
   integer, intent(in) :: iformat   
   integer :: ndims
   integer :: d
   integer, dimension(3) :: idims
   type(edtype),target :: cgrid

   ! First check to see if there is lat/lon data in this dataset
   ! if the data exists, load it
   
   if(.not.no_ll(iformat)) then
         
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

