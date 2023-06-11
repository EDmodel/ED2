!==========================================================================================!
!==========================================================================================!
!     Module ed_met_driver.  This module contains routines that read, interpolate, and     !
! update meteorological information.                                                       !
!------------------------------------------------------------------------------------------!
module ed_met_driver
  contains

   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine is the main procedure that will take care of reading the meteoro- !
   ! logic data to be used during the run.                                                 !
   !---------------------------------------------------------------------------------------!
   subroutine read_met_driver_head()
      use ed_max_dims    , only : max_met_vars     ! ! intent(in)
      use met_driver_coms, only : nformats         & ! intent(in)
                                , ed_met_driver_db & ! intent(in)
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
                                , met_ll_header    & ! intent(out)
                                , met_land_mask    ! ! intent(out)
      implicit none  
      !----- Local variables --------------------------------------------------------------!
      logical :: l1
      logical :: yes_lon     ! Variable lon is provided in the HDF5 file
      logical :: yes_lat     ! Variable lat is provided in the HDF5 file
      integer :: iformat
      integer :: n
      !------------------------------------------------------------------------------------!


      !----- First thing, let's check whether the meterological data metafile exists. -----!
      inquire(file=trim(ed_met_driver_db),exist=l1)
      if (.not. l1) then
         write (unit=*,fmt='(a)') 'File '//trim(ed_met_driver_db)//' not found!'
         write (unit=*,fmt='(a)') 'Specify ED_MET_DRIVER_DB properly in ED namelist.'
         call fatal_error('Ed_met_driver_db not found!','read_met_driver_head'             &
                         &,'ed_met_driver.f90')
      end if
      !------------------------------------------------------------------------------------!


      !----- Loading the meterorological data metafile information. -----------------------!
      open(unit=12,file=trim(ed_met_driver_db),form='formatted',status='old')
      read(unit=12,fmt=*)  ! skip header

      !------ Read the number of different file formats. ----------------------------------!
      read(unit=12,fmt=*) nformats
      !------ Allocate the header information for each format -----------------------------!
      allocate(met_names    (nformats)              )
      allocate(met_nlon     (nformats)              )
      allocate(met_nlat     (nformats)              )
      allocate(met_dx       (nformats)              )
      allocate(met_dy       (nformats)              )
      allocate(met_xmin     (nformats)              )
      allocate(met_ymin     (nformats)              )
      allocate(met_nv       (nformats)              )
      allocate(met_vars     (nformats, max_met_vars))
      allocate(met_frq      (nformats, max_met_vars))
      allocate(met_interp   (nformats, max_met_vars))
      allocate(met_ll_header(nformats)              )
      allocate(met_land_mask(nformats)              )
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     By default this is true.  This means that ED will use the header information   !
      ! to obtain a regular longitude and latitude grid.                                   !
      !------------------------------------------------------------------------------------!
      met_ll_header(:) = .true.    
      !------------------------------------------------------------------------------------!

      !----- Go through all formats and read information. ---------------------------------!
      do iformat = 1,nformats
         read(unit=12,fmt='(a)')  met_names(iformat)
         read(unit=12,fmt=*)      met_nlon(iformat), met_nlat(iformat), met_dx(iformat)    &
                                , met_dy(iformat)  , met_xmin(iformat), met_ymin(iformat)
         read(unit=12,fmt=*)      met_nv(iformat)
         read(unit=12,fmt=*)      (met_vars(iformat,n)  ,n=1,met_nv(iformat))
         read(unit=12,fmt=*)      (met_frq(iformat,n)   ,n=1,met_nv(iformat))
         read(unit=12,fmt=*)      (met_interp(iformat,n),n=1,met_nv(iformat))
         
         !----- Make variable list case-insensitive. --------------------------------------!
         call tolower(met_vars(iformat,1:met_nv(iformat)),met_nv(iformat))
         !---------------------------------------------------------------------------------!

         
         !----- First check - find out whether lon/lat data are there. --------------------!
         yes_lon = any(met_vars(iformat,1:met_nv(iformat)) == 'lon')
         yes_lat = any(met_vars(iformat,1:met_nv(iformat)) == 'lat')
         !---------------------------------------------------------------------------------!


         !----- Check that either both lon and lat are provided, or none are provided. ----!
         if (yes_lat .and. yes_lon) then
            met_ll_header (iformat) = .false.
         elseif (yes_lat .neqv. yes_lon) then
            write (unit=*,fmt='(a,1x,l1)') ' Longitude in HDF5 : ',yes_lon
            write (unit=*,fmt='(a,1x,l1)') ' Longitude in HDF5 : ',yes_lat
            call fatal_error('Found only one coordinate. Either none or both should be'//  &
                             'present','read_met_driver_head','ed_met_driver.f90')
         end if
         !---------------------------------------------------------------------------------!


         !----- Find out whether a land/sea mask is provided. -----------------------------!
         met_land_mask(iformat) = any(met_vars(iformat,1:met_nv(iformat)) == 'land')
         !---------------------------------------------------------------------------------!

      end do
      !------------------------------------------------------------------------------------!

      close (unit=12,status='keep')

      return
   end subroutine read_met_driver_head
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine will initialise the arrays that will temporarily receive the       !
   ! meteorological dataset. It will nullify the pointers, allocate them, then assign a    !
   ! very strange number so in case the variable is not properly initialised the user will !
   ! notice something went wrong right away.                                               !
   !---------------------------------------------------------------------------------------!
   subroutine init_met_drivers
      use ed_max_dims     , only : max_met_vars
      use met_driver_coms , only : nformats          & ! intent(in)
                                 , met_nv            & ! intent(out)
                                 , met_vars          & ! intent(out)
                                 , met_frq           & ! intent(out)
                                 , met_interp        & ! intent(out)
                                 , has_co2           & ! intent(out)
                                 , has_ustar         ! ! intent(out)
      use canopy_air_coms , only : isfclyrm          ! ! intent(in)
      use ed_state_vars   , only : edgrid_g          & ! structure
                                 , edtype            & ! structure
                                 , polygontype       ! ! structure
      use grid_coms       , only : ngrids            ! ! intent(in)
      use consts_coms     , only : day_sec           ! ! intent(in)
      use lapse           , only : setLapseParms     ! ! sub-routine
      implicit none
      !----- Local variables --------------------------------------------------------------!
      type(polygontype), pointer :: cpoly
      type(edtype)     , pointer :: cgrid
      integer                    :: ipy,igr
      integer                    :: iformat
      integer                    :: iv
      integer                    :: mem_size
      !------------------------------------------------------------------------------------!


      !----- Set the lapse rates. ---------------------------------------------------------!
      do igr=1,ngrids
         call setLapseParms(edgrid_g(igr))
      end do

      !----- Read the information for each format. ----------------------------------------!
      has_co2   = .false.
      has_ustar = .false.
      formloop: do iformat = 1,nformats
         !----- Finding the met driver boundaries. ----------------------------------------!
         gridloop: do igr = 1,ngrids
            cgrid => edgrid_g(igr)

            polyloop: do ipy = 1,cgrid%npolygons

               cpoly => cgrid%polygon(ipy)
               
               !----- Loop over variables. ------------------------------------------------!
               varloop: do iv = 1,met_nv(iformat)
                  
                  !------------------------------------------------------------------------!
                  !  Either this is a constant or variable in time.                        !
                  !------------------------------------------------------------------------!
                  select case(met_interp(iformat,iv))
                  case(0,3)    !----- Variable. -------------------------------------------!
                     mem_size = nint(day_sec / met_frq(iformat,iv)) * 31
                  case(1,5)    !----- Interpolated.   Read two months in. -----------------!
                     mem_size = 2 * nint(day_sec / met_frq(iformat,iv)) * 31
                  case (2,4)  !----- Constant in time. ------------------------------------!
                     mem_size = 1
                  case default
                     call fatal_error('Invalid met_interp! It should be between 0 and 5.'  &
                                     ,'init_met_drivers','ed_met_driver.f90')
                  end select

                  !------------------------------------------------------------------------!
                  !     Allocate variables accordingly, and initialise them with a crazy   !
                  ! number to make the model crash in case it attempts to use them without !
                  ! assigning an appropriate number.                                       !
                  !------------------------------------------------------------------------!
                  select case (trim(met_vars(iformat,iv)))
                  case ('nbdsf')   !- Near IR beam downward shortwave flux. [   W/m2] -----!
                     nullify(cgrid%metinput(ipy)%nbdsf)
                     allocate(cgrid%metinput(ipy)%nbdsf(mem_size))
                     cgrid%metinput(ipy)%nbdsf = huge(1.)

                  case ('nddsf')   !- Near IR diffuse downward shortwave flux. [   W/m2] --!
                     nullify(cgrid%metinput(ipy)%nddsf)
                     allocate(cgrid%metinput(ipy)%nddsf(mem_size))
                     cgrid%metinput(ipy)%nddsf = huge(1.)

                  case ('vbdsf')   !- Visible beam downward shortwave flux. [   W/m2] -----!
                     nullify(cgrid%metinput(ipy)%vbdsf)
                     allocate(cgrid%metinput(ipy)%vbdsf(mem_size))
                     cgrid%metinput(ipy)%vbdsf = huge(1.)

                  case ('vddsf')   !- Visible diffuse downward shortwave flux. [   W/m2] --!
                     nullify(cgrid%metinput(ipy)%vddsf)
                     allocate(cgrid%metinput(ipy)%vddsf(mem_size))
                     cgrid%metinput(ipy)%vddsf = huge(1.)

                  case ('prate')   !- Precipitation rate. [kg/m2/s] -----------------------!
                     nullify(cgrid%metinput(ipy)%prate)
                     allocate(cgrid%metinput(ipy)%prate(mem_size))
                     cgrid%metinput(ipy)%prate = huge(1.)

                  case ('dlwrf')   !- Downwelling longwave radiation. [   W/m2] -----------!
                     nullify(cgrid%metinput(ipy)%dlwrf)
                     allocate(cgrid%metinput(ipy)%dlwrf(mem_size))
                     cgrid%metinput(ipy)%dlwrf = huge(1.)

                  case ('pres')   !----- Air pressure. [   N/m2] --------------------------!
                     nullify(cgrid%metinput(ipy)%pres)
                     allocate(cgrid%metinput(ipy)%pres(mem_size))
                     cgrid%metinput(ipy)%pres = huge(1.)

                  case ('hgt')     !----- Reference height [      m] ----------------------!
                     nullify(cgrid%metinput(ipy)%hgt)
                     allocate(cgrid%metinput(ipy)%hgt(mem_size))
                     cgrid%metinput(ipy)%hgt = huge(1.)

                  case ('ugrd')    !----- Zonal wind. [    m/s] ---------------------------!
                     nullify(cgrid%metinput(ipy)%ugrd)
                     allocate(cgrid%metinput(ipy)%ugrd(mem_size))
                     cgrid%metinput(ipy)%ugrd = huge(1.)

                  case ('vgrd')    !----- Meridional wind.[    m/s]  ----------------------!
                     nullify(cgrid%metinput(ipy)%vgrd)
                     allocate(cgrid%metinput(ipy)%vgrd(mem_size))
                     cgrid%metinput(ipy)%vgrd = huge(1.)

                  case ('sh')      !----- Specific humidity. [kg/kg_a] --------------------!
                     nullify(cgrid%metinput(ipy)%sh)
                     allocate(cgrid%metinput(ipy)%sh(mem_size))
                     cgrid%metinput(ipy)%sh = huge(1.)

                  case ('tmp')     !----- Air temperature. [      K] ----------------------!
                     nullify(cgrid%metinput(ipy)%tmp)
                     allocate(cgrid%metinput(ipy)%tmp(mem_size))
                     cgrid%metinput(ipy)%tmp = huge(1.)

                  case ('co2')     !----- CO2 mixing ratio. [umol/mol] --------------------!
                     has_co2 = .true.
                     nullify(cgrid%metinput(ipy)%co2)
                     allocate(cgrid%metinput(ipy)%co2(mem_size))
                     cgrid%metinput(ipy)%co2 = huge(1.)

                  case ('ustar')   !----- Friction velocity [    m/s]. --------------------!
                     has_ustar = .true.
                     nullify(cgrid%metinput(ipy)%atm_ustar)
                     allocate(cgrid%metinput(ipy)%atm_ustar(mem_size))
                     cgrid%metinput(ipy)%atm_ustar = huge(1.)

                  case ('lon','lat','land') !---- Skip coordinates and land/sea mask. -----!
                  case default
                     call fatal_error('Invalid met variable'//trim(met_vars(iformat,iv))   &
                                     //'!','init_met_drivers','ed_met_driver.f90')
                  end select
               end do varloop
            end do polyloop
         end do gridloop
      end do formloop


      !------------------------------------------------------------------------------------!
      !     Friction velocity is not necessary (and is deprecated except when developing   !
      ! or testing the code). In case the user wants to run the model with prescribed u*,  !
      ! then they must provide it.                                                         !
      !------------------------------------------------------------------------------------!
      if ( isfclyrm == 0 .and. .not. has_ustar) then
         call fatal_error('You must provide u* if you want to run with prescribed u*!!!'   &
                         ,'init_met_drivers','ed_met_driver.f90')
      end if 
      !------------------------------------------------------------------------------------!

      return
   end subroutine init_met_drivers
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine will control the HDF5 file reading.                               !
   !---------------------------------------------------------------------------------------!
   subroutine read_met_drivers_init
      use ed_max_dims       , only : str_len          ! ! intent(in)
      use ed_state_vars     , only : edgrid_g         & ! structure
                                   , edtype           & ! structure
                                   , polygontype      ! ! structure
      use hdf5_utils        , only : shdf5_open_f     & ! subroutine
                                   , shdf5_close_f    ! ! subroutine
      use met_driver_coms   , only : nformats         & ! intent(in)
                                   , ishuffle         & ! intent(in)
                                   , met_names        & ! intent(in)
                                   , met_vars         & ! intent(in)
                                   , met_nv           & ! intent(in)
                                   , met_interp       & ! intent(in)
                                   , met_frq          & ! intent(in)
                                   , metcyc1          & ! intent(in)
                                   , metcycf          & ! intent(in)
                                   , metyears         ! ! intent(inout)
      use ed_misc_coms      , only : current_time     & ! intent(in)
                                   , iyeara           & ! intent(in)
                                   , iyearz           & ! intent(in)
                                   , imonthz          ! ! intent(in)
      use grid_coms         , only : ngrids           ! ! intent(in)
      use consts_coms       , only : day_sec          ! ! intent(in)
      use ed_node_coms      , only : mynum            ! ! intent(in)
      use random_utils      , only : init_random_seed ! ! subroutine
      implicit none
      !----- Local variables --------------------------------------------------------------!
      type(edtype)          , pointer :: cgrid
      character(len=str_len)          :: infile
      integer                         :: igr
      integer                         :: year_cyc
      integer                         :: year_cyc_2
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
      real                            :: rndnow
      logical                         :: exans
      logical                         :: not_cycle_co2
      !----- Local constants --------------------------------------------------------------!
      character(len=3), dimension(12), parameter :: mname = (/ 'JAN', 'FEB', 'MAR', 'APR'  &
                                                             , 'MAY', 'JUN', 'JUL', 'AUG'  &
                                                             , 'SEP', 'OCT', 'NOV', 'DEC' /) 
      !----- External functions -----------------------------------------------------------!
      logical         , external  :: isleap
      !----- Variables to be saved. -------------------------------------------------------!
      logical         , save      :: first_time=.true.
      !------------------------------------------------------------------------------------!
    


      !------------------------------------------------------------------------------------!
      !    If this is the first time this subroutine is called, we build the meteoro-      !
      ! logical driver sequence of years.  Often is the case that the user has a meteoro-  !
      ! logical data set that is shorter than the time span he or she is intending to      !
      ! run. In this case we need to fill the other years with some data.                  !
      !------------------------------------------------------------------------------------!
      if (first_time) then
         first_time = .false.

         !---------------------------------------------------------------------------------!
         !     Define the number of years we have meteorological information available.    !
         ! In case the simulation is supposed to end in December, we must add an extra     !
         ! year, because the meteorological forcing may need to access the following       !
         ! month.                                                                          !
         !---------------------------------------------------------------------------------!
         ncyc   = metcycf - metcyc1 + 1
         if (imonthz == 12) then
            nyears = iyearz  - iyeara  + 2
         else
            nyears = iyearz  - iyeara  + 1
         end if

         !------ Allocating the sequence of years. ----------------------------------------!
         allocate(metyears(nyears))      
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !     Now we must decide how we are going to select the meteorological driver     !
         ! years, in particular, the years outside the meteorological driver range.        !
         !---------------------------------------------------------------------------------!
         select case (ishuffle)
         case (0)
            !------------------------------------------------------------------------------!
            !     Make a loop over the years, repeating the loop as many times as we need. !
            ! First we determine how many full cycles fit between iyeara and metcyc1. Note !
            ! that this can be negative in case metcyc1 is less than iyeara, but this will !
            ! work too.                                                                    !
            !------------------------------------------------------------------------------!
            nfullcyc = floor(real(metcyc1-iyeara)/real(ncyc))
            i1stfull = metcyc1 - nfullcyc * ncyc
            
            year_use   = metcycf - (i1stfull - iyeara)
            if (mynum == 1) then
               write (unit=*,fmt='(a)'       )  '------------------------------'
               write (unit=*,fmt='(a,1x,i12)')  ' - METCYC1 =',metcyc1
               write (unit=*,fmt='(a,1x,i12)')  ' - METCYCF =',metcycf
               write (unit=*,fmt='(a,1x,i12)')  ' - NYEARS  =',nyears
               write (unit=*,fmt='(a)'       )  '------------------------------'
               write (unit=*,fmt='(2(1x,a))') '       IYEAR','    YEAR_USE'
            end if
            do iyear=1,nyears
               if (year_use == metcycf) then
                  year_use = metcyc1
               else
                  year_use = year_use + 1
               end if
               metyears(iyear) = year_use
               if (mynum == 1) then
                  write (unit=*,fmt='(2(1x,i12))') iyear,metyears(iyear)
               end if
            end do
            if (mynum == 1) then
               write (unit=*,fmt='(a)'       )  '------------------------------'
               write (unit=*,fmt='(a)'       )  ' '
            end if

         !---------------------------------------------------------------------------------!
         !    For years outside the met driver range, we pick up a random year.  Because   !
         ! we use the actual sequence of random years given by the code random generator,  !
         ! this sequence will be always the same for a given met driver and first year.    !
         !---------------------------------------------------------------------------------!
         case (1)
            !----- Reset the seed for random number generation. ---------------------------!
            call random_seed()
            !------------------------------------------------------------------------------!

            if (mynum == 1) then
               write (unit=*,fmt='(a)'       )  '------------------------------'
               write (unit=*,fmt='(a,1x,i12)')  ' - METCYC1 =',metcyc1
               write (unit=*,fmt='(a,1x,i12)')  ' - METCYCF =',metcycf
               write (unit=*,fmt='(a,1x,i12)')  ' - NYEARS  =',nyears
               write (unit=*,fmt='(a)'       )  '------------------------------'
               write (unit=*,fmt='(2(1x,a))') '       IYEAR','    YEAR_USE'
            end if
            do year_use=iyeara,iyearz
               iyear=year_use-iyeara+1
               if (year_use < metcyc1 .or. year_use > metcycf) then
                  call random_number(rndnow)
                  metyears(iyear) = metcyc1 + mod(floor(real(ncyc)*rndnow),ncyc)
               else
                  metyears(iyear) = year_use
               end if
               if (mynum == 1) then
                  write (unit=*,fmt='(2(1x,i12))') iyear,metyears(iyear)
               end if
            end do
            if (mynum == 1) then
               write (unit=*,fmt='(a)'       )  '------------------------------'
               write (unit=*,fmt='(a)'       )  ' '
            end if
            !------------------------------------------------------------------------------!

            !------------------------------------------------------------------------------!
            !    Call the random seed initialisation that will make things really random.  !
            ! This is to avoid other routines that may use random numbers to get the same  !
            ! results every time the model runs.                                           !
            !------------------------------------------------------------------------------!
            call init_random_seed()
            !------------------------------------------------------------------------------!
         !---------------------------------------------------------------------------------!

         !---------------------------------------------------------------------------------!
         !    For years outside the met driver range, we pick up a random year.  Unlike    !
         ! option 1, this will use different random numbers every time the model is called !
         ! because we no longer reset random_seed.                                         !
         !---------------------------------------------------------------------------------!
         case (2)

            if (mynum == 1) then
               write (unit=*,fmt='(a)'       )  '------------------------------'
               write (unit=*,fmt='(a,1x,i12)')  ' - METCYC1 =',metcyc1
               write (unit=*,fmt='(a,1x,i12)')  ' - METCYCF =',metcycf
               write (unit=*,fmt='(a,1x,i12)')  ' - NYEARS  =',nyears
               write (unit=*,fmt='(a)'       )  '------------------------------'
               write (unit=*,fmt='(2(1x,a))') '       IYEAR','    YEAR_USE'
            end if
            do year_use=iyeara,iyearz
               iyear=year_use-iyeara+1
               if (year_use < metcyc1 .or. year_use > metcycf) then
                  call date_and_time(values=seedtime)
                  call random_number(rndnow)
                  metyears(iyear) = metcyc1 + mod(floor(real(ncyc)*rndnow),ncyc)
               else
                  metyears(iyear) = year_use
               end if
               if (mynum == 1) then
                  write (unit=*,fmt='(2(1x,i12))') iyear,metyears(iyear)
               end if
            end do
            if (mynum == 1) then
               write (unit=*,fmt='(a)'       )  '------------------------------'
               write (unit=*,fmt='(a)'       )  ' '
            end if
         end select
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     We now retrieve the met driver year based on the stored sequence.              !
      !------------------------------------------------------------------------------------!
      !----- If we need to recycle over years, find the appropriate year to apply. --------!
      year_cyc = current_time%year
      ncyc = metcycf - metcyc1 + 1
      !----- If we are after the last year... ---------------------------------------------!
      do while(year_cyc > metcycf)
         year_cyc = year_cyc - ncyc
      end do
      !----- If we are before the first year... -------------------------------------------!
      do while(year_cyc < metcyc1)
         year_cyc = year_cyc + ncyc
      end do
      !------------------------------------------------------------------------------------!

      gridloop: do igr = 1,ngrids

         cgrid => edgrid_g(igr)

         !----- Loop over the different file formats --------------------------------------!
         formloop: do iformat = 1, nformats

            !------------------------------------------------------------------------------!
            !   SPECIAL CASE FOR CO2:                                                      !
            !     Usually we do not want to cycle CO2 but only the other meteorology.      !
            !     Here, we overwrite year_use with current_time%year if CO2 is provided    !
            !     independently.                                                           !
            !------------------------------------------------------------------------------!
            not_cycle_co2 = (met_nv(iformat) == 1 .and. trim(met_vars(iformat,1)) == 'co2')
            if (not_cycle_co2) then
               !---------------------------------------------------------------------------!
               !      Make sure that the special CO2 file exists.  If not, then use the    !
               ! default met cycle years.                                                  !
               !---------------------------------------------------------------------------!
               write(infile,fmt='(a,i4.4,a,a)') trim(met_names(iformat)),current_time%year &
                                               ,mname(current_time%month),'.h5'
               inquire(file=trim(infile),exist=exans)
               if (exans) then
                  !----- Allow CO2 outside met cycle. -------------------------------------!
                  year_use = current_time%year
                  !------------------------------------------------------------------------!
               else
                  !----- Typical case. Pool data from the meteorological cycle. -----------!
                  year_use = year_cyc
                  !------------------------------------------------------------------------!
               end if
               !---------------------------------------------------------------------------!
            else
               !----- Typical case. Pool data from the meteorological cycle. --------------!
               year_use = year_cyc
               !---------------------------------------------------------------------------!
            end if
            !------------------------------------------------------------------------------!

            !----- Create the file name and check whether it exists. ----------------------!
            write(infile,fmt='(a,i4.4,a,a)')   trim(met_names(iformat)), year_use          &
                                              ,mname(current_time%month),'.h5'
            inquire(file=trim(infile),exist=exans)
            if (exans) then
               call shdf5_open_f(trim(infile),'R')
            else
               write (unit=*,fmt='(a)'       )  '------------------------------'
               write (unit=*,fmt='(a,1x,a)'  )  ' - MET_VARS    =',trim(met_vars(iformat,1))
               write (unit=*,fmt='(a,1x,l1)' )  ' - CYCLE_CO2   =',.not. not_cycle_co2
               write (unit=*,fmt='(a,1x,i12)')  ' - METCYC1     =',metcyc1
               write (unit=*,fmt='(a,1x,i12)')  ' - METCYCF     =',metcycf
               write (unit=*,fmt='(a,1x,i12)')  ' - NCYC        =',ncyc
               write (unit=*,fmt='(a,1x,i12)')  ' - NYEARS      =',nyears
               write (unit=*,fmt='(a,1x,i12)')  ' - IYEAR       =',iyear
               write (unit=*,fmt='(a,1x,i12)')  ' - MONTH_CURR  =',current_time%month
               write (unit=*,fmt='(a,1x,i12)')  ' - YEAR_CURR   =',current_time%year
               write (unit=*,fmt='(a,1x,i12)')  ' - YEAR_CYC    =',year_cyc
               write (unit=*,fmt='(a,1x,i12)')  ' - YEAR_USE    =',year_use
               write (unit=*,fmt='(a)'       )  '------------------------------'
               call fatal_error('Cannot open met driver input file '//trim(infile)//'!'    &
                               ,'read_met_drivers_init','ed_met_driver.f90')
            end if
            
            !------------------------------------------------------------------------------!
            !     Find the closest met driver.                                             !
            !------------------------------------------------------------------------------!
            call get_lonlat_land(cgrid,iformat)
            !------------------------------------------------------------------------------!

            !----- Loop over variables. and read the data. --------------------------------!
            do iv = 1, met_nv(iformat)
               offset = 0
               call read_ol_file(infile,iformat, iv, mname(current_time%month)             &
                                ,current_time%year, offset, cgrid)
            end do

            !-----  Close the HDF5 file. --------------------------------------------------!
            call shdf5_close_f()
            
            !------------------------------------------------------------------------------!
            !------------------------------------------------------------------------------!
            !      For all interpolated variables, we also need the next time.             !
            !------------------------------------------------------------------------------!

            !------------------------------------------------------------------------------!
            !     Find next month and year. If this takes us into the next year, increment !
            ! year and reset month to January.                                             !
            !------------------------------------------------------------------------------!
            m2 = current_time%month + 1
            select case (m2)
            case (13)
               m2         = 1
               y2         = current_time%year + 1
               year_cyc_2 = y2
               
               !----- If we are now after the last year... --------------------------------!
               do while(year_cyc_2 > metcycf)
                  year_cyc_2 = year_cyc_2 - ncyc
               end do
               !---------------------------------------------------------------------------!

               !----- If we are now before the first year... ------------------------------!
               do while(year_cyc_2 < metcyc1)
                  year_cyc_2 = year_cyc_2 + ncyc
               end do
               !---------------------------------------------------------------------------!
            case default
               !---- Same year as the previous month. -------------------------------------!
               y2         = current_time%year
               year_cyc_2 = year_cyc
               !---------------------------------------------------------------------------!
            end select

            !------------------------------------------------------------------------------#
            !   Again consider the special case of not cycling co2
            !------------------------------------------------------------------------------#
            if (not_cycle_co2) then
               !---------------------------------------------------------------------------!
               !      Make sure that the special CO2 file exists.  If not, then use the    !
               ! default met cycle years.                                                  !
               !---------------------------------------------------------------------------!
               write(infile,fmt='(a,i4.4,a,a)') trim(met_names(iformat)),y2,mname(m2),'.h5'
               inquire(file=trim(infile),exist=exans)
               if (exans) then
                  !----- Allow CO2 outside met cycle. -------------------------------------!
                  year_use_2 = y2
                  !------------------------------------------------------------------------!
               else
                  !----- Typical case. Pool data from the meteorological cycle. -----------!
                  year_use_2 = year_cyc_2
                  !------------------------------------------------------------------------!
               end if
               !---------------------------------------------------------------------------!
            else
               !----- Typical case. Pool data from the meteorological cycle. --------------!
               year_use_2 = year_cyc_2
               !---------------------------------------------------------------------------!
            end if
            !----- Now, open the file once. -----------------------------------------------!
            write(infile,fmt='(a,i4.4,a,a)')  trim(met_names(iformat)), year_use_2         &
                                             ,mname(m2),'.h5'
            inquire(file=trim(infile),exist=exans)
            if (exans) then
               call shdf5_open_f(trim(infile),'R')
            else
               write (unit=*,fmt='(a)'       )  '------------------------------'
               write (unit=*,fmt='(a,1x,a)'  )  ' - MET_VARS    =',trim(met_vars(iformat,1))
               write (unit=*,fmt='(a,1x,l1)' )  ' - CYCLE_CO2   =',.not. not_cycle_co2
               write (unit=*,fmt='(a,1x,i12)')  ' - METCYC1     =',metcyc1
               write (unit=*,fmt='(a,1x,i12)')  ' - METCYCF     =',metcycf
               write (unit=*,fmt='(a,1x,i12)')  ' - NCYC        =',ncyc
               write (unit=*,fmt='(a,1x,i12)')  ' - NYEARS      =',nyears
               write (unit=*,fmt='(a,1x,i12)')  ' - IYEAR       =',iyear
               write (unit=*,fmt='(a,1x,i12)')  ' - MONTH_CURR  =',current_time%month
               write (unit=*,fmt='(a,1x,i12)')  ' - YEAR_CURR   =',current_time%year
               write (unit=*,fmt='(a,1x,i12)')  ' - YEAR_CYC    =',year_cyc
               write (unit=*,fmt='(a,1x,i12)')  ' - YEAR_USE    =',year_use
               write (unit=*,fmt='(a,1x,i12)')  ' - MONTH_CYC_2 =',m2
               write (unit=*,fmt='(a,1x,i12)')  ' - YEAR_CYC_2  =',year_cyc_2
               write (unit=*,fmt='(a,1x,i12)')  ' - YEAR_USE_2  =',year_use_2
               write (unit=*,fmt='(a)'       )  '------------------------------'
               call fatal_error('Cannot open met driver input file '//trim(infile)//'!'    &
                               ,'read_met_drivers_init','ed_met_driver.f90')
            end if

            !----- Loop over variables. ---------------------------------------------------!
            varloop: do iv = 1, met_nv(iformat)

               select case (met_interp(iformat,iv))
               case (1,5)
                  offset = nint(day_sec / met_frq(iformat,iv)) * 31

                  !----- Read the file. ---------------------------------------------------!
                  call read_ol_file(infile,iformat,iv,mname(m2),y2,offset,cgrid)
               end select

            end do varloop
            
            !----- Close the HDF5 file. ---------------------------------------------------!
            call shdf5_close_f()

         end do formloop
      end do gridloop


      return
   end subroutine read_met_drivers_init
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   subroutine read_met_drivers
      use ed_max_dims    , only : str_len       ! ! intent(in)
      use ed_state_vars  , only : edgrid_g      & ! structure
                                , edtype        & ! structure
                                , polygontype   ! ! structure
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
      !----- Local variables --------------------------------------------------------------!
      type(edtype)          , pointer :: cgrid
      character(len=str_len)          :: infile
      integer                         :: igr
      integer                         :: year_cyc
      integer                         :: year_cyc_2
      integer                         :: year_use
      integer                         :: ncyc
      integer                         :: iformat
      integer                         :: iv
      integer                         :: offset
      integer                         :: m2
      integer                         :: y2
      integer                         :: year_use_2
      logical                         :: exans
      logical                         :: not_cycle_co2
      !----- Local constants --------------------------------------------------------------!
      character(len=3), dimension(12), parameter :: mname = (/ 'JAN', 'FEB', 'MAR', 'APR'  &
                                                             , 'MAY', 'JUN', 'JUL', 'AUG'  &
                                                             , 'SEP', 'OCT', 'NOV', 'DEC' /)
      !------------------------------------------------------------------------------------!




      !----- If we need to recycle over years, find the appropriate year to apply. --------!
      year_cyc = current_time%year
      ncyc     = metcycf - metcyc1 + 1

      !----- If we are after the last year... ---------------------------------------------!
      do while(year_cyc > metcycf)
         year_cyc = year_cyc - ncyc
      end do

      !----- If we are before the first year... -------------------------------------------!
      do while(year_cyc < metcyc1)
         year_cyc = year_cyc + ncyc
      end do

      gridloop: do igr=1,ngrids

         cgrid => edgrid_g(igr)

         !----- Loop over the different file formats --------------------------------------!
         formloop: do iformat = 1, nformats

            !------------------------------------------------------------------------------!
            !   SPECIAL CASE FOR CO2:                                                      !
            !     Usually we do not want to cycle CO2 but only the other meteorology.      !
            !     Here, we overwrite year_use with current_time%year if CO2 is provided    !
            !     independently.                                                           !
            !------------------------------------------------------------------------------!
            not_cycle_co2 = (met_nv(iformat) == 1 .and. trim(met_vars(iformat,1)) == 'co2')
            if (not_cycle_co2) then
               !---------------------------------------------------------------------------!
               !      Make sure that the special CO2 file exists.  If not, then use the    !
               ! default met cycle years.                                                  !
               !---------------------------------------------------------------------------!
               write(infile,fmt='(a,i4.4,a,a)') trim(met_names(iformat)),current_time%year &
                                               ,mname(current_time%month),'.h5'
               inquire(file=trim(infile),exist=exans)
               if (exans) then
                  !----- Allow CO2 outside met cycle. -------------------------------------!
                  year_use = current_time%year
                  !------------------------------------------------------------------------!
               else
                  !----- Typical case. Pool data from the meteorological cycle. -----------!
                  year_use = year_cyc
                  !------------------------------------------------------------------------!
               end if
               !---------------------------------------------------------------------------!
            else
               !----- Typical case. Pool data from the meteorological cycle. --------------!
               year_use = year_cyc
               !---------------------------------------------------------------------------!
            end if
            !------------------------------------------------------------------------------!



            !----- Create the file name and check whether it exists. ----------------------!
            write(infile,'(a,i4.4,a,a)')trim(met_names(iformat)), year_use,   &
                 mname(current_time%month),'.h5'

            inquire(file=trim(infile),exist=exans)
            if(exans)then
               call shdf5_open_f(trim(infile),'R')
            else
               write (unit=*,fmt='(a)'       )  '------------------------------'
               write (unit=*,fmt='(a,1x,a)'  )  ' - MET_VARS    =',trim(met_vars(iformat,1))
               write (unit=*,fmt='(a,1x,l1)' )  ' - CYCLE_CO2   =',.not. not_cycle_co2
               write (unit=*,fmt='(a,1x,i12)')  ' - METCYC1     =',metcyc1
               write (unit=*,fmt='(a,1x,i12)')  ' - METCYCF     =',metcycf
               write (unit=*,fmt='(a,1x,i12)')  ' - NCYC        =',ncyc
               write (unit=*,fmt='(a,1x,i12)')  ' - NYEARS      =',nyears
               write (unit=*,fmt='(a,1x,i12)')  ' - IYEAR       =',iyear
               write (unit=*,fmt='(a,1x,i12)')  ' - MONTH_CURR  =',current_time%month
               write (unit=*,fmt='(a,1x,i12)')  ' - YEAR_CURR   =',current_time%year
               write (unit=*,fmt='(a,1x,i12)')  ' - YEAR_CYC    =',year_cyc
               write (unit=*,fmt='(a,1x,i12)')  ' - YEAR_USE    =',year_use
               write (unit=*,fmt='(a)'       )  '------------------------------'
               call fatal_error('Cannot open met driver input file '//trim(infile)//'!'    &
                               ,'read_met_drivers','ed_met_driver.f90')
            end if
            
            !------------------------------------------------------------------------------!
            !     Find the closest met driver.                                             !
            !------------------------------------------------------------------------------!
            call get_lonlat_land(cgrid,iformat)
            !------------------------------------------------------------------------------!
            
            !----- Loop over variables. and read the data. --------------------------------!
            do iv = 1, met_nv(iformat)
               
               !----- Check whether this is an interpolation variable. --------------------!
               select case (met_interp(iformat,iv))
               case (0,2,3,4)
                  !----- If not, things are simple.  Just read in the month. --------------!
                  offset = 0
                  call read_ol_file(infile,iformat,iv,mname(current_time%month)      &
                                   ,current_time%year,offset,cgrid)
               case (1,5)
                  !----- Here, just transfer future to current month.  --------------------!
                  call transfer_ol_month(trim(met_vars(iformat,iv)),met_frq(iformat,iv)    &
                                        ,cgrid)
               end select
            end do

            !----- Close the HDF5 file. ---------------------------------------------------!
            call shdf5_close_f()

            !------------------------------------------------------------------------------!
            !------------------------------------------------------------------------------!
            !      For all interpolated variables, we also need the next time.             !
            !------------------------------------------------------------------------------!

            !------------------------------------------------------------------------------!
            !     Find next month and year. If this takes us into the next year, increment !
            ! year and reset month to January.                                             !
            !------------------------------------------------------------------------------!
            m2 = current_time%month + 1
            select case (m2)
            case (13)
               m2 = 1
               y2 = current_time%year + 1
               year_cyc_2 = y2
               
               !----- If we are now after the last year... --------------------------------!
               do while(year_cyc_2 > metcycf)
                  year_cyc_2 = year_cyc_2 - ncyc
               end do
               !---------------------------------------------------------------------------!

               !----- If we are now before the first year... ------------------------------!
               do while(year_cyc_2 < metcyc1)
                  year_cyc_2 = year_cyc_2 + ncyc
               end do
               !---------------------------------------------------------------------------!
            case default
               !---- Same year as the previous month. -------------------------------------!
               y2         = current_time%year
               year_cyc_2 = year_cyc
               !---------------------------------------------------------------------------!
            end select

            !---- Again, consider the special case for not_cycle_co2. ---------------------!
            if (not_cycle_co2) then
               !---------------------------------------------------------------------------!
               !      Make sure that the special CO2 file exists.  If not, then use the    !
               ! default met cycle years.                                                  !
               !---------------------------------------------------------------------------!
               write(infile,fmt='(a,i4.4,a,a)') trim(met_names(iformat)),y2                &
                                               ,mname(current_time%month),'.h5'
               inquire(file=trim(infile),exist=exans)
               if (exans) then
                  !----- Allow CO2 outside met cycle. -------------------------------------!
                  year_use_2 = y2
                  !------------------------------------------------------------------------!
               else
                  !----- Typical case. Pool data from the meteorological cycle. -----------!
                  year_use_2 = year_cyc_2
                  !------------------------------------------------------------------------!
               end if
               !---------------------------------------------------------------------------!
            else
               !----- Typical case. Pool data from the meteorological cycle. --------------!
               year_use = year_cyc_2
               !---------------------------------------------------------------------------!
            end if
            !------------------------------------------------------------------------------!


            !----- Now, open the file once. -----------------------------------------------!
            write(infile,fmt='(a,i4.4,a,a)')  trim(met_names(iformat)), year_use_2         &
                                             ,mname(m2),'.h5'
            inquire(file=trim(infile),exist=exans) 
            if (exans) then
               call shdf5_open_f(trim(infile),'R')
            else
               write (unit=*,fmt='(a)'       )  '------------------------------'
               write (unit=*,fmt='(a,1x,a)'  )  ' - MET_VARS    =',trim(met_vars(iformat,1))
               write (unit=*,fmt='(a,1x,l1)' )  ' - CYCLE_CO2   =',.not. not_cycle_co2
               write (unit=*,fmt='(a,1x,i12)')  ' - METCYC1     =',metcyc1
               write (unit=*,fmt='(a,1x,i12)')  ' - METCYCF     =',metcycf
               write (unit=*,fmt='(a,1x,i12)')  ' - NCYC        =',ncyc
               write (unit=*,fmt='(a,1x,i12)')  ' - NYEARS      =',nyears
               write (unit=*,fmt='(a,1x,i12)')  ' - IYEAR       =',iyear
               write (unit=*,fmt='(a,1x,i12)')  ' - MONTH_CURR  =',current_time%month
               write (unit=*,fmt='(a,1x,i12)')  ' - YEAR_CURR   =',current_time%year
               write (unit=*,fmt='(a,1x,i12)')  ' - YEAR_CYC    =',year_cyc
               write (unit=*,fmt='(a,1x,i12)')  ' - YEAR_USE    =',year_use
               write (unit=*,fmt='(a,1x,i12)')  ' - MONTH_CYC_2 =',m2
               write (unit=*,fmt='(a,1x,i12)')  ' - YEAR_CYC_2  =',year_cyc_2
               write (unit=*,fmt='(a,1x,i12)')  ' - YEAR_USE_2  =',year_use_2
               write (unit=*,fmt='(a)'       )  '------------------------------'
               call fatal_error ('Cannot open met driver input file '//trim(infile)//'!'   &
                                ,'read_met_drivers','ed_met_driver.f90')
            end if
         
            !----- Loop over variables. ---------------------------------------------------!
            do iv = 1, met_nv(iformat)
               select case (met_interp(iformat,iv))
               case (1,5)
                  offset = nint(day_sec / met_frq(iformat,iv)) * 31
                  !----- Read the file. ---------------------------------------------------!
                  call read_ol_file(infile,iformat,iv,mname(m2),y2,offset,cgrid)
               end select
            end do
            
            !----- Close the HDF5 file. ---------------------------------------------------!
            call shdf5_close_f()
         end do formloop
      end do gridloop


      return
   end subroutine read_met_drivers
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   subroutine update_met_drivers(cgrid)
      use ed_state_vars        , only : edtype                    & ! structure
                                      , polygontype               ! ! structure
      use met_driver_coms      , only : met_frq                   & ! intent(in)
                                      , nformats                  & ! intent(in)
                                      , imetavg                   & ! intent(in)
                                      , met_nv                    & ! intent(in)
                                      , met_interp                & ! intent(in)
                                      , met_vars                  & ! intent(in)
                                      , has_co2                   & ! intent(in)
                                      , has_ustar                 & ! intent(in)
                                      , initial_co2               & ! intent(in)
                                      , dt_radinterp              & ! intent(in)
                                      , atm_tmp_intercept         & ! intent(in)
                                      , atm_tmp_slope             & ! intent(in)
                                      , prec_intercept            & ! intent(in)
                                      , prec_slope                & ! intent(in)
                                      , humid_scenario            & ! intent(in)
                                      , atm_rhv_min               & ! intent(in)
                                      , atm_rhv_max               & ! intent(in)
                                      , print_radinterp           ! ! intent(in)
      use ed_misc_coms         , only : simtime                   & ! intent(in)
                                      , current_time              ! ! intent(in)
      use canopy_air_coms      , only : ubmin                     & ! intent(in)
                                      , ustmin                    ! ! intent(in)
      use canopy_radiation_coms, only : cosz_min                  ! ! intent(in)
      use consts_coms          , only : day_sec                   & ! intent(in)
                                      , t00                       & ! intent(in)
                                      , t3ple                     & ! intent(in)
                                      , wdnsi                     & ! intent(in)
                                      , toodry                    & ! intent(in)
                                      , pio180                    & ! intent(in)
                                      , tiny_num                  ! ! intent(in)
      use pft_coms             , only : include_pft               & ! intent(in)
                                      , hgt_max                   ! ! intent(in)
      use ed_max_dims          , only : n_pft                     ! ! intent(in)
      use therm_lib            , only : tl2uint                   & ! function
                                      , ptrh2rvapil               & ! function
                                      , press2exner               & ! function
                                      , extemp2theta              & ! function
                                      , thetaeiv                  & ! function
                                      , vpdefil                   & ! function
                                      , rehuil                    ! ! function
      use lapse                , only : calc_met_lapse            ! ! sub-routine
      use update_derived_utils , only : update_model_time_dm      ! ! sub-routine
      use radiate_utils        , only : ed_zen                    & ! function
                                      , mean_daysecz              & ! function
                                      , solar_radiation_breakdown & ! sub-routine
                                      , dump_radinfo              ! ! sub-routine

      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(edtype)     , target  :: cgrid
      !----- Local variables --------------------------------------------------------------!
      type(polygontype), pointer :: cpoly
      type(simtime)              :: prevmet_timea
      type(simtime)              :: nextmet_timea
      integer                    :: points_per_day
      integer                    :: nday
      integer                    :: np
      integer                    :: ndays_elapsed
      integer                    :: nseconds_elapsed
      integer                    :: mnext
      integer                    :: mprev
      integer                    :: iformat
      integer                    :: iv
      integer                    :: ipy
      integer                    :: isi
      integer                    :: ipft
      real(kind=4)               :: wnext
      real(kind=4)               :: wprev
      real(kind=4)               :: dtnext
      real(kind=4)               :: dtprev
      real(kind=4)               :: rvaux
      real(kind=4)               :: temp0
      real(kind=4)               :: relhum
      real(kind=4)               :: snden            ! snow density (kg/m3)
      real(kind=4)               :: fice             ! Ice fraction precipication
      real(kind=4)               :: fliq             ! Liquid fraction precipitation
      real(kind=4)               :: secz_prev        ! Mean of sec(zenith angle) - previous
      real(kind=4)               :: secz_next        ! Mean of sec(zenith angle) - next
      real(kind=4)               :: fperp_prev       ! Perpendicular flux - previous
      real(kind=4)               :: fperp_next       ! Perpendicular flux - next
      logical                    :: night_next
      logical                    :: night_prev
      !----- External functions -----------------------------------------------------------!
      logical         , external :: isleap       ! Check whether the year is leap or not.
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Flush vels to zero, so we can integrate by adding zonal and meridional winds.  !
      ! Also, make sure the cosine of the zenith angle is up to date.                      !
      !------------------------------------------------------------------------------------!
      do ipy=1,cgrid%npolygons
         cgrid%met(ipy)%vels = 0.0
         cgrid%cosz(ipy)     = ed_zen(cgrid%lon(ipy),cgrid%lat(ipy),current_time)
      end do
      !------------------------------------------------------------------------------------!


      formloop: do iformat = 1, nformats
         
         !----- Loop over variables -------------------------------------------------------!
         varloop: do iv = 1, met_nv(iformat)
            
            select case (met_interp(iformat,iv))
            case (0,2,3,4)
               !---------------------------------------------------------------------------!
               !      If we fall in this part of the if block, this means that the vari-   !
               ! able doesn't require time interpolation.                                  !
               !---------------------------------------------------------------------------!
               select case(met_interp(iformat,iv))
               case (2,4)
                  !------------------------------------------------------------------------!
                  !     Time independent variables, set mprev to 1 and the previous time   !
                  ! to the current time.                                                   !
                  !------------------------------------------------------------------------!
                  mprev = 1
                  prevmet_timea = current_time

               case (0,3)
                  
                  !------------------------------------------------------------------------!
                  !     Time dependent variables, but with no time interpolation.          !
                  !------------------------------------------------------------------------!
                  ndays_elapsed    = current_time%date - 1
                  nseconds_elapsed = nint(current_time%time) + ndays_elapsed * nint(day_sec)
                  mprev            = int(float(nseconds_elapsed) / met_frq(iformat,iv)) + 1
                  dtprev           = floor(real(nseconds_elapsed)/met_frq(iformat,iv))     &
                                   * met_frq(iformat,iv) - real(nseconds_elapsed)



                  !------------------------------------------------------------------------!
                  !     Find the actual time of the previous and the upcoming met          !
                  ! information.  We need this to find the integral of the secant of the   !
                  ! zenith angle, which is used to make a time interpolation that pre-     !
                  ! serves the diurnal cycle consistent.  This is particularly important   !
                  ! for tropical regions because the zenith angle changes very quickly and !
                  ! the radiation could be too skewed even with hourly data.               !
                  !------------------------------------------------------------------------!
                  prevmet_timea = current_time
                  nextmet_timea = current_time ! we won't use it, though.
                  select case (imetavg)
                  case (0)
                     !----- Instantaneous averages. ---------------------------------------!
                     call update_model_time_dm(prevmet_timea,dtprev)
                     !---------------------------------------------------------------------!

                  case (1)
                     !----- Averages ending at the reference time. ------------------------!
                     call update_model_time_dm(prevmet_timea,dtprev - met_frq(iformat,iv))
                     !---------------------------------------------------------------------!

                  case (2)
                     !----- Averages beginning at the reference time. ---------------------!
                     call update_model_time_dm(prevmet_timea,dtprev)
                     !---------------------------------------------------------------------!

                  case (3)
                     !----- Averages centred at the reference time. -----------------------!
                     call update_model_time_dm(prevmet_timea                               &
                                              ,dtprev - 0.5 * met_frq(iformat,iv))
                     !---------------------------------------------------------------------!

                  end select
                  !------------------------------------------------------------------------!
                  
               end select

               !----- Find which variable it is, and then fill the polygons. --------------!
               select case (trim(met_vars(iformat,iv)))
               case('nbdsf')   !----- Near IR beam downward shortwave flux. [   W/m2] -----!

                  !------------------------------------------------------------------------!
                  !    Decide whether to interpolate or not.                               !
                  !------------------------------------------------------------------------!
                  if ( imetavg == -1 .or. met_interp(iformat,iv) == 2 .or.                 &
                                          met_interp(iformat,iv) == 4 ) then
                     !----- Constant, do not attempt to do any interpolation. -------------!
                     do ipy = 1,cgrid%npolygons
                        cgrid%met(ipy)%nir_beam = cgrid%metinput(ipy)%nbdsf(mprev)
                     end do
                  else
                     do ipy = 1,cgrid%npolygons
                        !------------------------------------------------------------------!
                        !     Check whether this is day time or night time.  We should     !
                        ! only do the full interpolation thing during the day, the night   !
                        ! should have no incoming radiation until we incorporate fireflies !
                        ! (or a good twilight scheme) in the model.                        !
                        !------------------------------------------------------------------!
                        if (cgrid%cosz(ipy) > cosz_min) then
                           !---------------------------------------------------------------!
                           !     Define the normalisation factors for the previous and the !
                           ! next time.                                                    !
                           !---------------------------------------------------------------!
                           select case (imetavg)
                           case (0)
                              secz_prev = mean_daysecz(cgrid%lon(ipy),cgrid%lat(ipy)       &
                                                      ,prevmet_timea,dt_radinterp,0.)
                           case default
                              secz_prev = mean_daysecz(cgrid%lon(ipy),cgrid%lat(ipy)       &
                                                      ,prevmet_timea,dt_radinterp          &
                                                      ,met_frq(iformat,iv))
                           end select
                           !---------------------------------------------------------------!



                           !---------------------------------------------------------------!
                           !     Decide whether we want to use the previous or the next    !
                           ! data based on both the zenith angle.                          !
                           !---------------------------------------------------------------!
                           night_prev = secz_prev == 0.
                           !---------------------------------------------------------------!




                           !---------------------------------------------------------------!
                           !     Decide what to do based on the weighting factors.  This   !
                           ! is to avoid putting too little or too much energy at noon,    !
                           ! dawn, and dusk.                                               !
                           !---------------------------------------------------------------!
                           if (night_prev) then
                              !----- Middle of the night, set the value to zero. ----------!
                              fperp_prev = 0.
                              fperp_next = 0.
                              cgrid%met(ipy)%nir_beam = 0.
                           else
                              !------------------------------------------------------------!
                              !     Daytime: we scale the current time using the previous  !
                              ! secant and scaling back with the current cosine of zenith  !
                              ! angle.                                                     !
                              !------------------------------------------------------------!
                              fperp_prev = cgrid%metinput(ipy)%nbdsf(mprev) * secz_prev
                              cgrid%met(ipy)%nir_beam = fperp_prev * cgrid%cosz(ipy)
                           end if
                           !---------------------------------------------------------------!



                           !----- Assume all future stuff to be zero, we won't use them. --!
                           secz_next               = 0.
                           night_next              = .true.
                           fperp_next              = 0.
                           wprev                   = 1.0
                           wnext                   = 0.0
                           mnext                   = mprev
                           !---------------------------------------------------------------!

                        else
                           !----- Night time, assign it zero. -----------------------------!
                           secz_prev               = 0.
                           secz_next               = 0.
                           night_prev              = .true.
                           night_next              = .true.
                           fperp_prev              = 0.
                           fperp_next              = 0.
                           cgrid%met(ipy)%nir_beam = 0.0
                           wprev                   = 0.0
                           wnext                   = 0.0
                           mnext                   = mprev
                           !---------------------------------------------------------------!
                        end if
                        !------------------------------------------------------------------!



                        !------------------------------------------------------------------!
                        !     Print the debugging information if the user wants so.        !
                        !------------------------------------------------------------------!
                        if (print_radinterp) then
                           call dump_radinfo(cgrid,ipy,mprev,mnext,prevmet_timea           &
                                            ,nextmet_timea,met_frq(iformat,iv),wprev,wnext &
                                            ,secz_prev,secz_next,fperp_prev,fperp_next     &
                                            ,night_prev,night_next,'nbdsf')
                        end if
                        !------------------------------------------------------------------!
                     end do
                  end if
                  !------------------------------------------------------------------------!

               case('nddsf')   !----- Near IR diffuse downward shortwave flux. [   W/m2] --!

                  !------------------------------------------------------------------------!
                  !    Decide whether to interpolate or not.                               !
                  !------------------------------------------------------------------------!
                  if ( imetavg == -1 .or. met_interp(iformat,iv) == 2 .or.                 &
                                          met_interp(iformat,iv) == 4 ) then
                     !----- Constant, do not attempt to do any interpolation. -------------!
                     do ipy = 1,cgrid%npolygons
                        cgrid%met(ipy)%nir_diffuse = cgrid%metinput(ipy)%nddsf(mprev)
                     end do
                  else
                     do ipy = 1,cgrid%npolygons

                        !------------------------------------------------------------------!
                        !     Check whether this is day time or night time.  We should     !
                        ! only do the full interpolation thing during the day, the night   !
                        ! should have no incoming radiation until we incorporate fireflies !
                        ! (or a good twilight scheme) in the model.                        !
                        !------------------------------------------------------------------!
                        if (cgrid%cosz(ipy) > cosz_min) then
                           !---------------------------------------------------------------!
                           !     Define the normalisation factors for the previous and the !
                           ! next time.                                                    !
                           !---------------------------------------------------------------!
                           select case (imetavg)
                           case (0)
                              secz_prev = mean_daysecz(cgrid%lon(ipy),cgrid%lat(ipy)       &
                                                      ,prevmet_timea,dt_radinterp,0.)
                           case default
                              secz_prev = mean_daysecz(cgrid%lon(ipy),cgrid%lat(ipy)       &
                                                      ,prevmet_timea,dt_radinterp          &
                                                      ,met_frq(iformat,iv))
                           end select
                           !---------------------------------------------------------------!



                           !---------------------------------------------------------------!
                           !     Decide whether we want to use the previous or the next    !
                           ! data based on both the zenith angle.                          !
                           !---------------------------------------------------------------!
                           night_prev = secz_prev == 0.
                           !---------------------------------------------------------------!




                           !---------------------------------------------------------------!
                           !     Decide what to do based on the weighting factors.  This   !
                           ! is to avoid putting too little or too much energy at noon,    !
                           ! dawn, and dusk.                                               !
                           !---------------------------------------------------------------!
                           if (night_prev) then
                              !----- Middle of the night, set the value to zero. ----------!
                              fperp_prev = 0.
                              fperp_next = 0.
                              cgrid%met(ipy)%nir_diffuse = 0.
                           else
                              !------------------------------------------------------------!
                              !     Daytime: we scale the current time using the previous  !
                              ! secant and scaling back with the current cosine of zenith  !
                              ! angle.                                                     !
                              !------------------------------------------------------------!
                              fperp_prev = cgrid%metinput(ipy)%nddsf(mprev) * secz_prev
                              cgrid%met(ipy)%nir_diffuse = fperp_prev * cgrid%cosz(ipy)
                           end if
                           !---------------------------------------------------------------!



                           !----- Assume all future stuff to be zero, we won't use them. --!
                           secz_next               = 0.
                           night_next              = .true.
                           fperp_next              = 0.
                           wprev                   = 1.0
                           wnext                   = 0.0
                           mnext                   = mprev
                           !---------------------------------------------------------------!

                        else
                           !----- Night time, assign it zero. -----------------------------!
                           secz_prev                  = 0.
                           secz_next                  = 0.
                           night_prev                 = .true.
                           night_next                 = .true.
                           fperp_prev                 = 0.
                           fperp_next                 = 0.
                           cgrid%met(ipy)%nir_diffuse = 0.0
                           wprev                      = 0.0
                           wnext                      = 0.0
                           mnext                      = mprev
                           !---------------------------------------------------------------!
                        end if
                        !------------------------------------------------------------------!



                        !------------------------------------------------------------------!
                        !     Print the debugging information if the user wants so.        !
                        !------------------------------------------------------------------!
                        if (print_radinterp) then
                           call dump_radinfo(cgrid,ipy,mprev,mnext,prevmet_timea           &
                                            ,nextmet_timea,met_frq(iformat,iv),wprev,wnext &
                                            ,secz_prev,secz_next,fperp_prev,fperp_next     &
                                            ,night_prev,night_next,'nddsf')
                        end if
                        !------------------------------------------------------------------!
                     end do
                  end if
                  !------------------------------------------------------------------------!

               case('vbdsf')   !----- Visible beam downward shortwave flux. [   W/m2] -----!
                  !------------------------------------------------------------------------!
                  !    Decide whether to interpolate or not.                               !
                  !------------------------------------------------------------------------!
                  if ( imetavg == -1 .or. met_interp(iformat,iv) == 2 .or.                 &
                                          met_interp(iformat,iv) == 4 ) then
                     !----- Constant, do not attempt to do any interpolation. -------------!
                     do ipy = 1,cgrid%npolygons
                        cgrid%met(ipy)%par_beam = cgrid%metinput(ipy)%vbdsf(mprev)
                     end do
                  else
                     do ipy = 1,cgrid%npolygons

                        !------------------------------------------------------------------!
                        !     Check whether this is day time or night time.  We should     !
                        ! only do the full interpolation thing during the day, the night   !
                        ! should have no incoming radiation until we incorporate fireflies !
                        ! (or a good twilight scheme) in the model.                        !
                        !------------------------------------------------------------------!
                        if (cgrid%cosz(ipy) > cosz_min) then
                           !---------------------------------------------------------------!
                           !     Define the normalisation factors for the previous and the !
                           ! next time.                                                    !
                           !---------------------------------------------------------------!
                           select case (imetavg)
                           case (0)
                              secz_prev = mean_daysecz(cgrid%lon(ipy),cgrid%lat(ipy)       &
                                                      ,prevmet_timea,dt_radinterp,0.)
                           case default
                              secz_prev = mean_daysecz(cgrid%lon(ipy),cgrid%lat(ipy)       &
                                                      ,prevmet_timea,dt_radinterp          &
                                                      ,met_frq(iformat,iv))
                           end select
                           !---------------------------------------------------------------!



                           !---------------------------------------------------------------!
                           !     Decide whether we want to use the previous or the next    !
                           ! data based on both the zenith angle.                          !
                           !---------------------------------------------------------------!
                           night_prev = secz_prev == 0.
                           !---------------------------------------------------------------!




                           !---------------------------------------------------------------!
                           !     Decide what to do based on the weighting factors.  This   !
                           ! is to avoid putting too little or too much energy at noon,    !
                           ! dawn,and dusk.                                                !
                           !---------------------------------------------------------------!
                           if (night_prev) then
                              !----- Middle of the night, set the value to zero. ----------!
                              fperp_prev = 0.
                              fperp_next = 0.
                              cgrid%met(ipy)%par_beam = 0.
                           else
                              !------------------------------------------------------------!
                              !     Daytime: we scale the current time using the previous  !
                              ! secant and scaling back with the current cosine of zenith  !
                              ! angle.                                                     !
                              !------------------------------------------------------------!
                              fperp_prev = cgrid%metinput(ipy)%vbdsf(mprev) * secz_prev
                              cgrid%met(ipy)%par_beam = fperp_prev * cgrid%cosz(ipy)
                           end if
                           !---------------------------------------------------------------!



                           !----- Assume all future stuff to be zero, we won't use them. --!
                           secz_next               = 0.
                           night_next              = .true.
                           fperp_next              = 0.
                           wprev                   = 1.0
                           wnext                   = 0.0
                           mnext                   = mprev
                           !---------------------------------------------------------------!

                        else
                           !----- Night time, assign it zero. -----------------------------!
                           secz_prev               = 0.
                           secz_next               = 0.
                           night_prev              = .true.
                           night_next              = .true.
                           fperp_prev              = 0.
                           fperp_next              = 0.
                           cgrid%met(ipy)%par_beam = 0.0
                           wprev                   = 0.0
                           wnext                   = 0.0
                           mnext                   = mprev
                           !---------------------------------------------------------------!
                        end if
                        !------------------------------------------------------------------!



                        !------------------------------------------------------------------!
                        !     Print the debugging information if the user wants so.        !
                        !------------------------------------------------------------------!
                        if (print_radinterp) then
                           call dump_radinfo(cgrid,ipy,mprev,mnext,prevmet_timea           &
                                            ,nextmet_timea,met_frq(iformat,iv),wprev,wnext &
                                            ,secz_prev,secz_next,fperp_prev,fperp_next     &
                                            ,night_prev,night_next,'vbdsf')
                        end if
                        !------------------------------------------------------------------!
                     end do
                  end if
                  !------------------------------------------------------------------------!

               case('vddsf')   !----- Visible diffuse downward shortwave flux. [   W/m2] --!

                  !------------------------------------------------------------------------!
                  !    Decide whether to interpolate or not.                               !
                  !------------------------------------------------------------------------!
                  if ( imetavg == -1 .or. met_interp(iformat,iv) == 2 .or.                 &
                                          met_interp(iformat,iv) == 4 ) then
                     !----- Constant, do not attempt to do any interpolation. -------------!
                     do ipy = 1,cgrid%npolygons
                        cgrid%met(ipy)%par_diffuse = cgrid%metinput(ipy)%vddsf(mprev)
                     end do
                  else
                     do ipy = 1,cgrid%npolygons

                        !------------------------------------------------------------------!
                        !     Check whether this is day time or night time.  We should     !
                        ! only do the full interpolation thing during the day, the night   !
                        ! should have no incoming radiation until we incorporate fireflies !
                        ! (or a good twilight scheme) in the model.                        !
                        !------------------------------------------------------------------!
                        if (cgrid%cosz(ipy) > cosz_min) then
                           !---------------------------------------------------------------!
                           !     Define the normalisation factors for the previous and the !
                           ! next time.                                                    !
                           !---------------------------------------------------------------!
                           select case (imetavg)
                           case (0)
                              secz_prev = mean_daysecz(cgrid%lon(ipy),cgrid%lat(ipy)       &
                                                      ,prevmet_timea,dt_radinterp,0.)
                           case default
                              secz_prev = mean_daysecz(cgrid%lon(ipy),cgrid%lat(ipy)       &
                                                      ,prevmet_timea,dt_radinterp          &
                                                      ,met_frq(iformat,iv))
                           end select
                           !---------------------------------------------------------------!



                           !---------------------------------------------------------------!
                           !     Decide whether we want to use the previous or the next    !
                           ! data based on both the zenith angle.                          !
                           !---------------------------------------------------------------!
                           night_prev = secz_prev == 0.
                           !---------------------------------------------------------------!




                           !---------------------------------------------------------------!
                           !     Decide what to do based on the weighting factors.  This   !
                           ! is to avoid putting too little or too much energy at noon,    !
                           ! dawn, and dusk.                                               !
                           !---------------------------------------------------------------!
                           if (night_prev) then
                              !----- Middle of the night, set the value to zero. ----------!
                              fperp_prev = 0.
                              fperp_next = 0.
                              cgrid%met(ipy)%par_diffuse = 0.
                           else
                              !------------------------------------------------------------!
                              !     Daytime: we scale the current time using the previous  !
                              ! secant and scaling back with the current cosine of zenith  !
                              ! angle.                                                     !
                              !------------------------------------------------------------!
                              fperp_prev = cgrid%metinput(ipy)%vddsf(mprev) * secz_prev
                              cgrid%met(ipy)%par_diffuse = fperp_prev * cgrid%cosz(ipy)
                           end if
                           !---------------------------------------------------------------!



                           !----- Assume all future stuff to be zero, we won't use them. --!
                           secz_next               = 0.
                           night_next              = .true.
                           fperp_next              = 0.
                           wprev                   = 1.0
                           wnext                   = 0.0
                           mnext                   = mprev
                           !---------------------------------------------------------------!

                        else
                           !----- Night time, assign it zero. -----------------------------!
                           secz_prev                  = 0.
                           secz_next                  = 0.
                           night_prev                 = .true.
                           night_next                 = .true.
                           fperp_prev                 = 0.
                           fperp_next                 = 0.
                           cgrid%met(ipy)%par_diffuse = 0.0
                           wprev                      = 0.0
                           wnext                      = 0.0
                           mnext                      = mprev
                           !---------------------------------------------------------------!
                        end if
                        !------------------------------------------------------------------!



                        !------------------------------------------------------------------!
                        !     Print the debugging information if the user wants so.        !
                        !------------------------------------------------------------------!
                        if (print_radinterp) then
                           call dump_radinfo(cgrid,ipy,mprev,mnext,prevmet_timea           &
                                            ,nextmet_timea,met_frq(iformat,iv),wprev,wnext &
                                            ,secz_prev,secz_next,fperp_prev,fperp_next     &
                                            ,night_prev,night_next,'vddsf')
                        end if
                        !------------------------------------------------------------------!
                     end do
                  end if
                  !------------------------------------------------------------------------!

               case('prate')   !----- Precipitation rate. [kg/m2/s] -----------------------!
                  do ipy = 1,cgrid%npolygons
                     cgrid%met(ipy)%pcpg = cgrid%metinput(ipy)%prate(mprev)
                  end do

               case('dlwrf')   !----- Downwelling longwave radiation. [   W/m2] -----------!
                  do ipy = 1,cgrid%npolygons
                     cgrid%met(ipy)%rlong = cgrid%metinput(ipy)%dlwrf(mprev)
                  end do

               case('pres')    
                  !------------------------------------------------------------------------!
                  !     Air pressure [   N/m2] and Exner function [J/kg/K].                !
                  !------------------------------------------------------------------------!
                  do ipy = 1,cgrid%npolygons
                     cgrid%met(ipy)%prss  = cgrid%metinput(ipy)%pres(mprev)
                  end do

               case('hgt')     !----- Air pressure. [      m] -----------------------------!
                  do ipy = 1,cgrid%npolygons
                     cgrid%met(ipy)%geoht = cgrid%metinput(ipy)%hgt(mprev)

                     !---------------------------------------------------------------------!
                     !    The reference height must be always above the top cohort.  We    !
                     ! never know when this may happen, but if there is a risk of the      !
                     ! reference height becoming lower than the top cohort, stop the run   !
                     ! and explain it to the user.                                         !
                     !---------------------------------------------------------------------!
                     do ipft = 1,n_pft
                        if (include_pft(ipft) .and. cgrid%met(ipy)%geoht <= hgt_max(ipft)) &
                        then
                           write (unit=*,fmt='(a)') '-------------------------------------'
                           write(unit=*,fmt='(a,1x,i12)')                                  &
                                                 ' - Polygon        =',ipy
                           write(unit=*,fmt='(a,1x,i12)')                                  &
                                                 ' - PFT            =',ipft
                           write (unit=*,fmt='(a,1x,es12.5)')                              &
                                                 ' - Ref. Height    =',cgrid%met(ipy)%geoht
                           write (unit=*,fmt='(a,1x,es12.5)')                              &
                                                 ' - PFT max height =',hgt_max(ipft)
                           write (unit=*,fmt='(a)') '-------------------------------------'
                           call fatal_error('Reference height must be higher than PFT '//  &
                                            'maximum height!!','update_met_drivers'        &
                                           ,'ed_met_driver.f90')
                        end if
                     end do
                     !---------------------------------------------------------------------!
                  end do

               case('ugrd')    !----- Zonal wind. [    m/s] -------------------------------!
                  !------------------------------------------------------------------------!
                  !     Here we add the square of the zonal wind, so in the end vels will  !
                  ! have (twice) the kinetic energy.  This way if we interpolate in the    !
                  ! vertical for the multi-site case, we interpolate energy.               !
                  !------------------------------------------------------------------------!
                  do ipy = 1,cgrid%npolygons
                     cgrid%met(ipy)%vels = cgrid%met(ipy)%vels                             &
                                         + cgrid%metinput(ipy)%ugrd(mprev)**2
                  end do

               case('vgrd')    !----- Meridional wind. [    m/s] --------------------------!
                  !------------------------------------------------------------------------!
                  !     Here we add the square of the meridional wind, so in the end vels  !
                  ! will have (twice) the kinetic energy.  This way if we interpolate in   !
                  ! the vertical for the multi-site case, we interpolate energy.           !
                  !------------------------------------------------------------------------!
                  do ipy = 1,cgrid%npolygons
                     cgrid%met(ipy)%vels = cgrid%met(ipy)%vels                             &
                                         + cgrid%metinput(ipy)%vgrd(mprev)**2
                  end do

               case('ustar')   !----- Prescribed u*. [    m/s] ----------------------------!
                  !------------------------------------------------------------------------!
                  !     This is used only when the user wants u* to be prescribed instead  !
                  ! of calculated.  We add the square of the u* so we integrate and        !
                  ! interpolate momentum flux.                                             !
                  !------------------------------------------------------------------------!
                  do ipy = 1,cgrid%npolygons
                     cgrid%met(ipy)%atm_ustar = cgrid%metinput(ipy)%atm_ustar(mprev)**2
                  end do

               case('sh')      !----- Specific humidity. [kg/kg_a] ------------------------!
                  !------------------------------------------------------------------------!
                  !     Here we will just get the specific humidity. But soon we will      !
                  ! check whether this value is consistent with our thermodynamics (i.e.,  !
                  ! it will not be supersaturated).  If not we will impose the maximum.    !
                  !------------------------------------------------------------------------!
                  do ipy = 1,cgrid%npolygons
                     cgrid%met(ipy)%atm_shv = cgrid%metinput(ipy)%sh(mprev)
                  end do

               case('tmp')    !------ Air temperature. [      K] --------------------------!
                  do ipy = 1,cgrid%npolygons
                     cgrid%met(ipy)%atm_tmp = cgrid%metinput(ipy)%tmp (mprev)
                  end do

               case('co2')     !----- CO2 mixing ratio. [umol/mol] ------------------------!
                  do ipy = 1,cgrid%npolygons
                     cgrid%met(ipy)%atm_co2 = cgrid%metinput(ipy)%co2(mprev)
                  end do
               end select
               
            !----- In this case, we need to interpolate. ----------------------------------!
            case (1,5)
               !---------------------------------------------------------------------------!
               !      If we fall in this part of the if block, this means that the vari-   !
               ! able must be interpolated over time.                                      !
               !---------------------------------------------------------------------------!

               !----- First, get the number of points per day and per month. --------------!
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
               !---------------------------------------------------------------------------!



               !----- Find the number of seconds since the beginning of this month. -------!
               ndays_elapsed    = current_time%date - 1
               nseconds_elapsed = nint(current_time%time) + ndays_elapsed * nint(day_sec)
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               ! mprev  - Index of the element of the previous metinput data.              !
               ! mnext  - Index of the element of the next metinput data.                  !
               ! dtprev = Time in seconds relative to the previous input time.             !
               ! dtnext = Time in seconds relative to the next input time.                 !
               !---------------------------------------------------------------------------!
               mprev  = int(real(nseconds_elapsed)/met_frq(iformat,iv)) + 1
               dtprev = floor(real(nseconds_elapsed)/met_frq(iformat,iv))                  &
                      * met_frq(iformat,iv) - real(nseconds_elapsed)
               mnext  = mprev + 1
               dtnext = dtprev + met_frq(iformat,iv)
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !     Find the actual time of the previous and the upcoming met             !
               ! information.  We need this to find the integral of the secant of the      !
               ! zenith angle, which is used to make a time interpolation that preserves   !
               ! the diurnal cycle consistent.  This is particularly important for         !
               ! tropical regions because the zenith angle changes very quickly and the    !
               ! radiation could be too skewed even with hourly data.                      !
               !---------------------------------------------------------------------------!
               prevmet_timea = current_time
               nextmet_timea = current_time
               select case (imetavg)
               case (0)
                  !----- Instantaneous averages. ------------------------------------------!
                  call update_model_time_dm(prevmet_timea,dtprev)
                  call update_model_time_dm(nextmet_timea,dtnext)
                  !------------------------------------------------------------------------!

               case (1)
                  !----- Averages ending at the reference time. ---------------------------!
                  call update_model_time_dm(prevmet_timea,dtprev - met_frq(iformat,iv))
                  call update_model_time_dm(nextmet_timea,dtnext - met_frq(iformat,iv))
                  !------------------------------------------------------------------------!

               case (2)
                  !----- Averages beginning at the reference time. ------------------------!
                  call update_model_time_dm(prevmet_timea,dtprev)
                  call update_model_time_dm(nextmet_timea,dtnext)
                  !------------------------------------------------------------------------!

               case (3)
                  !----- Averages centred at the reference time. --------------------------!
                  call update_model_time_dm( prevmet_timea                                 &
                                           , dtprev - 0.5 * met_frq(iformat,iv) )
                  call update_model_time_dm( nextmet_timea                                 &
                                           , dtnext - 0.5 * met_frq(iformat,iv) )
                  !------------------------------------------------------------------------!

               end select
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !     For the particular case where we at the last day of the month, after  !
               ! the last met driver data of that month, the next time is the first day of !
               ! the next month.  This will always be in the next position after day 31,   !
               ! even when the current month is shorter.                                   !
               !---------------------------------------------------------------------------!
               if(mnext > np)then
                  mnext = 1 + nint(day_sec / met_frq(iformat,iv)) * 31
               endif
               
               !---------------------------------------------------------------------------!
               !      Get interpolation factors.                                           !
               !---------------------------------------------------------------------------!
               wnext = max(0.0,min(1.0,(-dtprev)/met_frq(iformat,iv)))
               wprev = 1.0 - wnext

               !----- Find the variable and fill the sites. -------------------------------!
               select case (trim(met_vars(iformat,iv)))
               case('prate')   !----- Precipitation rate. [kg/m2/s] -----------------------!
                  do ipy = 1,cgrid%npolygons
                     cgrid%met(ipy)%pcpg = cgrid%metinput(ipy)%prate(mnext) * wnext        &
                                         + cgrid%metinput(ipy)%prate(mprev) * wprev
                  end do

               case('dlwrf')   !----- Downwelling longwave radiation. [   W/m2] -----------!
                  do ipy = 1,cgrid%npolygons
                     cgrid%met(ipy)%rlong = cgrid%metinput(ipy)%dlwrf(mnext) * wnext       &
                                          + cgrid%metinput(ipy)%dlwrf(mprev) * wprev
                  end do

               case('pres')
                  !------------------------------------------------------------------------!
                  !     Air pressure [   N/m2] and Exner function [J/kg/K].                !
                  !------------------------------------------------------------------------!
                  do ipy = 1,cgrid%npolygons
                     cgrid%met(ipy)%prss  = cgrid%metinput(ipy)%pres(mnext) * wnext        &
                                          + cgrid%metinput(ipy)%pres(mprev) * wprev
                  end do

               case('hgt')     !----- Reference height. [      m] -------------------------!
                  do ipy = 1,cgrid%npolygons
                     cgrid%met(ipy)%geoht = cgrid%metinput(ipy)%hgt(mnext) * wnext         &
                                          + cgrid%metinput(ipy)%hgt(mprev) * wprev
                     !---------------------------------------------------------------------!
                     !    The reference height must be always above the top cohort.  We    !
                     ! never know when this may happen, but if there is a risk of the      !
                     ! reference height becoming lower than the top cohort, stop the run   !
                     ! and explain it to the user.                                         !
                     !---------------------------------------------------------------------!
                     do ipft = 1,n_pft
                        if (include_pft(ipft) .and. cgrid%met(ipy)%geoht <= hgt_max(ipft)) &
                        then
                           write (unit=*,fmt='(a)') '-------------------------------------'
                           write(unit=*,fmt='(a,1x,i12)')                                  &
                                                 ' - Polygon        =',ipy
                           write(unit=*,fmt='(a,1x,i12)')                                  &
                                                 ' - PFT            =',ipft
                           write (unit=*,fmt='(a,1x,es12.5)')                              &
                                                 ' - Ref. Height    =',cgrid%met(ipy)%geoht
                           write (unit=*,fmt='(a,1x,es12.5)')                              &
                                                 ' - PFT max height =',hgt_max(ipft)
                           write (unit=*,fmt='(a)') '-------------------------------------'
                           call fatal_error('Reference height must be higher than PFT '//  &
                                            'maximum height!!','update_met_drivers'        &
                                           ,'ed_met_driver.f90')
                        end if
                     end do
                  end do
                  !------------------------------------------------------------------------!

               case('ugrd')    !----- Zonal wind. [    m/s] -------------------------------!
                  !------------------------------------------------------------------------!
                  !     Here we add the square of the zonal wind, so in the end vels will  !
                  ! have (twice) the kinetic energy.  This way if we interpolate in the    !
                  ! vertical for the multi-site case, we interpolate energy.               !
                  !------------------------------------------------------------------------!
                  do ipy = 1,cgrid%npolygons
                     cgrid%met(ipy)%vels =  cgrid%met(ipy)%vels                            &
                                         +  cgrid%metinput(ipy)%ugrd(mnext) ** 2 * wnext   &
                                         +  cgrid%metinput(ipy)%ugrd(mprev) ** 2 * wprev
                  end do

               case('vgrd')    !----- Meridional wind. [    m/s] --------------------------!
                  !------------------------------------------------------------------------!
                  !     Here we add the square of the meridional wind, so in the end vels  !
                  ! will have (twice) the kinetic energy.  This way if we interpolate in   !
                  ! the vertical for the multi-site case, we interpolate energy.           !
                  !------------------------------------------------------------------------!
                  do ipy = 1,cgrid%npolygons
                     cgrid%met(ipy)%vels = cgrid%met(ipy)%vels                             &
                                         + cgrid%metinput(ipy)%vgrd(mnext) **2 * wnext     &
                                         + cgrid%metinput(ipy)%vgrd(mprev) **2 * wprev
                  enddo

               case('ustar')   !----- Friction velocity. [    m/s] ------------------------!
                  !------------------------------------------------------------------------!
                  !     Here we add the square of the friction velocity, so in the end     !
                  ! atm_ustar will have (twice) the turbulent kinetic energy.  This way if !
                  ! we interpolate in the vertical for the multi-site case, we interpolate !
                  ! energy.                                                                !
                  !------------------------------------------------------------------------!
                  do ipy = 1,cgrid%npolygons
                     cgrid%met(ipy)%atm_ustar =                                            &
                                           cgrid%metinput(ipy)%atm_ustar(mnext)**2 * wnext &
                                         + cgrid%metinput(ipy)%atm_ustar(mprev)**2 * wprev
                  enddo
               case('sh')      !----- Specific humidity. [kg/kg_a] ------------------------!
                  do ipy = 1,cgrid%npolygons
                     cgrid%met(ipy)%atm_shv = cgrid%metinput(ipy)%sh(mnext) * wnext        &
                                            + cgrid%metinput(ipy)%sh(mprev) * wprev
                  end do

               case('tmp')     !----- Air temperature. [      K]---------------------------!
                  do ipy = 1,cgrid%npolygons
                     cgrid%met(ipy)%atm_tmp = cgrid%metinput(ipy)%tmp(mnext) * wnext       &
                                            + cgrid%metinput(ipy)%tmp(mprev) * wprev
                     !---------------------------------------------------------------------!
                  end do

               case('nbdsf')   !----- Near IR beam downward shortwave flux. [   W/m2] -----!

                  select case (imetavg)
                  case (-1)
                     do ipy=1,cgrid%npolygons
                        cgrid%met(ipy)%nir_beam = cgrid%metinput(ipy)%nbdsf(mnext) * wnext &
                                                + cgrid%metinput(ipy)%nbdsf(mprev) * wprev
                     end do
                  case default


                     do ipy= 1,cgrid%npolygons
                        !------------------------------------------------------------------!
                        !     Check whether this is day time or night time.  We should     !
                        ! only do the full interpolation thing during the day, the night   !
                        ! should have no incoming radiation until we incorporate fireflies !
                        ! (or a good twilight scheme) in the model.                        !
                        !------------------------------------------------------------------!
                        if (cgrid%cosz(ipy) > cosz_min) then

                           !---------------------------------------------------------------!
                           !     Define the normalisation factors for the previous and the !
                           ! next time.                                                    !
                           !---------------------------------------------------------------!
                           select case (imetavg)
                           case (0)
                              secz_prev = mean_daysecz(cgrid%lon(ipy),cgrid%lat(ipy)       &
                                                      ,prevmet_timea,dt_radinterp,0.)
                              secz_next = mean_daysecz(cgrid%lon(ipy),cgrid%lat(ipy)       &
                                                      ,nextmet_timea,dt_radinterp,0.)
                           case default
                              secz_prev = mean_daysecz(cgrid%lon(ipy),cgrid%lat(ipy)       &
                                                      ,prevmet_timea,dt_radinterp          &
                                                      ,met_frq(iformat,iv))
                              secz_next = mean_daysecz(cgrid%lon(ipy),cgrid%lat(ipy)       &
                                                      ,nextmet_timea,dt_radinterp          &
                                                      ,met_frq(iformat,iv))
                           end select
                           !---------------------------------------------------------------!



                           !---------------------------------------------------------------!
                           !     Decide whether we want to use the previous or the next    !
                           ! data based on the zenith angle.                               !
                           !---------------------------------------------------------------!
                           night_prev = secz_prev == 0.
                           night_next = secz_next == 0.
                           !---------------------------------------------------------------!




                           !---------------------------------------------------------------!
                           !     Decide what to do based on the weighting factors.  This   !
                           ! is to avoid putting too little or too much energy at noon,    !
                           ! dawn, and dusk.                                               !
                           !---------------------------------------------------------------!
                           if (night_prev .and. night_next) then
                              !----- Middle of the night, set the value to zero. ----------!
                              fperp_prev = 0.
                              fperp_next = 0.
                              cgrid%met(ipy)%nir_beam = 0.
                           elseif (night_prev) then
                              !------------------------------------------------------------!
                              !     Dawn time, the previous is zero so we only have        !
                              ! meaningful information for the next time: we scale the     !
                              ! next time using the next secant and scaling back with the  !
                              ! current cosine of zenith angle.                            !
                              !------------------------------------------------------------!
                              fperp_prev = 0.
                              fperp_next = cgrid%metinput(ipy)%nbdsf(mnext) * secz_next
                              cgrid%met(ipy)%nir_beam = fperp_next * cgrid%cosz(ipy)
                           elseif (night_next) then
                              !-----------------------------------------------------------!
                              !     Dusk time, the next time is zero so we only have      !
                              ! meaningful information for the previous time: we scale    !
                              ! the next time using the previous secant and scaling back  !
                              ! with the current cosine of zenith angle.                  !
                              !-----------------------------------------------------------!
                              fperp_prev = cgrid%metinput(ipy)%nbdsf(mprev) * secz_prev
                              fperp_next = 0.
                              cgrid%met(ipy)%nir_beam = fperp_prev * cgrid%cosz(ipy)
                           else
                              !----- Daytime, use both previous and next values. ---------!
                              fperp_next = cgrid%metinput(ipy)%nbdsf(mnext) * secz_next
                              fperp_prev = cgrid%metinput(ipy)%nbdsf(mprev) * secz_prev
                              cgrid%met(ipy)%nir_beam = cgrid%cosz(ipy)                    &
                                                      * ( fperp_next * wnext               &
                                                        + fperp_prev * wprev )
                           end if
                           !---------------------------------------------------------------!
                        else
                           !----- Night time, assign it zero. -----------------------------!
                           secz_prev               = 0.
                           secz_next               = 0.
                           night_prev              = .true.
                           night_next              = .true.
                           fperp_prev              = 0.
                           fperp_next              = 0.
                           cgrid%met(ipy)%nir_beam = 0.0
                           
                           !---------------------------------------------------------------!
                        end if
                        !------------------------------------------------------------------!



                        !------------------------------------------------------------------!
                        !     Print the debugging information if the user wants so.        !
                        !------------------------------------------------------------------!
                        if (print_radinterp) then
                           call dump_radinfo(cgrid,ipy,mprev,mnext,prevmet_timea           &
                                            ,nextmet_timea,met_frq(iformat,iv),wprev,wnext &
                                            ,secz_prev,secz_next,fperp_prev,fperp_next     &
                                            ,night_prev,night_next,'nbdsf')
                        end if
                        !------------------------------------------------------------------!
                     end do
                     !---------------------------------------------------------------------!
                  end select
                  !------------------------------------------------------------------------!


               case('nddsf')   !----- Near IR diffuse downward shortwave flux. [   W/m2] --!

                  select case (imetavg)
                  case (-1)
                     do ipy=1,cgrid%npolygons
                        cgrid%met(ipy)%nir_diffuse =                                       &
                                         cgrid%metinput(ipy)%nddsf(mnext) * wnext          &
                                       + cgrid%metinput(ipy)%nddsf(mprev) * wprev
                     end do
                  case default


                     do ipy= 1,cgrid%npolygons
                        !------------------------------------------------------------------!
                        !     Check whether this is day time or night time.  We should     !
                        ! only do the full interpolation thing during the day, the night   !
                        ! should have no incoming radiation until we incorporate fireflies !
                        ! (or a good twilight scheme) in the model.                        !
                        !------------------------------------------------------------------!
                        if (cgrid%cosz(ipy) > cosz_min) then

                           !---------------------------------------------------------------!
                           !     Define the normalisation factors for the previous and the !
                           ! next time.                                                    !
                           !---------------------------------------------------------------!
                           select case (imetavg)
                           case (0)
                              secz_prev = mean_daysecz(cgrid%lon(ipy),cgrid%lat(ipy)       &
                                                      ,prevmet_timea,dt_radinterp,0.)
                              secz_next = mean_daysecz(cgrid%lon(ipy),cgrid%lat(ipy)       &
                                                      ,nextmet_timea,dt_radinterp,0.)
                           case default
                              secz_prev = mean_daysecz(cgrid%lon(ipy),cgrid%lat(ipy)       &
                                                      ,prevmet_timea,dt_radinterp          &
                                                      ,met_frq(iformat,iv))
                              secz_next = mean_daysecz(cgrid%lon(ipy),cgrid%lat(ipy)       &
                                                      ,nextmet_timea,dt_radinterp          &
                                                      ,met_frq(iformat,iv))
                           end select
                           !---------------------------------------------------------------!



                           !---------------------------------------------------------------!
                           !     Decide whether we want to use the previous or the next    !
                           ! data based on the zenith angle.                               !
                           !---------------------------------------------------------------!
                           night_prev = secz_prev == 0.
                           night_next = secz_next == 0.
                           !---------------------------------------------------------------!




                           !---------------------------------------------------------------!
                           !     Decide what to do based on the weighting factors.  This   !
                           ! is to avoid putting too little or too much energy at noon,    !
                           ! dawn, and dusk.                                               !
                           !---------------------------------------------------------------!
                           if (night_prev .and. night_next) then
                              !----- Middle of the night, set the value to zero. ----------!
                              fperp_prev = 0.
                              fperp_next = 0.
                              cgrid%met(ipy)%nir_diffuse = 0.
                           elseif (night_prev) then
                              !------------------------------------------------------------!
                              !     Dawn time, the previous is zero so we only have        !
                              ! meaningful information for the next time: we scale the     !
                              ! next time using the next secant and scaling back with the  !
                              ! current cosine of zenith angle.                            !
                              !------------------------------------------------------------!
                              fperp_prev = 0.
                              fperp_next = cgrid%metinput(ipy)%nddsf(mnext) * secz_next
                              cgrid%met(ipy)%nir_diffuse = fperp_next * cgrid%cosz(ipy)
                           elseif (night_next) then
                              !------------------------------------------------------------!
                              !     Dusk time, the next time is zero so we only have       !
                              ! meaningful information for the previous time: we scale the !
                              ! next time using the previous secant and scaling back with  !
                              ! the current cosine of zenith angle.                        !
                              !------------------------------------------------------------!
                              fperp_prev = cgrid%metinput(ipy)%nddsf(mprev) * secz_prev
                              fperp_next = 0.
                              cgrid%met(ipy)%nir_diffuse = fperp_prev * cgrid%cosz(ipy)
                           else
                              !----- Daytime, use both previous and next values. ----------!
                              fperp_next = cgrid%metinput(ipy)%nddsf(mnext) * secz_next
                              fperp_prev = cgrid%metinput(ipy)%nddsf(mprev) * secz_prev
                              cgrid%met(ipy)%nir_diffuse = cgrid%cosz(ipy)                 &
                                                         * ( fperp_next * wnext            &
                                                           + fperp_prev * wprev )
                           end if
                           !---------------------------------------------------------------!
                        else
                           !----- Night time, assign it zero. -----------------------------!
                           secz_prev                  = 0.
                           secz_next                  = 0.
                           night_prev                 = .true.
                           night_next                 = .true.
                           fperp_prev                 = 0.
                           fperp_next                 = 0.
                           cgrid%met(ipy)%nir_diffuse = 0.0
                           
                           !---------------------------------------------------------------!
                        end if
                        !------------------------------------------------------------------!



                        !------------------------------------------------------------------!
                        !     Print the debugging information if the user wants so.        !
                        !------------------------------------------------------------------!
                        if (print_radinterp) then
                           call dump_radinfo(cgrid,ipy,mprev,mnext,prevmet_timea           &
                                            ,nextmet_timea,met_frq(iformat,iv),wprev,wnext &
                                            ,secz_prev,secz_next,fperp_prev,fperp_next     &
                                            ,night_prev,night_next,'nddsf')
                        end if
                        !------------------------------------------------------------------!
                     end do
                     !---------------------------------------------------------------------!
                  end select
                  !------------------------------------------------------------------------!

               case('vbdsf')   !----- Visible beam downward shortwave flux. [   W/m2] -----!

                  select case (imetavg)
                  case (-1)
                     do ipy=1,cgrid%npolygons
                        cgrid%met(ipy)%par_beam = cgrid%metinput(ipy)%vbdsf(mnext) * wnext &
                                                + cgrid%metinput(ipy)%vbdsf(mprev) * wprev
                     end do
                  case default


                     do ipy= 1,cgrid%npolygons
                        !------------------------------------------------------------------!
                        !     Check whether this is day time or night time.  We should     !
                        ! only do the full interpolation thing during the day, the night   !
                        ! should have no incoming radiation until we incorporate fireflies !
                        ! (or a good twilight scheme) in the model.                        !
                        !------------------------------------------------------------------!
                        if (cgrid%cosz(ipy) > cosz_min) then

                           !---------------------------------------------------------------!
                           !     Define the normalisation factors for the previous and the !
                           ! next time.                                                    !
                           !---------------------------------------------------------------!
                           select case (imetavg)
                           case (0)
                              secz_prev = mean_daysecz(cgrid%lon(ipy),cgrid%lat(ipy)       &
                                                      ,prevmet_timea,dt_radinterp,0.)
                              secz_next = mean_daysecz(cgrid%lon(ipy),cgrid%lat(ipy)       &
                                                      ,nextmet_timea,dt_radinterp,0.)
                           case default
                              secz_prev = mean_daysecz(cgrid%lon(ipy),cgrid%lat(ipy)       &
                                                      ,prevmet_timea,dt_radinterp          &
                                                      ,met_frq(iformat,iv))
                              secz_next = mean_daysecz(cgrid%lon(ipy),cgrid%lat(ipy)       &
                                                      ,nextmet_timea,dt_radinterp          &
                                                      ,met_frq(iformat,iv))
                           end select
                           !---------------------------------------------------------------!



                           !---------------------------------------------------------------!
                           !     Decide whether we want to use the previous or the next    !
                           ! data based on the zenith angle.                               !
                           !---------------------------------------------------------------!
                           night_prev = secz_prev == 0.
                           night_next = secz_next == 0.
                           !---------------------------------------------------------------!




                           !---------------------------------------------------------------!
                           !     Decide what to do based on the weighting factors.  This   !
                           ! is to avoid putting too little or too much energy at noon,    !
                           ! dawn, and dusk.                                               !
                           !---------------------------------------------------------------!
                           if (night_prev .and. night_next) then
                              !----- Middle of the night, set the value to zero. ----------!
                              fperp_prev = 0.
                              fperp_next = 0.
                              cgrid%met(ipy)%par_beam = 0.
                           elseif (night_prev) then
                              !------------------------------------------------------------!
                              !     Dawn time, the previous is zero so we only have        !
                              ! meaningful information for the next time: we scale the     !
                              ! next time using the next secant and scaling back with the  !
                              ! current cosine of zenith angle.                            !
                              !------------------------------------------------------------!
                              fperp_prev = 0.
                              fperp_next = cgrid%metinput(ipy)%vbdsf(mnext) * secz_next
                              cgrid%met(ipy)%par_beam = fperp_next * cgrid%cosz(ipy)
                           elseif (night_next) then
                              !------------------------------------------------------------!
                              !     Dusk time, the next time is zero so we only have       !
                              ! meaningful information for the previous time: we scale the !
                              ! next time using the previous secant and scaling back with  !
                              ! the current cosine of zenith angle.                        !
                              !------------------------------------------------------------!
                              fperp_prev = cgrid%metinput(ipy)%vbdsf(mprev) * secz_prev
                              fperp_next = 0.
                              cgrid%met(ipy)%par_beam = fperp_prev * cgrid%cosz(ipy)
                           else
                              !----- Daytime, use both previous and next values. ----------!
                              fperp_next = cgrid%metinput(ipy)%vbdsf(mnext) * secz_next
                              fperp_prev = cgrid%metinput(ipy)%vbdsf(mprev) * secz_prev
                             cgrid%met(ipy)%par_beam = cgrid%cosz(ipy)                     &
                                                     * ( fperp_next * wnext                &
                                                        + fperp_prev * wprev )
                           end if
                           !---------------------------------------------------------------!
                        else
                           !----- Night time, assign it zero. -----------------------------!
                           secz_prev               = 0.
                           secz_next               = 0.
                           night_prev              = .true.
                           night_next              = .true.
                           fperp_prev              = 0.
                           fperp_next              = 0.
                           cgrid%met(ipy)%par_beam = 0.0
                           
                           !---------------------------------------------------------------!
                        end if
                        !------------------------------------------------------------------!



                        !------------------------------------------------------------------!
                        !     Print the debugging information if the user wants so.        !
                        !------------------------------------------------------------------!
                        if (print_radinterp) then
                           call dump_radinfo(cgrid,ipy,mprev,mnext,prevmet_timea           &
                                            ,nextmet_timea,met_frq(iformat,iv),wprev,wnext &
                                            ,secz_prev,secz_next,fperp_prev,fperp_next     &
                                            ,night_prev,night_next,'vbdsf')
                        end if
                        !------------------------------------------------------------------!
                     end do
                     !---------------------------------------------------------------------!
                  end select
                  !------------------------------------------------------------------------!


               case('vddsf')   !----- Visible diffuse downward shortwave flux. [   W/m2] --!

                  select case (imetavg)
                  case (-1)
                     do ipy=1,cgrid%npolygons
                        cgrid%met(ipy)%par_diffuse =                                       &
                                           cgrid%metinput(ipy)%vddsf(mnext) * wnext        &
                                         + cgrid%metinput(ipy)%vddsf(mprev) * wprev
                     end do
                  case default

                     do ipy= 1,cgrid%npolygons
                        !------------------------------------------------------------------!
                        !     Check whether this is day time or night time.  We should     !
                        ! only do the full interpolation thing during the day, the night   !
                        ! should have no incoming radiation until we incorporate fireflies !
                        ! (or a good twilight scheme) in the model.                        !
                        !------------------------------------------------------------------!
                        if (cgrid%cosz(ipy) > cosz_min) then

                           !---------------------------------------------------------------!
                           !     Define the normalisation factors for the previous and the !
                           ! next time.                                                    !
                           !---------------------------------------------------------------!
                           select case (imetavg)
                           case (0)
                              secz_prev = mean_daysecz(cgrid%lon(ipy),cgrid%lat(ipy)       &
                                                      ,prevmet_timea,dt_radinterp,0.)
                              secz_next = mean_daysecz(cgrid%lon(ipy),cgrid%lat(ipy)       &
                                                      ,nextmet_timea,dt_radinterp,0.)
                           case default
                              secz_prev = mean_daysecz(cgrid%lon(ipy),cgrid%lat(ipy)       &
                                                      ,prevmet_timea,dt_radinterp          &
                                                      ,met_frq(iformat,iv))
                              secz_next = mean_daysecz(cgrid%lon(ipy),cgrid%lat(ipy)       &
                                                      ,nextmet_timea,dt_radinterp          &
                                                      ,met_frq(iformat,iv))
                           end select
                           !---------------------------------------------------------------!



                           !---------------------------------------------------------------!
                           !     Decide whether we want to use the previous or the next    !
                           ! data based on the zenith angle.                               !
                           !---------------------------------------------------------------!
                           night_prev = secz_prev == 0.
                           night_next = secz_next == 0.
                           !---------------------------------------------------------------!




                           !---------------------------------------------------------------!
                           !     Decide what to do based on the weighting factors.  This   !
                           ! is to avoid putting too little or too much energy at noon,    !
                           ! dawn, and dusk.                                               !
                           !---------------------------------------------------------------!
                           if (night_prev .and. night_next) then
                              !----- Middle of the night, set the value to zero. ----------!
                              fperp_prev = 0.
                              fperp_next = 0.
                              cgrid%met(ipy)%par_diffuse = 0.
                           elseif (night_prev) then
                              !------------------------------------------------------------!
                              !     Dawn time, the previous is zero so we only have        !
                              ! meaningful information for the next time: we scale the     !
                              ! next time using the next secant and scaling back with the  !
                              ! current cosine of zenith angle.                            !
                              !------------------------------------------------------------!
                              fperp_prev = 0.
                              fperp_next = cgrid%metinput(ipy)%vddsf(mnext) * secz_next
                              cgrid%met(ipy)%par_diffuse = fperp_next * cgrid%cosz(ipy)
                           elseif (night_next) then
                              !------------------------------------------------------------!
                              !     Dusk time, the next time is zero so we only have       !
                              ! meaningful information for the previous time: we scale the !
                              ! next time using the previous secant and scaling back with  !
                              ! the current cosine of zenith angle.                        !
                              !------------------------------------------------------------!
                              fperp_prev = cgrid%metinput(ipy)%vddsf(mprev) * secz_prev
                              fperp_next = 0.
                              cgrid%met(ipy)%par_diffuse = fperp_prev * cgrid%cosz(ipy)
                           else
                              !----- Daytime, use both previous and next values. ----------!
                              fperp_next = cgrid%metinput(ipy)%vddsf(mnext) * secz_next
                              fperp_prev = cgrid%metinput(ipy)%vddsf(mprev) * secz_prev
                              cgrid%met(ipy)%par_diffuse = cgrid%cosz(ipy)                 &
                                                         * ( fperp_next * wnext            &
                                                           + fperp_prev * wprev )
                           end if
                           !---------------------------------------------------------------!
                        else
                           !----- Night time, assign it zero. -----------------------------!
                           secz_prev                  = 0.
                           secz_next                  = 0.
                           night_prev                 = .true.
                           night_next                 = .true.
                           fperp_prev                 = 0.
                           fperp_next                 = 0.
                           cgrid%met(ipy)%par_diffuse = 0.0
                           
                           !---------------------------------------------------------------!
                        end if
                        !------------------------------------------------------------------!



                        !------------------------------------------------------------------!
                        !     Print the debugging information if the user wants so.        !
                        !------------------------------------------------------------------!
                        if (print_radinterp) then
                           call dump_radinfo(cgrid,ipy,mprev,mnext,prevmet_timea           &
                                            ,nextmet_timea,met_frq(iformat,iv),wprev,wnext &
                                            ,secz_prev,secz_next,fperp_prev,fperp_next     &
                                            ,night_prev,night_next,'vddsf')
                        end if
                        !------------------------------------------------------------------!
                     end do
                     !---------------------------------------------------------------------!
                  end select
                  !------------------------------------------------------------------------!

               case('co2')     !----- CO2 mixing ratio. [umol/mol] ------------------------!
                  do ipy = 1,cgrid%npolygons
                     cgrid%met(ipy)%atm_co2 = cgrid%metinput(ipy)%co2(mnext) * wnext       &
                                            + cgrid%metinput(ipy)%co2(mprev) * wprev
                  end do

               end select
            end select
         end do varloop
      end do formloop
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !     Now that all variables from the met driver were read and updated, we update    !
      ! the derived variables.                                                             !
      !------------------------------------------------------------------------------------!
      polyloop: do ipy = 1,cgrid%npolygons
         !----- CO2 (only if it hasn't been read). ----------------------------------------!
         if (.not. has_co2) cgrid%met(ipy)%atm_co2 = initial_co2
         !---------------------------------------------------------------------------------!

         !----- Set the default Exner function from pressure. -----------------------------!
         cgrid%met(ipy)%exner = press2exner(cgrid%met(ipy)%prss)
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !     Set default potential temperature from Exner function and air temperature.  !
         !---------------------------------------------------------------------------------!
         cgrid%met(ipy)%atm_theta = extemp2theta( cgrid%met(ipy)%exner                     &
                                                , cgrid%met(ipy)%atm_tmp )
         temp0                    = cgrid%met(ipy)%atm_tmp
         !---------------------------------------------------------------------------------!

         if (atm_tmp_intercept /= 0.0 .or. atm_tmp_slope /= 1.0) then
            cgrid%met(ipy)%atm_tmp = atm_tmp_intercept                                     &
                                   + cgrid%met(ipy)%atm_tmp * atm_tmp_slope
            !----- We must update potential temperature too. ------------------------------!
            cgrid%met(ipy)%atm_theta = extemp2theta( cgrid%met(ipy)%exner                  &
                                                   , cgrid%met(ipy)%atm_tmp )
            !------------------------------------------------------------------------------!
         end if
         cgrid%met(ipy)%pcpg = max(0.0,prec_intercept + cgrid%met(ipy)%pcpg * prec_slope)
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !      Modify the short-wave radiation splitting in case the user wants so.       !
         !---------------------------------------------------------------------------------!
         call solar_radiation_breakdown(cgrid,ipy)
         !---------------------------------------------------------------------------------!


         if (humid_scenario == 1) then 
            !------------------------------------------------------------------------------!
            !     Update atm_shv so the relative humidity remains the same.  We use the    !
            ! functions from therm_lib.f90 so it is consistent with the rest of the code.  !
            !------------------------------------------------------------------------------!
            !----- 1. Find relative humidity. ---------------------------------------------!
            relhum                 = rehuil( cgrid%met(ipy)%prss                           &
                                           , temp0                                         &
                                           , cgrid%met(ipy)%atm_shv                        &
                                           , .true.                                        )
            !----- 2. Find the equivalent relative humidity with the new temperature. -----!
            cgrid%met(ipy)%atm_shv = ptrh2rvapil( relhum                                   &
                                                , cgrid%met(ipy)%prss                      &
                                                , cgrid%met(ipy)%atm_tmp                   &
                                                , .true.                                   )
            !------------------------------------------------------------------------------!
         end if


         !---------------------------------------------------------------------------------!
         !     Check the relative humidity associated with the current pressure, temper-   !
         ! ature, and specific humidity.  Impose lower and upper bounds as prescribed by   !
         ! the variables atm_rhv_min and atm_rhv_max (from met_driver_coms.f90, and        !
         ! defined at the init_met_params subroutine in ed_params.f90).                    !
         !---------------------------------------------------------------------------------!
         relhum = rehuil(cgrid%met(ipy)%prss,cgrid%met(ipy)%atm_tmp,cgrid%met(ipy)%atm_shv &
                        ,.true.)
         !---------------------------------------------------------------------------------!
         !      Check whether the relative humidity is off-bounds.  If it is, then we re-  !
         ! calculate the mixing ratio and convert to specific humidity.                    !
         !---------------------------------------------------------------------------------!
         if (relhum < atm_rhv_min) then
            relhum                 = atm_rhv_min
            cgrid%met(ipy)%atm_shv = ptrh2rvapil( relhum                                   &
                                                , cgrid%met(ipy)%prss                      &
                                                , cgrid%met(ipy)%atm_tmp                   &
                                                , .true.                                   )
         elseif (relhum > atm_rhv_max) then
            relhum                 = atm_rhv_max
            cgrid%met(ipy)%atm_shv = ptrh2rvapil( relhum                                   &
                                                , cgrid%met(ipy)%prss                      &
                                                , cgrid%met(ipy)%atm_tmp                   &
                                                , .true.                                   )
         end if

         !---------------------------------------------------------------------------------!
         !    We now find the equivalent potential temperature.                            !
         !---------------------------------------------------------------------------------!
         rvaux                    = cgrid%met(ipy)%atm_shv / (1.0 - cgrid%met(ipy)%atm_shv)
         cgrid%met(ipy)%atm_theiv = thetaeiv(cgrid%met(ipy)%atm_theta,cgrid%met(ipy)%prss  &
                                            ,cgrid%met(ipy)%atm_tmp,rvaux,rvaux)
         !---------------------------------------------------------------------------------!

         !---------------------------------------------------------------------------------!
         !    And the vapour pressure deficit.                                             !
         !---------------------------------------------------------------------------------!
         cgrid%met(ipy)%atm_vpdef = vpdefil(cgrid%met(ipy)%prss,cgrid%met(ipy)%atm_tmp     &
                                           ,cgrid%met(ipy)%atm_shv,.true.)
         !---------------------------------------------------------------------------------!


         !------ Apply met to sites, and adjust met variables for topography. -------------!
         call calc_met_lapse(cgrid,ipy)

         !----- Vels.  At this point vels is 2*Kinetic Energy, take the square root. ------!
         cgrid%met(ipy)%vels = sqrt(max(0.0,cgrid%met(ipy)%vels))


         !---------------------------------------------------------------------------------!
         !     Normally u* shouldn't be part of the met driver, so we check whether this   !
         ! is a test run.                                                                  !
         !---------------------------------------------------------------------------------!
         if (has_ustar) then
            !----- u*.  At this point u* is u*^2, take the square root. -------------------!
            cgrid%met(ipy)%atm_ustar = sqrt(max(ustmin*ustmin,cgrid%met(ipy)%atm_ustar))
            !------------------------------------------------------------------------------!
         else
            !----- No u* from met driver, assign a dummy value. ---------------------------!
            cgrid%met(ipy)%atm_ustar = ustmin
            !------------------------------------------------------------------------------!
         end if
         !---------------------------------------------------------------------------------!

         cpoly => cgrid%polygon(ipy)
         siteloop: do isi = 1,cpoly%nsites

            !----- Vels. The site level is also still in kinetic energy form. -------------!
            cpoly%met(isi)%vels        = sqrt(max(ubmin*ubmin,cpoly%met(isi)%vels))
            cpoly%met(isi)%vels_stab   = cpoly%met(isi)%vels
            cpoly%met(isi)%vels_unstab = cpoly%met(isi)%vels
            !------------------------------------------------------------------------------!



            !----- CO2.  In case we used the namelist, use that value. --------------------!
            if (.not. has_co2) cpoly%met(isi)%atm_co2 = initial_co2
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !     Normally u* shouldn't be part of the met driver, so we check whether     !
            ! this is a test run.                                                          !
            !------------------------------------------------------------------------------!
            if (has_ustar) then
               !----- u*.  At this point u* is u*^2, take the square root. ----------------!
               cpoly%met(isi)%atm_ustar = sqrt(max(ustmin*ustmin,cpoly%met(isi)%atm_ustar))
               !---------------------------------------------------------------------------!
            else
               !----- No u* from met driver, assign a dummy value. ------------------------!
               cpoly%met(isi)%atm_ustar = ustmin
               !---------------------------------------------------------------------------!
            end if
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !     We now find some derived properties.  In case several sites exist, the   !
            ! lapse rate was applied to pressure, temperature, and mixing ratio.  Then we  !
            ! calculate the Exner function, potential temperature and equivalent potential !
            ! temperature, so it will respect the ideal gas law and first law of thermo-   !
            ! dynamics.                                                                    !
            !------------------------------------------------------------------------------!
            cpoly%met(isi)%exner     = press2exner(cpoly%met(isi)%prss)
            cpoly%met(isi)%atm_theta = extemp2theta( cpoly%met(isi)%exner                  &
                                                   , cpoly%met(isi)%atm_tmp )
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !     Check the relative humidity associated with the current pressure,        !
            ! temperature, and specific humidity.  Impose lower and upper bounds as        !
            ! prescribed by the variables atm_rhv_min and atm_rhv_max (from                !
            ! met_driver_coms.f90, and defined at the init_met_params subroutine in        !
            ! ed_params.f90).                                                              !
            !------------------------------------------------------------------------------!
            relhum = rehuil(cpoly%met(isi)%prss,cpoly%met(isi)%atm_tmp                     &
                           ,cpoly%met(isi)%atm_shv,.true.)
            !------------------------------------------------------------------------------!
            !      Check whether the relative humidity is off-bounds.  If it is, then we   !
            ! re-calculate the mixing ratio and convert to specific humidity.              !
            !------------------------------------------------------------------------------!
            if (relhum < atm_rhv_min) then
               relhum = atm_rhv_min
               cpoly%met(isi)%atm_shv = ptrh2rvapil( relhum                                &
                                                   , cgrid%met(ipy)%prss                   &
                                                   , cgrid%met(ipy)%atm_tmp                &
                                                   , .true.                                )
            elseif (relhum > atm_rhv_max) then
               relhum                 = atm_rhv_max
               cpoly%met(isi)%atm_shv = ptrh2rvapil( relhum                                &
                                                   , cgrid%met(ipy)%prss                   &
                                                   , cgrid%met(ipy)%atm_tmp                &
                                                   , .true.                                )
            end if

            !----- Find the atmospheric equivalent potential temperature. -----------------!
            rvaux                    = cgrid%met(ipy)%atm_shv                              &
                                     / (1.0 - cgrid%met(ipy)%atm_shv)
            cpoly%met(isi)%atm_theiv = thetaeiv( cpoly%met(isi)%atm_theta                  &
                                               , cpoly%met(isi)%prss                       &
                                               , cpoly%met(isi)%atm_tmp,rvaux,rvaux)
            !------------------------------------------------------------------------------!

            !----- Find the vapour pressure deficit. --------------------------------------!
            cpoly%met(isi)%atm_vpdef = vpdefil (cpoly%met(isi)%prss,cpoly%met(isi)%atm_tmp &
                                               ,cpoly%met(isi)%atm_shv,.true.)
            !------------------------------------------------------------------------------!

            !----- Solar radiation --------------------------------------------------------!
            cpoly%met(isi)%rshort_diffuse = cpoly%met(isi)%par_diffuse                     &
                                          + cpoly%met(isi)%nir_diffuse
            cpoly%met(isi)%rshort         = cpoly%met(isi)%rshort_diffuse                  &
                                          + cpoly%met(isi)%par_beam                        &
                                          + cpoly%met(isi)%nir_beam
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !  Precipitation internal energy and "depth".                                  !
            !  fcns derived from                                                           !
            !  Jin et al 1999 Hydrol Process. 13:2467-2482 Table 2                         !
            !  [[modified 11/16/09 by MCD]]                                                !
            !------------------------------------------------------------------------------!
            if (cpoly%met(isi)%atm_tmp > (t3ple + 2.5)) then
               !----- Rain only. ----------------------------------------------------------!
               fliq = 1.0
               fice = 0.0
               cpoly%met(isi)%dpcpg = max(0.0, cpoly%met(isi)%pcpg) *wdnsi

            elseif (cpoly%met(isi)%atm_tmp <= (t3ple + 2.5) .and.                          &
                    cpoly%met(isi)%atm_tmp  > (t3ple + 2.0) ) then
               !---------------------------------------------------------------------------!
               !     Slight modification so the liquid fraction becomes continuous.        !
               !---------------------------------------------------------------------------!
               fliq  = min(1.0, 0.4 + 1.2 * (cpoly%met(isi)%atm_tmp - t3ple - 2.0))
               fice  = 1.0 - fliq
               snden = 189.0
               cpoly%met(isi)%dpcpg = max(0.0, cpoly%met(isi)%pcpg)                        &
                                    * ((1.0-fice) * wdnsi + fice / snden)

            elseif (cpoly%met(isi)%atm_tmp <= (t3ple + 2.0) .and.                          &
                    cpoly%met(isi)%atm_tmp > t3ple          ) then
               !---------------------------------------------------------------------------!
               !     Increasing the fraction of snow. (N.B. May not be appropriate for     !
               ! sub-tropical regions where the level of the melting layer is higher...).  !
               !---------------------------------------------------------------------------!
               fliq  = max(0.0, 0.2 * (cpoly%met(isi)%atm_tmp - t3ple))
               fice  = 1.0 - fliq
               snden = (50.0+1.7*(cpoly%met(isi)%atm_tmp-258.15)**1.5 )
               cpoly%met(isi)%dpcpg = max(0.0, cpoly%met(isi)%pcpg)                        &
                                    * ((1.0-fice) * wdnsi + fice / snden)

            elseif (cpoly%met(isi)%atm_tmp <= t3ple         .and.                          &
                    cpoly%met(isi)%atm_tmp > (t3ple - 15.0) ) then
               !----- Below freezing point, snow only. ------------------------------------!
               fliq  = 0.0
               fice  = 1.0
               snden = (50.0+1.7*(cpoly%met(isi)%atm_tmp-258.15)**1.5 )
               cpoly%met(isi)%dpcpg = max(0.0, cpoly%met(isi)%pcpg) / snden

            else ! if (copy%met(isi)%atm_tmp < (t3ple - 15.0)) then
               !----- Below freezing point, snow only. ------------------------------------!
               fliq  = 0.0
               fice  = 1.0
               snden = 50.
               cpoly%met(isi)%dpcpg = max(0.0, cpoly%met(isi)%pcpg) / snden

            end if
            !------------------------------------------------------------------------------!

            !------------------------------------------------------------------------------!
            !     Set internal energy.  This will be the precipitation times the specific  !
            ! internal energy of water (above or at triple point) multiplied by the liquid !
            ! fraction plus the specific internal energy of ice (below or at the triple    !
            ! point) multiplied by the ice fraction.                                       !
            !------------------------------------------------------------------------------!
            cpoly%met(isi)%qpcpg = max(0.0, cpoly%met(isi)%pcpg)                           &
                                 * ( fliq * tl2uint(max(t3ple,cpoly%met(isi)%atm_tmp),1.0) &
                                   + fice * tl2uint(min(t3ple,cpoly%met(isi)%atm_tmp),0.0) )
            !------------------------------------------------------------------------------!
         end do siteloop
         !---------------------------------------------------------------------------------!
      end do polyloop
      !------------------------------------------------------------------------------------!

      return
   end subroutine update_met_drivers
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !      This subroutine reads in the meteorological data.  Subroutine shdf5_irec_f uses  !
   ! HDF5 routines to extract the data from the archive.  The rest of this code manages    !
   ! the allocation and storage of that data.  The extracted data may be from a year,      !
   ! other than the current simulation year.  This happens when we use a range of years to !
   ! loop data over, typically during very long simulations.  In such cases, the number of !
   ! days in February will change during leap years, and the number of days in February    !
   ! during the model year may not match that of the simulation year, therefore this must  !
   ! be taken into account.                                                                !
   !=======================================================================================!
   subroutine read_ol_file(infile,iformat, iv, mname, year, offset, cgrid)
      use ed_max_dims    , only : str_len       ! ! intent(in)
      use ed_state_vars  , only : edtype        & ! structure
                                , polygontype   ! ! structure
      use met_driver_coms, only : met_frq       & ! intent(in)
                                , met_nlon      & ! intent(in)
                                , met_nlat      & ! intent(in)
                                , met_vars      & ! intent(in)
                                , met_interp    ! ! intent(in)
      use hdf5_utils     , only : shdf5_irec_f  & ! subroutine
                                , shdf5_info_f  ! ! subroutine
      use consts_coms    , only : day_sec       ! ! intent(in)
     
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      character(len=str_len), intent(in)  :: infile
      type(edtype)          , target      :: cgrid
      integer               , intent(in)  :: iformat
      integer               , intent(in)  :: iv
      character(len=3)      , intent(in)  :: mname
      integer               , intent(in)  :: year
      integer               , intent(in)  :: offset
      !----- Local variables. -------------------------------------------------------------!
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
      !----- External functions. ----------------------------------------------------------!
      logical, external                   :: isleap
      !------------------------------------------------------------------------------------!


      !----- Geoposition variables do not apply. ------------------------------------------!
      if (met_vars(iformat,iv).eq.'lat' .or. met_vars(iformat,iv).eq.'lon') return
     

      !------------------------------------------------------------------------------------!
      !  Determine how many times to read.                                                 !
      !------------------------------------------------------------------------------------!
     
      select case (met_interp(iformat,iv))
      case (2,4)
         !---------------------------------------------------------------------------------!
         !     In cases 2 and 4 we are dealing with constant time, only 1 point.           !
         !---------------------------------------------------------------------------------!
         np = 1
      case default

         !---------------------------------------------------------------------------------!
         !      Allocate for the number of points in the current year. In case this month  !
         ! is February, we must make sure that we account for the possibility of having a  !
         ! mismatch in the number of days in the model and in the dataset.                 !
         ! 1.  If the model year is leap but the dataset is not, we repeat February 28     !
         !     twice.                                                                      !
         ! 2.  If the model year is not leap but the dataset is, we drop February 29.      !
         !---------------------------------------------------------------------------------!

         points_per_day = nint(day_sec/met_frq(iformat,iv))

         !----- Determine the number of points necessary for the model year. --------------!
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

            !----- Find out the number of days in the input data set. ---------------------!
            call shdf5_info_f(met_vars(iformat,iv),ndims,idims)
            nday_dset=idims(1)/points_per_day
            !------------------------------------------------------------------------------!
         end select 

         np      = nday      * points_per_day
         np_dset = nday_dset * points_per_day

      end select

      !----- Allocate the buffer space. ---------------------------------------------------!
     
      select case (met_interp(iformat,iv))
      case (3,5) 
         !---------------------------------------------------------------------------------!
         !     Cases 3 and 5, reading in one value representing the entire domain, though  !
         ! they change over time.                                                          !
         !---------------------------------------------------------------------------------!
         ndims = 1
         idims(1) = np_dset
         
         allocate(metvar(np,met_nlon(iformat),met_nlat(iformat)))
         allocate(metscalar(np_dset))
         call shdf5_irec_f(ndims,idims,trim(met_vars(iformat,iv)),rvara = metscalar)
         
         if (np < np_dset) then     
            !----- No interpolation is needed, but we need two buffers. -------------------!
            do ip=1,np
               metvar(ip,:,:) = metscalar(ip)
            end do

         elseif(np == np_dset) then
            !----- Datasets are the same, no double allocation. ---------------------------!
            do ip=1,np
               metvar(ip,:,:) = metscalar(ip)
            end do

         else
            !------------------------------------------------------------------------------!
            !     The dataset buffer is smaller than the models buffer, so we recycle the  !
            ! last day.                                                                    !
            !------------------------------------------------------------------------------!
            do ip=1,np_dset
               metvar(ip,:,:) = metscalar(ip)
            end do

            do ip=np_dset+1,np
               !----- Then reuse the 28th for the 29th day. -------------------------------!
               metvar(ip,:,:) = metscalar(ip-points_per_day+1)
            end do
         end if

         deallocate(metscalar)

      case (4)
         !---------------------------------------------------------------------------------!
         !      Case 4, we don't read, instead we use the value stored in met_freq as the  !
         ! scalar value.                                                                   !
         !---------------------------------------------------------------------------------!
         allocate (metvar(1,met_nlon(iformat),met_nlat(iformat)))
         metvar(1,:,:) = met_frq(iformat,iv)
         
      case (2)
         !----- Case 2, 2D grid is read in, ie. np = 1. -----------------------------------!

         allocate(metvar(1,met_nlon(iformat),met_nlat(iformat)))

         call shdf5_info_f(trim(met_vars(iformat,iv)),ndims,idims)
         
         !----- Check to see if the dimensions are consistent. ----------------------------!
         if (ndims /= 2 .or. idims(1) /= met_nlon(iformat)                                 &
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
            call fatal_error('Mismatch between dataset and specified input'                &
                            ,'read_ol_file','ed_met_driver.f90')
         end if

         call shdf5_irec_f(ndims,idims,trim(met_vars(iformat,iv)),rvara = metvar)

      case default
         !----- Case 1 or 0, we must read time-and-space variable array. ------------------!
         call shdf5_info_f(trim(met_vars(iformat,iv)),ndims,idims)
         
         ! Check to see if the dimensions are consistent
         if (ndims /= 3 .or. idims(1).ne.np_dset .or. idims(2).ne.met_nlon(iformat)        &
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
            !----- No interpolation is needed, but we need two buffers. -------------------!
            allocate(metvar(np,met_nlon(iformat),met_nlat(iformat)))
            allocate(metvar_dset(np_dset,met_nlon(iformat),met_nlat(iformat)))
            
            call shdf5_irec_f(ndims, idims, trim(met_vars(iformat,iv)), rvara = metvar_dset)
            metvar(1:np,:,:) = metvar_dset(1:np,:,:)
            deallocate(metvar_dset)
            
            
         elseif (np == np_dset) then
            !----- Datasets are the same, no double allocation needed. --------------------!
            allocate(metvar(np,met_nlon(iformat),met_nlat(iformat)))
            call shdf5_irec_f(ndims, idims, trim(met_vars(iformat,iv)), rvara = metvar)

         else
            !------------------------------------------------------------------------------!
            !     The dataset buffer is smaller than the models buffer, so we recycle the  !
            ! last day.                                                                    !
            !------------------------------------------------------------------------------!
            allocate(metvar(np,met_nlon(iformat),met_nlat(iformat)))
            allocate(metvar_dset(np_dset,met_nlon(iformat),met_nlat(iformat)))

            call shdf5_irec_f(ndims,idims,trim(met_vars(iformat,iv)),rvara=metvar_dset)
            metvar(1:np_dset,:,:) = metvar_dset(1:np_dset,:,:)
            
            !----- Then reuse the 28th for the 29th day. ----------------------------------!
            metvar(np_dset+1:np,:,:) = metvar_dset(np_dset-points_per_day+1:np_dset,:,:)
            deallocate(metvar_dset)
         end if
      end select

      !----- Define some aliases for the indices. -----------------------------------------!
      ioa = offset + 1
      ioz = offset + np

      !------------------------------------------------------------------------------------!
      !      The dataset for that variable has been read into an array.  Now, the polygons !
      ! are cycled through, and they are associated with an index in the grid.  If this is !
      ! a lat-lon grid, the procedure is a simple interpolation between the domain end-    !
      ! points, if this is a polar stereo-graphic projection, we make use of the           !
      ! associated indices saved in memory, ie. (cgrid%ilon(ipy)).                         !
      !------------------------------------------------------------------------------------!
      do ipy = 1,cgrid%npolygons
         
         !----- Alias for longitude/latitude indices.  ------------------------------------!
         ilon = cgrid%ilon(ipy)
         ilat = cgrid%ilat(ipy)
         !---------------------------------------------------------------------------------!

         !----- Get the time series. ------------------------------------------------------!
         select case (trim(met_vars(iformat,iv)))
         case('nbdsf')
            cgrid%metinput(ipy)%nbdsf(ioa:ioz)     = metvar(1:np,ilon,ilat)
         case('nddsf')
            cgrid%metinput(ipy)%nddsf(ioa:ioz)     = metvar(1:np,ilon,ilat)
         case('vbdsf')
            cgrid%metinput(ipy)%vbdsf(ioa:ioz)     = metvar(1:np,ilon,ilat)
         case('vddsf')
            cgrid%metinput(ipy)%vddsf(ioa:ioz)     = metvar(1:np,ilon,ilat)
         case('prate')
            cgrid%metinput(ipy)%prate(ioa:ioz)     = metvar(1:np,ilon,ilat)
         case('dlwrf')
            cgrid%metinput(ipy)%dlwrf(ioa:ioz)     = metvar(1:np,ilon,ilat)
         case('pres')
            cgrid%metinput(ipy)%pres(ioa:ioz)      = metvar(1:np,ilon,ilat)
         case('hgt')
            cgrid%metinput(ipy)%hgt(ioa:ioz)       = metvar(1:np,ilon,ilat)
         case('ugrd')
            cgrid%metinput(ipy)%ugrd(ioa:ioz)      = metvar(1:np,ilon,ilat)
         case('vgrd')
            cgrid%metinput(ipy)%vgrd(ioa:ioz)      = metvar(1:np,ilon,ilat)
         case('ustar')
            cgrid%metinput(ipy)%atm_ustar(ioa:ioz) = metvar(1:np,ilon,ilat)
         case('sh')
            cgrid%metinput(ipy)%sh(ioa:ioz)        = metvar(1:np,ilon,ilat)
         case('tmp')
            cgrid%metinput(ipy)%tmp(ioa:ioz)       = metvar(1:np,ilon,ilat)
         case('co2')
            cgrid%metinput(ipy)%co2(ioa:ioz)       = metvar(1:np,ilon,ilat)
         end select
         
      end do
     
      !----- Deallocate the buffer. -------------------------------------------------------!
      deallocate(metvar)
     
      return
   end subroutine read_ol_file
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine transfers the dataset stored in the following month to the        !
   ! current month.                                                                        !
   !---------------------------------------------------------------------------------------!
   subroutine transfer_ol_month(vname, frq, cgrid)
     
      use ed_state_vars, only : edtype   ! ! variable type
      use consts_coms  , only : day_sec  ! ! intent(in)

      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(edtype)    , target     :: cgrid     ! Current grid
      character(len=*), intent(in) :: vname     ! Variable name
      real            , intent(in) :: frq       ! Frequency of update
      !----- Local variables. -------------------------------------------------------------!
      integer                      :: mem_size  ! Memory size
      integer                      :: ipy       ! Polygon counter
      integer                      :: ica       ! First element of current   month
      integer                      :: icz       ! Last  element of current   month
      integer                      :: ifa       ! First element of following month
      integer                      :: ifz       ! Last  element of following month
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !     The memory size is the size of each month.  The number of days is always 31    !
      ! because we allocate the maximum possible number of days.                           !
      !------------------------------------------------------------------------------------!
      mem_size = nint(day_sec / frq) * 31

      !------------------------------------------------------------------------------------!
      !     Compute the indices.                                                           !
      !------------------------------------------------------------------------------------!
      ica = 1
      icz = mem_size
      ifa = icz + 1
      ifz = 2 * mem_size

      !------------------------------------------------------------------------------------!
      !      Transfer the time series from following month to current month.               !
      !------------------------------------------------------------------------------------!
      select case (trim(vname))
      case ('nbdsf')
         do ipy=1,cgrid%npolygons
            cgrid%metinput(ipy)%nbdsf(ica:icz)      = cgrid%metinput(ipy)%nbdsf(ifa:ifz)
         end do

      case ('nddsf')
         do ipy=1,cgrid%npolygons
            cgrid%metinput(ipy)%nddsf(ica:icz)      = cgrid%metinput(ipy)%nddsf(ifa:ifz)
         end do

      case ('vbdsf')
         do ipy=1,cgrid%npolygons
            cgrid%metinput(ipy)%vbdsf(ica:icz)      = cgrid%metinput(ipy)%vbdsf(ifa:ifz)
         end do

      case ('vddsf')
         do ipy=1,cgrid%npolygons
            cgrid%metinput(ipy)%vddsf(ica:icz)      = cgrid%metinput(ipy)%vddsf(ifa:ifz)
         end do

      case ('prate')
         do ipy=1,cgrid%npolygons
            cgrid%metinput(ipy)%prate(ica:icz)      = cgrid%metinput(ipy)%prate(ifa:ifz)
         end do

      case ('dlwrf')
         do ipy=1,cgrid%npolygons
            cgrid%metinput(ipy)%dlwrf(ica:icz)      = cgrid%metinput(ipy)%dlwrf(ifa:ifz)
         end do

      case ('pres')
         do ipy=1,cgrid%npolygons
            cgrid%metinput(ipy)%pres(ica:icz)       = cgrid%metinput(ipy)%pres(ifa:ifz)
         end do

      case ('hgt')
         do ipy=1,cgrid%npolygons
            cgrid%metinput(ipy)%hgt(ica:icz)        = cgrid%metinput(ipy)%hgt(ifa:ifz)
         end do

      case ('ugrd')
         do ipy=1,cgrid%npolygons
            cgrid%metinput(ipy)%ugrd(ica:icz)       = cgrid%metinput(ipy)%ugrd(ifa:ifz)
         end do

      case ('vgrd')
         do ipy=1,cgrid%npolygons
            cgrid%metinput(ipy)%vgrd(ica:icz)       = cgrid%metinput(ipy)%vgrd(ifa:ifz)
         end do

      case ('ustar')
         do ipy=1,cgrid%npolygons
            cgrid%metinput(ipy)%atm_ustar(ica:icz)  = cgrid%metinput(ipy)%atm_ustar(ifa:ifz)
         end do

      case ('sh')
         do ipy=1,cgrid%npolygons
            cgrid%metinput(ipy)%sh(ica:icz)         = cgrid%metinput(ipy)%sh(ifa:ifz)
         end do

      case ('tmp')
         do ipy=1,cgrid%npolygons
            cgrid%metinput(ipy)%tmp(ica:icz)        = cgrid%metinput(ipy)%tmp(ifa:ifz)
         end do

      case ('co2')
         do ipy=1,cgrid%npolygons
            cgrid%metinput(ipy)%co2(ica:icz)        = cgrid%metinput(ipy)%co2(ifa:ifz)
         end do

      end select

      return
   end subroutine transfer_ol_month
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This sub-routine finds the met driver grid point that is the closest to the       !
   ! polygon of interest and is valid (i.e., met driver grid point is on land).            !
   !---------------------------------------------------------------------------------------!
   subroutine match_poly_grid(cgrid,nlon,nlat,lon,lat,land)
      use ed_state_vars  , only : edtype        ! ! structure
      use met_driver_coms, only : met_land_min  ! ! intent(in)
      

      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(edtype)              , target        :: cgrid
      integer                   , intent(in)    :: nlon
      integer                   , intent(in)    :: nlat
      real, dimension(nlon,nlat), intent(inout) :: lon
      real, dimension(nlon,nlat), intent(inout) :: lat
      real, dimension(nlon,nlat), intent(inout) :: land
      !----- Local variables --------------------------------------------------------------!
      integer                                   :: ilon
      integer                                   :: ilat
      integer                                   :: ipy
      logical                                   :: any_land
      real                                      :: min_dist
      real                                      :: this_dist
      !----- External function ------------------------------------------------------------!
      real                       , external     :: dist_gc
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !    Check that at least one point in the met driver is land.                        !
      !------------------------------------------------------------------------------------!
      any_land = any(land(:,:) >= met_land_min)
      if (.not. any_land) then
         call fatal_error("Met driver problem -- All grid points are over water."          &
                         ,"match_poly_grid","ed_met_driver.f90")
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Go through all polygons and grid points, and find the closest site that is not !
      ! over land.                                                                         !
      !------------------------------------------------------------------------------------!
      polyloop: do ipy = 1,cgrid%npolygons
         min_dist = huge (1.)

         lonloop: do ilon = 1,nlon
            latloop: do ilat = 1,nlat

               !------ Skip point in case it has not enough land in it. -------------------!
               if (land(ilon,ilat) < met_land_min) cycle latloop
               !---------------------------------------------------------------------------!

               !------ Make sure longitude goes from -180 to 180. -------------------------!
               if (lon(ilon,ilat) > 180.0) lon(ilon,ilat) = lon(ilon,ilat) - 360.0
               !---------------------------------------------------------------------------!

               !---------------------------------------------------------------------------!
               !    Find distance from met driver grid cell and polygon of interest.  In   !
               ! case this is the closest point so far, select it as the candidate.        !
               !---------------------------------------------------------------------------!
               this_dist = dist_gc(cgrid%lon(ipy), lon(ilon,ilat)                          &
                                  ,cgrid%lat(ipy), lat(ilon,ilat) )
               if (this_dist < min_dist) then
                  cgrid%ilon(ipy) = ilon
                  cgrid%ilat(ipy) = ilat
                  min_dist        = this_dist
               end if
               !---------------------------------------------------------------------------!
            end do latloop
            !------------------------------------------------------------------------------!
         end do lonloop
         !---------------------------------------------------------------------------------!
      end do polyloop
      !------------------------------------------------------------------------------------!

      return
   end subroutine match_poly_grid
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This sub-routine finds the indices of the input meteorological forcing that are   !
   ! the closest to each polygon, skipping those met driver points that are not over land  !
   ! in case a land/sea  mask is provided.                                                 !
   !---------------------------------------------------------------------------------------!
   subroutine get_lonlat_land(cgrid,iformat)
      use ed_state_vars  , only : edtype        ! ! structure
      use met_driver_coms, only : met_ll_header & ! intent(in)
                                , met_land_mask & ! intent(in)
                                , met_xmin      & ! intent(in)
                                , met_ymin      & ! intent(in)
                                , met_dx        & ! intent(in)
                                , met_dy        & ! intent(in)
                                , met_nlon      & ! intent(inout)
                                , met_nlat      & ! intent(inout)
                                , lon2d         & ! intent(out)
                                , lat2d         & ! intent(out)
                                , land2d        ! ! intent(out)
      use hdf5_utils     , only : shdf5_info_f  & ! sub-routine
                                , shdf5_irec_f  ! ! sub-routine
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(edtype)         , target         :: cgrid
      integer              , intent(in)     :: iformat
      !----- Local variables. -------------------------------------------------------------!
      integer                               :: ndims
      integer                               :: d
      integer                               :: ilon
      integer                               :: ilat
      integer              , dimension(3)   :: idims
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Check whether this data set contains longitude and latitude.  In case it does, !
      ! load it, otherwise assign the matrix.                                              !
      !------------------------------------------------------------------------------------!
      if (met_ll_header(iformat)) then


         !----- Allocate 2d arrays. -------------------------------------------------------!
         allocate (lon2d(met_nlon(iformat),met_nlat(iformat)))
         allocate (lat2d(met_nlon(iformat),met_nlat(iformat)))
         !---------------------------------------------------------------------------------!




         !---------------------------------------------------------------------------------!
         !    Coordinates are not provided in the HDF5 file.  Make the 2D matrices, keep-  !
         ! ing in mind that latitude is flipped.                                           !
         !---------------------------------------------------------------------------------!
         do ilat=1,met_nlat(iformat)
            do ilon=1,met_nlon(iformat)
               lon2d(ilon,ilat) = met_xmin(iformat) + real(ilon -  1) * met_dx(iformat)
               lat2d(ilon,ilat) = met_ymin(iformat)                                        &
                                + real(met_nlat(iformat) - ilat)      * met_dy(iformat)
            end do
         end do
         !---------------------------------------------------------------------------------!
      else
         !----- Get the dimension information from longitude and check it. ----------------!
         call shdf5_info_f('lon',ndims,idims)
         if(ndims /= 2) then
            write(unit=*,fmt="(a)") "Incorrect # of dimensions of variable ""lon""."
            write(unit=*,fmt="(a,1x,i5)") "NDIMS=",ndims
            do d=1,ndims
               write(unit=*,fmt="(a,1x,i5)") "---> ",d,": DIM=",idims(d)
            end do
            call fatal_error ("ED-2 cannot handle time-varying longitude."                 &
                             ,"get_lonlat_land","ed_met_driver.f90")
         end if
         !---------------------------------------------------------------------------------!



         !----- Allocate the longitude array, then read it. -------------------------------!
         allocate(lon2d(idims(1),idims(2)))
         call shdf5_irec_f(ndims, idims, 'lon',rvara = lon2d )
         !---------------------------------------------------------------------------------!

         !----- Save dimensions. ----------------------------------------------------------!
         met_nlon(iformat) = idims(1)
         met_nlat(iformat) = idims(2)
         !---------------------------------------------------------------------------------!



         !----- Get the dimension information from latitude and check it. -----------------!
         call shdf5_info_f('lat',ndims,idims)
         if(ndims /= 2) then
            write(unit=*,fmt="(a)") "Incorrect # of dimensions of variable ""lat""."
            write(unit=*,fmt="(a,1x,i5)") "NDIMS=",ndims
            do d=1,ndims
               write(unit=*,fmt='(a,1x,i5)') "---> ",d,": DIM=",idims(d)
            end do
            call fatal_error ("ED-2 cannot handle time-varying latitude."                  &
                             ,"get_lonlat_land","ed_met_driver.f90")
         end if
         !---------------------------------------------------------------------------------!



         !----- Stop in case longitude and latitude dimensions don't match. ---------------!
         if (met_nlon(iformat) /= idims(1) .or. met_nlat(iformat) /= idims(2)) then
            write(unit=*,fmt='(a,2(1x,i5))') " Dimensions of variable ""lon"": "           &
                                            ,met_nlon(iformat),met_nlat(iformat)
            write(unit=*,fmt='(a,2(1x,i5))') " Dimensions of variable ""lat"" : "          &
                                            ,idims(1),idims(2)
            call fatal_error ("Longitude and latitude dimensions should match."            &
                             ,"get_lonlat_land","ed_met_driver.f90")
         end if
         !---------------------------------------------------------------------------------!



         !----- Allocate the longitude array, then read it. -------------------------------!
         allocate(lat2d(idims(1),idims(2)))
         call shdf5_irec_f(ndims, idims, 'lat',rvara = lat2d )
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Check whether this data set contains land/sea mask.  In case it does, load it, !
      ! otherwise assume that everything is land.                                          !
      !------------------------------------------------------------------------------------!
      if (met_land_mask(iformat)) then
         !----- Get the dimension information from land/sea mask and check it. ------------!
         call shdf5_info_f('land',ndims,idims)
         if(ndims /= 2) then
            write(unit=*,fmt="(a)") "Incorrect # of dimensions of variable ""land""."
            write(unit=*,fmt="(a,1x,i5)") "NDIMS=",ndims
            do d=1,ndims
               write(unit=*,fmt="(a,1x,i5)") "---> ",d,": DIM=",idims(d)
            end do
            call fatal_error ("ED-2 cannot handle time-varying land/sea mask."             &
                             ,"get_lonlat_land","ed_met_driver.f90")
         end if
         !---------------------------------------------------------------------------------!



         !----- Stop in case longitude and latitude dimensions don't match. ---------------!
         if (met_nlon(iformat) /= idims(1) .or. met_nlat(iformat) /= idims(2)) then
            write(unit=*,fmt='(a,2(1x,i5))') " Dimensions of variables ""lon/lat"": "      &
                                            ,met_nlon(iformat),met_nlat(iformat)
            write(unit=*,fmt='(a,2(1x,i5))') " Dimensions of variable ""land"" : "         &
                                            ,idims(1),idims(2)
            call fatal_error ("Land mask and lon/lat dimensions should match."             &
                             ,"get_lonlat_land","ed_met_driver.f90")
         end if
         !---------------------------------------------------------------------------------!



         !----- Allocate the longitude array, then read it. -------------------------------!
         allocate(land2d(idims(1),idims(2)))
         call shdf5_irec_f(ndims, idims, 'land',rvara = land2d )
         !---------------------------------------------------------------------------------!

      else

         !----- Allocate 2d arrays. -------------------------------------------------------!
         allocate (land2d(met_nlon(iformat),met_nlat(iformat)))
         land2d(:,:) = 1.0
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !     Look for the met driver point that is on land and is the closest to the        !
      ! polygon of interest.                                                               !
      !------------------------------------------------------------------------------------!
      call match_poly_grid(cgrid,met_nlon(iformat),met_nlat(iformat),lon2d,lat2d,land2d)
      !------------------------------------------------------------------------------------!

      !------ Deallocate matrices. --------------------------------------------------------!
      deallocate(lat2d,lon2d,land2d)
      !------------------------------------------------------------------------------------!

      return
   end subroutine get_lonlat_land
   !=======================================================================================!
   !=======================================================================================!

end module ed_met_driver
!==========================================================================================!
!==========================================================================================!
