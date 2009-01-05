subroutine read_met_driver_head()
  use max_dims, only: max_met_vars
  use met_driver_coms, only: nformats, met_names, met_nlon,   &
       met_nlat, met_dx, met_dy, met_xmin, met_ymin, met_nv,   &
       met_vars, met_frq, met_interp, ed_met_driver_db, no_ll
  implicit none
  logical :: l1
  logical :: yes_lat     ! Logical for determining if latitude grids are present
  logical :: yes_lon     ! Logical for determining if longitude grids are present
  integer :: iformat,n

  inquire(file=trim(ed_met_driver_db),exist=l1)
  if(.not.l1)then
     write (unit=*,fmt='(a)') 'File '//trim(ed_met_driver_db)//' not found!'
     write (unit=*,fmt='(a)') 'Specify ED_MET_DRIVER_DB properly in ED namelist.'
     call fatal_error('Ed_met_driver_db not found!' &
                     ,'read_met_driver_head'        &
                     ,'ed_met_driver.f90')
  endif

  open(unit=12,file=trim(ed_met_driver_db),form='formatted',status='old')
  read(unit=12,fmt=*)  ! skip header

  ! Read the number of different file formats
  read(unit=12,fmt=*) nformats

  ! Allocate the header information for each format
  allocate(met_names(nformats))
  allocate(met_nlon(nformats))
  allocate(met_nlat(nformats))
  allocate(met_dx(nformats))
  allocate(met_dy(nformats))
  allocate(met_xmin(nformats))
  allocate(met_ymin(nformats))
  allocate(met_nv(nformats))
  allocate(met_vars(nformats, max_met_vars))
  allocate(met_frq(nformats, max_met_vars))
  allocate(met_interp(nformats, max_met_vars))
  
  no_ll = .true. ! Just to initialize, if lon/lat are both found, it will become .false.
  
  ! Read the information for each format
  do iformat = 1,nformats
     read(unit=12,fmt='(a)')  met_names(iformat)
     read(unit=12,fmt=*)      met_nlon(iformat), met_nlat(iformat), met_dx(iformat)        &
                            , met_dy(iformat)  , met_xmin(iformat), met_ymin(iformat)
     read(unit=12,fmt=*)      met_nv(iformat)
     read(unit=12,fmt=*)      (met_vars(iformat,n)  ,n=1,met_nv(iformat))
     read(unit=12,fmt=*)      (met_frq(iformat,n)   ,n=1,met_nv(iformat))
     read(unit=12,fmt=*)      (met_interp(iformat,n),n=1,met_nv(iformat))
     
     !Just making sure that the variable list is case insensitive.
     call tolower(met_vars(iformat,1:met_nv(iformat)),met_nv(iformat))
     
     ! First check - see if lat/lon data are there
     yes_lon = any(met_vars(iformat,1:met_nv(iformat)) == 'lon')
     yes_lat = any(met_vars(iformat,1:met_nv(iformat)) == 'lat')
     

     ! Check to see if you have both, none or one of each
     if (yes_lat .and. yes_lon) then
        no_ll = .false.
     elseif (yes_lat .neqv. yes_lon) then
        call fatal_error("You are missing a lat or a lon variable in the met nl"           &
                        ,'read_met_driver_head'                                            &
                        ,'ed_met_driver.f90')
     endif
  end do

  close (unit=12,status='keep')

  return
end subroutine read_met_driver_head
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine init_met_drivers_array

  use max_dims, only: max_met_vars
  use met_driver_coms, only: nformats, met_names, met_nlon,   &
       met_nlat, met_dx, met_dy, met_xmin, met_ymin, met_nv,   &
       met_vars, met_frq, met_interp, ed_met_driver_db,no_ll,have_co2
  use ed_state_vars,only : edgrid_g,edtype,polygontype
  use grid_coms,only:ngrids

  implicit none

  type(polygontype), pointer :: cpoly
  type(edtype), pointer      :: cgrid

  integer :: ipy,igr
  integer :: iformat
  integer :: iv,n
  integer :: mem_size
  logical :: l1

  ! Set the lapse rates
  
  do igr=1,ngrids
     call setLapseParms_ar(edgrid_g(igr))
  enddo

  ! Read the information for each format
  have_co2=.false.
  do iformat = 1,nformats

     do igr = 1,ngrids

        cgrid => edgrid_g(igr)

        do ipy = 1,cgrid%npolygons

           cpoly => cgrid%polygon(ipy)

           ! Make sure site falls within file domain
           if( no_ll .and. &
                (cgrid%lon(ipy) < (met_xmin(iformat) - 0.5 * met_dx(iformat)) .or.  &
                cgrid%lat(ipy) < (met_ymin(iformat) - 0.5 * met_dy(iformat)) .or.  &
                cgrid%lon(ipy) > (met_xmin(iformat) + (met_nlon(iformat)-1) *  &
                met_dx(iformat) + 0.5 * met_dx(iformat)) .or.  &
                cgrid%lat(ipy) > (met_ymin(iformat) + (met_nlat(iformat)-1) *  &
                met_dy(iformat) + 0.5 * met_dy(iformat)) ))then
              print*
              print*,'========================================================'
              print*,'Site is not within the domain of the meteorological drivers'
              print*,'========================================================'
              print*,iformat,trim(met_names(iformat))
              print*,cgrid%lon(ipy), met_xmin(iformat) - 0.5 * met_dx(iformat)
              print*,cgrid%lat(ipy), met_ymin(iformat) - 0.5 * met_dy(iformat)
              print*,cgrid%lon(ipy), met_xmin(iformat) + (met_nlon(iformat)-1) *  &
                   met_dx(iformat) + 0.5 * met_dx(iformat)
              print*,cgrid%lat(ipy), met_ymin(iformat) + (met_nlat(iformat)-1) *  &
                   met_dy(iformat) + 0.5 * met_dy(iformat)
              stop
           endif
           
           ! Allocate memory
           do iv = 1,met_nv(iformat)
              
              !  Either this is a constant or variable in time
              !  Constant - type 2 and 4
              if ( met_interp(iformat,iv) == 2 .or. &
                   met_interp(iformat,iv) == 4 ) then
                 
                 mem_size = 1
              
                 !  Variable - type 0,1,3
              else
                 ! Calculate number of points in a month
                 mem_size = nint(86400.0 / met_frq(iformat,iv)) * 31
              endif
              
              ! If this is an interpolated variable, read two months in.
              if(met_interp(iformat,iv) == 1)mem_size = 2 * mem_size
              if (trim(met_vars(iformat,iv)) == 'co2') have_co2=.true.
              select case (trim(met_vars(iformat,iv)))
              case ('nbdsf')
                 nullify(cgrid%metinput(ipy)%nbdsf)
                 allocate(cgrid%metinput(ipy)%nbdsf(mem_size))
                 cgrid%metinput(ipy)%nbdsf = huge(1.)
              case ('nddsf')
                 nullify(cgrid%metinput(ipy)%nddsf)
                 allocate(cgrid%metinput(ipy)%nddsf(mem_size))
                 cgrid%metinput(ipy)%nddsf = huge(1.)
              case ('vbdsf')
                 nullify(cgrid%metinput(ipy)%vbdsf)
                 allocate(cgrid%metinput(ipy)%vbdsf(mem_size))
                 cgrid%metinput(ipy)%vbdsf = huge(1.)
              case ('vddsf')
                 nullify(cgrid%metinput(ipy)%vddsf)
                 allocate(cgrid%metinput(ipy)%vddsf(mem_size))
                 cgrid%metinput(ipy)%vddsf = huge(1.)
              case ('prate')
                 nullify(cgrid%metinput(ipy)%prate)
                 allocate(cgrid%metinput(ipy)%prate(mem_size))
                 cgrid%metinput(ipy)%prate = huge(1.)
              case ('dlwrf')
                 nullify(cgrid%metinput(ipy)%dlwrf)
                 allocate(cgrid%metinput(ipy)%dlwrf(mem_size))
                 cgrid%metinput(ipy)%dlwrf = huge(1.)
              case ('pres')
                 nullify(cgrid%metinput(ipy)%pres)
                 allocate(cgrid%metinput(ipy)%pres(mem_size))
                 cgrid%metinput(ipy)%pres = huge(1.)
              case ('hgt')
                 nullify(cgrid%metinput(ipy)%hgt)
                 allocate(cgrid%metinput(ipy)%hgt(mem_size))
                 cgrid%metinput(ipy)%hgt = huge(1.)
              case ('ugrd')
                 nullify(cgrid%metinput(ipy)%ugrd)
                 allocate(cgrid%metinput(ipy)%ugrd(mem_size))
                 cgrid%metinput(ipy)%ugrd = huge(1.)
              case ('vgrd')
                 nullify(cgrid%metinput(ipy)%vgrd)
                 allocate(cgrid%metinput(ipy)%vgrd(mem_size))
                 cgrid%metinput(ipy)%vgrd = huge(1.)
              case ('sh')
                 nullify(cgrid%metinput(ipy)%sh)
                 allocate(cgrid%metinput(ipy)%sh(mem_size))
                 cgrid%metinput(ipy)%sh = huge(1.)
              case ('tmp')
                 nullify(cgrid%metinput(ipy)%tmp)
                 allocate(cgrid%metinput(ipy)%tmp(mem_size))
                 cgrid%metinput(ipy)%tmp = huge(1.)
              case ('co2')
                 have_co2=.true.
                 nullify(cgrid%metinput(ipy)%co2)
                 allocate(cgrid%metinput(ipy)%co2(mem_size))
                 cgrid%metinput(ipy)%co2 = huge(1.)
              case ('lat','lon')
                 ! Do nothing
              case default
                 call fatal_error('Invalid met variable'//trim(met_vars(iformat,iv))//'!'  &
                                 ,'init_met_drivers_array'                                 &
                                 ,'ed_met_driver_array.f90')
              end select

           end do
        end do
     end do
  end do
  return
end subroutine init_met_drivers_array
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine read_met_drivers_init_array

  use ed_state_vars,only:edgrid_g,edtype,polygontype
  use hdf5_utils, only: shdf5_open_f, shdf5_close_f
  use met_driver_coms, only: nformats, met_names, met_nv, met_interp,  &
       met_frq, metcyc1, metcycf
  use mem_sites, only: grid_type
  use misc_coms, only: current_time
  use grid_coms,only: ngrids

  implicit none

  type(polygontype), pointer :: cpoly
  type(edtype), pointer      :: cgrid

  integer :: igr,ipy,isi
  integer :: year_use
  integer :: ncyc
  integer :: iformat
  character(len=256) :: infile
  character(len=3), dimension(12), parameter :: mname = (/'JAN', 'FEB',   &
       'MAR', 'APR', 'MAY', 'JUN', 'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC'/) 
  logical :: exans
  integer :: iv
  integer :: offset
  integer :: m2
  integer :: y2
  integer :: year_use_2
  logical,external :: isleap
  integer :: nday

  ! If we need to recycle over years, find the appropriate year to apply.
  year_use = current_time%year
  ncyc = metcycf - metcyc1 + 1

  ! If we are after the last year...
  do while(year_use > metcycf)
     year_use = year_use - ncyc
  enddo

  ! If we are before the first year...
  do while(year_use < metcyc1)
     year_use = year_use + ncyc
  enddo


  do igr = 1,ngrids

     cgrid => edgrid_g(igr)
     
     ! Loop over the different file formats
     do iformat = 1, nformats
        
        ! Open the file
        write(infile,'(a,i4.4,a,a)')trim(met_names(iformat)), year_use,   &
             mname(current_time%month),'.h5'
        inquire(file=trim(infile),exist=exans)
        if(exans)then
           call shdf5_open_f(trim(infile),'R')
        else
           print*,'Cannot open met driver input file',trim(infile)
           stop
        endif
        
        ! The following subroutine determines grid indices
        ! of each polygon's match to the met data
        
        call getll_array(cgrid,iformat)
        
        ! Loop over variables
        do iv = 1, met_nv(iformat)

           offset = 0
           call read_ol_file_ar(iformat, iv, year_use, mname(current_time%month),  &
                current_time%year, offset, cgrid)

        enddo

        ! Close the HDF5 file.
        call shdf5_close_f()
        
        ! For all interpolated variables, we also need the next time.
        
        ! Find next month and year
        m2 = current_time%month + 1
        y2 = current_time%year
        year_use_2 = year_use
        
        ! If this takes us into the next year, increment year and 
        ! reset month to January.
        if(m2 == 13)then
           m2 = 1
           y2 = current_time%year + 1
           year_use_2 = y2
           
           ! If we are now after the last year...
           do while(year_use_2 > metcycf)
              year_use_2 = year_use_2 - ncyc
           enddo
           
           ! If we are now before the first year...
           do while(year_use_2 < metcyc1)
              year_use_2 = year_use_2 + ncyc
           enddo
        endif
        
        ! Now, open the file once.
        write(infile,'(a,i4.4,a,a)')trim(met_names(iformat)), year_use_2,   &
             mname(m2),'.h5'
        inquire(file=trim(infile),exist=exans)
        if(exans)then
           call shdf5_open_f(trim(infile),'R')
        else
           print*,'Cannot open met driver input file',trim(infile)
           stop
        endif
        
        ! Loop over variables
        do iv = 1, met_nv(iformat)
           
           if(met_interp(iformat,iv) == 1)then

              offset = nint(86400.0 / met_frq(iformat,iv)) * 31

              ! Read the file.

              call read_ol_file_ar(iformat, iv, year_use_2, mname(m2),  &
                   y2, offset, cgrid)
              
           endif
        enddo
        
        ! Close the HDF5 file.
        call shdf5_close_f()
        
     enddo

  enddo
  

!  stop

  return
end subroutine read_met_drivers_init_array
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine read_met_drivers_array

  use ed_state_vars,only:edgrid_g,edtype,polygontype
  use mem_sites, only: grid_type
  use met_driver_coms, only: nformats, met_names, met_nv, met_interp,  &
       met_frq, met_vars, metcyc1, metcycf
  use misc_coms, only: current_time
  use hdf5_utils, only: shdf5_open_f, shdf5_close_f
  use ed_state_vars,only:edtype,edgrid_g
  use grid_coms,only : ngrids

  implicit none
  type(edtype),pointer :: cgrid
  integer :: igr
  integer :: year_use
  integer :: ncyc
  integer :: iformat
  character(len=256) :: infile
  character(len=3), dimension(12), parameter :: mname = (/'JAN', 'FEB',   &
       'MAR', 'APR', 'MAY', 'JUN', 'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC'/) 
  logical :: exans
  integer :: iv
  integer :: offset
  integer :: m2
  integer :: y2
  integer :: year_use_2
  integer :: nday
  logical, external :: isleap

  ! If we need to recycle over years, find the appropriate year to apply.
  year_use = current_time%year
  ncyc = metcycf - metcyc1 + 1

  ! If we are after the last year...
  do while(year_use > metcycf)
     year_use = year_use - ncyc
  enddo

  ! If we are before the first year...
  do while(year_use < metcyc1)
     year_use = year_use + ncyc
  enddo

  do igr=1,ngrids

     cgrid => edgrid_g(igr)

     ! Loop over the different file formats
     do iformat = 1, nformats
        
        ! Open the file
        write(infile,'(a,i4.4,a,a)')trim(met_names(iformat)), year_use,   &
             mname(current_time%month),'.h5'

        inquire(file=trim(infile),exist=exans)
        if(exans)then
           call shdf5_open_f(trim(infile),'R')
        else
           print*,'Cannot open met driver input file',trim(infile)
           stop
        endif
        
        !  Get the mapping between the polygons and the gridded variables
        call getll_array(cgrid,iformat)
        
        ! Loop over variables
        do iv = 1, met_nv(iformat)
           
           ! See if it is an interpolation variable.
           if(met_interp(iformat,iv) /= 1)then
              
              ! If not, things are simple.  Just read in the month.
              offset = 0
              call read_ol_file_ar(iformat, iv, year_use,   &
                   mname(current_time%month), current_time%year, offset,cgrid)

           else
              
              ! Here, just transfer future to current month.
              call transfer_ol_month_ar(trim(met_vars(iformat,iv)),   &
                   met_frq(iformat,iv),cgrid)
              
           endif
           
        enddo

        ! Close the HDF5 file.
        call shdf5_close_f()
        
        
        ! For all interpolated variables, get the future month.
        
        ! Find next month and year
        m2 = current_time%month + 1
        y2 = current_time%year
        year_use_2 = year_use
        
        ! If this takes us into the next year, increment year and 
        ! reset month to January.
        if(m2 == 13)then
           m2 = 1
           y2 = current_time%year + 1
           year_use_2 = y2
           
           ! If we are now after the last year...
           do while(year_use_2 > metcycf)
              year_use_2 = year_use_2 - ncyc
           enddo
           
           ! If we are now before the first year...
           do while(year_use_2 < metcyc1)
              year_use_2 = year_use_2 + ncyc
           enddo
        endif
        
        ! Now, open the file once.
        write(infile,'(a,i4.4,a,a)')trim(met_names(iformat)), year_use_2,   &
             mname(m2),'.h5'
        inquire(file=trim(infile),exist=exans) 
        if(exans)then
           call shdf5_open_f(trim(infile),'R')
        else
           print*,'Cannot open met driver input file',trim(infile)
           stop
        endif
     
        ! Loop over variables
        do iv = 1, met_nv(iformat)
           
           if(met_interp(iformat,iv) == 1)then

              offset = nint(86400.0 / met_frq(iformat,iv)) * 31
              call read_ol_file_ar(iformat, iv, year_use_2, mname(m2),  &
                   y2, offset,cgrid)
              
           endif
        enddo
        
        ! Close the HDF5 file.
        call shdf5_close_f()

     enddo
  enddo

  return
end subroutine read_met_drivers_array
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine update_met_drivers_array(cgrid)
  
  use ed_state_vars,only:edtype,polygontype
  use met_driver_coms, only: met_frq, nformats, met_nv, met_interp, met_vars &
                            ,have_co2,initial_co2
  use misc_coms, only: current_time
  use consts_coms, only: rdry, cice, cliq, alli, rocp, p00, cp,day_sec,t3ple
  use therm_lib, only: virtt

  implicit none

  type(edtype),target :: cgrid
  type(polygontype), pointer :: cpoly
  integer :: points_per_day
  integer :: nday
  integer :: np
  integer :: ndays_elapsed
  integer :: nseconds_elapsed
  integer :: mlo
  integer :: mhi
  real :: t1
  real :: t2
  integer :: iformat
  integer :: iv
  real    :: hillshade
  integer :: ipy,isi
  logical, external :: isleap

  ! Initialize vels
  do ipy=1,cgrid%npolygons
     cgrid%met(ipy)%vels = 0.0
  enddo
  


  do iformat = 1, nformats
     
     ! Loop over variables
     do iv = 1, met_nv(iformat)
        
        if(met_interp(iformat,iv) /= 1)then
           
           ! In these cases there is no time dynamic, simply set mlo to 1
           if(met_interp(iformat,iv) == 2 .or. met_interp(iformat,iv) == 4 ) then
              
              mlo = 1
              
              ! In these cases there is no interpolation, but there are
              ! time varying data.
           elseif(met_interp(iformat,iv) == 0 .or. met_interp(iformat,iv) == 3) then
              
              ! If this is not an interpolation variable, just find the time
              ! point.
              ndays_elapsed = current_time%date - 1
              nseconds_elapsed = nint(current_time%time) + ndays_elapsed * 86400
              mlo = int(float(nseconds_elapsed) / met_frq(iformat,iv)) + 1
              
           endif
           
           ! Find which variable it is, and then fill the sites.
           select case (trim(met_vars(iformat,iv)))
           case('nbdsf')
              do ipy = 1,cgrid%npolygons
                 cgrid%met(ipy)%nir_beam = cgrid%metinput(ipy)%nbdsf(mlo)
              enddo

           case('nddsf')

              do ipy = 1,cgrid%npolygons
                 cgrid%met(ipy)%nir_diffuse = cgrid%metinput(ipy)%nddsf(mlo)
              enddo

           case('vbdsf')

              do ipy = 1,cgrid%npolygons
                 cgrid%met(ipy)%par_beam = cgrid%metinput(ipy)%vbdsf(mlo)
              enddo

           case('vddsf')
              
              do ipy = 1,cgrid%npolygons
                 cgrid%met(ipy)%par_diffuse = cgrid%metinput(ipy)%vddsf(mlo)
              enddo

           case('prate')
              
              do ipy = 1,cgrid%npolygons
                 cgrid%met(ipy)%pcpg = cgrid%metinput(ipy)%prate(mlo)
              enddo

           case('dlwrf')
              
              do ipy = 1,cgrid%npolygons
                 cgrid%met(ipy)%rlong = cgrid%metinput(ipy)%dlwrf(mlo)
              enddo

           case('pres')
              
              do ipy = 1,cgrid%npolygons
                 cgrid%met(ipy)%prss = cgrid%metinput(ipy)%pres(mlo)
              enddo

           case('hgt')
              
              do ipy = 1,cgrid%npolygons
                 cgrid%met(ipy)%geoht = cgrid%metinput(ipy)%hgt(mlo)
              enddo

           case('ugrd')
              
              do ipy = 1,cgrid%npolygons
                 cgrid%met(ipy)%vels = cgrid%met(ipy)%vels +   &
                      cgrid%metinput(ipy)%ugrd(mlo)**2
              enddo

           case('vgrd')
              
              do ipy = 1,cgrid%npolygons
                 cgrid%met(ipy)%vels = cgrid%met(ipy)%vels +   &
                      cgrid%metinput(ipy)%vgrd(mlo)**2
              enddo

           case('sh')
              
              do ipy = 1,cgrid%npolygons
                 cgrid%met(ipy)%atm_shv = cgrid%metinput(ipy)%sh(mlo)
              enddo

           case('tmp')
              
              do ipy = 1,cgrid%npolygons
                 cgrid%met(ipy)%atm_tmp = cgrid%metinput(ipy)%tmp(mlo)
              enddo
              
           case('co2')
              
              do ipy = 1,cgrid%npolygons
                 cgrid%met(ipy)%atm_co2 = cgrid%metinput(ipy)%co2(mlo)
              enddo
           end select
           
        else
           
           ! In this case, we need to interpolate.
           ndays_elapsed = current_time%date - 1
           nseconds_elapsed = nint(current_time%time) + ndays_elapsed * day_sec
           
           ! First, get the number of points per day and per month.
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
           
           
           ! If this is not an interpolation variable, just find the time
           ! point.
           ndays_elapsed = current_time%date - 1
           nseconds_elapsed = nint(current_time%time) + ndays_elapsed * day_sec
           mlo = int(float(nseconds_elapsed) / met_frq(iformat,iv)) + 1
           
           
           ! Get indices
           mlo = int(float(nseconds_elapsed)/met_frq(iformat,iv)) + 1
           mhi = mlo + 1
           
           if(mhi > np)then
              mhi = 1 + nint(86400.0 / met_frq(iformat,iv)) * 31
           endif
           
           ! Get interpolation factors
           t1 = mod(float(nseconds_elapsed), met_frq(iformat,iv)) /   &
                met_frq(iformat,iv)
           t2 = 1.0 - t1

           ! Find the variable and fill the sites.
           select case (trim(met_vars(iformat,iv)))
           case('prate')

              do ipy = 1,cgrid%npolygons
                 cgrid%met(ipy)%pcpg = cgrid%metinput(ipy)%prate(mhi) *   &
                      t1 + cgrid%metinput(ipy)%prate(mlo) * t2
              enddo

           case('dlwrf')

              do ipy = 1,cgrid%npolygons
                 cgrid%met(ipy)%rlong = cgrid%metinput(ipy)%dlwrf(mhi) * t1 +  &
                      cgrid%metinput(ipy)%dlwrf(mlo) * t2
              enddo

           case('pres')

              do ipy = 1,cgrid%npolygons
                 cgrid%met(ipy)%prss = cgrid%metinput(ipy)%pres(mhi) * t1 +   &
                      cgrid%metinput(ipy)%pres(mlo) * t2
              enddo

           case('hgt')

              do ipy = 1,cgrid%npolygons
                 cgrid%met(ipy)%geoht = cgrid%metinput(ipy)%hgt(mhi) * t1 +  &
                      cgrid%metinput(ipy)%hgt(mlo) * t2
              enddo

           case('ugrd')

              do ipy = 1,cgrid%npolygons
                 cgrid%met(ipy)%vels = cgrid%met(ipy)%vels +   &
                      (cgrid%metinput(ipy)%ugrd(mhi) * t1 +   &
                      cgrid%metinput(ipy)%ugrd(mlo) * t2)**2
              enddo

           case('vgrd')

              do ipy = 1,cgrid%npolygons
                 cgrid%met(ipy)%vels = cgrid%met(ipy)%vels +   &
                      (cgrid%metinput(ipy)%vgrd(mhi) * t1 +  &
                      cgrid%metinput(ipy)%vgrd(mlo) * t2)**2
              enddo
           case('sh')

              do ipy = 1,cgrid%npolygons
                 cgrid%met(ipy)%atm_shv = cgrid%metinput(ipy)%sh(mhi) *   &
                      t1 + cgrid%metinput(ipy)%sh(mlo) * t2
              enddo

           case('tmp')

              do ipy = 1,cgrid%npolygons
                 cgrid%met(ipy)%atm_tmp = cgrid%metinput(ipy)%tmp(mhi) *   &
                      t1 + cgrid%metinput(ipy)%tmp(mlo) * t2
              enddo

           case('nbdsf')

              do ipy = 1,cgrid%npolygons
                 cgrid%met(ipy)%nir_beam = cgrid%metinput(ipy)%nbdsf(mhi) *  &
                      t1 + cgrid%metinput(ipy)%nbdsf(mlo) * t2
              enddo

           case('nddsf')

              do ipy = 1,cgrid%npolygons
                 cgrid%met(ipy)%nir_diffuse =   &
                      cgrid%metinput(ipy)%nddsf(mhi) * t1 +  &
                      cgrid%metinput(ipy)%nddsf(mlo) * t2
              enddo

           case('vbdsf')
              
              do ipy = 1,cgrid%npolygons
                 cgrid%met(ipy)%par_beam = cgrid%metinput(ipy)%vbdsf(mhi) * &
                      t1 + cgrid%metinput(ipy)%vbdsf(mlo) * t2
              enddo

           case('vddsf')

              do ipy = 1,cgrid%npolygons
                 cgrid%met(ipy)%par_diffuse =   &
                      cgrid%metinput(ipy)%vddsf(mhi) * t1 +  &
                      cgrid%metinput(ipy)%vddsf(mlo) * t2
              enddo

           case('co2')
              
              do ipy = 1,cgrid%npolygons
                 cgrid%met(ipy)%atm_co2 = cgrid%metinput(ipy)%co2(mhi) * t1 +  &
                      cgrid%metinput(ipy)%co2(mlo) * t2
              enddo

           end select
           
        endif
        
     enddo
  enddo
  
  ! Change from velocity squared to velocity, get the rhos, and compute
  ! qpcpg and dpcpg.

  do ipy = 1,cgrid%npolygons
        
     ! co2
     if (.not.have_co2) then
        cgrid%met(ipy)%atm_co2 = initial_co2
     end if

     !! APPLY MET TO SITES
     !! AND ADJUST MET VARIABLES FOR TOPOGRAPHY
     call calc_met_lapse_ar(cgrid,ipy)

     cpoly => cgrid%polygon(ipy)

     do isi = 1,cpoly%nsites
        
        ! vels
        cpoly%met(isi)%vels = sqrt(max(0.0,cpoly%met(isi)%vels))
        cpoly%met(isi)%vels_stab = max(0.1,cpoly%met(isi)%vels)
        cpoly%met(isi)%vels_unstab = max(1.0,cpoly%met(isi)%vels_stab)
        
        ! co2
        if (.not.have_co2) cpoly%met(isi)%atm_co2 = initial_co2
        
        ! exner
        cpoly%met(isi)%exner = cp * (cpoly%met(isi)%prss / p00)**rocp
        ! solar radiation
        cpoly%met(isi)%rshort_diffuse = cpoly%met(isi)%par_diffuse +   &
             cpoly%met(isi)%nir_diffuse
        
        cpoly%met(isi)%rshort = cpoly%met(isi)%rshort_diffuse +   &
             cpoly%met(isi)%par_beam + cpoly%met(isi)%nir_beam
        
        
        ! rho
        cpoly%met(isi)%rhos = cpoly%met(isi)%prss &
                            / (rdry * virtt(cpoly%met(isi)%atm_tmp,cpoly%met(isi)%atm_shv))
        
        ! qpcpg, dpcpg
        if(cpoly%met(isi)%atm_tmp > t3ple)then
           cpoly%met(isi)%qpcpg = (cliq * cpoly%met(isi)%atm_tmp + alli) * cpoly%met(isi)%pcpg
           cpoly%met(isi)%dpcpg = max(0.0, cpoly%met(isi)%pcpg * 0.001)
        else
           cpoly%met(isi)%qpcpg = cice * cpoly%met(isi)%atm_tmp * cpoly%met(isi)%pcpg
           cpoly%met(isi)%dpcpg = max(0.0, cpoly%met(isi)%pcpg * 0.01)
        endif
        !! note: dpcpg currently never gets used
        !! snow density is calculated in the integrator (grep snowdens)
        
     enddo
          
  enddo

  return
end subroutine update_met_drivers_array
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine read_ol_file_ar(iformat, iv, year_use, mname, year, offset, cgrid)

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
  integer :: ipy,isi
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
  integer :: hdferr
 


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
        stop
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
        print*,ndims,np_dset,met_nlon(iformat),met_nlat(iformat),idims
        stop
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
end subroutine read_ol_file_ar
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine transfer_ol_month_ar(vname, frq, cgrid)
  
  use ed_state_vars,only : edtype

  implicit none

  type(edtype),target :: cgrid
  character*(*) :: vname
  real, intent(in) :: frq
  integer :: mem_size
  integer :: ipy

  mem_size = nint(86400.0 / frq) * 31

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
end subroutine transfer_ol_month_ar
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine match_poly_grid_array(cgrid,nlon,nlat,lon,lat)

  use ed_state_vars,only:edtype

  implicit none
  
  type(edtype),target :: cgrid

  integer :: nlon,nlat,ilat,ilon,offlat,offlon
  integer :: located
  integer :: ipy
  real,dimension(nlon,nlat) :: lon,lat
  real :: deltax,deltay


  do ipy = 1,cgrid%npolygons
     
     located = 0
     do ilon = 1,nlon

        do ilat = 1,nlat
           
           if (ilat.ne.nlat) then
              offlat = 1
           else
              offlat = -1
           endif
           if (ilon.ne.nlon) then
              offlon = 1
           else
              offlon = -1
           endif

           if (lon(ilon,ilat)>180) lon(ilon,ilat) = lon(ilon,ilat) - 360.0
           if (lon(ilon+offlon,ilat)>180) lon(ilon+offlon,ilat) = lon(ilon+offlon,ilat) - 360.0
           if (lon(ilon,ilat+offlat)>180) lon(ilon,ilat+offlat) = lon(ilon,ilat+offlat) - 360.0

           deltay = 0.50*abs( lat(ilon,ilat+offlat) - lat(ilon,ilat) )
           deltax = 0.50*abs( lon(ilon+offlon,ilat) - lon(ilon,ilat) )
           
           if(  abs(cgrid%lat(ipy) - lat(ilon,ilat)).lt.deltay .and. &
                abs(cgrid%lon(ipy) - lon(ilon,ilat)).lt.deltax  ) then

              if(located.eq.1) then
                 print*,"OVERLAP PROBLEM"
                 stop
              endif
              
              cgrid%ilon(ipy) = ilon
              cgrid%ilat(ipy) = ilat
              located = 1
              
           endif
           
        enddo
     enddo

     if (located.eq.0) then
        print*,"Could not find a corresponding grid space"
        print*,"ixn the met data for polygon:",ipy
        print*,cgrid%lat(ipy),cgrid%lon(ipy)

        print*,lat
        print*,""
        print*,lon

        stop
     endif
     
  enddo

  return
end subroutine match_poly_grid_array
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine getll_array(cgrid,iformat)

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
  integer :: iv,i,j
  logical :: yes_lat,yes_lon  
  integer :: ndims
  integer, dimension(3) :: idims
  type(edtype),target :: cgrid

  ! First check to see if their is lat/lon data in this dataset
  ! if the data exists, load it
  
  if(.not.no_ll) then
        
     !  Get the dimensioning information on latitude
     call shdf5_info_f('lat',ndims,idims)
     
     if(ndims /= 2) then
        print*,"NOT SET UP TO HAVE TIME VARYING LAT/LON"
        stop
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
        print*,"NOT SET UP TO HAVE TIME VARYING LAT/LON"
        stop
     endif
     
     !  Allocate the latitude array
     allocate(lon2d(idims(1),idims(2)))
     ndims = 2
     
     !  Read in the latitude array
     call shdf5_irec_f(ndims, idims, 'lon',  &
          rvara = lon2d )
     
     !  Determine the indices of the grid that each polygon sees
     !  returns poly%ilon and poly%ilat
     
     call match_poly_grid_array(cgrid,met_nlon(iformat),met_nlat(iformat),lon2d,lat2d)
     
     ! Deallocate the lat-lon arrays
     deallocate(lat2d,lon2d)

  endif
  
  return
  
end subroutine getll_array
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine calc_met_lapse_ar(cgrid,ipy)
  
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

  cpoly => cgrid%polygon(ipy)

  if (lapse_scheme == 0) then  !!!!! bypass
     do isi = 1,cpoly%nsites
        cpoly%met(isi)%geoht       = cgrid%met(ipy)%geoht
        cpoly%met(isi)%atm_tmp     = cgrid%met(ipy)%atm_tmp
        cpoly%met(isi)%atm_shv     = cgrid%met(ipy)%atm_shv
        cpoly%met(isi)%prss        = cgrid%met(ipy)%prss
        cpoly%met(isi)%pcpg        = cgrid%met(ipy)%pcpg
        cpoly%met(isi)%par_diffuse = cgrid%met(ipy)%par_diffuse
        cpoly%met(isi)%atm_co2     = cgrid%met(ipy)%atm_co2
        cpoly%met(isi)%rlong       = cgrid%met(ipy)%rlong
        cpoly%met(isi)%par_beam    = cgrid%met(ipy)%par_beam
        cpoly%met(isi)%nir_diffuse = cgrid%met(ipy)%nir_diffuse
        cpoly%met(isi)%nir_beam    = cgrid%met(ipy)%nir_beam
        cpoly%met(isi)%vels        = cgrid%met(ipy)%vels
      
!        if ( cpoly%met(isi)%rlong < rlong_min) then
!           print*,"Problems with RLONG A"  
!           print*,cpoly%met(isi)%rlong,cgrid%met(ipy)%rlong,rlong_min
!           call fatal_error('Problems with RLONG A','calc_met_lapse_ar','ed_met_driver.f90')
!        end if
!        if ( cpoly%met(isi)%atm_tmp < 150.0) then
!           print*,cpoly%met(isi)%atm_tmp,cgrid%met(ipy)%atm_tmp
!           call fatal_error('Problems with ATM TEMP A','calc_met_lapse_ar','ed_met_driver.f90')
!        endif
        if ( cpoly%met(isi)%atm_shv < 0.01e-2) then
            cpoly%met(isi)%atm_shv = 0.01e-2
!           print*,cpoly%met(isi)%atm_shv,cgrid%met(ipy)%atm_shv
!           call fatal_error('Problems with ATM MOISTURE A','calc_met_lapse_ar','ed_met_driver.f90')
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
        cpoly%met(isi)%pcpg    = cgrid%met(ipy)%pcpg    + cgrid%lapse(ipy)%pcpg*delE
        cpoly%met(isi)%atm_co2 = cgrid%met(ipy)%atm_co2 + cgrid%lapse(ipy)%atm_co2*delE
        cpoly%met(isi)%rlong   = cgrid%met(ipy)%rlong   + cgrid%lapse(ipy)%rlong*delE
        cpoly%met(isi)%par_diffuse = cgrid%met(ipy)%par_diffuse + cgrid%lapse(ipy)%par_diffuse*delE
        cpoly%met(isi)%par_beam    = cgrid%met(ipy)%par_beam    + cgrid%lapse(ipy)%par_beam*delE
        cpoly%met(isi)%nir_diffuse = cgrid%met(ipy)%nir_diffuse + cgrid%lapse(ipy)%nir_diffuse*delE
        cpoly%met(isi)%nir_beam    = cgrid%met(ipy)%nir_beam    + cgrid%lapse(ipy)%nir_beam*delE
        cpoly%met(isi)%vels    = cgrid%met(ipy)%vels    + cgrid%lapse(ipy)%vels*delE
        !! note: at this point VELS is vel^2.  Thus this lapse preserves mean wind ENERGY
        !! not wind SPEED
!        if ( cpoly%met(isi)%rlong < 200.0) then
!           print*,"Problems with RLONG A"  
!           print*,cpoly%met(isi)%rlong,cgrid%met(ipy)%rlong
!           call fatal_error('Problems with RLONG A','calc_met_lapse_ar','ed_met_driver.f90')
!        end if
!        if ( cpoly%met(isi)%atm_tmp < 150.0) then
!           print*,cpoly%met(isi)%atm_tmp,cgrid%met(ipy)%atm_tmp
!           call fatal_error('Problems with ATM TEMP A','calc_met_lapse_ar','ed_met_driver.f90')
!        endif
        if ( cpoly%met(isi)%atm_shv < 0.01e-2) then
            cpoly%met(isi)%atm_shv = 0.01e-2
!           print*,cpoly%met(isi)%atm_shv,cgrid%met(ipy)%atm_shv
!           call fatal_error('Problems with ATM MOISTURE A','calc_met_lapse_ar','ed_met_driver.f90')
        endif

        call MetDiagnostics(cpoly,ipy,isi)
        
     enddo
  endif
  return
end subroutine calc_met_lapse_ar
!==========================================================================================!
!==========================================================================================!


subroutine MetDiagnostics(cpoly,ipy,isi)
  use ed_state_vars,only    : edtype,polygontype
  use canopy_radiation_coms,only : rlong_min
  
  implicit none
  integer, intent(in) :: ipy
  integer, intent(in) :: isi
  type(polygontype),target :: cpoly
  
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
  if(cpoly%met(isi)%pcpg .lt. 0. .or. cpoly%met(isi)%pcpg .gt. 0.01)  then
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
     print*,cpoly%met(isi)%par_diffuse
     call fatal_error('Problems with DIFFUSE PAR','MetDiagnostics','ed_met_driver.f90')
  endif
  if(cpoly%met(isi)%par_beam .lt. 0. .or. cpoly%met(isi)%par_beam .gt. 1400.)  then
     print*,cpoly%met(isi)%par_beam
     call fatal_error('Problems with DIRECT BEAM PAR','MetDiagnostics','ed_met_driver.f90')
  endif
  if(cpoly%met(isi)%nir_diffuse .lt. 0. .or. cpoly%met(isi)%nir_diffuse .gt. 1400.)  then
     print*,cpoly%met(isi)%nir_diffuse
     call fatal_error('Problems with DIFUSE NIR','MetDiagnostics','ed_met_driver.f90')
  endif
  if(cpoly%met(isi)%nir_beam .lt. 0. .or. cpoly%met(isi)%nir_beam .gt. 1400.)  then
     print*,cpoly%met(isi)%nir_beam
     call fatal_error('Problems with DIRECT BEAM NIR','MetDiagnostics','ed_met_driver.f90')
  endif
  if(cpoly%met(isi)%vels .lt. 0.)  then
     print*,cpoly%met(isi)%vels
     call fatal_error('Problems with WIND VELOCITY','MetDiagnostics','ed_met_driver.f90')
  endif

  return
end subroutine MetDiagnostics




!==========================================================================================!
!==========================================================================================!
subroutine setLapseParms_ar(cgrid)
  
  use ed_state_vars,only:edtype
  use met_driver_coms, only: lapse

  implicit none

  
  integer :: init
  integer :: ipy
  type(edtype), target :: cgrid

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
  enddo

  return
end subroutine setLapseParms_ar
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine int_met_avg(cgrid)
  
  ! Increment the time averaged polygon met-forcing
  ! variables. THese will be normalized by the output
  ! period to give time averages of each quanity. The
  ! polygon level variables are derived from the
  ! weighted spatial average from the site level quantities.


  use ed_state_vars,only:edtype,polygontype,sitetype,patchtype
  use misc_coms,only : dtlsm,frqfast
  
  implicit none

  type(edtype), target :: cgrid
  type(polygontype),pointer :: cpoly
  type(sitetype),pointer :: csite
  type(patchtype),pointer :: cpatch

  integer :: ipy,isi,ipa,ico
  real :: frqfasti,tfact

  frqfasti = 1.0 / frqfast
  tfact = dtlsm * frqfasti

  do ipy = 1,cgrid%npolygons

     cpoly => cgrid%polygon(ipy)

     do isi = 1,cpoly%nsites
        
        cgrid%avg_nir_beam(ipy)    = cgrid%avg_nir_beam(ipy) + &
             cpoly%met(isi)%nir_beam * cpoly%area(isi) * tfact
        cgrid%avg_nir_diffuse(ipy) = cgrid%avg_nir_diffuse(ipy) + &
             cpoly%met(isi)%nir_diffuse * cpoly%area(isi) * tfact
        cgrid%avg_par_beam(ipy)    = cgrid%avg_par_beam(ipy) + &
             cpoly%met(isi)%par_beam * cpoly%area(isi) * tfact
        cgrid%avg_par_diffuse(ipy) = cgrid%avg_par_diffuse(ipy) + &
             cpoly%met(isi)%par_diffuse * cpoly%area(isi) * tfact
        cgrid%avg_atm_tmp(ipy)     = cgrid%avg_atm_tmp(ipy) + &
             cpoly%met(isi)%atm_tmp * cpoly%area(isi) * tfact
        cgrid%avg_atm_shv(ipy)     = cgrid%avg_atm_shv(ipy) + &
             cpoly%met(isi)%atm_shv * cpoly%area(isi) * tfact
        cgrid%avg_rhos(ipy)        = cgrid%avg_rhos(ipy) + &
             cpoly%met(isi)%rhos * cpoly%area(isi) * tfact
        cgrid%avg_rshort(ipy)      = cgrid%avg_rshort(ipy) + &
             cpoly%met(isi)%rshort * cpoly%area(isi) * tfact
        cgrid%avg_rshort_diffuse(ipy) = cgrid%avg_rshort_diffuse(ipy) + &
             cpoly%met(isi)%rshort_diffuse * cpoly%area(isi) * tfact
        cgrid%avg_rlong(ipy)       = cgrid%avg_rlong(ipy) + &
             cpoly%met(isi)%rlong * cpoly%area(isi) * tfact
        cgrid%avg_pcpg(ipy)        = cgrid%avg_pcpg(ipy) + &
             cpoly%met(isi)%pcpg * cpoly%area(isi) * tfact
        cgrid%avg_qpcpg(ipy)       = cgrid%avg_qpcpg(ipy) + &
             cpoly%met(isi)%qpcpg * cpoly%area(isi) * tfact
        cgrid%avg_dpcpg(ipy)       = cgrid%avg_dpcpg(ipy) + &
             cpoly%met(isi)%dpcpg * cpoly%area(isi) * tfact
        cgrid%avg_vels(ipy)        = cgrid%avg_vels(ipy) + &
             cpoly%met(isi)%vels * cpoly%area(isi) * tfact
        cgrid%avg_prss(ipy)        = cgrid%avg_prss(ipy) + &
             cpoly%met(isi)%prss * cpoly%area(isi) * tfact
        cgrid%avg_exner(ipy)       = cgrid%avg_exner(ipy) + &
             cpoly%met(isi)%exner * cpoly%area(isi) * tfact
        cgrid%avg_geoht(ipy)       = cgrid%avg_geoht(ipy) + &
             cpoly%met(isi)%geoht * cpoly%area(isi) * tfact
        cgrid%avg_atm_co2(ipy)     = cgrid%avg_atm_co2(ipy) + &
             cpoly%met(isi)%atm_co2 * cpoly%area(isi) * tfact
        cgrid%avg_albedt(ipy)      = cgrid%avg_albedt(ipy) + &
             0.5*(cpoly%albedo_beam(isi)+cpoly%albedo_diffuse(isi))*cpoly%area(isi) * tfact
        cgrid%avg_rlongup(ipy)     = cgrid%avg_rlongup(ipy) + &
             cpoly%rlongup(isi) * cpoly%area(isi) * tfact


     enddo
        
  enddo
  
  return
end subroutine int_met_avg
