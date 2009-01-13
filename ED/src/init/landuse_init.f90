!===================================================================
subroutine landuse_init_array

  use ed_state_vars, only: edtype,polygontype,sitetype,edgrid_g
  use consts_coms, only: erad, pio180
  use disturb_coms, only: lutime, max_lu_years,num_lu_trans, ianth_disturb
  use misc_coms, only: iyeara
  use grid_coms,only: ngrids

  implicit none

  type(edtype),pointer      :: cgrid
  type(polygontype),pointer :: cpoly
  type(sitetype),pointer    :: csite

  real :: file_lat
  real :: file_lon
  character(len=256) :: fname
  real :: lu_area
  logical :: exans
  integer :: iyear,igr,ipy,isi

  do igr = 1,ngrids

     cgrid=>edgrid_g(igr)

     do ipy = 1,cgrid%npolygons

        cpoly => cgrid%polygon(ipy)

        ! Generate the landuse file name
        call landuse_file_name(cgrid%lat(ipy), cgrid%lon(ipy), file_lat,   &
             file_lon, fname)
        
        ! Use file_lat to compute the physical area sampled by the file
        lu_area = (erad * pio180)**2 * abs(cos(pio180 * file_lat))
        
        ! If we are doing anthropogenic disturbance, check for file existence.
        inquire(file=trim(fname),exist=exans)
        
        if(exans .and. ianth_disturb==1)then
           
           allocate(cpoly%clutimes(max_lu_years,cpoly%nsites))
           
           do isi = 1,cpoly%nsites
              
              csite => cpoly%site(isi)
              
              ! Land use file exists
              open(12,file=trim(fname),form='formatted',status='old')
              
              ! Skip header
              read(12,*)
              
              ! Each GLU file has 300 years, 1700--1999.
              cpoly%num_landuse_years(isi) = max_lu_years
              
              ! Loop over years
              do iyear = 1,cpoly%num_landuse_years(isi)
                 
                 read(12,*) &
                      cpoly%clutimes(iyear,isi)%landuse_year, &
                      cpoly%clutimes(iyear,isi)%landuse(1:num_lu_trans)
                 
                 ! Normalize by the area
                 cpoly%clutimes(iyear,isi)%landuse(12) = &
                      cpoly%clutimes(iyear,isi)%landuse(12) / lu_area
                 cpoly%clutimes(iyear,isi)%landuse(14) = &
                      cpoly%clutimes(iyear,isi)%landuse(14) / lu_area
                 cpoly%clutimes(iyear,isi)%landuse(16) = &
                      cpoly%clutimes(iyear,isi)%landuse(16) / lu_area
                 cpoly%clutimes(iyear,isi)%landuse(18) = &
                      cpoly%clutimes(iyear,isi)%landuse(18) / lu_area
                 
              enddo
              close(12)

           enddo
           
        else
           
           ! no GLU data for this site.  probably water, or anthropogenic
           ! disturbance is turned off.
           
           if (ianth_disturb==1) print*,'file: ',trim(fname),' not found.  assigning 0s for landuse'
           
           ! Allocate 1 landuse year
           allocate(cpoly%clutimes(1,cpoly%nsites))

           do isi = 1,cpoly%nsites
              cpoly%num_landuse_years(isi) = 1
              cpoly%clutimes(1,isi)%landuse_year = iyeara
              cpoly%clutimes(1,isi)%landuse(1:num_lu_trans) = 0.0
           enddo

        endif
 
        cpoly%plantation(:) = 0
        call read_plantation_fractions_array(cpoly, file_lat, file_lon)
        
     enddo
  enddo

  return
end subroutine landuse_init_array
   
! ==============================================================

subroutine landuse_file_name(lat, lon, file_lat, file_lon, fname)
  
  use misc_coms, only: ed_inputs_dir

  implicit none

  real, intent(in) :: lat
  real, intent(in) :: lon
  real, intent(out) :: file_lat
  character(len=256), intent(out) :: fname
  character(len=5) :: file_lat_string
  real, intent(out) :: file_lon
  character(len=6) :: file_lon_string
  
  if(lat > 0.0)then
     file_lat = 0.5 + real(int(lat))
     if(file_lat < 10.0)then
        write(file_lat_string,'(f3.1)')file_lat
     else
        write(file_lat_string,'(f4.1)')file_lat
     endif
  else
     file_lat = - (0.5 + real(int(-lat)))
     if(file_lat > -10.0)then
        write(file_lat_string,'(f4.1)')file_lat
     else
        write(file_lat_string,'(f5.1)')file_lat
     endif
  endif

  if(lon > 0.0)then
     file_lon = 0.5 + real(int(lon))
     if(file_lon < 10.0)then
        write(file_lon_string,'(f3.1)')file_lon
     elseif(lon < 100.0)then
        write(file_lon_string,'(f4.1)')file_lon
     else
        write(file_lon_string,'(f5.1)')file_lon
     endif
  else
     file_lon = - (0.5 + real(int(-lon)))
     if(lon > -10.0)then
        write(file_lon_string,'(f4.1)')file_lon
     elseif(lon > -100.0)then
        write(file_lon_string,'(f5.1)')file_lon
     else
        write(file_lon_string,'(f6.1)')file_lon
     endif
  endif
  
  fname = trim(ed_inputs_dir)//'glu/lat'//trim(file_lat_string)
  fname = trim(fname)//'lon'//trim(file_lon_string)//'.lu'

  return
end subroutine landuse_file_name

!===========================================================================

subroutine read_plantation_fractions_array(cpoly, file_lat, file_lon)

  use ed_state_vars,only:polygontype
  use misc_coms, only: ed_inputs_dir

  implicit none

  type(polygontype),target :: cpoly
  real, intent(in) :: file_lat
  real, intent(in) :: file_lon
  character(len=256) :: fname
  logical :: exans
  integer :: irec,isi
  real :: lat
  real :: lon
  real :: fracplant

  ! Set plantation fraction
  fname = trim(ed_inputs_dir)//'fraction.plantation'
  inquire(file=trim(fname),exist=exans)
  if(.not.exans)then
     print*,'There is no plantation file.  Exiting.'
     print*,'Assigning no plantations'
     print*,trim(fname)
     return
!     stop
  endif

  open(12, file=trim(fname), form='formatted', status='old')
  read_plantation: do
     read(12,*,iostat=irec)lat, lon, fracplant
     if(irec /= 0)exit read_plantation
     if(lat == file_lat .and. lon == file_lon)then
        if(fracplant > 0.125) then
           do isi = 1,cpoly%nsites
              cpoly%plantation(isi) = 1
           enddo
        endif
        exit read_plantation
     endif
  enddo read_plantation
  close(12)

  return
end subroutine read_plantation_fractions_array
