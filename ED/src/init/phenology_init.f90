module phenology_startup
  
contains
    
  subroutine phenology_init
    
    use phenology_coms, only: iphen_scheme
    use ed_misc_coms, only: ied_init_mode

    implicit none

    logical,parameter :: bypass=.true.

    ! Initialize the Botta et al. scheme.
    select case (iphen_scheme)
    case (1)

       ! Initialize from satellite.  This subroutine gives ALL SITES the
       ! phenological parameters from Harvard Forest.
       
!!       call read_harvard_phenology
       call read_prescribed_phenology

    case default
       ! Only initialize thermal sums here if no thermal sums information is
       ! available from the restart file, or if this is a run with bare 
       ! ground initialization.
       
       if(ied_init_mode /= 4)then
          
          print'(/,a)','    Reading thermal sums'
          call read_thermal_sums
          
       end if

    end select


    return
  end subroutine phenology_init
    
  ! =============================================

  subroutine read_thermal_sums

    use ed_state_vars,only:edtype,polygontype,sitetype,edgrid_g
    use grid_coms,only:ngrids
    use ed_misc_coms, only: iyeara, imontha, idatea, itimea, ed_inputs_dir

    implicit none

    character(len=256) :: fname
    character(len=3), dimension(2), parameter :: ftype=(/'chd','dgd'/)
    logical :: exans
    integer :: npoints
    integer :: ierr
    real :: tlat
    real :: tlon
    real, dimension(2,12) :: var_current_year
    real, dimension(2,12) :: var_past_year
    real :: partial_month_fraction
    real, dimension(:), allocatable :: flat
    real, dimension(:), allocatable :: flon
    real, dimension(:,:,:), allocatable :: varc
    real, dimension(:,:,:), allocatable :: varp
    real, dimension(2,12) :: var_out
    integer :: i,j
    integer :: itype
    integer :: igr,ipy,isi
    logical, external :: isleap

    type(edtype),pointer :: cgrid
    type(polygontype),pointer :: cpoly
    type(sitetype),pointer :: csite


    ! How far are we into the current month?  Use this to weight the 
    ! current and past year data for this month.
    select case (imontha)
    case(1,3,5,7,8,10,12)
       partial_month_fraction = (real(idatea-1) + real(itimea) * 0.01 /   &
            24.0) / 31.0

    case(4,6,9,11)
       partial_month_fraction = (real(idatea-1) + real(itimea) * 0.01 /   &
            24.0) / 30.0

    case(2)
       if (isleap(iyeara)) then
          partial_month_fraction = (real(idatea-1) + real(itimea) * 0.01 /   &
               24.0) / 29.0
      else
          partial_month_fraction = (real(idatea-1) + real(itimea) * 0.01 /   &
               24.0) / 28.0
       end if

    end select
    
    ! The two types are chilling days and degree days.
    do itype = 1, 2
       
       ! See if there is a file for this year
       write(fname,'(5a,i4.4,a)')trim(ed_inputs_dir),ftype(itype),'/temp.',  &
            ftype(itype),'.y',iyeara, '.dat'
       inquire(file=trim(fname),exist=exans)
       
       ! If not, use average data
       if(.not.exans)write(fname,'(2a)') trim(ed_inputs_dir)//ftype(itype),  &
            '/temp.'//ftype(itype)//'.avg.dat'
       inquire(file=trim(fname),exist=exans)
       if(.not.exans)then
          print*,'File ',trim(fname),' does not exist.'
          print*,'Stopping in subroutine read_thermal_sums'
          stop
       endif
       
       ! Open file
       open(12,file=trim(fname),form='formatted',status='old')
       
       if(itype == 1)then
          ! Count number of points
          npoints = 0
          count_points: do
             read(12,*,iostat=ierr)tlat,tlon,(var_current_year(itype,j),j=1,12)
             if(ierr == 0)then
                npoints = npoints + 1
             else
                exit count_points
             endif
          enddo count_points
          
          ! Allocate arrays
          allocate(flat(npoints))
          allocate(flon(npoints))
          allocate(varc(2,npoints,12))
          allocate(varp(2,npoints,12))
          varc(:,:,:) = 0.0
          varp(:,:,:) = 0.0
          
          rewind(12)
       endif
       
       ! Read in file
       do i = 1, npoints
          read(12,*)flat(i),flon(i),(varc(itype,i,j),j=1,12)
       enddo
       close(12)
       
       ! Now do the previous year.  First, get file name.
       write(fname,'(2a,i4.4,a)')trim(ed_inputs_dir)//ftype(itype),  &
            '/temp.'//ftype(itype)//'.y',iyeara-1, '.dat'
       inquire(file=trim(fname),exist=exans)
       if(.not.exans)write(fname,'(2a)') trim(ed_inputs_dir)//ftype(itype),  &
            '/temp.'//ftype(itype)//'.avg.dat'
       inquire(file=trim(fname),exist=exans)
       if(.not.exans)then
          print*,'File ',trim(fname),' does not exist.'
          print*,'Stopping in ed_database.f90'
          stop
       endif
       
       ! Open file.
       open(12,file=trim(fname),form='formatted',status='old')
       do i = 1, npoints
          read(12,*)tlat,tlon,(varp(itype,i,j),j=1,12)
       enddo
       close(12)
       
    enddo
    
    do igr=1,ngrids

       cgrid => edgrid_g(igr)

       do ipy=1,cgrid%npolygons

          cpoly => cgrid%polygon(ipy)

          do isi = 1,cpoly%nsites
             
             csite => cpoly%site(isi)
             
             ! Find the right index
             exans = .false.
             var_out(:,:) = 0.0
             find_index: do i=1,npoints
                if(abs(flat(i)-cgrid%lat(ipy)) < 1.0 .and.   &
                     abs(flon(i)-cgrid%lon(ipy)) < 1.0)then
                   exans = .true.
                   exit find_index
                endif
             enddo find_index
             
             ! Fill this site's info
             if(exans)then
                do itype = 1,2
                   var_current_year(itype,1:12) = varc(itype,i,1:12)
                   var_past_year(itype,1:12) = varp(itype,i,1:12)
                   ! Fill contribution from current month
                   var_out(itype,imontha) = partial_month_fraction *   &
                        var_current_year(itype,imontha)
                   ! Fill contribution from previous months in this year
                   var_out(itype,1:(imontha-1)) =   &
                        var_current_year(itype,1:(imontha-1))
                   ! Fill contribution from previous year
                   var_out(itype,imontha) = (1.0 - partial_month_fraction) *   &
                        var_past_year(itype,imontha)
                   ! Fill remaining info from past year
                   var_out(itype,(imontha+1):12) =   &
                        var_past_year(itype,(imontha+1):12)
                enddo
             endif
             
             ! Fill the patch-level degree and chill days
             call fill_thermal_sums(csite, cgrid%lat(ipy), imontha, var_out)
             
          enddo
       enddo
    enddo
       
    deallocate(flat)
    deallocate(flon)
    deallocate(varc)
    deallocate(varp)
    
    return
  end subroutine read_thermal_sums

  !==================================================================
  
  subroutine fill_thermal_sums(csite, lat, imontha, therm_sums)
    
    use ed_state_vars,only : sitetype

    implicit none
    
    type(sitetype),target :: csite
    real, intent(in) :: lat
    integer, intent(in) :: imontha
    real, dimension(2,12), intent(in) :: therm_sums
    real :: dgd
    real :: chd
    integer :: ipa

    if(lat >= 0.0)then
       
       if(imontha <= 8)then
          dgd = sum(therm_sums(2,1:imontha))
       else
          dgd = 0.0
       endif
       
       if(imontha >= 11)then
          chd = sum(therm_sums(1,imontha:12))
       elseif(imontha <= 6)then
          chd = sum(therm_sums(1,11:12)) + sum(therm_sums(1,1:imontha))
       else
          chd = 0.0
       endif
       
    else
       
       if(imontha <= 2)then
          dgd = sum(therm_sums(2,7:12)) + sum(therm_sums(2,1:imontha))
       elseif(imontha >= 7)then
          dgd = sum(therm_sums(2,7:imontha))
       else
          dgd = 0.0
       endif
       
       if(imontha >= 5)then
          chd = sum(therm_sums(1,5:imontha))
       else
          chd = 0.0
       endif
       
    endif
    
    ! Loop over patches
    
    do ipa = 1,csite%npatches
       csite%sum_chd(ipa) = chd
       csite%sum_dgd(ipa) = dgd
    enddo
    
    return
  end subroutine fill_thermal_sums

    ! ==============================================

  subroutine read_harvard_phenology

    use ed_state_vars,only:edgrid_g,edtype,polygontype,sitetype
    use ed_misc_coms, only: imontha, idatea, iyeara
    use grid_coms,only:ngrids

    implicit none

    type(edtype),pointer :: cgrid
    type(polygontype),pointer :: cpoly
    integer :: igr,isi,ipy
    integer :: doy
    integer, external :: julday

    ! Give ALL SITES the phenology parameters of Harvard Forest.

    do igr = 1,ngrids
       
       cgrid => edgrid_g(igr)
       
       do ipy = 1,cgrid%npolygons
          
          cpoly => cgrid%polygon(ipy)
          
          do isi = 1,cpoly%nsites
             
             ! Allocate memory for all years having data.

             ! THIS IS SERIOUSLY HARD-CODED, AND SHOULD BE AVOIDED BY ALL MEANS AT
             ! THIS STAGE!
             
             cpoly%phen_pars(isi)%nyears = 12
             allocate(cpoly%phen_pars(isi)%years(cpoly%phen_pars(isi)%nyears))
             allocate(cpoly%phen_pars(isi)%flush_a(cpoly%phen_pars(isi)%nyears))
             allocate(cpoly%phen_pars(isi)%flush_b(cpoly%phen_pars(isi)%nyears))
             allocate(cpoly%phen_pars(isi)%color_a(cpoly%phen_pars(isi)%nyears))
             allocate(cpoly%phen_pars(isi)%color_b(cpoly%phen_pars(isi)%nyears))
                   
             ! 1992
             cpoly%phen_pars(isi)%years(1) = 1992
             cpoly%phen_pars(isi)%flush_a(1) = 0.00621
             cpoly%phen_pars(isi)%flush_b(1) = -15.1
             cpoly%phen_pars(isi)%color_a(1) = 0.00360
             cpoly%phen_pars(isi)%color_b(1) = 50.1
             
             ! 1993
             cpoly%phen_pars(isi)%years(2) = 1993
             cpoly%phen_pars(isi)%flush_a(2) = 0.00709
             cpoly%phen_pars(isi)%flush_b(2) = -19.1
             cpoly%phen_pars(isi)%color_a(2) = 0.00359
             cpoly%phen_pars(isi)%color_b(2) = 48.4
             
             ! 1994
             cpoly%phen_pars(isi)%years(3) = 1994
             cpoly%phen_pars(isi)%flush_a(3) = 0.00688
             cpoly%phen_pars(isi)%flush_b(3) = -19.5
             cpoly%phen_pars(isi)%color_a(3) = 0.00363
             cpoly%phen_pars(isi)%color_b(3) = 37.1
             
             ! 1995
             cpoly%phen_pars(isi)%years(4) = 1995
             cpoly%phen_pars(isi)%flush_a(4) = 0.00680
             cpoly%phen_pars(isi)%flush_b(4) = -24.2
             cpoly%phen_pars(isi)%color_a(4) = 0.00365
             cpoly%phen_pars(isi)%color_b(4) = 31.4
             
             ! 1996
             cpoly%phen_pars(isi)%years(5) = 1996
             cpoly%phen_pars(isi)%flush_a(5) = 0.00673
             cpoly%phen_pars(isi)%flush_b(5) = -20.0
             cpoly%phen_pars(isi)%color_a(5) = 0.00357
             cpoly%phen_pars(isi)%color_b(5) = 47.3
             
             ! 1997
             cpoly%phen_pars(isi)%years(6) = 1997
             cpoly%phen_pars(isi)%flush_a(6) = 0.00653
             cpoly%phen_pars(isi)%flush_b(6) = -20.9
             cpoly%phen_pars(isi)%color_a(6) = 0.00358
             cpoly%phen_pars(isi)%color_b(6) = 46.9
             
             ! 1998
             cpoly%phen_pars(isi)%years(7) = 1998
             cpoly%phen_pars(isi)%flush_a(7) = 0.00726
             cpoly%phen_pars(isi)%flush_b(7) = -14.8
             cpoly%phen_pars(isi)%color_a(7) = 0.00363
             cpoly%phen_pars(isi)%color_b(7) = 34.8
             
             ! 1999
             cpoly%phen_pars(isi)%years(8) = 1999
             cpoly%phen_pars(isi)%flush_a(8) = 0.00702
             cpoly%phen_pars(isi)%flush_b(8) = -17.3
             cpoly%phen_pars(isi)%color_a(8) = 0.00355
             cpoly%phen_pars(isi)%color_b(8) = 45.1
             
             ! 2000
             cpoly%phen_pars(isi)%years(9) = 2000
             cpoly%phen_pars(isi)%flush_a(9) = 0.00683
             cpoly%phen_pars(isi)%flush_b(9) = -14.8
             cpoly%phen_pars(isi)%color_a(9) = 0.00360
             cpoly%phen_pars(isi)%color_b(9) = 40.5
             
             ! 2001
             cpoly%phen_pars(isi)%years(10) = 2001
             cpoly%phen_pars(isi)%flush_a(10) = 0.00705
             cpoly%phen_pars(isi)%flush_b(10) = -13.7
             cpoly%phen_pars(isi)%color_a(10) = 0.00358
             cpoly%phen_pars(isi)%color_b(10) = 36.1
             
             ! 2002
             cpoly%phen_pars(isi)%years(11) = 2002
             cpoly%phen_pars(isi)%flush_a(11) = 0.00682
             cpoly%phen_pars(isi)%flush_b(11) = -19.6
             cpoly%phen_pars(isi)%color_a(11) = 0.00352
             cpoly%phen_pars(isi)%color_b(11) = 36.3
             
             ! 2003
             cpoly%phen_pars(isi)%years(12) = 2003
             cpoly%phen_pars(isi)%flush_a(12) = 0.0066
             cpoly%phen_pars(isi)%flush_b(12) = -23.0
             cpoly%phen_pars(isi)%color_a(12) = 0.00356
             cpoly%phen_pars(isi)%color_b(12) = 48.3
             
             ! Initialize green_leaf_factor and leaf_aging_factor.
             doy = julday(imontha,idatea,iyeara)
             call prescribed_leaf_state(cgrid%lat(ipy), imontha, iyeara, doy,   &
                  cpoly%green_leaf_factor(:,isi), cpoly%leaf_aging_factor(:,ipy), &
                  cpoly%phen_pars(isi))
             
          enddo
       enddo
    enddo
    
    return
  end subroutine read_harvard_phenology
  
    ! ==============================================

  subroutine read_prescribed_phenology

    use ed_state_vars,only:edgrid_g,edtype,polygontype,sitetype
    use ed_misc_coms, only: imontha, idatea, iyeara
    use grid_coms,only:ngrids
    use phenology_coms, only: prescribed_phen,phenpath
    use ed_max_dims, only: str_len

    implicit none

    type(edtype),pointer :: cgrid
    type(polygontype),pointer :: cpoly
    type(prescribed_phen) :: phen_temp
    integer :: igr,isi,ipy,iyr
    integer :: doy
    integer, external :: julday
    character(len=str_len) :: phen_name,clat,clon
    logical :: phenology_exist
    
    do igr = 1,ngrids
       
       cgrid => edgrid_g(igr)
       
       !!Open up phenology folder

       do ipy = 1,cgrid%npolygons
          
          cpoly => cgrid%polygon(ipy)
          
          !! build phenology file name
          !!lat
          if(cgrid%lat(ipy) <= -10.0)then
             write(clat,'(f5.1)')cgrid%lat(ipy)
          elseif(cgrid%lat(ipy) < 0.0 .or. cgrid%lat(ipy) >= 10.0) then
             write(clat,'(f4.1)')cgrid%lat(ipy)
          else
             write(clat,'(f3.1)')cgrid%lat(ipy)
          endif
          !!lon
          if(cgrid%lon(ipy) <= -100.0)then
             write(clon,'(f6.1)')cgrid%lon(ipy)
          elseif(cgrid%lon(ipy) <= -10.0 .or. cgrid%lon(ipy) >= 100.0) then
             write(clon,'(f5.1)')cgrid%lon(ipy)
          elseif(cgrid%lon(ipy) < 0.0) then
             write(clon,'(f4.1)')cgrid%lon(ipy)
          elseif(cgrid%lon(ipy) < 10.) then
             write(clon,'(f3.1)')cgrid%lon(ipy)
          else
             write(clon,'(f4.1)')cgrid%lon(ipy)
          endif          

          !!open phenology file
          write(phen_name,'(6a)')trim(phenpath),'.lat',trim(clat),'lon',trim(clon),'.txt'
          inquire(file=trim(phen_name),exist=phenology_exist)
          if(.not.phenology_exist) then
             print*,"phenology file not found :: ",phen_name
             stop
          endif

          !! read phenology file
          open(12,file=trim(phen_name),form='formatted',status='old')
          read(12,*) phen_temp%nyears
          allocate(phen_temp%years(phen_temp%nyears))
          allocate(phen_temp%flush_a(phen_temp%nyears))
          allocate(phen_temp%flush_b(phen_temp%nyears))
          allocate(phen_temp%color_a(phen_temp%nyears))
          allocate(phen_temp%color_b(phen_temp%nyears))

          do iyr = 1,phen_temp%nyears
             read(12,*) phen_temp%years(iyr),phen_temp%flush_a(iyr),phen_temp%flush_b(iyr), &
                  phen_temp%color_a(iyr),phen_temp%color_b(iyr)
          enddo

          !! write phenology to each site
          do isi = 1,cpoly%nsites
             
             ! Allocate memory for all years having data.
             cpoly%phen_pars(isi)%nyears = phen_temp%nyears
             allocate(cpoly%phen_pars(isi)%years(cpoly%phen_pars(isi)%nyears))
             allocate(cpoly%phen_pars(isi)%flush_a(cpoly%phen_pars(isi)%nyears))
             allocate(cpoly%phen_pars(isi)%flush_b(cpoly%phen_pars(isi)%nyears))
             allocate(cpoly%phen_pars(isi)%color_a(cpoly%phen_pars(isi)%nyears))
             allocate(cpoly%phen_pars(isi)%color_b(cpoly%phen_pars(isi)%nyears))
                   
             do iyr = 1,phen_temp%nyears
                cpoly%phen_pars(isi)%years(iyr)   = phen_temp%years(iyr)
                cpoly%phen_pars(isi)%flush_a(iyr) = phen_temp%flush_a(iyr)
                cpoly%phen_pars(isi)%flush_b(iyr) = phen_temp%flush_b(iyr)
                cpoly%phen_pars(isi)%color_a(iyr) = phen_temp%color_a(iyr)
                cpoly%phen_pars(isi)%color_b(iyr) = phen_temp%color_b(iyr)
             enddo
             
             ! Initialize green_leaf_factor and leaf_aging_factor.
             doy = julday(imontha,idatea,iyeara)
             call prescribed_leaf_state(cgrid%lat(ipy), imontha, iyeara, doy,   &
                  cpoly%green_leaf_factor(:,isi), cpoly%leaf_aging_factor(:,ipy), &
                  cpoly%phen_pars(isi))
             
          enddo
       enddo
    enddo
    
    deallocate(phen_temp%years,phen_temp%flush_a,phen_temp%flush_b,phen_temp%color_a,phen_temp%color_b)
    
    return
  end subroutine read_prescribed_phenology


end module phenology_startup
