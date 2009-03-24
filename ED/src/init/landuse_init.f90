!==========================================================================================!
!==========================================================================================!
!     This file will read the landuse files and assign the anthropogenic disturbance       !
! matrices.                                                                                !
!  MLO.  Added the possibility to start the run before the initial anthropogenic           !
!        disturbance, by padding zeroes in the years before and after the dataset.         !
!------------------------------------------------------------------------------------------!
subroutine landuse_init_array

   use ed_state_vars , only : edtype         & ! structure
                            , polygontype    & ! structure
                            , sitetype       & ! structure
                            , edgrid_g       ! ! structure
   use consts_coms   , only : erad           & ! intent(in)
                            , pio180         ! ! intent(in)
   use disturb_coms  , only : lutime         & ! intent(in)
                            , max_lu_years   & ! intent(in)
                            , num_lu_trans   & ! intent(in)
                            , ianth_disturb  ! ! intent(in)
   use misc_coms     , only : iyeara         & ! intent(in)
                            , iyearz         ! ! intent(in)
   use grid_coms     , only : ngrids         ! ! intent(in)
   use max_dims      , only : str_len        ! ! intent(in)
   implicit none
   !----- Local variables -----------------------------------------------------------------!
   type(edtype)     , pointer    :: cgrid
   type(polygontype), pointer    :: cpoly
   type(sitetype)   , pointer    :: csite
   real                          :: file_lat
   real                          :: file_lon
   character(len=256)            :: fname
   real                          :: lu_area_i
   integer                       :: ierr
   logical                       :: exans
   integer                       :: iyear
   integer                       :: igr
   integer                       :: ipy
   integer                       :: isi
   integer                       :: sim_years
   integer                       :: lu_years
   integer                       :: yd_1st
   integer                       :: yd_this
   integer                       :: yd_last
   integer                       :: yd_tot
   !---------------------------------------------------------------------------------------!

   !----- Finding number of simulation years ----------------------------------------------!
   sim_years = iyearz-iyeara+1

   !----- Crashing the run if the user set up a very long run... --------------------------!
   if (ianth_disturb == 1 .and. sim_years > max_lu_years) then
      write (unit=*,fmt='(a,1x,i5)') 'IYEARA       (From namelist)        :',iyeara
      write (unit=*,fmt='(a,1x,i5)') 'IYEARZ       (From namelist)        : ',iyearz
      write (unit=*,fmt='(a,1x,i5)') 'MAX_LU_YEARS (From disturb_coms.f90): ',max_lu_years
      write (unit=*,fmt='(a)') ' Your run is too long.  Try increasing max_lu_years,'
      write (unit=*,fmt='(a)') ' so MAX_LU_YEARS >= IYEARZ-IYEARA+1, then recompile'
      write (unit=*,fmt='(a)') ' your model.'
      call fatal_error ('Simulation is too long for anthropogenic disturbance.'            &
                      &,'landuse_init_array','landuse_init.f90')
   end if

   do igr = 1,ngrids
      cgrid=>edgrid_g(igr)

      do ipy = 1,cgrid%npolygons
         cpoly => cgrid%polygon(ipy)

         !----- Generate the landuse file name. -------------------------------------------!
         call landuse_file_name(cgrid%lat(ipy),cgrid%lon(ipy),file_lat,file_lon,fname)

         !----- Use file_lat to compute the physical area sampled by the file. ------------!
         lu_area_i = 1. / ((erad * pio180)*(erad * pio180) * abs(cos(pio180 * file_lat)))

         !----- If we are doing anthropogenic disturbance, check for file existence. ------!
         inquire(file=trim(fname),exist=exans)
         
         if (exans .and. ianth_disturb == 1) then
            !----- File exists, allocate the maximum number of years ----------------------!
            allocate(cpoly%clutimes(max_lu_years,cpoly%nsites))
            
            siteloop: do isi = 1,cpoly%nsites
               csite => cpoly%site(isi)
               
               !----- Open the file -------------------------------------------------------!
               open(unit=12,file=trim(fname),form='formatted',status='old')
               
               !---------------------------------------------------------------------------!
               !    Here we pretend we read the file, but only read the first and last     !
               ! year.  We will zero-pad the years before and after the land use data set  !
               ! if needed.                                                                !
               !---------------------------------------------------------------------------!
               read (unit=12,fmt=*) ! Header
               read(unit=12,fmt=*,iostat=ierr) yd_1st
               if (ierr /= 0) then
                  call fatal_error('Error reading landuse file :'//trim(fname)//'...'      &
                                  &,'landuse_init_array','landuse_init.f90')
               end if
               yd_last = yd_1st
               !----- Read years until the end to find out first and last year. -----------!
               rangeread: do 
                  read (unit=12,fmt=*,iostat=ierr) yd_this
                  if (ierr == 0) then
                     !----- Found another year. -------------------------------------------!
                     yd_last = yd_this
                  else
                     !----- End of file. --------------------------------------------------!
                     exit rangeread
                  end if
               end do rangeread
               !----- Rewind the file and skip the header. --------------------------------!
               rewind(unit=12)
               read(unit=12,fmt=*)
               !---------------------------------------------------------------------------!


               !----- Determine the number of disturbance years. --------------------------!
               cpoly%num_landuse_years(isi) = max(yd_last,iyearz)-min(yd_1st,iyeara) + 1


               !----- Padding disturbances with zero before first available lu year. ------!
               iyear = 0
               do yd_this = iyeara,(yd_1st-1)
                  iyear = iyear + 1
                  cpoly%clutimes(iyear,isi)%landuse_year            = yd_this
                  cpoly%clutimes(iyear,isi)%landuse(1:num_lu_trans) = 0.0
               end do

               !---- Reading the years that have data -------------------------------------!
               do yd_this = yd_1st,yd_last
                  iyear = iyear + 1
                  read(unit=12,fmt=*) cpoly%clutimes(iyear,isi)%landuse_year               &
                                     ,cpoly%clutimes(iyear,isi)%landuse(1:num_lu_trans)
                  
                  !----- Normalize by the area --------------------------------------------!
                  cpoly%clutimes(iyear,isi)%landuse(12) = lu_area_i                        &
                                           * cpoly%clutimes(iyear,isi)%landuse(12) 

                  cpoly%clutimes(iyear,isi)%landuse(14) = lu_area_i                        &
                                           * cpoly%clutimes(iyear,isi)%landuse(14) 
                                           
                  cpoly%clutimes(iyear,isi)%landuse(16) = lu_area_i                        &
                                           * cpoly%clutimes(iyear,isi)%landuse(16)

                  cpoly%clutimes(iyear,isi)%landuse(18) = lu_area_i                        &
                                           * cpoly%clutimes(iyear,isi)%landuse(18)
                  
               end do
               close(unit=12,status='keep')

               !----- Padding disturbances with zero after last available lu year. --------!
               do yd_this = (yd_last+1),iyearz
                  iyear = iyear + 1
                  cpoly%clutimes(iyear,isi)%landuse_year            = yd_this
                  cpoly%clutimes(iyear,isi)%landuse(1:num_lu_trans) = 0.0
               end do

            end do siteloop
         else
            !------------------------------------------------------------------------------!
            !      No GLU data for this site.  Probably water, or anthropogenic            !
            ! disturbance is turned off.                                                   !
            !------------------------------------------------------------------------------!
            
            if (ianth_disturb==1) then
               write(unit=*,fmt='(3(a,1x))') 'Land use file:',trim(fname)                  &
                                           &,' not found.  assigning 0s for landuse'
            end if

            !----- Allocate just 1 landuse year. ------------------------------------------!
            allocate(cpoly%clutimes(1,cpoly%nsites))
            do isi = 1,cpoly%nsites
               cpoly%num_landuse_years(isi) = 1
               cpoly%clutimes(1,isi)%landuse_year = iyeara
               cpoly%clutimes(1,isi)%landuse(1:num_lu_trans) = 0.0
            end do

         end if
 
         cpoly%plantation(:) = 0
         call read_plantation_fractions_array(cpoly, file_lat, file_lon)
         
      enddo
   enddo

   return
end subroutine landuse_init_array
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
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
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine reads the plantation fraction.                                       !
!------------------------------------------------------------------------------------------!
subroutine read_plantation_fractions_array(cpoly, file_lat, file_lon)
   use ed_state_vars , only : polygontype   ! ! structure
   use misc_coms     , only : ed_inputs_dir ! ! intent(in)
   use max_dims      , only : str_len       ! ! intent(in)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(polygontype), target     :: cpoly
   real             , intent(in) :: file_lat
   real             , intent(in) :: file_lon
   !----- Local variables -----------------------------------------------------------------!
   character(len=str_len)        :: fname
   logical                       :: exans
   integer                       :: ierr
   integer                       :: isi
   real                          :: lat
   real                          :: lon
   real                          :: fracplant
   !---------------------------------------------------------------------------------------!


   !----- Set plantation fraction file and checking whether it exists. --------------------!
   fname = trim(ed_inputs_dir)//'fraction.plantation'
   inquire(file=trim(fname),exist=exans)
   if (.not.exans)then
      write(unit=*,fmt='(a)') 'There is no plantation file.  Exiting.'
      write(unit=*,fmt='(a)') 'Assigning no plantations'
      write(unit=*,fmt='(a)') 'File :'//trim(fname)//' not found...'
      return
   end if

   open(unit=12, file=trim(fname), form='formatted', status='old')
   read_plantation: do
      read (unit=12,fmt=*,iostat=ierr) lat, lon, fracplant
      !----- End of file. -----------------------------------------------------------------!
      if(ierr /= 0) exit read_plantation

      if (lat == file_lat .and. lon == file_lon) then
         if(fracplant > 0.125) then
            do isi = 1,cpoly%nsites
               cpoly%plantation(isi) = 1
            end do
         end if
         exit read_plantation
      end if
   end do read_plantation
   close(unit=12)

   return
end subroutine read_plantation_fractions_array
!==========================================================================================!
!==========================================================================================!
