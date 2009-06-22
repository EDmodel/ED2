!==========================================================================================!
!==========================================================================================!
!    This subroutine will generate the meteorological driver header at the same directory  !
! as the HDF5 files.                                                                       !
!------------------------------------------------------------------------------------------!
subroutine ed_metd_header()
   use mod_maxdims , only : maxstr     ! ! intent(in)
   use mod_ncep    , only : nvars_ol1  & ! intent(in)
                          , nvars_ol2  & ! intent(in)
                          , vars_ol1   & ! intent(in)
                          , vars_ol2   & ! intent(in)
                          , vtype_ol1  & ! intent(in)
                          , vtype_ol2  ! ! intent(in)
   use mod_grid    , only : grid_g     & ! intent(in)
                          , ssxp       & ! intent(in)
                          , ssyp       & ! intent(in)
                          , sstp       ! ! intent(in)
   use mod_ioopts  , only : intype     & ! intent(in)
                          , inpfrq     & ! intent(in)
                          , outpref    & ! intent(in)
                          , edgeoff    ! ! intent(in)
   implicit none
   !----- Local variables. ----------------------------------------------------------------!
   character(len=maxstr)                  :: headname
   character(len=maxstr)                  :: outname
   character(len=maxstr)                  :: varistring
   character(len=maxstr)                  :: timestring
   character(len=maxstr)                  :: typestring
   integer                                :: mxo
   integer                                :: myo
   integer                                :: xa
   integer                                :: xz
   integer                                :: ya
   integer                                :: yz
   integer                                :: nv
   real                                   :: dlon
   real                                   :: dlat
   real                                   :: lon0
   real                                   :: lat0
   !---------------------------------------------------------------------------------------!

   !----- Creating the header file name and open it. --------------------------------------!
   write(headname,fmt='(2a)') trim(outpref),'_HEADER'
   open(unit=59,file=trim(headname),status='replace',form='formatted',access='sequential')

   !----- Write the first line. -----------------------------------------------------------!
   write(unit=59,fmt='(a)') '# See README at the bottom of this file.'

   !----- Now we write the intype-dependent part. -----------------------------------------!
   select case (trim(intype))
   case ('ncep')

      !------------------------------------------------------------------------------------!
      !     Fill some variables that will be written in the output.                        !
      !------------------------------------------------------------------------------------!
      !----- Domain variables. ------------------------------------------------------------!
      xa   = 1+edgeoff
      xz   = ssxp(1)-edgeoff
      ya   = 1+edgeoff
      yz   = ssyp(1)-edgeoff
      mxo  = ssxp(1)-2*edgeoff
      myo  = ssyp(1)-2*edgeoff
      dlon = grid_g(1)%lon(xa+1,1)-grid_g(1)%lon(xa,1)
      dlat = grid_g(1)%lat(1,ya)-grid_g(1)%lat(1,ya+1)
      lon0 = grid_g(1)%lon(xa,1)
      lat0 = grid_g(1)%lat(1,yz)
      !------------------------------------------------------------------------------------!

      write (unit=59,fmt='(i1)') 2

      !------------------------------------------------------------------------------------!
      !     Filling the information about the first file set.                              !
      !------------------------------------------------------------------------------------!
      !----- Output file name. ------------------------------------------------------------!
      write(outname,fmt='(2a)') trim(outpref),'_OL1_'
      write (unit=59,fmt='(a)')  trim(outname)
      !----- Dimension information. -------------------------------------------------------!
      write (unit=59,fmt='(2(i5,1x),4(f10.4,1x))') mxo,myo,dlon,dlat,lon0,lat0
      write (unit=59,fmt='(i2)') nvars_ol1
      !----- Building then writing the strings with variables, times, and types. ----------!
      varistring=''
      timestring=''
      typestring=''
      do nv=1,nvars_ol1
         write(varistring,fmt='(3a)'       ) trim(varistring)//' '''//trim(vars_ol1(nv))
         write(timestring,fmt='(a,1x,f6.0)') trim(timestring),inpfrq
         write(typestring,fmt='(a,1x,i2)'  ) trim(typestring),vtype_ol1(nv)
      end do
      write (unit=59,fmt='(a)')  trim(varistring)
      write (unit=59,fmt='(a)')  trim(timestring)
      write (unit=59,fmt='(a)')  trim(typestring)
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !     Filling the information about the second file set.                             !
      !------------------------------------------------------------------------------------!
      !----- Output file name. ------------------------------------------------------------!
      write(outname,fmt='(2a)') trim(outpref),'_OL2_'
      write (unit=59,fmt='(a)')  trim(outname)
      !----- File and dimension information. ----------------------------------------------!
      write (unit=59,fmt='(2(i5,1x),4(f10.4,1x))') mxo,myo,dlon,dlat,lon0,lat0
      write (unit=59,fmt='(i2)') nvars_ol2
      !----- Building then writing the strings with variables, times, and types. ----------!
      varistring=''
      timestring=''
      typestring=''
      do nv=1,nvars_ol2
         write(varistring,fmt='(3a)'       ) trim(varistring)//' '''//trim(vars_ol2(nv))
         write(timestring,fmt='(a,1x,f6.0)') trim(timestring),inpfrq
         write(typestring,fmt='(a,1x,i2)'  ) trim(typestring),vtype_ol2(nv)
      end do
      write (unit=59,fmt='(a)')  trim(varistring)
      write (unit=59,fmt='(a)')  trim(timestring)
      write (unit=59,fmt='(a)')  trim(typestring)
      !------------------------------------------------------------------------------------!

   end select

   !----- Printing the README. ------------------------------------------------------------!
   write (unit=59,fmt=*)
   write (unit=59,fmt=*)
   write (unit=59,fmt='(a)') '!===========================================================!'
   write (unit=59,fmt='(a)') '! README                                                    !'
   write (unit=59,fmt='(a)') '!===========================================================!'
   write (unit=59,fmt='(a)') '!     The header of the meteorological driver must contain  !'
   write (unit=59,fmt='(a)') '! the following lines:                                      !'
   write (unit=59,fmt='(a)') '!                                                           !'
   write (unit=59,fmt='(a)') '! Line  1 : Banner, it will not be read;                    !'
   write (unit=59,fmt='(a)') '! Line  2 : Number of file formats, hereafter N;            !'
   write (unit=59,fmt='(a)') '! Lines 3+: For each of the N formats, add the following    !'
   write (unit=59,fmt='(a)') '!           lines, going through a-f for the first format,  !'
   write (unit=59,fmt='(a)') '!           then through a-f for the second format and so   !'
   write (unit=59,fmt='(a)') '!            on:                                            !'
   write (unit=59,fmt='(a)') '!    a. Prefixes of the file format;                        !'
   write (unit=59,fmt='(a)') '!    b. nlon, nlat, deltalon, deltalat, lon0, lat0.  If     !'
   write (unit=59,fmt='(a)') '!       lon and lat are also variables, only nlon and nlat  !'
   write (unit=59,fmt='(a)') '!       will be used;                                       !'
   write (unit=59,fmt='(a)') '!    c. Number of variables contained in this format;       !'
   write (unit=59,fmt='(a)') '!    d. List of variables for each format (see Table 1);    !'
   write (unit=59,fmt='(a)') '!    e. Frequency at which vares are updated, or the        !'
   write (unit=59,fmt='(a)') '!       constant value if the variable type is 4;           !'
   write (unit=59,fmt='(a)') '!    f. Variable type (see Table 2);                        !'
   write (unit=59,fmt='(a)') '!                                                           !'
   write (unit=59,fmt='(a)') '!===========================================================!'
   write (unit=59,fmt='(a)') '! Table 1. Variable names recognized by ED.                 !'
   write (unit=59,fmt='(a)') '!===========================================================!'
   write (unit=59,fmt='(a)') '! -> lon    -  Longitude                        [    deg]   !'
   write (unit=59,fmt='(a)') '! -> lat    -  Latitude                         [    deg]   !'
   write (unit=59,fmt='(a)') '! -> hgt    -  Reference height                 [  m AGL]   !'
   write (unit=59,fmt='(a)') '! -> tmp    -  Air temperature                  [      K]   !'
   write (unit=59,fmt='(a)') '! -> pres   -  Pressure                         [     Pa]   !'
   write (unit=59,fmt='(a)') '! -> sh     -  Specific humidity                [  kg/kg]   !'
   write (unit=59,fmt='(a)') '! -> ugrd   -  Zonal wind                       [    m/s]   !'
   write (unit=59,fmt='(a)') '! -> vgrd   -  Zonal wind                       [    m/s]   !'
   write (unit=59,fmt='(a)') '! -> prate  -  Precipitation rate               [kg/m2/s]   !'
   write (unit=59,fmt='(a)') '! -> dlwrf  -  Downward long wave radiation     [   W/m2]   !'
   write (unit=59,fmt='(a)') '! -> nbdsf  -  Near-IR beam radiation           [   W/m2]   !'
   write (unit=59,fmt='(a)') '! -> nddsf  -  Near-IR diffuse radiation        [   W/m2]   !'
   write (unit=59,fmt='(a)') '! -> vbdsf  -  Visible beam radiation           [   W/m2]   !'
   write (unit=59,fmt='(a)') '! -> vddsf  -  Visible beam radiation           [   W/m2]   !'
   write (unit=59,fmt='(a)') '!===========================================================!'
   write (unit=59,fmt='(a)') '!                                                           !'
   write (unit=59,fmt='(a)') '!===========================================================!'
   write (unit=59,fmt='(a)') '! Table 2. Variable types recognized by ED.                 !'
   write (unit=59,fmt='(a)') '!===========================================================!'
   write (unit=59,fmt='(a)') '!                                                           !'
   write (unit=59,fmt='(a)') '! 0. Read gridded data - no time interpolation;             !'
   write (unit=59,fmt='(a)') '! 1. Read gridded data - with time interpolatation;         !'
   write (unit=59,fmt='(a)') '! 2. Read gridded data that is constant in time.            !'
   write (unit=59,fmt='(a)') '!    If any of this is lon or lat, then deltalon, deltalat  !'
   write (unit=59,fmt='(a)') '!    lon0, and lat0 will be ignored;                        !'
   write (unit=59,fmt='(a)') '! 3. Read one value representing the whole grid, no time    !'
   write (unit=59,fmt='(a)') '!   interpolation;                                          !'
   write (unit=59,fmt='(a)') '! 4. Specify a constant for all polygons, constant in time. !'
   write (unit=59,fmt='(a)') '!    In this case, give the constant value at line "e"      !'
   write (unit=59,fmt='(a)') '!    instead of the frequency.                              !'
   write (unit=59,fmt='(a)') '!===========================================================!'

   !----- Close the file. -----------------------------------------------------------------!
   close(unit=59,status='keep')

   return
end subroutine ed_metd_header
!==========================================================================================!
!==========================================================================================!
