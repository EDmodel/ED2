!==========================================================================================!
!==========================================================================================!
!  RAPP. Subroutine ncep_fill_infotable                                                    !
!                                                                                          !
!     This subroutine will loop through the NCEP reanalysis files, and fill the informa-   !
! tion table with the information retrieved from the netcdf files.  Here we will follow    !
! RAMS notation as sort of standard.  After the NCEP variables are read, they will be in-  !
! line with RAMS to ease the output step.  Note that this subroutine will do something     !
! only if netcdf libraries are loaded, otherwise this subroutine will be a dummy one.      !
!------------------------------------------------------------------------------------------!
subroutine ncep_fill_infotable(year)
   use mod_maxdims    , only : maxstr           & ! intent(in)
                             , maxfiles         & ! intent(in)
                             , maxrank          ! ! intent(in)
   use mod_model      , only : ngrids           & ! intent(out)
                             , polelon          & ! intent(out)
                             , polelat          & ! intent(out)
                             , centlon          & ! intent(out)
                             , centlat          & ! intent(out)
                             , nnxp             & ! intent(out)
                             , nnyp             & ! intent(out)
                             , nnzp             & ! intent(out)
                             , nntp             & ! intent(out)
                             , xtn              & ! intent(out)
                             , ytn              & ! intent(out)
                             , this_time        & ! intent(out)
                             , zero_time        ! ! intent(out)
   use an_header      , only : nfiles           & ! intent(out)
                             , info_table       & ! intent(out)
                             , alloc_anheader   & ! subroutine
                             , nullify_anheader ! ! subroutine
   use mod_ioopts     , only : inpath           & ! intent(in)
                             , radratio         & ! intent(in)
                             , outtimes         & ! intent(out)
                             , nouttimes        & ! intent(out)
                             , missflg_int      ! ! intent(in)
   use mod_ncep       , only : ngrids_ncep      & ! intent(in)
                             , nvars_ncep       & ! intent(in)
                             , prefvars_ncep    & ! intent(in)
                             , grids_ncep       ! ! intent(in)
#if USE_NCDF
   use mod_netcdf     , only : ncid             & ! intent(out)
                             , ndimensions      & ! intent(out)
                             , nvariables       & ! intent(out)
                             , nglobals         & ! intent(out)
                             , unlimiteddimid   & ! intent(out)
                             , xtype            & ! intent(out)
                             , timeid           & ! intent(out)
                             , dummy_vname      & ! intent(out)
                             , ndims            & ! intent(out)
                             , dimids           & ! intent(out)
                             , natts            & ! intent(out)
                             , dimglobal        & ! intent(out)
                             , globalid         ! ! intent(out)
   use netcdf
   use mod_ncdf_globio, only : ncdf_load_err
#endif
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer               , intent(in)          :: year       ! Year that we are processing
#if USE_NCDF
   !----- Local variables (only if NetCDF is available). ----------------------------------!
   character(len=maxstr)                       :: flnm_full     ! Full file name
   character(len=maxstr) , dimension(maxfiles) :: fnames        ! File name list
   integer                                     :: nf            ! Various counters
   integer                                     :: nv,nt,tcnt    ! Various counters
   integer                                     :: ierr          ! Error flag.
   integer                                     :: ngrid         ! Error flag.
   logical                                     :: file_is_there ! Is the 
   !---------------------------------------------------------------------------------------!

   !----- Initialising the number of times and vertical levels ----------------------------!
   nouttimes = 0
   nnzp      = 0

   !---------------------------------------------------------------------------------------!
   !     NCEP dataset will always have three grids.                                        !
   ! 1. Gaussian, same time resolution as input, and state variables for output;           !
   ! 2. Gaussian, higher time resolution for radiation output;                             !
   ! 3. Lon/Lat, state variables that are read here.                                       !
   !---------------------------------------------------------------------------------------!
   ngrids    = ngrids_ncep

   !---------------------------------------------------------------------------------------!
   !     Making the list of files for this year:                                           !
   !---------------------------------------------------------------------------------------!
   nfiles = nvars_ncep
   do nf=1,nfiles
      write(flnm_full,fmt='(4a,i4.4,a)')                                                   &
                      trim(inpath),'/',trim(prefvars_ncep(nf)),'.',year,'.nc'
      inquire(file=flnm_full,exist=file_is_there)
      if (file_is_there) then
         fnames(nf) = flnm_full
      else
         call fatal_error('File '//trim(flnm_full)//' is missing!!!'                       &
                         ,'ncep_fill_infotable','ncep_fill_infotable.F90')
      end if
   end do

   !----- Now I will allocate the variable table accordingly ------------------------------!
   allocate (info_table(nfiles))

   !----- Loop through the files and store the variable info ------------------------------!
   fileloop: do nf=1,nfiles
      write (unit=*,fmt='(92a)') ('-',nv=1,92)
      write (unit=*,fmt='(a,1x,2a)') ' [+] Opening file :',trim(fnames(nf)),'...'

      !----- Assigning the grid to this file. ---------------------------------------------!
      ngrid = grids_ncep(nf)

      !----- Opening the file -------------------------------------------------------------!
      ierr = nf90_open(fnames(nf),NF90_NOWRITE,ncid)
      if (ierr /= NF90_NOERR) then
         call ncdf_load_err(ierr)
         call fatal_error ('Error opening the file '//trim(fnames(nf))//'...'              &
                          ,'ncep_fill_infotable','ncep_fill_infotable.F90'                   ) 
      end if


      !----- Reading the header and extracting the file main dimensions -------------------!
      ierr = nf90_inquire(ncid,ndimensions,nvariables,nglobals,unlimiteddimid)
      info_table(nf)%filename = fnames(nf)
      info_table(nf)%nvars    = nvariables
      
      !----- Getting some global arguments ------------------------------------------------!
      write (unit=*,fmt='(a)') '     - Loading global attributes and dimensions...'
      call commio_ncep(ngrid)
      info_table(nf)%ngrids    = ngrids
      info_table(nf)%ntimes    = nntp(ngrid)
      info_table(nf)%init_time = this_time(1,ngrid)

      !----- Initialising the analysis structure for this file ----------------------------!
      call nullify_anheader(info_table(nf))
      call alloc_anheader(info_table(nf))
      info_table(nf)%avail_grid(:)            = .false.
      info_table(nf)%avail_grid(ngrid)        = .true.
      info_table(nf)%file_time(1:nntp(ngrid)) = this_time(1:nntp(ngrid),ngrid)


      !----- Adding the new times into the outtimes array ---------------------------------!
      tcnt=0
      addtimeloop: do nt=1,nntp(ngrid)
         !----- I will only add times that were not there before --------------------------!
         if (nouttimes > 0) then
           if (any(outtimes(1:nouttimes)%elapsed == this_time(nt,ngrid)%elapsed))          &
              cycle addtimeloop
         end if
         tcnt=tcnt+1
         outtimes(nouttimes+tcnt) = this_time(nt,ngrid)
      end do addtimeloop
      !----- Updating time and sorting them up --------------------------------------------!
      nouttimes = nouttimes + tcnt
      call sort_time(nouttimes,outtimes(1:nouttimes))

      !----- Getting the variable information----------------------------------------------!
      write (unit=*,fmt='(a)') '     - Loading variable table...'

      !------------------------------------------------------------------------------------!
      !     We now build the NCEP variable table.                                          !
      !------------------------------------------------------------------------------------!
      do nv=1,nvariables
         call ncep_load_var_table(nv,ngrid                                                 &
                             ,info_table(nf)%varname(nv)  , info_table(nf)%npointer(nv)    &
                             ,info_table(nf)%idim_type(nv), info_table(nf)%ngrid(nv)       &
                             ,info_table(nf)%nvalues(nv)  , info_table(nf)%rank(nv)        &
                             ,info_table(nf)%dims(:,nv)   , info_table(nf)%stagger(:,nv)   )
         write (unit=*,fmt='(a,1x,a,1x,2(a,1x,i6,1x))') &
            '         [|] Retrieving:',info_table(nf)%varname(nv)                          &
           , 'idim_type=',info_table(nf)%idim_type(nv),'rank=',info_table(nf)%rank(nv)
      end do
      
      !----- Closing the file. It will be opened again later when we load the data. -------!
      ierr = nf90_close(ncid)
      if (ierr /= NF90_NOERR) then
         call ncdf_load_err(ierr)
         call fatal_error ('Error closing the file '//trim(fnames(nf))//'...'              &
                          ,'ncep_fill_infotable','ncep_fill_infotable.F90'                 ) 
      end if

      write (unit=*,fmt='(92a)') ('-',nv=1,92)
      write (unit=*,fmt=*)
   end do fileloop
   !---------------------------------------------------------------------------------------!
   !     Here we need to fix the number of vertical levels. In case all that was provided  !
   ! was surface variables, nnzp must be 1 then.                                           !
   !---------------------------------------------------------------------------------------!
   if (nnzp(1) == 0) nnzp(1) = 1
   if (nnzp(3) == 0) nnzp(3) = 1

   !---------------------------------------------------------------------------------------!
   !    Here we just build the grid 2 information. It has the same LON/LAT size as grid 1, !
   ! but with a larger time dimension.                                                     !
   !---------------------------------------------------------------------------------------!
   nnxp(2) = nnxp(1)
   nnyp(2) = nnyp(1)
   nnzp(2) = nnzp(1)
   nntp(2) = radratio*(nntp(1)-1) + 1

   return
#else
   call fatal_error('You can''t use ncep input without compiling with netcdf libraries...' &
                   ,'ncep_fill_infotable','ncep_fill_infotable.F90')
#endif

   return
end subroutine ncep_fill_infotable
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine ncep_load_var_table(nv,current_grid,varname,npointer,idim_type,ngrid,nvalues    &
                             ,rank,dims,stagger)
   use mod_maxdims    , only : maxrank       & ! intent(in)
                             , maxstr        ! ! intent(in)
#if USE_NCDF
   use netcdf
   use mod_netcdf     , only : ncid          & ! intent(in)
                             , nvariables    & ! intent(in)
                             , varid         & ! intent(in)
                             , xtype         & ! intent(in)
                             , dummy_vname   & ! intent(in)
                             , ndims         & ! intent(in)
                             , dimids        & ! intent(in)
                             , natts         & ! intent(in)
                             , idnnxp        & ! intent(in)
                             , idnnyp        & ! intent(in)
                             , idnnzp        & ! intent(in)
                             , idnntp        & ! intent(in)
                             , idnnxpst      & ! intent(in)
                             , idnnypst      & ! intent(in)
                             , idnnzpst      ! ! intent(in)
   use mod_ioopts     , only : missflg_int   & ! intent(in)
                             , missflg_real  & ! intent(in)
                             , missflg_char  ! ! intent(in)
   use mod_model      , only : nnxp          & ! intent(in)
                             , nnyp          & ! intent(in)
                             , nnzp          & ! intent(in)
                             , nntp          ! ! intent(in)
   use mod_ncdf_globio, only : ncdf_load_err ! ! intent(in)
#endif
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                    , intent(in)  :: nv, current_grid
   character(len=*)           , intent(out) :: varname
   integer                    , intent(out) :: npointer,idim_type,ngrid,nvalues,rank
   integer, dimension(maxrank), intent(out) :: dims
   logical, dimension(maxrank), intent(out) :: stagger
   !----- Local variables -----------------------------------------------------------------!
   integer                                  :: ierr
   !---------------------------------------------------------------------------------------!  
#if USE_NCDF
   ierr = nf90_inquire_variable(ncid,nv,varname,xtype,ndims,dimids,natts)
   
   !----- Crash in case something went wrong ----------------------------------------------!
   if (ierr /= NF90_NOERR) then
      call ncdf_load_err(ierr)
      call fatal_error (                                                                   &
          'Unable to get variable information for variable '//trim(varname)//'...'         &
         ,'ncep_load_var_table','ncep_fill_infotable.F90')
   end if

   !---- The following variables are just for BRAMS files, fill it with anything ----------! 
   npointer  = missflg_int
   nvalues   = missflg_int
   ngrid     = current_grid
   
   !---- Initialising dims and stagger with "default" values ------------------------------!
   stagger   = .false.
   dims      = 1
   idim_type = 99

   !---------------------------------------------------------------------------------------!
   !     Now I decide which category the variable falls in. Currently only real vectors    !
   ! and arrays and all scalars are stored with special dimension flag. Otherwise, we save !
   ! the variable with the default 99 category, which currently means "ignore me".         !
   !     Rank is the rank of each array (0 being a scalar, 1 a 1-D vector, 2 a 2-D, and so !
   ! on). The value assigned to idim_type is standardised for all models and the look-up   !
   ! table is available at ${RAPP_ROOT}/doc/dimension_table.txt.                           !
   !---------------------------------------------------------------------------------------!
   !----- Character, only scalars are considered ------------------------------------------!
   if (xtype == NF90_CHAR .and. ndims == 1) then
      idim_type = 90
      rank      = 0 
   !----- Integer, only scalars are considered --------------------------------------------!
   elseif (xtype == NF90_INT .and. ndims == 1) then
      idim_type = 80
      rank      = 0
   !----- Reals, only 1D vectors are considered -------------------------------------------!
   elseif (xtype == NF90_FLOAT .and. ndims == 1) then
      rank      = 1
      if (dimids(1) == idnnzp) then
         idim_type = 13
         rank      = 1
         dims(1)   = nnzp(ngrid)
      elseif (dimids(1) == idnnxp) then
         idim_type = 11
         rank      = 1
         dims(1)   = nnxp(ngrid)
      elseif (dimids(1) == idnnyp) then
         idim_type = 12
         rank      = 1
         dims(1)   = nnyp(ngrid)
      elseif (dimids(1) == idnntp) then
         idim_type = 14
      end if
   !----- Double precision, only 1D vectors are considered for the time being. ------------!
   elseif (xtype == NF90_DOUBLE .and. ndims == 1) then
      rank      = 1
      if (dimids(1) == idnnzp) then
         idim_type = 73
         rank      = 1
         dims(1)   = nnzp(ngrid)
      elseif (dimids(1) == idnnxp) then
         idim_type = 71
         rank      = 1
         dims(1)   = nnxp(ngrid)
      elseif (dimids(1) == idnnyp) then
         idim_type = 72
         rank      = 1
         dims(1)   = nnyp(ngrid)
      elseif (dimids(1) == idnntp) then
         idim_type = 74
      end if

   !----- Real, scalar, vectors and higher-rank arrays are considered, check everything ---!
   elseif (xtype == NF90_SHORT) then
      select case(ndims)

      !------------------------------------------------------------------------------------!
      !   3D fields.                                                                       !
      !------------------------------------------------------------------------------------!
      case (3)
         !----- Currently the only three-dimension variable known is XYT, check that. -----!
         rank = 3

         !----- X dimension ---------------------------------------------------------------!
         if (dimids(1) == idnnxp) then
            dims(1) = nnxp(ngrid)
         else 
            write (unit=*,fmt='(3(a,1x,i5,1x))')                                           &
               'Dimids(1)=',dimids(1),'IDNNXP=',idnnxp,'IDNNXPST=',idnnxpst
            call fatal_error ('X Dimension is wrong for '//trim(varname)//'!!!'            &
                             ,'ncep_load_var_table','ncep_fill_infotable.F90')
         end if

         !----- Y dimension ---------------------------------------------------------------!
         if (dimids(2) == idnnyp) then
            dims(2) = nnyp(ngrid)
         else 
            call fatal_error ('Y Dimension is wrong for '//trim(varname)//'!!!'            &
                             ,'ncep_load_var_table','ncep_fill_infotable.F90')
         end if

         !----- T dimension ---------------------------------------------------------------!
         if (dimids(3) == idnntp) then
            idim_type = 62
            dims(3) = nntp(ngrid)
         else 
            call fatal_error ('T Dimension is wrong for '//trim(varname)//'!!!'            &
                             ,'ncep_load_var_table','ncep_fill_infotable.F90')
         end if


      !------------------------------------------------------------------------------------!
      !   4D fields.                                                                       !
      !------------------------------------------------------------------------------------!
      case (4)
         rank = 4
         !----- X dimension ---------------------------------------------------------------!
         if (dimids(1) == idnnxp) then
            dims(1) = nnxp(ngrid)
         else 
            call fatal_error ('X Dimension is wrong for '//trim(varname)//'!!!'            &
                             ,'ncep_load_var_table','ncep_fill_infotable.F90')
         end if
         !----- Y dimension ---------------------------------------------------------------!
         if (dimids(2) == idnnyp) then
            dims(2) = nnyp(ngrid)
         else 
            call fatal_error ('Y Dimension is wrong for '//trim(varname)//'!!!'            &
                             ,'ncep_load_var_table','ncep_fill_infotable.F90')
         end if

         !----- 3rd dimension -------------------------------------------------------------!
         if (dimids(3) == idnnzp) then
            dims(3) = nnzp(ngrid)
         else 
            call fatal_error ('Z Dimension is wrong for '//trim(varname)//'!!!'            &
                             ,'ncep_load_var_table','ncep_fill_infotable.F90')
         end if

         !----- 4th dimension -------------------------------------------------------------!
         if (dimids(4) == idnntp) then
            idim_type = 63
            dims(4) = nntp(ngrid)
         else 
            call fatal_error ('T Dimension is wrong for '//trim(varname)//'!!!'            &
                             ,'ncep_load_var_table','ncep_fill_infotable.F90')
         end if
      end select
   end if
#else
    varname   = missflg_char
    npointer  = missflg_int
    idim_type = missflg_int
    ngrid     = missflg_int
    nvalues   = missflg_int
    rank      = missflg_int
    dims      = missflg_int
    stagger   = .false.
#endif


   return
end subroutine ncep_load_var_table
!==========================================================================================!
!==========================================================================================!
