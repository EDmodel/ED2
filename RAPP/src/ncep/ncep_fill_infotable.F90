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
subroutine ncep_fill_infotable(year,flnm)
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
                             , this_time        & ! intent(out)
                             , zero_time        ! ! intent(out)
   use an_header      , only : nfiles           & ! intent(out)
                             , info_table       & ! intent(out)
                             , alloc_anheader   & ! subroutine
                             , nullify_anheader ! ! subroutine
   use mod_ioopts     , only : outtimes         & ! intent(out)
                             , nouttimes        & ! intent(out)
                             , missflg_int      ! ! intent(in)
   use mod_ncep       , only : nvars_ncep       & ! intent(in)
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
   character(len=maxstr) , intent(in)          :: flnm       ! File prefix
#if USE_NCDF
   !----- Local variables (only if NetCDF is available). ----------------------------------!
   character(len=maxstr)                       :: flnm_full     ! Full file name
   character(len=maxstr) , dimension(maxfiles) :: fnames        ! File name list
   integer                                     :: na,nd,nf      ! Various counters
   integer                                     :: nv,nt,tcnt    ! Various counters
   integer               , pointer             :: ifm           ! Shortcut to current grid
   integer                                     :: nvbtab        ! Scratch for var count
   integer                                     :: itrim         ! String length
   integer                                     :: ierr          ! Error flag.
   integer                                     :: ngrid         ! Error flag.
   integer                                     :: ntimes        ! number of times/file
   integer                                     :: avant         ! Cleaning handle
   logical                                     :: file_is_there ! Is the 
   !---------------------------------------------------------------------------------------!

   !----- Deallocating table if necessary, it will be allocated soon. ---------------------!
   if(allocated(info_table)) deallocate(info_table)

   !----- Initialising the number of times and vertical levels ----------------------------!
   nouttimes = 0
   nnzp      = 0
   ngrids    = 2 ! NCEP dataset will always have two grids (1=Gaussian; 2=Lon/Lat) 

   !---------------------------------------------------------------------------------------!
   !     Making the list of files for this year:                                           !
   !---------------------------------------------------------------------------------------!
   nfiles = nvars_ncep
   do nf=1,nfiles
      write(flnm_full,fmt='(4a,i4.4,a)')                                                   &
                      trim(flnm),'/',trim(prefvars_ncep(nf)),'.',year,'.nc'
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
      call commio_ncep(ngrid,ntimes)
      info_table(nf)%ngrids    = 2
      info_table(nf)%ntimes    = ntimes
      info_table(nf)%init_time = this_time(1)

      !----- Initialising the analysis structure for this file ----------------------------!
      call nullify_anheader(info_table(nf))
      call alloc_anheader(info_table(nf))
      info_table(nf)%avail_grid(:)       = .false.
      info_table(nf)%avail_grid(ngrid)   = .true.
      info_table(nf)%file_time(1:ntimes) = this_time(1:ntimes)


      !----- Adding the new times into the outtimes array ---------------------------------!
      tcnt=0
      addtimeloop: do nt=1,ntimes
         !----- I will only add times that were not there before --------------------------!
         if (nouttimes > 0) then
           if (any(outtimes(1:nouttimes)%elapsed == this_time(nt)%elapsed))                &
              cycle addtimeloop
         end if
         tcnt=tcnt+1
         outtimes(nouttimes+tcnt) = this_time(nt)
      end do addtimeloop
      !----- Updating time and sorting them up --------------------------------------------!
      nouttimes = nouttimes + tcnt
      call sort_time(nouttimes,outtimes(1:nouttimes))

      !----- Getting the variable information----------------------------------------------!
      write (unit=*,fmt='(a)') ' '
      write (unit=*,fmt='(a)') '-----------------------------------------------------------'
      write (unit=*,fmt='(a)') ' - Loading variable table...'

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
   end do fileloop
   !---------------------------------------------------------------------------------------!
   !     Here we need to fix the number of vertical levels. In case all that was provided  !
   ! was surface variables, nnzp by 1 and set up any junk for ztn(1,:).                    !
   !---------------------------------------------------------------------------------------!
   where (nnzp(1:ngrids) == 0)
      nnzp(1:ngrids)  = 1
   end where
   
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
                             , idtimes       & ! intent(in)
                             , idnpatch      & ! intent(in)
                             , idnnxpst      & ! intent(in)
                             , idnnypst      & ! intent(in)
                             , idnnzpst      ! ! intent(in)
   use mod_ioopts     , only : missflg_int   & ! intent(in)
                             , missflg_real  & ! intent(in)
                             , missflg_char  ! ! intent(in)
   use mod_model      , only : nnxp          & ! intent(in)
                             , nnyp          & ! intent(in)
                             , nnzp          ! ! intent(in)
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
   character(len=maxstr)                    :: memorder,stagdim
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
   ! on). The value assigned to idim_type is standardized for all models and the look-up   !
   ! table is available at ${MEVI_ROOT}/doc/dimension_table.txt.                           !
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
         idim_type = 70
         if (dimids(1) == idnnzp) then
            idim_type = 13
            dims(1) = nnzp(ngrid)
         elseif (dimids(1) == idnnxp) then
            dims(1) = nnxp(ngrid)
         elseif (dimids(1) == idnnyp) then
            dims(1) = nnyp(ngrid)
         end if

         rank      = 1
   elseif (xtype == NF90_DOUBLE .and. ndims == 1) then
         idim_type = 60
         rank      = 0
   !----- Real, scalar, vectors and higher-rank arrays are considered, check everything ---!
   elseif (xtype == NF90_SHORT) then
      select case(ndims)
      !------------------------------------------------------------------------------------!
      !   Time-dependent scalar or time-independent vector, just set up idim_type and we   !
      ! are all set.                                                                       !
      !------------------------------------------------------------------------------------!
      case (1)
         idim_type = 70
         rank      = 0

      !------------------------------------------------------------------------------------!
      !   Single rank vector, checking which dimension it is associated. This will not     !
      ! crash if the vector has a non-spatial dimension, we will simply ignore it and      !
      ! leave the dimension to be 99.                                                      !
      !------------------------------------------------------------------------------------!
      case (2) !----- The variable has one dimension. We will figure out which one... -----!
         rank      = 1
         if (dimids(1) == idnnzp) then
            idim_type = 13
            dims(1) = nnzp(ngrid)
         elseif (dimids(1) == idnnxp) then
            idim_type = 11
            dims(1) = nnxp(ngrid)
         elseif (dimids(1) == idnnyp) then
            idim_type = 12
            dims(1) = nnyp(ngrid)
         end if

      !------------------------------------------------------------------------------------!
      !   2D fields. Currently only XY variables are known...                              !
      !------------------------------------------------------------------------------------!
      case (3)
         rank = 2
         !----- Currently the only 2D dimension is X and Y. -------------------------------!
         idim_type = 62

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


      !------------------------------------------------------------------------------------!
      !   3D fields. There are some possible idim_type, but all of them  have X and Y,     !
      !              what changes is the third dimension.                                  !
      !------------------------------------------------------------------------------------!
      case (4)
         rank = 3
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
            idim_type = 63
            dims(3) = nnzp(ngrid)
         else 
            call fatal_error ('Z Dimension is wrong for '//trim(varname)//'!!!'            &
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
