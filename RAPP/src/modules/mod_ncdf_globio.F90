!==========================================================================================!
!==========================================================================================!
! RAPP. Module mod_ncdf_globio. This module contains two functions that handle reading     !
!       global dimensions and global attributes. It will only exist when netcdf libraries  !
!       are available.                                                                     !
!==========================================================================================!
!==========================================================================================!
module mod_ncdf_globio
#if USE_NCDF
   contains
!==========================================================================================!
!==========================================================================================!
!   This function retrieves an array of values from the global variable attribute list.    !
!------------------------------------------------------------------------------------------!
   integer function ncio_glo_sca(varname,crashmissing,intout,realout,charout)
      use netcdf
      use mod_ioopts, only: missflg_int, missflg_real, missflg_char
      use mod_netcdf, only: ncid,xtype,dimglobal,globalid
      implicit none
      
      !----- Input arguments --------------------------------------------------------------!
      character(len=*)               , intent(in)  :: varname
      logical                        , intent(in)  :: crashmissing
      !----- Optional output, only one of them should appear there ------------------------!
      integer              , optional, intent(out) :: intout
      real                 , optional, intent(out) :: realout
      character(len=*)     , optional, intent(out) :: charout
      !----- Local variable, assigning the correct output type. ---------------------------!
      integer, dimension(3)     , parameter   :: NF90_ICI=(/NF90_INT,NF90_FLOAT,NF90_CHAR/)
      integer                                 :: ici, ierr
      !------------------------------------------------------------------------------------!



      !----- Checking which variable was provided -----------------------------------------!
      ici = 0
      if (present(intout)) ici=1
      if (present(realout) .and. ici == 0) then
         ici = 2
      elseif (present(realout)) then
          call fatal_error('Variable '//trim(varname)//' has too many arguments...'        &
                           ,'ncio_glo','mod_ncdf_globio.F90'          ) 
      end if
      if (present(charout) .and. ici == 0) then
         ici = 3
      elseif (present(charout)) then
          call fatal_error('Variable '//trim(varname)//' has too many arguments...'        &
                           ,'ncio_glo','mod_ncdf_globio.F90'          ) 
      end if


      !----- Checking the variable information --------------------------------------------!
      ierr = NF90_inquire_attribute(ncid,NF90_GLOBAL,varname,xtype,dimglobal,globalid)

      !------------------------------------------------------------------------------------!
      !     Now I will check for possible common errors, and decide how to handle them.    !
      !------------------------------------------------------------------------------------!
      !----- 1. Variable is non-optional and it is not there. -----------------------------!
      if (ierr /= NF90_NOERR .and. crashmissing) then
         call fatal_error ('Global attribute '//trim(varname)//' wasn''t found!!!'         &
                          ,'ncio_glo','mod_ncdf_globio.F90')
      !----- 2. That's fine that the variable is not here, fill with junk and return... ---!
      elseif (ierr /= NF90_NOERR) then
         select case (ici)
           case (1) ; intout  = missflg_int
           case (2) ; realout = missflg_real
           case (3) ; charout = missflg_char
         end select
      !----- 3. Attribute does not have the same type as the output, crash it! ------------!
      elseif (xtype /= NF90_ICI(ici))  then
         call fatal_error(' Global attribute '//trim(varname)//' has a different format!'  &
                         ,'ncio_glo','mod_ncdf_globio.F90') 
      !----- 4. We are dealing with scalars now, which means that dimglobal should be 2. --!
      elseif (dimglobal /= 1) then
         call fatal_error(' Global attribute '//trim(varname)//' is not scalar!'  &
                         ,'ncio_glo','mod_ncdf_globio.F90') 
      !----- 5. Guess it's okay, I will take it... ----------------------------------------!
      else
         select case (ici)
            case (1) ; ierr = NF90_get_att(ncid,NF90_GLOBAL,varname,intout)
            case (2) ; ierr = NF90_get_att(ncid,NF90_GLOBAL,varname,realout)
            case (3) ; ierr = NF90_get_att(ncid,NF90_GLOBAL,varname,charout)
         end select 

         if (ierr /= NF90_NOERR) then
            call ncdf_load_err(ierr)
            call  fatal_error('Failed retrieving global attribute '//trim(varname)//'!!!'  &
                             ,'ncio_i_glo','mod_ncdf_globio.F90')
         else 
            select case (ici)
            case (1)
               write (unit=*,fmt='(3(a,1x),(1x,i5),a)')                                    &
                     '         [|] Retrieving :',varname,'=',intout,'...'
            case (2)
               write (unit=*,fmt='(3(a,1x),(1x,es12.5),a)')                                &
                     '         [|] Retrieving :',varname,'=',realout,'...'
            case (3)
               write (unit=*,fmt='(3(a,1x),2a)')                                           &
                     '         [|] Retrieving :',varname,'=',trim(charout),'...'
            end select 
         end if
      end if

      ncio_glo_sca = ierr

      return
   end function ncio_glo_sca
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !   This function retrieves a value from the global variable attribute list.            !
   !---------------------------------------------------------------------------------------!
   integer function ncio_glo(varname,crashmissing,intout,realout,charout)
      use netcdf
      use mod_ioopts, only: missflg_int, missflg_real, missflg_char
      use mod_netcdf, only: ncid,xtype,dimglobal,globalid
      implicit none
      
      !----- Input arguments --------------------------------------------------------------!
      character(len=*)               , intent(in)  :: varname
      logical                        , intent(in)  :: crashmissing
      !----- Optional output, only one of them should appear there ------------------------!
      integer, dimension(*), optional, intent(out) :: intout
      real,    dimension(*), optional, intent(out) :: realout
      character(len=*)     , optional, intent(out) :: charout
      !----- Local variable, assigning the correct output type. ---------------------------!
      integer, dimension(3)     , parameter   :: NF90_ICI=(/NF90_INT,NF90_FLOAT,NF90_CHAR/)
      integer                                 :: ici, ierr
      character(len=32)                       :: fmtdump
      !------------------------------------------------------------------------------------!



      !----- Checking which variable was provided -----------------------------------------!
      ici = 0
      if (present(intout)) ici=1
      if (present(realout) .and. ici == 0) then
         ici = 2
      elseif (present(realout)) then
          call fatal_error('Variable '//trim(varname)//' has too many arguments...'        &
                           ,'ncio_glo','mod_ncdf_globio.F90'          ) 
      end if
      if (present(charout) .and. ici == 0) then
         ici = 3
      elseif (present(charout)) then
          call fatal_error('Variable '//trim(varname)//' has too many arguments...'        &
                           ,'ncio_glo','mod_ncdf_globio.F90'          ) 
      end if


      !----- Checking the variable information --------------------------------------------!
      ierr = NF90_inquire_attribute(ncid,NF90_GLOBAL,varname,xtype,dimglobal,globalid)

      !------------------------------------------------------------------------------------!
      !     Now I will check for possible common errors, and decide how to handle them.    !
      !------------------------------------------------------------------------------------!
      !----- 1. Variable is non-optional and it is not there. -----------------------------!
      if (ierr /= NF90_NOERR .and. crashmissing) then
         call fatal_error ('Global attribute '//trim(varname)//' wasn''t found!!!'         &
                          ,'ncio_glo','mod_ncdf_globio.F90')
      !----- 2. That's fine that the variable is not here, fill with junk and return... ---!
      elseif (ierr /= NF90_NOERR) then
         select case (ici)
           case (1) ; intout(1:dimglobal)  = missflg_int
           case (2) ; realout(1:dimglobal) = missflg_real
           case (3) ; charout              = missflg_char
         end select
      !----- 3. Attribute does not have the same type as the output, crash it! ------------!
      elseif (xtype /= NF90_ICI(ici))  then
         call fatal_error(' Global attribute '//trim(varname)//' is not integer!'          &
                         ,'ncio_glo','mod_ncdf_globio.F90') 
      !----- 5. Guess it's okay, I will take it... ----------------------------------------!
      else
         select case (ici)
            case (1) ; ierr = NF90_get_att(ncid,NF90_GLOBAL,varname,intout(1:dimglobal))
            case (2) ; ierr = NF90_get_att(ncid,NF90_GLOBAL,varname,realout(1:dimglobal))
            case (3) ; ierr = NF90_get_att(ncid,NF90_GLOBAL,varname,charout)
         end select 

         if (ierr /= NF90_NOERR) then
            call ncdf_load_err(ierr)
            call  fatal_error('Failed retrieving global attribute '//trim(varname)//'!!!'  &
                             ,'ncio_i_glo','mod_ncdf_globio.F90')
         else 
            select case (ici)
            case (1)
               write (fmtdump,fmt='(a,i6.6,a)') '(3(a,1x),',dimglobal,'(1x,i5),a)'
               write (unit=*,fmt=trim(fmtdump))                                            &
                     '         [|] Retrieving :',varname,'=',intout(1:dimglobal),'...'
            case (2)
               write (fmtdump,fmt='(a,i6.6,a)') '(3(a,1x),',dimglobal,'(1x,es12.5),a)'
               write (unit=*,fmt=trim(fmtdump))                                            &
                     '         [|] Retrieving :',varname,'=',realout(1:dimglobal),'...'
            case (3)
               write (unit=*,fmt='(3(a,1x),2a)')                                           &
                     '         [|] Retrieving :',varname,'=',trim(charout),'...'
            end select 
         end if
      end if

      ncio_glo = ierr

      return
   end function ncio_glo
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !   This function retrieves an integer value from the dimension list. Integer only...   !
   !---------------------------------------------------------------------------------------!
   integer function ncio_dim(varname,crashmissing,intout,thisid)
      use netcdf
      use mod_ioopts, only: missflg_int
      use mod_netcdf, only: ncid,dummy_vname,dimid
      implicit none
      
      !----- Input arguments --------------------------------------------------------------!
      character(len=*), intent(in)            :: varname
      logical         , intent(in)            :: crashmissing
      !----- Outputs: the dimension value and the corresponding ID (optional) -------------!
      integer         , intent(out)           :: intout
      integer         , intent(out), optional :: thisid
      !----- Local variables --------------------------------------------------------------!
      integer                                 :: ierr
      !------------------------------------------------------------------------------------!

      !----- Checking where the variable I want to retrieve is ----------------------------!
      ierr = NF90_inq_dimid(ncid,varname,dimid)

      !------------------------------------------------------------------------------------!
      !     Now I will check for possible common errors, and decide how to handle them.    !
      !------------------------------------------------------------------------------------!
      !----- 1. Variable is non-optional and it is not there. -----------------------------!
      if (ierr /= NF90_NOERR .and. crashmissing) then
         call fatal_error ('Dimension '//trim(varname)//' wasn''t found!!!'                &
                          ,'ncio_dim','mod_ncdf_globio.F90')
      !----- 2. That's fine that the variable is not here, fill with junk and return... ---!
      elseif (ierr /= NF90_NOERR) then
         intout  = missflg_int
      !----- 3. Guess it's okay, I will take it... ----------------------------------------!
      else
         ierr = NF90_inquire_dimension(ncid,dimid,dummy_vname,intout)

         if (ierr /= NF90_NOERR) then
            call ncdf_load_err(ierr)
            call  fatal_error('Failed retrieving dimension '//trim(varname)//'!!!'         &
                             ,'ncio_dim','mod_ncdf_globio.F90')
         else 
            write (unit=*,fmt='(3(a,1x),i5,a)')                                            &
                  '         [|] Retrieving :',varname,'=',intout,'...'
         end if
      end if

      if (present(thisid)) thisid = dimid

      ncio_dim = ierr

      return
   end function ncio_dim
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !   This function retrieves the time variable with a single dimension. This should be   !
   ! used to retrieve WRF time based on the "time" variable only.                          !
   !---------------------------------------------------------------------------------------!
   integer function ncio_wrf_time(varname,crashmissing,ntimes,timelength,timeout)
      use netcdf
      use mod_time,   only: time_stt
      use mod_netcdf, only: ncid,varid,dummy_vname,dimid,xtype,ndims,dimids,natts
      use mod_ioopts, only: missflg_int,missflg_char,missflg_real,missflg_dble
      implicit none
      
      !----- Input arguments --------------------------------------------------------------!
      character(len=*)                 , intent(in)  :: varname
      logical                          , intent(in)  :: crashmissing
      integer                          , intent(in)  :: ntimes
      integer                          , intent(in)  :: timelength
      !----- Optional output, only one of them should appear there ------------------------!
      type(time_stt), dimension(ntimes), intent(out) :: timeout
      !----- External functions for time handling -----------------------------------------!
      integer          , external         :: julday,v5d_datestamp,v5d_timestamp,julday1000
      character(len=3) , external         :: monchar
      real             , external         :: day_fraction
      character(len=17), external         :: grads_dtstamp
      character(len=19), external         :: rapp_dtstamp
      !----- Local variables --------------------------------------------------------------!
      character(len=timelength), dimension(ntimes) :: tmpvar
      integer                                      :: ierr
      integer                                      :: t
      !------------------------------------------------------------------------------------!
      
      
      !----- Setting a default, non-sense value. ------------------------------------------!
      timeout(:)%year       = missflg_int
      timeout(:)%month      = missflg_int
      timeout(:)%day        = missflg_int
      timeout(:)%hour       = missflg_int
      timeout(:)%minu       = missflg_int
      timeout(:)%seco       = missflg_int
      timeout(:)%fracday    = missflg_real
      timeout(:)%doy        = missflg_int
      timeout(:)%mmm        = missflg_char
      timeout(:)%yyyyddd    = missflg_int
      timeout(:)%hhmmss     = missflg_int
      timeout(:)%gradsstamp = missflg_char
      timeout(:)%timestr    = missflg_char
      timeout(:)%elapsed    = missflg_dble



      !----- Checking where the variable I want to retrieve is ----------------------------!
      ierr = nf90_inq_varid(ncid,varname,varid)

      !------------------------------------------------------------------------------------!
      !     Now I will check for possible common errors, and decide how to handle them.    !
      !------------------------------------------------------------------------------------!
      !----- 1. Variable is non-optional and it is not there. -----------------------------!
      if (ierr /= NF90_NOERR .and. crashmissing) then
         call fatal_error ('Time variable  '//trim(varname)//' wasn''t found!!!'           &
                          ,'ncio_wrf_time','mod_ncdf_globio.F90')
      !----- 2. Guess it's okay, now I need to check whether other problems exist. --------!
      else
         ierr = nf90_inquire_variable(ncid,varid,dummy_vname,xtype,ndims,dimids,natts)
         !----- Now I check several possible problems that will keep me to get the data ---!
         if (ierr /= NF90_NOERR .and. crashmissing) then
            write (unit=*,fmt='(a,1x,2a,1x,i5)')                                           &
                  'Unable to load dimension',trim(varname),'. Error # :',ierr
            call  fatal_error('Failed retrieving dimension '//trim(varname)//'!!!'         &
                             ,'ncio_wrf_time','mod_ncdf_globio.F90')
         elseif (xtype /= NF90_CHAR) then
            call fatal_error('Variable '//trim(varname)//' is not a character!'            &
                            ,'ncio_wrf_time','mod_ncdf_globio.F90')
         elseif (ndims /= 2) then !--- 2 dimensions because one is time... ----------------!
            call fatal_error(                                                              &
              'ncio_wrf_time downloads only 2D variables and '// trim(varname)//' is not!' &
             ,'ncio_wrf_time','mod_ncdf_globio.F90')
         !----- Guess it's okay to get the variable now -----------------------------------!
         else
            write (unit=*,fmt='(2(a,1x),a)')                                               &
                  '         [|] Retrieving :',varname,'...'

            ierr = nf90_get_var(ncid,varid,tmpvar)
            if (ierr /= NF90_NOERR) then
               call ncdf_load_err(ierr)
               call fatal_error ('Error retrieving the variable '//trim(varname)//'...'    &
                                ,'ncio_wrf_time','mod_ncdf_globio.F90')
            end if


            do t=1,ntimes
               read(tmpvar(t),fmt='(i4.4,5(1x,i2.2))')                                     &
                       timeout(t)%year, timeout(t)%month, timeout(t)%day                   &
                      ,timeout(t)%hour, timeout(t)%minu , timeout(t)%seco
               timeout(t)%doy     = julday (timeout(t)%month,timeout(t)%day                &
                                           ,timeout(t)%year)
               timeout(t)%mmm     = monchar(timeout(t)%month)
               timeout(t)%fracday = day_fraction(timeout(t)%hour,timeout(t)%minu           &
                                                ,timeout(t)%seco)
               timeout(t)%yyyyddd    = v5d_datestamp(timeout(t)%year,timeout(t)%doy)
               timeout(t)%hhmmss     = v5d_timestamp(timeout(t)%hour,timeout(t)%minu       &
                                                    ,timeout(t)%seco)
               timeout(t)%gradsstamp = grads_dtstamp(timeout(t)%year,timeout(t)%mmm        &
                                                    ,timeout(t)%day,timeout(t)%hour        &
                                                    ,timeout(t)%minu)
               timeout(t)%timestr    = rapp_dtstamp(timeout(t)%year,timeout(t)%month       &
                                                   ,timeout(t)%day,timeout(t)%hour         &
                                                   ,timeout(t)%minu,timeout(t)%seco)
               timeout(t)%elapsed   = dble(julday1000(timeout(t)%month,timeout(t)%day      &
                                                     ,timeout(t)%year))                    &
                                    + dble(timeout(t)%fracday)
            end do

         end if
      end if

      ncio_wrf_time = ierr

      return
   end function ncio_wrf_time
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !   This function retrieves the time variable with a single dimension. This should be   !
   ! used to retrieve ECMWF time based on the "time" variable only.                        !
   !---------------------------------------------------------------------------------------!
   integer function ncio_ecmwf_time(varname,crashmissing,ntimes,timeout)
      use netcdf
      use mod_time,   only: time_stt
      use mod_netcdf, only: ncid,varid,dummy_vname,dimid,xtype,ndims,dimids,natts
      use mod_ioopts, only: missflg_int,missflg_char,missflg_real,missflg_dble
      use rconstants, only: hr_sec
      implicit none
      
      !----- Input arguments --------------------------------------------------------------!
      character(len=*)                 , intent(in)  :: varname
      logical                          , intent(in)  :: crashmissing
      integer                          , intent(in)  :: ntimes
      !----- Optional output, only one of them should appear there ------------------------!
      type(time_stt), dimension(ntimes), intent(out) :: timeout
      !----- External functions for time handling -----------------------------------------!
      integer          , external         :: julday,v5d_datestamp,v5d_timestamp,julday1000
      character(len=3) , external         :: monchar
      real             , external         :: day_fraction
      character(len=17), external         :: grads_dtstamp
      character(len=19), external         :: rapp_dtstamp
      !----- Local variables --------------------------------------------------------------!
      integer     , dimension(ntimes)          :: tmpvar
      real(kind=8), dimension(ntimes)          :: tdble
      integer                                  :: ierr
      integer                                  :: t
      !------------------------------------------------------------------------------------!
      
      
      !----- Setting a default, non-sense value. ------------------------------------------!
      timeout(:)%year       = missflg_int
      timeout(:)%month      = missflg_int
      timeout(:)%day        = missflg_int
      timeout(:)%hour       = missflg_int
      timeout(:)%minu       = missflg_int
      timeout(:)%seco       = missflg_int
      timeout(:)%fracday    = missflg_real
      timeout(:)%doy        = missflg_int
      timeout(:)%mmm        = missflg_char
      timeout(:)%yyyyddd    = missflg_int
      timeout(:)%hhmmss     = missflg_int
      timeout(:)%gradsstamp = missflg_char
      timeout(:)%timestr    = missflg_char
      timeout(:)%elapsed    = missflg_dble



      !----- Checking where the variable I want to retrieve is ----------------------------!
      ierr = nf90_inq_varid(ncid,varname,varid)

      !------------------------------------------------------------------------------------!
      !     Now I will check for possible common errors, and decide how to handle them.    !
      !------------------------------------------------------------------------------------!
      !----- 1. Variable is non-optional and it is not there. -----------------------------!
      if (ierr /= NF90_NOERR .and. crashmissing) then
         call fatal_error ('Time variable '//trim(varname)//' wasn''t found!!!'            &
                          ,'ncio_ecmwf_time','mod_ncdf_globio.F90')
      !----- 2. Guess it's okay, now I need to check whether other problems exist. --------!
      else
         ierr = nf90_inquire_variable(ncid,varid,dummy_vname,xtype,ndims,dimids,natts)
         !----- Now I check several possible problems that will keep me to get the data ---!
         if (ierr /= NF90_NOERR .and. crashmissing) then
            write (unit=*,fmt='(a,1x,2a,1x,i5)')                                           &
                  'Unable to load dimension',trim(varname),'. Error # :',ierr
            call  fatal_error('Failed retrieving dimension '//trim(varname)//'!!!'         &
                             ,'ncio_ecmwf_time','mod_ncdf_globio.F90')
         elseif (xtype /= NF90_INT) then
            call fatal_error('Variable '//trim(varname)//' is not integer!'                &
                            ,'ncio_ecmwf_time','mod_ncdf_globio.F90')
         elseif (ndims /= 1) then !--- Time should be the only time -----------------------!
            call fatal_error(                                                              &
              'ncio_ecmwf_time downloads only 1D variables; '//trim(varname)//' is not!'   &
             ,'ncio_ecmwf_time','mod_ncdf_globio.F90')
         !----- Guess it's okay to get the variable now -----------------------------------!
         else
            write (unit=*,fmt='(2(a,1x),a)')                                               &
                  '         [|] Retrieving :',varname,'...'

            ierr = nf90_get_var(ncid,varid,tmpvar)
            if (ierr /= NF90_NOERR) then
               call ncdf_load_err(ierr)
               call fatal_error ('Error retrieving the variable '//trim(varname)//'...'    &
                                ,'ncio_ecmwf_time','mod_ncdf_globio.F90')
            end if

            !----- ECMWF gives time in hours ----------------------------------------------!
            tdble = dble(tmpvar) * dble(hr_sec)

            do t=1,ntimes
               !---------------------------------------------------------------------------!
               !    ECMWF time origin is Jan 1, 1900, 000000GMT...                         !
               !---------------------------------------------------------------------------!
               call date_add_to(1900,1,1,0,tdble(t),'s',timeout(t)%year                    &
                               ,timeout(t)%month,timeout(t)%day,timeout(t)%hhmmss)
                               
               timeout(t)%hour = timeout(t)%hhmmss/10000
               timeout(t)%minu = mod(timeout(t)%hhmmss/100,100)
               timeout(t)%seco = mod(timeout(t)%hhmmss,100)
               timeout(t)%doy     = julday (timeout(t)%month,timeout(t)%day                &
                                           ,timeout(t)%year)
               timeout(t)%mmm     = monchar(timeout(t)%month)
               timeout(t)%fracday = day_fraction(timeout(t)%hour,timeout(t)%minu           &
                                                ,timeout(t)%seco)
               timeout(t)%yyyyddd    = v5d_datestamp(timeout(t)%year,timeout(t)%doy)
               timeout(t)%gradsstamp = grads_dtstamp(timeout(t)%year,timeout(t)%mmm        &
                                                    ,timeout(t)%day,timeout(t)%hour        &
                                                    ,timeout(t)%minu)
               timeout(t)%timestr    = rapp_dtstamp(timeout(t)%year,timeout(t)%month       &
                                                   ,timeout(t)%day,timeout(t)%hour         &
                                                   ,timeout(t)%minu,timeout(t)%seco)
               timeout(t)%elapsed   = dble(julday1000(timeout(t)%month,timeout(t)%day      &
                                                     ,timeout(t)%year))                    &
                                    + dble(timeout(t)%fracday)
            end do
         end if
      end if

      ncio_ecmwf_time = ierr

      return
   end function ncio_ecmwf_time
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !   This function retrieves the time variable with a single dimension. This should be   !
   ! used to retrieve NCEP time based on the "time" variable only.                         !
   !---------------------------------------------------------------------------------------!
   integer function ncio_ncep_time(varname,crashmissing,ntimes,timeout)
      use netcdf
      use mod_time,   only: time_stt
      use mod_netcdf, only: ncid,varid,dummy_vname,dimid,xtype,ndims,dimids,natts
      use mod_ioopts, only: missflg_int,missflg_char,missflg_real,missflg_dble
      use rconstants, only: hr_sec
      implicit none
      
      !----- Input arguments --------------------------------------------------------------!
      character(len=*)                 , intent(in)  :: varname
      logical                          , intent(in)  :: crashmissing
      integer                          , intent(in)  :: ntimes
      !----- Optional output, only one of them should appear there ------------------------!
      type(time_stt), dimension(ntimes), intent(out) :: timeout
      !----- External functions for time handling -----------------------------------------!
      integer          , external         :: julday,v5d_datestamp,v5d_timestamp,julday1000
      character(len=3) , external         :: monchar
      real             , external         :: day_fraction
      character(len=17), external         :: grads_dtstamp
      character(len=19), external         :: rapp_dtstamp
      !----- Local variables --------------------------------------------------------------!
      real(kind=8), dimension(ntimes)              :: tmpvar
      integer                                      :: ierr
      integer                                      :: t
      !------------------------------------------------------------------------------------!
      
      
      !----- Setting a default, non-sense value. ------------------------------------------!
      timeout(:)%year       = missflg_int
      timeout(:)%month      = missflg_int
      timeout(:)%day        = missflg_int
      timeout(:)%hour       = missflg_int
      timeout(:)%minu       = missflg_int
      timeout(:)%seco       = missflg_int
      timeout(:)%fracday    = missflg_real
      timeout(:)%doy        = missflg_int
      timeout(:)%mmm        = missflg_char
      timeout(:)%yyyyddd    = missflg_int
      timeout(:)%hhmmss     = missflg_int
      timeout(:)%gradsstamp = missflg_char
      timeout(:)%timestr    = missflg_char
      timeout(:)%elapsed    = missflg_dble



      !----- Checking where the variable I want to retrieve is ----------------------------!
      ierr = nf90_inq_varid(ncid,varname,varid)

      !------------------------------------------------------------------------------------!
      !     Now I will check for possible common errors, and decide how to handle them.    !
      !------------------------------------------------------------------------------------!
      !----- 1. Variable is non-optional and it is not there. -----------------------------!
      if (ierr /= NF90_NOERR .and. crashmissing) then
         call fatal_error ('Time variable '//trim(varname)//' wasn''t found!!!'            &
                          ,'ncio_ncep_time','mod_ncdf_globio.F90')
      !----- 2. Guess it's okay, now I need to check whether other problems exist. --------!
      else
         ierr = nf90_inquire_variable(ncid,varid,dummy_vname,xtype,ndims,dimids,natts)
         !----- Now I check several possible problems that will keep me to get the data ---!
         if (ierr /= NF90_NOERR .and. crashmissing) then
            write (unit=*,fmt='(a,1x,2a,1x,i5)')                                           &
                  'Unable to load dimension',trim(varname),'. Error # :',ierr
            call  fatal_error('Failed retrieving dimension '//trim(varname)//'!!!'         &
                             ,'ncio_ncep_time','mod_ncdf_globio.F90')
         elseif (xtype /= NF90_DOUBLE) then
            call fatal_error('Variable '//trim(varname)//' is not double precision real!'  &
                            ,'ncio_ncep_time','mod_ncdf_globio.F90')
         elseif (ndims /= 1) then !--- Time should be the only time -----------------------!
            call fatal_error(                                                              &
              'ncio_ncep_time downloads only 1D variables and '//trim(varname)//' is not!' &
             ,'ncio_ncep_time','mod_ncdf_globio.F90')
         !----- Guess it's okay to get the variable now -----------------------------------!
         else
            write (unit=*,fmt='(2(a,1x),a)')                                               &
                  '         [|] Retrieving :',varname,'...'

            ierr = nf90_get_var(ncid,varid,tmpvar)
            if (ierr /= NF90_NOERR) then
               call ncdf_load_err(ierr)
               call fatal_error ('Error retrieving the variable '//trim(varname)//'...'    &
                                ,'ncio_ncep_time','mod_ncdf_globio.F90')
            end if

            !----- NCEP gives time in hours -----------------------------------------------!
            tmpvar = tmpvar * hr_sec

            do t=1,ntimes
               !---------------------------------------------------------------------------!
               !    NCEP time origin is Jan 1, 01, 000000GMT. If I use this, it will give  !
               ! the actual date 1 days off because RAMS time origin is 1 and NCEP origin  !
               ! is zero. A quick way to fix it is by offsetting NCEP origin by one day... !
               ! Quick and dirty solution, but it works...                                 !
               !---------------------------------------------------------------------------!
               call date_add_to(0,12,31,0,tmpvar(t),'s',timeout(t)%year,timeout(t)%month   &
                               ,timeout(t)%day,timeout(t)%hhmmss)
                               
               timeout(t)%hour = timeout(t)%hhmmss/10000
               timeout(t)%minu = mod(timeout(t)%hhmmss/100,100)
               timeout(t)%seco = mod(timeout(t)%hhmmss,100)
               timeout(t)%doy     = julday (timeout(t)%month,timeout(t)%day                &
                                           ,timeout(t)%year)
               timeout(t)%mmm     = monchar(timeout(t)%month)
               timeout(t)%fracday = day_fraction(timeout(t)%hour,timeout(t)%minu           &
                                                ,timeout(t)%seco)
               timeout(t)%yyyyddd    = v5d_datestamp(timeout(t)%year,timeout(t)%doy)
               timeout(t)%gradsstamp = grads_dtstamp(timeout(t)%year,timeout(t)%mmm        &
                                                    ,timeout(t)%day,timeout(t)%hour        &
                                                    ,timeout(t)%minu)
               timeout(t)%timestr    = rapp_dtstamp(timeout(t)%year,timeout(t)%month       &
                                                   ,timeout(t)%day,timeout(t)%hour         &
                                                   ,timeout(t)%minu,timeout(t)%seco)
               timeout(t)%elapsed   = dble(julday1000(timeout(t)%month,timeout(t)%day      &
                                                     ,timeout(t)%year))                    &
                                    + dble(timeout(t)%fracday)
            end do

         end if
      end if

      ncio_ncep_time = ierr

      return
   end function ncio_ncep_time
   !=======================================================================================!
   !=======================================================================================!












   !=======================================================================================!
   !=======================================================================================!
   !   This function retrieves a scalar.                                                   !
   !---------------------------------------------------------------------------------------!
   integer function ncio_0dvar(varname,crashmissing,it,realout)
      use netcdf
      use mod_ioopts, only: missflg_real
      use mod_netcdf, only: ncid,varid,dummy_vname,dimid,xtype,ndims,dimids,natts
      implicit none
      
      !----- Input arguments --------------------------------------------------------------!
      character(len=*)      , intent(in)  :: varname
      logical               , intent(in)  :: crashmissing
      integer               , intent(in)  :: it
      !----- Optional output, only one of them should appear there ------------------------!
      real   ,  intent(out)               :: realout
      !----- Local variables --------------------------------------------------------------!
      integer                             :: ierr
      real,    dimension(1)               :: tmpvar
      integer, dimension(1)               :: start,count
      !------------------------------------------------------------------------------------!

      !----- Default is non-sense value in case it crashes --------------------------------!
      realout = missflg_real
      
      !----- Filling the "vectors" with the time we want to download ----------------------!
      start = (/ it /)
      count = (/  1 /)
      !----- Checking where the variable I want to retrieve is ----------------------------!
      ierr = nf90_inq_varid(ncid,varname,varid)

      !------------------------------------------------------------------------------------!
      !     Now I will check for possible common errors, and decide how to handle them.    !
      !------------------------------------------------------------------------------------!
      !----- 1. Variable is non-optional and it is not there. -----------------------------!
      if (ierr /= NF90_NOERR .and. crashmissing) then
         call ncdf_load_err(ierr)
         call fatal_error ('Scalar Variable  '//trim(varname)//' wasn''t found!!!'         &
                          ,'ncio_0dvar','mod_ncdf_globio.F90')
      !----- 2. Guess it's okay, now I need to check whether other problems exist. --------!
      else if (ierr /= NF90_NOERR .and. (.not. crashmissing)) then
         continue
      else
         ierr = nf90_inquire_variable(ncid,varid,dummy_vname,xtype,ndims,dimids,natts)
         !----- Now I check several possible problems that will keep me to get the data ---!
         if (ierr /= NF90_NOERR .and. crashmissing) then
            call ncdf_load_err(ierr)
            call  fatal_error('Failed retrieving dimension '//trim(varname)//'!!!'         &
                             ,'ncio_0dvar','mod_ncdf_globio.F90')
         elseif (xtype /= NF90_FLOAT) then
            call fatal_error('Variable '//trim(varname)//' is not real!'                   &
                            ,'ncio_0dvar','mod_ncdf_globio.F90')
         elseif (ndims /= 1) then !--- 1 dimension because one is time... -----------------!
            call fatal_error('ncio_0dvar downloads only scalar variables and '//           &
                             trim(varname)//' is not!','ncio_0dvar','mod_ncdf_globio.F90')
         !----- Guess it's okay to get the variable now -----------------------------------!
         else
            ierr = nf90_get_var(ncid,varid,tmpvar,start=start,count=count)
            realout = tmpvar(1)
            write (unit=*,fmt='(2(a,1x),a,1x,es12.5)')                                     &
                  '         [|] Retrieving :',varname,'...',realout
            if (ierr /= NF90_NOERR) then
               call ncdf_load_err(ierr)
               call fatal_error ('Error retrieving the variable '//trim(varname)//'...'    &
                                ,'ncio_0dvar','mod_ncdf_globio.F90')
            end if
         end if
      end if

      ncio_0dvar = ierr

      return
   end function ncio_0dvar
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This function retrieves a real or integer variable with a single dimension that is !
   ! time-dependent. IMPORTANT: even if the data is stored as integer, the output will be  !
   ! stored as real.                                                                       !
   !---------------------------------------------------------------------------------------!
   integer function ncio_1dnotime(varname,crashmissing,dimin,realout,offset,strin)
      use netcdf
      use mod_ioopts, only: missflg_real
      use mod_netcdf, only: ncid,varid,dummy_vname,dimid,xtype,ndims,dimids,natts
      implicit none
      
      !----- Input arguments --------------------------------------------------------------!
      character(len=*)      , intent(in)           :: varname
      logical               , intent(in)           :: crashmissing
      integer               , intent(in)           :: dimin
      integer               , intent(in), optional :: offset
      integer               , intent(in), optional :: strin
      !----- Optional output, only one of them should appear there ------------------------!
      real   , dimension(dimin), intent(out) :: realout
      !----- Local variables --------------------------------------------------------------!
      integer, dimension(1)                  :: start,count,stride
      integer                                :: ierr,off1,str1
      !------------------------------------------------------------------------------------!

      !----- Default is non-sense value in case it crashes --------------------------------!
      realout = missflg_real
      
      !----- Defining the time we are getting and the subset of data if requested. --------!
      off1=0
      if (present(offset)) off1=offset
      str1=1
      if (present(strin)) str1=strin 
      start  = (/ 1+off1 /)
      count  = (/  dimin /)
      stride = (/   str1 /)
      
      
      !----- Checking where the variable I want to retrieve is ----------------------------!
      ierr = nf90_inq_varid(ncid,varname,varid)

      !------------------------------------------------------------------------------------!
      !     Now I will check for possible common errors, and decide how to handle them.    !
      !------------------------------------------------------------------------------------!
      !----- 1. Variable is non-optional and it is not there. -----------------------------!
      if (ierr /= NF90_NOERR .and. crashmissing) then
         call ncdf_load_err(ierr)
         call fatal_error ('1D-Variable  '//trim(varname)//' wasn''t found!!!'             &
                          ,'ncio_1dvar','mod_ncdf_globio.F90')
      !----- 2. Guess it's okay, now I need to check whether other problems exist. --------!
      else if (ierr /= NF90_NOERR .and. (.not. crashmissing)) then
         continue
      else
         ierr = nf90_inquire_variable(ncid,varid,dummy_vname,xtype,ndims,dimids,natts)
         !----- Now I check several possible problems that will keep me to get the data ---!
         if (ierr /= NF90_NOERR .and. crashmissing) then
            call ncdf_load_err(ierr)
            call  fatal_error('Failed retrieving dimension '//trim(varname)//'!!!'         &
                             ,'ncio_1dvar','mod_ncdf_globio.F90')
         elseif (xtype /= NF90_FLOAT .and. xtype /= NF90_INT) then
            call fatal_error('Variable '//trim(varname)//' is not real or integer!'        &
                            ,'ncio_1dnotime','mod_ncdf_globio.F90')
         elseif (ndims /= 1) then !--- 1 because this is not time-dependent ---------------!
            call fatal_error('ncio_1dnotime downloads only time independent 1D vars and '//&
                            trim(varname)//' is not!','ncio_1dnotime','mod_ncdf_globio.F90')
         !----- Guess it's okay to get the variable now -----------------------------------!
         else
            ierr = nf90_get_var(ncid,varid,realout,start=start,count=count,stride=stride)
            write (unit=*,fmt='(2(a,1x),2a,2(1x,es14.7))')                                 &
                  '         [|] Retrieving :',trim(varname),'...'                          &
                             ,' range=',minval(realout),maxval(realout)
            if (ierr /= NF90_NOERR) then
               call ncdf_load_err(ierr)
               call fatal_error ('Error retrieving the variable '//trim(varname)//'...'    &
                                ,'ncio_1dnotime','mod_ncdf_globio.F90')
            end if
         end if
      end if

      ncio_1dnotime = ierr

      return
   end function ncio_1dnotime
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !   This function retrieves a real variable with a single dimension that is time-depen- !
   ! dent.                                                                                 !
   !---------------------------------------------------------------------------------------!
   integer function ncio_1dvar(varname,crashmissing,it,dimin,realout,offset,strin)
      use netcdf
      use mod_ioopts, only: missflg_real
      use mod_netcdf, only: ncid,varid,dummy_vname,dimid,xtype,ndims,dimids,natts
      implicit none
      
      !----- Input arguments --------------------------------------------------------------!
      character(len=*)      , intent(in)           :: varname
      logical               , intent(in)           :: crashmissing
      integer               , intent(in)           :: it
      integer               , intent(in)           :: dimin
      integer               , intent(in), optional :: offset
      integer               , intent(in), optional :: strin
      !----- Optional output, only one of them should appear there ------------------------!
      real   , dimension(dimin), intent(out) :: realout
      !----- Local variables --------------------------------------------------------------!
      integer, dimension(2)                  :: start,count,stride
      integer                                :: ierr,off1,str1
      !------------------------------------------------------------------------------------!

      !----- Default is non-sense value in case it crashes --------------------------------!
      realout = missflg_real
      
      !----- Defining the time we are getting and the subset of data if requested. --------!
      off1=0
      if (present(offset)) off1=offset
      str1=1
      if (present(strin)) str1=strin 
      start  = (/ 1+off1, it /)
      count  = (/  dimin,  1 /)
      stride = (/   str1,  1 /)
      
      
      !----- Checking where the variable I want to retrieve is ----------------------------!
      ierr = nf90_inq_varid(ncid,varname,varid)

      !------------------------------------------------------------------------------------!
      !     Now I will check for possible common errors, and decide how to handle them.    !
      !------------------------------------------------------------------------------------!
      !----- 1. Variable is non-optional and it is not there. -----------------------------!
      if (ierr /= NF90_NOERR .and. crashmissing) then
         call ncdf_load_err(ierr)
         call fatal_error ('1D-Variable  '//trim(varname)//' wasn''t found!!!'             &
                          ,'ncio_1dvar','mod_ncdf_globio.F90')
      !----- 2. Guess it's okay, now I need to check whether other problems exist. --------!
      else if (ierr /= NF90_NOERR .and. (.not. crashmissing)) then
         continue
      else
         ierr = nf90_inquire_variable(ncid,varid,dummy_vname,xtype,ndims,dimids,natts)
         !----- Now I check several possible problems that will keep me to get the data ---!
         if (ierr /= NF90_NOERR .and. crashmissing) then
            call ncdf_load_err(ierr)
            call  fatal_error('Failed retrieving dimension '//trim(varname)//'!!!'         &
                             ,'ncio_1dvar','mod_ncdf_globio.F90')
         elseif (xtype /= NF90_FLOAT) then
            call fatal_error('Variable '//trim(varname)//' is not real!'                   &
                            ,'ncio_1dvar','mod_ncdf_globio.F90')
         elseif (ndims /= 2) then !--- 2 dimensions because one is time... ----------------!
            call fatal_error('ncio_1dvar downloads only 1D variables and '//               &
                             trim(varname)//' is not!','ncio_1dvar','mod_ncdf_globio.F90')
         !----- Guess it's okay to get the variable now -----------------------------------!
         else
            write (unit=*,fmt='(2(a,1x),2a,2(1x,es14.7))')                                 &
                  '         [|] Retrieving :',trim(varname),'...'                          &
                             ,' range=',minval(realout),maxval(realout)
            ierr = nf90_get_var(ncid,varid,realout,start=start,count=count,stride=stride)
            if (ierr /= NF90_NOERR) then
               call ncdf_load_err(ierr)
               call fatal_error ('Error retrieving the variable '//trim(varname)//'...'    &
                                ,'ncio_1dvar','mod_ncdf_globio.F90')
            end if
         end if
      end if

      ncio_1dvar = ierr

      return
   end function ncio_1dvar
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !   This function retrieves a real variable with rank 3. The 2D in the name is because  !
   ! one dimension is time, and we are downloading one time per call.                      !
   !---------------------------------------------------------------------------------------!
   integer function ncio_2dvar(varname,crashmissing,it,dim1,dim2,realout,offin1,offin2     &
                             ,strin1,strin2)
      use netcdf
      use mod_ioopts, only: missflg_real 
      use mod_netcdf, only: ncid,varid,dummy_vname,dimid,xtype,ndims,dimids,natts
      implicit none
      
      !----- Input arguments --------------------------------------------------------------!
      character(len=*)          , intent(in)           :: varname
      logical                   , intent(in)           :: crashmissing
      integer                   , intent(in)           :: it
      integer                   , intent(in)           :: dim1, dim2
      integer                   , intent(in), optional :: offin1, offin2
      integer                   , intent(in), optional :: strin1, strin2
      !----- Optional output, only one of them should appear there ------------------------!
      real, dimension(dim1,dim2), intent(out) :: realout
      !----- Local variables --------------------------------------------------------------!
      integer                                 :: ierr,off1,off2,str1,str2
      integer, dimension(3)                   :: start,count,stride
      !------------------------------------------------------------------------------------!

      !----- Default is non-sense value in case it crashes --------------------------------!
      realout = missflg_real

      !----- Checking whether the optional input variables are present --------------------!
      off1 = 0
      off2 = 0
      if (present(offin1)) off1 = offin1
      if (present(offin2)) off2 = offin2

      str1 = 1
      str2 = 1
      if (present(strin1)) str1 = strin1
      if (present(strin2)) str2 = strin2
      
      !----- Defining the dimensions so it gets only 1 time -------------------------------!
      start   = (/ 1+off1, 1+off2,   it /)
      count   = (/   dim1,   dim2,    1 /)
      stride  = (/   str1,   str2,    1 /)


      !----- Checking where the variable I want to retrieve is ----------------------------!
      ierr = nf90_inq_varid(ncid,varname,varid)

      !------------------------------------------------------------------------------------!
      !     Now I will check for possible common errors, and decide how to handle them.    !
      !------------------------------------------------------------------------------------!
      !----- 1. Variable is non-optional and it is not there. -----------------------------!
      if (ierr /= NF90_NOERR .and. crashmissing) then
         call ncdf_load_err(ierr)
         call fatal_error ('2D-Variable  '//trim(varname)//' wasn''t found!!!'             &
                          ,'ncio_2dvar','mod_ncdf_globio.F90')
      !----- 2. Guess it's okay, now I need to check whether other problems exist. --------!
      else if (ierr /= NF90_NOERR .and. (.not. crashmissing)) then
         continue
      else
         ierr = nf90_inquire_variable(ncid,varid,dummy_vname,xtype,ndims,dimids,natts)
         !----- Now I check several possible problems that will keep me to get the data ---!
         if (ierr /= NF90_NOERR .and. crashmissing) then
            call ncdf_load_err(ierr)
            call  fatal_error('Failed retrieving dimension '//trim(varname)//'!!!'         &
                             ,'ncio_2dvar','mod_ncdf_globio.F90')
         elseif (xtype /= NF90_FLOAT) then
            call fatal_error('Variable '//trim(varname)//' is not real!'                   &
                            ,'ncio_2dvar','mod_ncdf_globio.F90')
         elseif (ndims /= 3) then !--- 3 dimensions because one is time... ----------------!
            call fatal_error('ncio_2dvar downloads only 2D variables and '//               &
                             trim(varname)//' is not!','ncio_2dvar','mod_ncdf_globio.F90')
         !----- Guess it's okay to get the variable now -----------------------------------!
         else
            ierr = nf90_get_var(ncid,varid,realout,start=start,count=count,stride=stride)
            write (unit=*,fmt='(2(a,1x),2a,2(1x,es14.7))')                                 &
                  '         [|] Retrieving :',trim(varname),'...'                          &
                             ,' range=',minval(realout),maxval(realout)
            if (ierr /= NF90_NOERR) then
               write (unit=*,fmt=*) 'dimensions=',dim1,dim2
               write (unit=*,fmt=*) 'start=',start
               write (unit=*,fmt=*) 'count=',count
               write (unit=*,fmt=*) 'stride=',stride
               call ncdf_load_err(ierr)
               call fatal_error ('Error retrieving the variable '//trim(varname)//'...'    &
                                ,'ncio_2dvar','mod_ncdf_globio.F90')
            end if
         end if
      end if

      ncio_2dvar = ierr

      return
   end function ncio_2dvar
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !   This function retrieves a real variable with rank 4. The 3D in the name is because  !
   ! one dimension is time, and we are downloading one time per call.                      !
   !---------------------------------------------------------------------------------------!
   integer function ncio_3dvar(varname,crashmissing,it,dim1,dim2,dim3,realout,offin1       &
                              ,offin2,offin3,strin1,strin2,strin3)
      use netcdf
      use mod_ioopts, only: missflg_real
      use mod_netcdf, only: ncid,varid,dummy_vname,dimid,xtype,ndims,dimids,natts
      implicit none
      
      !----- Input arguments --------------------------------------------------------------!
      character(len=*)               , intent(in)           :: varname
      logical                        , intent(in)           :: crashmissing
      integer                        , intent(in)           :: it
      integer                        , intent(in)           :: dim1,dim2,dim3
      integer                        , intent(in), optional :: offin1, offin2, offin3
      integer                        , intent(in), optional :: strin1, strin2, strin3
      !----- Optional output, only one of them should appear there ------------------------!
      real, dimension(dim1,dim2,dim3), intent(out) :: realout
      !----- Local variables --------------------------------------------------------------!
      integer                                      :: ierr,off1,off2,off3,str1,str2,str3
      integer, dimension(4)                        :: start,count,stride
      !------------------------------------------------------------------------------------!

      !----- Default is non-sense value in case it crashes --------------------------------!
      realout = missflg_real

      !----- Checking whether the optional input variables are present --------------------!
      off1 = 0
      off2 = 0
      off3 = 0
      if (present(offin1)) off1 = offin1
      if (present(offin2)) off2 = offin2
      if (present(offin3)) off3 = offin3

      str1 = 1
      str2 = 1
      str3 = 1
      if (present(strin1)) str1 = strin1
      if (present(strin2)) str2 = strin2
      if (present(strin3)) str3 = strin3

      !----- Defining the dimensions so it gets only 1 time -------------------------------!
      start   = (/ 1+off1, 1+off2, 1+off3,   it /)
      count   = (/   dim1,   dim2,   dim3,    1 /)
      stride  = (/   str1,   str2,   str3,    1 /)
      !----- Checking where the variable I want to retrieve is ----------------------------!
      ierr = nf90_inq_varid(ncid,varname,varid)

      !------------------------------------------------------------------------------------!
      !     Now I will check for possible common errors, and decide how to handle them.    !
      !------------------------------------------------------------------------------------!
      !----- 1. Variable is non-optional and it is not there. -----------------------------!
      if (ierr /= NF90_NOERR .and. crashmissing) then
         call ncdf_load_err(ierr)
         call fatal_error ('3D-Variable  '//trim(varname)//' wasn''t found!!!'             &
                          ,'ncio_3dvar','mod_ncdf_globio.F90')
      !----- 2. Guess it's okay, now I need to check whether other problems exist. --------!
      else if (ierr /= NF90_NOERR .and. (.not. crashmissing)) then
         continue
      else
         ierr = nf90_inquire_variable(ncid,varid,dummy_vname,xtype,ndims,dimids,natts)
         !----- Now I check several possible problems that will keep me to get the data ---!
         if (ierr /= NF90_NOERR .and. crashmissing) then
            call ncdf_load_err(ierr)
            call  fatal_error('Failed retrieving dimension '//trim(varname)//'!!!'         &
                             ,'ncio_3dvar','mod_ncdf_globio.F90')
         elseif (xtype /= NF90_FLOAT) then
            call fatal_error('Variable '//trim(varname)//' is not real!'                   &
                            ,'ncio_3dvar','mod_ncdf_globio.F90')
         elseif (ndims /= 4) then !--- 4 dimensions because one is time... ----------------!
            call fatal_error('ncio_3dvar downloads only 3D variables and '//               &
                             trim(varname)//' is not!','ncio_3dvar','mod_ncdf_globio.F90')
         !----- Guess it's okay to get the variable now -----------------------------------!
         else
            ierr = nf90_get_var(ncid,varid,realout,start=start,count=count,stride=stride)
            write (unit=*,fmt='(2(a,1x),2a,2(1x,es14.7))')                                 &
                  '         [|] Retrieving :',trim(varname),'...'                          &
                             ,' range=',minval(realout),maxval(realout)
            if (ierr /= NF90_NOERR) then
               write (unit=*,fmt=*) 'dimensions=',dim1,dim2
               write (unit=*,fmt=*) 'start=',start
               write (unit=*,fmt=*) 'count=',count
               write (unit=*,fmt=*) 'stride=',stride
               call ncdf_load_err(ierr)
               call fatal_error ('Error retrieving the variable '//trim(varname)//'...'    &
                                ,'ncio_3dvar','mod_ncdf_globio.F90')
            end if
         end if
      end if

      ncio_3dvar = ierr

      return
   end function ncio_3dvar
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !   This function retrieves a real variable with rank 5. The 4D in the name is because  !
   ! one dimension is time, and we are downloading one time per call.                      !
   !---------------------------------------------------------------------------------------!
   integer function ncio_4dvar(varname,crashmissing,it,dim1,dim2,dim3,dim4,realout,offin1  &
                             ,offin2,offin3,offin4,strin1,strin2,strin3,strin4)
      use netcdf
      use mod_ioopts, only: missflg_real
      use mod_netcdf, only: ncid,varid,dummy_vname,dimid,xtype,ndims,dimids,natts
      implicit none
      
      !----- Input arguments --------------------------------------------------------------!
      character(len=*)               , intent(in)           :: varname
      logical                        , intent(in)           :: crashmissing
      integer                        , intent(in)           :: it
      integer                        , intent(in)           :: dim1,dim2,dim3,dim4
      integer                        , intent(in), optional :: offin1,offin2,offin3,offin4
      integer                        , intent(in), optional :: strin1,strin2,strin3,strin4
      !----- Optional output, only one of them should appear there ------------------------!
      real, dimension(dim1,dim2,dim3,dim4), intent(out)     :: realout
      !----- Local variables --------------------------------------------------------------!
      integer                                               :: ierr,off1,off2,off3,off4
      integer                                               :: str1,str2,str3,str4
      integer, dimension(5)                                 :: start,count,stride
      !------------------------------------------------------------------------------------!

      !----- Default is non-sense value in case it crashes --------------------------------!
      realout = missflg_real

      !----- Checking whether the optional input variables are present --------------------!
      off1 = 0
      off2 = 0
      off3 = 0
      off4 = 0
      if (present(offin1)) off1 = offin1
      if (present(offin2)) off2 = offin2
      if (present(offin3)) off3 = offin3
      if (present(offin4)) off4 = offin4

      str1 = 1
      str2 = 1
      str3 = 1
      str4 = 1
      if (present(strin1)) str1 = strin1
      if (present(strin2)) str2 = strin2
      if (present(strin3)) str3 = strin3
      if (present(strin4)) str4 = strin4

      !----- Defining the dimensions so it gets only 1 time -------------------------------!
      start   = (/ 1+off1, 1+off2, 1+off3, 1+off4,   it /)
      count   = (/   dim1,   dim2,   dim3,   dim4,    1 /)
      stride  = (/   str1,   str2,   str3,   str4,    1 /)
      !----- Checking where the variable I want to retrieve is ----------------------------!
      ierr = nf90_inq_varid(ncid,varname,varid)

      !------------------------------------------------------------------------------------!
      !     Now I will check for possible common errors, and decide how to handle them.    !
      !------------------------------------------------------------------------------------!
      !----- 1. Variable is non-optional and it is not there. -----------------------------!
      if (ierr /= NF90_NOERR .and. crashmissing) then
         call ncdf_load_err(ierr)
         call fatal_error ('4D-Variable  '//trim(varname)//' wasn''t found!!!'             &
                          ,'ncio_4dvar','mod_ncdf_globio.F90')
      !----- 2. Guess it's okay, now I need to check whether other problems exist. --------!
      else if (ierr /= NF90_NOERR .and. (.not. crashmissing)) then
         continue
      else
         ierr = nf90_inquire_variable(ncid,varid,dummy_vname,xtype,ndims,dimids,natts)
         !----- Now I check several possible problems that will keep me to get the data ---!
         if (ierr /= NF90_NOERR .and. crashmissing) then
            call ncdf_load_err(ierr)
            call  fatal_error('Failed retrieving dimension '//trim(varname)//'!!!'         &
                             ,'ncio_4dvar','mod_ncdf_globio.F90')
         elseif (xtype /= NF90_FLOAT) then
            call fatal_error('Variable '//trim(varname)//' is not real!'                   &
                            ,'ncio_4dvar','mod_ncdf_globio.F90')
         elseif (ndims /= 5) then !--- 5 dimensions because one is time... ----------------!
            call fatal_error('ncio_4dvar downloads only 3D variables and '//               &
                             trim(varname)//' is not!','ncio_4dvar','mod_ncdf_globio.F90')
         !----- Guess it's okay to get the variable now -----------------------------------!
         else
            write (unit=*,fmt='(2(a,1x),a)')                                               &
                  '         [|] Retrieving :',trim(varname),'...'
            ierr = nf90_get_var(ncid,varid,realout,start=start,count=count,stride=stride)
            if (ierr /= NF90_NOERR) then
               call ncdf_load_err(ierr)
               call fatal_error ('Error retrieving the variable '//trim(varname)//'...'    &
                                ,'ncio_4dvar','mod_ncdf_globio.F90')
            end if
         end if
      end if

      ncio_4dvar = ierr

      return
   end function ncio_4dvar
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !   This function retrieves a real variable with rank 3. The 2D in the name is because  !
   ! one dimension is time, and we are downloading one time per call. This should be used  !
   ! to extract fields from the ECMWF netcdf (and NCEP too).                               !
   !---------------------------------------------------------------------------------------!
   integer function ncio_2dshort(varname,crashmissing,it,dim1,dim2,realout,offin1,offin2   &
                                ,strin1,strin2)
      use netcdf
      use mod_ioopts, only: missflg_real 
      use mod_netcdf, only: ncid,varid,dummy_vname,dimid,xtype,ndims,dimids,natts
      implicit none
      
      !----- Input arguments --------------------------------------------------------------!
      character(len=*)          , intent(in)           :: varname
      logical                   , intent(in)           :: crashmissing
      integer                   , intent(in)           :: it
      integer                   , intent(in)           :: dim1, dim2
      integer                   , intent(in), optional :: offin1, offin2
      integer                   , intent(in), optional :: strin1, strin2
      !----- Optional output, only one of them should appear there ------------------------!
      real, dimension(dim1,dim2), intent(out) :: realout
      !----- Local variables --------------------------------------------------------------!
      integer(kind=2), dimension(dim1,dim2)   :: tmpvar
      integer                                 :: ierr,off1,off2,str1,str2
      integer, dimension(3)                   :: start,count,stride
      integer(kind=2)                         :: missshort
      real                                    :: refval,scalefac
      !------------------------------------------------------------------------------------!

      !----- Default is non-sense value in case it crashes --------------------------------!
      realout = missflg_real

      !----- Checking whether the optional input variables are present --------------------!
      off1 = 0
      off2 = 0
      if (present(offin1)) off1 = offin1
      if (present(offin2)) off2 = offin2

      str1 = 1
      str2 = 1
      if (present(strin1)) str1 = strin1
      if (present(strin2)) str2 = strin2
      
      !----- Defining the dimensions so it gets only 1 time -------------------------------!
      start   = (/ 1+off1, 1+off2,   it /)
      count   = (/   dim1,   dim2,    1 /)
      stride  = (/   str1,   str2,    1 /)


      !----- Checking where the variable I want to retrieve is ----------------------------!
      ierr = nf90_inq_varid(ncid,varname,varid)

      !------------------------------------------------------------------------------------!
      !     Now I will check for possible common errors, and decide how to handle them.    !
      !------------------------------------------------------------------------------------!
      !----- 1. Variable is non-optional and it is not there. -----------------------------!
      if (ierr /= NF90_NOERR .and. crashmissing) then
         call ncdf_load_err(ierr)
         call fatal_error ('2D-Variable  '//trim(varname)//' wasn''t found!!!'             &
                          ,'ncio_2dshort','mod_ncdf_globio.F90')
      !----- 2. Guess it's okay, now I need to check whether other problems exist. --------!
      else if (ierr /= NF90_NOERR .and. (.not. crashmissing)) then
         continue
      else
         ierr = nf90_inquire_variable(ncid,varid,dummy_vname,xtype,ndims,dimids,natts)
         !----- Now I check several possible problems that will keep me to get the data ---!
         if (ierr /= NF90_NOERR .and. crashmissing) then
            call ncdf_load_err(ierr)
            call  fatal_error('Failed retrieving dimension '//trim(varname)//'!!!'         &
                             ,'ncio_2dshort','mod_ncdf_globio.F90')
         elseif (xtype /= NF90_SHORT) then
            call fatal_error('Variable '//trim(varname)//' is not short integer!'          &
                            ,'ncio_2dshort','mod_ncdf_globio.F90')
         elseif (ndims /= 3) then !--- 3 dimensions because one is time... ----------------!
            call fatal_error('ncio_2dshort downloads only 2D variables and '//             &
                            trim(varname)//' is not!','ncio_2dshort','mod_ncdf_globio.F90')
         !----- Guess it's okay to get the variable now -----------------------------------!
         else
            write (unit=*,fmt='(2(a,1x),a)')                                               &
                  '         [|] Retrieving :',trim(varname),'...'
            ierr = nf90_get_var(ncid,varid,tmpvar,start=start,count=count,stride=stride)
            if (ierr /= NF90_NOERR) then
               call ncdf_load_err(ierr)
               call fatal_error ('Error retrieving the variable '//trim(varname)//'...'    &
                                ,'ncio_2dshort','mod_ncdf_globio.F90')
            end if

            !----- Getting the rescaling attributes ---------------------------------------!
            ierr = nf90_get_att(ncid,varid,'missing_value',missshort)
            ierr = nf90_get_att(ncid,varid,'add_offset'   ,refval   )
            ierr = nf90_get_att(ncid,varid,'scale_factor' ,scalefac )
            
            !----- Converting the values to the real form ---------------------------------!
            where (tmpvar /= missshort)
               realout = refval + scalefac*real(tmpvar)
            end where


         end if
      end if

      ncio_2dshort = ierr

      return
   end function ncio_2dshort
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !   This function retrieves a real variable with rank 4. The 3D in the name is because  !
   ! one dimension is time, and we are downloading one time per call. This will assume     !
   ! that the variable is stored as short integer, like ECMWF reanalysis (NCEP too...).    !
   !---------------------------------------------------------------------------------------!
   integer function ncio_3dshort(varname,crashmissing,it,dim1,dim2,dim3,realout,offin1     &
                              ,offin2,offin3,strin1,strin2,strin3)
      use netcdf
      use mod_ioopts, only: missflg_real
      use mod_netcdf, only: ncid,varid,dummy_vname,dimid,xtype,ndims,dimids,natts
      implicit none
      
      !----- Input arguments --------------------------------------------------------------!
      character(len=*)               , intent(in)           :: varname
      logical                        , intent(in)           :: crashmissing
      integer                        , intent(in)           :: it
      integer                        , intent(in)           :: dim1,dim2,dim3
      integer                        , intent(in), optional :: offin1, offin2, offin3
      integer                        , intent(in), optional :: strin1, strin2, strin3
      !----- Optional output, only one of them should appear there ------------------------!
      real, dimension(dim1,dim2,dim3), intent(out) :: realout
      !----- Local variables --------------------------------------------------------------!
      integer(kind=2), dimension(dim1,dim2,dim3)   :: tmpvar
      integer                                      :: ierr,off1,off2,off3,str1,str2,str3
      integer, dimension(4)                        :: start,count,stride
      integer(kind=2)                              :: missshort
      real                                         :: refval,scalefac
      !------------------------------------------------------------------------------------!

      !----- Default is non-sense value in case it crashes --------------------------------!
      realout = missflg_real

      !----- Checking whether the optional input variables are present --------------------!
      off1 = 0
      off2 = 0
      off3 = 0
      if (present(offin1)) off1 = offin1
      if (present(offin2)) off2 = offin2
      if (present(offin3)) off3 = offin3

      str1 = 1
      str2 = 1
      str3 = 1
      if (present(strin1)) str1 = strin1
      if (present(strin2)) str2 = strin2
      if (present(strin3)) str3 = strin3

      !----- Defining the dimensions so it gets only 1 time -------------------------------!
      start   = (/ 1+off1, 1+off2, 1+off3,   it /)
      count   = (/   dim1,   dim2,   dim3,    1 /)
      stride  = (/   str1,   str2,   str3,    1 /)
      !----- Checking where the variable I want to retrieve is ----------------------------!
      ierr = nf90_inq_varid(ncid,varname,varid)

      !------------------------------------------------------------------------------------!
      !     Now I will check for possible common errors, and decide how to handle them.    !
      !------------------------------------------------------------------------------------!
      !----- 1. Variable is non-optional and it is not there. -----------------------------!
      if (ierr /= NF90_NOERR .and. crashmissing) then
         call ncdf_load_err(ierr)
         call fatal_error ('3D-Variable  '//trim(varname)//' wasn''t found!!!'             &
                          ,'ncio_3dshort','mod_ncdf_globio.F90')
      !----- 2. Guess it's okay, now I need to check whether other problems exist. --------!
      else if (ierr /= NF90_NOERR .and. (.not. crashmissing)) then
         continue
      else
         ierr = nf90_inquire_variable(ncid,varid,dummy_vname,xtype,ndims,dimids,natts)
         !----- Now I check several possible problems that will keep me to get the data ---!
         if (ierr /= NF90_NOERR .and. crashmissing) then
            call ncdf_load_err(ierr)
            call  fatal_error('Failed retrieving dimension '//trim(varname)//'!!!'         &
                             ,'ncio_3dshort','mod_ncdf_globio.F90')
         elseif (xtype /= NF90_SHORT) then
            call fatal_error('Variable '//trim(varname)//' is not short integer!'          &
                            ,'ncio_3dshort','mod_ncdf_globio.F90')
         elseif (ndims /= 4) then !--- 4 dimensions because one is time... ----------------!
            call fatal_error('ncio_3dshort downloads only 3D variables and '//             &
                            trim(varname)//' is not!','ncio_3dshort','mod_ncdf_globio.F90')
         !----- Guess it's okay to get the variable now -----------------------------------!
         else
            write (unit=*,fmt='(2(a,1x),a)')                                               &
                  '         [|] Retrieving :',trim(varname),'...'
            ierr = nf90_get_var(ncid,varid,tmpvar,start=start,count=count,stride=stride)
            if (ierr /= NF90_NOERR) then
               call ncdf_load_err(ierr)
               call fatal_error ('Error retrieving the variable '//trim(varname)//'...'    &
                                ,'ncio_3dshort','mod_ncdf_globio.F90')
            end if

            !----- Getting the rescaling attributes ---------------------------------------!
            ierr = nf90_get_att(ncid,varid,'missing_value',missshort)
            ierr = nf90_get_att(ncid,varid,'add_offset'   ,refval   )
            ierr = nf90_get_att(ncid,varid,'scale_factor' ,scalefac )
            
            !----- Converting the values to the real form ---------------------------------!
            where (tmpvar /= missshort)
               realout = refval + scalefac*real(tmpvar)
            end where

         end if
      end if

      ncio_3dshort = ierr

      return
   end function ncio_3dshort
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine checks the error and print the most suitable message for the       !
   ! unfortunate user.                                                                     !
   !---------------------------------------------------------------------------------------!
   subroutine ncdf_load_err(ierr)
      use netcdf
      implicit none

      !----- Argument. --------------------------------------------------------------------!
      integer, intent(in) :: ierr
      !------------------------------------------------------------------------------------!

      write (unit=*,fmt='(a)')    '           [!] netCDF error reading the variable!'
      write (unit=*,fmt='(a,1x,2a)')                                                       &
                                  '           [!] Message:',trim(nf90_strerror(ierr)),'...'
      return
   end subroutine ncdf_load_err
   !=======================================================================================!
   !=======================================================================================!
#endif
end module mod_ncdf_globio
!==========================================================================================!
!==========================================================================================!
