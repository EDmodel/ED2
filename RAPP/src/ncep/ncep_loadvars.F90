!==========================================================================================!
!==========================================================================================!
!    Subroutine ncep_loadvars.                                                             !
!------------------------------------------------------------------------------------------!
logical function ncep_loadvars(month,year)
   use an_header       , only : info_table    ! ! intent(in)
   use mod_maxdims     , only : maxstr        ! ! intent(in)
   use mod_grid        , only : ssxp          & ! intent(in)
                              , ssyp          & ! intent(in)
                              , sstp          & ! intent(in)
                              , x_1st         & ! intent(in)
                              , y_1st         & ! intent(in)
                              , t_1st         & ! intent(in)
                              , tlast         ! ! intent(in)
   use mod_model       , only : ngrids        ! ! intent(in)
   use mod_ncep        , only : ncep_g        & ! intent(inout)
                              , nvars_ncep    & ! intent(in)
                              , vars_ncep     & ! intent(in)
                              , grids_ncep    & ! intent(in)
                              , prefvars_ncep ! ! intent(in)
   use mod_ioopts      , only : inpath        ! ! intent(in)

   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer              , intent(in) :: month
   integer              , intent(in) :: year
#if USE_NCDF
   !----- Local variables. ----------------------------------------------------------------!
   character(len=maxstr)             :: filename_foye ! File name (following year).
   integer                           :: x,y,t         ! X, Y, and T Counters
   integer                           :: ifm           ! Grid counter
   integer                           :: nv            ! Variable counter
   integer                           :: ierr          ! Error flag (currently not used)
   logical                           :: file_is_there ! Flag for file existence.
   !---------------------------------------------------------------------------------------!
   

   write (unit=*,fmt='(a)') '     - Loading the NCEP variables...'

   !---------------------------------------------------------------------------------------!
   !    Initialise with error flag.  Only if the function reaches the end we assign the    !
   ! success flag.                                                                         !
   !---------------------------------------------------------------------------------------!
   ncep_loadvars = .false.
   
   select case (month)
   case (1:11)
      !------------------------------------------------------------------------------------!
      !    January-November, use the information already available.                        !
      !------------------------------------------------------------------------------------!

      !----- This loop is kind of silly, but it avoids some logical tests. ----------------!
      varloop: do nv=1,nvars_ncep

         !----- Setting some short cuts for useful variables. -----------------------------!
         ifm = grids_ncep(nv)

         !----- This case selection is to send the variable to the right pointer. ---------!
         select case(trim(vars_ncep(nv)))
         case ('air')
            call loadshort_jannov(info_table(nv)%filename,vars_ncep(nv),ssxp(ifm)          &
                                 ,ssyp(ifm),sstp(ifm),x_1st(ifm),y_1st(ifm),t_1st(ifm)     &
                                 ,tlast(ifm),ncep_g(ifm)%temp,ierr)

         case ('pres')
            call loadshort_jannov(info_table(nv)%filename,vars_ncep(nv),ssxp(ifm)          &
                                 ,ssyp(ifm),sstp(ifm),x_1st(ifm),y_1st(ifm),t_1st(ifm)     &
                                 ,tlast(ifm),ncep_g(ifm)%pres,ierr)

         case ('rhum')
            call loadshort_jannov(info_table(nv)%filename,vars_ncep(nv),ssxp(ifm)          &
                                 ,ssyp(ifm),sstp(ifm),x_1st(ifm),y_1st(ifm),t_1st(ifm)     &
                                 ,tlast(ifm),ncep_g(ifm)%rhum,ierr)

         case ('uwnd')
            call loadshort_jannov(info_table(nv)%filename,vars_ncep(nv),ssxp(ifm)          &
                                 ,ssyp(ifm),sstp(ifm),x_1st(ifm),y_1st(ifm),t_1st(ifm)     &
                                 ,tlast(ifm),ncep_g(ifm)%uwnd,ierr)

         case ('vwnd')
            call loadshort_jannov(info_table(nv)%filename,vars_ncep(nv),ssxp(ifm)          &
                                 ,ssyp(ifm),sstp(ifm),x_1st(ifm),y_1st(ifm),t_1st(ifm)     &
                                 ,tlast(ifm),ncep_g(ifm)%vwnd,ierr)

         case ('dlwrf')
            call loadshort_jannov(info_table(nv)%filename,vars_ncep(nv),ssxp(ifm)          &
                                 ,ssyp(ifm),sstp(ifm),x_1st(ifm),y_1st(ifm),t_1st(ifm)     &
                                 ,tlast(ifm),ncep_g(ifm)%dlwrf,ierr)

         case ('nbdsf')
            call loadshort_jannov(info_table(nv)%filename,vars_ncep(nv),ssxp(ifm)          &
                                 ,ssyp(ifm),sstp(ifm),x_1st(ifm),y_1st(ifm),t_1st(ifm)     &
                                 ,tlast(ifm),ncep_g(ifm)%nbdsf,ierr)

         case ('nddsf')
            call loadshort_jannov(info_table(nv)%filename,vars_ncep(nv),ssxp(ifm)          &
                                 ,ssyp(ifm),sstp(ifm),x_1st(ifm),y_1st(ifm),t_1st(ifm)     &
                                 ,tlast(ifm),ncep_g(ifm)%nddsf,ierr)

         case ('vbdsf')
            call loadshort_jannov(info_table(nv)%filename,vars_ncep(nv),ssxp(ifm)          &
                                 ,ssyp(ifm),sstp(ifm),x_1st(ifm),y_1st(ifm),t_1st(ifm)     &
                                 ,tlast(ifm),ncep_g(ifm)%vbdsf,ierr)

         case ('vddsf')
            call loadshort_jannov(info_table(nv)%filename,vars_ncep(nv),ssxp(ifm)          &
                                 ,ssyp(ifm),sstp(ifm),x_1st(ifm),y_1st(ifm),t_1st(ifm)     &
                                 ,tlast(ifm),ncep_g(ifm)%vddsf,ierr)

         case ('prate')
            call loadshort_jannov(info_table(nv)%filename,vars_ncep(nv),ssxp(ifm)          &
                                 ,ssyp(ifm),sstp(ifm),x_1st(ifm),y_1st(ifm),t_1st(ifm)     &
                                 ,tlast(ifm),ncep_g(ifm)%prate,ierr)

         end select
      end do varloop

   case (12)
      !------------------------------------------------------------------------------------!
      !    December is an especial case because the last point is in January the next      !
      ! year. And we need to make sure that we have at least one point in the year after.  !
      !------------------------------------------------------------------------------------!

      !----- This loop is kind of silly, but it avoids some logical tests. ----------------!
      dec_varloop: do nv=1,nvars_ncep

         !----- Setting some short cuts for useful variables. -----------------------------!
         ifm = grids_ncep(nv)

         write(filename_foye,fmt='(4a,i4.4,a)')                                            &
                      trim(inpath),'/',trim(prefvars_ncep(nv)),'.',year+1,'.nc'
         inquire(file=filename_foye,exist=file_is_there)

         !---------------------------------------------------------------------------------!
         !    Check whether the next year is available.  If not, we will exit the function !
         ! and prematurely quit the run...                                                 !
         !---------------------------------------------------------------------------------!
         if (.not. file_is_there) then
            write(unit=*,fmt='(3a)') '    File ',trim(filename_foye),' wasn''t found, and'
            write(unit=*,fmt='(3a)') 'I need one point from it for time interpolation...'
            write(unit=*,fmt='(a,1x,i5,1x,a)')                                             &
                               'Therefore, I won''t be able to process the ',year,'data...'
            return
         end if

         !----- This case selection is to send the variable to the right pointer. ---------!
         select case(trim(vars_ncep(nv)))
         case ('air')
            call loadshort_december(info_table(nv)%filename,filename_foye,vars_ncep(nv)    &
                                   ,ssxp(ifm),ssyp(ifm),sstp(ifm),x_1st(ifm),y_1st(ifm)    &
                                   ,t_1st(ifm),tlast(ifm),ncep_g(ifm)%temp,ierr)

         case ('pres')
            call loadshort_december(info_table(nv)%filename,filename_foye,vars_ncep(nv)    &
                                   ,ssxp(ifm),ssyp(ifm),sstp(ifm),x_1st(ifm),y_1st(ifm)    &
                                   ,t_1st(ifm),tlast(ifm),ncep_g(ifm)%pres,ierr)

         case ('rhum')
            call loadshort_december(info_table(nv)%filename,filename_foye,vars_ncep(nv)    &
                                   ,ssxp(ifm),ssyp(ifm),sstp(ifm),x_1st(ifm),y_1st(ifm)    &
                                   ,t_1st(ifm),tlast(ifm),ncep_g(ifm)%rhum,ierr)

         case ('uwnd')
            call loadshort_december(info_table(nv)%filename,filename_foye,vars_ncep(nv)    &
                                   ,ssxp(ifm),ssyp(ifm),sstp(ifm),x_1st(ifm),y_1st(ifm)    &
                                   ,t_1st(ifm),tlast(ifm),ncep_g(ifm)%uwnd,ierr)

         case ('vwnd')
            call loadshort_december(info_table(nv)%filename,filename_foye,vars_ncep(nv)    &
                                   ,ssxp(ifm),ssyp(ifm),sstp(ifm),x_1st(ifm),y_1st(ifm)    &
                                   ,t_1st(ifm),tlast(ifm),ncep_g(ifm)%vwnd,ierr)

         case ('dlwrf')
            call loadshort_december(info_table(nv)%filename,filename_foye,vars_ncep(nv)    &
                                   ,ssxp(ifm),ssyp(ifm),sstp(ifm),x_1st(ifm),y_1st(ifm)    &
                                   ,t_1st(ifm),tlast(ifm),ncep_g(ifm)%dlwrf,ierr)

         case ('nbdsf')
            call loadshort_december(info_table(nv)%filename,filename_foye,vars_ncep(nv)    &
                                   ,ssxp(ifm),ssyp(ifm),sstp(ifm),x_1st(ifm),y_1st(ifm)    &
                                   ,t_1st(ifm),tlast(ifm),ncep_g(ifm)%nbdsf,ierr)

         case ('nddsf')
            call loadshort_december(info_table(nv)%filename,filename_foye,vars_ncep(nv)    &
                                   ,ssxp(ifm),ssyp(ifm),sstp(ifm),x_1st(ifm),y_1st(ifm)    &
                                   ,t_1st(ifm),tlast(ifm),ncep_g(ifm)%nddsf,ierr)

         case ('vbdsf')
            call loadshort_december(info_table(nv)%filename,filename_foye,vars_ncep(nv)    &
                                   ,ssxp(ifm),ssyp(ifm),sstp(ifm),x_1st(ifm),y_1st(ifm)    &
                                   ,t_1st(ifm),tlast(ifm),ncep_g(ifm)%vbdsf,ierr)

         case ('vddsf')
            call loadshort_december(info_table(nv)%filename,filename_foye,vars_ncep(nv)    &
                                   ,ssxp(ifm),ssyp(ifm),sstp(ifm),x_1st(ifm),y_1st(ifm)    &
                                   ,t_1st(ifm),tlast(ifm),ncep_g(ifm)%vddsf,ierr)

         case ('prate')
            call loadshort_december(info_table(nv)%filename,filename_foye,vars_ncep(nv)    &
                                   ,ssxp(ifm),ssyp(ifm),sstp(ifm),x_1st(ifm),y_1st(ifm)    &
                                   ,t_1st(ifm),tlast(ifm),ncep_g(ifm)%prate,ierr)

         end select
      end do dec_varloop
   end select

   !---------------------------------------------------------------------------------------!
   !     Here we fill the ice-liquid potential temperature and specific humidity arrays.   !
   ! We will interpolate these variables rather than temperature and relative humidity.    !
   !---------------------------------------------------------------------------------------!
   ifm = 3 ! The grid with state variables that will be interpolated.
   do t=1,sstp(ifm)
      do y=1,ssyp(ifm)
         do x=1,ssxp(ifm)
         end do
      end do
   end do

   !----- If the function reached this point, switch the flag to success. -----------------!
   ncep_loadvars = .true. 

#else

   !----- If netCDF is not available, no chance for success... ----------------------------!
   ncep_loadvars = .false.

   call fatal_error('You can''t use ncep input without compiling with netcdf libraries...' &
                   ,'ncep_fill_infotable','ncep_fill_infotable.F90')
#endif
end function ncep_loadvars
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This function will load the dataset into the variable array, provided that the month  !
! is not December.                                                                         !
!------------------------------------------------------------------------------------------!
subroutine loadshort_jannov(filename,varname,mxp,myp,mtp,xa,ya,ta,tz,var,ierr)
   use mod_maxdims     , only : maxstr
   use mod_ioopts      , only : missflg_real
#if USE_NCDF
   use netcdf
   use mod_netcdf      , only : ncid          ! ! intent(in)
   use mod_ncdf_globio , only : ncio_2dvar    & ! function
                              , ncio_3dvar    & ! function
                              , ncio_2dshort  & ! function
                              , ncio_3dshort  ! ! function
#endif
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   character(len=maxstr)          , intent(in)    :: filename
   character(len=maxstr)          , intent(in)    :: varname
   integer                        , intent(in)    :: mxp
   integer                        , intent(in)    :: myp
   integer                        , intent(in)    :: mtp
   integer                        , intent(in)    :: xa
   integer                        , intent(in)    :: ya
   integer                        , intent(in)    :: ta
   integer                        , intent(in)    :: tz
   real   , dimension(mxp,myp,mtp), intent(inout) :: var
   integer                        , intent(out)   :: ierr
   !----- Local variables. ----------------------------------------------------------------!
   integer                                        :: xoff
   integer                                        :: yoff
   integer                                        :: tabs
   integer                                        :: tloc
   integer                                        :: misscnt
   real                                           :: varmin
   real                                           :: varmax
   !---------------------------------------------------------------------------------------!

#if USE_NCDF

   !----- Initialise some short cut variables. --------------------------------------------!
   xoff = xa -1
   yoff = ya -1

   !----- Open the file. ------------------------------------------------------------------!
   ierr = nf90_open(filename,NF90_NOWRITE,ncid)

   !----- Loop over times and load the variable into the array. ---------------------------!
   do tabs=ta,tz
      tloc = tabs - ta + 1
      ierr = ncio_2dshort(varname,.true.,tabs,mxp,myp,var(:,:,tloc)                        &
                         ,offin1=xoff,offin2=yoff)
   end do

   !----- Closing the file ----------------------------------------------------------------!
   ierr = nf90_close(ncid)

   !----- Quick statistics to entretain the user and warn possible problems. --------------!
   varmin  = minval(var,mask=var /= missflg_real)
   varmax  = maxval(var,mask=var /= missflg_real)
   misscnt = count(var == missflg_real)
   write (unit=*,fmt='(2(a,1x),2(a,1x,es14.7,1x),a,1x,i6,a)')                              &
      '         [|] Retrieved :',trim(varname),'. Range: [',varmin,':',varmax              &
                    ,'] . # missing: ',misscnt,'.'

#else
   !----- Putting anything, this run is about to crash anyway... --------------------------!
   ierr = -4321

   call fatal_error('You can''t use ncep input without compiling with netcdf libraries...' &
                   ,'loadshort_jannov','ncep_loadvars.F90')
#endif

end subroutine loadshort_jannov
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This function will load the dataset into the variable array for December.  This is an !
! especial case because we need one extra point, which is the first data in the following  !
! year file.                                                                               !
!------------------------------------------------------------------------------------------!
subroutine loadshort_december(filename,filename_foye,varname,mxp,myp,mtp,xa,ya,ta,tz       &
                             ,var,ierr)
   use mod_maxdims     , only : maxstr
   use mod_ioopts      , only : missflg_real
#if USE_NCDF
   use netcdf
   use mod_netcdf      , only : ncid          ! ! intent(in)
   use mod_ncdf_globio , only : ncio_2dvar    & ! function
                              , ncio_3dvar    & ! function
                              , ncio_2dshort  & ! function
                              , ncio_3dshort  ! ! function
#endif
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   character(len=maxstr)          , intent(in)    :: filename
   character(len=maxstr)          , intent(in)    :: filename_foye
   character(len=maxstr)          , intent(in)    :: varname
   integer                        , intent(in)    :: mxp
   integer                        , intent(in)    :: myp
   integer                        , intent(in)    :: mtp
   integer                        , intent(in)    :: xa
   integer                        , intent(in)    :: ya
   integer                        , intent(in)    :: ta
   integer                        , intent(in)    :: tz
   real   , dimension(mxp,myp,mtp), intent(inout) :: var
   integer                        , intent(out)   :: ierr
   !----- Local variables. ----------------------------------------------------------------!
   integer                                        :: xoff
   integer                                        :: yoff
   integer                                        :: tzm1
   integer                                        :: tabs
   integer                                        :: tloc
   integer                                        :: misscnt
   real                                           :: varmin
   real                                           :: varmax
   !---------------------------------------------------------------------------------------!

#if USE_NCDF

   !----- Initialise some short cut variables. --------------------------------------------!
   tzm1 = tz -1
   xoff = xa -1
   yoff = ya -1


   !----- Open this year file. ------------------------------------------------------------!
   ierr = nf90_open(filename,NF90_NOWRITE,ncid)
   !----- Loop over times (except the last time) and load the variable into the array. ----!
   do tabs=ta,tzm1
      tloc = tabs - ta + 1
      write(unit=*,fmt='(4(a,1x,i5,1x))') 'TA=',ta,'TZ=',tz,'TABS=',tabs,'TLOC=',tloc
      ierr = ncio_2dshort(varname,.true.,tabs,mxp,myp,var(:,:,tloc)                        &
                         ,offin1=xoff,offin2=yoff)
   end do
   !----- Closing the file ----------------------------------------------------------------!
   ierr = nf90_close(ncid)
   !---------------------------------------------------------------------------------------!



   !----- Open the following year file. ---------------------------------------------------!
   ierr = nf90_open(filename_foye,NF90_NOWRITE,ncid)
   !----- Load the last time (first time of following year). ------------------------------!
   ierr = ncio_2dshort(varname,.true.,1,mxp,myp,var(:,:,tz),offin1=xa,offin2=ya)
   !----- Closing the file ----------------------------------------------------------------!
   ierr = nf90_close(ncid)
   !---------------------------------------------------------------------------------------!

   !----- Quick statistics to entretain the user and warn possible problems. --------------!
   varmin  = minval(var,mask=var /= missflg_real)
   varmax  = maxval(var,mask=var /= missflg_real)
   misscnt = count(var == missflg_real)
   write (unit=*,fmt='(2(a,1x),2(a,1x,es14.7,1x),a,1x,i6,a)')                              &
      '         [|] Retrieved :',trim(varname),'. Range: [',varmin,':',varmax              &
                    ,'] . # missing: ',misscnt,'.'


#else
   !----- Putting anything, this run is about to crash anyway... --------------------------!
   ierr = -4321

   call fatal_error('You can''t use ncep input without compiling with netcdf libraries...' &
                   ,'loadshort_december','ncep_loadvars.F90')
#endif

end subroutine loadshort_december
!==========================================================================================!
!==========================================================================================!
