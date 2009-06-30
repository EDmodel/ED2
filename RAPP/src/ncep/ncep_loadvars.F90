!==========================================================================================!
!==========================================================================================!
!    Subroutine ncep_loadvars.                                                             !
!------------------------------------------------------------------------------------------!
logical function ncep_loadvars()
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
   use mod_ioopts      , only : inpath        & ! intent(in)
                              , missflg_real  ! ! intent(in)
   use therm_lib       , only : ptrh2rvapil   & ! function
                              , theta_iceliq  ! ! function

   implicit none
#if USE_NCDF
   !----- Local variables. ----------------------------------------------------------------!
   integer                           :: x,y,t         ! X, Y, and T Counters
   integer                           :: ifm           ! Grid counter
   integer                           :: nv            ! Variable counter
   integer                           :: ierr          ! Error flag (currently not used)
   real                              :: mixr          ! Mixing ratio               [ kg/kg]
   !---------------------------------------------------------------------------------------!
   

   write (unit=*,fmt='(a)') '     - Loading the NCEP variables...'

   !---------------------------------------------------------------------------------------!
   !    Initialise with error flag.  Only if the function reaches the end we assign the    !
   ! success flag.                                                                         !
   !---------------------------------------------------------------------------------------!
   ncep_loadvars = .false.


   !----- This loop is kind of silly, but it avoids some logical tests. ----------------!
   varloop: do nv=1,nvars_ncep

      !----- Setting some short cuts for useful variables. --------------------------------!
      ifm = grids_ncep(nv)

      !----- This case selection is to send the variable to the right pointer. ------------!
      select case(trim(vars_ncep(nv)))
      case ('air')
         call loadshort(info_table(nv)%filename,vars_ncep(nv),ssxp(ifm),ssyp(ifm)          &
                       ,sstp(ifm),x_1st(ifm),y_1st(ifm),t_1st(ifm),tlast(ifm)              &
                       ,ncep_g(ifm)%temp,ierr)

      case ('pres')
         call loadshort(info_table(nv)%filename,vars_ncep(nv),ssxp(ifm),ssyp(ifm)          &
                       ,sstp(ifm),x_1st(ifm),y_1st(ifm),t_1st(ifm),tlast(ifm)              &
                       ,ncep_g(ifm)%pres,ierr)

      case ('rhum')
         call loadshort(info_table(nv)%filename,vars_ncep(nv),ssxp(ifm),ssyp(ifm)          &
                       ,sstp(ifm),x_1st(ifm),y_1st(ifm),t_1st(ifm),tlast(ifm)              &
                       ,ncep_g(ifm)%rhum,ierr)

      case ('uwnd')
         call loadshort(info_table(nv)%filename,vars_ncep(nv),ssxp(ifm),ssyp(ifm)          &
                       ,sstp(ifm),x_1st(ifm),y_1st(ifm),t_1st(ifm),tlast(ifm)              &
                       ,ncep_g(ifm)%uwnd,ierr)

      case ('vwnd')
         call loadshort(info_table(nv)%filename,vars_ncep(nv),ssxp(ifm),ssyp(ifm)          &
                       ,sstp(ifm),x_1st(ifm),y_1st(ifm),t_1st(ifm),tlast(ifm)              &
                       ,ncep_g(ifm)%vwnd,ierr)

      case ('dlwrf')
         call loadshort(info_table(nv)%filename,vars_ncep(nv),ssxp(ifm),ssyp(ifm)          &
                       ,sstp(ifm),x_1st(ifm),y_1st(ifm),t_1st(ifm),tlast(ifm)              &
                       ,ncep_g(ifm)%dlwrf,ierr)

      case ('nbdsf')
         call loadshort(info_table(nv)%filename,vars_ncep(nv),ssxp(ifm),ssyp(ifm)          &
                       ,sstp(ifm),x_1st(ifm),y_1st(ifm),t_1st(ifm),tlast(ifm)              &
                       ,ncep_g(ifm)%nbdsf,ierr)

      case ('nddsf')
         call loadshort(info_table(nv)%filename,vars_ncep(nv),ssxp(ifm),ssyp(ifm)          &
                       ,sstp(ifm),x_1st(ifm),y_1st(ifm),t_1st(ifm),tlast(ifm)              &
                       ,ncep_g(ifm)%nddsf,ierr)

      case ('vbdsf')
         call loadshort(info_table(nv)%filename,vars_ncep(nv),ssxp(ifm),ssyp(ifm)          &
                       ,sstp(ifm),x_1st(ifm),y_1st(ifm),t_1st(ifm),tlast(ifm)              &
                       ,ncep_g(ifm)%vbdsf,ierr)

      case ('vddsf')
         call loadshort(info_table(nv)%filename,vars_ncep(nv),ssxp(ifm),ssyp(ifm)          &
                       ,sstp(ifm),x_1st(ifm),y_1st(ifm),t_1st(ifm),tlast(ifm)              &
                       ,ncep_g(ifm)%vddsf,ierr)

      case ('prate')
         call loadshort(info_table(nv)%filename,vars_ncep(nv),ssxp(ifm),ssyp(ifm)          &
                       ,sstp(ifm),x_1st(ifm),y_1st(ifm),t_1st(ifm),tlast(ifm)              &
                       ,ncep_g(ifm)%prate,ierr)

      end select
   end do varloop

   !---------------------------------------------------------------------------------------!
   !     Here we will find the specific humidity that will be written in the output        !
   ! instead of relative humidity.                                                         !
   ! 1. Relative humidity will be converted to fraction instead of percentage;             !
   ! 2. Find mixing ratio using the thermodynamic library;                                 !
   ! 3. Convert mixing ratio to specific humidity;                                         !
   !---------------------------------------------------------------------------------------!
   ifm = 3 ! The grid with state variables that will be interpolated.
   do t=1,sstp(ifm)
      do y=1,ssyp(ifm)
         do x=1,ssxp(ifm)
            !------------------------------------------------------------------------------!
            !     Converting relative humidity to fraction.  Forcing it to never be super- !
            ! saturated.                                                                   !
            !------------------------------------------------------------------------------!
            if (ncep_g(ifm)%rhum(x,y,t) /= missflg_real) then
               ncep_g(ifm)%rhum(x,y,t) = min(1.,max(0.,ncep_g(ifm)%rhum(x,y,t) * 0.01))
            end if
            
            !----- Finding mixing ratio and then specific humidity. -----------------------!
            mixr = ptrh2rvapil(ncep_g(ifm)%rhum(x,y,t),ncep_g(ifm)%pres(x,y,t)             &
                              ,ncep_g(ifm)%temp(x,y,t))
            ncep_g(ifm)%shum(x,y,t) = mixr / (1. + mixr)

         end do
      end do
   end do

   !----- Here we ensure that radiation and precipitation are never negative. -------------!
   where (ncep_g(1)%prate < 0.) ncep_g(1)%prate = 0.
   where (ncep_g(1)%dlwrf < 0.) ncep_g(1)%dlwrf = 0.
   where (ncep_g(1)%nbdsf < 0.) ncep_g(1)%nbdsf = 0.
   where (ncep_g(1)%nddsf < 0.) ncep_g(1)%nddsf = 0.
   where (ncep_g(1)%vbdsf < 0.) ncep_g(1)%vbdsf = 0.
   where (ncep_g(1)%vddsf < 0.) ncep_g(1)%vddsf = 0.

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
!    This function will load the dataset into the variable array.                          !
!------------------------------------------------------------------------------------------!
subroutine loadshort(filename,varname,mxp,myp,mtp,xa,ya,ta,tz,var,ierr)
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

end subroutine loadshort
!==========================================================================================!
!==========================================================================================!
