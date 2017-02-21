!==========================================================================================!
!    These subroutines are intended to perform all necessary checks on the parameters      !
! given by ED2IN. If there is any incompability, then the model shouldn't even start, and  !
! an error message returned to the standard outputs so the user knows what to fix.         !
!    Avoid leaving such set up checks in the model, it will be very likely to be inside    !
! the time loop and that means that the check will be performed many times...              !
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!   This subroutine checks the grid definition for ED, like grid size, dimensions and if   !
! the provided numbers are allowed.                                                        !
!------------------------------------------------------------------------------------------!
subroutine ed_opspec_grid
   use ed_max_dims , only : maxgrds             & ! intent(in)
                          , nxpmax              & ! intent(in)
                          , nypmax              & ! intent(in)
                          , nzpmax              & ! intent(in)
                          , nzgmax              & ! intent(in)
                          , nzsmax              & ! intent(in)
                          , max_poi             & ! intent(in)
                          , max_ed_regions      ! ! intent(in)
   use grid_coms   , only : nnxp                & ! intent(in)
                          , nnyp                & ! intent(in)
                          , ngrids              & ! intent(in)
                          , polelat             & ! intent(in)
                          , polelon             & ! intent(in)
                          , centlat             & ! intent(in)
                          , centlon             & ! intent(in)
                          , deltax              & ! intent(in)
                          , deltay              & ! intent(in)
                          , nstratx             & ! intent(in)
                          , nstraty             & ! intent(in)
                          , nzg                 & ! intent(in)
                          , nzs                 ! ! intent(in)
   use mem_polygons, only : n_ed_region         & ! intent(in)
                          , n_poi               & ! intent(in)
                          , grid_type           & ! intent(in)
                          , grid_res            & ! intent(in)
                          , poi_lat             & ! intent(in)
                          , poi_lon             & ! intent(in)
                          , poi_res             & ! intent(in)
                          , ed_reg_latmin       & ! intent(in)
                          , ed_reg_latmax       & ! intent(in)
                          , ed_reg_lonmin       & ! intent(in)
                          , ed_reg_lonmax       ! ! intent(in)
   use soil_coms   , only : slz                 ! ! intent(in)

   implicit none
   !----- Local variables. ----------------------------------------------------------------!
   integer            :: ifaterr
   integer            :: ifm
   integer            :: k
   character(len=222) :: reason
   !---------------------------------------------------------------------------------------!



   !----- Start assuming everything will be fine. -----------------------------------------!
   ifaterr = 0
   !---------------------------------------------------------------------------------------!


   !----- Check whether the number of grids is reasonable. --------------------------------!
   if (ngrids < 1) then
      write(reason,'(a,1x,i4,a)') &
           'Number of grids needs to be positive. (Yours is ',ngrids,').'
      call opspec_fatal(reason,'opspec_grid')
      ifaterr=ifaterr+1
   elseif (ngrids > maxgrds) then
      write(reason,'(2(a,1x,i4,a))')  &
           'Number of grids cannot be larger than ',maxgrds,'.',' (Yours is ',ngrids,').'
      call opspec_fatal(reason,'opspec_grid')  
      ifaterr=ifaterr+1
   end if
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Check whether we are trying to run POI and regional simultaneously.  This is     !
   ! currently forbidden.                                                                  !
   !---------------------------------------------------------------------------------------!
   if ((.not. (n_poi > 0 .and. n_ed_region == 0)) .and.                                    &
       (.not. (n_poi == 0 .and. n_ed_region > 0))       ) then
      write(reason,'(a,1x,i4,1x,a,1x,i4,1x,a)')                                            &
            'One of n_poi or n_ed_region needs to be zero.( Yours: n_poi='                 &
           ,n_poi,', n_ed_region=',n_ed_region,').'
      call opspec_fatal(reason,'opspec_grid')  
      ifaterr=ifaterr+1

   elseif (n_poi < 0) then
      write(reason,'(a,1x,i4,a)')                                                          &
         'N_POI needs to be non-negative. Yours is currently set to ',n_poi,'.'
      call opspec_fatal(reason,'opspec_grid')  
      ifaterr=ifaterr+1

   elseif (n_ed_region < 0) then
      write(reason,'(a,1x,i4,a)')                                                          &
         'N_ED_REGION needs to be non-negative. Yours is currently set to ',n_ed_region,'.'
      call opspec_fatal(reason,'opspec_grid')  
      ifaterr=ifaterr+1

   elseif (n_poi > max_poi) then
      write(reason,'(a,1x,i4,a,1x,i4,a)')                                                  &
         'N_POI cannot be greater than',max_poi,'. Yours is currently set to ',n_poi,'.'
      call opspec_fatal(reason,'opspec_grid')  
      ifaterr=ifaterr+1

   elseif (n_ed_region> max_ed_regions) then
      write(reason,'(a,1x,i4,a,1x,i4,a)')                                                  &
          'N_ED_REGION cannot be greater than',max_ed_regions                              &
         ,'. Yours is currently set to ',n_ed_region,'.'
      call opspec_fatal(reason,'opspec_grid')  
      ifaterr=ifaterr+1
   end if
   !---------------------------------------------------------------------------------------!




   !----- Check whether the number of grid points makes sense for ED regions. -------------!
   do ifm=1,n_ed_region
      if (nnxp(ifm) > nxpmax) then
         write(reason,'(3(a,1x,i4,a))')                                                    &
              'Number of x-points cannot be greater than ',nxpmax,'.',' (Yours is'         &
             ,nnxp(ifm),' ','in grid',ifm,').'
         call opspec_fatal(reason,'opspec_grid')  
         ifaterr=ifaterr+1
      end if
      if (nnyp(ifm) > nypmax) then
         write(reason,'(3(a,1x,i4,a))')                                                    &
              'Number of y-points cannot be greater than ',nypmax,'.',' (Yours is'         &
             ,nnyp(ifm),' ','in grid',ifm,').'
         call opspec_fatal(reason,'opspec_grid')  
         ifaterr=ifaterr+1
      end if
      if (nstratx(ifm) < 1) then
         write (reason,'(a,1x,i4,1x,a,1x,i4,a)')                                           &
               'Nest x ratio should be at least 1. Your nstratx for grid',ifm              &
              ,'is currently set to',nstratx(ifm),'.'
         call opspec_fatal(reason,'opspec_grid')  
         ifaterr=ifaterr+1        
      end if
      if (nstraty(ifm) < 1) then
         write (reason,'(a,1x,i4,1x,a,1x,i4,a)')                                           &
               'Nest y ratio should be at least 1. Your nstraty for grid',ifm              &
              ,'is currently set to',nstratx(ifm),'.'
         call opspec_fatal(reason,'opspec_grid')  
         ifaterr=ifaterr+1        
      end if
   end do
   !---------------------------------------------------------------------------------------!



   !-----  Check whether ED has reasonable values to define the grids. --------------------!
   if (n_ed_region > 0) then
      if (grid_type < 0 .or. grid_type > 1) then
         write (reason,'(a,2x,a,1x,i4,a)')                                                 &
               'Grid_type should be either 0 or 1 for regional runs.'                      &
              ,'Yours is currently set to',grid_type,'.'
         call opspec_fatal(reason,'opspec_grid')  
         ifaterr=ifaterr+1
      end if
      if (grid_type == 0) then
         do ifm=1,n_ed_region
            if (ed_reg_latmin(ifm) < -90. .or. ed_reg_latmin(ifm) > 90. ) then
               write (reason,'(a,1x,i4,a,1x,es14.7,a)')                                    &
                     'ED_REG_LATMIN is outside [-90.;90.] for region ',ifm                 &
                    ,'. Yours is currently set to',ed_reg_latmin(ifm),'...'
               call opspec_fatal(reason,'opspec_grid')  
               ifaterr=ifaterr+1
            end if
            if (ed_reg_latmax(ifm) < -90. .or. ed_reg_latmax(ifm) > 90. ) then
               write (reason,'(a,1x,i4,a,1x,es14.7,a)')                                    &
                     'ED_REG_LATMAX is outside [-90.;90.] for region ',ifm                 &
                    ,'. Yours is currently set to',ed_reg_latmax(ifm),'...'
               call opspec_fatal(reason,'opspec_grid')  
               ifaterr=ifaterr+1
            end if
            if (ed_reg_lonmin(ifm) < -180. .or. ed_reg_lonmin(ifm) > 180. ) then
               write (reason,'(a,1x,i4,a,1x,es14.7,a)')                                    &
                     'ED_REG_LONMIN is outside [-180.;180.] for region ',ifm               &
                    ,'. Yours is currently set to',ed_reg_latmin(ifm),'...'
               call opspec_fatal(reason,'opspec_grid')  
               ifaterr=ifaterr+1
            end if
            if (ed_reg_lonmax(ifm) < -180. .or. ed_reg_lonmax(ifm) > 180. ) then
               write (reason,'(a,1x,i4,a,1x,es14.7,a)')                                    &
                     'ED_REG_LONMAX is outside [-180.;180.] for region ',ifm               &
                    ,'. Yours is currently set to',ed_reg_latmax(ifm),'...'
               call opspec_fatal(reason,'opspec_grid')  
               ifaterr=ifaterr+1
            end if
         end do
         if (grid_res <= 0.) then
            write (reason,'(a,1x,es14.7,a)')                                               &
                  'GRID_RES must be positive. Yours is currently set to',grid_res,'...'
            call opspec_fatal(reason,'opspec_grid')  
            ifaterr=ifaterr+1
         end if
      elseif (grid_type == 1) then
         if (deltax <= 0.) then
            write (reason,'(a,1x,es14.7,a)')                                               &
                  'DELTAX must be positive. Yours is currently set to',deltax,'...'
            call opspec_fatal(reason,'opspec_grid')  
            ifaterr=ifaterr+1
         end if
         if (deltay <= 0.) then
            write (reason,'(a,1x,es14.7,a)')                                               &
                  'DELTAY must be positive. Yours is currently set to',deltax,'...'
            call opspec_fatal(reason,'opspec_grid')  
            ifaterr=ifaterr+1
         end if
         if (polelat < -90. .or. polelat > 90.) then
            write (reason,'(a,1x,i4,a)')                                                   &
                  'POLELAT is outside [-90.;90.].',' Yours is currently set to'            &
                 ,polelat,'...'
            call opspec_fatal(reason,'opspec_grid')  
            ifaterr=ifaterr+1
         end if
         if (polelon < -180. .or. polelon > 180.) then
            write (reason,'(a,1x,i4,a)')                                                   &
                  'POLELON is outside [-180.;180.].',' Yours is currently set to'          &
                 ,polelon,'...'
            call opspec_fatal(reason,'opspec_grid')  
            ifaterr=ifaterr+1
         end if
         do ifm=1,n_ed_region
            if (centlat(ifm) < -90. .or. centlat(ifm) > 90. ) then
               write (reason,'(a,1x,i4,a,1x,es14.7,a)')                                    &
                     'CENTLAT is outside [-90.;90.] for region ',ifm                       &
                    ,'. Yours is currently set to',centlat(ifm),'...'
               call opspec_fatal(reason,'opspec_grid')  
               ifaterr=ifaterr+1
            end if
            if (centlon(ifm) < -180. .or. centlon(ifm) > 180. ) then
               write (reason,'(a,1x,i4,a,1x,es14.7,a)')                                    &
                     'CENTLON is outside [-180.;180.] for region ',ifm                     &
                    ,'. Yours is currently set to',centlon(ifm),'...'
               call opspec_fatal(reason,'opspec_grid')  
               ifaterr=ifaterr+1
            end if
         end do
      end if
   end if
   !---------------------------------------------------------------------------------------!




   !----- Check whether the grid specifications for POI runs are good. --------------------!
   if (n_poi > 0)then
      do ifm=1,n_poi
         if (poi_lat(ifm) < -90. .or. poi_lat(ifm) > 90. ) then
            write (reason,'(a,1x,i4,a,1x,es14.7,a)')                                       &
                  'POI_LAT is outside [-90.;90.] for POI #',ifm                            &
                 ,'. Yours is currently set to',poi_lat(ifm),'...'
            call opspec_fatal(reason,'opspec_grid')  
            ifaterr=ifaterr+1
         end if
         if (poi_lon(ifm) < -180. .or. poi_lon(ifm) > 180. ) then
            write (reason,'(a,1x,i4,a,1x,es14.7,a)')                                       &
                  'POI_LON is outside [-180.;180.] for POI #',ifm                          &
                  ,'. Yours is currently set to',poi_lon(ifm),'...'
            call opspec_fatal(reason,'opspec_grid')  
            ifaterr=ifaterr+1
         end if
         if (poi_res(ifm) < 0.001 .or. poi_res(ifm) > 5.) then
            write (reason,'(a,1x,i4,a,1x,es14.7,a)')                                       &
                  'POI_RES is outside [0.001;5.] range for POI #'                          &
                 ,ifm,'. Yours is currently set to',poi_res(ifm),'...'
            call opspec_fatal(reason,'opspec_grid')  
            ifaterr=ifaterr+1
         end if
      end do

   end if

   !---------------------------------------------------------------------------------------!
   !      Check whether ED soil layers are reasonable, i.e, enough layers, sorted from the !
   ! deepest to the shallowest.                                                            !
   !---------------------------------------------------------------------------------------!
   if (nzg < 2) then
      write (reason,'(a,1x,i4,a)')                                                         &
            'Too few soil layers.  Set it to at least 2. Your nzg is currently set to'     &
           ,nzg,'...'
      call opspec_fatal(reason,'opspec_grid')  
      ifaterr=ifaterr+1        
   elseif (nzg > nzgmax) then 
      write (reason,'(2(a,1x,i5,a))')                                                      &
            'The number of soil layers cannot be greater than ',nzgmax,'.'                 &
           ,' Your nzg is currently set to',nzg,'.'
      call opspec_fatal(reason,'opspec_grid') 
      ifaterr=ifaterr+1 
   end if
   do k=1,nzg
      if (slz(k) > -.001) then
         write (reason,'(a,1x,i4,1x,a,1x,es14.7,a)')                                       &
               'Your soil level #',k,'is not enough below ground. It is currently set to'  &
               ,slz(k),', make it deeper than -0.001...'
         call opspec_fatal(reason,'opspec_grid')  
         ifaterr=ifaterr+1        
      end if
   end do

   do k=1,nzg-1
      if (slz(k)-slz(k+1) > .001) then
         write (reason,'(2(a,1x,i4,1x),a,2x,a,1x,es14.7,1x,a,1x,es14.7,a)')                &
               'Soil layers #',k,'and',k+1,'are not enough apart (i.e. > 0.001).'          &
              ,'They are currently set as ',slz(k),'and',slz(k+1),'...'
         call opspec_fatal(reason,'opspec_grid')  
         ifaterr=ifaterr+1        
      end if
   end do


   !---------------------------------------------------------------------------------------!
   !     Check whether ED snow layers are well set, i.e., the number of soil levels is     !
   ! within the allowed range.                                                             !
   !---------------------------------------------------------------------------------------!
   if (nzs < 1) then
      write (reason,'(a,2x,a,1x,i4,a)')                                                    &
            'Too few maximum # of snow layers. Set it to at least 1.'                      &
           ,'Your nzs is currently set to',nzs,'.'
      call opspec_fatal(reason,'opspec_grid')  
      ifaterr=ifaterr+1        
   elseif (nzs > nzsmax) then 
      write (reason,'(2(a,1x,i5,a))')                                                      &
            'The number of snow layers cannot be greater than ',nzsmax,'.'                 &
           ,' Your nzs is currently set to',nzs,'.'
      call opspec_fatal(reason,'opspec_grid') 
      ifaterr=ifaterr+1 
   end if
   !---------------------------------------------------------------------------------------!


   !----- Stop the run in case there is any fatal error. ----------------------------------!
   if (ifaterr > 0) then
      write (unit=*,fmt='(a)')       ' -----------ED_OPSPEC_GRID -------------------------'
      write (unit=*,fmt='(a,1x,i5)') ' Fatal errors:',ifaterr
      write (unit=*,fmt='(a)')       ' ---------------------------------------------------'
      call fatal_error('Fatal errors at namelist - Grid settings '                         &
                      ,'ed_opspec_grid','ed_opspec.f90')
   end if
   !---------------------------------------------------------------------------------------!
   return
end subroutine ed_opspec_grid
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine ed_opspec_par
!------------------------------------------------------------------------------------------!
!   This subroutine checks the grid definition for ED, like grid size, dimensions and if   !
! the provided numbers are allowed.                                                        !
!------------------------------------------------------------------------------------------!
   use ed_max_dims  , only : maxmach
   use grid_coms    , only : nnxp,nnyp,ngrids
   use mem_polygons , only : n_ed_region,n_poi
   use ed_para_coms , only : nmachs
   implicit none
   integer :: ifaterr,ifm
   character(len=222)           :: reason
   ifaterr=0

   ! Check whether the number of slaves is more than what is allowed.
   if (nmachs > maxmach-1) then
      write(reason,'(2(a,1x,i4,a))')  &
           'Number of nodes cannot be larger than ',maxmach,'.',' (Yours is ',nmachs+1,').'
      call opspec_fatal(reason,'opspec_par')  
      ifaterr=ifaterr+1
   end if
   do ifm=1,ngrids
      if (n_ed_region > 0 .and. nnxp(ifm)*nnyp(ifm) < (nmachs +1)) then
         write (reason,'(3(a,1x,i5,1x),a)')                                                 &
               'Region number',ifm,'is too small. You have only',nnxp(ifm)*nnyp(ifm)        &
              ,'points for',nmachs+1,'nodes.'
         call opspec_fatal(reason,'opspec_par')  
         ifaterr=ifaterr+1
      end if
   end do
   ! Checking if the user is trying to run POI in parallel 
   ! (this will be allowed at some point).
   if (n_poi > 0 .and. nmachs > 0) then
      write (reason,'(3(a,1x,i5,1x),a)')                                                    &
        'You are attempting to run POI runs in parallel, this is not supported currently...'
      call opspec_fatal(reason,'opspec_par')  
      ifaterr=ifaterr+1
   end if

   ! stop the run if there are any fatal errors.
   if (ifaterr > 0) then
      write (unit=*,fmt='(a)')       ' -----------ED_OPSPEC_GRID --------------------------'
      write (unit=*,fmt='(a,1x,i5)') ' Fatal errors:',ifaterr
      write (unit=*,fmt='(a)')       ' ----------------------------------------------------'
      call fatal_error('Fatal errors at namelist - Grid settings '&
                     & ,'ed_opspec_grid','ed_opspec.f90')
   end if

   return
end subroutine ed_opspec_par
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine checks time related settings from ED2IN.                              !
!------------------------------------------------------------------------------------------!
subroutine ed_opspec_times
   use ed_misc_coms , only : frqfast          & ! intent(in)
                           , frqstate         & ! intent(in)
                           , frqsum           & ! intent(in)
                           , imontha          & ! intent(in)
                           , idatea           & ! intent(in)
                           , iyeara           & ! intent(in)
                           , itimea           & ! intent(in)
                           , imonthz          & ! intent(in)
                           , idatez           & ! intent(in)
                           , iyearz           & ! intent(in)
                           , itimez           & ! intent(in)
                           , dtlsm            & ! intent(in)
                           , radfrq           & ! intent(in)
                           , ifoutput         & ! intent(in)
                           , isoutput         & ! intent(in)
                           , idoutput         & ! intent(in)
                           , imoutput         & ! intent(in)
                           , iyoutput         & ! intent(in)
                           , iqoutput         & ! intent(in)
                           , itoutput         & ! intent(in)
                           , nrec_fast        & ! intent(in)
                           , nrec_state       & ! intent(in)
                           , outfast          & ! intent(in)
                           , outstate         & ! intent(in)
                           , unitfast         & ! intent(in)
                           , unitstate        ! ! intent(in)
   use consts_coms  , only : day_sec          & ! intent(in)
                           , hr_sec           ! ! intent(in)
   use grid_coms    , only : timmax           ! ! intent(in)
   use ed_misc_coms , only : fast_diagnostics ! ! intent(in)
   use ed_max_dims  , only : str_len          ! ! intent(in)

   implicit none
   !------ Local variables. ---------------------------------------------------------------!
   character(len=str_len) :: reason
   integer                :: ifaterr
   !---------------------------------------------------------------------------------------!
   ifaterr=0

   !---------------------------------------------------------------------------------------!
   !     Check whether the user provided a valid combination of unitfast, frqfast, and     !
   ! outfast.                                                                              !
   !---------------------------------------------------------------------------------------!
   select case(unitfast)
   !---------------------------------------------------------------------------------------!
   !    Seconds                                                                            !
   !---------------------------------------------------------------------------------------!
   case (0)
      if (ifoutput == 0 .and. iqoutput == 0) then
         !---- Useless, no output will be created. ----------------------------------------!
         nrec_fast = 1
      elseif (iqoutput /= 0 .and. mod(day_sec,frqfast) /= 0.) then
         !---------------------------------------------------------------------------------!
         !    If mean diurnal cycle is on, then frqfast must be a divisor of one day.      !
         !---------------------------------------------------------------------------------!
         write(reason,fmt='(a,1x,f8.2,a,2x,1x,es14.7,a)')                                  &
            'FRQFAST must be a divisor of ',day_sec,'sec when mean diurnal cycle is on.'   &
           ,'Yours is set to ',frqfast,'sec...'
         call opspec_fatal(reason,'opspec_times')

      elseif ((imoutput /= 0 .or. idoutput /= 0 .or. outfast == -1. .or. outfast == -2. )  &
             .and. (mod(day_sec,frqfast) /= 0. .and. mod(frqfast,day_sec) /= 0.)) then
         !---------------------------------------------------------------------------------!
         !    If unitfrq is 0, then frqfast must be either a divisor or a multiple of one  !
         ! day in case there is daily and/or monthly analysis.                             !
         !---------------------------------------------------------------------------------!
         write(reason,fmt='(a,1x,f8.2,1x,a,1x,es14.7)')  &
             'FRQFAST must be a divisor or an integer multiple of ',day_sec,               &
             ' sec when daily and/or monthly analysis are on.  Yours is set to ',frqfast
         call opspec_fatal(reason,'opspec_times')
         ifaterr = ifaterr + 1
      !------------------------------------------------------------------------------------!
      !   Check outfast setting and adjusting it if needed. Depending on the               !
      ! configuration, this will cause the run to crash.                                   !
      !------------------------------------------------------------------------------------!
      elseif (frqfast < 10*60) then
         write(reason,fmt='(a,1x,a,1x,es14.7)')                                    &
             'FRQFAST must be greater than 10min (600 sec) when daily and/or monthly'    &
            ,'analysis are on or you will create a memory leak. Yours is set to ',frqfast
         call opspec_fatal(reason,'opspec_misc')
         ifaterr = ifaterr +1
      elseif(outfast == 0.) then
         !----- User didn't specify any outfast, use frqfast ------------------------------!
         outfast    = frqfast
         nrec_fast  = 1
      elseif (outfast == -1.) then
         outfast = day_sec
         nrec_fast  = int(outfast/frqfast)
      elseif (outfast == -2.) then
         nrec_fast = 0 !---- This must be reset every month.
      !----- User set outfast < frqfast. Resetting it and printing a warning --------------!
      elseif(outfast > 0. .and. outfast < frqfast) then
         outfast = frqfast
         nrec_fast  = 1 
         write (unit=*,fmt='(a)') ' '
         write (unit=*,fmt='(a)') '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         write (unit=*,fmt='(a)') '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         write (unit=*,fmt='(a)') '!!!  WARNING! WARNING! WARNING! WARNING! WARNING!   !!!'
         write (unit=*,fmt='(a)') '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         write (unit=*,fmt='(a)') '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         write (unit=*,fmt='(a)') ' -> Outfast cannot be less than frqfast.'
         write (unit=*,fmt='(a,1x,f7.0,1x,a)')                                             &
                                  '    Oufast was redefined to ',outfast,'seconds.'
         write (unit=*,fmt='(a)') '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         write (unit=*,fmt='(a)') '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         write (unit=*,fmt='(a)') ' '
      elseif(mod(outfast,frqfast) /= 0.0) then
         call opspec_fatal('OUTFAST must be a multiple of FRQFAST','opspec_times')
         ifaterr = ifaterr + 1
      else
         nrec_fast = int(outfast/frqfast)
      end if

   !---------------------------------------------------------------------------------------!
   !    Days.                                                                              !
   !---------------------------------------------------------------------------------------!
   case (1)
      if (ifoutput == 0) then
         nrec_fast = 1 ! Useless, no output will be created.
      elseif (iqoutput /= 0) then
         !---------------------------------------------------------------------------------!
         !    If mean diurnal cycle is on, frqfast must be given in seconds.               !
         !---------------------------------------------------------------------------------!
         write(reason,fmt='(a,2x,a)')                                                      &
            'FRQFAST must be given in sec when IQOUTPUT (mean diurnal cycle) is not zero.' &
           ,'Your UNITFAST is set to 2, which means months, check your namelist...'
         call opspec_fatal(reason,'opspec_times')

      elseif ((imoutput /= 0 .or. idoutput /= 0) .and.                                     &
              (mod(day_sec,day_sec*frqfast) /= 0. .and. mod(frqfast,1.) /= 0.)) then
         !---------------------------------------------------------------------------------!
         !    If unitfrq is 0, then frqfast must be either a divisor or a multiple of one  !
         ! day in case there is daily and/or monthly analysis.                             !
         !---------------------------------------------------------------------------------!
         write(reason,fmt='(a,1x,f8.2,1x,a,1x,es14.7)')                                    &
             'FRQFAST must be a divisor or an integer multiple of ',1,                     &
             ' day when daily and/or monthly analysis are on. Yours is set to ',frqfast
         call opspec_fatal(reason,'opspec_times')

      !------------------------------------------------------------------------------------!
      !   Checking outfast setting and adjusting it if needed. Depending on the            !
      ! configuration, this will cause the run to crash.                                   !
      !------------------------------------------------------------------------------------!
      !----- User didn't specify any outfast, use frqfast ---------------------------------!
      elseif(outfast == 0.) then
         outfast    = frqfast
         nrec_fast  = 1
      elseif (outfast == -1.) then
         outfast    = 1.
         nrec_fast  = int(outfast/frqfast)
      elseif (outfast == -2.) then
         nrec_fast = 0 !---- This must be reset every month.
      !----- User set outfast < frqfast. Resetting it and printing a warning --------------!
      elseif(outfast > 0. .and. outfast < frqfast) then
         outfast = frqfast
         nrec_fast  = 1 
         write (unit=*,fmt='(a)') ' '
         write (unit=*,fmt='(a)') '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         write (unit=*,fmt='(a)') '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         write (unit=*,fmt='(a)') '!!!  WARNING! WARNING! WARNING! WARNING! WARNING!   !!!'
         write (unit=*,fmt='(a)') '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         write (unit=*,fmt='(a)') '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         write (unit=*,fmt='(a)') ' -> Outfast cannot be less than frqfast.'
         write (unit=*,fmt='(a,1x,f7.0,1x,a)')                                             &
                                  '    Oufast was redefined to ',outfast,'days.'
         write (unit=*,fmt='(a)') '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         write (unit=*,fmt='(a)') '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         write (unit=*,fmt='(a)') ' '
      elseif(mod(outfast,frqfast) /= 0.0) then
         call opspec_fatal('OUTFAST must be a multiple of FRQFAST','opspec_times')
         ifaterr = ifaterr + 1
      else
         nrec_fast = int(outfast/frqfast)
      end if

      !----- Since days are always 86,400 sec, use seconds instead. -----------------------!
      frqfast = frqfast * day_sec
      if (outfast /= -2. ) outfast = outfast * day_sec


   !---------------------------------------------------------------------------------------!
   !    Months.                                                                            !
   !---------------------------------------------------------------------------------------!
   case (2)
      if (ifoutput == 0) then
         nrec_fast = 1 ! Useless, no output will be created.
      elseif (iqoutput /= 0) then
         !---------------------------------------------------------------------------------!
         !    If mean diurnal cycle is on, frqfast must be given in seconds.               !
         !---------------------------------------------------------------------------------!
         write(reason,fmt='(a,2x,a)')                                                      &
            'FRQFAST must be given in sec when IQOUTPUT (mean diurnal cycle) is not zero.' &
           ,'Your UNITFAST is set to 2, which means months, check your namelist...'
         call opspec_fatal(reason,'opspec_times')
      !------------------------------------------------------------------------------------!
      !    If unifrq is monthly, it needs to be a round number.                            !
      !------------------------------------------------------------------------------------!
      elseif (mod(frqfast,1.) /= 0.0) then
         write(reason,fmt='(a,1x,f8.2,1x,a,1x,es14.7)')                                    &
             'FRQFAST must be a round number when unitfast is ',unitfast,                  &
             '(months). Yours is currently set to ',frqfast
         call opspec_fatal(reason,'opspec_times')  
         ifaterr = ifaterr + 1
      !------------------------------------------------------------------------------------!
      !    If unifrq is monthly, it needs to be a divisor of 12, so it closes a year cycle.!
      ! This could be softened at some point, but that would require more calculations.    !
      !  Therefore, the only values accepted now are 1,2,3,4,6, and 12.                    !
      !------------------------------------------------------------------------------------!
      elseif (mod(12.,frqfast) /= 0.0) then
         write(reason,fmt='(a,1x,f8.2,1x,a,1x,es14.7)')                                    &
             'FRQFAST must be a divisor of 12 (one year) when unitfast is ',unitfast,      &
             '(months). Yours is currently set to ',frqfast
         call opspec_fatal(reason,'opspec_times')  
         ifaterr = ifaterr + 1
      !------------------------------------------------------------------------------------!
      !    This is fine but now outfast must be set exactly as frqfast. If the user wasn't !
      ! aware of this, print an informative banner.                                        !
      !------------------------------------------------------------------------------------!
      elseif (outfast /= 0. .or. outfast > frqfast) then
         outfast = frqfast
         nrec_fast = 1
         write (unit=*,fmt='(a)') ' '
         write (unit=*,fmt='(a)') '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         write (unit=*,fmt='(a)') '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         write (unit=*,fmt='(a)') '!!!  WARNING! WARNING! WARNING! WARNING! WARNING!   !!!'
         write (unit=*,fmt='(a)') '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         write (unit=*,fmt='(a)') '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         write (unit=*,fmt='(a)') ' -> Outfast cannot be different than frqfast when'
         write (unit=*,fmt='(a)') '    unitfast is set to 2 (months).'
         write (unit=*,fmt='(a,1x,f7.0,1x,a)')                                             &
                                  '    Oufast was redefined to ',outfast,'months.'
         write (unit=*,fmt='(a)') '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         write (unit=*,fmt='(a)') '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         write (unit=*,fmt='(a)') ' '
      else !----- The user either knew or was lucky, don't print the banner ---------------!
         outfast = frqfast
         nrec_fast = 1
      end if


   !---------------------------------------------------------------------------------------!
   !    Years.                                                                             !
   !---------------------------------------------------------------------------------------!
   case (3)
      if (ifoutput == 0) then
         nrec_fast = 1 ! Useless, no output will be created.
      ! If unifrq is yearly, then frqfast must be integer ---------------------------------!
      elseif (mod(frqfast,1.) /= 0.0) then
         write(reason,fmt='(a,1x,f8.2,1x,a,1x,es14.7)')                                    &
             'FRQFAST must be a round number when unitfast is ',unitfast,                  &
             '(years). Yours is currently set to ',frqfast
         call opspec_fatal(reason,'opspec_times')  
         ifaterr = ifaterr + 1
      elseif (iqoutput /= 0) then
         !---------------------------------------------------------------------------------!
         !    If mean diurnal cycle is on, frqfast must be given in seconds.               !
         !---------------------------------------------------------------------------------!
         write(reason,fmt='(a,2x,a)')                                                      &
            'FRQFAST must be given in sec when IQOUTPUT (mean diurnal cycle) is not zero.' &
           ,'Your UNITFAST is set to 3, which means years, check your namelist...'
         call opspec_fatal(reason,'opspec_times')
      !------------------------------------------------------------------------------------!
      !    This is fine but now outfast must be set exactly as frqfast. If the user wasn't !
      ! aware of this, print an informative banner.                                        !
      !------------------------------------------------------------------------------------!
      elseif (outfast /= 0. .or. outfast > frqfast) then
         outfast = frqfast
         nrec_fast = 1
         write (unit=*,fmt='(a)') ' '
         write (unit=*,fmt='(a)') '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         write (unit=*,fmt='(a)') '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         write (unit=*,fmt='(a)') '!!!  WARNING! WARNING! WARNING! WARNING! WARNING!   !!!'
         write (unit=*,fmt='(a)') '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         write (unit=*,fmt='(a)') '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         write (unit=*,fmt='(a)') ' -> Outfast cannot be different than frqfast when'
         write (unit=*,fmt='(a)') '    unitfast is set to 3 (years).'
         write (unit=*,fmt='(a,1x,f7.0,1x,a)')                                             &
                                  '    Oufast was redefined to ',outfast,'years.'
         write (unit=*,fmt='(a)') '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         write (unit=*,fmt='(a)') '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         write (unit=*,fmt='(a)') ' '
      else !----- The user either knew or was lucky, don't print the banner ---------------!
         outfast = frqfast
         nrec_fast = 1
      end if


   !---------------------------------------------------------------------------------------!
   !    Invalid unit, stopping the run.                                                    !
   !---------------------------------------------------------------------------------------!
   case default
      write (reason,fmt='(a,1x,i4,a)') &
        'Invalid UNITFAST, it must be between 0 and 3. Yours is set to',unitfast,'...'
      call opspec_fatal(reason,'opspec_times')  
      ifaterr = ifaterr +1
   end select
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !     Checking whether the user provided a valid combination of unitstate, frqstate,    !
   ! and outstate.                                                                         !
   !---------------------------------------------------------------------------------------!
   select case(unitstate)
   !---------------------------------------------------------------------------------------!
   !    Seconds                                                                            !
   !---------------------------------------------------------------------------------------!
   case (0)
      if (isoutput == 0) then
         nrec_state = 1 ! Useless, no output will be created.
      !------------------------------------------------------------------------------------!
      !    If unitstate is 0, then frqstate must be either a divisor or a multiple of one   !
      ! day in case there is daily and/or monthly analysis.                                !
      !------------------------------------------------------------------------------------!
      elseif ((imoutput /= 0 .or. idoutput /= 0 .or. outstate == -1. .or. outstate == -2.) &
             .and. (mod(day_sec,frqstate) /= 0. .and. mod(frqstate,day_sec) /= 0.)) then
         write(reason,fmt='(a,1x,f8.2,1x,a,1x,es14.7)')  &
             'FRQSTATE must be a divisor or an integer multiple of ',day_sec,              &
             ' sec when daily and/or monthly analysis are on. Yours is set to ',frqstate
         call opspec_fatal(reason,'opspec_times')
         ifaterr = ifaterr + 1
      !------------------------------------------------------------------------------------!
      !   Checking outstate setting and adjusting it if needed. Depending on the           !
      ! configuration, this will cause the run to crash.                                   !
      !------------------------------------------------------------------------------------!
      !----- User didn't specify any outstate, use frqstate -------------------------------!
      elseif(outstate == 0.) then
         outstate    = frqstate
         nrec_state  = 1
      elseif (outstate == -1.) then
         outstate = day_sec
         nrec_state  = int(outstate/frqstate)
      elseif (outstate == -2.) then
         nrec_state = 0 !---- This must be reset every month.
      !----- User set outstate < frqstate. Resetting it and printing a warning ------------!
      elseif(outstate > 0. .and. outstate < frqstate) then
         outstate = frqstate
         nrec_state  = 1 
         write (unit=*,fmt='(a)') ' '
         write (unit=*,fmt='(a)') '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         write (unit=*,fmt='(a)') '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         write (unit=*,fmt='(a)') '!!!  WARNING! WARNING! WARNING! WARNING! WARNING!   !!!'
         write (unit=*,fmt='(a)') '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         write (unit=*,fmt='(a)') '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         write (unit=*,fmt='(a)') ' -> Outstate cannot be less than frqstate.'
         write (unit=*,fmt='(a,1x,f7.0,1x,a)')                                             &
                                  '    Oustate was redefined to ',outstate,'seconds.'
         write (unit=*,fmt='(a)') '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         write (unit=*,fmt='(a)') '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         write (unit=*,fmt='(a)') ' '
      elseif(mod(outstate,frqstate) /= 0.0) then
         call opspec_fatal('OUTSTATE must be a multiple of FRQSTATE','opspec_times')
         ifaterr = ifaterr + 1
      else
         nrec_state = int(outstate/frqstate)
      end if

   !---------------------------------------------------------------------------------------!
   !    Days.                                                                              !
   !---------------------------------------------------------------------------------------!
   case (1)
      if (isoutput == 0) then
         nrec_state = 1 ! Useless, no output will be created.
      !------------------------------------------------------------------------------------!
      !    If unitfrq is 0, then frqstate must be either a divisor or a multiple of one    !
      ! day in case there is daily and/or monthly analysis.                                !
      !------------------------------------------------------------------------------------!
      elseif ((imoutput /= 0 .or. idoutput /= 0) .and.                                     &
              (mod(day_sec,day_sec*frqstate) /= 0. .and. mod(frqstate,1.) /= 0.)) then
         write(reason,fmt='(a,1x,f8.2,1x,a,1x,es14.7)')                                    &
             'FRQSTATE must be a divisor or an integer multiple of ',1,                     &
             ' day when daily and/or monthly analysis are on. Yours is set to ',frqstate
         call opspec_fatal(reason,'opspec_times')

      !------------------------------------------------------------------------------------!
      !   Checking outstate setting and adjusting it if needed. Depending on the           !
      ! configuration, this will cause the run to crash.                                   !
      !------------------------------------------------------------------------------------!
      !----- User didn't specify any outstate, use frqstate ---------------------------------!
      elseif(outstate == 0.) then
         outstate    = frqstate
         nrec_state  = 1
      elseif (outstate == -1.) then
         outstate    = 1.
         nrec_state  = int(outstate/frqstate)
      elseif (outstate == -2.) then
         nrec_state = 0 !---- This must be reset every month.
      !----- User set outstate < frqstate. Resetting it and printing a warning --------------!
      elseif(outstate > 0. .and. outstate < frqstate) then
         outstate = frqstate
         nrec_state  = 1 
         write (unit=*,fmt='(a)') ' '
         write (unit=*,fmt='(a)') '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         write (unit=*,fmt='(a)') '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         write (unit=*,fmt='(a)') '!!!  WARNING! WARNING! WARNING! WARNING! WARNING!   !!!'
         write (unit=*,fmt='(a)') '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         write (unit=*,fmt='(a)') '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         write (unit=*,fmt='(a)') ' -> Outstate cannot be less than frqstate.'
         write (unit=*,fmt='(a,1x,f7.0,1x,a)')                                             &
                                  '    Oustate was redefined to ',outstate,'days.'
         write (unit=*,fmt='(a)') '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         write (unit=*,fmt='(a)') '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         write (unit=*,fmt='(a)') ' '
      elseif(mod(outstate,frqstate) /= 0.0) then
         call opspec_fatal('OUTSTATE must be a multiple of FRQSTATE','opspec_times')
         ifaterr = ifaterr + 1
      else
         nrec_state = int(outstate/frqstate)
      end if

      !----- Since days are always 86,400 sec, use seconds instead. -----------------------!
      frqstate = frqstate * day_sec
      if (outstate /= -2. ) outstate = outstate * day_sec


   !---------------------------------------------------------------------------------------!
   !    Months.                                                                            !
   !---------------------------------------------------------------------------------------!
   case (2)
      if (isoutput == 0) then
         nrec_state = 1 ! Useless, no output will be created.
      !------------------------------------------------------------------------------------!
      !    If unifrq is monthly, it needs to be a round number.                            !
      !------------------------------------------------------------------------------------!
      elseif (mod(frqstate,1.) /= 0.) then
         write(reason,fmt='(a,1x,f8.2,1x,a,1x,es14.7)')                                    &
             'FRQSTATE must be a round number when unitstate is ',unitstate,               &
             '(months). Yours is currently set to ',frqstate
         call opspec_fatal(reason,'opspec_times')  
         ifaterr = ifaterr + 1
      !------------------------------------------------------------------------------------!
      !    If unistate is monthly, it needs to be a divisor of 12, so it closes a year     !
      ! cycle. This could be softened at some point, but that would require more           !
      ! calculations.  Therefore, the only values accepted now are 1,2,3,4,6, and 12.      !
      !------------------------------------------------------------------------------------!
      elseif (mod(12.,frqstate) /= 0.) then
         write(reason,fmt='(a,1x,f8.2,1x,a,1x,es14.7)')                                    &
             'FRQSTATE must be a divisor of 12 (one year) when unitstate is ',unitstate,   &
             '(months). Yours is currently set to ',frqstate
         call opspec_fatal(reason,'opspec_times')  
         ifaterr = ifaterr + 1
      !------------------------------------------------------------------------------------!
      !    This is fine but now outstate must be set exactly as frqstate. If the user wasn't !
      ! aware of this, print an informative banner.                                        !
      !------------------------------------------------------------------------------------!
      elseif (outstate /= 0. .and. outstate > frqstate) then
         outstate = frqstate
         nrec_state = 1
         write (unit=*,fmt='(a)') ' '
         write (unit=*,fmt='(a)') '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         write (unit=*,fmt='(a)') '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         write (unit=*,fmt='(a)') '!!!  WARNING! WARNING! WARNING! WARNING! WARNING!   !!!'
         write (unit=*,fmt='(a)') '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         write (unit=*,fmt='(a)') '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         write (unit=*,fmt='(a)') ' -> Outstate cannot be different than frqstate when'
         write (unit=*,fmt='(a)') '    unitstate is set to 2 (months).'
         write (unit=*,fmt='(a,1x,f7.0,1x,a)')                                             &
                                  '    Oustate was redefined to ',outstate,'months.'
         write (unit=*,fmt='(a)') '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         write (unit=*,fmt='(a)') '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         write (unit=*,fmt='(a)') ' '
      else !----- The user either knew or was lucky, don't print the banner ---------------!
         outstate = frqstate
         nrec_state = 1
      end if


   !---------------------------------------------------------------------------------------!
   !    Years.                                                                             !
   !---------------------------------------------------------------------------------------!
   case (3)
      if (isoutput == 0) then
         nrec_state = 1 ! Useless, no output will be created.
      ! If unifrq is yearly, then frqstate must be integer ---------------------------------!
      elseif (mod(frqstate,1.) /= 0.0) then
         write(reason,fmt='(a,1x,f8.2,1x,a,1x,es14.7)')                                    &
             'FRQSTATE must be a round number when unitstate is ',unitstate,                  &
             '(years). Yours is currently set to ',frqstate
         call opspec_fatal(reason,'opspec_times')  
         ifaterr = ifaterr + 1
      !------------------------------------------------------------------------------------!
      !    This is fine but now outstate must be set exactly as frqstate. If the user      !
      ! wasn't aware of this, print an informative banner.                                 !
      !------------------------------------------------------------------------------------!
      elseif (outstate /= 0. .or. outstate > frqstate) then

         write (unit=*,fmt='(a)') ' '
         write (unit=*,fmt='(a)') '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         write (unit=*,fmt='(a)') '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         write (unit=*,fmt='(a)') '!!!  WARNING! WARNING! WARNING! WARNING! WARNING!   !!!'
         write (unit=*,fmt='(a)') '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         write (unit=*,fmt='(a)') '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         write (unit=*,fmt='(a)') ' -> Outstate cannot be different than frqstate when'
         write (unit=*,fmt='(a)') '    unitstate is set to 3 (years).'
         write (unit=*,fmt='(a,1x,f7.0,1x,a)')                                             &
                                  '    Oustate was set to ',outstate,'years.'
         write (unit=*,fmt='(a,1x,f7.0,1x,a)')                                             &
                                  '    Oustate was redefined to ',frqstate,'years.'
         write (unit=*,fmt='(a)') '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         write (unit=*,fmt='(a)') '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         write (unit=*,fmt='(a)') ' '
         outstate = frqstate
         nrec_state = 1
      else !----- The user either knew or was lucky, don't print the banner ---------------!
         outstate = frqstate
         nrec_state = 1
      end if


   !---------------------------------------------------------------------------------------!
   !    Invalid unit, stopping the run.                                                    !
   !---------------------------------------------------------------------------------------!
   case default
      write (reason,fmt='(a,1x,i4,a)') &
        'Invalid UNITSTATE, it must be between 0 and 3. Yours is set to',unitstate,'...'
      call opspec_fatal(reason,'opspec_times')  
      ifaterr = ifaterr +1
   end select
   !---------------------------------------------------------------------------------------!





   !---------------------------------------------------------------------------------------!
   !    Frequency of the analysis must be a divisor of the frequency of history output     !
   ! The history files contain variables that are integrated over the analysis frequency.  !
   ! If the analysis and history times are not synchronized, then the integrated values    !
   ! will be improperly normalized and give strange/biased results. Notes: If it is abso-  !
   ! lutely necessary to have them non-divisible, then do so, but realize any integrated   !
   ! fast-time variables are meaningless, or scale them yourself. Note 2: If there are no  !
   ! analysis outputs, only history, then the analysis time will be set to the history     !
   ! time, so that the integrated variables reported in the history file are integrated    !
   ! over the history period.                                                              !
   !---------------------------------------------------------------------------------------!
   if (isoutput /= 0 .and. (ifoutput /= 0 .or. iqoutput /= 0)) then
      if (unitfast == unitstate .and. mod(frqstate,frqfast) /= 0.0 ) then
         write(reason,fmt='(a,1x,2(a,1x,f10.2,1x))')                                       &
              'FRQFAST must be a divisor of FRQSTATE if both are outputing data.',         &
              'Yours is set to FRQFAST=',frqfast,' sec and FRQSTATE=',frqstate,'days!'
         call opspec_fatal(reason,'opspec_times')  
         ifaterr=ifaterr+1
      end if
      !----- Now we check other combinations ----------------------------------------------!
      select case (unitfast)
      case (0) !----- Seconds -------------------------------------------------------------!
         select case(unitstate)
         case (1) !---- state in days -----------------------------------------------------!
            if (mod(frqstate*day_sec,frqfast) /= 0.0) then
               write(reason,fmt='(a,1x,2(a,1x,f10.2,1x))')                                 &
                    'FRQFAST must be a divisor of FRQSTATE if both are outputing data.',   &
                    'Yours is set to FRQFAST=',frqfast,'sec and FRQSTATE=',frqstate,'days!'
               ifaterr=ifaterr+1
            end if
         case (2,3) !---- State in months or years, frqfast must be divisor of 1 day ------!
            if (mod(frqstate*day_sec,frqfast) /= 0.0) then
               write(reason,fmt='(a,1x,a,1x,f10.2,1x,a,1x,f10.2)')                         &
                    'FRQFAST must be a divisor of one day when frqfast is in seconds and', &
                    'FRQSTATE is in months or years. Your FRQFAST=',frqfast,'sec.'
               ifaterr=ifaterr+1
            end if
         end select

      case (1) !----- Days ----------------------------------------------------------------!
         select case(unitstate)
         case (0) !---- state in seconds --------------------------------------------------!
            if (mod(frqstate,frqfast*day_sec) /= 0.0) then
               write(reason,fmt='(a,1x,a,1x,f10.2,1x,a,1x,f10.2,1x,a)')                    &
                    'FRQFAST must be a divisor of FRQSTATE if both are outputing data.',   &
                    'Yours is set to FRQFAST=',frqfast,'days and FRQSTATE=',frqstate,'sec!'
               ifaterr=ifaterr+1
            end if

         case (2,3) !---- State in months or years, frqfast must be divisor of 1 day ------!
            if (frqfast /= 1.0) then
               write(reason,fmt='(a,1x,a,1x,f10.2,1x,a,1x,f10.2,1x,a)')                    &
                    'FRQFAST must be a divisor of one day or 1 day when UNITFAST is in',   &
                    'days and FRQSTATE is in months or years. Your FRQFAST=',frqfast,'days.'
               ifaterr=ifaterr+1
            end if

         end select

      case (2) !----- Months --------------------------------------------------------------!
         select case(unitstate)
         case (0,1) !---- state in seconds or days. Can't be since months are irregular. --!
            write(reason,fmt='(a,1x,a,1x,i5)')                                             &
               'If UNITFAST is in months, UNITSTATE must be either months or year, and',   &
               'currently UNITSTATE=',unitstate
            ifaterr=ifaterr+1
      
         case (3) !---- State years, frqfast must be divisor ------------------------------!
            if (mod(frqstate*12.,frqfast) /= 0.) then
               write(reason,fmt='(a,1x,2(a,1x,f10.2,1x),a)')                               &
                    'FRQFAST must be a divisor of FRQSTATE if both are outputing data.',   &
                    'Yours is set to FRQFAST=',frqfast,'mon and FRQSTATE=',frqstate,'yrs!'
            end if
            ifaterr=ifaterr+1

         end select

      case (3) !----- Years ---------------------------------------------------------------!
         select case(unitstate)
         case (0,1) !---- state in seconds or days. Can't be since months are irregular. --!
            write(reason,fmt='(a,1x,a,1x,i5)')                                             &
               'If UNITFAST is in months, UNITSTATE must be either months or year, and',   &
               'currently UNITSTATE=',unitstate
            ifaterr=ifaterr+1
      
         case (2) !---- State months, frqfast must be divisor -----------------------------!
            if (mod(frqstate,frqfast*12.) /= 0.) then
               write(reason,fmt='(a,1x,2(a,1x,f10.2,1x),a)')                               &
                    'FRQFAST must be a divisor of FRQSTATE if both are outputing data.',   &
                    'Yours is set to FRQFAST=',frqfast,'yrs and FRQSTATE=',frqstate,'mon!'
               ifaterr=ifaterr+1
            end if
         end select
      end select
   end if

   !----- 


   !Check if this simulation has a positive timmax.
   if (timmax < 0.0d0) then
      write(reason,fmt='(a,2(a,1x,2(i2.2,a),i4.4,1x,i4.4,a,1x))')     &
         'Your end time is before the initial time.'                  &
         ,'Initial:',imontha,'/',idatea,'/',iyeara,itimea,'GMT'       &
         ,'Final  :',imonthz,'/',idatez,'/',iyearz,itimez,'GMT'
      call opspec_fatal(reason,'opspec_times')  
      ifaterr=ifaterr+1
   end if

   !---------------------------------------------------------------------------------------!
   !    If we made up to this point, find frqsum so we make sure that the main time step   !
   ! and the radiation are both divisors of the frequency of output.                       !
   !---------------------------------------------------------------------------------------!
   if (ifaterr == 0) then
      call find_frqsum()

      !----- Check whether the time steps make sense. -------------------------------------!
      if (mod(frqsum,dtlsm) /= 0.0) then
         write(reason,fmt='(a,1x,f8.2,1x,a,1x,es14.7)')  &
             'DTLSM must be a divisor of ',frqsum,' sec. Yours is set to ',dtlsm
         call opspec_fatal(reason,'opspec_times')  
         ifaterr=ifaterr+1
      end if

      if (mod(frqsum,radfrq) /= 0.0) then
         write(reason,fmt='(a,1x,f8.2,1x,a,1x,es14.7)')  &
             'RADFRQ must be a divisor of ',frqsum,' sec. Yours is set to ',radfrq
         call opspec_fatal(reason,'opspec_times')  
         ifaterr=ifaterr+1
      end if
      !------------------------------------------------------------------------------------!
   end if
   !---------------------------------------------------------------------------------------!





   !----- DTLSM must be an integer divisor of RADFRQ so integrals make sense. -------------!
   if (mod(radfrq,dtlsm) /= 0.0) then
      write(reason,fmt='(a,1x,f8.2,1x,a,1x,f8.2,a)')  &
          'DTLSM must be a divisor of RADFRQ. Your DTLSM is set to',dtlsm, &
          'and your RADFRQ is set to',radfrq,'...'
      call opspec_fatal(reason,'opspec_times')  
      ifaterr=ifaterr+1
   end if
   !---------------------------------------------------------------------------------------!




   !------ Stop the run if there are any fatal errors. ------------------------------------!
   if (ifaterr > 0) then
      write (unit=*,fmt='(a)')       ' -----------ED_OPSPEC_TIMES ------------------------'
      write (unit=*,fmt='(a,1x,i5)') ' Fatal errors:',ifaterr
      write (unit=*,fmt='(a)')       ' ---------------------------------------------------'
      call fatal_error('Fatal errors at namelist - Time settings '&
                     & ,'ed_opspec_times','ed_opspec.f90')
   end if
   !---------------------------------------------------------------------------------------!




   return
end subroutine ed_opspec_times
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine performs miscellaneous tests over the options, like values outside    !
! the allowed ranges and conflicting dynamic settings.                                     !
!------------------------------------------------------------------------------------------!
subroutine ed_opspec_misc
   use ed_max_dims           , only : n_pft                        & ! intent(in)
                                    , str_len                      ! ! intent(in)
   use ed_misc_coms          , only : ifoutput                     & ! intent(in)
                                    , idoutput                     & ! intent(in)
                                    , iqoutput                     & ! intent(in)
                                    , imoutput                     & ! intent(in)
                                    , iyoutput                     & ! intent(in)
                                    , itoutput                     & ! intent(in)
                                    , isoutput                     & ! intent(in)
                                    , igoutput                     & ! intent(in)
                                    , iadd_site_means              & ! intent(in)
                                    , iadd_patch_means             & ! intent(in)
                                    , iadd_cohort_means            & ! intent(in)
                                    , iclobber                     & ! intent(in)
                                    , runtype                      & ! intent(in)
                                    , ied_init_mode                & ! intent(in)
                                    , ivegt_dynamics               & ! intent(in)
                                    , ibigleaf                     & ! intent(in)
                                    , integration_scheme           & ! intent(in)
                                    , iallom                       & ! intent(in)
                                    , igrass                       & ! intent(in)
                                    , growth_resp_scheme           & ! intent(in)
                                    , storage_resp_scheme          & ! intent(in)
                                    , min_site_area                ! ! intent(in)
   use canopy_air_coms       , only : icanturb                     & ! intent(in)
                                    , isfclyrm                     & ! intent(in)
                                    , ied_grndvap                  & ! intent(in)
                                    , ubmin                        & ! intent(in)
                                    , ugbmin                       & ! intent(in)
                                    , ustmin                       & ! intent(in)
                                    , gamm                         & ! intent(in)
                                    , gamh                         & ! intent(in)
                                    , tprandtl                     & ! intent(in)
                                    , ribmax                       & ! intent(in)
                                    , leaf_maxwhc                  ! ! intent(in)
   use soil_coms             , only : ed_nstyp                     & ! intent(in)
                                    , ed_nscol                     & ! intent(in)
                                    , isoilflg                     & ! intent(in)
                                    , nslcon                       & ! intent(in)
                                    , isoilcol                     & ! intent(in)
                                    , slxclay                      & ! intent(in)
                                    , slxsand                      & ! intent(in)
                                    , isoilstateinit               & ! intent(in)
                                    , isoildepthflg                & ! intent(in)
                                    , isoilbc                      & ! intent(in)
                                    , sldrain                      & ! intent(in)
                                    , zrough                       & ! intent(in)
                                    , runoff_time                  ! ! intent(in)
   use mem_polygons          , only : maxsite                      & ! intent(in)
                                    , maxpatch                     ! ! intent(in)
   use grid_coms             , only : ngrids                       ! ! intent(in)
   use physiology_coms       , only : iphysiol                     & ! intent(in)
                                    , h2o_plant_lim                & ! intent(in)
                                    , iddmort_scheme               & ! intent(in)
                                    , cbr_scheme                   & ! intent(in)
                                    , ddmort_const                 & ! intent(in)
                                    , n_plant_lim                  & ! intent(in)
                                    , vmfact_c3                    & ! intent(in)
                                    , vmfact_c4                    & ! intent(in)
                                    , mphoto_trc3                  & ! intent(in)
                                    , mphoto_tec3                  & ! intent(in)
                                    , mphoto_c4                    & ! intent(in)
                                    , bphoto_blc3                  & ! intent(in)
                                    , bphoto_nlc3                  & ! intent(in)
                                    , bphoto_c4                    & ! intent(in)
                                    , kw_grass                     & ! intent(in)
                                    , kw_tree                      & ! intent(in)
                                    , gamma_c3                     & ! intent(in)
                                    , gamma_c4                     & ! intent(in)
                                    , d0_grass                     & ! intent(in)
                                    , d0_tree                      & ! intent(in)
                                    , alpha_c3                     & ! intent(in)
                                    , alpha_c4                     & ! intent(in)
                                    , klowco2in                    & ! intent(in)
                                    , rrffact                      & ! intent(in)
                                    , growthresp                   & ! intent(in)
                                    , lwidth_grass                 & ! intent(in)
                                    , lwidth_bltree                & ! intent(in)
                                    , lwidth_nltree                & ! intent(in)
                                    , q10_c3                       & ! intent(in)
                                    , q10_c4                       & ! intent(in)
                                    , quantum_efficiency_T         ! ! intent(in)
   use decomp_coms           , only : n_decomp_lim                 & ! intent(in)
                                    , decomp_scheme                ! ! intent(in)
   use disturb_coms          , only : include_fire                 & ! intent(in)
                                    , fire_parameter               & ! intent(in)
                                    , ianth_disturb                & ! intent(in)
                                    , sm_fire                      & ! intent(in)
                                    , time2canopy                  & ! intent(in)
                                    , treefall_disturbance_rate    & ! intent(in)
                                    , min_patch_area               ! ! intent(in)
   use phenology_coms        , only : iphen_scheme                 & ! intent(in)
                                    , repro_scheme                 & ! intent(in)
                                    , radint                       & ! intent(in)
                                    , radslp                       & ! intent(in)
                                    , thetacrit                    ! ! intent(in)
   use pft_coms              , only : include_these_pft            & ! intent(in)
                                    , pft_1st_check                & ! intent(in)
                                    , agri_stock                   & ! intent(in)
                                    , plantation_stock             ! ! intent(in)
   use canopy_layer_coms     , only : crown_mod                    ! ! intent(in)
   use canopy_radiation_coms , only : icanrad                      & ! intent(in)
                                    , ihrzrad                      & ! intent(in)
                                    , ltrans_vis                   & ! intent(in)
                                    , ltrans_nir                   & ! intent(in)
                                    , lreflect_vis                 & ! intent(in)
                                    , lreflect_nir                 & ! intent(in)
                                    , orient_tree                  & ! intent(in)
                                    , orient_grass                 & ! intent(in)
                                    , clump_tree                   & ! intent(in)
                                    , clump_grass                  ! ! intent(in)
   use rk4_coms              , only : ibranch_thermo               & ! intent(in)
                                    , ipercol                      & ! intent(in)
                                    , rk4_tolerance                ! ! intent(in)
   use mem_polygons          , only : n_ed_region                  & ! intent(in)
                                    , n_poi                        ! ! intent(in)
   use detailed_coms         , only : dt_census                    & ! intent(in)
                                    , yr1st_census                 & ! intent(in)
                                    , mon1st_census                & ! intent(in)
                                    , min_recruit_dbh              & ! intent(in)
                                    , idetailed                    & ! intent(in)
                                    , patch_keep                   ! ! intent(in)

   use met_driver_coms       , only : imetrad                      ! ! intent(in)
#if defined(COUPLED)
#else
   use met_driver_coms       , only : ishuffle                     & ! intent(in)
                                    , imetavg                      ! ! intent(in)
#endif

   implicit none
   !----- Local variables -----------------------------------------------------------------!
   character(len=str_len) :: reason
   integer                :: ifaterr
   integer                :: ifm
   integer                :: ipft
   logical                :: agri_ok
   logical                :: plantation_ok
   logical                :: patch_detailed
   !----- Local constants -----------------------------------------------------------------!
   integer, parameter :: skip=huge(6)
   !---------------------------------------------------------------------------------------!

   !----- IFATERR will count the number of bad set ups. -----------------------------------!
   ifaterr=0


   if (maxsite < 1 .and. maxsite > ed_nstyp) then
      write (reason,fmt='(a,1x,i5,a,2x,a,1x,i5,a)')                                        &
         'Invalid MAXSITE, it must be between 1 and ed_nstyp (',ed_nstyp,').'              &
        ,'Yours is set to ',maxsite,'...'
   end if

   if (min_site_area < 0.0001 .or. min_site_area > 0.10) then
      write (reason,fmt='(a,2x,a,1x,es12.5,a)')                                            &
         'Invalid MIN_SITE_AREA, it must be between 0.0001 and 0.10.'                      &
        ,'Yours is set to ',min_site_area,'...'
   end if

   if (min_patch_area < 0.000001 .or. min_patch_area > 0.02) then
      write (reason,fmt='(a,2x,a,1x,es12.5,a)')                                            &
         'Invalid MIN_PATCH_AREA, it must be between 0.000001 and 0.02.'                   &
        ,'Yours is set to ',min_patch_area,'...'
   end if

   if (ifoutput /= 0 .and. ifoutput /= 3) then
      write (reason,fmt='(a,1x,i4,a)')                                                     &
        'Invalid IFOUTPUT, it must be 0 (none) or 3 (HDF5). Yours is set to',ifoutput,'...'
      call opspec_fatal(reason,'opspec_misc')
      ifaterr = ifaterr +1
   end if
   if (idoutput /= 0 .and. idoutput /= 3) then
      write (reason,fmt='(a,1x,i4,a)')                                                     &
        'Invalid IDOUTPUT, it must be 0 (none) or 3 (HDF5). Yours is set to',idoutput,'...'
      call opspec_fatal(reason,'opspec_misc')
      ifaterr = ifaterr +1
   end if
   if (imoutput /= 0 .and. imoutput /= 3) then
      write (reason,fmt='(a,1x,i4,a)')                                                     &
        'Invalid IMOUTPUT, it must be 0 (none) or 3 (HDF5). Yours is set to',imoutput,'...'
      call opspec_fatal(reason,'opspec_misc')
      ifaterr = ifaterr +1
   end if
   if (iqoutput /= 0 .and. iqoutput /= 3) then
      write (reason,fmt='(a,1x,i4,a)')                                                     &
        'Invalid IQOUTPUT, it must be 0 (none) or 3 (HDF5). Yours is set to',iqoutput,'...'
      call opspec_fatal(reason,'opspec_misc')
      ifaterr = ifaterr +1
   end if
   if (iyoutput /= 0 .and. iyoutput /= 3) then
      write (reason,fmt='(a,1x,i4,a)')                                                     &
        'Invalid IYOUTPUT, it must be 0 (none) or 3 (HDF5). Yours is set to',iyoutput,'...'
      call opspec_fatal(reason,'opspec_misc')  
      ifaterr = ifaterr +1
   end if
   if (itoutput /= 0 .and. itoutput /= 3) then
      write (reason,fmt='(a,1x,i4,a)')                                                     &
        'Invalid ITOUTPUT, it must be 0 (none) or 3 (HDF5). Yours is set to',itoutput,'...'
      call opspec_fatal(reason,'opspec_misc')  
      ifaterr = ifaterr +1
   end if
   if (isoutput /= 0 .and. isoutput /= 3) then
      write (reason,fmt='(a,1x,i4,a)')                                                     &
        'Invalid ISOUTPUT, it must be 0 (none) or 3 (HDF5). Yours is set to',isoutput,'...'
      call opspec_fatal(reason,'opspec_misc')
      ifaterr = ifaterr +1
   end if
   if (iadd_site_means < 0 .or. iadd_site_means > 1) then
      write (reason,fmt='(a,1x,i4,a)')                                                     &
        'Invalid IADD_SITE_MEANS, it must be 0 (no) or 1 (yes).  Yours is set to'          &
       ,iadd_site_means,'...'
      call opspec_fatal(reason,'opspec_misc')
      ifaterr = ifaterr +1
   end if
   if (iadd_patch_means < 0 .or. iadd_patch_means > 1) then
      write (reason,fmt='(a,1x,i4,a)')                                                     &
        'Invalid IADD_PATCH_MEANS, it must be 0 (no) or 1 (yes).  Yours is set to'         &
       ,iadd_patch_means,'...'
      call opspec_fatal(reason,'opspec_misc')
      ifaterr = ifaterr +1
   end if
   if (iadd_cohort_means < 0 .or. iadd_cohort_means > 1) then
      write (reason,fmt='(a,1x,i4,a)')                                                     &
        'Invalid IADD_COHORT_MEANS, it must be 0 (no) or 1 (yes).  Yours is set to'        &
       ,iadd_cohort_means,'...'
      call opspec_fatal(reason,'opspec_misc')
      ifaterr = ifaterr +1
   end if
   if (iclobber < 0 .or. iclobber > 1) then
      write (reason,fmt='(a,1x,i4,a)')                                                     &
        'Invalid ICLOBBER, it must be 0 or 1. Yours is set to',iclobber,'...'
      call opspec_fatal(reason,'opspec_misc')
      ifaterr = ifaterr +1
   end if
   if (trim(runtype) /= 'INITIAL' .and. trim(runtype) /= 'HISTORY') then
      write (reason,fmt='(a,1x,2a)')                                                       &
                         'Invalid RUNTYPE, it must be INITIAL or HISTORY. Yours is set to' &
                         ,trim(runtype),'...'
      call opspec_fatal(reason,'opspec_misc')
      ifaterr = ifaterr +1
   end if

   if (ied_init_mode == -8) then 
      !------------------------------------------------------------------------------------!
      !     The special 8-layer model works only in size- and age-structured runs.         !
      !------------------------------------------------------------------------------------!
      if (ibigleaf == 1) then
         write (reason,fmt='(a)')                                                          &
                            'IED_INIT_MODE can''t be -8 when running big leaf mode.'       &
                            ,trim(runtype),'...'
         call opspec_fatal(reason,'opspec_misc')
         ifaterr = ifaterr +1
      end if
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     This is just for idealised test runs and shouldn't be used as a regular        !
      ! option.                                                                            !
      !------------------------------------------------------------------------------------!
      write (unit=*,fmt='(a)') '==========================================================='
      write (unit=*,fmt='(a)') '   WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!   '
      write (unit=*,fmt='(a)') '   WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!   '
      write (unit=*,fmt='(a)') '   WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!   '
      write (unit=*,fmt='(a)') '   WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!   '
      write (unit=*,fmt='(a)') '   WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!   '
      write (unit=*,fmt='(a)') '   WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!   '
      write (unit=*,fmt='(a)') '==========================================================='
      write (unit=*,fmt='(a)') '    You have chosen to run a prescribed 8-layer run with   '
      write (unit=*,fmt='(a)') ' no ecosystem dynamics.  This should be used for TEST      '
      write (unit=*,fmt='(a)') ' simulations only.  If that''s not what you wanted, change '
      write (unit=*,fmt='(a)') ' your IED_INIT_MODE variable on your ED2IN.                '
      write (unit=*,fmt='(a)') '==========================================================='
   elseif ((ied_init_mode < -1 .or. ied_init_mode > 7) .and. &
           (ied_init_mode /= 99 )) then
      write (reason,fmt='(a,1x,i4,a)')                                                     &
                     'Invalid IED_INIT_MODE, it must be between -1 and 7. Yours is set to' &
                    ,ied_init_mode,'...'
      call opspec_fatal(reason,'opspec_misc')
      ifaterr = ifaterr +1
   end if



   if (ied_init_mode == 7 .and. isoilstateinit>0 ) then
      write (reason,fmt='(a)')                                                   &
           'Please set ISOILSTATEINIT=0 if using IED_INIT_MODE=7'
      call opspec_fatal(reason,'opspec_misc')
      ifaterr = ifaterr +1
   end if

      
#if defined(COUPLED)
   do ifm=1,ngrids
      if (isoilflg(ifm) < 0 .or. isoilflg(ifm) > 3) then
         write (reason,fmt='(a,1x,i4,1x,a,1x,i4,a)')                                       &
                       'Invalid ISOILFLG, it must be between 0 and 3. Yours is set to'     &
                       ,isoilflg(ifm),'for grid',ifm,'...'
         call opspec_fatal(reason,'opspec_misc')
         ifaterr = ifaterr +1
      end if
   end do
#else
   do ifm=1,ngrids
      if (isoilflg(ifm) < 1 .or. isoilflg(ifm) > 2) then
         write (reason,fmt='(a,1x,i4,1x,a,1x,i4,a)')                                       &
                       'Invalid ISOILFLG, it must be between 1 and 3. Yours is set to'     &
                      ,isoilflg(ifm),'for grid',ifm,'...'
         call opspec_fatal(reason,'opspec_misc')
         ifaterr = ifaterr +1
      end if
   end do
#endif

   if (nslcon < 1 .or. nslcon > ed_nstyp) then
      write (reason,fmt='(2(a,1x,i4),a)')                                                  &
             'Invalid NSLCON, it must be between 1 and ',ed_nstyp,'. Yours is set to'      &
            ,nslcon,'...'
      call opspec_fatal(reason,'opspec_misc')
      ifaterr = ifaterr +1
   end if

   if (isoilcol < 1 .or. isoilcol > ed_nscol) then
      write (reason,fmt='(2(a,1x,i4),a)')                                                  &
             'Invalid ISOILCOL, it must be between 1 and ',ed_nscol,'. Yours is set to'    &
            ,isoilcol,'...'
      call opspec_fatal(reason,'opspec_misc')
      ifaterr = ifaterr +1
   end if

do ifm=1,ngrids
   if (isoilflg(ifm)==1 .and. slxclay>0. .and. slxclay<1. .and. slxsand>0. .and. slxsand<1.) then
      write (unit=*,fmt='(a)') '==========================================================='
      write (unit=*,fmt='(a)') '   WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!   '
      write (unit=*,fmt='(a)') '   WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!   '
      write (unit=*,fmt='(a)') '   WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!   '
      write (unit=*,fmt='(a)') '   WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!   '
      write (unit=*,fmt='(a)') '   WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!   '
      write (unit=*,fmt='(a)') '   WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!   '
      write (unit=*,fmt='(a)') '==========================================================='
      write (unit=*,fmt='(a)') '    You have set up the run to read in soil types based on '
      write (unit=*,fmt='(a)') ' Lat/Lon.  This overrides the option to use user defined   '
      write (unit=*,fmt='(a)') ' sand and clay fractions- the DEFAULT SOIL PARAMETERS WILL '
      write (unit=*,fmt='(a)') ' BE USED. To use the NL%SLXCLAY and NL%SLXSAND options, set'
      write (unit=*,fmt='(a)') ' NL%ISOILFLG to 2 and define the soil type through NL%NLSCON'
      write (unit=*,fmt='(a)') '==========================================================='
   end if


   if (isoilflg(ifm)==2 .and. slxclay<0.) then
      write (unit=*,fmt='(a)') '==========================================================='
      write (unit=*,fmt='(a)') '   WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!   '
      write (unit=*,fmt='(a)') '==========================================================='
      write (unit=*,fmt='(a)') 'You have defined a soil clay fraction < 0. The DEFAULT SOIL'
      write (unit=*,fmt='(a)') 'PARAMETERS WILL BE USED based on NL%NLSCON.                '
      write (unit=*,fmt='(a)') '==========================================================='
   end if
   
   if (isoilflg(ifm)==2 .and. slxclay>1.) then
      write (unit=*,fmt='(a)') '==========================================================='
      write (unit=*,fmt='(a)') '   WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!   '
      write (unit=*,fmt='(a)') '==========================================================='
      write (unit=*,fmt='(a)') 'You have defined a soil clay fraction > 1. The DEFAULT SOIL'
      write (unit=*,fmt='(a)') 'PARAMETERS WILL BE USED based on NL%NLSCON.                '
      write (unit=*,fmt='(a)') '==========================================================='
   end if

   if (isoilflg(ifm)==2 .and. slxsand<0.) then
      write (unit=*,fmt='(a)') '==========================================================='
      write (unit=*,fmt='(a)') '   WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!   '
      write (unit=*,fmt='(a)') '==========================================================='
      write (unit=*,fmt='(a)') 'You have defined a soil sand fraction < 0. The DEFAULT SOIL'
      write (unit=*,fmt='(a)') 'PARAMETERS WILL BE USED based on NL%NLSCON.                '
      write (unit=*,fmt='(a)') '==========================================================='
   end if
   
   if (isoilflg(ifm)==2 .and. slxsand>1.) then
      write (unit=*,fmt='(a)') '==========================================================='
      write (unit=*,fmt='(a)') '   WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!   '
      write (unit=*,fmt='(a)') '==========================================================='
      write (unit=*,fmt='(a)') 'You have defined a soil sand fraction > 1. The DEFAULT SOIL'
      write (unit=*,fmt='(a)') 'PARAMETERS WILL BE USED based on NL%NLSCON.                '
      write (unit=*,fmt='(a)') '==========================================================='
   end if
   

   if (isoilflg(ifm)==2 .and. (slxsand+slxclay)>1.) then
      write (unit=*,fmt='(a)') '==========================================================='
      write (unit=*,fmt='(a)') '   WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!   '
      write (unit=*,fmt='(a)') '==========================================================='
      write (unit=*,fmt='(a)') 'Your soil and clay fractions add up to more than 1!  The   '
      write (unit=*,fmt='(a)') 'DEFAULT SOIL PARAMETERS WILL BE USED based on NL%NLSCON.   '
      write (unit=*,fmt='(a)') '==========================================================='
   end if
end do
   
#if defined(COUPLED)
   if (isoilstateinit < 0 .or. isoilstateinit > 2) then
      write (reason,fmt='(a,1x,i4,a)')                                                     &
             'Invalid ISOILSTATEINIT, it must be between 0 and 2. Yours is set to'         &
             ,isoilstateinit,'...'
      call opspec_fatal(reason,'opspec_misc')
      ifaterr = ifaterr +1
   end if
#else
   if (isoilstateinit < 0 .or. isoilstateinit > 1) then
      write (reason,fmt='(a,1x,i4,a)')                                                     &
             'Invalid ISOILSTATEINIT, it must be between 0 and 1. Yours is set to'         &
             ,isoilstateinit,'...'
      call opspec_fatal(reason,'opspec_misc')
      ifaterr = ifaterr +1
   end if
#endif


   if (isoildepthflg < 0 .or. isoildepthflg > 1) then
      write (reason,fmt='(a,1x,i4,a)')                                                     &
              'Invalid ISOILDEPTHFLG, it must be between 0 and 1. Yours is set to'         &
              ,isoildepthflg,'...'
      call opspec_fatal(reason,'opspec_misc')
      ifaterr = ifaterr +1
   end if

   if (isoilbc < 0 .or. isoilbc > 3) then
      write (reason,fmt='(a,1x,i4,a)')                                                     &
        'Invalid ISOILBC, it must be between 0 and 3.  Yours is set to',isoilbc,'...'
      call opspec_fatal(reason,'opspec_misc')
      ifaterr = ifaterr +1
   else if(isoilbc == 2 .and. (sldrain < 0. .or. sldrain > 90.)) then
      write (reason,fmt='(a,1x,es12.5,a)')                                                 &
        'Invalid SLDRAIN, it must be between 0 and 90.  Yours is set to',sldrain,'...'
      call opspec_fatal(reason,'opspec_misc')
      ifaterr = ifaterr +1
   end if

   if (ivegt_dynamics < 0 .or. ivegt_dynamics > 1) then
      write (reason,fmt='(a,1x,i4,a)')                                                     &
         'Invalid IVEGT_DYNAMICS, it must be between 0 and 1. Yours is set to'              &
        ,ivegt_dynamics,'...'
      call opspec_fatal(reason,'opspec_misc')
      ifaterr = ifaterr +1
   end if
   
   if (growth_resp_scheme < 0 .or. growth_resp_scheme > 1) then
      write (reason,fmt='(a,1x,i4,a)')                                                     &
         'Invalid GROWTH_RESP_SCHEME, it must be 0 or 1. Yours is set to'                  &
        ,growth_resp_scheme,'...'
      call opspec_fatal(reason,'opspec_misc')
      ifaterr = ifaterr +1
   end if
   
   if (storage_resp_scheme < 0 .or. storage_resp_scheme > 1) then
      write (reason,fmt='(a,1x,i4,a)')                                                     &
         'Invalid STORAGE_RESP_SCHEME, it must be 0 or 1. Yours is set to'                 &
        ,storage_resp_scheme,'...'
      call opspec_fatal(reason,'opspec_misc')
      ifaterr = ifaterr +1
   end if

   if (ibigleaf < 0 .or. ibigleaf > 1) then
      write (reason,fmt='(a,1x,i4,a)')                                                     &
         'Invalid IBIGLEAF, it must be between 0 and 1. Yours is set to',ibigleaf,'...'
      call opspec_fatal(reason,'opspec_misc')
      ifaterr = ifaterr +1
   elseif (ibigleaf == 1 .and. ( crown_mod /= 0 .or. abs(maxpatch) == 1 )) then
      write (reason,fmt='(a,1x,i4,a)')                                                     &
         'CROWN_MOD must be 0 and MAXPATCH cannot be -1 or 1 when IBIGLEAF is set to 1...'
      call opspec_fatal(reason,'opspec_misc')
      ifaterr = ifaterr +1
   end if

   !---------------------------------------------------------------------------------------!
   !      Integration scheme can be only 0 (Euler) or 1 (4th order Runge-Kutta).  The      !
   ! branch thermodynamics is currently working only with Runge-Kutta, so we won't allow   !
   ! using it in case the user decides for Euler.                                          !
   !---------------------------------------------------------------------------------------!
   select case (integration_scheme)
   case (0:3)
      !------------------------------------------------------------------------------------!
      !   Check the branch thermodynamics.                                                 !
      !------------------------------------------------------------------------------------!
      if (ibranch_thermo < 0 .or. ibranch_thermo > 2) then
         write (reason,fmt='(a,1x,i4,a)')                                                  &
                   'Invalid IBRANCH_THERMO, it must be between 0 and 2. Yours is set to'   &
                   ,ibranch_thermo,'...'
         call opspec_fatal(reason,'opspec_misc')
         ifaterr = ifaterr +1
      end if
      
      !------------------------------------------------------------------------------------!
      !   Check the Runge-Kutta tolerance.                                                 !
      !------------------------------------------------------------------------------------!
      if (rk4_tolerance < 1.e-7 .or. rk4_tolerance > 1.e-1) then
         write (reason,fmt='(a,2x,a,1x,es12.5,a)')                                         & 
              'Invalid RK4_TOLERANCE, it must be between 1.e-7 and 1.e-1.'                 &
             ,'Yours is set to',rk4_tolerance,'...'
         call opspec_fatal(reason,'opspec_misc')
         ifaterr = ifaterr +1
      end if
   case default
      write (reason,fmt='(a,1x,i4,a)')                                                     &
               'Invalid INTEGRATION_SCHEME, it must be 0, 1, 2, or 3. Yours is set to'     &
               ,integration_scheme,'...'
      call opspec_fatal(reason,'opspec_misc')
      ifaterr = ifaterr +1
   end select
   !---------------------------------------------------------------------------------------!

   if (iphysiol < 0 .or. iphysiol > 3) then
      write (reason,fmt='(a,1x,i4,a)')                                                     &
                    'Invalid IPHYSIOL, it must be between 0 and 3. Yours is set to'        &
                    ,iphysiol,'...'
      call opspec_fatal(reason,'opspec_misc')
      ifaterr = ifaterr +1
   end if

   if (iallom < 0 .or. iallom > 4) then
      write (reason,fmt='(a,1x,i4,a)')                                                     &
                    'Invalid IALLOM, it must be between 0 and 4. Yours is set to'          &
                    ,iallom,'...'
      call opspec_fatal(reason,'opspec_misc')
      ifaterr = ifaterr +1
   end if

   if (igrass < 0 .or. igrass > 1) then
      write (reason,fmt='(a,1x,i4,a)')                                                     &
                    'Invalid IGRASS, it must be between 0 and 1. Yours is set to'          &
                    ,igrass,'...'
      call opspec_fatal(reason,'opspec_misc')
      ifaterr = ifaterr +1
   elseif (igrass == 1 .and. ibigleaf == 1) then
      write (reason,fmt='(a,1x,a)')                                                        &
                    'Invalid setting.  New grass allometry (IGRASS = 1) is not supported'  &
                   ,'in big leaf ED (IBIGLEAF = 1)...'
      call opspec_fatal(reason,'opspec_misc')
      ifaterr = ifaterr +1
   
   end if

   if (iphen_scheme < -1 .or. iphen_scheme > 3) then
      write (reason,fmt='(a,1x,i4,a)')                                                     &
                    'Invalid IPHEN_SCHEME, it must be between -1 and 3. Yours is set to'   &
                    ,iphen_scheme,'...'
      call opspec_fatal(reason,'opspec_misc')
      ifaterr = ifaterr +1
   end if

   if (repro_scheme < 0 .or. repro_scheme > 3) then
      write (reason,fmt='(a,1x,i4,a)')                                                     &
                    'Invalid REPRO_SCHEME, it must be between  0 and 3. Yours is set to'   &
                    ,repro_scheme,'...'
      call opspec_fatal(reason,'opspec_misc')
      ifaterr = ifaterr +1
   end if

   if (radint < -100.0 .or. radint > 100.0) then
      write (reason,fmt='(a,1x,es12.5,a)')                                                 &
                    'Invalid RADINT, it must be between -100 and 100. Yours is set to'     &
                    ,radint,'...'
      call opspec_fatal(reason,'opspec_misc')
      ifaterr = ifaterr +1
   end if
 
    if (radslp < 0.0 .or. radslp > 1.0) then
      write (reason,fmt='(a,1x,es12.5,a)')                                                 &
                    'Invalid RADSLP, it must be between 0 and 1. Yours is set to'          &
                    ,radslp,'...'
      call opspec_fatal(reason,'opspec_misc')
      ifaterr = ifaterr +1
   end if  
   if (h2o_plant_lim < 0 .or. h2o_plant_lim > 2) then
      write (reason,fmt='(a,1x,i4,a)')                                                     &
                    'Invalid H2O_PLANT_LIM, it must be between 0 and 2.  Yours is set to'  &
                    ,h2o_plant_lim,'...'
      call opspec_fatal(reason,'opspec_misc')
      ifaterr = ifaterr +1
   end if

   if (iddmort_scheme < 0 .or. iddmort_scheme > 1) then
      write (reason,fmt='(a,1x,i4,a)')                                                     &
                    'Invalid IDDMORT_SCHEME, it must be between 0 and 1.  Yours is set to' &
                    ,iddmort_scheme,'...'
      call opspec_fatal(reason,'opspec_misc')
      ifaterr = ifaterr +1
   end if

   if (cbr_scheme < 0 .or. cbr_scheme > 2) then
      write (reason,fmt='(a,1x,i4,a)')                                                     &
                    'Invalid CBR_SCHEME, it must range from 0 to 2.  Yours is set to' &
                    ,cbr_scheme,'...'
      call opspec_fatal(reason,'opspec_misc')
      ifaterr = ifaterr +1
   end if

   if (ddmort_const < 0. .or. ddmort_const > 1.) then
      write (reason,fmt='(a,1x,es12.5,a)')                                                 &
                    'Invalid DDMORT_CONST, it must be between 0 and 1.  Yours is set to'   &
                    ,ddmort_const,'...'
      call opspec_fatal(reason,'opspec_misc')
      ifaterr = ifaterr +1
   end if

   if (vmfact_c3 < 0.01 .or. vmfact_c3 > 100.) then
      write (reason,fmt='(a,1x,es12.5,a)')                                                 &
                    'Invalid VMFACT_C3, it must be between 0.01 and 100.  Yours is set to' &
                    ,vmfact_c3,'...'
      call opspec_fatal(reason,'opspec_misc')
      ifaterr = ifaterr +1
   end if

   if (vmfact_c4 < 0.01 .or. vmfact_c4 > 100.) then
      write (reason,fmt='(a,1x,es12.5,a)')                                                 &
                    'Invalid VMFACT_C4, it must be between 0.01 and 100.  Yours is set to' &
                    ,vmfact_c4,'...'
      call opspec_fatal(reason,'opspec_misc')
      ifaterr = ifaterr +1
   end if

   if (mphoto_trc3 < 0.1 .or. mphoto_trc3 > 20.) then
      write (reason,fmt='(a,1x,es12.5,a)')                                                 &
                    'Invalid MPHOTO_TRC3, it must be between 0.1 and 20.  Yours is set to' &
                    ,mphoto_trc3,'...'
      call opspec_fatal(reason,'opspec_misc')
      ifaterr = ifaterr +1
   end if

   if (mphoto_tec3 < 0.1 .or. mphoto_tec3 > 20.) then
      write (reason,fmt='(a,1x,es12.5,a)')                                                 &
                    'Invalid MPHOTO_TEC3, it must be between 0.1 and 20.  Yours is set to' &
                    ,mphoto_trc3,'...'
      call opspec_fatal(reason,'opspec_misc')
      ifaterr = ifaterr +1
   end if

   if (mphoto_c4 < 0.1 .or. mphoto_c4 > 20.) then
      write (reason,fmt='(a,1x,es12.5,a)')                                                 &
                    'Invalid MPHOTO_C4, it must be between 0.1 and 20.  Yours is set to'   &
                    ,mphoto_c4,'...'
      call opspec_fatal(reason,'opspec_misc')
      ifaterr = ifaterr +1
   end if
 
   if (bphoto_blc3 < 500. .or. bphoto_blc3 > 50000.) then
      write (reason,fmt='(a,1x,es12.5,a)')                                                 &
                'Invalid BPHOTO_BLC3, it must be between 500. and 50000.  Yours is set to' &
               ,bphoto_blc3,'...'
      call opspec_fatal(reason,'opspec_misc')
      ifaterr = ifaterr +1
   end if
 
   if (bphoto_nlc3 < 500. .or. bphoto_nlc3 > 50000.) then
      write (reason,fmt='(a,1x,es12.5,a)')                                                 &
                'Invalid BPHOTO_NLC3, it must be between 500. and 50000.  Yours is set to' &
               ,bphoto_nlc3,'...'
      call opspec_fatal(reason,'opspec_misc')
      ifaterr = ifaterr +1
   end if
 
   if (bphoto_c4 < 500. .or. bphoto_c4 > 50000.) then
      write (reason,fmt='(a,1x,es12.5,a)')                                                 &
                'Invalid BPHOTO_C4, it must be between 500. and 50000.  Yours is set to'   &
               ,bphoto_c4,'...'
      call opspec_fatal(reason,'opspec_misc')
      ifaterr = ifaterr +1
   end if
  
  if (kw_grass < 15. .or. kw_grass > 15000.) then
      write (reason,fmt='(a,1x,es12.5,a)')                                                 &
                    'Invalid KW_GRASS, it must be between 15 and 15000.  Yours is set to'  &
                    ,kw_grass,'...'
      call opspec_fatal(reason,'opspec_misc')
      ifaterr = ifaterr +1
   end if
   
  if (kw_tree < 15. .or. kw_tree > 15000.) then
      write (reason,fmt='(a,1x,es12.5,a)')                                                 &
                    'Invalid KW_TREE, it must be between 15 and 15000.  Yours is set to'   &
                    ,kw_tree,'...'
      call opspec_fatal(reason,'opspec_misc')
      ifaterr = ifaterr +1
   end if

   if (gamma_c3 < 0.001 .or. gamma_c3 > 0.10) then
      write (reason,fmt='(a,1x,es12.5,a)')                                                 &
                    'Invalid GAMMA_C3, it must be between 0.001 and 0.1.  Yours is set to' &
                    ,gamma_c3,'...'
      call opspec_fatal(reason,'opspec_misc')
      ifaterr = ifaterr +1
   end if

   if (gamma_c4 < 0.001 .or. gamma_c4 > 0.10) then
      write (reason,fmt='(a,1x,es12.5,a)')                                                 &
                    'Invalid GAMMA_C4, it must be between 0.001 and 0.1.  Yours is set to' &
                    ,gamma_c4,'...'
      call opspec_fatal(reason,'opspec_misc')
      ifaterr = ifaterr +1
   end if

   if (d0_grass < 0.01 .or. d0_grass > 1.) then
      write (reason,fmt='(a,1x,es12.5,a)')                                                 &
                    'Invalid D0_GRASS, it must be between 0.01 and 1.  Yours is set to'    &
                    ,d0_grass,'...'
      call opspec_fatal(reason,'opspec_misc')
      ifaterr = ifaterr +1
   end if

   if (d0_tree < 0.01 .or. d0_tree > 1.) then
      write (reason,fmt='(a,1x,es12.5,a)')                                                 &
                    'Invalid D0_TREE, it must be between 0.01 and 1. Yours is set to'      &
                    ,d0_tree,'...'
      call opspec_fatal(reason,'opspec_misc')
      ifaterr = ifaterr +1
   end if
   
   if (alpha_c3 < 0.001 .or. alpha_c3 > 1.) then
      write (reason,fmt='(a,1x,es12.5,a)')                                                 &
                    'Invalid ALPHA_C3, it must be between 0.001 and 1.  Yours is set to'   &
                    ,alpha_c3,'...'
      call opspec_fatal(reason,'opspec_misc')
      ifaterr = ifaterr +1
   end if

   if (alpha_c4 < 0.001 .or. alpha_c4 > 1.) then
      write (reason,fmt='(a,1x,es12.5,a)')                                                 &
                    'Invalid ALPHA_C4, it must be between 0.001 and 1.  Yours is set to'   &
                    ,alpha_c4,'...'
      call opspec_fatal(reason,'opspec_misc')
      ifaterr = ifaterr +1
   end if

   if (klowco2in < 10. .or. klowco2in > 1000000.) then
      write (reason,fmt='(a,1x,es12.5,a)')                                                 &
                    'Invalid KLOWCO2IN, it must be between 10. and 1.e6.  Yours is set to' &
                    ,klowco2in,'...'
      call opspec_fatal(reason,'opspec_misc')
      ifaterr = ifaterr +1
   end if

   if (rrffact < 0.1 .or. rrffact > 10.) then
      write (reason,fmt='(a,1x,es12.5,a)')                                                 &
                    'Invalid RRFFACT, it must be between 0.1 and 10. Yours is set to'      &
                    ,rrffact,'...'
      call opspec_fatal(reason,'opspec_misc')
      ifaterr = ifaterr +1
   end if
     
   if (growthresp < 0.0 .or. growthresp > 1.0) then
      write (reason,fmt='(a,1x,es12.5,a)')                                                 &
                    'Invalid GROWTHRESP, it must be between 0 and 1. Yours is set to'      &
                    ,growthresp,'...'
      call opspec_fatal(reason,'opspec_misc')
      ifaterr = ifaterr +1
   end if

   if (lwidth_grass < 0.01 .or. lwidth_grass > 0.30) then
      write (reason,fmt='(a,1x,es12.5,a)')                                                 &
            'Invalid LWIDTH_GRASS, it must be between 0.01 and 0.30 Yours is set to'       &
           ,lwidth_grass,'...'
      call opspec_fatal(reason,'opspec_misc')
      ifaterr = ifaterr +1
   end if

   if (lwidth_bltree < 0.01 .or. lwidth_bltree > 0.30) then
      write (reason,fmt='(a,1x,es12.5,a)')                                                 &
            'Invalid LWIDTH_BLTREE, it must be between 0.01 and 0.30 Yours is set to'      &
           ,lwidth_bltree,'...'
      call opspec_fatal(reason,'opspec_misc')
      ifaterr = ifaterr +1
   end if

   if (lwidth_nltree < 0.01 .or. lwidth_nltree > 0.30) then
      write (reason,fmt='(a,1x,es12.5,a)')                                                 &
            'Invalid LWIDTH_NLTREE, it must be between 0.01 and 0.30 Yours is set to'      &
           ,lwidth_nltree,'...'
      call opspec_fatal(reason,'opspec_misc')
      ifaterr = ifaterr +1
   end if
   
   if (q10_c3 < 1.0 .or. q10_c3 > 10.) then
      write (reason,fmt='(a,1x,es12.5,a)')                                                 &
                    'Invalid Q10_C3, it must be between 1.0 and 10.  Yours is set to'      &
                    ,q10_c3,'...'
      call opspec_fatal(reason,'opspec_misc')
      ifaterr = ifaterr +1
   end if
   
   if (q10_c4 < 1.0 .or. q10_c4 > 10.) then
      write (reason,fmt='(a,1x,es12.5,a)')                                                 &
                    'Invalid Q10_C4, it must be between 1.0 and 10.  Yours is set to'      &
                    ,q10_c4,'...'
      call opspec_fatal(reason,'opspec_misc')
      ifaterr = ifaterr +1
   end if

   if (thetacrit < -1.49 .or. thetacrit > 1.) then
      write (reason,fmt='(a,1x,es12.5,a)')                                                 &
                    'Invalid THETACRIT, it must be between -1.49 and 1. Yours is set to'   &
                    ,thetacrit,'...'
      call opspec_fatal(reason,'opspec_misc')
      ifaterr = ifaterr +1
   end if

   if (quantum_efficiency_T < 0 .or. quantum_efficiency_T > 1) then
      write (reason,fmt='(a,1x,i4,a)')                                                     &
                    'Invalid QUANTUM_EFFICIENCY_T, it must be either 0 or 1. Yours is set to'         &
                    ,quantum_efficiency_T,'...'
      call opspec_fatal(reason,'opspec_misc')
      ifaterr = ifaterr +1
   end if

   if (n_plant_lim < 0 .or. n_plant_lim > 1) then
      write (reason,fmt='(a,1x,i4,a)')                                                     &
                    'Invalid N_PLANT_LIM, it must be between 0 and 1. Yours is set to'     &
                    ,n_plant_lim,'...'
      call opspec_fatal(reason,'opspec_misc')
      ifaterr = ifaterr +1
   end if

   if (n_decomp_lim < 0 .or. n_decomp_lim > 1) then
      write (reason,fmt='(a,1x,i4,a)')                                                     &
                    'Invalid N_DECOMP_LIM, it must be between 0 and 1. Yours is set to'    &
                    ,n_decomp_lim,'...'
      call opspec_fatal(reason,'opspec_misc')
      ifaterr = ifaterr +1
   end if

   if (decomp_scheme < 0 .or. decomp_scheme > 2) then
      write (reason,fmt='(a,1x,i4,a)')                                                     &
                    'Invalid DECOMP_SCHEME, it must be between 0 and 2. Yours is set to'   &
                    ,decomp_scheme,'...'
      call opspec_fatal(reason,'opspec_misc')
      ifaterr = ifaterr +1
   end if

   if (include_fire < 0 .or. include_fire > 3) then
      write (reason,fmt='(a,1x,i4,a)')                                                     &
                    'Invalid INCLUDE_FIRE, it must be between 0 and 3. Yours is set to'    &
                    ,include_fire,'...'
      call opspec_fatal(reason,'opspec_misc')
      ifaterr = ifaterr +1
   else if (include_fire /= 0) then
      if (fire_parameter < 0.0 .or. fire_parameter > 100.) then
         write (reason,fmt='(a,1x,es12.5,a)')                                              &
               'Invalid FIRE_PARAMETER, it must be between 0 and 100.. Yours is set to'    &
             , fire_parameter,'...'
         call opspec_fatal(reason,'opspec_misc')
         ifaterr = ifaterr +1
      end if
   end if

   select case (include_fire)
   case(3)
      if (sm_fire < 0. .or. sm_fire > 2) then
         write (reason,fmt='(a,1x,a,1x,i4,a,1x,es12.5,a)')                                 &
                        'Invalid SM_FIRE, it must be between 0 and 2'                      &
                       ,'when INLCUDE_FIRE is ',include_fire,'.  Your SM_FIRE is set to'   &
                       ,sm_fire,'...'
         call opspec_fatal(reason,'opspec_misc')
         ifaterr = ifaterr +1
      end if
   case(0:2)
      if (sm_fire < -3.1 .or. sm_fire > 1.) then
         write (reason,fmt='(a,1x,a,1x,i4,a,1x,es12.5,a)')                                 &
                        'Invalid SM_FIRE, it must be between -3.1 and 1.0'                 &
                       ,'when INLCUDE_FIRE is ',include_fire,'.  Your SM_FIRE is set to'   &
                       ,sm_fire,'...'
         call opspec_fatal(reason,'opspec_misc')
         ifaterr = ifaterr +1
      end if
   end select

   if (ianth_disturb < 0 .or. ianth_disturb > 1) then
      write (reason,fmt='(a,1x,i4,a)') &
        'Invalid IANTH_DISTURB, it must be between 0 and 1. Yours is set to',ianth_disturb,'...'
      call opspec_fatal(reason,'opspec_misc')  
      ifaterr = ifaterr +1
   end if

   select case (icanturb)
   case (0:4)
      continue
   case default
      write (reason,fmt='(a,1x,i4,a)') &
        'Invalid ICANTURB, it must be between 0 and 4. Yours is set to',icanturb,'...'
      call opspec_fatal(reason,'opspec_misc')  
      ifaterr = ifaterr +1
   end select

   if (isfclyrm < 0 .or. isfclyrm > 4) then
      write (reason,fmt='(a,1x,i4,a)') &
        'Invalid ISFCLYRM, it must be between 0 and 4. Yours is set to',isfclyrm,'...'
      call opspec_fatal(reason,'opspec_misc')  
      ifaterr = ifaterr +1
   end if

   if (ied_grndvap < 0 .or. ied_grndvap > 5) then
      write (reason,fmt='(a,1x,i4,a)') &
        'Invalid IED_GRNDVAP, it must be between 0 and 5.  Yours is set to',ied_grndvap    &
       ,'...'
      call opspec_fatal(reason,'opspec_misc')  
      ifaterr = ifaterr +1
   end if

   if (ipercol < 0 .or. ipercol > 2) then
      write (reason,fmt='(a,1x,i4,a)')                                                     &
        'Invalid IPERCOL, it must be between 0 and 2. Yours is set to',ipercol,'...'
      call opspec_fatal(reason,'opspec_misc')  
      ifaterr = ifaterr +1
   end if
   
   ! Checking if I am attempting to include invalid pfts. I can leave the loop once I hit
   ! the first "huge" value
   pftloop: do ipft=1,n_pft
      if (include_these_pft(ipft) == skip) then
        if (ipft /= 1) then
           exit pftloop
        else
           write (reason,fmt='(a)') &
              'You did not specify any valid INCLUDE_THESE_PFT, the run must have at least one pft...'
           call opspec_fatal(reason,'opspec_misc')
           ifaterr=ifaterr+1
        end if
      elseif (include_these_pft(ipft) < 0 .or. include_these_pft(ipft) > n_pft) then
         write (reason,fmt='(a,1x,i4,a,1x,i4,a)') &
            'Invalid INCLUDE_THESE_PFT, it must be between 1 and ',n_pft, &
            '. One of yours is set to',include_these_pft(ipft),'...'
         call opspec_fatal(reason,'opspec_misc')  
         ifaterr=ifaterr+1
      end if
   end do pftloop
   
   !----- Checking whether the user choice for agriculture and plantation make sense. -----!
   if (ianth_disturb == 1) then
      !------ Checking the plantation PFT.  It must be a tree PFT. ------------------------!
      select case (plantation_stock)
      case (2,3,4,6,7,8,9,10,11,17)
         continue
      case default
         write(reason,fmt='(a,1x,i5,a)')                                                   &
                      'Invalid plantation_stock , it can''t be grass and yours is set to'  &
                      ,plantation_stock,'...'
         ifaterr = ifaterr +1
         call opspec_fatal(reason,'opspec_misc')
      end select
   
      !------ Checking the plantation PFT. It must be a grass PFT. ------------------------!
      select case (agri_stock)
      case (1,5,12,13,14,15,16)
         continue
      case default
         write(reason,fmt='(a,1x,i5,a)')                                                   &
            'Invalid AGRI_STOCK , it can''t be a tree and yours is set to',agri_stock,'...'
         ifaterr = ifaterr +1
         call opspec_fatal(reason,'opspec_misc')
      end select

      agri_ok       = .false.
      plantation_ok = .false.
      agriloop: do ipft=1,n_pft
         if (include_these_pft(ipft) == skip) exit agriloop
         if (agri_stock == include_these_pft(ipft)) agri_ok = .true.
         if (plantation_stock == include_these_pft(ipft)) plantation_ok=.true.
      end do agriloop
      
      if (.not. agri_ok) then
         write(reason,fmt='(a,1x,i5,a)')                                                   &
            'Invalid AGRI_STOCK (',agri_stock,'). The pft must be in INCLUDE_THESE_PFT.'
         ifaterr = ifaterr +1
         call opspec_fatal(reason,'opspec_misc')
      end if
      
      if (.not. plantation_ok) then
         write(reason,fmt='(a,1x,i5,a)')                                                   &
            'Invalid PLANTATION_STOCK (',plantation_stock,                                 &
             &'). The pft must be in INCLUDE_THESE_PFT.'
         ifaterr = ifaterr +1
         call opspec_fatal(reason,'opspec_misc')
      end if
   end if

   if (pft_1st_check < 0 .or. pft_1st_check > 2) then
      write (reason,fmt='(a,1x,i4,a)')                                                     &
                    'Invalid PFT_1ST_CHECK, it must be between 0 and 2.  Yours is set to'  &
                    ,pft_1st_check,'...'
      ifaterr = ifaterr +1
      call opspec_fatal(reason,'opspec_misc')
   end if

   if (crown_mod < 0 .or. crown_mod > 1) then
      write (reason,fmt='(a,1x,i4,a)')                                                     &
                    'Invalid CROWN_MOD, it must be either 0 or 1.  Yours is set to'        &
                    ,crown_mod,'...'
      ifaterr = ifaterr +1
      call opspec_fatal(reason,'opspec_misc')
   end if

   if  (icanrad < 1 .or. icanrad > 2) then
      write (reason,fmt='(a,1x,i4,a)')                                                     &
                    'Invalid ICANRAD, it must be between 1 and 2.  Yours is set to'        &
                    ,icanrad,'...'
      ifaterr = ifaterr +1
      call opspec_fatal(reason,'opspec_misc')
   elseif (icanrad /= 0 .and. crown_mod == 1) then
      write(reason,fmt='(a)') 'CROWN_MOD must be turned off when ICANRAD is 1 or 2...'
      ifaterr = ifaterr +1
      call opspec_fatal(reason,'opspec_misc')
   end if

   if  (ihrzrad < 0 .or. ihrzrad > 4) then
      write (reason,fmt='(a,1x,i4,a)')                                                     &
                    'Invalid IHRZRAD, it must be between 0 and 4.  Yours is set to'        &
                    ,ihrzrad,'...'
      ifaterr = ifaterr +1
      call opspec_fatal(reason,'opspec_misc')
   end if

   if  (ltrans_vis < 0.01 .or. ltrans_vis > 0.99) then
      write (reason,fmt='(a,2x,a,1x,es12.5,a)')                                            &
                    'Invalid LTRANS_VIS, it must be between 0.01 and 0.99.'                &
                   ,'Yours is set to',ltrans_vis,'...'
      ifaterr = ifaterr +1
      call opspec_fatal(reason,'opspec_misc')
   end if

   if  (ltrans_nir < 0.01 .or. ltrans_nir > 0.99) then
      write (reason,fmt='(a,2x,a,1x,es12.5,a)')                                            &
                    'Invalid LTRANS_NIR, it must be between 0.01 and 0.99.'                &
                   ,'Yours is set to',ltrans_nir,'...'
      ifaterr = ifaterr +1
      call opspec_fatal(reason,'opspec_misc')
   end if

   if  (lreflect_vis < 0.01 .or. lreflect_vis > 0.99) then
      write (reason,fmt='(a,2x,a,1x,es12.5,a)')                                            &
                    'Invalid LREFLECT_VIS, it must be between 0.01 and 0.99.'              &
                   ,'Yours is set to',lreflect_vis,'...'
      ifaterr = ifaterr +1
      call opspec_fatal(reason,'opspec_misc')
   end if

   if  (lreflect_nir < 0.01 .or. lreflect_nir > 0.99) then
      write (reason,fmt='(a,2x,a,1x,es12.5,a)')                                            &
                    'Invalid LREFLECT_NIR, it must be between 0.01 and 0.99.'              &
                   ,'Yours is set to',lreflect_nir,'...'
      ifaterr = ifaterr +1
      call opspec_fatal(reason,'opspec_misc')
   end if

   if  (lreflect_vis + ltrans_vis > 0.99) then
      write (unit=*,fmt='(a,1x,es12.5)') ' LTRANS_VIS   = ',ltrans_vis
      write (unit=*,fmt='(a,1x,es12.5)') ' LREFLECT_VIS = ',lreflect_vis
      write (unit=*,fmt='(a,1x,es12.5)') ' LABSORPT_VIS = ',1. - ltrans_vis - lreflect_vis
      write (reason,fmt='(a,2x,a)')                                                        &
                    'LTRANS_VIS + LREFLECT_VIS cannot exceed 0.99.'                        &
                   ,'This causes absorptance to be weird (a bad thing).'
      ifaterr = ifaterr +1
      call opspec_fatal(reason,'opspec_misc')
   end if

   if  (lreflect_nir + ltrans_nir > 0.99) then
      write (unit=*,fmt='(a,1x,es12.5)') ' LTRANS_NIR   = ',ltrans_nir
      write (unit=*,fmt='(a,1x,es12.5)') ' LREFLECT_NIR = ',lreflect_nir
      write (unit=*,fmt='(a,1x,es12.5)') ' LABSORPT_NIR = ',1. - ltrans_nir - lreflect_nir
      write (reason,fmt='(a,2x,a)')                                                        &
                    'LTRANS_NIR + LREFLECT_NIR cannot exceed 0.99.'                        &
                   ,'This causes absorptance to be weird (a bad thing).'
      ifaterr = ifaterr +1
      call opspec_fatal(reason,'opspec_misc')
   end if

   if  (orient_tree < -0.40 .or. orient_tree > 0.60) then
      write (reason,fmt='(a,2x,a,1x,es12.5,a)')                                            &
                    'Invalid ORIENT_TREE, it must be between -0.40 and 0.60.'              &
                   ,'Yours is set to',orient_tree,'...'
      ifaterr = ifaterr +1
      call opspec_fatal(reason,'opspec_misc')
   end if

   if  (orient_grass < -0.40 .or. orient_grass > 0.60) then
      write (reason,fmt='(a,2x,a,1x,es12.5,a)')                                            &
                    'Invalid ORIENT_GRASS, it must be between -0.40 and 0.60.'             &
                   ,'Yours is set to',orient_grass,'...'
      ifaterr = ifaterr +1
      call opspec_fatal(reason,'opspec_misc')
   end if

   if  (clump_tree < 0.01 .or. clump_tree > 1.00) then
      write (reason,fmt='(a,2x,a,1x,es12.5,a)')                                            &
                    'Invalid CLUMP_TREE, it must be between 0.01 and 1.00.'                &
                   ,'Yours is set to',clump_tree,'...'
      ifaterr = ifaterr +1
      call opspec_fatal(reason,'opspec_misc')
   end if

   if  (clump_grass < 0.01 .or. clump_grass > 1.00) then
      write (reason,fmt='(a,2x,a,1x,es12.5,a)')                                            &
                    'Invalid CLUMP_GRASS, it must be between 0.01 and 1.00.'               &
                   ,'Yours is set to',clump_grass,'...'
      ifaterr = ifaterr +1
      call opspec_fatal(reason,'opspec_misc')
   end if

   if (ihrzrad /= 0 .and. ihrzrad /= 4 .and. (igoutput < 0 .or. igoutput > 1)) then
      write (reason,fmt='(a,1x,i4,a)')                                                     &
                    'Invalid IGOUTPUT, it must be 0 or 1.  Yours is set to'                &
                    ,igoutput,'...'
      ifaterr = ifaterr +1
      call opspec_fatal(reason,'opspec_misc')
   end if


    
   if (zrough <= 0.0) then
      write (reason,fmt='(a,1x,es14.7,a)')                                                 &
                    'Invalid ZROUGH, it must be positive.  Yours is set to',zrough,'...'
      call opspec_fatal(reason,'opspec_misc')
      ifaterr = ifaterr +1
   end if

   if (treefall_disturbance_rate /= 0.0) then
      if (time2canopy < 0.0) then
         write (reason,fmt='(a,1x,es14.7,a)')                                              &
                'Invalid TIME2CANOPY, it can''t be negative.  Yours is set to'             &
                ,time2canopy,'...'
         call opspec_fatal(reason,'opspec_misc')
         ifaterr = ifaterr +1
      elseif (time2canopy >= 2.0 * (1. - epsilon(1.)) / abs(treefall_disturbance_rate))    &
      then
         !---------------------------------------------------------------------------------!
         !     We need this if statement because the effective treefall tends to infinity  !
         ! as time2canopy approaches 2. / treefall_disturbance_rate.                       !
         !---------------------------------------------------------------------------------!
         write (unit=*,fmt='(a,1x,es14.7)')   ' TREEFALL_DISTURBANCE_RATE ='               &
                                             ,treefall_disturbance_rate
         write (unit=*,fmt='(a,1x,es14.7)')   ' MAX(TIME2CANOPY)          ='               &
                                             , 2.0 * (1. - epsilon(1.))                    &
                                             / abs(treefall_disturbance_rate)
         write (reason,fmt='(a,2x,a,1x,es14.7,a)')                                         &
             ' Invalid TIME2CANOPY, it can''t be >= 2. / TREEFALL_DISTURBANCE_RATE.'       &
            ,'Yours is set to',time2canopy,'...'
         call opspec_fatal(reason,'opspec_misc')  
         ifaterr = ifaterr +1
      end if
   end if
    
   
    
   if (runoff_time < 0.0) then
      write (reason,fmt='(a,1x,es14.7,a)')                                                 &
            'Invalid RUNOFF_TIME, it can''t be negative. Yours is set to',runoff_time,'...'
      call opspec_fatal(reason,'opspec_misc')  
      ifaterr = ifaterr +1
   end if
    
   if (ubmin < 0.0001 .or. ubmin > 2.0) then
      write (reason,fmt='(a,1x,es14.7,a)')                                                 &
            'Invalid UBMIN, it must be between 0.0001 and 2.0.  Yours is set to'           &
           ,ustmin,'...'
      call opspec_fatal(reason,'opspec_misc')  
      ifaterr = ifaterr +1
   end if
    
   if (ustmin < 0.0001 .or. ustmin > 1.0) then
      write (reason,fmt='(a,1x,es14.7,a)')                                                 &
            'Invalid USTMIN, it must be between 0.0001 and 1.0. Yours is set to'           &
           ,ustmin,'...'
      call opspec_fatal(reason,'opspec_misc')  
      ifaterr = ifaterr +1
   elseif (ustmin > ubmin) then
      write (unit=*,fmt='(a,1x,es12.5)') ' UBMIN  = ',ubmin
      write (unit=*,fmt='(a,1x,es12.5)') ' UGBMIN = ',ugbmin
      write (unit=*,fmt='(a,1x,es12.5)') ' USTMIN = ',ustmin
      write (reason,fmt='(a)') 'Invalid USTMIN, it can''t be greater than UBMIN...'
      call opspec_fatal(reason,'opspec_misc')  
      ifaterr = ifaterr +1
   end if

   if (ugbmin < ustmin .or. ugbmin > ubmin) then
      write (unit=*,fmt='(a,1x,es12.5)') ' UBMIN  = ',ubmin
      write (unit=*,fmt='(a,1x,es12.5)') ' UGBMIN = ',ugbmin
      write (unit=*,fmt='(a,1x,es12.5)') ' USTMIN = ',ustmin
      write (reason,fmt='(a)') 'Invalid UGBMIN, it can''t be between USTMIN and UBMIN...'
      call opspec_fatal(reason,'opspec_misc')  
      ifaterr = ifaterr +1
   end if

   if (gamm < 0.1 .or. gamm > 100.0) then
      write (reason,fmt='(a,1x,es14.7,a)')                                                 &
            'Invalid GAMM, it must be between 0.1 and 100.0. Yours is set to'              &
           ,gamm,'...'
      call opspec_fatal(reason,'opspec_misc')  
      ifaterr = ifaterr +1
   end if

   if (gamh < 0.1 .or. gamh > 100.0) then
      write (reason,fmt='(a,1x,es14.7,a)')                                                 &
            'Invalid GAMH, it must be between 0.1 and 100.0. Yours is set to'              &
           ,gamh,'...'
      call opspec_fatal(reason,'opspec_misc')  
      ifaterr = ifaterr +1
   end if

   if (tprandtl < 0.01 .or. tprandtl > 100.0) then
      write (reason,fmt='(a,1x,es14.7,a)')                                                 &
            'Invalid TPRANDTL, it must be between 0.01 and 100.0. Yours is set to'         &
           ,tprandtl,'...'
      call opspec_fatal(reason,'opspec_misc')  
      ifaterr = ifaterr +1
   end if

   if (ribmax < 0.01 .or. ribmax > 1.0) then
      write (reason,fmt='(a,1x,es14.7,a)')                                                 &
            'Invalid RIBMAX, it must be between 0.01 and 1.0..  Yours is set to'           &
           ,ribmax,'...'
      call opspec_fatal(reason,'opspec_misc')  
      ifaterr = ifaterr +1
   end if

   if (leaf_maxwhc < 0.0 .or. leaf_maxwhc > 10.) then
      write (reason,fmt='(a,1x,es14.7,a)')                                                 &
            'Invalid LEAF_MAXWHC, it must be between 0.0 and 10..  Yours is set to'        &
           ,leaf_maxwhc,'...'
      call opspec_fatal(reason,'opspec_misc')  
      ifaterr = ifaterr +1
   end if

#if defined(COUPLED)
#else
   if (ishuffle < 0 .or. ishuffle > 2) then
      write (reason,fmt='(a,1x,i4,a)')                                                     &
                    'Invalid ISHUFFLE, it must be between 0 and 2.  Yours is set to'       &
                    ,ishuffle,'...'
      ifaterr = ifaterr +1
      call opspec_fatal(reason,'opspec_misc')
   end if
   if (imetavg < -1 .or. imetavg > 3) then
      write (reason,fmt='(a,1x,i4,a)')                                                     &
                    'Invalid IMETAVG, it must be between -1 and 3.  Yours is set to'       &
                    ,imetavg,'...'
      ifaterr = ifaterr +1
      call opspec_fatal(reason,'opspec_misc')
   end if

#endif

   if (imetrad < 0 .or. imetrad > 5) then
      write (reason,fmt='(a,1x,i4,a)')                                                     &
                    'Invalid IMETRAD, it must be between 0 and 5.  Yours is set to'        &
                    ,imetrad,'...'
      ifaterr = ifaterr +1
      call opspec_fatal(reason,'opspec_misc')
   end if

   if (dt_census < 1 .or. dt_census > 120) then
      write (reason,fmt='(a,1x,i4,a)')                                                     &
                    'Invalid DT_CENSUS, it must be between 1 and 120.  Yours is set to'    &
                    ,dt_census,'...'
      ifaterr = ifaterr +1
      call opspec_fatal(reason,'opspec_misc')
   end if

   if (yr1st_census < 1200 .or. yr1st_census > 3200) then
      write (reason,fmt='(a,1x,i4,a)')                                                     &
             'Invalid YR1ST_CENSUS, it must be between 1200 and 3200.  Yours is set to'    &
            ,yr1st_census,'...'
      ifaterr = ifaterr +1
      call opspec_fatal(reason,'opspec_misc')
   end if

   if (mon1st_census < 1 .or. mon1st_census > 12) then
      write (reason,fmt='(a,1x,i4,a)')                                                     &
             'Invalid MON1ST_CENSUS, it must be between 1 and 12.  Yours is set to'        &
            ,mon1st_census,'...'
      ifaterr = ifaterr +1
      call opspec_fatal(reason,'opspec_misc')
   end if

   if (min_recruit_dbh < 0 .or. min_recruit_dbh > 100.) then
      write (reason,fmt='(a,1x,es12.5,a)')                                                 &
             'Invalid MIN_RECRUIT_DBH, it must be between 0 and 100.  Yours is set to'     &
            ,min_recruit_dbh,'...'
      ifaterr = ifaterr +1
      call opspec_fatal(reason,'opspec_misc')
   end if


   if (idetailed < 0 .or. idetailed > 63) then
      write (reason,fmt='(a,1x,i4,a)')                                                     &
                    'Invalid IDETAILED, it must be between 0 and 63.  Yours is set to'     &
                    ,idetailed,'...'
      ifaterr = ifaterr +1
      call opspec_fatal(reason,'opspec_misc')
   elseif (idetailed > 0) then
      patch_detailed = ibclr(idetailed,5) > 0

      if (patch_detailed .and. (n_poi > 1 .or. n_ed_region > 0)) then
         write(unit=*,fmt='(a)')       '--------------------------------------------------'
         write(unit=*,fmt='(a,1x,i6)') ' IDETAILED   = ',idetailed
         write(unit=*,fmt='(a,1x,i6)') ' N_POI       = ',n_poi
         write(unit=*,fmt='(a,1x,i6)') ' N_ED_REGION = ',n_ed_region
         write(unit=*,fmt='(a)')       ' '
         write(unit=*,fmt='(a)')       ' The following should be all F in regional runs'
         write(unit=*,fmt='(a)')       '    or multiple polygon runs'
         write(unit=*,fmt='(a,1x,l1)') ' - BUDGET              [ 1] = ',btest(idetailed,0)
         write(unit=*,fmt='(a,1x,l1)') ' - PHOTOSYNTHESIS      [ 2] = ',btest(idetailed,1)
         write(unit=*,fmt='(a,1x,l1)') ' - DETAILED INTEGRATOR [ 4] = ',btest(idetailed,2)
         write(unit=*,fmt='(a,1x,l1)') ' - SANITY CHECK BOUNDS [ 8] = ',btest(idetailed,3)
         write(unit=*,fmt='(a,1x,l1)') ' - ERROR RECORDING     [16] = ',btest(idetailed,4)
         write(unit=*,fmt='(a)')       '--------------------------------------------------'
         write (reason,fmt='(2(a,1x))') 'Only single polygon runs are allowed'             &
                                       ,'with detailed patch-level output...'
         ifaterr = ifaterr +1
         call opspec_fatal(reason,'opspec_misc')
      end if

      if (patch_keep < -2) then
         write (reason,fmt='(a,2x,a,1x,i4,a)')                                             &
                       'Invalid PATCH_KEEP, it must be between -2 and number of patches.'  &
                      ,'Yours is set to',patch_keep,'...'
         ifaterr = ifaterr +1
         call opspec_fatal(reason,'opspec_misc')
      end if
      !------------------------------------------------------------------------------------!
   end if
   !---------------------------------------------------------------------------------------!




   !----- Stop the run if there are any fatal errors. -------------------------------------!
   if (ifaterr > 0) then
      write (unit=*,fmt='(a)')       ' -----------ED_OPSPEC_MISC --------------------------'
      write (unit=*,fmt='(a,1x,i5)') ' Fatal errors:',ifaterr
      write (unit=*,fmt='(a)')       ' ----------------------------------------------------'
      call fatal_error('Fatal errors at namelist - Misc settings '                         &
                      ,'ed_opspec_misc','ed_opspec.f90')
   end if

   return
end subroutine ed_opspec_misc
!==========================================================================================!
!==========================================================================================!
