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
subroutine ed_opspec_grid
!------------------------------------------------------------------------------------------!
!   This subroutine checks the grid definition for ED, like grid size, dimensions and if   !
! the provided numbers are allowed.                                                        !
!------------------------------------------------------------------------------------------!
  use max_dims, only : maxgrds,nxpmax,nypmax,nzpmax,nzgmax,nzsmax,max_soi,max_ed_regions
  use grid_coms, only : nnxp,nnyp,ngrids,polelat,polelon,centlat,centlon,deltax,deltay    &
                       ,nstratx,nstraty,nzg,nzs
  use mem_sites, only : n_ed_region,n_soi,grid_type,grid_res,soi_lat,soi_lon              &
                       ,ed_reg_latmin,ed_reg_latmax,ed_reg_lonmin,ed_reg_lonmax
  use soil_coms, only : slz

  implicit none

  integer :: ifaterr, ifm,k
  character(len=222)           :: reason

  ifaterr=0

  ! Check whether the number of grids is okay.
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
  endif

  ! Check whether you are trying to run SOI and regional simultaneously. 
  ! This is currently forbidden.
  if ((.not. (n_soi > 0 .and. n_ed_region == 0)) .and. (.not. (n_soi == 0 .and. n_ed_region > 0))) then
     write(reason,'(a,1x,i4,1x,a,1x,i4,1x,a)')  &
          'One of n_soi or n_ed_region needs to be zero.( Yours: n_soi='                   &
          ,n_soi,', n_ed_region=',n_ed_region,').'
     call opspec_fatal(reason,'opspec_grid')  
     ifaterr=ifaterr+1
  elseif (n_soi < 0) then
     write(reason,'(a,1x,i4,a)') &
        'N_SOI needs to be non-negative. Yours is currently set to ',n_soi,'.'
     call opspec_fatal(reason,'opspec_grid')  
     ifaterr=ifaterr+1
  elseif (n_ed_region < 0) then
     write(reason,'(a,1x,i4,a)') &
        'N_ED_REGION needs to be non-negative. Yours is currently set to ',n_ed_region,'.'
     call opspec_fatal(reason,'opspec_grid')  
     ifaterr=ifaterr+1
  elseif (n_soi > max_soi) then
     write(reason,'(a,1x,i4,a,1x,i4,a)') &
        'N_SOI cannot be greater than',max_soi,'. Yours is currently set to ',n_soi,'.'
     call opspec_fatal(reason,'opspec_grid')  
     ifaterr=ifaterr+1
  elseif (n_ed_region> max_ed_regions) then
     write(reason,'(a,1x,i4,a,1x,i4,a)') &
        'N_ED_REGION cannot be greater than',max_ed_regions,'. Yours is currently set to ',n_ed_region,'.'
     call opspec_fatal(reason,'opspec_grid')  
     ifaterr=ifaterr+1
  end if

  ! Check whether the number of grid points makes sense for ED regions
  do ifm=1,n_ed_region
     if (nnxp(ifm) > nxpmax) then
        write(reason,'(3(a,1x,i4,a))')  &
             'Number of x-points cannot be greater than ',nxpmax,'.',' (Yours is',nnxp(ifm),' ','in grid',ifm,').'
        call opspec_fatal(reason,'opspec_grid')  
        ifaterr=ifaterr+1
     end if
     if (nnyp(ifm) > nypmax) then
        write(reason,'(3(a,1x,i4,a))')  &
             'Number of y-points cannot be greater than ',nypmax,'.',' (Yours is',nnyp(ifm),' ','in grid',ifm,').'
        call opspec_fatal(reason,'opspec_grid')  
        ifaterr=ifaterr+1
     end if
     if (nstratx(ifm) < 1) then
        write (reason,'(a,1x,i4,1x,a,1x,i4,a)') &
              'Nest x ratio should be at least 1. Your nstratx for grid',ifm,'is currently set to',nstratx(ifm),'.'
        call opspec_fatal(reason,'opspec_grid')  
        ifaterr=ifaterr+1        
     endif
     if (nstraty(ifm) < 1) then
        write (reason,'(a,1x,i4,1x,a,1x,i4,a)') &
              'Nest y ratio should be at least 1. Your nstraty for grid',ifm,'is currently set to',nstratx(ifm),'.'
        call opspec_fatal(reason,'opspec_grid')  
        ifaterr=ifaterr+1        
     endif
  end do

  !  Check whether ED has allowable values to define the grids
  if (n_ed_region > 0) then
     if (grid_type < 0 .or. grid_type > 1) then
        write (reason,'(a,1x,i4,a)') &
              'Grid_type should be either 0 or 1 for regional runs. Yours is currently set to',grid_type,'.'
        call opspec_fatal(reason,'opspec_grid')  
        ifaterr=ifaterr+1
     end if
     if (grid_type == 0) then
        do ifm=1,n_ed_region
           if (ed_reg_latmin(ifm) < -90. .or. ed_reg_latmin(ifm) > 90. ) then
              write (reason,'(a,1x,i4,a,1x,es14.7,a)') &
                    'ED_REG_LATMIN is outside [-90.;90.] for region ',ifm,'. Yours is currently set to',ed_reg_latmin(ifm),'...'
              call opspec_fatal(reason,'opspec_grid')  
              ifaterr=ifaterr+1
           end if
           if (ed_reg_latmax(ifm) < -90. .or. ed_reg_latmax(ifm) > 90. ) then
              write (reason,'(a,1x,i4,a,1x,es14.7,a)') &
                    'ED_REG_LATMAX is outside [-90.;90.] for region ',ifm,'. Yours is currently set to',ed_reg_latmax(ifm),'...'
              call opspec_fatal(reason,'opspec_grid')  
              ifaterr=ifaterr+1
           end if
           if (ed_reg_lonmin(ifm) < -180. .or. ed_reg_lonmin(ifm) > 180. ) then
              write (reason,'(a,1x,i4,a,1x,es14.7,a)') &
                    'ED_REG_LONMIN is outside [-180.;180.] for region ',ifm,'. Yours is currently set to',ed_reg_latmin(ifm),'...'
              call opspec_fatal(reason,'opspec_grid')  
              ifaterr=ifaterr+1
           end if
           if (ed_reg_lonmax(ifm) < -180. .or. ed_reg_lonmax(ifm) > 180. ) then
              write (reason,'(a,1x,i4,a,1x,es14.7,a)') &
                    'ED_REG_LONMAX is outside [-180.;180.] for region ',ifm,'. Yours is currently set to',ed_reg_latmax(ifm),'...'
              call opspec_fatal(reason,'opspec_grid')  
              ifaterr=ifaterr+1
           end if
        end do
        if (grid_res <= 0.) then
           write (reason,'(a,1x,es14.7,a)') &
                 'GRID_RES must be positive. Yours is currently set to',grid_res,'...'
           call opspec_fatal(reason,'opspec_grid')  
           ifaterr=ifaterr+1
        end if
     elseif (grid_type == 1) then
        if (deltax <= 0.) then
           write (reason,'(a,1x,es14.7,a)') &
                 'DELTAX must be positive. Yours is currently set to',deltax,'...'
           call opspec_fatal(reason,'opspec_grid')  
           ifaterr=ifaterr+1
        end if
        if (deltay <= 0.) then
           write (reason,'(a,1x,es14.7,a)') &
                 'DELTAY must be positive. Yours is currently set to',deltax,'...'
           call opspec_fatal(reason,'opspec_grid')  
           ifaterr=ifaterr+1
        end if
        if (polelat < -90. .or. polelat > 90.) then
           write (reason,'(a,1x,i4,a)') &
                 'POLELAT is outside [-90.;90.].',' Yours is currently set to',polelat,'...'
           call opspec_fatal(reason,'opspec_grid')  
           ifaterr=ifaterr+1
        end if
        if (polelon < -180. .or. polelon > 180.) then
           write (reason,'(a,1x,i4,a)') &
                 'POLELON is outside [-180.;180.].',' Yours is currently set to',polelon,'...'
           call opspec_fatal(reason,'opspec_grid')  
           ifaterr=ifaterr+1
        end if
        do ifm=1,n_ed_region
           if (centlat(ifm) < -90. .or. centlat(ifm) > 90. ) then
              write (reason,'(a,1x,i4,a,1x,es14.7,a)') &
                    'CENTLAT is outside [-90.;90.] for region ',ifm,'. Yours is currently set to',centlat(ifm),'...'
              call opspec_fatal(reason,'opspec_grid')  
              ifaterr=ifaterr+1
           end if
           if (centlon(ifm) < -180. .or. centlon(ifm) > 180. ) then
              write (reason,'(a,1x,i4,a,1x,es14.7,a)') &
                    'CENTLON is outside [-180.;180.] for region ',ifm,'. Yours is currently set to',centlon(ifm),'...'
              call opspec_fatal(reason,'opspec_grid')  
              ifaterr=ifaterr+1
           end if
        end do
     end if
  end if

  
  if (n_soi > 0)then
     do ifm=1,n_soi
        if (soi_lat(ifm) < -90. .or. soi_lat(ifm) > 90. ) then
           write (reason,'(a,1x,i4,a,1x,es14.7,a)') &
                 'SOI_LAT is outside [-90.;90.] for SOI #',ifm,'. Yours is currently set to',soi_lat(ifm),'...'
           call opspec_fatal(reason,'opspec_grid')  
           ifaterr=ifaterr+1
        end if
        if (soi_lon(ifm) < -180. .or. soi_lon(ifm) > 180. ) then
           write (reason,'(a,1x,i4,a,1x,es14.7,a)') &
                 'SOI_LON is outside [-180.;180.] for SOI #',ifm,'. Yours is currently set to',soi_lon(ifm),'...'
           call opspec_fatal(reason,'opspec_grid')  
           ifaterr=ifaterr+1
        end if
     end do

  end if

  ! Check whether ED soil is okay
  ! check whether the number of soil levels is outside the allowed range
  if (nzg < 2) then
     write (reason,'(a,1x,i4,a)') &
           'Too few soil layers. Set it to at least 2. Your nzg is currently set to',nzg,'.'
     call opspec_fatal(reason,'opspec_grid')  
     ifaterr=ifaterr+1        
  elseif (nzg > nzgmax) then 
     write (reason,'(2(a,1x,i5,a))') &
           'The number of soil layers cannot be greater than ',nzgmax,'.',' Your nzg is currently set to',nzg,'.'
     call opspec_fatal(reason,'opspec_grid') 
     ifaterr=ifaterr+1 
  end if
  do k=1,nzg
     if (slz(k) > -.001) then
        write (reason,'(a,1x,i4,1x,a,1x,es14.7,a)') &
              'Your soil level #',k,'is not enough below ground. It is currently set to'   &
              ,slz(k),', make it deeper than -0.001...'
        call opspec_fatal(reason,'opspec_grid')  
        ifaterr=ifaterr+1        
     endif
  enddo

  do k=1,nzg-1
     if (slz(k)-slz(k+1) > .001) then
        write (reason,'(2(a,1x,i4,1x),a,1x,es14.7,1x,a,1x,es14.7,a)') &
              'Soil layers #',k,'and',k+1,'are not enough apart (i.e. > 0.001) . They are currently set as '   &
              ,slz(k),'and',slz(k+1),'...'
        call opspec_fatal(reason,'opspec_grid')  
        ifaterr=ifaterr+1        
     endif
  enddo


  ! Check whether ED snow layers is okay
  ! check whether the number of soil levels is outside the allowed range
  if (nzs < 1) then
     write (reason,'(a,1x,i4,a)') &
           'Too few maximum # of snow layers. Set it to at least 1. Your nzs is currently set to',nzs,'.'
     call opspec_fatal(reason,'opspec_grid')  
     ifaterr=ifaterr+1        
  elseif (nzs > nzsmax) then 
     write (reason,'(2(a,1x,i5,a))') &
           'The number of snow layers cannot be greater than ',nzsmax,'.',' Your nzs is currently set to',nzs,'.'
     call opspec_fatal(reason,'opspec_grid') 
     ifaterr=ifaterr+1 
  end if


  ! stop the run if there are any fatal errors.
  if (ifaterr > 0) then
     write (unit=*,fmt='(a)')       ' -----------ED_OPSPEC_GRID --------------------------'
     write (unit=*,fmt='(a,1x,i5)') ' Fatal errors:',ifaterr
     write (unit=*,fmt='(a)')       ' ----------------------------------------------------'
     call fatal_error('Fatal errors at namelist - Grid settings '&
                     ,'ed_opspec_grid','ed_opspec.f90')
  end if
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
   use max_dims     , only : maxmach
   use grid_coms    , only : nnxp,nnyp,ngrids
   use mem_sites    , only : n_ed_region,n_soi
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
   endif
   do ifm=1,ngrids
      if (n_ed_region > 0 .and. nnxp(ifm)*nnyp(ifm) < (nmachs +1)) then
         write (reason,'(3(a,1x,i5,1x),a)')                                                 &
               'Region number',ifm,'is too small. You have only',nnxp(ifm)*nnyp(ifm)        &
              ,'points for',nmachs+1,'nodes.'
         call opspec_fatal(reason,'opspec_par')  
         ifaterr=ifaterr+1
      end if
   end do
   ! Checking if the user is trying to run SOI in parallel 
   ! (this will be allowed at some point).
   if (n_soi > 0 .and. nmachs > 0) then
      write (reason,'(3(a,1x,i5,1x),a)')                                                    &
        'You are attempting to run SOI runs in parallel, this is not supported currently...'
      call opspec_fatal(reason,'opspec_par')  
      ifaterr=ifaterr+1
   end if

   ! stop the run if there are any fatal errors.
   if (ifaterr > 0) then
      write (unit=*,fmt='(a)')       ' -----------ED_OPSPEC_GRID --------------------------'
      write (unit=*,fmt='(a,1x,i5)') ' Fatal errors:',ifaterr
      write (unit=*,fmt='(a)')       ' ----------------------------------------------------'
      call fatal_error('Fatal errors at namelist - Grid settings '&
                      ,'ed_opspec_grid','ed_opspec.f90')
   end if

   return
end subroutine ed_opspec_par
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine ed_opspec_times
!------------------------------------------------------------------------------------------!
!    This subroutine checks time related settings from ED2IN.                              !
!------------------------------------------------------------------------------------------!
   use misc_coms , only : frqfast,frqstate,imontha,idatea,iyeara,itimea  &
                         ,imonthz,idatez,iyearz,itimez,dtlsm,radfrq      &
                         ,ifoutput,isoutput,idoutput,imoutput,iyoutput,  &
                         nrec_fast,nrec_state,outfast,outstate,unitfast,unitstate
   use consts_coms, only : day_sec,hr_sec
   use grid_coms , only : timmax


   implicit none
   character(len=222)           :: reason
   integer :: ifaterr,n
   ifaterr=0




   !---------------------------------------------------------------------------------------!
   !     Checking whether the user provided a valid combination of unitfast, frqfast, and  !
   ! outfast.                                                                              !
   !---------------------------------------------------------------------------------------!
   select case(unitfast)
   !---------------------------------------------------------------------------------------!
   !    Seconds                                                                            !
   !---------------------------------------------------------------------------------------!
   case (0)
      if (ifoutput == 0) then
         nrec_fast = 1 ! Useless, no output will be created.
      !------------------------------------------------------------------------------------!
      !    If unitfrq is 0, then frqfast must be either a divisor or a multiple of one day !
      ! in case there is daily and/or monthly analysis.                                    !
      !------------------------------------------------------------------------------------!
      elseif ((imoutput /= 0 .or. idoutput /= 0 .or. outfast == -1. .or. outfast == -2. )  &
             .and. (mod(day_sec,frqfast) /= 0. .and. mod(frqfast,day_sec) /= 0.)) then
         write(reason,fmt='(a,1x,f8.2,1x,a,1x,es14.7)')  &
             'FRQFAST must be a divisor or an integer multiple of ',day_sec,               &
             ' sec when daily and/or monthly analysis are on. Yours is set to ',frqfast
         call opspec_fatal(reason,'opspec_times')
         ifaterr = ifaterr + 1
      !------------------------------------------------------------------------------------!
      !   Checking outfast setting and adjusting it if needed. Depending on the            !
      ! configuration, this will cause the run to crash.                                   !
      !------------------------------------------------------------------------------------!
      !----- User didn't specify any outfast, use frqfast ---------------------------------!
      elseif(outfast == 0.) then
         outfast    = frqfast
         nrec_fast  = 1
      elseif (outfast == -1.) then
         outfast = day_sec
         nrec_fast  = outfast/frqfast
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
         nrec_fast = outfast/frqfast
      end if

   !---------------------------------------------------------------------------------------!
   !    Days.                                                                              !
   !---------------------------------------------------------------------------------------!
   case (1)
      if (ifoutput == 0) then
         nrec_fast = 1 ! Useless, no output will be created.
      !------------------------------------------------------------------------------------!
      !    If unitfrq is 0, then frqfast must be either a divisor or a multiple of one day !
      ! in case there is daily and/or monthly analysis.                                    !
      !------------------------------------------------------------------------------------!
      elseif ((imoutput /= 0 .or. idoutput /= 0) .and.                                     &
              (mod(day_sec,day_sec*frqfast) /= 0. .and. mod(frqfast,1.) /= 0.)) then
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
         nrec_fast  = outfast/frqfast
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
         nrec_fast = outfast/frqfast
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
      !------------------------------------------------------------------------------------!
      !    If unifrq is monthly, it needs to be a round number.                            !
      !------------------------------------------------------------------------------------!
      elseif (mod(frqfast,1.) /= 0) then
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
      elseif (mod(12.,frqfast) /= 0) then
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
      elseif (mod(frqfast,1.) /= 0) then
         write(reason,fmt='(a,1x,f8.2,1x,a,1x,es14.7)')                                    &
             'FRQFAST must be a round number when unitfast is ',unitfast,                  &
             '(years). Yours is currently set to ',frqfast
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
         nrec_state  = outstate/frqstate
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
         nrec_state = outstate/frqstate
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
         nrec_state  = outstate/frqstate
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
         nrec_state = outstate/frqstate
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
      elseif (mod(frqstate,1.) /= 0) then
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
      elseif (mod(12.,frqstate) /= 0) then
         write(reason,fmt='(a,1x,f8.2,1x,a,1x,es14.7)')                                    &
             'FRQSTATE must be a divisor of 12 (one year) when unitstate is ',unitstate,   &
             '(months). Yours is currently set to ',frqstate
         call opspec_fatal(reason,'opspec_times')  
         ifaterr = ifaterr + 1
      !------------------------------------------------------------------------------------!
      !    This is fine but now outstate must be set exactly as frqstate. If the user wasn't !
      ! aware of this, print an informative banner.                                        !
      !------------------------------------------------------------------------------------!
      elseif (outstate /= 0. .or. outstate > frqstate) then
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
      elseif (mod(frqstate,1.) /= 0) then
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
         outstate = frqstate
         nrec_state = 1
         write (unit=*,fmt='(a)') ' '
         write (unit=*,fmt='(a)') '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         write (unit=*,fmt='(a)') '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         write (unit=*,fmt='(a)') '!!!  WARNING! WARNING! WARNING! WARNING! WARNING!   !!!'
         write (unit=*,fmt='(a)') '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         write (unit=*,fmt='(a)') '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         write (unit=*,fmt='(a)') ' -> Outstate cannot be different than frqstate when'
         write (unit=*,fmt='(a)') '    unitstate is set to 3 (years).'
         write (unit=*,fmt='(a,1x,f7.0,1x,a)')                                             &
                                  '    Oustate was redefined to ',outstate,'years.'
         write (unit=*,fmt='(a)') '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         write (unit=*,fmt='(a)') '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         write (unit=*,fmt='(a)') ' '
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
   if (isoutput /= 0. .and. ifoutput /= 0.) then
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
            end if
         case (2,3) !---- State in months or years, frqfast must be divisor of 1 day ------!
            if (mod(frqstate*day_sec,frqfast) /= 0.0) then
               write(reason,fmt='(a,1x,a,1x,f10.2,1x,a,1x,f10.2)')                         &
                    'FRQFAST must be a divisor of one day when frqfast is in seconds and', &
                    'FRQSTATE is in months or years. Your FRQFAST=',frqfast,'sec.'
            end if
         end select

      case (1) !----- Days ----------------------------------------------------------------!
         select case(unitstate)
         case (0) !---- state in seconds --------------------------------------------------!
            if (mod(frqstate,frqfast*day_sec) /= 0.0) then
               write(reason,fmt='(a,1x,a,1x,f10.2,1x,a,1x,f10.2,1x,a)')                    &
                    'FRQFAST must be a divisor of FRQSTATE if both are outputing data.',   &
                    'Yours is set to FRQFAST=',frqfast,'days and FRQSTATE=',frqstate,'sec!'
            end if

         case (2,3) !---- State in months or years, frqfast must be divisor of 1 day ------!
            if (frqfast /= 1.0) then
               write(reason,fmt='(a,1x,a,1x,f10.2,1x,a,1x,f10.2,1x,a)')                    &
                    'FRQFAST must be a divisor of one day or 1 day when UNITFAST is in',   &
                    'days and FRQSTATE is in months or years. Your FRQFAST=',frqfast,'days.'
            end if

         end select

      case (2) !----- Months --------------------------------------------------------------!
         select case(unitstate)
         case (0,1) !---- state in seconds or days. Can't be since months are irregular. --!
            write(reason,fmt='(a,1x,a,1x,i5)')                                             &
               'If UNITFAST is in months, UNITSTATE must be either months or year, and',   &
               'currently UNITSTATE=',unitstate
      
         case (3) !---- State years, frqfast must be divisor ------------------------------!
            if (mod(frqstate*12,frqfast) /= 0.) then
               write(reason,fmt='(a,1x,2(a,1x,f10.2,1x),a)')                               &
                    'FRQFAST must be a divisor of FRQSTATE if both are outputing data.',   &
                    'Yours is set to FRQFAST=',frqfast,'mon and FRQSTATE=',frqstate,'yrs!'
            end if

         end select

      case (3) !----- Years ---------------------------------------------------------------!
         select case(unitstate)
         case (0,1) !---- state in seconds or days. Can't be since months are irregular. --!
            write(reason,fmt='(a,1x,a,1x,i5)')                                             &
               'If UNITFAST is in months, UNITSTATE must be either months or year, and',   &
               'currently UNITSTATE=',unitstate
      
         case (2) !---- State months, frqfast must be divisor -----------------------------!
            if (mod(frqstate,frqfast*12) /= 0.) then
               write(reason,fmt='(a,1x,2(a,1x,f10.2,1x),a)')                               &
                    'FRQFAST must be a divisor of FRQSTATE if both are outputing data.',   &
                    'Yours is set to FRQFAST=',frqfast,'yrs and FRQSTATE=',frqstate,'mon!'
            end if
         end select
      end select
   end if


   !Check if this simulation has a positive timmax.
   if (timmax < 0.0) then
      write(reason,fmt='(a,2(a,1x,2(i2.2,a),i4.4,1x,i4.4,a,1x))')     &
         'Your end time is before the initial time.'                  &
         ,'Initial:',imontha,'/',idatea,'/',iyeara,itimea,'GMT'       &
         ,'Final  :',imonthz,'/',idatez,'/',iyearz,itimez,'GMT'
      call opspec_fatal(reason,'opspec_times')  
      ifaterr=ifaterr+1
   end if

   if (mod(hr_sec,dtlsm) /= 0.0) then
      write(reason,fmt='(a,1x,f8.2,1x,a,1x,es14.7)')  &
          'DTLSM must be a divisor of ',hr_sec,' sec. Yours is set to ',dtlsm
      call opspec_fatal(reason,'opspec_times')  
      ifaterr=ifaterr+1
   end if

   if (mod(hr_sec,radfrq) /= 0.0) then
      write(reason,fmt='(a,1x,f8.2,1x,a,1x,es14.7)')  &
          'RADFRQ must be a divisor of ',hr_sec,' sec. Yours is set to ',radfrq
      call opspec_fatal(reason,'opspec_times')  
      ifaterr=ifaterr+1
   end if

   ! Not sure if this is really necessary. If not, please remove it...
   if (mod(radfrq,dtlsm) /= 0.0) then
      write(reason,fmt='(a,1x,f8.2,1x,a,1x,f8.2)')  &
          'DTLSM must be a divisor of RADFRQ. Your DTLSM is set to',dtlsm, &
          'and your RADFRQ is set to',radfrq,'...'
      call opspec_fatal(reason,'opspec_times')  
      ifaterr=ifaterr+1
   end if

   ! Stop the run if there are any fatal errors.
   if (ifaterr > 0) then
      write (unit=*,fmt='(a)')       ' -----------ED_OPSPEC_TIMES -------------------------'
      write (unit=*,fmt='(a,1x,i5)') ' Fatal errors:',ifaterr
      write (unit=*,fmt='(a)')       ' ----------------------------------------------------'
      call fatal_error('Fatal errors at namelist - Time settings '&
                      ,'ed_opspec_times','ed_opspec.f90')
   end if
   return
end subroutine ed_opspec_times
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine ed_opspec_misc
!------------------------------------------------------------------------------------------!
!   This subroutine performs miscellaneous tests over the options, like values outside the !
! allowed ranges and conflicting dynamic settings.
!------------------------------------------------------------------------------------------!
   use max_dims, only : n_pft
   use misc_coms, only : ifoutput,idoutput,imoutput,iyoutput,isoutput,iclobber,runtype,ied_init_mode &
                        ,integration_scheme
   use soil_coms, only : isoilflg, nslcon,isoilstateinit,isoildepthflg,zrough,runoff_time
   use mem_sites, only : n_soi,n_ed_region,maxpatch,maxcohort
   use grid_coms , only : ngrids
   use physiology_coms, only : istoma_scheme,n_plant_lim
   use decomp_coms, only : n_decomp_lim
   use disturb_coms, only : include_fire,ianth_disturb,treefall_disturbance_rate
   use phenology_coms, only : iphen_scheme
   use pft_coms, only  : include_these_pft,pft_1st_check

   implicit none
   character(len=222)           :: reason
   integer, parameter           :: skip=huge(6)
   integer                      :: ifaterr,ifm,ipft
   ifaterr=0

   if (ifoutput /= 0 .and. ifoutput /= 3) then
      write (reason,fmt='(a,1x,i4,a)') &
        'Invalid IFOUTPUT, it must be 0 (none) or 3 (HDF5). Yours is set to',ifoutput,'...'
      call opspec_fatal(reason,'opspec_misc')  
      ifaterr = ifaterr +1
   end if
   if (idoutput /= 0 .and. idoutput /= 3) then
      write (reason,fmt='(a,1x,i4,a)') &
        'Invalid IDOUTPUT, it must be 0 (none) or 3 (HDF5). Yours is set to',idoutput,'...'
      call opspec_fatal(reason,'opspec_misc')  
      ifaterr = ifaterr +1
   end if
   if (imoutput /= 0 .and. imoutput /= 3) then
      write (reason,fmt='(a,1x,i4,a)') &
        'Invalid IMOUTPUT, it must be 0 (none) or 3 (HDF5). Yours is set to',imoutput,'...'
      call opspec_fatal(reason,'opspec_misc')  
      ifaterr = ifaterr +1
   end if
   if (iyoutput /= 0 .and. iyoutput /= 3) then
      write (reason,fmt='(a,1x,i4,a)') &
        'Invalid IYOUTPUT, it must be 0 (none) or 3 (HDF5). Yours is set to',iyoutput,'...'
      call opspec_fatal(reason,'opspec_misc')  
      ifaterr = ifaterr +1
   end if
   if (isoutput /= 0 .and. isoutput /= 3) then
      write (reason,fmt='(a,1x,i4,a)') &
        'Invalid ISOUTPUT, it must be 0 (none) or 3 (HDF5). Yours is set to',isoutput,'...'
      call opspec_fatal(reason,'opspec_misc')  
      ifaterr = ifaterr +1
   end if
   if (iclobber < 0 .or. iclobber > 1) then
      write (reason,fmt='(a,1x,i4,a)') &
        'Invalid ICLOBBER, it must be 0 or 1. Yours is set to',iclobber,'...'
      call opspec_fatal(reason,'opspec_misc')  
      ifaterr = ifaterr +1
   end if
   if (trim(runtype) /= 'INITIAL' .and. trim(runtype) /= 'HISTORY') then
      write (reason,fmt='(a,1x,2a)') &
        'Invalid RUNTYPE, it must be INITIAL or HISTORY. Yours is set to',trim(runtype),'...'
      call opspec_fatal(reason,'opspec_misc')  
      ifaterr = ifaterr +1
   end if

   if (ied_init_mode == -1) then ! This should be avoided. Use as a last resort!
      write (unit=*,fmt='(a)') '==========================================================='
      write (unit=*,fmt='(a)') '   WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!   '
      write (unit=*,fmt='(a)') '   WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!   '
      write (unit=*,fmt='(a)') '   WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!   '
      write (unit=*,fmt='(a)') '   WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!   '
      write (unit=*,fmt='(a)') '   WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!   '
      write (unit=*,fmt='(a)') '   WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!   '
      write (unit=*,fmt='(a)') '==========================================================='
      write (unit=*,fmt='(a)') '    You have set up the run to read ED-1 and ED-2 files    '
      write (unit=*,fmt='(a)') ' mixed in the same run. This is against my beliefs and you '
      write (unit=*,fmt='(a)') ' should know that I''m only doing this because I am a very '
      write (unit=*,fmt='(a)') ' nice model and I don''t know how to say NO! to such a     '
      write (unit=*,fmt='(a)') ' desperate user. But don''t put high expectations on this  '
      write (unit=*,fmt='(a)') ' run, and in case it crashes, it is going to be all your   '
      write (unit=*,fmt='(a)') ' fault and I will remind you that!!!                       '
      write (unit=*,fmt='(a)') '==========================================================='
   elseif (ied_init_mode < 0 .or. ied_init_mode > 3) then
      write (reason,fmt='(a,1x,i4,a)') &
        'Invalid IED_INIT_MODE, it must be between 0 and 3. Yours is set to',ied_init_mode,'...'
      call opspec_fatal(reason,'opspec_misc')  
      ifaterr = ifaterr +1
   end if

   do ifm=1,ngrids
      if (isoilflg(ifm) < 1 .or. isoilflg(ifm) > 2) then
         write (reason,fmt='(a,1x,i4,1x,a,1x,i4,a)') &
           'Invalid ISOILFLG, it must be between 0 and 3. Yours is set to',isoilflg(ifm),'for grid',ifm,'...'
         call opspec_fatal(reason,'opspec_misc')  
         ifaterr = ifaterr +1
      end if
   end do

   if (nslcon < 1 .or. nslcon > 12) then
      write (reason,fmt='(a,1x,i4,a)') &
        'Invalid NSLCON, it must be between 1 and 12. Yours is set to',nslcon,'...'
      call opspec_fatal(reason,'opspec_misc')  
      ifaterr = ifaterr +1
   end if

   if (isoilstateinit < 0 .or. isoilstateinit > 1) then
      write (reason,fmt='(a,1x,i4,a)') &
        'Invalid ISOILSTATEINIT, it must be between 0 and 1. Yours is set to',isoilstateinit,'...'
      call opspec_fatal(reason,'opspec_misc')  
      ifaterr = ifaterr +1
   end if

   if (isoildepthflg < 0 .or. isoildepthflg > 1) then
      write (reason,fmt='(a,1x,i4,a)') &
        'Invalid ISOILDEPTHFLG, it must be between 0 and 1. Yours is set to',isoildepthflg,'...'
      call opspec_fatal(reason,'opspec_misc')  
      ifaterr = ifaterr +1
   end if

   if (integration_scheme < 0 .or. integration_scheme > 1) then
      write (reason,fmt='(a,1x,i4,a)') &
        'Invalid INTEGRATION_SCHEME, it must be between 0 and 1. Yours is set to',integration_scheme,'...'
      call opspec_fatal(reason,'opspec_misc')  
      ifaterr = ifaterr +1
   end if

   if (istoma_scheme < 0 .or. istoma_scheme > 1) then
      write (reason,fmt='(a,1x,i4,a)') &
        'Invalid ISTOMA_SCHEME, it must be between 0 and 1. Yours is set to',istoma_scheme,'...'
      call opspec_fatal(reason,'opspec_misc')  
      ifaterr = ifaterr +1
   end if

   if (iphen_scheme < 0 .or. iphen_scheme > 1) then
      write (reason,fmt='(a,1x,i4,a)') &
        'Invalid IPHEN_SCHEME, it must be between 0 and 1. Yours is set to',iphen_scheme,'...'
      call opspec_fatal(reason,'opspec_misc')  
      ifaterr = ifaterr +1
   end if

   if (n_plant_lim < 0 .or. n_plant_lim > 1) then
      write (reason,fmt='(a,1x,i4,a)') &
        'Invalid N_PLANT_LIM, it must be between 0 and 1. Yours is set to',n_plant_lim,'...'
      call opspec_fatal(reason,'opspec_misc')  
      ifaterr = ifaterr +1
   end if

   if (n_decomp_lim < 0 .or. n_decomp_lim > 1) then
      write (reason,fmt='(a,1x,i4,a)') &
        'Invalid N_DECOMP_LIM, it must be between 0 and 1. Yours is set to',n_decomp_lim,'...'
      call opspec_fatal(reason,'opspec_misc')  
      ifaterr = ifaterr +1
   end if

   if (include_fire < 0 .or. include_fire > 1) then
      write (reason,fmt='(a,1x,i4,a)') &
        'Invalid INCLUDE_FIRE, it must be between 0 and 1. Yours is set to',include_fire,'...'
      call opspec_fatal(reason,'opspec_misc')  
      ifaterr = ifaterr +1
   end if

   if (ianth_disturb < 0 .or. ianth_disturb > 1) then
      write (reason,fmt='(a,1x,i4,a)') &
        'Invalid IANTH_DISTURB, it must be between 0 and 1. Yours is set to',ianth_disturb,'...'
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
   
   if (pft_1st_check < 0 .or. pft_1st_check > 2) then
      write (reason,fmt='(a,1x,i4,a)') &
        'Invalid PFT_1ST_CHECK, it must be between 0 and 2. Yours is set to',pft_1st_check,'...'
      call opspec_fatal(reason,'opspec_misc')  
   end if
   
   !if (maxpatch < 0) then
   !   write (reason,fmt='(a,1x,i4,a)') &
   !     'Invalid MAXPATCH, it must be either 0 (no limit) or positive. Yours is set to',maxpatch,'...'
   !   call opspec_fatal(reason,'opspec_misc')  
   !   ifaterr = ifaterr +1
   !end if
   ! 
   !if (maxcohort < 0) then
   !   write (reason,fmt='(a,1x,i4,a)') &
   !     'Invalid MAXCOHORT, it must be either 0 (no limit) or positive. Yours is set to',maxcohort,'...'
   !   call opspec_fatal(reason,'opspec_misc')  
   !   ifaterr = ifaterr +1
   !end if
    
   if (zrough <= 0) then
      write (reason,fmt='(a,1x,es14.7,a)') &
        'Invalid ZROUGH, it must be positive. Yours is set to',zrough,'...'
      call opspec_fatal(reason,'opspec_misc')  
      ifaterr = ifaterr +1
   end if
    
   if (treefall_disturbance_rate < 0) then
      write (reason,fmt='(a,1x,es14.7,a)') &
        'Invalid TREEFALL_DISTURBANCE_RATE, it must be non-negative. Yours is set to' &
        ,treefall_disturbance_rate,'...'
      call opspec_fatal(reason,'opspec_misc')  
      ifaterr = ifaterr +1
   end if
    
   if (runoff_time < 0) then
      write (reason,fmt='(a,1x,es14.7,a)') &
        'Invalid RUNOFF_TIME, it must be non-negative. Yours is set to',runoff_time,'...'
      call opspec_fatal(reason,'opspec_misc')  
      ifaterr = ifaterr +1
   end if

   ! Stop the run if there are any fatal errors.
   if (ifaterr > 0) then
      write (unit=*,fmt='(a)')       ' -----------ED_OPSPEC_MISC --------------------------'
      write (unit=*,fmt='(a,1x,i5)') ' Fatal errors:',ifaterr
      write (unit=*,fmt='(a)')       ' ----------------------------------------------------'
      call fatal_error('Fatal errors at namelist - Misc settings '&
                      ,'ed_opspec_misc','ed_opspec.f90')
   end if
   return
end subroutine ed_opspec_misc
