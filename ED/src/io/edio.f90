module edio
   contains
   !=======================================================================================!
   !=======================================================================================!
   !     This is the main driver for file output in ED.                                    !
   !---------------------------------------------------------------------------------------!
   subroutine ed_output(analysis_time,new_day,new_year,dail_analy_time,mont_analy_time     &
                       ,dcyc_analy_time,annual_time,history_time,dcycle_time)
      use ed_state_vars, only : edgrid_g                & ! structure
                              , filltab_alltypes        & ! subroutine
                              , filltables              ! ! intent(inout)
      use grid_coms    , only : ngrids                  ! ! intent(in)
      use ed_misc_coms , only : isoutput                & ! intent(in)
                              , ifoutput                & ! intent(in)
                              , itoutput                & ! intent(in)
                              , writing_dail            & ! intent(in)
                              , writing_mont            & ! intent(in)
                              , writing_dcyc            & ! intent(in)
                              , iprintpolys             ! ! intent(in)
      use average_utils, only : aggregate_polygon_fmean & ! sub-routine
                              , normalize_ed_fmean_vars & ! sub-routine
                              , integrate_ed_dmean_vars & ! sub-routine
                              , integrate_ed_qmean_vars & ! sub-routine
                              , normalize_ed_dmean_vars & ! sub-routine
                              , integrate_ed_mmean_vars & ! sub-routine
                              , zero_ed_dmean_vars      & ! sub-routine
                              , normalize_ed_mmean_vars & ! sub-routine
                              , normalize_ed_qmean_vars & ! sub-routine
                              , zero_ed_mmean_vars      & ! sub-routine
                              , zero_ed_qmean_vars      & ! sub-routine
                              , zero_ed_qmean_vars      & ! sub-routine
                              , zero_ed_yearly_vars     ! ! sub-routine
      use ed_print     , only : print_fields            ! ! sub-routine
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      logical, intent(in)  :: analysis_time
      logical, intent(in)  :: dail_analy_time
      logical, intent(in)  :: mont_analy_time
      logical, intent(in)  :: dcyc_analy_time
      logical, intent(in)  :: history_time
      logical, intent(in)  :: dcycle_time
      logical, intent(in)  :: new_day
      logical, intent(in)  :: new_year
      logical, intent(in)  :: annual_time
      !----- Local variables. -------------------------------------------------------------!
      integer              :: ifm
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      If there is any IO, then we need to check whether the pointer tables are up   !
      ! to date or not.  They must be rehashed if there has been any change in the number  !
      ! of cohorts or patches (e.g, a cohort or a patch has been terminated or created).   !
      !------------------------------------------------------------------------------------!
      if (analysis_time   .or. history_time    .or.                                        &
          dail_analy_time .or. mont_analy_time .or. dcyc_analy_time .or. annual_time ) then

         if (filltables) then

            !----- Re-hash the tables. ----------------------------------------------------!
            call filltab_alltypes
            !----- Reset the rehash flag. -------------------------------------------------!
            filltables=.false.
         end if
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      If this is the time for an output, we shall call routines to prepare the      !
      ! variables for output.                                                              !
      !------------------------------------------------------------------------------------!
      if ( analysis_time .or.   history_time .or. dcycle_time  .or.                        &
          (new_day       .and. (writing_dail .or. writing_mont .or. writing_dcyc))) then
         do ifm=1,ngrids
            call normalize_ed_fmean_vars    (edgrid_g(ifm))
            call aggregate_polygon_fmean    (edgrid_g(ifm))
         end do

         if (writing_dail .or. writing_mont .or. writing_dcyc) then
            do ifm=1,ngrids
               call integrate_ed_dmean_vars(edgrid_g(ifm))
               if (writing_dcyc) call integrate_ed_qmean_vars(edgrid_g(ifm))
            end do
         end if
      end if
      !------------------------------------------------------------------------------------!



      !----- Instantaneous analysis. ------------------------------------------------------!
      if (analysis_time) then
         !----- Write out analysis fields - mostly polygon averages. ----------------------!
         if (ifoutput == 3) call h5_output('INST')
         if (itoutput == 3) call h5_output('OPTI')

         !----- If printpolys is on then print this info to the screen. -------------------!
         if (iprintpolys == 1) then
            do ifm=1,ngrids
               call print_fields(ifm,edgrid_g(ifm))
            end do
         end if
      end if
      !------------------------------------------------------------------------------------!



      !----- Daily analysis output and monthly integration. -------------------------------!
      if (new_day .and. (writing_dail .or. writing_mont .or. writing_dcyc)) then

         do ifm=1,ngrids
            call normalize_ed_dmean_vars(edgrid_g(ifm))
            if (writing_mont .or. writing_dcyc) then
               call integrate_ed_mmean_vars(edgrid_g(ifm))
            end if
         end do

         if (dail_analy_time) call h5_output('DAIL')

         do ifm=1,ngrids
            call zero_ed_dmean_vars(edgrid_g(ifm))
         end do
      end if
      !------------------------------------------------------------------------------------!




      !----- Monthly analysis and monthly mean diurnal cycle output. ----------------------!
      if (mont_analy_time .or. dcyc_analy_time) then
         do ifm=1,ngrids
            call normalize_ed_mmean_vars(edgrid_g(ifm))
            if (writing_dcyc) call normalize_ed_qmean_vars(edgrid_g(ifm))
         end do
         if (mont_analy_time) call h5_output('MONT')
         if (dcyc_analy_time) call h5_output('DCYC')
         do ifm=1,ngrids
            call zero_ed_mmean_vars(edgrid_g(ifm))
            if (writing_dcyc) call zero_ed_qmean_vars(edgrid_g(ifm))
         end do
      end if
      !------------------------------------------------------------------------------------!



      !----- Yearly analysis output. ------------------------------------------------------!
      if (annual_time) then
         call h5_output('YEAR')
      end if
      if (new_year) then
         do ifm=1,ngrids
            call zero_ed_yearly_vars(edgrid_g(ifm))
         end do
      end if
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !      History files should only be output at a frequency which divides by frqanl,   !
      ! thus the integrated fast-time variables are valid, but representative of the last  !
      ! frqanl period, not the last frqhist period.                                        !
      !------------------------------------------------------------------------------------!
      if (history_time .and. isoutput /= 0) then
         call h5_output('HIST')
      end if
      !------------------------------------------------------------------------------------!

      return
   end subroutine ed_output
   !=======================================================================================!
   !=======================================================================================!
end module edio
