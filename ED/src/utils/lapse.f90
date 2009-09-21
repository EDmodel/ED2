!==========================================================================================!
!==========================================================================================!
!     This subroutine calculates the site-level variables based on the polygon-level ones. !
! If needed, it will also apply the lapse rate on some variables, and recalculate the      !
! others in order to keep satisfying the ideal gas law and some thermodynamic properties.  !
!------------------------------------------------------------------------------------------!
subroutine calc_met_lapse(cgrid,ipy)
  
   use ed_state_vars         , only : edtype      & ! structure
                                    , polygontype ! ! structure
   use canopy_radiation_coms , only : rlong_min   ! ! intent(in)
   use consts_coms           , only : toodry      ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   type(edtype)     , target     :: cgrid
   integer          , intent(in) :: ipy
   !----- Local variables -----------------------------------------------------------------!
   type(polygontype), pointer    :: cpoly
   integer                       :: isi
   real                          :: ebar   !! mean elevation
   real                          :: delE   !! deviation from mean elevation
   real                          :: aterr  !! terrestrial area
   !----- Local constants -----------------------------------------------------------------!
   real             , parameter  :: offset=tiny(1.)/epsilon(1.) !! Tiny offset to avoid FPE
   logical          , parameter  :: bypass=.true.
   !---------------------------------------------------------------------------------------!

   !----- Pass over sites once to calc preliminary stats. ---------------------------------!
   cpoly => cgrid%polygon(ipy)

   ebar = 0.0
   aterr = 0.0

   do isi=1,cpoly%nsites
      ebar  = ebar + cpoly%area(isi)*cpoly%elevation(isi)
      aterr = aterr + cpoly%area(isi)
   end do
   ebar = ebar/aterr

   if (bypass) then
      do isi = 1,cpoly%nsites
         cpoly%met(isi)%geoht       = cgrid%met(ipy)%geoht
         cpoly%met(isi)%atm_tmp     = cgrid%met(ipy)%atm_tmp
         cpoly%met(isi)%atm_shv     = cgrid%met(ipy)%atm_shv
         cpoly%met(isi)%prss        = cgrid%met(ipy)%prss
         cpoly%met(isi)%pcpg        = cgrid%met(ipy)%pcpg
         cpoly%met(isi)%par_diffuse = cgrid%met(ipy)%par_diffuse
         cpoly%met(isi)%atm_co2     = cgrid%met(ipy)%atm_co2
         cpoly%met(isi)%rlong       = cgrid%met(ipy)%rlong
         cpoly%met(isi)%par_beam    = cgrid%met(ipy)%par_beam
         cpoly%met(isi)%nir_diffuse = cgrid%met(ipy)%nir_diffuse
         cpoly%met(isi)%nir_beam    = cgrid%met(ipy)%nir_beam
         cpoly%met(isi)%vels        = cgrid%met(ipy)%vels
       
         !---------------------------------------------------------------------------------!
         !     Sanity check:                                                               !
         !---------------------------------------------------------------------------------!
         if ( cpoly%met(isi)%rlong < rlong_min) then
            write(unit=*,fmt='(a)')           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            write(unit=*,fmt='(a)')           '+ RLONG is too low!!!'
            write(unit=*,fmt='(a,1x,es12.5)') '+ Longitude       : ',cgrid%lon(ipy)
            write(unit=*,fmt='(a,1x,es12.5)') '+ Latitude        : ',cgrid%lat(ipy)
            write(unit=*,fmt='(a,1x,es12.5)') '+ Site-level      : ',cpoly%met(isi)%rlong
            write(unit=*,fmt='(a,1x,es12.5)') '+ Polygon-level   : ',cgrid%met(ipy)%rlong
            write(unit=*,fmt='(a,1x,es12.5)') '+ Minimum OK value: ',rlong_min
            write(unit=*,fmt='(a)')           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            call fatal_error('Problems with RLONG A','calc_met_lapse'                      &
                            ,'edcp_met.f90')
         elseif ( cpoly%met(isi)%atm_tmp < 150.0) then
            write(unit=*,fmt='(a)')           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            write(unit=*,fmt='(a)')           '+ ATM_TMP is too low!!!'
            write(unit=*,fmt='(a,1x,es12.5)') '+ Longitude       : ',cgrid%lon(ipy)
            write(unit=*,fmt='(a,1x,es12.5)') '+ Latitude        : ',cgrid%lat(ipy)
            write(unit=*,fmt='(a,1x,es12.5)') '+ Site-level      : ',cpoly%met(isi)%atm_tmp
            write(unit=*,fmt='(a,1x,es12.5)') '+ Polygon-level   : ',cgrid%met(ipy)%atm_tmp
            write(unit=*,fmt='(a,1x,es12.5)') '+ Minimum OK value: ',150.0
            write(unit=*,fmt='(a)')           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            call fatal_error('Problems with ATM_TMP','calc_met_lapse'                      &
                            ,'edcp_met.f90')
         else if ( cpoly%met(isi)%atm_shv < toodry) then
            write(unit=*,fmt='(a)')           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            write(unit=*,fmt='(a)')           '+ ATM_SHV is too low!!!'
            write(unit=*,fmt='(a,1x,es12.5)') '+ Longitude       : ',cgrid%lon(ipy)
            write(unit=*,fmt='(a,1x,es12.5)') '+ Latitude        : ',cgrid%lat(ipy)
            write(unit=*,fmt='(a,1x,es12.5)') '+ Site-level      : ',cpoly%met(isi)%atm_shv
            write(unit=*,fmt='(a,1x,es12.5)') '+ Polygon-level   : ',cgrid%met(ipy)%atm_shv
            write(unit=*,fmt='(a,1x,es12.5)') '+ Minimum OK value: ',1.e-5
            write(unit=*,fmt='(a)')           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            call fatal_error('Problems with ATM_SHV','calc_met_lapse'                      &
                            ,'edcp_met.f90')
         elseif ( cpoly%met(isi)%rlong > 600.0) then
            write(unit=*,fmt='(a)')           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            write(unit=*,fmt='(a)')           '+ RLONG is too high!!!'
            write(unit=*,fmt='(a,1x,es12.5)') '+ Longitude       : ',cgrid%lon(ipy)
            write(unit=*,fmt='(a,1x,es12.5)') '+ Latitude        : ',cgrid%lat(ipy)
            write(unit=*,fmt='(a,1x,es12.5)') '+ Site-level      : ',cpoly%met(isi)%rlong
            write(unit=*,fmt='(a,1x,es12.5)') '+ Polygon-level   : ',cgrid%met(ipy)%rlong
            write(unit=*,fmt='(a,1x,es12.5)') '+ Maximum OK value: ',600.0
            write(unit=*,fmt='(a)')           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            call fatal_error('Problems with RLONG A','calc_met_lapse'                      &
                            ,'edcp_met.f90')
         elseif ( cpoly%met(isi)%atm_tmp > 317.0) then
            write(unit=*,fmt='(a)')           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            write(unit=*,fmt='(a)')           '+ ATM_TMP is too high!!!'
            write(unit=*,fmt='(a,1x,es12.5)') '+ Longitude       : ',cgrid%lon(ipy)
            write(unit=*,fmt='(a,1x,es12.5)') '+ Latitude        : ',cgrid%lat(ipy)
            write(unit=*,fmt='(a,1x,es12.5)') '+ Site-level      : ',cpoly%met(isi)%atm_tmp
            write(unit=*,fmt='(a,1x,es12.5)') '+ Polygon-level   : ',cgrid%met(ipy)%atm_tmp
            write(unit=*,fmt='(a,1x,es12.5)') '+ Maximum OK value: ',317.0
            write(unit=*,fmt='(a)')           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            call fatal_error('Problems with ATM_TMP','calc_met_lapse'                      &
                            ,'edcp_met.f90')
         elseif ( cpoly%met(isi)%atm_shv > 3.0e-2) then
            write(unit=*,fmt='(a)')           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            write(unit=*,fmt='(a)')           '+ ATM_SHV is too high!!!'
            write(unit=*,fmt='(a,1x,es12.5)') '+ Longitude       : ',cgrid%lon(ipy)
            write(unit=*,fmt='(a,1x,es12.5)') '+ Latitude        : ',cgrid%lat(ipy)
            write(unit=*,fmt='(a,1x,es12.5)') '+ Site-level      : ',cpoly%met(isi)%atm_shv
            write(unit=*,fmt='(a,1x,es12.5)') '+ Polygon-level   : ',cgrid%met(ipy)%atm_shv
            write(unit=*,fmt='(a,1x,es12.5)') '+ Maximum OK value: ',30.0e-3
            write(unit=*,fmt='(a)')           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            call fatal_error('Problems with ATM_SHV','calc_met_lapse'                      &
                            ,'edcp_met.f90')
         elseif ( cpoly%met(isi)%par_beam + cpoly%met(isi)%nir_beam &
                + cpoly%met(isi)%par_diffuse + cpoly%met(isi)%nir_diffuse > 1320.0 ) then
            write(unit=*,fmt='(a)')           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            write(unit=*,fmt='(a)')           '+ SOLAR RADIATION is non-sense !!!'
            write(unit=*,fmt='(a,1x,es12.5)') '+ Longitude       : ',cgrid%lon(ipy)
            write(unit=*,fmt='(a,1x,es12.5)') '+ Latitude        : ',cgrid%lat(ipy)
            write(unit=*,fmt='(a,1x,es12.5)') '+ PAR_BEAM        : ',cpoly%met(isi)%par_beam
            write(unit=*,fmt='(a,1x,es12.5)') '+ PAR_DIFFUSE     : '                       &
                                              , cpoly%met(isi)%par_diffuse
            write(unit=*,fmt='(a,1x,es12.5)') '+ NIR_BEAM        : ',cpoly%met(isi)%nir_beam
            write(unit=*,fmt='(a,1x,es12.5)') '+ NIR_DIFFUSE     : '                       &
                                              , cpoly%met(isi)%nir_diffuse
            write(unit=*,fmt='(a,1x,es12.5)') '+ Maximum OK value (sum): ',1320.0
            write(unit=*,fmt='(a)')           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            call fatal_error('Problems with solar radiation','calc_met_lapse'              &
                            ,'ed_met_driver.f90')
         end if
      end do
    
   else
      
      !----- Second pass, calculate lapse rate adjustment. --------------------------------!
      do isi = 1,cpoly%nsites
         
         delE = cpoly%elevation(isi) - ebar
         
         !----- Perform linear adjustments. -----------------------------------------------!
         cpoly%met(isi)%geoht   = cgrid%met(ipy)%geoht   + cgrid%lapse(ipy)%geoht   * delE
         cpoly%met(isi)%atm_tmp = cgrid%met(ipy)%atm_tmp + cgrid%lapse(ipy)%atm_tmp * delE
         cpoly%met(isi)%atm_shv = cgrid%met(ipy)%atm_shv + cgrid%lapse(ipy)%atm_shv * delE
         cpoly%met(isi)%prss    = cgrid%met(ipy)%prss    + cgrid%lapse(ipy)%prss    * delE
         cpoly%met(isi)%pcpg    = cgrid%met(ipy)%pcpg    + cgrid%lapse(ipy)%pcpg    * delE
         cpoly%met(isi)%atm_co2 = cgrid%met(ipy)%atm_co2 + cgrid%lapse(ipy)%atm_co2 * delE
         cpoly%met(isi)%rlong   = cgrid%met(ipy)%rlong   + cgrid%lapse(ipy)%rlong   * delE
         cpoly%met(isi)%par_diffuse = cgrid%met(ipy)%par_diffuse                           &
                                    + cgrid%lapse(ipy)%par_diffuse * delE
         cpoly%met(isi)%par_beam    = cgrid%met(ipy)%par_beam                              &
                                    + cgrid%lapse(ipy)%par_beam * delE
         cpoly%met(isi)%nir_diffuse = cgrid%met(ipy)%nir_diffuse                           &
                                    + cgrid%lapse(ipy)%nir_diffuse * delE
         cpoly%met(isi)%nir_beam    = cgrid%met(ipy)%nir_beam                              &
                                    + cgrid%lapse(ipy)%nir_beam * delE
         !---------------------------------------------------------------------------------!
         ! Note: at this point VELS is vel^2.  Thus this lapse preserves mean wind ENERGY  !
         !       not wind SPEED.                                                           !
         !---------------------------------------------------------------------------------!
         cpoly%met(isi)%vels    = cgrid%met(ipy)%vels    + cgrid%lapse(ipy)%vels*delE

         !---------------------------------------------------------------------------------!
         !     Sanity check:                                                               !
         !---------------------------------------------------------------------------------!
         if ( cpoly%met(isi)%rlong < rlong_min) then
            write(unit=*,fmt='(a)')           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            write(unit=*,fmt='(a)')           '+ RLONG is too low!!!'
            write(unit=*,fmt='(a,1x,es12.5)') '+ Site-level      : ',cpoly%met(isi)%rlong
            write(unit=*,fmt='(a,1x,es12.5)') '+ Polygon-level   : ',cgrid%met(ipy)%rlong
            write(unit=*,fmt='(a,1x,es12.5)') '+ Minimum OK value: ',rlong_min
            write(unit=*,fmt='(a)')           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            call fatal_error('Problems with RLONG A','calc_met_lapse'                      &
                            ,'edcp_met.f90')
         elseif ( cpoly%met(isi)%atm_tmp < 150.0) then
            write(unit=*,fmt='(a)')           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            write(unit=*,fmt='(a)')           '+ ATM_TMP is too low!!!'
            write(unit=*,fmt='(a,1x,es12.5)') '+ Site-level      : ',cpoly%met(isi)%atm_tmp
            write(unit=*,fmt='(a,1x,es12.5)') '+ Polygon-level   : ',cgrid%met(ipy)%atm_tmp
            write(unit=*,fmt='(a,1x,es12.5)') '+ Minimum OK value: ',150.0
            write(unit=*,fmt='(a)')           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            call fatal_error('Problems with ATM_TMP','calc_met_lapse'                      &
                            ,'edcp_met.f90')
         else if ( cpoly%met(isi)%atm_shv < 1.e-5) then
            write(unit=*,fmt='(a)')           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            write(unit=*,fmt='(a)')           '+ ATM_SHV is too low!!!'
            write(unit=*,fmt='(a,1x,es12.5)') '+ Site-level      : ',cpoly%met(isi)%atm_shv
            write(unit=*,fmt='(a,1x,es12.5)') '+ Polygon-level   : ',cgrid%met(ipy)%atm_shv
            write(unit=*,fmt='(a,1x,es12.5)') '+ Minimum OK value: ',1.e-5
            write(unit=*,fmt='(a)')           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            call fatal_error('Problems with ATM_SHV','calc_met_lapse'                      &
                            ,'edcp_met.f90')
         elseif ( cpoly%met(isi)%rlong > 600.0) then
            write(unit=*,fmt='(a)')           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            write(unit=*,fmt='(a)')           '+ RLONG is too high!!!'
            write(unit=*,fmt='(a,1x,es12.5)') '+ Site-level      : ',cpoly%met(isi)%rlong
            write(unit=*,fmt='(a,1x,es12.5)') '+ Polygon-level   : ',cgrid%met(ipy)%rlong
            write(unit=*,fmt='(a,1x,es12.5)') '+ Maximum OK value: ',600.0
            write(unit=*,fmt='(a)')           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            call fatal_error('Problems with RLONG A','calc_met_lapse'                      &
                            ,'edcp_met.f90')
         elseif ( cpoly%met(isi)%atm_tmp > 317.0) then
            write(unit=*,fmt='(a)')           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            write(unit=*,fmt='(a)')           '+ ATM_TMP is too high!!!'
            write(unit=*,fmt='(a,1x,es12.5)') '+ Site-level      : ',cpoly%met(isi)%atm_tmp
            write(unit=*,fmt='(a,1x,es12.5)') '+ Polygon-level   : ',cgrid%met(ipy)%atm_tmp
            write(unit=*,fmt='(a,1x,es12.5)') '+ Maximum OK value: ',317.0
            write(unit=*,fmt='(a)')           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            call fatal_error('Problems with ATM_TMP','calc_met_lapse'                      &
                            ,'edcp_met.f90')
         elseif ( cpoly%met(isi)%atm_shv > 30.0e-3) then
            write(unit=*,fmt='(a)')           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            write(unit=*,fmt='(a)')           '+ ATM_SHV is too high!!!'
            write(unit=*,fmt='(a,1x,es12.5)') '+ Site-level      : ',cpoly%met(isi)%atm_shv
            write(unit=*,fmt='(a,1x,es12.5)') '+ Polygon-level   : ',cgrid%met(ipy)%atm_shv
            write(unit=*,fmt='(a,1x,es12.5)') '+ Maximum OK value: ',30.0e-3
            write(unit=*,fmt='(a)')           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            call fatal_error('Problems with ATM_SHV','calc_met_lapse'                      &
                            ,'edcp_met.f90')
         elseif ( cpoly%met(isi)%par_beam + cpoly%met(isi)%nir_beam &
                + cpoly%met(isi)%par_diffuse + cpoly%met(isi)%nir_diffuse > 1320.0 ) then
            write(unit=*,fmt='(a)')           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            write(unit=*,fmt='(a)')           '+ SOLAR RADIATION is non-sense !!!'
            write(unit=*,fmt='(a,1x,es12.5)') '+ PAR_BEAM        : ',cpoly%met(isi)%par_beam
            write(unit=*,fmt='(a,1x,es12.5)') '+ PAR_DIFFUSE     : '                       &
                                              , cpoly%met(isi)%par_diffuse
            write(unit=*,fmt='(a,1x,es12.5)') '+ NIR_BEAM        : ',cpoly%met(isi)%nir_beam
            write(unit=*,fmt='(a,1x,es12.5)') '+ NIR_DIFFUSE     : '                       &
                                              , cpoly%met(isi)%nir_diffuse
            write(unit=*,fmt='(a,1x,es12.5)') '+ Maximum OK value (sum): ',1320.0
            write(unit=*,fmt='(a)')           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            call fatal_error('Problems with solar radiation','calc_met_lapse'              &
                            ,'ed_met_driver.f90')
         end if
      end do
   end if
   return
end subroutine calc_met_lapse
!==========================================================================================!
!==========================================================================================!
