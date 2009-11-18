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
   real                          :: hillshade
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
         hillshade = cpoly%slope(isi)/180.
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
            write(unit=*,fmt='(a,1x,i12)')    '+ X point         : ',cgrid%ilon(ipy)
            write(unit=*,fmt='(a,1x,i12)')    '+ Y point         : ',cgrid%ilat(ipy)
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
         cpoly%met(isi)%atm_shv = max(toodry, cgrid%met(ipy)%atm_shv                       &
                                            + cgrid%lapse(ipy)%atm_shv * delE)
         cpoly%met(isi)%prss    = cgrid%met(ipy)%prss    + cgrid%lapse(ipy)%prss    * delE
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
         ! Note: Precipitation adjustment is based on proportional change rather than a    !
         !       linear scaling.                                                           !
         !---------------------------------------------------------------------------------!
         cpoly%met(isi)%pcpg    = cgrid%met(ipy)%pcpg*cpoly%pptweight(isi)

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

         call MetDiagnostics(cpoly,ipy,isi)

      end do
   end if
   return
end subroutine calc_met_lapse
!==========================================================================================!
!==========================================================================================!



subroutine MetDiagnostics(cpoly,ipy,isi)
  use ed_state_vars,only    : edtype,polygontype
  use canopy_radiation_coms,only : rlong_min
  
  implicit none
  integer, intent(in) :: ipy
  integer, intent(in) :: isi
  type(polygontype),target :: cpoly

  logical, parameter :: fixRad = .true.


  if(cpoly%met(isi)%geoht .le. 0.)  then
     print*,cpoly%met(isi)%geoht
     call fatal_error('Problems with GEOHT','MetDiagnostics','ed_met_driver.f90')
  endif
  if(cpoly%met(isi)%atm_tmp .le. 200. .or. cpoly%met(isi)%atm_tmp .ge. 350.)  then
     print*,cpoly%met(isi)%atm_tmp
     call fatal_error('Problems with atm_tmp','MetDiagnostics','ed_met_driver.f90')
  endif
  if(cpoly%met(isi)%atm_shv .le. 0. .or. cpoly%met(isi)%atm_shv .ge. 1.)  then
     print*,cpoly%met(isi)%atm_shv
     call fatal_error('Problems with atm_shv','MetDiagnostics','ed_met_driver.f90')
  endif
  if(cpoly%met(isi)%prss .le. 0.)  then
     print*,cpoly%met(isi)%prss
     call fatal_error('Problems with prss','MetDiagnostics','ed_met_driver.f90')
  endif
  if(cpoly%met(isi)%pcpg .lt. 0. .or. cpoly%met(isi)%pcpg .gt. 0.1)  then
     print*,cpoly%met(isi)%pcpg
     call fatal_error('Problems with precipitation','MetDiagnostics','ed_met_driver.f90')
  endif
  if(cpoly%met(isi)%atm_co2 .le. 100. .or. cpoly%met(isi)%atm_co2 .gt. 2000.)  then
     print*,cpoly%met(isi)%atm_co2
     call fatal_error('Problems with ATM CO2','MetDiagnostics','ed_met_driver.f90')
  endif
  if(cpoly%met(isi)%rlong .le. rlong_min)  then
     print*,cpoly%met(isi)%rlong
     call fatal_error('Problems with LONGWAVE RADIATION','MetDiagnostics','ed_met_driver.f90')
  endif
  if(cpoly%met(isi)%par_diffuse .lt. 0. .or. cpoly%met(isi)%par_diffuse .gt. 1400.)  then
     if(cpoly%met(isi)%par_diffuse .lt. 0. .and. fixRad) then
        cpoly%met(isi)%par_diffuse = 0.0
     else
        print*,cpoly%met(isi)%par_diffuse
        call fatal_error('Problems with DIFFUSE PAR','MetDiagnostics','ed_met_driver.f90')
     endif
  endif
  if(cpoly%met(isi)%par_beam .lt. 0. .or. cpoly%met(isi)%par_beam .gt. 1400.)  then
     if(cpoly%met(isi)%par_beam .lt. 0. .and. fixRad) then
        cpoly%met(isi)%par_beam = 0.0
     else
        print*,cpoly%met(isi)%par_beam
        call fatal_error('Problems with DIRECT BEAM PAR','MetDiagnostics','ed_met_driver.f90')
     end if
  endif
  if(cpoly%met(isi)%nir_diffuse .lt. 0. .or. cpoly%met(isi)%nir_diffuse .gt. 1400.)  then
     if(cpoly%met(isi)%nir_diffuse .lt. 0. .and. fixRad) then
        cpoly%met(isi)%nir_diffuse = 0.0
     else
        print*,cpoly%met(isi)%nir_diffuse
        call fatal_error('Problems with DIFUSE NIR','MetDiagnostics','ed_met_driver.f90')
     endif
  end if

  if(cpoly%met(isi)%nir_beam .lt. 0. .or. cpoly%met(isi)%nir_beam .gt. 1400.)  then
     if(cpoly%met(isi)%nir_beam .lt. 0. .and. fixRad) then
        cpoly%met(isi)%nir_beam = 0.0
     else
        print*,cpoly%met(isi)%nir_beam
        call fatal_error('Problems with DIRECT BEAM NIR','MetDiagnostics','ed_met_driver.f90')
     end if
  endif
  if(cpoly%met(isi)%vels .lt. 0.)  then
     print*,cpoly%met(isi)%vels
     call fatal_error('Problems with WIND VELOCITY','MetDiagnostics','ed_met_driver.f90')
  endif

  return
end subroutine MetDiagnostics


!==========================================================================================!
!==========================================================================================!
subroutine setLapseParms(cgrid)
  
  use ed_state_vars,only:edtype,polygontype
  use met_driver_coms, only: lapse

   implicit none

  integer :: ipy,isi
  type(edtype), target :: cgrid
  type(polygontype),pointer :: cpoly
  real :: ebar, aterr

   !! right now, simply transfer lapse rates from ed_data
   !! in future, could set parms based on spatial maps of parms
  
  do ipy = 1,cgrid%npolygons
     
     cgrid%lapse(ipy)%geoht   = lapse%geoht
     cgrid%lapse(ipy)%vels    = lapse%vels
     cgrid%lapse(ipy)%atm_tmp = lapse%atm_tmp
     cgrid%lapse(ipy)%atm_shv = lapse%atm_shv
     cgrid%lapse(ipy)%prss    = lapse%prss
     cgrid%lapse(ipy)%pcpg    = lapse%pcpg
     cgrid%lapse(ipy)%atm_co2 = lapse%atm_co2
     cgrid%lapse(ipy)%rlong   = lapse%rlong
     cgrid%lapse(ipy)%nir_beam    = lapse%nir_beam
     cgrid%lapse(ipy)%nir_diffuse = lapse%nir_diffuse
     cgrid%lapse(ipy)%par_beam    = lapse%par_beam
     cgrid%lapse(ipy)%par_diffuse = lapse%par_diffuse
     cgrid%lapse(ipy)%pptnorm = lapse%pptnorm

     !!! PRECIP WEIGHTS !!!
     cpoly => cgrid%polygon(ipy)
     ebar = 0.0
     aterr = 0.0
     do isi = 1,cpoly%nsites
        ebar = ebar + cpoly%area(isi)*cpoly%elevation(isi)
        aterr = aterr + cpoly%area(isi)
     enddo
     ebar = ebar/aterr
     do isi = 1,cpoly%nsites
        cpoly%pptweight(isi) = (cgrid%lapse(ipy)%pptnorm + &
             cgrid%lapse(ipy)%pcpg*(cpoly%elevation(isi) - ebar))&
             /cgrid%lapse(ipy)%pptnorm
     enddo

  enddo

  return
end subroutine setLapseParms
!==========================================================================================!
!==========================================================================================!
