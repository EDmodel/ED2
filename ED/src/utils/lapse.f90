!==========================================================================================!
!==========================================================================================!
!     This subroutine calculates the site-level variables based on the polygon-level ones. !
! If needed, it will also apply the lapse rate on some variables, and recalculate the      !
! others in order to keep satisfying the ideal gas law and some thermodynamic properties.  !
!------------------------------------------------------------------------------------------!
subroutine calc_met_lapse(cgrid,ipy)
  
   use ed_state_vars         , only : edtype      & ! structure
                                    , polygontype ! ! structure
   use consts_coms           , only : toodry,pi1      ! ! intent(in)
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
         hillshade = sin(pi1*cpoly%slope(isi)/180.)/2.
         cpoly%met(isi)%geoht       = cgrid%met(ipy)%geoht
         cpoly%met(isi)%atm_tmp     = cgrid%met(ipy)%atm_tmp
         cpoly%met(isi)%atm_shv     = cgrid%met(ipy)%atm_shv
         cpoly%met(isi)%prss        = cgrid%met(ipy)%prss
         cpoly%met(isi)%pcpg        = cgrid%met(ipy)%pcpg
         cpoly%met(isi)%par_diffuse = cgrid%met(ipy)%par_diffuse*(1.0-hillshade)
         cpoly%met(isi)%atm_co2     = cgrid%met(ipy)%atm_co2
         cpoly%met(isi)%rlong       = cgrid%met(ipy)%rlong
         cpoly%met(isi)%par_beam    = cgrid%met(ipy)%par_beam
         cpoly%met(isi)%nir_diffuse = cgrid%met(ipy)%nir_diffuse*(1.0-hillshade)
         cpoly%met(isi)%nir_beam    = cgrid%met(ipy)%nir_beam
         cpoly%met(isi)%vels        = cgrid%met(ipy)%vels
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
         cpoly%met(isi)%par_diffuse = (cgrid%met(ipy)%par_diffuse                           &
                                    + cgrid%lapse(ipy)%par_diffuse * delE)*(1.0-hillshade)
         cpoly%met(isi)%par_beam    = cgrid%met(ipy)%par_beam                              &
                                    + cgrid%lapse(ipy)%par_beam * delE
         cpoly%met(isi)%nir_diffuse = (cgrid%met(ipy)%nir_diffuse                           &
                                    + cgrid%lapse(ipy)%nir_diffuse * delE)*(1.0-hillshade)
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
         cpoly%met(isi)%pcpg    = cgrid%met(ipy)%pcpg * cpoly%pptweight(isi)

      end do
   end if

   !----- Check whether the radiation terms make sense... ---------------------------------!
   call met_sanity_check(cgrid,ipy)

   return
end subroutine calc_met_lapse
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!      Right now, this subroutine simply transfers lapse rates from ed_data.  In the       !
! future, it could set parmameters based on spatial maps of parameters.                    !
!------------------------------------------------------------------------------------------!
subroutine setLapseParms(cgrid)
   use ed_state_vars  , only : edtype      & ! structure
                             , polygontype ! ! structure
   use met_driver_coms, only : lapse       ! ! structure

   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   type(edtype)     , target  :: cgrid
   !----- Local variables. ----------------------------------------------------------------!
   type(polygontype), pointer :: cpoly
   real                       :: ebar
   real                       :: aterr
   integer                    :: ipy
   integer                    :: isi
   !---------------------------------------------------------------------------------------!

  
   do ipy = 1,cgrid%npolygons
      
      cgrid%lapse(ipy)%geoht       = lapse%geoht
      cgrid%lapse(ipy)%vels        = lapse%vels
      cgrid%lapse(ipy)%atm_tmp     = lapse%atm_tmp
      cgrid%lapse(ipy)%atm_shv     = lapse%atm_shv
      cgrid%lapse(ipy)%prss        = lapse%prss
      cgrid%lapse(ipy)%pcpg        = lapse%pcpg
      cgrid%lapse(ipy)%atm_co2     = lapse%atm_co2
      cgrid%lapse(ipy)%rlong       = lapse%rlong
      cgrid%lapse(ipy)%nir_beam    = lapse%nir_beam
      cgrid%lapse(ipy)%nir_diffuse = lapse%nir_diffuse
      cgrid%lapse(ipy)%par_beam    = lapse%par_beam
      cgrid%lapse(ipy)%par_diffuse = lapse%par_diffuse
      cgrid%lapse(ipy)%pptnorm     = lapse%pptnorm

      !----- Precipitation weight. --------------------------------------------------------!
      cpoly => cgrid%polygon(ipy)
      ebar  = 0.0
      aterr = 0.0
      do isi = 1,cpoly%nsites
         ebar  = ebar  + cpoly%area(isi)*cpoly%elevation(isi)
         aterr = aterr + cpoly%area(isi)
      end do

      ebar = ebar/aterr
      do isi = 1,cpoly%nsites
         if (cgrid%lapse(ipy)%pptnorm /= 0.) then
            cpoly%pptweight(isi) = ( cgrid%lapse(ipy)%pptnorm                              &
                                   + cgrid%lapse(ipy)%pcpg*(cpoly%elevation(isi) - ebar) ) &
                                 / cgrid%lapse(ipy)%pptnorm
         else
            cpoly%pptweight(isi) = 1.0
         end if
      end do

   end do

   return
end subroutine setLapseParms
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine will check the meteorological drivers, and print some information    !
! on the standard output, then print the banner that makes any user shiver...              !
!------------------------------------------------------------------------------------------!
subroutine met_sanity_check(cgrid,ipy)
   use ed_state_vars  , only : edtype       & ! structure
                             , polygontype  ! ! structure
   use met_driver_coms, only : rshort_min   & ! intent(in)
                             , rshort_max   & ! intent(in)
                             , rlong_min    & ! intent(in)
                             , rlong_max    & ! intent(in)
                             , atm_tmp_min  & ! intent(in)
                             , atm_tmp_max  & ! intent(in)
                             , atm_shv_min  & ! intent(in)
                             , atm_shv_max  & ! intent(in)
                             , atm_rhv_min  & ! intent(in)
                             , atm_rhv_max  & ! intent(in)
                             , atm_co2_min  & ! intent(in)
                             , atm_co2_max  & ! intent(in)
                             , prss_min     & ! intent(in)
                             , prss_max     & ! intent(in)
                             , pcpg_min     & ! intent(in)
                             , pcpg_max     & ! intent(in)
                             , vels_min     & ! intent(in)
                             , vels_max     & ! intent(in)
                             , geoht_min    & ! intent(in)
                             , geoht_max    ! ! intent(in)
   use ed_misc_coms   , only : current_time ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   type(edtype)     , target     :: cgrid
   integer          , intent(in) :: ipy
   !----- Local variables. ----------------------------------------------------------------!
   type(polygontype), pointer    :: cpoly
   integer                       :: isi
   integer                       :: ifaterr
   character(len=10)             :: now_date
   character(len=12)             :: now_time
   !----- Local constants. ----------------------------------------------------------------!
   logical          , parameter  :: fixRad = .true.
   character(len=3) , parameter  :: fmtc = '(a)'
   character(len=9) , parameter  :: fmti = '(a,1x,i6)'
   character(len=13), parameter  :: fmtf = '(a,1x,es12.5)'
   character(len=8) , parameter  :: fmtt = '(a,1x,a)'
   !---------------------------------------------------------------------------------------!


   !----- Reset the error counter. --------------------------------------------------------!
   ifaterr = 0
   cpoly => cgrid%polygon(ipy)

   !----- Write the current date and time. ------------------------------------------------!
   write(now_date,fmt='(2(i2.2,a),i4.4)') current_time%month,'/',current_time%date ,'/'    &
                                         ,current_time%year
   write(now_time,fmt='(3(i2.2,a))'     ) current_time%hour ,':',current_time%min  ,':'    &
                                         ,current_time%sec,' UTC'
   siteloop: do isi=1,cpoly%nsites
      if (cpoly%met(isi)%geoht < geoht_min .or.                                            &
          cpoly%met(isi)%geoht > geoht_max) then
         write (unit=*,fmt=fmtc) '---------------------------------------------------'
         write (unit=*,fmt=fmtc) ' Reference height doesn''t make sense... '
         write (unit=*,fmt=fmtc) ' '
         write (unit=*,fmt=fmti) ' - Polygon        :',ipy
         write (unit=*,fmt=fmti) ' - Site           :',isi
         write (unit=*,fmt=fmtt) ' - Date           :',now_date
         write (unit=*,fmt=fmtt) ' - Time           :',now_time
         write (unit=*,fmt=fmtf) ' - Longitude      :',cgrid%lon(ipy)
         write (unit=*,fmt=fmtf) ' - Latitude       :',cgrid%lat(ipy)
         write (unit=*,fmt=fmtf) ' - Site height    :',cpoly%met(isi)%geoht
         write (unit=*,fmt=fmtf) ' - Polygon height :',cgrid%met(ipy)%geoht
         write (unit=*,fmt=fmtf) ' - Minimum OK     :',geoht_min
         write (unit=*,fmt=fmtf) ' - Maximum OK     :',geoht_max
         write (unit=*,fmt=fmtc) '---------------------------------------------------'
         write (unit=*,fmt=fmtc) ' '
         ifaterr = ifaterr + 1
      end if

      if (cpoly%met(isi)%atm_tmp < atm_tmp_min .or.                                        &
          cpoly%met(isi)%atm_tmp > atm_tmp_max) then
         write (unit=*,fmt=fmtc) '---------------------------------------------------'
         write (unit=*,fmt=fmtc) ' Air temperature doesn''t make sense... '
         write (unit=*,fmt=fmtc) ' '
         write (unit=*,fmt=fmti) ' - Polygon        :',ipy
         write (unit=*,fmt=fmti) ' - Site           :',isi
         write (unit=*,fmt=fmtt) ' - Date           :',now_date
         write (unit=*,fmt=fmtt) ' - Time           :',now_time
         write (unit=*,fmt=fmtf) ' - Longitude      :',cgrid%lon(ipy)
         write (unit=*,fmt=fmtf) ' - Latitude       :',cgrid%lat(ipy)
         write (unit=*,fmt=fmtf) ' - Site temp.     :',cpoly%met(isi)%atm_tmp
         write (unit=*,fmt=fmtf) ' - Polygon temp.  :',cgrid%met(ipy)%atm_tmp
         write (unit=*,fmt=fmtf) ' - Minimum OK     :',atm_tmp_min
         write (unit=*,fmt=fmtf) ' - Maximum OK     :',atm_tmp_max
         write (unit=*,fmt=fmtc) '---------------------------------------------------'
         write (unit=*,fmt=fmtc) ' '
         ifaterr = ifaterr + 1
      end if

      if (cpoly%met(isi)%atm_shv < atm_shv_min .or.                                        &
          cpoly%met(isi)%atm_shv > atm_shv_max) then
         write (unit=*,fmt=fmtc) '---------------------------------------------------'
         write (unit=*,fmt=fmtc) ' Air specific humidity doesn''t make sense... '
         write (unit=*,fmt=fmtc) ' '
         write (unit=*,fmt=fmti) ' - Polygon        :',ipy
         write (unit=*,fmt=fmti) ' - Site           :',isi
         write (unit=*,fmt=fmtt) ' - Date           :',now_date
         write (unit=*,fmt=fmtt) ' - Time           :',now_time
         write (unit=*,fmt=fmtf) ' - Longitude      :',cgrid%lon(ipy)
         write (unit=*,fmt=fmtf) ' - Latitude       :',cgrid%lat(ipy)
         write (unit=*,fmt=fmtf) ' - Site spec. hum.:',cpoly%met(isi)%atm_shv
         write (unit=*,fmt=fmtf) ' - Polygon sp. h. :',cgrid%met(ipy)%atm_shv
         write (unit=*,fmt=fmtf) ' - Minimum OK     :',atm_shv_min
         write (unit=*,fmt=fmtf) ' - Maximum OK     :',atm_shv_max
         write (unit=*,fmt=fmtc) '---------------------------------------------------'
         write (unit=*,fmt=fmtc) ' '
         ifaterr = ifaterr + 1
      end if

      if (cpoly%met(isi)%atm_co2 < atm_co2_min .or.                                        &
          cpoly%met(isi)%atm_co2 > atm_co2_max) then
         write (unit=*,fmt=fmtc) '---------------------------------------------------'
         write (unit=*,fmt=fmtc) ' Air carbon dioxide doesn''t make sense... '
         write (unit=*,fmt=fmtc) ' '
         write (unit=*,fmt=fmti) ' - Polygon        :',ipy
         write (unit=*,fmt=fmti) ' - Site           :',isi
         write (unit=*,fmt=fmtt) ' - Date           :',now_date
         write (unit=*,fmt=fmtt) ' - Time           :',now_time
         write (unit=*,fmt=fmtf) ' - Longitude      :',cgrid%lon(ipy)
         write (unit=*,fmt=fmtf) ' - Latitude       :',cgrid%lat(ipy)
         write (unit=*,fmt=fmtf) ' - Site CO2       :',cpoly%met(isi)%atm_co2
         write (unit=*,fmt=fmtf) ' - Polygon CO2    :',cgrid%met(ipy)%atm_co2
         write (unit=*,fmt=fmtf) ' - Minimum OK     :',atm_co2_min
         write (unit=*,fmt=fmtf) ' - Maximum OK     :',atm_co2_max
         write (unit=*,fmt=fmtc) '---------------------------------------------------'
         write (unit=*,fmt=fmtc) ' '
         ifaterr = ifaterr + 1
      end if

      if (cpoly%met(isi)%prss < prss_min .or.                                              &
          cpoly%met(isi)%prss > prss_max) then
         write (unit=*,fmt=fmtc) '---------------------------------------------------'
         write (unit=*,fmt=fmtc) ' Atmospheric pressure doesn''t make sense... '
         write (unit=*,fmt=fmtc) ' '
         write (unit=*,fmt=fmti) ' - Polygon        :',ipy
         write (unit=*,fmt=fmti) ' - Site           :',isi
         write (unit=*,fmt=fmtt) ' - Date           :',now_date
         write (unit=*,fmt=fmtt) ' - Time           :',now_time
         write (unit=*,fmt=fmtf) ' - Longitude      :',cgrid%lon(ipy)
         write (unit=*,fmt=fmtf) ' - Latitude       :',cgrid%lat(ipy)
         write (unit=*,fmt=fmtf) ' - Site press.    :',cpoly%met(isi)%prss
         write (unit=*,fmt=fmtf) ' - Polygon press. :',cgrid%met(ipy)%prss
         write (unit=*,fmt=fmtf) ' - Minimum OK     :',prss_min
         write (unit=*,fmt=fmtf) ' - Maximum OK     :',prss_max
         write (unit=*,fmt=fmtc) '---------------------------------------------------'
         write (unit=*,fmt=fmtc) ' '
         ifaterr = ifaterr + 1
      end if

      if (cpoly%met(isi)%pcpg < pcpg_min .or.                                              &
          cpoly%met(isi)%pcpg > pcpg_max) then
         write (unit=*,fmt=fmtc) '---------------------------------------------------'
         write (unit=*,fmt=fmtc) ' Precipitation rate doesn''t make sense... '
         write (unit=*,fmt=fmtc) ' '
         write (unit=*,fmt=fmti) ' - Polygon        :',ipy
         write (unit=*,fmt=fmti) ' - Site           :',isi
         write (unit=*,fmt=fmtt) ' - Date           :',now_date
         write (unit=*,fmt=fmtt) ' - Time           :',now_time
         write (unit=*,fmt=fmtf) ' - Longitude      :',cgrid%lon(ipy)
         write (unit=*,fmt=fmtf) ' - Latitude       :',cgrid%lat(ipy)
         write (unit=*,fmt=fmtf) ' - Site Precip    :',cpoly%met(isi)%pcpg
         write (unit=*,fmt=fmtf) ' - Polygon Precip :',cgrid%met(ipy)%pcpg
         write (unit=*,fmt=fmtf) ' - Minimum OK     :',pcpg_min
         write (unit=*,fmt=fmtf) ' - Maximum OK     :',pcpg_max
         write (unit=*,fmt=fmtc) '---------------------------------------------------'
         write (unit=*,fmt=fmtc) ' '
         ifaterr = ifaterr + 1
      end if

      !------------------------------------------------------------------------------------!
      !     At this point, vels is twice the kinetic energy.  Check the square root of     !
      ! vels, so we are comparing apples to apples.                                        !
      !------------------------------------------------------------------------------------!
      if (cpoly%met(isi)%vels < vels_min * vels_min) then
         write (unit=*,fmt=fmtc) '---------------------------------------------------'
         write (unit=*,fmt=fmtc) ' Kinetic energy is negative...'
         write (unit=*,fmt=fmtc) ' '
         write (unit=*,fmt=fmti) ' - Polygon                :',ipy
         write (unit=*,fmt=fmti) ' - Site                   :',isi
         write (unit=*,fmt=fmtt) ' - Date                   :',now_date
         write (unit=*,fmt=fmtt) ' - Time                   :',now_time
         write (unit=*,fmt=fmtf) ' - Longitude              :',cgrid%lon(ipy)
         write (unit=*,fmt=fmtf) ' - Latitude               :',cgrid%lat(ipy)
         write (unit=*,fmt=fmtf) ' - Site kinetic energy    :',0.5 * cpoly%met(isi)%vels
         write (unit=*,fmt=fmtf) ' - Polygon kinetic energy :',0.5 * cgrid%met(ipy)%vels
         write (unit=*,fmt=fmtf) ' - Minimum OK             :',0.5 * vels_min * vels_min
         write (unit=*,fmt=fmtc) '---------------------------------------------------'
         write (unit=*,fmt=fmtc) ' '
         ifaterr = ifaterr + 1
      elseif (cpoly%met(isi)%vels > vels_max * vels_max) then
         write (unit=*,fmt=fmtc) '---------------------------------------------------'
         write (unit=*,fmt=fmtc) ' Wind speed doesn''t make sense... '
         write (unit=*,fmt=fmtc) ' '
         write (unit=*,fmt=fmti) ' - Polygon        :',ipy
         write (unit=*,fmt=fmti) ' - Site           :',isi
         write (unit=*,fmt=fmtt) ' - Date           :',now_date
         write (unit=*,fmt=fmtt) ' - Time           :',now_time
         write (unit=*,fmt=fmtf) ' - Longitude      :',cgrid%lon(ipy)
         write (unit=*,fmt=fmtf) ' - Latitude       :',cgrid%lat(ipy)
         write (unit=*,fmt=fmtf) ' - Site wind      :',sqrt(cpoly%met(isi)%vels)
         write (unit=*,fmt=fmtf) ' - Polygon wind   :',sqrt(max(0.,cgrid%met(ipy)%vels))
         write (unit=*,fmt=fmtf) ' - Minimum OK     :',vels_min
         write (unit=*,fmt=fmtf) ' - Maximum OK     :',vels_max
         write (unit=*,fmt=fmtc) '---------------------------------------------------'
         write (unit=*,fmt=fmtc) ' '
         ifaterr = ifaterr + 1
      end if
      !------------------------------------------------------------------------------------!

      if (cpoly%met(isi)%rlong < rlong_min .or.                                            &
          cpoly%met(isi)%rlong > rlong_max) then
         write (unit=*,fmt=fmtc) '---------------------------------------------------'
         write (unit=*,fmt=fmtc) ' Longwave radiation doesn''t make sense... '
         write (unit=*,fmt=fmtc) ' '
         write (unit=*,fmt=fmti) ' - Polygon          :',ipy
         write (unit=*,fmt=fmti) ' - Site             :',isi
         write (unit=*,fmt=fmtt) ' - Date             :',now_date
         write (unit=*,fmt=fmtt) ' - Time             :',now_time
         write (unit=*,fmt=fmtf) ' - Longitude        :',cgrid%lon(ipy)
         write (unit=*,fmt=fmtf) ' - Latitude         :',cgrid%lat(ipy)
         write (unit=*,fmt=fmtf) ' - Site Longwave    :',cpoly%met(isi)%rlong
         write (unit=*,fmt=fmtf) ' - Polygon Longwave :',cgrid%met(ipy)%rlong
         write (unit=*,fmt=fmtf) ' - Site Temp        :',cpoly%met(isi)%atm_tmp
         write (unit=*,fmt=fmtf) ' - Site SH          :',cpoly%met(isi)%atm_shv
         write (unit=*,fmt=fmtf) ' - Minimum OK       :',rlong_min
         write (unit=*,fmt=fmtf) ' - Maximum OK       :',rlong_max
         write (unit=*,fmt=fmtc) '---------------------------------------------------'
         write (unit=*,fmt=fmtc) ' '
         ifaterr = ifaterr + 1
      end if

      if ((cpoly%met(isi)%par_diffuse < rshort_min .and. (.not. fixRad)) .or.              &
           cpoly%met(isi)%par_diffuse > rshort_max) then
         write (unit=*,fmt=fmtc) '---------------------------------------------------'
         write (unit=*,fmt=fmtc) ' Diffuse PAR radiation doesn''t make sense... '
         write (unit=*,fmt=fmtc) ' '
         write (unit=*,fmt=fmti) ' - Polygon        :',ipy
         write (unit=*,fmt=fmti) ' - Site           :',isi
         write (unit=*,fmt=fmtt) ' - Date           :',now_date
         write (unit=*,fmt=fmtt) ' - Time           :',now_time
         write (unit=*,fmt=fmtf) ' - Longitude      :',cgrid%lon(ipy)
         write (unit=*,fmt=fmtf) ' - Latitude       :',cgrid%lat(ipy)
         write (unit=*,fmt=fmtf) ' - Site rshort    :',cpoly%met(isi)%par_diffuse
         write (unit=*,fmt=fmtf) ' - Polygon rshort :',cgrid%met(ipy)%par_diffuse
         write (unit=*,fmt=fmtf) ' - Minimum OK     :',rshort_min
         write (unit=*,fmt=fmtf) ' - Maximum OK     :',rshort_max
         write (unit=*,fmt=fmtc) '---------------------------------------------------'
         write (unit=*,fmt=fmtc) ' '
         ifaterr = ifaterr + 1
      elseif (cpoly%met(isi)%par_diffuse < rshort_min) then
         cpoly%met(isi)%par_diffuse = rshort_min
      end if

      if ((cpoly%met(isi)%par_beam < rshort_min .and. (.not. fixRad)) .or.              &
           cpoly%met(isi)%par_beam > rshort_max) then
         write (unit=*,fmt=fmtc) '---------------------------------------------------'
         write (unit=*,fmt=fmtc) ' Direct PAR radiation doesn''t make sense... '
         write (unit=*,fmt=fmtc) ' '
         write (unit=*,fmt=fmti) ' - Polygon        :',ipy
         write (unit=*,fmt=fmti) ' - Site           :',isi
         write (unit=*,fmt=fmtt) ' - Date           :',now_date
         write (unit=*,fmt=fmtt) ' - Time           :',now_time
         write (unit=*,fmt=fmtf) ' - Longitude      :',cgrid%lon(ipy)
         write (unit=*,fmt=fmtf) ' - Latitude       :',cgrid%lat(ipy)
         write (unit=*,fmt=fmtf) ' - Site PAR       :',cpoly%met(isi)%par_beam
         write (unit=*,fmt=fmtf) ' - Polygon PAR    :',cgrid%met(ipy)%par_beam
         write (unit=*,fmt=fmtf) ' - Minimum OK     :',rshort_min
         write (unit=*,fmt=fmtf) ' - Maximum OK     :',rshort_max
         write (unit=*,fmt=fmtc) '---------------------------------------------------'
         write (unit=*,fmt=fmtc) ' '
         ifaterr = ifaterr + 1
      elseif (cpoly%met(isi)%par_beam < rshort_min) then
         cpoly%met(isi)%par_beam = rshort_min
      end if

      if ((cpoly%met(isi)%nir_diffuse < rshort_min .and. (.not. fixRad)) .or.              &
           cpoly%met(isi)%nir_diffuse > rshort_max) then
         write (unit=*,fmt=fmtc) '---------------------------------------------------'
         write (unit=*,fmt=fmtc) ' Diffuse near IR radiation doesn''t make sense... '
         write (unit=*,fmt=fmtc) ' '
         write (unit=*,fmt=fmti) ' - Polygon        :',ipy
         write (unit=*,fmt=fmti) ' - Site           :',isi
         write (unit=*,fmt=fmtt) ' - Date           :',now_date
         write (unit=*,fmt=fmtt) ' - Time           :',now_time
         write (unit=*,fmt=fmtf) ' - Longitude      :',cgrid%lon(ipy)
         write (unit=*,fmt=fmtf) ' - Latitude       :',cgrid%lat(ipy)
         write (unit=*,fmt=fmtf) ' - Site radNIR    :',cpoly%met(isi)%nir_diffuse
         write (unit=*,fmt=fmtf) ' - Polygon radNIR :',cgrid%met(ipy)%nir_diffuse
         write (unit=*,fmt=fmtf) ' - Minimum OK     :',rshort_min
         write (unit=*,fmt=fmtf) ' - Maximum OK     :',rshort_max
         write (unit=*,fmt=fmtc) '---------------------------------------------------'
         write (unit=*,fmt=fmtc) ' '
         ifaterr = ifaterr + 1
      elseif (cpoly%met(isi)%nir_diffuse < rshort_min) then
         cpoly%met(isi)%nir_diffuse = rshort_min
      end if

      if ((cpoly%met(isi)%nir_beam < rshort_min .and. (.not. fixRad)) .or.                 &
           cpoly%met(isi)%nir_beam > rshort_max) then
         write (unit=*,fmt=fmtc) '---------------------------------------------------'
         write (unit=*,fmt=fmtc) ' Direct near IR radiation doesn''t make sense... '
         write (unit=*,fmt=fmtc) ' '
         write (unit=*,fmt=fmti) ' - Polygon        :',ipy
         write (unit=*,fmt=fmti) ' - Site           :',isi
         write (unit=*,fmt=fmtt) ' - Date           :',now_date
         write (unit=*,fmt=fmtt) ' - Time           :',now_time
         write (unit=*,fmt=fmtf) ' - Longitude      :',cgrid%lon(ipy)
         write (unit=*,fmt=fmtf) ' - Latitude       :',cgrid%lat(ipy)
         write (unit=*,fmt=fmtf) ' - Site radNIR    :',cpoly%met(isi)%nir_beam
         write (unit=*,fmt=fmtf) ' - Polygon radNIR :',cgrid%met(ipy)%nir_beam
         write (unit=*,fmt=fmtf) ' - Minimum OK     :',rshort_min
         write (unit=*,fmt=fmtf) ' - Maximum OK     :',rshort_max
         write (unit=*,fmt=fmtc) '---------------------------------------------------'
         write (unit=*,fmt=fmtc) ' '
         ifaterr = ifaterr + 1
      elseif (cpoly%met(isi)%nir_beam < rshort_min) then
         cpoly%met(isi)%nir_beam = rshort_min
      end if
      
      !------------------------------------------------------------------------------------!
      !     Also, check that the sum of all four components do not exceed the maximum      !
      ! shortwave radiation. Otherwise, recalculate the radiation to account for some      !
      ! of the components that may have been affected.                                     !
      !------------------------------------------------------------------------------------!
      if (fixRad) then
         cpoly%met(isi)%rshort_diffuse = cpoly%met(isi)%par_diffuse                        &
                                       + cpoly%met(isi)%nir_diffuse
         cpoly%met(isi)%rshort         = cpoly%met(isi)%rshort_diffuse                     &
                                       + cpoly%met(isi)%par_beam                           &
                                       + cpoly%met(isi)%nir_beam
      end if
      !----- Then check whether the full terms do not exceed maximum... -------------------!
      if (cpoly%met(isi)%rshort_diffuse > rshort_max) then
         write (unit=*,fmt=fmtc) '---------------------------------------------------'
         write (unit=*,fmt=fmtc) ' Diffuse SW radiation doesn''t make sense... '
         write (unit=*,fmt=fmtc) ' '
         write (unit=*,fmt=fmti) ' - Polygon        :',ipy
         write (unit=*,fmt=fmti) ' - Site           :',isi
         write (unit=*,fmt=fmtt) ' - Date           :',now_date
         write (unit=*,fmt=fmtt) ' - Time           :',now_time
         write (unit=*,fmt=fmtf) ' - Longitude      :',cgrid%lon(ipy)
         write (unit=*,fmt=fmtf) ' - Latitude       :',cgrid%lat(ipy)
         write (unit=*,fmt=fmtf) ' - PAR diffuse    :',cpoly%met(isi)%par_diffuse
         write (unit=*,fmt=fmtf) ' - NIR diffuse    :',cpoly%met(isi)%nir_diffuse
         write (unit=*,fmt=fmtf) ' - Site rshortd   :',cpoly%met(isi)%rshort_diffuse
         write (unit=*,fmt=fmtf) ' - Polygon rshortd:',cgrid%met(ipy)%rshort_diffuse
         write (unit=*,fmt=fmtf) ' - Minimum OK     :',rshort_min
         write (unit=*,fmt=fmtf) ' - Maximum OK     :',rshort_max
         write (unit=*,fmt=fmtc) '---------------------------------------------------'
         write (unit=*,fmt=fmtc) ' '
         ifaterr = ifaterr + 1
      end if

      if (cpoly%met(isi)%rshort > rshort_max) then
         write (unit=*,fmt=fmtc) '---------------------------------------------------'
         write (unit=*,fmt=fmtc) ' Total SW radiation doesn''t make sense... '
         write (unit=*,fmt=fmtc) ' '
         write (unit=*,fmt=fmti) ' - Polygon        :',ipy
         write (unit=*,fmt=fmti) ' - Site           :',isi
         write (unit=*,fmt=fmtt) ' - Date           :',now_date
         write (unit=*,fmt=fmtt) ' - Time           :',now_time
         write (unit=*,fmt=fmtf) ' - Longitude      :',cgrid%lon(ipy)
         write (unit=*,fmt=fmtf) ' - Latitude       :',cgrid%lat(ipy)
         write (unit=*,fmt=fmtf) ' - PAR direct     :',cpoly%met(isi)%par_beam
         write (unit=*,fmt=fmtf) ' - PAR diffuse    :',cpoly%met(isi)%par_diffuse
         write (unit=*,fmt=fmtf) ' - NIR direct     :',cpoly%met(isi)%nir_beam
         write (unit=*,fmt=fmtf) ' - NIR diffuse    :',cpoly%met(isi)%nir_diffuse
         write (unit=*,fmt=fmtf) ' - Site rshort    :',cpoly%met(isi)%rshort
         write (unit=*,fmt=fmtf) ' - Polygon rshort :',cgrid%met(ipy)%rshort
         write (unit=*,fmt=fmtf) ' - Minimum OK     :',rshort_min
         write (unit=*,fmt=fmtf) ' - Maximum OK     :',rshort_max
         write (unit=*,fmt=fmtc) '---------------------------------------------------'
         write (unit=*,fmt=fmtc) ' '
         ifaterr = ifaterr + 1
      end if
      
      if (ifaterr > 0) then
         call fatal_error('Meteorological forcing has issues (check message).'             &
                         ,'met_sanity_check','lapse.f90')
      end if

   end do siteloop 

   return
end subroutine met_sanity_check
!==========================================================================================!
!==========================================================================================!

