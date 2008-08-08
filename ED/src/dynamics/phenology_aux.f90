subroutine prescribed_leaf_state(lat, imonth, iyear, doy,   &
     green_leaf_factor, leaf_aging_factor, phen_pars)
  
  ! Calculate phenology factors for prescribed phenology schemes.

  use phenology_coms, only: iphenys1, iphenysf, iphenyf1, iphenyff,  &
       prescribed_phen
  use max_dims, only: n_pft
  use pft_coms, only: phenology
  
  implicit none
  
  type(prescribed_phen), intent(in) :: phen_pars
  real, intent(in) :: lat
  integer, intent(in) :: imonth
  integer :: n_recycle_years
  integer :: my_year
  integer, intent(in) :: iyear
  integer, intent(in) :: doy
  real :: elongf
  real :: delay
  integer :: pft
  real, dimension(n_pft), intent(out) :: green_leaf_factor
  real, dimension(n_pft), intent(out) :: leaf_aging_factor
  
  if( (lat >= 0.0 .and. imonth <= 7) .or.   & ! in northern hemisphere, 
       ! this assumes dropping between August 1 and December 31 and flushing
       ! between January 1 and July 31.
       (lat < 0.0 .and. (imonth > 7 .or. imonth == 1)) )then ! in the 
     ! southern hemisphere, this assumes dropping between February 1 and
     ! July 31 and flushing between August 1 and January 31.
     
     ! get the year
     n_recycle_years = iphenysf - iphenys1 + 1
     if(iyear > iphenysf)then
        my_year = mod(iyear-iphenys1,n_recycle_years) + 1
     elseif(iyear < iphenys1)then
        my_year = n_recycle_years - mod(iphenysf-iyear,n_recycle_years)
     else
        my_year = iyear - iphenys1 + 1
     endif
     
     ! calculate the factors
     elongf = 1.0 / (1.0 +   &
          (phen_pars%flush_a(my_year) * doy)**phen_pars%flush_b(my_year))
     delay = elongf

!     print*,phen_pars%flush_a(my_year),phen_pars%flush_b(my_year),elongf
     
  else
     ! leaves turning color
     
     ! get the year
     n_recycle_years = iphenyff - iphenyf1 + 1
     if(iyear > iphenyff)then
        my_year = mod(iyear-iphenyf1,n_recycle_years) + 1
     elseif(iyear < iphenyf1)then
        my_year = n_recycle_years - mod(iphenyff-iyear,n_recycle_years)
     else
        my_year = iyear - iphenyf1 + 1
     endif
     
     ! calculate the factors
     elongf = 1.0 / (1.0 +   &
          (phen_pars%color_a(my_year) * doy)**phen_pars%color_b(my_year))
     delay = 1.0 / (1.0 +   &
          (phen_pars%color_a(my_year) *   &
          doy * 1.095)**phen_pars%color_b(my_year))
 !    print*,phen_pars%color_a(my_year),phen_pars%color_b(my_year),elongf,doy
  endif
  
  ! load the values for each PFT
  do pft = 1, n_pft
     if(phenology(pft) == 2)then
        green_leaf_factor(pft) = elongf
        leaf_aging_factor(pft) = delay
     elseif(phenology(pft) == 0 .or. phenology(pft) == 1)then
        green_leaf_factor(pft) = 1.0
        leaf_aging_factor(pft) = 1.0
     endif
  enddo
  
  return
end subroutine prescribed_leaf_state
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine update_thermal_sums_ar(month, cpoly, isi, lat)
  
  use ed_state_vars,only:polygontype,sitetype

  implicit none

  type(polygontype),target :: cpoly
  type(sitetype),pointer   :: csite
  integer :: isi,ipa
  integer, intent(in) :: month
  real, intent(in) :: lat


  ! Chill days --- number of days with average temperatures below 278.15 K
  ! Degree days --- sum of daily average temperatures above 278.15 K 

  ! loop over patches

  csite => cpoly%site(isi)
  
  do ipa = 1,csite%npatches

     ! Minimum monthly temperature of the site
     cpoly%min_monthly_temp(isi) = min(cpoly%min_monthly_temp(isi), csite%avg_daily_temp(ipa))
     
     if(csite%avg_daily_temp(ipa) > 278.15)then  
        
        ! update dgd
        
        if(lat >= 0.0)then  
           
           ! northern hemisphere
           
           if(month <= 8)then 
              
              ! update only for Jan-Aug.
              
              csite%sum_dgd(ipa) = csite%sum_dgd(ipa)   &
                   + (csite%avg_daily_temp(ipa)-278.15)
              
           else 
              
              ! not accumulating degree days in Sep-Dec
              
              csite%sum_dgd(ipa) = 0.0
              
           endif
           
        else 
           
           ! in southern hemisphere
           
           if(month <= 2 .or. month >= 7)then 
              
              ! Accumulating only Jul-Feb
              
              csite%sum_dgd(ipa) = csite%sum_dgd(ipa) +   &
                   csite%avg_daily_temp(ipa) - 278.15
              
           else 
              
              ! not accumulating degree days Mar-Jun
              
              csite%sum_dgd(ipa) = 0.0
              
           endif
        endif
     else 
        
        ! update chilling days
        
        if(lat >= 0.0)then  
           
           ! northern hemisphere
           
           if(month >= 11 .or. month <= 6)then 
              
              ! accumulate only Nov-Jun
              
              csite%sum_chd(ipa) = csite%sum_chd(ipa) + 1.0
              
           else 
              
              ! not accumulating chilling days, Jul-Oct
              
              csite%sum_chd(ipa) = 0.0
              
           endif
        else 
           
           ! in southern hemisphere
           
           if(month >= 5)then 
              
              ! accumulate only in May-Dec
              
              csite%sum_chd(ipa) = csite%sum_chd(ipa) + 1.0
              
           else 
              
              ! not accumulating chilling days Jan-Apr
              
              csite%sum_chd(ipa) = 0.0
              
           endif
        endif
     endif
     
  enddo
  
  return
end subroutine update_thermal_sums_ar
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
real function daylength(lat,day)

  use consts_coms, only: pio180

  implicit none
  
  real :: lat
  real :: arg
  integer :: day
  
  arg = -tan(lat*pio180)*tan(-23.5*pio180*cos(6.283/365.0*(float(day)+9.0)))
  if( abs(arg) < 1.0 )then
     daylength = 120.0 * acos(arg)/(15.0*pio180)
  elseif(arg >= 1.0)then
     daylength = 0.0
  elseif(arg <= -1.0)then
     daylength = 1440.0
  endif

  return
end function daylength
