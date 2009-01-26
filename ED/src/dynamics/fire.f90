subroutine fire_frequency_ar(month, cgrid)

  use ed_state_vars,only:edtype,polygontype,sitetype,patchtype
  use pft_coms, only: agf_bs, qsw, q
  use grid_coms, only: nzg
  use soil_coms, only: dslz
  use disturb_coms, only: fire_dryness_threshold, fire_parameter
  use allometry, only: ed_biomass

  implicit none

  type(edtype),target :: cgrid
  type(polygontype),pointer :: cpoly
  type(sitetype),pointer    :: csite
  type(patchtype),pointer   :: cpatch
  integer :: ipy,isi,ipa,ico
  integer, intent(in) :: month
  real :: ignition_rate
  real :: babove
  real :: patch_water_depth
  integer :: k
  real :: fuel

  ! Loop over polygons and sites.
        
        
  do ipy = 1,cgrid%npolygons
     
     cpoly => cgrid%polygon(ipy)
     
     do isi = 1,cpoly%nsites
        
        csite => cpoly%site(isi)
        
        ! Initialize
        ignition_rate = 0.0
        
        do ipa=1,csite%npatches
           
           cpatch => csite%patch(ipa)
           
           ! Initialize patch fuel.
           fuel = 0.0
           
           ! Add up fuel from all the cohorts
           
           do ico = 1,cpatch%ncohorts
              babove = ed_biomass(cpatch%bdead(ico), cpatch%balive(ico), cpatch%bleaf(ico),  &
                   cpatch%pft(ico), cpatch%hite(ico), cpatch%bstorage(ico)) * cpatch%nplant(ico)
              fuel = fuel + babove
           enddo
     
           ! Calculate patch water in meters
           patch_water_depth = sum(csite%sfcwater_depth(1:csite%nlev_sfcwater(ipa),ipa)) *  &
                0.001
           do k = cpoly%lsl(isi), nzg
              patch_water_depth = patch_water_depth + real(csite%soil_water(k,ipa) *   &
                   dble(dslz(k)))
           enddo
           
           ! calculate patch contribution to the ignition rate
           if(patch_water_depth < fire_dryness_threshold)then
              ignition_rate = ignition_rate + fuel * csite%area(ipa)
           endif
           
        enddo
        
        ! calculate fire dist rate [1/month]
        cpoly%lambda_fire(month,isi) = fire_parameter * ignition_rate
        cpoly%ignition_rate(isi) = ignition_rate
        cpoly%fuel(isi) = fuel
       

     enddo

     
  enddo

  return
end subroutine fire_frequency_ar
