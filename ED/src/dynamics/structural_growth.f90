subroutine structural_growth_ar(cgrid, month)

  ! Do not change the order of operations below unless you know what you
  ! are doing.  Changing the order can affect the C/N budgets.

  use ed_state_vars,only: edtype,polygontype,sitetype,patchtype
  use pft_coms, only: q, qsw, seedling_mortality, c2n_leaf, c2n_storage,   &
       c2n_recruit, c2n_stem, l2n_stem
  use decomp_coms, only: f_labile
  use max_dims, only: n_pft, n_dbh
  use therm_lib,only: update_veg_energy_cweh

  implicit none

  type(edtype),target       :: cgrid
  type(polygontype),pointer :: cpoly
  type(sitetype),pointer    :: csite
  type(patchtype),pointer   :: cpatch
  integer :: ipy,isi,ipa,ico
  integer, intent(in) :: month

  real :: salloc
  real :: salloci
  real :: balive_in
  real :: bdead_in
  real :: hite_in
  real :: dbh_in
  real :: nplant_in
  real :: bstorage_in
  real :: f_bseeds
  real :: f_bdead
  real :: balive_mort_litter
  real :: bstorage_mort_litter
  real :: struct_litter
  real :: mort_litter
  real :: seed_litter
  real :: net_seed_N_uptake
  real :: net_stem_N_uptake
  integer :: update_month
  real :: cb_act
  real :: cb_max
  integer :: imonth

 

  do ipy = 1,cgrid%npolygons
     
     cpoly => cgrid%polygon(ipy)
             
     ! Initialization
     cpoly%basal_area(:,:,:) = 0.0
     cpoly%agb(:,:,:) = 0.0
     ! cs%basal_area_growth  = 0.0
     ! cpoly%agb_growth  = 0.0
     ! cpoly%basal_area_mort = 0.0
     ! cpoly%agb_mort = 0.0

     do isi = 1,cpoly%nsites
        
        csite => cpoly%site(isi)

        do ipa=1,csite%npatches

           cpatch => csite%patch(ipa)

           do ico = 1,cpatch%ncohorts

              salloc = 1.0 + q(cpatch%pft(ico)) + qsw(cpatch%pft(ico)) * cpatch%hite(ico)
              salloci = 1.0 / salloc

              ! Remember inputs in order to calculate increments later on.
              balive_in = cpatch%balive(ico)
              bdead_in = cpatch%bdead(ico)
              hite_in = cpatch%hite(ico)
              dbh_in = cpatch%dbh(ico)
              nplant_in = cpatch%nplant(ico)
              bstorage_in = cpatch%bstorage(ico)

              ! Apply mortality.  Do not allow nplant < 3.0e-8.  Such a 
              ! sparse cohort will soon be terminated, anyway.  
              ! NB: monthly_dndt is negative.
              cpatch%monthly_dndt(ico) = max(cpatch%monthly_dndt(ico), 3.0e-8 - cpatch%nplant(ico))

              cpatch%nplant(ico) = cpatch%nplant(ico) + cpatch%monthly_dndt(ico)

              ! Calculate litter owing to mortality.
              balive_mort_litter = - cpatch%balive(ico) * cpatch%monthly_dndt(ico)
              bstorage_mort_litter = - cpatch%bstorage(ico) * cpatch%monthly_dndt(ico)
              struct_litter = - cpatch%bdead(ico) * cpatch%monthly_dndt(ico)
              mort_litter = balive_mort_litter + bstorage_mort_litter +  &
                   struct_litter

              ! Reset monthly_dndt.
              cpatch%monthly_dndt(ico) = 0.0

              ! Determine how to distribute what is in bstorage.
              call plant_structural_allocation_ar(cpatch,ico, month, f_bseeds,   &
                   f_bdead, cgrid%lat(ipy))

              ! Grow plants; bdead gets fraction f_bdead of bstorage.
              cpatch%bdead(ico) = cpatch%bdead(ico) + f_bdead * cpatch%bstorage(ico)
              
              ! Rebalance the plant nitrogen uptake considering the actual 
              ! allocation to structural growth.  This is necessary because
              ! c2n_stem does not necessarily equal c2n_storage.
              net_stem_N_uptake = (cpatch%bdead(ico) - bdead_in) * ( 1.0 / c2n_stem &
                   - 1.0 / c2n_storage) * cpatch%nplant(ico)
              
              ! Calculate total seed production and seed litter.  The seed
              ! pool gets a fraction f_bseeds of bstorage.

              cpatch%bseeds(ico) = f_bseeds * cpatch%bstorage(ico)
              seed_litter = cpatch%bseeds(ico) * cpatch%nplant(ico) * seedling_mortality(cpatch%pft(ico))
              
              ! Rebalance the plant nitrogen uptake considering the actual 
              ! allocation to seeds.  This is necessary because c2n_recruit
              ! does not have to equal c2n_storage.
              net_seed_N_uptake = cpatch%bseeds(ico) * (1.0 / c2n_recruit(cpatch%pft(ico)) &
                   - 1.0 / c2n_storage) * cpatch%nplant(ico)
              
              ! Decrement the storage pool.
              cpatch%bstorage(ico) = cpatch%bstorage(ico) * (1.0 - f_bdead - f_bseeds)

              ! Finalize litter inputs
              csite%fsc_in(ipa) = csite%fsc_in(ipa) + f_labile(cpatch%pft(ico)) *   &
                   balive_mort_litter + bstorage_mort_litter + seed_litter
              csite%fsn_in(ipa) = csite%fsn_in(ipa) + f_labile(cpatch%pft(ico)) *   &
                   balive_mort_litter / c2n_leaf(cpatch%pft(ico)) +   &
                   bstorage_mort_litter/ c2n_storage + seed_litter /  &
                   c2n_recruit(cpatch%pft(ico))
              csite%ssc_in(ipa) = csite%ssc_in(ipa) + (1.0 - f_labile(cpatch%pft(ico))) *  &
                   balive_mort_litter + struct_litter
              csite%ssl_in(ipa) = csite%ssl_in(ipa) + ((1.0 - f_labile(cpatch%pft(ico))) * &
                   balive_mort_litter + struct_litter ) * l2n_stem / c2n_stem
              csite%total_plant_nitrogen_uptake(ipa) = &
                   csite%total_plant_nitrogen_uptake(ipa)   &
                   + net_seed_N_uptake + net_stem_N_uptake

              ! Calculate the derived cohort properties
              call update_derived_cohort_props_ar(cpatch,ico,   &
                   cpoly%green_leaf_factor(cpatch%pft(ico),isi), cpoly%lsl(isi))
              

              ! Update the vegetation internal energy and heat capacity
              call update_veg_energy_cweh(cpatch,ico)


              ! Update annual average carbon balances for mortality
              update_month = month - 1
              if(update_month == 0)update_month = 12
              cpatch%cb(update_month,ico) = cpatch%cb(13,ico)
              cpatch%cb_max(update_month,ico) = cpatch%cb_max(13,ico)
              cpatch%cb(13,ico) = 0.0
              cpatch%cb_max(13,ico) = 0.0
              cb_act = 0.0
              cb_max = 0.0
              do imonth = 1,12
                 cb_act = cb_act + cpatch%cb(imonth,ico)
                 cb_max = cb_max + cpatch%cb_max(imonth,ico)
              enddo
              if(cb_max > 0.0)then
                 cpatch%cbr_bar(ico) = cb_act / cb_max
              else
                 cpatch%cbr_bar(ico) = 0.0
              endif

              ! Update interesting output quantities
              call update_vital_rates_ar(cpatch,ico, dbh_in, bdead_in,   &
                   balive_in, hite_in, bstorage_in, nplant_in, mort_litter,  &
                   csite%area(ipa), cpoly%basal_area(:,:,isi), cpoly%agb(:,:,isi),  &
                   cpoly%basal_area_growth(:,:,isi), cpoly%agb_growth(:,:,isi),   &
                   cpoly%basal_area_mort(:,:,isi), cpoly%agb_mort(:,:,isi))

           enddo
           
           ! Age the patch if this is not agriculture
           if(csite%dist_type(ipa) /= 1)csite%age(ipa) = csite%age(ipa) + 1.0/12.0
           

        enddo
        
     enddo
     
  enddo


  return
end subroutine structural_growth_ar

!===============================================================

subroutine plant_structural_allocation_ar(cpatch,ico, month, f_bseeds, f_bdead, lat)
  
  use ed_state_vars,only:patchtype
  use pft_coms, only: phenology, repro_min_h, r_fract

  implicit none

  type(patchtype),target :: cpatch
  integer :: ico

  integer, intent(in) :: month
  real, intent(out) :: f_bseeds
  real, intent(out) :: f_bdead
  real, intent(in) :: lat

  ! Calculate fraction of bstorage going to bdead and reproduction
  if(phenology(cpatch%pft(ico)) /= 2   .or.  &  ! for NOT broad leaf deciduous
       (lat >= 0.0 .and. month == 6) .or.  &  ! or Jun in north
       (lat < 0.0 .and. month == 12) )then    ! or Dec in south
     
     ! For all PFTs except broadleaf deciduous
     
     if(cpatch%hite(ico) <= repro_min_h(cpatch%pft(ico)))then
        f_bseeds = 0.0
     else

        f_bseeds = r_fract(cpatch%pft(ico))
     endif
     f_bdead = 1.0 - f_bseeds
     
  else
     
     f_bdead = 0.0
     f_bseeds = 0.0
     
  endif
        
  return
end subroutine plant_structural_allocation_ar

!===================================================================

subroutine update_derived_cohort_props_ar(cpatch,ico, green_leaf_factor, lsl)

  use ed_state_vars,only:patchtype
  use pft_coms, only: phenology, sla, q, qsw

  implicit none

  type(patchtype),target :: cpatch
  integer :: ico
  real, intent(in) :: green_leaf_factor
  real, external :: bd2dbh
  real, external :: dbh2h
  real :: bl
  real :: bl_max
  real, external :: dbh2bl
  real :: rootdepth
  real :: calc_root_depth
  integer, external :: assign_root_depth
  integer, intent(in) :: lsl

  cpatch%dbh(ico) = bd2dbh(cpatch%pft(ico), cpatch%bdead(ico)) 
  cpatch%hite(ico) = dbh2h(cpatch%pft(ico), cpatch%dbh(ico))
     
  if(cpatch%phenology_status(ico) /= 2)then

     ! Update status
     bl = cpatch%balive(ico) * green_leaf_factor / (1.0 +   &
          q(cpatch%pft(ico)) + qsw(cpatch%pft(ico)) * cpatch%hite(ico))
     bl_max = dbh2bl(cpatch%dbh(ico),cpatch%pft(ico)) * green_leaf_factor
     if(bl.lt.bl_max)then
        cpatch%phenology_status(ico) = 1
     else
        cpatch%phenology_status(ico) = 0
     endif
     
     ! Update LAI
     cpatch%lai(ico) = bl * cpatch%nplant(ico) * sla(cpatch%pft(ico))
     cpatch%bleaf(ico) = bl

  endif

  ! Update rooting depth
  rootdepth = calc_root_depth(cpatch%hite(ico), cpatch%dbh(ico), cpatch%pft(ico))
  
  ! See which discrete soil level this corresponds to
  cpatch%krdepth(ico) = assign_root_depth(rootdepth, lsl)
  
  return
end subroutine update_derived_cohort_props_ar

!=====================================================================

subroutine update_vital_rates_ar(cpatch,ico, dbh_in, bdead_in, balive_in, hite_in,  &
     bstorage_in, nplant_in, mort_litter, area, basal_area, agb,  &
     basal_area_growth, agb_growth, basal_area_mort, agb_mort)
   
  use ed_state_vars,only:patchtype
  use max_dims, only: n_pft, n_dbh
  use consts_coms, only: pi1
  use pft_coms, only: agf_bs,q,qsw
  
  implicit none

  real, intent(in) :: dbh_in
  real, intent(in) :: bdead_in
  real, intent(in) :: balive_in
  real, intent(in) :: hite_in
  real, intent(in) :: bstorage_in
  real, intent(in) :: nplant_in
  real, intent(in) :: mort_litter  
  
  type(patchtype),target :: cpatch
  integer :: ico

  integer :: bdbh
  real, external :: ed_biomass
  real, intent(in) :: area
  real, dimension(n_pft, n_dbh) :: basal_area
  real, dimension(n_pft, n_dbh) :: agb
  real, dimension(n_pft, n_dbh) :: basal_area_growth
  real, dimension(n_pft, n_dbh) :: agb_growth
  real, dimension(n_pft, n_dbh) :: basal_area_mort
  real, dimension(n_pft, n_dbh) :: agb_mort

  ! Get dbh bin
  bdbh = min(int(dbh_in*0.1),10)+1

  ! Update current basal area, agb
  basal_area(cpatch%pft(ico), bdbh) = basal_area(cpatch%pft(ico), bdbh) + area * cpatch%nplant(ico) *  &
       pi1 * 0.25 * cpatch%dbh(ico)**2

  agb(cpatch%pft(ico), bdbh) = agb(cpatch%pft(ico), bdbh) + area * cpatch%nplant(ico) * 10.0 * &
       ed_biomass(cpatch%bdead(ico), cpatch%balive(ico), cpatch%bleaf(ico), cpatch%pft(ico), &
       cpatch%hite(ico), cpatch%bstorage(ico)) 

  ! Only update rates for cohorts on the first census.
  if(cpatch%first_census(ico) /= 1)return

  ! Computed for plants alive both at past census and current census
  basal_area_growth(cpatch%pft(ico),bdbh) = basal_area_growth(cpatch%pft(ico),bdbh) +  &
       area * cpatch%nplant(ico) * pi1 * 0.25 * (cpatch%dbh(ico)**2 - dbh_in**2)
  agb_growth(cpatch%pft(ico),bdbh) = agb_growth(cpatch%pft(ico),bdbh) +  &
       area * cpatch%nplant(ico) * 10.0 * (ed_biomass(cpatch%bdead(ico), cpatch%balive(ico),   &
       cpatch%bleaf(ico), cpatch%pft(ico), cpatch%hite(ico), cpatch%bstorage(ico)) - ed_biomass(bdead_in,   &
       balive_in, cpatch%bleaf(ico), cpatch%pft(ico), hite_in, bstorage_in))
  ! note bleaf is unchanged.

  ! Computed for plants alive at past census but dead at current census
  basal_area_mort(cpatch%pft(ico),bdbh) = basal_area_mort(cpatch%pft(ico),bdbh) +  &
       area * (nplant_in - cpatch%nplant(ico)) * pi1 * 0.25 * dbh_in**2
        
  agb_mort(cpatch%pft(ico),bdbh) = agb_mort(cpatch%pft(ico),bdbh) + area * mort_litter * 10.0

  return
end subroutine update_vital_rates_ar
        

!==========================================================================================!
!==========================================================================================!
subroutine print_C_and_N_budgets(cgrid)
  
  use ed_state_vars,only:edtype,polygontype,sitetype,patchtype
  
  implicit none

  type(edtype),target :: cgrid
  type(polygontype),pointer :: cpoly
  type(sitetype),pointer    :: csite
  type(patchtype),pointer   :: cpatch
  integer :: ipy
  real :: soil_C
  real :: soil_N
  real :: veg_C
  real :: veg_N
  logical,parameter :: print_on = .false.

  do ipy = 1,cgrid%npolygons

     call compute_C_and_N_storage(cgrid,ipy, soil_C, soil_N, veg_C, veg_N)

     if (print_on) then
        
        print*,"================================================="
        print*,"CBUDGET-INITIAL STORAGE, SOIL_C+VEG_C-NEP, SOIL_C+VEG_C"
        print*,cgrid%cbudget_initialstorage(ipy), soil_C+veg_C-cgrid%cbudget_nep(ipy), soil_C+veg_C
        print*,""
        print*,"SOIL_C+VEG_C,   NEP,   N-BUDGET_INITIAL_STORAGE, SOIL_N + VEG_N"
        print*,soil_C+veg_C,cgrid%cbudget_nep(ipy), cgrid%nbudget_initialstorage(ipy), soil_N+ veg_N
        print*,""
        print*,"================================================="
        
     endif
     
  enddo
  
  return
end subroutine print_C_and_N_budgets

!=====================================================================

subroutine compute_C_and_N_storage(cgrid,ipy, soil_C, soil_N, veg_C, veg_N)

  use ed_state_vars,only:edtype,polygontype,sitetype,patchtype

  use max_dims, only: n_pft
  use pft_coms, only: include_pft, c2n_recruit, c2n_stem, c2n_leaf,  &
       c2n_storage, c2n_slow, c2n_structural

  implicit none

  type(edtype),target :: cgrid
  type(polygontype),pointer :: cpoly
  type(sitetype),pointer :: csite
  type(patchtype),pointer :: cpatch

  integer :: isi,ipa,ico
  integer, intent(in) :: ipy
  real, intent(out) :: soil_C
  real, intent(out) :: soil_N
  real, intent(out) :: veg_C
  real, intent(out) :: veg_N

  integer :: ipft
  real(kind=8) :: area_factor, this_carbon, this_nitrogen
  real(kind=8) :: soil_C8, soil_N8, veg_C8, veg_N8
  
  real(kind=8), parameter :: almostnothing=1.e-30

  ! Initialize C and N pools
  soil_C8 = 0.0
  soil_N8 = 0.0
  veg_C8 = 0.0
  veg_N8 = 0.0

  cpoly => cgrid%polygon(ipy)

  do isi = 1,cpoly%nsites

     csite => cpoly%site(isi)

     do ipa = 1,csite%npatches
          
        cpatch => csite%patch(ipa)

        ! site area times patch area
        area_factor   = dble(cpoly%area(isi)) * dble(csite%area(ipa))
        
        this_carbon   = dble(csite%fast_soil_C(ipa)) + dble(csite%slow_soil_C(ipa)) &
                      + dble(csite%structural_soil_C(ipa))
        this_nitrogen = dble(csite%fast_soil_N(ipa)) + dble(csite%mineralized_soil_N(ipa))  &
                      + dble(csite%slow_soil_C(ipa)) / dble(c2n_slow)                       &
                      + dble(csite%structural_soil_C(ipa)) / dble(c2n_structural)
        
        ! Get soil nitrogen and carbon
        soil_C8 = soil_C8 + area_factor * this_carbon
        soil_N8 = soil_N8 + area_factor * this_nitrogen

        ! Account for carbon/nitrogen in repro arrays
        do ipft = 1, n_pft
           if(include_pft(ipft) == 1)then
              veg_C8 = veg_C8 + dble(csite%repro(ipft,ipa)) * area_factor
              veg_N8 = veg_N8 + dble(csite%repro(ipft,ipa)) / dble(c2n_recruit(ipft)) * area_factor
           endif
        enddo
        
        do ico = 1,cpatch%ncohorts
           
           ! Get the carbon and nitrogen in vegetation.
           veg_C8 = veg_C8 + area_factor * (dble(cpatch%balive(ico)) +   &
                dble(cpatch%bdead(ico)) + dble(cpatch%bstorage(ico))) * dble(cpatch%nplant(ico))
           
           veg_N8 = veg_N8 + area_factor * (dble(cpatch%balive(ico)) /   &
                dble(c2n_leaf(cpatch%pft(ico))) + dble(cpatch%bdead(ico)) / dble(c2n_stem) + dble(cpatch%bstorage(ico)) / &
                dble(c2n_storage)) * dble(cpatch%nplant(ico))
        enddo
     enddo
  enddo
  
  if (abs(soil_C8) < almostnothing) soil_C8 = sign(almostnothing,soil_C8)
  if (abs(soil_N8) < almostnothing) soil_N8 = sign(almostnothing,soil_N8)
  if (abs(veg_C8) < almostnothing)  veg_C8  = sign(almostnothing,veg_C8)
  if (abs(veg_N8) < almostnothing)  veg_N8  = sign(almostnothing,veg_N8)

  soil_C = sngl(soil_C8)
  soil_N = sngl(soil_N8)
  veg_C  = sngl(veg_C8)
  veg_N  = sngl(veg_N8)
  return
end subroutine compute_C_and_N_storage
