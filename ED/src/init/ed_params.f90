!==========================================================================================!
!==========================================================================================!
subroutine load_ed_ecosystem_params()

   use max_dims, only: n_pft
   use pft_coms, only: include_pft, include_pft_ag,C2B,frost_mort,include_these_pft,grass_pft
   use disturb_coms,only:min_new_patch_area,num_lu_trans


   implicit none
   integer :: p

   call init_decomp_params()
   call init_ff_coms()
   call init_disturb_params()
   call init_lapse_params()
   call init_can_rad_params()
   call init_hydro_coms()
   call init_soil_coms()
   call init_phen_coms()
   call init_ed_misc_coms()

   ! PLANT FUNCTIONAL TYPES (PFTs):
   ! 1 - C4 grass
   ! 2 - early tropical
   ! 3 - mid tropical
   ! 4 - late tropical
   ! 5 - C3 grass
   ! 6 - northern pines
   ! 7 - southern pines
   ! 8 - late conifers
   ! 9 - early temperate deciduous
   ! 10 - mid temperate deciduous
   ! 11 - late temperate deciduous 
   ! 12 - c3 pasture
   ! 13 - c3 crop (e.g.,wheat, rice, soybean) 
   ! 14 - c4 pasture
   ! 15 - c4 crop (e.g.,corn/maize)

   ! grass PFTs
   grass_pft=huge(1)
   grass_pft(1)=1
   grass_pft(2)=5

   grass_pft(3)=12
   grass_pft(4)=13
   grass_pft(5)=14
   grass_pft(6)=15

   ! Include_pft: flag specifying to whether you want to include a plant functional 
   ! type (1) or whether you want it excluded (0) from the simulation.
   include_pft = 0
   include_pft_ag = 0
   do p=1,n_pft
      if (include_these_pft(p) > 0 .and. include_these_pft(p) <= n_pft) then
         include_pft(include_these_pft(p)) = 1
      end if
   end do

   ! Grasses can grow anywhere, including agricultural patches!
   p=1
   do while (grass_pft(p) > 0 .and. grass_pft(p) <= n_pft)
      if (include_pft(grass_pft(p)) == 1) include_pft_ag(grass_pft(p)) = 1
      p = p+1
   end do
   if (sum(include_pft_ag) == 0) then
!      WHY DO WE REQUIRE THERE TO BE AT LEAST ONE GRASS??  (MCD)
!      MLO - Because pft 1 and 5 used to be the only PFTs allowed in agricultural patches.
!            So when a new agricultural patch is created not having one was causing memory
!            allocation issues. Yeonjoo also mentioned that this may be a problem at the
!            reproduction.
!      call fatal_error ('No grass included in include_these_pft, you should have at least one kind of grass...' &
!                       ,'load_ecosystem_params','ed_params.f90')
      call warning ('No grass included in include_these_pft, you should have at least one kind of grass...' &
                       ,'load_ecosystem_params','ed_params.f90')
   end if

   call init_pft_photo_params()
   call init_pft_resp_params()
   call init_pft_mort_params()
   call init_pft_alloc_params()
   call init_pft_nitro_params()
   call init_pft_leaf_params()
   call init_pft_repro_params()

   call init_pft_derived_params()

   return
end subroutine load_ed_ecosystem_params
!==========================================================================================!
!==========================================================================================!

subroutine init_ed_misc_coms

  use ed_misc_coms,only:burnin,outputMonth, &
       restart_target_year,use_target_year

  burnin = 0 !! number of years to ignore demography when starting a run

  outputMonth = 6 !! month to output annual files

  restart_target_year = 2000 !! year to read when parsing pss/css with multiple years

  use_target_year = 0    !! flag specifying whether to search for a target year in pss/css
  
  return
end subroutine init_ed_misc_coms




!==========================================================================================!
!==========================================================================================!
subroutine init_lapse_params()
!! define defaults for lapse rates

  use met_driver_coms, only: lapse

  lapse%geoht   = 0.0
  lapse%vels    = 0.0
  lapse%atm_tmp = 0.0
  lapse%atm_shv = 0.0
  lapse%prss    = 0.0
  lapse%pcpg    = 0.0
  lapse%atm_co2 = 0.0
  lapse%rlong   = 0.0
  lapse%nir_beam    = 0.0
  lapse%nir_diffuse = 0.0
  lapse%par_beam    = 0.0
  lapse%par_diffuse = 0.0

end subroutine init_lapse_params
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine init_can_rad_params()

use canopy_radiation_coms, only: leaf_reflect_nir,leaf_trans_nir,  &
     leaf_scatter_nir,leaf_reflect_vis_temperate,leaf_trans_vis_temperate, &
     leaf_scatter_vis,leaf_reflect_vis_tropics, leaf_trans_vis_tropics,  &
     diffuse_backscatter_vis, diffuse_backscatter_nir, emis_v, &
     mubar,visible_fraction,visible_fraction_dir,visible_fraction_dif, &
     leaf_reflect_nir,leaf_trans_nir,lai_min,rlong_min,veg_temp_min

use max_dims, only: n_pft
use pft_coms, only: phenology

implicit none

integer :: ipft
real :: leaf_scatter_vis_temperate
real :: leaf_scatter_vis_tropics
real :: diffuse_bscat_vis_temp
real :: diffuse_bscat_vis_trop

mubar     = 1.0 

visible_fraction = 0.45
visible_fraction_dir = 0.43
visible_fraction_dif = 0.52
leaf_reflect_nir = 0.577
leaf_trans_nir = 0.248
lai_min = 1.0e-5


leaf_scatter_nir = leaf_reflect_nir + leaf_trans_nir

leaf_scatter_vis_temperate = leaf_reflect_vis_temperate +   &
     leaf_trans_vis_temperate

leaf_scatter_vis_tropics = leaf_reflect_vis_tropics +   &
     leaf_trans_vis_tropics

diffuse_bscat_vis_temp= (2.0 * leaf_reflect_vis_temperate -   &
     leaf_trans_vis_temperate) / ( 3.0 * leaf_scatter_vis_temperate )

diffuse_bscat_vis_trop = (2.0 * leaf_reflect_vis_tropics -   &
     leaf_trans_vis_tropics) / ( 3.0 * leaf_scatter_vis_tropics )

diffuse_backscatter_nir = (2.0 * leaf_reflect_nir -   &
     leaf_trans_nir) / ( 3.0 * leaf_scatter_nir )

leaf_scatter_vis(1:4) = leaf_scatter_vis_tropics
leaf_scatter_vis(5:11) = leaf_scatter_vis_temperate
leaf_scatter_vis(12:13) = leaf_scatter_vis_temperate
leaf_scatter_vis(14:15) = leaf_scatter_vis_tropics

diffuse_backscatter_vis(1:4) =  diffuse_bscat_vis_trop
diffuse_backscatter_vis(5:11) = diffuse_bscat_vis_temp
diffuse_backscatter_vis(12:13) = diffuse_bscat_vis_temp
diffuse_backscatter_vis(14:15) = diffuse_bscat_vis_trop

emis_v(1) = 0.96d0
emis_v(2:4) = 0.95d0
emis_v(5) = 0.96d0
emis_v(6:8) = 0.97d0
emis_v(9:11) = 0.95d0
emis_v(12:15) = 0.96d0

rlong_min = 50.0
veg_temp_min = 150.0


return
end subroutine init_can_rad_params
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine init_pft_photo_params()

use max_dims,only : n_pft
use pft_coms, only: D0, Vm_low_temp, Vm0, stomatal_slope, cuticular_cond, &
     quantum_efficiency, photosyn_pathway

implicit none


D0 = 0.01 ! same for all PFTs

Vm_low_temp(1:4) = 5.0 ! tropical PFTs
Vm_low_temp(5:13) = 4.7137 ! temperate PFTs
Vm_low_temp(14:15) = 5.0 

Vm0(1) = 12.5
Vm0(2) = 18.8
Vm0(3) = 12.5
Vm0(4) = 6.25
Vm0(5) = 18.3
Vm0(6) = 15.625 * 0.7264
Vm0(7) = 15.625 * 0.7264
Vm0(8) = 6.25 * 0.7264
Vm0(9) = 18.25 * 1.1171
Vm0(10) = 15.625 * 1.1171
Vm0(11) = 6.25 * 1.1171
Vm0(12:13) = 18.3
Vm0(14:15) = 12.5

stomatal_slope(1) = 10.0
stomatal_slope(2:4) = 8.0
stomatal_slope(5:13) = 6.3949
stomatal_slope(14:15) = 10.0 

cuticular_cond(1:5) = 10000.0
cuticular_cond(6:8) = 1000.0
cuticular_cond(9:11) = 20000.0
cuticular_cond(12:15) = 10000.0

quantum_efficiency(1) = 0.06
quantum_efficiency(2:13) = 0.08
quantum_efficiency(14:15) = 0.06

photosyn_pathway(1) = 4
photosyn_pathway(2:13) = 3
photosyn_pathway(14:15) = 4

return
end subroutine init_pft_photo_params
!==========================================================================================!
!==========================================================================================!

subroutine init_decomp_params()
  
  use decomp_coms,only:resp_opt_water,resp_water_below_opt, &
       resp_water_above_opt,resp_temperature_increase, &
       N_immobil_supply_scale,cwd_frac,r_fsc,r_stsc,r_ssc,K1,K2,K3

  resp_opt_water = 0.8938
  resp_water_below_opt = 5.0786
  resp_water_above_opt = 4.5139
  resp_temperature_increase = 0.0757
  N_immobil_supply_scale = 40.0 / 365.2425
  cwd_frac = 0.2
  r_fsc=1.0  
  r_stsc=0.3 
  r_ssc=1.0  
  K1=4.5 / 365.2425
  K2=11.0 / 365.2425
  K3=100.2 / 365.2425

  return

end subroutine init_decomp_params




!==========================================================================================!
!==========================================================================================!
subroutine init_pft_resp_params()

use pft_coms, only: growth_resp_factor, leaf_turnover_rate,   &
     dark_respiration_factor, storage_turnover_rate,  &
     root_respiration_factor
use decomp_coms, only: f_labile

implicit none

growth_resp_factor(1:5) = 0.333
growth_resp_factor(6:8) = 0.4503
growth_resp_factor(9:11) = 0.0
growth_resp_factor(12:15) = 0.333

leaf_turnover_rate(1) = 2.0
leaf_turnover_rate(2) = 1.0
leaf_turnover_rate(3) = 0.5
leaf_turnover_rate(4) = 0.333
leaf_turnover_rate(5) = 2.0
leaf_turnover_rate(6:8) = 0.333
leaf_turnover_rate(9:11) = 0.0
leaf_turnover_rate(12:15) = 2.0

dark_respiration_factor(1) = 0.04
dark_respiration_factor(2:4) = 0.02
dark_respiration_factor(5) = 0.04
dark_respiration_factor(6:11) = 0.02
dark_respiration_factor(12:15) = 0.04

storage_turnover_rate(1:8) = 0.0
storage_turnover_rate(9:11) = 0.6243
storage_turnover_rate(12:15) = 0.0

root_respiration_factor = 0.528

!!!! Added f_labile  [[MCD]]
f_labile(1:5) = 1.0
f_labile(6:11) = 0.79
f_labile(12:15) = 1.0

return
end subroutine init_pft_resp_params
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine init_pft_mort_params()

use pft_coms, only: mort1, mort2, mort3, seedling_mortality, treefall_s_gtht,  &
     treefall_s_ltht, plant_min_temp,frost_mort

implicit none


frost_mort = 3.0


mort1(1:4) = 10.0
mort1(5:11) = 1.0
mort1(12:13) = 1.0
mort1(14:15) = 10.0

mort2 = 20.0

mort3(1) = 0.037
mort3(2) = 0.037
mort3(3) = 0.019
mort3(4) = 0.0
mort3(5) = 0.066
mort3(6) = 0.0033928
mort3(7) = 0.0043
mort3(8) = 0.0023568
mort3(9) = 0.006144
mort3(10) = 0.003808
mort3(11) = 0.00428
mort3(12:13) = 0.066
mort3(14:15) = 0.037

seedling_mortality = 0.95

treefall_s_gtht = 0.0

treefall_s_ltht(1) = 0.25
treefall_s_ltht(2:4) = 0.1
treefall_s_ltht(5) = 0.25
treefall_s_ltht(6:11) = 0.1
treefall_s_ltht(12:15) = 0.25

plant_min_temp(1:4) = 0.0
plant_min_temp(5:6) = -80.0
plant_min_temp(7) = -10.0
plant_min_temp(8) = -60.0
plant_min_temp(9) = -80.0
plant_min_temp(10:11) = -20.0
plant_min_temp(12:13) = -80.0
plant_min_temp(14:15) = 0.0

return
end subroutine init_pft_mort_params
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine init_pft_alloc_params()

use pft_coms, only: rho, SLA, q, qsw, hgt_min, b1Ht, b2Ht, b1Bs,  &
     b2Bs, b1Bl, b2Bl, C2B, leaf_turnover_rate, hgt_ref

implicit none

C2B    = 2.0               !  Carbon-to-biomass ratio of plant tissues.

rho(1:2) = 0.53
rho(3) = 0.71
rho(4) = 0.9
rho(5:11) = 0.0
rho(12:13) = 0.0
rho(14:15) = 0.53


SLA(1:4) = 10.0**((2.4-0.46*log10(12.0/leaf_turnover_rate(1:4)))) * C2B * 0.1
SLA(5) = 22.0
SLA(6) = 6.0
SLA(7) = 9.0
SLA(8) = 10.0
SLA(9) = 30.0
SLA(10) = 24.2
SLA(11) = 60.0
SLA(12:13) = 22.0
SLA(14:15) =  10.0**((2.4-0.46*log10(12.0/leaf_turnover_rate(13)))) * C2B * 0.1

q(1:5) = 1.0
q(6:8) = 0.3463
q(9:11) = 1.1274
q(12:15) = 1.0


!    Using the wrong qsw for mid-latitude PFTs (5-13), since the other parameters need to 
! be optimized with the fixed version. I assumed 12 and 13 to be mid-latitudes, and 
! 14 and 15 to be tropical, is that correct?
qsw(1:4)    = SLA(1:4)   / (3900.0*2.0/1000.0)
qsw(14:15)  = SLA(14:15) / (3900.0*2.0/1000.0)
qsw(5:13)    = SLA(5:13) / 3900.0 !KIM - ED1/ED2 codes and Moorcroft et al.'re wrong!

hgt_min = 1.5
hgt_min(5) = 0.2
hgt_min(12:13) = 0.2

hgt_ref = 0.0
hgt_ref(6:11) = 1.3

b1Ht(1:4) = 0.0
b1Ht(5) = 0.4778
b1Ht(6) = 27.14
b1Ht(7) = 27.14
b1Ht(8) = 22.79
b1Ht(9) = 22.6799
b1Ht(10) = 25.18
b1Ht(11) = 23.3874
b1Ht(12:13) = 0.4778
b1Ht(14:15) = 0.0

b2Ht(1:4) = 0.0
b2Ht(5) = -0.75
b2Ht(6) = -0.03884
b2Ht(7) = -0.03884
b2Ht(8) = -0.04445 
b2Ht(9) = -0.06534
b2Ht(10) = -0.04964
b2Ht(11) = -0.05404
b2Ht(12:13) = -0.75
b2Ht(14:15) = 0.0

b1Bl(1:4) = 0.0
b1Bl(5) = 0.08
b1Bl(6) = 0.024
b1Bl(7) = 0.024
b1Bl(8) = 0.0454
b1Bl(9) = 0.0129
b1Bl(10) = 0.048
b1Bl(11) = 0.017
b1Bl(12:13) = 0.08
b1Bl(14:15) = 0.0

b2Bl(1:4) = 0.0
b2Bl(5) = 1.0
b2Bl(6) = 1.899
b2Bl(7) = 1.899
b2Bl(8) = 1.6829
b2Bl(9) = 1.7477
b2Bl(10) = 1.455
b2Bl(11) = 1.731
b2Bl(12:13) = 1.0
b2Bl(14:15) = 0.0

b1Bs(1:4) = 0.0 
b1Bs(5) = 1.0e-5
b1Bs(6) = 0.147
b1Bs(7) = 0.147
b1Bs(8) = 0.1617
b1Bs(9) = 0.02648
b1Bs(10) = 0.1617
b1Bs(11) = 0.235
b1Bs(12:13) = 1.0e-5
b1Bs(14:15) = 0.0 

b2Bs(1:4) = 0.0
b2Bs(5) = 1.0
b2Bs(6) = 2.238
b2Bs(7) = 2.238
b2Bs(8) = 2.1536
b2Bs(9) = 2.95954
b2Bs(10) = 2.4572
b2Bs(11) = 2.2518
b2Bs(12:13) = 1.0
b2Bs(14:15) = 0.0

return
end subroutine init_pft_alloc_params
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine init_pft_nitro_params()

use pft_coms, only: c2n_leaf, Vm0, SLA, water_conductance, &
     c2n_slow,c2n_structural,c2n_storage,c2n_stem,l2n_stem, &
     C2B,agf_bs,plant_N_supply_scale

implicit none

c2n_slow       = 10.0 ! Carbon to Nitrogen ratio, slow pool.
c2n_structural = 150.0 ! Carbon to Nitrogen ratio, structural pool.
c2n_storage    = 150.0 ! Carbon to Nitrogen ratio, storage pool.
c2n_stem       = 150.0 ! Carbon to Nitrogen ratio, structural stem.
l2n_stem       = 150.0 ! Carbon to Nitrogen ratio, structural stem.


agf_bs = 0.7  ! fraction of structural stem that is assumed to be above ground.
plant_N_supply_scale = 0.5 

c2n_leaf = 1000.0 / ((0.11289 + 0.12947 * Vm0) * SLA)
c2n_leaf(6) = 1000.0 / ((0.11289 + 0.12947 * 15.625) * SLA(6))
c2n_leaf(7) = 1000.0 / ((0.11289 + 0.12947 * 15.625) * SLA(7))
c2n_leaf(8) = 1000.0 / ((0.11289 + 0.12947 * 6.25) * SLA(8))
c2n_leaf(9) = 1000.0 / ((0.11289 + 0.12947 * 18.25) * SLA(9))
c2n_leaf(10) = 1000.0 / ((0.11289 + 0.12947 * 15.625) * SLA(10))
c2n_leaf(11) = 1000.0 / ((0.11289 + 0.12947 * 6.25) * SLA(11))



water_conductance(1:4) = 0.001904
water_conductance(5:11) = 0.00476
water_conductance(12:13) = 0.00476
water_conductance(14:15) = 0.001904

return
end subroutine init_pft_nitro_params
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine init_pft_leaf_params()

use pft_coms, only: phenology, clumping_factor, leaf_width

implicit none

phenology(1:5) = 1
phenology(6:8) = 0
phenology(9:11) = 2
phenology(12:15) = 1

clumping_factor(1) = 1.0
clumping_factor(2:4) = 0.735
clumping_factor(5) = 0.84
clumping_factor(6:8) = 0.735
clumping_factor(9:11) = 0.84
clumping_factor(12:13) = 0.84
clumping_factor(14:15) = 1.0

leaf_width(1:4) = 0.2
leaf_width(5:11) = 0.05
leaf_width(12:13) = 0.05
leaf_width(14:15) = 0.2

return
end subroutine init_pft_leaf_params
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine init_pft_repro_params()

use pft_coms, only: r_fract, seed_rain, nonlocal_dispersal, repro_min_h

implicit none

r_fract(1) = 1.0
r_fract(2:4) = 0.3
r_fract(5) = 1.0
r_fract(6:11) = 0.3
r_fract(12:15) = 1.0

seed_rain = 0.01

nonlocal_dispersal(1:5) = 1.0
nonlocal_dispersal(6:7) = 0.766
nonlocal_dispersal(8) = 0.001
nonlocal_dispersal(9) = 1.0
nonlocal_dispersal(10) = 0.325
nonlocal_dispersal(11) = 0.074
nonlocal_dispersal(12:15) = 1.0

repro_min_h(1:5) = 0.0
repro_min_h(6:11) = 5.0
repro_min_h(12:15) = 0.0

return
end subroutine init_pft_repro_params
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine init_pft_derived_params()

use pft_coms, only: root_turnover_rate, c2n_leaf, max_dbh, b1Ht, b2Ht,  &
     hgt_min, q, qsw, c2n_recruit, c2n_stem,hgt_ref
use decomp_coms, only: f_labile
use max_dims, only: n_pft

implicit none

integer :: ipft
real :: dbh,balive,bleaf,bdead
real, external :: h2dbh,dbh2bl,dbh2bd

root_turnover_rate(1) = 2.0
root_turnover_rate(2) = 1.0
root_turnover_rate(3) = 0.5
root_turnover_rate(4:5) = 0.333
!root_turnover_rate(6:11) = 5.1* c2n_leaf(10) / c2n_leaf(6:11)  ! leaves and fine roots have the same c2n.
root_turnover_rate(6) = 3.927218
root_turnover_rate(7) = 4.117847
root_turnover_rate(8) = 3.800132
root_turnover_rate(9) = 5.772506
root_turnover_rate(10) = 5.083700
root_turnover_rate(11) = 5.070992
root_turnover_rate(12:13) = 0.333
root_turnover_rate(14:15) = 2.0

max_dbh(1) = 0.498
max_dbh(2:4) = 68.31
max_dbh(5) = 0.498
max_dbh(6:11) = log(1.0-(0.999*b1Ht(6:11)-hgt_ref(6:11))/b1Ht(6:11))/b2Ht(6:11)
max_dbh(12:15) = 0.498

do ipft = 1,n_pft
   dbh = h2dbh(hgt_min(ipft),ipft)
   bleaf = dbh2bl(dbh,ipft)
   bdead = dbh2bd(dbh,hgt_min(ipft),ipft)
   balive = bleaf * (1.0 + q(ipft) + qsw(ipft) * hgt_min(ipft))
   c2n_recruit(ipft) = (balive + bdead) / (balive * (f_labile(ipft) /   &
        c2n_leaf(ipft) + (1.0 - f_labile(ipft)) / c2n_stem) +   &
        bdead/c2n_stem)
enddo

end subroutine init_pft_derived_params
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine init_disturb_params

  use disturb_coms,only:patch_dynamics,min_new_patch_area, &
       treefall_hite_threshold,treefall_age_threshold, &
       forestry_on,agriculture_on,plantation_year,plantation_rotation, &
       mature_harvest_age,fire_dryness_threshold,fire_parameter

  implicit none
  
  patch_dynamics = 1
  
  min_new_patch_area = 0.005
  
  treefall_hite_threshold = 10.0  !  Only trees above this height create a gap when they fall.
  
  treefall_age_threshold = 0.0  !  Minimum patch age for treefall disturbance.  
  
  forestry_on = 0  ! Set to 1 if to do forest harvesting.
  
  agriculture_on = 0  ! Set to 1 if to do agriculture.
  
  plantation_year = 1960 ! Earliest year at which plantations occur
  
  plantation_rotation = 25.0 ! Number of years that a plantation requires to reach maturity
  
  mature_harvest_age = 50.0 ! Years that a non-plantation patch requires to reach maturity
  
  fire_dryness_threshold = 0.2  !  (meters) Fire can occur if total soil water falls below this threshold.
  
  fire_parameter = 1.0  ! Dimensionless parameter controlling speed of fire spread.
  
  return

end subroutine init_disturb_params
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine init_hydro_coms

  use hydrology_coms,only:useTOPMODEL,useRUNOFF,HydroOutputPeriod, &
       MoistRateTuning,MoistSatThresh,Moist_dWT,FracLiqRunoff, &
       GrassLAIMax,inverse_runoff_time
  use misc_coms, only: ied_init_mode

  implicit none

  if(ied_init_mode == 3)then
     ! Signifies a restart from an ED2 history file
     useTOPMODEL = 1
     useRUNOFF   = 0
  else
     useTOPMODEL = 0
     useRUNOFF = 0
  endif

  HydroOutputPeriod = 96 !! multiples of dtlsm

  MoistRateTuning = 1.0    

  MoistSatThresh = 0.95     
  
  Moist_dWT = 2.0 

  FracLiqRunoff = 0.5
  
  GrassLAImax = 4.0

  inverse_runoff_time = 0.1

  return
end subroutine init_hydro_coms
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine init_soil_coms


  use soil_coms,only:water_stab_thresh,dewmax

  implicit none

  water_stab_thresh = 60.0

  dewmax = 3.0e-5

  return

end subroutine init_soil_coms
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine init_phen_coms

  use phenology_coms,only: retained_carbon_fraction, &
       theta_crit,dl_tr,st_tr1,st_tr2,phen_a,phen_b,phen_c

  implicit none

  retained_carbon_fraction = 0.5
  theta_crit = 0.2
  dl_tr = 655.0
  st_tr1 = 284.3
  st_tr2 = 275.15
  
  phen_a = -68.0   
  phen_b = 638.0    
  phen_c = -0.01    


  return
end subroutine init_phen_coms
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine init_ff_coms

  use fusion_fission_coms,only:min_recruit_size,min_dbh_class,maxdbh,min_hgt_class, &
       fusetol,fusetol_h,lai_fuse_tol,lai_tol,ntol,profile_tol,max_patch_age

  implicit none

  min_recruit_size  = 1.0e-3
  
  min_dbh_class = 0.0  
  
  maxdbh = 200.0 
  
  min_hgt_class = 0.0
  
  fusetol = 0.4
  
  fusetol_h = 0.5
  
  lai_fuse_tol = 0.8
  
  lai_tol = 1.0
  
  ntol = 0.001
  
  profile_tol = 0.2

  max_patch_age = 500.0
  

  return

end subroutine init_ff_coms
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine overwrite_with_xml_config(thisnode)
   !!! PARSE XML FILE
   use max_dims, only: n_pft
   use misc_coms, only: iedcnfgf
   
   implicit none
   integer, intent(in) :: thisnode
   integer             :: max_pft_xml
   logical             :: iamhere

   if (iedcnfgf /= '') then
      inquire (file=iedcnfgf,exist=iamhere)
      if (iamhere) then

         !! FIRST, determine number of pft's defined in xml file
         call count_pft_xml_config(trim(iedcnfgf),max_pft_xml)
         if(max_pft_xml > n_pft) then

            write(unit=*,fmt='(a)') '*********************************************'
            write(unit=*,fmt='(a)') '**                                         **'
            write(unit=*,fmt='(a)') '**  Number of PFTs required by XML Config  **'
            write(unit=*,fmt='(a)') '**  exceeds the memory available           **'
            write(unit=*,fmt='(a)') '**                                         **'
            write(unit=*,fmt='(a)') '**  Please change n_pft in Module max_dims **'
            write(unit=*,fmt='(a)') '**  and recompile                          **'
            write(unit=*,fmt='(a)') '**                                         **'
            write(unit=*,fmt='(a)') '*********************************************'
            call fatal_error('Too many PFTs','overwrite_with_xml_config','ed_params.f90')
         end if

         !! SECOND, update parameter defaults from XML
         call read_ed_xml_config(trim(iedcnfgf))

         !! THIRD, reset any values based on xml

         !! FINALLY, write out copy of settings
         call write_ed_xml_config()
      elseif (thisnode == 1) then
            write(unit=*,fmt='(a)') '*********************************************'
            write(unit=*,fmt='(a)') '**               WARNING!                  **'
            write(unit=*,fmt='(a)') '**                                         **'
            write(unit=*,fmt='(a)') '**    XML file wasn''t found. Using default **'
            write(unit=*,fmt='(a)') '** parameters in ED.                       **'
            write(unit=*,fmt='(a)') '** (You provided '//trim(iedcnfgf)//').'
            write(unit=*,fmt='(a)') '**                                         **'
            write(unit=*,fmt='(a)') '*********************************************'
      end if
   end if  !! end XML
   return
end subroutine overwrite_with_xml_config
!==========================================================================================!
!==========================================================================================!
