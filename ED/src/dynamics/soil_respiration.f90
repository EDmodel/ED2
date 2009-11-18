subroutine soil_respiration(csite,ipa)

  use ed_state_vars,only:sitetype,patchtype
  use soil_coms, only: soil
  use grid_coms, only: nzg
  use pft_coms, only: root_respiration_factor, q, qsw

  implicit none

  type(sitetype),target :: csite
  type(patchtype),pointer :: cpatch
  integer :: ipa,ico
  real :: r_resp_temp_fac
  real :: Lc
  real :: r_resp

  ! This is the temperature dependence of root respiration.  Same for all 
  ! cohorts.

  cpatch => csite%patch(ipa)

  r_resp_temp_fac = 1.0 / (1.0 + exp(0.4   &
       * ( 278.15 - csite%soil_tempk(nzg,ipa) ) ) ) &
       / (1.0 + exp(0.4 * ( csite%soil_tempk(nzg,ipa) - 318.15 ) ) )  &
       * exp( 10.41 - 3000.0/csite%soil_tempk(nzg,ipa) )
  
  do ico = 1,cpatch%ncohorts

     r_resp = root_respiration_factor(cpatch%pft(ico)) * r_resp_temp_fac *   &
          cpatch%balive(ico) * q(cpatch%pft(ico)) / (1.0 + q(cpatch%pft(ico)) + qsw(cpatch%pft(ico)) *   &
          cpatch%hite(ico)) * cpatch%nplant(ico)
                        
     cpatch%root_respiration(ico) = r_resp
     cpatch%mean_root_resp(ico)  = cpatch%mean_root_resp(ico)  + r_resp
     cpatch%today_root_resp(ico) = cpatch%today_root_resp(ico) + r_resp
       
  enddo

  ! Compute soil/temperature modulation of heterotrophic respiration
  call resp_index(csite%ntext_soil(nzg,ipa),csite%soil_tempk(nzg,ipa),  &
       csite%soil_water(nzg,ipa),soil(csite%ntext_soil(nzg,ipa))%slmsts,  &
       csite%A_decomp(ipa))

  ! Compute nitrogen immobilization factor
  call resp_f_decomp(csite,ipa, Lc)

  ! Compute heterotrophic respiration
  call resp_rh(csite,ipa, Lc)

  ! Update averaged variables
  csite%today_A_decomp(ipa) = csite%today_A_decomp(ipa) + csite%A_decomp(ipa)
  csite%today_Af_decomp(ipa) = csite%today_Af_decomp(ipa) +   &
       csite%A_decomp(ipa) * csite%f_decomp(ipa)
  csite%mean_rh(ipa) = csite%mean_rh(ipa) + csite%rh(ipa)

  return
end subroutine soil_respiration

!=====================================================================

subroutine resp_index(nsoil,tempk,theta,slmsts,resp_weight)

  use decomp_coms, only: resp_temperature_increase, resp_opt_water,  &
       resp_water_below_opt, resp_water_above_opt

  implicit none
  
  logical :: LloydTaylor = .true.
  integer :: nsoil
  real :: tempk
  real(kind=8) :: theta
  real :: slmsts
  real :: resp_weight
  real :: temperature_limitation
  real :: water_limitation
  real :: Ws

  ! temperature dependence
  if(LloydTaylor)then 
     !! Use Lloyd and Taylor 1994 temperature dependence
     temperature_limitation = min(1.0, &
          resp_temperature_increase*exp(308.56*(1./56.02-1./(tempk-227.15))))
  else !! use original exponential temperature dependence
     temperature_limitation = min(1.0,exp(resp_temperature_increase *   &
          (tempk-318.15)))
  end if
  ! moisture dependence
  Ws = real(theta/dble(slmsts))



  if(Ws.le.resp_opt_water)then
     water_limitation = exp( (Ws - resp_opt_water) * resp_water_below_opt)
  else
     water_limitation = exp((resp_opt_water-Ws) * resp_water_above_opt)
  endif
  
  ! compute the weight
  resp_weight = temperature_limitation * water_limitation
     
  return
end subroutine resp_index


!===================================================================

subroutine resp_f_decomp(csite,ipa,Lc)

  use ed_state_vars,only: sitetype
  use decomp_coms, only: r_stsc, N_immobil_supply_scale, K1,   &
       n_decomp_lim
  use pft_coms, only: c2n_structural, c2n_slow

  implicit none

  type(sitetype),target :: csite
  integer :: ipa
  real :: N_immobilization_demand
  real, intent(out) :: Lc

  if(csite%structural_soil_C(ipa) > 0.0)then
     if(csite%structural_soil_L(ipa) == csite%structural_soil_C(ipa))then
        Lc = 0.049787 ! = exp(-3.0)
     else
        Lc = exp(-3.0*csite%structural_soil_L(ipa)/csite%structural_soil_C(ipa))
     endif
  else
     Lc=0.0
  endif
  
  if(n_decomp_lim == 1)then
     N_immobilization_demand = csite%A_decomp(ipa) * Lc * K1   &
          * csite%structural_soil_C(ipa)   &
          *((1.0 - r_stsc) / c2n_slow - 1.0 / c2n_structural) 
     
     csite%f_decomp(ipa) = N_immobil_supply_scale * csite%mineralized_soil_N(ipa)   &
          / (N_immobilization_demand + N_immobil_supply_scale *   &
          csite%mineralized_soil_N(ipa))
  else
     ! Option for no plant N limitation
     csite%f_decomp(ipa) = 1.0
  endif

  return
end subroutine resp_f_decomp

!=============================================================

subroutine resp_rh(csite,ipa,Lc)

  use decomp_coms, only: K1, K2, K3, r_fsc, r_ssc, r_stsc, cwd_frac
  use ed_state_vars,only:sitetype

  implicit none

  type(sitetype), target :: csite
  integer :: ipa
  real, intent(in) :: Lc
  real :: fast_C_loss
  real :: structural_C_loss
  real :: slow_C_loss

  ! These have units [kgC/m2/day]
  fast_C_loss = csite%A_decomp(ipa) * K2 * csite%fast_soil_C(ipa)
  structural_C_loss = csite%A_decomp(ipa) * Lc * K1   &
       * csite%structural_soil_C(ipa) * csite%f_decomp(ipa)
  slow_C_loss = csite%A_decomp(ipa) * K3 * csite%slow_soil_C(ipa)

  ! Unit conversion is (kg C)/day  to (umol CO2)/s

  csite%rh(ipa) = 964.5062 * (r_fsc*fast_C_loss + r_stsc*structural_C_loss    &
       + r_ssc*slow_C_loss)
  csite%cwd_rh(ipa) = 964.5062 * (r_stsc*structural_C_loss + r_ssc*slow_C_loss) *   &
       cwd_frac

  return
end subroutine resp_rh

!==========================================================================

subroutine update_C_and_N_pools(cgrid)
  
  use ed_state_vars,only:edtype,polygontype,sitetype
  use decomp_coms, only: K1, K2, K3, r_stsc
  use pft_coms, only: c2n_slow, c2n_structural
  
  implicit none
  
  type(edtype),target :: cgrid
  type(polygontype),pointer :: cpoly
  type(sitetype),pointer    :: csite
  integer :: ipy,isi,ipa
  real :: Lc
  real :: fast_C_loss
  real :: fast_N_loss
  real :: structural_C_loss
  real :: structural_L_loss
  real :: slow_C_input
  real :: slow_C_loss

  do ipy = 1,cgrid%npolygons

     cpoly => cgrid%polygon(ipy)

     do isi = 1,cpoly%nsites
        
        csite => cpoly%site(isi)

        do ipa = 1,csite%npatches

           if(csite%structural_soil_C(ipa) > 0.0)then
              if(csite%structural_soil_L(ipa) == csite%structural_soil_C(ipa))then
                 Lc = 0.049787 ! = exp(-3.0)
              else
                 Lc = exp(-3.0*csite%structural_soil_L(ipa) /   &
                      csite%structural_soil_C(ipa))
              endif
           else
              Lc=0.0
           endif
     
           ! fast pools
           fast_C_loss = csite%today_A_decomp(ipa) * K2 * csite%fast_soil_C(ipa)
           fast_N_loss = csite%today_A_decomp(ipa) * K2 * csite%fast_soil_N(ipa)
!fast_C_loss = 0.0
           
           ! structural pools
           structural_C_loss = csite%today_Af_decomp(ipa) * Lc * K1 *   &
                csite%structural_soil_C(ipa)
           structural_L_loss = csite%today_Af_decomp(ipa) * Lc * K1 *   &
                csite%structural_soil_L(ipa)
!structural_C_loss = 0.0
           
           ! slow pools
           slow_C_input = (1.0 - r_stsc) * structural_C_loss
           slow_C_loss = csite%today_A_decomp(ipa) * K3 * csite%slow_soil_C(ipa)
!slow_C_loss = 0.0
           
           ! mineralized pool
           csite%mineralized_N_input = fast_N_loss + slow_C_loss / c2n_slow
           csite%mineralized_N_loss = csite%total_plant_nitrogen_uptake(ipa) +   &
                csite%today_Af_decomp(ipa) * Lc * K1 *   &
                csite%structural_soil_C(ipa)   &
                * ( (1.0 - r_stsc) / c2n_slow - 1.0 / c2n_structural)
     
           ! all C fluxes have units kgC/m2/day, and we are updating on 
           ! the daily time step.
           csite%fast_soil_C(ipa) = csite%fast_soil_C(ipa) + csite%fsc_in(ipa) -   &
                fast_C_loss
           csite%structural_soil_C(ipa) = csite%structural_soil_C(ipa) +   &
                csite%ssc_in(ipa) - structural_C_loss
           csite%structural_soil_L(ipa) = csite%structural_soil_L(ipa) +   &
                csite%ssl_in(ipa) - structural_L_loss
           csite%slow_soil_C(ipa) = csite%slow_soil_C(ipa) + slow_C_input -  &
                slow_C_loss
           
           ! all N fluxes have units kgN/m2/day, and we are updating on
           ! the daily time step
           csite%fast_soil_N(ipa) = csite%fast_soil_N(ipa) + csite%fsn_in(ipa) -   &
                fast_N_loss
           csite%mineralized_soil_N(ipa) = csite%mineralized_soil_N(ipa) +   &
                csite%mineralized_N_input(ipa) - csite%mineralized_N_loss(ipa)
     
           ! reset average variables
           csite%dmean_A_decomp(ipa) = 0.0
           csite%dmean_Af_decomp(ipa) = 0.0
           
           ! require pools to be >= 0.0
           csite%fast_soil_C(ipa) = max(0.0,csite%fast_soil_C(ipa))
           csite%structural_soil_C(ipa) = max(0.0,csite%structural_soil_C(ipa))
           csite%structural_soil_L(ipa) = max(0.0,csite%structural_soil_L(ipa))
           csite%slow_soil_C(ipa) = max(0.0,csite%slow_soil_C(ipa))
           csite%fast_soil_N(ipa) = max(0.0,csite%fast_soil_N(ipa))
           csite%mineralized_soil_N(ipa) = max(0.0,csite%mineralized_soil_N(ipa))
           
        enddo
     enddo
  enddo
  
  return
end subroutine update_C_and_N_pools
