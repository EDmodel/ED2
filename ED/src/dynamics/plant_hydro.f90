!==========================================================================================!
!==========================================================================================!
! MODULE: ED_STATE_VARS
!> \brief Calculations of plant hydrodynamics at ED2 timestep, including various
!> utils for plant hydraulics
!> \details Util functions to perform conversions between psi, rwc and water_int
!> \author Xiangtao Xu, 26 Jan. 2018
!==========================================================================================!
!==========================================================================================!
module plant_hydro

   contains
!   !=======================================================================================!
!   !=======================================================================================!
!   !      Driver of plant hydrodynamics. This subroutine tracks changes of leaf/stem water !
!   ! water potential as well as water flows from root to stem and from stem to canopy.     !
!   ! This subroutine currently works at DTLSM level, should be integrated into
!   ! RK4 scheme in the future.   [XXT]
!   !---------------------------------------------------------------------------------------!
!   subroutine update_plant_hydrodynamics(csite,ipa,ntext_soil) 
!      use ed_state_vars   , only : sitetype             & ! structure
!                                 , patchtype            ! ! structure
!      use ed_misc_coms    , only : dtlsm                & ! intent(in)
!                                 , frqsum               & ! intent(in)        
!                                 , current_time         ! ! intent(in) 
!      use soil_coms       , only : dslzi8               &
!                                 , dslz8                &
!                                 , nzg                  &
!                                 , soil8                &
!                                 , slz                  &
!                                 , dslz                 &
!                                 , dslzi                &
!                                 , soil                    ! ! intent(in)
!      use consts_coms     , only : wdns8                &
!                                 , wdnsi8               &
!                                 , wdns                 &
!                                 , wdnsi                &
!                                 , pi1                  
!      use rk4_coms        , only : rk4patchtype         &
!                                 , rk4site              &
!                                 , tiny_offset
!      use therm_lib      , only  : uextcm2tl            & ! function
!                                 , tl2uint              ! ! function
!      use therm_lib8     , only  : tl2uint8             ! ! function
!      use ed_therm_lib   , only  : calc_veg_hcap        ! ! subroutine
!      use pft_coms       , only  : Ks_stem              &
!                                 , Ks_stem_b            &
!                                 , psi50                &
!                                 , Cap_leaf             &
!                                 , Cap_stem             &
!                                 , SRA                  &
!                                 , rho                  &
!                                 , root_beta            &
!                                 , xylem_fraction       &
!                                 , wat_dry_ratio_grn    &
!                                 , wat_dry_ratio_ngrn   &
!                                 , leaf_rwc_min         &
!                                 , wood_rwc_min         &
!                                 , vessel_curl_factor   &
!                                 , C2B                  &
!                                 , is_grass             
!      implicit none
!      !----- Arguments --------------------------------------------------------------------!
!      type(sitetype)            , target      :: csite
!      integer                   , intent(in)  :: ipa
!      integer    ,dimension(nzg), intent(in)  :: ntext_soil
!
!      !----- Local Vars  --------------------------------------------------------------------!
!      type(patchtype)       , pointer       :: cpatch
!      !----- temporary state variables for water flow within plant
!      real(kind=8)                          :: transp               !transpiration
!      real(kind=8)                          :: J_sr                 !soil-root flow
!      real(kind=8)                          :: J_rl                 !root-leaf flow
!      real(kind=8)                          :: org_psi_leaf         !original leaf water potential
!      real(kind=8)                          :: org_psi_stem         !original stem water potential
!      real(kind=8)                          :: exp_term             !exponent term
!      !----- variables necessary to calcualte water flow
!      real(kind=8)                          :: ap
!      real(kind=8)                          :: bp
!      real(kind=8)                          :: stem_cond            !stem conductance
!      real(kind=8)                          :: c_leaf               !leaf water capacitance
!      real(kind=8)                          :: c_stem               !stem water capacitance
!      real(kind=8)                          :: RAI                  !root area index
!      real(kind=8)                          :: cohort_crown_area
!      real(kind=8)                          :: wgpfrac              !relative soil water
!      real(kind=8)                          :: wflx                 !water flux
!      real(kind=8)                          :: qflx                 !internal energy flux
!      real(kind=8)                          :: soil_wflx            !soil water flux
!      real(kind=8)                          :: soil_qflx            !soil internal energy flux
!      real(kind=8)                          :: soil_water_cond
!      real(kind=8),dimension(nzg)           :: layer_psi
!      real(kind=8)                          :: weighted_soil_psi
!      real(kind=8)                          :: weighted_soil_cond
!      real(kind=8)                          :: above_layer_depth
!      real(kind=8)                          :: current_layer_depth
!      real(kind=8)                          :: total_water_supply
!      real(kind=8)      ,dimension(nzg)     :: layer_water_supply
!      !----- variables to deal with energy budget of plant water absorption
!      real(kind=8)                          :: dsw                  !change of soil water
!      real(kind=8)                          :: org_soil_tempk 
!      real(kind=4)                          :: rwc_err
!      real(kind=4)                          :: psi_pot              ! potential psi
!      !----- variables for loops
!      integer                               :: ico
!      integer                               :: k
!      integer                               :: nsoil
!      integer                               :: ipft
!      !--------------- other
!      logical,parameter                     :: quality_check = .false.
!      logical                               :: update_psi
!      integer,parameter                     :: check_ico = 1
!      logical                               :: small_tree_flag
!      real(kind=4), parameter               :: soilcp_MPa = -3.1 ! keep consistent with ed_params.f90
!      !----- External function ------------------------------------------------------------!
!      real(kind=4), external                :: sngloff           ! Safe dble 2 single precision
!      !----- Locally saved variables. -----------------------------------------------------!
!      real        , save              :: dtlsm_o_frqsum
!      logical     , save              :: first_time = .true.
!      !------------------------------------------------------------------------------------!
!
!      !----- Assign the constant scaling factor. ------------------------------------------!
!      if (first_time) then
!         first_time     = .false.
!         dtlsm_o_frqsum = dtlsm / frqsum
!      end if
!      !------------------------------------------------------------------------------------!
!      !------------ Check whether it is time to updated psi-----------------------------------------
!      ! We update the dmin/dmax psi, one timestep after new_day.
!      ! Because we need to conserve the value for phenology_driv and output
!      !------------------------------------------------------------
!      update_psi       = current_time%time == dtlsm
!
!      cpatch => csite%patch(ipa)
!
!      !----------------------------------------------------------------------
!      ! Preparation for plant hydrodynamics
!      !-----------------------------------------------------------------------
!      
!      ! Calculate layer_psi, soil water potential of each soil layer
!      do k = 1,nzg
!        nsoil = ntext_soil(k)
!        
!        !get relative soil moisture
!        wgpfrac = min(1.0,dble(csite%soil_water(k,ipa)) &
!                        * dble(csite%soil_fracliq(k,ipa)) / soil8(nsoil)%slmsts)
!
!        !To ensure numerical stability, we assign a very large negative number
!        !to soil layer water potential if soil moisture is very low.
!        if (wgpfrac < 1e-3) then
!            layer_psi(k) = -1e6
!        else
!            layer_psi(k) = soil8(nsoil)%slpots / wgpfrac ** soil8(nsoil)%slbs
!        endif
!
!      enddo
!      
!
!      !----------------------------------------------------------------------
!      ! Update plant hydrodynamics
!      ! The regular slover assumes stem water pool is way larger than the leaf water
!      ! pool. In cases where leaf water pool is of similar magnitude to stem
!      ! water pool (small seedlings and grasses), leaf water potential is forced
!      ! to be the same as stem water potential to maintain numerical stability
!      !-----------------------------------------------------------------------
!
!
!      !Loop over all the cohorts
!      cohortloop:do ico = 1,cpatch%ncohorts
!        ipft = cpatch%pft(ico)
!      
!        !1. get transpiration rate of the last time step
!        if (cpatch%leaf_resolvable(ico)) then
!            transp = (cpatch%fs_open(ico) * cpatch%psi_open(ico)              &
!                    + (1.0d0 - cpatch%fs_open(ico)) *                               &
!                    cpatch%psi_closed(ico)) * cpatch%lai(ico) / cpatch%nplant(ico)
!        else
!            transp = 0.
!        endif
!
!        ! the unit of transp is kg/sec
!          
!
!        !2. update leaf, stem psi and water flows
!        
!        ! First, check the relative magnitude of leaf and stem capacitance
!        ! If it is a small tree/grass, force psi_leaf to be the same as psi_stem
!        c_leaf = Cap_leaf(ipft) * cpatch%bleaf(ico) * C2B
!        c_stem = Cap_stem(ipft) * & ! kg H2O kg m-1
!                 (cpatch%broot(ico) + cpatch%bdead(ico) * xylem_fraction(ipft)) * C2B
!        small_tree_flag = c_leaf > (c_stem / 2.)
!
!
!        !2.1 update leaf psi
!        if (small_tree_flag) then
!            ! 2.1.1
!            ! small tree, leaf_psi is the same as stem_psi
!            ! use c_leaf+c_stem as the total capacitance to solve stem_psi
!            c_stem = c_stem + c_leaf
!
!            org_psi_leaf = cpatch%psi_leaf(ico)
!            J_rl = transp -                                         &
!                   (cpatch%psi_leaf(ico) - cpatch%psi_stem(ico))    &
!                   * c_leaf / dble(dtlsm)                           
!                   ! water necessary to bring psi_leaf to be the same as
!                   ! psi_stem
!            ! For the convenience of calculation later, we attribute J_rl as
!            ! transp here...
!        else
!
!            !2.1.2 Regular case, big trees
!            !    update leaf psi while assuming stem psi is constant....
!            !    only calcualte the leaf psi if it is not grass. For grass, we treat
!            !    leaf_psi the same as stem_psi
!        
!
!            stem_cond =  Ks_stem(ipft)   &                                                    
!                         *  1 / (1 + (cpatch%psi_stem(ico) / psi50(ipft)) ** Ks_stem_b(ipft))  &  !cavitation effect
!                         * xylem_fraction(ipft) * (pi1 * (cpatch%dbh(ico) / 200.) ** 2)         &  !conducting area     m2
!                         / (cpatch%hite(ico) * vessel_curl_factor(ipft))                          !conducting length   m
!            
!            if (c_leaf > 0. .and. stem_cond > 0.) then
!                ! if there are leaves and the tree has at least some level of
!                ! conductivity
!                ap = - stem_cond / c_leaf
!                ! the unit of ap is s-1
!
!                bp = ((cpatch%psi_stem(ico) - cpatch%hite(ico)) * stem_cond -   &                       
!                         transp)                                                &
!                         / c_leaf                                                           
!                ! the unit of bp is m s-1
!
!                ! calculate new psi_leaf
!                org_psi_leaf = cpatch%psi_leaf(ico)
!                
!                exp_term = ap * dtlsm
!                if (ap * dtlsm < - 40.) then
!                    exp_term = 0.0d0
!                else
!                    exp_term = exp(ap * dble(dtlsm))
!                endif
!
!                cpatch%psi_leaf(ico) =                   &
!                    sngloff(((ap * org_psi_leaf + bp) *  &!  m s-1
!                            exp_term - bp) / ap,tiny_offset)
!
!                ! calculate the average sapflow rate from stem to leaf within the
!                ! time step
!                J_rl = (cpatch%psi_leaf(ico) - org_psi_leaf) * c_leaf / dtlsm + &
!                        transp ! kgH2O s-1
!            elseif (c_leaf > 0.) then
!                ! If there is no conductance in the trunk but the cohort still
!                ! has some leaves...
!
!                ! changes of psi_leaf is solely due to transpiration
!                    
!                org_psi_leaf = cpatch%psi_leaf(ico)
!                cpatch%psi_leaf(ico) = sngloff(org_psi_leaf -              &
!                                       transp * dtlsm  / c_leaf,tiny_offset)
!                J_rl = 0.
!
!            else
!                ! If there is no leaves...
!                ! psi_leaf is set to be equal to psi_stem minus the gravitational
!                ! effect
!                org_psi_leaf = cpatch%psi_leaf(ico)
!                ! only include gravitational effect
!                cpatch%psi_leaf(ico) = cpatch%psi_stem(ico) - cpatch%hite(ico) 
!                J_rl = 0.
!
!            endif
!
!!            ! J_rl can overchrage wood storage when leaf water storage is smiliar to stem water
!!            ! storage. This scenario can happen for small seedlings
!!            ! In this case, we set J_rl to 0.
!!
!            if (J_rl < 0.) then
!                ! downward sapflow
!                wflx = -J_rl * dble(dtlsm)
!
!                psi_pot = cpatch%psi_stem(ico) +             &
!                          sngloff(wflx / c_stem, tiny_offset)
!
!                if (psi_pot > -1e-3) then
!                    J_rl = 0.
!                
!                    cpatch%psi_leaf(ico) = sngloff(org_psi_leaf                       &
!                                                   +(J_rl - transp) * dtlsm  / c_leaf &
!                                                   ,tiny_offset)
!                endif
!            endif
!!
!            if (J_rl > 0. .and. cpatch%wood_rwc(ico) < wood_rwc_min(ipft)) then
!                ! wood_rwc is too low to support upward sapflow
!                J_rl = 0.
!                ! modify psi_leaf
!                cpatch%psi_leaf(ico) = sngloff(org_psi_leaf -              &
!                                       transp * dtlsm  / c_leaf,tiny_offset)
!            elseif (J_rl < 0. .and. cpatch%leaf_rwc(ico) < leaf_rwc_min(ipft)) then
!                ! leaf_rwc is too low to support downward sapflow
!                J_rl = 0.
!                ! modify psi_leaf
!                cpatch%psi_leaf(ico) = sngloff(org_psi_leaf -              &
!                                       transp * dtlsm  / c_leaf,tiny_offset)
!                
!            endif
! 
!        endif
!
!        ! For debugging purpose
!        if(isnan(J_rl)) then
!            write (unit=*,fmt='(80a)')         ('=',k=1,80)
!            write (unit=*,fmt='(a)')           'Sapflow is NaN!!'
!            write (unit=*,fmt='(a,1x,i9)')   ' + ico:                 ',ico
!            write (unit=*,fmt='(a,1x,es12.4)')   ' + LEAF CAPACITANCE:    ',c_leaf
!            write (unit=*,fmt='(a,1x,es12.4)')   ' + LAI:                 ',cpatch%lai(ico)
!            write (unit=*,fmt='(a,1x,es12.4)')   ' + NPLANT:              ',cpatch%nplant(ico)
!            write (unit=*,fmt='(a,1x,es12.4)')   ' + TRANSPIRATION:       ',transp
!            write (unit=*,fmt='(a,1x,es12.4)')   ' + STEM CONDUCTANCE:    ',stem_cond
!            write (unit=*,fmt='(a,1x,es12.4)')   ' + ORG LEAF PSI:        ',org_psi_leaf
!            write (unit=*,fmt='(a,1x,es12.4)')   ' + ORG STEM PSI:        ',cpatch%psi_stem(ico)
!           
!            call fatal_error('Sapflow is wrong','update_plant_hydrodynamics'               &
!                            &,'plant_hydrodynamics.f90')
!        endif
!
!        
!        !2.2 update stem psi while assuming soil psi is constant....
!
!        weighted_soil_psi  = 0.
!        weighted_soil_cond = 0.
!        layer_water_supply = 0.
!        total_water_supply = 0.
!
!        ! loop over all soil layers to get the aggregated water conductance
!        do k = cpatch%krdepth(ico),nzg
!            current_layer_depth = -slz(k)
!            if (k+1 .le. nzg) then
!                above_layer_depth = -slz(k+1)
!            else
!                above_layer_depth = 0.0
!            endif
!                                       
!            !  Calculate RAI in each layer
!            !  Based on Katul et al. 2003 PCE
!            cohort_crown_area = cpatch%crown_area(ico) / cpatch%nplant(ico)
!
!            if (cohort_crown_area == 0.) then
!                RAI = 0.
!            else
!                RAI = cpatch%broot(ico) * & !kgC
!                    (root_beta(ipft) ** (above_layer_depth / -slz(cpatch%krdepth(ico))) - &
!                    root_beta(ipft) ** (current_layer_depth / -slz(cpatch%krdepth(ico))) ) * &
!                    SRA(ipft) / &  !    m2/kgC
!                    (2.0 * cohort_crown_area)     ! m2
!            endif
!
!            ! The unit of RAI is m2/m2
!
!            !  Calculate soil water conductance
!            nsoil = ntext_soil(k)
!            wgpfrac = min(1.0,dble(csite%soil_water(k,ipa)) &
!                        * dble(csite%soil_fracliq(k,ipa)) / soil8(nsoil)%slmsts)
!            
!            soil_water_cond = soil8(nsoil)%slcons * wgpfrac ** (2.0 * soil(nsoil)%slbs + 3.0) * 1.e3 &   ! kgH2O m-2 s-1
!                                  * sqrt(RAI) / (pi1 * dslz(k))  &  ! m-1
!                                  * 2.0 * cohort_crown_area   ! m2
!
!
!            ! disable hydraulic redistribution
!            ! soil water conductivity is 0. if soil water potential is smaller
!            ! than air dry soil water potential (soilcp)
!            if ((layer_psi(k) <= cpatch%psi_stem(ico)) &
!                .or. (layer_psi(k) < (soilcp_MPa * 102.))) then
!                soil_water_cond = 0.0
!            endif
!            ! The unit of soil_water_cond is kgH2O m-1 s-1
!
!            
!            ! Calculate weighted conductance, weighted psi, and
!            ! water_supply_layer_frac
!            weighted_soil_cond = weighted_soil_cond + soil_water_cond
!            weighted_soil_psi  = weighted_soil_psi + soil_water_cond * layer_psi(k) ! kgH2O s-1
!
!
!            layer_water_supply(k) = &
!                    soil_water_cond * (layer_psi(k) - cpatch%psi_stem(ico))!kgH2O s-1 
!        enddo
!
!        ! Update psi_stem
!        org_psi_stem = cpatch%psi_stem(ico)
!        
!        if (c_stem > 0. .and. weighted_soil_cond > 0.) then
!            ! If the tree stem has some kind of capacity
!
!            ap = - weighted_soil_cond  & !kgH2O m-1 s-1
!                / c_stem  ! kgH2O m-1
!            !the unit of ap is s-1
!
!            bp = (weighted_soil_psi - J_rl) & ! kgH2O s-1
!                / c_stem                    ! kgH2O m-1
!            !the unit of bp is m s-1
!            exp_term = ap * dtlsm
!            if (ap * dtlsm < - 40.) then
!                exp_term = 0.0d0
!            else
!                exp_term = exp(ap * dble(dtlsm))
!            endif
!            cpatch%psi_stem(ico) = sngloff(((ap * org_psi_stem + bp) * exp_term- bp) &
!                                / ap,tiny_offset)
!            J_sr = (cpatch%psi_stem(ico) - org_psi_stem) * c_stem  / dtlsm + &! kgH2O s-1
!                    J_rl
!
!            if (J_sr > 0.) then
!                wflx = (J_sr - J_rl) * dble(dtlsm)
!
!                psi_pot = org_psi_stem +             &
!                          sngloff(wflx / c_stem, tiny_offset)
!
!                if (psi_pot > -1e-6) then
!                    J_sr = 0.
!                
!                    cpatch%psi_stem(ico) = sngloff(org_psi_stem                       &
!                                                   +(J_sr - J_rl) * dtlsm  / c_stem   &
!                                                   ,tiny_offset)
!                endif
!            endif
!
!
!        elseif (c_stem > 0.) then
!                    
!            ! the plants cannot take up water
!            ! change of psi_stem is solely due to J_rl
!            cpatch%psi_stem(ico) = sngloff(org_psi_stem - J_rl * dtlsm/c_stem, &
!                                           tiny_offset)
!            J_sr = 0.0d0
!                
!        else
!            ! This is when the tree stem has no capacity
!            ! This means the all the sapwood and fine roots are dead
!            ! psi_stem cannot change in this case
!            cpatch%psi_stem(ico) = sngloff(org_psi_stem,tiny_offset)
!            J_sr = 0.0d0
!
!        endif
!
!
!
!
!        ! For debugging purpose
!        if(isnan(J_sr)) then
!            write (unit=*,fmt='(80a)')         ('=',k=1,80)
!            write (unit=*,fmt='(a)')           'Soil to root water flow is NaN!!'
!            write (unit=*,fmt='(a,1x,i9)')   ' + ico:                 ',ico
!            write (unit=*,fmt='(a,1x,es12.4)')   ' + STEM CAPACITANCE:    ',c_stem
!            write (unit=*,fmt='(a,1x,es12.4)')   ' + RAI:                 ',RAI
!            write (unit=*,fmt='(a,1x,es12.4)')   ' + NPLANT:              ',cpatch%nplant(ico)
!            write (unit=*,fmt='(a,1x,es12.4)')   ' + SAPFLOW:             ',J_rl
!            write (unit=*,fmt='(a,1x,es12.4)')   ' + SOIL CONDUCTANCE:    ',weighted_soil_cond
!            write (unit=*,fmt='(a,1x,es12.4)')   ' + ORG LEAF PSI:        ',org_psi_leaf
!            write (unit=*,fmt='(a,1x,es12.4)')   ' + ORG STEM PSI:        ',cpatch%psi_stem(ico)
!           
!            call fatal_error('soil to root water flow is wrong','update_plant_hydrodynamics'               &
!                            &,'plant_hydrodynamics.f90')
!        endif
!
!
!        ! for small trees, force psi_leaf to be the same as psi_stem
!        ! adjust fluxes to make sure water budget is conserved
!        if (small_tree_flag) then
!            cpatch%psi_leaf(ico) = cpatch%psi_stem(ico)
!            J_rl = (cpatch%psi_leaf(ico) - org_psi_leaf) * c_leaf  / dtlsm + &! kgH2O s-1
!                    transp
!        endif
!
!        cpatch%water_flux_rl(ico) = sngloff(J_rl,tiny_offset)
!        cpatch%water_flux_sr(ico) = sngloff(J_sr,tiny_offset)
!        !----------------------------------------------------------------------
!        ! Record extracted water from soil layers
!        !-----------------------------------------------------------------------
!  
!        ! update layer_water supply according to the final J_sr value
!  
!        if (sum(layer_water_supply) /= 0.) then
!          layer_water_supply = layer_water_supply / sum(layer_water_supply) * &
!              J_sr
!        else
!          layer_water_supply = 0.d0
!        endif
!
!        ! record layer_water_supply
!        do k = 1,nzg
!            cpatch%water_flux_sr_layer(k,ico) = sngloff(layer_water_supply(k),tiny_offset)
!            cpatch%fmean_water_flux_sr_layer(k,ico) = cpatch%fmean_water_flux_sr_layer(k,ico) &
!                + cpatch%water_flux_sr_layer(k,ico) * dtlsm_o_frqsum
!        enddo
!
!        !----------------------------------------------------------------------!
!        !     Update relative water content and energy                         !
!        !----------------------------------------------------------------------!
!
!        ! Leaf
!        ! Sapflow was assume to be equal to transp during the integration
!        ! Now, correct it with the actual sapflow
!        wflx = (J_rl - transp) * dble(dtlsm)            ! kg H2O
!        qflx = wflx * tl2uint8(dble(cpatch%wood_temp(ico)),1.0d0)            &
!              * dble(cpatch%nplant(ico))                ! J/m2
!       
!        if (cpatch%bleaf(ico) > 0.) then
!            cpatch%leaf_rwc(ico) = cpatch%leaf_rwc(ico)                       &
!                                 + sngloff(wflx                               & ! kg H2O
!                                 / (dble( wat_dry_ratio_grn(ipft))            & ! kg H2O/kg biomass
!                                       * dble(cpatch%bleaf(ico)) * dble(C2B)) &  ! kg biomass
!                                 ,tiny_offset)      
!        else
!            ! no leaves, use psi_leaf to calculate virtual leaf_rwc
!            cpatch%leaf_rwc(ico) = cpatch%psi_leaf(ico) * Cap_leaf(ipft)      &
!                                 / wat_dry_ratio_grn(ipft)                    &
!                                 + 1.0
!        endif
!
!        cpatch%leaf_energy(ico) = cpatch%leaf_energy(ico)                 &
!                                + sngloff(qflx,tiny_offset)
!
!
!        ! Wood and soil
!        wflx = -J_rl
!        qflx = wflx * tl2uint8(dble(cpatch%wood_temp(ico)),1.0d0)
!        ! Loop over soil layer to account for soil water extraction
!        do k = 1,nzg
!            soil_wflx = layer_water_supply(k)                              ! kg/s
!            soil_qflx = layer_water_supply(k)                              &
!                      * tl2uint8(dble(csite%soil_tempk(k,ipa)),1.0d0)
!            !wflx      = wflx + soil_wflx
!            qflx      = qflx + soil_qflx
!
!            ! update soil water and energy
!            csite%soil_water(k,ipa) = csite%soil_water(k,ipa)               &
!                    - sngloff(soil_wflx                                     & ! kg/s
!                    * dble(cpatch%nplant(ico) * wdnsi * dslzi(k) * dtlsm)   &! s/kg
!                    , tiny_offset)
!            csite%soil_energy(k,ipa) = csite%soil_energy(k,ipa)             &
!                    - sngloff(soil_qflx                                     & ! J/s
!                    * dble(cpatch%nplant(ico) * dslzi(k) * dtlsm)           & ! s/m3
!                    , tiny_offset)
!        enddo
!        wflx = (wflx + J_sr) * dble(dtlsm)
!        qflx = qflx * dble(cpatch%nplant(ico) * dtlsm)
!
!        cpatch%wood_rwc(ico) = cpatch%wood_rwc(ico)                         &
!                             + sngloff(wflx                                 & ! kg H2O
!                             / (dble( wat_dry_ratio_ngrn(ipft))             & ! kg H2O/kg biomass
!                                   * ( dble(cpatch%broot(ico))              & ! kg
!                                     + dble(cpatch%bdead(ico))              & ! kg
!                                     * dble(xylem_fraction(ipft)))          &
!                                     * dble(C2B))                           &
!                             , tiny_offset)
!
!        cpatch%wood_energy(ico) = cpatch%wood_energy(ico)                   &
!                                + sngloff(qflx,tiny_offset)
!
! 
!
!        !----------------------------------------------------------!
!        ! Check the relative error btween psi and rwc              !
!        ! Make sure they are consistent                            !
!        !----------------------------------------------------------!
!        ! leaf
!        rwc_err = abs(cpatch%leaf_rwc(ico)                              &
!                - (1.0 + cpatch%psi_leaf(ico)  * Cap_leaf(ipft)       &
!                  / wat_dry_ratio_grn(ipft)))
!        if (rwc_err > 1.e-2 .or. &
!            cpatch%leaf_rwc(ico) > 1.01 .or. &
!            cpatch%leaf_rwc(ico) < -0.01) then
!            print*, 'ico',ico
!            print*, 'pft',ipft
!            print*, 'leaf_rwc',cpatch%leaf_rwc(ico)
!            print*, 'leaf_psi',cpatch%psi_leaf(ico)
!            print*, 'org_leaf_psi',org_psi_leaf
!            print*, 'stem_psi',cpatch%psi_stem(ico)
!            print*, 'org_stem_psi',org_psi_stem
!            print*, 'rwc_err',rwc_err
!            print*, 'transp',transp
!            print*, 'sapflow',J_rl
!            print*, 'water uptake', sum(layer_water_supply)
!
!
!            call fatal_error('Error of Leaf RWC is too large','update_plant_hydrodynamics'               &
!                            &,'plant_hydrodynamics.f90')
!        else
!            ! update rwc using leaf_psi
!!            cpatch%psi_leaf(ico) = (cpatch%leaf_rwc(ico) - 1.0)         &
!!                                 * wat_dry_ratio_grn(ipft)              &
!!                                 / Cap_leaf(ipft)
!            cpatch%leaf_rwc(ico) = (1. + cpatch%psi_leaf(ico) * Cap_leaf(ipft) &
!                                   / wat_dry_ratio_grn(ipft))
!        endif
!
!        ! wood
!        rwc_err = abs(cpatch%wood_rwc(ico)                              &
!                - (1.0 + cpatch%psi_stem(ico)  * Cap_stem(ipft)       &
!                  / wat_dry_ratio_ngrn(ipft)))
!        if (rwc_err > 1.e-2 .or. &
!            cpatch%wood_rwc(ico) > 1.01 .or. &
!            cpatch%wood_rwc(ico) < -0.01) then
!            print*, 'ico',ico
!            print*, 'pft',ipft
!            print*, 'wood_rwc',cpatch%wood_rwc(ico)
!            print*, 'leaf_psi',cpatch%psi_leaf(ico)
!            print*, 'org_leaf_psi',org_psi_leaf
!            print*, 'stem_psi',cpatch%psi_stem(ico)
!            print*, 'org_stem_psi',org_psi_stem
!            print*, 'rwc_err',rwc_err
!            print*, 'transp',transp
!            print*, 'sapflow',J_rl
!            print*, 'water uptake', sum(layer_water_supply)
!            call fatal_error('Error of Wood RWC is too large','update_plant_hydrodynamics'               &
!                            &,'plant_hydrodynamics.f90')
!        else
!            ! update psi_stem using rwc
!  !          cpatch%psi_stem(ico) = (cpatch%wood_rwc(ico) - 1.0)         &
!  !                               * wat_dry_ratio_ngrn(ipft)              &
!  !                               / Cap_stem(ipft)
!            ! update rwc using psi_stem
!            cpatch%wood_rwc(ico) = (1. + cpatch%psi_stem(ico) * Cap_stem(ipft) &
!                                    / wat_dry_ratio_ngrn(ipft))
!        endif
!
!
!        !now update plant interal water
!        call update_veg_water_int(cpatch,ico)
!
!        !----------------------------------------------------------!
!        ! Heat capacity changes when rwc changes
!        ! Update Hcap and update temperature
!        !----------------------------------------------------------!
!        call calc_veg_hcap(cpatch%bleaf(ico),cpatch%broot(ico)                              &
!                            ,cpatch%bdead(ico),cpatch%bsapwooda(ico)                        &
!                            ,cpatch%nplant(ico),cpatch%pft(ico)                             &
!                            ,cpatch%leaf_rwc(ico),cpatch%wood_rwc(ico)                      &
!                            ,cpatch%leaf_hcap(ico),cpatch%wood_hcap(ico) )
!        ! energy is correct, update temperature
!        call uextcm2tl(cpatch%leaf_energy(ico),cpatch%leaf_water(ico)            &
!                       ,cpatch%leaf_hcap(ico),cpatch%leaf_temp(ico)               &
!                       ,cpatch%leaf_fliq(ico))
!        call uextcm2tl(cpatch%wood_energy(ico),cpatch%wood_water(ico)            &
!                       ,cpatch%wood_hcap(ico),cpatch%wood_temp(ico)               &
!                       ,cpatch%wood_fliq(ico))
!
!        
!        !----------------------------------------------------------------------
!        ! integrate fast-analysis and dmax/dmin variables
!        !----------------------------------------------------------------------
!        cpatch%fmean_psi_leaf(ico)      = cpatch%fmean_psi_leaf(ico)        &
!                                        + cpatch%psi_leaf(ico) * dtlsm_o_frqsum
!        cpatch%fmean_psi_stem(ico)      = cpatch%fmean_psi_stem(ico)        &
!                                        + cpatch%psi_stem(ico) * dtlsm_o_frqsum
!        cpatch%fmean_water_flux_rl(ico) = cpatch%fmean_water_flux_rl(ico)   &
!                                        + cpatch%water_flux_rl(ico) * dtlsm_o_frqsum
!        cpatch%fmean_water_flux_sr(ico) = cpatch%fmean_water_flux_sr(ico)   &
!                                        + cpatch%water_flux_sr(ico) * dtlsm_o_frqsum
!        cpatch%fmean_leaf_rwc(ico)      = cpatch%fmean_leaf_rwc(ico)        &
!                                        + cpatch%leaf_rwc(ico) * dtlsm_o_frqsum
!        cpatch%fmean_wood_rwc(ico)      = cpatch%fmean_wood_rwc(ico)        &
!                                        + cpatch%wood_rwc(ico) * dtlsm_o_frqsum
!        cpatch%fmean_leaf_water_int(ico)= cpatch%fmean_leaf_water_int(ico)  &
!                                        + cpatch%leaf_water_int(ico) * dtlsm_o_frqsum
!        cpatch%fmean_wood_water_int(ico)= cpatch%fmean_wood_water_int(ico)  &
!                                        + cpatch%wood_water_int(ico) * dtlsm_o_frqsum
!
!        if (cpatch%dmax_psi_leaf(ico) == 0. .or. update_psi) then
!            cpatch%dmax_psi_leaf(ico) = cpatch%psi_leaf(ico)
!        else
!            cpatch%dmax_psi_leaf(ico) = max(cpatch%dmax_psi_leaf(ico),      &
!                                            cpatch%psi_leaf(ico))
!        endif
!
!        if (cpatch%dmin_psi_leaf(ico) == 0. .or. update_psi) then
!            cpatch%dmin_psi_leaf(ico) = cpatch%psi_leaf(ico)
!        else
!            cpatch%dmin_psi_leaf(ico) = min(cpatch%dmin_psi_leaf(ico),      &
!                                            cpatch%psi_leaf(ico))
!        endif
!
!        if (cpatch%dmax_psi_stem(ico) == 0. .or. update_psi) then
!            cpatch%dmax_psi_stem(ico) = cpatch%psi_stem(ico)
!        else
!            cpatch%dmax_psi_stem(ico) = max(cpatch%dmax_psi_stem(ico),      &
!                                            cpatch%psi_stem(ico))
!        endif
!
!        if (cpatch%dmin_psi_stem(ico) == 0. .or. update_psi) then
!            cpatch%dmin_psi_stem(ico) = cpatch%psi_stem(ico)
!        else
!            cpatch%dmin_psi_stem(ico) = min(cpatch%dmin_psi_stem(ico),      &
!                                            cpatch%psi_stem(ico))
!        endif
!
!
!  
!        if ((quality_check .and. &
!            ((check_ico == 0) .or. (ico == check_ico))) &
!            .or. cpatch%psi_leaf(ico) > 0. .or. cpatch%psi_stem(ico) > 0.) then
!           write (unit=*,fmt='(80a)')         ('=',k=1,80)
!           write (unit=*,fmt='(a)')           'Plant Hydrodynamics Quality Check:'
!           write (unit=*,fmt='(a,1x,i9)')   ' + HOUR:                ',current_time%hour
!           write (unit=*,fmt='(a,1x,i9)')   ' + ICO:                 ',ico
!           write (unit=*,fmt='(a,1x,i9)')   ' + PFT:                 ',cpatch%pft(ico)
!           write (unit=*,fmt='(a,1x,f9.4)')   ' + DBH:                 ',cpatch%dbh(ico)
!           write (unit=*,fmt='(a,1x,f9.4)')   ' + ELONGF:              ',cpatch%elongf(ico)
!           write (unit=*,fmt='(a,1x,f9.4)')   ' + BLEAF:               ',cpatch%bleaf(ico) 
!           write (unit=*,fmt='(a,1x,i9)')   ' + KRDEPTH:             ',cpatch%krdepth(ico)
!           write (unit=*,fmt='(a,1x,es12.4)')   ' + PSI_STEM:            ',cpatch%psi_stem(ico)
!           write (unit=*,fmt='(a,1x,es12.4)')   ' + PSI_LEAF:            ',cpatch%psi_leaf(ico)
!           write (unit=*,fmt='(a,1x,es12.4)')   ' + ORG_PSI_STEM:        ',org_psi_stem
!           write (unit=*,fmt='(a,1x,es12.4)')   ' + ORG_PSI_LEAF         ',org_psi_leaf
!           write (unit=*,fmt='(a,1x,es12.4)')   ' + LEAF_RWC:            ',cpatch%leaf_rwc(ico)
!           write (unit=*,fmt='(a,1x,es12.4)')   ' + WOOD_RWC:            ',cpatch%wood_rwc(ico)
!           write (unit=*,fmt='(a,1x,es12.4)')   ' + LEAF_TEMP:           ',cpatch%leaf_temp(ico)
!           write (unit=*,fmt='(a,1x,es12.4)')   ' + WOOD_TEMP:           ',cpatch%wood_temp(ico)
!           write (unit=*,fmt='(a,1x,es12.4)')   ' + TRANSP:              ',transp
!           write (unit=*,fmt='(a,1x,es12.4)')   ' + WATER_FLUX_RL:       ',cpatch%water_flux_rl(ico)
!           write (unit=*,fmt='(a,1x,es12.4)')   ' + WATER_FLUX_SR:       ',cpatch%water_flux_sr(ico)
!           write (unit=*,fmt='(a,1x,3f9.4)')  ' + SOIL_WATER (top 3 layer):       ',csite%soil_water(nzg-3:nzg,ipa)
!           write (unit=*,fmt='(a,1x,3f9.4)')  ' + SOIL_TEMPK (top 3 layer):       ',csite%soil_tempk(nzg-3:nzg,ipa)
!
!           if (cpatch%psi_leaf(ico) > 0. .or. cpatch%psi_stem(ico) > 0.) then
!                call fatal_error('Water potential is positive','update_plant_hydrodynamics' &
!                                ,'plant_hydrodynamics.f90')
!           endif
!        endif
!
!    enddo cohortloop
!  
!    ! update soil temperature, it should not change
!    do k = 1,nzg
!      nsoil = ntext_soil(k)
!      call uextcm2tl(csite%soil_energy(k,ipa), csite%soil_water(k,ipa)*wdns                &
!                    ,soil(nsoil)%slcpd, csite%soil_tempk(k,ipa), csite%soil_fracliq(k,ipa))
!    enddo
!
!
!    
!   end subroutine update_plant_hydrodynamics

   !=======================================================================================!
   !=======================================================================================!
   !   Util functions/subroutines for plant hydrodynamic calculations                      !
   !=======================================================================================!
   !=======================================================================================!


   !=======================================================================================!
   !=======================================================================================!
   !  SUBROUTINE: PSI2RWC           
   !> \breif Convert water potential of leaf and wood to relative water content
   !> \details Here we assume a constant hydraulic capacitance for both leaf and
   !> wood. From the definition of hydraulic capacitance we have \n
   !>       hydro_cap = delta_water_content / delta_psi \n
   !> Since psi = 0. when water_content is at saturation, we have \n
   !>       hydro_cap = (1. - rwc) * water_content_at_saturation / (0. - psi) \n
   !> Reorganize the equation above, we can get \n
   !>       rwc = 1. + psi * hydro_cap / water_content_at_saturation
   !=======================================================================================!
   subroutine psi2rwc(leaf_psi,wood_psi,ipft,leaf_rwc,wood_rwc)
      use pft_coms          ,   only : leaf_hydro_cap       & ! intent(in)
                                     , wood_hydro_cap       & ! intent(in)
                                     , leaf_water_sat       & ! intent(in)
                                     , wood_water_sat       ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real      , intent(in)    ::  leaf_psi    ! Water potential of leaves            [  m]
      real      , intent(in)    ::  wood_psi    ! Water potential of wood              [  m]
      integer   , intent(in)    ::  ipft        ! plant functional type                [  -]
      real      , intent(out)   ::  leaf_rwc    ! Relative water content of leaves     [0-1]
      real      , intent(out)   ::  wood_rwc    ! Relative water content of wood       [0-1]

      ! first caculate for leaf
      leaf_rwc  =   1.  + leaf_psi                  & ! [m]
                        * leaf_hydro_cap(ipft)      & ! [kg H2O/kg biomass/m]
                        / leaf_water_sat(ipft)      ! ! [kg H2O/kg biomass]

      ! same for wood
      wood_rwc  =   1.  + wood_psi                  & ! [m]
                        * wood_hydro_cap(ipft)      & ! [kg H2O/kg biomass/m]
                        / wood_water_sat(ipft)      ! ! [kg H2O/kg biomass]

   end subroutine psi2rwc

   !=======================================================================================!
   !=======================================================================================!
   !  SUBROUTINE: RWC2PSI
   !> \brief Convert relative water content to water potential
   !> \details The inverse of psi2rwc
   !=======================================================================================!
   subroutine rwc2psi(leaf_rwc,wood_rwc,ipft,leaf_psi,wood_psi)
      use pft_coms          ,   only : leaf_hydro_cap       & ! intent(in)
                                     , wood_hydro_cap       & ! intent(in)
                                     , leaf_water_sat       & ! intent(in)
                                     , wood_water_sat       ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real      , intent(in)    ::  leaf_rwc    ! Relative water content of leaves     [0-1]
      real      , intent(in)    ::  wood_rwc    ! Relative water content of wood       [0-1]
      integer   , intent(in)    ::  ipft        ! plant functional type                [  -]
      real      , intent(out)   ::  leaf_psi    ! Water potential of leaves            [  m]
      real      , intent(out)   ::  wood_psi    ! Water potential of wood              [  m]

      ! first caculate for leaf
      leaf_psi  =   (leaf_rwc - 1.) * leaf_water_sat(ipft) / leaf_hydro_cap(ipft)
      ! same for wood
      wood_psi  =   (wood_rwc - 1.) * wood_water_sat(ipft) / wood_hydro_cap(ipft)

   end subroutine rwc2psi


   !=======================================================================================!
   !=======================================================================================!
   !  SUBROUTINE: RWC2TW            
   !> \breif Convert relative water content to total water for both leaf and wood
   !> \details total water content = \n
   !>  relative water content * water_content_at_saturation * biomass
   !> \warning In the hydro version, bsapwood is set as 0 and bdead is assumedto contain
   !> both sapwood and heart wood. Root is counted as wood.
   !=======================================================================================!
   subroutine rwc2tw(leaf_rwc,wood_rwc,bleaf,bdead,broot,ipft                    &
                    ,leaf_water_int,wood_water_int)
      use pft_coms      ,   only : xylem_fraction       & ! intent(in)
                                 , wood_rwc_min         & ! intent(in)
                                 , C2B                  & ! intent(in)
                                 , leaf_water_sat       & ! intent(in)
                                 , wood_water_sat       ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real      , intent(in)    ::  leaf_rwc    ! Relative water content of leaves     [0-1]
      real      , intent(in)    ::  wood_rwc    ! Relative water content of wood       [0-1]
      real      , intent(in)    ::  bleaf       ! Biomass of leaf                      [kgC]
      real      , intent(in)    ::  bdead       ! Biomass of wood                      [kgC]
      real      , intent(in)    ::  broot       ! Biomass of fine root                 [kgC]
      integer   , intent(in)    ::  ipft        ! Plant functional type                [  -]
      real      , intent(out)   ::  leaf_water_int ! Total internal water of leaf      [ kg]
      real      , intent(out)   ::  wood_water_int ! Total internal water of wood      [ kg]

      ! leaf
      leaf_water_int    =   leaf_rwc * leaf_water_sat(ipft) * bleaf * C2B 

      ! wood
      wood_water_int    =   wood_rwc * wood_water_sat(ipft)                 &
                        *   (broot + bdead * xylem_fraction(ipft)) * C2B    & ! live biomass
                        +   wood_rwc_min(ipft) * wood_water_sat(ipft)       &
                        *   bdead * (1. - xylem_fraction(ipft)) * C2B         ! dead biomass

   end subroutine rwc2tw

   !=======================================================================================!
   !=======================================================================================!
   !  SUBROUTINE: TW2RWC            
   !> \breif Convert total water to relative water content for both leaf and wood
   !> \details the inverse of rwc2tw \n
   !=======================================================================================!
   subroutine tw2rwc(leaf_water_int,wood_water_int,bleaf,bdead,broot,ipft        &
                    ,leaf_rwc,wood_rwc)
      use pft_coms      ,   only : xylem_fraction       & ! intent(in)
                                 , wood_rwc_min         & ! intent(in)
                                 , C2B                  & ! intent(in)
                                 , leaf_water_sat       & ! intent(in)
                                 , wood_water_sat       ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real      , intent(in)    ::  leaf_water_int ! Total internal water of leaf      [ kg]
      real      , intent(in)    ::  wood_water_int ! Total internal water of wood      [ kg]
      real      , intent(in)    ::  bleaf       ! Biomass of leaf                      [kgC]
      real      , intent(in)    ::  bdead       ! Biomass of wood                      [kgC]
      real      , intent(in)    ::  broot       ! Biomass of fine root                 [kgC]
      integer   , intent(in)    ::  ipft        ! Plant functional type                [  -]
      real      , intent(out)   ::  leaf_rwc    ! Relative water content of leaves     [0-1]
      real      , intent(out)   ::  wood_rwc    ! Relative water content of wood       [0-1]

      ! leaf
      leaf_rwc          =   leaf_water_int / (leaf_water_sat(ipft) * bleaf * C2B)
      
      ! wood
      wood_rwc          =   (wood_water_int                                 &
                        -   wood_rwc_min(ipft) * wood_water_sat(ipft)       &
                        *   bdead * (1. - xylem_fraction(ipft)) * C2B)      &
                        /   (wood_water_sat(ipft)                           &
                        *   (broot + bdead * xylem_fraction(ipft)) * C2B)
   end subroutine tw2rwc

   !=======================================================================================!
   !=======================================================================================!
   !  SUBROUTINE: PSI2TW            
   !> \breif Convert water potential to total water for both leaf and wood
   !=======================================================================================!
   subroutine psi2tw(leaf_psi,wood_psi,bleaf,bdead,broot,ipft                           &
                    ,leaf_water_int,wood_water_int)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real      , intent(in)    ::  leaf_psi    ! Water potential of leaves            [  m]
      real      , intent(in)    ::  wood_psi    ! Water potential of wood              [  m]
      real      , intent(in)    ::  bleaf       ! Biomass of leaf                      [kgC]
      real      , intent(in)    ::  bdead       ! Biomass of wood                      [kgC]
      real      , intent(in)    ::  broot       ! Biomass of fine root                 [kgC]
      integer   , intent(in)    ::  ipft        ! Plant functional type                [  -]
      real      , intent(out)   ::  leaf_water_int ! Total internal water of leaf      [ kg]
      real      , intent(out)   ::  wood_water_int ! Total internal water of wood      [ kg]
      !----- Local Variables --------------------------------------------------------------!
      real                      ::  leaf_rwc    ! Relative water content of leaf       [  -]
      real                      ::  wood_rwc    ! Relative water content of wood       [  -]

      ! first convert to rwc
      call psi2rwc(leaf_psi,wood_psi,ipft,leaf_rwc,wood_rwc)
      ! second convert to tw
      call rwc2tw(leaf_rwc,wood_rwc,bleaf,bdead,broot,ipft,leaf_water_int,wood_water_int)

   end subroutine psi2tw


   !=======================================================================================!
   !=======================================================================================!
   !  SUBROUTINE: TW2PSI            
   !> \breif Convert total water to water potential for both leaf and wood
   !> \details the inverse of psi2tw \n
   !=======================================================================================!
   subroutine tw2psi(leaf_water_int,wood_water_int,bleaf,bdead,broot,ipft        &
                    ,leaf_psi,wood_psi)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real      , intent(in)    ::  leaf_water_int ! Total internal water of leaf      [ kg]
      real      , intent(in)    ::  wood_water_int ! Total internal water of wood      [ kg]
      real      , intent(in)    ::  bleaf       ! Biomass of leaf                      [kgC]
      real      , intent(in)    ::  bdead       ! Biomass of wood                      [kgC]
      real      , intent(in)    ::  broot       ! Biomass of fine root                 [kgC]
      integer   , intent(in)    ::  ipft        ! Plant functional type                [  -]
      real      , intent(out)   ::  leaf_psi    ! Water potential of leaves            [  m]
      real      , intent(out)   ::  wood_psi    ! Water potential of wood              [  m]
      !----- Local Variables --------------------------------------------------------------!
      real                      ::  leaf_rwc    ! Relative water content of leaf       [  -]
      real                      ::  wood_rwc    ! Relative water content of wood       [  -]

      ! first convert to rwc
      call tw2rwc(leaf_water_int,wood_water_int,bleaf,bdead,broot,ipft,leaf_rwc,wood_rwc)
      ! second convert to psi
      call rwc2psi(leaf_rwc,wood_rwc,ipft,leaf_psi,wood_psi)

   end subroutine tw2psi

end module plant_hydro

!==========================================================================================!
!==========================================================================================!
