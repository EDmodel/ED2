!==========================================================================================!
!==========================================================================================!
! MODULE: PLANT_HYDRO
!> \brief Calculations of plant hydrodynamics at ED2 timestep, including various
!> utils for plant hydraulics
!> \details Util functions to perform conversions between psi, rwc and water_int
!> \author Xiangtao Xu, 26 Jan. 2018
!==========================================================================================!
!==========================================================================================!
module plant_hydro

   contains
    
   !=======================================================================================!
   !=======================================================================================!
   ! SUBROUTINE PLANT_HYDRO_DRIVER
   !> \brief   Main driver to calculate plant hydrodynamics within a site.
   !> \details This subroutine works at DTLSM scale, similar to photosyn_driv.
   !> Keep in mind that this subroutine uses average transpiration from last timestep.
   !> Therefore, we should also use water potential at the start of the last
   !> timestep. \n
   !>   Alternatively, one can directly call calc_plant_water_flux
   !> subroutine within the rk4_derivs.f90, which will give an estimate of water
   !> flux within each integration timestep using the transpiration/water
   !> potential at the start of the integration timestep. However, this can
   !> incur extra computational cost.
   !>
   !> \author Xiangtao Xu, 30 Jan. 2018
   !---------------------------------------------------------------------------------------!
   subroutine plant_hydro_driver(csite,ipa,ntext_soil)
      use ed_state_vars        , only : sitetype               & ! structure
                                      , patchtype              ! ! structure
      use ed_misc_coms         , only : dtlsm                  & ! intent(in)
                                      , frqsum                 & ! intent(in)
                                      , current_time           ! ! intent(in)
      use soil_coms            , only : soil                   ! ! intent(in)
      use grid_coms            , only : nzg                    ! ! intent(in)
      use consts_coms          , only : pio4                   ! ! intent(in)
      use allometry            , only : dbh2sf                 ! ! function
      use physiology_coms      , only : plant_hydro_scheme     ! ! intent(in)
      use pft_coms             , only : C2B                    & ! intent(in)
                                      , leaf_water_cap         ! ! intent(in)

      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(sitetype)            , target      :: csite
      integer                   , intent(in)  :: ipa
      integer    ,dimension(nzg), intent(in)  :: ntext_soil
      !----- Local Vars  --------------------------------------------------------------------!
      type(patchtype)           , pointer     :: cpatch     ! patch strcture
      real                                    :: wgpfrac    ! relative soil moisture
      integer                                 :: nsoil      ! soil type for soil
      integer                                 :: k          ! iterator for soil layer
      integer                                 :: ico        ! iterator for cohort
      real       ,dimension(nzg)              :: soil_psi   ! soil water potential    [m]
      real       ,dimension(nzg)              :: soil_cond  ! soil water conductance  [kg/m2/s]
      real                                    :: sap_frac   ! sapwood fraction      
      real                                    :: sap_area   ! sapwood area            [m2]
      real                                    :: bsap       ! sapwood biomass         [kgC]
      real                                    :: crown_area ! crown area              [m2]
      real                                    :: transp     ! transpiration rate      [kg/s]
      real                                    :: c_leaf     ! leaf capacitance        [kg/m]
      !----- Variables for debugging purposes ---------------------------------------------!
      integer,parameter                       :: dco        = 0 ! the cohort to debug
      logical                                 :: error_flag
      logical,parameter                       :: debug_flag = .false.
      character(len=13) , parameter           :: efmt       = '(a,1x,es12.5)'
      character(len=9)  , parameter           :: ifmt       = '(a,1x,i5)'
      character(len=9)  , parameter           :: lfmt       = '(a,1x,l1)'
      !----- Locally saved variables. -----------------------------------------------------!
      real        , save                      :: dtlsm_o_frqsum
      logical     , save                      :: first_time = .true.
      !----- Assign the constant scaling factor. ------------------------------------------!
      if (first_time) then
         first_time     = .false.
         dtlsm_o_frqsum = dtlsm / frqsum
      end if
      !------------------------------------------------------------------------------------!


      !-- Point to the cohort structures --------------------------------------------------!
      cpatch => csite%patch(ipa)
      !------------------------------------------------------------------------------------!

      select case (plant_hydro_scheme)
      case (0)
        ! Compatible with original ED-2.2, do not track plant hydraulics
        do ico = 1, cpatch%ncohorts
            cpatch%wflux_wl        (ico)    = 0.
            cpatch%wflux_gw        (ico)    = 0.
            cpatch%wflux_gw_layer(:,ico)    = 0.
            
            cpatch%leaf_rwc        (ico)    = 0.
            cpatch%wood_rwc        (ico)    = 0.
            cpatch%leaf_psi        (ico)    = 0.
            cpatch%wood_psi        (ico)    = 0.
        end do
      case default
        ! track plant hydraulics

        !--------------------------------------------------------------------------
        ! Calculate water potential and conductance in each soil layer in preparation for
        ! later calculations
        !--------------------------------------------------------------------------
        do k = 1,nzg
          nsoil = ntext_soil(k)
          
          !get relative soil moisture

          wgpfrac = max(soil(nsoil)%soilcp,min(1.0,                                &
                    csite%soil_water(k,ipa) * csite%soil_fracliq(k,ipa)           &
                    / soil(nsoil)%slmsts))

          ! Clapp & Horn curves
          soil_psi(k)  = soil(nsoil)%slpots / wgpfrac ** soil(nsoil)%slbs ! m


          !In the model, soil can't get drier than Dry soil capacity (-3.1MPa, can
          !be changed in ed_params.f90)
          !So, we turn off soil-wood water conductance when that's the case
          if (csite%soil_water(k,ipa) <= soil(nsoil)%soilcp) then
              soil_cond(k) = 0.
          else
              soil_cond(k) = soil(nsoil)%slcons                                    &
                           * wgpfrac ** (2.0 * soil(nsoil)%slbs + 3.0) * 1.e3
              ! kgH2O m-2 s-1
          end if

        end do

        ! Loop over cohorts, calculate plant hydraulic fluxes
        cohortloop: do ico = 1, cpatch%ncohorts

            ! track the plant hydraulics only when leaf is resolvable
            ! also need to track plant hydro when leaf is not resolvable
            ! When plants shed all the leaves during the dry season, leaf_hcap
            ! becomes 0 and the cohort becomes non-resolvable. But if we do not
            ! track water potential in this case, we can never allow soil water
            ! to refill wood

            if (cpatch%leaf_resolvable(ico)) then
                !---- prepare input for plant water flux calculations
                sap_frac    = dbh2sf(cpatch%dbh(ico),cpatch%pft(ico))         ! m2
                sap_area    = sap_frac * pio4 * (cpatch%dbh(ico) / 100.) ** 2 ! m2
                bsap        = (cpatch%bdeada(ico) + cpatch%bdeadb(ico)) * sap_frac
                crown_area  = cpatch%crown_area(ico) / cpatch%nplant(ico)  ! m2
                transp      = ( cpatch%fs_open(ico) * cpatch%psi_open(ico)               &
                              + (1. - cpatch%fs_open(ico)) * cpatch%psi_closed(ico)      &
                              ) * cpatch%lai(ico) / cpatch%nplant(ico)     ! kg / s


                ! Please notice that the current leaf_water_int has included the
                ! transpirational lost of this time step but not the sapflow. So
                ! leaf psi can be very low.
                ! To get meaningful leaf_psi, I decoupled leaf_psi and
                ! leaf_water_int by adding the transpirational cost back to
                ! leaf_water_int
                ! In this case, leaf_psi represents the water potential at the
                ! START of the timestep.
                call rwc2psi(cpatch%leaf_rwc(ico),cpatch%wood_rwc(ico),cpatch%pft(ico)     &
                            ,cpatch%leaf_psi(ico),cpatch%wood_psi(ico))

                c_leaf = leaf_water_cap(cpatch%pft(ico)) * C2B * cpatch%bleaf(ico)
                if (c_leaf > 0.) then
            
                    cpatch%leaf_psi(ico) = cpatch%leaf_psi(ico)                  &
                                         + transp * dtlsm                        &  ! kg/H2O
                                         / c_leaf                                   ! kgH2O/m
                else
                    ! no leaves, set leaf_psi the same as wood_psi - hite
                    cpatch%leaf_psi(ico) = cpatch%wood_psi(ico) - cpatch%hite(ico)
!                    ! need to reset rwc
!                    call psi2rwc(cpatch%leaf_psi(ico),cpatch%wood_psi(ico),cpatch%pft(ico) &
!                                 cpatch%leaf_rwc(ico),cpatch%wood_rwc(ico))
                end if
                
                !--------------------------------------------------------------------------
                ! Handling Potential Errors and Help Debugging
                !--------------------------------------------------------------------------
                error_flag = (isnan(cpatch%leaf_psi(ico)) )           & ! NaN values
                         .or.(cpatch%leaf_psi(ico) > 0.)                ! psi is positive
  
                ! I copy the error printing from rk4_misc.f90
                if((debug_flag .and. (dco == 0 .or. ico == dco)) .or. error_flag) then
                   write (unit=*,fmt='(a)') ' '
                   write (unit=*,fmt='(92a)') ('=',k=1,92)
                   write (unit=*,fmt='(92a)') ('=',k=1,92)
                   write (unit=*,fmt='(a)'  ) ' Input leaf_psi is too high for plant hydodynamics'
                   write (unit=*,fmt='(92a)') ('-',k=1,92)
                   write (unit=*,fmt='(a,i4.4,2(1x,i2.2),1x,f6.0)') ' TIME           : '   &
                                                   ,current_time%year,current_time%month   &
                                                   ,current_time%date,current_time%time
                   write (unit=*,fmt='(a)'  ) ' '
                   write (unit=*,fmt=ifmt   ) ' + IPA              =',ipa
                   write (unit=*,fmt=ifmt   ) ' + ICO              =',ico
                   write (unit=*,fmt=ifmt   ) ' + PFT              =',cpatch%pft(ico)
                   write (unit=*,fmt=ifmt   ) ' + KRDEPTH          =',cpatch%krdepth(ico)
                   write (unit=*,fmt=efmt   ) ' + HEIGHT           =',cpatch%hite(ico)
  
                   write (unit=*,fmt='(a)'  ) ' '
                   write (unit=*,fmt=efmt   ) ' + BLEAF            =',cpatch%bleaf(ico)
                   write (unit=*,fmt=efmt   ) ' + LAI              =',cpatch%lai(ico) 
                   write (unit=*,fmt=efmt   ) ' + NPLANT           =',cpatch%nplant(ico) 
                   write (unit=*,fmt=efmt   ) ' + BSAPWOOD         =',bsap 
                   write (unit=*,fmt=efmt   ) ' + BROOT            =',cpatch%broot(ico)
                   write (unit=*,fmt=efmt   ) ' + SAPWOOD_AREA     =',sap_area
                   write (unit=*,fmt=efmt   ) ' + CROWN_AERA       =',crown_area
  
                   write (unit=*,fmt='(a)'  ) ' '
                   write (unit=*,fmt=efmt   ) ' + TRANSP           =',transp
                   write (unit=*,fmt=efmt   ) ' + C_LEAF           =',c_leaf
                   write (unit=*,fmt=efmt   ) ' + PSI_OPEN         =',cpatch%psi_open(ico)
                   write (unit=*,fmt=efmt   ) ' + PSI_CLOSED       =',cpatch%psi_closed(ico)
                   write (unit=*,fmt=efmt   ) ' + FS_OPEN          =',cpatch%fs_open(ico)
                   write (unit=*,fmt='(a)'  ) ' '
                   write (unit=*,fmt=efmt   ) ' + LEAF_PSI         =',cpatch%leaf_psi(ico)
                   write (unit=*,fmt=efmt   ) ' + LEAF_RWC         =',cpatch%leaf_rwc(ico)
                   write (unit=*,fmt=efmt   ) ' + LEAF_WATER_INT   =',cpatch%leaf_water_int(ico)
                   write (unit=*,fmt=efmt   ) ' + LEAF_WATER_IM2   =',cpatch%leaf_water_im2(ico)
                   write (unit=*,fmt='(a)'  ) ' '
                   write (unit=*,fmt=efmt   ) ' + WOOD_PSI         =',cpatch%wood_psi(ico)
                   write (unit=*,fmt=efmt   ) ' + WOOD_RWC         =',cpatch%wood_rwc(ico)
                   write (unit=*,fmt=efmt   ) ' + WOOD_WATER_INT   =',cpatch%wood_water_int(ico)
                   write (unit=*,fmt=efmt   ) ' + WOOD_WATER_IM2   =',cpatch%wood_water_im2(ico)
                   write (unit=*,fmt='(a)'  ) ' '
                   write (unit=*,fmt=efmt   ) ' + WFLUX_GW (LAST)  =',cpatch%wflux_gw(ico) 
                   write (unit=*,fmt=efmt   ) ' + WFLUX_WL (LAST)  =',cpatch%wflux_wl(ico)
                   write (unit=*,fmt='(a)'  ) ' + WFLUX_GW_LAYER   ='            
                   do k = 1, nzg
                       write (unit=*,fmt='(i5,1x,es12.5)') k, cpatch%wflux_gw_layer(k,ico)
                   end do
  
                   write (unit=*,fmt='(a)'  ) ' '
                   write (unit=*,fmt='(a)'  ) ' + SOIL_PSI         ='            
                   do k = 1, nzg
                       write (unit=*,fmt='(i5,1x,es12.5)') k, soil_psi(k)                 
                   end do
                   write (unit=*,fmt='(92a)') ('=',k=1,92)
                   write (unit=*,fmt='(92a)') ('=',k=1,92)
                   write (unit=*,fmt='(a)') ' '
  
                   if (error_flag) then 
                       call fatal_error('Plant Hydrodynamics is wrong',&
                                        'plant_hydro_driver'       ,&
                                        'plant_hydro.f90')
                   end if
                end if

                ! note here, transp is from last timestep's psi_open and psi_closed
                call calc_plant_water_flux(                           &
                        dtlsm                                         &!input
                       ,sap_area,crown_area,cpatch%krdepth(ico)       &!input
                       ,cpatch%bleaf(ico),bsap,cpatch%broot(ico)      &!input
                       ,cpatch%hite(ico),cpatch%pft(ico),transp       &!input
                       ,cpatch%leaf_psi(ico),cpatch%wood_psi(ico)     &!input
                       ,soil_psi,soil_cond,ipa,ico                    &!input
                       ,cpatch%wflux_wl(ico),cpatch%wflux_gw(ico)     &!output
                       ,cpatch%wflux_gw_layer(:,ico))                 !!output

            else
                cpatch%wflux_wl(ico) = 0.
                cpatch%wflux_gw(ico) = 0.
                cpatch%wflux_gw_layer(:,ico)  = 0.
            end if
        end do cohortloop

      end select

      
      
      
      
      
      ! Update Fast timescale output
      do ico = 1, cpatch%ncohorts
         cpatch%fmean_leaf_psi   (ico) = cpatch%fmean_leaf_psi   (ico)                     &
                                       + cpatch%leaf_psi         (ico) * dtlsm_o_frqsum
         cpatch%fmean_wood_psi   (ico) = cpatch%fmean_wood_psi   (ico)                     &
                                       + cpatch%wood_psi         (ico) * dtlsm_o_frqsum
         cpatch%fmean_leaf_water_int(ico) = cpatch%fmean_leaf_water_int(ico)               &
                                          + cpatch%leaf_water_int(ico) * dtlsm_o_frqsum
         cpatch%fmean_wood_water_int(ico) = cpatch%fmean_wood_water_int(ico)               &
                                          + cpatch%wood_water_int(ico) * dtlsm_o_frqsum
         if (cpatch%dmax_leaf_psi(ico) == 0.) then
             cpatch%dmax_leaf_psi(ico) =  cpatch%leaf_psi(ico)
         else
             cpatch%dmax_leaf_psi(ico) =  max(cpatch%dmax_leaf_psi(ico),cpatch%leaf_psi(ico))
         end if

         if (cpatch%dmin_leaf_psi(ico) == 0.) then
             cpatch%dmin_leaf_psi(ico) =  cpatch%leaf_psi(ico)
         else
             cpatch%dmin_leaf_psi(ico) =  min(cpatch%dmin_leaf_psi(ico),cpatch%leaf_psi(ico))
         end if

         if (cpatch%dmax_wood_psi(ico) == 0.) then
             cpatch%dmax_wood_psi(ico) =  cpatch%wood_psi(ico)
         else
             cpatch%dmax_wood_psi(ico) =  max(cpatch%dmax_wood_psi(ico),cpatch%wood_psi(ico))
         end if

         if (cpatch%dmin_wood_psi(ico) == 0.) then
             cpatch%dmin_wood_psi(ico) =  cpatch%wood_psi(ico)
         else
             cpatch%dmin_wood_psi(ico) =  min(cpatch%dmin_wood_psi(ico),cpatch%wood_psi(ico))
         end if

          cpatch%fmean_wflux_wl             (ico) = cpatch%fmean_wflux_wl(ico)         &
                                                  + cpatch%wflux_wl(ico)               &
                                                  * dtlsm_o_frqsum
          cpatch%fmean_wflux_gw             (ico) = cpatch%fmean_wflux_gw(ico)         &
                                                  + cpatch%wflux_gw(ico)               &
                                                  * dtlsm_o_frqsum
          do k = 1, nzg

              cpatch%fmean_wflux_gw_layer (k,ico) = cpatch%fmean_wflux_gw_layer(k,ico) &
                                                  + cpatch%wflux_gw_layer(k,ico)       &
                                                  * dtlsm_o_frqsum
          end do

       end do

      return

   end subroutine plant_hydro_driver
   !=======================================================================================!
   !=======================================================================================!


   !=======================================================================================!
   !=======================================================================================!
   ! SUBROUTINE CALC_PLANT_WATER_FLUX
   !> \brief Calculate water flow within plants driven by hydraulic laws
   !> \details This subroutine calculates ground to wood/root, and wood/root to
   !> leaf water flow from initial plant water potential and transpiration rates
   !> based on hydraulic laws. To simplify calculation and reduce numerial error,
   !> It is assumed that within one timestep, leaf water pool << wood water pool
   !> << soil water pool. Special cases are handled when these
   !> assumptions are incorrect such as grasses and very small seedlings.\n
   !>   Note that this subroutine only update fluxes but not the water potential
   !> or water content, which will be updated within the RK4 Integrator.
   !>   Currently, this subroutine works at DTLSM level, but it is
   !> straightforwad to implement the plant hydrodynamics within the RK4
   !> integration scheme. Yet, it is not yet tested how much extra computational
   !> cost would it incur\n
   !> References:\n
   !> [X16]   Diversity in plant hydraulic traits explains seasonal and inter--annual
   !> variations of vegetation dynamics in seasonally dry tropical forests
   !> X Xu, D Medvigy, JS Powers, JM Becknell, K Guan - New Phytologist, 2016
   !>
   !> \author Xiangtao Xu, 29 Jan. 2018
   !---------------------------------------------------------------------------------------!
   subroutine calc_plant_water_flux(dt                                  & !timestep
               ,sap_area,crown_area,krdepth,bleaf,bsap,broot,hite,ipft  & !plant input
               ,transp,leaf_psi,wood_psi                                & !plant input
               ,soil_psi,soil_cond                                      & !soil  input
               ,ipa,ico                                                 & !for debugging
               ,wflux_wl,wflux_gw,wflux_gw_layer)                       ! !flux  output
      use soil_coms       , only : slz8                 & ! intent(in)
                                 , dslz8                ! ! intent(in)
      use grid_coms       , only : nzg                  ! ! intent(in)
      use consts_coms     , only : pi18                 & ! intent(in)
                                 , lnexp_min8           ! ! intent(in)
      use rk4_coms        , only : tiny_offset          ! ! intent(in)
      use pft_coms        , only : leaf_water_cap       & ! intent(in) 
                                 , wood_water_cap       & ! intent(in)
                                 , leaf_psi_min         & ! intent(in)
                                 , wood_psi_min         & ! intent(in)
                                 , wood_psi50           & ! intent(in)
                                 , wood_Kmax            & ! intent(in)
                                 , wood_Kexp            & ! intent(in)
                                 , vessel_curl_factor   & ! intent(in)
                                 , root_beta            & ! intent(in)
                                 , SRA                  & ! intent(in)
                                 , C2B                  & ! intent(in)
                                 , hgt_min              ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real   ,                 intent(in)  :: dt              !time step [s]
      real   ,                 intent(in)  :: sap_area        !sapwood_area [m2]
      real   ,                 intent(in)  :: crown_area      !crown_area [m2]
      integer,                 intent(in)  :: krdepth         !maximum rooting depth
      real   ,                 intent(in)  :: bleaf           !leaf biomass [kgC]
      real   ,                 intent(in)  :: bsap            !sapwood biomass [kgC]
      real   ,                 intent(in)  :: broot           !fine root biomass [kgC]
      real   ,                 intent(in)  :: hite            !plant height [m]       
      integer,                 intent(in)  :: ipft            !plant functional type
      real   ,                 intent(in)  :: transp          !transpiration [kg/s]
      real   ,                 intent(in)  :: leaf_psi        !leaf water pot. [m]
      real   ,                 intent(in)  :: wood_psi        !wood water pot. [m]
      real   , dimension(nzg), intent(in)  :: soil_psi        !soil water pot. [m]
      real   , dimension(nzg), intent(in)  :: soil_cond       !soil water cond.[kg/m2/s]
      integer,                 intent(in)  :: ipa             !Patch index          
      integer,                 intent(in)  :: ico             !Cohort index         
      real   ,                 intent(out) :: wflux_wl        !wood-leaf flux  [kg/s]
      real   ,                 intent(out) :: wflux_gw        !ground-wood flux [kg/s]
      real   , dimension(nzg), intent(out) :: wflux_gw_layer  !wflux_gw for each soil layer

      !----- Local Vars  --------------------------------------------------------------------!
      !----- temporary double precision variables for input and output variables
      real(kind=8)                 :: dt_d              
      real(kind=8)                 :: sap_area_d        
      real(kind=8)                 :: crown_area_d      
      real(kind=8)                 :: bleaf_d           
      real(kind=8)                 :: bsap_d            
      real(kind=8)                 :: broot_d           
      real(kind=8)                 :: hite_d            
      real(kind=8)                 :: transp_d          
      real(kind=8)                 :: leaf_psi_d        
      real(kind=8)                 :: wood_psi_d        
      real(kind=8), dimension(nzg) :: soil_psi_d        
      real(kind=8), dimension(nzg) :: soil_cond_d        
      real(kind=8)                 :: wflux_wl_d        
      real(kind=8)                 :: wflux_gw_d        
      real(kind=8), dimension(nzg) :: wflux_gw_layer_d  
      !----- Auxiliary variables
      real(kind=8)                          :: exp_term             !exponent term
      real(kind=8)                          :: ap                   ![s-1]
      real(kind=8)                          :: bp                   ![m s-1]
      real(kind=8)                          :: stem_cond            !stem conductance
      real(kind=8)                          :: plc                  !plant loss of conductance
      real(kind=8)                          :: c_leaf               !leaf water capacitance
      real(kind=8)                          :: c_stem               !stem water capacitance
      real(kind=8)                          :: RAI                  !root area index
      real(kind=8)                          :: root_frac            !fraction of roots
      real(kind=8)                          :: proj_leaf_psi        !projected leaf water pot.
      real(kind=8)                          :: proj_wood_psi        !projected wood water pot. 
      real(kind=8)                          :: gw_cond              !g->w water conductivity
      real(kind=8)                          :: org_wood_psi         !used for small tree
      real(kind=8)                          :: org_leaf_psi         !used for small tree
      real(kind=8)                          :: weighted_soil_psi
      real(kind=8)                          :: weighted_gw_cond
      real(kind=8)                          :: above_layer_depth
      real(kind=8)                          :: current_layer_depth
      real(kind=8)                          :: total_water_supply
      real(kind=8)      , dimension(nzg)    :: layer_water_supply
      character(len=13) , parameter         :: efmt       = '(a,1x,es12.5)'
      character(len=9)  , parameter         :: ifmt       = '(a,1x,i5)'
      character(len=9)  , parameter         :: lfmt       = '(a,1x,l1)'
      !----- variables for loops
      integer                               :: k
      integer,parameter                     :: dco        = 0
      !--------------- Flags
      logical,parameter                     :: debug_flag = .false.
      logical                               :: small_tree_flag
      logical                               :: zero_flow_flag
      logical                               :: error_flag
      !----- External function ------------------------------------------------------------!
      real(kind=4), external                :: sngloff           ! Safe dble 2 single precision
      !------------------------------------------------------------------------------------!

   
      !--------------------------------------------------------------------------
      ! Convert all input state vars to double precision               
      !--------------------------------------------------------------------------
      dt_d          = dble(dt)
      sap_area_d    = dble(sap_area)
      crown_area_d  = dble(crown_area)
      bleaf_d       = dble(bleaf)
      bsap_d        = dble(bsap)
      broot_d       = dble(broot)
      hite_d        = dble(hite) 
      transp_d      = dble(transp)
      leaf_psi_d    = dble(leaf_psi)                            
      wood_psi_d    = dble(wood_psi)
      soil_psi_d    = dble(soil_psi)
      soil_cond_d   = dble(soil_cond)
   



      !----------------------------------------------------------------------
      ! Update plant hydrodynamics
      ! The regular slover assumes stem water pool is way larger than the leaf water
      ! pool. In cases where leaf water pool is of similar magnitude to stem
      ! water pool (small seedlings and grasses), leaf water potential is forced
      ! to be the same as stem water potential to maintain numerical stability
      ! in sacrifice of baises
      !
      ! Water flow is calculated from canopy to roots
      ! Positive flow means upward (g->w, w->l, l->air)
      !-----------------------------------------------------------------------

      ! initiate proj_psi as the starting psi
      proj_leaf_psi = leaf_psi_d 
      proj_wood_psi = wood_psi_d


      org_wood_psi = wood_psi_d
      org_leaf_psi = leaf_psi_d



      !--------------------------------------------------------------------------
      ! 1. Calculate wood/stem/root to leaf water flow
      !--------------------------------------------------------------------------
        
      ! First, check the relative magnitude of leaf and sapwood water pool
      ! If it is a small tree/grass, force psi_leaf to be the same as psi_stem
      c_leaf = dble(leaf_water_cap(ipft) * C2B) * bleaf_d        ! kg H2O / m
      c_stem = dble(wood_water_cap(ipft) * C2B) * (broot_d + bsap_d) ! kg H2O / m

      ! Consider a tree is too small if c_leaf is larger than half of c_stem
      ! This is an arbitrary threshold. Users are welcomed to modify this term
      ! if leaf_psi has strong oscillations from each timestep to another

      ! we also assume it is a small tree if the tree is too short
      small_tree_flag = ((c_leaf > (c_stem / 2.d0)) .or. (hite_d == hgt_min(ipft)))





      if (small_tree_flag) then
          ! 1.1 small tree, force leaf_psi to be the same as wood_psi
          
          ! calculate the new veg_psi of mixing leaf and wood
          wood_psi_d   = (c_leaf * leaf_psi_d + c_stem * wood_psi_d) / (c_leaf+c_stem)
          leaf_psi_d   = wood_psi_d


          ! use c_leaf+c_stem as the total capacitance to solve stem_psi later
          c_stem = c_stem + c_leaf

          ! in this case, we temporally assign transpiration as wflux_wl_d since
          ! the leaf and wood are treated as a whole. The value will be
          ! recalculated once we got the projected water potential

          wflux_wl_d   = transp_d

      else
          ! 1.2 Regular case, big trees
  
          ! Special case handling.... Maybe there are better ways to avoid
          ! these traps in a more systematic way. Currently, I just do it mannually
  
          ! Special case 0 if negative wflux_wl overchange wood storage?
          ! Should not be a problem if initiate leaf waterpotential a little
          ! lower than the stem water potential.....
  
          ! Special case 1, if there are no leaves, we zero the flow
  
          ! Sepcial case 2, if both wood and leaves are very dry, either wood
          ! cannot support upward sapflow or leaf cannot support downward
          ! flow... This could happen for dying trees experiencing extreme
          ! drought. In this case, we zero the flow

          ! Special case 3, if the cohort just grows out of 'small tree status'.
          ! Their leaves can be over-charged with water because gravitational
          ! effect was not considered for leaf water potential of small trees.
          ! As a result, this can lead to a down-ward sapflow, and potentially
          ! over-charging the sapwood. We need to zero the flow in this case as
          ! well, until leaf_psi_d drops below wood_psi_d - hite_d
  
          zero_flow_flag = (c_leaf == 0.d0)                          .or.      & ! Case 1
                           (leaf_psi_d >= (wood_psi_d - hite_d) .and.          &
                            leaf_psi_d <= dble(leaf_psi_min(ipft)))  .or.      & ! Case 2
                           (leaf_psi_d <= (wood_psi_d - hite_d) .and.          &
                            wood_psi_d <= dble(wood_psi_min(ipft)))  .or.      & ! Case 2
                           (leaf_psi_d >  (wood_psi_d - hite_d))               ! ! Case 3
  
                          
          if (zero_flow_flag) then
              ! 1.2.1 Special case. No need to calculate sapflow
              wflux_wl_d = 0.d0

              ! proj_leaf_psi is only depdent on transp
              if (c_leaf > 0.) then
                  proj_leaf_psi = leaf_psi_d - transp_d * dt_d / c_leaf
              else
                  proj_leaf_psi = leaf_psi_d
              end if
  
          else
              ! We do need to calculate sapflow

              ! calculate plant loss of conductivity [dimensionless]
              plc = 1.d0 / (1.d0 +                                                  &
                    (wood_psi_d / dble(wood_psi50(ipft))) ** dble(wood_Kexp(ipft)))
              ! calculate stem conductance [kg / s]
              stem_cond = dble(wood_Kmax(ipft)) * plc                 & ! kg/m/s
                        * sap_area_d                                  & ! conducting area m2
                        / (hite_d * dble(vessel_curl_factor(ipft)))   ! ! conducting length m
  
              if (stem_cond == 0.) then
              ! 1.2.2 Special case when there are no stem conductivity
                  wflux_wl_d = 0.d0
              else
                  ! 1.2.3 Normal case with positive c_leaf and positive stem_cond
                  ! Check ref X16 for derivation of the equations
  
                  ap = - stem_cond / c_leaf
                  ! the unit of ap is s-1
  
                  bp = ((wood_psi_d - hite_d) * stem_cond - transp_d) &
                     / c_leaf                                                           
                  ! the unit of bp is m s-1
                  exp_term = exp(max(ap * dt_d,lnexp_min8))
  
                  ! project the final leaf psi
                  proj_leaf_psi = ((ap * leaf_psi_d + bp) * exp_term - bp) / ap
  
                  ! calculate the average sapflow rate  within the time step [kg H2O /s]
                  wflux_wl_d = (proj_leaf_psi - leaf_psi_d) * c_leaf / dt_d + transp_d
              end if
   
          end if
      end if

        
      !--------------------------------------------------------------------------
      ! 2. Calculate ground -> wood/stem/root water flow
      !--------------------------------------------------------------------------
        weighted_soil_psi  = 0.d0
        weighted_gw_cond = 0.d0
        layer_water_supply = 0.d0
        total_water_supply = 0.d0

        ! loop over all soil layers to get the aggregated water conductance
        do k = krdepth,nzg
            current_layer_depth = -slz8(k)
            if (k+1 .le. nzg) then
                above_layer_depth = -slz8(k+1)
            else
                above_layer_depth = 0.d0
            end if

            ! calcualte the root fraction of this layer
            root_frac = &
                (dble(root_beta(ipft)) ** (above_layer_depth   / (-slz8(krdepth)))  &
                -dble(root_beta(ipft)) ** (current_layer_depth / (-slz8(krdepth)))  &
                )
                                       
            !  Calculate RAI in each layer
            !  Assume root can extent to an area 4 times of crown area (twice as
            !  much as crown radius)
            if (crown_area_d == 0.d0) then
                RAI = 0.d0
            else
                RAI = broot_d * dble(SRA(ipft)) * root_frac     & !m2
                    / (4.d0 * crown_area_d)                         !m2
            end if

            !  Calculate soil-root water conductance kg H2O / m / s
            !  Based on Katul et al. 2003 PCE
            gw_cond = soil_cond_d(k) * sqrt(RAI) / (pi18 * dslz8(k))  & ! kg H2O / m3 / s
                    * (4.d0 * crown_area_d)           ! ! conducting area  m2
            
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! disable hydraulic redistribution
            ! assume roots will shut down if they are going to lose water to
            ! soil
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            if (soil_psi_d(k) <= wood_psi_d) then
                gw_cond = 0.d0
            end if
            
            ! Calculate weighted conductance, weighted psi, and
            ! water_supply_layer_frac
            weighted_gw_cond   = weighted_gw_cond + gw_cond ! kg H2O /m /s
            weighted_soil_psi  = weighted_soil_psi + gw_cond * soil_psi_d(k) ! kgH2O/s

            layer_water_supply(k) = gw_cond * (soil_psi_d(k) - wood_psi_d)!kgH2O s-1 
        end do

        ! Now we can calculate ground->wood water flow
        ! First we handle special cases
        zero_flow_flag = (c_stem == 0.d0)       .or.      & ! No sapwood or fine  roots
                         (weighted_gw_cond == 0.d0)     ! ! soil is drier than wood

        if (zero_flow_flag) then
            ! no need to calculate water flow
            ! wood psi is only dependent on sapflow
            wflux_gw_d    = 0.d0
            if (c_stem > 0.) then
                proj_wood_psi = wood_psi_d - wflux_wl_d * dt_d / c_stem
            else
                proj_wood_psi = wood_psi_d
            end if
        else
            ! calculate the average soil water uptake
            ap = - weighted_gw_cond  / c_stem 
            !the unit of ap is s-1

            bp = (weighted_soil_psi - wflux_wl_d) / c_stem
            !the unit of bp is m s-1

            exp_term        = exp(max(ap * dt_d,lnexp_min8))
            proj_wood_psi   = ((ap * wood_psi_d + bp) * exp_term - bp) / ap
            wflux_gw_d     = (proj_wood_psi - wood_psi_d) * c_stem  / dt_d + wflux_wl_d
        end if

        ! We need to re-calculate the water fluxes for small tree scenarios
        if (small_tree_flag) then
            ! for small tree

            ! first wflux_gw is correct
            ! no need to update

            ! but we need to update wflux_wl_d
            proj_leaf_psi = proj_wood_psi
            wflux_wl_d    = (proj_leaf_psi - org_leaf_psi) &
                          * c_leaf / dt_d + transp_d

        end if

       


        !----------------------------------------------------------------------
        ! Now estimate the water uptake from each layer based on
        ! layer_water_supply
        !-----------------------------------------------------------------------
        if (sum(layer_water_supply) == 0.d0) then
            wflux_gw_layer_d = 0.d0
        else
            wflux_gw_layer_d = layer_water_supply / sum(layer_water_supply) * &
                               wflux_gw_d
        end if




      !--------------------------------------------------------------------------
      ! Handling Potential Errors and Help Debugging
      !--------------------------------------------------------------------------
      error_flag = (isnan(wflux_wl_d) .or. isnan(wflux_gw_d))       & ! NaN values
               .or.(proj_leaf_psi > 0. .or. proj_wood_psi > 0.)     & ! psi is positive
               .or.(leaf_psi_d > 0. .or. wood_psi_d > 0.)

      ! I copy the error printing from rk4_misc.f90
      if((debug_flag .and. (dco == 0 .or. ico == dco)) .or. error_flag) then
         write (unit=*,fmt='(a)') ' '
         write (unit=*,fmt='(92a)') ('=',k=1,92)
         write (unit=*,fmt='(92a)') ('=',k=1,92)
         write (unit=*,fmt='(a)'  ) ' plant hydrodynamics inconsistency detected!!'
         write (unit=*,fmt='(92a)') ('-',k=1,92)
         write (unit=*,fmt=ifmt   ) ' + IPA              =',ipa

         write (unit=*,fmt='(a)'  ) ' '
         write (unit=*,fmt=ifmt   ) ' + ICO              =',ico
         write (unit=*,fmt=ifmt   ) ' + PFT              =',ipft
         write (unit=*,fmt=ifmt   ) ' + KRDEPTH          =',krdepth
         write (unit=*,fmt=efmt   ) ' + HEIGHT           =',hite
         write (unit=*,fmt=lfmt   ) ' + SMALL_TREE_FLAG  =',small_tree_flag

         write (unit=*,fmt='(a)'  ) ' '
         write (unit=*,fmt=efmt   ) ' + BLEAF            =',bleaf                      
         write (unit=*,fmt=efmt   ) ' + BSAPWOOD         =',bsap                       
         write (unit=*,fmt=efmt   ) ' + BROOT            =',broot                      
         write (unit=*,fmt=efmt   ) ' + SAPWOOD_AREA     =',sap_area                   
         write (unit=*,fmt=efmt   ) ' + CROWN_AERA       =',crown_area

         write (unit=*,fmt='(a)'  ) ' '
         write (unit=*,fmt=efmt   ) ' + TRANSP           =',transp                     
         write (unit=*,fmt=efmt   ) ' + LEAF_PSI (INPUT) =',leaf_psi                   
         write (unit=*,fmt=efmt   ) ' + WOOD_PSI (INPUT) =',wood_psi                   
         write (unit=*,fmt=efmt   ) ' + LEAF_PSI (PROJ.) =',proj_leaf_psi              
         write (unit=*,fmt=efmt   ) ' + WOOD_PSI (PROJ.) =',proj_wood_psi              
         write (unit=*,fmt=efmt   ) ' + WFLUX_GW         =',wflux_gw_d                 
         write (unit=*,fmt=efmt   ) ' + WFLUX_WL         =',wflux_wl_d                 
         write (unit=*,fmt='(a)'  ) ' + WFLUX_GW_LAYER   ='            
         do k = 1, nzg
             write (unit=*,fmt='(i5,1x,es12.5)') k, wflux_gw_layer_d(k)                 
         end do

         write (unit=*,fmt='(a)'  ) ' '
         write (unit=*,fmt='(a)'  ) ' + SOIL_PSI         ='            
         do k = 1, nzg
             write (unit=*,fmt='(i5,1x,es12.5)') k, soil_psi(k)                 
         end do
         write (unit=*,fmt='(92a)') ('=',k=1,92)
         write (unit=*,fmt='(92a)') ('=',k=1,92)
         write (unit=*,fmt='(a)') ' '

         if (error_flag) then 
             call fatal_error('Plant Hydrodynamics is wrong',&
                              'calc_plant_water_flux'       ,&
                              'plant_hydro.f90')
         end if
      end if
 


      !--------------------------------------------------------------------------
      ! Copy all the results to output variables
      !--------------------------------------------------------------------------
      wflux_wl = sngloff(wflux_wl_d,tiny_offset)
      ! wflux_gw = sngloff(wflux_gw_d,tiny_offset)
      do k = 1, nzg
        wflux_gw_layer(k) = sngloff(wflux_gw_layer_d(k),tiny_offset)
      end do
      ! MLO -> XX. I imposed wflux_gw to be the same as the sum of layers to ensure
      !            that water would be conserved.  Please check that this is ok.
      wflux_gw = sum(wflux_gw_layer)


      return
   end subroutine calc_plant_water_flux
   !=======================================================================================!
   !=======================================================================================!







   !=======================================================================================!
   !=======================================================================================!
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
      use pft_coms          ,   only : leaf_water_cap       & ! intent(in)
                                     , wood_water_cap       & ! intent(in)
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
                        * leaf_water_cap(ipft)      & ! [kg H2O/kg biomass/m]
                        / leaf_water_sat(ipft)      ! ! [kg H2O/kg biomass]

      ! same for wood
      wood_rwc  =   1.  + wood_psi                  & ! [m]
                        * wood_water_cap(ipft)      & ! [kg H2O/kg biomass/m]
                        / wood_water_sat(ipft)      ! ! [kg H2O/kg biomass]

      return
   end subroutine psi2rwc
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !  SUBROUTINE: RWC2PSI
   !> \brief Convert relative water content to water potential
   !> \details The inverse of psi2rwc
   !=======================================================================================!
   subroutine rwc2psi(leaf_rwc,wood_rwc,ipft,leaf_psi,wood_psi)
      use pft_coms          ,   only : leaf_water_cap       & ! intent(in)
                                     , wood_water_cap       & ! intent(in)
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
      leaf_psi  =   (leaf_rwc - 1.) * leaf_water_sat(ipft) / leaf_water_cap(ipft)
      ! same for wood
      wood_psi  =   (wood_rwc - 1.) * wood_water_sat(ipft) / wood_water_cap(ipft)

      return
   end subroutine rwc2psi
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !  SUBROUTINE: RWC2TW            
   !> \brief Convert relative water content to total water for both leaf and wood
   !> \details total water content = \n
   !>  relative water content * water_content_at_saturation * biomass
   !> \warning In the hydro version, bsapwood is set as 0 and bdead is assumedto contain
   !> both sapwood and heart wood. Root is counted as wood.  When dynamic hydraulics is
   !> turned off, we account for water in sapwood using the sapwood biomass definition.
   !=======================================================================================!
   subroutine rwc2tw(leaf_rwc,wood_rwc,bleaf,bsapwooda,bsapwoodb,bdeada,bdeadb,broot,dbh   &
                    ,ipft,leaf_water_int,wood_water_int)
      use pft_coms       , only : leaf_water_sat      & ! intent(in)
                                , wood_water_sat      & ! intent(in)
                                , C2B                 ! ! intent(in)
      use allometry      , only : dbh2sf              ! ! function
      use physiology_coms, only : plant_hydro_scheme  ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real      , intent(in)    ::  leaf_rwc       ! Relative water content of leaves  [0-1]
      real      , intent(in)    ::  wood_rwc       ! Relative water content of wood    [0-1]
      real      , intent(in)    ::  bleaf          ! Biomass of leaf                   [kgC]
      real      , intent(in)    ::  bsapwooda      ! Aboveground sapwood biomass       [kgC]
      real      , intent(in)    ::  bsapwoodb      ! Belowground sapwood biomass       [kgC]
      real      , intent(in)    ::  bdeada         ! Aboveground (heart)wood biomass   [kgC]
      real      , intent(in)    ::  bdeadb         ! Belowground (heart)wood biomass   [kgC]
      real      , intent(in)    ::  broot          ! Biomass of fine root              [kgC]
      real      , intent(in)    ::  dbh            ! Diameter at breast height         [ cm]
      integer   , intent(in)    ::  ipft           ! Plant functional type             [  -]
      real      , intent(out)   ::  leaf_water_int ! Total internal water of leaf      [ kg]
      real      , intent(out)   ::  wood_water_int ! Total internal water of wood      [ kg]
      !----- Local variables. -------------------------------------------------------------!
      real                      ::  sap_frac    ! Fraction of sapwood to basal area    [0-1]
      !------------------------------------------------------------------------------------!


      !----- Leaf.  This is the same, regardless of the plant hydraulic scheme. -----------!
      leaf_water_int    =   leaf_rwc * leaf_water_sat(ipft) * bleaf * C2B 
      !------------------------------------------------------------------------------------!


      !----- Wood.  Check the scheme to decide between sapwood fraction or biomass. -------!
      select case (plant_hydro_scheme)
      case (0)
         !----- Use sapwood biomass to obtain sapwood internal water. ---------------------!
         wood_water_int = wood_rwc * wood_water_sat(ipft) * C2B                            &
                        * (broot + bsapwooda + bsapwoodb )
         !---------------------------------------------------------------------------------!
      case default
         !----- Find the sapwood fraction. ------------------------------------------------!
         sap_frac       = dbh2sf(dbh,ipft)
         !---------------------------------------------------------------------------------!

         !----- Total water only includes live biomass (fine roots and sapwood). ----------!
         wood_water_int = wood_rwc * wood_water_sat(ipft) * C2B                            &
                        * ( broot + (bdeada + bdeadb + bsapwooda + bsapwoodb) * sap_frac )
         !---------------------------------------------------------------------------------!
      end select
      !------------------------------------------------------------------------------------!


      return
   end subroutine rwc2tw
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !  SUBROUTINE: TW2RWC            
   !> \brief Convert total water to relative water content for both leaf and wood
   !> \details the inverse of rwc2tw \n
   !=======================================================================================!
   subroutine tw2rwc(leaf_water_int,wood_water_int,bleaf,bsapwooda,bsapwoodb,bdeada,bdeadb &
                    ,broot,dbh,ipft,leaf_rwc,wood_rwc)
      use pft_coms       , only : leaf_water_sat     & ! intent(in)
                                , wood_water_sat     & ! intent(in)
                                , C2B                ! ! intent(in)
      use allometry      , only : dbh2sf             ! ! function
      use physiology_coms, only : plant_hydro_scheme ! ! intent(in)
      use consts_coms    , only : tiny_num           ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real      , intent(in)    ::  leaf_water_int ! Total internal water of leaf      [ kg]
      real      , intent(in)    ::  wood_water_int ! Total internal water of wood      [ kg]
      real      , intent(in)    ::  bleaf          ! Biomass of leaf                   [kgC]
      real      , intent(in)    ::  bsapwooda      ! Aboveground sapwood biomass       [kgC]
      real      , intent(in)    ::  bsapwoodb      ! Belowground sapwood biomass       [kgC]
      real      , intent(in)    ::  bdeada         ! Aboveground (heart)wood biomass   [kgC]
      real      , intent(in)    ::  bdeadb         ! Belowground (heart)wood biomass   [kgC]
      real      , intent(in)    ::  broot          ! Biomass of fine root              [kgC]
      real      , intent(in)    ::  dbh            ! Diameter at breast height         [ cm]
      integer   , intent(in)    ::  ipft           ! Plant functional type             [  -]
      real      , intent(out)   ::  leaf_rwc       ! Relative water content of leaves  [0-1]
      real      , intent(out)   ::  wood_rwc       ! Relative water content of wood    [0-1]
      !----- Local variables --------------------------------------------------------------!
      real                      :: tot_water_sat
      real                      :: sap_frac        ! Fraction of sapwood to basal area [0-1]
      !------------------------------------------------------------------------------------!


      !----- Leaf.  This is the same, regardless of the plant hydraulic scheme. -----------!
      tot_water_sat = leaf_water_sat(ipft) * C2B * bleaf
      if (tot_water_sat > tiny_num) then
         leaf_rwc  = leaf_water_int / tot_water_sat
      else
         leaf_rwc  = 0.
      end if
      !------------------------------------------------------------------------------------!


      !----- Wood.  Check the scheme to decide between sapwood fraction or biomass. -------!
      select case (plant_hydro_scheme)
      case (0)
         !----- Use sapwood biomass to obtain sapwood internal water. ---------------------!
         tot_water_sat = wood_water_sat(ipft) * C2B                                        &
                       * (broot + bsapwooda + bsapwoodb)
         !---------------------------------------------------------------------------------!
      case default
         !----- Find the sapwood fraction. ------------------------------------------------!
         sap_frac = dbh2sf(dbh,ipft)
         !---------------------------------------------------------------------------------!


         !----- Find the sapwood fraction. ------------------------------------------------!
         tot_water_sat = wood_water_sat(ipft) * C2B                                        &
                       * (broot + (bdeada + bdeadb + bsapwooda + bsapwoodb) * sap_frac)
         !---------------------------------------------------------------------------------!
      end select
      !------------------------------------------------------------------------------------!


      !------ Make sure the denominator is not zero. --------------------------------------!
      if (tot_water_sat > tiny_num) then
         wood_rwc = wood_water_int / tot_water_sat
      else
         wood_rwc = 0.
      end if
      !------------------------------------------------------------------------------------!

      return
   end subroutine tw2rwc
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !  SUBROUTINE: PSI2TW            
   !> \breif Convert water potential to total water for both leaf and wood
   !=======================================================================================!
   subroutine psi2tw(leaf_psi,wood_psi,bleaf,bsapwooda,bsapwoodb,bdeada,bdeadb,broot,dbh   &
                    ,ipft,leaf_water_int,wood_water_int)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real      , intent(in)    ::  leaf_psi       ! Water potential of leaves         [  m]
      real      , intent(in)    ::  wood_psi       ! Water potential of wood           [  m]
      real      , intent(in)    ::  bleaf          ! Biomass of leaf                   [kgC]
      real      , intent(in)    ::  bsapwooda      ! Aboveground sapwood biomass       [kgC]
      real      , intent(in)    ::  bsapwoodb      ! Belowground sapwood biomass       [kgC]
      real      , intent(in)    ::  bdeada         ! Aboveground (heart)wood biomass   [kgC]
      real      , intent(in)    ::  bdeadb         ! Belowground (heart)wood biomass   [kgC]
      real      , intent(in)    ::  broot          ! Biomass of fine root              [kgC]
      real      , intent(in)    ::  dbh            ! Diameter at breast height         [ cm]
      integer   , intent(in)    ::  ipft           ! Plant functional type             [  -]
      real      , intent(out)   ::  leaf_water_int ! Total internal water of leaf      [ kg]
      real      , intent(out)   ::  wood_water_int ! Total internal water of wood      [ kg]
      !----- Local Variables --------------------------------------------------------------!
      real                      ::  leaf_rwc       ! Relative water content of leaf    [  -]
      real                      ::  wood_rwc       ! Relative water content of wood    [  -]
      !------------------------------------------------------------------------------------!

      ! first convert to rwc
      call psi2rwc(leaf_psi,wood_psi,ipft,leaf_rwc,wood_rwc)
      ! second convert to tw
      call rwc2tw(leaf_rwc,wood_rwc,bleaf,bsapwooda,bsapwoodb,bdeada,bdeadb,broot,dbh,ipft &
                 ,leaf_water_int,wood_water_int)

      return

   end subroutine psi2tw
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !  SUBROUTINE: TW2PSI            
   !> \brief Convert total water to water potential for both leaf and wood
   !> \details the inverse of psi2tw \n
   !=======================================================================================!
   subroutine tw2psi(leaf_water_int,wood_water_int,bleaf,bsapwooda,bsapwoodb,bdeada,bdeadb &
                    ,broot,dbh,ipft,leaf_psi,wood_psi)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real      , intent(in)    ::  leaf_water_int ! Total internal water of leaf      [ kg]
      real      , intent(in)    ::  wood_water_int ! Total internal water of wood      [ kg]
      real      , intent(in)    ::  bleaf          ! Biomass of leaf                   [kgC]
      real      , intent(in)    ::  bsapwooda      ! Aboveground sapwood biomass       [kgC]
      real      , intent(in)    ::  bsapwoodb      ! Belowground sapwood biomass       [kgC]
      real      , intent(in)    ::  bdeada         ! Aboveground (heart)wood biomass   [kgC]
      real      , intent(in)    ::  bdeadb         ! Belowground (heart)wood biomass   [kgC]
      real      , intent(in)    ::  broot          ! Biomass of fine root              [kgC]
      real      , intent(in)    ::  dbh            ! Diameter at breast height         [ cm]
      integer   , intent(in)    ::  ipft           ! Plant functional type             [  -]
      real      , intent(out)   ::  leaf_psi       ! Water potential of leaves         [  m]
      real      , intent(out)   ::  wood_psi       ! Water potential of wood           [  m]
      !----- Local Variables --------------------------------------------------------------!
      real                      ::  leaf_rwc       ! Relative water content of leaf    [  -]
      real                      ::  wood_rwc       ! Relative water content of wood    [  -]
      !------------------------------------------------------------------------------------!

      ! first convert to rwc
      call tw2rwc(leaf_water_int,wood_water_int,bleaf,bsapwooda,bsapwoodb,bdeada,bdeadb    &
                 ,broot,dbh,ipft,leaf_rwc,wood_rwc)
      ! second convert to psi
      call rwc2psi(leaf_rwc,wood_rwc,ipft,leaf_psi,wood_psi)

      return

   end subroutine tw2psi
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !  SUBROUTINE twi2twe
   !> \brief  Intensive to extensive internal water converter.
   !> \details This subroutine converts intensive internal water (kg/plant) to extensive
   !>          water content (kg/m2).  To avoid energy leaks, we assume that all water
   !>          stored in sapwood is included in the heat capacity.  This is not the most
   !>          elegant solution, and in the future, we should make it only the fraction
   !>          associated with branches, but this requires additional changes in the
   !>          budget checks.
   !> \author Marcos Longo 08 Sep 2019
   !---------------------------------------------------------------------------------------!
   subroutine twi2twe(leaf_water_int,wood_water_int,nplant,leaf_water_im2,wood_water_im2)
      use physiology_coms      , only : plant_hydro_scheme  ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real      , intent(in)    ::  leaf_water_int ! Total internal water of leaf  [ kg/pl]
      real      , intent(in)    ::  wood_water_int ! Total internal water of wood  [ kg/pl]
      real      , intent(in)    ::  nplant         ! Stem density                  [ pl/m2]
      real      , intent(out)   ::  leaf_water_im2 ! Extensive leaf internal water [ kg/m2]
      real      , intent(out)   ::  wood_water_im2 ! Water potential of wood       [ kg/m2]
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Convert water from intensive to extensive.                                    !
      !------------------------------------------------------------------------------------!
      select case (plant_hydro_scheme)
      case (0)
         !----- Plant hydraulics is disabled.  Make this number zero. ---------------------!
         leaf_water_im2 = 0.
         wood_water_im2 = 0.
         !---------------------------------------------------------------------------------!
      case default
         !----- Plant hydraulics is enabled.  Scale the internal water. -------------------!
         leaf_water_im2 = nplant * leaf_water_int
         wood_water_im2 = nplant * wood_water_int
         !---------------------------------------------------------------------------------!
      end select
      !------------------------------------------------------------------------------------!

      return
   end subroutine twi2twe
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !  SUBROUTINE twe2twi
   !> \brief  Extensive to extensive internal water converter.
   !> \details This subroutine converts extensive internal water (kg/m2) to extensive 
   !>          water content (kg/plant).  To avoid energy leaks, we assume that all water
   !>          stored in wood is included in the heat capacity.  In the future, we should
   !>          make it only the fraction associated with branches, but this requires 
   !>          additional changes in the budget checks.
   !> \author Marcos Longo 08 Sep 2019
   !---------------------------------------------------------------------------------------!
   subroutine twe2twi(leaf_water_im2,wood_water_im2,nplant,ipft,bleaf,broot,bsapwooda      &
                     ,bsapwoodb,leaf_water_int,wood_water_int)
      use physiology_coms      , only : plant_hydro_scheme  ! ! intent(in)
      use pft_coms             , only : C2B                 & ! intent(in)
                                      , leaf_water_sat      & ! intent(in)
                                      , wood_water_sat      ! ! intent(in)
      use consts_coms          , only : tiny_num            ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real      , intent(in)    ::  leaf_water_im2 ! Extensive leaf internal water [ kg/m2]
      real      , intent(in)    ::  wood_water_im2 ! Water potential of wood       [ kg/m2]
      real      , intent(in)    ::  nplant         ! Stem density                  [ pl/m2]
      integer   , intent(in)    ::  ipft           ! Plant functional type         [     -]
      real      , intent(in)    ::  bleaf          ! Leaf biomass                  [kgC/pl]
      real      , intent(in)    ::  broot          ! Fine-root biomass             [kgC/pl]
      real      , intent(in)    ::  bsapwooda      ! Aboveground sapwood biomass   [kgC/pl]
      real      , intent(in)    ::  bsapwoodb      ! Belowground sapwood biomass   [kgC/pl]
      real      , intent(inout) ::  leaf_water_int ! Total internal water of leaf  [ kg/pl]
      real      , intent(inout) ::  wood_water_int ! Total internal water of wood  [ kg/pl]
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Convert water from intensive to extensive.                                    !
      !------------------------------------------------------------------------------------!
      select case (plant_hydro_scheme)
      case (0)
         !----- Plant hydraulics is disabled.  Assume saturation. -------------------------!
         leaf_water_int = leaf_water_sat(ipft) * C2B * bleaf
         wood_water_int = wood_water_sat(ipft) * C2B * (broot + bsapwooda + bsapwoodb)
         !---------------------------------------------------------------------------------!
      case default
         !----- Plant hydraulics is enabled.  Scale the internal water. -------------------!
         leaf_water_int = leaf_water_im2 / nplant
         wood_water_int = wood_water_im2 / nplant
         !---------------------------------------------------------------------------------!
      end select
      !------------------------------------------------------------------------------------!

      return
   end subroutine twe2twi
   !=======================================================================================!
   !=======================================================================================!
end module plant_hydro

!==========================================================================================!
!==========================================================================================!
