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

   !------ Tolerance for minimum checks.  -------------------------------------------------!
   real(kind=4), parameter :: tol_buff   = 1.e-4
   real(kind=8), parameter :: tol_buff_d = dble(tol_buff)

   real(kind=4), parameter :: om_buff   = 1. - tol_buff
   real(kind=8), parameter :: om_buff_d = dble(om_buff)

   real(kind=4), parameter :: op_buff   = 1. + tol_buff
   real(kind=8), parameter :: op_buff_d = dble(op_buff)

   real(kind=4), parameter :: mg_safe   = 0.05
   real(kind=4), parameter :: op_safe   = 1. + mg_safe
   real(kind=4), parameter :: om_safe   = 1. - mg_safe
   !---------------------------------------------------------------------------------------!

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
                                      , dtlsm_o_frqsum         & ! intent(in)
                                      , current_time           ! ! intent(in)
      use soil_coms            , only : soil                   & ! intent(in)
                                      , matric_potential       & ! function
                                      , hydr_conduct           ! ! function
      use grid_coms            , only : nzg                    ! ! intent(in)
      use consts_coms          , only : pio4                   & ! intent(in)
                                      , wdns                   ! ! intent(in)
      use allometry            , only : dbh2sf                 ! ! function
      use physiology_coms      , only : plant_hydro_scheme     ! ! intent(in)
      use pft_coms             , only : C2B                    & ! intent(in)
                                      , leaf_water_cap         & ! intent(in)
                                      , leaf_psi_min           & ! intent(in)
                                      , small_psi_min          ! ! intent(in)

      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(sitetype)        , target      :: csite
      integer               , intent(in)  :: ipa
      integer,dimension(nzg), intent(in)  :: ntext_soil
      !----- Local Vars  ------------------------------------------------------------------!
      type(patchtype)       , pointer     :: cpatch      !< patch strcture
      real                                :: swater_min  !< Min. soil moisture for condct.
      real                                :: swater_max  !< Max. soil moisture for condct.
      real                                :: swater_use  !< soil moisture
      integer                             :: nsoil       !< soil type for soil
      integer                             :: k           !< iterator for soil lyr
      integer                             :: ico         !< iterator for cohort
      integer                             :: ipft        !< PFT index
      real ,dimension(nzg)                :: soil_psi    !< soil water potential   [      m]
      real ,dimension(nzg)                :: soil_cond   !< soil water conductance [kg/m2/s]
      real                                :: sap_frac    !< sapwood fraction       [    ---]
      real                                :: sap_area    !< sapwood area           [     m2]
      real                                :: bsap        !< sapwood biomass        [    kgC]
      real                                :: transp      !< transpiration rate     [   kg/s]
      real                                :: c_leaf      !< leaf capacitance       [   kg/m]
      logical                             :: track_hydraulics !< whether track hydraulics
      !----- Variables for debugging purposes ---------------------------------------------!
      integer, parameter                  :: dco        = 0 ! the cohort to debug
      logical, dimension(3)               :: error_flag
      logical, parameter                  :: debug_flag = .false.
      character(len=13)     , parameter   :: efmt       = '(a,1x,es12.5)'
      character(len=9)      , parameter   :: ifmt       = '(a,1x,i5)'
      character(len=9)      , parameter   :: lfmt       = '(a,1x,l1)'
      !----- External functions. ----------------------------------------------------------!
      logical               , external    :: isnan_real
      !------------------------------------------------------------------------------------!


      !-- Point to the cohort structures --------------------------------------------------!
      cpatch => csite%patch(ipa)
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !      Decide whether or not to solve dynamic plant hydraulics.                      !
      !------------------------------------------------------------------------------------!
      select case (plant_hydro_scheme)
      case (0)
         !------ Compatible with original ED-2.2, do not track plant hydraulics. ----------!
         do ico = 1, cpatch%ncohorts
             ipft = cpatch%pft(ico)

             cpatch%wflux_wl        (ico)    = 0.
             cpatch%wflux_gw        (ico)    = 0.
             cpatch%wflux_gw_layer(:,ico)    = 0.

             cpatch%leaf_rwc        (ico)    = 1.0
             cpatch%wood_rwc        (ico)    = 1.0
             cpatch%leaf_psi        (ico)    = 0.
             cpatch%wood_psi        (ico)    = 0.
         end do
        !----------------------------------------------------------------------------------!
      case default
         !---------------------------------------------------------------------------------!
         !    Dynamic plant hydraulics.                                                    !
         !---------------------------------------------------------------------------------!

         !---------------------------------------------------------------------------------!
         !     Calculate water potential and conductance in each soil layer in preparation
         ! for later calculations.
         !---------------------------------------------------------------------------------!
         do k = 1,nzg
            nsoil = ntext_soil(k)

            !------------------------------------------------------------------------------!
            !      Get bounded soil moisture.                                              !
            !  MLO.  The lower bound used to be air-dry soil moisture.  This causes issues !
            !  in the RK4 integrator if the soil moisture is just slightly above air-dry   !
            !  and dtlsm is long.  For the time being, I am assuming that soil             !
            !  conductivity is halted just below the permanent wilting point.  Similarly,  !
            !  I am assuming that matric potential cannot exceed a value slightly less     !
            !  than the bubbling point.                                                    !
            !------------------------------------------------------------------------------!
            swater_min = mg_safe * soil(nsoil)%soilcp  + om_safe * soil(nsoil)%soilwp
            swater_max = mg_safe * soil(nsoil)%sfldcap + om_safe * soil(nsoil)%slmsts
            swater_use = max( swater_min                                                   &
                            , min(swater_max                                               &
                                 ,csite%soil_water(k,ipa) * csite%soil_fracliq(k,ipa) ) )
            !------------------------------------------------------------------------------!


            !----- Clapp & Hornberger curves. ---------------------------------------------!
            soil_psi(k)  = matric_potential(nsoil,swater_use)
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !    In the model, soil can't get drier than residual soil moisture.  Ensure   !
            ! that hydraulic conductivity is effectively zero in case soil moisture        !
            ! reaches this level or drier.                                                 !
            !------------------------------------------------------------------------------!
            if (csite%soil_water(k,ipa) < swater_min) then
               soil_cond(k) = 0.
            else
               soil_cond(k) = wdns * hydr_conduct(k,nsoil,csite%soil_water(k,ipa)          &
                                                 ,csite%soil_fracliq(k,ipa))
            end if
            !------------------------------------------------------------------------------!
         end do
         !---------------------------------------------------------------------------------!




         !---------------------------------------------------------------------------------!
         !      Loop over cohorts, calculate plant hydraulic fluxes.                       !
         !---------------------------------------------------------------------------------!
         cohortloop: do ico = 1, cpatch%ncohorts
            ipft = cpatch%pft(ico)

            !------------------------------------------------------------------------------!
            !     Track the plant hydraulics when either leaf or wood are resolvable. Leaf !
            ! become un-resolvable in dry season when all leaves are shed. In this scena-  !
            ! rio, we still nedd to track plant hydraulics to update wood_psi so that the  !
            ! model knows when to reflush the leaves. Otherwise, soil water will never re- !
            ! fill wood water pool.                                                        !
            !     Special case: leaf is not resolvable when bleaf is on allometry while    !
            ! wood can be still resolvable. Transpiration is always zero, whereas soil     !
            ! water can still flow into the stem. This will ultimately leads to unreal-    !
            ! istically high wood_psi and even positive psi due to numerical erros. There- !
            ! fore, we do not track plant hydraulics in this case.                         !
            !------------------------------------------------------------------------------!

            track_hydraulics = cpatch%leaf_resolvable(ico) .or.                            &
                              (cpatch%wood_resolvable(ico) .and.                           & 
                               .not. (cpatch%elongf(ico) == 1.0 .and.                      &
                                      .not. cpatch%leaf_resolvable(ico) ))


            if (track_hydraulics) then
               !----- Prepare input for plant water flux calculations. --------------------!
               sap_frac    = dbh2sf(cpatch%dbh(ico),ipft)                    ! m2
               sap_area    = sap_frac * pio4 * (cpatch%dbh(ico) / 100.) ** 2 ! m2
               bsap        = ( cpatch%bdeada   (ico) + cpatch%bdeadb   (ico)               &
                             + cpatch%bsapwooda(ico) + cpatch%bsapwoodb(ico) ) * sap_frac
               transp      = ( cpatch%fs_open(ico) * cpatch%psi_open(ico)                  &
                             + (1. - cpatch%fs_open(ico)) * cpatch%psi_closed(ico) )       &
                           * cpatch%lai(ico) / cpatch%nplant(ico)            ! kg / s
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !     Note  that the current leaf_water_int has already deducted losses     !
               ! through transpiration (but not sapflow).  Consequently, leaf psi can be   !
               ! very low.  To get meaningful leaf_psi, leaf_psi and leaf_water_int        !
               ! become temporarily decoupled: transpiration is added back to              !
               ! leaf_water_int.  Therefore, leaf_psi represents the water potential at    !
               ! the START of the timestep.                                                !
               !---------------------------------------------------------------------------!
               call rwc2psi(cpatch%leaf_rwc(ico),cpatch%wood_rwc(ico),ipft                 &
                           ,cpatch%leaf_psi(ico),cpatch%wood_psi(ico))
               c_leaf = leaf_water_cap(ipft) * C2B * cpatch%bleaf(ico)
               if (c_leaf > 0.) then
                  cpatch%leaf_psi(ico) = cpatch%leaf_psi(ico)  & ! m
                                       + transp * dtlsm        & ! kgH2O
                                       / c_leaf                ! ! kgH2O/m
                  !------------------------------------------------------------------------!
               else
                  !----- No leaves, set leaf_psi the same as wood_psi - hite. -------------!
                  cpatch%leaf_psi(ico) = cpatch%wood_psi(ico) - cpatch%hite(ico)
                  !------------------------------------------------------------------------!
               end if
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !      Run sanity check.  The code will crash if any of these happen.       !
               !                                                                           !
               ! 1.  If leaf_psi is invalid (run the debugger, the problem may be else-    !
               !     where)                                                                !
               ! 2.  If leaf_psi is positive (non-sensical)                                !
               ! 3.  If leaf_psi is too negative (also non-sensical)                       !
               !---------------------------------------------------------------------------!
               error_flag(1) = isnan_real(cpatch%leaf_psi(ico)) ! NaN values
               error_flag(2) = cpatch%leaf_psi(ico) > 0.        ! Positive potential
               error_flag(3) = merge( cpatch%leaf_psi(ico) < small_psi_min(ipft)           &
                                    , cpatch%leaf_psi(ico) < leaf_psi_min (ipft)           &
                                    , cpatch%is_small(ico)                        )
               if ((debug_flag .and. (dco == 0 .or. ico == dco)) .or. any(error_flag)) then
                  write (unit=*,fmt='(a)') ' '
                  write (unit=*,fmt='(92a)') ('=',k=1,92)
                  write (unit=*,fmt='(92a)') ('=',k=1,92)
                  write (unit=*,fmt='(a)'  )                                               &
                     ' Invalid leaf_psi detected.'
                  write (unit=*,fmt='(92a)') ('-',k=1,92)
                  write (unit=*,fmt='(a,i4.4,2(1x,i2.2),1x,f6.0)') ' TIME           : '    &
                                                  ,current_time%year,current_time%month    &
                                                  ,current_time%date,current_time%time
                  write (unit=*,fmt='(a)'  ) ' '
                  write (unit=*,fmt=ifmt   ) ' + IPA              =',ipa
                  write (unit=*,fmt=ifmt   ) ' + ICO              =',ico
                  write (unit=*,fmt=ifmt   ) ' + PFT              =',ipft
                  write (unit=*,fmt=ifmt   ) ' + KRDEPTH          =',cpatch%krdepth(ico)
                  write (unit=*,fmt=efmt   ) ' + HEIGHT           =',cpatch%hite(ico)
                  write (unit=*,fmt=lfmt   ) ' + SMALL            =',cpatch%is_small(ico)

                  write (unit=*,fmt='(a)'  ) ' '
                  write (unit=*,fmt=lfmt   ) ' + FINITE           =',.not. error_flag(1)
                  write (unit=*,fmt=lfmt   ) ' + NEGATIVE         =',.not. error_flag(2)
                  write (unit=*,fmt=lfmt   ) ' + BOUNDED          =',.not. error_flag(3)

                  write (unit=*,fmt='(a)'  ) ' '
                  write (unit=*,fmt=efmt   ) ' + LEAF_PSI_MIN     =',leaf_psi_min (ipft)
                  write (unit=*,fmt=efmt   ) ' + SMALL_PSI_MIN    =',small_psi_min(ipft)

                  write (unit=*,fmt='(a)'  ) ' '
                  write (unit=*,fmt=efmt   ) ' + BLEAF            =',cpatch%bleaf(ico)
                  write (unit=*,fmt=efmt   ) ' + LAI              =',cpatch%lai(ico) 
                  write (unit=*,fmt=efmt   ) ' + NPLANT           =',cpatch%nplant(ico) 
                  write (unit=*,fmt=efmt   ) ' + BSAPWOOD (Hydro) =',bsap 
                  write (unit=*,fmt=efmt   ) ' + BSAPWOOD (Allom) ='                       &
                                                                  , cpatch%bsapwooda(ico)  &
                                                                  + cpatch%bsapwoodb(ico)

                  write (unit=*,fmt=efmt   ) ' + BROOT            =',cpatch%broot(ico)
                  write (unit=*,fmt=efmt   ) ' + SAPWOOD_AREA     =',sap_area
  
                  write (unit=*,fmt='(a)'  ) ' '
                  write (unit=*,fmt=efmt   ) ' + TRANSP           =',transp
                  write (unit=*,fmt=efmt   ) ' + C_LEAF           =',c_leaf
                  write (unit=*,fmt=efmt   ) ' + PSI_OPEN         =',cpatch%psi_open(ico)
                  write (unit=*,fmt=efmt   ) ' + PSI_CLOSED       =',cpatch%psi_closed(ico)
                  write (unit=*,fmt=efmt   ) ' + FS_OPEN          =',cpatch%fs_open(ico)
                  write (unit=*,fmt='(a)'  ) ' '
                  write (unit=*,fmt=efmt   ) ' + LEAF_PSI         =',cpatch%leaf_psi(ico)
                  write (unit=*,fmt=efmt   ) ' + LEAF_RWC         =',cpatch%leaf_rwc(ico)
                  write (unit=*,fmt=efmt   ) ' + LEAF_WATER_INT   ='                       &
                                                               ,cpatch%leaf_water_int(ico)
                  write (unit=*,fmt=efmt   ) ' + LEAF_WATER_IM2   ='                       &
                                                               ,cpatch%leaf_water_im2(ico)
                  write (unit=*,fmt='(a)'  ) ' '
                  write (unit=*,fmt=efmt   ) ' + WOOD_PSI         =',cpatch%wood_psi(ico)
                  write (unit=*,fmt=efmt   ) ' + WOOD_RWC         =',cpatch%wood_rwc(ico)
                  write (unit=*,fmt=efmt   ) ' + WOOD_WATER_INT   ='                       &
                                                               ,cpatch%wood_water_int(ico)
                  write (unit=*,fmt=efmt   ) ' + WOOD_WATER_IM2   ='                       &
                                                               ,cpatch%wood_water_im2(ico)
                  write (unit=*,fmt='(a)'  ) ' '
                  write (unit=*,fmt=efmt   ) ' + WFLUX_GW (LAST)  =',cpatch%wflux_gw(ico) 
                  write (unit=*,fmt=efmt   ) ' + WFLUX_WL (LAST)  =',cpatch%wflux_wl(ico)


                  write (unit=*,fmt='(a)'        ) ' '
                  write (unit=*,fmt='(92a)'      ) ('-',k=1,92)
                  write (unit=*,fmt='(a,2(1x,a))') '    K','    SOIL_PSI','WFLUX_GW_LYR'
                  write (unit=*,fmt='(92a)'      ) ('-',k=1,92)
                  do k = 1, nzg
                     write (unit=*,fmt='(i5,2(1x,es12.5))')                                &
                                                k,soil_psi(k),cpatch%wflux_gw_layer(k,ico)
                  end do
                  write (unit=*,fmt='(92a)'   ) ('-',k=1,92)
                  write (unit=*,fmt='(a)'     ) ' '
                  write (unit=*,fmt='(92a)'   ) ('=',k=1,92)
                  write (unit=*,fmt='(92a)'   ) ('=',k=1,92)
                  write (unit=*,fmt='(a)'     ) ' '

                  if (any(error_flag)) then 
                     call fatal_error('Plant Hydrodynamics is off-track.'                  &
                                     ,'plant_hydro_driver','plant_hydro.f90')
                  end if
                  !------------------------------------------------------------------------!
               end if
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !    Find water fluxes.  Note that transp is from last timestep's psi_open  !
               ! and psi_closed.                                                           !
               !---------------------------------------------------------------------------!
               call calc_plant_water_flux(                            &
                        dtlsm                                            & ! input
                       ,sap_area,cpatch%nplant(ico),ipft                 & ! input
                       ,cpatch%is_small(ico),cpatch%krdepth(ico)         & ! input
                       ,cpatch%bleaf(ico),bsap,cpatch%broot(ico)         & ! input
                       ,cpatch%hite(ico),cpatch%root_frac(:,ico)         & ! input
                       ,transp,cpatch%leaf_psi(ico),cpatch%wood_psi(ico) & ! input
                       ,soil_psi,soil_cond,ipa,ico                       & ! input
                       ,cpatch%wflux_wl(ico),cpatch%wflux_gw(ico)        & ! output
                       ,cpatch%wflux_gw_layer(:,ico))                    ! ! output
               !---------------------------------------------------------------------------!
            else
               !----- Neither leaves nor wood are resolvable.  Assume zero flow. ----------!
               cpatch%wflux_wl(ico) = 0.
               cpatch%wflux_gw(ico) = 0.
               cpatch%wflux_gw_layer(:,ico)  = 0.
               !---------------------------------------------------------------------------!
            end if
            !------------------------------------------------------------------------------!
         end do cohortloop
         !---------------------------------------------------------------------------------!
      end select
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Update the most frequent timescale averages.                                   !
      !------------------------------------------------------------------------------------!
      do ico = 1, cpatch%ncohorts
         cpatch%fmean_leaf_psi      (ico) = cpatch%fmean_leaf_psi      (ico)               &
                                          + cpatch%leaf_psi            (ico)               &
                                          * dtlsm_o_frqsum
         cpatch%fmean_wood_psi      (ico) = cpatch%fmean_wood_psi      (ico)               &
                                          + cpatch%wood_psi            (ico)               &
                                          * dtlsm_o_frqsum
         cpatch%fmean_leaf_water_int(ico) = cpatch%fmean_leaf_water_int(ico)               &
                                          + cpatch%leaf_water_int      (ico)               &
                                          * dtlsm_o_frqsum
         cpatch%fmean_wood_water_int(ico) = cpatch%fmean_wood_water_int(ico)               &
                                          + cpatch%wood_water_int      (ico)               &
                                          * dtlsm_o_frqsum
         if (cpatch%dmax_leaf_psi(ico) == 0.) then
             cpatch%dmax_leaf_psi(ico) =  cpatch%leaf_psi(ico)
         else
             cpatch%dmax_leaf_psi(ico) =  max( cpatch%dmax_leaf_psi(ico)                   &
                                             , cpatch%leaf_psi     (ico) )
         end if
         if (cpatch%dmin_leaf_psi(ico) == 0.) then
             cpatch%dmin_leaf_psi(ico) =  cpatch%leaf_psi(ico)
         else
             cpatch%dmin_leaf_psi(ico) =  min( cpatch%dmin_leaf_psi(ico)                   &
                                             , cpatch%leaf_psi     (ico) )
         end if
         if (cpatch%dmax_wood_psi(ico) == 0.) then
             cpatch%dmax_wood_psi(ico) =  cpatch%wood_psi(ico)
         else
             cpatch%dmax_wood_psi(ico) =  max( cpatch%dmax_wood_psi(ico)                   &
                                             , cpatch%wood_psi     (ico) )
         end if
         if (cpatch%dmin_wood_psi(ico) == 0.) then
             cpatch%dmin_wood_psi(ico) =  cpatch%wood_psi(ico)
         else
             cpatch%dmin_wood_psi(ico) =  min( cpatch%dmin_wood_psi(ico)                   &
                                             , cpatch%wood_psi     (ico) )
         end if
         !---------------------------------------------------------------------------------!
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
   !> [X16] Xu X, Medvigy D, Powers JS, Becknell JM , Guan K. 2016. Diversity in plant 
   !>       hydraulic traits explains seasonal and inter-annual variations of vegetation
   !>       dynamics in seasonally dry tropical forests. New Phytol. 212: 80-95. 
   !>       doi:10.1111/nph.14009.
   !>
   !> [K03] Katul G, Leuning R , Oren R. 2003. Relationship between plant hydraulic and
   !>       biochemical properties derived from a steady-state coupled water and carbon 
   !>       transport model. Plant Cell Environ. 26: 339-350. 
   !>       doi:10.1046/j.1365-3040.2003.00965.x.
   !>
   !> \author Xiangtao Xu, 29 Jan. 2018
   !---------------------------------------------------------------------------------------!
   subroutine calc_plant_water_flux(dt                                  & !timestep
               ,sap_area,nplant,ipft,is_small,krdepth                   & !plant input
               ,bleaf,bsap,broot,hite ,root_frac                        & !plant input
               ,transp,leaf_psi,wood_psi                                & !plant input
               ,soil_psi,soil_cond                                      & !soil  input
               ,ipa,ico                                                 & !debug input
               ,wflux_wl,wflux_gw,wflux_gw_layer)                       ! !flux  output
      use soil_coms       , only : dslz8                ! ! intent(in)
      use grid_coms       , only : nzg                  ! ! intent(in)
      use consts_coms     , only : pi18                 & ! intent(in)
                                 , lnexp_min8           ! ! intent(in)
      use rk4_coms        , only : tiny_offset          ! ! intent(in)
      use pft_coms        , only : leaf_water_cap       & ! intent(in) 
                                 , wood_water_cap       & ! intent(in)
                                 , leaf_psi_min         & ! intent(in)
                                 , wood_psi_min         & ! intent(in)
                                 , small_psi_min        & ! intent(in)
                                 , wood_psi50           & ! intent(in)
                                 , wood_Kmax            & ! intent(in)
                                 , wood_Kexp            & ! intent(in)
                                 , vessel_curl_factor   & ! intent(in)
                                 , SRA                  & ! intent(in)
                                 , C2B                  ! ! intent(in)
      use ed_misc_coms    , only : current_time         ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real   ,                 intent(in)  :: dt             !time step           [      s]
      real   ,                 intent(in)  :: sap_area       !sapwood_area        [     m2]
      real   ,                 intent(in)  :: nplant         !plant density       [  pl/m2]
      integer,                 intent(in)  :: ipft           !plant funct. type   [    ---]
      integer,                 intent(in)  :: krdepth        !Max. rooting depth  [    ---]
      logical,                 intent(in)  :: is_small       !Small cohort?       [    T|F]
      real   ,                 intent(in)  :: bleaf          !leaf biomass        [    kgC]
      real   ,                 intent(in)  :: bsap           !sapwood biomass     [ kgC/pl]
      real   ,                 intent(in)  :: broot          !fine root biomass   [ kgC/pl]
      real   ,                 intent(in)  :: hite           !plant height        [      m]
      real   , dimension(nzg), intent(in)  :: root_frac      !Root fraction       [      m]
      real   ,                 intent(in)  :: transp         !transpiration       [   kg/s]
      real   ,                 intent(in)  :: leaf_psi       !leaf water pot.     [      m]
      real   ,                 intent(in)  :: wood_psi       !wood water pot.     [      m]
      real   , dimension(nzg), intent(in)  :: soil_psi       !soil water pot.     [      m]
      real   , dimension(nzg), intent(in)  :: soil_cond      !soil water cond.    [kg/m2/s]
      integer,                 intent(in)  :: ipa            !Patch index         [    ---]
      integer,                 intent(in)  :: ico            !Cohort index        [    ---]
      real   ,                 intent(out) :: wflux_wl       !wood-leaf flux      [   kg/s]
      real   ,                 intent(out) :: wflux_gw       !ground-wood flux    [   kg/s]
      real   , dimension(nzg), intent(out) :: wflux_gw_layer !wflux_gw for each soil layer
      !----- Temporary double precision variables (input/output). -------------------------!
      real(kind=8)                 :: dt_d
      real(kind=8)                 :: sap_area_d
      real(kind=8)                 :: bleaf_d
      real(kind=8)                 :: bsap_d
      real(kind=8)                 :: broot_d
      real(kind=8)                 :: nplant_d
      real(kind=8)                 :: hite_d
      real(kind=8)                 :: transp_d
      real(kind=8)                 :: leaf_psi_d
      real(kind=8)                 :: wood_psi_d
      real(kind=8), dimension(nzg) :: soil_psi_d
      real(kind=8), dimension(nzg) :: soil_cond_d
      real(kind=8)                 :: wflux_wl_d
      real(kind=8)                 :: wflux_gw_d
      real(kind=8), dimension(nzg) :: wflux_gw_layer_d
      !----- Temporary double precision variables (PFT parameters). -----------------------!
      real(kind=8)                 :: leaf_psi_min_d
      real(kind=8)                 :: wood_psi_min_d
      real(kind=8)                 :: leaf_psi_lwr_d
      real(kind=8)                 :: wood_psi_lwr_d
      real(kind=8)                 :: SRA_d
      real(kind=8)                 :: wood_psi50_d
      real(kind=8)                 :: wood_Kexp_d
      real(kind=8)                 :: wood_Kmax_d
      real(kind=8)                 :: vessel_curl_factor_d
      real(kind=8)                 :: root_frac_d          !fraction of roots
      !----- Auxiliary variables. ---------------------------------------------------------!
      real(kind=8)                          :: exp_term             !exponent term
      real(kind=8)                          :: ap                   ![s-1]
      real(kind=8)                          :: bp                   ![m s-1]
      real(kind=8)                          :: stem_cond            !stem conductance
      real(kind=8)                          :: plc                  !plant loss of conductance
      real(kind=8)                          :: c_leaf               !leaf water capacitance
      real(kind=8)                          :: c_stem               !stem water capacitance
      real(kind=8)                          :: RAI                  !root area index
      real(kind=8)                          :: proj_leaf_psi        !projected leaf water pot.
      real(kind=8)                          :: proj_wood_psi        !projected wood water pot. 
      real(kind=8)                          :: gw_cond              !g->w water conductivity
      real(kind=8)                          :: org_wood_psi         !used for small tree
      real(kind=8)                          :: org_leaf_psi         !used for small tree
      real(kind=8)                          :: weighted_soil_psi
      real(kind=8)                          :: weighted_gw_cond
      real(kind=8)                          :: total_water_supply
      real(kind=8)      , dimension(nzg)    :: layer_water_supply
      !----- Counters. --------------------------------------------------------------------!
      integer                               :: k
      !----- Boolean flags. ---------------------------------------------------------------!
      logical                               :: zero_flow_wl
      logical                               :: zero_flow_gw
      logical           , dimension(5)      :: error_flag
      !----- Local constants. -------------------------------------------------------------!
      character(len=13) , parameter         :: efmt       = '(a,1x,es12.5)'
      character(len=9)  , parameter         :: ifmt       = '(a,1x,i5)'
      character(len=9)  , parameter         :: lfmt       = '(a,1x,l1)'
      integer           , parameter         :: dco        = 0
      logical           , parameter         :: debug_flag = .false.
      !----- External function ------------------------------------------------------------!
      real(kind=4)      , external          :: sngloff       ! Safe dble 2 single precision
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !    Convert all input state vars and some PFT-dependent parameters to double 
      ! precision.
      !------------------------------------------------------------------------------------!
      dt_d                 = dble(dt                      )
      sap_area_d           = dble(sap_area                )
      bleaf_d              = dble(bleaf                   )
      bsap_d               = dble(bsap                    )
      broot_d              = dble(broot                   )
      nplant_d             = dble(nplant                  )
      hite_d               = dble(hite                    )
      transp_d             = dble(transp                  )
      leaf_psi_d           = dble(leaf_psi                )
      wood_psi_d           = dble(wood_psi                )
      soil_psi_d           = dble(soil_psi                )
      soil_cond_d          = dble(soil_cond               )
      SRA_d                = dble(SRA               (ipft))
      !----- Minimum threshold depends on whether the plant is small or large. ------------!
      if (is_small) then
         leaf_psi_min_d    = dble(small_psi_min     (ipft))
         wood_psi_min_d    = dble(small_psi_min     (ipft))
      else
         leaf_psi_min_d    = dble(leaf_psi_min      (ipft))
         wood_psi_min_d    = dble(wood_psi_min      (ipft))
      end if 
      leaf_psi_lwr_d       = om_buff_d * leaf_psi_min_d
      wood_psi_lwr_d       = om_buff_d * wood_psi_min_d
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Update plant hydrodynamics.
      !
      !     The regular solver assumes stem water pool is way larger than the leaf water
      ! pool. In cases where leaf water pool is of similar magnitude to stem water pool 
      ! (seedlings and grasses), leaf water potential is forced to be the same as stem 
      ! water potential, to maintain numerical stability.  This, however, may bias the 
      ! water potential estimates for these plants.
      !
      ! Water flow is calculated from canopy to roots
      ! Positive flux means upward flow (g->w, w->l, l->air)
      !------------------------------------------------------------------------------------!


      !----- Initialise proj_psi as the starting psi. Also save the initial psi values. ---!
      proj_leaf_psi = leaf_psi_d 
      proj_wood_psi = wood_psi_d
      org_wood_psi  = wood_psi_d
      org_leaf_psi  = leaf_psi_d
      !------------------------------------------------------------------------------------!



      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      ! 1. Calculate wood/stem/root to leaf water flow
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      ! First, check the relative magnitude of leaf and sapwood water pool
      ! If it is a small tree/grass, force psi_leaf to be the same as psi_stem
      !------------------------------------------------------------------------------------!
      c_leaf = dble(leaf_water_cap(ipft) * C2B) * bleaf_d            ! kg H2O / m
      c_stem = dble(wood_water_cap(ipft) * C2B) * (broot_d + bsap_d) ! kg H2O / m
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     If cohort is considered small (see stable_cohorts.f90 for flag), then we do
      ! not distinguish between leaves and wood when calculating hydrodynamics.
      !------------------------------------------------------------------------------------!
      if (is_small) then
         !---------------------------------------------------------------------------------!
         !   1.1.  Small tree, force leaf_psi to be the same as wood_psi.  Calculate the 
         !         new veg_psi of mixing leaf and wood
         !---------------------------------------------------------------------------------!
         wood_psi_d   = (c_leaf * leaf_psi_d + c_stem * wood_psi_d) / (c_leaf+c_stem)
         leaf_psi_d   = wood_psi_d
         !---------------------------------------------------------------------------------!


         !----- Use c_leaf+c_stem as the total capacitance to solve stem_psi later. -------!
         c_stem = c_stem + c_leaf
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !    In this case, we temporarily assign transpiration as wflux_wl_d since leaves 
         ! and wood are treated as a single entity.  The value will be recalculated once 
         ! we obtain the projected water potential.
         !---------------------------------------------------------------------------------!
         wflux_wl_d   = transp_d
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !    Set zero_flow_wl to .false., so it appears correctly in the error message in !
         ! case the model crashes.                                                         !
         !---------------------------------------------------------------------------------!
         zero_flow_wl = .false.
         !---------------------------------------------------------------------------------!

      else
         !---------------------------------------------------------------------------------!
         ! 1.2.  Regular case, big trees.
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !    Special cases in which flow between leaves and wood should be zero.  Perhaps
         ! there are better ways to systematically avoid these traps, but currently, this
         ! is done in a case-by-case manner.
         !
         ! Case 0.  Negative wflux_wl overchange wood storage?  This should never occur
         !          if we initialise leaf water potential slightly lower than the stem
         !          water potential.
         !
         ! Case 1.  Cohort has no leaves.
         !
         ! Case 2.  Both wood and leaves are very dry, and either (a) wood cannot support 
         !          upward sapflow or (b) leaf cannot support downward flow.  This could 
         !          happen for dying trees experiencing extreme drought. 
         !
         ! Case 3.  The cohort just grows out of 'small tree status'.  Their leaves can be 
         !          over-charged with water because gravitational effect was not 
         !          considered for leaf water potential of small trees.  As a result, this 
         !          can lead to a down-ward sapflow, and potentially over-charging the 
         !          sapwood. We need to zero the flow in this case as well, until 
         !          leaf_psi_d drops below wood_psi_d - hite_d.
         !---------------------------------------------------------------------------------!
         zero_flow_wl = ( c_leaf == 0.d0                            ) .or.  & ! Case 1
                        ( leaf_psi_d >= (wood_psi_d - hite_d) .and.         &
                          leaf_psi_d <= leaf_psi_lwr_d              ) .or.  & ! Case 2a
                        ( leaf_psi_d <= (wood_psi_d - hite_d) .and.         &
                          wood_psi_d <= wood_psi_lwr_d              ) .or.  & ! Case 2b
                        ( leaf_psi_d >  (wood_psi_d - hite_d)       )       ! ! Case 3
         !---------------------------------------------------------------------------------!




         !---------------------------------------------------------------------------------!
         !    Decide whether or not to calculate sapflow.
         !---------------------------------------------------------------------------------!
         if (zero_flow_wl) then
            !------------------------------------------------------------------------------!
            ! 1.2.1. No need to calculate sapflow
            !------------------------------------------------------------------------------!
            wflux_wl_d = 0.d0

            !------ Proj_leaf_psi is only dependent upon transpiration. -------------------!
            if (c_leaf > 0.) then
                proj_leaf_psi = leaf_psi_d - transp_d * dt_d / c_leaf
            else
                proj_leaf_psi = leaf_psi_d
            end if
            !------------------------------------------------------------------------------!

         else
            !------------------------------------------------------------------------------!
            !     We do need to calculate sapflow.  First convert some PFT-dependent 
            ! parameters to double precision.
            !------------------------------------------------------------------------------!
            wood_psi50_d         = dble(wood_psi50        (ipft))
            wood_Kexp_d          = dble(wood_Kexp         (ipft))
            wood_Kmax_d          = dble(wood_Kmax         (ipft))
            vessel_curl_factor_d = dble(vessel_curl_factor(ipft))
            !------------------------------------------------------------------------------!

            !----- Calculate plant loss of conductivity [dimensionless]. ------------------!
            plc = 1.d0 / (1.d0 + (wood_psi_d / wood_psi50_d) ** wood_Kexp_d)
            !------------------------------------------------------------------------------!



            !----- Calculate stem conductance [kg / s]. -----------------------------------!
            stem_cond = wood_Kmax_d * plc                 & ! kg/m/s
                      * sap_area_d                        & ! conducting area m2
                      / (hite_d * vessel_curl_factor_d)   ! ! conducting length m
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !     Find sapflow.
            !------------------------------------------------------------------------------!
            if (stem_cond == 0.) then
               !---- 1.2.2. Zero flux because stem conductivity is also zero. -------------!
               wflux_wl_d = 0.d0
               !---------------------------------------------------------------------------!
            else
               !---------------------------------------------------------------------------!
               ! 1.2.3. "Normal case", with positive c_leaf and positive stem_cond.  Check
               !        reference X16 for derivation of the equations.
               !---------------------------------------------------------------------------!
               ap = - stem_cond / c_leaf                                            ! [1/s]
               bp = ((wood_psi_d - hite_d) * stem_cond - transp_d) / c_leaf         ! [m/s]

               !----- Project the final leaf psi. -----------------------------------------!
               exp_term      = exp(max(ap * dt_d,lnexp_min8))
               proj_leaf_psi = max( leaf_psi_lwr_d                                         &
                                  , ((ap * leaf_psi_d + bp) * exp_term - bp) / ap )
               !---------------------------------------------------------------------------!


               !----- Calculate the average sapflow rate within the time step [kgH2O/s]. --!
               wflux_wl_d = (proj_leaf_psi - leaf_psi_d) * c_leaf / dt_d + transp_d
               !---------------------------------------------------------------------------!
            end if
            !------------------------------------------------------------------------------!
         end if
         !---------------------------------------------------------------------------------!
      end if
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!


      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      ! 2.  Calculate ground -> wood/stem/root water flow
      !------------------------------------------------------------------------------------!
      weighted_soil_psi  = 0.d0
      weighted_gw_cond   = 0.d0
      layer_water_supply = 0.d0
      total_water_supply = 0.d0

      !----- Loop over all soil layers to get the aggregated water conductance. -----------!
      do k = krdepth,nzg

         !----- Retrieve the root fraction of this layer. ---------------------------------!
         root_frac_d = dble(root_frac(k))
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !  Calculate RAI in each layer.                                                   !
         !---------------------------------------------------------------------------------!
         RAI = broot_d * SRA_d * root_frac_d * nplant_d  ! m2/m2
         !---------------------------------------------------------------------------------!

         !---------------------------------------------------------------------------------!
         !    Calculate soil-root water conductance kg H2O/m/s based on reference [K03].
         !---------------------------------------------------------------------------------!
         gw_cond = soil_cond_d(k) * sqrt(RAI) / (pi18 * dslz8(k))  & ! kg H2O / m3 / s
                 / nplant_d                                        ! ! conducting area  m2
         !---------------------------------------------------------------------------------!




         !---------------------------------------------------------------------------------!
         !      Disable hydraulic redistribution.  Assume roots will shut down if they are 
         ! going to lose water to soil.
         !---------------------------------------------------------------------------------!
         if (soil_psi_d(k) <= wood_psi_d) then
            gw_cond = 0.d0
         end if
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !    Calculate weighted conductance, weighted psi, and water_supply_layer_frac.
         !---------------------------------------------------------------------------------!
         weighted_gw_cond      = weighted_gw_cond + gw_cond                  ! kgH2O/m/s
         weighted_soil_psi     = weighted_soil_psi + gw_cond * soil_psi_d(k) ! kgH2O/s
         layer_water_supply(k) = gw_cond * (soil_psi_d(k) - wood_psi_d)      ! kgH2O/s
         !---------------------------------------------------------------------------------!
      end do
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !    Now we can calculate ground->wood water flow.
      ! First we handle special cases
      !------------------------------------------------------------------------------------!
      zero_flow_gw = (c_stem           == 0.d0) .or.  & ! No sapwood or fine  roots
                     (weighted_gw_cond == 0.d0)       ! ! soil is drier than wood

      if (zero_flow_gw) then
         !---------------------------------------------------------------------------------!
         !     No need to calculate water flow: wood psi is only dependent upon sapflow.
         !---------------------------------------------------------------------------------!
         wflux_gw_d    = 0.d0
         if (c_stem > 0.) then
            !----- Make sure that projected wood psi will be bounded. ---------------------!
            wflux_wl_d    = min(wflux_wl_d, (wood_psi_d - wood_psi_lwr_d) * c_stem / dt_d )
            proj_wood_psi = wood_psi_d - wflux_wl_d * dt_d / c_stem
            !------------------------------------------------------------------------------!
         else
            proj_wood_psi = wood_psi_d
         end if
         !---------------------------------------------------------------------------------!
      else
         !---------------------------------------------------------------------------------!
         !     Calculate the average soil water uptake. Check reference X16 for derivation
         ! of the equations.
         !---------------------------------------------------------------------------------!
         ap = - weighted_gw_cond  / c_stem  ! ! 1/s
         bp = (weighted_soil_psi - wflux_wl_d) / c_stem ! m/s
         !---------------------------------------------------------------------------------!

         !----- Project the final wood psi, but ensure it will be bounded. ----------------!
         exp_term        = exp(max(ap * dt_d,lnexp_min8))
         proj_wood_psi   = max( wood_psi_lwr_d                                             &
                              , ((ap * wood_psi_d + bp) * exp_term - bp) / ap )
         !---------------------------------------------------------------------------------!


         !----- Calculate the average root extraction within the time step [kgH2O/s]. -----!
         wflux_gw_d     = (proj_wood_psi - wood_psi_d) * c_stem  / dt_d + wflux_wl_d
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !     Re-calculate the water fluxes in the cases of small cohorts.
      !------------------------------------------------------------------------------------!
      if (is_small) then
         !---------------------------------------------------------------------------------!
         !     Ground->wood flux (wflux_gw_d) is correct, no need to update.  However, we
         ! do need to update wood->leaf (wflux_wl_d).
         !---------------------------------------------------------------------------------!
         proj_leaf_psi = proj_wood_psi
         wflux_wl_d    = (proj_leaf_psi - org_leaf_psi)  * c_leaf / dt_d + transp_d
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Now estimate the water uptake from each layer based on layer_water_supply.
      !------------------------------------------------------------------------------------!
      if (sum(layer_water_supply) == 0.d0) then
         wflux_gw_layer_d = 0.d0
      else
         wflux_gw_layer_d = layer_water_supply / sum(layer_water_supply) * wflux_gw_d
      end if
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!




      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      ! 3.  Sanity check.  Stop the simulation in case anything went wrong.   Sympotms for
      !     things going wrong:
      !     a.  NaN values --- Run the debugger
      !     b.  Projected leaf/wood potential is positive
      !     c.  Current leaf/wood potential is positive
      !     d.  Projected leaf/wood potential is less than minimum acceptable
      !     e.  Current leaf/wood potential is less than minimum acceptable
      !------------------------------------------------------------------------------------!
      error_flag(1) = isnan(wflux_wl_d)              .or. isnan(wflux_gw_d)
      error_flag(2) = proj_leaf_psi > 0.             .or. proj_wood_psi > 0.
      error_flag(3) = leaf_psi_d    > 0.             .or. wood_psi_d    > 0.
      error_flag(4) = proj_leaf_psi < leaf_psi_min_d .or. proj_wood_psi < wood_psi_min_d
      error_flag(5) = leaf_psi_d    < leaf_psi_min_d .or. wood_psi_d    < wood_psi_min_d

      if ( (debug_flag .and. (dco == 0 .or. ico == dco)) .or. any(error_flag)) then
         write (unit=*,fmt='(a)') ' '
         write (unit=*,fmt='(92a)') ('=',k=1,92)
         write (unit=*,fmt='(92a)') ('=',k=1,92)
         write (unit=*,fmt='(a)'  ) ' Plant hydrodynamics inconsistency detected!!'
         write (unit=*,fmt='(92a)') ('-',k=1,92)
         write (unit=*,fmt='(a,i4.4,2(1x,i2.2),1x,f6.0)') ' TIME           : '             &
                                                     ,current_time%year,current_time%month &
                                                     ,current_time%date,current_time%time
         write (unit=*,fmt='(a)'  ) ' '
         write (unit=*,fmt=ifmt   ) ' + IPA              =',ipa
         write (unit=*,fmt=ifmt   ) ' + ICO              =',ico
         write (unit=*,fmt=ifmt   ) ' + PFT              =',ipft
         write (unit=*,fmt=ifmt   ) ' + KRDEPTH          =',krdepth
         write (unit=*,fmt=efmt   ) ' + HEIGHT           =',hite

         write (unit=*,fmt='(a)'  ) ' '
         write (unit=*,fmt=lfmt   ) ' + IS_SMALL         =',is_small
         write (unit=*,fmt=lfmt   ) ' + ZERO_FLOW_WL     =',zero_flow_wl
         write (unit=*,fmt=lfmt   ) ' + ZERO_FLOW_GW     =',zero_flow_gw

         write (unit=*,fmt='(a)'  ) ' '
         write (unit=*,fmt=efmt   ) ' + BLEAF            =',bleaf
         write (unit=*,fmt=efmt   ) ' + BSAPWOOD         =',bsap
         write (unit=*,fmt=efmt   ) ' + BROOT            =',broot
         write (unit=*,fmt=efmt   ) ' + SAPWOOD_AREA     =',sap_area

         write (unit=*,fmt='(a)'  ) ' '
         write (unit=*,fmt=lfmt   ) ' + Finite fluxes     =',.not. error_flag(1)
         write (unit=*,fmt=lfmt   ) ' + Negative Proj Psi =',.not. error_flag(2)
         write (unit=*,fmt=lfmt   ) ' + Negative Curr Psi =',.not. error_flag(3)
         write (unit=*,fmt=lfmt   ) ' + Bounded Proj Psi  =',.not. error_flag(4)
         write (unit=*,fmt=lfmt   ) ' + Bounded Curr Psi  =',.not. error_flag(5)

         write (unit=*,fmt='(a)'  ) ' '
         write (unit=*,fmt=efmt   ) ' + LEAF_PSI_MIN      =',leaf_psi_min (ipft)
         write (unit=*,fmt=efmt   ) ' + SMALL_PSI_MIN     =',small_psi_min(ipft)

         write (unit=*,fmt='(a)'  ) ' '
         write (unit=*,fmt=efmt   ) ' + TRANSP           =',transp
         write (unit=*,fmt=efmt   ) ' + LEAF_PSI (INPUT) =',leaf_psi
         write (unit=*,fmt=efmt   ) ' + WOOD_PSI (INPUT) =',wood_psi
         write (unit=*,fmt=efmt   ) ' + LEAF_PSI (PROJ.) =',proj_leaf_psi
         write (unit=*,fmt=efmt   ) ' + WOOD_PSI (PROJ.) =',proj_wood_psi
         write (unit=*,fmt=efmt   ) ' + WFLUX_GW         =',wflux_gw_d
         write (unit=*,fmt=efmt   ) ' + WFLUX_WL         =',wflux_wl_d


         write (unit=*,fmt='(a)'        ) ' '
         write (unit=*,fmt='(92a)'      ) ('-',k=1,92)
         write (unit=*,fmt='(a,2(1x,a))') '    K','    SOIL_PSI','WFLUX_GW_LYR'
         write (unit=*,fmt='(92a)'      ) ('-',k=1,92)
         do k = 1, nzg
            write (unit=*,fmt='(i5,2(1x,es12.5))') k,soil_psi(k),wflux_gw_layer_d(k)
         end do
         write (unit=*,fmt='(92a)') ('=',k=1,92)
         write (unit=*,fmt='(92a)') ('=',k=1,92)
         write (unit=*,fmt='(a)'  ) ' '

         if (any(error_flag)) then 
            call fatal_error('Plant Hydrodynamics is off-track.'                           &
                            ,'calc_plant_water_flux','plant_hydro.f90')
         end if
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!
 


      !------------------------------------------------------------------------------------!
      !     Copy all the results to output variables.
      !------------------------------------------------------------------------------------!
      wflux_wl = sngloff(wflux_wl_d,tiny_offset)
      do k = 1, nzg
         wflux_gw_layer(k) = sngloff(wflux_gw_layer_d(k),tiny_offset)
      end do
      wflux_gw = sum(wflux_gw_layer)
      !------------------------------------------------------------------------------------!


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
   subroutine tw2rwc(leaf_water_int,wood_water_int,is_small,bleaf,bsapwooda,bsapwoodb      &
                    ,bdeada,bdeadb,broot,dbh,ipft,leaf_rwc,wood_rwc)
      use pft_coms       , only : leaf_water_sat     & ! intent(in)
                                , wood_water_sat     & ! intent(in)
                                , leaf_rwc_min       & ! intent(in)
                                , wood_rwc_min       & ! intent(in)
                                , small_rwc_min      & ! intent(in)
                                , C2B                ! ! intent(in)
      use allometry      , only : dbh2sf             ! ! function
      use physiology_coms, only : plant_hydro_scheme ! ! intent(in)
      use consts_coms    , only : tiny_num           ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real      , intent(in)    ::  leaf_water_int ! Total internal water of leaf      [ kg]
      real      , intent(in)    ::  wood_water_int ! Total internal water of wood      [ kg]
      logical   , intent(in)    ::  is_small       ! Small/large plant flag            [T|F]
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
      elseif (is_small) then
         leaf_rwc  = op_buff * small_rwc_min(ipft)
      else
         leaf_rwc  = op_buff * leaf_rwc_min(ipft)
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
      elseif (is_small) then
         wood_rwc = op_buff * small_rwc_min(ipft)
      elseif (is_small) then
         wood_rwc = op_buff * wood_rwc_min(ipft)
      end if
      !------------------------------------------------------------------------------------!

      return
   end subroutine tw2rwc
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !  SUBROUTINE: PSI2TW            
   !> \brief Convert water potential to total water for both leaf and wood
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
   subroutine tw2psi(leaf_water_int,wood_water_int,is_small,bleaf,bsapwooda,bsapwoodb      &
                    ,bdeada,bdeadb,broot,dbh,ipft,leaf_psi,wood_psi)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real      , intent(in)    ::  leaf_water_int ! Total internal water of leaf      [ kg]
      real      , intent(in)    ::  wood_water_int ! Total internal water of wood      [ kg]
      logical   , intent(in)    ::  is_small       ! Small/large plant flag            [T|F]
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
      call tw2rwc(leaf_water_int,wood_water_int,is_small,bleaf,bsapwooda,bsapwoodb,bdeada  &
                 ,bdeadb,broot,dbh,ipft,leaf_rwc,wood_rwc)
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
      leaf_water_im2 = nplant * leaf_water_int
      wood_water_im2 = nplant * wood_water_int
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
   subroutine twe2twi(leaf_water_im2,wood_water_im2,nplant,leaf_water_int,wood_water_int)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real      , intent(in)  ::  leaf_water_im2 ! Extensive leaf internal water   [ kg/m2]
      real      , intent(in)  ::  wood_water_im2 ! Water potential of wood         [ kg/m2]
      real      , intent(in)  ::  nplant         ! Stem density                    [ pl/m2]
      real      , intent(out) ::  leaf_water_int ! Total internal water of leaf    [ kg/pl]
      real      , intent(out) ::  wood_water_int ! Total internal water of wood    [ kg/pl]
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Convert water from intensive to extensive.                                    !
      !------------------------------------------------------------------------------------!
      leaf_water_int = leaf_water_im2 / nplant
      wood_water_int = wood_water_im2 / nplant
      !------------------------------------------------------------------------------------!

      return
   end subroutine twe2twi
   !=======================================================================================!
   !=======================================================================================!





   !=======================================================================================!
   !=======================================================================================!
   !  SUBROUTINE: PSI2TWE
   !> \breif Convert water potential of leaf and wood to extensive water storage
   !> \details This sub-routine is useful when we need to go from water potential to 
   !>          storage, but we don't need the intermediate quantities (relative and
   !>          individual water contents).
   !=======================================================================================!
   subroutine psi2twe(leaf_psi,wood_psi,ipft,nplant,bleaf,bsapwooda,bsapwoodb,bdeada       &
                     ,bdeadb,broot,dbh,leaf_water_im2,wood_water_im2)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real      , intent(in)    ::  leaf_psi       ! Water potential of leaves     [     m]
      real      , intent(in)    ::  wood_psi       ! Water potential of wood       [     m]
      integer   , intent(in)    ::  ipft           ! Plant functional type         [     -]
      real      , intent(in)    ::  nplant         ! Stem density                  [ pl/m2]
      real      , intent(in)    ::  bleaf          ! Biomass of leaf               [kgC/pl]
      real      , intent(in)    ::  bsapwooda      ! Aboveground sapwood biomass   [kgC/pl]
      real      , intent(in)    ::  bsapwoodb      ! Belowground sapwood biomass   [kgC/pl]
      real      , intent(in)    ::  bdeada         ! Aboveground heartwood biomass [kgC/pl]
      real      , intent(in)    ::  bdeadb         ! Belowground heartwood biomass [kgC/pl]
      real      , intent(in)    ::  broot          ! Biomass of fine root          [kgC/pl]
      real      , intent(in)    ::  dbh            ! Diameter at breast height     [    cm]
      real      , intent(out)   ::  leaf_water_im2 ! Extensive leaf internal water [ kg/m2]
      real      , intent(out)   ::  wood_water_im2 ! Extensive wood internal water [ kg/m2]
      !----- Local variables. -------------------------------------------------------------!
      real                      ::  leaf_rwc       ! Relative leaf water content   [    --]
      real                      ::  wood_rwc       ! Relative wood water content   [    --]
      real                      ::  leaf_water_int ! Intensive leaf internal water [ kg/pl]
      real                      ::  wood_water_int ! Intensive wood internal water [ kg/pl]
      !------------------------------------------------------------------------------------!

      !----- 1. Potential -> relative water content. --------------------------------------!
      call psi2rwc(leaf_psi,wood_psi,ipft,leaf_rwc,wood_rwc)
      !----- 2. Relative water content -> Intensive internal water. -----------------------!
      call rwc2tw(leaf_rwc,wood_rwc,bleaf,bsapwooda,bsapwoodb,bdeada,bdeadb,broot,dbh,ipft &
                 ,leaf_water_int,wood_water_int)
      !----- 3. Intensive internal water -> Extensive internal water. ---------------------!
      call twi2twe(leaf_water_int,wood_water_int,nplant,leaf_water_im2,wood_water_im2)
      !------------------------------------------------------------------------------------!


      return
   end subroutine psi2twe
   !=======================================================================================!
   !=======================================================================================!

   !=======================================================================================!
   !  SUBROUTINE: UPDATE_PLC
   !> \breif update percentage loss of xylem conductance using daily minimum leaf psi
   !> \details This subroutine is called at daily time scale in growth_balive.f90
   !> Daily minimum leaf psi is used because upper branch, which has similar
   !> water potential as leaf, should be the most vulnerable section along the hydraulic
   !> pathway
   !=======================================================================================!
   subroutine update_plc(cpatch,ico)
      use ed_state_vars,  only : patchtype             ! ! structure
     use ed_misc_coms   , only : current_time          & ! intent(in)
                               , simtime               ! ! structure
      use physiology_coms,only : plant_hydro_scheme    ! ! intent(in)
      use pft_coms,       only : wood_psi50            & ! intent(in)
                               , wood_Kexp             ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(patchtype), target       :: cpatch
      integer        , intent(in)   :: ico
      !----- Locals    --------------------------------------------------------------------!
      real                          :: plc_today
      integer                       :: ipft
      real                          :: ndaysi
      type(simtime)                 :: lastmonth
      !------------------------------------------------------------------------------------!
      
      ! No need to update PLC if we are not tracking hydro-dynamics
      if (plant_hydro_scheme == 0) return

      call lastmonthdate(current_time,lastmonth,ndaysi)
      ipft = cpatch%pft(ico)
      
      plc_today  =  max(0., 1. - 1. /                                         &
                        (1. + (cpatch%dmin_leaf_psi(ico)                      &
                        / wood_psi50(ipft)) ** wood_Kexp(ipft)))
      cpatch%plc_monthly   (13,ico) = cpatch%plc_monthly   (13,ico) + plc_today * ndaysi

    return

   end subroutine update_plc
   !=======================================================================================!
   !=======================================================================================!



end module plant_hydro

!==========================================================================================!
!==========================================================================================!
