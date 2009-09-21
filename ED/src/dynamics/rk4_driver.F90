!==========================================================================================!
!==========================================================================================!
!     This module contains the wrappers for the Runge-Kutta integration scheme.            !
!==========================================================================================!
!==========================================================================================!
module rk4_driver

   contains
   !=======================================================================================!
   !=======================================================================================!
   !      Main driver of short-time scale dynamics of the Runge-Kutta integrator           !
   !      for the land surface model.                                                      !
   !---------------------------------------------------------------------------------------!
   subroutine rk4_timestep(cgrid,ifm)
      use rk4_coms               , only : integration_vars     & ! structure
                                        , rk4patchtype         & ! structure
                                        , zero_rk4_patch       & ! subroutine
                                        , zero_rk4_cohort      & ! subroutine
                                        , integration_buff     & ! intent(out)
                                        , rk4met               ! ! intent(out)
      use ed_state_vars          , only : edtype               & ! structure
                                        , polygontype          & ! structure
                                        , sitetype             & ! structure
                                        , patchtype            ! ! structure
      use grid_coms              , only : nzg                  ! ! intent(in)
      use canopy_struct_dynamics , only : canopy_turbulence8   ! ! subroutine
      use ed_misc_coms           , only : current_time         ! ! intent(in)
      implicit none

      !----------- Use MPI timing calls, need declarations --------------------------------!
      include 'mpif.h'
      !----- Arguments --------------------------------------------------------------------!
      type(edtype)              , target      :: cgrid
      integer                   , intent (in) :: ifm
      !----- Local variables --------------------------------------------------------------!
      type(polygontype)         , pointer     :: cpoly
      type(sitetype)            , pointer     :: csite
      type(patchtype)           , pointer     :: cpatch
      integer                                 :: ipy,isi,ipa
      integer                                 :: iun
      integer, dimension(nzg)                 :: ed_ktrans
      real                                    :: sum_lai_rbi
      real                                    :: wcurr_loss2atm
      real                                    :: ecurr_loss2atm
      real                                    :: co2curr_loss2atm
      real                                    :: wcurr_loss2drainage
      real                                    :: ecurr_loss2drainage
      real                                    :: wcurr_loss2runoff
      real                                    :: ecurr_loss2runoff
      real                                    :: ecurr_latent
      !----- Variables declared differently depending on the user's compilation options. --!
#if USE_MPIWTIME
      real(kind=8)                            :: time_py_start
      real(kind=8)                            :: time_py_end
#else
      real                                    :: time_py_start
      real                                    :: time_py_spent
#endif
      !----- Locally saved variables. -----------------------------------------------------!
      logical                   , save        :: first_time=.true.
      !----- Local constants. -------------------------------------------------------------!
      logical                   , parameter   :: print_fields = .false.
      !----- Functions --------------------------------------------------------------------!
      real                      , external    :: walltime
      !------------------------------------------------------------------------------------!
      
      polygonloop: do ipy = 1,cgrid%npolygons
         cpoly => cgrid%polygon(ipy)

#if USE_MPIWTIME
         time_py_start = MPI_Wtime() 
#else
         time_py_start = walltime(0.)
#endif

         siteloop: do isi = 1,cpoly%nsites
            csite => cpoly%site(isi)

            patchloop: do ipa = 1,csite%npatches
               cpatch => csite%patch(ipa)

               if (print_fields) then
                  iun = 30+ipa
                  if (first_time) then
                     write (unit=iun,fmt='(19(a,1x))') 'YEAR','MONTH','  DAY','  TIME'     &
                          ,'  VEG.HEIGHT','   CAN.DEPTH','       SPEED','       PRESS'     &
                          ,'    ATM.TEMP','    CAN.TEMP','   SFCW.TEMP','    SOIL.TMP'     &
                          ,'    ATM.SHV ','     CAN.SHV','   SFCW.MASS','  SOIL.MOIST'     &
                          ,'    ATM.CO2 ','     CAN.CO2','    CAN.DENS'
                  end if
                  write (unit=iun,fmt='(i4.4,2(4x,i2.2),1x,f6.0,15(1x,es12.5))')           &
                     current_time%year,current_time%month,current_time%date                &
                    ,current_time%time,csite%veg_height(ipa),csite%can_depth(ipa)          &
                    ,cpoly%met(isi)%vels,cpoly%met(isi)%prss,cpoly%met(isi)%atm_tmp        &
                    ,csite%can_temp(ipa),csite%sfcwater_tempk(1,ipa)                       &
                    ,csite%soil_tempk(nzg,ipa),cpoly%met(isi)%atm_shv,csite%can_shv(ipa)   &
                    ,csite%sfcwater_mass(1,ipa),csite%soil_water(nzg,ipa)                  &
                    ,cpoly%met(isi)%atm_co2,csite%can_co2(ipa),csite%can_rhos(ipa)
                   
               end if

               
               !----- Reset all buffers to zero, as a safety measure. ---------------------!
               call zero_rk4_patch(integration_buff%initp)
               call zero_rk4_patch(integration_buff%dinitp)
               call zero_rk4_patch(integration_buff%yscal)
               call zero_rk4_patch(integration_buff%y)
               call zero_rk4_patch(integration_buff%dydx)
               call zero_rk4_patch(integration_buff%yerr)
               call zero_rk4_patch(integration_buff%ytemp)
               call zero_rk4_patch(integration_buff%ak2)
               call zero_rk4_patch(integration_buff%ak3)
               call zero_rk4_patch(integration_buff%ak4)
               call zero_rk4_patch(integration_buff%ak5)
               call zero_rk4_patch(integration_buff%ak6)
               call zero_rk4_patch(integration_buff%ak7)
               call zero_rk4_cohort(integration_buff%initp)
               call zero_rk4_cohort(integration_buff%dinitp)
               call zero_rk4_cohort(integration_buff%yscal)
               call zero_rk4_cohort(integration_buff%y)
               call zero_rk4_cohort(integration_buff%dydx)
               call zero_rk4_cohort(integration_buff%yerr)
               call zero_rk4_cohort(integration_buff%ytemp)
               call zero_rk4_cohort(integration_buff%ak2)
               call zero_rk4_cohort(integration_buff%ak3)
               call zero_rk4_cohort(integration_buff%ak4)
               call zero_rk4_cohort(integration_buff%ak5)
               call zero_rk4_cohort(integration_buff%ak6)
               call zero_rk4_cohort(integration_buff%ak7)

               !----- Get velocity for aerodynamic resistance. ----------------------------!
               if (csite%can_theta(ipa) < cpoly%met(isi)%atm_theta) then
                  cpoly%met(isi)%vels = cpoly%met(isi)%vels_stab
               else
                  cpoly%met(isi)%vels = cpoly%met(isi)%vels_unstab
               end if

               !---------------------------------------------------------------------------!
               !    Copy the meteorological variables to the rk4met structure.             !
               !---------------------------------------------------------------------------!
               call copy_met_2_rk4met(cpoly%met(isi)%vels,cpoly%met(isi)%atm_enthalpy      &
                                     ,cpoly%met(isi)%atm_theta,cpoly%met(isi)%atm_tmp      &
                                     ,cpoly%met(isi)%atm_shv,cpoly%met(isi)%atm_co2        &
                                     ,cpoly%met(isi)%geoht,cpoly%met(isi)%exner            &
                                     ,cpoly%met(isi)%pcpg,cpoly%met(isi)%qpcpg             &
                                     ,cpoly%met(isi)%dpcpg,cpoly%met(isi)%prss             &
                                     ,cpoly%met(isi)%geoht,cpoly%lsl(isi)                  &
                                     ,cgrid%lon(ipy),cgrid%lat(ipy))


               !---------------------------------------------------------------------------!
               !     Set up the integration patch.                                         !
               !---------------------------------------------------------------------------!
               call copy_patch_init(csite,ipa,integration_buff%initp)



               !---------------------------------------------------------------------------!
               !     Calculate the canopy geometry, and the scalar transport coefficients. !
               !---------------------------------------------------------------------------!
               call canopy_turbulence8(csite,integration_buff%initp,isi,ipa,.true.)



               !----- Get photosynthesis, stomatal conductance, and transpiration. --------!
               call canopy_photosynthesis(csite,ipa,cpoly%met(isi)%vels                    &
                          ,cpoly%met(isi)%atm_tmp,cpoly%met(isi)%prss,ed_ktrans            &
                          ,csite%ntext_soil(:,ipa),csite%soil_water(:,ipa)                 &
                          ,csite%soil_fracliq(:,ipa),cpoly%lsl(isi),sum_lai_rbi            &
                          ,cpoly%leaf_aging_factor(:,isi),cpoly%green_leaf_factor(:,isi) )

               !----- Compute root and heterotrophic respiration. -------------------------!
               call soil_respiration(csite,ipa)

               !---------------------------------------------------------------------------!
               !    This is the driver for the integration process...                      !
               !---------------------------------------------------------------------------!
               call integrate_patch(csite,integration_buff%initp,ipa,isi,ipy,ifm           &
                                   ,wcurr_loss2atm,ecurr_loss2atm,co2curr_loss2atm         &
                                   ,wcurr_loss2drainage,ecurr_loss2drainage                &
                                   ,wcurr_loss2runoff,ecurr_loss2runoff,ecurr_latent)

               !---------------------------------------------------------------------------!
               !    Update the minimum monthly temperature, based on canopy temperature.   !
               !---------------------------------------------------------------------------!
               if (cpoly%site(isi)%can_temp(ipa) < cpoly%min_monthly_temp(isi)) then
                  cpoly%min_monthly_temp(isi) = cpoly%site(isi)%can_temp(ipa)
               end if
               
               !-------------------------------------------------------------!
               !     Compute the residuals.                                  !
               !-------------------------------------------------------------!
               call compute_budget(csite,cpoly%lsl(isi)                      &
                                  ,cpoly%met(isi)%pcpg,cpoly%met(isi)%qpcpg  &
                                  ,ipa,wcurr_loss2atm,ecurr_loss2atm         &
                                  ,co2curr_loss2atm,wcurr_loss2drainage      &
                                  ,ecurr_loss2drainage,wcurr_loss2runoff     &
                                  ,ecurr_loss2runoff,ecurr_latent            &
                                  ,cpoly%area(isi),cgrid%cbudget_nep(ipy))
               

            end do patchloop
         end do siteloop

#if USE_MPIWTIME
         time_py_end            = MPI_Wtime() 
         cgrid%walltime_py(ipy) = cgrid%walltime_py(ipy) + (time_py_end-time_py_start)
#else
         ! You will get funky results unless you use MPI_Wtime, best to flag as
         ! nonesense so the analysis is not misleading
!         time_py_spent          = walltime(time_py_start)
         cgrid%walltime_py(ipy) = cgrid%walltime_py(ipy) -9.9 !+ dble(time_py_spent)
#endif

      end do polygonloop

      first_time = .false.
      return
   end subroutine rk4_timestep
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine will drive the integration process.                               !
   !---------------------------------------------------------------------------------------!
   subroutine integrate_patch(csite,initp,ipa,isi,ipy,ifm,wcurr_loss2atm,ecurr_loss2atm    &
                             ,co2curr_loss2atm,wcurr_loss2drainage,ecurr_loss2drainage     &
                             ,wcurr_loss2runoff,ecurr_loss2runoff,ecurr_latent)
      use ed_state_vars   , only : sitetype             & ! structure
                                 , patchtype            ! ! structure
      use ed_misc_coms    , only : dtlsm                ! ! intent(in)
      use soil_coms       , only : soil_rough           & ! intent(in)
                                 , snow_rough           ! ! intent(in)
      use canopy_air_coms , only : exar8                ! ! intent(in)
      use consts_coms     , only : vonk8                & ! intent(in)
                                 , cp8                  & ! intent(in)
                                 , cpi8                 ! ! intent(in)
      use rk4_coms        , only : integration_vars     & ! structure
                                 , rk4patchtype         & ! structure
                                 , rk4met               & ! intent(inout)
                                 , zero_rk4_patch       & ! subroutine
                                 , zero_rk4_cohort      & ! subroutine
                                 , tbeg                 & ! intent(inout)
                                 , tend                 & ! intent(inout)
                                 , dtrk4                & ! intent(inout)
                                 , dtrk4i               & ! intent(inout)
                                 , ibranch_thermo       & ! intent(in)
                                 , effarea_water        & ! intent(out)
                                 , effarea_heat         ! ! intent(out)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(sitetype)        , target      :: csite   
      type(rk4patchtype)    , target      :: initp
      integer               , intent(in)  :: ifm     
      integer               , intent(in)  :: ipy     
      integer               , intent(in)  :: isi     
      integer               , intent(in)  :: ipa     
      real                  , intent(out) :: wcurr_loss2atm
      real                  , intent(out) :: ecurr_loss2atm
      real                  , intent(out) :: co2curr_loss2atm
      real                  , intent(out) :: wcurr_loss2drainage
      real                  , intent(out) :: ecurr_loss2drainage
      real                  , intent(out) :: wcurr_loss2runoff
      real                  , intent(out) :: ecurr_loss2runoff
      real                  , intent(out) :: ecurr_latent
      !----- Local variables --------------------------------------------------------------!
      real(kind=8)                          :: hbeg
      !----- Locally saved variable -------------------------------------------------------!
      logical                  , save       :: first_time=.true.
      !------------------------------------------------------------------------------------!

      !----- Assigning some constants which will remain the same throughout the run. ------!
      if (first_time) then
         first_time = .false.
         tbeg   = 0.d0
         tend   = dble(dtlsm)
         dtrk4  = tend - tbeg
         dtrk4i = 1.d0/dtrk4
         
         !---------------------------------------------------------------------------------!
         !    The area factor for heat and water exchange between canopy and vegetation is !
         ! applied only on LAI, and it depends on how we are considering the branches and  !
         ! twigs.  If their area isn't explicitly defined, we add a 0.2 factor to the      !
         ! area because WPA will be 0.  Otherwise, we don't add anything to the LAI, and   !
         ! let WPA to do the job.                                                          !
         !---------------------------------------------------------------------------------!
         select case (ibranch_thermo)
         case (0)
            effarea_water = 1.2d0
            effarea_heat  = 2.2d0
         case default
            effarea_water = 1.0d0
            effarea_heat  = 2.0d0
         end select
      end if

      !------------------------------------------------------------------------------------!
      !      Initial step size.  Experience has shown that giving this too large a value   !
      ! causes the integrator to fail (e.g., soil layers become supersaturated).           !
      !------------------------------------------------------------------------------------!
      hbeg = dble(csite%htry(ipa))


      !------------------------------------------------------------------------------------!
      !     Zero the canopy-atmosphere flux values.  These values are updated every dtlsm, !
      ! so they must be zeroed at each call.                                               !
      !------------------------------------------------------------------------------------!
      initp%upwp = 0.d0
      initp%tpwp = 0.d0
      initp%qpwp = 0.d0
      initp%cpwp = 0.d0
      initp%wpwp = 0.d0

      !----- Go into the ODE integrator. --------------------------------------------------!
      call odeint(hbeg,csite,ipa,isi,ipy,ifm)

      !------------------------------------------------------------------------------------!
      !      Normalize canopy-atmosphere flux values.  These values are updated every      !
      ! dtlsm, so they must be normalized every time.                                      !
      !------------------------------------------------------------------------------------!
      initp%upwp = initp%can_rhos * initp%upwp * dtrk4i
      initp%tpwp = initp%can_rhos * initp%tpwp * dtrk4i
      initp%qpwp = initp%can_rhos * initp%qpwp * dtrk4i
      initp%cpwp = initp%can_rhos * initp%cpwp * dtrk4i
      initp%wpwp = initp%can_rhos * initp%wpwp * dtrk4i
      
      
      !------------------------------------------------------------------------------------!
      ! Move the state variables from the integrated patch to the model patch.             !
      !------------------------------------------------------------------------------------!
      call initp2modelp(tend-tbeg,initp,csite,ipa,isi,ipy,wcurr_loss2atm,ecurr_loss2atm    &
                       ,co2curr_loss2atm,wcurr_loss2drainage,ecurr_loss2drainage           &
                       ,wcurr_loss2runoff,ecurr_loss2runoff,ecurr_latent)

      return
   end subroutine integrate_patch
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine will copy the variables from the integration buffer to the state  !
   ! patch and cohorts.                                                                    !
   !---------------------------------------------------------------------------------------!
   subroutine initp2modelp(hdid,initp,csite,ipa,isi,ipy,wbudget_loss2atm,ebudget_loss2atm  &
                          ,co2budget_loss2atm,wbudget_loss2drainage,ebudget_loss2drainage  &
                          ,wbudget_loss2runoff,ebudget_loss2runoff,ebudget_latent)
      use rk4_coms             , only : rk4patchtype         & ! structure
                                      , rk4met               & ! intent(in)
                                      , rk4min_veg_temp      & ! intent(in)
                                      , rk4max_veg_temp      & ! intent(in)
                                      , tiny_offset          & ! intent(in) 
                                      , checkbudget          ! ! intent(in)
      use ed_state_vars        , only : sitetype             & ! structure
                                      , patchtype            & ! structure
                                      , edgrid_g             ! ! structure
      use consts_coms          , only : day_sec              & ! intent(in)
                                      , t3ple8               ! ! intent(in)
      use ed_misc_coms         , only : fast_diagnostics     & ! intent(in)
                                      , dtlsm                ! ! intent(in)
      use soil_coms            , only : soil8                & ! intent(in)
                                      , slz8                 ! ! intent(in)
      use grid_coms            , only : nzg                  & ! intent(in)
                                      , nzs                  ! ! intent(in)
      use therm_lib            , only : qwtk                 ! ! subroutine
      use ed_therm_lib         , only : ed_grndvap           ! ! subroutine
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(rk4patchtype), target      :: initp
      type(sitetype)    , target      :: csite
      real(kind=8)      , intent(in)  :: hdid
      integer           , intent(in)  :: ipa
      integer           , intent(in)  :: ipy
      integer           , intent(in)  :: isi
      real              , intent(out) :: wbudget_loss2atm
      real              , intent(out) :: ebudget_loss2atm
      real              , intent(out) :: co2budget_loss2atm
      real              , intent(out) :: wbudget_loss2drainage
      real              , intent(out) :: ebudget_loss2drainage
      real              , intent(out) :: wbudget_loss2runoff
      real              , intent(out) :: ebudget_loss2runoff
      real              , intent(out) :: ebudget_latent
      !----- Local variables --------------------------------------------------------------!
      type(patchtype)   , pointer     :: cpatch
      integer                         :: mould
      integer                         :: ico
      integer                         :: k
      integer                         :: kclosest
      integer                         :: ksn
      integer                         :: nsoil
      integer                         :: nlsw1
      real(kind=8)                    :: available_water
      real(kind=8)                    :: tmp_energy
      real                            :: surface_temp
      real                            :: surface_fliq
      !----- Local contants ---------------------------------------------------------------!
      real        , parameter         :: tendays_sec=10.*day_sec
      !----- External function ------------------------------------------------------------!
      real        , external          :: sngloff
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Most variables require just a simple copy.  More comments will be made next to !
      ! those in which this is not true.  All floating point variables are converted back  !
      ! to single precision.                                                               !
      !------------------------------------------------------------------------------------!
      csite%can_enthalpy(ipa) = sngloff(initp%can_enthalpy  ,tiny_offset)
      csite%can_theta(ipa)    = sngloff(initp%can_theta     ,tiny_offset)
      csite%can_prss(ipa)     = sngloff(initp%can_prss      ,tiny_offset)
      csite%can_temp(ipa)     = sngloff(initp%can_temp      ,tiny_offset)
      csite%can_shv(ipa)      = sngloff(initp%can_shv       ,tiny_offset)
      csite%can_co2(ipa)      = sngloff(initp%can_co2       ,tiny_offset)
      csite%can_rhos(ipa)     = sngloff(initp%can_rhos      ,tiny_offset)
      csite%can_depth(ipa)    = sngloff(initp%can_depth     ,tiny_offset)

      csite%ustar(ipa)    = sngloff(initp%ustar   ,tiny_offset)
      csite%tstar(ipa)    = sngloff(initp%tstar   ,tiny_offset)
      csite%qstar(ipa)    = sngloff(initp%qstar   ,tiny_offset)
      csite%cstar(ipa)    = sngloff(initp%cstar   ,tiny_offset)

      csite%upwp(ipa)     = sngloff(initp%upwp    ,tiny_offset)
      csite%wpwp(ipa)     = sngloff(initp%wpwp    ,tiny_offset)
      csite%tpwp(ipa)     = sngloff(initp%tpwp    ,tiny_offset)
      csite%qpwp(ipa)     = sngloff(initp%qpwp    ,tiny_offset)
      csite%cpwp(ipa)     = sngloff(initp%cpwp    ,tiny_offset)

      !------------------------------------------------------------------------------------!
      !    These variables are fast scale fluxes, and they may not be allocated, so just   !
      ! check this before copying.                                                         !
      !------------------------------------------------------------------------------------!
      if(fast_diagnostics) then

         csite%avg_vapor_vc(ipa)         =sngloff(initp%avg_vapor_vc      ,tiny_offset)
         csite%avg_dew_cg(ipa)           =sngloff(initp%avg_dew_cg        ,tiny_offset)
         csite%avg_vapor_gc(ipa)         =sngloff(initp%avg_vapor_gc      ,tiny_offset)
         csite%avg_wshed_vg(ipa)         =sngloff(initp%avg_wshed_vg      ,tiny_offset)
         csite%avg_vapor_ac(ipa)         =sngloff(initp%avg_vapor_ac      ,tiny_offset)
         csite%avg_transp(ipa)           =sngloff(initp%avg_transp        ,tiny_offset)
         csite%avg_evap(ipa)             =sngloff(initp%avg_evap          ,tiny_offset)
         csite%avg_drainage(ipa)         =sngloff(initp%avg_drainage      ,tiny_offset)
         csite%avg_netrad(ipa)           =sngloff(initp%avg_netrad        ,tiny_offset)
         csite%avg_sensible_vc(ipa)      =sngloff(initp%avg_sensible_vc   ,tiny_offset)
         csite%avg_qwshed_vg(ipa)        =sngloff(initp%avg_qwshed_vg     ,tiny_offset)
         csite%avg_sensible_gc(ipa)      =sngloff(initp%avg_sensible_gc   ,tiny_offset)
         csite%avg_sensible_ac(ipa)      =sngloff(initp%avg_sensible_ac   ,tiny_offset)
         csite%avg_carbon_ac(ipa)        =sngloff(initp%avg_carbon_ac     ,tiny_offset)
         do k = rk4met%lsl, nzg
            csite%avg_sensible_gg(k,ipa) =sngloff(initp%avg_sensible_gg(k)   ,tiny_offset)
            csite%avg_smoist_gg(k,ipa)   =sngloff(initp%avg_smoist_gg(k)     ,tiny_offset)
            csite%avg_smoist_gc(k,ipa)   =sngloff(initp%avg_smoist_gc(k)     ,tiny_offset)
         end do
      end if

      if(checkbudget) then
         co2budget_loss2atm    = sngloff(initp%co2budget_loss2atm   ,tiny_offset)
         ebudget_loss2atm      = sngloff(initp%ebudget_loss2atm     ,tiny_offset)
         ebudget_loss2drainage = sngloff(initp%ebudget_loss2drainage,tiny_offset)
         ebudget_loss2runoff   = sngloff(initp%ebudget_loss2runoff  ,tiny_offset)
         ebudget_latent        = sngloff(initp%ebudget_latent       ,tiny_offset)
         wbudget_loss2atm      = sngloff(initp%wbudget_loss2atm     ,tiny_offset)
         wbudget_loss2drainage = sngloff(initp%wbudget_loss2drainage,tiny_offset)
         wbudget_loss2runoff   = sngloff(initp%wbudget_loss2runoff  ,tiny_offset)
      else
         co2budget_loss2atm             = 0.
         ebudget_loss2atm               = 0.
         ebudget_loss2drainage          = 0.
         ebudget_loss2runoff            = 0.
         ebudget_latent                 = 0.
         wbudget_loss2atm               = 0.
         wbudget_loss2drainage          = 0.
         wbudget_loss2runoff            = 0.
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     The following is not a pure diagnostic, it is used for phenology and mortality !
      ! functions, preserve this variable and its dependencies in all contexts.            !
      !------------------------------------------------------------------------------------!
      csite%avg_daily_temp(ipa) = csite%avg_daily_temp(ipa) + csite%can_temp(ipa)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      ! paw_avg - 10-day average of plant available water.                              !
      !     I don't think this is currently used, but it may be turned on for drought-     !
      ! -related  phenology.                                                               !
      !------------------------------------------------------------------------------------!
      cpatch => csite%patch(ipa)
      do ico = 1,cpatch%ncohorts
         available_water = 0.d0
         do k = cpatch%krdepth(ico), nzg - 1
            nsoil = csite%ntext_soil(k,ipa)
            available_water = available_water                                              &
                            + (initp%soil_water(k) - soil8(nsoil)%soilcp)                  &
                            * (slz8(k+1)-slz8(k))                                          &
                            / (soil8(nsoil)%slmsts - soil8(nsoil)%soilcp)
         end do
         nsoil = csite%ntext_soil(nzg,ipa)
         available_water = available_water                                                 &
                         + (initp%soil_water(nzg) - soil8(nsoil)%soilcp)                   &
                         * (-1.d0*slz8(nzg))                                               &
                         / (soil8(nsoil)%slmsts -soil8(nsoil)%soilcp) 
         available_water = available_water / (-1.d0*slz8(cpatch%krdepth(ico)))


         cpatch%paw_avg(ico) = cpatch%paw_avg(ico)*(1.0-sngl(hdid)/tendays_sec)            &
                             + sngl(available_water)*sngl(hdid)/tendays_sec
      end do

      
      do k = rk4met%lsl, nzg
         csite%soil_water(k,ipa)   = sngloff(initp%soil_water(k)  ,tiny_offset)
         csite%soil_energy(k,ipa)  = sngloff(initp%soil_energy(k) ,tiny_offset)
         csite%soil_tempk(k,ipa)   = sngloff(initp%soil_tempk(k)  ,tiny_offset)
         csite%soil_fracliq(k,ipa) = sngloff(initp%soil_fracliq(k),tiny_offset)
      end do
      

      !------------------------------------------------------------------------------------!
      !    Surface water energy is computed in J/m² inside the integrator. Converting it   !
      ! back to J/kg in the layers that surface water/snow still exists.                   !
      !------------------------------------------------------------------------------------!
      csite%nlev_sfcwater(ipa) = initp%nlev_sfcwater
      do k = 1, csite%nlev_sfcwater(ipa)
         csite%sfcwater_depth(k,ipa)   = sngloff(initp%sfcwater_depth(k)   ,tiny_offset)
         csite%sfcwater_mass(k,ipa)    = sngloff(initp%sfcwater_mass(k)    ,tiny_offset)
         csite%sfcwater_tempk(k,ipa)   = sngloff(initp%sfcwater_tempk(k)   ,tiny_offset)
         csite%sfcwater_fracliq(k,ipa) = sngloff(initp%sfcwater_fracliq(k) ,tiny_offset)
         tmp_energy                    = initp%sfcwater_energy(k)/initp%sfcwater_mass(k)
         csite%sfcwater_energy(k,ipa)  = sngloff(tmp_energy                ,tiny_offset)
      end do
      !------------------------------------------------------------------------------------!
      !    For the layers that no longer exist, assign zeroes for prognostic variables,    !
      ! and something for temperature and liquid fraction (just to avoid singularities,    !
      ! and funny numbers in the output, but these values are meaningless and should never !
      ! be used).                                                                          !
      !------------------------------------------------------------------------------------!
      do k = csite%nlev_sfcwater(ipa)+1,nzs
         csite%sfcwater_energy(k,ipa)  = 0.
         csite%sfcwater_mass(k,ipa)    = 0.
         csite%sfcwater_depth(k,ipa)   = 0.
         if (k == 1) then
            csite%sfcwater_fracliq(k,ipa) = csite%soil_fracliq(nzg,ipa)
            csite%sfcwater_tempk(k,ipa)   = csite%soil_tempk(nzg,ipa)
         else
            csite%sfcwater_fracliq(k,ipa) = csite%sfcwater_fracliq(k-1,ipa)
            csite%sfcwater_tempk(k,ipa)   = csite%sfcwater_tempk(k-1,ipa)
         end if
      end do
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !     Cohort variables.  Here we must check whether the cohort was really solved or  !
      ! it was skipped after being flagged as "unsafe".  Here the reason why it was flag-  !
      ! ged as such matters.                                                               !
      !------------------------------------------------------------------------------------!
      do ico = 1,cpatch%ncohorts

         if (initp%solvable(ico)) then
            !------------------------------------------------------------------------------!
            !    The cohort was solved, update internal energy and water, and re-calculate !
            ! temperature.  Note that energy may need to be scaled back.                   !
            !------------------------------------------------------------------------------!
            cpatch%veg_water(ico)  = sngloff(initp%veg_water(ico),tiny_offset)
            tmp_energy             = initp%veg_energy(ico)                                 &
                                   + (dble(cpatch%hcapveg(ico))-initp%hcapveg(ico))        &
                                   * initp%veg_temp(ico)

            cpatch%veg_energy(ico) = sngloff(tmp_energy,tiny_offset)
            call qwtk(cpatch%veg_energy(ico),cpatch%veg_water(ico),cpatch%hcapveg(ico)     &
                     ,cpatch%veg_temp(ico),cpatch%veg_fliq(ico))
         elseif (cpatch%hite(ico) <=  csite%total_snow_depth(ipa)) then
            !------------------------------------------------------------------------------!
            !    For plants buried in snow, fix the leaf temperature to the snow temper-   !
            ! ature of the layer that is the closest to the leaves.                        !
            !------------------------------------------------------------------------------!
            kclosest = 1
            do k = csite%nlev_sfcwater(ipa), 1, -1
               if (sum(csite%sfcwater_depth(1:k,ipa)) > cpatch%hite(ico)) kclosest = k
            end do
            cpatch%veg_temp(ico)   = csite%sfcwater_tempk(kclosest,ipa)
            cpatch%veg_fliq(ico)   = 0.
            cpatch%veg_water(ico)  = 0.
            cpatch%veg_energy(ico) = cpatch%hcapveg(ico) * cpatch%veg_temp(ico)
         else
            !------------------------------------------------------------------------------!
            !     For plants with minimal foliage or very sparse patches, fix the leaf     !
            ! temperature to the canopy air space and force veg_water to be zero.          !
            !------------------------------------------------------------------------------!
            cpatch%veg_temp(ico)   = csite%can_temp(ipa)
            cpatch%veg_fliq(ico)   = 0.
            cpatch%veg_water(ico)  = 0. 
            cpatch%veg_energy(ico) = cpatch%hcapveg(ico) * cpatch%veg_temp(ico)
         end if

         !---------------------------------------------------------------------------------!
         !     Final sanity check.  This should be removed soon, since it should never     !
         ! happen (well, if this still happens, then it's a bug, and we should remove the  !
         ! bug first...).                                                                  !
         !---------------------------------------------------------------------------------!
         if (cpatch%veg_temp(ico) < sngl(rk4min_veg_temp) .or.                             &
             cpatch%veg_temp(ico) > sngl(rk4max_veg_temp)   ) then
            write (unit=*,fmt='(80a)')         ('=',k=1,80)
            write (unit=*,fmt='(a)')           'FINAL VEG_TEMP IS WRONG IN INITP2MODELP'
            write (unit=*,fmt='(80a)')         ('-',k=1,80)
            write (unit=*,fmt='(a,1x,f9.4)')   ' + LONGITUDE:    ',rk4met%lon
            write (unit=*,fmt='(a,1x,f9.4)')   ' + LATITUDE:     ',rk4met%lat
            write (unit=*,fmt='(a,1x,i6)')     ' + POLYGON:      ',ipy
            write (unit=*,fmt='(a,1x,i6)')     ' + SITE:         ',isi
            write (unit=*,fmt='(a,1x,i6)')     ' + PATCH:        ',ipa
            write (unit=*,fmt='(a,1x,i6)')     ' + COHORT:       ',ico
            write (unit=*,fmt='(a)')           ' + PATCH AGE:    ',csite%age(ipa)
            write (unit=*,fmt='(a,1x,es12.4)') '   - AGE:        ',csite%age(ipa)
            write (unit=*,fmt='(a,1x,i6)')     '   - DIST_TYPE:  ',csite%dist_type(ipa)
            write (unit=*,fmt='(a)')           ' + BUFFER_COHORT (initp):'
            write (unit=*,fmt='(a,1x,es12.4)') '   - ENERGY:     ',initp%veg_energy(ico)
            write (unit=*,fmt='(a,1x,es12.4)') '   - WATER:      ',initp%veg_water(ico)
            write (unit=*,fmt='(a,1x,es12.4)') '   - TEMPERATURE:',initp%veg_temp(ico)
            write (unit=*,fmt='(a,1x,es12.4)') '   - FRACLIQ:    ',initp%veg_fliq(ico)
            write (unit=*,fmt='(a,1x,es12.4)') '   - HEAT_CAP:   ',initp%hcapveg(ico)
            write (unit=*,fmt='(a)')           ' + STATE_COHORT (cpatch):'
            write (unit=*,fmt='(a,1x,es12.4)') '   - ENERGY:     ',cpatch%veg_energy(ico)
            write (unit=*,fmt='(a,1x,es12.4)') '   - WATER:      ',cpatch%veg_water(ico)
            write (unit=*,fmt='(a,1x,es12.4)') '   - TEMPERATURE:',cpatch%veg_temp(ico)
            write (unit=*,fmt='(a,1x,es12.4)') '   - FRACLIQ:    ',cpatch%veg_fliq(ico)
            write (unit=*,fmt='(a,1x,es12.4)') '   - HEAT_CAP:   ',cpatch%hcapveg(ico)
            write (unit=*,fmt='(80a)') ('-',k=1,80)
            call print_rk4patch(initp, csite,ipa)
            call fatal_error('extreme vegetation temperature','initp2modelp'               &
                            &,'rk4_driver.f90')
         end if
         !---------------------------------------------------------------------------------!
      end do


      ksn   = csite%nlev_sfcwater(ipa)
      nsoil = csite%ntext_soil(nzg,ipa)
      nlsw1 = max(1, ksn)
      call ed_grndvap(ksn,nsoil,csite%soil_water(nzg,ipa),csite%soil_energy(nzg,ipa)       &
                     ,csite%sfcwater_energy(nlsw1,ipa),csite%can_rhos(ipa)                 &
                     ,csite%can_shv(ipa),csite%ground_shv(ipa),csite%surface_ssh(ipa)      &
                     ,surface_temp,surface_fliq)
      return
   end subroutine initp2modelp
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    Currently not in use.                                                              !
   !---------------------------------------------------------------------------------------!
   subroutine canopy_atm_fluxes(csite,cpoly,ipa,isi)
      use ed_state_vars , only : polygontype  & ! structure
                               , sitetype     & ! structure
                               , patchtype    ! ! structure
      use consts_coms   , only : cpi          ! ! intent(in)
      implicit none
    
      !----- Arguments --------------------------------------------------------------------!
      type(polygontype) , target     :: cpoly
      type(sitetype)    , target     :: csite
      integer           , intent(in) :: ipa
      integer           , intent(in) :: isi
      !------------------------------------------------------------------------------------!

      call fatal_error('Decide how to set vels in canopy_atm_fluxes.'                      &
                     &,'canopy_atm_fluxes','rk4_driver.f90')

      !----- Calculate turbulent fluxes between atmosphere and canopy. --------------------!
      !    pis = cpoly%pi0 * cpi
      !    thetacan = pss%can_temp / pis
      !    if(thetacan.lt.cpoly%theta)then
      !       cpoly%vels = cpoly%vels_stab
      !    else
      !       cpoly%vels = cpoly%vels_unstab
      !    endif

      return
   end subroutine canopy_atm_fluxes
   !=======================================================================================!
   !=======================================================================================! 
end module rk4_driver
!==========================================================================================!
!==========================================================================================!
