!==========================================================================================!
!==========================================================================================!
! grell_cupar_driver.f90                                                                   !
!                                                                                          !
!    This file contains the driver that will compute the Cumulus based on Grell-Dévényi    !
! parameterization, and send the information needed by the other modules.                  !
!------------------------------------------------------------------------------------------!
subroutine grell_cupar_driver(banneron,cldd,clds)

   use catt_start        , only: &
           catt                  ! ! intent(in) - flag for CATT. 

   use extras            , only: &
           extra3d               ! ! intent(inout) - Extra scratch for CATT

   use grell_coms        , only: & ! All variables here are intent(in)
           closure_type          & ! Flag for closure to be used on dyn. control
          ,comp_down             & ! I will compute downdrafts.              [T/F]
          ,comp_noforc_cldwork   & ! I will compute no forced cloud work     [T/F]
          ,comp_modif_thermo     & ! I will compute the x_* variables        [T/F]
          ,checkmass             & ! I will check mass balance               [T/F]
          ,maxens_cap            & ! Ensemble size on level of capping inversion
          ,maxens_dyn            & ! Ensemble size on dynamic control
          ,maxens_eff            & ! Ensemble size on precipitation efficiency
          ,maxens_lsf            & ! Ensemble size on large scale perturbations
          ,mgmzp                 & ! Vertical grid size
          ,cap_maxs              & ! Maximum depth of capping inversion     [ hPa]
          ,cap_max_increment     & ! Extra cap_maxs due to upstream conv.   [ hPa]
          ,depth_min             & ! Minimum cloud depth to qualify it      [ hPa]
          ,depth_max             & ! Maximum cloud depth to qualify it      [ hPa]
          ,edtmax                & ! Maximum Grell's epsilon (dnmf/upmf)
          ,edtmin                & ! Minimum Grell's epsilon (dnmf/upmf)
          ,inv_ensdim            & ! Inverse of ensemble dimension size
          ,iupmethod             & ! Method to define updraft originating level.
          ,masstol               & ! Maximum mass leak allowed to happen    [ ---]
          ,max_heat              & ! Maximum heating scale                  [K/dy]
          ,pmass_left            & ! Fraction of mass left at the ground    [ ---]
          ,radius                & ! Radius, for entrainment rate.          [   m]
          ,relheight_down        & ! Relative height for downdraft origin   [ ---]
          ,wnorm_max             & ! Normalised trigger vertical velocity
          ,wnorm_increment       & ! Increments on wnorm_max for ensemble
          ,zkbmax                & ! Top height for updrafts to originate   [   m]
          ,zcutdown              & ! Top height for downdrafts to originate [   m]
          ,z_detr                ! ! Top height for downdraft detrainment   [   m]

   use io_params         , only: &
           frqanl                ! ! intent(in) - Frequency of analysis.

   use mem_basic         , only: &
           co2_on                & ! intent(in) - Flag for CO2 presence.
          ,co2con                & ! intent(in) - CO2 mixing ratio if constant.
          ,basic_g               ! ! intent(in) - Basic variables structure

   use mem_cuparm        , only: &
           confrq                & ! intent(in)    - Convective frequency 
          ,cptime                & ! intent(in)    - Time to start computing the cloud.
          ,cuparm_g              & ! intent(inout) - Structure with convection
          ,nclouds               ! ! intent(in)    - # of clouds available.

   use mem_ensemble, only:       & !
           ensemble_e            & ! intent(inout) - Ensemble structure
          ,zero_ensemble         ! ! subroutine to flush scratch arrays to zero.

   use mem_grid          , only: & ! All variables are intent(in)
           deltax                & ! Current grid resolution in x (m)
          ,deltay                & ! Current grid resolution in y (m)
          ,dtlt                  & ! current grid time step
          ,dtlongn               & ! all grids time step
          ,grid_g                & ! grid structure
          ,initial               & ! flag for initial run
          ,jdim                  & ! j-dimension (usually 1 except for rare 2-D runs)
          ,naddsc                & ! number of additional scalars
          ,ngrid                 & ! current grid ID
          ,ngrids                & ! number of grids in this run. Used @ initialization.
          ,time                  & ! simulation elapsed time
          ,zt                    & ! Vertical thermodynamic levels for current grid
          ,zm                    ! ! Vertical momentum levels for current grid

   use mem_mass          , only: &
           imassflx              & ! intent(in) - flag for mass flux outout
          ,mass_g                ! ! intent(inout) - mass structure

   use mem_micro         , only: &
           micro_g               ! ! intent(in) - microphysics structure
           
   use mem_scalar        , only: &
           scalar_g      ! ! intent(inout) - This is used by CATT only.

   use mem_scratch       , only: &
           scratch       & ! intent(out) - Scratch array, to save old info.
          ,vctr6         & ! intent(out) - Scratch, contains the liquid water mix. ratio.
          ,vctr7         & ! intent(out) - Scratch, contains the ice mixing ratio.
          ,vctr8         & ! intent(out) - Scratch, contains the column CO2 mixing ratio.
          ,vctr9         & ! intent(out) - Scratch, contains the vertical velocity sigma.
          ,vctr28        ! ! intent(out) - Scratch, contains the convective CO2 forcing.
   
   use mem_scratch_grell , only: &
           zero_scratch_grell    ! ! Subroutine - Flushes scratch variables to zero.

   use mem_tend          , only: &
           tend                  ! ! intent(inout) - Tendency structure

   use mem_turb          , only: &
           turb_g                & ! intent(in) - turbulence structure
          ,idiffk                ! ! intent(in) - turbulence closure flag

   use micphys           , only: &
           availcat              ! ! intent(in) - Flag: hydrometeor is available [T|F] 

   use node_mod          , only: & ! All variables are intent(in)
           mynum                 & ! This node ID
          ,mxp                   & ! # of x-points for current grid in this node
          ,myp                   & ! # of y-points for current grid in this node
          ,mzp                   & ! # of z-points for current grid in this node
          ,ia                    & ! westernmost point for current grid in this node
          ,iz                    & ! easternmost point for current grid in this node
          ,ja                    & ! southernmost point for current grid in this node
          ,jz                    ! ! northernmost point for current grid in this node

   use therm_lib           , only: &
           level                 ! ! intent(in) - Phase complexity level

   implicit none


   !----- Arguments. ----------------------------------------------------------------------!
   logical, intent(in)               :: banneron         ! Flag for banner printing.
   integer, intent(in)               :: cldd             ! Deepest Grell cloud.
   integer, intent(in)               :: clds             ! Shallowest Grell cloud.
   !----- Local variables. ----------------------------------------------------------------!
   integer                           :: icld             ! Cloud type counter.
   integer                           :: i,j              ! Counters for x and y
   integer                           :: iscl             ! Scalar counter
   integer                           :: ccc              ! For initialization cloud loop
   real                              :: dti              ! confrq/frqanl
   integer                           :: iedt,imbp,idync  ! Counters - debugging
   !----- Local constants. ----------------------------------------------------------------!
   logical, parameter                :: printing=.false. ! Should I print debugging stuff?
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !                           ----- Main convection block -----                           !
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
  
   !---------------------------------------------------------------------------------------!
   ! 1.  Flushing the feedback to zero, so the default is no convection happened.          !
   !---------------------------------------------------------------------------------------!
   do icld=cldd,clds
      call azero(mzp*mxp*myp,cuparm_g(ngrid)%thsrc  (:,:,:,icld))
      call azero(mzp*mxp*myp,cuparm_g(ngrid)%rtsrc  (:,:,:,icld))
      if (co2_on) call azero(mzp*mxp*myp,cuparm_g(ngrid)%co2src (:,:,:,icld))
      call azero(mzp*mxp*myp,cuparm_g(ngrid)%cuprliq(:,:,:,icld))
      call azero(mzp*mxp*myp,cuparm_g(ngrid)%cuprice(:,:,:,icld))

      call azero(    mxp*myp,cuparm_g(ngrid)%conprr (  :,:,icld))
      call azero(    mxp*myp,cuparm_g(ngrid)%aadn   (  :,:,icld))
      call azero(    mxp*myp,cuparm_g(ngrid)%aaup   (  :,:,icld))
      call azero(    mxp*myp,cuparm_g(ngrid)%areadn (  :,:,icld))
      call azero(    mxp*myp,cuparm_g(ngrid)%areaup (  :,:,icld))
      call azero(    mxp*myp,cuparm_g(ngrid)%dnmf   (  :,:,icld))
      call azero(    mxp*myp,cuparm_g(ngrid)%edt    (  :,:,icld))
      call azero(    mxp*myp,cuparm_g(ngrid)%upmf   (  :,:,icld))
      call aone (    mxp*myp,cuparm_g(ngrid)%xierr  (  :,:,icld))
      call azero(    mxp*myp,cuparm_g(ngrid)%zjmin  (  :,:,icld))
      call azero(    mxp*myp,cuparm_g(ngrid)%zk22   (  :,:,icld))
      call azero(    mxp*myp,cuparm_g(ngrid)%zkdt   (  :,:,icld))
      call azero(    mxp*myp,cuparm_g(ngrid)%zkbcon (  :,:,icld))
      call azero(    mxp*myp,cuparm_g(ngrid)%zktop  (  :,:,icld))
   end do
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   ! 2. Flushing a scratch array to zero, and depending on whether CO2 is prognosed, this  !
   !    will store the large-scale tendency.                                               !
   !    scratch%vt3dj => large-scale CO2 mixing ratio forcing.                             !
   !---------------------------------------------------------------------------------------!
   if (co2_on) then
      call atob(mxp*myp*mzp,tend%co2t,scratch%vt3dj)
   else
      call azero(mzp*mxp*myp,scratch%vt3dj)
   end if
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   ! 3. Flush CATT-related variables to zero.                                              !
   !---------------------------------------------------------------------------------------!
   if (catt == 1) then
       dti = confrq(icld) / frqanl
       call azero(mzp*mxp*myp,scratch%scr1)
       call azero(mzp*mxp*myp,scratch%vt3dp)
       call azero(mzp*mxp*myp,extra3d(1,ngrid)%d3)
       call azero(mzp*mxp*myp,extra3d(2,ngrid)%d3)
       if (mod(time,dble(frqanl)) < dtlongn(1)) then
          call azero(mxp*myp*mzp,extra3d(4,ngrid)%d3)
       end if
   end if
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Big loop accross the horizontal grid points. Now we call the model column by      !
   ! column.                                                                               !
   !---------------------------------------------------------------------------------------!
   jloop: do j=ja,jz
      iloop: do i=ia,iz

         ![[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[!
         !---------------------------------------------------------------------------------!
         ! 4. We now initialise some variables that don't depend on the cloud spectral     !
         !    size because they must be done only once.                                    !
         !---------------------------------------------------------------------------------!

         !---------------------------------------------------------------------------------!
         ! 4a. Finding the ice and liquid mixing ratio.  Since ice and liquid may not be   !
         !     solved for some user's configuration (level = 1), we must check that first. !
         !---------------------------------------------------------------------------------!
         if (banneron) write (unit=60+mynum,fmt='(2(a,1x,i5,1x))')                         &
                             '       [~] Maybe calling integ_liq_ice... i=',i,'j=',j
         select case (level)
         !------ No condensed phase, simply set up both to zero. --------------------------!
         case (0,1)
            call azero(mzp,vctr6(1:mzp))
            call azero(mzp,vctr7(1:mzp))
         !----- Liquid condensation only, use saturation adjustment -----------------------!
         case (2)
             call atob(mzp,micro_g(ngrid)%rcp(:,i,j),vctr6(1:mzp))
             call azero(mzp,vctr7(1:mzp))
         case (3)
            call integ_liq_ice(mzp,availcat                                                &
              , micro_g(ngrid)%rcp           (:,i,j), micro_g(ngrid)%rrp           (:,i,j) &
              , micro_g(ngrid)%rpp           (:,i,j), micro_g(ngrid)%rsp           (:,i,j) &
              , micro_g(ngrid)%rap           (:,i,j), micro_g(ngrid)%rgp           (:,i,j) &
              , micro_g(ngrid)%rhp           (:,i,j), micro_g(ngrid)%q6            (:,i,j) &
              , micro_g(ngrid)%q7            (:,i,j), vctr6(1:mzp)                         &
              , vctr7(1:mzp)                        )
         end select
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         ! 4b. Copy CO2 array to scratch variable, if CO2 is actively prognosed.  Other-   !
         !     wise, put the constant CO2 mixing ratio, although it will not be really     !
         !     used.                                                                       !
         !---------------------------------------------------------------------------------!
         if (co2_on) then
            call atob(mzp,basic_g(ngrid)%co2p(:,i,j),vctr8(1:mzp))
            call azero(mzp,vctr28(1:mzp))
         else
            call ae0(mzp,vctr8(1:mzp),co2con(1))
            call azero(mzp,vctr28(1:mzp))
         end if
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         ! 4c. Copying sigma-w to a scratch array. This is because it is available only    !
         !     for idiffk = 1 or idiffk = 7. It's really unlikely that one would use the   !
         !     other TKE-related schemes with cumulus parameterization though, because     !
         !     they are LES schemes. If that happens, use special flag (sigmaw=0).         !
         !---------------------------------------------------------------------------------!
         select case (idiffk(ngrid))
         case (1,7)
            call atob(mzp,turb_g(ngrid)%sigw(:,i,j),vctr9(1:mzp))
         case default
            call azero(mzp,vctr9(1:mzp))
         end select
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         ! 4d. Flushes all scratch variables to zero.                                      !
         !---------------------------------------------------------------------------------!
         call zero_scratch_grell(3)


         !---------------------------------------------------------------------------------!
         ! 4e. Initialise grid-related variables (how many levels, offset, etc.)           !
         !---------------------------------------------------------------------------------!
         if (banneron) write (unit=60+mynum,fmt='(2(a,1x,i5,1x))')                         &
                             '       [~] Calling initial_grid_grell... i=',i,'j=',j
         call initial_grid_grell(mzp,deltax,deltay,zt(1:mzp),zm(1:mzp)                     &
              , grid_g(ngrid)%flpw             (i,j), grid_g(ngrid)%rtgt             (i,j) &
              , confrq                        (cldd), turb_g(ngrid)%kpbl             (i,j) )


         !---------------------------------------------------------------------------------!
         ! 4f. Copying the tendencies to the scratch array.                                !
         !---------------------------------------------------------------------------------!
         if (banneron) write (unit=60+mynum,fmt='(2(a,1x,i5,1x))')                         &
                             '       [~] Calling initial_tend_grell... i=',i,'j=',j
         call initial_tend_grell(mzp,mxp,myp,i,j,tend%tht,tend%tket,tend%rtt,scratch%vt3dj)


         !---------------------------------------------------------------------------------!
         ! 4g. Initialising pressure, temperature, and mixing ratio. This will include the !
         !     effect of previously called shallow convection in case it happened. If this !
         !     is the shallowest cloud, or if the shallower clouds didn't happen, this     !
         !     will simply copy the BRAMS fields.                                          !
         !---------------------------------------------------------------------------------!
         if (banneron) write (unit=60+mynum,fmt='(2(a,1x,i5,1x))')                         &
                             '       [~] Calling initial_thermo_grell... i=',i,'j=',j
         call initial_thermo_grell(mzp,dtlt         , basic_g(ngrid)%thp           (:,i,j) &
              , basic_g(ngrid)%theta         (:,i,j), basic_g(ngrid)%rtp           (:,i,j) &
              , vctr8                        (1:mzp), basic_g(ngrid)%pi0           (:,i,j) &
              , basic_g(ngrid)%pp            (:,i,j), basic_g(ngrid)%pc            (:,i,j) &
              , basic_g(ngrid)%wp            (:,i,j), basic_g(ngrid)%dn0           (:,i,j) &
              , turb_g(ngrid)%tkep           (:,i,j), vctr6                        (1:mzp) &
              , vctr7                        (1:mzp), vctr9                        (1:mzp))


         !---------------------------------------------------------------------------------!
         ! 5. We will now go through the cloud sizes for the first time, in order to solve !
         !    the static control.                                                          !
         !---------------------------------------------------------------------------------!
         staticloop: do icld = cldd,clds

            !------------------------------------------------------------------------------!
            ! 5a. Reset the entire ensemble structure for this cloud.                      !
            !------------------------------------------------------------------------------!
            call zero_scratch_grell(1)
            call zero_ensemble(ensemble_e(icld))


            !------------------------------------------------------------------------------!
            ! 5b. Initialise the remainder Grell's scratch variables. This is done only if !
            !     the user asked for this check, otherwise I will pretend that convection  !
            !     never happen at the surroundings. This requires full grid information,   !
            !     and needs some old information as well, thus the scratch structures.     !
            !------------------------------------------------------------------------------!
            if (banneron) write (unit=60+mynum,fmt='(3(a,1x,i5,1x))')                      &
                   '       [~] Calling initial_upstream_grell... i=',i,'j=',j,'icld=',icld
            call initial_winds_grell(comp_down(icld),mzp,mxp,myp,i,j,jdim                  &
                                    , cuparm_g(ngrid)%dnmf (:,:,icld)                      &
                                    , basic_g(ngrid)%up                                    &
                                    , basic_g(ngrid)%vp                                    &
                                    , ensemble_e(icld)%prev_dnmf                           )
         
            !------------------------------------------------------------------------------!
            ! 5c. Call the subroutine which will deal with the static control.             !
            !------------------------------------------------------------------------------!
            if (banneron) write (unit=60+mynum,fmt='(3(a,1x,i5,1x))')                      &
                       '       [~] Calling grell_cupar_static... i=',i,'j=',j,'icld=',icld
            call grell_cupar_static(comp_down(icld),comp_noforc_cldwork,checkmass          &
                                   ,iupmethod,maxens_cap,maxens_eff,mgmzp,cap_maxs(icld)   &
                                   ,cap_max_increment,wnorm_max(icld),wnorm_increment      &
                                   ,depth_min(icld),depth_max(icld),edtmax,edtmin,masstol  &
                                   ,pmass_left,radius(icld),relheight_down,zkbmax(icld)    &
                                   ,zcutdown(icld),z_detr(icld)                            &
                                   ,cuparm_g(ngrid)%aadn(i,j,icld)                         &
                                   ,cuparm_g(ngrid)%aaup(i,j,icld)                         &
                                   ,ensemble_e(icld)%dellatheiv_eff                        &
                                   ,ensemble_e(icld)%dellathil_eff                         &
                                   ,ensemble_e(icld)%dellaqtot_eff                         &
                                   ,ensemble_e(icld)%dellaco2_eff                          &
                                   ,ensemble_e(icld)%pw_eff                                &
                                   ,ensemble_e(icld)%edt_eff                               &
                                   ,ensemble_e(icld)%aatot0_eff                            &
                                   ,ensemble_e(icld)%aatot_eff                             &
                                   ,ensemble_e(icld)%ierr_cap                              &
                                   ,ensemble_e(icld)%jmin_cap                              &
                                   ,ensemble_e(icld)%k22_cap                               &
                                   ,ensemble_e(icld)%kbcon_cap                             &
                                   ,ensemble_e(icld)%kdet_cap                              &
                                   ,ensemble_e(icld)%kstabi_cap                            &
                                   ,ensemble_e(icld)%kstabm_cap                            &
                                   ,ensemble_e(icld)%ktop_cap                              &
                                   ,ensemble_e(icld)%pwav_cap                              &
                                   ,ensemble_e(icld)%pwev_cap                              &
                                   ,ensemble_e(icld)%wbuoymin_cap                          &
                                   ,ensemble_e(icld)%cdd_cap                               &
                                   ,ensemble_e(icld)%cdu_cap                               &
                                   ,ensemble_e(icld)%mentrd_rate_cap                       &
                                   ,ensemble_e(icld)%mentru_rate_cap                       &
                                   ,ensemble_e(icld)%dbyd_cap                              &
                                   ,ensemble_e(icld)%dbyu_cap                              &
                                   ,ensemble_e(icld)%etad_cld_cap                          &
                                   ,ensemble_e(icld)%etau_cld_cap                          &
                                   ,ensemble_e(icld)%rhod_cld_cap                          &
                                   ,ensemble_e(icld)%rhou_cld_cap                          &
                                   ,ensemble_e(icld)%qliqd_cld_cap                         &
                                   ,ensemble_e(icld)%qliqu_cld_cap                         &
                                   ,ensemble_e(icld)%qiced_cld_cap                         &
                                   ,ensemble_e(icld)%qiceu_cld_cap                         &
                                   ,i,j,icld,mynum                                         )
         end do staticloop

         !---------------------------------------------------------------------------------!
         ! 6. We now compute the dynamic control, which will determine the characteristic  !
         !    mass flux for each Grell cumulus cloud.                                      !
         !---------------------------------------------------------------------------------!
         if (banneron) write (unit=60+mynum,fmt='(2(a,1x,i5,1x))')                         &
                             '       [~] Calling grell_cupar_dynamic... i=',i,'j=',j
         call grell_cupar_dynamic(cldd,clds,nclouds,dtlt,maxens_cap,maxens_eff,maxens_lsf  &
                                 ,maxens_dyn,mgmzp,closure_type,comp_modif_thermo          &
                                 ,comp_down,mynum,i,j,printing)

         !---------------------------------------------------------------------------------!
         ! 7. We now go through the cloud sizes again, to compute the feedback to the      !
         !    large scale and any other cloud-dependent variable that needs to be          !
         !    computed.
         !---------------------------------------------------------------------------------!
         fdbkloop: do icld=cldd,clds
            !------------------------------------------------------------------------------!
            ! 7a. Computing the feedback to the large-scale per se.  This will determine   !
            !     the forcing of large-scale variables, namely ice-liquid potential        !
            !     temperature, total water vapour mixing ratio, and carbon dioxide mixing  !
            !     ratio.                                                                   !
            !------------------------------------------------------------------------------!
            if (banneron) write (unit=60+mynum,fmt='(3(a,1x,i5,1x))')                      &
                   '       [~] Calling grell_cupar_feedback... i=',i,'j=',j,'icld=',icld
            call grell_cupar_feedback( comp_down(icld),mgmzp,maxens_cap,maxens_eff         &
                                     , maxens_lsf,maxens_dyn,inv_ensdim,max_heat(icld)     &
                                     , ensemble_e(icld)%dnmf_ens                           &
                                     , ensemble_e(icld)%upmf_ens                           &
                                     , ensemble_e(icld)%dellathil_eff                      &
                                     , ensemble_e(icld)%dellaqtot_eff                      &
                                     , ensemble_e(icld)%dellaco2_eff                       &
                                     , ensemble_e(icld)%pw_eff                             &
                                     , ensemble_e(icld)%ierr_cap                           &
                                     , cuparm_g(ngrid)%upmf(i,j,icld)                      &
                                     , cuparm_g(ngrid)%upmf(i,j,icld)                      &
                                     , cuparm_g(ngrid)%edt (i,j,icld)                      &
                                     , i,j,icld,mynum)

            !------------------------------------------------------------------------------!
            ! 7b. We now compute the other output variables such as mass fluxes and inter- !
            !     esting levels.                                                           !
            !------------------------------------------------------------------------------!
            if (banneron) write (unit=60+mynum,fmt='(3(a,1x,i5,1x))')                      &
                     '       [~] Calling grell_cupar_output... i=',i,'j=',j,'icld=',icld
            call grell_cupar_output(comp_down(icld),mzp,mgmzp,maxens_cap                   &
              , grid_g(ngrid)%rtgt             (i,j), zt(1:mzp)                            &
              , zm(1:mzp)                           , cuparm_g(ngrid)%dnmf      (i,j,icld) &
              , cuparm_g(ngrid)%upmf      (i,j,icld), ensemble_e(icld)%ierr_cap            &
              , ensemble_e(icld)%kdet_cap           , ensemble_e(icld)%k22_cap             &
              , ensemble_e(icld)%kbcon_cap          , ensemble_e(icld)%jmin_cap            &
              , ensemble_e(icld)%ktop_cap           , ensemble_e(icld)%wbuoymin_cap        &
              , ensemble_e(icld)%etad_cld_cap       , ensemble_e(icld)%mentrd_rate_cap     &
              , ensemble_e(icld)%cdd_cap            , ensemble_e(icld)%dbyd_cap            &
              , ensemble_e(icld)%rhod_cld_cap       , ensemble_e(icld)%etau_cld_cap        &
              , ensemble_e(icld)%mentru_rate_cap    , ensemble_e(icld)%cdu_cap             &
              , ensemble_e(icld)%dbyu_cap           , ensemble_e(icld)%rhou_cld_cap        &
              , ensemble_e(icld)%qliqd_cld_cap      , ensemble_e(icld)%qliqu_cld_cap       &
              , ensemble_e(icld)%qiced_cld_cap      , ensemble_e(icld)%qiceu_cld_cap       &
              , cuparm_g(ngrid)%xierr     (i,j,icld), cuparm_g(ngrid)%zjmin     (i,j,icld) &
              , cuparm_g(ngrid)%zk22      (i,j,icld), cuparm_g(ngrid)%zkbcon    (i,j,icld) &
              , cuparm_g(ngrid)%zkdt      (i,j,icld), cuparm_g(ngrid)%zktop     (i,j,icld) &
              , cuparm_g(ngrid)%conprr    (i,j,icld), cuparm_g(ngrid)%thsrc   (:,i,j,icld) &
              , cuparm_g(ngrid)%rtsrc   (:,i,j,icld), vctr28                       (1:mzp) &
              , cuparm_g(ngrid)%areadn    (i,j,icld), cuparm_g(ngrid)%areaup    (i,j,icld) &
              , cuparm_g(ngrid)%cuprliq (:,i,j,icld), cuparm_g(ngrid)%cuprice (:,i,j,icld))

            !------------------------------------------------------------------------------!
            ! 7c. Copy the CO2 forcing if CO2 is available.                                !
            !------------------------------------------------------------------------------!
            if (co2_on) call atob(mzp,vctr28(1:mzp),cuparm_g(ngrid)%co2src(:,i,j,icld))

            !------------------------------------------------------------------------------!
            ! 7d. If the user needs the full mass flux information to drive Lagrangian     !
            !     models do it now.                                                        !
            !------------------------------------------------------------------------------!
            if (banneron) write (unit=60+mynum,fmt='(3(a,1x,i5,1x))')                      &
               '       [~] Maybe calling prep_convflx_to_mass... i=',i,'j=',j,'icld=',icld
            if (imassflx == 1) then
               call prep_convflx_to_mass(mzp,mgmzp,maxens_cap                              &
                  , cuparm_g(ngrid)%dnmf    (i,j,icld), cuparm_g(ngrid)%upmf    (i,j,icld) &
                  , ensemble_e(icld)%ierr_cap         , ensemble_e(icld)%jmin_cap          &
                  , ensemble_e(icld)%k22_cap          , ensemble_e(icld)%kbcon_cap         &
                  , ensemble_e(icld)%kdet_cap         , ensemble_e(icld)%ktop_cap          &
                  , ensemble_e(icld)%cdd_cap          , ensemble_e(icld)%cdu_cap           &
                  , ensemble_e(icld)%mentrd_rate_cap  , ensemble_e(icld)%mentru_rate_cap   &
                  , ensemble_e(icld)%etad_cld_cap     , ensemble_e(icld)%etau_cld_cap      &           
                  , mass_g(ngrid)%cfxdn   (:,i,j,icld), mass_g(ngrid)%cfxup   (:,i,j,icld) &
                  , mass_g(ngrid)%dfxdn   (:,i,j,icld), mass_g(ngrid)%dfxup   (:,i,j,icld) &
                  , mass_g(ngrid)%efxdn   (:,i,j,icld), mass_g(ngrid)%efxup   (:,i,j,icld) )
            end if

            !------------------------------------------------------------------------------!
            ! 7e. Here are CATT-related procedures. It does pretty much the same as it     !
            !     used to, except that I took the statistics out from inside               !
            !     grell_cupar_main and brought trans_conv_mflx inside the i/j loop.        !
            !------------------------------------------------------------------------------!
            if (banneron) write (unit=60+mynum,fmt='(2(a,1x,i5,1x))')                      &
                '       [~] Maybe calling grell_massflx_stats... i=',i,'j=',j,'icld=',icld
            if (catt == 1) then
               if (icld == 1) then
                  call grell_massflx_stats(mzp,icld,.true.,dti,maxens_dyn,maxens_lsf       &
                                          , maxens_eff,maxens_cap,inv_ensdim,closure_type  &
                                          , ensemble_e(icld)%ierr_cap                      &
                                          , ensemble_e(icld)%upmf_ens                      &
                                          , extra3d(1,ngrid)%d3        (:,i,j)             &
                                          , extra3d(4,ngrid)%d3        (:,i,j)             )
                  call grell_massflx_stats(mzp,icld,.false.,dti,maxens_dyn,maxens_lsf      &
                                          , maxens_eff,maxens_cap,inv_ensdim,closure_type  &
                                          , ensemble_e(icld)%ierr_cap                      &
                                          , ensemble_e(icld)%upmf_ens                      &
                                          , extra3d(1,ngrid)%d3        (:,i,j)             &
                                          , extra3d(4,ngrid)%d3        (:,i,j)             )
               elseif (icld == 2) then
                  call grell_massflx_stats(mzp,icld,.true.,dti,maxens_dyn,maxens_lsf       &
                                          , maxens_eff,maxens_cap,inv_ensdim,closure_type  &
                                          , ensemble_e(icld)%ierr_cap                      &
                                          , ensemble_e(icld)%upmf_ens                      &
                                          , extra3d(2,ngrid)%d3        (:,i,j)             &
                                          , scratch%vt3dp                                  )
               end if
               do iscl=1,naddsc
                  call trans_conv_mflx(icld,iscl,mzp,mxp,myp,maxens_cap,i,j                &
                                      ,ensemble_e(icld)%ierr_cap                           &
                                      ,ensemble_e(icld)%jmin_cap                           &
                                      ,ensemble_e(icld)%k22_cap                            &
                                      ,ensemble_e(icld)%kbcon_cap                          &
                                      ,ensemble_e(icld)%kdet_cap                           &
                                      ,ensemble_e(icld)%kstabi_cap                         &
                                      ,ensemble_e(icld)%kstabm_cap                         &
                                      ,ensemble_e(icld)%ktop_cap                           &
                                      ,ensemble_e(icld)%cdd_cap                            &
                                      ,ensemble_e(icld)%cdu_cap                            &
                                      ,ensemble_e(icld)%mentrd_rate_cap                    &
                                      ,ensemble_e(icld)%mentru_rate_cap                    &
                                      ,ensemble_e(icld)%etad_cld_cap                       &
                                      ,ensemble_e(icld)%etau_cld_cap                       &
                                      ,cuparm_g(ngrid)%edt     (i,j,icld)                  &
                                      ,cuparm_g(ngrid)%upmf    (i,j,icld)                  &
                                      ,scalar_g(iscl,ngrid)%sclt                           )
               end do
            end if
         end do fdbkloop


      end do iloop
   end do jloop
   !---------------------------------------------------------------------------------------!

   return 
end subroutine grell_cupar_driver
!==========================================================================================!
!==========================================================================================!
