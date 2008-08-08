!==========================================================================================!
! grell_cupar_driver.f90                                                                   !
!                                                                                          !
!    This file contains the driver that will compute the Cumulus based on Grell-Dévényi    !
! parameterization, and send the information needed by the other modules.                  !
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine grell_cupar_driver(icld)

   use catt_start        , only: &
           catt                  ! ! intent(in) - flag for CATT. 

   use extras            , only: &
           extra3d               ! ! intent(inout) - Extra scratch for CATT

   use grell_coms        , only: & ! All variables here are intent(in)
           closure_type          & ! Flag for closure to be used on dyn. control
          ,comp_down             & ! I will compute downdrafts.           [T/F]
          ,comp_noforc_cldwork   & ! I will compute no forced cloud work  [T/F]
          ,comp_modif_thermo     & ! I will compute the x_* variables     [T/F]
          ,checkmass             & ! I will check mass balance            [T/F]
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
          ,iupstrm               & ! Check for upstream convection?
          ,masstol               & ! Maximum mass leak allowed to happen    [ ---]
          ,max_heat              & ! Maximum heating scale                  [K/dy]
          ,pmass_left            & ! Fraction of mass left at the ground    [ ---]
          ,radius                & ! Radius, for entrainment rate.          [   m]
          ,relheight_down        & ! Relative height for downdraft origin   [ ---]
          ,zkbmax                & ! Top height for updrafts to originate   [   m]
          ,zcutdown              & ! Top height for downdrafts to originate [   m]
          ,z_detr                ! ! Top height for downdraft detrainment   [   m]

   use io_params         , only: &
           frqanl                ! ! intent(in) - Frequency of analysis.

   use mem_basic         , only: &
           basic_g               ! ! intent(in) - Basic variables structure

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
           scratch       ! ! intent(inout) - Scratch array, to save old info.
   
   use mem_scratch_grell , only: &
           zero_scratch_grell    ! ! Subroutine - Flushes scratch variables to zero.

   use mem_tend          , only: &
           tend                  ! ! intent(inout) - Tendency structure

   use mem_turb          , only: &
           turb_g                ! ! intent(in) - turbulence structure

   use node_mod          , only: & ! All variables are intent(in)
           mynum                 & ! This node ID
          ,mxp                   & ! # of x-points for current grid in this node
          ,myp                   & ! # of y-points for current grid in this node
          ,mzp                   & ! # of z-points for current grid in this node
          ,ia                    & ! westernmost point for current grid in this node
          ,iz                    & ! easternmost point for current grid in this node
          ,ja                    & ! southernmost point for current grid in this node
          ,jz                    ! ! northernmost point for current grid in this node

   implicit none


   !------ Variable declaration -----------------------------------------------------------!
   integer, intent(in)               :: icld                 ! Current cloud type:
                                                             ! 1. shallow; 2. deep.
   integer                           :: i,j                  ! Counters for x and y
   integer                           :: iscl                 ! Scalar counter
   integer                           :: ccc                  ! For initialization cloud loop
   real                              :: dti                  ! confrq/frqanl
   integer                           :: iedt,imbp,idync      ! Counters - debugging
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   ! Initial settings - moved to mem_cuparm. Now the values are initialized there.         !
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !    If you reached this point, you are ready to run the parameterization.              !
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
     
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !                           ----- Main convection block -----                           !
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
  
   !---------------------------------------------------------------------------------------!
   ! 1.  Flushing the feedback to zero, so the default is no convection happened.          !
   !---------------------------------------------------------------------------------------!
   cuparm_g(ngrid)%thsrc(1:mzp,ia:iz,ja:jz,icld)  = 0.
   cuparm_g(ngrid)%rtsrc(1:mzp,ia:iz,ja:jz,icld)  = 0.
   cuparm_g(ngrid)%conprr(ia:iz,ja:jz,icld)       = 0.
   call azero (mxp*myp,scratch%vt2de(1))
   call azero (mxp*myp,scratch%vt2df(1))

   !---------------------------------------------------------------------------------------!
   ! 2. Accumulate all convection sources into scratch arrays. This is to remove           !
   !    instability as we move to deeper convection.                                       !
   !    scratch%vt3di => potential temperature forcing.                                    !
   !    scratch%vt3dj => total mixing ratio forcing.                                       !
   !---------------------------------------------------------------------------------------!
   call azero(mzp*mxp*myp,scratch%vt3di(1))
   call azero(mzp*mxp*myp,scratch%vt3dj(1))
   do ccc=icld+1,nclouds
      call accum(mxp*myp*mzp, scratch%vt3di(1), cuparm_g(ngrid)%thsrc(1,1,1,ccc))
      call accum(mxp*myp*mzp, scratch%vt3dj(1), cuparm_g(ngrid)%rtsrc(1,1,1,ccc))
   end do

   !---------------------------------------------------------------------------------------!
   ! 3. Initialize variables for CATT                                                      !
   !---------------------------------------------------------------------------------------!
   if (catt == 1) then
       dti = confrq(icld) / frqanl
       call azero(mzp*mxp*myp,scratch%scr1(1))
       select case (icld)
       case (1)
          call azero(mzp*mxp*myp,extra3d(1,ngrid)%d3(1,1,1))
          if(mod(time,dble(frqanl)) < dtlongn(1) ) then
             call azero(mxp*myp*mzp,extra3d(4,ngrid)%d3(1,1,1)  )  
          end if
       case (2)
          call azero(mzp*mxp*myp,extra3d(2,ngrid)%d3(1,1,1))
          call azero(mzp*mxp*myp,scratch%vt3dp(1))
       end select
   end if
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   ! 4. If upstream is on, save the error flag and previous downdraft in some scratch      !
   !    array before the grid loop, so they are safely stored.                             !
   !---------------------------------------------------------------------------------------!
   if (iupstrm == 1) then
      call atob(mxp*myp,cuparm_g(ngrid)%dnmf(1,1,icld),scratch%vt2de(1))
      call atob(mxp*myp,cuparm_g(ngrid)%xierr(1,1,icld),scratch%vt2df(1))
   end if
  
   !---------------------------------------------------------------------------------------!
   ! 5. Big loop accross the horizontal grid points. Now I will call the model column by   !
   !    column.                                                                            !
   !---------------------------------------------------------------------------------------!
   jloop: do j=ja,jz
      iloop: do i=ia,iz

         !---------------------------------------------------------------------------------!
         ! 5a. Flushes scratch variables to zero.                                          !
         !---------------------------------------------------------------------------------!
         call zero_scratch_grell()
         call zero_ensemble(ensemble_e(icld))
         
         !---------------------------------------------------------------------------------!
         ! 5b. Initialize grid-related variables (how many levels, offset, etc.)           !
         !---------------------------------------------------------------------------------!
         call initialize_grid_grell(mzp,deltax,deltay,zt(1:mzp),zm(1:mzp)                  &
              , grid_g(ngrid)%flpw             (i,j), grid_g(ngrid)%rtgt             (i,j) &
              , turb_g(ngrid)%kpbl             (i,j))

         !---------------------------------------------------------------------------------!
         ! 5c. Copying the tendencies to the scratch array. This happens apart from the    !
         !     other variables because tend is a scratch structure, so I need to provide   !
         !     i and j.                                                                    !
         !---------------------------------------------------------------------------------!
         call initialize_tend_grell(mzp,mxp,myp,i,j,tend%tht(1),tend%rtt(1),tend%pt(1)     &
                                   ,scratch%vt3di(1),scratch%vt3dj(1))


         !---------------------------------------------------------------------------------!
         ! 5d. Initializing pressure, temperature, and mixing ratio. This will include the !
         !     effect of previously called shallow convection in case it happened. If this !
         !     is a shallow cumulus call, I just reset thsrc and rtsrc, so it will add     !
         !     zero, silly, but harmless.                                                  !
         !---------------------------------------------------------------------------------!
         call initialize_thermo_grell(mzp,dtlt                                             &
              , basic_g(ngrid)%theta         (1,i,j), basic_g(ngrid)%rv            (1,i,j) &
              , basic_g(ngrid)%pi0           (1,i,j), basic_g(ngrid)%pp            (1,i,j) &
              , basic_g(ngrid)%wp            (1,i,j), basic_g(ngrid)%dn0           (1,i,j) &
              , turb_g(ngrid)%tkep           (1,i,j), micro_g(ngrid)%rcp           (1,i,j) )


         !---------------------------------------------------------------------------------!
         ! 5e. Initialize the remainder Grell's scratch variables. This is done only if    !
         !     the user asked for this check, otherwise I will pretend that convection     !
         !     never happen at the surroundings. This requires full grid information, and  !
         !     needs some old information as well, thus the scratch structures.            !
         !---------------------------------------------------------------------------------!
         call initialize_upstream_grell(comp_down(icld),iupstrm,mzp,mxp,myp,i,j,jdim       &
              , scratch%vt2de                    (1), scratch%vt2df                    (1) &
              , basic_g(ngrid)%up            (1,1,1), basic_g(ngrid)%vp            (1,1,1) )
         
         !---------------------------------------------------------------------------------!
         ! 5f. Call the main cumulus parameterization subroutine, which will deal with the !
         !     static control, dynamic control and feedback. This sends back the reference !
         !     upward and downward mass flux in kg/m²/s.                                   !
         !---------------------------------------------------------------------------------!
         call grell_cupar_main(closure_type(icld),comp_down(icld)                          &
              , comp_noforc_cldwork(icld),comp_modif_thermo(icld),checkmass,iupmethod      &
              , maxens_cap(icld),maxens_dyn(icld),maxens_eff(icld),maxens_lsf(icld)        &
              , mgmzp,cap_maxs(icld),cap_max_increment,dtlt,depth_min(icld)                &
              , depth_max(icld),edtmax,edtmin,inv_ensdim(icld),masstol,max_heat(icld)      &
              , pmass_left,radius(icld),relheight_down,zkbmax(icld),zcutdown(icld)         &
              , z_detr(icld)                                                               &
              , ensemble_e(icld)%edt_eff       (1,1), ensemble_e(icld)%dellahe_eff (1,1,1) &
              , ensemble_e(icld)%dellaq_eff  (1,1,1), ensemble_e(icld)%dellaqc_eff (1,1,1) &
              , ensemble_e(icld)%dellat_eff  (1,1,1), ensemble_e(icld)%pw_eff      (1,1,1) &
              , ensemble_e(icld)%dnmf_ens  (1,1,1,1), ensemble_e(icld)%upmf_ens  (1,1,1,1) &
              , cuparm_g(ngrid)%aadn      (i,j,icld), cuparm_g(ngrid)%aaup      (i,j,icld) &
              , cuparm_g(ngrid)%edt       (i,j,icld), cuparm_g(ngrid)%dnmf      (i,j,icld) &
              , cuparm_g(ngrid)%upmf      (i,j,icld), mynum                                )
         !---------------------------------------------------------------------------------!
         ! 5g. Compute the other output variables                                          !
         !---------------------------------------------------------------------------------!
         call grell_output(comp_down(icld),mzp,mgmzp, grid_g(ngrid)%rtgt             (i,j) &
              , zt(1:mzp)                           , zm(1:mzp)                            &
              , basic_g(ngrid)%pi0           (1,i,j), basic_g(ngrid)%pp            (1,i,j) &
              , cuparm_g(ngrid)%dnmf      (i,j,icld), cuparm_g(ngrid)%upmf      (i,j,icld) &
              , cuparm_g(ngrid)%xierr     (i,j,icld), cuparm_g(ngrid)%zjmin     (i,j,icld) &
              , cuparm_g(ngrid)%zk22      (i,j,icld), cuparm_g(ngrid)%zkbcon    (i,j,icld) &
              , cuparm_g(ngrid)%zkdt      (i,j,icld), cuparm_g(ngrid)%zktop     (i,j,icld) &
              , cuparm_g(ngrid)%conprr    (i,j,icld), cuparm_g(ngrid)%thsrc   (1,i,j,icld) &
              , cuparm_g(ngrid)%rtsrc   (1,i,j,icld), cuparm_g(ngrid)%areadn  (1,i,j,icld) &
              , cuparm_g(ngrid)%areaup  (1,i,j,icld), cuparm_g(ngrid)%cupcond (1,i,j,icld) )

         !---------------------------------------------------------------------------------!
         ! 5h. If the user needs the full mass flux information to drive Lagrangian models !
         !     do it now.                                                                  !
         !---------------------------------------------------------------------------------!
         if (imassflx == 1) then
            call prep_convflx_to_mass(mzp                                                  &
                  , cuparm_g(ngrid)%dnmf    (i,j,icld), cuparm_g(ngrid)%upmf    (i,j,icld) &
                  , mass_g(ngrid)%cfxdn   (1,i,j,icld), mass_g(ngrid)%cfxup   (1,i,j,icld) &
                  , mass_g(ngrid)%dfxdn   (1,i,j,icld), mass_g(ngrid)%dfxup   (1,i,j,icld) &
                  , mass_g(ngrid)%efxdn   (1,i,j,icld), mass_g(ngrid)%efxup   (1,i,j,icld) )
         end if

         !---------------------------------------------------------------------------------!
         ! 5i. Here are CATT-related procedures. It does pretty much the same as it used   !
         !     to, except that I took the statistics out from inside grell_cupar_main and  !
         !     brought trans_conv_mflx inside the i/j loop.                                !
         !---------------------------------------------------------------------------------!
         if (catt == 1) then
            if (icld == 1) then
               call grell_massflx_stats(mzp,icld,.true.,dti,maxens_dyn(icld)               &
                  , maxens_lsf(icld),maxens_eff(icld),maxens_cap(icld),inv_ensdim(icld)    &
                  , closure_type(icld)                , ensemble_e(icld)%upmf_ens(1,1,1,1) &
                  , extra3d(1,ngrid)%d3        (1,i,j), extra3d(4,ngrid)%d3        (1,i,j) )
               call grell_massflx_stats(mzp,icld,.false.,dti,maxens_dyn(icld)              &
                  , maxens_lsf(icld),maxens_eff(icld),maxens_cap(icld),inv_ensdim(icld)    &
                  , closure_type(icld)                , ensemble_e(icld)%upmf_ens(1,1,1,1) &
                  , extra3d(1,ngrid)%d3        (1,i,j), extra3d(4,ngrid)%d3        (1,i,j) )
            elseif (icld == 2) then
               call grell_massflx_stats(mzp,icld,.true.,dti,maxens_dyn(icld)               &
                  , maxens_lsf(icld), maxens_eff(icld),maxens_cap(icld),inv_ensdim(icld)   &
                  , closure_type(icld)                , ensemble_e(icld)%upmf_ens(1,1,1,1) &
                  , extra3d(2,ngrid)%d3        (1,i,j), scratch%vt3dp                  (1) )
            end if
            do iscl=1,naddsc
               call trans_conv_mflx(icld,iscl,mzp,mxp,myp,i,j                              &
                  , cuparm_g(ngrid)%edt     (i,j,icld), cuparm_g(ngrid)%upmf    (i,j,icld) &
                  , scalar_g(iscl,ngrid)%sclt(1)                                           )
            end do
         end if
      end do iloop
   end do jloop

   !---------------------------------------------------------------------------------------!
   ! 6. Now I will update the tendencies due to this type of cloud, as well as the         !
   !    precipitation                                                                      !
   !---------------------------------------------------------------------------------------!
   call accum(mxp*myp*mzp, tend%tht(1), cuparm_g(ngrid)%thsrc(1,1,1,icld))
   call accum(mxp*myp*mzp, tend%rtt(1), cuparm_g(ngrid)%rtsrc(1,1,1,icld))
   call update(mxp*myp,cuparm_g(ngrid)%aconpr(1,1),cuparm_g(ngrid)%conprr(1,1,icld)   &
              ,dtlt)

   return 
end subroutine grell_cupar_driver
!==========================================================================================!
!==========================================================================================!
