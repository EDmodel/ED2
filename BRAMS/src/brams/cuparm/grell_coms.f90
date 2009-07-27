!==========================================================================================!
! grell_coms.f90                                                                           !
!                                                                                          !
!    This is the module containing all parameters needed to define scratch arrays in       !
! Grell's parameterization and all parameters that define each cloud type. Variables that  !
! are read from RAMSIN are also stored here.                                               !
!==========================================================================================!
!==========================================================================================!

module grell_coms
  use grid_dims, only: maxclouds
  !----------------------------------------------------------------------------------------!
  !   Dimension related variables                                                          !
  !----------------------------------------------------------------------------------------!
  !----- Maximum number of layers among the grids calling Grell ---------------------------!
  integer                                       :: mgmzp
  
  !----- Number of large-scale forcing members --------------------------------------------!
  integer     , dimension(maxclouds)            :: maxens_lsf
  
  !----- Number of precipitation efficiency members ---------------------------------------!
  integer     , dimension(maxclouds)            :: maxens_eff
  
  !----- Number of dynamic control members, it depends on the closure type ----------------!
  integer     , dimension(maxclouds)            :: maxens_dyn
  
  !----- Number of control members on cap_maxs --------------------------------------------!
  integer     , dimension(maxclouds)            :: maxens_cap
  
  !----- 1./ (maxens_lsf * maxens_eff * maxens_dyn) ---------------------------------------!
  real        , dimension(maxclouds)            :: inv_ensdim
   
  !----------------------------------------------------------------------------------------!
  !   Flags to bypass part of the cumulus parameterization. Depending on the type of cloud !
  ! and closure_type, we can skip parts of the code and speed up the run.                  !
  !----------------------------------------------------------------------------------------!
  !----- I will compute cloud work function with no forcing -------------------------------! 
  logical     , dimension(maxclouds)            :: comp_noforc_cldwork
  !----- I will compute the cloud work function with the arbitrary mass flux --------------!
  logical     , dimension(maxclouds)            :: comp_modif_thermo
  !----- I will compute downdrafts --------------------------------------------------------!
  logical     , dimension(maxclouds)            :: comp_down

  !----------------------------------------------------------------------------------------!
  !    These variables are obtained from RAMSIN.                                           !
  !----------------------------------------------------------------------------------------!
  !----- Method to choose updraft originating level (1. Max. Moist static energy; 2. PBL) -!
  integer                 :: iupmethod
  
  !------ Check for upstream convection (0 - no; 1- check only; 2- check and use it) ------!
  integer                 :: iupstrm
  
  !----- Minimum depth for the cloud to be considered [m] ---------------------------------!
  real              , dimension(maxclouds) :: depth_min
  
  !----------------------------------------------------------------------------------------!
  !   Maximum depth of inversion capping [hPa]. If negative, then the absolute value is    !
  ! the percentage of area that needs to be have enough buoyancy to reach the LFC.         !
  !----------------------------------------------------------------------------------------!
  real              , dimension(maxclouds) :: cap_maxs
  
  !----- Closure type to use for dynamic control ------------------------------------------!
  character (len=2) , dimension(maxclouds) :: closure_type
  
  !------ Cloud radius, which will affect entrainment -------------------------------------!
  real              , dimension(maxclouds) :: radius
  
  !------ Maximum height in which updrafts can originate [m] ------------------------------!
  real              , dimension(maxclouds) :: zkbmax
  
  !------ Maximum heat rate allowed for feedback [K/day] ----------------------------------!
  real              , dimension(maxclouds)  :: max_heat
  
  !------ Height above which downdrafts cannot originate ----------------------------------!
  real              , dimension(maxclouds) :: zcutdown
  
  !------ Top of downdraft detrainment layer ----------------------------------------------!
  real              , dimension(maxclouds) :: z_detr

  !----------------------------------------------------------------------------------------!

  !----------------------------------------------------------------------------------------!
  !    Variables that should stay outside reach in RAMSIN. These encompasses variables     !
  ! that are too much detail for the user to handle, a handful of parameters and also some !
  ! debugging options.                                                                     !
  !----------------------------------------------------------------------------------------!
  !------ This is the maximum normalized vertical velocity in which convection is called --!
  real              , dimension(maxclouds) :: wnorm_max

  !------ Maximum depth that the cloud can have to still belong to this class -------------!
  real              , dimension(maxclouds) :: depth_max
  !------ Perform rigorous mass balance check every call ----------------------------------!
  logical , parameter  :: checkmass=.true.
  
  !----- Increment on maximum depth of inversion capping [Pa] -----------------------------!
  real    , parameter  :: cap_maxs_increment  = 1500.
  !----- Increment on maximum normalized vertical velocity in which convection is called --!
  real    , parameter  :: wnorm_increment     = 0.25
  !----------------------------------------------------------------------------------------!


  !----------------------------------------------------------------------------------------!
  !  These variables are parameters for various Grell's computation                        !
  !----------------------------------------------------------------------------------------!
  !------ Minimum diameter for clouds to develop downdrafts and rain ----------------------!
  real                              , parameter  :: min_down_radius = 900.
  
  !------ Bottom "height" for wind mean [Pa] ----------------------------------------------!
  real                              , parameter  :: pbotmean = 15000.
  
  !------ Pressure at the top for wind mean [Pa] ------------------------------------------!
  real                              , parameter  :: ptopmean = 30000.
  
  !------ Minimum wind speed to bother to compute direction. [m/s] ------------------------!
  real                              , parameter  :: vspeedmin = 5.
  
  !------ Epsilon is the ratio between reference downdraft and updraft mass fluxes --------!
  real                              , parameter  :: edtmax = .95  ! Upper bound
  real                              , parameter  :: edtmin = .20  ! Lower bound
  
  !------ Maximum acceptable PBL height ---------------------------------------------------!
  real                              , parameter ::  pblhmax = 3000.
  
  !------ Minimum cloud mixing ratio to consider the layer wet ----------------------------!
  real                              , parameter ::  rcpmin  = 1.e-5
  
  !------ Height relative to the top above which no downdrafts can occur ------------------!
  real                              , parameter ::  relheight_down = 0.6
  
  !------ Percentage of mass left when hitting the ground ---------------------------------!
  real                              , parameter ::  pmass_left     = 0.03
  
  !------ Increment cap_max by this amount for different elements -------------------------!
  real                              , parameter  :: cap_max_increment=20.
  
  !------ Maximum "leakage" of mass allowed (normalized) ----------------------------------!
  real, parameter ::   masstol        = 1.e-6

  !----- Maximum height that a cloud can ever possibly reach [m] --------------------------!
  real, parameter ::   zmaxtpse       = 18000.

  !----------------------------------------------------------------------------------------!
  !    Ensemble related variables. acrit and acritt are a look-up table for climatological !
  ! cloud work function                                                                    !
  !----------------------------------------------------------------------------------------!
  !----- Maximum number of levels with climatological cloud work function -----------------!
  integer                 , parameter :: maxcrit=30
  !----- Levels in which cloud work function is available ---------------------------------!
  real, dimension(maxcrit), parameter :: pclim=                                            &
                         (/ 85000., 83750., 82500., 81250., 80000., 78750., 77500., 76250. &
                          , 75000., 73750., 72500., 71250., 70000., 68750., 67500., 66250. &
                          , 65000., 63750., 62500., 61250., 60000., 55000., 50000., 45000. &
                          , 40000., 35000., 30000., 25000., 20000., 15000.                /)
  !----- Cloud work function from the first dataset ---------------------------------------!
  real, dimension(maxcrit), parameter :: aclim1=                                           &
                         (/  63.23,  57.95,  53.90,  52.36,  44.50,  49.65,  50.00,  49.83 &
                          ,  55.30,  52.89,  60.80,  58.83,  66.40,  67.66,  70.70,  79.37 &
                          ,  75.00,  93.86, 108.0 , 111.0 , 130.0 , 152.0 , 221.0 , 315.0  &
                          , 368.0 , 410.0 , 525.5 , 766.3 ,1168.6 ,1685.1                 /)
  !----- Cloud work function from the second dataset --------------------------------------!
  real, dimension(maxcrit), parameter :: aclim2=(/                                         &
                            203.0 , 299.0 , 359.0 , 403.0 , 515.0 , 478.0 , 518.0 , 530.0  &
                          , 521.0 , 565.0 , 543.0 , 588.0 , 566.0 , 602.0 , 596.0 , 611.0  &
                          , 625.0 , 619.0 , 645.0 , 627.0 , 665.0 , 659.0 , 688.0 , 743.0  &
                          , 813.0 , 886.0 , 947.0 ,1138.0 ,1377.0 ,1896.0                 /)
  !----------------------------------------------------------------------------------------!

  contains
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine simply assign proper values for some Grell's cumulus parameter-      !
! ization. This will try to assign as little memory as needed and turn options on and off  !
! or assign proper values depending on the user settings.                                  !
!------------------------------------------------------------------------------------------!
   subroutine define_grell_coms(ngrids,nclouds,mmzp,nnqparm,grell_1st,grell_last)

      implicit none

      integer                   , intent(in) :: ngrids     ! Number of grids (nested)
      integer                   , intent(in) :: nclouds    ! Number of clouds
      integer, dimension(ngrids), intent(in) :: mmzp       ! Number of points in Z direction
      integer, dimension(ngrids), intent(in) :: nnqparm    ! Flag for running cumulus
      integer, dimension(ngrids), intent(in) :: grell_1st  ! First cloud for Grell to solve
      integer, dimension(ngrids), intent(in) :: grell_last ! Last cloud for Grell to solve
      integer                                :: icld       ! Cloud counter
      real                                   :: mycdf      ! User-based CDF
      real   , external                      :: cdf2normal ! Cumulative distribution func.
      !----- Getting the maximum number of layers among the grids -------------------------!
      mgmzp = maxval(mmzp)

      !------------------------------------------------------------------------------------!
      !     Now I check whether I am running an ensemble of a simple dynamic control, and  !
      ! assign maxens_dyn accordingly. Also, only Grell (1993) and Kain-Fritsch (1990)     !
      ! require the use of non-forced cloud work function, so I skip it if it is           !
      ! unecessary.                                                                        !
      !------------------------------------------------------------------------------------!
      cloudloop: do icld=1,nclouds


         !---------------------------------------------------------------------------------!
         !   First I decide whether this cloud will be allowed to have precipitation and   !
         ! downdrafts. If not, reduce the dimension of precipitation efficiency since it   !
         ! won't matter. For now precipitating clouds come together with downdrafts (both  !
         ! or nothing). And this decision is currently defined by a single number          !
         ! min_down_radius...                                                              !
         !---------------------------------------------------------------------------------!
         comp_down(icld) = radius(icld) > min_down_radius
         if (.not. comp_down(icld)) maxens_eff(icld) = 1


         !---------------------------------------------------------------------------------!
         !    Then I define the dynamic control dimension and the step bypasses. I will    !
         ! also check whether the cloud has a radius too small to have precipitation and   !
         ! downdrafts, in which case some ensembles will be forbidden.                     !
         !---------------------------------------------------------------------------------!
         select case (closure_type(icld))
         case ('en')
            if (.not. comp_down(icld)) then
               write(unit=*,fmt='(a)') '--------------------------------------------------'
               write(unit=*,fmt='(a,1x,i4)')   ' - For cloud #',icld
               write(unit=*,fmt='(a,1x,f8.2)') ' - Radius is ',radius(icld)
               write(unit=*,fmt='(a,1x,a)')    ' - Closure type is ',closure_type(icld)
               write(unit=*,fmt='(a,1x,f8.2)') ' - Minimum radius for this closure ',      &
                                                   min_down_radius
               call abort_run('Radius is too small for the chosen clousure type',          &
                              'define_grell_coms','grell_coms.f90')
            end if
            maxens_dyn(icld)          = 16
            comp_noforc_cldwork(icld) = .true.
            comp_modif_thermo(icld)   = .true.

         case ('nc')
            maxens_dyn(icld)          = 10
            comp_noforc_cldwork(icld) = .true.
            comp_modif_thermo(icld)   = .true.

         case ('as','kf','gr')
            maxens_dyn(icld)          = 1
            comp_noforc_cldwork(icld) = .true.
            comp_modif_thermo(icld)   = .true.

         case ('pp') !AS should be here, moved up for debugging only.
            maxens_dyn(icld)          = 1
            comp_noforc_cldwork(icld) = .false.
            comp_modif_thermo(icld)   = .true.

         case default
            if (.not. comp_down(icld)) then
               write(unit=*,fmt='(a)') '--------------------------------------------------'
               write(unit=*,fmt='(a,1x,i4)')   ' - For cloud #',icld
               write(unit=*,fmt='(a,1x,f8.2)') ' - Radius is ',radius(icld)
               write(unit=*,fmt='(a,1x,a)')    ' - Closure type is ',closure_type(icld)
               write(unit=*,fmt='(a,1x,f8.2)') ' - Minimum radius for this closure ',      &
                                                   min_down_radius
               call abort_run('Radius is too small for the chosen clousure type',          &
                              'define_grell_coms','grell_coms.f90')
            end if
            maxens_dyn(icld)          = 1
            comp_noforc_cldwork(icld) = .false.
            comp_modif_thermo(icld)   = .false.

         end select
         !---------------------------------------------------------------------------------!


         !----- Finding the inverse of ensemble dimension ---------------------------------!
         inv_ensdim(icld) = 1./ real(maxens_dyn(icld)*maxens_lsf(icld)                     &
                                    *maxens_eff(icld)*maxens_cap(icld))

         !----- Finding the maximum depth. For the deepest cloud, set no bounds -----------!
         if (icld /= 1) then
            depth_max(icld)=depth_min(icld-1)*1.00
         else
            depth_max(icld)=huge(1.)
         end if

         !---------------------------------------------------------------------------------!
         !  Configure the maximum normalized wind based on cap_maxs                        !
         !---------------------------------------------------------------------------------!
         if (cap_maxs(icld) < 0.) then
            mycdf = (100.+cap_maxs(icld))/100.
            wnorm_max(icld) = cdf2normal(mycdf)
         else
            cap_maxs(icld) = cap_maxs(icld) * 100.
            wnorm_max(icld) = 0.
         end if

      end do cloudloop
      !------------------------------------------------------------------------------------!
      !    Here we overwrite the ensemble dimensions to a very small number in case the    !
      ! user opted for another cumulus parameterization for those grids.                   !
      !------------------------------------------------------------------------------------!
      if (.not. any(grell_1st == 1)) then
         maxens_dyn(1) = 1
         maxens_eff(1) = 1
         maxens_lsf(1) = 1
         maxens_cap(1) = 1
         inv_ensdim(1) = 1.
      end if
      if (.not. any(grell_last == nclouds)) then
         maxens_dyn(nclouds) = 1
         maxens_eff(nclouds) = 1
         maxens_lsf(nclouds) = 1
         maxens_cap(nclouds) = 1
         inv_ensdim(nclouds) = 1.
      end if
      !------------------------------------------------------------------------------------!

      return
   end subroutine define_grell_coms
!==========================================================================================!
!==========================================================================================!
end module grell_coms
