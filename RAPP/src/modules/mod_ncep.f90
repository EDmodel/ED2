!==========================================================================================!
!==========================================================================================!
! Module ncep: This module contains the variables that we will need to save for NCEP input !
!              variables.                                                                  !
!------------------------------------------------------------------------------------------!
module mod_ncep
   use mod_maxdims , only : maxstr
   use mod_time    , only : time_stt

   !----- NCEP variable structure. --------------------------------------------------------!
   type ncep_vars
      type(time_stt), dimension(:)    , pointer :: when
      real          , dimension(:,:,:), pointer :: pres
      real          , dimension(:,:,:), pointer :: temp
      real          , dimension(:,:,:), pointer :: rhum
      real          , dimension(:,:,:), pointer :: uwnd
      real          , dimension(:,:,:), pointer :: vwnd
      real          , dimension(:,:,:), pointer :: prate
      real          , dimension(    :), pointer :: prate1d
      real          , dimension(:,:,:), pointer :: dlwrf
      real          , dimension(:,:,:), pointer :: nbdsf
      real          , dimension(:,:,:), pointer :: nddsf
      real          , dimension(:,:,:), pointer :: vbdsf
      real          , dimension(:,:,:), pointer :: vddsf
      real          , dimension(:,:,:), pointer :: shum
   end type ncep_vars
   
   !----- NCEP variable array, it may have more than one grid in the input. ---------------!
   type(ncep_vars), dimension(:), allocatable :: ncep_g
   !---------------------------------------------------------------------------------------!


   
   !---------------------------------------------------------------------------------------!
   !   Number of grids we will work with:                                                  !
   ! 1. Gaussian grid, with input data frequency + the output grid for state variables;    !
   ! 2. Gaussian grid, with higher time resolution, the output grid for flux variables;    !
   ! 3. Lon/lat grid, the input grid for state variables.                                  !
   !---------------------------------------------------------------------------------------!
   integer                          , parameter :: ngrids_ncep  = 3
   logical, dimension(ngrids_ncep+1), parameter :: flux_g  = (/  .true.,  .true.           &
                                                              , .false., .false. /)
   logical, dimension(ngrids_ncep+1), parameter :: state_g = (/ .true.,  .false.           &
                                                             ,  .true.,  .false. /)
   logical, dimension(ngrids_ncep+1), parameter :: rain_g  = (/ .true.,  .false.           &
                                                             , .false.,   .true. /)
   !---------------------------------------------------------------------------------------!

   
   !----- Number of variables. ------------------------------------------------------------!
   integer, parameter :: nvars_ncep = 12 ! # of input variables from NCEP
   integer, parameter :: nvars_ol1  = 10 ! # of input variables for "state" grid output
   integer, parameter :: nvars_ol2  =  6 ! # of input variables for flux grid output.
   integer, parameter :: nvars_ol4  =  3 ! # of input variables for flux grid output.
   !---------------------------------------------------------------------------------------!


   !----- Variable names ------------------------------------------------------------------!
   character(len=maxstr), parameter, dimension(nvars_ncep) :: vars_ncep =                  &
                       (/ 'air  '            & ! Air temperature                  [      K]
                        , 'pres '            & ! Pressure                         [     Pa]
                        , 'rhum '            & ! Relative humidity                [      %]
                        , 'uwnd '            & ! Zonal wind                       [    m/s]
                        , 'vwnd '            & ! Zonal wind                       [    m/s]
                        , 'pres '            & ! Pressure                         [     Pa]
                        , 'dlwrf'            & ! Downward long wave radiation     [   W/m2]
                        , 'nbdsf'            & ! Near-IR beam radiation           [   W/m2]
                        , 'nddsf'            & ! Near-IR diffuse radiation        [   W/m2]
                        , 'vbdsf'            & ! Visible beam radiation           [   W/m2]
                        , 'vddsf'            & ! Visible beam radiation           [   W/m2]
                        , 'prate'           /) ! Precipitation rate               [kg/m2/s]
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Prefix of files containing the variables, one file per variable, not with full    !
   ! path, but only what is after INPATH.  This    
   !---------------------------------------------------------------------------------------!
   character(len=maxstr), parameter, dimension(nvars_ncep) :: prefvars_ncep =              &
                (/ 'air.sig995/air.sig995          '                                       &
                 , 'pres.sfc/pres.sfc              '                                       &
                 , 'rhum.sig995/rhum.sig995        '                                       &
                 , 'uwnd.sig995/uwnd.sig995        '                                       &
                 , 'vwnd.sig995/vwnd.sig995        '                                       &
                 , 'pres.sfc.gauss/pres.sfc.gauss  '                                       &
                 , 'dlwrf.sfc.gauss/dlwrf.sfc.gauss'                                       &
                 , 'nbdsf.sfc.gauss/nbdsf.sfc.gauss'                                       &
                 , 'nddsf.sfc.gauss/nddsf.sfc.gauss'                                       &
                 , 'vbdsf.sfc.gauss/vbdsf.sfc.gauss'                                       &
                 , 'vddsf.sfc.gauss/vddsf.sfc.gauss'                                       &
                 , 'prate.sfc.gauss/prate.sfc.gauss'                                      /)  
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Grid where the variable should be loaded:                                         !
   ! + Gaussian variables go to grid 1 (the one that we will actually use).                !
   ! + Lon/lat variables go to grid 3  (the one to be interpolated).                       !
   !---------------------------------------------------------------------------------------!
   integer, parameter, dimension(nvars_ncep) :: grids_ncep = (/3,3,3,3,3,1,1,1,1,1,1,1/)
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !     Table with the output variables for both the state and the flux grids.            !
   !---------------------------------------------------------------------------------------!
   !----- State grid. ---------------------------------------------------------------------!
   character(len=maxstr), parameter, dimension(nvars_ol1) :: vars_ol1 =                    &
                       (/ 'lon  '            & ! Longitude                        [    deg]
                        , 'lat  '            & ! Latitude                         [    deg]
                        , 'hgt  '            & ! Reference height                 [  m AGL]
                        , 'tmp  '            & ! Air temperature                  [      K]
                        , 'pres '            & ! Pressure                         [     Pa]
                        , 'sh   '            & ! Specific humidity                [  kg/kg]
                        , 'ugrd '            & ! Zonal wind                       [    m/s]
                        , 'vgrd '            & ! Zonal wind                       [    m/s]
                        , 'prate'            & ! Precipitation rate               [kg/m2/s]
                        , 'dlwrf'           /) ! Downward long wave radiation     [   W/m2]
   !----- Flux grid. ----------------------------------------------------------------------!
   character(len=maxstr), parameter, dimension(nvars_ol2) :: vars_ol2 =                   &
                       (/ 'lon  '            & ! Longitude                        [    deg]
                        , 'lat  '            & ! Latitude                         [    deg]
                        , 'nbdsf'            & ! Near-IR beam radiation           [   W/m2]
                        , 'nddsf'            & ! Near-IR diffuse radiation        [   W/m2]
                        , 'vbdsf'            & ! Visible beam radiation           [   W/m2]
                        , 'vddsf'           /) ! Visible beam radiation           [   W/m2]
   !----- State grid. ---------------------------------------------------------------------!
   character(len=maxstr), parameter, dimension(nvars_ol4) :: vars_ol4 =                    &
                       (/ 'lon  '            & ! Longitude                        [    deg]
                        , 'lat  '            & ! Latitude                         [    deg]
                        , 'prate'           /) ! Precipitation rate               [kg/m2/s]
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !     Variable type.  Acceptable types include:                                         !
   ! 0. Read gridded data - no time interpolation;                                         !
   ! 1. Read gridded data - with time interpolatation;                                     !
   ! 2. Read gridded data - constant in time, not changing (if this is lat/lon, it will    !
   !    overwrite line 3 information);                                                     !
   ! 3. Read one value representing the whole grid - no time interpolation;                !
   ! 4. Specify a constant for all polygons, constant in time (most likely reference       !
   !    height).                                                                           !
   !---------------------------------------------------------------------------------------!
   integer, parameter, dimension(nvars_ol1) :: vtype_ol1=(/ 2, 2, 2, 1, 1, 1, 1, 1, 1, 0/)
   integer, parameter, dimension(nvars_ol2) :: vtype_ol2=(/ 2, 2, 0, 0, 0, 0/)
   integer, parameter, dimension(nvars_ol4) :: vtype_ol4=(/ 2, 2, 0 /)
   !=======================================================================================!
   !=======================================================================================!



   contains


   !=======================================================================================!
   !=======================================================================================!
   subroutine alloc_ncep(ncep,nx,ny,nt,state,flux,rain)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(ncep_vars), intent(inout) :: ncep
      integer        , intent(in)    :: nx
      integer        , intent(in)    :: ny
      integer        , intent(in)    :: nt
      logical        , intent(in)    :: state
      logical        , intent(in)    :: flux
      logical        , intent(in)    :: rain
      !------------------------------------------------------------------------------------!


      !------ Time is always allocated. ---------------------------------------------------!
      allocate(ncep%when(nt))

      if (state) then
         allocate (ncep%pres    (nx,ny,nt))
         allocate (ncep%temp    (nx,ny,nt))
         allocate (ncep%rhum    (nx,ny,nt))
         allocate (ncep%uwnd    (nx,ny,nt))
         allocate (ncep%vwnd    (nx,ny,nt))
         allocate (ncep%shum    (nx,ny,nt))
         allocate (ncep%dlwrf   (nx,ny,nt))
      end if

      if (flux) then
         allocate (ncep%nbdsf   (nx,ny,nt))
         allocate (ncep%nddsf   (nx,ny,nt))
         allocate (ncep%vbdsf   (nx,ny,nt))
         allocate (ncep%vddsf   (nx,ny,nt))
      end if
      
      if (rain) then
         allocate (ncep%prate   (nx,ny,nt))
         allocate (ncep%prate1d       (nt))
      end if

      return
   end subroutine alloc_ncep
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   subroutine nullify_ncep(ncep)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(ncep_vars), intent(inout) :: ncep
      !------------------------------------------------------------------------------------!
      if (associated(ncep%when    )) nullify (ncep%when    )
      if (associated(ncep%pres    )) nullify (ncep%pres    )
      if (associated(ncep%temp    )) nullify (ncep%temp    )
      if (associated(ncep%rhum    )) nullify (ncep%rhum    )
      if (associated(ncep%uwnd    )) nullify (ncep%uwnd    )
      if (associated(ncep%vwnd    )) nullify (ncep%vwnd    )
      if (associated(ncep%shum    )) nullify (ncep%shum    )
      if (associated(ncep%prate   )) nullify (ncep%prate   )
      if (associated(ncep%prate1d )) nullify (ncep%prate1d )
      if (associated(ncep%dlwrf   )) nullify (ncep%dlwrf   )
      if (associated(ncep%nbdsf   )) nullify (ncep%nbdsf   )
      if (associated(ncep%nddsf   )) nullify (ncep%nddsf   )
      if (associated(ncep%vbdsf   )) nullify (ncep%vbdsf   )
      if (associated(ncep%vddsf   )) nullify (ncep%vddsf   )

      return
   end subroutine nullify_ncep
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine will assign an initial value to the arrays.  The default value     !
   ! is the missing number flag.                                                           !
   !---------------------------------------------------------------------------------------!
   subroutine init_ncep(ncep,nt)
      use mod_time   , only : missflg_time ! ! subroutine
      use mod_ioopts , only : missflg_real & ! intent(in)
                            , missflg_int  & ! intent(in)
                            , missflg_char & ! intent(in)
                            , missflg_dble ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(ncep_vars), intent(inout) :: ncep
      integer        , intent(in)    :: nt
      !------------------------------------------------------------------------------------!
      if (associated(ncep%when    )) call missflg_time(nt,missflg_int,missflg_real         &
                                                      ,missflg_dble,missflg_char,ncep%when)
      if (associated(ncep%pres    )) ncep%pres    = missflg_real
      if (associated(ncep%temp    )) ncep%temp    = missflg_real
      if (associated(ncep%rhum    )) ncep%rhum    = missflg_real
      if (associated(ncep%uwnd    )) ncep%uwnd    = missflg_real
      if (associated(ncep%vwnd    )) ncep%vwnd    = missflg_real
      if (associated(ncep%shum    )) ncep%shum    = missflg_real
      if (associated(ncep%prate   )) ncep%prate   = missflg_real
      if (associated(ncep%prate1d )) ncep%prate1d = missflg_real
      if (associated(ncep%dlwrf   )) ncep%dlwrf   = missflg_real
      if (associated(ncep%nbdsf   )) ncep%nbdsf   = missflg_real
      if (associated(ncep%nddsf   )) ncep%nddsf   = missflg_real
      if (associated(ncep%vbdsf   )) ncep%vbdsf   = missflg_real
      if (associated(ncep%vddsf   )) ncep%vddsf   = missflg_real

      return
   end subroutine init_ncep
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   subroutine dealloc_ncep(ncep)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(ncep_vars), intent(inout) :: ncep
      !------------------------------------------------------------------------------------!
      if (associated(ncep%when    )) deallocate (ncep%when    )
      if (associated(ncep%pres    )) deallocate (ncep%pres    )
      if (associated(ncep%temp    )) deallocate (ncep%temp    )
      if (associated(ncep%rhum    )) deallocate (ncep%rhum    )
      if (associated(ncep%uwnd    )) deallocate (ncep%uwnd    )
      if (associated(ncep%vwnd    )) deallocate (ncep%vwnd    )
      if (associated(ncep%shum    )) deallocate (ncep%shum    )
      if (associated(ncep%prate   )) deallocate (ncep%prate   )
      if (associated(ncep%prate1d )) deallocate (ncep%prate1d )
      if (associated(ncep%dlwrf   )) deallocate (ncep%dlwrf   )
      if (associated(ncep%nbdsf   )) deallocate (ncep%nbdsf   )
      if (associated(ncep%nddsf   )) deallocate (ncep%nddsf   )
      if (associated(ncep%vbdsf   )) deallocate (ncep%vbdsf   )
      if (associated(ncep%vddsf   )) deallocate (ncep%vddsf   )

      return
   end subroutine dealloc_ncep
   !=======================================================================================!
   !=======================================================================================!
end module mod_ncep
!==========================================================================================!
!==========================================================================================!
