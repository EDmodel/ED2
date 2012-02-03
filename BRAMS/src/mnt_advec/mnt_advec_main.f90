!==========================================================================================!
!==========================================================================================!
!    Subroutine radvc_mnt.f90.  This subroutine is the main driver for the monotic         !!
! advection scheme.  This was originally implemented by Saulo R. Freitas                   !
! saulo.freitas@cptec.inpe.br) in June 2009, and it was originally parallelised by         !
! Luiz Flavio and Jairo Panneta.  This advection scheme is highly conservative, monotonic, !
! and it keeps the mass mixing ratio positive.                                             !
!                                                                                          !
!    The current version was adapted to BRAMS-4.0.6 by Marcos Longo in December 2011.  In  !
! this version we don't use the additional memory, instead we try to follow the original   !
! wrapper functions.                                                                       !
!                                                                                          !
!    References:                                                                           !
!                                                                                          !
!    Walcek, C. J., N. M. Aleksic, 1998: A simple but accurate mass conservative, peak-    !
!        -preserving, mixing ratio bounded advection algorithm with Fortran code.  Atmos.  !
!        Environ., 32, 3863-3880.                                                          !
!    Walcek, C. J., 2000: Minor flux adjustment near mixing ratio extremes for simplified  !
!        yet highly accurate monotonic calculation of tracer advection.  J. Geophys. Res., !
!        105(D7), 9335-9348.                                                               !
!------------------------------------------------------------------------------------------!
subroutine radvc_mnt_driver(m1,m2,m3,ia,iz,ja,jz,mynum)

   use grid_dims    , only : maxgrds        ! ! intent(in)
   use mem_grid     , only : ngrid          & ! intent(in)
                           , grid_g         & ! intent(in)
                           , dzt            & ! intent(in)
                           , dtlt           ! ! intent(in)
   use mem_basic    , only : basic_g        ! ! intent(in)
   use mem_mnt_advec, only : advec_g        ! ! intent(in)
   use var_tables   , only : num_scalar     & ! intent(in)
                           , scalar_tab     ! ! intent(in)
   use mem_scratch  , only : vctr11         & ! Scratch for scal_in
                           , vctr12         & ! Scratch for uavg, vavg, wavg
                           , vctr13         & ! Scratch for den(i-1)_wal
                           , vctr14         & ! Scratch for den(i)_wal
                           , vctr15         & ! Scratch for  densu, densv, densw
                           , vctr16         & ! Scratch for dxtw
                           , vctr17         ! ! Scratch for scal_out
   implicit none

   !----- Arguments. ----------------------------------------------------------------------!
   integer                    , intent(in) :: m1
   integer                    , intent(in) :: m2
   integer                    , intent(in) :: m3
   integer                    , intent(in) :: ia
   integer                    , intent(in) :: iz
   integer                    , intent(in) :: ja
   integer                    , intent(in) :: jz
   integer                    , intent(in) :: mynum
   !----- Local variables. ----------------------------------------------------------------!
   integer                                 :: nv
   integer                                 :: ka
   integer                                 :: kz
   integer                                 :: k
   integer                                 :: i
   integer                                 :: j
   integer                                 :: mzxyp
   real   , dimension(:)      , pointer    :: scalarp
   real   , dimension(:)      , pointer    :: scalart
   !----- Locally saved variables. --------------------------------------------------------!
   logical, dimension(maxgrds), save       :: first_time = .true.
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     If this is the first time we call this subroutine, we must initialise the grid    !
   ! spacing.                                                                              !
   !---------------------------------------------------------------------------------------!
   if (first_time(ngrid)) then
      call init_grid_spacing( m1, m2, m3, grid_g (ngrid)%dxt  , grid_g (ngrid)%dyt         &
                                        , dzt                 , grid_g (ngrid)%fmapt       &
                                        , grid_g (ngrid)%rtgt , advec_g(ngrid)%dxtw        &
                                        , advec_g(ngrid)%dytw , advec_g(ngrid)%dztw        )
      first_time(ngrid) = .false.
   end if
   !---------------------------------------------------------------------------------------!



   !----- Alias for number of volume points in this node's sub-domain. --------------------!
   mzxyp = m1 * m2 * m3
   !---------------------------------------------------------------------------------------!



   !----- Alias for first and last levels to be computed. ---------------------------------!
   ka = 2
   kz = m1-1
   !---------------------------------------------------------------------------------------!



   !----- Find actual air densities. ------------------------------------------------------!
   call find_actual_densities( m1, m2, m3, basic_g(ngrid)%rtp   , basic_g(ngrid)%rv        &
                                         , basic_g(ngrid)%pp    , basic_g(ngrid)%pi0       &
                                         , basic_g(ngrid)%theta , advec_g(ngrid)%denst     &
                                         , advec_g(ngrid)%densu , advec_g(ngrid)%densv     &
                                         , advec_g(ngrid)%densw )
   !---------------------------------------------------------------------------------------!



   !----- Find the winds that will be used by the advection scheme. -----------------------!
   call find_avg_winds( m1, m2, m3, ia, iz, ja, jz, ka, kz                                 &
                      , basic_g(ngrid)%uc    , basic_g(ngrid)%up                           &
                      , basic_g(ngrid)%vc    , basic_g(ngrid)%vp                           &
                      , basic_g(ngrid)%wc    , basic_g(ngrid)%wp                           &
                      , grid_g (ngrid)%fmapui, grid_g (ngrid)%fmapvi                       &
                      , grid_g (ngrid)%rtgt  , grid_g (ngrid)%rtgu                         &
                      , grid_g (ngrid)%rtgv  , grid_g (ngrid)%f13t                         &
                      , grid_g (ngrid)%f23t  , advec_g(ngrid)%uavg                         &
                      , advec_g(ngrid)%vavg  , advec_g(ngrid)%wavg                         )
   !---------------------------------------------------------------------------------------!



   !----- Find the Walcek's density terms. ------------------------------------------------!
   call find_walcek_densities( dtlt, m1, m2, m3                                            &
                             , advec_g(ngrid)%uavg    , advec_g(ngrid)%vavg                &
                             , advec_g(ngrid)%wavg    , advec_g(ngrid)%denst               &
                             , advec_g(ngrid)%densu   , advec_g(ngrid)%densv               &
                             , advec_g(ngrid)%densw   , advec_g(ngrid)%den0_wal            &
                             , advec_g(ngrid)%den1_wal, advec_g(ngrid)%den2_wal            &
                             , advec_g(ngrid)%den3_wal, advec_g(ngrid)%dxtw                &
                             , advec_g(ngrid)%dytw    , advec_g(ngrid)%dztw    )
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Exchange the boundary conditions for most fields.                                 !
   !---------------------------------------------------------------------------------------!
   call mpilbc_driver('fulladv',0)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Loop over all scalars.                                                            !
   !---------------------------------------------------------------------------------------!
   do nv=1,num_scalar(ngrid)
      !----- Use some pointers to make it easier to read. ---------------------------------!
      scalarp => scalar_tab(nv,ngrid)%var_p
      scalart => scalar_tab(nv,ngrid)%var_t
      !------------------------------------------------------------------------------------!



      !----- Copy the scalar to the advection scratch array and update boundaries. --------!
      call atob(mzxyp,scalarp,advec_g(ngrid)%scal_in )
      call mpilbc_driver('fulladv',5)
      call atob(mzxyp,advec_g(ngrid)%scal_in,advec_g(ngrid)%scal_out)
      !------------------------------------------------------------------------------------!



      !----- Make the advection for the X direction. --------------------------------------!
      do j=ja,jz
         do k=2,m1-1
            !----- Copy vectors to the scratch vectors. -----------------------------------!
            call array2xcol(m1,m2,m3,k,j,advec_g(ngrid)%scal_in ,vctr11)
            call array2xcol(m1,m2,m3,k,j,advec_g(ngrid)%uavg    ,vctr12)
            call array2xcol(m1,m2,m3,k,j,advec_g(ngrid)%den0_wal,vctr13)
            call array2xcol(m1,m2,m3,k,j,advec_g(ngrid)%den1_wal,vctr14)
            call array2xcol(m1,m2,m3,k,j,advec_g(ngrid)%densu   ,vctr15)
            call array2xcol(m1,m2,m3,k,j,advec_g(ngrid)%dxtw    ,vctr16)
            !------------------------------------------------------------------------------!



            !----- Solve the monotonic advection. -----------------------------------------!
            call monotonic_advec(m2,ia,iz,dtlt,vctr11,vctr12,vctr13,vctr14,vctr15,vctr16   &
                                ,vctr17)
            !------------------------------------------------------------------------------!



            !----- Copy solution back to the output array. --------------------------------!
            call xcol2array(m1,m2,m3,k,j,vctr17,advec_g(ngrid)%scal_out)
            !------------------------------------------------------------------------------!
         end do
      end do
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Copy the output to the input and exchange B.C. before we run the advection for !
      ! the Y direction.                                                                   !
      !------------------------------------------------------------------------------------!
      call atob(mzxyp,advec_g(ngrid)%scal_out,advec_g(ngrid)%scal_in)
      call mpilbc_driver('fulladv',5)
      call atob(mzxyp,advec_g(ngrid)%scal_in,advec_g(ngrid)%scal_out)
      !------------------------------------------------------------------------------------!



      !----- Make the advection for the Y direction. --------------------------------------!
      do i=ia,iz
         do k=2,m1-1
            !----- Copy vectors to the scratch vectors. -----------------------------------!
            call array2ycol(m1,m2,m3,k,i,advec_g(ngrid)%scal_in ,vctr11)
            call array2ycol(m1,m2,m3,k,i,advec_g(ngrid)%vavg    ,vctr12)
            call array2ycol(m1,m2,m3,k,i,advec_g(ngrid)%den1_wal,vctr13)
            call array2ycol(m1,m2,m3,k,i,advec_g(ngrid)%den2_wal,vctr14)
            call array2ycol(m1,m2,m3,k,i,advec_g(ngrid)%densv   ,vctr15)
            call array2ycol(m1,m2,m3,k,i,advec_g(ngrid)%dytw    ,vctr16)
            !------------------------------------------------------------------------------!



            !----- Solve the monotonic advection. -----------------------------------------!
            call monotonic_advec(m3,ja,jz,dtlt,vctr11,vctr12,vctr13,vctr14,vctr15,vctr16   &
                                ,vctr17)
            !------------------------------------------------------------------------------!



            !----- Copy solution back to the output array. --------------------------------!
            call ycol2array(m1,m2,m3,k,i,vctr17,advec_g(ngrid)%scal_out)
            !------------------------------------------------------------------------------!
         end do
      end do
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Copy the output to the input and exchange B.C. before we run the advection for !
      ! the Z direction.                                                                   !
      !------------------------------------------------------------------------------------!
      call atob(mzxyp,advec_g(ngrid)%scal_out,advec_g(ngrid)%scal_in)
      call mpilbc_driver('fulladv',5)
      call atob(mzxyp,advec_g(ngrid)%scal_in,advec_g(ngrid)%scal_out)
      !------------------------------------------------------------------------------------!



      !----- Make the advection for the Y direction. --------------------------------------!
      do j=ja,jz
         do i=ia,iz
            !----- Copy vectors to the scratch vectors. -----------------------------------!
            call array2zcol(m1,m2,m3,i,j,advec_g(ngrid)%scal_in ,vctr11)
            call array2zcol(m1,m2,m3,i,j,advec_g(ngrid)%wavg    ,vctr12)
            call array2zcol(m1,m2,m3,i,j,advec_g(ngrid)%den2_wal,vctr13)
            call array2zcol(m1,m2,m3,i,j,advec_g(ngrid)%den3_wal,vctr14)
            call array2zcol(m1,m2,m3,i,j,advec_g(ngrid)%densw   ,vctr15)
            call array2zcol(m1,m2,m3,i,j,advec_g(ngrid)%dztw    ,vctr16)
            !------------------------------------------------------------------------------!



            !----- Solve the monotonic advection. -----------------------------------------!
            call monotonic_advec(m1,ka,kz,dtlt,vctr11,vctr12,vctr13,vctr14,vctr15,vctr16   &
                                ,vctr17)
            !------------------------------------------------------------------------------!



            !----- Copy solution back to the output array. --------------------------------!
            call zcol2array(m1,m2,m3,i,j,vctr17,advec_g(ngrid)%scal_out)
            !------------------------------------------------------------------------------!
         end do
      end do
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Update the tendency due to advection.                                         !
      !------------------------------------------------------------------------------------!
      call advtndc(m1,m2,m3,ia,iz,ja,jz,scalarp,advec_g(ngrid)%scal_out,scalart,dtlt,mynum)
      !------------------------------------------------------------------------------------!
   end do
   !---------------------------------------------------------------------------------------!

   return
end subroutine radvc_mnt_driver
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine calculates change in the analysed property qp during time step       !
! dtime due to advection along a grid that has npts points in a single direction. The      !
! property is loaded into the qp array, and the result is updated to qc.  Velocities       !
! "wind" and fluxes "flux" are specified at the staggered grid points.  ddim0 is the       !
! dimensionally dependent density before the advection of this term, whilst ddimm1 is the  !
! dimensionally dependent density including this direction. densp is the air density       !
! before the any advection.                                                                !
!                                                                                          !
! 1D grid   ->    |  na     |  na+1   |    ...     |  nz-1   |  nz     |                   !
! Wind      -> u(na-1)   u(na)     u(na+1) ...  u(nz-1)   u(nz)     u(nz+1)                !
! Property  ->    | q(na)   | q(na+1) |    ...     | q(nz-1) | q(nz)   |                   !
! Density   -> d(na-1)   d(na)     d(na+1) ...  d(nz-1)   d(nz)     d(nz+1)                !
!                                                                                          !
!     Boundary conditions for the Q arrays are stored at cells na-1 and nz+1.              !
!                                                                                          !
!     For this subroutine and comments, we use the some generic notation.  Their actual    !
! meaning depends on the direction, and this table is the reference:                       !
!                                                                                          !
!    /---------------------------------------------------------------------------------\   !
!    |  Direction  |  Wind  |  npts  |    na  |    nz  |  delta_n | Left    |  Right   |   !
!    |-------------+--------+--------+--------+--------+----------+---------+----------|   !
!    |          X  |  uavg  |    m2  |    ia  |    iz  |  dxtw    | West    |  East    |   !
!    |          Y  |  vavg  |    m3  |    ja  |    jz  |  dytw    | South   |  North   |   !
!    |          Z  |  wavg  |    m1  |     2  |  m1-1  |  dztw    | Nadir   |  Zenith  |   !
!    \---------------------------------------------------------------------------------/   !
!                                                                                          !
!------------------------------------------------------------------------------------------!
subroutine monotonic_advec(npts,na,nz,dtime,qp,wind,ddm1,ddp0,densn,delta_n,qc)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                 , intent(in)    :: npts
   integer                 , intent(in)    :: na
   integer                 , intent(in)    :: nz
   real                    , intent(in)    :: dtime
   real   , dimension(npts), intent(in)    :: qp
   real   , dimension(npts), intent(in)    :: wind
   real   , dimension(npts), intent(in)    :: ddm1
   real   , dimension(npts), intent(in)    :: ddp0
   real   , dimension(npts), intent(in)    :: densn
   real   , dimension(npts), intent(in)    :: delta_n
   real   , dimension(npts), intent(inout) :: qc
   !----- Local variables. ----------------------------------------------------------------!
   integer                                 :: n
   integer                                 :: nam1
   integer                                 :: nzp1
   integer                                 :: nm1
   integer                                 :: np1
   logical, dimension(npts)                :: extreme
   logical                                 :: locmax
   logical                                 :: locmin
   real   , dimension(npts)                :: flux
   real   , dimension(npts)                :: qcmin
   real   , dimension(npts)                :: qcmax
   real                                    :: qleft
   real                                    :: qright
   real                                    :: qhalf
   real                                    :: qguess
   real                                    :: courant
   real                                    :: alpha
   !----- Local constants. ----------------------------------------------------------------!
   real                    , parameter     :: eps = epsilon(1.)
   !---------------------------------------------------------------------------------------!


   !----- Edge indices for winds. ---------------------------------------------------------!
   nam1 = na - 1
   nzp1 = nz + 1
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Initialise the extreme indices as false so only the local extremes will be true.  !
   !---------------------------------------------------------------------------------------!
   extreme(:) = .false.
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Initialise the "current" q with "previous" q so the indices outside the na:nz     !
   ! range will have something.                                                            !
   !---------------------------------------------------------------------------------------!
   qc(:) = qp(:)
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     Flag every local extreme.                                                         !
   !---------------------------------------------------------------------------------------!
   do n=na,nz
      !----- Auxiliary indices for neighbour cells. ---------------------------------------!
      nm1 = n - 1
      np1 = n + 1
      !------------------------------------------------------------------------------------!


      !----- Flag the point as an extreme if it is one. -----------------------------------!
      locmax     = qp(n) >= ( max(qp(nm1),qp(np1)) - eps )
      locmin     = qp(n) <= ( max(qp(nm1),qp(np1)) + eps )
      extreme(n) = locmax .or. locmin
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     qcmin and qcmax are the absolute physical minimum limits to the property       !
      ! qc ar time t + dtime.  If these limits are ever violated, then a non-monotonic     !
      ! (i.e. oscillatory) behaviour will happen.                                          !
      !------------------------------------------------------------------------------------!
      if (wind(nm1) >= 0.) then
         qleft  = qp(nm1)
      else
         qleft  = qp(n)
      end if
      if (wind(n)   <  0.) then
         qright = qp(np1)
      else
         qright = qp(n)
      end if
      qcmin(n) = min(qp(n),qleft,qright)
      qcmax(n) = max(qp(n),qleft,qright)
      !------------------------------------------------------------------------------------!
   end do
   !---------------------------------------------------------------------------------------!



   !>->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->!
   !>->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->!
   !     We first solve the "left-to-right" flux.                                          !
   !>->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->!
   !>->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->!
   !------ Leftmost boundary condition. ---------------------------------------------------!
   if (wind(nam1) >= 0.0) flux(nam1) = qp(nam1) * wind(nam1) * dtime * densn(nam1)
   !------ Update values for those points that are experiencing westerly advection. -------!
   lrloop: do n=na,nz
      nm1 = n - 1
      np1 = n + 1
      !----- Solve this point only if the winds are coming from the west. -----------------!
      if (wind(n) >= 0.0) then
         !---------------------------------------------------------------------------------!
         !    Check whether there is only outflow from this grid point, or there is        !
         ! inflow.                                                                         !
         !---------------------------------------------------------------------------------!
         if (wind(nm1) < 0.0) then
            !----- Outflow only. ----------------------------------------------------------!
            flux(n) = qp(n) * wind(n) * dtime * densn(n)
            !------------------------------------------------------------------------------!
         else
            !----- Outflow and inflow, check stability. -----------------------------------!



            !----- Find the Courant number. -----------------------------------------------!
            courant  = dtime * wind(n) / delta_n(n)
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !     Decide which alpha to use.                                               !
            !------------------------------------------------------------------------------!
            if (extreme(np1)) then
               !----- The cell downwind is a local extreme. -------------------------------!
               alpha = 1.75 - 0.45 * courant
            elseif (extreme(nm1)) then
               !----- The cell two cells upwind is a local extreme. -----------------------!
               alpha = max(1.5, 1.2 + 0.6 * courant)
            else
               !----- Default is 1.0. -----------------------------------------------------!
               alpha = 1.0
            end if
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !     Find the q(i+1/2) term.                                                  !
            !------------------------------------------------------------------------------!
            qhalf = qp(n) + 0.25 * (qp(np1) - qp(nm1)) * (1. - courant) * alpha
            qhalf = min(max(qhalf,min(qp(n),qp(np1))),max(qp(n),qp(np1)))
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !     Find the first guess for this qc.                                        !
            !------------------------------------------------------------------------------!
            qguess = (qp(n) * ddm1(n) - courant * qhalf * densn(n) + flux(nm1)/delta_n(n)) &
                   / ddp0(n)
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !      Final bounded answer and the flux.                                      !
            !------------------------------------------------------------------------------!
            qc(n)   = max(qcmin(n),min(qcmax(n),qguess))
            flux(n) = delta_n(n) * (qp(n)*ddm1(n) - qc(n)*ddp0(n)) + flux(nm1)
            !------------------------------------------------------------------------------!
         end if
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!
   end do lrloop
   !>->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->!
   !>->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->!




   !<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<!
   !<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<!
   !     Now we solve the "right-to-left" advection.                                       !
   !<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<!
   !<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<!
   !------ Rightmost boundary condition. --------------------------------------------------!
   if (wind(nz) < 0.) flux(nz) = qp(nzp1) * wind(nz) * dtime * densn(nz)
   !------ Update values for those points that are experiencing easterly advection. -------!
   rlloop: do n=nz,na,-1
      nm1 = n - 1
      np1 = n + 1
      if (wind(nm1) >= 0.) then
         if (wind(n) < 0.) then
            !------ Inflow-only cell. -----------------------------------------------------!
            qguess = (qp(n)*ddm1(n) - flux(n)/delta_n(n) + flux(nm1)/delta_n(n)) / ddp0(n)
            qc(n)  = max(qcmin(n),min(qcmax(n),qguess))
         end if
      else
         !------ Outflow and inflow, check stability. -------------------------------------!

         !----- Find the Courant number. --------------------------------------------------!
         courant  = dtime * abs(wind(nm1)) / delta_n(n)
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Decide which alpha to use.                                                  !
         !---------------------------------------------------------------------------------!
         if (extreme(nm1)) then
            !----- The cell downwind is a local extreme. ----------------------------------!
            alpha = 1.75 - 0.45 * courant
         elseif (extreme(np1)) then
            !----- The cell two cells upwind is a local extreme. --------------------------!
            alpha = max(1.5, 1.2 + 0.6 * courant)
         else
            !----- Default is 1.0. --------------------------------------------------------!
            alpha = 1.0
         end if
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !     Find the q(i-1/2) term.                                                     !
         !---------------------------------------------------------------------------------!
         qhalf = qp(n) + 0.25 * (qp(nm1) - qp(np1)) * (1. - courant) * alpha
         qhalf = min(max(qhalf,min(qp(n),qp(nm1))),max(qp(n),qp(nm1)))
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !     Find the first guess for this qc.                                           !
         !---------------------------------------------------------------------------------!
         qguess = (qp(n)*ddm1(n) - flux(n)/delta_n(n) - courant*qhalf*densn(nm1)) / ddp0(n)
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !      Final bounded answer and the flux.                                         !
         !---------------------------------------------------------------------------------!
         qc(n)     = max(qcmin(n),min(qcmax(n),qguess))
         flux(nm1) = delta_n(n) * (qc(n)*ddp0(n) - qp(n)*ddp0(n)) + flux(n)
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!
   end do rlloop
   !<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<!
   !<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<!

   return
end subroutine monotonic_advec
!==========================================================================================!
!==========================================================================================!
