!==========================================================================================!
! rconv_driver.f90                                                                         !
!                                                                                          !
!    This file contains the main convective parameterization driver. This is actually just !
! a wrapper to decide whether to call Grell, Kuo, or Souza parameterizations for this      !
! grid.                                                                                    !
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine rconv_driver()
   use mem_cuparm , only : confrq            & ! How often I should call this cloud
                         , cptime            & ! Time to start calling this cloud
                         , nclouds           & ! Number of clouds to use here 
                         , nnqparm           & ! Flag for cumulus parameterization
                         , ndeepest          & ! Flag for deepest cloud
                         , nshallowest       & ! Flag for shallowest cloud
                         , grell_1st         & ! First cloud to use Grell's method
                         , grell_last        & ! Last cloud to use Grell's scheme
                         , cuparm_g          & ! The cumulus scheme structure
                         , initialize_cuparm ! ! Routine that resets structure to zero
          
   use mem_grid   , only : initial      & ! Flag for "initial" run
                         , ngrid        & ! Current grid
                         , ngrids       & ! Total # of grids
                         , time         & ! Current time
                         , dtlt         & ! Current grid delta-t
                         , dtlongn      & ! Delta-t array
                         , grid_g       ! ! Grid structure

   use node_mod   , only : mynum        & ! This node number
                         , mzp          & ! Current grid # of vertical levels for this node
                         , myp          & ! Current grid # of meridional pts for this node
                         , mxp          & ! Current grid # of zonal points for this node
                         , ia           & ! Westernmost boundary
                         , iz           & ! Easternmost boundary
                         , ja           & ! Southernmost boundary
                         , jz           ! ! Northermost boundary
   
   use mem_turb   , only : turb_g       & ! Turbulence structure
                         , idiffk       ! ! Turbulence closure flag
   
   use mem_scratch, only : scratch      ! ! Scratch structure
   use mem_tend   , only : tend_g       ! ! Tendency structure
   use mem_basic  , only : co2_on       & ! Flag, whether CO2 is prognosed in this run.
                         , basic_g      ! ! Basic variable structure

   implicit none
   !----- Local variables. ----------------------------------------------------------------!
   integer             :: icld              ! Cloud counter
   integer             :: i                 ! x-position 
   integer             :: j                 ! y-position
   integer             :: k                 ! z-position
   integer             :: ifm               ! Grid counter
   integer             :: iun               ! File unit for debugging
   logical             :: cumulus_now       ! Is now a good time to call cumulus? [T|F]
   !----- Local constants. ----------------------------------------------------------------!
   logical          , parameter      :: print_debug=.false. ! Print debugging stuff [T|F]
   character(len=9) , parameter      :: fmti='(a,1x,i6)'
   character(len=13), parameter      :: fmtf='(a,1x,es14.7)'   
   !----- External functions. -------------------------------------------------------------!
   logical, external   :: cumulus_time      ! Function to check whether it's cumulus time.
   !---------------------------------------------------------------------------------------!

   !----- Nothing to be done in here for this grid ----------------------------------------!
   if (nnqparm(ngrid) == 0) return

   cumulus_now = cumulus_time(initial,time,cptime,confrq,dtlt)
   
   
   !---------------------------------------------------------------------------------------!
   !     Initialise all cumulus variables. Before we do that, though, we must copy aconpr  !
   ! to a scratch array and copy back.                                                     !
   !---------------------------------------------------------------------------------------!
   if (cumulus_now) then
      call atob(mxp*myp,cuparm_g(ngrid)%aconpr,scratch%vt2de)
      call initialize_cuparm(cuparm_g(ngrid))
      call atob(mxp*myp,scratch%vt2de,cuparm_g(ngrid)%aconpr)
   end if


   !---------------------------------------------------------------------------------------!
   !   Checking whether I should use Souza's or old Grell's scheme                         !
   !---------------------------------------------------------------------------------------!
   if ((nshallowest(ngrid) == 1 .or. nshallowest(ngrid) == 3) .and. cumulus_now) then
      select case(nshallowest(ngrid))
      case (1) !
         call souza_cupar_driver()
      case (3)
         call old_grell_cupar_driver(nclouds)
      end select
      
      !----- If I have TKE available, I may use it to define PBL top ----------------------!
      select case (idiffk(ngrid))
      case (1,4,5,6,7,8)
         call atob(mzp*mxp*myp,turb_g(ngrid)%tkep,scratch%vt3de)
      case default
         call azero(mzp*mxp*myp,scratch%vt3de)
      end select
   end if


   !---------------------------------------------------------------------------------------!
   !    Solve the cloud sizes that should be determined by Grell's scheme, if any of such  !
   ! exists and if this is the time to call Grell's parametrisation.                       !
   !---------------------------------------------------------------------------------------!
   if (grell_1st(ngrid) <= grell_last(ngrid) .and. cumulus_now) then
      call grell_cupar_driver(grell_1st(ngrid),grell_last(ngrid))
   end if

   !---------------------------------------------------------------------------------------!
   !   Checking whether I was supposed to run Kuo for this grid.                           !
   !---------------------------------------------------------------------------------------!
   if ((ndeepest(ngrid) == 1 .or. ndeepest(ngrid) == 3) .and. cumulus_now) then
       select case (ndeepest(ngrid))
       case (1)
          call kuo_cupar_driver()
       case (3)
          call old_grell_cupar_driver(1)
       end select


       !----- If I have TKE available, I may use it to define PBL top ---------------------!
       select case (idiffk(ngrid))
       case (1,4,5,6,7,8)
          call atob(mzp*mxp*myp,turb_g(ngrid)%tkep,scratch%vt3de)
       case default
          call azero(mzp*mxp*myp,scratch%vt3de)
       end select
   end if

   !---------------------------------------------------------------------------------------!
   !   Update the heating, moistening, and CO2 exchange rates, and the convective          !
   ! precipitation.                                                                        !
   !---------------------------------------------------------------------------------------!
   do icld=1,nclouds
      do j=ja,jz
         do i=ia,iz
            !----- Update the tendencies. -------------------------------------------------!
            do k=1,mzp
               tend_g(ngrid)%tht(k,i,j) = tend_g(ngrid)%tht(k,i,j)                         &
                                        + cuparm_g(ngrid)%thsrc(k,i,j,icld)
               tend_g(ngrid)%rtt(k,i,j) = tend_g(ngrid)%rtt(k,i,j)                         &
                                        + cuparm_g(ngrid)%rtsrc(k,i,j,icld)
               if (co2_on) then
                  tend_g(ngrid)%co2t(k,i,j) = tend_g(ngrid)%co2t(k,i,j)                    &
                                            + cuparm_g(ngrid)%co2src(k,i,j,icld)
               end if
            end do
            !----- Update the total precipitation. ----------------------------------------!
            cuparm_g(ngrid)%aconpr(i,j) = cuparm_g(ngrid)%aconpr(i,j)                      &
                                        + cuparm_g(ngrid)%conprr(i,j,icld) * dtlt
            !------------------------------------------------------------------------------!
         end do
      end do
   end do

   if (print_debug) then
      iun=20+mynum
      do icld=1,nclouds
         do j=1,myp
            do i=1,mxp
               write(unit=iun,fmt='(a)') '------------------------------------------------'
               write(unit=iun,fmt=fmti)  'I    = ',i
               write(unit=iun,fmt=fmti)  'J    = ',j
               write(unit=iun,fmt=fmti)  'ICLD = ',icld
               write(unit=iun,fmt='(4(a,1x))') 'LEVEL','         THP','         THT'       &
                                                      ,'       THSRC'
               do k=1,mzp
                  write (unit=iun,fmt='(i5,(3(1x,es12.5)))')                               &
                       k,basic_g(ngrid)%thp(k,i,j),cuparm_g(ngrid)%thsrc(k,i,j,icld)       &
                        ,tend_g(ngrid)%tht(k,i,j)
               end do
               write(unit=iun,fmt='(a)') '------------------------------------------------'
               write(unit=iun,fmt='(a)') ' '
            end do
         end do
      end do
   end if



   return
end subroutine rconv_driver
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!   This function simply tests whether this is time to run the cumulus parameterization.   !
!------------------------------------------------------------------------------------------!
logical function cumulus_time(initial,time,cptime,confrq,deltat)
   implicit none
   integer    , intent(in) :: initial ! Flag to state whether this is initial run
   real(kind=8),intent(in) :: time    ! Elapsed time
   real(kind=8),intent(in) :: cptime  ! Time to start computing cumulus
   real        ,intent(in) :: confrq  ! How often should this cloud be called?
   real        ,intent(in) :: deltat  ! model time step

 
   cumulus_time = (.not. (initial == 2 .and. time < cptime - dble(deltat))) .and. &
                  mod(time,confrq) < deltat 

   return
end function cumulus_time
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine will create a cloud based on the moistening and precipitation rate   !
! for Kuo and Souza cumulus parameterization. Ideally, this variable should come directly  !
! from the cumulus parameterization, this is a sketch just to not leave it empty.          !
!------------------------------------------------------------------------------------------!
subroutine cloud_sketch(m1,m2,m3,ia,iz,ja,jz,deltat,flpw,rtgt,kpbl,tke,upmf,rtsrc,conprr   &
                       ,cupcond)

   use mem_grid, only: zt,zm
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                     , intent(in)  :: m1,m2,m3    ! Matrix dimensions
   integer                     , intent(in)  :: ia,iz,ja,jz ! Node boundaries
   real                        , intent(in)  :: deltat      ! Time step
   real   , dimension   (m2,m3), intent(in)  :: flpw        ! Lowest point above ground
   real   , dimension   (m2,m3), intent(in)  :: rtgt        ! Correction
   integer, dimension   (m2,m3), intent(in)  :: kpbl        ! PBL top
   real   , dimension(m1,m2,m3), intent(in)  :: tke         ! TKE
   real   , dimension(m1,m2,m3), intent(in)  :: rtsrc       ! Moistening rate
   real   , dimension   (m2,m3), intent(in)  :: conprr      ! Precipitation rate
   real   , dimension   (m2,m3), intent(in)  :: upmf        ! Updraft reference 
   real   , dimension(m1,m2,m3), intent(out) :: cupcond     ! My Frankeinstein cloud
   !----- Local variables. ----------------------------------------------------------------!
   integer                                   :: i,j,k,lpw,kcldbase,kzdown
   real                                      :: dnmf
   !----- Local constants. ----------------------------------------------------------------!
   real, parameter :: zcldbase=1600. ! Fixed height for cloud base
   real, parameter :: zdown=3000.    ! Fixed height for downdraft origin
   real, parameter :: epsilon= 0.3   ! Just an arbitrary epsilon, dnmf/upmf
   real, parameter :: dnmfdef= 0.1   ! Just an arbitrary defaul dnmf
   !---------------------------------------------------------------------------------------!

   call azero(m1*m2*m3,cupcond)

   do j=ja,jz
      do i=ia,iz
         lpw=nint(flpw(i,j))
         if (conprr(i,j) > 0.) then
            !----- Finding cloud base. If PBL top was computed, use it, otherwise guess. --!
            if (kpbl(i,j) > 0) then
                kcldbase = kpbl(i,j)
            elseif (any(tke(lpw:m1,i,j) > 0)) then
               kcldbase=(lpw-1)+maxloc(tke(lpw:m1,i,j),dim=1)
            else
               baseloop:do kcldbase=lpw,m1-1
                  if ((zt(kcldbase+1)-zm(lpw-1))*rtgt(i,j) > zcldbase) exit baseloop
               end do baseloop
            end if
            !----- Finding downdraft origin -----------------------------------------------!
            downloop:do kzdown=lpw,m1-1
               if ((zt(kzdown+1)-zm(lpw-1))*rtgt(i,j) > zdown) exit downloop
            end do downloop
            !------------------------------------------------------------------------------!
            !    If there was rain, put it back in the cloud, between the cloud base and   !
            ! the level in which downdrafts begin.                                         !
            !------------------------------------------------------------------------------!
            if (upmf(i,j) > 0.) then
              dnmf = epsilon*upmf(i,j)
            else
              dnmf = dnmfdef
            end if
            !----- Putting the rain as a constant between the cloud base and dndraft top --!
            do k=kcldbase,kzdown
               cupcond(k,i,j) = conprr(i,j) / dnmf
            end do
         end if

         !----- Now I assume that the positive tendency on moisture = cloud ---------------!
         do k=lpw,m1
            if (rtsrc(k,i,j) > 0.) cupcond(k,i,j) = cupcond(k,i,j) + rtsrc(k,i,j) * deltat
         end do

      end do
   end do
   return
end subroutine cloud_sketch
!==========================================================================================!
!==========================================================================================!
