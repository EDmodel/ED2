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
subroutine rconv_driver(banneron)
   use mem_cuparm, only : &
           confrq       &  ! How often I should call this cloud
          ,cptime       &  ! Time to start calling this cloud
          ,nclouds      &  ! Number of clouds to use here 
          ,nnqparm      &  ! Flag for cumulus parameterization
          ,ndeepest     &  ! Flag for deepest cloud
          ,nshallowest  &  ! Flag for shallowest cloud
          ,grell_1st    &  ! First cloud to use Grell's parameterization
          ,grell_last   &  ! Last cloud to use Grell's parameterization
          ,cuparm_g     !
          
   use mem_grid, only:  &
           initial      &  ! Flag for "initial" run
          ,ngrid        &  ! Current grid
          ,ngrids       &  ! Total # of grids
          ,time         &  ! Current time
          ,dtlt         &  ! Current grid delta-t
          ,grid_g       !  ! Grid structure

   use node_mod, only:  &
           mzp          &  ! Current grid # of vertical levels for this node
          ,myp          &  ! Current grid # of meridional points for this node
          ,mxp          &  ! Current grid # of zonal points for this node
          ,ia           &  ! Westernmost boundary
          ,iz           &  ! Easternmost boundary
          ,ja           &  ! Southernmost boundary
          ,jz           !  ! Northermost boundary
   
   use mem_turb, only:  &
           turb_g       &  ! Turbulence structure
          ,idiffk       !  ! Turbulence closure flag
   
   use mem_scratch, only : & 
           scratch      !  ! Scratch structure
   use mem_tend, only:  &
           tend         !  ! Tendency structure

   implicit none
   logical, intent(in) :: banneron          ! Flag to print the banner
   integer             :: icld              ! Cloud counter
   integer             :: ifm               ! Grid counter
   logical, external   :: cumulus_time      ! Function to check whether it's cumulus time.
   logical, save       :: first_time=.true. ! Flag to check whether this is the 1st call

   !----- Nothing to be done in here for this grid ----------------------------------------!
   if (nnqparm(ngrid) == 0) return

   !---------------------------------------------------------------------------------------!
   !    Checking whether this is the first time I access the subroutine. If so, perform a  !
   ! handful of initialisations. This should be unecessary because the variables are       !
   ! initialised at the rams_mem_alloc time, but apparently it is reset again.             !
   !---------------------------------------------------------------------------------------!
   if (first_time) then
      do ifm=1,ngrids
         cuparm_g(ifm)%rtsrc  = 0.
         cuparm_g(ifm)%thsrc  = 0.
         cuparm_g(ifm)%conprr = 0.
         if (associated(cuparm_g(ifm)%dnmf))   cuparm_g(ifm)%dnmf  = 0.
         if (associated(cuparm_g(ifm)%xierr))  cuparm_g(ifm)%xierr = 1.
      end do
      first_time=.false.
   end if


   !---------------------------------------------------------------------------------------!
   !   Checking whether I should use Souza's or old Grell's scheme                         !
   !---------------------------------------------------------------------------------------!
   if ((nshallowest(ngrid) == 1 .or. nshallowest(ngrid) == 3) .and.                        &
        cumulus_time(initial,time,cptime(nclouds),confrq(nclouds),dtlt)) then
      select case(nshallowest(ngrid))
      case (1) !
         call souza_cupar_driver()
      case (3)
         call old_grell_cupar_driver(nclouds)
      end select
      
      !----- If I have TKE available, I may use it to define PBL top ----------------------!
      select case (idiffk(ngrid))
      case (1,4,5,6,7)
         call atob(mzp*mxp*myp,turb_g(ngrid)%tkep,scratch%vt3de)
      case default
         call azero(mzp*mxp*myp,scratch%vt3de)
      end select
      
      !----- This is just temporary, it is not working properly ---------------------------!
      call cloud_sketch(mzp,mxp,myp,ia,iz,ja,jz,dtlt                                       &
         , grid_g(ngrid)%flpw                   , grid_g(ngrid)%rtgt                       &
         , turb_g(ngrid)%kpbl                   , scratch%vt3de                            &
         , cuparm_g(ngrid)%upmf    (:,:,nclouds), cuparm_g(ngrid)%rtsrc    (:,:,:,nclouds) &
         , cuparm_g(ngrid)%conprr  (:,:,nclouds), cuparm_g(ngrid)%cuprliq  (:,:,:,nclouds) )
   end if


   !---------------------------------------------------------------------------------------!
   !   Looping accross the cloud sizes, from last to first                                 !
   !---------------------------------------------------------------------------------------!
   do icld = grell_last(ngrid),grell_1st(ngrid),-1
      !----- Checking whether this is time to call Grell's parameterization ---------------!
      if (cumulus_time(initial,time,cptime(icld),confrq(icld),dtlt)) then
         call grell_cupar_driver(banneron,icld)
      end if
   end do

   !---------------------------------------------------------------------------------------!
   !   Checking whether I was supposed to run Kuo for this grid.                           !
   !---------------------------------------------------------------------------------------!
   if ((ndeepest(ngrid) == 1 .or. ndeepest(ngrid) == 3) .and.                              &
       cumulus_time(initial,time,cptime(1),confrq(1),dtlt)) then
       select case (ndeepest(ngrid))
       case (1)
          call kuo_cupar_driver()
       case (3)
          call old_grell_cupar_driver(1)
       end select


       !----- If I have TKE available, I may use it to define PBL top ---------------------!
       select case (idiffk(ngrid))
       case (1,4,5,6,7)
          call atob(mzp*mxp*myp,turb_g(ngrid)%tkep,scratch%vt3de)
       case default
          call azero(mzp*mxp*myp,scratch%vt3de)
       end select
       
       !----- This is temporary, it is not working properly -------------------------------!
       call cloud_sketch(mzp,mxp,myp,ia,iz,ja,jz,dtlt                                      &
         , grid_g(ngrid)%flpw                     , grid_g(ngrid)%rtgt                     &
         , turb_g(ngrid)%kpbl                     , scratch%vt3de                          &
         , cuparm_g(ngrid)%upmf           (:,:,1) , cuparm_g(ngrid)%rtsrc        (:,:,:,1) &
         , cuparm_g(ngrid)%conprr         (:,:,1) , cuparm_g(ngrid)%cuprliq      (:,:,:,1) )
   end if

   !---------------------------------------------------------------------------------------!
   !   Updating the heating and moistening rates, and the convective precipitation.        !
   !---------------------------------------------------------------------------------------!
   do icld=1,nclouds
      !----- Accumulating the tendencies --------------------------------------------------!
      call accum(mxp*myp*mzp,tend%tht,cuparm_g(ngrid)%thsrc(:,:,:,icld)    )
      call accum(mxp*myp*mzp,tend%rtt,cuparm_g(ngrid)%rtsrc(:,:,:,icld)    )
      !----- Updating total precipitation -------------------------------------------------!
      call update(mxp*myp,cuparm_g(ngrid)%aconpr,cuparm_g(ngrid)%conprr(:,:,icld),dtlt)
   end do


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

   integer                  , intent(in)  :: m1,m2,m3    ! Matrix dimensions
   integer                  , intent(in)  :: ia,iz,ja,jz ! Node boundaries
   real                     , intent(in)  :: deltat      ! Time step
   real, dimension   (m2,m3), intent(in)  :: flpw        ! Lowest point above ground
   real, dimension   (m2,m3), intent(in)  :: rtgt        ! Correction
   integer, dimension(m2,m3), intent(in)  :: kpbl        ! PBL top
   real, dimension(m1,m2,m3), intent(in)  :: tke         ! TKE
   real, dimension(m1,m2,m3), intent(in)  :: rtsrc       ! Moistening rate
   real, dimension   (m2,m3), intent(in)  :: conprr      ! Precipitation rate
   real, dimension   (m2,m3), intent(in)  :: upmf        ! Updraft reference 
   real, dimension(m1,m2,m3), intent(out) :: cupcond     ! My Frankeinstein cloud
   
   integer :: i,j,k,lpw,kcldbase,kzdown
   real :: dnmf
   
   real, parameter :: zcldbase=1600. ! Fixed height for cloud base
   real, parameter :: zdown=3000.    ! Fixed height for downdraft origin
   real, parameter :: epsilon= 0.3   ! Just an arbitrary epsilon, dnmf/upmf
   real, parameter :: dnmfdef= 0.1   ! Just an arbitrary defaul dnmf

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
