!############################# Change Log #################################################!
! 5.0.0                                                                                    !
!                                                                                          !
!##########################################################################################!
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved                       !
!  Regional Atmospheric Modeling System - RAMS                                             !
!##########################################################################################!


module mem_cuparm

   use grid_dims

   type cuparm_vars

      !------------------------------------------------------------------------------------!
      !    Variables to be dimensioned by (nzp,nxp,nyp,nclouds). These are filled for each !
      ! cloud separatedly, because it is for output or cloud-dependent processes only. The !
      ! actual tendency will be added together with all other non-convective processes.    !
      !------------------------------------------------------------------------------------!
      real, pointer, dimension(:,:,:,:) :: &
            thsrc                          & ! Heating rate due to this cloud
           ,rtsrc                          & ! Moistening rate due to this cloud
           ,co2src                         & ! CO2 change rate due to this cloud
           ,cuprliq                        & ! Param. cumulus liquid water mixing ratio
           ,cuprice                        ! ! Param. cumulus ice mixing ratio

      !------------------------------------------------------------------------------------!
      !    Variables to be dimensioned by (nxp,nyp,nclouds). These are filled for each     !
      ! cloud separatedly, because it is for output or cloud-dependent processes only.     !
      !------------------------------------------------------------------------------------!
      real, pointer, dimension(:,:,:) :: &
            aadn                         & ! Downdraft cloud work (not always computed)
           ,aaup                         & ! Updraft cloud work (not always computed)
           ,areadn                       & ! Downdraft relative area.
           ,areaup                       & ! Updraft relative area.
           ,conprr                       & ! Convective precipitation rate
           ,dnmf                         & ! Reference downdraft mass flux
           ,edt                          & ! dnmf/upmf (cloud work related)
           ,upmf                         & ! Reference updraft mass flux
           ,xierr                        & ! Error flag
           ,zjmin                        & ! Downdraft originating level
           ,zk22                         & ! Updraft originating level
           ,zkdt                         & ! Top of detrainment level
           ,zkbcon                       & ! Level of free convection
           ,zktop                        & ! Cloud top
           ,xiact_c                      & ! For old Grell
           ,xiact_p                      ! ! For old Grell

      !------------------------------------------------------------------------------------!
      !     Variables to be dimensioned by (nzp,nxp,nyp), these are used for cumulus       !
      ! inversion                                                                          !
      !------------------------------------------------------------------------------------!
      real, pointer, dimension(:,:,:) :: &
            thsrcf                       &
           ,rtsrcf                       &
           ,thsrcp                       &
           ,rtsrcp                       !
      !------------------------------------------------------------------------------------!
      !    Variables to be dimensioned by (nxp,nyp). Aconpr and conprr should contain the  !
      ! precipitation rate due to all precipitating clouds                                 !
      !------------------------------------------------------------------------------------!
      real, pointer, dimension(:,:) :: &
            aconpr                     &
           ,conprrp                    &
           ,conprrf                    !

   end type cuparm_vars

   type (cuparm_vars), allocatable :: cuparm_g(:), cuparmm_g(:)


   !---------------------------------------------------------------------------------------!
   !  These variables will control all cumulus parameterizations (Grell, Souza, Kuo)       !
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !    Number of clouds I want to use (always reserve the first for deep, even when the   !
   ! deepest is off. This is to make sure we won't mess with the original Kuo parameter-   !
   ! ization or inverse Kuo                                                                !
   !---------------------------------------------------------------------------------------!
   integer                            :: nclouds
   !------ Flag to turn on the cumulus parameterization -----------------------------------!
   integer, dimension(maxgrds)        :: nnqparm
   !------ Flag to decide which parameterization to use for deepest cloud -----------------!
   integer, dimension(maxgrds)        :: ndeepest
   !------ Flag to decide which parameterization to use for shallowest cloud --------------!
   integer, dimension(maxgrds)        :: nshallowest
   !------ Index to tell which one is the first cloud for Grell to solve ------------------!
   integer, dimension(maxgrds)        :: grell_1st
   !------ Index to tell which one is the first cloud for Grell to solve ------------------!
   integer, dimension(maxgrds)        :: grell_last
   !------ How often the cumulus parameterization should be called? -----------------------!
   real                               :: confrq
   !------ Time to start running the parameterization for this cloud type -----------------!
   real(kind=8)                       :: cptime
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    Miscellaneous parameters                                                           !
   !---------------------------------------------------------------------------------------!
   !----- Parameters for inverse cumulus --------------------------------------------------!
   integer           , parameter             :: maxcufiles=100, maxcugrids=10
   integer                                   :: if_cuinv
   real(kind=8)                              :: tcu_beg,tcu_end
   real                                      :: cu_til,cu_tel,tnudcu
   real              , dimension(maxcugrids) :: wt_cu_grid
   character(len=128)                        :: cu_prefix
   character(len=128), dimension(maxcufiles) :: fnames_cu
   character(len=14) , dimension(maxcufiles) :: itotdate_cu
   real(kind=8)      , dimension(maxcufiles) :: cu_times
   integer                                   :: ncufiles,ncufl
   real(kind=8)                              :: cutime1,cutime2
   !----- Parameter for Kuo parameterization ----------------------------------------------!
   real                                      :: wcldbs ! Minimum vertical velocity trigger.
   !---------------------------------------------------------------------------------------!


   !=======================================================================================!
   !=======================================================================================!


   contains



   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine adjust the parameters given by the user to make sure the          !
   ! configuration is consistent. For example, it will set up the number of clouds and the !
   ! first and last cloud spectrum that will be computed by Grell parameterization, and    !
   ! turn off any clouds other than the deepest when running cumulus inversion.            !
   !---------------------------------------------------------------------------------------!
   subroutine define_cumulus_dimensions(ngrids)
      implicit none
      integer, intent(in) :: ngrids
      integer             :: ngr
      
      !------------------------------------------------------------------------------------!
      !    Cumulus inversion. By using it you cannot run any Grell or Souza shallow        !
      ! cumulus. Your only options are Kuo or nothing. Resetting all other variables to    !
      ! reflect this.                                                                      !
      !------------------------------------------------------------------------------------!
      if (if_cuinv == 1) then
         nclouds          = 1
         ndeepest         = 1
         nshallowest      = 0
         grell_1st        = maxclouds ! 1st greater than last will keep Greel from running.
         grell_last       = 0
      end if

      !------------------------------------------------------------------------------------!
      !  Now loop over the grids, to adjust or maybe turn of cumulus inversion.            !
      !------------------------------------------------------------------------------------!
      do ngr=1,ngrids
         !----- No convection, nothing should happen. -------------------------------------!
         if (nnqparm(ngr) == 0) then
            ndeepest(ngr)    = 0
            nshallowest(ngr) = 0
            grell_1st (ngr)  = maxclouds ! 1st greater than last = no Grell 
            grell_last(ngr)  = 0
         else
            if (ndeepest(ngr) == 2) then
               grell_1st(ngr) = 1
            else
               grell_1st(ngr) = 2
            end if
            if (nshallowest(ngr) == 2) then
               grell_last(ngr) = nclouds
            else
               grell_last(ngr) = nclouds-1
            end if
         end if
      end do
      return
   end subroutine define_cumulus_dimensions
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   subroutine alloc_cuparm(cuparm,n1,n2,n3,ng,co2_on)
      implicit none
      type (cuparm_vars) :: cuparm
      integer, intent(in) :: n1,n2,n3,ng
      logical, intent(in) :: co2_on

      !------------------------------------------------------------------------------------!
      !     Do not allocate anything if cumulus won't be used in this grid.                !
      !------------------------------------------------------------------------------------!
      if(nnqparm(ng) == 0) return


      !------------------------------------------------------------------------------------!
      !     These variables should always be allocated if any cumulus parameterization was !
      ! on.                                                                                !
      !------------------------------------------------------------------------------------!
      allocate (cuparm%thsrc   (n1,n2,n3,nclouds))
      allocate (cuparm%rtsrc   (n1,n2,n3,nclouds))
      allocate (cuparm%cuprliq (n1,n2,n3,nclouds))
      allocate (cuparm%cuprice (n1,n2,n3,nclouds))
      allocate (cuparm%aconpr     (n2,n3)        )
      allocate (cuparm%upmf       (n2,n3,nclouds))
      allocate (cuparm%conprr     (n2,n3,nclouds))
      allocate (cuparm%xiact_c    (n2,n3,nclouds))
      allocate (cuparm%xiact_p    (n2,n3,nclouds))
      allocate (cuparm%areadn     (n2,n3,nclouds))
      allocate (cuparm%areaup     (n2,n3,nclouds))

      !----- If CO2 is on, allocate the CO2 source. ---------------------------------------!
      if (co2_on) allocate (cuparm%co2src  (n1,n2,n3,nclouds))

      !----- If cumulus inversion is on, include extra variables. -------------------------!
      if (if_cuinv == 1) then
         allocate (cuparm%thsrcp  (n1,n2,n3)        )
         allocate (cuparm%rtsrcp  (n1,n2,n3)        )
         allocate (cuparm%thsrcf  (n1,n2,n3)        )
         allocate (cuparm%rtsrcf  (n1,n2,n3)        )
         allocate (cuparm%conprrp    (n2,n3)        )
         allocate (cuparm%conprrf    (n2,n3)        )
      endif

      !----- If Grell will be used in this grid, allocate the other variables -------------!
      if (grell_1st(ng) <= grell_last(ng) .or. ndeepest(ng) == 3 .or.                      &
          nshallowest(ng) == 3) then
         allocate (cuparm%aadn       (n2,n3,nclouds))
         allocate (cuparm%aaup       (n2,n3,nclouds))
         allocate (cuparm%dnmf       (n2,n3,nclouds))
         allocate (cuparm%edt        (n2,n3,nclouds))
         allocate (cuparm%upmf       (n2,n3,nclouds))
         allocate (cuparm%xierr      (n2,n3,nclouds))
         allocate (cuparm%zjmin      (n2,n3,nclouds))
         allocate (cuparm%zk22       (n2,n3,nclouds))
         allocate (cuparm%zkdt       (n2,n3,nclouds))
         allocate (cuparm%zkbcon     (n2,n3,nclouds))
         allocate (cuparm%zktop      (n2,n3,nclouds))
      end if

      return
   end subroutine alloc_cuparm
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   subroutine nullify_cuparm(cuparm)

      implicit none
      type (cuparm_vars) :: cuparm

      if(associated(cuparm%thsrc      ))  nullify (cuparm%thsrc      )
      if(associated(cuparm%rtsrc      ))  nullify (cuparm%rtsrc      )
      if(associated(cuparm%co2src     ))  nullify (cuparm%co2src     )
      if(associated(cuparm%areadn     ))  nullify (cuparm%areadn     )
      if(associated(cuparm%areaup     ))  nullify (cuparm%areaup     )
      if(associated(cuparm%cuprliq    ))  nullify (cuparm%cuprliq    )
      if(associated(cuparm%cuprice    ))  nullify (cuparm%cuprice    )
      if(associated(cuparm%aconpr     ))  nullify (cuparm%aconpr     )
      if(associated(cuparm%conprr     ))  nullify (cuparm%conprr     )
      if(associated(cuparm%thsrcp     ))  nullify (cuparm%thsrcp     )
      if(associated(cuparm%rtsrcp     ))  nullify (cuparm%rtsrcp     )
      if(associated(cuparm%thsrcf     ))  nullify (cuparm%thsrcf     )
      if(associated(cuparm%rtsrcf     ))  nullify (cuparm%rtsrcf     )
      if(associated(cuparm%conprrp    ))  nullify (cuparm%conprrp    )
      if(associated(cuparm%conprrf    ))  nullify (cuparm%conprrf    )
      if(associated(cuparm%aadn       ))  nullify (cuparm%aadn       )
      if(associated(cuparm%aaup       ))  nullify (cuparm%aaup       )
      if(associated(cuparm%dnmf       ))  nullify (cuparm%dnmf       )
      if(associated(cuparm%edt        ))  nullify (cuparm%edt        )
      if(associated(cuparm%upmf       ))  nullify (cuparm%upmf       )
      if(associated(cuparm%xierr      ))  nullify (cuparm%xierr      )
      if(associated(cuparm%zjmin      ))  nullify (cuparm%zjmin      )
      if(associated(cuparm%zk22       ))  nullify (cuparm%zk22       )
      if(associated(cuparm%zkdt       ))  nullify (cuparm%zkdt       )
      if(associated(cuparm%zkbcon     ))  nullify (cuparm%zkbcon     )
      if(associated(cuparm%zktop      ))  nullify (cuparm%zktop      )
      if(associated(cuparm%xiact_c    ))  nullify (cuparm%xiact_c    )
      if(associated(cuparm%xiact_p    ))  nullify (cuparm%xiact_p    )

      return
   end subroutine nullify_cuparm
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   subroutine dealloc_cuparm(cuparm)

      implicit none

      type (cuparm_vars) :: cuparm
      if(associated(cuparm%thsrc      ))  deallocate (cuparm%thsrc      )
      if(associated(cuparm%rtsrc      ))  deallocate (cuparm%rtsrc      )
      if(associated(cuparm%co2src     ))  deallocate (cuparm%co2src     )
      if(associated(cuparm%areadn     ))  deallocate (cuparm%areadn     )
      if(associated(cuparm%areaup     ))  deallocate (cuparm%areaup     )
      if(associated(cuparm%cuprliq    ))  deallocate (cuparm%cuprliq    )
      if(associated(cuparm%cuprice    ))  deallocate (cuparm%cuprice    )
      if(associated(cuparm%aconpr     ))  deallocate (cuparm%aconpr     )
      if(associated(cuparm%conprr     ))  deallocate (cuparm%conprr     )
      if(associated(cuparm%thsrcp     ))  deallocate (cuparm%thsrcp     )
      if(associated(cuparm%rtsrcp     ))  deallocate (cuparm%rtsrcp     )
      if(associated(cuparm%thsrcf     ))  deallocate (cuparm%thsrcf     )
      if(associated(cuparm%rtsrcf     ))  deallocate (cuparm%rtsrcf     )
      if(associated(cuparm%conprrp    ))  deallocate (cuparm%conprrp    )
      if(associated(cuparm%conprrf    ))  deallocate (cuparm%conprrf    )
      if(associated(cuparm%aadn       ))  deallocate (cuparm%aadn       )
      if(associated(cuparm%aaup       ))  deallocate (cuparm%aaup       )
      if(associated(cuparm%dnmf       ))  deallocate (cuparm%dnmf       )
      if(associated(cuparm%edt        ))  deallocate (cuparm%edt        )
      if(associated(cuparm%upmf       ))  deallocate (cuparm%upmf       )
      if(associated(cuparm%xierr      ))  deallocate (cuparm%xierr      )
      if(associated(cuparm%zjmin      ))  deallocate (cuparm%zjmin      )
      if(associated(cuparm%zk22       ))  deallocate (cuparm%zk22       )
      if(associated(cuparm%zkdt       ))  deallocate (cuparm%zkdt       )
      if(associated(cuparm%zkbcon     ))  deallocate (cuparm%zkbcon     )
      if(associated(cuparm%zktop      ))  deallocate (cuparm%zktop      )
      if(associated(cuparm%xiact_c    ))  deallocate (cuparm%xiact_c    )
      if(associated(cuparm%xiact_p    ))  deallocate (cuparm%xiact_p    )

      return
   end subroutine dealloc_cuparm
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   subroutine initialize_cuparm(cuparm)

      implicit none

      type (cuparm_vars) :: cuparm
      if(associated(cuparm%thsrc      ))  cuparm%thsrc      = 0.
      if(associated(cuparm%rtsrc      ))  cuparm%rtsrc      = 0.
      if(associated(cuparm%co2src     ))  cuparm%co2src     = 0.
      if(associated(cuparm%upmf       ))  cuparm%upmf       = 0.
      if(associated(cuparm%areadn     ))  cuparm%areadn     = 0.
      if(associated(cuparm%areaup     ))  cuparm%areaup     = 0.
      if(associated(cuparm%cuprliq    ))  cuparm%cuprliq    = 0.
      if(associated(cuparm%cuprice    ))  cuparm%cuprice    = 0.
      if(associated(cuparm%aconpr     ))  cuparm%aconpr     = 0.
      if(associated(cuparm%conprr     ))  cuparm%conprr     = 0.
      if(associated(cuparm%thsrcp     ))  cuparm%thsrcp     = 0.
      if(associated(cuparm%rtsrcp     ))  cuparm%rtsrcp     = 0.
      if(associated(cuparm%thsrcf     ))  cuparm%thsrcf     = 0.
      if(associated(cuparm%rtsrcf     ))  cuparm%rtsrcf     = 0.
      if(associated(cuparm%conprrp    ))  cuparm%conprrp    = 0.
      if(associated(cuparm%conprrf    ))  cuparm%conprrf    = 0.
      if(associated(cuparm%aadn       ))  cuparm%aadn       = 0.
      if(associated(cuparm%aaup       ))  cuparm%aaup       = 0.
      if(associated(cuparm%dnmf       ))  cuparm%dnmf       = 0.
      if(associated(cuparm%edt        ))  cuparm%edt        = 0.
      if(associated(cuparm%upmf       ))  cuparm%upmf       = 0.
      if(associated(cuparm%xierr      ))  cuparm%xierr      = 1. !---- No convection
      if(associated(cuparm%zjmin      ))  cuparm%zjmin      = 0.
      if(associated(cuparm%zk22       ))  cuparm%zk22       = 0.
      if(associated(cuparm%zkdt       ))  cuparm%zkdt       = 0.
      if(associated(cuparm%zkbcon     ))  cuparm%zkbcon     = 0.
      if(associated(cuparm%zktop      ))  cuparm%zktop      = 0.
      if(associated(cuparm%xiact_c    ))  cuparm%xiact_c    = 0.
      if(associated(cuparm%xiact_p    ))  cuparm%xiact_p    = 0.

      return
   end subroutine initialize_cuparm
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   subroutine filltab_cuparm(cuparm,cuparmm,imean,n1,n2,n3,ng)
      use var_tables

      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type (cuparm_vars) , intent(in) :: cuparm
      type (cuparm_vars) , intent(in) :: cuparmm
      integer            , intent(in) :: imean
      integer            , intent(in) :: n1
      integer            , intent(in) :: n2
      integer            , intent(in) :: n3
      integer            , intent(in) :: ng
      !----- Local variables. -------------------------------------------------------------!
      integer :: npts
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !   Four-dimensional, space + cloud spectral size                                    ! 
      !------------------------------------------------------------------------------------!
      npts=n1*n2*n3*nclouds
      if (associated(cuparm%thsrc))                                                        &
         call vtables2(cuparm%thsrc(1,1,1,1),cuparmm%thsrc(1,1,1,1),ng,npts,imean          &
                       ,'THSRC :8:hist:anal:mpti:mpt3')

      if (associated(cuparm%rtsrc))                                                        &
         call vtables2(cuparm%rtsrc(1,1,1,1),cuparmm%rtsrc(1,1,1,1),ng,npts,imean          &
                      ,'RTSRC :8:hist:anal:mpti:mpt3')

      if (associated(cuparm%co2src))                                                       &
           call vtables2(cuparm%co2src(1,1,1,1),cuparmm%co2src(1,1,1,1),ng,npts,imean      &
                        ,'CO2SRC :8:hist:anal:mpti:mpt3')

      if (associated(cuparm%cuprliq))                                                      &
           call vtables2(cuparm%cuprliq(1,1,1,1),cuparmm%cuprliq(1,1,1,1),ng,npts,imean    &
                        ,'CUPRLIQ :8:hist:anal:mpti:mpt3')                              

      if (associated(cuparm%cuprice))                                                      &
           call vtables2(cuparm%cuprice(1,1,1,1),cuparmm%cuprice(1,1,1,1),ng,npts,imean    &
                        ,'CUPRICE :8:hist:anal:mpti:mpt3')
      !--------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !   Three-dimensional, space only                                                    ! 
      !------------------------------------------------------------------------------------!
      npts=n1*n2*n3
      if (associated(cuparm%thsrcp))                                                       &
         call vtables2(cuparm%thsrcp(1,1,1),cuparmm%thsrcp(1,1,1),ng,npts,imean            &
                      ,'THSRCP :3:mpti:mpt3')
      if (associated(cuparm%rtsrcp))                                                       &
         call vtables2(cuparm%rtsrcp(1,1,1),cuparmm%rtsrcp(1,1,1),ng,npts,imean            &
                      ,'RTSRCP :3:mpti:mpt3')
      if (associated(cuparm%thsrcf))                                                       &
         call vtables2(cuparm%thsrcf(1,1,1),cuparmm%thsrcf(1,1,1),ng,npts,imean            &
                      ,'THSRCF :3:mpti:mpt3')
      if (associated(cuparm%rtsrcf))                                                       &
         call vtables2(cuparm%rtsrcf(1,1,1),cuparmm%rtsrcf(1,1,1),ng,npts,imean            &
                      ,'RTSRCF :3:mpti:mpt3')
      !--------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !   Three-dimensional, horizontal + cloud spectrum                                   ! 
      !------------------------------------------------------------------------------------!
      npts=n2*n3*nclouds

      if (associated(cuparm%aadn))                                                         &
         call vtables2(cuparm%aadn(1,1,1),cuparmm%aadn(1,1,1),ng,npts,imean                &
                      ,'AADN :9:hist:anal:mpti:mpt3')

      if (associated(cuparm%aaup))                                                         &
         call vtables2(cuparm%aaup(1,1,1),cuparmm%aaup(1,1,1),ng, npts, imean              &
                      ,'AAUP :9:hist:anal:mpti:mpt3')

      if (associated(cuparm%areadn))                                                       &
         call vtables2(cuparm%areadn(1,1,1),cuparmm%areadn(1,1,1),ng,npts,imean            &
                      ,'AREADN :9:hist:anal:mpti:mpt3')

      if (associated(cuparm%areaup))                                                       &
         call vtables2(cuparm%areaup(1,1,1),cuparmm%areaup(1,1,1),ng,npts,imean            &
                      ,'AREAUP :9:hist:anal:mpti:mpt3')

      if (associated(cuparm%conprr))                                                       &
         call vtables2(cuparm%conprr(1,1,1),cuparmm%conprr(1,1,1),ng,npts,imean            &
                      ,'CONPRR :9:hist:anal:mpti:mpt3')

      if (associated(cuparm%dnmf))                                                         &
         call vtables2(cuparm%dnmf(1,1,1),cuparmm%dnmf(1,1,1),ng,npts,imean                &
                      ,'DNMF :9:hist:anal:mpti:mpt3')

      if (associated(cuparm%edt))                                                          &
         call vtables2(cuparm%edt(1,1,1),cuparmm%edt(1,1,1),ng,npts,imean                  &
                      ,'EDT :9:hist:anal:mpti:mpt3')

      if (associated(cuparm%upmf))                                                         &
         call vtables2(cuparm%upmf(1,1,1),cuparmm%upmf(1,1,1),ng,npts,imean                &
                      ,'UPMF :9:hist:anal:mpti:mpt3')

      if (associated(cuparm%xierr))                                                        &
         call vtables2(cuparm%xierr(1,1,1),cuparmm%xierr(1,1,1),ng,npts,imean              &
                      ,'XIERR :9:hist:anal:mpti:mpt3')

      if (associated(cuparm%zjmin))                                                        &
         call vtables2(cuparm%zjmin(1,1,1),cuparmm%zjmin(1,1,1),ng,npts,imean              &
                      ,'ZJMIN :9:hist:anal:mpti:mpt3')

      if (associated(cuparm%zk22))                                                         &
         call vtables2(cuparm%zk22(1,1,1),cuparmm%zk22(1,1,1),ng,npts,imean                &
                      ,'ZK22 :9:hist:anal:mpti:mpt3')

      if (associated(cuparm%zkdt))                                                         &
         call vtables2(cuparm%zkdt(1,1,1),cuparmm%zkdt(1,1,1),ng,npts,imean                &
                      ,'ZKDT :9:hist:anal:mpti:mpt3')

      if (associated(cuparm%zkbcon))                                                       &
         call vtables2(cuparm%zkbcon(1,1,1),cuparmm%zkbcon(1,1,1),ng,npts,imean            &
                      ,'ZKBCON :9:hist:anal:mpti:mpt3')

      if (associated(cuparm%zktop))                                                        &
         call vtables2(cuparm%zktop(1,1,1),cuparmm%zktop(1,1,1),ng,npts,imean              &
                      ,'ZKTOP :9:hist:anal:mpti:mpt3')

      if (associated(cuparm%xiact_c))                                                      &
         call vtables2(cuparm%xiact_c(1,1,1),cuparmm%xiact_c(1,1,1),ng,npts,imean          &
                      ,'XIACT_C :9:hist:anal:mpti:mpt3')

      if (associated(cuparm%xiact_p))                                                      &
         call vtables2(cuparm%xiact_p(1,1,1),cuparmm%xiact_p(1,1,1),ng,npts,imean          &
                      ,'XIACT_P :9:hist:anal:mpti:mpt3')
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !   Two-dimensional variables, in space only                                         !
      !------------------------------------------------------------------------------------!
      npts=n2*n3
      if (associated(cuparm%aconpr))                                                       &
         call vtables2(cuparm%aconpr(1,1),cuparmm%aconpr(1,1),ng,npts,imean                &
                      ,'ACONPR :2:hist:anal:mpti:mpt3')
      if (associated(cuparm%conprrp))                                                      &
         call vtables2(cuparm%conprrp(1,1),cuparmm%conprrp(1,1),ng,npts,imean              &
                      ,'CONPRRP :2:mpti')
      if (associated(cuparm%conprrf))                                                      &
         call vtables2(cuparm%conprrf(1,1),cuparmm%conprrf(1,1),ng,npts,imean              &
                      ,'CONPRRF :2:mpti')
      !------------------------------------------------------------------------------------!

      return
   end subroutine filltab_cuparm
   !=======================================================================================!
   !=======================================================================================!
end module mem_cuparm
!==========================================================================================!
!==========================================================================================!
