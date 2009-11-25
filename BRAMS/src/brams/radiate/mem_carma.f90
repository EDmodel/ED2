!==================================== Change Log ==========================================!
! 5.0.0                                                                                    !
!                                                                                          !
!==========================================================================================!
!==========================================================================================!
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved                       ! 
!  Regional Atmospheric Modeling System - RAMS                                             ! 
!==========================================================================================!
module mem_carma

   use grid_dims

   implicit none

   type carma_v
      real, pointer, dimension(:,:,:) :: aot
   end type carma_v

   type (carma_v), dimension(:), allocatable :: carma, carma_m

   !----- Used in radcomp_carma (radcomp) - radriv. ---------------------------------------!
   real :: tdec
   real :: sdec
   real :: cdec
   real :: declin
   real :: rvr(nzpmax)
   real :: rtr(nzpmax)
   real :: dn0r(nzpmax)
   real :: pird(nzpmax)
   real :: prd(nzpmax)
   real :: temprd(nzpmax+1)  
   real :: fthrl(nzpmax)
   real :: dzmr(nzpmax)
   real :: dztr(nzpmax)
   real :: fthrs(nzpmax) 

   !---------------------------------------------------------------------------------------!
   !    These variables are scratch variables/arrays used by CARMA.                        !
   !---------------------------------------------------------------------------------------!
   real                                :: p_surf 
   real                                :: p_top  
   real                                :: t_surf 
   real                                :: tabove_aerad 
   real                                :: solnet
   real                                :: xirdown
   real                                :: RAIN_aerad
   real                                :: xland_aerad

   logical                             :: isl_aerad 
   logical                             :: ir_aerad 
   integer                             :: lla
   integer                             :: lls
     
   real, allocatable, dimension(:)     :: p
   real, allocatable, dimension(:)     :: t
   real, allocatable, dimension(:)     :: t_aerad
   real, allocatable, dimension(:)     :: p_aerad
   real, allocatable, dimension(:)     :: qv_aerad

   real, allocatable, dimension(:)     :: LWL_aerad
   real, allocatable, dimension(:)     :: IWL_aerad
   real, allocatable, dimension(:)     :: LWP_aerad
   real, allocatable, dimension(:)     :: IWP_aerad


   real, allocatable, dimension(:)     :: press   
   real, allocatable, dimension(:)     :: dpg     
   real, allocatable, dimension(:)     :: tt      
   real, allocatable, dimension(:)     :: rhoa  
   real, allocatable, dimension(:)     :: rsfx   
   real, allocatable, dimension(:)     :: heats_aerad 
   real, allocatable, dimension(:)     :: heati_aerad 
   real, allocatable, dimension(:)     :: rdh2o   
   real, allocatable, dimension(:)     :: ptempg  
   real, allocatable, dimension(:)     :: ptempt  
   real, allocatable, dimension(:)     :: u1i

   real, allocatable, dimension(:,:)   :: gc     
   real, allocatable, dimension(:,:)   :: pah2o  
   real, allocatable, dimension(:,:)   :: paco2  
   real, allocatable, dimension(:,:)   :: pao2   
   real, allocatable, dimension(:,:)   :: pao3   
   real, allocatable, dimension(:,:)   :: taugas 
   real, allocatable, dimension(:,:)   :: paray  
   real, allocatable, dimension(:,:)   :: tauaer 
   real, allocatable, dimension(:,:)   :: taul   
   real, allocatable, dimension(:,:)   :: taucld 
   real, allocatable, dimension(:,:)   :: wcld   
   real, allocatable, dimension(:,:)   :: gcld   
   real, allocatable, dimension(:,:)   :: w0   
   real, allocatable, dimension(:,:)   :: g0

   real, allocatable, dimension(:,:)   :: taucldlw 
   real, allocatable, dimension(:,:)   :: wolc   
   real, allocatable, dimension(:,:)   :: gl 
   real, allocatable, dimension(:,:)   :: taucldice 
   real, allocatable, dimension(:,:)   :: woice   
   real, allocatable, dimension(:,:)   :: gice 

   real, allocatable, dimension(:,:)   :: opd   
   real, allocatable, dimension(:,:)   :: uopd   
   real, allocatable, dimension(:,:)   :: ptemp  
   real, allocatable, dimension(:,:)   :: slope  
   real, allocatable, dimension(:,:)   :: b1
   real, allocatable, dimension(:,:)   :: b2
   real, allocatable, dimension(:,:)   :: b3
   real, allocatable, dimension(:,:)   :: ak
   real, allocatable, dimension(:,:)   :: gami
   real, allocatable, dimension(:,:)   :: ee1
   real, allocatable, dimension(:,:)   :: el1
   real, allocatable, dimension(:,:)   :: em1
   real, allocatable, dimension(:,:)   :: el2
   real, allocatable, dimension(:,:)   :: em2
   real, allocatable, dimension(:,:)   :: af
   real, allocatable, dimension(:,:)   :: bf
   real, allocatable, dimension(:,:)   :: ef
   real, allocatable, dimension(:,:)   :: cp
   real, allocatable, dimension(:,:)   :: cpb
   real, allocatable, dimension(:,:)   :: cmb
   real, allocatable, dimension(:,:)   :: ck1
   real, allocatable, dimension(:,:)   :: ck2
   real, allocatable, dimension(:,:)   :: fnet
   real, allocatable, dimension(:,:)   :: tmi
   real, allocatable, dimension(:,:)   :: tmid
   real, allocatable, dimension(:,:)   :: tmiu
   real, allocatable, dimension(:,:)   :: direc
   real, allocatable, dimension(:,:)   :: directu
     
   real, allocatable, dimension(:,:,:) :: pc     
   real, allocatable, dimension(:,:,:) :: pc_aerad 
   real, allocatable, dimension(:,:,:) :: caer    
   real, allocatable, dimension(:,:,:) :: y3    

  
   real,parameter :: emisir_aerad = 1.0
   real,parameter :: h2ocol_aerad = 0.01

   !=======================================================================================!
   !=======================================================================================!


   contains



   !=======================================================================================!
   !=======================================================================================!
   subroutine alloc_carma(car,ng,nmxp,nmyp,nw)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type (carma_v), dimension(:), intent(inout) :: car
      integer                     , intent(in)    :: ng
      integer                     , intent(in)    :: nw
      integer                     , intent(in)    :: nmxp
      integer                     , intent(in)    :: nmyp
      !------------------------------------------------------------------------------------!

      allocate (car(ng)%aot(nmxp,nmyp,nw))

      return
   end subroutine alloc_carma
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   subroutine nullify_carma(car,ng)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type (carma_v), dimension(:), intent(inout) :: car
      integer                     , intent(in)    :: ng
      !------------------------------------------------------------------------------------!

      if (associated(car(ng)%aot ))  nullify (car(ng)%aot )

  
      return
   end subroutine nullify_carma
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   subroutine dealloc_carma(car,ng) 
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type (carma_v), dimension(:), intent(inout) :: car
      integer                     , intent(in)    :: ng
      !------------------------------------------------------------------------------------!
  
      if (associated(car(ng)%aot ))  deallocate (car(ng)%aot )

      return
   end subroutine dealloc_carma
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   subroutine zero_carma(car,ng)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type (carma_v), dimension(:), intent(inout) :: car
      integer                     , intent(in)    :: ng
      !------------------------------------------------------------------------------------!
      if (associated(car(ng)%aot )) car(ng)%aot(:,:,:)=0. 

      return
   end subroutine zero_carma
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   subroutine filltab_carma(cv,cvm,ng,imean,n1,n2,n3)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type (carma_v), intent(inout) :: cv
      type (carma_v), intent(inout) :: cvm
      integer       , intent(in)    :: ng
      integer       , intent(in)    :: n1
      integer       , intent(in)    :: n2
      integer       , intent(in)    :: n3
      integer       , intent(in)    :: imean
      !----- Local variables. -------------------------------------------------------------!
      integer                       :: npts
      character(len=7)              :: sname
      !------------------------------------------------------------------------------------!

      if (associated(cv%aot)) then
         npts=n1*n2*n3
         write(sname,fmt='(a4)') 'AOT'
         call vtables2(cv%aot(1,1,1),cvm%aot(1,1,1),ng,npts,imean                          &
                      ,sname//' :7:hist:anal:mpti:mpt3')
      end if
      return
   end subroutine filltab_carma
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   subroutine init_carma(m1)
      use mem_globrad , only : ntotal   & ! intent(in)
                             , nlayer   & ! intent(in)
                             , ndbl     & ! intent(in)
                             , ngauss   & ! intent(in)
                             , nrad     ! ! intent(in)
      use mem_aerad   , only : nz_rad   & ! intent(in)
                             , ngas     & ! intent(in)
                             , nelem    & ! intent(in)
                             , ngroup   & ! intent(in)
                             , nbin     ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      integer, intent(in) :: m1
      !------------------------------------------------------------------------------------!

      !----- 1D variables. ----------------------------------------------------------------!
      if(.not. allocated(p          ))    allocate(p           (m1))
      if(.not. allocated(t          ))    allocate(t           (m1))
      if(.not. allocated(t_aerad    ))    allocate(t_aerad     (m1))
      if(.not. allocated(p_aerad    ))    allocate(p_aerad     (m1))
      if(.not. allocated(qv_aerad   ))    allocate(qv_aerad    (m1))
      if(.not. allocated(LWL_aerad  ))    allocate(LWL_aerad   (m1))
      if(.not. allocated(IWL_aerad  ))    allocate(IWL_aerad   (m1))
      if(.not. allocated(LWP_aerad  ))    allocate(LWP_aerad   (m1))
      if(.not. allocated(IWP_aerad  ))    allocate(IWP_aerad   (m1))
      if(.not. allocated(rhoa       ))    allocate(rhoa        (m1))

      if(.not. allocated(press      )) allocate(press      (nlayer))
      if(.not. allocated(dpg        )) allocate(dpg        (nlayer))
      if(.not. allocated(tt         )) allocate(tt         (nlayer))
      if(.not. allocated(rdh2o      )) allocate(rdh2o      (nlayer))
      if(.not. allocated(rsfx       )) allocate(rsfx       (ntotal))
      if(.not. allocated(heats_aerad)) allocate(heats_aerad(nz_rad))
      if(.not. allocated(heati_aerad)) allocate(heati_aerad(nz_rad))
      if(.not. allocated(ptempg     )) allocate(ptempg     (ntotal))
      if(.not. allocated(ptempt     )) allocate(ptempt     (ntotal))
      if(.not. allocated(u1i        )) allocate(u1i        (ntotal))

      !----- 2D variables. ----------------------------------------------------------------!
      if(.not. allocated(gc         )) allocate(gc         (    m1,ngas  ))
      if(.not. allocated(pah2o      )) allocate(pah2o      (ntotal,nlayer))
      if(.not. allocated(paco2      )) allocate(paco2      (ntotal,nlayer))
      if(.not. allocated(pao2       )) allocate(pao2       (ntotal,nlayer))
      if(.not. allocated(pao3       )) allocate(pao3       (ntotal,nlayer))
      if(.not. allocated(taugas     )) allocate(taugas     (ntotal,nlayer))
      if(.not. allocated(paray      )) allocate(paray      (ntotal,nlayer))
      if(.not. allocated(tauaer     )) allocate(tauaer     (ntotal,nlayer))
      if(.not. allocated(taul       )) allocate(taul       (ntotal,nlayer))
      if(.not. allocated(taucld     )) allocate(taucld     (ntotal,nlayer))
      if(.not. allocated(wcld       )) allocate(wcld       (ntotal,nlayer))
      if(.not. allocated(gcld       )) allocate(gcld       (ntotal,nlayer))

      if(.not. allocated(taucldlw   )) allocate(taucldlw   (ntotal,nlayer))
      if(.not. allocated(wolc       )) allocate(wolc       (ntotal,nlayer))
      if(.not. allocated(gl         )) allocate(gl         (ntotal,nlayer))
      if(.not. allocated(taucldice  )) allocate(taucldice  (ntotal,nlayer))
      if(.not. allocated(woice      )) allocate(woice      (ntotal,nlayer))
      if(.not. allocated(gice       )) allocate(gice       (ntotal,nlayer))

      if(.not. allocated(w0         )) allocate(w0         (ntotal,nlayer))
      if(.not. allocated(g0         )) allocate(g0         (ntotal,nlayer))
      if(.not. allocated(opd        )) allocate(opd        (ntotal,nlayer))
      if(.not. allocated(uopd       )) allocate(uopd       (ntotal,nlayer))
      if(.not. allocated(ptemp      )) allocate(ptemp      (ntotal,nlayer))
      if(.not. allocated(slope      )) allocate(slope      (ntotal,nlayer))
      if(.not. allocated(b1         )) allocate(b1         (ntotal,nlayer))
      if(.not. allocated(b2         )) allocate(b2         (ntotal,nlayer))
      if(.not. allocated(b3         )) allocate(b3         (ntotal,nlayer))
      if(.not. allocated(ak         )) allocate(ak         (ntotal,nlayer))
      if(.not. allocated(gami       )) allocate(gami       (ntotal,nlayer))
      if(.not. allocated(ee1        )) allocate(ee1        (ntotal,nlayer))
      if(.not. allocated(el1        )) allocate(el1        (ntotal,nlayer))
      if(.not. allocated(em1        )) allocate(em1        (ntotal,nlayer))
      if(.not. allocated(el2        )) allocate(el2        (ntotal,nlayer))
      if(.not. allocated(em2        )) allocate(em2        (ntotal,nlayer))
      if(.not. allocated(af         )) allocate(af         (ntotal,ndbl  ))
      if(.not. allocated(bf         )) allocate(bf         (ntotal,ndbl  ))
      if(.not. allocated(ef         )) allocate(ef         (ntotal,ndbl  ))
      if(.not. allocated(cp         )) allocate(cp         (ntotal,nlayer))
      if(.not. allocated(cpb        )) allocate(cpb        (ntotal,nlayer))
      if(.not. allocated(cmb        )) allocate(cmb        (ntotal,nlayer))
      if(.not. allocated(ck1        )) allocate(ck1        (ntotal,nlayer))
      if(.not. allocated(ck2        )) allocate(ck2        (ntotal,nlayer))
      if(.not. allocated(fnet       )) allocate(fnet       (ntotal,nlayer))
      if(.not. allocated(tmi        )) allocate(tmi        (ntotal,nlayer))
      if(.not. allocated(tmid       )) allocate(tmid       (ntotal,nlayer))
      if(.not. allocated(tmiu       )) allocate(tmiu       (ntotal,nlayer))
      if(.not. allocated(direc      )) allocate(direc      (ntotal,nlayer))
      if(.not. allocated(directu    )) allocate(directu    (ntotal,nlayer))

      !----- 4D variables. ----------------------------------------------------------------!
      if(.not. allocated(pc         )) allocate(pc      (m1,nbin,nelem       ))
      if(.not. allocated(pc_aerad   )) allocate(pc_aerad(nz_rad,nbin,ngroup  ))
      if(.not. allocated(caer       )) allocate(caer    (nlayer,nrad,ngroup  ))
      if(.not. allocated(y3         )) allocate(y3      (ntotal,ngauss,nlayer))

      !----- Initialise the variables recently allocated. ---------------------------------!
      call flush_carma()

      return
   end subroutine init_carma
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine will initialise all scratch variables with zero.                  !
   !---------------------------------------------------------------------------------------!
   subroutine flush_carma()
      implicit none

      p_surf       = 0.0
      p_top        = 0.0
      t_surf       = 0.0
      tabove_aerad = 0.0
      lla          = 0
      lls          = 0
      solnet       = 0.0
      xirdown      = 0.0

      !----- 2D variables. ----------------------------------------------------------------!
      p            = 0.0
      t            = 0.0
      t_aerad      = 0.0
      p_aerad      = 0.0
      qv_aerad     = 0.0

      LWL_aerad    = 0.0
      IWL_aerad    = 0.0
      LWP_aerad    = 0.0
      IWP_aerad    = 0.0
      xland_aerad  = 0.0
      RAIN_aerad   = 0.0

      press        = 0.0
      dpg          = 0.0
      tt           = 0.0
      rdh2o        = 0.0
      rhoa         = 0.0
      rsfx         = 0.0
      heats_aerad  = 0.0
      heati_aerad  = 0.0
      ptempg       = 0.0
      ptempt       = 0.0
      u1i          = 0.0

      !----- 3D variables. ----------------------------------------------------------------!
      gc           = 0.0
      pah2o        = 0.0
      paco2        = 0.0
      pao2         = 0.0
      pao3         = 0.0
      taugas       = 0.0
      paray        = 0.0
      tauaer       = 0.0
      taul         = 0.0
      taucld       = 0.0
      taucldlw     = 0.0
      wolc         = 0.0
      gl           = 0.0
      taucldice    = 0.0
      woice        = 0.0
      gice         = 0.0
      wcld         = 0.0
      gcld         = 0.0
      w0           = 0.0
      g0           = 0.0
      opd          = 0.0
      uopd         = 0.0
      ptemp        = 0.0
      slope        = 0.0
      b1           = 0.0
      b2           = 0.0
      b3           = 0.0
      ak           = 0.0
      gami         = 0.0
      ee1          = 0.0
      el1          = 0.0
      em1          = 0.0
      el2          = 0.0
      em2          = 0.0
      af           = 0.0
      bf           = 0.0
      ef           = 0.0
      cp           = 0.0
      cpb          = 0.0
      cmb          = 0.0
      ck1          = 0.0
      ck2          = 0.0
      fnet         = 0.0
      tmi          = 0.0
      tmid         = 0.0
      tmiu         = 0.0
      direc        = 0.0
      directu      = 0.0
      !----- 4D Variables. ----------------------------------------------------------------!
      pc           = 0.0
      pc_aerad     = 0.0
      caer         = 0.0
      y3           = 0.0

      return 
   end subroutine flush_carma
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   subroutine end_carma()
      implicit none


      !----- 3D variables. ----------------------------------------------------------------!
      if (allocated(p           )) deallocate(p           )
      if (allocated(t           )) deallocate(t           )
      if (allocated(t_aerad     )) deallocate(t_aerad     )
      if (allocated(p_aerad     )) deallocate(p_aerad     )
      if (allocated(qv_aerad    )) deallocate(qv_aerad    )
      if (allocated(LWL_aerad   )) deallocate(LWL_aerad   )
      if (allocated(IWL_aerad   )) deallocate(IWL_aerad   )
      if (allocated(LWP_aerad   )) deallocate(LWP_aerad   )
      if (allocated(IWP_aerad   )) deallocate(IWP_aerad   )
      if (allocated(rhoa        )) deallocate(rhoa        )

      if (allocated(press       )) deallocate(press       )
      if (allocated(dpg         )) deallocate(dpg         )
      if (allocated(tt          )) deallocate(tt          )
      if (allocated(rdh2o       )) deallocate(rdh2o       )
      if (allocated(rsfx        )) deallocate(rsfx        )
      if (allocated(heats_aerad )) deallocate(heats_aerad )
      if (allocated(heati_aerad )) deallocate(heati_aerad )
      if (allocated(ptempg      )) deallocate(ptempg      )
      if (allocated(ptempt      )) deallocate(ptempt      )
      if (allocated(u1i         )) deallocate(u1i         )

      !----- 4D variables. ----------------------------------------------------------------!
      if (allocated(gc          )) deallocate(gc          )
      if (allocated(pah2o       )) deallocate(pah2o       )
      if (allocated(paco2       )) deallocate(paco2       )
      if (allocated(pao2        )) deallocate(pao2        )
      if (allocated(pao3        )) deallocate(pao3        )
      if (allocated(taugas      )) deallocate(taugas      )
      if (allocated(paray       )) deallocate(paray       )
      if (allocated(tauaer      )) deallocate(tauaer      )
      if (allocated(taul        )) deallocate(taul        )
      if (allocated(taucld      )) deallocate(taucld      )
      if (allocated(wcld        )) deallocate(wcld        )
      if (allocated(gcld        )) deallocate(gcld        )

      if (allocated(taucldlw    )) deallocate(taucldlw    )
      if (allocated(wolc        )) deallocate(wolc        )
      if (allocated(gl          )) deallocate(gl          )
      if (allocated(taucldice   )) deallocate(taucldice   )
      if (allocated(woice       )) deallocate(woice       )
      if (allocated(gice        )) deallocate(gice        )

      if (allocated(w0          )) deallocate(w0          )
      if (allocated(g0          )) deallocate(g0          )
      if (allocated(opd         )) deallocate(opd         )
      if (allocated(uopd        )) deallocate(uopd        )
      if (allocated(ptemp       )) deallocate(ptemp       )
      if (allocated(slope       )) deallocate(slope       )
      if (allocated(b1          )) deallocate(b1          )
      if (allocated(b2          )) deallocate(b2          )
      if (allocated(b3          )) deallocate(b3          )
      if (allocated(ak          )) deallocate(ak          )
      if (allocated(gami        )) deallocate(gami        )
      if (allocated(ee1         )) deallocate(ee1         )
      if (allocated(el1         )) deallocate(el1         )
      if (allocated(em1         )) deallocate(em1         )
      if (allocated(el2         )) deallocate(el2         )
      if (allocated(em2         )) deallocate(em2         )
      if (allocated(af          )) deallocate(af          )
      if (allocated(bf          )) deallocate(bf          )
      if (allocated(ef          )) deallocate(ef          )
      if (allocated(cp          )) deallocate(cp          )
      if (allocated(cpb         )) deallocate(cpb         )
      if (allocated(cmb         )) deallocate(cmb         )
      if (allocated(ck1         )) deallocate(ck1         )
      if (allocated(ck2         )) deallocate(ck2         )
      if (allocated(fnet        )) deallocate(fnet        )
      if (allocated(tmi         )) deallocate(tmi         )
      if (allocated(tmid        )) deallocate(tmid        )
      if (allocated(tmiu        )) deallocate(tmiu        )
      if (allocated(direc       )) deallocate(direc       )
      if (allocated(directu     )) deallocate(directu     )

      !----- 5D variables. ----------------------------------------------------------------!
      if (allocated(pc          )) deallocate(pc          )
      if (allocated(pc_aerad    )) deallocate(pc_aerad    )
      if (allocated(caer        )) deallocate(caer        )
      if (allocated(y3          )) deallocate(y3          )
      
      return
   end subroutine end_carma
   !=======================================================================================!
   !=======================================================================================!
end module mem_carma
!==========================================================================================!
!==========================================================================================!
