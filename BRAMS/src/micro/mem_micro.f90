!====================================== Change Log ========================================!
! 5.0.0                                                                                    !
!                                                                                          !
!==========================================================================================!
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved                       !
!  Regional Atmospheric Modeling System - RAMS                                             !
!==========================================================================================!
!==========================================================================================!
module mem_micro
   !=======================================================================================!
   !=======================================================================================!
   type micro_vars
      !----- Variables to be dimensioned by (nzp,nxp,nyp) ---------------------------------!
      real, pointer, dimension(:,:,:) :: rcp,rrp,rpp,rsp,rap,rgp,rhp                       &
                                        ,ccp,crp,cpp,csp,cap,cgp,chp                       &
                                        ,cccnp,cifnp,q2,q6,q7

      !----- Variables to be dimensioned by (nxp,nyp) -------------------------------------!
      real, pointer, dimension(  :,:) :: accpr,accpp,accps,accpa,accpg,accph               &
                                        ,pcprr,pcprp,pcprs,pcpra,pcprg,pcprh               &
                                        ,pcpg,qpcpg,dpcpg
   end type micro_vars
   
   !----- Structure -----------------------------------------------------------------------!
   type (micro_vars), allocatable :: micro_g(:), microm_g(:)
                          
   
   contains
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine performs the allocation of the micro_vars components for each      !
   ! grid, avoiding unecessary allocations based on the user input.                        !
   !---------------------------------------------------------------------------------------!
   subroutine alloc_micro(micro,n1,n2,n3,ng,isfcl)
      use therm_lib, only : cloud_on,bulk_on
      use micphys  , only : icloud,irain,ipris,isnow,iaggr,igraup,ihail

      implicit none          
      type (micro_vars) :: micro
      integer, intent(in) :: n1,n2,n3,ng,isfcl

      !----- Cloud should be allocated whenever cloud is on -------------------------------!
      if (cloud_on) allocate (micro%rcp(n1,n2,n3))

      !------------------------------------------------------------------------------------!
      !    Other hydrometeors will allocated only when the bulk microphysics is activated. !
      !------------------------------------------------------------------------------------!
      if (bulk_on) then

         !----- Rain ----------------------------------------------------------------------!
         if(irain >= 1)  then
            allocate (micro%rrp(n1,n2,n3))
            allocate (micro%accpr(n2,n3))
            allocate (micro%pcprr(n2,n3))
            allocate (micro%q2(n1,n2,n3))
         end if
         !---------------------------------------------------------------------------------!

         !----- Pristine ice --------------------------------------------------------------!
         if(ipris >= 1)  then
            allocate (micro%rpp(n1,n2,n3))
            allocate (micro%accpp(n2,n3))
            allocate (micro%pcprp(n2,n3))
         end if
         !---------------------------------------------------------------------------------!

         !----- Snow ----------------------------------------------------------------------!
         if(isnow >= 1)  then
            allocate (micro%rsp(n1,n2,n3))
            allocate (micro%accps(n2,n3))
            allocate (micro%pcprs(n2,n3))
         end if
         !---------------------------------------------------------------------------------!

         !----- Aggregates ----------------------------------------------------------------!
         if(iaggr >= 1)  then
            allocate (micro%rap(n1,n2,n3))
            allocate (micro%accpa(n2,n3))
            allocate (micro%pcpra(n2,n3))
         end if
         !---------------------------------------------------------------------------------!

         !----- Graupel -------------------------------------------------------------------!
         if(igraup >= 1) then
            allocate (micro%rgp(n1,n2,n3))
            allocate (micro%accpg(n2,n3))
            allocate (micro%pcprg(n2,n3))
            allocate (micro%q6(n1,n2,n3))
         end if
         !---------------------------------------------------------------------------------!

         !----- Hail ----------------------------------------------------------------------!
         if(ihail >= 1)  then
            allocate (micro%rhp(n1,n2,n3))
            allocate (micro%accph(n2,n3))
            allocate (micro%pcprh(n2,n3))
            allocate (micro%q7(n1,n2,n3))
         end if
         !---------------------------------------------------------------------------------!

         !----- Allocating concentration for prognostic hydrometeors ----------------------!
         if(icloud == 5) allocate (micro%ccp(n1,n2,n3)) !----- Cloud        ---------------!
         if(irain  == 5) allocate (micro%crp(n1,n2,n3)) !----- Rain         ---------------!
         if(ipris  == 5) allocate (micro%cpp(n1,n2,n3)) !----- Pristine ice ---------------!
         if(isnow  == 5) allocate (micro%csp(n1,n2,n3)) !----- Snow         ---------------!
         if(iaggr  == 5) allocate (micro%cap(n1,n2,n3)) !----- Aggregates   ---------------!
         if(igraup == 5) allocate (micro%cgp(n1,n2,n3)) !----- Graupel      ---------------!
         if(ihail  == 5) allocate (micro%chp(n1,n2,n3)) !----- Hail         ---------------!

         !----- Allocating CCN concentrations for prognostic CCNs -------------------------!
         if (icloud == 7) allocate (micro%cccnp(n1,n2,n3))
         if (ipris  == 7) allocate (micro%cifnp(n1,n2,n3))
         !---------------------------------------------------------------------------------!

         !----- Allocating the precipitation rates and alikes -----------------------------!
         allocate (micro%pcpg (n2,n3))
         allocate (micro%qpcpg(n2,n3))
         allocate (micro%dpcpg(n2,n3))
         !---------------------------------------------------------------------------------!

      end if

      return
   end subroutine alloc_micro
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine nullifies all pointers for a safe allocation.                      !
   !---------------------------------------------------------------------------------------!
   subroutine nullify_micro(micro)

      implicit none
      !----- Argument ---------------------------------------------------------------------!
      type (micro_vars) :: micro
      !------------------------------------------------------------------------------------!
   
      if (associated(micro%rcp   )) nullify (micro%rcp   )
      if (associated(micro%rrp   )) nullify (micro%rrp   )
      if (associated(micro%rpp   )) nullify (micro%rpp   )
      if (associated(micro%rsp   )) nullify (micro%rsp   )
      if (associated(micro%rap   )) nullify (micro%rap   )
      if (associated(micro%rgp   )) nullify (micro%rgp   )
      if (associated(micro%rhp   )) nullify (micro%rhp   )
      if (associated(micro%ccp   )) nullify (micro%ccp   )
      if (associated(micro%crp   )) nullify (micro%crp   )
      if (associated(micro%cpp   )) nullify (micro%cpp   )
      if (associated(micro%csp   )) nullify (micro%csp   )
      if (associated(micro%cap   )) nullify (micro%cap   )
      if (associated(micro%cgp   )) nullify (micro%cgp   )
      if (associated(micro%chp   )) nullify (micro%chp   )
      if (associated(micro%cccnp )) nullify (micro%cccnp )
      if (associated(micro%cifnp )) nullify (micro%cifnp )
      if (associated(micro%q2    )) nullify (micro%q2    )
      if (associated(micro%q6    )) nullify (micro%q6    )
      if (associated(micro%q7    )) nullify (micro%q7    )

      if (associated(micro%accpr )) nullify (micro%accpr )
      if (associated(micro%accpp )) nullify (micro%accpp )
      if (associated(micro%accps )) nullify (micro%accps )
      if (associated(micro%accpa )) nullify (micro%accpa )
      if (associated(micro%accpg )) nullify (micro%accpg )
      if (associated(micro%accph )) nullify (micro%accph )
      if (associated(micro%pcprr )) nullify (micro%pcprr )
      if (associated(micro%pcprp )) nullify (micro%pcprp )
      if (associated(micro%pcprs )) nullify (micro%pcprs )
      if (associated(micro%pcpra )) nullify (micro%pcpra )
      if (associated(micro%pcprg )) nullify (micro%pcprg )
      if (associated(micro%pcprh )) nullify (micro%pcprh )
      if (associated(micro%pcpg  )) nullify (micro%pcpg  )
      if (associated(micro%qpcpg )) nullify (micro%qpcpg )
      if (associated(micro%dpcpg )) nullify (micro%dpcpg )

      return
   end subroutine nullify_micro
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine deallocates all pointers when they are no longer needed.           !
   !---------------------------------------------------------------------------------------!
   subroutine dealloc_micro(micro)

      implicit none
      !----- Argument ---------------------------------------------------------------------!
      type (micro_vars) :: micro
      !------------------------------------------------------------------------------------!
   
      if (associated(micro%rcp   )) deallocate (micro%rcp   )
      if (associated(micro%rrp   )) deallocate (micro%rrp   )
      if (associated(micro%rpp   )) deallocate (micro%rpp   )
      if (associated(micro%rsp   )) deallocate (micro%rsp   )
      if (associated(micro%rap   )) deallocate (micro%rap   )
      if (associated(micro%rgp   )) deallocate (micro%rgp   )
      if (associated(micro%rhp   )) deallocate (micro%rhp   )
      if (associated(micro%ccp   )) deallocate (micro%ccp   )
      if (associated(micro%crp   )) deallocate (micro%crp   )
      if (associated(micro%cpp   )) deallocate (micro%cpp   )
      if (associated(micro%csp   )) deallocate (micro%csp   )
      if (associated(micro%cap   )) deallocate (micro%cap   )
      if (associated(micro%cgp   )) deallocate (micro%cgp   )
      if (associated(micro%chp   )) deallocate (micro%chp   )
      if (associated(micro%cccnp )) deallocate (micro%cccnp )
      if (associated(micro%cifnp )) deallocate (micro%cifnp )
      if (associated(micro%q2    )) deallocate (micro%q2    )
      if (associated(micro%q6    )) deallocate (micro%q6    )
      if (associated(micro%q7    )) deallocate (micro%q7    )

      if (associated(micro%accpr )) deallocate (micro%accpr )
      if (associated(micro%accpp )) deallocate (micro%accpp )
      if (associated(micro%accps )) deallocate (micro%accps )
      if (associated(micro%accpa )) deallocate (micro%accpa )
      if (associated(micro%accpg )) deallocate (micro%accpg )
      if (associated(micro%accph )) deallocate (micro%accph )
      if (associated(micro%pcprr )) deallocate (micro%pcprr )
      if (associated(micro%pcprp )) deallocate (micro%pcprp )
      if (associated(micro%pcprs )) deallocate (micro%pcprs )
      if (associated(micro%pcpra )) deallocate (micro%pcpra )
      if (associated(micro%pcprg )) deallocate (micro%pcprg )
      if (associated(micro%pcprh )) deallocate (micro%pcprh )
      if (associated(micro%pcpg  )) deallocate (micro%pcpg  )
      if (associated(micro%qpcpg )) deallocate (micro%qpcpg )
      if (associated(micro%dpcpg )) deallocate (micro%dpcpg )

      return
   end subroutine dealloc_micro
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine deallocates all pointers when they are no longer needed.           !
   !---------------------------------------------------------------------------------------!
   subroutine zero_micro(micro)

      implicit none
      !----- Argument ---------------------------------------------------------------------!
      type (micro_vars) :: micro
      !------------------------------------------------------------------------------------!
   
      if (associated(micro%rcp   )) micro%rcp   = 0.0
      if (associated(micro%rrp   )) micro%rrp   = 0.0
      if (associated(micro%rpp   )) micro%rpp   = 0.0
      if (associated(micro%rsp   )) micro%rsp   = 0.0
      if (associated(micro%rap   )) micro%rap   = 0.0
      if (associated(micro%rgp   )) micro%rgp   = 0.0
      if (associated(micro%rhp   )) micro%rhp   = 0.0
      if (associated(micro%ccp   )) micro%ccp   = 0.0
      if (associated(micro%crp   )) micro%crp   = 0.0
      if (associated(micro%cpp   )) micro%cpp   = 0.0
      if (associated(micro%csp   )) micro%csp   = 0.0
      if (associated(micro%cap   )) micro%cap   = 0.0
      if (associated(micro%cgp   )) micro%cgp   = 0.0
      if (associated(micro%chp   )) micro%chp   = 0.0
      if (associated(micro%cccnp )) micro%cccnp = 0.0
      if (associated(micro%cifnp )) micro%cifnp = 0.0
      if (associated(micro%q2    )) micro%q2    = 0.0
      if (associated(micro%q6    )) micro%q6    = 0.0
      if (associated(micro%q7    )) micro%q7    = 0.0

      if (associated(micro%accpr )) micro%accpr = 0.0
      if (associated(micro%accpp )) micro%accpp = 0.0
      if (associated(micro%accps )) micro%accps = 0.0
      if (associated(micro%accpa )) micro%accpa = 0.0
      if (associated(micro%accpg )) micro%accpg = 0.0
      if (associated(micro%accph )) micro%accph = 0.0
      if (associated(micro%pcprr )) micro%pcprr = 0.0
      if (associated(micro%pcprp )) micro%pcprp = 0.0
      if (associated(micro%pcprs )) micro%pcprs = 0.0
      if (associated(micro%pcpra )) micro%pcpra = 0.0
      if (associated(micro%pcprg )) micro%pcprg = 0.0
      if (associated(micro%pcprh )) micro%pcprh = 0.0
      if (associated(micro%pcpg  )) micro%pcpg  = 0.0
      if (associated(micro%qpcpg )) micro%qpcpg = 0.0
      if (associated(micro%dpcpg )) micro%dpcpg = 0.0

      return
   end subroutine zero_micro
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine fills the variable table information for each microphysics vari-   !
   ! able.                                                                                 !
   !---------------------------------------------------------------------------------------!
   subroutine filltab_micro(micro,microm,imean,n1,n2,n3,ng)

      use var_tables
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type (micro_vars), intent(inout) :: micro,microm
      integer          , intent(in)    :: imean,n1,n2,n3,ng
      !----- Local variables --------------------------------------------------------------!
      integer                          :: npts
      real, pointer                    :: var,varm
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !                      ****** 3-D variables (NZP,NXP,NYP) ******                     !
      !------------------------------------------------------------------------------------!
      npts=n1*n2*n3
      if (associated(micro%rcp))                                                           &
         call vtables2 (micro%rcp,microm%rcp,ng, npts, imean                               &
                       ,'RCP :3:hist:anal:mpti:mpt3:mpt1')

      if (associated(micro%rrp))                                                           &
         call vtables2 (micro%rrp,microm%rrp,ng, npts, imean                               &
                       ,'RRP :3:hist:anal:mpti:mpt3:mpt1')

      if (associated(micro%rpp))                                                           &
         call vtables2 (micro%rpp,microm%rpp,ng, npts, imean                               &
                       ,'RPP :3:hist:anal:mpti:mpt3:mpt1')

      if (associated(micro%rsp))                                                           &
         call vtables2 (micro%rsp,microm%rsp,ng, npts, imean                               &
                       ,'RSP :3:hist:anal:mpti:mpt3:mpt1')

      if (associated(micro%rap))                                                           &
         call vtables2 (micro%rap,microm%rap,ng, npts, imean                               &
                       ,'RAP :3:hist:anal:mpti:mpt3:mpt1')

      if (associated(micro%rgp))                                                           &
         call vtables2 (micro%rgp,microm%rgp,ng, npts, imean                               &
                       ,'RGP :3:hist:anal:mpti:mpt3:mpt1')

      if (associated(micro%rhp))                                                           &
         call vtables2 (micro%rhp,microm%rhp,ng, npts, imean                               &
                       ,'RHP :3:hist:anal:mpti:mpt3:mpt1')

      if (associated(micro%ccp))                                                           &
         call vtables2 (micro%ccp,microm%ccp,ng, npts, imean                               &
                       ,'CCP :3:hist:anal:mpti:mpt3:mpt1')

      if (associated(micro%crp))                                                           &
         call vtables2 (micro%crp,microm%crp,ng, npts, imean                               &
                       ,'CRP :3:hist:anal:mpti:mpt3:mpt1')

      if (associated(micro%cpp))                                                           &
         call vtables2 (micro%cpp,microm%cpp,ng, npts, imean                               &
                       ,'CPP :3:hist:anal:mpti:mpt3:mpt1')

      if (associated(micro%csp))                                                           &
         call vtables2 (micro%csp,microm%csp,ng, npts, imean                               &
                       ,'CSP :3:hist:anal:mpti:mpt3:mpt1')

      if (associated(micro%cap))                                                           &
         call vtables2 (micro%cap,microm%cap,ng, npts, imean                               &
                       ,'CAP :3:hist:anal:mpti:mpt3:mpt1')

      if (associated(micro%cgp))                                                           &
         call vtables2 (micro%cgp,microm%cgp,ng, npts, imean                               &
                       ,'CGP :3:hist:anal:mpti:mpt3:mpt1')

      if (associated(micro%chp))                                                           &
         call vtables2 (micro%chp,microm%chp,ng, npts, imean                               &
                       ,'CHP :3:hist:anal:mpti:mpt3:mpt1')

      if (associated(micro%cccnp))                                                         &
         call vtables2 (micro%cccnp,microm%cccnp,ng, npts, imean                           &
                       ,'CCCNP :3:hist:anal:mpti:mpt3:mpt1')

      if (associated(micro%cifnp))                                                         &
         call vtables2 (micro%cifnp,microm%cifnp,ng, npts, imean                           &
                       ,'CIFNP :3:hist:anal:mpti:mpt3:mpt1')

      if (associated(micro%q2))                                                            &
         call vtables2 (micro%q2,microm%q2,ng, npts, imean                                 &
                       ,'Q2 :3:hist:anal:mpti:mpt3:mpt1')

      if (associated(micro%q6))                                                            &
         call vtables2 (micro%q6,microm%q6,ng, npts, imean                                 &
                       ,'Q6 :3:hist:anal:mpti:mpt3:mpt1')

      if (associated(micro%q7))                                                            &
         call vtables2 (micro%q7,microm%q7,ng, npts, imean                                 &
                       ,'Q7 :3:hist:anal:mpti:mpt3:mpt1')
      !------------------------------------------------------------------------------------!


                    
      !------------------------------------------------------------------------------------!
      !                      ****** 2-D variables (NZP,NXP,NYP) ******                     !
      !------------------------------------------------------------------------------------!
      npts=n2*n3

      if (associated(micro%accpr))                                                         &
         call vtables2 (micro%accpr,microm%accpr,ng, npts, imean                           &
                       ,'ACCPR :2:hist:anal:mpti:mpt3')

      if (associated(micro%accpp))                                                         &
         call vtables2 (micro%accpp,microm%accpp,ng, npts, imean                           &
                       ,'ACCPP :2:hist:anal:mpti:mpt3')

      if (associated(micro%accps))                                                         &
         call vtables2 (micro%accps,microm%accps,ng, npts, imean                           &
                       ,'ACCPS :2:hist:anal:mpti:mpt3')

      if (associated(micro%accpa))                                                         &
         call vtables2 (micro%accpa,microm%accpa,ng, npts, imean                           &
                       ,'ACCPA :2:hist:anal:mpti:mpt3')

      if (associated(micro%accpg))                                                         &
         call vtables2 (micro%accpg,microm%accpg,ng, npts, imean                           &
                       ,'ACCPG :2:hist:anal:mpti:mpt3')

      if (associated(micro%accph))                                                         &
         call vtables2 (micro%accph,microm%accph,ng, npts, imean                           &
                       ,'ACCPH :2:hist:anal:mpti:mpt3')

      if (associated(micro%pcprr))                                                         &
         call vtables2 (micro%pcprr,microm%pcprr,ng, npts, imean                           &
                       ,'PCPRR :2:anal:mpt3')

      if (associated(micro%pcprp))                                                         &
         call vtables2 (micro%pcprp,microm%pcprp,ng, npts, imean                           &
                       ,'PCPRP :2:anal:mpt3')

      if (associated(micro%pcprs))                                                         &
         call vtables2 (micro%pcprs,microm%pcprs,ng, npts, imean                           &
                       ,'PCPRS :2:anal:mpt3')

      if (associated(micro%pcpra))                                                         &
         call vtables2 (micro%pcpra,microm%pcpra,ng, npts, imean                           &
                       ,'PCPRA :2:anal:mpt3')

      if (associated(micro%pcprg))                                                         &
         call vtables2 (micro%pcprg,microm%pcprg,ng, npts, imean                           &
                       ,'PCPRG :2:anal:mpt3')

      if (associated(micro%pcprh))                                                         &
         call vtables2 (micro%pcprh,microm%pcprh,ng, npts, imean                           &
                       ,'PCPRH :2:anal:mpt3')

      if (associated(micro%pcpg))                                                          &
         call vtables2 (micro%pcpg,microm%pcpg,ng, npts, imean                             &
                       ,'PCPG :2:hist:mpti:mpt3')

      if (associated(micro%qpcpg))                                                         &
         call vtables2 (micro%qpcpg,microm%qpcpg,ng, npts, imean                           &
                       ,'QPCPG :2:hist:mpti:mpt3')

      if (associated(micro%dpcpg))                                                         &
         call vtables2 (micro%dpcpg,microm%dpcpg,ng, npts, imean                           &
                       ,'DPCPG :2:hist:mpti:mpt3')
      return
   end subroutine filltab_micro
   !=======================================================================================!
   !=======================================================================================!
end module mem_micro
!==========================================================================================!
!==========================================================================================!
