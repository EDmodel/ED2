!===================================== Change Log =========================================!
! 5.0.0                                                                                    !
!                                                                                          !
! MLO - 08/31/08 Rearranging the cumulus structure, eliminating grell and shcu structures. !
!                                                                                          !
!==========================================================================================!
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved                       !
!  Regional Atmospheric Modeling System - RAMS                                             !
!==========================================================================================!
!    This module contains various variables to be used as scratch arrays.                  !
!------------------------------------------------------------------------------------------!
module mem_scratch

   use grid_dims, only: maxdimp

   type scratch_vars
      !----- The largest arrays available, they should hold the largest possible variable -!
      real, pointer, dimension(:) ::  scr1,scr2,scr3
      !----- 2-D variables, they should hold any (X,Y) variable for any grid. -------------!
      real, pointer, dimension(:) ::  vt2da,vt2db,vt2dc,vt2dd,vt2de,vt2df
      !----- 3-D variables, for (Z,X,Y), (X,Y,P), (X,Y,C), and (X,Y,W) variables. ---------!
      real, pointer, dimension(:) ::  vt3da,vt3db,vt3dc,vt3dd,vt3de,vt3df,vt3dg,vt3dh
      real, pointer, dimension(:) ::  vt3di,vt3dj,vt3dk,vt3dl,vt3dm,vt3dn,vt3do,vt3dp,vt3dq
      real, pointer, dimension(:) ::  vt3dr,vt3ds
      !----- 4-D variables, for (Z,X,Y,C), (G,X,Y,P), and (S,X,Y,P) variables -------------!
      real, pointer, dimension(:) ::  vt4da,vt4db,vt4dc
   end type scratch_vars

   type (scratch_vars) :: scratch

   !------A bunch of 1-D vectors ----------------------------------------------------------!
   real   , dimension(maxdimp) ::  vctr1,  vctr2,  vctr3,  vctr4,  vctr5,  vctr6
   real   , dimension(maxdimp) ::  vctr7,  vctr8,  vctr9, vctr10, vctr11, vctr12 
   real   , dimension(maxdimp) :: vctr13, vctr14, vctr15, vctr16, vctr17, vctr18
   real   , dimension(maxdimp) :: vctr19, vctr20, vctr21, vctr22, vctr23, vctr24
   real   , dimension(maxdimp) :: vctr25, vctr26, vctr27, vctr28, vctr29, vctr30
    
   real   , dimension(maxdimp) :: vctr31, vctr32, vctr33, vctr34, vctr35, vctr36
   real   , dimension(maxdimp) :: vctr37, vctr38, vctr39, vctr40, vctr41
   integer, dimension(maxdimp) :: ivctr
   !---------------------------------------------------------------------------------------!

   contains
   !=======================================================================================!
   !=======================================================================================!
   subroutine alloc_scratch(ngrs,nmzp,nmxp,nmyp,nnzp,nnxp,nnyp,nzg,nzs,npatch,nclouds      &
                           ,maxx,maxy,maxz)

      use mem_aerad,    only : nwave    ! ! intent(in)
      use mem_radiate , only : ilwrtyp  & ! intent(in)
                             , iswrtyp  ! ! intent(in)
      use catt_start  , only : CATT     ! ! intent(in)
      implicit none
      !------ Arguments -------------------------------------------------------------------!
      integer                  , intent(in)  :: ngrs,nzg,nzs,npatch,nclouds
      integer, dimension (ngrs), intent(in)  :: nmzp,nmxp,nmyp,nnzp,nnxp,nnyp
      integer                  , intent(out) :: maxx,maxy,maxz
      !------ Local variables -------------------------------------------------------------!
      integer                                :: ng,ntpts,ntpts2,ntpts4,ntptsx
      integer                                :: ntpts_catt
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !    Find the maximum number of grid points needed for any grid. The max points in   !
      ! each direction are passed back for use by various nesting things.                  !
      !------------------------------------------------------------------------------------!
      maxx       = 0
      maxy       = 0
      maxz       = 0
      ntpts      = 0
      ntpts2     = 0
      ntpts_catt = 0
      ntpts4     = 0
      do ng=1,ngrs
         maxx   = max(maxx,nnxp(ng))
         maxy   = max(maxy,nnyp(ng))
         maxz   = max(maxz,nnzp(ng))
         ntpts  = max(nmxp(ng)*nmyp(ng)*nmzp(ng),nmxp(ng)*nmyp(ng)*nclouds                 &
                     ,nmxp(ng)*nmyp(ng)*npatch,ntpts )
         ntpts2 = max(nmxp(ng)*nmyp(ng),ntpts2 )
         ntpts4 = max(nmxp(ng)*nmyp(ng)*max(nmzp(ng),nzg,nzs,nwave)*max(npatch,nclouds)    &
                     ,ntpts4)
      enddo
      
      !----- Setting the maximum possible number for scr1 and scr2. -----------------------!
      ntptsx=max(maxx*maxy*maxz,ntpts4,ntpts2*nzg*npatch,ntpts2*nzs*npatch,maxz*40)+1000

      !----- Making it even larger in case it's a CATT run --------------------------------!
      if (ilwrtyp==4 .or. iswrtyp==4) then
         ntpts_catt = max(ntptsx,maxx*maxy*nwave)+1000
      else
         ntpts_catt = ntptsx
      end if

      !------------------------------------------------------------------------------------!
      !    Allocate arrays based on options (if necessary).                                !
      !------------------------------------------------------------------------------------!
      allocate (scratch%scr1 (ntpts_catt))
      allocate (scratch%scr2 (ntpts_catt))

      if (ilwrtyp==4 .or. iswrtyp==4) then
         allocate (scratch%scr3 (ntpts_catt))
      else
         allocate (scratch%scr3 (1))  ! Not used in this case
      endif

      allocate (scratch%vt3da(ntpts ))
      allocate (scratch%vt3db(ntpts ))

      allocate (scratch%vt3dc(ntpts ))
      allocate (scratch%vt3dd(ntpts ))
      allocate (scratch%vt3de(ntpts ))
      allocate (scratch%vt3df(ntpts ))
      allocate (scratch%vt3dg(ntpts ))
      allocate (scratch%vt3dh(ntpts ))
      allocate (scratch%vt3di(ntpts ))
      allocate (scratch%vt3dj(ntpts ))
      allocate (scratch%vt3dk(ntpts ))
      allocate (scratch%vt3dl(ntpts ))
      allocate (scratch%vt3dm(ntpts ))
      allocate (scratch%vt3dn(ntpts ))
      allocate (scratch%vt3do(ntpts ))
      allocate (scratch%vt3dp(ntpts ))
      allocate (scratch%vt3dq(ntpts ))
      allocate (scratch%vt3dr(ntpts ))
      allocate (scratch%vt3ds(ntpts ))

      allocate (scratch%vt2da(ntpts2))
      allocate (scratch%vt2db(ntpts2))
      allocate (scratch%vt2dc(ntpts2))
      allocate (scratch%vt2dd(ntpts2))
      allocate (scratch%vt2de(ntpts2))
      allocate (scratch%vt2df(ntpts2))


      allocate (scratch%vt4da(ntpts4))
      allocate (scratch%vt4db(ntpts4))
      allocate (scratch%vt4dc(ntpts4))

      ! ALF - Putting zero in all variables
      ! For CATT
      call azero(ntpts_catt, scratch%scr1)
      call azero(ntpts_catt, scratch%scr2)
      if (CATT == 1) then
         call azero(ntpts_catt, scratch%scr3)
      else
         call azero(1, scratch%scr3)  ! Not used in this case
      endif

      call azero(ntpts , scratch%vt3da)
      call azero(ntpts , scratch%vt3db)
      call azero(ntpts , scratch%vt3dc)
      call azero(ntpts , scratch%vt3dd)
      call azero(ntpts , scratch%vt3de)
      call azero(ntpts , scratch%vt3df)
      call azero(ntpts , scratch%vt3dg)
      call azero(ntpts , scratch%vt3dh)
      call azero(ntpts , scratch%vt3di)
      call azero(ntpts , scratch%vt3dj)
      call azero(ntpts , scratch%vt3dk)
      call azero(ntpts , scratch%vt3dl)
      call azero(ntpts , scratch%vt3dm)
      call azero(ntpts , scratch%vt3dn)
      call azero(ntpts , scratch%vt3do)
      call azero(ntpts , scratch%vt3dp)
      call azero(ntpts , scratch%vt3dq)
      call azero(ntpts , scratch%vt3dr)
      call azero(ntpts , scratch%vt3ds)
      call azero(ntpts2, scratch%vt2da)
      call azero(ntpts2, scratch%vt2db)
      call azero(ntpts2, scratch%vt2dc)
      call azero(ntpts2, scratch%vt2dd)
      call azero(ntpts2, scratch%vt2de)
      call azero(ntpts2, scratch%vt2df)
      call azero(ntpts4, scratch%vt4da)
      call azero(ntpts4, scratch%vt4db)
      call azero(ntpts4, scratch%vt4dc)

      return
   end subroutine alloc_scratch
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   subroutine nullify_scratch()
      implicit none

      !----- Deallocate all scratch arrays ------------------------------------------------!
      if (associated(scratch%scr1 ))  nullify (scratch%scr1 )
      if (associated(scratch%scr2 ))  nullify (scratch%scr2 )
      if (associated(scratch%scr3 ))  nullify (scratch%scr3 )
      if (associated(scratch%vt3da))  nullify (scratch%vt3da)
      if (associated(scratch%vt3db))  nullify (scratch%vt3db)
      if (associated(scratch%vt3dc))  nullify (scratch%vt3dc)
      if (associated(scratch%vt3dd))  nullify (scratch%vt3dd)
      if (associated(scratch%vt3de))  nullify (scratch%vt3de)
      if (associated(scratch%vt3df))  nullify (scratch%vt3df)
      if (associated(scratch%vt3dg))  nullify (scratch%vt3dg)
      if (associated(scratch%vt3dh))  nullify (scratch%vt3dh)
      if (associated(scratch%vt3di))  nullify (scratch%vt3di)
      if (associated(scratch%vt3dj))  nullify (scratch%vt3dj)
      if (associated(scratch%vt3dk))  nullify (scratch%vt3dk)
      if (associated(scratch%vt3dl))  nullify (scratch%vt3dl)
      if (associated(scratch%vt3dm))  nullify (scratch%vt3dm)
      if (associated(scratch%vt3dn))  nullify (scratch%vt3dn)
      if (associated(scratch%vt3do))  nullify (scratch%vt3do)
      if (associated(scratch%vt3dp))  nullify (scratch%vt3dp)
      if (associated(scratch%vt3dq))  nullify (scratch%vt3dq)
      if (associated(scratch%vt3dr))  nullify (scratch%vt3dr)
      if (associated(scratch%vt3ds))  nullify (scratch%vt3ds)
      if (associated(scratch%vt2da))  nullify (scratch%vt2da)
      if (associated(scratch%vt2db))  nullify (scratch%vt2db)
      if (associated(scratch%vt2dc))  nullify (scratch%vt2dc)
      if (associated(scratch%vt2dd))  nullify (scratch%vt2dd)
      if (associated(scratch%vt2de))  nullify (scratch%vt2de)
      if (associated(scratch%vt2df))  nullify (scratch%vt2df)
      if (associated(scratch%vt4da))  nullify (scratch%vt4da)
      if (associated(scratch%vt4db))  nullify (scratch%vt4db)
      if (associated(scratch%vt4dc))  nullify (scratch%vt4dc)

      return
   end subroutine nullify_scratch
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   subroutine dealloc_scratch()
      implicit none

      !----- Deallocate all scratch arrays ------------------------------------------------!
      if (associated(scratch%scr1 ))  deallocate (scratch%scr1 )
      if (associated(scratch%scr2 ))  deallocate (scratch%scr2 )
      if (associated(scratch%scr3 ))  deallocate (scratch%scr3 )
      if (associated(scratch%vt3da))  deallocate (scratch%vt3da)
      if (associated(scratch%vt3db))  deallocate (scratch%vt3db)
      if (associated(scratch%vt3dc))  deallocate (scratch%vt3dc)
      if (associated(scratch%vt3dd))  deallocate (scratch%vt3dd)
      if (associated(scratch%vt3de))  deallocate (scratch%vt3de)
      if (associated(scratch%vt3df))  deallocate (scratch%vt3df)
      if (associated(scratch%vt3dg))  deallocate (scratch%vt3dg)
      if (associated(scratch%vt3dh))  deallocate (scratch%vt3dh)
      if (associated(scratch%vt3di))  deallocate (scratch%vt3di)
      if (associated(scratch%vt3dj))  deallocate (scratch%vt3dj)
      if (associated(scratch%vt3dk))  deallocate (scratch%vt3dk)
      if (associated(scratch%vt3dl))  deallocate (scratch%vt3dl)
      if (associated(scratch%vt3dm))  deallocate (scratch%vt3dm)
      if (associated(scratch%vt3dn))  deallocate (scratch%vt3dn)
      if (associated(scratch%vt3do))  deallocate (scratch%vt3do)
      if (associated(scratch%vt3dp))  deallocate (scratch%vt3dp)
      if (associated(scratch%vt3dq))  deallocate (scratch%vt3dq)
      if (associated(scratch%vt3dr))  deallocate (scratch%vt3dr)
      if (associated(scratch%vt3ds))  deallocate (scratch%vt3ds)
      if (associated(scratch%vt2da))  deallocate (scratch%vt2da)
      if (associated(scratch%vt2db))  deallocate (scratch%vt2db)
      if (associated(scratch%vt2dc))  deallocate (scratch%vt2dc)
      if (associated(scratch%vt2dd))  deallocate (scratch%vt2dd)
      if (associated(scratch%vt2de))  deallocate (scratch%vt2de)
      if (associated(scratch%vt2df))  deallocate (scratch%vt2df)
      if (associated(scratch%vt4da))  deallocate (scratch%vt4da)
      if (associated(scratch%vt4db))  deallocate (scratch%vt4db)
      if (associated(scratch%vt4dc))  deallocate (scratch%vt4dc)
      return
   end subroutine dealloc_scratch
   !=======================================================================================!
   !=======================================================================================!
end module mem_scratch
!==========================================================================================!
!==========================================================================================!
