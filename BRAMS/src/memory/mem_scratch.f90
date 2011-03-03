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
      real, pointer, dimension(:) ::  scr1,scr2,scr3,scr4,scr5,scr6
      !----- 2-D variables, they should hold any (X,Y) variable for any grid. -------------!
      real, pointer, dimension(:) ::  vt2da,vt2db,vt2dc,vt2dd,vt2de,vt2df,vt2dg,vt2dh
      real, pointer, dimension(:) ::  vt2di,vt2dj,vt2dk,vt2dl,vt2dm,vt2dn,vt2do,vt2dp
      real, pointer, dimension(:) ::  vt2dq,vt2dr,vt2ds
      !----- 3-D variables, for (Z,X,Y), (X,Y,P), (X,Y,C), and (X,Y,W) variables. ---------!
      real, pointer, dimension(:) ::  vt3da,vt3db,vt3dc,vt3dd,vt3de,vt3df,vt3dg,vt3dh
      real, pointer, dimension(:) ::  vt3di,vt3dj,vt3dk,vt3dl,vt3dm,vt3dn,vt3do,vt3dp
      real, pointer, dimension(:) ::  vt3dq,vt3dr,vt3ds
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

      use mem_aerad,    only : nwave       ! ! intent(in)
      use mem_radiate , only : ilwrtyp     & ! intent(in)
                             , iswrtyp     ! ! intent(in)
      use catt_start  , only : CATT        ! ! intent(in)
      use grid_dims   , only : ndim_types  & ! intent(in)
                             , number_dims ! ! intent(in)
      implicit none
      !------ Arguments -------------------------------------------------------------------!
      integer                        , intent(in)  :: ngrs
      integer                        , intent(in)  :: nzg
      integer                        , intent(in)  :: nzs
      integer                        , intent(in)  :: npatch
      integer                        , intent(in)  :: nclouds
      integer, dimension (ngrs)      , intent(in)  :: nmzp
      integer, dimension (ngrs)      , intent(in)  :: nmxp
      integer, dimension (ngrs)      , intent(in)  :: nmyp
      integer, dimension (ngrs)      , intent(in)  :: nnzp
      integer, dimension (ngrs)      , intent(in)  :: nnxp
      integer, dimension (ngrs)      , intent(in)  :: nnyp
      integer                        , intent(out) :: maxx
      integer                        , intent(out) :: maxy
      integer                        , intent(out) :: maxz
      !------ Local variables -------------------------------------------------------------!
      integer, dimension (ndim_types)              :: npts_local
      integer, dimension (ndim_types)              :: npts_global
      integer                                      :: nptsmax
      integer                                      :: mpts1d
      integer                                      :: mpts2d
      integer                                      :: mpts3d
      integer                                      :: mpts4d
      integer                                      :: ng
      integer                                      :: idim
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !    Find the maximum number of grid points needed for any grid. The max points in   !
      ! each direction are passed back for use by various nesting things.                  !
      !------------------------------------------------------------------------------------!
      maxx       = 0
      maxy       = 0
      maxz       = 0
      do ng=1,ngrs
         maxx   = max(maxx,nnxp(ng))
         maxy   = max(maxy,nnyp(ng))
         maxz   = max(maxz,nnzp(ng))
      end do

      !----- Initialise all dimensions with zero. -----------------------------------------!
      npts_local (:) = 0
      npts_global(:) = 0
      do ng = 1,ngrs
         !----- 1-D. ----------------------------------------------------------------------!
         npts_local (1) = max(npts_local(1),nmxp(ng),nmyp(ng),nmzp(ng),nclouds,npatch      &
                             ,nwave,nzg,nzs)
         npts_global(1) = max(npts_global(1),nnxp(ng),nnyp(ng),nnzp(ng),nclouds,npatch     &
                             ,nwave,nzg,nzs)
         !----- 2-D (nxp,nyp) -------------------------------------------------------------!
         npts_local (2) = max(npts_local(2),nmxp(ng)*nmyp(ng))
         npts_global(2) = max(npts_global(2),nnxp(ng)*nnyp(ng))
         !----- 3-D (nzp,nxp,nyp) ---------------------------------------------------------!
         npts_local (3) = max(npts_local(3),nmzp(ng)*nmxp(ng)*nmyp(ng))
         npts_global(3) = max(npts_global(3),nnzp(ng)*nnxp(ng)*nnyp(ng))
         !----- 4-D (nzg,nxp,nyp,npatch) --------------------------------------------------!
         npts_local (4) = max(npts_local(4),nzg*nmxp(ng)*nmyp(ng)*npatch)
         npts_global(4) = max(npts_global(4),nzg*nnxp(ng)*nnyp(ng)*npatch)
         !----- 4-D (nzs,nxp,nyp,npatch) --------------------------------------------------!
         npts_local (5) = max(npts_local(5),nzs*nmxp(ng)*nmyp(ng)*npatch)
         npts_global(5) = max(npts_global(5),nzs*nnxp(ng)*nnyp(ng)*npatch)
         !----- 3-D (nxp,nyp,npatch) ------------------------------------------------------!
         npts_local (6) = max(npts_local(6),nmxp(ng)*nmyp(ng)*npatch)
         npts_global(6) = max(npts_global(6),nnxp(ng)*nnyp(ng)*npatch)
         !----- 3-D (nxp,nyp,nwave) -------------------------------------------------------!
         npts_local (7) = max(npts_local(7),nmxp(ng)*nmyp(ng)*nwave)
         npts_global(7) = max(npts_global(7),nnxp(ng)*nnyp(ng)*nwave)
         !----- 4-D (nzp,nxp,nyp,nclouds) -------------------------------------------------!
         npts_local (8) = max(npts_local(8),nmzp(ng)*nmxp(ng)*nmyp(ng)*nclouds)
         npts_global(8) = max(npts_global(8),nnzp(ng)*nnxp(ng)*nnyp(ng)*nclouds)
         !----- 3-D (nxp,nyp,nclouds) -----------------------------------------------------!
         npts_local (9) = max(npts_local(9),nmxp(ng)*nmyp(ng)*nclouds)
         npts_global(9) = max(npts_global(9),nnxp(ng)*nnyp(ng)*nclouds)
      end do
      
      !----- Set the maximum possible number for scr1 and scr2. ---------------------------!
      nptsmax = maxval(npts_global)
      mpts1d  = maxval(npts_local,mask = number_dims == 1)
      mpts2d  = maxval(npts_local,mask = number_dims == 2)
      mpts3d  = maxval(npts_local,mask = number_dims == 3)
      mpts4d  = maxval(npts_local,mask = number_dims == 4)
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !    Allocate arrays based on options (if necessary).                                !
      !------------------------------------------------------------------------------------!
      allocate (scratch%scr1 (nptsmax))
      allocate (scratch%scr2 (nptsmax))
      allocate (scratch%scr3 (nptsmax))
      allocate (scratch%scr4 (nptsmax))
      allocate (scratch%scr5 (nptsmax))
      allocate (scratch%scr6 (nptsmax))

      allocate (scratch%vt2da(mpts2d))
      allocate (scratch%vt2db(mpts2d))
      allocate (scratch%vt2dc(mpts2d))
      allocate (scratch%vt2dd(mpts2d))
      allocate (scratch%vt2de(mpts2d))
      allocate (scratch%vt2df(mpts2d))
      allocate (scratch%vt2dg(mpts2d))
      allocate (scratch%vt2dh(mpts2d))
      allocate (scratch%vt2di(mpts2d))
      allocate (scratch%vt2dj(mpts2d))
      allocate (scratch%vt2dk(mpts2d))
      allocate (scratch%vt2dl(mpts2d))
      allocate (scratch%vt2dm(mpts2d))
      allocate (scratch%vt2dn(mpts2d))
      allocate (scratch%vt2do(mpts2d))
      allocate (scratch%vt2dp(mpts2d))
      allocate (scratch%vt2dq(mpts2d))
      allocate (scratch%vt2dr(mpts2d))
      allocate (scratch%vt2ds(mpts2d))

      allocate (scratch%vt3da(mpts3d))
      allocate (scratch%vt3db(mpts3d))
      allocate (scratch%vt3dc(mpts3d))
      allocate (scratch%vt3dd(mpts3d))
      allocate (scratch%vt3de(mpts3d))
      allocate (scratch%vt3df(mpts3d))
      allocate (scratch%vt3dg(mpts3d))
      allocate (scratch%vt3dh(mpts3d))
      allocate (scratch%vt3di(mpts3d))
      allocate (scratch%vt3dj(mpts3d))
      allocate (scratch%vt3dk(mpts3d))
      allocate (scratch%vt3dl(mpts3d))
      allocate (scratch%vt3dm(mpts3d))
      allocate (scratch%vt3dn(mpts3d))
      allocate (scratch%vt3do(mpts3d))
      allocate (scratch%vt3dp(mpts3d))
      allocate (scratch%vt3dq(mpts3d))
      allocate (scratch%vt3dr(mpts3d))
      allocate (scratch%vt3ds(mpts3d))


      allocate (scratch%vt4da(mpts4d))
      allocate (scratch%vt4db(mpts4d))
      allocate (scratch%vt4dc(mpts4d))

      !----- ALF - Put zero in every variable. --------------------------------------------!
      call azero(nptsmax, scratch%scr1)
      call azero(nptsmax, scratch%scr2)
      call azero(nptsmax, scratch%scr3)
      call azero(nptsmax, scratch%scr4)
      call azero(nptsmax, scratch%scr5)
      call azero(nptsmax, scratch%scr6)

      call azero(mpts2d , scratch%vt2da)
      call azero(mpts2d , scratch%vt2db)
      call azero(mpts2d , scratch%vt2dc)
      call azero(mpts2d , scratch%vt2dd)
      call azero(mpts2d , scratch%vt2de)
      call azero(mpts2d , scratch%vt2df)
      call azero(mpts2d , scratch%vt2dg)
      call azero(mpts2d , scratch%vt2dh)
      call azero(mpts2d , scratch%vt2di)
      call azero(mpts2d , scratch%vt2dj)
      call azero(mpts2d , scratch%vt2dk)
      call azero(mpts2d , scratch%vt2dl)
      call azero(mpts2d , scratch%vt2dm)
      call azero(mpts2d , scratch%vt2dn)
      call azero(mpts2d , scratch%vt2do)
      call azero(mpts2d , scratch%vt2dp)
      call azero(mpts2d , scratch%vt2dq)
      call azero(mpts2d , scratch%vt2dr)
      call azero(mpts2d , scratch%vt2ds)

      call azero(mpts3d , scratch%vt3da)
      call azero(mpts3d , scratch%vt3db)
      call azero(mpts3d , scratch%vt3dc)
      call azero(mpts3d , scratch%vt3dd)
      call azero(mpts3d , scratch%vt3de)
      call azero(mpts3d , scratch%vt3df)
      call azero(mpts3d , scratch%vt3dg)
      call azero(mpts3d , scratch%vt3dh)
      call azero(mpts3d , scratch%vt3di)
      call azero(mpts3d , scratch%vt3dj)
      call azero(mpts3d , scratch%vt3dk)
      call azero(mpts3d , scratch%vt3dl)
      call azero(mpts3d , scratch%vt3dm)
      call azero(mpts3d , scratch%vt3dn)
      call azero(mpts3d , scratch%vt3do)
      call azero(mpts3d , scratch%vt3dp)
      call azero(mpts3d , scratch%vt3dq)
      call azero(mpts3d , scratch%vt3dr)
      call azero(mpts3d , scratch%vt3ds)

      call azero(mpts4d , scratch%vt4da)
      call azero(mpts4d , scratch%vt4db)
      call azero(mpts4d , scratch%vt4dc)

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
      if (associated(scratch%scr4 ))  nullify (scratch%scr4 )
      if (associated(scratch%scr5 ))  nullify (scratch%scr5 )
      if (associated(scratch%scr6 ))  nullify (scratch%scr6 )

      if (associated(scratch%vt2da))  nullify (scratch%vt2da)
      if (associated(scratch%vt2db))  nullify (scratch%vt2db)
      if (associated(scratch%vt2dc))  nullify (scratch%vt2dc)
      if (associated(scratch%vt2dd))  nullify (scratch%vt2dd)
      if (associated(scratch%vt2de))  nullify (scratch%vt2de)
      if (associated(scratch%vt2df))  nullify (scratch%vt2df)
      if (associated(scratch%vt2dg))  nullify (scratch%vt2dg)
      if (associated(scratch%vt2dh))  nullify (scratch%vt2dh)
      if (associated(scratch%vt2di))  nullify (scratch%vt2di)
      if (associated(scratch%vt2dj))  nullify (scratch%vt2dj)
      if (associated(scratch%vt2dk))  nullify (scratch%vt2dk)
      if (associated(scratch%vt2dl))  nullify (scratch%vt2dl)
      if (associated(scratch%vt2dm))  nullify (scratch%vt2dm)
      if (associated(scratch%vt2dn))  nullify (scratch%vt2dn)
      if (associated(scratch%vt2do))  nullify (scratch%vt2do)
      if (associated(scratch%vt2dp))  nullify (scratch%vt2dp)
      if (associated(scratch%vt2dq))  nullify (scratch%vt2dq)
      if (associated(scratch%vt2dr))  nullify (scratch%vt2dr)
      if (associated(scratch%vt2ds))  nullify (scratch%vt2ds)

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
      if (associated(scratch%scr4 ))  deallocate (scratch%scr4 )
      if (associated(scratch%scr5 ))  deallocate (scratch%scr5 )
      if (associated(scratch%scr6 ))  deallocate (scratch%scr6 )

      if (associated(scratch%vt2da))  deallocate (scratch%vt2da)
      if (associated(scratch%vt2db))  deallocate (scratch%vt2db)
      if (associated(scratch%vt2dc))  deallocate (scratch%vt2dc)
      if (associated(scratch%vt2dd))  deallocate (scratch%vt2dd)
      if (associated(scratch%vt2de))  deallocate (scratch%vt2de)
      if (associated(scratch%vt2df))  deallocate (scratch%vt2df)
      if (associated(scratch%vt2dg))  deallocate (scratch%vt2dg)
      if (associated(scratch%vt2dh))  deallocate (scratch%vt2dh)
      if (associated(scratch%vt2di))  deallocate (scratch%vt2di)
      if (associated(scratch%vt2dj))  deallocate (scratch%vt2dj)
      if (associated(scratch%vt2dk))  deallocate (scratch%vt2dk)
      if (associated(scratch%vt2dl))  deallocate (scratch%vt2dl)
      if (associated(scratch%vt2dm))  deallocate (scratch%vt2dm)
      if (associated(scratch%vt2dn))  deallocate (scratch%vt2dn)
      if (associated(scratch%vt2do))  deallocate (scratch%vt2do)
      if (associated(scratch%vt2dp))  deallocate (scratch%vt2dp)
      if (associated(scratch%vt2dq))  deallocate (scratch%vt2dq)
      if (associated(scratch%vt2dr))  deallocate (scratch%vt2dr)
      if (associated(scratch%vt2ds))  deallocate (scratch%vt2ds)

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
