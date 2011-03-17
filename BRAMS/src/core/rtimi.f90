!============================= Change Log =================================================!
! 5.0.0                                                                                    !
!                                                                                          !
!==========================================================================================!
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved                       !
!  Regional Atmospheric Modeling System - RAMS                                             !
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!      This routine simply sets all tendency arrays to zero.                               !
!------------------------------------------------------------------------------------------!
subroutine tend0()

   use mem_grid
   use mem_tend
   use var_tables
   use node_mod

   implicit none
   !----- Local variables. ----------------------------------------------------------------!
   integer     :: n
   integer     :: mxyzp
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     First we flush the "acoustic variables" tendencies to zero.                       !
   !---------------------------------------------------------------------------------------!
   mxyzp = mxp * myp * mzp
   call azero(mxyzp,tend_g(ngrid)%ut)
   call azero(mxyzp,tend_g(ngrid)%vt)
   call azero(mxyzp,tend_g(ngrid)%wt)
   call azero(mxyzp,tend_g(ngrid)%pt)
   !---------------------------------------------------------------------------------------!



   !----- Now we flush the other scalar variables. ----------------------------------------!
   do n = 1,num_scalar(ngrid)
      call azero(mxyzp,scalar_tab(n,ngrid)%var_t)
   end do
   !---------------------------------------------------------------------------------------!

   return
end subroutine tend0
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This sub-routine applies the Asselin filter.  For the velocities and pressure, this  !
! must be done in two stages, the first when IAC=1 and the second when IAC=2.              !
!------------------------------------------------------------------------------------------!
subroutine hadvance(iac)
   use mem_grid
   use mem_tend
   use mem_basic
   use mem_scratch
   use node_mod
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer, intent(in) :: iac
   !----- Local variables. ----------------------------------------------------------------!
   integer :: mxyzp
   !---------------------------------------------------------------------------------------!

   mxyzp = mxp * myp * mzp
   eps   = .2

   !----- For both IAC=1 and IAC=2, call PREDICT for U, V, W, and P. ----------------------!
   call predict(mxyzp,basic_g(ngrid)%uc,basic_g(ngrid)%up,tend_g(ngrid)%ut,scratch%vt3da   &
               ,iac,dtlv)

   if (icorflg == 1 .or. jdim == 1) then
      call predict(mxyzp,basic_g(ngrid)%vc,basic_g(ngrid)%vp,tend_g(ngrid)%vt              &
                  ,scratch%vt3da,iac,dtlv)
   end if

   call predict(mxyzp,basic_g(ngrid)%wc,basic_g(ngrid)%wp,tend_g(ngrid)%wt,scratch%vt3da   &
               ,iac,dtlv)
   call predict(mxyzp,basic_g(ngrid)%pc,basic_g(ngrid)%pp,tend_g(ngrid)%pt,scratch%vt3da   &
               ,iac,dtlv)

   return
end subroutine hadvance
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     For IAC=3, this routine moves the arrays AC and AP forward by 1 time level by adding !
! in the prescribed tendency. It also applies the Asselin filter given by:                 !
!                                                                                          !
!              {AC} = AC + EPS * (AP - 2 * AC + AF)                                        !
!                                                                                          !
! where AP,AC,AF are the past, current and future time levels of A.  All IAC=1 does is to  !
! perform the {AC} calculation without the AF term present.  IAC=2 completes the           !
! calculation of {AC} by adding the AF term only, and advances AC by filling it with input !
! AP values which were already updated in ACOUSTC.                                         !
!------------------------------------------------------------------------------------------!
subroutine predict(npts,ac,ap,fa,af,iac,dtlp)

   use mem_grid
   use node_mod

   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                 , intent(in)    :: npts
   integer                 , intent(in)    :: iac
   real                    , intent(in)    :: dtlp
   real   , dimension(npts), intent(in)    :: fa
   real   , dimension(npts), intent(inout) :: ac
   real   , dimension(npts), intent(inout) :: ap
   real   , dimension(npts), intent(inout) :: af
   !----- Local variables. ----------------------------------------------------------------!
   integer                                 :: m
   real                                    :: epsu
   !---------------------------------------------------------------------------------------!



   !----- Set up the epsilon factor. ------------------------------------------------------!
   epsu = eps
   if (ngbegun(ngrid) == 0) epsu = 0.5
   !---------------------------------------------------------------------------------------!


   select case (iac)
   case (1)
      do m = 1,npts
         ac(m) = ac(m) + epsu * (ap(m) - 2. * ac(m))
      end do
      return

   case (2)
      do m = 1,npts
         af(m) = ap(m)
         ap(m) = ac(m) + epsu * af(m)
      end do

   case (3)
      do m = 1,npts
         af(m) = ap(m) + dtlp * fa(m)
      end do
      if (ngrid == 1 .and. ipara == 0) call cyclic_set(nzp,nxp,nyp,af,'T')
      do m = 1,npts
         ap(m) = ac(m) + epsu * (ap(m) - 2. * ac(m) + af(m))
      end do
   end select

   do m = 1,npts
     ac(m) = af(m)
   end do

   return
end subroutine predict
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!   This sub-routine updates the thermodynamic fields.                                     !
!------------------------------------------------------------------------------------------!
subroutine predtr()
   use mem_grid
   use var_tables
   use node_mod

   implicit none
   !----- Local variables. ----------------------------------------------------------------!
   integer :: mxyzp
   integer :: n
   !---------------------------------------------------------------------------------------!


   mxyzp = mxp * myp * mzp

   do n = 1,num_scalar(ngrid)
      call update(mxyzp,scalar_tab(n,ngrid)%var_p,scalar_tab(n,ngrid)%var_t, dtlt)
   end do

   return
end subroutine predtr
!==========================================================================================!
!==========================================================================================!
