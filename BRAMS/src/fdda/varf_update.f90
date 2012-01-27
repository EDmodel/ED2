!===================================== Change Log =========================================!
! 5.0.0                                                                                    !
!                                                                                          !
!==========================================================================================!
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved                       !
!  Regional Atmospheric Modeling System - RAMS                                             !
!==========================================================================================!
subroutine varf_update(iswap,initflag,ifileok)
   use mem_leaf
   use mem_varinit
   use mem_basic
   use mem_grid
   use mem_scratch
   use therm_lib  , only : level
   use rconstants , only : toodry
   use grid_dims  , only : str_len
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer               , intent(in)  :: initflag
   integer               , intent(in)  :: iswap
   logical               , intent(out) :: ifileok
   !----- Local variables. ----------------------------------------------------------------!
   character(len=str_len)              :: flnm
   character(len=7)                    :: cgrid
   integer                             :: iver_var
   integer                             :: nc
   integer                             :: iyearx
   integer                             :: imonthx
   integer                             :: idatex
   integer                             :: ihourx
   integer                             :: nxpx
   integer                             :: nypx
   integer                             :: nzpx
   integer                             :: imarker
   integer                             :: i
   integer                             :: j
   integer                             :: k
   real                                :: rlatx
   real                                :: wlon1x
   real                                :: deltaxx
   real                                :: deltayx
   real                                :: deltazx
   real                                :: dzratx
   real                                :: dzmaxx
   !----- Local constants. ----------------------------------------------------------------!
   integer               , parameter   :: iun = 22
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !    Check and see what we are doing. If it is initial time, read fields into regular   !
   ! arrays.  If not, see if nudging will be done on this grid if it is a nested grid.     !
   !---------------------------------------------------------------------------------------!
   if (ngrid > 1 .and. tnudcent+tnudtop < .001 .and. initflag == 0) return


   !---------------------------------------------------------------------------------------!
   !     Put new fields into varinit future arrays. If iswap is 1, then swap future into   !
   ! past first.                                                                           !
   !---------------------------------------------------------------------------------------!
   if (iswap == 1) then
      call atob(nxyzp,varinit_g(ngrid)%varuf,varinit_g(ngrid)%varup)
      call atob(nxyzp,varinit_g(ngrid)%varvf,varinit_g(ngrid)%varvp)
      call atob(nxyzp,varinit_g(ngrid)%varpf,varinit_g(ngrid)%varpp)
      call atob(nxyzp,varinit_g(ngrid)%vartf,varinit_g(ngrid)%vartp)
      call atob(nxyzp,varinit_g(ngrid)%varrf,varinit_g(ngrid)%varrp)
      if (co2_on) then
         call atob(nxyzp,varinit_g(ngrid)%varof,varinit_g(ngrid)%varop)
      end if
   end if

   !----- Make data file name from tag file name. -----------------------------------------!
   write(cgrid,fmt='(a2,i1,a4)') '-g',ngrid,'.vfm'
   nc=len_trim(fnames_varf(nvarffl))
   flnm=fnames_varf(nvarffl)(1:nc-4)//trim(cgrid)

   !----- Check for existence. ------------------------------------------------------------!
   inquire(file=trim(flnm),exist=ifileok)

   !----- It must have at least one grid... -----------------------------------------------!
   if (.not. ifileok .and. ngrid == 1) then
      call abort_run ('Missing grid 1 varfile :'//trim(flnm)//'...'                        &
                     ,'varf_update','varf_update.f90')
   elseif (.not. ifileok) then
      !----- File is not here, but it has a parent grid, so nothing to be done here... ----!
      return
   end if


   !---------------------------------------------------------------------------------------!
   !      Read the varfile fields into the "future" varinit arrays.  These will be swapped !
   ! to the past arrays when needed.                                                       !
   !---------------------------------------------------------------------------------------!
   call rams_f_open(iun,flnm,'FORMATTED','OLD','READ',0)

   !----- Find varfile "version". ---------------------------------------------------------!
   read (unit=iun,fmt=*) imarker
   rewind (unit=iun)

   if(imarker == 999999) then
      read(unit=iun,fmt=*) imarker,iver_var
   else
      iver_var=1
   end if

   read(unit=iun,fmt=*) iyearx,imonthx,idatex,ihourx,nxpx,nypx,nzpx,rlatx,wlon1x           &
                       ,deltaxx,deltayx,deltazx,dzratx,dzmaxx

   if(nxp /= nxpx .or. nyp /= nypx .or. nzp /= nzpx .or. abs(deltax-deltaxx) > .001 .or.   &
      abs(deltay-deltayx) > .001 .or. abs(deltaz-deltazx) > .001 .or.                      &
      abs(dzrat-dzratx) > .001 .or. abs(dzmax-dzmaxx) > .001 .or.                          &
      abs(platn(ngrid)-rlatx) > .001 .or. abs(plonn(ngrid)-wlon1x) > .001) then
      
      write(unit=*,fmt='(a)')            '================================================'
      write(unit=*,fmt='(a)')            '!! Grid mismatch between varfile and namelist !!'
      write(unit=*,fmt='(a,1x,a)')       ' File:',trim(flnm)
      write(unit=*,fmt='(a,1x,i5,1x,a)') ' Check values for grid:',ngrid,'(File,Namelist).'
      write(unit=*,fmt='(a,2(1x,i5))')     ' NXP     :',nxpx,nxp
      write(unit=*,fmt='(a,2(1x,i5))')     ' NYP     :',nypx,nyp
      write(unit=*,fmt='(a,2(1x,i5))')     ' NZP     :',nzpx,nzp
      write(unit=*,fmt='(a,2(1x,es12.5))') ' DELTAX  :',deltaxx,deltax
      write(unit=*,fmt='(a,2(1x,es12.5))') ' DELTAY  :',deltayx,deltay
      write(unit=*,fmt='(a,2(1x,es12.5))') ' DELTAZ  :',deltazx,deltaz
      write(unit=*,fmt='(a,2(1x,es12.5))') ' DZRAT   :',dzratx ,dzrat
      write(unit=*,fmt='(a,2(1x,es12.5))') ' DZMAX   :',dzmaxx ,dzmax
      write(unit=*,fmt='(a,2(1x,es12.5))') ' POLELAT :',rlatx  ,platn(ngrid)
      write(unit=*,fmt='(a,2(1x,es12.5))') ' POLELON :',wlon1x ,plonn(ngrid)
      write(unit=*,fmt='(a)')            '================================================'

      call abort_run('Grid mismatch between varfile and namelist...'                       &
                    ,'varf_update','varf_update.f90')
   end if

   call vfirec(iun,varinit_g(ngrid)%varuf,nxyzp,'LIN')
   call vfirec(iun,varinit_g(ngrid)%varvf,nxyzp,'LIN')
   call vfirec(iun,varinit_g(ngrid)%varpf,nxyzp,'LIN')
   call vfirec(iun,varinit_g(ngrid)%vartf,nxyzp,'LIN')
   call vfirec(iun,varinit_g(ngrid)%varrf,nxyzp,'LIN')

   !---------------------------------------------------------------------------------------! 
   !     Mixing ratio should never be zero or negative, making it at least "toodry"...     !
   !---------------------------------------------------------------------------------------!
   where (varinit_g(ngrid)%varrf(:,:,:) < toodry)
      varinit_g(ngrid)%varrf(:,:,:) = toodry
   end where
   
   !---------------------------------------------------------------------------------------!
   !     Decide what to do with CO2 based on what the user asked for.  In any case, we     !
   ! copy CO2 to scratch array vt3do so it will work even when CO2 is not used.            !
   !---------------------------------------------------------------------------------------!
   select case (ico2)
   case (0) !----- Skip, nothing needs to be done in this case. ---------------------------!
      call ae0(nzp*nxp*nyp,scratch%vt3do,co2con(1))
   case (1) !----- Initialise with constant value. ----------------------------------------!
      call ae0(nzp*nxp*nyp,varinit_g(ngrid)%varof,co2con(1))
      call atob(nzp*nxp*nyp,varinit_g(ngrid)%varof,scratch%vt3do)
   case (2) !----- Initialise with constant profile. --------------------------------------!
      call a3e1(nzp,nxp,nyp,1,nxp,1,nyp,1,nzp,varinit_g(ngrid)%varof,co2con(1:nzp))
      call atob(nzp*nxp*nyp,varinit_g(ngrid)%varof,scratch%vt3do)
   case (3) !----- Initialise with data read from ISAN files. -----------------------------!
      call vfirec(iun,varinit_g(ngrid)%varof,nxyzp,'LIN')
      call atob(nxyzp,varinit_g(ngrid)%varof,scratch%vt3do)
   end select
   !---------------------------------------------------------------------------------------!

            
   if(initflag == 1 .and. iver_var == 2) then
      ! Extract snow depth from the varfile. Ignore other 2D fields for now.
      call vfirec(iun,scratch%vt2da,nxyp,'LIN')
      call vfirec(iun,scratch%vt2da,nxyp,'LIN')
      call vfirec(iun,scratch%vt2da,nxyp,'LIN')
      call vfirec(iun,leaf_g(ngrid)%snow_mass,nxyp,'LIN')
      call vfirec(iun,scratch%vt2da,nxyp,'LIN')
   endif

   close(unit=iun, status='keep')

   !----- If running ADAP coord, do interpolation to Cartesian levels. --------------------!

   if (if_adap == 1) then
      call varf_adap(nzp,nxp,nyp,varinit_g(ngrid)%varuf            &
                    ,varinit_g(ngrid)%varvf,varinit_g(ngrid)%varpf,varinit_g(ngrid)%vartf  &
                    ,varinit_g(ngrid)%varrf,scratch%vt3do,grid_g(ngrid)%topta)
      if (co2_on) then
         call atob(nxyzp,scratch%vt3do,varinit_g(ngrid)%varof)
      end if
   end if

   !----- Find the reference state. -------------------------------------------------------!

   if (initflag == 1 .and. ngrid == 1)  then
        call varref(nzp, nxp, nyp, varinit_g(ngrid)%vartf     , varinit_g(ngrid)%varpf     &
                                 , basic_g(ngrid)%pi0         , basic_g(ngrid)%th0         &
                                 , varinit_g(ngrid)%varrf     , scratch%vt3do              &
                                 , basic_g(ngrid)%dn0         , basic_g(ngrid)%dn0u        &
                                 , basic_g(ngrid)%dn0v        , varinit_g(ngrid)%varuf     &
                                 , varinit_g(ngrid)%varvf     , grid_g(ngrid)%topt         &
                                 , grid_g(ngrid)%topu         , grid_g(ngrid)%topv         &
                                 , grid_g(ngrid)%rtgt         , grid_g(ngrid)%rtgu         &
                                 , grid_g(ngrid)%rtgv         , grid_g(ngrid)%topta )
   end if

   !----- We don't really deal with the Exner function, but the perturbation. -------------!
   do j=1,nyp
      do i=1,nxp
         do k=1,nzp
            varinit_g(ngrid)%varpf(k,i,j) = varinit_g(ngrid)%varpf(k,i,j)                  &
                                          - basic_g(ngrid)%pi0(k,i,j)
         end do
      end do
   end do

   !----- If this is an initialization, put data into regular arrays. ---------------------!
   if (initflag == 1 ) then
      call atob(nxyzp,varinit_g(ngrid)%varuf,basic_g(ngrid)%uc  )
      call atob(nxyzp,varinit_g(ngrid)%varvf,basic_g(ngrid)%vc  )
      call atob(nxyzp,varinit_g(ngrid)%varpf,basic_g(ngrid)%pc  )
      call atob(nxyzp,varinit_g(ngrid)%vartf,basic_g(ngrid)%thp )
      call atob(nxyzp,varinit_g(ngrid)%varrf,basic_g(ngrid)%rtp )
      if (co2_on) call atob(nxyzp,varinit_g(ngrid)%varof,basic_g(ngrid)%co2p )
   end if

   return
end subroutine varf_update
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!      This subroutine will interpolate from sigma-z varfile vertical coords to ADAP grid. !
!------------------------------------------------------------------------------------------!
subroutine varf_adap(n1,n2,n3,varu,varv,varp,vart,varr,varc,topta)
   use mem_scratch
   use mem_grid
   use rconstants
   use therm_lib, only: virtt
   implicit none

   !----- Arguments. ----------------------------------------------------------------------!
   integer                  , intent(in)    :: n1
   integer                  , intent(in)    :: n2
   integer                  , intent(in)    :: n3
   !----- Local variables. ----------------------------------------------------------------!
   real, dimension(n1,n2,n3), intent(inout) :: varu
   real, dimension(n1,n2,n3), intent(inout) :: varv
   real, dimension(n1,n2,n3), intent(inout) :: varp
   real, dimension(n1,n2,n3), intent(inout) :: vart
   real, dimension(n1,n2,n3), intent(inout) :: varr
   real, dimension(n1,n2,n3), intent(inout) :: varc
   real, dimension(n2,n3)   , intent(in)    :: topta
   !----- Local constants. ----------------------------------------------------------------!
   integer                                  :: i
   integer                                  :: j
   integer                                  :: k
   !---------------------------------------------------------------------------------------!
   

   do j=1,n3
      do i=1,n2
      
         do k=1,n1
            vctr10(k) = topta(i,j) + (1. - topta(i,j)/ztop) * ztn(k,ngrid)
         end do

         !----- Copy fields to temporary arrays. ------------------------------------------!
         vctr1(1:n1)=varu(1:n1,i,j)
         vctr2(1:n1)=varv(1:n1,i,j)
         vctr3(1:n1)=vart(1:n1,i,j)
         vctr4(1:n1)=varr(1:n1,i,j)
         vctr5(1:n1)=varc(1:n1,i,j)

         call htint2(n1,vctr1,vctr10,n1,vctr11,ztn(:,ngrid))
         call htint2(n1,vctr2,vctr10,n1,vctr12,ztn(:,ngrid))
         call htint2(n1,vctr3,vctr10,n1,vctr13,ztn(:,ngrid))
         call htint2(n1,vctr4,vctr10,n1,vctr14,ztn(:,ngrid))
         call htint2(n1,vctr5,vctr10,n1,vctr15,ztn(:,ngrid))

         !----- Do a hydrostatic balance. -------------------------------------------------!
         do k=1,n1
            vctr27(k) = virtt(vctr13(k),vctr14(k))
         end do

         vctr28(n1)= varp(n1,i,j) + grav * (ztn(n1,ngrid) - vctr10(n1)) / vctr27(n1)
         do k = n1-1,1,-1
            vctr28(k) = vctr28(k+1) + grav * (ztn(k+1,ngrid)-ztn(k,ngrid))                 &
                                    / ((vctr27(k)+vctr27(k+1))*.5)
         end do
         
         !----- Copying back the variables. -----------------------------------------------!
         varu(1:n1,i,j) = vctr11(1:n1)
         varv(1:n1,i,j) = vctr12(1:n1)
         vart(1:n1,i,j) = vctr13(1:n1)
         varr(1:n1,i,j) = vctr14(1:n1)
         varc(1:n1,i,j) = vctr15(1:n1)
         varp(1:n1,i,j) = vctr28(1:n1)
         
      end do
   end do

   return
end subroutine varf_adap
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine will compute the reference sounding.                                 !
!------------------------------------------------------------------------------------------!
subroutine varref(n1,n2,n3,thp,pc,pi0,th0,rtp,co2p,dn0,dn0u,dn0v,uc,vc,topt,topu,topv,rtgt &
                 ,rtgu,rtgv,topta)

   use mem_grid
   use ref_sounding
   use mem_scratch
   use rconstants
   use therm_lib   , only : virtt     & ! intent(in)
                          , vapour_on ! ! intent(in)
   use mem_basic   , only : co2_on    & ! intent(in)
                          , co2con    ! ! intent(in)

   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                  , intent(in)    :: n1
   integer                  , intent(in)    :: n2
   integer                  , intent(in)    :: n3
   !----- Local variables. ----------------------------------------------------------------!
   real, dimension(n1,n2,n3), intent(in)    :: thp
   real, dimension(n1,n2,n3), intent(in)    :: pc
   real, dimension(n1,n2,n3), intent(in)    :: rtp
   real, dimension(n1,n2,n3), intent(in)    :: co2p
   real, dimension(n1,n2,n3), intent(in)    :: uc
   real, dimension(n1,n2,n3), intent(in)    :: vc
   real, dimension(n2,n3)   , intent(in)    :: topt
   real, dimension(n2,n3)   , intent(in)    :: topu
   real, dimension(n2,n3)   , intent(in)    :: topv
   real, dimension(n2,n3)   , intent(in)    :: rtgt
   real, dimension(n2,n3)   , intent(in)    :: rtgu
   real, dimension(n2,n3)   , intent(in)    :: rtgv
   real, dimension(n2,n3)   , intent(in)    :: topta
   real, dimension(n1,n2,n3), intent(inout) :: pi0
   real, dimension(n1,n2,n3), intent(inout) :: th0
   real, dimension(n1,n2,n3), intent(inout) :: dn0
   real, dimension(n1,n2,n3), intent(inout) :: dn0u
   real, dimension(n1,n2,n3), intent(inout) :: dn0v
   !----- Local constants. ----------------------------------------------------------------!
   integer                                  :: i
   integer                                  :: j
   integer                                  :: k
   !---------------------------------------------------------------------------------------!

   !----- Reference sounding is point with lowest topography. -----------------------------!
   topref = 1.e10
   do j=1,nyp
      do i=1,nxp
         if (topta(i,j) < topref) then
            iref   = i
            jref   = j
            topref = topta(i,j)
         end if
      end do
   end do

   !-----  Set up 1-D reference state, depending on the coordinate. -----------------------!
   select case (if_adap)
   case (0) !----- Sigma-z, terrain following coordinate. ---------------------------------!
      do k=1,nzp
         vctr2(k) = ztn(k,ngrid) * (1. - topref/ztop) + topref
      end do
      call htint2(nzp,thp(:,iref,jref),vctr2,nzp,vctr1         ,zt)
      call htint2(nzp,uc (:,iref,jref),vctr2,nzp,u01dn(:,ngrid),zt)
      call htint2(nzp,vc (:,iref,jref),vctr2,nzp,v01dn(:,ngrid),zt)
      if (vapour_on) then
         call htint2(nzp,rtp(:,iref,jref),vctr2,nzp,rt01dn(:,ngrid),zt)
      else
         rt01dn(1:nzp,ngrid) = 0.
      end if
      if (co2_on) then
         call htint2(nzp,co2p(:,iref,jref),vctr2,nzp,co201dn(:,ngrid),zt)
      else
         co201dn(1:nzp,ngrid) = co2con(1)
      end if
   case (1)
      vctr2(1:nzp)       = ztn(1:nzp,ngrid)
      vctr1(1:nzp)       = thp(1:nzp,iref,jref)
      u01dn(1:nzp,ngrid) = uc(1:nzp,iref,jref)
      v01dn(1:nzp,ngrid) = vc(1:nzp,iref,jref)
      if (vapour_on) then
         rt01dn(1:nzp,ngrid) = rtp(1:nzp,iref,jref)
      else
         rt01dn(1:nzp,ngrid) = 0.
      end if
      if (co2_on) then
         co201dn(1:nzp,ngrid) = co2p(1:nzp,iref,jref)
      else
         co201dn(1:nzp,ngrid) = co2con(1)
      end if
   end select

   !------ We now compute the reference theta, which will indeed be theta_v. --------------!
   do k = 1,nzp
      th01dn(k,ngrid) = virtt(vctr1(k),rt01dn(k,ngrid))
   end do

   !----- Boundary condition. -------------------------------------------------------------!
   u01dn(1,ngrid)   = u01dn(2,ngrid)
   v01dn(1,ngrid)   = v01dn(2,ngrid)
   rt01dn(1,ngrid)  = rt01dn(2,ngrid)
   th01dn(1,ngrid)  = th01dn(2,ngrid)
   co201dn(1,ngrid) = co201dn(1,ngrid)

   !----- Hydrostatic adjustment for pressure. --------------------------------------------!
   pi01dn(1,ngrid) = pc(1,iref,jref) + grav * (vctr2(1) - zt(1))                           &
                                     / (.5 * ( th01dn(1,ngrid)                             &
                                             + virtt(thp(1,iref,jref),rtp(1,iref,jref)) ))
   do k = 2,nzp
     pi01dn(k,ngrid) = pi01dn(k-1,ngrid)                                                   &
                     - grav / (dzm(k-1) * .5 * (th01dn(k,ngrid) + th01dn(k-1,ngrid)))
   end do

   !----- Finding the density. ------------------------------------------------------------!
   do k = 1,nzp
     vctr4(k)        = (pi01dn(k,ngrid) / cp) ** cpor * p00
     dn01dn(k,ngrid) = cp * vctr4(k) / (rdry * th01dn(k,ngrid) * pi01dn(k,ngrid))
   end do

   !------ Compute 3-D reference state from 1-D reference state. --------------------------!
   call refs3d(nzp,nxp,nyp,pi0,dn0,dn0u,dn0v,th0,topt,rtgt)

   return
end subroutine varref
!==========================================================================================!
!==========================================================================================!
