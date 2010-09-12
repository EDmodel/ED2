!==========================================================================================!
!==========================================================================================!
!     This module contains the subroutines that drive CARMA radiation.                     !
!------------------------------------------------------------------------------------------!
module rad_carma
  
   use mem_carma

   !=======================================================================================!
   !=======================================================================================!


   contains



   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine is an interface between BRAMS and the actual CARMA driver.         !
   !---------------------------------------------------------------------------------------!
   subroutine radcomp_carma(m1,m2,m3,npat,nclouds,ncrad,ia,iz,ja,jz,mynum,iswrtyp,ilwrtyp  &
                           ,icumfdbk,solfac,theta,pi0,pp,rv,rain,lwl,iwl,cuprliq,cuprice   &
                           ,cuparea,cupierr,dn0,rtp,fthrd,rtgt,f13t,f23t,glat,glon,rshort  &
                           ,rshort_top,rshortup_top,rlong,rlongup_top,albedt,cosz,rlongup  &
                           ,fmapt,pm,patch_area)
      use catt_start  , only : catt               ! ! intent(in)
      use mem_grid    , only : ngrid              ! ! intent(in)
      use grid_dims   , only : nzpmax             ! ! intent(in)
      use mem_aerad   , only : nwave              ! ! intent(in)
      use rconstants  , only : day_sec            & ! intent(in)
                             , p00                & ! intent(in)
                             , t00                & ! intent(in)
                             , cpor               & ! intent(in)
                             , cpi                ! ! intent(in)
         
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      integer                          , intent(in)    :: m1
      integer                          , intent(in)    :: m2
      integer                          , intent(in)    :: m3
      integer                          , intent(in)    :: ia
      integer                          , intent(in)    :: iz
      integer                          , intent(in)    :: ja
      integer                          , intent(in)    :: jz
      integer                          , intent(in)    :: mynum
      integer                          , intent(in)    :: npat
      integer                          , intent(in)    :: nclouds
      integer                          , intent(inout) :: ncrad
      integer                          , intent(in)    :: iswrtyp
      integer                          , intent(in)    :: ilwrtyp
      integer                          , intent(in)    :: icumfdbk
      real                             , intent(in)    :: solfac
      real, dimension(m1,m2,m3)        , intent(in)    :: theta
      real, dimension(m1,m2,m3)        , intent(in)    :: pi0
      real, dimension(m1,m2,m3)        , intent(in)    :: pp
      real, dimension(m1,m2,m3)        , intent(in)    :: rv
      real, dimension(m1,m2,m3)        , intent(inout) :: lwl
      real, dimension(m1,m2,m3)        , intent(inout) :: iwl
      real, dimension(m1,m2,m3,nclouds), intent(in)    :: cuprliq
      real, dimension(m1,m2,m3,nclouds), intent(in)    :: cuprice
      real, dimension(   m2,m3,nclouds), intent(in)    :: cuparea
      real, dimension(   m2,m3,nclouds), intent(in)    :: cupierr
      real, dimension(   m2,m3)        , intent(in)    :: rain
      real, dimension(   m2,m3,npat)   , intent(in)    :: patch_area
      real, dimension(m1,m2,m3)        , intent(in)    :: dn0
      real, dimension(m1,m2,m3)        , intent(in)    :: rtp
      real, dimension(m2,m3)           , intent(in)    :: rtgt
      real, dimension(m2,m3)           , intent(in)    :: f13t
      real, dimension(m2,m3)           , intent(in)    :: f23t
      real, dimension(m2,m3)           , intent(in)    :: glat
      real, dimension(m2,m3)           , intent(in)    :: glon
      real, dimension(m2,m3)           , intent(in)    :: cosz
      real, dimension(m2,m3)           , intent(in)    :: albedt
      real, dimension(m2,m3)           , intent(in)    :: fmapt
      real, dimension(m1,m2,m3)        , intent(in)    :: pm
      real, dimension(m2,m3)           , intent(inout) :: rshort
      real, dimension(m2,m3)           , intent(inout) :: rshort_top
      real, dimension(m2,m3)           , intent(inout) :: rshortup_top
      real, dimension(m2,m3)           , intent(inout) :: rlong
      real, dimension(m2,m3)           , intent(inout) :: rlongup_top
      real, dimension(m1,m2,m3)        , intent(inout) :: fthrd
      real, dimension(m2,m3)           , intent(in)    :: rlongup
      !----- Local variables. -------------------------------------------------------------!
      integer                                          :: i
      integer                                          :: j
      integer                                          :: k
      integer                                          :: ipat
      integer                                          :: icld
      real                                             :: xland
      real, dimension(m1)                              :: lwl_cld
      real, dimension(m1)                              :: iwl_cld
      real                                             :: area_csky
      real                                             :: area_cld
      real                                             :: rain_eff
      real                                             :: rshort_cld
      real                                             :: rshort_top_cld
      real                                             :: rshortup_top_cld
      real                                             :: rlong_cld
      real                                             :: rlongup_top_cld
      real, dimension(m1)                              :: fthrd_cld
      real, dimension(nwave)                           :: aotl_cld
      real, dimension(nwave)                           :: aotl
      real                                             :: press
      real                                             :: tempc
      real                                             :: rvcgs
      !----- Locally saved variables. -----------------------------------------------------!
      logical                          , save          :: first_time = .true.
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !     If this is the first time we are calling CARMA, we must define some flags      !
      ! and allocate scratch variables.                                                    !
      !------------------------------------------------------------------------------------!
      if (first_time) then
         if (icumfdbk == 1) then
            ncrad = nclouds
         else
            ncrad = 0
         end if
         first_time = .false.
      end if

      !------------------------------------------------------------------------------------!
      !     Here we initialise some variables.                                             !
      !------------------------------------------------------------------------------------!
      call init_carma(m1)
      call setupbins

      jloop: do j=ja,jz
         iloop: do i=ia,iz
            
            !------------------------------------------------------------------------------!
            !      These variables turn the LW and SW radiation on or off.  The default is !
            ! to call longwave if ILWRTYP is 4, and shortwave is turned off unless ISWRTYP !
            ! is 4 and the sun is sufficiently above the horizon (no twilight).            !
            !------------------------------------------------------------------------------!
            ir_aerad  = ilwrtyp == 4
            isl_aerad = iswrtyp == 4 .and. cosz(i,j) > .03

            !----- Finding the amount of land. --------------------------------------------!
            xland = 1. - patch_area(i,j,1)

            !----- Precipitation is used only if the full cumulus feedback is off. --------!
            if (icumfdbk == 0) then
               rain_eff = rain(i,j)
            else
               rain_eff = 0.
            end if

            !----- Make sure that the hydrometeor mixing ratio is non-negative. -----------!
            do k = 1, m1
               lwl(k,i,j) = max(0.,lwl(k,i,j))
               iwl(k,i,j) = max(0.,iwl(k,i,j))
            end do

            !----- Initialise the clear sky area and AOT. ---------------------------------!
            area_csky = 1.
            aotl(:) = 0.
            cldloop: do icld = ncrad,0,-1

               !---------------------------------------------------------------------------!
               !     Copy cumulus cloud (or environment) mixing ratio to the scratch       !
               ! variables.  The fluxes will be an average of the area covered by cumulus  !
               ! clouds and cumulus-free.                                                  !
               !---------------------------------------------------------------------------!
               if (icld == 0) then
                  area_cld = area_csky
                  do k = 1, m1
                     lwl_cld(k) = 0.
                     iwl_cld(k) = 0.
                  end do
               else
                  if (cupierr(i,j,icld) /= 0.) cycle cldloop

                  area_cld  = cuparea(i,j,icld)
                  area_csky = area_csky - area_cld

                  do k = 1, m1
                     lwl_cld(k) = max(0.,cuprliq(k,i,j,icld))
                     iwl_cld(k) = max(0.,cuprice(k,i,j,icld))
                  end do
               end if
               !---------------------------------------------------------------------------!

               !----- Reset the scratch variables. ----------------------------------------!
               call flush_carma()
               rshort_cld       = 0.
               rshort_top_cld   = 0.
               rshortup_top_cld = 0.
               rlong_cld        = 0.
               rlongup_top_cld  = 0.
               fthrd_cld (:)    = 0.

               !----- Copying the aerosol optical transmission to a scratch array. --------!
               call ci_3d_1d(m2,m3,nwave,carma(ngrid)%aot,aotl_cld,i,j)
               !---------------------------------------------------------------------------!

               !----- Calling the main CARMA driver. --------------------------------------!
               call radcarma(nzpmax,m1,solfac,theta(1:m1,i,j),pi0(1:m1,i,j),pp(1:m1,i,j)   &
                            ,rv(1:m1,i,j),rain_eff,lwl(1:m1,i,j),iwl(1:m1,i,j)             &
                            ,lwl_cld,iwl_cld,dn0(1:m1,i,j),rtp(1:m1,i,j),fthrd_cld         &
                            ,rtgt(i,j),f13t(i,j),f23t(i,j),glat(i,j),glon(i,j),rshort_cld  &
                            ,rshort_top_cld,rshortup_top_cld,rlong_cld,rlongup_top_cld     &
                            ,albedt(i,j),cosz(i,j),rlongup(i,j),mynum,fmapt(i,j)           &
                            ,pm(1:m1,i,j),aotl_cld,xland)

               !----- Integrating the fluxes and optical depth. ---------------------------!
               if (isl_aerad) then
                   rshort      (i,j) = rshort      (i,j) + rshort_cld       * area_cld
                   rshort_top  (i,j) = rshort_top  (i,j) + rshort_top_cld   * area_cld
                   rshortup_top(i,j) = rshortup_top(i,j) + rshortup_top_cld * area_cld
               end if
               if (ir_aerad ) then
                  rlong        (i,j) = rlong       (i,j) + rlong_cld        * area_cld
                  rlongup_top  (i,j) = rlongup_top (i,j) + rlongup_top_cld  * area_cld
               end if
               do k = 1, m1
                  fthrd(k,i,j) = fthrd(k,i,j) + fthrd_cld(k) * area_cld
               end do
               aotl(:) = aotl(:) + aotl_cld(:) * area_cld
               
               !---------------------------------------------------------------------------!
               !     Print the output if the radiation heating rate is screwy.             !
               !---------------------------------------------------------------------------!
               if (any(abs(fthrd_cld) > 300./day_sec)) then
                  do k=1,m1
                     fthrd_cld(k) = fthrd_cld(k) * day_sec
                     lwl(k,i,j)   = lwl(k,i,j)   * 1000.
                     iwl(k,i,j)   = iwl(k,i,j)   * 1000.
                     lwl_cld(k)   = lwl_cld(k)   * 1000.
                     iwl_cld(k)   = iwl_cld(k)   * 1000.
                  end do

                  write (unit=*,fmt='(123a)') ('=',k=1,123)
                  write (unit=*,fmt='(a)')     ' Radiative heating rate is screwy!!! '
                  write (unit=*,fmt='(123a)') ('=',k=1,123)
                  write (unit=*,fmt='(a,1x,i4)')     ' - I       = ',i
                  write (unit=*,fmt='(a,1x,i4)')     ' - J       = ',j
                  write (unit=*,fmt='(a,1x,i4)')     ' - ICLD    = ',icld
                  write (unit=*,fmt='(a,1x,es12.5)') ' - RSHORT  = ',rshort_cld
                  write (unit=*,fmt='(a,1x,es12.5)') ' - RLONG   = ',rlong_cld
                  write (unit=*,fmt='(a,1x,es12.5)') ' - RLONGUP = ',rlongup(i,j)
                  write (unit=*,fmt='(a,1x,es12.5)') ' - COSZ    = ',cosz(i,j)
                  write (unit=*,fmt='(a,1x,es12.5)') ' - ALBEDT  = ',albedt(i,j)
                  write (unit=*,fmt='(a,1x,es12.5)') ' - RAIN    = ',rain(i,j)
                  write (unit=*,fmt='(a,1x,es12.5)') ' - AREA    = ',area_cld
                  write (unit=*,fmt='(123a)') ('-',k=1,123)
                  write (unit=*,fmt='(10(a,1x))') '    K','  TEMP [degC]','  PRESS [hPa]'  &
                                                         ,'  DN0 [kg/m3]','    RV [g/kg]'  &
                                                         ,'   LWL [g/kg]','   IWL [g/kg]'  &
                                                         ,'LWLCLD [g/kg]','IWLCLD [g/kg]'  &
                                                         ,'FTHRD [K/day]'
                  write (unit=*,fmt='(123a)') ('-',k=1,123)
                  do k=1,m1
                     press = ((pp(k,i,j) + pi0(k,i,j)) * cpi) ** cpor * p00 * 0.01
                     tempc = (theta(k,i,j) * (pp(k,i,j) + pi0(k,i,j)) * cpi) - t00
                     rvcgs = rv(k,i,j) * 1000.
                     write (unit=*,fmt='(i5,1x,8(f13.6,1x),es13.6,1x)')                    &
                                                                k, tempc     , press       &
                                                                 , dn0(k,i,j), rvcgs       &
                                                                 , lwl(k,i,j), iwl(k,i,j)  &
                                                                 , lwl_cld(k), iwl_cld(k)  &
                                                                 , fthrd_cld(k)
                  end do
                  write (unit=*,fmt='(123a)') ('=',k=1,123)
                  write (unit=*,fmt='(a)')     ' '
                  call abort_run('Radiative heating rate is weird...','radcomp_carma'      &
                                ,'rad_carma.f90')
               end if
               !---------------------------------------------------------------------------!
            end do cldloop

            !----- Copying the scratch array back to the CARMA structure. -----------------!
            call ci_1d_3d(m2,m3,nwave,aotl,carma(ngrid)%aot,i,j)
            
            if (ir_aerad .and. rlong(i,j) /= rlong(i,j)) then
               write (unit=*,fmt='(a)')       '==========================================='
               write (unit=*,fmt='(a,1x,i5)')     ' MYNUM   = ',mynum
               write (unit=*,fmt='(a,1x,i5)')     ' IA      = ',ia
               write (unit=*,fmt='(a,1x,i5)')     ' IZ      = ',iz
               write (unit=*,fmt='(a,1x,i5)')     ' JA      = ',ja
               write (unit=*,fmt='(a,1x,i5)')     ' JZ      = ',jz
               write (unit=*,fmt='(a,1x,i5)')     ' I       = ',i
               write (unit=*,fmt='(a,1x,i5)')     ' J       = ',j
               write (unit=*,fmt='(a,1x,es12.5)') ' RSHORT  = ',rshort(i,j)
               write (unit=*,fmt='(a,1x,es12.5)') ' RLONG   = ',rlong(i,j)
               write (unit=*,fmt='(a,1x,es12.5)') ' RLONGUP = ',rlongup(i,j)
               write (unit=*,fmt='(a)')       '==========================================='
               call abort_run('Weird RLONG on CARMA...','radcomp_carma','rad_carma.f90')
            end if
         end do iloop
      end do jloop

      return
   end subroutine radcomp_carma
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   subroutine radcarma(nzpmax,m1,solfac,theta,pi0,pp,rv,rain,lwl,iwl,lwl_cld,iwl_cld,dn0   &
                      ,rtp,fthrd,rtgt,f13t,f23t,glat,glon,rshort,rshort_top,rshortup_top   &
                      ,rlong,rlongup_top,albedt,cosz,rlongup,mynum,fmapt,pm,aotl,xland)
      use catt_start  , only : catt       ! ! intent(in)
      use mem_grid    , only : centlon    & ! intent(in)
                             , dzm        & ! intent(in)
                             , dzt        & ! intent(in)
                             , idatea     & ! intent(in)
                             , imontha    & ! intent(in)
                             , itimea     & ! intent(in)
                             , itopo      & ! intent(in)
                             , iyeara     & ! intent(in)
                             , ngrid      & ! intent(in)
                             , nzp        & ! intent(in)
                             , plonn      & ! intent(in)
                             , time       ! ! intent(in)
      use mem_radiate , only : lonrad     & ! intent(in)
                             , iswrtyp    & ! intent(in)
                             , ilwrtyp    ! ! intent(in)
      use rconstants  , only : cp         & ! intent(in)
                             , cpor       & ! intent(in)
                             , p00        & ! intent(in)
                             , pio180     & ! intent(in)
                             , halfpi     & ! intent(in)
                             , wdns       & ! intent(in)
                             , stefan     & ! intent(in)
                             , toodry     ! ! intent(in)
      use mem_aerad   , only : ngas       & ! intent(in)
                             , nwave      & ! intent(in)
                             , lprocopio  ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      integer                  , intent(in)    :: nzpmax
      integer                  , intent(in)    :: m1
      integer                  , intent(in)    :: mynum
      real                     , intent(in)    :: solfac
      real                     , intent(in)    :: rtgt
      real                     , intent(in)    :: f13t
      real                     , intent(in)    :: f23t
      real                     , intent(in)    :: glat
      real                     , intent(in)    :: glon
      real                     , intent(in)    :: cosz
      real                     , intent(in)    :: albedt
      real                     , intent(in)    :: fmapt
      real                     , intent(in)    :: rain
      real                     , intent(in)    :: xland
      real, dimension(m1)      , intent(in)    :: pm     ! part. mat. (kg[pm]/kg[air])
      real, dimension(m1)      , intent(in)    :: theta
      real, dimension(m1)      , intent(in)    :: pi0
      real, dimension(m1)      , intent(in)    :: pp
      real, dimension(m1)      , intent(in)    :: rv
      real, dimension(m1)      , intent(in)    :: dn0
      real, dimension(m1)      , intent(in)    :: rtp
      real, dimension(m1)      , intent(in)    :: lwl
      real, dimension(m1)      , intent(in)    :: iwl
      real, dimension(m1)      , intent(in)    :: lwl_cld
      real, dimension(m1)      , intent(in)    :: iwl_cld
      real                     , intent(in)    :: rlongup
      real, dimension(m1)      , intent(inout) :: fthrd
      real, dimension(nwave)   , intent(inout) :: aotl
      real                     , intent(out)   :: rshort
      real                     , intent(out)   :: rshort_top
      real                     , intent(out)   :: rshortup_top
      real                     , intent(out)   :: rlong
      real                     , intent(out)   :: rlongup_top
      !----- Local variables. -------------------------------------------------------------!
      real, dimension(nzpmax)                  :: prd
      real, dimension(nzpmax+1)                :: temprd
      real, dimension(nzpmax)                  :: dn0r
      real, dimension(nzpmax)                  :: dztr
      real, dimension(nzpmax)                  :: pmr
      real, dimension(nzpmax)                  :: rvr
      real                                     :: rainr
      real, dimension(nzpmax)                  :: lwlr
      real, dimension(nzpmax)                  :: iwlr
      real                                     :: xlandr
      real, dimension(nzpmax)                  :: fthrl
      real, dimension(nzpmax)                  :: fthrs
      real                                     :: pird
      real                                     :: dzsdx
      real                                     :: dzsdy
      real                                     :: dlon
      real                                     :: a1
      real                                     :: a2
      real                                     :: dayhr
      real                                     :: gglon
      real                                     :: dztri
      real                                     :: dayhrr
      real                                     :: hrangl
      real                                     :: sinz
      real                                     :: sazmut
      real                                     :: slazim
      real                                     :: slangl
      real                                     :: cosi
      integer                                  :: igas
      integer                                  :: kk
      integer                                  :: ik
      integer                                  :: i
      integer                                  :: j
      integer                                  :: k
      integer                                  :: nzz
      !----- Local constants. -------------------------------------------------------------!
      real                     , parameter     :: fcui = 1.e-6 ! mg/kg => kg/kg
      !------------------------------------------------------------------------------------!


      !----- Copying environment variables to some scratch arrays. ------------------------!
      nzz = m1 - 1
      do k = 1,m1
         pird      = (pp(k) + pi0(k)) / cp
         temprd(k) = theta(k) * pird ! air temperature (k)
         rvr(k)    = max(toodry,rv(k))

         !----- Convert the next 7 variables to cgs for the time being. ----------------!
         prd(k)  = pird ** cpor * p00 * 10.  ! pressure
         dn0r(k) = dn0(k) * 1.e-3         ! air density
         dztr(k) = dzt(k) / rtgt * 1.e-2

         if (catt == 1) then
            pmr(k) = pm(k) * fcui
         else
            pmr(k) = 0.
         end if
      end do

      !----- Finding the lower boundary condition. ----------------------------------------!
      temprd(1) = (rlongup / stefan) ** 0.25
      temprd(nzp+1) = temprd(nzp)
      !----- Initialise the surface atmospheric state. ---------------------------------!
      p_surf = 0.5*(prd(1) + prd(2))
      p_top  = prd(m1)
      t_surf = temprd(1)

      !------------------------------------------------------------------------------------!
      !      K level on CARMA grid corresponds to k+1 level on BRAMS grid.  Here we        !
      ! transfer values from BRAMS grid to CARMA grid
      !------------------------------------------------------------------------------------!
      do k=1,m1-1
         p(k)    =    prd(k+1)
         t(k)    = temprd(k+1)
         rhoa(k) =   dn0r(k+1)
         lwlr(k) = (lwl(k+1)+lwl_cld(k+1)) * dn0r(k+1) * wdns  ![kg/m3]
         iwlr(k) = (iwl(k+1)+iwl_cld(k+1)) * dn0r(k+1) * wdns  ![kg/m3]
      end do


      !----- Revert the vertical index when cartesian coordinates are used. ---------------!
      do ik = 1,nzz
         kk = nzz + 1 - ik
         t_aerad(kk)   = t(ik)
         p_aerad(kk)   = p(ik)
         lwl_aerad(kk) = lwlr(ik)
         iwl_aerad(kk) = iwlr(ik)
      end do
      do ik = 1,nzz
         dztri=1./(dztr(ik) * 1.e+2)
         lwp_aerad(ik) = lwl_aerad(ik) * dztri   ![kg/m2]
         iwp_aerad(ik) = iwl_aerad(ik) * dztri   ![kg/m2]
      end do

      xland_aerad=xland
      rain_aerad=rain
      tabove_aerad  = t(nzz)

      !----- Initialise gas concentrations. -----------------------------------------------!
      do igas = 1,ngas
         do k = 1,m1-1
            !----- K level on CARMA grid corresponds to k+1 level on BRAMS grid. ----------!
            if (igas == 1) gc(k,igas) = rvr(k+1) * dn0r(k+1)
         end do
      end do

 
      call initaer(m1,pmr,dn0r,nzpmax)

      !----- Initialise radiation arrays. -------------------------------------------------!
      call initrad(imontha,idatea,iyeara,itimea,time,m1)


      call prerad(m1,nzpmax,dztr,fmapt)
      call radtran(albedt,cosz,m1,aotl(11))
      call radtran_to_rams(m1,fthrl,rlong,rlongup_top,fthrs,rshort,rshort_top,rshortup_top &
                          ,aotl,mynum)

      !------------------------------------------------------------------------------------!
      !    Modify the downward surface shortwave flux by considering the slope of the      !
      ! topography.                                                                        !
      !------------------------------------------------------------------------------------!
      if (itopo == 1 .and. isl_aerad) then
         dzsdx = f13t * rtgt
         dzsdy = f23t * rtgt
  
         !---------------------------------------------------------------------------------!
         !     The y- and x-directions must be true north and east for this correction.    !
         ! The following lines will rotate the model y/x to the true north/east.           !
         !---------------------------------------------------------------------------------!
         dlon =   (plonn(ngrid) - glon) * pio180
         a1   =   dzsdx * cos(dlon) + dzsdy * sin(dlon)
         a2   = - dzsdx * sin(dlon) + dzsdy * cos(dlon)
         dzsdx = a1
         dzsdy = a2
  
         dayhr = real(time / 3600.) + real(itimea/100)  + real(mod(itimea,100)) / 60.
         gglon = glon
         if (lonrad == 0) gglon = centlon(1)
         dayhrr = mod(dayhr+gglon/15.+24.,24.)
         hrangl = 15. * (dayhrr - 12.) * pio180
         sinz   = max(0.000001,sqrt(max(0., (1. - cosz * cosz ) ) ))
  
         sazmut = asin(max(-1.,min(1.,cdec*sin(hrangl)/sinz)))
         !---------------------------------------------------------------------------------!
         !    Imposing a lower bound for dzsdx and dzsdy, they will be squared soon after  !
         ! and we don't want it to cause underflow.                                        !
         !---------------------------------------------------------------------------------!
         if (abs(dzsdx) < 1e-16) dzsdx = 1.e-16
         if (abs(dzsdy) < 1e-16) dzsdy = 1.e-16
         slazim       = halfpi - atan2(dzsdy,dzsdx)
         slangl       = atan(sqrt(dzsdx*dzsdx+dzsdy*dzsdy))
         cosi         = cos(slangl) * cosz + sin(slangl) * sinz * cos(sazmut-slazim)
         rshort       = max(0.,rshort       * cosi / cosz)
      end if
    
      do k = 2,m1-1
         fthrd(k) = fthrl(k) + fthrs(k)
      end do
  
      rshort  = rshort / (1. - albedt)
      fthrd(1) = fthrd(2)
      
      return
   end subroutine radcarma
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This routine evaluates the derived mapping arrays and sets up the particle size   !
   ! bins.                                                                                 !
   !---------------------------------------------------------------------------------------!
   subroutine setupbins
      use rconstants , only : pi4o3         & ! intent(in)
                            , pi1           & ! intent(in)
                            , onethird      ! ! intent(in)
      use mem_aerad  , only : lunoprt       & ! intent(something)
                            , nbin          ! ! intent(something)
      use mem_globaer, only : ngroup        & ! intent(something)
                            , nelem         & ! intent(something)
                            , itype         & ! intent(something)
                            , i_involatile  & ! intent(something)
                            , i_volatile    & ! intent(something)
                            , ienconc       & ! intent(something)
                            , igelem        & ! intent(something)
                            , ncore         & ! intent(something)
                            , nelemg        & ! intent(something)
                            ,i_coremass     & ! intent(something)
                            ,i_volcore      & ! intent(something)
                            ,i_core2mom     & ! intent(something)
                            ,ixyz           & ! intent(something)
                            ,nxyz           & ! intent(something)
                            ,rhop3          & ! intent(something)
                            ,rhoelem        & ! intent(something)
                            ,rhopcore3      & ! intent(something)
                            ,rhocore        & ! intent(something)
                            ,rmassmin       & ! intent(something)
                            ,rmin           & ! intent(something)
                            ,rmrat          & ! intent(something)
                            ,rmass          & ! intent(something)
                            ,rmasscore      & ! intent(something)
                            ,pcore          & ! intent(something)
                            ,rmassup        & ! intent(something)
                            ,rmasscoreup    & ! intent(something)
                            ,dm             & ! intent(something)
                            ,vol            & ! intent(something)
                            ,r              & ! intent(something)
                            ,rcore          & ! intent(something)
                            ,rup            & ! intent(something)
                            ,rcoreup        & ! intent(something)
                            ,dr             & ! intent(something)
                            ,rlow           & ! intent(something)
                            ,diffmass       ! ! intent(something)
      implicit none
     
      !----- Local variables. -------------------------------------------------------------!
      integer :: igrp
      integer :: ielem
      integer :: j
      integer :: ie
      integer :: ig
      integer :: ibin
      real    :: vrfact
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Determine which elements are particle number concentrations.                  !
      ! <ienconc(igroup)> is the element corresponding to particle number concentration in !
      ! group <igroup>.                                                                    !
      !------------------------------------------------------------------------------------!
      igrp = 0
      do ielem = 1, nelem
         if (itype(ielem) == i_involatile .or. itype(ielem) .eq. i_volatile) then
            igrp = igrp + 1
            ienconc(igrp) = ielem
         end if
      end do

      !------------------------------------------------------------------------------------!
      !      Determine which group each element belongs to, i.e., <igelem(ielem)> is the   !
      ! group to which element <ielem> belongs.                                            !
      !------------------------------------------------------------------------------------!
      igrp = 0
      do ielem = 1, nelem
         if (itype(ielem) == i_involatile .or. itype(ielem) .eq. i_volatile) then
            igrp = igrp + 1
         end if
         igelem(ielem) = igrp
      end do

      !------------------------------------------------------------------------------------!
      !       Particle mass densities (nxyz*nbin for each group) -- the user might want    !
      !  to modify this (this code segment does not appear in setupaer subroutine          !
      !  because <igelem> is not defined until this subroutine).                           !
      !------------------------------------------------------------------------------------!

      do ie = 1,nelem
         ig = igelem(ie)
         do ibin = 1,nbin
            do ixyz = 1,nxyz
               rhop3(ixyz,ibin,ig)     = rhoelem(ie)
               rhopcore3(ixyz,ibin,ig) = rhocore(ie)
            end do
         end do
      end do


      !------------------------------------------------------------------------------------!
      !      Set up the particle bins.  For each particle group, the mass of a particle in !
      !  bin j is <rmrat> times that in bin j-1.                                           !
      !                                                                                    !
      !    rmass(nbin,ngroup)     =  bin center mass [g]                                   !
      !    r(nbin,ngroup)         =  bin mean (volume-weighted) radius [cm]                !
      !    vol(nbin,ngroup)       =  bin center volume [cm^3]                              !
      !    dr(nbin,ngroup)        =  bin width in radius space [cm]                        !
      !    dv(nbin,ngroup)        =  bin width in volume space [cm^3]                      !
      !    dm(nbin,ngroup)        =  bin width in mass space [g]                           !
      !------------------------------------------------------------------------------------!
      do igrp = 1, ngroup

         rmassmin(igrp) = pi4o3 * rhop3(1,1,igrp) * rmin(igrp)**3
         vrfact = ( (3./2./pi1/(rmrat(igrp)+1.))**onethird )* (rmrat(igrp)**onethird - 1.)

         do j = 1, nbin
            rmass(j,igrp)     = rmassmin(igrp) * rmrat(igrp)**(j-1)
            rmasscore(j,igrp) = pcore/100.     * rmass(j,igrp)
   
            rmassup(j,igrp)     = 2. * rmrat(igrp) / (rmrat(igrp)+1.)*rmass(j,igrp)
            rmasscoreup(j,igrp) = pcore/100. * rmassup(j,igrp)
   
            dm(j,igrp)          = 2.*(rmrat(igrp)-1.)/(rmrat(igrp)+1.)*rmass(j,igrp)
            vol(j,igrp)         = rmass(j,igrp) / rhop3(1,1,igrp)
   
            r(j,igrp)           = (rmass(j,igrp)/rhop3(1,1,igrp)/pi4o3)**onethird
            rcore(j,igrp)       = (rmasscore(j,igrp)/rhopcore3(1,1,igrp)/pi4o3 )**onethird
   
            rup(j,igrp)         = (rmassup(j,igrp)/rhop3(1,1,igrp)/pi4o3)**onethird
            rcoreup(j,igrp)     =  (rmasscoreup(j,igrp)/rhopcore3(1,1,igrp)/pi4o3 )        &
                                ** onethird

            dr(j,igrp)          = vrfact*(rmass(j,igrp)/rhop3(1,1,igrp))**onethird
            rlow(j,igrp)        = rup(j,igrp) - dr(j,igrp)
         end do
      end do

      return
   end subroutine setupbins
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   subroutine initaer(m1,pmr,dn0r,nzpmax)
      use mem_aerad   , only : nbin        ! !
      use mem_globaer , only : nelem       & !
                             , igelem      & !
                             , ienconc     & !
                             , small_pc    & !
                             , itype       & !
                             , i_coremass  & !
                             , rmass       & !
                             , fix_coref   & !
                             , i_core2mom  & !
                             , rhop3       & !
                             , dr          & !
                             , r           ! !
      use rconstants  , only : pi1         & ! intent(in)
                             , twopi       ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      integer                   , intent(in) :: m1,nzpmax
      real   , dimension(nzpmax), intent(in) :: pmr
      real   , dimension(nzpmax), intent(in) :: dn0r
      !----- Local variables. -------------------------------------------------------------!
      real   , dimension(m1)                 :: totm
      real   , dimension(m1)                 :: r0
      real   , dimension(m1)                 :: rsig
      real                                   :: sum
      real                                   :: totn
      integer                                :: ie,ix,iy
      integer                                :: ielem
      integer                                :: ig
      integer                                :: ip
      integer                                :: j
      integer                                :: k
      integer                                :: kr
      integer                                :: nzz
      real                                   :: arg1
      real                                   :: arg2
      !------------------------------------------------------------------------------------!

      nzz = m1 - 1

      !----- Transfer variables from BRAMS grid to CARMA. ---------------------------------!
      do k = 1,nzz
         !---------------------------------------------------------------------------------!
         !     k level in CARMA corresponds to the k+1 level in BRAMS.  This is true for   !
         ! terrain-following coordinate, not for the adaptive one...                       !
         !---------------------------------------------------------------------------------!
         kr = k + 1
         !---- totm  is the total mass particle concentration (g/cm3). --------------------!
         totm(k)  = pmr(kr) *  dn0r(kr) 
      end do

      !------------------------------------------------------------------------------------!
      !      Initialize particle number densities.  Core mass is assumed to be 100% of     !
      ! particle mass.                                                                     !
      !------------------------------------------------------------------------------------!
      do ielem = 1,nelem
         ig = igelem(ielem)
         ip = ienconc(ig)
         do j = 1,nbin
            do k = 1,nzz
               if (ielem == ip) then
                  !----- Particle number concentration [#/cm3]. ---------------------------!
                  pc(k,j,ielem) = small_pc
               elseif (itype(ielem) == i_coremass) then
                  !----- Core mass concentration [g/cm3]. ---------------------------------!
                  pc(k,j,ielem) = pc(k,j,ip) * rmass(j,ig) * fix_coref
               else if( itype(ielem) .eq. i_core2mom )then
                  !----- Second moment of core mass distribution [g2/cm6]. ----------------!
                  pc(k,j,ielem) = pc(k,j,ip) * (rmass(j,ig)*fix_coref)**2
               end if
            end do
         end do
      end do

      !------------------------------------------------------------------------------------!
      !      Initial particle distribution: log-normal size distribution for first         !
      ! particle group (which has only one particle element) in a single column.           !
      !------------------------------------------------------------------------------------!
      ig = 1
      ie = ienconc(ig)
  
      do k = 1,nzz
         !---------------------------------------------------------------------------------!
         !    Log-normal parameters:                                                       !
         !                                                                                 !
         !    r0   = number mode radius.                                                   !
         !    rsig = geometric standard deviation.                                         !
         !    totm = total mass particle concentration (g/cm3).                            !
         !---------------------------------------------------------------------------------!
         r0(k)   = 1.95e-5
         rsig(k) = 1.62
         totn    = (6. * totm(k) / (rhop3(1,1,1)*pi1*r0(k)**3) )                           &
                 * exp((-9./2)*log(rsig(k))**2)
         !----- Adjust prefactor to yield particle number concentration <ntot>. -----------!
         sum = 0.
         do j = 1,nbin
            arg1 = dr(j,ig) / ( sqrt(twopi) * r(j,ig) * log(rsig(k)) ) 
            arg2 = - log( r(j,ig) / r0(k) )**2 / ( 2.*log(rsig(k))**2 )
            sum  = sum + arg1 * exp( arg2 )
         end do

         totn = totn / sum

         do j = 1,nbin
            arg1       = totn * dr(j,ig) / ( sqrt(twopi) * r(j,ig) * log(rsig(k)) ) 
            arg2       = -log( r(j,ig) / r0(k) )**2 / ( 2.*log(rsig(k))**2 )
            pc(k,j,ie) = max( arg1 * exp( arg2 ), small_pc )
         end do
      end do

      return
   end subroutine initaer
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   subroutine initrad(imonth1,idate1,iyear1,itime1,time_rams,m1)
      use mem_aerad   , only : is_grp_ice_aerad  & ! intent(something)
                             , r_aerad           & ! intent(something)
                             , rup_aerad         & ! intent(something)
                             , rcore_aerad       & ! intent(something)
                             , rcoreup_aerad     & ! intent(something)
                             , ptop_aerad        & ! intent(something)
                             , pbot_aerad        & ! intent(something)
                             , u0_aerad          & ! intent(something)
                             , sfc_alb_aerad     & ! intent(something)
                             , emisir_aerad      & ! intent(something)
                             , tsfc_aerad        & ! intent(something)
                             , tabove_aerad      & ! intent(something)
                             , wave_aerad        & ! intent(something)
                             , lprocopio         & ! intent(something)
                             , nx                & ! intent(something)
                             , ny                & ! intent(something)
                             , nbin              & ! intent(something)
                             , nwave             & ! intent(something)
                             , nsol              ! ! intent(something)
      use mem_globaer , only : time              & ! intent(something)
                             , do_solar          & ! intent(something)
                             , do_ir             & ! intent(something)
                             , rad_start         & ! intent(something)
                             , ix                & ! intent(something)
                             , iy                & ! intent(something)
                             , rlat              & ! intent(something)
                             , u0                & ! intent(something)
                             , ngroup            & ! intent(something)
                             , is_grp_ice        & ! intent(something)
                             , ienconc           & ! intent(something)
                             , r                 & ! intent(something)
                             , rup               & ! intent(something)
                             , rcore             & ! intent(something)
                             , rcoreup           & ! intent(something)
                             , t_surf            & ! intent(something)
                             , wave              & ! intent(something)
                             , z_sin             & ! intent(something)
                             , z_cos             ! ! intent(something)
      use mem_globrad , only : lmie              ! ! intent(something)
      use mem_carma   , only : declin            ! ! intent(something)
      use rconstants  , only : day_sec           & ! intent(in)
                             , pi1               ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      integer     , intent(in) :: m1
      integer     , intent(in) :: imonth1
      integer     , intent(in) :: idate1
      integer     , intent(in) :: iyear1
      integer     , intent(in) :: itime1
      real(kind=8), intent(in) :: time_rams
      !----- Local variables. -------------------------------------------------------------!
      integer                  :: iday
      integer                  :: iwave
      integer                  :: julday
      real                     :: saz
      real                     :: wavetemp
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !    Define values needed for calculation of solar zenith angle:                     !
      !    - iday      : is day of year                                                    !
      !    - rad_start : is solar time corresponding to <time> = 0, in seconds             !
      !                = 0 means <time> = 0 corresponds to midnight,                       !
      !                = 6 * 3600 means <time> = 0 corresponds to 6 am                     !
      !    Note: all times are local standard time.                                        !
      !------------------------------------------------------------------------------------!
      iday      = julday(imonth1,idate1,iyear1)
      iday      = iday + nint(time_rams/day_sec)
      rad_start = (float(itime1/100) + float(mod(itime1,100)) / 60.)*day_sec


      !------------------------------------------------------------------------------------!
      !    Calculate terms in solar zenith angle computation (adapted from original Toon's !
      ! model).                                                                            !
      !  - saz: is solar azimuth angle [rad]                                               !
      !  - declin: is solar declination [rad]                                              !
      !  - z_sin: is sin term in precalculation                                            !
      !  - z_cos: is cos term in precalculation                                            !
      !------------------------------------------------------------------------------------!
  
      saz = 2. * pi1 / 365. * iday 
  
      declin = 0.006918 - 0.399912 * cos(     saz) + 0.070257 * sin(     saz)              &
                        - 0.006758 * cos(2. * saz) + 0.000907 * sin(2. * saz)              &
                        - 0.002697 * cos(3. * saz) + 0.001480 * sin(3. * saz)

      !------------------------------------------------------------------------------------!
      !    Initialize the radiative transfer model.                                        !
      !------------------------------------------------------------------------------------!
      call setuprad(m1)
     
       if (.not. lmie) call calcproperties

      !---- Get radiative wavelengths. ----------------------------------------------------!
      do iwave = 1,nwave
         !---------------------------------------------------------------------------------!
         !      Solar wavelengths in radiative transfer model are bin centers, infrared    !
         ! are bin edges.                                                                  !
         !---------------------------------------------------------------------------------!
         if (iwave <= nsol) then
            wave(iwave) = wave_aerad(iwave)
         else
            wave(iwave) = 0.5*( wave_aerad(iwave) + wave_aerad(iwave+1) )
         end if
      end do

      !---- Switch bins 11 and 12 (stored at 17 and 18???). -------------------------------!
      wavetemp = wave(17)
      wave(17) = wave(18)
      wave(18) = wavetemp

      return
   end subroutine initrad
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine defines various constants, and calculates the pressure-averaged   !
   ! absorption coefficients.                                                              !
   !---------------------------------------------------------------------------------------!
   subroutine setuprad(m1)
      use rconstants , only : gcgs            & ! intent(in)
                            , mmo2cgs         & ! intent(in)
                            , mmo2cgs         & ! intent(in)
                            , mmo3cgs         & ! intent(in)
                            , mmco2cgs        & ! intent(in)
                            , mmdrycgs        & ! intent(in)
                            , avogrado        & ! intent(in)
                            , loschcgs        ! ! intent(in)
      use mem_aerad  , only : wave_aerad      & ! intent(something)
                            , u0_aerad        & ! intent(something)
                            , sfc_alb_aerad   & ! intent(something)
                            , emisir_aerad    & ! intent(something)
                            , ptop_aerad      & ! intent(something)
                            , pbot_aerad      & ! intent(something)
                            , tsfc_aerad      ! ! intent(something)
      use mem_globrad, only : nlayer          & ! intent(something)
                            , ntotal          & ! intent(something)
                            , nprob           & ! intent(something)
                            , wave            & ! intent(something)
                            , nvert           & ! intent(something)
                            , p               & ! intent(something)
                            , t               & ! intent(something)
                            , o3c             & ! intent(something)
                            , wol             & ! intent(something)
                            , gol             & ! intent(something)
                            , tauray          & ! intent(something)
                            , nsolp           & ! intent(something)
                            , ltemp           & ! intent(something)
                            , nsol            & ! intent(something)
                            , aco2            & ! intent(something)
                            , xaco2           & ! intent(something)
                            , ao2             & ! intent(something)
                            , xao2            & ! intent(something)
                            , ao3             & ! intent(something)
                            , xao3            & ! intent(something)
                            , xah2o           & ! intent(something)
                            , psh2o           & ! intent(something)
                            , psco2           & ! intent(something)
                            , pso2            & ! intent(something)
                            , pso3            & ! intent(something)
                            , pj              & ! intent(something)
                            , o3mixp          & ! intent(something)
                            , akh2o           & ! intent(something)
                            , ako3            & ! intent(something)
                            , akco2           & ! intent(something)
                            , nirp            & ! intent(something)
                            , lmie            & ! intent(something)
                            , corereal        & ! intent(something)
                            , coreimag        ! ! intent(something)
       
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      integer               ,intent(in) :: m1
      !----- Local variables. -------------------------------------------------------------!
      integer, dimension(nlayer)        :: ii
      integer, dimension(nlayer)        :: ik
      real   , dimension(m1)            :: pbar
      real   , dimension(m1)            :: o3mix
      real   , dimension(nlayer)        :: dp
      real   , dimension(nlayer)        :: ps
      integer                           :: i
      integer                           :: i1
      integer                           :: j1
      integer                           :: j
      integer                           :: k
      integer                           :: l
      real                              :: co2mix
      real                              :: o2mix
      real                              :: o3mix2
      real                              :: pm
      real                              :: wvo
      real                              :: x
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !    Compute the layer average pressure (pbar), pressure at the edge of the layer    !
      ! (press), and the layer total mass (dpg). Units are bars, dyne/cm2, and g/cm2,      !
      ! respectively.                                                                      !
      !------------------------------------------------------------------------------------!
      pbar (1) = p_top * 5.0e-7
      press(1) = p_top
      do k  = 2,nvert
         pbar(k)  = p_aerad(k-1) * 1.0e-6
         press(k) = (p_aerad(k-1) + p_aerad(k)) * 0.5
         dpg(k-1) = (press(k)-press(k-1)) / gcgs
      end do
      pbar(nlayer)  = p_aerad(nvert) * 1.0e-6
      press(nlayer) = p_surf
      dpg(nvert)    = (press(nlayer)-press(nvert)) / gcgs
      !----- Skin temperature. ------------------------------------------------------------!

      tt(nlayer) = t_surf

      !------------------------------------------------------------------------------------!
      !    This is the amount of water vapor above model domain (gm/cm2), based on the     !
      ! 1976 US standard atmosphere, mid-latitude sounding.                                !
      !------------------------------------------------------------------------------------!
      rdh2o(1)   = h2ocol_aerad
      !----- Interpolate temperature from layer centre (t) to layer edge (tt). ------------!
      tt(1) = t_aerad(1)
      do k = 2, nvert
         tt(k) =  t_aerad(k-1) * (press(k)/p_aerad(k-1))                                   &
               ** (log(t_aerad(k)/t_aerad(k-1))/log(p_aerad(k)/p_aerad(k-1)))
      end do
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Define mass mixing ratios.  O3mix taken from US standard atmosphere,           !
      ! mid-latitude sounding.  Likewise, ozone column abundance o3c (#/cm2) above ptop    !
      ! was calculated from this same profile.                                             !
      !------------------------------------------------------------------------------------!
      o2mix          = 0.22   * mmo2cgs  / mmdrycgs
      co2mix         = 3.5e-4 * mmco2cgs / mmdrycgs
      o3c            = 9.02e18
      o3mix2         = o3c * mmo3cgs * gcgs / (p_top * avogrado)

      do l = nsolp+1,ntotal
         ltemp(l-nsolp) = nprob(l) - nsol
      end do

      x =   loschcgs / avogrado

      !----- Convert solar absorption coefficients to cm**2/gm. ---------------------------!
      do l = 1,nsolp
        aco2(l) = xaco2(l) / (x * mmco2cgs)
        ao2(l)  = xao2(l)  / (x * mmo2cgs )
        ao3(l)  = xao3(l)  / (x * mmo3cgs )
      end do

      !----- Calculate absorption coefficients. -------------------------------------------!
      pah2o = 0.0
      paco2 = 0.0
      pao2  = 0.0
      pao3  = 0.0

      do j = 1,nlayer
         do l = 1,nsolp
            pah2o(l,j) = xah2o(l) * pbar(j) ** psh2o(l)
            paco2(l,j) =  aco2(l) * pbar(j) ** psco2(l)
            pao2(l,j)  =  ao2(l)  * pbar(j) ** pso2(l)
            pao3(l,j)  =  ao3(l)  * pbar(j) ** pso3(l)
         end do
      end do

      do j = 1,nlayer
         iloop: do i = 1,6
            ii(j) = i
            if(pbar(j) > pj(i)) exit iloop
         end do iloop
      end do
     
      ps = 0.0
      do j = 1,nlayer
         if (ii(j) == 1) ii(j) = 2
         dp(j) = log(pj(ii(j)-1)/pj(ii(j)))

         if (pbar(j) > pj(6)) then
            ik(j) = ii(j) - 1
         else
            ik(j) = ii(j)
         end if

         ps(j) = pbar(j)/pj(ik(j))
         if (j /= 1) then
            o3mix(j) =o3mixp(ik(j)) * ps(j)**(log(o3mixp(ii(j)-1) / o3mixp(ii(j))) / dp(j))
         end if
      end do
     
      do j = 1,nlayer

         do l = 1,31
            pah2o(nsolp+l,j) = akh2o(l,ik(j))                                              &
                             * ps(j)**(log(akh2o(l,ii(j)-1) / akh2o(l,ii(j))) / dp(j))
         end do

         do l = 32,35
            pah2o(nsolp+l,j) = akh2o(32,ik(j))                                             &
                             * ps(j)**(log(akh2o(32,ii(j)-1) / akh2o(32,ii(j))) / dp(j))
            pao3(nsolp+l,j)  = ako3(l-31,ik(j))                                            &
                             * ps(j)**(log(ako3(l-31,ii(j)-1) / ako3(l-31,ii(j))) / dp(j))
         end do

         pah2o(nsolp+36,j)  = akh2o(33,ik(j))                                              &
                            * ps(j)**(log(akh2o(33,ii(j)-1) / akh2o(33,ii(j))) / dp(j))

         do l = 37,40
            paco2(nsolp+l,j) = akco2(1,ik(j))                                              &
                             * ps(j)**(log(akco2(1,ii(j)-1) / akco2(1,ii(j))) / dp(j))
            pah2o(nsolp+l,j) = akh2o(l-3,ik(j))                                            &
                             * ps(j)**(log(akh2o(l-3,ii(j)-1) / akh2o(l-3,ii(j))) / dp(j))
         end do

         do l = 41,44
            paco2(nsolp+l,j) = akco2(2,ik(j))                                              &
                             * ps(j)**(log(akco2(2,ii(j)-1) / akco2(2,ii(j))) / dp(j))
            pah2o(nsolp+l,j) = akh2o(l-7,ik(j))                                            &
                             * ps(j)**(log(akh2o(l-7,ii(j)-1) / akh2o(l-7,ii(j))) / dp(j))
         end do

         do l = 45,48
            paco2(nsolp+l,j) = akco2(3,ik(j))                                              &
                             * ps(j)**(log(akco2(3,ii(j)-1) / akco2(3,ii(j))) / dp(j))
            pah2o(nsolp+l,j) = akh2o(l-11,ik(j))                                           &
                             * ps(j)**(log(akh2o(l-11,ii(j)-1)/akh2o(l-11,ii(j))) / dp(j))
         end do

         do l = 49,51
            paco2(nsolp+l,j) = akco2(4,ik(j))                                              &
                             * ps(j)**(log(akco2(4,ii(j)-1) / akco2(4,ii(j))) / dp(j))
            pah2o(nsolp+l,j) = akh2o(l-11,ik(j))                                           &
                             * ps(j)**(log(akh2o(l-11,ii(j)-1)/akh2o(l-11,ii(j))) / dp(j))
         end do

         do l = 52,54
            paco2(nsolp+l,j) = akco2(5,ik(j))                                              &
                             * ps(j)**(log (akco2(5,ii(j)-1) / akco2(5,ii(j))) / dp(j))
            pah2o(nsolp+l,j) = akh2o(l-14,ik(j))                                           &
                             * ps(j)**(log(akh2o(l-14,ii(j)-1)/akh2o(l-14,ii(j))) / dp(j))
         end do

         do l = 55,57
            paco2(nsolp+l,j) = akco2(6,ik(j))                                              &
                             * ps(j)**(log (akco2(6,ii(j)-1) / akco2(6,ii(j))) / dp(j))
            pah2o(nsolp+l,j) = akh2o(l-17,ik(j))                                           &
                             * ps(j)**(log(akh2o(l-17,ii(j)-1)/akh2o(l-17,ii(j))) / dp(j))
         end do

         do l = 58,nirp
            pah2o(nsolp+l,j) = akh2o(l-17,ik(j))                                           &
                             * ps(j)**(log(akh2o(l-17,ii(j)-1)/akh2o(l-17,ii(j))) / dp(j))
         end do
      end do

      !----- Store o3mix2 in o3mix(1). ----------------------------------------------------!
      o3mix(1) = o3mix2
     
      !----- Here we find taugas, which is the sum of tauco2, tauo2, and tauo3. -----------!
      do l=1,ntotal
         pm = p_top/gcgs
         taugas(l,1) = pm * (o2mix*pao2(l,1) + co2mix*paco2(l,1) + o3mix(1)*pao3(l,1)) 
      end do
      do j = 2, nlayer
         do l = 1,ntotal
            pm = dpg(j-1)
            taugas(l,j) = pm * (o2mix*pao2(l,j) + co2mix*paco2(l,j) + o3mix(j)*pao3(l,j)) 
         end do
      end do
   
      !------------------------------------------------------------------------------------!
      !     Wave length must be in microns.  Here we compute the Rayleigh optical depth    !
      ! parameters.                                                                        !
      !------------------------------------------------------------------------------------!
      do l = 1,ntotal
        wvo       = wave(nprob(l))
        tauray(l) = (8.46e-9/wvo**4) * ( 1.+0.0113/wvo**2+0.00013/wvo**4 )
      end do
   
      !------------------------------------------------------------------------------------!
      !     We don't include Rayleigh scattering in infrared.                              !
      !------------------------------------------------------------------------------------!
      do j = 1,nvert
         do l= 1,ntotal
            if (l <= nsolp) then
               paray(l,j+1) = tauray(l) * dpg(j) * gcgs
            else
               paray(l,j+1) = 0.
            end if
         end do
         do l   =   1,ntotal
            if (l <= nsolp) then
               paray(l,1) = tauray(l) * p_top
            else
               paray(l,1) = 0.
            end if
         end do
      end do

      return
   end subroutine setuprad
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
  
  SUBROUTINE calcproperties
    ! **********************************************************************
    !
    !		 CALCULATE THE AEROSOL EXTINCTION CROSS SECTIONS
    !
    ! **********************************************************************
    !
    !	  Get <is_grp_ice> and radius grid from interface common block
    !	  and calculate cross-sectional area for each bin.
    use rconstants, only : pii, pi1, twopi, stefan
    USE mem_aerad, ONLY: is_grp_ice_aerad,r_aerad,rcore_aerad,rup_aerad, &
         rcoreup_aerad,lunmie,lunoprt, tabove_aerad
  
    USE mem_globrad, ONLY: ngroup,nrad,core_rad,xsecta, &
         coreup_rad,l_write,l_read,nwave,rmin,rdqext, &
         qscat,qbrqs,nsol,wave,corereal,coreimag,nsolp, &
         weight,ntotal,sol,solfx,nprob,iblackbody_above, &
         t_above,ncount,nlow,planck
  
    USE mem_globaer, ONLY: r,rcore,rup,rcoreup,is_grp_ice
  
    IMPLICIT NONE
    
    INTEGER :: i
    INTEGER :: ibeyond_spectrum
    INTEGER :: ig
    INTEGER :: irefr
    INTEGER :: j
    INTEGER :: jj
    INTEGER :: k
    INTEGER :: l
    INTEGER :: mgroup
    INTEGER :: mrad
    INTEGER :: mwave
    LOGICAL :: all_ok
    REAL    :: awave
    REAL    :: corerad
    REAL    :: ctbrqs
    REAL    :: ddr
    REAL    :: ddrc
    REAL    :: qextd
    REAL    :: qscatd
    REAL    :: r_real
    REAL    :: rr
    REAL    :: sum
    REAL    :: sum1
    REAL    :: sum2
    REAL    :: t1
    REAL    :: thetd
    REAL    :: tmag
    REAL    :: v
    REAL    :: wvno
  
    CHARACTER(LEN=*),PARAMETER :: &
         lab355='(//,"setuprad: error in weights ",/," ' &
         //'sum of weights for solar =",1pe15.5,/,"' &
         //' sum of weights for ir = ",1pe15.5,/," ' &
         //'total sum =  ",1pe15.5)'
    
    DO ig = 1, ngroup
       DO I = 1, nrad
          xsecta(i,ig) = pi1 * r(i,ig)**2.
       END DO
    END DO
    !
    !	  Set <i_mie> = L_READ to WRITE the mie coefficients to a data file,
    !		      = L_WRITE to READ them
    !
    
    ! ------------------------------------------------------------------
    ! "mie.data" Parameters definitions
    mwave=50
    mrad=30
    mgroup=1
    rmin=(/0.1000000E-05/)
    !(1:30,1:1,1:50)
    rdqext(1:30,1:1,1:6)=reshape( (/                                 &
         0.1238996E-05, 0.9044850E-03, 0.1058328E-02, 0.1310231E-02, &
         0.1623058E-02, 0.2011839E-02, 0.2492962E-02, 0.3089402E-02, &
         0.3829632E-02, 0.4747202E-02, 0.5886665E-02, 0.7303554E-02, &
         0.9068833E-02, 0.1127584E-01, 0.1405085E-01, 0.1757244E-01, &
         0.2210986E-01, 0.2809986E-01, 0.3631130E-01, 0.4820201E-01, &
         0.6671599E-01, 0.9805208E-01, 0.1554057E+00, 0.2653104E+00, &
         0.4692753E+00, 0.7888026E+00, 0.1204756E+01, 0.1752262E+01, &
         0.2199486E+01, 0.2739768E+01, 0.2356571E-05, 0.8107370E-03, &
         0.9485792E-03, 0.1175354E-02, 0.1457698E-02, 0.1806319E-02, &
         0.2237227E-02, 0.2773230E-02, 0.3436975E-02, 0.4259836E-02, &
         0.5281698E-02, 0.6551373E-02, 0.8131125E-02, 0.1010238E-01, &
         0.1257303E-01, 0.1569186E-01, 0.1967621E-01, 0.2486489E-01, &
         0.3183057E-01, 0.4162014E-01, 0.5628673E-01, 0.8007501E-01, &
         0.1220185E+00, 0.2010520E+00, 0.3521606E+00, 0.6149174E+00, &
         0.9806738E+00, 0.1471468E+01, 0.1975219E+01, 0.2499236E+01, &
         0.4482199E-05, 0.6878105E-03, 0.8052588E-03, 0.9958311E-03, &
         0.1234672E-02, 0.1529180E-02, 0.1895103E-02, 0.2348306E-02, &
         0.2910469E-02, 0.3606677E-02, 0.4470860E-02, 0.5543770E-02, &
         0.6876979E-02, 0.8536757E-02, 0.1060989E-01, 0.1321133E-01, &
         0.1650270E-01, 0.2072231E-01, 0.2624893E-01, 0.3373386E-01, &
         0.4438854E-01, 0.6062073E-01, 0.8745370E-01, 0.1355955E+00, &
         0.2271680E+00, 0.4006258E+00, 0.6897491E+00, 0.1075740E+01, &
         0.1600417E+01, 0.2066453E+01, 0.8525142E-05, 0.5920068E-03, &
         0.6939178E-03, 0.8606558E-03, 0.1066172E-02, 0.1321362E-02, &
         0.1637412E-02, 0.2028739E-02, 0.2514813E-02, 0.3116922E-02, &
         0.3863455E-02, 0.4789023E-02, 0.5938904E-02, 0.7368284E-02, &
         0.9149896E-02, 0.1137746E-01, 0.1417908E-01, 0.1773632E-01, &
         0.2232307E-01, 0.2838602E-01, 0.3671290E-01, 0.4880296E-01, &
         0.6768916E-01, 0.9976597E-01, 0.1586369E+00, 0.2715197E+00, &
         0.4801895E+00, 0.8039033E+00, 0.1225015E+01, 0.1773117E+01, &
         0.1621482E-04, 0.5126877E-03, 0.5990996E-03, 0.7428692E-03, &
         0.9220786E-03, 0.1141519E-02, 0.1414444E-02, 0.1753329E-02, &
         0.2172156E-02, 0.2691440E-02, 0.3335611E-02, 0.4134388E-02, &
         0.5125767E-02, 0.6357079E-02, 0.7889097E-02, 0.9800074E-02, &
         0.1219303E-01, 0.1521039E-01, 0.1905701E-01, 0.2405036E-01, &
         0.3071964E-01, 0.4002312E-01, 0.5382273E-01, 0.7594469E-01, &
         0.1145147E+00, 0.1866735E+00, 0.3250482E+00, 0.5710782E+00, &
         0.9245852E+00, 0.1392398E+01, 0.3084059E-04, 0.4278156E-03, &
         0.4990920E-03, 0.6195636E-03, 0.7684358E-03, 0.9515302E-03, &
         0.1178236E-02, 0.1461325E-02, 0.1809853E-02, 0.2242170E-02, &
         0.2778809E-02, 0.3443943E-02, 0.4269208E-02, 0.5292946E-02, &
         0.6565067E-02, 0.8148359E-02, 0.1012359E-01, 0.1259963E-01, &
         0.1572569E-01, 0.1971983E-01, 0.2492236E-01, 0.3190922E-01, &
         0.4173374E-01, 0.5646301E-01, 0.8037232E-01, 0.1225615E+00, &
         0.2020946E+00, 0.3541155E+00, 0.6180223E+00, 0.9846209E+00/),&
         (/30,1,6/))

    rdqext(1:30,1:1,7:12)=reshape( (/&
         0.5865879E-04, 0.3386404E-03, 0.3943375E-03, 0.4902230E-03, &
         0.6069456E-03, 0.7545453E-03, 0.9325470E-03, 0.1156643E-02, &
         0.1432099E-02, 0.1774964E-02, 0.2199463E-02, 0.2725592E-02, &
         0.3378012E-02, 0.4187341E-02, 0.5191328E-02, 0.6438720E-02, &
         0.7990645E-02, 0.9926975E-02, 0.1235256E-01, 0.1541248E-01, &
         0.1931667E-01, 0.2439146E-01, 0.3118407E-01, 0.4068905E-01, &
         0.5484687E-01, 0.7765548E-01, 0.1176138E+00, 0.1926048E+00, &
         0.3362654E+00, 0.5893933E+00, 0.1115690E-03, 0.2516592E-03, &
         0.2934801E-03, 0.3651093E-03, 0.4550150E-03, 0.5632643E-03, &
         0.6970990E-03, 0.8643277E-03, 0.1070622E-02, 0.1324547E-02, &
         0.1642295E-02, 0.2034783E-02, 0.2521489E-02, 0.3124935E-02, &
         0.3872860E-02, 0.4801260E-02, 0.5953945E-02, 0.7387332E-02, &
         0.9173343E-02, 0.1140706E-01, 0.1421631E-01, 0.1778384E-01, &
         0.2238511E-01, 0.2846927E-01, 0.3683005E-01, 0.4897850E-01, &
         0.6797409E-01, 0.1002687E+00, 0.1595861E+00, 0.2733425E+00, &
         0.2122042E-03, 0.1571589E-03, 0.1840258E-03, 0.2313218E-03, &
         0.2901404E-03, 0.3512961E-03, 0.4390195E-03, 0.5448884E-03, &
         0.6724349E-03, 0.8375491E-03, 0.1036835E-02, 0.1285574E-02, &
         0.1592136E-02, 0.1972729E-02, 0.2444628E-02, 0.3029664E-02, &
         0.3754555E-02, 0.4654815E-02, 0.5771693E-02, 0.7160614E-02, &
         0.8890757E-02, 0.1105248E-01, 0.1376909E-01, 0.1721303E-01, &
         0.2164283E-01, 0.2747480E-01, 0.3543673E-01, 0.4689975E-01, &
         0.6461873E-01, 0.9437821E-01, 0.4036125E-03, 0.4093335E-01, &
         0.5522372E-01, 0.7828710E-01, 0.1187611E+00, 0.1948035E+00, &
         0.3404121E+00, 0.5961001E+00, 0.9566938E+00, 0.1437752E+01, &
         0.1951394E+01, 0.2464083E+01, 0.2949664E+01, 0.3406605E+01, &
         0.3629449E+01, 0.3356592E+01, 0.2565616E+01, 0.1934376E+01, &
         0.2454620E+01, 0.2495687E+01, 0.2086732E+01, 0.2298010E+01, &
         0.2203053E+01, 0.2162956E+01, 0.2135190E+01, 0.2162642E+01, &
         0.2130731E+01, 0.2120723E+01, 0.2111154E+01, 0.2073417E+01, &
         0.3014119E-01, 0.3649874E-01, 0.4848227E-01, 0.6716944E-01, &
         0.9884995E-01, 0.1569090E+00, 0.2682001E+00, 0.4743640E+00, &
         0.7958630E+00, 0.1214206E+01, 0.1762091E+01, 0.2209994E+01, &
         0.2747663E+01, 0.3190050E+01, 0.3513948E+01, 0.3566287E+01, &
         0.2919923E+01, 0.2062144E+01, 0.2168095E+01, 0.2654343E+01, &
         0.2093068E+01, 0.2346398E+01, 0.2183643E+01, 0.2237450E+01, &
         0.2201363E+01, 0.2169828E+01, 0.2122523E+01, 0.2131855E+01, &
         0.2098403E+01, 0.2083777E+01, 0.2825850E-01, 0.3407567E-01, &
         0.4489003E-01, 0.6141419E-01, 0.8881936E-01, 0.1381295E+00, &
         0.2320493E+00, 0.4095595E+00, 0.7030742E+00, 0.1092780E+01, &
         0.1622257E+01, 0.2083007E+01, 0.2637337E+01, 0.3099045E+01, &
         0.3466475E+01, 0.3576084E+01, 0.3138925E+01, 0.2289058E+01, &
         0.2024652E+01, 0.2619606E+01, 0.2271519E+01, 0.2340897E+01, &
         0.2215942E+01, 0.2211942E+01, 0.2123458E+01, 0.2147658E+01, &
         0.2145206E+01, 0.2089916E+01, 0.2105316E+01, 0.2094704E+01/),&
         (/30,1,6/))

    rdqext(1:30,1:1,13:18)=reshape( (/&
         0.2595173E-01, 0.3114211E-01, 0.4062882E-01, 0.5475407E-01, &
         0.7750003E-01, 0.1173315E+00, 0.1920643E+00, 0.3352450E+00, &
         0.5877376E+00, 0.9459963E+00, 0.1422653E+01, 0.1940509E+01, &
         0.2447865E+01, 0.2935348E+01, 0.3393993E+01, 0.3634520E+01, &
         0.3377689E+01, 0.2585614E+01, 0.1928602E+01, 0.2453908E+01, &
         0.2537595E+01, 0.2120584E+01, 0.2355100E+01, 0.2236805E+01, &
         0.2206528E+01, 0.2137048E+01, 0.2110582E+01, 0.2129888E+01, &
         0.2112728E+01, 0.2111520E+01, 0.2456848E-01, 0.2940105E-01, &
         0.3814531E-01, 0.5096052E-01, 0.7121034E-01, 0.1060126E+00, &
         0.1704660E+00, 0.2941915E+00, 0.5193896E+00, 0.8569373E+00, &
         0.1297529E+01, 0.1841395E+01, 0.2305595E+01, 0.2818348E+01, &
         0.3263292E+01, 0.3597307E+01, 0.3494461E+01, 0.2745187E+01, &
         0.2013707E+01, 0.2287942E+01, 0.2655969E+01, 0.2131976E+01, &
         0.2384275E+01, 0.2205176E+01, 0.2206340E+01, 0.2215608E+01, &
         0.2134920E+01, 0.2115479E+01, 0.2095288E+01, 0.2098952E+01, &
         0.2235298E-01, 0.2663950E-01, 0.3427472E-01, 0.4518270E-01, &
         0.6187829E-01, 0.8962012E-01, 0.1396183E+00, 0.2349173E+00, &
         0.4147888E+00, 0.7108101E+00, 0.1102706E+01, 0.1634750E+01, &
         0.2092761E+01, 0.2647524E+01, 0.3107261E+01, 0.3467156E+01, &
         0.3579113E+01, 0.3125170E+01, 0.2247379E+01, 0.2026859E+01, &
         0.2657282E+01, 0.2265156E+01, 0.2332558E+01, 0.2186288E+01, &
         0.2227778E+01, 0.2146957E+01, 0.2169617E+01, 0.2140218E+01, &
         0.2092722E+01, 0.2088475E+01, 0.1921333E-01, 0.2277928E-01, &
         0.2899945E-01, 0.3757729E-01, 0.5010217E-01, 0.6980440E-01, &
         0.1035102E+00, 0.1657181E+00, 0.2851043E+00, 0.5038010E+00, &
         0.8360577E+00, 0.1268737E+01, 0.1815424E+01, 0.2272147E+01, &
         0.2793396E+01, 0.3235432E+01, 0.3570348E+01, 0.3524326E+01, &
         0.2809205E+01, 0.2151345E+01, 0.2281744E+01, 0.2666198E+01, &
         0.2126539E+01, 0.2388358E+01, 0.2219515E+01, 0.2215051E+01, &
         0.2223578E+01, 0.2160036E+01, 0.2126196E+01, 0.2106740E+01, &
         0.1828648E-01, 0.2165078E-01, 0.2748548E-01, 0.3545161E-01, &
         0.4692184E-01, 0.6465408E-01, 0.9444001E-01, 0.1486191E+00, &
         0.2522533E+00, 0.4460807E+00, 0.7561581E+00, 0.1161547E+01, &
         0.1705030E+01, 0.2152755E+01, 0.2702824E+01, 0.3151649E+01, &
         0.3478999E+01, 0.3582580E+01, 0.3043263E+01, 0.2252439E+01, &
         0.2120131E+01, 0.2719469E+01, 0.2260474E+01, 0.2383116E+01, &
         0.2185176E+01, 0.2205244E+01, 0.2167593E+01, 0.2190734E+01, &
         0.2141221E+01, 0.2117435E+01, 0.1705720E-01, 0.2016090E-01, &
         0.2550504E-01, 0.3270831E-01, 0.4289185E-01, 0.5826783E-01, &
         0.8343055E-01, 0.1281675E+00, 0.2128689E+00, 0.3742239E+00, &
         0.6495301E+00, 0.1024606E+01, 0.1532292E+01, 0.2017566E+01, &
         0.2558653E+01, 0.3033171E+01, 0.3455450E+01, 0.3578496E+01, &
         0.3254357E+01, 0.2367896E+01, 0.1937363E+01, 0.2580973E+01, &
         0.2382489E+01, 0.2224710E+01, 0.2240751E+01, 0.2242849E+01, &
         0.2187788E+01, 0.2117765E+01, 0.2155853E+01, 0.2127444E+01/),&
         (/30,1,6/))

    rdqext(1:30,1:1,19:24)=reshape( (/&
         0.1661168E-01, 0.1962307E-01, 0.2479483E-01, 0.3173472E-01, &
         0.4148188E-01, 0.5607228E-01, 0.7971362E-01, 0.1213590E+00, &
         0.1997862E+00, 0.3497850E+00, 0.6111342E+00, 0.9758610E+00, &
         0.1464721E+01, 0.1970491E+01, 0.2492324E+01, 0.2974746E+01, &
         0.3425845E+01, 0.3616955E+01, 0.3324818E+01, 0.2534858E+01, &
         0.2057498E+01, 0.2525669E+01, 0.2476210E+01, 0.2140753E+01, &
         0.2303478E+01, 0.2236078E+01, 0.2172018E+01, 0.2146353E+01, &
         0.2130648E+01, 0.2142246E+01, 0.1598668E-01, 0.1886983E-01, &
         0.2380455E-01, 0.3038575E-01, 0.3954578E-01, 0.5309152E-01, &
         0.7472860E-01, 0.1123202E+00, 0.1824803E+00, 0.3170920E+00, &
         0.5579348E+00, 0.9075783E+00, 0.1368394E+01, 0.1899847E+01, &
         0.2387506E+01, 0.2883463E+01, 0.3339896E+01, 0.3636674E+01, &
         0.3422941E+01, 0.2665677E+01, 0.2101645E+01, 0.2395268E+01, &
         0.2593362E+01, 0.2087997E+01, 0.2348550E+01, 0.2197394E+01, &
         0.2205405E+01, 0.2186869E+01, 0.2118561E+01, 0.2120739E+01, &
         0.1447862E-01, 0.1705971E-01, 0.2144386E-01, 0.2720915E-01, &
         0.3506626E-01, 0.4635068E-01, 0.6373934E-01, 0.9284600E-01, &
         0.1456350E+00, 0.2465072E+00, 0.4357705E+00, 0.7413899E+00, &
         0.1142247E+01, 0.1682719E+01, 0.2132625E+01, 0.2685511E+01, &
         0.3137660E+01, 0.3472792E+01, 0.3583071E+01, 0.3077484E+01, &
         0.2215237E+01, 0.2114642E+01, 0.2624185E+01, 0.2220664E+01, &
         0.2292433E+01, 0.2218276E+01, 0.2214263E+01, 0.2160881E+01, &
         0.2179611E+01, 0.2127588E+01, 0.1328048E-01, 0.1562869E-01, &
         0.1959484E-01, 0.2475764E-01, 0.3168387E-01, 0.4140854E-01, &
         0.5595865E-01, 0.7952227E-01, 0.1210101E+00, 0.1991165E+00, &
         0.3485271E+00, 0.6091268E+00, 0.9733056E+00, 0.1461134E+01, &
         0.1967971E+01, 0.2488624E+01, 0.2971455E+01, 0.3423538E+01, &
         0.3618834E+01, 0.3328263E+01, 0.2539237E+01, 0.2039369E+01, &
         0.2503954E+01, 0.2453909E+01, 0.2115271E+01, 0.2290610E+01, &
         0.2230195E+01, 0.2159497E+01, 0.2120082E+01, 0.2116015E+01, &
         0.1200182E-01, 0.1410724E-01, 0.1764445E-01, 0.2220361E-01, &
         0.2822565E-01, 0.3648769E-01, 0.4846574E-01, 0.6714261E-01, &
         0.9880274E-01, 0.1568201E+00, 0.2680290E+00, 0.4740632E+00, &
         0.7954467E+00, 0.1213648E+01, 0.1761516E+01, 0.2209371E+01, &
         0.2747197E+01, 0.3189625E+01, 0.3513435E+01, 0.3566643E+01, &
         0.2920861E+01, 0.2063354E+01, 0.2168825E+01, 0.2651909E+01, &
         0.2093758E+01, 0.2346431E+01, 0.2182174E+01, 0.2229954E+01, &
         0.2199746E+01, 0.2167265E+01, 0.1181864E-01, 0.1388969E-01, &
         0.1736677E-01, 0.2184242E-01, 0.2774175E-01, 0.3580971E-01, &
         0.4745414E-01, 0.6550952E-01, 0.9593535E-01, 0.1514246E+00, &
         0.2576532E+00, 0.4557129E+00, 0.7698092E+00, 0.1179520E+01, &
         0.1725135E+01, 0.2171919E+01, 0.2718426E+01, 0.3164545E+01, &
         0.3487779E+01, 0.3580781E+01, 0.2997886E+01, 0.2219576E+01, &
         0.2137670E+01, 0.2651932E+01, 0.2223345E+01, 0.2425306E+01, &
         0.2227556E+01, 0.2204142E+01, 0.2194215E+01, 0.2185364E+01/),&
         (/30,1,6/))

    rdqext(1:30,1:1,25:30)=reshape( (/&
         0.1146874E-01, 0.1347448E-01, 0.1683753E-01, 0.2115573E-01, &
         0.2682521E-01, 0.3453235E-01, 0.4556208E-01, 0.6248134E-01, &
         0.9066294E-01, 0.1415598E+00, 0.2386575E+00, 0.4215866E+00, &
         0.7207971E+00, 0.1115563E+01, 0.1650671E+01, 0.2105534E+01, &
         0.2660306E+01, 0.3117507E+01, 0.3468289E+01, 0.3581810E+01, &
         0.3107697E+01, 0.2204360E+01, 0.2083869E+01, 0.2692604E+01, &
         0.2246163E+01, 0.2310285E+01, 0.2205017E+01, 0.2212293E+01, &
         0.2154623E+01, 0.2157726E+01, 0.1036541E-01, 0.1216810E-01, &
         0.1517874E-01, 0.1901641E-01, 0.2399688E-01, 0.3064706E-01, &
         0.3991923E-01, 0.5366334E-01, 0.7567919E-01, 0.1140351E+00, &
         0.1857564E+00, 0.3233098E+00, 0.5682176E+00, 0.9208928E+00, &
         0.1387183E+01, 0.1914258E+01, 0.2408734E+01, 0.2901393E+01, &
         0.3359867E+01, 0.3639189E+01, 0.3413174E+01, 0.2648792E+01, &
         0.2071824E+01, 0.2370969E+01, 0.2552859E+01, 0.2186358E+01, &
         0.2357531E+01, 0.2219193E+01, 0.2180621E+01, 0.2176074E+01, &
         0.1099698E-01, 0.1291546E-01, 0.1612651E-01, 0.2023638E-01, &
         0.2560486E-01, 0.3284562E-01, 0.4309161E-01, 0.5858051E-01, &
         0.8396278E-01, 0.1291468E+00, 0.2147529E+00, 0.3777217E+00, &
         0.6549321E+00, 0.1031459E+01, 0.1541618E+01, 0.2024099E+01, &
         0.2567300E+01, 0.3040632E+01, 0.3457704E+01, 0.3574554E+01, &
         0.3238655E+01, 0.2353804E+01, 0.1952011E+01, 0.2580399E+01, &
         0.2390371E+01, 0.2200474E+01, 0.2244411E+01, 0.2220684E+01, &
         0.2163651E+01, 0.2148918E+01, 0.9705304E-02, 0.1138785E-01, &
         0.1419229E-01, 0.1775317E-01, 0.2234502E-01, 0.2841550E-01, &
         0.3675437E-01, 0.4886510E-01, 0.6778988E-01, 0.9994365E-01, &
         0.1589722E+00, 0.2721637E+00, 0.4813170E+00, 0.8054546E+00, &
         0.1227107E+01, 0.1775224E+01, 0.2224483E+01, 0.2758411E+01, &
         0.3200041E+01, 0.3526455E+01, 0.3557312E+01, 0.2899145E+01, &
         0.2049803E+01, 0.2180376E+01, 0.2742225E+01, 0.2126863E+01, &
         0.2353465E+01, 0.2187709E+01, 0.2224508E+01, 0.2222815E+01, &
         0.9146382E-02, 0.1072809E-01, 0.1336029E-01, 0.1669216E-01, &
         0.2096744E-01, 0.2657457E-01, 0.3418472E-01, 0.4505027E-01, &
         0.6166824E-01, 0.8925748E-01, 0.1389438E+00, 0.2336179E+00, &
         0.4124213E+00, 0.7073137E+00, 0.1098216E+01, 0.1629121E+01, &
         0.2088338E+01, 0.2642952E+01, 0.3103581E+01, 0.3466837E+01, &
         0.3577819E+01, 0.3131474E+01, 0.2266839E+01, 0.2022272E+01, &
         0.2655272E+01, 0.2286961E+01, 0.2383289E+01, 0.2214466E+01, &
         0.2205394E+01, 0.2140238E+01, 0.8506889E-02, 0.9974504E-02, &
         0.1241214E-01, 0.1548789E-01, 0.1941376E-01, 0.2451917E-01, &
         0.3135827E-01, 0.4093947E-01, 0.5523325E-01, 0.7830309E-01, &
         0.1187902E+00, 0.1948592E+00, 0.3405170E+00, 0.5962692E+00, &
         0.9569100E+00, 0.1438057E+01, 0.1951612E+01, 0.2464408E+01, &
         0.2949951E+01, 0.3406850E+01, 0.3629331E+01, 0.3356170E+01, &
         0.2565260E+01, 0.1935007E+01, 0.2454862E+01, 0.2496692E+01, &
         0.2086216E+01, 0.2314208E+01, 0.2219513E+01, 0.2162504E+01/),&
         (/30,1,6/))

    rdqext(1:30,1:1,31:36)=reshape( (/&
         0.7831769E-02, 0.9179345E-02, 0.1141457E-01, 0.1422596E-01, &
         0.1779609E-01, 0.2240106E-01, 0.2849078E-01, 0.3686021E-01, &
         0.4902377E-01, 0.6804754E-01, 0.1003985E+00, 0.1598312E+00, &
         0.2738129E+00, 0.4842010E+00, 0.8094147E+00, 0.1232453E+01, &
         0.1780571E+01, 0.2230534E+01, 0.2762867E+01, 0.3204298E+01, &
         0.3531933E+01, 0.3553329E+01, 0.2890255E+01, 0.2053158E+01, &
         0.2194848E+01, 0.2679376E+01, 0.2134710E+01, 0.2353978E+01, &
         0.2197762E+01, 0.2211876E+01, 0.7078757E-02, 0.8293798E-02, &
         0.1030588E-01, 0.1282875E-01, 0.1601650E-01, 0.2009444E-01, &
         0.2541709E-01, 0.3258747E-01, 0.4271624E-01, 0.5799330E-01, &
         0.8296370E-01, 0.1273095E+00, 0.2112185E+00, 0.3711553E+00, &
         0.6447718E+00, 0.1018571E+01, 0.1524037E+01, 0.2011806E+01, &
         0.2550897E+01, 0.3026439E+01, 0.3453130E+01, 0.3582561E+01, &
         0.3267390E+01, 0.2385726E+01, 0.1939468E+01, 0.2592190E+01, &
         0.2432119E+01, 0.2266550E+01, 0.2245551E+01, 0.2216399E+01, &
         0.5893333E-02, 0.6901805E-02, 0.8567906E-02, 0.1064833E-01, &
         0.1325984E-01, 0.1656435E-01, 0.2080209E-01, 0.2635475E-01, &
         0.3388031E-01, 0.4460318E-01, 0.6096004E-01, 0.8803720E-01, &
         0.1366774E+00, 0.2292520E+00, 0.4044449E+00, 0.6954624E+00, &
         0.1083038E+01, 0.1609830E+01, 0.2073517E+01, 0.2627049E+01, &
         0.3090688E+01, 0.3465830E+01, 0.3572824E+01, 0.3151445E+01, &
         0.2316983E+01, 0.2038591E+01, 0.2610828E+01, 0.2285346E+01, &
         0.2321426E+01, 0.2202826E+01, 0.5025937E-02, 0.5884357E-02, &
         0.7300625E-02, 0.9065186E-02, 0.1127127E-01, 0.1404497E-01, &
         0.1756500E-01, 0.2210026E-01, 0.2808706E-01, 0.3629318E-01, &
         0.4817502E-01, 0.6667232E-01, 0.9797529E-01, 0.1552611E+00, &
         0.2650325E+00, 0.4687851E+00, 0.7881204E+00, 0.1203845E+01, &
         0.1751305E+01, 0.2198479E+01, 0.2739005E+01, 0.3182265E+01, &
         0.3504886E+01, 0.3572367E+01, 0.2939144E+01, 0.2093931E+01, &
         0.2184590E+01, 0.2653599E+01, 0.2120970E+01, 0.2362548E+01, &
         0.4433422E-02, 0.5189727E-02, 0.6437071E-02, 0.7988810E-02, &
         0.9924201E-02, 0.1234918E-01, 0.1540817E-01, 0.1931125E-01, &
         0.2438432E-01, 0.3117430E-01, 0.4067514E-01, 0.5482540E-01, &
         0.7761944E-01, 0.1175483E+00, 0.1924794E+00, 0.3360287E+00, &
         0.5890094E+00, 0.9476252E+00, 0.1424953E+01, 0.1942187E+01, &
         0.2450353E+01, 0.2937539E+01, 0.3395996E+01, 0.3633851E+01, &
         0.3374526E+01, 0.2582175E+01, 0.1926354E+01, 0.2447965E+01, &
         0.2536666E+01, 0.2114620E+01, 0.3825215E-02, 0.4476689E-02, &
         0.5550712E-02, 0.6885700E-02, 0.8548012E-02, 0.1062361E-01, &
         0.1322856E-01, 0.1652461E-01, 0.2075062E-01, 0.2628650E-01, &
         0.3378580E-01, 0.4446473E-01, 0.6074107E-01, 0.8766055E-01, &
         0.1359789E+00, 0.2279065E+00, 0.4019802E+00, 0.6917781E+00, &
         0.1078330E+01, 0.1603768E+01, 0.2068956E+01, 0.2621971E+01, &
         0.3086536E+01, 0.3465486E+01, 0.3571346E+01, 0.3157148E+01, &
         0.2325003E+01, 0.2047926E+01, 0.2618621E+01, 0.2334253E+01/),&
         (/30,1,6/))

    rdqext(1:30,1:1,37:42)=reshape( (/&
         0.3549440E-02, 0.4153976E-02, 0.5149985E-02, 0.6387233E-02, &
         0.7926833E-02, 0.9847006E-02, 0.1225196E-01, 0.1528491E-01, &
         0.1915287E-01, 0.2417619E-01, 0.3089087E-01, 0.4026837E-01, &
         0.5419936E-01, 0.7657288E-01, 0.1156511E+00, 0.1888472E+00, &
         0.3291641E+00, 0.5778280E+00, 0.9332777E+00, 0.1404682E+01, &
         0.1927353E+01, 0.2428208E+01, 0.2918166E+01, 0.3377411E+01, &
         0.3638303E+01, 0.3399381E+01, 0.2617698E+01, 0.1981718E+01, &
         0.2419561E+01, 0.2593332E+01, 0.3034851E-02, 0.3551537E-02, &
         0.4402093E-02, 0.5458310E-02, 0.6770817E-02, 0.8404701E-02, &
         0.1044440E-01, 0.1300307E-01, 0.1623797E-01, 0.2038019E-01, &
         0.2579521E-01, 0.3310764E-01, 0.4347322E-01, 0.5917915E-01, &
         0.8498399E-01, 0.1310288E+00, 0.2183748E+00, 0.3844299E+00, &
         0.6652285E+00, 0.1044529E+01, 0.1559250E+01, 0.2036551E+01, &
         0.2583322E+01, 0.3054318E+01, 0.3461036E+01, 0.3569443E+01, &
         0.3209359E+01, 0.2341575E+01, 0.2014180E+01, 0.2588365E+01, &
         0.2387623E-02, 0.2794191E-02, 0.3462401E-02, 0.4291734E-02, &
         0.5320973E-02, 0.6599757E-02, 0.8191557E-02, 0.1017778E-01, &
         0.1266782E-01, 0.1581217E-01, 0.1983106E-01, 0.2506921E-01, &
         0.3211028E-01, 0.4202449E-01, 0.5691476E-01, 0.8113549E-01, &
         0.1239570E+00, 0.2047747E+00, 0.3591337E+00, 0.6259582E+00, &
         0.9947014E+00, 0.1491055E+01, 0.1988876E+01, 0.2518942E+01, &
         0.2998372E+01, 0.3440423E+01, 0.3601781E+01, 0.3304095E+01, &
         0.2485376E+01, 0.2048215E+01, 0.2057811E-02, 0.2408417E-02, &
         0.2984623E-02, 0.3698977E-02, 0.4585142E-02, 0.5685347E-02, &
         0.7053476E-02, 0.8756823E-02, 0.1088493E-01, 0.1355791E-01, &
         0.1694377E-01, 0.2129343E-01, 0.2700860E-01, 0.3478724E-01, &
         0.4593816E-01, 0.6308043E-01, 0.9170132E-01, 0.1434965E+00, &
         0.2423881E+00, 0.4283417E+00, 0.7306451E+00, 0.1128293E+01, &
         0.1666130E+01, 0.2118356E+01, 0.2672537E+01, 0.3127278E+01, &
         0.3469966E+01, 0.3582943E+01, 0.3093040E+01, 0.2194497E+01, &
         0.1826921E-02, 0.2137052E-02, 0.2648383E-02, 0.3281895E-02, &
         0.4067792E-02, 0.5043134E-02, 0.6254602E-02, 0.7761572E-02, &
         0.9640561E-02, 0.1199300E-01, 0.1495711E-01, 0.1873215E-01, &
         0.2362413E-01, 0.3014096E-01, 0.3919665E-01, 0.5255830E-01, &
         0.7384456E-01, 0.1107295E+00, 0.1794446E+00, 0.3113199E+00, &
         0.5483193E+00, 0.8950609E+00, 0.1350772E+01, 0.1885957E+01, &
         0.2367349E+01, 0.2866825E+01, 0.3320561E+01, 0.3631001E+01, &
         0.3433996E+01, 0.2677806E+01, 0.1663156E-02, 0.1945878E-02, &
         0.2411773E-02, 0.2989146E-02, 0.3704706E-02, 0.4592277E-02, &
         0.5694410E-02, 0.7064310E-02, 0.8770572E-02, 0.1090218E-01, &
         0.1357968E-01, 0.1697154E-01, 0.2132941E-01, 0.2705653E-01, &
         0.3485396E-01, 0.4603667E-01, 0.6323758E-01, 0.9197413E-01, &
         0.1440057E+00, 0.2433691E+00, 0.4301138E+00, 0.7332162E+00, &
         0.1131626E+01, 0.1670127E+01, 0.2121743E+01, 0.2675674E+01, &
         0.3129785E+01, 0.3470570E+01, 0.3583046E+01, 0.3089532E+01/),&
         (/30,1,6/))

    rdqext(1:30,1:1,43:48)=reshape( (/&
         0.1567061E-02, 0.1833278E-02, 0.2271720E-02, 0.2815047E-02, &
         0.3488794E-02, 0.4324771E-02, 0.5361873E-02, 0.6650588E-02, &
         0.8254934E-02, 0.1025699E-01, 0.1276743E-01, 0.1593864E-01, &
         0.1999407E-01, 0.2528441E-01, 0.3240535E-01, 0.4245194E-01, &
         0.5758059E-01, 0.8226308E-01, 0.1260231E+00, 0.2087452E+00, &
         0.3665486E+00, 0.6375951E+00, 0.1009468E+01, 0.1511518E+01, &
         0.2003096E+01, 0.2538950E+01, 0.3016001E+01, 0.3448988E+01, &
         0.3589492E+01, 0.3284389E+01, 0.1488439E-02, 0.1741797E-02, &
         0.2157624E-02, 0.2673772E-02, 0.3314134E-02, 0.4107578E-02, &
         0.5092440E-02, 0.6315934E-02, 0.7837981E-02, 0.9735893E-02, &
         0.1211273E-01, 0.1510853E-01, 0.1892637E-01, 0.2387865E-01, &
         0.3048642E-01, 0.3968957E-01, 0.5331150E-01, 0.7509397E-01, &
         0.1129788E+00, 0.1837380E+00, 0.3194806E+00, 0.5618942E+00, &
         0.9127132E+00, 0.1375636E+01, 0.1905448E+01, 0.2395723E+01, &
         0.2890356E+01, 0.3347702E+01, 0.3638081E+01, 0.3419334E+01, &
         0.1409372E-02, 0.1649454E-02, 0.2044782E-02, 0.2532639E-02, &
         0.3138920E-02, 0.3890525E-02, 0.4823158E-02, 0.5981331E-02, &
         0.7421338E-02, 0.9215627E-02, 0.1146004E-01, 0.1428342E-01, &
         0.1786957E-01, 0.2249676E-01, 0.2861935E-01, 0.3704118E-01, &
         0.4929540E-01, 0.6848914E-01, 0.1011789E+00, 0.1613057E+00, &
         0.2766435E+00, 0.4891382E+00, 0.8161697E+00, 0.1241602E+01, &
         0.1789592E+01, 0.2240942E+01, 0.2770504E+01, 0.3211761E+01, &
         0.3541587E+01, 0.3546342E+01, 0.1330705E-02, 0.1558003E-02, &
         0.1930830E-02, 0.2391852E-02, 0.2964355E-02, 0.3673612E-02, &
         0.4553958E-02, 0.5646588E-02, 0.7004765E-02, 0.8696110E-02, &
         0.1080916E-01, 0.1346239E-01, 0.1682217E-01, 0.2113581E-01, &
         0.2679852E-01, 0.3449527E-01, 0.4550745E-01, 0.6239447E-01, &
         0.9051254E-01, 0.1412796E+00, 0.2381176E+00, 0.4206070E+00, &
         0.7193625E+00, 0.1113713E+01, 0.1648399E+01, 0.2103685E+01, &
         0.2658494E+01, 0.3116058E+01, 0.3468100E+01, 0.3581523E+01, &
         0.1251933E-02, 0.1466046E-02, 0.1816172E-02, 0.2251175E-02, &
         0.2788802E-02, 0.3456600E-02, 0.4284468E-02, 0.5312049E-02, &
         0.6588972E-02, 0.8177907E-02, 0.1016081E-01, 0.1264632E-01, &
         0.1578495E-01, 0.1979593E-01, 0.2502291E-01, 0.3204691E-01, &
         0.4193268E-01, 0.5677203E-01, 0.8089415E-01, 0.1235155E+00, &
         0.2039266E+00, 0.3575470E+00, 0.6234539E+00, 0.9915216E+00, &
         0.1486626E+01, 0.1985795E+01, 0.2514534E+01, 0.2994468E+01, &
         0.3438259E+01, 0.3604452E+01, 0.1174132E-02, 0.1373839E-02, &
         0.1703073E-02, 0.2110218E-02, 0.2614515E-02, 0.3239858E-02, &
         0.4015790E-02, 0.4978259E-02, 0.6173875E-02, 0.7660876E-02, &
         0.9514965E-02, 0.1183530E-01, 0.1475762E-01, 0.1847636E-01, &
         0.2328914E-01, 0.2968762E-01, 0.3855173E-01, 0.5157678E-01, &
         0.7222367E-01, 0.1078229E+00, 0.1739076E+00, 0.3007662E+00, &
         0.5305643E+00, 0.8717529E+00, 0.1318129E+01, 0.1859136E+01, &
         0.2329573E+01, 0.2836743E+01, 0.3284825E+01, 0.3613000E+01/),&
         (/30,1,6/))

    rdqext(1:30,1:1,49:50)=reshape( (/&
         0.1095844E-02, 0.1281458E-02, 0.1588486E-02, 0.1968793E-02, &
         0.2439237E-02, 0.3023065E-02, 0.3747140E-02, 0.4644954E-02, &
         0.5759797E-02, 0.7145516E-02, 0.8871601E-02, 0.1102877E-01, &
         0.1373922E-01, 0.1717497E-01, 0.2159321E-01, 0.2740864E-01, &
         0.3534436E-01, 0.4676276E-01, 0.6439901E-01, 0.9399497E-01, &
         0.1477852E+00, 0.2506478E+00, 0.4432062E+00, 0.7520572E+00, &
         0.1156173E+01, 0.1698893E+01, 0.2147104E+01, 0.2698069E+01, &
         0.3147782E+01, 0.3476962E+01, 0.9974022E-03, 0.1165836E-02, &
         0.1443403E-02, 0.1788278E-02, 0.2216780E-02, 0.2746567E-02, &
         0.3404155E-02, 0.4219844E-02, 0.5231818E-02, 0.6489299E-02, &
         0.8053846E-02, 0.1000574E-01, 0.1245149E-01, 0.1553765E-01, &
         0.1947781E-01, 0.2460348E-01, 0.3147328E-01, 0.4110512E-01, &
         0.5548911E-01, 0.7873266E-01, 0.1195716E+00, 0.1963573E+00, &
         0.3433384E+00, 0.6008126E+00, 0.9627100E+00, 0.1446232E+01, &
         0.1957438E+01, 0.2473066E+01, 0.2957629E+01, 0.3413126E+01/),&
         (/30,1,2/))

    qscat(1:30,1:1,1:6)=reshape( (/&
         0.1535112E-05, 0.8998886E-09, 0.1729664E-08, 0.4070829E-08, &
         0.9589681E-08, 0.2261375E-07, 0.5328723E-07, 0.1255802E-06, &
         0.2960315E-06, 0.6976022E-06, 0.1644110E-05, 0.3875383E-05, &
         0.9135672E-05, 0.2154004E-04, 0.5080069E-04, 0.1198578E-03, &
         0.2829565E-03, 0.6685623E-03, 0.1581542E-02, 0.3746952E-02, &
         0.8891535E-02, 0.2111296E-01, 0.4993847E-01, 0.1159231E+00, &
         0.2536170E+00, 0.4850497E+00, 0.7988749E+00, 0.1241824E+01, &
         0.1650353E+01, 0.2180239E+01, 0.2919784E-05, 0.5829020E-09, &
         0.1120681E-08, 0.2640323E-08, 0.6234001E-08, 0.1468913E-07, &
         0.3459422E-07, 0.8157038E-07, 0.1922175E-06, 0.4529455E-06, &
         0.1067609E-05, 0.2516517E-05, 0.5931983E-05, 0.1398515E-04, &
         0.3297836E-04, 0.7779063E-04, 0.1835863E-03, 0.4335718E-03, &
         0.1024991E-02, 0.2426435E-02, 0.5753419E-02, 0.1366061E-01, &
         0.3240611E-01, 0.7616633E-01, 0.1729649E+00, 0.3576158E+00, &
         0.6282287E+00, 0.1009327E+01, 0.1438690E+01, 0.1939538E+01, &
         0.5553430E-05, 0.3008072E-09, 0.5785822E-09, 0.1359781E-08, &
         0.3207448E-08, 0.7555678E-08, 0.1780665E-07, 0.4196275E-07, &
         0.9891541E-07, 0.2330545E-06, 0.5492852E-06, 0.1294664E-05, &
         0.3051624E-05, 0.7193325E-05, 0.1696003E-04, 0.3999544E-04, &
         0.9435236E-04, 0.2227023E-03, 0.5260573E-03, 0.1243972E-02, &
         0.2945860E-02, 0.6987582E-02, 0.1659290E-01, 0.3932323E-01, &
         0.9199918E-01, 0.2059242E+00, 0.4121182E+00, 0.7000768E+00, &
         0.1114519E+01, 0.1523778E+01, 0.1056262E-04, 0.1667335E-09, &
         0.3215872E-09, 0.7584442E-09, 0.1788107E-08, 0.4214065E-08, &
         0.9933288E-08, 0.2340783E-07, 0.5519188E-07, 0.1300859E-06, &
         0.3066088E-06, 0.7225132E-06, 0.1702958E-05, 0.4013857E-05, &
         0.9462706E-05, 0.2231118E-04, 0.5262020E-04, 0.1241540E-03, &
         0.2931054E-03, 0.6925743E-03, 0.1638427E-02, 0.3881993E-02, &
         0.9212516E-02, 0.2187444E-01, 0.5172065E-01, 0.1198937E+00, &
         0.2612906E+00, 0.4962334E+00, 0.8145505E+00, 0.1259660E+01, &
         0.2009011E-04, 0.9312614E-10, 0.1788817E-09, 0.4216057E-09, &
         0.9970077E-09, 0.2345042E-08, 0.5529286E-08, 0.1303760E-07, &
         0.3071025E-07, 0.7236788E-07, 0.1705414E-06, 0.4019309E-06, &
         0.9472516E-06, 0.2232579E-05, 0.5262621E-05, 0.1240705E-04, &
         0.2925507E-04, 0.6900538E-04, 0.1628382E-03, 0.3845306E-03, &
         0.9089075E-03, 0.2151185E-02, 0.5099573E-02, 0.1210664E-01, &
         0.2873163E-01, 0.6768199E-01, 0.1548432E+00, 0.3259844E+00, &
         0.5861368E+00, 0.9459786E+00, 0.3821139E-04, 0.4500313E-10, &
         0.8613462E-10, 0.2036637E-09, 0.4809451E-09, 0.1131263E-08, &
         0.2664907E-08, 0.6288710E-08, 0.1480967E-07, 0.3488467E-07, &
         0.8223534E-07, 0.1937976E-06, 0.4568121E-06, 0.1076551E-05, &
         0.2537450E-05, 0.5981515E-05, 0.1410117E-04, 0.3325193E-04, &
         0.7843683E-04, 0.1851134E-03, 0.4371810E-03, 0.1033536E-02, &
         0.2446699E-02, 0.5801552E-02, 0.1377500E-01, 0.3267637E-01, &
         0.7678847E-01, 0.1742809E+00, 0.3598650E+00, 0.6311989E+00/),&
         (/30,1,6/))

    qscat(1:30,1:1,7:12)=reshape( (/&
         0.7267806E-04, 0.1765787E-10, 0.3378414E-10, 0.7984936E-10, &
         0.1884185E-09, 0.4455544E-09, 0.1046103E-08, 0.2469826E-08, &
         0.5812942E-08, 0.1370527E-07, 0.3230256E-07, 0.7611979E-07, &
         0.1794039E-06, 0.4228477E-06, 0.9965112E-06, 0.2348792E-05, &
         0.5536366E-05, 0.1305259E-04, 0.3077852E-04, 0.7260025E-04, &
         0.1713272E-03, 0.4045929E-03, 0.9563949E-03, 0.2263770E-02, &
         0.5366995E-02, 0.1274226E-01, 0.3023528E-01, 0.7115984E-01, &
         0.1623103E+00, 0.3391695E+00, 0.1382336E-03, 0.5461978E-11, &
         0.1045634E-10, 0.2480865E-10, 0.5877250E-10, 0.1386328E-09, &
         0.3260322E-09, 0.7689554E-09, 0.1809900E-08, 0.4257048E-08, &
         0.1004131E-07, 0.2366112E-07, 0.5576113E-07, 0.1314203E-06, &
         0.3096752E-06, 0.7299032E-06, 0.1720274E-05, 0.4054981E-05, &
         0.9559084E-05, 0.2253919E-04, 0.5315715E-04, 0.1254218E-03, &
         0.2961037E-03, 0.6996637E-03, 0.1655233E-02, 0.3921884E-02, &
         0.9307324E-02, 0.2209932E-01, 0.5224660E-01, 0.1210630E+00, &
         0.2629204E-03, 0.8523221E-12, 0.1655710E-11, 0.3974937E-11, &
         0.9475196E-11, 0.2160493E-10, 0.5173064E-10, 0.1217817E-09, &
         0.2857532E-09, 0.6785524E-09, 0.1597406E-08, 0.3767277E-08, &
         0.8871726E-08, 0.2090630E-07, 0.4927249E-07, 0.1161264E-06, &
         0.2736318E-06, 0.6449878E-06, 0.1520047E-05, 0.3583024E-05, &
         0.8446612E-05, 0.1991419E-04, 0.4696518E-04, 0.1108043E-03, &
         0.2615690E-03, 0.6179706E-03, 0.1461678E-02, 0.3462458E-02, &
         0.8215324E-02, 0.1950821E-01, 0.5000746E-03, 0.2305856E-02, &
         0.5466960E-02, 0.1297981E-01, 0.3079700E-01, 0.7245710E-01, &
         0.1650817E+00, 0.3440087E+00, 0.6102081E+00, 0.9822114E+00, &
         0.1416896E+01, 0.1905245E+01, 0.2396360E+01, 0.2865409E+01, &
         0.3078600E+01, 0.2821023E+01, 0.2079308E+01, 0.1451826E+01, &
         0.2001947E+01, 0.2059541E+01, 0.1667592E+01, 0.1893693E+01, &
         0.1813272E+01, 0.1784861E+01, 0.1766699E+01, 0.1802252E+01, &
         0.1777354E+01, 0.1773361E+01, 0.1769136E+01, 0.1735761E+01, &
         0.8336438E-03, 0.1607949E-02, 0.3809640E-02, 0.9040533E-02, &
         0.2146645E-01, 0.5076601E-01, 0.1177685E+00, 0.2571919E+00, &
         0.4902766E+00, 0.8061812E+00, 0.1250220E+01, 0.1660402E+01, &
         0.2188359E+01, 0.2646668E+01, 0.2970078E+01, 0.3017053E+01, &
         0.2416580E+01, 0.1578704E+01, 0.1703342E+01, 0.2210760E+01, &
         0.1668151E+01, 0.1938131E+01, 0.1788607E+01, 0.1854164E+01, &
         0.1828766E+01, 0.1805977E+01, 0.1766851E+01, 0.1782201E+01, &
         0.1754907E+01, 0.1744322E+01, 0.6669299E-03, 0.1286015E-02, &
         0.3045615E-02, 0.7224638E-02, 0.1715594E-01, 0.4064883E-01, &
         0.9501187E-01, 0.2120690E+00, 0.4218818E+00, 0.7130325E+00, &
         0.1132582E+01, 0.1539426E+01, 0.2076132E+01, 0.2551370E+01, &
         0.2925028E+01, 0.3023993E+01, 0.2621071E+01, 0.1807222E+01, &
         0.1550507E+01, 0.2171909E+01, 0.1841270E+01, 0.1927457E+01, &
         0.1817936E+01, 0.1826131E+01, 0.1748519E+01, 0.1781892E+01, &
         0.1787518E+01, 0.1738967E+01, 0.1760810E+01, 0.1754529E+01/),&
         (/30,1,6/))

    qscat(1:30,1:1,13:18)=reshape( (/&
         0.4939059E-03, 0.9520435E-03, 0.2253458E-02, 0.5342496E-02, &
         0.1268401E-01, 0.3009751E-01, 0.7084160E-01, 0.1616292E+00, &
         0.3379757E+00, 0.6021811E+00, 0.9701180E+00, 0.1407054E+01, &
         0.1889471E+01, 0.2381627E+01, 0.2852786E+01, 0.3083985E+01, &
         0.2840933E+01, 0.2098689E+01, 0.1445787E+01, 0.2000638E+01, &
         0.2100897E+01, 0.1701578E+01, 0.1950433E+01, 0.1846692E+01, &
         0.1828247E+01, 0.1768301E+01, 0.1749795E+01, 0.1776199E+01, &
         0.1766105E+01, 0.1769362E+01, 0.4058536E-03, 0.7821483E-03, &
         0.1850704E-02, 0.4385977E-02, 0.1041047E-01, 0.2471486E-01, &
         0.5835165E-01, 0.1345477E+00, 0.2890299E+00, 0.5356299E+00, &
         0.8710567E+00, 0.1318777E+01, 0.1752060E+01, 0.2261234E+01, &
         0.2721528E+01, 0.3050559E+01, 0.2949300E+01, 0.2250151E+01, &
         0.1530140E+01, 0.1829173E+01, 0.2215410E+01, 0.1710420E+01, &
         0.1978602E+01, 0.1811628E+01, 0.1824998E+01, 0.1844743E+01, &
         0.1773051E+01, 0.1760529E+01, 0.1747295E+01, 0.1756304E+01, &
         0.2877247E-03, 0.5543064E-03, 0.1310881E-02, 0.3104613E-02, &
         0.7364837E-02, 0.1748890E-01, 0.4143228E-01, 0.9678915E-01, &
         0.2156751E+00, 0.4275572E+00, 0.7205915E+00, 0.1142951E+01, &
         0.1548665E+01, 0.2086365E+01, 0.2559992E+01, 0.2925613E+01, &
         0.3027108E+01, 0.2608634E+01, 0.1765504E+01, 0.1553489E+01, &
         0.2209895E+01, 0.1835341E+01, 0.1919490E+01, 0.1788675E+01, &
         0.1842055E+01, 0.1772597E+01, 0.1804116E+01, 0.1782525E+01, &
         0.1742764E+01, 0.1744599E+01, 0.1639157E-03, 0.3156405E-03, &
         0.7458854E-03, 0.1764762E-02, 0.4181915E-02, 0.9925411E-02, &
         0.2356508E-01, 0.5567071E-01, 0.1286458E+00, 0.2779670E+00, &
         0.5200978E+00, 0.8485439E+00, 0.1296149E+01, 0.1719943E+01, &
         0.2235500E+01, 0.2693254E+01, 0.3024626E+01, 0.2977588E+01, &
         0.2311464E+01, 0.1667730E+01, 0.1821186E+01, 0.2224497E+01, &
         0.1704216E+01, 0.1982319E+01, 0.1825384E+01, 0.1832962E+01, &
         0.1851820E+01, 0.1797334E+01, 0.1770782E+01, 0.1758385E+01, &
         0.1360369E-03, 0.2619239E-03, 0.6188128E-03, 0.1463670E-02, &
         0.3467186E-02, 0.8226552E-02, 0.1953486E-01, 0.4623934E-01, &
         0.1076391E+00, 0.2373873E+00, 0.4609284E+00, 0.7655962E+00, &
         0.1201772E+01, 0.1605729E+01, 0.2142393E+01, 0.2606626E+01, &
         0.2936520E+01, 0.3031624E+01, 0.2534039E+01, 0.1769609E+01, &
         0.1651198E+01, 0.2273975E+01, 0.1833097E+01, 0.1972575E+01, &
         0.1789315E+01, 0.1821071E+01, 0.1794554E+01, 0.1826451E+01, &
         0.1784653E+01, 0.1767444E+01, 0.1044451E-03, 0.2010619E-03, &
         0.4748902E-03, 0.1122807E-02, 0.2658433E-02, 0.6304602E-02, &
         0.1497037E-01, 0.3549854E-01, 0.8326766E-01, 0.1878861E+00, &
         0.3827521E+00, 0.6613459E+00, 0.1058644E+01, 0.1477897E+01, &
         0.1997881E+01, 0.2482627E+01, 0.2914294E+01, 0.3026468E+01, &
         0.2727181E+01, 0.1884788E+01, 0.1458210E+01, 0.2130828E+01, &
         0.1948996E+01, 0.1807852E+01, 0.1839294E+01, 0.1855008E+01, &
         0.1810482E+01, 0.1750667E+01, 0.1797086E+01, 0.1775596E+01/),&
         (/30,1,6/))

    qscat(1:30,1:1,19:24)=reshape( (/&
         0.9441071E-04, 0.1817370E-03, 0.4291996E-03, 0.1014639E-02, &
         0.2401892E-02, 0.5695102E-02, 0.1352202E-01, 0.3207864E-01, &
         0.7541227E-01, 0.1713672E+00, 0.3548766E+00, 0.6246089E+00, &
         0.1003889E+01, 0.1434349E+01, 0.1932783E+01, 0.2422200E+01, &
         0.2884667E+01, 0.3065633E+01, 0.2791466E+01, 0.2049569E+01, &
         0.1575594E+01, 0.2073899E+01, 0.2041062E+01, 0.1722158E+01, &
         0.1899726E+01, 0.1846775E+01, 0.1793869E+01, 0.1778048E+01, &
         0.1770816E+01, 0.1789440E+01, 0.8151573E-04, 0.1569007E-03, &
         0.3704934E-03, 0.8756850E-03, 0.2072421E-02, 0.4912498E-02, &
         0.1166200E-01, 0.2767925E-01, 0.6524293E-01, 0.1495746E+00, &
         0.3165503E+00, 0.5734129E+00, 0.9269114E+00, 0.1370421E+01, &
         0.1830986E+01, 0.2328271E+01, 0.2798606E+01, 0.3087561E+01, &
         0.2882272E+01, 0.2176012E+01, 0.1618124E+01, 0.1939663E+01, &
         0.2155108E+01, 0.1668577E+01, 0.1943499E+01, 0.1805721E+01, &
         0.1826187E+01, 0.1817263E+01, 0.1757304E+01, 0.1767201E+01, &
         0.5565022E-04, 0.1070932E-03, 0.2528019E-03, 0.5972348E-03, &
         0.1412555E-02, 0.3345877E-02, 0.7938247E-02, 0.1885045E-01, &
         0.4463263E-01, 0.1040232E+00, 0.2302090E+00, 0.4500423E+00, &
         0.7507968E+00, 0.1183003E+01, 0.1586547E+01, 0.2124766E+01, &
         0.2591938E+01, 0.2930691E+01, 0.3031670E+01, 0.2565965E+01, &
         0.1732758E+01, 0.1644230E+01, 0.2178035E+01, 0.1792361E+01, &
         0.1880221E+01, 0.1821903E+01, 0.1829716E+01, 0.1787137E+01, &
         0.1815011E+01, 0.1770356E+01, 0.3980661E-04, 0.7659413E-04, &
         0.1807602E-03, 0.4268896E-03, 0.1009173E-02, 0.2388930E-02, &
         0.5664303E-02, 0.1344884E-01, 0.3190569E-01, 0.7501386E-01, &
         0.1705220E+00, 0.3534238E+00, 0.6226878E+00, 0.1001000E+01, &
         0.1432039E+01, 0.1929170E+01, 0.2418807E+01, 0.2882358E+01, &
         0.3067569E+01, 0.2794610E+01, 0.2053816E+01, 0.1557362E+01, &
         0.2052078E+01, 0.2018643E+01, 0.1696542E+01, 0.1886784E+01, &
         0.1840827E+01, 0.1781331E+01, 0.1751741E+01, 0.1756093E+01, &
         0.2681827E-04, 0.5159516E-04, 0.1217334E-03, 0.2873892E-03, &
         0.6790494E-03, 0.1606384E-02, 0.3805922E-02, 0.9031690E-02, &
         0.2144548E-01, 0.5071694E-01, 0.1176591E+00, 0.2569804E+00, &
         0.4899683E+00, 0.8057491E+00, 0.1249726E+01, 0.1659806E+01, &
         0.2187880E+01, 0.2646229E+01, 0.2969582E+01, 0.3017386E+01, &
         0.2417457E+01, 0.1579919E+01, 0.1704027E+01, 0.2208308E+01, &
         0.1668819E+01, 0.1938163E+01, 0.1787129E+01, 0.1846662E+01, &
         0.1827163E+01, 0.1803410E+01, 0.2525160E-04, 0.4857985E-04, &
         0.1146153E-03, 0.2705707E-03, 0.6392652E-03, 0.1512128E-02, &
         0.3582197E-02, 0.8499914E-02, 0.2018366E-01, 0.4776109E-01, &
         0.1110539E+00, 0.2441137E+00, 0.4710058E+00, 0.7794141E+00, &
         0.1218766E+01, 0.1624015E+01, 0.2158341E+01, 0.2620128E+01, &
         0.2944883E+01, 0.3030327E+01, 0.2490729E+01, 0.1736488E+01, &
         0.1670145E+01, 0.2207124E+01, 0.1796909E+01, 0.2016083E+01, &
         0.1831616E+01, 0.1820178E+01, 0.1820695E+01, 0.1821183E+01/),&
         (/30,1,6/))

    qscat(1:30,1:1,25:30)=reshape( (/&
         0.2244690E-04, 0.4318224E-04, 0.1018737E-03, 0.2404698E-03, &
         0.5680734E-03, 0.1343483E-02, 0.3181970E-02, 0.7548691E-02, &
         0.1792551E-01, 0.4245912E-01, 0.9911473E-01, 0.2203728E+00, &
         0.4348921E+00, 0.7303964E+00, 0.1156200E+01, 0.1560783E+01, &
         0.2099243E+01, 0.2570755E+01, 0.2926596E+01, 0.3029954E+01, &
         0.2592838E+01, 0.1722351E+01, 0.1611468E+01, 0.2245607E+01, &
         0.1816836E+01, 0.1897303E+01, 0.1807774E+01, 0.1827002E+01, &
         0.1780776E+01, 0.1792670E+01, 0.1508714E-04, 0.2902153E-04, &
         0.6845398E-04, 0.1615359E-03, 0.3814493E-03, 0.9016190E-03, &
         0.2133902E-02, 0.5058516E-02, 0.1200905E-01, 0.2850071E-01, &
         0.6714708E-01, 0.1536901E+00, 0.3239290E+00, 0.5833729E+00, &
         0.9418291E+00, 0.1383321E+01, 0.1851519E+01, 0.2346708E+01, &
         0.2818617E+01, 0.3089535E+01, 0.2873804E+01, 0.2160190E+01, &
         0.1588495E+01, 0.1916183E+01, 0.2115090E+01, 0.1767522E+01, &
         0.1952622E+01, 0.1827966E+01, 0.1801739E+01, 0.1806674E+01, &
         0.1903640E-04, 0.3661983E-04, 0.8638486E-04, 0.2038823E-03, &
         0.4815583E-03, 0.1138599E-02, 0.2695894E-02, 0.6393607E-02, &
         0.1518182E-01, 0.3599739E-01, 0.8440971E-01, 0.1902647E+00, &
         0.3866871E+00, 0.6665245E+00, 0.1066252E+01, 0.1483992E+01, &
         0.2006419E+01, 0.2490376E+01, 0.2916542E+01, 0.3022475E+01, &
         0.2712411E+01, 0.1870931E+01, 0.1473295E+01, 0.2130482E+01, &
         0.1957074E+01, 0.1783501E+01, 0.1843272E+01, 0.1833044E+01, &
         0.1786526E+01, 0.1781775E+01, 0.1164130E-04, 0.2239145E-04, &
         0.5281003E-04, 0.1246017E-03, 0.2941640E-03, 0.6950784E-03, &
         0.1644360E-02, 0.3896082E-02, 0.9245978E-02, 0.2195383E-01, &
         0.5190634E-01, 0.1203065E+00, 0.2620847E+00, 0.4973831E+00, &
         0.8161713E+00, 0.1261464E+01, 0.1674268E+01, 0.2199425E+01, &
         0.2656993E+01, 0.2982174E+01, 0.3008611E+01, 0.2397170E+01, &
         0.1566271E+01, 0.1716662E+01, 0.2299051E+01, 0.1702490E+01, &
         0.1945370E+01, 0.1792890E+01, 0.1841396E+01, 0.1849848E+01, &
         0.9210702E-05, 0.1771512E-04, 0.4177722E-04, 0.9855754E-04, &
         0.2326356E-03, 0.5495443E-03, 0.1299600E-02, 0.3077847E-02, &
         0.7301237E-02, 0.1733785E-01, 0.4107690E-01, 0.9598324E-01, &
         0.2140418E+00, 0.4249914E+00, 0.7171717E+00, 0.1138275E+01, &
         0.1544475E+01, 0.2081769E+01, 0.2556129E+01, 0.2925340E+01, &
         0.3025771E+01, 0.2614344E+01, 0.1784989E+01, 0.1548554E+01, &
         0.2207918E+01, 0.1856935E+01, 0.1970010E+01, 0.1816692E+01, &
         0.1819728E+01, 0.1765595E+01, 0.6915328E-05, 0.1330029E-04, &
         0.3136207E-04, 0.7397675E-04, 0.1745795E-03, 0.4122800E-03, &
         0.9745910E-03, 0.2306919E-02, 0.5469492E-02, 0.1298585E-01, &
         0.3081126E-01, 0.7249000E-01, 0.1651520E+00, 0.3441308E+00, &
         0.6103705E+00, 0.9824560E+00, 0.1417094E+01, 0.1905561E+01, &
         0.2396655E+01, 0.2865654E+01, 0.3078475E+01, 0.2820626E+01, &
         0.2078964E+01, 0.1452464E+01, 0.2002201E+01, 0.2060558E+01, &
         0.1667076E+01, 0.1909250E+01, 0.1829093E+01, 0.1784407E+01/),&
         (/30,1,6/))

    qscat(1:30,1:1,31:36)=reshape( (/&
         0.4983547E-05, 0.9583965E-05, 0.2259792E-04, 0.5329637E-04, &
         0.1257496E-03, 0.2968784E-03, 0.7014993E-03, 0.1659575E-02, &
         0.3932194E-02, 0.9331833E-02, 0.2215747E-01, 0.5238259E-01, &
         0.1213650E+00, 0.2641168E+00, 0.5003189E+00, 0.8203176E+00, &
         0.1266057E+01, 0.1680060E+01, 0.2204016E+01, 0.2661379E+01, &
         0.2987473E+01, 0.3004859E+01, 0.2388839E+01, 0.1569596E+01, &
         0.1731562E+01, 0.2236392E+01, 0.1710583E+01, 0.1946052E+01, &
         0.1803016E+01, 0.1828850E+01, 0.3336525E-05, 0.6416181E-05, &
         0.1512720E-04, 0.3567183E-04, 0.8414762E-04, 0.1985995E-03, &
         0.4690658E-03, 0.1109019E-02, 0.2625725E-02, 0.6226893E-02, &
         0.1478571E-01, 0.3506286E-01, 0.8226940E-01, 0.1858021E+00, &
         0.3792889E+00, 0.6567879E+00, 0.1051920E+01, 0.1472531E+01, &
         0.1990234E+01, 0.2475642E+01, 0.2911975E+01, 0.3030586E+01, &
         0.2739404E+01, 0.1902393E+01, 0.1459937E+01, 0.2141851E+01, &
         0.1998446E+01, 0.1849737E+01, 0.1843773E+01, 0.1828411E+01, &
         0.1609733E-05, 0.3095464E-05, 0.7296939E-05, 0.1720301E-04, &
         0.4056917E-04, 0.9570643E-04, 0.2259008E-03, 0.5336171E-03, &
         0.1261881E-02, 0.2988344E-02, 0.7088546E-02, 0.1683272E-01, &
         0.3988795E-01, 0.9328347E-01, 0.2085486E+00, 0.4163024E+00, &
         0.7056225E+00, 0.1122294E+01, 0.1530449E+01, 0.2065823E+01, &
         0.2542608E+01, 0.2924463E+01, 0.3020677E+01, 0.2632277E+01, &
         0.1835126E+01, 0.1563670E+01, 0.2162810E+01, 0.1854607E+01, &
         0.1907334E+01, 0.1804355E+01, 0.8536749E-06, 0.1641531E-05, &
         0.3869212E-05, 0.9121159E-05, 0.2150566E-04, 0.5071913E-04, &
         0.1196660E-03, 0.2825053E-03, 0.6674958E-03, 0.1579006E-02, &
         0.3740940E-02, 0.8877230E-02, 0.2107904E-01, 0.4985903E-01, &
         0.1157458E+00, 0.2532728E+00, 0.4845448E+00, 0.7981708E+00, &
         0.1241009E+01, 0.1649390E+01, 0.2179454E+01, 0.2638595E+01, &
         0.2961328E+01, 0.3022736E+01, 0.2434684E+01, 0.1610582E+01, &
         0.1719000E+01, 0.2209665E+01, 0.1695634E+01, 0.1954261E+01, &
         0.5176838E-06, 0.9954088E-06, 0.2346362E-05, 0.5530849E-05, &
         0.1303859E-04, 0.3074600E-04, 0.7252304E-04, 0.1711466E-03, &
         0.4041660E-03, 0.9553821E-03, 0.2261380E-02, 0.5361313E-02, &
         0.1272874E-01, 0.3020329E-01, 0.7108603E-01, 0.1621523E+00, &
         0.3388926E+00, 0.6034029E+00, 0.9719585E+00, 0.1408524E+01, &
         0.1891890E+01, 0.2383881E+01, 0.2854789E+01, 0.3083266E+01, &
         0.2837948E+01, 0.2095346E+01, 0.1443578E+01, 0.1994793E+01, &
         0.2100053E+01, 0.1695583E+01, 0.2871854E-06, 0.5520937E-06, &
         0.1301202E-05, 0.3067068E-05, 0.7230215E-05, 0.1704622E-04, &
         0.4019847E-04, 0.9483170E-04, 0.2238343E-03, 0.5287331E-03, &
         0.1250309E-02, 0.2960899E-02, 0.7023325E-02, 0.1667778E-01, &
         0.3952316E-01, 0.9245398E-01, 0.2068544E+00, 0.4136037E+00, &
         0.7020445E+00, 0.1117285E+01, 0.1526141E+01, 0.2060742E+01, &
         0.2538260E+01, 0.2924151E+01, 0.3019181E+01, 0.2637346E+01, &
         0.1843114E+01, 0.1572628E+01, 0.2170437E+01, 0.1903235E+01/),&
         (/30,1,6/))

    qscat(1:30,1:1,37:42)=reshape( (/&
         0.2130258E-06, 0.4095630E-06, 0.9652463E-06, 0.2275018E-05, &
         0.5362792E-05, 0.1264280E-04, 0.2981126E-04, 0.7031697E-04, &
         0.1659372E-03, 0.3918542E-03, 0.9262424E-03, 0.2192283E-02, &
         0.5197192E-02, 0.1233868E-01, 0.2928065E-01, 0.6895282E-01, &
         0.1575780E+00, 0.3308383E+00, 0.5926470E+00, 0.9557645E+00, &
         0.1395105E+01, 0.1870389E+01, 0.2363955E+01, 0.2836186E+01, &
         0.3088190E+01, 0.2861268E+01, 0.2129971E+01, 0.1498628E+01, &
         0.1965536E+01, 0.2156023E+01, 0.1139642E-06, 0.2191180E-06, &
         0.5163641E-06, 0.1217040E-05, 0.2868713E-05, 0.6762461E-05, &
         0.1594324E-04, 0.3759716E-04, 0.8869241E-04, 0.2093319E-03, &
         0.4944418E-03, 0.1169102E-02, 0.2768247E-02, 0.6565538E-02, &
         0.1559031E-01, 0.3696075E-01, 0.8661259E-01, 0.1948363E+00, &
         0.3941952E+00, 0.6764117E+00, 0.1080671E+01, 0.1495647E+01, &
         0.2022275E+01, 0.2504614E+01, 0.2919853E+01, 0.3017293E+01, &
         0.2684913E+01, 0.1859088E+01, 0.1536330E+01, 0.2138908E+01, &
         0.4368971E-07, 0.8401571E-07, 0.1979443E-06, 0.4665041E-06, &
         0.1099423E-05, 0.2591267E-05, 0.6108387E-05, 0.1440077E-04, &
         0.3395877E-04, 0.8010484E-04, 0.1890504E-03, 0.4464910E-03, &
         0.1055577E-02, 0.2498972E-02, 0.5925733E-02, 0.1407011E-01, &
         0.3337342E-01, 0.7839159E-01, 0.1776647E+00, 0.3656191E+00, &
         0.6387883E+00, 0.1025151E+01, 0.1451272E+01, 0.1958829E+01, &
         0.2446588E+01, 0.2899261E+01, 0.3050108E+01, 0.2773031E+01, &
         0.2001010E+01, 0.1567242E+01, 0.2413002E-07, 0.4641096E-07, &
         0.1093743E-06, 0.2577491E-06, 0.6074086E-06, 0.1431569E-05, &
         0.3374563E-05, 0.7954735E-05, 0.1875470E-04, 0.4422979E-04, &
         0.1043463E-03, 0.2463109E-03, 0.5818865E-03, 0.1376202E-02, &
         0.3259612E-02, 0.7733199E-02, 0.1836365E-01, 0.4348901E-01, &
         0.1014429E+00, 0.2250517E+00, 0.4421329E+00, 0.7401199E+00, &
         0.1169107E+01, 0.1572968E+01, 0.2111608E+01, 0.2581027E+01, &
         0.2928094E+01, 0.3031281E+01, 0.2579799E+01, 0.1712282E+01, &
         0.1497814E-07, 0.2878358E-07, 0.6784342E-07, 0.1598607E-06, &
         0.3767407E-06, 0.8879427E-06, 0.2092868E-05, 0.4933270E-05, &
         0.1163001E-04, 0.2742337E-04, 0.6468265E-04, 0.1526333E-03, &
         0.3604091E-03, 0.8518165E-03, 0.2015840E-02, 0.4778121E-02, &
         0.1134260E-01, 0.2692300E-01, 0.6348782E-01, 0.1457673E+00, &
         0.3096642E+00, 0.5640594E+00, 0.9129618E+00, 0.1358052E+01, &
         0.1811523E+01, 0.2311157E+01, 0.2779208E+01, 0.3082438E+01, &
         0.2892156E+01, 0.2187006E+01, 0.1030113E-07, 0.1980520E-07, &
         0.4668839E-07, 0.1100485E-06, 0.2593487E-06, 0.6112119E-06, &
         0.1440549E-05, 0.3395485E-05, 0.8004393E-05, 0.1887181E-04, &
         0.4450613E-04, 0.1049992E-03, 0.2478533E-03, 0.5855319E-03, &
         0.1384842E-02, 0.3280106E-02, 0.7781909E-02, 0.1847929E-01, &
         0.4376074E-01, 0.1020566E+00, 0.2262809E+00, 0.4440247E+00, &
         0.7426682E+00, 0.1172451E+01, 0.1576190E+01, 0.2114785E+01, &
         0.2583662E+01, 0.2928657E+01, 0.3031442E+01, 0.2576708E+01/),&
         (/30,1,6/))

    qscat(1:30,1:1,43:48)=reshape( (/&
         0.8111209E-08, 0.1559276E-07, 0.3674471E-07, 0.8658950E-07, &
         0.2040660E-06, 0.4809785E-06, 0.1133488E-05, 0.2671544E-05, &
         0.6297743E-05, 0.1484709E-04, 0.3501190E-04, 0.8259039E-04, &
         0.1949215E-03, 0.4603695E-03, 0.1088433E-02, 0.2576899E-02, &
         0.6110877E-02, 0.1451004E-01, 0.3441217E-01, 0.8077721E-01, &
         0.1826787E+00, 0.3740697E+00, 0.6499188E+00, 0.1041741E+01, &
         0.1464437E+01, 0.1978474E+01, 0.2464825E+01, 0.2907833E+01, &
         0.3037615E+01, 0.2755203E+01, 0.6604543E-08, 0.1270030E-07, &
         0.2991808E-07, 0.7050917E-07, 0.1661954E-06, 0.3916392E-06, &
         0.9230437E-06, 0.2175630E-05, 0.5128434E-05, 0.1209001E-04, &
         0.2850855E-04, 0.6724245E-04, 0.1586770E-03, 0.3746897E-03, &
         0.8856207E-03, 0.2095975E-02, 0.4968437E-02, 0.1179495E-01, &
         0.2799396E-01, 0.6597276E-01, 0.1511538E+00, 0.3193898E+00, &
         0.5772527E+00, 0.9326557E+00, 0.1375426E+01, 0.1838931E+01, &
         0.2335359E+01, 0.2806430E+01, 0.3088753E+01, 0.2879159E+01, &
         0.5314572E-08, 0.1022444E-07, 0.2410669E-07, 0.5677523E-07, &
         0.1338119E-06, 0.3153663E-06, 0.7432756E-06, 0.1751925E-05, &
         0.4129471E-05, 0.9734469E-05, 0.2295285E-04, 0.5413433E-04, &
         0.1277282E-03, 0.3015520E-03, 0.7125568E-03, 0.1685775E-02, &
         0.3994397E-02, 0.9479688E-02, 0.2250814E-01, 0.5320236E-01, &
         0.1231845E+00, 0.2675992E+00, 0.5053291E+00, 0.8274198E+00, &
         0.1273818E+01, 0.1690028E+01, 0.2211889E+01, 0.2669053E+01, &
         0.2996819E+01, 0.2998283E+01, 0.4225172E-08, 0.8133153E-08, &
         0.1916917E-07, 0.4515045E-07, 0.1064213E-06, 0.2507710E-06, &
         0.5910236E-06, 0.1392893E-05, 0.3283100E-05, 0.7739130E-05, &
         0.1824681E-04, 0.4303140E-04, 0.1015181E-03, 0.2396298E-03, &
         0.5660809E-03, 0.1338762E-02, 0.3170775E-02, 0.7522085E-02, &
         0.1786233E-01, 0.4231057E-01, 0.9877851E-01, 0.2196952E+00, &
         0.4338378E+00, 0.7289849E+00, 0.1154307E+01, 0.1559031E+01, &
         0.2097416E+01, 0.2569233E+01, 0.2926430E+01, 0.3029643E+01, &
         0.3311696E-08, 0.6376909E-08, 0.1502053E-07, 0.3541976E-07, &
         0.8341913E-07, 0.1966386E-06, 0.4633977E-06, 0.1092155E-05, &
         0.2574293E-05, 0.6068137E-05, 0.1430610E-04, 0.3373469E-04, &
         0.7957646E-04, 0.1878019E-03, 0.4435409E-03, 0.1048594E-02, &
         0.2482403E-02, 0.5886369E-02, 0.1397656E-01, 0.3315251E-01, &
         0.7788367E-01, 0.1765939E+00, 0.3638026E+00, 0.6363935E+00, &
         0.1021569E+01, 0.1448428E+01, 0.1954509E+01, 0.2442554E+01, &
         0.2897095E+01, 0.3052830E+01, 0.2558936E-08, 0.4921040E-08, &
         0.1160372E-07, 0.2734706E-07, 0.6442724E-07, 0.1518185E-06, &
         0.3578024E-06, 0.8431982E-06, 0.1987328E-05, 0.4684338E-05, &
         0.1104337E-04, 0.2603939E-04, 0.6141687E-04, 0.1449215E-03, &
         0.3421815E-03, 0.8086877E-03, 0.1913604E-02, 0.4535340E-02, &
         0.1076550E-01, 0.2555612E-01, 0.6031040E-01, 0.1388400E+00, &
         0.2969850E+00, 0.5466675E+00, 0.8872278E+00, 0.1334345E+01, &
         0.1775121E+01, 0.2280190E+01, 0.2743268E+01, 0.3065531E+01/),&
         (/30,1,6/))

    qscat(1:30,1:1,49:50)=reshape( (/&
         0.1941416E-08, 0.3729086E-08, 0.8796222E-08, 0.2073586E-07, &
         0.4884795E-07, 0.1151409E-06, 0.2713934E-06, 0.6395728E-06, &
         0.1507417E-05, 0.3552973E-05, 0.8375334E-05, 0.1974695E-04, &
         0.4657068E-04, 0.1098731E-03, 0.2593660E-03, 0.6127629E-03, &
         0.1449341E-02, 0.3433179E-02, 0.8145735E-02, 0.1934303E-01, &
         0.4578911E-01, 0.1066269E+00, 0.2353837E+00, 0.4579038E+00, &
         0.7614718E+00, 0.1196600E+01, 0.1600340E+01, 0.2137544E+01, &
         0.2602569E+01, 0.2934595E+01, 0.1327624E-08, 0.2548038E-08, &
         0.5997629E-08, 0.1413148E-07, 0.3332353E-07, 0.7850970E-07, &
         0.1850432E-06, 0.4361592E-06, 0.1027941E-05, 0.2423005E-05, &
         0.5711535E-05, 0.1346519E-04, 0.3175162E-04, 0.7489492E-04, &
         0.1767484E-03, 0.4174066E-03, 0.9867256E-03, 0.2335701E-02, &
         0.5537855E-02, 0.1314833E-01, 0.3119540E-01, 0.7337637E-01, &
         0.1670410E+00, 0.3474123E+00, 0.6147258E+00, 0.9890165E+00, &
         0.1422405E+01, 0.1913994E+01, 0.2404561E+01, 0.2871936E+01/),&
         (/30,1,2/))

    qbrqs(1:30,1:1,1:6)=reshape( (/&
         0.1901999E-05, 0.4716377E-12, 0.7370087E-12, 0.3116265E-12, &
         -0.2044577E-11, 0.1364757E-11, 0.2750097E-11, 0.1206867E-10, &
         0.7804748E-10, 0.2244802E-09, 0.7640885E-09, 0.2786042E-08, &
         0.1003014E-07, 0.3626222E-07, 0.1312045E-06, 0.4747086E-06, &
         0.1717792E-05, 0.6216306E-05, 0.2250494E-04, 0.8146375E-04, &
         0.2948020E-03, 0.1065178E-02, 0.3828139E-02, 0.1352974E-01, &
         0.4564070E-01, 0.1391708E+00, 0.3677734E+00, 0.7288581E+00, &
         0.1024653E+01, 0.1386942E+01, 0.3617603E-05, 0.2991549E-12, &
         -0.5700891E-12, 0.1586302E-11, 0.9922366E-12, 0.1795592E-11, &
         -0.3836203E-11, 0.1123241E-10, 0.2974527E-10, 0.8357062E-10, &
         0.3698330E-09, 0.1441129E-08, 0.5219182E-08, 0.1896117E-07, &
         0.6871815E-07, 0.2483693E-06, 0.8989706E-06, 0.3253124E-05, &
         0.1177489E-04, 0.4262225E-04, 0.1542823E-03, 0.5580288E-03, &
         0.2012467E-02, 0.7188866E-02, 0.2496183E-01, 0.8068376E-01, &
         0.2299847E+00, 0.5466217E+00, 0.8696541E+00, 0.1232296E+01, &
         0.6880682E-05, 0.1916953E-12, 0.4707819E-12, 0.3414841E-12, &
         0.4034255E-12, 0.1418572E-12, 0.3161054E-12, 0.1475615E-11, &
         0.1669863E-10, 0.3005890E-10, 0.1336998E-09, 0.5276825E-09, &
         0.1935742E-08, 0.6909711E-08, 0.2535454E-07, 0.9163686E-07, &
         0.3316923E-06, 0.1200273E-05, 0.4344654E-05, 0.1572453E-04, &
         0.5692060E-04, 0.2060147E-03, 0.7448617E-03, 0.2682767E-02, &
         0.9544524E-02, 0.3276870E-01, 0.1033715E+00, 0.2852662E+00, &
         0.6315490E+00, 0.9307482E+00, 0.1308705E-04, 0.3488645E-12, &
         -0.2643488E-12, 0.3843259E-12, 0.5917927E-12, 0.1269802E-11, &
         -0.1910060E-11, 0.4475100E-11, 0.2189990E-11, 0.1775183E-10, &
         0.7611258E-10, 0.2165611E-09, 0.8112211E-09, 0.2849104E-08, &
         0.1058281E-07, 0.3819969E-07, 0.1382846E-06, 0.5006514E-06, &
         0.1810602E-05, 0.6554102E-05, 0.2372228E-04, 0.8587552E-04, &
         0.3107595E-03, 0.1122709E-02, 0.4033369E-02, 0.1423876E-01, &
         0.4789057E-01, 0.1452604E+00, 0.3812047E+00, 0.7420156E+00, &
         0.2489158E-04, 0.4701444E-13, 0.1535593E-12, 0.3381911E-12, &
         0.5717167E-12, 0.4095736E-12, 0.2617854E-12, 0.3672750E-11, &
         0.3647954E-11, 0.8980333E-11, 0.2740179E-10, 0.1001760E-09, &
         0.3395771E-09, 0.1195295E-08, 0.4350772E-08, 0.1593059E-07, &
         0.5727330E-07, 0.2076147E-06, 0.7506055E-06, 0.2718510E-05, &
         0.9838780E-05, 0.3561608E-04, 0.1289225E-03, 0.4663866E-03, &
         0.1683064E-02, 0.6024596E-02, 0.2104279E-01, 0.6894487E-01, &
         0.2004059E+00, 0.4936804E+00, 0.4734378E-04, 0.1351300E-13, &
         -0.1268058E-12, 0.2985778E-13, 0.3237699E-12, 0.6953172E-13, &
         -0.5454742E-12, 0.2173147E-11, 0.1272297E-11, 0.2537878E-11, &
         0.6615754E-11, 0.2416349E-10, 0.1269510E-09, 0.4061239E-09, &
         0.1475770E-08, 0.5395762E-08, 0.1915448E-07, 0.6942866E-07, &
         0.2513579E-06, 0.9103257E-06, 0.3293702E-05, 0.1192200E-04, &
         0.4315420E-04, 0.1562069E-03, 0.5649814E-03, 0.2037445E-02, &
         0.7276999E-02, 0.2525676E-01, 0.8155691E-01, 0.2321569E+00/),&
         (/30,1,6/))

    qbrqs(1:30,1:1,7:12)=reshape( (/&
         0.9004787E-04, 0.2465319E-14, 0.6956497E-13, 0.7135020E-13, &
         -0.7006778E-13, 0.4401354E-12, 0.4894050E-12, 0.6742820E-12, &
         -0.1286480E-11, 0.2726573E-12, 0.1833192E-11, 0.4639781E-11, &
         0.2582518E-10, 0.1105154E-09, 0.3573153E-09, 0.1309629E-08, &
         0.4652057E-08, 0.1709559E-07, 0.6191839E-07, 0.2242348E-06, &
         0.8105023E-06, 0.2933317E-05, 0.1061740E-04, 0.3843134E-04, &
         0.1391129E-03, 0.5032207E-03, 0.1815523E-02, 0.6493290E-02, &
         0.2262557E-01, 0.7371690E-01, 0.1712710E-03, 0.1129701E-13, &
         -0.3926947E-13, 0.2638894E-13, 0.5981981E-13, 0.1928346E-12, &
         0.1923591E-12, 0.7062969E-12, 0.8035807E-12, 0.1282500E-11, &
         0.7179135E-12, 0.1216088E-11, 0.4708282E-11, 0.2094201E-10, &
         0.5307086E-10, 0.2274950E-09, 0.8008884E-09, 0.2959651E-08, &
         0.1066964E-07, 0.3889165E-07, 0.1402413E-06, 0.5080287E-06, &
         0.1838745E-05, 0.6654002E-05, 0.2408813E-04, 0.8719498E-04, &
         0.3155253E-03, 0.1139886E-02, 0.4094611E-02, 0.1445007E-01, &
         0.3257575E-03, 0.8448067E-14, 0.9466076E-14, 0.6598221E-14, &
         0.5819855E-13, 0.1483693E-12, 0.2723606E-13, 0.1134216E-12, &
         -0.7536919E-12, 0.2585718E-12, 0.3713785E-13, 0.1066277E-11, &
         0.1599707E-12, 0.3999822E-12, 0.3919878E-11, 0.1730125E-10, &
         0.4113725E-10, 0.1999410E-09, 0.6582198E-09, 0.2469145E-08, &
         0.8993994E-08, 0.3211079E-07, 0.1164509E-06, 0.4218663E-06, &
         0.1527693E-05, 0.5526857E-05, 0.2000606E-04, 0.7242112E-04, &
         0.2620910E-03, 0.9472192E-03, 0.3353513E-01, 0.3950379E-04, &
         0.1429898E-03, 0.5172197E-03, 0.1865839E-02, 0.6671186E-02, &
         0.2322445E-01, 0.7551139E-01, 0.2170407E+00, 0.5241100E+00, &
         0.8542780E+00, 0.1208818E+01, 0.1538107E+01, 0.1891108E+01, &
         0.2099686E+01, 0.1850941E+01, 0.1153786E+01, 0.7584199E+00, &
         0.1408991E+01, 0.1537233E+01, 0.1201549E+01, 0.1408215E+01, &
         0.1379510E+01, 0.1359040E+01, 0.1355852E+01, 0.1382048E+01, &
         0.1385486E+01, 0.1393122E+01, 0.1399668E+01, 0.1368301E+01, &
         0.8438055E-05, 0.2306761E-04, 0.8350273E-04, 0.3021766E-03, &
         0.1091762E-02, 0.3922987E-02, 0.1385760E-01, 0.4668255E-01, &
         0.1419962E+00, 0.3740283E+00, 0.7350657E+00, 0.1032138E+01, &
         0.1392090E+01, 0.1722258E+01, 0.1978830E+01, 0.2038919E+01, &
         0.1466292E+01, 0.8217813E+00, 0.1068222E+01, 0.1678929E+01, &
         0.1179439E+01, 0.1465021E+01, 0.1322734E+01, 0.1388140E+01, &
         0.1398959E+01, 0.1402826E+01, 0.1374126E+01, 0.1401267E+01, &
         0.1384118E+01, 0.1380084E+01, 0.6045238E-05, 0.1652407E-04, &
         0.5981760E-04, 0.2164966E-03, 0.7826978E-03, 0.2818324E-02, &
         0.1001874E-01, 0.3432159E-01, 0.1077876E+00, 0.2957576E+00, &
         0.6457256E+00, 0.9422171E+00, 0.1321832E+01, 0.1654439E+01, &
         0.1947359E+01, 0.2037713E+01, 0.1652677E+01, 0.9366399E+00, &
         0.8927257E+00, 0.1617906E+01, 0.1334269E+01, 0.1372657E+01, &
         0.1323403E+01, 0.1372922E+01, 0.1323401E+01, 0.1366046E+01, &
         0.1377430E+01, 0.1352949E+01, 0.1383479E+01, 0.1394892E+01/),&
         (/30,1,6/))

    qbrqs(1:30,1:1,13:18)=reshape( (/&
         0.3857907E-05, 0.1054521E-04, 0.3817201E-04, 0.1381718E-03, &
         0.4998088E-03, 0.1803249E-02, 0.6449907E-02, 0.2247932E-01, &
         0.7327781E-01, 0.2114086E+00, 0.5139958E+00, 0.8473598E+00, &
         0.1197895E+01, 0.1527089E+01, 0.1880320E+01, 0.2101947E+01, &
         0.1867903E+01, 0.1170279E+01, 0.7712408E+00, 0.1393570E+01, &
         0.1562685E+01, 0.1217631E+01, 0.1419034E+01, 0.1369893E+01, &
         0.1385373E+01, 0.1366546E+01, 0.1356363E+01, 0.1383033E+01, &
         0.1387707E+01, 0.1393639E+01, 0.2875744E-05, 0.7860280E-05, &
         0.2845365E-04, 0.1030003E-03, 0.3726821E-03, 0.1345828E-02, &
         0.4827840E-02, 0.1696944E-01, 0.5645249E-01, 0.1680575E+00, &
         0.4297750E+00, 0.7848734E+00, 0.1099823E+01, 0.1440201E+01, &
         0.1776021E+01, 0.2050056E+01, 0.1986646E+01, 0.1322586E+01, &
         0.7947940E+00, 0.1190645E+01, 0.1659948E+01, 0.1198948E+01, &
         0.1479658E+01, 0.1347349E+01, 0.1393301E+01, 0.1416820E+01, &
         0.1365937E+01, 0.1370853E+01, 0.1368950E+01, 0.1383935E+01, &
         0.1718626E-05, 0.4697151E-05, 0.1700411E-04, 0.6155527E-04, &
         0.2227757E-03, 0.8053585E-03, 0.2899469E-02, 0.1030230E-01, &
         0.3524727E-01, 0.1104061E+00, 0.3019377E+00, 0.6538075E+00, &
         0.9490179E+00, 0.1328286E+01, 0.1660758E+01, 0.1947923E+01, &
         0.2039193E+01, 0.1638955E+01, 0.9184783E+00, 0.8948992E+00, &
         0.1630617E+01, 0.1313450E+01, 0.1395140E+01, 0.1328123E+01, &
         0.1377631E+01, 0.1357848E+01, 0.1388879E+01, 0.1391962E+01, &
         0.1361340E+01, 0.1369491E+01, 0.7404256E-06, 0.2023248E-05, &
         0.7321936E-05, 0.2650530E-04, 0.9594424E-04, 0.3471650E-03, &
         0.1253919E-02, 0.4500857E-02, 0.1584823E-01, 0.5295629E-01, &
         0.1588167E+00, 0.4104221E+00, 0.7685978E+00, 0.1076261E+01, &
         0.1422771E+01, 0.1755290E+01, 0.2024872E+01, 0.2010120E+01, &
         0.1377721E+01, 0.8087896E+00, 0.1148468E+01, 0.1692759E+01, &
         0.1168727E+01, 0.1450756E+01, 0.1345646E+01, 0.1390995E+01, &
         0.1421059E+01, 0.1388422E+01, 0.1380926E+01, 0.1383566E+01, &
         0.5598029E-06, 0.1530578E-05, 0.5538715E-05, 0.2004728E-04, &
         0.7256944E-04, 0.2626222E-03, 0.9491406E-03, 0.3413769E-02, &
         0.1209360E-01, 0.4104747E-01, 0.1265943E+00, 0.3394618E+00, &
         0.6988942E+00, 0.9913729E+00, 0.1363254E+01, 0.1694228E+01, &
         0.1954880E+01, 0.2044281E+01, 0.1563253E+01, 0.8779636E+00, &
         0.9851181E+00, 0.1671286E+01, 0.1245029E+01, 0.1393614E+01, &
         0.1312447E+01, 0.1383697E+01, 0.1368982E+01, 0.1410644E+01, &
         0.1386400E+01, 0.1377081E+01, 0.3768961E-06, 0.1029766E-05, &
         0.3728270E-05, 0.1349196E-04, 0.4884070E-04, 0.1767822E-03, &
         0.6393095E-03, 0.2304285E-02, 0.8216626E-02, 0.2838739E-01, &
         0.9074473E-01, 0.2547902E+00, 0.5869321E+00, 0.8975460E+00, &
         0.1271339E+01, 0.1603128E+01, 0.1935510E+01, 0.2057565E+01, &
         0.1761917E+01, 0.1026894E+01, 0.8074358E+00, 0.1538698E+01, &
         0.1437847E+01, 0.1310425E+01, 0.1362680E+01, 0.1381976E+01, &
         0.1356858E+01, 0.1333898E+01, 0.1381830E+01, 0.1388422E+01/),&
         (/30,1,6/))

    qbrqs(1:30,1:1,19:24)=reshape( (/&
         0.3236136E-06, 0.8855458E-06, 0.3204303E-05, 0.1159740E-04, &
         0.4198187E-04, 0.1519609E-03, 0.5496361E-03, 0.1982331E-02, &
         0.7082549E-02, 0.2460569E-01, 0.7962730E-01, 0.2273520E+00, &
         0.5421261E+00, 0.8665863E+00, 0.1227701E+01, 0.1557529E+01, &
         0.1907946E+01, 0.2091759E+01, 0.1820815E+01, 0.1128173E+01, &
         0.7792270E+00, 0.1441805E+01, 0.1500771E+01, 0.1177203E+01, &
         0.1397181E+01, 0.1378007E+01, 0.1360600E+01, 0.1342718E+01, &
         0.1381574E+01, 0.1394904E+01, 0.2598837E-06, 0.7105884E-06, &
         0.2571407E-05, 0.9306488E-05, 0.3368797E-04, 0.1219437E-03, &
         0.4411651E-03, 0.1592339E-02, 0.5703115E-02, 0.1995317E-01, &
         0.6563381E-01, 0.1919264E+00, 0.4775247E+00, 0.8215518E+00, &
         0.1156758E+01, 0.1487721E+01, 0.1835645E+01, 0.2094382E+01, &
         0.1909802E+01, 0.1222000E+01, 0.7729526E+00, 0.1325953E+01, &
         0.1618831E+01, 0.1129606E+01, 0.1441841E+01, 0.1368787E+01, &
         0.1390773E+01, 0.1396055E+01, 0.1353166E+01, 0.1377240E+01, &
         0.1468962E-06, 0.4008892E-06, 0.1451859E-05, 0.5252375E-05, &
         0.1901069E-04, 0.6881658E-04, 0.2490552E-03, 0.9001884E-03, &
         0.3238772E-02, 0.1148521E-01, 0.3908658E-01, 0.1211622E+00, &
         0.3270030E+00, 0.6846489E+00, 0.9770864E+00, 0.1352290E+01, &
         0.1683818E+01, 0.1951200E+01, 0.2043353E+01, 0.1586950E+01, &
         0.8647165E+00, 0.9542888E+00, 0.1668071E+01, 0.1251961E+01, &
         0.1387784E+01, 0.1315097E+01, 0.1377233E+01, 0.1364871E+01, &
         0.1410496E+01, 0.1380941E+01, 0.8878270E-07, 0.2427971E-06, &
         0.8784166E-06, 0.3178303E-05, 0.1150401E-04, 0.4164483E-04, &
         0.1507380E-03, 0.5452277E-03, 0.1966480E-02, 0.7026603E-02, &
         0.2441811E-01, 0.7907017E-01, 0.2259613E+00, 0.5397338E+00, &
         0.8649546E+00, 0.1225237E+01, 0.1554972E+01, 0.1905898E+01, &
         0.2093079E+01, 0.1824446E+01, 0.1131443E+01, 0.7718529E+00, &
         0.1433445E+01, 0.1505972E+01, 0.1179354E+01, 0.1397423E+01, &
         0.1375723E+01, 0.1354594E+01, 0.1339999E+01, 0.1368543E+01, &
         0.4909048E-07, 0.1342343E-06, 0.4857623E-06, 0.1758437E-05, &
         0.6363671E-05, 0.2303496E-04, 0.8338189E-04, 0.3017350E-03, &
         0.1090177E-02, 0.3917343E-02, 0.1383809E-01, 0.4662061E-01, &
         0.1418286E+00, 0.3736582E+00, 0.7347008E+00, 0.1031695E+01, &
         0.1391785E+01, 0.1721950E+01, 0.1978449E+01, 0.2039113E+01, &
         0.1467247E+01, 0.8228753E+00, 0.1067643E+01, 0.1679000E+01, &
         0.1180979E+01, 0.1465477E+01, 0.1323149E+01, 0.1388206E+01, &
         0.1398749E+01, 0.1403252E+01, 0.4491110E-07, 0.1226495E-06, &
         0.4438543E-06, 0.1606513E-05, 0.5813939E-05, 0.2104629E-04, &
         0.7618610E-04, 0.2757049E-03, 0.9963178E-03, 0.3582313E-02, &
         0.1267850E-01, 0.4292409E-01, 0.1317570E+00, 0.3511752E+00, &
         0.7116793E+00, 0.1005011E+01, 0.1373192E+01, 0.1703715E+01, &
         0.1960462E+01, 0.2044763E+01, 0.1533077E+01, 0.8851303E+00, &
         0.1014989E+01, 0.1669717E+01, 0.1245077E+01, 0.1428272E+01, &
         0.1312307E+01, 0.1375859E+01, 0.1382066E+01, 0.1413976E+01/),&
         (/30,1,6/))

    qbrqs(1:30,1:1,25:30)=reshape( (/&
         0.3771462E-07, 0.1028282E-06, 0.3719264E-06, 0.1346401E-05, &
         0.4873332E-05, 0.1763927E-04, 0.6385347E-04, 0.2310964E-03, &
         0.8353932E-03, 0.3006974E-02, 0.1067756E-01, 0.3646917E-01, &
         0.1138473E+00, 0.3100131E+00, 0.6640764E+00, 0.9579683E+00, &
         0.1336368E+01, 0.1668591E+01, 0.1948658E+01, 0.2041109E+01, &
         0.1620370E+01, 0.8906084E+00, 0.9080491E+00, 0.1650102E+01, &
         0.1287251E+01, 0.1404940E+01, 0.1312431E+01, 0.1369794E+01, &
         0.1347291E+01, 0.1383594E+01, 0.2066275E-07, 0.5672286E-07, &
         0.2053315E-06, 0.7420299E-06, 0.2685598E-05, 0.9722087E-05, &
         0.3519086E-04, 0.1273800E-03, 0.4608098E-03, 0.1663016E-02, &
         0.5953574E-02, 0.2080238E-01, 0.6821626E-01, 0.1985453E+00, &
         0.4901721E+00, 0.8306453E+00, 0.1171311E+01, 0.1501209E+01, &
         0.1851867E+01, 0.2099972E+01, 0.1896659E+01, 0.1208096E+01, &
         0.8045296E+00, 0.1344672E+01, 0.1607270E+01, 0.1197700E+01, &
         0.1446757E+01, 0.1368162E+01, 0.1369020E+01, 0.1384346E+01, &
         0.2943546E-07, 0.8037221E-07, 0.2906546E-06, 0.1051397E-05, &
         0.3806524E-05, 0.1377661E-04, 0.4987255E-04, 0.1805111E-03, &
         0.6527653E-03, 0.2352573E-02, 0.8386328E-02, 0.2895008E-01, &
         0.9238102E-01, 0.2587797E+00, 0.5930740E+00, 0.9019188E+00, &
         0.1276954E+01, 0.1608966E+01, 0.1937796E+01, 0.2052794E+01, &
         0.1750075E+01, 0.1011107E+01, 0.8058437E+00, 0.1549847E+01, &
         0.1425072E+01, 0.1307907E+01, 0.1361917E+01, 0.1381650E+01, &
         0.1353116E+01, 0.1336342E+01, 0.1401500E-07, 0.3835394E-07, &
         0.1390823E-06, 0.5033391E-06, 0.1820643E-05, 0.6589878E-05, &
         0.2385175E-04, 0.8634376E-04, 0.3124391E-03, 0.1128766E-02, &
         0.4054963E-02, 0.1431323E-01, 0.4812629E-01, 0.1458958E+00, &
         0.3825957E+00, 0.7433397E+00, 0.1042454E+01, 0.1399161E+01, &
         0.1729491E+01, 0.1988325E+01, 0.2033609E+01, 0.1445695E+01, &
         0.8003992E+00, 0.1082630E+01, 0.1693802E+01, 0.1152685E+01, &
         0.1457584E+01, 0.1329207E+01, 0.1388109E+01, 0.1413331E+01, &
         0.9938845E-08, 0.2700967E-07, 0.9789189E-07, 0.3541029E-06, &
         0.1281403E-05, 0.4637087E-05, 0.1678571E-04, 0.6076441E-04, &
         0.2199208E-03, 0.7950518E-03, 0.2862566E-02, 0.1017336E-01, &
         0.3482668E-01, 0.1092176E+00, 0.2991365E+00, 0.6501684E+00, &
         0.9459308E+00, 0.1325391E+01, 0.1657931E+01, 0.1947675E+01, &
         0.2038492E+01, 0.1645274E+01, 0.9275196E+00, 0.8941895E+00, &
         0.1624768E+01, 0.1322001E+01, 0.1387358E+01, 0.1331443E+01, &
         0.1374707E+01, 0.1341997E+01, 0.6412746E-08, 0.1769780E-07, &
         0.6366894E-07, 0.2303611E-06, 0.8336778E-06, 0.3016920E-05, &
         0.1092051E-04, 0.3952926E-04, 0.1430871E-03, 0.5175775E-03, &
         0.1867124E-02, 0.6675717E-02, 0.2323973E-01, 0.7555697E-01, &
         0.2171555E+00, 0.5243142E+00, 0.8544180E+00, 0.1209036E+01, &
         0.1538329E+01, 0.1891320E+01, 0.2099623E+01, 0.1850582E+01, &
         0.1153475E+01, 0.7582853E+00, 0.1409339E+01, 0.1536727E+01, &
         0.1201235E+01, 0.1408159E+01, 0.1379772E+01, 0.1358475E+01/),&
         (/30,1,6/))

    qbrqs(1:30,1:1,31:36)=reshape( (/&
         0.3923926E-08, 0.1067798E-07, 0.3897862E-07, 0.1409177E-06, &
         0.5099691E-06, 0.1845879E-05, 0.6680979E-05, 0.2418180E-04, &
         0.8753642E-04, 0.3167605E-03, 0.1144346E-02, 0.4110507E-02, &
         0.1450484E-01, 0.4873199E-01, 0.1475266E+00, 0.3861555E+00, &
         0.7467078E+00, 0.1046755E+01, 0.1402116E+01, 0.1732574E+01, &
         0.1992608E+01, 0.2031040E+01, 0.1437798E+01, 0.7941419E+00, &
         0.1091892E+01, 0.1695350E+01, 0.1146589E+01, 0.1452048E+01, &
         0.1327341E+01, 0.1386919E+01, 0.2159660E-08, 0.5853871E-08, &
         0.2139184E-07, 0.7725549E-07, 0.2794547E-06, 0.1011294E-05, &
         0.3659937E-05, 0.1324533E-04, 0.4794678E-04, 0.1735502E-03, &
         0.6276303E-03, 0.2262394E-02, 0.8069280E-02, 0.2789819E-01, &
         0.8931855E-01, 0.2513025E+00, 0.5814862E+00, 0.8937050E+00, &
         0.1266289E+01, 0.1597861E+01, 0.1933207E+01, 0.2062025E+01, &
         0.1771295E+01, 0.1042512E+01, 0.8124586E+00, 0.1524010E+01, &
         0.1449887E+01, 0.1308147E+01, 0.1367361E+01, 0.1380090E+01, &
         0.7260568E-09, 0.2020068E-08, 0.7262182E-08, 0.2583743E-07, &
         0.9367815E-07, 0.3389539E-06, 0.1226493E-05, 0.4437584E-05, &
         0.1606419E-04, 0.5814762E-04, 0.2104592E-03, 0.7609050E-03, &
         0.2740245E-02, 0.9745683E-02, 0.3342814E-01, 0.1052505E+00, &
         0.2897407E+00, 0.6376665E+00, 0.9356289E+00, 0.1315299E+01, &
         0.1647984E+01, 0.1946704E+01, 0.2036762E+01, 0.1665433E+01, &
         0.9483317E+00, 0.8870602E+00, 0.1611704E+01, 0.1349355E+01, &
         0.1338266E+01, 0.1332024E+01, 0.2729792E-09, 0.7788261E-09, &
         0.2781656E-08, 0.1002817E-07, 0.3617201E-07, 0.1306958E-06, &
         0.4733878E-06, 0.1714013E-05, 0.6202667E-05, 0.2245028E-04, &
         0.8127029E-04, 0.2940956E-03, 0.1062638E-02, 0.3819079E-02, &
         0.1349840E-01, 0.4554095E-01, 0.1388998E+00, 0.3671713E+00, &
         0.7282541E+00, 0.1023936E+01, 0.1386446E+01, 0.1716618E+01, &
         0.1972223E+01, 0.2042006E+01, 0.1485218E+01, 0.8435971E+00, &
         0.1052239E+01, 0.1682687E+01, 0.1208819E+01, 0.1468988E+01, &
         0.1138483E-09, 0.3352924E-09, 0.1323467E-08, 0.4765051E-08, &
         0.1696672E-07, 0.6176972E-07, 0.2235700E-06, 0.8091876E-06, &
         0.2928622E-05, 0.1059961E-04, 0.3837255E-04, 0.1388960E-03, &
         0.5024281E-03, 0.1812669E-02, 0.6483227E-02, 0.2259160E-01, &
         0.7361497E-01, 0.2122605E+00, 0.5155379E+00, 0.8483830E+00, &
         0.1199574E+01, 0.1528771E+01, 0.1882021E+01, 0.2101714E+01, &
         0.1865457E+01, 0.1167623E+01, 0.7686450E+00, 0.1395737E+01, &
         0.1559300E+01, 0.1215185E+01, 0.6871950E-10, 0.1435970E-09, &
         0.5198381E-09, 0.1935983E-08, 0.7118665E-08, 0.2559486E-07, &
         0.9233558E-07, 0.3342500E-06, 0.1209471E-05, 0.4377391E-05, &
         0.1584374E-04, 0.5735440E-04, 0.2075856E-03, 0.7505263E-03, &
         0.2703073E-02, 0.9615588E-02, 0.3300178E-01, 0.1040363E+00, &
         0.2868509E+00, 0.6337283E+00, 0.9324754E+00, 0.1312067E+01, &
         0.1644767E+01, 0.1946319E+01, 0.2036625E+01, 0.1671548E+01, &
         0.9529492E+00, 0.8828400E+00, 0.1606683E+01, 0.1356629E+01/),&
         (/30,1,6/))

    qbrqs(1:30,1:1,37:42)=reshape( (/&
         0.4279296E-10, 0.1019830E-09, 0.3495865E-09, 0.1240186E-08, &
         0.4548052E-08, 0.1637553E-07, 0.5894112E-07, 0.2132708E-06, &
         0.7721861E-06, 0.2795869E-05, 0.1012025E-04, 0.3663430E-04, &
         0.1326112E-03, 0.4797276E-03, 0.1731045E-02, 0.6194473E-02, &
         0.2161726E-01, 0.7068184E-01, 0.2048290E+00, 0.5019357E+00, &
         0.8389409E+00, 0.1184580E+01, 0.1513944E+01, 0.1866371E+01, &
         0.2102310E+01, 0.1884416E+01, 0.1191447E+01, 0.7944575E+00, &
         0.1358286E+01, 0.1585731E+01, 0.1320506E-10, 0.3754272E-10, &
         0.1172757E-09, 0.4610952E-09, 0.1747573E-08, 0.6396914E-08, &
         0.2307197E-07, 0.8355009E-07, 0.3025998E-06, 0.1094051E-05, &
         0.3960026E-05, 0.1433091E-04, 0.5187821E-04, 0.1877786E-03, &
         0.6790176E-03, 0.2446735E-02, 0.8717073E-02, 0.3004434E-01, &
         0.9555029E-01, 0.2664733E+00, 0.6046580E+00, 0.9103130E+00, &
         0.1287310E+01, 0.1619669E+01, 0.1941297E+01, 0.2045045E+01, &
         0.1725672E+01, 0.9880461E+00, 0.8178003E+00, 0.1567536E+01, &
         0.5899305E-11, 0.2042627E-10, 0.3764592E-10, 0.1280842E-09, &
         0.4243975E-09, 0.1496634E-08, 0.5502207E-08, 0.1978670E-07, &
         0.7175582E-07, 0.2595964E-06, 0.9389952E-06, 0.3399184E-05, &
         0.1230400E-04, 0.4453738E-04, 0.1612071E-03, 0.5830508E-03, &
         0.2102342E-02, 0.7505785E-02, 0.2602144E-01, 0.8381479E-01, &
         0.2377558E+00, 0.5596477E+00, 0.8785649E+00, 0.1245335E+01, &
         0.1575926E+01, 0.1921143E+01, 0.2079966E+01, 0.1798919E+01, &
         0.1101398E+01, 0.8241425E+00, 0.1203991E-11, 0.5900688E-11, &
         0.1824277E-10, 0.5389456E-10, 0.1671786E-09, 0.5994772E-09, &
         0.2299403E-08, 0.8175137E-08, 0.2941320E-07, 0.1066256E-06, &
         0.3856113E-06, 0.1395808E-05, 0.5051475E-05, 0.1828509E-04, &
         0.6619136E-04, 0.2395459E-03, 0.8658932E-03, 0.3116101E-02, &
         0.1105803E-01, 0.3770426E-01, 0.1173082E+00, 0.3180817E+00, &
         0.6740167E+00, 0.9669973E+00, 0.1344094E+01, 0.1676006E+01, &
         0.1949607E+01, 0.2042465E+01, 0.1603115E+01, 0.8707204E+00, &
         0.3155373E-11, 0.1247106E-11, 0.9425561E-11, 0.2065244E-10, &
         0.7427725E-10, 0.2928719E-09, 0.1081926E-08, 0.3947230E-08, &
         0.1429128E-07, 0.5198186E-07, 0.1883016E-06, 0.6817618E-06, &
         0.2467739E-05, 0.8929537E-05, 0.3232451E-04, 0.1170097E-03, &
         0.4233310E-03, 0.1528165E-02, 0.5475521E-02, 0.1917978E-01, &
         0.6327049E-01, 0.1858341E+00, 0.4656517E+00, 0.8128141E+00, &
         0.1142859E+01, 0.1475338E+01, 0.1820214E+01, 0.2086410E+01, &
         0.1925889E+01, 0.1240782E+01, 0.1531637E-11, 0.2570500E-11, &
         0.1066159E-11, 0.1496788E-10, 0.5043575E-10, 0.1727316E-09, &
         0.6206500E-09, 0.2259410E-08, 0.8238917E-08, 0.2965713E-07, &
         0.1075832E-06, 0.3893423E-06, 0.1409170E-05, 0.5098709E-05, &
         0.1845820E-04, 0.6681245E-04, 0.2417935E-03, 0.8740014E-03, &
         0.3145104E-02, 0.1115909E-01, 0.3803171E-01, 0.1182230E+00, &
         0.3202052E+00, 0.6765814E+00, 0.9693897E+00, 0.1346075E+01, &
         0.1677899E+01, 0.1949971E+01, 0.2042721E+01, 0.1599026E+01/),&
         (/30,1,6/))

    qbrqs(1:30,1:1,43:48)=reshape( (/&
         0.9917659E-12, 0.1381379E-11, 0.3523106E-11, 0.9577097E-11, &
         0.3459986E-10, 0.1411091E-09, 0.4419279E-09, 0.1553517E-08, &
         0.5765216E-08, 0.2065135E-07, 0.7509527E-07, 0.2717656E-06, &
         0.9833307E-06, 0.3558670E-05, 0.1288058E-04, 0.4662503E-04, &
         0.1687613E-03, 0.6103279E-03, 0.2200291E-02, 0.7850773E-02, &
         0.2717153E-01, 0.8719374E-01, 0.2460894E+00, 0.5732117E+00, &
         0.8879268E+00, 0.1258482E+01, 0.1589697E+01, 0.1929189E+01, &
         0.2068972E+01, 0.1783306E+01, 0.3368454E-12, 0.1753294E-11, &
         0.1704585E-12, 0.3925133E-11, 0.2875746E-10, 0.8277545E-10, &
         0.3136528E-09, 0.1162799E-08, 0.4243529E-08, 0.1521952E-07, &
         0.5533235E-07, 0.1997696E-06, 0.7227579E-06, 0.2614614E-05, &
         0.9465122E-05, 0.3426161E-04, 0.1240174E-03, 0.4486567E-03, &
         0.1619298E-02, 0.5798678E-02, 0.2027742E-01, 0.6662144E-01, &
         0.1944622E+00, 0.4824011E+00, 0.8250818E+00, 0.1162402E+01, &
         0.1492890E+01, 0.1841950E+01, 0.2096908E+01, 0.1904410E+01, &
         -0.1283701E-11, 0.2700828E-12, 0.5694661E-11, 0.1144088E-11, &
         0.1337860E-10, 0.5726562E-10, 0.2234361E-09, 0.8505050E-09, &
         0.3057868E-08, 0.1090547E-07, 0.3982935E-07, 0.1442852E-06, &
         0.5222154E-06, 0.1889433E-05, 0.6839028E-05, 0.2475330E-04, &
         0.8960722E-04, 0.3242505E-03, 0.1171335E-02, 0.4206709E-02, &
         0.1483641E-01, 0.4977831E-01, 0.1503363E+00, 0.3922572E+00, &
         0.7523814E+00, 0.1054152E+01, 0.1407214E+01, 0.1737990E+01, &
         0.2000353E+01, 0.2026293E+01, 0.1193975E-11, 0.1670941E-11, &
         0.5259805E-11, 0.4074115E-11, 0.2072773E-10, 0.5170789E-10, &
         0.1838727E-09, 0.6078937E-09, 0.2173212E-08, 0.7777988E-08, &
         0.2828418E-07, 0.1023848E-06, 0.3703800E-06, 0.1340259E-05, &
         0.4847633E-05, 0.1754610E-04, 0.6351911E-04, 0.2298873E-03, &
         0.8310285E-03, 0.2991342E-02, 0.1062301E-01, 0.3629179E-01, &
         0.1133487E+00, 0.3088466E+00, 0.6626130E+00, 0.9566727E+00, &
         0.1335223E+01, 0.1667487E+01, 0.1948545E+01, 0.2040859E+01, &
         -0.1518614E-11, 0.1025309E-11, 0.2377287E-12, 0.8970014E-11, &
         0.3335786E-11, 0.3682587E-10, 0.1141844E-09, 0.4113278E-09, &
         0.1523402E-08, 0.5454212E-08, 0.1970126E-07, 0.7100462E-07, &
         0.2571013E-06, 0.9296294E-06, 0.3366042E-05, 0.1218411E-04, &
         0.4409765E-04, 0.1596165E-03, 0.5773018E-03, 0.2081709E-02, &
         0.7433035E-02, 0.2577847E-01, 0.8309824E-01, 0.2359817E+00, &
         0.5567063E+00, 0.8765479E+00, 0.1242425E+01, 0.1572880E+01, &
         0.1919152E+01, 0.2082179E+01, 0.3449489E-12, 0.1323594E-12, &
         0.2524434E-11, 0.7129730E-11, 0.1099827E-10, 0.2720793E-10, &
         0.9374397E-10, 0.2881573E-09, 0.1014721E-08, 0.3635337E-08, &
         0.1330892E-07, 0.4817363E-07, 0.1743618E-06, 0.6308808E-06, &
         0.2282745E-05, 0.8262761E-05, 0.2990859E-04, 0.1082676E-03, &
         0.3917242E-03, 0.1414399E-02, 0.5071576E-02, 0.1780289E-01, &
         0.5903469E-01, 0.1748266E+00, 0.4436502E+00, 0.7959911E+00, &
         0.1116611E+01, 0.1453352E+01, 0.1792378E+01, 0.2065975E+01/),&
         (/30,1,6/))

    qbrqs(1:30,1:1,49:50)=reshape( (/&
         -0.9806707E-13, 0.1529088E-11, 0.2389664E-12, 0.2416812E-11, &
         -0.3707256E-12, 0.1301249E-10, 0.6405868E-10, 0.1982637E-09, &
         0.7039382E-09, 0.2456551E-08, 0.8758652E-08, 0.3174901E-07, &
         0.1151133E-06, 0.4168501E-06, 0.1507989E-05, 0.5457498E-05, &
         0.1975455E-04, 0.7150998E-04, 0.2587977E-03, 0.9353363E-03, &
         0.3364424E-02, 0.1192217E-01, 0.4049582E-01, 0.1250701E+00, &
         0.3359799E+00, 0.6949816E+00, 0.9873564E+00, 0.1360237E+01, &
         0.1691364E+01, 0.1953644E+01, 0.1310020E-11, 0.7826037E-12, &
         -0.8811858E-12, 0.2869265E-11, 0.2279223E-11, 0.5720607E-12, &
         0.1676336E-10, 0.9887384E-10, 0.3575822E-09, 0.1381762E-08, &
         0.4984234E-08, 0.1798714E-07, 0.6499247E-07, 0.2345360E-06, &
         0.8492439E-06, 0.3073006E-05, 0.1112344E-04, 0.4026809E-04, &
         0.1457573E-03, 0.5272286E-03, 0.1901805E-02, 0.6798243E-02, &
         0.2365158E-01, 0.7678756E-01, 0.2202472E+00, 0.5297817E+00, &
         0.8581598E+00, 0.1214843E+01, 0.1544259E+01, 0.1896759E+01/),&
         (/30,1,2/))

    ! ------------------------------------------------------------------

   
    IF (l_read ) THEN
       all_ok = (mwave==nwave) .and. (mrad==nrad) .and. (mgroup==ngroup)
       DO ig = 1, ngroup
          !...dbg:
          all_ok = all_ok .and. ( r(1,ig) .eq. rmin(ig) )
       END DO
       !
       IF ( .not. all_ok )THEN
          WRITE(lunoprt,*) ' setuprad: mie.data grid(s) bad: '
          WRITE(lunoprt,*) ' in mie.data, mwave, mrad, mgroup = ', &
               mwave, mrad, mgroup
          WRITE(lunoprt,*) ' and rmin = ', rmin
          call abort_run ('Bad set up data...','calcproperties','rad_carma.f90')
       END IF
    ELSE
       ! calculate extinction and scattering coefficients
       !
       DO ig = 1,ngroup
          !
          !	select ice/liquid index of refractive index array
          !
          IF( is_grp_ice(ig) )THEN
             irefr = 2
          ELSE
             irefr = 1
          END IF
          !
          !	<thetd> is angle between incident and scattered radiation
          !	<j_thetd> is number of <thetd> values to consider
          !
          thetd = 0.0
          DO  l=1,nwave
             !
             !kml for biomass burning particles:
             r_real = 1.495
             tmag = 1.0e-6
  
             !	  real = 1.52
             !	  tmag = 0.015
             
             !
             !	 calculate the center of the wavelength interval of an ir interval
             !
             IF( l .le. nsol ) THEN
                awave = wave(l)
             ELSE
                awave = 0.5*(wave(l)+wave(l+1))
             END IF
             wvno     =	twopi/(awave*1.0e-4)
             !
             DO i=1,nrad
                IF(i .eq. 1) THEN
                   ddr     = 0.2*(rup(1,ig)-r(1,ig))
                   rr      = r(1,ig)
                   corerad = rcore(1,ig)
                   ddrc    = 0.2*(rcoreup(i,ig)- rcore(1,ig))
                ELSE
                   ddr     = 0.2*(rup(i,ig)-rup(i-1,ig))
                   rr      = rup(i-1,ig)
                   corerad = rcoreup(i-1,ig)
                   ddrc    = 0.2*(rcoreup(i,ig)-rcoreup(i-1,ig))
                END IF
                !
                rdqext(i,ig,l) = 0.0
                qscat(i,ig,l)  = 0.0
                qbrqs(i,ig,l)  = 0.0
                !
                DO j=1,6
                   !
                   !
                   ! limit x=twopi/wave to no larger 1000 to avoid anguish in the mie code.
                   !
                   IF( wvno*rr .gt. 1000. ) rr = 1000./wvno

                   CALL miess(rr,r_real,tmag,thetd,qextd,qscatd,ctbrqs,&
                        corerad,corereal,coreimag,wvno)
                   !     IF(l.eq.7 .or. l.eq.8)
                   ! &      print*,'qex qsc=',qextd,qscatd
                   !
                   rdqext(i,ig,l) = rdqext(i,ig,l)+qextd/6.
                   qscat(i,ig,l)  = qscat(i,ig,l)+qscatd/6.
                   qbrqs(i,ig,l)  = qbrqs(i,ig,l)+ctbrqs/6.
                   rr             = rr+ddr
                   corerad        = corerad + ddrc
                   !
                END DO
             END DO
          END DO
          !
          !	 stop
       END DO	  ! ig=1,ngroup
       !
    END IF
    !
    ! *********************************************************************
    !
    !				 check sum of weights
    !
    ! **********************************************************************
    !
    sum  = 0.0
    sum1 = 0.0
    sum2 = 0.0
    DO l = 1,nsolp
       sum  = sum+weight(l)
    END DO
    DO l = nsolp+1,ntotal
       sum1 = sum1+weight(l)
    END DO
    sum2	=   sum+sum1
  
    !
    IF ( abs(nwave-sum2) .gt. 1.e-3 ) WRITE(lunoprt,FMT=lab355) sum,sum1,sum2
    !
    DO l = 1,nsolp
       sol(l) = solfx(nprob(l)) * weight(l)
    END DO
  
    !
    ! *********************************************************************
    !
    !	compute planck function table. wave is in units of microns.
    !
    ! **********************************************************************
    !
    !	set <iblackbody_above> = .true. to include a source of radiation
    !	at the top of the radiative transfer model DOmain
    
    iblackbody_above = .false.
    t_above = tabove_aerad
    !
    !	set <ibeyond_spectrum> = 1 to include blackbody radiation at
    !	wavelengths longer than wave(nwave+1) in planck(nwave+1-nsol) and
    !	at wavelengths shorter than wave(nsol+1) in planck(1)
    !
    ibeyond_spectrum = 1
    !
    IF( ibeyond_spectrum .eq. 1 )THEN 
       DO j = 1,ncount
          planck(nwave+1-nsol,j) = (0.01*float(nlow+j))**4
       END DO
       DO i = nsol+2,nwave
          DO j = 1,ncount
             k = i-nsol
             v = 1.438e4 / wave(i)
             CALL plnk(v,(0.01*float(nlow+j)),planck(k,j))
          END DO
       END DO
    ELSE 
       DO i = nsol+1,nwave+1
          DO j = 1,ncount
             k = i-nsol
             v = 1.438e4 / wave(i)
             CALL plnk(v,(0.01*float(nlow+j)),planck(k,j))
          END DO
       END DO
    END IF
    !
    DO j = 1,ncount
  
       IF( ibeyond_spectrum .eq. 1 )THEN
  
          planck(1,j) = planck(2,j) * stefan * pii
          DO l = nsol+2,nwave
             k = l-nsol
             planck(k,j) = (planck(k+1,j)-planck(k,j)) * stefan * pii
          END DO
  
       ELSE
  
          DO l = nsol+1,nwave
             k = l-nsol
             planck(k,j) = (planck(k+1,j)-planck(k,j))*stefan * pii
          END DO
        
       END IF
  
    END DO
    !
  
  END SUBROUTINE calcproperties
   !=======================================================================================!
   !=======================================================================================!







   !=======================================================================================!
   !=======================================================================================!
   !    This subroutines loads profiles of temperature [k], water vapour [g/g], and        !
   ! aerosol particle concentrations [#/cm^2].                                             !
   !---------------------------------------------------------------------------------------!
   subroutine prerad(m1,nzpmax,dztr,fmapt)
      use mem_aerad  , only : nbin    ! !
      use mem_globaer, only : ngroup  & !
                            , ienconc & !
                            , nelem   ! ! 
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      integer                   , intent(in) :: m1
      integer                   , intent(in) :: nzpmax
      real   , dimension(nzpmax), intent(in) :: dztr
      real                      , intent(in) :: fmapt
      !----- Local variables. -------------------------------------------------------------!
      integer                                :: ibin
      integer                                :: iep
      integer                                :: igas
      integer                                :: igroup
      integer                                :: k,kk
      integer                                :: nzz
      real                                   :: xymet
      !------------------------------------------------------------------------------------!
      !
      igas = 1
      nzz = m1 - 1
  
      do k = 1,nzz
         !----- Reverse the vertical index when in cartesian coordinates. -----------------!
         kk = nzz + 1 - k
         !----- qv-aerad has g[h20]/g[ar] units in the radiation code... ------------------!
         qv_aerad(kk) = gc(k,igas) / rhoa(k)
         do igroup = 1,ngroup
            iep = ienconc(igroup)
            do ibin = 1,nbin
               xymet = fmapt*fmapt
               pc_aerad(kk,ibin,igroup) = pc(k,ibin,iep) * (1./dztr(k)) / xymet
            end do
         end do
      end do
      
      return
   end subroutine prerad
   !=======================================================================================!
   !=======================================================================================!







   !=======================================================================================!
   !=======================================================================================!
   subroutine radtran(albedt,cosz,m1,aot11)
      use rconstants  , only : pi1             & ! intent(in)
                             , gcgs            & ! intent(in)
                             , day_sec         ! ! intent(in)
      use mem_aerad   , only : u0_aerad        & ! intent()
                             , qrad_aerad      & ! intent()
                             , alb_toai_aerad  & ! intent()
                             , alb_tomi_aerad  & ! intent()
                             , alb_toa_aerad   & ! intent()
                             , fsl_up_aerad    & ! intent()
                             , fsl_dn_aerad    & ! intent()
                             , fir_up_aerad    & ! intent()
                             , fir_dn_aerad    & ! intent()
                             , nir             ! ! intent()
      use mem_globrad , only : isl             & ! intent()
                             , nvert           & ! intent()
                             , nlayer          & ! intent()
                             , ngroup          & ! intent()
                             , nrad            & ! intent()
                             , u0              & ! intent()
                             , nsolp           & ! intent()
                             , albedo_sfc      & ! intent()
                             , emis            & ! intent()
                             , ntotal          & ! intent()
                             , emisir          & ! intent()
                             , ir              & ! intent()
                             , irs             & ! intent()
                             , fdegday         & ! intent()
                             , qrad            & ! intent()
                             , roundoff        & ! intent()
                             , xsecta          & ! intent()
                             , rdqext          & ! intent()
                             , nprob           & ! intent()
                             , qscat           & ! intent()
                             , nwave           & ! intent()
                             , tslu            & ! intent()
                             , tsld            & ! intent()
                             , fupbs           & ! intent()
                             , fdownbs         & ! intent()
                             , fnetbs          & ! intent()
                             , nsol            & ! intent()
                             , fslu            & ! intent()
                             , fsld            & ! intent()
                             , alb_toa         & ! intent()
                             , alb_tomi        & ! intent()
                             , alb_toai        & ! intent()
                             , solfx           & ! intent()
                             , tiru            & ! intent()
                             , fupbi           & ! intent()
                             , fdownbi         & ! intent()
                             , fnetbi          & ! intent()
                             , firu            ! ! intent()
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      integer                    , intent(in) :: m1
      real                       , intent(in) :: albedt
      real                       , intent(in) :: cosz
      real                       , intent(in) :: aot11
      !----- Local variables. -------------------------------------------------------------!
      real    , dimension(nlayer)             :: heati
      real    , dimension(nlayer)             :: heats
      real    , dimension(nlayer)             :: heat 
      integer                                 :: i,i1,j1
      integer                                 :: ig
      integer                                 :: j
      integer                                 :: l,k
      real                                    :: term1
      integer                                 :: count=0
      !------------------------------------------------------------------------------------!
  
      heats   =  0.0
      heati   =  0.0
     
      !----- Interpolate temperature from layer centre (t) to layer edge (tt). ------------!
      tt(1) = t_aerad(1)
      do j = 2, nvert
         tt(j) =  t_aerad(j-1) * (press(j)/p_aerad(j-1))                                   &
               ** (log(t_aerad(j)/t_aerad(j-1)) / log(p_aerad(j)/p_aerad(j-1)))
      end do
  
      !----- Water vapour (g/cm2). --------------------------------------------------------!
      do j = 2, nlayer
         rdh2o(j)   = qv_aerad(j-1) * dpg(j-1)
      end do
  
      !----- Aerosol concentrations (#/cm2). ----------------------------------------------!
      do ig = 1, ngroup
         do j = 2, nvert
            do i = 1, nrad
               caer(j,i,ig) = pc_aerad(j-1,i,ig)
            end do
         end do
      end do  
  
      !----- Surface reflectivity and emissivity. -----------------------------------------!
      do l = 1,nsolp
         rsfx(l) = albedt
         emis(l) = 0.0
      end do
      do l = nsolp+1,ntotal
         emis(l) = emisir_aerad
         rsfx(l) = 1.0 - emis(l)
      end do
     
      !----- Set wavelength limits lla and lls based on values of isl and ir. -------------!
      lla = ntotal
      lls = 1
      if (.not. isl_aerad) lls  =  nsolp+1
      if (.not. ir_aerad)  lla  =  nsolp
       
             
      !----- Compute some optical properties. ---------------------------------------------!
      call oppr(m1,aot11)

      !------------------------------------------------------------------------------------!
      !     In case infrared calculations are required, then we also compute the Planck's  !
      ! function.                                                                          !
      !------------------------------------------------------------------------------------!
      if (ir_aerad) call oppr1(m1)

      !------------------------------------------------------------------------------------!
      !     If no infrared scattering exists, then we set the index to the number of       !
      ! solar intervals.                                                                   !
      !------------------------------------------------------------------------------------!
      if(irs == 0) lla = nsolp

      !------------------------------------------------------------------------------------!
      !     If either solar or infrared scattering calculations are required, we call the  ! 
      ! two stream code and find the solution.                                             !
      !------------------------------------------------------------------------------------!
      call twostr(m1)
      call add(m1,cosz)

      !------------------------------------------------------------------------------------!
      !     If infrared calculations are required then call newflux1 for a more accurate   !
      ! solution.                                                                          !
      !------------------------------------------------------------------------------------!
      if (ir_aerad) call newflux1(m1)

      !----- Calculate infrafred and solar heating rates (deg/day). -----------------------!
      if (isl_aerad) then
         do j = 1,nvert
            term1 = fdegday/(dpg(j) * gcgs)
            do l =  1,nsolp
               heats(j) = heats(j) + (fnet(l,j+1) - fnet(l,j)) * term1
            end do
         end do
      end if

      do j = 1,nvert
         if (ir_aerad) term1 = fdegday / (dpg(j) * gcgs)
         do l = nsolp+1,ntotal
            if (ir_aerad) then
               heati(j)  = heati(j)                                                        &
                         + (directu(l,j+1) - direc(l,j+1) - (directu(l,j) - direc(l,j)))   &
                         * term1
            end if
         end do
      end do

      do j = 1,nvert
         !----- Load heating rates [deg_k/s] into interface common block. -----------------!
         heat(j)        =  heats(j) + heati(j)
         heats_aerad(j) =  heats(j) / day_sec
         heati_aerad(j) =  heati(j) / day_sec
      end do

      if (any(abs(heats) > 200.)) then
         write (unit=*,fmt='(77a)') ('-',j=1,77)
         write (unit=*,fmt='(7(a,1x))') '    J','    L','       HEATS','        FNET'      &
                                                       ,'         DPG','     FDEGDAY'      &
                                                       ,'       TERM1'
         write (unit=*,fmt='(77a)') ('-',j=1,77)
         do l = 1, nsolp
            do j = 1, nvert
               term1 = fdegday/(dpg(j) * gcgs)
               write (unit=*,fmt='(2(i5,1x),5(es12.5,1x))') j,l,heats(j),fnet(l,j),dpg(j)  &
                                                               ,fdegday,term1
            end do
         end do
         write (unit=*,fmt='(77a)') ('-',j=1,77)
      end if

      solnet      = 0.0
      soldowntoa  = 0.0
      soluptoa    = 0.0
      if (isl_aerad) then
         do l = 1,nsolp
            solnet     = solnet     - fnet(l,nlayer)
            soldowntoa = soldowntoa + direc(l,1)
            soluptoa   = soluptoa   - fnet(l,1) - direc(l,1)
         end do
      end if
    
      xirdown = 0.0
      xirup   = 0.0
      if (ir_aerad) then
         do l = nsolp+1,ntotal
            xirdown = xirdown + direc(l,nlayer)
            xirup   = xirup   + directu(l,1)
         end do
      end if

      return
   end subroutine radtran
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !      This subroutine calculates optical properties such as single scattering albedo,  !
   ! asymmetry parameter, etc.                                                             !
   !                                                                                       !
   !                                                                                       !
   !    w0(nwave,nlayer)   : single scattering albedo *** delta scaled ***                 !
   !    g0(nwave,nlayer)   : asymmetry parameter *** delta scaled ***                      !
   !    opd(nwave,nlayer)  : cumulative optical depth *** delta scaled ***                 !
   !    sfl(nwave)         : solar flux                                                    !
   !    uw0(nwave,nlayer)  : unscaled single scattering albedo                             !
   !    ug0(nwave,nlayer)  : unscaled asymmetry parameter                                  !
   !    utaul(nwave,nlayer): unscaled optical depth of layer                               !
   !                                                                                       !
   !     Assume that p is same on all sigma levels. if psurface varies a lot, then we may  !
   ! need to calculate tauo3, tauco2, tauo2 for each point.                                !
   !                                                                                       !
   !     NOTE : The top layer is a dummy.  It contains a defined gas amount.  Different    !
   ! models will require different treatment of this.                                      !
   !                                                                                       !
   !     Calculate total optical depth including gases.  Then given the aerosol optical    !
   ! depths and cloud optical depths, calculate final optical properties. we use a delta   !
   ! two stream approach to find w0, single scattering albedo, g0, asymmmetry parameter,   !
   ! taul, layer optical depth, opd, cumulative optical depth to base of layer.            !
   !---------------------------------------------------------------------------------------!
   subroutine oppr(m1,aot11)
      use rconstants , only : gcgs          ! ! intent(in)
      use mem_globrad, only : nlayer        & ! intent(something)
                            , nwave         & ! intent(something)
                            , ngroup        & ! intent(something)
                            , nrad          & ! intent(something)
                            , nirp          & ! intent(something)
                            , xsecta        & ! intent(something)
                            , rdqext        & ! intent(something)
                            , qscat         & ! intent(something)
                            , qbrqs         & ! intent(something)
                            , ntotal        & ! intent(something)
                            , nprob         & ! intent(something)
                            , roundoff      & ! intent(something)
                            , roundoff8     & ! intent(something)
                            , ptop          & ! intent(something)
                            , p             & ! intent(something)
                            , q             & ! intent(something)
                            , nsolp         & ! intent(something)
                            , contnm        & ! intent(something)
                            , uw0           & ! intent(something)
                            , ug0           & ! intent(something)
                            , ir            & ! intent(something)
                            , ngauss        & ! intent(something)
                            , gangle        & ! intent(something)
                            , ta            & ! intent(something)
                            , tb            & ! intent(something)
                            , wa            & ! intent(something)
                            , wb            & ! intent(something)
                            , ga            & ! intent(something)
                            , gb            & ! intent(something)
                            , tia           & ! intent(something)
                            , tib           & ! intent(something)
                            , wia           & ! intent(something)
                            , wib           & ! intent(something)
                            , gia           & ! intent(something)
                            , gib           & ! intent(something)
                            , alpha         & ! intent(something)
                            , gama          & ! intent(something)
                            , caseE         & ! intent(something)
                            , caseW         & ! intent(something)
                            , caseG         & ! intent(something)
                            , wave          & ! intent(something)
                            , lmie          ! ! intent(something)
      use mem_aerad  , only : lprocopio     ! ! intent(in)
      use rconstants , only : t008          ! ! intent(in)

      implicit none

      !----- Arguments. -------------------------------------------------------------------!
      integer                                    , intent(in) :: m1
      real                                       , intent(in) :: aot11
      !----- Local variables. -------------------------------------------------------------!
      real(kind=8), dimension(ntotal,nlayer)                  :: taua8
      real(kind=8), dimension(ntotal,nlayer)                  :: taus8
      real(kind=8), dimension(ntotal,nlayer)                  :: g018
      real(kind=8), dimension(ntotal,nlayer)                  :: wol8
      real(kind=8), dimension(ntotal,nlayer)                  :: gol8
      real(kind=8), dimension(ntotal,nlayer)                  :: tauh2o8
      real(kind=8), dimension(ntotal,nlayer)                  :: utaul8
      real(kind=8), dimension(ntotal,nlayer)                  :: taucld8
      real(kind=8), dimension(ntotal,nlayer)                  :: wcld8
      real(kind=8), dimension(ntotal,nlayer)                  :: gcld8
      real(kind=8), dimension(ntotal,nlayer)                  :: wolc8
      real(kind=8), dimension(ntotal,nlayer)                  :: woice8
      real(kind=8), dimension(ntotal,nlayer)                  :: worain8
      real(kind=8), dimension(ntotal,nlayer)                  :: gl8
      real(kind=8), dimension(ntotal,nlayer)                  :: gice8
      real(kind=8), dimension(ntotal,nlayer)                  :: grain8
      real(kind=8), dimension(ntotal,nlayer)                  :: denc8
      real(kind=8), dimension(ntotal,nlayer)                  :: taucldlw8
      real(kind=8), dimension(ntotal,nlayer)                  :: taucldice8
      real(kind=8), dimension(ntotal,nlayer)                  :: taurain8
      real(kind=8), dimension(ntotal)                         :: wot8
      real(kind=8), dimension(ntotal)                         :: got8
      real(kind=8), dimension(nlayer)                         :: corr8
      real(kind=8), dimension(nlayer)                         :: reffi8
      real(kind=8), dimension(nwave)                          :: rdqextnew8
      real(kind=8), dimension(nwave)                          :: wonew8
      real(kind=8), dimension(nwave)                          :: gonew8
      integer                                                 :: i,i1,j1,kk
      integer                                                 :: ig
      integer                                                 :: iradgas
      integer                                                 :: j
      integer                                                 :: l
      real(kind=8)                                            :: cco8
      real(kind=8)                                            :: den8
      real(kind=8)                                            :: denom8
      real(kind=8)                                            :: fo8
      real(kind=8)                                            :: pcorr8
      real(kind=8)                                            :: qcorr8
      real(kind=8)                                            :: ttas8
      real(kind=8)                                            :: x_teste8
      integer                                                 :: idaot
      integer                                                 :: jjj
      !----- Double precision version of most variables in this subroutine. ---------------!
      real(kind=8)                                            :: rain_aerad8
      real(kind=8)                                            :: p_top8
      real(kind=8), dimension(nrad,ngroup)                    :: xsecta8
      real(kind=8), dimension(nlayer,nrad,ngroup)             :: caer8
      real(kind=8), dimension(ntotal,nlayer)                  :: tauaer8
      real(kind=8), dimension(nrad,ngroup,nwave)              :: rdqext8
      real(kind=8), dimension(nrad,ngroup,nwave)              :: qscat8
      real(kind=8), dimension(nrad,ngroup,nwave)              :: qbrqs8
      real(kind=8), dimension(ntotal)                         :: ta8
      real(kind=8), dimension(ntotal)                         :: tb8
      real(kind=8), dimension(ntotal)                         :: wa8
      real(kind=8), dimension(ntotal)                         :: wb8
      real(kind=8), dimension(ntotal)                         :: ga8
      real(kind=8), dimension(ntotal)                         :: gb8
      real(kind=8), dimension(ntotal)                         :: tia8
      real(kind=8), dimension(ntotal)                         :: tib8
      real(kind=8), dimension(ntotal)                         :: wia8
      real(kind=8), dimension(ntotal)                         :: wib8
      real(kind=8), dimension(ntotal)                         :: gia8
      real(kind=8), dimension(ntotal)                         :: gib8
      real(kind=8), dimension(ntotal)                         :: alpha8
      real(kind=8), dimension(ntotal)                         :: gama8
      real(kind=8), dimension(ngauss)                         :: gangle8
      real(kind=8), dimension(ntotal,nlayer)                  :: paray8
      real(kind=8), dimension(ntotal,nlayer)                  :: taugas8
      real(kind=8), dimension(ntotal,nlayer)                  :: pah2o8
      real(kind=8), dimension(nlayer)                         :: rdh2o8
      real(kind=8), dimension(nlayer)                         :: tt8
      real(kind=8), dimension(nirp)                           :: contnm8
      real(kind=8), dimension(m1)                             :: p_aerad8
      real(kind=8), dimension(m1)                             :: qv_aerad8
      real(kind=8), dimension(m1)                             :: t_aerad8
      real(kind=8), dimension(m1)                             :: lwl_aerad8
      real(kind=8), dimension(m1)                             :: iwl_aerad8
      real(kind=8), dimension(m1)                             :: lwp_aerad8
      real(kind=8), dimension(m1)                             :: iwp_aerad8
      real(kind=8), dimension(ntotal,nlayer)                  :: taul8
      real(kind=8), dimension(ntotal,nlayer)                  :: w08
      real(kind=8), dimension(ntotal,nlayer)                  :: g08
      real(kind=8), dimension(ntotal,nlayer)                  :: opd8
      real(kind=8), dimension(ntotal,nlayer)                  :: uopd8
      real(kind=8)                                            :: y38
      logical                                                 :: alright
      !----- External functions. ----------------------------------------------------------!
      logical, external                                       :: is_finite
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !     Initialise all local arrays.                                                   !
      !------------------------------------------------------------------------------------!
      
      !----- Copy the following values to the scratch double precision variables. ---------!
      xsecta8      = dble(xsecta    )
      
      rdqext8      = dble(rdqext    )
      qscat8       = dble(qscat     )
      qbrqs8       = dble(qbrqs     )
      ta8          = dble(ta        )
      tb8          = dble(tb        )
      tt8          = dble(tt        )
      wa8          = dble(wa        )
      wb8          = dble(wb        )
      ga8          = dble(ga        )
      gb8          = dble(gb        )
      tia8         = dble(tia       )
      tib8         = dble(tib       )
      wia8         = dble(wia       )
      wib8         = dble(wib       )
      gia8         = dble(gia       )
      gib8         = dble(gib       )
      alpha8       = dble(alpha     )
      gama8        = dble(gama      )
      gangle8      = dble(gangle    )
      paray8       = dble(paray     )
      taugas8      = dble(taugas    )
      pah2o8       = dble(pah2o     )
      rdh2o8       = dble(rdh2o     )
      p_aerad8     = dble(p_aerad   )
      qv_aerad8    = dble(qv_aerad  )
      t_aerad8     = dble(t_aerad   )
      lwl_aerad8   = dble(lwl_aerad )
      iwl_aerad8   = dble(iwl_aerad )
      lwp_aerad8   = dble(lwp_aerad )
      iwp_aerad8   = dble(iwp_aerad )
      contnm8      = dble(contnm    )
      rain_aerad8  = dble(rain_aerad)
      p_top8       = dble(p_top     )
      !------------------------------------------------------------------------------------!
      
      !------ Copy over the following global variables from mem_carma, to dp --------------!

      tauaer8      = dble(tauaer    )
      taucldlw8    = dble(taucldlw  )
      taucldice8   = dble(taucldice )
      taucld8      = dble(taucld    )
      wcld8        = dble(wcld      )
      gcld8        = dble(gcld      )
      gl8          = dble(gl        )
      gice8        = dble(gice      )
      wolc8        = dble(wolc      )
      caer8        = dble(caer      )
      woice8       = dble(woice     )
      gice8        = dble(gice      )
      taul8        = dble(taul      )
      opd8         = dble(opd       )
      uopd8        = dble(uopd      )

      corr8      = 0.d0
      denc8      = 0.d0
      reffi8     = 0.d0
      taucldlw8  = 0.d0
      taucldice8 = 0.d0
      taurain8   = 0.d0
      wolc8      = 0.d0
      woice8     = 0.d0
      worain8    = 0.d0
      gl8        = 0.d0
      gice8      = 0.d0
      grain8     = 0.d0
      taucld8    = 0.d0
      wcld8      = 0.d0
      gcld8      = 0.d0
      rdqextnew8 = 0.d0
      wonew8     = 0.d0
      gonew8     = 0.d0


      if (lprocopio .and. lmie) then
 
         idaot = max(min(int(10*((anint(10.*max(aot11,1e-20))/10.)+0.1)/2.),9),1)

         do l = 1,nwave
            rdqextnew8(l) = dble(casee(idaot,l)) 
            wonew8(l)     = dble(casew(idaot,l)) 
            gonew8(l)     = dble(caseg(idaot,l))
         end do

         taua8=0.d0
         do j=1,nlayer
            do ig = 1,ngroup
               do i = 1,nrad
                 do l = 1,ntotal
                    taua8(l,j)   = taua8(l,j)                                              &
                                 + rdqextnew8(nprob(l)) * xsecta8(i,ig) * caer8(j,i,ig)
                    tauaer8(l,j) = max(taua8(nprob(l),j),roundoff8)
                    wol8(l,j)    = wonew8(nprob(l))
                    gol8(l,j)    = gonew8(nprob(l))
                    tauaer(l,j)  = sngl(tauaer8(l,j))
                 end do
               end do
            end do
         end do

      else
         taua8 = 0.0
         taus8 = 0.0
         g018  = 0.0

         do j=1,nlayer
            do ig = 1,ngroup
               do  i = 1,nrad
                  do  l = 1,nwave
                      taua8(l,j) = taua8(l,j)                                              &
                                 + rdqext8(i,ig,l) * xsecta8(i,ig) * caer8(j,i,ig)
                      taus8(l,j) = taus8(l,j)                                              &
                                 + qscat8(i,ig,l)  * xsecta8(i,ig) * caer8(j,i,ig)
                      g018(l,j)  = g018(l,j)                                               &
                                 + qbrqs8(i,ig,l)  * xsecta8(i,ig) * caer8(j,i,ig)
                  end do
               end do
            end do
            

            do l= 1,ntotal
               tauaer8(l,j) = max(taua8(nprob(l),j),roundoff8)
               tauaer(l,j)  = sngl(tauaer8(l,j))
               wol8(l,j)    = taus8(nprob(l),j)/tauaer8(l,j)
               if (wol8(l,j) /= 0.d0) then
                  ttas8 = taus8(nprob(l),j)
               else
                  ttas8 = 1.d0
               end if
               gol8(l,j)    = g018(nprob(l),j)/ttas8
            end do
         end do
         lmie = .true.

      end if


      do j=1,nlayer
         if (xland_aerad >= .009) then
           reffi8(j) =  7.0d0 * 1.d3 * lwl_aerad8(j) + 5.5d0 
         else 
           reffi8(j) =  9.5d0 * 1.d3 * lwl_aerad8(j) + 4.0d0
         end if
         
         corr8(j) = 1.047d0 - 9.13d-5 * (tt8(j)-t008) + 2.03d-4 * (tt8(j)-t008) **2        &
                            - 1.06d-5 * (tt8(j)-t008) **3
         corr8(j) = max(corr8(j),roundoff8)  
      end do

      do j=1,nlayer
         do l= 1,ntotal        
            if (j == 1) taurain8(l,j) = 1.8d-4 * rain_aerad8 * 2.0d3
            
            taucldlw8  (l,j) = 1.d3 * lwp_aerad8(j) *(ta8(l)/reffi8(j)+tb8(l)/reffi8(j)**2)

            worain8    (l,j) = 5.5d-1
            wolc8      (l,j) = (1.d0 - wa8(l)) + wb8(l) * reffi8(j)
        
            gl8        (l,j) = ga8(l) + gb8(l) * reffi8(j)
        
            grain8     (l,j) = 9.5d-1
            denc8      (l,j) = 1.d0 / (tia8(l) + tib8(l) * 1.d3 * iwl_aerad8(j))
            taucldice8 (l,j) = corr8(j) * 1.0d3 * iwp_aerad8(j) * denc8(l,j)
            
            if (wib8(l) < roundoff8) then
               x_teste8 = 0.d0
            else
               x_teste8 = (1.0d3 * iwl_aerad8(j))**wib8(l)
            end if

            woice8(l,j) = ( 1.d0 -  wia8(l) * x_teste8)                                    &
                        * ( 1.d0 - gama8(l) * (corr8(j) - 1.d0)/ corr8(j) )
            
            if(gib8(l) < roundoff8) then
               x_teste8 = 1.d0
            else
               x_teste8 = (1.0d3 * iwl_aerad8(j))**gib8(l)
            end if
            gice8(l,j)  = (gia8(l) * x_teste8)                                             &
                        * ( 1.d0 + alpha8(l) * (corr8(j) - 1.d0)/ corr8(j) )
                           
            
            select case (l)
            case (91:113)
               taucldlw8(l,j)= 1.0d3 * lwp_aerad8(j) * ta8(l) * exp(tb8(l)* reffi8(j))

            case (114:154)
               taucldlw8(l,j) = 1.0d3 * lwp_aerad8(j) * ( ta8(l) + tb8(l)* reffi8(j))
               gl8(l,j)       = 1.d0 - ga8(l) * exp( gb8(l) * reffi8(j))

            end select
            
            taucld8(l,j) = taucldlw8(l,j) + taucldice8(l,j) + taurain8(l,j)

            if ( taucld8(l,j) > roundoff8) then
               wcld8(l,j) = ( wolc8(l,j) * taucldlw8(l,j)  + woice8(l,j) * taucldice8(l,j) &
                            + worain8(l,j) * taurain8(l,j)) / taucld8(l,j)
               gcld8(l,j) = ( wolc8(l,j) *  taucldlw8(l,j) * gl8(l,j)                      &
                            + woice8(l,j) * taucldice8(l,j)* gice8(l,j)                    &
                            + worain8(l,j) * taurain8(l,j) * grain8(l,j))                  &
                            / (wcld8(l,j) * taucld8(l,j))            
            else 
              wcld8(l,j)  = 1.d0
              gcld8(l,j)  = 0.d0
            end if
         end do
      end do

      iradgas = 1
      do  j = 1,nlayer
         kk = max( 1, j-1 )
         !---------------------------------------------------------------------------------!
         !   Bergstrom water vapor continuum fix:                                          !
         !                                                                                 !
         !    <qcorr8> is layer average water vapor mixing ratio                           !
         !    <pcorr8> is layer average pressure                                           !
         !                                                                                 !
         !   For layer 0, calculate mixing ratio [g/g] from vapor column [g/cm^2]          !
         !   and average pressure [dyne/cm^2].                                             !
         !---------------------------------------------------------------------------------!
         if (j == 1)then
            qcorr8 = rdh2o8(1) * dble(gcgs) / dble(p_top)
            pcorr8 = 5.d-1 * p_aerad8(1)
            cco8   = exp(1.8d3/t_aerad8(kk))*(qcorr8*pcorr8/2.87d0 + pcorr8/4.61d3)
         else
            qcorr8 = qv_aerad8(kk)
            pcorr8 = p_aerad8(kk)
            cco8 = exp(1.8d3/t_aerad8(kk))*(qcorr8*pcorr8/2.87d0 + pcorr8/4.61d3)
         end if
  
         do l = lls, lla
            !---- Bergstrom water vapor continuum. ----------------------------------------!
            select case (l-nsolp)
            case (:30)
               tauh2o8(l,j) = rdh2o8(j) * pah2o8(l,j)
            case (31:36)
               tauh2o8(l,j) = rdh2o8(j) * pah2o8(l,j) * cco8
            case (37:)
               tauh2o8(l,j) = rdh2o8(j) * pah2o8(l,j) + cco8 * rdh2o8(j) * contnm8(l-nsolp)
            end select

            if (iradgas == 0) then
               taul8(l,j) = max(roundoff8,tauaer8(l,j))
            else
               taul8(l,j) = max( roundoff8, tauh2o8(l,j) + taugas8(l,j) + paray8(l,j)      &
                                          + tauaer8(l,j) + taucld8(l,j))
            end if
         end do

         do l = lls,lla
            utaul8(l,j)  = taul8(l,j)
            if (iradgas == 0) then
               wot8(l) = wol8(l,j)
            else
               wot8(l)      = ( paray8(l,j) + tauaer8(l,j) * wol8(l,j)                     &
                              + taucld8(l,j) * wcld8(l,j) ) / taul8(l,j)
            end if
            wot8(l)      = min(1.d0 - roundoff8,wot8(l))

            if (iradgas == 0) then
               got8(l) = gol8(l,j)
            else
               denom8       = paray8(l,j) + taucld8(l,j) * wcld8(l,j)                      &
                                          + tauaer8(l,j) * wol8(l,j)
               if (denom8 > roundoff8 ) then
                  got8(l) = ( wcld8(l,j) * gcld8(l,j) * taucld8(l,j)                       &
                            + gol8(l,j)  * wol8(l,j)  * tauaer8(l,j) ) / denom8
                  got8(l) = max(roundoff8, got8(l))
               else
                  got8(l) = roundoff8 !----- Used to be zero... ---------------------------!
               end if
            end if
         end do

         do l = lls, lla
            fo8        = got8(l) * got8(l)
            den8       = 1.d0 - wot8(l) * fo8
            taul8(l,j) = taul8(l,j) * den8
            w08  (l,j) = (1.d0 - fo8) * wot8(l) / den8
            g08  (l,j) = got8(l) / (1.d0 + got8(l))
            opd8 (l,j) = 0.d0
            opd8 (l,j) = opd8(l,kk) + taul8(l,j)
            uopd8(l,j) = 0.d0
            uopd8(l,j) = uopd8 (l,kk) + utaul8(l,j)
            
            !----- Copying the double precision arrays to the single precision. -----------!
            taul (l,j) = sngl(taul8(l,j))
            w0   (l,j) = sngl(w08  (l,j))
            g0   (l,j) = sngl(g08  (l,j))
            opd  (l,j) = sngl(opd8 (l,j))
            uopd (l,j) = sngl(uopd8(l,j))
         end do
      end do
  
      if (ir_aerad) then
         do j = 1,nlayer
            do i = 1,ngauss
               do l = lls, lla
                  y38        = exp(-min(4.6d1,taul8(l,j)/gangle8(i)))
                  y3 (l,i,j) = sngl(y38)
               end do
            end do
         end do
      end if
      
      !----- Sanity check. ----------------------------------------------------------------!
      alright = .true.
      okloop1: do j=1,nlayer
         okloop2: do l=lls,lla
            alright = is_finite(taul(l,j))
            if (.not. alright) exit okloop1
         end do okloop2
      end do okloop1

      if (.not. alright) then
         write (unit=*,fmt='(12(a,1x))') '   KK','   IB','        TAUL','       TAUL8'     &
                                                        ,'    TAUCLDI8','       CORR8'     &
                                                        ,'        TIA8','        TIB8'     &
                                                        ,'       DENC8','         TT8'     &
                                                        ,'        IWP8','        IWL8'
         do j=1,nlayer
            do l=lls,lla
               write (unit=*,fmt='(2(i5,1x),10(es12.5,1x))')                               &
                  j,l,taul(l,j),taul8(l,j),taucldice8(l,j),corr8(j)                        &
                     ,tia8(l),tib8(l),denc8(l,j),tt8(j),iwl_aerad8(j),iwp_aerad8(j)
            end do
         end do
      end if
      !------------------------------------------------------------------------------------!

      !------ Copy the dp values back to their global singles -----------------------------!

      tauaer      = sngl(tauaer8    )
      taucldlw    = sngl(taucldlw8  )
      taucldice   = sngl(taucldice8 )
      taucld      = sngl(taucld8    )
      wcld        = sngl(wcld8      )
      gcld        = sngl(gcld8      )
      gl          = sngl(gl8        )
      gice        = sngl(gice8      )
      wolc        = sngl(wolc8      )
      woice       = sngl(woice8     )
      gice        = sngl(gice8      )
      taul        = sngl(taul8      )
      opd         = sngl(opd8       )
      uopd        = sngl(uopd8      )
      !------------------------------------------------------------------------------------!

      return
   end subroutine oppr 
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine computes the Planck function and its derivative, at the ground as !
   ! well as the other levels.                                                             !
   !---------------------------------------------------------------------------------------!
   subroutine oppr1(m1)
      use mem_globrad, only : ntotal            & ! intent(in)
                            , tgrnd             & ! intent(in)
                            , nlow              & ! intent(in)
                            , nirp              & ! intent(in)
                            , planck            & ! intent(in)
                            , ltemp             & ! intent(in)
                            , nsolp             & ! intent(in)
                            , weight            & ! intent(in)
                            , iblackbody_above  & ! intent(in)
                            , t_above           & ! intent(in)
                            , nlayer            & ! intent(in)
                            , ncount            & ! intent(in)
                            , roundoff          ! ! intent(in)

      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      integer, intent(in) :: m1
      !----- Local variables. -------------------------------------------------------------!
      integer, dimension(nlayer)        :: it1
      integer                           :: itg
      integer                           :: itp
      integer                           :: i
      integer                           :: j
      integer                           :: kindex
      integer                           :: l
      integer                           :: i1
      integer                           :: j1
      real   , dimension(ntotal)        :: pltemp1
      real   , dimension(ntotal,nlayer) :: ptemp2
      !------------------------------------------------------------------------------------!

      !----- Calculate the wavelength dependent Planck's function at the ground. ----------!
      itg = anint(100. * t_surf) - nlow
      do i = 1,nirp
         pltemp1(i) = planck(ltemp(i),itg)
      end do
      do l = nsolp+1,ntotal
         ptempg(l) = pltemp1(l-nsolp)*weight(l)
      end do
      if (iblackbody_above) then
         !---------------------------------------------------------------------------------!
         !     Calculate the wavelength dependent planck function at the top of the model. !
         !---------------------------------------------------------------------------------!
         itp = anint(100. * tabove_aerad) - nlow
         do i = 1,nirp
            pltemp1(i) = planck(ltemp(i),itp)
         end do
         do l = nsolp+1,ntotal
            ptempt(l) = pltemp1(l-nsolp)*weight(l)
         end do
      end if

      do j = 1,nlayer
         it1(j) = anint(100.*tt(j)) - nlow
      end do
      do j = 1,nlayer
         do i = 1,nirp
            ptemp2(i,j)=planck(ltemp(i),it1(j))
         end do
      end do

      !------------------------------------------------------------------------------------!
      !      kindex makes the top layer isothermal.  Using kindex, find Planck's function  !
      ! at the bottom of each layer.  Please notice that when slope is set to 0, then      !
      ! there will be isothermal layers with tt(j) corresponding to average temperature    !
      ! of layer.  In this case, tt(nlayer) must be set to tgrnd.                          !
      !------------------------------------------------------------------------------------!
      do j = 1,nlayer
         kindex = max(1,j-1)
         do l = nsolp+1, ntotal
            ptemp(l,j) = ptemp2(l-nsolp,j) * weight(l)
            slope(l,j) = (ptemp(l,j)-ptemp(l,kindex)) / taul(l,j)
            if (taul(l,j) <= roundoff) then
               slope(l,j) = 0.
            else 
               slope(l,j) = (ptemp(l,j)-ptemp(l,kindex)) / taul(l,j)
            end if
         end do
      end do

      return
   end subroutine oppr1
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
  SUBROUTINE twostr(m1)
    !
    !	 ******************************************************************
    !	 *  Purpose		:  Defines matrix properties and sets up  *
    !	 *			   matrix coefficients that do not depend *
    !	 *			   on zenith angle or temperature.	  *
    !	 *  Subroutines Called  :  None 				  *
    !	 *  Input		:  W0, G0				  *
    !	 *  Output		:  B1, B2, GAMI, ACON, EL1, AF, ETC	  *
    !	 * ****************************************************************
    !
    use rconstants , only: srthree, twopi
    USE mem_globrad, ONLY: nsolp,nlayer,jn,jdble,irs,ntotal
  
    IMPLICIT NONE
  
    INTEGER,INTENT(IN) :: m1
    INTEGER	   :: j
    INTEGER	   :: jd
    INTEGER	   :: l
    REAL,PARAMETER :: two = 2.d0
    
      
    DO  l    =  1,ntotal !lls(i1,j1),lla(i1,j1)
  	IF(isl_aerad .OR. irs .NE. 0 ) THEN
  	  IF( l>=lls .AND. l<=lla) THEN
  	    IF(l <= nsolp ) THEN
  	      u1i(l) = srthree
  	    ELSE
  	      u1i(l) = two
  	    END IF
  	    !u1s(l)  =  twopi/u1i(l)
  	  END IF
  	END IF
    END DO
    !
    !	   here we define layer properties following general scheme
    !	   of meador and weavor. then we set up layer properties
    !	   needed for matrix.
    !
    DO  j =  1,nlayer
      DO  l=  1,ntotal
  	  IF(isl_aerad  .OR. irs .NE. 0 ) THEN
  	    IF( l>=lls .AND. l<=lla) THEN
  	      !these are for two stream and hemispheric means
  	      b1(l,j)   =  0.5*u1i(l)*(2.-w0(l,j)*(1. + g0(l,j)))
  	      b2(l,j)   =  0.5*u1i(l)*w0(l,j)*(1. - g0(l,j))
  	      ak(l,j)   =  SQRT(ABS(b1(l,j)**2 - b2(l,j)**2))
  	      gami(l,j)  =  b2(l,j)/(b1(l,j) + ak(l,j))
  	      ee1(l,j)   =  EXP(-ak(l,j)*taul(l,j))
  	      el1(l,j)   =  1.0 + gami(l,j) *ee1(l,j)
  	      em1(l,j)   =  1.0 - gami(l,j) * ee1(l,j)
  	      el2(l,j)   =  gami(l,j) + ee1(l,j)
  	      em2(l,j)   =  gami(l,j) - ee1(l,j)
  	    END IF
  	  END IF
      END DO
    END DO
    !
    !	  we seek to solve ax(l-1)+bx(l)+ex(l+1) = d.
    !	  l=2n for even l, l=n+1 for odd l. the mean intensity (tmi/4pi)
    !	  and the net flux (fnet) are related to x's as noted in add.
    !	  first we set up the coefficients that are independent of solar
    !	  angle or temparature: a(i),b(i),e(i). d(i) is defined in add.
    !
    j=  0
    DO  jd=  2,jn,2
      j=  j + 1
      DO  l=  1,ntotal
  	  IF(isl_aerad .OR. irs .NE. 0 ) THEN
  	    IF( l>=lls .AND. l<=lla) THEN
  	      !here are the even matrix elements
  	      af(l,jd)   =  em1(l,j+1)*el1(l,j)- &
  				  em2(l,j+1)*el2(l,j)
  	      bf(l,jd)   =  em1(l,j+1)* em1(l,j)- &
  				  em2(l,j+1)*em2(l,j)
  	      ef(l,jd)   =  el1(l,j+1)*em2(l,j+1) - &
  				  el2(l,j+1)*em1(l,j+1)
  	      !here are the odd matrix elements except for the top.
  	      af(l,jd+1) =  em1(l,j)*el2(l,j)- &
  				  el1(l,j)*em2(l,j)
  	      bf(l,jd+1) =  el1(l,j+1)*el1(l,j) - &
  				  el2(l,j+1)*el2(l,j)
  	      ef(l,jd+1) =  el2(l,j)*em2(l,j+1)- &
  				el1(l,j)*em1(l,j+1)
  	    END IF
  	  END IF
      END DO
    END DO
    !
    !	  HERE ARE THE TOP AND BOTTOM BOUNDARY CONDITIONS AS WELL AS THE
    !	  BEGINNING OF THE TRIDIAGONAL SOLUTION DEFINITIONS. I ASSUME
    !	  NO DIFFUSE RADIATION IS INCIDENT AT UPPER BOUNDARY.
    !
    DO  l=  1,ntotal
  	IF(isl_aerad .OR. irs .NE. 0 ) THEN
  	  IF( l>=lls .AND. l<=lla) THEN
  	    af(l,1)    = 0.0
  	    bf(l,1)    = el1(l,1)
  	    ef(l,1)    = -em1(l,1)
  	    af(l,jdble) = el1(l,nlayer)-rsfx(l)*el2(l,nlayer)
  	    bf(l,jdble) = em1(l,nlayer)-rsfx(l)*em2(l,nlayer)
  	    ef(l,jdble) = 0.0
  	  END IF
  	END IF
    END DO
  
  END SUBROUTINE twostr
  
  SUBROUTINE add(m1,cosz)
    use rconstants , only: srthree, twopi, pi1
    USE mem_globrad, ONLY: isl,u0,nlayer,nsolp,sol,roundoff, &
  			   irs,ntotal,u1s,emis,jn, &
  			   jdble,ndbl
    
    IMPLICIT NONE 
   
    !	  THIS SUBROUTINE FORMS THE MATRIX FOR THE MULTIPLE LAYERS AND
    !	  USES A TRIDIAGONAL ROUTINE TO FIND RADIATION IN THE ENTIRE
    !	  ATMOSPHERE.
   
    !	  ******************************
    !	  *   CALCULATIONS FOR SOLAR   *
    !	  ******************************
    INTEGER,INTENT(IN) :: m1
    REAL,INTENT(IN) :: cosz  
    INTEGER :: j,kk
    INTEGER :: jd
    INTEGER :: kindex
    INTEGER :: l
    REAL    :: b4
    REAL    :: c1
    REAL    :: c2
    REAL    :: cm1
    REAL    :: cp1
    REAL    :: du0
    REAL    :: x
    REAL    :: x2
    REAL    :: x3
    REAL    :: x4
    REAL,DIMENSION(nsolp,nlayer) :: direct,el3,ee3,cm
    REAL,DIMENSION(nsolp) :: sfcs
    REAL,DIMENSION(ntotal,ndbl) :: df,as,ds,xk
    INTEGER :: i1,j1
  
    DO  j	 =  1,nlayer
      kk = MAX( 1, j-1 )
      DO  l    =  1,nsolp
  	  du0=1./cosz
  	  IF(isl_aerad)  THEN
  	    b3(l,j)     =  0.5*(1.-srthree*g0(l,j)*cosz)
  	    b4         =  1. - b3(l,j)
  	    x2         =  taul(l,j)*du0
  	    ee3(l,j)   =  EXP(-x2)
  	    x3         =  opd(l,j)*du0
  	    el3(l,j)   =  EXP(-x3)*sol(l)
  	    direct(l,j) =  cosz*el3(l,j)
  	    c1         =  b1(l,j) - du0
  	    IF( ABS(c1) < roundoff ) c1 = sign(roundoff,c1)
  	    c2         =  ak(l,j)*ak(l,j) - du0*du0
  	    IF( ABS(c2) <= roundoff ) c2 = roundoff
  	    cp1        =  w0(l,j)*(b3(l,j)*c1+b4*b2(l,j))/c2
  	    cpb(l,j)    =  cp1 * el3(l,j)
  	    IF( j /= 1 ) THEN
  	      x4 = el3(l,kk)
  	    ELSE
  	      x4 = sol(l)
  	    END IF
  	    cp(l,j)     =  cp1 * x4
  	    cm1        =  ( cp1*b2(l,j) + w0(l,j)*b4 )/c1
  	    cmb(l,j)    =  cm1 * el3(l,j)
  	    cm(l,j)    =  cm1 * x4
  	  END IF
      END DO
    END DO
    !	     CALCULATE SFCS, THE SOURCE AT THE BOTTOM.
    IF(isl_aerad)  THEN
       DO l=  1,nsolp
  	  sfcs(l)=  direct(l,nlayer) * rsfx(l)
       END DO
    END IF
   
    !	  ******************************
    !	  * CALCULATIONS FOR INFRARED. *
    !	  ******************************
    DO  j= 1,nlayer
      DO  l = nsolp+1,ntotal
  	  IF(irs /= 0)  THEN
  	      kindex = MAX(1,j-1)
  	      b3(l,j)     = 1.0/(b1(l,j)+b2(l,j))
  	      cp(l,j)     = (ptemp(l,kindex)+slope(l,j)* &
  				   b3(l,j))*(twopi/u1i(l))
  	      cpb(l,j)    = cp(l,j) + slope(l,j)* &
  				  taul(l,j)*(twopi/u1i(l))
  	      cm(l,j)     = (ptemp(l,kindex)-slope(l,j)* &
  				   b3(l,j))*(twopi/u1i(l))
  	      cmb(l,j)    = cm(l,j) + slope(l,j)* &
  				  taul(l,j)*(twopi/u1i(l))
  	      el3(l,j)    = 0.0
  	      direct(l,j) = 0.0
  	      ee3(l,j)    = 0.0
  	  END IF
      END DO
    END DO
    
    IF (irs /= 0)  THEN
       DO  l= nsolp+1,ntotal
          sfcs(l)= emis(l)*ptempg(l)*pi1
       END DO
    END IF
   
    j=  0
    DO  jd=  2,jn,2
     j=  j + 1
     DO  l=1,ntotal
  	 IF(isl_aerad .OR. irs .NE. 0 ) THEN
  	   IF(l>=lls .AND. l<=lla) THEN
  	 !	    HERE ARE THE EVEN MATRIX ELEMENTS
  	   df(l,jd) = (cp(l,j+1) - cpb(l,j))*em1(l,j+1) -  &
  		(cm(l,j+1) - cmb(l,j))*em2(l,j+1)
  	 !	    HERE ARE THE ODD MATRIX ELEMENTS EXCEPT FOR THE TOP.
  	   df(l,jd+1) =  el2(l,j) * (cp(l,j+1)-cpb(l,j)) +  &
  		el1(l,j) * (cmb(l,j) - cm(l,j+1))
  	    END IF
  	  END IF
      END DO
    END DO
   
    !	  HERE ARE THE TOP AND BOTTOM BOUNDARY CONDITIONS AS WELL AS THE
    !	  BEGINNING OF THE TRIDIAGONAL SOLUTION DEFINITIONS. I ASSUME NO
    !	  DIFFUSE RADIATION IS INCIDENT AT THE TOP.
    DO  l=1,ntotal
  	IF(isl_aerad .OR. irs .NE. 0 ) THEN
  	  IF(l>=lls.AND. l<=lla) THEN
  	    df(l,1)   = -cm(l,1)
  	    df(l,jdble) = sfcs(l)+rsfx(l)*cmb(l,nlayer)- &
  			  cpb(l,nlayer)
  	    ds(l,jdble) = df(l,jdble)/bf(l,jdble)
  	    as(l,jdble) = af(l,jdble)/bf(l,jdble)
  	  END IF
  	END IF
    END DO
   
    !	  ********************************************
    !	  *	WE SOLVE THE TRIDIAGONAL EQUATIONS   *
    !	  ********************************************
   
    DO  j = 2, jdble
      DO  l=1,ntotal
  	  IF(isl_aerad .OR. irs .NE. 0 ) THEN
  	    IF(l>=lls .AND. l<=lla) THEN
  	      x  = 1./(bf(l,jdble+1-j) - ef(l,jdble+1-j)* &
  			  as(l,jdble+2-j))
  	      as(l,jdble+1-j) = af(l,jdble+1-j)*x
  	      ds(l,jdble+1-j) = (df(l,jdble+1-j) - &
  			  ef(l,jdble+1-j) *ds(l,jdble+2-j))*x
  	    END IF
  	  END IF
      END DO
    END DO
   
    DO  l=1,ntotal
  	IF(isl_aerad .OR. irs .NE. 0 ) THEN
  	  IF(l>=lls .AND. l<=lla) THEN
  	    xk(l,1)    = ds(l,1)
  	  END IF 
  	END IF
    END DO
  
    DO  j	= 2, jdble
      DO  l=1,ntotal
  	  IF(isl_aerad .OR. irs .NE. 0 ) THEN
  	    IF(l>=lls .AND. l<=lla) THEN
  	      xk(l,j) = ds(l,j) - as(l,j)*xk(l,j-1)
  	    END IF 
  	  END IF
      END DO
    END DO
   
    !  ***************************************************************
    !	  CALCULATE LAYER COEFFICIENTS, NET FLUX AND MEAN INTENSITY
    !  ***************************************************************
      
     DO j = 1,nlayer
       DO  l=1,ntotal
  	   IF(isl_aerad .OR. irs .NE. 0 ) THEN
  	     IF(l>=lls .AND. l<=lla) THEN
  	       ck1(l,j)   = xk(l,2*j-1)
  	       ck2(l,j)   = xk(l,2*j)
	       
  	       fnet(l,j)  = ck1(l,j)  *( el1(l,j) &
  				 -el2(l,j)) + ck2(l,j) * &
  				    ( em1(l,j)-em2(l,j) ) + &
  				    cpb(l,j) - cmb(l,j) - direct(l,j)
 	       tmi(l,j)   =  el3(l,j) + u1i(l) *(ck1(l,j)  *  &
  			  ( el1(l,j) + el2(l,j))   + ck2(l,j) * &
  			  ( em1(l,j)+em2(l,j) ) +  cpb(l,j) + &
  			  cmb(l,j) )
  	     END IF 
  	   END IF
       END DO
     END DO
   
    END SUBROUTINE add
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine calculates upward and downward intensities and fluxes using Gauss !
   ! quadrature angles and weights.                                                        !
   !---------------------------------------------------------------------------------------!
   subroutine newflux1(m1)
      use rconstants  , only : twopi             ! ! intent(in)
      use mem_globrad , only : ntotal            & ! intent(in)
                             , ngauss            & ! intent(in)
                             , nlayer            & ! intent(in)
                             , nsolp             & ! intent(in)
                             , irs               & ! intent(in)
                             , gangle            & ! intent(in)
                             , iblackbody_above  & ! intent(in)
                             , gratio            & ! intent(in)
                             , gweight           & ! intent(in)
                             , emis              ! ! intent(in)
      use node_mod    , only : mynum             ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      integer                                  , intent(in) :: m1
      real    , dimension(ntotal,ngauss,nlayer)             :: y1
      real    , dimension(ntotal,ngauss,nlayer)             :: y2
      real    , dimension(ntotal,ngauss,nlayer)             :: y4
      real    , dimension(ntotal,ngauss,nlayer)             :: y8
      real    , dimension(ntotal,ngauss,nlayer)             :: dintent
      real    , dimension(ntotal,ngauss,nlayer)             :: uintent
      real    , dimension(ntotal,nlayer)                    :: a1
      real    , dimension(ntotal,nlayer)                    :: a2
      real    , dimension(ntotal,nlayer)                    :: a3
      real    , dimension(ntotal,nlayer)                    :: a4
      real    , dimension(ntotal,nlayer)                    :: a7
      real    , dimension(ntotal,nlayer)                    :: y5
      integer                                               :: i
      integer                                               :: j
      integer                                               :: kindex
      integer                                               :: l
      integer                                               :: m
      real                                                  :: ckm
      real                                                  :: ckp
      real                                                  :: x4
      real                                                  :: ya
      real                                                  :: yb
      integer                                               :: i1
      integer                                               :: j1
      logical                                               :: alright
      !----- External functions. ----------------------------------------------------------!
      logical                                  , external   :: is_finite
      !------------------------------------------------------------------------------------!

      do j =  1,nlayer
         kindex = max(1, j-1)
         !----- First the non-scattering coefficients. ------------------------------------!
         do l = nsolp+1,ntotal
            a3(l,j) = ptemp(l,kindex) * twopi
            a4(l,j) = twopi * slope(l,j)
            a7(l,j) = a3(l,j)
            y5(l,j) = a4(l,j) * taul(l,j)
         end do
         !----- Then the scattering coefficients. -----------------------------------------!
         do l = nsolp+1,ntotal
            if (irs /= 0) then
               x4      = slope(l,j) * (twopi*b3(l,j)-(twopi/u1i(l)))
               a1(l,j) = u1i(l) - ak(l,j)
               a2(l,j) = gami(l,j) * (ak(l,j)+u1i(l))
               a3(l,j) = a3(l,j) + x4
               a7(l,j) = a7(l,j) - x4
            end if
         end do
      end do

      !----- Calculations for all Gaussian points. ----------------------------------------!
      do j=  1,nlayer
         do i=  1,ngauss
            !----- No scattering. ---------------------------------------------------------!
            do l =  nsolp+1,ntotal
               y1(l,i,j)  =  0.0
               y2(l,i,j)  =  0.0
               y4(l,i,j)  =  a7(l,j) - a4(l,j) * gangle(i)
               y8(l,i,j)  =  a3(l,j) + a4(l,j) * gangle(i)
            end do
            !----- Scattering. ------------------------------------------------------------!
            do l = nsolp+1, ntotal
               if (irs /= 0) then
                  ya  = a1(l,j) * (y3(l,i,j)-ee1(l,j)) / (ak(l,j)*gangle(i)-1.)
                  yb  = a2(l,j) * (1.- ee1(l,j)*y3(l,i,j)) / (ak(l,j)*gangle(i)+1.)
                  ckp = ck1(l,j) + ck2(l,j)
                  ckm = ck1(l,j) - ck2(l,j)
                  y1(l,i,j) = ckp * yb + ckm * ya
                  y2(l,i,j) = ckp * ya + ckm * yb
               end if
            end do
         end do
      end do

      do j = 1,nlayer
         do l = nsolp+1,ntotal
            direc(l,j)   = 0.0
            directu(l,j) = 0.0
         end do
      end do

      !------------------------------------------------------------------------------------!
      !     Direc is downward flux, and directu is upward flux.  Here we calculate dintent !
      ! the downward intensity, and direc the downward flux.                               !
      !------------------------------------------------------------------------------------!
      do i = 1,ngauss
         do l = nsolp+1,ntotal
            if (iblackbody_above) then
               dintent(l,i,1) = ptempt(l) * y3(l,i,1) * twopi                              &
                              + y1(l,i,1) + (1.-y3(l,i,1)) * y4(l,i,1)
            else
               dintent(l,i,1) = (1. - y3(l,i,1)) * y4(l,i,1) + y1(l,i,1)
            end if
            direc(l,1) = direc(l,1) + dintent(l,i,1) * gweight(i)
         end do
      end do

      !------------------------------------------------------------------------------------!
      !    Dintent is downward intensity * twopi.  Direc is the downward flux.             !
      !------------------------------------------------------------------------------------!
      do j= 2,nlayer
         do i = 1,ngauss
            do l = nsolp+1,ntotal
               dintent(l,i,j) = dintent(l,i,j-1) * y3(l,i,j) + y1(l,i,j) + y5(l,j)         &
                              + (1.-y3(l,i,j)) * y4(l,i,j)
               direc(l,j)     = direc(l,j) + dintent(l,i,j) * gweight(i)
            end do
         end do
      end do

      !------------------------------------------------------------------------------------!
      !    Uintent is the upward intensity * twopi.  Directu is the upward flux.           !
      !    Assume that the reflectivity is Lambert.                                        !
      !------------------------------------------------------------------------------------!
      do i = 1,ngauss
         do l = nsolp+1,ntotal
            uintent(l,i,nlayer) = ptempg(l) * emis(l) * twopi                              &
                                + 2. * rsfx(l) * direc(l,nlayer)
            directu(l,nlayer)   = directu(l,nlayer) + uintent(l,i,nlayer)*gweight(i)
         end do
      end do

      do m= 2,nlayer
         j = nlayer-m+1
         do i = 1,ngauss
            do l = nsolp+1,ntotal
               uintent(l,i,j) = (uintent(l,i,j+1)-y5(l,j+1)) * y3(l,i,j+1) + y2(l,i,j+1)   &
                              + (1.-y3(l,i,j+1))*y8(l,i,j+1)
               directu(l,j)   = directu(l,j) + gweight(i) * uintent(l,i,j)
            end do
         end do
      end do
      
      !----- Sanity check. ----------------------------------------------------------------!
      alright = .true.
      okloop1: do j=1, nlayer
         okloop2: do l=nsolp+1,ntotal
            alright = is_finite(direc(l,j)) .and. is_finite(taul(l,j))
            if (.not. alright) exit okloop1
         end do okloop2
      end do okloop1
      if (.not. alright) then
         write (unit=*,fmt='(13(a,1x))')       '   kk',       '   ig',       '   ib'       &
                                       ,'      ptempt','          y1','          y3'       &
                                       ,'          y4','          y5','       slope'       &
                                       ,'          a4','        taul','     dintent'       &
                                       ,'       direc'
         write (unit=*,fmt='(148a)') ('-',i=1,148)
         do j=1,nlayer
            do i=1,ngauss
               do l=nsolp+1,ntotal
                  write (unit=*,fmt='(3(i5,1x),10(es12.5,1x))')                            &
                         j,i,l,ptempt(l),y1(l,i,j),y3(l,i,j),y4(l,i,j),y5(l,j)             &
                        ,slope(l,j),a4(l,j),taul(l,j),dintent(l,i,j),direc(l,j)
               end do
            end do
         end do
         write (unit=*,fmt='(148a)') ('-',i=1,148)
      end if

      return
   end subroutine newflux1
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!


  SUBROUTINE plnk(e,t1,d)
    
    !	  ******************************************************
    !	  *  Purpose		 :  Calculate Planck Function  *
    !	  *  Subroutines Called  :  None		       *
    !	  *  Input		 :  WAVE, NCOUNT	       *
    !	  *  Output		 :  PLANCK		       *
    !	  * ****************************************************
   
    !  THIS SUBROUTINE COMPUTES THE INTEGRAL OF THE PLANCK FUNCTION BETWEEN
    !  ZERO AND THE SPECIFIED VALUE OF LAMBDA.  THUS (USING XL AS LAMBDA)
    !  WE WANT TO INTEGRATE
    !  R = INTEGRAL(XL=0 TO XL=XLSPEC) ( C1*XL**-5* / (EXP(C2/XL*T)-1) )*DXL
    !  SUBSTITUTING U=C2/(XL*T), THE INTEGRAL BECOMES
    !  R = A CONSTANT TIMES INTEGRAL (USPEC TO INFINITY) OF
    !		 ( U**3 / (EXP(U) - 1) )*DU
    !  THE APPROXIMATIONS SHOWN HERE ARE ON PAGE 998 OF ABRAMOWITZ AND SEGUN
    !  UNDER THE HEADING OF DEBYE FUNCTIONS.  C2 IS THE PRODUCT OF PLANCK'S
    !  CONSTANT AND THE SPEED OF LIGHT DIVIDED BY BOLTZMANN'S CONSTANT.
    !  C2 = 14390 WHEN LAMBDA IS IN MICRONS.
    !  THE FACTOR 0.15399 IS THE RECIPROCAL OF SIX TIMES
    !  THE SUM OF (1/N**2) FOR ALL N FROM ONE TO INFINITY.  IT IS CHOSEN TO
    !  NORMALIZE THE INTEGRAL TO A MAXIMUM VALUE OF UNITY.
    !  RADIATION IN REAL UNITS IS OBTAINED BY MULTIPLYING THE INTEGRAL BY
    !  THE STEFAN-BOLTZMANN CONSTANT TIMES T**4.
    IMPLICIT NONE
   
    REAL				     :: e
    REAL, INTENT(IN)			     :: t1
    REAL, INTENT(OUT)			     :: d
    REAL :: am(5)
    REAL :: v1,a
    INTEGER :: m
   
    d		 =   0.0
    v1  	 =   e/t1
   
    IF (v1 <= 1.) THEN
      d 	=  1.0 - 0.15399*v1**3 *  &
  	  (1./3.-v1/8. + v1**2/60. - v1**4/5040. +  &
  	  v1**6/272160. - v1**8/13305600	 )
    END IF
   
    IF ( v1 > 1. .AND. v1 <= 50.) THEN
      DO  m   =  1,5
  	a	=  FLOAT(m)*v1
  	am(m)	=  0.15399 * EXP(-a)/m**4 * (((a+3.)*a+6.)*a+6.)
      END DO
   
      d 	 =  am(1)+am(2)+am(3)+am(4)+am(5)
    END IF
   
    d		  =  d*t1**4
   
  END SUBROUTINE plnk
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine computes mie scattering by a stratified sphere, i.e. a particle    !
   ! consisting of a spherical core surrounded by a spherical shell.  The basic code used  !
   ! was that described in the report: " subroutines for computing the parameters of the   !
   ! electromagnetic radiation scattered by a sphere ", J.V. Dave, I B M scientific        !
   ! center, Palo Alto , California. Report no. 320 - 3236 .. May 1968.                    !
   !    The modifications for stratified spheres are described in Toon and Ackerman, Appl. !
   ! Optics, in press, 1981                                                                !
   !    The parameters in the calling statement are defined as follows :                   !
   !    - ro is the outer (shell) radius;                                                  !
   !    - r  is the core radius;                                                           !
   !    - rfr, rfi  are the real and imaginary parts of the shell index of refraction in   !
   !          the form (rfr - i* rfi);                                                     !
   !    - re2, tmag2  are the index parts for the core;                                    !
   !          ( we assume space has unit index. )                                          !
   !    - thetd: angle in degrees between the directions of the incident and the           !
   !          scattered radiation.  thetd is< or= 90.0 if thetd should happen to be        !
   !          greater than 90.0, enter with supplementary value, see comments below on     !
   !          eltrmx;                                                                      !
   !                                                                                       !
   !    The definitions for the following symbols can be found in "Light scattering by     !
   ! small particles,H.C.van de Hulst, John Wiley' sons, inc., New York, 1957" .           !
   !    - qext: efficiency factor for extinction,  van de Hulst, p.14 ' 127.               !
   !    - qscat: efficiency factor for scattering, van de Hulst, p.14 ' 127.               !
   !    - ctbrqs: average(cosine theta) * qscat,   van de Hulst, p.128                     !
   !    - eltrmx(i,j,k): elements of the transformation matrix f, van de Hulst,            !
   !          p.34,45 ' 125. i=1: element m sub 2..i=2: element m sub 1..                  !
   !          i = 3: element s sub 21.. i = 4: element d sub 21..                          !
   !    - eltrmx(i,j,1) represents the ith element of the matrix for  the angle thetd(j).. !
   !          eltrmx(i,j,2) represents the ith element of the matrix for the angle         !
   !          180.0 - thetd(j) ..                                                          !
   !    - qbs is the back scatter cross section.                                           !
   !                                                                                       !
   !    - iacap is the dimension of acap in the original program the dimension of acap was !
   !          7000. for conserving space this should be not much higher than               !
   !          the value, n=1.1*(nreal**2 + nimag**2)**.5 * x + 1                           !
   !    - wvno: 2*pi / wavelength                                                          !
   !                                                                                       !
   !    Also, the subroutine computes the capital a function by making use of downward     !
   ! recurrence relationship.                                                              !
   !                                                                                       !
   !    - ta(1): real part of wfn(1).  ta(2): imaginary part of wfn(1).                    !
   !    - ta(3): real part of wfn(2).  ta(4): imaginary part of wfn(2).                    !
   !    - tb(1): real part of fna.     tb(2): imaginary part of fna.                       !
   !    - tc(1): real part of fnb.     tc(2): imaginary part of fnb.                       !
   !    - td(1): real part of fnap.    td(2): imaginary part of fnap.                      !
   !    - te(1): real part of fnbp.    te(2): imaginary part of fnbp.                      !
   !    - fnap, fnbp  are the preceding values of fna, fnb respectively.                   !
   !---------------------------------------------------------------------------------------!
   subroutine miess(ro,rfr,rfi,thetd,qext,qscat,ctbrqs,r,re2,tmag2,wvno)
      use rconstants, only: pio1808 ! ! intent(in)
      implicit none
      !----- Local parameters. ------------------------------------------------------------!
      integer                             , parameter        :: iacap = 20000
      real(kind=8)                        , parameter        :: epsilon_mie = 1.d-14
      !----- Arguments. -------------------------------------------------------------------!
      real                                , intent(in)       :: ro
      real                                , intent(in)       :: rfr
      real                                , intent(in)       :: rfi
      real                                , intent(inout)    :: thetd
      real                                , intent(out)      :: qext
      real                                , intent(out)      :: qscat
      real                                , intent(out)      :: ctbrqs
      real                                , intent(in)       :: r
      real                                , intent(in)       :: re2
      real                                , intent(in)       :: tmag2
      real                                , intent(in)       :: wvno
      !----- Local variables. -------------------------------------------------------------!
      complex(kind=16), dimension(iacap)                     :: acap
      complex(kind=16), dimension(2)                         :: wfn
      complex(kind=16), dimension(4)                         :: z
      complex(kind=16), dimension(3,iacap)                   :: w
      complex(kind=16), dimension(8)                         :: u
      complex(kind=16)                                       :: fnap
      complex(kind=16)                                       :: fnbp
      complex(kind=16)                                       :: fna
      complex(kind=16)                                       :: fnb
      complex(kind=16)                                       :: rf
      complex(kind=16)                                       :: rrf
      complex(kind=16)                                       :: rrfx
      complex(kind=16)                                       :: wm1
      complex(kind=16)                                       :: fn1
      complex(kind=16)                                       :: fn2
      complex(kind=16)                                       :: tc1
      complex(kind=16)                                       :: tc2
      complex(kind=16)                                       :: k1
      complex(kind=16)                                       :: k2
      complex(kind=16)                                       :: k3
      complex(kind=16)                                       :: rc
      complex(kind=16)                                       :: dh1
      complex(kind=16)                                       :: dh2
      complex(kind=16)                                       :: dh4
      complex(kind=16)                                       :: p24h24
      complex(kind=16)                                       :: p24h21
      complex(kind=16)                                       :: pstore
      complex(kind=16)                                       :: hstore
      complex(kind=16)                                       :: dummy
      complex(kind=16)                                       :: dumsq
      real(kind=8)    , dimension(5)                         :: t
      real(kind=8)    , dimension(4)                         :: ta
      real(kind=8)    , dimension(2)                         :: tb
      real(kind=8)    , dimension(2)                         :: tc
      real(kind=8)    , dimension(2)                         :: td
      real(kind=8)    , dimension(2)                         :: te
      real(kind=8)    , dimension(3)                         :: piz
      real(kind=8)    , dimension(3)                         :: tau
      real(kind=8)                                           :: cstht
      real(kind=8)                                           :: si2tht
      real(kind=8)    , dimension(4,2)                       :: eltrmx
      real(kind=8)                                           :: x
      real(kind=8)                                           :: x1
      real(kind=8)                                           :: x4
      real(kind=8)                                           :: y1
      real(kind=8)                                           :: y4
      real(kind=8)                                           :: rx
      real(kind=8)                                           :: sinx1
      real(kind=8)                                           :: sinx4
      real(kind=8)                                           :: cosx1
      real(kind=8)                                           :: cosx4
      real(kind=8)                                           :: ey1
      real(kind=8)                                           :: e2y1
      real(kind=8)                                           :: ey4
      real(kind=8)                                           :: ey1my4
      real(kind=8)                                           :: ey1py4
      real(kind=8)                                           :: aa
      real(kind=8)                                           :: bb
      real(kind=8)                                           :: cc
      real(kind=8)                                           :: dd
      real(kind=8)                                           :: denom
      real(kind=8)                                           :: realp
      real(kind=8)                                           :: amagp
      real(kind=8)                                           :: qbsr
      real(kind=8)                                           :: qbsi
      real(kind=8)                                           :: rmm
      logical                                                :: iflag
      integer                                                :: nmx1
      integer                                                :: nmx2
      integer                                                :: n
      integer                                                :: nn
      integer                                                :: m
      integer                                                :: k
      integer                                                :: i
      !------------------------------------------------------------------------------------!


      !----- Iflag is a logical test, used in several places. -----------------------------!
      iflag = r/ro < 1.d-6

      rf   =  cmplx( rfr,  -rfi )
      rc   =  cmplx( re2, -tmag2 )
      x    =  ro * wvno
      k1   =  rc * wvno
      k2   =  rf * wvno
      k3   =  cmplx( wvno, 0.d0 )
      z(1) =  k2 * ro
      z(2) =  k3 * ro
      z(3) =  k1 * r
      z(4) =  k2 * r
      x1   =  dble (z(1))
      x4   =  dble (z(4))
      y1   =  aimag(z(1))
      y4   =  aimag(z(4))
      rrf  =  1.d0 / rf
      rx   =  1.d0 / x
      rrfx =  rrf * rx
      t(1) =  ( x*x ) * ( rfr*rfr + rfi*rfi )
      t(1) =  sqrt( t(1) )
      nmx1 =  1.1d0 * t(1)

      if ( nmx1 > iacap-1 )  then
         write (unit=*,fmt='(a,1x,i)') ' NMX1  = ',nmx1
         write (unit=*,fmt='(a,1x,i)') ' IACAP = ',iacap
         call abort_run('The upper limit for acap is not enough.'                          &
                       ,'miess','rad_carma.f90')
      end if

      nmx2 = t(1)

      if ( nmx1 <= 150 ) then
         nmx1 = 150
         nmx2 = 135
      end if
      acap(nmx1+1 ) = dcmplx(0.d0,0.d0)

      if (.not. iflag) then
         do n = 1,3
           w(n,nmx1+1) = dcmplx(0.d0,0.d0)
         end do
      end if

      loop1: do n = 1,nmx1
         nn       = nmx1 - n + 1
         acap(nn) = (nn+1) * rrfx - 1.d0 / ( (nn+1) * rrfx + acap(nn+1))
         if (iflag) cycle loop1
         do m = 1,3
           w(m,nn) = (nn+1) / z(m+1) - 1.d0 / ((nn+1)/z(m+1) + w(m,nn+1))
         end do
      end do loop1

      loop2: do
         if ( thetd < 0.0 )  thetd = abs(thetd)

         if ( thetd <= 0.0 ) then
            cstht  = 1.0
            si2tht = 0.0
            exit loop2
         end if

         if ( thetd < 90. )  then
            t(1)   =  pio1808 * thetd
            cstht  =  cos( t(1) )
            si2tht =  1.d0 - cstht*cstht
            exit loop2
         end if

         if ( thetd <= 90. ) then
            cstht  =  0.0
            si2tht =  1.0
            exit loop2
         end if

         write (unit=*,fmt='(a,1x,es12.5)') ' THETD = ',thetd
         call abort_run('The value of the scattering angle is greater than 90.0 degrees!'  &
                       ,'miess','rad_carma.f90')
      end do loop2

         piz(1)  =  0.0
         piz(2)  =  1.0
         tau(1)  =  0.0
         tau(2)  =  cstht

      !----- Initialisation of homogeneous sphere. ----------------------------------------!
      t(1)   =  cos(x)
      t(2)   =  sin(x)
      wm1    =  cmplx( t(1),-t(2) )
      wfn(1) =  cmplx( t(2), t(1) )
      ta(1)  =  t(2)
      ta(2)  =  t(1)
      wfn(2) =  rx * wfn(1) - wm1
      ta(3)  =  dble(wfn(2))
      ta(4)  =  aimag(wfn(2))

      if (iflag)  then
         tc1   = acap(1) * rrf  +  rx
         tc2   = acap(1) * rf   +  rx
         fna   = (tc1 * ta(3) - ta(1)) / (tc1 * wfn(2) - wfn(1))
         fnb   = (tc2 * ta(3) - ta(1)) / (tc2 * wfn(2) - wfn(1))
         tb(1) = dble(fna)
         tb(2) = aimag(fna)
         tc(1) = dble(fnb)
         tc(2) = aimag(fnb)
      else
         n = 1
         !----- Initialisation procedure for stratified sphere begins here. ---------------!
         sinx1   =  sin(x1)
         sinx4   =  sin(x4)
         cosx1   =  cos(x1)
         cosx4   =  cos(x4)
         ey1     =  exp(y1)
         e2y1    =  ey1 * ey1
         ey4     =  exp(y4)
         ey1my4  =  exp(y1 - y4)
         ey1py4  =  ey1 * ey4
         ey1my4  =  exp(y1 - y4)
         aa      =  sinx4 * (ey1py4 + ey1my4)
         bb      =  cosx4 * (ey1py4 - ey1my4)
         cc      =  sinx1 * (e2y1 + 1.d0)
         dd      =  cosx1 * (e2y1 - 1.d0)
         denom   =  1.d0  +  e2y1 * ( 4.d0 * sinx1 * sinx1 - 2.d0 + e2y1 )
         realp   =  (aa * cc + bb * dd) / denom
         amagp   =  (bb * cc - aa * dd) / denom
         dummy   =  cmplx( realp, amagp )
         aa      =  sinx4 * sinx4 - 0.5
         bb      =  cosx4 * sinx4
         p24h24  =  5.d-1 + cmplx( aa,bb ) * ey4 * ey4
         aa      =  sinx1 * sinx4  -  cosx1 * cosx4
         bb      =  sinx1 * cosx4  +  cosx1 * sinx4
         cc      =  sinx1 * sinx4  +  cosx1 * cosx4
         dd      = -sinx1 * cosx4  +  cosx1 * sinx4
         p24h21  =  0.5 * cmplx( aa,bb ) * ey1 * ey4  + 0.5 * cmplx( cc,dd ) * ey1my4
         dh4     =  z(4) / ( 1.0 + dcmplx( 0.d0,1.d0 ) * z(4) )  -  1.d0 / z(4)
         dh1     =  z(1) / ( 1.0 + dcmplx( 0.d0,1.d0 ) * z(1) )  -  1.d0 / z(1)
         dh2     =  z(2) / ( 1.0 + dcmplx( 0.d0,1.d0 ) * z(2) )  -  1.d0 / z(2)
         pstore  =  ( dh4 + n / z(4) )  *  ( w(3,n) + n / z(4) )
         p24h24  =  p24h24 / pstore
         hstore  =  ( dh1 + n / z(1) )  *  ( w(3,n) + n / z(4) )
         p24h21  =  p24h21 / hstore
         pstore  =  ( acap(n) + n / z(1) )  /  ( w(3,n) + n / z(4) )
         dummy   =  dummy * pstore
         dumsq   =  dummy * dummy
         !---------------------------------------------------------------------------------!
         ! Note:  the definitions of u(i) in this program are not the same as the usubi    !
         !        defined in the article by Toon and Ackerman.  The corresponding terms    !
         !        are:                                                                     !
         !          usub1 = u(1)                       usub2 = u(5)                        !
         !          usub3 = u(7)                       usub4 = dumsq                       !
         !          usub5 = u(2)                       usub6 = u(3)                        !
         !          usub7 = u(6)                       usub8 = u(4)                        !
         !          ratio of spherical bessel ftn to spherical henkal ftn = u(8)           !
         !---------------------------------------------------------------------------------!
         u(1)  = k3 * acap(n)  -  k2 * w(1,n)
         u(2)  = k3 * acap(n)  -  k2 * dh2
         u(3)  = k2 * acap(n)  -  k3 * w(1,n)
         u(4)  = k2 * acap(n)  -  k3 * dh2
         u(5)  = k1 *  w(3,n)  -  k2 * w(2,n)
         u(6)  = k2 *  w(3,n)  -  k1 * w(2,n)
         u(7)  = dcmplx( 0.d0,-1.d0 )  *  ( dummy * p24h21 - p24h24 )
         u(8)  = ta(3) / wfn(2)
         fna   = u(8) * (u(1)*u(5)*u(7) + k1*u(1) - dumsq*k3*u(5))                         &
                      / (u(2)*u(5)*u(7) + k1*u(2) - dumsq*k3*u(5))
         fnb   = u(8) * (u(3)*u(6)*u(7) + k2*u(3) - dumsq*k2*u(6))                         &
                      / (u(4)*u(6)*u(7) + k2*u(4) - dumsq*k2*u(6))
         tb(1) = dble(fna)
         tb(2) = aimag(fna)
         tc(1) = dble(fnb)
         tc(2) = aimag(fnb)
      end if

      fnap  = fna
      fnbp  = fnb
      td(1) = dble(fnap)
      td(2) = aimag(fnap)
      te(1) = dble(fnbp)
      te(2) = aimag(fnbp)
      t(1)  = 1.50
      !------------------------------------------------------------------------------------!
      !    From here on eltrmx(i,j,k) has the following meaning:                           !
      !    - eltrmx(1,j,k): real part of the first complex amplitude.                      !
      !    - eltrmx(2,j,k): imaginary part of the first complex amplitude.                 !
      !    - eltrmx(3,j,k): real part of the second complex amplitude.                     !
      !    - eltrmx(4,j,k): imaginary part of the second complex amplitude.                !
      !    - k = 1 : for thetd(j) and k = 2 : for 180.0 - thetd(j)                         !
      !    - definition of the complex amplitude: van de Hulst,p.125.                      !
      !------------------------------------------------------------------------------------!
      tb(1) = t(1) * tb(1)
      tb(2) = t(1) * tb(2)
      tc(1) = t(1) * tc(1)
      tc(2) = t(1) * tc(2)
      fna   = dcmplx(tb(1),tb(2))
      fnb   = dcmplx(tc(1),tc(2))

      eltrmx(1,1) = tb(1) * piz(2) + tc(1) * tau(2)
      eltrmx(2,1) = tb(2) * piz(2) + tc(2) * tau(2)
      eltrmx(3,1) = tc(1) * piz(2) + tb(1) * tau(2)
      eltrmx(4,1) = tc(2) * piz(2) + tb(2) * tau(2)
      eltrmx(1,2) = tb(1) * piz(2) - tc(1) * tau(2)
      eltrmx(2,2) = tb(2) * piz(2) - tc(2) * tau(2)
      eltrmx(3,2) = tc(1) * piz(2) - tb(1) * tau(2)
      eltrmx(4,2) = tc(2) * piz(2) - tb(2) * tau(2)

      qext   = 2.0 * sngl( tb(1) + tc(1))
      qscat  = sngl( tb(1)* tb(1) + tb(2)* tb(2) + tc(1)*tc(1) + tc(2)*tc(2) ) / 0.75
      ctbrqs = 0.0
      qbsr   = -2.d0*(tc(1) - tb(1))
      qbsi   = -2.d0*(tc(2) - tb(2))
      rmm    = -1.d0
      n = 2

      bigloop: do
         t(1)   = 2*n - 1
         t(2)   =   n - 1
         t(3)   = 2*n + 1
         piz(3) = (t(1)*piz(2)*cstht - n*piz(1)) / t(2)
         tau(3) = cstht * (piz(3)-piz(1)) - t(1)*si2tht*piz(2) + tau(1)
         !----- Here we set up the homogeneous sphere. ------------------------------------!
         wm1    =  wfn(1)
         wfn(1) =  wfn(2)
         ta(1)  =  dble(wfn(1))
         ta(2)  =  aimag(wfn(1))
         ta(4)  =  aimag(wfn(2))
         wfn(2) =  t(1) * rx * wfn(1)  -  wm1
         ta(3)  =  dble(wfn(2))

         if (.not. iflag ) then
            !
            !----- Here we set up the stratified sphere. ----------------------------------!
            dh2     =  - n / z(2) + 1.0 / ( n / z(2) - dh2 )
            dh4     =  - n / z(4) + 1.0 / ( n / z(4) - dh4 )
            dh1     =  - n / z(1) + 1.0 / ( n / z(1) - dh1 )
            pstore  =  ( dh4 + n / z(4) )  *  ( w(3,n) + n / z(4) )
            p24h24  =  p24h24 / pstore
            hstore  =  ( dh1 + n / z(1) )  *  ( w(3,n) + n / z(4) )
            p24h21  =  p24h21 / hstore
            pstore  =  ( acap(n) + n / z(1) )  /  ( w(3,n) + n / z(4) )
            dummy   =  dummy * pstore
            dumsq   =  dummy * dummy
            u(1) =  k3 * acap(n)  -  k2 * w(1,n)
            u(2) =  k3 * acap(n)  -  k2 * dh2
            u(3) =  k2 * acap(n)  -  k3 * w(1,n)
            u(4) =  k2 * acap(n)  -  k3 * dh2
            u(5) =  k1 *  w(3,n)  -  k2 * w(2,n)
            u(6) =  k2 *  w(3,n)  -  k1 * w(2,n)
            u(7) =  ( 0.0,-1.0 )  *  ( dummy * p24h21 - p24h24 )
            u(8) =  ta(3) / wfn(2)
            fna  =  u(8) * (u(1)*u(5)*u(7) + k1*u(1) - dumsq*k3*u(5))                      &
                         / (u(2)*u(5)*u(7) + k1*u(2) - dumsq*k3*u(5))
            fnb  =  u(8) * (u(3)*u(6)*u(7) + k2*u(3) - dumsq*k2*u(6))                      &
                         / (u(4)*u(6)*u(7) + k2*u(4) - dumsq*k2*u(6))
            tb(1) = dble(fna)
            tb(2) = aimag(fna)
            tc(1) = dble(fnb)
            tc(2) = aimag(fnb)
         end if

         tc1 = acap(n) * rrf + n * rx
         tc2 = acap(n) * rf  + n * rx
         fn1 = (tc1 * ta(3) - ta(1)) / (tc1 * wfn(2) - wfn(1))
         fn2 = (tc2 * ta(3) - ta(1)) / (tc2 * wfn(2) - wfn(1))
         m   = wvno * r
         if ( n >= m ) then
            if (.not. iflag) then
               iflag = abs((fn1-fna)/fn1) < epsilon_mie .and.                              &
                       abs(  ( fn2-fnb ) / fn2  ) < epsilon_mie  
            end if
            if (iflag) then
               fna   =  fn1
               fnb   =  fn2
               tb(1) = dble(fna)
               tb(2) = aimag(fna)
               tc(1) = dble(fnb)
               tc(2) = aimag(fnb)
            end if
         end if
         t(5)    = n
         t(4)    = t(1) / (t(5) * t(2))
         t(2)    = (t(2) * (t(5) + 1.0)) / t(5)
         ctbrqs  = ctbrqs + t(2) * (td(1)*tb(1) + td(2)*tb(2) + te(1)*tc(1) + te(2)*tc(2)) &
                          + t(4) * (td(1)*te(1) + td(2)*te(2))
         qext    = qext + t(3) * (tb(1) + tc(1))
         t(4)    =  tb(1)*tb(1) + tb(2)*tb(2) + tc(1)*tc(1) + tc(2)*tc(2)
         qscat   =  qscat  +  t(3) * t(4)
         rmm     =  -rmm
         qbsr    =  qbsr + t(3)*rmm*(tc(1) - tb(1))
         qbsi    =  qbsi + t(3)*rmm*(tc(2) - tb(2))
         t(2)    =  n * (n+1)
         t(1)    =  t(3) / t(2)
         k       = (n/2)*2

         eltrmx(1,1) = eltrmx(1,1)+t(1)*(tb(1)*piz(3)+tc(1)*tau(3))
         eltrmx(2,1) = eltrmx(2,1)+t(1)*(tb(2)*piz(3)+tc(2)*tau(3))
         eltrmx(3,1) = eltrmx(3,1)+t(1)*(tc(1)*piz(3)+tb(1)*tau(3))
         eltrmx(4,1) = eltrmx(4,1)+t(1)*(tc(2)*piz(3)+tb(2)*tau(3))
         if ( k == n )  then
           eltrmx(1,2) = eltrmx(1,2)+t(1)*(-tb(1)*piz(3)+tc(1)*tau(3)) 
           eltrmx(2,2) = eltrmx(2,2)+t(1)*(-tb(2)*piz(3)+tc(2)*tau(3)) 
           eltrmx(3,2) = eltrmx(3,2)+t(1)*(-tc(1)*piz(3)+tb(1)*tau(3)) 
           eltrmx(4,2) = eltrmx(4,2)+t(1)*(-tc(2)*piz(3)+tb(2)*tau(3)) 
         else
           eltrmx(1,2) = eltrmx(1,2)+t(1)*(tb(1)*piz(3)-tc(1)*tau(3))
           eltrmx(2,2) = eltrmx(2,2)+t(1)*(tb(2)*piz(3)-tc(2)*tau(3))
           eltrmx(3,2) = eltrmx(3,2)+t(1)*(tc(1)*piz(3)-tb(1)*tau(3))
           eltrmx(4,2) = eltrmx(4,2)+t(1)*(tc(2)*piz(3)-tb(2)*tau(3))
         end if

         if ( t(4) < epsilon_mie ) exit bigloop

         n = n + 1
         piz(1) = piz(2)
         piz(2) = piz(3)
         tau(1) = tau(2)
         tau(2) = tau(3)
         fnap   = fna
         fnbp   = fnb
         td(1)  = dble(fnap)
         td(2)  = aimag(fnap)
         te(1)  = dble(fnbp)
         te(2)  = aimag(fnbp)
         if ( n > nmx2 ) then
            write (unit=*,fmt='(a,1x,i6)') 'N    = ',n
            write (unit=*,fmt='(a,1x,i6)') 'NMX2 = ',nmx2
            call abort_run('The upper limit for acap is not large enough.'                 &
                          ,'miess','rad_carma.f90')
         end if
      end do bigloop

      do k = 1,2
         do i= 1,4
            t(i) = eltrmx(i,k)
         end do
         eltrmx(2,k) = t(1)*t(1) + t(2)*t(2)
         eltrmx(1,k) = t(3)*t(3) + t(4)*t(4)
         eltrmx(3,k) = t(1)*t(3) + t(2)*t(4)
         eltrmx(4,k) = t(2)*t(3) - t(4)*t(1)
      end do
      t(1)   =   2.0 * rx * rx
      qext   =  qext * t(1)
      qscat  = qscat * t(1)
      ctbrqs =   2.0 * ctbrqs * t(1)

      return
   end subroutine miess
   !=======================================================================================!
   !=======================================================================================!







   !=======================================================================================!
   !=======================================================================================!
   subroutine radtran_to_rams(m1,fthrl,rlong,rlongup_top,fthrs,rshort,rshort_top           &
                             ,rshortup_top,aotr,mynum)
      use mem_grid   , only : nzpmax  ! ! intent(in)
      use mem_globrad, only : nwave   & ! intent(in)
                            , nsolp   & ! intent(in)
                            , ntotal  & ! intent(in)
                            , nlayer  & ! intent(in)
                            , ngauss  & ! intent(in)
                            , gweight & ! intent(in)
                            , nprob   ! ! intent(in)
      use mem_radiate, only : ilwrtyp ! ! intent(in)
      implicit none
  
      !----- Arguments. -------------------------------------------------------------------!
      integer                   , intent(in)    :: m1,mynum
      real                      , intent(out)   :: rlong
      real                      , intent(out)   :: rlongup_top
      real                      , intent(out)   :: rshort
      real                      , intent(out)   :: rshort_top
      real                      , intent(out)   :: rshortup_top
      real   , dimension(nwave) , intent(out)   :: aotr
      real   , dimension(nzpmax), intent(inout) :: fthrl
      real   , dimension(nzpmax), intent(inout) :: fthrs
      !----- Local variables. -------------------------------------------------------------!
      real   , dimension(nzpmax)                :: duml
      real   , dimension(nzpmax)                :: dums
      real   , dimension(nwave)                 :: dumaot
      real   , dimension(ntotal,nzpmax)         :: dum2aot
      integer                                   :: k
      integer                                   :: j1
      integer                                   :: i1
      integer                                   :: k1
      integer                                   :: kr
      integer                                   :: ib
      integer                                   :: ig
      integer                                   :: nzz
      integer                                   :: l
      !------------------------------------------------------------------------------------!

      nzz  = m1 - 1
      aotr = 0.0

      !----- Reverse the vertical and transfer values from CARMA grid to BRAMS grid. ------!
      do k=2,m1
         kr = nzz+2- k 
         fthrl(k) = heati_aerad(kr)
         fthrs(k) = heats_aerad(kr)
      end do


      rshort       = solnet      ! Total short wave absorbed by the surface.
      rshort_top   = soldowntoa  ! Total incoming short wave at the top layer.
      rshortup_top = soluptoa    ! Total outgoing short wave at the top layer.
      rlong        = xirdown ! Surface long wave radiation.
      rlongup_top  = xirup   ! Emerging long wave radiation at the top layer.
      if (rlong /= rlong .and. ilwrtyp == 4) then 
         write(unit=*,fmt='(a)') '-------------------------------------------------------'
         call abort_run('Weird RLONG... ','radtrans_to_rams','rad_carma.f90')
      end if
      dumaot=0.0

      do l=1,ntotal
         do k=1,m1                                  
            dum2aot(nprob(l),k)=tauaer(nprob(l),k)
         end do                                     
      end do

      do l=1,nwave
         do k=1,m1                                  
            dumaot(l) = dumaot(l) + dum2aot(l,k)   
         end do                                     
      end do
     
      do l=1,nwave
         aotr(l)= dumaot(l)
      end do

      return
   end subroutine radtran_to_rams
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine converts the 3-D, wavenumber-like arrays into 2-D (X*Y,WN).        !
   !---------------------------------------------------------------------------------------!
   subroutine ci_3d_1d(m2,m3,m4,a3d,a1d,i,j)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      integer                     , intent(in)  :: m2
      integer                     , intent(in)  :: m3
      integer                     , intent(in)  :: m4
      integer                     , intent(in)  :: i
      integer                     , intent(in)  :: j
      real   , dimension(m2,m3,m4), intent(in)  :: a3d
      real   , dimension(m4)      , intent(out) :: a1d
      !----- Local variables. -------------------------------------------------------------!
      integer                                   :: w
      !------------------------------------------------------------------------------------!

      do w=1,m4
         a1d(w)=a3d(i,j,w)
      end do

      return
   end subroutine ci_3d_1d
   !=======================================================================================!
   !=======================================================================================!







   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine converts back the 2-D, wavenumber-like arrays into 3-D (X,Y,WN).   !
   !---------------------------------------------------------------------------------------!
   subroutine ci_1d_3d(m2,m3,m4,a1d,a3d,i,j)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      integer                     , intent(in)    :: m2
      integer                     , intent(in)    :: m3
      integer                     , intent(in)    :: m4
      integer                     , intent(in)    :: i
      integer                     , intent(in)    :: j
      real   , dimension(m2,m3,m4), intent(inout) :: a3d
      real   , dimension(m4)      , intent(in)    :: a1d
      !----- Local variables. -------------------------------------------------------------!
      integer                                     :: w
      !------------------------------------------------------------------------------------!

      do w=1,m4
         a3d(i,j,w)=a1d(w)
      end do

      return
   end subroutine ci_1d_3d
   !=======================================================================================!
   !=======================================================================================!
end module rad_carma
!==========================================================================================!
!==========================================================================================!
