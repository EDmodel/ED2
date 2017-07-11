!===================================== Change Log =========================================!
! 5.0.0                                                                                    !
!                                                                                          !
!==========================================================================================!
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved                       !
!  Regional Atmospheric Modeling System - RAMS                                             !
!==========================================================================================!
!    This subroutine will drive the initialization of model grid, topography, and fields   !
! for horizontally homogeneous fields.                                                     !
!------------------------------------------------------------------------------------------!
subroutine inithh()
   use mem_basic  , only : basic_g  & ! intent(inout)
                         , co2_on   & ! intent(in)
                         , co2con   ! ! intent(in)
   use mem_grid   , only : nxtnest  & ! intent(in)
                         , grid_g   & ! intent(in)
                         , ngrids   & ! intent(in)
                         , nzp      & ! intent(in)
                         , nxp      & ! intent(in)
                         , nyp      & ! intent(in)
                         , if_adap  ! ! intent(in)
   use mem_scratch, only : scratch  ! ! intent(inout)
   implicit none
   !----- Local variables. ----------------------------------------------------------------!
   integer :: ifm
   integer :: icm
   !---------------------------------------------------------------------------------------!

   !----- Arrange the input sounding. -----------------------------------------------------!
   call arrsnd(co2_on,co2con(1))


   !---------------------------------------------------------------------------------------!
   !     For GRID 1, and the second hemispheric grid if applicable, compute the 1-D        !
   ! reference state variables, the 3-D reference state variables, the 3-D model fields,   !
   ! the surface layer parameters, and the initial soil model fields.                      !
   !---------------------------------------------------------------------------------------!
   do ifm = 1,ngrids
      icm = nxtnest(ifm)
      if (icm == 0) then

         call newgrid(ifm)
         call refs1d(co2_on,co2con(1))

         call refs3d( nzp,nxp,nyp           , basic_g(ifm)%pi0      , basic_g(ifm)%dn0     &
                    , basic_g(ifm)%dn0u     , basic_g(ifm)%dn0v     , basic_g(ifm)%th0     &
                    , grid_g(ifm)%topt      , grid_g(ifm)%rtgt      )


         select case (if_adap)
         case (0)
            call flds3d( nzp,nxp,nyp          , basic_g(ifm)%uc      , basic_g(ifm)%vc     &
                       , basic_g(ifm)%pi0     , basic_g(ifm)%theta   , basic_g(ifm)%thp    &
                       , basic_g(ifm)%rtp     , basic_g(ifm)%pc      , basic_g(ifm)%rv     &
                       , scratch%vt3do        , grid_g(ifm)%topt     , grid_g(ifm)%topu    &
                       , grid_g(ifm)%topv     , grid_g(ifm)%rtgt     , grid_g(ifm)%rtgu    &
                       , grid_g(ifm)%rtgv     )

         case (1)
            call flds3d_adap( nzp,nxp,nyp        , grid_g(ifm)%flpu   , grid_g(ifm)%flpv   &
                            , grid_g(ifm)%flpw   , basic_g(ifm)%uc    , basic_g(ifm)%vc    &
                            , basic_g(ifm)%pi0   , basic_g(ifm)%theta , basic_g(ifm)%thp   &
                            , basic_g(ifm)%rtp   , basic_g(ifm)%pc    , basic_g(ifm)%rv    &
                            , scratch%vt3do      )
         end select

         !----- Copy the results for CO2 only if there is a CO2 array allocated. ----------!
         if (co2_on) call atob(nzp*nxp*nyp,scratch%vt3do,basic_g(ifm)%co2p)

      end if
   end do
   return
end subroutine inithh
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine will arrange the input sounding so everything will be in standard    !
! units.                                                                                   !
!------------------------------------------------------------------------------------------!
subroutine arrsnd(co2_on,co2con)
   use mem_grid
   use mem_scratch
   use ref_sounding
   use rconstants
   use therm_lib    , only : rslf       & ! function
                           , ptrh2rvapl & ! function
                           , virtt      ! ! function
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   logical           , intent(in)   :: co2_on
   real              , intent(in)   :: co2con
   !----- Local variables. ----------------------------------------------------------------!
   integer                          :: nnns
   integer                          :: k
   integer                          :: here
   integer                          :: beneath
   integer                          :: above
   integer                          :: ierr
   real                             :: toffset
   real                             :: dir
   real                             :: spd
   real                             :: zold1
   real                             :: zold2
   real                             :: tavg
   real                             :: wt
   real, dimension(1)               :: rtss
   real              , parameter    :: undefflg = 9999. 
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      First of all, check whether ps(1) is zero.  This is the flag to tell the model   !
   ! to skip the namelist variables and read the SOUND_IN file.                            !
   !---------------------------------------------------------------------------------------!
   if (ps(1) == 0.) then

      do nsndg=1,maxsndg
         ps  (nsndg) = 0.
         ts  (nsndg) = 0.
         rts (nsndg) = 0.
         us  (nsndg) = 0.
         vs  (nsndg) = 0.
         co2s(nsndg) = 0.
      end do

      !----- Open the SOUND_IN file. ------------------------------------------------------!
      open (unit=19,file='SOUND_IN',status='old',form='formatted')

      !------------------------------------------------------------------------------------!
      !     Now we will read the fields, acknowledging that CO2S may not be there if the   !
      ! user is not running with CO2.                                                      !
      !------------------------------------------------------------------------------------!
      if (co2_on) then
         readloop1: do nsndg=1,maxsndg
            read(unit=19,fmt=*,iostat=ierr) ps(nsndg),ts(nsndg),rts(nsndg)                 &
                                           ,us(nsndg),vs(nsndg),co2s(nsndg)
            if(ps(nsndg) <= 0. .or. ierr /= 0 ) exit readloop1
         end do readloop1
      else
         readloop2: do nsndg=1,maxsndg
            read(unit=19,fmt=*,iostat=ierr) ps(nsndg),ts(nsndg),rts(nsndg)                 &
                                           ,us(nsndg),vs(nsndg)
            if(ps(nsndg) <= 0. .or. ierr /= 0 ) exit readloop2
            co2s(nsndg) = co2con
         end do readloop2
      end if
      !----- Close the SOUND_IN file. -----------------------------------------------------!
      close(unit=19,status='keep')
   end if

   !----- Deciding the temperature offset based on the units the user provided. -----------!
   select case (itsflg)
   case (0) !----- Degrees Celsius. -------------------------------------------------------!
      toffset = t00
   case default !----- Kelvin. ------------------------------------------------------------!
      toffset = 0.
   end select

   mainloop: do nsndg=1,maxsndg
      nnns=nsndg
      !----- If we've reached the last level with actual data, exit the loop... -----------!
      if(ps(nsndg) == 0.) exit mainloop

      !----- If wind is defined and the user provided speed and direction, decompose it. --!
      if (us(nsndg) /= undefflg .and. iusflg /= 0 )then
         dir       =  us(nsndg)
         spd       =  vs(nsndg)
         us(nsndg) = -spd * sin(pio180*dir)
         vs(nsndg) = -spd * cos(pio180*dir)
      end if

      !----- Converting pressure to Pascals. ----------------------------------------------!
      select case (ipsflg) 
      case (0) !----- Pressure given in millibars. ----------------------------------------!
         ps(nsndg) = ps(nsndg) * 1.e2
      case (1) !----- Pressure is height in meters with PS(1)=surface pressure. -----------!

         !---------------------------------------------------------------------------------!
         !     If sounding moisture is expressed as a mixing ratio (IRTSFLG=2), take       !
         ! advantage of knowing virtual temperature effect when integrating hydrostatic-   !
         ! ally to get sounding pressure.                                                  !
         !---------------------------------------------------------------------------------!
         select case (irtsflg)
         case (2)
            vctr4(nsndg) = virtt(ts(nsndg)+toffset, rts(nsndg)*1.e-3)
         case default
            vctr4(nsndg) = ts(nsndg)+toffset
         end select

         if(nsndg == 1)then
            ps(nsndg) = ps(nsndg)*100.
            zold2 = 0.
         else
            zold1 = zold2
            zold2 = ps(nsndg)

            select case (itsflg)
            case (0,1) !----- Temperature. ------------------------------------------------!
               tavg      = 0.5 * (vctr4(nsndg) + vctr4(nsndg-1))
               ps(nsndg) = ps(nsndg-1) * exp(-grav * (zold2-zold1) / (rdry*tavg))
            case (2) !----- Potential temperature. ----------------------------------------!
               tavg      = (vctr4(nsndg) + vctr4(nsndg-1)*p00k/ps(nsndg-1)**rocp) * .5
               ps(nsndg) = (ps(nsndg-1)**rocp                                              &
                         - grav * (zold2-zold1) * p00k/(cpdry*tavg)) ** cpor
            end select
         end if
      case default
         write(unit=*,fmt='(a,1x,i2,a)') 'Invalid pressure type (IPSFLG): ',ipsflg,'...'
         call abort_run('Incorrect pressure unit in MODEL_SOUND, check your namelist!'     &
                       ,'arrsnd','rhhi.f90')
      end select

      !----- Deciding what to do with the temperature array. ------------------------------!
      select case (itsflg)
      case (0) !----- Temperature in degrees Celsius. -------------------------------------!
         ts(nsndg) = ts(nsndg) + t00
      case (1) !----- Temperature in Kelvin, no need to do anything. ----------------------!
         continue
      case (2) !----- Temperature is potential temperature in Kelvin. ---------------------!
         ts(nsndg) = (ps(nsndg)*p00i)**rocp * ts(nsndg)
      case default
         write(unit=*,fmt='(a,1x,i2,a)') 'Invalid temperature type (ITSFLG): ',itsflg,'...'
         call abort_run('Incorrect temperature unit in MODEL_SOUND, check your namelist!'  &
                       ,'arrsnd','rhhi.f90')
      end select

      !------------------------------------------------------------------------------------!
      !     Deciding what to do with the humidity.  In the end, it must become water       !
      ! vapour mixing ratio in kg/kg.                                                      !
      !------------------------------------------------------------------------------------!
      select case (irtsflg)
      case (0) !----- Humidity given as dew point in degrees Celsius. ---------------------!
         rts(nsndg) = rslf(ps(nsndg),rts(nsndg)+t00)

      case (1) !----- Humidity given as dew point in Kelvin. ------------------------------!
         rts(nsndg) = rslf(ps(nsndg),rts(nsndg))

      case (2) !----- Humidity given as mixing ratio in g/kg. -----------------------------!
         rts(nsndg) = rts(nsndg) * 1.e-3
      case (3) !----- Humidity given as relative humidity in percent. ---------------------!
         rts(nsndg) = ptrh2rvapl(rts(nsndg)*.01,ps(nsndg),ts(nsndg),.false.)
      case (4) !----- Humidity given as dew point depression in Kelvin. -------------------!
         rts(nsndg) = rslf(ps(nsndg),ts(nsndg)-rts(nsndg))
      case default
         write(unit=*,fmt='(a,1x,i2,a)') 'Invalid humidity type (IRTSFLG): ',irtsflg,'...'
         call abort_run('Incorrect humidity unit in MODEL_SOUND, check your namelist!'     &
                       ,'arrsnd','rhhi.f90')
      end select
   end do mainloop

   !---------------------------------------------------------------------------------------!
   !    We now look for levels with missing wind information, and perform an interpolation !
   ! between the last level with data and the next one.                                    !
   !---------------------------------------------------------------------------------------!
   nsndg=nnns-1
   do here=1,nsndg
      if (us(here) == undefflg) then
         !----- Find the level beneath the current that has data. -------------------------!
         beneath_loop: do beneath = here,1,-1
            if (us(beneath) /= undefflg) exit beneath_loop
         end do beneath_loop

         above_loop: do above = here,nsndg,1
            if(us(above) /= undefflg) exit above_loop
         end do above_loop

         wt       = (ps(here) - ps(beneath)) / (ps(above) - ps(beneath))
         us(here) = us(beneath) + wt * (us(above) - us(beneath))
         vs(here) = vs(beneath) + wt * (vs(above) - vs(beneath))

      end if
   end do

   !----- Compute height levels of input sounding. ----------------------------------------!
   do k=2,nsndg
      hs(k) = hs(k-1) - rdry * .5 * (virtt(ts(k),rts(k)) + virtt(ts(k-1),rts(k-1)))        &
                                  * (log(ps(k)) - log(ps(k-1))) / grav
   end do

   !----- Check whether the provided sounding goes high enough. ---------------------------!
   if (hs(nsndg) < zt(nzp)) then
      write (unit=*,fmt='(a,1x,es12.5)') 'Sounding top [  m] : ',hs(nsndg)
      write (unit=*,fmt='(a,1x,es12.5)') 'Sounding top [hPa] : ',0.01*ps(nsndg)
      write (unit=*,fmt='(a,1x,es12.5)') 'Model top    [  m] : ',zt(nzp)
      call abort_run('Input sounding needs to go higher!!!','arrsnd','rhhi.f90')
   end if

   !----- Find the theta profile. ---------------------------------------------------------!
   do k=1,nsndg
      thds(k) = ts(k) * (p00/ps(k))**rocp
   end do
   
   !---------------------------------------------------------------------------------------!
   !   CO2 has no option for units, so it must be given in ppm [µmol_CO2/mol_air].  BRAMS  !
   ! also uses CO2 in ppm, so we don't need to do any unit conversion.                     !
   !---------------------------------------------------------------------------------------!

   return
end subroutine arrsnd
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This routine computes the reference state sounding on the model sigma-z levels from   !
! input sounding defined on pressure levels.                                               !
!------------------------------------------------------------------------------------------!
subroutine refs1d(co2_on,co2con)
   use mem_grid
   use mem_scratch
   use ref_sounding
   use rconstants
   use therm_lib   , only : virtt        & ! function
                          , vapour_on    & ! intent(in)
                          , exner2press  & ! function
                          , extheta2temp ! ! function
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   logical           , intent(in)   :: co2_on
   real              , intent(in)   :: co2con
   !----- Local variables. ----------------------------------------------------------------!
   integer                          :: k
   !---------------------------------------------------------------------------------------!

   if (ztn(nnzp(ngrid),ngrid) > hs(nsndg)) then
      write (unit=*,fmt='(a,1x,es12.5)') ' - Sounding top [m]  : ',hs(nsndg)
      write (unit=*,fmt='(a,1x,es12.5)') ' - Model top    [m]  : ',ztn(nnzp(ngrid),ngrid)
      call abort_run('Input sounding is not high enough!!','refs1d','rhhi.f90')
   end if

   call htint(nsndg,thds,hs,nnzp(ngrid),           vctr1,ztn(1,ngrid))
   call htint(nsndg,  us,hs,nnzp(ngrid),  u01dn(1,ngrid),ztn(1,ngrid))
   call htint(nsndg,  vs,hs,nnzp(ngrid),  v01dn(1,ngrid),ztn(1,ngrid))

   !----- Fill water vapour reference state only if this run solves water. ----------------!
   if (vapour_on) then
      call htint(nsndg, rts,hs,nnzp(ngrid), rt01dn(1,ngrid),ztn(1,ngrid))
   else
      do k = 1,nnzp(ngrid)
         rt01dn(k,ngrid) = 0.
      enddo
   endif

   !----- Fill carbon dioxide reference state only if this run solves CO2. ----------------!
   if (co2_on) then
      call htint(nsndg,co2s,hs,nnzp(ngrid),co201dn(1,ngrid),ztn(1,ngrid))
   else
      do k = 1,nnzp(ngrid)
         co201dn(k,ngrid) = co2con
      end do
   end if

   !----- Convert the reference state to virtual potential temperature. -------------------!
   do k = 1,nnzp(ngrid)
      th01dn(k,ngrid) = virtt(vctr1(k),rt01dn(k,ngrid))
   end do

   !----- Boundary conditions, copying second level into the first. -----------------------!
   u01dn(1,ngrid)   = u01dn(2,ngrid)
   v01dn(1,ngrid)   = v01dn(2,ngrid)
   rt01dn(1,ngrid)  = rt01dn(2,ngrid)
   th01dn(1,ngrid)  = th01dn(2,ngrid)
   co201dn(1,ngrid) = co201dn(2,ngrid)

   !----- Finding the reference Exner function, using pressure and hydrostatic assumption. !
   pi01dn(1,ngrid) = cpdry * (ps(1) * p00i) ** rocp                                        &
                   + grav * (hs(1) - ztn(1,ngrid))                                         &
                   / (.5 * (th01dn(1,ngrid) + virtt(thds(1),rts(1)) ) )
   do k = 2,nnzp(ngrid)
      pi01dn(k,ngrid) = pi01dn(k-1,ngrid)                                                  &
                      - grav/(dzmn(k-1,ngrid) * .5 * (th01dn(k,ngrid) + th01dn(k-1,ngrid)))
   end do

   !----- Finding the reference density. --------------------------------------------------!
   do k = 1,nnzp(ngrid)
      vctr4(k) = exner2press(pi01dn(k,ngrid))
      dn01dn(k,ngrid) = vctr4(k) / (rdry * extheta2temp(pi01dn(k,ngrid),th01dn(k,ngrid)))
   end do

   return
end subroutine refs1d
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This routine initializes the 3-D velocity and thermodynamic fields from the 1-D       !
! reference state sounding when the terrain-following coordinate is used.                  !
!------------------------------------------------------------------------------------------!
subroutine flds3d(n1,n2,n3,uc,vc,pi0,theta,thp,rtp,pc,rv,co2p,topt,topu,topv,rtgt,rtgu,rtgv)
   use mem_grid
   use mem_scratch
   use ref_sounding
   use rconstants
   use mem_basic   , only : co2_on       ! ! intent(in)
   use therm_lib   , only : rslf         & ! function
                          , theta_iceliq & ! function
                          , tv2temp      & ! function
                          , vapour_on    & ! intent(in)
                          , cloud_on     & ! intent(in)
                          , exner2press  & ! function
                          , extheta2temp ! ! function
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer , intent(in)                     :: n1
   integer , intent(in)                     :: n2
   integer , intent(in)                     :: n3
   real, dimension(n1,n2,n3), intent(in)    :: pi0
   real, dimension(n1,n2,n3), intent(out)   :: thp
   real, dimension(n1,n2,n3), intent(out)   :: theta
   real, dimension(n1,n2,n3), intent(out)   :: rtp
   real, dimension(n1,n2,n3), intent(out)   :: pc
   real, dimension(n1,n2,n3), intent(out)   :: rv
   real, dimension(n1,n2,n3), intent(out)   :: uc
   real, dimension(n1,n2,n3), intent(out)   :: vc
   real, dimension(n1,n2,n3), intent(out)   :: co2p
   real, dimension(   n2,n3), intent(in)    :: topt 
   real, dimension(   n2,n3), intent(in)    :: topu 
   real, dimension(   n2,n3), intent(in)    :: topv 
   real, dimension(   n2,n3), intent(in)    :: rtgt
   real, dimension(   n2,n3), intent(in)    :: rtgu
   real, dimension(   n2,n3), intent(in)    :: rtgv
   !----- Local variables. ----------------------------------------------------------------!
   integer                                  :: i
   integer                                  :: j
   integer                                  :: k
   real                                     :: qlatu
   real                                     :: qlonu
   real                                     :: qlatv
   real                                     :: qlonv
   real                                     :: dummy
   real, dimension(nzpmax)                  :: p0
   real, dimension(nzpmax)                  :: temp
   real, dimension(nzpmax)                  :: rvls
   real, dimension(nzpmax)                  :: rc
   !---------------------------------------------------------------------------------------!


   jloop: do j=1,nyp
      iloop: do i=1,nxp
         !----- Find the actual height using the terrain-following coordinate. ------------!
         do k=1,nzp
            vctr11(k) = zt(k) * rtgt(i,j) + topt(i,j)
            vctr12(k) = zt(k) * rtgu(i,j) + topu(i,j)
            vctr13(k) = zt(k) * rtgv(i,j) + topv(i,j)
         enddo

         !----- Interpolate the fields using the reference state. -------------------------!
         call htint(nzp, u01dn(:,ngrid),zt,nzp,vctr5,vctr12)
         call htint(nzp, v01dn(:,ngrid),zt,nzp,vctr6,vctr13)
         call htint(nzp,th01dn(:,ngrid),zt,nzp,vctr3,vctr11)
         if(vapour_on) call htint(nzp, rt01dn(:,ngrid),zt,nzp, rtp(:,i,j),vctr11)
         if(co2_on)    call htint(nzp,co201dn(:,ngrid),zt,nzp,co2p(:,i,j),vctr11)

         !---------------------------------------------------------------------------------!
         !     If sounding winds are to be interpreted as eastward (U) and northward (V)   !
         ! components, rotate winds from geographic to polar stereographic orientation.    !
         !---------------------------------------------------------------------------------!
         select case (ihtran)
         case (0) !----- Cartesian. -------------------------------------------------------!
            do k = 1,nzp
               uc(k,i,j) = vctr5(k)
               vc(k,i,j) = vctr6(k)
            end do

         case (1) !----- Polar stereographic. ---------------------------------------------!
            call xy_ll(qlatu,qlonu,platn(ngrid),plonn(ngrid),xm(i),yt(j))
            call xy_ll(qlatv,qlonv,platn(ngrid),plonn(ngrid),xt(i),ym(j))

            do k = 1,nzp
               call uevetouv(uc(k,i,j),dummy,vctr5(k),vctr6(k),qlatu,qlonu                 &
                            ,platn(ngrid),plonn(ngrid))
               call uevetouv(dummy,vc(k,i,j),vctr5(k),vctr6(k),qlatv,qlonv                 &
                            ,platn(ngrid),plonn(ngrid))
            end do
         end select
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Convert virtual potential temperature into ice-liquid potential temper-     !
         ! ature.  Since the initial state has no condensates, this will be potential      !
         ! temperature.                                                                    !
         !---------------------------------------------------------------------------------!
         if (vapour_on) then
            do k=1,nzp
               thp(k,i,j) = tv2temp(vctr3(k),rtp(k,i,j))
            end do
         else
            do k=1,nzp
               thp(k,i,j) = vctr3(k)
            end do
         end if
         !---------------------------------------------------------------------------------!

         !---------------------------------------------------------------------------------!
         !    Potential temperature will be the same as thp and pc will be zero because    !
         ! there is no disturbance right at the beginning.                                 !
         !---------------------------------------------------------------------------------!
         do k=1,nzp
            theta(k,i,j) = thp(k,i,j)
            pc(k,i,j)    = 0.
         end do
         !---------------------------------------------------------------------------------!

         !---------------------------------------------------------------------------------!
         !    If there is chance for condensation, check whether this initial state is     !
         ! super-saturated or not. In case it is, put the excess as cloud droplets, and    !
         ! correct both the ice-liquid potential temperature and the water vapour mixing   !
         ! ratio.  If no condensation is allowed but water vapour is allowed, then make    !
         ! it the same as the total water mixing ratio.                                    !
         !---------------------------------------------------------------------------------!
         if (cloud_on) then
            do k=1,nzp
               p0(k)      = exner2press(pi0(k,i,j))
               temp(k)    = extheta2temp(pi0(k,i,j),theta(k,i,j))
               rvls(k)    = rslf(p0(k),temp(k))
               rc(k)      = max(0.,rtp(k,i,j)-rvls(k))
               thp(k,i,j) = theta_iceliq(pi0(k,i,j),temp(k),rc(k),0.)
               rv(k,i,j)  = rtp(k,i,j)-rc(k)
            enddo
         elseif (vapour_on) then
            do k=1,nzp
               rv(k,i,j)=rtp(k,i,j)
            end do
         end if
         !---------------------------------------------------------------------------------!

      end do iloop
   end do jloop

   return
end subroutine flds3d
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This routine initializes the 3-D velocity and thermodynamic fields from the 1-D       !
! reference state sounding when the adaptive (shaved-eta) coordinate is used.              !
!------------------------------------------------------------------------------------------!
subroutine flds3d_adap(n1,n2,n3,flpu,flpv,flpw,uc,vc,pi0,theta,thp,rtp,pc,rv,co2p)
   use mem_grid
   use ref_sounding
   use rconstants
   use mem_basic   , only : co2_on       ! ! intent(in)
   use therm_lib   , only : rslf         & ! function
                          , theta_iceliq & ! function
                          , tv2temp      & ! function
                          , exner2press  & ! function
                          , extheta2temp & ! function
                          , vapour_on    & ! intent(in)
                          , cloud_on     ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer , intent(in)                     :: n1
   integer , intent(in)                     :: n2
   integer , intent(in)                     :: n3
   real, dimension(n1,n2,n3), intent(in)    :: pi0
   real, dimension(n1,n2,n3), intent(inout) :: thp
   real, dimension(n1,n2,n3), intent(inout) :: theta
   real, dimension(n1,n2,n3), intent(out)   :: rtp
   real, dimension(n1,n2,n3), intent(out)   :: pc
   real, dimension(n1,n2,n3), intent(out)   :: rv
   real, dimension(n1,n2,n3), intent(out)   :: uc
   real, dimension(n1,n2,n3), intent(out)   :: vc
   real, dimension(n1,n2,n3), intent(out)   :: co2p
   real, dimension(   n2,n3), intent(in)    :: flpu
   real, dimension(   n2,n3), intent(in)    :: flpv 
   real, dimension(   n2,n3), intent(in)    :: flpw 
   !----- Local variables. ----------------------------------------------------------------!
   integer                                  :: i
   integer                                  :: j
   integer                                  :: k
   real                                     :: qlatu
   real                                     :: qlonu
   real                                     :: qlatv
   real                                     :: qlonv
   real                                     :: dummy
   real, dimension(nzpmax)                  :: p0
   real, dimension(nzpmax)                  :: temp
   real, dimension(nzpmax)                  :: rvls
   real, dimension(nzpmax)                  :: rc
   !---------------------------------------------------------------------------------------!

   jloop: do j = 1,n3
      iloop: do i = 1,n2

         !---------------------------------------------------------------------------------!
         !     If sounding winds are to be interpreted as eastward (U) and northward (V)   !
         ! components, rotate winds from geographic to polar stereographic orientation.    !
         !---------------------------------------------------------------------------------!
         select case (ihtran)
         case (0) !----- Cartesian. -------------------------------------------------------!
            do k = 1,n1
               uc(k,i,j) = u01dn(k,ngrid)
            end do
            do k = 1,n1
               vc(k,i,j) = v01dn(k,ngrid)
            end do

         case (1) !----- Polar stereographic. ---------------------------------------------!
            call xy_ll(qlatu,qlonu,platn(ngrid),plonn(ngrid),xm(i),yt(j))
            call xy_ll(qlatv,qlonv,platn(ngrid),plonn(ngrid),xt(i),ym(j))

            do k = 1,n1
               call uevetouv(uc(k,i,j),dummy,u01dn(k,ngrid),v01dn(k,ngrid),qlatu,qlonu     &
                            ,platn(ngrid),plonn(ngrid))
            end do
            do k = 1,n1
               call uevetouv(dummy,vc(k,i,j),u01dn(k,ngrid),v01dn(k,ngrid),qlatv,qlonv     &
                            ,platn(ngrid),plonn(ngrid))
            end do
         end select
         !---------------------------------------------------------------------------------!

         !---------------------------------------------------------------------------------!
         !    Copy CO2 mixing ratio only if CO2 is being solved.                           !
         !---------------------------------------------------------------------------------!
         if (co2_on) then
            do k=1,n1
               co2p(k,i,j) = co201dn(k,ngrid)
            end do
         end if
         !---------------------------------------------------------------------------------!

         !---------------------------------------------------------------------------------!
         !     Convert virtual potential temperature into ice-liquid potential temper-     !
         ! ature.  Since the initial state has no condensates, this will be potential      !
         ! temperature.  Also, load total mixing ratio only if water is solved.            !
         !---------------------------------------------------------------------------------!
         if (vapour_on) then
            do k = 1,n1
               rtp(k,i,j) = rt01dn(k,ngrid)
               thp(k,i,j) = tv2temp(th01dn(k,ngrid),rtp(k,i,j))
            end do
         else
            do k = 1,n1
               thp(k,i,j) = th01dn(k,ngrid)
            end do
         end if
         !---------------------------------------------------------------------------------!

         !---------------------------------------------------------------------------------!
         !    Potential temperature will be the same as thp and pc will be zero because    !
         ! there is no disturbance right at the beginning.                                 !
         !---------------------------------------------------------------------------------!
         do k = 1,n1
            theta(k,i,j) = thp(k,i,j)
            pc(k,i,j)    = 0.
         end do
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !    If there is chance for condensation, check whether this initial state is     !
         ! super-saturated or not. In case it is, put the excess as cloud droplets, and    !
         ! correct both the ice-liquid potential temperature and the water vapour mixing   !
         ! ratio.  If no condensation is allowed but water vapour is allowed, then make    !
         ! it the same as the total water mixing ratio.                                    !
         !---------------------------------------------------------------------------------!
         if (cloud_on) then
            do k = 1,n1
               p0(k)      = exner2press(pi0(k,i,j))
               temp(k)    = extheta2temp(theta(k,i,j),pi0(k,i,j))
               rvls(k)    = rslf(p0(k),temp(k))
               rc(k)      = max(0.,rtp(k,i,j) - rvls(k))
               thp(k,i,j) = theta_iceliq(pi0(k,i,j),temp(k),rc(k),0.)
               rv(k,i,j)  = rtp(k,i,j) - rc(k)
            end do
         elseif (vapour_on) then
            do k = 1,n1
               rv(k,i,j) = rtp(k,i,j)
            end do
         end if
         !---------------------------------------------------------------------------------!

      end do iloop
   end do jloop

   return
end  subroutine flds3d_adap
!==========================================================================================!
!==========================================================================================!






