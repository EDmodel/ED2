!==========================================================================================!
!==========================================================================================!
subroutine simple_lake_model(time,dtlongest)


  use node_mod,only:ja,jz,ia,iz

  use consts_coms,only:stefan,cpi,vonk,cp,grav,p00,rocp,cpor,alvl
  use canopy_air_coms, only : ubmin
  use mem_edcp,only:wgrid_g

  !------- Transfer these arrays to the polygons --------------!
  use io_params,  only: ssttime1,ssttime2,iupdsst
  use mem_leaf,   only: leaf_g
  use mem_basic,  only: basic_g
  use mem_radiate,only: radiate_g
  use mem_cuparm, only: cuparm_g
  use mem_micro,  only: micro_g
  use mem_grid,   only: zt,grid_g,dzt,zm,if_adap,jdim,ngrid
  use therm_lib,  only: rslif
  !------------------------------------------------------------!

  implicit none
  real, intent(in) :: dtlongest
  real(kind=8),intent(in) :: time
  real :: cosz
  real :: prss
  real :: ustar,tstar,rstar,thetacan,water_rsat,zts,water_rough
  real :: vels_pat
  real :: bot,top,last_rv,last_th
  real :: hcapcan,wcapcan
  real :: z0fac_water,pis,rdi,idt
  real,external :: vertical_vel_flux


  integer :: i,j,n
  integer :: k1w,k2w,k3w,k2u,k2u_1,k2v,k2v_1
  real :: up_mean,vp_mean,exner_mean,dn0_mean
  real :: rv_mean,theta_mean
  real :: topma_t,wtw,wtu1,wtu2,wtv1,wtv2
  real :: canopy_water_vapor
  real :: canopy_tempk
  real :: sflux_u,sflux_v,sflux_w,sflux_t,sflux_r
  integer :: niter_leaf
  
  real :: dtll_factor,dtll,dtlc_factor,dtlc
  real :: dtllohcc,dtllowcc,dtlcohcc,dtlcowcc
  real :: gzotheta, patarea,vels
  real :: timefac_sst,seatc
  real    :: ustaro,delz,d_vel,d_veln,vel_new
  integer :: ifixu
  real, parameter   :: co2_mean   = 370.00
  real, parameter   :: canopy_co2 = 370.00
  real, parameter   :: emiss_w  = 0.97  ! emissivity of water (super gross approximation!)

  real              :: cstar ! Dummy variables for now, if you have 
                             ! Actual data feel free to compute co2_mean
                             ! and canopy_co2, then cstar will have meaning.

   ! Define leaf3 and canopy time-split timesteps here.  This ensures that leaf3
   ! will not use a timestep longer than about 40 seconds, and canopy will not
   ! use a timestep longer than about 15 seconds.  This allows values of
   ! hcapcan = 2.e4, wcapcan = 2.e1, and hcapveg = 3.e4 as are now defined below.
   niter_leaf  = max(1,nint(dtlongest/40.+.4))
   dtll_factor = 1. / float(niter_leaf)
   dtll        = dtlongest * dtll_factor
  
   
   
   ! Update the sea-surface temperature
   
   if (iupdsst == 0) then
      timefac_sst = 0.
   else
      timefac_sst = sngl((time - ssttime1(ngrid)) / (ssttime2(ngrid) - ssttime1(ngrid)))
   endif

   


  ! Note: Water albedo from Atwater and Bell (1981), excerpted from the LEAF3 scheme
  
  ! The target step size for the water body suface flux is about 30 seconds

!  hcapcan = 2.0e4
!  wcapcan = 2.0e1

!  hcapcan = 2.0e4
!  wcapcan = 2.0e1

   



!  dtllohcc = dtll / hcapcan
!  dtllowcc = dtll / wcapcan

  
  z0fac_water = .016 / grav

  ! Transfer atmospheric information to the ed lsm and fill those
  ! values.  Requires some mean calculations.  Also perform a 
  ! sanity check on the pressure; exit if it is unacceptable.
  !------------------------------------------------------------


  ! Prepare the atm data fields
  
  
  do j=ja,jz
     do i=ia,iz

        seatc = leaf_g(ngrid)%seatp(i,j) + timefac_sst*(leaf_g(ngrid)%seatf(i,j) - leaf_g(ngrid)%seatp(i,j))

        if (if_adap == 1) then
           
           ! Shaved Eta coordinate system
           !--------------------------------------
           
           k2w = nint(grid_g(ngrid)%flpw(i,j))
           k1w = k2w - 1
           k3w = k2w + 1
           
           k2u   = nint(grid_g(ngrid)%flpu(i,j))
           k2u_1 = nint(grid_g(ngrid)%flpu(i-1,j))
           
           k2v   = nint(grid_g(ngrid)%flpv(i,j))
           k2v_1 = nint(grid_g(ngrid)%flpv(i,j-jdim))
           
           topma_t = .25 * (grid_g(ngrid)%topma(i,j) + grid_g(ngrid)%topma(i-1,j)  &
                + grid_g(ngrid)%topma(i,j-jdim) + grid_g(ngrid)%topma(i-1,j-jdim))
           
           ! weights for lowest predicted points, relative to points above them
           
           wtw = (zm(k2w) - topma_t) * dzt(k2w)
           wtu1 = grid_g(ngrid)%aru(k2u_1,i-1,j)   / grid_g(ngrid)%aru(k2u_1+1,i-1,j)
           wtu2 = grid_g(ngrid)%aru(k2u,i,j)       / grid_g(ngrid)%aru(k2u+1,i,j)
           wtv1 = grid_g(ngrid)%arv(k2v_1,i,j-jdim) / grid_g(ngrid)%arv(k2v_1+1,i,j-jdim)
           wtv2 = grid_g(ngrid)%arv(k2v,i,j)       / grid_g(ngrid)%arv(k2v+1,i,j)
           
           theta_mean   =  wtw * basic_g(ngrid)%theta(k2w,i,j) + (1. - wtw)  * basic_g(ngrid)%theta(k3w,i,j)
           
           rv_mean      =  wtw * basic_g(ngrid)%rv(k2w,i,j)    + (1. - wtw)  * basic_g(ngrid)%rv(k3w,i,j)
           
           up_mean      = (wtu1        * basic_g(ngrid)%up(k2u_1,i-1,j)    &
                +  (1. - wtu1) * basic_g(ngrid)%up(k2u_1+1,i-1,j)  &
                +  wtu2        * basic_g(ngrid)%up(k2u,i,j)        &
                +  (1. - wtu2) * basic_g(ngrid)%up(k2u+1,i,j)) * .5
        
           vp_mean      = (wtv1        * basic_g(ngrid)%vp(k2v_1,i,j-jdim)    &
                +  (1. - wtv1) * basic_g(ngrid)%vp(k2v_1+1,i,j-jdim)  &
                +  wtv2        * basic_g(ngrid)%vp(k2v,i,j)          &
                +  (1. - wtv2) * basic_g(ngrid)%vp(k2v+1,i,j)) * .5
           
           if (wtw >= .5) then
              exner_mean   = ((wtw - .5) * (basic_g(ngrid)%pp(k1w,i,j) + basic_g(ngrid)%pi0(k1w,i,j))  &
                   + (1.5 - wtw) * (basic_g(ngrid)%pp(k2w,i,j) + basic_g(ngrid)%pi0(k2w,i,j)))
              dn0_mean     = (wtw - .5)  * basic_g(ngrid)%dn0(k1w,i,j)  &
                   + (1.5 - wtw) * basic_g(ngrid)%dn0(k2w,i,j)
           else
              exner_mean  = ((wtw + .5) * (basic_g(ngrid)%pp(k2w,i,j) + basic_g(ngrid)%pi0(k2w,i,j))  &
                   + (.5 - wtw) * (basic_g(ngrid)%pp(k3w,i,j) + basic_g(ngrid)%pi0(k3w,i,j)))
              dn0_mean    = (wtw + .5) * basic_g(ngrid)%dn0(k2w,i,j)  &
                   + (.5 - wtw) * basic_g(ngrid)%dn0(k3w,i,j)
           endif
           
           
        else
           
           ! Terrain following coordinate system
           !--------------------------------------       
           
           theta_mean   = basic_g(ngrid)%theta(2,i,j)
           rv_mean      = basic_g(ngrid)%rv(2,i,j)
           
           up_mean       = (basic_g(ngrid)%up(2,i,j) + basic_g(ngrid)%up(2,i-1,j))     * 0.5
           vp_mean       = (basic_g(ngrid)%vp(2,i,j) + basic_g(ngrid)% vp(2,i,j-jdim)) * 0.5
           exner_mean    = (basic_g(ngrid)%pp(1,i,j) + basic_g(ngrid)%pp(2,i,j) & 
                         +  basic_g(ngrid)%pi0(1,i,j)+basic_g(ngrid)%pi0(2,i,j))       * 0.5
           
           
           dn0_mean     = (basic_g(ngrid)%dn0(1,i,j) + basic_g(ngrid)%dn0(2,i,j)) * 0.5
           
        endif


        cosz  = radiate_g(ngrid)%cosz(i,j)

        ustar              = wgrid_g(ngrid)%ustar(i,j)
        canopy_tempk       = leaf_g(ngrid)%can_temp(i,j,1)
        canopy_water_vapor = leaf_g(ngrid)%can_rvap(i,j,1)

        sflux_u = 0.0
        sflux_v = 0.0
        sflux_w = 0.0
        sflux_r = 0.0
        sflux_t = 0.0

        idt = 0
        do n = 1,niter_leaf

           pis =  exner_mean * cpi
           
           prss = pis ** cpor * p00

           water_rsat  = rslif(prss,seatc)
           
           water_rough = max(z0fac_water * ustar ** 2,.0001)
           
           thetacan = canopy_tempk / pis
           
           zts = zt(2) + grid_g(ngrid)%rtgt(i,j)

           patarea  = leaf_g(ngrid)%patch_area(i,j,1)
           vels     = sqrt(up_mean**2 + vp_mean**2)
           vels_pat = max(vels,ubmin)
           
           !---- This is the LEAF stars subroutine. ---------------------------------------!
           call ed_stars(theta_mean,rv_mean,co2_mean,thetacan,zts,vels,water_rough &
                        ,ustar,rstar,tstar,cstar,canopy_water_vapor,canopy_co2)
           
           !----- This part is on LEAF, but not on ED, don't know how necessary this is. --!
!           ifixu=0
!           ustaro=ustar
!           delz = 2.*zts
!           d_vel =  - ustar * ustar *dtlongest / delz
!           vel_new = vels_pat + d_vel
!           if (vel_new < .5 * vels_pat) then
!              ifixu=1
!              d_veln = .5 * vels_pat
!              ustar=sqrt(d_veln*delz/dtlongest)
!           end if

           ! Calculate the heat,moisture and momentum fluxes
           ! -----------------------------------------------

           sflux_u = sflux_u - ustar*ustar*up_mean/vels_pat
           sflux_v = sflux_v - ustar*ustar*vp_mean/vels_pat
           sflux_t = sflux_t - ustar*tstar
           sflux_r = sflux_r - ustar*rstar
           
           gzotheta = grav * zts / theta_mean
           
           sflux_w = sflux_w + vertical_vel_flux(gzotheta,tstar,ustar)
           
           ! Update the sea surface air temperature and water vapor mixing ratio
           ! -------------------------------------------------------------------

           ! In calculating the water capacity and heat capacity of the
           ! sea surface air space
           ! We will assume a layer that is 20 meters thick
           
           wcapcan = 20.0 * dn0_mean
           hcapcan = cp * 20.0 * dn0_mean

           dtllohcc = dtll/hcapcan
           dtllowcc = dtll/wcapcan


           rdi = .2 * ustar
           
           last_th = canopy_tempk

           canopy_tempk  = canopy_tempk        &
                + dtllohcc * dn0_mean * cp                    &
                * ( (seatc -  canopy_tempk) * rdi    &
                + ustar * tstar * pis)
           
           bot = ( water_rsat - canopy_water_vapor) * rdi
           top = ustar * rstar
           last_rv = canopy_water_vapor
           
           canopy_water_vapor = canopy_water_vapor &
                + dtllowcc * dn0_mean * (( water_rsat-canopy_water_vapor) * rdi  &
                + ustar * rstar )
           
           if(canopy_water_vapor < 0.001 .or. canopy_water_vapor /= canopy_water_vapor .or. &
              canopy_tempk /= canopy_tempk ) then
              write(unit=*,fmt='(a)') '======= WATER VAPOR IN CAS ABOVE WATER-BODIES IS SCREWY! ======'
              write(unit=*,fmt='(3(a,1x,i5,1x))') 'i=',i,'j=',j,'n=',n
              write(unit=*,fmt='(2(a,1x,f8.2,1x))') 'Lon: ',grid_g(ngrid)%glon(i,j) &
                                                   ,'Lat: ',grid_g(ngrid)%glat(i,j)
              write(unit=*,fmt=*) 'EXNER (PIO)        : ',exner_mean
              write(unit=*,fmt=*) 'DN0_MEAN           : ',dn0_mean
              write(unit=*,fmt=*) 'THETA_MEAN         : ',theta_mean
              write(unit=*,fmt=*) 'RV_MEAN            : ',rv_mean
              write(unit=*,fmt=*) 'SST                : ',seatc
              write(unit=*,fmt=*) 'CANOPY_TEMPK       : ',canopy_tempk
              write(unit=*,fmt=*) 'CANOPY_RVAP/RSAT   : ',canopy_water_vapor,water_rsat
              write(unit=*,fmt=*) 'RDI                : ',rdi
              write(unit=*,fmt=*) 'DTLLOWCC           : ',dtllowcc
              write(unit=*,fmt=*) 'USTAR              : ',ustar
              write(unit=*,fmt=*) 'TSTAR              : ',tstar
              write(unit=*,fmt=*) 'RSTAR              : ',rstar
              write(unit=*,fmt=*) 'BOTTOM             : ',bot
              write(unit=*,fmt=*) 'TOP                : ',top
              write(unit=*,fmt=*) 'LAST_RV            : ',last_rv
              write(unit=*,fmt=*) 'LAST_TH            : ',last_th

              write(unit=*,fmt='(a)') '=======FLUXES======='
              
              write(unit=*,fmt=*) 'U momentum flux:    ',wgrid_g(ngrid)%sflux_u(i,j)
              write(unit=*,fmt=*) 'V momentum flux:    ',wgrid_g(ngrid)%sflux_v(i,j)
              write(unit=*,fmt=*) 'W momentum flux:    ',wgrid_g(ngrid)%sflux_w(i,j)
              write(unit=*,fmt=*) 'Sensible Heat flux: ',wgrid_g(ngrid)%sflux_t(i,j)
              write(unit=*,fmt=*) 'Latent Heat flux:   ',alvl*wgrid_g(ngrid)%sflux_r(i,j)

              call fatal_error('Lake model failed','simple_lake_model','edcp_water.f90')
           end if
     
        enddo

        ! Transfer model scalars back to global arrays
        ! --------------------------------------------

        wgrid_g(ngrid)%ustar(i,j) = ustar
        wgrid_g(ngrid)%rstar(i,j) = rstar
        wgrid_g(ngrid)%tstar(i,j) = tstar

        wgrid_g(ngrid)%sflux_u(i,j) = dn0_mean*sflux_u/real(niter_leaf)
        wgrid_g(ngrid)%sflux_v(i,j) = dn0_mean*sflux_v/real(niter_leaf)
        wgrid_g(ngrid)%sflux_w(i,j) = dn0_mean*sflux_w/real(niter_leaf)
        wgrid_g(ngrid)%sflux_t(i,j) = dn0_mean*sflux_t/real(niter_leaf)
        wgrid_g(ngrid)%sflux_r(i,j) = dn0_mean*sflux_r/real(niter_leaf)

        wgrid_g(ngrid)%albedt(i,j) = min(max(-.0139 + .0467 * tan(acos(cosz)),.03),.999)
        wgrid_g(ngrid)%rlongup(i,j) = emiss_w * stefan * seatc**4


        leaf_g(ngrid)%can_temp(i,j,1) = canopy_tempk
        leaf_g(ngrid)%can_rvap(i,j,1) = canopy_water_vapor
        

     end do
  end do

!  print*,"LAKE MODEL",ja,ja,ia,iz
!  print*,"ustar:", wgridf_g(ngrid)%ustar(2,2)," tstar: " &
!       , wgridf_g(ngrid)%tstar(2,2)," rstar:", wgridf_g(ngrid)%rstar(2,2)
!  print*,"CANOPY TEMP:",wgrids_g(ngrid)%canopy_tempk(2,2)
!  print*,"SEA TEMP:",leaf_g(ngrid)%seatp(2,2)
!  print*,"WATER VAPOR ",wgrids_g(ngrid)%canopy_water_vapor(ia:iz,ja:jz)
!  print*,"U MOMENTUM FLUX",wgridf_g(ngrid)%sflux_u(2,2)
!  print*,"HEAT FLUX:",wgridf_g(ngrid)%sflux_t(2,2)
!  print*,"MOISTURE FLUX:",wgridf_g(ngrid)%sflux_r(2,2)

  return
end subroutine simple_lake_model
!==========================================================================================!
!==========================================================================================!
