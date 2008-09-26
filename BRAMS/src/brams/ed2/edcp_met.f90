!==========================================================================

subroutine copy_atm2lsm(ifm,init)

  ! This subroutine produces the following meteorologic variables
  use node_mod, only : &
       master_num,mmzp,mmxp,mmyp,  &
       ia,iz,ja,jz,ia_1,iz1,ja_1,jz1
  
  use rconstants,only:cpi,cp,p00,rocp,rgas,cliq,alli,cice,t3ple
  use met_driver_coms, only: have_co2,initial_co2
  use ed_state_vars,only: edgrid_g,edtype,polygontype
  use ed_node_coms,only:mynum
  use canopy_radiation_coms,only : rlong_min

  !------- Transfer these arrays to the polygons ----!
  use mem_basic,   only: basic_g
  use mem_radiate, only: radiate_g
  use mem_cuparm,  only: cuparm_g
  use mem_micro,   only: micro_g
  use mem_grid,    only: zt,grid_g,dzt,zm,if_adap,jdim
  use therm_lib,   only: virtt
  !--------------------------------------------------!




  implicit none

  type(edtype),pointer :: cgrid
  type(polygontype),pointer :: cpoly
  
  integer :: ipy,isi
  logical :: init
  integer :: m1,m2,m3,i,j,m1max,ifm,n
  integer :: k1w,k2w,k3w,k2u,k2u_1,k2v,k2v_1
  real, dimension(mmxp(ifm),mmyp(ifm)) :: rshortd
  real, dimension(mmxp(ifm),mmyp(ifm)) :: up_mean,vp_mean,pi0_mean,dn0_mean
  real, dimension(mmxp(ifm),mmyp(ifm)) :: rv_mean,theta_mean
  real, dimension(mmxp(ifm),mmyp(ifm)) :: map_2d_lsm
  real :: rshort1,cosz1,rshortd1,scalar1
  real :: topma_t,wtw,wtu1,wtu2,wtv1,wtv2
  integer :: ix,iy

  m2=mmxp(ifm)
  m3=mmyp(ifm)

  cgrid => edgrid_g(ifm)

  up_mean = 0.0
  vp_mean = 0.0
  pi0_mean = 0.0
  dn0_mean = 0.0
  rv_mean = 0.0
  theta_mean = 0.0


  ! Transfer atmospheric information to the ed lsm and fill those
  ! values.  Requires some mean calculations.  Also perform a 
  ! sanity check on the pressure; exit if it is unacceptable.
  !------------------------------------------------------------

  ! Use sib2 routine raddrv to calculate diffuse shortwave radiation
  !-----------------------------------------------------------------

  ! Prepare the atm data fields

  do j=ja,jz
     do i=ia,iz
        if (if_adap == 1) then
           
        ! Shaved Eta coordinate system
        !--------------------------------------
           
        k2w = nint(grid_g(ifm)%flpw(i,j))
        k1w = k2w - 1
        k3w = k2w + 1

        k2u   = nint(grid_g(ifm)%flpu(i,j))
        k2u_1 = nint(grid_g(ifm)%flpu(i-1,j))
        
        k2v   = nint(grid_g(ifm)%flpv(i,j))
        k2v_1 = nint(grid_g(ifm)%flpv(i,j-jdim))

        topma_t = .25 * (grid_g(ifm)%topma(i,j) + grid_g(ifm)%topma(i-1,j)  &
             + grid_g(ifm)%topma(i,j-jdim) + grid_g(ifm)%topma(i-1,j-jdim))
        
        ! weights for lowest predicted points, relative to points above them

        wtw = (zm(k2w) - topma_t) * dzt(k2w)
        wtu1 = grid_g(ifm)%aru(k2u_1,i-1,j)   / grid_g(ifm)%aru(k2u_1+1,i-1,j)
        wtu2 = grid_g(ifm)%aru(k2u,i,j)       / grid_g(ifm)%aru(k2u+1,i,j)
        wtv1 = grid_g(ifm)%arv(k2v_1,i,j-jdim) / grid_g(ifm)%arv(k2v_1+1,i,j-jdim)
        wtv2 = grid_g(ifm)%arv(k2v,i,j)       / grid_g(ifm)%arv(k2v+1,i,j)

        theta_mean(i,j)=  wtw * basic_g(ifm)%theta(k2w,i,j) + (1. - wtw)  * basic_g(ifm)%theta(k3w,i,j)
        
        rv_mean(i,j) =  wtw * basic_g(ifm)%rv(k2w,i,j)    + (1. - wtw)  * basic_g(ifm)%rv(k3w,i,j)
        
        up_mean(i,j) = (wtu1        * basic_g(ifm)%up(k2u_1,i-1,j)    &
             +  (1. - wtu1) * basic_g(ifm)%up(k2u_1+1,i-1,j)  &
             +  wtu2        * basic_g(ifm)%up(k2u,i,j)        &
             +  (1. - wtu2) * basic_g(ifm)%up(k2u+1,i,j)) * .5
        
        vp_mean(i,j) = (wtv1        * basic_g(ifm)%vp(k2v_1,i,j-jdim)    &
             +  (1. - wtv1) * basic_g(ifm)%vp(k2v_1+1,i,j-jdim)  &
             +  wtv2        * basic_g(ifm)%vp(k2v,i,j)          &
             +  (1. - wtv2) * basic_g(ifm)%vp(k2v+1,i,j)) * .5
        
        if (wtw >= .5) then
           pi0_mean(i,j)  = ((wtw - .5) * (basic_g(ifm)%pp(k1w,i,j) + basic_g(ifm)%pi0(k1w,i,j))  &
                + (1.5 - wtw) * (basic_g(ifm)%pp(k2w,i,j) + basic_g(ifm)%pi0(k2w,i,j)))
           dn0_mean(i,j) = (wtw - .5)  * basic_g(ifm)%dn0(k1w,i,j)  &
                + (1.5 - wtw) * basic_g(ifm)%dn0(k2w,i,j)
        else
           pi0_mean(i,j)  = ((wtw + .5) * (basic_g(ifm)%pp(k2w,i,j) + basic_g(ifm)%pi0(k2w,i,j))  &
                + (.5 - wtw) * (basic_g(ifm)%pp(k3w,i,j) + basic_g(ifm)%pi0(k3w,i,j)))
           dn0_mean(i,j) = (wtw + .5) * basic_g(ifm)%dn0(k2w,i,j)  &
                + (.5 - wtw) * basic_g(ifm)%dn0(k3w,i,j)
        endif
        

      else

        ! Terrain following coordinate system
        !--------------------------------------       
        
        theta_mean(i,j) = basic_g(ifm)%theta(2,i,j)
        rv_mean(i,j) = basic_g(ifm)%rv(2,i,j)
        up_mean(i,j) = (basic_g(ifm)%up(2,i,j) + basic_g(ifm)%up(2,i-1,j)) * 0.5
        vp_mean(i,j) = (basic_g(ifm)%vp(2,i,j) +basic_g(ifm)% vp(2,i,j-jdim)) * 0.5
        pi0_mean(i,j) = (basic_g(ifm)%pp(1,i,j) + basic_g(ifm)%pp(2,i,j) & 
             +basic_g(ifm)%pi0(1,i,j)+basic_g(ifm)%pi0(2,i,j)) * 0.5
        
        if(pi0_mean(i,j) .ne. pi0_mean(i,j))then
           print*,'bad pi0 mean'
           print*,i,j,basic_g(ifm)%pp(1,i,j),basic_g(ifm)%pp(2,i,j), &
                basic_g(ifm)%pi0(1,i,j),basic_g(ifm)%pi0(2,i,j)
           stop
        endif
        
        dn0_mean(i,j) = (basic_g(ifm)%dn0(1,i,j) + basic_g(ifm)%dn0(2,i,j)) * 0.5

      endif
     
      rshort1 = radiate_g(ifm)%rshort(i,j)
      cosz1 = radiate_g(ifm)%cosz(i,j)
     
     ! The following subroutine calculates diffuse shortwave
     ! radiation.  This is a gross approximation that uses
     ! constant scattering coefficients.  A more appropriate
     ! way of doing this was probably already designed by
     ! DMM.  Valeriy Ivanov is also investigating methods of 
     ! parameterizing diffuse shortwave radiation.
      call short2diff(rshort1,cosz1,rshortd1)
      rshortd(i,j) = rshortd1
     
    enddo
  enddo

  ! Project them to the lan-surface sites

  do ipy=1,cgrid%npolygons
     
     ix =  cgrid%ilon(ipy)
     iy =  cgrid%ilat(ipy)

     cgrid%met(ipy)%nir_beam     =  0.45*(radiate_g(ifm)%rshort(ix,iy)-rshortd(ix,iy))
     cgrid%met(ipy)%nir_diffuse  =  0.45*(rshortd(ix,iy))
     cgrid%met(ipy)%par_beam     =  0.55*(radiate_g(ifm)%rshort(ix,iy)-rshortd(ix,iy))
     cgrid%met(ipy)%par_diffuse  =  0.55*(rshortd(ix,iy))

     cgrid%met(ipy)%rlong    = radiate_g(ifm)%rlong(ix,iy)
     cgrid%met(ipy)%prss     = pi0_mean(ix,iy)*100.0    ! Convert from millibars -> Pascals
     cgrid%met(ipy)%geoht    = zt(2) + grid_g(ifm)%rtgt(ix,iy)
     cgrid%met(ipy)%vels     = sqrt( up_mean(ix,iy)**2 + vp_mean(ix,iy)**2)

     cgrid%met(ipy)%atm_shv  = rv_mean(ix,iy)
     cgrid%met(ipy)%atm_tmp  = theta_mean(ix,iy)*pi0_mean(ix,iy) * cpi
     cgrid%met(ipy)%atm_co2  = initial_co2
  enddo

  call fill_site_precip(ifm,cgrid,m2,m3,ia,iz,ja,jz,pi0_mean,theta_mean)

  
  do ipy = 1,cgrid%npolygons
     
     ! co2
     if (.not.have_co2) then
        cgrid%met(ipy)%atm_co2 = initial_co2
     end if

     ! if this is the first call, it is only looking for atm_tmp and atm_shv, so we can help
     ! bypass the error check
     
     if (init) then
           cgrid%met(ipy)%rlong=rlong_min
     endif


     !! APPLY MET TO SITES
     !! AND ADJUST MET VARIABLES FOR TOPOGRAPHY
     call calc_met_lapse_ar(cgrid,ipy)

     cpoly => cgrid%polygon(ipy)

     do isi = 1,cpoly%nsites
        
        ! vels
        cpoly%met(isi)%vels = sqrt(max(0.0,cpoly%met(isi)%vels))
        cpoly%met(isi)%vels_stab = max(0.1,cpoly%met(isi)%vels)
        cpoly%met(isi)%vels_unstab = max(1.0,cpoly%met(isi)%vels_stab)
        
        ! co2
        if (.not.have_co2) cpoly%met(isi)%atm_co2 = initial_co2
        
        ! exner
        cpoly%met(isi)%exner = cp * (cpoly%met(isi)%prss / p00)**rocp
        ! solar radiation
        cpoly%met(isi)%rshort_diffuse = cpoly%met(isi)%par_diffuse +   &
             cpoly%met(isi)%nir_diffuse
        
        cpoly%met(isi)%rshort = cpoly%met(isi)%rshort_diffuse +   &
             cpoly%met(isi)%par_beam + cpoly%met(isi)%nir_beam
        
        
        ! rho
        cpoly%met(isi)%rhos = cpoly%met(isi)%prss / (rgas *   &
             virtt(cpoly%met(isi)%atm_tmp,cpoly%met(isi)%atm_shv))
        
        ! qpcpg, dpcpg
        if(cpoly%met(isi)%atm_tmp > t3ple)then
           cpoly%met(isi)%qpcpg = (cliq * (cpoly%met(isi)%atm_tmp - t3ple) +   &
                alli) * cpoly%met(isi)%pcpg
           cpoly%met(isi)%dpcpg = max(0.0, cpoly%met(isi)%pcpg * 0.001)
        else
           cpoly%met(isi)%qpcpg = cice * (cpoly%met(isi)%atm_tmp - t3ple) *  &
                cpoly%met(isi)%pcpg
           cpoly%met(isi)%dpcpg = max(0.0, cpoly%met(isi)%pcpg * 0.01)
        endif
        !! note: dpcpg currently never gets used
        !! snow density is calculated in the integrator (grep snowdens)
        
     enddo
     
  enddo

  return
end subroutine copy_atm2lsm

!============================================================================

subroutine fill_site_precip(ifm,cgrid,m2,m3,ia,iz,ja,jz,pi0_mean,theta_mean)
  
  use mem_cuparm, only: cuparm_g,nnqparm
  use mem_micro,  only: micro_g
  use therm_lib,  only: bulk_on
  use mem_basic,  only: basic_g
  use rconstants, only: cpi, cliq
  use ed_state_vars,only: edtype
  use mem_edcp, only : ed_precip_g
  use misc_coms, only: dtlsm
  
 
  implicit none
  type(edtype),target :: cgrid
  integer :: ifm,m2,m3,i,j,ia,iz,ja,jz,n,lsmnest
  real,dimension(m2,m3) :: pi0_mean,theta_mean,conprr_mean,bulkprr_mean
  real    :: totpcp
  integer :: ix,iy
  integer :: ipy
  logical :: cumulus_on, bulkmicro_on
  real    :: dtlsmi

  cumulus_on   = nnqparm(ifm) > 0
  dtlsmi = 1./dtlsm

  ! Zero the precipitation structures
  ! --------------------------------
  do ipy=1,cgrid%npolygons
     cgrid%met(ipy)%pcpg = 0
  enddo

  do i=ia,iz
     do j=ja,jz
        conprr_mean(i,j)  = 0.
        bulkprr_mean(i,j) = 0.
     end do
  end do

  ! Depending on the microphysics being implemented
  ! Some precipitation arrays may or may not exist
  ! -----------------------------------------------

  ! 1) Check for convective precipitation
  ! -----------------------------------------------

  if (nnqparm(ifm) > 0) then

     ! We need to use difference in aconpr instead of conprr because cumulus parametrisation is usually called 
     ! less often than ED. Therefore, using conprr would repeat the same precipitation rate until cumulus is 
     ! called again, which would make ED overestimate the total convective precipitation. Using aconpr and 
     ! comparing to the previous one will give the average precipitation rate between two ED calls, so it will 
     ! take the right total amount of convective precipitation. 
     do i=ia,iz
        do j=ja,jz
           conprr_mean(i,j) = (cuparm_g(ifm)%aconpr(i,j) - ed_precip_g(ifm)%prev_aconpr(i,j)) * dtlsmi
          ed_precip_g(ifm)%prev_aconpr(i,j) = cuparm_g(ifm)%aconpr(i,j)
        end do
     end do

     ! 2) Check for resolved precipitation
     ! -----------------------------------
     if (bulk_on) then

        ! Again, we need to use the accumulated precipitation because the bulk microphysics is called every
        ! BRAMS time step. Therefore, if we use pcpg we will capture the precipitation that happened only at
        ! one time step before the current ED call, which would underestimate the resolved precipitation unless
        ! dtlsm = dtlong(ifm). Using the accumulated precipitation difference will give the amount of resolved
        ! precipitation since the last ED call, thus taking the right amount of resolved precipitation.
        do i=ia,iz
           do j=ja,jz
              !------------------------------------------------------------------------------------!
              !    Here I add all forms of precipitation, previously converted into the equivalent !
              ! amount of liquid water. This can probably be done better, especially if snow is    !
              ! falling, by averaging pcpg inside the microphysics scheme, but I will leave like   !
              ! this for now even because I don't think ED currently distinguishes between rain    !
              ! and snow...                                                                        !
              !------------------------------------------------------------------------------------!
              totpcp = micro_g(ifm)%accpr(i,j) + micro_g(ifm)%accpp(i,j) &
                     + micro_g(ifm)%accps(i,j) + micro_g(ifm)%accpa(i,j) &
                     + micro_g(ifm)%accpg(i,j) + micro_g(ifm)%accph(i,j)
              bulkprr_mean(i,j) = (totpcp-ed_precip_g(ifm)%prev_abulkpr(i,j)) * dtlsmi
              ed_precip_g(ifm)%prev_abulkpr(i,j) = totpcp
           end do
        end do
     end if
  end if

  do ipy=1,cgrid%npolygons
     ix =  cgrid%ilon(ipy)
     iy =  cgrid%ilat(ipy)
     cgrid%met(ipy)%pcpg = conprr_mean(ix,iy) + bulkprr_mean(ix,iy)
  end do

  return
end subroutine fill_site_precip

!----------------------------------------------------------------------------------------!
! Subroutine to copy the old "future" flux field to the "past"                           !
!----------------------------------------------------------------------------------------!

subroutine copy_fluxes_future_2_past(ifm)

  use mem_edcp,only: &
       ed_fluxf_g,   &
       ed_fluxp_g,   &
       wgridp_g,     &
       wgridf_g
                          
  implicit none
  integer :: ifm

  ! Do terrestrial fluxes
  ed_fluxp_g(ifm)%ustar   = ed_fluxf_g(ifm)%ustar
  ed_fluxp_g(ifm)%tstar   = ed_fluxf_g(ifm)%tstar
  ed_fluxp_g(ifm)%rstar   = ed_fluxf_g(ifm)%rstar
  ed_fluxp_g(ifm)%albedt  = ed_fluxf_g(ifm)%albedt
  ed_fluxp_g(ifm)%rlongup = ed_fluxf_g(ifm)%rlongup
  ed_fluxp_g(ifm)%sflux_u = ed_fluxf_g(ifm)%sflux_u
  ed_fluxp_g(ifm)%sflux_v = ed_fluxf_g(ifm)%sflux_v
  ed_fluxp_g(ifm)%sflux_w = ed_fluxf_g(ifm)%sflux_w
  ed_fluxp_g(ifm)%sflux_t = ed_fluxf_g(ifm)%sflux_t
  ed_fluxp_g(ifm)%sflux_r = ed_fluxf_g(ifm)%sflux_r

  ! Do water body fluxes
  wgridp_g(ifm)%ustar   = wgridf_g(ifm)%ustar
  wgridp_g(ifm)%tstar   = wgridf_g(ifm)%tstar
  wgridp_g(ifm)%rstar   = wgridf_g(ifm)%rstar
  wgridp_g(ifm)%albedt  = wgridf_g(ifm)%albedt
  wgridp_g(ifm)%rlongup = wgridf_g(ifm)%rlongup
  wgridp_g(ifm)%sflux_u = wgridf_g(ifm)%sflux_u
  wgridp_g(ifm)%sflux_v = wgridf_g(ifm)%sflux_v
  wgridp_g(ifm)%sflux_w = wgridf_g(ifm)%sflux_w
  wgridp_g(ifm)%sflux_t = wgridf_g(ifm)%sflux_t
  wgridp_g(ifm)%sflux_r = wgridf_g(ifm)%sflux_r

  return
end subroutine copy_fluxes_future_2_past
              
! ===============================================================================


!============================================================
subroutine copy_fluxes_lsm2atm(ifm)

  use ed_state_vars,only: edgrid_g,edtype,polygontype,sitetype
  use mem_edcp,only:ed_fluxf_g,ed_flux
  use mem_grid, only: zt,grid_g,dzt,zm,if_adap,jdim
  use mem_basic,only:basic_g
  use mem_leaf,only:leaf_g
  implicit none

  type(edtype),pointer      :: cgrid
  type(polygontype),pointer :: cpoly
  type(sitetype),pointer    :: csite
  type(ed_flux),pointer     :: fluxp
  integer :: ipy,isi,ifm
  integer :: ix,iy
  real    :: up_mean,vp_mean,cosine,sine

  !average and copy values
  
  ! Set the pointers, the state variable and the future flux grid
  cgrid => edgrid_g(ifm)
  fluxp => ed_fluxf_g(ifm)
  
  do ipy=1,cgrid%npolygons
     cpoly => cgrid%polygon(ipy)

     ix = cgrid%ilon(ipy)
     iy = cgrid%ilat(ipy)

     fluxp%ustar(ix,iy)   = 0.0
     fluxp%tstar(ix,iy)   = 0.0
     fluxp%rstar(ix,iy)   = 0.0
     fluxp%sflux_u(ix,iy)   = 0.0
     fluxp%sflux_v(ix,iy)   = 0.0
     fluxp%sflux_w(ix,iy)   = 0.0
     fluxp%sflux_t(ix,iy)   = 0.0
     fluxp%sflux_r(ix,iy)   = 0.0
     fluxp%rlongup(ix,iy) = 0.0
     fluxp%albedt(ix,iy)  = 0.0
     
     ! IGNORING ADAP CONVERSION FOR WINDSPEED, THIS CALCULATION RESULTS
     ! ONLY IN GETTING A WIND DIRECTION, THE END RESULT, THE MOMENTUM FLUX
     ! IN THE U AND V DIRECTION SHOULD NOT BE SO SENSITIVE TO THIS.
     
     up_mean = (basic_g(ifm)%up(2,ix,iy) + basic_g(ifm)%up(2,ix-1,iy)) * 0.5
     vp_mean = (basic_g(ifm)%vp(2,ix,iy) + basic_g(ifm)% vp(2,ix,iy-jdim)) * 0.5
     cosine = up_mean/sqrt(up_mean**2 + vp_mean**2)
     sine   = vp_mean/sqrt(up_mean**2 + vp_mean**2)

     do isi=1,cpoly%nsites
        csite => cpoly%site(isi)
        
        fluxp%ustar(ix,iy) = fluxp%ustar(ix,iy) + cpoly%area(isi)*sum(csite%area * csite%ustar)

        fluxp%tstar(ix,iy) = fluxp%tstar(ix,iy) + cpoly%area(isi)*sum(csite%area * csite%tstar)
        
        fluxp%rstar(ix,iy) = fluxp%rstar(ix,iy) + cpoly%area(isi)*sum(csite%area * csite%rstar)

        fluxp%sflux_u(ix,iy) = fluxp%sflux_u(ix,iy) + cosine*cpoly%area(isi)*sum(csite%area * csite%upwp)

        fluxp%sflux_v(ix,iy) = fluxp%sflux_v(ix,iy) + sine*cpoly%area(isi)*sum(csite%area * csite%upwp)

        fluxp%sflux_w(ix,iy) = fluxp%sflux_w(ix,iy) + cpoly%area(isi)*sum(csite%area * csite%wpwp)

        fluxp%sflux_t(ix,iy) = fluxp%sflux_t(ix,iy) + cpoly%area(isi)*sum(csite%area * csite%tpwp)

        fluxp%sflux_r(ix,iy) = fluxp%sflux_r(ix,iy) + cpoly%area(isi)*sum(csite%area * csite%rpwp)

        fluxp%rlongup(ix,iy) = fluxp%rlongup(ix,iy) + cpoly%area(isi)*cpoly%rlongup(isi)

        fluxp%albedt(ix,iy)  = fluxp%albedt(ix,iy)  + &
             cpoly%area(isi)*0.5*(cpoly%albedo_beam(isi) + cpoly%albedo_diffuse(isi))

     enddo
  enddo

  return
end subroutine copy_fluxes_lsm2atm

!===================================================================================!

subroutine initialize_ed2leaf(ifm,mxp,myp)
  
  use mem_edcp,only:    &
       ed_fluxf_g,      &
       ed_fluxp_g,      &
       wgridp_g,        &
       wgridf_g,        &
       wgrids_g,        &
       ed_precip_g,     &
       alloc_edprecip,  &
       zero_edprecip,   &
       alloc_edflux,    &
       zero_edflux,     &
       alloc_wgrid,     &
       zero_wgrid,      &
       alloc_wgrid_s,   &
       zero_wgrid_s

  use node_mod,only:mmxp,mmyp,ia,iz,ja,jz,mynum
  use ed_work_vars,only:work_e
  use mem_leaf,only:leaf_g
  use mem_basic,only: basic_g
  use mem_grid,only: grid_g,zt,dzt,zm,if_adap,jdim
  use rconstants,only:cpi
 

  implicit none
  
  integer :: ix,iy,ifm,i,j
  integer :: mxp,myp      ! This is the size of the local work_e array
  integer :: k2w,k3w,k1w
  real    :: topma_t,wtw
  real,dimension(mmxp(ifm),mmyp(ifm)) :: theta_mean,pi0_mean,rv_mean
  

  ! Note that the mmxp and mmyp are the size of the leaf_g arrays

  ! Step 1 - Set the patch are fraction in the leaf array to the land
  ! fraction in the work array. The leaf array has the border cells and 
  ! is 2 cells larger in each dimension than the work array

  leaf_g(ifm)%patch_area(:,:,1)=1.0
  leaf_g(ifm)%patch_area(:,:,2)=0.0

  do i=1,mxp
     do j=1,myp
        ix = work_e(ifm)%xatm(i,j)
        iy = work_e(ifm)%yatm(i,j)
        leaf_g(ifm)%patch_area(ix,iy,1) = 1.0-work_e(ifm)%landfrac(i,j)
        leaf_g(ifm)%patch_area(ix,iy,2) = work_e(ifm)%landfrac(i,j)
     enddo
  enddo

  ! Step 2 - allocate the flux arrays from terrestrial and water bodies
  ! These arrays are the same size as the leaf arrays, and 2 greater than
  ! the work arrays.



  call alloc_edflux(ed_fluxf_g(ifm),mmxp(ifm),mmyp(ifm))
  call alloc_edflux(ed_fluxp_g(ifm),mmxp(ifm),mmyp(ifm))

  call zero_edflux(ed_fluxf_g(ifm))
  call zero_edflux(ed_fluxp_g(ifm))

  call alloc_wgrid(wgridp_g(ifm),mmxp(ifm),mmyp(ifm))
  call alloc_wgrid(wgridf_g(ifm),mmxp(ifm),mmyp(ifm))
  
  call zero_wgrid(wgridp_g(ifm))
  call zero_wgrid(wgridf_g(ifm))

  call alloc_wgrid_s(wgrids_g(ifm),mmxp(ifm),mmyp(ifm))
  call zero_wgrid_s(wgrids_g(ifm))
  

  call alloc_edprecip(ed_precip_g(ifm),mmxp(ifm),mmyp(ifm))
  call zero_edprecip (ed_precip_g(ifm))

  do j=ja,jz
     do i=ia,iz
        if (if_adap == 1) then
           
        ! Shaved Eta coordinate system
        !--------------------------------------
        k2w = nint(grid_g(ifm)%flpw(i,j))   !CHECK
        k1w = k2w - 1
        k3w = k2w + 1                       !CHECK
        topma_t = .25 * (grid_g(ifm)%topma(i,j) + grid_g(ifm)%topma(i-1,j)  &
             + grid_g(ifm)%topma(i,j-jdim) + grid_g(ifm)%topma(i-1,j-jdim))
        wtw = (zm(k2w) - topma_t) * dzt(k2w)   !CHECK

        theta_mean(i,j)=  wtw * basic_g(ifm)%theta(k2w,i,j) + (1. - wtw)  * basic_g(ifm)%theta(k3w,i,j)
        rv_mean(i,j) =  wtw * basic_g(ifm)%rv(k2w,i,j)    + (1. - wtw)  * basic_g(ifm)%rv(k3w,i,j)

        if (wtw >= .5) then
           pi0_mean(i,j)  = ((wtw - .5) * (basic_g(ifm)%pp(k1w,i,j) + basic_g(ifm)%pi0(k1w,i,j))  &
                + (1.5 - wtw) * (basic_g(ifm)%pp(k2w,i,j) + basic_g(ifm)%pi0(k2w,i,j)))
        else
           pi0_mean(i,j)  = ((wtw + .5) * (basic_g(ifm)%pp(k2w,i,j) + basic_g(ifm)%pi0(k2w,i,j))  &
                + (.5 - wtw) * (basic_g(ifm)%pp(k3w,i,j) + basic_g(ifm)%pi0(k3w,i,j)))
        endif

      else
         
         ! Terrain following coordinate system
         !--------------------------------------       
         
         theta_mean(i,j) = basic_g(ifm)%theta(2,i,j)
         rv_mean(i,j) = basic_g(ifm)%rv(2,i,j)
         pi0_mean(i,j) = (basic_g(ifm)%pp(1,i,j) + basic_g(ifm)%pp(2,i,j) & 
              +basic_g(ifm)%pi0(1,i,j)+basic_g(ifm)%pi0(2,i,j)) * 0.5
      endif
      
      wgrids_g(ifm)%canopy_tempk(i,j)       =  theta_mean(i,j) * pi0_mean(i,j) * cpi
      wgrids_g(ifm)%canopy_water_vapor(i,j) =  rv_mean(i,j)
    enddo
  enddo

!  wgrids_g(ifm)%canopy_tempk       =  theta_mean * pi0_mean * cpi
  
!  wgrids_g(ifm)%canopy_water_vapor =  rv_mean


  return
end subroutine initialize_ed2leaf

!===================================================================================!

subroutine transfer_ed2leaf(ifm,timel)
  
  use mem_edcp,only:  &
       ed_fluxf_g,    &
       ed_fluxp_g,    &
       wgridp_g,      &
       wgridf_g,      &
       edtime1,       &
       edtime2
  use mem_turb,only: turb_g
  use rconstants,only:stefan
  use mem_leaf,only:leaf_g
  use mem_radiate,only:radiate_g
  use node_mod, only : &
       master_num,mmzp,mmxp,mmyp,  &
       ia,iz,ja,jz,ia_1,iz1,ja_1,jz1
  use mem_grid,only:jdim

  implicit none

  integer :: ifm
  real(kind=8) :: timel
  real :: tfact
  logical,save :: first=.true.
  integer :: ic,jc,ici,jci,i,j
  integer :: m2,m3


  m2 = mmxp(ifm)
  m3 = mmyp(ifm)

  
  ! There are effectively two patches land and water
  ! We could also assign leaf_g(ifm)%zrough based on the surface
  ! characteristics of ED, and discretize the sites or ed-patches
  ! into leaf patches, but for now they are left homogeneous as land
  
  tfact = (timel-edtime1)/(edtime2-edtime1)
  
  ! Interpolate the terrestrial fluxes
  leaf_g(ifm)%ustar(ia:iz,ja:jz,2) = (1-tfact)*ed_fluxp_g(ifm)%ustar(ia:iz,ja:jz) + tfact*ed_fluxf_g(ifm)%ustar(ia:iz,ja:jz)
  leaf_g(ifm)%tstar(ia:iz,ja:jz,2) = (1-tfact)*ed_fluxp_g(ifm)%tstar(ia:iz,ja:jz) + tfact*ed_fluxf_g(ifm)%tstar(ia:iz,ja:jz)
  leaf_g(ifm)%rstar(ia:iz,ja:jz,2) = (1-tfact)*ed_fluxp_g(ifm)%rstar(ia:iz,ja:jz) + tfact*ed_fluxf_g(ifm)%rstar(ia:iz,ja:jz)     

  ! Interpolate the water-body fluxes
  leaf_g(ifm)%ustar(ia:iz,ja:jz,1) = (1-tfact)*wgridp_g(ifm)%ustar(ia:iz,ja:jz) + tfact*wgridf_g(ifm)%ustar(ia:iz,ja:jz)
  leaf_g(ifm)%tstar(ia:iz,ja:jz,1) = (1-tfact)*wgridp_g(ifm)%tstar(ia:iz,ja:jz) + tfact*wgridf_g(ifm)%tstar(ia:iz,ja:jz)
  leaf_g(ifm)%rstar(ia:iz,ja:jz,1) = (1-tfact)*wgridp_g(ifm)%rstar(ia:iz,ja:jz) + tfact*wgridf_g(ifm)%rstar(ia:iz,ja:jz)


  ! Interpolate and blend the albedo, upwelling longwave, and turbulent fluxes
  
  radiate_g(ifm)%albedt(ia:iz,ja:jz) = &
       leaf_g(ifm)%patch_area(ia:iz,ja:jz,2)* &
       ((1-tfact)*ed_fluxp_g(ifm)%albedt(ia:iz,ja:jz) + tfact*ed_fluxf_g(ifm)%albedt(ia:iz,ja:jz)) + &
       leaf_g(ifm)%patch_area(ia:iz,ja:jz,1)* &
       ((1-tfact)*wgridp_g(ifm)%albedt(ia:iz,ja:jz) + tfact*wgridf_g(ifm)%albedt(ia:iz,ja:jz))
  
  radiate_g(ifm)%rlongup(ia:iz,ja:jz) = &
       leaf_g(ifm)%patch_area(ia:iz,ja:jz,2)* &
       ((1-tfact)*ed_fluxp_g(ifm)%rlongup(ia:iz,ja:jz) + tfact*ed_fluxf_g(ifm)%rlongup(ia:iz,ja:jz)) + &
       leaf_g(ifm)%patch_area(ia:iz,ja:jz,1)* &
       ((1-tfact)*wgridp_g(ifm)%rlongup(ia:iz,ja:jz) + tfact*wgridf_g(ifm)%rlongup(ia:iz,ja:jz))
  
  turb_g(ifm)%sflux_u(ia:iz,ja:jz) = &
       leaf_g(ifm)%patch_area(ia:iz,ja:jz,2)* &
       ((1-tfact)*ed_fluxp_g(ifm)%sflux_u(ia:iz,ja:jz) + tfact*ed_fluxf_g(ifm)%sflux_u(ia:iz,ja:jz)) + &
       leaf_g(ifm)%patch_area(ia:iz,ja:jz,1)* &
       ((1-tfact)*wgridp_g(ifm)%sflux_u(ia:iz,ja:jz) + tfact*wgridf_g(ifm)%sflux_u(ia:iz,ja:jz))

  turb_g(ifm)%sflux_v(ia:iz,ja:jz) = &
       leaf_g(ifm)%patch_area(ia:iz,ja:jz,2)* &
       ((1-tfact)*ed_fluxp_g(ifm)%sflux_v(ia:iz,ja:jz) + tfact*ed_fluxf_g(ifm)%sflux_v(ia:iz,ja:jz)) + &
       leaf_g(ifm)%patch_area(ia:iz,ja:jz,1)* &
       ((1-tfact)*wgridp_g(ifm)%sflux_v(ia:iz,ja:jz) + tfact*wgridf_g(ifm)%sflux_v(ia:iz,ja:jz))

  turb_g(ifm)%sflux_w(ia:iz,ja:jz) = &
       leaf_g(ifm)%patch_area(ia:iz,ja:jz,2)* &
       ((1-tfact)*ed_fluxp_g(ifm)%sflux_w(ia:iz,ja:jz) + tfact*ed_fluxf_g(ifm)%sflux_w(ia:iz,ja:jz)) + &
       leaf_g(ifm)%patch_area(ia:iz,ja:jz,1)* &
       ((1-tfact)*wgridp_g(ifm)%sflux_w(ia:iz,ja:jz) + tfact*wgridf_g(ifm)%sflux_w(ia:iz,ja:jz))

  turb_g(ifm)%sflux_t(ia:iz,ja:jz) = &
       leaf_g(ifm)%patch_area(ia:iz,ja:jz,2)* &
       ((1-tfact)*ed_fluxp_g(ifm)%sflux_t(ia:iz,ja:jz) + tfact*ed_fluxf_g(ifm)%sflux_t(ia:iz,ja:jz)) + &
       leaf_g(ifm)%patch_area(ia:iz,ja:jz,1)* &
       ((1-tfact)*wgridp_g(ifm)%sflux_t(ia:iz,ja:jz) + tfact*wgridf_g(ifm)%sflux_t(ia:iz,ja:jz))
  
  turb_g(ifm)%sflux_r(ia:iz,ja:jz) = &
       leaf_g(ifm)%patch_area(ia:iz,ja:jz,2)* &
       ((1-tfact)*ed_fluxp_g(ifm)%sflux_r(ia:iz,ja:jz) + tfact*ed_fluxf_g(ifm)%sflux_r(ia:iz,ja:jz)) + &
       leaf_g(ifm)%patch_area(ia:iz,ja:jz,1)* &
       ((1-tfact)*wgridp_g(ifm)%sflux_r(ia:iz,ja:jz) + tfact*wgridf_g(ifm)%sflux_r(ia:iz,ja:jz))


  ! The boundary cells of these arrays have not been filled.  These must be filled by the adjacent
  ! cells, likely the 2nd or 2nd to last cell in each row or column
  ! ----------------------------------------------------------------------------------------------


  do j = 1,m3
     ! Top Row
     radiate_g(ifm)%albedt(1,j) = radiate_g(ifm)%albedt(2,j)
     radiate_g(ifm)%rlongup(1,j)= radiate_g(ifm)%rlongup(2,j)
     turb_g(ifm)%sflux_u(1,j)= turb_g(ifm)%sflux_u(2,j)
     turb_g(ifm)%sflux_v(1,j)= turb_g(ifm)%sflux_v(2,j)
     turb_g(ifm)%sflux_w(1,j)= turb_g(ifm)%sflux_w(2,j)
     turb_g(ifm)%sflux_t(1,j)= turb_g(ifm)%sflux_t(2,j)
     turb_g(ifm)%sflux_r(1,j)= turb_g(ifm)%sflux_r(2,j)
     
     leaf_g(ifm)%ustar(1,j,1)=leaf_g(ifm)%ustar(2,j,1)
     leaf_g(ifm)%rstar(1,j,1)=leaf_g(ifm)%rstar(2,j,1)
     leaf_g(ifm)%tstar(1,j,1)=leaf_g(ifm)%tstar(2,j,1)
     leaf_g(ifm)%ustar(1,j,2)=leaf_g(ifm)%ustar(2,j,2)
     leaf_g(ifm)%rstar(1,j,2)=leaf_g(ifm)%rstar(2,j,2)
     leaf_g(ifm)%tstar(1,j,2)=leaf_g(ifm)%tstar(2,j,2)

     !Bottom Row
     radiate_g(ifm)%albedt(m2,j) = radiate_g(ifm)%albedt(m2-1,j)
     radiate_g(ifm)%rlongup(m2,j)= radiate_g(ifm)%rlongup(m2-1,j)
     turb_g(ifm)%sflux_u(m2,j)= turb_g(ifm)%sflux_u(m2-1,j)
     turb_g(ifm)%sflux_v(m2,j)= turb_g(ifm)%sflux_v(m2-1,j)
     turb_g(ifm)%sflux_w(m2,j)= turb_g(ifm)%sflux_w(m2-1,j)
     turb_g(ifm)%sflux_t(m2,j)= turb_g(ifm)%sflux_t(m2-1,j)
     turb_g(ifm)%sflux_r(m2,j)= turb_g(ifm)%sflux_r(m2-1,j)
     
     leaf_g(ifm)%ustar(m2,j,1)=leaf_g(ifm)%ustar(m2-1,j,1)
     leaf_g(ifm)%rstar(m2,j,1)=leaf_g(ifm)%rstar(m2-1,j,1)
     leaf_g(ifm)%tstar(m2,j,1)=leaf_g(ifm)%tstar(m2-1,j,1)
     leaf_g(ifm)%ustar(m2,j,2)=leaf_g(ifm)%ustar(m2-1,j,2)
     leaf_g(ifm)%rstar(m2,j,2)=leaf_g(ifm)%rstar(m2-1,j,2)
     leaf_g(ifm)%tstar(m2,j,2)=leaf_g(ifm)%tstar(m2-1,j,2)

  enddo
  
  if (jdim == 1) then
     do i = 1,m2

        radiate_g(ifm)%albedt(i,1) = radiate_g(ifm)%albedt(i,2)
        radiate_g(ifm)%rlongup(i,1)= radiate_g(ifm)%rlongup(i,2)
        turb_g(ifm)%sflux_u(i,1)= turb_g(ifm)%sflux_u(i,2)
        turb_g(ifm)%sflux_v(i,1)= turb_g(ifm)%sflux_v(i,2)
        turb_g(ifm)%sflux_w(i,1)= turb_g(ifm)%sflux_w(i,2)
        turb_g(ifm)%sflux_t(i,1)= turb_g(ifm)%sflux_t(i,2)
        turb_g(ifm)%sflux_r(i,1)= turb_g(ifm)%sflux_r(i,2)
        
        leaf_g(ifm)%ustar(i,1,1)=leaf_g(ifm)%ustar(i,2,1)
        leaf_g(ifm)%rstar(i,1,1)=leaf_g(ifm)%rstar(i,2,1)
        leaf_g(ifm)%tstar(i,1,1)=leaf_g(ifm)%tstar(i,2,1)
        leaf_g(ifm)%ustar(i,1,2)=leaf_g(ifm)%ustar(i,2,2)
        leaf_g(ifm)%rstar(i,1,2)=leaf_g(ifm)%rstar(i,2,2)
        leaf_g(ifm)%tstar(i,1,2)=leaf_g(ifm)%tstar(i,2,2)


        radiate_g(ifm)%albedt(i,m3) = radiate_g(ifm)%albedt(i,m3-1)
        radiate_g(ifm)%rlongup(i,m3)= radiate_g(ifm)%rlongup(i,m3-1)
        turb_g(ifm)%sflux_u(i,m3)= turb_g(ifm)%sflux_u(i,m3-1)
        turb_g(ifm)%sflux_v(i,m3)= turb_g(ifm)%sflux_v(i,m3-1)
        turb_g(ifm)%sflux_w(i,m3)= turb_g(ifm)%sflux_w(i,m3-1)
        turb_g(ifm)%sflux_t(i,m3)= turb_g(ifm)%sflux_t(i,m3-1)
        turb_g(ifm)%sflux_r(i,m3)= turb_g(ifm)%sflux_r(i,m3-1)
        
        leaf_g(ifm)%ustar(i,m3,1)=leaf_g(ifm)%ustar(i,m3-1,1)
        leaf_g(ifm)%rstar(i,m3,1)=leaf_g(ifm)%rstar(i,m3-1,1)
        leaf_g(ifm)%tstar(i,m3,1)=leaf_g(ifm)%tstar(i,m3-1,1)
        leaf_g(ifm)%ustar(i,m3,2)=leaf_g(ifm)%ustar(i,m3-1,2)
        leaf_g(ifm)%rstar(i,m3,2)=leaf_g(ifm)%rstar(i,m3-1,2)
        leaf_g(ifm)%tstar(i,m3,2)=leaf_g(ifm)%tstar(i,m3-1,2)

     enddo
  endif


  ! Fill the four corners of each grid
  ! ----------------------------------

  ic=1
  jc=1
  ici=2
  jci=2
  radiate_g(ifm)%albedt(ic,jc) = radiate_g(ifm)%albedt(ici,jci)
  radiate_g(ifm)%rlongup(ic,jc)= radiate_g(ifm)%rlongup(ici,jci)
  turb_g(ifm)%sflux_u(ic,jc)= turb_g(ifm)%sflux_u(ici,jci)
  turb_g(ifm)%sflux_v(ic,jc)= turb_g(ifm)%sflux_v(ici,jci)
  turb_g(ifm)%sflux_w(ic,jc)= turb_g(ifm)%sflux_w(ici,jci)
  turb_g(ifm)%sflux_t(ic,jc)= turb_g(ifm)%sflux_t(ici,jci)
  turb_g(ifm)%sflux_r(ic,jc)= turb_g(ifm)%sflux_r(ici,jci)
  leaf_g(ifm)%ustar(ic,jc,1)=leaf_g(ifm)%ustar(ici,jci,1)
  leaf_g(ifm)%rstar(ic,jc,1)=leaf_g(ifm)%rstar(ici,jci,1)
  leaf_g(ifm)%tstar(ic,jc,1)=leaf_g(ifm)%tstar(ici,jci,1)
  leaf_g(ifm)%ustar(ic,jc,2)=leaf_g(ifm)%ustar(ici,jci,2)
  leaf_g(ifm)%rstar(ic,jc,2)=leaf_g(ifm)%rstar(ici,jci,2)
  leaf_g(ifm)%tstar(ic,jc,2)=leaf_g(ifm)%tstar(ici,jci,2)

  ic=m2
  jc=1
  ici=m2-1
  jci=2
  radiate_g(ifm)%albedt(ic,jc) = radiate_g(ifm)%albedt(ici,jci)
  radiate_g(ifm)%rlongup(ic,jc)= radiate_g(ifm)%rlongup(ici,jci)
  turb_g(ifm)%sflux_u(ic,jc)= turb_g(ifm)%sflux_u(ici,jci)
  turb_g(ifm)%sflux_v(ic,jc)= turb_g(ifm)%sflux_v(ici,jci)
  turb_g(ifm)%sflux_w(ic,jc)= turb_g(ifm)%sflux_w(ici,jci)
  turb_g(ifm)%sflux_t(ic,jc)= turb_g(ifm)%sflux_t(ici,jci)
  turb_g(ifm)%sflux_r(ic,jc)= turb_g(ifm)%sflux_r(ici,jci)
  leaf_g(ifm)%ustar(ic,jc,1)=leaf_g(ifm)%ustar(ici,jci,1)
  leaf_g(ifm)%rstar(ic,jc,1)=leaf_g(ifm)%rstar(ici,jci,1)
  leaf_g(ifm)%tstar(ic,jc,1)=leaf_g(ifm)%tstar(ici,jci,1)
  leaf_g(ifm)%ustar(ic,jc,2)=leaf_g(ifm)%ustar(ici,jci,2)
  leaf_g(ifm)%rstar(ic,jc,2)=leaf_g(ifm)%rstar(ici,jci,2)
  leaf_g(ifm)%tstar(ic,jc,2)=leaf_g(ifm)%tstar(ici,jci,2)
  
  ic=m2
  jc=m3
  ici=m2-1
  jci=m3-1
  radiate_g(ifm)%albedt(ic,jc) = radiate_g(ifm)%albedt(ici,jci)
  radiate_g(ifm)%rlongup(ic,jc)= radiate_g(ifm)%rlongup(ici,jci)
  turb_g(ifm)%sflux_u(ic,jc)= turb_g(ifm)%sflux_u(ici,jci)
  turb_g(ifm)%sflux_v(ic,jc)= turb_g(ifm)%sflux_v(ici,jci)
  turb_g(ifm)%sflux_w(ic,jc)= turb_g(ifm)%sflux_w(ici,jci)
  turb_g(ifm)%sflux_t(ic,jc)= turb_g(ifm)%sflux_t(ici,jci)
  turb_g(ifm)%sflux_r(ic,jc)= turb_g(ifm)%sflux_r(ici,jci)
  leaf_g(ifm)%ustar(ic,jc,1)=leaf_g(ifm)%ustar(ici,jci,1)
  leaf_g(ifm)%rstar(ic,jc,1)=leaf_g(ifm)%rstar(ici,jci,1)
  leaf_g(ifm)%tstar(ic,jc,1)=leaf_g(ifm)%tstar(ici,jci,1)
  leaf_g(ifm)%ustar(ic,jc,2)=leaf_g(ifm)%ustar(ici,jci,2)
  leaf_g(ifm)%rstar(ic,jc,2)=leaf_g(ifm)%rstar(ici,jci,2)
  leaf_g(ifm)%tstar(ic,jc,2)=leaf_g(ifm)%tstar(ici,jci,2)
  
  ic=1
  jc=m3
  ici=2
  jci=m3-1
  radiate_g(ifm)%albedt(ic,jc) = radiate_g(ifm)%albedt(ici,jci)
  radiate_g(ifm)%rlongup(ic,jc)= radiate_g(ifm)%rlongup(ici,jci)
  turb_g(ifm)%sflux_u(ic,jc)= turb_g(ifm)%sflux_u(ici,jci)
  turb_g(ifm)%sflux_v(ic,jc)= turb_g(ifm)%sflux_v(ici,jci)
  turb_g(ifm)%sflux_w(ic,jc)= turb_g(ifm)%sflux_w(ici,jci)
  turb_g(ifm)%sflux_t(ic,jc)= turb_g(ifm)%sflux_t(ici,jci)
  turb_g(ifm)%sflux_r(ic,jc)= turb_g(ifm)%sflux_r(ici,jci)
  leaf_g(ifm)%ustar(ic,jc,1)=leaf_g(ifm)%ustar(ici,jci,1)
  leaf_g(ifm)%rstar(ic,jc,1)=leaf_g(ifm)%rstar(ici,jci,1)
  leaf_g(ifm)%tstar(ic,jc,1)=leaf_g(ifm)%tstar(ici,jci,1)
  leaf_g(ifm)%ustar(ic,jc,2)=leaf_g(ifm)%ustar(ici,jci,2)
  leaf_g(ifm)%rstar(ic,jc,2)=leaf_g(ifm)%rstar(ici,jci,2)
  leaf_g(ifm)%tstar(ic,jc,2)=leaf_g(ifm)%tstar(ici,jci,2)


  return
end subroutine transfer_ed2leaf

!===================================================================================!


subroutine calc_met_lapse_ar(cgrid,ipy)
  
  use ed_state_vars,only    : edtype,polygontype
  use canopy_radiation_coms,only : rlong_min
  
  implicit none
  integer, intent(in) :: ipy
  logical,parameter :: bypass=.true.
  integer :: isi
  type(edtype),target       :: cgrid
  type(polygontype),pointer :: cpoly
  
  real            :: ebar   !! mean elevation
  real            :: delE   !! deviation from mean elevation
  real            :: aterr  !! terrestrial area
  real, parameter :: offset=tiny(1.)/epsilon(1.) !! Tiny offset to avoid FPE

  !pass over sites once to calc preliminary stats

  cpoly => cgrid%polygon(ipy)

  ebar = 0.0
  aterr = 0.0

  do isi=1,cpoly%nsites
     ebar = ebar + cpoly%area(isi)*cpoly%elevation(isi)
     aterr = aterr + cpoly%area(isi)
  enddo
  ebar = ebar/aterr

  if (bypass) then
     do isi = 1,cpoly%nsites
        cpoly%met(isi)%geoht       = cgrid%met(ipy)%geoht
        cpoly%met(isi)%atm_tmp     = cgrid%met(ipy)%atm_tmp
        cpoly%met(isi)%atm_shv     = cgrid%met(ipy)%atm_shv
        cpoly%met(isi)%prss        = cgrid%met(ipy)%prss
        cpoly%met(isi)%pcpg        = cgrid%met(ipy)%pcpg
        cpoly%met(isi)%par_diffuse = cgrid%met(ipy)%par_diffuse
        cpoly%met(isi)%atm_co2     = cgrid%met(ipy)%atm_co2
        cpoly%met(isi)%rlong       = cgrid%met(ipy)%rlong
        cpoly%met(isi)%par_beam    = cgrid%met(ipy)%par_beam
        cpoly%met(isi)%nir_diffuse = cgrid%met(ipy)%nir_diffuse
        cpoly%met(isi)%nir_beam    = cgrid%met(ipy)%nir_beam
        cpoly%met(isi)%vels        = cgrid%met(ipy)%vels
      
        if ( cpoly%met(isi)%rlong < rlong_min) then
           print*,"Problems with RLONG A"  
           print*,cpoly%met(isi)%rlong,cgrid%met(ipy)%rlong,rlong_min
           call fatal_error('Problems with RLONG A','calc_met_lapse_ar','ed_met_driver.f90')
        end if
        if ( cpoly%met(isi)%atm_tmp < 150.0) then
           print*,cpoly%met(isi)%atm_tmp,cgrid%met(ipy)%atm_tmp
           call fatal_error('Problems with ATM TEMP A','calc_met_lapse_ar','ed_met_driver.f90')
        endif
        if ( cpoly%met(isi)%atm_shv < 0.01e-2) then
           print*,cpoly%met(isi)%atm_shv,cgrid%met(ipy)%atm_shv
           call fatal_error('Problems with ATM MOISTURE A','calc_met_lapse_ar','ed_met_driver.f90')
        endif
        if ( cpoly%met(isi)%rlong > 600.0) then
           print*,"Problems with RLONG too high"  
           print*,cpoly%met(isi)%rlong,cgrid%met(ipy)%rlong
           call fatal_error('Problems with RLONG A','calc_met_lapse_ar','ed_met_driver.f90')
        end if
        if ( cpoly%met(isi)%atm_tmp > 317.0) then
           print*,"Problems with atm temp"  
           print*,cpoly%met(isi)%atm_tmp,cgrid%met(ipy)%atm_tmp
           call fatal_error('Problems with atm temp','calc_met_lapse_ar','ed_met_driver.f90')
        end if
        if ( cpoly%met(isi)%atm_tmp < 130.0) then
           print*,"Problems with atm temp,too low"  
           print*,cpoly%met(isi)%atm_tmp,cgrid%met(ipy)%atm_tmp
           call fatal_error('Problems with atm temp','calc_met_lapse_ar','ed_met_driver.f90')
        end if
        if ( cpoly%met(isi)%par_beam   +cpoly%met(isi)%nir_beam+ &
             cpoly%met(isi)%par_diffuse+cpoly%met(isi)%nir_diffuse > 1320.0 ) then
           print*,"Problems with solar radiation, higher than solar const."
           print*,cpoly%met(isi)%par_beam,cpoly%met(isi)%nir_beam
           print*,cpoly%met(isi)%par_diffuse,cpoly%met(isi)%nir_diffuse
           call fatal_error('Problems with solar radiation','calc_met_lapse_ar','ed_met_driver.f90')
        end if
        if ( cpoly%met(isi)%atm_shv < 0.001e-3) then
           print*,"Atmospheric moisture at surface too low"  
           print*,cpoly%met(isi)%atm_shv
           call fatal_error('Problems with atm shv','calc_met_lapse_ar','ed_met_driver.f90')
        end if

        if ( cpoly%met(isi)%atm_shv > 30.0e-3) then
           print*,"Atmospheric moisture at surface too high"  
           print*,cpoly%met(isi)%atm_shv
           call fatal_error('Problems with atm shv','calc_met_lapse_ar','ed_met_driver.f90')
        end if
     enddo
   
  else
     
     !!second pass, calc lapse rate adjustment
     do isi = 1,cpoly%nsites
        
        delE = cpoly%elevation(isi) - ebar
        
        !! perform linear adjustments
        cpoly%met(isi)%geoht   = cgrid%met(ipy)%geoht   + cgrid%lapse(ipy)%geoht*delE
        cpoly%met(isi)%atm_tmp = cgrid%met(ipy)%atm_tmp + cgrid%lapse(ipy)%atm_tmp*delE
        cpoly%met(isi)%atm_shv = cgrid%met(ipy)%atm_shv + cgrid%lapse(ipy)%atm_shv*delE
        cpoly%met(isi)%prss    = cgrid%met(ipy)%prss    + cgrid%lapse(ipy)%prss*delE
        cpoly%met(isi)%pcpg    = cgrid%met(ipy)%pcpg    + cgrid%lapse(ipy)%pcpg*delE
        cpoly%met(isi)%atm_co2 = cgrid%met(ipy)%atm_co2 + cgrid%lapse(ipy)%atm_co2*delE
        cpoly%met(isi)%rlong   = cgrid%met(ipy)%rlong   + cgrid%lapse(ipy)%rlong*delE
        cpoly%met(isi)%par_diffuse = cgrid%met(ipy)%par_diffuse + cgrid%lapse(ipy)%par_diffuse*delE
        cpoly%met(isi)%par_beam    = cgrid%met(ipy)%par_beam    + cgrid%lapse(ipy)%par_beam*delE
        cpoly%met(isi)%nir_diffuse = cgrid%met(ipy)%nir_diffuse + cgrid%lapse(ipy)%nir_diffuse*delE
        cpoly%met(isi)%nir_beam    = cgrid%met(ipy)%nir_beam    + cgrid%lapse(ipy)%nir_beam*delE
        cpoly%met(isi)%vels    = cgrid%met(ipy)%vels    + cgrid%lapse(ipy)%vels*delE
        !! note: at this point VELS is vel^2.  Thus this lapse preserves mean wind ENERGY
        !! not wind SPEED
        if ( cpoly%met(isi)%rlong < 200.0) then
           print*,"Problems with RLONG A"  
           print*,cpoly%met(isi)%rlong,cgrid%met(ipy)%rlong
           call fatal_error('Problems with RLONG A','calc_met_lapse_ar','ed_met_driver.f90')
        end if
        if ( cpoly%met(isi)%atm_tmp < 150.0) then
           print*,cpoly%met(isi)%atm_tmp,cgrid%met(ipy)%atm_tmp
           call fatal_error('Problems with ATM TEMP A','calc_met_lapse_ar','ed_met_driver.f90')
        endif
        if ( cpoly%met(isi)%atm_shv < 0.01e-2) then
           print*,cpoly%met(isi)%atm_shv,cgrid%met(ipy)%atm_shv
           call fatal_error('Problems with ATM MOISTURE A','calc_met_lapse_ar','ed_met_driver.f90')
        endif
        
     enddo
  endif
  return
end subroutine calc_met_lapse_ar

!===========================================================================!
!===========================================================================!

subroutine int_met_avg(cgrid)
  
  ! Increment the time averaged polygon met-forcing
  ! variables. THese will be normalized by the output
  ! period to give time averages of each quanity. The
  ! polygon level variables are derived from the
  ! weighted spatial average from the site level quantities.


  use ed_state_vars,only:edtype,polygontype,sitetype,patchtype
  use misc_coms,only : dtlsm,frqfast
  
  implicit none

  type(edtype), target :: cgrid
  type(polygontype),pointer :: cpoly
  type(sitetype),pointer :: csite
  type(patchtype),pointer :: cpatch

  integer :: ipy,isi,ipa,ico
  real :: frqfasti,tfact

  frqfasti = 1.0 / frqfast
  tfact = dtlsm * frqfasti

  do ipy = 1,cgrid%npolygons

     cpoly => cgrid%polygon(ipy)

     do isi = 1,cpoly%nsites
        
        cgrid%avg_nir_beam(ipy)    = cgrid%avg_nir_beam(ipy) + &
             cpoly%met(isi)%nir_beam * cpoly%area(isi) * tfact
        cgrid%avg_nir_diffuse(ipy) = cgrid%avg_nir_diffuse(ipy) + &
             cpoly%met(isi)%nir_diffuse * cpoly%area(isi) * tfact
        cgrid%avg_par_beam(ipy)    = cgrid%avg_par_beam(ipy) + &
             cpoly%met(isi)%par_beam * cpoly%area(isi) * tfact
        cgrid%avg_par_diffuse(ipy) = cgrid%avg_par_diffuse(ipy) + &
             cpoly%met(isi)%par_diffuse * cpoly%area(isi) * tfact
        cgrid%avg_atm_tmp(ipy)     = cgrid%avg_atm_tmp(ipy) + &
             cpoly%met(isi)%atm_tmp * cpoly%area(isi) * tfact
        cgrid%avg_atm_shv(ipy)     = cgrid%avg_atm_shv(ipy) + &
             cpoly%met(isi)%atm_shv * cpoly%area(isi) * tfact
        cgrid%avg_rhos(ipy)        = cgrid%avg_rhos(ipy) + &
             cpoly%met(isi)%rhos * cpoly%area(isi) * tfact
        cgrid%avg_theta(ipy)       = cgrid%avg_theta(ipy) + &
             cpoly%met(isi)%theta * cpoly%area(isi) * tfact
        cgrid%avg_rshort(ipy)      = cgrid%avg_rshort(ipy) + &
             cpoly%met(isi)%rshort * cpoly%area(isi) * tfact
        cgrid%avg_rshort_diffuse(ipy) = cgrid%avg_rshort_diffuse(ipy) + &
             cpoly%met(isi)%rshort_diffuse * cpoly%area(isi) * tfact
        cgrid%avg_rlong(ipy)       = cgrid%avg_rlong(ipy) + &
             cpoly%met(isi)%rlong * cpoly%area(isi) * tfact
        cgrid%avg_pcpg(ipy)        = cgrid%avg_pcpg(ipy) + &
             cpoly%met(isi)%pcpg * cpoly%area(isi) * tfact
        cgrid%avg_qpcpg(ipy)       = cgrid%avg_qpcpg(ipy) + &
             cpoly%met(isi)%qpcpg * cpoly%area(isi) * tfact
        cgrid%avg_dpcpg(ipy)       = cgrid%avg_dpcpg(ipy) + &
             cpoly%met(isi)%dpcpg * cpoly%area(isi) * tfact
        cgrid%avg_vels(ipy)        = cgrid%avg_vels(ipy) + &
             cpoly%met(isi)%vels * cpoly%area(isi) * tfact
        cgrid%avg_prss(ipy)        = cgrid%avg_prss(ipy) + &
             cpoly%met(isi)%prss * cpoly%area(isi) * tfact
        cgrid%avg_exner(ipy)       = cgrid%avg_exner(ipy) + &
             cpoly%met(isi)%exner * cpoly%area(isi) * tfact
        cgrid%avg_geoht(ipy)       = cgrid%avg_geoht(ipy) + &
             cpoly%met(isi)%geoht * cpoly%area(isi) * tfact
        cgrid%avg_atm_co2(ipy)     = cgrid%avg_atm_co2(ipy) + &
             cpoly%met(isi)%atm_co2 * cpoly%area(isi) * tfact
        cgrid%avg_albedt(ipy)      = cgrid%avg_albedt(ipy) + &
             0.5*(cpoly%albedo_beam(isi)+cpoly%albedo_diffuse(isi))*cpoly%area(isi) * tfact
        cgrid%avg_rlongup(ipy)     = cgrid%avg_rlongup(ipy) + &
             cpoly%rlongup(isi) * cpoly%area(isi) * tfact


     enddo
        
  enddo
  
  return
end subroutine int_met_avg
