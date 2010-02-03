! Adapted in 17/07/2002 by Alvaro L.Fazenda for version 5.0
! Adapted in 15/04/2003 by Alvaro L.Fazenda for version 5.04
! Adapted in 06/03/2007 by Alvaro L.Fazenda for BRAMS V.4.0
!-------------------------

subroutine old_grell_cupar_driver(icld)

  use mem_basic         , only: basic_g & ! intent(in)
                              , co2_on  ! ! intent(in)

  use mem_tend          , only: tend
![MLO - Adding shallow cumulus output
  use mem_cuparm        , only: nnqparm, confrq, cuparm_g,nclouds
!MLO]
  use node_mod          , only: &
       mynum, &   ! INTENT(IN)
       mxp,   &   ! INTENT(IN)
       myp,   &   ! INTENT(IN)
       mzp,   &   ! INTENT(IN)
       ia,    &   ! INTENT(IN)
       iz,    &   ! INTENT(IN)
       ja,    &   ! INTENT(IN)
       jz,    &   ! INTENT(IN)
       i0,    &   ! INTENT(IN)  ! Rever função
       j0         ! INTENT(IN)  ! Rever função

  ! USE Modules for Grell Parameterization
  use mem_grell_param   , only: &
       mgmxp,     & ! INTENT(IN)
       mgmyp,     & ! INTENT(IN)
       mgmzp,     & ! INTENT(IN)
       maxiens,   & ! INTENT(IN)
       ngrids_cp  ! ! INTENT(IN)

  use grell_coms, only:  &
       iupmethod,    & ! INTENT(IN)
       depth_min,    & ! INTENT(IN)
       cap_maxs,     & ! INTENT(IN)
       radius,       & !
       zkbmax,       & !
       zcutdown,     & !
       z_detr        ! !

  use mem_scratch1_grell    , only: sc1_grell_g
  use mem_scratch2_grell    , only: zero_scratch2_grell
  use mem_scratch2_grell_sh , only: zero_scratch2_grell_sh
  use mem_scratch3_grell    , only: zero_scratch3_grell
  use mem_scratch3_grell_sh , only: zero_scratch3_grell_sh

  use mem_grid          , only: time,   & ! INTENT(IN)
       initial,                & ! INTENT(IN)
       dtlt,                   & ! INTENT(IN)
       ngrid,                  &
       grid_g,                 &
       naddsc,                 &
       itimea                    ! INTENT(IN)

  use rconstants        , only: tkmin

  use mem_turb          , only: turb_g,idiffk

  use mem_micro         , only: micro_g

  use mem_scratch       , only: scratch

  use mem_scalar        , only: scalar_g
![MLO - Extra variables for output
  use io_params         , only: frqanl
  use mem_leaf          , only: leaf_g
  use mem_mass          , only: imassflx,mass_g
!MLO]

  implicit none
  integer, intent(in) :: icld

  !MLO. This is to use the PBL height that was found when running Nakanishi and Niino.
  !     If we didn't run Nakanishi and Niino, then use the default.
  
  

  if (icld == 1) then

     !----- Remove part of the instability due to shallow cumulus action ------------------!
     call include_shal_effect(mzp,mxp,myp,ia,iz,ja,jz,dtlt                                 &
                   ,sc1_grell_g(ngrid)%thetasta , sc1_grell_g(ngrid)%rvsta      &
                   ,basic_g(ngrid)%theta, basic_g(ngrid)%rv                     &
                   ,basic_g(ngrid)%pi0  , basic_g(ngrid)%pp                     &
                   ,cuparm_g(ngrid)%thsrc(:,:,:,nclouds)                        &
                   ,cuparm_g(ngrid)%rtsrc(:,:,:,nclouds))

     !-------------------------------------------------------------------------------------!
     !     Zero out CO2 tendency if CO2 is prognosed.  The old Grell scheme won't compute  !
     ! the transport of CO2 through updrafts and downdrafts, feel free to add this.  It    !
     ! should be similar to the water transport, except that CO2 doesn't change phase.     !
     !-------------------------------------------------------------------------------------!
     if (co2_on) call azero(mxp*myp*mzp,cuparm_g(ngrid)%co2src(:,:,:,icld))

     call azero(mxp*myp*mzp,cuparm_g(ngrid)%thsrc(:,:,:,icld))
     call azero(mxp*myp*mzp,cuparm_g(ngrid)%rtsrc(:,:,:,icld))
     call azero(mxp*myp,cuparm_g(ngrid)%conprr(:,:,icld))

     call cuparth(mynum,mgmxp,mgmyp,mgmzp,mzp,mxp,myp,ia,iz,ja,jz,i0,j0,maxiens         &
                 ,icld,iupmethod,depth_min(icld),cap_maxs,radius(icld)                  &
                 ,zkbmax,zcutdown,z_detr,DTLT,time                                      &
                 ,basic_g(ngrid)%UP              , basic_g(ngrid)%VP                    &
                 ,basic_g(ngrid)%WP              , sc1_grell_g(ngrid)%thetasta          &
                 ,basic_g(ngrid)%PP              , basic_g(ngrid)%PI0                   &
                 ,basic_g(ngrid)%DN0             , sc1_grell_g(ngrid)%rvsta             &
                 ,turb_g(ngrid)%kpbl             , turb_g(ngrid)%tkep  ,tkmin           &
                 ,micro_g(ngrid)%rcp             , grid_g(ngrid)%topt                   &
                 ,grid_g(ngrid)%RTGT             , tend%THT                             &
                 ,tend%RTT                       , tend%PT                              &
                 ,cuparm_g(ngrid)%THSRC(:,:,:,icld), cuparm_g(ngrid)%RTSRC(:,:,:,icld)  &
                 ,cuparm_g(ngrid)%CONPRR(:,:,icld) , sc1_grell_g(ngrid)%ierr4d          &
                 ,sc1_grell_g(ngrid)%jmin4d      , sc1_grell_g(ngrid)%kdet4d            &
                 ,sc1_grell_g(ngrid)%k224d       , sc1_grell_g(ngrid)%kbcon4d           &
                 ,sc1_grell_g(ngrid)%ktop4d      , sc1_grell_g(ngrid)%kpbl4d            &
                 ,sc1_grell_g(ngrid)%kstabi4d    , sc1_grell_g(ngrid)%kstabm4d          &
                 ,sc1_grell_g(ngrid)%xmb4d       , sc1_grell_g(ngrid)%edt4d             &
                 ,sc1_grell_g(ngrid)%zcup5d      , sc1_grell_g(ngrid)%pcup5d            &
                 ,sc1_grell_g(ngrid)%enup5d      , sc1_grell_g(ngrid)%endn5d            &
                 ,sc1_grell_g(ngrid)%deup5d      , sc1_grell_g(ngrid)%dedn5d            &
                 ,sc1_grell_g(ngrid)%zup5d       , sc1_grell_g(ngrid)%zdn5d             &
                 ,sc1_grell_g(ngrid)%prup5d      , sc1_grell_g(ngrid)%clwup5d           &
                 ,sc1_grell_g(ngrid)%tup5d                                              &
                 ,cuparm_g(ngrid)%upmf(:,:,icld)   , cuparm_g(ngrid)%dnmf(:,:,icld)     &
                 ,cuparm_g(ngrid)%xierr(:,:,icld)  , cuparm_g(ngrid)%zklnb(:,:,icld)    &
                 ,cuparm_g(ngrid)%zklfc(:,:,icld)  , cuparm_g(ngrid)%zklou(:,:,icld)    &
                 ,cuparm_g(ngrid)%zklod(:,:,icld)  , cuparm_g(ngrid)%zkdet(:,:,icld)    &
                 ,cuparm_g(ngrid)%xiact_p(:,:,icld), cuparm_g(ngrid)%xiact_c(:,:,icld)  )
! [MLO ------------- Stilt - RAMS coupling  ------------------
     if (imassflx == 1) then
        call old_prep_convflx_to_mass(mzp,mxp,myp,ia,iz,ja,jz,maxiens,ngrid             &
           ,sc1_grell_g(ngrid)%ierr4d      ,sc1_grell_g(ngrid)%jmin4d                   &
           ,sc1_grell_g(ngrid)%kdet4d      ,sc1_grell_g(ngrid)%k224d                    &
           ,sc1_grell_g(ngrid)%kbcon4d     ,sc1_grell_g(ngrid)%ktop4d                   &
           ,sc1_grell_g(ngrid)%kpbl4d      ,sc1_grell_g(ngrid)%kstabi4d                 &
           ,sc1_grell_g(ngrid)%kstabm4d    ,sc1_grell_g(ngrid)%xmb4d                    &
           ,sc1_grell_g(ngrid)%edt4d       ,sc1_grell_g(ngrid)%zcup5d                   &
           ,sc1_grell_g(ngrid)%pcup5d      ,sc1_grell_g(ngrid)%enup5d                   &
           ,sc1_grell_g(ngrid)%endn5d      ,sc1_grell_g(ngrid)%deup5d                   &
           ,sc1_grell_g(ngrid)%dedn5d      ,sc1_grell_g(ngrid)%zup5d                    &
           ,sc1_grell_g(ngrid)%zdn5d       ,icld                                        )
     end if ! (imassflx == 1)
! ------------- Stilt - RAMS coupling  ------------------ MLO]

!----- Shallow cumulus -------------------
  else if (icld == 2) then


     call azero(mxp*myp*mzp,cuparm_g(ngrid)%thsrc(:,:,:,icld))
     call azero(mxp*myp*mzp,cuparm_g(ngrid)%rtsrc(:,:,:,icld))

     !-------------------------------------------------------------------------------------!
     !     Zero out CO2 tendency if CO2 is prognosed.  The old Grell scheme won't compute  !
     ! the transport of CO2 through updrafts and downdrafts, feel free to add this.  It    !
     ! should be similar to the water transport, except that CO2 doesn't change phase.     !
     !-------------------------------------------------------------------------------------!
     if (co2_on) call azero(mxp*myp*mzp,cuparm_g(ngrid)%co2src(:,:,:,icld))

     call cuparth_shal(mynum,mgmxp,mgmyp,mgmzp,mzp,mxp,myp,ia,iz,ja,jz,i0,j0,maxiens,icld &
                      ,iupmethod,depth_min(icld),cap_maxs,radius(icld),zkbmax,zcutdown    &
                      ,z_detr,dtlt,time                                                   &
                      ,basic_g(ngrid)%up            ,basic_g(ngrid)%vp                    &
                      ,basic_g(ngrid)%wp            ,basic_g(ngrid)%theta                 &
                      ,basic_g(ngrid)%pp            ,basic_g(ngrid)%pi0                   &
                      ,basic_g(ngrid)%dn0           ,basic_g(ngrid)%rv                    &
                      ,turb_g(ngrid)%kpbl           ,turb_g(ngrid)%tkep,  tkmin           &
                      ,micro_g(ngrid)%rcp           ,grid_g(ngrid)%topt                   &
                      ,grid_g(ngrid)%rtgt           ,tend%tht                             &
                      ,tend%rtt                     ,tend%pt                              &
                      ,cuparm_g(ngrid)%thsrc(:,:,:,icld)                                  &
                      ,cuparm_g(ngrid)%rtsrc(:,:,:,icld)                                  &
                      ,sc1_grell_g(ngrid)%ierr4d    ,sc1_grell_g(ngrid)%jmin4d            &
                      ,sc1_grell_g(ngrid)%kdet4d    ,sc1_grell_g(ngrid)%k224d             &
                      ,sc1_grell_g(ngrid)%kbcon4d   ,sc1_grell_g(ngrid)%ktop4d            &
                      ,sc1_grell_g(ngrid)%kpbl4d    ,sc1_grell_g(ngrid)%kstabi4d          &
                      ,sc1_grell_g(ngrid)%kstabm4d  ,sc1_grell_g(ngrid)%xmb4d             &
                      ,sc1_grell_g(ngrid)%edt4d     ,sc1_grell_g(ngrid)%zcup5d            &
                      ,sc1_grell_g(ngrid)%pcup5d    ,sc1_grell_g(ngrid)%enup5d            &
                      ,sc1_grell_g(ngrid)%endn5d    ,sc1_grell_g(ngrid)%deup5d            &
                      ,sc1_grell_g(ngrid)%dedn5d    ,sc1_grell_g(ngrid)%zup5d             &
                      ,sc1_grell_g(ngrid)%zdn5d     ,sc1_grell_g(ngrid)%prup5d            &
                      ,sc1_grell_g(ngrid)%clwup5d   ,sc1_grell_g(ngrid)%tup5d             &
                      ,cuparm_g(ngrid)%upmf(:,:,icld)  ,cuparm_g(ngrid)%xierr(:,:,icld)   &
                      ,cuparm_g(ngrid)%zklnb(:,:,icld) ,cuparm_g(ngrid)%zklfc(:,:,icld)  &
                      ,cuparm_g(ngrid)%zklou(:,:,icld)  ,cuparm_g(ngrid)%xierr(:,:,1)       &
                      )
! [MLO ------------- Stilt - RAMS coupling  ------------------
     if (imassflx == 1) then
       call old_prep_convflx_to_mass(mzp,mxp,myp,ia,iz,ja,jz,maxiens,ngrid                &
           ,sc1_grell_g(ngrid)%ierr4d      ,sc1_grell_g(ngrid)%jmin4d                     &
           ,sc1_grell_g(ngrid)%kdet4d      ,sc1_grell_g(ngrid)%k224d                      &
           ,sc1_grell_g(ngrid)%kbcon4d     ,sc1_grell_g(ngrid)%ktop4d                     &
           ,sc1_grell_g(ngrid)%kpbl4d      ,sc1_grell_g(ngrid)%kstabi4d                   &
           ,sc1_grell_g(ngrid)%kstabm4d    ,sc1_grell_g(ngrid)%xmb4d                      &
           ,sc1_grell_g(ngrid)%edt4d       ,sc1_grell_g(ngrid)%zcup5d                     &
           ,sc1_grell_g(ngrid)%pcup5d      ,sc1_grell_g(ngrid)%enup5d                     &
           ,sc1_grell_g(ngrid)%endn5d      ,sc1_grell_g(ngrid)%deup5d                     &
           ,sc1_grell_g(ngrid)%dedn5d      ,sc1_grell_g(ngrid)%zup5d                      &
           ,sc1_grell_g(ngrid)%zdn5d       ,icld                                          )
     end if
! ------------- Stilt - RAMS coupling  ------------------ MLO]

  end if
  !MLO- Flushing scratch variables to zero, otherwise nested grids can get left-overs
  !     from the previous grid.
  call zero_scratch2_grell()
  call zero_scratch3_grell()
  call zero_scratch2_grell_sh() !LFR
  call zero_scratch3_grell_sh() !LFR
  return

end subroutine old_grell_cupar_driver
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine include_shal_effect(m1,m2,m3,ia,iz,ja,jz,dtlt                                   &
                              ,thetasta,rvsta,theta,rv,pi0,pp,thsrc,rtsrc)
   use therm_lib , only : rslf
   use rconstants, only : cpi,cpor,p00
   implicit none
   integer, intent(in)                                           :: ia,iz,ja,jz
   integer, intent(in)                                           :: m1,m2,m3
   real   , intent(in)                                           :: dtlt
   real   , intent(out) , dimension(m1,m2,m3)                    :: thetasta,rvsta
   real   , intent(in)  , dimension(m1,m2,m3)                    :: theta,rv,pi0,pp
   real   , intent(in)  , dimension(m1,m2,m3)                    :: thsrc,rtsrc
   integer                                                       :: i,j,k
   real                                                          :: press,rsat,tempk
   do j=ja,jz
      do i=ia,iz
         do k=2,m1
            ! Updating the potential temperature
            thetasta(k,i,j) = theta(k,i,j)+dtlt*thsrc(k,i,j)
            ! Finding the vapour mixing ratio after the shallow cumulus call
            press=p00*(cpi*(pi0(k,i,j)+pp(k,i,j)))**cpor
            tempk=cpi*theta(k,i,j)*(pi0(k,i,j)+pp(k,i,j))
            rsat =rslf(press,tempk)
            
            rvsta(k,i,j) = max(epsilon(1.),min(rsat,rv(k,i,j)+dtlt*rtsrc(k,i,j)))
         end do
      end do
   end do
   return
end subroutine include_shal_effect








!------------------------------------------------------------------------------------------!
! Subroutine prep_convflx_to_mass                                                          !
! Developed by Saulo R. Freitas (CPTEC/INPE)                                               !
!                                                                                          !
!   This subroutine prepares the convective fluxes to be used in STILT (or other           !
! Lagrangian models).                                                                      !
!------------------------------------------------------------------------------------------!
subroutine old_prep_convflx_to_mass(m1,m2,m3,ia,iz,ja,jz,maxiens,ifm                       &
                                ,ierr4d,jmin4d,kdet4d,k224d,kbcon4d,ktop4d                 &
                                ,kpbl4d,kstabi4d,kstabm4d,xmb4d,edt4d,zcup5d,pcup5d,enup5d &
                                ,endn5d,deup5d,dedn5d,zup5d,zdn5d,iens)

use mem_mass

implicit none

integer, intent(in) :: ifm,iens, maxiens,m1,m2,m3,ia,iz,ja,jz

integer, intent(in),dimension(m2,m3,maxiens) ::                                            &
                         ierr4d,jmin4d,kdet4d,k224d,kbcon4d,ktop4d,kpbl4d,kstabi4d,kstabm4d
                       
real, intent(in), dimension(m2,m3,maxiens) :: xmb4d,edt4d
                        
real, intent(in), dimension(m1,m2,m3,maxiens) ::                        &
                                      enup5d,endn5d,deup5d,dedn5d,zup5d,zdn5d,zcup5d,pcup5d

integer :: i, j


!zero out all convflx
  call azero(m1*m2*m3,mass_g(ifm)%cfxup(:,:,:,iens))
  call azero(m1*m2*m3,mass_g(ifm)%cfxdn(:,:,:,iens))
  call azero(m1*m2*m3,mass_g(ifm)%dfxup(:,:,:,iens))
  call azero(m1*m2*m3,mass_g(ifm)%efxup(:,:,:,iens))
  call azero(m1*m2*m3,mass_g(ifm)%dfxdn(:,:,:,iens))
  call azero(m1*m2*m3,mass_g(ifm)%efxdn(:,:,:,iens))

do j=ja,jz
  do i=ia,iz
!    if((iens == 1 .and. ierr4d(i,j,iens,ngrid) == 0) .or. iens == 2) then
    if(ierr4d(i,j,iens) == 0) then
      call get_convflx(iens,i,j,m1,m2,m3                                       &
                  ,   xmb4d(i,j,iens),   edt4d(i,j,iens)                       &
                  ,  jmin4d(i,j,iens),  kdet4d(i,j,iens)                       &
                  ,   k224d(i,j,iens), kbcon4d(i,j,iens)                       &
                  ,  ktop4d(i,j,iens),  kpbl4d(i,j,iens)                       &
                  ,kstabi4d(i,j,iens),kstabm4d(i,j,iens)                       &
                  ,zcup5d(1:m1,i,j,iens),pcup5d(1:m1,i,j,iens)                 &
                  ,deup5d(1:m1,i,j,iens),enup5d(1:m1,i,j,iens)                 &
                  ,dedn5d(1:m1,i,j,iens),endn5d(1:m1,i,j,iens)                 &
                  , zup5d(1:m1,i,j,iens), zdn5d(1:m1,i,j,iens)                 &
                  ,mass_g(ifm)%cfxup(:,:,:,iens),mass_g(ifm)%cfxdn(:,:,:,iens) &
                  ,mass_g(ifm)%dfxup(:,:,:,iens),mass_g(ifm)%efxup(:,:,:,iens) &
                  ,mass_g(ifm)%dfxdn(:,:,:,iens),mass_g(ifm)%efxdn(:,:,:,iens))
    endif 
  enddo 
enddo  

return
end subroutine old_prep_convflx_to_mass
!------------------------------------------------------------------------------------------!

!------------------------------------------------------------------------------------------!
! Subroutine get_convflx                                                                   !
! Developed by Saulo R. Freitas (CPTEC/INPE)                                               !
!                                                                                          !
!   This subroutine aims at getting the convective fluxes from Grell shallow and           !
! deep convective parameterizations.                                                       !
!------------------------------------------------------------------------------------------!
subroutine get_convflx(iens,i,j,m1,m2,m3,xmb,edt,jmin,kdet,k22,kbcon,ktop,kpbl             &
                      ,kstabi,kstabm,z_cup,p_cup,cd,entr,cdd,entrd,zu,zd,cfxup,cfxdn       &
                      ,dfxup,efxup,dfxdn,efxdn)

implicit none

integer, intent(in)                      ::  iens,i,j,m1,m2,m3,jmin,kdet,k22,kbcon   &
                                            ,ktop,kpbl,kstabi,kstabm
                                            
real, intent(in)                         ::  xmb,edt

real, intent(in),    dimension(m1)       ::  z_cup,p_cup,cd,entr,cdd,entrd,zu,zd

real, intent(inout), dimension(m1,m2,m3) ::  cfxup,cfxdn,dfxup,efxup,dfxdn,efxdn
                                             
integer                                  ::  k,kr
real                                     ::  dz,totmas,entup,detup,entdoj,entupk,detupk    &
                                            ,detdo,entdo,subdown,subin,detdo1,detdo2

do k=2,ktop+1
        
  kr= k + 1   ! level K of conv grid  corresponds to level K + 1 of RAMS grid
  dz =  z_cup(kr) - z_cup(k)
   
  entup   = 0.
  detup   = 0.
  entdoj  = 0.
  entupk  = 0.
  detupk  = 0.
  detdo   = edt*cdd(k)*dz*zd(kr)
  entdo   = edt*entrd(k)*dz*zd(kr)
  subdown = zu(k ) - edt*zd(k )
  subin   = zu(kr) - edt*zd(kr)

  if (k >= kbcon .and. k < ktop) then
    entup  = entr(k)   *dz*zu(k)
    detup  =   cd(kr) *dz*zu(k)
  end if

  if(k == jmin)  entdoj = edt*zd(k)
  if(k == k22-1) entupk = zu(kpbl)
  if(k == ktop)  detupk = zu(ktop)
  if(k > kdet)   detdo  = 0.
  if(k == ktop)  subin  = 0.
  if(k < kbcon)  detup  = 0.
  
  if(iens == 1) then ! Deep convection
      cfxup(k ,i,j) =     xmb* zu(k)
      cfxdn(k ,i,j) =-edt*xmb* zd(k)
      dfxup(kr,i,j) =     xmb*(detup + detupk)
      efxup(kr,i,j) =    -xmb*(entup + entupk)
      dfxdn(kr,i,j) =     xmb*(detdo         ) !edt already is at detdo
      efxdn(kr,i,j) =    -xmb*(entdo + entdoj) !edt already is at entdo,entdoj
  elseif(iens == 2)  then ! Shallow convection
      cfxup(k ,i,j) =     xmb* zu(k)
      dfxup(kr,i,j) =     xmb*(detup + detupk)
      efxup(kr,i,j) =    -xmb*(entup + entupk)
  end if
!------------------------------------------------------------------------------------------!
! Checking the mass conservation                                                           !
!------------------------------------------------------------------------------------------!
  totmas=subin-subdown+detup-entup-entdo+detdo-entupk-entdoj+detupk
  if(abs(totmas) > 1.e-6) then
    write (unit=*,fmt='(a)')                 '----------- Subroutine Get_convflx ----------'
    write(unit=*, fmt='(4(a,1x,i3,1x))')     '  K= ',k,'   I=',i,'   J=',j,'   IENS=',iens
    write(unit=*, fmt='(2(a,1x,es10.3,1x))') '  subin=  ',    subin,'subdown=',subdown
    write(unit=*, fmt='(2(a,1x,es10.3,1x))') '  detup=  ',    detup,'entup=  ',entup
    write(unit=*, fmt='(2(a,1x,es10.3,1x))') '  entdo=  ',    entdo,'detdo=  ',detdo
    write(unit=*, fmt='(2(a,1x,es10.3,1x))') '  entupk= ',   entupk,'detupk= ',detupk
    write(unit=*, fmt='(2(a,1x,es10.3,1x))') '  zu(k)=  ',    zu(k),'zd(k)=  ',zd(k)
    write(unit=*, fmt='(2(a,1x,es10.3,1x))') '  zu(kr)= ',   zu(kr),'zd(kr)= ',zd(kr)
    write(unit=*, fmt='(2(a,1x,es10.3,1x))') '  entdoj= ',   entdoj,'edt=    ',edt
    write(unit=*, fmt='(1(a,1x,es10.3,1x))') '  totmas= ',   totmas
    write(unit=*, fmt='(a)')                 '---------------------------------------------'
    call abort_run('The model will stop since it is not conserving mass...' &
                  ,'get_convfx','rmass.f90')
  end if

end do

!------------------------------------------------------------------------------------------!
! Bottom layer                                                                             !
!------------------------------------------------------------------------------------------!
if(iens == 1) then  ! Deep convection
  k = 1
  kr= k + 1  ! the K-level of Grell is equivalent to the BRAMS K+1-level

  dz        =  z_cup(2)-z_cup(1)

  detdo1    = edt*zd(2)*  cdd(1)*dz
  detdo2    = edt*zd(1)
  entdo     = edt*zd(2)*entrd(1)*dz
  subin     =-edt*zd(2)           

  cfxup(kr,i,j) = 0.
  cfxdn(k ,i,j) =-edt*xmb* zd(1)
  dfxup(kr,i,j) = 0.
  efxup(kr,i,j) = 0.
  dfxdn(kr,i,j) = xmb*(detdo1+detdo2) !edt already is at detdo1,2
  efxdn(kr,i,j) =-xmb* entdo          !edt already is at entdo
 

!------------------------------------------------------------------------------------------!
! Checking the mass conservation                                                           !
!------------------------------------------------------------------------------------------!
  totmas = detdo1+detdo2-entdo+subin
  if(abs(totmas) > 1.e-6) then
    write (unit=*,fmt='(a)')                 '----------- Subroutine Get_convflx ----------'
    write(unit=*, fmt='(4(a,1x,i3,1x))')     '  K= ',k,'   I=',i,'   J=',j,'   IENS=',iens
    write(unit=*, fmt='(2(a,1x,es10.3,1x))') '  subin=  ',    subin,'entdo=  ',entdo
    write(unit=*, fmt='(2(a,1x,es10.3,1x))') '  detdo1= ',   detdo1,'detdo2= ',detdo2
    write(unit=*, fmt='(1(a,1x,es10.3,1x))') '  totmas= ',   totmas
    write(unit=*, fmt='(a)')                 '---------------------------------------------'
    call abort_run('The model will stop since it is not conserving mass...' &
                  ,'get_convfx','rmass.f90')
  end if
end if

return
end subroutine get_convflx
!------------------------------------------------------------------------------------------!
