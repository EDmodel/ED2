!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

subroutine radiate(mzp,mxp,myp,ia,iz,ja,jz,mynum)

  use mem_tend   ,  only: tend
  use mem_grid   ,  only: ngrid,time,dtlt,itimea,nzg,nzs,npatch,grid_g,nnzp, &
                          if_adap,zm,zt
  use mem_leaf   ,  only: leaf_g
  use mem_radiate,  only: ilwrtyp,iswrtyp,icumfdbk,radiate_g,radfrq, &
                           ncall_i,prsnz,prsnzp,jday,solfac,nadd_rad
  use mem_basic  ,  only: basic_g
  use mem_scratch,  only: vctr1,vctr2,vctr3,vctr4,vctr5,vctr6,vctr7,vctr8, &
                          vctr9,vctr10,vctr11,vctr12,scratch
  use mem_micro  ,  only: micro_g
  use rconstants ,  only: cp,cpor,p00
  use mem_harr   ,  only: ng,nb,nsolb,npsb,nuum,prf,alpha,trf,beta,&
                          xp,wght,wlenlo,wlenhi,solar0,ralcs,a0,a1,a2,a3, &
		 	  exptabc,ulim,npartob,npartg,ncog,ncb,ocoef,bcoef,gcoef
  use ref_sounding, only: pi01dn

  ! CATT
  !kmlnew
  use therm_lib   , only: qtk, level, cloud_on, bulk_on
  use micphys     , only: gnu,icloud,irain,ipris,isnow,iaggr,igraup,ihail
  use mem_cuparm  , only: cuparm_g, cuparm_g_sh, nnqparm
  !kmlnew
  use rad_carma   , only: radcomp_carma
  use catt_start  , only: catt           ! intent(in)
  use mem_scalar  , only: scalar_g

  ! teb_spm
  use teb_spm_start, only: teb_spm !intent(in)
  use mem_teb_common, only: &
       tebc_g
       
  !MLO - For the new deal with parameterized clouds in Harrington
  use mem_grell,       only : grell_g
  use mem_shcu,        only : shcu_g
  use shcu_vars_const, only : nnshcu

  implicit none

  integer, intent(in) :: mzp,mxp,myp,ia,iz,ja,jz,mynum
  integer :: koff,ka,nrad

  real :: solc

  !kmlnew
  real, dimension(mzp,mxp,myp) :: lwl,iwl,fracao_liq
  real, dimension(mxp,myp)     :: rain
  real :: dummy
  integer :: ncall = 0
  integer :: i,j,k
  real :: max_albedt,max_rlongup
  !kmlnew
  real                 :: time_rfrq
  ! teb_spm
  real, pointer :: emis_town(:,:), alb_town(:,:), ts_town(:,:), g_urban(:,:,:)
  !

  !MLO - Pointers to ease the Harrington radiation call:
  real, pointer, dimension(:,:,:) :: ha_rcp,ha_rrp,ha_rpp,ha_rsp,ha_rap,ha_rgp,ha_rhp &
                                    ,ha_ccp,ha_crp,ha_cpp,ha_csp,ha_cap,ha_cgp,ha_chp
  real, pointer, dimension(:,:)   :: parm_rain,parm_xkbcon,parm_xjmin,parm_dnmf
  real, pointer, dimension(:,:,:) :: parm_rtsrc,parm_rtsrc_sh
  real, allocatable, target, dimension(:,:,:) :: dumzero3d
  real, allocatable, target, dimension(:,:)   :: dumzero2d

  logical :: grell_on

  if (ilwrtyp + iswrtyp == 0) return

  !MLO - Nullifying the pointers for Harrington/Bulk microphysics interface
  nullify(ha_rcp,ha_rrp,ha_rpp,ha_rsp,ha_rap,ha_rgp,ha_rhp &
         ,ha_ccp,ha_crp,ha_cpp,ha_csp,ha_cap,ha_cgp,ha_chp &
         ,parm_rain,parm_xkbcon,parm_xjmin,parm_dnmf       &
         ,parm_rtsrc,parm_rtsrc_sh,emis_town,alb_town,ts_town,g_urban)
  allocate(dumzero3d(mzp,mxp,myp),dumzero2d(mxp,myp))
  call azero(mzp*mxp*myp,dumzero3d)
  call azero(mxp*myp,dumzero2d)

  ! teb_spm
  if (teb_spm==1) then
     emis_town => tebc_g(ngrid)%emis_town
     alb_town  => tebc_g(ngrid)%alb_town
     ts_town   => tebc_g(ngrid)%ts_town
     g_urban   => leaf_g(ngrid)%g_urban
  else
     emis_town => dumzero2d
     alb_town  => dumzero2d
     ts_town   => dumzero2d
     g_urban   => dumzero3d
  endif
  !

  call tend_accum(mzp,mxp,myp,ia,iz,ja,jz,tend%tht(1),radiate_g(ngrid)%fthrd(1,1,1))
  
  if ((iswrtyp == 3 .or. ilwrtyp == 3) .and. ncall_i == 0) then
     ! if first call for this node, initialize several quantities & mclatchy
     ! sounding data.
     write (unit=*,fmt=*) '----> Initializing Harrington on node ',mynum
     call harr_radinit(nnzp(1))
     ncall_i = ncall_i + 1
  end if


  time_rfrq = real(dmod(time + 0.001,dble(radfrq)))

  if ( time_rfrq  < dtlt .or. time < 0.001) then
     ! compute solar zenith angle, multiplier for solar constant, sfc albedo,
     ! and surface upward longwave radiation.
     call radprep(mxp,myp,nzg,nzs,npatch,ia,iz,ja,jz   &
          ,iswrtyp,ilwrtyp                             &
          ,leaf_g(ngrid)%soil_water      (:,:,:,:)     &
          ,leaf_g(ngrid)%soil_energy     (:,:,:,:)     &
          ,leaf_g(ngrid)%soil_text       (:,:,:,:)     &
          ,leaf_g(ngrid)%sfcwater_energy (:,:,:,:)     &
          ,leaf_g(ngrid)%sfcwater_depth  (:,:,:,:)     &
          ,leaf_g(ngrid)%leaf_class      (:,:,:)       &
          ,leaf_g(ngrid)%veg_fracarea    (:,:,:)       &
          ,leaf_g(ngrid)%veg_height      (:,:,:)       &
          ,leaf_g(ngrid)%veg_albedo      (:,:,:)       &
          ,leaf_g(ngrid)%patch_area      (:,:,:)       &
          ,leaf_g(ngrid)%sfcwater_nlev   (:,:,:)       &
          ,leaf_g(ngrid)%veg_temp        (:,:,:)       &
          ,leaf_g(ngrid)%can_temp        (:,:,:)       &
          ! teb_spm                                   
          ,emis_town                                   &
          ,alb_town                                    &
          ,ts_town                                     &
          ,g_urban                                     &
          !                                           
          ,grid_g(ngrid)%glat             (:,:)        &
          ,grid_g(ngrid)%glon             (:,:)        &
          ,radiate_g(ngrid)%rlongup       (:,:)        &
          ,radiate_g(ngrid)%rlong_albedo  (:,:)        &
          ,radiate_g(ngrid)%albedt        (:,:)        &
          ,radiate_g(ngrid)%rshort        (:,:)        &
          ,radiate_g(ngrid)%rlong         (:,:)        &
          ,radiate_g(ngrid)%rshort_top    (:,:)        &
          ,radiate_g(ngrid)%rshortup_top  (:,:)        &
          ,radiate_g(ngrid)%cosz          (:,:)        )
                                                      
     call azero(mzp*mxp*myp,radiate_g(ngrid)%fthrd(1,1,1))
     call azero(mzp*mxp*myp,radiate_g(ngrid)%fthrd_lw(1,1,1))

     ! CARMA radiation
     if (ilwrtyp==4 .or. iswrtyp==4) then

        if (level == 1) then

           ! not used with level 1 - putting zeros
           call azero(mzp*mxp*myp,lwl)
           call azero(mzp*mxp*myp,iwl)
           call azero(mxp*myp,rain)

        elseif (level == 2) then
           call atob (mxp * myp * mzp,micro_g(ngrid)%rcp(:,:,:),lwl)

           call azero(mzp*mxp*myp,iwl) ! not used with level 2 - putting zeros

           if (nnqparm(ngrid) > 0) then
              do i=ia,iz
                 do j=ja,jz
                    ! conv  precip at mm/h
                    rain(i,j)= cuparm_g(ngrid)%conprr(i,j)
                 enddo
              enddo
           endif

        elseif (bulk_on) then

           if (nnqparm(ngrid) > 0) then
              do i=ia,iz
                 do j=ja,jz
                    ! conv + explic  precip at mm/h
                    rain(i,j)= (cuparm_g(ngrid)%conprr(i,j) + &
                         micro_g(ngrid)%pcpg(i,j)) * 3600.
                 enddo
              enddo
           endif
           
           call azero2(mzp*mxp*myp,lwl,iwl)
           if(icloud>0) &
                call accum (mxp * myp * mzp,lwl,micro_g(ngrid)%rcp(:,:,:))

           if(igraup>0) then
              do k=1,mzp
                 
                 do i=ia,iz
                    do j=ja,jz
                       call qtk(micro_g(ngrid)%q6(k,i,j), &
                            dummy,fracao_liq(k,i,j))
                    enddo
                 enddo
              enddo
              call ae1t1p1(mxp * myp * mzp,lwl,fracao_liq,  &
                   micro_g(ngrid)%rgp(:,:,:),lwl)

              fracao_liq(:,:,:) = 1. - fracao_liq(:,:,:)

              call ae1t1p1(mxp * myp * mzp,iwl,fracao_liq,&
                      micro_g(ngrid)%rgp(:,:,:),iwl)

           endif

           if(ihail>0) then
              do k=1,mzp

                 do i=ia,iz
                    do j=ja,jz
                       call qtk(micro_g(ngrid)%q7(k,i,j),dummy,  &
                            fracao_liq(k,i,j))
                    enddo
                 enddo
              enddo
              call ae1t1p1(mxp*myp*mzp,lwl,fracao_liq, &
                   micro_g(ngrid)%rhp(:,:,:),lwl)
              
              fracao_liq(:,:,:) = 1. - fracao_liq(:,:,:)
              
              call ae1t1p1(mxp*myp*mzp,iwl,  &
                   fracao_liq,micro_g(ngrid)%rhp(:,:,:),iwl)

           endif

           if(iaggr>0) call accum (mxp*myp*mzp,iwl,micro_g(ngrid)%rap(:,:,:))
           if(isnow>0) call accum (mxp*myp*mzp,iwl,micro_g(ngrid)%rsp(:,:,:))
           if(ipris>0) call accum (mxp*myp*mzp,iwl,micro_g(ngrid)%rpp(:,:,:))

        endif
        !kmln

        call radcomp_carma(mzp,mxp,myp,ia,iz,ja,jz,solfac  &
             ,basic_g(ngrid)%theta       &
             ,basic_g(ngrid)%pi0         &
             ,basic_g(ngrid)%pp          &
             ,basic_g(ngrid)%rv          &
             !kmlnew
             ,rain,lwl,iwl               &
             !kmlnew
             ,basic_g(ngrid)%dn0         &
             ,basic_g(ngrid)%rtp         &
             ,radiate_g(ngrid)%fthrd     &
             ,grid_g(ngrid)%rtgt         &
             ,grid_g(ngrid)%f13t         &
             ,grid_g(ngrid)%f23t         &
             ,grid_g(ngrid)%glat         &
             ,grid_g(ngrid)%glon         &
             ,radiate_g(ngrid)%rshort    &
             ,radiate_g(ngrid)%rlong     &
             ,radiate_g(ngrid)%albedt    &
             ,radiate_g(ngrid)%cosz      &
             ,radiate_g(ngrid)%rlongup   &
             ,mynum                      &
             !srf - carma arrays
             ,grid_g(ngrid)%fmapt        &
             ,scalar_g(3,ngrid)%sclp     &
             !kmlnew
             ,leaf_g(ngrid)%patch_area   &
             ,npatch                     &
             !lmlnew
             )
     end if

     if (ilwrtyp <= 2 .or. iswrtyp <= 2) then

        ! if using mahrer-pielke and/or chen-cotton radiation, call radcomp.

        call radcomp(mzp,mxp,myp,ngrid,ia,iz,ja,jz    &
             ,basic_g(ngrid)%theta      (1,1,1)  &
             ,basic_g(ngrid)%pi0        (1,1,1)  &
             ,basic_g(ngrid)%pp         (1,1,1)  &
             ,basic_g(ngrid)%rv         (1,1,1)  &
             ,basic_g(ngrid)%dn0        (1,1,1)  &
             ,basic_g(ngrid)%rtp        (1,1,1)  &
             ,radiate_g(ngrid)%fthrd    (1,1,1)  &
             ,grid_g(ngrid)%rtgt        (1,1)    &
             ,grid_g(ngrid)%f13t        (1,1)    &
             ,grid_g(ngrid)%f23t        (1,1)    &
             ,grid_g(ngrid)%glon        (1,1)    &
             ,radiate_g(ngrid)%rshort   (1,1)    &
             ,radiate_g(ngrid)%rlong    (1,1)    &
             ,radiate_g(ngrid)%albedt   (1,1)    &
             ,radiate_g(ngrid)%cosz     (1,1)    &
             ,radiate_g(ngrid)%rlongup  (1,1)    &
             ,radiate_g(ngrid)%fthrd_lw (1,1,1)  &
             ,mynum)

     end if

     ! using Harrington radiation
     if (iswrtyp == 3 .or. ilwrtyp == 3) then

        ![MLO - Associating pointers to the arrays or to the dummy array depending on the 
        !       microphysics and convection set up. First I assume no condensed hydrometeors
        !       and no cumulus parameterization, so everything is pointing to the dummy arrays
        ha_rcp        => dumzero3d
        ha_rrp        => dumzero3d
        ha_rpp        => dumzero3d
        ha_rsp        => dumzero3d
        ha_rap        => dumzero3d
        ha_rgp        => dumzero3d
        ha_rhp        => dumzero3d
        
        ha_ccp        => dumzero3d
        ha_crp        => dumzero3d
        ha_cpp        => dumzero3d
        ha_csp        => dumzero3d
        ha_cap        => dumzero3d
        ha_cgp        => dumzero3d
        ha_chp        => dumzero3d

        parm_rain     => dumzero2d
        parm_xkbcon   => dumzero2d
        parm_xjmin    => dumzero2d
        parm_dnmf     => dumzero2d
        
        parm_rtsrc    => dumzero3d
        parm_rtsrc_sh => dumzero3d

        ! Now I check whether condensation is allowed. If so, I repoint to the actual arrays.
        if (cloud_on) then
           ha_rcp => micro_g(ngrid)%rcp
           if (icloud == 5) ha_ccp => micro_g(ngrid)%ccp
        end if
        ! Now I check whether the bulk microphysics is activated. If so, I repoint to the actual arrays.
        if (level == 3) then
           ha_rrp => micro_g(ngrid)%rrp
           ha_rpp => micro_g(ngrid)%rpp
           ha_rsp => micro_g(ngrid)%rsp
           ha_rap => micro_g(ngrid)%rap
           ha_rgp => micro_g(ngrid)%rgp
           ha_rhp => micro_g(ngrid)%rhp
           if (irain  == 5) ha_crp => micro_g(ngrid)%crp
           if (ipris  == 5) ha_cpp => micro_g(ngrid)%cpp
           if (isnow  == 5) ha_csp => micro_g(ngrid)%csp
           if (iaggr  == 5) ha_cap => micro_g(ngrid)%cap
           if (igraup == 5) ha_cgp => micro_g(ngrid)%cgp
           if (ihail  == 5) ha_chp => micro_g(ngrid)%chp
        end if
! MLO - Now I'm going to add the effects of parameterized rain in the radiation.
!       The way I'm going to do this is by drawing a very stupid cloud, which
!       contains rain drops between the cloud base and level of origin of downdrafts, and 
!       cloud droplets wherever rtsrc is greater than zero. Yes, you are right, these clouds 
!       don't even have ice particles, but, come on, they didn't even exist until the last revision...
!       Here is the shopping section. The values were already set up as zeroes, so they will only be 
!       reset if the feature is available.
! 1. Getting the precipitation rate (which will be magically transformed into rain drops);
! 2. Getting the moistening rate (which will be magically transformed into cloud droplets);
        if (nnqparm(ngrid) > 0) then
          parm_rain  => cuparm_g(ngrid)%conprr
          parm_rtsrc => cuparm_g(ngrid)%rtsrc
        end if
        if (nnshcu(ngrid) == 1) then
          parm_rtsrc_sh => shcu_g(ngrid)%rtsrcsh
        elseif (nnshcu(ngrid) == 2) then
          parm_rtsrc_sh => cuparm_g_sh(ngrid)%rtsrc
        end if
! 3. If I am using Grell scheme, I also provide the cloud base, level of origin of downdrafts and
!    the downdraft mass flux (just to find a typical value for convective downdrafts)
        grell_on = nnqparm(ngrid) == 2
        if (grell_on) then
          parm_xkbcon => grell_g(ngrid)%xkbcon
          parm_xjmin  => grell_g(ngrid)%xjmin
          parm_dnmf   => grell_g(ngrid)%dnmf
        end if
        call harr_raddriv(mzp,mxp,myp,ngrid,if_adap,ia,iz,ja,jz,nadd_rad,iswrtyp,ilwrtyp,icumfdbk    &
             ,grid_g(ngrid)%flpw            (1,1  )  ,grid_g(ngrid)%topt            (1,1)    &
             ,grid_g(ngrid)%glat            (1,1  )  ,grid_g(ngrid)%rtgt            (1,1  )  &
             ,basic_g(ngrid)%pi0            (1,1,1)  ,basic_g(ngrid)%pp             (1,1,1)  &
             ,basic_g(ngrid)%dn0            (1,1,1)  ,basic_g(ngrid)%theta          (1,1,1)  &
             ,basic_g(ngrid)%rv             (1,1,1)  ,radiate_g(ngrid)%rshort       (1,1  )  &
             ,radiate_g(ngrid)%rlong        (1,1  )  ,radiate_g(ngrid)%fthrd        (1,1,1)  &
             ,radiate_g(ngrid)%rlongup      (1,1  )  ,radiate_g(ngrid)%cosz         (1,1  )  &
             ,radiate_g(ngrid)%albedt       (1,1  )  ,radiate_g(ngrid)%rshort_top   (1,1  )  &
             ,radiate_g(ngrid)%rshortup_top (1,1  )  ,radiate_g(ngrid)%rlongup_top  (1,1  )  &
             ,radiate_g(ngrid)%fthrd_lw     (1,1,1)  ,ha_rcp                        (1,1,1)  &
             ,ha_rrp                        (1,1,1)  ,ha_rpp                        (1,1,1)  &
             ,ha_rsp                        (1,1,1)  ,ha_rap                        (1,1,1)  &
             ,ha_rgp                        (1,1,1)  ,ha_rhp                        (1,1,1)  &
             ,ha_ccp                        (1,1,1)  ,ha_crp                        (1,1,1)  &
             ,ha_cpp                        (1,1,1)  ,ha_csp                        (1,1,1)  &
             ,ha_cap                        (1,1,1)  ,ha_cgp                        (1,1,1)  &
             ,ha_chp                        (1,1,1)  ,parm_rain                     (1,1)    &
             ,parm_xkbcon                   (1,1)    ,parm_xjmin                    (1,1)    &
             ,parm_dnmf                     (1,1)    ,parm_rtsrc                    (1,1,1)  &
             ,parm_rtsrc_sh                 (1,1,1)  &
             ,grell_on                               ,mynum                                  )
     end if
  end if
  deallocate(dumzero2d,dumzero3d)
  return
end subroutine radiate

!*****************************************************************************

subroutine tend_accum(m1,m2,m3,ia,iz,ja,jz,at,at2)

  implicit none
  integer :: m1,m2,m3,ia,iz,ja,jz,i,j,k
  real, dimension(m1,m2,m3) :: at,at2

  do j = ja,jz
     do i = ia,iz
        do k = 1,m1
           at(k,i,j) = at(k,i,j) + at2(k,i,j)
        end do
     end do
  end do

  return
end subroutine tend_accum

!*****************************************************************************

subroutine radprep(m2,m3,mzg,mzs,np,ia,iz,ja,jz,iswrtyp,ilwrtyp,             &
                   soil_water,soil_energy, &
                   soil_text,sfcwater_energy,sfcwater_depth,leaf_class,      &
                   veg_fracarea,veg_height,veg_albedo,patch_area,            &
                   sfcwater_nlev,veg_temp,can_temp,                          &
                   ! TEB_SPM
                   emis_town,alb_town,ts_town,g_urban,                       &
                   !
                   glat,glon,rlongup,rlong_albedo,albedt,rshort,             &
                   rlong,rshort_top,rshortup_top,cosz                        )

  ! TEB_SPM
  use teb_spm_start, only: TEB_SPM !INTENT(IN)
  use mem_leaf, only: isfcl   !RGK

  implicit none

  ! arguments:
  integer                       :: m2,m3,mzg,mzs,np,ia,iz,ja,jz,iswrtyp,ilwrtyp
  !dimension(mzg,m2,m3,np)
  real, dimension(mzg,m2,m3,np) :: soil_water,soil_energy,soil_text
  !dimension(mzs,m2,m3,np)
  real, dimension(mzs,m2,m3,np) :: sfcwater_energy,sfcwater_depth
  !dimension(m2,m3,np)
  real, dimension(m2,m3,np)     :: leaf_class,veg_fracarea,veg_height, &
                                   veg_albedo,patch_area,sfcwater_nlev,veg_temp,can_temp
  !TEB_SPM
  !dimension(m2,m3)
  real, dimension(m2,m3) :: emis_town,alb_town,ts_town
  !dimension(m2,m3,np)
  real, dimension(m2,m3,np) :: g_urban
  !dimension(m2,m3)
  real, dimension(m2,m3) :: glat,glon,rlongup,rlong_albedo,albedt,rshort,     &
                   rlong,rshort_top,rshortup_top,cosz


  INTEGER :: ip,i,j
  ! REAL :: c1,c2

  ! Interface necessary to use pointer as argument - TEB_SPM

  ! Compute solar zenith angle [cosz(i,j)] & solar constant factr [solfac].

  call zen(m2,m3,ia,iz,ja,jz,iswrtyp,ilwrtyp,glon,glat,cosz)

  ! Compute patch-averaged surface albeDO [albedt(i,j)] and up longwave
  ! radiative flux [rlongup(i,j)].

  if (isfcl /= 5) then   !RGK
     ! Zero out rlongup, albedt, rshort, rlong, and fthrd prior to summing over 
     ! land/sea flux cells

     call azero2(m2*m3,rlongup,albedt)
     do ip = 1,np
        do j = 1,jz
           do i = 1,iz
                            
              call sfcrad(mzg,mzs,ip               &
                   ,soil_energy    (1:mzg,i,j,ip) ,soil_water      (1:mzg,i,j,ip)  &
                   ,soil_text      (1:mzg,i,j,ip) ,sfcwater_energy (1:mzs,i,j,ip)  &
                   ,sfcwater_depth (1:mzs,i,j,ip) ,patch_area      (i,j,ip)    &
                   ,can_temp       (i,j,ip)   ,veg_temp        (i,j,ip)    &
                   ,leaf_class     (i,j,ip)   ,veg_height      (i,j,ip)    &
                   ,veg_fracarea   (i,j,ip)   ,veg_albeDO      (i,j,ip)    &
                   ,sfcwater_nlev  (i,j,ip)                                &
                   ,rshort         (i,j)      ,rlong           (i,j)       &
                   ,albedt         (i,j)      ,rlongup         (i,j)       &
                   ,cosz           (i,j)                                   &
                   ! TEB_SPM
                   ,g_urban        (i,j,ip)   ,emis_town       (i,j)       &
                   ,alb_town       (i,j)      ,ts_town         (i,j)       &
                   !
                   )
           end do
        end do
     end do
  endif   !RGK
  return
end subroutine radprep

!*****************************************************************************

subroutine radcomp(m1,m2,m3,ifm,ia,iz,ja,jz  &
     ,theta,pi0,pp,rv,dn0,rtp,fthrd  &
     ,rtgt,f13t,f23t,glon,rshort,rlong,albedt,cosz,rlongup,fthrd_lw  &
     ,mynum)

  use mem_grid   , only : dzm,dzt,itopo,plonn,ngrid,time,itimea,centlon
  use mem_scratch, only : scratch
  use mem_radiate, only : ilwrtyp,iswrtyp,lonrad,cdec,jday,solfac,sun_longitude
  use rconstants , only : cp,cpor,p00,stefan,solar,pio180,pi1,halfpi

  implicit none
  !----- Arguments --------------------------------------------------------------------------!
  integer , intent(in)                        :: m1,m2,m3,ifm,ia,iz,ja,jz,mynum
  real    , intent(in)  , dimension(m1,m2,m3) :: theta,pi0,pp,rv,dn0,rtp
  real    , intent(in)  , dimension(m2,m3)    :: rtgt,f13t,f23t,glon,cosz,albedt,rlongup
  real    , intent(out) , dimension(m1,m2,m3) :: fthrd, fthrd_lw
  real    , intent(out) , dimension(m2,m3)    :: rshort,rlong
  !----- Local variables --------------------------------------------------------------------!
  integer                                     :: i,j,k,kk
  real(kind=8)                                :: dzsdx,dzsdy,dlon,a1,a2,hrangl,sinz
  real(kind=8)                                :: sazmut,slazim,slangl,cosi,gglon 
  real                  , dimension(m1)       :: rvr,rtr,dn0r,pird,prd,dzmr,dztr,temprd
  real                  , dimension(m1)       :: fthrl,fthrs
  !----- Constants --------------------------------------------------------------------------!
  real(kind=8), parameter :: offset=1.d-20

  ! MLO - 2007-11-03. Used the new variables computed on zen subroutine to find the hour angle, and 
  !                   switched most sin/cos operations to double precision.

  !Constant longitude if lonrad is set to 0
  gglon=dble(centlon(1))

  do j = ja,jz
     do i = ia,iz

        do k = 1,m1
           pird(k) = (pp(k,i,j) + pi0(k,i,j)) / cp
           temprd(k) = theta(k,i,j) * pird(k)
           rvr(k) = max(0.,rv(k,i,j))
           rtr(k) = max(rvr(k),rtp(k,i,j))
           ! convert the next 7 variables to cgs for now.
           prd(k) = pird(k) ** cpor * p00 * 10.
           dn0r(k) = dn0(k,i,j) * 1.e-3
           dzmr(k) = dzm(k) / rtgt(i,j) * 1.e-2
           dztr(k) = dzt(k) / rtgt(i,j) * 1.e-2
        end do
        temprd(1) = (rlongup(i,j) / stefan) ** 0.25

        do k=1,m1
           if (prd(k) <  0. .or. dn0r(k) <   0. .or.  &
               rtr(k) <  0. .or. temprd(k) < 160.) then   
               ! tl(k) < 160.: This is -113 C, which is much colder than the 
               ! Vostok, Antarctica world record and should also
               ! be colder than any atmospheric temperature

              print*, 'Temperature too low or negative value of'
              print*, 'density, vapor, or pressure!'
              print*, 'before calling Chen-Cotton/Mahrer-Pielke radiation'
              print*, 'at grid: ',ifm,' at (k,i,j) = ',k,i,j
              print*, 'stopping model'
              print*, 'rad: k, rl(k), dl(k), pl(k), o3l(k), tl(k)'
              do kk=1,m1
                print'(i3,1x,5(es12.3,1x))', kk, rtr(kk), dn0r(kk), prd(kk), temprd(kk)
              enddo
              stop '»»»» STOP! radcomp (rad_driv.f90)'
           endif
        end do

        ! call the longwave parameterizations.

        if (ilwrtyp == 2) then     !>-- Mahrer-Pielke --<!
           call lwradp(m1,temprd,rvr,dn0r,dztr,pird,scratch%vt3dq(1),fthrl,rlong(i,j))
        elseif (ilwrtyp == 1) then !>-- Chen-Cotton --<!
!           call lwradc(m1,rvr,rtr,dn0r,temprd,prd,dztr,fthrl,rlong(i,j))
           call lwradc(m1,rvr,rvr,dn0r,temprd,prd,dztr,fthrl,rlong(i,j))
        end if

        ! the shortwave parameterizations are only valid if the cosine
        !    of the zenith angle is greater than .03 .

        if (cosz(i,j) > .03) then

           if (iswrtyp == 2) then     !>-- Mahrer-Pielke --<!
              call shradp(m1,rvr,dn0r,dzmr,scratch%vt3dq(1),pird,cosz(i,j)  &
                   ,albedt(i,j),solar*1e3*solfac,fthrs,rshort(i,j))
           elseif (iswrtyp == 1) then !>-- Chen-Cotton --<!
              call shradc(m1,rvr,rtr,dn0r,dztr,prd  &
                   ,albedt(i,j),solar*1.e3*solfac,cosz(i,j),fthrs,rshort(i,j))
           end if

           ! modify the downward surface shortwave flux by considering
           !    the slope of the topography.

           if (itopo == 1) then
              dzsdx = dble(f13t(i,j) * rtgt(i,j))
              dzsdy = dble(f23t(i,j) * rtgt(i,j))

              ! the y- and x-directions must be true north and east for
              ! this correction. the following rotates the model y/x
              ! to the true north/east.

              ! the following rotation seems to be incorrect, so call this instead:

              dlon = (dble(plonn(ngrid)) - dble(glon(i,j))) * pio180
              a1 = dzsdx*dcos(dlon) + dzsdy * dsin(dlon)
              a2 = -dzsdx*dsin(dlon) + dzsdy * dcos(dlon)
              dzsdx = a1
              dzsdy = a2

              if (lonrad == 1) gglon = dble(glon(i,j))
              hrangl = (sun_longitude-gglon)*pio180
              sinz = sqrt(dble(1.) - dble(cosz(i,j)) ** 2)
              sazmut = dasin(max(dble(-1.),min(dble(1.),cdec*dsin(hrangl)/dble(sinz))))
              
              ! Offset but preserving the sign...
              if (abs(dzsdx) < offset) dzsdx = sign(offset,dzsdx)
              if (abs(dzsdy) < offset) dzsdy = sign(offset,dzsdy)

              slazim = halfpi - datan2(dzsdy,dzsdx)
              slangl = datan(sqrt(dzsdx*dzsdx+dzsdy*dzsdy))
              cosi = dcos(slangl) * dble(cosz(i,j)) + dsin(slangl) * sinz  &
                   * dcos(sazmut-slazim)
              rshort(i,j) = rshort(i,j) * sngl(cosi) / cosz(i,j)
           end if

        else
           do k = 1,m1
              fthrs(k) = 0.
           end do
           rshort(i,j) = 0.
        end if

        do k = 2,m1-1
           fthrd(k,i,j) = fthrl(k) + fthrs(k)
           fthrd_lw(k,i,j) = fthrl(k)
        end do

        ! convert the downward flux at the ground to SI.

        rshort(i,j)     = rshort(i,j) * 1.e-3 / (1. - albedt(i,j))
        rlong(i,j)      = rlong(i,j) * 1.e-3
        fthrd(1,i,j)    = fthrd(2,i,j)
        fthrd_lw(1,i,j) = fthrd_lw(2,i,j)
     end do
  end do
  return
end subroutine radcomp

!******************************************************************************

subroutine zen(m2,m3,ia,iz,ja,jz,iswrtyp,ilwrtyp,glon,glat,cosz)

  use mem_grid   , only: nzpmax,imontha,idatea,iyeara,time,itimea,centlat, &
                         centlon
  use mem_radiate, only: lonrad, jday, solfac, sdec, cdec, declin, sun_longitude
  use rconstants , only: pio180,twopi,day_sec
  use mem_mclat,   only: mclat_spline
  use mem_harr,    only: nsolb, solar0,solar1
  
  implicit none

  !----- Arguments -----------------------------------------------------------------------!
  integer, intent(in)                    :: m2,m3           ! Grid dimensions
  integer, intent(in)                    :: ia,iz,ja,jz     ! Node dimensions
  integer, intent(in)                    :: iswrtyp,ilwrtyp ! Radiation scheme flag
  real   , intent(in) , dimension(m2,m3) :: glon,glat       ! Grid coordinates
  real   , intent(out), dimension(m2,m3) :: cosz
  !----- Local variables -----------------------------------------------------------------!
  integer             :: i,j            ! Grid counters
  integer             :: is             ! counter over solar bands in Harrington radiation
  integer, external   :: julday         ! Function to compute Julian day
  real(kind=8)        :: t1             ! 2 pi times fraction of year elapsed
  real(kind=8)        :: t2             ! 2 pi times fraction of year elapsed with offset
  real(kind=8)        :: eqn_of_time    ! equation of time solution [s]
  real(kind=8)        :: d0             ! coefficient for solfac computation
  real(kind=8)        :: d02            ! coefficient for solfac computation
  real(kind=8)        :: utc_sec        ! seconds elapsed in current simulation day (UTC)
  real(kind=8)        :: radlat         ! latitude in radians
  real(kind=8)        :: hrangl         ! Hour angle in radians

  ! Find current Julian day
  jday = julday(imontha,idatea,iyeara)
  jday = jday + floor(time/day_sec)

  ! Solfac is a multiplier of the solar constant to correct for Earth's
  ! varying distance to the sun.
  d0 =  twopi * dble(jday-1) / 365.
  d02 = d0 * 2.
  solfac = sngl(1.000110 + 0.034221 * dcos (d0) + 0.001280 * dsin(d0)  &
                         + 0.000719 * dcos(d02) + 0.000077 * dsin(d02))

  ! Check whether Harrington shortwave or longwave radiation is used
  if (iswrtyp == 3 .or. ilwrtyp == 3) then
   ! Adjust solar fluxes at top of atmosphere for current Earth-Sun distance      
   ! for Harrington shortwave radiation                                           
    do is = 1,nsolb                                                             
       solar1(is) = solar0(is) * solfac                                         
    enddo                                                                       
    ! Interpolate Mclatchy soundings between summer and winter values, and prepare 
    ! spline coefficients for interpolation by latitude.                           
    call mclat_spline(jday)                                                     
  end if


  ! Declin is the solar latitude in degrees   
  t1 = twopi * dble(jday) / 366.
  declin =  .322003                                                               &
           - 22.971  * dcos(t1) - .357898 * dcos(t1 * 2.) - .14398  * dcos(t1 * 3.)  &
           + 3.94638 * dsin(t1) + .019334 * dsin(t1 * 2.) + .05928  * dsin(t1 * 3.)
  t2 = (279.134 + .985647 * dble(jday)) * pio180
  !      sdec - sine of declination, cdec - cosine of declination
  sdec = dsin(declin*pio180)
  cdec = dcos(declin*pio180)

! The equation of time gives the number of seconds by which sundial time
! leads clock time
  eqn_of_time = 5.0323                  &
              - 100.976 * dsin(t2)       - 430.847 * dcos(t2)       &
              + 595.275 * dsin(t2 * 2.)  + 12.5024 * dcos(t2 * 2.)  &
              + 3.6858  * dsin(t2 * 3.)  + 18.25   * dcos(t2 * 3.)  &
              - 12.47   * dsin(t2 * 4.)

  ! Find the UTC time in seconds
  utc_sec = dmod(time+3600.*dble(itimea/100)+60.*dble(mod(itimea,100)),dble(day_sec))

  ! Find the longitude where the sun is at zenith
  sun_longitude = 180. - 360. * (utc_sec + eqn_of_time) / day_sec
  
  ! Compute the cosine of zenith angle
  if (lonrad == 0) then
    hrangl=(sun_longitude-dble(centlon(1)))*pio180
    do j=ja,jz
      do i=ia,iz
        radlat=dble(glat(i,j))*pio180
        cosz(i,j) = sngl(dsin(radlat)*sdec+dcos(radlat)*cdec*dcos(hrangl))
        ! Making sure that it is bounded
        cosz(i,j) = max(-1.,min(1.,cosz(i,j)))
      end do
    end do
  else
    do j=ja,jz
      do i=ia,iz
        radlat=dble(glat(i,j))*pio180
        hrangl=(sun_longitude-dble(glon(i,j)))*pio180
        cosz(i,j) = sngl(dsin(radlat)*sdec+dcos(radlat)*cdec*dcos(hrangl))
        ! Making sure that it is bounded
        cosz(i,j) = max(-1.,min(1.,cosz(i,j)))
      end do
    end do
  end if

 return
end subroutine zen

!****************************************************************************
