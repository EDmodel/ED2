! Adapted in 17/07/2002 by Alvaro L.Fazenda for version 5.0
! Adapted in 15/04/2003 by Alvaro L.Fazenda for version 5.04
! Adapted in 06/03/2007 by Alvaro L.Fazenda for BRAMS V.4.0
!-------------------------

subroutine cuparm_grell(iens)

  use mem_basic         , only: basic_g

  use mem_tend          , only: tend
![MLO - Adding shallow cumulus output
  use mem_cuparm        , only: nnqparm, confrq, cuparm_g,cuparm_g_sh
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
       ngrids_cp, & ! INTENT(IN)
       icbase,    & ! INTENT(IN)
       depth_min, & ! INTENT(IN)
       cap_maxs   ! ! INTENT(IN)

  use mem_scratch1_grell    , only: sc1_grell_g
  use mem_scratch2_grell    , only: zero_scratch2_grell
  use mem_scratch2_grell_sh , only: zero_scratch2_grell_sh
  use mem_scratch3_grell    , only: zero_scratch3_grell
  use mem_scratch3_grell_sh , only: zero_scratch3_grell_sh

  use mem_grell         , only: grell_g,grell_g_sh
![MLO - Adding shallow cumulus
  use shcu_vars_const   , only: nnshcu,shcufrq
!MLO]
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
  use mem_mass          , only: imassflx
!MLO]

  implicit none
  integer, intent(in) :: iens
! Two dimensions to start calling deep and shallow convection at different times
  real(kind=8), dimension(2), parameter :: CPTIME = (/7200.,0./) 

  !MLO. This is to use the PBL height that was found when running Nakanishi and Niino.
  !     If we didn't run Nakanishi and Niino, then use the default.
  
  

  if(iens == 1) then

     if (time == 0.) then
        call azero(mxp*myp*mzp,cuparm_g(ngrid)%THSRC(1,1,1))
        call azero(mxp*myp*mzp,cuparm_g(ngrid)%RTSRC(1,1,1))
        if (nnshcu(ngrid) == 0) then
           call azero(mxp*myp*mzp,cuparm_g_sh(ngrid)%THSRC(1,1,1))
           call azero(mxp*myp*mzp,cuparm_g_sh(ngrid)%RTSRC(1,1,1))
        end if
     end if

     if(initial == 2 .and. TIME < cptime(iens)-dble(dtlt)) return

     if(mod(time,dble(confrq)) < dtlt .or. time < .01 .or. abs(time-cptime(iens)) < 0.01)  &
     then
        !----- Remove part of the instability due to shallow cumulus action ---------------------!
        call include_shal_effect(mzp,mxp,myp,ia,iz,ja,jz,dtlt  &
                                ,sc1_grell_g(ngrid)%thetasta(1,1,1),sc1_grell_g(ngrid)%rvsta(1,1,1) &
                                ,basic_g(ngrid)%theta(1,1,1), basic_g(ngrid)%rv(1,1,1)           &
                                ,basic_g(ngrid)%pi0(1,1,1)  , basic_g(ngrid)%pp(1,1,1)           &
                                ,cuparm_g_sh(ngrid)%thsrc(1,1,1), cuparm_g_sh(ngrid)%rtsrc(1,1,1))

        call azero(mxp*myp*mzp,cuparm_g(ngrid)%thsrc(1,1,1))
        call azero(mxp*myp*mzp,cuparm_g(ngrid)%rtsrc(1,1,1))
        call azero(mxp*myp,cuparm_g(ngrid)%conprr(1,1))

        call cuparth(mynum,mgmxp,mgmyp,mgmzp,mzp,mxp,myp,ia,iz,ja,jz,i0,j0,maxiens         &
                    ,iens,icbase,depth_min(iens),cap_maxs(iens),DTLT,time                  &
                    ,basic_g(ngrid)%UP              , basic_g(ngrid)%VP                    &
!                    ,basic_g(ngrid)%WP              , basic_g(ngrid)%theta                 &
                    ,basic_g(ngrid)%WP              , sc1_grell_g(ngrid)%thetasta          &
                    ,basic_g(ngrid)%PP              , basic_g(ngrid)%PI0                   &
!                    ,basic_g(ngrid)%DN0             , basic_g(ngrid)%rv                    &
                    ,basic_g(ngrid)%DN0             , sc1_grell_g(ngrid)%rvsta             &
                    ,turb_g(ngrid)%kpbl             , turb_g(ngrid)%tkep  ,tkmin           &
                    ,micro_g(ngrid)%rcp             , grid_g(ngrid)%topt                   &
                    ,grid_g(ngrid)%RTGT             , tend%THT                             &
                    ,tend%RTT                       , tend%PT                              &
                    ,cuparm_g(ngrid)%THSRC          , cuparm_g(ngrid)%RTSRC                &
                    ,cuparm_g(ngrid)%CONPRR         , sc1_grell_g(ngrid)%ierr4d            &
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
                    ,grell_g(ngrid)%upmf            , grell_g(ngrid)%dnmf                  &
                    ,grell_g(ngrid)%xierr           , grell_g(ngrid)%xktop                 &
                    ,grell_g(ngrid)%xkbcon          , grell_g(ngrid)%xk22                  &
                    ,grell_g(ngrid)%xjmin           , grell_g(ngrid)%xkdt                  &
                    ,grell_g(ngrid)%xiact_p         , grell_g(ngrid)%xiact_c               )
! [MLO ------------- Stilt - RAMS coupling  ------------------
       if (imassflx == 1) then
         call prep_convflx_to_mass(mzp,mxp,myp,ia,iz,ja,jz,maxiens,ngrid                  &
            ,sc1_grell_g(ngrid)%ierr4d      ,sc1_grell_g(ngrid)%jmin4d                     &
            ,sc1_grell_g(ngrid)%kdet4d      ,sc1_grell_g(ngrid)%k224d                      &
            ,sc1_grell_g(ngrid)%kbcon4d     ,sc1_grell_g(ngrid)%ktop4d                     &
            ,sc1_grell_g(ngrid)%kpbl4d      ,sc1_grell_g(ngrid)%kstabi4d                   &
            ,sc1_grell_g(ngrid)%kstabm4d    ,sc1_grell_g(ngrid)%xmb4d                      &
            ,sc1_grell_g(ngrid)%edt4d       ,sc1_grell_g(ngrid)%zcup5d                     &
            ,sc1_grell_g(ngrid)%pcup5d      ,sc1_grell_g(ngrid)%enup5d                     &
            ,sc1_grell_g(ngrid)%endn5d      ,sc1_grell_g(ngrid)%deup5d                     &
            ,sc1_grell_g(ngrid)%dedn5d      ,sc1_grell_g(ngrid)%zup5d                      &
            ,sc1_grell_g(ngrid)%zdn5d       ,iens                                          )
       end if ! (imassflx == 1)
! ------------- Stilt - RAMS coupling  ------------------ MLO]

     end if !(mod(time,dble(confrq)) < dtlt .or. time < .01 .or. abs(time-cptime) < 0.01)

     call ACCUM(mxp*myp*mzp, tend%tht(1), cuparm_g(ngrid)%THSRC(1,1,1))
     call ACCUM(mxp*myp*mzp, tend%rtt(1), cuparm_g(ngrid)%RTSRC(1,1,1))

     call UPDATE(mxp*myp, cuparm_g(ngrid)%ACONPR(1,1),cuparm_g(ngrid)%CONPRR(1,1),DTLT)

!----- Shallow cumulus -------------------
  else if (iens == 2) then

     if(time == 0.) then !004
        call azero(mxp*myp*mzp,cuparm_g_sh(ngrid)%THSRC(1,1,1))
        call azero(mxp*myp*mzp,cuparm_g_sh(ngrid)%RTSRC(1,1,1))
        if (nnqparm(ngrid) /= 2 .or. time < cptime(1)-dble(dtlt) ) then
          call aone(mxp*myp,grell_g(ngrid)%xierr(1,1))
        end if
     end if 

     if(initial == 2 .and. time < cptime(iens)-dble(dtlt)) return

     if(mod(TIME,dble(shcufrq)) < dtlt .or. time < .01 .or. abs(time-cptime(iens)) < 0.01) &
     then

        call azero(mxp*myp*mzp,cuparm_g_sh(ngrid)%thsrc(1,1,1))
        call azero(mxp*myp*mzp,cuparm_g_sh(ngrid)%rtsrc(1,1,1))           

        call cuparth_shal(mynum,mgmxp,mgmyp,mgmzp,mzp,mxp,myp,ia,iz,ja,jz,i0,j0,maxiens,iens &
                         ,icbase,depth_min(iens),cap_maxs(iens),dtlt,time                    &
                         ,basic_g(ngrid)%up            ,basic_g(ngrid)%vp                    &
                         ,basic_g(ngrid)%wp            ,basic_g(ngrid)%theta                 &
                         ,basic_g(ngrid)%pp            ,basic_g(ngrid)%pi0                   &
                         ,basic_g(ngrid)%dn0           ,basic_g(ngrid)%rv                    &
                         ,turb_g(ngrid)%kpbl           ,turb_g(ngrid)%tkep,  tkmin           &
                         ,micro_g(ngrid)%rcp           ,grid_g(ngrid)%topt                   &
                         ,grid_g(ngrid)%rtgt           ,tend%tht                             &
                         ,tend%rtt                     ,tend%pt                              &
                         ,cuparm_g_sh(ngrid)%thsrc     ,cuparm_g_sh(ngrid)%rtsrc             &
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
                         ,grell_g_sh(ngrid)%upmf       ,grell_g_sh(ngrid)%xierr              &
                         ,grell_g_sh(ngrid)%xktop      ,grell_g_sh(ngrid)%xkbcon             &
                         ,grell_g_sh(ngrid)%xk22       ,grell_g(ngrid)%xierr                 &
                         )
! [MLO ------------- Stilt - RAMS coupling  ------------------
       if (imassflx == 1) then
         call prep_convflx_to_mass(mzp,mxp,myp,ia,iz,ja,jz,maxiens,ngrid                  &
            ,sc1_grell_g(ngrid)%ierr4d      ,sc1_grell_g(ngrid)%jmin4d                     &
            ,sc1_grell_g(ngrid)%kdet4d      ,sc1_grell_g(ngrid)%k224d                      &
            ,sc1_grell_g(ngrid)%kbcon4d     ,sc1_grell_g(ngrid)%ktop4d                     &
            ,sc1_grell_g(ngrid)%kpbl4d      ,sc1_grell_g(ngrid)%kstabi4d                   &
            ,sc1_grell_g(ngrid)%kstabm4d    ,sc1_grell_g(ngrid)%xmb4d                      &
            ,sc1_grell_g(ngrid)%edt4d       ,sc1_grell_g(ngrid)%zcup5d                     &
            ,sc1_grell_g(ngrid)%pcup5d      ,sc1_grell_g(ngrid)%enup5d                     &
            ,sc1_grell_g(ngrid)%endn5d      ,sc1_grell_g(ngrid)%deup5d                     &
            ,sc1_grell_g(ngrid)%dedn5d      ,sc1_grell_g(ngrid)%zup5d                      &
            ,sc1_grell_g(ngrid)%zdn5d       ,iens                                          )
       end if
! ------------- Stilt - RAMS coupling  ------------------ MLO]
     end if

     call accum(mxp*myp*mzp, tend%tht(1), cuparm_g_sh(ngrid)%thsrc(1,1,1))
     call accum(mxp*myp*mzp, tend%rtt(1), cuparm_g_sh(ngrid)%rtsrc(1,1,1))

  end if
  !MLO- Flushing scratch variables to zero, otherwise nested grids can get left-overs
  !     from the previous grid.
  call zero_scratch2_grell()
  call zero_scratch3_grell()
  call zero_scratch2_grell_sh() !LFR
  call zero_scratch3_grell_sh() !LFR
  return

end subroutine cuparm_grell
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
