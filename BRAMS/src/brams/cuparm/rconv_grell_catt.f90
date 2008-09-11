!------------------------------------------------------------------------------
! Cumulus Parameterization by G. Grell
!-Convective transport for smoke aerosols and non-hygroscopic gases
!-Developed by Saulo Freitas (sfreitas@cptec.inpe.br)
!-ref: Freitas, S.R., et al.: Monitoring the transport of biomass burning
!      emissions in South America. Environmental Fluid Mechanics,
!      Kluwer Academic Publishers, 2005.
!------------------------------------------------------------------------------

subroutine CUPARM_GRELL_CATT(iens)

  use mem_basic         , only: basic_g
  use mem_tend          , only: tend
  use mem_cuparm        , only: nnqparm,confrq,cuparm_g,cuparm_g_sh

  use node_mod          , only: mynum,   &   ! INTENT(IN)
             		        mxp,     &   ! INTENT(IN)
             		        myp,     &   ! INTENT(IN)
             		        mzp,     &   ! INTENT(IN)
             		        ia,      &   ! INTENT(IN)
             		        iz,      &   ! INTENT(IN)
             		        ja,      &   ! INTENT(IN)
             		        jz,      &   ! INTENT(IN)
             		        i0,      &   ! INTENT(IN)  
             		        j0	     ! INTENT(IN) 
  use mem_grid          , only: time,    &   ! INTENT(IN)
            		   	initial, &   ! INTENT(IN)
            		   	dtlt,	 &   ! INTENT(IN)
            		   	itimea,  &   ! INTENT(IN)
            		   	ngrid,   &   ! INTENT(IN)
            		   	grid_g,  &   ! INTENT(IN)  	  
            		   	dtlongn, &   ! INTENT(IN)
           		   	deltaxn, &   ! INTENT(IN)
           		   	deltayn, &
				npatch
				

  use rconstants        , only: tkmin 
  use extras            , only: extra3d,extra2d 
  use mem_turb          , only: turb_g,idiffk
  use mem_micro         , only: micro_g
  use mem_scratch       , only: scratch
  use mem_scalar        , only: scalar_g
!srf
  use io_params         , only: frqanl
  use mem_leaf          , only: leaf_g
  use micphys           , only: level
  
![MLO -> To allow writing convective mass fluxes and call shallow cumulus in a different frequency
  use mem_mass         , only: imassflx
  use shcu_vars_const   , only: nnshcu,shcufrq
! MLO]


!- use modules for Grell Parameterization
  use mem_grell_param       , only: mgmxp,mgmyp,mgmzp,maxiens,ngrids_cp,depth_min,icbase,cap_maxs
  use mem_scratch1_grell    , only: sc1_grell_g,iruncon
  use mem_scratch1_grell    , only: sc1_grell_g
  use mem_scratch2_grell    , only: zero_scratch2_grell
  use mem_scratch2_grell_sh , only: zero_scratch2_grell_sh
  use mem_scratch3_grell    , only: zero_scratch3_grell
  use mem_scratch3_grell_sh , only: zero_scratch3_grell_sh
  use mem_grell             , only: grell_g,grell_g_sh
  
  implicit none

  integer, intent(in) :: iens
  real(kind=8),parameter,dimension(2) :: cptime = (/7200.,0./) !orig: CPTIME = 7200.
!!$  integer,parameter :: i_forcing = 1


  !------------------------ deep convection -----------------------------------
  if(iens == 1) then 
  
  !
  !        Zero out tendencies initially
     if (time == 0.) then 
        call azero(mxp*myp*mzp,cuparm_g(ngrid)%thsrc(1,1,1))
        call azero(mxp*myp*mzp,cuparm_g(ngrid)%rtsrc(1,1,1))
        call azero(mxp*myp*mzp,extra3d(1,ngrid)%d3(1,1,1)  )
        call azero(mxp*myp    ,extra2d(3,ngrid)%d2(1,1)    )  
        if (nnshcu(ngrid) == 0) then
           call azero(mxp*myp*mzp,cuparm_g_sh(ngrid)%THSRC(1,1,1))
           call azero(mxp*myp*mzp,cuparm_g_sh(ngrid)%RTSRC(1,1,1))
           call azero(mxp*myp    ,extra2d(2,ngrid)%d2(1,1)       )  
        end if
     end if 

     if(initial == 2 .and. TIME < cptime(iens)-dble(dtlt)) return
     if(mod(time,dble(confrq)) < dtlt .or. time < .01 .or. abs(time-cptime(iens)) < 0.01) then

        iruncon=1
        !----- Remove part of the instability due to shallow cumulus action ---------------------!
        call include_shal_effect(mzp,mxp,myp,ia,iz,ja,jz,dtlt                                    &
                                ,sc1_grell_g(ngrid)%thetasta(1,1,1),sc1_grell_g(ngrid)%rvsta(1,1,1) &
                                ,basic_g(ngrid)%theta(1,1,1), basic_g(ngrid)%rv(1,1,1)           &
                                ,basic_g(ngrid)%pi0(1,1,1)  , basic_g(ngrid)%pp(1,1,1)           &
                                ,cuparm_g_sh(ngrid)%thsrc(1,1,1), cuparm_g_sh(ngrid)%rtsrc(1,1,1))

        call azero(mxp*myp*mzp,cuparm_g(ngrid)%thsrc(1,1,1))
        call azero(mxp*myp*mzp,cuparm_g(ngrid)%rtsrc(1,1,1))
        call azero(mxp*myp    ,cuparm_g(ngrid)%conprr(1,1) )
        call azero(mxp*myp*mzp,extra3d(1,ngrid)%d3(1,1,1)  )  
        call azero(mxp*myp*mzp,extra3d(5,ngrid)%d3(1,1,1)  )!cloud/ice tendency
        call azero(mxp*myp    ,extra2d(3,ngrid)%d2(1,1)    )  

        if(mod(time,dble(frqanl)) < dtlongn(1) ) then
             call azero(mxp*myp*mzp,extra3d(4,ngrid)%d3(1,1,1)  )  
	endif

        call cuparth_catt(                 &
             mynum,                        & !01
             mgmxp,                        & !02
             mgmyp,                        & !03
             mgmzp,                        & !04
             mzp,                          & !05
             mxp,                          & !06
             myp,                          & !07
             ia,                           & !08
             iz,                           & !09
             ja,                           & !10
             jz,                           & !11
             i0,                           & !12
             j0,                           & !13
             maxiens,                      & !14
             iens,                         & !15
             ngrid,                        & !16
             ngrids_cp,                    & !17
             icbase,                       & !17
             depth_min(iens),              & !17
             cap_maxs(iens),               & !17
             DTLT,                         & !18
             time,                         & !19
             !
             basic_g(ngrid)%up   (1,1,1),    & !20
             basic_g(ngrid)%vp   (1,1,1),    & !21
             basic_g(ngrid)%wp   (1,1,1),    & !22
             sc1_grell_g(ngrid)%thetasta(1,1,1),& !23
             basic_g(ngrid)%pp   (1,1,1),    & !24
             basic_g(ngrid)%pi0  (1,1,1),    & !25
             basic_g(ngrid)%dn0  (1,1,1),    & !26
             sc1_grell_g(ngrid)%rvsta(1,1,1),   & !27
             turb_g(ngrid)%kpbl(1,1),        & !27½
             turb_g(ngrid)%tkep  (1,1,1),    & !28 *
             tkmin,                          & !29 *
             micro_g(ngrid)%rcp  (1,1,1),    & !30 *
	     grid_g(ngrid)%topt  (1,1),      & !31
             grid_g(ngrid)%RTGT  (1,1),      & !32
             !
             !grell_g(ngrid)%lsfth(1,1,1)  ,  & !33 
             !grell_g(ngrid)%lsfrt(1,1,1)  ,  & !34 
             !
             tend%THT(1),                    & !33 
             tend%RTT(1),                    & !34 
             !
             tend%PT(1),                     & !35
             cuparm_g(ngrid)%THSRC (1,1,1),  & !36 
             cuparm_g(ngrid)%RTSRC (1,1,1),  & !37 
             cuparm_g(ngrid)%CONPRR(1,1),    & !38	
             extra3d(1,ngrid)%d3   (1,1,1),  & !39
             extra3d(4,ngrid)%d3   (1,1,1),  & !40
             extra3d(5,ngrid)%d3   (1,1,1),  & !41
             extra2d(3,ngrid)%d2   (1,1),    & !42
             !
             sc1_grell_g(ngrid)%ierr4d(1,1,1),     & !43
             sc1_grell_g(ngrid)%jmin4d(1,1,1),     & !44
             sc1_grell_g(ngrid)%kdet4d(1,1,1),     & !45
             sc1_grell_g(ngrid)%k224d(1,1,1),      & !46
             sc1_grell_g(ngrid)%kbcon4d(1,1,1),    & !47
             sc1_grell_g(ngrid)%ktop4d(1,1,1),     & !48
             sc1_grell_g(ngrid)%kpbl4d(1,1,1),     & !49
             sc1_grell_g(ngrid)%kstabi4d(1,1,1),   & !50
             sc1_grell_g(ngrid)%kstabm4d(1,1,1),   & !51
             sc1_grell_g(ngrid)%xmb4d(1,1,1),      & !52
             sc1_grell_g(ngrid)%edt4d(1,1,1),      & !53
             !
             sc1_grell_g(ngrid)%zcup5d(1,1,1,1),   & !54 
             sc1_grell_g(ngrid)%pcup5d(1,1,1,1),   & !55 
             sc1_grell_g(ngrid)%enup5d(1,1,1,1),   & !56
             sc1_grell_g(ngrid)%endn5d(1,1,1,1),   & !57
             sc1_grell_g(ngrid)%deup5d(1,1,1,1),   & !58 
             sc1_grell_g(ngrid)%dedn5d(1,1,1,1),   & !59
             sc1_grell_g(ngrid)%zup5d(1,1,1,1),    & !60
             sc1_grell_g(ngrid)%zdn5d(1,1,1,1),    & !61
             sc1_grell_g(ngrid)%prup5d(1,1,1,1),   & !62
             sc1_grell_g(ngrid)%clwup5d(1,1,1,1),  & !63
             sc1_grell_g(ngrid)%tup5d(1,1,1,1),    & !64
             !
             grell_g(ngrid)%upmf   (1,1),  & !65
             grell_g(ngrid)%dnmf   (1,1),  & !66
             grell_g(ngrid)%xierr  (1,1),  & !67
             grell_g(ngrid)%xktop  (1,1),  & !68
             grell_g(ngrid)%xkbcon (1,1),  & !69
             grell_g(ngrid)%xk22   (1,1),  & !70
             grell_g(ngrid)%xjmin  (1,1),  & !71
             grell_g(ngrid)%xkdt   (1,1),  & !72
             grell_g(ngrid)%xiact_p(1,1),  & !73
             grell_g(ngrid)%xiact_c(1,1),  & !74
	     confrq,                       & !75
             frqanl,                       & !76
	     deltaxn(ngrid)*deltayn(ngrid),& !77
             leaf_g(ngrid)%patch_area(1,1,1), &!78
	     npatch,                       & !79
	     level)                          !80

     end if 
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

     call accum(mxp*myp*mzp, tend%tht(1), cuparm_g(ngrid)%thsrc(1,1,1))
     call accum(mxp*myp*mzp, tend%rtt(1), cuparm_g(ngrid)%rtsrc(1,1,1))
     call update(mxp*myp, cuparm_g(ngrid)%aconpr(1,1),cuparm_g(ngrid)%conprr(1,1),dtlt)

  !---------------   Shallow cumulus scheme -----------------------------------
  else if(iens == 2) then
     if(TIME == 0.) then
        call azero(mxp*myp*mzp,cuparm_g_sh(ngrid)%thsrc(1,1,1))
        call azero(mxp*myp*mzp,cuparm_g_sh(ngrid)%rtsrc(1,1,1))
! [MLO -> this is to allow the run to call the shallow cumulus if the deep wasn't called
        if (nnqparm(ngrid) /= 2 .or. time < cptime(1)-dble(dtlt) ) then
          call aone(mxp*myp,grell_g(ngrid)%xierr(1,1))
        end if
!MLO]
     end if 

     if(initial == 2 .and. time < cptime(iens)-dble(dtlt)) return
![MLO -> Changed to shcufrq to allow different scales
     if(mod(TIME,dble(shcufrq)) < dtlt .or. time < .01 .or. abs(time-cptime(iens)) <  0.01) then
! MLO]
        iruncon=1

        call azero(mxp*myp*mzp,cuparm_g_sh(ngrid)%thsrc(1,1,1))
        call azero(mxp*myp*mzp,cuparm_g_sh(ngrid)%rtsrc(1,1,1))           
        call azero(mxp*myp    ,extra2d(2,ngrid)%d2(1,1)       )  


        call cuparth_shal_catt(       &
             mynum,           &       !01
             mgmxp,           &       !02
             mgmyp,           &       !03
             mgmzp,           &       !04
             mzp,             &       !05
             mxp,             &       !06
             myp,             &       !07
             ia,              &       !08
             iz,              &       !09
             !
             ja,              &       !10
             jz,              &       !11
             i0,              &       !12
             j0,              &       !13
             maxiens,         &       !15
             iens,            &       !16
             ngrid,           &       !17
             ngrids_cp,       &       !18
             icbase,          &       !18¼
             depth_min(iens), &       !18½
             cap_maxs(iens),  &       !18¾
             dtlt,            &       !19
             time,            &       !20
             !
             basic_g(ngrid)%up   (1,1,1),     &   !21
             basic_g(ngrid)%vp   (1,1,1),     &   !22
             basic_g(ngrid)%wp   (1,1,1),     &   !23
             basic_g(ngrid)%theta(1,1,1),     &   !24
             basic_g(ngrid)%pp   (1,1,1),     &   !25
             basic_g(ngrid)%pi0  (1,1,1),     &   !26
             basic_g(ngrid)%dn0  (1,1,1),     &   !27
             basic_g(ngrid)%rv   (1,1,1),     &   !28
             turb_g(ngrid)%kpbl(1,1),         &   !28½
             turb_g(ngrid)%tkep  (1,1,1),     &   !29
             !
             tkmin,                          &   !30
             micro_g(ngrid)%rcp     (1,1,1), &   !31
             grid_g(ngrid)%topt     (1,1),   &   !32
             grid_g(ngrid)%RTGT     (1,1),   &   !33
             !
             !grell_g_sh(ngrid)%lsfth(1,1,1), &   !34
             !grell_g_sh(ngrid)%lsfrt(1,1,1), &   !35
             !
             tend%THT(1),                    &   !34
             tend%RTT(1),                    &   !35
             !
             tend%PT(1),		     &   !36
             cuparm_g_sh(ngrid)%thsrc(1,1,1),&   !37
             cuparm_g_sh(ngrid)%rtsrc(1,1,1),&   !38
             extra3d(2,ngrid)%d3     (1,1,1),&   !39 !<< usando extra3d(2)
             extra2d(2,ngrid)%d2     (1,1),  &   !39 !<< usando extra2d(2)
             !
             sc1_grell_g(ngrid)%ierr4d(1,1,1),      & !40
             sc1_grell_g(ngrid)%jmin4d(1,1,1),      & !41
             sc1_grell_g(ngrid)%kdet4d(1,1,1),      & !42
             sc1_grell_g(ngrid)%k224d(1,1,1),       & !43
             sc1_grell_g(ngrid)%kbcon4d(1,1,1),     & !44
             sc1_grell_g(ngrid)%ktop4d(1,1,1),      & !45
             sc1_grell_g(ngrid)%kpbl4d(1,1,1),      & !46
             sc1_grell_g(ngrid)%kstabi4d(1,1,1),    & !47
             sc1_grell_g(ngrid)%kstabm4d(1,1,1),    & !48
             !
             sc1_grell_g(ngrid)%xmb4d(1,1,1),       & !49
             sc1_grell_g(ngrid)%edt4d(1,1,1),       & !50
             sc1_grell_g(ngrid)%zcup5d(1,1,1,1),    & !51
             sc1_grell_g(ngrid)%pcup5d(1,1,1,1),    & !52
             sc1_grell_g(ngrid)%enup5d(1,1,1,1),    & !53
             sc1_grell_g(ngrid)%endn5d(1,1,1,1),    & !54
             sc1_grell_g(ngrid)%deup5d(1,1,1,1),    & !55
             sc1_grell_g(ngrid)%dedn5d(1,1,1,1),    & !56
             sc1_grell_g(ngrid)%zup5d(1,1,1,1),     & !57
             sc1_grell_g(ngrid)%zdn5d(1,1,1,1),     & !58
             sc1_grell_g(ngrid)%prup5d(1,1,1,1),    & !59
             sc1_grell_g(ngrid)%clwup5d(1,1,1,1),   & !60
             sc1_grell_g(ngrid)%tup5d(1,1,1,1),     & !61
             !
             grell_g_sh(ngrid)%upmf  (1,1),  & !62
             grell_g_sh(ngrid)%xierr (1,1),  & !63
             grell_g_sh(ngrid)%xktop (1,1),  & !64
             grell_g_sh(ngrid)%xkbcon(1,1),  & !65
             grell_g_sh(ngrid)%xk22  (1,1),  & !66
             grell_g   (ngrid)%xierr (1,1),  & !67
	     confrq,shcufrq,                 &
	     deltaxn(ngrid)*deltayn(ngrid),  &
             leaf_g(ngrid)%patch_area(1,1,1),&
	     npatch,                         &
	     level)

     end if !005
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
     call accum(mxp*myp*mzp, tend%tht(1), cuparm_g_sh(ngrid)%thsrc(1,1,1))
     call accum(mxp*myp*mzp, tend%rtt(1), cuparm_g_sh(ngrid)%rtsrc(1,1,1))

  end if

  !--------- Convective Transport based on mass flux scheme -------------------

  if(iruncon == 1) then
 
     call azero(mxp*myp*mzp,scratch%scr1(1)) 

     call trans_conv_mflx(iens,scratch%scr1(1))

     if(iens == 2) call AE1M1(mxp*myp*mzp,extra3d(3,ngrid)%d3(1,1,1), &
                              scalar_g(1,ngrid)%sclt(1),              &
			      extra3d(3,ngrid)%d3(1,1,1)) 
  end if 

  !MLO- Flushing scratch variables to zero, otherwise nested grids can get left-overs
  !     from the previous grid.
  call zero_scratch2_grell()
  call zero_scratch3_grell()
  call zero_scratch2_grell_sh() !LFR
  call zero_scratch3_grell_sh() !LFR

  return

end subroutine CUPARM_GRELL_CATT
