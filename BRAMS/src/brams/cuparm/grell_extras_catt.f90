subroutine trans_conv_mflx(icld,iscl,m1,m2,m3,maxens_cap,i,j,ierr_cap,jmin_cap,k22_cap     &
                          ,kbcon_cap,kdet_cap,kstabi_cap,kstabm_cap,ktop_cap,cdd_cap       &
                          ,cdu_cap,mentrd_rate_cap,mentru_rate_cap,etad_cld_cap            &
                          ,etau_cld_cap,edt,upmf,sttend)
  !------------------------------------------------------------------
  !-Convective transport for smoke aerosols and non-hygroscopic gases
  !-Developed by Saulo Freitas (sfreitas@cptec.inpe.br)
  !-ref: Freitas, S.R., et al.: Monitoring the transport of biomass burning
  !      emissions in South America. Environmental Fluid Mechanics,
  !      Kluwer Academic Publishers, 2005.
  !------------------------------------------------------------------
  use mem_tconv         ,only: stcum1d,dn01d,se,se_cup,sc_up,sc_dn,sc_up_c,sc_dn_c, &
                               henry_coef,pw_up,pw_dn,&
                               trans_conv_alloc, &    !Control alloc variable
                               alloc_trans_conv, &    !alloc subroutine
                               zero_tconv

  use mem_scalar        ,only : scalar_g
  use mem_grid          ,only : dtlt,ngrid,naddsc
  use mem_scratch       ,only : scratch
  use mem_basic         ,only : basic_g
  use grell_coms        ,only : mgmzp
  use mem_ensemble      ,only : ensemble_vars
  use mem_scratch_grell ,only :  z_cup,kgoff
  

  implicit none

  integer, intent(in)                      :: m1,m2,m3,icld,iscl,i,j,maxens_cap
  real, dimension(m1,m2,m3), intent(inout) :: sttend
  real, intent(in)                         :: edt,upmf
  integer                                  :: icap

   integer, dimension      (maxens_cap), intent(in) :: ierr_cap
   integer, dimension      (maxens_cap), intent(in) :: jmin_cap
   integer, dimension      (maxens_cap), intent(in) :: k22_cap
   integer, dimension      (maxens_cap), intent(in) :: kbcon_cap
   integer, dimension      (maxens_cap), intent(in) :: kdet_cap
   integer, dimension      (maxens_cap), intent(in) :: kstabi_cap
   integer, dimension      (maxens_cap), intent(in) :: kstabm_cap
   integer, dimension      (maxens_cap), intent(in) :: ktop_cap
   real   , dimension(mgmzp,maxens_cap), intent(in) :: cdd_cap
   real   , dimension(mgmzp,maxens_cap), intent(in) :: cdu_cap
   real   , dimension(mgmzp,maxens_cap), intent(in) :: mentrd_rate_cap
   real   , dimension(mgmzp,maxens_cap), intent(in) :: mentru_rate_cap
   real   , dimension(mgmzp,maxens_cap), intent(in) :: etad_cld_cap
   real   , dimension(mgmzp,maxens_cap), intent(in) :: etau_cld_cap



  integer                                  :: k,kr,iconv,iwet
  real, dimension(2) :: c0

  c0(1) = 0.002
  c0(2) = 0.000
  
  
  if (upmf == 0. .or. all(ierr_cap /= 0)) return

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!! THIS IS NOT RIGHT... Yet... The right way to do this is averaging all non-zero !!!!
  !!!!!! members.                                                                       !!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  stacloop: do icap=1,maxens_cap
     if (ierr_cap(icap) == 0) exit stacloop
  end do stacloop

!  data (c0(i),i=1,2)  &
!       /0.002 &  
!       ,0.000 /  

  if(.not. trans_conv_alloc) then
     call alloc_trans_conv(mgmzp)
     call zero_tconv()
     !PRINT*,'======================================================='
     !PRINT*,'convective transport: memory allocation'
     !PRINT*,'======================================================='
  end if


  iwet = 0
  !---
  !02/05/2005: COstc nao tera apenas plumerise
  if(iscl == 2) return    ! sem transporte convectivo
  !---
  if(iscl == 3 .and. icld == 1) iwet = 1   ! com deposicao umida

  call azero(mgmzp,stcum1d)
  iconv = 0       
  if(ierr_cap(icap) == 0) then
     
     iconv = 1
     call get_dn01d(mgmzp,m1,dn01d,basic_g(ngrid)%dn0(:,i,j))
     call get_se(mgmzp,m1,m2,m3,i,j,scalar_g(iscl,ngrid)%sclp,se,se_cup)

     if(iwet == 0) then
        call get_sc_up( mgmzp,m1,se,se_cup,sc_up,    &
             k22_cap(icap),kbcon_cap(icap),cdu_cap(1:m1,icap), &
             mentru_rate_cap(1:m1,icap),etau_cld_cap(1:m1,icap))
     else
        ! print*,'I J PREC:',i,j
        ! print*,'PM25=',scalar_g(iscl,ngrid)%sclp(:,I,J)
                       
        call get_sc_up_wet( mgmzp,m1,se,se_cup,sc_up,   &
             k22_cap(icap), kbcon_cap(icap), ktop_cap(icap), &
             cdu_cap(1:m1,icap),mentru_rate_cap(1:m1,icap),  &
             etau_cld_cap(1:m1,icap),sc_up_c,henry_coef,pw_up,dn01d, &
             c0(icld),etau_cld_cap(1:m1,icap) )

     endif

     call get_sc_dn(mgmzp,m1,se,se_cup,sc_dn,jmin_cap(icap),kdet_cap(icap), &
                    cdd_cap(:,icap),mentrd_rate_cap(:,icap), &
                    etad_cld_cap(:,icap)  )

     call get_stcum(mgmzp,m1,dn01d,stcum1d,se,se_cup,sc_up,sc_dn,                      &
          upmf  , edt, jmin_cap(icap)  , kdet_cap(icap), k22_cap(icap)   ,    &
          kbcon_cap(icap), ktop_cap(icap),  kstabi_cap(icap),                 &
          kstabm_cap(icap), z_cup(1:m1), cdu_cap(1:m1,icap), &
          mentru_rate_cap(1:m1,icap), cdd_cap(1:m1,icap),    &
          mentrd_rate_cap(1:m1,icap), etau_cld_cap(1:m1,icap), &
          etad_cld_cap(1:m1,icap))

     if(iwet == 1) then
        !------------- WET REMOVAL SCHEMES -------------
        !-- contribuicao devido `a concentracao aquosa desentranhada
        !-- junto com a agua liquida
        call get_stcum_detrain(mgmzp,m1,dn01d,stcum1d,sc_up_c,sc_dn_c,upmf, edt,          &
                               jmin_cap(icap), kdet_cap(icap), k22_cap(icap),    &
                               kbcon_cap(icap),ktop_cap(icap), kstabi_cap(icap), &
                               kstabm_cap(icap), z_cup(1:m1), cdu_cap(:,icap),      &
                               mentru_rate_cap(1:m1,icap), cdd_cap(:,icap),         &
                               mentrd_rate_cap(:,icap), etau_cld_cap(:,icap),       &
                               etad_cld_cap(:,icap))

        !-- calcula a massa removida e depositada na superficie
        call get_wet_deposited_mass(mgmzp,m1,dtlt,dn01d,pw_up,pw_dn,                   &
             scalar_g(iscl,ngrid)%wetdep(i,j),                                         &
             upmf,  edt,kbcon_cap(icap),ktop_cap(icap))
        !print*,'I J WM:',i,j, scalar_g(iscl,ngrid)%wetdep(i,j)


     endif

  endif    ! endif of test if is there or not convection
   if(iconv .eq. 1) then
     do k=1,m1-1
        kr= k + kgoff !nivel k da grade do grell corresponde ao nivel k+kgoff do rams
        sttend(kr,i,j)= sttend(kr,i,j) + stcum1d(k)
     enddo
  endif  !endif do iconv

  return
end subroutine trans_conv_mflx

!--------------------------------------------------

subroutine get_se(mgmzp,n1,n2,n3,i,j,sp,se,se_cup)

  implicit none

  ! Arguments
  integer :: mgmzp, n1, n2, n3, i, j
  real    :: sp(n1,n2,n3), se(mgmzp), se_cup(mgmzp)

  !local Variables
  integer :: k, kr

  do k=2,n1-1
     kr = k+1   ! nivel K da grade DO Grell corresponde ao nivel K + 1 DO RAMS
     se(k) = sp(kr,i,j)                         ! concentration on Z levels
     se_cup(k) = .5*( sp(kr-1,i,j) + sp(kr,i,j) )   ! concentration on Z_CUP levels
  enddo
  se(1)     = sp(2,i,j)
  se_cup(1) = sp(2,i,j)

end subroutine get_se

!----------------------------------------------------------------------

subroutine get_sc_up( mgmzp,n1,se,se_cup,sc_up,k22,kbcon,cd,entr,z_cup)

  implicit none
  integer :: mgmzp,n1,k22,kbcon
  real, dimension(n1)    :: cd,entr,z_cup
  real, dimension(mgmzp) :: se,se_cup,sc_up
  real ::dz
  integer :: k

  ! Scalar concentration in-cloud - updraft

  do k=1,k22-1
     sc_up(k) = se_cup(k)
  enddo

  do k=k22,kbcon
     sc_up(k) = se_cup(k22)
  enddo

  do k=kbcon+1,n1-1

     dz=z_cup(k)-z_cup(k-1)
     sc_up(k)= ((1.-.5*cd(k)*dz)*sc_up(k-1)+entr(k)*dz*se(k-1)) / &
                (1. + entr(k)*dz - .5*cd(k)*dz )
     !expression II based on dH_c/dz = entr*(H_env - H_c)
     !  sc_up(k)= ( (1.- .5*entr(k)*dz)*sc_up(k-1) + entr(k)*dz*se(k-1) ) &
     !              / ( 1. + .5*entr(k)*dz )
  enddo

end subroutine get_sc_up

!-------------------------------------------------------------------------

subroutine get_sc_dn( mgmzp,n1,se,se_cup,sc_dn,jmin,kdet,cdd,entrd,z_cup )

  implicit none

  ! Arguments
  integer                :: mgmzp, n1, jmin, kdet
  real, dimension(n1)    :: cdd, entrd, z_cup
  real, dimension(mgmzp) ::  se, se_cup, sc_dn

  ! Local variables
  integer                :: k
  real                   :: dz


  !Scalar concentration in-cloud - DOwndraft
  do k=jmin+1,n1-1
     sc_dn(k) = se_cup(k)
  enddo

  !k=jmim
  sc_dn(jmin) = se_cup(jmin)

  !added fev-2003 for shallow scheme
  if(jmin == 1 ) return

  do k=jmin-1,1,-1
     dz=z_cup(k+1)-z_cup(k)
     sc_dn(k) = ((1.-.5*cdd(k)*dz)*sc_dn(K+1) + entrd(k)*dz*se(k)) /&
          (1. + entrd(k)*dz - .5*cdd(k)*dz )
  enddo

end subroutine get_sc_dn

!------------------------------------------------------------------------


subroutine get_stcum(mgmzp,n1,dn0,stcum1d,se,se_cup,sc_up,sc_dn,xmb,edt, &
     jmin,kdet,k22,kbcon,ktop,kstabi,kstabm,z_cup,                  &
     cd,entr,cdd,entrd,zu,zd                                       )

  use rconstants, only: g

  implicit none

  ! Arguments
  integer                :: mgmzp, n1, jmin, kdet, k22, &
       kbcon, ktop, kstabi, kstabm
  real                   :: xmb,edt
  real, dimension(mgmzp) :: se,se_cup,sc_up,sc_dn,stcum1d,dn0
  real, dimension(n1) :: z_cup,cd,entr,cdd,entrd,zu,zd

  ! Local variables
  integer                :: k
  real                   :: dz, dp, entup, detup, entdoj, entupk, detupk, &
       detdo, entdo, subin, subdown, detdo1, detdo2

  do k=2,ktop

     dz =   z_cup(k+1) - z_cup(k)

     dp      = g*dn0(k)*dz  
     entup   = 0.
     detup   = 0.
     entdoj  = 0.
     entupk  = 0.
     detupk  = 0.
     detdo   = edt*  cdd(k)*dz*zd(k+1)
     entdo   = edt*entrd(k)*dz*zd(k+1)
     subin   = zu(k+1) - edt*zd(k+1)
     subdown = zu(k  ) - edt*zd(k  )

     if(k.ge.kbcon .and. k.lt.ktop)   then
        detup =   cd(k+1) *dz*zu(k)
        entup = entr(k)   *dz*zu(k)
     endif

     if(k.eq.jmin)          entdoj = edt*zd(k)
     if(k.eq.k22-1)         entupk = zu(k22)
     if(k.gt.kdet)          detdo  = 0.
     if(k.lt.kbcon)         detup  = 0.
     if(k.eq.ktop) then
        detupk  = zu(ktop)
        subin   = 0.
     endif

     !tendency due cumulus transport ( k>=2)
     stcum1d(k) = stcum1d(k) +                         & !stcum1d comparece
                  xmb*( 			       & !no lado direito
                  subin  *se_cup(k+1)		       & !para computar
                  - subdown*se_cup(k)		       & !a soma sobre o
                  + detup*( sc_up(k+1)+ sc_up(k) )*.5  & !espectro de nuvens
                  + detdo*( sc_dn(k+1)+ sc_dn(k) )*.5  & !(maxiens)
                  - entup*se(k) 		       &
                  - entdo*se(k) 		       &
                  - entupk*se_cup(k22)  	       &
                  - entdoj*se_cup(jmin) 	       &
                  + detupk* sc_up(ktop) 	       &
                  )*g/dp

  end do
  !
  !tendency due cumulus transport (at bottom, k = 1)
  dz        =       z_cup(2)-z_cup(1)
  dp =   g*dn0(1)*dz ! ja' esta' na grade do rams

  detdo1    = edt*zd(2)*cdd(1)*dz
  detdo2    = edt*zd(1)
  entdo     = edt*zd(2)*entrd(1)*dz
  subin     =-edt*zd(2)

  stcum1d(1) = stcum1d(1) +                      &
               xmb*(				 &
               detdo1*(sc_dn(1)+sc_dn(2))*.5	 &
               + detdo2*    sc_dn(1)		 &
               + subin *    se_cup(2)		 &
               - entdo *    se(1)		 &
               )*g/dp

  return
end subroutine get_stcum

!-----------------------------------------------

subroutine get_dn01d(mgmzp,m1,dn01d,dn0)

  implicit none
  integer mgmzp,m1,k
  real, dimension(mgmzp) :: dn01d
  real, dimension(m1)    :: dn0

  do k=1,m1-1
     dn01d(k) = dn0(k+1)  ! nivel K da grade DO Grell corresponde ao nivel K + 1 DO RAMS
  enddo
  dn01d(m1)=dn01d(m1-1)

end subroutine get_dn01d

!---------------------------------------------------------------------------

subroutine get_sc_up_wet( mgmzp,n1,se,se_cup,sc_up,k22,kbcon,ktop,cd, &
     entr,z_cup,sc_up_c,henry_coef,pw_up,dn0,c0,zu)

 
  implicit none
  integer mgmzp,n1,k,k22,kbcon,ktop,iwd,iall
  real dz,conc_equi,conc_mxr,qrch,c0,scav_eff
  real, dimension(n1) :: cd,entr,z_cup, &
       dn0,        & ! basic state air density [kg[air]/m^3]
       zu          ! ! norm mass flux updrat
  
  real, dimension(mgmzp) :: se,  &
       se_cup ,    &! gas-phase in environment mixing ratio [kg(gas phase)/kg(air)]
       sc_up        ! gas-phase in updraft     mixing ratio [kg(gas phase)/kg(air)]

  real, dimension(mgmzp) :: &
       henry_coef, & ! Henry's constant [(kg(aq)/kg(water))/(kg(gas phase)/kg(air))]
       sc_up_c,    & ! aqueous-phase in updraft mixing ratio [kg(aq)/kg(air)]
       pw_up         ! precitable gas/aer amount [kg[..]/kg[air]
  data scav_eff/0.6/ ! for smoke : Chuang et al. (1992) J. Atmos. Sci.

  iwd =1
  iall=0

  call azero3(mgmzp,sc_up_c,pw_up,henry_coef)


  do k=1,k22-1
     sc_up(k) = se_cup(k)
  enddo

  do k=k22,kbcon
     sc_up(k) = se_cup(k22)
  enddo

  do k=kbcon+1,ktop
     dz=z_cup(k)-z_cup(k-1)
     !
     !
     sc_up(k)= ((1.-.5*cd(k)*dz)*sc_up(k-1)+entr(k)*dz*se(k-1))/ &
               ( 1. + entr(k)*dz - .5*cd(k)*dz )

     if(iwd.eq.1)then
        conc_mxr  =  scav_eff* sc_up(k)    
     else
        conc_mxr=0.       
     endif

     !
     qrch = sc_up(k) - conc_mxr          

     sc_up_c(k) = (sc_up(k)-qrch)/(1.+c0*dz*zu(k)) 
                                                   
     !
     pw_up(k) = c0*dz*sc_up_c(k)*zu(k)      
     !
     if(sc_up_c(k).lt.0.) sc_up_c(k) = 0.
     !
     !
     if(iall.eq.1)then
        sc_up_c(k) = 0.
        pw_up(k)	 = max( (sc_up(k)-QRCH)*zu(k), 0. )
     endif
     !
     !----- set next level
     !
     sc_up(k) = sc_up_c(k) + qrch 
  enddo
  !-----

  do k=ktop+1,n1-1
     sc_up(k) = se_cup(k)
  enddo
  !
   do k=kbcon+1,ktop
     sc_up(k) = sc_up(k) - sc_up_c(k)
  enddo

end subroutine get_sc_up_wet

!---------------------------------------------------------------------

subroutine get_stcum_detrain(mgmzp,n1,dn0,stcum1d,sc_up_c,sc_dn_c,  &
     	  		     xmb,edt,jmin,kdet,k22,kbcon,ktop, &
     	  		     kstabi,kstabm,z_cup,cd,entr,cdd, &
     	  		     entrd,zu,zd)

  use rconstants, only: g  			

  implicit none

  ! Arguments
  integer                :: mgmzp,n1,jmin,kdet,k22,kbcon,ktop,kstabi,kstabm
  real                   :: xmb,edt
  real, dimension(mgmzp) :: sc_up_c,sc_dn_c,stcum1d,dn0
  real, dimension(n1)    :: z_cup,cd,entr,cdd,entrd,zu,zd

  ! Local Variables
  integer :: k
  real    :: dz, dp,detup, detupk

  do k=kbcon+1,ktop
     dz =   z_cup(k+1) - z_cup(k)
     dp =   g*dn0(k)*dz  
     detup  = 0.
     if(k.lt.ktop)  detup  = cd(k+1) *dz*zu(k)
     detupk = 0.
     if(k.eq.ktop)  detupk = zu(ktop)
     !tendency due cumulus transport ( k>=2)
     stcum1d(k) = stcum1d(k) +                             & ! stcum1d do lado direito corresponde 
                  xmb*( 				   & !         `a contribuicao da fase 
                  + detup*( sc_up_c(k+1)+ sc_up_c(k) )*.5  & !  	gasosa.
                  + detupk* sc_up_c(ktop)		   & ! stcum1d final = fase gasosa +
                  )*g/dp				     !  	       fase liquida

  end do
end subroutine get_stcum_detrain

!---------------------------------------------------------------------

subroutine get_wet_deposited_mass(mgmzp,m1,dt,dn01d,pw_up,pw_dn, &
     wetdep,xmb,edt,kbcon,ktop)

  implicit none
  integer k,mgmzp,m1,kbcon,ktop
  real dt,wetdep,xmb,edt
  real, dimension(mgmzp) :: dn01d,pw_up,pw_dn

  !--- wetdep units kg m^-2
  !!wetdep = 0.  ! use for instantaneous rate (not integrated over the time)
  do  k=1,ktop+1
     wetdep = wetdep + dt*(pw_up(k)+edt*pw_dn(k))*xmb
  !   print *,'pw_up,pw_dn,WM:',k,pw_up(k),pw_dn(k),wetdep
  enddo

end subroutine get_wet_deposited_mass
