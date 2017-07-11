!========================================================================
! Dry deposition and sedimentation of particles
! Coded and implemented by Saulo Freitas and Karla Longo
! (CPTEC/INPE: sfreitas@cptec.inpe.br ; longo@cptec.inpe.br)
! Ref.: Longo, K. M.; Freitas, S. R.; Silva Dias, M.A.F. Dias, P. Silva Dias. Numerical
! modeling developments towards a system suitable to a real time air quality forecast
! and climate changes studies in South America. Newsletter of the International Global
! Atmospheric Chemistry Project, Taiwan, v. 33, p. 12-16, 2006.
!========================================================================

subroutine drydep_driver(m1,m2,m3,ia,iz,ja,jz)
  
  USE mem_grid         ,  ONLY :  npatch,naddsc,MAXGRDS,nzpmax,dzt,zt,ngrids  &
                                 ,ngrid,nzg,nzs,jdim,dtlt,grid_g
  USE mem_basic        ,  ONLY :  basic_g
  USE mem_turb         ,  ONLY :  turb_g
  USE mem_micro        ,  ONLY :  micro_g
  USE mem_scratch      ,  ONLY :  scratch
  USE mem_leaf         ,  ONLY :  leaf_g 
  USE extras           ,  ONLY :  extra2d 
    
  IMPLICIT NONE

  ! Arguments:
  INTEGER,INTENT(IN) :: ia,iz,ja,jz,m1,m2,m3
  
  ! Local Variables:
  ! Automatic arrays:
  real :: v_dep(m2, m3, naddsc) ! dry deposition velocity
  real :: r_lsl(m2, m3, npatch) ! laminar sub-layer resistance
  real :: v_sed(m2,m3)          ! particle fall velocity

  !- If no tracers then go back
  if(naddsc == 0) return

  call dry_dep(ngrid,m1,m2,m3,npatch,ia,iz,ja,jz,jdim,dtlt &
       ,v_dep  		                                   &
       ,leaf_g(ngrid)%r_aer                                &
       ,r_lsl  		                                   &
       ,v_sed  		                                   &
       ,basic_g(ngrid)%theta	                           &
       ,basic_g(ngrid)%rv	                           &
       ,basic_g(ngrid)%pp	                           &
       ,basic_g(ngrid)%dn0	                           &
       ,basic_g(ngrid)%pi0	                           &
       ,basic_g(ngrid)%up	                           &
       ,basic_g(ngrid)%vp	                           &
       ,turb_g(ngrid)%tkep	                           &
       ,turb_g(ngrid)%sflux_t                              &
       ,micro_g(ngrid)%rcp	                           &
       ,grid_g(ngrid)%rtgt	                           &
       ,scratch%vt2da  	                                   &
       ,scratch%vt2db  	                                   &
       ,scratch%vt2dc  	                                   &
       ,scratch%vt2dd  	                                   &
       ,scratch%vt2de  	                                   &
       ,scratch%vt2df  	                                   &
       ,leaf_g(ngrid)%ustar	                           &
       ,leaf_g(ngrid)%tstar	                           &
       ,leaf_g(ngrid)%patch_area                           &
       ,leaf_g(ngrid)%leaf_class                           &
       ,leaf_g(ngrid)%patch_rough                          &
       ,extra2d(4,ngrid)%d2                                &  !save v_dep
       ,MAXGRDS		                                   &
       ,dzt,zt,nzpmax,naddsc	                           )

  return
end subroutine drydep_driver

!========================================================================

subroutine dry_dep(ngrid,m1,m2,m3,npatch,ia,iz,ja,jz,jdim,dt       &
     ,v_dep,r_aer,r_lsl,v_sed                                      &
     ,theta,rv,pp,dn0,pi0,up,vp                                    &
     ,tke,sfl,rcp,rtgt                                             &
     ,temps,prss,dens,vels,rvs,Zi                                  &
     ,ustar,tstar,patch_area,veg,Z0m	                           &
     ,v_dep_pm25, maxgrds,dzt,zt,nzpmax,naddsc                     )

  USE rconstants         ,  ONLY :  cpdryi,cpor,p00
  USE mem_scalar         ,  ONLY :  scalar_g

  IMPLICIT NONE
  INTEGER, INTENT(IN)  ::  ngrid,m1,m2,m3,npatch,ia,iz,ja,jz,jdim
  INTEGER, INTENT(IN)  ::  maxgrds,nzpmax,naddsc
  REAL, INTENT (IN)    ::  dt
  REAL, PARAMETER      ::  ubmin = 0.25
  
  REAL, INTENT (INOUT)    ::  v_dep(m2,m3,naddsc)
  REAL, INTENT (IN)       ::  r_aer(m2,m3,npatch)
  REAL, INTENT (INOUT)    ::  r_lsl(m2,m3,npatch)
  REAL, INTENT (INOUT)    ::       v_sed(m2,m3)
  REAL, INTENT (OUT)      ::  v_dep_pm25(m2,m3)		
  
  REAL, INTENT (IN)    :: theta(m1,m2,m3)
  REAL, INTENT (IN)    ::    rv(m1,m2,m3)
  REAL, INTENT (IN)    ::    pp(m1,m2,m3)
  REAL, INTENT (IN)    ::   dn0(m1,m2,m3)
  REAL, INTENT (IN)    ::   pi0(m1,m2,m3)
  REAL, INTENT (IN)    ::    up(m1,m2,m3)
  REAL, INTENT (IN)    ::    vp(m1,m2,m3)
  REAL, INTENT (IN)    ::   tke(m1,m2,m3)
  REAL, INTENT (IN)    ::   rcp(m1,m2,m3)
  
  REAL, INTENT (IN)    ::   sfl(m2,m3)
  REAL, INTENT (IN)    ::  rtgt(m2,m3)
  REAL, INTENT (INOUT) :: temps(m2,m3)   
  REAL, INTENT (INOUT) ::  prss(m2,m3)    
  REAL, INTENT (INOUT) ::  dens(m2,m3)    
  REAL, INTENT (INOUT) ::  vels(m2,m3)    
  REAL, INTENT (INOUT) ::   rvs(m2,m3)     
  REAL, INTENT (IN)    ::    Zi(m2,m3)
  
  REAL, INTENT (IN)    ::       ustar(m2,m3,npatch)
  REAL, INTENT (IN)    ::       tstar(m2,m3,npatch)
  REAL, INTENT (IN)    ::  patch_area(m2,m3,npatch)
  REAL, INTENT (IN)    ::         veg(m2,m3,npatch)
  REAL, INTENT (IN)    ::         Z0m(m2,m3,npatch)
  
  REAL, INTENT (IN)    :: dzt(nzpmax)
  REAL, INTENT (IN)    ::  zt(nzpmax)
  
  integer :: i,j,ipatch,idry_part,iscl
  real    :: ups,vps,pis
  
  v_dep = 0.0
  r_lsl = 0.0
  v_sed = 0.0

  !-aux variables

  do j = ja,jz
     do i = ia,iz
        rvs  (i,j) = rv(2,i,j)
        pis        = (pp(1,i,j) + pp(2,i,j) + pi0(1,i,j) + pi0(2,i,j))*.5 * cpdryi
        prss (i,j) = pis ** cpor * p00                                               
        dens (i,j) = ( dn0(1,i,j) + dn0(2,i,j) ) * .5
        temps(i,j) = theta(2,i,j) * pis        ! temps=theta*Exner/CP
        ups        = (up(2,i,j) + up(2,i-1,j   )) * .5
        vps        = (vp(2,i,j) + vp(2,i,j-jdim)) * .5
        vels (i,j) = max(ubmin,sqrt(ups** 2 + vps** 2))
     enddo
  enddo
  
  call define_PBL_height(m1,m2,m3,npatch,ia,iz,ja,jz,zt,tke,sfl,rcp,rtgt,zi)

  !print*,'!- loop  gases/aerossois'
  do iscl = 1,naddsc         

     idry_part = 0
     if(iscl == 3) idry_part = 1     ! whit dry deposition (pm25)

     if(idry_part == 0) then
        !--  dry deposition for gases
        call dry_dep_gaseous(m1,m2,m3,naddsc,npatch,ia,iz,ja,jz  &
             ,V_dep(:,:,iscl),patch_area)
     else
        !--  dry deposition for particles
        call dry_dep_particles(iscl,m1,m2,m3,naddsc,npatch,ia,iz,ja,jz  &
             ,V_dep(:,:,iscl)                             &
             ,r_aer,r_lsl,v_sed                           &
             ,temps,dens,vels,rvs,Zi 		     &
             ,ustar,tstar,patch_area,veg,Z0m	             )
      
        !for output - v_edp of PM25
        v_dep_pm25(:,:)=V_dep(:,:,3) 
      
        !----------  
      
     endif

     !-apply dry deposition on the tracers concentration and get the 
     ! deposited mass on surface
     ! **(JP)** Pass to subroutine first pointer element,copy are not necessary
     call apply_drydep(m1,m2,m3,ia,iz,ja,jz,V_dep(:,:,iscl), &
          scalar_g(iscl,ngrid)%drydep,                       &
          scalar_g(iscl,ngrid)%sclp,                         &
          scalar_g(iscl,ngrid)%sclt,                         &
          dens,rtgt,dzt,dt                                   )
   
  enddo

  return
end subroutine dry_dep

!========================================================================

subroutine dry_dep_particles(iscl,m1,m2,m3,naddsc,npatch,ia,iz,ja,jz &
     ,V_dep,r_aer,r_lsl,v_sed                                        &
     ,temps,dens,vels,rvs,Zi                                         &
     ,ustar,tstar                                                    &
     ,patch_area,veg,Z0m                                             )
  use leaf_coms, only : min_patch_area
  implicit none
  integer :: m1,m2,m3,naddsc,npatch,ia,iz,ja,jz,i,j,ipatch,iscl
  real    :: vdtmp
  real, dimension(m2,m3)        :: temps,dens,vels,rvs,Zi
  real, dimension(m2,m3,npatch) :: ustar,tstar,patch_area,veg,Z0m
  
  real, dimension(m2,m3)     :: V_dep
  real, dimension(m2,m3,npatch) :: r_aer,r_lsl
  real, dimension(m2,m3)     :: v_sed

  !- sedimentation  parametrization 

  call sedim_particles(m2,m3,npatch,ia,iz,ja,jz,temps,dens,v_sed)

  !- laminar sub-layer resistance
  call lsl_particles(m2,m3,npatch,ia,iz,ja,jz                  &
       ,temps,dens,vels,rvs,Zi,ustar,tstar,patch_area,veg,Z0m  &
       ,v_sed,r_lsl)

  !print*,'!- Particles Deposition velocity (cm/s)'

  do j = ja,jz
     do i = ia,iz
        do ipatch = 1,npatch

           if (patch_area(i,j,ipatch) >= min_patch_area) then

              vdtmp = v_sed(i,j) + 1./(r_aer(i,j,ipatch) + r_lsl(i,j,ipatch) +&
                   r_aer(i,j,ipatch) * r_lsl(i,j,ipatch) * v_sed(i,j))
              
              V_dep(i,j) = V_dep(i,j) + patch_area(i,j,ipatch)*vdtmp

           endif
        enddo
     enddo
  enddo

  return
end subroutine dry_dep_particles

!========================================================================

subroutine sedim_particles(m2,m3,npatch,ia,iz,ja,jz,temps,dens,v_sed)

  implicit none
  
  integer :: m2,m3,naddsc,npatch,ia,iz,ja,jz,i,j,ipatch
  real, dimension(m2,m3) :: temps,dens
  REAL,PARAMETER :: ASP = 1.257          ! 1.249
  REAL,PARAMETER :: BSP = 0.4            ! 0.42
  REAL,PARAMETER :: CSP = 1.1            ! 0.87 

  !
  !For small particles diameter (<20 micrometer) (low Reynolds number)
  !the fall velocity is given by (Seinfeld & Pandis, Jacobson)
  ! V_s = 2 r**2 (rho_p - rho_air) * g * Gi / 9 n_air
  ! where:
  ! r == radius of particle (m)
  ! rho_p,a = density of particle, air (kg/m^3)
  ! g = gravity accel (9.8 m/s^2)
  ! Gi = Kn(A' + B' + C' exp(-C'/Kn) , Kn = Knudsen number
  ! n_air = dynamic viscosity of air
  !

  real, dimension(m2,m3)    :: v_sed

  real :: A,B,C,Kn,n_air,Gi,v_air,mfp,nu,rey,part_dens,part_radius

  !- constantes
  !parameter (g=9.80622, pi=3.141593)
  !parameter (ASP=1.257,BSP=0.4,CSP=1.1) ! Seinfeld & Pandis 
  
  !parameter (kB = 1.3807e-23)     ! const Boltzmann - kg m^2 s^-2 K^-1 molecule^-1
  !parameter (M_AVEG = 4.8096e-26) ! average mass of one molecure - kg  molecule^-1
  
  !- particle properties ! ver com a Karla:
  !kml **part_dens ja esta definido na rotina setupaer com o nome rhoelem, mas a
  !kml **ordem de grandeza e essa mesmo.
  !kml **para o part_radius idealmente deveria ser usado o r0(ix,iy,ik) definido
  !kml **na initaer

  part_radius = 1.95e-5 * 1.e-2       ! effective particle radius (meter)
  part_dens   = 1.35 * 1.e+3    ! particle density kg/m^3

  do j = ja,jz
     do i = ia,iz

        !- several particle/environment properties

        !- mean speed of air molecules (m/s)
        !  v_air = sqrt( 8. * kB   * temps(i,j) / (pi * M_AVEG) )
        v_air = sqrt( 7.3102e+2 * temps(i,j)                        )
        
        !-dynamic viscosity of air (kg m^-1 s^-1)
        !  n_air = 1.8325e-5*(416.16/(temps(i,j)+120.))*(temps(i,j)/296.16)**1.5
        !optimized version
        n_air = 1.8325e-5*(416.16/(temps(i,j)+120.))*(temps(i,j)/296.16) * & 
             sqrt(temps(i,j)/296.16)   
        !-- kinematic viscosity of air ()
        nu = n_air/dens(i,j)
        !- mean free path of an air molecule (m)
        mfp = 2.* n_air /(dens(i,j)*v_air)
        
        !- Knudsen number
        Kn = mfp/part_radius

        !- Gi determination
        Gi = 1. + Kn*( ASP + BSP*exp(-CSP/Kn) )

        !- particle sedimentation velocity (m/s)
        !  v_sed(i,j) = (2./9.)*g*part_radius**2  *(part_dens-dens(i,j))*Gi/n_air
        v_sed(i,j) =     2.18    *part_radius**2. *(part_dens-dens(i,j))*Gi/n_air

        ! Reynolds number
        !   rey = 2*part_radius*v_sed(i,j)/nu
        
        !   print*,'Sedimentation velocity cm/s'   
        !   print*,'VSED',i,j,v_sed(i,j)*100.
        !   print*,'VSED',n_air,Gi,temps(i,j)

     enddo
  enddo

  return
end subroutine sedim_particles

!========================================================================

subroutine lsl_particles(m2,m3,npatch,ia,iz,ja,jz           &
     ,temps,dens,vels,rvs,Zi,ustar,tstar,patch_area,veg,Z0m &
     ,v_sed,r_lsl)
  use rconstants, only: t00, vonk,cpdry,pi1,grav,boltzmann
  use leaf_coms, only : min_patch_area
  implicit none
  REAL,PARAMETER :: ASP = 1.257          ! 1.249
  REAL,PARAMETER :: BSP = 0.4            ! 0.42
  REAL,PARAMETER :: CSP = 1.1            ! 0.87 
  
  REAL,PARAMETER :: M_AVEG = 4.8096e-26  ! average mass of one molecure - kg  molecule^-1
  REAL,PARAMETER :: em23 = -2./3., ep23 = +2./3., ep13=1./3.    ! exponents 2./3. 1/3
 
  integer :: m2,m3,npatch,ia,iz,ja,jz,i,j,ipatch
  real, dimension(m2,m3) :: temps,dens,vels,rvs,Zi
  real, dimension(m2,m3,npatch) :: ustar,tstar,patch_area,veg,Z0m
  real :: wptp,wstar
  !
  real, dimension(m2,m3)    :: v_sed       ! particle sedimentation velocity   (m/s)
  real, dimension(m2,m3,npatch) :: r_lsl   ! resistance to molecular diffusion (s/m)
  
  real :: Kn,n_air,part_dens,part_radius,Gi,v_air &
       ,mfp,D,nu,Sc,St,Kd,Dh,Z0h,Pr    

  real :: limite !ALF
  
  limite = -30
  
  part_radius = 1.95e-5 * 1.e-2       ! effective particle radius (meter)
  part_dens   = 1.35 * 1.e+3    ! particle density kg/m^3

  !        print*,'DENTRO DA LSL_PARTICLES'

  do j = ja,jz
     do i = ia,iz

        !- several particle/environment properties
        
        !- mean speed of air molecules (m/s)
        !  v_air = sqrt( 8. * boltzmann   * temps(i,j) / (pi1 * M_AVEG) )
        v_air = sqrt( 7.3102e+2 * temps(i,j)                        )
        
        !-dynamic viscosity of air (kg m^-1 s^-1)
        !  n_air = 1.8325e-5*(416.16/(temps(i,j)+120.))*(temps(i,j)/296.16)**1.5
        !optimized version
        n_air = 1.8325e-5*(416.16/(temps(i,j)+120.))*(temps(i,j)/296.16) * & 
             sqrt(temps(i,j)/296.16)
        
        !- mean free path of an air molecule (m)
        mfp = 2.* n_air /(dens(i,j)*v_air)
        
        !- Knudsen number
        Kn = mfp/part_radius
        
        !- Gi determination
        Gi = 1. + Kn*( ASP + BSP*exp(-CSP/Kn) )
        
        !- Schmidt number determination (Jacobson)
        !-- molecular diffusivity (Brownian diffusivity coeficient)   
        !XXXXXXXXXXXXXXXXXXXXXXXXXXXX
        D =  boltzmann * temps(i,j) * Gi / (6.*pi1*part_radius*n_air)
        
        !x    D =  (boltzmann/M_AVEG) * temps(i,j) * Gi / (6.*pi1*part_radius*n_air)
        
        !-  kinematic viscosity of air ()
        nu = n_air/dens(i,j)
        !-  Schmidt number
        Sc = nu/D
        
        !- laminar sub-layer resistance for particles (s m^-1)
        do ipatch=1,npatch

           !-  Stokes number determination (Slinn, 1980)
           !-  St =     ustar^2 * V_sedim /( g * kinematic viscosity of air)   
           St = max(ustar(i,j,ipatch)**2. * v_sed(i,j) / (grav * nu), 0.01) 
    
           !-  laminar sub-layer resistance for particles (s m^-1) 

           if(ipatch == 1) then !- water

              !-from Slinn 1980 (Atmos. Env.)

              r_lsl(i,j,ipatch) = ( vonK * vels(i,j)/ustar(i,j,ipatch)**2. ) /&
                   ( 1./sqrt(Sc) + 10.**max((-3./St),limite) )
			  
              !print*,'LSL-OCEAN',i,j,ipatch,nint(veg(i,j,ipatch)),patch_area(i,j,ipatch)
              !PRINT*,ustar(i,j,ipatch),r_lsl(i,j,ipatch)

           else

              if(patch_area(i,j,ipatch) >= min_patch_area) then
                 
                 !-- for smooth land surfaces !- bare ground (3), desert(3) or ice (2)
                 !-- Seinfeld & Pandis (1998)
                 if(nint(veg(i,j,ipatch)) == 3 ) then 
            
                    r_lsl(i,j,ipatch) = 1./(ustar(i,j,ipatch) * &
                         (Sc**em23 + 10.**max((-3./St),limite)))
         
                 else

                    !- expression from Jacobson(1999)
!
                    !- thermal conductivity of dry air (Kd)
                    Kd = 0.023807 + 7.1128e-5*(temps(i,j) - t00) !- Eq.(2.3) 
                    !- Prandt number
                    Pr =  n_air*cpdry*(1.+0.859*rvs(i,j))/Kd           !- Eq.(17.32)  
                    
                    !- energy moisture roughness lengths (Z0h)                
                    !- Eq.(8.10)
                    !-- molecular thermal diffusion coeff. (m^2 s^-1)
                    Dh=Kd/(dens(i,j)*cpdry)
                    !- Z0h	   
                    Z0h=Dh/(vonK*ustar(i,j,ipatch))

                    !- lsl resistance according Jacobson Eq. (20.14)
                    !         r_lsl(i,j,ipatch) =  log( Z0m(i,j,ipatch)/Z0h )* ( (Sc/Pr)**ep23 ) &
                    !	                     /  ( vonK*ustar(i,j,ipatch) )
                    
                    !       print*,'================'
                    !       print*,'LSL-JAC',i,j,ipatch,nint(veg(i,j,ipatch)),patch_area(i,j,ipatch)
                    !	PRINT*,Sc,Pr,Z0m(i,j,ipatch),Z0h
                    !	PRINT*,ustar(i,j,ipatch),r_lsl(i,j,ipatch)

                    !---------	             
                    !- lsl resistance according :
                    !- from Wesely et al. (1985) following Slinn (1982)
                    !- also Binkowski & Shankar, JGR, 100,D12,26191-26209, 1995
                    !- Rb= (u* (1+ (w*/u*)^2)) (Sc^2./3 + 10^(-3/St))    

                    wptp  = - ustar(i,j,ipatch) * tstar(i,j,ipatch) ! sensible heat flux
                    wstar = ( max (0., grav* Zi(i,j)*  wptp/temps(i,j) ) )**ep13  
                    !	 
                    r_lsl(i,j,ipatch) = 1./(                         &
                         ustar(i,j,ipatch)*                          &
                         (1. + 0.24*(wstar/ustar(i,j,ipatch))**2.)*  &
                         ( (Sc**em23 + 10.**max((-3./St),limite) ) ) )

                    !print*,'================'
                    !print*,'LSL_WE',i,j,ipatch,nint(veg(i,j,ipatch)),patch_area(i,j,ipatch)
                    !PRINT*,wstar,ustar(i,j,ipatch),r_lsl(i,j,ipatch)

                 endif
              endif
           endif

           !print*,'laminar sub-layer resistance'
           !print*,'LSL_J',i,j,ipatch,r_lsl(i,j,ipatch)
           !print*,'ZOM',Z0m(i,j,ipatch)
           !print*,'Z0H',Z0h
           !print*,'SC PR',Sc,Pr
           !print*,'u*',ustar(i,j,ipatch)

        enddo

     enddo
  enddo

  return
end subroutine lsl_particles

!========================================================================

subroutine dry_dep_gaseous(m1,m2,m3,naddsc,npatch,ia,iz,ja,jz &
     ,V_dep,patch_area)
  use leaf_coms, only : min_patch_area
  implicit none
  integer :: m1,m2,m3,naddsc,npatch,ia,iz,ja,jz,i,j,ipatch
  real    :: vdtmp
  real, dimension(m2,m3,npatch) :: patch_area

  real, dimension(m2,m3)     :: V_dep
  real, dimension (2) :: V_dep_CO
  data (V_dep_CO(i),i=1,2)  &
       /0.        &  ! m/s (ocean)
       ,0.0003 /     !     (continent)      

  !print*,'!- Gaseous Deposition Velocity (m/s)'

  do j = ja,jz
     do i = ia,iz
       
        ! ipatch = 1
        V_dep(i,j)=  V_dep_CO(1)*patch_area(i,j,1)
        
        do ipatch = 2,npatch

           if (patch_area(i,j,ipatch) >= min_patch_area) &
                V_dep(i,j) = V_dep(i,j) + patch_area(i,j,ipatch)*V_dep_CO(2)

        enddo
     enddo
  enddo

return
end subroutine dry_dep_gaseous

!========================================================================

subroutine apply_drydep(m1,m2,m3,ia,iz,ja,jz,V_dep,M_dep,sclp,sclt &
     ,dens,rtgt,dzt,dt)

  implicit none
  integer :: m1,m2,m3,ia,iz,ja,jz,i,j

  real, dimension(m2,m3)    :: V_dep &! dry deposition velocity (m/s)
       ,M_dep  ! accumulated mass on surface due dry dep 
  ! process (kg m^-2)
  real, dimension(m1,m2,m3) :: sclp & ! tracer concentration (kg/kg)
       ,sclt   ! tendency (kg[tracer] kg[air]^-1 s^-1)
  real, dimension(m2,m3) :: rtgt,dens  
  real, dimension(m1) :: dzt  
  real :: dt,dz,tend_drydep			    

  do j=ja,jz
     do i=ia,iz

        !-1st vertical layer thickness
        dz = rtgt(i,j)/dzt(2) ! dzt=1/(z(k)-z(k-1))
        !- tendency to dry deposition process
        ! tend_drydep  = - V_dep(i,j)*sclp(2,i,j)/(dz*dens(i,j))
        tend_drydep  = - V_dep(i,j)*sclp(2,i,j)/(dz)
        !- update total tendency kg[tracer]/kg[air]/s
        sclt(2,i,j) = sclt(2,i,j) + tend_drydep
        !- accumulate the surface deposited mass of the tracer by this process
        ! M_dep(i,j) = M_dep(i,j) + dt*tend_drydep*dens(i,j) ! kg[tracer] m^-2
        M_dep(i,j) = M_dep(i,j) - dt*tend_drydep*dens(i,j)*dz ! kg[tracer] m^-2

     enddo
  enddo

return
end subroutine apply_drydep

!========================================================================

subroutine define_PBL_height(m1,m2,m3,npatch,ia,iz,ja,jz,zt,tke,shf,rcp,rtgt,Zi)

  implicit none
  integer :: m1,m2,m3,npatch,ia,iz,ja,jz,i,j,k
  real, dimension(m1)    ::    zt
  real, dimension(m2,m3) ::    shf,rtgt,Zi
  real, dimension(m1,m2,m3) :: tke,rcp
  real :: pblht
  REAL,PARAMETER :: tkethrsh=0.001  !   tke threshold for PBL height in m2/s2
                                    !   tkmin    = 5.e-4   minimum TKE in RAMS 
  REAL,PARAMETER :: rcmin=1.e-4     !   min liq water = 0.1 g/kg

  do j=ja,jz
     do i=ia,iz

        Zi(i,j) = 0.
        !- convective layer
        if(shf(i,j) >= 1.e-8) then 
      
           pblht=0.
           do k=2,m1-1
              pblht=zt(k)*rtgt(i,j) 
              !if(i.ge.10.and.i.le.25.and.j.ge.13.and.j.le.25) &
              !     print*,'i,j,k,z,pbl=',i,j,k,ztn(k,ngrd),pblht
              if( rcp(k,i,j) .gt. rcmin     ) goto 10 ! dry convective layer
              if( tke(k,i,j) .le. tkethrsh  ) goto 10 
           enddo
10         continue

           Zi(i,j)=pblht
        endif
        !print*,'PBLh',Zi(i,j)

     enddo
  enddo

  return
end subroutine define_PBL_height

