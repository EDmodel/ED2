 !############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

!     *****************************************************************

subroutine tkescl(m1,m2,m3,m4,ia,iz,ja,jz &
          ,tkep,tket,epsp &
          ,vt3da,vt3dc,vt3dh,vt3di,vt3dj &
          ,scr1,scr2,rtgt &
          ,vt3dd,vt3de,dxt &
          ,ustar,patch_area,flpw,dn0)     
!_STC  +-------------------------------------------------------------------+
!_STC  \  Authors: S. Trini Castelli, E. Ferrero                           \
!_STC  \  First version:        1998                                       \
!_STC  \  Last updated version: 2000                                       \
!_STC  \  Version for RAMS4.3 : 2001                                       \
!_STC  +-------------------------------------------------------------------+
!_STC  \  This routine calculates the mixing coefficients based on the     \
!_STC  \    prognostic TKE and diagnostic length scale, according to a     \
!_STC  \    standard E-l closure model                                     \
!_STC  +-------------------------------------------------------------------+

use mem_grid, only: &
              zt, &       !INTENT(IN)
              zm, &       !INTENT(IN)
              dzm, &      !INTENT(IN)
              npatch, &   !INTENT(IN)
              dtlt        !INTENT(IN)
        
use mem_scratch, only: &
                 vctr30, &   !INTENT(OUT)
                 vctr31, &   !INTENT(OUT)
                 vctr1, &    !INTENT(OUT)
                 vctr33, &   !INTENT(OUT)
                 vctr32, &   !INTENT(OUT)
                 vctr25, &   !INTENT(OUT)
                 vctr27, &   !INTENT(OUT)
                 vctr28      !INTENT(OUT)                
                 
use ke_coms, only: &
             c_eps, &       !INTENT(IN)
             iopzl, &       !INTENT(IN)
             al0_const, &   !INTENT(IN)
             alf_tht, &     !INTENT(IN)
             c_mi, &        !INTENT(IN)
             alf_tke        !INTENT(IN)
             
             
use rconstants, only: &
                vonk, &        !INTENT(IN)
                grav           !INTENT(IN)

implicit none

integer, intent(in) :: &
                       m1, &
                       m2, &
                       m3, &
                       m4, &
                       ia, &
                       iz, &
                       ja, &
                       jz  

real, intent(in),    dimension(m1, m2, m3) :: &
                                              tkep,  &
                                              vt3da, &
                                              vt3dc, &
                                              vt3dj, &
                                              vt3dd, &
                                              vt3de, &
                                              dn0
                                              

real, intent(inout), dimension(m1, m2, m3) :: & 
                                              tket,  &
                                              vt3dh, &
                                              scr2 
                                              

real, intent(out),   dimension(m1, m2, m3) :: &
                                              epsp,  &
                                              vt3di, &
                                              scr1
                                              
                                            
    
real, intent(in),    dimension(m2,m3)      :: &
                                              rtgt, &
                                              dxt

real, intent(in),    dimension(m2,m3,m4)   :: &
                                              ustar, &
                                              patch_area

real, intent(in), dimension(m2,m3)      :: flpw

real :: coef2
real :: dudz_1,dudz_2,dvdz_1,dvdz_2,dvelmin_1,dvelmin_2 &
        ,psizil_1q,psizil_2q,psizil_q
real :: al0,al0_zil,scl,sumtkz,sumtk,tkep_k2,dzloc,dpsi2dz,scl_max
real, external :: ssum
integer :: k,i,j,np,k2

coef2= 1./(c_eps**(2./3.))

!_STC.................................................................
!STC
!  - C_EPS,C_MI,ALF_TKE,ALF_EPS,ALF_THT assigned in rconstants.h
!  - aymptotic mixing length option assigned in rconstants.h
!  - Constant value (Ying, 1992) for the asymptotic mixing length assigned 
!  in rconstants.h
!STC
!_STC.................................................................

      do j=ja,jz
         do i=ia,iz
            k2=nint(flpw(i,j))

!_STC.........................................................
!_STC   Calculation of the asymptotic lenght scale: OPTIONS 
!_STC.........................................................
!
!_STC-----> MELLOR-YAMADA 1974: pre-calculation 
            if(iopzl.eq.2) then	

               do k=k2,m1-1
                  vctr30(k) = sqrt(2.0 * tkep(k,i,j))
                  vctr31(k)=.75*vctr30(k)/sqrt(max(1.e-20,vt3dj(k,i,j)))
                  vctr1(k)=(zt(k)-zm(k2-1))*rtgt(i,j)
                  dzloc=(zm(k)-zm(k-1))*rtgt(i,j)
                  vctr33(k)=vctr30(k)*dzloc
                  vctr32(k)=vctr33(k)*vctr1(k)
               enddo
!
               sumtkz=ssum(m1-k2,vctr32(k2),1)
               sumtk =ssum(m1-k2,vctr33(k2),1)
               al0=.1*sumtkz/sumtk

            endif


            do k=k2,m1-1


!_STC....................................
!_STC-----> Ying's value
!_STC....................................

               if(iopzl == 1) then    ! YING's length scale
                  scl=vonk*vctr1(k)/(1.+vonk*vctr1(k)/al0_const) 

!_STC....................................
!_STC-----> Mellor-Yamada 1974 formulation
!_STC....................................

               elseif(iopzl == 2) then
                  scl=min(vonk*vctr1(k)/(1.+vonk*vctr1(k)/al0),vctr31(k))

!_STC................................................
!_STC-----> Zilitinkevich-Laikhtman 1965 formulation
!_STC................................................
!_5Mar02     trasformo la formulazione di Zilitinkevich da funzione di
!_5Mar02     psi a funzione di psi**2

               elseif(iopzl == 3) then	! opzione Zilitinkevich

!
!... calcolo della derivata in z e del valore della lunghezza l0


                  dudz_1=vt3dd(k,i,j)
                  dudz_2=vt3dd(k+1,i,j)
                  dvdz_1=vt3de(k,i,j)
                  dvdz_2=vt3de(k+1,i,j)
   
!.... valore minimo per le derivate in z di u e v
   
                  dvelmin_1=0.001/dzm(k)
                  dvelmin_2=0.001/dzm(k+1)
   
                  if(abs(dudz_1)<dvelmin_1 &
                            .and. dudz_1>=0.) dudz_1=dvelmin_1
                  if(abs(dudz_1)<dvelmin_1 &
                            .and. dudz_1<0.) dudz_1=-dvelmin_1
                  if(abs(dvdz_1)<dvelmin_1 &
                            .and. dvdz_1>=0.) dvdz_1=dvelmin_1
                  if(abs(dvdz_1)<dvelmin_1 &
                            .and. dvdz_1<0.) dvdz_1=-dvelmin_1
                  if(abs(dudz_2)<dvelmin_2 &
                            .and. dudz_2>=0.) dudz_2=dvelmin_2
                  if(abs(dudz_2)<dvelmin_2 &
                            .and. dudz_2<0.) dudz_2=-dvelmin_2
                  if(abs(dvdz_2)<dvelmin_2 &
                            .and. dvdz_2>=0.) dvdz_2=dvelmin_2
                  if(abs(dvdz_2)<dvelmin_2 &
                            .and. dvdz_2<0.) dvdz_2=-dvelmin_2
   
                  psizil_1q=dudz_1*dudz_1 + dvdz_1*dvdz_1  &
                                - alf_tht*vt3dj(k,i,j)
                  psizil_2q=dudz_2*dudz_2 + dvdz_2*dvdz_2  &
                                - alf_tht*vt3dj(k+1,i,j)

                  psizil_q=0.5*(psizil_1q+psizil_2q)

! 'special' unlucky case control !
                  if(psizil_2q == psizil_1q) then
                   psizil_1q=2*dvelmin_1*dvelmin_1 - alf_tht*vt3dj(k,i,j)
                   psizil_2q=2*dvelmin_2*dvelmin_2 - alf_tht*vt3dj(k+1,i,j)
                  endif
   
                  dpsi2dz=(psizil_2q-psizil_1q)*dzm(k)/rtgt(i,j)
!
                  if(dpsi2dz.eq.0.) stop 'dpsi2dz = 0. chissa perche!'

                  if(dpsi2dz.ne.0.) then
 
                  al0_zil = -4.* vonk *(psizil_q/dpsi2dz)    
                  if(al0_zil.lt.0.) al0_zil=abs(al0_zil)
               endif
   
               scl=vonk*vctr1(k)/(1. + vonk*vctr1(k)/al0_zil ) ! Zily

!_8Mar02   maximum limit for SCL: copied from Deardoff scheme

               scl_max = 1.5*(rtgt(i,j)*(zm(k) - zm(k-1)))**(1./3.) &
                   /dxt(i,j) ** (2./3.)
               if (scl .gt. scl_max) scl = scl_max

   
            endif     !fine variante Zilitinkevic
!----------------------------------------------------------------------------
!_STC................................................
!_STC	T.K.E. equation terms
!_STC................................................
!
!_STC................................................
!_STC......vertical diffusion coefficient assignment
!_STC................................................
           scr1(k,i,j)=c_mi*scl*sqrt(tkep(k,i,j))    ! Km

!_STC................................................
!_STC   ......production term
!_STC................................................

           vctr25(k)=scr1(k,i,j)*vt3dh(k,i,j)-0.333333*tkep(k,i,j)  &
                *(vt3da(k,i,j)+vt3dc(k,i,j)+scr2(k,i,j))

!_STC................................................
!_STC	......buoyancy term
!_STC................................................

           vctr27(k)=-(scr1(k,i,j)*alf_tht)*vt3dj(k,i,j)


!_STC................................................
!_STC	......dissipation term (diagnostic)
!_STC................................................

           vctr28(k)=c_eps*tkep(k,i,j)*sqrt(max(1.e-20,tkep(k,i,j)))/scl
             
           epsp(k,i,j)= vctr28(k)

!_STC................................................
!_STC	T.K.E. equation
!_STC................................................
           tket(k,i,j)=tket(k,i,j)+vctr25(k)+vctr27(k)-vctr28(k)

!_STC....................................................................
!_STC   Assignment of diffusion coefficients to the correspondent arrays
!_STC	Km(oriz) = Km(vert) : isotropy
!_STC....................................................................
!
           scr1(k,i,j)=scr1(k,i,j)*dn0(k,i,j)   ! account for density
           scr2(k,i,j)=scr1(k,i,j)   ! used in DIFFVEL as horiz. diff. coef.
           vt3dh(k,i,j)=scr1(k,i,j)*alf_tht   ! HEAT diffusion coefficient
           vt3di(k,i,j)=scr1(k,i,j)*alf_tke   ! TKE diffusion coefficient Ke

        enddo
        
! Bottom boundary condition        

      tkep_k2=0.
      do np=1,npatch
         tkep_k2=tkep_k2+(ustar(i,j,np)**2)*coef2  &
                         *patch_area(i,j,np)
      enddo
      tket(k2,i,j)=(tkep_k2-tkep(k2,i,j))/dtlt
   
     enddo
  enddo

return
end subroutine tkescl



!**********************************************************************

subroutine tkeeps(m1,m2,m3,m4,ia,iz,ja,jz  &
          ,tkep,tket,epsp,epst &
          ,vt3da,vt3dc,vt3dh,vt3di,vt3dj,scr1,scr2 &
          ,rtgt,ustar,patch_area,flpw,dn0)
!_STC  +-------------------------------------------------------------------+
!_STC  \  Authors: S. Trini Castelli, E. Ferrero                           \
!_STC  \  First version:        1998                                       \
!_STC  \  Last updated version: 2000                                       \
!_STC  \  Version for RAMS4.3 : 2001                                       \
!_STC  +-------------------------------------------------------------------+
!_STC  \  This routine calculates the mixing coefficients based on the     \
!_STC  \  prognostic TKE and prognostic TKE dissipation, according         \
!_STC  \  for instance, to the RODI (1980)                                 \
!_STC  \  (see also DETERING-ETLING's paper (BLM, 1985, n. 33) )           \
!_STC  +-------------------------------------------------------------------+

use mem_grid, only: &
              zt, &           !INTENT(IN)
              zm, &           !INTENT(IN)
              ngrid, &        !INTENT(IN)
              ngrids, &       !INTENT(IN)
              dtlt, &         !INTENT(IN)
              time            !INTENT(IN)
              
use mem_scratch, only: &
                 vctr25, &      !INTENT(OUT)
                 vctr27, &      !INTENT(OUT)
                 vctr28         !INTENT(OUT)

use ke_coms, only: &
             coef_km, &         !INTENT(IN)
             alf_tht, &         !INTENT(IN)
             c1_eps, &          !INTENT(IN)
             c2_eps, &          !INTENT(IN)
             alf_tke            !INTENT(IN)

use rconstants, only: &
                vonk            !INTENT(IN)

implicit none

integer, intent(in)  ::  &
                         m1, &
                         m2, &
                         m3, &
                         m4, &
                         ia, &
                         iz, &
                         ja, &
                         jz

real, intent(in),    dimension(m1,m2,m3) :: &
                                            tkep,  &
                                            epsp,  &
                                            vt3da, &
                                            vt3dc, &
                                            vt3dj, &
                                            dn0
                                            

real, intent(inout), dimension(m1,m2,m3) :: &
                                            tket,  &
                                            epst,  &
                                            vt3dh, &
                                            scr2
                                            
real, intent(out), dimension(m1,m2,m3)   :: &
                                            vt3di, &
                                            scr1


real, intent(in), dimension(m2,m3)    :: rtgt

real, intent(in), dimension(m2,m3,m4) :: &
                                         ustar, &
                                         patch_area

real, intent(in), dimension(m2,m3) :: flpw

integer :: k,i,j,np,k2
real :: epsp_k2,tkep_k2
real :: coef1,coef_km_sqr

coef_km_sqr = sqrt(coef_km)
!
!_STC................................................................
!     Empirical constants and coefficients for the E-eps closure
!           C_MI,C1_EPS,C2_EPS,ALF_TKE,ALF_EPS,ALF_THT
!     assigned in ke_coms
!_STC................................................................


do j=ja,jz
   do i=ia,iz
      k2=nint(flpw(i,j))
      coef1=1./(vonk*(zt(k2)-zm(k2-1))*rtgt(i,j))

      do k=k2,m1-1         !! loop extended to level 2

!_STC................................................
!_STC	......vertical diffusion coefficient
!_STC................................................
!
         scr1(k,i,j)=coef_km*(tkep(k,i,j)**2.)/epsp(k,i,j)
!
!_STC................................................
!_STC   wind shear production term: P_s
!_STC................................................
!
         vctr25(k)=scr1(k,i,j)*vt3dh(k,i,j)-0.333333*tkep(k,i,j)  &
                  *(vt3da(k,i,j)+vt3dc(k,i,j)+scr2(k,i,j))

!_STC................................................
!_STC  buoyancy term
!_STC    VT3DJ(K,I,J) = Brunt-Vaisala frequency from subroutine BRUVAIS
!_STC................................................

         vctr27(k)=-(scr1(k,i,j)*alf_tht)*vt3dj(k,i,j)

!_STC................................................
!_STC  dissipation term (from prognostic equation at a former time step)
!_STC................................................

         vctr28(k)= epsp(k,i,j)

!_STC................................................
!_STC     equation for the TKE 
!_STC................................................

         tket(k,i,j)=tket(k,i,j)+vctr25(k)+vctr27(k)-vctr28(k)
!_STC.....................................................................
!_STC    equation for the TKE dissipation: RODI
!_STC
!_STC in the "vertical layers" the lateral component of the buoyancy 
!_STC is orthogonal to gravity, so that there is not latera contribution
!_STC from the buoyancy and the correction with c3_eps vanishes. (see Rodi)
!_STC.....................................................................

         epst(k,i,j)= epst(k,i,j) &
           + c1_eps*epsp(k,i,j)*(vctr25(k)+vctr27(k))/tkep(k,i,j) &
           - c2_eps*epsp(k,i,j)*epsp(k,i,j)/tkep(k,i,j)
!
!
!_STC.....................................................................
!_STC  Array assignment of the diffusion coefficients for the
!_STC  scalar diffusion. In SCR2 is assigned the diffusion coefficient
!_STC  for the dissipation tkeps.
!_STC.....................................................................

         scr1(k,i,j)=scr1(k,i,j)*dn0(k,i,j)   ! account for density
         scr2(k,i,j)=scr1(k,i,j)            ! used in DIFFVEL as horizontal
                                            !      diffusion coefficient
         vt3dh(k,i,j)=scr1(k,i,j)*alf_tht   ! HEAT diffusion coefficient
         vt3di(k,i,j)=scr1(k,i,j)*alf_tke   ! TKE diffusion coefficient Ke

if(i==25.and.j==20.and.ngrid==ngrids.and.k<10) print '(a,i3,12e14.5)','tkeeps:',k  &
,tkep(k,i,j),epsp(k,i,j),scr1(k,i,j),tket(k,i,j),vctr25(k),vctr27(k),-vctr28(k)

      enddo
      
! bottom boundary conditions
      epsp_k2=0.
      tkep_k2=0.
      do np=1,m4
         tkep_k2=tkep_k2+(ustar(i,j,np)**2)*coef_km_sqr  &
             *patch_area(i,j,np)
!altern       tkep_k2=tkep_k2+(ustar(i,j,np)**2.)*patch_area(i,j,np)/C_ust_tk
         epsp_k2=epsp_k2+(ustar(i,j,np)**3)*coef1  &
                *patch_area(i,j,np)
!alternative  epsp_k2=c_eps*(tkep_k2**1.5)*coef1
      enddo
      tket(k2,i,j)=(tkep_k2-tkep(k2,i,j))/dtlt
      epst(k2,i,j)=(epsp_k2-epsp(k2,i,j))/dtlt

   enddo
enddo

if(ngrid==ngrids)print*,'-------------------------------',time,coef_km

return
end subroutine tkeeps



!     ******************************************************************

subroutine mxtked(m1,m2,m3,ia,iz,ja,jz,ibcon,jd  &
   ,tkep,tket,up,vp,wp,rtp,rv,theta,vt3da  &
   ,vt3dc,vt3dh,vt3dj,scr1,scr2,sflux_u,sflux_v,sflux_w,sflux_t,dxt,rtgt,flpw)
!  +-------------------------------------------------------------------+
!  \  this routine calculates the mixing coefficients based on the     \
!  \    prognostic tke and diagnostic length scale, according to the   \
!  \    deardorff scheme.                                              \
!  +-------------------------------------------------------------------+

use mem_grid, only: &
              zm                 !INTENT(IN)

use mem_scratch, only: &
                 vctr25, &       !INTENT(OUT)
                 vctr26, &       !INTENT(OUT)
                 vctr27, &       !INTENT(OUT)
                 vctr28          !INTENT(OUT)

!use rconstants, only:

implicit none

integer, intent(in) :: m1,m2,m3,ia,iz,ja,jz,ibcon,jd

real, intent(in),    dimension(m1, m2, m3) :: &
                                              tkep, &
                                              up, &     !not used in this subroutine
                                              vp, &     !not used in this subroutine
                                              wp, &     !not used in this subroutine
                                              rtp, &    !not used in this subroutine
                                              rv, &     !not used in this subroutine
                                              theta, &  !not used in this subroutine
                                              vt3da, & 
                                              vt3dc, &
                                              vt3dj



real, intent(inout), dimension(m1, m2, m3) :: &
                                              tket,  &
                                              vt3dh, &
                                              scr2
                                              
real, intent(out), dimension(m1, m2, m3)   :: &
                                              scr1

real, intent(in), dimension(m2, m3)        :: &
                                              dxt, &
                                              sflux_u, &  !not used in this subroutine
                                              sflux_v, &  !not used in this subroutine
                                              sflux_w, &  !not used in this subroutine
                                              sflux_t, &  !not used in this subroutine
                                              rtgt
                                              

real, intent(in), dimension(m2,m3)      :: flpw

integer :: i,j,k,lpw

real :: sqrttkep,tket2,c1,sclu,scl

do j = ja,jz
   do i = ia,iz
      lpw = nint(flpw(i,j))
      tket2 = tket(lpw,i,j)
      c1 = 3.375 * rtgt(i,j) / dxt(i,j) ** 2  ! 3.375 = 1.5 ** 3
      do k = lpw,m1-1
         sqrttkep = sqrt(tkep(k,i,j))
         sclu = (c1 * (zm(k) - zm(k-1))) ** .333333
         if (vt3dj(k,i,j) .lt. 1.e-10) then
            scl = sclu                                   ! unstable case
         else
            scl = .76 * sqrt(tkep(k,i,j) / vt3dj(k,i,j))   !stable case
            if (scl .gt. sclu) scl = sclu
         endif

! now, vt3dh is strain rate and scr2 is dw/dz

         scr1(k,i,j) = .1 * scl * sqrttkep
         vctr25(k) = scr1(k,i,j) * vt3dh(k,i,j)  &
            - .333333 * tkep(k,i,j)  &
            * (vt3da(k,i,j) + vt3dc(k,i,j) + scr2(k,i,j))
         vt3dh(k,i,j) = scr1(k,i,j) * (1. + 2. * scl / sclu)
         scr2(k,i,j) = scr1(k,i,j)

! now, vt3dh is vkh and scr2 is hkm

         vctr26(k) = .19 + .51 * scl / sclu
         vctr27(k) = - vt3dh(k,i,j) * vt3dj(k,i,j)
         vctr28(k) = vctr26(k) * tkep(k,i,j) * sqrttkep / scl
         tket(k,i,j) = tket(k,i,j) + vctr25(k) + vctr27(k)  &
            - vctr28(k)
      enddo
      tket( lpw ,i,j) = tket2 + vctr25(3) + vctr27(2) - vctr28(2)

   enddo
enddo

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     call friclyr(nzp,nxp,nyp,a(iscr1),a(iustarl),a(itstarl)
!    +    ,a(iustarw),a(itstarw),a(ipctlnd),a(itheta),a(irtgt))
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

return
end subroutine mxtked



!     *****************************************************************

subroutine tkemy(m1,m2,m3,ia,iz,ja,jz,ibcon,jd,i0,j0  &
   ,tkep,tket,vt3dh,vt3di,vt3dj,scr1,rtgt  &
   ,theta,dn0,up,vp,wp,sflux_u,sflux_v,sflux_w,sflux_t,tkep2,flpw,flpu,flpv,sigw)

use mem_grid, only:  &
              zt, &       ! intent(in)
              zm, &       ! intent(in)
              nstbot      ! intent(in)

use mem_scratch, only: &
                 vctr30, &     ! intent(out)
                 vctr31, &     ! intent(out)
                 vctr1, &      ! intent(out)
                 vctr32, &     ! intent(out)
                 vctr33, &     ! intent(out)
                 vctr9, &      ! intent(out)
                 vctr5         ! intent(out)

use rconstants, only: &
                onethird, &    ! intent(in)
                vonk, &        ! intent(in)
                tkmin, &       ! intent(in)
                grav,  &       ! intent(in)
                sigwmin        ! intent(in)

implicit none

integer, intent(in) :: m1,m2,m3,ia,iz,ja,jz,ibcon,jd,i0,j0

real, intent(in)   , dimension(m1,m2,m3) :: tkep,  &
                                            vt3dj, &
                                            theta, &
                                            dn0,   &
                                            up,    &
                                            vp,    &
                                            wp


real, intent(out)  , dimension(m1,m2,m3) :: vt3dh, &
                                            scr1,  &
                                            sigw
                                            
                                         
real, intent(inout), dimension(m1,m2,m3) :: tket,  &
                                            vt3di


real, intent(in), dimension(m2,m3) :: rtgt,     & 
                                      sflux_u,  &
                                      sflux_v,  &
                                      sflux_w,  &
                                      sflux_t

real, intent(out), dimension(m1)      :: tkep2

real, intent(in), dimension(m2,m3) :: flpw, &
                                      flpu, &
                                      flpv

integer :: i,j,k,k2
integer :: lpu,lpu1,lpv,lpv1

real :: a1,a2,b1,b2,c1,aux1,aux2,rf1,rf2,rf3,rf4,wght1,wght3,sumtkz,sumtk  &
   ,al0,tket2,ri,rf,shr,smr,tker,qq,ssmf,shf,sh0,ssm,aux,gm,gh,sm1,sm2  &
   ,sh1,sh2,dzloc,ssum


data a1,a2,b1,b2,c1/0.92,0.74,16.6,10.1,0.08/
data aux1,aux2/0.758964199,2.58286747/
data rf1,rf2,rf3,rf4/1.,0.191232309,0.223117196,0.234067819/

!        1 - mellor and yamada (after andre et al, 1978)

wght3=1.
!c      wght3=zt(nint(flpw(2))/zt(3)
wght1=1.0-wght3

do j=ja,jz
   do i=ia,iz
     k2=nint(flpw(i,j))
!     wght3=zt(k2)/zt(k2+1)
     do k=k2,m1-1
      !if(i+i0==28 .and. j+j0==28) then
      !   tkep(k,i,j)=1.
      !   vt3di(k,i,j)=.286e-3
      !   vt3dj(k,i,j)=-.1e-5
      !endif
         tkep2(k) = 2.0 * tkep(k,i,j)
         vctr30(k) = sqrt(tkep2(k))
         vctr31(k)=.75*vctr30(k)/sqrt(max(1.e-20,vt3dj(k,i,j)))
         vctr1(k)=(zt(k)-zm(k2-1))*rtgt(i,j)
         dzloc=(zm(k)-zm(k-1))*rtgt(i,j)
         vctr33(k)=vctr30(k)*dzloc
         vctr32(k)=vctr33(k)*vctr1(k)
      enddo

      sumtkz=ssum(m1-k2,vctr32(k2),1)
      sumtk =ssum(m1-k2,vctr33(k2),1)
      al0=.1*sumtkz/sumtk
      tket2=tket(k2,i,j)

      do k=k2,m1-1
         vctr9(k)=min(vonk*vctr1(k)/(1.+vonk*vctr1(k)/al0)  &
              ,vctr31(k))

! --- for growing turbulence use helfand and labraga's modified sm and sh (sh0)

         ri=min(vt3dj(k,i,j)/max(vt3di(k,i,j),1.e-11),0.190)
         rf=min(0.6588*(ri+0.1776-sqrt(ri*(ri-0.3221)  &
              +0.03156)),0.16)
         shr=aux2*(rf-rf2)/(rf-rf1)
         smr=aux1*(rf-rf4)/(rf-rf3)*shr

         aux=vctr9(k)*vctr9(k)/tkep2(k)
         gm=aux*vt3di(k,i,j)
         gh=-aux*vt3dj(k,i,j)

         tker=max(16.6*vctr9(k)*vctr9(k)*(smr*vt3di(k,i,j)  &
            -shr*vt3dj(k,i,j)),2.*tkmin)

         if(tker.gt.tkep2(k) )then
            qq=sqrt(tkep2(k)/tker)
            ssmf=qq*smr
            shf=qq*shr
            sh0=shf
            ssm=ssmf
         else

! --- for decaying turbulence use mellor and yamada's sm and sh (sh0)
            sm1=0.6992-9.33948672*gh
            sm2=1.-(36.7188-187.4408515*gh+88.83949824*gm)*gh  &
                 +5.0784*gm
            ssm=sm1/max(sm2,1e-10)
            sh1=0.74-4.0848*ssm*gm
            sh2=1.-30.5916*gh
            sh0=sh1/max(sh2,1e-10)
         endif
         
         !-------- MLO. Finding sigma-w for cumulus parameterization ----------------------!
         sigw(k,i,j) = sqrt(max(sigwmin*sigwmin,                 &
                       tkep2(k)*(onethird - 2.*a1*ssm*gm + 4.*a1*sh0*gh)))
         !---------------------------------------------------------------------------------!
         
         
         scr1(k,i,j)=vctr9(k)*vctr30(k)*ssm
         vt3dh(k,i,j)=vctr9(k)*vctr30(k)*sh0
         vctr5(k)=scr1(k,i,j)*vt3di(k,i,j)  &
            -vt3dh(k,i,j)*vt3dj(k,i,j)
         tket(k,i,j)=tket(k,i,j) + 0.5 * vctr5(k)  &
            -tkep2(k) * sqrt(max(1.e-20,tkep2(k)))/(vctr9(k)*16.6)
     ! if(i+i0==28 .and. j+j0==28 .and. k<20) then
     !  print '(i3,12e13.4)',k,scr1(k,i,j),vctr9(k),vctr30(k),ssm

     !  print '(i3,12e13.4)',k,tkep(k,i,j),tket(k,i,j)  &
     !                      ,.5*scr1(k,i,j)*vt3di(k,i,j)  &
     !                      ,scr1(k,i,j),vt3di(k,i,j)  &
     !                      ,-.5*vt3dh(k,i,j)*vt3dj(k,i,j) &
     !                      ,vt3dh(k,i,j),vt3dj(k,i,j) &
     !       ,-tkep2(k) * sqrt(max(1.e-20,tkep2(k)))/(vctr9(k)*16.6)
      ! print '(i3,12f12.7)',k,sflux_u(i,j)/dn0(k2,i,j),sflux_v(i,j)/dn0(k2,i,j)  &
      !                       ,sflux_w(i,j)/dn0(k2,i,j),sflux_t(i,j)/dn0(k2,i,j)
      !   enddo
      !endif
         scr1(k,i,j)=scr1(k,i,j)*dn0(k,i,j)
         vt3dh(k,i,j)=vt3dh(k,i,j)*dn0(k,i,j)
         vt3di(k,i,j)=0.2*vctr9(k)*vctr30(k)*dn0(k,i,j)
      enddo
      if(nstbot.eq.1)then
!c
!c  experimental diagnostic equilibrium form - may try this someday
!c
!   CJT - 2 July 2001 - this seems to be way too small for larger grid spacings

    !     tket(k2,i,j) = (3.25 * sqrt(sflux_u(i,j)**2+sflux_v(i,j)**2) &
    !                     -tkep(k2,i,j))/dtlt
!c


         lpu  = nint(flpu(i,j))
         lpu1 = nint(flpu(i-1,j))
         lpv  = nint(flpv(i,j))
         lpv1 = nint(flpv(i,j-jd))


         tket(k2,i,j) =tket2  &
            + 0.5 * wght3 * vctr5(k2+1)   &
            + wght1 *  &
             (-(sflux_u(i,j)*(up( lpu ,i,j)+up( lpu1 ,i-1,j))  &
               +sflux_v(i,j)*(vp( lpv ,i,j)+vp( lpv1 ,i,j-jd))  &
               +sflux_w(i,j)*wp(k2,i,j))/vctr1(k2)  &
               +sflux_t(i,j)*grav/theta(k2,i,j))/dn0(k2,i,j) &
            
            - tkep2(k2)*sqrt(max(1.e-20,tkep2(k2)))/(vctr9(k2)*16.6)

      endif
      
!      if(i+i0==28 .and. j+j0==28) then
!         do k=2,10
!            print '(i3,12f10.5)',k,tkep(k,i,j),tket(k,i,j),scr1(k,i,j) &
!                                  ,vt3dh(k,i,j),vt3di(k,i,j)
!         enddo
!       print*,'-------------------------------',time
!      endif
      
   enddo
enddo

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     call friclyr(nzp,nxp,nyp,a(iscr1),a(iustarl),a(itstarl)
!    +    ,a(iustarw),a(itstarw),a(ipctlnd),a(itheta),a(irtgt))
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

return
end subroutine tkemy

!     *****************************************************************

!     subroutine friclyr(m1,m2,m3,scr1,ustarl,tstarl,ustarw,tstarw
!    +    ,pctlnd,theta,rtgt)
!     include 'rcommons.h'
!     dimension scr1(m1,m2,m3),ustarl(m2,m3),ustarw(m2,m3)
!    +   ,tstarl(m2,m3),tstarw(m2,m3),pctlnd(m2,m3),theta(m1,m2,m3)
!    +   ,rtgt(m2,m3)
!
!     do j=1,n3
!        do i=1,n2
!           ust=ustarl(i,j)*pctlnd(i,j)+ustarw(i,j)*(1.-pctlnd(i,j))
!           tst=tstarl(i,j)*pctlnd(i,j)+tstarw(i,j)*(1.-pctlnd(i,j))
!           if(abs(tst).lt.1.e-10)tst=1.e-10
!           obuklen=theta(2,i,j)*ust*ust/(3.92*tst)
!           do k=2,m1-1
!              zst=zt(k)*rtgt(i,j)
!              xi=zst/obuklen
!              if(xi.le.0.)then
!                 phim=(1.-15.*xi)**-0.25
!              else
!                 phim=1.+4.7*xi
!              endif
!              scr1(k,i,j)=min(scr1(k,i,j),0.4*ust*zst/phim)
!           enddo
!        enddo
!     enddo
!     return
!     end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!     *****************************************************************

subroutine tkeinit(n1,n2,n3,ia,iz,ja,jz)

use mem_grid, only: &
              ngrid, &         !INTENT(IN)
              zt               !INTENT(IN)
              
use mem_turb, only: &
              turb_g           !INTENT(INOUT)

use ke_coms, only: &
             c_eps             !INTENT(IN)

use rconstants, only: &
                tkmin, &       !INTENT(IN)
                vonk           !INTENT(IN)
                
implicit none
integer, intent(in) :: n1,n2,n3,ia,iz,ja,jz

integer :: i,j,k
real :: epsmin

!        limit the values to the minimum

if( associated(turb_g(ngrid)%tkep) ) then
   do j = ja,jz
      do i = ia,iz
         do k = 1,n1
            turb_g(ngrid)%tkep(k,i,j) = max(tkmin,turb_g(ngrid)%tkep(k,i,j))
         enddo
      enddo
   enddo
endif
if( associated(turb_g(ngrid)%epsp) ) then
   epsmin=(c_eps*(tkmin**1.5))/(vonk*zt(2))
   do j = ja,jz
      do i = ia,iz
         do k = 1,n1
            turb_g(ngrid)%epsp(k,i,j) = max(epsmin,turb_g(ngrid)%epsp(k,i,j))
         enddo
      enddo
   enddo
endif

return
end subroutine tkeinit
