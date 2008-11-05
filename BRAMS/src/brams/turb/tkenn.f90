!==========================================================================================!
!==========================================================================================!
! Subroutine nakanishi                                                                     !
! Developed by Marcos Longo (Lab. MASTER - Univ. São Paulo/EPS Harvard University)         !
! Sao Paulo, April 29, 2005                                                                !
!                                                                                          !
! This routine is intended to find the TKE tendency by using a 2.5-level model             !
! based on Mellor-Yamada scheme. It's actually an update of the existent M-Y code          !
! in RAMS, outputing also the Sig-W (vertical velocity standard-deviation) and             !
! the Lagrangian time scale, Obukhov Length and PBL depth, using the following references: !
!                                                                                          !
! JANJIC, Z. I. Nonsingular implementation of the Mellor-Yamada level 2.5 scheme in the    !
!   NCEP meso model, office note # 437, National Centers for Environmental Prediction,     !
!   2001, 61 pp.                                                                           !
!                                                                                          !
! NAKANISHI, M. Improvement of the Mellor-Yamada turbulence closure model based on         !
!   large-eddy simulation data. Boundary-Layer Meteor., v. 99, p. 349-378, 2001.           !
!                                                                                          !
! NAKANISHI, M.; NIINO, H. An improved Mellor-Yamada level-3 model with condens-           !
!   ation physics: its design and verification. Boundary-Layer Meteor., v.112,             !
!   p. 1-31, 2004.                                                                         !
!                                                                                          !
! NAKANISHI, M.; NIINO, H. An improved Mellor-Yamada level-3 model with condensation       !
!   physics: its numerical stability and application to a regional prediction of advection !
!   fog. Boundary-Layer Meteor., vol. 119, p. 397-407, 2006.                               !
!                                                                                          !
! HANNA, S. R. Application in air pollution modeling. In: NIEUSWSTADT, F. M. T.;           !
!   VAN DOP, H. Atmospheric turbulence and air pollution modelling. Dordrecht: D.          !
!   Reidel Publishing Company, 1982, chap. 7, p. 275-310.                                  !
!                                                                                          !
! VOGEZELANG, D. H. P.; HOLTSLAG, A. M. Evaluation and model impacts of alternati-         !
!   ve boundary-layer height formulations. Boundary-Layer Meteor., v. 81, p. 245-          !
!   269, 1996.                                                                             !
!==========================================================================================!
!==========================================================================================!
subroutine nakanishi(m1,m2,m3,m4,ia,iz,ja,jz,jd,tkep,tket,vt3dd,vt3de,vt3dh,vt3di          &
                    ,vt3dj,scr1,rtgt,theta,rv,rtp,dn0,up,vp,vegz0,patchz0,tstar,ustar      &
                    ,patch_area,sflux_u,sflux_v,sflux_t,flpu,flpv,flpw,kpbl,pblhgt,lmo     &
                    ,tl,sigw)
   !---------------------------------------------------------------------------------------!
   !   Reference for variables, following Nakanishi (2001) and Nakanishi and Niino         !
   ! (2004,2006), which are based on Helfand and Labraga (1988) and Mellor and Yamada      !
   ! (1974,1982).                                                                          !
   !                                                                                       !
   ! dzloc  -> dz for integrals                                                            !
   ! scr1   -> Km×rho0                                                                     !
   ! sumtkz -> int_0^H z sqrt(e) dz                                                        !
   ! sumtk  -> int_0^H sqrt(e) dz                                                          !
   ! tkep   -> 1/2 q² = e (TKE)                                                            !
   ! vctr1  -> z                                                                           !
   ! vctr5  -> zeta
   ! vctr9  -> L                                                                           !
   ! vctr19 -> Gm                                                                          !
   ! vctr20 -> Gh                                                                          !
   ! vctr21 -> THETAV0                                                                     !
   ! vctr22 -> THETA V', used for PBL in the unstable case                                 !
   ! vctr23 -> Sm                                                                          !
   ! vctr24 -> Sh                                                                          !
   ! vctr29 -> zonal wind in thermodynamic grid                                            !
   ! vctr30 -> sqrt(2e) = q                                                                !
   ! vctr31 -> meridional wind in thermodynamic grid                                       !
   ! vctr32 -> q dz                                                                        !
   ! vctr33 -> q z dz                                                                      !
   ! vctr35 -> q³                                                                          !
   ! vctr36 -> Ps+Pb=(q³/L)×(Sm×Gm+Sh×Gh)                                                  !
   ! vctr37 -> epsilon=(q³/B1×L)                                                           !
   ! vctr38 -> 2e = q²´                                                                    !
   ! vt3dd  -> dU/dz at the u points                                                       !
   ! vt3de  -> dV/dz at the v points                                                       !
   ! vt3dh  -> Kh*rho0                                                                     !
   ! vt3di  -> (dU/dz)²+(dV/dz)² (M²) at the beginning, then Kq×rho0                       !
   ! vt3dj  -> g/theta_v d(THETA_v)/dz (N²)                                                   !
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   ! Loading modules                                                                       !
   !---------------------------------------------------------------------------------------!
   use mem_grid,     only: nstbot,zm,zt
   use mem_scratch,  only:  vctr1, vctr5, vctr9,vctr19,vctr20,vctr21,vctr22,vctr23,vctr24  &
                           ,vctr29,vctr30,vctr31,vctr32,vctr33,vctr38,vctr35,vctr36,vctr37

   use rconstants,     only: abslmomin,abswltlmin,cp,g,ltscalemax,sigwmin,vonk,lturbmin    &
                            ,tkmin,onethird

   use turb_constants, only:  nna1,nna2,nnb1,nnb2,nnc1,nnc2,nnc3,nnc4,nnc5                 &
                             ,nnq1,nnq2,nnq3,nnq4,nnq5,nnq6                                &
                             ,nngama1,nngama2,nnf1,nnf2,nnrf1,nnrf2,nnri1,nnri2,nnri3      &
                             ,nnrfc,nnce1a,nnce1b,nnce2,nnce3,nnce4,nncr1                  &
                             ,nno1,nno2,nno3,nno4, nno5, nno6, nno7, nno8, nnreq,nnrsl     &
                             ,nnmacheps

   use therm_lib     , only: virtt
   implicit none

   !---------------------------------------------------------------------------------------!
   !   Variable declaration section.                                                       !
   !---------------------------------------------------------------------------------------!
   !----- Arguments (Input/Output/Both) ---------------------------------------------------!
   integer  , intent(in)                          :: m1,m2,m3,m4,ia,iz,ja,jz,jd
   real     , intent(in)    , dimension(m2,m3)    :: flpu,flpv,flpw
   real     , intent(in)    , dimension(m2,m3)    :: rtgt, sflux_u, sflux_v, sflux_t
   real     , intent(in)    , dimension(m1,m2,m3) :: tkep, vt3dd, vt3de, theta
   real     , intent(in)    , dimension(m1,m2,m3) :: vt3dj, rv, rtp, dn0, up, vp
   real     , intent(in)    , dimension(m2,m3,m4) :: vegz0, patchz0, ustar
   real     , intent(in)    , dimension(m2,m3,m4) :: tstar, patch_area 
   integer  , intent(out)   , dimension(m2,m3)    :: kpbl
   real     , intent(out)   , dimension(m2,m3)    :: pblhgt, lmo
   real     , intent(out)   , dimension(m1,m2,m3) :: vt3dh, scr1, tl, sigw
   real     , intent(inout) , dimension(m1,m2,m3) :: vt3di,tket
   !----- Internal variables --------------------------------------------------------------!
   integer                                        :: i,j,k,p,k2u,k2um1,k2v,k2vm1,k2w 
   real                                           :: weightsurf,du0dz,dv0dz,sumtkz,sumtk       
   real                                           :: tket2,aux,dzloc,rf,ri,sh2,sm2,gm2
   real                                           :: qlevel2,oneminusalpha,oneminusalpha2
   real                                           :: e1,e2,e3,e4,r1,r2,wltl0,lt,ls,lb
   real                                           :: suminv,z0w,ustarw,tstarw,wstarw,qc
   real                                           :: dzloc1up,dzloc1dn,dzloc2,dq3dz
   real                                           :: dlsmdz,d2q3dz2,janjc,janjd,janjg
   real                                           :: janjh,janji,janjp1,janjt1
   logical                                        :: stable,neutral
   !----- Function ------------------------------------------------------------------------!
   real, external                                 :: cbrt
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !   Big loops on horizontal.                                                            !
   !---------------------------------------------------------------------------------------!
   do j=ja,jz
      do i=ia,iz
         k2u   = nint(flpu(i,j))
         k2um1 = nint(flpu(i-1,j))
         k2v   = nint(flpv(i,j))
         k2vm1 = nint(flpv(i,j-jd))
         k2w   = nint(flpw(i,j))
         !---------------------------------------------------------------------------------!
         !  k2w is the first useful level                                                  !
         !---------------------------------------------------------------------------------!
         do k=k2w,m1-1
            vctr30(k) = max(sqrt(2.0 * tkep(k,i,j)),sqrt(2*tkmin))    ! q 
            vctr38(k) = max(2.0 * tkep(k,i,j),2*tkmin)                ! q²
            vctr35(k) = vctr38(k)*vctr30(k)                           ! q³

            !------------------------------------------------------------------------------!
            !  Restrictions in M², just to avoid singularities (following Janjic 2001)     !
            !------------------------------------------------------------------------------!
            if (vt3di(k,i,j) < nnmacheps)  vt3di(k,i,j)=(1.+nnmacheps)*nnmacheps*nnreq

            !------------------------------------------------------------------------------!
            !    Determining some integrated variables which will be necessary for length  !
            ! scale derivation.                                                            !
            !------------------------------------------------------------------------------!
            vctr1(k)=(zt(k)-zm(k2w-1))*rtgt(i,j) 
            dzloc=(zm(k)-zm(k-1))*rtgt(i,j)
            vctr33(k)=vctr30(k)*dzloc
            vctr32(k)=vctr33(k)*vctr1(k)
            !----- Deriving some variables that are needed for PBL depth estimation -------!
            vctr21(k)=virtt(theta(k,i,j),rv(k,i,j),rtp(k,i,j))
            vctr29(k)=up(k,i,j)
            vctr31(k)=vp(k,i,j)
            vctr29(k)=0.5*(up(k,i,j)+up(k,i-1,j)) 
            vctr31(k)=0.5*(vp(k,i,j)+vp(k,i,j-jd))
         end do !----- k=k2,m1-1 !
    
         !----- The sum is between k2w and m1-1: k1 and m1 are just boundaries ------------!
         sumtkz=sum(vctr32(k2w:(m1-1)))
         sumtk =sum(vctr33(k2w:(m1-1)))
         !---------------------------------------------------------------------------------!

         !----- Finding the average dU/dz and dV/dz at the thermodynamic level k2 ---------!
         du0dz=0.25 * ( vt3dd(k2um1,i-1,j)    + vt3dd(k2u,i,j)                             &
                      + vt3dd(k2um1-1,i-1,j)  + vt3dd(k2u-1,i,j) )
         dv0dz=0.25 * ( vt3de(k2vm1,i,j-jd)   + vt3de(k2v,i,j)                             &
                      + vt3de(k2vm1-1,i,j-jd) + vt3de(k2v-1,i,j) )
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !    Here I'll add a weighted average to z0, ustar and tstar, to avoid patch      !
         ! dependance.                                                                     !
         !---------------------------------------------------------------------------------!
         z0w=0.                                                      
         ustarw=0.
         tstarw=0.
         do p=1,m4
            z0w    =    z0w + max(vegz0(i,j,p),patchz0(i,j,p)) * patch_area(i,j,p)
            ustarw = ustarw + ustar(i,j,p)                     * patch_area(i,j,p)
            tstarw = tstarw + tstar(i,j,p)                     * patch_area(i,j,p)
         end do
         !---------------------------------------------------------------------------------!
         ! Here I'll find the <w'Theta'>g and truncate for small values                    !
         !---------------------------------------------------------------------------------!
         if (abs(sflux_t(i,j)/dn0(k2w,i,j)) < abswltlmin) then
            wltl0=sign(abswltlmin,sflux_t(i,j)/dn0(k2w,i,j))
         else
            wltl0=sflux_t(i,j)/dn0(k2w,i,j)
         end if

         !---------------------------------------------------------------------------------!
         ! Finding the Obukhov length (LMO)                                                !
         ! From Nakanishi (2001), Eq. 4:                                                   !
         !              THETA0 u*³                                                         !
         ! LMO = - --------------------                                                    !
         !           k g <w theta>(g)                                                      !
         !                                                                                 !
         !   Since LMO can be close to zero, but either positive or negative, it must be   !
         ! truncated to a small value, but keeping the same sign                           !
         !---------------------------------------------------------------------------------!
         lmo(i,j)= - vctr21(k2w) * ustarw * ustarw * ustarw / (vonk * g * wltl0)
         if (abs(lmo(i,j)) < abslmomin) lmo(i,j)=sign(abslmomin,lmo(i,j))

         do k=k2w,m1-1
            !------------------------------------------------------------------------------!
            !    Calculating the length scale L, based on Nakanishi (2001), and GH (nega-  !
            ! tive of dimensionless square of Brunt Väisälä frequency) and GM (dimension-  !
            ! less square of mean shear)                                                   !
            !------------------------------------------------------------------------------!
            vctr5(k) = vctr1(k)/lmo(i,j) ! zeta
            !------------------------------------------------------------------------------!
            ! Finding Ls, following equation A2 of Nakanishi and Niino (2004):             !
            !------------------------------------------------------------------------------!
            if ( vctr5(k) >= 1 ) then
               ls = vonk * vctr1(k) / nnq1                       
            elseif (vctr5(k) >= 0 ) then                       
               ls = vonk * vctr1(k) / (1 + nnq2 * vctr5(k))       
            else !if (vctr5(k) < 0) then                       
               ls = vonk * vctr1(k) * (1 - nnq3 * vctr5(k))**nnq4   
            end if                                            
            !------------------------------------------------------------------------------!
            ! Finding Lt, following equation A3 of Nakanishi and Niino (2004):             !
            !------------------------------------------------------------------------------!
            lt = nnq5 * sumtkz / sumtk
            !------------------------------------------------------------------------------!
            ! Finding Lb, following equation A4 of Nakanishi and Niino (2004):             !
            !------------------------------------------------------------------------------!
            if (vt3dj(k,i,j) > 0 .and. vctr5(k) >= 0) then
               lb = vctr30(k) / max(sqrt(vt3dj(k,i,j)),1.e-10)
            elseif (vt3dj(k,i,j) > 0 .and. vctr5(k) < 0) then
               qc= cbrt((g/vctr21(k2w))*wltl0*lt)
               lb = (1. + nnq6 *sqrt(qc/max((lt*sqrt(vt3dj(k,i,j))),1.e-20)))*             &
                          vctr30(k)/max(sqrt(vt3dj(k,i,j)),1.e-10)
            else
               lb = 1.e20 !----- Making it very large but without risking overflow. -------!
            end if
            suminv=(1./ls)+(1./lt)+(1./lb) !-----Equation A1, Nakanishi and Niino (2004) --!
            vctr9(k)=1./suminv
            !------------------------------------------------------------------------------!
            !    Restriction on L: from Nakanishi and Niino (2006), based on Janjic(2001). !
            ! This is to guarantee that Cw, defined as <w²>/q² < Rsl (about 0.14)          !
            !------------------------------------------------------------------------------!
            if (vt3dj(k,i,j) > 0) then
               janjc=(nno1*vt3dj(k,i,j)+nno2*vt3di(k,i,j))*vt3dj(k,i,j)
               janjd=(nno3*vt3di(k,i,j)+nno4*vt3dj(k,i,j))
               janjg=(nno5*vt3dj(k,i,j)+nno6*vt3di(k,i,j))*vt3dj(k,i,j)-3.*nnrsl*janjc
               janjh=nno7*vt3di(k,i,j)+nno8*vt3dj(k,i,j)-3.*nnrsl*janjd
               janji=1-3.*nnrsl
               janjt1=(-janjh+sqrt(janjh*janjh-4.*janjg*janji))/(2.*janji)
               if (janjt1 > nnmacheps*nnmacheps) then
                  vctr9(k)=min(vctr9(k),vctr30(k)/sqrt(janjt1))
               end if
            end if
         end do !----- k=k2,m1-1 !

         !---------------------------------------------------------------------------------!
         !   Obtaining an estimation of the PBL height. First, I need to check whether the !
         ! PBL is stable, neutral, or convective. This is obtained by verifying whether    !
         ! z/LMO is > 0, =0 or <0                                                          !
         !---------------------------------------------------------------------------------!
         stable = vctr5(k2w)  >=  1.e-4
         neutral = abs(vctr5(k2w)) < 1.e-4 
         if (stable .or. neutral) then 
            !------------------------------------------------------------------------------!
            !   If the boundary layer is STABLE or NEUTRAL, then we follow Vogelezang and  !
            ! Holtslag (1996) equation (3); (SBL is defined when (w'theta')g < 0). The PBL !
            ! depth is defined as the first level where average Ri > 0.25.                 !
            !------------------------------------------------------------------------------!
            kpbl(i,j) = k2w
            pblhgt(i,j)=vctr1(k2w)
            sboundlay: do k=k2w+1,m1-1
               ri = g * (vctr21(k)-vctr21(k2w)) * (vctr1(k)-vctr1(k2w))                    &
                  / ( vctr21(k2w) * ( (vctr29(k)-vctr29(k2w)) * (vctr29(k)-vctr29(k2w))    &
                                    + (vctr31(k)-vctr31(k2w)) * (vctr31(k)-vctr31(k2w))    &
                                    + 100.*ustarw*ustarw ) )
               if (ri >= 0.25) then
                  kpbl(i,j) = k
                  pblhgt(i,j) = 0.5 * (vctr1(k-1)+vctr1(k))
                  exit sboundlay
               end if
            end do sboundlay
         else
            !------------------------------------------------------------------------------!
            !    -> Or, if the PBL is convective, then the PBL is defined as the level     !
            !       where the first minimum of w'theta'                                    !
            !------------------------------------------------------------------------------!
            kpbl(i,j) = k2w
            pblhgt(i,j)=vctr1(k2w)
            aux= wltl0*g/(cp * vctr21(k2w))
            convmixlay: do k=k2w+1,m1-1
               kpbl(i,j) = k
               pblhgt(i,j)=0.5*(vctr1(k)+vctr1(k-1))
               wstarw=cbrt(aux * pblhgt(i,j) )
               vctr22(k)=vctr21(k2w)+8.5*wltl0/(wstarw*cp)
               ri = g * (vctr21(k)-vctr22(k)) * (vctr1(k)-vctr1(k2w))                      &
                  / ( vctr22(k) * ((vctr29(k)-vctr29(k2w)) * (vctr29(k)-vctr29(k2w))       &
                    + (vctr31(k)-vctr31(k2w)) * (vctr31(k)-vctr31(k2w))                    &
                    + 100.*ustarw*ustarw) )
               if (ri >= 0.25) exit convmixlay
            end do convmixlay
         end if
         !---------------------------------------------------------------------------------!

         do k=k2w,m1-1

            !------------------------------------------------------------------------------!
            !    Calculating both gradient and flux Richardson numbers (ri and rf, respec- !
            ! tively):                                                                     !
            !------------------------------------------------------------------------------!
            ri=vt3dj(k,i,j)/vt3di(k,i,j)
            rf=nnri1*(ri+nnri2-sqrt(ri*ri-nnri3*ri+nnri2*nnri2))
            !------------------------------------------------------------------------------!
            ! Finding the SH and SM of the Level 2 model (Nakanish, 2001):                 !
            !------------------------------------------------------------------------------!
            if (rf < nnrfc) then
               sh2=3.*nna2*(nngama1+nngama2)*(nnrfc-rf)/(1.-rf)
               sm2=nna1*nnf1*(nnrf1-rf)*sh2/(nna2*nnf2*(nnrf2-rf))
               qlevel2=max(vctr9(k)*sqrt(nnb1*vt3di(k,i,j)*sm2*(1.-rf)),sqrt(2.*tkmin))
            else
               qlevel2=sqrt(2.*tkmin)
            end if
            !------------------------------------------------------------------------------!
            ! Finding Gm and Gh                                                            !
            !------------------------------------------------------------------------------!
            aux=vctr9(k)*vctr9(k)/vctr38(k)
            vctr19(k)=aux*vt3di(k,i,j)       !----- Gm ------------------------------------!
            vctr20(k)=-aux*vt3dj(k,i,j)      !----- Gh ------------------------------------!

            !------------------------------------------------------------------------------!
            !     Case of growing turbulence, it must use the modified Level 2½ model      !
            ! (Nakanishi and Niino, 2004), otherwise it just uses the original Level 2½.   !
            !------------------------------------------------------------------------------!
            if (vctr30(k) < qlevel2) then
               oneminusalpha=vctr30(k)/max(qlevel2,1.e-10)
               oneminusalpha2=oneminusalpha*oneminusalpha
            else
               oneminusalpha=1.
               oneminusalpha2=1.
            end if
            e1=1.+oneminusalpha2*(nnce1a*vctr19(k)-nnce1b*vctr20(k))
            e2=-oneminusalpha2*nnce2*vctr20(k)
            e3=oneminusalpha2*nnce3*vctr19(k)
            e4=1.-oneminusalpha2*nnce4*vctr20(k)
            r1=oneminusalpha*nncr1
            r2=oneminusalpha*nna2
            vctr23(k)=(r2*e2-r1*e4)/(e2*e3-e1*e4) !----- Sm -------------------------------!
            vctr24(k)=(r1*e3-r2*e1)/(e2*e3-e1*e4) !----- Sh -------------------------------!

            !------------------------------------------------------------------------------!
            !    Deriving the convective velocity standard-deviation.                      !
            !    The max was inserted because there is no physical limitation which effec- !
            ! tively prevents the variance to be negative (although in preliminary tests   !
            ! the negative values were significantly smaller than the positive ones.       !
            !------------------------------------------------------------------------------!
            sigw(k,i,j) = sqrt(max(vctr38(k)*(onethird-2*nna1*vctr23(k)*vctr19(k)          &
                        + 4.*nna1*(1. - nnc2)*vctr24(k)*vctr20(k)),sigwmin*sigwmin))

         end do
         !---------------------------------------------------------------------------------!

         !----- Finding the surface weight ------------------------------------------------!
         weightsurf=vctr1(k2w)/vctr1(k2w+1)

         !----- Saving the TKE tendency at the bottom for later ---------------------------!
         tket2=tket(k2w,i,j)
         
         do k=k2w,m1-1
            !------------------------------------------------------------------------------!
            ! The vertical Lagrangian timescale is found following Hanna (1982)            !
            !------------------------------------------------------------------------------!
            if (stable) then
            !------------------------------------------------------------------------------!
            ! Stable PBL, so we shall use equation 7.24                                    !
            !------------------------------------------------------------------------------!
               tl(k,i,j)= 0.10 * (pblhgt(i,j)**0.20) * (vctr1(k)**0.80) / sigw(k,i,j)
            elseif (neutral) then
            !------------------------------------------------------------------------------!
            !    This is the formula proposed by Hanna (1982), but without the Coriolis    !
            ! term: later papers, such as Wilson (200, JAM), also disconsider this.        !
            !------------------------------------------------------------------------------!
               tl(k,i,j)=0.5*vctr1(k)/sigw(k,i,j)
            !------------------------------------------------------------------------------!
            !   For unstable PBLs we should follow equation 7.17, considering all possible !
            ! conditions                                                                   !
            !------------------------------------------------------------------------------!
            elseif (vctr1(k) >= 0.1* pblhgt(i,j)) then ! if z/h > 0.1
               ! If PBL is too low, this number can cause underflow in exp. ---------------!
               aux       = max(40.,5.0*vctr1(k)/pblhgt(i,j))
               tl(k,i,j) = 0.15*pblhgt(i,j)*(1-exp(-aux))/sigw(k,i,j)
            elseif ((z0w - vctr1(k)) <= lmo(i,j)) then ! -(z-z0)/Lmo > 1 and z/h < 0.1
               tl(k,i,j)= 0.59*vctr1(k)/sigw(k,i,j)
            else                                       ! -(z-z0)/Lmo < 1 and z/h < 0.1
               tl(k,i,j)= 0.10*vctr1(k)/(sigw(k,i,j)*(0.55+0.38*(vctr1(k)-z0w)/lmo(i,j)))
            end if
            tl(k,i,j)=min(tl(k,i,j),ltscalemax) ! Avoiding exaggerated values...
            !------------------------------------------------------------------------------!

            !------------------------------------------------------------------------------!
            ! Finding the TKE tendency.                                                    !
            !------------------------------------------------------------------------------!
            !----- Finding Ps+Pb-epsilon=(q³/L)×(SmGm+ShGh-1/B1) --------------------------!
            vctr36(k)=vctr35(k)*(vctr23(k)*vctr19(k)+vctr24(k)*vctr20(k))/vctr9(k)
            vctr37(k)=vctr35(k)/(nnb1*vctr9(k))
            !------------------------------------------------------------------------------!
            !    Finding the tendency term, based on Nakanishi and Niino equation 5. Note  !
            ! that all the 2 factors vanish because here it's the TKE tendency, and        !
            ! TKE=q²/2.                                                                    !
            !------------------------------------------------------------------------------!
            tket(k,i,j)=tket(k,i,j)+vctr36(k)-vctr37(k)

            !----- Km=L×q×Sm (scr1) -------------------------------------------------------!
            scr1(k,i,j)=vctr9(k)*vctr30(k)*vctr23(k)*dn0(k,i,j)
            !----- Kh=L×q×Sh (vt3dh) ------------------------------------------------------!
            vt3dh(k,i,j)=vctr9(k)*vctr30(k)*vctr24(k)*dn0(k,i,j)
            !-----  Kq=L×q×Sq (vt3di). N&N (2004), proposed Sq=2*Sm instead of Sq=0.2 -----!
            vt3di(k,i,j)=2.*scr1(k,i,j) ! Since Km=L×q×Sm... 
         end do

         !---------------------------------------------------------------------------------!
         !    Different closure for the surface, using surface fluxes. If weightsurf =0    !
         ! it'll ignore.                                                                   !
         !---------------------------------------------------------------------------------!
         if (nstbot == 1 .and. weightsurf > 0) then
            tket(k2w,i,j) = tket2 -vctr37(k2w) + (1.-weightsurf)*vctr36(k2w+1)             &
                          + weightsurf*(-sflux_u(i,j)*du0dz-sflux_v(i,j)*dv0dz             &
                          + g*sflux_t(i,j)/theta(k2w,i,j))/dn0(k2w,i,j)
         end if
      end do
   end do
   return
end subroutine nakanishi
!==========================================================================================!
!==========================================================================================!
