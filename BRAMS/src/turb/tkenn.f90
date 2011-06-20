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
   use mem_grid,       only : nstbot               & ! intent(in)
                            , zm                   & ! intent(in)
                            , zt                   ! ! intent(in)
   use mem_scratch,    only : zagl     => vctr1    & ! intent(out)
                            , zeta     => vctr5    & ! intent(out)
                            , lscal    => vctr9    & ! intent(out)
                            , gm       => vctr19   & ! intent(out)
                            , gh       => vctr20   & ! intent(out)
                            , thetav0  => vctr21   & ! intent(out)
                            , thetavp  => vctr22   & ! intent(out)
                            , ssm      => vctr23   & ! intent(out)
                            , ssh      => vctr24   & ! intent(out)
                            , uspd     => vctr29   & ! intent(out)
                            , q        => vctr30   & ! intent(out)
                            , vspd     => vctr31   & ! intent(out)
                            , qdz      => vctr32   & ! intent(out)
                            , qzdz     => vctr33   & ! intent(out)
                            , qqq      => vctr35   & ! intent(out)
                            , pspluspb => vctr36   & ! intent(out)
                            , epsq     => vctr37   & ! intent(out)
                            , qq       => vctr38   ! ! intent(out)
   use rconstants,     only : abslmomin            & ! intent(in)
                            , abswltlmin           & ! intent(in)
                            , cp                   & ! intent(in)
                            , grav                 & ! intent(in)
                            , ltscalemax           & ! intent(in)
                            , sigwmin              & ! intent(in)
                            , vonk                 & ! intent(in)
                            , lturbmin             & ! intent(in)
                            , tkmin                & ! intent(in)
                            , onethird             & ! intent(in)
                            , lnexp_max            ! ! intent(in)
   use leaf_coms     , only : ustmin               & ! intent(in)
                            , min_patch_area       ! ! intent(in)
   use turb_coms     , only : a1      => nna1      & ! intent(in)
                            , a2      => nna2      & ! intent(in)
                            , b1      => nnb1      & ! intent(in)
                            , b2      => nnb2      & ! intent(in)
                            , c1      => nnc1      & ! intent(in)
                            , c2      => nnc2      & ! intent(in)
                            , c3      => nnc3      & ! intent(in)
                            , c4      => nnc4      & ! intent(in)
                            , c5      => nnc5      & ! intent(in)
                            , s1      => nns1      & ! intent(in)
                            , s2      => nns2      & ! intent(in)
                            , s3      => nns3      & ! intent(in)
                            , s4      => nns4      & ! intent(in)
                            , s5      => nns5      & ! intent(in)
                            , s6      => nns6      & ! intent(in)
                            , gama1   => nngama1   & ! intent(in)
                            , gama2   => nngama2   & ! intent(in)
                            , f1      => nnf1      & ! intent(in)
                            , f2      => nnf2      & ! intent(in)
                            , rf1     => nnrf1     & ! intent(in)
                            , rf2     => nnrf2     & ! intent(in)
                            , ri1     => nnri1     & ! intent(in)
                            , ri2     => nnri2     & ! intent(in)
                            , ri3     => nnri3     & ! intent(in)
                            , rfc     => nnrfc     & ! intent(in)
                            , ce1a    => nnce1a    & ! intent(in)
                            , ce1b    => nnce1b    & ! intent(in)
                            , ce2     => nnce2     & ! intent(in)
                            , ce3     => nnce3     & ! intent(in)
                            , ce4     => nnce4     & ! intent(in)
                            , cr1     => nncr1     & ! intent(in)
                            , o1      => nno1      & ! intent(in)
                            , o2      => nno2      & ! intent(in)
                            , o3      => nno3      & ! intent(in)
                            , o4      => nno4      & ! intent(in)
                            , o5      => nno5      & ! intent(in)
                            , o6      => nno6      & ! intent(in)
                            , o7      => nno7      & ! intent(in)
                            , o8      => nno8      & ! intent(in)
                            , req     => nnreq     & ! intent(in)
                            , rsl     => nnrsl     & ! intent(in)
                            , macheps => nnmacheps ! ! intent(in)
   use therm_lib,      only : virtt                ! ! function
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
   integer  , intent(inout) , dimension(m2,m3)    :: kpbl
   real     , intent(inout) , dimension(m2,m3)    :: pblhgt, lmo
   real     , intent(inout) , dimension(m1,m2,m3) :: vt3dh, scr1, tl, sigw
   real     , intent(inout) , dimension(m1,m2,m3) :: vt3di,tket
   !----- Internal variables --------------------------------------------------------------!
   integer                                        :: i,j,k,p,k2u,k2um1,k2v,k2vm1,k2w 
   real                                           :: weightsurf,du0dz,dv0dz,sumtkz,sumtk       
   real                                           :: tket2,aux,dzloc,rf,ri,sh2,sm2,gm2
   real                                           :: qlevel2,oneminusalpha,oneminusalpha2
   real                                           :: e1,e2,e3,e4,r1,r2,wltl0,lt,ls,lb
   real                                           :: suminv,z0w,ustarw,tstarw,wstarw,qc
   real                                           :: janjc,janjd,janjg
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
            q(k)   = max(sqrt(2.0 * tkep(k,i,j)),sqrt(2*tkmin))    ! q 
            qq(k)  = max(2.0 * tkep(k,i,j),2*tkmin)                ! q²
            qqq(k) = qq(k)*q(k)                                    ! q³

            !------------------------------------------------------------------------------!
            !  Restrictions in M², just to avoid singularities (following Janjic 2001)     !
            !------------------------------------------------------------------------------!
            if (vt3di(k,i,j) < macheps)  vt3di(k,i,j)=(1.+macheps)*macheps*req

            !------------------------------------------------------------------------------!
            !    Determining some integrated variables which will be necessary for length  !
            ! scale derivation.                                                            !
            !------------------------------------------------------------------------------!
            zagl(k) = (zt(k)-zm(k2w-1))*rtgt(i,j)
            dzloc   = (zm(k)-zm(k-1))*rtgt(i,j)
            qzdz(k) = q(k)*dzloc
            qdz(k)  = qzdz(k)*zagl(k)
            !----- Deriving some variables that are needed for PBL depth estimation -------!
            thetav0(k) = virtt(theta(k,i,j),rv(k,i,j),rtp(k,i,j))
            uspd(k)    = 0.5*(up(k,i,j)+up(k,i-1,j)) 
            vspd(k)    = 0.5*(vp(k,i,j)+vp(k,i,j-jd))
         end do !----- k=k2,m1-1 !

         !----- The sum is between k2w and m1-1: k1 and m1 are just boundaries ------------!
         sumtkz = sum( qdz(k2w:(m1-1)))
         sumtk  = sum(qzdz(k2w:(m1-1)))
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
         patchloop: do p=1,m4
            if (patch_area(i,j,p) < min_patch_area) cycle patchloop
            z0w    =    z0w + patchz0(i,j,p) * patch_area(i,j,p)
            ustarw = ustarw + ustar(i,j,p)   * patch_area(i,j,p)
            tstarw = tstarw + tstar(i,j,p)   * patch_area(i,j,p)
         end do patchloop
         ustarw = max(ustarw,ustmin)
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
         lmo(i,j)= - thetav0(k2w) * ustarw * ustarw * ustarw / (vonk * grav * wltl0)
         if (abs(lmo(i,j)) < abslmomin) lmo(i,j)=sign(abslmomin,lmo(i,j))

         do k=k2w,m1-1
            !------------------------------------------------------------------------------!
            !    Calculating the length scale L, based on Nakanishi (2001), and GH (nega-  !
            ! tive of dimensionless square of Brunt Väisälä frequency) and GM (dimension-  !
            ! less square of mean shear)                                                   !
            !------------------------------------------------------------------------------!
            zeta(k) = zagl(k)/lmo(i,j) ! zeta
            !------------------------------------------------------------------------------!
            ! Finding Ls, following equation A2 of Nakanishi and Niino (2004):             !
            !------------------------------------------------------------------------------!
            if ( zeta(k) >= 1 ) then
               ls = vonk * zagl(k) / s1
            elseif (zeta(k) >= 0 ) then
               ls = vonk * zagl(k) / (1 + s2 * zeta(k))
            else !if (zeta(k) < 0) then
               ls = vonk * zagl(k) * (1 - s3 * zeta(k))**s4
            end if
            !------------------------------------------------------------------------------!
            ! Finding Lt, following equation A3 of Nakanishi and Niino (2004):             !
            !------------------------------------------------------------------------------!
            lt = s5 * sumtkz / sumtk
            !------------------------------------------------------------------------------!
            ! Finding Lb, following equation A4 of Nakanishi and Niino (2004):             !
            !------------------------------------------------------------------------------!
            if (vt3dj(k,i,j) > 0 .and. zeta(k) >= 0) then
               lb = q(k) / max(sqrt(vt3dj(k,i,j)),1.e-10)
            elseif (vt3dj(k,i,j) > 0 .and. zeta(k) < 0) then
               qc= cbrt((grav/thetav0(k2w))*wltl0*lt)
               lb = (1. + s6 *sqrt(qc/max((lt*sqrt(vt3dj(k,i,j))),1.e-20)))*             &
                          q(k)/max(sqrt(vt3dj(k,i,j)),1.e-10)
            else
               lb = 1.e20 !----- Making it very large but without risking overflow. -------!
            end if
            suminv=(1./ls)+(1./lt)+(1./lb) !-----Equation A1, Nakanishi and Niino (2004) --!
            lscal(k)=1./suminv
            !------------------------------------------------------------------------------!
            !    Restriction on L: from Nakanishi and Niino (2006), based on Janjic(2001). !
            ! This is to guarantee that Cw, defined as <w²>/q² < Rsl (about 0.14)          !
            !------------------------------------------------------------------------------!
            if (vt3dj(k,i,j) > 0) then
               janjc=(o1*vt3dj(k,i,j)+o2*vt3di(k,i,j))*vt3dj(k,i,j)
               janjd=(o3*vt3di(k,i,j)+o4*vt3dj(k,i,j))
               janjg=(o5*vt3dj(k,i,j)+o6*vt3di(k,i,j))*vt3dj(k,i,j)-3.*rsl*janjc
               janjh=o7*vt3di(k,i,j)+o8*vt3dj(k,i,j)-3.*rsl*janjd
               janji=1-3.*rsl
               janjt1=(-janjh+sqrt(janjh*janjh-4.*janjg*janji))/(2.*janji)
               if (janjt1 > macheps*macheps) then
                  lscal(k)=min(lscal(k),q(k)/sqrt(janjt1))
               end if
            end if
         end do !----- k=k2,m1-1 !

         !---------------------------------------------------------------------------------!
         !   Obtaining an estimation of the PBL height. First, I need to check whether the !
         ! PBL is stable, neutral, or convective. This is obtained by verifying whether    !
         ! z/LMO is > 0, =0 or <0                                                          !
         !---------------------------------------------------------------------------------!
         stable = zeta(k2w)  >=  1.e-4
         neutral = abs(zeta(k2w)) < 1.e-4 
         if (stable .or. neutral) then 
            !------------------------------------------------------------------------------!
            !   If the boundary layer is STABLE or NEUTRAL, then we follow Vogelezang and  !
            ! Holtslag (1996) equation (3); (SBL is defined when (w'theta')g < 0). The PBL !
            ! depth is defined as the first level where average Ri > 0.25.                 !
            !------------------------------------------------------------------------------!
            kpbl(i,j) = k2w
            pblhgt(i,j)=zagl(k2w)
            sboundlay: do k=k2w+1,m1-1
               ri = grav * (thetav0(k)-thetav0(k2w)) * (zagl(k)-zagl(k2w))                 &
                  / ( thetav0(k2w) * ( (uspd(k)-uspd(k2w)) * (uspd(k)-uspd(k2w))    &
                                    + (vspd(k)-vspd(k2w)) * (vspd(k)-vspd(k2w))    &
                                    + 100.*ustarw*ustarw ) )
               if (ri >= 0.25) then
                  kpbl(i,j) = k
                  pblhgt(i,j) = 0.5 * (zagl(k-1)+zagl(k))
                  exit sboundlay
               end if
            end do sboundlay
         else
            !------------------------------------------------------------------------------!
            !    -> Or, if the PBL is convective, then the PBL is defined as the level     !
            !       where the first minimum of w'theta'                                    !
            !------------------------------------------------------------------------------!
            kpbl(i,j) = k2w
            pblhgt(i,j)=zagl(k2w)
            aux= wltl0 * grav / (cp * thetav0(k2w))
            convmixlay: do k=k2w+1,m1-1
               kpbl(i,j) = k
               pblhgt(i,j)=0.5*(zagl(k)+zagl(k-1))
               wstarw=cbrt(aux * pblhgt(i,j) )
               thetavp(k)=thetav0(k2w)+8.5*wltl0/(wstarw*cp)
               ri = grav * (thetav0(k)-thetavp(k)) * (zagl(k)-zagl(k2w))                   &
                  / ( thetavp(k) * ((uspd(k)-uspd(k2w)) * (uspd(k)-uspd(k2w))       &
                    + (vspd(k)-vspd(k2w)) * (vspd(k)-vspd(k2w))                    &
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
            rf=ri1*(ri+ri2-sqrt(ri*ri-ri3*ri+ri2*ri2))
            !------------------------------------------------------------------------------!
            ! Finding the SH and SM of the Level 2 model (Nakanish, 2001):                 !
            !------------------------------------------------------------------------------!
            if (rf < rfc) then
               sh2=3.*a2*(gama1+gama2)*(rfc-rf)/(1.-rf)
               sm2=a1*f1*(rf1-rf)*sh2/(a2*f2*(rf2-rf))
               qlevel2=max(lscal(k)*sqrt(b1*vt3di(k,i,j)*sm2*(1.-rf)),sqrt(2.*tkmin))
            else
               qlevel2=sqrt(2.*tkmin)
            end if
            !------------------------------------------------------------------------------!
            ! Finding Gm and Gh                                                            !
            !------------------------------------------------------------------------------!
            aux=lscal(k)*lscal(k)/qq(k)
            gm(k)=aux*vt3di(k,i,j)       !----- Gm ------------------------------------!
            gh(k)=-aux*vt3dj(k,i,j)      !----- Gh ------------------------------------!

            !------------------------------------------------------------------------------!
            !     Case of growing turbulence, it must use the modified Level 2½ model      !
            ! (Nakanishi and Niino, 2004), otherwise it just uses the original Level 2½.   !
            !------------------------------------------------------------------------------!
            if (q(k) < qlevel2) then
               oneminusalpha=q(k)/max(qlevel2,1.e-10)
               oneminusalpha2=oneminusalpha*oneminusalpha
            else
               oneminusalpha=1.
               oneminusalpha2=1.
            end if
            e1=1.+oneminusalpha2*(ce1a*gm(k)-ce1b*gh(k))
            e2=-oneminusalpha2*ce2*gh(k)
            e3=oneminusalpha2*ce3*gm(k)
            e4=1.-oneminusalpha2*ce4*gh(k)
            r1=oneminusalpha*cr1
            r2=oneminusalpha*a2
            ssm(k)=(r2*e2-r1*e4)/(e2*e3-e1*e4) !----- Sm -------------------------------!
            ssh(k)=(r1*e3-r2*e1)/(e2*e3-e1*e4) !----- Sh -------------------------------!

            !------------------------------------------------------------------------------!
            !    Deriving the convective velocity standard-deviation.                      !
            !    The max was inserted because there is no physical limitation which effec- !
            ! tively prevents the variance to be negative (although in preliminary tests   !
            ! the negative values were significantly smaller than the positive ones.       !
            !------------------------------------------------------------------------------!
            sigw(k,i,j) = sqrt(max(qq(k)*(onethird-2*a1*ssm(k)*gm(k)          &
                        + 4.*a1*(1. - c2)*ssh(k)*gh(k)),sigwmin*sigwmin))

         end do
         !---------------------------------------------------------------------------------!

         !----- Finding the surface weight ------------------------------------------------!
         weightsurf=zagl(k2w)/zagl(k2w+1)

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
               tl(k,i,j)= 0.10 * (pblhgt(i,j)**0.20) * (zagl(k)**0.80) / sigw(k,i,j)
            elseif (neutral) then
            !------------------------------------------------------------------------------!
            !    This is the formula proposed by Hanna (1982), but without the Coriolis    !
            ! term: later papers, such as Wilson (200, JAM), also disconsider this.        !
            !------------------------------------------------------------------------------!
               tl(k,i,j)=0.5*zagl(k)/sigw(k,i,j)
            !------------------------------------------------------------------------------!
            !   For unstable PBLs we should follow equation 7.17, considering all possible !
            ! conditions                                                                   !
            !------------------------------------------------------------------------------!
            elseif (zagl(k) >= 0.1* pblhgt(i,j)) then ! if z/h > 0.1
               ! If PBL is too low, this number can cause underflow in exp. ---------------!

               !----- CHANGED MAX TO MIN (RGK 03-2011). -----------------------------------!
               aux       = min(lnexp_max,5.0*zagl(k)/pblhgt(i,j))
               tl(k,i,j) = 0.15*pblhgt(i,j)*(1-exp(-aux))/sigw(k,i,j)
            elseif ((z0w - zagl(k)) <= lmo(i,j)) then ! -(z-z0)/Lmo > 1 and z/h < 0.1
               tl(k,i,j)= 0.59*zagl(k)/sigw(k,i,j)
            else                                       ! -(z-z0)/Lmo < 1 and z/h < 0.1
               tl(k,i,j)= 0.10*zagl(k)/(sigw(k,i,j)*(0.55+0.38*(zagl(k)-z0w)/lmo(i,j)))
            end if
            tl(k,i,j)=min(tl(k,i,j),ltscalemax) ! Avoiding exaggerated values...
            !------------------------------------------------------------------------------!

            !------------------------------------------------------------------------------!
            ! Finding the TKE tendency.                                                    !
            !------------------------------------------------------------------------------!
            !----- Finding Ps+Pb-epsilon=(q³/L)×(SmGm+ShGh-1/B1) --------------------------!
            pspluspb(k)=qqq(k)*(ssm(k)*gm(k)+ssh(k)*gh(k))/lscal(k)
            epsq(k)=qqq(k)/(b1*lscal(k))
            !------------------------------------------------------------------------------!
            !    Finding the tendency term, based on Nakanishi and Niino equation 5. Note  !
            ! that all the 2 factors vanish because here it's the TKE tendency, and        !
            ! TKE=q²/2.                                                                    !
            !------------------------------------------------------------------------------!
            tket(k,i,j)=tket(k,i,j)+pspluspb(k)-epsq(k)

            !----- Km=L×q×Sm (scr1) -------------------------------------------------------!
            scr1(k,i,j)=lscal(k)*q(k)*ssm(k)*dn0(k,i,j)
            !----- Kh=L×q×Sh (vt3dh) ------------------------------------------------------!
            vt3dh(k,i,j)=lscal(k)*q(k)*ssh(k)*dn0(k,i,j)
            !-----  Kq=L×q×Sq (vt3di). N&N (2004), proposed Sq=2*Sm instead of Sq=0.2 -----!
            vt3di(k,i,j)=2.*scr1(k,i,j) ! Since Km=L×q×Sm... 
         end do

         !---------------------------------------------------------------------------------!
         !    Different closure for the surface, using surface fluxes. If weightsurf =0    !
         ! it'll ignore.                                                                   !
         !---------------------------------------------------------------------------------!
         if (nstbot == 1 .and. weightsurf > 0) then
            tket(k2w,i,j) = tket2 -epsq(k2w) + (1.-weightsurf)*pspluspb(k2w+1)             &
                          + weightsurf*(-sflux_u(i,j)*du0dz-sflux_v(i,j)*dv0dz             &
                          + grav * sflux_t(i,j) / theta(k2w,i,j)) / dn0(k2w,i,j)
         end if
      end do
   end do
   return
end subroutine nakanishi
!==========================================================================================!
!==========================================================================================!
