!===============================================================================
! OLAM version 2.12  -

! Copyright (C) 2002-2006; All Rights Reserved; 
! Duke University, Durham, North Carolina, USA 

! Portions of this software are copied or derived from the RAMS software
! package.  The following copyright notice pertains to RAMS and its derivatives,
! including OLAM:  

   !----------------------------------------------------------------------------
   ! Copyright (C) 1991-2006  ; All Rights Reserved ; Colorado State University; 
   ! Colorado State University Research Foundation ; ATMET, LLC 

   ! This software is free software; you can redistribute it and/or modify it 
   ! under the terms of the GNU General Public License as published by the Free
   ! Software Foundation; either version 2 of the License, or (at your option)
   ! any later version. 

   ! This software is distributed in the hope that it will be useful, but
   ! WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
   ! or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
   ! for more details.
 
   ! You should have received a copy of the GNU General Public License along
   ! with this program; if not, write to the Free Software Foundation, Inc.,
   ! 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA 
   ! (http://www.gnu.org/licenses/gpl.html) 
   !----------------------------------------------------------------------------

! It is requested that any scientific publications based on application of OLAM
! include the following acknowledgment:  "OLAM was developed at the 
! Edmund T. Pratt Jr. School of Engineering, Duke University."

! For additional information, including published references, please contact
! the software authors, Robert L. Walko (robert.walko@duke.edu)
! or Roni Avissar (avissar@duke.edu).
!===============================================================================
subroutine harr_swrad(nrad,albedt,cosz,time,   &
                      u,pl,tl,dzl,vp,tp,omgp,gp,fu,fd,flxus,flxds,ngass,mynum)

  ! Bob's interface subroutine

  use mem_harr, only: mg, ng, mb, nb, nsolb, npsb, xp, alpha, beta, wght,  &
                      prf,trf,ralcs, solar1, ulim

  implicit none

  integer, intent(in) :: nrad
  integer, intent(in) :: mynum
  integer, intent(in) :: ngass(mg)

  real, intent(in) :: albedt
  real, intent(in) :: cosz
  real, intent(in) :: time
  real, intent(in) :: u   (nrad,3)
  real, intent(in) :: pl  (nrad)
  real, intent(in) :: tl  (nrad)
  real, intent(in) :: dzl (nrad)
  real, intent(in) :: vp  (nrad)
  real, intent(in) :: tp  (nrad,mb)
  real, intent(in) :: omgp(nrad,mb)
  real, intent(in) :: gp  (nrad,mb)
  real, intent(in) :: fu  (nrad,6)
  real, intent(in) :: fd  (nrad,6)

  real, intent(inout) :: flxus(nrad)
  real, intent(inout) :: flxds(nrad)

  call swrad(nrad,ng,nb,nsolb,npsb,     &
     u,pl,tl,dzl,vp,                    &
     xp,alpha,beta,wght,prf,trf,ralcs,  &
     solar1,ngass,                      &
     albedt,cosz,                       &
     tp,omgp,gp,fu,fd,flxus,flxds,ulim,time,mynum)

  return
end subroutine harr_swrad

!-----------------------------------------------------------------------------------------!
subroutine harr_lwrad(nrad,u,pl,tl,dzl,vp,tp,omgp,gp,fu,fd,flxul,flxdl,ngast,mynum)

  ! Bob's interface subroutine

  use mem_harr, only: mg, ng, mb, nb, nsolb, npsb, nuum, xp, alpha, beta, wght,  &
                      prf, trf, ralcs, a0, a1, a2, a3, exptabc, ulim

  implicit none

  integer, intent(in) :: nrad
  integer, intent(in) :: mynum
  integer, intent(in) :: ngast(mg)

  real, intent(in) :: u   (nrad,3)
  real, intent(in) :: pl  (nrad)
  real, intent(in) :: tl  (nrad)
  real, intent(in) :: dzl (nrad)
  real, intent(in) :: vp  (nrad)
  real, intent(in) :: tp  (nrad,mb)
  real, intent(in) :: omgp(nrad,mb)
  real, intent(in) :: gp  (nrad,mb)
  real, intent(in) :: fu  (nrad,6)
  real, intent(in) :: fd  (nrad,6)

  real, intent(inout) :: flxul(nrad)
  real, intent(inout) :: flxdl(nrad)

  call lwrad(nrad,ng,nb,nsolb,npsb,nuum,   &
     u,pl,tl,dzl,vp,                       &
     xp,alpha,beta,wght,prf,trf,ralcs,     &
     a0,a1,a2,a3,                          &
     exptabc,ngast,                        &
     tp,omgp,gp,fu,fd,flxul,flxdl,ulim,mynum)
  return
end subroutine harr_lwrad


!
! ---------------------------------------------------------------------
!
! * New swrad (Jerry, March 8, 1996)

subroutine swrad(nz,ng,nb,ns,npsb,               &
           u,pl,tl,dz,vp,                        &
           xp,alpha,beta,wght,prf,trf,ral,       &
           solar,ngas,                           &
           alb,amu0,                             &
           tp,omgp,asym,fu,fd,flxsu,flxsd,ulim,time,mynum)

  use rconstants, only : cp
  use mem_harr, only : mb,mg,mk
  implicit none
  integer :: nz,ng,nb,ns,mynum
  real :: time

  !
  !     two-stream radiative transfer code.
  !     This version will attempt the use of FESFT
  !     as outlined in R&G (MWR 1991)
  !      JH 5-95
  !
  real, parameter :: top=1800.,tm=1800./296.,gma=0.002

  !     input

  real    :: u(nz,mg),pl(nz),tl(nz),dz(nz),vp(nz)
  real    :: ulim(mg,mb)
  integer :: npsb(mg,mb),na(mg),ngas(mg)

  !     input parameters

  real    :: ral(mb),wlenlo(mb),wlenhi(mb)
  real    :: xp(mg,mk,mb),alpha(mg,mk,mb),beta(mg,mk,mb)
  real    :: wght(mg,mk,mb),prf(mg,mb),trf(mg,mb)
  real    :: solar(mb)

  !     arrays used in fluxes

  real    :: tg(nz),tp(nz,nb),tcr(nz),omgp(nz,nb)
  real    :: alb,amu0,asym(nz,nb)
  real    :: t(nz),r(nz),tc(nz),sigu(nz),sigd(nz)
  real    :: re(nz),vd(nz),td(nz),vu(nz)
  real    :: fu(nz,6),fd(nz,6),flxsu(nz),flxsd(nz)

  integer ib,iz,mxa,ig,ik,ig1,ig2,ik1,ik2,ig3,ik3
  real :: uu,fact,ufac,dfac,tu1,tu2,td1,td2,fbu,fbd,tu3,td3

  !     output

  !     remember, for this code it is assumed that u, the gaseous
  !     absorber amounts, are specified in Pa
  !----------------------------------------------------------------------------
  !     zero out flxsu and flxsd

  !     do iz = 1,nz
  !      flxsu(iz) = 0.0
  !      flxsd(iz) = 0.0
  !     enddo

  !     loop through each band

  bandloop: do ib=1,ns

  !        calculate here the properties that are considered grey,
  !        i.e. averaged values are used across the band...
  !        rayleigh scatter...cloud is now done outside of this routine

  !           get rayleigh scattering AND
  !           zero out the local flux arrays

       do iz=1,nz
           tcr(iz) = ral(ib)*pl(iz)*dz(iz)/tl(iz)
           fu(iz,2) = 0.
           fu(iz,3) = 0.
           fu(iz,4) = 0.
           fd(iz,2) = 0.
           fd(iz,3) = 0.
           fd(iz,4) = 0.
       enddo

  !        determine if, and how many overlaps...also check
  !        to see if the gas is used

     mxa = 0
     do ig=1,mg
        if (npsb(ig,ib) > 0 .and. ngas(ig) == 1)then
           mxa = mxa+1
           na(mxa) = ig
        endif
     enddo

     if (mxa == 0) cycle bandloop
     if (mxa == 1) then

  !-------------------------------------------------------------------------
  !           no overlapping gasses, single aborber

        ig = na(1)

        do ik=1,npsb(ig,ib)

            do iz=nz,2,-1
              uu=min(ulim(ig,ib),u(iz,ig))
              tg(iz) = xp(ig,ik,ib)*uu*  &
                 (pl(iz)/prf(ig,ib))**alpha(ig,ik,ib)*  &
                 (trf(ig,ib)/tl(iz))**beta(ig,ik,ib)
            enddo

  !               now do rest of stuff in subroutine

            call flxsw(nz,tg,tp(1:nz,ib),tcr,omgp(1:nz,ib),  &
                    alb,solar(ib),amu0,t,r,tc,  &
                    sigu,sigd,re,vd,td,vu,  &
                    fu(1:nz,6),fd(1:nz,6),asym(1:nz,ib),time,mynum)

  !               add pseudo-band fluxes to the total flux

            do iz=1,nz
              flxsu(iz)=flxsu(iz)+wght(ig,ik,ib)*fu(iz,6)
              flxsd(iz)=flxsd(iz)+wght(ig,ik,ib)*fd(iz,6)
            enddo

        enddo
     else if (mxa == 2) then

  !--------------------------------------------------------------------------
  !           overlap of two gasses using the FESFT

        ig1 = na(1)
        ig2 = na(2)

  !           do the gray fluxes first

            tg(1:nz) = 0.

            call flxsw(nz,tg,tp(1:nz,ib),tcr,omgp(1:nz,ib),  &
                    alb,solar(ib),amu0,t,r,tc,  &
                    sigu,sigd,re,vd,td,vu,  &
                    fu(1:nz,1),fd(1:nz,1),asym(1:nz,ib),time,mynum)

  !           do the 1st gas


        do ik1=1,npsb(ig1,ib)

            do iz=2,nz
              uu=min(ulim(ig1,ib),u(iz,ig1))
              tg(iz) = xp(ig1,ik1,ib)*uu*  &
                 (pl(iz)/prf(ig1,ib))**alpha(ig1,ik1,ib)*  &
                 (trf(ig1,ib)/tl(iz))**beta(ig1,ik1,ib)
            enddo

            call flxsw(nz,tg,tp(1:nz,ib),tcr,omgp(1:nz,ib),  &
                    alb,solar(ib),amu0,t,r,tc,  &
                    sigu,sigd,re,vd,td,vu,  &
                    fu(1:nz,6),fd(1:nz,6),asym(1:nz,ib),time,mynum)

            fact = wght(ig1,ik1,ib)
            do iz=1,nz
              fu(iz,2) = fu(iz,2)+fact*fu(iz,6)
              fd(iz,2) = fd(iz,2)+fact*fd(iz,6)
            enddo

         enddo

  !            do the 2nd gas

         do ik2=1,npsb(ig2,ib)

            do iz=1,nz
              uu=min(ulim(ig2,ib),u(iz,ig2))
              tg(iz) = xp(ig2,ik2,ib)*uu*  &
                 (pl(iz)/prf(ig2,ib))**alpha(ig2,ik2,ib)*  &
                 (trf(ig2,ib)/tl(iz))**beta(ig2,ik2,ib)
            enddo

            call flxsw(nz,tg,tp(1:nz,ib),tcr,omgp(1:nz,ib),  &
                    alb,solar(ib),amu0,t,r,tc,  &
                    sigu,sigd,re,vd,td,vu,  &
                    fu(1:nz,6),fd(1:nz,6),asym(1:nz,ib),time,mynum)

            fact = wght(ig2,ik2,ib)
            do iz=1,nz
              fu(iz,3) = fu(iz,3)+fact*fu(iz,6)
              fd(iz,3) = fd(iz,3)+fact*fd(iz,6)
            enddo
        enddo

  ! this is the stuff that is fixed (old stuff deleted)

        do iz=1,nz
          ufac=max(1.e-14,fu(iz,1))
          dfac=max(1.e-14,fd(iz,1))
          tu1 = max(0.0,min(1.1,fu(iz,2)/ufac))
          tu2 = max(0.0,min(1.1,fu(iz,3)/ufac))
          td1 = max(0.0,min(1.1,fd(iz,2)/dfac))
          td2 = max(0.0,min(1.1,fd(iz,3)/dfac))
          fbu = tu1*tu2*fu(iz,1)
          fbd = td1*td2*fd(iz,1)

          flxsu(iz)=flxsu(iz)+fbu
          flxsd(iz)=flxsd(iz)+fbd
        enddo


        else if (mxa == 3) then

  !--------------------------------------------------------------------------
  !           overlap of three gasses

        ig1 = na(1)
        ig2 = na(2)
        ig3 = na(3)

  !           do the gray fluxes first

            tg(1:nz) = 0.

            call flxsw(nz,tg,tp(1:nz,ib),tcr,omgp(1:nz,ib),  &
                    alb,solar(ib),amu0,t,r,tc,  &
                    sigu,sigd,re,vd,td,vu,  &
                    fu(1:nz,1),fd(1:nz,1),asym(1:nz,ib),time,mynum)

  !           do the 1st gas


        do ik1=1,npsb(ig1,ib)

            do iz=2,nz
              uu=min(ulim(ig1,ib),u(iz,ig1))
              tg(iz) = xp(ig1,ik1,ib)*uu*  &
                 (pl(iz)/prf(ig1,ib))**alpha(ig1,ik1,ib)*  &
                 (trf(ig1,ib)/tl(iz))**beta(ig1,ik1,ib)
!              write(unit=60+mynum,fmt='(3(a,1x,i3,1x),9(a,1x,es10.3,1x))') '1st gas... iz=',iz,'ik1=',ik1,'ib=',ib &
!                                         ,'xp=',xp(ig1,ik1,ib),'uu=',uu,'pl=',pl(iz) &
!                                         ,'prf=',prf(ig1,ib),'alpha',alpha(ig1,ik1,ib) & 
!                                         ,'trf=',trf(ig1,ib),'tl=',tl(iz),'beta=',beta(ig1,ik1,ib),'tg=',tg(iz)
            enddo

            call flxsw(nz,tg,tp(1:nz,ib),tcr,omgp(1:nz,ib),  &
                    alb,solar(ib),amu0,t,r,tc,  &
                    sigu,sigd,re,vd,td,vu,  &
                    fu(1:nz,6),fd(1:nz,6),asym(1:nz,ib),time,mynum)

            fact = wght(ig1,ik1,ib)
            do iz=1,nz
              fu(iz,2) = fu(iz,2)+fact*fu(iz,6)
              fd(iz,2) = fd(iz,2)+fact*fd(iz,6)
            enddo
        enddo

  !           do the 2nd gas

        do ik2=1,npsb(ig2,ib)

            do iz=2,nz
              uu=min(ulim(ig2,ib),u(iz,ig2))
              tg(iz) = xp(ig2,ik2,ib)*uu*  &
                 (pl(iz)/prf(ig2,ib))**alpha(ig2,ik2,ib)*  &
                 (trf(ig2,ib)/tl(iz))**beta(ig2,ik2,ib)
!              write(unit=70+mynum,fmt='(3(a,1x,i3,1x),9(a,1x,es10.3,1x))') '2nd gas... iz=',iz,'ik1=',ik1,'ib=',ib &
!                                         ,'xp=',xp(ig1,ik1,ib),'uu=',uu,'pl=',pl(iz) &
!                                         ,'prf=',prf(ig1,ib),'alpha',alpha(ig1,ik1,ib) & 
!                                         ,'trf=',trf(ig1,ib),'tl=',tl(iz),'beta=',beta(ig1,ik1,ib),'tg=',tg(iz)
            enddo

            call flxsw(nz,tg,tp(1:nz,ib),tcr,omgp(1:nz,ib),  &
                    alb,solar(ib),amu0,t,r,tc,  &
                    sigu,sigd,re,vd,td,vu,  &
                    fu(1:nz,6),fd(1:nz,6),asym(1:nz,ib),time,mynum)

            fact = wght(ig2,ik2,ib)
            do iz=1,nz
              fu(iz,3) = fu(iz,3)+fact*fu(iz,6)
              fd(iz,3) = fd(iz,3)+fact*fd(iz,6)
            enddo
          enddo

  !          do the 3rd gas

        do ik3=1,npsb(ig3,ib)

            do iz=2,nz
              uu=min(ulim(ig3,ib),u(iz,ig3))
              tg(iz) = xp(ig3,ik3,ib)*uu*  &
                 (pl(iz)/prf(ig3,ib))**alpha(ig3,ik3,ib)*  &
                 (trf(ig3,ib)/tl(iz))**beta(ig3,ik3,ib)
!              write(unit=80+mynum,fmt='(3(a,1x,i3,1x),9(a,1x,es10.3,1x))') '3rd gas... iz=',iz,'ik1=',ik1,'ib=',ib &
!                                         ,'xp=',xp(ig1,ik1,ib),'uu=',uu,'pl=',pl(iz) &
!                                         ,'prf=',prf(ig1,ib),'alpha',alpha(ig1,ik1,ib) & 
!                                         ,'trf=',trf(ig1,ib),'tl=',tl(iz),'beta=',beta(ig1,ik1,ib),'tg=',tg(iz)
            enddo

  !               now do rest of stuff in subroutine

            call flxsw(nz,tg,tp(1:nz,ib),tcr,omgp(1:nz,ib),  &
                    alb,solar(ib),amu0,t,r,tc,  &
                    sigu,sigd,re,vd,td,vu,  &
                    fu(1:nz,6),fd(1:nz,6),asym(1:nz,ib),time,mynum)

  !               sum the pseudo-band fluxes to get total flux

            fact = wght(ig3,ik3,ib)
            do iz=1,nz
              fu(iz,4) = fu(iz,4)+fact*fu(iz,6)
              fd(iz,4) = fd(iz,4)+fact*fd(iz,6)
            enddo
        enddo

  ! this is the stuff that is fixed (old stuff deleted)

        do iz=1,nz
          ufac=max(1.e-14,fu(iz,1))
          dfac=max(1.e-14,fd(iz,1))
          tu1 = max(0.0,min(1.1,fu(iz,2)/ufac))
          tu2 = max(0.0,min(1.1,fu(iz,3)/ufac))
          tu3 = max(0.0,min(1.1,fu(iz,4)/ufac))
          td1 = max(0.0,min(1.1,fd(iz,2)/dfac))
          td2 = max(0.0,min(1.1,fd(iz,3)/dfac))
          td3 = max(0.0,min(1.1,fd(iz,4)/dfac))

          fbu = tu1*tu2*tu3*fu(iz,1)
          fbd = td1*td2*td3*fd(iz,1)
          flxsu(iz)=flxsu(iz)+fbu
          flxsd(iz)=flxsd(iz)+fbd
        end do
  !-------------------------------------------------------------------------

     else
          write(*,'(a)') '╗╗╗ ERROR ллллллллллллллллллллллллллллллллллллллл'
          write(*,'(a)') 'Two is the maximum amount of overlapping'
          write(*,'(a)') 'gasses allowed: if you really want to have '
          write(*,'(a)') 'more gasses overlap, they better have few'
          write(*,'(a)') 'pseudo-bands for all of them, and in any'
          write(*,'(a)') 'case it will cost you. To do that come into'
          write(*,'(a)') 'the code and modify it here, look at previous'
          write(*,'(a)') 'structures to see how it is done. Its easy'
          write(*,'(a)') 'but BEWARE: ITS BOUND TO BE HORRIBLY EXPENSIVE'
          stop '╗╗╗ swrad (harr_rad.f90)'
     endif
!  111     continue

  end do bandloop
  !
  return
end subroutine swrad


!---------------------------------------------------------------------------------------!
! New lwrad (March 8, 1996)

subroutine lwrad(nz,ng,nb,ns,npsb,nuum,            &
           u,pl,tl,dz,vp,                          &
           xp,alpha,beta,wght,prf,trf,ral,         &
           a0,a1,a2,a3,                            &
           exptabc,ngas,                           &
           tp,omgp,asym,fu,fd,flxlu,flxld,ulim,mynum)

  use rconstants, only: cp
  use mem_harr, only: mb,mg,mk
  implicit none

  real, parameter :: top=1800.,tm=1800./296.,gma=0.002,tr=296.,ccc=6.08
  integer :: nz,ng,nb,ns,mynum
  !
  !     two-stream radiative transfer code.
  !     This version will attempt the use of FESFT
  !     as outlined in R&G (MWR 1991)
  !      JH 5-95
  !
  !     input

  real    u(nz,mg),pl(nz),tl(nz),dz(nz),vp(nz)
  real    exptabc(150),ulim(mg,mb)
  integer npsb(mg,mb),na(mg),nuum(mb),ngas(mg)

  !     input parameters

  real    ral(mb),wlenlo(mb),wlenhi(mb)
  real    xp(mg,mk,mb),alpha(mg,mk,mb),beta(mg,mk,mb)
  real    wght(mg,mk,mb),prf(mg,mb),trf(mg,mb)
  real    a0(mb),a1(mb),a2(mb),a3(mb)

  !     arrays used in fluxes

  real    tg(nz),tp(nz,nb),tcr(nz),omgp(nz,nb),src(nz)
  real    t(nz),r(nz),tc(nz),sigu(nz),sigd(nz)
  real    re(nz),vd(nz),td(nz),vu(nz),asym(nz,nb)
  real    fu(nz,6),fd(nz,6),flxlu(nz),flxld(nz)

  integer :: ib,iflag,iz,mxa,ig,ik,ig1,ig2,ik1,ik2,nir,ii,ig3,ik3
  real :: tf,ewght,expp,uu,chck,fact,ufac,tu1,tu2,fbu,dx         &
         ,dxc,tn1c,tn2c,tn1,tn2,difflxb,fbd,dfac,tu3,tn3c,tn3

  !     remember, for this code it is assumed that u, the gaseous
  !     absorber amounts, are specified in g/m^3
  !----------------------------------------------------------------------------
  !     zero out flxsu and flxsd

  !     do iz = 1,nz
  !      flxlu(iz) = 0.0
  !      flxld(iz) = 0.0
  !     enddo

  !     loop through each band

  nir=ns+1
  bandloop: do ib=nir,nb

  !        calculate the properties that are grey across the band
  !        Planck function at the interfaces, continuum absorption,
  !        ...cloud is now done outside of this routine

        iflag=0

  !           do surface and top of atmosphere first

       src(1) = a0(ib)+tl(1)*(a1(ib)+tl(1)*  &
                       (a2(ib)+a3(ib)*tl(1)))
       src(2) = a0(ib)+tl(2)*(a1(ib)+tl(2)*  &
                       (a2(ib)+a3(ib)*tl(2)))
       src(nz) =  a0(ib)+tl(nz)*(a1(ib)+tl(nz)*  &
                           (a2(ib)+a3(ib)*tl(nz)))

       tcr(1)=0.
       tcr(2)=0.
       tcr(nz)=0.

       do iz=nz-1,3,-1

  !                 get sources at the interface, temp first

              tf = 0.5*(tl(iz)+tl(iz-1))
              src(iz) = a0(ib)+tf*(a1(ib)+tf*  &
                           (a2(ib)+a3(ib)*tf))
              tcr(iz) = 0.
        enddo


        if (nuum(ib).eq.1) then

  !              this band has continuum

           do iz=2,nz
               ii=int(tl(iz)-179.)
               ewght=tl(iz)-(float(ii)+179.)
               if(ii > 0)then
                  expp=(1.-ewght)*exptabc(ii)+ewght*exptabc(ii+1)
               else
                  expp=exp(1800.0/tl(iz)-1800.0/296.0)
               endif
!MLO - Switching by what is given in J. Harrington's dissertation
!               expp=exp(ccc*(tr/tl(iz) - 1))
               tcr(iz) = ral(ib)*u(iz,1)*  &
                      expp*(vp(iz)+  &
                      gma*(pl(iz)-vp(iz)))
!               write (unit=60+mynum,fmt='(2(a,1x,i5,1x),7(a,1x,es10.3,1x))') 'iz=',iz,'ib=',ib &
!                          ,'ral=',ral(ib),'u=',u(iz,1),'expp=',expp &
!                          ,'tl=',tl(iz),'vp=',vp(iz),'pl=',pl(iz),'tcr=',tcr(iz)
           enddo

        endif


     do iz=1,nz
           fu(iz,2) = 0.
           fu(iz,3) = 0.
           fu(iz,4) = 0.
           fd(iz,2) = 0.
           fd(iz,3) = 0.
           fd(iz,4) = 0.
           tg(iz)=0.
     enddo

  !        determine if, and how many overlaps

     mxa = 0
     do ig=1,mg
        if (npsb(ig,ib).gt.0.and.ngas(ig).eq.1)then
           mxa = mxa+1
           na(mxa) = ig
        endif
     enddo

     if (mxa.eq.0) cycle bandloop
     if (mxa.eq.1) then

  !-------------------------------------------------------------------------
  !           no overlapping gasses, single aborber

        ig = na(1)

        do ik=1,npsb(ig,ib)

            do iz=nz,2,-1
              uu=min(ulim(ig,ib),u(iz,ig))
              tg(iz) = xp(ig,ik,ib)*uu*  &
                 (pl(iz)/prf(ig,ib))**alpha(ig,ik,ib)*  &
                 (trf(ig,ib)/tl(iz))**beta(ig,ik,ib)
            enddo

  !               now do rest of stuff in subroutine

            call flxlw(nz,tg,tp(1:nz,ib),tcr,omgp(1:nz,ib),  &
                    src,t,r,tc,  &
                    sigu,sigd,re,vd,td,vu,  &
                    fu(1:nz,6),fd(1:nz,6),1.,asym(1:nz,ib),mynum)

  !               add pseudo-band fluxes to the total flux

            do iz=1,nz
              flxlu(iz)=flxlu(iz)+wght(ig,ik,ib)*fu(iz,6)
              flxld(iz)=flxld(iz)+wght(ig,ik,ib)*fd(iz,6)
            enddo
        enddo
     else if (mxa.eq.2) then

  !--------------------------------------------------------------------------
  !           overlap of two gasses using the FESFT

        ig1 = na(1)
        ig2 = na(2)

  !           do the gray fluxes first

        chck=float(iflag)
        
            tg(1:nz) = 0.
        
            call flxlw(nz,tg,tp(1:nz,ib),tcr,omgp(1:nz,ib),  &
                    src,t,r,tc,  &
                    sigu,sigd,re,vd,td,vu,  &
                    fu(1:nz,1),fd(1:nz,1),chck,asym(1:nz,ib),mynum)

  !           if there is continuum abs. do it now since tg=0
  !
       if(nuum(ib).eq.1)then

            call flxlw(nz,tg,tp(1:nz,ib),tcr,omgp(1:nz,ib),  &
                    src,t,r,tc,  &
                    sigu,sigd,re,vd,td,vu,  &
                    fu(1:nz,5),fd(1:nz,5),1.,asym(1:nz,ib),mynum)
       endif

  !           do the 1st gas


        do ik1=1,npsb(ig1,ib)

            do iz=2,nz
              uu=min(ulim(ig1,ib),u(iz,ig1))
              tg(iz) = xp(ig1,ik1,ib)*uu*  &
                 (pl(iz)/prf(ig1,ib))**alpha(ig1,ik1,ib)*  &
                 (trf(ig1,ib)/tl(iz))**beta(ig1,ik1,ib)
            enddo

            call flxlw(nz,tg,tp(1:nz,ib),tcr,omgp(1:nz,ib),  &
                    src,t,r,tc,  &
                    sigu,sigd,re,vd,td,vu,  &
                    fu(1:nz,6),fd(1:nz,6),chck,asym(1:nz,ib),mynum)

            fact = wght(ig1,ik1,ib)
            do iz=1,nz
              fu(iz,2) = fu(iz,2)+fact*fu(iz,6)
              fd(iz,2) = fd(iz,2)+fact*fd(iz,6)
            enddo
         enddo

  !            do the 2nd gas

         do ik2=1,npsb(ig2,ib)

            do iz=1,nz
              uu=min(ulim(ig2,ib),u(iz,ig2))
              tg(iz) = xp(ig2,ik2,ib)*uu*  &
                 (pl(iz)/prf(ig2,ib))**alpha(ig2,ik2,ib)*  &
                 (trf(ig2,ib)/tl(iz))**beta(ig2,ik2,ib)
            enddo

            call flxlw(nz,tg,tp(1:nz,ib),tcr,omgp(1:nz,ib),  &
                    src,t,r,tc,  &
                    sigu,sigd,re,vd,td,vu,  &
                    fu(1:nz,6),fd(1:nz,6),chck,asym(1:nz,ib),mynum)

            fact = wght(ig2,ik2,ib)
            do iz=1,nz
              fu(iz,3) = fu(iz,3)+fact*fu(iz,6)
              fd(iz,3) = fd(iz,3)+fact*fd(iz,6)
            enddo
        enddo


  !  above (now wiped out...) is old section, below is the new bounded part:

        do iz=1,nz
          ufac=max(1.e-14,fu(iz,1))
          tu1 = max(0.0,min(1.1,fu(iz,2)/ufac))
          tu2 = max(0.0,min(1.1,fu(iz,3)/ufac))
          fbu = tu1*tu2*fu(iz,1)

          dx=fd(iz,1)-fu(iz,1)
          dxc = sign(max(1.e-14,abs(dx)),dx)
          tn1c = (fd(iz,2)-fu(iz,2))/dxc + 1.
          tn2c = (fd(iz,3)-fu(iz,3))/dxc + 1.
          tn1 = min(max(tn1c,0.),2.) - 1.
          tn2 = min(max(tn2c,0.),2.) - 1.
          difflxb = tn1*tn2*dx
          fbd=difflxb+fbu

          flxlu(iz)=flxlu(iz)+fbu
          flxld(iz)=flxld(iz)+fbd
        enddo


        else if (mxa.eq.3) then

  !--------------------------------------------------------------------------
  !           overlap of three gasses

        ig1 = na(1)
        ig2 = na(2)
        ig3 = na(3)

  !           do the gray fluxes first

        chck=float(iflag)

            tg(1:nz) = 0.
        
            call flxlw(nz,tg,tp(1:nz,ib),tcr,omgp(1:nz,ib),  &
                    src,t,r,tc,  &
                    sigu,sigd,re,vd,td,vu,  &
                    fu(1:nz,1),fd(1:nz,1),chck,asym(1:nz,ib),mynum)

  !           if there is continuum abs. do it now since tg=0
  !
       if(nuum(ib).eq.1)then

            call flxlw(nz,tg,tp(1:nz,ib),tcr,omgp(1:nz,ib),  &
                    src,t,r,tc,  &
                    sigu,sigd,re,vd,td,vu,  &
                    fu(1:nz,5),fd(1:nz,5),1.,asym(1:nz,ib),mynum)
       endif

  !           do the 1st gas


        do ik1=1,npsb(ig1,ib)

            do iz=2,nz
              uu=min(ulim(ig1,ib),u(iz,ig1))
              tg(iz) = xp(ig1,ik1,ib)*uu*  &
                 (pl(iz)/prf(ig1,ib))**alpha(ig1,ik1,ib)*  &
                 (trf(ig1,ib)/tl(iz))**beta(ig1,ik1,ib)
!              write(unit=60+mynum,fmt='(3(a,1x,i3,1x),9(a,1x,es10.3,1x))') '1st gas... iz=',iz,'ik1=',ik1,'ib=',ib &
!                                         ,'xp=',xp(ig1,ik1,ib),'uu=',uu,'pl=',pl(iz) &
!                                         ,'prf=',prf(ig1,ib),'alpha',alpha(ig1,ik1,ib) & 
!                                         ,'trf=',trf(ig1,ib),'tl=',tl(iz),'beta=',beta(ig1,ik1,ib),'tg=',tg(iz)
            enddo

            call flxlw(nz,tg,tp(1:nz,ib),tcr,omgp(1:nz,ib),  &
                    src,t,r,tc,  &
                    sigu,sigd,re,vd,td,vu,  &
                    fu(1:nz,6),fd(1:nz,6),chck,asym(1:nz,ib),mynum)

            fact = wght(ig1,ik1,ib)
            do iz=1,nz
              fu(iz,2) = fu(iz,2)+fact*fu(iz,6)
              fd(iz,2) = fd(iz,2)+fact*fd(iz,6)
            enddo
          enddo

  !           do the 2nd gas

        do ik2=1,npsb(ig2,ib)

            do iz=2,nz
              uu=min(ulim(ig2,ib),u(iz,ig2))
              tg(iz) = xp(ig2,ik2,ib)*uu*  &
                 (pl(iz)/prf(ig2,ib))**alpha(ig2,ik2,ib)*  &
                 (trf(ig2,ib)/tl(iz))**beta(ig2,ik2,ib)
!              write(unit=70+mynum,fmt='(3(a,1x,i3,1x),9(a,1x,es10.3,1x))') '2nd gas... iz=',iz,'ik1=',ik1,'ib=',ib &
!                                         ,'xp=',xp(ig1,ik1,ib),'uu=',uu,'pl=',pl(iz) &
!                                         ,'prf=',prf(ig1,ib),'alpha',alpha(ig1,ik1,ib) & 
!                                         ,'trf=',trf(ig1,ib),'tl=',tl(iz),'beta=',beta(ig1,ik1,ib),'tg=',tg(iz)
            enddo

            call flxlw(nz,tg,tp(1:nz,ib),tcr,omgp(1:nz,ib),  &
                    src,t,r,tc,  &
                    sigu,sigd,re,vd,td,vu,  &
                    fu(1:nz,6),fd(1:nz,6),chck,asym(1:nz,ib),mynum)

            fact = wght(ig2,ik2,ib)
            do iz=1,nz
              fu(iz,3) = fu(iz,3)+fact*fu(iz,6)
              fd(iz,3) = fd(iz,3)+fact*fd(iz,6)
            enddo
          enddo

  !          do the 3rd gas

        do ik3=1,npsb(ig3,ib)

            do iz=2,nz
              uu=min(ulim(ig3,ib),u(iz,ig3))
              tg(iz) = xp(ig3,ik3,ib)*uu*  &
                 (pl(iz)/prf(ig3,ib))**alpha(ig3,ik3,ib)*  &
                 (trf(ig3,ib)/tl(iz))**beta(ig3,ik3,ib)
!              write(unit=80+mynum,fmt='(3(a,1x,i3,1x),9(a,1x,es10.3,1x))') '3rd gas... iz=',iz,'ik1=',ik1,'ib=',ib &
!                                         ,'xp=',xp(ig1,ik1,ib),'uu=',uu,'pl=',pl(iz) &
!                                         ,'prf=',prf(ig1,ib),'alpha',alpha(ig1,ik1,ib) & 
!                                         ,'trf=',trf(ig1,ib),'tl=',tl(iz),'beta=',beta(ig1,ik1,ib),'tg=',tg(iz)
             enddo

  !               now do rest of stuff in subroutine

            call flxlw(nz,tg,tp(1:nz,ib),tcr,omgp(1:nz,ib),  &
                    src,t,r,tc,  &
                    sigu,sigd,re,vd,td,vu,  &
                    fu(1:nz,6),fd(1:nz,6),chck,asym(1:nz,ib),mynum)

  !               sum the pseudo-band fluxes to get total flux

            fact = wght(ig3,ik3,ib)
            do iz=1,nz
              fu(iz,4) = fu(iz,4)+fact*fu(iz,6)
              fd(iz,4) = fd(iz,4)+fact*fd(iz,6)
            enddo
        enddo


        do iz=1,nz
          ufac=max(1.e-14,fu(iz,1))
          dfac=max(1.e-14,fd(iz,1))
          tu1 = max(0.0,min(1.1,fu(iz,2)/ufac))
          tu2 = max(0.0,min(1.1,fu(iz,3)/ufac))
          tu3 = max(0.0,min(1.1,fu(iz,4)/ufac))
          fbu = tu1*tu2*tu3*fu(iz,1)

          dx=fd(iz,1)-fu(iz,1)
          dxc = sign(max(1.e-14,abs(dx)),dx)
          tn1c = (fd(iz,2)-fu(iz,2))/dxc + 1.
          tn2c = (fd(iz,3)-fu(iz,3))/dxc + 1.
          tn3c = (fd(iz,4)-fu(iz,4))/dxc + 1.
          tn1 = min(max(tn1c,0.),2.) - 1.
          tn2 = min(max(tn2c,0.),2.) - 1.
          tn3 = min(max(tn3c,0.),2.) - 1.
          difflxb = tn1*tn2*tn3*dx
          fbd=difflxb+fbu

          tn1 = max(0.0,min(1.1,fd(iz,2)/dfac))
          tn2 = max(0.0,min(1.1,fd(iz,3)/dfac))
          tn3 = max(0.0,min(1.1,fd(iz,4)/dfac))
          fbd = tn1*tn2*tn3*fd(iz,1)


          flxlu(iz)=flxlu(iz)+fbu
          flxld(iz)=flxld(iz)+fbd
        enddo



  !-------------------------------------------------------------------------

     else
          write(*,'(a)') '╗╗╗ ERROR ллллллллллллллллллллллллллллллллллллллл'
          write(*,'(a)') 'Two is the maximum amount of overlapping'
          write(*,'(a)') 'gasses allowed: if you really want to have '
          write(*,'(a)') 'more gasses overlap, they better have few'
          write(*,'(a)') 'pseudo-bands for all of them, and in any'
          write(*,'(a)') 'case it will cost you. To do that come into'
          write(*,'(a)') 'the code and modify it here, look at previous'
          write(*,'(a)') 'structures to see how it is done. Its easy'
          write(*,'(a)') 'but BEWARE: ITS BOUND TO BE HORRIBLY EXPENSIVE'
          stop '╗╗╗ lwrad (harr_rad.f90)'
     endif
  end do bandloop

  return
end subroutine lwrad



!
! ---------------------------------------------------------------------
!
subroutine flxsw(nz,tg,tp,tcr,omgp,alb,slr,amu0,t,r,tc,sigu,sigd,re,vd  &
                 ,td,vu,fu,fd,asym,time,mynum)
  implicit none
  integer :: nz,mynum
  real, intent(in) , dimension(nz) :: tg,tp,tcr,omgp,asym
  real, intent(in)                 :: alb,slr,amu0,time
  real, intent(out), dimension(nz) :: t,r,tc,sigu,sigd,re,vd,td,vu
  real, intent(out), dimension(nz) :: fu,fd

  ! Local variables
  integer :: iz
  ! Double precision version of I/O variables
  real(kind=8)     ,  dimension(nz) :: tg8,tp8,tcr8,omgp8,asym8
  real(kind=8)                      :: alb8,slr8,amu08
  real(kind=8)     ,  dimension(nz) :: t8,r8,tc8,sigu8,sigd8,re8,vd8,td8,vu8
  real(kind=8)     ,  dimension(nz) :: fu8,fd8
  !
  real(kind=8)                      :: expt1
  real(kind=8) :: tau,omg0,af,fact,beta0,g1,g2,gg,rinf,ggtau,expp1,expp2,denomi  &
                 ,cc,g3,g4,aa,bb,tcm2,exp1,exp2,amu0i
  real(kind=8), parameter :: diffac=1.66, tinyreal=1.d-20,eps=1.d-6,flush=46.
  real        , external  :: dble2sngl



  !     at this stage we have gotten rid of all the stuff
  !     that are specific absorber gas dependent. The effects
  !     of all the absorbing gasses in this specific pseudo-
  !     band overlap (or no-overlap) is accumulated in tg.

  ! Initializing dble precision version of output variables 
   t8    = dble(0.)
   r8    = dble(0.)
   tc8   = dble(0.)
   sigu8 = dble(0.)
   sigd8 = dble(0.)
   re8   = dble(0.)
   vd8   = dble(0.)
   td8   = dble(0.)
   vu8   = dble(0.)
  ! Transferring data from single to double precision
   tg8   = dble(tg)
   tp8   = dble(tp)
   tcr8  = dble(tcr)
   omgp8 = dble(omgp)
   asym8 = dble(asym)
   alb8  = dble(alb)
   slr8  = dble(slr)
   amu08 = dble(amu0)

  !     get total optical depth, single scattering albedo
  !     and assymetry parameter
  amu0i = 1.0 / amu08

  re8(nz) = 0.
  vd8(nz) = 0.
  expt1 = 1.
  fd8(nz) = amu08*slr8
  tc8(nz) = 0.

  do iz=nz,2,-1
     tau = tg8(iz) + tp8(iz) + tcr8(iz)
     omg0 = min(.999999,(tcr8(iz) + omgp8(iz) * tp8(iz)) / max(1.e-20, tau))
     af = asym8(iz) * omgp8(iz) * tp8(iz) / (omg0 * max(1.e-20, tau))
  !           do delta-m scaling (wiscombe)
     fact = af * af
     tau = (1.0 - omg0 * fact) * tau
     omg0 = ( (1.0 - fact) * omg0) / (1.0 - omg0 * fact)

  !           determine the ODE matrix coefficients (Ritter and Geleyn)

     beta0 = (4.0 + af) / (8.0 * (1.0 + af))
     g1 = diffac * (1.0 - omg0 * (1.0 - beta0))
     g2 = diffac * omg0 * beta0
     gg = sqrt(g1**2 - g2**2)

  !           determine the local (true) reflection and transmission coefficients

     rinf = g2 / (gg + g1)
     ggtau = tau * gg
     if (ggtau > flush) then
        expp1=0.
     else
        expp1=dexp(-ggtau)
     end if
     expp2 = expp1**2
     denomi = 1./(1.0 - rinf**2 * expp2)
     t8(iz) = ((1.0 - rinf**2) * expp1) * denomi
     r8(iz) = rinf * (1.0 - expp2) * denomi
    
  !        get the source functions, go from top down to accomodate solar terms

     if ((gg - amu0i) < eps .and. gg > amu0i) then
        fact = 1.0 / (gg**2 - 1.0 / (amu08 + eps)**2)
     elseif ( (amu0i - gg) < eps .and. gg <= amu0i) then
        fact = 1.0 / (gg**2 - 1.0 / (amu08 - eps)**2)
     else
        fact = 1.0 / (gg**2 - amu0i**2)
     endif

     cc = omg0 * slr8 * fact
     g3 = 0.5 - 0.75 * af * amu08 / (1.0 + af)
     g4 = 1.0 - g3
     aa = g3 * (g1 - amu0i) + g4 * g2
     bb = g4 * (g1 + amu0i) + g3 * g2
     tc8(iz-1) = tc8(iz) + tau
     tcm2 = tc8(iz-1) * amu0i
     if (tcm2 > flush) then
        exp2 = 0.
     else
        exp2 = exp(-tcm2)
     end if
     exp1 = expt1
     expt1 = exp2
     sigu8(iz) = cc * ( (aa - r8(iz) * bb) * exp1 - aa * t8(iz) * exp2 )
     sigd8(iz) = cc * ( -bb * t8(iz) * exp1 + (bb - r8(iz) * aa) * exp2)
     fd8(iz-1) = amu08 * slr8 * exp2
     td8(iz)   = 1.0 - re8(iz) * r8(iz)
     re8(iz-1) = r8(iz) + t8(iz)**2 * re8(iz) / td8(iz)

     vd8(iz-1) = max(tinyreal,sigd8(iz) + ( t8(iz) * vd8(iz) + t8(iz) * re8(iz) * sigu8(iz) ) / td8(iz))
     vu8(iz-1) = ( r8(iz) * vd8(iz) + sigu8(iz) ) / td8(iz)
  enddo

  !     specify the boundary conditions
  if (tc8(1)*amu0i > flush) then
     fu8(1) = alb8 * vd8(1)  / (1.0 - alb8 * re8(1))
  else
     fu8(1) = alb8 * (vd8(1) + slr8 * amu08 * exp(-tc8(1) * amu0i)) / (1.0 - alb8 * re8(1))
  end if
  !     do adding, going from top down
  !     calculate fluxes going up through the layers

  do iz = 2, nz
     fd8(iz-1) = re8(iz-1) * fu8(iz-1) + vd8(iz-1) + fd8(iz-1)
     fu8(iz) = t8(iz) * fu8(iz-1) / td8(iz) + vu8(iz-1)
  enddo

  do iz=2,nz
     fd(iz)   = dble2sngl(tinyreal,fd8(iz)  )
     fu(iz)   = dble2sngl(tinyreal,fu8(iz)  )
     t(iz)    = dble2sngl(tinyreal,t8(iz)   )
     r(iz)    = dble2sngl(tinyreal,r8(iz)   )
     tc(iz)   = dble2sngl(tinyreal,tc8(iz)  )
     sigu(iz) = dble2sngl(tinyreal,sigu8(iz))
     sigd(iz) = dble2sngl(tinyreal,sigd8(iz))
     re(iz)   = dble2sngl(tinyreal,re8(iz)  )
     vd(iz)   = dble2sngl(tinyreal,vd8(iz)  )
     td(iz)   = dble2sngl(tinyreal,td8(iz)  )
     vu(iz)   = dble2sngl(tinyreal,vu8(iz)  )
  end do
  fd(1)   = dble2sngl(tinyreal,fd8(1)  )
  fu(1)   = dble2sngl(tinyreal,fu8(1)  )

  return
end subroutine flxsw
!
! ---------------------------------------------------------------------
!
subroutine flxlw(nz,tg,tp,tcr,omgp,src,t,r,tc,sigu,sigd,re,vd,td,vu   &
                ,fu,fd,chck,asym,mynum)
   use rconstants, only: halfpi,pi1
   implicit none
   integer            , intent(in)  :: nz,mynum
   real, dimension(nz), intent(in)  :: tg,tp,tcr,omgp,src,asym
   real               , intent(in)  :: chck
   ! Output variables
   real, dimension(nz), intent(out) :: fu,fd
   real, dimension(nz), intent(out) :: t,r,tc,sigu,sigd
   real, dimension(nz), intent(out) :: re,vd,td,vu

   !Local variables
   real(kind=8), dimension(nz)      :: tg8,tp8,tcr8,omgp8,src8,asym8,fu8,fd8
   real(kind=8), dimension(nz)      :: t8,r8,tc8,sigu8,sigd8
   real(kind=8), dimension(nz)      :: re8,vd8,td8,vu8
   real(kind=8)                     :: chck8

   real(kind=8)                     :: tau,omg0,af,fact,beta0,g1,g2,gg
   real(kind=8)                     :: rinf,ggtau,expp1,expp2,aa,bb,cc

   real(kind=8), parameter :: diffac=1.66, tinyreal=1.d-20,flush=46.
   real        , external  :: dble2sngl
   integer :: iz

   !Initialize double precision variables
   ! Output, initialized with zero
   t8    = dble(0.)
   r8    = dble(0.)
   tc8   = dble(0.)
   sigu8 = dble(0.)
   sigd8 = dble(0.)
   re8   = dble(0.)
   vd8   = dble(0.)
   td8   = dble(0.)
   vu8   = dble(0.)
   ! Input, data transformed into double precision
   tg8   = dble(tg)
   tp8   = dble(tp)
   tcr8  = dble(tcr)
   omgp8 = dble(omgp)
   src8  = dble(src)
   asym8 = dble(asym)
   chck8 = dble(chck)

   !     get total optical depth, single scattering albedo
   !     and assymetry parameter
   do iz=nz,2,-1

      tau = tg8(iz) + tp8(iz) + tcr8(iz) * chck8

      if ( tau == 0.0 ) then
         omg0 = 0.
         af = 0.
      else
         omg0 = min(.999999, omgp8(iz) * tp8(iz) / max(1.d-20, tau))
         af = asym8(iz)
      end if
      !           do delta-m scaling (wiscombe)
      fact = af * af
      tau = (1.0 - omg0 * fact) * tau
      omg0 = ((1.0 - fact) * omg0) / (1.0 - omg0 * fact)

      ! Determine the ODE matrix coefficients (Ritter and Geleyn)
      beta0 = (4.+af)/(8.*(1.+af))
      g1 = diffac*(1.-omg0*(1.-beta0))
      g2 = diffac*omg0*beta0
      gg = sqrt(g1**2-g2**2)
   
      ! Determine the local (true) reflection and transmission coefficients

      rinf  = g2/(gg+g1)
      ggtau = gg*tau
      if (ggtau > flush) then
         expp1=0.
      else
         expp1=dexp(-ggtau)
      end if
      expp2=expp1**2
      t8(iz) = (1.0-rinf**2)*expp1 / (1.0-rinf**2*expp2)
      r8(iz) = rinf*(1.-expp2) / (1.-rinf**2*expp2)
      ! get the source functions, go from top down to accomodate solar terms

      if (tau < 4.d-2) then      !changed June 12 after Jerry's recom.
         sigu8(iz) = dble(halfpi) * (src8(iz) + src8(iz-1)) * tau * diffac
         sigd8(iz) = sigu8(iz)
      else
         aa =  (g1 + g2) * (1.0 - r8(iz)) - (1.0 + r8(iz) - t8(iz)) / tau
         bb = -(g1 + g2) * t8(iz) + (1.0 + r8(iz) - t8(iz)) / tau
         cc = diffac * pi1 * (1.0 - omg0) / gg**2
         sigu8(iz) = cc*(aa*src8(iz)+bb*src8(iz-1))
         sigd8(iz) = cc*(bb*src8(iz)+aa*src8(iz-1))
         if (sigu8(iz) < 0.0 .or. sigd8(iz) < 0.0) then
            print*,'negative source in flxlw: iz',iz
            print*,'aa,bb,cc:',aa,bb,cc
            print*,'src(iz), src(iz-1):',src8(iz),src8(iz-1)
            print*,'g1, g2, tau:',g1,g2,tau
            print*,'r(iz), t(iz):',r8(iz),t8(iz)
            print*,'sigu(iz), sigd(iz):',sigu8(iz),sigd8(iz)
            call abort_run('Negative source','flxlw','harr_rad.f90')
         end if
      end if
   
   end do

   !     do adding
   !     initialize

   re8(nz)  = 0.
   vd8(nz)  = 0.
   fd8(nz) = 0.
   fu8(1)  = dble(pi1)*src8(1)

   do iz=nz,2,-1
      td8(iz)   = 1. - re8(iz)*r8(iz)
      re8(iz-1) = max(tinyreal,r8(iz) + t8(iz)**2*re8(iz) / td8(iz))
      vd8(iz-1) = sigd8(iz) + ( t8(iz)*vd8(iz) + t8(iz)*re8(iz)*sigu8(iz) ) / td8(iz)
      vu8(iz-1) = ( r8(iz)*vd8(iz) + sigu8(iz) ) / td8(iz)
!     write (unit=*,fmt='(a,1x,i3,1x,8(a,1x,es10.3,1x))') &
!        '     ┐┐┐┐ iz=',iz,'tau=',tau,'tg=',tg(iz),'tp=',tp(iz),'tcr=',tcr(iz),'omgp=',omgp(iz),'asym=',asym(iz),'omg0=',omg0,'fact=',fact
   enddo

   !     calculate fluxes going up through the layers
   do iz=2,nz
      fd8(iz-1) = re8(iz-1)*fu8(iz-1) + vd8(iz-1)
      fu8(iz)   = t8(iz)*fu8(iz-1)/td8(iz)+vu8(iz-1)
   enddo
   
   !Transfer values back to single precision
   do iz=2,nz
      fd(iz)   = dble2sngl(tinyreal,fd8(iz)  )
      fu(iz)   = dble2sngl(tinyreal,fu8(iz)  )
      t(iz)    = dble2sngl(tinyreal,t8(iz)   )
      r(iz)    = dble2sngl(tinyreal,r8(iz)   )
      tc(iz)   = dble2sngl(tinyreal,tc8(iz)  )
      sigu(iz) = dble2sngl(tinyreal,sigu8(iz))
      sigd(iz) = dble2sngl(tinyreal,sigd8(iz))
      re(iz)   = dble2sngl(tinyreal,re8(iz)  )
      vd(iz)   = dble2sngl(tinyreal,vd8(iz)  )
      td(iz)   = dble2sngl(tinyreal,td8(iz)  )
      vu(iz)   = dble2sngl(tinyreal,vu8(iz)  )
   end do
   fd(1)   = dble2sngl(tinyreal,fd8(1)  )
   fu(1)   = dble2sngl(tinyreal,fu8(1)  )


   return
end subroutine flxlw
!
! ---------------------------------------------------------------------

subroutine rayleigh(wlnlo,wlnhi,rayavg)
  use mem_harr, only: sun
  implicit none

  real, intent(in)  :: wlnlo,wlnhi
  real, intent(out) :: rayavg
  real :: an0,ssum,h1,sum,wl1,f1,fac,f2,h2,wl2,al
  integer :: num,i

  !     this stuff comes from Houghton 1985 (book), see also Slingo
  !     and Schrecker 1982 (QJ)
  !     calculate constant part (flux weighted) of rayleigh scattering
  !     rayleigh scattering = rayavg*delta z*press/temp
  !     written JV May 1993

  real, parameter :: an01=1.000064328,an02=0.0294981,an03=.2554e-3
  real, parameter :: an0d1=146.,an0d2=41.

  an0(al) = an01+an02/(an0d1-al**2)+an03/(an0d2-al**2)

  num = int(wlnhi-wlnlo)+1
  num = min(max(25,num),5000)
  sum = 0.
  ssum= 0.
  wl1 = wlnlo
  f1  = (an0(wl1)**2-1.)**2/wl1**4*sun(1.e4/wl1)
  h1  = sun(1.e4/wl1)
  do i=2,num
     fac = (real(i)-1)/(real(num)-1.)
     wl2 = wlnlo+(wlnhi-wlnlo)*fac
     f2  = (an0(wl2)**2-1.)**2/wl2**4*sun(1.e4/wl2)
     h2  = sun(1.e4/wl2)
     sum = sum + 0.5*(f1+f2)*(wl2-wl1)
     ssum= ssum + 0.5*(h1+h2)*(wl2-wl1)
     wl1 = wl2
     f1  = f2
     h1  = h2
  enddo
  rayavg = 0.945319e-2*sum/ssum

  return
end subroutine rayleigh
!
! ---------------------------------------------------------------------
!
subroutine csband(wlnlo,wlnhi,csavg)
  implicit none
  real :: wlnlo,wlnhi,csavg
  real :: cs,wnhi,wnlo,sum,wn1,f1,fac,f2,wn,wn2,plkavg
  integer :: num,i
  !
  !     calculate the blackbody flux (t=296.) weighted self-
  !     broadening coefficient for water vapor. See Kniezys
  !     et al 1980 (LOWTRAN 5). Units are converted to have
  !     the water vapor content input a Pascal.

  cs(wn)  = 4.18+5578.*exp(-7.87e-3*wn)

  wnhi = 1.e4/wlnlo
  wnlo = 1.e4/wlnhi
  num  = int(wnhi-wnlo)+1
  num  = min(max(25,num),5000)
  sum  = 0.
  wn1  = wnlo
  f1   = cs(wn1)
  do i = 2,num
     fac = (real(i)-1)/(real(num)-1.)
     wn2 = wnlo+(wnhi-wnlo)*fac
     f2  = cs(wn2)
     sum = sum + 0.5*(f1+f2)*plkavg(wn1,wn2,296.)
     wn1 = wn2
     f1  = f2
  enddo
  csavg = sum/(plkavg(wnlo,wnhi,296.)*1013250.*9.81)
  return
end subroutine csband
!
! ---------------------------------------------------------------------
!
real function  plkavg( wnumlo, wnumhi, t )

  implicit none

  !
  ! Computes Planck function integrated between two wavenumbers
  !
  ! INPUT:
  !
  !   wnumlo : lower wavenumber ( inv cm ) of spectral interval
  !   wnumhi : upper wavenumber
  !   t       : temperature (k)
  !
  ! OUTPUT        
  !
  !   plkavg : integrated Planck function ( watts/sq m )
  !
  ! References-- (1) Houghton,physics of atmospheres,appendix 7
  !              (2) Specifications of the physical world: new value
  !                  of the fundamental constants, dimensions/n.b.s.
  !                  Jan. 1974
  !
  ! Method-- Houghton's exponential series is used for v.gt.vcut
  !          ( 'v' is Houghton's notation ) and his power series
  !          in v for v.le.vcut.  more terms are taken in the
  !          exponential series, the larger v is.  ( note that
  !          Houghton's assessment that the power series is useful
  !          for  v.lt.2*pi  is incorrect--vcut must be less than
  !          2 just in order to get 4-5 significant digits. )
  !
  ! Accuracy-- 6 significant digits
  !
  ! ARGUMENTS
  real ::    t, wnumlo, wnumhi
  ! LOCAL CONSTANTS

  real, parameter :: c2     = 1.438786     ! second radiation constant
  real, parameter :: cona   = .11573303e-7 ! 15 / pi**4 / 13305600
  real, parameter :: conb   = 48.888889    ! 440 / 9
  real, parameter :: f15pi4 = 0.15398973   ! 15 / pi**4
  real, parameter :: sigdpi = 1.804919e-8  ! stefan-boltzmann constant divided by pi
  real, parameter :: vcut   = 1.5          ! power-series cutoff point
  ! vcp      :  exponential series cutoff points
  real, dimension(7), parameter :: vcp=(/10.25, 5.7, 3.9, 2.9, 2.3, 1.9, 0.0 /)

  ! LOCAL VARIABLES
  real               :: ex ! exp( - v )
  real               :: mv ! multiples of *v*
  real, dimension(2) :: d  ! integral of normalized planck function
                           !   from 0 to current wavenumber
  real               :: v  ! h*c*nu / (k*t), where h=plancks constant,
                           !   c=speed of light, k=boltzmann constant,
                           !   nu=wavenumber
  integer :: mmax          ! no. of terms to take in exponential series
  real :: vsq
  integer :: i,m
  !
  !
  if( t <= 0.0 .or. wnumhi <= wnumlo .or. wnumlo <= 0. ) stop '╗╗╗ ahah: plkavg (harr_rad.f90)'

  do i = 1, 2
     if(i == 1) then
        v = ( c2 / t ) * wnumlo
     elseif(i == 2)  then 
        v = ( c2 / t ) * wnumhi
     end if

     if( v < vcut )  then
  !
  !            *** Houghton's power series (factored using Horner's rule)
        vsq = v**2
        d(i) = 1.0 - cona *vsq * v * ( 4435200. + v * ( -1663200. &
                 + v * ( 221760. + vsq * ( -2640. + vsq * ( conb  &
                 - vsq ) ) ) ) )
     else
  !
  !*** Houghton's exponential series
  !*** set upper limit of series depending on value of v
        mmax=0
        uplimit: do
          mmax=mmax+1
          if (v >= vcp(mmax)) exit uplimit
        end do uplimit

        ex = exp(-v)
        d(i)= ex * (6. + v * ( 6. + v * ( 3. + v ) ) )

        do m = 2, mmax
           mv = m * v
           d(i) = d(i) + ex**m * ( 6. + mv * ( 6. + mv *  &
                                 ( 3. + mv ) ) ) / m**4
        end do
        d(i) = f15pi4 * d(i)
     end if
  end do
  plkavg = sigdpi * t**4 * (d(1) - d(2))
  return
end function plkavg
!
! ---------------------------------------------------------------------
!
real function sunavg(wnumlo, wnumhi, solcon)
  use mem_harr, only: sun
  implicit none

  real :: wnumlo, wnumhi, solcon
  real :: scln=1372.844

  integer :: num,i
  real :: v1,s1,v2,s2,fac2,sum
  ! Spectral solar flux  (after LOWTRAN7)
  !
  ! On input:
  ! wnumlo --- lower wavelength of a spectral interval (wavenumbers)
  ! wnumhi --- upper wavelength of a spectral interval (wavenumbers)
  ! solcon - value of solar "constant" in Watts/ sq m
  !
  ! On output:
  ! Total extraterrestrial solar flux (W/sq m) between
  ! wavelengths wnumlo, wnumhi
  ! Trapezoidal integration is used with the resolution 1 cm^-1
  ! or such that there is at least 25 points between wnumhi and
  ! wnumlo but no more than 5000 points.
  !
  num = int(wnumhi-wnumlo)+1
  num = min(max(25, num),5000)
  sum = 0.
  v1 = wnumlo
  s1 = sun(v1)
  do i=2, num
    fac2 = (real(i)-1.)/(real(num)-1.)
    v2 =  wnumlo+(wnumhi-wnumlo)*fac2
    s2 = sun(v2)
    sum = sum + (s1+s2)/2.*(10000./v1-10000./v2)
    v1 = v2
    s1 = s2
  enddo
  sunavg = (solcon/scln)*sum

  return
end function sunavg
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
real function dble2sngl(lowerthreshold,dblevar)
   implicit none
   real(kind=8), intent(in) :: lowerthreshold,dblevar
   
   if (dblevar == 0.) then
      dble2sngl = 0.
   else
      ! This will offset the value by a minimum threshold, but it will preserv the sign.
      dble2sngl = sngl(sign(max(lowerthreshold,abs(dblevar)),dblevar))
   end if
   
   return
end function dble2sngl
!==========================================================================================!
!==========================================================================================!
