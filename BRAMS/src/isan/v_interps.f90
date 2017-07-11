!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

SUBROUTINE obs_isen (m1,m2,tsnd,psnd,zsnd,rsnd,usndz,vsndz  &
                    ,zsndz,stlt,stln,topsnd,lpst,lzst  &
                    ,obsu,obsv,obsp,obss,obsr,no,chstid)
use isan_coms
use rconstants

implicit none
integer :: m1,m2,no
real :: tsnd(m1,m2),psnd(m1,m2),zsnd(m1,m2)  &
         ,rsnd(m1,m2),stlt(m1),stln(m1),topsnd(m1),usndz(m1,m2)  &
         ,vsndz(m1,m2),zsndz(m1,m2)  &
         ,obsu(no,*),obsv(no,*),obsp(no,*),obss(no,*),obsr(no,*)
integer, dimension(m1) :: lpst,lzst
character(len=8), dimension(*) :: chstid

real,dimension(maxlev) :: pk,pth,zi
integer ns,lp,lz,k,lmbot,nglev,lmtop,l,lbchyd,lbc,lbcp
real :: th,wt,pkn,syo,po,tho
logical :: failure

staloop: do ns=1,nsta
   lp=lpst(ns)
   lz=lzst(ns)

   do k=1,nisn
      obsp(ns,k)=1E30
      obsr(ns,k)=1E30
      obsu(ns,k)=1E30
      obsv(ns,k)=1E30
      obss(ns,k)=1E30
   enddo

   lmbot=0
   nglev=0
   do k=1,lp
      pk(k)=1.e30
      if(psnd(ns,k) < 1.e29) pk(k)=psnd(ns,k)**rocp
      pth(k)=1.e30
      if(tsnd(ns,k) < 1.e29 .and. psnd(ns,k) < 1.e29)  &
           pth(k)=tsnd(ns,k)*(p00/psnd(ns,k))**rocp
      if(pth(k) < 1.e30 .and. lmbot == 0) lmbot=k
      if(pth(k) < 1.e29) nglev=nglev+1
   enddo
   lmtop=0
   do k=lp,1,-1
      if(pth(k) < 1.e30 .and. lmtop == 0) lmtop=k
   enddo
   if(nglev >= 1) cycle staloop

   l=2
   isnloop1: do k=1,nisn
      th=levth(k)
      if(th < pth(lmbot) .or. th > pth(lmtop)) then
         obsp(ns,k)=1e30
         obsr(ns,k)=1e30
         cycle isnloop1
      endif

      do while (.not. (th >= pth(l-1) .and. th.lt.pth(l)) )
        l=l+1
        if(l.gt.lp) exit isnloop1
      end do
      wt=(th-pth(l-1))/(pth(l)-pth(l-1))
      pkn=(pk(l)-pk(l-1))*wt+pk(l-1)
      obsp(ns,k)=pkn**cpor
      if(rsnd(ns,l) < 1e19 .and. rsnd(ns,l-1) < 1e19) then
         obsr(ns,k)=rsnd(ns,l-1)+(rsnd(ns,l)-rsnd(ns,l-1))*wt
      else
         obsr(ns,k)=1e30
      endif

   end do isnloop1

   ! Find a boundary condition as first level at or below 360K

   lbchyd=360
   failure=.true.
   isnloop2: do k=nisn,1,-1
      !print*,'mjb p',k,ns,levth(k),lbchyd,obsp(ns,k)
      if(levth(k).le.lbchyd.and.obsp(ns,k).lt.1e19) then
         lbc=k
         failure=.false.
         exit isnloop2
      endif
   end do isnloop2
   !stop 'obs_isen'
   !MJW Seems like just set everything to missing and continue rather than stop
   if (failure) then
     print*,'---Could not find 1st obs isentropic boundary level ',ns,' ',chstid(ns)
     do k=1,nisn
        obss(ns,k)=1.e30
        obsu(ns,k)=1.e30
        obsv(ns,k)=1.e30
     enddo
     cycle staloop
   end if

   failure=.true.
   lploop: do k=lp,1,-1
      if(pth(k).le.float(levth(lbc))) then
         lbcp=k
         failure=.false.
         exit lploop
      end if
   end do lploop
   if (failure) then
     print*,'---Could not find 2nd obs isentropic boundary level ',ns,' ',chstid(ns)
   !stop 'obs_isen'
   !MJW Seems like just set everything to missing and continue rather than stop
     do k=1,nisn
        obss(ns,k)=1.e30
        obsu(ns,k)=1.e30
        obsv(ns,k)=1.e30
     end do
     cycle staloop
   end if

   syo=cpdry*pth(lbcp)*pk(lbcp)/p00**rocp+grav*zsnd(ns,lbcp)
   obss(ns,lbc)=syo+cpdry*(pk(lbcp)+obsp(ns,lbc)**rocp)  &
        *.5/p00**rocp *(levth(lbc)-pth(lbcp))
   po=obsp(ns,lbc)
   syo=obss(ns,lbc)
   tho=levth(lbc)
   do k=lbc+1,nisn
      obss(ns,k)=1e30
      if(obsp(ns,k) < 1e19) then
         obss(ns,k)=syo+cpdry*(po**rocp+obsp(ns,k)**rocp)  &
              *.5/p00**rocp *(levth(k)-tho)
         syo=obss(ns,k)
         po=obsp(ns,k)
         tho=levth(k)
      endif
   enddo

   syo=obss(ns,lbc)
   po=obsp(ns,lbc)
   tho=levth(lbc)
   do k=lbc-1,1,-1
      obss(ns,k)=1e30
      if(obsp(ns,k).lt.1e19) then
         obss(ns,k)=syo+cpdry*(po**rocp+obsp(ns,k)**rocp)*.5  &
              /p00**rocp*(levth(k)-tho)
         syo=obss(ns,k)
         po=obsp(ns,k)
         tho=levth(k)
      endif
   enddo

   do k=1,nisn
      if(obss(ns,k).lt.1e19) then
         zi(k)=(obss(ns,k)-cpdry*levth(k)*(obsp(ns,k)*p00i)**rocp)/grav
      else
         zi(k)=1e30
      endif
   enddo

![MLO - Changing if to avoid crashing with Check bounds (and removing some goto statements...)
   l=2
   isnloop3: do k=1,nisn

      if(lz==0) then
         obsu(ns,k)=1e30
         obsv(ns,k)=1e30
         cycle isnloop3
      elseif(zi(k) < zsndz(ns,1) .or. zi(k) > zsndz(ns,lz)) then
         obsu(ns,k)=1e30
         obsv(ns,k)=1e30
         cycle isnloop3
      endif

      do while (.not. (zi(k) >= zsndz(ns,l-1) .and. zi(k) < zsndz(ns,l)) )
        l=l+1
        if(l.gt.lz) cycle staloop
      end do
      wt=(zi(k)-zsndz(ns,l-1))/(zsndz(ns,l)-zsndz(ns,l-1))
      if(usndz(ns,l).lt.1e19.and.usndz(ns,l-1).lt.1e19  &
           .and.vsndz(ns,l).lt.1e19.and.vsndz(ns,l-1).lt.1e19) then
         obsu(ns,k)=usndz(ns,l-1)+(usndz(ns,l)-usndz(ns,l-1))*wt
         obsv(ns,k)=vsndz(ns,l-1)+(vsndz(ns,l)-vsndz(ns,l-1))*wt
      else
         obsu(ns,k)=1e30
         obsv(ns,k)=1e30
      endif
      
   end do isnloop3

end do staloop

return
end

!***************************************************************************

SUBROUTINE obs_sigz (m1,m2,tsnd,psnd,zsnd,rsnd,usndz,vsndz  &
                    ,zsndz,stlt,stln,topsnd,lpst,lzst,sndtopg  &
                    ,obsu,obsv,obsp,obst,obsr,no,ztop)
use isan_coms
use rconstants
use therm_lib,only: ptrh2rvapl,virtt

implicit none
integer :: m1,m2,no,lpst(m1),lzst(m1)
real ::   tsnd(m1,m2),psnd(m1,m2),zsnd(m1,m2)  &
         ,rsnd(m1,m2),stlt(m1),stln(m1),topsnd(m1),usndz(m1,m2)  &
         ,vsndz(m1,m2),zsndz(m1,m2),sndtopg(m1)  &
         ,obsu(no,*),obsv(no,*),obsp(no,*),obst(no,*),obsr(no,*)
         
real :: pk(maxlev),pth(maxlev),zi(maxlev),sigzr(maxsigz)  &
         ,tsz(maxsigz),rsz(maxsigz),psz(maxsigz),sdat(2)

integer :: ns,k,lmbot,lp,lz,nglev,lmtop,l,lbc,lbcp,kbcu
real :: ztop,wt,bchyd,pio,zso,tho,rs
logical :: failure

staloop: do ns=1,nsta
   lp=lpst(ns)
   lz=lzst(ns)
   
   !DO K=1,nsigz
   !   print*,'mjb a5',ns,nsigz,k,zsnd(ns,k)
   !enddo

   ! Fill sigzr from interpolated grid topography height

   do K=1,nsigz
      sigzr(K)=sndtopg(ns)+sigz(k)*(1.-sndtopg(ns)/ztop)
   end do

   do k=1,nsigz
      obsp(ns,k)=1E30
      obsr(ns,k)=1E30
      obst(ns,k)=1E30
   enddo

   lmbot=0
   nglev=0
   do k=1,lp
      pk(k)=1.e30
      if(psnd(ns,k) < 1.e29) pk(k)=psnd(ns,k)**rocp
      pth(k)=1.e30
      if(tsnd(ns,k) < 1.e29 .and. psnd(ns,k) < 1.e29)  &
         pth(k)=tsnd(ns,k)*(p00/psnd(ns,k))**rocp
      if(pth(k) < 1.e30 .and. lmbot == 0) lmbot=k
      if(pth(k) < 1.e29) nglev=nglev+1
   enddo
   lmtop=0
   do k=lp,1,-1
      if(pth(k) < 1.e30 .and. lmtop == 0) lmtop=k
   enddo

   if(nglev > 1) then

     l=2

     !print*,'mjb a',ns,zsnd(ns,1)
     sigzloop1: do k=1,nsigz
        if(sigzr(k) < zsnd(ns,1)) then
           obst(ns,k)=pth(1)
           obsr(ns,k)=rsnd(ns,1)
           cycle sigzloop1
        endif
        if(sigzr(k) > zsnd(ns,lp)) then
           obst(ns,k)=1e30
           obsr(ns,k)=1e30
           cycle sigzloop1
        endif

        do while (.not. (sigzr(k) >= zsnd(ns,l-1) .and. sigzr(k) < zsnd(ns,l)) )
          l=l+1
          if(l.gt.lp) exit sigzloop1
        end do
        wt=(sigzr(k)-zsnd(ns,l-1))/(zsnd(ns,l)-zsnd(ns,l-1))
        obst(ns,k)=(pth(l)-pth(l-1))*wt+pth(l-1)

        if(rsnd(ns,l) < 1e19 .and. rsnd(ns,l-1) < 1e19) then
           obsr(ns,k)=rsnd(ns,l-1)+(rsnd(ns,l)-rsnd(ns,l-1))*wt
        else
           obsr(ns,k)=1e30
        endif
        !print*,'mjb b',k,ns,lp,obst(ns,k),obsr(ns,k)
     end do sigzloop1

     ! Find a boundary condition as first level at or below 10000m

     bchyd=10000.
     failure=.true.
     sigzloop2: do k=nsigz,1,-1
        !print*,'mjb c',ns,k,sigzr(k),bchyd,obst(ns,k)
        if(sigzr(k) <= bchyd .and. obst(ns,k) < 1e19) then
           lbc=k
           failure=.false.
           exit sigzloop2
        end if
     end do sigzloop2
     if (failure) then 
       print*,'---Could not find 1st obs sigma-z boundary level',ns
     !stop 'obs_sigz'
     !MJW Seems like just set everything to missing and continue rather than stop
       do k=1,nsigz
          obsp(ns,k)=1.e30
          obst(ns,k)=1.e30
          obsu(ns,k)=1.e30
          obsv(ns,k)=1.e30
       enddo
       cycle staloop
     end if

     failure=.true.
     lploop: do k=lp,1,-1
        if(zsnd(ns,k) <= sigzr(lbc)) then
           lbcp=k
           failure=.false.
           exit lploop
        endif
     end do lploop
     if (failure) then
       print*,'---Could not find 2nd obs sigma-z boundary level',ns
     !stop 'obs_sigz'
     !MJW Seems like just set everything to missing and continue rather than stop
       do k=1,nsigz
          obsp(ns,k)=1.e30
          obst(ns,k)=1.e30
          obsu(ns,k)=1.e30
          obsv(ns,k)=1.e30
       enddo
       cycle staloop
     end if

     pio=cpdry*(psnd(ns,lbcp)/p00)**rocp
     obsp(ns,lbc)=pio-(sigzr(lbc)-zsnd(ns,lbcp))*grav/((pth(lbcp)+obst(ns,lbc))*.5)
     pio=obsp(ns,lbc)
     zso=sigzr(lbc)
     tho=obst(ns,lbc)
     do k=lbc+1,nsigz
        obsp(ns,k)=1e30
        if(obst(ns,k) < 1e19) then
           obsp(ns,k)=pio-(sigzr(k)-zso)*grav/((obst(ns,k)+tho)*.5)
           zso=sigzr(k)
           pio=obsp(ns,k)
           tho=obst(ns,k)
        endif
     enddo

     zso=sigzr(lbc)
     pio=obsp(ns,lbc)
     tho=obst(ns,lbc)
     do k=lbc-1,1,-1
        obsp(ns,k)=1e30
        if(obst(ns,k) < 1e19) then
           obsp(ns,k)=pio-(sigzr(k)-zso)*grav/((obst(ns,k)+tho)*.5)
           zso=sigzr(k)
           pio=obsp(ns,k)
           tho=obst(ns,k)
        endif
     enddo

     ! Compute virtual temp

     DO K=1,nsigz
        IF(obsp(ns,K)+obst(ns,k) < 1E19) THEN
           tsz(k)=obst(ns,k)*obsp(ns,k)/cpdry
           psz(K)=(obsp(ns,k)/cpdry)**cpor*p00
           IF(obsr(ns,k).LT.1E19) THEN
              rsz(k)=ptrh2rvapl(obsr(ns,k),psz(k),tsz(k),.false.)
              tsz(k)=virtt(obst(ns,k),rsz(k))
           else
              tsz(k)=obst(ns,k)
           endif
        endif
     enddo

     ! Recompute pressure profile with virtual temp
     failure=.true.
     sigzloop3: do k=nsigz,1,-1
        if(obsp(ns,k).lt.1e19) then
           kbcu=k
           failure=.false.
           exit sigzloop3
        endif
     end do sigzloop3
     if (failure) then
       print*, '---Could not find good obs sigma-z boundary level 3',ns
       !stop 'obs_sigz'
       !MJW Seems like just set everything to missing and continue rather than stop
       do k=1,nsigz
          obsp(ns,k)=1.e30
          obst(ns,k)=1.e30
          obsu(ns,k)=1.e30
          obsv(ns,k)=1.e30
       enddo
       cycle staloop
     end if

     zso=sigzr(kbcu)
     pio=obsp(ns,kbcu)
     tho=tsz(kbcu)
     do k=kbcu-1,1,-1
        obsp(ns,k)=1e30
        if(obst(ns,k).lt.1e19) then
           obsp(ns,k)=pio-(sigzr(k)-zso)*grav/((tsz(k)+tho)*.5)
           zso=sigzr(k)
           pio=obsp(ns,k)
           tho=tsz(k)
        endif
     enddo
     do k=1,nsigz
        if(obsp(ns,k) < 1e19) obsp(ns,k)=(obsp(ns,k)/cpdry)**cpor*p00
     enddo

   ! Vertically interpolate winds in height
   end if

![MLO - Some changes to remove some goto and to avoid crashing when using check bounds
   l=2
   sigzloop4: do k=1,nsigz

![MLO - Changed here so if lz == 0 zsndz doesn't 
      if(lz==0) then
         obsu(ns,k)=1e30
         obsv(ns,k)=1e30
         cycle sigzloop4
      elseif(sigzr(k) < zsndz(ns,1) .or. sigzr(k) > zsndz(ns,lz)) then
         obsu(ns,k)=1e30
         obsv(ns,k)=1e30
         cycle sigzloop4
      endif

      do while (.not. (sigzr(k) >= zsndz(ns,l-1) .and. sigzr(k) < zsndz(ns,l)))
        l=l+1
        if (l > lz) exit sigzloop4
      end do
      wt=(sigzr(k)-zsndz(ns,l-1))/(zsndz(ns,l)-zsndz(ns,l-1))
      if(usndz(ns,l).lt.1e19.and.usndz(ns,l-1).lt.1e19  &
           .and.vsndz(ns,l).lt.1e19.and.vsndz(ns,l-1).lt.1e19) then
         obsu(ns,k)=usndz(ns,l-1)+(usndz(ns,l)-usndz(ns,l-1))*wt
         obsv(ns,k)=vsndz(ns,l-1)+(vsndz(ns,l)-vsndz(ns,l-1))*wt
      else
         obsu(ns,k)=1e30
         obsv(ns,k)=1e30
      endif
      
   end do sigzloop4

end do staloop
!MLO]

return
end
!***************************************************************************

subroutine vterpp_i (np1,np2,np3,npi3,un,vn,tn,zn,rn,ui2,vi2,pi2,si2,ri2)
     
use isan_coms
use rconstants

implicit none
integer :: np1,np2,np3,npi3
real :: un(np1,np2,np3),vn(np1,np2,np3),tn(np1,np2,np3)  &
         ,zn(np1,np2,np3),rn(np1,np2,np3)  &
         ,ui2(np1,np2,npi3),vi2(np1,np2,npi3),pi2(np1,np2,npi3)  &
         ,si2(np1,np2,npi3),ri2(np1,np2,npi3)

integer, parameter :: npr=maxpr+2
real :: ppd(npr),thd(npr),pkd(npr),ud(npr),vd(npr),zd(npr),rd(npr)

integer :: i,j,k,mcnt,npd,lp,kpbc,kibc
real :: thl,wt,pkn,pbc,sy

do j=1,np2
   do i=1,np1

      DO K=1,NISN
         UI2(I,J,K)=1.E30
         VI2(I,J,K)=1.E30
         RI2(I,J,K)=1.E30
         PI2(I,J,K)=1.E30
         SI2(I,J,K)=1.E30
      ENDDO

      ! determine if this column is all missing. if so, just leave
      !   isentropic data as missing

      mcnt=0
      DO K=1,NPRZ
         if(tn(i,j,k).gt.1000.) mcnt=mcnt+1
      ENDDO
      if(mcnt.eq.nprz) goto 4500

      DO K=1,NPRZ
         pnpr(k)=levpr(k)*100.
         UD(K+2)=UN(I,J,K)
         VD(K+2)=VN(I,J,K)
         RD(K+2)=RN(I,J,K)
         PPD(K+2)=PNPR(K)
         THD(K+2)=TN(I,J,K)*(P00/PNPR(K))**ROCP
      ENDDO

      ! Define two phony levels for isentropes underground

      UD(1)=UD(3)
      UD(2)=UD(3)
      VD(1)=VD(3)
      VD(2)=VD(3)
      RD(1)=RD(3)
      RD(2)=RD(3)
      PPD(1)=120000.
      PPD(2)=110000.
      THD(2)=(TN(I,J,1)-2.25)*(P00/PNPR(1))**ROCP
      THD(1)=180.
      npd=nprz+2
      DO K=1,NPD
         PKD(K)=PPD(K)**ROCP
      ENDDO

      lp = npd

      do k = nisn,1,-1
         35 continue
         thl = levth(k)

         if (thl .gt. thd(npd)) then

            pi2(i,j,k) = (thd(npd) / thl) ** cpor * ppd(npd)
            ui2(i,j,k) = ud(npd)
            vi2(i,j,k) = vd(npd)
            ri2(i,j,k) = 0.

         elseif (thl .le. thd(lp) .and. thl .gt. thd(lp-1)) then

            wt = (thl - thd(lp-1)) / (thd(lp) - thd(lp-1))
            pkn = pkd(lp-1) + (pkd(lp) - pkd(lp-1)) * wt
            pi2(i,j,k) = pkn ** cpor

            ui2(i,j,k) = ud(lp-1) + (ud(lp) - ud(lp-1)) * wt
            vi2(i,j,k) = vd(lp-1) + (vd(lp) - vd(lp-1)) * wt
            ri2(i,j,k) = rd(lp-1) + (rd(lp) - rd(lp-1)) * wt

         else

            lp = lp - 1
            if (lp .le. 1) then
               print*, 'vterpp_i interpolation tried to go below'
               print*, 'lowest pressure level.'
               stop 'vterpp_i'
            endif
            goto 35

         endif
      enddo

      ! Find a pressure b.c. as second isentrope above 400 mb

      pbc=40000.
      do k=1,npd
         if(ppd(k).lt.pbc) then
            kpbc=k
            exit
         endif
      enddo
     
      do k=2,nisn
         if(pi2(i,j,k).lt.ppd(kpbc)) then
            kibc=k
            exit
         endif
      enddo

      sy=cpdry*thd(kpbc)*pkd(kpbc)/p00k+grav*zn(i,j,kpbc-2)

      si2(i,j,kibc-1)=sy-cpdry*(pi2(i,j,kibc-1)**rocp  &
           +ppd(kpbc)**rocp)/(2.*p00k)*(thd(kpbc)-levth(kibc-1))
      do k=kibc-2,1,-1
         si2(i,j,k)=si2(i,j,k+1)+cpdry*(pi2(i,j,k+1)**rocp  &
              +pi2(i,j,k)**rocp)/(2.*p00k)  &
              *(levth(k)-levth(k+1))
      enddo

      si2(i,j,kibc)=sy+cpdry*(pi2(i,j,kibc)**rocp  &
           +ppd(kpbc)**rocp)/(2.*p00k)*(levth(kibc)-thd(kpbc))
      do k=kibc+1,nisn
         si2(i,j,k)=si2(i,j,k-1)+cpdry*(pi2(i,j,k-1)**rocp  &
              +pi2(i,j,k)**rocp)/(2.*p00k)  &
              *(levth(k)-levth(k-1))
      enddo

      4500 continue
   enddo
enddo

! If any level is above 100 mb, compute geostrophic winds

!GDKM=2.*SPCON
!DO J=2,NP2-1
!   GDLAT=(XSWLAT+(J-1)*GDATDY)*PI180
!   FCORI=1./(2.*7.292E-5*SIN(GDLAT))
!   DO I=2,NP1-1
!      DO K=NISN,1,-1
!         IF(PI2(I,J,K).LT.10000.) THEN
!            UI2(I,J,K)=-FCORI*(SI2(I,J+1,K)-SI2(I,J-1,K))/(GDKM*GDATDY)
!            VI2(I,J,K)= FCORI*(SI2(I+1,J,K)-SI2(I-1,J,K))  &
!                 /(GDKM*GDATDX*COS(GDLAT))
!         ENDIF
!      ENDDO
!   ENDDO
!ENDDO

return
end

!***************************************************************************

subroutine vterpp_s (np1,np2,np3,npi3,un,vn,tn,zn,rn  &
                    ,ui2,vi2,pi2,ti2,ri2,topt,rtgt)
     
use isan_coms
use rconstants
use therm_lib, only: ptrh2rvapil,virtt

implicit none
integer :: np1,np2,np3,npi3
real :: un(np1,np2,np3),vn(np1,np2,np3),tn(np1,np2,np3)  &
         ,zn(np1,np2,np3),rn(np1,np2,np3)  &
         ,ui2(np1,np2,npi3),vi2(np1,np2,npi3),pi2(np1,np2,npi3)  &
         ,ti2(np1,np2,npi3),ri2(np1,np2,npi3)  &
         ,topt(np1,np2),rtgt(np1,np2)

integer, parameter :: npr=maxpr+2
real :: ppd(npr),thetd(npr),pkd(npr),ud(npr),vd(npr),zd(npr)  &
         ,rd(npr),pid(npr),tempd(npr),rtd(npr),thvd(npr)  &
         ,sigzr(maxsigz),vvv(maxsigz)

integer :: i,j,k,mcnt,npd,kl,kpbc,kibc
real :: pbc,thvp,piibc,raux

do j=1,np2
   do i=1,np1

      DO K=1,npi3
         UI2(I,J,K)=1.E30
         VI2(I,J,K)=1.E30
         RI2(I,J,K)=1.E30
         PI2(I,J,K)=1.E30
         TI2(I,J,K)=1.E30
      ENDDO

      ! determine if this column is all missing.
      ! if so, just leave data as missing

      mcnt=0
      DO K=1,NPRZ
         if(tn(i,j,k).gt.1000.) mcnt=mcnt+1
      ENDDO
      if(mcnt.eq.nprz) goto 4500

      DO K=1,NPRZ
         PNPR(K)=levpr(k)*100.
         UD(K+2)=UN(I,J,K)
         VD(K+2)=VN(I,J,K)
         RD(K+2)=RN(I,J,K)
         PPD(K+2)=PNPR(K)
         THETD(K+2)=TN(I,J,K)*(P00/PNPR(K))**ROCP
         ZD(K+2)=ZN(I,J,K)
      ENDDO

      ! Define two phony levels for isentropes underground

      UD(1)=UD(3)
      UD(2)=UD(3)
      VD(1)=VD(3)
      VD(2)=VD(3)
      RD(1)=RD(3)
      RD(2)=RD(3)
      PPD(1)=120000.
      PPD(2)=110000.
      THETD(2)=(TN(I,J,1)-2.25)*(P00/PNPR(1))**ROCP
      THETD(1)=220.
      npd=nprz+2
      DO K=1,NPD
         PKD(K)=PPD(K)**ROCP
         pid(k)=cpdry*(ppd(k)/p00)**rocp
         tempd(k)=thetd(k)*pid(k)/cpdry
         rtd(k)=ptrh2rvapil(rd(k),ppd(k),tempd(k),.false.)
         thvd(k)=virtt(thetd(k),rtd(k))
      ENDDO

      zd(2)=zd(3)+(thvd(3)+thvd(2))*.5*(pid(3)-pid(2))/grav
      zd(1)=zd(2)+(thvd(2)+thvd(1))*.5*(pid(2)-pid(1))/grav

      do k=1,npi3
         sigzr(k)=topt(i,j)+sigz(k)*rtgt(i,j)
      enddo

      call htint(npd,ud,zd,npi3,vvv,sigzr)
      call psfill(npi3,vvv,ui2,np1,np2,i,j)
      call htint(npd,vd,zd,npi3,vvv,sigzr)
      call psfill(npi3,vvv,vi2,np1,np2,i,j)
      call htint(npd,thetd,zd,npi3,vvv,sigzr)
      call psfill(npi3,vvv,ti2,np1,np2,i,j)
      call htint(npd,rd,zd,npi3,vvv,sigzr)
      call psfill(npi3,vvv,ri2,np1,np2,i,j)

      call htint(npd,pid,zd,npi3,vvv,sigzr)
      call psfill(npi3,vvv,pi2,np1,np2,i,j)
      do k=1,npi3
         pi2(i,j,k)=(pi2(i,j,k)/cpdry)**cpor*p00
      enddo

      DO KL=Npi3,1,-1
         IF(TI2(I,J,KL).LT.1E19) GOTO 402
      ENDDO
      STOP 'ST2-MISS'
      402 CONTINUE
      KL=KL+1
      DO K=KL,Npi3
         TI2(I,J,K)=tempd(npd)*(p00/pi2(i,j,k))**rocp
         UI2(I,J,K)=UD(NPD)
         VI2(I,J,K)=VD(NPD)
         RI2(I,J,K)=0.
      ENDDO

      ! Find a pressure b.c. as second level above 700 mb

      pbc=70000.
      DO K=1,npd
         if(ppd(k).lt.pbc) then
            kpbc=k
            goto 320
         endif
      ENDDO
      320 continue
      DO K=2,npi3
         if(ti2(i,j,k).gt.thetd(kpbc)) then
            kibc=k
            goto 321
         endif
      ENDDO
      print*,'ISAN error: domain top not high enough'
      stop 'vterpp_s'
      321 continue

      do k=1,npi3
         vvv(k)=ti2(i,j,k)*(pi2(i,j,k)/p00)**rocp
         raux  =ptrh2rvapil(ri2(i,j,k),pi2(i,j,k),vvv(k),.false.)
         vvv(k)=virtt(ti2(i,j,k),raux)
      enddo
      raux = ptrh2rvapil(rd(kpbc),ppd(kpbc),tempd(kpbc),.false.)
      thvp=virtt(thetd(kpbc),raux)


      piibc=cpdry*pkd(kpbc)/p00**rocp
      pi2(i,j,kibc-1)=piibc+(zd(kpbc)-sigzr(kibc-1))*grav/(.5*(thvp+vvv(kibc-1)))
      do k=kibc-2,1,-1
         pi2(i,j,k)=pi2(i,j,k+1)+(sigzr(k+1)-sigzr(k))*grav/(.5*(vvv(k+1)+vvv(k)))
      enddo

      do k=kibc,npi3
         pi2(i,j,k)=pi2(i,j,k-1)-(sigzr(k)-sigzr(k-1))*grav/(.5*(vvv(k-1)+vvv(k)))
      enddo
      do k=1,npi3
         pi2(i,j,k)=(pi2(i,j,k)/cpdry)**cpor*p00
      enddo

      4500 continue
   ENDDO
ENDDO

! If any level is above 100 mb, compute geostrophic winds

!GDKM=2.*SPCON
!DO J=2,NPRY-1
!   GDLAT=(XSWLAT+(J-1)*GDATDY)*PI180
!   FCORI=1./(2.*7.292E-5*SIN(GDLAT))
!   DO I=2,NPRX-1
!      DO K=NISN,1,-1
!         IF(PI2(I,J,K).LT.10000.) THEN
!            UI2(I,J,K)=-FCORI*(SI2(I,J+1,K)-SI2(I,J-1,K))/(GDKM*GDATDY)
!            VI2(I,J,K)= FCORI*(SI2(I+1,J,K)-SI2(I-1,J,K))  &
!                        /(GDKM*GDATDX*COS(GDLAT))
!         ENDIF
!      ENDDO
!   ENDDO
!ENDDO


return
end

!***************************************************************************

subroutine psfill(nz,vvv,aaa,np1,np2,i,j)
implicit none
integer :: nz,np1,np2,i,j,k
real :: vvv(nz),aaa(np1,np2,nz)

do k=1,nz
   aaa(i,j,k)=vvv(k)
enddo

return
end
