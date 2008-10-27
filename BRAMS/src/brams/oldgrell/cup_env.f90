!--------------------------------------------------------------------
subroutine cup_env(j,z,qes,he,hes,t,q,p,z1,mix,mgmxp,mkx,mgmzp,istart,iend     &
                  ,psur,ierr,tcrit,itest)

  use rconstants, only : rgas, cp, alvl, aklv, akiv, ep, g, rocp
  use therm_lib, only: virtt

  implicit none
  integer                         :: mix,mgmxp,mkx,mgmzp,i,k,istart,iend,iph,j &
                                    ,m,itest
  real,    dimension(mgmxp,mgmzp) :: z,p,qes,he,hes,t,q,tv
  real,    dimension(mgmxp)       :: psur,z1
  real,    dimension(2)           :: ae, be
  real                            :: tcrit, e, tvbar
  integer, dimension (mgmxp)      :: ierr

  BE(1)=ep*aklv/rocp
  AE(1)=BE(1)/273.+ALOG(610.71)
  BE(2)=ep*akiv/rocp
  AE(2)=BE(2)/273.+ALOG(610.71)
  do K=1,MKX
     do I=ISTART,IEND
        if(ierr(i).eq.0)then
           !sgb - IPH is for phase, dependent on TCRIT (water or ice)
           IPH = 1
           if (T(I,K).le.TCRIT) IPH = 2
           E = exp(AE(IPH)-BE(IPH)/T(I,K))
           QES(I,K) = ep*E/(100.*P(I,K)-E)
           if (QES(I,K) .le. 1.E-08)   QES(I,K)=1.E-08
           if (Q(I,K)   .gt. QES(I,K)) Q(I,K)=QES(I,K)
           TV(I,K) = virtt(T(I,K),Q(I,K))
        endif
     enddo
  enddo

  !--- z's are calculated with changed h's and q's and t's
  !--- if itest=2

  if (itest.ne.2) then
     do I=ISTART,IEND
        if (ierr(i).eq.0) then
           Z(I,1)=max(0.,Z1(I))-(ALOG(P(I,1))-ALOG(PSUR(I)))*rgas*TV(I,1)/g
        endif
     enddo

     !--- Calculate heights
     do K=2,MKX
        do I=ISTART,IEND
           if (ierr(i).eq.0) then
              TVBAR =  .5*TV(I,K)+.5*TV(I,K-1)
              Z(I,K) = Z(I,K-1)-(ALOG(P(I,K))-ALOG(P(I,K-1)))*rgas*TVBAR/g
           endif
        enddo
     enddo
  else
     do k=1,mkx
        do i=istart,iend
           if (ierr(i).eq.0) then
              z(i,k) = (he(i,k)-cp*t(i,k)-alvl*q(i,k))/g
              z(i,k) = max(1.e-3,z(i,k))
           endif
        enddo
     enddo
  endif

  !--- calculate moist static energy - HE
  !--- Saturated moist static energy - HES

  do K=1,MKX
     do I=ISTART,IEND
        if (ierr(i).eq.0) then
           if (itest.eq.0) HE(I,K) = g*Z(I,K)+cp*T(I,K)+alvl*Q(I,K)
           HES(I,K) = g*Z(I,K)+cp*T(I,K)+alvl*QES(I,K)

           if (HE(I,K).ge.HES(I,K)) HE(I,K) = HES(I,K)

        endif
     enddo
  enddo

  return
end subroutine cup_env

!--------------------------------------------------------------------
subroutine cup_env_clev(j,t,qes,q,he,hes,z,p,qes_cup,q_cup,he_cup,hes_cup      &
                       ,z_cup,p_cup,gamma_cup,t_cup,psur,mix,mgmxp,mkx,mgmzp   &
                       ,istart,iend,ierr,z1)

  use rconstants, only : alvl, rm, aklv
  implicit none

  integer                         :: i, j, k, mix, mgmxp, mkx, mgmzp, istart   &
                                    ,iend, m
  integer, dimension(mgmxp)       :: ierr
  real,    dimension(mgmxp)       :: psur,z1
  real,    dimension(mgmxp,mgmzp) :: qes_cup,q_cup,he_cup,hes_cup,z_cup,p_cup  &
                                    ,gamma_cup,t_cup,qes,q,he,hes,z,p,t

  do k=2,mkx
     do i=istart,iend
        if (ierr(i).eq.0)then
           qes_cup(i,k) = .5*(qes(i,k-1) + qes(i,k))
           q_cup(i,k)   = .5*(  q(i,k-1) +   q(i,k))
           hes_cup(i,k) = .5*(hes(i,k-1) + hes(i,k))
           he_cup(i,k)  = .5*( he(i,k-1) +  he(i,k))
           if (he_cup(i,k).gt.hes_cup(i,k)) he_cup(i,k) = hes_cup(i,k)

           z_cup(i,k) = .5*(z(i,k-1) + z(i,k))
           p_cup(i,k) = .5*(p(i,k-1) + p(i,k))
           t_cup(i,k) = .5*(t(i,k-1) + t(i,k))

           gamma_cup(i,k) =aklv*(alvl/(rm*t_cup(i,k)*t_cup(i,k)))*qes_cup(i,k)
        endif
     enddo
  enddo
  do i=istart,iend
     if (ierr(i).eq.0)then
        qes_cup(i,1) = qes(i,1)
        q_cup(i,1)   = q(i,1)
        hes_cup(i,1) = hes(i,1)
        he_cup(i,1)  = he(i,1)

        !srf
        z_cup(i,1)   = z1(i)
        p_cup(i,1)   = psur(i)
        t_cup(i,1)   = t(i,1)
        !srf	
        gamma_cup(i,1) = aklv*(alvl/(rm*t_cup(i,1)*t_cup(i,1)))*qes_cup(i,1)
     endif
  enddo

  return
end subroutine cup_env_clev

!--------------------------------------------------------------------
subroutine cup_direction2(i,j,dir,id,mix,mjx,mgmxp,mgmyp,massflx,iresult,num   &
                         ,imass,nall,maxens3,massfld)
  implicit none
  integer                         :: mix,mjx,mgmxp,mgmyp,i,j,k,num,iresult     &
                                    ,maxens3,imass,nall,ia,ib,ja,jb
  integer, dimension(mgmxp,mgmyp) :: id
  real                            :: massfld,diff
  real,    dimension(mgmxp)       :: dir
  real,    dimension(mgmxp,mgmyp) :: massflx
  
  if (imass.eq.1) then
     massfld = massflx(i,j)
  endif
  iresult=0
  !      return

  diff = 22.5
  if (dir(i).lt.22.5) dir(i)=360.+dir(i)
  if (id(i,j).eq.1)   iresult=1
  ja=j-1
  ia=i-1
  jb=j+1
  ib=i+1
  if (dir(i).gt.90.-diff.and.dir(i).le.90.+diff) then

     !--- Steering flow from the east
     if (id(ib,j).eq.1) then
        iresult=1
        if (imass.eq.1) then
           massfld = max(massflx(ib,j),massflx(i,j))
        endif
        return
     endif
  else if (dir(i).gt.135.-diff.and.dir(i).le.135.+diff)then

     !--- Steering flow from the south-east
     if (id(ib,ja).eq.1) then
        iresult=1
        if (imass.eq.1) then
           massfld = max(massflx(ib,ja),massflx(i,j))
        endif
        return
     endif

     !--- Steering flow from the south
  else if (dir(i).gt.180.-diff.and.dir(i).le.180.+diff) then
     if (id(i,ja).eq.1) then
        iresult=1
        if (imass.eq.1) then
           massfld = max(massflx(i,ja),massflx(i,j))
        endif
        return
     endif

     !--- Steering flow from the south west
  else if (dir(i).gt.225.-diff.and.dir(i).le.225.+diff) then
     if (id(ia,ja).eq.1) then
        iresult=1
        if (imass.eq.1) then
           massfld = max(massflx(ia,ja),massflx(i,j))
        endif
        return
     endif

     !--- Steering flow from the west
  else if (dir(i).gt.270.-diff.and.dir(i).le.270.+diff) then
     if (id(ia,j).eq.1) then
        iresult=1
        if (imass.eq.1)then
           massfld = max(massflx(ia,j),massflx(i,j))
        endif
        return
     endif

     !--- Steering flow from the north-west
  else if (dir(i).gt.305.-diff.and.dir(i).le.305.+diff) then
     if (id(ia,jb).eq.1) then
        iresult=1
        if (imass.eq.1) then
           massfld = max(massflx(ia,jb),massflx(i,j))
        endif
        return
     endif

     !--- Steering flow from the north
  else if (dir(i).gt.360.-diff.and.dir(i).le.360.+diff) then
     if (id(i,jb).eq.1) then
        iresult=1
        if (imass.eq.1) then
           massfld = max(massflx(i,jb),massflx(i,j))
        endif
        return
     endif

     !--- Steering flow from the north-east
  else if (dir(i).gt.45.-diff.and.dir(i).le.45.+diff) then
     if (id(ib,jb).eq.1) then
        iresult=1
        if (imass.eq.1) then
           massfld = max(massflx(ib,jb),massflx(i,j))
        endif
        return
     endif
  endif

  return
end subroutine cup_direction2

!--------------------------------------------------------------------
subroutine cup_ktop(ilo,dby,kbcon,ktop,mix,mgmxp,mkx,mgmzp,istart,iend,ierr)
  implicit none
  integer                         :: mix,mgmxp,mkx,mgmzp,i,k,istart,iend,ilo
  integer, dimension(mgmxp)       :: ierr,kbcon,ktop
  real,    dimension(mgmxp,mgmzp) :: dby

  do I=ISTART,IEND
     ktop(i)=1
     if (ierr(I).eq.0) then

        do K=KBCON(I)+1,MKX-2
           if (DBY(I,K).le.0.) then
              KTOP(I) = K-1
              GO TO 41
           endif
        enddo
        if (ilo.eq.1) ierr(i)=5
        if (ilo.eq.2) ierr(i)=998
        GO TO 42
41      continue
        do k=ktop(i)+1,mkx
           dby(i,k)=0.
        enddo
     endif
42   continue
  enddo
  return
end subroutine cup_ktop

!--------------------------------------------------------------------
subroutine cup_kbcon(iloop,k22,kbcon,he_cup,hes_cup,mix,mgmxp,mkx,mgmzp,istart &
                    ,iend,ierr,kbmax,p_cup,cap_max)
  implicit none
  integer                         :: i,mix,mgmxp,mkx,mgmzp,istart,iend,iloop
  integer, dimension(mgmxp)       :: kbcon,k22,ierr,kbmax
  real                            :: pbcdif
  real,    dimension(mgmxp)       :: cap_max
  real,    dimension(mgmxp,mgmzp) :: he_cup,hes_cup,p_cup

  !--- Determine the level of convective cloud base  - KBCON

  mainloop: do I=ISTART,IEND
     kbcon(i)=1
     if (ierr(I).ne.0) cycle mainloop
     KBCON(I)=K22(I)
     GO TO 32
31   continue
     KBCON(I)=KBCON(I)+1
     if (KBCON(I).gt.KBMAX(i)+2) then
        if(iloop.eq.1)ierr(i)=3
        if(iloop.eq.2)ierr(i)=997
        cycle mainloop
     endif
32   continue
     if (HE_cup(I,K22(I)).lt.HES_cup(I,KBCON(I))) GO TO 31

     !     Cloud base pressure and max moist static energy pressure
     !     i.e., the depth (in mb) of the layer of negative buoyancy

     if (KBCON(I)-K22(I).eq.1) then
        cycle mainloop
     endif

     PBCDIF = -P_cup(I,KBCON(I)) + P_cup(I,K22(I))
     if (PBCDIF.gt.cap_max(i)) then
        K22(I)   = K22(I)+1
        KBCON(I) = K22(I)
        GO TO 32
     endif
  enddo mainloop

  return
end subroutine cup_kbcon


!--------------------------------------------------------------------
subroutine cup_kbcon_cin(iloop,k22,kbcon,he_cup,hes_cup,z,tmean,qes,mix,mgmxp  &
                        ,mkx,mgmzp,istart,iend,ierr,kbmax,p_cup,cap_max)
  use rconstants, only : alvl,aklv,rm,cp,g
  implicit none
  integer                         :: i,mix,mgmxp,mkx,mgmzp,istart,iend,iloop
  integer, dimension(mgmxp)       :: kbcon,k22,ierr,kbmax
  real                            :: pbcdif,cin,cin_max,dh,tprim,gamma
  real,    dimension(mgmxp)       :: cap_max
  real,    dimension(mgmxp,mgmzp) :: he_cup,hes_cup,p_cup,z,tmean,qes

  !
  !--- DETERMINE THE LEVEL OF CONVECTIVE CLOUD BASE  - KBCON
  !

  do I=ISTART,IEND
     cin_max=-cap_max(i)

     kbcon(i)=1
     cin = 0.
     if (ierr(I).ne.0) cycle
     KBCON(I)=K22(I)
     GO TO 32
31   continue
     KBCON(I)=KBCON(I)+1
     if (KBCON(I).gt.KBMAX(i)+2) then
        if (iloop.eq.1) ierr(i)=3
        if (iloop.eq.2) ierr(i)=997
        cycle
     endif
32   continue
     dh = HE_cup(I,K22(I)) - HES_cup(I,KBCON(I))
     if (dh.lt. 0.) then
        GAMMA = aklv*(alvl/(rm*(Tmean(I,K22(i))**2)))*QES(I,K22(i))
        tprim = dh/(cp*(1.+gamma))

        cin   = cin + g*tprim*(z(i,k22(i))-z(i,k22(i)-1))/tmean(i,k22(i))
        go to 31
     end if

     !     If negative energy in negatively buoyant layer
     !       exceeds convective inhibition (CIN) threshold,
     !       then set K22 level one level up and see if that
     !       will work.

     if (cin.lt.cin_max) then
        K22(I)=K22(I)+1
        KBCON(I)=K22(I)
        GO TO 32
     endif
  enddo

  return
end subroutine cup_kbcon_cin


!--------------------------------------------------------------------
subroutine MINIMI(ARRAY,MXX,mgmxp,MZX,mgmzp,KS,KEND,KT,ISTART,IEND,ierr)

  implicit none

  integer                          :: MXX,MZX,ISTART,IEND,KSTOP,I,K,mgmxp,mgmzp
  integer, dimension(mgmxp)        :: KT,KS,KEND,ierr
  real,    dimension(mgmxp)        :: X
  real,    dimension(mgmxp,mgmzp)  :: ARRAY

  do I=ISTART,IEND
     KT(I)=KS(I)
     if (ierr(i).eq.0) then
        X(I)=ARRAY(I,KS(I))
        KSTOP=max(KS(I)+1,KEND(I))

        do K=KS(I)+1,KSTOP
           if(ARRAY(I,K).lt.X(I)) then
              X(I)=ARRAY(I,K)
              KT(I)=K
           endif
        enddo
     endif
  enddo

  return
end subroutine MINIMI

!--------------------------------------------------------------------
subroutine MAXIMI(ARRAY,MXX,mgmxp,MZX,mgmzp,KS,KE,MAXX,ISTART,IEND,ierr)

  implicit none

  integer                         :: MXX,MZX,KS,ISTART,IEND,I,K,mgmxp,mgmzp
  integer, dimension(mgmxp)       :: MAXX,ierr(mgmxp),KE(mgmxp)
  real,    dimension(mgmxp)       :: X
  real,    dimension(mgmxp,mgmzp) :: ARRAY
  real XAR

  do I=ISTART,IEND
     MAXX(I)=KS
     if(ierr(i).eq.0)then
        X(I)=ARRAY(I,KS)

        do K=KS,KE(i)
           XAR=ARRAY(I,K)
           if(XAR.ge.X(I)) then
              X(I)=XAR
              MAXX(I)=K
           endif
        enddo
     endif
  enddo

  return
end subroutine MAXIMI


!-------------------------------------------------------------------

subroutine get_zi(mix,mgmxp,mkx,mgmzp,istart,iend,j,ierr,kzi,tkeg, &
                  rcpg,z,ztop,tkmin)

  implicit none
  integer mix,mgmxp,mkx,mgmzp,i,k,istart,iend,j,kzimax,ktke_max
  real tkmin,rcpmin,pblhmax,tke_tmp
  real,    dimension(mgmxp,mgmzp) :: tkeg,rcpg,z
  real,    dimension(mgmxp)	  :: ztop
  integer, dimension(mgmxp)	  :: kzi,ierr

  data rcpmin/1.e-5/, pblhmax/3000./ 
  !print*,j,mgmxp,mgmzp,mix,istart,iend
  kzimax=2
  do i=istart,iend
    kzi(i)  = 2

    if(ierr(i).eq.0)then
     tke_tmp = 0.
     ktke_max= 1
     !---  max level for kzi
     do k=1,mkx
       if(z(i,k).ge. pblhmax+ztop(i)) then	
          kzimax = k
          !print*,z(i,k), pblhmax,ztop(i),kzimax
          exit
       endif
     enddo
     !---
     !       go to 201
     !level of max tke  below kzimax and w/out clouds
     do  k=1,kzimax
       !print*,k,tkeg(i,k), tke_tmp,ktke_max,kzimax
       if(rcpg(i,k) .lt. rcpmin) then
	 if( tkeg(i,k) .ge. tke_tmp) then
	   tke_tmp = tkeg(i,k)
	   cycle
	 else
	   ktke_max= max(1,k-1)
	   exit
	 endif
       endif	   
     enddo	     
     !201    continue
!         print*,ktke_max

     do k=ktke_max,kzimax+1
!  	print*,rcpg(i,k),tkeg(i,k),k,kzi(i),i
	if(rcpg(i,k) .lt. rcpmin) then
          if(tkeg(i,k) .gt. 1.1*tkmin)  then
	    kzi(i) = k
	    cycle
	  endif
        else
	   kzi(i) = k
	   exit
	endif
     enddo
     kzi(i) = max(2     ,kzi(i))
     kzi(i) = min(kzimax,kzi(i))

   endif
 enddo
end subroutine get_zi
