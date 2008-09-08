!--------------------------------------------------------------------
!--------------------------------------------------------------------
subroutine cup_kbcon_catt(cap_inc,iloop, k22, kbcon, he_cup, hes_cup, &
     hkb,kzi,mix, mgmxp, mkx, mgmzp, istart, iend, ierr, kbmax, p_cup,&
     cap_max,j)
  implicit none
  integer i,j,k,kp, mix, mgmxp, mkx, mgmzp, istart, iend, iloop
  integer kbcon(mgmxp), k22(mgmxp), ierr(mgmxp), kbmax(mgmxp),kzi(mgmxp)
  real he_cup(mgmxp,mgmzp), hes_cup(mgmxp,mgmzp), p_cup(mgmxp,mgmzp)
  real pbcdif, cap_max(mgmxp),plus,cap_inc(mgmxp), hkb(mgmxp),hkbpbl
!  real pbcdif, cap_max(mgmxp),plus,cap_inc, hkb(mgmxp),hkbpbl

  !--- Determine the level of convective cloud base  - KBCON
  ! kbcon  = LFC of parcel from k22
  ! k22    = updraft originating level
  ! hkb    = moist static energy at originating level

  do I=ISTART,IEND
     kbcon(i)=1
     if (ierr(I).ne.0) cycle
     KBCON(I)=K22(I)
     
     GO TO 32
31   continue
     KBCON(I)=KBCON(I)+1
          
     if (KBCON(I).gt.KBMAX(i)+2) then
         if(iloop.lt.4)ierr(i)=3
         if(iloop.eq.4)ierr(i)=997

        cycle
     endif
32   continue

     if (HE_cup(I,K22(I)).lt.HES_cup(I,KBCON(I))) then
      
      GO TO 31
     else
    endif

     !     Cloud base pressure and max moist static energy pressure
     !     i.e., the depth (in mb) of the layer of negative buoyancy

     if (KBCON(I)-K22(I).eq.1) then
        !GO TO 27
        cycle
     endif

     PBCDIF = -P_cup(I,KBCON(I)) + P_cup(I,K22(I))
     plus=cap_max(i)-float(iloop-1)*cap_inc(i)
     plus=max(25.,cap_max(i)-float(iloop-1)*cap_inc(i))

     if (PBCDIF.gt.plus) then
        K22(I)   = K22(I)+1
        KBCON(I) = K22(I)         
        GO TO 32
     endif
  enddo

  return
end subroutine cup_kbcon_catt

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
