
!*******************************************************************************
!############################# Change Log ##################################
! 2.3.0.1
!
! 000830 CJT rams_comp ##
!            Corrected reference to "cpr" and "r" to "cpor" and "rgas" ##
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!  Mission Research Corporation / *ASTeR Division
!###########################################################################

subroutine RAMS_comp(n1,n2,n3,n4,n5,n6)

use somevars
use rconstants
use rpost_coms
use therm_lib, only : qtk, dewfrostpoint, rslif, virtt, thetaeiv

dimension a(n1,n2,n3),b(n1,n2,n3),c(n1,n2,n3),d(n1,n2,n3),e(n1,n2,n3),o(n1,n2,n3),topt(n1,n2)
dimension a2(n1,n2,n4,n5),a6(n1,n2,n3,n6)
real f1(n1,n2,n3),f2(n1,n2,n3)
dimension theta(n1,n2,n3),pp(n1,n2,n3),slp(n1,n2),z(n1,n2,n3)
real :: temptemp,fracliq
integer :: nsoil

entry RAMS_transf_ppb_day(n1,n2,n3,a)
!  TRANSFORMACAO DE kg[CO]/kg[AR] para ppb (PARTE POR BILHAO)
     do k=1,n3
 	do j=1,n2
 	   do i=1,n1

     a(i,j,k)=a(i,j,k)*(mmdry/mmco2)*1.E+9*1.e-6*day_sec

 	   enddo
 	enddo
     enddo
return

entry RAMS_transf_ppb(n1,n2,n3,a)
!  TRANSFORMACAO DE kg[CO]/kg[AR] para ppb (PARTE POR BILHAO)
     do k=1,n3
 	do j=1,n2
 	   do i=1,n1

     a(i,j,k)=a(i,j,k)*(mmdry/mmco)*1.E+9

 	   enddo
 	enddo
     enddo
return

entry RAMS_transf_ppm(n1,n2,n3,a)
!  TRANSFORMACAO DE kg[CO2]/kg[AR] para ppm (PARTE POR MILHAO)
     do k=1,n3
 	do j=1,n2
 	   do i=1,n1
 	      a(i,j,k)=a(i,j,k)*(mmdry/mmco2)*1.E+6
 	   enddo
 	enddo
     enddo
return

entry RAMS_transf_ugm3(n1,n2,n3,a,b)
!   Transformacao de conc  de kg/kg para ug/m3
    do k=1,n3
      do j=1,n2
    	do i=1,n1
        !print*,a(i,j,k),b(i,j,k)
	a(i,j,k)=a(i,j,k)*b(i,j,k)*1.e+9
    
    	enddo
      enddo
    enddo
return

entry get_ZI(n1,n2,n3,a,c,ngrd)
tkemin =1.E-2
      
      do i=1,n1
         do j=1,n2
          a(i,j,1) = (ztn(2,ngrd)+ ztn(1,ngrd))
         enddo
      enddo

      do i=1,n1
        do j=1,n2
         do k = 2,n3-1
!
!         print*,k,c(i,j,k),ztn(k,ngrd)
         if(c(i,j,k).lt.tkemin) then
            kzi = k
            a(i,j,1)=0.5*(ztn(kzi,ngrd)+ ztn(kzi-1,ngrd))
	    go to 500
	 endif
	 enddo
 500     continue
!         print*,i,j,a(i,j,1)	  
        enddo
      enddo
return

entry RAMS_comp_nebulosidade(n1,n2,n3,a,b,c,e,ngrd)
!
! b=cloud
! c=dn0
! e=topo
       do j=1,n2
           do i=1,n1
         rnebu=0.
         rodzint=0
         do k=2,n3-1
          dz=(ztn(k,ngrd)-ztn(k-1,ngrd))*(1.-e(i,j,1)/zmn(nnzp(1)-1,1))
          rnebu   = rnebu   + b(i,j,k)*c(i,j,k)*dz	  
          rodzint = rodzint +	       c(i,j,k)*dz	  
         enddo
              a(i,j,1) = 1000.*rnebu/(rodzint+1.e-5)
	 enddo
       enddo
return

!entry RAMS_get_leaf(n1,n2,n5,a)
!       do j=1,n2
!        do i=1,n1
!         do ip=1,n5
!            a(i,j,ip)=a(i,j,ip)
!         enddo
!        enddo
!       enddo
!return

entry RAMS_comp_vegclass2(n1,n2,n3,a)
      do i=1,n1
      a(i,1,1) = a(i,1,1) + .5
      enddo
return

!entry RAMS_comp_vegclass(n1,n2,n5,a)
!       do j=1,n2
!        do i=1,n1
!         do ip=1,n5
!!            print*,i,j,ip,' veg=',a(i,j,ip)
!            a(i,j,ip)=float(int(a(i,j,ip)+0.5))
!!            print*,i,j,ip,' veg=',a(i,j,ip)
!         enddo
!        enddo
!       enddo
!return

entry D3toD2(n1,n2,n3,klevel,a,c)
  do j=1,n2
      do i=1,n1
	a(i,j,1)=c(i,j,klevel)
!	print*,i,j,klevel,c(i,j,klevel)
      enddo
  enddo
return


entry up_to_tp(n1,n2,n3,a,b)
         do k=1,n3
            do j=1,n2
               a(1,j,k)=b(1,j,k)
               do i=2,n1
               a(i,j,k)=0.5* (b(i,j,k) + b(i-1,j,k))
               enddo
            enddo
         enddo
return

entry vp_to_tp(n1,n2,n3,a,b)
         do k=1,n3
            do i=1,n1
               a(i,1,k)=b(i,1,k)
               do j=2,n2
               a(i,j,k)=0.5* (b(i,j,k) + b(i,j-1,k))
               enddo
            enddo
         enddo
return

entry wp_to_tp(n1,n2,n3,a,b)
         do j=1,n2
           do i=1,n1
              a(i,j,1)=b(i,j,1)
              do k=2,n3
               a(i,j,k)=0.5* (b(i,j,k) + b(i,j,k-1))
              enddo
            enddo
         enddo
return

entry set_undef(n1,n2,n3,a,c)
  do j=1,n2
      do i=1,n1
	a(i,j,1)=c(i,j,1)
!        if(a(i,j,1) .lt. 0.0001) a(i,j,1) = -9.99e33
!	print*,i,j,klevel,c(i,j,klevel)
      enddo
  enddo
return
!SRF
!***********************************************************

entry RAMS_comp_vals(n1,n2,n3,a)
   print*,'==================== values =========================='
   do k=1,n3
      do j=1,n2
         do i=1,n1
            if(a(i,j,k).ne.0.) write(*,'(3i3,e14.6)') i,j,k  &
                                                     ,a(i,j,k)
         enddo
      enddo
   enddo
   print*,'======================================================'
return

entry RAMS_comp_tot(n1,n2,n3,a)
   tot=0.
   do k=1,n3
      do j=1,n2
         do i=1,n1
            tot=tot+a(i,j,k)
         enddo
      enddo
   enddo
   write(*,'(a,e12.6)') '-> total- ',tot
return

entry RAMS_comp_maxval(n1,n2,n3,a)
   zmax=-1.0e30
   do k=1,n3
      do j=1,n2
         do i=1,n1
            if(a(i,j,k).gt.zmax) then
               zmax=a(i,j,k)
               maxx=i
               maxy=j
               maxz=k
            endif
         enddo
      enddo
   enddo
   write(*,'(a,e12.6,a,3i3)') '-> max- ',zmax,' at i,j,k-',maxx  &
                              ,maxy,maxz
return

entry RAMS_comp_minval(n1,n2,n3,a)
   zmin=1.0e30
   do k=1,n3
      do j=1,n2
         do i=1,n1
            if(a(i,j,k).lt.zmin) then
               zmin=a(i,j,k)
               minx=i
               miny=j
               minz=k
            endif
         enddo
      enddo
   enddo
   write(*,'(a,e12.6,a,3i3)') '-> min- ',zmin,' at i,j,k-',minx  &
                              ,miny,minz
return

entry RAMS_comp_zero(n1,n2,n3,a)
   do k=1,n3
      do j=1,n2
         do i=1,n1
            a(i,j,k)=0.
         enddo
      enddo
   enddo
return

entry RAMS_comp_1minus(n1,n2,n3,a)
   do k=1,n3
      do j=1,n2
         do i=1,n1
            a(i,j,k)=1.-a(i,j,k)
         enddo
      enddo
   enddo
return

entry RAMS_comp_mults(n1,n2,n3,a,s)
   do k=1,n3
      do j=1,n2
         do i=1,n1
            a(i,j,k)=a(i,j,k) * s
         enddo
      enddo
   enddo
return


entry RAMS_comp_accum(n1,n2,n3,a,b)
   do k=1,n3
      do j=1,n2
         do i=1,n1
            a(i,j,k)=a(i,j,k)+b(i,j,k)
         enddo
      enddo
   enddo
return

entry RAMS_comp_noneg(n1,n2,n3,a)
   do k=1,n3
      do j=1,n2
         do i=1,n1
            a(i,j,k)=max(a(i,j,k),0.)
         enddo
      enddo
   enddo
return

entry RAMS_comp_nopos(n1,n2,n3,a)
   do k=1,n3
      do j=1,n2
         do i=1,n1
            a(i,j,k)=min(a(i,j,k),0.)
         enddo
      enddo
   enddo
return

entry RAMS_comp_subt(n1,n2,n3,a,b)
   do k=1,n3
      do j=1,n2
         do i=1,n1
            a(i,j,k)=a(i,j,k)-b(i,j,k)
         enddo
      enddo
   enddo
return

entry RAMS_comp_mult(n1,n2,n3,a,b)
   do k=1,n3
      do j=1,n2
         do i=1,n1
            a(i,j,k)=a(i,j,k)*b(i,j,k)
         enddo
      enddo
   enddo
return

entry RAMS_comp_z(n1,n2,n3,a,c,ngrd)
   do k=1,n3
      do j=1,n2
         do i=1,n1
            a(i,j,k)=c(i,j,1)  &
                 +ztn(k,ngrd)*(1.-c(i,j,1)/zmn(nnzp(1)-1,1))
         enddo
      enddo
   enddo
return

entry RAMS_comp_rotate(n1,n2,n3,a,b,ngrd)
   do k=1,n3
      do j=1,n2
         do i=1,n1
            call xy_ll(qlat,qlon,platn(ngrd),plonn(ngrd)  &
               ,xtn(i,ngrd),ytn(j,ngrd))
            u=a(i,j,k)
            v=b(i,j,k)
            call uvtoueve(u,v,a(i,j,k),b(i,j,k)  &
                         ,qlat,qlon,platn(ngrd),plonn(ngrd))
         enddo
      enddo
   enddo
return

entry RAMS_comp_tempK(n1,n2,n3,a,b)
   do k=1,n3
      do j=1,n2
         do i=1,n1
            a(i,j,k)=a(i,j,k)*b(i,j,k)/cp
         enddo
      enddo
   enddo
return

entry RAMS_comp_press(n1,n2,n3,a)
   do k=1,n3
      do j=1,n2
         do i=1,n1
            a(i,j,k)=(a(i,j,k)/cp)**cpor*p00*.01
         enddo
      enddo
   enddo
return

entry RAMS_comp_wcms(n1,n2,n3,a)
   do k=1,n3
      do j=1,n2
         do i=1,n1
            a(i,j,k)=a(i,j,k)*100.
         enddo
      enddo
   enddo
return

entry RAMS_comp_avgw(n1,n2,n3,a)
   do k=n3,2,-1
      do j=1,n2
         do i=1,n1
            a(i,j,k)=0.5*(a(i,j,k)+a(i,j,k-1))
         enddo
      enddo
   enddo
return

entry RAMS_comp_avgu(n1,n2,n3,a)
   do k=1,n3
      do j=1,n2
         do i=n1,2,-1
            a(i,j,k)=0.5*(a(i,j,k)+a(i-1,j,k))
         enddo
      enddo
   enddo
return

entry RAMS_comp_avgv(n1,n2,n3,a)
   do k=1,n3
      do j=n2,2,-1
         do i=1,n1
            a(i,j,k)=0.5*(a(i,j,k)+a(i,j-1,k))
         enddo
      enddo
   enddo
return

entry RAMS_comp_sfcdiv(n1,n2,n3,a,ngrd)
   do j=1,n2
      do i=1,n1
         a(i,j,1)=-(a(i,j,2)-a(i,j,1))*dztn(2,ngrd)
      enddo
   enddo
return

entry RAMS_comp_rt(n1,n2,n3,a)
   do k=1,n3
      do j=1,n2
         do i=1,n1
            a(i,j,k)=max(0.,a(i,j,k))*1000.
         enddo
      enddo
   enddo
return

entry RAMS_comp_hr_pcprate(n1,n2,n3,a,b,c)
   do k=1,n3
      do j=1,n2
         do i=1,n1
            a(i,j,k)=(b(i,j,k)-c(i,j,k))
         enddo
      enddo
   enddo
return

entry RAMS_comp_speed(n1,n2,n3,a,b)
   do k=1,n3   
      do j=1,n2
         do i=1,n1
            a(i,j,k)=sqrt(a(i,j,k)**2+b(i,j,k)**2)
         enddo
      enddo
   enddo
return

entry RAMS_comp_dir(n1,n2,n3,a,b,ngrd)
   do k=1,n3
      do j=1,n2
         do i=1,n1
            call xy_ll(qlat,qlon,platn(ngrd),plonn(ngrd)  &
               ,xtn(i,ngrd),ytn(j,ngrd))
            u=a(i,j,k)
            v=b(i,j,k)
            call uvtoueve(u,v,a(i,j,k),b(i,j,k)  &
                         ,qlat,qlon,platn(ngrd),plonn(ngrd))
            call winddf(a(i,j,k),ff,a(i,j,k),b(i,j,k))
         enddo
      enddo
   enddo
return

entry RAMS_comp_dewK(n1,n2,n3,a,b,c)
   do k=1,n3
      do j=1,n2
         do i=1,n1
            xpress=(b(i,j,k)/cp)**cpor*p00
            xtemp=c(i,j,k)*b(i,j,k)/cp
            xwatsat=rslif(xpress,xtemp)
            a(i,j,k)=dewfrostpoint(xpress,min(a(i,j,k),xwatsat) )
         enddo
      enddo
   enddo
return

!<Demerval
entry RAMS_comp_thete(n1,n2,n3,a,e,f1,f2)
   do k=1,n3
      do j=1,n2
         do i=1,n1
         TL=55.+1./( 1./(f1(i,j,k)-55.) -  &
            log(f2(i,j,k)/100.)/2840.)
            a(i,j,k)=a(i,j,k)*exp((alvl*e(i,j,k))/(cp*TL))
         enddo
      enddo
   enddo
return

!entry RAMS_comp_thete(n1,n2,n3,a,b,c)
!   do k=1,n3
!      do j=1,n2
!         do i=1,n1
!            xpress=(b(i,j,k)/cp)**cpor*p00
!            xtemp=c(i,j,k)*b(i,j,k)/cp
!            xwatsat=rslif(xpress,xtemp)
!            a(i,j,k)=c(i,j,k)*exp( alvl*xwatsat  &
!                 /(cp*dewfrostpoint(xpress,min(a(i,j,k),xwatsat) )) )
!         enddo
!      enddo
!   enddo
!return

!Demerval>

entry RAMS_comp_thetv(n1,n2,n3,a,b,c)
   do k=1,n3
      do j=1,n2
         do i=1,n1
            a(i,j,k)=virtt(a(i,j,k),b(i,j,k),c(i,j,k))
         enddo
      enddo
   enddo
return

entry RAMS_comp_bowen(n1,n2,n3,a,b)
   do k=1,n3
      do j=1,n2
         do i=1,n1
            a(i,j,k)=a(i,j,k)/max(1.e-12,b(i,j,k))*1004./2.5e6
         enddo
      enddo
   enddo
return

entry RAMS_comp_rh(n1,n2,n3,a,b,c)
   do k=1,n3
      do j=1,n2
         do i=1,n1
            xtemp=c(i,j,k)*b(i,j,k)/cp
            xpress=(b(i,j,k)/cp)**cpor*p00
            a(i,j,k)=100.*min(1.  &
                 ,max(0.,a(i,j,k)/rslif(xpress,xtemp)))
         enddo
      enddo
   enddo
return

entry RAMS_comp_watsat(n1,n2,n3,a,b,c)
   do k=1,n3
      do j=1,n2
         do i=1,n1
            c(i,j,k)=c(i,j,k)*a(i,j,k)/cp
            b(i,j,k)=(b(i,j,k)/cp)**cpor*p00
            a(i,j,k)=rslif(b(i,j,k),c(i,j,k))
         enddo
      enddo
   enddo
return


entry RAMS_comp_vegclass(n1,n2,n3,a)
   do i=1,n1
      a(i,1,1) = float(nint(a(i,1,1)))
!      print*,i,int(a(i,1,1)) 
   enddo
   return

entry RAMS_comp_horizdiv(n1,n2,n3,a)
   do k=2,n3
      do j=1,n2
         do i=1,n1
            a(i,j,k)=-(a(i,j,k)-a(i,j,k-1))*dztn(k,ngrd)
         enddo
      enddo
   enddo
return

entry RAMS_comp_vertint(n1,n2,n3,a,topt,ngrd)
   ztop = zmn(nnzp(1)-1,1)
   do j = 1,n2
      do i = 1,n1
         rtgt = 1. - topt(i,j) / ztop
         a(i,j,1) = 0.
         do k = 2,n3-1
            a(i,j,1) = a(i,j,1) + a(i,j,k) * (zmn(k,ngrd)-zmn(k-1,ngrd)) * rtgt
         enddo
      enddo
   enddo
return

entry RAMS_comp_ppress(n1,n2,n3,a,c)
   do k=1,n3
      do j=1,n2
         do i=1,n1
            a(i,j,k) = 1000. * (a(i,j,k)/cp) ** cpor  &
                     - 1000. * (c(i,j,k)/cp) ** cpor
         enddo
      enddo
   enddo
return

entry RAMS_comp_raintemp(n1,n2,n3,a)
   do k=1,n3
      do j=1,n2
         do i=1,n1
            if (a(i,j,k) > 0.) then
               a(i,j,k) = tsupercool + a(i,j,k) * cliqi - t00
            else
               a(i,j,k) = -9.99e33
            end if
         enddo
      enddo
   enddo
return

entry RAMS_comp_qtcpcp(n1,n2,n3,a)
   do k=1,n3
      do j=1,n2
         do i=1,n1
            if (a(i,j,k) > 0.) then
               call qtk(a(i,j,k),temptemp,fracliq)
               a(i,j,k) = temptemp - t00
            else
               a(i,j,k) = -9.99e33
            end if
         end do
      end do
   end do
return

entry RAMS_comp_fracliq(n1,n2,n3,a)
   do k=1,n3
      do j=1,n2
         do i=1,n1
            call qtk(a(i,j,k),temptemp,fracliq)
            a(i,j,k) = fracliq
         end do
      end do
   enddo
return

entry RAMS_comp_fracice(n1,n2,n3,a)
   do k=1,n3
      do j=1,n2
         do i=1,n1
            call qtk(a(i,j,k),temptemp,fracliq)
            a(i,j,k) = 1.0 - fracliq
         enddo
      enddo
   enddo
return

entry RAMS_comp_hydrodiam(n1,n2,n3,a,c,ccfmas,ppwmas)
   rpwmas = 1. / ppwmas
   do k=1,n3
      do j=1,n2
         do i=1,n1
            if(a(i,j,k) .gt. 1.e-10 .and. c(i,j,k).gt.1.e-10)then
               a(i,j,k) = (a(i,j,k) / (c(i,j,k) * ccfmas))**rpwmas
            else
               a(i,j,k) = 0.
            endif
         enddo
      enddo
   enddo
return

entry rams_sum_snowlayers(n1,n2,n3,a)
   do ip=1,n3
      do k=2,n2
         do ij=1,n1
            a(ij,1,ip) = a(ij,1,ip) + a(ij,k,ip)
         enddo
      enddo
   enddo
return

entry rams_fill_sst(n1,n2,n3,kp,a,c)
   do j=1,n2
      do i = 1,n1
         call qtk(c(i,j,kp),temptemp,fracliq)
         a(i,j,1) = temptemp-t00
      enddo
   enddo
return

entry rams_comp_pcpgnorm(n1,n2,n3,a,c)
return

entry rams_comp_vapnorm(n1,n2,n3,a,c)
return

entry rams_comp_snownorm(n1,n2,n3,a,c)
return

entry rams_comp_vegnorm(n1,n2,n3,a,c)
return

entry rams_comp_cannorm(n1,n2,n3,a,c)
return


entry acha_ztropop(n1,n2,n3,a,c,ngrd)
  a(1:n1,1:n2,1)=0.
  estratosfera=0.018 ! 15 K/Km
  do i=1,n1
    do j=1,n2
      malhakz: do k=3,n3
        dtheta  = (c(i,j,k)-c(i,j,k-1))/(ztn(k,ngrd)-ztn(k-1,ngrd))
        if (dtheta.gt.estratosfera) exit malhakz
      end do malhakz
      a(i,j,1)=0.5*(ztn(k,ngrd)+ztn(k-1,ngrd))
    end do
  end do
return

entry acha_ttropop(n1,n2,n3,a,c,e,ngrd)
  a(1:n1,1:n2,1)=0.
  estratosfera=0.018 ! 15 K/Km
  do i=1,n1
    do j=1,n2
      malhakt: do k=3,n3
        dtheta  = (c(i,j,k)-c(i,j,k-1))/(ztn(k,ngrd)-ztn(k-1,ngrd))
        if (dtheta.gt.estratosfera) exit malhakt
      end do malhakt
      a(i,j,1)=0.5*(e(i,j,k)+e(i,j,k-1))
    end do
  end do
return

entry acha_ptropop(n1,n2,n3,a,c,e,ngrd)
  a(1:n1,1:n2,1)=0.
  estratosfera=0.018 ! 15 K/Km
  do i=1,n1
    do j=1,n2
      malhakp: do k=3,n3
        dtheta  = (c(i,j,k)-c(i,j,k-1))/(ztn(k,ngrd)-ztn(k-1,ngrd))
        if (dtheta.gt.estratosfera) exit malhakp
      end do malhakp
      a(i,j,1)=exp(0.5*(log(e(i,j,k))+log(e(i,j,k-1))))
    end do
  end do
return


entry get_ZItheta(n1,n2,n3,a,c,e,ngrd)
      do i=1,n1
        do j=1,n2
          a(i,j,1) = 0.
        enddo
      enddo

      do i=1,n1
        do j=1,n2

          do k=3,n3-5            
            dtheta=c(i,j,k)-c(i,j,k-1)
            x=0.6
            if(dtheta.gt.x.or.abs(e(i,j,k)).gt.1.e-5) go to 1881
          enddo
 1881     continue
          kzi = k
          if(abs(e(i,j,k)).gt.1.e-5) then

            a(i,j,1) =0.5*(ztn(kzi,ngrd)+ ztn(kzi-1,ngrd))

          else
        
            a(i,j,1) =0.5*(ztn(kzi,ngrd)+ ztn(kzi-1,ngrd))
          endif

          if(a(i,j,1).lt.0..or.kzi.le.2) a(i,j,1)=0.
       enddo
     enddo

   return



entry RAMS_comp_pbl (n1,n2,n3,a,c,ngrd)

tkethrsh=0.001   ! tke threshold for PBL height in m2/s2
do j=1,n2
   do i=1,n1
      do k=1,n3
        a(i,j,k)=0.
      enddo
   enddo
enddo

do j=1,n2
   do i=1,n1
      pblht=0.
      do k=2,n3
         pblht=ztn(k,ngrd)*(1.-c(i,j,1)/zmn(nnzp(1)-1,1))
!DSM         if(a(i,j,k).le.tkethrsh) goto 10
      enddo
      10 continue
      do k=1,n3
        a(i,j,k)=pblht
      enddo
   enddo
enddo

return

entry RAMS_comp_etrans(n1,n2,n3,a,b,a2d)
   do j=1,n2
      do i=1,n1
         temp1=a(i,j,1)*b(i,j,1)/cp
         press1=(b(i,j,1)/cp)**cpor*p00
         dens=press1/(rgas*temp1)
         if(i.eq.5.and.j.eq.5) then
            print*,'============++++++'
            print*,temp1,press1,dens,a2d(i,j)
         endif
         a(i,j,1)=a2d(i,j)*dens*1.e-3*39.37*3600.
         do k=2,n3
            a(i,j,k)=a(i,j,1)
         enddo
      enddo
   enddo
return

entry RAMS_comp_slpress(n1,n2,n3,theta,pp,z,slp)
!
!     This subroutine calculates the pressure at level zlev. it
!       is hardwired here to calculate mean sea level pressure,
!       but can be easily changed to calculate pressure at any level
!       by just changing zlev.
!     a standard atmosphere lapse rate of 6.5 c/km is used to
!       interpolate temperature down from level 2 in the model.
!

!c      rlap=-.0065  ! standard temp lapse rate
   rlap=.0025     ! approx standard theta lapse rate
   zlev=0.

   do j=1,n2
      do i=1,n1
         levloop: do k=2,n3
            if(z(i,j,k).ge.zlev) then
               ktop=k
               kbot=k-1
               exit levloop
            endif
         end do levloop

         ddz=zlev-z(i,j,kbot)
         if(zlev.lt.z(i,j,kbot))then
            thbar=(theta(i,j,kbot)-.5*ddz*rlap)
         else
            thbar=.5*(theta(i,j,kbot)+theta(i,j,ktop))
         endif
         slp(i,j)=pp(i,j,kbot)-ddz*sl_g/thbar
         slp(i,j)=(slp(i,j) * cpi)**cpor*p00
      end do
   end do
return

entry RAMS_comp_ctprof(n1,n2,n3,a,b,ngrd)
   do i=1,n1
      do j=1,n2
         kmax=0
         do k=1,n3
            if(a(i,j,k).ge.0.0001.and.b(i,j,k).ge.0.99)kmax=k
         enddo
         if(kmax.gt.2)then
            a(i,j,1)=ztn(kmax,ngrd)
         else
            a(i,j,1)=0.0
         endif
      enddo
   enddo
return

end


subroutine RAMS_comp_multap(n1,n2,n3,n4,a,b)
dimension a(n1,n2,n4),b(n1,n2,n3)
   do k=1,n4
      do j=1,n2
         do i=1,n1
            a(i,j,k)=a(i,j,k)*b(i,j,1)
         enddo
      enddo
   enddo
return
end
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine get_leaf_soil(n1,n2,n4,n5,a,a2)
   implicit none
   integer, intent(in) :: n1,n2,n4,n5
   real, dimension(n1,n2,n4,n5), intent(out) :: a2
   real, dimension(n1,n2,n4*n5), intent(in)  :: a
   integer :: kip, k,i,j,ip
   kip=0
   do ip=1,n5
      do k=1,n4
         kip=kip+1
         do j=1,n2
            do i=1,n1
               a2(i,j,k,ip)=a(i,j,kip)
            end do
         end do
      end do
   end do
   return
end subroutine get_leaf_soil
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine get_cumulus(n1,n2,n3,n6,a,a6)
   implicit none
   integer, intent(in) :: n1,n2,n3,n6
   real, dimension(n1,n2,n3,n6), intent(out) :: a6
   real, dimension(n1,n2,n3*n6), intent(in) :: a
   integer :: kip, k,i,j,ip

   kip=0
   do ip=1,n6
      do k=1,n3
         kip=kip+1
         do j=1,n2
            do i=1,n1
               a6(i,j,k,ip)=a(i,j,kip)
            end do
         end do
      end do
   end do
   return
end subroutine get_cumulus
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine RAMS_comp_richardson(n1,n2,n3,np,rib,z0,speed,thetav_atm,thetav_can,topt,ngrd)

   use somevars
   use rconstants
   use therm_lib, only : qtk, qwtk, dewfrostpoint, rslif, virtt, thetaeiv
   implicit none
   integer                     , intent(in)    :: n1,n2,n3,np
   real   , dimension(n1,n2,np), intent(inout) :: rib
   real   , dimension(n1,n2,np), intent(in)    :: z0,thetav_can
   real   , dimension(n1,n2,n3), intent(in)    :: speed,thetav_atm
   real   , dimension(n1,n2)   , intent(in)    :: topt
   integer                     , intent(in)    :: ngrd
   integer                                     :: i,j,p
   real                                        :: zedtop,zagl,spd

   zedtop = myzmn(mynnzp(1)-1,1)
   do j=1,n2
      do i=1,n1
         zagl=myztn(2,ngrd)*(1.-topt(i,j)/zedtop)
         do p=1,np
            spd = max(speed(i,j,2),0.65)
            rib(i,j,p) = grav * (zagl-z0(i,j,p)) * (thetav_atm(i,j,2)-thetav_can(i,j,p))   &
                       / (0.5 * (thetav_can(i,j,2)+thetav_can(i,j,p)) * spd**2)
         end do
      end do
   end do
   return
end subroutine RAMS_comp_richardson
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine RAMS_comp_dn0(n1,n2,n3,a,b,c,topt,ngrd)

   use somevars
   use rconstants
   use therm_lib, only : qtk, qwtk, dewfrostpoint, rslif, virtt, thetaeiv
   implicit none
   integer                     , intent(in)    :: n1,n2,n3
   real   , dimension(n1,n2,n3), intent(inout) :: a,b,c
   real   , dimension(n1,n2)   , intent(in)    :: topt
   integer                     , intent(in)    :: ngrd
   integer                                     :: i,j,k
   real   , dimension(n3)                      :: scratch2,scratch11,scratch12
   real                                        :: zedtop,c1,c2,c3

   zedtop = myzmn(mynnzp(1)-1,1)
   do j=1,n2
      do i=1,n1
         do k=1,n3
           scratch2(k)=myztn(k,ngrd)*(1.-topt(i,j)/zedtop)+topt(i,j)
         enddo
         call htint(n3,mypi01dn(1,ngrd),myztn(1,ngrd),n3,scratch11,scratch2)
         call htint(n3,myth01dn(1,ngrd),myztn(1,ngrd),n3,scratch12,scratch2)       
         do k=1,n3
            b(i,j,k)=scratch12(k)
         enddo
         a(i,j,n3) = scratch11(n3)

         c1=grav*2.*(1.-topt(i,j)/zedtop)
         c2=(1-cpor)
         c3=cp**c2
         do k=n3-1,1,-1
            a(i,j,k)=a(i,j,k+1) +c1/((b(i,j,k)+b(i,j,k+1))*mydzmn(k,ngrd))
         enddo
         do k=1,n3
            c(i,j,k)=(c3*p00)/(rdry*b(i,j,k)*a(i,j,k)**c2)
         enddo

      enddo
   enddo
   return
end subroutine RAMS_comp_dn0
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine RAMS_comp_relvortx(n1,n2,n3,a,b,c,topt,ngrd)

   use somevars
   use rconstants
   use therm_lib, only : qtk, qwtk, dewfrostpoint, rslif, virtt, thetaeiv
   implicit none
   integer                     , intent(in)    :: n1,n2,n3
   real   , dimension(n1,n2,n3), intent(inout) :: a,b,c
   real   , dimension(n1,n2)   , intent(in)    :: topt
   integer                     , intent(in)    :: ngrd
   integer                                     :: i,j,k,j1,j2,k1,k2
   real                                        :: factor
   real   , dimension(n1+n2+n3)                :: dum1,dum2

   factor = myztn(1,ngrd) / myztn(2,ngrd)
   do j=1,n2
      do i=1,n1
         a(i,j,1) = a(i,j,2) * factor
      enddo
   enddo

   call gradr(n1,n2,n3,2,n1-1,1,n2-1,b,c,'ydir','wpnt',topt  &
      ,myxmn(:,ngrd),myxtn(:,ngrd),myymn(:,ngrd),myytn(:,ngrd)  &
      ,myzmn(:,ngrd),myztn(:,ngrd),mydeltayn(ngrd)  &
      ,mydzmn(:,ngrd),mydztn(:,ngrd),dum1,dum2,myzmn(mynnzp(1)-1,1)  &
      ,myjdim,myihtran)
   call gradr(n1,n2,n3,2,n1-1,1,n2-1,a,b,'zdir','vpnt',topt  &
      ,myxmn(:,ngrd),myxtn(:,ngrd),myymn(:,ngrd),myytn(:,ngrd)  &
      ,myzmn(:,ngrd),myztn(:,ngrd),mydeltayn(ngrd)  &
      ,mydzmn(:,ngrd),mydztn(:,ngrd),dum1,dum2,myzmn(mynnzp(1)-1,1)  &
      ,myjdim,myihtran)

   do k=1,n3
      do j=2,n2-1
         do i=1,n1
            b(i,j,k) = c(i,j,k) - b(i,j,k)
         enddo
      enddo
   enddo
   do j = 1,n2
      do i =1,n1
         do k =1,n3
            a(i,j,k) = 0.
         end do
      end do
   end do
   do j = 2,n2-1
      j1 = max(j-1,2)
      j2 = min(j,n2-2)
      do i = 2,n1-1
         do k = 1,n3
            k1 = max(k-1,1)
            k2 = min(k,n3-1)
            a(i,j,k) =0.25 * (b(i,j1,k1) + b(i,j1,k2)  &
                            + b(i,j2,k1) + b(i,j2,k2))
         enddo
      enddo
   enddo
   
  return
end subroutine RAMS_comp_relvortx
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine RAMS_comp_relvorty(n1,n2,n3,a,b,c,topt,ngrd)

   use somevars
   use rconstants
   use therm_lib, only : qtk, qwtk, dewfrostpoint, rslif, virtt, thetaeiv
   implicit none
   integer                     , intent(in)    :: n1,n2,n3
   real   , dimension(n1,n2,n3), intent(inout) :: a,b,c
   real   , dimension(n1,n2)   , intent(in)    :: topt
   integer                     , intent(in)    :: ngrd
   integer                                     :: i,j,k,i1,i2,k1,k2
   real                                        :: factor
   real   , dimension(n1+n2+n3)                :: dum1,dum2

   factor = myztn(2,ngrd) / myztn(3,ngrd)
   do j=1,n2
      do i=1,n1
         a(i,j,1) = a(i,j,2) * factor
      enddo
   enddo

   call gradr(n1,n2,n3,1,n1-1,2,n2-1,b,c,'xdir','wpnt',topt  &
      ,myxmn(:,ngrd),myxtn(:,ngrd),myymn(:,ngrd),myytn(:,ngrd)  &
      ,myzmn(:,ngrd),myztn(:,ngrd),mydeltayn(ngrd)  &
      ,mydzmn(:,ngrd),mydztn(:,ngrd),dum1,dum2,myzmn(mynnzp(1)-1,1)  &
      ,myjdim,myihtran)
   call gradr(n1,n2,n3,1,n1-1,2,n2-1,a,b,'zdir','upnt',topt  &
      ,myxmn(:,ngrd),myxtn(:,ngrd),myymn(:,ngrd),myytn(:,ngrd)  &
      ,myzmn(:,ngrd),myztn(:,ngrd),mydeltayn(ngrd)  &
      ,mydzmn(:,ngrd),mydztn(:,ngrd),dum1,dum2,myzmn(mynnzp(1)-1,1)  &
      ,myjdim,myihtran)

   do k=1,n3
      do j=1,n2
         do i=1,n1
            b(i,j,k) = b(i,j,k) - c(i,j,k)
         enddo
      enddo
   enddo

   do j = 1,n2
      do i =1,n1
         do k =1,n3
            a(i,j,k) = 0.
         end do
      end do
   end do

   do j = 2,n2-1
      do i = 2,n1-1
         i1 = max(i-1,2)
         i2 = min(i,n1-2)
         do k = 1,n3
            k1 = max(k-1,1)
            k2 = min(k,n3-1)
            a(i,j,k) = 0.25 * (b(i1,j,k1) + b(i1,j,k2)  &
                             + b(i2,j,k1) + b(i2,j,k2))
         enddo
      enddo
   enddo


   return
end subroutine RAMS_comp_relvorty
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine RAMS_comp_relvortz(n1,n2,n3,a,b,c,topt,ngrd)

   use somevars
   use rconstants
   use therm_lib, only : qtk, qwtk, dewfrostpoint, rslif, virtt, thetaeiv
   implicit none
   integer                     , intent(in)    :: n1,n2,n3
   real   , dimension(n1,n2,n3), intent(inout) :: a,b,c
   real   , dimension(n1,n2)   , intent(in)    :: topt
   integer                     , intent(in)    :: ngrd
   integer                                     :: i,j,k,i1,i2,j1,j2
   real                                        :: factor
   real   , dimension(n1+n2+n3)                :: dum1,dum2

   factor = myztn(1,ngrd) / myztn(2,ngrd)
   do j=1,n2
      do i=1,n1
         a(i,j,1) = a(i,j,2) * factor
         b(i,j,1) = b(i,j,2) * factor
      enddo
   enddo

   call gradr(n1,n2,n3,1,n1-1,1,n2-1,b,c,'xdir','vpnt',topt  &
      ,myxmn(:,ngrd),myxtn(:,ngrd),myymn(:,ngrd),myytn(:,ngrd)  &
      ,myzmn(:,ngrd),myztn(:,ngrd),mydeltayn(ngrd)  &
      ,mydzmn(:,ngrd),mydztn(:,ngrd),dum1,dum2,myzmn(mynnzp(1)-1,1)  &
      ,myjdim,myihtran)

   call gradr(n1,n2,n3,1,n1-1,1,n2-1,a,b,'ydir','upnt',topt  &
      ,myxmn(:,ngrd),myxtn(:,ngrd),myymn(:,ngrd),myytn(:,ngrd)  &
      ,myzmn(:,ngrd),myztn(:,ngrd),mydeltayn(ngrd)  &
      ,mydzmn(:,ngrd),mydztn(:,ngrd),dum1,dum2,myzmn(mynnzp(1)-1,1)  &
      ,myjdim,myihtran)

   do k=1,n3
      do j=1,n2
         do i=1,n1
            b(i,j,k) = c(i,j,k) - b(i,j,k)
         enddo
      enddo
   enddo

   do j = 1,n2
      j1 = max(j-1,1)
      j2 = min(j,n2-1)
      do i = 1,n1
         i1 = max(i-1,1)
         i2 = min(i,n1-1)
         do k = 1,n3
            a(i,j,k) = 0.25 * (b(i1,j1,k) + b(i1,j2,k)  &
                             + b(i2,j1,k) + b(i2,j2,k))
         enddo
      enddo
   enddo
   return
end subroutine RAMS_comp_relvortz
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine RAMS_comp_totvortz(n1,n2,n3,a,b,c,topt,ngrd)

   use somevars
   use rconstants
   use therm_lib, only : qtk, qwtk, dewfrostpoint, rslif, virtt, thetaeiv
   implicit none
   integer                     , intent(in)    :: n1,n2,n3
   real   , dimension(n1,n2,n3), intent(inout) :: a,b,c
   real   , dimension(n1,n2)   , intent(in)    :: topt
   integer                     , intent(in)    :: ngrd
   integer                                     :: i,j,k,i1,i2,j1,j2
   real                                        :: factor,omega2,fcor,xlon,xlat
   real   , dimension(n1+n2+n3)                :: dum1,dum2

   factor = myztn(1,ngrd) / myztn(2,ngrd)
   do j=1,n2
      do i=1,n1
         a(i,j,1) = a(i,j,2) * factor
         b(i,j,1) = b(i,j,2) * factor
      enddo
   enddo

   call gradr(n1,n2,n3,1,n1-1,1,n2-1,b,c,'xdir','vpnt',topt  &
      ,myxmn(:,ngrd),myxtn(:,ngrd),myymn(:,ngrd),myytn(:,ngrd)  &
      ,myzmn(:,ngrd),myztn(:,ngrd),mydeltayn(ngrd)  &
      ,mydzmn(:,ngrd),mydztn(:,ngrd),dum1,dum2,myzmn(mynnzp(1)-1,1)  &
      ,myjdim,myihtran)
   call gradr(n1,n2,n3,1,n1-1,1,n2-1,a,b,'ydir','upnt',topt  &
      ,myxmn(:,ngrd),myxtn(:,ngrd),myymn(:,ngrd),myytn(:,ngrd)  &
      ,myzmn(:,ngrd),myztn(:,ngrd),mydeltayn(ngrd)  &
      ,mydzmn(:,ngrd),mydztn(:,ngrd),dum1,dum2,myzmn(mynnzp(1)-1,1)  &
      ,myjdim,myihtran)

   do k=1,n3
      do j=1,n2
         do i=1,n1
            b(i,j,k) = c(i,j,k) - b(i,j,k)
         enddo
      enddo
   enddo

   do j = 1,n2
      j1 = max(j-1,1)
      j2 = min(j,n2-1)
      do i = 1,n1
         i1 = max(i-1,1)
         i2 = min(i,n1-1)
         do k = 1,n3
            a(i,j,k) = 0.25 * (b(i1,j1,k) + b(i1,j2,k)  &
                             + b(i2,j1,k) + b(i2,j2,k))
         enddo
      enddo
   enddo

   omega2 = 2. * omega
   do j = 1,n2
      do i = 1,n1
         call xy_ll(xlat,xlon,myplatn(ngrd),myplonn(ngrd),myxtn(i,ngrd),myytn(j,ngrd))
         fcor = omega2 * sin(xlat * pio180)
         do k = 1,n3
            a(i,j,k) = a(i,j,k) + fcor
         enddo
      enddo
   enddo
   return
end subroutine RAMS_comp_totvortz
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine RAMS_comp_potvortz(n1,n2,n3,a,b,c,e,topt,ngrd)
   use somevars
   use rconstants
   use therm_lib, only : qtk, qwtk, dewfrostpoint, rslif, virtt, thetaeiv
   implicit none
   integer                     , intent(in)    :: n1,n2,n3
   real   , dimension(n1,n2,n3), intent(inout) :: a,b,c,e
   real   , dimension(n1,n2)   , intent(in)    :: topt
   integer                     , intent(in)    :: ngrd
   integer                                     :: i,j,k
   real   , dimension(n1+n2+n3)                :: dum1,dum2

   call gradr(n1,n2,n3,1,n1-1,1,n2-1,b,e,'zdir','tpnt',topt  &
      ,myxmn(:,ngrd),myxtn(:,ngrd),myymn(:,ngrd),myytn(:,ngrd)  &
      ,myzmn(:,ngrd),myztn(:,ngrd),mydeltayn(ngrd)  &
      ,mydzmn(:,ngrd),mydztn(:,ngrd),dum1,dum2,myzmn(mynnzp(1)-1,1)  &
      ,myjdim,myihtran)

   do k=1,n3
      do j=1,n2
         do i=1,n1
            a(i,j,k) = a(i,j,k) * e(i,j,k) / (grav * c(i,j,k))
         enddo
      enddo
   enddo
   return
end subroutine RAMS_comp_potvortz
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine RAMS_comp_solenoidx(nx,ny,nz,alpha,press,solex,topt,ngrd)
   use somevars
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                     , intent(in)    :: nx
   integer                     , intent(in)    :: ny
   integer                     , intent(in)    :: nz
   real   , dimension(nx,ny,nz), intent(inout) :: alpha ! Specific Volume
   real   , dimension(nx,ny,nz), intent(inout) :: press
   real   , dimension(nx,ny,nz), intent(inout) :: solex
   real   , dimension(nx,ny)   , intent(in)    :: topt
   integer                     , intent(in)    :: ngrd
   !----- Local variables. ----------------------------------------------------------------!
   real   , dimension(nx,ny,nz)                :: dady
   real   , dimension(nx,ny,nz)                :: dadz
   real   , dimension(nx,ny,nz)                :: dpdy
   real   , dimension(nx,ny,nz)                :: dpdz
   integer                                     :: x
   integer                                     :: y
   integer                                     :: z
   real   , dimension(nx+ny+nz)                :: dum1
   real   , dimension(nx+ny+nz)                :: dum2
   !---------------------------------------------------------------------------------------!
   do z=1,nz
      do y=1,ny
         do x=1,nx
            solex(x,y,z) = 0.
         end do
      end do
   end do

   call gradr(nx,ny,nz,2,nx-1,2,ny-1,alpha,dadz,'zdir','tpnt',topt,myxmn(:,ngrd)           &
             ,myxtn(:,ngrd),myymn(:,ngrd),myytn(:,ngrd),myzmn(:,ngrd),myztn(:,ngrd)        &
             ,mydeltayn(ngrd),mydzmn(:,ngrd),mydztn(:,ngrd),dum1,dum2,myzmn(mynnzp(1)-1,1) &
             ,myjdim,myihtran)
   call gradr(nx,ny,nz,2,nx-1,2,ny-1,press,dpdy,'ydir','tpnt',topt,myxmn(:,ngrd)           &
             ,myxtn(:,ngrd),myymn(:,ngrd),myytn(:,ngrd),myzmn(:,ngrd),myztn(:,ngrd)        &
             ,mydeltayn(ngrd),mydzmn(:,ngrd),mydztn(:,ngrd),dum1,dum2,myzmn(mynnzp(1)-1,1) &
             ,myjdim,myihtran)

   do z=1,nz-1
      do y=2,ny-1
         do x=2,nx-1
            solex(x,y,z) = dadz(x,y,z) * dpdy(x,y,z)
         end do
      end do
   end do


   call gradr(nx,ny,nz,2,nx-1,2,ny-1,press,dpdz,'zdir','tpnt',topt,myxmn(:,ngrd)           &
             ,myxtn(:,ngrd),myymn(:,ngrd),myytn(:,ngrd),myzmn(:,ngrd),myztn(:,ngrd)        &
             ,mydeltayn(ngrd),mydzmn(:,ngrd),mydztn(:,ngrd),dum1,dum2,myzmn(mynnzp(1)-1,1) &
             ,myjdim,myihtran)
   call gradr(nx,ny,nz,2,nx-1,2,ny-1,alpha,dady,'ydir','tpnt',topt,myxmn(:,ngrd)           &
             ,myxtn(:,ngrd),myymn(:,ngrd),myytn(:,ngrd),myzmn(:,ngrd),myztn(:,ngrd)        &
             ,mydeltayn(ngrd),mydzmn(:,ngrd),mydztn(:,ngrd),dum1,dum2,myzmn(mynnzp(1)-1,1) &
             ,myjdim,myihtran)


   do z=1,nz-1
      do y=2,ny-1
         do x=2,nx-1
            solex(x,y,z) = solex(x,y,z) - dpdz(x,y,z) * dady(x,y,z)
         end do
      end do
   end do

  return
end subroutine RAMS_comp_solenoidx
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine RAMS_comp_solenoidy(nx,ny,nz,alpha,press,soley,topt,ngrd)
   use somevars
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                     , intent(in)    :: nx
   integer                     , intent(in)    :: ny
   integer                     , intent(in)    :: nz
   real   , dimension(nx,ny,nz), intent(inout) :: alpha ! Specific Volume
   real   , dimension(nx,ny,nz), intent(inout) :: press
   real   , dimension(nx,ny,nz), intent(inout) :: soley
   real   , dimension(nx,ny)   , intent(in)    :: topt
   integer                     , intent(in)    :: ngrd
   !----- Local variables. ----------------------------------------------------------------!
   real   , dimension(nx,ny,nz)                :: dadx
   real   , dimension(nx,ny,nz)                :: dadz
   real   , dimension(nx,ny,nz)                :: dpdx
   real   , dimension(nx,ny,nz)                :: dpdz
   integer                                     :: x
   integer                                     :: y
   integer                                     :: z
   real   , dimension(nx+ny+nz)                :: dum1
   real   , dimension(nx+ny+nz)                :: dum2
   !----- Local constants. ----------------------------------------------------------------!
   logical                     , parameter     :: printdbg = .false.
   !---------------------------------------------------------------------------------------!
   do z=1,nz
      do y=1,ny
         do x=1,nx
            soley(x,y,z) = 0.
         end do
      end do
   end do

   call gradr(nx,ny,nz,2,nx-1,2,ny-1,alpha,dadx,'xdir','tpnt',topt,myxmn(:,ngrd)           &
             ,myxtn(:,ngrd),myymn(:,ngrd),myytn(:,ngrd),myzmn(:,ngrd),myztn(:,ngrd)        &
             ,mydeltayn(ngrd),mydzmn(:,ngrd),mydztn(:,ngrd),dum1,dum2,myzmn(mynnzp(1)-1,1) &
             ,myjdim,myihtran)
   call gradr(nx,ny,nz,2,nx-1,2,ny-1,press,dpdz,'zdir','tpnt',topt,myxmn(:,ngrd)           &
             ,myxtn(:,ngrd),myymn(:,ngrd),myytn(:,ngrd),myzmn(:,ngrd),myztn(:,ngrd)        &
             ,mydeltayn(ngrd),mydzmn(:,ngrd),mydztn(:,ngrd),dum1,dum2,myzmn(mynnzp(1)-1,1)  &
              ,myjdim,myihtran)

   do z=1,nz-1
      do y=2,ny-1
         do x=2,nx-1
            soley(x,y,z) = dadx(x,y,z) * dpdz(x,y,z)
         end do
      end do
   end do

   call gradr(nx,ny,nz,2,nx-1,2,ny-1,press,dpdx,'xdir','tpnt',topt,myxmn(:,ngrd)           &
             ,myxtn(:,ngrd),myymn(:,ngrd),myytn(:,ngrd),myzmn(:,ngrd),myztn(:,ngrd)        &
             ,mydeltayn(ngrd),mydzmn(:,ngrd),mydztn(:,ngrd),dum1,dum2,myzmn(mynnzp(1)-1,1) &
             ,myjdim,myihtran)
   call gradr(nx,ny,nz,2,nx-1,2,ny-1,alpha,dadz,'zdir','tpnt',topt,myxmn(:,ngrd)           &
             ,myxtn(:,ngrd),myymn(:,ngrd),myytn(:,ngrd),myzmn(:,ngrd),myztn(:,ngrd)        &
             ,mydeltayn(ngrd),mydzmn(:,ngrd),mydztn(:,ngrd),dum1,dum2,myzmn(mynnzp(1)-1,1) &
             ,myjdim,myihtran)

   do z=1,nz-1
      do y=2,ny-1
         do x=2,nx-1
            soley(x,y,z) = soley(x,y,z) - dadz(x,y,z) * dpdx(x,y,z)
         end do
      end do
   end do
   
   if (printdbg) then
      write (unit=61,fmt='(3(a5,1x),7(a12,1x))') '    X','    Y','    Z','       ALPHA'    &
                                                ,'       PRESS','        DADX'             &
                                                ,'        DADZ','        DPDX'             &
                                                ,'        DPDZ','       SOLEY'

      do z=1,nz-1
         do y=2,ny-1
            do x=2,nx-1
               write(unit=61,fmt='(3(i5,1x),7(es12.5,1x))') x,y,z,alpha(x,y,z)             &
                                                                 ,press(x,y,z)             &
                                                                 , dadx(x,y,z)             &
                                                                 , dadz(x,y,z)             &
                                                                 , dpdx(x,y,z)             &
                                                                 , dpdz(x,y,z)             &
                                                                 ,soley(x,y,z)
            end do
         end do
      end do
   end if

   return
end subroutine RAMS_comp_solenoidy
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine RAMS_comp_sfcwmeantemp(n1,n2,ns,np,a,b,c,d,e)
   use rconstants
   use therm_lib, only: qtk
   implicit none 
   integer :: n1,n2,ns,np,nlev,i,j,ip,k
   real, dimension(n1,n2,np)    :: a,d,e
   real, dimension(n1,n2,ns,np) :: b,c
   real :: temptemp,fracliq,snowarea,xmasstot
   real, parameter :: undef=-9.99e33
   !a  area
   !b energy
   !c  mass
   !d  nlev
   !e integrated value
   do j=1,n2
     do i=1,n1
        e(i,j,1) = 0.
        if (a(i,j,1) <= 0.99) then
           do ip=2,np
              xmasstot = 0.
              nlev=nint(d(i,j,ip))
              e(i,j,ip) = 0.
              if (nlev > 0 ) then
                 do k=1,nlev
                    if (c(i,j,k,ip) > 1.e-6) then
                       xmasstot  = xmasstot + c(i,j,k,ip)
                       call qtk(b(i,j,k,ip),temptemp,fracliq)
                       e(i,j,ip) = e(i,j,ip) + temptemp*c(i,j,k,ip)
                    end if
                 end do
                 if (xmasstot > 1.e-6) e(i,j,ip) = e(i,j,ip) / xmasstot - t00
              else
                 e(i,j,ip) = 0.
              end if
           end do

           snowarea = 0.
           do ip=2,np
              if( nint(d(i,j,ip)) > 0) then
                 snowarea = snowarea + a(i,j,ip)
                 e(i,j,1) = e(i,j,1) + e(i,j,ip)*a(i,j,ip)
              end if
           end do
           if (snowarea > 0.) then
              e(i,j,1) = e(i,j,1) / snowarea
           else
              e(i,j,1) = undef
           end if
        else
           e(i,j,1) = undef
        end if
     end do
   end do
   return
end subroutine RAMS_comp_sfcwmeantemp
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine RAMS_comp_thetaeiv(n1,n2,n3,xxx,temp,pres,rv,rtp)
   use therm_lib, only: thetaeiv
   implicit none
   integer, intent(in) :: n1, n2, n3
   real, dimension(n1,n2,n3), intent(inout) :: xxx,temp,pres,rv,rtp
   integer :: i,j,k
   !----- xxx comes as thil, gets out as theta_e_iv. --------------------------------------!
   do k=1,n3
      do j=1,n2
         do i=1,n1
            if (rtp(i,j,k) < rv(i,j,k)) rtp(i,j,k) = rv(i,j,k)
            xxx(i,j,k)=thetaeiv(xxx(i,j,k),pres(i,j,k),temp(i,j,k),rv(i,j,k),rtp(i,j,k)    &
                               ,12,.true.)
         end do
      end do
   end do
   return
end subroutine RAMS_comp_thetaeiv
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine RAMS_comp_sfcwinteg(n1,n2,ns,np,a,c,d,e)
   implicit none 
   integer :: n1,n2,ns,np,nlev,i,j,ip,k
   real, dimension(n1,n2,np)    :: a,d,e
   real, dimension(n1,n2,ns,np) :: c
   real, parameter :: undef=-9.99e33
   !a  area                                                   
   !c mass/depth                                             
   !d  nlev                                                   
   !e  integrated value                                       
      do j=1,n2                                               
        do i=1,n1                                             
           if (a(i,j,1) > 0.99) then                          
              e(i,j,1) = -9.99e33                                
           else
              e(i,j,1) = 0.                                  
              do ip=2,np                                     
                 nlev=nint(d(i,j,ip))
                 e(i,j,ip) = 0.                               
                 do k=1,nlev                                  
                    e(i,j,ip) = e(i,j,ip) + c(i,j,k,ip)      
                 end do                                       
                 e(i,j,1) = e(i,j,1) + e(i,j,ip)*a(i,j,ip)    
              end do                                          
              e(i,j,1) = e(i,j,1) / (1.-a(i,j,1))            
           end if                                             
        end do                                                
      end do                                                  
   return        
end subroutine RAMS_comp_sfcwinteg
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This routine is for quantities defined for all patches.                               !
!------------------------------------------------------------------------------------------!
subroutine RAMS_comp_patchsum(nx,ny,nz,np,iovar)
   use leaf_coms, only : min_patch_area
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                         , intent(in)    :: nx
   integer                         , intent(in)    :: ny
   integer                         , intent(in)    :: nz
   integer                         , intent(in)    :: np
   real   , dimension(*)           , intent(inout) :: iovar
   !----- Local variables. ----------------------------------------------------------------!
   integer                                         :: x
   integer                                         :: y
   integer                                         :: z
   integer                                         :: p
   integer                                         :: iareaa
   integer                                         :: iareaz
   integer                                         :: ivala
   integer                                         :: ivalz
   integer                                         :: iouta
   integer                                         :: ioutz
   integer                                         :: narea
   integer                                         :: nvals
   integer                                         :: nout
   real                                            :: totarea
   real   , dimension(nx,ny,nz,np)                 :: patval
   real   , dimension(nx,ny,nz)                    :: psum
   real   , dimension(nx,ny,np)                    :: pfarea
   !---------------------------------------------------------------------------------------!

   !----- Define the indices for copying in and out. --------------------------------------!
   narea  = nx * ny * np
   nvals  = nx * ny * nz * np
   nout   = nx * ny * nz
   iareaa = 1
   iareaz = narea
   ivala  = iareaz + 1
   ivalz  = iareaz + nvals
   iouta  = 1
   ioutz  = nout

   !----- Copy the patch fraction area to a scratch array. --------------------------------!
   call atob(narea,iovar(iareaa:iareaz),pfarea)
   
   !----- Copy the patch-structured data to a scratch array. ------------------------------!
   call atob(nvals,iovar(ivala:ivalz),patval)

   !----- Compute the patch weigthed average, including water. ----------------------------!
   do z = 1,nz
      do y = 1,ny
         do x = 1,nx
            psum(x,y,z) = 0.
            totarea     = 0.
            ploop: do p = 1,np
               if (pfarea(x,y,p) < min_patch_area) cycle ploop
               
               psum(x,y,z) = psum(x,y,z) + pfarea(x,y,p) * patval(x,y,z,p)
               totarea     = totarea + pfarea(x,y,p)
            end do ploop
            psum(x,y,z) = psum(x,y,z) / totarea
         end do
      end do
   end do

   !----- Copy psum into iovar, which will be used for output. ----------------------------!
   call atob(nout,psum,iovar(iouta:ioutz))

   return
end subroutine RAMS_comp_patchsum
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!      This routine is for quantities that are not defined for water patches.              !
!------------------------------------------------------------------------------------------!
subroutine RAMS_comp_patchsum_l(nx,ny,nz,np,iovar)
   use leaf_coms, only : min_patch_area
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                         , intent(in)    :: nx
   integer                         , intent(in)    :: ny
   integer                         , intent(in)    :: nz
   integer                         , intent(in)    :: np
   real   , dimension(*)           , intent(inout) :: iovar
   !----- Local variables. ----------------------------------------------------------------!
   integer                                         :: x
   integer                                         :: y
   integer                                         :: z
   integer                                         :: p
   integer                                         :: iareaa
   integer                                         :: iareaz
   integer                                         :: ivala
   integer                                         :: ivalz
   integer                                         :: iouta
   integer                                         :: ioutz
   integer                                         :: narea
   integer                                         :: nvals
   integer                                         :: nout
   real                                            :: landarea
   real   , dimension(nx,ny,nz,np)                 :: patval
   real   , dimension(nx,ny,nz)                    :: psum
   real   , dimension(nx,ny,np)                    :: pfarea
   !---------------------------------------------------------------------------------------!


   !----- Define the indices for copying in and out. --------------------------------------!
   narea  = nx * ny * np
   nvals  = nx * ny * nz * np
   nout   = nx * ny * nz
   iareaa = 1
   iareaz = narea
   ivala  = iareaz + 1
   ivalz  = iareaz + nvals
   iouta  = 1
   ioutz  = nout

   !----- Copy the patch fraction area to a scratch array. --------------------------------!
   call atob(narea,iovar(iareaa:iareaz),pfarea)
   
   !----- Copy the patch-structured data to a scratch array. ------------------------------!
   call atob(nvals,iovar(ivala:ivalz),patval)

   !----- Compute the patch weigthed average, excluding water. ----------------------------!
   do z = 1,nz
      do y = 1,ny
         do x = 1,nx
            if (pfarea(x,y,1) < 1.-min_patch_area) then
               psum(x,y,z) = 0.
               landarea    = 0.
               do p = 2,np
                  psum(x,y,z) = psum(x,y,z) + pfarea(x,y,p) * patval(x,y,z,p)
                  landarea    = landarea    + pfarea(x,y,p)
               end do
               psum(x,y,z)    = psum(x,y,z) / landarea
            else
               psum(x,y,z)    = -9.99e33
            end if
         end do
      end do
   end do

   !----- Copy psum into iovar, which will be used for output. ----------------------------!
   call atob(nout,psum,iovar(iouta:ioutz))

   return
end subroutine RAMS_comp_patchsum_l
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!      Extract value from largest patch.                                                   !
!------------------------------------------------------------------------------------------!
subroutine RAMS_comp_bigpatch(n1,n2,n3,n4,a,f,bpat)
   use somevars, only : mynbig
   implicit none
   integer, intent(in) :: n1,n2,n3,n4
   real   , dimension(n1,n2,n3,n4) , intent(in)    :: a
   real   , dimension(n1,n2,mynbig), intent(inout) :: f
   real   , dimension(n1,n2,n3)    , intent(out)   :: bpat
   integer                                         :: i,j,k,ip
   do k = 1,n3
      do j = 1,n2
         do i = 1,n1
            if (f(i,j,2) >= f(i,j,1)) then
               bpat(i,j,k) = a(i,j,k,2)
            else
               bpat(i,j,k) = a(i,j,k,1)
            end if
         end do
      end do
   end do

   !---------------------------------------------------------------------------------------!
   !      Copy bpat into f, which was passed in as a(1).                                   !
   !---------------------------------------------------------------------------------------!

   do k = 1,n3
      do j = 1,n2
         do i = 1,n1
            f(i,j,k) = bpat(i,j,k)
         end do
      end do
   end do
   return
end subroutine RAMS_comp_bigpatch
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine RAMS_comp_tvegc(n1,n2,n3,a,b,c,e)
   use rconstants, only : t00
   use therm_lib , only : qwtk
   implicit none
   integer, intent(in) :: n1,n2,n3
   real, dimension(n1,n2,n3), intent(in)    :: c,e
   real, dimension(n1,n2,n3), intent(inout) :: a,b
   integer :: i,j,k
   real :: temptemp,fracliq
!MLO: Input:  a = Veg. Energy        in     J/m2
!             b = Veg. Water         in    kg/m2
!             c = Veg. Heat capacity in   J/m2/K
!             e = Canopy Theta       in  Celsius
!     Output: a = Temperature in Celsius
!             b = Fraction in liquid phase

   do k=1,n3
      do j=1,n2
         do i=1,n1
            !------------------------------------------------------------------------ 
            !   Compute tveg only if there is enough heat capacity, otherwise assign
            ! canopy temperature.
            !------------------------------------------------------------------------ 
            if (c(i,j,k) > 10.) then
               call qwtk(a(i,j,k),b(i,j,k),c(i,j,k),temptemp,fracliq)
               a(i,j,k) = temptemp-t00
               b(i,j,k) = fracliq
            else
               a(i,j,k) = e(i,j,k)
               b(i,j,k) = 0.5
            end if
         end do
      end do
   end do
   return
end subroutine RAMS_comp_tvegc

subroutine RAMS_comp_5050(n1,n2,n3,a,d)
real a(n1,n2),d(n1,n2,n3)

do j = 1,n2
   do i = 1,n1
      a(i,j) = .5 * (a(i,j) + d(i,j,2))
   enddo
enddo

return
end subroutine RAMS_comp_5050
!------------------subroutine to calculate cloud fraction
subroutine cldfraction(n1,n2,n3,frac,pi,rh)
use rconstants
implicit none

integer :: i,j,k,kmax,n1,n2,n3
real :: frac(n1,n2),pi(n1,n2,n3),rh(n1,n2,n3)
real, allocatable::rhc(:), cs(:)
real :: kappai,c_1,c_2,c_junk,pop2,csmax

c_1     = 2.
c_junk  = 3.
c_2     = c_junk**0.5
kappai = (1./.286)


allocate (rhc(n3),cs(n3) )
      print*,'+++++++:',n1,n2,n3

do j=1,n2
   do i=1,n1
      frac(i,j) = 0.
      csmax = 0.
      kmax  = 0
      do k = 1, n3
         rhc(k)= 0.
         cs(k)= 0.
      enddo

      do k = 1, n3
         pop2 = (pi(i,j,k)/pi(i,j,2))**kappai

         rhc(k) = 100. - (100.*c_1*pop2)*  &
               (1.-pop2)*(1.+c_2*(pop2-0.5))

         if(rh(i,j,k) .ge. rhc(k))then
            if(rhc(k).eq.100.)rhc(k)=rhc(k)+0.0000001
            cs(k) = ( (rh(i,j,k)-rhc(k))/(100.-rhc(k)) ) **2. 
         else
            cs(k) = 0.
         endif
         if(cs(k).gt.csmax)then
            csmax=cs(k)
            kmax = k
         endif
         frac(i,j) = frac(i,j) + cs(k)*(1./float(k))
      if(i==20.and.j==20) print*,'+++++++:',k,pi(i,j,k),rh(i,j,k),frac(i,j)
      enddo
      
      csmax=max(csmax,0.)

!      frac(i,j) = 1.-min(1.,max(0.,csmax))
      frac(i,j) = 1.-min(1.,max(0.,frac(i,j)))  ! actually returns 
                                                ! clear sky fraction

   enddo
enddo

deallocate (rhc,cs)

return
end
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine will compute the mixing ratio (water vapour or CO2 or any other      !
! tracer), using the similarity theory.  It currently uses either the Louis (1979), the    !
! Oncley and Dudhia (1995) or the Beljaars-Holtslag (1991) methods, and the choice is done !
! based on the ISTAR used in the run.                                                      !
!                                                                                          !
! 1. Based on L79;                                                                         !
! 2. Based on: OD95, but with some terms computed as in L79 and B71 to avoid singular-     !
!    ities.                                                                                !
! 3. Based on BH91, using an iterative method to find zeta, and using the modified         !
!    equation for stable layers.                                                           !
!                                                                                          !
! References:                                                                              !
! B71.  BUSINGER, J.A, et. al; Flux-Profile relationships in the atmospheric surface       !
!           layer. J. Atmos. Sci., 28, 181-189, 1971.                                      !
! L79.  LOUIS, J.F.; Parametric Model of vertical eddy fluxes in the atmosphere.           !
!           Boundary-Layer Meteor., 17, 187-202, 1979.                                     !
! BH91. BELJAARS, A.C.M.; HOLTSLAG, A.A.M.; Flux parameterization over land surfaces for   !
!           atmospheric models. J. Appl. Meteor., 30, 327-341, 1991.                       !
! OD95. ONCLEY, S.P.; DUDHIA, J.; Evaluation of surface fluxes from MM5 using observa-     !
!           tions.  Mon. Wea. Rev., 123, 3344-3357, 1995.                                  !
!------------------------------------------------------------------------------------------!
subroutine RAMS_reduced_prop(nx,ny,nz,np,ng,which,topt,theta_atm,rvap_atm,co2_atm,uspd_atm &
                            ,theta_can,rvap_can,co2_can,prss_can,zout,rough,rib,zeta,parea &
                            ,ustar,tstar,rstar,cstar,varred)
   use rpost_coms, only : isfcl          ! ! intent(in)
   use somevars  , only : myztn          & ! intent(in)
                        , myzmn          & ! intent(in)
                        , mynnzp         & ! intent(in)
                        , myistar        ! ! intent(in)
   use rconstants, only : grav           & ! intent(in)
                        , p00i           & ! intent(in)
                        , rocp           & ! intent(in)
                        , vonk           & ! intent(in)
                        , ep             & ! intent(in)
                        , toodry         ! ! intent(in)
   use therm_lib , only : virtt          & ! function
                        , eslif          & ! function
                        , tslif          ! ! function
   use leaf_coms , only : ustmin         & ! intent(in)
                        , ubmin          & ! intent(in)
                        , bl79           & ! intent(in)
                        , csm            & ! intent(in)
                        , csh            & ! intent(in)
                        , dl79           & ! intent(in)
                        , ribmax         & ! intent(in)
                        , tprandtl       & ! intent(in)
                        , z0moz0h        & ! intent(in)
                        , z0hoz0m        & ! intent(in)
                        , min_patch_area & ! intent(in)
                        , psim           & ! function
                        , psih           & ! function
                        , zoobukhov      ! ! function
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                  , intent(in)    :: nx
   integer                  , intent(in)    :: ny
   integer                  , intent(in)    :: nz
   integer                  , intent(in)    :: np
   integer                  , intent(in)    :: ng
   character(len=4)         , intent(in)    :: which
   real, dimension(nx,ny,nz), intent(in)    :: theta_atm
   real, dimension(nx,ny,nz), intent(in)    :: rvap_atm
   real, dimension(nx,ny,nz), intent(in)    :: co2_atm
   real, dimension(nx,ny,nz), intent(in)    :: uspd_atm
   real, dimension(nx,ny,np), intent(in)    :: theta_can
   real, dimension(nx,ny,np), intent(in)    :: rvap_can
   real, dimension(nx,ny,np), intent(in)    :: co2_can
   real, dimension(nx,ny,np), intent(in)    :: prss_can
   real, dimension(nx,ny)   , intent(in)    :: topt
   real                     , intent(in)    :: zout
   real, dimension(nx,ny,np), intent(in)    :: rough
   real, dimension(nx,ny,np), intent(inout) :: rib
   real, dimension(nx,ny,np), intent(inout) :: zeta
   real, dimension(nx,ny,np), intent(in)    :: parea
   real, dimension(nx,ny,np), intent(inout) :: ustar
   real, dimension(nx,ny,np), intent(inout) :: tstar
   real, dimension(nx,ny,np), intent(inout) :: rstar
   real, dimension(nx,ny,np), intent(inout) :: cstar
   real, dimension(nx,ny)   , intent(inout) :: varred
   !----- Local variables. ----------------------------------------------------------------!
   integer           :: x            ! Longitude counter
   integer           :: y            ! Latitude counter
   integer           :: p            ! Patch counter
   real              :: zgrd         ! Grid bottom
   real              :: ztop         ! Grid top
   real              :: zref         ! Reference height
   real              :: rtgt         ! Terrain-following coordinate correction factor. 
   real              :: thetav_atm   ! Atmospheric virtual potential temperature
   real              :: thetav_can   ! Canopy air space virtual potential temperature
   logical           :: stable       ! Stable state
   logical           :: is_ed2       ! This is an ED-2 run
   real              :: zroz0m       ! zref/rough(momentum)
   real              :: lnzroz0m     ! ln[zref/rough(momentum)]
   real              :: zroz0h       ! zref/rough(heat)
   real              :: lnzroz0h     ! ln[zref/rough(heat)]
   real              :: zooz0m       ! zout/rough(momentum)
   real              :: lnzooz0m     ! ln[zout/rough(momentum)]
   real              :: zooz0h       ! zout/rough(heat)
   real              :: lnzooz0h     ! ln[zout/rough(heat)]
   real              :: uref         ! Reference wind speed.
   real              :: ured         ! Wind reduced to the level of interest.
   real              :: redp         ! Output variable for this patch.
   real              :: validarea    ! Total area where we have results.
   real              :: stabcorr     ! Correction for very stable cases.
   !----- Local variables, used by L79. ---------------------------------------------------!
   real              :: a2r          ! Drag coefficient in neutral conditions
   real              :: a2o          ! Drag coefficient in neutral conditions
   real              :: fhr          ! Stability parameter for heat at z = zref
   real              :: fmr          ! Stability parameter for momentum at z = zref
   real              :: fho          ! Stability parameter for heat at z = zout
   real              :: fmo          ! Stability parameter for momentum at z = zout
   real              :: c2           ! Part of the c coefficient common to momentum & heat.
   real              :: c3           ! Another auxiliary variable.
   real              :: multh        ! Factor to be multiplied to get the heat/water.
   real              :: cm           ! c coefficient times |Rib|^1/2 for momentum.
   real              :: ch           ! c coefficient times |Rib|^1/2 for heat.
   real              :: ee           ! (z/z0)^1/3 -1. for eqn. 20 w/o assuming z/z0 >> 1.
   !----- Local variables, used by OD95 and/or BH91. --------------------------------------!
   real              :: zetaom       ! (zout + roughness(momentum))/(Obukhov length).
   real              :: zetaoh       ! (zout + roughness(heat)    )/(Obukhov length).
   real              :: zeta0m       ! roughness(momentum)/(Obukhov length).
   real              :: zeta0h       ! roughness(heat)/(Obukhov length).
   real              :: ribold       ! Bulk richardson number.
   !----- External functions. -------------------------------------------------------------!
   real, external    :: cbrt         ! Cubic root
   !---------------------------------------------------------------------------------------!

   !----- Decide whether this is an ED-2 run or not. --------------------------------------!
   is_ed2 = isfcl == 5


   !----- Define grid bottom and top. -----------------------------------------------------!
   zgrd = myztn(2,ng)
   ztop = myzmn(mynnzp(1)-1,1)

   yloop: do y = 1,ny
      xloop: do x = 1,nx
         rtgt = 1. - topt(x,y) / ztop
         zref = zgrd * rtgt
         
         !----- Compute the virtual potential temperature at the model first level. -------!
         thetav_atm = virtt(theta_atm(x,y,2),rvap_atm(x,y,2),rvap_atm(x,y,2))

         !----- Initialise the output variable. -------------------------------------------!
         varred(x,y) = 0.
         validarea   = 0.

         ploop: do p = 1,np
            !----- Skip patch if the area is tiny. ----------------------------------------!
            if (parea(x,y,p) < min_patch_area) cycle ploop

            !----- Compute the virtual pot. temperature at the canopy air space (CAS). ----!
            thetav_can = virtt(theta_can(x,y,p),rvap_can(x,y,p),rvap_can(x,y,p))

            !----- Compute the reference wind speed. --------------------------------------!
            uref = max(uspd_atm(x,y,2),ubmin)

            !----- Re-compute the Richardson number if this is an ED-2 run. ---------------!
            if (is_ed2) then
               rib(x,y,p) = 2.0 * grav * (zref-rough(x,y,p)) * (thetav_atm-thetav_can)     &
                          / ( (thetav_atm+thetav_can) * uref * uref)
            end if
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !     Find the bulk Richardson number and determine whether the layer is       !
            ! stable or not.                                                               !
            !------------------------------------------------------------------------------!
            stable     = rib(x,y,p) > 0.0
            !------------------------------------------------------------------------------!


            !------ Check whether we must correct the bulk Richardson number and stars. ---!
            if (rib(x,y,p) > ribmax .and. myistar /= 1 .and. is_ed2) then
               stabcorr   = ribmax / rib(x,y,p)
               rib(x,y,p) = ribmax
            else
               stabcorr   = 1.0
            end if


            !------------------------------------------------------------------------------!
            !     Find some variables common to all methods.  Notice that, unlike          !
            ! leaf_stars, we here use the output height, not the reference height.         !
            !------------------------------------------------------------------------------!
            zroz0m      = (zref)/rough(x,y,p)
            lnzroz0m    = log(zroz0m)
            zroz0h      = z0moz0h * zroz0m
            lnzroz0h    = log(zroz0h)
            zooz0m      = (zout+rough(x,y,p))/rough(x,y,p)
            lnzooz0m    = log(zooz0m)
            zooz0h      = z0moz0h * zooz0m
            lnzooz0h    = log(zooz0h)
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !     Here we find the standard functions of heat and momentum to integrate    !
            ! the result.                                                                  !
            !------------------------------------------------------------------------------!
            select case (myistar)
            case (1)
               !---------------------------------------------------------------------------!
               !     Here we will use L79 model, the BRAMS default.                        !
               !---------------------------------------------------------------------------!

               !----- Compute the a-square factor and the coefficient to find theta*. -----!
               a2r  = (vonk / lnzroz0m) ** 2.
               a2o  = (vonk / lnzooz0m) ** 2.

               if (stable) then
                  !------------------------------------------------------------------------!
                  !     Stable case.                                                       !
                  !------------------------------------------------------------------------!
                  fmr = 1.0 / (1.0 + (2.0*bl79 * rib(x,y,p) / sqrt(1.0 + dl79*rib(x,y,p))))
                  fhr = 1.0 / (1.0 + (3.0*bl79 * rib(x,y,p) * sqrt(1.0 + dl79*rib(x,y,p))))
                  fmo = fmr
                  fho = fhr

               else
                  !------------------------------------------------------------------------!
                  !     Unstable case.  The only difference from the original method is    !
                  ! that we no longer assume z >> z0, so the "c" coefficient uses the full !
                  ! z/z0 term.                                                             !
                  !------------------------------------------------------------------------!
                  ee  = cbrt(zroz0m) - 1.
                  c2  = bl79 * a2r * ee * sqrt(ee * abs(rib(x,y,p)))
                  cm  = csm * c2
                  ch  = csh * c2
                  fmr = (1.0 - 2.0 * bl79 * rib(x,y,p) / (1.0 + 2.0 * cm))
                  fhr = (1.0 - 3.0 * bl79 * rib(x,y,p) / (1.0 + 3.0 * ch))
                  ee  = cbrt(zooz0m) - 1.
                  c2  = bl79 * a2o * ee * sqrt(ee * abs(rib(x,y,p)))
                  cm  = csm * c2
                  ch  = csh * c2
                  fmo = (1.0 - 2.0 * bl79 * rib(x,y,p) / (1.0 + 2.0 * cm))
                  fho = (1.0 - 3.0 * bl79 * rib(x,y,p) / (1.0 + 3.0 * ch))
               end if

               if (is_ed2) then
                  !----- Re-compute the stars if this is an ED-2 run. ---------------------!
                  ustar(x,y,p) = max(ustmin,uref * sqrt(a2r * fmr))
                  !----- Finding the coefficient to scale the other stars. ----------------!
                  c3 = a2r * uref * fhr / ustar(x,y,p)
                  !----- Computing the other scales. --------------------------------------!
                  rstar(x,y,p) = c3 * (rvap_atm (x,y,2) - rvap_can (x,y,p)   ) * stabcorr
                  tstar(x,y,p) = c3 * (theta_atm(x,y,2) - theta_can(x,y,p)   ) * stabcorr
                  cstar(x,y,p) = c3 * (co2_atm  (x,y,2) - co2_can  (x,y,p)   ) * stabcorr

                  !----- Compute zeta from u* and T* --------------------------------------!
                  zeta(x,y,p)  = grav * vonk * c3 * (thetav_atm - thetav_can)              &
                               / (thetav_atm * ustar(x,y,p) * ustar(x,y,p))
               end if

               ured  = max(0., ustar(x,y,p) * lnzooz0m / (vonk * sqrt(fmo)))
               multh = tprandtl * ustar(x,y,p) * lnzooz0m / (vonk * ured * fho)

            case (2,4)
               !---------------------------------------------------------------------------!
               ! 2. Here we use the model proposed by OD95, the standard for MM5, but with !
               !    some terms that were computed in B71 (namely, the "0" terms). which    !
               !    prevent singularities.  Since we use OD95 to estimate zeta, which      !
               !    avoids the computation of the Obukhov length L , we can't compute      !
               !    zeta0 by its definition(z0/L). However we know zeta, so zeta0 can be   !
               !    written as z0/z * zeta.                                                !
               !                                                                           !
               ! 4. We use the model proposed by BH91, but we find zeta using the          !
               !    approximation given by OD95.                                           !
               !---------------------------------------------------------------------------!
               if (is_ed2) then 

                  !----- We now compute the stability correction functions. ---------------!
                  if (stable) then
                     !----- Stable case. --------------------------------------------------!
                     zeta(x,y,p) = rib(x,y,p) * lnzroz0m / (1.1 - 5.0 * rib(x,y,p))
                  else
                     !----- Unstable case. ------------------------------------------------!
                     zeta(x,y,p) = rib(x,y,p) * lnzroz0m
                  end if
               end if

               zetaom = (zout + rough(x,y,p)) * zeta(x,y,p) / zref
               zetaoh = zetaom
               zeta0m = rough(x,y,p) * zeta(x,y,p) / zref
               zeta0h = zeta0m
               !---------------------------------------------------------------------------!


               !----- Re-compute the stars if this is an ED-2 run. ------------------------!
               if (is_ed2) then
                  ustar(x,y,p) = max (ustmin, vonk * uref                                  &
                                            / (lnzroz0m - psim(zeta(x,y,p),stable,myistar) &
                                                       + psim(zeta0m,stable,myistar)     ))

                  !----- Finding the coefficient to scale the other stars. ----------------!
                  c3    = vonk / (tprandtl * (lnzroz0m - psih(zeta(x,y,p),stable,myistar)  &
                                                      + psih(zeta0m,stable,myistar)     ))
                  !----- Computing the other scales. --------------------------------------!
                  rstar(x,y,p) = c3 * (rvap_atm (x,y,2) - rvap_can (x,y,p)   ) * stabcorr
                  tstar(x,y,p) = c3 * (theta_atm(x,y,2) - theta_can(x,y,p)   ) * stabcorr
                  cstar(x,y,p) = c3 * (co2_atm  (x,y,2) - co2_can  (x,y,p)   ) * stabcorr
               end if
               !---------------------------------------------------------------------------!


               ured  = ustar(x,y,p) * ( lnzooz0m - psim(zetaom,stable,myistar)             &
                                      + psim(zeta0m,stable,myistar) ) / vonk
               multh = tprandtl     * ( lnzooz0m - psih(zetaoh,stable,myistar)             &
                                      + psih(zeta0h,stable,myistar) ) / vonk

            case (3,5)
               !---------------------------------------------------------------------------!
               ! 3. Here we use the model proposed by BH91, which is almost the same as    !
               !    the OD95 method, with the two following (important) differences.       !
               !    a. Zeta (z/L) is actually found using the iterative method.            !
               !    b. Stable functions are computed in a more generic way.  BH91 claim    !
               !       that the oft-used approximation (-beta*zeta) can cause poor         !
               !       ventilation of the stable layer, leading to decoupling between the  !
               !       atmosphere and the canopy air space and excessive cooling.          !
               ! 5. Similar as 3, but we compute the stable functions the same way as      !
               !    OD95.                                                                  !
               !---------------------------------------------------------------------------!
               if (is_ed2) then 
                  !----- Make sure that the bulk Richardson number is not above ribmax. ---!
                  zeta(x,y,p) = zoobukhov(rib(x,y,p),zref,rough(x,y,p),zroz0m,lnzroz0m     &
                                         ,zroz0h,lnzroz0h,stable,myistar)
               end if

               zetaom = (zout + rough(x,y,p)) * zeta(x,y,p) / zref
               zetaoh = zetaom
               zeta0m = rough(x,y,p) * zeta(x,y,p) / zref
               zeta0h = zeta0m
               !---------------------------------------------------------------------------!


               !----- Re-compute the stars if this is an ED-2 run. ------------------------!
               if (is_ed2) then
                  ustar(x,y,p) = max (ustmin, vonk * uref                                  &
                                            / (lnzroz0m - psim(zeta(x,y,p),stable,myistar) &
                                                        + psim(zeta0m,stable,myistar)    ))

                  !----- Finding the coefficient to scale the other stars. ----------------!
                  c3    = vonk / (tprandtl * (lnzroz0m - psih(zeta(x,y,p),stable,myistar)  &
                                                       + psih(zeta0m,stable,myistar)     ))
                  !----- Computing the other scales. --------------------------------------!
                  rstar(x,y,p) = c3 * (rvap_atm (x,y,2) - rvap_can (x,y,p)   ) * stabcorr
                  tstar(x,y,p) = c3 * (theta_atm(x,y,2) - theta_can(x,y,p)   ) * stabcorr
                  cstar(x,y,p) = c3 * (co2_atm  (x,y,2) - co2_can  (x,y,p)   ) * stabcorr
               end if
               !---------------------------------------------------------------------------!


               ured  = ustar(x,y,p) * ( lnzooz0m - psim(zetaom,stable,myistar)             &
                                      + psim(zeta0m,stable,myistar) ) / vonk
               multh = tprandtl     * ( lnzooz0m - psih(zetaoh,stable,myistar)             &
                                      + psih(zeta0h,stable,myistar) ) / vonk

            end select
            !------------------------------------------------------------------------------!


            !----- We now compute the reference variable for this patch. ------------------!
            select case (which)
            case('WIND')
               redp = ured
            case('THET')
               redp = theta_can(x,y,p) + tstar(x,y,p) * multh
            case('TEMP')
               !----- Find the potential temperature. -------------------------------------!
               redp = theta_can(x,y,p) + tstar(x,y,p) * multh
               !----- Convert it to temperature. ------------------------------------------!
               redp = redp * (p00i * prss_can(x,y,p)) ** rocp
            case('RVAP')
               redp = max(toodry,rvap_can(x,y,p)  + rstar(x,y,p) * multh)
            case('CO_2')
               redp = max(toodry,co2_can(x,y,p)   + cstar(x,y,p) * multh)
            case('TDEW')
               !----- Find the mixing ratio. ----------------------------------------------!
               redp = max(toodry,rvap_can(x,y,p)  + rstar(x,y,p) * multh)
               !----- Convert it to vapour partial pressure. ------------------------------!
               redp = prss_can(x,y,p) * redp / (ep + redp)
               !----- Find the dew/frost point. -------------------------------------------!
               redp = tslif(redp)
            case('ZETA')
               redp = zeta(x,y,p)
            case('RICH')
               redp = rib(x,y,p)
            case('USTR')
               redp = ustar(x,y,p)
            case('TSTR')
               redp = tstar(x,y,p)
            case('RSTR')
               redp = rstar(x,y,p)
            case('CSTR')
               redp = cstar(x,y,p)
            end select

            validarea   = validarea   + parea(x,y,p)
            varred(x,y) = varred(x,y) + parea(x,y,p) * redp
         end do ploop

         varred(x,y) = varred(x,y) / validarea
      end do xloop
   end do yloop

   return
end subroutine RAMS_reduced_prop
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine computes the relative humidity.                                      !
!------------------------------------------------------------------------------------------!
subroutine RAMS_comp_relhum(nx,ny,no,tdinrhout,temp)
   use therm_lib, only : eslif
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                     , intent(in)    :: nx
   integer                     , intent(in)    :: ny
   integer                     , intent(in)    :: no
   real, dimension(nx,ny,no)   , intent(inout) :: tdinrhout
   real, dimension(nx,ny,no)   , intent(in)    :: temp
   !----- Local variables. ----------------------------------------------------------------!
   integer                                     :: x
   integer                                     :: y
   integer                                     :: o
   real                                        :: e
   real                                        :: es
   !---------------------------------------------------------------------------------------!

   oloop: do o=1,no
      yloop: do y=1,ny
         xloop: do x=1,nx
            e  = eslif(tdinrhout(x,y,o))
            es = eslif(temp(x,y,o))
            tdinrhout(x,y,o) = max(0.,min(1.,e / es))
         end do xloop
      end do yloop
   end do oloop
   return
end subroutine RAMS_comp_relhum
!==========================================================================================!
!==========================================================================================!




!==========================================================================================!
!==========================================================================================!
!     This subroutine avoids using tiny numbers to some variables, by flushing them to     !
! zero if they are too small.  This avoids FPE, especially when patch integration is about !
! to happen.                                                                               !
!------------------------------------------------------------------------------------------!
subroutine RAMS_flush_to_zero(nx,ny,nz,np,myvar,threshold)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                     , intent(in)    :: nx
   integer                     , intent(in)    :: ny
   integer                     , intent(in)    :: nz
   integer                     , intent(in)    :: np
   real, dimension(nx,ny,nz,np), intent(inout) :: myvar
   real                        , intent(in)    :: threshold
   !----- Local variables. ----------------------------------------------------------------!
   integer                                     :: x
   integer                                     :: y
   integer                                     :: z
   integer                                     :: p
   !---------------------------------------------------------------------------------------!
   ploop: do p=1,np
      zloop: do z=1,nz
         yloop: do y=1,ny
            xloop: do x=1,nx
               if (abs(myvar(x,y,z,p)) < threshold) myvar(x,y,z,p) = 0.
            end do xloop
         end do yloop
      end do zloop
   end do ploop
   return
end subroutine RAMS_flush_to_zero
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!       Will Cheng's code for calculating slp with mm5's GRAPH method.  Added for          !
!  calculating SLP from MM5 algorithm.                                                     !
!    The subroutine calculates SLP from an algorithm taken from  GRAPH, a post-processing  !
! package of MM5 V3.3                                                                      !
!                                                                                          !
!    Input: theta - potential temperature (K)         3D                                   !
!           pp    - Exner function        (J/kg K)    3D                                   !
!           z     - terrain               (m)         2D                                   !
!                                                                                          !
!    Ouput: SLP   - sea-level pressure    (hPa)       2D                                   !
!------------------------------------------------------------------------------------------!
subroutine RAMS_comp_slpmm5(n1,n2,n3,theta,pp,z,slp)
   use rconstants, only : cp   & ! intent(in)
                        , cpi  & ! intent(in)
                        , rdry & ! intent(in)
                        , cpor & ! intent(in)
                        , p00  ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                    , intent(in)    :: n1
   integer                    , intent(in)    :: n2
   integer                    , intent(in)    :: n3
   real, dimension(n1,n2,n3)  , intent(in)    :: theta
   real, dimension(n1,n2,n3)  , intent(in)    :: pp
   real, dimension(n1,n2)     , intent(in)    :: z
   real, dimension(n1,n2)     , intent(inout) :: slp
   !----- Local variables. ----------------------------------------------------------------!
   integer                                    :: i
   integer                                    :: j
   integer                                    :: k
   integer                                    :: kk
   real, dimension(n1,n2)                     :: sfp
   real, dimension(n1,n2)                     :: ts
   real, dimension(n1,n2,n3-1)                :: t_mm5
   real, dimension(n1,n2,n3-1)                :: p_mm5
   !---------------------------------------------------------------------------------------!
   
   do j = 1,n2
      do i = 1,n1
         !----- Calculate surface pressure. -----------------------------------------------!
         sfp(i,j) = (0.5*(pp(i,j,1)+pp(i,j,2))/cp)**cpor*p00*.01
         !----- Calculate surface temp. ---------------------------------------------------!
         ts(i,j)  = 0.5 * cpi * (theta(i,j,1)*pp(i,j,1) + theta(i,j,2)*pp(j,j,2))
      end do
   end do

   do k = 2,n3
      kk = n3-k+1
      do j = 1,n2
         do i = 1,n1
            !----- Flip arrays upside down for input to GRAPH subroutine. -----------------!
            t_mm5(i,j,kk) = theta(i,j,k) * pp(i,j,k) * cpi
            p_mm5(i,j,kk) = (pp(i,j,k) * cpi)**cpor * p00 * .01
         end do
      end do
   end do

   call seaprs_0(t_mm5,p_mm5,z,sfp,ts,n1,n2,n3-1,slp)

   return
end subroutine RAMS_comp_slpmm5
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine computes sea level pressure from the rule                            !
!              t1/t2=(p1/p2)**(gamma*r/g).                                                 !
!                                                                                          !
!     *** Levels go from top-down ***                                                      !
!                                                                                          !
!     Input       t        temperature (Kelvin)                3D                          !
!                 ter      terrain     (m)                     2D                          !
!                 sfp      surface pressure (hPa)              2D                          !
!                 imx      dot point dimension n-s                                         !
!                 jmx      dot point dimension e-w                                         !
!                 kx       number of vertical levels                                       !
!                                                                                          !
!     Output      slp      sea level pressure (hPa)            2D                          !
!                                                                                          !
!------------------------------------------------------------------------------------------!
subroutine seaprs_0(t,pp,ter,sfp,ts,imx,jmx,kx,slp)
   use rconstants, only : rdry & ! intent(in)
                        , grav & ! intent(in)
                        , t00  ! ! intent(in)
   implicit none
   !----- Local constants. ----------------------------------------------------------------!
   real                          , parameter     :: gamma  = 6.5e-3
   real                          , parameter     :: tcrit  = t00+17.5
   real                          , parameter     :: pconst = 100.
   real                          , parameter     :: xterm  = gamma * rdry / grav
   !----- Arguments. ----------------------------------------------------------------------!
   integer                       , intent(in)    :: imx
   integer                       , intent(in)    :: jmx
   integer                       , intent(in)    :: kx
   real   , dimension(imx,jmx,kx), intent(in)    :: t
   real   , dimension(imx,jmx,kx), intent(in)    :: pp
   real   , dimension(imx,jmx)   , intent(in)    :: ter
   real   , dimension(imx,jmx)   , intent(in)    :: sfp
   real   , dimension(imx,jmx)   , intent(inout) :: ts
   real   , dimension(imx,jmx)   , intent(inout) :: slp
   !----- Local variables. ----------------------------------------------------------------!
   integer                                       :: i
   integer                                       :: j
   integer                                       :: k
   integer                                       :: kupto
   integer                                       :: klo
   integer                                       :: khi
   logical                                       :: l1
   logical                                       :: l2
   logical                                       :: l3
   real   , dimension(imx,jmx)                   :: ps
   real   , dimension(imx,jmx)                   :: pl
   real   , dimension(imx,jmx)                   :: t0
   real   , dimension(imx,jmx)                   :: xklev
   real                                          :: xk
   real                                          :: xkhold
   real                                          :: plo
   real                                          :: phi
   real                                          :: tlo
   real                                          :: thi
   real                                          :: tl
   real                                          :: tbar
   real                                          :: t0hold
   real                                          :: hl
   !---------------------------------------------------------------------------------------!

   !------ Compute pressure at pconst mb above surface (pl). ------------------------------!
   kupto=kx/2


   mainloop: do
      do j=1, jmx
         do i=1,imx
            pl(i,j)=sfp(i,j)-pconst
            xklev(i,j)=0.
         end do
      end do

      !----- Find 2 levels on sigma surfaces surrounding pl at each i,j. ------------------!
      jloop: do j=1,jmx
         iloop: do i=1,imx
            kloop: do k=kx-1,kupto,-1
               xk     = real(k)
               xkhold = xklev(i,j)
               xklev(i,j) = merge(xk,xkhold                                                &
                                 ,((pp(i,j,k)   <  pl(i,j)) .and.                          &
                                   (pp(i,j,k+1) >= pl(i,j))       ))
            end do kloop

            if (xklev(i,j) < 1.) then
               write (unit=*,fmt='(a,1x,es12.5,1x,a)')                                     &
                  ' Error finding pressure level ',pconst,' mb above the surface!'
               write (unit=*,fmt='(a,1x,i5,a)') ' Last k level =',kupto,'...'

               if (kupto /= 1) then
                 write (unit=*,fmt='(a)') ' Trying again with kupto=1...'
                 kupto = 1
                 cycle mainloop
               else
                  write(unit=*,fmt='(a,1x,i5,1x)')     ' - I    =',i
                  write(unit=*,fmt='(a,1x,i5,1x)')     ' - J    =',i
                  write(unit=*,fmt='(a,1x,es12.5,1x)') ' - PL   =',pl(i,j)
                  write(unit=*,fmt='(a,1x,es12.5,1x)') ' - PSFC =',sfp(i,j)
                  stop
               end if
            end if
         end do iloop
      end do jloop
      !---- The default is to leave the loop... -------------------------------------------!
      exit mainloop
   end do mainloop

   !---------------------------------------------------------------------------------------!
   !      Get temperature at pl (tl), extrapolate t at surface (ts) and T at sea level     !
   ! (t0) with 6.5 k/km lapse rate.                                                        !
   !---------------------------------------------------------------------------------------!
   jloop2: do j=1,jmx
      iloop2: do i=1,imx
         klo     = nint(xklev(i,j))+1
         khi     = nint(xklev(i,j))
         plo     = pp(i,j,klo)
         phi     = pp(i,j,khi)
         tlo     = t(i,j,klo)
         thi     = t(i,j,khi)
         tl      = thi-(thi-tlo)*alog(pl(i,j)/phi)/alog(plo/phi)
         ts(i,j) = tl*(sfp(i,j)/pl(i,j))**xterm
         tbar    = (ts(i,j)+tl)*0.5
         hl      = ter(i,j)-rdry/grav*alog(pl(i,j)/sfp(i,j))*tbar
         t0(i,j) = tl+gamma*hl
      end do iloop2
   end do jloop2
   !---------------------------------------------------------------------------------------!



   !----- Correct sea level temperature if too hot. ---------------------------------------!
   jloop3: do j=1,jmx
      iloop3: do i=1,imx
         l1      = t0(i,j) <  tcrit
         l2      = ts(i,j) >= tcrit
         l3      = .not. l1
         t0hold  = t0(i,j)

         t0(i,j) = merge(t0hold,merge(tcrit,tcrit-0.005*(ts(i,j)-tcrit)**2,l2.and.l3)      &
                        ,l1.and.l2)
      end do iloop3
   end do jloop3
   !---------------------------------------------------------------------------------------!



   !----- Compute sea level pressure. -----------------------------------------------------!
   jloop4: do j=1,jmx
      iloop4: do i=1,imx
         slp(i,j)=sfp(i,j)*exp(2.*grav*ter(i,j)/(rdry*(ts(i,j)+t0(i,j))))
      end do iloop4
   end do jloop4
   !---------------------------------------------------------------------------------------!

   return
end subroutine seaprs_0
!==========================================================================================!
!==========================================================================================!





!==========================================================================================!
!==========================================================================================!
!     This subroutine computes the relative soil moisture.                                 !
!------------------------------------------------------------------------------------------!
subroutine RAMS_comp_slmstf(nx,ny,nz,np,inwoutf,soil_text)
   use soil_coms, only : soil ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                     , intent(in)    :: nx
   integer                     , intent(in)    :: ny
   integer                     , intent(in)    :: nz
   integer                     , intent(in)    :: np
   real, dimension(nx,ny,nz,np), intent(inout) :: inwoutf
   real, dimension(nx,ny,nz,np), intent(in)    :: soil_text
   !----- Local variables. ----------------------------------------------------------------!
   integer                                     :: x
   integer                                     :: y
   integer                                     :: z
   integer                                     :: p
   integer                                     :: nsoil
   real                                        :: soil_water
   real                                        :: soil_rmois
   !---------------------------------------------------------------------------------------!
   
   do p=2,np
      do z=1,nz
         do y=1,ny
            do x=1,nx
               !----- Copy point to temporary variable. -----------------------------------!
               soil_water = inwoutf(x,y,z,p)

               !----- Find the soil class. ------------------------------------------------!
               nsoil = nint(soil_text(x,y,z,p))

               select case (nsoil)
               case (0) ! Water
                  soil_rmois = 1.
               case default ! Soil
                  soil_rmois = (soil_water         - soil(nsoil)%soilcp)                   &
                             / (soil(nsoil)%slmsts - soil(nsoil)%soilcp)
               end select

               !----- Copy the result back to the array. ----------------------------------!
               inwoutf(x,y,z,p) = soil_rmois
            end do
         end do
      end do
   end do

   return
end subroutine RAMS_comp_slmstf
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine computes the temperature and liquid fraction given the internal      !
! energy and water content.                                                                !
!------------------------------------------------------------------------------------------!
subroutine RAMS_comp_qwtk(nx,ny,nz,np,inqoutt,inwoutl,soil_text)
   use rconstants, only : wdns ! ! intent(in)
   use soil_coms , only : soil ! ! intent(in)
   use therm_lib , only : qwtk ! ! subroutine
   !----- Arguments. ----------------------------------------------------------------------!
   integer                     , intent(in)    :: nx
   integer                     , intent(in)    :: ny
   integer                     , intent(in)    :: nz
   integer                     , intent(in)    :: np
   real, dimension(nx,ny,nz,np), intent(inout) :: inqoutt
   real, dimension(nx,ny,nz,np), intent(inout) :: inwoutl
   real, dimension(nx,ny,nz,np), intent(in)    :: soil_text
   !----- Local variables. ----------------------------------------------------------------!
   integer                                     :: x
   integer                                     :: y
   integer                                     :: z
   integer                                     :: p
   integer                                     :: nsoil
   real                                        :: energy
   real                                        :: water
   real                                        :: dryhcap
   real                                        :: temperature
   real                                        :: fracliq
   !---------------------------------------------------------------------------------------!

   !----- Loop through all grid points, skipping the water patch. -------------------------!
   do p=2,np
      do z=1,nz
         do y=1,ny
            do x=1,nx
               !----- Save variables into temporary places. -------------------------------!
               energy  = inqoutt(x,y,z,p)
               water   = inwoutl(x,y,z,p) * wdns
               nsoil   = nint(soil_text(x,y,z,p))
               dryhcap = soil(nsoil)%slcpd
               !----- Compute temperature and liquid water fraction. ----------------------!
               call qwtk(energy,water,dryhcap,temperature,fracliq)
               !----- Save in the variables that will be returned. ------------------------!
               inqoutt(x,y,z,p) = temperature
               inwoutl(x,y,z,p) = fracliq
            end do
         end do
      end do
   end do

   return
end subroutine RAMS_comp_qwtk
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine RAMS_comp_copysst(nx,ny,nz,inqoutt)
   use therm_lib , only : qtk  ! ! subroutine
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                    , intent(in)    :: nx
   integer                    , intent(in)    :: ny
   integer                    , intent(in)    :: nz
   real, dimension(nx,ny,nz)  , intent(inout) :: inqoutt
   !----- Local variables. ----------------------------------------------------------------!
   integer                                    :: x
   integer                                    :: y
   integer                                    :: z
   real                                       :: energy
   real, dimension(nx,ny)                     :: temperature
   real                                       :: fracliq
   !---------------------------------------------------------------------------------------!

   !----- Copy energy and compute the temperature, saving it into a 2-D array. ------------!
   do y=1,ny
      do x=1,nx
         energy = inqoutt(x,y,nz)
         call qtk(energy,temperature(x,y),fracliq)
      end do
   end do

   !----- Copy temperature to the output array. -------------------------------------------!
   do z=1,nz
      do y=1,ny
         do x=1,nx
            inqoutt(x,y,z) = temperature(x,y)
         end do
      end do
   end do

   return
end subroutine RAMS_comp_copysst
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine RAMS_comp_tempC(nx,ny,nz,np,temp)
   use rconstants, only : t00 ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                     , intent(in)    :: nx
   integer                     , intent(in)    :: ny
   integer                     , intent(in)    :: nz
   integer                     , intent(in)    :: np
   real, dimension(nx,ny,nz,np), intent(inout) :: temp
   !----- Local variables. ----------------------------------------------------------------!
   integer                                    :: x
   integer                                    :: y
   integer                                    :: z
   integer                                    :: p
   !---------------------------------------------------------------------------------------!

   do p=1,np
      do z=1,nz
         do y=1,ny
            do x=1,nx
               temp(x,y,z,p) = temp(x,y,z,p) - t00
            end do
         end do
      end do
   end do
   return
end subroutine RAMS_comp_tempC
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine RAMS_comp_pvap(nx,ny,nz,pres,inroute)
   use rconstants, only : ep ! ! intent(in) 
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                     , intent(in)    :: nx
   integer                     , intent(in)    :: ny
   integer                     , intent(in)    :: nz
   real   , dimension(nx,ny,nz), intent(in)    :: pres
   real   , dimension(nx,ny,nz), intent(inout) :: inroute
   !----- Local variables. ----------------------------------------------------------------!
   integer                                     :: x
   integer                                     :: y
   integer                                     :: z
   real                                        :: rvap
   real                                        :: pvap
   !---------------------------------------------------------------------------------------!

   do x=1,nx
      do y=1,ny
         do z=1,nz
            rvap           = inroute(x,y,z)
            pvap           = pres(x,y,z) * rvap / (ep + rvap)
            inroute(x,y,z) = pvap
         end do
      end do
   end do

   return
end subroutine RAMS_comp_pvap
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine RAMS_comp_spvol(nx,ny,nz,intouta,pvap,pres)
   use rconstants, only : ep   & ! intent(in) 
                        , rdry ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                     , intent(in)    :: nx
   integer                     , intent(in)    :: ny
   integer                     , intent(in)    :: nz
   real   , dimension(nx,ny,nz), intent(inout) :: intouta
   real   , dimension(nx,ny,nz), intent(in)    :: pvap
   real   , dimension(nx,ny,nz), intent(in)    :: pres
   !----- Local variables. ----------------------------------------------------------------!
   integer                                     :: x
   integer                                     :: y
   integer                                     :: z
   real                                        :: temp
   real                                        :: alpha
   !---------------------------------------------------------------------------------------!

   do z=1,nz
      do y=1,ny
         do x=1,nx
            temp           = intouta(x,y,z)
            alpha          = rdry * temp /(pres(x,y,z) - (1.-ep) * pvap(x,y,z))
            intouta(x,y,z) = alpha
         end do
      end do
   end do

   return
end subroutine RAMS_comp_spvol
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine RAMS_comp_theta2temp(nx,ny,nz,inthoutt,press)
   use rconstants, only : p00i & ! intent(in)
                        , rocp ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                     , intent(in)    :: nx
   integer                     , intent(in)    :: ny
   integer                     , intent(in)    :: nz
   real   , dimension(nx,ny,nz), intent(inout) :: inthoutt
   real   , dimension(nx,ny,nz), intent(in)    :: press
   !----- Local variables. ----------------------------------------------------------------!
   integer                                     :: x
   integer                                     :: y
   integer                                     :: z
   real                                        :: theta
   real                                        :: temp
   !---------------------------------------------------------------------------------------!
   do z=1,nz
      do y=1,ny
         do x=1,nx
            theta           = inthoutt(x,y,z)
            temp            = theta * (p00i * press(x,y,z)) ** rocp
            inthoutt(x,y,z) = temp
         end do
      end do
   end do
   return
end subroutine RAMS_comp_theta2temp
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine RAMS_comp_zenith(nx,ny,cosz,zenith)
   use rconstants, only : onerad ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                  , intent(in)    :: nx
   integer                  , intent(in)    :: ny
   real   , dimension(nx,ny), intent(in)    :: cosz
   real   , dimension(nx,ny), intent(inout) :: zenith
   !----- Local variables. ----------------------------------------------------------------!
   integer                                     :: x
   integer                                     :: y
   integer                                     :: z
   real                                        :: theta
   real                                        :: temp
   !---------------------------------------------------------------------------------------!

   do y=1,ny
      do x=1,nx
         zenith(x,y) = acos(cosz(x,y)) * onerad
         if (zenith(x,y) > 90. ) zenith(x,y) = 90.
      end do
   end do

   return
end subroutine RAMS_comp_zenith
!==========================================================================================!
!==========================================================================================!

