!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################


subroutine oda_nudge()

use mem_oda
use mem_grid
use mem_scratch
use mem_tend
use mem_basic
use io_params

use node_mod

implicit none

integer :: ncall=0, nobs, ng,i,j,k
real::valugp
real, allocatable, dimension(:,:) :: plt

!!!!print*,'start nud'
!!!!!!!!!!!!!!!!! Modifications needed for ADAP


! Namelist variables:

! - Flag to turn on oda                     (IF_ODA)
! - Model start and end time for oda        (TODABEG,TODAEND)
! - Frequency of tendency updates           (FRQODA) (0. for every coarse grid dt)
! - For each grid:
!     - Nudging timescale                   (TNUDODA)
!     - Number of vertical levels to apply  (no)
!     - Flag to turn on for this grid       (IG_ODA)


! Covariance stuff:

! For each grid:
! - Surface radii - 1) e-2  2) zero         (RODA_SFCE, RODA_SFC0) (SFC0 <= UPA0)
! - "mid-level" radii - same                (RODA_UPAE, RODA_UPA0)
! - "mid-level" radii height                (RODA_HGT)
!         (linear change from surface to mid, constant above)
!
! - vertical "100" factor                   (RODA_ZFACT)

! Obs processing:
! - Time interpolate limit (TIL)- if the future-past obs time 
!    is > this limit, do not use to interpolate (ODA_SFC_TIL,ODA_UPA_TIL)
!
! - Time extrapolate limit (TEL)- if past/future obs is greater than TIL,
!    but less than TEL, use obs                 (ODA_SFC_TEL,ODA_UPA_TEL)



if (ncall == 0) then
   call oda_nudge_init(ngrids,mmzp(1),mmxp(1),mmyp(1))
   ncall = 1
endif

!print*,'time,frqoda:',mynum,time,frqoda,mod(time,frqoda),dtlongn(1)

if (frqoda > 0. .and. mod(time,dble(frqoda)) >= dtlongn(1) ) go to 20


if ( (ngrid == 1 .and. time >= todabeg .and. time <= todaend)  &
      .or. time == timstr ) then
      
   ! Compute new krigged fields and variances

   do ng=1,ngrids
   
      if (wt_oda_grid(ng) > 0.0) then
   
      ! Process observations and call krigging 

      call oda_proc_obs(mmzp(ng),mmxp(ng),mmyp(ng),mi0(ng),mj0(ng)  &
            ,basic_g(ng)%pp,basic_g(ng)%pi0,scratch%scr1,ng,nobs)
            
            !!!!print*,'oda-sfc nobs:',ng,nobs

      !do i=1,nobs
      !print*,'uobs:',i,ukobs(i),oda_sfc_obs(i)%u(1:oda_sfc_info(i)%ntimes)
      !enddo
      oda_g(ng)%uk(1:mmzp(ng),1:mmxp(ng),1:mmyp(ng))=0.
      oda_g(ng)%ukv(1:mmzp(ng),1:mmxp(ng),1:mmyp(ng))=0.
      call krig(mmzp(ng),mmxp(ng),mmyp(ng)  &
               ,xtn(1+mi0(ng),ng),ytn(1+mj0(ng),ng),ztn(1,ng)  &
               ,nobs,xkobs,ykobs,zkobs,ekobs,ukobs  &
               ,ng,nnzp(ng),grid_g(ng)%topt  &
               ,oda_g(ng)%uk,oda_g(ng)%ukv,1.)

      oda_g(ng)%vk(1:mmzp(ng),1:mmxp(ng),1:mmyp(ng))=0.
      oda_g(ng)%vkv(1:mmzp(ng),1:mmxp(ng),1:mmyp(ng))=0.
      call krig(mmzp(ng),mmxp(ng),mmyp(ng)  &
               ,xtn(1+mi0(ng),ng),ytn(1+mj0(ng),ng),ztn(1,ng)  &
               ,nobs,xkobs,ykobs,zkobs,ekobs,vkobs  &
               ,ng,nnzp(ng),grid_g(ng)%topt,oda_g(ng)%vk,oda_g(ng)%vkv,1.)

      !do i=1,nobs
      !if(tkobs(i) > 200. .and. tkobs(i) < 273.) print*,'tobs:',ng,i,tkobs(i)
      !enddo
      oda_g(ng)%tk(1:mmzp(ng),1:mmxp(ng),1:mmyp(ng))=0.
      oda_g(ng)%tkv(1:mmzp(ng),1:mmxp(ng),1:mmyp(ng))=0.
      call krig(mmzp(ng),mmxp(ng),mmyp(ng)  &
               ,xtn(1+mi0(ng),ng),ytn(1+mj0(ng),ng),ztn(1,ng)  &
               ,nobs,xkobs,ykobs,zkobs,ekobs,tkobs  &
               ,ng,nnzp(ng),grid_g(ng)%topt  &
               ,oda_g(ng)%tk,oda_g(ng)%tkv,1.)
    
      !if(mynum == 2) then
      !print*,'nobs:',nobs,xtn(1+mi0(ng),ng),ytn(1+mj0(ng),ng)
      !do i=1,nobs
      !   if(xkobs(i) > xtn(1+i0,ng) .and. tkobs(i) > 0.) &
      !      print*,'obs:',i,tkobs(i),xkobs(i)
      !enddo

      !allocate(plt(mmxp(ng),mmyp(ng)))
      !call opngks
      !do k=1,nnzp(ng)
      !print*,'plotting:',mynum,k
      !do j=1,mmyp(ng)
      !   do i=1,mmxp(ng)
      !       !if(k==2)print*,'t,tv:',k,i,j,oda_g(ng)%tk(k,i,j),oda_g(ng)%tkv(k,i,j)
      !      plt(i,j)=oda_g(ng)%tk(k,i,j)
      !   enddo
      !enddo
      !plt(1:mmxp(ng),1:mmyp(ng))=oda_g(ng)%tk(k,1:mmxp(ng),1:mmyp(ng))
      !call ezcntr(plt,mmxp(ng),mmyp(ng))
      !do j=1,mmyp(ng)
      !   do i=1,mmxp(ng)
      !      plt(i,j)=oda_g(ng)%tkv(k,i,j)
      !   enddo
      !enddo
      !plt(1:mmxp(ng),1:mmyp(ng))=oda_g(ng)%tkv(k,1:mmxp(ng),1:mmyp(ng))
      !call ezcntr(plt,mmxp(ng),mmyp(ng))
      !enddo
      !call clsgks
      !endif
    
      oda_g(ng)%rk(1:mmzp(ng),1:mmxp(ng),1:mmyp(ng))=0.
      oda_g(ng)%rkv(1:mmzp(ng),1:mmxp(ng),1:mmyp(ng))=0.
      call krig(mmzp(ng),mmxp(ng),mmyp(ng)  &
               ,xtn(1+mi0(ng),ng),ytn(1+mj0(ng),ng),ztn(1,ng)  &
               ,nobs,xkobs,ykobs,zkobs,ekobs,rkobs  &
               ,ng,nnzp(ng),grid_g(ng)%topt  &
               ,oda_g(ng)%rk,oda_g(ng)%rkv,1.)
      endif
      
   enddo

endif

20 continue
!print*,'here1'
! Compute and apply tendencies

ng=ngrid
if (wt_oda_grid(ng) > 0.0 .and. time >= todabeg .and. time <= todaend) then
      if(allocated(plt)) deallocate(plt);allocate(plt(nnxp(ng),nnyp(ng)))

call oda_tendency(mmzp(ng),mmxp(ng),mmyp(ng),mia(ng),miz(ng),mja(ng),mjz(ng)  &
                 ,basic_g(ng)%up,tend%ut,oda_g(ng)%uk,oda_g(ng)%ukv  &
                 ,wt_oda_uv * wt_oda_grid(ng)/tnudoda,time,mi0(ng),mj0(ng))
!      do j=1,nnyp(ng)
!      do i=1,nnxp(ng)
      
!      plt(i,j)=valugp(nnzp(ng),nnxp(ng),nnyp(ng),22,i,j,tend%ut(1))
!      enddo
!      enddo
    !  call ezcntr(plt,nnxp(ng),nnyp(ng))
call oda_tendency(mmzp(ng),mmxp(ng),mmyp(ng),mia(ng),miz(ng),mja(ng),mjz(ng)  &
                 ,basic_g(ng)%vp,tend%vt,oda_g(ng)%vk,oda_g(ng)%vkv  &
                 ,wt_oda_uv * wt_oda_grid(ng)/tnudoda,time,mi0(ng),mj0(ng))
!      do j=1,nnyp(ng)
!      do i=1,nnxp(ng)
      
!      plt(i,j)=valugp(nnzp(ng),nnxp(ng),nnyp(ng),22,i,j,tend%vt(1))
!      enddo
!      enddo
   !   call ezcntr(plt,nnxp(ng),nnyp(ng))
call oda_tendency(mmzp(ng),mmxp(ng),mmyp(ng),mia(ng),miz(ng),mja(ng),mjz(ng)  &
                 ,basic_g(ng)%theta,tend%tht,oda_g(ng)%tk,oda_g(ng)%tkv  &
                 ,wt_oda_th * wt_oda_grid(ng)/tnudoda,time,mi0(ng),mj0(ng))
!print*,'ttttt tend:',time,minval(tend%tht(1:mmzp(ng)*mmxp(ng)*mmyp(ng))) &
!                     ,maxval(tend%tht(1:mmzp(ng)*mmxp(ng)*mmyp(ng)))
      !do j=1,nnyp(ng)
      !do i=1,nnxp(ng)
     
      !plt(i,j)=valugp(nnzp(ng),nnxp(ng),nnyp(ng),2,i,j,tend%tht(1))
      !enddo
      !enddo
      !call ezcntr(plt,nnxp(ng),nnyp(ng))
call oda_tendency(mmzp(ng),mmxp(ng),mmyp(ng),mia(ng),miz(ng),mja(ng),mjz(ng)  &
                 ,basic_g(ng)%rtp,tend%rtt,oda_g(ng)%rk,oda_g(ng)%rkv         &
                 ,wt_oda_rt * wt_oda_grid(ng)/tnudoda,time,mi0(ng),mj0(ng))


endif

return
end

!========================================================================

subroutine oda_nudge_init(ngg,m1m,m2m,m3m)

use mem_oda
use mem_grid

implicit none

integer :: ngg
integer, dimension(ngg) :: m1m,m2m,m3m

integer ::ng,k

! Turn namelist parameters into krigging routine parameters
print*,'nud_init:'

do ng=1,ngrids
   if (wt_oda_grid(ng) > 0.0) then
      nstkrg(ng)=1
      ckrg(1:nstkrg(ng),ng) = 1
      do k=1,nnzp(ng)
         if(ztn(k,ng) < roda_hgt(ng)) then
            rmaxkrg(k,ng)=   &
               roda_sfc0(ng)+(roda_upa0(ng)-roda_sfc0(ng))  &
                             * ztn(k,ng)/roda_hgt(ng)
            ! Kriging routine needs the following to be negative...
            akrg(k,ng)   = -  &
               (roda_sfce(ng)+(roda_upae(ng)-roda_sfce(ng))  &
                             * ztn(k,ng)/roda_hgt(ng)) 
         else
            rmaxkrg(k,ng)=roda_upa0(ng)
            akrg(k,ng) = -roda_upae(ng)
         endif
         print*,k,ztn(k,ng), roda_hgt(ng),rmaxkrg(k,ng),akrg(k,ng)
      enddo
      
      caxkrg(1,ng)= 1. ; caykrg(1,ng)= 0. ; cazkrg(1,ng)= 0.
      caxkrg(2,ng)= 0. ; caykrg(2,ng)= 1. ; cazkrg(2,ng)= 0.
      caxkrg(3,ng)= 0. ; caykrg(3,ng)= 0. ; cazkrg(3,ng)= roda_zfact(ng)
      
      ! Initialize analysis and variance fields so zero tendency will be computed

      oda_g(ng)%uk(1:m1m(ng),1:m2m(ng),1:m3m(ng))=0.
      oda_g(ng)%vk(1:m1m(ng),1:m2m(ng),1:m3m(ng))=0.
      oda_g(ng)%tk(1:m1m(ng),1:m2m(ng),1:m3m(ng))=0.
      oda_g(ng)%rk(1:m1m(ng),1:m2m(ng),1:m3m(ng))=0.
      oda_g(ng)%ukv(1:m1m(ng),1:m2m(ng),1:m3m(ng))=2.
      oda_g(ng)%vkv(1:m1m(ng),1:m2m(ng),1:m3m(ng))=2.
      oda_g(ng)%tkv(1:m1m(ng),1:m2m(ng),1:m3m(ng))=2.
      oda_g(ng)%rkv(1:m1m(ng),1:m2m(ng),1:m3m(ng))=2.
      
   endif
enddo

return
end

!========================================================================

subroutine oda_tendency(n1,n2,n3,ia,iz,ja,jz,ap,at,ak,akv,tscalei,time,i0,j0)

use mem_oda

implicit none

integer :: n1,n2,n3,ia,iz,ja,jz,i0,j0
real, dimension(n1,n2,n3) :: ap,at,ak,akv
real :: tscalei
real(kind=8) :: time
real :: ttmin,ttmax

integer :: i,j,k,iiix,jjjx,kkkx,iiim,jjjm,kkkm
real :: fnna

! Compute and apply tendencies based on krigged field/variance 

!print*,'oda_tend:',ia,iz,ja,jz,n1,n2,n3

ttmax=-1e20; ttmin=1e20
do j=ja,jz
   do i=ia,iz
      do k=1,n1
         fnna=(1.0-akv(k,i,j)/2.0)*(ak(k,i,j)-ap(k,i,j))*tscalei
         at(k,i,j)=at(k,i,j)+fnna
         if(fnna > ttmax) then
            ttmax=fnna
            iiix=i+i0
            jjjx=j+j0
            kkkx=k
         endif
         if(fnna < ttmin) then
            ttmin=fnna
            iiim=i+i0
            jjjm=j+j0
            kkkm=k
         endif
      enddo
   enddo
enddo

! print*,'fffffff max:',time,ttmax,iiix,jjjx,kkkx,i0,j0
! print*,'fffffff min:',time,ttmin,iiim,jjjm,kkkm,i0,j0

return
end
