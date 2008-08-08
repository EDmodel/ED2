!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

subroutine oda_sta_input (plat,plon,ngds)

use mem_oda
use obs_input
use mem_grid

implicit none

real :: plat,plon
integer :: ngds

integer :: ierr,ns,nf,iexist,ngd,ifile,nloc,idsta,ntimes,ng,nlevs,k
real :: stx,sty,xy_ij
real(kind=8) :: secs_init,secs_obs

integer, parameter :: num_convars=2
character(len=8) :: varn(num_convars)
data varn/'ue','ve'/
real :: vars(num_convars),varp1(1000),varp2(1000)
  
! Read all observation data

! Transfer station id's to structure and assign integer id

do ns=1,num_oda_sfc
   oda_sfc_info(ns)%id=staid_sfc(ns)
   oda_sfc_info(ns)%intid=ns
   oda_sfc_info(ns)%ntimes=0
enddo

! Upper air
do ns=1,num_oda_upa
   oda_upa_info(ns)%id=staid_upa(ns)
   oda_upa_info(ns)%intid=ns
   oda_upa_info(ns)%ntimes=0
enddo

! Get abs seconds of run start

call date_abs_secs2(iyeara,imontha,idatea,itimea*100,secs_init)

! Process surface files

do nf=1,nsfcfiles

   print*,'Opening sfc obs file: ',trim(fnames_sfc(nf))
   
   open(unit=31,file=fnames_sfc(nf),status='old')
      
   ifile=1
   header(ifile)%iun=31
   call rr_sfc_ver (ifile)
   
   do while (.TRUE.)
      call rr_sfc_obs (ifile,'no',ierr)
      if(ierr==1) exit
                 
         ! Find station integer id
      do ns=1,num_oda_sfc
      
         if(trim(rsfc_obs%id) == trim(staid_sfc(ns))) then

            idsta=ns
            oda_sfc_info(ns)%ntimes=oda_sfc_info(ns)%ntimes+1
            ntimes=oda_sfc_info(ns)%ntimes
            
            if(ntimes > maxtimes_sfc) then
               print*,'Surface obs exceeds memory allocation, stopping'
               print*,'idsta,ntimes:',ns,idsta,ntimes
               stop 'too many times for sfc obs'
            endif
            
            if(ntimes == 1) then
               ! Compute locations relative to grids
               call ll_xy (rsfc_obs%lat,rsfc_obs%lon,plat,plon,stx,sty)
               oda_sfc_info(ns)%xsta=stx
               oda_sfc_info(ns)%ysta=sty
               oda_sfc_info(ns)%xlat=rsfc_obs%lat
               oda_sfc_info(ns)%xlon=rsfc_obs%lon
               oda_sfc_info(ns)%stopo=rsfc_obs%elev
               call findgrid (stx,sty,ngd)
               
               ! Assume telescoping grids
               oda_sfc_info(ns)%iactive(1:ngds)=0
               oda_sfc_info(ns)%iactive(1:ngd)=1
               
               ! Find i,j on all active grids
               do ng=1,ngds
                  if (oda_sfc_info(ns)%iactive(ng) == 1) then
                     oda_sfc_info(ns)%xista(ng)=  &
                                      xy_ij (nnxp(ng),xtn(1,ng),stx)
                     oda_sfc_info(ns)%xjsta(ng)=  &
                                      xy_ij (nnyp(ng),ytn(1,ng),sty)
                  else
                     oda_sfc_info(ns)%xista(ng)= 0.
                     oda_sfc_info(ns)%xjsta(ng)= 0.
                  endif
               enddo
            endif
            
            ! Find time in seconds relative to run start         
            call date_abs_secs2(rsfc_obs%jyear,rsfc_obs%jmonth  &
                               ,rsfc_obs%jdate,rsfc_obs%jtime*100,secs_obs)
            oda_sfc_obs(ns)%time(ntimes)=secs_obs - secs_init

            ! Fill data arrays
            
            call sfc_data_convert(vars,varn,num_convars)
            
            oda_sfc_obs(ns)%temp(ntimes)=rsfc_obs%t
            oda_sfc_obs(ns)%dewpt(ntimes)=rsfc_obs%td
            oda_sfc_obs(ns)%us(ntimes)=vars(1)
            oda_sfc_obs(ns)%vs(ntimes)=vars(2)
            oda_sfc_obs(ns)%ps(ntimes)=rsfc_obs%p
            
            oda_sfc_obs(ns)%u(ntimes)=-999.
            oda_sfc_obs(ns)%v(ntimes)=-999.
            if(vars(1) > -998. .and. vars(2) > -998.) then
               call uevetouv(oda_sfc_obs(ns)%u(ntimes)  &
                         ,oda_sfc_obs(ns)%v(ntimes)  &
                         ,vars(1),vars(2)  &
                         ,rsfc_obs%lat,rsfc_obs%lon,plat,plon)
            endif
            if(ns==461) print*,'obs_in:',ntimes,trim(rsfc_obs%id)  &
           ,oda_sfc_info(ns)%xista(1),oda_sfc_info(ns)%xjsta(1) &
           ,oda_sfc_obs(ns)%time(ntimes),oda_sfc_obs(ns)%us(ntimes)  &
           ,oda_sfc_obs(ns)%vs(ntimes),oda_sfc_obs(ns)%temp(ntimes) &
           ,oda_sfc_obs(ns)%u(ntimes),oda_sfc_obs(ns)%v(ntimes) &
           ,oda_sfc_obs(ns)%dewpt(ntimes),rsfc_obs%ff,rsfc_obs%dd
            exit
         endif  
     enddo   
   enddo
   
         
   close(31)

enddo
   
! Process upper air files

   print*,'Opening upa obs file: ',nupafiles
do nf=1,nupafiles

   print*,'Opening upa obs file: ',fnames_upa(nf)(1:len_trim(fnames_upa(nf)))
   
   open(unit=31,file=fnames_upa(nf),status='old')
      
   ifile=1
   header(ifile)%iun=31
   call rr_upa_ver (ifile)
   
   do while (.TRUE.)
      call rr_upa_obs (ifile,'no',ierr)
      if(ierr==1) exit
                 
         ! Find station integer id
      do ns=1,num_oda_upa
         if(trim(rupa_obs%id) == trim(staid_upa(ns))) then

            idsta=ns
            oda_upa_info(ns)%ntimes = oda_upa_info(ns)%ntimes + 1
            ntimes=oda_upa_info(ns)%ntimes
            
            if(ntimes > maxtimes_upa) then
               print*,'Upa obs exceeds memory allocation, stopping'
               print*,'idsta,ntimes:',ns,idsta,ntimes
               stop 'too many times for obs'
            endif
            
            if(ntimes == 1) then
               ! Compute locations relative to grids
               call ll_xy (rupa_obs%lat,rupa_obs%lon,plat,plon,stx,sty)
               oda_upa_info(ns)%xsta=stx
               oda_upa_info(ns)%ysta=sty
               oda_upa_info(ns)%xlat=rupa_obs%lat
               oda_upa_info(ns)%xlon=rupa_obs%lon
               oda_upa_info(ns)%stopo=rupa_obs%elev
               call findgrid (stx,sty,ngd)
               
               ! Assume telescoping grids
               oda_upa_info(ns)%iactive(1:ngds)=0
               oda_upa_info(ns)%iactive(1:ngd)=1
               
               
               ! Find i,j on all active grids
               do ng=1,ngds
                  if (oda_upa_info(ns)%iactive(ng) == 1) then
                     oda_upa_info(ns)%xista(ng)=  &
                                      xy_ij (nnxp(ng),xtn(1,ng),stx)
                     oda_upa_info(ns)%xjsta(ng)=  &
                                      xy_ij (nnyp(ng),ytn(1,ng),sty)
                  else
                     oda_upa_info(ns)%xista(ng)= 0.
                     oda_upa_info(ns)%xjsta(ng)= 0.
                  endif
               enddo
            endif
            
            ! Find time in seconds relative to run start         
            call date_abs_secs2(rupa_obs%jyear,rupa_obs%jmonth  &
                               ,rupa_obs%jdate,rupa_obs%jtime*100,secs_obs)
            oda_upa_obs(ns)%time(ntimes)=secs_obs - secs_init

            ! Fill data arrays
            
            call upa_get_profile(varp1,nlevs,'ue','z')
            call upa_get_profile(varp2,nlevs,'ve','z')
            if (nlevs > maxupalevs) then
               print*, 'ODA error: maxupalevs exceeded for winds:',nlevs
               stop 'ODA maxupalevs winds'
            endif
            do k=1,nlevs
               if(varp1(k) > -998. .and. varp2(k) > -998.) then
                  call uevetouv(oda_upa_obs(ns)%u(k,ntimes)  &
                               ,oda_upa_obs(ns)%v(k,ntimes)  &
                               ,varp1(k),varp2(k)  &
                               ,rupa_obs%lat,rupa_obs%lon,plat,plon)
               endif
            enddo

            call upa_get_profile(oda_upa_obs(ns)%zz(1,ntimes),nlevs,'zz','z')
            oda_upa_obs(ns)%us(1:nlevs,ntimes)=varp1(1:nlevs)
            oda_upa_obs(ns)%vs(1:nlevs,ntimes)=varp2(1:nlevs)
            oda_upa_obs(ns)%lz(ntimes)=nlevs
            
            ! For thermo variables, unlike surface obs, we usually have 
            !    pressure reported in standard raobs. 
            !    Convert temp,prs,rh => theta,pi,r
            !
            ! First get geopotential height, 
            call upa_get_profile(varp2,nlevs,'geo','p')
            if (nlevs > maxupalevs) then
               print*, 'ODA error: maxupalevs exceeded for prs:',nlevs
               stop 'ODA maxupalevs press'
            endif
            oda_upa_obs(ns)%lp(ntimes)=nlevs
            oda_upa_obs(ns)%zgeo(1:nlevs,ntimes)=varp2(1:nlevs)
            
            call upa_get_profile(varp1,nlevs,'theta','p')
            oda_upa_obs(ns)%theta(1:nlevs,ntimes)=varp1(1:nlevs)
            
            call upa_get_profile(varp1,nlevs,'pi','p')
            oda_upa_obs(ns)%pi(1:nlevs,ntimes)=varp1(1:nlevs)
            
            call upa_get_profile(varp1,nlevs,'mixrat','p')
            oda_upa_obs(ns)%rv(1:nlevs,ntimes)=varp1(1:nlevs)
            
!print*,plat,plon,ns,ntimes
!do k=1,oda_upa_obs(ns)%lp(ntimes)
!   print '(i3,5f12.3)',k,oda_upa_obs(ns)%pi(k,ntimes)  &
!               ,oda_upa_obs(ns)%zgeo(k,ntimes) &
!               ,oda_upa_obs(ns)%theta(k,ntimes) &
!               ,oda_upa_obs(ns)%rv(k,ntimes)
!enddo
!do k=1,oda_upa_obs(ns)%lz(ntimes)
!   print '(i3,5f12.3)',k,oda_upa_obs(ns)%u(k,ntimes)  &
!               ,oda_upa_obs(ns)%v(k,ntimes) &
!               ,oda_upa_obs(ns)%us(k,ntimes) &
!               ,oda_upa_obs(ns)%vs(k,ntimes) &
!               ,oda_upa_obs(ns)%zz(k,ntimes)
!enddo

            exit
         endif  
      enddo   
   enddo
        
   close(31)


enddo

return
end



real function xy_ij (nih,xh,x)

! gets i (or j) grid point index from x (or y) value in m

implicit none

integer :: nih
real :: xh(nih),x

integer :: inti,m
real :: reali

if(x < xh(1) .or. x > xh(nih))then
   print*, 'x, y, or z value exceeds grid limits'
   stop 'xy_ij'
endif

!print*,x
!print*,xh

do m=2,nih
   if (x >= xh(m-1) .and. x < xh(m)) then
      inti=m-1
      reali=(x-xh(m-1))/(xh(m)-xh(m-1))
      xy_ij=float(inti)+reali
    !  print*,'i,j:',inti,reali,xy_ij
      exit
   endif
enddo

return
end

