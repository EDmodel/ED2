!==========================================================================================!
! rmass.f90                                                                                !
!                                                                                          !
!     This file contains the subroutines to compute mass-flux related stuff, as well as    !
! the averaging for advection and turbulence variables                                     !
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
! Subroutine prep_advflx_to_mass                                                           !
! Developed by Saulo R. Freitas (CPTEC/INPE)                                               !
!                                                                                          !
!   This subroutine prepares the advective fluxes to be used in Lagrangian models).        !
!------------------------------------------------------------------------------------------!
subroutine prep_advflx_to_mass(mzp,mxp,myp,ia,iz,ja,jz,ng)

   use mem_grid   ,  only: dtlt
   use mem_scratch,  only: scratch
   use mem_mass   ,  only: mass_g,frqmassave

   implicit none

   integer, intent(in) :: mzp
   integer, intent(in) :: mxp
   integer, intent(in) :: myp
   integer, intent(in) :: ia
   integer, intent(in) :: iz
   integer, intent(in) :: ja
   integer, intent(in) :: jz
   integer, intent(in) :: ng

   real                :: dtlti
   real                :: frqmassi

   dtlti    = 1./dtlt
   frqmassi = 1./frqmassave

   call compute_mass_flux(mzp,mxp,myp,ia,iz,ja,jz,dtlti,frqmassi                           &
           ,scratch%vt3da(1)         , scratch%vt3db(1)        , scratch%vt3dc(1)          &
           ,mass_g(ng)%afxu  (1,1,1) , mass_g(ng)%afxv  (1,1,1), mass_g(ng)%afxw  (1,1,1)  &
           ,mass_g(ng)%afxub (1,1,1) , mass_g(ng)%afxvb (1,1,1), mass_g(ng)%afxwb (1,1,1)  )

   return
end subroutine prep_advflx_to_mass
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
! Subroutine compute_mass_flux                                                             !
! Based on original Saulo R. Freitas (CPTEC/INPE) subroutine                               !
!                                                                                          !
! This subroutine compute the integrated mass flux from advection.                         !
!------------------------------------------------------------------------------------------!
subroutine compute_mass_flux(mzp,mxp,myp,ia,iz,ja,jz,dtlti,frqmassi,vt3da,vt3db,vt3dc,afxu &
                            ,afxv,afxw,afxub,afxvb,afxwb)
   implicit none
   integer                        , intent(in)    :: mzp,mxp,myp
   integer                        , intent(in)    :: ia,iz,ja,jz
   real                           , intent(in)    :: dtlti,frqmassi
   real   , dimension(mzp,mxp,myp), intent(in)    :: vt3da,vt3db,vt3dc
   real   , dimension(mzp,mxp,myp), intent(out)   :: afxu,afxv,afxw
   real   , dimension(mzp,mxp,myp), intent(inout) :: afxub,afxvb,afxwb

   integer                                         :: i,j,k
   do k=1,mzp
      do i=ia,iz
         do j=ja,jz
            afxu(k,i,j)  =                vt3da(k,i,j) *dtlti
            afxv(k,i,j)  =                vt3db(k,i,j) *dtlti
            afxw(k,i,j)  =                vt3dc(k,i,j) *dtlti
            afxub(k,i,j) = afxub(k,i,j) + vt3da(k,i,j) * frqmassi
            afxvb(k,i,j) = afxvb(k,i,j) + vt3db(k,i,j) * frqmassi
            afxwb(k,i,j) = afxwb(k,i,j) + vt3dc(k,i,j) * frqmassi
         end do
      end do
   end do

   return
end subroutine compute_mass_flux
!==========================================================================================!
!==========================================================================================!








!------------------------------------------------------------------------------------------!
! Subroutine prep_convflx_to_mass                                                          !
! Developed by Saulo R. Freitas (CPTEC/INPE)                                               !
!                                                                                          !
!   This subroutine prepares the convective fluxes to be used in STILT (or other           !
! Lagrangian models).                                                                      !
!------------------------------------------------------------------------------------------!
subroutine prep_convflx_to_mass(m1,m2,m3,ia,iz,ja,jz,maxiens,ifm                          &
                                ,ierr4d,jmin4d,kdet4d,k224d,kbcon4d,ktop4d                 &
                                ,kpbl4d,kstabi4d,kstabm4d,xmb4d,edt4d,zcup5d,pcup5d,enup5d &
                                ,endn5d,deup5d,dedn5d,zup5d,zdn5d,iens)

use mem_mass

implicit none

integer, intent(in) :: ifm,iens, maxiens,m1,m2,m3,ia,iz,ja,jz

integer, intent(in),dimension(m2,m3,maxiens) ::                                            &
                         ierr4d,jmin4d,kdet4d,k224d,kbcon4d,ktop4d,kpbl4d,kstabi4d,kstabm4d
                       
real, intent(in), dimension(m2,m3,maxiens) :: xmb4d,edt4d
                        
real, intent(in), dimension(m1,m2,m3,maxiens) ::                        &
                                      enup5d,endn5d,deup5d,dedn5d,zup5d,zdn5d,zcup5d,pcup5d

integer :: i, j


!zero out all convflx
if (iens == 1) then
!Deep conv
  call azero(m1*m2*m3,mass_g(ifm)%cfxup1(1,1,1))
  call azero(m1*m2*m3,mass_g(ifm)%cfxdn1(1,1,1))
  call azero(m1*m2*m3,mass_g(ifm)%dfxup1(1,1,1))
  call azero(m1*m2*m3,mass_g(ifm)%efxup1(1,1,1))
  call azero(m1*m2*m3,mass_g(ifm)%dfxdn1(1,1,1))
  call azero(m1*m2*m3,mass_g(ifm)%efxdn1(1,1,1))
!shallow conv
elseif (iens == 2) then
 call azero(m1*m2*m3,mass_g(ifm)%cfxup2(1,1,1))
 call azero(m1*m2*m3,mass_g(ifm)%dfxup2(1,1,1))       
 call azero(m1*m2*m3,mass_g(ifm)%efxup2(1,1,1))
end if

do j=ja,jz
  do i=ia,iz
!    if((iens == 1 .and. ierr4d(i,j,iens,ngrid) == 0) .or. iens == 2) then
    if(ierr4d(i,j,iens) == 0) then
      call get_convflx(iens,i,j,m1,m2,m3                                       &
                  ,   xmb4d(i,j,iens),   edt4d(i,j,iens)                       &
                  ,  jmin4d(i,j,iens),  kdet4d(i,j,iens)                       &
                  ,   k224d(i,j,iens), kbcon4d(i,j,iens)                       &
                  ,  ktop4d(i,j,iens),  kpbl4d(i,j,iens)                       &
                  ,kstabi4d(i,j,iens),kstabm4d(i,j,iens)                       &
                  ,zcup5d(1:m1,i,j,iens),pcup5d(1:m1,i,j,iens)                 &
                  ,deup5d(1:m1,i,j,iens),enup5d(1:m1,i,j,iens)                 &
                  ,dedn5d(1:m1,i,j,iens),endn5d(1:m1,i,j,iens)                 &
                  , zup5d(1:m1,i,j,iens), zdn5d(1:m1,i,j,iens)                 &
                  ,mass_g(ifm)%cfxup1(1,1,1),mass_g(ifm)%cfxdn1(1,1,1)   &
                  ,mass_g(ifm)%dfxup1(1,1,1),mass_g(ifm)%efxup1(1,1,1)   &
                  ,mass_g(ifm)%dfxdn1(1,1,1),mass_g(ifm)%efxdn1(1,1,1)   &
                  ,mass_g(ifm)%cfxup2(1,1,1),mass_g(ifm)%dfxup2(1,1,1)   &
                  ,mass_g(ifm)%efxup2(1,1,1))
    endif 
  enddo 
enddo  

return
end subroutine prep_convflx_to_mass
!------------------------------------------------------------------------------------------!

!------------------------------------------------------------------------------------------!
! Subroutine get_convflx                                                                   !
! Developed by Saulo R. Freitas (CPTEC/INPE)                                               !
!                                                                                          !
!   This subroutine aims at getting the convective fluxes from Grell shallow and           !
! deep convective parameterizations.                                                       !
!------------------------------------------------------------------------------------------!
subroutine get_convflx(iens,i,j,m1,m2,m3,xmb,edt,jmin,kdet,k22,kbcon,ktop,kpbl             &
                      ,kstabi,kstabm,z_cup,p_cup,cd,entr,cdd,entrd,zu,zd,cfxup1,cfxdn1     &
                      ,dfxup1,efxup1,dfxdn1,efxdn1,cfxup2,dfxup2,efxup2)

implicit none

integer, intent(in)                      ::  iens,i,j,m1,m2,m3,jmin,kdet,k22,kbcon   &
                                            ,ktop,kpbl,kstabi,kstabm
                                            
real, intent(in)                         ::  xmb,edt

real, intent(in),    dimension(m1)       ::  z_cup,p_cup,cd,entr,cdd,entrd,zu,zd

real, intent(inout), dimension(m1,m2,m3) ::  cfxup1,cfxdn1,dfxup1,efxup1,dfxdn1,efxdn1     &
			                    ,cfxup2,dfxup2,efxup2  
                                             
integer                                  ::  k,kr
real                                     ::  dz,totmas,entup,detup,entdoj,entupk,detupk    &
                                            ,detdo,entdo,subdown,subin,detdo1,detdo2

do k=2,ktop+1
        
  kr= k + 1   ! level K of conv grid  corresponds to level K + 1 of RAMS grid
  dz =  z_cup(kr) - z_cup(k)
   
  entup   = 0.
  detup   = 0.
  entdoj  = 0.
  entupk  = 0.
  detupk  = 0.
  detdo   = edt*cdd(k)*dz*zd(kr)
  entdo   = edt*entrd(k)*dz*zd(kr)
  subdown = zu(k ) - edt*zd(k )
  subin   = zu(kr) - edt*zd(kr)

  if (k >= kbcon .and. k < ktop) then
    entup  = entr(k)   *dz*zu(k)
    detup  =   cd(kr) *dz*zu(k)
  end if

  if(k == jmin)  entdoj = edt*zd(k)
  if(k == k22-1) entupk = zu(kpbl)
  if(k == ktop)  detupk = zu(ktop)
  if(k > kdet)   detdo  = 0.
  if(k == ktop)  subin  = 0.
  if(k < kbcon)  detup  = 0.
  
  if(iens == 1) then ! Deep convection
      cfxup1(k ,i,j) =     xmb* zu(k)
      cfxdn1(k ,i,j) =-edt*xmb* zd(k)
      dfxup1(kr,i,j) =     xmb*(detup + detupk)
      efxup1(kr,i,j) =    -xmb*(entup + entupk)
      dfxdn1(kr,i,j) =     xmb*(detdo         ) !edt already is at detdo
      efxdn1(kr,i,j) =    -xmb*(entdo + entdoj) !edt already is at entdo,entdoj
  elseif(iens == 2)  then ! Shallow convection
      cfxup2(k ,i,j) =     xmb* zu(k)
      dfxup2(kr,i,j) =     xmb*(detup + detupk)
      efxup2(kr,i,j) =    -xmb*(entup + entupk)
  end if
!------------------------------------------------------------------------------------------!
! Checking the mass conservation                                                           !
!------------------------------------------------------------------------------------------!
  totmas=subin-subdown+detup-entup-entdo+detdo-entupk-entdoj+detupk
  if(abs(totmas) > 1.e-6) then
    write (unit=*,fmt='(a)')                 '----------- Subroutine Get_convflx ----------'
    write(unit=*, fmt='(4(a,1x,i3,1x))')     '  K= ',k,'   I=',i,'   J=',j,'   IENS=',iens
    write(unit=*, fmt='(2(a,1x,es10.3,1x))') '  subin=  ',    subin,'subdown=',subdown
    write(unit=*, fmt='(2(a,1x,es10.3,1x))') '  detup=  ',    detup,'entup=  ',entup
    write(unit=*, fmt='(2(a,1x,es10.3,1x))') '  entdo=  ',    entdo,'detdo=  ',detdo
    write(unit=*, fmt='(2(a,1x,es10.3,1x))') '  entupk= ',   entupk,'detupk= ',detupk
    write(unit=*, fmt='(2(a,1x,es10.3,1x))') '  zu(k)=  ',    zu(k),'zd(k)=  ',zd(k)
    write(unit=*, fmt='(2(a,1x,es10.3,1x))') '  zu(kr)= ',   zu(kr),'zd(kr)= ',zd(kr)
    write(unit=*, fmt='(2(a,1x,es10.3,1x))') '  entdoj= ',   entdoj,'edt=    ',edt
    write(unit=*, fmt='(1(a,1x,es10.3,1x))') '  totmas= ',   totmas
    write(unit=*, fmt='(a)')                 '---------------------------------------------'
    call abort_run('The model will stop since it is not conserving mass...' &
                  ,'get_convfx','rmass.f90')
  end if

end do

!------------------------------------------------------------------------------------------!
! Bottom layer                                                                             !
!------------------------------------------------------------------------------------------!
if(iens == 1) then  ! Deep convection
  k = 1
  kr= k + 1  ! the K-level of Grell is equivalent to the BRAMS K+1-level

  dz        =  z_cup(2)-z_cup(1)

  detdo1    = edt*zd(2)*  cdd(1)*dz
  detdo2    = edt*zd(1)
  entdo     = edt*zd(2)*entrd(1)*dz
  subin     =-edt*zd(2)           

  cfxup1(kr,i,j) = 0.
  cfxdn1(k,i,j) =-edt*xmb* zd(1)
  dfxup1(kr,i,j) = 0.
  efxup1(kr,i,j) = 0.
  dfxdn1(kr,i,j) = xmb*(detdo1+detdo2) !edt already is at detdo1,2
  efxdn1(kr,i,j) =-xmb* entdo          !edt already is at entdo
 

!------------------------------------------------------------------------------------------!
! Checking the mass conservation                                                           !
!------------------------------------------------------------------------------------------!
  totmas = detdo1+detdo2-entdo+subin
  if(abs(totmas) > 1.e-6) then
    write (unit=*,fmt='(a)')                 '----------- Subroutine Get_convflx ----------'
    write(unit=*, fmt='(4(a,1x,i3,1x))')     '  K= ',k,'   I=',i,'   J=',j,'   IENS=',iens
    write(unit=*, fmt='(2(a,1x,es10.3,1x))') '  subin=  ',    subin,'entdo=  ',entdo
    write(unit=*, fmt='(2(a,1x,es10.3,1x))') '  detdo1= ',   detdo1,'detdo2= ',detdo2
    write(unit=*, fmt='(1(a,1x,es10.3,1x))') '  totmas= ',   totmas
    write(unit=*, fmt='(a)')                 '---------------------------------------------'
    call abort_run('The model will stop since it is not conserving mass...' &
                  ,'get_convfx','rmass.f90')
  end if
end if

return
end subroutine get_convflx
!------------------------------------------------------------------------------------------!






!==========================================================================================!
!==========================================================================================!
! Subroutine prepare_tke_to_mass                                                           !
!                                                                                          !
!   The aim of this subroutine is simply save the mean TKE value at the regular and lite   !
! analysis.                                                                                !
!------------------------------------------------------------------------------------------!
subroutine prepare_tke_to_mass(m1,m2,m3,ia,iz,ja,jz,dtlt,tkep,tkepb)
   use mem_mass, only : frqmassave
   implicit none
   integer, intent(in)                      :: m1,m2,m3
   integer, intent(in)                      :: ia,iz,ja,jz
   real                     , intent(in)    :: dtlt
   real, dimension(m1,m2,m3), intent(in)    :: tkep
   real, dimension(m1,m2,m3), intent(inout) :: tkepb
   integer                                  :: i, j, k
   real                                     :: timefac

   timefac = dtlt/frqmassave

   do k=1, m1
     do i= ia, iz
       do j= ja, jz
         tkepb(k,i,j)= tkepb(k,i,j) + tkep(k,i,j) * timefac
       end do
     end do
   end do

   return
end subroutine prepare_tke_to_mass
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
! Subroutine prep_turb_to_mass                                                             !
!    Saving the standard deviation of vertical velocity and vertical Lagrangian time scale !
! to the average structures.                                                               !
!------------------------------------------------------------------------------------------!
subroutine prepare_turb_to_mass(m1,m2,m3,ia,iz,ja,jz,dtlt,sigw,tl,sigwb,tlb)
   use mem_mass, only : frqmassave
   implicit none
   integer                  , intent(in)    :: m1,m2,m3
   integer                  , intent(in)    :: ia,iz,ja,jz
   real                     , intent(in)    :: dtlt
   real, dimension(m1,m2,m3), intent(in)    :: sigw,tl
   real, dimension(m1,m2,m3), intent(inout) :: sigwb,tlb
   integer                                  :: i, j, k
   real                                     :: timefac

   timefac = dtlt/frqmassave

   do k=1, m1
      do i= ia, iz
         do j= ja, jz
            sigwb(k,i,j) = sigwb(k,i,j) + sigw(k,i,j) * timefac
            tlb(k,i,j)   = tlb(k,i,j)   + tl(k,i,j)   * timefac
         end do
      end do
   end do

   return
end subroutine prepare_turb_to_mass
!==========================================================================================!
!==========================================================================================!


