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






!==========================================================================================!
!==========================================================================================!
! Subroutine save_convective_mass                                                          !
!                                                                                          !
!   This subroutine prepares the convective fluxes to be used in Lagrangian models.        !
!------------------------------------------------------------------------------------------!
subroutine prep_convflx_to_mass(m1,dnmf,upmf,cfxdn,cfxup,dfxdn,dfxup,efxdn,efxup)

   use mem_scratch_grell, only: &
          cdd                   & ! Downdraft detrainment function;
         ,cdu                   & ! Updraft detrainment function;
         ,dzd_cld               & ! Delta-z for downdrafts;
         ,etad_cld              & ! Normalized dndraft mass flux;
         ,etau_cld              & ! Normalized updraft mass flux;
         ,ierr                  & ! Error flag;
         ,jmin                  & ! Downdraft origin level
         ,k22                   & ! Updraft origin level
         ,kbcon                 & ! Cloud base level
         ,kdet                  & ! Detrainment level
         ,kgoff                 & ! BRAMS grid offset
         ,ktop                  & ! Cloud top level
         ,mentrd_rate           & ! Downdraft entrainment rate
         ,mentru_rate           ! ! Updraft entrainment rate

   implicit none
   integer            , intent(in)    :: m1      ! BRAMS vertical dimension
   real               , intent(in)    :: dnmf    ! Reference downdraft mass flux   [kg/m²/s]
   real               , intent(in)    :: upmf    ! Reference updraft mass flux     [kg/m²/s]

   !---------------------------------------------------------------------------------------!
   !    All fluxes are given in kg/m²/s.                                                   !
   !---------------------------------------------------------------------------------------!
   !----- Deep convection variables -------------------------------------------------------!
   real, dimension(m1), intent(out) :: cfxdn  ! Convective downdraft flux
   real, dimension(m1), intent(out) :: cfxup  ! Convective updraft flux
   real, dimension(m1), intent(out) :: dfxdn  ! Detrainment associated with downdraft
   real, dimension(m1), intent(out) :: dfxup  ! Detrainment associated with updraft
   real, dimension(m1), intent(out) :: efxdn  ! Entrainment associated with downdraft
   real, dimension(m1), intent(out) :: efxup  ! Entrainment associated with updraft
   !---------------------------------------------------------------------------------------!

   integer                            :: k       ! Grell's grid counter
   integer                            :: kr      ! BRAMS's grid counter
   real                               :: subin   ! Subsidence from level aloft;
   real                               :: subout  ! Subsidence to level below;
   real                               :: detdo   ! Downdraft detrainment term
   real                               :: detdo1  ! Downdraft detrainment term @ 1st level
   real                               :: detdo2  ! Downdraft detrainment term @ 2nd level
   real                               :: entdo   ! Downdraft entrainment term
   real                               :: detup   ! Updraft detrainment term
   real                               :: entup   ! Updraft entrainment term
   real                               :: totmass ! Total mass balance

   cfxdn = 0.
   cfxup = 0.
   dfxdn = 0.
   dfxup = 0.
   efxdn = 0.
   efxup = 0.

   vertloop: do k=2,ktop
      kr = k + kgoff ! Output variables should use BRAMS grid.

      !------------------------------------------------------------------------------------!
      ! Computing updraft terms, depending on where I am.                                  !
      !------------------------------------------------------------------------------------!
      !----- Below the cloud base, no entrainment or detrainment --------------------------!
      if (k < kbcon .and. k /= k22-1) then 
         entup = 0.
         detup = 0.

      !------------------------------------------------------------------------------------!
      !    Where the updrafts begin, entrainment only. You may ask yourself why only at    !
      ! k22-1, and the levels between k22 and kbcon-1 are all zero? This is because the    !
      ! net value is zero, since the rates cancel out in this layer.                       !
      !------------------------------------------------------------------------------------!
      elseif (k == k22-1) then
         entup = etau_cld(k22)
         detup = 0.

      !----- In-cloud, both entrainment and detrainment -----------------------------------!
      elseif (k >= kbcon .and. k < ktop) then
         entup = mentru_rate(k) * dzd_cld(k) * etau_cld(k)
         detup =       cdu(k+1) * dzd_cld(k) * etau_cld(k)
      
      !----- At the cloud top, detrainment only -------------------------------------------!
      else
         entup = 0.
         subin = 0.          !---- NOTHING enters through the top, not even from above. ---!
         detup = etau_cld(k)
      end if
      !------------------------------------------------------------------------------------!
      
      
      !------------------------------------------------------------------------------------!
      !    Compute downdraft terms, depending on where I am. Note that it's safe to use    !
      ! this for shallow clouds, because it kdet and jmin will be both zero, so it will    !
      ! always fall in the "nothing happens" case.                                         !
      !------------------------------------------------------------------------------------!
      !----- Below detrainment level, both entrainment and detrainment happen -------------!
      if (k <= kdet) then
         detdo  = cdd(k) * dzd_cld(k) * etad_cld(k+1)
         entdo  = mentrd_rate(k) * dzd_cld(k) * etad_cld(k+1)
      !----- Within the downdraft layer, but above kdet, only entrainment happens ---------!
      elseif (k > kdet .and. k < jmin) then
         detdo  = 0.
         entdo  = mentrd_rate(k) * dzd_cld(k) * etad_cld(k+1)
      !----- Jmin requires special assumption otherwise the entrainment would be zero -----!
      elseif (k == jmin) then 
         detdo  = 0.
         entdo  = etad_cld(k)
      !----- Outside the downdraft layer, nothing happens ---------------------------------!
      else 
         detdo  = 0.
         entdo  = 0.
      end if
      !------------------------------------------------------------------------------------!
      
      cfxup(kr) =  upmf * etau_cld(k)
      cfxdn(kr) = -dnmf * etad_cld(k)
      dfxup(kr) =  upmf*detup
      efxup(kr) = -upmf*entup
      dfxdn(kr) =  dnmf*detdo
      efxdn(kr) = -dnmf*entdo

   end do vertloop

   !---------------------------------------------------------------------------------------!
   ! Bottom layer                                                                          !
   !---------------------------------------------------------------------------------------!
   kr= 1 + kgoff  ! the K-level of Grell is equivalent to the BRAMS K+1-level


   detdo1  = etad_cld(2)*cdd(1)*dzd_cld(1)
   detdo2  = etad_cld(1)
   entdo   = etad_cld(2)*mentrd_rate(1)*dzd_cld(1)

   cfxup(kr) = 0.
   cfxdn(kr) =-dnmf * etad_cld(1)
   dfxup(kr) = 0.
   efxup(kr) = 0.
   dfxdn(kr) = dnmf*(detdo1+detdo2) !edt already is at detdo1,2
   efxdn(kr) =-dnmf*entdo           !edt already is at entdo
   return
end subroutine prep_convflx_to_mass
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
! Subroutine prepare_timeavg_to_mass                                                       !
!                                                                                          !
!   The aim of this subroutine is simply save the mean value of any variable at the        !
! regular and lite analysis.                                                               !
!------------------------------------------------------------------------------------------!
subroutine prepare_timeavg_to_mass(m1,m2,m3,ia,iz,ja,jz,dtlt,var,avgvar)
   use mem_mass, only : frqmassave
   implicit none
   integer, intent(in)                      :: m1,m2,m3
   integer, intent(in)                      :: ia,iz,ja,jz
   real                     , intent(in)    :: dtlt
   real, dimension(m1,m2,m3), intent(in)    :: var
   real, dimension(m1,m2,m3), intent(inout) :: avgvar
   integer                                  :: i, j, k
   real                                     :: timefac

   timefac = dtlt/frqmassave

   do k=1, m1
     do i= ia, iz
       do j= ja, jz
         avgvar(k,i,j)= avgvar(k,i,j) + var(k,i,j) * timefac
       end do
     end do
   end do

   return
end subroutine prepare_timeavg_to_mass
!==========================================================================================!
!==========================================================================================!
