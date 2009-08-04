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
           ,scratch%vt3da            , scratch%vt3db           , scratch%vt3dc             &
           ,mass_g(ng)%afxu          , mass_g(ng)%afxv         , mass_g(ng)%afxw           &
           ,mass_g(ng)%afxub         , mass_g(ng)%afxvb        , mass_g(ng)%afxwb          )

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
subroutine prep_convflx_to_mass(m1,mgmzp,maxens_cap,dnmf,upmf,ee,cfxdn,cfxup,dfxdn,dfxup   &
                               ,efxdn,efxup)
   use mem_ensemble     , only : ensemble_vars ! ! Ensemble structure
   use mem_scratch_grell, only : dzd_cld       & ! Delta-z for downdrafts;
                               , kgoff         ! ! BRAMS grid offset
   implicit none
   !---------------------------------------------------------------------------------------!
   !    Arguments.                                                                         !
   !---------------------------------------------------------------------------------------!
   !----- Input variables. ----------------------------------------------------------------!
   integer            , intent(in)    :: m1          ! BRAMS vertical dimension
   integer            , intent(in)    :: mgmzp       ! Grell scheme vertical dimension
   integer            , intent(in)    :: maxens_cap  ! # of static control members
   real               , intent(in)    :: dnmf        ! Reference downdraft mass flux
   real               , intent(in)    :: upmf        ! Reference updraft mass flux
   type(ensemble_vars), intent(in)    :: ee          ! Ensemble structure
   !----- Output variables.  All fluxes are given in kg/m²/s. -----------------------------!
   real, dimension(m1), intent(out)   :: cfxdn       ! Convective downdraft flux
   real, dimension(m1), intent(out)   :: cfxup       ! Convective updraft flux
   real, dimension(m1), intent(out)   :: dfxdn       ! Detrainment associated w/ downdraft
   real, dimension(m1), intent(out)   :: dfxup       ! Detrainment associated w/ updraft
   real, dimension(m1), intent(out)   :: efxdn       ! Entrainment associated w/ downdraft
   real, dimension(m1), intent(out)   :: efxup       ! Entrainment associated w/ updraft
   !---------------------------------------------------------------------------------------!
   !    Local variables.                                                                   !
   !---------------------------------------------------------------------------------------!
   integer                            :: k           ! Grell's grid counter
   integer                            :: kr          ! BRAMS's grid counter
   integer                            :: icap        ! Dynamic control counter
   integer                            :: nmok        ! # of members that were okay.
   integer                            :: jmin        ! Downdraft origin
   integer                            :: k22         ! Updraft origin
   integer                            :: kbcon       ! Level of free convection
   integer                            :: kdet        ! Origin of downdraft detrainment
   integer                            :: ktop        ! Cloud top
   real                               :: nmoki       ! 1./nmok
   real                               :: subin       ! Subsidence from level aloft;
   real                               :: subout      ! Subsidence to level below;
   real                               :: detdo       ! Downdraft detrainment term
   real                               :: detdo1      ! Downdraft detrainm. term @ 1st level
   real                               :: detdo2      ! Downdraft detrainm. term @ 2nd level
   real                               :: entdo       ! Downdraft entrainm. term
   real                               :: detup       ! Updraft detrainment term
   real                               :: entup       ! Updraft entrainment term
   real                               :: totmass     ! Total mass balance
   real, dimension(mgmzp)             :: cdd         ! Norm. downdraft detrainment rate
   real, dimension(mgmzp)             :: cdu         ! Norm. updraft detrainment rate
   real, dimension(mgmzp)             :: mentrd_rate ! Norm. downdraft entrainment rate
   real, dimension(mgmzp)             :: mentru_rate ! Norm. updraft entrainment rate
   real, dimension(mgmzp)             :: dbyd        ! Buoyancy associated with downdrafts
   real, dimension(mgmzp)             :: dbyu        ! Buoyancy associated with updrafts
   real, dimension(mgmzp)             :: etad_cld    ! normalised downdraft mass flux
   real, dimension(mgmzp)             :: etau_cld    ! normalised updraft mass flux
   real, dimension(mgmzp)             :: rhod_cld    ! Downdraft density
   real, dimension(mgmzp)             :: rhou_cld    ! Updraft density
   real, dimension(mgmzp)             :: qliqd_cld   ! Liquid water mix. ratio at downdraft
   real, dimension(mgmzp)             :: qliqu_cld   ! Liquid water mix. ratio at updraft
   real, dimension(mgmzp)             :: qiced_cld   ! Ice mixing ratio at downdraft
   real, dimension(mgmzp)             :: qiceu_cld   ! Ice mixing ratio at updraft
   !----- Intermediate flux variables.  They will be averaged later. ----------------------!
   real, dimension(m1,maxens_cap)     :: cfxdn_cap   ! Convective downdraft flux
   real, dimension(m1,maxens_cap)     :: cfxup_cap   ! Convective updraft flux
   real, dimension(m1,maxens_cap)     :: dfxdn_cap   ! Detrainment associated w/ downdraft
   real, dimension(m1,maxens_cap)     :: dfxup_cap   ! Detrainment associated w/ updraft
   real, dimension(m1,maxens_cap)     :: efxdn_cap   ! Entrainment associated w/ downdraft
   real, dimension(m1,maxens_cap)     :: efxup_cap   ! Entrainment associated w/ updraft
   !---------------------------------------------------------------------------------------!

   nmok      = 0
   cfxdn     = 0.
   cfxup     = 0.
   dfxdn     = 0.
   dfxup     = 0.
   efxdn     = 0.
   efxup     = 0.
   cfxdn_cap = 0.
   cfxup_cap = 0.
   dfxdn_cap = 0.
   dfxup_cap = 0.
   efxdn_cap = 0.
   efxup_cap = 0.

   stacloop: do icap=1,maxens_cap
   
      !----- If this member failed to produce a cloud, skip it. ---------------------------!
      if (ee%ierr_cap(icap) /= 0) cycle stacloop
      
      !----- Otherwise we copy the array values to local variables. -----------------------!
      jmin        = ee%jmin_cap (icap)
      k22         = ee%k22_cap  (icap)
      kbcon       = ee%kbcon_cap(icap)
      kdet        = ee%kdet_cap (icap)
      ktop        = ee%ktop_cap (icap)

      cdd         = ee%cdd_cap(:,icap)
      cdu         = ee%cdu_cap(:,icap)
      mentrd_rate = ee%mentrd_rate_cap(:,icap)
      mentru_rate = ee%mentru_rate_cap(:,icap)
      etad_cld    = ee%etad_cld_cap(:,icap)
      etau_cld    = ee%etau_cld_cap(:,icap)

      !----- We include this member in the final average. ---------------------------------!
      nmok = nmok + 1

      !----- Find the fluxes for this member. ---------------------------------------------!
      vertloop: do k=2,ktop
         kr = k + kgoff ! Output variables should use BRAMS grid.

         !---------------------------------------------------------------------------------!
         ! Computing updraft terms, depending on where I am.                               !
         !---------------------------------------------------------------------------------!
         !----- Below the cloud base, no entrainment or detrainment -----------------------!
         if (k < kbcon .and. k /= k22-1) then 
            entup = 0.
            detup = 0.

         !---------------------------------------------------------------------------------!
         !    Where the updrafts begin, entrainment only. You may ask yourself why only at !
         ! k22-1, and the levels between k22 and kbcon-1 are all zero? This is because the !
         ! net value is zero, since the rates cancel out in this layer.                    !
         !---------------------------------------------------------------------------------!
         elseif (k == k22-1) then
            entup = etau_cld(k22)
            detup = 0.

         !----- In-cloud, both entrainment and detrainment --------------------------------!
         elseif (k >= kbcon .and. k < ktop) then
            entup = mentru_rate(k) * dzd_cld(k) * etau_cld(k)
            detup =       cdu(k+1) * dzd_cld(k) * etau_cld(k)
         
         !----- At the cloud top, detrainment only ----------------------------------------!
         else
            entup = 0.
            subin = 0.       !---- NOTHING enters through the top, not even from above. ---!
            detup = etau_cld(k)
         end if
         !---------------------------------------------------------------------------------!
         
         
         !---------------------------------------------------------------------------------!
         !    Compute downdraft terms, depending on where I am. Note that it's safe to use !
         ! this for shallow clouds, because it kdet and jmin will be both zero, so it will !
         ! always fall in the "nothing happens" case.                                      !
         !---------------------------------------------------------------------------------!
         !----- Below detrainment level, both entrainment and detrainment happen ----------!
         if (k <= kdet) then
            detdo  =         cdd(k) * dzd_cld(k) * etad_cld(k+1)
            entdo  = mentrd_rate(k) * dzd_cld(k) * etad_cld(k+1)
         !----- Within the downdraft layer, but above kdet, only entrainment happens ------!
         elseif (k > kdet .and. k < jmin) then
            detdo  = 0.
            entdo  = mentrd_rate(k) * dzd_cld(k) * etad_cld(k+1)
         !----- Jmin requires special assumption otherwise the entrainment would be zero --!
         elseif (k == jmin) then 
            detdo  = 0.
            entdo  = etad_cld(k)
         !----- Outside the downdraft layer, nothing happens ------------------------------!
         else 
            detdo  = 0.
            entdo  = 0.
         end if
         !---------------------------------------------------------------------------------!
         
         cfxup_cap(kr,icap) =  upmf * etau_cld(k)
         cfxdn_cap(kr,icap) = -dnmf * etad_cld(k)
         dfxup_cap(kr,icap) =  upmf * detup
         efxup_cap(kr,icap) = -upmf * entup
         dfxdn_cap(kr,icap) =  dnmf * detdo
         efxdn_cap(kr,icap) = -dnmf * entdo

      end do vertloop

      !------------------------------------------------------------------------------------!
      ! Bottom layer                                                                       !
      !------------------------------------------------------------------------------------!
      kr= 1 + kgoff  ! the K-level of Grell is equivalent to the BRAMS K+1-level


      detdo1  = etad_cld(2) * cdd(1)         * dzd_cld(1)
      detdo2  = etad_cld(1)
      entdo   = etad_cld(2) * mentrd_rate(1) * dzd_cld(1)

      cfxup_cap(kr,icap) = 0.
      cfxdn_cap(kr,icap) =-dnmf * etad_cld(1)
      dfxup_cap(kr,icap) = 0.
      efxup_cap(kr,icap) = 0.
      dfxdn_cap(kr,icap) = dnmf * (detdo1+detdo2) ! edt already is at detdo1,2
      efxdn_cap(kr,icap) =-dnmf * entdo           ! edt already is at entdo
   end do stacloop
   
   !----- If there was at least one non-zero member, find the average. --------------------!
   if (nmok /= 0) then
      nmoki = 1. / real(nmok)
      do kr=1,m1
         cfxup(kr) = sum(cfxup_cap(kr,:)) * nmoki
         cfxdn(kr) = sum(cfxdn_cap(kr,:)) * nmoki
         dfxup(kr) = sum(dfxup_cap(kr,:)) * nmoki
         dfxdn(kr) = sum(dfxdn_cap(kr,:)) * nmoki
         efxup(kr) = sum(efxup_cap(kr,:)) * nmoki
         efxdn(kr) = sum(efxdn_cap(kr,:)) * nmoki
      end do
   end if

   return
end subroutine prep_convflx_to_mass
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine prepare_timeavg_driver(mzp,mxp,myp,ia,iz,ja,jz,dtlt,ifm,idiffkk)
   use mem_mass , only : imassflx & ! intent(in)
                       , mass_g   ! ! intent(inout)
   use mem_turb , only : turb_g   ! ! intent(in)
   implicit none
   !------ Arguments ----------------------------------------------------------------------!
   integer, intent(in) :: mzp, mxp, myp, ia, iz, ja, jz, ifm, idiffkk
   real   , intent(in) :: dtlt
   !---------------------------------------------------------------------------------------!
   
   !---------------------------------------------------------------------------------------!
   !     Skip this if the user doesn't want to store time averages or if the turbulence    !
   ! closure is not TKE-based.                                                             !
   !---------------------------------------------------------------------------------------!
   if (imassflx /= 1 .or. idiffkk == 2 .or. idiffkk == 3) return

   !----- Checking which closure ----------------------------------------------------------!
   select case (idiffkk)
   case (1) !----- Helfand-Labraga, only TKE and sig-W are available. ---------------------!
      call prepare_timeavg_to_mass(mzp,mxp,myp,ia,iz,ja,jz,dtlt                            &
                                  ,turb_g(ifm)%tkep,mass_g(ifm)%tkepb                      )

      call prepare_timeavg_to_mass(mzp,mxp,myp,ia,iz,ja,jz,dtlt                            &
                                  ,turb_g(ifm)%sigw     ,mass_g(ifm)%sigwb                 )

   case (4,5,6) !----- LES closures, only TKE is available --------------------------------!
      call prepare_timeavg_to_mass(mzp,mxp,myp,ia,iz,ja,jz,dtlt                            &
                                  ,turb_g(ifm)%tkep     ,mass_g(ifm)%tkepb                 )

   case (7) !----- Nakanishi-Niino closure, TKE, sig-W and Lagrangian time scale exist. ---!
      call prepare_timeavg_to_mass(mzp,mxp,myp,ia,iz,ja,jz,dtlt                            &
                                  ,turb_g(ifm)%tkep     ,mass_g(ifm)%tkepb                 )

      call prepare_timeavg_to_mass(mzp,mxp,myp,ia,iz,ja,jz,dtlt                            &
                                  ,turb_g(ifm)%sigw     ,mass_g(ifm)%sigwb                 )

      call prepare_timeavg_to_mass(mzp,mxp,myp,ia,iz,ja,jz,dtlt                            &
                                  ,turb_g(ifm)%ltscale  ,mass_g(ifm)%ltscaleb              )
   end select

   return
end subroutine prepare_timeavg_driver
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
