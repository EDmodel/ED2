!==========================================================================================!
!==========================================================================================!
!    This subroutine computes takes the ensemble averages, integrate them and then compute !
! should be called inside the other loops, so the ensemble variables that should have      !
! the permutation of all ensemble dimensions is fully filled.                              !
!------------------------------------------------------------------------------------------!
subroutine grell_cupar_feedback(comp_down,mgmzp,maxens_cap,maxens_eff,maxens_lsf           &
                               ,maxens_dyn,inv_ensdim,max_heat,dnmf_ens,upmf_ens           &
                               ,dellathil_eff,dellaqtot_eff,dellaco2_eff,pw_eff,ierr_cap   &
                               ,upmf,dnmf,edt,i,j,icld,mynum)

   use rconstants        , only : day_sec       ! ! intent(in)
   use mem_ensemble      , only : ensemble_vars ! ! type
   use mem_scratch_grell , only : mkx           & ! intent(out)
                                , outco2        & ! intent(out)
                                , outqtot       & ! intent(out)
                                , outthil       & ! intent(out)
                                , precip        ! ! intent(out)
   implicit none
   !----- Arguments, input variables. -----------------------------------------------------!
   logical               , intent(in)    :: comp_down         ! I am computing downdrafts
   integer               , intent(in)    :: mgmzp             ! # of levels
   integer               , intent(in)    :: maxens_cap        ! # of static controls
   integer               , intent(in)    :: maxens_eff        ! # of prec. effic. members
   integer               , intent(in)    :: maxens_lsf        ! # of mb members
   integer               , intent(in)    :: maxens_dyn        ! # of dyn. control members
   real                  , intent(in)    :: inv_ensdim        ! 1 / # of members 
   real                  , intent(in)    :: max_heat          ! Maximum heat allowed
   integer               , intent(in)    :: i                 ! X coordinate
   integer               , intent(in)    :: j                 ! Y coordinate
   integer               , intent(in)    :: icld              ! Cloud "coordinate"
   integer               , intent(in)    :: mynum             ! My number, for debugging
   !----- Arguments, input/output variables (ensemble structure). -------------------------!
   real, dimension(maxens_dyn,maxens_lsf,maxens_eff,maxens_cap), intent(inout) ::          &
            dnmf_ens       & ! Reference downdraft mass flux                      [kg/m²/s]
           ,upmf_ens       ! ! Reference updraft mass flux                        [kg/m²/s]
   real   , dimension(mgmzp,maxens_eff,maxens_cap), intent(inout) ::                       &
            dellathil_eff  & ! Change in ice-liquid potential temperature                
           ,dellaqtot_eff  & ! Change in total mixing ratio                              
           ,dellaco2_eff   & ! Change in CO2 mixing ratio                                
           ,pw_eff         ! ! Water that doesn't evaporate (aka rain).                  
   integer, dimension(maxens_cap), intent(inout) ::                                        &
            ierr_cap       ! ! Convection failure flag.
   !----- Arguments, output variables. ----------------------------------------------------!
   real                  , intent(out)   :: upmf              ! Ref. upward mass flx (mb)
   real                  , intent(out)   :: dnmf              ! Ref. dnward mass flx (m0)
   real                  , intent(out)   :: edt               ! m0/mb, Grell's epsilon
   !----- Local variables. ----------------------------------------------------------------!
   integer                               :: k                 ! Counter
   integer                               :: l                 ! Counter
   integer                               :: kmin              ! Minimum location index
   integer                               :: kmax              ! Maximum location index
   integer                               :: iedt              ! efficiency counter
   integer                               :: icap              ! Cap_max counter
   integer                               :: iun
   real                                  :: rescale           ! Rescaling factor
   real                                  :: inv_maxens_effcap ! 1/( maxens_eff*maxens_cap)
   !----- Local constant, controlling debugging information. ------------------------------!
   logical               , parameter     :: print_debug = .false.
   !---------------------------------------------------------------------------------------!

   !----- Assigning the inverse of part of the ensemble dimension. ------------------------!
   inv_maxens_effcap = 1./(maxens_eff*maxens_cap)

   !---------------------------------------------------------------------------------------!
   !    Initialise all output variables.  They may become the actual values in case con-   !
   ! vection didn't happen.                                                                !
   !---------------------------------------------------------------------------------------!
   upmf    = 0.
   dnmf    = 0.
   edt     = 0.
   precip  = 0.
   do k=1,mgmzp
      outthil (k) = 0.
      outqtot (k) = 0.
      outco2  (k) = 0.
   end do
   !---------------------------------------------------------------------------------------!
   !     Before we average, we just need to make sure we don't have negative reference     !
   ! mass fluxes.  If we do, then we flush those to zero.                                  !
   !---------------------------------------------------------------------------------------!
   where (upmf_ens < 0.)
      upmf_ens = 0.
   end where
   
   where (dnmf_ens < 0.)
      dnmf_ens = 0.
   end where

   !---------------------------------------------------------------------------------------!
   !    We will average over all members.  Here we are also doing something slightly       !
   ! different from the original code: upmf (the mb term using Grell's notation) is simply !
   ! the average accross all members, not only the positive ones.  This should avoid bias- !
   ! ing the convection towards the strong convective members.  The argument here is that  !
   ! if most terms are zero, then the environment is unfavourable for convection, so       !
   ! little, if any, convection should happen.                                             !
   !---------------------------------------------------------------------------------------!
   upmf=sum(upmf_ens) * inv_ensdim
   if (comp_down .and. upmf > 0) then
      !------------------------------------------------------------------------------------!
      !     Convection happened and it this cloud supports downdrafts.  Find the downdraft !
      ! reference and the ratio between downdrafts and updrafts.                           !
      !------------------------------------------------------------------------------------!
      dnmf = sum(dnmf_ens) * inv_ensdim
   elseif (upmf == 0.) then
      !------------------------------------------------------------------------------------!
      !     Unlikely, but if the reference upward mass flux is zero, that means that there !
      ! is no cloud...                                                                     !
      !------------------------------------------------------------------------------------!
      where (ierr_cap == 0)
         ierr_cap = 10
      end where
      return
   end if

   !---------------------------------------------------------------------------------------!
   !   Average the temperature tendency among the precipitation efficiency ensemble. If it !
   ! is heating/cooling too much, rescale the reference upward mass flux. Since max_heat   !
   ! is given in K/day, I will compute outt in K/day just to test the value and check      !
   ! whether it is outside the allowed range or not. After I'm done, I will rescale it     !
   ! back to K/s.                                                                          !
   !---------------------------------------------------------------------------------------!
   do k=1,mkx
      do icap = 1, maxens_cap
         do iedt = 1, maxens_eff
            outthil(k) = outthil(k) + dellathil_eff(k,iedt,icap)
         end do
      end do
      outthil(k) = upmf * outthil(k) * inv_maxens_effcap * day_sec
   end do
   !----- Get minimum and maximum outt, and where they happen -----------------------------!
   kmin = minloc(outthil,dim=1)
   kmax = maxloc(outthil,dim=1)
   
   !----- If excessive heat happens, scale down both updrafts and downdrafts --------------!
   if (kmax > 2 .and. outthil(kmax) > max_heat) then
      rescale = max_heat / outthil(kmax)
      upmf = upmf * rescale
      dnmf = dnmf * rescale
      do k=1,mkx
         outthil(k) = outthil(k) * rescale
      end do
   end if
   !----- If excessive cooling happens, scale down both updrafts and downdrafts. ----------!
   if (outthil(kmin)  < - max_heat) then
      rescale = - max_heat/ outthil(kmin)
      upmf = upmf * rescale
      dnmf = dnmf * rescale
      do k=1,mkx
         outthil(k) = outthil(k) * rescale
      end do
   end if
   !----- Heating close to the surface needs to be smaller, being strict there ------------!
   do k=1,2
      if (outthil(k) > 0.5 * max_heat) then
         rescale = 0.5 * max_heat / outthil(k)
         upmf = upmf * rescale
         dnmf = dnmf * rescale
         do l=1,mkx
            outthil(l) = outthil(l) * rescale
         end do
      end if
   end do
   !----- Converting outt to K/s ----------------------------------------------------------!
   do k=1,mkx
      outthil(k) = outthil(k) / day_sec
   end do
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !    With the mass flux and heating on the track, compute other sources/sinks           !
   !---------------------------------------------------------------------------------------!
   do k=1,mkx
      do icap = 1, maxens_cap
         do iedt= 1, maxens_eff
            outqtot(k) = outqtot(k) + dellaqtot_eff(k,iedt,icap)
            outco2 (k) = outco2 (k) + dellaco2_eff (k,iedt,icap)
         end do
      end do
      outqtot(k) = upmf * outqtot(k) * inv_maxens_effcap
      outco2(k)  = upmf * outco2(k)  * inv_maxens_effcap
   end do

   !---------------------------------------------------------------------------------------!
   !   Computing precipitation. It should never be negative, so making sure that this      !
   ! never happens. I will skip this in case this cloud is too shallow.                    !
   !---------------------------------------------------------------------------------------!
   precip = 0.
   if (comp_down) then
      do icap=1,maxens_cap
         do iedt=1,maxens_eff
            precip        = precip + max(0.,upmf*sum(pw_eff(1:mkx,iedt,icap)))
         end do
      end do
      precip = precip * inv_maxens_effcap
   end if
   
   !---------------------------------------------------------------------------------------!
   !    Redefining epsilon.                                                                !
   !---------------------------------------------------------------------------------------!
   if (comp_down .and. upmf > 0) then
      edt  = dnmf/upmf
   end if

   !----- If the user wants print outs, now that is the time... ---------------------------!
   if (print_debug) then
      iun=mynum+50
      write(unit=iun,fmt='(a)') '---------------------------------------------------------'
      write(unit=iun,fmt='(3(a,1x,i5,1x))') ' I=',i,'J=',j,'ICLD=',icld
      write(unit=iun,fmt='(4(a,1x,f10.4,1x))') ' PRECIP   =',precip                        &
                                              ,' EDT      =',edt                           &
                                              ,' DNMF     =',dnmf                          &
                                              ,' UPMF     =',upmf
      write(unit=iun,fmt='(a,3(1x,a))')  'Level','    DTHILDT','    DQTOTDT','     DCO2DT'
      do k=1,mkx
         write(unit=iun,fmt='(i5,3(1x,es12.5))') outthil(k),outqtot(k),outco2(k)
      end do
      write(unit=iun,fmt='(a)') '---------------------------------------------------'
      write(unit=iun,fmt='(a)') ' '
      
   end if

   return
end subroutine grell_cupar_feedback
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine organises the variables that go to the output, namely the heating and !
! moistening rates due to convection.                                                      !
!------------------------------------------------------------------------------------------!
subroutine grell_cupar_output(comp_down,m1,mgmzp,maxens_cap,rtgt,zt,zm,dnmf,upmf,ierr_cap  &
                             ,kdet_cap,k22_cap,kbcon_cap,jmin_cap,ktop_cap,wbuoymin_cap    &
                             ,etad_cld_cap,mentrd_rate_cap,cdd_cap,dbyd_cap,rhod_cld_cap   &
                             ,etau_cld_cap,mentru_rate_cap,cdu_cap,dbyu_cap,rhou_cld_cap   &
                             ,qliqd_cld_cap,qliqu_cld_cap,qiced_cld_cap,qiceu_cld_cap      &
                             ,xierr,zjmin,zk22,zkbcon,zkdt,zktop,conprr,thsrc,rtsrc,co2src &
                             ,areadn,areaup,cuprliq,cuprice)
   use mem_ensemble     , only: &
           ensemble_vars        ! ! type
   use mem_scratch_grell, only: &
           dzd_cld              & ! intent(in) - Top-down layer thickness         [      m]
          ,dzu_cld              & ! intent(in) - Bottom-up layer thickness        [      m]
          ,kgoff                & ! intent(in) - BRAMS grid offset
          ,mkx                  & ! intent(in) - # of cloud grid levels
          ,outco2               & ! intent(in) - Total CO2 mixing ratio forcing   [  ppm/s]
          ,outqtot              & ! intent(in) - Total mixing ratio forcing       [kg/kg/s]
          ,outthil              & ! intent(in) - Theta_il forcing                 [    K/s]
          ,precip               & ! intent(in) - Precipitation rate               [kg/m²/s]
          ,sigw                 & ! intent(in) - Vertical velocity std. deviation [    m/s]
          ,tke                  & ! intent(in) - Turbulent kinetic energy         [   J/kg]
          ,wwind                ! ! intent(in) - Vertical velocity                [    m/s]

   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   logical            , intent(in)  :: comp_down   ! I will compute downdraft stuff
   integer            , intent(in)  :: m1          ! Number of levels
   integer            , intent(in)  :: mgmzp       ! Number of Grell's levels
   integer            , intent(in)  :: maxens_cap  ! Number of static control realisations
   real               , intent(in)  :: rtgt        ! Corr. to get the heights     [   ----]
   real, dimension(m1), intent(in)  :: zt          ! Height at thermodyn. levels  [      m]
   real, dimension(m1), intent(in)  :: zm          ! Height at momentum levels    [      m]
   real               , intent(in)  :: dnmf        ! Reference downdraft mass flux[kg/m²/s]
   real               , intent(in)  :: upmf        ! Reference updraft mass flux  [kg/m²/s]
   !----- Input, ensemble structure variables. --------------------------------------------!
   integer, dimension      (maxens_cap), intent(in) :: ierr_cap
   integer, dimension      (maxens_cap), intent(in) :: kdet_cap
   integer, dimension      (maxens_cap), intent(in) :: k22_cap
   integer, dimension      (maxens_cap), intent(in) :: kbcon_cap
   integer, dimension      (maxens_cap), intent(in) :: jmin_cap
   integer, dimension      (maxens_cap), intent(in) :: ktop_cap
   real   , dimension      (maxens_cap), intent(in) :: wbuoymin_cap
   real   , dimension(mgmzp,maxens_cap), intent(in) :: etad_cld_cap
   real   , dimension(mgmzp,maxens_cap), intent(in) :: mentrd_rate_cap
   real   , dimension(mgmzp,maxens_cap), intent(in) :: cdd_cap
   real   , dimension(mgmzp,maxens_cap), intent(in) :: dbyd_cap
   real   , dimension(mgmzp,maxens_cap), intent(in) :: rhod_cld_cap
   real   , dimension(mgmzp,maxens_cap), intent(in) :: etau_cld_cap
   real   , dimension(mgmzp,maxens_cap), intent(in) :: mentru_rate_cap
   real   , dimension(mgmzp,maxens_cap), intent(in) :: cdu_cap
   real   , dimension(mgmzp,maxens_cap), intent(in) :: dbyu_cap
   real   , dimension(mgmzp,maxens_cap), intent(in) :: rhou_cld_cap
   real   , dimension(mgmzp,maxens_cap), intent(in) :: qliqd_cld_cap
   real   , dimension(mgmzp,maxens_cap), intent(in) :: qliqu_cld_cap
   real   , dimension(mgmzp,maxens_cap), intent(in) :: qiced_cld_cap
   real   , dimension(mgmzp,maxens_cap), intent(in) :: qiceu_cld_cap
   !----- Output variables. ---------------------------------------------------------------!
   real, dimension(m1), intent(out) :: thsrc       ! Potential temperature fdbck  [    K/s]
   real, dimension(m1), intent(out) :: rtsrc       ! Total mixing ratio feedback  [kg/kg/s]
   real, dimension(m1), intent(out) :: co2src      ! Total CO2 mixing ratio fdbk  [  ppm/s]
   real, dimension(m1), intent(out) :: cuprliq     ! Cumulus water mixing ratio   [  kg/kg]
   real, dimension(m1), intent(out) :: cuprice     ! Cumulus ice mixing ratio     [  kg/kg]
   real               , intent(out) :: areadn      ! Fractional downdraft area    [    ---]
   real               , intent(out) :: areaup      ! Fractional updraft area      [    ---]
   real               , intent(out) :: conprr      ! Rate of convective precip.   [kg/m²/s]
   real               , intent(out) :: xierr       ! Error flag
   real               , intent(out) :: zjmin       ! Downdraft originating level  [      m]
   real               , intent(out) :: zk22        ! Updraft origin               [      m]
   real               , intent(out) :: zkbcon      ! Level of free convection     [      m]
   real               , intent(out) :: zkdt        ! Top of the dndraft detr.     [      m]
   real               , intent(out) :: zktop       ! Cloud top                    [      m]
   !----- Local variables. ----------------------------------------------------------------!
   integer                          :: icap        ! Static control counter.
   integer                          :: nmok        ! # of members that had clouds.
   integer                          :: k           ! Cloud level counter
   integer                          :: kr          ! BRAMS level counter
   real                             :: exner       ! Exner fctn. for tend. conv.  [ J/kg/K]
   real                             :: nmoki       ! 1/nmok
   integer                          :: jmin        ! Downdraft origin
   integer                          :: k22         ! Updraft origin
   integer                          :: kbcon       ! Level of free convection
   integer                          :: kdet        ! Origin of downdraft detrainment
   integer                          :: ktop        ! Cloud top
   !----- Aux. variables for fractional area. ---------------------------------------------! 
   real, dimension(maxens_cap)      :: areadn_cap  ! Fractional downdraft area    [    ---]
   real, dimension(maxens_cap)      :: areaup_cap  ! Fractional updraft area      [    ---]
   real, dimension(m1,maxens_cap)   :: cuprliq_cap ! Cumulus water mixing ratio   [  kg/kg]
   real, dimension(m1,maxens_cap)   :: cuprice_cap ! Cumulus ice mixing ratio     [  kg/kg]
   
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !    Flushing all variables to zero in case convection didn't happen.                   !
   !---------------------------------------------------------------------------------------!
   do k=1,m1
      thsrc  (k) = 0.
      rtsrc  (k) = 0.
      co2src (k) = 0.
      cuprliq(k) = 0.
      cuprice(k) = 0.
   end do

   areadn  = 0.
   areaup  = 0.
   conprr  = 0.
   zkdt    = 0.
   zk22    = 0.
   zkbcon  = 0.
   zjmin   = 0.
   zktop   = 0.
   
   !---------------------------------------------------------------------------------------!
   !   Copying the error flag. This should not be zero in case convection failed.          !
   !---------------------------------------------------------------------------------------!
   if (upmf > 0. .and. any(ierr_cap == 0)) then
      xierr = 0.
   else
      xierr = real(ierr_cap(1))
   end if


   if (xierr /= 0.) return !----- No cloud, no need to return anything --------------------!

   !---------------------------------------------------------------------------------------!
   !     We choose the level to be the most level that will encompass all non-zero         !
   ! members.                                                                              !
   !---------------------------------------------------------------------------------------!
   kdet  = maxval(kdet_cap ,mask = ierr_cap == 0)
   k22   = minval(k22_cap  ,mask = ierr_cap == 0)
   kbcon = minval(kbcon_cap,mask = ierr_cap == 0)
   jmin  = maxval(jmin_cap ,mask = ierr_cap == 0)
   ktop  = maxval(ktop_cap ,mask = ierr_cap == 0)

   !---------------------------------------------------------------------------------------!
   !    Fixing the levels, here I will add back the offset so the output will be consist-  !
   ! ent. I will return these variables even when no cloud developed for debugging         !
   ! purposes. When the code is running fine, then I should return them only when          !
   ! convection happens.                                                                   !
   !---------------------------------------------------------------------------------------!
   zkdt   = (zt(kdet  + kgoff)-zm(kgoff))*rtgt
   zk22   = (zt(k22   + kgoff)-zm(kgoff))*rtgt
   zkbcon = (zt(kbcon + kgoff)-zm(kgoff))*rtgt
   zjmin  = (zt(jmin  + kgoff)-zm(kgoff))*rtgt
   zktop  = (zt(ktop  + kgoff)-zm(kgoff))*rtgt

   !---------------------------------------------------------------------------------------!
   !    Precipitation is simply copied, it could even be output directly from the main     !
   ! subroutine, brought here just to be together with the other source terms.             !
   !---------------------------------------------------------------------------------------!
   conprr = precip
   
   do k=1,mkx
      kr    = k + kgoff
      !----- Here we are simply copying including the offset back. ------------------------!
      thsrc(kr)   = outthil(k)
      rtsrc(kr)   = outqtot(k)
      co2src(kr)  = outco2(k)
   end do


   !---------------------------------------------------------------------------------------!
   !   Computing the relative area covered by downdrafts and updrafts.                     !
   !---------------------------------------------------------------------------------------!
   nmok = 0
   stacloop: do icap=1,maxens_cap
      if (ierr_cap(icap) /= 0) cycle stacloop
      nmok = nmok + 1
      call grell_draft_area(comp_down,m1,mgmzp,kgoff,jmin_cap(icap),k22_cap(icap)          &
                           ,kbcon_cap(icap),ktop_cap(icap),dzu_cld,wwind,tke,sigw          &
                           ,wbuoymin_cap(icap),etad_cld_cap(1:mgmzp,icap)                  &
                           ,mentrd_rate_cap(1:mgmzp,icap),cdd_cap(1:mgmzp,icap)            &
                           ,dbyd_cap(1:mgmzp,icap),rhod_cld_cap(1:mgmzp,icap),dnmf         &
                           ,etau_cld_cap(1:mgmzp,icap),mentru_rate_cap(1:mgmzp,icap)       &
                           ,cdu_cap(1:mgmzp,icap),dbyu_cap(1:mgmzp,icap)                   &
                           ,rhou_cld_cap(1:mgmzp,icap),upmf,areadn_cap(icap)               &
                           ,areaup_cap(icap))

      !------------------------------------------------------------------------------------!
      !   I compute the cloud condensed mixing ratio for this realisation.  This is used   !
      ! by Harrington when cumulus feedback is requested, so we rescale the liquid water   !
      ! at the downdrafts and updrafts by their area.                                      !
      !------------------------------------------------------------------------------------!
      do k=1,ktop_cap(icap)
         kr=k+kgoff
         cuprliq_cap(kr,icap) = max(0., qliqd_cld_cap(k,icap) * areadn_cap(icap)           &
                              + qliqu_cld_cap(k,icap) * areaup_cap(icap) )
         cuprice_cap(kr,icap) = max(0., qiced_cld_cap(k,icap) * areadn_cap(icap)           &
                              + qiceu_cld_cap(k,icap) * areaup_cap(icap) )
      end do
   end do stacloop
   !----- Find the averaged area. ---------------------------------------------------------! 
   if (nmok /= 0) then
      nmoki = 1./real(nmok)
      areadn = sum(areadn_cap) * nmoki
      areaup = sum(areaup_cap) * nmoki
      do kr=1,m1
         cuprliq(kr) = sum(cuprliq_cap(kr,1:maxens_cap)) * nmoki
         cuprice(kr) = sum(cuprice_cap(kr,1:maxens_cap)) * nmoki
      end do
   end if


   return
end subroutine grell_cupar_output
!==========================================================================================!
!==========================================================================================!
