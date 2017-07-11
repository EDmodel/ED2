!==========================================================================================!
!==========================================================================================!
!    This subroutine computes takes the ensemble averages, integrate them and then compute !
! should be called inside the other loops, so the ensemble variables that should have      !
! the permutation of all ensemble dimensions is fully filled.                              !
!------------------------------------------------------------------------------------------!
subroutine grell_cupar_feedback(mgmzp,maxens_cap,maxens_eff,maxens_lsf,maxens_dyn          &
                               ,inv_ensdim,max_heat,dnmf_ens,upmf_ens,dnmx_ens             &
                               ,upmx_ens,dnmf_cap,upmf_cap,dellathil_eff,dellaqtot_eff     &
                               ,dellaco2_eff,pw_eff,outco2,outqtot,outthil,precip,ierr_cap &
                               ,comp_down_cap,dnmf,upmf,dnmx,upmx,edt,i,j,icld,mynum)

   use rconstants        , only : day_sec       ! ! intent(in)
   use mem_ensemble      , only : ensemble_vars ! ! type
   use mem_scratch_grell , only : mkx           & ! intent(in)
                                , z_cup         ! ! intent(in)
   implicit none
   !----- Arguments, input variables. -----------------------------------------------------!
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
           ,upmf_ens       & ! Reference updraft mass flux                        [kg/m²/s]
           ,dnmx_ens       & ! Potential downdraft mass flux                      [kg/m²/s]
           ,upmx_ens       ! ! Potential updraft mass flux                        [kg/m²/s]
   real   , dimension(mgmzp,maxens_eff,maxens_cap), intent(in) ::                          &
            dellathil_eff  & ! Change in ice-liquid potential temperature                
           ,dellaqtot_eff  & ! Change in total mixing ratio                              
           ,dellaco2_eff   & ! Change in CO2 mixing ratio                                
           ,pw_eff         ! ! Water that doesn't evaporate (aka rain).                  
   integer, dimension(maxens_cap), intent(inout) ::                                        &
            ierr_cap       ! ! Convection failure flag.                           [    ---]
   logical, dimension(maxens_cap), intent(inout) ::                                        &
            comp_down_cap  ! ! Computed downdraft for this stat. control ensemble.[    ---]
   real   , dimension(maxens_cap), intent(inout) ::                                        &
            dnmf_cap       & ! Reference downdraft mass flux                      [kg/m²/s]
           ,upmf_cap       ! ! Reference updraft mass flux                        [kg/m²/s]
   real   , dimension(mgmzp)     , intent(inout) ::                                        &
            outco2         & ! CO2 mixing ratio tendency                         [   ppm/s]
           ,outqtot        & ! Total water mixing ratio tendency                 [ kg/kg/s]
           ,outthil        ! ! Ice-liquid pot. temperature tendency              [     K/s]
   real   , dimension(1)         , intent(inout) ::                                        &
            precip         ! ! Precipitation rate                                 [kg/m²/s]
   !----- Arguments, output variables. ----------------------------------------------------!
   real                  , intent(inout) :: dnmf              ! Ref. dnward mass flx (m0)
   real                  , intent(inout) :: upmf              ! Ref. upward mass flx (mb)
   real                  , intent(inout) :: dnmx              ! Pot. dnward mass flx (m0p)
   real                  , intent(inout) :: upmx              ! Pot. upward mass flx (mbp)
   real                  , intent(inout) :: edt               ! m0/mb, Grell's epsilon
   !----- Local variables. ----------------------------------------------------------------!
   integer                               :: k              ! Counter
   integer                               :: l              ! Counter
   integer                               :: kmin           ! Minimum location index
   integer                               :: kmax           ! Maximum location index
   integer                               :: iedt           ! efficiency counter
   integer                               :: idyn           ! dynamic control counter
   integer                               :: icap           ! Cap_max counter
   integer                               :: imbp           ! Large-scale forcing counter
   integer                               :: iun            ! Unit
   real                                  :: rescale        ! Rescaling factor
   real                                  :: inv_maxens_ec  ! 1 / ( maxens_eff * maxens_cap)
   real                                  :: inv_maxens_eld ! 1 / ( maxens_eff * maxens_cap &
                                                           !     * maxens_dyn )
   real                                  :: max_heat_si    ! Maximum heating rate in K/s
   !----- Local constant, controlling debugging information. ------------------------------!
   logical               , parameter     :: print_debug = .false.
   !---------------------------------------------------------------------------------------!



   !----- Assign the inverse of part of the ensemble dimension. ---------------------------!
   inv_maxens_ec  = 1. / (maxens_eff * maxens_cap )
   inv_maxens_eld = 1. / (maxens_eff * maxens_lsf * maxens_dyn)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Find the maximum heating rate in K/s.                                             !
   !---------------------------------------------------------------------------------------!
   max_heat_si = max_heat / day_sec
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    Initialise all output variables.  They may become the actual values in case con-   !
   ! vection didn't happen.                                                                !
   !---------------------------------------------------------------------------------------!
   dnmf     = 0.
   upmf     = 0.
   dnmx     = 0.
   upmx     = 0.
   edt      = 0.
   precip   = 0.
   outthil  = 0.
   outqtot  = 0.
   outco2   = 0.
   !---------------------------------------------------------------------------------------!


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

   where (upmx_ens < 0.)
      upmx_ens = 0.
   end where
   
   where (dnmx_ens < 0.)
      dnmx_ens = 0.
   end where
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Find the averaged mass fluxes for each static control.  We will average over all  !
   ! members.  Here we are also doing something slightly different from the original code: !
   ! upmf (the mb term using Grell's notation) is simply the average accross all members,  !
   ! not only the positive ones.  This should avoid biasing the convection towards the     !
   ! strong convective members.  The argument here is that if most terms are zero, then    !
   ! the environment is unfavourable for convection, so little, if any, convection should  !
   ! happen.                                                                               !
   !---------------------------------------------------------------------------------------!
   dnmf = 0.
   dnmx = 0.
   upmf = 0.
   upmx = 0.
   do icap=1,maxens_cap

      dnmf_cap(icap) = 0.
      upmf_cap(icap) = 0.

      do iedt=1,maxens_eff
         do imbp=1,maxens_lsf
            do idyn=1,maxens_dyn
               dnmf_cap(icap) = dnmf_cap(icap) + max(0.,dnmf_ens(idyn,imbp,iedt,icap))
               upmf_cap(icap) = upmf_cap(icap) + max(0.,upmf_ens(idyn,imbp,iedt,icap))
               dnmf           = dnmf           + max(0.,dnmf_ens(idyn,imbp,iedt,icap))
               dnmx           = dnmx           + max(0.,dnmx_ens(idyn,imbp,iedt,icap))
               upmf           = upmf           + max(0.,upmf_ens(idyn,imbp,iedt,icap))
               upmx           = upmx           + max(0.,upmx_ens(idyn,imbp,iedt,icap))
            end do
         end do
      end do

      !----- Normalise the mass fluxes. ---------------------------------------------------!
      dnmf_cap(icap) = dnmf_cap(icap) * inv_maxens_eld
      upmf_cap(icap) = upmf_cap(icap) * inv_maxens_eld
      !------------------------------------------------------------------------------------!

   end do
   !---------------------------------------------------------------------------------------!



   !----- Normalise the mass fluxes. ------------------------------------------------------!
   dnmf = dnmf * inv_ensdim
   dnmx = dnmx * inv_ensdim
   upmf = upmf * inv_ensdim
   upmx = upmx * inv_ensdim
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   if (upmf == 0. .and. any(ierr_cap == 0)) then
      !------------------------------------------------------------------------------------!
      !     Unlikely, but if the reference upward mass flux is zero, that means that there !
      ! is no cloud...                                                                     !
      !------------------------------------------------------------------------------------!


      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
      if (print_debug) then
         iun = 40 + icld
         write (unit=iun,fmt='(143a)') ('=',k=1,143)
         write(unit=iun,fmt='(1(a,1x,i12,1x))'  ) ' I        =',i
         write(unit=iun,fmt='(1(a,1x,i12,1x))'  ) ' J        =',j
         write(unit=iun,fmt='(1(a,1x,f12.5,1x))') ' DNMF     =',dnmf
         write(unit=iun,fmt='(1(a,1x,f12.5,1x))') ' DNMX     =',dnmx
         write(unit=iun,fmt='(1(a,1x,f12.5,1x))') ' UPMF     =',upmf
         write(unit=iun,fmt='(1(a,1x,f12.5,1x))') ' UPMX     =',upmx
         write (unit=iun,fmt='(143a)') ('-',k=1,143)
         write (unit=iun,fmt='(11(a,1x))') '        ICAP','        IEDT','        IMBP'    &
                                          ,'        IDYN','    IERR_CAP','    DNMF_ENS'    &
                                          ,'    DNMX_ENS','    DNMF_CAP','    UPMF_ENS'    &
                                          ,'    UPMX_ENS','    DNMF_CAP'
         write (unit=iun,fmt='(143a)') ('-',k=1,143)
         do icap=1,maxens_cap
            do iedt=1,maxens_eff
               do imbp=1,maxens_lsf
                  do idyn=1,maxens_dyn
                     write (unit=iun,fmt='(5(i12,1x),6(f12.5,1x))')                        &
                                                      icap,iedt,imbp,idyn,ierr_cap(icap)   &
                                                     ,dnmf_ens(idyn,imbp,iedt,icap)        &
                                                     ,dnmx_ens(idyn,imbp,iedt,icap)        &
                                                     ,dnmf_cap(icap)                       &
                                                     ,upmf_ens(idyn,imbp,iedt,icap)        &
                                                     ,upmx_ens(idyn,imbp,iedt,icap)        &
                                                     ,upmf_cap(icap)
                  end do
               end do
            end do
         end do
         write (unit=iun,fmt='(143a)') ('=',k=1,143)
         write (unit=iun,fmt='(a)'   ) ' '
      end if
      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!


      where (ierr_cap == 0)
         ierr_cap = 13
      end where
      return
   end if
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !   Average the temperature tendency among the precipitation efficiency ensemble. If it !
   ! is heating/cooling too much, rescale the reference upward mass flux.                  !
   !---------------------------------------------------------------------------------------!
   do k=1,mkx
      do icap=1,maxens_cap
         do iedt=1,maxens_eff
            outthil(k) = outthil(k) + upmf * dellathil_eff(k,iedt,icap)
         end do
      end do
      outthil(k) = outthil(k) * inv_maxens_ec
   end do
   !----- Get minimum and maximum outt, and where they happen -----------------------------!
   kmin = minloc(outthil,dim=1)
   kmax = maxloc(outthil,dim=1)
   
   !----- If excessive heat happens, scale down both updrafts and downdrafts --------------!
   if (kmax > 2 .and. outthil(kmax) > max_heat_si) then
      rescale = max_heat_si / outthil(kmax)
      upmf        = upmf        * rescale
      dnmf        = dnmf        * rescale
      upmx        = upmx        * rescale
      dnmx        = dnmx        * rescale
      dnmf_cap(:) = dnmf_cap(:) * rescale
      upmf_cap(:) = upmf_cap(:) * rescale
      do k=1,mkx
         outthil(k) = outthil(k) * rescale
      end do
   end if
   !----- If excessive cooling happens, scale down both updrafts and downdrafts. ----------!
   if (outthil(kmin)  < - max_heat_si) then
      rescale = - max_heat_si/ outthil(kmin)
      upmf        = upmf        * rescale
      dnmf        = dnmf        * rescale
      upmx        = upmx        * rescale
      dnmx        = dnmx        * rescale
      dnmf_cap(:) = dnmf_cap(:) * rescale
      upmf_cap(:) = upmf_cap(:) * rescale
      do k=1,mkx
         outthil(k) = outthil(k) * rescale
      end do
   end if
   !----- Heating close to the surface needs to be smaller, being strict there ------------!
   do k=1,2
      if (outthil(k) > 0.5 * max_heat_si) then
         rescale = 0.5 * max_heat_si / outthil(k)
         upmf        = upmf        * rescale
         dnmf        = dnmf        * rescale
         upmx        = upmx        * rescale
         dnmx        = dnmx        * rescale
         dnmf_cap(:) = dnmf_cap(:) * rescale
         upmf_cap(:) = upmf_cap(:) * rescale
         do l=1,mkx
            outthil(l) = outthil(l) * rescale
         end do
      end if
   end do
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !    With the mass flux and heating on the track, compute other sources/sinks           !
   !---------------------------------------------------------------------------------------!
   do k=1,mkx
      do icap=1,maxens_cap
         do iedt=1,maxens_eff
            outqtot(k) = outqtot(k) + upmf * dellaqtot_eff(k,iedt,icap)
            outco2 (k) = outco2 (k) + upmf * dellaco2_eff (k,iedt,icap)
         end do
      end do
      outqtot(k) = outqtot(k) * inv_maxens_ec
      outco2 (k) = outco2 (k) * inv_maxens_ec
   end do
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !   Compute precipitation.  It should never be negative, so we check whether this ever  !
   ! happens.                                                                              !
   !---------------------------------------------------------------------------------------!
   if (any(comp_down_cap)) then
      do icap=1,maxens_cap
         do iedt=1,maxens_eff
            do k=1,mkx
               precip(1) = precip(1) + upmf * pw_eff(k,iedt,icap)
            end do
         end do
      end do
      precip(1) = max(0.,precip(1)) * inv_maxens_ec
   end if
   
   !---------------------------------------------------------------------------------------!
   !    Redefine epsilon.                                                                  !
   !---------------------------------------------------------------------------------------!
   if (any(comp_down_cap) .and. upmf > 0.) then
      edt  = dnmf/upmf
   end if
   !---------------------------------------------------------------------------------------!




   !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
   !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
   !      If the user wants print outs, now that is the time...                            !
   !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
   if (print_debug .and. upmf /= 0.) then
      iun=mynum+50
      write(unit=iun,fmt='(70a)'       ) ('=',k=1,70)
      write(unit=iun,fmt='(a,1x,i5)'   ) ' I        =',i
      write(unit=iun,fmt='(a,1x,i5)'   ) ' J        =',j
      write(unit=iun,fmt='(a,1x,i5)'   ) ' ICLD     =',icld
      write(unit=iun,fmt='(a,1x,f10.4)') ' PRECIP   =',precip(1)
      write(unit=iun,fmt='(a,1x,f10.4)') ' EDT      =',edt
      write(unit=iun,fmt='(a,1x,f10.4)') ' DNMF     =',dnmf
      write(unit=iun,fmt='(a,1x,f10.4)') ' UPMF     =',upmf
      write(unit=iun,fmt='(70a)'       ) ('-',k=1,70)
      write(unit=iun,fmt='(5(a,1x))'   )  '       LEVEL','      HEIGHT'                    &
                                         ,'     DTHILDT','     DQTOTDT','      DCO2DT'
      write(unit=iun,fmt='(5(a,1x))'   )  '            ','           m'                    &
                                         ,'       K/day','    g/kg/day','umol/mol/day'
      write(unit=iun,fmt='(70a)'       ) ('-',k=1,70)
      do k=mkx,1,-1
         write (unit=iun,fmt='(i13,1x,4(f13.5,1x))') k,z_cup  (k)                          &
                                                      ,outthil(k) * day_sec                &
                                                      ,outqtot(k) * 1000. * day_sec        &
                                                      ,outco2 (k) * day_sec
      end do
      write(unit=iun,fmt='(70a)'       ) ('=',k=1,70)
      write(unit=iun,fmt='(a)') ' '
   end if
   !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
   !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!

   return
end subroutine grell_cupar_feedback
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine finds the area covered by each cloud.  Because overlap of potential  !
! areas may exist, we must scale the actual area by the maximum area, using the actual     !
! areas as weighting factors.                                                              !
!------------------------------------------------------------------------------------------!
subroutine grell_cupar_area_scaler(cldd,clds,m1,mgmzp,maxens_cap)
   use mem_ensemble     , only : &
           ensemble_e          ! ! structure
   use mem_scratch_grell, only : &
           dzd_cld             & ! intent(in) - Top-down layer thickness          [      m]
          ,dzu_cld             & ! intent(in) - Bottom-up layer thickness         [      m]
          ,kgoff               & ! intent(in) - BRAMS grid offset               
          ,mkx                 & ! intent(in) - # of cloud grid levels          
          ,sigw                & ! intent(in) - Vertical velocity std. dev.       [    m/s]
          ,tke                 & ! intent(in) - Turbulent kinetic energy          [   J/kg]
          ,wwind               ! ! intent(in) - Vertical velocity                 [    m/s]
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer            , intent(in) :: m1         ! Number of levels
   integer            , intent(in) :: mgmzp      ! Number of Grell's levels
   integer            , intent(in) :: cldd       ! Deepest cloud
   integer            , intent(in) :: clds       ! Shallowest cloud
   integer            , intent(in) :: maxens_cap ! Number of static control realisations
   !----- Local variables. ----------------------------------------------------------------!
   integer                         :: icld       ! Cloud size counter
   integer                         :: icap       ! Static control counter.
   integer                         :: klod       ! Level of origin of downdrafts
   integer                         :: klou       ! Level of origin of updrafts
   real                            :: sumdnarea  ! Total area for a SC
   real                            :: sumdnmf    ! Sum of downdraft MF for a given SC
   real                            :: sumuparea  ! Total area for a SC
   real                            :: sumupmf    ! Sum of updraft MF for a given SC
   real                            :: scalfac    ! Scaling factor
   !----- Local constants. ----------------------------------------------------------------!
   logical,parameter               :: do_rescale = .false.
   !---------------------------------------------------------------------------------------!
   
   cloudloop1: do icld = cldd, clds
      !----- Set the areas and drafts to zero, in case convection has failed. -------------!
      call azero(maxens_cap,ensemble_e(icld)%areadn_cap)
      call azero(maxens_cap,ensemble_e(icld)%areaup_cap)
      call azero(maxens_cap,ensemble_e(icld)%wdndraft_cap)
      call azero(maxens_cap,ensemble_e(icld)%wupdraft_cap)


      !------------------------------------------------------------------------------------!
      !     Compute the maximum relative area covered by downdrafts and updrafts for each  !
      ! static control if no other cloud were present.                                     !
      !------------------------------------------------------------------------------------!
      stacloop1: do icap=1,maxens_cap
         if (ensemble_e(icld)%ierr_cap(icap) == 0) then
            call grell_draft_area(m1,mgmzp,kgoff                                           &
                                 ,ensemble_e(icld)%comp_down_cap(icap)                     &
                                 ,ensemble_e(icld)%klod_cap(icap)                          &
                                 ,ensemble_e(icld)%klou_cap(icap)                          &
                                 ,ensemble_e(icld)%klfc_cap(icap)                          &
                                 ,ensemble_e(icld)%klnb_cap(icap)                          &
                                 ,ensemble_e(icld)%ktop_cap(icap)                          &
                                 ,dzu_cld,wwind,tke,sigw                                   &
                                 ,ensemble_e(icld)%wbuoymin_cap(icap)                      &
                                 ,ensemble_e(icld)%etad_cld_cap(1:mgmzp,icap)              &
                                 ,ensemble_e(icld)%mentrd_rate_cap(1:mgmzp,icap)           &
                                 ,ensemble_e(icld)%cdd_cap(1:mgmzp,icap)                   &
                                 ,ensemble_e(icld)%dbyd_cap(1:mgmzp,icap)                  &
                                 ,ensemble_e(icld)%rhod_cld_cap(1:mgmzp,icap)              &
                                 ,ensemble_e(icld)%dnmf_cap(icap)                          &
                                 ,ensemble_e(icld)%etau_cld_cap(1:mgmzp,icap)              &
                                 ,ensemble_e(icld)%mentru_rate_cap(1:mgmzp,icap)           &
                                 ,ensemble_e(icld)%cdu_cap(1:mgmzp,icap)                   &
                                 ,ensemble_e(icld)%dbyu_cap(1:mgmzp,icap)                  &
                                 ,ensemble_e(icld)%rhou_cld_cap(1:mgmzp,icap)              &
                                 ,ensemble_e(icld)%upmf_cap(icap)                          &
                                 ,ensemble_e(icld)%areadn_cap(icap)                        &
                                 ,ensemble_e(icld)%areaup_cap(icap)                        &
                                 ,ensemble_e(icld)%wdndraft_cap(icap)                      &
                                 ,ensemble_e(icld)%wupdraft_cap(icap))
         end if
      end do stacloop1
   end do cloudloop1

   !---------------------------------------------------------------------------------------!
   !    Now we loop again, this time in the reverse order, and rescale the areas because   !
   ! the areas overlap.  I.e., if the updraft that generates a cloud of size 1 has 20%     !
   ! of the area that may support, and cloud of size 2 has 30%, that does not mean that    !
   ! their areas are 20% and 30% of the grid respectively.  That would be true if they     !
   ! were the only clouds.  Because there is overlap then we must scale the areas.         !
   !---------------------------------------------------------------------------------------!
   stacloop2: do icap = 1, maxens_cap
      sumdnarea = 0.
      sumuparea = 0.
      sumdnmf   = 0.
      sumupmf   = 0.
      !----- Find the total area and the sum of all areas. --------------------------------!
      cloudloop2: do icld = cldd,clds
         if (ensemble_e(icld)%ierr_cap(icap) == 0) then
            sumdnmf   = sumdnmf + ensemble_e(icld)%dnmf_cap(icap)
            sumupmf   = sumupmf + ensemble_e(icld)%upmf_cap(icap)
         end if
      end do cloudloop2
      sumdnarea = 0.
      sumuparea = 0.
      !----- Scale the areas. -------------------------------------------------------------!
      cloudloop3: do icld = cldd,clds
         klod = ensemble_e(icld)%klod_cap(icap)
         klou = ensemble_e(icld)%klou_cap(icap)
         !----- Find the scaled area and reference downdraft. -----------------------------!
         if (do_rescale) then
            if ( ensemble_e(icld)%comp_down_cap(icap) .and.                                &
                 ensemble_e(icld)%ierr_cap(icap) == 0 .and. sumdnmf > 0. ) then

               ensemble_e(icld)%areadn_cap(icap)   = ensemble_e(icld)%areadn_cap(icap)     &
                                                   * ensemble_e(icld)%dnmf_cap(icap)       &
                                                   / sumdnmf
            end if

            !----- Find the scaled area and reference updraft. ----------------------------!
            if (ensemble_e(icld)%ierr_cap(icap) == 0 .and. sumupmf > 0.) then
               ensemble_e(icld)%areaup_cap(icap)   = ensemble_e(icld)%areaup_cap(icap)     &
                                                   * ensemble_e(icld)%upmf_cap(icap)       &
                                                   / sumupmf
            end if
         end if
         sumdnarea = sumdnarea + ensemble_e(icld)%areadn_cap(icap)
         sumuparea = sumuparea + ensemble_e(icld)%areaup_cap(icap)
      end do cloudloop3

      !----- Avoid the total area of this static control to exceed 99% of the area. -------!
      if (sumdnarea + sumuparea > 0.99) then
         scalfac = 0.99 / (sumdnarea + sumuparea) 
         cloudloop4: do icld = cldd,clds
            ensemble_e(icld)%areadn_cap(icap) = ensemble_e(icld)%areadn_cap(icap) * scalfac
            ensemble_e(icld)%areaup_cap(icap) = ensemble_e(icld)%areaup_cap(icap) * scalfac
         end do cloudloop4
      end if
   end do stacloop2

   !---------------------------------------------------------------------------------------!
   !    Final loop, this time we find the characteristic downdraft and updraft velocities, !
   ! which are currently simple diagnostic variables.                                      !
   !---------------------------------------------------------------------------------------!
   cloudloop5: do icld = cldd,clds
      stacloop3: do icap = 1, maxens_cap
         klod = ensemble_e(icld)%klod_cap(icap)
         klou = ensemble_e(icld)%klou_cap(icap)
         !----- Find the scaled area and reference downdraft. -----------------------------!
         if ( ensemble_e(icld)%comp_down_cap(icap) .and.                                   &
              ensemble_e(icld)%ierr_cap(icap) == 0 .and. sumdnmf > 0. ) then
            ensemble_e(icld)%wdndraft_cap(icap) = ensemble_e(icld)%dnmf_cap(icap)          &
                                                * ensemble_e(icld)%etad_cld_cap(klod,icap) &
                                                / (ensemble_e(icld)%rhod_cld_cap(klod,icap)&
                                                  *ensemble_e(icld)%areadn_cap(icap) )
         end if

         !----- Find the scaled area and reference updraft. -------------------------------!
         if (ensemble_e(icld)%ierr_cap(icap) == 0 .and. sumupmf > 0.) then
            ensemble_e(icld)%wupdraft_cap(icap) = ensemble_e(icld)%upmf_cap(icap)          &
                                                * ensemble_e(icld)%etau_cld_cap(klou,icap) &
                                                / (ensemble_e(icld)%rhou_cld_cap(klou,icap)&
                                                  *ensemble_e(icld)%areaup_cap(icap) )
         end if
      end do stacloop3
   end do cloudloop5


   return
end subroutine grell_cupar_area_scaler
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine organises the variables that go to the output, namely the heating and !
! moistening rates due to convection.                                                      !
!------------------------------------------------------------------------------------------!
subroutine grell_cupar_output(m1,mgmzp,maxens_cap,rtgt,zm,zt,dnmf,upmf,dnmx,upmx           &
                             ,ierr_cap,dnmf_cap,upmf_cap,kdet_cap,klou_cap,klcl_cap        &
                             ,klfc_cap,klod_cap,klnb_cap,ktop_cap,areadn_cap,areaup_cap    &
                             ,wdndraft_cap,wupdraft_cap,wbuoymin_cap,etad_cld_cap          &
                             ,mentrd_rate_cap,cdd_cap,dbyd_cap,rhod_cld_cap,etau_cld_cap   &
                             ,mentru_rate_cap,cdu_cap,dbyu_cap,rhou_cld_cap,qliqd_cld_cap  &
                             ,qliqu_cld_cap,qiced_cld_cap,qiceu_cld_cap,outco2,outqtot     &
                             ,outthil,precip,xierr,zklod,zklou,zklcl,zklfc,zkdt,zklnb      &
                             ,zktop,conprr,thsrc,rtsrc,co2src,areadn,areaup,wdndraft       &
                             ,wupdraft,wbuoymin,cuprliq,cuprice,i,j,icld,mynum)
   use mem_ensemble     , only : ensemble_vars ! ! type
   use mem_scratch_grell, only : kgoff         & ! intent(in)
                               , mkx           ! ! intent(in)
   use rconstants       , only : toodry        & ! intent(in)
                               , day_sec       ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer            , intent(in)  :: m1          ! Number of levels
   integer            , intent(in)  :: mgmzp       ! Number of Grell's levels
   integer            , intent(in)  :: maxens_cap  ! Number of static control realisations
   real               , intent(in)  :: rtgt        ! Corr. to get the heights     [   ----]
   real, dimension(m1), intent(in)  :: zm          ! Height at thermodyn. levels  [      m]
   real, dimension(m1), intent(in)  :: zt          ! Height at momentum levels    [      m]
   real               , intent(in)  :: dnmf        ! Reference downdraft mass flux[kg/m²/s]
   real               , intent(in)  :: upmf        ! Reference updraft mass flux  [kg/m²/s]
   real               , intent(in)  :: dnmx        ! Potential downdraft mass flux[kg/m²/s]
   real               , intent(in)  :: upmx        ! Potential updraft mass flux  [kg/m²/s]
   integer            , intent(in)  :: i           ! X-position
   integer            , intent(in)  :: j           ! Y-position
   integer            , intent(in)  :: icld        ! Cloud size
   integer            , intent(in)  :: mynum       ! Node number
   !----- Input, ensemble structure variables. --------------------------------------------!
   integer, dimension      (maxens_cap), intent(in) :: ierr_cap
   integer, dimension      (maxens_cap), intent(in) :: kdet_cap
   integer, dimension      (maxens_cap), intent(in) :: klou_cap
   integer, dimension      (maxens_cap), intent(in) :: klcl_cap
   integer, dimension      (maxens_cap), intent(in) :: klfc_cap
   integer, dimension      (maxens_cap), intent(in) :: klod_cap
   integer, dimension      (maxens_cap), intent(in) :: klnb_cap
   integer, dimension      (maxens_cap), intent(in) :: ktop_cap
   real   , dimension      (maxens_cap), intent(in) :: dnmf_cap
   real   , dimension      (maxens_cap), intent(in) :: upmf_cap
   real   , dimension      (maxens_cap), intent(in) :: areadn_cap
   real   , dimension      (maxens_cap), intent(in) :: areaup_cap
   real   , dimension      (maxens_cap), intent(in) :: wdndraft_cap
   real   , dimension      (maxens_cap), intent(in) :: wupdraft_cap
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
   real   , dimension(mgmzp)           , intent(in) :: outco2
   real   , dimension(mgmzp)           , intent(in) :: outqtot
   real   , dimension(mgmzp)           , intent(in) :: outthil
   real   , dimension(    1)           , intent(in) :: precip
   !----- Output variables. ---------------------------------------------------------------!
   real, dimension(m1), intent(inout) :: thsrc     ! Potential temperature fdbck  [    K/s]
   real, dimension(m1), intent(inout) :: rtsrc     ! Total mixing ratio feedback  [kg/kg/s]
   real, dimension(m1), intent(inout) :: co2src    ! Total CO2 mixing ratio fdbk  [  ppm/s]
   real, dimension(m1), intent(inout) :: cuprliq   ! Cumulus water mixing ratio   [  kg/kg]
   real, dimension(m1), intent(inout) :: cuprice   ! Cumulus ice mixing ratio     [  kg/kg]
   real               , intent(inout) :: areadn    ! Fractional downdraft area    [    ---]
   real               , intent(inout) :: areaup    ! Fractional updraft area      [    ---]
   real               , intent(inout) :: wdndraft  ! Expected downdraft at LOD    [    m/s]
   real               , intent(inout) :: wupdraft  ! Expected updraft at LOU      [    m/s]
   real               , intent(inout) :: wbuoymin  ! Minimum buoyant velocity     [    m/s]
   real               , intent(inout) :: conprr    ! Rate of convective precip.   [kg/m²/s]
   real               , intent(inout) :: xierr     ! Error flag                   [    ---]
   real               , intent(inout) :: zklod     ! Level of origin of downdrafts[      m]
   real               , intent(inout) :: zklou     ! Level of origin of updrafts  [      m]
   real               , intent(inout) :: zklcl     ! Lifting condensation level   [      m]
   real               , intent(inout) :: zklfc     ! Level of free convection     [      m]
   real               , intent(inout) :: zkdt      ! Top of the dndraft detr.     [      m]
   real               , intent(inout) :: zklnb     ! Cloud top                    [      m]
   real               , intent(inout) :: zktop     ! Cloud top                    [      m]
   !----- Local variables. ----------------------------------------------------------------!
   logical, dimension(maxens_cap)     :: is_cloud  ! This static control has a cloud
   integer                            :: icap      ! Static control counter.
   real                               :: nmok      ! # of members that had clouds.
   integer                            :: k         ! Cloud level counter
   integer                            :: kr        ! BRAMS level counter
   real                               :: exner     ! Exner fctn. for tend. conv.  [ J/kg/K]
   real                               :: nmoki     ! 1/nmok
   real                               :: zhgt      ! Height
   real                               :: f_dn      ! Fraction for downdraft
   real                               :: f_up      ! Fraction for updraft
   integer                            :: klod      ! Downdraft origin
   integer                            :: klou      ! Updraft origin
   integer                            :: klcl      ! Lifting condensation level
   integer                            :: klfc      ! Level of free convection
   integer                            :: kdet      ! Origin of downdraft detrainment
   integer                            :: klnb      ! Cloud top
   integer                            :: ktop      ! Cloud top
   real, dimension(m1)                :: cuprliqd  ! Cumulus water mixing ratio   [  kg/kg]
   real, dimension(m1)                :: cupriced  ! Cumulus ice mixing ratio     [  kg/kg]
   real, dimension(m1)                :: cuprliqu  ! Cumulus water mixing ratio   [  kg/kg]
   real, dimension(m1)                :: cupriceu  ! Cumulus ice mixing ratio     [  kg/kg]
   !----- Aux. variables for fractional area. ---------------------------------------------! 
   real, dimension(m1,maxens_cap)   :: cuprliq_cap  ! Cumulus water mixing ratio  [  kg/kg]
   real, dimension(m1,maxens_cap)   :: cuprice_cap  ! Cumulus ice mixing ratio    [  kg/kg]
   real, dimension(m1,maxens_cap)   :: cuprliqd_cap ! Cumulus water mixing ratio  [  kg/kg]
   real, dimension(m1,maxens_cap)   :: cupriced_cap ! Cumulus ice mixing ratio    [  kg/kg]
   real, dimension(m1,maxens_cap)   :: cuprliqu_cap ! Cumulus water mixing ratio  [  kg/kg]
   real, dimension(m1,maxens_cap)   :: cupriceu_cap ! Cumulus ice mixing ratio    [  kg/kg]
   !----- Local constants, for debugging. -------------------------------------------------!
   integer                          :: iun
   logical          , parameter     :: print_debug=.false.
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !    Flushing all variables to zero in case convection didn't happen.                   !
   !---------------------------------------------------------------------------------------!
   thsrc         = 0.
   rtsrc         = 0.
   co2src        = 0.
   cuprliq       = 0.
   cuprice       = 0.
   cuprliqd      = 0.
   cupriced      = 0.
   cuprliqu      = 0.
   cupriceu      = 0.
   cuprliq_cap   = 0.
   cuprice_cap   = 0.
   cuprliqd_cap  = 0.
   cupriced_cap  = 0.
   cuprliqu_cap  = 0.
   cupriceu_cap  = 0.

   areadn        = 0.
   areaup        = 0.
   conprr        = 0.
   zkdt          = 0.
   zklou         = 0.
   zklcl         = 0.
   zklfc         = 0.
   zklod         = 0.
   zklnb         = 0.
   zktop         = 0.
   wdndraft      = 0.
   wupdraft      = 0.
   wbuoymin      = 0.

   !---------------------------------------------------------------------------------------!
   !    Find the number of clouds and the cloud mask.                                      !
   !---------------------------------------------------------------------------------------!
   is_cloud(:) = ierr_cap(:) == 0
   nmok        = real(count(is_cloud))

   !---------------------------------------------------------------------------------------!
   !   Copying the error flag. This should not be zero in case convection failed.          !
   !---------------------------------------------------------------------------------------!
   if (upmf > 0. .and. nmok > 0. ) then
      xierr = 0.
      nmoki = 1./nmok
   else
      xierr = real(ierr_cap(1))
      return !----- No cloud, we don't return any variable. -------------------------------!
   end if
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     We choose the level to be the most level that will encompass all non-zero         !
   ! members.                                                                              !
   !---------------------------------------------------------------------------------------!
   kdet  = max(1,min(mkx,nint(real(sum(kdet_cap ,mask = is_cloud))/nmok)))
   klou  = max(1,min(mkx,nint(real(sum(klou_cap ,mask = is_cloud))/nmok)))
   klcl  = max(1,min(mkx,nint(real(sum(klcl_cap ,mask = is_cloud))/nmok)))
   klfc  = max(1,min(mkx,nint(real(sum(klfc_cap ,mask = is_cloud))/nmok)))
   klod  = max(1,min(mkx,nint(real(sum(klod_cap ,mask = is_cloud))/nmok)))
   klnb  = max(1,min(mkx,nint(real(sum(klnb_cap ,mask = is_cloud))/nmok)))
   ktop  = max(1,min(mkx,nint(real(sum(ktop_cap ,mask = is_cloud))/nmok)))
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    Fix the levels, here I will add back the offset so the output will be consistent.  !
   ! I shall return these variables even when no cloud developed for debugging purposes.   !
   !---------------------------------------------------------------------------------------!
   zkdt   = (zt(kdet  + kgoff)-zm(kgoff))*rtgt
   zklou  = (zt(klou  + kgoff)-zm(kgoff))*rtgt
   zklcl  = (zt(klcl  + kgoff)-zm(kgoff))*rtgt
   zklfc  = (zt(klfc  + kgoff)-zm(kgoff))*rtgt
   zklod  = (zt(klod  + kgoff)-zm(kgoff))*rtgt
   zklnb  = (zt(klnb  + kgoff)-zm(kgoff))*rtgt
   zktop  = (zt(ktop  + kgoff)-zm(kgoff))*rtgt
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    Precipitation is simply copied, it could even be output directly from the main     !
   ! subroutine, brought here just to be together with the other source terms.             !
   !---------------------------------------------------------------------------------------!
   conprr = precip(1)
   
   do k=1,mkx
      kr    = k + kgoff
      !----- Here we are simply copying including the offset back. ------------------------!
      thsrc(kr)   = outthil(k)
      rtsrc(kr)   = outqtot(k)
      co2src(kr)  = outco2(k)
   end do
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !   Compute the relative area covered by downdrafts and updrafts.                       !
   !---------------------------------------------------------------------------------------!
   stacloop: do icap=1,maxens_cap
      if (ierr_cap(icap) == 0) then

         !---------------------------------------------------------------------------------!
         !     Find the fraction due to downdrafts and updrafts.                           !
         !---------------------------------------------------------------------------------!
         f_dn = areadn_cap(icap) / (areadn_cap(icap) + areaup_cap(icap))
         f_up = areaup_cap(icap) / (areadn_cap(icap) + areaup_cap(icap))
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !   Compute the cloud condensed mixing ratio for this realisation.  This is used  !
         ! by Harrington when cumulus feedback is requested, so we rescale the liquid      !
         ! water at the downdrafts and updrafts by their area.  Levels outside the cloud   !
         ! range may contain environment values, and we don't want to double account.      !
         !---------------------------------------------------------------------------------!
         do k=1,klod_cap(icap)
            kr=k+kgoff
            cuprliq_cap (kr,icap) = cuprliq_cap (kr,icap) + qliqd_cld_cap(k,icap) * f_dn
            cuprice_cap (kr,icap) = cuprice_cap (kr,icap) + qiced_cld_cap(k,icap) * f_dn
            cuprliqd_cap(kr,icap) = qliqd_cld_cap(k,icap)
            cupriced_cap(kr,icap) = qiced_cld_cap(k,icap)
         end do
         do k=klou_cap(icap),ktop_cap(icap)
            kr=k+kgoff
            cuprliq_cap (kr,icap) = cuprliq_cap (kr,icap) + qliqu_cld_cap(k,icap) * f_up
            cuprice_cap (kr,icap) = cuprice_cap (kr,icap) + qiceu_cld_cap(k,icap) * f_up
            cuprliqu_cap(kr,icap) = qliqu_cld_cap(k,icap)
            cupriceu_cap(kr,icap) = qiceu_cld_cap(k,icap)
         end do
      end if
   end do stacloop
   !---------------------------------------------------------------------------------------!



   !----- Find the averaged area. ---------------------------------------------------------! 
   areadn   = sum(areadn_cap)   * nmoki
   areaup   = sum(areaup_cap)   * nmoki
   wdndraft = sum(wdndraft_cap) * nmoki
   wupdraft = sum(wupdraft_cap) * nmoki
   wbuoymin = sum(wbuoymin_cap) * nmoki
   do icap=1,maxens_cap
      do kr=1,m1
         cuprliq (kr) = cuprliq (kr) + cuprliq_cap (kr,icap) * nmoki
         cuprice (kr) = cuprice (kr) + cuprice_cap (kr,icap) * nmoki
         cuprliqd(kr) = cuprliqd(kr) + cuprliqd_cap(kr,icap) * nmoki
         cupriced(kr) = cupriced(kr) + cupriced_cap(kr,icap) * nmoki
         cuprliqu(kr) = cuprliqu(kr) + cuprliqu_cap(kr,icap) * nmoki
         cupriceu(kr) = cupriceu(kr) + cupriceu_cap(kr,icap) * nmoki
      end do
   end do
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     If printing debug, check whether the cloud happened and print the cloud           !
   ! characteristics.                                                                      !
   !---------------------------------------------------------------------------------------!
   if (print_debug) then
      write(unit=20+icld,fmt='(143a)'            ) ('-',k=1,143)
      write(unit=20+icld,fmt='(a)'               ) ''
      write(unit=20+icld,fmt='(a,1x,i5)'         ) '  I        =',i
      write(unit=20+icld,fmt='(a,1x,i5)'         ) '  J        =',j
      write(unit=20+icld,fmt='(a,1x,i5)'         ) '  NMOK     =',nint(nmok)
      write(unit=20+icld,fmt='(a)'               ) ''
      write(unit=20+icld,fmt='(a,1x,i5,1x,f10.2)') '  DET      =',kdet,zkdt
      write(unit=20+icld,fmt='(a,1x,i5,1x,f10.2)') '  LOU      =',klou,zklou
      write(unit=20+icld,fmt='(a,1x,i5,1x,f10.2)') '  LCL      =',klcl,zklcl
      write(unit=20+icld,fmt='(a,1x,i5,1x,f10.2)') '  LFC      =',klfc,zklfc
      write(unit=20+icld,fmt='(a,1x,i5,1x,f10.2)') '  LOD      =',klod,zklod
      write(unit=20+icld,fmt='(a,1x,i5,1x,f10.2)') '  LNB      =',klnb,zklnb
      write(unit=20+icld,fmt='(a,1x,i5,1x,f10.2)') '  TOP      =',ktop,zktop
      write(unit=20+icld,fmt='(a)'               ) ''
      write(unit=20+icld,fmt='(a,1x,es12.5)'     ) '  DNMF     =',dnmf
      write(unit=20+icld,fmt='(a,1x,es12.5)'     ) '  UPMF     =',upmf
      write(unit=20+icld,fmt='(a,1x,es12.5)'     ) '  DNMX     =',dnmx
      write(unit=20+icld,fmt='(a,1x,es12.5)'     ) '  UPMX     =',upmx
      write(unit=20+icld,fmt='(a,1x,es12.5)'     ) '  AREADN   =',areadn
      write(unit=20+icld,fmt='(a,1x,es12.5)'     ) '  AREAUP   =',areaup
      write(unit=20+icld,fmt='(a,1x,es12.5)'     ) '  WDNDRAFT =',wdndraft
      write(unit=20+icld,fmt='(a,1x,es12.5)'     ) '  WUPDRAFT =',wupdraft
      write(unit=20+icld,fmt='(a,1x,es12.5)'     ) '  WBUOYMIN =',wbuoymin
      write(unit=20+icld,fmt='(a,1x,es12.5)'     ) '  CONPRR   =',conprr * day_sec
      write(unit=20+icld,fmt='(a)'               ) ''
      write(unit=20+icld,fmt='(11(a,1x))')  '           K','      HEIGHT','       THSRC'   &
                                           ,'       RTSRC','      CO2SRC','     CUPRICE'   &
                                           ,'     CUPRLIQ','  CUPRICE_DN','  CUPRLIQ_DN'   &
                                           ,'  CUPRICE_UP','  CUPRLIQ_UP'
      write(unit=20+icld,fmt='(143(a))') ('-',k=1,143)
      do k=m1,1,-1
         zhgt = ( zt(k+kgoff) - zm(kgoff) ) * rtgt
         write (unit=20+icld,fmt='(i12,1x,f12.2,1x,9(es12.5,1x))')                         &
                                                              k                            &
                                                            , zhgt                         &
                                                            , thsrc  (k)  * day_sec        &
                                                            , rtsrc  (k)  * day_sec*1000.  &
                                                            , co2src (k)  * day_sec        &
                                                            , cuprice(k)  * 1000.          &
                                                            , cuprliq(k)  * 1000.          &
                                                            , cupriced(k) * 1000.          &
                                                            , cuprliqd(k) * 1000.          &
                                                            , cupriceu(k) * 1000.          &
                                                            , cuprliqu(k) * 1000.
      end do
      write(unit=20+icld,fmt='(143(a))') ('-',k=1,143)
      write(unit=20+icld,fmt='(a)'               ) ''
   end if
   !---------------------------------------------------------------------------------------!


   return
end subroutine grell_cupar_output
!==========================================================================================!
!==========================================================================================!
