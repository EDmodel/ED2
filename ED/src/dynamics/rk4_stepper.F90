!==========================================================================================!
!==========================================================================================!
!    This module contains some subroutines used to compute the stages for advancing the    !
! Runge-Kutta step.                                                                        !
!------------------------------------------------------------------------------------------!
module rk4_stepper_ar

   contains 
   !=======================================================================================!
   !=======================================================================================!
   !   This subroutine is the main Runge-Kutta step driver.                                !
   !---------------------------------------------------------------------------------------!
   subroutine rkqs_ar(integration_buff, x,htry,hdid,hnext,csite,ipa,isi,ipy,ifm,rhos,vels  &
                     ,atm_tmp,atm_shv,atm_co2,geoht,exner,pcpg,qpcpg,dpcpg, prss,lsl)

      use ed_state_vars , only : sitetype            & ! structure
                               , patchtype           & ! structure
                               , edtype              & ! structure
                               , rk4patchtype        & ! structure
                               , integration_vars_ar & ! structure
                               , edgrid_g            ! ! structure
      use rk4_coms      , only : hmin                & ! intent(in)
                               , rk4eps              & ! intent(in)
                               , rk4epsi             & ! intent(in)
                               , safety              & ! intent(in)
                               , pgrow               & ! intent(in)
                               , pshrnk              & ! intent(in)
                               , errcon              ! ! intent(in)

      implicit none
      !----- Arguments --------------------------------------------------------------------!
      integer                  , intent(in)    :: lsl, ipa,isi,ipy,ifm
      type(sitetype)           , target        :: csite
      type(integration_vars_ar), target        :: integration_buff
      real                     , intent(in)    :: rhos,vels,atm_tmp,atm_shv,atm_co2,geoht
      real                     , intent(in)    :: exner,pcpg,qpcpg,dpcpg,prss,htry
      real                     , intent(inout) :: x
      real                     , intent(out)   :: hdid,hnext
      !----- Local variables --------------------------------------------------------------!
      type(edtype)             , pointer       :: cgrid
      real                                     :: h,errmax,xnew,newh,oldh
      logical                                  :: reject_step,minstep,stuck,test_reject
      integer                                  :: k
      !------------------------------------------------------------------------------------!

      cgrid       => edgrid_g(ifm)
      h           =  htry
      reject_step =  .false.
      hstep:   do

         !---------------------------------------------------------------------------------!
         ! 1. Try a step of varying size.                                                  !
         !---------------------------------------------------------------------------------!
         call rkck_ar(integration_buff%y,integration_buff%dydx,integration_buff%ytemp      &
                     ,integration_buff%yerr,integration_buff%ak2,integration_buff%ak3      &
                     ,integration_buff%ak4,integration_buff%ak5,integration_buff%ak6       &
                     ,integration_buff%ak7,x,h,csite,ipa,isi,ipy,reject_step,rhos          &
                     ,vels,atm_tmp,atm_shv,atm_co2,geoht,exner,pcpg,qpcpg,dpcpg,prss,lsl)


         !---------------------------------------------------------------------------------!
         ! 2. Check to see how accurate the step was.  Errors were calculated by integrat- !
         !    ing the derivative of that last step.                                        !
         !---------------------------------------------------------------------------------!
         if (reject_step) then
            !------------------------------------------------------------------------------!
            !    If step was already rejected, that means the step had finished premature- !
            ! ly, so we assign a standard large error (10.0).                              !
            !------------------------------------------------------------------------------!
            errmax = 10.0
         else
            call get_errmax_ar(errmax, integration_buff%yerr,integration_buff%yscal        &
                              ,csite%patch(ipa),lsl,integration_buff%y                     &
                              ,integration_buff%ytemp)
            errmax = errmax * rk4epsi
         end if

         !---------------------------------------------------------------------------------!
         ! 3. If that error was large, then calculate a new step size to try.  There are   !
         !    two types of new tries.  If step failed to be minimally reasonable (reject-  !
         !    ed) we have assigned a standard large error (10.0).  Otherwise a new step is !
         !    calculated based on the size of that error.  Hopefully, those new steps      !
         !    should be less then the previous h.  If the error was small, i.e. less then  !
         !    rk4eps, then we are done with this step, and we can move forward             !
         !    time: x = x + h                                                              !
         !---------------------------------------------------------------------------------!
         if (errmax > 1.0) then
            !----- Defining new step and checking if it can be ----------------------------!
            oldh    = h
            newh    = safety * h * errmax**pshrnk
            minstep = (newh == h)

            !----- Defining next time, and checking if it really added something ----------!
            h       = max(0.1*h, newh)
            xnew    = x + h
            stuck   = xnew == x

         !---------------------------------------------------------------------------------!
         ! 4. Here is the moment of truth... If we reached a tiny step and yet the model   !
         !    didn't converge, then we print various values to inform the user and abort   !
         !    the run.  Please, don't hate the messenger.                                  !
         !---------------------------------------------------------------------------------!
            if (minstep .or. stuck) then

               write (unit=*,fmt='(80a)')         ('=',k=1,80)
               write (unit=*,fmt='(a)')           '   STEPSIZE UNDERFLOW IN RKQS'
               write (unit=*,fmt='(80a)')         ('-',k=1,80)
               write (unit=*,fmt='(a,1x,f9.4)')   ' + LONGITUDE:   ',cgrid%lon(ipy)
               write (unit=*,fmt='(a,1x,f9.4)')   ' + LATITUDE:    ',cgrid%lat(ipy)
               write (unit=*,fmt='(a,1x,i6)')     ' + POLYGON:     ',ipy
               write (unit=*,fmt='(a)')           ' + PATCH AGE:   ',csite%age(ipa)
               write (unit=*,fmt='(a,1x,es12.5)') '   - AGE:       ',csite%age(ipa)
               write (unit=*,fmt='(a,1x,i6)')     '   - DIST_TYPE: ',csite%dist_type(ipa)
               write (unit=*,fmt='(a,1x,l1)')     ' + REJECT_STEP: ',reject_step
               write (unit=*,fmt='(a,1x,l1)')     ' + MINSTEP:     ',minstep
               write (unit=*,fmt='(a,1x,l1)')     ' + STUCK:       ',stuck
               write (unit=*,fmt='(a,1x,es12.5)') ' + ERRMAX:      ',errmax
               write (unit=*,fmt='(a,1x,es12.5)') ' + X:           ',x
               write (unit=*,fmt='(a,1x,es12.5)') ' + H:           ',h
               write (unit=*,fmt='(a,1x,es12.5)') ' + OLDH:        ',oldh
               write (unit=*,fmt='(a,1x,es12.5)') ' + NEWH:        ',newh
               write (unit=*,fmt='(a,1x,es12.5)') ' + SAFETY:      ',safety
               write (unit=*,fmt='(a,1x,es12.5)') ' + ERRMAX:      ',errmax
               write (unit=*,fmt='(80a)') ('-',k=1,80)
               if(reject_step)then
                  write (unit=*,fmt='(a)') '   Likely to be a rejected step problem.'
                  write (unit=*,fmt='(80a)') ('=',k=1,80)
               else
                  write (unit=*,fmt='(a)') '   Likely to be an errmax problem.'
                  write (unit=*,fmt='(80a)') ('=',k=1,80)
               endif

               if (reject_step) then
                  !----- Run the LSM sanity check but this time we force the print. -------!
                  call lsm_sanity_check_ar(integration_buff%ytemp,test_reject,csite,ipa    &
                                          ,lsl,integration_buff%dydx,h,atm_tmp,atm_shv     &
                                          ,atm_co2,prss,exner,rhos,vels,geoht,pcpg,qpcpg   &
                                          ,dpcpg,.true.)
                  if (.not. test_reject) then
                     call lsm_sanity_check_ar(integration_buff%ak7,test_reject,csite,ipa   &
                                             ,lsl,integration_buff%dydx,h,atm_tmp,atm_shv  &
                                             ,atm_co2,prss,exner,rhos,vels,geoht,pcpg      &
                                             ,qpcpg,dpcpg,.true.)
                  end if
                  call print_sanity_check_ar(integration_buff%y,csite,ipa,lsl)
               else
                  call print_errmax_ar(errmax,integration_buff%yerr,integration_buff%yscal &
                                      ,csite%patch(ipa),lsl,integration_buff%y             &
                                      ,integration_buff%ytemp)
                  write (unit=*,fmt='(80a)') ('=',k=1,80)
                  write (unit=*,fmt='(a,1x,es12.5)') ' - Rel. errmax:',errmax*rk4epsi
                  write (unit=*,fmt='(a,1x,es12.5)') ' - Raw errmax: ',errmax
                  write (unit=*,fmt='(a,1x,es12.5)') ' - Epsilon:',rk4eps
                  write (unit=*,fmt='(80a)') ('=',k=1,80)
               end if
               call print_patch_ar(integration_buff%y, csite,ipa, lsl,atm_tmp,atm_shv      &
                                  ,atm_co2,prss,exner,rhos,vels,geoht,pcpg,qpcpg,dpcpg)
            endif
         
         else
            !------------------------------------------------------------------------------!
            !   Great, it worked, so now we can advance to the next step, copy the stuff   !
            ! on the integration buffer to the actual patch/cohort structures, ans we will !
            ! also set up h for the next time.  And here we can relax h for the next step, !
            ! and try something faster.                                                    !
            !------------------------------------------------------------------------------!
            if (errmax > errcon) then
               hnext = safety * h * errmax**pgrow
            else
               hnext = 5.0 * h
            endif
            hnext = max(2.0*hmin,hnext)

            !----- Updating time ----------------------------------------------------------!
            x    = x + h
            hdid = h
            call copy_rk4_patch_ar(integration_buff%ytemp,integration_buff%y               &
                                  ,csite%patch(ipa),lsl)
            exit hstep
         end if
      end do hstep

      return
   end subroutine rkqs_ar
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine will update the variables and perform the actual time stepping.    !
   !---------------------------------------------------------------------------------------!
   subroutine rkck_ar(y,dydx,yout,yerr,ak2,ak3,ak4,ak5,ak6,ak7,x,h,csite,ipa,isi,ipy       &
                     ,reject_step,rhos,vels,atm_tmp,atm_shv,atm_co2,geoht,exner,pcpg,qpcpg &
                     ,dpcpg,prss,lsl)

      use ed_state_vars , only : sitetype            & ! structure
                               , patchtype           & ! structure
                               , rk4patchtype        & ! structure
                               , integration_vars_ar ! ! structure
      use rk4_coms      , only : print_diags         & ! intent(in)
                               , a2                  & ! intent(in)
                               , a3                  & ! intent(in)
                               , a4                  & ! intent(in)
                               , a5                  & ! intent(in)
                               , a6                  & ! intent(in)
                               , b21                 & ! intent(in)
                               , b31                 & ! intent(in)
                               , b32                 & ! intent(in)
                               , b41                 & ! intent(in)
                               , b42                 & ! intent(in)
                               , b43                 & ! intent(in)
                               , b51                 & ! intent(in)
                               , b52                 & ! intent(in)
                               , b53                 & ! intent(in)
                               , b54                 & ! intent(in)
                               , b61                 & ! intent(in)
                               , b62                 & ! intent(in)
                               , b63                 & ! intent(in)
                               , b64                 & ! intent(in)
                               , b65                 & ! intent(in)
                               , c1                  & ! intent(in)
                               , c3                  & ! intent(in)
                               , c4                  & ! intent(in)
                               , c6                  & ! intent(in)
                               , dc5                 & ! intent(in)
                               , dc1                 & ! intent(in)
                               , dc3                 & ! intent(in)
                               , dc4                 & ! intent(in)
                               , dc6                 ! ! intent(in)
      implicit none

      !----- Arguments --------------------------------------------------------------------!
      integer           , intent(in)  :: lsl,ipa,isi,ipy
      real              , intent(in)  :: x,h,rhos,vels,atm_tmp,atm_shv,atm_co2
      real              , intent(in)  :: geoht,exner,pcpg,qpcpg,dpcpg,prss
      type(rk4patchtype), target      :: y,dydx,yout,yerr
      type(rk4patchtype), target      :: ak2,ak3,ak4,ak5,ak6,ak7
      type(sitetype)    , target      :: csite
      logical           , intent(out) :: reject_step
      !----- Local variables --------------------------------------------------------------!
      type(patchtype)   , pointer     :: cpatch
      !------------------------------------------------------------------------------------!


      !----- Interfaces, in case the model is compiled without forcing them. --------------!
#if USE_INTERF
      interface
         subroutine leaf_derivs_ar(initp,dydx,csite,ipa,isi,ipy,rhos,prss,pcpg,qpcpg,dpcpg &
                                  ,atm_tmp,exner,geoht,vels,atm_shv,atm_co2,lsl)
           
            use ed_state_vars ,only : sitetype,rk4patchtype,patchtype
            implicit none
            integer             , intent(in) :: lsl,ipa,isi,ipy
            real                , intent(in) :: rhos,prss,pcpg,qpcpg,dpcpg
            real                , intent(in) :: atm_tmp,exner,geoht,vels,atm_shv,atm_co2
            type (rk4patchtype) , target     :: initp
            type (rk4patchtype) , target     :: dydx
            type (sitetype)     , target     :: csite
         end subroutine leaf_derivs_ar
      end interface
#endif

      !------------------------------------------------------------------------------------!
      !     Start and assume that nothing went wrong up to this point... If we find any    !
      ! seriously bad step, quit and reduce the time step without even bothering to try    !
      ! further.                                                                           !
      !------------------------------------------------------------------------------------!
      reject_step = .false.
      cpatch => csite%patch(ipa)

      call copy_rk4_patch_ar(y, ak7, cpatch, lsl)
      call inc_rk4_patch_ar(ak7, dydx, b21*h, cpatch, lsl)
      !call adjust_veg_properties(ak7,dydx,csite,ipa,rhos,b21*h)
      call update_veg_properties(ak7,csite,ipa)
      call stabilize_snow_layers_ar(ak7, csite, ipa, lsl)
      call lsm_sanity_check_ar(ak7, reject_step, csite, ipa, lsl ,dydx,h,atm_tmp,atm_shv   &
                              ,atm_co2,prss,exner,rhos,vels,geoht,pcpg,qpcpg,dpcpg         &
                              ,print_diags)
      if (reject_step) return



      call leaf_derivs_ar(ak7, ak2, csite, ipa,isi,ipy, rhos, prss, pcpg, qpcpg, dpcpg     &
                         ,atm_tmp, exner, geoht, vels, atm_shv, atm_co2, lsl)
      call copy_rk4_patch_ar(y, ak7, cpatch, lsl)
      call inc_rk4_patch_ar(ak7, dydx, b31*h, cpatch, lsl)
      !call adjust_veg_properties(ak7,dydx,csite,ipa,rhos,b31*h)
      call inc_rk4_patch_ar(ak7, ak2, b32*h, cpatch, lsl)
      !call adjust_veg_properties(ak7,ak2,csite,ipa,rhos,b32*h)
      call update_veg_properties(ak7,csite,ipa)
      call stabilize_snow_layers_ar(ak7, csite,ipa, lsl)
      call lsm_sanity_check_ar(ak7,reject_step,csite,ipa, lsl,dydx,h,atm_tmp,atm_shv       &
                              ,atm_co2,prss,exner,rhos,vels,geoht,pcpg,qpcpg,dpcpg         &
                              ,print_diags)
      if (reject_step) return



      call leaf_derivs_ar(ak7, ak3, csite,ipa,isi,ipy, rhos, prss, pcpg, qpcpg, dpcpg      &
                         ,atm_tmp, exner, geoht, vels, atm_shv, atm_co2, lsl)
      call copy_rk4_patch_ar(y, ak7, cpatch, lsl)
      call inc_rk4_patch_ar(ak7, dydx, b41*h, cpatch, lsl)
      !call adjust_veg_properties(ak7,dydx,csite,ipa,rhos,b41*h)
      call inc_rk4_patch_ar(ak7, ak2, b42*h, cpatch, lsl)
      !call adjust_veg_properties(ak7,dydx,csite,ipa,rhos,b42*h)
      call inc_rk4_patch_ar(ak7, ak3, b43*h, cpatch, lsl)
      !call adjust_veg_properties(ak7,dydx,csite,ipa,rhos,b43*h)
      call update_veg_properties(ak7,csite,ipa)
      call stabilize_snow_layers_ar(ak7, csite,ipa, lsl)
      call lsm_sanity_check_ar(ak7, reject_step, csite,ipa, lsl,dydx,h,atm_tmp,atm_shv     &
                              ,atm_co2,prss,exner,rhos,vels,geoht,pcpg,qpcpg,dpcpg         &
                              ,print_diags )
      if (reject_step) return



      call leaf_derivs_ar(ak7, ak4, csite, ipa,isi,ipy, rhos, prss, pcpg, qpcpg, dpcpg     &
                         ,atm_tmp, exner, geoht, vels, atm_shv, atm_co2, lsl)
      call copy_rk4_patch_ar(y, ak7, cpatch, lsl)
      call inc_rk4_patch_ar(ak7, dydx, b51*h, cpatch, lsl)
      !call adjust_veg_properties(ak7,dydx,csite,ipa,rhos,b51*h)
      call inc_rk4_patch_ar(ak7, ak2, b52*h, cpatch, lsl)
      !call adjust_veg_properties(ak7,dydx,csite,ipa,rhos,b52*h)
      call inc_rk4_patch_ar(ak7, ak3, b53*h, cpatch, lsl)
      !call adjust_veg_properties(ak7,dydx,csite,ipa,rhos,b53*h)
      call inc_rk4_patch_ar(ak7, ak4, b54*h, cpatch, lsl)
      !call adjust_veg_properties(ak7,dydx,csite,ipa,rhos,b54*h)
      call update_veg_properties(ak7,csite,ipa)
      call stabilize_snow_layers_ar(ak7, csite, ipa, lsl)
      call lsm_sanity_check_ar(ak7,reject_step,csite,ipa,lsl,dydx,h,atm_tmp,atm_shv        &
                              ,atm_co2,prss,exner,rhos,vels,geoht,pcpg,qpcpg,dpcpg         &
                              ,print_diags)
      if (reject_step) return



      call leaf_derivs_ar(ak7, ak5, csite, ipa,isi,ipy, rhos, prss, pcpg, qpcpg, dpcpg     &
                         ,atm_tmp, exner, geoht, vels, atm_shv, atm_co2, lsl)
      call copy_rk4_patch_ar(y, ak7, cpatch, lsl)
      call inc_rk4_patch_ar(ak7, dydx, b61*h, cpatch, lsl)
      !call adjust_veg_properties(ak7,dydx,csite,ipa,rhos,b61*h)
      call inc_rk4_patch_ar(ak7, ak2, b62*h, cpatch, lsl)
      !call adjust_veg_properties(ak7,dydx,csite,ipa,rhos,b62*h)
      call inc_rk4_patch_ar(ak7, ak3, b63*h, cpatch, lsl)
      !call adjust_veg_properties(ak7,dydx,csite,ipa,rhos,b63*h)
      call inc_rk4_patch_ar(ak7, ak4, b64*h, cpatch, lsl)
      !call adjust_veg_properties(ak7,dydx,csite,ipa,rhos,b64*h)
      call inc_rk4_patch_ar(ak7, ak5, b65*h, cpatch, lsl)
      !call adjust_veg_properties(ak7,dydx,csite,ipa,rhos,b65*h)
      call update_veg_properties(ak7,csite,ipa)
      call stabilize_snow_layers_ar(ak7, csite,ipa, lsl)
      call lsm_sanity_check_ar(ak7, reject_step, csite,ipa, lsl,dydx,h,atm_tmp,atm_shv     &
                              ,atm_co2,prss,exner,rhos,vels,geoht,pcpg,qpcpg,dpcpg         &
                              ,print_diags)
      if(reject_step)return

      call leaf_derivs_ar(ak7, ak6, csite,ipa,isi,ipy, rhos, prss, pcpg, qpcpg, dpcpg      &
                         ,atm_tmp, exner, geoht, vels, atm_shv, atm_co2, lsl)
      call copy_rk4_patch_ar(y, yout, cpatch, lsl)
      call inc_rk4_patch_ar(yout, dydx, c1*h, cpatch, lsl)
      !call adjust_veg_properties(ak7,dydx,csite,ipa,rhos,c1*h)
      call inc_rk4_patch_ar(yout, ak3, c3*h, cpatch, lsl)
      !call adjust_veg_properties(ak7,dydx,csite,ipa,rhos,c3*h)
      call inc_rk4_patch_ar(yout, ak4, c4*h, cpatch, lsl)
      !call adjust_veg_properties(ak7,dydx,csite,ipa,rhos,c4*h)
      call inc_rk4_patch_ar(yout, ak6, c6*h, cpatch, lsl)
      !call adjust_veg_properties(ak7,dydx,csite,ipa,rhos,c6*h)
      call update_veg_properties(yout,csite,ipa)
      call stabilize_snow_layers_ar(yout, csite,ipa, lsl)
      call lsm_sanity_check_ar(yout, reject_step, csite,ipa, lsl,dydx,h,atm_tmp,atm_shv    &
                              ,atm_co2,prss,exner,rhos,vels,geoht,pcpg,qpcpg,dpcpg         &
                              ,print_diags)
      if(reject_step)return

      call copy_rk4_patch_ar(dydx, yerr, cpatch, lsl)
      call inc_rk4_patch_ar(yerr, dydx, dc1*h-1.0, cpatch, lsl)
      call inc_rk4_patch_ar(yerr, ak3, dc3*h, cpatch, lsl)
      call inc_rk4_patch_ar(yerr, ak4, dc4*h, cpatch, lsl)
      call inc_rk4_patch_ar(yerr, ak5, dc5*h, cpatch, lsl)
      call inc_rk4_patch_ar(yerr, ak6, dc6*h, cpatch, lsl )

      return
   end subroutine rkck_ar
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine will check for potentially serious problems.  Note that the upper  !
   ! and lower bound are defined in rk4_coms.f90, so if you need to change any limit for   !
   ! some reason, you can adjust there. (Only exception is veg_temp_min, which is in       !
   ! canopy_radiation_coms).                                                               !
   !---------------------------------------------------------------------------------------!
   subroutine lsm_sanity_check_ar(y,reject_step, csite,ipa, lsl,dydx,h,atm_tmp,atm_shv     &
                                 ,atm_co2,prss,exner,rhos,vels,geoht,pcpg,qpcpg,dpcpg      &
                                 ,print_problems)

      use ed_state_vars         , only : sitetype             & ! structure
                                       , patchtype            & ! structure
                                       , rk4patchtype         & ! structure
                                       , integration_vars_ar  ! ! structure
      use grid_coms             , only : nzg                  ! ! intent(in)
      use soil_coms             , only : soil                 ! ! intent(in), lookup table
      use canopy_radiation_coms , only : lai_min              & ! intent(in)
                                       , veg_temp_min         ! ! intent(in)
      use consts_coms           , only : t3ple                ! ! intent(in)
      use misc_coms             , only : integ_err            & ! intent(inout)
                                       , record_err           ! ! intent(inout)
      use rk4_coms              , only : rk4min_can_temp      & ! intent(in)
                                       , rk4max_can_temp      & ! intent(in)
                                       , rk4max_can_shv       & ! intent(in)
                                       , rk4min_can_shv       & ! intent(in)
                                       , rk4max_veg_temp      & ! intent(in)
                                       , rk4min_veg_water     & ! intent(in)
                                       , rk4min_sfcw_temp     & ! intent(in)
                                       , rk4max_soil_temp     & ! intent(in)
                                       , rk4min_soil_temp     & ! intent(in)
                                       , rk4min_sfcw_mass     & ! intent(in)
                                       , rk4min_virt_water    ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(rk4patchtype) , target      :: y,dydx
      type(sitetype)     , target      :: csite
      integer            , intent(in)  :: lsl
      real               , intent(in)  :: atm_tmp,atm_shv,atm_co2
      real               , intent(in)  :: prss,exner,rhos,vels,geoht,pcpg,qpcpg,dpcpg
      logical            , intent(in)  :: print_problems
      logical            , intent(out) :: reject_step
      !----- Local variables --------------------------------------------------------------!
      type(patchtype)    , pointer     :: cpatch
      integer                          :: k
      real                             :: h
      integer                          :: ipa,ico
      logical                          :: cflag1,cflag2,introuble
      !------------------------------------------------------------------------------------!
      
      !if (print_problems) then
      !   write (unit=*,fmt='(a)') '//////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\'
      !   write (unit=*,fmt='(a)') '            VERBOSE LSM_SANITY_CHECK...             '
      !   write (unit=*,fmt='(a)') '\\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////'
      !   introuble=.false.
      !end if
      
      !----- Assigning initial values for flags. ------------------------------------------!
      cflag1 = .false.
      cflag2 = .false.

      !----- First check, if top soil temperature is NaN, end of story... -----------------!
      if (y%soil_tempk(nzg) /= y%soil_tempk(nzg)) then
         write (unit=*,fmt='(a)') 'Top soil temperature is NaN!!!'
         call print_patch_ar(y,csite,ipa,lsl,atm_tmp,atm_shv,atm_co2,prss,exner,rhos,vels  &
                            ,geoht,pcpg,qpcpg,dpcpg)
      end if

      !----- Being optimistic and assuming things are fine --------------------------------!
      reject_step = .false.


      do k = lsl, nzg
         !----- Checking for water in completely frozen soil ------------------------------!
         if ((y%soil_tempk(k) < (t3ple-0.01) .and. y%soil_fracliq(k) > 0.001) .or.         &
             (y%soil_tempk(k) > (t3ple+0.01) .and. y%soil_fracliq(k) < 0.999) ) then
            reject_step = .true.
            if(record_err) cflag1 = .true.      
            if (print_problems) then
               write(unit=*,fmt='(a)')           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
               write(unit=*,fmt='(a)')           ' + Fracliq doesn''t match temperature...'
               write(unit=*,fmt='(a)')           '----------------------------------------'
               write(unit=*,fmt='(a,1x,i6)')     ' Level:       ',k
               write(unit=*,fmt='(a,1x,es12.5)') ' SOIL_TEMPK:  ',y%soil_tempk(k)
               write(unit=*,fmt='(a,1x,es12.5)') ' SOIL_FLIQ :  ',y%soil_fracliq(k)
               write(unit=*,fmt='(a,1x,es12.5)') ' SOIL_ENERGY: ',y%soil_energy(k)
               write(unit=*,fmt='(a,1x,es12.5)') ' SOIL_WATER:  ',y%soil_water(k)
               write(unit=*,fmt='(a)')           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            elseif (.not. record_err) then
               return
            end if
         end if
         !----- Checking for absurd liquid fraction ---------------------------------------!
         if (y%soil_fracliq(k) > 1.0 .or. y%soil_fracliq(k) < 0.0) then
            reject_step = .true.
            if(record_err) cflag2 = .true.   
            if (print_problems) then
               write(unit=*,fmt='(a)')           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
               write(unit=*,fmt='(a)')           ' + Soil fracliq makes no sense...'
               write(unit=*,fmt='(a)')           '----------------------------------------'
               write(unit=*,fmt='(a,1x,i6)')     ' Level:       ',k
               write(unit=*,fmt='(a,1x,es12.5)') ' SOIL_TEMPK:  ',y%soil_tempk(k)
               write(unit=*,fmt='(a,1x,es12.5)') ' SOIL_FLIQ :  ',y%soil_fracliq(k)
               write(unit=*,fmt='(a,1x,es12.5)') ' SOIL_ENERGY: ',y%soil_energy(k)
               write(unit=*,fmt='(a,1x,es12.5)') ' SOIL_WATER:  ',y%soil_water(k)
               write(unit=*,fmt='(a)')           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            elseif (.not. record_err) then
               return
            end if
         end if
      end do
      if (record_err .and. cflag1) integ_err(42,2) = integ_err(42,2) + 1_8
      if (record_err .and. cflag2) integ_err(43,2) = integ_err(43,2) + 1_8
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !   Checking whether the canopy temperature is too hot or too cold.                  !
      !------------------------------------------------------------------------------------! 
      if (y%can_temp > rk4max_can_temp .or. y%can_temp < rk4min_can_temp) then
         reject_step = .true.
         if(record_err) integ_err(1,2) = integ_err(1,2) + 1_8
         if (print_problems) then
            write(unit=*,fmt='(a)')           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            write(unit=*,fmt='(a)')           ' + Canopy air is off-track...'
            write(unit=*,fmt='(a)')           '----------------------------------------'
            write(unit=*,fmt='(a,1x,es12.5)') ' CAN_TEMP:      ',y%can_temp
            write(unit=*,fmt='(a,1x,es12.5)') ' D(CAN_TEMP)/Dt:',dydx%can_temp
            write(unit=*,fmt='(a,1x,es12.5)') ' H:             ',h
            write(unit=*,fmt='(a)')           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         elseif (.not. record_err) then
            return
         end if
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !   Checking whether the canopy air is too dry or too humid.                         !
      !------------------------------------------------------------------------------------!
      if(y%can_shv > rk4max_can_shv .or. y%can_shv <= rk4min_can_shv)then
         reject_step = .true.
         if(record_err) integ_err(2,2) = integ_err(2,2) + 1_8
         if (print_problems) then
            write(unit=*,fmt='(a)')           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            write(unit=*,fmt='(a)')           ' + Canopy air is off-track...'
            write(unit=*,fmt='(a)')           '----------------------------------------'
            write(unit=*,fmt='(a,1x,es12.5)') ' CAN_SHV:      ',y%can_shv
            write(unit=*,fmt='(a,1x,es12.5)') ' D(CAN_SHV)/Dt:',dydx%can_shv
            write(unit=*,fmt='(a,1x,es12.5)') ' H:            ',h
            write(unit=*,fmt='(a)')           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         elseif (.not. record_err) then
            return
         end if
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Checking whether the temporary snow/water layer(s) has(ve) reasonable values.   !
      !------------------------------------------------------------------------------------!
      k = y%nlev_sfcwater
      if (k >= 1) then
         !----- Temperature ---------------------------------------------------------------!
         if (y%sfcwater_tempk(1) < rk4min_sfcw_temp) then
            reject_step = .true.
            if(record_err) integ_err(44,2) = integ_err(44,2) + 1_8
            if (print_problems) then
               write(unit=*,fmt='(a)')           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
               write(unit=*,fmt='(a)')           ' + Only snow layer is way too cold...'
               write(unit=*,fmt='(a)')           '----------------------------------------'
               write(unit=*,fmt='(a,1x,es12.5)') ' SFCW_TEMP:        ',y%sfcwater_tempk(1)
               write(unit=*,fmt='(a,1x,es12.5)') ' SFCW_ENERGY:      ',y%sfcwater_energy(1)
               write(unit=*,fmt='(a,1x,es12.5)') ' SFCW_MASS:        ',y%sfcwater_mass(1)
               write(unit=*,fmt='(a,1x,es12.5)') ' D(SFCW_ENERGY)/Dt:'                     &
                                                 ,dydx%sfcwater_energy(1)
               write(unit=*,fmt='(a,1x,es12.5)') ' D(SFCW_MASS)/Dt:  '                     &
                                                 ,dydx%sfcwater_mass(1)
               write(unit=*,fmt='(a,1x,es12.5)') ' H:                ',h
               write(unit=*,fmt='(a)')           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            elseif (.not. record_err) then
               return
            end if
         end if
         !----- Mass ----------------------------------------------------------------------!
        if (y%sfcwater_mass(k) < rk4min_sfcw_mass) then
            reject_step = .true.
            if(record_err) integ_err(32+k,2) = integ_err(32+k,2) + 1_8
            if (print_problems) then
               write(unit=*,fmt='(a)')           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
               write(unit=*,fmt='(a)')           ' + Top snow/pond layer has weird mass...'
               write(unit=*,fmt='(a)')           '----------------------------------------'
               write(unit=*,fmt='(a,1x,i6)')     ' # of layer        ',k
               write(unit=*,fmt='(a,1x,es12.5)') ' SFCW_TEMP:        ',y%sfcwater_tempk(k)
               write(unit=*,fmt='(a,1x,es12.5)') ' SFCW_ENERGY:      ',y%sfcwater_energy(k)
               write(unit=*,fmt='(a,1x,es12.5)') ' SFCW_MASS:        ',y%sfcwater_mass(k)
               write(unit=*,fmt='(a,1x,es12.5)') ' D(SFCW_ENERGY)/Dt:'                     &
                                                 ,dydx%sfcwater_energy(k)
               write(unit=*,fmt='(a,1x,es12.5)') ' D(SFCW_MASS)/Dt:  '                     &
                                                 ,dydx%sfcwater_mass(k)
               write(unit=*,fmt='(a,1x,es12.5)') ' H:                ',h
               write(unit=*,fmt='(a)')           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            elseif (.not. record_err) then
               return
            end if
         end if
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Checking whether the soil layers have decent moisture and temperatures.         !
      !------------------------------------------------------------------------------------!
      cflag1 = .false.
      do k=lsl,nzg
         !----- Soil moisture -------------------------------------------------------------!
         if (y%soil_water(k) < 0.95*dble(soil(csite%ntext_soil(k,ipa))%soilcp) .or.        &
             y%soil_water(k) > 1.05*dble(soil(csite%ntext_soil(k,ipa))%slmsts)  ) then
            reject_step = .true.
            if(record_err) integ_err(3+k,2) = integ_err(3+k,2) + 1_8
            if (print_problems) then
               write(unit=*,fmt='(a)')           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
               write(unit=*,fmt='(a)')           ' + Soil layer is off-track...'
               write(unit=*,fmt='(a)')           '----------------------------------------'
               write(unit=*,fmt='(a,1x,i6)')     ' Level:       ',k
               write(unit=*,fmt='(a,1x,es12.5)') ' SOIL_TEMPK:  ',y%soil_tempk(k)
               write(unit=*,fmt='(a,1x,es12.5)') ' SOIL_FLIQ :  ',y%soil_fracliq(k)
               write(unit=*,fmt='(a,1x,es12.5)') ' SOIL_ENERGY: ',y%soil_energy(k)
               write(unit=*,fmt='(a,1x,es12.5)') ' SOIL_WATER:  ',y%soil_water(k)
               write(unit=*,fmt='(a)')           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            elseif (.not. record_err) then
               return
            end if
         end if

         !----- Soil temperature ----------------------------------------------------------!
         if (y%soil_tempk(k) > rk4max_soil_temp .or.                                       &
             y%soil_tempk(k) < rk4min_soil_temp ) then
            reject_step = .true.
            if(record_err) cflag1 = .true.
            if (print_problems) then
               write(unit=*,fmt='(a)')           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
               write(unit=*,fmt='(a)')           ' + Soil layer is off-track...'
               write(unit=*,fmt='(a)')           '----------------------------------------'
               write(unit=*,fmt='(a,1x,i6)')     ' Level:       ',k
               write(unit=*,fmt='(a,1x,es12.5)') ' SOIL_TEMPK:  ',y%soil_tempk(k)
               write(unit=*,fmt='(a,1x,es12.5)') ' SOIL_FLIQ :  ',y%soil_fracliq(k)
               write(unit=*,fmt='(a,1x,es12.5)') ' SOIL_ENERGY: ',y%soil_energy(k)
               write(unit=*,fmt='(a,1x,es12.5)') ' SOIL_WATER:  ',y%soil_water(k)
               write(unit=*,fmt='(a)')           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            elseif (.not. record_err) then
               return
            end if
         end if
      end do
      if(record_err .and. cflag1) integ_err(45,2) = integ_err(45,2) + 1_8
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Negative rk4 increment factors can make this value significantly negative, but !
      ! it should not slow down the integration.
      !------------------------------------------------------------------------------------!
      if (y%virtual_water < rk4min_virt_water) then
         reject_step = .true.
         if(record_err) integ_err(39,2) = integ_err(39,2) + 1_8
         if (print_problems) then
            write(unit=*,fmt='(a)')           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            write(unit=*,fmt='(a)')           ' + Virtual layer mass is off-track...'
            write(unit=*,fmt='(a)')           '----------------------------------------'
            write(unit=*,fmt='(a,1x,i6)')     ' Level:           ',k
            write(unit=*,fmt='(a,1x,es12.5)') ' VIRT_WATER:      ',y%virtual_water
            write(unit=*,fmt='(a,1x,es12.5)') ' VIRT_HEAT :      ',y%virtual_heat
            write(unit=*,fmt='(a,1x,es12.5)') ' D(VIRT_WATER)/Dt:',dydx%virtual_water
            write(unit=*,fmt='(a,1x,es12.5)') ' D(VIRT_HEAT)/Dt :',dydx%virtual_heat
            write(unit=*,fmt='(a)')           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         elseif (.not. record_err) then
            return
         end if
         return
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Checking leaf temperature, but only for those cohorts with sufficient LAI.     !
      !------------------------------------------------------------------------------------!
      cpatch => csite%patch(ipa)
      cflag1 = .false.
      do ico = 1,cpatch%ncohorts
         if (cpatch%lai(ico) > lai_min) then
            if (y%veg_temp(ico) > rk4max_veg_temp .or. y%veg_temp(ico) < veg_temp_min ) then
               reject_step = .true.
               if(record_err) cflag1 = .true.
               if (print_problems) then
                  write(unit=*,fmt='(a)')           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                  write(unit=*,fmt='(a)')           ' + Leaf temperature is off-track...'
                  write(unit=*,fmt='(a)')           '--------------------------------------'
                  write(unit=*,fmt='(a,1x,i6)')     ' PFT:          ',cpatch%pft(ico)
                  write(unit=*,fmt='(a,1x,es12.5)') ' LAI:          ',cpatch%lai(ico)
                  write(unit=*,fmt='(a,1x,es12.5)') ' HCAPVEG:      ',y%hcapveg(ico)
                  write(unit=*,fmt='(a,1x,es12.5)') ' VEG_TEMP:     ',y%veg_temp(ico)
                  write(unit=*,fmt='(a,1x,es12.5)') ' VEG_FRACLIQ:  ',y%veg_fliq(ico)
                  write(unit=*,fmt='(a,1x,es12.5)') ' VEG_ENERGY:   ',y%veg_energy(ico)
                  write(unit=*,fmt='(a,1x,es12.5)') ' VEG_WATER:    ',y%veg_water(ico)
                  write(unit=*,fmt='(a,1x,es12.5)') ' D(VEG_EN)/Dt: ',dydx%veg_energy(ico)
                  write(unit=*,fmt='(a,1x,es12.5)') ' D(VEG_WAT)/Dt:',dydx%veg_water(ico)
                  write(unit=*,fmt='(a,1x,es12.5)') ' H:            ',h
                  write(unit=*,fmt='(a)')           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
               elseif (.not. record_err) then
                  return
               end if
            end if
            if (y%veg_water(ico) < rk4min_veg_water) then
               reject_step = .true.
               if(record_err) cflag1 = .true.
               if (print_problems) then
                  write(unit=*,fmt='(a)')           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                  write(unit=*,fmt='(a)')           ' + Leaf water is off-track...'
                  write(unit=*,fmt='(a)')           '--------------------------------------'
                  write(unit=*,fmt='(a,1x,i6)')     ' PFT:          ',cpatch%pft(ico)
                  write(unit=*,fmt='(a,1x,es12.5)') ' LAI:          ',cpatch%lai(ico)
                  write(unit=*,fmt='(a,1x,es12.5)') ' HCAPVEG:      ',y%hcapveg(ico)
                  write(unit=*,fmt='(a,1x,es12.5)') ' VEG_TEMP:     ',y%veg_temp(ico)
                  write(unit=*,fmt='(a,1x,es12.5)') ' VEG_FRACLIQ:  ',y%veg_fliq(ico)
                  write(unit=*,fmt='(a,1x,es12.5)') ' VEG_ENERGY:   ',y%veg_energy(ico)
                  write(unit=*,fmt='(a,1x,es12.5)') ' VEG_WATER:    ',y%veg_water(ico)
                  write(unit=*,fmt='(a,1x,es12.5)') ' D(VEG_EN)/Dt: ',dydx%veg_energy(ico)
                  write(unit=*,fmt='(a,1x,es12.5)') ' D(VEG_WAT)/Dt:',dydx%veg_water(ico)
                  write(unit=*,fmt='(a,1x,es12.5)') ' H:            ',h
                  write(unit=*,fmt='(a)')           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
               elseif (.not. record_err) then
                  return
               end if
            end if
         end if      

      end do
      if(record_err .and. cflag1) integ_err(46,2) = integ_err(46,2) + 1_8
      !------------------------------------------------------------------------------------!

      !if (print_problems .and. (.not. introuble)) then
      !   write (unit=*,fmt='(a)') '    No serious problem spotted             '
      !   write (unit=*,fmt='(a)') '\\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////'
      !end if


      return
   end subroutine lsm_sanity_check_ar
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This will print the values whenever the step didn't converge due to crazy values  !
   ! no matter how small the steps were (reject_step=.true.).                              !
   !---------------------------------------------------------------------------------------!
   subroutine print_sanity_check_ar(y, csite, ipa, lsl)

      use ed_state_vars         , only : sitetype      & ! structure
                                       , patchtype     & ! structure
                                       , rk4patchtype  ! ! structure
      use grid_coms             , only : nzg           ! ! intent(in)
      use soil_coms             , only : soil          ! ! intent(in), look-up table
      use canopy_radiation_coms , only : lai_min       ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(rk4patchtype) , target     :: y
      type(sitetype)     , target     :: csite
      integer            , intent(in) :: ipa,lsl
      !----- Local variables --------------------------------------------------------------!
      type(patchtype)    , pointer    :: cpatch
      integer                         :: ico,k
      !------------------------------------------------------------------------------------!

      write(unit=*,fmt='(78a)') ('=',k=1,78)
      write(unit=*,fmt='(78a)') ('=',k=1,78)
      write(unit=*,fmt='(a,20x,a,20x,a)') '======','SANITY CHECK','======'
      write(unit=*,fmt='(78a)') ('=',k=1,78)

      write(unit=*,fmt='(a)') ' '
      write(unit=*,fmt='(78a)') ('-',k=1,78)
      write(unit=*,fmt='(a5,3(1x,a12))') 'LEVEL','  SOIL_TEMPK','SOIL_FRACLIQ'             &
                                               &,'  SOIL_WATER'
      do k=lsl,nzg
         write(unit=*,fmt='(i5,3(1x,es12.5))') &
              k, y%soil_tempk(k), y%soil_fracliq(k), y%soil_water(k)
      end do
      write(unit=*,fmt='(78a)') ('-',k=1,78)

      write(unit=*,fmt='(a)') ' '
      write(unit=*,fmt='(78a)') ('-',k=1,78)
      write(unit=*,fmt='(a5,3(1x,a12))') 'LEVEL','  OLD_SOIL_T','OLD_SOIL_FLQ'             &
                                               &,'OLD_SOIL_H2O'
      do k=lsl,nzg
         write(unit=*,fmt='(i5,3(1x,es12.5))')                                             &
              k, csite%soil_tempk(k,ipa), csite%soil_fracliq(k,ipa)                        &
               , csite%soil_water(k,ipa)
      end do
      write(unit=*,fmt='(78a)') ('-',k=1,78)

      write(unit=*,fmt='(a)') ' '
      write(unit=*,fmt='(78a)') ('-',k=1,78)
      write (unit=*,fmt='(a,1x,es12.5)') ' CAN_TEMP=     ',y%can_temp
      write (unit=*,fmt='(a,1x,es12.5)') ' OLD_CAN_TEMP= ',csite%can_temp(ipa)
      write (unit=*,fmt='(a,1x,es12.5)') ' CAN_VAPOR=    ',y%can_shv
      write (unit=*,fmt='(a,1x,es12.5)') ' OLD_CAN_VAP=  ',csite%can_shv(ipa)
      write (unit=*,fmt='(a,1x,i12)')    ' #LEV_SFCH2O=  ',y%nlev_sfcwater
      write (unit=*,fmt='(a,1x,i12)')    ' OLD_#_SFCH2O= ',csite%nlev_sfcwater(ipa)
      if(y%nlev_sfcwater == 1) then
         write(unit=*,fmt='(a,1x,es12.5)') ,'SFCWATER_TEMPK=',y%sfcwater_tempk(1)
      end if
      write(unit=*,fmt='(78a)') ('-',k=1,78)

      write(unit=*,fmt='(a)') ' '
      write(unit=*,fmt='(78a)') ('-',k=1,78)
      cpatch => csite%patch(ipa)
      write (unit=*,fmt='(2(a5,1x),5(a12,1x))')                                            &
         '  COH','  PFT','         LAI','  VEG_ENERGY','OLD_VEG_ENER','    VEG_TEMP'&
         &,'OLD_VEG_TEMP'
      do ico = 1,cpatch%ncohorts
         if(cpatch%lai(ico) > lai_min) then
            write(unit=*,fmt='(2(i5,1x),5(es12.5,1x))')                                    &
               ico,cpatch%pft(ico),cpatch%lai(ico),y%veg_energy(ico)                       &
                  ,cpatch%veg_energy(ico),y%veg_temp(ico),cpatch%veg_temp(ico)
         end if
      end do
      write(unit=*,fmt='(78a)') ('-',k=1,78)

      write(unit=*,fmt='(a)') ' '
      write(unit=*,fmt='(78a)') ('-',k=1,78)
      write (unit=*,fmt='(2(a5,1x),6(a12,1x))') &
         '  COH','  PFT','         LAI','   VEG_WATER',' OLD_VEG_H2O','    HEAT_CAP'&
        &,'RK4_HEAT_CAP' ,'     FRACLIQ'
      do ico = 1,cpatch%ncohorts
         if(cpatch%lai(ico) > lai_min) then
            write(unit=*,fmt='(2(i5,1x),6(es12.5,1x))') &
               ico,cpatch%pft(ico),cpatch%lai(ico),y%veg_water(ico),cpatch%veg_water(ico)  &
                  ,cpatch%hcapveg(ico),y%hcapveg(ico),y%veg_fliq(ico)
         end if
      end do
      write(unit=*,fmt='(78a)') ('-',k=1,78)
      write(unit=*,fmt='(a)') ' '
     
      write(unit=*,fmt='(78a)') ('=',k=1,78)
      write(unit=*,fmt='(78a)') ('=',k=1,78)
      write(unit=*,fmt='(a)') ' '

      return
   end subroutine print_sanity_check_ar
   !=======================================================================================!
   !=======================================================================================!
end module rk4_stepper_ar
!==========================================================================================!
!==========================================================================================!
