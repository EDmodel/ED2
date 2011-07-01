!==========================================================================================!
!==========================================================================================!
!    This module contains some subroutines used to compute the stages for advancing the    !
! Runge-Kutta step.                                                                        !
!------------------------------------------------------------------------------------------!
module rk4_stepper

   contains 
   !=======================================================================================!
   !=======================================================================================!
   !   This subroutine is the main Runge-Kutta step driver.                                !
   !---------------------------------------------------------------------------------------!
   subroutine rkqs(x,htry,hdid,hnext,csite,ipa)

      use rk4_coms      , only : rk4patchtype        & ! structure
                               , integration_buff    & ! intent(inout)
                               , rk4site             & ! intent(in)
                               , hmin                & ! intent(in)
                               , rk4eps              & ! intent(in)
                               , rk4epsi             & ! intent(in)
                               , safety              & ! intent(in)
                               , pgrow               & ! intent(in)
                               , pshrnk              & ! intent(in)
                               , errcon              & ! intent(in)
                               , print_diags         & ! intent(in)
                               , print_detailed      & ! intent(in)
                               , norm_rk4_fluxes     & ! intent(in)
                               , reset_rk4_fluxes    ! ! intent(in)
      use ed_state_vars , only : sitetype            & ! structure
                               , patchtype           ! ! structure
      use grid_coms     , only : nzg                 & ! intent(in)
                               , nzs                 ! ! intent(in)

      implicit none
      !----- Arguments --------------------------------------------------------------------!
      integer                  , intent(in)    :: ipa
      type(sitetype)           , target        :: csite
      real(kind=8)             , intent(in)    :: htry
      real(kind=8)             , intent(inout) :: x
      real(kind=8)             , intent(out)   :: hdid,hnext
      !----- Local variables --------------------------------------------------------------!
      real(kind=8)                             :: h,errmax,xnew,newh,oldh
      logical                                  :: reject_step,reject_result
      logical                                  :: minstep,stuck,test_reject,pdo
      integer                                  :: k
      !------------------------------------------------------------------------------------!

      h           =  htry
      reject_step =  .false.
      hstep:   do

         !---------------------------------------------------------------------------------!
         ! 1. Try a step of varying size.                                                  !
         !---------------------------------------------------------------------------------!
         call rkck(integration_buff%y,integration_buff%dydx,integration_buff%ytemp         &
                     ,integration_buff%yerr,integration_buff%ak2,integration_buff%ak3      &
                     ,integration_buff%ak4,integration_buff%ak5,integration_buff%ak6       &
                     ,integration_buff%ak7,x,h,csite,ipa,reject_step,reject_result)

         !---------------------------------------------------------------------------------!
         ! 2. Check to see how accurate the step was.  Errors were calculated by integrat- !
         !    ing the derivative of that last step.                                        !
         !---------------------------------------------------------------------------------!
         if (reject_step .or. reject_result) then
            !------------------------------------------------------------------------------!
            !    If step was already rejected, that means the step had finished premature- !
            ! ly, so we assign a standard large error (10.0).                              !
            !------------------------------------------------------------------------------!
            errmax = 1.d1
         else
            call get_errmax(errmax, integration_buff%yerr,integration_buff%yscal           &
                           ,csite%patch(ipa),integration_buff%y,integration_buff%ytemp)
            errmax = errmax * rk4epsi
         end if

         !---------------------------------------------------------------------------------!
         ! 3. If that error was large, then calculate a new step size to try.  There are   !
         !    two types of new tries.  If step failed to be minimally reasonable (reject-  !
         !    ed) we have assigned a standard large error (10.0).  Otherwise a new step is !
         !    calculated based on the size of that error.  Hopefully, those new steps      !
         !    should be less than the previous h.  If the error was small, i.e. less then  !
         !    rk4eps, then we are done with this step, and we can move forward             !
         !    time: x = x + h                                                              !
         !---------------------------------------------------------------------------------!
         if (errmax > 1.d0) then
            !----- Defining new step and checking if it can be. ---------------------------!
            oldh    = h
            newh    = safety * h * errmax**pshrnk
            minstep = (newh == h) .or. newh < hmin

            !----- Defining next time, and checking if it really added something. ---------!
            h       = max(1.d-1*h, newh)
            xnew    = x + h
            stuck   = xnew == x

            !------------------------------------------------------------------------------!
            ! 3a. Here is the moment of truth... If we reached a tiny step and yet the     !
            !     model didn't converge, then we print various values to inform the user   !
            !     and abort the run.  Please, don't hate the messenger.                    !
            !------------------------------------------------------------------------------!
            if (minstep .or. stuck) then

               write (unit=*,fmt='(80a)')         ('=',k=1,80)
               write (unit=*,fmt='(a)')           '   STEPSIZE UNDERFLOW IN RKQS'
               write (unit=*,fmt='(80a)')         ('-',k=1,80)
               write (unit=*,fmt='(a,1x,f9.4)')   ' + LONGITUDE:     ',rk4site%lon
               write (unit=*,fmt='(a,1x,f9.4)')   ' + LATITUDE:      ',rk4site%lat
               write (unit=*,fmt='(a)')           ' + PATCH INFO:    '
               write (unit=*,fmt='(a,1x,i6)')     '   - NUMBER:      ',ipa
               write (unit=*,fmt='(a,1x,es12.4)') '   - AGE:         ',csite%age(ipa)
               write (unit=*,fmt='(a,1x,i6)')     '   - DIST_TYPE:   ',csite%dist_type(ipa)
               write (unit=*,fmt='(a,1x,l1)')     ' + REJECT_STEP:   ',reject_step
               write (unit=*,fmt='(a,1x,l1)')     ' + REJECT_RESULT: ',reject_result
               write (unit=*,fmt='(a,1x,l1)')     ' + MINSTEP:       ',minstep
               write (unit=*,fmt='(a,1x,l1)')     ' + STUCK:         ',stuck
               write (unit=*,fmt='(a,1x,es12.4)') ' + ERRMAX:        ',errmax
               write (unit=*,fmt='(a,1x,es12.4)') ' + X:             ',x
               write (unit=*,fmt='(a,1x,es12.4)') ' + H:             ',h
               write (unit=*,fmt='(a,1x,es12.4)') ' + OLDH:          ',oldh
               write (unit=*,fmt='(a,1x,es12.4)') ' + NEWH:          ',newh
               write (unit=*,fmt='(a,1x,es12.4)') ' + SAFETY:        ',safety
               write (unit=*,fmt='(80a)') ('-',k=1,80)
               if (reject_step .or. reject_result) then
                  write (unit=*,fmt='(a)') '   Likely to be a rejected step problem.'
                  write (unit=*,fmt='(80a)') ('=',k=1,80)
               else
                  write (unit=*,fmt='(a)') '   Likely to be an errmax problem.'
                  write (unit=*,fmt='(80a)') ('=',k=1,80)
               end if

               if (reject_result) then
                  !----- Run the LSM sanity check but this time we force the print. -------!
                  call rk4_sanity_check(integration_buff%ytemp,test_reject,csite,ipa       &
                                       ,integration_buff%dydx,h,.true.)
                  call print_sanity_check(integration_buff%y,csite,ipa)
               elseif (reject_step) then
                  call rk4_sanity_check(integration_buff%ak7,test_reject,csite,ipa         &
                                       ,integration_buff%dydx,h,.true.)
                  call print_sanity_check(integration_buff%y,csite,ipa)
               else
                  call print_errmax(errmax,integration_buff%yerr,integration_buff%yscal    &
                                      ,csite%patch(ipa),integration_buff%y                 &
                                      ,integration_buff%ytemp)
                  write (unit=*,fmt='(80a)') ('=',k=1,80)
                  write (unit=*,fmt='(a,1x,es12.4)') ' - Rel. errmax:',errmax
                  write (unit=*,fmt='(a,1x,es12.4)') ' - Raw errmax: ',errmax*rk4eps
                  write (unit=*,fmt='(a,1x,es12.4)') ' - Epsilon:',rk4eps
                  write (unit=*,fmt='(80a)') ('=',k=1,80)
               end if
               call print_rk4patch(integration_buff%y, csite,ipa)
            endif
         
         else

            !------------------------------------------------------------------------------!
            ! 3b.  Great, it worked, so now we can advance to the next step.  We just need !
            !      to do some minor adjustments before...                                  !
            !------------------------------------------------------------------------------!
            !----- i.   Final update of leaf properties to avoid negative water. ----------!
            call adjust_veg_properties(integration_buff%ytemp,h,csite,ipa)
            !----- ii.  Final update of top soil properties to avoid off-bounds moisture. -!
            call adjust_topsoil_properties(integration_buff%ytemp,h,csite,ipa)
            !----- iii. Make temporary surface water stable and positively defined. -------!
            call adjust_sfcw_properties(nzg,nzs,integration_buff%ytemp, h, csite,ipa)
            !----- iv.  Update the diagnostic variables. ----------------------------------!
            call update_diagnostic_vars(integration_buff%ytemp, csite,ipa)
            !------------------------------------------------------------------------------!

            !------------------------------------------------------------------------------!
            ! 3c. Set up h for the next time.  And here we can relax h for the next step,  !
            !    and try something faster.                                                 !
            !------------------------------------------------------------------------------!
            if (errmax > errcon) then
               hnext = safety * h * errmax**pgrow
            else
               hnext = 5.d0 * h
            endif
            hnext = max(2.d0*hmin,hnext)

            !------ 3d. Normalise the fluxes if the user wants detailed debugging. --------!
            if (print_detailed) then
               call norm_rk4_fluxes(integration_buff%ytemp,h)
               call print_rk4_state(integration_buff%y,integration_buff%ytemp,csite,ipa,x,h)
            end if

            !------------------------------------------------------------------------------!
            !    3e. Copy the temporary structure to the intermediate state.               !
            !------------------------------------------------------------------------------!
            call copy_rk4_patch(integration_buff%ytemp,integration_buff%y                  &
                               ,csite%patch(ipa))

            !------------------------------------------------------------------------------!
            !    3f. Flush step-by-step fluxes to zero if the user wants detailed          !
            !        debugging.                                                            !
            !------------------------------------------------------------------------------!
            if (print_detailed) then
               call reset_rk4_fluxes(integration_buff%y)
            end if

            !----- 3g. Update time. -------------------------------------------------------!
            x    = x + h
            hdid = h

            exit hstep
         end if
      end do hstep

      return
   end subroutine rkqs
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine will update the variables and perform the actual time stepping.    !
   !---------------------------------------------------------------------------------------!
   subroutine rkck(y,dydx,yout,yerr,ak2,ak3,ak4,ak5,ak6,ak7,x,h,csite,ipa                  &
                  ,reject_step,reject_result)

      use rk4_coms      , only : rk4patchtype        & ! structure
                               , integration_vars    & ! structure
                               , rk4site             & ! intent(in)
                               , print_diags         & ! intent(in)
                               , rk4_a2              & ! intent(in)
                               , rk4_a3              & ! intent(in)
                               , rk4_a4              & ! intent(in)
                               , rk4_a5              & ! intent(in)
                               , rk4_a6              & ! intent(in)
                               , rk4_b21             & ! intent(in)
                               , rk4_b31             & ! intent(in)
                               , rk4_b32             & ! intent(in)
                               , rk4_b41             & ! intent(in)
                               , rk4_b42             & ! intent(in)
                               , rk4_b43             & ! intent(in)
                               , rk4_b51             & ! intent(in)
                               , rk4_b52             & ! intent(in)
                               , rk4_b53             & ! intent(in)
                               , rk4_b54             & ! intent(in)
                               , rk4_b61             & ! intent(in)
                               , rk4_b62             & ! intent(in)
                               , rk4_b63             & ! intent(in)
                               , rk4_b64             & ! intent(in)
                               , rk4_b65             & ! intent(in)
                               , rk4_c1              & ! intent(in)
                               , rk4_c3              & ! intent(in)
                               , rk4_c4              & ! intent(in)
                               , rk4_c6              & ! intent(in)
                               , rk4_dc5             & ! intent(in)
                               , rk4_dc1             & ! intent(in)
                               , rk4_dc3             & ! intent(in)
                               , rk4_dc4             & ! intent(in)
                               , rk4_dc6             & ! intent(in)
                               , zero_rk4_patch      & ! intent(in)
                               , zero_rk4_cohort     ! ! intent(in)
      use ed_state_vars , only : sitetype            & ! structure
                               , patchtype           ! ! structure
      use grid_coms     , only : nzg                 & ! intent(in)
                               , nzs                 ! ! intent(in)
      implicit none

      !----- Arguments --------------------------------------------------------------------!
      integer           , intent(in)  :: ipa
      real(kind=8)      , intent(in)  :: x,h
      type(rk4patchtype), target      :: y,dydx,yout,yerr
      type(rk4patchtype), target      :: ak2,ak3,ak4,ak5,ak6,ak7
      type(sitetype)    , target      :: csite
      logical           , intent(out) :: reject_step
      logical           , intent(out) :: reject_result
      !----- Local variables --------------------------------------------------------------!
      type(patchtype)   , pointer     :: cpatch
      real(kind=8)                    :: combh
      real(kind=8)                    :: dpdt
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Start and assume that nothing went wrong up to this point... If we find any    !
      ! seriously bad step, quit and reduce the time step without even bothering to try    !
      ! further.                                                                           !
      !------------------------------------------------------------------------------------!
      reject_step   = .false.
      reject_result = .false.
      cpatch => csite%patch(ipa)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     For each stage we checked the sanity, so we avoid using non-sense values to    !
      ! advance.  Also, we estimate the derivative of pressure after each stage, and in    !
      ! the end we estimate the full step derivative as the weighted average for each of   !
      ! these partial steps taken.                                                         !
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Second stage (the first was the Euler, whose derivative is already computed    !
      ! and saved in dydx.                                                                 !
      !------------------------------------------------------------------------------------!
      call copy_rk4_patch(y, ak7, cpatch)
      call inc_rk4_patch(ak7, dydx, rk4_b21*h, cpatch)
      combh = rk4_b21*h
      call update_diagnostic_vars(ak7, csite,ipa)
      call rk4_sanity_check(ak7, reject_step, csite, ipa,dydx,h,print_diags)
      if (reject_step) return
      !------------------------------------------------------------------------------------!


      !------ Estimate the derivative of canopy pressure. ---------------------------------!
      dpdt          = (ak7%can_prss - y%can_prss) / combh
      dydx%can_prss = dpdt * rk4_b21
      !------------------------------------------------------------------------------------!


      !------ Get the new derivative evaluation. ------------------------------------------!
      call leaf_derivs(ak7, ak2, csite, ipa)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Third stage.                                                                    !
      !------------------------------------------------------------------------------------!
      call copy_rk4_patch(y, ak7, cpatch)
      call inc_rk4_patch(ak7, dydx, rk4_b31*h, cpatch)
      call inc_rk4_patch(ak7,  ak2, rk4_b32*h, cpatch)
      combh = (rk4_b31+rk4_b32)*h
      call update_diagnostic_vars(ak7, csite,ipa)
      call rk4_sanity_check(ak7,reject_step,csite,ipa,dydx,h,print_diags)
      if (reject_step) return
      !------------------------------------------------------------------------------------!


      !------ Estimate the derivative of canopy pressure. ---------------------------------!
      dpdt          = (ak7%can_prss - y%can_prss) / combh
      dydx%can_prss = dydx%can_prss + dpdt * rk4_b31
      ak2%can_prss  =                 dpdt * rk4_b32
      !------------------------------------------------------------------------------------!


      !------ Get the new derivative evaluation. ------------------------------------------!
      call leaf_derivs(ak7, ak3, csite,ipa)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Fourth stage.                                                                   !
      !------------------------------------------------------------------------------------!
      call copy_rk4_patch(y, ak7, cpatch)
      call inc_rk4_patch(ak7, dydx, rk4_b41*h, cpatch)
      call inc_rk4_patch(ak7,  ak2, rk4_b42*h, cpatch)
      call inc_rk4_patch(ak7,  ak3, rk4_b43*h, cpatch)
      combh = (rk4_b41+rk4_b42+rk4_b43)*h
      call update_diagnostic_vars(ak7, csite,ipa)
      call rk4_sanity_check(ak7, reject_step, csite,ipa,dydx,h,print_diags)
      if (reject_step) return
      !------------------------------------------------------------------------------------!


      !------ Estimate the derivative of canopy pressure. ---------------------------------!
      dpdt          = (ak7%can_prss - y%can_prss) / combh
      dydx%can_prss = dydx%can_prss + dpdt * rk4_b41
      ak2%can_prss  = ak2%can_prss  + dpdt * rk4_b42
      ak3%can_prss  =                 dpdt * rk4_b43
      !------------------------------------------------------------------------------------!


      !------ Get the new derivative evaluation. ------------------------------------------!
      call leaf_derivs(ak7, ak4, csite, ipa)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Fifth stage.                                                                    !
      !------------------------------------------------------------------------------------!
      call copy_rk4_patch(y, ak7, cpatch)
      call inc_rk4_patch(ak7, dydx, rk4_b51*h, cpatch)
      call inc_rk4_patch(ak7,  ak2, rk4_b52*h, cpatch)
      call inc_rk4_patch(ak7,  ak3, rk4_b53*h, cpatch)
      call inc_rk4_patch(ak7,  ak4, rk4_b54*h, cpatch)
      combh = (rk4_b51+rk4_b52+rk4_b53+rk4_b54)*h
      call update_diagnostic_vars(ak7, csite,ipa)
      call rk4_sanity_check(ak7,reject_step,csite,ipa,dydx,h,print_diags)
      if (reject_step) return
      !------------------------------------------------------------------------------------!


      !------ Estimate the derivative of canopy pressure. ---------------------------------!
      dpdt          = (ak7%can_prss - y%can_prss) / combh
      dydx%can_prss = dydx%can_prss + dpdt * rk4_b51
      ak2%can_prss  = ak2%can_prss  + dpdt * rk4_b52
      ak3%can_prss  = ak3%can_prss  + dpdt * rk4_b53
      ak4%can_prss  =                 dpdt * rk4_b54
      !------------------------------------------------------------------------------------!


      !------ Get the new derivative evaluation. ------------------------------------------!
      call leaf_derivs(ak7, ak5, csite, ipa)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Sixth stage.                                                                    !
      !------------------------------------------------------------------------------------!
      call copy_rk4_patch(y, ak7, cpatch)
      call inc_rk4_patch(ak7, dydx, rk4_b61*h, cpatch)
      call inc_rk4_patch(ak7,  ak2, rk4_b62*h, cpatch)
      call inc_rk4_patch(ak7,  ak3, rk4_b63*h, cpatch)
      call inc_rk4_patch(ak7,  ak4, rk4_b64*h, cpatch)
      call inc_rk4_patch(ak7,  ak5, rk4_b65*h, cpatch)
      combh = (rk4_b61+rk4_b62+rk4_b63+rk4_b64+rk4_b65)*h
      call update_diagnostic_vars(ak7, csite,ipa)
      call rk4_sanity_check(ak7, reject_step, csite,ipa,dydx,h,print_diags)
      if(reject_step)return
      !------------------------------------------------------------------------------------!


      !------ Estimate the derivative of canopy pressure. ---------------------------------!
      dpdt          = (ak7%can_prss - y%can_prss) / combh
      dydx%can_prss = dydx%can_prss + dpdt * rk4_b61
      ak2%can_prss  = ak2%can_prss  + dpdt * rk4_b62
      ak3%can_prss  = ak3%can_prss  + dpdt * rk4_b63
      ak4%can_prss  = ak4%can_prss  + dpdt * rk4_b64
      ak5%can_prss  =                 dpdt * rk4_b65
      !------------------------------------------------------------------------------------!


      !------ Get the new derivative evaluation. ------------------------------------------!
      call leaf_derivs(ak7, ak6, csite,ipa)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Aggregate all derivatives to make the new guess.                                !
      !------------------------------------------------------------------------------------!
      call copy_rk4_patch(y, yout, cpatch)
      call inc_rk4_patch(yout, dydx, rk4_c1*h, cpatch)
      call inc_rk4_patch(yout,  ak3, rk4_c3*h, cpatch)
      call inc_rk4_patch(yout,  ak4, rk4_c4*h, cpatch)
      call inc_rk4_patch(yout,  ak6, rk4_c6*h, cpatch)
      combh = (rk4_c1+rk4_c3+rk4_c4+rk4_c6)*h
      call update_diagnostic_vars   (yout, csite,ipa)
      call rk4_sanity_check(yout, reject_result, csite,ipa,dydx,h,print_diags)
      !------------------------------------------------------------------------------------!
      if(reject_result)return
      !------------------------------------------------------------------------------------!


      !------ Estimate the derivative of canopy pressure. ---------------------------------!
      dpdt          = (ak7%can_prss - y%can_prss) / combh
      dydx%can_prss = dydx%can_prss + dpdt * rk4_c1
      ak3%can_prss  = ak3%can_prss  + dpdt * rk4_c3
      ak4%can_prss  = ak4%can_prss  + dpdt * rk4_c4
      ak6%can_prss  =                 dpdt * rk4_c6
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Average the pressure derivative estimates.                                    !
      !------------------------------------------------------------------------------------!
      dydx%can_prss = dydx%can_prss                                                        &
                    / (rk4_b21 + rk4_b31 + rk4_b41 + rk4_b51 + rk4_b61 + rk4_c1)
      ak2%can_prss  = ak2%can_prss                                                         &
                    / (          rk4_b32 + rk4_b42 + rk4_b52 + rk4_b62         )
      ak3%can_prss  = ak3%can_prss                                                         &
                    / (                    rk4_b43 + rk4_b53 + rk4_b63 + rk4_c3)
      ak4%can_prss  = ak4%can_prss                                                         &
                    / (                              rk4_b54 + rk4_b64 + rk4_c4)
      ak5%can_prss  = ak5%can_prss                                                         &
                    / (                                        rk4_b65         )
      ak6%can_prss  = ak6%can_prss                                                         &
                    / (                                                  rk4_c6)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Estimate the error for this step.                                               !
      !------------------------------------------------------------------------------------!
      call zero_rk4_patch (yerr)
      call zero_rk4_cohort(yerr)
      call inc_rk4_patch(yerr, dydx, rk4_dc1*h, cpatch)
      call inc_rk4_patch(yerr, ak3,  rk4_dc3*h, cpatch)
      call inc_rk4_patch(yerr, ak4,  rk4_dc4*h, cpatch)
      call inc_rk4_patch(yerr, ak5,  rk4_dc5*h, cpatch)
      call inc_rk4_patch(yerr, ak6,  rk4_dc6*h, cpatch)
      !------------------------------------------------------------------------------------!

      return
   end subroutine rkck
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine will check for potentially serious problems.  Note that the upper  !
   ! and lower bound are defined in rk4_coms.f90, so if you need to change any limit for   !
   ! some reason, you can adjust there.                                                    !
   !---------------------------------------------------------------------------------------!
   subroutine rk4_sanity_check(y,reject_step, csite,ipa,dydx,h,print_problems)
      use rk4_coms              , only : rk4patchtype          & ! structure
                                       , integration_vars      & ! structure
                                       , rk4site               & ! intent(in)
                                       , rk4eps                & ! intent(in)
                                       , toocold               & ! intent(in)
                                       , rk4min_can_theiv      & ! intent(in)
                                       , rk4max_can_theiv      & ! intent(in)
                                       , rk4min_can_theta      & ! intent(in)
                                       , rk4max_can_theta      & ! intent(in)
                                       , rk4max_can_shv        & ! intent(in)
                                       , rk4min_can_shv        & ! intent(in)
                                       , rk4min_can_rhv        & ! intent(in)
                                       , rk4max_can_rhv        & ! intent(in)
                                       , rk4min_can_temp       & ! intent(in)
                                       , rk4max_can_temp       & ! intent(in)
                                       , rk4min_can_theiv      & ! intent(in)
                                       , rk4max_can_theiv      & ! intent(in)
                                       , rk4min_can_prss       & ! intent(in)
                                       , rk4max_can_prss       & ! intent(in)
                                       , rk4min_can_co2        & ! intent(in)
                                       , rk4max_can_co2        & ! intent(in)
                                       , rk4max_veg_temp       & ! intent(in)
                                       , rk4min_veg_temp       & ! intent(in)
                                       , rk4min_veg_lwater     & ! intent(in)
                                       , rk4min_sfcw_temp      & ! intent(in)
                                       , rk4max_sfcw_temp      & ! intent(in)
                                       , rk4max_soil_temp      & ! intent(in)
                                       , rk4min_soil_temp      & ! intent(in)
                                       , rk4max_soil_water     & ! intent(in)
                                       , rk4min_soil_water     & ! intent(in)
                                       , rk4min_sfcw_mass      & ! intent(in)
                                       , rk4min_virt_water     & ! intent(in)
                                       , rk4tiny_sfcw_mass     & ! intent(in)
                                       , rk4water_stab_thresh  & ! intent(in)
                                       , integ_err             & ! intent(inout)
                                       , record_err            & ! intent(in)
                                       , osow                  & ! intent(in)
                                       , osoe                  & ! intent(in)
                                       , oswe                  & ! intent(in)
                                       , oswm                  ! ! intent(in)
      use ed_state_vars         , only : sitetype              & ! structure
                                       , patchtype             ! ! structure
      use grid_coms             , only : nzg                   ! ! intent(in)

      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(rk4patchtype) , target      :: y
      type(rk4patchtype) , target      :: dydx
      type(sitetype)     , target      :: csite
      logical            , intent(in)  :: print_problems
      logical            , intent(out) :: reject_step
      real(kind=8)       , intent(in)  :: h
      !----- Local variables --------------------------------------------------------------!
      type(patchtype)    , pointer     :: cpatch
      integer                          :: k
      integer                          :: ksn
      real(kind=8)                     :: rk4min_leaf_water
      real(kind=8)                     :: rk4min_wood_water
      integer                          :: ipa
      integer                          :: ico
      logical                          :: cflag7
      logical                          :: cflag8
      logical                          :: cflag9
      logical                          :: cflag10
      !------------------------------------------------------------------------------------!


      !----- Be optimistic and start assuming that things are fine. -----------------------!
      reject_step = .false.
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !   Check whether the canopy air equivalent potential temperature is off.            !
      !------------------------------------------------------------------------------------!
      if (y%can_theiv > rk4max_can_theiv .or. y%can_theiv < rk4min_can_theiv ) then
         reject_step = .true.
         if(record_err) integ_err(1,2) = integ_err(1,2) + 1_8
         if (print_problems) then
            write(unit=*,fmt='(a)')           '==========================================='
            write(unit=*,fmt='(a)')           ' + Canopy air theta_Eiv is off-track...'
            write(unit=*,fmt='(a)')           '-------------------------------------------'
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_THEIV:         ',y%can_theiv
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_THETA:         ',y%can_theta
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_SHV:           ',y%can_shv
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_RHV:           ',y%can_rhv
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_TEMP:          ',y%can_temp
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_RHOS:          ',y%can_rhos
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_CO2:           ',y%can_co2
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_DEPTH:         ',y%can_depth
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_PRSS:          ',y%can_prss
            write(unit=*,fmt='(a,1x,es12.4)') ' D(CAN_LNTHETA )/Dt:',dydx%can_lntheta
            write(unit=*,fmt='(a,1x,es12.4)') ' D(CAN_SHV     )/Dt:',dydx%can_shv
            write(unit=*,fmt='(a,1x,es12.4)') ' D(CAN_CO2     )/Dt:',dydx%can_co2
            write(unit=*,fmt='(a)')           '==========================================='
            write(unit=*,fmt='(a)')           ' '
         elseif (.not. record_err) then
            return
         end if
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !   Check whether the canopy air potential temperature is off.                       !
      !------------------------------------------------------------------------------------!
      if (y%can_theta > rk4max_can_theta .or. y%can_theta < rk4min_can_theta ) then
         reject_step = .true.
         if(record_err) integ_err(2,2) = integ_err(2,2) + 1_8
         if (print_problems) then
            write(unit=*,fmt='(a)')           '==========================================='
            write(unit=*,fmt='(a)')           ' + Canopy air pot. temp. is off-track...'
            write(unit=*,fmt='(a)')           '-------------------------------------------'
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_THEIV:         ',y%can_theiv
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_THETA:         ',y%can_theta
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_SHV:           ',y%can_shv
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_RHV:           ',y%can_rhv
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_TEMP:          ',y%can_temp
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_RHOS:          ',y%can_rhos
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_CO2:           ',y%can_co2
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_DEPTH:         ',y%can_depth
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_PRSS:          ',y%can_prss
            write(unit=*,fmt='(a,1x,es12.4)') ' D(CAN_LNTHETA )/Dt:',dydx%can_lntheta
            write(unit=*,fmt='(a,1x,es12.4)') ' D(CAN_SHV     )/Dt:',dydx%can_shv
            write(unit=*,fmt='(a,1x,es12.4)') ' D(CAN_CO2     )/Dt:',dydx%can_co2
            write(unit=*,fmt='(a)')           '==========================================='
            write(unit=*,fmt='(a)')           ' '
         elseif (.not. record_err) then
            return
         end if
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !   Check whether the canopy air equivalent potential temperature is off.            !
      !------------------------------------------------------------------------------------!
      if ( y%can_shv > rk4max_can_shv .or. y%can_shv < rk4min_can_shv  ) then
         reject_step = .true.
         if(record_err) integ_err(3,2) = integ_err(3,2) + 1_8
         if (print_problems) then
            write(unit=*,fmt='(a)')           '==========================================='
            write(unit=*,fmt='(a)')           ' + Canopy air sp. humidity is off-track...'
            write(unit=*,fmt='(a)')           '-------------------------------------------'
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_THEIV:         ',y%can_theiv
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_THETA:         ',y%can_theta
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_SHV:           ',y%can_shv
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_RHV:           ',y%can_rhv
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_TEMP:          ',y%can_temp
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_RHOS:          ',y%can_rhos
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_CO2:           ',y%can_co2
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_DEPTH:         ',y%can_depth
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_PRSS:          ',y%can_prss
            write(unit=*,fmt='(a,1x,es12.4)') ' D(CAN_LNTHETA )/Dt:',dydx%can_lntheta
            write(unit=*,fmt='(a,1x,es12.4)') ' D(CAN_SHV     )/Dt:',dydx%can_shv
            write(unit=*,fmt='(a,1x,es12.4)') ' D(CAN_CO2     )/Dt:',dydx%can_co2
            write(unit=*,fmt='(a)')           '==========================================='
            write(unit=*,fmt='(a)')           ' '
         elseif (.not. record_err) then
            return
         end if
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !   Check whether the canopy air equivalent potential temperature is off.            !
      !------------------------------------------------------------------------------------!
      if (y%can_temp > rk4max_can_temp .or. y%can_temp < rk4min_can_temp) then
         reject_step = .true.
         if(record_err) integ_err(4,2) = integ_err(4,2) + 1_8
         if (print_problems) then
            write(unit=*,fmt='(a)')           '==========================================='
            write(unit=*,fmt='(a)')           ' + Canopy air temperature is off-track...'
            write(unit=*,fmt='(a)')           '-------------------------------------------'
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_THEIV:         ',y%can_theiv
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_THETA:         ',y%can_theta
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_SHV:           ',y%can_shv
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_RHV:           ',y%can_rhv
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_TEMP:          ',y%can_temp
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_RHOS:          ',y%can_rhos
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_CO2:           ',y%can_co2
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_DEPTH:         ',y%can_depth
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_PRSS:          ',y%can_prss
            write(unit=*,fmt='(a,1x,es12.4)') ' D(CAN_LNTHETA )/Dt:',dydx%can_lntheta
            write(unit=*,fmt='(a,1x,es12.4)') ' D(CAN_SHV     )/Dt:',dydx%can_shv
            write(unit=*,fmt='(a,1x,es12.4)') ' D(CAN_CO2     )/Dt:',dydx%can_co2
            write(unit=*,fmt='(a)')           '==========================================='
            write(unit=*,fmt='(a)')           ' '
         elseif (.not. record_err) then
            return
         end if
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !   Check whether the canopy air equivalent potential temperature is off.            !
      !------------------------------------------------------------------------------------!
      if (y%can_prss > rk4max_can_prss .or. y%can_prss < rk4min_can_prss) then
         reject_step = .true.
         if(record_err) integ_err(5,2) = integ_err(5,2) + 1_8
         if (print_problems) then
            write(unit=*,fmt='(a)')           '==========================================='
            write(unit=*,fmt='(a)')           ' + Canopy air pressure is off-track...'
            write(unit=*,fmt='(a)')           '-------------------------------------------'
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_THEIV:         ',y%can_theiv
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_THETA:         ',y%can_theta
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_SHV:           ',y%can_shv
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_RHV:           ',y%can_rhv
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_TEMP:          ',y%can_temp
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_RHOS:          ',y%can_rhos
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_CO2:           ',y%can_co2
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_DEPTH:         ',y%can_depth
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_PRSS:          ',y%can_prss
            write(unit=*,fmt='(a,1x,es12.4)') ' D(CAN_LNTHETA )/Dt:',dydx%can_lntheta
            write(unit=*,fmt='(a,1x,es12.4)') ' D(CAN_SHV     )/Dt:',dydx%can_shv
            write(unit=*,fmt='(a,1x,es12.4)') ' D(CAN_CO2     )/Dt:',dydx%can_co2
            write(unit=*,fmt='(a)')           '==========================================='
            write(unit=*,fmt='(a)')           ' '
         elseif (.not. record_err) then
            return
         end if
      end if
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !   Check whether the canopy air equivalent potential temperature is off.            !
      !------------------------------------------------------------------------------------!
      if (y%can_co2 > rk4max_can_co2 .or. y%can_co2 < rk4min_can_co2) then
         reject_step = .true.
         if(record_err) integ_err(6,2) = integ_err(6,2) + 1_8
         if (print_problems) then
            write(unit=*,fmt='(a)')           '==========================================='
            write(unit=*,fmt='(a)')           ' + Canopy air CO2 is off-track...'
            write(unit=*,fmt='(a)')           '-------------------------------------------'
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_THEIV:         ',y%can_theiv
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_THETA:         ',y%can_theta
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_SHV:           ',y%can_shv
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_RHV:           ',y%can_rhv
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_TEMP:          ',y%can_temp
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_RHOS:          ',y%can_rhos
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_CO2:           ',y%can_co2
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_DEPTH:         ',y%can_depth
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_PRSS:          ',y%can_prss
            write(unit=*,fmt='(a,1x,es12.4)') ' D(CAN_LNTHETA )/Dt:',dydx%can_lntheta
            write(unit=*,fmt='(a,1x,es12.4)') ' D(CAN_SHV     )/Dt:',dydx%can_shv
            write(unit=*,fmt='(a,1x,es12.4)') ' D(CAN_CO2     )/Dt:',dydx%can_co2
            write(unit=*,fmt='(a)')           '==========================================='
            write(unit=*,fmt='(a)')           ' '
         elseif (.not. record_err) then
            return
         end if
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Check leaf properties, but only for those cohorts with sufficient LAI.         !
      !------------------------------------------------------------------------------------!
      cpatch => csite%patch(ipa)
      cflag7 = .false.
      cflag8 = .false.
      leafloop: do ico = 1,cpatch%ncohorts
         if (.not. y%leaf_resolvable(ico)) cycle leafloop

         !----- Find the minimum leaf surface water. --------------------------------------!
         rk4min_leaf_water = rk4min_veg_lwater * y%lai(ico)

         !----- Check leaf surface water. -------------------------------------------------!
         if (y%leaf_water(ico) < rk4min_leaf_water) then
            reject_step = .true.
            if(record_err) cflag7 = .true.
            if (print_problems) then
               write(unit=*,fmt='(a)')           '========================================'
               write(unit=*,fmt='(a)')           ' + Leaf surface water is off-track...'
               write(unit=*,fmt='(a)')           '========================================'
               write(unit=*,fmt='(a,1x,i6)')     ' PFT:           ',cpatch%pft(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' HEIGHT:        ',cpatch%hite(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' LAI:           ',y%lai(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' WAI:           ',y%wai(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' WPA:           ',y%wpa(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' TAI:           ',y%tai(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' NPLANT:        ',y%nplant(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' CROWN_AREA:    ',y%crown_area(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' LEAF_HCAP:     ',y%leaf_hcap(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' LEAF_TEMP:     ',y%leaf_temp(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' LEAF_FRACLIQ:  ',y%leaf_fliq(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' LEAF_ENERGY:   ',y%leaf_energy(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' LEAF_WATER:    ',y%leaf_water(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' VEG_WIND:      ',y%veg_wind(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' LINT_SHV:      ',y%lint_shv(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' MIN_LEAF_WATER:',rk4min_leaf_water
               write(unit=*,fmt='(a,1x,es12.4)') ' LEAF_GBH:      ',y%leaf_gbh(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' LEAF_GBW:      ',y%leaf_gbw(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' LEAF_REYNOLDS: ',y%leaf_reynolds(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' LEAF_GRASHOF:  ',y%leaf_grashof(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' LEAF_NUFREE:   ',y%leaf_nussfree(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' LEAF_NUFORC:   ',y%leaf_nussforc(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' D(LEAF_EN)/Dt: ',dydx%leaf_energy(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' D(LEAF_WAT)/Dt:',dydx%leaf_water(ico)
               write(unit=*,fmt='(a)')           '========================================'
            elseif (.not. record_err) then
               return
            end if
         end if

         !----- Check leaf temperature. ---------------------------------------------------!
         if (y%leaf_temp(ico) > rk4max_veg_temp .or.                                       &
             y%leaf_temp(ico) < rk4min_veg_temp      ) then
            reject_step = .true.
            if(record_err) cflag8 = .true.
            if (print_problems) then
               write(unit=*,fmt='(a)')           '========================================'
               write(unit=*,fmt='(a)')           ' + Leaf temperature is off-track...'
               write(unit=*,fmt='(a)')           '========================================'
               write(unit=*,fmt='(a,1x,i6)')     ' PFT:           ',cpatch%pft(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' HEIGHT:        ',cpatch%hite(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' LAI:           ',y%lai(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' WAI:           ',y%wai(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' WPA:           ',y%wpa(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' TAI:           ',y%tai(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' NPLANT:        ',y%nplant(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' CROWN_AREA:    ',y%crown_area(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' LEAF_HCAP:     ',y%leaf_hcap(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' LEAF_TEMP:     ',y%leaf_temp(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' LEAF_FRACLIQ:  ',y%leaf_fliq(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' LEAF_ENERGY:   ',y%leaf_energy(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' LEAF_WATER:    ',y%leaf_water(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' VEG_WIND:      ',y%veg_wind(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' LINT_SHV:      ',y%lint_shv(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' MIN_LEAF_WATER:',rk4min_leaf_water
               write(unit=*,fmt='(a,1x,es12.4)') ' LEAF_GBH:      ',y%leaf_gbh(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' LEAF_GBW:      ',y%leaf_gbw(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' LEAF_REYNOLDS: ',y%leaf_reynolds(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' LEAF_GRASHOF:  ',y%leaf_grashof(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' LEAF_NUFREE:   ',y%leaf_nussfree(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' LEAF_NUFORC:   ',y%leaf_nussforc(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' D(LEAF_EN)/Dt: ',dydx%leaf_energy(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' D(LEAF_WAT)/Dt:',dydx%leaf_water(ico)
               write(unit=*,fmt='(a)')           '========================================'
            elseif (.not. record_err) then
               return
            end if
         end if
      end do leafloop
      if(record_err .and. cflag7) integ_err(7,2) = integ_err(7,2) + 1_8
      if(record_err .and. cflag8) integ_err(8,2) = integ_err(8,2) + 1_8
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Check wood properties, but only for those cohorts with sufficient LAI.         !
      !------------------------------------------------------------------------------------!
      cpatch => csite%patch(ipa)
      cflag9  = .false.
      cflag10 = .false.
      woodloop: do ico = 1,cpatch%ncohorts
         if (.not. y%wood_resolvable(ico)) cycle woodloop

         !----- Find the minimum wood surface water. --------------------------------------!
         rk4min_wood_water = rk4min_veg_lwater * y%wai(ico)

         !----- Check wood surface water. -------------------------------------------------!
         if (y%wood_water(ico) < rk4min_wood_water) then
            reject_step = .true.
            if(record_err) cflag9 = .true.
            if (print_problems) then
               write(unit=*,fmt='(a)')           '========================================'
               write(unit=*,fmt='(a)')           ' + Wood surface water is off-track...'
               write(unit=*,fmt='(a)')           '========================================'
               write(unit=*,fmt='(a,1x,i6)')     ' PFT:           ',cpatch%pft(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' HEIGHT:        ',cpatch%hite(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' LAI:           ',y%lai(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' WAI:           ',y%wai(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' WPA:           ',y%wpa(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' TAI:           ',y%tai(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' NPLANT:        ',y%nplant(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' CROWN_AREA:    ',y%crown_area(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' WOOD_HCAP:     ',y%wood_hcap(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' WOOD_TEMP:     ',y%wood_temp(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' WOOD_FRACLIQ:  ',y%wood_fliq(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' WOOD_ENERGY:   ',y%wood_energy(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' WOOD_WATER:    ',y%wood_water(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' VEG_WIND:      ',y%veg_wind(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' LINT_SHV:      ',y%lint_shv(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' MIN_WOOD_WATER:',rk4min_wood_water
               write(unit=*,fmt='(a,1x,es12.4)') ' WOOD_GBH:      ',y%wood_gbh(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' WOOD_GBW:      ',y%wood_gbw(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' WOOD_REYNOLDS: ',y%wood_reynolds(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' WOOD_GRASHOF:  ',y%wood_grashof(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' WOOD_NUFREE:   ',y%wood_nussfree(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' WOOD_NUFORC:   ',y%wood_nussforc(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' D(WOOD_EN)/Dt: ',dydx%wood_energy(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' D(WOOD_WAT)/Dt:',dydx%wood_water(ico)
               write(unit=*,fmt='(a)')           '========================================'
            elseif (.not. record_err) then
               return
            end if
         end if

         !----- Check wood temperature. ---------------------------------------------------!
         if (y%wood_temp(ico) > rk4max_veg_temp .or.                                       &
             y%wood_temp(ico) < rk4min_veg_temp      ) then
            reject_step = .true.
            if(record_err) cflag10 = .true.
            if (print_problems) then
               write(unit=*,fmt='(a)')           '========================================'
               write(unit=*,fmt='(a)')           ' + Wood temperature is off-track...'
               write(unit=*,fmt='(a)')           '========================================'
               write(unit=*,fmt='(a,1x,i6)')     ' PFT:           ',cpatch%pft(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' HEIGHT:        ',cpatch%hite(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' LAI:           ',y%lai(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' WAI:           ',y%wai(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' WPA:           ',y%wpa(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' TAI:           ',y%tai(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' NPLANT:        ',y%nplant(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' CROWN_AREA:    ',y%crown_area(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' WOOD_HCAP:     ',y%wood_hcap(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' WOOD_TEMP:     ',y%wood_temp(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' WOOD_FRACLIQ:  ',y%wood_fliq(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' WOOD_ENERGY:   ',y%wood_energy(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' WOOD_WATER:    ',y%wood_water(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' VEG_WIND:      ',y%veg_wind(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' LINT_SHV:      ',y%lint_shv(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' MIN_WOOD_WATER:',rk4min_wood_water
               write(unit=*,fmt='(a,1x,es12.4)') ' WOOD_GBH:      ',y%wood_gbh(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' WOOD_GBW:      ',y%wood_gbw(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' WOOD_REYNOLDS: ',y%wood_reynolds(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' WOOD_GRASHOF:  ',y%wood_grashof(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' WOOD_NUFREE:   ',y%wood_nussfree(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' WOOD_NUFORC:   ',y%wood_nussforc(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' D(WOOD_EN)/Dt: ',dydx%wood_energy(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' D(WOOD_WAT)/Dt:',dydx%wood_water(ico)
               write(unit=*,fmt='(a)')           '========================================'
            elseif (.not. record_err) then
               return
            end if
         end if
      end do woodloop
      if(record_err .and. cflag9 ) integ_err( 9,2) = integ_err( 9,2) + 1_8
      if(record_err .and. cflag10) integ_err(10,2) = integ_err(10,2) + 1_8
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !     Check the water mass of the virtual pool.  The energy is checked only when     !
      ! there is enough mass.                                                              !
      !------------------------------------------------------------------------------------!
      if (y%virtual_water < rk4min_virt_water) then
         reject_step = .true.
         if(record_err) integ_err(12,2) = integ_err(12,2) + 1_8
         if (print_problems) then
            write(unit=*,fmt='(a)')           '==========================================='
            write(unit=*,fmt='(a)')           ' + Virtual layer mass is off-track...'
            write(unit=*,fmt='(a)')           '-------------------------------------------'
            write(unit=*,fmt='(a,1x,es12.4)') ' VIRTUAL_ENERGY:   ',y%virtual_energy
            write(unit=*,fmt='(a,1x,es12.4)') ' VIRTUAL_WATER:    ',y%virtual_water
            write(unit=*,fmt='(a,1x,es12.4)') ' VIRTUAL_DEPTH:    ',y%virtual_depth
            write(unit=*,fmt='(a,1x,es12.4)') ' VIRTUAL_TEMPK:    ',y%virtual_tempk
            write(unit=*,fmt='(a,1x,es12.4)') ' VIRTUAL_FLIQ :    ',y%virtual_fracliq
            write(unit=*,fmt='(a,1x,es12.4)') ' D(VIRT_WATER)/Dt: ',dydx%virtual_water
            write(unit=*,fmt='(a,1x,es12.4)') ' D(VIRT_ENERGY)/Dt:',dydx%virtual_energy
            write(unit=*,fmt='(a)')           '==========================================='
         elseif (.not. record_err) then
            return
         end if
      elseif (y%virtual_water > 5.d-1 * rk4water_stab_thresh .and.                         &
             (y%virtual_tempk < rk4min_sfcw_temp .or. y%virtual_tempk > rk4max_sfcw_temp)) &
      then
         reject_step = .true.
         if(record_err) integ_err(11,2) = integ_err(11,2) + 1_8
         if (print_problems) then
            write(unit=*,fmt='(a)')           '==========================================='
            write(unit=*,fmt='(a)')           ' + Virtual layer temp. is off-track...'
            write(unit=*,fmt='(a)')           '-------------------------------------------'
            write(unit=*,fmt='(a,1x,es12.4)') ' VIRTUAL_ENERGY:   ',y%virtual_energy
            write(unit=*,fmt='(a,1x,es12.4)') ' VIRTUAL_WATER:    ',y%virtual_water
            write(unit=*,fmt='(a,1x,es12.4)') ' VIRTUAL_DEPTH:    ',y%virtual_depth
            write(unit=*,fmt='(a,1x,es12.4)') ' VIRTUAL_TEMPK:    ',y%virtual_tempk
            write(unit=*,fmt='(a,1x,es12.4)') ' VIRTUAL_FLIQ :    ',y%virtual_fracliq
            write(unit=*,fmt='(a,1x,es12.4)') ' D(VIRT_WATER)/Dt: ',dydx%virtual_water
            write(unit=*,fmt='(a,1x,es12.4)') ' D(VIRT_ENERGY)/Dt:',dydx%virtual_energy
            write(unit=*,fmt='(a)')           '==========================================='
         elseif (.not. record_err) then
            return
         end if
         return
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Checking whether the soil layers have decent moisture and temperatures.         !
      !------------------------------------------------------------------------------------!
      do k=rk4site%lsl,nzg
         !----- Soil moisture -------------------------------------------------------------!
         if (y%soil_water(k)< rk4min_soil_water(k) .or.                                    &
             y%soil_water(k)> rk4max_soil_water(k) ) then
            reject_step = .true.
            if(record_err) integ_err(osow+k,2) = integ_err(osow+k,2) + 1_8
            if (print_problems) then
               write(unit=*,fmt='(a)')           '========================================'
               write(unit=*,fmt='(a)')           ' + Soil layer water is off-track...'
               write(unit=*,fmt='(a)')           '----------------------------------------'
               write(unit=*,fmt='(a,1x,i6)')     ' Level:       ',k
               write(unit=*,fmt='(a,1x,es12.4)') ' H:           ',h
               write(unit=*,fmt='(a,1x,es12.4)') ' SOIL_TEMPK:  ',y%soil_tempk(k)
               write(unit=*,fmt='(a,1x,es12.4)') ' SOIL_FLIQ :  ',y%soil_fracliq(k)
               write(unit=*,fmt='(a,1x,es12.4)') ' SOIL_ENERGY: ',y%soil_energy(k)
               write(unit=*,fmt='(a,1x,es12.4)') ' SOIL_WATER:  ',y%soil_water(k)
               write(unit=*,fmt='(a,1x,es12.4)') ' D(SOIL_E)/Dt:',dydx%soil_energy(k)
               write(unit=*,fmt='(a,1x,es12.4)') ' D(SOIL_M)/Dt:',dydx%soil_water(k)
               if (k == nzg .and. y%nlev_sfcwater > 0) then
                  write(unit=*,fmt='(a,1x,es12.4)') ' SFCW_TEMP:   ',y%sfcwater_tempk(1)
                  write(unit=*,fmt='(a,1x,es12.4)') ' SFCW_ENERGY: ',y%sfcwater_energy(1)
                  write(unit=*,fmt='(a,1x,es12.4)') ' SFCW_MASS:   ',y%sfcwater_mass(1)
                  write(unit=*,fmt='(a,1x,es12.4)') ' SFCW_DEPTH:  ',y%sfcwater_depth(1)
                  write(unit=*,fmt='(a,1x,es12.4)') ' D(SFCW_E)/Dt:',dydx%sfcwater_energy(1)
                  write(unit=*,fmt='(a,1x,es12.4)') ' D(SFCW_M)/Dt:',dydx%sfcwater_mass(1)
               end if
               write(unit=*,fmt='(a)')           '========================================'
            elseif (.not. record_err) then
               return
            end if
         end if

         !----- Soil temperature ----------------------------------------------------------!
         if (y%soil_tempk(k) > rk4max_soil_temp .or. y%soil_tempk(k) < rk4min_soil_temp )  &
         then
            reject_step = .true.
            if(record_err) integ_err(osoe+k,2) = integ_err(osoe+k,2) + 1_8
            if (print_problems) then
               write(unit=*,fmt='(a)')           '========================================'
               write(unit=*,fmt='(a)')           ' + Soil layer temp is off-track...'
               write(unit=*,fmt='(a)')           '----------------------------------------'
               write(unit=*,fmt='(a,1x,i6)')     ' Level:       ',k
               write(unit=*,fmt='(a,1x,es12.4)') ' H:           ',h
               write(unit=*,fmt='(a,1x,es12.4)') ' SOIL_TEMPK:  ',y%soil_tempk(k)
               write(unit=*,fmt='(a,1x,es12.4)') ' SOIL_FLIQ :  ',y%soil_fracliq(k)
               write(unit=*,fmt='(a,1x,es12.4)') ' SOIL_ENERGY: ',y%soil_energy(k)
               write(unit=*,fmt='(a,1x,es12.4)') ' SOIL_WATER:  ',y%soil_water(k)
               write(unit=*,fmt='(a,1x,es12.4)') ' D(SOIL_E)/Dt:',dydx%soil_energy(k)
               write(unit=*,fmt='(a,1x,es12.4)') ' D(SOIL_M)/Dt:',dydx%soil_water(k)
               if (k == nzg .and. y%nlev_sfcwater > 0) then
                  write(unit=*,fmt='(a,1x,es12.4)') ' SFCW_TEMP:   ',y%sfcwater_tempk(1)
                  write(unit=*,fmt='(a,1x,es12.4)') ' SFCW_ENERGY: ',y%sfcwater_energy(1)
                  write(unit=*,fmt='(a,1x,es12.4)') ' SFCW_MASS:   ',y%sfcwater_mass(1)
                  write(unit=*,fmt='(a,1x,es12.4)') ' SFCW_DEPTH:  ',y%sfcwater_depth(1)
                  write(unit=*,fmt='(a,1x,es12.4)') ' D(SFCW_E)/Dt:',dydx%sfcwater_energy(1)
                  write(unit=*,fmt='(a,1x,es12.4)') ' D(SFCW_M)/Dt:',dydx%sfcwater_mass(1)
               end if
               write(unit=*,fmt='(a)')           '========================================'
            elseif (.not. record_err) then
               return
            end if
         end if
      end do
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !    Check whether the temporary snow/water layer(s) has(ve) reasonable values.      !
      !------------------------------------------------------------------------------------!
      ksn = y%nlev_sfcwater

      do k=1, ksn
         !----- Temperature ---------------------------------------------------------------!
         if (y%sfcwater_tempk(k) < rk4min_sfcw_temp .or.                                   &
             y%sfcwater_tempk(k) > rk4max_sfcw_temp      ) then
            reject_step = .true.
            if(record_err) integ_err(oswe+ksn,2) = integ_err(oswe+ksn,2) + 1_8
            if (print_problems) then
               write(unit=*,fmt='(a)')           '========================================'
               write(unit=*,fmt='(a)')           ' + Snow/pond temperature is off...'
               write(unit=*,fmt='(a)')           '----------------------------------------'
               write(unit=*,fmt='(a,1x,i6)')     ' This layer:    ',k
               write(unit=*,fmt='(a,1x,i6)')     ' # of layers:   ',y%nlev_sfcwater
               write(unit=*,fmt='(a,1x,i6)')     ' Stability flag:',y%flag_sfcwater
               write(unit=*,fmt='(a,1x,es12.4)') ' SFCW_TEMP:     ',y%sfcwater_tempk(k)
               write(unit=*,fmt='(a,1x,es12.4)') ' SFCW_ENERGY:   ',y%sfcwater_energy(k)
               write(unit=*,fmt='(a,1x,es12.4)') ' SFCW_MASS:     ',y%sfcwater_mass(k)
               write(unit=*,fmt='(a,1x,es12.4)') ' SFCW_DEPTH:    ',y%sfcwater_depth(k)
               write(unit=*,fmt='(a,1x,es12.4)') ' D(SFCW_E)/Dt:  ',dydx%sfcwater_energy(k)
               write(unit=*,fmt='(a,1x,es12.4)') ' D(SFCW_M)/Dt:  ',dydx%sfcwater_mass(k)
               write(unit=*,fmt='(a)')           '========================================'
            elseif (.not. record_err) then
               return
            end if
         end if

         !----- Mass ----------------------------------------------------------------------!
         if (y%sfcwater_mass(k) < rk4min_sfcw_mass) then
            reject_step = .true.
            if(record_err) integ_err(oswm+ksn,2) = integ_err(oswm+ksn,2) + 1_8
            if (print_problems) then
               write(unit=*,fmt='(a)')           '========================================'
               write(unit=*,fmt='(a)')           ' + Snow/pond mass is off...'
               write(unit=*,fmt='(a)')           '----------------------------------------'
               write(unit=*,fmt='(a,1x,i6)')     ' This layer:    ',k
               write(unit=*,fmt='(a,1x,i6)')     ' # of layers:   ',y%nlev_sfcwater
               write(unit=*,fmt='(a,1x,i6)')     ' Stability flag:',y%flag_sfcwater
               write(unit=*,fmt='(a,1x,es12.4)') ' SFCW_TEMP:     ',y%sfcwater_tempk(k)
               write(unit=*,fmt='(a,1x,es12.4)') ' SFCW_ENERGY:   ',y%sfcwater_energy(k)
               write(unit=*,fmt='(a,1x,es12.4)') ' SFCW_MASS:     ',y%sfcwater_mass(k)
               write(unit=*,fmt='(a,1x,es12.4)') ' SFCW_DEPTH:    ',y%sfcwater_depth(k)
               write(unit=*,fmt='(a,1x,es12.4)') ' D(SFCW_E)/Dt:  ',dydx%sfcwater_energy(k)
               write(unit=*,fmt='(a,1x,es12.4)') ' D(SFCW_M)/Dt:  ',dydx%sfcwater_mass(k)
               write(unit=*,fmt='(a)')           '========================================'
            elseif (.not. record_err) then
               return
            end if
         end if
      end do
      !------------------------------------------------------------------------------------!



      if (reject_step .and. print_problems) then
         write(unit=*,fmt='(a)')           ' '
         write(unit=*,fmt='(78a)')         ('=',k=1,78)
         write(unit=*,fmt='(a,1x,es12.4)') ' TIMESTEP:          ',h
         write(unit=*,fmt='(a)')           ' '
         write(unit=*,fmt='(a)')           '         ---- SANITY CHECK BOUNDS ----'
         write(unit=*,fmt='(a)')           ' '
         write(unit=*,fmt='(a)')           ' 1. CANOPY AIR SPACE: '
         write(unit=*,fmt='(a)')           ' '
         write(unit=*,fmt='(6(a,1x))')     '   MIN_THEIV','   MAX_THEIV','     MIN_SHV'    &
                                          ,'     MAX_SHV','     MIN_RHV','     MAX_RHV'
         write(unit=*,fmt='(6(es12.5,1x))')  rk4min_can_theiv,rk4max_can_theiv             &
                                            ,rk4min_can_shv  ,rk4max_can_shv               &
                                            ,rk4min_can_rhv  ,rk4max_can_rhv
         write(unit=*,fmt='(a)') ' '
         write(unit=*,fmt='(4(a,1x))')     '    MIN_TEMP','    MAX_TEMP','   MIN_THETA'    &
                                          ,'   MAX_THETA'
         write(unit=*,fmt='(4(es12.5,1x))') rk4min_can_temp ,rk4max_can_temp               &
                                           ,rk4min_can_theta,rk4max_can_theta
         write(unit=*,fmt='(a)') ' '
         write(unit=*,fmt='(4(a,1x))')     '    MIN_PRSS','    MAX_PRSS','     MIN_CO2'    &
                                          ,'     MAX_CO2'
         write(unit=*,fmt='(4(es12.5,1x))') rk4min_can_prss ,rk4max_can_prss               &
                                           ,rk4min_can_co2  ,rk4max_can_co2
         write(unit=*,fmt='(a)') ' '
         write(unit=*,fmt='(78a)')         ('-',k=1,78)
         write(unit=*,fmt='(a)')           ' '
         write(unit=*,fmt='(a)')           ' 2. LEAF PROPERTIES: '
         write(unit=*,fmt='(3(a,1x))')     '    MIN_TEMP','    MAX_TEMP','  MIN_LWATER'
         write(unit=*,fmt='(3(es12.5,1x))') rk4min_veg_temp ,rk4max_veg_temp               &
                                           ,rk4min_veg_lwater
         write(unit=*,fmt='(a)')           ' '
         write(unit=*,fmt='(78a)')         ('-',k=1,78)
         write(unit=*,fmt='(a)')           ' '
         write(unit=*,fmt='(a)')           ' 3. SURFACE WATER / VIRTUAL POOL PROPERTIES: '
         write(unit=*,fmt='(3(a,1x))')     '    MIN_TEMP','    MAX_TEMP','   MIN_WMASS'
         write(unit=*,fmt='(3(es12.5,1x))') rk4min_sfcw_temp ,rk4max_sfcw_temp             &
                                           ,rk4min_sfcw_mass
         write(unit=*,fmt='(a)')           ' '
         write(unit=*,fmt='(78a)')         ('-',k=1,78)
         write(unit=*,fmt='(a)')           ' '
         write(unit=*,fmt='(a)')           ' 4. SOIL (TEXTURE CLASS AT TOP LAYER): '
         write(unit=*,fmt='(4(a,1x))')     '   MIN_WATER','   MAX_WATER','    MIN_TEMP'    &
                                          ,'    MAX_TEMP'
         write(unit=*,fmt='(4(es12.5,1x))') rk4min_soil_water(nzg),rk4max_soil_water(nzg)  &
                                           ,rk4min_soil_temp      ,rk4max_soil_temp
         write(unit=*,fmt='(a)')           ' '
         write(unit=*,fmt='(78a)')         ('=',k=1,78)
         write(unit=*,fmt='(a)')           ' '
      end if

      return
   end subroutine rk4_sanity_check
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This will print the values whenever the step didn't converge due to crazy values  !
   ! no matter how small the steps were (reject_step=.true.).                              !
   !---------------------------------------------------------------------------------------!
   subroutine print_sanity_check(y, csite, ipa)

      use rk4_coms              , only : rk4patchtype  & ! structure
                                       , rk4site       ! ! intent(in)
      use ed_state_vars         , only : sitetype      & ! structure
                                       , patchtype     ! ! structure
      use grid_coms             , only : nzg           ! ! intent(in)
      use soil_coms             , only : soil8         ! ! intent(in), look-up table
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(rk4patchtype) , target     :: y
      type(sitetype)     , target     :: csite
      integer            , intent(in) :: ipa
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
      do k=rk4site%lsl,nzg
         write(unit=*,fmt='(i5,3(1x,es12.4))') &
              k, y%soil_tempk(k), y%soil_fracliq(k), y%soil_water(k)
      end do
      write(unit=*,fmt='(78a)') ('-',k=1,78)

      write(unit=*,fmt='(a)') ' '
      write(unit=*,fmt='(78a)') ('-',k=1,78)
      write(unit=*,fmt='(a5,3(1x,a12))') 'LEVEL','  OLD_SOIL_T','OLD_SOIL_FLQ'             &
                                               &,'OLD_SOIL_H2O'
      do k=rk4site%lsl,nzg
         write(unit=*,fmt='(i5,3(1x,es12.4))')                                             &
              k, csite%soil_tempk(k,ipa), csite%soil_fracliq(k,ipa)                        &
               , csite%soil_water(k,ipa)
      end do
      write(unit=*,fmt='(78a)') ('-',k=1,78)

      write(unit=*,fmt='(a)') ' '
      write(unit=*,fmt='(78a)') ('-',k=1,78)
      write (unit=*,fmt='(a,1x,es12.4)') ' CAN_TEMP=     ',y%can_temp
      write (unit=*,fmt='(a,1x,es12.4)') ' OLD_CAN_TEMP= ',csite%can_temp(ipa)
      write (unit=*,fmt='(a,1x,es12.4)') ' CAN_VAPOR=    ',y%can_shv
      write (unit=*,fmt='(a,1x,es12.4)') ' OLD_CAN_VAP=  ',csite%can_shv(ipa)
      write (unit=*,fmt='(a,1x,i12)')    ' #LEV_SFCH2O=  ',y%nlev_sfcwater
      write (unit=*,fmt='(a,1x,i12)')    ' OLD_#_SFCH2O= ',csite%nlev_sfcwater(ipa)
      if(y%nlev_sfcwater == 1) then
         write(unit=*,fmt='(a,1x,es12.4)') ,'SFCWATER_TEMPK=',y%sfcwater_tempk(1)
      end if
      write(unit=*,fmt='(78a)') ('-',k=1,78)

      write(unit=*,fmt='(a)') ' '
      write(unit=*,fmt='(78a)') ('-',k=1,78)
      cpatch => csite%patch(ipa)
      write (unit=*,fmt='(2(a5,1x),8(a12,1x))')                                            &
         '  COH','  PFT','         LAI','         WAI','         WPA','         TAI'       &
                        ,' LEAF_ENERGY',' OLD_LEAF_EN','   LEAF_TEMP','OLD_LEAF_TMP'
      do ico = 1,cpatch%ncohorts
         if(y%leaf_resolvable(ico)) then
            write(unit=*,fmt='(2(i5,1x),8(es12.4,1x))')                                    &
               ico,cpatch%pft(ico),y%lai(ico),y%wai(ico),y%wpa(ico),y%tai(ico)             &
                  ,y%leaf_energy(ico),cpatch%leaf_energy(ico),y%leaf_temp(ico)             &
                  ,cpatch%leaf_temp(ico)
         end if
      end do
      write(unit=*,fmt='(78a)') ('-',k=1,78)

      write(unit=*,fmt='(a)') ' '
      write(unit=*,fmt='(78a)') ('-',k=1,78)
      write (unit=*,fmt='(2(a5,1x),8(a12,1x))') &
         '  COH','  PFT','         LAI','         WAI','         WPA','         TAI'       &
                        ,'  LEAF_WATER','OLD_LEAF_H2O','   LEAF_HCAP','   LEAF_FLIQ'
      do ico = 1,cpatch%ncohorts
         if(y%leaf_resolvable(ico)) then
            write(unit=*,fmt='(2(i5,1x),8(es12.4,1x))')                                    &
               ico,cpatch%pft(ico),y%lai(ico),y%wai(ico),y%wpa(ico),y%tai(ico)             &
                  ,y%leaf_water(ico),cpatch%leaf_water(ico),cpatch%leaf_hcap(ico)          &
                  ,y%leaf_hcap(ico)
         end if
      end do
      write(unit=*,fmt='(78a)') ('-',k=1,78)
      write(unit=*,fmt='(a)') ' '

      write(unit=*,fmt='(a)') ' '
      write(unit=*,fmt='(78a)') ('-',k=1,78)
      cpatch => csite%patch(ipa)
      write (unit=*,fmt='(2(a5,1x),8(a12,1x))')                                            &
         '  COH','  PFT','         LAI','         WAI','         WPA','         TAI'       &
                        ,' WOOD_ENERGY',' OLD_WOOD_EN','   WOOD_TEMP','OLD_WOOD_TMP'
      do ico = 1,cpatch%ncohorts
         if(y%wood_resolvable(ico)) then
            write(unit=*,fmt='(2(i5,1x),8(es12.4,1x))')                                    &
               ico,cpatch%pft(ico),y%lai(ico),y%wai(ico),y%wpa(ico),y%tai(ico)             &
                  ,y%wood_energy(ico),cpatch%wood_energy(ico),y%wood_temp(ico)             &
                  ,cpatch%wood_temp(ico)
         end if
      end do
      write(unit=*,fmt='(78a)') ('-',k=1,78)

      write(unit=*,fmt='(a)') ' '
      write(unit=*,fmt='(78a)') ('-',k=1,78)
      write (unit=*,fmt='(2(a5,1x),8(a12,1x))') &
         '  COH','  PFT','         LAI','         WAI','         WPA','         TAI'       &
                        ,'  WOOD_WATER','OLD_WOOD_H2O','   WOOD_HCAP','   WOOD_FLIQ'
      do ico = 1,cpatch%ncohorts
         if(y%wood_resolvable(ico)) then
            write(unit=*,fmt='(2(i5,1x),8(es12.4,1x))')                                    &
               ico,cpatch%pft(ico),y%lai(ico),y%wai(ico),y%wpa(ico),y%tai(ico)             &
                  ,y%wood_water(ico),cpatch%wood_water(ico),cpatch%wood_hcap(ico)          &
                  ,y%wood_hcap(ico)
         end if
      end do
      write(unit=*,fmt='(78a)') ('-',k=1,78)
      write(unit=*,fmt='(a)') ' '


      write(unit=*,fmt='(78a)') ('=',k=1,78)
      write(unit=*,fmt='(78a)') ('=',k=1,78)
      write(unit=*,fmt='(a)') ' '

      return
   end subroutine print_sanity_check
   !=======================================================================================!
   !=======================================================================================!
end module rk4_stepper
!==========================================================================================!
!==========================================================================================!
