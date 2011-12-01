!==========================================================================================!
!==========================================================================================!
!     This subroutine is the main step driver for the simple lake model.  This subroutine  !
! will be shared by the three integration schemes, and just the actual step will be called !
! differently.                                                                             !
!------------------------------------------------------------------------------------------!
subroutine integrate_lake(dtfull,htryio)
   use lake_coms    , only : lakemet            & ! intent(in)
                           , tlbeg              & ! intent(inout)
                           , tlend              & ! intent(inout)
                           , dtlake             & ! intent(inout)
                           , dtlakei            & ! intent(inout)
                           , lake_buff          & ! intent(inout)
                           , clone_lakesite     & ! subroutine
                           , integ_lakesite     & ! subroutine
                           , normal_lakesite    & ! subroutine
                           , lake_yscal         & ! subroutine
                           , lake_errmax        ! ! subroutine
   use ed_misc_coms , only : integration_scheme ! ! intent(in)
   use rk4_coms     , only : hmin               & ! intent(in)
                           , maxstp             & ! intent(in)
                           , rk4eps             & ! intent(in)
                           , rk4epsi            & ! intent(in)
                           , safety             & ! intent(in)
                           , pgrow              & ! intent(in)
                           , pshrnk             & ! intent(in)
                           , errcon             ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   real(kind=4), intent(in)    :: dtfull
   real(kind=4), intent(inout) :: htryio
   !----- Local variables. ----------------------------------------------------------------!
   logical                     :: reject_step      ! Should I reject the step?
   logical                     :: reject_result    ! Should I reject the result?
   logical                     :: minstep          ! Minimum time step reached
   logical                     :: stuck            ! Tiny step, it won't advance
   logical                     :: test_reject      ! Reject the test
   integer                     :: i                ! Step counter
   integer                     :: k                ! Banner counter
   real(kind=8)                :: x                ! Elapsed time
   real(kind=8)                :: xnew             ! Elapsed time + h
   real(kind=8)                :: newh             ! New time step suggested
   real(kind=8)                :: oldh             ! Old time step
   real(kind=8)                :: h                ! Current delta-t attempt
   real(kind=8)                :: hnext            ! Next delta-t
   real(kind=8)                :: hdid             ! delta-t that worked (???)
   real(kind=8)                :: errmax           ! Maximum error of this step
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Assign time step dimensions, which will remain the same to the end of this this   !
   ! integration call.                                                                     !
   !---------------------------------------------------------------------------------------!
   tlbeg   = 0.d0
   tlend   = dble(dtfull)
   dtlake  = tlend - tlbeg
   dtlakei = 1.d0/dtlake
   !---------------------------------------------------------------------------------------!


   !----- Initial step size.  We use the previous step as our first guess. ----------------!
   x    = tlbeg
   if (dtlake >= 0.d0) then
      h =  dble(htryio)
   else
      h = -dble(htryio)
   end if
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Copy temporary patches.                                                           !
   !---------------------------------------------------------------------------------------!
   call clone_lakesite(lake_buff%initp,lake_buff%y)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   ! Begin timestep loop                                                                   !
   !---------------------------------------------------------------------------------------!
   timesteploop: do i=1,maxstp

      !----- Get initial derivatives ------------------------------------------------------!
      call lake_derivs(lake_buff%y,lake_buff%dydx)

      !----- Get scalings used to determine stability -------------------------------------!
      call lake_yscal(lake_buff%y,lake_buff%dydx,h,lake_buff%yscal)

      !----- Be sure not to overstep ------------------------------------------------------!
      if((x+h-tlend)*(x+h-tlbeg) > 0.d0) h=tlend-x

      !------------------------------------------------------------------------------------!
      !     Here we will perform the Heun's integration using the time step.  As in Runge- !
      ! Kutta, we also check whether the integration is going well and if needed we shrink !
      ! the intermediate time steps.                                                       !
      !------------------------------------------------------------------------------------!
      reject_step =  .false.
      hstep:   do

         !---------------------------------------------------------------------------------!
         !     Try a step of varying size.                                                 !
         !---------------------------------------------------------------------------------!
         select case (integration_scheme)
         case (0)
            !------------------------------------------------------------------------------!
            !    Euler scheme.  This is very simple so it won't have a routine by itself.  !
            ! Integrate, then update and correct diagnostic variables to avoid overshoot-  !
            ! ing, provided that the overshooting is small.                                !
            !------------------------------------------------------------------------------!
            call clone_lakesite  (lake_buff%y    ,lake_buff%ytemp)
            call integ_lakesite  (lake_buff%ytemp,lake_buff%dydx,h)
            call lake_diagnostics(lake_buff%ytemp)

            !----- Perform a sanity check. ------------------------------------------------!
            call lake_sanity_check(lake_buff%ytemp,reject_step,lake_buff%dydx,h,.false.)
            reject_result = reject_step
            !------------------------------------------------------------------------------!

         case (1)
            call lake_heun(lake_buff%y,lake_buff%dydx,lake_buff%ytemp,lake_buff%yerr       &
                          ,lake_buff%ak2,lake_buff%ak3,x,h,reject_step,reject_result)
         case (2)
            call lake_rk4(lake_buff%y,lake_buff%dydx,lake_buff%ytemp,lake_buff%yerr        &
                         ,lake_buff%ak2,lake_buff%ak3,lake_buff%ak4,lake_buff%ak5          &
                         ,lake_buff%ak6,lake_buff%ak7,x,h,reject_step,reject_result)
         end select
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Here we check the error of this step.  Three outcomes are possible:         !
         ! 1.  The updated values make no sense.  Reject step, assign a large error and    !
         !     try again with a smaller time step;                                         !
         ! 2.  The updated values are reasonable, but the error is large.  Reject step and !
         !     try again with a smaller time step;                                         !
         ! 3.  The updated values are reasonable, and the error is small.  Accept step and !
         !     try again with a larger time step.                                          !
         !---------------------------------------------------------------------------------!
         if (reject_step .or. reject_result) then
            !------------------------------------------------------------------------------!
            !    If step was already rejected, that means the step had finished premature- !
            ! ly, so we assign a standard large error (10.0).                              !
            !------------------------------------------------------------------------------!
            errmax = 1.d1
         elseif (integration_scheme == 0 .or. integration_scheme == 3) then
            !------ Euler scheme, we can't estimate the error, assume it's fine. ----------!
            errmax = 1.d-1
         else
            call lake_errmax(errmax,lake_buff%yscal,lake_buff%yerr)
            !----- Scale the error based on the prescribed tolerance. ---------------------!
            errmax = errmax * rk4epsi
         end if
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !    If that error was large, then calculate a new step size to try.  There are   !
         ! two types of new tries.  If step failed to be minimally reasonable (rejected)   !
         ! we have assigned a standard large error (10.0).  Otherwise a new step is        !
         ! calculated based on the size of that error.  Hopefully, those new steps should  !
         ! be less than the previous h.  If the error was small, i.e. less than rk4eps,    !
         ! then we are done with this step, and we can move forward                        !
         ! time: x = x + h                                                                 !
         !---------------------------------------------------------------------------------!
         if (errmax > 1.d0) then
            !----- Define new step and checking if it can be. -----------------------------!
            oldh    = h
            newh    = safety * h * errmax**pshrnk
            minstep = (newh == h) .or. newh < hmin

            !----- Find next time, and check whether it really added something. -----------!
            h       = max(1.d-1*h, newh)
            xnew    = x + h
            stuck   = xnew == x

            !------------------------------------------------------------------------------!
            !     Here is the moment of truth... If we reached a tiny step and yet the     !
            ! model didn't converge, then we print various values to inform the user and   !
            ! abort the run.  Please, don't hate the messenger...                          !
            !------------------------------------------------------------------------------!
            if (minstep .or. stuck) then

               write (unit=*,fmt='(80a)')         ('=',k=1,80)
               write (unit=*,fmt='(a)')           '   STEPSIZE UNDERFLOW IN INTEGRATE_LAKE'
               write (unit=*,fmt='(80a)')         ('-',k=1,80)
               write (unit=*,fmt='(a,1x,f9.4)')   ' + LONGITUDE:     ',lakemet%lon
               write (unit=*,fmt='(a,1x,f9.4)')   ' + LATITUDE:      ',lakemet%lat
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
               if (reject_step) then
                  write (unit=*,fmt='(a)') '   Likely to be a rejected step problem.'
                  write (unit=*,fmt='(80a)') ('=',k=1,80)
               else
                  write (unit=*,fmt='(a)') '   Likely to be an errmax problem.'
                  write (unit=*,fmt='(80a)') ('=',k=1,80)
               end if

               if (reject_result) then
                  !----- Run the LSM sanity check but this time we force the print. -------!
                  call lake_sanity_check(lake_buff%ytemp,test_reject,lake_buff%dydx,h      &
                                        ,.true.)
               elseif (reject_step) then
                  select case (integration_scheme)
                  case (1)
                     call lake_sanity_check(lake_buff%ak7,test_reject,lake_buff%dydx,h     &
                                           ,.true.)
                  case (2)
                     call lake_sanity_check(lake_buff%ak3,test_reject,lake_buff%dydx,h     &
                                           ,.true.)
                  end select
               else
                  call print_lake_errmax(errmax,lake_buff%yerr,lake_buff%yscal             &
                                        ,lake_buff%y,lake_buff%ytemp)
                  write (unit=*,fmt='(80a)') ('=',k=1,80)
                  write (unit=*,fmt='(a,1x,es12.4)') ' - Rel. errmax:',errmax
                  write (unit=*,fmt='(a,1x,es12.4)') ' - Raw errmax: ',errmax*rk4eps
                  write (unit=*,fmt='(a,1x,es12.4)') ' - Epsilon:',rk4eps
                  write (unit=*,fmt='(80a)') ('=',k=1,80)
               end if
               call print_lakesite(lake_buff%y,lake_buff%initp,h)
            end if

         else
            !------------------------------------------------------------------------------!
            !     Great, it worked, so now we can advance to the next step.  Update the    !
            ! diagnostic variables.                                                        !
            !------------------------------------------------------------------------------!
            call lake_diagnostics(lake_buff%ytemp)
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !     Set up h for the next time.  And here we can relax h for the next step,  !
            ! and try something faster.                                                    !
            !------------------------------------------------------------------------------!
            if (errmax > errcon) then
               hnext = safety * h * errmax**pgrow
            else
               hnext = 5.d0 * h
            endif
            hnext = max(2.d0*hmin,hnext)
            !------------------------------------------------------------------------------!



            !----- Copy the temporary structure to the intermediate state. ----------------!
            call clone_lakesite(lake_buff%ytemp,lake_buff%y)
            !------------------------------------------------------------------------------!


            !----- Update time. -----------------------------------------------------------!
            x        = x + h
            hdid     = h
            h        = hnext
            exit hstep
            !------------------------------------------------------------------------------!
         end if
      end do hstep

      !----- If the integration reached the next step, make some final adjustments --------!
      if((x-tlend)*dtlake >= 0.d0)then

         !------ Copy the temporary patch to the next intermediate step -------------------!
         call clone_lakesite(lake_buff%y,lake_buff%initp)
         !---------------------------------------------------------------------------------!


         !------ Normalise the fluxes so they are back in flux units. ---------------------!
         call normal_lakesite(lake_buff%initp,dtlake)
         !---------------------------------------------------------------------------------!


         !------ Update the substep for next time and leave -------------------------------!
         htryio = sngl(hnext)
         return
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!


      !----- Use hnext as the next substep ------------------------------------------------!
      h = hnext
      !------------------------------------------------------------------------------------!
   end do timesteploop
   !---------------------------------------------------------------------------------------!



   !----- If it reached this point, that is really bad news... ----------------------------!
   write (unit=*,fmt='(a)') ' ==> Too many steps in routine integrate_lake'
   call print_lakesite(lake_buff%y,lake_buff%initp,h)
   !---------------------------------------------------------------------------------------!


   return
end subroutine integrate_lake
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!      This will run the integration of one Heun step for the simple lake model.           !
! Here we use several buffers, so a quick guide of what each one means...                  !
!                                                                                          !
! 1. y      => This is the state y(t)                                                      !
! 2. ytemp  => This is the state y(t+h)                                                    !
! 3. yerr   => This is the error e(t+h)                                                    !
! 4. dydx   => This is the Runge-Kutta term K1                                             !
! 5. ak2    => This is the Runge-Kutta term K2                                             !
! 6. ak3    => This is the next step using Euler: ye(t+h) = y(t) + K1 * h                  !
!                                                                                          !
!------------------------------------------------------------------------------------------!
subroutine lake_heun(y,dydx,ytemp,yerr,ak2,ak3,x,h,reject_step,reject_result)
   use rk4_coms      , only : heun_a2             & ! intent(in)
                            , heun_b21            & ! intent(in)
                            , heun_c1             & ! intent(in)
                            , heun_c2             & ! intent(in)
                            , heun_dc1            & ! intent(in)
                            , heun_dc2            ! ! intent(in)
   use lake_coms     , only : lakesitetype        & ! structure
                            , lakemet             & ! intent(in)
                            , zero_lakesite       & ! subroutine
                            , clone_lakesite      & ! subroutine
                            , integ_lakesite      ! ! subroutine
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(lakesitetype), target      :: y
   type(lakesitetype), target      :: dydx
   type(lakesitetype), target      :: ytemp
   type(lakesitetype), target      :: yerr
   type(lakesitetype), target      :: ak2
   type(lakesitetype), target      :: ak3
   logical           , intent(out) :: reject_step
   logical           , intent(out) :: reject_result
   real(kind=8)      , intent(in)  :: x
   real(kind=8)      , intent(in)  :: h
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Start and assume that nothing went wrong up to this point... If we find any       !
   ! seriously bad step, quit and reduce the time step without even bothering to try       !
   ! further.                                                                              !
   !---------------------------------------------------------------------------------------!
   reject_step   = .false.
   reject_result = .false.
   !---------------------------------------------------------------------------------------!




   !----- ak3 is the temporary array with the Euler step with no correction. --------------!
   call clone_lakesite  (y,ak3)
   call integ_lakesite  (ak3,dydx,heun_b21*h)
   call lake_diagnostics(ak3)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    Check to see if the Euler result makes sense.  Since we will compute the           !
   ! derivative correction using it, the Euler step must be bounded.  If not, reject the   !
   ! step and try a smaller step size.                                                     !
   !---------------------------------------------------------------------------------------!
   call lake_sanity_check(ak3,reject_step,dydx,h,.false.)
   if (reject_step) return
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Compute the second term (correction) of the derivative, using the Euler's         !
   ! predicted state.                                                                      !
   !---------------------------------------------------------------------------------------!
   call lake_derivs(ak3,ak2)
   !---------------------------------------------------------------------------------------!



   !----- We now combine both derivatives and find the potential next step. ---------------!
   call clone_lakesite(y,ytemp)
   call integ_lakesite(ytemp,dydx, heun_c1*h)
   call integ_lakesite(ytemp, ak2, heun_c2*h)

   !---------------------------------------------------------------------------------------!
   !      Update the diagnostic properties.                                                !
   !---------------------------------------------------------------------------------------!
   call lake_diagnostics(ytemp)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    Check to see if this attempt of advancing one time step makes sense.  If not,      !
   ! reject the result and try a smaller step size.                                        !
   !---------------------------------------------------------------------------------------!
   call lake_sanity_check(ytemp,reject_result,ak2,h,.false.)
   if(reject_result)return
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Compute the estimate of the error associated with the step.                       !
   !---------------------------------------------------------------------------------------!
   call zero_lakesite (yerr)
   call integ_lakesite(yerr,dydx,heun_dc1*h)
   call integ_lakesite(yerr, ak2,heun_dc2*h)
   !---------------------------------------------------------------------------------------!

   return
end subroutine lake_heun
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine will update the variables and perform the actual time stepping.       !
!------------------------------------------------------------------------------------------!
subroutine lake_rk4(y,dydx,yout,yerr,ak2,ak3,ak4,ak5,ak6,ak7,x,h,reject_step,reject_result)
   use rk4_coms      , only : rk4_a2              & ! intent(in)
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
                            , rk4_dc6             ! ! intent(in)
   use lake_coms     , only : lakesitetype        & ! structure
                            , lakemet             & ! intent(in)
                            , zero_lakesite       & ! subroutine
                            , integ_lakesite      & ! subroutine
                            , clone_lakesite      ! ! subroutine
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(lakesitetype), target      :: y
   type(lakesitetype), target      :: dydx
   type(lakesitetype), target      :: yout
   type(lakesitetype), target      :: yerr
   type(lakesitetype), target      :: ak2
   type(lakesitetype), target      :: ak3
   type(lakesitetype), target      :: ak4
   type(lakesitetype), target      :: ak5
   type(lakesitetype), target      :: ak6
   type(lakesitetype), target      :: ak7
   logical           , intent(out) :: reject_step
   logical           , intent(out) :: reject_result
   real(kind=8)      , intent(in)  :: x
   real(kind=8)      , intent(in)  :: h
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Start and assume that nothing went wrong up to this point... If we find any       !
   ! seriously bad step, quit and reduce the time step without even bothering to try       !
   ! further.                                                                              !
   !---------------------------------------------------------------------------------------!
   reject_step   = .false.
   reject_result = .false.
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     For each stage we checked the sanity, so we avoid using non-sense values to       !
   ! advance.  Also, we estimate the derivative of pressure after each stage, and in       !
   ! the end we estimate the full step derivative as the weighted average for each of      !
   ! these partial steps taken.                                                            !
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Second stage (the first was the Euler, whose derivative is already computed       !
   ! and saved in dydx.                                                                    !
   !---------------------------------------------------------------------------------------!
   call clone_lakesite   (y, ak7)
   call integ_lakesite   (ak7, dydx, rk4_b21*h)
   call lake_diagnostics (ak7)
   call lake_sanity_check(ak7, reject_step,dydx,h,.false.)
   if (reject_step) return
   !---------------------------------------------------------------------------------------!


   !------ Get the new derivative evaluation. ---------------------------------------------!
   call lake_derivs(ak7, ak2)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    Third stage.                                                                       !
   !---------------------------------------------------------------------------------------!
   call clone_lakesite   (  y,  ak7)
   call integ_lakesite   (ak7, dydx, rk4_b31*h)
   call integ_lakesite   (ak7,  ak2, rk4_b32*h)
   call lake_diagnostics (ak7)
   call lake_sanity_check(ak7,reject_step,dydx,h,.false.)
   if (reject_step) return
   !---------------------------------------------------------------------------------------!


   !------ Get the new derivative evaluation. ---------------------------------------------!
   call lake_derivs(ak7, ak3)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    Fourth stage.                                                                      !
   !---------------------------------------------------------------------------------------!
   call clone_lakesite   (  y,  ak7)
   call integ_lakesite   (ak7, dydx, rk4_b41*h)
   call integ_lakesite   (ak7,  ak2, rk4_b42*h)
   call integ_lakesite   (ak7,  ak3, rk4_b43*h)
   call lake_diagnostics (ak7)
   call lake_sanity_check(ak7, reject_step,dydx,h,.false.)
   if (reject_step) return
   !---------------------------------------------------------------------------------------!


   !------ Get the new derivative evaluation. ---------------------------------------------!
   call lake_derivs(ak7, ak4)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    Fifth stage.                                                                       !
   !---------------------------------------------------------------------------------------!
   call clone_lakesite   (  y,  ak7)
   call integ_lakesite   (ak7, dydx, rk4_b51*h)
   call integ_lakesite   (ak7,  ak2, rk4_b52*h)
   call integ_lakesite   (ak7,  ak3, rk4_b53*h)
   call integ_lakesite   (ak7,  ak4, rk4_b54*h)
   call lake_diagnostics (ak7)
   call lake_sanity_check(ak7,reject_step,dydx,h,.false.)
   if (reject_step) return
   !---------------------------------------------------------------------------------------!


   !------ Get the new derivative evaluation. ---------------------------------------------!
   call lake_derivs(ak7, ak5)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    Sixth stage.                                                                       !
   !---------------------------------------------------------------------------------------!
   call clone_lakesite   (  y,  ak7)
   call integ_lakesite   (ak7, dydx, rk4_b61*h)
   call integ_lakesite   (ak7,  ak2, rk4_b62*h)
   call integ_lakesite   (ak7,  ak3, rk4_b63*h)
   call integ_lakesite   (ak7,  ak4, rk4_b64*h)
   call integ_lakesite   (ak7,  ak5, rk4_b65*h)
   call lake_diagnostics (ak7)
   call lake_sanity_check(ak7, reject_step,dydx,h,.false.)
   if(reject_step)return
   !---------------------------------------------------------------------------------------!


   !------ Get the new derivative evaluation. ---------------------------------------------!
   call lake_derivs(ak7, ak6)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    Aggregate all derivatives to make the new guess.                                   !
   !---------------------------------------------------------------------------------------!
   call clone_lakesite   (   y, yout)
   call integ_lakesite   (yout, dydx, rk4_c1*h)
   call integ_lakesite   (yout,  ak3, rk4_c3*h)
   call integ_lakesite   (yout,  ak4, rk4_c4*h)
   call integ_lakesite   (yout,  ak6, rk4_c6*h)
   call lake_diagnostics (yout)
   call lake_sanity_check(yout, reject_result,dydx,h,.false.)
   !---------------------------------------------------------------------------------------!
   if(reject_result)return
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    Estimate the error for this step.                                                  !
   !---------------------------------------------------------------------------------------!
   call zero_lakesite (yerr)
   call integ_lakesite(yerr, dydx, rk4_dc1*h)
   call integ_lakesite(yerr, ak3,  rk4_dc3*h)
   call integ_lakesite(yerr, ak4,  rk4_dc4*h)
   call integ_lakesite(yerr, ak5,  rk4_dc5*h)
   call integ_lakesite(yerr, ak6,  rk4_dc6*h)
   !---------------------------------------------------------------------------------------!

   return
end subroutine lake_rk4
!==========================================================================================!
!==========================================================================================!
