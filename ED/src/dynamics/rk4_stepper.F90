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
   subroutine rkqs(integration_buff, x,htry,hdid,hnext,csite,ipa,isi,ipy,ifm)

      use rk4_coms      , only : rk4patchtype        & ! structure
                               , integration_vars & ! structure
                               , rk4met              & ! intent(in)
                               , hmin                & ! intent(in)
                               , rk4eps              & ! intent(in)
                               , rk4epsi             & ! intent(in)
                               , safety              & ! intent(in)
                               , pgrow               & ! intent(in)
                               , pshrnk              & ! intent(in)
                               , errcon              ! ! intent(in)
      use ed_state_vars , only : sitetype            & ! structure
                               , patchtype           & ! structure
                               , edtype              & ! structure
                               , edgrid_g            ! ! structure

      implicit none
      !----- Arguments --------------------------------------------------------------------!
      integer                  , intent(in)    :: ipa,isi,ipy,ifm
      type(sitetype)           , target        :: csite
      type(integration_vars), target        :: integration_buff
      real(kind=8)             , intent(in)    :: htry
      real(kind=8)             , intent(inout) :: x
      real(kind=8)             , intent(out)   :: hdid,hnext
      !----- Local variables --------------------------------------------------------------!
      type(edtype)             , pointer       :: cgrid
      real(kind=8)                             :: h,errmax,xnew,newh,oldh
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
         call rkck(integration_buff%y,integration_buff%dydx,integration_buff%ytemp         &
                     ,integration_buff%yerr,integration_buff%ak2,integration_buff%ak3      &
                     ,integration_buff%ak4,integration_buff%ak5,integration_buff%ak6       &
                     ,integration_buff%ak7,x,h,csite,ipa,isi,ipy,reject_step)


         !---------------------------------------------------------------------------------!
         ! 2. Check to see how accurate the step was.  Errors were calculated by integrat- !
         !    ing the derivative of that last step.                                        !
         !---------------------------------------------------------------------------------!
         if (reject_step) then
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
         !    should be less then the previous h.  If the error was small, i.e. less then  !
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
               write (unit=*,fmt='(a)')           ' + PATCH INFO:  '
               write (unit=*,fmt='(a,1x,es12.4)') '   - AGE:       ',csite%age(ipa)
               write (unit=*,fmt='(a,1x,i6)')     '   - DIST_TYPE: ',csite%dist_type(ipa)
               write (unit=*,fmt='(a,1x,l1)')     ' + REJECT_STEP: ',reject_step
               write (unit=*,fmt='(a,1x,l1)')     ' + MINSTEP:     ',minstep
               write (unit=*,fmt='(a,1x,l1)')     ' + STUCK:       ',stuck
               write (unit=*,fmt='(a,1x,es12.4)') ' + ERRMAX:      ',errmax
               write (unit=*,fmt='(a,1x,es12.4)') ' + X:           ',x
               write (unit=*,fmt='(a,1x,es12.4)') ' + H:           ',h
               write (unit=*,fmt='(a,1x,es12.4)') ' + OLDH:        ',oldh
               write (unit=*,fmt='(a,1x,es12.4)') ' + NEWH:        ',newh
               write (unit=*,fmt='(a,1x,es12.4)') ' + SAFETY:      ',safety
               write (unit=*,fmt='(a,1x,es12.4)') ' + ERRMAX:      ',errmax
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
                  call rk4_sanity_check(integration_buff%ytemp,test_reject,csite,ipa       &
                                       ,integration_buff%dydx,h,.true.)
                  if (.not. test_reject) then
                     call rk4_sanity_check(integration_buff%ak7,test_reject,csite,ipa      &
                                          ,integration_buff%dydx,h,.true.)
                  end if
                  call print_sanity_check(integration_buff%y,csite,ipa)
               else
                  call print_errmax(errmax,integration_buff%yerr,integration_buff%yscal    &
                                      ,csite%patch(ipa),integration_buff%y                 &
                                      ,integration_buff%ytemp)
                  write (unit=*,fmt='(80a)') ('=',k=1,80)
                  write (unit=*,fmt='(a,1x,es12.4)') ' - Rel. errmax:',errmax*rk4epsi
                  write (unit=*,fmt='(a,1x,es12.4)') ' - Raw errmax: ',errmax
                  write (unit=*,fmt='(a,1x,es12.4)') ' - Epsilon:',rk4eps
                  write (unit=*,fmt='(80a)') ('=',k=1,80)
               end if
               call print_rk4patch(integration_buff%y, csite,ipa)
            endif
         
         else
            !------------------------------------------------------------------------------!
            !   Great, it worked, so now we can advance to the next step.  We just need to !
            ! do some minor adjustments before...                                          !
            !------------------------------------------------------------------------------!
            !----- 1. Final update of leaf properties to avoid negative water. ------------!
            call adjust_veg_properties(integration_buff%ytemp,h,csite,ipa)
            !----- 2. Make snow layers stable and positively defined. ---------------------!
            call redistribute_snow(integration_buff%ytemp, csite,ipa)
            !----- 3. Update the diagnostic variables. ------------------------------------!
            call update_diagnostic_vars(integration_buff%ytemp, csite,ipa)
            
            !------------------------------------------------------------------------------!
            ! 4. Set up h for the next time.  And here we can relax h for the next step,   !
            !    and try something faster.                                                 !
            !------------------------------------------------------------------------------!
            if (errmax > errcon) then
               hnext = safety * h * errmax**pgrow
            else
               hnext = 5.d0 * h
            endif
            hnext = max(2.d0*hmin,hnext)

            !----- 5. Updating time. ------------------------------------------------------!
            x    = x + h
            hdid = h
            
            !----- 6. Copying the temporary structure to the intermediate state. ----------!
            call copy_rk4_patch(integration_buff%ytemp,integration_buff%y                  &
                                  ,csite%patch(ipa))
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
   subroutine rkck(y,dydx,yout,yerr,ak2,ak3,ak4,ak5,ak6,ak7,x,h,csite,ipa,isi,ipy          &
                     ,reject_step)

      use rk4_coms      , only : rk4patchtype        & ! structure
                               , integration_vars    & ! structure
                               , rk4met              & ! intent(in)
                               , print_diags         & ! intent(in)
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
      use ed_state_vars , only : sitetype            & ! structure
                               , patchtype           ! ! structure
      implicit none

      !----- Arguments --------------------------------------------------------------------!
      integer           , intent(in)  :: ipa,isi,ipy
      real(kind=8)      , intent(in)  :: x,h
      type(rk4patchtype), target      :: y,dydx,yout,yerr
      type(rk4patchtype), target      :: ak2,ak3,ak4,ak5,ak6,ak7
      type(sitetype)    , target      :: csite
      logical           , intent(out) :: reject_step
      !----- Local variables --------------------------------------------------------------!
      type(patchtype)   , pointer     :: cpatch
      real(kind=8)                    :: combh
      !------------------------------------------------------------------------------------!


      !----- Interfaces, in case the model is compiled without forcing them. --------------!
#if USE_INTERF
      interface
         subroutine leaf_derivs(initp,dydx,csite,ipa,isi,ipy)
            use rk4_coms      ,only : rk4patchtype
            use ed_state_vars ,only : sitetype
            implicit none
            integer             , intent(in) :: ipa,isi,ipy
            type (rk4patchtype) , target     :: initp
            type (rk4patchtype) , target     :: dydx
            type (sitetype)     , target     :: csite
         end subroutine leaf_derivs
      end interface
#endif

      !------------------------------------------------------------------------------------!
      !     Start and assume that nothing went wrong up to this point... If we find any    !
      ! seriously bad step, quit and reduce the time step without even bothering to try    !
      ! further.                                                                           !
      !------------------------------------------------------------------------------------!
      reject_step = .false.
      cpatch => csite%patch(ipa)

      call copy_rk4_patch(y, ak7, cpatch)
      call inc_rk4_patch(ak7, dydx, b21*h, cpatch)
      combh = b21*h
      call adjust_veg_properties(ak7,combh,csite,ipa)
      call redistribute_snow(ak7, csite,ipa)
      call update_diagnostic_vars(ak7, csite,ipa)
      call rk4_sanity_check(ak7, reject_step, csite, ipa,dydx,h,print_diags)
      if (reject_step) return



      call leaf_derivs(ak7, ak2, csite, ipa,isi,ipy)
      call copy_rk4_patch(y, ak7, cpatch)
      call inc_rk4_patch(ak7, dydx, b31*h, cpatch)
      call inc_rk4_patch(ak7, ak2, b32*h, cpatch)
      combh = (b31+b32)*h
      call adjust_veg_properties(ak7,combh,csite,ipa)
      call redistribute_snow(ak7, csite,ipa)
      call update_diagnostic_vars(ak7, csite,ipa)
      call rk4_sanity_check(ak7,reject_step,csite,ipa,dydx,h,print_diags)
      if (reject_step) return



      call leaf_derivs(ak7, ak3, csite,ipa,isi,ipy)
      call copy_rk4_patch(y, ak7, cpatch)
      call inc_rk4_patch(ak7, dydx, b41*h, cpatch)
      call inc_rk4_patch(ak7,  ak2, b42*h, cpatch)
      call inc_rk4_patch(ak7,  ak3, b43*h, cpatch)
      combh = (b41+b42+b43)*h
      call adjust_veg_properties(ak7,combh,csite,ipa)
      call redistribute_snow(ak7, csite,ipa)
      call update_diagnostic_vars(ak7, csite,ipa)
      call rk4_sanity_check(ak7, reject_step, csite,ipa,dydx,h,print_diags)
      if (reject_step) return



      call leaf_derivs(ak7, ak4, csite, ipa,isi,ipy)
      call copy_rk4_patch(y, ak7, cpatch)
      call inc_rk4_patch(ak7, dydx, b51*h, cpatch)
      call inc_rk4_patch(ak7,  ak2, b52*h, cpatch)
      call inc_rk4_patch(ak7,  ak3, b53*h, cpatch)
      call inc_rk4_patch(ak7,  ak4, b54*h, cpatch)
      combh = (b51+b52+b53+b54)*h
      call adjust_veg_properties(ak7,combh,csite,ipa)
      call redistribute_snow(ak7, csite,ipa)
      call update_diagnostic_vars(ak7, csite,ipa)
      call rk4_sanity_check(ak7,reject_step,csite,ipa,dydx,h,print_diags)
      if (reject_step) return



      call leaf_derivs(ak7, ak5, csite, ipa,isi,ipy)
      call copy_rk4_patch(y, ak7, cpatch)
      call inc_rk4_patch(ak7, dydx, b61*h, cpatch)
      call inc_rk4_patch(ak7,  ak2, b62*h, cpatch)
      call inc_rk4_patch(ak7,  ak3, b63*h, cpatch)
      call inc_rk4_patch(ak7,  ak4, b64*h, cpatch)
      call inc_rk4_patch(ak7,  ak5, b65*h, cpatch)
      combh = (b61+b62+b63+b64+b65)*h
      call adjust_veg_properties(ak7,combh,csite,ipa)
      call redistribute_snow(ak7, csite,ipa)
      call update_diagnostic_vars(ak7, csite,ipa)
      call rk4_sanity_check(ak7, reject_step, csite,ipa,dydx,h,print_diags)
      if(reject_step)return

      call leaf_derivs(ak7, ak6, csite,ipa,isi,ipy)
      call copy_rk4_patch(y, yout, cpatch)
      call inc_rk4_patch(yout, dydx, c1*h, cpatch)
      call inc_rk4_patch(yout,  ak3, c3*h, cpatch)
      call inc_rk4_patch(yout,  ak4, c4*h, cpatch)
      call inc_rk4_patch(yout,  ak6, c6*h, cpatch)
      combh = (c1+c3+c4+c6)*h
      call adjust_veg_properties(yout,combh,csite,ipa)
      call redistribute_snow(yout, csite,ipa)
      call update_diagnostic_vars(yout, csite,ipa)
      call rk4_sanity_check(yout, reject_step, csite,ipa,dydx,h,print_diags)
      if(reject_step)return

      call copy_rk4_patch(dydx, yerr, cpatch)
      call inc_rk4_patch(yerr, dydx, dc1*h-1.d0, cpatch)
      call inc_rk4_patch(yerr, ak3,  dc3*h     , cpatch)
      call inc_rk4_patch(yerr, ak4,  dc4*h     , cpatch)
      call inc_rk4_patch(yerr, ak5,  dc5*h     , cpatch)
      call inc_rk4_patch(yerr, ak6,  dc6*h     , cpatch)

      return
   end subroutine rkck
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine will check for potentially serious problems.  Note that the upper  !
   ! and lower bound are defined in rk4_coms.f90, so if you need to change any limit for   !
   ! some reason, you can adjust there. (Only exception is veg_temp_min, which is in       !
   ! canopy_radiation_coms).                                                               !
   !---------------------------------------------------------------------------------------!
   subroutine rk4_sanity_check(y,reject_step, csite,ipa,dydx,h,print_problems)
      use rk4_coms              , only : rk4patchtype         & ! structure
                                       , integration_vars     & ! structure
                                       , rk4met               & ! intent(in)
                                       , rk4eps               & ! intent(in)
                                       , toocold              & ! intent(in)
                                       , rk4min_can_temp      & ! intent(in)
                                       , rk4max_can_temp      & ! intent(in)
                                       , rk4max_can_rhv       & ! intent(in)
                                       , rk4max_can_shv       & ! intent(in)
                                       , rk4min_can_shv       & ! intent(in)
                                       , rk4min_can_co2       & ! intent(in)
                                       , rk4max_can_co2       & ! intent(in)
                                       , rk4max_veg_temp      & ! intent(in)
                                       , rk4min_veg_temp      & ! intent(in)
                                       , rk4min_veg_lwater    & ! intent(in)
                                       , rk4min_sfcw_temp     & ! intent(in)
                                       , rk4max_sfcw_temp     & ! intent(in)
                                       , rk4max_soil_temp     & ! intent(in)
                                       , rk4min_soil_temp     & ! intent(in)
                                       , rk4min_sfcw_mass     & ! intent(in)
                                       , rk4min_virt_water    & ! intent(in)
                                       , rk4min_sfcwater_mass ! ! intent(in)
      use ed_state_vars         , only : sitetype             & ! structure
                                       , patchtype            ! ! structure
      use grid_coms             , only : nzg                  ! ! intent(in)
      use soil_coms             , only : soil8                ! ! intent(in), lookup table
      use consts_coms           , only : t3ple8               ! ! intent(in)
      use ed_misc_coms          , only : integ_err            & ! intent(inout)
                                       , record_err           ! ! intent(inout)
      use therm_lib             , only : rehuil               & ! function
                                       , qtk8                 ! ! subroutine

      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(rk4patchtype) , target      :: y,dydx
      type(sitetype)     , target      :: csite
      logical            , intent(in)  :: print_problems
      logical            , intent(out) :: reject_step
      !----- Local variables --------------------------------------------------------------!
      type(patchtype)    , pointer     :: cpatch
      integer                          :: k
      integer                          :: ksn
      real(kind=8)                     :: h
      real(kind=8)                     :: can_rhv
      real(kind=8)                     :: total_sfcw_energy
      real(kind=8)                     :: total_sfcw_mass
      real(kind=8)                     :: mean_sfcw_tempk
      real(kind=8)                     :: mean_sfcw_fracliq
      integer                          :: ipa,ico
      logical                          :: cflag1,cflag2,introuble
      !------------------------------------------------------------------------------------!


      !----- Assigning initial values for flags. ------------------------------------------!
      cflag1 = .false.
      cflag2 = .false.

      !----- First check, if top soil temperature is NaN, end of story... -----------------!
      if (y%soil_tempk(nzg) /= y%soil_tempk(nzg)) then
         write (unit=*,fmt='(a)') 'Top soil temperature is NaN!!!'
         call print_rk4patch(y,csite,ipa)
      end if

      !----- Being optimistic and assuming things are fine --------------------------------!
      reject_step = .false.
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !  MLO.  The sanity check on liquid fraction against temperature and values of       !
      !        liquid water fraction were removed, since the only way they would fail is   !
      !        if the thermodynamic functions were incorrect. The soil properties may be   !
      !        absurd, but the liquid fraction should be consistent with that.  Same       !
      !        argument for the temporary surface water.                                   !
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !   Checking whether the canopy temperature is too hot or too cold.                  !
      !------------------------------------------------------------------------------------! 
      if (y%can_temp > rk4max_can_temp .or. y%can_temp < rk4min_can_temp) then
         reject_step = .true.
         if(record_err) integ_err(1,2) = integ_err(1,2) + 1_8
         if (print_problems) then
            write(unit=*,fmt='(a)')           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            write(unit=*,fmt='(a)')           ' + Canopy air temperature is off-track...'
            write(unit=*,fmt='(a)')           '-----------------------------------------'
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_SHV:       ',y%can_shv
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_RHV:       ',can_rhv
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_TEMP:      ',y%can_temp
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_CO2:       ',y%can_co2
            write(unit=*,fmt='(a,1x,es12.4)') ' PRESSURE:      ',rk4met%prss
            write(unit=*,fmt='(a,1x,es12.4)') ' D(CAN_TEMP)/Dt:',dydx%can_temp
            write(unit=*,fmt='(a,1x,es12.4)') ' D(CAN_SHV)/Dt: ',dydx%can_shv
            write(unit=*,fmt='(a,1x,es12.4)') ' D(CAN_CO2)/Dt: ',dydx%can_co2
            write(unit=*,fmt='(a,1x,es12.4)') ' H:            ',h
            write(unit=*,fmt='(a)')           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         elseif (.not. record_err) then
            return
         end if
      end if
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     The check of canopy humidity is done only when temperature makes sense, to     !
      ! avoid floating point exceptions when temperature is too cold or too hot.           !
      !------------------------------------------------------------------------------------!
      if (.not. reject_step) then
         !----- Finding canopy relative humidity ------------------------------------------!
         can_rhv = dble(rehuil(sngl(rk4met%prss),sngl(max(y%can_temp,toocold))             &
                              ,sngl(y%can_shv)))

         !----- Checking whether the canopy air is too dry or too humid. ------------------!
         if ((can_rhv > rk4max_can_rhv .and. y%can_shv > rk4max_can_shv) .or.              &
             y%can_shv < rk4min_can_shv )then
            reject_step = .true.
            if(record_err) integ_err(2,2) = integ_err(2,2) + 1_8
            if (print_problems) then
               write(unit=*,fmt='(a)')           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
               write(unit=*,fmt='(a)')           ' + Canopy air spec. hum. is off-track...'
               write(unit=*,fmt='(a)')           '----------------------------------------'
               write(unit=*,fmt='(a,1x,es12.4)') ' CAN_SHV:       ',y%can_shv
               write(unit=*,fmt='(a,1x,es12.4)') ' CAN_RHV:       ',can_rhv
               write(unit=*,fmt='(a,1x,es12.4)') ' CAN_TEMP:      ',y%can_temp
               write(unit=*,fmt='(a,1x,es12.4)') ' CAN_CO2:       ',y%can_co2
               write(unit=*,fmt='(a,1x,es12.4)') ' PRESSURE:      ',rk4met%prss
               write(unit=*,fmt='(a,1x,es12.4)') ' D(CAN_TEMP)/Dt:',dydx%can_temp
               write(unit=*,fmt='(a,1x,es12.4)') ' D(CAN_SHV)/Dt: ',dydx%can_shv
               write(unit=*,fmt='(a,1x,es12.4)') ' D(CAN_CO2)/Dt: ',dydx%can_co2
               write(unit=*,fmt='(a,1x,es12.4)') ' H:            ',h
               write(unit=*,fmt='(a)')           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            elseif (.not. record_err) then
               return
            end if
         end if
      end if
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !   Checking whether the canopy CO2 is reasonable.                                   !
      !------------------------------------------------------------------------------------! 
      !if (y%can_co2 > rk4max_can_co2 .or. y%can_co2 < rk4min_can_co2) then
      !   reject_step = .true.
      !   if(record_err) integ_err(3,2) = integ_err(3,2) + 1_8
      !   if (print_problems) then
      !      write(unit=*,fmt='(a)')           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      !      write(unit=*,fmt='(a)')           ' + Canopy air CO2  is off-track...       '
      !      write(unit=*,fmt='(a)')           '-----------------------------------------'
      !      write(unit=*,fmt='(a,1x,es12.4)') ' CAN_SHV:       ',y%can_shv
      !      write(unit=*,fmt='(a,1x,es12.4)') ' CAN_RHV:       ',can_rhv
      !      write(unit=*,fmt='(a,1x,es12.4)') ' CAN_TEMP:      ',y%can_temp
      !      write(unit=*,fmt='(a,1x,es12.4)') ' CAN_CO2:       ',y%can_co2
      !      write(unit=*,fmt='(a,1x,es12.4)') ' PRESSURE:      ',rk4met%prss
      !      write(unit=*,fmt='(a,1x,es12.4)') ' D(CAN_TEMP)/Dt:',dydx%can_temp
      !      write(unit=*,fmt='(a,1x,es12.4)') ' D(CAN_SHV)/Dt: ',dydx%can_shv
      !      write(unit=*,fmt='(a,1x,es12.4)') ' D(CAN_CO2)/Dt: ',dydx%can_co2
      !      write(unit=*,fmt='(a,1x,es12.4)') ' H:            ',h
      !      write(unit=*,fmt='(a)')           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      !   elseif (.not. record_err) then
      !      return
      !   end if
      !end if
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !    Checking whether the temporary snow/water layer(s) has(ve) reasonable values.   !
      !------------------------------------------------------------------------------------!
      ksn = y%nlev_sfcwater
      if (ksn >= 1) then
         total_sfcw_energy = sum(y%sfcwater_energy(1:ksn))
         total_sfcw_mass   = sum(y%sfcwater_mass(1:ksn))
         if (abs(total_sfcw_mass) > rk4min_sfcwater_mass) then
            call qtk8(total_sfcw_energy/total_sfcw_mass,mean_sfcw_tempk,mean_sfcw_fracliq)
         else
            mean_sfcw_tempk   = t3ple8
            mean_sfcw_fracliq = 5.d-1
         end if

         !----- Temperature ---------------------------------------------------------------!
         if (mean_sfcw_tempk < rk4min_sfcw_temp .or. mean_sfcw_tempk > rk4max_sfcw_temp)   &
         then
            reject_step = .true.
            if(record_err) integ_err(44,2) = integ_err(44,2) + 1_8
            if (print_problems) then
               write(unit=*,fmt='(a)')              '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
               write(unit=*,fmt='(a)')              ' + Snow/pond temperature is off...'
               write(unit=*,fmt='(a)')              '-------------------------------------'
               write(unit=*,fmt='(a,1x,es12.4)')    ' TOTAL_MASS:  ',total_sfcw_mass
               write(unit=*,fmt='(a,1x,es12.4)')    ' TOTAL_ENERGY:',total_sfcw_energy
               write(unit=*,fmt='(a,1x,es12.4)')    ' MEAN_TEMPK:  ',mean_sfcw_tempk
               write(unit=*,fmt='(a,1x,es12.4)')    ' MEAN_FLIQ:   ',mean_sfcw_fracliq
               write(unit=*,fmt='(a)')              '-------------------------------------'
               do k=1,ksn
                  write(unit=*,fmt='(a,1x,i6)')     ' # of layer   ',k
                  write(unit=*,fmt='(a,1x,es12.4)') ' SFCW_TEMP:   ',y%sfcwater_tempk(k)
                  write(unit=*,fmt='(a,1x,es12.4)') ' SFCW_ENERGY: ',y%sfcwater_energy(k)
                  write(unit=*,fmt='(a,1x,es12.4)') ' SFCW_MASS:   ',y%sfcwater_mass(k)
                  write(unit=*,fmt='(a,1x,es12.4)') ' SFCW_DEPTH:  ',y%sfcwater_depth(k)
                  write(unit=*,fmt='(a,1x,es12.4)') ' D(SFCW_E)/Dt:',dydx%sfcwater_energy(k)
                  write(unit=*,fmt='(a,1x,es12.4)') ' D(SFCW_M)/Dt:',dydx%sfcwater_mass(k)
                  write(unit=*,fmt='(a)')           '-------------------------------------'
               end do
               write(unit=*,fmt='(a,1x,es12.4)')    ' H:           ',h
               write(unit=*,fmt='(a)')              '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            elseif (.not. record_err) then
               return
            end if
         end if

         !----- Mass ----------------------------------------------------------------------!
         if (total_sfcw_mass < rk4min_sfcw_mass) then
            reject_step = .true.
            if(record_err) integ_err(32+ksn,2) = integ_err(32+ksn,2) + 1_8
            if (print_problems) then
               write(unit=*,fmt='(a)')              '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
               write(unit=*,fmt='(a)')              ' + Snow/pond layer has weird mass...'
               write(unit=*,fmt='(a)')              '-------------------------------------'
               write(unit=*,fmt='(a,1x,es12.4)')    ' TOTAL_MASS:  ',total_sfcw_mass
               write(unit=*,fmt='(a,1x,es12.4)')    ' TOTAL_ENERGY:',total_sfcw_energy
               write(unit=*,fmt='(a,1x,es12.4)')    ' MEAN_TEMPK:  ',mean_sfcw_tempk
               write(unit=*,fmt='(a,1x,es12.4)')    ' MEAN_FLIQ:   ',mean_sfcw_fracliq
               write(unit=*,fmt='(a)')              '-------------------------------------'
               do k=1,ksn
                  write(unit=*,fmt='(a,1x,i6)')     ' # of layer   ',k
                  write(unit=*,fmt='(a,1x,es12.4)') ' SFCW_TEMP:   ',y%sfcwater_tempk(k)
                  write(unit=*,fmt='(a,1x,es12.4)') ' SFCW_ENERGY: ',y%sfcwater_energy(k)
                  write(unit=*,fmt='(a,1x,es12.4)') ' SFCW_MASS:   ',y%sfcwater_mass(k)
                  write(unit=*,fmt='(a,1x,es12.4)') ' SFCW_DEPTH:  ',y%sfcwater_depth(k)
                  write(unit=*,fmt='(a,1x,es12.4)') ' D(SFCW_E)/Dt:',dydx%sfcwater_energy(k)
                  write(unit=*,fmt='(a,1x,es12.4)') ' D(SFCW_M)/Dt:',dydx%sfcwater_mass(k)
                  write(unit=*,fmt='(a)')           '-------------------------------------'
               end do
               write(unit=*,fmt='(a,1x,es12.4)')    ' H:           ',h
               write(unit=*,fmt='(a)')              '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
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
      do k=rk4met%lsl,nzg
         !----- Soil moisture -------------------------------------------------------------!
         if (y%soil_water(k)< (1.d0-rk4eps)*soil8(csite%ntext_soil(k,ipa))%soilcp .or.     &
             y%soil_water(k)> (1.d0+rk4eps)*soil8(csite%ntext_soil(k,ipa))%slmsts     )    &
         then
            reject_step = .true.
            if(record_err) integ_err(3+k,2) = integ_err(3+k,2) + 1_8
            if (print_problems) then
               write(unit=*,fmt='(a)')           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
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
               write(unit=*,fmt='(a)')           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            elseif (.not. record_err) then
               return
            end if
         end if

         !----- Soil temperature ----------------------------------------------------------!
         if (y%soil_tempk(k) > rk4max_soil_temp .or. y%soil_tempk(k) < rk4min_soil_temp )  &
         then
            reject_step = .true.
            if(record_err) cflag1 = .true.
            if (print_problems) then
               write(unit=*,fmt='(a)')           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
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
            write(unit=*,fmt='(a,1x,es12.4)') ' VIRT_WATER:      ',y%virtual_water
            write(unit=*,fmt='(a,1x,es12.4)') ' VIRT_HEAT :      ',y%virtual_heat
            write(unit=*,fmt='(a,1x,es12.4)') ' D(VIRT_WATER)/Dt:',dydx%virtual_water
            write(unit=*,fmt='(a,1x,es12.4)') ' D(VIRT_HEAT)/Dt :',dydx%virtual_heat
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
         if (y%solvable(ico)) then
            if (y%veg_temp(ico) > rk4max_veg_temp .or. y%veg_temp(ico) < rk4min_veg_temp)  &
            then
               reject_step = .true.
               if(record_err) cflag1 = .true.
               if (print_problems) then
                  write(unit=*,fmt='(a)')           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                  write(unit=*,fmt='(a)')           ' + Leaf temperature is off-track...'
                  write(unit=*,fmt='(a)')           '--------------------------------------'
                  write(unit=*,fmt='(a,1x,i6)')     ' PFT:          ',cpatch%pft(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' LAI:          ',y%lai(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' WPA:          ',y%wpa(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' TAI:          ',y%tai(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' HCAPVEG:      ',y%hcapveg(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' VEG_TEMP:     ',y%veg_temp(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' VEG_FRACLIQ:  ',y%veg_fliq(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' VEG_ENERGY:   ',y%veg_energy(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' VEG_WATER:    ',y%veg_water(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' D(VEG_EN)/Dt: ',dydx%veg_energy(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' D(VEG_WAT)/Dt:',dydx%veg_water(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' H:            ',h
                  write(unit=*,fmt='(a)')           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
               elseif (.not. record_err) then
                  return
               end if
            end if
            if (y%veg_water(ico) < rk4min_veg_lwater * y%tai(ico)) then
               reject_step = .true.
               if(record_err) cflag1 = .true.
               if (print_problems) then
                  write(unit=*,fmt='(a)')           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                  write(unit=*,fmt='(a)')           ' + Leaf water is off-track...'
                  write(unit=*,fmt='(a)')           '--------------------------------------'
                  write(unit=*,fmt='(a,1x,i6)')     ' PFT:          ',cpatch%pft(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' LAI:          ',y%lai(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' WPA:          ',y%wpa(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' TAI:          ',y%tai(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' HCAPVEG:      ',y%hcapveg(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' VEG_TEMP:     ',y%veg_temp(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' VEG_FRACLIQ:  ',y%veg_fliq(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' VEG_ENERGY:   ',y%veg_energy(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' VEG_WATER:    ',y%veg_water(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' D(VEG_EN)/Dt: ',dydx%veg_energy(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' D(VEG_WAT)/Dt:',dydx%veg_water(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' H:            ',h
                  write(unit=*,fmt='(a)')           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
               elseif (.not. record_err) then
                  return
               end if
            end if
         end if      

      end do
      if(record_err .and. cflag1) integ_err(46,2) = integ_err(46,2) + 1_8
      !------------------------------------------------------------------------------------!


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
                                       , rk4met        ! ! intent(in)
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
      do k=rk4met%lsl,nzg
         write(unit=*,fmt='(i5,3(1x,es12.4))') &
              k, y%soil_tempk(k), y%soil_fracliq(k), y%soil_water(k)
      end do
      write(unit=*,fmt='(78a)') ('-',k=1,78)

      write(unit=*,fmt='(a)') ' '
      write(unit=*,fmt='(78a)') ('-',k=1,78)
      write(unit=*,fmt='(a5,3(1x,a12))') 'LEVEL','  OLD_SOIL_T','OLD_SOIL_FLQ'             &
                                               &,'OLD_SOIL_H2O'
      do k=rk4met%lsl,nzg
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
      write (unit=*,fmt='(2(a5,1x),7(a12,1x))')                                            &
         '  COH','  PFT','         LAI','         WPA','         TAI','  VEG_ENERGY'       &
                        ,'OLD_VEG_ENER','    VEG_TEMP','OLD_VEG_TEMP'
      do ico = 1,cpatch%ncohorts
         if(y%solvable(ico)) then
            write(unit=*,fmt='(2(i5,1x),7(es12.4,1x))')                                    &
               ico,cpatch%pft(ico),y%lai(ico),y%wpa(ico),y%tai(ico)                        &
                  ,y%veg_energy(ico),cpatch%veg_energy(ico),y%veg_temp(ico)                &
                  ,cpatch%veg_temp(ico)
         end if
      end do
      write(unit=*,fmt='(78a)') ('-',k=1,78)

      write(unit=*,fmt='(a)') ' '
      write(unit=*,fmt='(78a)') ('-',k=1,78)
      write (unit=*,fmt='(2(a5,1x),8(a12,1x))') &
         '  COH','  PFT','         LAI','         WPA','         TAI'                      &
                        ,'   VEG_WATER',' OLD_VEG_H2O','    HEAT_CAP','RK4_HEAT_CAP'       &
                        ,'     FRACLIQ'
      do ico = 1,cpatch%ncohorts
         if(y%solvable(ico)) then
            write(unit=*,fmt='(2(i5,1x),8(es12.4,1x))')                                    &
               ico,cpatch%pft(ico),y%lai(ico),y%wpa(ico),y%tai(ico)                        &
                  ,y%veg_water(ico),cpatch%veg_water(ico),cpatch%hcapveg(ico)              &
                  ,y%hcapveg(ico),y%veg_fliq(ico)
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
