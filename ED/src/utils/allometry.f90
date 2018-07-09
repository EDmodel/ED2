!==========================================================================================!
!==========================================================================================!
!    This module is a library with several allometric relationships.                       !
!------------------------------------------------------------------------------------------!
module allometry

   implicit none

   contains
   !=======================================================================================!
   !=======================================================================================!
   real function h2dbh(h,ipft)

      use pft_coms    , only : is_tropical & ! intent(in)
                             , is_liana    & ! intent(in)
                             , b1Ht        & ! intent(in), lookup table
                             , b2Ht        & ! intent(in), lookup table
                             , hgt_ref     ! ! intent(in)
      use ed_misc_coms, only : iallom      ! ! intent(in)
      implicit none

      !----- Arguments --------------------------------------------------------------------!
      real   , intent(in) :: h
      integer, intent(in) :: ipft
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !      Select which type of model we are running.                                    !
      !------------------------------------------------------------------------------------!
      !----- Size- and age-structure (typical ED model). ----------------------------------!
      if (is_liana(ipft)) then
         !----- Inverse Weibull function, similar to Poorter et al. (2006). ---------------!
         h2dbh =  ( log(hgt_ref(ipft) / ( hgt_ref(ipft) - h ) ) / b1Ht(ipft) )             &
               ** ( 1.0 / b2Ht(ipft) )
         !---------------------------------------------------------------------------------!
      elseif (is_tropical(ipft)) then
         select case (iallom)
            case (0,1)
               !----- Default ED-2.1 allometry. -------------------------------------------!
               h2dbh = exp((log(h)-b1Ht(ipft))/b2Ht(ipft))
            case default
               !----- Poorter et al. (2006) allometry. ------------------------------------!
               h2dbh =  ( log(hgt_ref(ipft) / ( hgt_ref(ipft) - h ) ) / b1Ht(ipft) )       &
                  ** ( 1.0 / b2Ht(ipft) )
         end select
      else ! Temperate
         h2dbh = log(1.0-(h-hgt_ref(ipft))/b1Ht(ipft))/b2Ht(ipft)
      end if
      !------------------------------------------------------------------------------------!

      return
   end function h2dbh
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   real function dbh2h(ipft, dbh)
      use pft_coms     , only : is_tropical & ! intent(in)
                              , is_liana    & ! intent(in)
                              , dbh_crit    & ! intent(in)
                              , b1Ht        & ! intent(in)
                              , b2Ht        & ! intent(in)
                              , hgt_ref     & ! intent(in)
                              , hgt_max     ! ! intent(in)
      use ed_misc_coms , only : iallom      & ! intent(in)
                              , ibigleaf    ! ! intent(in)
      use ed_state_vars, only : patchtype   ! ! structure
      use consts_coms  , only : lnexp_max   & ! intent(in)
                              , lnexp_min   ! ! intent(in)
      implicit none

      !----- Arguments --------------------------------------------------------------------!
      integer       , intent(in) :: ipft
      real          , intent(in) :: dbh
      !----- Local variables --------------------------------------------------------------!
      real                :: mdbh
      real                :: lnexp
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !      Select which type of model we are running.                                    !
      !------------------------------------------------------------------------------------!
      select case (ibigleaf)
         case (0)
            !----- Size- and age-structure (typical ED model). ----------------------------!
            if (is_liana(ipft)) then
               !----- Weibull function, similar to Poorter et al. (2006). -----------------!
               mdbh = min(dbh,dbh_crit(ipft))
               dbh2h = hgt_ref(ipft) * (1. - exp(-b1Ht(ipft) * mdbh ** b2Ht(ipft)))
            elseif (is_tropical(ipft)) then
               mdbh = min(dbh,dbh_crit(ipft))
               select case (iallom)
               case (0,1)
                  !----- Default ED-2.1 allometry. ----------------------------------------!
                  dbh2h = exp (b1Ht(ipft) + b2Ht(ipft) * log(mdbh) )
               case default
                  !----- Weibull function. ------------------------------------------------!
                  lnexp = max(lnexp_min,min(lnexp_max,b1Ht(ipft) * mdbh ** b2Ht(ipft)))
                  dbh2h = hgt_ref(ipft) * (1. - exp(-lnexp))
               end select
            else !----- Temperate PFT allometry. ------------------------------------------!
               lnexp = max(lnexp_min,min(lnexp_max,b2Ht(ipft) * dbh))
               dbh2h = hgt_ref(ipft) + b1Ht(ipft) * (1.0 - exp(lnexp))
            end if

         case (1)
            !------------------------------------------------------------------------------!
            !     Big-leaf version of ED. DBH is not really meaningful, but in the big-leaf!
            ! model the typical allometry doesn't really make sense so we impose maximum   !
            ! height.                                                                      !
            !------------------------------------------------------------------------------!
            dbh2h = hgt_max(ipft)
         !---------------------------------------------------------------------------------!
      end select
      !------------------------------------------------------------------------------------!

      return
   end function dbh2h
   !=======================================================================================!
   !=======================================================================================!




   !=======================================================================================!
   !=======================================================================================!
   !     Function that finds Bdead from DBH.                                               !
   !---------------------------------------------------------------------------------------!
   real function size2bd(dbh,hite,ipft)

      use pft_coms    , only : C2B         & ! intent(in)
                             , dbh_crit    & ! intent(in)
                             , b1Bs_small  & ! intent(in)
                             , b2Bs_small  & ! intent(in)
                             , b1Bs_large  & ! intent(in)
                             , b2Bs_large  & ! intent(in)
                             , is_grass    & ! intent(in)
                             , is_tropical & ! intent(in)
                             , is_liana    ! ! intent(in)
      use ed_misc_coms, only : igrass      & ! intent(in)
                             , iallom      ! ! intent(in)
      implicit none

      !----- Arguments --------------------------------------------------------------------!
      real   , intent(in) :: dbh
      real   , intent(in) :: hite
      integer, intent(in) :: ipft
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Structural biomass depends on the allometry and the size.                      !
      !------------------------------------------------------------------------------------!
      if (igrass == 1 .and. is_grass(ipft)   ) then
         size2bd = 0.0
      else if (iallom == 3 .and. is_tropical(ipft) .and. (.not. is_liana(ipft))) then
         !----- Decide parameters based on seedling/adult size. ---------------------------!
         if (dbh <= dbh_crit(ipft)) then
            size2bd = b1Bs_small(ipft) / C2B * (dbh*dbh*hite) ** b2Bs_small(ipft)
         else
            size2bd = b1Bs_large(ipft) / C2B * (dbh*dbh*hite) ** b2Bs_large(ipft)
         end if
         !---------------------------------------------------------------------------------!
      else 
         !----- Decide parameters baded on the maximum height. ----------------------------!
         if (dbh <= dbh_crit(ipft)) then
            size2bd = b1Bs_small(ipft) / C2B * dbh ** b2Bs_small(ipft)
         else
            size2bd = b1Bs_large(ipft) / C2B * dbh ** b2Bs_large(ipft)
         end if
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!

      return
   end function size2bd
   !=======================================================================================!
   !=======================================================================================!




   !=======================================================================================!
   !=======================================================================================!
   !     This function computes the diameter at the breast height given the structural     !
   ! (dead) biomass and the PFT type.  In addition to these inputs, it is also required to !
   ! give a first guess for DBH, which can be the previous value, for example.   This is   !
   ! only used for the new allometry, because we must solve an iterative root-finding      !
   ! method.                                                                               !
   !---------------------------------------------------------------------------------------!
   real function bd2dbh(ipft, bdead)
      use pft_coms    , only : b1Bs_small  & ! intent(in)
                             , b2Bs_small  & ! intent(in)
                             , b1Bs_large  & ! intent(in)
                             , b2Bs_large  & ! intent(in)
                             , bdead_crit  & ! intent(in)
                             , hgt_max     & ! intent(in)
                             , dbh_lut     & ! intent(in)
                             , bdead_lut   & ! intent(in)
                             , le_mask_lut & ! intent(out)
                             , ge_mask_lut & ! intent(out)
                             , is_tropical & ! intent(in)
                             , is_liana    & ! intent(in)
                             , C2B         ! ! intent(in)
      use ed_misc_coms, only : iallom      ! ! intent(in)
      use consts_coms , only : lnexp_min   & ! intent(in)
                             , lnexp_max   ! ! intent(in)
      use therm_lib   , only : toler       & ! intent(in)
                             , maxfpo      ! ! intent(in)
      implicit none

      !----- Arguments --------------------------------------------------------------------!
      integer, intent(in)  :: ipft      ! PFT type                           [         ---]
      real   , intent(in)  :: bdead     ! Structural (dead) biomass          [   kgC/plant]
      !----- Local variables. -------------------------------------------------------------!
      integer              :: ilwr      ! Lower index of the lookup table
      integer              :: iupr      ! Upper index of the lookup table
      integer              :: it        ! Iteration counter
      real                 :: dbha      ! Lower guess for DBH
      real                 :: dbhz      ! Upper guess for DBH
      real                 :: hgta      ! Lower guess for height
      real                 :: hite      ! New   guess for height
      real                 :: hgtz      ! Upper guess for height
      real                 :: funa      ! Function evaluation for lower guess
      real                 :: funz      ! Function evaluation for upper guess
      real                 :: fun       ! Function evaluation for new   guess
      logical              :: zside     ! Last update was on zside
      logical              :: converged ! Iterative method converged
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Decide which coefficients to use based on the critical bdead.                   !
      !------------------------------------------------------------------------------------!
      if (iallom == 3 .and. is_tropical(ipft) .and. (.not. is_liana(ipft)) ) then
         if (bdead <= bdead_lut(1,ipft)) then
            !----- Use the look-up table to find the best dbh. ----------------------------!
            bd2dbh = dbh_lut(1,ipft) * bdead / bdead_lut(1,ipft)
            !------------------------------------------------------------------------------!
         else if (bdead >= bdead_crit(ipft)) then
            !------------------------------------------------------------------------------!
            !     Bdead is above critical value, height is known.                          !
            !------------------------------------------------------------------------------!
            bd2dbh =  ( C2B * bdead                                                        &
                      / ( b1Bs_large(ipft) * hgt_max(ipft)**b2Bs_large(ipft) ))            &
                   ** ( 1. / (2. * b2Bs_large(ipft) ) )
            !------------------------------------------------------------------------------!
         else
            !----- Use the look-up table to find the best dbh. ----------------------------!
            le_mask_lut(:) = bdead <= bdead_lut(:,ipft)
            ge_mask_lut(:) = bdead >= bdead_lut(:,ipft)
            ilwr  = maxloc (bdead_lut(:,ipft),dim=1,mask=le_mask_lut)
            iupr  = minloc (bdead_lut(:,ipft),dim=1,mask=ge_mask_lut)
            !------------------------------------------------------------------------------!

            !------------------------------------------------------------------------------!
            !      In case ilwr and iupr are the same, we have an exact estimate.  Other-  !
            ! wise, use log-linear interpolation.                                          !
            !------------------------------------------------------------------------------!
            if (ilwr == iupr) then
               bd2dbh = dbh_lut(ilwr,ipft)
            else
               !------ Define the first guess for Regula Falsi (Illinois) method. ---------!
               dbha    = dbh_lut(ilwr,ipft)
               dbhz    = dbh_lut(iupr,ipft)
               hgta    = dbh2h(ipft,dbha)
               hgtz    = dbh2h(ipft,dbhz)
               funa    = size2bd(dbha,hgta,ipft) - bdead
               funz    = size2bd(dbhz,hgtz,ipft) - bdead
               zside   = funa * funz < 0.
               !---------------------------------------------------------------------------!


               !----- This should not happen, but check it anyway. ------------------------!
               if (.not. zside) then
                  write (unit=*,fmt='(a)')           '------------------------------------'
                  write (unit=*,fmt='(a)')           ' Ill-posed first guesses:           '
                  write (unit=*,fmt='(a)')           '------------------------------------'
                  write (unit=*,fmt='(a)')           ' '
                  write (unit=*,fmt='(a)')           ' Input:    '
                  write (unit=*,fmt='(a,1x,i14)'   ) ' + ipft    =',ipft
                  write (unit=*,fmt='(a,1x,es14.7)') ' + bdead   =',bdead
                  write (unit=*,fmt='(a)')           ' '
                  write (unit=*,fmt='(a)')           ' Guesses:  '
                  write (unit=*,fmt='(a,1x,es14.7)') ' + dbha    =',dbha
                  write (unit=*,fmt='(a,1x,es14.7)') ' + dbhz    =',dbhz
                  write (unit=*,fmt='(a,1x,es14.7)') ' + hgta    =',hgta
                  write (unit=*,fmt='(a,1x,es14.7)') ' + hgtz    =',hgtz
                  write (unit=*,fmt='(a,1x,es14.7)') ' + bdeada  =',size2bd(dbha,hgta,ipft)
                  write (unit=*,fmt='(a,1x,es14.7)') ' + bdeadz  =',size2bd(dbhz,hgtz,ipft)
                  write (unit=*,fmt='(a,1x,es14.7)') ' + funa    =',funa
                  write (unit=*,fmt='(a,1x,es14.7)') ' + funz    =',funz
                  write (unit=*,fmt='(a)')           ' '
                  write (unit=*,fmt='(a)')           '-----------------------------------'
                  call fatal_error('Ill-posed first guesses for Regula Falsi'              &
                                  ,'bd2dbh','allometry.f90')
               end if
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !     Loop until convergence.                                               !
               !---------------------------------------------------------------------------!
               rfaloop: do it=1,maxfpo
                  bd2dbh = (funz * dbha - funa * dbhz) / ( funz - funa)
                  hite   = dbh2h(ipft,bd2dbh)

                  !------------------------------------------------------------------------!
                  !     Now that we updated the guess, check whether they are really       !
                  ! close.  If so, it converged within tolerance.                          !
                  !------------------------------------------------------------------------!
                  converged = abs(bd2dbh-dbha) < toler * bd2dbh
                  if (converged) exit rfaloop
                  !------------------------------------------------------------------------!


                  !------ Find the new function evaluation. -------------------------------!
                  fun  =  size2bd(bd2dbh,hite,ipft) - bdead
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !     Define the new interval based on the intermediate value theorem.   !
                  ! We check whether the same side is being picked repeatedly, and apply   !
                  ! the Regula Falsi step to try to speed up convergence.                  !
                  !------------------------------------------------------------------------!
                  if (fun*funa < 0. ) then
                     dbhz  = bd2dbh
                     funz  = fun
                     if (zside) funa = funa * 0.5
                     zside = .true.
                  else
                     dbha  = bd2dbh
                     funa  = fun
                     if (.not. zside) funz = funz * 0.5
                     zside = .false.
                  end if
                  !------------------------------------------------------------------------!
               end do rfaloop
               !---------------------------------------------------------------------------!




               !---------------------------------------------------------------------------!
               !      This method should always congerge in case the function is           !
               ! continuous, but just in case.                                             !
               !---------------------------------------------------------------------------!
               if (.not. converged) then
                  write (unit=*,fmt='(a)')           '------------------------------------'
                  write (unit=*,fmt='(a)')           ' Function did not converge:         '
                  write (unit=*,fmt='(a)')           '------------------------------------'
                  write (unit=*,fmt='(a)')           ' '
                  write (unit=*,fmt='(a,1x,i14)'   ) ' Iterations =',maxfpo
                  write (unit=*,fmt='(a)')           ' '
                  write (unit=*,fmt='(a)')           ' Input:     '
                  write (unit=*,fmt='(a,1x,i14)'   ) ' + ipft    =',ipft
                  write (unit=*,fmt='(a,1x,es14.7)') ' + bdead   =',bdead
                  write (unit=*,fmt='(a)')           ' '
                  write (unit=*,fmt='(a)')           ' Guesses:  '
                  write (unit=*,fmt='(a,1x,es14.7)') ' + dbha    =',dbha
                  write (unit=*,fmt='(a,1x,es14.7)') ' + dbh     =',bd2dbh
                  write (unit=*,fmt='(a,1x,es14.7)') ' + dbhz    =',dbhz
                  write (unit=*,fmt='(a,1x,es14.7)') ' + hgta    =',hgta
                  write (unit=*,fmt='(a,1x,es14.7)') ' + hite    =',hite
                  write (unit=*,fmt='(a,1x,es14.7)') ' + hgtz    =',hgtz
                  write (unit=*,fmt='(a,1x,es14.7)') ' + bdeada  =',size2bd(dbha,hgta,ipft)
                  write (unit=*,fmt='(a,1x,es14.7)') ' + bdeadz  =',size2bd(dbhz,hgtz,ipft)
                  write (unit=*,fmt='(a,1x,es14.7)') ' + funa    =',funa
                  write (unit=*,fmt='(a,1x,es14.7)') ' + fun     =',fun
                  write (unit=*,fmt='(a,1x,es14.7)') ' + funz    =',funz
                  write (unit=*,fmt='(a)')           ' '
                  write (unit=*,fmt='(a)')           '-------------------------------------'
                  call fatal_error('Ill-posed first guess for Regula Falsi'                &
                                  ,'bd2dbh','allometry.f90')
               end if
               !------------------------------------------------------------------------------!
            end if
            !------------------------------------------------------------------------------!
         end if
         !---------------------------------------------------------------------------------!
      else
         !----- Bdead is not a direct function of height, simple inversion is sufficient. -!
         if (bdead <= bdead_crit(ipft)) then
            bd2dbh = (bdead / b1Bs_small(ipft) * C2B)**(1.0/b2Bs_small(ipft))
         else
            bd2dbh = (bdead / b1Bs_large(ipft) * C2B)**(1.0/b2Bs_large(ipft))
         end if
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!

      return
   end function bd2dbh
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !      This function computes the maximum leaf biomass [kgC/m2] given the size (either  !
   ! height or DBH).  Trees would always use DBH as the starting point, but grasses may    !
   ! use DBH (old style) or height (new style).  This replaces dbh2bl and h2bl with a      !
   ! single generic function that should be used by all plants.                            !
   !---------------------------------------------------------------------------------------!
   real function size2bl(dbh,hite,ipft)
      use pft_coms     , only : dbh_crit       & ! intent(in)
                              , C2B            & ! intent(in)
                              , b1Bl           & ! intent(in)
                              , b2Bl           & ! intent(in)
                              , is_liana       & ! intent(in)
                              , is_grass       & ! intent(in)
                              , is_tropical    & ! intent(in)
                              , liana_dbh_crit ! ! intent(in)
      use ed_misc_coms , only : igrass         & ! intent(in)
                              , iallom         ! ! intent(in)
      use ed_state_vars, only : patchtype      ! ! structure
      implicit none

      !----- Arguments --------------------------------------------------------------------!
      real          , intent(in) :: dbh
      real          , intent(in) :: hite
      integer       , intent(in) :: ipft
      !----- Local variables --------------------------------------------------------------!
      real                       :: mdbh
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Get the DBH, or potential DBH in case of grasses.                              !
      !------------------------------------------------------------------------------------!
      if (igrass == 1 .and. is_grass(ipft)) then 
         !---- Use height for new grasses. ------------------------------------------------!
         mdbh   = min(dbh,h2dbh(hite,ipft))
      elseif (is_liana(ipft)) then
         mdbh   = min(dbh,liana_dbh_crit)
      else
         !---- Use dbh for trees. ---------------------------------------------------------!
         mdbh   = min(dbh,dbh_crit(ipft))
      end if
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !     Find leaf biomass depending on the allometry.  The new allometry uses dbh and  !
      ! height, whereas the old allometry uses dbh only.                                   !
      !------------------------------------------------------------------------------------!
      if (iallom == 3 .and. is_tropical(ipft) .and. (.not. is_liana(ipft))) then
         size2bl = b1Bl(ipft) / C2B * (mdbh*mdbh*hite) ** b2Bl(ipft)
      else 
         size2bl = b1Bl(ipft) / C2B * mdbh ** b2Bl(ipft)
      end if
      !------------------------------------------------------------------------------------!


      return
   end function size2bl
   !=======================================================================================!
   !=======================================================================================!





   !=======================================================================================!
   !=======================================================================================!
   !      This subroutine finds height given the biomass of living tissues, assuming       !
   ! minimum sapwood biomass given DBH.                                                    !
   !---------------------------------------------------------------------------------------!
   real function ba2h(balive,ipft)
      use pft_coms    , only : hgt_min     & ! intent(in)
                             , hgt_max     & ! intent(in)
                             , balive_crit & ! intent(in)
                             , dbh_lut     & ! intent(in)
                             , balive_lut  & ! intent(in)
                             , le_mask_lut & ! intent(out)
                             , ge_mask_lut ! ! intent(out)
      use consts_coms , only : lnexp_min   & ! intent(in)
                             , lnexp_max   ! ! intent(in)
      implicit none

      !----- Arguments --------------------------------------------------------------------!
      real   , intent(in) :: balive      ! Live biomass
      integer, intent(in) :: ipft        ! PFT
      !----- Local variables. -------------------------------------------------------------!
      integer             :: ilwr        ! Lower index of the lookup table
      integer             :: iupr        ! Upper index of the lookup table
      real                :: dbh         ! Best guess for DBH.
      real                :: finterp     ! Interpolation factor
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !    In case there is more biomass than balive_crit, assume that the plant is at     !
      ! maximum height.                                                                    !
      !------------------------------------------------------------------------------------!
      if (balive >= balive_crit(ipft)) then
         ba2h        = hgt_max(ipft)
      elseif (balive <= balive_lut(1,ipft)) then
         !----- Do not let height be less than minimum. -----------------------------------!
         ba2h        = hgt_min(ipft)
         !---------------------------------------------------------------------------------!
      else
         !----- Use the look-up table to find the best dbh. -------------------------------!
         le_mask_lut(:) = balive <= balive_lut(:,ipft)
         ge_mask_lut(:) = balive >= balive_lut(:,ipft)
         ilwr  = maxloc (balive_lut(:,ipft),dim=1,mask=le_mask_lut)
         iupr  = minloc (balive_lut(:,ipft),dim=1,mask=ge_mask_lut)
         !---------------------------------------------------------------------------------!

         !---------------------------------------------------------------------------------!
         !      In case ilwr and iupr are the same, we have an exact estimate.  Otherwise, !
         ! use log-linear interpolation.                                                   !
         !---------------------------------------------------------------------------------!
         if (ilwr == iupr) then
            dbh     = dbh_lut(ilwr,ipft)
         else
            finterp = log( dbh_lut   (iupr,ipft) / dbh_lut   (ilwr,ipft))                  &
                    * log( balive                / balive_lut(ilwr,ipft))                  &
                    / log( balive_lut(iupr,ipft) / balive_lut(ilwr,ipft))
            finterp = max(lnexp_min,min(lnexp_max,finterp))
            dbh     = dbh_lut(ilwr,ipft) * exp(finterp)
         end if
         ba2h = dbh2h(ipft,dbh)
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!

      return
    end function ba2h
   !=======================================================================================!
   !=======================================================================================!




   !=======================================================================================!
   !=======================================================================================!
   !      This function determines an effective DBH for grasses given their leaf biomass.  !
   ! DBH has no real meaning for grasses with the new allometry.                           !
   !---------------------------------------------------------------------------------------!
   real function bl2dbh(bleaf,ipft)
      use pft_coms,     only : dbh_crit    & ! intent(in)
                             , hgt_max     & ! intent(in)
                             , is_grass    & ! intent(in)
                             , C2B         & ! intent(in)
                             , bleaf_crit  & ! intent(in)
                             , b1Bl        & ! intent(in)
                             , b2Bl        & ! intent(in)
                             , dbh_lut     & ! intent(in)
                             , bleaf_lut   & ! intent(in)
                             , le_mask_lut & ! intent(out)
                             , ge_mask_lut & ! intent(out)
                             , is_tropical & ! intent(in)
                             , is_liana    ! ! intent(in)
      use ed_misc_coms, only : igrass      & ! intent(in)
                             , iallom      ! ! intent(in)
      use consts_coms , only : lnexp_min   & ! intent(in)
                             , lnexp_max   ! ! intent(in)
      implicit none

      !----- Arguments --------------------------------------------------------------------!
      real   , intent(in) :: bleaf
      integer, intent(in) :: ipft
      !----- Local variables --------------------------------------------------------------!
      real                :: mdbh
      integer             :: ilwr        ! Lower index of the lookup table
      integer             :: iupr        ! Upper index of the lookup table
      real                :: finterp     ! Interpolation factor
      !------------------------------------------------------------------------------------!


      if ( iallom == 3 .and. is_tropical(ipft) .and. (.not. is_liana(ipft)) ) then
         if (bleaf <= bleaf_lut(1,ipft)) then
            !----- Use the look-up table to find the best dbh. ----------------------------!
            bl2dbh = dbh_lut(1,ipft) * bleaf / bleaf_lut(1,ipft)
            !------------------------------------------------------------------------------!
         else if (bleaf >= bleaf_crit(ipft)) then
            !----- Bleaf should not exceed bleaf_crit.  Set dbh to maximum. ---------------!
            bl2dbh = dbh_crit(ipft)
            !------------------------------------------------------------------------------!
         else
            !----- Use the look-up table to find the best dbh. ----------------------------!
            le_mask_lut(:) = bleaf <= bleaf_lut(:,ipft)
            ge_mask_lut(:) = bleaf >= bleaf_lut(:,ipft)
            ilwr  = maxloc (bleaf_lut(:,ipft),dim=1,mask=le_mask_lut)
            iupr  = minloc (bleaf_lut(:,ipft),dim=1,mask=ge_mask_lut)
            !------------------------------------------------------------------------------!

            !------------------------------------------------------------------------------!
            !      In case ilwr and iupr are the same, we have an exact estimate.  Other-  !
            ! wise, use log-linear interpolation.                                          !
            !------------------------------------------------------------------------------!
            if (ilwr == iupr) then
               bl2dbh = dbh_lut(ilwr,ipft)
            else
               finterp = log( dbh_lut  (iupr,ipft) / dbh_lut  (ilwr,ipft))                 &
                       * log( bleaf                / bleaf_lut(ilwr,ipft))                 &
                       / log( bleaf_lut(iupr,ipft) / bleaf_lut(ilwr,ipft))
               finterp = max(lnexp_min,min(lnexp_max,finterp))
               bl2dbh  = dbh_lut(ilwr,ipft) * exp(finterp)
            end if
            !------------------------------------------------------------------------------!
         end if
         !---------------------------------------------------------------------------------!

      else
         !----- Find out whether this is an adult tree or a sapling/grass. ----------------!
         mdbh = (bleaf * C2B / b1Bl(ipft) ) ** (1./b2Bl(ipft))
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     For grasses, limit maximum effective dbh by maximum height.                 !
         !---------------------------------------------------------------------------------!
         if (is_grass(ipft) .and. igrass == 1) then
            bl2dbh = min(mdbh, h2dbh(hgt_max(ipft),ipft))
         else
            bl2dbh = min(mdbh, dbh_crit(ipft))
         end if
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!

      return
   end function bl2dbh
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !      This function determines the height for grasses given their leaf biomass.        !
   !---------------------------------------------------------------------------------------!
   real function bl2h(bleaf,ipft)
      use pft_coms,      only:  hgt_max    ! ! intent(in), lookup table
      implicit none
      
      !----- Arguments --------------------------------------------------------------------!
      real   , intent(in)      :: bleaf
      integer, intent(in)      :: ipft
      !------------------------------------------------------------------------------------!

      !----- Use existing allometric equations to convert leaves to height. ---------------!
      bl2h = min(hgt_max(ipft),dbh2h(ipft,bl2dbh(bleaf,ipft)))

      return
   end function bl2h
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    Canopy Area allometry from Dietze and Clark (2008).                                !
   !---------------------------------------------------------------------------------------!
   real function size2ca(dbh,hite,sla,ipft,cap_crit)
      use pft_coms     , only : dbh_crit       & ! intent(in)
                              , hgt_max        & ! intent(in)
                              , is_grass       & ! intent(in)
                              , is_liana       & ! intent(in)
                              , is_tropical    & ! intent(in)
                              , b1Ca           & ! intent(in)
                              , b2Ca           & ! intent(in)
                              , liana_dbh_crit ! ! intent(in)
      use ed_misc_coms , only : igrass         & ! intent(in)
                              , iallom         ! ! intent(in)
      use ed_state_vars, only : patchtype      ! ! structure
      use consts_coms  , only : tiny_num       ! ! intent(in)
      implicit none

      !----- Arguments --------------------------------------------------------------------!
      real   , intent(in)           :: dbh       !> Diameter at breast height     [     cm]
      real   , intent(in)           :: hite      !> Height                        [      m]
      real   , intent(in)           :: sla       !> Specific leaf area            [ m2/kgC]
      integer, intent(in)           :: ipft      !> Current PFT                   [     --]
      logical, intent(in), optional :: cap_crit  !> Force cap at dbh_crit         [    T|F]
      !----- Internal variables -----------------------------------------------------------!
      real                          :: loclai    !> The maximum local LAI         [  m2/m2]
      real                          :: mdbh      !> The maximum DBH               [     cm]
      logical                       :: upr_loose !> Cap at dbh_crit               [    T|F]
      !------------------------------------------------------------------------------------!


      !----- Decide whether to cap crown area at critical size. ---------------------------!
      if (present(cap_crit)) then
         upr_loose = .not. cap_crit
      else
         upr_loose = .false.
      end if
      !------------------------------------------------------------------------------------!

      if (dbh < tiny_num) then
         loclai  = 0.0
         size2ca = 0.0
      else

         !----- make this function generic to size, not just dbh. -------------------------!
         loclai = sla * size2bl(dbh,hite,ipft) 
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Check whether to impose a maximum height.                                   !
         !---------------------------------------------------------------------------------!
         select case (iallom)
         case (0)
            !----- No upper bound in the allometry. ---------------------------------------!
            mdbh = dbh
            !------------------------------------------------------------------------------!

         case default
            !----- Check whether this is grass or tree, impose maximum DBH. ---------------!
            if (is_grass(ipft) .and. igrass == 1) then
               mdbh = min(dbh,h2dbh(hgt_max(ipft),ipft) )
            elseif (upr_loose) then
               mdbh = dbh
            elseif (is_liana(ipft)) then
               mdbh = min(dbh,liana_dbh_crit)
            else
               mdbh = min(dbh,dbh_crit(ipft))
            end if
            !------------------------------------------------------------------------------!
         end select
         !---------------------------------------------------------------------------------!


         !----- Find the nominal crown area. ----------------------------------------------!
         if (iallom == 3 .and. is_tropical(ipft) .and. (.not. is_liana(ipft))) then
            size2ca = b1Ca(ipft) * ( mdbh * mdbh * hite ) ** b2Ca(ipft)
         else
            size2ca = b1Ca(ipft) * mdbh ** b2Ca(ipft)
         end if
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!



      !----- Local LAI / Crown area should never be less than one. ------------------------!
      size2ca = min(loclai, size2ca)
      !------------------------------------------------------------------------------------!

      return
   end function size2ca
   !=======================================================================================!
   !=======================================================================================!




   !=======================================================================================!
   !=======================================================================================!
   !    This function computes the commercial volume of a tree (stem/bole volume, in m3).  !
   ! The current equation is a re-fit from Nogueira et al. (2008) so a single set of       !
   ! parameters can be used.  Their equation is for tropical trees only, so temperate and  !
   ! boreal forests may need a different set of parameters.  Grasses are assumed to have   !
   ! no commercial volume.                                                                 !
   !                                                                                       !
   ! Nogueira, E. M., et al.  Estimates of forest biomass in the Brazilian Amazon: new     !
   !    allometric equations and adjustments to biomass from wood-volume inventories.      !
   !    Forest Ecol. Manag., 256(11), 1853-1867, Nov. 2008,                                !
   !    doi:10.1016/j.foreco.2008.07.022.                                                  !
   !---------------------------------------------------------------------------------------!
   real function size2vol(dbh,hgt,ipft)
      use pft_coms    , only : b1Vol    & ! intent(in)
                             , b2Vol    ! ! intent(in)
      implicit none

      !----- Arguments --------------------------------------------------------------------!
      real   , intent(in) :: dbh
      real   , intent(in) :: hgt
      integer, intent(in) :: ipft
      !------------------------------------------------------------------------------------!


      size2vol = b1Vol(ipft) * (dbh * dbh * hgt) ** b2Vol(ipft)

      return
   end function size2vol
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This function computes the commercial timber biomass for different PFTs.           !
   !---------------------------------------------------------------------------------------!
   real function size2bt(dbh,hgt,bdead,bsapwooda,bbark,ipft)
      use pft_coms    , only : is_grass    & ! intent(in)
                             , is_tropical & ! intent(in)
                             , rho         & ! intent(in)
                             , agf_bs      & ! intent(in)
                             , C2B         ! ! intent(in)

      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real(kind=4), intent(in) :: dbh
      real(kind=4), intent(in) :: hgt
      real(kind=4), intent(in) :: bdead
      real(kind=4), intent(in) :: bsapwooda
      real(kind=4), intent(in) :: bbark
      integer     , intent(in) :: ipft
      !----- Local variables. -------------------------------------------------------------!
      real(kind=4)             :: btimber
      real(kind=4)             :: bagwood
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !      Calculate commercial timber biomass based on the life form and habitat.       !
      !------------------------------------------------------------------------------------!
      if (is_grass(ipft)) then
         !---------------------------------------------------------------------------------!
         !     Grasses don't have commercial timber.                                       !
         !---------------------------------------------------------------------------------!
         size2bt = 0.0
         !---------------------------------------------------------------------------------!
      else if (is_tropical(ipft)) then
         !---------------------------------------------------------------------------------!
         !     Tropical trees.  We use the tree volume scaled with typical wood density.   !
         !---------------------------------------------------------------------------------!
         btimber  = 1000. * rho(ipft) * size2vol(dbh,hgt,ipft) / C2B
         bagwood  = agf_bs(ipft) * (bdead + bbark) + bsapwooda
         size2bt  = min(btimber,bagwood)
         !---------------------------------------------------------------------------------!
      else
         !---------------------------------------------------------------------------------!
         !     Temperate trees.  For the time being we use total aboveground wood biomass, !
         ! so results are consistent with former implementation.                           !
         !---------------------------------------------------------------------------------!
         size2bt  = agf_bs(ipft) * (bdead + bbark) + bsapwooda
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!

      return
   end function size2bt
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This function computes bark thickness.  To obtain the actual thickness, we compare !
   ! the actual bark biomass with on-allometry bark biomass.                               !
   !---------------------------------------------------------------------------------------!
   real function size2xb(dbh,hgt,bbark,ipft)
      use pft_coms    , only : is_grass       & ! intent(in)
                             , is_liana       & ! intent(in)
                             , qbark          & ! intent(in)
                             , dbh_crit       & ! intent(in)
                             , hgt_max        & ! intent(in)
                             , liana_dbh_crit & ! intent(in)
                             , b1Xb           ! ! intent(in)
      use consts_coms , only : tiny_num       ! ! intent(in)
      use ed_misc_coms, only : igrass         ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real(kind=4), intent(in) :: dbh
      real(kind=4), intent(in) :: hgt
      real(kind=4), intent(in) :: bbark
      integer     , intent(in) :: ipft
      !----- Local variables. -------------------------------------------------------------!
      real(kind=4)             :: mdbh
      real(kind=4)             :: bleaf_max
      real(kind=4)             :: bbark_max
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !      Calculate commercial timber biomass based on the life form and habitat.       !
      !------------------------------------------------------------------------------------!
      if (qbark(ipft) < tiny_num) then
         size2xb = 0.0
      else
         !----- Check whether this is grass, liana or tree, impose maximum DBH. -----------!
         if (is_grass(ipft) .and. igrass == 1) then
            mdbh = min(dbh,h2dbh(hgt_max(ipft),ipft) )
         elseif (is_liana(ipft)) then
            mdbh = min(dbh,liana_dbh_crit)
         else
            mdbh = min(dbh,dbh_crit(ipft))
         end if
         !---------------------------------------------------------------------------------!


         !----- Find on-allometry bark biomass. -------------------------------------------!
         bleaf_max = size2bl(mdbh,hgt,ipft) 
         bbark_max = qbark(ipft) * hgt * bleaf_max
         !---------------------------------------------------------------------------------!

         !----- Scale bark thickness with the ratio between actual and on-allometry. ------!
         size2xb   = b1Xb(ipft) * mdbh * bbark / bbark_max
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!

      return
   end function size2xb
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This function find the potential rooting depth (i.e. based only on allometry, and  !
   ! ignoring soil depth.).                                                                !
   !---------------------------------------------------------------------------------------!
   real function size2prd(hite,dbh,ipft)
      use ed_misc_coms, only : iallom          & ! intent(in)
                             , use_efrd_trtree ! ! intent(in)
      use pft_coms    , only : is_tropical     & ! intent(in)
                             , is_grass        & ! intent(in)
                             , is_liana        & ! intent(in)
                             , dbh_crit        & ! intent(in)
                             , b1Rd            & ! intent(in)
                             , b2Rd            & ! intent(in)
                             , d18O_ref        & ! intent(in)
                             , b1d18O          & ! intent(in)
                             , b2d18O          & ! intent(in)
                             , b1Efrd          & ! intent(in)
                             , b2Efrd          ! ! intent(in)
      implicit none

      !----- Arguments --------------------------------------------------------------------!
      real   , intent(in) :: hite
      real   , intent(in) :: dbh
      integer, intent(in) :: ipft
      !----- Local variables --------------------------------------------------------------!
      logical             :: is_efrd
      real                :: dbhuse
      real                :: size
      real                :: d18O
      !------------------------------------------------------------------------------------!



      !----- Decide whether to use the delta-18O based allometry or size-dependent only. --!
      is_efrd = use_efrd_trtree        .and. is_tropical(ipft)      .and.                  &
                (.not. is_liana(ipft)) .and. (.not. is_grass(ipft))
      !------------------------------------------------------------------------------------!



      if (is_efrd) then
         !---------------------------------------------------------------------------------!
         !    For tropical trees, use the allometric model to obtain the Effective         !
         ! Functional Rooting Depth based on B18.  We made a slight modification in their  !
         ! equation relating delta 18O and depth:                                          !
         !                                                                                 !
         !    depth = exp(a + b * d18O^2)                                                  !
         !                                                                                 !
         ! because it fits the data better than the original equation without the square,  !
         ! and it avoids extremely shallow soils for small trees.  We also use a           !
         ! heteroscedastic least squares, using the algorithm developed by L16.            !
         !                                                                                 !
         ! References:                                                                     !
         !                                                                                 !
         ! Brum M, Vadeboncoeur MA, Ivanov V, Asbjornsen H, Saleska S, Alves LF, Penha D,  !
         !    Dias JD, Aragao LEOC, Barros F, Bittencourt P, Pereira L, Oliveira RS, 2018. !
         !    Hydrological niche segregation defines forest structure and drought          !
         !    tolerance strategies in a seasonal Amazonian forest. J. Ecol., in press.     !
         !    doi:10.1111/1365-2745.13022 (B18).                                           !
         !                                                                                 !
         ! Longo M, Keller M, dos-Santos MN, Leitold V, Pinage ER, Baccini A, Saatchi S,   !
         !    Nogueira EM, Batistella M , Morton DC. 2016. Aboveground biomass variability !
         !    across intact and degraded forests in the Brazilian Amazon.                  !
         !    Global Biogeochem. Cycles, 30(11):1639-1660. doi:10.1002/2016GB005465 (L16). !
         !---------------------------------------------------------------------------------!
         dbhuse   = min(dbh_crit(ipft),dbh)
         d18O     = d18O_ref(ipft) * (1. - exp(- b1d18O(ipft) * dbhuse ** b2d18O(ipft)))
         size2prd = - exp(b1Efrd(ipft) + b2Efrd(ipft) * d18O * d18O)
         !---------------------------------------------------------------------------------!
      else
         !---------------------------------------------------------------------------------!
         !     Standard size-dependent rooting depth.                                      !
         !---------------------------------------------------------------------------------!
         select case (iallom)
         case (0)
            !---------------------------------------------------------------------------------!
            !    Original ED-2.1 (I don't know the source for this equation, though).         !
            ! Grasses get a fixed rooting depth of 70cm.                                      !
            !---------------------------------------------------------------------------------!
            size     = dbh * dbh * hite
            !---------------------------------------------------------------------------------!
         case (1,2)
            !---------------------------------------------------------------------------------!
            !    This is just a test allometry, that imposes root depth to be 0.5 m for       !
            ! plants that are 0.15-m tall, and 5.0 m for plants that are 35-m tall.           !
            !---------------------------------------------------------------------------------!
            size     = hite
            !---------------------------------------------------------------------------------!
         case (3)
            !---------------------------------------------------------------------------------!
            !    For tropical trees, use size-based parameters loosely based on B18, applying !
            ! a simple log-linear interpolation as a function of size (dbh*dbh*hgt) that      !
            ! makes the minimum and maximum rooting depths match the values predicted with    !
            ! the modified B18 equations.  For grasses, we also use a log-linear size inter-  !
            ! polation that makes their seedling rooting depth the same as seedling trees,    !
            ! and puts the depth at 1.5m when they are at maximum height.  We use the same    !
            ! coefficients and functional form as iallom=2 for non-tropical plants and for tropical        !
            ! lianas.                                                                         !
            !                                                                                 !
            !    depth = exp(a + b * d18O^2)                                                  !
            !                                                                                 !
            ! because it fits the data better than the original (no square for d18O), and     !
            ! avoids extremely shallow soils for small trees.  We also use a heteroscedastic  !
            ! least squares, using the algorithm developed by L16.  For grasses, we applied   !
            ! a simple log-linear interpolation as a function of size (dbh*dbh*hgt) that      !
            ! makes their seedling rooting depth the same as seedling trees, and puts the     !
            ! depth at 1.5m when they are at maximum height.  We use the same coefficients    !
            ! and functional form as iallom=2 for non-tropical plants and for tropical        !
            ! lianas.                                                                         !
            !                                                                                 !
            ! References:                                                                     !
            !                                                                                 !
            ! Brum M, Vadeboncoeur MA, Ivanov V, Asbjornsen H, Saleska S, Alves LF, Penha D,  !
            !    Dias JD, Aragao LEOC, Barros F, Bittencourt P, Pereira L, Oliveira RS, 2018. !
            !    Hydrological niche segregation defines forest structure and drought          !
            !    tolerance strategies in a seasonal Amazonian forest. J. Ecol., in press.     !
            !    doi:10.1111/1365-2745.13022 (B18).                                           !
            !                                                                                 !
            ! Longo M, Keller M, dos-Santos MN, Leitold V, Pinage ER, Baccini A, Saatchi S,   !
            !    Nogueira EM, Batistella M , Morton DC. 2016. Aboveground biomass variability !
            !    across intact and degraded forests in the Brazilian Amazon.                  !
            !    Global Biogeochem. Cycles, 30(11):1639-1660. doi:10.1002/2016GB005465 (L16). !
            !---------------------------------------------------------------------------------!
            if ( is_tropical(ipft) .and. (.not. is_liana(ipft)) ) then
               dbhuse   = min(dbh_crit(ipft),dbh)
               size     = dbhuse * dbhuse * hite
            else
               size     = hite
            end if
            !---------------------------------------------------------------------------------!
         end select
         !------------------------------------------------------------------------------------!



         !----- Find rooting depth based on selected size variable. -----------------------!
         size2prd = b1Rd(ipft)  * size ** b2Rd(ipft)
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!

      return
   end function size2prd
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function finds the actual rooting depth, which mlimited by soil depth.       !
   !---------------------------------------------------------------------------------------!
   integer function size2krdepth(hite,dbh,ipft,lsl)
      use ed_misc_coms, only : iallom      ! ! intent(in)
      use grid_coms   , only : nzg         ! ! intent(in)
      use soil_coms   , only : slz         ! ! intent(in)
      implicit none

      !----- Arguments --------------------------------------------------------------------!
      real   , intent(in) :: hite
      real   , intent(in) :: dbh
      integer, intent(in) :: ipft
      integer, intent(in) :: lsl
      !----- Local variables --------------------------------------------------------------!
      real                :: pot_root_depth
      integer             :: k
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Find the potential rooting depth, which is only based on allometric equations.  !
      !------------------------------------------------------------------------------------!
      pot_root_depth = size2prd(hite,dbh,ipft)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Root depth is the maximum root depth if the soil is that deep.  Find what is   !
      ! the deepest soil layer this root can go.                                           !
      !------------------------------------------------------------------------------------!
      size2krdepth = nzg
      do k=nzg,lsl+1,-1
         if (pot_root_depth < slz(k)) size2krdepth = k-1
      end do
      !------------------------------------------------------------------------------------!

      return
   end function size2krdepth
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This function finds the trunk height.  Currently this is based on the following    !
   ! reference, which is for a site in Bolivia:                                            !
   !                                                                                       !
   ! Poorter L., L. Bongers, F. Bongers, 2006: Architecture of 54 moist-forest tree        !
   !     species: traits, trade-offs, and functional groups. Ecology, 87, 1289-1301.       !
   !                                                                                       !
   ! Crown length is defined as the height of the tree minus the height of the first leaf  !
   !---------------------------------------------------------------------------------------!
   real function h2crownbh(height,ipft)
      use pft_coms, only : b1Cl & ! intent(in)
                         , b2Cl ! ! intent(in)
      implicit none 

      !----- Arguments --------------------------------------------------------------------!
      real   , intent(in) :: height
      integer, intent(in) :: ipft
      !----- Local variables. -------------------------------------------------------------!
      real                :: crown_length
      !------------------------------------------------------------------------------------!
      
      crown_length = b1Cl(ipft) * height ** b2Cl(ipft)
      h2crownbh    = max(0.05,height - crown_length)

      return
   end function h2crownbh
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine finds the total above ground biomass (wood + leaves)              !
   !---------------------------------------------------------------------------------------!
   real function ed_biomass(cpatch,ico)
      use pft_coms,      only : agf_bs     ! ! intent(in)
      use ed_state_vars, only : patchtype  ! ! Structure
      implicit none 

      !----- Arguments --------------------------------------------------------------------!
      type(patchtype), target :: cpatch
      integer, intent(in)     :: ico
      !------------------------------------------------------------------------------------!

      ed_biomass = cpatch%bleaf(ico) + cpatch%bsapwooda(ico)                               &
                 + ( cpatch%bdead(ico) + cpatch%bbark(ico) ) * agf_bs(cpatch%pft(ico))

      return
   end function ed_biomass
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !      This function computes the biomass of everything but storage ("bevery") given    !
   ! size (dbh and height) and the PFT.  This assumes that cohort is in perfect allometry. !
   !                                                                                       !
   !---------------------------------------------------------------------------------------!
   real function size2be(dbh,hite,ipft)
      use pft_coms, only : q     & ! intent(in)
                         , qsw   & ! intent(in)
                         , qbark ! ! intent(in)
      implicit none

      !----- Arguments --------------------------------------------------------------------!
      real          , intent(in) :: dbh
      real          , intent(in) :: hite
      integer       , intent(in) :: ipft
      !----- Local variables --------------------------------------------------------------!
      real                       :: bleaf
      real                       :: bdead
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Find potential leaf and heartwood biomass.                                     !
      !------------------------------------------------------------------------------------!
      bleaf   = size2bl(dbh,hite,ipft)
      bdead   = size2bd(dbh,hite,ipft)
      size2be = bleaf * (1. + q(ipft) + (qsw(ipft)+qbark(ipft)) * hite) + bdead
      !------------------------------------------------------------------------------------!

      return
   end function size2be
   !=======================================================================================!
   !=======================================================================================!




   !=======================================================================================!
   !=======================================================================================!
   !     This function decomposes total biomass (except for storage) into biomass of each  !
   ! tissue, plus the dbh and height.                                                      !
   !---------------------------------------------------------------------------------------!
   subroutine expand_bevery(ipft,bevery,dbh,hite,bleaf,broot,bsapa,bsapb,bbark,balive,bdead)
      use pft_coms    , only : bevery_crit & ! intent(in)
                             , balive_crit & ! intent(in)
                             , agf_bs      & ! intent(in)
                             , q           & ! intent(in)
                             , qsw         & ! intent(in)
                             , qbark       & ! intent(in)
                             , hgt_max     & ! intent(in)
                             , dbh_lut     & ! intent(in)
                             , bevery_lut  & ! intent(in)
                             , balive_lut  & ! intent(in)
                             , bdead_lut   & ! intent(in)
                             , le_mask_lut & ! intent(out)
                             , ge_mask_lut ! ! intent(in)
      use consts_coms , only : lnexp_min   & ! intent(in)
                             , lnexp_max   ! ! intent(in)
      use therm_lib   , only : toler       & ! intent(in)
                             , maxfpo      ! ! intent(in)
      implicit none

      !----- Arguments --------------------------------------------------------------------!
      integer, intent(in)  :: ipft      ! PFT type                            [        ---]
      real   , intent(in)  :: bevery    ! Biomass (Everything but storage)    [  kgC/plant]
      real   , intent(out) :: dbh       ! Diameter at breast height           [         cm]
      real   , intent(out) :: hite      ! Cohort height                       [          m]
      real   , intent(out) :: bleaf     ! Leaf biomass                        [  kgC/plant]
      real   , intent(out) :: broot     ! Root biomass                        [  kgC/plant]
      real   , intent(out) :: bsapa     ! Above-ground sapwood biomass        [  kgC/plant]
      real   , intent(out) :: bsapb     ! Below-ground sapwood biomass        [  kgC/plant]
      real   , intent(out) :: bbark     ! Bark biomass                        [  kgC/plant]
      real   , intent(out) :: balive    ! Live biomass                        [  kgC/plant]
      real   , intent(out) :: bdead     ! Heartwood biomass                   [  kgC/plant]
      !----- Local variables. -------------------------------------------------------------!
      integer              :: ilwr      ! Lower index of the lookup table
      integer              :: iupr      ! Upper index of the lookup table
      integer              :: it        ! Iteration counter
      real                 :: finterp   ! Interpolation factor
      real                 :: salloci   ! Allocation parameters
      real                 :: dbha      ! Lower guess for DBH
      real                 :: dbhz      ! Upper guess for DBH
      real                 :: hgta      ! Lower guess for height
      real                 :: hgtz      ! Upper guess for height
      real                 :: funa      ! Function evaluation for lower guess
      real                 :: funz      ! Function evaluation for upper guess
      real                 :: fun       ! Function evaluation for new   guess
      logical              :: zside     ! Last update was on zside
      logical              :: converged ! Iterative method converged
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Decide which coefficients to use based on the critical bdead.                   !
      !------------------------------------------------------------------------------------!
      if (bevery <= bevery_lut(1,ipft)) then
         !----- Use the look-up table to find the best dbh. -------------------------------!
         finterp = bevery / bevery_lut(1,ipft)
         dbh     = dbh_lut(1,ipft)    * bevery / bevery_lut(1,ipft)
         hite    = dbh2h(ipft,dbh)
         bdead   = bdead_lut(1,ipft)  * bevery / bevery_lut(1,ipft)
         balive  = balive_lut(1,ipft) * bevery / bevery_lut(1,ipft)
         salloci = 1. / (1. + q(ipft) + (qsw(ipft)+qbark(ipft)) * hite )
         bleaf   =                                            salloci * balive
         broot   =                       q    (ipft)        * salloci * balive
         bsapa   =       agf_bs(ipft)  * qsw  (ipft) * hite * salloci * balive
         bsapb   = (1. - agf_bs(ipft)) * qsw  (ipft) * hite * salloci * balive
         bbark   =                       qbark(ipft) * hite * salloci * balive
         !---------------------------------------------------------------------------------!
      else if (bevery >= bevery_crit(ipft)) then
         !---------------------------------------------------------------------------------!
         !     Bdead is above critical value, height is known.                             !
         !---------------------------------------------------------------------------------!
         balive  = balive_crit(ipft)
         bdead   = bevery - balive
         dbh     = bd2dbh(ipft,bdead)
         hite    = hgt_max(ipft)
         salloci = 1. / (1. + q(ipft) + (qsw(ipft)+qbark(ipft)) * hite )
         bleaf   =                                            salloci * balive
         broot   =                       q    (ipft)        * salloci * balive
         bsapa   =       agf_bs(ipft)  * qsw  (ipft) * hite * salloci * balive
         bsapb   = (1. - agf_bs(ipft)) * qsw  (ipft) * hite * salloci * balive
         bbark   =                       qbark(ipft) * hite * salloci * balive

         !---------------------------------------------------------------------------------!
      else
         !----- Use the look-up table to find the best dbh. -------------------------------!
         le_mask_lut(:) = bevery <= bevery_lut(:,ipft)
         ge_mask_lut(:) = bevery >= bevery_lut(:,ipft)
         ilwr           = maxloc (bevery_lut(:,ipft),dim=1,mask=le_mask_lut)
         iupr           = minloc (bevery_lut(:,ipft),dim=1,mask=ge_mask_lut)
         !---------------------------------------------------------------------------------!

         !---------------------------------------------------------------------------------!
         !      In case ilwr and iupr are the same, we have an exact estimate.  Other-     !
         ! wise, use log-linear interpolation.                                             !
         !---------------------------------------------------------------------------------!
         if (ilwr == iupr) then
            dbh     = dbh_lut(ilwr,ipft)
            hite    = dbh2h(ipft,dbh)
         else
            !------ Define the first guess for Regula Falsi (Illinois) method. ------------!
            dbha    = dbh_lut(ilwr,ipft)
            dbhz    = dbh_lut(iupr,ipft)
            hgta    = dbh2h(ipft,dbha)
            hgtz    = dbh2h(ipft,dbhz)
            funa    = size2be(dbha,hgta,ipft) - bevery
            funz    = size2be(dbhz,hgtz,ipft) - bevery
            zside   = funa * funz < 0.
            !------------------------------------------------------------------------------!


            !----- This should not happen, but check it anyway. ---------------------------!
            if (.not. zside) then
               write (unit=*,fmt='(a)')           '---------------------------------------'
               write (unit=*,fmt='(a)')           ' Function does not have proper guess:  '
               write (unit=*,fmt='(a)')           '---------------------------------------'
               write (unit=*,fmt='(a)')           ' '
               write (unit=*,fmt='(a)')           ' Input:    '
               write (unit=*,fmt='(a,1x,i14)'   ) ' + ipft   =',ipft
               write (unit=*,fmt='(a,1x,es14.7)') ' + bevery =',bevery
               write (unit=*,fmt='(a)')           ' '
               write (unit=*,fmt='(a)')           ' Guesses:  '
               write (unit=*,fmt='(a,1x,es14.7)') ' + dbha    =',dbha
               write (unit=*,fmt='(a,1x,es14.7)') ' + dbhz    =',dbhz
               write (unit=*,fmt='(a,1x,es14.7)') ' + hgta    =',hgta
               write (unit=*,fmt='(a,1x,es14.7)') ' + hgtz    =',hgtz
               write (unit=*,fmt='(a,1x,es14.7)') ' + beverya =',size2be(dbha,hgta,ipft)
               write (unit=*,fmt='(a,1x,es14.7)') ' + beveryz =',size2be(dbhz,hgtz,ipft)
               write (unit=*,fmt='(a,1x,es14.7)') ' + funa    =',funa
               write (unit=*,fmt='(a,1x,es14.7)') ' + funz    =',funz
               write (unit=*,fmt='(a)')           ' '
               write (unit=*,fmt='(a)')           '--------------------------------------'
               call fatal_error('Ill-posed first guesses for Regula Falsi'                 &
                               ,'expand_bevery','allometry.f90')
               
            end if
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !     Loop until convergence.                                                  !
            !------------------------------------------------------------------------------!
            rfaloop: do it=1,maxfpo
               dbh  = (funz * dbha - funa * dbhz) / ( funz - funa)
               hite = dbh2h(ipft,dbh)

               !---------------------------------------------------------------------------!
               !     Now that we updated the guess, check whether they are really close.   !
               ! If so, it converged within tolerance.                                     !
               !---------------------------------------------------------------------------!
               converged = abs(dbh-dbha) < toler * dbh
               if (converged) exit rfaloop
               !---------------------------------------------------------------------------!


               !------ Find the new function evaluation. ----------------------------------!
               fun  =  size2be(dbh,hite,ipft) - bevery
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !     Define the new interval based on the intermediate value theorem.  We  !
               ! check whether the same side is being picked repeatedly, and apply the     !
               ! Regula Falsi step to try to speed up convergence.                         !
               !---------------------------------------------------------------------------!
               if (fun*funa < 0. ) then
                  dbhz  = dbh
                  funz  = fun
                  if (zside) funa = funa * 0.5
                  zside = .true.
               else
                  dbha  = dbh
                  funa  = fun
                  if (.not. zside) funz = funz * 0.5
                  zside = .false.
               end if
               !---------------------------------------------------------------------------!
            end do rfaloop
            !------------------------------------------------------------------------------!




            !------------------------------------------------------------------------------!
            !      This method should always congerge in case the function is continuous,  !
            ! but just in case.                                                            !
            !------------------------------------------------------------------------------!
            if (.not. converged) then
               write (unit=*,fmt='(a)')           '---------------------------------------'
               write (unit=*,fmt='(a)')           ' Function did not converge:            '
               write (unit=*,fmt='(a)')           '---------------------------------------'
               write (unit=*,fmt='(a)')           ' '
               write (unit=*,fmt='(a,1x,i14)'   ) ' Iterations =',maxfpo
               write (unit=*,fmt='(a)')           ' '
               write (unit=*,fmt='(a)')           ' Input:      '
               write (unit=*,fmt='(a,1x,i14)'   ) ' + ipft     =',ipft
               write (unit=*,fmt='(a,1x,es14.7)') ' + bevery   =',bevery
               write (unit=*,fmt='(a)')           ' '
               write (unit=*,fmt='(a)')           ' Guesses:    '
               write (unit=*,fmt='(a,1x,es14.7)') ' + dbha      =',dbha
               write (unit=*,fmt='(a,1x,es14.7)') ' + dbh       =',dbh
               write (unit=*,fmt='(a,1x,es14.7)') ' + dbhz      =',dbhz
               write (unit=*,fmt='(a,1x,es14.7)') ' + hgta      =',hgta
               write (unit=*,fmt='(a,1x,es14.7)') ' + hite      =',hite
               write (unit=*,fmt='(a,1x,es14.7)') ' + hgtz      =',hgtz
               write (unit=*,fmt='(a,1x,es14.7)') ' + beverya   =',size2be(dbha,hgta,ipft)
               write (unit=*,fmt='(a,1x,es14.7)') ' + beveryz   =',size2be(dbhz,hgtz,ipft)
               write (unit=*,fmt='(a,1x,es14.7)') ' + funa      =',funa
               write (unit=*,fmt='(a,1x,es14.7)') ' + fun       =',fun
               write (unit=*,fmt='(a,1x,es14.7)') ' + funz      =',funz
               write (unit=*,fmt='(a)')           ' '
               write (unit=*,fmt='(a)')           '--------------------------------------'
               call fatal_error('Ill-posed first guess for Regula Falsi'                   &
                               ,'expand_bevery','allometry.f90')
               
            end if
            !------------------------------------------------------------------------------!
         end if
         !---------------------------------------------------------------------------------!


         !------ Solution for dbh was determined, derive tissue biomass. ------------------!
         bdead   = size2bd(dbh,hite,ipft)
         balive  = bevery - bdead
         salloci = 1. / (1. + q(ipft) + (qsw(ipft)+qbark(ipft)) * hite )
         bleaf   =                                            salloci * balive
         broot   =                       q    (ipft)        * salloci * balive
         bsapa   =       agf_bs(ipft)  * qsw  (ipft) * hite * salloci * balive
         bsapb   = (1. - agf_bs(ipft)) * qsw  (ipft) * hite * salloci * balive
         bbark   =                       qbark(ipft) * hite * salloci * balive
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!

      return
   end subroutine expand_bevery
   !=======================================================================================!
   !=======================================================================================!


   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine estimates the tree area indices, namely leaf, branch(plus twigs), !
   ! and stem.  For the leaf area index (LAI), we use the specific leaf area (SLA), a      !
   ! constant.  The wood area index WAI is found either from Hormann et al. (2003) or      !
   ! based on the average WAI/LAI ratio from Olivas et al. (2013).                         !
   !                                                                                       !
   ! G. Hormann, S. Irrgan, H. Jochheim, M. Lukes, H. Meesenburg, J. Muller, B. Scheler,   !
   !    J. Scherzer, G. Schuler, B. Schultze, B. Strohbach, F. Suckow, M. Wegehenkel, and  !
   !    G. Wessolek.   Wasserhaushalt von waldokosystemen: methodenleitfaden zur           !
   !    bestimmung der wasserhaushaltskomponenten auf level II-flachen. Technical note,    !
   !    Bundesministerium fur Verbraucherschutz, Ernahrung und Landwirtschaft (BMVEL),     !
   !    Bonn, Germany, 2003. URL http://www.wasklim.de/download/Methodenband.pdf.          !
   !                                                                                       !
   ! P. C. Olivas, S. F. Oberbauer, D. B. Clark, D. A. Clark, M. G. Ryan, J. J. O'Brien,   !
   !    and H. Ordonez.  Comparison of direct and indirect methods for assessing leaf area !
   !    index across a tropical rain forest landscape.  Agric. For. Meteorol., 177,        !
   !    110--116, 2013.                                                                    !
   !---------------------------------------------------------------------------------------!
   subroutine area_indices(cpatch, ico)
      use ed_state_vars, only : patchtype       ! ! Structure
      use pft_coms     , only : dbh_crit        & ! intent(in)
                              , is_liana        & ! intent(in)
                              , is_tropical     & ! intent(in)
                              , is_grass        & ! intent(in)
                              , b1WAI           & ! intent(in)
                              , b2WAI           & ! intent(in)
                              , liana_dbh_crit  ! ! intent(in)
      use rk4_coms     , only : ibranch_thermo  ! ! intent(in)
      use ed_misc_coms , only : igrass          & ! intent(in)
                              , iallom          ! ! intent(in)
      implicit none

      !----- Arguments --------------------------------------------------------------------!
      type(patchtype), target :: cpatch
      integer, intent(in)     :: ico
      !----- Local variables --------------------------------------------------------------!
      real                    :: loccai
      real                    :: mdbh
      integer                 :: ipft
      !------------------------------------------------------------------------------------!


      ipft   = cpatch%pft(ico)


      !------------------------------------------------------------------------------------!
      !     First, we compute the LAI.  Because SLA may not be constant, we must use the   !
      ! cohort-level variable as opposed to the pft look-up table.  Likewise, we use       !
      ! actual bleaf because leaves may not be fully flushed.                              !
      !------------------------------------------------------------------------------------!
      cpatch%lai(ico) = cpatch%bleaf(ico) * cpatch%nplant(ico) * cpatch%sla(ico)
      !------------------------------------------------------------------------------------!

      !----- Find the crown area. ---------------------------------------------------------!
      loccai = size2ca(cpatch%dbh(ico),cpatch%hite(ico),cpatch%sla(ico),ipft)
      cpatch%crown_area(ico) = min(1.0, cpatch%nplant(ico) * loccai)
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !     Here we check whether we need to compute the branch, stem, and effective       !
      ! branch area indices.  These are only needed when branch thermodynamics is used,    !
      ! otherwise, simply assign zeroes to them.                                           !
      !------------------------------------------------------------------------------------!
      select case (ibranch_thermo)
      case (0) 
         !----- Ignore branches and trunk. ------------------------------------------------!
         cpatch%wai(ico)  = 0.
         !---------------------------------------------------------------------------------!

      case (1,2)
         !---------------------------------------------------------------------------------!
         !     Get the DBH, or potential DBH in case of grasses.                           !
         !---------------------------------------------------------------------------------!
         if (is_grass(ipft) .and. igrass == 1) then 
            !---- Use height for new grasses. ---------------------------------------------!
            mdbh   = min(cpatch%dbh(ico),h2dbh(cpatch%hite(ico),ipft))
         elseif (is_liana(ipft)) then
            mdbh   = min(cpatch%dbh(ico),liana_dbh_crit)
         else
            !---- Use dbh for trees. ------------------------------------------------------!
            mdbh   = min(cpatch%dbh(ico),dbh_crit(ipft))
         end if
         !---------------------------------------------------------------------------------!


         !-----Find WAI. ------------------------------------------------------------------!
         if (iallom == 3 .and. is_tropical(ipft) .and. (.not. is_liana(ipft))) then
            cpatch%wai(ico) = cpatch%nplant(ico) * b1WAI(ipft)                             &
                            * ( mdbh * mdbh * cpatch%hite(ico) ) ** b2WAI(ipft)
         else
            cpatch%wai(ico) = cpatch%nplant(ico) * b1WAI(ipft) * mdbh ** b2WAI(ipft)
         end if
         !---------------------------------------------------------------------------------!
      end select
      !------------------------------------------------------------------------------------!

      return
   end subroutine area_indices
   !=======================================================================================!
   !=======================================================================================!
end module allometry
!==========================================================================================!
!==========================================================================================!


