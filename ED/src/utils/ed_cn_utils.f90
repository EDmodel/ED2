!==========================================================================================!
!==========================================================================================!
! MODULE: ED_CN_UTILS
!
!> \brief   Sub-routines that handle carbon:nitrogen budgets and sthoichiometry.
!> \details This module contains some sub-routines that relate carbon and nitrogen and 
!>          assess their budgets.  These routines used to be in structural growth, however
!>          they were migrated to a separate module to avoid circular dependency.
!> \author  David Medvigy and Michael Dietze.
!> \author  20 Oct 2017.  Migrated routines from structural growth, Marcos Longo
!------------------------------------------------------------------------------------------!
module ed_cn_utils
   contains
   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine will display the carbon and nitrogen budgets on screen.            !
   !---------------------------------------------------------------------------------------!
   subroutine print_C_and_N_budgets(cgrid)
      use ed_state_vars , only : edtype       ! ! structure
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(edtype)                 , target    :: cgrid
      !----- Local variables --------------------------------------------------------------!
      integer                                  :: ipy
      real                                     :: soil_C
      real                                     :: soil_N
      real                                     :: veg_C
      real                                     :: veg_N
      !----- Local constants --------------------------------------------------------------!
      logical                      , parameter :: print_on = .false.
      !------------------------------------------------------------------------------------!



      do ipy = 1,cgrid%npolygons

         call compute_C_and_N_storage(cgrid,ipy, soil_C, soil_N, veg_C, veg_N)

         if (print_on) then

            write (unit=*,fmt='(a)') '-----------------------------------------------------'
            write (unit=*,fmt='(a,1x,i6)')     'POLYGON           =',ipy
            write (unit=*,fmt='(a,1x,f9.3)')   'LON               =',cgrid%lon(ipy)
            write (unit=*,fmt='(a,1x,f9.3)')   'LAT               =',cgrid%lat(ipy)
            write (unit=*,fmt='(a,1x,es12.5)') 'C_INITIAL_STORAGE ='                       &
                                                       ,cgrid%cbudget_initialstorage(ipy)
            write (unit=*,fmt='(a,1x,es12.5)') 'SOIL_C+VEG_C-NCEP ='                       &
                                                       ,soil_C+veg_C-cgrid%cbudget_nep(ipy)
            write (unit=*,fmt='(a,1x,es12.5)') 'SOIL_C+VEG_C      = ',soil_C + veg_C
            write (unit=*,fmt='(a,1x,es12.5)') 'NEP               = ',cgrid%cbudget_nep(ipy)
            write (unit=*,fmt='(a,1x,es12.5)') 'C_REMOVED_STORAGE = '                      &
                                                         ,cgrid%cbudget_removedstorage(ipy)
            write (unit=*,fmt='(a,1x,es12.5)') 'N_INITIAL_STORAGE ='                       &
                                                       ,cgrid%nbudget_initialstorage(ipy)
            write (unit=*,fmt='(a,1x,es12.5)') 'SOIL_N+VEG_N      = ',soil_N + veg_N
            write (unit=*,fmt='(a)') '-----------------------------------------------------'
            write (unit=*,fmt='(a)') ' '
         end if
      end do
      return
   end subroutine print_C_and_N_budgets
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine will compute the carbon and nitrogen pools.                        !
   !---------------------------------------------------------------------------------------!
   subroutine compute_C_and_N_storage(cgrid,ipy, soil_C, soil_N, veg_C, veg_N)

      use ed_state_vars , only : edtype         & ! structure
                               , polygontype    & ! structure
                               , sitetype       & ! structure
                               , patchtype      ! ! structure
      use ed_max_dims   , only : n_pft          ! ! intent(in)
      use pft_coms      , only : include_pft    & ! intent(in)
                               , c2n_recruit    & ! intent(in)
                               , c2n_stem       & ! intent(in)
                               , c2n_leaf       & ! intent(in)
                               , c2n_storage    & ! intent(in)
                               , f_labile_leaf  & ! intent(in)
                               , f_labile_stem  ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(edtype)      , target      :: cgrid
      integer           , intent(in)  :: ipy
      real              , intent(out) :: soil_C
      real              , intent(out) :: soil_N
      real              , intent(out) :: veg_C
      real              , intent(out) :: veg_N
      !----- Local variables --------------------------------------------------------------!
      type(polygontype) , pointer     :: cpoly
      type(sitetype)    , pointer     :: csite
      type(patchtype)   , pointer     :: cpatch
      integer                         :: isi
      integer                         :: ipa
      integer                         :: ico
      integer                         :: ipft
      real                            :: bfast_C
      real                            :: bstruct_C
      real(kind=8)                    :: area_factor
      real(kind=8)                    :: this_carbon
      real(kind=8)                    :: this_nitrogen
      real(kind=8)                    :: soil_C8
      real(kind=8)                    :: soil_N8
      real(kind=8)                    :: veg_C8
      real(kind=8)                    :: veg_N8
      !----- Local constants --------------------------------------------------------------!
      real(kind=8)      , parameter   :: almostnothing=1.d-30
      !----- External functions. ----------------------------------------------------------!
      real              , external    :: sngloff
      !------------------------------------------------------------------------------------!

      !----- Initialize C and N pools. ----------------------------------------------------!
      soil_C8 = 0.0d0
      soil_N8 = 0.0d0
      veg_C8  = 0.0d0
      veg_N8  = 0.0d0

      cpoly => cgrid%polygon(ipy)
      siteloop: do isi = 1,cpoly%nsites

         csite => cpoly%site(isi)
         patchloop: do ipa = 1,csite%npatches
            cpatch => csite%patch(ipa)

            !----- Site area times patch area. --------------------------------------------!
            area_factor   = dble(cpoly%area(isi)) * dble(csite%area(ipa))

            !----- Find carbon and nitrogen soil pools for this patch. --------------------!
            this_carbon   = dble(csite%fast_grnd_C       (ipa))                            &
                          + dble(csite%fast_soil_C       (ipa))                            &
                          + dble(csite%structural_grnd_C (ipa))                            &
                          + dble(csite%structural_soil_C (ipa))                            &
                          + dble(csite%microbial_soil_C  (ipa))                            &
                          + dble(csite%slow_soil_C       (ipa))                            &
                          + dble(csite%passive_soil_C    (ipa))
            this_nitrogen = dble(csite%fast_grnd_N       (ipa))                            &
                          + dble(csite%fast_soil_N       (ipa))                            &
                          + dble(csite%structural_grnd_N (ipa))                            &
                          + dble(csite%structural_grnd_N (ipa))                            &
                          + dble(csite%mineralized_soil_N(ipa))

            !----- Add to the full counter. -----------------------------------------------!
            soil_C8 = soil_C8 + area_factor * this_carbon
            soil_N8 = soil_N8 + area_factor * this_nitrogen

            !----- Loop over PFT so we account for veg carbon/nitrogen in repro arrays. ---!
            pftloop: do ipft = 1, n_pft
               if (include_pft(ipft)) then
                  veg_C8 = veg_C8 + dble(csite%repro(ipft,ipa)) * area_factor
                  veg_N8 = veg_N8 + dble(csite%repro(ipft,ipa))                            &
                                  / dble(c2n_recruit(ipft)) * area_factor
               end if
            end do pftloop

            cohortloop: do ico = 1,cpatch%ncohorts
               ipft = cpatch%pft(ico)

               !----- Get the carbon and nitrogen in vegetation. --------------------------!
               veg_C8 = veg_C8 + area_factor * dble(cpatch%nplant(ico))                    &
                               * ( dble(cpatch%balive  (ico))                              &
                                 + dble(cpatch%bdeada  (ico))                              &
                                 + dble(cpatch%bdeadb  (ico))                              &
                                 + dble(cpatch%bstorage(ico)) ) 
                               
               bfast_C   = f_labile_leaf(ipft)                                             &
                         * ( cpatch%bleaf(ico) + cpatch%broot(ico))                        &
                         + f_labile_stem(ipft)                                             &
                         * ( cpatch%bsapwooda(ico)+cpatch%bbarka(ico)+cpatch%bdeada(ico)   &
                           + cpatch%bsapwoodb(ico)+cpatch%bbarkb(ico)+cpatch%bdeadb(ico))
               bstruct_C = (1.0 - f_labile_leaf(ipft))                                     &
                         * ( cpatch%bleaf(ico) + cpatch%broot(ico))                        &
                         + (1.0 - f_labile_stem(ipft))                                     &
                         * ( cpatch%bsapwooda(ico)+cpatch%bbarka(ico)+cpatch%bdeada(ico)   &
                           + cpatch%bsapwoodb(ico)+cpatch%bbarkb(ico)+cpatch%bdeadb(ico))

               veg_N8 = veg_N8 + area_factor                                               &
                      * ( dble(bfast_C  ) / dble(c2n_leaf(cpatch%pft(ico)))                &
                        + dble(bstruct_C) / dble(c2n_stem(cpatch%pft(ico)))                &
                        + dble(cpatch%bstorage (ico)) / dble(c2n_storage))                 &
                      * dble(cpatch%nplant(ico))
            end do cohortloop
         end do patchloop
      end do siteloop

      soil_C = sngloff(soil_C8,almostnothing)
      soil_N = sngloff(soil_N8,almostnothing)
      veg_C  = sngloff(veg_C8 ,almostnothing)
      veg_N  = sngloff(veg_N8 ,almostnothing)

      return
   end subroutine compute_C_and_N_storage
   !=======================================================================================!
   !=======================================================================================!
end module ed_cn_utils
!==========================================================================================!
!==========================================================================================!
