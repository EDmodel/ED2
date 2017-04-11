module patch_pft_size_profile_mod
  contains

   !=======================================================================================!
   !=======================================================================================!
   subroutine patch_pft_size_profile(csite,ipa)
      use ed_state_vars       , only : sitetype   & ! structure
                                     , patchtype  ! ! structure
      use fusion_fission_coms , only : ff_nhgt    & ! intent(in)
                                     , hgt_class  ! ! intent(in)
      use allometry           , only : size2bl    ! ! intent(in)
      use ed_max_dims         , only : n_pft      ! ! intent(in)
      use pft_coms            , only : hgt_min    & ! intent(in)
                                     , is_grass   ! ! intent(in)
      use ed_misc_coms        , only : igrass     ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(sitetype)         , target     :: csite     ! Current site
      integer                , intent(in) :: ipa       ! Current patch index
      !----- Local variables --------------------------------------------------------------!
      type(patchtype)        , pointer    :: cpatch    ! Current patch
      integer                             :: ipft      ! PFT index
      integer                             :: ihgt      ! Height class index
      integer                             :: ico       ! Counters
      real                                :: lai_pot   ! Potential LAI
      !------------------------------------------------------------------------------------!


      !----- Reset all bins to zero. ------------------------------------------------------!
      do ipft=1,n_pft
         do ihgt=1,ff_nhgt
            csite%cumlai_profile(ipft,ihgt,ipa)=0.0
         end do
      end do
      !------------------------------------------------------------------------------------!



      !----- Update bins ------------------------------------------------------------------!
      cpatch => csite%patch(ipa)
      cohortloop: do ico = 1,cpatch%ncohorts

         !----- Find the PFT class. -------------------------------------------------------!
         ipft    = cpatch%pft(ico)
         ihgt    = min(ff_nhgt,max(1,count(hgt_class < cpatch%hite(ico))))
         !---------------------------------------------------------------------------------!

         !---------------------------------------------------------------------------------!
         !     Check whether this cohort is almost at the minimum height given its PFT.    !
         ! If it is, then we will skip it.                                                 !
         !---------------------------------------------------------------------------------!
         if (cpatch%hite(ico) < hgt_min(ipft) + 0.2) cycle cohortloop
         !---------------------------------------------------------------------------------!


         !----- Find the height class. ----------------------------------------------------!
         ihgt    = min(ff_nhgt,max(1,count(hgt_class < cpatch%hite(ico))))
         !---------------------------------------------------------------------------------!


         !----- Find the potential (on-allometry) leaf area index. ------------------------!
         if (is_grass(ipft) .and. igrass==1) then
             !--use actual bleaf for grass
             lai_pot = cpatch%nplant(ico) * cpatch%sla(ico) * cpatch%bleaf(ico)
         else
             !--use dbh for trees
             lai_pot = cpatch%nplant(ico) * cpatch%sla(ico)                                &
                     * size2bl(cpatch%dbh(ico),cpatch%hite(ico),ipft)
         end if
         !---------------------------------------------------------------------------------!


         !----- Add the potential LAI to the bin. -----------------------------------------!
         csite%cumlai_profile(ipft,ihgt,ipa) = lai_pot                                     &
                                             + csite%cumlai_profile(ipft,ihgt,ipa)
         !---------------------------------------------------------------------------------!
      end do cohortloop
      !------------------------------------------------------------------------------------!



      !----- Integrate the leaf area index from top to bottom. ----------------------------!
      do ihgt=ff_nhgt-1,1,-1
         do ipft=1,n_pft
            csite%cumlai_profile(ipft,ihgt,ipa) = csite%cumlai_profile(ipft,ihgt  ,ipa)    &
                                                + csite%cumlai_profile(ipft,ihgt+1,ipa)
         end do
      end do
      !------------------------------------------------------------------------------------!

      return
   end subroutine patch_pft_size_profile
   !=======================================================================================!
   !=======================================================================================!


end module patch_pft_size_profile_mod