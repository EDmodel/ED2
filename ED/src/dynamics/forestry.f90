!==========================================================================================!
! Forestry.f90. These subroutines will control the creation of new patches based on        !
!               anthropogenic disturbances such as conversion to agriculture, land         !
!               abandonment and harvesting.                                                !
!                                                                                          !
! NOTICE:  These subroutines have not been thoroughly tested. Please                       !
!          report problems to David Medvigy, medvigy@post.harvard.edu.                     !
!==========================================================================================!
subroutine apply_forestry(cpoly, isi, year)
   use ed_state_vars        , only : polygontype                & ! structure
                                   , sitetype                   & ! structure
                                   , patchtype                  & ! structure
                                   , allocate_sitetype          & ! subroutine
                                   , deallocate_sitetype        & ! subroutine
                                   , copy_sitetype_mask         ! ! subroutine
   use disturb_coms         , only : ianth_disturb              & ! intent(in)
                                   , lutime                     & ! intent(in)
                                   , min_new_patch_area         & ! intent(in)
                                   , plantation_year            & ! intent(in)
                                   , plantation_rotation        & ! intent(in)
                                   , mature_harvest_age         ! ! intent(in)
   use disturbance_utils    , only : initialize_disturbed_patch & ! subroutine
                                   , plant_patch                ! ! subroutine
   use fuse_fiss_utils      , only : terminate_patches          ! ! subroutine
   use ed_max_dims          , only : n_pft                      & ! intent(in)
                                   , n_dbh                      ! ! intent(in)
   use grid_coms            , only : nzg                        & ! intent(in)
                                   , nzs                        ! ! intent(in)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(polygontype)             , target      :: cpoly
   integer                       , intent(in)  :: year
   integer                       , intent(in)  :: isi
   !----- Local variables -----------------------------------------------------------------!
   type(sitetype)                , pointer     :: csite
   type(patchtype)               , pointer     :: cpatch
   type(sitetype)                , pointer     :: tempsite
   type(lutime)                  , pointer     :: clutime
   logical        ,  dimension(:), allocatable :: mask
   logical                                     :: mature_plantation
   logical                                     :: mature_primary
   logical                                     :: mature_secondary
   integer                                     :: ipft
   integer                                     :: idbh
   integer                                     :: newp
   integer                                     :: iyear
   integer                                     :: useyear
   integer                                     :: ipa
   integer                                     :: ico
   real                                        :: primary_harvest_target
   real                                        :: secondary_harvest_target
   real                                        :: total_harvest_target
   real                                        :: total_site_biomass
   real                                        :: area_mature_primary
   real                                        :: agb_mature_primary
   real                                        :: area_mature_secondary
   real                                        :: agb_mature_secondary
   real                                        :: area_mature_plantation
   real                                        :: agb_mature_plantation
   real                                        :: total_harvested_area
   real                                        :: lambda_mature_primary
   real                                        :: lambda_mature_secondary
   real                                        :: lambda_mature_plantation
   real                                        :: harvest_deficit
   !---------------------------------------------------------------------------------------!

   !----- This subroutine will be skipped if anthropogenic disturbance is turned off. -----!
   if (ianth_disturb == 0) return

   csite => cpoly%site(isi)

   !---------------------------------------------------------------------------------------!
   !    The site is harvested up to a target biomass.  Patches are harvested first from    !
   ! those above the harvest_age, with equal rates.  If the biomass target is not met, the !
   ! remaining biomass is harvested from the patches below the minimum age, starting with  !
   ! the oldest.  Harvest rates are taken from George Hurtt's GLU, Global landuse files.   !
   ! Elements 12 and 16 are secondary harvesting and 14 and 18 are primary.                !
   !---------------------------------------------------------------------------------------!


   !---- First, find the right year in the clutimes vector. -------------------------------!
   useyear = cpoly%num_landuse_years(isi)

   if (year >= cpoly%clutimes(1,isi)%landuse_year) then
      !----- Loop over years. -------------------------------------------------------------!
      find_lu_year: do iyear = 1,cpoly%num_landuse_years(isi)
         if (year == cpoly%clutimes(iyear,isi)%landuse_year) then
            useyear = iyear
            exit find_lu_year
         end if
      end do find_lu_year
   else
      call fatal_error('Invalid initial year when using land-use disturbance'              &
                      &,'apply_forestry','forestry.f90')
   end if

   clutime => cpoly%clutimes(useyear,isi)
   
   !---------------------------------------------------------------------------------------!
   !      Set primary target based on current rates and unapplied harvest from previous    !
   ! years (memory).  These are in kgC/m2.  At this point we must check whether a PFT-     !
   ! blind or a selective logging must be applied, and if it is selective logging, we do   !
   ! not log anything here.                                                                !
   !---------------------------------------------------------------------------------------!
   if (clutime%landuse(14) > 0.) then
      primary_harvest_target   = clutime%landuse(14) + clutime%landuse(18)                 &
                               + cpoly%primary_harvest_memory(isi)
   else
      primary_harvest_target   = 0.
   end if
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Set primary and secondary targets based on current rates and unapplied harvest   !
   ! from previous years (memory).  These are in kgC/m2.  At this point we must check      !
   ! whether a PFT-blind or a selective logging must be applied, and if it is selective    !
   ! logging, we do not log anything here.                                                 !
   !---------------------------------------------------------------------------------------!
   if (clutime%landuse(12) > 0.) then
      secondary_harvest_target   = clutime%landuse(12) + clutime%landuse(16)               &
                                 + cpoly%primary_harvest_memory(isi)
   else
      secondary_harvest_target   = 0.
   end if
   !---------------------------------------------------------------------------------------!



   !----- The total target is the sum of both targets. ------------------------------------!
   total_harvest_target     = primary_harvest_target + secondary_harvest_target
   !---------------------------------------------------------------------------------------!



   !=======================================================================================!
   !=======================================================================================!
   !     Decide whether or not to create a harvest patch. A new patch must meet the        !
   ! following criteria:                                                                   !
   ! a. It must have site agb > 0;                                                         !
   ! b. Harvest must exceed some minimum threshold.                                        !
   !    Note: this is not necessarily the total harvest area.                              !
   !---------------------------------------------------------------------------------------!
   

   !---------------------------------------------------------------------------------------!
   !    Finding total biomass density in kgC/m2.                                           !
   !---------------------------------------------------------------------------------------!
   total_site_biomass = 0.
   pftagbloop: do ipft=1,n_pft
      pftdbhloop: do idbh=1,n_dbh
         total_site_biomass = total_site_biomass + cpoly%agb(ipft, idbh, isi)
      end do pftdbhloop
   end do pftagbloop
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     If site has no biomass, or if the area to be harvested is less than the minimum   !
   ! area for a new patch, do not harvest, and update the memory for the next year.        !
   !---------------------------------------------------------------------------------------!
   if (total_site_biomass == 0.0 .or.                                                      &
       total_harvest_target <= total_site_biomass * min_new_patch_area) then
      cpoly%primary_harvest_memory(isi)   = primary_harvest_target
      cpoly%secondary_harvest_memory(isi) = secondary_harvest_target
      return
   end if
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    Conditions were met, making a new patch.  The following routines extends the       !
   ! allocation of the current site by one patch, and preserves the original patches       !
   !---------------------------------------------------------------------------------------!
   nullify(tempsite)
   allocate(tempsite)
   allocate(mask(csite%npatches))
   call allocate_sitetype(tempsite,csite%npatches)
   mask(:) = .true.
   call copy_sitetype_mask(csite,tempsite,mask,csite%npatches,csite%npatches)
   call deallocate_sitetype(csite)
   call allocate_sitetype(csite,tempsite%npatches + 1)
   call copy_sitetype_mask(tempsite,csite,mask,tempsite%npatches,tempsite%npatches)
   call deallocate_sitetype(tempsite)
   deallocate(tempsite)
   deallocate(mask)
   !---------------------------------------------------------------------------------------!


   
   newp = csite%npatches

   !------ Initialize the new patch (newp) in the last position. --------------------------!
   csite%dist_type(newp) = 2
   ! Area is not initialized.
   call initialize_disturbed_patch(csite,cpoly%met(isi)%atm_tmp,newp,1,cpoly%lsl(isi))

   !------ Compute current stocks of agb in mature forests. -------------------------------!
   call inventory_mat_forests(cpoly,isi,area_mature_primary,agb_mature_primary             &
                             ,area_mature_secondary,agb_mature_secondary                   &
                             ,area_mature_plantation,agb_mature_plantation)     

   !------ Compute the mature-forest harvest rates. ---------------------------------------!
   call mat_forest_harv_rates(agb_mature_primary,agb_mature_secondary                      &
                             ,agb_mature_plantation,primary_harvest_target                 &
                             ,secondary_harvest_target,lambda_mature_primary               &
                             ,lambda_mature_secondary,lambda_mature_plantation             &
                             ,harvest_deficit)                                    

   !------ Apply harvesting to the mature stands. -----------------------------------------!
   call harv_mat_patches(cpoly,isi,newp,lambda_mature_primary                              &
                        ,lambda_mature_secondary,lambda_mature_plantation)
   
   !---------------------------------------------------------------------------------------!
   !     Compute harvested area from mature patches.  This is also updated in              !
   ! harvest_immature_patches().                                                           !
   !---------------------------------------------------------------------------------------!
   total_harvested_area = lambda_mature_primary    * area_mature_primary                   &
                        + lambda_mature_secondary  * area_mature_secondary                 &
                        + lambda_mature_plantation * area_mature_plantation
   call harv_immat_patches(cpoly,isi,newp,harvest_deficit,total_harvested_area)

   !---------------------------------------------------------------------------------------!
   !     Now we know the area of the new patch, and can normalize the averaged patch       !
   ! quantities. But this is done only when the new patch area is significant otherwise,   ! 
   ! just terminate it.                                                                    !
   !---------------------------------------------------------------------------------------!
   csite%area(newp) = total_harvested_area
   if (total_harvested_area > min_new_patch_area) then
      write(unit=*,fmt='(a,1x,i5)')     'LANDUSE YEAR          =',clutime%landuse_year
      write(unit=*,fmt='(a,1x,es12.5)') 'LANDUSE 14            =',clutime%landuse(14)
      write(unit=*,fmt='(a,1x,es12.5)') 'LANDUSE 18            =',clutime%landuse(18)
      write(unit=*,fmt='(a,1x,es12.5)') 'PRIMARY_HARVEST_MEM   ='                          &
                                                         ,cpoly%primary_harvest_memory(isi)
      write(unit=*,fmt='(a,1x,es12.5)') 'LAMBDA_MATURE_PRIMARY =',lambda_mature_primary
      write(unit=*,fmt='(a,1x,es12.5)') 'LAMBDA_MATURE_2NDARY  =',lambda_mature_secondary
      write(unit=*,fmt='(a,1x,es12.5)') 'LAMBDA_MATURE_PLANT   =',lambda_mature_plantation
      write(unit=*,fmt='(a,1x,es12.5)') 'AREA_MATURE_PRIMARY   =',area_mature_primary
      write(unit=*,fmt='(a,1x,es12.5)') 'AREA_MATURE_2NDARY    =',area_mature_secondary
      write(unit=*,fmt='(a,1x,es12.5)') 'AREA_MATURE_PLANT     =',area_mature_plantation
      write(unit=*,fmt='(a,1x,es12.5,1x,a,1x,i5)')                                         &
          ' ---> Making new patch (harvesting), with area=',csite%area(newp)               &
              ,' for dist_type=',2
      call norm_harv_patch(csite,newp)
      
      !----- Update temperature and density. ----------------------------------------------!
      call update_patch_thermo_props(csite,newp,newp,nzg,nzs,cpoly%ntext_soil(:,isi))

      !----- Plant the patch if it is a plantation. ---------------------------------------!
      if (cpoly%plantation(isi) == 1 .and. year > plantation_year) then
         call plant_patch(csite,newp,nzg,cpoly%plantation_stocking_pft(isi)                &
                         ,cpoly%plantation_stocking_density(isi)                           &
                         ,cpoly%ntext_soil(:,isi),cpoly%green_leaf_factor(:,isi),2.0       &
                         ,cpoly%lsl(isi))
         csite%plantation(newp) = 1
      end if
      call update_patch_derived_props(csite,cpoly%lsl(isi),cpoly%met(isi)%prss,newp)
      call new_patch_sfc_props(csite,newp,nzg,nzs,cpoly%ntext_soil(:,isi))
      call update_budget(csite,cpoly%lsl(isi),newp,newp)
   end if

   !----- Eliminate those patches with small area. ----------------------------------------!
   call terminate_patches(csite)

   !----- Clear out the primary harvest memory. -------------------------------------------!
   cpoly%primary_harvest_memory(isi) = 0.0

   !----- There still may be a deficit if we have harvested all of the patch agb. ---------!
   cpoly%secondary_harvest_memory(isi) = harvest_deficit
   
   return
end subroutine apply_forestry
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine inventory_mat_forests(cpoly,isi,area_mature_primary,agb_mature_primary          &
                                   ,area_mature_secondary, agb_mature_secondary            &
                                   ,area_mature_plantation, agb_mature_plantation)
   use ed_state_vars , only : polygontype         & ! structure
                            , sitetype            & ! structure
                            , patchtype           ! ! structure
   use disturb_coms  , only : plantation_rotation & ! intent(in)
                            , mature_harvest_age  ! ! intent(in)
   use allometry     , only : ed_biomass          ! ! function
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(polygontype) , target      :: cpoly
   integer           , intent(in)  :: isi
   real              , intent(out) :: area_mature_primary
   real              , intent(out) :: agb_mature_primary
   real              , intent(out) :: area_mature_secondary
   real              , intent(out) :: agb_mature_secondary
   real              , intent(out) :: area_mature_plantation
   real              , intent(out) :: agb_mature_plantation
   !----- Local variables -----------------------------------------------------------------!
   type(sitetype)    , pointer     :: csite
   type(patchtype)   , pointer     :: cpatch
   logical                         :: mature_plantation
   logical                         :: mature_primary
   logical                         :: mature_secondary
   integer                         :: ipa
   integer                         :: ico
   integer                         :: ipft
   !---------------------------------------------------------------------------------------!



   !----- Initialize inventory. -----------------------------------------------------------!
   area_mature_primary    = 0.0
   area_mature_secondary  = 0.0
   area_mature_plantation = 0.0
   agb_mature_primary     = 0.0
   agb_mature_secondary   = 0.0
   agb_mature_plantation  = 0.0

   csite => cpoly%site(isi)
   patchloop: do ipa=1,csite%npatches

      cpatch => csite%patch(ipa)

      !----- Compute the patch AGB --------------------------------------------------------!
      csite%plant_ag_biomass(ipa) = 0.0
      do ico=1,cpatch%ncohorts
         csite%plant_ag_biomass(ipa) = csite%plant_ag_biomass(ipa)                         &
                                     + cpatch%agb(ico) * cpatch%nplant(ico)
      end do

      !----- Skip the patch if the biomass is low. ----------------------------------------!
      if (csite%plant_ag_biomass(ipa) < 0.01) cycle patchloop


      !----- Run some checks to determine if this patch could be possibly harvested. ------!
      mature_primary    = csite%dist_type(ipa)  == 3                   .and.               &
                          csite%age(ipa)         > mature_harvest_age
      mature_plantation = csite%plantation(ipa) == 1                   .and.               &
                          csite%age(ipa)         > plantation_rotation
      mature_secondary  = csite%dist_type(ipa)  == 2                   .and.               &
                          csite%plantation(ipa) /= 1                   .and.               &
                          csite%age(ipa)         > mature_harvest_age


      !------------------------------------------------------------------------------------!
      !    Add the biomass and area to the appropriate type of vegetation.                 !
      !------------------------------------------------------------------------------------!
      if (mature_plantation) then
         area_mature_plantation = area_mature_plantation      + csite%area(ipa)
         agb_mature_plantation  = agb_mature_plantation                                    &
                                + csite%plant_ag_biomass(ipa) * csite%area(ipa)

      elseif (mature_secondary) then
         area_mature_secondary = area_mature_secondary + csite%area(ipa)
         agb_mature_secondary  = agb_mature_secondary                                      &
                               + csite%plant_ag_biomass(ipa) * csite%area(ipa)

      elseif (mature_primary) then
         area_mature_primary = area_mature_primary + csite%area(ipa)
         agb_mature_primary  = agb_mature_primary                                          &
                             + csite%plant_ag_biomass(ipa) * csite%area(ipa)

      end if

   end do patchloop
   return
end subroutine inventory_mat_forests
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine mat_forest_harv_rates(agb_mature_primary,agb_mature_secondary                   &
                                ,agb_mature_plantation, primary_harvest_target             &
                                ,secondary_harvest_target,lambda_mature_primary            &
                                ,lambda_mature_secondary,lambda_mature_plantation          &
                                ,harvest_deficit)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   real, intent(in)    :: agb_mature_primary
   real, intent(in)    :: agb_mature_secondary
   real, intent(in)    :: agb_mature_plantation
   real, intent(in)    :: primary_harvest_target
   real, intent(inout) :: secondary_harvest_target
   real, intent(out)   :: lambda_mature_primary
   real, intent(out)   :: lambda_mature_plantation
   real, intent(out)   :: lambda_mature_secondary
   real, intent(out)   :: harvest_deficit
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !    Compute harvesting rate in mature primary forest.  If there is not enough biomass  !
   ! to harvest, harvest what is possible from primary and attempt to harvest the          !
   ! remainder of the target from secondary.                                               !
   !---------------------------------------------------------------------------------------!
   if (agb_mature_primary > primary_harvest_target) then
      lambda_mature_primary = primary_harvest_target / agb_mature_primary
   else
      lambda_mature_primary    = 1.0
      harvest_deficit          = primary_harvest_target   - agb_mature_primary
      secondary_harvest_target = secondary_harvest_target + harvest_deficit
   end if

   !---------------------------------------------------------------------------------------!
   !      Compute harvesting rate in mature plantations and mature secondary forests.      !
   ! First try to remove all biomass from plantations.  If this is not possible, remove    !
   ! what you could and try to remove the remainder from mature secondary forests.  If     !
   ! this is again not possible, store what is left in harvest_deficit.                    !
   !---------------------------------------------------------------------------------------!
   if (agb_mature_plantation > secondary_harvest_target) then
      lambda_mature_plantation = secondary_harvest_target / agb_mature_plantation
      lambda_mature_secondary  = 0.0
      harvest_deficit          = 0.0
   else
      lambda_mature_plantation = 1.0
      harvest_deficit          = secondary_harvest_target - agb_mature_plantation

      if(agb_mature_secondary > harvest_deficit)then
         lambda_mature_secondary = harvest_deficit / agb_mature_secondary
         harvest_deficit         = 0.0
      else
         lambda_mature_secondary = 1.0
         harvest_deficit         = harvest_deficit - agb_mature_secondary
      end if
   endif

   return
end subroutine mat_forest_harv_rates
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine harv_mat_patches(cpoly,isi,newp,lambda_mature_primary                           &
                           ,lambda_mature_secondary,lambda_mature_plantation)
   use ed_state_vars    , only : polygontype          & ! structure
                               , sitetype             & ! structure
                               , patchtype            ! ! structure
   use disturb_coms     , only : mature_harvest_age   & ! intent(in)
                               , plantation_rotation  ! ! intent(in)
   use disturbance_utils, only : accum_dist_litt      & ! subroutine
                               , insert_survivors     & ! subroutine
                               , increment_patch_vars ! ! subroutine
   use ed_max_dims      , only : n_pft                ! ! intent(in)

   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(polygontype)                  , target     :: cpoly
   integer                            , intent(in) :: isi
   integer                            , intent(in) :: newp
   real                               , intent(in) :: lambda_mature_plantation
   real                               , intent(in) :: lambda_mature_secondary
   real                               , intent(in) :: lambda_mature_primary
   !----- Local variables -----------------------------------------------------------------!
   type(sitetype)                     , pointer    :: csite
   type(patchtype)                    , pointer    :: cpatch
   real             , dimension(n_pft)             :: mindbh_harvest
   integer                                         :: ipa
   logical                                         :: mature_plantation
   logical                                         :: mature_primary
   logical                                         :: mature_secondary
   real                                            :: dA
   !----- Local constants. ----------------------------------------------------------------!
   integer                            , parameter  :: new_lu         = 2
   integer                            , parameter  :: poly_dist_type = 1
   !---------------------------------------------------------------------------------------!



   !----- Loop over patches. --------------------------------------------------------------!
   csite => cpoly%site(isi)

   do ipa=1,csite%npatches

      cpatch => csite%patch(ipa)

      !----- Run some checks to determine if this patch could be possibly harvested. ------!
      mature_primary    = csite%dist_type(ipa)  == 3                   .and.               &
                          csite%age(ipa)         > mature_harvest_age
      mature_plantation = csite%plantation(ipa) == 1                   .and.               &
                          csite%age(ipa)         > plantation_rotation
      mature_secondary  = csite%dist_type(ipa)  == 2                   .and.               &
                          csite%plantation(ipa) /= 1                   .and.               &
                          csite%age(ipa)         > mature_harvest_age

      if (mature_plantation) then
         !----- Harvest mature plantations. -----------------------------------------------!
         dA                      = csite%area(ipa) * lambda_mature_plantation
         mindbh_harvest(1:n_pft) = cpoly%mindbh_secondary(1:n_pft,isi)
      elseif(mature_secondary)then
         !----- Harvest mature secondary. -------------------------------------------------!
         dA                      = csite%area(ipa) * lambda_mature_secondary
         mindbh_harvest(1:n_pft) = cpoly%mindbh_secondary(1:n_pft,isi)
      elseif (mature_primary) then
         !----- Harvest mature primary. ---------------------------------------------------!
         dA                      = csite%area(ipa) * lambda_mature_primary
         mindbh_harvest(1:n_pft) = cpoly%mindbh_primary(1:n_pft,isi)
      else
         !----- Immature patches not harvested here. --------------------------------------!
         dA                      = 0.0  
         mindbh_harvest(1:n_pft) = huge(1.)
      end if


     
      !------ Found a patch that is contributing to the new patch. ------------------------!
      if (dA > 0.0 .and. csite%plant_ag_biomass(ipa) >= 0.01) then
         csite%area(ipa) = csite%area(ipa) - dA
         call increment_patch_vars(csite,newp,ipa,dA)
         !---------------------------------------------------------------------------------!
         !     The destination patch disturbance type was previously set to 1 (agri-       !
         ! culture), but I think it should be 2 (secondary forest) - MLO.  Added the       !
         ! insert survivors subroutine here just to generalise, with the target biomass    !
         ! the survivorship should be 0.                                                   !
         !---------------------------------------------------------------------------------!
         call accum_dist_litt(csite,newp,ipa,new_lu,dA,cpoly%loss_fraction(new_lu,isi)     &
                             ,poly_dist_type,mindbh_harvest)
      end if
   end do

   return
end subroutine harv_mat_patches
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine harv_immat_patches(cpoly,isi, newp, harvest_deficit,total_harvest_area)
   use ed_state_vars     , only : polygontype          & ! structure
                                , sitetype             & ! structure
                                , patchtype            ! ! structure
   use disturb_coms      , only : plantation_rotation  & ! intent(in)
                                , mature_harvest_age   ! ! intent(in)
   use disturbance_utils , only : accum_dist_litt      & ! subroutine
                                , insert_survivors     & ! subroutine
                                , increment_patch_vars ! ! subroutine
   use ed_max_dims       , only : n_pft                ! ! intent(in)
   implicit none

   !----- Arguments -----------------------------------------------------------------------!
   type(polygontype)                  , target        :: cpoly
   integer                            , intent(in)    :: isi
   integer                            , intent(in)    :: newp
   real                               , intent(inout) :: harvest_deficit
   real                               , intent(inout) :: total_harvest_area
   !----- Local variables -----------------------------------------------------------------!
   type(sitetype)                     , pointer       :: csite
   real             , dimension(n_pft)                :: mindbh_harvest
   integer                                            :: ipa
   real                                               :: lambda
   real                                               :: dA
   !----- Local constants. ----------------------------------------------------------------!
   integer                            , parameter     :: new_lu         = 2
   integer                            , parameter     :: poly_dist_type = 1
   !---------------------------------------------------------------------------------------!


   !----- Loop over patches ---------------------------------------------------------------!
 
   csite => cpoly%site(isi)

   !---------------------------------------------------------------------------------------!
   !     David, csite%patch(npatches) is the patch that has just been created due to       !
   !  harvesting.  It doesn't have area assigned yet.  In any case, I don't think it falls !
   !  into the "immature secondary" or "immature plantation" categories since it is al-    !
   !  ready harvested... so I switched the loop to stop at csite%npatches-1 rather than    !
   !  csite%npatches. Is that correct?                                                     !
   !---------------------------------------------------------------------------------------!
   patchloop1: do ipa=1,csite%npatches-1

      !----- Skip the patch with too little biomass ---------------------------------------!
      if (csite%plant_ag_biomass(ipa) < 0.01) cycle patchloop1

      !------------------------------------------------------------------------------------!
      !    First harvest the immature secondary. This will happen only when the following  !
      ! conditions are met:                                                                !
      ! 1. There is still a deficit;                                                       !
      ! 2. Secondary forest                                                                !
      ! 3. Either immature plantation or immature secondary                                !
      !------------------------------------------------------------------------------------!
      if ( harvest_deficit > 0.0 .and. csite%dist_type(ipa) == 2 .and.                     &
          ((csite%plantation(ipa) == 1 .and. csite%age(ipa) < plantation_rotation) .or.    &
           (csite%plantation(ipa) /= 1 .and. csite%age(ipa) < mature_harvest_age) )) then

         !----- Dump the minimum DBH here.  It is an argument, but not really used here. --!
         mindbh_harvest(1:n_pft) = cpoly%mindbh_secondary(1:n_pft,isi)

         if( (csite%area(ipa) * csite%plant_ag_biomass(ipa)) > harvest_deficit)then
            !----- Patch is not totally harvested. ----------------------------------------!
            lambda = harvest_deficit / (csite%area(ipa) * csite%plant_ag_biomass(ipa))
            dA = csite%area(ipa) * lambda
            harvest_deficit = 0.0
         else
            !----- Patch is totally harvested. --------------------------------------------!
            dA = csite%area(ipa)
            harvest_deficit = harvest_deficit                                              &
                            - csite%area(ipa) * csite%plant_ag_biomass(ipa)
         end if
         total_harvest_area = total_harvest_area + dA
         csite%area(ipa)    = csite%area(ipa) - dA
         call increment_patch_vars(csite,newp,ipa, dA)
         call accum_dist_litt(csite,newp,ipa,new_lu,dA,cpoly%loss_fraction(new_lu,isi)     &
                             ,poly_dist_type,mindbh_harvest)
      end if
   end do patchloop1

   !----- Return if we have reached our harvest target. -----------------------------------!
   if(harvest_deficit <= 0.0)return

   !---------------------------------------------------------------------------------------!
   !    If we did not reach our target, loop again through patches, this time harvesting   !
   ! from immature primary.                                                                !
   !---------------------------------------------------------------------------------------!
   patchloop2: do ipa=1,csite%npatches
   
      !----- Skip the patch with too little biomass ---------------------------------------!
      if (csite%plant_ag_biomass(ipa) < 0.01) cycle patchloop2

      !------------------------------------------------------------------------------------!
      !    If necessary, harvest the immature primary. This will happen only when the      !
      ! following conditions are met:                                                      !
      ! 1. There is still a deficit;                                                       !
      ! 2. Primary forest;                                                                 !
      ! 3. It is immature.                                                                 !
      !------------------------------------------------------------------------------------!
      if (harvest_deficit > 0.0 .and. csite%dist_type(ipa) == 3 .and.                      &
          csite%age(ipa) < mature_harvest_age) then

         !----- Dump the minimum DBH here.  It is an argument, but not really used here. --!
         mindbh_harvest(1:n_pft) = cpoly%mindbh_primary(1:n_pft,isi)

         if ((csite%area(ipa) * csite%plant_ag_biomass(ipa)) > harvest_deficit) then
            
            !----- Patch is not totally harvested. ----------------------------------------!
            lambda = harvest_deficit / (csite%area(ipa) * csite%plant_ag_biomass(ipa))
            dA = csite%area(ipa) * lambda
            harvest_deficit = 0.0
         else
            !----- Patch is totally harvested. --------------------------------------------!
            dA = csite%area(ipa)
            harvest_deficit = harvest_deficit                                              &
                            - csite%area(ipa) * csite%plant_ag_biomass(ipa)
         end if

         total_harvest_area = total_harvest_area + dA
         csite%area(ipa)    = csite%area(ipa) - dA

         call increment_patch_vars(csite,newp,ipa,dA)
         call accum_dist_litt(csite,newp,ipa,new_lu,dA,cpoly%loss_fraction(new_lu,isi)     &
                             ,poly_dist_type,mindbh_harvest)
      end if
   end do patchloop2

   return
end subroutine harv_immat_patches
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine norm_harv_patch(csite,newp)

   use ed_state_vars , only : sitetype            & ! structure
                            , patchtype           ! ! structure
   use disturb_coms  , only : min_new_patch_area  ! ! intent(in)
   use ed_max_dims   , only : n_pft               ! ! intent(in)
   use grid_coms     , only : nzg                 & ! intent(in)
                            , nzs                 ! ! intent(in)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(sitetype), target     :: csite
   integer       , intent(in) :: newp
   !----- Local variables -----------------------------------------------------------------!
   real                       :: area_fac
   integer                    :: k
   !---------------------------------------------------------------------------------------!

   !----- Skip normalization when the patch is too small.  It will be terminated soon. ----!
   if (csite%area(newp) < min_new_patch_area) then
      return
   else
      !----- To make the values the weighted average of all contributing patches. ---------!
      area_fac = 1.0 / csite%area(newp)
   end if

   csite%fast_soil_C(newp)                 = csite%fast_soil_C(newp)         * area_fac
   csite%slow_soil_C(newp)                 = csite%slow_soil_C(newp)         * area_fac
   csite%structural_soil_C(newp)           = csite%structural_soil_C(newp)   * area_fac
   csite%structural_soil_L(newp)           = csite%structural_soil_L(newp)   * area_fac
   csite%mineralized_soil_N(newp)          = csite%mineralized_soil_N(newp)  * area_fac
   csite%fast_soil_N(newp)                 = csite%fast_soil_N(newp)         * area_fac
   csite%sum_dgd(newp)                     = csite%sum_dgd(newp)             * area_fac
   csite%sum_chd(newp)                     = csite%sum_chd(newp)             * area_fac
   csite%can_theta(newp)                   = csite%can_theta(newp)           * area_fac
   csite%can_theiv(newp)                   = csite%can_theiv(newp)           * area_fac
   csite%can_prss(newp)                    = csite%can_prss(newp)            * area_fac
   csite%can_co2(newp)                     = csite%can_co2(newp)             * area_fac
   csite%can_shv(newp)                     = csite%can_shv(newp)             * area_fac
   csite%can_depth(newp)                   = csite%can_depth(newp)           * area_fac
   csite%ggbare(newp)                      = csite%ggbare(newp)              * area_fac
   csite%ggveg(newp)                       = csite%ggveg(newp)               * area_fac
   csite%rough(newp)                       = csite%rough(newp)               * area_fac
   csite%mean_rh(newp)                     = csite%mean_rh(newp)             * area_fac
   csite%today_A_decomp(newp)              = csite%today_A_decomp(newp)      * area_fac
   csite%today_Af_decomp(newp)             = csite%today_Af_decomp(newp)     * area_fac
   csite%repro(1:n_pft,newp)               = csite%repro(1:n_pft,newp)       * area_fac
   csite%fsc_in(newp)                      = csite%fsc_in(newp)              * area_fac
   csite%ssc_in(newp)                      = csite%ssc_in(newp)              * area_fac
   csite%ssl_in(newp)                      = csite%ssl_in(newp)              * area_fac
   csite%fsn_in(newp)                      = csite%fsn_in(newp)              * area_fac
   csite%total_plant_nitrogen_uptake(newp) = csite%total_plant_nitrogen_uptake(newp)       &
                                           * area_fac

   do k = 1, nzs
      csite%sfcwater_mass(k,newp)   = csite%sfcwater_mass(k,newp)   * area_fac
      csite%sfcwater_energy(k,newp) = csite%sfcwater_energy(k,newp) * area_fac
      csite%sfcwater_depth(k,newp)  = csite%sfcwater_depth(k,newp)  * area_fac
   end do
   do k = 1, nzg
      csite%soil_energy(k,newp) = csite%soil_energy(k,newp) * area_fac
      csite%soil_water(k,newp)  = csite%soil_water(k,newp)  * area_fac
   end do
   return
end subroutine norm_harv_patch
!==========================================================================================!
!==========================================================================================!
