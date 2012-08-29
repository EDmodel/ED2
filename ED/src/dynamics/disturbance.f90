!==========================================================================================!
!==========================================================================================!
!    This module contains subroutines and functions that will apply disturbances to        !
! patches.  This is usually done once a year, and the main disturbance driver will decide  !
! which kind of disturbance should be applied.                                             !
!------------------------------------------------------------------------------------------!
module disturbance_utils

   use ed_state_vars   , only : allocate_patchtype    & ! subroutine
                              , copy_patchtype        & ! subroutine
                              , deallocate_patchtype  & ! subroutine
                              , allocate_sitetype     & ! subroutine
                              , deallocate_sitetype   & ! subroutine
                              , copy_sitetype_mask    ! ! subroutine
   use fuse_fiss_utils , only : fuse_cohorts          & ! subroutine
                              , terminate_cohorts     & ! subroutine
                              , split_cohorts         ! ! subroutine
   !=======================================================================================!
   !=======================================================================================!


   contains



   !=======================================================================================!
   !=======================================================================================!
   !     This is the main disturbance driver.  It will be called every New Year day, and   !
   ! it will decide whether a new patch should be created.  Disturbances can be natural    !
   ! and anthropogenic, and the new patch will be always assigned an identification flag   !
   ! that will tell how that patch was created.  Three categories are currently possible:  !
   ! 1 - agriculture: conversion to agriculture by land clearing;                          !
   ! 2 - secondary forest: logging, land abandonment, and harvest create this patch;       !
   ! 3 - primary forest: natural disturbances (treefall or fire).                          !
   !---------------------------------------------------------------------------------------!
   subroutine apply_disturbances(cgrid)

      use ed_state_vars, only : edtype                  & ! structure
                              , polygontype             & ! structure
                              , sitetype                & ! structure
                              , patchtype               ! ! structure
      use ed_misc_coms , only : current_time            & ! intent(in)
                              , ibigleaf                ! ! intent(in)
      use disturb_coms , only : min_patch_area          & ! intent(in)
                              , mature_harvest_age      & ! intent(in)
                              , plantation_rotation     & ! intent(in)
                              , ianth_disturb           & ! intent(in)
                              , time2canopy             ! ! intent(in)
      use ed_max_dims  , only : n_dist_types            & ! intent(in)
                              , n_pft                   & ! intent(in)
                              , n_dbh                   ! ! intent(in)
      use mem_polygons , only : maxcohort               ! ! intent(in)
      use grid_coms    , only : nzg                     & ! intent(in)
                              , nzs                     ! ! intent(in)
      use pft_coms     , only : include_pft             ! ! intent(in)
      use allometry    , only : area_indices            ! ! function
      use mortality    , only : disturbance_mortality   ! ! subroutine
      use consts_coms  , only : lnexp_max               ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(edtype)                   , target      :: cgrid
      !----- Local variables. -------------------------------------------------------------!
      type(polygontype)              , pointer     :: cpoly
      type(sitetype)                 , pointer     :: csite
      type(sitetype)                 , pointer     :: tsite
      type(patchtype)                , pointer     :: cpatch
      type(patchtype)                , pointer     :: qpatch
      integer                                      :: ipy
      integer                                      :: isi
      integer                                      :: ipa
      integer                                      :: ico
      integer                                      :: ipft
      integer                                      :: i
      integer                                      :: i1
      integer                                      :: i2
      integer                                      :: mypfts
      integer                                      :: onsp
      integer                                      :: old_lu
      integer                                      :: new_lu
      integer                                      :: new_add
      integer                                      :: poly_dest_type
      integer, dimension(:)          , allocatable :: pfts
      logical, dimension(:)          , allocatable :: disturb_mask
      logical, dimension(:)          , allocatable :: acceptor_mask
      logical                                      :: have_acceptor
      logical                                      :: ploughed
      logical                                      :: abandoned
      logical                                      :: natural
      logical                                      :: logged
      logical                                      :: is_plantation
      logical                                      :: mature_plantation
      logical                                      :: mature_primary
      logical                                      :: mature_secondary
      real   , dimension(n_pft,n_dbh)              :: initial_agb
      real   , dimension(n_pft,n_dbh)              :: initial_basal_area
      real   , dimension(n_pft)                    :: mindbh_harvest
      real                                         :: area
      real                                         :: area_fac
      real                                         :: dA
      real                                         :: total_distrate
      real                                         :: elim_nplant
      real                                         :: elim_lai
      real   , dimension(3)                        :: area_old
      real   , dimension(3)                        :: area_new
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Select which type of vegetation structure this run is solving.                  !
      !------------------------------------------------------------------------------------!
      select case (ibigleaf)
      case (0)
         !----- Allocating the temporary site that will host the original patches. --------!
         nullify(tsite)
         allocate(tsite)

         polyloop: do ipy = 1,cgrid%npolygons
            
            cpoly => cgrid%polygon(ipy)
            siteloop: do isi = 1,cpoly%nsites

               csite => cpoly%site(isi)

               !----- Store AGB, basal area profiles in memory. ---------------------------!
               call update_site_derived_props(cpoly, 1,isi)
               initial_agb(1:n_pft,1:n_dbh)        = cpoly%agb(1:n_pft,1:n_dbh,isi)
               initial_basal_area(1:n_pft,1:n_dbh) = cpoly%basal_area(1:n_pft,1:n_dbh,isi)

               !---------------------------------------------------------------------------!
               !      First take care of harvesting, i.e., secondary -> secondary and      !
               ! primary -> secondary.                                                     !
               !---------------------------------------------------------------------------!
               call apply_forestry(cpoly,isi, current_time%year)

               !----- Update the cut output variables. ------------------------------------!
               call update_site_derived_props(cpoly, 1,isi)
               cpoly%agb_cut(1:n_pft,1:n_dbh,isi) = cpoly%agb_cut(1:n_pft, 1:n_dbh,isi)    &
                                                  + initial_agb(1:n_pft, 1:n_dbh)          &
                                                  - cpoly%agb(1:n_pft, 1:n_dbh,isi)

               cpoly%basal_area_cut(1:n_pft,1:n_dbh,isi) =                                 &
                                                 cpoly%basal_area_cut(1:n_pft,1:n_dbh,isi) &
                                               + initial_basal_area(1:n_pft,1:n_dbh)       &
                                               - cpoly%basal_area(1:n_pft,1:n_dbh,isi)

               !----- Save the Original Number (of) Site Patches, onsp... -----------------!
               onsp = csite%npatches


               !---------------------------------------------------------------------------!
               !     Create a temporary site with vectors containing all current patches   !
               ! as well as n_dist_types patches.  Create the newly disturbed patches in   !
               ! here, and depending on how many are created, repopulate the existing      !
               ! site's patch vectors.                                                     !
               !---------------------------------------------------------------------------!
               call allocate_sitetype(tsite,onsp)

               allocate(disturb_mask(onsp + n_dist_types))
               disturb_mask         = .false.
               disturb_mask(1:onsp) = .true.

               !---------------------------------------------------------------------------!
               !     Transfer the origial patch values into the front end of the temp's    !
               ! space.                                                                    !
               !---------------------------------------------------------------------------!
               call copy_sitetype_mask(csite,tsite,disturb_mask(1:onsp)                    &
                                      ,count(disturb_mask),count(disturb_mask))

               !----- Reallocate and transfer them back. ----------------------------------!
               call deallocate_sitetype(csite)
               call allocate_sitetype(csite,onsp + n_dist_types)
               call copy_sitetype_mask(tsite,csite,disturb_mask(1:onsp)                    &
                                      ,count(disturb_mask),count(disturb_mask))
               call deallocate_sitetype(tsite)


               !---------------------------------------------------------------------------!
               !      Initialize all the potential as well as implemented disturbance      !
               ! patches.  n_dist_types new patches will be created, each one containing a !
               ! different patch type.  In case no conversion to that kind of patch has    !
               ! happened, or if the newly created patch is tiny, it will be removed soon. !
               !---------------------------------------------------------------------------!
               do new_lu = onsp+1, onsp+n_dist_types
                  call initialize_disturbed_patch(csite,cpoly%met(isi)%atm_tmp,new_lu,1    &
                                                 ,cpoly%lsl(isi))
               end do
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !       Reset mortality rates due to disturbance.                           !
               !---------------------------------------------------------------------------!
               resetdist: do ipa=1,onsp
                  cpatch => csite%patch(ipa)
                  do ico=1,cpatch%ncohorts
                     cpatch%mort_rate(5,ico) = 0.0
                  end do
               end do resetdist
               !---------------------------------------------------------------------------!

               !----- Loop over q, the *destination* landuse type. ------------------------!
               new_lu_loop: do new_lu = 1, n_dist_types
                  !----- Set up area to zero, in case no conversion happens. --------------!
                  area = 0.0

                  do ipa=1,onsp
                     !----- Save the old land use in a shorter variable for convenience. --!
                     old_lu        = csite%dist_type(ipa)
                     is_plantation = csite%plantation(ipa) == 1

                     !---------------------------------------------------------------------!
                     !    Now we add the area associated with each kind of possible        !
                     ! disturbance that can happen.  Types of conversion that are solved   !
                     ! here are:                                                           !
                     ! * ploughed  - conversion from primary/secondary land to agriculture.!
                     ! * abandoned - conversion from agriculture to secondary land.        !
                     ! * natural   - natural disturbance from primary/secondary land to    !
                     !               primary land (fires or tree fall)                     !
                     ! * logged    - conversion from primary/secondary land to secondary   !
                     !               land due to logging.                                  !
                     !---------------------------------------------------------------------!
                     ploughed  = new_lu == 1 .and. old_lu /= 1
                     abandoned = new_lu == 2 .and. old_lu == 1
                     !----- Natural disturbance, either trees are old or there is a fire. -!
                     natural   = new_lu == 3 .and. old_lu /= 1 .and.                       &
                                 ( csite%age(ipa) > time2canopy .or.                       &
                                   cpoly%nat_dist_type(isi) == 1)
                     !----- Check whether the patch is ready  be harvested. ---------------!
                     mature_primary    = old_lu == 3                               .and.   &
                                         csite%age(ipa)         > mature_harvest_age
                     mature_plantation = is_plantation                             .and.   &
                                         csite%age(ipa)         > plantation_rotation
                     mature_secondary  = old_lu == 2 .and. (.not. is_plantation)   .and.   &
                                         csite%age(ipa)         > mature_harvest_age
                     !---------------------------------------------------------------------!
                     logged    = new_lu == 2 .and.                                         &
                                 ( mature_primary .or. mature_plantation .or.              &
                                   mature_secondary)

                     !---------------------------------------------------------------------!
                     !    Add area if any of the disturbances that produce of type new_lu  !
                     ! has happened.  The ones that produce type other than new_lu will be !
                     ! always false for new_lu.                                            !
                     !---------------------------------------------------------------------!
                     if  (ploughed .or. abandoned .or. natural .or. logged) then
                        total_distrate = cpoly%disturbance_rates (new_lu,old_lu,isi)       &
                                       + cpoly%disturbance_memory(new_lu,old_lu,isi)
                        total_distrate = min(lnexp_max,max(0.,total_distrate))
                        dA             = csite%area(ipa) * ( 1.0 - exp(- total_distrate) )
                        area           = area + dA
                     end if
                     !---------------------------------------------------------------------!


                  end do
                  
                  if (area > min_patch_area) then
                     write(unit=*,fmt='(a,1x,es12.5,1x,a,1x,i5)')                          &
                         ' ---> Making new patch, with area=',area,' for dist_type=',new_lu

                     !---------------------------------------------------------------------!
                     !     Set the flag that this patch should be kept as a newly created  !
                     ! transition patch.                                                   !
                     !---------------------------------------------------------------------!
                     disturb_mask(onsp+new_lu)     = .true.
     
                     csite%dist_type(onsp+new_lu)  = new_lu
                     csite%plantation(onsp+new_lu) = 0
                     csite%area(onsp+new_lu)       = area
                     
                     !---------------------------------------------------------------------!
                     !     Initialize to zero the new trasitioned patches.                 !
                     !---------------------------------------------------------------------!
                     call initialize_disturbed_patch(csite,cpoly%met(isi)%atm_tmp          &
                                                    ,onsp+new_lu,1,cpoly%lsl(isi))

                     !---------------------------------------------------------------------!
                     !    Now go through patches, adding its contribution to the new       !
                     ! patch.                                                              !
                     !---------------------------------------------------------------------!
                     do ipa=1,onsp
                        !------------------------------------------------------------------!
                        !     Save the old land use in a shorter variable for convenience. !
                        !------------------------------------------------------------------!
                        old_lu        = csite%dist_type(ipa)
                        is_plantation = csite%plantation(ipa) == 1

                        !----- Check whether this patch can be disturbed. -----------------!
                        ploughed  = new_lu == 1 .and. old_lu /= 1
                        abandoned = new_lu == 2 .and. old_lu == 1
                        natural   = new_lu == 3 .and. old_lu /= 1 .and.                    &
                                    ( csite%age(ipa) > time2canopy .or.                    &
                                      cpoly%nat_dist_type(isi) == 1)
                        !----- Check whether the patch is ready  be harvested. ------------!
                        mature_primary    = old_lu == 3                              .and. &
                                            csite%age(ipa)         > mature_harvest_age
                        mature_plantation = is_plantation                            .and. &
                                            csite%age(ipa)         > plantation_rotation
                        mature_secondary  = old_lu == 2 .and. (.not. is_plantation)  .and. &
                                            csite%age(ipa)         > mature_harvest_age
                        !------------------------------------------------------------------!
                        logged    = new_lu == 2 .and.                                      &
                                    ( mature_primary .or. mature_plantation .or.           &
                                      mature_secondary)

                        !------------------------------------------------------------------!
                        !     Adjust some information to be sent to the disturbance        !
                        ! routine.                                                         !
                        !------------------------------------------------------------------!
                        if (natural) then
                           poly_dest_type          = cpoly%nat_dist_type(isi)
                           mindbh_harvest(1:n_pft) = huge(1.)
                        elseif (logged .and. mature_primary) then
                           poly_dest_type          = 2
                           mindbh_harvest(1:n_pft) = cpoly%mindbh_primary(1:n_pft,isi)
                        elseif (logged) then
                           poly_dest_type          = 2
                           mindbh_harvest(1:n_pft) = cpoly%mindbh_secondary(1:n_pft,isi)
                        else
                           poly_dest_type          = 0
                           mindbh_harvest(1:n_pft) = huge(1.)
                        end if
                        !------------------------------------------------------------------!



                        !------------------------------------------------------------------!
                        !    If the patch is going to be disturbed, compute the area of    !
                        ! the disturbed patch to be added to the new destination patch and !
                        ! update the litter layer.                                         !
                        !------------------------------------------------------------------!
                        if (ploughed .or. abandoned .or. natural .or. logged) then
                           total_distrate = cpoly%disturbance_rates (new_lu,old_lu,isi)    &
                                          + cpoly%disturbance_memory(new_lu,old_lu,isi)
                           total_distrate = min(lnexp_max,max(0.,total_distrate))
                           dA             = csite%area(ipa) * (1.0 - exp(- total_distrate))

                           area_fac = dA / csite%area(onsp+new_lu)
                           call disturbance_mortality(csite,ipa,total_distrate,new_lu      &
                                                     ,poly_dest_type,mindbh_harvest)
                           call increment_patch_vars(csite,new_lu+onsp,ipa,area_fac)
                           call insert_survivors(csite,new_lu+onsp,ipa,new_lu,area_fac     &
                                                ,poly_dest_type,mindbh_harvest)
                           call accum_dist_litt(csite,new_lu+onsp,ipa,new_lu,area_fac      &
                                               ,poly_dest_type,mindbh_harvest)

                           !----- Update patch area. --------------------------------------!
                           csite%area(ipa) = csite%area(ipa) - dA
                        end if
                     end do

                     !---------------------------------------------------------------------!
                     !      Update temperature and density.  This must be done before      !
                     ! planting, since the leaf temperature is initially assigned as the   !
                     ! canopy air temperature.                                             !
                     !---------------------------------------------------------------------!
                     call update_patch_thermo_props(csite,new_lu+onsp,new_lu+onsp,nzg,nzs  &
                                                   ,cpoly%ntext_soil(:,isi))

                     !----- If the new patch is agriculture, plant it with grasses. -------!
                     if (new_lu == 1) then 
                        call plant_patch(csite,new_lu+onsp,nzg                             &
                                        ,cpoly%agri_stocking_pft(isi)                      &
                                        ,cpoly%agri_stocking_density(isi)                  &
                                        ,cpoly%ntext_soil(:,isi)                           &
                                        ,cpoly%green_leaf_factor(:,isi), 1.0               &
                                        ,cpoly%lsl(isi))
                     end if

                     qpatch => csite%patch(new_lu+onsp)

                     !----- Fuse then terminate cohorts. ----------------------------------!
                     if (csite%patch(new_lu+onsp)%ncohorts > 0 .and. maxcohort >= 0) then
                        call fuse_cohorts(csite,new_lu+onsp,cpoly%green_leaf_factor(:,isi) &
                                         ,cpoly%lsl(isi))
                        call terminate_cohorts(csite,new_lu+onsp,elim_nplant,elim_lai)
                        call split_cohorts(qpatch,cpoly%green_leaf_factor(:,isi)           &
                                          ,cpoly%lsl(isi))
                     end if
                 
                     !----- Store AGB, basal area profiles in memory. ---------------------!
                     initial_agb(1:n_pft,1:n_dbh) = cpoly%agb(1:n_pft, 1:n_dbh,isi)
                     initial_basal_area(1:n_pft,1:n_dbh) =                                 &
                                                      cpoly%basal_area(1:n_pft,1:n_dbh,isi)

                     !---------------------------------------------------------------------!
                     !     Update the derived properties including veg_height, and patch-  !
                     ! -level LAI, WAI.                                                    !
                     !---------------------------------------------------------------------!
                     call update_patch_derived_props( csite,cpoly%lsl(isi)                 &
                                                    , cpoly%met(isi)%prss                  &
                                                    , new_lu+onsp)
                     !----- Update soil temperature, liquid fraction, etc. ----------------!
                     call new_patch_sfc_props(csite,new_lu+onsp,nzg,nzs                    &
                                             ,cpoly%ntext_soil(:,isi))
                     !----- Update budget properties. -------------------------------------!
                     call update_budget(csite,cpoly%lsl(isi),new_lu+onsp,new_lu+onsp)

                     !----- Update AGB, basal area. ---------------------------------------!
                     call update_site_derived_props(cpoly,1,isi)

                     !----- Update either cut or mortality. -------------------------------!
                     if (new_lu /= 3) then
                        cpoly%agb_cut(1:n_pft,1:n_dbh,isi) =                               &
                               cpoly%agb_cut(1:n_pft, 1:n_dbh,isi)                         &
                             + initial_agb(1:n_pft, 1:n_dbh)                               &
                             - cpoly%agb(1:n_pft, 1:n_dbh,isi)
                        cpoly%basal_area_cut(1:n_pft, 1:n_dbh,isi) =                       &
                               cpoly%basal_area_cut(1:n_pft, 1:n_dbh,isi)                  &
                             + initial_basal_area(1:n_pft, 1:n_dbh)                        &
                             - cpoly%basal_area(1:n_pft, 1:n_dbh,isi)
                     else
                        cpoly%agb_mort(1:n_pft,1:n_dbh,isi) =                              &
                               cpoly%agb_mort(1:n_pft,1:n_dbh,isi)                         &
                             + initial_agb(1:n_pft,1:n_dbh)                                &
                             - cpoly%agb(1:n_pft,1:n_dbh,isi)
                        cpoly%basal_area_mort(1:n_pft, 1:n_dbh,isi) =                      &
                               cpoly%basal_area_mort(1:n_pft, 1:n_dbh,isi)                 &
                             + initial_basal_area(1:n_pft, 1:n_dbh)                        &
                             - cpoly%basal_area(1:n_pft, 1:n_dbh,isi)
                     endif
                     
                     !----- Clear the disturbance memory for this disturbance type. -------!
                     cpoly%disturbance_memory(new_lu,1:n_dist_types,isi) = 0.0

                  elseif(area > 0.0)then
                     !---------------------------------------------------------------------!
                     !     The patch creation has been skipped because the area was too    !
                     ! small.  Put the current disturbance rates in memory to be added at  !
                     ! the next timestep.                                                  !
                     !---------------------------------------------------------------------!
                     cpoly%disturbance_memory(new_lu,1:n_dist_types,isi) =                 &
                            cpoly%disturbance_memory(new_lu,1:n_dist_types,isi)            &
                          + cpoly%disturbance_rates(new_lu,1:n_dist_types,isi)
                  end if
               end do new_lu_loop

               !---------------------------------------------------------------------------!
               !      Reallocate the current site to fit the original patches and whatever !
               ! was generated in disturbance (ie, make it non sparse).  Populate the      !
               ! original site with both the modified original patches, and the newly      !
               ! created patches.  The index of all of these are disturb_mask.  This mask  !
               ! should be one for all original patches, and sparse from there after.      !
               !---------------------------------------------------------------------------!
               call allocate_sitetype(tsite,count(disturb_mask))
               call copy_sitetype_mask(csite,tsite,disturb_mask,size(disturb_mask)         &
                                      ,count(disturb_mask))
               call deallocate_sitetype(csite)
               call allocate_sitetype(csite,count(disturb_mask))
             
               disturb_mask = .false.
               disturb_mask(1:csite%npatches) = .true.
               call copy_sitetype_mask(tsite,csite,disturb_mask(1:csite%npatches)          &
                                      ,count(disturb_mask),count(disturb_mask))
             
               call deallocate_sitetype(tsite)
               deallocate(disturb_mask)
               !---------------------------------------------------------------------------!

             
            end do siteloop
         end do polyloop

         !----- Free memory before leaving... ---------------------------------------------!
         deallocate(tsite)
      case (1)
         !---------------------------------------------------------------------------------!
         !      Big leaf.                                                                  !
         !---------------------------------------------------------------------------------!
         !----- Allocate the temporary site that will host the original patches. ----------!
         nullify (tsite)
         allocate(tsite)

         polyloop2: do ipy = 1,cgrid%npolygons

            cpoly => cgrid%polygon(ipy)
            siteloop2: do isi = 1,cpoly%nsites

               csite => cpoly%site(isi)

               !----- Store AGB, basal area profiles in memory. ---------------------------!
               call update_site_derived_props(cpoly, 1,isi)
               initial_agb(1:n_pft,1:n_dbh)        = cpoly%agb(1:n_pft,1:n_dbh,isi)
               initial_basal_area(1:n_pft,1:n_dbh) = cpoly%basal_area(1:n_pft,1:n_dbh,isi)

               !---------------------------------------------------------------------------!
               !      First take care of harvesting, i.e., secondary -> secondary and      !
               ! primary -> secondary.                                                     !
               !---------------------------------------------------------------------------!
               call apply_forestry(cpoly,isi, current_time%year)

               !----- Update the cut output variables. ------------------------------------!
               call update_site_derived_props(cpoly, 1,isi)
               cpoly%agb_cut(1:n_pft,1:n_dbh,isi) = cpoly%agb_cut(1:n_pft, 1:n_dbh,isi)    &
                                                  + initial_agb(1:n_pft, 1:n_dbh)          &
                                                  - cpoly%agb(1:n_pft, 1:n_dbh,isi)

               cpoly%basal_area_cut(1:n_pft,1:n_dbh,isi) =                                 &
                                                 cpoly%basal_area_cut(1:n_pft,1:n_dbh,isi) &
                                               + initial_basal_area(1:n_pft,1:n_dbh)       &
                                               - cpoly%basal_area(1:n_pft,1:n_dbh,isi)

               !----- Save the Original Number (of) Site Patches, onsp... -----------------!
               onsp = csite%npatches


               !---------------------------------------------------------------------------!
               !     Create a temporary site with vectors containing all current patches   !
               ! as well as n_dist_types * mypfts patches. Create the newly disturbed      !
               ! patches here and repopulate the existing site's patch vectors.            !
               !---------------------------------------------------------------------------!
               mypfts=count(include_pft)
               allocate(pfts(mypfts))
               i1=0
               do i=1,n_pft
                  if (include_pft(i)) then
                     i1=i1+1
                     pfts(i1)=i
                  end if
               end do
               call allocate_sitetype(tsite,onsp)

               allocate(disturb_mask(onsp + (n_dist_types-1) * mypfts +1 ))
               disturb_mask         = .false.
               disturb_mask(1:onsp) = .true.

               !---------------------------------------------------------------------------!
               !     Transfer the origial patch values into the front end of the temp's    !
               ! space.                                                                    !
               !---------------------------------------------------------------------------!
               call copy_sitetype_mask(csite,tsite,disturb_mask(1:onsp)                    &
                                      ,count(disturb_mask),count(disturb_mask))

               !----- Reallocate and transfer them back. ----------------------------------!
               call deallocate_sitetype(csite)
               call allocate_sitetype(csite,onsp + (n_dist_types-1) * mypfts +1)
               call copy_sitetype_mask(tsite,csite,disturb_mask(1:onsp)                    &
                                      ,count(disturb_mask),count(disturb_mask))
               call deallocate_sitetype(tsite)


               !---------------------------------------------------------------------------!
               !      Initialize all the potential as well as implemented disturbance      !
               ! patches.  n_dist_types new patches will be created, each one containing a !
               ! different patch type.  In case no conversion to that kind of patch has    !
               ! happened, or if the newly created patch is tiny, it will be removed soon. !
               !---------------------------------------------------------------------------!
               do new_lu = onsp+1, onsp+(n_dist_types-1)*mypfts+1
                  call initialize_disturbed_patch(csite,cpoly%met(isi)%atm_tmp,new_lu,1    &
                                                 ,cpoly%lsl(isi))
               end do



               !---------------------------------------------------------------------------!
               !    First round, determine areas disturbanced and change the areas of the  !
               ! existing patches.                                                         !
               !---------------------------------------------------------------------------!
               
               !----- Loop over q, the *destination* landuse type. ------------------------!
               new_lu_loop2: do new_lu = 1, n_dist_types
                  !----- Set up area to zero, in case no conversion happens. --------------!
                  area = 0.0

                  do ipa=1,onsp
                     !----- Save the old land use in a shorter variable for convenience. --!
                     old_lu        = csite%dist_type(ipa)
                     is_plantation = csite%plantation(ipa) == 1

                     !---------------------------------------------------------------------!
                     !    Now we add the area associated with each kind of possible        !
                     ! disturbance that can happen.  Types of conversion that are solved   !
                     ! here are:                                                           !
                     ! * ploughed  - conversion from primary/secondary land to agriculture.!
                     ! * abandoned - conversion from agriculture to secondary land.        !
                     ! * natural   - natural disturbance from primary/secondary land to    !
                     !               primary land (fires or tree fall)                     !
                     ! * logged    - conversion from primary/secondary land to secondary   !
                     !               land due to logging.                                  !
                     !---------------------------------------------------------------------!
                     ploughed  = new_lu == 1 .and. old_lu /= 1
                     abandoned = new_lu == 2 .and. old_lu == 1
                     !----- Natural disturbance, either trees are old or there is a fire. -!
                     natural   = new_lu == 3 .and. old_lu /= 1 .and.                       &
                                 ( csite%age(ipa) > time2canopy .or.                       &
                                   cpoly%nat_dist_type(isi) == 1)
                     !----- Check whether the patch is ready  be harvested. ---------------!
                     mature_primary    = old_lu == 3                                 .and. &
                                         csite%age(ipa)         > mature_harvest_age
                     mature_plantation = is_plantation                               .and. &
                                         csite%age(ipa)         > plantation_rotation
                     mature_secondary  = old_lu == 2 .and. (.not. is_plantation)     .and. &
                                         csite%age(ipa)         > mature_harvest_age
                     !---------------------------------------------------------------------!
                     logged    = new_lu == 2 .and.                                         &
                                 ( mature_primary .or. mature_plantation .or.              &
                                   mature_secondary)

                     !---------------------------------------------------------------------!
                     !    Add area if any of the disturbances that produce of type new_lu  !
                     ! has happened.  The ones that produce type other than new_lu will be !
                     ! always false for new_lu.                                            !
                     !---------------------------------------------------------------------!
                     if  (ploughed .or. abandoned .or. natural .or. logged) then
                        total_distrate = cpoly%disturbance_rates (new_lu,old_lu,isi)       &
                                       + cpoly%disturbance_memory(new_lu,old_lu,isi)
                        total_distrate = min(lnexp_max,max(0.,total_distrate))
                        dA             = csite%area(ipa) * (1.0 - exp(- total_distrate))
                        area = area + dA
                     end if
                     !---------------------------------------------------------------------!


                  end do

                  if (area > min_patch_area) then

                     !---------------------------------------------------------------------!
                     !     Set the flag that this patch should be kept as a newly created  !
                     ! transition patch.                                                   !
                     !---------------------------------------------------------------------!
                     if(new_lu == 1) then
                        disturb_mask(onsp+new_lu)      = .true.
                        csite%dist_type(onsp+new_lu)   = new_lu
                        csite%plantation(onsp+new_lu)  = 0
                        csite%area(onsp+new_lu)        = area
                        !----- Initialize to zero the new trasitioned patches. ------------!
                        call initialize_disturbed_patch(csite,cpoly%met(isi)%atm_tmp       &
                                                    ,onsp+new_lu,1,cpoly%lsl(isi))
                     else
                        do ipft=1,mypfts
                           i                           = (new_lu-2)*mypfts+1+ipft
                           disturb_mask(onsp+i)        = .true.
                           csite%dist_type(onsp+i)     = new_lu
                           csite%plantation(onsp+i)    = 0
                           csite%area(onsp+i)          = area / real(mypfts)
                           !----- Initialize to zero the new trasitioned patches. ---------!
                           call initialize_disturbed_patch(csite,cpoly%met(isi)%atm_tmp    &
                                                    ,onsp+i,1,cpoly%lsl(isi))
                        end do
                     end if
                     !---------------------------------------------------------------------!


                     !---------------------------------------------------------------------!
                     !    Now go through patches, adding its contribution to the new       !
                     ! patch.                                                              !
                     !---------------------------------------------------------------------!
                     do ipa=1,onsp
                        !------------------------------------------------------------------!
                        !     Save the old land use in a shorter variable for convenience. !
                        !------------------------------------------------------------------!
                        old_lu        = csite%dist_type(ipa)
                        is_plantation = csite%plantation(ipa) == 1

                        !----- Check whether this patch can be disturbed. -----------------!
                        ploughed  = new_lu == 1 .and. old_lu /= 1
                        abandoned = new_lu == 2 .and. old_lu == 1
                        natural   = new_lu == 3 .and. old_lu /= 1 .and.                    &
                                    ( csite%age(ipa) > time2canopy .or.                    &
                                      cpoly%nat_dist_type(isi) == 1)
                        !----- Check whether the patch is ready  be harvested. ------------!
                        mature_primary    = old_lu == 3                              .and. &
                                            csite%age(ipa)         > mature_harvest_age
                        mature_plantation = is_plantation                            .and. &
                                            csite%age(ipa)         > plantation_rotation
                        mature_secondary  = old_lu == 2 .and. (.not. is_plantation)  .and. &
                                            csite%age(ipa)         > mature_harvest_age
                        !------------------------------------------------------------------!
                        logged    = new_lu == 2 .and.                                      &
                                    ( mature_primary .or. mature_plantation .or.           &
                                      mature_secondary)

                        !------------------------------------------------------------------!
                        !     Adjust some information to be sent to the disturbance        !
                        !  routine.                                                        !
                        !------------------------------------------------------------------!
                        if (natural) then
                           poly_dest_type          = cpoly%nat_dist_type(isi)
                           mindbh_harvest(1:n_pft) = huge(1.)
                        elseif (logged .and. mature_primary) then
                           poly_dest_type          = 2
                           mindbh_harvest(1:n_pft) = cpoly%mindbh_primary(1:n_pft,isi)
                        elseif (logged) then
                           poly_dest_type          = 2
                           mindbh_harvest(1:n_pft) = cpoly%mindbh_secondary(1:n_pft,isi)
                        else
                           poly_dest_type          = 0
                           mindbh_harvest(1:n_pft) = huge(1.)
                        end if
                        !------------------------------------------------------------------!



                        !------------------------------------------------------------------!
                        !    If the patch is going to be disturbed, compute the area of    !
                        ! the disturbed patch to be added to the new destination patch and !
                        ! update the litter layer.                                         !
                        !------------------------------------------------------------------!
                        if (ploughed) then
                           total_distrate = cpoly%disturbance_rates (new_lu,old_lu,isi)    &
                                          + cpoly%disturbance_memory(new_lu,old_lu,isi)
                           total_distrate = min(lnexp_max,max(0.,total_distrate))
                           dA             = csite%area(ipa) * (1.0 - exp(- total_distrate))

                           area_fac = dA / csite%area(onsp+new_lu)
                           call disturbance_mortality(csite,ipa,total_distrate,new_lu      &
                                                     ,poly_dest_type,mindbh_harvest)
                           call increment_patch_vars(csite,new_lu+onsp,ipa,area_fac)
                           call accum_dist_litt(csite,new_lu+onsp,ipa,new_lu,area_fac      &
                                               ,poly_dest_type,mindbh_harvest)

                           !----- Update patch area. --------------------------------------!
                           csite%area(ipa) = csite%area(ipa) - dA
                        else if(abandoned .or. natural .or. logged)then
                           do ipft=1,mypfts
                              i   = (new_lu-2)*mypfts+1+ipft
                              total_distrate = cpoly%disturbance_rates (new_lu,old_lu,isi) &
                                             + cpoly%disturbance_memory(new_lu,old_lu,isi)
                              total_distrate = min(lnexp_max,max(0.,total_distrate))
                              dA             = csite%area(ipa)                             &
                                             * (1.0 - exp(- total_distrate))
                              dA             = dA / real(mypfts)

                              area_fac = dA / csite%area(onsp+i)
                              call disturbance_mortality(csite,ipa,total_distrate,new_lu   &
                                                        ,poly_dest_type,mindbh_harvest)
                              call increment_patch_vars(csite,onsp+i,ipa,area_fac)
                              call accum_dist_litt(csite,onsp+i,ipa,new_lu,area_fac        &
                                                  ,poly_dest_type,mindbh_harvest)
                           end do
                           !----- Update patch area. --------------------------------------!
                           csite%area(ipa) = csite%area(ipa) - dA
                        end if
                     end do

                     !----- Clear the disturbance memory for this disturbance type. -------!
                     cpoly%disturbance_memory(new_lu,1:n_dist_types,isi) = 0.0

                  elseif(area > 0.0)then
                     !---------------------------------------------------------------------!
                     !     The patch creation has been skipped because the area was too    !
                     ! small.  Put the current disturbance rates in memory to be added at  !
                     ! the next timestep.                                                  !
                     !---------------------------------------------------------------------!
                     cpoly%disturbance_memory(new_lu,1:n_dist_types,isi) =                 &
                            cpoly%disturbance_memory(new_lu,1:n_dist_types,isi)            &
                          + cpoly%disturbance_rates(new_lu,1:n_dist_types,isi)
                  end if
               end do new_lu_loop2



               !---------------------------------------------------------------------------!
               !    Now we try to find the acceptor patch(es). If there is no one,         !
               ! we create them.                                                           !
               !---------------------------------------------------------------------------!
               do new_lu = 1, n_dist_types
                  !----- Save the old land use in a shorter variable for convenience. -----!
                  old_lu        = csite%dist_type(ipa)
                  is_plantation = csite%plantation(ipa) == 1

                  !----- Check whether this patch can be disturbed. -----------------------!
                  ploughed  = new_lu == 1 .and. old_lu /= 1
                  abandoned = new_lu == 2 .and. old_lu == 1
                  natural   = new_lu == 3 .and. old_lu /= 1 .and.                          &
                              ( csite%age(ipa) > time2canopy .or.                          &
                                cpoly%nat_dist_type(isi) == 1)
                  !----- Check whether the patch is ready  be harvested. ------------------!
                  mature_primary    = old_lu == 3                                 .and.    &
                                      csite%age(ipa)         > mature_harvest_age
                  mature_plantation = is_plantation                               .and.    &
                                      csite%age(ipa)         > plantation_rotation
                  mature_secondary  = old_lu == 2 .and. (.not. is_plantation)     .and.    &
                                      csite%age(ipa)         > mature_harvest_age
                  !------------------------------------------------------------------------!
                  logged    = new_lu == 2 .and.                                            &
                              ( mature_primary .or. mature_plantation .or.                 &
                                mature_secondary)

                  !------------------------------------------------------------------------!
                  !     Adjust some information to be sent to the disturbance routine.     !
                  !------------------------------------------------------------------------!
                  if (natural) then
                     poly_dest_type          = cpoly%nat_dist_type(isi)
                     mindbh_harvest(1:n_pft) = huge(1.)
                  elseif (logged .and. mature_primary) then
                     poly_dest_type          = 2
                     mindbh_harvest(1:n_pft) = cpoly%mindbh_primary(1:n_pft,isi)
                  elseif (logged) then
                     poly_dest_type          = 2
                     mindbh_harvest(1:n_pft) = cpoly%mindbh_secondary(1:n_pft,isi)
                  else
                     poly_dest_type          = 0
                     mindbh_harvest(1:n_pft) = huge(1.)
                  end if
                  !------------------------------------------------------------------------!
                  
                  if(new_lu ==1)then
                     i1=1
                     i2=1
                  else
                     i1=(new_lu-2)*mypfts+1+1
                     i2=(new_lu-2)*mypfts+1+mypfts
                  end if
                  
                  do i=i1,i2
                     if (disturb_mask(onsp+i)) then
                        area_old = 0.0
                        allocate(acceptor_mask(onsp))
                        acceptor_mask = .false.
                        do ipa = 1,onsp
                           area_old(csite%dist_type(ipa)) = area_old(csite%dist_type(ipa)) &
                                                          + csite%area(ipa)
                           if ( csite%dist_type(ipa) == new_lu ) then
                              acceptor_mask(ipa) = .true.
                           end if
                        end do

                        if (count(acceptor_mask) > 0) then
                           do ipa =1,onsp
                              if (acceptor_mask(ipa)) then
                                 dA       = csite%area(ipa) * csite%area(onsp+i)           &
                                          / area_old(new_lu)
                                 area_fac = csite%area(ipa) / (csite%area(ipa) + dA)
                                 call normal_patch_vars(csite,ipa,area_fac)

                                 area_fac = dA / (csite%area(ipa) + dA)
                                 call increment_patch_vars(csite,ipa,i+onsp,area_fac)
                                 call accum_dist_litt(csite,ipa,i+onsp,new_lu,area_fac     &
                                                ,poly_dest_type,mindbh_harvest)
                                 area_fac = csite%area(ipa) / (csite%area(ipa) + dA)
                                 cpatch=>csite%patch(ipa)
                                 do ico=1,cpatch%ncohorts
                                    cpatch%nplant(ico)= cpatch%nplant(ico) * area_fac

                                    !------ Compute all area indices needed. --------------!
                                    call area_indices(cpatch%nplant(ico),cpatch%bleaf(ico) &
                                                     ,cpatch%bdead(ico),cpatch%balive(ico) &
                                                     ,cpatch%dbh(ico),cpatch%hite(ico)     &
                                                     ,cpatch%pft(ico),cpatch%sla(ico)      &
                                                     ,cpatch%lai(ico),cpatch%wai(ico)      &
                                                     ,cpatch%crown_area(ico)               &
                                                     ,cpatch%bsapwooda(ico))

                                 end do
                                 csite%area(ipa)     = csite%area(ipa)     / area_fac
                                 csite%age(ipa)      = csite%age(ipa)      * area_fac

                                 disturb_mask(onsp+i) = .false.
                               end if
                           end do
                        else
                           !---------------------------------------------------------------!
                           !      Update temperature and density.                          !
                           !---------------------------------------------------------------!
                           call update_patch_thermo_props(csite,i+onsp,i+onsp,nzg          &
                                                         ,nzs,cpoly%ntext_soil(:,isi))

                           if (new_lu == 1) then
                              !------------------------------------------------------------!
                              !    If the new patch is agriculture, plant it with grasses. !
                              !------------------------------------------------------------!
                              call plant_patch(csite,i+onsp,nzg                            &
                                              ,cpoly%agri_stocking_pft(isi)                &
                                              ,cpoly%agri_stocking_density(isi)            &
                                              ,cpoly%ntext_soil(:,isi)                     &
                                              ,cpoly%green_leaf_factor(:,isi), 1.0         &
                                              ,cpoly%lsl(isi))
                              write(unit=*,fmt='(a,1x,es12.5,1x,a,1x,i5,1x,a,1x,i5)')      &
                              ' ---> Making new patch, with area=',csite%area(i+onsp)      &
                                  ,' for dist_type=',new_lu                                &
                                  ,' with pft=',cpoly%agri_stocking_pft(isi)
                           else
                              call plant_patch(csite,i+onsp,nzg,pfts(i)                    &
                                              ,cpoly%agri_stocking_density(isi)            &
                                              ,cpoly%ntext_soil(:,isi)                     &
                                              ,cpoly%green_leaf_factor(:,isi), 1.0         &
                                              ,cpoly%lsl(isi))
                              write(unit=*,fmt='(a,1x,es12.5,1x,a,1x,i5,1x,a,1x,i5)')      &
                              ' ---> Making new patch, with area=',csite%area(i+onsp),     &
                              ' for dist_type=',new_lu,' with pft=',pfts(i)
                           end if
                           !---------------------------------------------------------------!
                           !     Update the derived properties including veg_height, and   !
                           ! patch--level LAI, WAI.                                         !
                           !---------------------------------------------------------------!
                           call update_patch_derived_props(csite,cpoly%lsl(isi)            &
                                                          ,cpoly%met(isi)%prss             &
                                                          ,i+onsp)
                           !----- Update soil temperature, liquid fraction, etc. ----------!
                           call new_patch_sfc_props(csite,i+onsp,nzg,nzs                   &
                                                   ,cpoly%ntext_soil(:,isi))
                           !----- Update budget properties. -------------------------------!
                           call update_budget(csite,cpoly%lsl(isi),i+onsp,i+onsp)

                        end if
                        !----- Store AGB, basal area profiles in memory. ------------------!
                        initial_agb(1:n_pft,1:n_dbh) = cpoly%agb(1:n_pft, 1:n_dbh,isi)
                        initial_basal_area(1:n_pft,1:n_dbh) =                              &
                                                     cpoly%basal_area(1:n_pft,1:n_dbh,isi)


                        !----- Update AGB, basal area. ------------------------------------!
                        call update_site_derived_props(cpoly,1,isi)

                        !----- Update either cut or mortality. ----------------------------!
                        if (new_lu /= 3) then
                           cpoly%agb_cut(1:n_pft,1:n_dbh,isi) =                            &
                                  cpoly%agb_cut(1:n_pft, 1:n_dbh,isi)                      &
                                + initial_agb(1:n_pft, 1:n_dbh)                            &
                                - cpoly%agb(1:n_pft, 1:n_dbh,isi)
                           cpoly%basal_area_cut(1:n_pft, 1:n_dbh,isi) =                    &
                                  cpoly%basal_area_cut(1:n_pft, 1:n_dbh,isi)               &
                                + initial_basal_area(1:n_pft, 1:n_dbh)                     &
                                - cpoly%basal_area(1:n_pft, 1:n_dbh,isi)
                        else
                           cpoly%agb_mort(1:n_pft,1:n_dbh,isi) =                           &
                                  cpoly%agb_mort(1:n_pft,1:n_dbh,isi)                      &
                                + initial_agb(1:n_pft,1:n_dbh)                             &
                                - cpoly%agb(1:n_pft,1:n_dbh,isi)
                           cpoly%basal_area_mort(1:n_pft, 1:n_dbh,isi) =                   &
                                  cpoly%basal_area_mort(1:n_pft, 1:n_dbh,isi)              &
                                + initial_basal_area(1:n_pft, 1:n_dbh)                     &
                                - cpoly%basal_area(1:n_pft, 1:n_dbh,isi)
                        end if
                     
                        deallocate(acceptor_mask)
                     end if
                  end do
               end do
               !---------------------------------------------------------------------------!
               !      Reallocate the current site to fit the original patches and whatever !
               ! was generated in disturbance (ie, make it non sparse).  Populate the      !
               ! original site with both the modified original patches, and the newly      !
               ! created patches.  The index of all of these are disturb_mask.  This mask  !
               ! should be ones for all original patches, and sparse from there after.     !
               !---------------------------------------------------------------------------!
               call allocate_sitetype(tsite,count(disturb_mask))
               call copy_sitetype_mask(csite,tsite,disturb_mask,size(disturb_mask)         &
                                      ,count(disturb_mask))
               call deallocate_sitetype(csite)
               call allocate_sitetype(csite,count(disturb_mask))

               disturb_mask = .false.
               disturb_mask(1:csite%npatches) = .true.
               call copy_sitetype_mask(tsite,csite,disturb_mask(1:csite%npatches)          &
                                      ,count(disturb_mask),count(disturb_mask))

               call deallocate_sitetype(tsite)
               deallocate(disturb_mask)
               deallocate(pfts)
               !---------------------------------------------------------------------------!


            end do siteloop2
         end do polyloop2

         !----- Free memory before leaving... ---------------------------------------------!
         deallocate(tsite)

      end select
      !------------------------------------------------------------------------------------!

      return
   end subroutine apply_disturbances
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   subroutine site_disturbance_rates(month, year, cgrid)

      use ed_state_vars, only : edtype                    & ! structure
                              , polygontype               & ! structure
                              , sitetype                  & ! structure
                              , patchtype                 ! ! structure
      use disturb_coms , only : treefall_disturbance_rate & ! intent(in)
                              , lutime                    & ! structure
                              , ianth_disturb             & ! intent(in)
                              , include_fire              & ! intent(in)
                              , plantation_rotation       & ! intent(in)
                              , mature_harvest_age        ! ! intent(in)
      use ed_max_dims  , only : n_pft                     ! ! intent(in)

      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(edtype)     , target           :: cgrid
      integer          , intent(in)       :: month
      integer          , intent(in)       :: year
      !----- Local variables. -------------------------------------------------------------!
      type(polygontype), pointer          :: cpoly
      type(sitetype)   , pointer          :: csite
      type(patchtype)  , pointer          :: cpatch
      type(lutime)     , pointer          :: clutime
      integer                             :: ipy
      integer                             :: isi
      integer                             :: ipa
      integer                             :: ico
      integer                             :: ipft
      integer                             :: iyear
      integer                             :: useyear
      logical                             :: mature_plantation
      logical                             :: mature_secondary
      real                                :: weight
      real                                :: sumweight
      real                                :: pharvest
      real                                :: fire_dist_rate
      !----- Local constants. -------------------------------------------------------------!
      real             , parameter        :: max_pharvest = 1.-2.*epsilon(1.)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !      Loop over sites and polygons.                                                 !
      !------------------------------------------------------------------------------------!
      polyloop: do ipy = 1,cgrid%npolygons
         cpoly => cgrid%polygon(ipy)

         siteloop: do isi = 1,cpoly%nsites
            csite => cpoly%site(isi)

            !------ Calculate fire disturbance rates only when fire is on. ----------------!
            select case (include_fire)
            case (0) 
               fire_dist_rate = 0.0
            case default
               fire_dist_rate = sum(cpoly%lambda_fire(1:12,isi)) / 12.0
            end select
            cpoly%fire_disturbance_rate(isi) = fire_dist_rate
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !      For natural disturbance, use largest disturbance mode.                  !
            !------------------------------------------------------------------------------!
            if (fire_dist_rate > treefall_disturbance_rate) then
               cpoly%nat_disturbance_rate(isi) = fire_dist_rate
               cpoly%nat_dist_type(isi)        = 1
            else
               cpoly%nat_disturbance_rate(isi) = treefall_disturbance_rate
               cpoly%nat_dist_type(isi)        = 0
            end if

            !----- Set disturbance rates assuming only natural disturbance. ---------------!
            cpoly%disturbance_rates(1:2,1:3,isi) = 0.0
            cpoly%disturbance_rates(3,1,isi)     = 0.0
            cpoly%disturbance_rates(3,2:3,isi)   = cpoly%nat_disturbance_rate(isi)

            !------------------------------------------------------------------------------!
            !      Now it is time for anthropogenic disturbance rates.  We no longer need  !
            ! to check the year because all years in this simulation are assigned          !
            ! prescribed disturbance rates (even if that means zero disturbance).          !
            !------------------------------------------------------------------------------!
            useyear = cpoly%num_landuse_years(isi)
            !----- Loop over years. -------------------------------------------------------!
            find_lu_year: do iyear = 1,cpoly%num_landuse_years(isi)
               
               if (year == cpoly%clutimes(iyear,isi)%landuse_year) then
                  useyear = iyear
                  exit find_lu_year
               end if
            end do find_lu_year

            !------------------------------------------------------------------------------!
            !    Update land-use transition matrix.  The matrix (and array for each        !
            ! site) has two indices.  The leftmost one is the destination, and the         !
            ! middle one is the source.  Land use types are:                               !
            !    - 1. Agriculture                                                          !
            !    - 2. Secondary forest                                                     !
            !    - 3. Primary forest                                                       !
            !------------------------------------------------------------------------------!
            clutime => cpoly%clutimes(useyear,isi)

            !----- Agriculture to agriculture (1 => 1). -----------------------------------!
            cpoly%disturbance_rates(1,1,isi) = 0.0
            !----- Secondary forest to agriculture (2 => 1). ------------------------------!
            cpoly%disturbance_rates(1,2,isi) = clutime%landuse(7) + clutime%landuse(9)
            !----- Primary forest to agriculture (3 => 1). --------------------------------!
            cpoly%disturbance_rates(1,3,isi) = clutime%landuse(4) + clutime%landuse(5)

            !------------------------------------------------------------------------------!
            !     Agriculture to secondary land (1 => 2).   Here we also account for       !
            ! transitions to primary land due to abandonment.  In ED definition, a land    !
            ! is considered primary only when the last disturbance was natural.            !
            !------------------------------------------------------------------------------!
            cpoly%disturbance_rates(2,1,isi) = clutime%landuse(8) + clutime%landuse(10)    &
                                             + clutime%landuse(3) + clutime%landuse(6)

            !------------------------------------------------------------------------------!
            !     Secondary forest to secondary forest (2 => 2).  Here we must convert the !
            ! harvest probability of being cut given that the DBH exceeds the minimum DBH. !
            ! This is done only when anthropogenic disturbance is on and we are not seek-  !
            ! ing the biomass target, otherwise we set it to zero.                         !
            !------------------------------------------------------------------------------!
            sumweight = 0.
            pharvest  = 0.
            if (ianth_disturb == 1 .and. clutime%landuse(12) < 0) then
               patchloop: do ipa=1,csite%npatches
                  cpatch => csite%patch(ipa)

                  !------------------------------------------------------------------------!
                  !    Check whether anything is going to be harvested from this patch.    !
                  !------------------------------------------------------------------------!
                  mature_plantation = csite%plantation(ipa) == 1                   .and.   &
                                      csite%age(ipa)         > plantation_rotation
                  mature_secondary  = csite%dist_type(ipa)  == 2                   .and.   &
                                      csite%plantation(ipa) /= 1                   .and.   &
                                      csite%age(ipa)         > mature_harvest_age
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !     If this patch is secondary and mature, compute the weighted aver-  !
                  ! age of the probability of harvest, using the basal area as the weight. !
                  !------------------------------------------------------------------------!
                  if (mature_plantation .or. mature_secondary) then
                     cohortloop: do ico=1,cpatch%ncohorts
                        ipft = cpatch%pft(ico)
                        if (cpatch%dbh(ico) >= cpoly%mindbh_secondary(ipft,isi)) then
                           weight    = cpatch%nplant(ico) * cpatch%basarea(ico)            &
                                     * csite%area(ipa)
                           pharvest  = pharvest                                            &
                                     + cpoly%probharv_secondary(ipft,isi) * weight
                           sumweight = sumweight + weight
                        end if  
                     end do cohortloop
                  end if
                  !----- Normalise the probability, unless it's zero. ---------------------!
                  if (sumweight > 0.) then
                     pharvest = min(max_pharvest,pharvest/sumweight)
                  else
                     pharvest = 0.
                  end if
               end do patchloop
            end if
            !----- Convert the probability into disturbance rate. -------------------------#
            cpoly%disturbance_rates(2,2,isi) = - log(1.0 - pharvest)
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !     Primary forest to secondary forest (3 => 2).  Here we must decide        !
            ! whether we want a biomass target harvesting, or a selective logging.  When   !
            ! harvesting is biomass-based, then we set the disturbance rate to zero,       !
            ! otherwise we read the disturbance rate.                                      !
            !------------------------------------------------------------------------------!
            if (clutime%landuse(12) > 0.) then
               cpoly%disturbance_rates(2,3,isi) = 0.0
            else
               cpoly%disturbance_rates(2,3,isi) = clutime%landuse(11)
            end if
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !    Agriculture to primary forest (1 => 3).  This must be zero, even though   !
            ! the landuse(3) and landuse(6) may not be (these will be added to the         !
            ! secondary instead).  This is because abandoned land is always considered     !
            ! secondary forest since the "disturbance" that creates this patch is not      !
            ! natural.                                                                     !
            !------------------------------------------------------------------------------!
            cpoly%disturbance_rates(3,1,isi) = 0.0


            !------------------------------------------------------------------------------!
            !    Secondary forest to primary forest (2 => 3).  This is the natural         !
            ! disturbance rate, because the definition of primary or secondary forest is   !
            ! defined based on the last disturbance.                                       !
            !------------------------------------------------------------------------------!
            cpoly%disturbance_rates(3,2,isi)= cpoly%nat_disturbance_rate(isi)
            !----- Primary forest to primary forest (3 => 3).  Natural disturbance. -------!
            cpoly%disturbance_rates(3,3,isi)= cpoly%nat_disturbance_rate(isi)
            !------------------------------------------------------------------------------!

         end do siteloop
      end do polyloop

      return
   end subroutine site_disturbance_rates
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine assigns initial conditions to a patch that has been disturbed     !
   !---------------------------------------------------------------------------------------!
   subroutine initialize_disturbed_patch(csite,atm_tmp,np,dp,lsl)

      use ed_state_vars, only : sitetype  & ! structure
                              , patchtype ! ! structure
      use consts_coms  , only : t3ple     ! ! intent(in)
      use grid_coms    , only : nzs       & ! intent(in)
                              , nzg       ! ! intent(in)

      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(sitetype), target      :: csite
      real          , intent(in)  :: atm_tmp
      integer       , intent(in)  :: np
      integer       , intent(in)  :: dp
      integer       , intent(in)  :: lsl
      !----- Local variables. -------------------------------------------------------------!
      integer                     :: k
      !------------------------------------------------------------------------------------!


      !----- Start with an empty patch.  Surviving cohorts will be added later. -----------!
      csite%patch(np)%ncohorts = 0

      !----- Initialise the surface water diagnostic variables with some defaults. --------!
      do k=1,nzs
         csite%sfcwater_tempk(k,np) = atm_tmp
         if (atm_tmp >= t3ple) then
            csite%sfcwater_fracliq(k,np) = 1.0
         else
            csite%sfcwater_fracliq(k,np) = 1.0
         end if
      end do

      !------------------------------------------------------------------------------------!
      !     Initialise most variables, except dist_type, plantation, and area, which will  !
      ! be defined outside this subroutine.  Most of the following variables will receive  !
      ! properties from the donor patches.                                                 !
      !------------------------------------------------------------------------------------!
      csite%age                        (np) = 0.0
      csite%fast_soil_C                (np) = 0.0
      csite%slow_soil_C                (np) = 0.0
      csite%structural_soil_C          (np) = 0.0
      csite%structural_soil_L          (np) = 0.0
      csite%mineralized_soil_N         (np) = 0.0
      csite%fast_soil_N                (np) = 0.0
      csite%sum_dgd                    (np) = 0.0
      csite%sum_chd                    (np) = 0.0
      csite%can_depth                  (np) = 0.0
      csite%can_theta                  (np) = 0.0
      csite%can_theiv                  (np) = 0.0
      csite%can_vpdef                  (np) = 0.0
      csite%can_prss                   (np) = 0.0
      csite%can_shv                    (np) = 0.0
      csite%can_co2                    (np) = 0.0
      csite%ggbare                     (np) = 0.0
      csite%ggveg                      (np) = 0.0
      csite%soil_energy          (1:nzg,np) = 0.0
      csite%soil_water           (1:nzg,np) = 0.0
      csite%sfcwater_mass        (1:nzs,np) = 0.0
      csite%sfcwater_energy      (1:nzs,np) = 0.0
      csite%sfcwater_depth       (1:nzs,np) = 0.0
      csite%rough                      (np) = 0.0
      csite%fsc_in                     (np) = 0.0
      csite%ssc_in                     (np) = 0.0
      csite%ssl_in                     (np) = 0.0
      csite%fsn_in                     (np) = 0.0
      csite%total_plant_nitrogen_uptake(np) = 0.0
      !------------------------------------------------------------------------------------!

      !----- Initialise all fast variables. -----------------------------------------------!
      call init_ed_patch_vars(csite,np,np,lsl)

      return
   end subroutine initialize_disturbed_patch
   !=======================================================================================!
   !=======================================================================================!





   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine will re-scale some patch variables using new area fraction.       !
   !---------------------------------------------------------------------------------------!
   subroutine normal_patch_vars(csite,ipa, area_fac)
      use ed_state_vars, only : sitetype  & ! structure
                              , patchtype ! ! structure
      use ed_max_dims  , only : n_pft     ! ! intent(in)
      use grid_coms    , only : nzg       ! ! intent(in)


      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(sitetype), target      :: csite
      integer       , intent(in)  :: ipa
      real          , intent(in)  :: area_fac
      !----- Local variables. -------------------------------------------------------------!
      integer                     :: k
      !------------------------------------------------------------------------------------!

      csite%fast_soil_C       (ipa) = csite%fast_soil_C       (ipa) * area_fac
      csite%slow_soil_C       (ipa) = csite%slow_soil_C       (ipa) * area_fac
      csite%structural_soil_C (ipa) = csite%structural_soil_C (ipa) * area_fac
      csite%structural_soil_L (ipa) = csite%structural_soil_L (ipa) * area_fac
      csite%mineralized_soil_N(ipa) = csite%mineralized_soil_N(ipa) * area_fac
      csite%fast_soil_N       (ipa) = csite%fast_soil_N       (ipa) * area_fac
      csite%sum_dgd           (ipa) = csite%sum_dgd           (ipa) * area_fac
      csite%sum_chd           (ipa) = csite%sum_chd           (ipa) * area_fac
      csite%can_theta         (ipa) = csite%can_theta         (ipa) * area_fac
      csite%can_theiv         (ipa) = csite%can_theiv         (ipa) * area_fac
      csite%can_vpdef         (ipa) = csite%can_vpdef         (ipa) * area_fac
      csite%can_prss          (ipa) = csite%can_prss          (ipa) * area_fac
      csite%can_shv           (ipa) = csite%can_shv           (ipa) * area_fac
      csite%can_co2           (ipa) = csite%can_co2           (ipa) * area_fac
      csite%can_depth         (ipa) = csite%can_depth         (ipa) * area_fac
      csite%ggbare            (ipa) = csite%ggbare            (ipa) * area_fac
      csite%ggveg             (ipa) = csite%ggveg             (ipa) * area_fac
      csite%rough             (ipa) = csite%rough             (ipa) * area_fac
      csite%mean_rh           (ipa) = csite%mean_rh           (ipa) * area_fac
      csite%mean_cwd_rh       (ipa) = csite%mean_cwd_rh       (ipa) * area_fac
      csite%today_A_decomp    (ipa) = csite%today_A_decomp    (ipa) * area_fac
      csite%today_Af_decomp   (ipa) = csite%today_Af_decomp   (ipa) * area_fac
      csite%fsc_in            (ipa) = csite%fsc_in            (ipa) * area_fac
      csite%ssc_in            (ipa) = csite%ssc_in            (ipa) * area_fac
      csite%ssl_in            (ipa) = csite%ssl_in            (ipa) * area_fac
      csite%fsn_in            (ipa) = csite%fsn_in            (ipa) * area_fac
      csite%total_plant_nitrogen_uptake(ipa) = csite%total_plant_nitrogen_uptake(ipa)      &
                                             * area_fac


      !----- Do the same thing for the multiple-level variables. --------------------------!
      do k=1,n_pft
         csite%repro                 (k,ipa) = csite%repro          (k,ipa) * area_fac
      end do
      do k = 1, csite%nlev_sfcwater(ipa)
         csite%sfcwater_mass         (k,ipa) = csite%sfcwater_mass  (k,ipa) * area_fac
         csite%sfcwater_energy       (k,ipa) = csite%sfcwater_energy(k,ipa) * area_fac
         csite%sfcwater_depth        (k,ipa) = csite%sfcwater_depth (k,ipa) * area_fac
      end do
      do k = 1, nzg
         csite%soil_energy           (k,ipa) = csite%soil_energy    (k,ipa) * area_fac
         csite%soil_water(k,ipa)             = csite%soil_water     (k,ipa) * area_fac
      end do
      
      csite%fast_soil_C(ipa)                 = csite%fast_soil_C(ipa)       * area_fac
      csite%structural_soil_C(ipa)           = csite%structural_soil_C(ipa) * area_fac
      csite%structural_soil_L(ipa)           = csite%structural_soil_L(ipa) * area_fac
      csite%fast_soil_N(ipa)                 = csite%fast_soil_N(ipa)       * area_fac
      return
   end subroutine normal_patch_vars
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine will include the contribution of each contributing patch (cp) to  !
   ! the new, disturbed patch (np).                                                        !
   !---------------------------------------------------------------------------------------!
   subroutine increment_patch_vars(csite,np, cp, area_fac)
      use ed_state_vars, only : sitetype  & ! structure
                              , patchtype ! ! structure
      use ed_max_dims  , only : n_pft     ! ! intent(in)
      use grid_coms    , only : nzg       ! ! intent(in)

      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(sitetype), target      :: csite
      integer       , intent(in)  :: np
      integer       , intent(in)  :: cp
      real          , intent(in)  :: area_fac
      !----- Local variables. -------------------------------------------------------------!
      integer                     :: k
      !------------------------------------------------------------------------------------!

      csite%fast_soil_C                (np) = csite%fast_soil_C                (np)        &
                                            + csite%fast_soil_C                (cp)        &
                                            * area_fac
      csite%slow_soil_C                (np) = csite%slow_soil_C                (np)        &
                                            + csite%slow_soil_C                (cp)        &
                                            * area_fac
      csite%structural_soil_C          (np) = csite%structural_soil_C          (np)        &
                                            + csite%structural_soil_C          (cp)        &
                                            * area_fac
      csite%structural_soil_L          (np) = csite%structural_soil_L          (np)        &
                                            + csite%structural_soil_L          (cp)        &
                                            * area_fac
      csite%mineralized_soil_N         (np) = csite%mineralized_soil_N         (np)        &
                                            + csite%mineralized_soil_N         (cp)        &
                                            * area_fac
      csite%fast_soil_N                (np) = csite%fast_soil_N                (np)        &
                                            + csite%fast_soil_N                (cp)        &
                                            * area_fac
      csite%sum_dgd                    (np) = csite%sum_dgd                    (np)        &
                                            + csite%sum_dgd                    (cp)        &
                                            * area_fac
      csite%sum_chd                    (np) = csite%sum_chd                    (np)        &
                                            + csite%sum_chd                    (cp)        &
                                            * area_fac
      csite%can_theta                  (np) = csite%can_theta                  (np)        &
                                            + csite%can_theta                  (cp)        &
                                            * area_fac
      csite%can_temp                   (np) = csite%can_temp                   (np)        &
                                            + csite%can_temp                   (cp)        &
                                            * area_fac
      csite%can_temp_pv                (np) = csite%can_temp_pv                (np)        &
                                            + csite%can_temp_pv                (cp)        &
                                            * area_fac
      csite%htry                       (np) = csite%htry                       (np)        &
                                            + csite%htry                       (cp)        &
                                            * area_fac
      csite%hprev                      (np) = csite%hprev                      (np)        &
                                            + csite%hprev                      (cp)        &
                                            * area_fac
      csite%can_theiv                  (np) = csite%can_theiv                  (np)        &
                                            + csite%can_theiv                  (cp)        &
                                            * area_fac
      csite%can_vpdef                  (np) = csite%can_vpdef                  (np)        &
                                            + csite%can_vpdef                  (cp)        &
                                            * area_fac
      csite%can_prss                   (np) = csite%can_prss                   (np)        &
                                            + csite%can_prss                   (cp)        &
                                            * area_fac
      csite%can_shv                    (np) = csite%can_shv                    (np)        &
                                            + csite%can_shv                    (cp)        &
                                            * area_fac
      csite%can_co2                    (np) = csite%can_co2                    (np)        &
                                            + csite%can_co2                    (cp)        &
                                            * area_fac
      csite%can_depth                  (np) = csite%can_depth                  (np)        &
                                            + csite%can_depth                  (cp)        &
                                            * area_fac
      csite%ggbare                     (np) = csite%ggbare                     (np)        &
                                            + csite%ggbare                     (cp)        &
                                            * area_fac
      csite%ggveg                      (np) = csite%ggveg                      (np)        &
                                            + csite%ggveg                      (cp)        &
                                            * area_fac
      csite%rough                      (np) = csite%rough                      (np)        &
                                            + csite%rough                      (cp)        &
                                            * area_fac
      csite%mean_rh                    (np) = csite%mean_rh                    (np)        &
                                            + csite%mean_rh                    (cp)        &
                                            * area_fac
      csite%mean_cwd_rh                (np) = csite%mean_cwd_rh                (np)        &
                                            + csite%mean_cwd_rh                (cp)        &
                                            * area_fac
      csite%today_A_decomp             (np) = csite%today_A_decomp             (np)        &
                                            + csite%today_A_decomp             (cp)        &
                                            * area_fac
      csite%today_Af_decomp            (np) = csite%today_Af_decomp            (np)        &
                                            + csite%today_Af_decomp            (cp)        &
                                            * area_fac
      csite%fsc_in                     (np) = csite%fsc_in                     (np)        &
                                            + csite%fsc_in                     (cp)        &
                                            * area_fac
      csite%ssc_in                     (np) = csite%ssc_in                     (np)        &
                                            + csite%ssc_in                     (cp)        &
                                            * area_fac
      csite%ssl_in                     (np) = csite%ssl_in                     (np)        &
                                            + csite%ssl_in                     (cp)        &
                                            * area_fac
      csite%fsn_in                     (np) = csite%fsn_in                     (np)        &
                                            + csite%fsn_in                     (cp)        &
                                            * area_fac
      csite%total_plant_nitrogen_uptake(np) = csite%total_plant_nitrogen_uptake(np)        &
                                            + csite%total_plant_nitrogen_uptake(cp)        &
                                            * area_fac

      !----- Do the same thing for the multiple-level variables. --------------------------!
      do k=1,n_pft
         csite%repro                 (k,np) = csite%repro                    (k,np)        &
                                            + csite%repro                    (k,cp)        &
                                            * area_fac
      end do
      do k = 1, csite%nlev_sfcwater(cp)
         csite%sfcwater_mass         (k,np) = csite%sfcwater_mass            (k,np)        &
                                            + csite%sfcwater_mass            (k,cp)        &
                                            * area_fac
         csite%sfcwater_energy       (k,np) = csite%sfcwater_energy          (k,np)        &
                                            + csite%sfcwater_energy          (k,cp)        &
                                            * csite%sfcwater_mass            (k,cp)        &
                                            * area_fac
         csite%sfcwater_depth        (k,np) = csite%sfcwater_depth           (k,np)        &
                                            + csite%sfcwater_depth           (k,cp)        &
                                            * area_fac
      end do
      do k = 1, nzg
         csite%soil_energy           (k,np) = csite%soil_energy              (k,np)        &
                                            + csite%soil_energy              (k,cp)        &
                                            * area_fac
         csite%soil_water(k,np)             = csite%soil_water               (k,np)        &
                                            + csite%soil_water               (k,cp)        &
                                            * area_fac
      end do
      return
   end subroutine increment_patch_vars
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine will populate the disturbed patch with the cohorts that were      !
   ! disturbed but did not go extinct.                                                     !
   !---------------------------------------------------------------------------------------!
   subroutine insert_survivors(csite, np, cp, q, area_fac,poly_dest_type,mindbh_harvest)

      use ed_state_vars, only : sitetype     & ! structure
                              , patchtype    ! ! structure
      use ed_misc_coms , only : idoutput     & ! intent(in)
                              , iqoutput     & ! intent(in)
                              , imoutput     ! ! intent(in)
      use ed_max_dims  , only : n_pft        ! ! intent(in)
      use mortality    , only : survivorship ! ! function
    
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(sitetype)                  , target      :: csite
      integer                         , intent(in)  :: q
      integer                         , intent(in)  :: poly_dest_type
      integer                         , intent(in)  :: np
      integer                         , intent(in)  :: cp
      real          , dimension(n_pft), intent(in)  :: mindbh_harvest
      real                            , intent(in)  :: area_fac
      !----- Local variables. -------------------------------------------------------------!
      type(patchtype)                 , pointer     :: cpatch
      type(patchtype)                 , pointer     :: npatch
      type(patchtype)                 , pointer     :: tpatch
      logical        , dimension(:)   , allocatable :: mask
      integer                                       :: ico
      integer                                       :: nco
      real                                          :: n_survivors
      real                                          :: survival_fac
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      ! cpatch => Old patches we are looping                                               !
      ! npatch => New transition patch                                                     !
      ! tpatch => temporary patch                                                          !
      !------------------------------------------------------------------------------------!
      cpatch => csite%patch(cp)
      npatch => csite%patch(np)
      nullify(tpatch)
      allocate(tpatch)

      !----- Mask: flag to decide whether the cohort survived or not. ---------------------!
      if (cpatch%ncohorts > 0) then
         allocate(mask(cpatch%ncohorts))
         mask(:) = .false.
    
         survivalloop: do ico = 1,cpatch%ncohorts
            survival_fac = survivorship(q,poly_dest_type, mindbh_harvest, csite, cp, ico)  &
                         * area_fac
            n_survivors     = cpatch%nplant(ico) * survival_fac

            !----- If something survived, make a new cohort. ------------------------------!
            mask(ico) = n_survivors > 0.0
         end do survivalloop
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     If the new patch has received survivors from a donor already, then it should   !
      ! have cohorts.  So the temporary patch vector will be the sum of the new cohorts    !
      ! found here, plus those already applied previously in the loop calling this sub-    !
      ! routine.                                                                           !
      !------------------------------------------------------------------------------------!
      if (npatch%ncohorts > 0) then
         nco = npatch%ncohorts
         call allocate_patchtype(tpatch,count(mask) + npatch%ncohorts)
         call copy_patchtype(npatch,tpatch,1,npatch%ncohorts,1,npatch%ncohorts)
         call deallocate_patchtype(npatch)
      else
         nco = 0
         call allocate_patchtype(tpatch,count(mask))
      end if


      cohortloop: do ico = 1,cpatch%ncohorts
         
         survival_fac = survivorship(q,poly_dest_type, mindbh_harvest, csite, cp, ico)     &
                      * area_fac
         n_survivors  = cpatch%nplant(ico) * survival_fac

         !----- If mask is true, at least some of this cohort survived. -------------------!
         if (mask(ico)) then
            nco = nco + 1
            call copy_patchtype(cpatch,tpatch,ico,ico,nco,nco)

            !------------------------------------------------------------------------------!
            !    Scale the total area based on the new population density and new area.    !
            ! We must also rescale all extensive properties from cohorts, since they are   !
            ! per unit area and we are effectively changing the area.                      !
            ! IMPORTANT: Only cohort-level variables that have units per area (m2) should  !
            !            be rescaled.  Variables whose units are per plant should _NOT_ be !
            !            included here.                                                    !
            !------------------------------------------------------------------------------!
            tpatch%lai                (nco) = tpatch%lai               (nco) * survival_fac
            tpatch%wai                (nco) = tpatch%wai               (nco) * survival_fac
            tpatch%nplant             (nco) = tpatch%nplant            (nco) * survival_fac
            tpatch%mean_gpp           (nco) = tpatch%mean_gpp          (nco) * survival_fac
            tpatch%mean_leaf_resp     (nco) = tpatch%mean_leaf_resp    (nco) * survival_fac
            tpatch%mean_root_resp     (nco) = tpatch%mean_root_resp    (nco) * survival_fac
            tpatch%mean_growth_resp   (nco) = tpatch%mean_growth_resp  (nco) * survival_fac
            tpatch%mean_storage_resp  (nco) = tpatch%mean_storage_resp (nco) * survival_fac
            tpatch%mean_vleaf_resp    (nco) = tpatch%mean_vleaf_resp   (nco) * survival_fac
            tpatch%today_gpp          (nco) = tpatch%today_gpp         (nco) * survival_fac
            tpatch%today_nppleaf      (nco) = tpatch%today_nppleaf     (nco) * survival_fac
            tpatch%today_nppfroot     (nco) = tpatch%today_nppfroot    (nco) * survival_fac
            tpatch%today_nppsapwood   (nco) = tpatch%today_nppsapwood  (nco) * survival_fac
            tpatch%today_nppcroot     (nco) = tpatch%today_nppcroot    (nco) * survival_fac
            tpatch%today_nppseeds     (nco) = tpatch%today_nppseeds    (nco) * survival_fac
            tpatch%today_nppwood      (nco) = tpatch%today_nppwood     (nco) * survival_fac
            tpatch%today_nppdaily     (nco) = tpatch%today_nppdaily    (nco) * survival_fac
            tpatch%today_gpp_pot      (nco) = tpatch%today_gpp_pot     (nco) * survival_fac
            tpatch%today_gpp_lightmax (nco) = tpatch%today_gpp_lightmax(nco) * survival_fac
            tpatch%today_gpp_moistmax (nco) = tpatch%today_gpp_moistmax(nco) * survival_fac
            tpatch%today_leaf_resp    (nco) = tpatch%today_leaf_resp   (nco) * survival_fac
            tpatch%today_root_resp    (nco) = tpatch%today_root_resp   (nco) * survival_fac
            tpatch%gpp                (nco) = tpatch%gpp               (nco) * survival_fac
            tpatch%leaf_respiration   (nco) = tpatch%leaf_respiration  (nco) * survival_fac
            tpatch%root_respiration   (nco) = tpatch%root_respiration  (nco) * survival_fac
            tpatch%monthly_dndt       (nco) = tpatch%monthly_dndt      (nco) * survival_fac
            tpatch%leaf_water         (nco) = tpatch%leaf_water        (nco) * survival_fac
            tpatch%leaf_hcap          (nco) = tpatch%leaf_hcap         (nco) * survival_fac
            tpatch%leaf_energy        (nco) = tpatch%leaf_energy       (nco) * survival_fac
            tpatch%wood_water         (nco) = tpatch%wood_water        (nco) * survival_fac
            tpatch%wood_hcap          (nco) = tpatch%wood_hcap         (nco) * survival_fac
            tpatch%wood_energy        (nco) = tpatch%wood_energy       (nco) * survival_fac
            !----- Crown area shall not exceed 1. -----------------------------------------!
            tpatch%crown_area         (nco) = min(1.,tpatch%crown_area(nco) * survival_fac)
            !----- Carbon flux monthly means are extensive, we must convert them. ---------!
            if (idoutput > 0 .or. imoutput > 0 .or. iqoutput > 0) then
               tpatch%dmean_par_l     (nco) = tpatch%dmean_par_l      (nco) * survival_fac
               tpatch%dmean_par_l_beam(nco) = tpatch%dmean_par_l_beam (nco) * survival_fac
               tpatch%dmean_par_l_diff(nco) = tpatch%dmean_par_l_diff (nco) * survival_fac
            end if
            if (imoutput > 0 .or. iqoutput > 0) then
               tpatch%mmean_par_l     (nco) = tpatch%mmean_par_l      (nco) * survival_fac
               tpatch%mmean_par_l_beam(nco) = tpatch%mmean_par_l_beam (nco) * survival_fac
               tpatch%mmean_par_l_diff(nco) = tpatch%mmean_par_l_diff (nco) * survival_fac
            end if
            if (iqoutput > 0) then
               tpatch%qmean_par_l     (:,nco) = tpatch%qmean_par_l      (:,nco)            &
                                              * survival_fac
               tpatch%qmean_par_l_beam(:,nco) = tpatch%qmean_par_l_beam (:,nco)            &
                                              * survival_fac
               tpatch%qmean_par_l_diff(:,nco) = tpatch%qmean_par_l_diff (:,nco)            &
                                              * survival_fac
            end if
            
            !----- Make mortality rate due to disturbance zero to avoid double counting. --!
            tpatch%mort_rate(5,nco) = 0.0
         end if
      end do cohortloop

      !----- Copy the temporary patch into the newpatch. ----------------------------------!
      call allocate_patchtype(npatch,tpatch%ncohorts)
      call copy_patchtype(tpatch,npatch,1,tpatch%ncohorts,1,tpatch%ncohorts)
      call deallocate_patchtype(tpatch)

      deallocate(tpatch)
      if (allocated(mask)) deallocate(mask)


      return
   end subroutine insert_survivors
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine updates the litter pools after a disturbance takes place.         !
   !---------------------------------------------------------------------------------------!
   subroutine accum_dist_litt(csite,np,cp,q,area_fac,poly_dest_type,mindbh_harvest)
      use ed_state_vars, only : sitetype     & ! structure
                              , patchtype    & ! structure
                              , polygontype  ! ! structure
      use decomp_coms  , only : f_labile     ! ! intent(in)
      use ed_max_dims  , only : n_pft        ! ! intent(in)
      use pft_coms     , only : c2n_storage  & ! intent(in)
                              , c2n_leaf     & ! intent(in)
                              , c2n_recruit  & ! intent(in)
                              , c2n_stem     & ! intent(in)
                              , l2n_stem     ! ! intent(in)
      use pft_coms     , only : agf_bs       ! ! intent(in)
      use grid_coms    , only : nzg          ! ! intent(in)
      use mortality    , only : survivorship ! ! function

      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(sitetype)                   , target     :: csite
      integer                          , intent(in) :: np
      integer                          , intent(in) :: cp
      real           , dimension(n_pft), intent(in) :: mindbh_harvest
      integer                          , intent(in) :: q
      real                             , intent(in) :: area_fac
      integer                          , intent(in) :: poly_dest_type
      !----- Local variables. -------------------------------------------------------------!
      type(patchtype)                  , pointer    :: cpatch
      type(patchtype)                  , pointer    :: npatch
      integer                                       :: ico
      integer                                       :: ipft
      real                                          :: loss_fraction
      real                                          :: fast_litter
      real                                          :: struct_litter
      real                                          :: struct_lignin
      real                                          :: fast_litter_n
      real                                          :: struct_cohort
      real                                          :: survival_fac
      !------------------------------------------------------------------------------------!

      !---- Initialise the non-scaled litter pools. ---------------------------------------!
      fast_litter   = 0.0
      struct_litter = 0.0
      struct_lignin = 0.0
      fast_litter_n = 0.0

      !------------------------------------------------------------------------------------!
      ! cpatch => contributing patch                                                       !
      ! npatch => new patch.                                                               !
      !------------------------------------------------------------------------------------!
      cpatch => csite%patch(cp)
      npatch => csite%patch(np)

      do ico = 1,cpatch%ncohorts
         ipft = cpatch%pft(ico)

         !---------------------------------------------------------------------------------!
         !     Find the loss fraction, which normally corresponds to the above-ground bio- !
         ! mass in case the patch was harvest/logged, or nothing in case it was a natural  !
         ! disturbance.                                                                    !
         !---------------------------------------------------------------------------------!
         select case(q)
         case (1,2)
            loss_fraction = agf_bs(ipft)
         case (3)
            loss_fraction = 0.
         end select
         !---------------------------------------------------------------------------------!

         !----- Find survivorship only once. ----------------------------------------------!
         survival_fac  = survivorship(q,poly_dest_type,mindbh_harvest,csite,cp,ico)


         fast_litter   = fast_litter                                                       &
                       + (1. - survival_fac)                                               &
                       * ( f_labile(ipft) * cpatch%balive(ico) + cpatch%bstorage(ico))     &
                       * cpatch%nplant(ico)
         fast_litter_n = fast_litter_n                                                     &
                       + (1. - survival_fac)                                               &
                       * ( f_labile(ipft) * cpatch%balive(ico) / c2n_leaf(ipft)            &
                         + cpatch%bstorage(ico) / c2n_storage )                            &
                       * cpatch%nplant(ico)

         struct_cohort = cpatch%nplant(ico)                                                &
                       * (1. - survival_fac)                                               &
                       * ( (1. - loss_fraction ) * cpatch%bdead(ico)                       &
                         + (1. - f_labile(ipft)) * cpatch%balive(ico) )

         struct_litter = struct_litter + struct_cohort
         struct_lignin = struct_lignin + struct_cohort * l2n_stem / c2n_stem(ipft)
      end do

      !----- Load disturbance litter directly into carbon and N pools. --------------------!
      csite%fast_soil_C(np)       = csite%fast_soil_C(np)       + fast_litter   * area_fac
      csite%structural_soil_C(np) = csite%structural_soil_C(np) + struct_litter * area_fac
      csite%structural_soil_L(np) = csite%structural_soil_L(np) + struct_lignin * area_fac
      csite%fast_soil_N(np)       = csite%fast_soil_N(np)       + fast_litter_n * area_fac

      return
   end subroutine accum_dist_litt
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    Add a cohort of the appropriate PFT type to populate a plantation/cropland/pasture !
   ! patch.                                                                                !
   !---------------------------------------------------------------------------------------!
   subroutine plant_patch(csite,np,mzg,pft,density,ntext_soil,green_leaf_factor            &
                         ,height_factor,lsl)
      use ed_state_vars , only  : sitetype                 & ! structure
                                , patchtype                ! ! structure
      use pft_coms       , only : q                        & ! intent(in)
                                , qsw                      & ! intent(in)
                                , sla                      & ! intent(in)
                                , hgt_min                  & ! intent(in)
                                , hgt_max                  & ! intent(in)
                                , dbh_bigleaf              ! ! intent(in)
      use ed_misc_coms   , only : dtlsm                    & ! intent(in)
                                , ibigleaf                 ! ! intent(in)
      use fuse_fiss_utils, only : sort_cohorts             ! ! sub-routine
      use ed_therm_lib   , only : calc_veg_hcap            ! ! function
      use consts_coms    , only : t3ple                    & ! intent(in)
                                , pio4                     ! ! intent(in)
      use allometry      , only : h2dbh                    & ! function
                                , dbh2bd                   & ! function
                                , area_indices             & ! function
                                , ed_biomass               ! ! function
      use ed_max_dims    , only : n_pft                    ! ! intent(in)
      use phenology_coms , only : retained_carbon_fraction ! ! intent(in)
      use phenology_aux  , only : pheninit_balive_bstorage ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(sitetype)                  , target     :: csite
      integer                         , intent(in) :: mzg
      integer                         , intent(in) :: np
      integer                         , intent(in) :: pft
      integer                         , intent(in) :: lsl
      integer       , dimension(mzg)  , intent(in) :: ntext_soil
      real          , dimension(n_pft), intent(in) :: green_leaf_factor
      real                            , intent(in) :: density
      real                            , intent(in) :: height_factor
      !----- Local variables. -------------------------------------------------------------!
      type(patchtype)                 , pointer    :: cpatch
      type(patchtype)                 , pointer    :: tpatch
      integer                                      :: nc
      real                                         :: salloc
      real                                         :: salloci
      real                                         :: bleaf_max
      real                                         :: balive_max
      !------------------------------------------------------------------------------------!



      cpatch => csite%patch(np)


      !------------------------------------------------------------------------------------!
      !      Reallocate the current cohort with an extra space for the planted plantation  !
      ! cohort.  If the patch was previously empty, we simply create the first cohort.     !
      !------------------------------------------------------------------------------------!
      nc = cpatch%ncohorts + 1
      if (cpatch%ncohorts > 0) then
         nullify(tpatch)
         allocate(tpatch)
         call allocate_patchtype(tpatch,cpatch%ncohorts)
         call copy_patchtype(cpatch,tpatch,1,cpatch%ncohorts,1,cpatch%ncohorts)
         call deallocate_patchtype(cpatch)
         call allocate_patchtype(cpatch,tpatch%ncohorts + 1)
         call copy_patchtype(tpatch,cpatch,1,tpatch%ncohorts,1,tpatch%ncohorts)
         call deallocate_patchtype(tpatch)
         deallocate(tpatch)
      else
         call allocate_patchtype(cpatch,1)
      end if
         
      cpatch%ncohorts = nc
      csite%paco_n(np)= nc

      !------------------------------------------------------------------------------------!
      !      We can add all the biomass needed in a single cohort, even if it exceeds the  !
      ! safe range (i.e., with maximum LAI exceeding the unity).  In case it is a large    !
      ! cohort, it will soon be split at the splitting call of apply_disturbances().  This !
      ! new cohort will always be the last here, they will be sorted afterwards too.       !
      !------------------------------------------------------------------------------------!
      cpatch%pft(nc)    = pft
      cpatch%nplant(nc) = density
      select case (ibigleaf)
      case (0)
         !---------------------------------------------------------------------------------!
         !    SAS approximation, assign height and use it to find DBH and the structural   !
         ! (dead) biomass.                                                                 !
         !---------------------------------------------------------------------------------!
         cpatch%hite (nc) = hgt_min(cpatch%pft(nc)) * min(1.0,height_factor)
         cpatch%dbh  (nc) = h2dbh(cpatch%hite(nc),cpatch%pft(nc))
         cpatch%bdead(nc) = dbh2bd(cpatch%dbh(nc),cpatch%pft(nc))
         !---------------------------------------------------------------------------------!

      case (1)
         !---------------------------------------------------------------------------------!
         !    Big leaf approximation, assign the typical DBH and height and use them to    !
         ! find height and the structural (dead) biomass.                                  !
         !---------------------------------------------------------------------------------!
         cpatch%hite (nc) = hgt_max(cpatch%pft(nc))
         cpatch%dbh  (nc) = dbh_bigleaf(cpatch%pft(nc))
         cpatch%bdead(nc) = dbh2bd(cpatch%dbh(nc),cpatch%pft(nc))
         !---------------------------------------------------------------------------------!
      end select
      !------------------------------------------------------------------------------------!




      !----- Initialise other cohort-level variables. -------------------------------------!
      call init_ed_cohort_vars(cpatch, nc, lsl)
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !      Initialise the active and storage biomass scaled by the leaf drought          !
      ! phenology (or start with 1.0 if the plant doesn't shed their leaves due to water   !
      ! stress.                                                                            !
      !------------------------------------------------------------------------------------!
      call pheninit_balive_bstorage(mzg,cpatch%pft(nc),cpatch%krdepth(nc),cpatch%hite(nc)  &
                                   ,cpatch%dbh(nc),csite%soil_water(:,np),ntext_soil       &
                                   ,green_leaf_factor,cpatch%paw_avg(nc),cpatch%elongf(nc) &
                                   ,cpatch%phenology_status(nc),cpatch%bleaf(nc)           &
                                   ,cpatch%broot(nc),cpatch%bsapwooda(nc)                  &
                                   ,cpatch%bsapwoodb(nc),cpatch%balive(nc)                 &
                                   ,cpatch%bstorage(nc))
      !------------------------------------------------------------------------------------!



      !----- Compute all area indices needed. ---------------------------------------------!
      call area_indices(cpatch%nplant(nc),cpatch%bleaf(nc),cpatch%bdead(nc)                &
                       ,cpatch%balive(nc),cpatch%dbh(nc),cpatch%hite(nc),cpatch%pft(nc)    &
                       ,cpatch%sla(nc),cpatch%lai(nc),cpatch%wai(nc),cpatch%crown_area(nc) &
                       ,cpatch%bsapwooda(nc))


      !----- Find the new basal area and above-ground biomass. ----------------------------!
      cpatch%basarea(nc)= pio4 * cpatch%dbh(nc) * cpatch%dbh(nc)
      cpatch%agb(nc)    = ed_biomass(cpatch%bdead(nc),cpatch%bleaf(nc)                     &
                                    ,cpatch%bsapwooda(nc),cpatch%pft(nc))

      cpatch%leaf_temp    (nc) = csite%can_temp  (np)
      cpatch%leaf_temp_pv (nc) = csite%can_temp  (np)
      cpatch%leaf_water   (nc) = 0.0
      cpatch%leaf_vpdef   (nc) = csite%can_vpdef (np)
      cpatch%leaf_fliq    (nc) = 0.0
      cpatch%wood_temp    (nc) = csite%can_temp  (np)
      cpatch%wood_temp_pv (nc) = csite%can_temp  (np)
      cpatch%wood_water   (nc) = 0.0
      cpatch%wood_fliq    (nc) = 0.0

      !----- Because we assigned no water, the internal energy is simply hcap*T. ----------!
      call calc_veg_hcap(cpatch%bleaf(nc),cpatch%bdead(nc),cpatch%bsapwooda(nc)            &
                        ,cpatch%nplant(nc),cpatch%pft(nc)                                  &
                        ,cpatch%leaf_hcap(nc),cpatch%wood_hcap(nc))

      cpatch%leaf_energy(nc) = cpatch%leaf_hcap(nc) * cpatch%leaf_temp(nc)
      cpatch%wood_energy(nc) = cpatch%wood_hcap(nc) * cpatch%wood_temp(nc)

      call is_resolvable(csite,np,nc)

      !----- Should plantations be considered recruits? -----------------------------------!
      cpatch%new_recruit_flag(nc) = 1

      !----- Sort the cohorts so that the new cohort is at the correct height bin. --------!
      call sort_cohorts(cpatch)

      return
   end subroutine plant_patch
   !=======================================================================================!
   !=======================================================================================!
end module disturbance_utils
!==========================================================================================!
!==========================================================================================!
