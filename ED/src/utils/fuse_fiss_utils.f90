module fuse_fiss_utils_ar

  use ed_state_vars,only : copy_patchtype,deallocate_patchtype,allocate_patchtype, &
       allocate_sitetype,deallocate_sitetype,copy_sitetype_mask,copy_sitetype, &
       copy_patchtype_mask

   contains

   !===================================================
   subroutine sort_cohorts_ar(cpatch)

     use ed_state_vars,only : patchtype,patchswap_g

     implicit none
     type(patchtype),target :: cpatch
     type(patchtype),pointer :: temppatch
     integer :: iico
     integer,dimension(1):: tallid
     
     allocate(temppatch)
     call allocate_patchtype(temppatch,cpatch%ncohorts)
     
     iico = 1
     do while(iico.le.cpatch%ncohorts)
        
        tallid = maxloc(cpatch%hite)

        call copy_patchtype(cpatch,temppatch,tallid(1),tallid(1),iico,iico)
        cpatch%hite(tallid) = -1.0
        iico = iico + 1
     enddo

     call copy_patchtype(temppatch,cpatch,1,cpatch%ncohorts,1,cpatch%ncohorts)

     call deallocate_patchtype(temppatch)
     deallocate(temppatch)

     return

   end subroutine sort_cohorts_ar

   !---------------------------------------------------------
   subroutine terminate_cohorts_ar(csite,ipa)

     use fusion_fission_coms, only : min_recruit_size
     use pft_coms, only: l2n_stem,c2n_stem,c2n_storage, c2n_leaf
     use decomp_coms, only : f_labile
     use ed_state_vars,only : patchtype,sitetype
     implicit none
     
     type(sitetype),target   :: csite
     type(patchtype),pointer :: cpatch,temppatch
     logical,allocatable :: remain_table(:)
     integer :: ico,ipa
     real :: csize

     cpatch => csite%patch(ipa)

     allocate(temppatch)

     allocate(remain_table(cpatch%ncohorts))
     remain_table = .true.
     
     do ico = 1,cpatch%ncohorts
     
        csize = cpatch%nplant(ico) * (cpatch%balive(ico) + cpatch%bdead(ico) + cpatch%bstorage(ico))
        
        if(csize < (0.1 * min_recruit_size) )then
           
           remain_table(ico) = .false.

           ! Update litter pools
           
           csite%fsc_in(ipa) = csite%fsc_in(ipa) + cpatch%nplant(ico) * &
                (f_labile(cpatch%pft(ico)) * cpatch%balive(ico) + cpatch%bstorage(ico))                                           
           
           csite%fsn_in(ipa) = csite%fsn_in(ipa) + cpatch%nplant(ico) * &
                (f_labile(cpatch%pft(ico)) * cpatch%balive(ico) /       &
                c2n_leaf(cpatch%pft(ico)) + cpatch%bstorage(ico) / c2n_storage)                      
           
           csite%ssc_in(ipa) = csite%ssc_in(ipa) + cpatch%nplant(ico) * &
                ((1.0 - f_labile(cpatch%pft(ico))) * cpatch%balive(ico) + cpatch%bdead(ico))
           
           csite%ssl_in(ipa) = csite%ssl_in(ipa) + cpatch%nplant(ico) * &
                ( (1.0 - f_labile(cpatch%pft(ico))) *  &     
                cpatch%balive(ico) + cpatch%bdead(ico) ) * l2n_stem / c2n_stem                       
           
        endif
        
     enddo

     ! Reallocate the patch via the remain table

     call allocate_patchtype(temppatch,count(remain_table))
     
     call copy_patchtype_mask(cpatch,temppatch,remain_table,size(remain_table),count(remain_table))
     
     call deallocate_patchtype(cpatch)
     
     call allocate_patchtype(cpatch,count(remain_table))
     
     call copy_patchtype(temppatch,cpatch,1,cpatch%ncohorts,1,cpatch%ncohorts)
     
     ! Remove the temporary patch
     
     call deallocate_patchtype(temppatch)
     deallocate(temppatch)
     deallocate(remain_table)
     csite%cohort_count(ipa) = cpatch%ncohorts

  return
end subroutine terminate_cohorts_ar
!====================================================================

   subroutine terminate_patches_ar(csite)

     use ed_state_vars,only:polygontype,sitetype,patchtype
     use disturb_coms, only : min_new_patch_area
     
     implicit none
     
     type(sitetype),target   :: csite
     type(sitetype),pointer  :: tsite
     type(patchtype),pointer :: cpatch
     integer :: ipa
     integer,dimension(csite%npatches) :: mask
     real :: epsilon
     

     ! Loop through all the patches in this site and determine
     ! which of these patches is too small in area to be valid
     ! Remove these patches via the mask function. Realocate a 
     ! new site with only the valid patches, and normalize
     ! their areas and plant densities to reflect the area loss

     mask = 1
     epsilon = 0.0
     
     do ipa = 1,csite%npatches
        if(csite%area(ipa) .lt. min_new_patch_area .or. csite%patch(ipa)%ncohorts == 0) then
           epsilon = epsilon + csite%area(ipa)
           mask(ipa) = 0
        endif
        
     enddo

     ! use the mask to resize the patch vectors in the current site
     ! ------------------------------------------------------------
     allocate(tsite)
     call allocate_sitetype(tsite,sum(mask))
     
     call copy_sitetype_mask(csite,tsite,mask,size(mask),sum(mask))

     call deallocate_sitetype(csite)

     call allocate_sitetype(csite,sum(mask))

     mask = 0
     mask(1:tsite%npatches) = 1
     call copy_sitetype_mask(tsite,csite,mask(1:tsite%npatches),sum(mask),sum(mask))

     call deallocate_sitetype(tsite)
     deallocate(tsite)
     ! ------------------------------------------------------------


     ! Renormalize the number of plants and the total area to reflect
     ! the deletions

     do ipa = 1,csite%npatches
        csite%area(:) = csite%area(:) / (1.0-epsilon)
        cpatch => csite%patch(ipa)
        cpatch%nplant(:) = cpatch%nplant(:) / (1.0-epsilon)
     enddo

     return
   end subroutine terminate_patches_ar


   !----------------------------------------------------------

   subroutine fuse_cohorts_ar(csite,ipa, green_leaf_factor, lsl)

     use ed_state_vars,only:sitetype,patchtype
     use pft_coms            , only: rho, b1Ht, max_dbh, sla,hgt_ref
     use fusion_fission_coms , only: fusetol_h, fusetol, lai_fuse_tol
     use max_dims            , only: n_pft
     use mem_sites           , only: maxcohort

     implicit none

     type(sitetype),target :: csite
     type(patchtype),pointer :: cpatch
     integer :: ipa,ico1,ico2
     real         , dimension(n_pft) , intent(in) :: green_leaf_factor
     integer      ,                    intent(in) :: lsl

     type(patchtype),pointer :: temppatch

     logical :: fusion_test
     real    :: hite_threshold
     real    :: newn
     real    :: total_lai
     real :: tolerance_mult
     logical , allocatable, dimension(:) :: fuse_table
     real    , external :: dbh2h
     real    , external :: dbh2bl
     real, parameter :: tolerance_max = 2.5
     integer :: ncohorts_old

     if(csite%cohort_count(ipa) == 0)return ! return if there aren't any cohorts
     
     tolerance_mult = 1.0

     cpatch => csite%patch(ipa)

     allocate(fuse_table(cpatch%ncohorts))
     fuse_table(:) = .true.

     force_fusion: do
     ncohorts_old =  count(fuse_table)
     donloop:do ico1 = 1,cpatch%ncohorts-1
        if (.not. fuse_table(ico1)) cycle donloop !This one is gone
         
        recloop: do ico2 = ico1+1,cpatch%ncohorts
           if (.not. fuse_table(ico2)) cycle recloop !This one is gone

           ! get fusion height threshold
           if(rho(cpatch%pft(ico1)) == 0.0)then
              hite_threshold = b1Ht(cpatch%pft(ico1)) + hgt_ref(cpatch%pft(ico1))
           else
              hite_threshold = dbh2h(cpatch%pft(ico1), max_dbh(cpatch%pft(ico1)))
           end if

           ! test for similarity
           if(cpatch%hite(ico1) < (0.95 * hite_threshold ))then
              fusion_test = (abs(cpatch%hite(ico1) - cpatch%hite(ico2)) < fusetol_h * tolerance_mult)
           else
              fusion_test = (abs(cpatch%dbh(ico1) - cpatch%dbh(ico2)) /   &
                   (0.5*(cpatch%dbh(ico1) + cpatch%dbh(ico2))) < fusetol * tolerance_mult)
           end if

           if(fusion_test)then

           ! Cohorts have a similar size
              
              newn = cpatch%nplant(ico1) + cpatch%nplant(ico2)

              total_lai = ( cpatch%nplant(ico2) * dbh2bl(cpatch%dbh(ico2),cpatch%pft(ico2))  &
                   + cpatch%nplant(ico1) * dbh2bl(cpatch%dbh(ico1),cpatch%pft(ico1))) * sla(cpatch%pft(ico2))

              if(  cpatch%pft(ico1)      ==     cpatch%pft(ico2)     .and. &  ! cohorts are the same PFT and
                   total_lai             <       lai_fuse_tol*tolerance_mult     .and. &  ! LAI won't be too big and
                   cpatch%first_census(ico1)     == cpatch%first_census(ico2)     .and. &  ! won't mess
                   cpatch%new_recruit_flag(ico1) == cpatch%new_recruit_flag(ico2) )then    ! up output

                 call fuse_2_cohorts_ar(cpatch,ico1,ico2, newn,green_leaf_factor(cpatch%pft(ico1)), lsl)

                 fuse_table(ico1) = .false.
                 cycle donloop
              end if
           end if
        end do recloop
     end do donloop

     
     if( count(fuse_table) <= maxcohort)exit force_fusion
     if( (count(fuse_table) == ncohorts_old) .and. (tolerance_mult > tolerance_max) ) exit force_fusion
     tolerance_mult = tolerance_mult * 1.01
     
  enddo force_fusion

  ! If fusion didn't happen at all, return
  if (all(fuse_table)) then
     deallocate(fuse_table)
     return
  endif


  ! Now copy the merged patch to a temporary patch using the fuse_table
  ! as a mask.  Then deallocate the patch, copy the values from the temp
  ! and deallocate the temp.

  allocate(temppatch)
  call allocate_patchtype(temppatch,cpatch%ncohorts)

  call copy_patchtype_mask(cpatch,temppatch,fuse_table,size(fuse_table),count(fuse_table))

  call deallocate_patchtype(cpatch)
  
  call allocate_patchtype(cpatch,count(fuse_table))

  fuse_table = .false.
  fuse_table(1:cpatch%ncohorts) = .true.

  call copy_patchtype_mask(temppatch,cpatch,fuse_table,size(fuse_table),count(fuse_table))

  call deallocate_patchtype(temppatch)
  deallocate(temppatch)  
  
  call sort_cohorts_ar(cpatch)
  csite%cohort_count(ipa) = count(fuse_table)

  deallocate(fuse_table)

  return
end subroutine fuse_cohorts_ar

!=============================================================

subroutine split_cohorts_ar(cpatch, green_leaf_factor, lsl)

     use ed_state_vars,only : patchtype
     use pft_coms             , only: q, qsw, sla
     use fusion_fission_coms  , only: lai_tol
     use max_dims             , only: n_pft
     use therm_lib            , only: update_veg_energy_ct

     implicit none
     
     type(patchtype),target :: cpatch
     type(patchtype),pointer :: temppatch
     
     integer :: ico,inew
     real        , dimension(n_pft), intent(in) :: green_leaf_factor

     real         , parameter :: epsilon=0.0001
     integer                  :: ncohorts_new
     real    :: slai

     real, external :: dbh2h
     real, external :: bd2dbh
     integer, intent(in) :: lsl
     integer,allocatable :: split_mask(:)

     allocate(temppatch)
     allocate(split_mask(cpatch%ncohorts))

     split_mask = 0

     do ico = 1,cpatch%ncohorts

        slai =  cpatch%nplant(ico) * cpatch%balive(ico) * green_leaf_factor(cpatch%pft(ico))   &
             /(1.0 + q(cpatch%pft(ico)) + qsw(cpatch%pft(ico)) * &
             cpatch%hite(ico)) * sla(cpatch%pft(ico))

        if(slai > lai_tol)then
           ! Determine the split list
           split_mask(ico) = 1

        endif
     enddo

     ncohorts_new = cpatch%ncohorts + sum(split_mask)

     ! Allocate the temppatch
     
     call allocate_patchtype(temppatch,cpatch%ncohorts)
     
     ! Fill the temp space with the current patches
     
     call copy_patchtype(cpatch,temppatch,1,cpatch%ncohorts,1,cpatch%ncohorts)
     
     ! Deallocate the current patch
     
     call deallocate_patchtype(cpatch)
     
     ! Re-allocate the current patch
     
     call allocate_patchtype(cpatch,ncohorts_new)
     
     ! Transfer the temp values back in
     
     call copy_patchtype(temppatch,cpatch,1,temppatch%ncohorts,1,temppatch%ncohorts)
     
     ! Remove the temporary patch
     
     call deallocate_patchtype(temppatch)
     
     inew = size(split_mask)
     do ico = 1,size(split_mask)

        if (split_mask(ico).eq.1) then

           inew = inew+1
           
           ! Half the densities of the oringal cohort

           cpatch%nplant(ico)           = cpatch%nplant(ico) * 0.5
           cpatch%lai(ico)              = cpatch%lai(ico) * 0.5   
           cpatch%veg_water(ico)        = cpatch%veg_water(ico) * 0.5
           cpatch%mean_gpp(ico)         = cpatch%mean_gpp(ico) * 0.5
           cpatch%mean_leaf_resp(ico)   = cpatch%mean_leaf_resp(ico) * 0.5
           cpatch%mean_root_resp(ico)   = cpatch%mean_root_resp(ico) * 0.5
           cpatch%dmean_leaf_resp(ico)  = cpatch%dmean_leaf_resp(ico) * 0.5
           cpatch%dmean_root_resp(ico)  = cpatch%dmean_root_resp(ico) * 0.5
           cpatch%dmean_gpp(ico)        = cpatch%dmean_gpp(ico) * 0.5
           cpatch%dmean_gpp_pot(ico)    = cpatch%dmean_gpp_pot(ico) * 0.5
           cpatch%dmean_gpp_max(ico)    = cpatch%dmean_gpp_max(ico) * 0.5
           cpatch%growth_respiration(ico)  = cpatch%growth_respiration(ico) * 0.5
           cpatch%storage_respiration(ico) = cpatch%storage_respiration(ico) * 0.5
           cpatch%vleaf_respiration(ico)   = cpatch%vleaf_respiration(ico) * 0.5
           cpatch%monthly_dndt(ico)     = cpatch%monthly_dndt(ico) * 0.5
           cpatch%Psi_open(ico)         = cpatch%Psi_open(ico) * 0.5

           cpatch%gpp(ico)               = cpatch%gpp(ico)              * 0.5
           cpatch%leaf_respiration(ico)  = cpatch%leaf_respiration(ico) * 0.5
           cpatch%root_respiration(ico)  = cpatch%root_respiration(ico) * 0.5

           ! Update the heat capacity and the vegetation energy
           call update_veg_energy_ct(cpatch,ico)


           ! Apply those values to the new cohort


!           cpatch%maker(inew) = 1
           
           call copy_cohort_ar(cpatch,ico,inew)

           ! Tweak the heights and DBHs
           cpatch%bdead(ico) = cpatch%bdead(ico)*(1.0 - epsilon)
           cpatch%dbh(ico) = bd2dbh(cpatch%pft(ico), cpatch%bdead(ico))
           cpatch%hite(ico)  = dbh2h(cpatch%pft(ico), cpatch%dbh(ico))

           cpatch%bdead(inew) = 2.0*cpatch%bdead(inew)-cpatch%bdead(ico)
           cpatch%dbh(inew) = bd2dbh(cpatch%pft(inew), cpatch%bdead(inew))
           cpatch%hite(inew) = dbh2h(cpatch%pft(inew), cpatch%dbh(inew))


           ! Update the vegetation energy again, due to tweaks

           call update_veg_energy_ct(cpatch,inew)


        endif

        !! SANITY CHECK
        if(cpatch%balive(ico) < 0.0) then
           print*," -- ERROR IN SPLIT COHORTS -- "
           stop
        endif

     enddo
     deallocate(split_mask)
     deallocate(temppatch)
     return
   end subroutine split_cohorts_ar

   !===================================================================

   subroutine copy_cohort_ar(cpatch,isc,idt)

     use ed_state_vars,only:patchtype
      
     implicit none

     type(patchtype),target :: cpatch
     integer :: isc,idt
     integer :: imonth

     cpatch%pft(idt)    = cpatch%pft(isc)
     cpatch%nplant(idt) = cpatch%nplant(isc)
     cpatch%hite(idt)   = cpatch%hite(isc)
     cpatch%dbh(idt)    = cpatch%dbh(isc)
     cpatch%bdead(idt)  = cpatch%bdead(isc)
     cpatch%bleaf(idt)  = cpatch%bleaf(isc)
     cpatch%phenology_status(idt) = cpatch%phenology_status(isc)
     cpatch%balive(idt) = cpatch%balive(isc)
     cpatch%lai(idt)    = cpatch%lai(isc)
     cpatch%bstorage(idt) = cpatch%bstorage(isc)
     cpatch%veg_energy(idt) = cpatch%veg_energy(isc)
     cpatch%hcapveg(idt) = cpatch%hcapveg(isc)
     cpatch%veg_temp(idt) = cpatch%veg_temp(isc)
     cpatch%veg_water(idt) = cpatch%veg_water(isc)
     cpatch%mean_gpp(idt) = cpatch%mean_gpp(isc)
     cpatch%mean_leaf_resp(idt) = cpatch%mean_leaf_resp(isc)
     cpatch%mean_root_resp(idt) = cpatch%mean_root_resp(isc)
     cpatch%dmean_leaf_resp(idt) = cpatch%dmean_leaf_resp(isc)
     cpatch%dmean_root_resp(idt) = cpatch%dmean_root_resp(isc)
     cpatch%dmean_gpp(idt) = cpatch%dmean_gpp(isc)
     cpatch%dmean_gpp_pot(idt) = cpatch%dmean_gpp_pot(isc)
     cpatch%dmean_gpp_max(idt) = cpatch%dmean_gpp_max(isc)
     cpatch%growth_respiration(idt) = cpatch%growth_respiration(isc)
     cpatch%storage_respiration(idt) = cpatch%storage_respiration(isc)
     cpatch%vleaf_respiration(idt) = cpatch%vleaf_respiration(isc)
     cpatch%fsn(idt) = cpatch%fsn(isc)
     cpatch%monthly_dndt(idt) = cpatch%monthly_dndt(isc)
     cpatch%Psi_open(idt) = cpatch%Psi_open(isc)
     do imonth = 1,13
        cpatch%cb(imonth,idt) = cpatch%cb(imonth,isc)
        cpatch%cb_max(imonth,idt) = cpatch%cb_max(imonth,isc)
     enddo
     cpatch%cbr_bar(idt) = cpatch%cbr_bar(isc)
     cpatch%krdepth(idt) = cpatch%krdepth(isc)
     cpatch%first_census(idt) = cpatch%first_census(isc)
     cpatch%new_recruit_flag(idt) = cpatch%new_recruit_flag(isc)
!================================================================================

     cpatch%par_v(idt) = cpatch%par_v(isc)
     cpatch%par_v_beam(idt) = cpatch%par_v_beam(isc)
     cpatch%par_v_diffuse(idt) = cpatch%par_v_diffuse(isc)
     cpatch%rshort_v(idt) = cpatch%rshort_v(isc)
     cpatch%rshort_v_beam(idt) = cpatch%rshort_v_beam(isc)
     cpatch%rshort_v_diffuse(idt) = cpatch%rshort_v_diffuse(isc)
     cpatch%rlong_v(idt) = cpatch%rlong_v(isc)
     cpatch%rlong_v_surf(idt) = cpatch%rlong_v_surf(isc)
     cpatch%rlong_v_incid(idt) = cpatch%rlong_v_incid(isc)
     cpatch%rb(idt) = cpatch%rb(isc)
     cpatch%A_open(idt) = cpatch%A_open(isc)
!     cpatch%status(idt) = cpatch%status(isc)

     cpatch%gpp(idt)               = cpatch%gpp(isc)
     cpatch%leaf_respiration(idt)  = cpatch%leaf_respiration(isc)
     cpatch%root_respiration(idt)  = cpatch%root_respiration(isc)
     
     return
   end subroutine copy_cohort_ar



   !===========================================================================

   subroutine  fuse_2_cohorts_ar(cpatch,ico1,ico2, newn,green_leaf_factor, lsl)
 
     use ed_state_vars,only:patchtype
     use pft_coms, only: q, qsw, sla
     use consts_coms,only:t3ple,alli,cliq,cice
     use therm_lib,only:calc_hcapveg,update_veg_energy_ct

     implicit none
     type(patchtype),target :: cpatch
     integer :: ico1,ico2
     real, intent(in) :: newn
     real, intent(in) :: green_leaf_factor
     integer, intent(in) :: lsl
     
     real :: newni
     real :: cb_max
     real :: root_depth
     
     real    , external :: calc_root_depth
     integer , external :: assign_root_depth
     real    , external :: bd2dbh
     real    , external :: dbh2h

     real :: laiold
     
     newni = 1.0 / newn
     

     ! Conserve carbon by calculating bdead first.
     cpatch%bdead(ico2) = (cpatch%nplant(ico2) * cpatch%bdead(ico2) + cpatch%nplant(ico1) * cpatch%bdead(ico1)) * newni

     ! Then get dbh and hite from bdead.
     cpatch%dbh(ico2) = bd2dbh(cpatch%pft(ico2), cpatch%bdead(ico2))
     cpatch%hite(ico2) = dbh2h(cpatch%pft(ico2), cpatch%dbh(ico2))

     ! keep the phenology_status of cc and conserve carbon to get balive.
     cpatch%balive(ico2) = (cpatch%nplant(ico2) * cpatch%balive(ico2) + cpatch%nplant(ico1) * cpatch%balive(ico1)) *newni

     ! Update bleaf and lai.
     if(cpatch%phenology_status(ico2) < 2)then
        cpatch%bleaf(ico2) = green_leaf_factor * cpatch%balive(ico2) /   &
             (1.0 + q(cpatch%pft(ico2)) + cpatch%hite(ico2) * qsw(cpatch%pft(ico2)))
        laiold=cpatch%lai(ico2)
        cpatch%lai(ico2) = cpatch%bleaf(ico2) * sla(cpatch%pft(ico2)) * newn
        ! write (unit=62,fmt='(a,i5,1x,a,1x,3(f13.2,1x))') 'PFT=',cpatch%pft(ico2),'(LAI1,LAI2,LAIM)=',cpatch%lai(ico1),laiold,cpatch%lai(ico2)
     else
        cpatch%lai(ico2) = 0.0
     endif

     cpatch%bstorage(ico2) = (cpatch%nplant(ico2) * cpatch%bstorage(ico2) + cpatch%nplant(ico1) *  &
          cpatch%bstorage(ico1)) * newni

     ! Do our best to conserve energy.  May be problematic
     ! if cohorts have different phenology_status codes.
     if(cpatch%phenology_status(ico2) < 2 .and. cpatch%phenology_status(ico1) < 2)then
        
        cpatch%veg_water(ico2) = (cpatch%veg_water(ico2) * cpatch%nplant(ico2) +   &
             cpatch%veg_water(ico1) * cpatch%nplant(ico1)) * newni

     elseif(cpatch%phenology_status(ico2) < 2)then

        cpatch%veg_water(ico2) = cpatch%veg_water(ico2) * cpatch%nplant(ico2) * newni

     endif

     ! No need to modify hcapveg until the new implementation.
     
     cpatch%mean_gpp(ico2) = (cpatch%mean_gpp(ico2) * cpatch%nplant(ico2) + cpatch%mean_gpp(ico1) *   &
          cpatch%nplant(ico1)) * newni

     cpatch%mean_leaf_resp(ico2) = (cpatch%mean_leaf_resp(ico2) * cpatch%nplant(ico2) +   &
          cpatch%mean_leaf_resp(ico1) * cpatch%nplant(ico1)) * newni

     cpatch%mean_root_resp(ico2) = (cpatch%mean_root_resp(ico2) * cpatch%nplant(ico2) +   &
          cpatch%mean_root_resp(ico1) * cpatch%nplant(ico1)) * newni
     
     cpatch%growth_respiration(ico2) = (cpatch%growth_respiration(ico2) * cpatch%nplant(ico2) +  &
          cpatch%growth_respiration(ico1) * cpatch%nplant(ico1)) * newni
     
     cpatch%storage_respiration(ico2) = (cpatch%storage_respiration(ico2) * cpatch%nplant(ico2) +  &
          cpatch%storage_respiration(ico1) * cpatch%nplant(ico1)) * newni
     
     cpatch%vleaf_respiration(ico2) = (cpatch%vleaf_respiration(ico2) * cpatch%nplant(ico2) +  &
          cpatch%vleaf_respiration(ico1) * cpatch%nplant(ico1)) * newni
     
     cpatch%fsn(ico2) = (cpatch%fsn(ico2) * cpatch%nplant(ico2) + cpatch%fsn(ico1) * cpatch%nplant(ico1)) * newni

     cpatch%Psi_open(ico2) = (cpatch%Psi_open(ico2) * cpatch%nplant(ico2) + cpatch%Psi_open(ico1) *   &
          cpatch%nplant(ico1)) * newni

     cpatch%cb(1:13,ico2) = (cpatch%cb(1:13,ico2) * cpatch%nplant(ico2) + cpatch%nplant(ico1) *   &
          cpatch%cb(1:13,ico1)) * newni

     cpatch%cb_max(1:13,ico2) = (cpatch%cb_max(1:13,ico2) * cpatch%nplant(ico2) +   &
          cpatch%nplant(ico1) * cpatch%cb_max(1:13,ico1)) * newni

     cpatch%gpp(ico2)              = ( cpatch%gpp(ico2)              * cpatch%nplant(ico2) &
                                   +   cpatch%gpp(ico1)              * cpatch%nplant(ico1) )* newni
     cpatch%leaf_respiration(ico2) = ( cpatch%leaf_respiration(ico2) * cpatch%nplant(ico2) &
                                   +   cpatch%leaf_respiration(ico1) * cpatch%nplant(ico1) )* newni
     cpatch%root_respiration(ico2) = ( cpatch%root_respiration(ico2) * cpatch%nplant(ico2) &
                                   +   cpatch%root_respiration(ico1) * cpatch%nplant(ico1) )* newni


     cb_max = sum(cpatch%cb_max(1:12,ico2))
     if(cb_max > 0.0)then
        cpatch%cbr_bar(ico2) = sum(cpatch%cb(1:12,ico2)) / cb_max
     else
        cpatch%cbr_bar(ico2) = 0.0
     endif
     
     root_depth = calc_root_depth(cpatch%hite(ico2), cpatch%dbh(ico2), cpatch%pft(ico2))
     cpatch%krdepth(ico2) = assign_root_depth(root_depth, lsl)

     cpatch%nplant(ico2) = newn

     !---- Will this conserve energy? Is the total LAI conserved? I hope so...
     ! I think this is fine because temperature is no longer prognostic, so temperature will be 
     ! adjusted at the next fast step...
     !cpatch%veg_energy(ico2) = (cpatch%veg_energy(ico2) * cpatch%nplant(ico2) +   &
     !     cpatch%veg_energy(ico1) * cpatch%nplant(ico1)) * newni
     
     ! Think it is safer to give the fused cohort an energy that is based off of it's
     ! new temperature
     !=======================================================================================

     call update_veg_energy_ct(cpatch,ico2)
     


     return
   end subroutine fuse_2_cohorts_ar

!-----------------------------------------------------------------------------

   subroutine fuse_patches_ar(cgrid)

     use ed_state_vars,only : edtype,polygontype,sitetype,patchtype
     use fusion_fission_coms, only : ff_ndbh, ntol, profile_tol
     use max_dims, only: n_pft
     use mem_sites, only: maxpatch,maxcohort

     implicit none

     type(edtype),target :: cgrid
     type(polygontype),pointer :: cpoly
     type(sitetype),pointer :: csite
     type(patchtype),pointer :: cpatch
     type(sitetype),pointer  :: tempsite
     integer :: ipy,isi,ipa,ipa_next,ipa_tp
     integer :: i,j,istop,npatches
     real :: norm,tolerance_mult
     integer,allocatable :: fuse_table(:)
     real, parameter :: tolerance_max=2.0
     integer :: npatches_old

     ! Allocate the swapper patches in the site type

     allocate(tempsite)

     do ipy = 1,cgrid%npolygons

     cpoly => cgrid%polygon(ipy)
        
     do isi = 1,cpoly%nsites

     csite => cpoly%site(isi)

     call allocate_sitetype(tempsite, csite%npatches )
     allocate(fuse_table(csite%npatches))
     fuse_table(:) = 1

     ! ALGORITHM
     ! set all fusion flags to true
     ! create size profiles
     ! goto every patch
     ! find next older patch with same dist_type
     ! check fusion criterion
     ! if within criterion, fuse, otherwise, skip


     ! Loop from the youngest to oldest patch

     do ipa = csite%npatches,1,-1
        call patch_pft_size_profile_ar(csite,ipa,ff_ndbh,cpoly%green_leaf_factor(:,isi))
     end do
     
 
     ! loop over sites
     tolerance_mult = 1.0
     max_patch: do
        npatches_old = sum(fuse_table)
        csite%fuse_flag(:) = 1
        
        ! Loop from youngest to the second oldest patch
        do ipa = csite%npatches,2,-1
           
           cpatch => csite%patch(ipa)

           if (fuse_table(ipa).ne.0) then
              
              ! Look at the next patch which is immediately older and 
              ! of the same disturbance type
              ipa_next = ipa - 1
              istop = 0
              ipa_tp = -1
              do while ((ipa_next > 0) .and. (istop == 0))

                 if( csite%dist_type(ipa) == csite%dist_type(ipa_next) .and. &
                    fuse_table(ipa_next).ne.0 ) then
                    istop  = 1
                    ipa_tp = ipa_next
                 endif
                 ipa_next = ipa_next - 1
              enddo


              ! Once we have identified the patch with the same disturbance type
              ! and closest age (ipa_tp), determine if it is similar enough to average (fuse)
              ! the two together.
              
              if (ipa_tp > 0 ) then ! found a fusion candidate?
                 
                 !  fusion criterion
                 do i=1,n_pft          ! loop over pft 
                    do j=1,ff_ndbh !      loop over hgt bins
                       
                       if(csite%pft_density_profile(i,j,ipa) > tolerance_mult * ntol  &
                            .or. csite%pft_density_profile(i,j,ipa_tp) >   &
                            tolerance_mult * ntol)then
                          
                          ! This is the normalized difference in their biodensity profiles
                          ! If the normalized difference is greater than the tolerance
                          ! for any of the pfts and dbh classes, then reject them as similar
                          norm = abs(csite%pft_density_profile(i,j,ipa)  &
                                  -csite%pft_density_profile(i,j,ipa_tp))  &
                                  /(0.5*(csite%pft_density_profile(i,j,ipa)  &
                                  +csite%pft_density_profile(i,j,ipa_tp)))

                          if(norm > profile_tol) csite%fuse_flag(ipa)=0   ! reject
                          
                       endif
                    enddo
                 enddo

          
                               
                 ! Create a mapping of the patches that fuse together
                 if(csite%fuse_flag(ipa) == 1)then
                                    

                    ! Take an average of the patch properties at index ipa and ipa_tp
                    ! assign the average to index ipa_tp

                    call fuse_2_patches_ar(csite,ipa,ipa_tp,cpoly%met(isi)%rhos, &
                         cpoly%lsl(isi),cpoly%green_leaf_factor(:,isi))

                    ! Recalculate the pft size profile for the averaged patch at ipa_tp
                    call patch_pft_size_profile_ar(csite,ipa_tp,ff_ndbh,cpoly%green_leaf_factor(:,isi))

                    ! The patch at index ipa is no longer valid, it should be flagged as such
                    fuse_table(ipa) = 0
                    
                 endif

              endif
                    
           endif
           
           
           
        enddo
        
        
        npatches = sum(fuse_table)

        if(npatches <= maxpatch)exit max_patch
        if((npatches == npatches_old) .and. (tolerance_mult > tolerance_max)) exit max_patch
        tolerance_mult = tolerance_mult * 1.01
     enddo max_patch
     
     
     ! Set the number of patches in the site to "npatches"
     tempsite%npatches = npatches
     
     ! Copy the selected data into the temporary space, args 1 and 3 must be dimension of arg 4
     ! argument 2 must be the dimension of the sum of the 3rd argument.
     
     call copy_sitetype_mask(csite,tempsite,fuse_table,size(fuse_table),npatches)
     
     call deallocate_sitetype(csite)
     
     ! Reallocate the current site
     call allocate_sitetype(csite,npatches)

     ! Copy the selected temporary data into the orignal site vectors
     fuse_table = 0
     fuse_table(1:npatches) = 1
     
     call copy_sitetype_mask(tempsite,csite,fuse_table,size(fuse_table),npatches)

     ! The new and fused csite is now complete, clean up the temporary data
     
     deallocate(fuse_table)

     call deallocate_sitetype(tempsite)
     

  enddo
  
enddo
deallocate(tempsite)


  return
end subroutine fuse_patches_ar
!-----------------------------------------------------------------------------

   !===============================================================

   subroutine fuse_2_patches_ar(csite,dp,rp,rhos,lsl,green_leaf_factor)

     use ed_state_vars,only:sitetype,patchtype
     use therm_lib, only: update_veg_energy_ct
     use soil_coms, only: soil
     use grid_coms, only: nzg, nzs
     use fusion_fission_coms, only: ff_ndbh
     use max_dims, only: n_pft,n_dbh
     use consts_coms, only: cpi, cpor, p00
     
     implicit none
     
     integer :: dp,rp
     real, intent(in) :: rhos
     integer :: ico,ndc,nrc
     real, dimension(n_pft), intent(in) :: green_leaf_factor
     integer :: lsl
     type(sitetype),target :: csite
     type(patchtype),pointer :: cpatch
     type(patchtype), pointer :: newpatch
     
     !  This function fuses the two patches specified in the argument.
     !  It fuses the first patch in the argument (the "donor" = dp ) into the second
     !  patch in the argument (the "recipient" = rp ), and frees the memory 
     !  associated with the second patch

     real :: newarea,newareai
     
     ! new area
     newarea = csite%area(dp) + csite%area(rp)
     newareai = 1.0/newarea

     csite%age(rp)                   = (csite%age(dp) * csite%area(dp)                          &
                                     +  csite%age(rp) * csite%area(rp)) * newareai
     csite%fast_soil_C(rp)           = (csite%fast_soil_C(dp)*csite%area(dp)                    &
                                     +  csite%fast_soil_C(rp)*csite%area(rp))  * newareai  
     csite%slow_soil_C(rp)           = (csite%slow_soil_C(dp)*csite%area(dp)                    &
                                     +  csite%slow_soil_C(rp)*csite%area(rp))  * newareai  
     csite%structural_soil_C(rp)     = (csite%structural_soil_C(dp)*csite%area(dp)              &
                                     +  csite%structural_soil_C(rp)*csite%area(rp)) * newareai  
     csite%structural_soil_L(rp)     = (csite%structural_soil_L(dp)*csite%area(dp)              &
                                     +  csite%structural_soil_L(rp)*csite%area(rp)) * newareai  
     csite%mineralized_soil_N(rp)    = (csite%mineralized_soil_N(dp)*csite%area(dp)             &
                                     +  csite%mineralized_soil_N(rp)*csite%area(rp)) * newareai 
     csite%fast_soil_N(rp)           = (csite%fast_soil_N(dp)*csite%area(dp)                    &
                                     +  csite%fast_soil_N(rp)*csite%area(rp)) * newareai  
     csite%sum_dgd(rp)               = (csite%sum_dgd(dp) * csite%area(dp)                      &
                                     +  csite%sum_dgd(rp) * csite%area(rp)) * newareai

     csite%sum_chd(rp)               = (csite%sum_chd(dp) * csite%area(dp)                      &
                                     +  csite%sum_chd(rp) * csite%area(rp)) * newareai
     csite%can_co2(rp)               = (csite%can_co2(dp) * csite%area(dp)                      &
                                     +  csite%can_co2(rp) * csite%area(rp)) * newareai

     csite%can_temp(rp)              = (csite%can_temp(dp) * csite%area(dp)                     &
                                     +  csite%can_temp(rp) * csite%area(rp)) * newareai
     csite%can_shv(rp)               = (csite%can_shv(dp) * csite%area(dp)                      &
                                     +  csite%can_shv(rp) * csite%area(rp))                     &
                                     *  newareai

     csite%sfcwater_energy(1:nzs,rp) &
            = (csite%sfcwater_energy(1:nzs,rp) * csite%sfcwater_mass(1:nzs,rp) * csite%area(rp) &
            +  csite%sfcwater_energy(1:nzs,dp) * csite%sfcwater_mass(1:nzs,dp) * csite%area(dp))&
            * newareai
     csite%sfcwater_mass(1:nzs,rp)   = (csite%sfcwater_mass(1:nzs,rp) * csite%area(rp)          &
                                     +  csite%sfcwater_mass(1:nzs,dp) * csite%area(dp))         &
                                     * newareai
     csite%sfcwater_depth(1:nzs,rp)  = (csite%sfcwater_depth(1:nzs,rp) * csite%area(rp)         &
                                     +  csite%sfcwater_depth(1:nzs,dp) * csite%area(dp))        &
                                     * newareai

     csite%soil_energy(1:nzg,rp)     = (csite%soil_energy(1:nzg,dp) * csite%area(dp)            &
                                     +  csite%soil_energy(1:nzg,rp) * csite%area(rp)) * newareai

     csite%soil_water(1:nzg,rp)      = (csite%soil_water(1:nzg,rp) * dble(csite%area(rp))             &
                                     +  csite%soil_water(1:nzg,dp) * dble(csite%area(dp))) * dble(newareai)

     !-----------------------------------------------------!
     ! This subroutine takes care of filling:              !
     !                                                     !
     ! + csite%ground_shv(rp)                              !
     ! + csite%surface_ssh(rp)                             !
     ! + csite%soil_tempk(k,rp)                            !
     ! + csite%soil_fracliq(k,rp)                          !
     ! + csite%nlev_sfcwater(rp)                           !
     ! + csite%sfcwater_energy(k,rp)                       !
     ! + csite%csite%sfcwater_tempk(k,rp)                  !
     ! + csite%sfcwater_fracliq(k,rp)                      !
     call new_patch_sfc_props_ar(csite,rp,rhos)            !
     !-----------------------------------------------------!

     csite%mean_rh(rp) = (csite%mean_rh(rp) * csite%area(rp)                                   &
                       +  csite%mean_rh(dp) * csite%area(dp)) * newareai
     
     csite%dmean_A_decomp(rp) = (csite%dmean_A_decomp(rp) * csite%area(rp)                     &
                              +  csite%dmean_A_decomp(dp) * csite%area(dp)) * newareai
     
     csite%dmean_Af_decomp(rp) = (csite%dmean_Af_decomp(rp) * csite%area(rp)                   &
                               + csite%dmean_Af_decomp(dp) * csite%area(dp)) * newareai
     
     csite%repro(1:n_pft,rp) = (csite%repro(1:n_pft,rp) * csite%area(rp)                       &
                             +  csite%repro(1:n_pft,dp) * csite%area(dp)) * newareai
     
     csite%watertable(rp) = (csite%watertable(rp) * csite%area(rp)                             &
                          +  csite%watertable(dp)*csite%area(dp)) *newareai
     
     ! Even though these variables are not prognostic, they need to be copied so the output will have the values.
     ! Other variables will probably be scaled here as well
     csite%avg_carbon_ac(rp)         = (csite%avg_carbon_ac(rp)         * csite%area(rp)        &
                                     +  csite%avg_carbon_ac(dp)         * csite%area(dp))       &
                                     *  newareai
     csite%avg_vapor_vc(rp)          = (csite%avg_vapor_vc(rp)          * csite%area(rp)        &
                                     +  csite%avg_vapor_vc(dp)          * csite%area(dp))       &
                                     *  newareai
     csite%avg_dew_cg(rp)            = (csite%avg_dew_cg(rp)            * csite%area(rp)        &
                                     +  csite%avg_dew_cg(dp)            * csite%area(dp))       &
                                     *  newareai
     csite%avg_vapor_gc(rp)          = (csite%avg_vapor_gc(rp)          * csite%area(rp)        &
                                     +  csite%avg_vapor_gc(dp)          * csite%area(dp))       &
                                     *  newareai
     csite%avg_wshed_vg(rp)          = (csite%avg_wshed_vg(rp)          * csite%area(rp)        &
                                     +  csite%avg_wshed_vg(dp)          * csite%area(dp))       &
                                     *  newareai
     csite%avg_vapor_ac(rp)          = (csite%avg_vapor_ac(rp)          * csite%area(rp)        &
                                     +  csite%avg_vapor_ac(dp)          * csite%area(dp))       &
                                     *  newareai
     csite%avg_transp(rp)            = (csite%avg_transp(rp)            * csite%area(rp)        &
                                     +  csite%avg_transp(dp)            * csite%area(dp))       &
                                     *  newareai
     csite%avg_evap(rp)              = (csite%avg_evap(rp)              * csite%area(rp)        &
                                     +  csite%avg_evap(dp)              * csite%area(dp))       &
                                     *  newareai
     csite%avg_smoist_gg(1:nzg,rp)   = (csite%avg_smoist_gg(1:nzg,rp)   * csite%area(rp)        &
                                     +  csite%avg_smoist_gg(1:nzg,dp)   * csite%area(dp))       &
                                     *  newareai
     csite%avg_smoist_gc(1:nzg,rp)   = (csite%avg_smoist_gc(1:nzg,rp)   * csite%area(rp)        &
                                     +  csite%avg_smoist_gc(1:nzg,dp)   * csite%area(dp))       &
                                     *  newareai
     csite%avg_runoff(rp)            = (csite%avg_runoff(rp)            * csite%area(rp)        &
                                     +  csite%avg_runoff(dp)            * csite%area(dp))       &
                                     *  newareai
     csite%aux(rp)                   = (csite%aux(rp)                   * csite%area(rp)        &
                                     +  csite%aux(dp)                   * csite%area(dp))       &
                                     *  newareai
     csite%aux_s(1:nzg,rp)           = (csite%aux_s(1:nzg,rp)           * csite%area(rp)        &
                                     +  csite%aux_s(1:nzg,dp)           * csite%area(dp))       &
                                     *  newareai
     csite%avg_sensible_vc(rp)       = (csite%avg_sensible_vc(rp)       * csite%area(rp)        &
                                     +  csite%avg_sensible_vc(dp)       * csite%area(dp))       &
                                     *  newareai
     csite%avg_sensible_2cas(rp)     = (csite%avg_sensible_2cas(rp)     * csite%area(rp)        &
                                     +  csite%avg_sensible_2cas(dp)     * csite%area(dp))       &
                                     *  newareai
     csite%avg_qwshed_vg(rp)         = (csite%avg_qwshed_vg(rp)         * csite%area(rp)        &
                                     +  csite%avg_qwshed_vg(dp)         * csite%area(dp))       &
                                     *  newareai
     csite%avg_sensible_gc(rp)       = (csite%avg_sensible_gc(rp)       * csite%area(rp)        &
                                     +  csite%avg_sensible_gc(dp)       * csite%area(dp))       &
                                     *  newareai
     csite%avg_sensible_ac(rp)       = (csite%avg_sensible_ac(rp)       * csite%area(rp)        &
                                     +  csite%avg_sensible_ac(dp)       * csite%area(dp))       &
                                     *  newareai
     csite%avg_sensible_tot(rp)      = (csite%avg_sensible_tot(rp)      * csite%area(rp)        &
                                     +  csite%avg_sensible_tot(dp)      * csite%area(dp))       &
                                     *  newareai
     csite%avg_sensible_gg(1:nzg,rp) = (csite%avg_sensible_gg(1:nzg,rp) * csite%area(rp)        &
                                     +  csite%avg_sensible_gg(1:nzg,dp) * csite%area(dp))       &
                                     *  newareai
     csite%avg_runoff_heat(rp)       = (csite%avg_runoff_heat(rp)       * csite%area(rp)        &
                                     +  csite%avg_runoff_heat(dp)       * csite%area(dp))       &
                                     *  newareai
     csite%avg_heatstor_veg(rp)      = (csite%avg_heatstor_veg(rp)      * csite%area(rp)        &
                                     +  csite%avg_heatstor_veg(dp)      * csite%area(dp))       &
                                     *  newareai
     csite%avg_veg_energy(rp)        = (csite%avg_veg_energy(rp)        * csite%area(rp)        &
                                     +  csite%avg_veg_energy(dp)        * csite%area(dp))       &
                                     *  newareai
     csite%avg_veg_temp(rp)          = (csite%avg_veg_temp(rp)          * csite%area(rp)        &
                                     +  csite%avg_veg_temp(dp)          * csite%area(dp))       &
                                     *  newareai
     csite%avg_veg_water(rp)         = (csite%avg_veg_water(rp)         * csite%area(rp)        &
                                     +  csite%avg_veg_water(dp)         * csite%area(dp))       &
                                     *  newareai

     csite%co2budget_gpp(rp)             = (csite%co2budget_gpp(rp)              * csite%area(rp)  &
                                         +  csite%co2budget_gpp(dp)              * csite%area(dp)) &
                                         *  newareai
     csite%co2budget_gpp_dbh(1:n_dbh,rp) = (csite%co2budget_gpp_dbh(1:n_dbh,rp)  * csite%area(rp)  &
                                         +  csite%co2budget_gpp_dbh(1:n_dbh,dp)  * csite%area(dp)) &
                                         *  newareai
     csite%co2budget_plresp(rp)          = (csite%co2budget_plresp(rp)           * csite%area(rp)  &
                                         +  csite%co2budget_plresp(dp)           * csite%area(dp)) &
                                         *  newareai
     csite%co2budget_rh(rp)              = (csite%co2budget_rh(rp)               * csite%area(rp)  &
                                         +  csite%co2budget_rh(dp)               * csite%area(dp)) &
                                         *  newareai

     ! adjust densities of cohorts in recipient patch

     cpatch => csite%patch(rp)
     nrc = cpatch%ncohorts
     do ico = 1,nrc
        cpatch%lai(ico) = cpatch%lai(ico) * csite%area(rp) * newareai
        cpatch%nplant(ico) = cpatch%nplant(ico) * csite%area(rp) * newareai
        cpatch%veg_water(ico) = cpatch%veg_water(ico) * csite%area(rp) * newareai
        cpatch%mean_gpp(ico) = cpatch%mean_gpp(ico) * csite%area(rp) * newareai
        cpatch%mean_leaf_resp(ico) = cpatch%mean_leaf_resp(ico) * csite%area(rp) * newareai
        cpatch%mean_root_resp(ico) = cpatch%mean_root_resp(ico) * csite%area(rp) * newareai
        cpatch%growth_respiration(ico) = cpatch%growth_respiration(ico) * csite%area(rp) * newareai
        cpatch%storage_respiration(ico) = cpatch%storage_respiration(ico) * csite%area(rp) * newareai
        cpatch%vleaf_respiration(ico) = cpatch%vleaf_respiration(ico) * csite%area(rp) * newareai
        cpatch%Psi_open(ico) = cpatch%Psi_open(ico) * csite%area(rp) * newareai
        
        ! These were giving problem at new year
        cpatch%gpp(ico)              = cpatch%gpp(ico)              * csite%area(rp) * newareai
        cpatch%leaf_respiration(ico) = cpatch%leaf_respiration(ico) * csite%area(rp) * newareai
        cpatch%root_respiration(ico) = cpatch%root_respiration(ico) * csite%area(rp) * newareai

        ! Now that the plant density and the amount of water has changed, there will
        ! be an inconsistency in the energy. If the energy stays the same, but temperature
        ! is diagnosed with a sudden drop in biomass or water, we will have sky-rocketing values
        ! so we have to adjust the energy accordingly also.

        call update_veg_energy_ct(cpatch,ico)


     enddo

     ! adjust densities of cohorts in donor patch 
     cpatch => csite%patch(dp)
     ndc = cpatch%ncohorts
     do ico = 1,ndc

        cpatch%nplant(ico) = cpatch%nplant(ico) * csite%area(dp) * newareai
        cpatch%lai(ico) = cpatch%lai(ico) * csite%area(dp) * newareai
        cpatch%veg_water(ico) = cpatch%veg_water(ico) * csite%area(dp) * newareai
        cpatch%mean_gpp(ico) = cpatch%mean_gpp(ico) * csite%area(dp) * newareai
        cpatch%mean_leaf_resp(ico) = cpatch%mean_leaf_resp(ico) * csite%area(dp) * newareai
        cpatch%mean_root_resp(ico) = cpatch%mean_root_resp(ico) * csite%area(dp) * newareai
        cpatch%growth_respiration(ico) = cpatch%growth_respiration(ico) * csite%area(dp) * newareai
        cpatch%storage_respiration(ico) = cpatch%storage_respiration(ico) * csite%area(dp) * newareai
        cpatch%vleaf_respiration(ico) = cpatch%vleaf_respiration(ico) * csite%area(dp) * newareai
        cpatch%Psi_open(ico) = cpatch%Psi_open(ico) * csite%area(dp) * newareai

        cpatch%gpp(ico)      = cpatch%gpp(ico)                      * csite%area(dp) * newareai
        cpatch%leaf_respiration(ico) = cpatch%leaf_respiration(ico) * csite%area(dp) * newareai
        cpatch%root_respiration(ico) = cpatch%root_respiration(ico) * csite%area(dp) * newareai

        
        ! Now that the plant density and the amount of water has changed, there will
        ! be an inconsistency in the energy. If the energy stays the same, but temperature
        ! is diagnosed with a sudden drop in biomass or water, we will have sky-rocketing values
        ! so we have to adjust the energy accordingly also.  This linear scaling should
        ! be suitable, because vegetation energy is a linear combination of biomass and water
        ! mutliplied by temperature.

        call update_veg_energy_ct(cpatch,ico)

     enddo
     
     ! Fill a new patch with the donor and recipient cohort vectors


     ! Allocate the new patch and its cohorts
     nullify(newpatch)
     allocate(newpatch)
     call allocate_patchtype(newpatch,ndc + nrc )


     call copy_patchtype(csite%patch(rp),newpatch,1,nrc,1,nrc)
 
     call copy_patchtype(csite%patch(dp),newpatch,1,ndc,nrc+1,nrc+ndc)

!     call deallocate_patchtype(csite%patch(dp))

     call deallocate_patchtype(csite%patch(rp))

     call allocate_patchtype(csite%patch(rp),ndc+nrc)

     call copy_patchtype(newpatch,csite%patch(rp),1,nrc+ndc,1,nrc+ndc)


     call deallocate_patchtype(newpatch)
     deallocate(newpatch)


     !-----------------------------------------------------!
     ! This subroutine takes care of filling:              !
     !                                                     !
     ! + csite%veg_height(rp)                              !
     ! + csite%lai(rp)                                     !
     ! + csite%veg_rough(rp)                               !
     ! + csite%wbudget_initialstorage(rp)                  !
     ! + csite%ebudget_initialstorage(rp)                  !
     ! + csite%co2budget_initialstorage(rp)                !
     call update_patch_derived_props_ar(csite,lsl, rhos,rp)
     !-----------------------------------------------------!

     !-----------------------------------------------------!
     !    This subroutine takes care of updating size      !
     ! profile within patch.                               !
     !                                                     !
     ! + csite%pft_density_profile(:,:,rp)                 !
     call patch_pft_size_profile_ar(csite,rp,ff_ndbh,green_leaf_factor)
     !-----------------------------------------------------!

     ! update patch area
     csite%area(rp) = newarea


     return
   end subroutine fuse_2_patches_ar
   !=====================================================================
   subroutine patch_pft_size_profile_ar(csite,ipa,nbins,green_leaf_factor)

     use ed_state_vars,only:sitetype,patchtype

     use fusion_fission_coms, only : maxdbh
     use max_dims, only : n_pft
      
     implicit none

     type(sitetype),target :: csite
     type(patchtype),pointer :: cpatch

     integer :: ipa
     real, dimension(n_pft), intent(in) :: green_leaf_factor
     integer :: nbins
     
     real    :: rmin,rmax,dh
     integer :: i,j,ico
     real    :: elong

     dh = maxdbh/real(nbins)

     ! initialize bins
     do i=1,n_pft
        do j=1,nbins
           csite%pft_density_profile(i,j,ipa)=0.0
        enddo
     enddo

     ! update bins
     cpatch => csite%patch(ipa)
     do ico = 1,cpatch%ncohorts

        elong = green_leaf_factor(cpatch%pft(ico))
        do j=1,nbins
           if (j == 1)then
              rmin = 0.0
           else
              rmin = real((j-1))*dh
           endif
           rmax = real(j)*dh

           if(cpatch%dbh(ico) > rmin .and. cpatch%dbh(ico) <= rmax)then
              csite%pft_density_profile(cpatch%pft(ico),j,ipa) =   &
                   csite%pft_density_profile(cpatch%pft(ico),j,ipa) + cpatch%nplant(ico)
           endif
        enddo
        
        ! deal with largest dbh bin
        j = nbins
        rmin = real(j)*dh
        if(cpatch%dbh(ico) > rmin)then
           csite%pft_density_profile(cpatch%pft(ico),j,ipa) =   &
                csite%pft_density_profile(cpatch%pft(ico),j,ipa) + cpatch%nplant(ico)
        endif
        
     enddo
     
     return
     
   end subroutine patch_pft_size_profile_ar
   

 end module fuse_fiss_utils_ar
