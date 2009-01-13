subroutine bare_ground_init(cgrid)
  ! This subroutine initializes a near-bare ground polygon.
  use ed_state_vars, only: edtype,polygontype,sitetype &
                          ,allocate_sitetype

  implicit none

  type(edtype)      , target  :: cgrid
  type(polygontype) , pointer :: cpoly
  type(sitetype)    , pointer :: csite
  integer :: ipy,isi

  ! Loop over all polygons
  do ipy=1,cgrid%npolygons
     cpoly => cgrid%polygon(ipy)
    
     do isi=1,cpoly%nsites
        csite => cpoly%site(isi)

        csite%npatches = 1
        call allocate_sitetype(csite,1)
        csite%dist_type          (1) = 3
        csite%age                (1) = 0.0
        csite%area               (1) = 1.0
        csite%fast_soil_C        (1) = 0.2
        csite%slow_soil_C        (1) = 0.01
        csite%structural_soil_C  (1) = 10.0
        csite%structural_soil_L  (1) = csite%structural_soil_C (1)
        csite%mineralized_soil_N (1) = 1.0
        csite%fast_soil_N        (1) = 1.0
        csite%sum_dgd            (1) = 0.0
        csite%sum_chd            (1) = 0.0
        csite%plantation         (1) = 0
        csite%plant_ag_biomass   (1) = 0.

        call init_bare_ground_patchtype(.true.,csite,cpoly%lsl(isi),cpoly%met(isi)%atm_tmp,1,csite%npatches)

        call init_ed_patch_vars_array(csite,1,csite%npatches)
     enddo
     call init_ed_site_vars_array(cpoly,cgrid%lat(ipy))
  enddo
  call init_ed_poly_vars_array(cgrid)

  return
end subroutine bare_ground_init
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine init_bare_ground_patchtype(zero_time,csite,lsl,atm_tmp,ipa_a,ipa_z)
   !
   ! This subroutine assigns a near-bare ground state for the given patch
   !
   use ed_state_vars, only: edtype,polygontype,sitetype,patchtype &
                           ,allocate_sitetype,allocate_patchtype
   use max_dims, only: n_pft
   use pft_coms, only: SLA, q, qsw, hgt_min, include_pft, include_these_pft,include_pft_ag
   use ed_therm_lib, only: calc_hcapveg

   ! This subroutine assigns a near-bare ground state for the given patch
   use consts_coms, only : t3ple
   implicit none
   logical, intent(in)      :: zero_time   ! This is a call done at initialisation.
   type(sitetype), target   :: csite
   integer, intent(in)      :: lsl         ! Lowest soil level
   real   , intent(in)      :: atm_tmp     ! Atmospheric temperature
   integer, intent(in)      :: ipa_a,ipa_z ! Patch indexes. 1st and last patches to
                                           ! initialize with near bare-ground state


   type(patchtype), pointer :: cpatch
   real, external :: h2dbh,dbh2bd,dbh2bl,ed_biomass
   integer :: ipa,ico,mypfts,ipft
   real :: laisum

   do ipa=ipa_a,ipa_z
      cpatch => csite%patch(ipa)

      ! Decide how many cohorts to allocate
      select case (csite%dist_type(ipa))
      case (1)   ! Agriculture
         mypfts = sum(include_pft_ag)
      case (2,3) ! Secondary or primary forest
         mypfts= sum(include_pft)
      end select
      laisum = 0.
      call allocate_patchtype(cpatch,mypfts)
      do ico = 1,mypfts
         ! Include_these_pft is sorted, so just the first elements will be used 
         select case (csite%dist_type(ipa))
         case (1)
            ipft=include_these_pft(ico)
         case (2,3)
            ipft=include_these_pft(ico)
         end select

         cpatch%pft(ico)              = ipft
   
         ! Define the near-bare ground
         cpatch%hite(ico)             = hgt_min(ipft)
         cpatch%dbh(ico)              = h2dbh(cpatch%hite(ico),ipft)
         cpatch%bdead(ico)            = dbh2bd(cpatch%dbh(ico),cpatch%hite(ico),ipft)
         cpatch%bleaf(ico)            = dbh2bl(cpatch%dbh(ico),ipft)
         if (ipft==1 .or. ipft==5) then   ! Grasses need more plants to start or their
                                          ! seed mass will never be large enough to reproduce
            cpatch%nplant(ico)           = 0.6
         else
            cpatch%nplant(ico)           = 0.1
         endif
         cpatch%phenology_status(ico) = 0
         cpatch%balive(ico)           = cpatch%bleaf(ico) * ( 1.0 + q(ipft) +  &
                                        qsw(ipft) * cpatch%hite(ico) )
         cpatch%lai(ico)              = cpatch%bleaf(ico) * cpatch%nplant(ico) * SLA(ipft)
         cpatch%bstorage(ico)         = 0.0   
         csite%plant_ag_biomass(ipa)  = csite%plant_ag_biomass(ipa) +                    &
               ed_biomass(cpatch%bdead(ico),cpatch%balive(ico), cpatch%bleaf(ico)        &
                         ,cpatch%pft(ico), cpatch%hite(ico),cpatch%bstorage(ico))        &
              * cpatch%nplant(ico)           
         
         ! Initialize cohort-level variables
         call init_ed_cohort_vars_array(cpatch,ico,lsl)
         
         laisum = laisum + cpatch%lai(ico)
      end do


      ! If this is not the initial time, set heat capacity for stability.
      do ico = 1,mypfts
         if (.not. zero_time) then
            
            cpatch%veg_temp(ico)  = atm_tmp
            cpatch%veg_water(ico) = 0.0
            
            cpatch%hcapveg(ico) = calc_hcapveg(cpatch%bleaf(ico),cpatch%bdead(ico), &
                 cpatch%nplant(ico),cpatch%pft(ico))
            
            cpatch%veg_energy(ico) = cpatch%hcapveg(ico) * (cpatch%veg_temp(ico)-t3ple)

         end if
      end do
   end do
   return
end subroutine init_bare_ground_patchtype
!==========================================================================================!
!==========================================================================================!
