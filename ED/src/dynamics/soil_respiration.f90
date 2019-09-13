module soil_respiration_module
  contains

!==========================================================================================!
!==========================================================================================!
!     This subroutine computes the soil respiration terms (root and heterotrophic).        !
!------------------------------------------------------------------------------------------!
subroutine soil_respiration(csite,ipa,mzg,ntext_soil)

   use ed_state_vars, only : sitetype                 & ! structure
                           , patchtype                ! ! structure
   use soil_coms    , only : soil                     & ! intent(in)
                           , dslz                     & ! intent(in)
                           , slz                      & ! intent(in)
                           , dolz                     & ! intent(in)
                           , olz                      ! ! intent(in)
   use decomp_coms  , only : k_rh_active              ! ! intent(in)
   use consts_coms  , only : wdns                     & ! intent(in)
                           , umols_2_kgCyr            ! ! intent(in)
   use therm_lib    , only : uextcm2tl                ! ! function
   use ed_misc_coms , only : dtlsm                    & ! intent(in)
                           , frqsum                   & ! intent(in)
                           , ivertresp                ! ! intent(in)
   use grid_coms    , only : nzl

   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   type(sitetype)                , target     :: csite
   integer                       , intent(in) :: ipa
   integer                       , intent(in) :: mzg
   integer       , dimension(mzg), intent(in) :: ntext_soil
   !----- Local variables. ----------------------------------------------------------------!
   type(patchtype)               , pointer    :: cpatch
   integer                                    :: ico
   integer                                    :: ipft
   integer                                    :: k
   integer                                    :: kroot
   integer                                    :: nsoil
   real                                       :: Lc   !!lignin factor
   real                                       :: rel_soil_moist
   real                                       :: sum_soil_energy
   real                                       :: sum_soil_water
   real                                       :: sum_soil_hcap
   real                                       :: sum_soil_slmsts
   real                                       :: sum_soil_soilcp
   real                                       :: avg_soil_temp
   real                                       :: avg_soil_fliq
   real                                       :: layer_soil_energy
   real                                       :: layer_soil_water
   real                                       :: layer_soil_hcap
   real                                       :: layer_soil_slmsts
   real                                       :: layer_soil_soilcp
   !----- External functions. -------------------------------------------------------------!
   real                          , external   :: het_resp_weight
   real                          , external   :: root_resp_norm
   !----- Locally saved variables. --------------------------------------------------------!
   real                          , save       :: dtlsm_o_frqsum
   logical                       , save       :: first_time = .true.
   !---------------------------------------------------------------------------------------!


   !----- Assign the constant scaling factor. ---------------------------------------------!
   if (first_time) then
      first_time     = .false.
      dtlsm_o_frqsum = dtlsm / frqsum
   end if
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Find the current root respiration.  This is done cohort by cohort because the    !
   ! parameters may be different depending on the PFT, and also because each layer has a   !
   ! different temperature.                                                                !
   !---------------------------------------------------------------------------------------!
   cpatch => csite%patch(ipa)
   do ico = 1,cpatch%ncohorts
      ipft  = cpatch%pft(ico)
      kroot = cpatch%krdepth(ico)

      !------------------------------------------------------------------------------------!
      !    Add "intensive" contribution of each layer, assuming that the roots are equally !
      ! spread throughout the entire depth.                                                !
      !------------------------------------------------------------------------------------!
      cpatch%root_respiration(ico) = 0.0
      do k = kroot,mzg
         cpatch%root_respiration(ico) = cpatch%root_respiration(ico)                       &
                                      + root_resp_norm(ipft,csite%soil_tempk(k,ipa))       &
                                      * dslz(k)
      end do
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !      Now we make the value in umol/m2/s, by dividing by the total depth and        !
      ! multiplying by the total root biomass.  The minus sign is because slz is negative. !
      !------------------------------------------------------------------------------------!
      cpatch%root_respiration(ico) = - cpatch%root_respiration(ico) * cpatch%broot(ico)    &
                                     * cpatch%nplant(ico) / slz(kroot)
      !------------------------------------------------------------------------------------!


      !----- Add this time step to the daily mean root respiration. -----------------------!
      cpatch%today_root_resp(ico)  = cpatch%today_root_resp(ico)                           &
                                   + cpatch%root_respiration(ico)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     The following is for output only, we switch the units to kgC/plant/yr.         !
      !------------------------------------------------------------------------------------!
      cpatch%fmean_root_resp(ico)   = cpatch%fmean_root_resp (ico)                         &
                                    + cpatch%root_respiration(ico)                         &
                                    * umols_2_kgCyr * dtlsm_o_frqsum                       &
                                    / cpatch%nplant          (ico)
      !------------------------------------------------------------------------------------!
   end do
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Integrate the soil extensive properties, plus the minimum and maximum possible    !
   ! soil water content of the active layer.                                               !
   !---------------------------------------------------------------------------------------!
   sum_soil_energy = 0.0
   sum_soil_hcap   = 0.0
   sum_soil_water  = 0.0
   sum_soil_slmsts = 0.0
   sum_soil_soilcp = 0.0
   do k = k_rh_active,mzg
      nsoil = ntext_soil(k)

      !------------------------------------------------------------------------------------!
      !    Convert the units so energy is in J/m2, heat capacity in J/m2/K, and water in   !
      ! kg/m2.                                                                             !
      !------------------------------------------------------------------------------------!
      sum_soil_energy = sum_soil_energy + csite%soil_energy(k,ipa)        * dslz(k)
      sum_soil_hcap   = sum_soil_hcap   + soil(nsoil)%slcpd               * dslz(k)
      sum_soil_water  = sum_soil_water  + csite%soil_water (k,ipa) * wdns * dslz(k)
      sum_soil_slmsts = sum_soil_slmsts + soil(nsoil)%slmsts       * wdns * dslz(k)
      sum_soil_soilcp = sum_soil_soilcp + soil(nsoil)%soilcp       * wdns * dslz(k)
      !------------------------------------------------------------------------------------!
   end do
   !---------------------------------------------------------------------------------------!


   !EJL This is the section that calculates respiration. Need to do this for
   !each vertical layer
   do k=1,nzl 
     nsoil = ntext_soil(k)
     
     !---Get values from above, soil energy, heat capacity, moisture, etc. for each layer   !
     !   instead of for the whole root zone                                                 !
     layer_soil_hcap = soil(nsoil)%slcpd               * dslz(k)
     layer_soil_water = csite%soil_water(k,ipa) * wdns * dslz(k)
     layer_soil_energy = csite%soil_energy(k,ipa)        * dslz(k)
     layer_soil_slmsts = soil(nsoil)%slmsts       * wdns * dslz(k)
     layer_soil_soilcp = soil(nsoil)%soilcp       * wdns * dslz(k)

     !----- Find the average temperature and the relative soil moisture. -------------------!
     select case (ivertresp)
       case (0)
       call uextcm2tl(sum_soil_energy,sum_soil_water,layer_soil_hcap,avg_soil_temp,avg_soil_fliq)
       case (1)
       call uextcm2tl(layer_soil_energy,layer_soil_water,layer_soil_hcap,avg_soil_temp,avg_soil_fliq)
     end select



     rel_soil_moist = min( 1.0, max(0.0, ( layer_soil_water  - layer_soil_soilcp )              &
                                     / ( layer_soil_slmsts - layer_soil_soilcp ) ) )
     !--------------------------------------------------------------------------------------!
     !----- Compute soil/temperature modulation of heterotrophic respiration. --------------!
     csite%A_decomp(k,ipa) = het_resp_weight(avg_soil_temp,rel_soil_moist)
     !--------------------------------------------------------------------------------------!

     !----- Compute nitrogen immobilization factor. ----------------------------------------!
     call resp_f_decomp(csite,ipa,k,Lc)
     !--------------------------------------------------------------------------------------!

     !----- Compute heterotrophic respiration. ---------------------------------------------!
     call resp_rh(csite,ipa,k,Lc)
     !--------------------------------------------------------------------------------------!

     !----- Update averaged variables. -----------------------------------------------------!
     csite%today_A_decomp (k,ipa) = csite%today_A_decomp(k,ipa) + csite%A_decomp(k,ipa)
     csite%today_Af_decomp(k,ipa) = csite%today_Af_decomp(k,ipa)                            &
                              + csite%A_decomp(k,ipa) * csite%f_decomp(k,ipa)
     !--------------------------------------------------------------------------------------!

     !----- The output is converted to kgC/m2/yr. ------------------------------------------!
     csite%fmean_rh(k,ipa) = csite%fmean_rh(k,ipa)                                              &
                           + csite%rh(k,ipa) * umols_2_kgCyr * dtlsm_o_frqsum
     csite%fmean_cwd_rh(ipa) = csite%fmean_cwd_rh(ipa)                                      &
                           + csite%cwd_rh(ipa) * umols_2_kgCyr * dtlsm_o_frqsum
     !--------------------------------------------------------------------------------------!
   end do ! vertical loop


   return
end subroutine soil_respiration
!==========================================================================================!
!==========================================================================================!







!==========================================================================================!
!==========================================================================================!
!     This subroutine computes the Nitrogen immobilization factor.                         !
!------------------------------------------------------------------------------------------!
subroutine resp_f_decomp(csite,ipa,k,Lc)

   use ed_state_vars, only : sitetype               ! ! structure
   use decomp_coms  , only : r_stsc                 & ! intent(in)
                           , N_immobil_supply_scale & ! intent(in)
                           , decay_rate_stsc        & ! intent(in)
                           , n_decomp_lim           ! ! intent(in)
   use pft_coms     , only : c2n_structural         & ! intent(in)
                           , c2n_slow               ! ! intent(in)

   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   type(sitetype), target      :: csite
   integer       , intent(in)  :: ipa
   integer       , intent(in)  :: k
   real          , intent(out) :: Lc
   !----- Local variables. ----------------------------------------------------------------!
   real                        :: N_immobilization_demand
   !---------------------------------------------------------------------------------------!

 
   if (csite%structural_soil_C(k,ipa) > 0.0) then
      if (csite%structural_soil_L(k,ipa) == csite%structural_soil_C(k,ipa)) then
         Lc = 0.049787 ! = exp(-3.0)
      else
         Lc = exp(-3.0 * csite%structural_soil_L(k,ipa)/csite%structural_soil_C(k,ipa))
      end if
   else
      Lc=0.0
   end if
   
   if (n_decomp_lim == 1) then
      N_immobilization_demand = csite%A_decomp(k,ipa) * Lc * decay_rate_stsc               &
                              * csite%structural_soil_C(k,ipa)                             &
                              * ((1.0 - r_stsc) / c2n_slow - 1.0 / c2n_structural)
      
      csite%f_decomp(k,ipa)     = N_immobil_supply_scale * csite%mineralized_soil_N(k,ipa) &
                              / ( N_immobilization_demand                                  &
                                + N_immobil_supply_scale  * csite%mineralized_soil_N(k,ipa))
   else
      !----- Option for no plant N limitation. --------------------------------------------!
      csite%f_decomp(k,ipa)     = 1.0
   end if

   return
end subroutine resp_f_decomp
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine computes the heterotrophic respiration.                              !
!------------------------------------------------------------------------------------------!
subroutine resp_rh(csite,ipa,k,Lc)

   use ed_state_vars, only : sitetype        ! ! structure
   use consts_coms  , only : kgCday_2_umols  ! ! intent(in)
   use decomp_coms  , only : decay_rate_stsc & ! intent(in)
                           , decay_rate_fsc  & ! intent(in)
                           , decay_rate_ssc  & ! intent(in)
                           , r_fsc           & ! intent(in)
                           , r_ssc           & ! intent(in)
                           , r_stsc          & ! intent(in)
                           , cwd_frac        ! ! intent(in)

   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   type(sitetype), target       :: csite
   integer       , intent(in)   :: ipa
   integer       , intent(in)   :: k
   real          , intent(in)   :: Lc
   !----- Local variables. ----------------------------------------------------------------!
   real                         :: fast_C_loss
   real                         :: structural_C_loss
   real                         :: slow_C_loss
   !---------------------------------------------------------------------------------------!



   !----- The following variables have units of [umol_CO2/m2/s]. --------------------------!
   ! EJL - make rh a function of depth
   fast_C_loss       = kgCday_2_umols * csite%A_decomp(k,ipa)                                &
                     * decay_rate_fsc * csite%fast_soil_C(k,ipa)
   structural_C_loss = kgCday_2_umols * csite%A_decomp(k,ipa) * Lc * decay_rate_stsc         &
                     * csite%structural_soil_C(k,ipa)* csite%f_decomp(k,ipa)
   slow_C_loss       = kgCday_2_umols * csite%A_decomp(k,ipa)                                &
                     * decay_rate_ssc * csite%slow_soil_C(k,ipa)
   !---------------------------------------------------------------------------------------!

   !----- Find the heterotrophic respiration and the fraction due to CWD. -----------------!
   csite%rh(k,ipa)     = r_fsc * fast_C_loss + r_stsc * structural_C_loss                    &
                     + r_ssc * slow_C_loss
   csite%cwd_rh(ipa) = cwd_frac * (r_stsc * structural_C_loss + r_ssc * slow_C_loss)
   ! csite%cwd_rh(ipa) = r_stsc * structural_C_loss
   !---------------------------------------------------------------------------------------!

   return
end subroutine resp_rh
!==========================================================================================!
!==========================================================================================!


!==========================================================================================!
!==========================================================================================!
!     This subroutine computes the depth of the organic layers and redistributes
!     the soil pools vertically. Should be called daily from veg dynamics. Maybe
!     move to update_soil_CN? 
!     !EJL 
!------------------------------------------------------------------------------------------!
subroutine organic_layer_depth(cgrid)

   use ed_state_vars, only : edtype        & ! structure
                           , polygontype   &
                           , sitetype      !
   use grid_coms,     only : nzg, nzl      !
   use soil_coms,     only : olz, dolz     ! 
   use ed_misc_coms,  only : icarbdyn      & ! intent(in) 
                           , isoiltext     ! intent(in) 
   use decomp_coms,   only : organic_soil_texture ! intent(in)

   implicit none
   !----- Arguments. ---------------------------------------------------------------------!
   type(edtype), target       :: cgrid
   !----- Local variables.    ------------------------------------------------------------!
   type(polygontype), pointer   :: cpoly
   type(sitetype)   , pointer   :: csite
   real                         :: fast_c_den
   real                         :: slow_c_den
   real                         :: struct_c_den
   real                         :: fillfrac
   real                         :: fillcheck
   real                         :: ld
   real                         :: nextld
   real                         :: edep
   real                         :: mdep
   real                         :: extc1
   real                         :: extc2
   real                         :: extc3
   real                         :: ld2 
   real                         :: c1d
   real                         :: c2d
   real                         :: c3d
   real                         :: n1d
   real                         :: n2d
   real                         :: n3d
   real    , dimension(nzl)     :: oldc1
   real    , dimension(nzl)     :: oldc2
   real    , dimension(nzl)     :: oldc3
   real    , dimension(nzl)     :: newc1
   real    , dimension(nzl)     :: newc2
   real    , dimension(nzl)     :: newc3
   integer                      :: k
   integer                      :: ipy
   integer                      :: isi
   integer                      :: ipa
   integer                      :: count
   real                         :: total_peat
   !---------------------------------------------------------------------------------------!

   !Carbon density of different soil carbon pools (kg/m3). Bulk density * carbon
   !content from Ise et al. (their source is Yi et al. 2009)
   fast_c_den = 70. * 0.4152
   slow_c_den = 179. * 0.3278
   struct_c_den = 124. * 0.37  ! no values, so just used mean of lit and hum

   polygonloop: do ipy = 1,cgrid%npolygons

      cpoly => cgrid%polygon(ipy)

      siteloop: do isi = 1,cpoly%nsites

         csite => cpoly%site(isi)

         patchloop: do ipa = 1,csite%npatches


           total_peat = 0.0
           fillfrac = 0.0
           ! Find total depth of organic layer
           do k=1,nzl
             total_peat = total_peat + csite%slow_soil_C(k,ipa) / slow_c_den   &
             + csite%fast_soil_C(k,ipa) /fast_c_den                            &
             + csite%structural_soil_c(k,ipa) / struct_c_den
           end do   
           csite%peat_depth(ipa) = total_peat

           ! Find the fill fraction of the layers 
           csite%litter_depth(:,ipa) = 0.0
           do k = 1, nzl
!              print*, olz(k), olz(k+1), dolz(k),fillfrac
             if (csite%peat_depth(ipa) .ge. (-1.0 * olz(k))) then 
               csite%litter_depth(k,ipa) = 1.0
             else if (csite%peat_depth(ipa) .lt. (-1.0*olz(k)) .and. &
               csite%peat_depth(ipa) .gt. (-1.0*olz(k+1))) then 
               csite%litter_depth(k,ipa) = (csite%peat_depth(ipa) + olz(k+1)) / dolz(k)
             endif
             fillfrac=fillfrac+csite%litter_depth(k,ipa)

             select case (isoiltext)
               case (1)
                 if (csite%litter_depth(k,ipa) .gt. 0.51) then 
                   cpoly%ntext_soil(k,isi) = organic_soil_texture
                 else
                 cpoly%ntext_soil(k,isi) = cpoly%ntext_soil(1,isi)
                 end if
               case (0)
             end select
           end do


           !Redistribute soil pools based on new depth starting from surface (k=nzl).
           !Option 1) We will maintain fraciton of each pool in each layer when pushing 
           ! up or down. i.e. if top layer is 70% met and 30% struct and has 10%
           ! extra carbon, then 10% of depth of met and depth of struct get
           ! pushed down.
           !Option 2) Metabolic carbon gets preferentially moved up while humic
           !then structural gets preferentially moved down. 

            fillcheck=0.0 ! Zero out a bunch of values and arrays.
!           ld = 0.    
!           nextld = 0.
!           extc1=0.
            oldc1(:)=csite%fast_soil_C(:,ipa)
            oldc2(:)=csite%structural_soil_C(:,ipa)
            oldc3(:)=csite%slow_soil_C(:,ipa)
            newc1(:)=0.0
            newc2(:)=0.0
            newc3(:)=0.0
          
           !debugging print statements
!           print*, 'dolz', dolz
!           print*, 'patch,carbon pools',ipa, oldc1, oldc2, oldc3 
   
            count=0 
            select case (icarbdyn)
!!!!!!!!!!!OPTION 1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!
             case (1)  
              do while (abs(fillfrac-fillcheck) > 0.01)
                 fillcheck=0.0
                 if (count > 99) then
                     print*,'fillfrac,fillcheck ' 
                     print*,fillfrac,fillcheck
                     call fatal_error('organic_layer_depth (opt 1) did not converge' &
                     ,'organic_layer_depth','soil_respiration.f90')
                 end if
!                 print*, 'while fillfrac ne fillcheck',fillfrac, fillcheck
                 do k=nzl,2,-1 
                    ld = oldc1(k) / fast_c_den &
                       + oldc2(k) / struct_c_den &
                       + oldc3(k) / slow_c_den
                    nextld = oldc1(k-1) / fast_c_den &
                       + oldc2(k-1) / struct_c_den &
                       + oldc3(k-1) /slow_c_den
                    if (ld .ge. dolz(k)) then ! If extra carbon, then new levels
                    ! equal to layer thickness and extra is difference
                       newc1(k)=oldc1(k)*dolz(k) / ld
                       newc2(k)=oldc2(k)*dolz(k) / ld
                       newc3(k)=oldc3(k)*dolz(k) / ld

                       extc1 = oldc1(k) - newc1(k)                       
                       extc2 = oldc2(k) - newc2(k)                       
                       extc3 = oldc3(k) - newc3(k)                       

                       oldc1(k-1)=oldc1(k-1) + extc1
                       oldc2(k-1)=oldc2(k-1) + extc2
                       oldc3(k-1)=oldc3(k-1) + extc3

                    else ! missing carbon
                       mdep = min(dolz(k)-ld, nextld)
                  
                       if (nextld .gt. 0.0) then
                         extc1 = max(oldc1(k-1)*mdep/nextld, 0.0)
                         extc2 = max(oldc2(k-1)*mdep/nextld, 0.0)
                         extc3 = max(oldc3(k-1)*mdep/nextld, 0.0)
                       else
                         extc1 = 0.0
                         extc2 = 0.0
                         extc3 = 0.0
                       endif

                       newc1(k) = oldc1(k) + extc1
                       newc2(k) = oldc2(k) + extc2
                       newc3(k) = oldc3(k) + extc3

                       oldc1(k-1) = oldc1(k-1) - extc1
                       oldc2(k-1) = oldc2(k-1) - extc2
                       oldc3(k-1) = oldc3(k-1) - extc3
                    endif !extra carbon?
                    
                    !Test that the algorithm created correct layer thickness
                    ld2 = newc1(k) / fast_c_den + newc2(k) / struct_c_den &
                        + newc3(k) / slow_c_den

                    fillcheck = fillcheck+ld2/dolz(k)       
                 end do ! vert loop

                 !!! ASSUME DEEPEST LAYER NEVER FILLS UP. IF THIS FAILS, CREATED
                 !!! THICKER LAYERS           
                 newc1(1) = oldc1(1)
                 newc2(1) = oldc2(1)
                 newc3(1) = oldc3(1)

                 ld2 = newc1(1) / fast_c_den + newc2(1) / struct_c_den &
                     + newc3(1) / slow_c_den

                 fillcheck = fillcheck+ld2/dolz(1)       

                 oldc1(:)=newc1(:)
                 oldc2(:)=newc2(:)
                 oldc3(:)=newc3(:)
                 count=count+1
              end do ! while loop 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! OPTION 2!!!!!!!!!!!!!!!!!!!!!!!!
             case (2)
              do while (abs(fillfrac-fillcheck) .gt. 0.01)
                 fillcheck=0.0
                 if (count .gt. 99) then 
                    print*,'fillfrac,fillcheck '         
                    print*,fillfrac,fillcheck
                    call fatal_error('organic_layer_depth (opt2) did not converge' &
                    ,'organic_layer_depth','soil_respiration.f90')
                 end if
                 do k=nzl,2,-1 
                    ld = oldc1(k) / fast_c_den &
                       + oldc2(k) / struct_c_den &
                       + oldc3(k) / slow_c_den
                    nextld = oldc1(k-1) / fast_c_den &
                       + oldc2(k-1) / struct_c_den &
                       + oldc3(k-1) /slow_c_den

                    c1d = oldc1(k) / fast_c_den
                    c2d = oldc2(k) / struct_c_den
                    c3d = oldc3(k) / slow_c_den
                    n1d = oldc1(k-1) / fast_c_den
                    n2d = oldc2(k-1) / struct_c_den
                    n3d = oldc3(k-1) / slow_c_den

                    if (ld .ge. dolz(k)) then ! if extra carbon
                      edep = ld - dolz(k) ! extra depth to be pushed down
                      if (edep .ge. c3d) then
                         extc3 = oldc3(k)
                         if (edep .ge. c3d+c2d) then
                            extc2 = oldc2(k)
                            extc1 = (edep - c3d -c2d) * fast_c_den
                         else
                            extc1 = 0.0
                            extc2 = (edep-c3d) * struct_c_den
                         end if
                      else
                         extc1 = 0.0
                         extc2 = 0.0
                         extc3 = edep*slow_c_den
                      end if                 

                       newc1(k)=oldc1(k)-extc1
                       newc2(k)=oldc2(k)-extc2
                       newc3(k)=oldc3(k)-extc3

                       oldc1(k-1)=oldc1(k-1) + extc1
                       oldc2(k-1)=oldc2(k-1) + extc2
                       oldc3(k-1)=oldc3(k-1) + extc3

                    else ! missing carbon
                       mdep = min(dolz(k)-ld, nextld)
                       if (mdep .ge. n1d) then
                          extc1 = oldc1(k-1)
                          if (mdep .ge. n1d+n2d) then
                             extc2 = oldc2(k-1)
                             extc3 = (mdep - n1d - n2d) * slow_c_den
                          else
                             extc3 = 0.0
                             extc2 = (mdep-n1d) * struct_c_den
                          endif 
                       else 
                          extc1 = mdep * fast_c_den
                          extc2 = 0.0
                          extc3 = 0.0
                       endif
 
                       newc1(k) = oldc1(k) + extc1
                       newc2(k) = oldc2(k) + extc2
                       newc3(k) = oldc3(k) + extc3

                       oldc1(k-1) = oldc1(k-1) - extc1
                       oldc2(k-1) = oldc2(k-1) - extc2
                       oldc3(k-1) = oldc3(k-1) - extc3
                    endif !extra carbon?
                    
                    !Test that the algorithm created correct layer thickness
                    ld2 = newc1(k) / fast_c_den + newc2(k) / struct_c_den &
                        + newc3(k) / slow_c_den

                    fillcheck = fillcheck+ld2/dolz(k)       
                 end do ! vert loop

                 !!! ASSUME DEEPEST LAYER NEVER FILLS UP. IF THIS FAILS, CREATE
                 !!! THICKER LAYERS           
                 newc1(1) = oldc1(1)
                 newc2(1) = oldc2(1)
                 newc3(1) = oldc3(1)

                 ld2 = newc1(1) / fast_c_den + newc2(1) / struct_c_den &
                     + newc3(1) / slow_c_den

                 fillcheck = fillcheck+ld2/dolz(1)       

                 oldc1(:)=newc1(:)
                 oldc2(:)=newc2(:)
                 oldc3(:)=newc3(:)

                 print*, count, fillcheck, fillfrac
                 count=count+1
              end do ! while loop
 
           end select !option of redistribution



           do k=1,nzl 
             csite%fast_soil_C(k,ipa) = max(oldc1(k), 0.0)
             csite%structural_soil_C(k,ipa) = max(oldc2(k),0.0)
             csite%slow_soil_C(k,ipa) = max(oldc3(k),0.0)
           end do

         end do patchloop
      end do siteloop
   end do polygonloop


   return
end subroutine organic_layer_depth
!==========================================================================================!
!==========================================================================================!




!==========================================================================================!
!==========================================================================================!
!     This subroutine will update the soil carbon and nitrogen pools.                      !
!------------------------------------------------------------------------------------------!
subroutine update_C_and_N_pools(cgrid)
   
   use ed_state_vars, only : edtype          & ! structure
                           , polygontype     & ! structure
                           , sitetype        ! ! structure
   use decomp_coms  , only : decay_rate_fsc  & ! intent(in)
                           , decay_rate_stsc & ! intent(in)
                           , decay_rate_ssc  & ! intent(in)
                           , r_stsc          ! ! intent(in)
   use pft_coms     , only : c2n_slow        & ! intent(in)
                           , c2n_structural  ! ! intent(in)
   use grid_coms    , only : nzl             ! ! intent(in)

   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   type(edtype)     , target   :: cgrid
   !----- Local variables. ----------------------------------------------------------------!
   type(polygontype), pointer  :: cpoly
   type(sitetype)   , pointer  :: csite
   integer                     :: ipy
   integer                     :: isi
   integer                     :: ipa
   integer                     :: k
   real     , dimension(nzl)   :: Lc
   real     , dimension(nzl)   :: fast_C_loss
   real     , dimension(nzl)   :: fast_N_loss
   real     , dimension(nzl)   :: structural_C_loss
   real     , dimension(nzl)   :: structural_L_loss
   real     , dimension(nzl)   :: slow_C_input
   real     , dimension(nzl)   :: slow_C_loss
   !---------------------------------------------------------------------------------------!

   polygonloop: do ipy = 1,cgrid%npolygons

      cpoly => cgrid%polygon(ipy)

      siteloop: do isi = 1,cpoly%nsites
         
         csite => cpoly%site(isi)

         patchloop: do ipa = 1,csite%npatches

            do k=1,nzl

              if (csite%structural_soil_C(k,ipa) > 0.0) then
                 if (csite%structural_soil_L(k,ipa) == csite%structural_soil_C(k,ipa)) then
                    Lc(k) = 0.049787 ! = exp(-3.0)
                 else
                    Lc(k) = exp( -3.0 * csite%structural_soil_L(k,ipa)                     &
                          /  csite%structural_soil_C(k,ipa))
                 end if
              else
                 Lc(k)=0.0
              end if
      
              !----- Fast pools. ----------------------------------------------------------!
              fast_C_loss(k) = csite%today_A_decomp(k,ipa) * decay_rate_fsc                &
                        * csite%fast_soil_C(k,ipa)
              fast_N_loss(k) = csite%today_A_decomp(k,ipa) * decay_rate_fsc                &
                        * csite%fast_soil_N(k,ipa)

              !----- Structural pools. ----------------------------------------------------!
             structural_C_loss(k) = csite%today_Af_decomp(k,ipa) * Lc(k) * decay_rate_stsc &
                              * csite%structural_soil_C(k,ipa)
             structural_L_loss(k) = csite%today_Af_decomp(k,ipa) * Lc(k) * decay_rate_stsc &
                              * csite%structural_soil_L(k,ipa)

              !----- Slow pools. ----------------------------------------------------------!
              slow_C_input(k) = (1.0 - r_stsc) * structural_C_loss(k)
              slow_C_loss(k)  = csite%today_A_decomp(k,ipa) * decay_rate_ssc               &
                         * csite%slow_soil_C(k,ipa)
            
              !----- Mineralized pool. ----------------------------------------------------!
              csite%mineralized_N_input(k,ipa) = fast_N_loss(k) + slow_C_loss(k) / c2n_slow
              csite%mineralized_N_loss(k,ipa)  = csite%total_plant_nitrogen_uptake(ipa)        &
                                 + csite%today_Af_decomp(k,ipa) * Lc(k) * decay_rate_stsc  &
                                 * csite%structural_soil_C(k,ipa)                          &
                                 * ( (1.0 - r_stsc) / c2n_slow - 1.0 / c2n_structural)


              !----------------------------------------------------------------------------!
              !      All carbon fluxes have units kgC/m2/day, and we are updating on the   !
              ! daily time step. Nitrogen has units kgN/m2/day. Updating with vertical res.!
              ! Only adding inputs to top organic layer (k==nzl). EJL                      !
              !----------------------------------------------------------------------------!
              if (k == nzl) then 
                csite%fast_soil_C(k,ipa)   = csite%fast_soil_C(k,ipa) + csite%fsc_in(ipa)  &
                                         - fast_C_loss(k)

                csite%fast_soil_N(k,ipa)   = csite%fast_soil_N(k,ipa) + csite%fsn_in(ipa)  &
                                         - fast_N_loss(k)

                csite%structural_soil_C(k,ipa) = csite%structural_soil_C(k,ipa)              &
                                         + csite%ssc_in(ipa) - structural_C_loss(k)
                csite%structural_soil_L(k,ipa) = csite%structural_soil_L(k,ipa)              &
                                         + csite%ssl_in(ipa) - structural_L_loss(k)
              else
                csite%fast_soil_C(k,ipa)   = csite%fast_soil_C(k,ipa)                      &
                                         - fast_C_loss(k)

                csite%fast_soil_N(k,ipa)   = csite%fast_soil_N(k,ipa)                      &
                                         - fast_N_loss(k)
                csite%structural_soil_C(k,ipa) = csite%structural_soil_C(k,ipa)              &
                                         - structural_C_loss(k)
                csite%structural_soil_L(k,ipa) = csite%structural_soil_L(k,ipa)              &
                                         - structural_L_loss(k)
              end if

              csite%slow_soil_C(k,ipa)       = csite%slow_soil_C(k,ipa) + slow_C_input(k)  &
                                         - slow_C_loss(k)
            
              csite%mineralized_soil_N(k,ipa) = csite%mineralized_soil_N(k,ipa)            &
                                          + csite%mineralized_N_input(k,ipa)               &
                                          - csite%mineralized_N_loss(k,ipa)

!              call organic_layer_depth(cgrid)

              !----------------------------------------------------------------------------!
              !      Force all pools to be either zero or positive.                        !
              !----------------------------------------------------------------------------!
              csite%fast_soil_C(k,ipa)        = max(0.0,csite%fast_soil_C(k,ipa))
              csite%structural_soil_C(k,ipa)  = max(0.0,csite%structural_soil_C(k,ipa))
              csite%structural_soil_L(k,ipa)  = max(0.0,csite%structural_soil_L(k,ipa))
              csite%slow_soil_C(k,ipa)        = max(0.0,csite%slow_soil_C(k,ipa))
              csite%fast_soil_N(k,ipa)        = max(0.0,csite%fast_soil_N(k,ipa))
              csite%mineralized_soil_N(k,ipa) = max(0.0,csite%mineralized_soil_N(k,ipa))
            
           end do !organic layer

         end do patchloop
      end do siteloop
   end do polygonloop
   
   return
end subroutine update_C_and_N_pools
!==========================================================================================!
!==========================================================================================!


!==========================================================================================!
!==========================================================================================!
!     This subroutine computes the anoxic heterotrophic respiration following Walter and   !
!  Heinmann (2000) and Ise et al. 2010 (Climate Change and Variability)                    !
!------------------------------------------------------------------------------------------!
subroutine anoxic_resp(csite,ipa,k,Lc)

   use ed_state_vars, only : sitetype        ! ! structure
   use consts_coms  , only : kgCday_2_umols  ! ! intent(in)
   use decomp_coms  , only : decay_rate_stsc & ! intent(in)
                           , decay_rate_fsc  & ! intent(in)
                           , decay_rate_ssc  & ! intent(in)
                           , r_fsc           & ! intent(in)
                           , r_ssc           & ! intent(in)
                           , r_stsc          & ! intent(in)
                           , cwd_frac        ! ! intent(in)

   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   type(sitetype), target       :: csite
   integer       , intent(in)   :: ipa
   integer       , intent(in)   :: k
   real          , intent(in)   :: Lc
   !----- Local variables. ----------------------------------------------------------------!
   real                         :: fast_C_loss
   real                         :: structural_C_loss
   real                         :: slow_C_loss
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   ! Production of CH4 fol


   return
end subroutine anoxic_resp
!==========================================================================================!
!==========================================================================================!

end module soil_respiration_module







!==========================================================================================!
!==========================================================================================!
!     This function determines the normalised root respiration (umol/kgC_fine_root/s) at a !
! given soil layer.                                                                        !
!------------------------------------------------------------------------------------------!
real function root_resp_norm(ipft,soil_temp)
   use pft_coms       , only : root_respiration_factor  & ! intent(in)
                             , rrf_low_temp             & ! intent(in)
                             , rrf_high_temp            & ! intent(in)
                             , rrf_decay_e              & ! intent(in)
                             , rrf_hor                  & ! intent(in)
                             , rrf_q10                  ! ! intent(in)
   use farq_leuning   , only : arrhenius                & ! function
                             , collatz                  ! ! function
   use rk4_coms       , only : tiny_offset              ! ! intent(in)
   use physiology_coms, only : iphysiol                 ! ! intent(in)
   use consts_coms    , only : lnexp_min8               & ! intent(in)
                             , lnexp_max8               & ! intent(in)
                             , t008                     ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer     , intent(in) :: ipft
   real(kind=4), intent(in) :: soil_temp
   !----- Local variables. ----------------------------------------------------------------!
   real(kind=8)             :: soil_temp8
   real(kind=8)             :: rrf08
   real(kind=8)             :: rrf_low_temp8
   real(kind=8)             :: rrf_high_temp8
   real(kind=8)             :: rrf_decay_e8
   real(kind=8)             :: rrf_hor8
   real(kind=8)             :: rrf_q108
   real(kind=8)             :: lnexplow
   real(kind=8)             :: lnexphigh
   real(kind=8)             :: tlow_fun
   real(kind=8)             :: thigh_fun
   real(kind=8)             :: rrf8
   !----- External functions. -------------------------------------------------------------!
   real(kind=4)             :: sngloff
   !---------------------------------------------------------------------------------------!

   !----- Copy some variables to double precision temporaries. ----------------------------!
   soil_temp8      = dble(soil_temp                    )
   rrf08           = dble(root_respiration_factor(ipft))
   rrf_low_temp8   = dble(rrf_low_temp           (ipft)) + t008
   rrf_high_temp8  = dble(rrf_high_temp          (ipft)) + t008
   rrf_decay_e8    = dble(rrf_decay_e            (ipft))
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    Compute the functions that will control the Rrf function for low and high temper-  !
   ! ature.  In order to avoid floating point exceptions, we check whether the temperature !
   ! will make the exponential too large or too small.                                     !
   !---------------------------------------------------------------------------------------!
   !----- Low temperature. ----------------------------------------------------------------!
   lnexplow  = rrf_decay_e8 * (rrf_low_temp8  - soil_temp8)
   lnexplow  = max(lnexp_min8,min(lnexp_max8,lnexplow))
   tlow_fun  = 1.d0 +  exp(lnexplow)
   !----- High temperature. ---------------------------------------------------------------!
   lnexphigh = rrf_decay_e8 * (soil_temp8 - rrf_high_temp8)
   lnexphigh = max(lnexp_min8,min(lnexp_max8,lnexphigh))
   thigh_fun = 1.d0 + exp(lnexphigh)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    Decide which functional form to use based on the physiology.  This is just to make !
   ! it look similar to the leaf respiration respiration.                                  !
   !---------------------------------------------------------------------------------------!
   select case (iphysiol)
   case (0,1)
      rrf_hor8 = dble(rrf_hor(ipft))
      rrf8     = arrhenius(soil_temp8,rrf08,rrf_hor8) / (tlow_fun * thigh_fun)
   case (2,3)
      rrf_q108 = dble(rrf_q10(ipft))
      rrf8     = collatz(soil_temp8,rrf08,rrf_q108)   / (tlow_fun * thigh_fun)
   end select
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Convert result to single precision.                                               !
   !---------------------------------------------------------------------------------------!
   root_resp_norm = sngloff(rrf8,tiny_offset)
   !---------------------------------------------------------------------------------------!


   return
end function root_resp_norm
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This function computes the heterotrophic respiration limitation factor, which        !
! includes limitations due to temperature and soil moisture.                               !
!------------------------------------------------------------------------------------------!
real function het_resp_weight(soil_tempk,rel_soil_moist)

   use decomp_coms, only : resp_temperature_increase  & ! intent(in)
                         , resp_opt_water             & ! intent(in)
                         , resp_water_below_opt       & ! intent(in)
                         , resp_water_above_opt       & ! intent(in)
                         , decomp_scheme              & ! intent(in)
                         , rh_lloyd_1                 & ! intent(in)
                         , rh_lloyd_2                 & ! intent(in)
                         , rh_lloyd_3                 & ! intent(in)
                         , rh_decay_low               & ! intent(in)
                         , rh_decay_high              & ! intent(in)
                         , rh_low_temp                & ! intent(in)
                         , rh_high_temp               & ! intent(in)
                         , rh_decay_dry               & ! intent(in)
                         , rh_decay_wet               & ! intent(in)
                         , rh_dry_smoist              & ! intent(in)
                         , rh_wet_smoist              ! ! intent(in)
   use consts_coms, only : lnexp_min                  & ! intent(in)
                         , lnexp_max                  ! ! intent(in)

   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   real(kind=4), intent(in) :: soil_tempk
   real(kind=4), intent(in) :: rel_soil_moist
   !----- Local variables. ----------------------------------------------------------------!
   real(kind=4)             :: temperature_limitation
   real(kind=4)             :: water_limitation
   real(kind=4)             :: lnexplloyd
   real(kind=4)             :: lnexplow
   real(kind=4)             :: lnexphigh
   real(kind=4)             :: tlow_fun
   real(kind=4)             :: thigh_fun
   real(kind=4)             :: lnexpdry
   real(kind=4)             :: lnexpwet
   real(kind=4)             :: smdry_fun
   real(kind=4)             :: smwet_fun
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Find the temperature dependence.                                                  !
   !---------------------------------------------------------------------------------------!
   select case(decomp_scheme)
   case (0,3)
      !----- Use original ED-2.1 exponential temperature dependence. ----------------------!
      temperature_limitation = min( 1.0                                                    &
                                  , exp( resp_temperature_increase * (soil_tempk-318.15)))
      !------------------------------------------------------------------------------------!
   case (1,4) 
      !----- Use Lloyd and Taylor (1994) temperature dependence. --------------------------!
      lnexplloyd             = rh_lloyd_1 * ( rh_lloyd_2 - 1. / (soil_tempk - rh_lloyd_3))
      lnexplloyd             = max(lnexp_min,min(lnexp_max,lnexplloyd))
      temperature_limitation = min( 1.0, resp_temperature_increase * exp(lnexplloyd) )
      !------------------------------------------------------------------------------------!
   case (2)
      !------------------------------------------------------------------------------------!
      !      Similar to the original ED-1.0 formulation, which is based on the CENTURY     !
      ! model.  The change in the functional form is to avoid power of negative numbers,   !
      ! but the coefficients were tuned to give a similar curve.                           !
      !------------------------------------------------------------------------------------!
      !----- Low temperature limitation. --------------------------------------------------!
      lnexplow               = rh_decay_low * (rh_low_temp - soil_tempk)
      lnexplow               = max(lnexp_min,min(lnexp_max,lnexplow))
      tlow_fun               = 1.0 + exp(lnexplow)
      !----- High temperature limitation. -------------------------------------------------!
      lnexphigh              = rh_decay_high * (soil_tempk - rh_high_temp)
      lnexphigh              = max(lnexp_min,min(lnexp_max,lnexphigh))
      thigh_fun              = 1.0 + exp(lnexphigh)
      !----- Temperature limitation is a combination of both. -----------------------------!
      temperature_limitation = 1.0 / (tlow_fun * thigh_fun )
      !------------------------------------------------------------------------------------!
   end select
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Find the limitation due to soil moisture.                                         !
   !---------------------------------------------------------------------------------------!
   select case (decomp_scheme)
   case (0,1)
      !----- ED-2.1 default, also used when decomp_scheme is 1. ---------------------------!
      if (rel_soil_moist <= resp_opt_water)then
         water_limitation = exp( (rel_soil_moist - resp_opt_water) * resp_water_below_opt)
      else
         water_limitation = exp( (resp_opt_water - rel_soil_moist) * resp_water_above_opt)
      end if
      !------------------------------------------------------------------------------------!
   case (2)
      !----- Dry soil limitation. ---------------------------------------------------------!
      lnexpdry         = rh_decay_dry * (rh_dry_smoist - rel_soil_moist)
      lnexpdry         = max(lnexp_min,min(lnexp_max,lnexpdry))
      smdry_fun        = 1.0 + exp(lnexpdry)
      !----- Wet soil limitation. ---------------------------------------------------------!
      lnexpwet         = rh_decay_wet * (rel_soil_moist - rh_wet_smoist)
      lnexpwet         = max(lnexp_min,min(lnexp_max,lnexpwet))
      smwet_fun        = 1.0 + exp(lnexpwet)
      !----- Soil moisture limitation is a combination of both. ---------------------------!
      water_limitation = 1.0 / (smdry_fun * smwet_fun)
      !------------------------------------------------------------------------------------!
   case (3,4)
      !----- From Jaclyn Matthes:
      !----- Empirical equation from meta-analysis in Moyano et al., Biogeosciences,2012 --!
      water_limitation = rel_soil_moist * 4.0893 - rel_soil_moist**2 * 3.1681 - 0.3195897
      !------------------------------------------------------------------------------------!
   end select
   !---------------------------------------------------------------------------------------!



   !----- Compute the weight, which is just the combination of both. ----------------------!
   het_resp_weight = temperature_limitation * water_limitation
   !---------------------------------------------------------------------------------------!

   return
end function het_resp_weight
!==========================================================================================!
!==========================================================================================!
