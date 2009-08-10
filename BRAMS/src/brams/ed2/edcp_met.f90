!==========================================================================================!
!==========================================================================================!
!     This subroutine copies various atmospheric fields that are needed by ED-2.           !
!------------------------------------------------------------------------------------------!
subroutine copy_atm2lsm(ifm,init)
   use ed_state_vars        , only : edgrid_g     & ! structure
                                   , edtype       & ! structure
                                   , polygontype  ! ! structure
   use mem_basic            , only : co2_on       & ! intent(in)
                                   , co2con       & ! intent(in)
                                   , basic_g      ! ! structure
   use mem_radiate          , only : radiate_g    ! ! structure
   use mem_cuparm           , only : cuparm_g     ! ! structure
   use mem_micro            , only : micro_g      ! ! structure
   use mem_grid             , only : grid_g       & ! structure
                                   , zt           & ! intent(in)
                                   , dzt          & ! intent(in)
                                   , zm           & ! intent(in)
                                   , if_adap      & ! intent(in)
                                   , jdim         ! ! intent(in)
   use node_mod             , only : master_num   & ! intent(in)
                                   , mmzp         & ! intent(in)
                                   , mmxp         & ! intent(in)
                                   , mmyp         & ! intent(in)
                                   , ia           & ! intent(in)
                                   , iz           & ! intent(in)
                                   , ja           & ! intent(in)
                                   , jz           & ! intent(in)
                                   , ia_1         & ! intent(in)
                                   , iz1          & ! intent(in)
                                   , ja_1         & ! intent(in)
                                   , jz1          ! ! intent(in)
   use rconstants           , only : cpi          & ! intent(in)
                                   , cp           & ! intent(in)
                                   , p00          & ! intent(in)
                                   , rocp         & ! intent(in)
                                   , rgas         & ! intent(in)
                                   , cliq         & ! intent(in)
                                   , alli         & ! intent(in)
                                   , cice         & ! intent(in)
                                   , t3ple        & ! intent(in)
                                   , cpor         & ! intent(in)
                                   , tsupercool   ! ! intent(in)
   use ed_node_coms         , only : mynum        ! ! intent(in)
   use canopy_radiation_coms, only : rlong_min    ! ! intent(in)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   integer                             , intent(in) :: ifm
   logical                             , intent(in) :: init
   !----- Local variables -----------------------------------------------------------------!
   type(edtype)                        , pointer    :: cgrid
   type(polygontype)                   , pointer    :: cpoly 
   integer                                          :: ipy,isi,ipa
   integer                                          :: m1,m2,m3,m1max
   integer                                          :: ix,iy,i,j,n
   integer                                          :: k1w,k2w,k3w
   integer                                          :: k2u,k3u,k2u_1,k3u_1
   integer                                          :: k2v,k3v,k2v_1,k3v_1
   real, dimension(mmxp(ifm),mmyp(ifm))             :: rshortd
   real, dimension(mmxp(ifm),mmyp(ifm))             :: up_mean,vp_mean,pi0_mean
   real, dimension(mmxp(ifm),mmyp(ifm))             :: rv_mean,rtp_mean,theta_mean
   real, dimension(mmxp(ifm),mmyp(ifm))             :: co2p_mean
   real, dimension(mmxp(ifm),mmyp(ifm))             :: map_2d_lsm
   real                                             :: rshort1,cosz1,rshortd1,scalar1
   real                                             :: topma_t,wtw,wtu1,wtu2,wtv1,wtv2
   !---------------------------------------------------------------------------------------!
  
   !----- Assigning some aliases. ---------------------------------------------------------!
   m2    =  mmxp(ifm)
   m3    =  mmyp(ifm)
   cgrid => edgrid_g(ifm)


   up_mean    = 0.0
   vp_mean    = 0.0
   pi0_mean   = 0.0
   rv_mean    = 0.0
   rtp_mean   = 0.0
   theta_mean = 0.0


   !---------------------------------------------------------------------------------------!
   !    Loop over the grid points, and prepare the atm data fields.  First, check which    !
   ! vertical coordinate we are using.                                                     !
   !---------------------------------------------------------------------------------------!
   select case (if_adap)
   case (0) 
      !------------------------------------------------------------------------------------!
      !      Terrain-following coordinate (sigma_z).  Here we simply use the lowest pre-   !
      ! dicted level for thermo and wind variables, and a simple average between levels    !
      ! one and two for Exner function and density.  I am not sure this is consistent with !
      ! what is done for density in adaptive coordinate...                                 !
      !------------------------------------------------------------------------------------!
      do j=ja,jz
         do i=ia,iz
            theta_mean(i,j) = basic_g(ifm)%theta(2,i,j)
            rv_mean(i,j)    = basic_g(ifm)%rv(2,i,j)
            rtp_mean(i,j)   = basic_g(ifm)%rtp(2,i,j)
            up_mean(i,j)    = (basic_g(ifm)%up(2,i,j) + basic_g(ifm)%up(2,i-1,j)) * 0.5
            vp_mean(i,j)    = (basic_g(ifm)%vp(2,i,j) + basic_g(ifm)% vp(2,i,j-jdim)) * 0.5
            pi0_mean(i,j)   = ( basic_g(ifm)%pp(1,i,j) + basic_g(ifm)%pp(2,i,j)            &
                              + basic_g(ifm)%pi0(1,i,j)+basic_g(ifm)%pi0(2,i,j)) * 0.5
            if (co2_on) then
               co2p_mean(i,j) = basic_g(ifm)%co2p(2,i,j)
            else
               co2p_mean(i,j) = co2con(1)
            end if
            !------------------------------------------------------------------------------!
            !      The following subroutine calculates diffuse shortwave radiation.  This  !
            ! is a gross approximation that uses constant scattering coefficients.  A more !
            ! appropriate way of doing this was probably already designed by DMM.  Valeriy !
            ! Ivanov is also investigating methods of parameterizing diffuse shortwave     !
            ! radiation.                                                                   !
            !------------------------------------------------------------------------------!
            rshort1 = radiate_g(ifm)%rshort(i,j)
            cosz1   = radiate_g(ifm)%cosz(i,j)
            call short2diff(rshort1,cosz1,rshortd1)
            rshortd(i,j) = rshortd1
         end do
      end do
   case (1)
      !------------------------------------------------------------------------------------!
      !     Adaptive height coordinate (shaved-eta).  Here we need to do a weighted aver-  !
      ! age using the lowest predicted points and the points avove them.                   !
      !------------------------------------------------------------------------------------!
      do j=ja,jz
         do i=ia,iz
            !------------------------------------------------------------------------------!
            !    Finding the lowest boundary levels, depending on which kind of variable   !
            ! we are averaging (staggered grid).                                           !
            !------------------------------------------------------------------------------!
            !----- U points, need to average between i-1 and i... -------------------------!
            k2u   = nint(grid_g(ifm)%flpu(i,j))
            k3u   = k2u + 1
            k2u_1 = nint(grid_g(ifm)%flpu(i-1,j))
            k3u_1 = k2u_1 + 1
            !----- V points, need to average between j-jdim and j... ----------------------!
            k2v   = nint(grid_g(ifm)%flpv(i,j))
            k3v   = k2v+1
            k2v_1 = nint(grid_g(ifm)%flpv(i,j-jdim))
            k3v_1 = k2v_1 + 1
            !----- W/T points, only the i,j points are needed...  -------------------------!
            k2w   = nint(grid_g(ifm)%flpw(i,j))
            k1w   = k2w - 1
            k3w   = k2w + 1
            !------------------------------------------------------------------------------!

            topma_t = .25 * ( grid_g(ifm)%topma(i,j) + grid_g(ifm)%topma(i-1,j)            &
                            + grid_g(ifm)%topma(i,j-jdim) + grid_g(ifm)%topma(i-1,j-jdim))

            !------------------------------------------------------------------------------!
            !     Computing the weights for lowest predicted points, relative to points    !
            ! above them.                                                                  !
            !------------------------------------------------------------------------------!
            wtw  = (zm(k2w) - topma_t) * dzt(k2w)
            wtu1 = grid_g(ifm)%aru(k2u_1,i-1,j)    / grid_g(ifm)%aru(k3u_1,i-1,j)
            wtu2 = grid_g(ifm)%aru(k2u,i,j)        / grid_g(ifm)%aru(k3u,i,j)
            wtv1 = grid_g(ifm)%arv(k2v_1,i,j-jdim) / grid_g(ifm)%arv(k3v_1,i,j-jdim)
            wtv2 = grid_g(ifm)%arv(k2v,i,j)        / grid_g(ifm)%arv(k3v,i,j)

            theta_mean(i,j)   =       wtw   * basic_g(ifm)%theta(k2w,i,j)                  &
                              + (1. - wtw)  * basic_g(ifm)%theta(k3w,i,j)

            rv_mean(i,j)      =       wtw   * basic_g(ifm)%rv(k2w,i,j)                     &
                              + (1. - wtw)  * basic_g(ifm)%rv(k3w,i,j)

            rtp_mean(i,j)     =       wtw   * basic_g(ifm)%rtp(k2w,i,j)                    &
                              + (1. - wtw)  * basic_g(ifm)%rtp(k3w,i,j)

            if (co2_on) then
               co2p_mean(i,j) =       wtw   * basic_g(ifm)%co2p(k2w,i,j)                   &
                              + (1. - wtw)  * basic_g(ifm)%co2p(k3w,i,j)
            else
               co2p_mean(i,j) = co2con(1)
            end if

            up_mean(i,j)      = (        wtu1  * basic_g(ifm)%up(k2u_1,i-1,j)              &
                                +  (1. - wtu1) * basic_g(ifm)%up(k3u_1,i-1,j)              &
                                +        wtu2  * basic_g(ifm)%up(k2u,i,j)                  &
                                +  (1. - wtu2) * basic_g(ifm)%up(k3u,i,j)       ) * .5
            vp_mean(i,j)      = (        wtv1  * basic_g(ifm)%vp(k2v_1,i,j-jdim)           &
                                +  (1. - wtv1) * basic_g(ifm)%vp(k3v_1,i,j-jdim)           &
                                +        wtv2  * basic_g(ifm)%vp(k2v,i,j)                  &
                                +  (1. - wtv2) * basic_g(ifm)%vp(k3v,i,j)       ) * .5
            
            !------------------------------------------------------------------------------!
            !     Exner function and density, need to consider the layer beneath the       !
            ! ground to make the profile.                                                  !
            !------------------------------------------------------------------------------!
            if (wtw >= .5) then
               pi0_mean(i,j)  = (wtw - .5)                                                 &
                              * (basic_g(ifm)%pp(k1w,i,j) + basic_g(ifm)%pi0(k1w,i,j))     &
                              + (1.5 - wtw)                                                &
                              * (basic_g(ifm)%pp(k2w,i,j) + basic_g(ifm)%pi0(k2w,i,j))
            else
               pi0_mean(i,j)  = (wtw + .5)                                                 &
                              * (basic_g(ifm)%pp(k2w,i,j) + basic_g(ifm)%pi0(k2w,i,j))     &
                              + (.5 - wtw)                                                 &
                              * (basic_g(ifm)%pp(k3w,i,j) + basic_g(ifm)%pi0(k3w,i,j))
            end if

            !------------------------------------------------------------------------------!
            !      The following subroutine calculates diffuse shortwave radiation.  This  !
            ! is a gross approximation that uses constant scattering coefficients.  A more !
            ! appropriate way of doing this was probably already designed by DMM.  Valeriy !
            ! Ivanov is also investigating methods of parameterizing diffuse shortwave     !
            ! radiation.                                                                   !
            !------------------------------------------------------------------------------!
            rshort1 = radiate_g(ifm)%rshort(i,j)
            cosz1   = radiate_g(ifm)%cosz(i,j)
            call short2diff(rshort1,cosz1,rshortd1)
            rshortd(i,j) = rshortd1
        end do
      end do
   end select

   !---------------------------------------------------------------------------------------!
   !     With these averages, we can start copying the values to the ED structures.        !
   !---------------------------------------------------------------------------------------!
   do ipy=1,cgrid%npolygons
      
      !---- Aliases for grid points of interest. ------------------------------------------!
      ix =  cgrid%ilon(ipy)
      iy =  cgrid%ilat(ipy)
      k2w   = nint(grid_g(ifm)%flpw(ix,iy))
      k1w   = k2w - 1

      !------------------------------------------------------------------------------------!
      !     Splitting the shortwave into direct and diffuse.  Also, assign 45% of the      !
      ! total as near infrared, and the remaining as photosynthetically active radiation.  !
      ! These could be probably better retrieved from the radiation model, need to check   !
      ! that.                                                                              !
      !------------------------------------------------------------------------------------!
      cgrid%met(ipy)%nir_beam     =  0.45*(radiate_g(ifm)%rshort(ix,iy)-rshortd(ix,iy))
      cgrid%met(ipy)%nir_diffuse  =  0.45*(rshortd(ix,iy))
      cgrid%met(ipy)%par_beam     =  0.55*(radiate_g(ifm)%rshort(ix,iy)-rshortd(ix,iy))
      cgrid%met(ipy)%par_diffuse  =  0.55*(rshortd(ix,iy))
      !------------------------------------------------------------------------------------!

      cgrid%met(ipy)%rlong    = radiate_g(ifm)%rlong(ix,iy)
      !----- Converting Exner function to pressure. ---------------------------------------!
      cgrid%met(ipy)%prss     = p00 * (cpi * pi0_mean(ix,iy))**cpor

      !----- Finding the actual height above ground for 2nd level. ------------------------!
      cgrid%met(ipy)%geoht    = (zt(k2w)-zm(k1w)) * grid_g(ifm)%rtgt(ix,iy)

      !------------------------------------------------------------------------------------!
      !   Finding the kinetic energy (twice its value, to be precise).  It will be con-    !
      ! verted to wind speed later, after calc_met_lapse.                                  !
      !------------------------------------------------------------------------------------!
      cgrid%met(ipy)%vels     = up_mean(ix,iy)**2 + vp_mean(ix,iy)**2
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!! WHY WERE WE FORCING WIND SPEED TO BE >= 0.65m/s???                       !!!!!!
      ! cgrid%met(ipy)%vels     = max(0.65,sqrt( up_mean(ix,iy)**2 + vp_mean(ix,iy)**2))
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !------------------------------------------------------------------------------------!
      !    ED needs specific humidity, and temperature, but BRAMS has mixing ratio and     !
      ! potential temperature, so we must convert before sending to ED.  CO2 is already in !
      ! µmol_CO2 / mol_air, no conversion needed...                                        !
      !------------------------------------------------------------------------------------!
      cgrid%met(ipy)%atm_shv  = rv_mean(ix,iy) / (1.+rtp_mean(ix,iy))
      cgrid%met(ipy)%atm_tmp  = theta_mean(ix,iy)*pi0_mean(ix,iy) * cpi
      cgrid%met(ipy)%atm_co2  = co2p_mean(ix,iy)
   end do

   !----- Filling the precipitation arrays. -----------------------------------------------!
   call fill_site_precip(ifm,cgrid,m2,m3,ia,iz,ja,jz,pi0_mean,theta_mean)

  
   do ipy = 1,cgrid%npolygons

      !------------------------------------------------------------------------------------!
      !      If this is the first call, it is only looking for atm_tmp and atm_shv, so we  !
      ! can bypass the error check.                                                        !
      !------------------------------------------------------------------------------------!
      if (init) cgrid%met(ipy)%rlong=rlong_min


      !------ Apply met to sites, and adjust met variables for topography. ----------------!
      call calc_met_lapse(cgrid,ipy)

      !----- Converting vels to wind (m/s) ------------------------------------------------!
      cgrid%met(ipy)%vels = sqrt(max(0.0,cgrid%met(ipy)%vels))


      cpoly => cgrid%polygon(ipy)
      do isi = 1,cpoly%nsites
         
         !----- Vels.  At this point vels is 2*Kinetic Energy, take the square root. ------!
         cpoly%met(isi)%vels        = sqrt(max(0.0,cpoly%met(isi)%vels))
         cpoly%met(isi)%vels_stab   = max(0.1,cpoly%met(isi)%vels)
         cpoly%met(isi)%vels_unstab = max(1.0,cpoly%met(isi)%vels_stab)

         !----- Exner function, simply copy from average ----------------------------------!
         cpoly%met(isi)%exner = pi0_mean(ix,iy)

         !----- Solar radiation -----------------------------------------------------------!
         cpoly%met(isi)%rshort_diffuse = cpoly%met(isi)%par_diffuse                        &
                                       + cpoly%met(isi)%nir_diffuse

         cpoly%met(isi)%rshort         = cpoly%met(isi)%rshort_diffuse                     &
                                       + cpoly%met(isi)%par_beam + cpoly%met(isi)%nir_beam
         
         !---------------------------------------------------------------------------------!
         !    Finding the mean energy and depth of precipitation rates: qpcpg and dpcpg,   !
         ! respectively.  This is just an approximation, because we are assuming all ice   !
         ! or all liquid, and a standard snow density.  We could do better by using the    !
         ! actual values of these things from the microphysics. For the cumulus parameter- !
         ! ization, that's the best we can do...                                           !
         !---------------------------------------------------------------------------------!
         if (cpoly%met(isi)%atm_tmp > t3ple) then
            cpoly%met(isi)%qpcpg = cliq * (cpoly%met(isi)%atm_tmp - tsupercool)            &
                                 * cpoly%met(isi)%pcpg
            cpoly%met(isi)%dpcpg = max(0.0, cpoly%met(isi)%pcpg * 0.001)
         else
            cpoly%met(isi)%qpcpg = cice * cpoly%met(isi)%atm_tmp * cpoly%met(isi)%pcpg
            cpoly%met(isi)%dpcpg = max(0.0, cpoly%met(isi)%pcpg * 0.01)
         end if
      end do
   end do

   return
end subroutine copy_atm2lsm
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!      This subroutine will consider both kinds of precipitation (cumulus scheme + bulk    !
! microphysics) to retrieve precipitation for ED. It is currently a simplified method,     !
! because it doesn't fully take advantage of the bulk microphysics details.                !
!------------------------------------------------------------------------------------------!
subroutine fill_site_precip(ifm,cgrid,m2,m3,ia,iz,ja,jz,pi0_mean,theta_mean)
   use mem_cuparm   , only : cuparm_g     & ! structure
                           , nnqparm      ! ! intent(in)
   use mem_micro    , only : micro_g      ! ! structure
   use micphys      , only : availcat     ! ! intent(in)
   use therm_lib    , only : bulk_on      ! ! intent(in)
   use mem_basic    , only : basic_g      ! ! structure
   use rconstants   , only : cpi          & ! intent(in)
                           , cliq         ! ! intent(in)
   use ed_state_vars, only : edtype       ! ! structure
   use mem_edcp     , only : ed_precip_g  ! ! structure
   use ed_misc_coms , only : dtlsm        ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   type(edtype)          , target     :: cgrid
   integer               , intent(in) :: ifm,m2,m3,ia,iz,ja,jz
   real, dimension(m2,m3), intent(in) :: pi0_mean,theta_mean
   !----- Local variables -----------------------------------------------------------------!
   integer                            :: i,j,n,ix,iy,ipy
   real, dimension(m2,m3)             :: conprr_mean,bulkprr_mean
   real                               :: totpcp, dtlsmi
   logical                            :: cumulus_on
   !---------------------------------------------------------------------------------------!


   !----- Some aliases.  ------------------------------------------------------------------!
   cumulus_on   = nnqparm(ifm) > 0
   dtlsmi = 1./dtlsm


   !----- Zero the precipitation structures. ----------------------------------------------!
   do ipy=1,cgrid%npolygons
      cgrid%met(ipy)%pcpg = 0
   end do
   do i=ia,iz
      do j=ja,jz
         conprr_mean(i,j)  = 0.
         bulkprr_mean(i,j) = 0.
      end do
   end do

   !---------------------------------------------------------------------------------------!
   !     Here we need to check which kind of precipitation we have available, since both   !
   ! cumulus and bulk microphysics may or may not exist.                                   !
   !---------------------------------------------------------------------------------------!
   !----- Cumulus clouds from cumulus scheme. ---------------------------------------------!
   if (cumulus_on) then
      !------------------------------------------------------------------------------------!
      !     We need to use difference in aconpr instead of conprr because cumulus parame-  !
      ! trisation is usually called less often than ED.  Therefore, using conprr would     !
      ! repeat the same precipitation rate until cumulus is called again, which would make !
      ! ED overestimate the total convective precipitation.  Using aconpr and comparing to !
      ! the previous one will give the average precipitation rate between two ED calls, so !
      ! it will take the right total amount of convective precipitation.                   !
      !------------------------------------------------------------------------------------!
      do i=ia,iz
         do j=ja,jz
            conprr_mean(i,j) = dtlsmi * ( cuparm_g(ifm)%aconpr(i,j)                        &
                                        - ed_precip_g(ifm)%prev_aconpr(i,j) )
            ed_precip_g(ifm)%prev_aconpr(i,j) = cuparm_g(ifm)%aconpr(i,j)
         end do
      end do
   end if
   !----- Bulk microphysics method. -------------------------------------------------------!
   if (bulk_on) then

      !------------------------------------------------------------------------------------!
      !     Again, we need to use the accumulated precipitation because the bulk micro-    !
      ! physics is called every BRAMS time step.  Therefore, if we use pcpg we will cap-   !
      ! ture the precipitation that happened only at one time step before the current ED   !
      ! call, which would underestimate the resolved precipitation unless                  !
      ! dtlsm = dtlong(ifm).  Using the accumulated precipitation difference will give the !
      ! amount of resolved precipitation since the last ED call, thus taking the right     !
      ! amount of resolved precipitation.                                                  !
      !------------------------------------------------------------------------------------!
      do i=ia,iz
         do j=ja,jz
            !------------------------------------------------------------------------------!
            !     Here I add all forms of precipitation available, previously converted    !
            ! into the equivalent amount of liquid water.  Even when the bulk microphysics !
            ! scheme is on, we may not have all forms of precipitation, so we check one by !
            ! one before adding them.                                                      !
            !------------------------------------------------------------------------------!
            totpcp = 0.
            if (availcat(2)) totpcp = totpcp + micro_g(ifm)%accpr(i,j) ! Rain
            if (availcat(3)) totpcp = totpcp + micro_g(ifm)%accpp(i,j) ! Pristine ice
            if (availcat(4)) totpcp = totpcp + micro_g(ifm)%accps(i,j) ! Snow
            if (availcat(5)) totpcp = totpcp + micro_g(ifm)%accpa(i,j) ! Aggregates
            if (availcat(6)) totpcp = totpcp + micro_g(ifm)%accpg(i,j) ! Graupel
            if (availcat(7)) totpcp = totpcp + micro_g(ifm)%accph(i,j) ! Hail
            !----- Finding the average rate. ----------------------------------------------!
            bulkprr_mean(i,j) = (totpcp-ed_precip_g(ifm)%prev_abulkpr(i,j)) * dtlsmi
            !----- Save current amount for next call. -------------------------------------!
            ed_precip_g(ifm)%prev_abulkpr(i,j) = totpcp
         end do
      end do
   end if

   !----- Now we combine both kinds of precipitation. -------------------------------------!
   do ipy=1,cgrid%npolygons
      ix =  cgrid%ilon(ipy)
      iy =  cgrid%ilat(ipy)
      cgrid%met(ipy)%pcpg = conprr_mean(ix,iy) + bulkprr_mean(ix,iy)
   end do

   return
end subroutine fill_site_precip
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!      Subroutine to copy the old "future" flux field to the "past".                       !
!------------------------------------------------------------------------------------------!
subroutine copy_fluxes_future_2_past(ifm)
   use mem_edcp,only: ed_fluxf_g  & ! structure
                    , ed_fluxp_g  ! ! structure
                           
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer               , intent(in) :: ifm
   !---------------------------------------------------------------------------------------!

   !----- Simple copy of fluxes over land. ------------------------------------------------!
   ed_fluxp_g(ifm)%ustar   = ed_fluxf_g(ifm)%ustar
   ed_fluxp_g(ifm)%tstar   = ed_fluxf_g(ifm)%tstar
   ed_fluxp_g(ifm)%rstar   = ed_fluxf_g(ifm)%rstar
   ed_fluxp_g(ifm)%cstar   = ed_fluxf_g(ifm)%cstar
   ed_fluxp_g(ifm)%albedt  = ed_fluxf_g(ifm)%albedt
   ed_fluxp_g(ifm)%rlongup = ed_fluxf_g(ifm)%rlongup
   ed_fluxp_g(ifm)%sflux_u = ed_fluxf_g(ifm)%sflux_u
   ed_fluxp_g(ifm)%sflux_v = ed_fluxf_g(ifm)%sflux_v
   ed_fluxp_g(ifm)%sflux_w = ed_fluxf_g(ifm)%sflux_w
   ed_fluxp_g(ifm)%sflux_t = ed_fluxf_g(ifm)%sflux_t
   ed_fluxp_g(ifm)%sflux_r = ed_fluxf_g(ifm)%sflux_r
   ed_fluxp_g(ifm)%sflux_c = ed_fluxf_g(ifm)%sflux_c

   return
end subroutine copy_fluxes_future_2_past
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine will copy the ED-2 variables that are used by BRAMS in other sub-    !
! routines, such as radiation, turbulence, etc.                                            !
!------------------------------------------------------------------------------------------! 
subroutine copy_fluxes_lsm2atm(ifm)
   use ed_state_vars , only : edgrid_g    & ! structure
                            , edtype      & ! structure
                            , polygontype & ! structure
                            , sitetype    ! ! structure
   use mem_edcp      , only : ed_fluxf_g  & ! structure
                            , ed_flux     ! ! structure
   use mem_grid      , only : zt          & ! intent(in)
                            , grid_g      & ! structure
                            , dzt         & ! intent(in)
                            , zm          & ! intent(in)
                            , if_adap     & ! intent(in)
                            , jdim        ! ! intent(in)
   use mem_basic     , only : basic_g     ! ! structure
   use mem_leaf      , only : leaf_g      ! ! structure
   use rconstants    , only : mmcod       ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer               , intent(in) :: ifm
   !----- Local variables -----------------------------------------------------------------!
   type(edtype)          , pointer    :: cgrid
   type(polygontype)     , pointer    :: cpoly
   type(sitetype)        , pointer    :: csite
   type(ed_flux)         , pointer    :: fluxp
   integer                            :: ipy,isi
   integer                            :: ix,iy
   integer                            :: k2u,k3u,k2u_1,k3u_1
   integer                            :: k2v,k3v,k2v_1,k3v_1
   real                               :: wtu1,wtu2,wtv1,wtv2
   real                               :: up_mean,vp_mean,cosine,sine
   real                               :: poly_area_i, site_area_i
   real(kind=8)                       :: angle
   !---------------------------------------------------------------------------------------!


   !----- Set the pointers, the state variable and the future flux grid. ------------------!
   cgrid => edgrid_g(ifm)
   fluxp => ed_fluxf_g(ifm)
  
   do ipy=1,cgrid%npolygons
      cpoly => cgrid%polygon(ipy)

      ix = cgrid%ilon(ipy)
      iy = cgrid%ilat(ipy)
      poly_area_i = 1./ sum(cpoly%area)

      !----- Initialising all fluxes. -----------------------------------------------------!
      fluxp%ustar(ix,iy)   = 0.0
      fluxp%tstar(ix,iy)   = 0.0
      fluxp%rstar(ix,iy)   = 0.0
      fluxp%cstar(ix,iy)   = 0.0
      fluxp%sflux_u(ix,iy)   = 0.0
      fluxp%sflux_v(ix,iy)   = 0.0
      fluxp%sflux_w(ix,iy)   = 0.0
      fluxp%sflux_t(ix,iy)   = 0.0
      fluxp%sflux_r(ix,iy)   = 0.0
      fluxp%sflux_c(ix,iy)   = 0.0
      fluxp%rlongup(ix,iy) = 0.0
      fluxp%albedt(ix,iy)  = 0.0

      !------------------------------------------------------------------------------------!
      !     Finding the wind average.  To do this, we first check which vertical coordi-   !
      ! nate we are using.                                                                 !
      !------------------------------------------------------------------------------------!
      select case(if_adap)
      case (0) !----- Terrain-following coordinate. ---------------------------------------!
         up_mean = (basic_g(ifm)%up(2,ix,iy) + basic_g(ifm)%up(2,ix-1,iy)) * 0.5
         vp_mean = (basic_g(ifm)%vp(2,ix,iy) + basic_g(ifm)% vp(2,ix,iy-jdim)) * 0.5
      case (1) !----- Adaptive coordinate. ------------------------------------------------!
         !---------------------------------------------------------------------------------!
         !    Finding the lowest boundary levels, depending on which kind of variable we   !
         ! are averaging (staggered grid).                                                 !
         !---------------------------------------------------------------------------------!
         !----- U points, need to average between i-1 and i... ----------------------------!
         k2u   = nint(grid_g(ifm)%flpu(ix,iy))
         k3u   = k2u + 1
         k2u_1 = nint(grid_g(ifm)%flpu(ix-1,iy))
         k3u_1 = k2u_1 + 1
         !----- V points, need to average between j-jdim and j... -------------------------!
         k2v   = nint(grid_g(ifm)%flpv(ix,iy))
         k3v   = k2v+1
         k2v_1 = nint(grid_g(ifm)%flpv(ix,iy-jdim))
         k3v_1 = k2v_1 + 1
         !---------------------------------------------------------------------------------!
         !     Computing the weights for lowest predicted points, relative to points above !
         ! them.                                                                           !
         !---------------------------------------------------------------------------------!
         wtu1 = grid_g(ifm)%aru(k2u_1,ix-1,iy)    / grid_g(ifm)%aru(k3u_1,ix-1,iy)
         wtu2 = grid_g(ifm)%aru(k2u,ix,iy)        / grid_g(ifm)%aru(k3u,ix,iy)
         wtv1 = grid_g(ifm)%arv(k2v_1,ix,iy-jdim) / grid_g(ifm)%arv(k3v_1,ix,iy-jdim)
         wtv2 = grid_g(ifm)%arv(k2v,ix,iy)        / grid_g(ifm)%arv(k3v,ix,iy)
         up_mean = (        wtu1  * basic_g(ifm)%up(k2u_1,ix-1,iy)                         &
                   +  (1. - wtu1) * basic_g(ifm)%up(k3u_1,ix-1,iy)                         &
                   +        wtu2  * basic_g(ifm)%up(k2u,ix,iy)                             &
                   +  (1. - wtu2) * basic_g(ifm)%up(k3u,ix,iy)       ) * .5
         vp_mean = (        wtv1  * basic_g(ifm)%vp(k2v_1,ix,iy-jdim)                      &
                   +  (1. - wtv1) * basic_g(ifm)%vp(k3v_1,ix,iy-jdim)                      &
                   +        wtv2  * basic_g(ifm)%vp(k2v,ix,iy)                             &
                   +  (1. - wtv2) * basic_g(ifm)%vp(k3v,ix,iy)       ) * .5
      end select
      !----- This is safer than using u/sqrt(u²+v²), when u and v are both zero... --------!
      angle  = datan2(dble(vp_mean),dble(up_mean))
      cosine = sngl(dcos(angle))
      sine   = sngl(dsin(angle))


      !------------------------------------------------------------------------------------!
      !    Finding the averaged flux variables.  For most variables, it will be simply     !
      ! the site and patch weighted average, in which the weights are the site and patch   !
      ! areas.  The exceptions are commented.                                              !
      !------------------------------------------------------------------------------------!
      do isi=1,cpoly%nsites
         csite => cpoly%site(isi)
         site_area_i = 1./sum(csite%area)

         fluxp%ustar(ix,iy) = fluxp%ustar(ix,iy)                                           &
                            + cpoly%area(isi)*sum(csite%area * csite%ustar)                &
                            * site_area_i * poly_area_i

         fluxp%tstar(ix,iy) = fluxp%tstar(ix,iy)                                           &
                            + cpoly%area(isi)*sum(csite%area * csite%tstar)                &
                            * site_area_i * poly_area_i

         fluxp%cstar(ix,iy) = fluxp%cstar(ix,iy)                                           &
                            + cpoly%area(isi)*sum(csite%area * csite%cstar)                &
                            * site_area_i * poly_area_i

         !---------------------------------------------------------------------------------!
         !     BRAMS needs rstar (mixing ratio), but ED computes qstar (specific humidity) !
         ! characteristic friction scales.  We must find the appropriate conversion.       !
         !---------------------------------------------------------------------------------!
         fluxp%rstar(ix,iy) = fluxp%rstar(ix,iy)                                           &
                            + cpoly%area(isi)* site_area_i * poly_area_i                   &
                            * sum(csite%area * csite%qstar / (1.-csite%can_shv))           &
                            / (1. - cpoly%met(isi)%atm_shv)
         !---------------------------------------------------------------------------------!

         !---------------------------------------------------------------------------------!
         !     The flux of momentum for the horizontal wind needs to be split with the     !
         ! appropriate direction.                                                          !
         !---------------------------------------------------------------------------------!
         fluxp%sflux_u(ix,iy) = fluxp%sflux_u(ix,iy)                                       &
                              + cosine*cpoly%area(isi)*sum(csite%area * csite%upwp)        &
                              * site_area_i * poly_area_i
         fluxp%sflux_v(ix,iy) = fluxp%sflux_v(ix,iy)                                       &
                              + sine*cpoly%area(isi)*sum(csite%area * csite%upwp)          &
                              * site_area_i * poly_area_i
         !---------------------------------------------------------------------------------!

         fluxp%sflux_w(ix,iy) = fluxp%sflux_w(ix,iy)                                       &
                              + cpoly%area(isi)*sum(csite%area * csite%wpwp)               &
                              * site_area_i * poly_area_i

         fluxp%sflux_t(ix,iy) = fluxp%sflux_t(ix,iy)                                       &
                              + cpoly%area(isi)*sum(csite%area * csite%tpwp)               &
                              * site_area_i * poly_area_i

         !---------------------------------------------------------------------------------!
         !     Same thing as in the rstar case: we must convert the ED-2 flux (specific    !
         ! volume) to mixing ratio, which is what BRAMS expects.                           !
         !---------------------------------------------------------------------------------!
         fluxp%sflux_r(ix,iy) = fluxp%sflux_r(ix,iy)                                       &
                              + cpoly%area(isi)* site_area_i * poly_area_i                 &
                              * sum(csite%area * csite%qpwp / (1.-csite%can_shv))          &
                              / (1. - cpoly%met(isi)%atm_shv)
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !     Same thing as the cstar case: we must convert the ED-2 flux to BRAMS units. !
         !---------------------------------------------------------------------------------!
         fluxp%sflux_c(ix,iy) = fluxp%sflux_c(ix,iy)                                       &
                              + cpoly%area(isi)*sum(csite%area * csite%cpwp)               &
                              * site_area_i * poly_area_i

         !---------------------------------------------------------------------------------!
         !   Include emission and reflected longwave in rlongup.                           !
         !---------------------------------------------------------------------------------!
         fluxp%rlongup(ix,iy) = fluxp%rlongup(ix,iy)                                       &
                              + ( cpoly%area(isi)*cpoly%rlongup(isi)                       &
                                + cgrid%met(ipy)%rlong *cpoly%area(isi)                    &
                                * cpoly%rlong_albedo(isi) ) * site_area_i * poly_area_i
         !---------------------------------------------------------------------------------!

         !---------------------------------------------------------------------------------!
         !    Total albedo is the average between albedo for direct (beam) and diffuse.    !
         !---------------------------------------------------------------------------------!
         fluxp%albedt(ix,iy)  = fluxp%albedt(ix,iy)                                        &
                              + cpoly%area(isi) * site_area_i * poly_area_i                &
                              * 0.5 * (cpoly%albedo_beam(isi) + cpoly%albedo_diffuse(isi) )
         !---------------------------------------------------------------------------------!
      end do
   end do

   return
end subroutine copy_fluxes_lsm2atm
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine will assign initial values for some of the state and flux variables   !
! in LEAF-3, which will get the data from ED-2.  Also, it will perform initialisation of   !
! some structures that are used in the ED-2/LEAF-3 interfacing.                            !
!------------------------------------------------------------------------------------------!
subroutine initialize_ed2leaf(ifm,mxp,myp)
  
   use mem_edcp     , only : ed_fluxf_g     & ! structure
                           , ed_fluxp_g     & ! structure
                           , wgrid_g        & ! structure
                           , ed_precip_g    & ! structure
                           , alloc_edprecip & ! structure
                           , zero_edprecip  & ! structure
                           , alloc_edflux   & ! structure
                           , zero_edflux    & ! structure
                           , alloc_wgrid    & ! structure
                           , zero_wgrid     ! ! structure
   use node_mod     , only : mmxp           & ! intent(in)
                           , mmyp           & ! intent(in)
                           , ia             & ! intent(in)
                           , iz             & ! intent(in)
                           , ja             & ! intent(in)
                           , jz             & ! intent(in)
                           , mynum          ! ! intent(in)
   use ed_work_vars , only : work_e         ! ! structure
   use mem_leaf     , only : leaf_g         ! ! structure
   use mem_basic    , only : basic_g        ! ! structure
   use mem_grid     , only : grid_g         & ! structure
                           , zt             & ! intent(in)
                           , dzt            & ! intent(in)
                           , zm             & ! intent(in)
                           , if_adap        & ! intent(in)
                           , jdim           ! ! intent(in)
   use rconstants   , only : cpi            ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer               , intent(in) :: ifm, mxp, myp
   !----- Local variables -----------------------------------------------------------------!
   integer                             :: ix,iy,i,j
   integer                             :: k2w,k3w,k1w
   real                                :: topma_t,wtw
   real,dimension(mmxp(ifm),mmyp(ifm)) :: theta_mean,pi0_mean,rv_mean
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   ! Step 1 - Set the patch are fraction in the leaf array to the land fraction in the     !
   !          work array. The leaf array has the border cells and is 2 cells larger in     !
   !          each dimension than the work array, thus the ix and iy instead of a simple   !
   !          loop.                                                                        !
   !---------------------------------------------------------------------------------------!
   leaf_g(ifm)%patch_area(:,:,1)=1.0
   leaf_g(ifm)%patch_area(:,:,2)=0.0
   do i=1,mxp
      do j=1,myp
         ix = work_e(ifm)%xatm(i,j)
         iy = work_e(ifm)%yatm(i,j)
         leaf_g(ifm)%patch_area(ix,iy,1) = 1.0-work_e(ifm)%landfrac(i,j)
         leaf_g(ifm)%patch_area(ix,iy,2) = work_e(ifm)%landfrac(i,j)
      end do
   end do

   !---------------------------------------------------------------------------------------!
   ! Step 2 - Allocate the flux arrays from terrestrial and water bodies.  These arrays    !
   !          are the same size as the leaf arrays, and 2 greater than the work arrays.    !
   !---------------------------------------------------------------------------------------!
   call alloc_edflux(ed_fluxf_g(ifm),mmxp(ifm),mmyp(ifm))
   call alloc_edflux(ed_fluxp_g(ifm),mmxp(ifm),mmyp(ifm))
   call alloc_wgrid(wgrid_g(ifm),mmxp(ifm),mmyp(ifm))
   call alloc_edprecip(ed_precip_g(ifm),mmxp(ifm),mmyp(ifm))
   !----- Assign zeroes to the newly allocated matrices. ----------------------------------!
   call zero_edflux(ed_fluxf_g(ifm))
   call zero_edflux(ed_fluxp_g(ifm))
   call zero_wgrid(wgrid_g(ifm))
   call zero_edprecip (ed_precip_g(ifm))
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !    Computing the bottom=most layer, then transfer the value to state variables.  The  !
   ! bottom-most layer definition depends on the vertical coordinate.                      !
   !---------------------------------------------------------------------------------------!
   select case (if_adap)
   case (0) !------ Terrain-following coordinates. ----------------------------------------!
      do j=ja,jz
         do i=ia,iz
            theta_mean(i,j) = basic_g(ifm)%theta(2,i,j)
            rv_mean(i,j)    = basic_g(ifm)%rv(2,i,j)
            pi0_mean(i,j)   = ( basic_g(ifm)%pp(1,i,j)  + basic_g(ifm)%pp(2,i,j)           &
                              + basic_g(ifm)%pi0(1,i,j) + basic_g(ifm)%pi0(2,i,j) ) * 0.5
         end do
      end do
   case (1)
      do j=ja,jz
         do i=ia,iz
            !----- Weighted average -------------------------------------------------------!
            k2w = nint(grid_g(ifm)%flpw(i,j)) 
            k1w = k2w - 1
            k3w = k2w + 1
            topma_t = .25 * ( grid_g(ifm)%topma(i,j) + grid_g(ifm)%topma(i-1,j)            &
                            + grid_g(ifm)%topma(i,j-jdim) + grid_g(ifm)%topma(i-1,j-jdim))
            wtw     = (zm(k2w) - topma_t) * dzt(k2w)

            theta_mean(i,j) =       wtw  * basic_g(ifm)%theta(k2w,i,j)                     &
                            + (1. - wtw) * basic_g(ifm)%theta(k3w,i,j)
            rv_mean(i,j)    =       wtw  * basic_g(ifm)%rv(k2w,i,j)                        &
                            + (1. - wtw) * basic_g(ifm)%rv(k3w,i,j)

            if (wtw >= .5) then
               pi0_mean(i,j)  = (wtw - .5)                                                 &
                              * (basic_g(ifm)%pp(k1w,i,j) + basic_g(ifm)%pi0(k1w,i,j))     &
                              + (1.5 - wtw)                                                &
                              * (basic_g(ifm)%pp(k2w,i,j) + basic_g(ifm)%pi0(k2w,i,j))
            else
               pi0_mean(i,j)  = (wtw + .5)                                                 &
                              * (basic_g(ifm)%pp(k2w,i,j) + basic_g(ifm)%pi0(k2w,i,j))     &
                              + (.5 - wtw)                                                 &
                              * (basic_g(ifm)%pp(k3w,i,j) + basic_g(ifm)%pi0(k3w,i,j))
            end if
         end do
      end do
   end select

   do j=ja,jz
      do i=ia,iz
         leaf_g(ifm)%can_temp(i,j,1) =  theta_mean(i,j) * pi0_mean(i,j) * cpi
         leaf_g(ifm)%can_rvap(i,j,1) =  rv_mean(i,j)
         leaf_g(ifm)%can_temp(i,j,2) =  theta_mean(i,j) * pi0_mean(i,j) * cpi
         leaf_g(ifm)%can_rvap(i,j,2) =  rv_mean(i,j)
     
         leaf_g(ifm)%gpp(i,j)        = 0.0
         leaf_g(ifm)%resphet(i,j)    = 0.0
         leaf_g(ifm)%plresp(i,j)     = 0.0
      end do
   end do
   return
end subroutine initialize_ed2leaf
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine will transfer the fluxes computed in ED-2 to the LEAF-3 arrays, so   !
! the other parts of the code will have the appropriate values.                            !
!------------------------------------------------------------------------------------------!
subroutine transfer_ed2leaf(ifm,timel)
   use mem_edcp      , only : ed_fluxf_g  & ! structure
                            , ed_fluxp_g  & ! structure
                            , wgrid_g     & ! structure
                            , edtime1     & ! intent(in)
                            , edtime2     ! ! intent(in)
   use mem_turb      , only : turb_g      ! ! structure
   use rconstants    , only : stefan      & ! intent(in)
                            , g           ! ! intent(in)
   use mem_leaf      , only : leaf_g      & ! structure
                            , zrough      ! ! intent(in)
   use mem_radiate   , only : radiate_g   ! ! structure
   use node_mod      , only : master_num  & ! intent(in)
                            , mmzp        & ! intent(in)
                            , mmxp        & ! intent(in)
                            , mmyp        & ! intent(in)
                            , ia          & ! intent(in)
                            , iz          & ! intent(in)
                            , ja          & ! intent(in)
                            , jz          & ! intent(in)
                            , ia_1        & ! intent(in)
                            , iz1         & ! intent(in)
                            , ja_1        & ! intent(in)
                            , jz1         & ! intent(in)
                            , ibcon       ! ! intent(in)
   use mem_grid      , only : jdim        & ! intent(in)
                            , grid_g      ! ! structure
   use ed_state_vars , only : edgrid_g    & ! structure
                            , edtype      & ! structure
                            , polygontype & ! structure
                            , sitetype    ! ! structure
   use soil_coms     , only : soil_rough  ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer     , intent(in):: ifm
   real(kind=8), intent(in) :: timel
   !----- Local variables -----------------------------------------------------------------!
   type(edtype)     , pointer   :: cgrid
   type(polygontype), pointer   :: cpoly
   type(sitetype)   , pointer   :: csite
   real                         :: tfact,la
   integer                      :: ic,jc,ici,jci,i,j,ix,iy
   integer                      :: m2,m3
   integer                      :: ipy, isi
   real                         :: site_area_i, polygon_area_i
   !----- Local constants -----------------------------------------------------------------!
   real             , parameter :: z0fac_water  = 0.016/g
   real             , parameter :: snowrough    = 0.001
   real             , parameter :: z0_min_water = 0.0001
   !---------------------------------------------------------------------------------------!

   !----- Assigning some aliases ----------------------------------------------------------!
   m2 = mmxp(ifm)
   m3 = mmyp(ifm)
   tfact = (timel-edtime1)/(edtime2-edtime1)

  
   !---------------------------------------------------------------------------------------!
   !     There are effectively two LEAF-3 patches: patch 1 is water, patch 2 is land (the  !
   ! ED polygon level.  We could also assign leaf_g(ifm)%zrough based on the surface       !
   ! characteristics of ED, and discretize the sites or ed-patches into leaf patches, but  !
   ! for now they are left homogeneous as land.                                            !
   !---------------------------------------------------------------------------------------!
  
  
   do j=ja,jz
      do i=ia,iz 
         !---------------------------------------------------------------------------------!
         !      For the stars, we use a simple linear interpolation on time.               !
         !---------------------------------------------------------------------------------!
         leaf_g(ifm)%ustar(i,j,2) = (1. - tfact) * ed_fluxp_g(ifm)%ustar(i,j)              &
                                  +       tfact  * ed_fluxf_g(ifm)%ustar(i,j)
         leaf_g(ifm)%tstar(i,j,2) = (1. - tfact) * ed_fluxp_g(ifm)%tstar(i,j)              &
                                  +       tfact  * ed_fluxf_g(ifm)%tstar(i,j)
         leaf_g(ifm)%rstar(i,j,2) = (1. - tfact) * ed_fluxp_g(ifm)%rstar(i,j)              &
                                  +       tfact  * ed_fluxf_g(ifm)%rstar(i,j)
         leaf_g(ifm)%cstar(i,j,2) = (1. - tfact) * ed_fluxp_g(ifm)%cstar(i,j)              &
                                  +       tfact  * ed_fluxf_g(ifm)%cstar(i,j)

         !---------------------------------------------------------------------------------!
         !      For the albedo and surface fluxes, we must blend land and water.  The land !
         ! (ED polygon) is done similarly to the stars, while the water uses the instant-  !
         ! aneous value.                                                                   !
         !---------------------------------------------------------------------------------!
         radiate_g(ifm)%albedt(i,j) = leaf_g(ifm)%patch_area(i,j,2)                        &
                                    * ( (1.-tfact) * ed_fluxp_g(ifm)%albedt(i,j)           &
                                      +     tfact  * ed_fluxf_g(ifm)%albedt(i,j) )         &
                                    + leaf_g(ifm)%patch_area(i,j,1)                        &
                                    * wgrid_g(ifm)%albedt(i,j)

         radiate_g(ifm)%rlongup(i,j) = leaf_g(ifm)%patch_area(i,j,2)                       &
                                     * ( (1.-tfact) * ed_fluxp_g(ifm)%rlongup(i,j)         &
                                       +     tfact  * ed_fluxf_g(ifm)%rlongup(i,j) )       &
                                     + leaf_g(ifm)%patch_area(i,j,1)                       &
                                     * wgrid_g(ifm)%rlongup(i,j)

         turb_g(ifm)%sflux_u(i,j)    = leaf_g(ifm)%patch_area(i,j,2)                       &
                                     * ( (1.-tfact) * ed_fluxp_g(ifm)%sflux_u(i,j)         &
                                       +     tfact  * ed_fluxf_g(ifm)%sflux_u(i,j))        &
                                     + leaf_g(ifm)%patch_area(i,j,1)                       &
                                     * wgrid_g(ifm)%sflux_u(i,j)

         turb_g(ifm)%sflux_v(i,j)    = leaf_g(ifm)%patch_area(i,j,2)                       &
                                     * ( (1.-tfact) * ed_fluxp_g(ifm)%sflux_v(i,j)         &
                                       +     tfact  * ed_fluxf_g(ifm)%sflux_v(i,j))        &
                                     + leaf_g(ifm)%patch_area(i,j,1)                       &
                                     * wgrid_g(ifm)%sflux_v(i,j)

         turb_g(ifm)%sflux_w(i,j)    = leaf_g(ifm)%patch_area(i,j,2)                       &
                                     * ( (1.-tfact) * ed_fluxp_g(ifm)%sflux_w(i,j)         &
                                       +     tfact  * ed_fluxf_g(ifm)%sflux_w(i,j) )       &
                                     + leaf_g(ifm)%patch_area(i,j,1)                       &
                                     * wgrid_g(ifm)%sflux_w(i,j)

         turb_g(ifm)%sflux_t(i,j)    = leaf_g(ifm)%patch_area(i,j,2)                       &
                                     * ( (1.-tfact) * ed_fluxp_g(ifm)%sflux_t(i,j)         &
                                       +     tfact  * ed_fluxf_g(ifm)%sflux_t(i,j) )       &
                                     + leaf_g(ifm)%patch_area(i,j,1)                       &
                                     * wgrid_g(ifm)%sflux_t(i,j)

         turb_g(ifm)%sflux_r(i,j)    = leaf_g(ifm)%patch_area(i,j,2)                       &
                                     * ( (1.-tfact) * ed_fluxp_g(ifm)%sflux_r(i,j)         &
                                       +     tfact  * ed_fluxf_g(ifm)%sflux_r(i,j))        &
                                     + leaf_g(ifm)%patch_area(i,j,1)                       &
                                     * wgrid_g(ifm)%sflux_r(i,j)

         turb_g(ifm)%sflux_c(i,j)    = leaf_g(ifm)%patch_area(i,j,2)                       &
                                     * ( (1.-tfact) * ed_fluxp_g(ifm)%sflux_c(i,j)         &
                                       +     tfact  * ed_fluxf_g(ifm)%sflux_c(i,j))        &
                                     + leaf_g(ifm)%patch_area(i,j,1)                       &
                                     * wgrid_g(ifm)%sflux_c(i,j)

         !----- Copy the water-body fluxes, they are synchronised with BRAMS --------------!
         leaf_g(ifm)%ustar(i,j,1) = wgrid_g(ifm)%ustar(i,j)
         leaf_g(ifm)%tstar(i,j,1) = wgrid_g(ifm)%tstar(i,j)
         leaf_g(ifm)%rstar(i,j,1) = wgrid_g(ifm)%rstar(i,j)
         leaf_g(ifm)%cstar(i,j,1) = wgrid_g(ifm)%cstar(i,j)

         !----- Find roughness scales for water bodies ------------------------------------!
         leaf_g(ifm)%patch_rough(i,j,1) = max(z0fac_water * leaf_g(ifm)%ustar(i,j,1) ** 2  &
                                             ,z0_min_water)
         leaf_g(ifm)%soil_rough(i,j,1)  = 0.0
         leaf_g(ifm)%veg_rough(i,j,1)   = 0.0
         
         !----- Initializing the land bodies. They will remain 0 over the ocean. ----------!
         leaf_g(ifm)%veg_rough(i,j,2)   = 0.0
         leaf_g(ifm)%soil_rough(i,j,2)  = soil_rough
         leaf_g(ifm)%patch_rough(i,j,2) = 0.0
      end do
   end do

   !---------------------------------------------------------------------------------------!
   !     Computing the vegetation and patch roughness for land.  We assign the polygon     !
   ! average in both cases.                                                                !
   !---------------------------------------------------------------------------------------!
   cgrid => edgrid_g(ifm)
   polyloop: do ipy=1,cgrid%npolygons
      cpoly => cgrid%polygon(ipy)
      ix     = cgrid%ilon(ipy)
      iy    = cgrid%ilat(ipy)     
      polygon_area_i = 1. / sum(cpoly%area(:))

      siteloop: do isi = 1, cpoly%nsites
         csite => cpoly%site(isi)
         site_area_i = 1. / sum(csite%area(:))

         leaf_g(ifm)%veg_rough(ix,iy,2)   = leaf_g(ifm)%veg_rough(ix,iy,2)                 &
                                          + sum(csite%veg_rough(:)*csite%area(:))          &
                                          * cpoly%area(isi)  * site_area_i * polygon_area_i
         leaf_g(ifm)%patch_rough(ix,iy,2) = leaf_g(ifm)%patch_rough(ix,iy,2)               &
                                          + sum(csite%rough(:)*csite%area(:))              &
                                          * cpoly%area(isi) * site_area_i * polygon_area_i
      end do siteloop
   end do polyloop
  
   do j=ja,jz
      do i=ia,iz
         !----- Not sure about this one, but this is consistent with LEAF-3 ---------------!
         leaf_g(ifm)%patch_rough(i,j,2) = max(leaf_g(ifm)%patch_rough(i,j,2)               &
                                             ,leaf_g(ifm)%soil_rough(i,j,2)                &
                                             ,grid_g(ifm)%topzo(i,j))
      end do
   end do
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !    The boundary cells of these arrays have not been filled.  These must be filled by  !
   ! the adjacen cells, likely the 2nd or 2nd to last cell in each row or column.  This    !
   ! will be done only for the true boundary cells (i.e. the ones at the edge of the total !
   ! domain).                                                                              !
   !---------------------------------------------------------------------------------------!



   !----- Western Boundary ----------------------------------------------------------------!
   if (iand(ibcon,1) /= 0) then
      do j=ja,jz
         radiate_g(ifm)%albedt(1,j) = radiate_g(ifm)%albedt(2,j)
         radiate_g(ifm)%rlongup(1,j)= radiate_g(ifm)%rlongup(2,j)

         turb_g(ifm)%sflux_u(1,j)   = turb_g(ifm)%sflux_u(2,j)
         turb_g(ifm)%sflux_v(1,j)   = turb_g(ifm)%sflux_v(2,j)
         turb_g(ifm)%sflux_w(1,j)   = turb_g(ifm)%sflux_w(2,j)
         turb_g(ifm)%sflux_t(1,j)   = turb_g(ifm)%sflux_t(2,j)
         turb_g(ifm)%sflux_r(1,j)   = turb_g(ifm)%sflux_r(2,j)
         turb_g(ifm)%sflux_c(1,j)   = turb_g(ifm)%sflux_c(2,j)
 
         leaf_g(ifm)%ustar(1,j,1)   = leaf_g(ifm)%ustar(2,j,1)
         leaf_g(ifm)%rstar(1,j,1)   = leaf_g(ifm)%rstar(2,j,1)
         leaf_g(ifm)%tstar(1,j,1)   = leaf_g(ifm)%tstar(2,j,1)
         leaf_g(ifm)%cstar(1,j,1)   = leaf_g(ifm)%cstar(2,j,1)
 
         leaf_g(ifm)%ustar(1,j,2)   = leaf_g(ifm)%ustar(2,j,2)
         leaf_g(ifm)%rstar(1,j,2)   = leaf_g(ifm)%rstar(2,j,2)
         leaf_g(ifm)%tstar(1,j,2)   = leaf_g(ifm)%tstar(2,j,2)
         leaf_g(ifm)%cstar(1,j,2)   = leaf_g(ifm)%cstar(2,j,2)
      end do
   end if

   !----- Eastern Boundary ----------------------------------------------------------------!
   if (iand(ibcon,2) /= 0) then
      do j = ja, jz
         radiate_g(ifm)%albedt(m2,j) = radiate_g(ifm)%albedt(m2-1,j)
         radiate_g(ifm)%rlongup(m2,j)= radiate_g(ifm)%rlongup(m2-1,j)

         turb_g(ifm)%sflux_u(m2,j)= turb_g(ifm)%sflux_u(m2-1,j)
         turb_g(ifm)%sflux_v(m2,j)= turb_g(ifm)%sflux_v(m2-1,j)
         turb_g(ifm)%sflux_w(m2,j)= turb_g(ifm)%sflux_w(m2-1,j)
         turb_g(ifm)%sflux_t(m2,j)= turb_g(ifm)%sflux_t(m2-1,j)
         turb_g(ifm)%sflux_r(m2,j)= turb_g(ifm)%sflux_r(m2-1,j)
         turb_g(ifm)%sflux_c(m2,j)= turb_g(ifm)%sflux_c(m2-1,j)
      
         leaf_g(ifm)%ustar(m2,j,1)=leaf_g(ifm)%ustar(m2-1,j,1)
         leaf_g(ifm)%rstar(m2,j,1)=leaf_g(ifm)%rstar(m2-1,j,1)
         leaf_g(ifm)%tstar(m2,j,1)=leaf_g(ifm)%tstar(m2-1,j,1)
         leaf_g(ifm)%cstar(m2,j,1)=leaf_g(ifm)%cstar(m2-1,j,1)

         leaf_g(ifm)%ustar(m2,j,2)=leaf_g(ifm)%ustar(m2-1,j,2)
         leaf_g(ifm)%rstar(m2,j,2)=leaf_g(ifm)%rstar(m2-1,j,2)
         leaf_g(ifm)%tstar(m2,j,2)=leaf_g(ifm)%tstar(m2-1,j,2)
         leaf_g(ifm)%cstar(m2,j,2)=leaf_g(ifm)%cstar(m2-1,j,2)
      end do
   end if
  
   !----- Southern Boundary ---------------------------------------------------------------!
   if (jdim == 1 .and. iand(ibcon,4) /= 0) then
      do i = ia,iz
         radiate_g(ifm)%albedt(i,1) = radiate_g(ifm)%albedt(i,2)
         radiate_g(ifm)%rlongup(i,1)= radiate_g(ifm)%rlongup(i,2)

         turb_g(ifm)%sflux_u(i,1)= turb_g(ifm)%sflux_u(i,2)
         turb_g(ifm)%sflux_v(i,1)= turb_g(ifm)%sflux_v(i,2)
         turb_g(ifm)%sflux_w(i,1)= turb_g(ifm)%sflux_w(i,2)
         turb_g(ifm)%sflux_t(i,1)= turb_g(ifm)%sflux_t(i,2)
         turb_g(ifm)%sflux_r(i,1)= turb_g(ifm)%sflux_r(i,2)
         turb_g(ifm)%sflux_c(i,1)= turb_g(ifm)%sflux_c(i,2)
         
         leaf_g(ifm)%ustar(i,1,1)=leaf_g(ifm)%ustar(i,2,1)
         leaf_g(ifm)%rstar(i,1,1)=leaf_g(ifm)%rstar(i,2,1)
         leaf_g(ifm)%tstar(i,1,1)=leaf_g(ifm)%tstar(i,2,1)
         leaf_g(ifm)%cstar(i,1,1)=leaf_g(ifm)%cstar(i,2,1)

         leaf_g(ifm)%ustar(i,1,2)=leaf_g(ifm)%ustar(i,2,2)
         leaf_g(ifm)%rstar(i,1,2)=leaf_g(ifm)%rstar(i,2,2)
         leaf_g(ifm)%tstar(i,1,2)=leaf_g(ifm)%tstar(i,2,2)
         leaf_g(ifm)%cstar(i,1,2)=leaf_g(ifm)%cstar(i,2,2)
      end do
   end if


   !----- Northern Boundary ---------------------------------------------------------------!
   if (jdim == 1 .and. iand(ibcon,8) /= 0) then
      do i = ia,iz
         radiate_g(ifm)%albedt(i,m3) = radiate_g(ifm)%albedt(i,m3-1)
         radiate_g(ifm)%rlongup(i,m3)= radiate_g(ifm)%rlongup(i,m3-1)

         turb_g(ifm)%sflux_u(i,m3)= turb_g(ifm)%sflux_u(i,m3-1)
         turb_g(ifm)%sflux_v(i,m3)= turb_g(ifm)%sflux_v(i,m3-1)
         turb_g(ifm)%sflux_w(i,m3)= turb_g(ifm)%sflux_w(i,m3-1)
         turb_g(ifm)%sflux_t(i,m3)= turb_g(ifm)%sflux_t(i,m3-1)
         turb_g(ifm)%sflux_r(i,m3)= turb_g(ifm)%sflux_r(i,m3-1)
         turb_g(ifm)%sflux_c(i,m3)= turb_g(ifm)%sflux_c(i,m3-1)
         
         leaf_g(ifm)%ustar(i,m3,1)=leaf_g(ifm)%ustar(i,m3-1,1)
         leaf_g(ifm)%rstar(i,m3,1)=leaf_g(ifm)%rstar(i,m3-1,1)
         leaf_g(ifm)%tstar(i,m3,1)=leaf_g(ifm)%tstar(i,m3-1,1)
         leaf_g(ifm)%cstar(i,m3,1)=leaf_g(ifm)%cstar(i,m3-1,1)

         leaf_g(ifm)%ustar(i,m3,2)=leaf_g(ifm)%ustar(i,m3-1,2)
         leaf_g(ifm)%rstar(i,m3,2)=leaf_g(ifm)%rstar(i,m3-1,2)
         leaf_g(ifm)%tstar(i,m3,2)=leaf_g(ifm)%tstar(i,m3-1,2)
         leaf_g(ifm)%cstar(i,m3,2)=leaf_g(ifm)%cstar(i,m3-1,2)

      end do
   end if

   !----- Southwestern corner -------------------------------------------------------------!
   if (iand(ibcon,5) /= 0 .or. (iand(ibcon,1) /= 0 .and. jdim == 0)) then
      ic=1
      jc=1
      ici=2
      jci=1+jdim

      radiate_g(ifm)%albedt(ic,jc) = radiate_g(ifm)%albedt(ici,jci)
      radiate_g(ifm)%rlongup(ic,jc)= radiate_g(ifm)%rlongup(ici,jci)

      turb_g(ifm)%sflux_u(ic,jc)= turb_g(ifm)%sflux_u(ici,jci)
      turb_g(ifm)%sflux_v(ic,jc)= turb_g(ifm)%sflux_v(ici,jci)
      turb_g(ifm)%sflux_w(ic,jc)= turb_g(ifm)%sflux_w(ici,jci)
      turb_g(ifm)%sflux_t(ic,jc)= turb_g(ifm)%sflux_t(ici,jci)
      turb_g(ifm)%sflux_r(ic,jc)= turb_g(ifm)%sflux_r(ici,jci)
      turb_g(ifm)%sflux_c(ic,jc)= turb_g(ifm)%sflux_c(ici,jci)

      leaf_g(ifm)%ustar(ic,jc,1)=leaf_g(ifm)%ustar(ici,jci,1)
      leaf_g(ifm)%rstar(ic,jc,1)=leaf_g(ifm)%rstar(ici,jci,1)
      leaf_g(ifm)%tstar(ic,jc,1)=leaf_g(ifm)%tstar(ici,jci,1)
      leaf_g(ifm)%cstar(ic,jc,1)=leaf_g(ifm)%cstar(ici,jci,1)

      leaf_g(ifm)%ustar(ic,jc,2)=leaf_g(ifm)%ustar(ici,jci,2)
      leaf_g(ifm)%rstar(ic,jc,2)=leaf_g(ifm)%rstar(ici,jci,2)
      leaf_g(ifm)%tstar(ic,jc,2)=leaf_g(ifm)%tstar(ici,jci,2)
      leaf_g(ifm)%cstar(ic,jc,2)=leaf_g(ifm)%cstar(ici,jci,2)
   end if

   !----- Southeastern corner -------------------------------------------------------------!
   if (iand(ibcon,6) /= 0 .or. (iand(ibcon,1) /= 0 .and. jdim == 0)) then
      ic=m2
      jc=1
      ici=m2-1
      jci=1+jdim

      radiate_g(ifm)%albedt(ic,jc) = radiate_g(ifm)%albedt(ici,jci)
      radiate_g(ifm)%rlongup(ic,jc)= radiate_g(ifm)%rlongup(ici,jci)

      turb_g(ifm)%sflux_u(ic,jc)= turb_g(ifm)%sflux_u(ici,jci)
      turb_g(ifm)%sflux_v(ic,jc)= turb_g(ifm)%sflux_v(ici,jci)
      turb_g(ifm)%sflux_w(ic,jc)= turb_g(ifm)%sflux_w(ici,jci)
      turb_g(ifm)%sflux_t(ic,jc)= turb_g(ifm)%sflux_t(ici,jci)
      turb_g(ifm)%sflux_r(ic,jc)= turb_g(ifm)%sflux_r(ici,jci)
      turb_g(ifm)%sflux_c(ic,jc)= turb_g(ifm)%sflux_c(ici,jci)

      leaf_g(ifm)%ustar(ic,jc,1)=leaf_g(ifm)%ustar(ici,jci,1)
      leaf_g(ifm)%rstar(ic,jc,1)=leaf_g(ifm)%rstar(ici,jci,1)
      leaf_g(ifm)%tstar(ic,jc,1)=leaf_g(ifm)%tstar(ici,jci,1)
      leaf_g(ifm)%cstar(ic,jc,1)=leaf_g(ifm)%cstar(ici,jci,1)

      leaf_g(ifm)%ustar(ic,jc,2)=leaf_g(ifm)%ustar(ici,jci,2)
      leaf_g(ifm)%rstar(ic,jc,2)=leaf_g(ifm)%rstar(ici,jci,2)
      leaf_g(ifm)%tstar(ic,jc,2)=leaf_g(ifm)%tstar(ici,jci,2)
      leaf_g(ifm)%cstar(ic,jc,2)=leaf_g(ifm)%cstar(ici,jci,2)
   end if
  
   !----- Northeastern corner -------------------------------------------------------------!
   if (iand(ibcon,9) /= 0 .or. (iand(ibcon,2) /= 0 .and. jdim == 0)) then
      ic=m2
      jc=m3
      ici=m2-1
      jci=m3-jdim

      radiate_g(ifm)%albedt(ic,jc) = radiate_g(ifm)%albedt(ici,jci)
      radiate_g(ifm)%rlongup(ic,jc)= radiate_g(ifm)%rlongup(ici,jci)

      turb_g(ifm)%sflux_u(ic,jc)= turb_g(ifm)%sflux_u(ici,jci)
      turb_g(ifm)%sflux_v(ic,jc)= turb_g(ifm)%sflux_v(ici,jci)
      turb_g(ifm)%sflux_w(ic,jc)= turb_g(ifm)%sflux_w(ici,jci)
      turb_g(ifm)%sflux_t(ic,jc)= turb_g(ifm)%sflux_t(ici,jci)
      turb_g(ifm)%sflux_r(ic,jc)= turb_g(ifm)%sflux_r(ici,jci)
      turb_g(ifm)%sflux_c(ic,jc)= turb_g(ifm)%sflux_c(ici,jci)

      leaf_g(ifm)%ustar(ic,jc,1)=leaf_g(ifm)%ustar(ici,jci,1)
      leaf_g(ifm)%rstar(ic,jc,1)=leaf_g(ifm)%rstar(ici,jci,1)
      leaf_g(ifm)%tstar(ic,jc,1)=leaf_g(ifm)%tstar(ici,jci,1)
      leaf_g(ifm)%cstar(ic,jc,1)=leaf_g(ifm)%cstar(ici,jci,1)

      leaf_g(ifm)%ustar(ic,jc,2)=leaf_g(ifm)%ustar(ici,jci,2)
      leaf_g(ifm)%rstar(ic,jc,2)=leaf_g(ifm)%rstar(ici,jci,2)
      leaf_g(ifm)%tstar(ic,jc,2)=leaf_g(ifm)%tstar(ici,jci,2)
      leaf_g(ifm)%cstar(ic,jc,2)=leaf_g(ifm)%cstar(ici,jci,2)
   end if
  
   !----- Northwestern corner -------------------------------------------------------------!
   if (iand(ibcon,10) /= 0 .or. (iand(ibcon,2) /= 0 .and. jdim == 0)) then
      ic=1
      jc=m3
      ici=2
      jci=m3-jdim

      radiate_g(ifm)%albedt(ic,jc) = radiate_g(ifm)%albedt(ici,jci)
      radiate_g(ifm)%rlongup(ic,jc)= radiate_g(ifm)%rlongup(ici,jci)

      turb_g(ifm)%sflux_u(ic,jc)= turb_g(ifm)%sflux_u(ici,jci)
      turb_g(ifm)%sflux_v(ic,jc)= turb_g(ifm)%sflux_v(ici,jci)
      turb_g(ifm)%sflux_w(ic,jc)= turb_g(ifm)%sflux_w(ici,jci)
      turb_g(ifm)%sflux_t(ic,jc)= turb_g(ifm)%sflux_t(ici,jci)
      turb_g(ifm)%sflux_r(ic,jc)= turb_g(ifm)%sflux_r(ici,jci)
      turb_g(ifm)%sflux_c(ic,jc)= turb_g(ifm)%sflux_c(ici,jci)

      leaf_g(ifm)%ustar(ic,jc,1)=leaf_g(ifm)%ustar(ici,jci,1)
      leaf_g(ifm)%rstar(ic,jc,1)=leaf_g(ifm)%rstar(ici,jci,1)
      leaf_g(ifm)%tstar(ic,jc,1)=leaf_g(ifm)%tstar(ici,jci,1)
      leaf_g(ifm)%cstar(ic,jc,1)=leaf_g(ifm)%cstar(ici,jci,1)

      leaf_g(ifm)%ustar(ic,jc,2)=leaf_g(ifm)%ustar(ici,jci,2)
      leaf_g(ifm)%rstar(ic,jc,2)=leaf_g(ifm)%rstar(ici,jci,2)
      leaf_g(ifm)%tstar(ic,jc,2)=leaf_g(ifm)%tstar(ici,jci,2)
      leaf_g(ifm)%cstar(ic,jc,2)=leaf_g(ifm)%cstar(ici,jci,2)
   end if


   return
end subroutine transfer_ed2leaf
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine calc_met_lapse(cgrid,ipy)
  
   use ed_state_vars         , only : edtype      & ! structure
                                    , polygontype ! ! structure
   use canopy_radiation_coms , only : rlong_min   ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   type(edtype)     , target     :: cgrid
   integer          , intent(in) :: ipy
   !----- Local variables -----------------------------------------------------------------!
   type(polygontype), pointer    :: cpoly
   integer                       :: isi
   real                          :: ebar   !! mean elevation
   real                          :: delE   !! deviation from mean elevation
   real                          :: aterr  !! terrestrial area
   !----- Local constants -----------------------------------------------------------------!
   real             , parameter  :: offset=tiny(1.)/epsilon(1.) !! Tiny offset to avoid FPE
   logical          , parameter  :: bypass=.true.
   !---------------------------------------------------------------------------------------!

   !----- Pass over sites once to calc preliminary stats. ---------------------------------!
   cpoly => cgrid%polygon(ipy)

   ebar = 0.0
   aterr = 0.0

   do isi=1,cpoly%nsites
      ebar  = ebar + cpoly%area(isi)*cpoly%elevation(isi)
      aterr = aterr + cpoly%area(isi)
   end do
   ebar = ebar/aterr

   if (bypass) then
      do isi = 1,cpoly%nsites
         cpoly%met(isi)%geoht       = cgrid%met(ipy)%geoht
         cpoly%met(isi)%atm_tmp     = cgrid%met(ipy)%atm_tmp
         cpoly%met(isi)%atm_shv     = cgrid%met(ipy)%atm_shv
         cpoly%met(isi)%prss        = cgrid%met(ipy)%prss
         cpoly%met(isi)%pcpg        = cgrid%met(ipy)%pcpg
         cpoly%met(isi)%par_diffuse = cgrid%met(ipy)%par_diffuse
         cpoly%met(isi)%atm_co2     = cgrid%met(ipy)%atm_co2
         cpoly%met(isi)%rlong       = cgrid%met(ipy)%rlong
         cpoly%met(isi)%par_beam    = cgrid%met(ipy)%par_beam
         cpoly%met(isi)%nir_diffuse = cgrid%met(ipy)%nir_diffuse
         cpoly%met(isi)%nir_beam    = cgrid%met(ipy)%nir_beam
         cpoly%met(isi)%vels        = cgrid%met(ipy)%vels
       
         !---------------------------------------------------------------------------------!
         !     Sanity check:                                                               !
         !---------------------------------------------------------------------------------!
         if ( cpoly%met(isi)%rlong < rlong_min) then
            write(unit=*,fmt='(a)')           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            write(unit=*,fmt='(a)')           '+ RLONG is too low!!!'
            write(unit=*,fmt='(a,1x,es12.5)') '+ Site-level      : ',cpoly%met(isi)%rlong
            write(unit=*,fmt='(a,1x,es12.5)') '+ Polygon-level   : ',cgrid%met(ipy)%rlong
            write(unit=*,fmt='(a,1x,es12.5)') '+ Minimum OK value: ',rlong_min
            write(unit=*,fmt='(a)')           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            call fatal_error('Problems with RLONG A','calc_met_lapse'                      &
                            ,'edcp_met.f90')
         elseif ( cpoly%met(isi)%atm_tmp < 150.0) then
            write(unit=*,fmt='(a)')           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            write(unit=*,fmt='(a)')           '+ ATM_TMP is too low!!!'
            write(unit=*,fmt='(a,1x,es12.5)') '+ Site-level      : ',cpoly%met(isi)%atm_tmp
            write(unit=*,fmt='(a,1x,es12.5)') '+ Polygon-level   : ',cgrid%met(ipy)%atm_tmp
            write(unit=*,fmt='(a,1x,es12.5)') '+ Minimum OK value: ',150.0
            write(unit=*,fmt='(a)')           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            call fatal_error('Problems with ATM_TMP','calc_met_lapse'                      &
                            ,'edcp_met.f90')
         else if ( cpoly%met(isi)%atm_shv < 1.e-5) then
            write(unit=*,fmt='(a)')           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            write(unit=*,fmt='(a)')           '+ ATM_SHV is too low!!!'
            write(unit=*,fmt='(a,1x,es12.5)') '+ Site-level      : ',cpoly%met(isi)%atm_shv
            write(unit=*,fmt='(a,1x,es12.5)') '+ Polygon-level   : ',cgrid%met(ipy)%atm_shv
            write(unit=*,fmt='(a,1x,es12.5)') '+ Minimum OK value: ',1.e-5
            write(unit=*,fmt='(a)')           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            call fatal_error('Problems with ATM_SHV','calc_met_lapse'                      &
                            ,'edcp_met.f90')
         elseif ( cpoly%met(isi)%rlong > 600.0) then
            write(unit=*,fmt='(a)')           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            write(unit=*,fmt='(a)')           '+ RLONG is too high!!!'
            write(unit=*,fmt='(a,1x,es12.5)') '+ Site-level      : ',cpoly%met(isi)%rlong
            write(unit=*,fmt='(a,1x,es12.5)') '+ Polygon-level   : ',cgrid%met(ipy)%rlong
            write(unit=*,fmt='(a,1x,es12.5)') '+ Maximum OK value: ',600.0
            write(unit=*,fmt='(a)')           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            call fatal_error('Problems with RLONG A','calc_met_lapse'                      &
                            ,'edcp_met.f90')
         elseif ( cpoly%met(isi)%atm_tmp > 317.0) then
            write(unit=*,fmt='(a)')           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            write(unit=*,fmt='(a)')           '+ ATM_TMP is too high!!!'
            write(unit=*,fmt='(a,1x,es12.5)') '+ Site-level      : ',cpoly%met(isi)%atm_tmp
            write(unit=*,fmt='(a,1x,es12.5)') '+ Polygon-level   : ',cgrid%met(ipy)%atm_tmp
            write(unit=*,fmt='(a,1x,es12.5)') '+ Maximum OK value: ',317.0
            write(unit=*,fmt='(a)')           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            call fatal_error('Problems with ATM_TMP','calc_met_lapse'                      &
                            ,'edcp_met.f90')
         elseif ( cpoly%met(isi)%atm_shv > 30.0e-3) then
            write(unit=*,fmt='(a)')           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            write(unit=*,fmt='(a)')           '+ ATM_SHV is too high!!!'
            write(unit=*,fmt='(a,1x,es12.5)') '+ Site-level      : ',cpoly%met(isi)%atm_shv
            write(unit=*,fmt='(a,1x,es12.5)') '+ Polygon-level   : ',cgrid%met(ipy)%atm_shv
            write(unit=*,fmt='(a,1x,es12.5)') '+ Maximum OK value: ',30.0e-3
            write(unit=*,fmt='(a)')           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            call fatal_error('Problems with ATM_SHV','calc_met_lapse'                      &
                            ,'edcp_met.f90')
         elseif ( cpoly%met(isi)%par_beam + cpoly%met(isi)%nir_beam &
                + cpoly%met(isi)%par_diffuse + cpoly%met(isi)%nir_diffuse > 1320.0 ) then
            write(unit=*,fmt='(a)')           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            write(unit=*,fmt='(a)')           '+ SOLAR RADIATION is non-sense !!!'
            write(unit=*,fmt='(a,1x,es12.5)') '+ PAR_BEAM        : ',cpoly%met(isi)%par_beam
            write(unit=*,fmt='(a,1x,es12.5)') '+ PAR_DIFFUSE     : '                       &
                                              , cpoly%met(isi)%par_diffuse
            write(unit=*,fmt='(a,1x,es12.5)') '+ NIR_BEAM        : ',cpoly%met(isi)%nir_beam
            write(unit=*,fmt='(a,1x,es12.5)') '+ NIR_DIFFUSE     : '                       &
                                              , cpoly%met(isi)%nir_diffuse
            write(unit=*,fmt='(a,1x,es12.5)') '+ Maximum OK value (sum): ',1320.0
            write(unit=*,fmt='(a)')           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            call fatal_error('Problems with solar radiation','calc_met_lapse'              &
                            ,'ed_met_driver.f90')
         end if
      end do
    
   else
      
      !----- Second pass, calculate lapse rate adjustment. --------------------------------!
      do isi = 1,cpoly%nsites
         
         delE = cpoly%elevation(isi) - ebar
         
         !----- Perform linear adjustments. -----------------------------------------------!
         cpoly%met(isi)%geoht   = cgrid%met(ipy)%geoht   + cgrid%lapse(ipy)%geoht   * delE
         cpoly%met(isi)%atm_tmp = cgrid%met(ipy)%atm_tmp + cgrid%lapse(ipy)%atm_tmp * delE
         cpoly%met(isi)%atm_shv = cgrid%met(ipy)%atm_shv + cgrid%lapse(ipy)%atm_shv * delE
         cpoly%met(isi)%prss    = cgrid%met(ipy)%prss    + cgrid%lapse(ipy)%prss    * delE
         cpoly%met(isi)%pcpg    = cgrid%met(ipy)%pcpg    + cgrid%lapse(ipy)%pcpg    * delE
         cpoly%met(isi)%atm_co2 = cgrid%met(ipy)%atm_co2 + cgrid%lapse(ipy)%atm_co2 * delE
         cpoly%met(isi)%rlong   = cgrid%met(ipy)%rlong   + cgrid%lapse(ipy)%rlong   * delE
         cpoly%met(isi)%par_diffuse = cgrid%met(ipy)%par_diffuse                           &
                                    + cgrid%lapse(ipy)%par_diffuse * delE
         cpoly%met(isi)%par_beam    = cgrid%met(ipy)%par_beam                              &
                                    + cgrid%lapse(ipy)%par_beam * delE
         cpoly%met(isi)%nir_diffuse = cgrid%met(ipy)%nir_diffuse                           &
                                    + cgrid%lapse(ipy)%nir_diffuse * delE
         cpoly%met(isi)%nir_beam    = cgrid%met(ipy)%nir_beam                              &
                                    + cgrid%lapse(ipy)%nir_beam * delE
         !---------------------------------------------------------------------------------!
         ! Note: at this point VELS is vel^2.  Thus this lapse preserves mean wind ENERGY  !
         !       not wind SPEED.                                                           !
         !---------------------------------------------------------------------------------!
         cpoly%met(isi)%vels    = cgrid%met(ipy)%vels    + cgrid%lapse(ipy)%vels*delE

       
         !---------------------------------------------------------------------------------!
         !     Sanity check:                                                               !
         !---------------------------------------------------------------------------------!
         if ( cpoly%met(isi)%rlong < rlong_min) then
            write(unit=*,fmt='(a)')           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            write(unit=*,fmt='(a)')           '+ RLONG is too low!!!'
            write(unit=*,fmt='(a,1x,es12.5)') '+ Site-level      : ',cpoly%met(isi)%rlong
            write(unit=*,fmt='(a,1x,es12.5)') '+ Polygon-level   : ',cgrid%met(ipy)%rlong
            write(unit=*,fmt='(a,1x,es12.5)') '+ Minimum OK value: ',rlong_min
            write(unit=*,fmt='(a)')           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            call fatal_error('Problems with RLONG A','calc_met_lapse'                      &
                            ,'edcp_met.f90')
         elseif ( cpoly%met(isi)%atm_tmp < 150.0) then
            write(unit=*,fmt='(a)')           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            write(unit=*,fmt='(a)')           '+ ATM_TMP is too low!!!'
            write(unit=*,fmt='(a,1x,es12.5)') '+ Site-level      : ',cpoly%met(isi)%atm_tmp
            write(unit=*,fmt='(a,1x,es12.5)') '+ Polygon-level   : ',cgrid%met(ipy)%atm_tmp
            write(unit=*,fmt='(a,1x,es12.5)') '+ Minimum OK value: ',150.0
            write(unit=*,fmt='(a)')           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            call fatal_error('Problems with ATM_TMP','calc_met_lapse'                      &
                            ,'edcp_met.f90')
         else if ( cpoly%met(isi)%atm_shv < 1.e-5) then
            write(unit=*,fmt='(a)')           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            write(unit=*,fmt='(a)')           '+ ATM_SHV is too low!!!'
            write(unit=*,fmt='(a,1x,es12.5)') '+ Site-level      : ',cpoly%met(isi)%atm_shv
            write(unit=*,fmt='(a,1x,es12.5)') '+ Polygon-level   : ',cgrid%met(ipy)%atm_shv
            write(unit=*,fmt='(a,1x,es12.5)') '+ Minimum OK value: ',1.e-5
            write(unit=*,fmt='(a)')           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            call fatal_error('Problems with ATM_SHV','calc_met_lapse'                      &
                            ,'edcp_met.f90')
         elseif ( cpoly%met(isi)%rlong > 600.0) then
            write(unit=*,fmt='(a)')           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            write(unit=*,fmt='(a)')           '+ RLONG is too high!!!'
            write(unit=*,fmt='(a,1x,es12.5)') '+ Site-level      : ',cpoly%met(isi)%rlong
            write(unit=*,fmt='(a,1x,es12.5)') '+ Polygon-level   : ',cgrid%met(ipy)%rlong
            write(unit=*,fmt='(a,1x,es12.5)') '+ Maximum OK value: ',600.0
            write(unit=*,fmt='(a)')           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            call fatal_error('Problems with RLONG A','calc_met_lapse'                      &
                            ,'edcp_met.f90')
         elseif ( cpoly%met(isi)%atm_tmp > 317.0) then
            write(unit=*,fmt='(a)')           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            write(unit=*,fmt='(a)')           '+ ATM_TMP is too high!!!'
            write(unit=*,fmt='(a,1x,es12.5)') '+ Site-level      : ',cpoly%met(isi)%atm_tmp
            write(unit=*,fmt='(a,1x,es12.5)') '+ Polygon-level   : ',cgrid%met(ipy)%atm_tmp
            write(unit=*,fmt='(a,1x,es12.5)') '+ Maximum OK value: ',317.0
            write(unit=*,fmt='(a)')           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            call fatal_error('Problems with ATM_TMP','calc_met_lapse'                      &
                            ,'edcp_met.f90')
         elseif ( cpoly%met(isi)%atm_shv > 30.0e-3) then
            write(unit=*,fmt='(a)')           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            write(unit=*,fmt='(a)')           '+ ATM_SHV is too high!!!'
            write(unit=*,fmt='(a,1x,es12.5)') '+ Site-level      : ',cpoly%met(isi)%atm_shv
            write(unit=*,fmt='(a,1x,es12.5)') '+ Polygon-level   : ',cgrid%met(ipy)%atm_shv
            write(unit=*,fmt='(a,1x,es12.5)') '+ Maximum OK value: ',30.0e-3
            write(unit=*,fmt='(a)')           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            call fatal_error('Problems with ATM_SHV','calc_met_lapse'                      &
                            ,'edcp_met.f90')
         elseif ( cpoly%met(isi)%par_beam + cpoly%met(isi)%nir_beam &
                + cpoly%met(isi)%par_diffuse + cpoly%met(isi)%nir_diffuse > 1320.0 ) then
            write(unit=*,fmt='(a)')           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            write(unit=*,fmt='(a)')           '+ SOLAR RADIATION is non-sense !!!'
            write(unit=*,fmt='(a,1x,es12.5)') '+ PAR_BEAM        : ',cpoly%met(isi)%par_beam
            write(unit=*,fmt='(a,1x,es12.5)') '+ PAR_DIFFUSE     : '                       &
                                              , cpoly%met(isi)%par_diffuse
            write(unit=*,fmt='(a,1x,es12.5)') '+ NIR_BEAM        : ',cpoly%met(isi)%nir_beam
            write(unit=*,fmt='(a,1x,es12.5)') '+ NIR_DIFFUSE     : '                       &
                                              , cpoly%met(isi)%nir_diffuse
            write(unit=*,fmt='(a,1x,es12.5)') '+ Maximum OK value (sum): ',1320.0
            write(unit=*,fmt='(a)')           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            call fatal_error('Problems with solar radiation','calc_met_lapse'              &
                            ,'ed_met_driver.f90')
         end if
      end do
   end if
   return
end subroutine calc_met_lapse
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutines increments the time averaged polygon met-forcing variables.  These  !
! will be normalized by the output period to give time averages of each quanity.  The      !
! polygon level variables are derived from the weighted spatial average from the site      !
! level quantities.                                                                        !
!------------------------------------------------------------------------------------------!
subroutine int_met_avg(cgrid)
   use ed_state_vars , only : edtype      & ! structure
                            , polygontype & ! structure
                            , sitetype    & ! structure
                            , patchtype   ! ! structure
   use ed_misc_coms  , only : dtlsm       & ! intent(in)
                            , frqfast     ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   type(edtype)      , target  :: cgrid
   !----- Local variables -----------------------------------------------------------------!
   type(polygontype) , pointer :: cpoly
   type(sitetype)    , pointer :: csite
   type(patchtype)   , pointer :: cpatch
   integer                     :: ipy,isi,ipa,ico
   real                        :: frqfasti,tfact
   !---------------------------------------------------------------------------------------!

   !----- Some aliases. -------------------------------------------------------------------!
   frqfasti = 1.0 / frqfast
   tfact = dtlsm * frqfasti

   do ipy = 1,cgrid%npolygons
      cpoly => cgrid%polygon(ipy)

      do isi = 1,cpoly%nsites
         
         cgrid%avg_nir_beam(ipy)    = cgrid%avg_nir_beam(ipy)                              &
                                    + cpoly%met(isi)%nir_beam * cpoly%area(isi) * tfact
         cgrid%avg_nir_diffuse(ipy) = cgrid%avg_nir_diffuse(ipy)                           &
                                    + cpoly%met(isi)%nir_diffuse * cpoly%area(isi) * tfact
         cgrid%avg_par_beam(ipy)    = cgrid%avg_par_beam(ipy)                              &
                                    + cpoly%met(isi)%par_beam * cpoly%area(isi) * tfact
         cgrid%avg_par_diffuse(ipy) = cgrid%avg_par_diffuse(ipy)                           &
                                    + cpoly%met(isi)%par_diffuse * cpoly%area(isi) * tfact
         cgrid%avg_atm_tmp(ipy)     = cgrid%avg_atm_tmp(ipy)                               &
                                    + cpoly%met(isi)%atm_tmp * cpoly%area(isi) * tfact
         cgrid%avg_atm_shv(ipy)     = cgrid%avg_atm_shv(ipy)                               &
                                    + cpoly%met(isi)%atm_shv * cpoly%area(isi) * tfact
         cgrid%avg_rshort(ipy)      = cgrid%avg_rshort(ipy)                                &
                                    + cpoly%met(isi)%rshort * cpoly%area(isi) * tfact
         cgrid%avg_rshort_diffuse(ipy) = cgrid%avg_rshort_diffuse(ipy)                     &
                                       + cpoly%met(isi)%rshort_diffuse * cpoly%area(isi)   &
                                       * tfact
         cgrid%avg_rlong(ipy)       = cgrid%avg_rlong(ipy)                                 &
                                    + cpoly%met(isi)%rlong * cpoly%area(isi) * tfact
         cgrid%avg_pcpg(ipy)        = cgrid%avg_pcpg(ipy)                                  &
                                    + cpoly%met(isi)%pcpg * cpoly%area(isi) * tfact
         cgrid%avg_qpcpg(ipy)       = cgrid%avg_qpcpg(ipy)                                 &
                                    + cpoly%met(isi)%qpcpg * cpoly%area(isi) * tfact
         cgrid%avg_dpcpg(ipy)       = cgrid%avg_dpcpg(ipy)                                 &
                                    + cpoly%met(isi)%dpcpg * cpoly%area(isi) * tfact
         cgrid%avg_vels(ipy)        = cgrid%avg_vels(ipy)                                  &
                                    + cpoly%met(isi)%vels * cpoly%area(isi) * tfact
         cgrid%avg_prss(ipy)        = cgrid%avg_prss(ipy)                                  &
                                    + cpoly%met(isi)%prss * cpoly%area(isi) * tfact
         cgrid%avg_exner(ipy)       = cgrid%avg_exner(ipy)                                 &
                                    + cpoly%met(isi)%exner * cpoly%area(isi) * tfact
         cgrid%avg_geoht(ipy)       = cgrid%avg_geoht(ipy)                                 &
                                    + cpoly%met(isi)%geoht * cpoly%area(isi) * tfact
         cgrid%avg_atm_co2(ipy)     = cgrid%avg_atm_co2(ipy)                               &
                                    + cpoly%met(isi)%atm_co2 * cpoly%area(isi) * tfact
         cgrid%avg_albedt(ipy)      = cgrid%avg_albedt(ipy)                                &
                                    + 0.5 * ( cpoly%albedo_beam(isi)                       &
                                            + cpoly%albedo_diffuse(isi) )                  &
                                    * cpoly%area(isi) * tfact
         cgrid%avg_rlongup(ipy)     = cgrid%avg_rlongup(ipy)                               &
                                    + cpoly%rlongup(isi) * cpoly%area(isi) * tfact
      end do
   end do
   return
end subroutine int_met_avg
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine copy_avgvars_to_leaf(ifm)

   use ed_state_vars , only: edgrid_g,edtype,polygontype,sitetype,patchtype
   use mem_leaf      , only: leaf_g
   use mem_grid      , only: nzg,nzs
   use rconstants    , only: t3ple,cliqvlme,cicevlme,allivlme
   use soil_coms     , only: soil
   use ed_misc_coms  , only: frqsum
   implicit none
   
   !----- Argument ------------------------------------------------------------------------!
   integer, intent(in)  :: ifm
   !----- Local variables -----------------------------------------------------------------!
   type(edtype)     , pointer :: cgrid
   type(polygontype), pointer :: cpoly
   type(sitetype)   , pointer :: csite
   type(patchtype)  , pointer :: cpatch
   integer                    :: ipy,isi,ipa,ico
   integer                    :: ix,iy,k
   real                       :: frqsumi,site_area_i,poly_area_i
   !---------------------------------------------------------------------------------------!

   !----- Set the pointers ----------------------------------------------------------------!
   cgrid => edgrid_g(ifm)

   frqsumi = 1.0 / frqsum
   do ipy=1,cgrid%npolygons
      cpoly => cgrid%polygon(ipy)

      ix = cgrid%ilon(ipy)
      iy = cgrid%ilat(ipy)
      
      do k=1,nzg
         leaf_g(ifm)%soil_text(k,ix,iy,2)   = cgrid%ntext_soil(k,ipy)
         leaf_g(ifm)%soil_energy(k,ix,iy,2) = cgrid%avg_soil_energy(k,ipy)
         leaf_g(ifm)%soil_water(k,ix,iy,2)  = cgrid%avg_soil_water(k,ipy)
      end do
      !----- Surface water is always 1, because we give the averaged value. ---------------!
      leaf_g(ifm)%sfcwater_nlev     (ix,iy,2) = 1.
      leaf_g(ifm)%sfcwater_energy (1,ix,iy,2) = cgrid%avg_snowenergy(ipy)
      leaf_g(ifm)%sfcwater_mass   (1,ix,iy,2) = cgrid%avg_snowmass  (ipy)
      leaf_g(ifm)%sfcwater_depth  (1,ix,iy,2) = cgrid%avg_snowdepth (ipy)
      do k=2,nzs
         leaf_g(ifm)%sfcwater_energy (k,ix,iy,2) = 0.
         leaf_g(ifm)%sfcwater_mass   (k,ix,iy,2) = 0.
         leaf_g(ifm)%sfcwater_depth  (k,ix,iy,2) = 0.
      end do
      
      leaf_g(ifm)%veg_water(ix,iy,2) = cgrid%avg_veg_water(ipy)
      leaf_g(ifm)%veg_temp(ix,iy,2)  = cgrid%avg_veg_temp(ipy)
      
      leaf_g(ifm)%can_temp(ix,iy,2)  = cgrid%avg_can_temp(ipy)
      leaf_g(ifm)%can_co2(ix,iy,2)  = cgrid%avg_can_co2(ipy)
      !----- ED uses specific humidity, converting it to mixing ratio. --------------------!
      leaf_g(ifm)%can_rvap(ix,iy,2)  = cgrid%avg_can_shv(ipy)/(1.-cgrid%avg_can_shv(ipy))
      
      
      leaf_g(ifm)%veg_lai(ix,iy,2)   = cgrid%lai(ipy)
      leaf_g(ifm)%veg_tai(ix,iy,2)   = cgrid%lai(ipy) + cgrid%wai(ipy)


      
      leaf_g(ifm)%gpp(ix,iy)         = 0.0
      leaf_g(ifm)%resphet(ix,iy)     = 0.0
      leaf_g(ifm)%plresp(ix,iy)      = 0.0

      poly_area_i = 1./sum(cpoly%area)
      do isi=1,cpoly%nsites
         csite => cpoly%site(isi)
         if (csite%npatches>0) then
            site_area_i=1./sum(csite%area)
            
            leaf_g(ifm)%gpp(ix,iy) = leaf_g(ifm)%gpp(ix,iy) + &
                 sum(csite%area*csite%co2budget_gpp)*cpoly%area(isi)   *frqsumi
            leaf_g(ifm)%plresp(ix,iy) = leaf_g(ifm)%plresp(ix,iy) + &
                 sum(csite%area*csite%co2budget_plresp)*cpoly%area(isi)*frqsumi
           
            leaf_g(ifm)%resphet(ix,iy) = leaf_g(ifm)%resphet(ix,iy) + &
                 sum(csite%area*csite%co2budget_rh) *cpoly%area(isi)   *frqsumi
         end if
      end do

   end do
   return

   
end subroutine copy_avgvars_to_leaf
!==========================================================================================!
!==========================================================================================!
