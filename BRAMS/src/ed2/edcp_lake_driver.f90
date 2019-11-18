!==========================================================================================!
!==========================================================================================!
!     This subroutine will compute the fluxes between water and the air.  ED solves only   !
! the land regions, so this subroutine is a simplified version of LEAF-3 without consider- !
! ing land, only water.  It's called lake model but it is used for oceans and rivers too   !
! (using prescribed surface temperature.)                                                  !
!------------------------------------------------------------------------------------------!
subroutine simple_lake_model()
   use node_mod               , only : ja                & ! intent(in)
                                     , jz                & ! intent(in)
                                     , ia                & ! intent(in)
                                     , iz                & ! intent(in)
                                     , mynum             ! ! intent(in)
   use consts_coms            , only : stefan            & ! intent(in)
                                     , vonk              & ! intent(in)
                                     , grav              & ! intent(in)
                                     , rdry              & ! intent(in)
                                     , t00               & ! intent(in)
                                     , epim1             & ! intent(in)
                                     , mmdryi            & ! intent(in)
                                     , mmdry             ! ! intent(in)
   use io_params              , only : ssttime1          & ! intent(in)
                                     , ssttime2          & ! intent(in)
                                     , iupdsst           ! ! intent(in)
   use mem_edcp               , only : ed_fluxf_g        & ! structure
                                     , ed_fluxp_g        & ! structure
                                     , edtime2           ! ! intent(in)
   use mem_leaf               , only : leaf_g            & ! structure
                                     , dtleaf            ! ! intent(in)
   use mem_radiate            , only : radiate_g         ! ! structure
   use mem_grid               , only : ngrid             & ! intent(in)
                                     , nzg               & ! intent(in)
                                     , nzs               & ! intent(in)
                                     , time              & ! intent(in)
                                     , dtlt              ! ! intent(in)
   use leaf_coms              , only : min_waterrough    & ! intent(in)
                                     , can_depth_min     ! ! intent(in)
   use rk4_coms               , only : rk4min_can_shv    & ! intent(in)
                                     , rk4max_can_shv    & ! intent(in)
                                     , rk4min_can_temp   & ! intent(in)
                                     , rk4max_can_temp   & ! intent(in)
                                     , rk4min_can_co2    & ! intent(in)
                                     , rk4max_can_co2    ! ! intent(in)
   use therm_lib              , only : rhovsil           & ! function
                                     , reducedpress      & ! function
                                     , idealdenssh       ! ! function
   use lake_coms              , only : lake_buff         & ! intent(out)
                                     , initial_lake_buff & ! subroutine
                                     , zero_lakesite     ! ! subroutine
   implicit none
   !----- Local variables -----------------------------------------------------------------!
   integer                  :: i
   integer                  :: j
   integer                  :: k
   real                     :: timefac_sst
   real(kind=8)             :: dsst_dt ! Derivative of sea surface temperature
   !----- Locally saved variables. --------------------------------------------------------!
   logical     , save       :: first_time = .true.
   !---------------------------------------------------------------------------------------!


   !----- If this is the first time, nullify the structure. -------------------------------!
   if (first_time) then
      call initial_lake_buff()
      first_time = .false.
   end if
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     If we have time-varying water surface temperature, we perform a linear interpol-  !
   ! ation over time now. Otherwise, we just use the initially prescribed value (which     !
   ! can be spatially heterogeneous.                                                       !
   !---------------------------------------------------------------------------------------!
   if (iupdsst == 0) then
      timefac_sst = 0. !----- This will force it to use the initial value only. -----------!
   else
      timefac_sst = sngl((edtime2 - ssttime1(ngrid)) / (ssttime2(ngrid) - ssttime1(ngrid)))
   end if
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    Big domain loops start here.                                                       !
   !---------------------------------------------------------------------------------------!
   jloop: do j=ja,jz
      iloop: do i=ia,iz


         !---------------------------------------------------------------------------------!
         !     Find the sea surface temperature, which will be used for the integration.   !
         !---------------------------------------------------------------------------------!
         leaf_g(ngrid)%ground_temp(i,j,1) = leaf_g(ngrid)%seatp(i,j)                       &
                                          + timefac_sst * ( leaf_g(ngrid)%seatf(i,j)       &
                                                          - leaf_g(ngrid)%seatp(i,j) )
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !     Find the time derivative of the sea surface temperature, which will be used !
         ! for the integration.                                                            !
         !---------------------------------------------------------------------------------!
         dsst_dt = (dble(leaf_g(ngrid)%seatf(i,j)) - dble(leaf_g(ngrid)%seatp(i,j)))       &
                 / (ssttime2(ngrid) - ssttime1(ngrid))
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !     Flush all buffers to zero, for a clean start.                               !
         !---------------------------------------------------------------------------------!
         call zero_lakesite(lake_buff%initp)
         call zero_lakesite(lake_buff%yscal)
         call zero_lakesite(lake_buff%y    )
         call zero_lakesite(lake_buff%dydx )
         call zero_lakesite(lake_buff%yerr )
         call zero_lakesite(lake_buff%ytemp)
         call zero_lakesite(lake_buff%ak2  )
         call zero_lakesite(lake_buff%ak3  )
         call zero_lakesite(lake_buff%ak4  )
         call zero_lakesite(lake_buff%ak5  )
         call zero_lakesite(lake_buff%ak6  )
         call zero_lakesite(lake_buff%ak7  )
         !---------------------------------------------------------------------------------!



         !----- Set up the meteorological forcing structure. ------------------------------!
         call copy_met_2_lake(i,j,ngrid,dsst_dt)
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !     Assign the initial values for "canopy" properties.                          !
         !---------------------------------------------------------------------------------!
         call copy_lake_init(i,j,ngrid,lake_buff%initp)
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Integrate the state variables and fluxes.                                   !
         !---------------------------------------------------------------------------------!
         call integrate_lake(dtlt,ed_fluxp_g(ngrid)%rk4step(i,j,1))
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !     Copy the state variables and the fluxes to the output arrays.               !
         !---------------------------------------------------------------------------------!
         call copy_lake_brams(i,j,ngrid,nzg,nzs,lake_buff%initp)
         !---------------------------------------------------------------------------------!
      end do iloop
   end do jloop

   return
end subroutine simple_lake_model
!==========================================================================================!
!==========================================================================================!





!==========================================================================================!
!==========================================================================================!
!     This subroutine will initialise the surface fields for the simple lake model.        !
!------------------------------------------------------------------------------------------!
subroutine copy_met_2_lake(i,j,ifm,dsst_dt)
   use canopy_air_coms        , only : ubmin             ! ! intent(in)
   use mem_basic              , only : co2_on            & ! intent(in)
                                     , co2con            & ! intent(in)
                                     , basic_g           ! ! structure
   use mem_leaf               , only : leaf_g            ! ! intent(in)
   use mem_radiate            , only : radiate_g         ! ! structure
   use mem_grid               , only : zt                & ! intent(in)
                                     , grid_g            & ! structure
                                     , dzt               & ! intent(in)
                                     , zm                & ! intent(in)
                                     , if_adap           & ! intent(in)
                                     , jdim              & ! intent(in)
                                     , ngrid             ! ! intent(in)
   use therm_lib8             , only : idealdenssh8      & ! function
                                     , rehuil8           & ! function
                                     , reducedpress8     & ! function
                                     , press2exner8      & ! function
                                     , exner2press8      & ! function
                                     , extheta2temp8     & ! function
                                     , tq2enthalpy8      ! ! function
   use lake_coms              , only : lakemet           ! ! intent(out)
   use canopy_air_coms        , only : ubmin8            ! ! intent(in)
   use leaf_coms              , only : can_depth_min     ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer     , intent(in)     :: i
   integer     , intent(in)     :: j
   integer     , intent(in)     :: ifm
   real(kind=8), intent(in)     :: dsst_dt
   !----- Local variables. ----------------------------------------------------------------!
   integer                      :: k1w
   integer                      :: k2w
   integer                      :: k3w
   integer                      :: k2u
   integer                      :: k3u
   integer                      :: k2u_1
   integer                      :: k3u_1
   integer                      :: k2v
   integer                      :: k3v
   integer                      :: k2v_1
   integer                      :: k3v_1
   logical                      :: ok_flpoint
   real(kind=4)                 :: topma_t
   real(kind=4)                 :: wtw
   real(kind=4)                 :: wtu1
   real(kind=4)                 :: wtu2
   real(kind=4)                 :: wtv1
   real(kind=4)                 :: wtv2
   real(kind=4)                 :: exner_mean
   real(kind=4)                 :: theta_mean
   real(kind=4)                 :: co2p_mean
   real(kind=4)                 :: up_mean
   real(kind=4)                 :: vp_mean
   real(kind=4)                 :: rv_mean
   real(kind=4)                 :: rtp_mean
   real(kind=4)                 :: zref_mean
   real(kind=8)                 :: can_theta8
   real(kind=8)                 :: can_shv8
   real(kind=8)                 :: can_depth8
   real(kind=8)                 :: can_exner8
   real(kind=8)                 :: can_prss8
   real(kind=8)                 :: angle
   !----- External functions. -------------------------------------------------------------!
   logical           , external :: is_finite
   logical           , external :: is_finite8
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !   We now retrieve some variables from the atmosphere.  We may need to do some inter-  !
   ! polation in the vertical, for example, when we use adaptive coordinate.  Therefore,   !
   ! we always check this.                                                                 !
   !---------------------------------------------------------------------------------------!
   select case (if_adap)
   case (0)
      !------------------------------------------------------------------------------------!
      !      Terrain-following coordinate (sigma_z).  Here we simply use the lowest pre-   !
      ! dicted level for thermo and wind variables, and a simple average between levels    !
      ! one and two for Exner function and density.  I am not sure this is consistent with !
      ! what is done for density in adaptive coordinate...                                 !
      !------------------------------------------------------------------------------------!
      !----- U points, need to average between i-1 and i... -------------------------------!
      k2u   = nint(grid_g(ifm)%flpu(i,j))
      k3u   = k2u + 1
      k2u_1 = nint(grid_g(ifm)%flpu(i-1,j))
      k3u_1 = k2u_1 + 1
      !----- V points, need to average between j-jdim and j... ----------------------------!
      k2v   = nint(grid_g(ifm)%flpv(i,j))
      k3v   = k2v+1
      k2v_1 = nint(grid_g(ifm)%flpv(i,j-jdim))
      k3v_1 = k2v_1 + 1
      !----- W/T points, only the i,j points are needed...  -------------------------------!
      k2w   = nint(grid_g(ifm)%flpw(i,j))
      k1w   = k2w - 1
      k3w   = k2w + 1
      !------------------------------------------------------------------------------------!


      theta_mean = basic_g(ifm)%theta(2,i,j)
      rv_mean    = basic_g(ifm)%rv(2,i,j)
      rtp_mean   = basic_g(ifm)%rtp(2,i,j)

      if (co2_on) then
         co2p_mean = basic_g(ifm)%co2p(2,i,j)
      else 
         co2p_mean = co2con(1)
      end if

      !----- Wind speed, we must average so it falls in the same grid point. --------------!
      up_mean   = (basic_g(ifm)%up(2,i,j) + basic_g(ifm)%up(2,i-1,j     ))     * 0.5
      vp_mean   = (basic_g(ifm)%vp(2,i,j) + basic_g(ifm)%vp(2,i  ,j-jdim))     * 0.5

      !------------------------------------------------------------------------------------!
      !    Average the Exner function.  Exner function is split into reference and         !
      ! perturbation, and it is vertically staggered.                                      !
      !------------------------------------------------------------------------------------!
      exner_mean = ( basic_g(ifm)%pp (1,i,j) + basic_g(ifm)%pp (2,i,j)                     &
                   + basic_g(ifm)%pi0(1,i,j) + basic_g(ifm)%pi0(2,i,j))  * 0.5
      !------------------------------------------------------------------------------------!

   case (1)
      !------------------------------------------------------------------------------------!
      !     Adaptive height coordinate (shaved-eta).  Here we need to do a weighted aver-  !
      ! age using the lowest predicted points and the points avove them.                   !
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !    Find the lowest boundary levels, depending on which kind of variable            !
      ! we are averaging (staggered grid).                                                 !
      !------------------------------------------------------------------------------------!
      !----- U points, need to average between i-1 and i... -------------------------------!
      k2u   = nint(grid_g(ifm)%flpu(i,j))
      k3u   = k2u + 1
      k2u_1 = nint(grid_g(ifm)%flpu(i-1,j))
      k3u_1 = k2u_1 + 1
      !----- V points, need to average between j-jdim and j... ----------------------------!
      k2v   = nint(grid_g(ifm)%flpv(i,j))
      k3v   = k2v+1
      k2v_1 = nint(grid_g(ifm)%flpv(i,j-jdim))
      k3v_1 = k2v_1 + 1
      !----- W/T points, only the i,j points are needed...  -------------------------------!
      k2w   = nint(grid_g(ifm)%flpw(i,j))
      k1w   = k2w - 1
      k3w   = k2w + 1
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     The topography above ground should be also averaged.                           !
      !------------------------------------------------------------------------------------!
      topma_t = .25 * ( grid_g(ifm)%topma(i  ,j     ) + grid_g(ifm)%topma(i-1,j     )      &
                      + grid_g(ifm)%topma(i  ,j-jdim) + grid_g(ifm)%topma(i-1,j-jdim))
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Find the weights for lowest predicted (above ground) points, relative to       !
      ! points above them.                                                                 !
      !------------------------------------------------------------------------------------!
      wtw  = (zm(k2w) - topma_t) * dzt(k2w)
      wtu1 = grid_g(ifm)%aru(k2u_1,i-1,     j)  / grid_g(ifm)%aru(k3u_1,i-1,     j)
      wtu2 = grid_g(ifm)%aru(k2u  ,  i,     j)  / grid_g(ifm)%aru(k3u  ,  i,     j)
      wtv1 = grid_g(ifm)%arv(k2v_1,  i,j-jdim)  / grid_g(ifm)%arv(k3v_1,  i,j-jdim)
      wtv2 = grid_g(ifm)%arv(k2v,    i,     j)  / grid_g(ifm)%arv(k3v  ,  i,     j)

      theta_mean    =       wtw   * basic_g(ifm)%theta(k2w,i,j)                            &
                    + (1. - wtw)  * basic_g(ifm)%theta(k3w,i,j)

      rv_mean       =       wtw   * basic_g(ifm)%rv(k2w,i,j)                               &
                    + (1. - wtw)  * basic_g(ifm)%rv(k3w,i,j)

      rtp_mean      =       wtw   * basic_g(ifm)%rtp(k2w,i,j)                              &
                    + (1. - wtw)  * basic_g(ifm)%rtp(k3w,i,j)

      if (co2_on) then
         co2p_mean  =       wtw   * basic_g(ifm)%co2p(k2w,i,j)                             &
                    + (1. - wtw)  * basic_g(ifm)%co2p(k3w,i,j)
      else 
         co2p_mean  = co2con(1)
      end if

      !------------------------------------------------------------------------------------!
      !     Wind speed must use the staggered grid, but it also must take the volume above !
      ! ground.                                                                            !
      !------------------------------------------------------------------------------------!
      up_mean      = (        wtu1  * basic_g(ifm)%up(k2u_1,i-1,   j)                      &
                     +  (1. - wtu1) * basic_g(ifm)%up(k3u_1,i-1,   j)                      &
                     +        wtu2  * basic_g(ifm)%up(k2u,    i,   j)                      &
                     +  (1. - wtu2) * basic_g(ifm)%up(k3u,    i,   j)) * .5
      vp_mean      = (        wtv1  * basic_g(ifm)%vp(k2v_1,i,j-jdim)                      &
                     +  (1. - wtv1) * basic_g(ifm)%vp(k3v_1,i,j-jdim)                      &
                     +        wtv2  * basic_g(ifm)%vp(k2v,  i,     j)                      &
                     +  (1. - wtv2) * basic_g(ifm)%vp(k3v,  i,     j)) * .5
      
      !------------------------------------------------------------------------------------!
      !     Exner function and density, need to consider the layer beneath the ground to   !
      ! make the profile.                                                                  !
      !------------------------------------------------------------------------------------!
      if (wtw >= .5) then
         exner_mean  = (wtw - .5)                                                          &
                     * (basic_g(ifm)%pp(k1w,i,j) + basic_g(ifm)%pi0(k1w,i,j))              &
                     + (1.5 - wtw)                                                         &
                     * (basic_g(ifm)%pp(k2w,i,j) + basic_g(ifm)%pi0(k2w,i,j))
      else
         exner_mean  = (wtw + .5)                                                          &
                     * (basic_g(ifm)%pp(k2w,i,j) + basic_g(ifm)%pi0(k2w,i,j))              &
                     + (.5 - wtw)                                                          &
                     * (basic_g(ifm)%pp(k3w,i,j) + basic_g(ifm)%pi0(k3w,i,j))
      end if
   end select
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !   Find the reference height.                                                          !
   !---------------------------------------------------------------------------------------!
   zref_mean = (zt(k2w)-zm(k1w))*grid_g(ifm)%rtgt(i,j)
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Copy 2-D variables, that don't depend on the type of vertical coordinate.         !
   !---------------------------------------------------------------------------------------!
   lakemet%rshort   = dble(radiate_g(ifm)%rshort (i,j))
   lakemet%rlong    = dble(radiate_g(ifm)%rlong  (i,j))
   lakemet%tanz     = tan(acos(dble(radiate_g(ifm)%cosz  (i,j))))
   lakemet%lon      = dble(grid_g   (ifm)%glon   (i,j))
   lakemet%lat      = dble(grid_g   (ifm)%glat   (i,j))
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Copy the values to the meteorological buffer.  Start with those that require only !
   ! a simple copy with the conversion from single to double precision.                    !
   !---------------------------------------------------------------------------------------!
   lakemet%geoht     = dble(zref_mean )
   lakemet%atm_theta = dble(theta_mean)
   lakemet%atm_co2   = dble(co2p_mean )
   lakemet%atm_exner = dble(exner_mean)
   lakemet%atm_shv   = dble(rtp_mean  ) / (1.d0 + dble(rtp_mean))
   !----- SST derivative is already in double precision, just copy it. --------------------!
   lakemet%dsst_dt   = dsst_dt
   !---------------------------------------------------------------------------------------!
   



   !---------------------------------------------------------------------------------------!
   !     Sanity check.                                                                     !
   !---------------------------------------------------------------------------------------!
   ok_flpoint = is_finite8(lakemet%rshort)    .and. is_finite8(lakemet%rlong)    .and.     &
                is_finite8(lakemet%tanz)      .and. is_finite8(lakemet%lon)      .and.     &
                is_finite8(lakemet%lat)       .and. is_finite8(lakemet%geoht)    .and.     &
                is_finite8(lakemet%atm_theta) .and. is_finite8(lakemet%atm_co2)  .and.     &
                is_finite8(lakemet%atm_exner) .and. is_finite8(lakemet%atm_shv)  .and.     &
                is_finite8(lakemet%dsst_dt)   .and. is_finite (exner_mean)       .and.     &
                is_finite (theta_mean)        .and. is_finite (co2p_mean)        .and.     &
                is_finite (up_mean)           .and. is_finite (vp_mean)          .and.     &
                is_finite (rv_mean)           .and. is_finite (rtp_mean)         .and.     &
                is_finite (zref_mean)
   if (.not. ok_flpoint) then
      write(unit=*,fmt='(a)'          ) '-------------------------------------------------'
      write(unit=*,fmt='(a)'          ) '  Something went wrong...                        '
      write(unit=*,fmt='(a)'          ) '-------------------------------------------------'
      write(unit=*,fmt='(a)'          ) ' Meteorological forcing.'
      write(unit=*,fmt='(a,1x,es12.5)') ' - Rshort         :',lakemet%rshort
      write(unit=*,fmt='(a,1x,es12.5)') ' - Rlong          :',lakemet%rlong
      write(unit=*,fmt='(a,1x,es12.5)') ' - Tanz           :',lakemet%tanz
      write(unit=*,fmt='(a,1x,es12.5)') ' - Lon            :',lakemet%lon
      write(unit=*,fmt='(a,1x,es12.5)') ' - Lat            :',lakemet%lat
      write(unit=*,fmt='(a,1x,es12.5)') ' - Geoht          :',lakemet%geoht
      write(unit=*,fmt='(a,1x,es12.5)') ' - Theta          :',lakemet%atm_theta
      write(unit=*,fmt='(a,1x,es12.5)') ' - CO2            :',lakemet%atm_co2
      write(unit=*,fmt='(a,1x,es12.5)') ' - Exner          :',lakemet%atm_exner
      write(unit=*,fmt='(a,1x,es12.5)') ' - Spec. hum.     :',lakemet%atm_shv
      write(unit=*,fmt='(a,1x,es12.5)') ' - d(SST)/dt      :',lakemet%dsst_dt
      write(unit=*,fmt='(a)'          ) ' Mean values.'
      write(unit=*,fmt='(a,1x,es12.5)') ' - Exner_mean     :',exner_mean
      write(unit=*,fmt='(a,1x,es12.5)') ' - Theta_mean     :',theta_mean
      write(unit=*,fmt='(a,1x,es12.5)') ' - CO2p_mean      :',co2p_mean
      write(unit=*,fmt='(a,1x,es12.5)') ' - Up_mean        :',up_mean
      write(unit=*,fmt='(a,1x,es12.5)') ' - Vp_mean        :',vp_mean
      write(unit=*,fmt='(a,1x,es12.5)') ' - Rv_mean        :',rv_mean
      write(unit=*,fmt='(a,1x,es12.5)') ' - Rtp_mean       :',rtp_mean
      write(unit=*,fmt='(a,1x,es12.5)') ' - Zref_mean      :',zref_mean
      write(unit=*,fmt='(a)'          ) ' BRAMS values.'
      write(unit=*,fmt='(a,1x,i12)')    ' - Grid           :',ifm
      write(unit=*,fmt='(a,1x,i12)')    ' - I              :',i
      write(unit=*,fmt='(a,1x,i12)')    ' - J              :',j
      write(unit=*,fmt='(a,1x,i12)')    ' - JDIM           :',jdim
      write(unit=*,fmt='(a,1x,i12)')    ' - K1W            :',k1w
      write(unit=*,fmt='(a,1x,i12)')    ' - K2W            :',k2w
      write(unit=*,fmt='(a,1x,i12)')    ' - K3W            :',k3w
      write(unit=*,fmt='(a,1x,es12.5)') ' - Glon           :',grid_g(ifm)%glon(i,j)
      write(unit=*,fmt='(a,1x,es12.5)') ' - Glat           :',grid_g(ifm)%glat(i,j)
      write(unit=*,fmt='(a,1x,es12.5)') ' - Theta          :',basic_g(ifm)%theta(2,i,j)
      write(unit=*,fmt='(a,1x,es12.5)') ' - Rv             :',basic_g(ifm)%rv(2,i,j)
      write(unit=*,fmt='(a,1x,es12.5)') ' - Rtp            :',basic_g(ifm)%rtp(2,i,j)
      if (co2_on) then
         write(unit=*,fmt='(a,1x,es12.5)') ' - CO2p           :',basic_g(ifm)%co2p(2,i,j)
      else 
         write(unit=*,fmt='(a,1x,es12.5)') ' - CO2p           :',co2con(1)
      end if
      write(unit=*,fmt='(a,1x,es12.5)') ' - Up(i,j)        :',basic_g(ifm)%up(2,i,j)
      write(unit=*,fmt='(a,1x,es12.5)') ' - Up(i-1,j)      :',basic_g(ifm)%up(2,i-1,j)
      write(unit=*,fmt='(a,1x,es12.5)') ' - Vp(i,j)        :',basic_g(ifm)%vp(2,i,j)
      write(unit=*,fmt='(a,1x,es12.5)') ' - Vp(i,j-jdim)   :',basic_g(ifm)%vp(2,i,j-jdim)
      write(unit=*,fmt='(a,1x,es12.5)') ' - Pi0(1,i,j)     :',basic_g(ifm)%pi0(1,i,j)
      write(unit=*,fmt='(a,1x,es12.5)') ' - Pi0(2,i,j)     :',basic_g(ifm)%pi0(2,i,j)
      write(unit=*,fmt='(a,1x,es12.5)') ' - Pp(1,i,j)      :',basic_g(ifm)%pp(1,i,j)
      write(unit=*,fmt='(a,1x,es12.5)') ' - Pp(2,i,j)      :',basic_g(ifm)%pp(2,i,j)
      write(unit=*,fmt='(a,1x,es12.5)') ' - zt(2)          :',zt(k2w)
      write(unit=*,fmt='(a,1x,es12.5)') ' - zm(1)          :',zm(k1w)
      write(unit=*,fmt='(a,1x,es12.5)') ' - Rtgt           :',grid_g(ifm)%rtgt(i,j)
      write(unit=*,fmt='(a,1x,es12.5)') ' - Rshort         :',radiate_g(ifm)%rshort(i,j)
      write(unit=*,fmt='(a,1x,es12.5)') ' - Rlong          :',radiate_g(ifm)%rlong(i,j)
      write(unit=*,fmt='(a,1x,es12.5)') ' - Cosz           :',radiate_g(ifm)%cosz(i,j)
      write(unit=*,fmt='(a)'          ) ' '
      write(unit=*,fmt='(a)'          ) '-------------------------------------------------'
      call abort_run('Non-resolvable values','copy_met_2_lake','edcp_lake_misc.f90')
   end if
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     Copy the canopy air space properties to double precision scratch variables.       !
   !---------------------------------------------------------------------------------------!
   can_theta8 = dble(leaf_g(ifm)%can_theta(i,j,1))
   can_shv8   = dble(leaf_g(ifm)%can_rvap (i,j,1))                                         &
              / (1.d0 +dble(leaf_g(ifm)%can_rvap (i,j,1)))
   can_depth8 = dble(can_depth_min)
   !---------------------------------------------------------------------------------------!



   !----- Pressure. -----------------------------------------------------------------------!
   lakemet%atm_prss  = exner2press8(lakemet%atm_exner)
   !---------------------------------------------------------------------------------------!


   !----- Air temperature. ----------------------------------------------------------------!
   lakemet%atm_tmp   = extheta2temp8(lakemet%atm_exner,lakemet%atm_theta)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Find the pressure and Exner functions at the canopy depth, find the temperature   !
   ! of the air above canopy at the canopy depth, and the specific enthalpy at that level. !
   !---------------------------------------------------------------------------------------!
   can_prss8            = reducedpress8(lakemet%atm_prss,lakemet%atm_theta,lakemet%atm_shv &
                                       ,lakemet%geoht,can_theta8,can_shv8,can_depth8)
   can_exner8           = press2exner8 (can_prss8)
   lakemet%atm_tmp_zcan = extheta2temp8(can_exner8,lakemet%atm_theta)
   lakemet%atm_enthalpy = tq2enthalpy8 (lakemet%atm_tmp_zcan,lakemet%atm_shv,.true.)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Update properties that need to use therm_lib8.                                    !
   !---------------------------------------------------------------------------------------!
   lakemet%atm_rhos  = idealdenssh8(lakemet%atm_prss,lakemet%atm_tmp,lakemet%atm_shv)
   lakemet%atm_rhv   = rehuil8(lakemet%atm_prss,lakemet%atm_tmp,lakemet%atm_shv,.true.)
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Find the mean wind speed and mean wind direction, so we can split the flux of     !
   ! momentum.                                                                             !
   !---------------------------------------------------------------------------------------!
   lakemet%atm_vels = max(ubmin8,sqrt(dble(up_mean)**2 + dble(vp_mean)**2))
   angle            = datan2(dble(vp_mean),dble(up_mean))
   lakemet%ucos     = dcos(angle)
   lakemet%usin     = dsin(angle)
   !---------------------------------------------------------------------------------------!




   return
end subroutine copy_met_2_lake
!==========================================================================================!
!==========================================================================================!





!==========================================================================================!
!==========================================================================================!
subroutine copy_lake_brams(i,j,ifm,mzg,mzs,initp)
   use mem_basic             , only : basic_g      ! ! structure
   use mem_radiate           , only : radiate_g    ! ! structure
   use mem_leaf              , only : leaf_g       ! ! structure
   use lake_coms             , only : lakemet      ! ! intent(out)
   use consts_coms           , only : grav         ! ! intent(in)
   use canopy_air_coms       , only : ubmin8       ! ! intent(in)
   use lake_coms             , only : lakesitetype & ! structure
                                    , lakemet      & ! intent(in)
                                    , tiny_lakeoff ! ! intent(in)
   use therm_lib             , only : thetaeiv     & ! function
                                    , vpdefil      & ! function
                                    , press2exner  & ! function
                                    , extheta2temp ! ! function
   use therm_lib8            , only : alvl8        & ! function
                                    , alvi8        & ! function
                                    , tq2enthalpy8 & ! function
                                    , tl2uint8     ! ! function
   use mem_edcp              , only : ed_fluxf_g   ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer           , intent(in) :: i
   integer           , intent(in) :: j
   integer           , intent(in) :: ifm
   integer           , intent(in) :: mzg
   integer           , intent(in) :: mzs
   type(lakesitetype), target     :: initp
   !----- Local variables. ----------------------------------------------------------------!
   integer                        :: k
   real                           :: can_shv
   real                           :: can_exner
   real                           :: can_temp
   !----- External functions. -------------------------------------------------------------!
   real              , external   :: sngloff
   !----- Local constants -----------------------------------------------------------------!
   real             , parameter   :: z0fac_water  = 0.016/grav
   real             , parameter   :: z0_min_water = 0.0001
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !     Transfer model scalars back to global arrays.                                     !
   !---------------------------------------------------------------------------------------!
   !----- Stars. --------------------------------------------------------------------------!
   ed_fluxf_g(ifm)%ustar (i,j,1) = sngloff(initp%avg_ustar,tiny_lakeoff)
   ed_fluxf_g(ifm)%tstar (i,j,1) = sngloff(initp%avg_tstar,tiny_lakeoff)
   ed_fluxf_g(ifm)%cstar (i,j,1) = sngloff(initp%avg_cstar,tiny_lakeoff)
   ed_fluxf_g(ifm)%zeta  (i,j,1) = sngloff(initp%zeta     ,tiny_lakeoff)
   ed_fluxf_g(ifm)%ribulk(i,j,1) = sngloff(initp%ribulk   ,tiny_lakeoff)
   !----- Water vapour:  we must convert specific humidity back to mixing ratio. ----------!
   ed_fluxf_g(ifm)%rstar (i,j,1) = sngloff( initp%avg_qstar                                &
                                          / ( (1.d0-lakemet%atm_shv)                       &
                                            * (1.d0 - initp%can_shv))                      &
                                         , tiny_lakeoff)
   !---------------------------------------------------------------------------------------!



   !----- Fluxes towards the free atmosphere. ---------------------------------------------!
   ed_fluxf_g(ifm)%sflux_u(i,j,1) = sngloff(initp%avg_sflux_u ,tiny_lakeoff)
   ed_fluxf_g(ifm)%sflux_v(i,j,1) = sngloff(initp%avg_sflux_v ,tiny_lakeoff)
   ed_fluxf_g(ifm)%sflux_w(i,j,1) = sngloff(initp%avg_sflux_w ,tiny_lakeoff)
   ed_fluxf_g(ifm)%sflux_t(i,j,1) = sngloff(initp%avg_sflux_t ,tiny_lakeoff)
   ed_fluxf_g(ifm)%sflux_r(i,j,1) = sngloff(initp%avg_sflux_r ,tiny_lakeoff)
   ed_fluxf_g(ifm)%sflux_c(i,j,1) = sngloff(initp%avg_sflux_c ,tiny_lakeoff)
   !---------------------------------------------------------------------------------------!



   !----- Radiation-related variables. ----------------------------------------------------!
   ed_fluxf_g(ifm)%albedt     (i,j,1) = sngloff(initp%avg_albedt     ,tiny_lakeoff)
   ed_fluxf_g(ifm)%rlongup    (i,j,1) = sngloff(initp%avg_rlongup    ,tiny_lakeoff)
   ed_fluxf_g(ifm)%rshort_gnd (i,j,1) = sngloff(initp%avg_rshort_gnd ,tiny_lakeoff)
   ed_fluxf_g(ifm)%rlong_gnd  (i,j,1) = 0.0
   !---------------------------------------------------------------------------------------!



   !----- Find the flux components of each patch. -----------------------------------------!
   leaf_g(ifm)%hflxac      (i,j,1) = sngloff( initp%avg_sensible_ac, tiny_lakeoff)
   leaf_g(ifm)%wflxac      (i,j,1) = sngloff( initp%avg_vapor_ac   , tiny_lakeoff)
   leaf_g(ifm)%qwflxac     (i,j,1) = sngloff( initp%avg_vapor_ac                           &
                                            * tq2enthalpy8(initp%can_temp,1.d0,.true.)     &
                                            , tiny_lakeoff )
   leaf_g(ifm)%eflxac      (i,j,1) = leaf_g(ifm)%hflxac (i,j,1)                            &
                                   + leaf_g(ifm)%qwflxac(i,j,1)
   leaf_g(ifm)%cflxac      (i,j,1) = sngloff( initp%avg_carbon_ac  , tiny_lakeoff)
   leaf_g(ifm)%hflxgc      (i,j,1) = sngloff( initp%avg_sensible_gc, tiny_lakeoff)
   leaf_g(ifm)%wflxgc      (i,j,1) = sngloff( initp%avg_vapor_gc   , tiny_lakeoff)
   leaf_g(ifm)%qwflxgc     (i,j,1) = sngloff( initp%avg_vapor_gc                           &
                                            * tq2enthalpy8(initp%lake_temp,1.d0,.true.)    &
                                            , tiny_lakeoff)
   leaf_g(ifm)%hflxvc      (i,j,1) = 0.0
   leaf_g(ifm)%wflxvc      (i,j,1) = 0.0
   leaf_g(ifm)%qwflxvc     (i,j,1) = 0.0
   leaf_g(ifm)%transp      (i,j,1) = 0.0
   leaf_g(ifm)%qtransp     (i,j,1) = 0.0
   leaf_g(ifm)%intercepted (i,j,1) = 0.0
   leaf_g(ifm)%qintercepted(i,j,1) = 0.0
   leaf_g(ifm)%wshed       (i,j,1) = 0.0
   leaf_g(ifm)%qwshed      (i,j,1) = 0.0
   leaf_g(ifm)%throughfall (i,j,1) = 0.0
   leaf_g(ifm)%qthroughfall(i,j,1) = 0.0
   leaf_g(ifm)%gpp         (i,j,1) = 0.0
   leaf_g(ifm)%plresp      (i,j,1) = 0.0
   leaf_g(ifm)%resphet     (i,j,1) = 0.0
   leaf_g(ifm)%growresp    (i,j,1) = 0.0
   leaf_g(ifm)%runoff      (i,j,1) = 0.0
   leaf_g(ifm)%qrunoff     (i,j,1) = 0.0
   leaf_g(ifm)%drainage    (i,j,1) = 0.0
   leaf_g(ifm)%qdrainage   (i,j,1) = 0.0
   !---------------------------------------------------------------------------------------!


   !----- Finding some canopy air properties. ---------------------------------------------!
   leaf_g(ifm)%can_theta(i,j,1)    = sngloff(initp%can_theta ,tiny_lakeoff)
   can_shv                         = sngloff(initp%can_shv   ,tiny_lakeoff)
   leaf_g(ifm)%can_rvap(i,j,1)     = can_shv / (1.0 - can_shv)
   leaf_g(ifm)%can_co2(i,j,1)      = sngloff(initp%can_co2   ,tiny_lakeoff)
   leaf_g(ifm)%can_prss(i,j,1)     = sngloff(initp%can_prss  ,tiny_lakeoff)
   can_exner                       = press2exner(leaf_g(ifm)%can_prss(i,j,1))
   can_temp                        = extheta2temp(can_exner,leaf_g(ifm)%can_theta(i,j,1))
   leaf_g(ifm)%can_theiv(i,j,1)    = thetaeiv( leaf_g(ifm)%can_theta(i,j,1)                &
                                             , leaf_g(ifm)%can_prss(i,j,1)                 &
                                             , can_temp                                    &
                                             , leaf_g(ifm)%can_rvap(i,j,1)                 &
                                             , leaf_g(ifm)%can_rvap(i,j,1)                 )
   leaf_g(ifm)%can_vpdef(i,j,1)    = vpdefil ( leaf_g(ifm)%can_prss(i,j,1)                 &
                                             , can_temp                                    &
                                             , leaf_g(ifm)%can_rvap(i,j,1)                 &
                                             , .false.                                     )
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Unless one includes a PFT for water lilies, there is no vegetation on water...    !
   !---------------------------------------------------------------------------------------!
   leaf_g(ifm)%veg_energy(i,j,1) = 0.
   leaf_g(ifm)%veg_water (i,j,1) = 0.
   leaf_g(ifm)%veg_lai   (i,j,1) = 0.
   leaf_g(ifm)%veg_tai   (i,j,1) = 0.
   leaf_g(ifm)%veg_agb   (i,j,1) = 0.
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !    Our simple lake model allows no ice, making sure the values are properly zeroed.   !
   !---------------------------------------------------------------------------------------!
   leaf_g(ifm)%sfcwater_nlev     (i,j,1) = 0.
   do k=1, mzs
      leaf_g(ifm)%sfcwater_energy (1,i,j,1) = 0.
      leaf_g(ifm)%sfcwater_mass   (1,i,j,1) = 0.
      leaf_g(ifm)%sfcwater_depth  (1,i,j,1) = 0.
   end do
   !---------------------------------------------------------------------------------------!


   !----- Find roughness scales for water bodies ------------------------------------------!
   leaf_g(ifm)%patch_rough(i,j,1) = max(z0fac_water * leaf_g(ifm)%ustar(i,j,1) ** 2        &
                                       ,z0_min_water)
   leaf_g(ifm)%soil_rough (i,j,1) = 0.0
   leaf_g(ifm)%veg_rough  (i,j,1) = 0.0
   !---------------------------------------------------------------------------------------!





   !---------------------------------------------------------------------------------------!
   !    "Soil" energy.  Because we can't store sea surface temperature, we store the       !
   ! internal energy
   !---------------------------------------------------------------------------------------!
   leaf_g(ifm)%soil_energy (mzg,i,j,1) = sngloff(tl2uint8(initp%lake_temp,initp%lake_fliq) &
                                                ,tiny_lakeoff)
   leaf_g(ifm)%soil_water  (mzg,i,j,1) = 0.
   do k=1, mzg-1
      leaf_g(ifm)%soil_energy (k,i,j,1) = leaf_g(ifm)%soil_energy (mzg,i,j,1)
      leaf_g(ifm)%soil_water  (k,i,j,1) = leaf_g(ifm)%soil_water  (mzg,i,j,1)
   end do
   leaf_g(ifm)%psibar_10d       (i,j,1) = 1.0
   !---------------------------------------------------------------------------------------!

   return
end subroutine copy_lake_brams
!==========================================================================================!
!==========================================================================================!
