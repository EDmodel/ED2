subroutine canopy_update_euler(csite, ipa, vels, atm_tmp, prss, pcpg, qpcpg,  &
     wshed_canopy, qwshed_canopy, canwcap, canhcap, dt_leaf, hxfergc,  &
     sxfer_t, wxfergc, hxfersc, wxfersc, sxfer_r, ed_transp,          &
     ntext_soil, soil_water, soil_fracliq, lsl,  &
     leaf_aging_factor, green_leaf_factor, sxfer_c)


  use ed_state_vars,only:sitetype,patchtype
  use grid_coms, only: nzg
  use consts_coms, only: alvi
  use ed_max_dims, only: n_pft

  implicit none

  type(sitetype),target   :: csite
  type(patchtype),pointer :: cpatch
  integer :: ipa,ico
  real, intent(in) :: vels
  real, intent(in) :: atm_tmp
  real, intent(in) :: prss
  real, intent(in) :: pcpg
  real, intent(in) :: qpcpg
  real, intent(inout) :: wshed_canopy
  real, intent(inout) :: qwshed_canopy
  real, intent(in) :: canwcap
  real, intent(in) :: canhcap
  real, intent(in) :: dt_leaf
  real, intent(in) :: hxfergc
  real, intent(in) :: sxfer_t
  real, intent(in) :: sxfer_c
  real, intent(in) :: wxfergc
  real, intent(in) :: hxfersc
  real, intent(in) :: wxfersc
  real, intent(in) :: sxfer_r
  real, intent(in), dimension(n_pft) :: leaf_aging_factor, green_leaf_factor
  real, dimension(nzg) :: ed_transp

  
  real :: sum_lai_rbi
  integer :: ndims
  integer, dimension(nzg) :: ed_ktrans
  integer, dimension(nzg), intent(in) :: ntext_soil
  real, dimension(nzg), intent(in) :: soil_water
  real, dimension(nzg), intent(in) :: soil_fracliq
  integer, intent(in) :: lsl

  ! Get photosynthesis, stomatal conductance, and transpiration
  call canopy_photosynthesis(csite,ipa, vels, atm_tmp, prss, &
       ed_ktrans, ntext_soil, soil_water, soil_fracliq, lsl, sum_lai_rbi, &
       leaf_aging_factor, green_leaf_factor)

  call canopy_precip_interception(csite,ipa, pcpg, qpcpg, wshed_canopy,  &
       qwshed_canopy)
  
  ndims = 2
  cpatch => csite%patch(ipa)
  do ico = 1,cpatch%ncohorts
     if(cpatch%solvable(ico) )then
        ndims = ndims + 2
     endif
  enddo

! IMPLICIT SCHEME STILL NEEDS UPDATING.
!  call canopy_implicit_driver(ed_patch, ndims, prss, canhcap, canair,   &
!       sum_lai_rbi, dt_leaf, hxfergc, hxferca, wxfergc, hxfersc, wxfersc,  &
  !       sxfer_r, ed_transp, ed_ktrans)
  call canopy_explicit_driver(csite, ipa, ndims, prss, canhcap, canwcap,   &
       sum_lai_rbi, dt_leaf, hxfergc, sxfer_t, wxfergc, hxfersc, wxfersc,  &
       sxfer_r, ed_transp, ed_ktrans, sxfer_c)

!  ed_patch%mean_wflux = ed_patch%mean_wflux + (sxfer_r) / dt_leaf
!  ed_patch%mean_latflux = ed_patch%mean_latflux + (wxfergc + wxfersc) /   &
!       dt_leaf * alvi
!  ed_patch%mean_hflux = ed_patch%mean_hflux + hxferca / dt_leaf

  return
end subroutine canopy_update_euler

!==============================================================

subroutine canopy_precip_interception(csite,ipa, pcpg, qpcpg, wshed_canopy,  &
     qwshed_canopy)

  use ed_state_vars,only:sitetype,patchtype
  use canopy_radiation_coms, only: lai_min
  use consts_coms, only: cice, cliq, alli, t3ple, tsupercool
  use therm_lib, only: qwtk

  implicit none

  type(sitetype),target :: csite
  type(patchtype),pointer :: cpatch
  integer :: ipa,ico
  real, intent(in) :: pcpg
  real, intent(in) :: qpcpg
  real, intent(inout) :: wshed_canopy
  real, intent(inout) :: qwshed_canopy

  real :: laii

  real :: tvegaux
  real :: qwtot
  real :: lai_fraction
  real :: wshed_layer
  real :: hcapveg_factor

  laii = 0.0
  cpatch => csite%patch(ipa)

  do ico = 1,cpatch%ncohorts
     if(cpatch%hite(ico) > csite%total_snow_depth(ipa))then
        laii = laii + cpatch%lai(ico)
     endif
  enddo

  if (laii > lai_min) then
     
     laii = 1.0 / laii
     
     do ico = 1,cpatch%ncohorts

        if(cpatch%solvable(ico))then
           tvegaux = cpatch%veg_temp(ico) - t3ple

           hcapveg_factor = cpatch%lai(ico) / csite%lai(ipa)

! If precipitation, add intercepted mass and energy to vegetation surface
!   (Ignore energy of water already on vegetation surface)

! Vegetation layers intercept precipitation in proportion to their LAI
           lai_fraction = cpatch%lai(ico) * laii
           !MLO - I guess we can use vegetation energy instead of recalculating qwtot
           !      from temperature. Please, let me know whether there is any problem 
           !      with this...
           
           !if(tvegaux > 0.0)then
           !   qwtot = cpatch%hcapveg(ico) * hcapveg_factor * tvegaux + qpcpg *   &
           !        lai_fraction + cpatch%veg_water(ico) * (cliq * tvegaux + alli)
           !else
           !   qwtot = cpatch%hcapveg(ico) * hcapveg_factor * tvegaux +   &
           !        qpcpg * lai_fraction + cpatch%veg_water(ico) * cice * tvegaux
           !endif
           cpatch%veg_energy(ico) = cpatch%veg_energy(ico) + qpcpg * lai_fraction
           cpatch%veg_water(ico)  = cpatch%veg_water(ico) + pcpg * lai_fraction

! Compute equilbrium temperature of veg + precipitation

           call qwtk(cpatch%veg_energy(ico), cpatch%veg_water(ico), cpatch%hcapveg(ico) * hcapveg_factor,   &
                cpatch%veg_temp(ico), cpatch%veg_fliq(ico))
      
! Shed any excess intercepted precipitation and its energy

           if (cpatch%veg_water(ico) > .22 * cpatch%lai(ico)) then
              wshed_layer = cpatch%veg_water(ico) - .22 * cpatch%lai(ico)
              wshed_canopy = wshed_canopy + wshed_layer

              if (cpatch%veg_fliq(ico) <= .0001) then
                 qwshed_canopy = qwshed_canopy + cice * cpatch%veg_temp(ico) * wshed_layer
              else
                 qwshed_canopy = qwshed_canopy + wshed_layer &
                               * ( cpatch%veg_fliq(ico)*cliq*(cpatch%veg_temp(ico)-tsupercool) &
                                 + (1.-cpatch%veg_fliq(ico))*cice*cpatch%veg_temp(ico))
              endif
              
              cpatch%veg_water(ico) = cpatch%veg_water(ico) - wshed_layer
           endif
           
        endif

     enddo
  else
     wshed_canopy = pcpg
     qwshed_canopy = qpcpg
  endif

  return
end subroutine canopy_precip_interception

!==============================================================

subroutine canopy_implicit_driver(csite,ipa, ndims, prss, canhcap, canwcap,  &
     sum_lai_rbi, dt_leaf, hxfergc, sxfer_t, wxfergc, hxfersc, wxfersc,    &
     sxfer_r, ed_transp, ed_ktrans, sxfer_c)

!-----------------------------------------------------------------------------
! Execute the implicit exchange of heat and moisture between vegetation and 
! canopy air
!-----------------------------------------------------------------------------

  use ed_state_vars,only: sitetype,patchtype
  use consts_coms, only: cp, alvl, alvi,rdry
  use grid_coms, only: nzg
  use therm_lib, only: rhovsil,rhovsilp,tq2enthalpy,idealdenssh
  implicit none

  type(sitetype), target  :: csite
  type(patchtype),pointer :: cpatch
  integer :: ipa,ico
  integer, intent(in) :: ndims
  real, intent(in) :: prss
  real, intent(in) :: canhcap
  real, intent(in) :: canwcap
  real, intent(in) :: sum_lai_rbi
  real, intent(in) :: dt_leaf
  real, intent(in) :: hxfergc
  real, intent(in) :: wxfergc
  real, intent(in) :: hxfersc
  real, intent(in) :: wxfersc
  real, intent(in) :: sxfer_c
  real, intent(in) :: sxfer_r
  real, intent(in) :: sxfer_t
  real, dimension(nzg), intent(out) :: ed_transp
  integer, dimension(nzg), intent(in) :: ed_ktrans

  ! Locals

  real, dimension(ndims, ndims) :: t_evolve_matrix
  integer :: ic
  real :: tvegaux
  real :: veg_rhovs
  real :: veg_rhovsp
  real :: vp_gradient
  real :: a1
  real :: sigmaw
  real :: dsigmaw_dW
  real :: a4
  real :: et_conductance
  real, dimension(ndims) :: explicit_deriv_portion
  real, dimension(ndims) :: original_state
  integer :: ind1
  integer :: ind2
  integer, dimension(ndims) :: indx
  real :: d
  real, dimension(ndims) :: implicit_new_state
  integer :: k
  real :: mult
  real :: a10
  integer :: idim


  cpatch => csite%patch(ipa)

  ! Initialize time evolution matrix and the explicit contribution
  t_evolve_matrix(1:ndims, 1:ndims) = 0.0
  explicit_deriv_portion(1:ndims) = 0.0
  ed_transp(:) = 0.0

  ! explicitly integrated contribution to the canopy air temperature
  explicit_deriv_portion(1) = (hxfergc + hxfersc - sxfer_t) /   &
       (dt_leaf * canhcap)

  ! Canopy air temperature  (d can_temp, Tcan)
  t_evolve_matrix(1,1) = - 2.2 * cp * csite%can_rhos(ipa) / canhcap * sum_lai_rbi

  ! explicitly integrated contribution to the canopy specific humidity
  explicit_deriv_portion(2) = (wxfergc + wxfersc - sxfer_r) /   &
       (dt_leaf * canwcap)

  ic = 0

  do ico = 1,cpatch%ncohorts
     if(cpatch%solvable(ico))then

        ! Set indices
        ic = ic + 1
        ind1 = 1 + 2 * ic
        ind2 = 2 + 2 * ic

        ! Compute heat transfer coefficient
        a10 = 2.2 * cp * csite%can_rhos(ipa) * cpatch%lai(ico) / cpatch%rb(ico)

        ! compute ET using variables at time n
        tvegaux = cpatch%veg_temp(ico)
        veg_rhovs= rhovsil(tvegaux)
        veg_rhovsp = rhovsilp(tvegaux)
        vp_gradient = veg_rhovs - csite%can_rhos(ipa) * csite%can_shv(ipa)
        a1 = 2.2 * cpatch%lai(ico) / cpatch%rb(ico)
        sigmaw = min(1.,(cpatch%veg_water(ico) / (.22 * cpatch%lai(ico)))**.66667)
        if(cpatch%veg_water(ico) > 0.0)then
           dsigmaw_dW = 1.0 / (0.33 * cpatch%lai(ico) * sqrt(sigmaw))
        else
           dsigmaw_dW = 0.0
        endif
        a4 = cpatch%lai(ico) / (cpatch%rb(ico) + cpatch%stomatal_resistance(ico))

        if(vp_gradient <= 0.0)then

           et_conductance = a1

           ! d veg_water, explicit
           explicit_deriv_portion(ind2) = - vp_gradient * a1

           ! d veg_water, can_shv
           mult = csite%can_rhos(ipa) * a1
           t_evolve_matrix(ind2,2) = mult
           explicit_deriv_portion(ind2) = explicit_deriv_portion(ind2) -   &
                mult * csite%can_shv(ipa)

           ! d veg_water, veg_temp
           mult = - veg_rhovsp * a1
           t_evolve_matrix(ind2,ind1) = mult
           explicit_deriv_portion(ind2) = explicit_deriv_portion(ind2) -   &
                mult * cpatch%veg_temp(ico)

        else

           et_conductance = a1 * sigmaw + a4 * (1.0 - sigmaw)

           ! d can_shv, veg_water
           mult = vp_gradient * (a1 - a4) * dsigmaw_dW / canwcap
           t_evolve_matrix(2,ind2) = t_evolve_matrix(2,ind2) + mult
           explicit_deriv_portion(2) = explicit_deriv_portion(2) - mult *   &
                cpatch%veg_water(ico)

           ! d veg_temp, veg_water
           mult = -alvl / cpatch%hcapveg(ico) * vp_gradient * (a1 - a4) * dsigmaw_dW
           t_evolve_matrix(ind1,ind2) = mult
           explicit_deriv_portion(ind1) = explicit_deriv_portion(ind1) -   &
                cpatch%veg_water(ico) * mult

           ! d veg_water, explicit
           explicit_deriv_portion(ind2) = - vp_gradient * sigmaw * a1

           ! d veg_water, can_shv
           mult = csite%can_rhos(ipa) * sigmaw * a1
           t_evolve_matrix(ind2,2) = mult
           explicit_deriv_portion(ind2) = explicit_deriv_portion(ind2) -   &
                mult * csite%can_shv(ipa)

           ! d veg_water, veg_temp
           mult = - veg_rhovsp * sigmaw * a1
           t_evolve_matrix(ind2,ind1) = mult
           explicit_deriv_portion(ind2) = explicit_deriv_portion(ind2) -   &
                mult * cpatch%veg_temp(ico)

           ! d veg_water, veg_water
           mult = - vp_gradient * a1 * dsigmaw_dW
           t_evolve_matrix(ind2,ind2) = mult
           explicit_deriv_portion(ind2) = explicit_deriv_portion(ind2) -   &
                mult * cpatch%veg_water(ico)
           
        endif

        ! d can_temp, veg_temp
        t_evolve_matrix(1,ind1) = a10 / canhcap

        ! dcan_shv, explicit
        explicit_deriv_portion(2) = explicit_deriv_portion(2) +   &
             vp_gradient * et_conductance / canwcap

        ! dcan_shv, can_shv
        mult = - csite%can_rhos(ipa) / canwcap * et_conductance
        t_evolve_matrix(2,2) = t_evolve_matrix(2,2) + mult
        explicit_deriv_portion(2) = explicit_deriv_portion(2) - mult *   &
             csite%can_shv(ipa)

        ! dcan_shv, veg_temp
        mult = veg_rhovsp * et_conductance / canwcap
        t_evolve_matrix(2,ind1) = t_evolve_matrix(2,ind1) + mult
        explicit_deriv_portion(2) = explicit_deriv_portion(2) - mult *   &
             cpatch%veg_temp(ico)

        ! dveg_temp, explicit
        explicit_deriv_portion(ind1) = (cpatch%rshort_v(ico) + cpatch%rlong_v(ico)) /   &
             cpatch%hcapveg(ico) - a10 / cpatch%hcapveg(ico) * (cpatch%veg_temp(ico) -   &
             csite%can_temp(ipa)) - alvl * et_conductance * vp_gradient /   &
             cpatch%hcapveg(ico)

        ! d veg_temp, can_temp
        mult = a10 / cpatch%hcapveg(ico)
        t_evolve_matrix(ind1,1) = mult
        explicit_deriv_portion(ind1) = explicit_deriv_portion(ind1) -   &
             csite%can_temp(ipa) * mult

        ! d veg_temp, can_shv
        mult = csite%can_rhos(ipa) * alvl / cpatch%hcapveg(ico) * et_conductance
        t_evolve_matrix(ind1,2) = mult
        explicit_deriv_portion(ind1) = explicit_deriv_portion(ind1) -   &
             csite%can_shv(ipa) * mult

        ! d veg_temp, veg_temp
        mult = - a10 / cpatch%hcapveg(ico) - alvl / cpatch%hcapveg(ico) * veg_rhovsp *   &
             et_conductance
        t_evolve_matrix(ind1,ind1) = mult
        explicit_deriv_portion(ind1) = explicit_deriv_portion(ind1) -   &
             cpatch%veg_temp(ico) * mult

        ! derivative matrix is done.  Now load initial state of the cohorts:

        original_state(ind1) = cpatch%veg_temp(ico)
        original_state(ind2) = cpatch%veg_water(ico)

     endif

  enddo

  ! original state of canopy
  original_state(1) = csite%can_temp(ipa)
  original_state(2) = csite%can_shv(ipa)
  
  if(ic > 0)then
     
     ! Compute (Identity matrix) - dt_leaf * (Derivs matrix)
     do ic=1,ndims
        do k=1,ndims
           if(ic == k)then
              t_evolve_matrix(ic,k) = 1.0 - dt_leaf * t_evolve_matrix(ic,k)
           else
              t_evolve_matrix(ic,k) = - dt_leaf * t_evolve_matrix(ic,k)
           endif
        enddo
     enddo
     
     ! Do the matrix inversion
     call ludcmp_dble(t_evolve_matrix,ndims,ndims,indx,d)
     
  endif
     
  ! This is the vector that the inverse matrix needs to multiply
  implicit_new_state(1:ndims) = original_state(1:ndims) +   &
       explicit_deriv_portion(1:ndims) * dt_leaf

  ! Do the multiplication
  if(ic > 0)call lubksb_dble(t_evolve_matrix,ndims,ndims,indx,  &
       implicit_new_state)
     
  ! Load the new state into the patch and cohort structures
  csite%can_temp(ipa) = implicit_new_state(1)
  csite%can_shv(ipa) = implicit_new_state(2)
  idim = 1
  ic = 0
  
  do ico = 1,cpatch%ncohorts

     if(cpatch%solvable(ico))then
        ic = ic + 1
        idim = idim + 2
        cpatch%veg_temp(ico) = implicit_new_state(idim)
        cpatch%veg_water(ico) = implicit_new_state(idim+1)
        
        veg_rhovs= rhovsil(original_state(idim))
        veg_rhovsp = rhovsilp(original_state(idim))
        vp_gradient = veg_rhovs - csite%can_rhos(ipa) * original_state(2)
        if(vp_gradient > 0.0)then
          ! Calculate the resulting transpiration
           sigmaw = min(1.,(original_state(idim+1) / (.22 * cpatch%lai(ico)))**.66667)
           if(original_state(idim+1) > 0.0)then
              dsigmaw_dW = 1.0 / (0.33 * cpatch%lai(ico) * sqrt(sigmaw))
           else
              dsigmaw_dW = 0.0
           endif
           a4 = cpatch%lai(ico) / (cpatch%rb(ico) + cpatch%stomatal_resistance(ico))
           ed_transp(ed_ktrans(cpatch%krdepth(ico))) =   &
                ed_transp(ed_ktrans(cpatch%krdepth(ico))) +  &
                vp_gradient * (1.0-sigmaw) * a4 - &
                csite%can_rhos(ipa) * (1.0 - sigmaw) * a4 * (csite%can_shv(ipa) -   &
                original_state(2)) + veg_rhovsp * (1.0 - sigmaw) * a4 *   &
                (cpatch%veg_temp(ico) - original_state(idim))  - vp_gradient *   &
                dsigmaw_dW * a4 * (cpatch%veg_water(ico) - original_state(idim+1))
        endif

     endif

  enddo

  ed_transp(:) = ed_transp(:) * dt_leaf

  csite%mean_latflux(ipa) = csite%mean_latflux(ipa) +   &
       ((csite%can_shv(ipa) - original_state(2)) * canwcap -  & 
       (wxfergc + wxfersc - sxfer_r)) / dt_leaf * alvl

  ! Require veg_water >= 0
  do ico = 1,cpatch%ncohorts
     if(cpatch%solvable(ico))then
        if(cpatch%veg_water(ico) < 0.0)then
           
           ! Take away from the canopy humidity...
           csite%can_shv(ipa) = csite%can_shv(ipa)+ cpatch%veg_water(ico) / canwcap
           
           ! ... and set the veg_water to zero
           cpatch%veg_water(ico) = 0.0
           
        endif
     endif
     
  enddo

  if(csite%can_temp(ipa) /= csite%can_temp(ipa) .or.   &
       csite%can_temp(ipa) > 400.0 .or.   &
       csite%can_temp(ipa) < 100.0 .or. csite%can_shv(ipa) < 0.0001)then
     print*,'can temp is bad'
     print*,csite%can_temp(ipa),original_state(1),explicit_deriv_portion(1)
     print*,csite%can_shv(ipa),original_state(2),explicit_deriv_portion(2)
     print*,original_state(1:ndims)
     print*,hxfergc/canhcap,hxfersc/canhcap,sxfer_t/canhcap,csite%ustar(ipa)
     stop
  endif
  
  csite%can_rhos(ipa)       = idealdenssh(prss,csite%can_temp(ipa),csite%can_shv(ipa))
  csite%can_enthalpy(ipa)   = tq2enthalpy(csite%can_temp(ipa),csite%can_shv(ipa))
  csite%avg_daily_temp(ipa) = csite%avg_daily_temp(ipa) + csite%can_temp(ipa)
  return
end subroutine canopy_implicit_driver
!==============================================================

subroutine canopy_explicit_driver(csite,ipa, ndims, prss, canhcap, canwcap, &
     sum_lai_rbi, dt_leaf, hxfergc, sxfer_t, wxfergc, hxfersc, wxfersc,     &
     sxfer_r, ed_transp, ed_ktrans, sxfer_c)

!-----------------------------------------------------------------------------
! Execute the explicit exchange of heat and moisture between vegetation and 
! canopy air
!-----------------------------------------------------------------------------

  use ed_state_vars,only:sitetype,patchtype
  use consts_coms, only: cp, alvl, alvi, cliq, cice, alli, t3ple,tsupercool,rdry
  use grid_coms, only: nzg
  use therm_lib, only: rhovsil,rhovsilp,tq2enthalpy,idealdenssh

  implicit none

  type(patchtype),pointer :: cpatch
  type(sitetype),target :: csite
  integer :: ipa,ico

  integer, intent(in) :: ndims
  real, intent(in) :: prss
  real, intent(in) :: canhcap
  real, intent(in) :: canwcap
  real, intent(in) :: sum_lai_rbi
  real, intent(in) :: dt_leaf
  real, intent(in) :: hxfergc
  real, intent(in) :: wxfergc
  real, intent(in) :: hxfersc
  real, intent(in) :: wxfersc
  real, intent(in) :: sxfer_t
  real, intent(in) :: sxfer_r
  real, intent(in) :: sxfer_c
  real, dimension(nzg), intent(out) :: ed_transp
  integer, dimension(nzg), intent(in) :: ed_ktrans

  ! Locals

  integer :: ic

  real :: tvegaux
  real :: veg_rhovs
  real :: veg_rhovsp
  real :: vp_gradient
  real :: a1
  real :: sigmaw
  real :: dsigmaw_dW
  real :: a4
  real :: et_conductance
  real, dimension(ndims) :: explicit_deriv_portion
  real, dimension(ndims) :: original_state
  integer :: ind1
  integer :: ind2
  real :: a10
  real, dimension(ndims) :: explicit_new_state
  integer :: idim
  real :: dQdt

  ! Initialize time evolution matrix and the explicit contribution
  explicit_deriv_portion(1:ndims) = 0.0
  ed_transp(:) = 0.0

  ! explicitly integrated contribution to the canopy air temperature
  explicit_deriv_portion(1) = (hxfergc + hxfersc - sxfer_t) /   &
       (dt_leaf * canhcap)

  ! explicitly integrated contribution to the canopy specific humidity
  explicit_deriv_portion(2) = (wxfergc + wxfersc - sxfer_r) /   &
       (dt_leaf * canwcap)

  cpatch => csite%patch(ipa)


  ic = 0
  do ico = 1,cpatch%ncohorts
     if(cpatch%solvable(ico))then

        ! Set indices
        ic = ic + 1
        ind1 = 1 + 2 * ic
        ind2 = 2 + 2 * ic

        ! Compute heat transfer coefficient
        a10 = 2.2 * cp * csite%can_rhos(ipa) * cpatch%lai(ico) / cpatch%rb(ico)

        ! compute ET using variables at time n
        tvegaux = cpatch%veg_temp(ico)
        veg_rhovs= rhovsil(tvegaux)
        veg_rhovsp = rhovsilp(tvegaux)
        vp_gradient = veg_rhovs - csite%can_rhos(ipa) * csite%can_shv(ipa)
        a1 = 2.2 * cpatch%lai(ico) / cpatch%rb(ico)
        sigmaw = min(1.,(cpatch%veg_water(ico) / (.22 * cpatch%lai(ico)))**.66667)
        if(cpatch%veg_water(ico) > 0.0)then
           dsigmaw_dW = 1.0 / (0.33 * cpatch%lai(ico) * sqrt(sigmaw))
        else
           dsigmaw_dW = 0.0
        endif
        a4 = cpatch%lai(ico) / (cpatch%rb(ico) + cpatch%stomatal_resistance(ico))

        if(vp_gradient <= 0.0)then

           et_conductance = a1

           ! d veg_water, explicit
           explicit_deriv_portion(ind2) = - vp_gradient * a1

        else

           et_conductance = a1 * sigmaw + a4 * (1.0 - sigmaw)

           ! d veg_water, explicit
           explicit_deriv_portion(ind2) = - vp_gradient * sigmaw * a1
           ed_transp(ed_ktrans(cpatch%krdepth(ico))) =   &
                ed_transp(ed_ktrans(cpatch%krdepth(ico))) +  &
                vp_gradient * (1.0 - sigmaw) * a4

        endif

        ! contribution to dcan_temp/dt, from this vegetation layer's veg_temp
        explicit_deriv_portion(1) = explicit_deriv_portion(1) + a10 /  &
             canhcap * (cpatch%veg_temp(ico) - csite%can_temp(ipa))

        ! contribution to dcan_shv/dt, from this vegetation layer's veg_water
        explicit_deriv_portion(2) = explicit_deriv_portion(2) +   &
             vp_gradient * et_conductance / canwcap

        ! dveg_temp, explicit
        dQdt = cpatch%rshort_v(ico) + cpatch%rlong_v(ico) -   &
             a10 * (cpatch%veg_temp(ico) - csite%can_temp(ipa)) -   &
             alvl * et_conductance * vp_gradient
        tvegaux = cpatch%veg_temp(ico) - t3ple
        if(tvegaux > 0.)then
           csite%mean_latflux(ipa) = csite%mean_latflux(ipa) +   &
                vp_gradient * et_conductance * alvl -   &
                (cliq *(cpatch%veg_temp(ico)  -tsupercool)) * explicit_deriv_portion(ind2)
           explicit_deriv_portion(ind1) = dQdt / (cpatch%hcapveg(ico) * cpatch%lai(ico) /   &
                csite%lai(ipa) + cliq *   &
                cpatch%veg_water(ico))
        else
           csite%mean_latflux(ipa) = csite%mean_latflux(ipa) +   &
                vp_gradient * et_conductance * alvl -   &
                cice * cpatch%veg_temp(ico)* explicit_deriv_portion(ind2)
           explicit_deriv_portion(ind1) = dQdt / (cpatch%hcapveg(ico) * cpatch%lai(ico) /   &
                csite%lai(ipa) + cice * cpatch%veg_water(ico))
        endif

        ! Load initial state of the cohorts:
        original_state(ind1) = cpatch%veg_temp(ico)
        original_state(ind2) = cpatch%veg_water(ico)

     endif
  enddo

  ! original state of canopy
  original_state(1) = csite%can_temp(ipa)
  original_state(2) = csite%can_shv(ipa)
  
  ! new state of canopy
  do ic=1,ndims
     explicit_new_state(ic) = original_state(ic) +   &
          explicit_deriv_portion(ic) * dt_leaf
  enddo

  ! Load the new state into the patch and cohort structures
  csite%can_temp(ipa) = explicit_new_state(1)
  csite%can_shv(ipa) = explicit_new_state(2)
  idim = 1
  ic = 0
  do ico = 1,cpatch%ncohorts
     if(cpatch%solvable(ico))then
        ic = ic + 1
        idim = idim + 2
        cpatch%veg_temp(ico) = explicit_new_state(idim)
        cpatch%veg_water(ico) = explicit_new_state(idim+1)
     endif
  enddo

  ed_transp(:) = ed_transp(:) * dt_leaf

  ! Require 0.22 LAI >= veg_water >= 0
  do ico = 1,cpatch%ncohorts
     if(cpatch%solvable(ico))then
        if(cpatch%veg_water(ico) < 0.0)then
           
           ! Take away from the canopy humidity...
           csite%can_shv(ipa) = csite%can_shv(ipa) + cpatch%veg_water(ico) / canwcap
           
           ! ... and set the veg_water to zero
           cpatch%veg_water(ico) = 0.0
           
        endif
     endif
  enddo
  csite%can_rhos(ipa)       = idealdenssh(prss,csite%can_temp(ipa),csite%can_shv(ipa))
  csite%can_enthalpy(ipa)   = tq2enthalpy(csite%can_temp(ipa),csite%can_shv(ipa))
  csite%avg_daily_temp(ipa) = csite%avg_daily_temp(ipa) + csite%can_temp(ipa)

  return
end subroutine canopy_explicit_driver

!======================================================================

subroutine ludcmp(a,n,np,indx,d)

  implicit none
  
  real, parameter :: tiny_offset=1.0e-20
  integer :: n
  integer :: np
  real :: d
  real, dimension(np,np) :: a
  integer, dimension(n) :: indx
  real, dimension(np) :: vv
  integer :: i
  integer :: j
  integer :: k
  integer :: imax
  real :: sum
  real :: aamax
  real :: dum

  d = 1.0
  do i=1,n
     aamax = 0.0
     do j=1,n
        if(abs(a(i,j)) > aamax)aamax = abs(a(i,j))
     enddo
     if(aamax == 0.0)then
        print*,'singular matrix in ludcmp'
        do j=1,n
           print*,i,j,abs(a(i,j)),aamax
        enddo
        stop
     endif
     vv(i) = 1.0 / aamax
  enddo

  do j=1,n
     if(j.gt.1)then
        do i=1,j-1
           sum = a(i,j)
           if(i > 1)then
              do k=1,i-1
                 sum = sum - a(i,k) * a(k,j)
              enddo
              a(i,j) = sum
           endif
        enddo
     endif
     aamax = 0.0
     do i=j,n
        sum = a(i,j)
        if (j > 1)then
           do k=1,j-1
              sum = sum - a(i,k) * a(k,j)
           enddo
           a(i,j) = sum
        endif
        dum = vv(i) * abs(sum)
        if(dum >= aamax)then
           imax = i
           aamax = dum
        endif
     enddo
     if(j /= imax)then
        do k=1,n
           dum = a(imax,k)
           a(imax,k) = a(j,k)
           a(j,k) = dum
        enddo
        d = -d
        vv(imax) = vv(j)
     endif
     indx(j) = imax
     if(j /= n)then
        if(a(j,j) == 0.0) a(j,j) = tiny_offset
        dum = 1.0 / a(j,j)
        do i=j+1,n
           a(i,j) = a(i,j)*dum
        enddo
     endif
  enddo
  if(a(n,n) == 0.0)a(n,n) = tiny_offset

  return
end subroutine ludcmp

!======================================================================

subroutine ludcmp_dble(a,n,np,indx,d)

  implicit none
  
  real(kind=8), parameter :: tiny_offset = 1.0d-20
  integer, intent(in) :: n
  integer, intent(in) :: np
  real, intent(out) :: d
  real, dimension(np,np), intent(inout) :: a
  real(kind=8), dimension(np,np) :: ad
  integer, dimension(n), intent(out) :: indx
  real(kind=8), dimension(np) :: vv
  integer :: i
  integer :: j
  integer :: k
  integer :: imax
  real(kind=8) :: sum
  real(kind=8) :: aamax
  real(kind=8) :: dum

  ad = dble(a)

  d = 1.0

  do i = 1, n
     aamax = 0.0d0
     do j = 1, n
        if(abs(ad(i,j)) > aamax)aamax = abs(ad(i,j))
     enddo
     if(aamax == 0.0d0)then
        print*,'singular matrix in ludcmp'
        do j=1,n
           print*,i,j,abs(a(i,j)),aamax
        enddo
        stop
     endif
     vv(i) = 1.0d0 / aamax
  enddo

  do j = 1, n
     if(j.gt.1)then
        do i = 1, j - 1
           sum = ad(i,j)
           if(i > 1)then
              do k = 1, i - 1
                 sum = sum - ad(i,k) * ad(k,j)
              enddo
              ad(i,j) = sum
           endif
        enddo
     endif
     aamax = 0.0d0
     do i=j,n
        sum = ad(i,j)
        if (j > 1)then
           do k = 1, j - 1
              sum = sum - ad(i,k) * ad(k,j)
           enddo
           ad(i,j) = sum
        endif
        dum = vv(i) * abs(sum)
        if(dum >= aamax)then
           imax = i
           aamax = dum
        endif
     enddo
     if(j /= imax)then
        do k = 1, n
           dum = ad(imax,k)
           ad(imax,k) = ad(j,k)
           ad(j,k) = dum
        enddo
        d = -d
        vv(imax) = vv(j)
     endif
     indx(j) = imax
     if(j /= n)then
        if(ad(j,j) == 0.0d0) ad(j,j) = tiny_offset
        dum = 1.0d0 / ad(j,j)
        do i = j + 1, n
           ad(i,j) = ad(i,j) * dum
        enddo
     endif
  enddo
  if(ad(n,n) == 0.0d0)ad(n,n) = tiny_offset

  a = real(ad)

  return
end subroutine ludcmp_dble

!==================================================================

subroutine lubksb(a,n,np,indx,b)
  implicit none
  integer :: n
  integer :: np
  integer, dimension(n) :: indx
  real, dimension(n) :: b
  real, dimension(np,np) :: a
  integer :: ii
  integer :: i
  integer :: ll
  real :: sum
  integer :: j

  ii = 0
  do i=1,n
     ll = indx(i)
     sum = b(ll)
     b(ll) = b(i)
     if(ii /= 0)then
        do j=ii,i-1
           sum = sum - a(i,j) * b(j)
        enddo
     elseif(sum /= 0.0)then
        ii = i
     endif
     b(i) = sum
  enddo
  do i=n,1,-1
     sum = b(i)
     if ( i < n )then
        do j=i+1,n
           sum = sum - a(i,j) * b(j)
        enddo
     endif
     b(i) = sum / a(i,i)
  enddo

  return
end subroutine lubksb

!==================================================================

subroutine lubksb_dble(a,n,np,indx,b)
  implicit none
  integer, intent(in) :: n
  integer, intent(in) :: np
  integer, dimension(n), intent(in) :: indx
  real, dimension(n), intent(inout) :: b
  real(kind=8), dimension(n) :: bd
  real, dimension(np,np), intent(in) :: a
  real(kind=8), dimension(np,np) :: ad
  integer :: ii
  integer :: i
  integer :: ll
  real(kind=8) :: sum
  integer :: j

  ad = dble(a)
  bd = dble(b)

  ii = 0
  do i=1,n
     ll = indx(i)
     sum = bd(ll)
     bd(ll) = bd(i)
     if(ii /= 0)then
        do j=ii,i-1
           sum = sum - ad(i,j) * bd(j)
        enddo
     elseif(sum /= 0.0d0)then
        ii = i
     endif
     bd(i) = sum
  enddo
  do i=n,1,-1
     sum = bd(i)
     if ( i < n )then
        do j=i+1,n
           sum = sum - ad(i,j) * bd(j)
        enddo
     endif
     bd(i) = sum / ad(i,i)
  enddo

  b = real(bd)

  return
end subroutine lubksb_dble

