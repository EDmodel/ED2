module rk4_stepper_ar

contains 

  subroutine rkqs_ar(integration_buff, x, htry, hmin, epsil, hdid, hnext, csite, &
       ipa,isi,ipy,ifm,rhos, vels, atm_tmp, atm_shv, atm_co2, geoht, exner, pcpg, qpcpg, &
       prss, lsl)

    use ed_state_vars,only:sitetype,patchtype,rk4patchtype,integration_vars_ar &
                          ,edgrid_g
!    use lsm_integ_utils, only: print_patch, copy_rk4_patch, get_errmax, &
!         print_errmax

    implicit none

    integer, intent(in) :: lsl
    type(sitetype),target :: csite
    integer :: ipa,isi,ipy,ifm
    type(integration_vars_ar), target :: integration_buff

    real, parameter :: safety = 0.9
    real, parameter :: pgrow = -0.2
    real, parameter :: pshrnk = -0.25
    real, parameter :: errcon = 1.89e-4
    ! Changed eps by epsil and iflag by iflag1 because ifort didn't like these names
    ! when I ran with full interfacing...
    real :: h,htry,x,epsil,hdid,hnext,errmax,xnew,newh
    integer :: iflag1,iflag2
    logical :: minstep
    real :: hmin
    real, intent(in) :: rhos

    real, intent(in) :: vels
    real, intent(in) :: atm_tmp
    real, intent(in) :: atm_shv
    real, intent(in) :: atm_co2
    real, intent(in) :: geoht
    real, intent(in) :: exner
    real, intent(in) :: pcpg
    real, intent(in) :: qpcpg
    real, intent(in) :: prss

    h = htry
    hstep:   do

       !1) Try a step of varying size.
       
       call rkck_ar(integration_buff%y,   &
            integration_buff%dydx,   &
            integration_buff%ytemp,  &
            integration_buff%yerr,  &
            integration_buff%ak2,  &
            integration_buff%ak3,  &
            integration_buff%ak4,  &
            integration_buff%ak5,  &
            integration_buff%ak6,  &
            integration_buff%ak7,  &
            x, h, csite, ipa,isi,ipy, iflag1,  &
            rhos, vels, atm_tmp, atm_shv, atm_co2, geoht, exner,   &
            pcpg, qpcpg, prss, lsl)

       !2) Check to see how accurate the step was.  Errors
       !   were calculated by integrating the derivative
       !   of that last step.  

       if(iflag1.eq.1)then
          call get_errmax_ar(errmax, integration_buff%yerr,   &
               integration_buff%yscal, csite%patch(ipa), lsl,  &
               integration_buff%y,integration_buff%ytemp)
          errmax = errmax/epsil
       else

          call get_errmax_ar(errmax, integration_buff%yerr,   &
               integration_buff%yscal, csite%patch(ipa), lsl,  &
               integration_buff%y,integration_buff%ytemp)
          errmax = errmax/epsil
          if ( errmax < 1) then
!             print*,"INTEGRATOR DID NOT GIVE SANE RESULTS"
!             print*,"YET IT PASSED THE ERROR CRITERIA"
!             print*,"THIS SHOULD NOT BE. STOPPING"
          endif

          errmax = 10.0
       endif

       !3) If that error was large, then calculate a new
       !   step size to try.  There are two types of
       !   new tries.  If iflag1 is not 1, that means
       !   the step had finished prematurely, so we
       !   assign a standard large error (10.0). Otherwise
       !   a new step is calculated based on the size of that
       !   error.  Hopefully, those new steps should be less
       !   then the previous h.
       !   If the error was small, ie less then epsil, then
       !   we are done with this step, and we can move forward
       !   time: x = x + h

       if(errmax > 1.0)then

          newh = safety * h * errmax**pshrnk

	  minstep = (newh .eq. h)

	  if(minstep) then	
            print*,h,newh,safety,errmax,safety * errmax**pshrnk
	  endif

          if(newh < 0.1*h)then
             h = h * 0.1
          else
             h = newh
          endif
          xnew = x + h

          if(xnew == x .or. minstep)then

             print*,'stepsize underflow in rkqs'
             print*,'Longitude:',edgrid_g(ifm)%lon(ipy)
             print*,'Latitude:',edgrid_g(ifm)%lat(ipy)
             print*,'Polygon:',ipy
             print*,'patch age: ',csite%age(ipa)
             print*,'patch dist_type: ',csite%dist_type(ipa)
             print*,'iflag1: ',iflag1
             print*,'errmax: ',errmax
             if(iflag1 == 1)then
                print*,'Likely to be an errmax problem.'
             else
                print*,'Likely to be an iflag1 problem.'
	        print*,'turn on print_diag in lsm_sanity_check'
		print*,'for additional information'
             endif

             if(iflag1 == 1)then
                call print_errmax_ar(errmax, integration_buff%yerr,  &
                     integration_buff%yscal, csite%patch(ipa), lsl,  &
                     integration_buff%y,integration_buff%ytemp,epsil)
                print*,'errmax',errmax/epsil,'raw',errmax,'epsilon',epsil
             else
                call print_sanity_check_ar(integration_buff%y,  &
                     iflag2, csite,ipa, lsl)
             endif
             call print_patch_ar(integration_buff%y, csite,ipa, lsl)
             stop                 ! TIME TO DEBUG
          endif

       else
          
          if(errmax > errcon)then
             hnext = safety * h * errmax**pgrow
          else
             hnext = 5.0 * h
          endif

!          if (hnext/h > 1.1) print*,hnext/h

          
          hnext = max(2*hmin,hnext)

          x = x + h
          hdid = h
          call copy_rk4_patch_ar(integration_buff%ytemp,  &
               integration_buff%y, csite%patch(ipa), lsl)
          exit hstep
       endif
    enddo hstep

    return
  end subroutine rkqs_ar
  
  !==================================================================
  subroutine rkck_ar(y, dydx, yout, yerr, ak2, ak3, ak4, ak5, ak6, ak7,  &
       x, h, csite, ipa,isi,ipy, iflag1,   &
       rhos, vels, atm_tmp, atm_shv, atm_co2, geoht, exner, pcpg, qpcpg, &
       prss, lsl)
    
    use ed_state_vars, only: sitetype,patchtype,rk4patchtype,integration_vars_ar
!    use lsm_integ_utils, only: stabilize_snow_layers, inc_rk4_patch, copy_rk4_patch

    implicit none

    integer, intent(in) :: lsl
    type(rk4patchtype), target :: y
    type(rk4patchtype), target :: dydx
    type(rk4patchtype), target :: yout
    type(rk4patchtype), target :: yerr
    type(rk4patchtype), target :: ak1
    type(rk4patchtype), target :: ak2
    type(rk4patchtype), target :: ak3
    type(rk4patchtype), target :: ak4
    type(rk4patchtype), target :: ak5
    type(rk4patchtype), target :: ak6
    type(rk4patchtype), target :: ak7
    real, intent(in) :: x
    real, intent(in) :: h
    type(sitetype),target :: csite
    type(patchtype),pointer :: cpatch
    integer :: ipa,isi,ipy
    integer, intent(out) :: iflag1
    real, intent(in) :: rhos
    real, intent(in) :: vels
    real, intent(in) :: atm_tmp
    real, intent(in) :: atm_shv
    real, intent(in) :: atm_co2
    real, intent(in) :: geoht
    real, intent(in) :: exner
    real, intent(in) :: pcpg
    real, intent(in) :: qpcpg
    real, intent(in) :: prss

    real, parameter :: a2=0.2,a3=0.3,a4=0.6,a5=1.0,a6=0.875,b21=0.2,  &
         b31=3.0/40.0,b32=9.0/40.0,b41=0.3,b42=-0.9,b43=1.2,  &
         b51=-11.0/54.0,b52=2.5,b53=-70.0/27.0,b54=35.0/27.0, &
         b61=1631.0/55296.0,b62=175.0/512.0,b63=575.0/13824.0, &
         b64=44275.0/110592.0,b65=253.0/4096.0,c1=37.0/378.0, &
         c3=250.0/621.0,c4=125.0/594.0,c6=512.0/1771.0, &
         dc5=-277.0/14336.0
    real, parameter :: dc1=c1-2825.0/27648.0,dc3=c3-18575.0/48384.0, &
         dc4=c4-13525.0/55296.0,dc6=c6-0.25


#if USE_INTERF
    interface
       subroutine leaf_derivs_ar(initp, dydx, csite, ipa,isi,ipy, rhos, prss, pcpg, qpcpg,   &
            atm_tmp, exner, geoht, vels, atm_shv, atm_co2, lsl)
         
         use ed_state_vars,only:sitetype,rk4patchtype,patchtype
         implicit none
         integer, intent(in) :: lsl
         real, intent(in) :: rhos
         type (rk4patchtype) ,target :: initp
         type (rk4patchtype) ,target :: dydx
         type (sitetype) ,target :: csite
         type (patchtype),pointer :: cpatch
         integer :: ipa,isi,ipy
         real, intent(in) :: prss
         real, intent(in) :: pcpg
         real, intent(in) :: qpcpg
         real, intent(in) :: atm_tmp
         real, intent(in) :: exner
         real, intent(in) :: geoht
         real, intent(in) :: vels
         real, intent(in) :: atm_shv
         real, intent(in) :: atm_co2
       end subroutine leaf_derivs_ar
       
    end interface
#endif

    iflag1 = 1

    cpatch => csite%patch(ipa)

    call copy_rk4_patch_ar(y, ak7, cpatch, lsl)
    call inc_rk4_patch_ar(ak7, dydx, b21*h, cpatch, lsl)
    call stabilize_snow_layers_ar(ak7, csite, ipa, b21*h, lsl)
    call lsm_sanity_check_ar(ak7, iflag1, csite, ipa, lsl ,dydx,h )

    if(iflag1 /= 1)return

    call leaf_derivs_ar(ak7, ak2, csite, ipa,isi,ipy, rhos, prss, pcpg, qpcpg, atm_tmp,   &
         exner, geoht, vels, atm_shv, atm_co2, lsl)
    call copy_rk4_patch_ar(y, ak7, cpatch, lsl)
    call inc_rk4_patch_ar(ak7, dydx, b31*h, cpatch, lsl)
    call inc_rk4_patch_ar(ak7, ak2, b32*h, cpatch, lsl)
    call stabilize_snow_layers_ar(ak7, csite,ipa,(b31+b32)*h, lsl)
    call lsm_sanity_check_ar(ak7, iflag1,csite,ipa, lsl,dydx,h )

    if(iflag1 /= 1)return

    call leaf_derivs_ar(ak7, ak3, csite,ipa,isi,ipy, rhos, prss, pcpg, qpcpg, atm_tmp,   &
         exner, geoht, vels, atm_shv, atm_co2, lsl)
    call copy_rk4_patch_ar(y, ak7, cpatch, lsl)
    call inc_rk4_patch_ar(ak7, dydx, b41*h, cpatch, lsl)
    call inc_rk4_patch_ar(ak7, ak2, b42*h, cpatch, lsl)
    call inc_rk4_patch_ar(ak7, ak3, b43*h, cpatch, lsl)
    call stabilize_snow_layers_ar(ak7, csite,ipa, (b41+b42+b43)*h, lsl)
    call lsm_sanity_check_ar(ak7, iflag1, csite,ipa, lsl,dydx,h )

    if(iflag1 /= 1)return

    call leaf_derivs_ar(ak7, ak4, csite, ipa,isi,ipy, rhos, prss, pcpg, qpcpg, atm_tmp,   &
         exner, geoht, vels, atm_shv, atm_co2, lsl)
    call copy_rk4_patch_ar(y, ak7, cpatch, lsl)
    call inc_rk4_patch_ar(ak7, dydx, b51*h, cpatch, lsl)
    call inc_rk4_patch_ar(ak7, ak2, b52*h, cpatch, lsl)
    call inc_rk4_patch_ar(ak7, ak3, b53*h, cpatch, lsl)
    call inc_rk4_patch_ar(ak7, ak4, b54*h, cpatch, lsl)
    call stabilize_snow_layers_ar(ak7, csite, ipa, (b51+b52+b53+b54)*h, lsl)
    call lsm_sanity_check_ar(ak7, iflag1, csite, ipa, lsl,dydx,h )

    if(iflag1 /= 1)return

    call leaf_derivs_ar(ak7, ak5, csite, ipa,isi,ipy, rhos, prss, pcpg, qpcpg, atm_tmp,   &
         exner, geoht, vels, atm_shv, atm_co2, lsl)
    call copy_rk4_patch_ar(y, ak7, cpatch, lsl)
    call inc_rk4_patch_ar(ak7, dydx, b61*h, cpatch, lsl)
    call inc_rk4_patch_ar(ak7, ak2, b62*h, cpatch, lsl)
    call inc_rk4_patch_ar(ak7, ak3, b63*h, cpatch, lsl)
    call inc_rk4_patch_ar(ak7, ak4, b64*h, cpatch, lsl)
    call inc_rk4_patch_ar(ak7, ak5, b65*h, cpatch, lsl)
    call stabilize_snow_layers_ar(ak7, csite,ipa, (b61+b62+b63+b64+b65)*h, lsl)
    call lsm_sanity_check_ar(ak7, iflag1, csite,ipa, lsl,dydx,h )
    
    if(iflag1 /= 1)return

    call leaf_derivs_ar(ak7, ak6, csite,ipa,isi,ipy, rhos, prss, pcpg, qpcpg, atm_tmp,   &
         exner, geoht, vels, atm_shv, atm_co2, lsl)
    call copy_rk4_patch_ar(y, yout, cpatch, lsl)
    call inc_rk4_patch_ar(yout, dydx, c1*h, cpatch, lsl)
    call inc_rk4_patch_ar(yout, ak3, c3*h, cpatch, lsl)
    call inc_rk4_patch_ar(yout, ak4, c4*h, cpatch, lsl)
    call inc_rk4_patch_ar(yout, ak6, c6*h, cpatch, lsl)
    call stabilize_snow_layers_ar(yout, csite,ipa, h, lsl)
    call lsm_sanity_check_ar(yout, iflag1, csite,ipa, lsl,dydx,h )

    if(iflag1 /= 1)return

    call copy_rk4_patch_ar(dydx, yerr, cpatch, lsl)
    call inc_rk4_patch_ar(yerr, dydx, dc1*h-1.0, cpatch, lsl)
    call inc_rk4_patch_ar(yerr, ak3, dc3*h, cpatch, lsl)
    call inc_rk4_patch_ar(yerr, ak4, dc4*h, cpatch, lsl)
    call inc_rk4_patch_ar(yerr, ak5, dc5*h, cpatch, lsl)
    call inc_rk4_patch_ar(yerr, ak6, dc6*h, cpatch, lsl )

    return
  end subroutine rkck_ar

  !=====================================================================
  subroutine lsm_sanity_check_ar(y, iflag1, csite,ipa, lsl,dydx,h )

    use ed_state_vars, only: sitetype,patchtype,rk4patchtype,integration_vars_ar
    use grid_coms, only: nzg
    use soil_coms, only: soil
    use canopy_radiation_coms, only: lai_min
    use consts_coms, only : t3ple
    use canopy_air_coms, only: hcapveg_ref,heathite_min
    use therm_lib, only: qwtk

    implicit none
    integer, intent(in) :: lsl
    type(sitetype), target :: csite
    type(patchtype),pointer :: cpatch
    type(rk4patchtype), target :: y,dydx
    integer iflag1,k
    real :: atm_tempk,h,hcapveg,veg_temp,fracliq
    integer :: ipa,ico
    integer, parameter :: print_diags=0

    if(y%soil_tempk(nzg) /= y%soil_tempk(nzg))then
       print*,'in the sanity check'
       call print_patch_ar(y, csite,ipa, lsl)
       stop
    endif

    iflag1 = 1

    do k = lsl, nzg
       if(y%soil_tempk(k) < (t3ple-0.01) .and. y%soil_fracliq(k).gt.0.001)then
          iflag1 = 0
          if(print_diags==1)print*,'too much liquid',iflag1,  &
               y%soil_tempk(k),y%soil_fracliq(k)    
          return
       endif
       if(y%soil_tempk(k) > (t3ple+0.01) .and. y%soil_fracliq(k).lt.0.999)then
          iflag1 = 0
          if(print_diags==1)print*,'too much ice',iflag1,  &
               y%soil_tempk(k),y%soil_fracliq(k)    
          return
       endif
       if(y%soil_fracliq(k).gt.1.0 .or. y%soil_fracliq(k).lt.0.0)then
          iflag1 = 0
          if(print_diags==1) print*,'bad fracliq',iflag1,y%soil_fracliq(k)    
          return
       endif
    enddo

    if(y%can_temp.gt.350.0)then
       iflag1 = 0
       if(print_diags==1)print*,'canopy tempk too high',y%can_temp,dydx%can_temp,h
       return
    endif

    if(y%can_temp.lt.200.0)then
       iflag1 = 0
       if(print_diags==1) print*,'canopy tempk too low',y%can_temp,dydx%can_temp,h
       return
    endif

    if(y%nlev_sfcwater.eq.1.and.y%sfcwater_tempk(1).lt.205.0)then
       iflag1 = 0
       if(print_diags==1) print*,'sfcwater_tempk too low',y%sfcwater_tempk(1)
       return
    endif

    if(y%can_shv > 1.0)then
       iflag1 = 0
       if(print_diags==1)print*,'canopy water vapor too high',y%can_shv,dydx%can_shv,h
       return
    endif

    if(y%can_shv <= 0.0)then
       iflag1 = 0
       if(print_diags==1)print*,'canopy water vapor too low',y%can_shv,dydx%can_shv,h
       return
    endif

      do k=lsl,nzg
         if(y%soil_water(k).lt.soil(csite%ntext_soil(k,ipa))%soilcp)then
            iflag1 = 0
            if(print_diags==1) print*,'soil water too low',k,  &
            y%soil_water(k),csite%ntext_soil(k,ipa)
         endif
       if(y%soil_water(k).gt.1.0)then
          iflag1 = 0
          if(print_diags==1) print*,'soil water too high',k,  &
          y%soil_water(k),csite%ntext_soil(k,ipa)
       endif
       if(y%soil_tempk(k) > 350.0)then
          iflag1 = 0
          if(print_diags==1) print*,'soil_tempk too high',k,y%soil_tempk(k)
       endif

    enddo
    
    if(y%nlev_sfcwater >= 1)then
       if(y%sfcwater_mass(y%nlev_sfcwater).lt.-1.0e-3)then
          iflag1 = 0
          if(print_diags==1) print*,'sfcwater_mass too low',  &
               y%nlev_sfcwater,y%sfcwater_mass(y%nlev_sfcwater)
          return
       endif
    endif

    ! Negative rk4 increment factors can make this value significantly
    ! negative, but it should not slow down the integration
    ! Changed from -1.0e-3 to -1.0e-1
    
    if(y%virtual_water < -1.0e-1)then
       iflag1 = 0
       if(print_diags==1)print*,'virtual water too low',y%virtual_water
       return
    endif

    cpatch => csite%patch(ipa)

    do ico = 1,cpatch%ncohorts
    
       if (cpatch%lai(ico) > lai_min) then
          hcapveg = hcapveg_ref * max(cpatch%hite(1),heathite_min) * cpatch%lai(ico) / csite%lai(ipa)
          call qwtk(y%veg_energy(ico),y%veg_water(ico),hcapveg,veg_temp,fracliq)

          if(veg_temp > 380.0)then
             iflag1 = 0
             if(print_diags==1) print*,'leaf temp too high',veg_temp,y%veg_energy(ico),  &
                  cpatch%lai(ico),cpatch%pft(ico),cpatch%veg_temp(ico)
             return
          end if
       end if
    end do

    
    return
  end subroutine lsm_sanity_check_ar

  !=====================================================================

  subroutine print_sanity_check_ar(y, iflag1, csite,ipa, lsl)

    use ed_state_vars,only: sitetype,patchtype,rk4patchtype
    use grid_coms, only: nzg
    use soil_coms, only: soil
    use canopy_radiation_coms, only: lai_min

    implicit none
    integer, intent(in) :: lsl
    type(sitetype),target :: csite
    type(patchtype),pointer :: cpatch
    type(rk4patchtype), target :: y
    integer :: ipa,ico
    integer iflag1,k
    real :: atm_tempk

    write(unit=*,fmt='(64a)') ('=',k=1,64)
    write(unit=*,fmt='(64a)') ('=',k=1,64)
    write(unit=*,fmt='(a,20x,a,20x,a)') '======','SANITY CHECK','======'
    write(unit=*,fmt='(64a)') ('=',k=1,64)

    write(unit=*,fmt='(a)') ' '
    write(unit=*,fmt='(64a)') ('-',k=1,64)
    write(unit=*,fmt='(a5,3(1x,a12))') 'LEVEL','  SOIL_TEMPK','SOIL_FRACLIQ','  SOIL_WATER'
    do k=lsl,nzg
       write(unit=*,fmt='(i5,3(1x,es12.5))') &
            k, y%soil_tempk(k), y%soil_fracliq(k), y%soil_water(k)
    end do
    write(unit=*,fmt='(64a)') ('-',k=1,64)

    write(unit=*,fmt='(a)') ' '
    write(unit=*,fmt='(64a)') ('-',k=1,64)
    write(unit=*,fmt='(a5,3(1x,a12))') 'LEVEL','  OLD_SOIL_T','OLD_SOIL_FLQ','OLD_SOIL_H2O'
    do k=lsl,nzg
       write(unit=*,fmt='(i5,3(1x,es12.5))') &
            k, csite%soil_tempk(k,ipa), csite%soil_fracliq(k,ipa), csite%soil_water(k,ipa)
    end do
    write(unit=*,fmt='(64a)') ('-',k=1,64)

    write(unit=*,fmt='(a)') ' '
    write(unit=*,fmt='(64a)') ('-',k=1,64)
    write (unit=*,fmt='(a,1x,es12.5)') ' CAN_TEMP=     ',y%can_temp
    write (unit=*,fmt='(a,1x,es12.5)') ' OLD_CAN_TEMP= ',csite%can_temp(ipa)
    write (unit=*,fmt='(a,1x,es12.5)') ' CAN_VAPOR=    ',y%can_shv
    write (unit=*,fmt='(a,1x,es12.5)') ' OLD_CAN_VAP=  ',csite%can_shv(ipa)
    write (unit=*,fmt='(a,1x,i12)')    ' #LEV_SFCH2O=  ',y%nlev_sfcwater
    write (unit=*,fmt='(a,1x,i12)')    ' OLD_#_SFCH2O= ',csite%nlev_sfcwater(ipa)
    if(y%nlev_sfcwater == 1) then
       write(unit=*,fmt='(a,1x,es12.5)') ,'SFCWATER_TEMPK=',y%sfcwater_tempk(1)
    end if
    write(unit=*,fmt='(64a)') ('-',k=1,64)

    write(unit=*,fmt='(a)') ' '
    write(unit=*,fmt='(64a)') ('-',k=1,64)
    cpatch => csite%patch(ipa)
    write (unit=*,fmt='(2(a5,1x),4(a12,1x))') &
       '  COH','  PFT','         LAI','  VEG_ENERGY','OLD_VEG_ENER','OLD_VEG_TEMP'
    do ico = 1,cpatch%ncohorts
       if(cpatch%lai(ico) > lai_min) then
          write(unit=*,fmt='(2(i5,1x),4(es12.5,1x))') &
             ico,cpatch%pft(ico),cpatch%lai(ico),y%veg_energy(ico),cpatch%veg_energy(ico)  &
                                ,cpatch%veg_temp(ico)
       end if
    end do
    write(unit=*,fmt='(64a)') ('-',k=1,64)

    write(unit=*,fmt='(a)') ' '
    write(unit=*,fmt='(64a)') ('-',k=1,64)
    write (unit=*,fmt='(2(a5,1x),3(a12,1x))') &
       '  COH','  PFT','         LAI','   VEG_WATER',' OLD_VEG_H2O'
    do ico = 1,cpatch%ncohorts
       if(cpatch%lai(ico) > lai_min) then
          write(unit=*,fmt='(2(i5,1x),3(es12.5,1x))') &
             ico,cpatch%pft(ico),cpatch%lai(ico),y%veg_water(ico),cpatch%veg_water(ico)
       end if
    end do
    write(unit=*,fmt='(64a)') ('-',k=1,64)
    write(unit=*,fmt='(a)') ' '
    
    write(unit=*,fmt='(64a)') ('=',k=1,64)
    write(unit=*,fmt='(64a)') ('=',k=1,64)
    write(unit=*,fmt='(a)') ' '

    return
  end subroutine print_sanity_check_ar

end module rk4_stepper_ar
