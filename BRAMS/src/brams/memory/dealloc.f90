!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################


subroutine dealloc_all()

  use mem_all

  use mem_opt    ! Needed for optimization - ALF

  use catt_start, only: &
       CATT                        ! intent(in)

  use mem_aerad, only: &
       nwave,          &         !INTENT(IN)
       final_definitions_aerad   !Subroutine

  use mem_globaer, only: &
       final_definitions_globaer !Subroutine

  use mem_globrad, only: &
       final_definitions_globrad !Subroutine

  use teb_spm_start, only: TEB_SPM ! INTENT(IN)

  use mem_teb, only:  &    !for urban parameterization
       teb_g, tebm_g, &
       dealloc_teb         ! Subroutine

  use mem_teb_common, only: &  !for TEB common use variables
       tebc_g, tebcm_g,     &
       dealloc_tebc            ! Subroutine

  use mem_gaspart, only:      & !for gas emission 
       gaspart_g, gaspartm_g, &
       dealloc_gaspart          ! Subroutine


  use mem_scratch_grell, only : dealloc_scratch_grell
  use mem_ensemble, only : ensemble_e,dealloc_ensemble
  use grell_coms, only : nclouds
  use mem_mass, only : mass_g, massm_g, dealloc_mass


  implicit none

  ! deallocate all model memory.  Used on dynamic balance

  integer :: ng,ne

  deallocate(num_var,vtab_r,scalar_tab,num_scalar)

  call dealloc_tend(naddsc)
  call dealloc_scratch()

  call dealloc_scratch_grell()
  if (allocated(ensemble_e)) then
     do ne=1,nclouds
        call dealloc_ensemble(ensemble_e(ne))
     end do
  end if
  deallocate(ensemble_e)

  call dealloc_opt_scratch() ! For optimization - ALF

  if (ilwrtyp==4 .or. iswrtyp==4) then ! For CARMA
     call final_definitions_aerad()
     call final_definitions_globrad()
     call final_definitions_globaer()
  endif

  do ng=1,ngrids
     call dealloc_basic(basic_g(ng)) 
     call dealloc_basic(basicm_g(ng))
     call dealloc_cuparm(cuparm_g(ng))  
     call dealloc_cuparm(cuparmm_g(ng))
     call dealloc_grid(grid_g(ng))    
     call dealloc_grid(gridm_g(ng))
     call dealloc_leaf(leaf_g(ng))    
     call dealloc_leaf(leafm_g(ng))
     call dealloc_micro(micro_g(ng)) 
     call dealloc_micro(microm_g(ng))
     call dealloc_radiate(radiate_g(ng)) 
     call dealloc_radiate(radiatem_g(ng))
     call dealloc_turb(turb_g(ng)) 
     call dealloc_turb(turbm_g(ng))
     call dealloc_varinit(varinit_g(ng)) 
     call dealloc_varinit(varinitm_g(ng))

     call dealloc_mass(mass_g(ng))
     call dealloc_mass(massm_g(ng))

     call dealloc_oda(oda_g(ng)) 
     call dealloc_oda(odam_g(ng))

     if (TEB_SPM==1) then
        if(allocated(tebc_g)) then
           call dealloc_tebc(tebc_g(ng))      !for teb common
           call dealloc_tebc(tebcm_g(ng))     !for teb common
        endif
        if(allocated(teb_g)) then
           call dealloc_teb(teb_g(ng))      !for teb
           call dealloc_teb(tebm_g(ng))     !for teb
        endif
        if(allocated(gaspart_g)) then
           call dealloc_gaspart(gaspart_g(ng))      !for gas/paticles
           call dealloc_gaspart(gaspartm_g(ng))     !for gas/paticles
        endif
     endif

  enddo
  if (allocated(basic_g  )) deallocate(basic_g   ,basicm_g   )
  if (allocated(cuparm_g )) deallocate(cuparm_g  ,cuparmm_g  )
  if (allocated(grid_g   )) deallocate(grid_g    ,gridm_g    )
  if (allocated(leaf_g   )) deallocate(leaf_g    ,leafm_g    )
  if (allocated(micro_g  )) deallocate(micro_g   ,microm_g   )
  if (allocated(radiate_g)) deallocate(radiate_g ,radiatem_g )
  if (allocated(turb_g   )) deallocate(turb_g    ,turbm_g    )
  if (allocated(varinit_g)) deallocate(varinit_g ,varinitm_g )
  if (allocated(oda_g    )) deallocate(oda_g     ,odam_g     )
  if (allocated(mass_g   )) deallocate(mass_g    ,massm_g    )

  if (TEB_SPM==1) then
     if(allocated(teb_g)) then
        deallocate(teb_g, tebm_g)         ! for urban parameterization
     endif
     if(allocated(tebc_g)) then
        deallocate(tebc_g, tebcm_g)         ! for urban parameterization
     endif
     if(allocated(gaspart_g)) then
        deallocate(gaspart_g, gaspartm_g)         ! for urban parameterization
     endif
  endif

  if(allocated(scalar_g)) then
     do ng=1,ngrids
        call dealloc_scalar(scalar_g(:,ng),naddsc)
     enddo
     deallocate(scalar_g,scalarm_g)
  endif

  return
end subroutine dealloc_all
