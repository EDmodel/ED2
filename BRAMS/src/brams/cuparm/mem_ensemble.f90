!==========================================================================================!
!    Module mem_ensemble - This module contains all Grell's scratch variables that do      !
! depend on ensemble dimensions. Since shallow and deep convection have different ensemble !
! sizes, two scratch arrays are allocated, one for each spectral size                      !
!==========================================================================================!
module mem_ensemble

   
   type ensemble_vars

      !------ 3D dependence (maxens_dyn,maxens_lsf,maxens_eff,maxens_cap) -----------------!
      real, pointer, dimension(:,:,:,:) :: &
            dnmf_ens,                      & ! Reference downdraft mass flux      [kg/m²/s]
            upmf_ens                       ! ! Reference updraft mass flux        [kg/m²/s]
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    These variables are associated with each member of the second kind of ensemble, !
      ! that varies the precipitation efficiency. As described by Grell and Dévényi, the   !
      ! precipitation efficiency has a huge impact on the final answer, so that is a good  !
      ! variable to test for sensibility.                                                  !
      !------------------------------------------------------------------------------------!
      !------ 2D dependence (maxens_eff,maxens_cap) ---------------------------------------!
      real, pointer, dimension(:,:) ::     &
            edt_eff                        ! ! Precipitation efficiency for each member

      !------ 3D dependence (mgmzp,maxens_eff,maxens_cap) ---------------------------------!
      real, pointer, dimension(:,:,:) ::   & ! These are changes per unit of mass.
            dellatheiv_eff                 & ! Change in ice-liquid potential temperature
           ,dellathil_eff                  & ! Change in ice-liquid potential temperature
           ,dellaqtot_eff                  & ! Change in total mixing ratio
           ,pw_eff                         ! ! Water that doesn't evaporate (aka rain).
      !------------------------------------------------------------------------------------!
   end type ensemble_vars

   type(ensemble_vars), allocatable, dimension(:) :: ensemble_e
   contains
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
   subroutine alloc_ensemble(ensemble,mgmzp,maxens_dyn,maxens_lsf,maxens_eff,maxens_cap)

      implicit none
      type(ensemble_vars) :: ensemble
      integer, intent(in) :: mgmzp,maxens_dyn,maxens_lsf,maxens_eff,maxens_cap

      allocate (ensemble%dnmf_ens            (maxens_dyn,maxens_lsf,maxens_eff,maxens_cap))
      allocate (ensemble%upmf_ens            (maxens_dyn,maxens_lsf,maxens_eff,maxens_cap))

      allocate (ensemble%edt_eff                                   (maxens_eff,maxens_cap))

      allocate (ensemble%dellathil_eff (mgmzp,                      maxens_eff,maxens_cap))
      allocate (ensemble%dellatheiv_eff(mgmzp,                      maxens_eff,maxens_cap))
      allocate (ensemble%dellaqtot_eff (mgmzp,                      maxens_eff,maxens_cap))

      allocate (ensemble%pw_eff        (mgmzp,                      maxens_eff,maxens_cap))

      return
   end subroutine alloc_ensemble
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
   subroutine nullify_ensemble(ensemble)

      implicit none
      type(ensemble_vars) :: ensemble

      if(associated(ensemble%dnmf_ens      ))  nullify(ensemble%dnmf_ens      )
      if(associated(ensemble%upmf_ens      ))  nullify(ensemble%upmf_ens      )

      if(associated(ensemble%edt_eff       ))  nullify(ensemble%edt_eff       )

      if(associated(ensemble%dellathil_eff ))  nullify(ensemble%dellathil_eff )
      if(associated(ensemble%dellatheiv_eff))  nullify(ensemble%dellatheiv_eff)
      if(associated(ensemble%dellaqtot_eff ))  nullify(ensemble%dellaqtot_eff )

      if(associated(ensemble%pw_eff        ))  nullify(ensemble%pw_eff        )

      return
   end subroutine nullify_ensemble
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
   subroutine dealloc_ensemble(ensemble)

      implicit none
      type(ensemble_vars) :: ensemble

      if(associated(ensemble%dnmf_ens      ))  deallocate(ensemble%dnmf_ens      )
      if(associated(ensemble%upmf_ens      ))  deallocate(ensemble%upmf_ens      )

      if(associated(ensemble%edt_eff       ))  deallocate(ensemble%edt_eff       )

      if(associated(ensemble%dellathil_eff ))  deallocate(ensemble%dellathil_eff )
      if(associated(ensemble%dellatheiv_eff))  deallocate(ensemble%dellatheiv_eff)
      if(associated(ensemble%dellaqtot_eff ))  deallocate(ensemble%dellaqtot_eff )

      if(associated(ensemble%pw_eff        ))  deallocate(ensemble%pw_eff        )

      return
   end subroutine dealloc_ensemble
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
   subroutine zero_ensemble(ensemble)

      implicit none
      type(ensemble_vars) :: ensemble

      if(associated(ensemble%dnmf_ens      ))  ensemble%dnmf_ens       = 0.
      if(associated(ensemble%upmf_ens      ))  ensemble%upmf_ens       = 0.

      if(associated(ensemble%edt_eff       ))  ensemble%edt_eff        = 0.

      if(associated(ensemble%dellathil_eff ))  ensemble%dellathil_eff  = 0.
      if(associated(ensemble%dellatheiv_eff))  ensemble%dellatheiv_eff = 0.
      if(associated(ensemble%dellaqtot_eff ))  ensemble%dellaqtot_eff  = 0.

      if(associated(ensemble%pw_eff        ))  ensemble%pw_eff         = 0.

      return
   end subroutine zero_ensemble
!==========================================================================================!
!==========================================================================================!
end module mem_ensemble
