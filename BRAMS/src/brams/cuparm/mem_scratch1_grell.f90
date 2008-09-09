! Module necessary to Grell Cumulus param.
! Scratch variables - module 1

MODULE mem_scratch1_grell

  !TYPE scratch1_grell_vars

     !srf- feb-05-2002 : Variables for cumulus transport scheme
     !                   adapted in july-15-2002 for 5.x version
     !
     !ngrids_cp = numero de grades onde a parametrizacao de cumulus e' usada
     !
     implicit none
     integer :: iruncon = 0

    type scratch1_grell_vars
    !3d dependence + grid (mmzp,mmxp,mmyp)
        real, pointer, dimension(:,:,:) :: &
             thetasta,                     & ! Temporary theta, with shallow cumulus effect
             rvsta                         ! ! Temporary rv, with shallow cumulus effect

    !3d dependence + grid (mmxp,mmyp,maxiens)
        integer, pointer, dimension(:,:,:) ::  &
             ierr4d,                             & ! ierr4d=0 there is convection
             jmin4d,                             & ! etl
             kdet4d,                             & ! detrainemnt level
             k224d,                              & !
             kbcon4d,                            & ! cloud base
             ktop4d,                             & ! cloud top
             kstabi4d,                           & ! cloud base
             kstabm4d,                           & ! cloud top
             kpbl4d                                ! pbl height

        real, pointer, dimension(:,:,:) :: &
             xmb4d,                              & ! updraft mass flux
             edt4d                                 ! mflx_down/mflx_up

        !4d dependence + grid, then  (mmzp,mmxp,mmyp,maxiens)
        real, pointer, dimension(:,:,:,:) ::   &
             prup5d,                             &
            clwup5d,                             &
              tup5d,                             &
             enup5d,                             & ! entrain up   rate
             endn5d,                             & ! entrain down rate
             deup5d,                             & ! detrain up   rate
             dedn5d,                             & ! detrain down rate
              zup5d,                             & ! norm mass flux up
              zdn5d,                             & ! norm mass flux down
             p_lw5d,                             & ! prec/clw up for wet deposition
             zcup5d,                             & ! z level
             pcup5d                                !p level

  end type scratch1_grell_vars

  type(scratch1_grell_vars), allocatable, dimension(:) :: sc1_grell_g

contains
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
  subroutine alloc_scratch1_grell(sc1_grell,m1,m2,m3,maxiens)

     implicit none
     
     type(scratch1_grell_vars) :: sc1_grell
     integer, intent(in) :: m1,m2,m3,maxiens

     allocate (sc1_grell%thetasta (m1,m2,m3)        )
     allocate (sc1_grell%rvsta    (m1,m2,m3)        )

     allocate (sc1_grell%ierr4d      (m2,m3,maxiens))
     allocate (sc1_grell%jmin4d      (m2,m3,maxiens))
     allocate (sc1_grell%kdet4d      (m2,m3,maxiens))
     allocate (sc1_grell%k224d       (m2,m3,maxiens))
     allocate (sc1_grell%kbcon4d     (m2,m3,maxiens))
     allocate (sc1_grell%ktop4d      (m2,m3,maxiens))
     allocate (sc1_grell%kstabi4d    (m2,m3,maxiens))
     allocate (sc1_grell%kstabm4d    (m2,m3,maxiens))
     allocate (sc1_grell%kpbl4d      (m2,m3,maxiens))

     allocate (sc1_grell%xmb4d       (m2,m3,maxiens))
     allocate (sc1_grell%edt4d       (m2,m3,maxiens))

     allocate (sc1_grell%zcup5d   (m1,m2,m3,maxiens))
     allocate (sc1_grell%pcup5d   (m1,m2,m3,maxiens))
     allocate (sc1_grell%prup5d   (m1,m2,m3,maxiens))
     allocate (sc1_grell%clwup5d  (m1,m2,m3,maxiens))
     allocate (sc1_grell%tup5d    (m1,m2,m3,maxiens))
     allocate (sc1_grell%p_lw5d   (m1,m2,m3,maxiens))

     allocate (sc1_grell%enup5d   (m1,m2,m3,maxiens))
     allocate (sc1_grell%endn5d   (m1,m2,m3,maxiens))
     allocate (sc1_grell%deup5d   (m1,m2,m3,maxiens))
     allocate (sc1_grell%dedn5d   (m1,m2,m3,maxiens))
     allocate (sc1_grell%zup5d    (m1,m2,m3,maxiens))
     allocate (sc1_grell%zdn5d    (m1,m2,m3,maxiens))

  end subroutine alloc_scratch1_grell
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
  subroutine nullify_scratch1_grell(sc1_grell)

     implicit none
     type (scratch1_grell_vars) :: sc1_grell
    
     if(associated(sc1_grell%thetasta )) nullify(sc1_grell%thetasta )
     if(associated(sc1_grell%rvsta    )) nullify(sc1_grell%rvsta    )

     if(associated(sc1_grell%ierr4d   )) nullify(sc1_grell%ierr4d   )
     if(associated(sc1_grell%jmin4d   )) nullify(sc1_grell%jmin4d   )
     if(associated(sc1_grell%kdet4d   )) nullify(sc1_grell%kdet4d   )
     if(associated(sc1_grell%k224d    )) nullify(sc1_grell%k224d    )
     if(associated(sc1_grell%kbcon4d  )) nullify(sc1_grell%kbcon4d  )
     if(associated(sc1_grell%ktop4d   )) nullify(sc1_grell%ktop4d   )
     if(associated(sc1_grell%kstabi4d )) nullify(sc1_grell%kstabi4d )
     if(associated(sc1_grell%kstabm4d )) nullify(sc1_grell%kstabm4d )
     if(associated(sc1_grell%kpbl4d   )) nullify(sc1_grell%kpbl4d   )

     if(associated(sc1_grell%xmb4d    )) nullify(sc1_grell%xmb4d    )
     if(associated(sc1_grell%edt4d    )) nullify(sc1_grell%edt4d    )

     if(associated(sc1_grell%zcup5d   )) nullify(sc1_grell%zcup5d   )
     if(associated(sc1_grell%pcup5d   )) nullify(sc1_grell%pcup5d   )
     if(associated(sc1_grell%prup5d   )) nullify(sc1_grell%prup5d   )
     if(associated(sc1_grell%clwup5d  )) nullify(sc1_grell%clwup5d  )
     if(associated(sc1_grell%tup5d    )) nullify(sc1_grell%tup5d    )
     if(associated(sc1_grell%p_lw5d   )) nullify(sc1_grell%p_lw5d   )

     if(associated(sc1_grell%enup5d   )) nullify(sc1_grell%enup5d   )
     if(associated(sc1_grell%endn5d   )) nullify(sc1_grell%endn5d   )
     if(associated(sc1_grell%deup5d   )) nullify(sc1_grell%deup5d   )
     if(associated(sc1_grell%dedn5d   )) nullify(sc1_grell%dedn5d   )
     if(associated(sc1_grell%zup5d    )) nullify(sc1_grell%zup5d    )
     if(associated(sc1_grell%zdn5d    )) nullify(sc1_grell%zdn5d    )
     
     return
  end subroutine nullify_scratch1_grell
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
  subroutine dealloc_scratch1_grell(sc1_grell)

     implicit none
     type (scratch1_grell_vars) :: sc1_grell
    
     if(associated(sc1_grell%thetasta )) deallocate(sc1_grell%thetasta )
     if(associated(sc1_grell%rvsta    )) deallocate(sc1_grell%rvsta    )

     if(associated(sc1_grell%ierr4d   )) deallocate(sc1_grell%ierr4d   )
     if(associated(sc1_grell%jmin4d   )) deallocate(sc1_grell%jmin4d   )
     if(associated(sc1_grell%kdet4d   )) deallocate(sc1_grell%kdet4d   )
     if(associated(sc1_grell%k224d    )) deallocate(sc1_grell%k224d    )
     if(associated(sc1_grell%kbcon4d  )) deallocate(sc1_grell%kbcon4d  )
     if(associated(sc1_grell%ktop4d   )) deallocate(sc1_grell%ktop4d   )
     if(associated(sc1_grell%kstabi4d )) deallocate(sc1_grell%kstabi4d )
     if(associated(sc1_grell%kstabm4d )) deallocate(sc1_grell%kstabm4d )
     if(associated(sc1_grell%kpbl4d   )) deallocate(sc1_grell%kpbl4d   )

     if(associated(sc1_grell%xmb4d    )) deallocate(sc1_grell%xmb4d    )
     if(associated(sc1_grell%edt4d    )) deallocate(sc1_grell%edt4d    )

     if(associated(sc1_grell%zcup5d   )) deallocate(sc1_grell%zcup5d   )
     if(associated(sc1_grell%pcup5d   )) deallocate(sc1_grell%pcup5d   )
     if(associated(sc1_grell%prup5d   )) deallocate(sc1_grell%prup5d   )
     if(associated(sc1_grell%clwup5d  )) deallocate(sc1_grell%clwup5d  )
     if(associated(sc1_grell%tup5d    )) deallocate(sc1_grell%tup5d    )
     if(associated(sc1_grell%p_lw5d   )) deallocate(sc1_grell%p_lw5d   )

     if(associated(sc1_grell%enup5d   )) deallocate(sc1_grell%enup5d   )
     if(associated(sc1_grell%endn5d   )) deallocate(sc1_grell%endn5d   )
     if(associated(sc1_grell%deup5d   )) deallocate(sc1_grell%deup5d   )
     if(associated(sc1_grell%dedn5d   )) deallocate(sc1_grell%dedn5d   )
     if(associated(sc1_grell%zup5d    )) deallocate(sc1_grell%zup5d    )
     if(associated(sc1_grell%zdn5d    )) deallocate(sc1_grell%zdn5d    )
     
     return
  end subroutine dealloc_scratch1_grell
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
 subroutine zero_scratch1_grell(sc1_grell)

   implicit none
   type (scratch1_grell_vars) :: sc1_grell

     if(associated(sc1_grell%thetasta )) sc1_grell%thetasta = 0.
     if(associated(sc1_grell%rvsta    )) sc1_grell%rvsta    = 0.

     if(associated(sc1_grell%ierr4d   )) sc1_grell%ierr4d   = 0
     if(associated(sc1_grell%jmin4d   )) sc1_grell%jmin4d   = 0
     if(associated(sc1_grell%kdet4d   )) sc1_grell%kdet4d   = 0
     if(associated(sc1_grell%k224d    )) sc1_grell%k224d    = 0
     if(associated(sc1_grell%kbcon4d  )) sc1_grell%kbcon4d  = 0
     if(associated(sc1_grell%ktop4d   )) sc1_grell%ktop4d   = 0
     if(associated(sc1_grell%kstabi4d )) sc1_grell%kstabi4d = 0
     if(associated(sc1_grell%kstabm4d )) sc1_grell%kstabm4d = 0
     if(associated(sc1_grell%kpbl4d   )) sc1_grell%kpbl4d   = 0

     if(associated(sc1_grell%xmb4d    )) sc1_grell%xmb4d    = 0.
     if(associated(sc1_grell%edt4d    )) sc1_grell%edt4d    = 0.

     if(associated(sc1_grell%zcup5d   )) sc1_grell%zcup5d   = 0.
     if(associated(sc1_grell%pcup5d   )) sc1_grell%pcup5d   = 0.
     if(associated(sc1_grell%prup5d   )) sc1_grell%prup5d   = 0.
     if(associated(sc1_grell%clwup5d  )) sc1_grell%clwup5d  = 0.
     if(associated(sc1_grell%tup5d    )) sc1_grell%tup5d    = 0.
     if(associated(sc1_grell%p_lw5d   )) sc1_grell%p_lw5d   = 0.

     if(associated(sc1_grell%enup5d   )) sc1_grell%enup5d   = 0.
     if(associated(sc1_grell%endn5d   )) sc1_grell%endn5d   = 0.
     if(associated(sc1_grell%deup5d   )) sc1_grell%deup5d   = 0.
     if(associated(sc1_grell%dedn5d   )) sc1_grell%dedn5d   = 0.
     if(associated(sc1_grell%zup5d    )) sc1_grell%zup5d    = 0.
     if(associated(sc1_grell%zdn5d    )) sc1_grell%zdn5d    = 0.


  end subroutine zero_scratch1_grell

end module mem_scratch1_grell
