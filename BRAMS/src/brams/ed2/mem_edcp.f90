module mem_edcp

  type cpl_vars
     integer,  pointer :: xatm(:,:)
     integer,  pointer :: yatm(:,:)
     real,     pointer :: glatll(:,:)
     real,     pointer :: glonll(:,:)

  end type cpl_vars

  type(cpl_vars),pointer :: cpl_g(:)

  type ed_flux
     real, pointer, dimension (:,:) :: ustar
     real, pointer, dimension (:,:) :: tstar
     real, pointer, dimension (:,:) :: rstar
     real, pointer, dimension (:,:) :: sflux_u
     real, pointer, dimension (:,:) :: sflux_v
     real, pointer, dimension (:,:) :: sflux_t
     real, pointer, dimension (:,:) :: sflux_r
     real, pointer, dimension (:,:) :: sflux_w
     real, pointer, dimension (:,:) :: rlongup
     real, pointer, dimension (:,:) :: albedt
  end type ed_flux
  
  type(ed_flux),pointer, dimension(:) :: ed_fluxp_g,ed_fluxf_g

  type water_grid
     real, pointer, dimension(:,:) :: ustar
     real, pointer, dimension(:,:) :: tstar
     real, pointer, dimension(:,:) :: rstar
     real, pointer, dimension(:,:) :: sflux_u
     real, pointer, dimension(:,:) :: sflux_v
     real, pointer, dimension(:,:) :: sflux_t
     real, pointer, dimension(:,:) :: sflux_r
     real, pointer, dimension(:,:) :: sflux_w
     real, pointer, dimension(:,:) :: rlongup
     real, pointer, dimension(:,:) :: albedt
  end type water_grid

  type(water_grid),pointer,dimension(:) :: wgrid_g

  !    These variables are needed because the microphysics, cumulus, and ED
  ! have different update frequency, so we store the total accumulated precipitation
  ! at the previous ED call.
  !------------------------------------------------------------------------------------
 
  type ed_precip
     real, pointer, dimension(:,:) :: prev_aconpr
     real, pointer, dimension(:,:) :: prev_abulkpr
  end type ed_precip
  type (ed_precip), pointer, dimension(:) :: ed_precip_g


  real(kind=8) :: edtime1,edtime2

contains
  
  subroutine alloc_edflux(ed,n2,n3)
    implicit none
    type(ed_flux) :: ed
    integer, intent(in) :: n2,n3
    
    call nullify_edflux(ed)
    
    allocate(ed%ustar   (n2,n3)   )
    allocate(ed%tstar   (n2,n3)   )
    allocate(ed%rstar   (n2,n3)   )
    allocate(ed%sflux_u   (n2,n3)   )
    allocate(ed%sflux_v   (n2,n3)   )
    allocate(ed%sflux_r   (n2,n3)   )
    allocate(ed%sflux_t   (n2,n3)   )
    allocate(ed%sflux_w   (n2,n3)   )
    allocate(ed%rlongup   (n2,n3)   )
    allocate(ed%albedt    (n2,n3)   )
    
    return
  end subroutine alloc_edflux
  
  !---------------------------------------------------------!
  
  subroutine nullify_edflux(ed)
    implicit none
    type(ed_flux) :: ed
    
    if (associated(ed%ustar     ))  nullify(ed%ustar   )
    if (associated(ed%tstar     ))  nullify(ed%tstar   )
    if (associated(ed%rstar     ))  nullify(ed%rstar   )
    if (associated(ed%sflux_u   ))  nullify(ed%sflux_u   )
    if (associated(ed%sflux_v   ))  nullify(ed%sflux_v   )
    if (associated(ed%sflux_r   ))  nullify(ed%sflux_r   )
    if (associated(ed%sflux_t   ))  nullify(ed%sflux_t   )
    if (associated(ed%sflux_w   ))  nullify(ed%sflux_w   )
    if (associated(ed%albedt    ))  nullify(ed%albedt   )
    if (associated(ed%rlongup   ))  nullify(ed%rlongup )
    return
  end subroutine nullify_edflux
  
  !---------------------------------------------------------!

  subroutine zero_edflux(ed)
    implicit none
    type(ed_flux) :: ed

    if (associated(ed%ustar     ))  ed%ustar   = 0.0
    if (associated(ed%tstar     ))  ed%tstar   = 0.0
    if (associated(ed%rstar     ))  ed%rstar   = 0.0
    if (associated(ed%sflux_u   ))  ed%sflux_u = 0.0
    if (associated(ed%sflux_v   ))  ed%sflux_v = 0.0
    if (associated(ed%sflux_r   ))  ed%sflux_r = 0.0
    if (associated(ed%sflux_t   ))  ed%sflux_t = 0.0
    if (associated(ed%sflux_w   ))  ed%sflux_w = 0.0
    if (associated(ed%albedt    ))  ed%albedt  = 0.0
    if (associated(ed%rlongup   ))  ed%rlongup = 0.0
    return
  end subroutine zero_edflux

  !---------------------------------------------------------!

  subroutine dealloc_edflux(ed)
    implicit none
    type(ed_flux) :: ed
    
    if (associated(ed%ustar  ))deallocate(ed%ustar  )
    if (associated(ed%tstar  ))deallocate(ed%tstar  )
    if (associated(ed%rstar  ))deallocate(ed%rstar  )
    if (associated(ed%rlongup  ))deallocate(ed%rlongup  )
    if (associated(ed%albedt   ))deallocate(ed%albedt   )
    if (associated(ed%sflux_u   ))  deallocate(ed%sflux_u   )
    if (associated(ed%sflux_v   ))  deallocate(ed%sflux_v   )
    if (associated(ed%sflux_r   ))  deallocate(ed%sflux_r   )
    if (associated(ed%sflux_t   ))  deallocate(ed%sflux_t   )
    if (associated(ed%sflux_w   ))  deallocate(ed%sflux_w   )
    return
  end subroutine dealloc_edflux

  ! =====================================================

  subroutine alloc_wgrid(wg,n2,n3)
    implicit none
    type(water_grid) :: wg
    integer, intent(in) :: n2,n3
    
    call nullify_wgrid(wg)
    
    allocate(wg%ustar     (n2,n3)   )
    allocate(wg%tstar     (n2,n3)   )
    allocate(wg%rstar     (n2,n3)   )
    allocate(wg%rlongup   (n2,n3)   )
    allocate(wg%albedt    (n2,n3)   )
    allocate(wg%sflux_u   (n2,n3)   )
    allocate(wg%sflux_v   (n2,n3)   )
    allocate(wg%sflux_r   (n2,n3)   )
    allocate(wg%sflux_t   (n2,n3)   )
    allocate(wg%sflux_w   (n2,n3)   )

    return
  end subroutine alloc_wgrid
  
  !---------------------------------------------------------!
  
  subroutine nullify_wgrid(wg)
    implicit none
    type(water_grid) :: wg
    
    if (associated(wg%ustar     ))  nullify(wg%ustar   )
    if (associated(wg%tstar     ))  nullify(wg%tstar   )
    if (associated(wg%rstar     ))  nullify(wg%rstar   )
    if (associated(wg%albedt    ))  nullify(wg%albedt   )
    if (associated(wg%rlongup   ))  nullify(wg%rlongup ) 
    if (associated(wg%sflux_u   ))  nullify(wg%sflux_u   )
    if (associated(wg%sflux_v   ))  nullify(wg%sflux_v   )
    if (associated(wg%sflux_r   ))  nullify(wg%sflux_r   )
    if (associated(wg%sflux_t   ))  nullify(wg%sflux_t   )
    if (associated(wg%sflux_w   ))  nullify(wg%sflux_w   )
    return
  end subroutine nullify_wgrid
  
  !---------------------------------------------------------!

  subroutine zero_wgrid(wg)
    implicit none
    type(water_grid) :: wg
    
    if (associated(wg%ustar     ))  wg%ustar   = 0.0
    if (associated(wg%tstar     ))  wg%tstar   = 0.0
    if (associated(wg%rstar     ))  wg%rstar   = 0.0
    if (associated(wg%albedt    ))  wg%albedt  = 0.0
    if (associated(wg%rlongup   ))  wg%rlongup = 0.0 
    if (associated(wg%sflux_u   ))  wg%sflux_u = 0.0
    if (associated(wg%sflux_v   ))  wg%sflux_v = 0.0
    if (associated(wg%sflux_r   ))  wg%sflux_r = 0.0
    if (associated(wg%sflux_t   ))  wg%sflux_t = 0.0
    if (associated(wg%sflux_w   ))  wg%sflux_w = 0.0
    return
  end subroutine zero_wgrid
  
  !---------------------------------------------------------!
  
  subroutine dealloc_wgrid(wg)
    implicit none
    
    type(water_grid) :: wg
    
    if (associated(wg%ustar  ))deallocate(wg%ustar  )
    if (associated(wg%tstar  ))deallocate(wg%tstar  )
    if (associated(wg%rstar  ))deallocate(wg%rstar  )
    if (associated(wg%rlongup  ))deallocate(wg%rlongup  )
    if (associated(wg%albedt   ))deallocate(wg%albedt   )
    if (associated(wg%sflux_u   ))  deallocate(wg%sflux_u   )
    if (associated(wg%sflux_v   ))  deallocate(wg%sflux_v   )
    if (associated(wg%sflux_r   ))  deallocate(wg%sflux_r   )
    if (associated(wg%sflux_t   ))  deallocate(wg%sflux_t   )
    if (associated(wg%sflux_w   ))  deallocate(wg%sflux_w   )
    return
  end subroutine dealloc_wgrid
  
  !---------------------------------------------------------!
  
  subroutine alloc_edprecip(edp,n2,n3)
     implicit none
     type(ed_precip), intent(inout) :: edp
     integer, intent(in)            :: n2,n3

     call nullify_edprecip(edp)

     allocate(edp%prev_aconpr (n2,n3))
     allocate(edp%prev_abulkpr(n2,n3))

     return
  end subroutine alloc_edprecip
  
  !---------------------------------------------------------!
  
  subroutine nullify_edprecip(edp)
     type(ed_precip), intent(inout) :: edp
     
     if (associated(edp%prev_aconpr )) nullify(edp%prev_aconpr )
     if (associated(edp%prev_abulkpr)) nullify(edp%prev_abulkpr)

     return
  end subroutine  nullify_edprecip
  
  !---------------------------------------------------------!
  
  subroutine zero_edprecip(edp)
     type(ed_precip), intent(inout) :: edp
     
     if (associated(edp%prev_aconpr )) edp%prev_aconpr  = 0.0
     if (associated(edp%prev_abulkpr)) edp%prev_abulkpr = 0.0

     return
  end subroutine  zero_edprecip
  
  !---------------------------------------------------------!
  
  subroutine dealloc_ed_precip(edp)
     type(ed_precip), intent(inout) :: edp
     
     if (associated(edp%prev_aconpr )) deallocate(edp%prev_aconpr )
     if (associated(edp%prev_abulkpr)) deallocate(edp%prev_abulkpr)

     return
  end subroutine  dealloc_ed_precip


end module mem_edcp
