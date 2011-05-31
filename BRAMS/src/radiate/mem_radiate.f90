!====================================== Change Log ========================================!
! 5.0.0                                                                                    !
!                                                                                          !
!==========================================================================================!
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved                       !
!  Regional Atmospheric Modeling System - RAMS                                             !
!==========================================================================================!
!==========================================================================================!
module mem_radiate

   type radiate_vars
      !----- Variables to be dimensioned by (nzp,nxp,nyp). --------------------------------!
      real, dimension(:,:,:), pointer :: fthrd          ! Heating rate - total     [   K/s]
      real, dimension(:,:,:), pointer :: fthrd_lw       ! Heating rate - Long wave [   K/s]
      !----- Variables to be dimensioned by (    nxp,nyp). --------------------------------!
      real, dimension(  :,:), pointer :: cosz           ! Cosine of zenith angle   [   ---]
      real, dimension(  :,:), pointer :: albedt         ! Total albedo             [   ---]
      real, dimension(  :,:), pointer :: rshort         ! Dnward SW sfc. radiation [  W/m²]
      real, dimension(  :,:), pointer :: rshort_top     ! Dnward SW TOA  radiation [  W/m²]
      real, dimension(  :,:), pointer :: rshortup_top   ! Upward SW TOA  radiation [  W/m²]
      real, dimension(  :,:), pointer :: rshort_diffuse ! Diffuse dnward SW rad.   [  W/m²]
      real, dimension(  :,:), pointer :: rlong          ! Dnward LW sfc. radiation [  W/m²]
      real, dimension(  :,:), pointer :: rlongup        ! Upward LW sfc. radiation [  W/m²]
      real, dimension(  :,:), pointer :: rlongup_top    ! Upward LW TOA  radiation [  W/m²]
      !------------------------------------------------------------------------------------!
   end Type
   
   type (radiate_vars), dimension(:), allocatable :: radiate_g  ! Instantaneous variables
   type (radiate_vars), dimension(:), allocatable :: radiatem_g ! Averaged variables

   !---------------------------------------------------------------------------------------!
   !     Other variables.                                                                  !
   !---------------------------------------------------------------------------------------!
   integer            :: lonrad
   integer            :: ilwrtyp
   integer            :: iswrtyp
   integer            :: icumfdbk
   real               :: radfrq
   integer            :: ncall_i 
   real               :: prsnz
   real               :: prsnzp  
   integer            :: ncrad                 ! Number of clouds that affect radiation
   real   , parameter :: rad_cosz_min   = 0.0009
   real   , parameter :: rad_rshort_min = 0.50
   !----- These are used for adding extra levels at the top with the Mclatchy soundings. --!
   integer, parameter :: maxadd_rad     = 10   ! max allowed # of added rad levels
   integer            :: nadd_rad              ! actual # of added radiation levels
   real               :: zmrad                 ! Top of radiation grid
   real               :: solfac                ! solar-constant coefficient for variable 
                                               !      Earth-Sun dist
   integer            :: jday                  ! Day of year
   real(kind=8)       :: cdec                  ! cosine of declination
   real(kind=8)       :: sdec                  ! sine of declination
   real(kind=8)       :: declin                ! Declination
   real(kind=8)       :: sun_longitude         ! longitude where Sun is at zenith


   !=======================================================================================!
   !=======================================================================================!


   contains



   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine will allocate the pointers of the radiate structure, based on     !
   ! options if necessary.                                                                 !
   !---------------------------------------------------------------------------------------!
   subroutine alloc_radiate(radiate,n1,n2,n3,ng)

      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(radiate_vars), intent(inout) :: radiate
      integer           , intent(in)    :: n1
      integer           , intent(in)    :: n2
      integer           , intent(in)    :: n3
      integer           , intent(in)    :: ng
      !------------------------------------------------------------------------------------!

      if( ilwrtyp == 0 .and. iswrtyp == 0) return

      allocate (radiate%fthrd          (n1,n2,n3))
      allocate (radiate%rshort            (n2,n3))
      allocate (radiate%rlong             (n2,n3))
      allocate (radiate%rlongup           (n2,n3))
      allocate (radiate%albedt            (n2,n3))
      allocate (radiate%cosz              (n2,n3))
      allocate (radiate%fthrd_lw       (n1,n2,n3))
      allocate (radiate%rshort_top        (n2,n3))
      allocate (radiate%rshortup_top      (n2,n3))
      allocate (radiate%rlongup_top       (n2,n3))
      allocate (radiate%rshort_diffuse    (n2,n3))
                         
      return
   end subroutine alloc_radiate
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine will nullify all pointers for a safe allocation.                   !
   !---------------------------------------------------------------------------------------!
   subroutine nullify_radiate(radiate)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type (radiate_vars), intent(inout) :: radiate
      !------------------------------------------------------------------------------------!

      if (associated(radiate%fthrd         )) nullify (radiate%fthrd         )
      if (associated(radiate%rshort        )) nullify (radiate%rshort        )
      if (associated(radiate%rlong         )) nullify (radiate%rlong         )
      if (associated(radiate%rlongup       )) nullify (radiate%rlongup       )
      if (associated(radiate%albedt        )) nullify (radiate%albedt        )
      if (associated(radiate%cosz          )) nullify (radiate%cosz          )
      if (associated(radiate%fthrd_lw      )) nullify (radiate%fthrd_lw      )
      if (associated(radiate%rshort_top    )) nullify (radiate%rshort_top    )
      if (associated(radiate%rshortup_top  )) nullify (radiate%rshortup_top  )
      if (associated(radiate%rlongup_top   )) nullify (radiate%rlongup_top   )
      if (associated(radiate%rshort_diffuse)) nullify (radiate%rshort_diffuse)

      return
   end subroutine nullify_radiate
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine will assign zeroes to all variables.  This is to avoid some bogus  !
   ! values to variables that are never used.                                              !
   !---------------------------------------------------------------------------------------!
   subroutine zero_radiate(radiate)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type (radiate_vars), intent(inout) :: radiate
      !------------------------------------------------------------------------------------!

      if(associated(radiate%fthrd         )) radiate%fthrd          =0.0
      if(associated(radiate%rshort        )) radiate%rshort         =0.0
      if(associated(radiate%rlong         )) radiate%rlong          =0.0
      if(associated(radiate%rlongup       )) radiate%rlongup        =0.0
      if(associated(radiate%albedt        )) radiate%albedt         =0.0
      if(associated(radiate%cosz          )) radiate%cosz           =0.0
      if(associated(radiate%fthrd_lw      )) radiate%fthrd_lw       =0.0
      if(associated(radiate%rshort_top    )) radiate%rshort_top     =0.0
      if(associated(radiate%rshortup_top  )) radiate%rshortup_top   =0.0
      if(associated(radiate%rlongup_top   )) radiate%rlongup_top    =0.0
      if(associated(radiate%rshort_diffuse)) radiate%rshort_diffuse =0.0

      return
   end subroutine zero_radiate
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine will deallocate all pointers before deallocating the structure.   !
   !---------------------------------------------------------------------------------------!
   subroutine dealloc_radiate(radiate)
      !----- Arguments. -------------------------------------------------------------------!
      type (radiate_vars), intent(inout) :: radiate
      !------------------------------------------------------------------------------------!

      if (associated(radiate%fthrd         )) deallocate (radiate%fthrd         )
      if (associated(radiate%rshort        )) deallocate (radiate%rshort        )
      if (associated(radiate%rlong         )) deallocate (radiate%rlong         )
      if (associated(radiate%rlongup       )) deallocate (radiate%rlongup       )
      if (associated(radiate%albedt        )) deallocate (radiate%albedt        )
      if (associated(radiate%cosz          )) deallocate (radiate%cosz          )
      if (associated(radiate%fthrd_lw      )) deallocate (radiate%fthrd_lw      )
      if (associated(radiate%rshort_top    )) deallocate (radiate%rshort_top    )
      if (associated(radiate%rshortup_top  )) deallocate (radiate%rshortup_top  )
      if (associated(radiate%rlongup_top   )) deallocate (radiate%rlongup_top   )
      if (associated(radiate%rshort_diffuse)) deallocate (radiate%rshort_diffuse)

      return
   end subroutine dealloc_radiate
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine will fill pointers to arrays into variable tables.                !
   !---------------------------------------------------------------------------------------!
   subroutine filltab_radiate(radiate,radiatem,imean,nmz,nmx,nmy,ng)
      use var_tables
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(radiate_vars), intent(inout) :: radiate,radiatem
      integer           , intent(in)    :: imean
      integer           , intent(in)    :: nmz
      integer           , intent(in)    :: nmx
      integer           , intent(in)    :: nmy
      integer           , intent(in)    :: ng
      !----- Local variables. -------------------------------------------------------------!
      integer                           :: npts
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     3-D variables, dimensioned by nmz*nmx*nmy.                                     !
      !------------------------------------------------------------------------------------!

      npts = nmz * nmx * nmy

      if (associated(radiate%fthrd))                                                       &
         call vtables2(radiate%fthrd,radiatem%fthrd,ng,npts,imean                          &
                      ,'FTHRD :3:hist:anal:mpti:mpt3')

      if (associated(radiate%fthrd_lw))                                                    &
         call vtables2(radiate%fthrd_lw,radiatem%fthrd_lw,ng,npts,imean                    &
                      ,'FTHRD_LW :3:hist:anal:mpti:mpt3')

      !------------------------------------------------------------------------------------!
      !     2-D variables, dimensioned by nx*ny.                                           !
      !------------------------------------------------------------------------------------!
      npts=nmx*nmy
      if (associated(radiate%rshort))                                                      &
         call vtables2(radiate%rshort,radiatem%rshort,ng,npts,imean                        &
                      ,'RSHORT :2:hist:anal:mpti:mpt3')

      if (associated(radiate%rshort_diffuse))                                              &
         call vtables2(radiate%rshort_diffuse,radiatem%rshort_diffuse,ng,npts,imean        &
                      ,'RSHORT_DIFFUSE :2:hist:anal:mpti:mpt3')

      if (associated(radiate%rshort_top))                                                  &
         call vtables2(radiate%rshort_top,radiatem%rshort_top,ng,npts,imean                &
                      ,'RSHORT_TOP :2:hist:anal:mpti:mpt3')

      if (associated(radiate%rshortup_top))                                                &
         call vtables2(radiate%rshortup_top,radiatem%rshortup_top,ng,npts,imean            &
                      ,'RSHORTUP_TOP :2:hist:anal:mpti:mpt3')

      if (associated(radiate%rlong))                                                       &
         call vtables2(radiate%rlong,radiatem%rlong,ng,npts,imean                          &
                      ,'RLONG :2:hist:anal:mpti:mpt3')

      if (associated(radiate%rlongup))                                                     &
         call vtables2(radiate%rlongup,radiatem%rlongup,ng,npts,imean                      &
                      ,'RLONGUP :2:hist:anal:mpti:mpt3')

      if (associated(radiate%rlongup_top))                                                 &
         call vtables2(radiate%rlongup_top,radiatem%rlongup_top,ng,npts,imean              &
                      ,'RLONGUP_TOP :2:hist:anal:mpti:mpt3')

      if (associated(radiate%albedt))                                                      &
         call vtables2(radiate%albedt,radiatem%albedt,ng,npts,imean                        &
                      ,'ALBEDT :2:hist:anal:mpti:mpt3')

      if (associated(radiate%cosz))                                                        &
         call vtables2(radiate%cosz,radiatem%cosz,ng,npts,imean                            &
                       ,'COSZ :2:hist:anal:mpti:mpt3')

      return
   end subroutine filltab_radiate
   !=======================================================================================!
   !=======================================================================================!
end module mem_radiate
!==========================================================================================!
!==========================================================================================!
