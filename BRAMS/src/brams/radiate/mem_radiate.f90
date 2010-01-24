!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################
!:DOC%BEGIN
!:DOC%TITLE Modulo Mem_Radiate

Module mem_radiate

   Type radiate_vars
   
      ! Variables to be dimensioned by (nzp,nxp,nyp)
   real, pointer, dimension(:,:,:) :: &
                          fthrd,fthrd_lw !DMM/MLO
                          
      ! Variables to be dimensioned by (nxp,nyp)
   real, pointer, dimension(:,:) :: &
                          rshort,rlong,rlongup,albedt,cosz &
!MLO: Adding fixed Harrington, adapted from David Medvigy's fix on OLAM
                         ,rshort_top,rshortup_top,rlongup_top &
                         ,albedt_beam,albedt_diffuse,rlong_albedo &
                         ,rshort_diffuse
   End Type
   
   type (radiate_vars), allocatable :: radiate_g(:), radiatem_g(:)
   
   integer :: lonrad,ilwrtyp,iswrtyp,icumfdbk
   real    :: radfrq
   integer :: ncall_i !Indica primeira chamada
   real    :: prsnz,prsnzp !Calculadas na primeira chamada

!DMM These are used for adding extra levels at the top with the Mclatchy soundings
   integer, parameter :: maxadd_rad = 10 ! max allowed # of added rad levels
   integer            :: nadd_rad        ! actual # of added radiation levels
   real               :: zmrad = 30.e3   ! top of radiation grid
   real               :: solfac          ! solar-constant coefficient for variable 
                                         !      Earth-Sun dist
!MLO Adding some variables that used to be in the radcom common block
   integer            :: jday            ! Day of year
   real(kind=8)       :: cdec            ! cosine of declination
   real(kind=8)       :: sdec            ! sine of declination
   real(kind=8)       :: declin          ! Declination
   real(kind=8)       :: sun_longitude   ! longitude where Sun is at zenith

!MLO Adding a variable to be used in the cumulus-radiation feedback.
   integer            :: ncrad     ! Number of clouds that affect radiation
Contains

   subroutine alloc_radiate(radiate,n1,n2,n3,ng)

   implicit none
   type (radiate_vars) :: radiate
   integer, intent(in) :: n1,n2,n3,ng

! Allocate arrays based on options (if necessary)
      
      if( ilwrtyp+iswrtyp > 0)  then
                         allocate (radiate%fthrd(n1,n2,n3))
                         allocate (radiate%rshort(n2,n3))
                         allocate (radiate%rlong(n2,n3))
                         allocate (radiate%rlongup(n2,n3))
                         allocate (radiate%albedt(n2,n3))
                         allocate (radiate%cosz(n2,n3))
![Harrington fix
                         allocate (radiate%fthrd_lw       (n1,n2,n3))
                         allocate (radiate%rshort_top     (n2,n3)   )
                         allocate (radiate%rshortup_top   (n2,n3)   )
                         allocate (radiate%rlongup_top    (n2,n3)   )
                         allocate (radiate%albedt_beam    (n2,n3)   )
                         allocate (radiate%albedt_diffuse (n2,n3)   )
                         allocate (radiate%rlong_albedo   (n2,n3)   )
                         allocate (radiate%rshort_diffuse (n2,n3)   )
!]
      endif
                         
   return
   end subroutine alloc_radiate


   subroutine nullify_radiate(radiate)

   implicit none
   type (radiate_vars) :: radiate
   

   if (associated(radiate%fthrd))          nullify (radiate%fthrd)
   if (associated(radiate%rshort))         nullify (radiate%rshort)
   if (associated(radiate%rlong))          nullify (radiate%rlong)
   if (associated(radiate%rlongup))        nullify (radiate%rlongup)
   if (associated(radiate%albedt))         nullify (radiate%albedt)
   if (associated(radiate%cosz))           nullify (radiate%cosz)
!Harrington fix
   if (associated(radiate%fthrd_lw      )) nullify (radiate%fthrd_lw      )
   if (associated(radiate%rshort_top    )) nullify (radiate%rshort_top    )
   if (associated(radiate%rshortup_top  )) nullify (radiate%rshortup_top  )
   if (associated(radiate%rlongup_top   )) nullify (radiate%rlongup_top   )
   if (associated(radiate%albedt_beam   )) nullify (radiate%albedt_beam   )
   if (associated(radiate%albedt_diffuse)) nullify (radiate%albedt_diffuse)
   if (associated(radiate%rlong_albedo  )) nullify (radiate%rlong_albedo  )
   if (associated(radiate%rshort_diffuse)) nullify (radiate%rshort_diffuse)

   return
   end subroutine nullify_radiate

   subroutine zero_radiate(radiate)

   implicit none
   type (radiate_vars) :: radiate
   

   if (associated(radiate%fthrd))          radiate%fthrd          =0.0
   if (associated(radiate%rshort))         radiate%rshort         =0.0
   if (associated(radiate%rlong))          radiate%rlong          =0.0
   if (associated(radiate%rlongup))        radiate%rlongup        =0.0
   if (associated(radiate%albedt))         radiate%albedt         =0.0
   if (associated(radiate%cosz))           radiate%cosz           =0.0
!Harrington fix
   if (associated(radiate%fthrd_lw      )) radiate%fthrd_lw       =0.0
   if (associated(radiate%rshort_top    )) radiate%rshort_top     =0.0
   if (associated(radiate%rshortup_top  )) radiate%rshortup_top   =0.0
   if (associated(radiate%rlongup_top   )) radiate%rlongup_top    =0.0
   if (associated(radiate%albedt_beam   )) radiate%albedt_beam    =0.0
   if (associated(radiate%albedt_diffuse)) radiate%albedt_diffuse =0.0
   if (associated(radiate%rlong_albedo  )) radiate%rlong_albedo   =0.0
   if (associated(radiate%rshort_diffuse)) radiate%rshort_diffuse =0.0

   return
   end subroutine zero_radiate

   subroutine dealloc_radiate(radiate)

   implicit none
   type (radiate_vars) :: radiate
   

   if (associated(radiate%fthrd))    deallocate (radiate%fthrd)
   if (associated(radiate%rshort))   deallocate (radiate%rshort)
   if (associated(radiate%rlong))    deallocate (radiate%rlong)
   if (associated(radiate%rlongup))  deallocate (radiate%rlongup)
   if (associated(radiate%albedt))   deallocate (radiate%albedt)
   if (associated(radiate%cosz))     deallocate (radiate%cosz)
!Harrington fix
   if (associated(radiate%fthrd_lw      )) deallocate (radiate%fthrd_lw      )
   if (associated(radiate%rshort_top    )) deallocate (radiate%rshort_top    )
   if (associated(radiate%rshortup_top  )) deallocate (radiate%rshortup_top  )
   if (associated(radiate%rlongup_top   )) deallocate (radiate%rlongup_top   )
   if (associated(radiate%albedt_beam   )) deallocate (radiate%albedt_beam   )
   if (associated(radiate%albedt_diffuse)) deallocate (radiate%albedt_diffuse)
   if (associated(radiate%rlong_albedo  )) deallocate (radiate%rlong_albedo  )
   if (associated(radiate%rshort_diffuse)) deallocate (radiate%rshort_diffuse)

   return
   end subroutine dealloc_radiate


   subroutine filltab_radiate(radiate,radiatem,imean,n1,n2,n3,ng)

   use var_tables

      implicit none
      type (radiate_vars) :: radiate,radiatem
      integer, intent(in) :: imean,n1,n2,n3,ng
      integer :: npts
      real, pointer :: var,varm

   ! Fill pointers to arrays into variable tables

      npts=n1*n2*n3

      if (associated(radiate%fthrd))  &
         call vtables2 (radiate%fthrd(1,1,1),radiatem%fthrd(1,1,1)  &
                    ,ng, npts, imean,  &
                    'FTHRD :3:hist:anal:mpti:mpt3')

      if (associated(radiate%fthrd_lw))  &
         call vtables2 (radiate%fthrd_lw(1,1,1),radiatem%fthrd_lw(1,1,1)  &
                    ,ng, npts, imean,  &
                    'FTHRD_LW :3:hist:anal:mpti:mpt3')

      npts=n2*n3
      if (associated(radiate%rshort))  &
         call vtables2 (radiate%rshort(1,1),radiatem%rshort(1,1)  &
                    ,ng, npts, imean,  &
                    'RSHORT :2:hist:anal:mpti:mpt3')
      if (associated(radiate%rlong))  &
         call vtables2 (radiate%rlong(1,1),radiatem%rlong(1,1)  &
                    ,ng, npts, imean,  &
                    'RLONG :2:hist:anal:mpti:mpt3')
      if (associated(radiate%rlongup))  &
         call vtables2 (radiate%rlongup(1,1),radiatem%rlongup(1,1)  &
                    ,ng, npts, imean,  &
                    'RLONGUP :2:hist:anal:mpti:mpt3')
      if (associated(radiate%albedt))  &
         call vtables2 (radiate%albedt(1,1),radiatem%albedt(1,1)  &
                    ,ng, npts, imean,  &
                    'ALBEDT :2:hist:anal:mpti:mpt3')
      if (associated(radiate%cosz))  &
         call vtables2 (radiate%cosz(1,1),radiatem%cosz(1,1)  &
                    ,ng, npts, imean,  &
                    'COSZ :2:anal:mpt3')

      return
      end subroutine filltab_radiate
!:DOC%END
End Module mem_radiate
