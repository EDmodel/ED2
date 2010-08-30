!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################


module mem_oda


  type oda_vars

     ! Variables to be dimensioned by (nzp,nxp,nyp)
     real, pointer, dimension(:,:,:) :: &
          uk,vk,tk,rk,ukv,vkv,tkv,rkv

  end type oda_vars

  type (oda_vars), allocatable :: oda_g(:), odam_g(:)



  integer, parameter :: maxodafiles=1000, maxodasta=2000, maxodagrids=10
  integer, parameter :: maxodanzp=200, maxodatimes=3*maxodafiles

  character(len=128), dimension(maxodafiles) :: fnames_upa, fnames_sfc
  character(len=14) , dimension(maxodafiles) :: itotdate_upa,itotdate_sfc
  character(len=8)  , dimension(maxodasta)   :: staid_sfc, staid_upa
  integer, dimension(maxodasta)              :: ntimes_sfc, ntimes_upa
  integer                                    :: maxtimes_sfc, maxtimes_upa

  ! Namelist inputs

  character(len=128) :: oda_upaprefix,oda_sfcprefix
  integer :: if_oda
  real :: frqoda,todabeg,todaend,tnudoda,wt_oda_grid(maxodagrids)  &
       ,oda_sfc_til,oda_sfc_tel,oda_upa_til,oda_upa_tel  &
       ,wt_oda_uv,wt_oda_th,wt_oda_pi,wt_oda_rt
  real, dimension(maxodagrids) :: roda_sfce, roda_sfc0  &
       , roda_upae, roda_upa0  &
       , roda_zfact, roda_hgt


  integer ::  nupafiles,nsfcfiles,num_oda_sfc,num_oda_upa

  ! Surface data 

  type oda_sfc_info_type
     character(len=8) :: id
     integer :: intid
     integer :: ntimes
     integer, dimension(maxodagrids) :: iactive
     real, dimension(maxodagrids) :: xista, xjsta
     real :: xlat,xlon,xsta,ysta,stopo
  end type oda_sfc_info_type


  type oda_sfc_type
     real, pointer, dimension(:) :: temp, dewpt, us, vs, ps,u,v 
     real(kind=8), pointer, dimension(:) :: time 
  end type oda_sfc_type

  type(oda_sfc_info_type), allocatable :: oda_sfc_info(:)
  type(oda_sfc_type)     , allocatable :: oda_sfc_obs(:)


  ! Upper air info

  type oda_upa_info_type
     character(len=8) :: id
     integer :: intid
     integer :: ntimes
     integer, dimension(maxodagrids) :: iactive
     real, dimension(maxodagrids) :: xista, xjsta
     real :: xlat,xlon,xsta,ysta,stopo
  end type oda_upa_info_type

  ! Upper air data
  type oda_upa_type
     character(len=14) :: ctotdate
     real, pointer, dimension(:,:) :: theta, rv, us, vs, u, v, zz, pi, zgeo
     real(kind=8), pointer, dimension(:) :: time 
     integer, pointer, dimension(:) :: lp,lz
  end type oda_upa_type

  type(oda_upa_info_type), allocatable :: oda_upa_info(:)
  type(oda_upa_type)     , allocatable :: oda_upa_obs(:)

  integer, parameter :: maxupalevs=6000

  ! Krigging routine info

  real :: rmaxkrg(maxodanzp,maxodagrids)  &
       ,ckrg(3,maxodagrids),akrg(maxodanzp,maxodagrids)  &
       ,caxkrg(9,maxodagrids),caykrg(9,maxodagrids),cazkrg(9,maxodagrids)
  integer :: nstkrg(maxodagrids)


  ! Filled obs arrays for an analysis time

  integer, parameter :: maxkobs=10000
  real, dimension(maxkobs) :: ukobs,vkobs,tkobs,rkobs,pkobs  &
       ,xkobs,ykobs,zkobs,ekobs  &
       ,ikobs,jkobs

contains

  subroutine alloc_oda(oda,n1,n2,n3,ng,proc_type)

    implicit none
    type (oda_vars) :: oda
    integer, intent(in) :: n1,n2,n3,ng,proc_type

    ! Allocate arrays based on options (if necessary)


    if( if_oda == 1 .and. proc_type /= 1) then
       allocate (oda%uk(n1,n2,n3))
       allocate (oda%vk(n1,n2,n3))
       allocate (oda%tk(n1,n2,n3))
       allocate (oda%rk(n1,n2,n3))
       allocate (oda%ukv(n1,n2,n3))
       allocate (oda%vkv(n1,n2,n3))
       allocate (oda%tkv(n1,n2,n3))
       allocate (oda%rkv(n1,n2,n3))
    endif

    return
  end subroutine alloc_oda


  subroutine nullify_oda(oda)

    implicit none
    type (oda_vars) :: oda

    if (associated(oda%uk))      nullify (oda%uk)
    if (associated(oda%vk))      nullify (oda%vk)
    if (associated(oda%tk))      nullify (oda%tk)
    if (associated(oda%rk))      nullify (oda%rk)
    if (associated(oda%ukv))     nullify (oda%ukv)
    if (associated(oda%vkv))     nullify (oda%vkv)
    if (associated(oda%tkv))     nullify (oda%tkv)
    if (associated(oda%rkv))     nullify (oda%rkv)

    return
  end subroutine nullify_oda

  subroutine dealloc_oda(oda)

    implicit none
    type (oda_vars) :: oda


    if (associated(oda%uk))      deallocate (oda%uk)
    if (associated(oda%vk))      deallocate (oda%vk)
    if (associated(oda%tk))      deallocate (oda%tk)
    if (associated(oda%rk))      deallocate (oda%rk)
    if (associated(oda%ukv))     deallocate (oda%ukv)
    if (associated(oda%vkv))     deallocate (oda%vkv)
    if (associated(oda%tkv))     deallocate (oda%tkv)
    if (associated(oda%rkv))     deallocate (oda%rkv)

    return
  end subroutine dealloc_oda


  subroutine filltab_oda(oda,odam,imean,n1,n2,n3,ng)

    use var_tables

    implicit none
    type (oda_vars) :: oda,odam
    integer, intent(in) :: imean,n1,n2,n3,ng
    integer :: npts
    real, pointer :: var,varm

    ! Fill pointers to arrays into variable tables

    npts=n1*n2*n3

    if (associated(oda%uk))  &
         call vtables2 (oda%uk(1,1,1),odam%uk(1,1,1)  &
         ,ng, npts, imean,  &
         'UKODA :3:')
    if (associated(oda%vk))  &
         call vtables2 (oda%vk(1,1,1),odam%vk(1,1,1)  &
         ,ng, npts, imean,  &
         'VKODA :3:')
    if (associated(oda%tk))  &
         call vtables2 (oda%tk(1,1,1),odam%tk(1,1,1)  &
         ,ng, npts, imean,  &
         'TKODA :3:')
    if (associated(oda%rk))  &
         call vtables2 (oda%rk(1,1,1),odam%rk(1,1,1)  &
         ,ng, npts, imean,  &
         'RKODA :3:')
    if (associated(oda%ukv))  &
         call vtables2 (oda%ukv(1,1,1),odam%ukv(1,1,1)  &
         ,ng, npts, imean,  &
         'UVODA :3:')
    if (associated(oda%vkv))  &
         call vtables2 (oda%vkv(1,1,1),odam%vkv(1,1,1)  &
         ,ng, npts, imean,  &
         'VVODA :3:')
    if (associated(oda%tkv))  &
         call vtables2 (oda%tkv(1,1,1),odam%tkv(1,1,1)  &
         ,ng, npts, imean,  &
         'TVODA :3:')
    if (associated(oda%rkv))  &
         call vtables2 (oda%rkv(1,1,1),odam%rkv(1,1,1)  &
         ,ng, npts, imean,  &
         'RVODA :3:')

    return
  end subroutine filltab_oda

end module mem_oda
