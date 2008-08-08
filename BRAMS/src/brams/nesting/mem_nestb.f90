!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################


module mem_nestb

  type nest_bounds

     real, pointer, dimension(:,:,:)   :: bux,buy,buz,bvx,bvy,bvz  &
          ,bwx,bwy,bwz,bpx,bpy,bpz
     real, pointer, dimension(:,:,:,:) :: bsx,bsy,bsz

  end type nest_bounds

!!!!!!   use max_dims, only : maxgrds,maxmach
  !  The following are here only until all include's are modules !!!!!!!
  integer, parameter, private :: maxgrds=8

  type (nest_bounds) :: nbounds(maxgrds)


contains

  subroutine alloc_nestb(ng,nx,ny,nz)

    use var_tables

    implicit none

    integer :: ng,nx,ny,nz

    !  Allocate "b" array components. All grids will be allocated,
    !     only to 1's if nesting isn't done.

    allocate( nbounds(ng)%bux(nz,ny,2) )
    allocate( nbounds(ng)%buy(nz,nx,2) )
    allocate( nbounds(ng)%buz(nx,ny,2) )

    allocate( nbounds(ng)%bvx(nz,ny,2) )
    allocate( nbounds(ng)%bvy(nz,nx,2) )
    allocate( nbounds(ng)%bvz(nx,ny,2) )

    allocate( nbounds(ng)%bwx(nz,ny,2) )
    allocate( nbounds(ng)%bwy(nz,nx,2) )
    allocate( nbounds(ng)%bwz(nx,ny,2) )

    allocate( nbounds(ng)%bpx(nz,ny,2) )
    allocate( nbounds(ng)%bpy(nz,nx,2) )
    allocate( nbounds(ng)%bpz(nx,ny,2) )

    allocate( nbounds(ng)%bsx(nz,ny,2,num_scalar(ng)) )
    allocate( nbounds(ng)%bsy(nz,nx,2,num_scalar(ng)) )
    allocate( nbounds(ng)%bsz(nx,ny,2,num_scalar(ng)) )

    ! ALF - Putting Zero on all nest variables

    call azero((nz*ny*2), nbounds(ng)%bux(1,1,1))
    call azero((nz*nx*2), nbounds(ng)%buy(1,1,1))
    call azero((nx*ny*2), nbounds(ng)%buz(1,1,1))

    call azero((nz*ny*2), nbounds(ng)%bvx(1,1,1))
    call azero((nz*nx*2), nbounds(ng)%bvy(1,1,1))
    call azero((nx*ny*2), nbounds(ng)%bvz(1,1,1))

    call azero((nz*ny*2), nbounds(ng)%bwx(1,1,1))
    call azero((nz*nx*2), nbounds(ng)%bwy(1,1,1))
    call azero((nx*ny*2), nbounds(ng)%bwz(1,1,1))

    call azero((nz*ny*2), nbounds(ng)%bpx(1,1,1))
    call azero((nz*nx*2), nbounds(ng)%bpy(1,1,1))
    call azero((nx*ny*2), nbounds(ng)%bpz(1,1,1))

    call azero((nz*ny*2*num_scalar(ng)), nbounds(ng)%bsx(1,1,1,1))
    call azero((nz*nx*2*num_scalar(ng)), nbounds(ng)%bsy(1,1,1,1))
    call azero((nx*ny*2*num_scalar(ng)), nbounds(ng)%bsz(1,1,1,1))

    return
  end subroutine alloc_nestb

end module mem_nestb
