!==========================================================================================!
!==========================================================================================!
! Copyright (C) 1991-2004  ; All Rights Reserved ; ATMET, LLC                              !
!                                                                                          !
! This file is free software; you can redistribute it and/or modify it under the           !
! terms of the GNU General Public License as published by the Free Software                !
! Foundation; either version 2 of the License, or (at your option) any later version.      !
!                                                                                          !
! This software is distributed in the hope that it will be useful, but WITHOUT ANY         !
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A          !
! PARTICULAR PURPOSE.  See the GNU General Public License for more details.                !
!                                                                                          !
! You should have received a copy of the GNU General Public License along with this        !
! program; if not, write to the Free Software Foundation, Inc.,                            !
! 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.                                 !
!==========================================================================================!
!==========================================================================================!
module an_header

   use mod_maxdims , only: maxstr,maxrank,maxrank
   use mod_time    , only: time_stt
   !----- Number of files (group of simultaneous times in BRAMS case) ---------------------!
   integer                                       :: nfiles
   type head_table

      !----- The following variables have a single value per file -------------------------!
      character(len=maxstr)                      :: filename   ! File name
      integer                                    :: nvars      ! Number of variables
      integer                                    :: ntimes     ! Number of times per file
      integer                                    :: ngrids     ! Number of grids in here 
      type(time_stt)                             :: init_time  ! Initial time
      !----- The following variables are dimensioned by nvars -----------------------------!
      character(len=16), pointer, dimension(:)   :: varname    ! Variable name
      integer          , pointer, dimension(:)   :: npointer   ! Pointer (offset) - (B)RAMS
      integer          , pointer, dimension(:)   :: idim_type  ! Dimension type   - (B)RAMS
      integer          , pointer, dimension(:)   :: ngrid      ! Current grid
      integer          , pointer, dimension(:)   :: nvalues    ! # of values      - (B)RAMS
      integer          , pointer, dimension(:)   :: rank       ! Rank of this dataset
      integer          , pointer, dimension(:,:) :: dims       ! Size of each dimension
      logical          , pointer, dimension(:,:) :: stagger    ! Is it staggered  - WRF
      !----- The following variables are dimensioned by ntimes ----------------------------!
      type(time_stt)   , pointer, dimension(:)   :: file_time  ! File time
      !----- The following variable is dimensioned by ngrids ------------------------------!
      logical          , pointer, dimension(:)   :: avail_grid ! Available grids flag.

   end type head_table

   type (head_table), allocatable, dimension(:)  :: info_table

   contains 
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine dynamically allocates the dimensions for this file. DO NOT CALL   !
   ! THIS SUBROUTINE BEFORE ASSIGNING NVARS AND NTIMES, OTHERWISE HORRIBLE THINGS WILL     !
   ! HAPPEN!!!                                                                             !
   !---------------------------------------------------------------------------------------!
   subroutine alloc_anheader(itable)
      implicit none

      type(head_table), intent(inout) :: itable
      integer                         :: invars,intimes,ingrids
      
      invars  = itable%nvars
      intimes = itable%ntimes
      ingrids = itable%ngrids

      allocate(itable%varname            (invars)         )
      allocate(itable%npointer           (invars)         )
      allocate(itable%idim_type          (invars)         )
      allocate(itable%ngrid              (invars)         )
      allocate(itable%nvalues            (invars)         )
      allocate(itable%rank               (invars)         )
      allocate(itable%dims       (maxrank,invars)         )
      allocate(itable%stagger    (maxrank,invars)         )
      
      allocate(itable%file_time                  (intimes))

      allocate(itable%avail_grid                 (ingrids))

      return
   end subroutine alloc_anheader
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine safely nullifies the an_table pointers before allocating anything. !
   !---------------------------------------------------------------------------------------!
   subroutine nullify_anheader(itable)
      type(head_table), intent(inout) :: itable 

      if(associated(itable%varname    )) nullify(itable%varname    )
      if(associated(itable%npointer   )) nullify(itable%npointer   )
      if(associated(itable%idim_type  )) nullify(itable%idim_type  )
      if(associated(itable%ngrid      )) nullify(itable%ngrid      )
      if(associated(itable%nvalues    )) nullify(itable%nvalues    )
      if(associated(itable%rank       )) nullify(itable%rank       )
      if(associated(itable%dims       )) nullify(itable%dims       )
      if(associated(itable%stagger    )) nullify(itable%stagger    )

      if(associated(itable%file_time  )) nullify(itable%file_time  )

      if(associated(itable%avail_grid )) nullify(itable%avail_grid )
      return
   end subroutine nullify_anheader
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine deallocates the header structure.                                  !
   !---------------------------------------------------------------------------------------!
   subroutine dealloc_anheader(itable)
      type(head_table), intent(inout) :: itable 

      if(associated(itable%varname    )) deallocate(itable%varname    )
      if(associated(itable%npointer   )) deallocate(itable%npointer   )
      if(associated(itable%idim_type  )) deallocate(itable%idim_type  )
      if(associated(itable%ngrid      )) deallocate(itable%ngrid      )
      if(associated(itable%nvalues    )) deallocate(itable%nvalues    )
      if(associated(itable%rank       )) deallocate(itable%rank       )
      if(associated(itable%dims       )) deallocate(itable%dims       )

      if(associated(itable%file_time  )) deallocate(itable%file_time  )

      if(associated(itable%avail_grid )) deallocate(itable%avail_grid )

      return
   end subroutine dealloc_anheader

end module an_header
!==========================================================================================!
!==========================================================================================!
