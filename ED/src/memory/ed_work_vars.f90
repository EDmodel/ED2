module ed_work_vars

   use ed_max_dims, only : maxgrds

   type work_vars
      ! Variables to be dimensioned by (nxp,nyp)
      real   , dimension(  :,:), pointer :: glon
      real   , dimension(  :,:), pointer :: glat
      real   , dimension(  :,:), pointer :: work
      logical, dimension(  :,:), pointer :: land
      real   , dimension(  :,:), pointer :: landfrac
      real   , dimension(:,:,:), pointer :: soilfrac
      integer, dimension(:,:,:), pointer :: ntext
      integer, dimension(  :,:), pointer :: xatm
      integer, dimension(  :,:), pointer :: yatm
   end type work_vars

   type work_vecs
      ! Polygon vectors
      real   , dimension(:  ), pointer :: glon
      real   , dimension(:  ), pointer :: glat
      real   , dimension(:  ), pointer :: work
      real   , dimension(:  ), pointer :: landfrac
      real   , dimension(:,:), pointer :: soilfrac
      integer, dimension(:,:), pointer :: ntext
      integer, dimension(:  ), pointer :: xid
      integer, dimension(:  ), pointer :: yid
   end type work_vecs


   type (work_vars), dimension(:), allocatable :: work_e
   type (work_vecs), dimension(:), allocatable :: work_v

   !----- Auxiliary variable to count number of polygons in this run. ---------------------!
   integer, dimension(maxgrds)  :: npolys_run

   !=======================================================================================!
   !=======================================================================================!


   contains



   !=======================================================================================!
   !=======================================================================================!
   subroutine ed_alloc_work(worke,n2,n3,nsite)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type (work_vars), intent(inout) :: worke
      integer         , intent(in)    :: n2
      integer         , intent(in)    :: n3
      integer         , intent(in)    :: nsite
      !------------------------------------------------------------------------------------!



      !----- Allocate arrays based on options (if necessary). -----------------------------!
      allocate (worke%glon    (      n2,n3))
      allocate (worke%glat    (      n2,n3))
      allocate (worke%xatm    (      n2,n3))
      allocate (worke%yatm    (      n2,n3))
      allocate (worke%work    (      n2,n3))
      allocate (worke%land    (      n2,n3))
      allocate (worke%landfrac(      n2,n3))
      allocate (worke%soilfrac(nsite,n2,n3))
      allocate (worke%ntext   (nsite,n2,n3))
      !------------------------------------------------------------------------------------!

      return
   end subroutine ed_alloc_work
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   subroutine ed_nullify_work(worke)
      implicit none

      !----- Arguments. -------------------------------------------------------------------!
      type (work_vars) :: worke
      !------------------------------------------------------------------------------------!

      if (associated(worke%glon     )) nullify (worke%glon    )
      if (associated(worke%glat     )) nullify (worke%glat    )
      if (associated(worke%xatm     )) nullify (worke%xatm    )
      if (associated(worke%yatm     )) nullify (worke%yatm    )
      if (associated(worke%work     )) nullify (worke%work    )
      if (associated(worke%land     )) nullify (worke%land    )
      if (associated(worke%landfrac )) nullify (worke%landfrac)
      if (associated(worke%soilfrac )) nullify (worke%soilfrac)
      if (associated(worke%ntext    )) nullify (worke%ntext   )

      return
   end subroutine ed_nullify_work
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   subroutine ed_dealloc_work(worke)
      implicit none
      type (work_vars) :: worke

      if (associated(worke%glon        )) deallocate (worke%glon        )
      if (associated(worke%glat        )) deallocate (worke%glat        )
      if (associated(worke%xatm        )) deallocate (worke%xatm        )
      if (associated(worke%yatm        )) deallocate (worke%yatm        )
      if (associated(worke%work        )) deallocate (worke%work        )
      if (associated(worke%land        )) deallocate (worke%land        )
      if (associated(worke%landfrac    )) deallocate (worke%landfrac    )
      if (associated(worke%soilfrac    )) deallocate (worke%soilfrac    )
      if (associated(worke%ntext       )) deallocate (worke%ntext       )
      return
   end subroutine ed_dealloc_work
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   subroutine ed_alloc_work_vec(workv,npolys,nsite)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type (work_vecs), intent(inout) :: workv
      integer         , intent(in)    :: npolys
      integer         , intent(in)    :: nsite
      !------------------------------------------------------------------------------------!



      !----- Allocate arrays based on options (if necessary). -----------------------------!
      allocate (workv%glon    (      npolys))
      allocate (workv%glat    (      npolys))
      allocate (workv%work    (      npolys))
      allocate (workv%landfrac(      npolys))
      allocate (workv%soilfrac(nsite,npolys))
      allocate (workv%ntext   (nsite,npolys))
      allocate (workv%xid     (      npolys))
      allocate (workv%yid     (      npolys))
      !------------------------------------------------------------------------------------!

      return
   end subroutine ed_alloc_work_vec
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   subroutine ed_nullify_work_vec(workv)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type (work_vecs) :: workv
      !------------------------------------------------------------------------------------!


      if (associated(workv%glon    ))   nullify(workv%glon    )
      if (associated(workv%glat    ))   nullify(workv%glat    )
      if (associated(workv%work    ))   nullify(workv%work    )
      if (associated(workv%landfrac))   nullify(workv%landfrac)
      if (associated(workv%soilfrac))   nullify(workv%soilfrac)
      if (associated(workv%ntext   ))   nullify(workv%ntext   )
      if (associated(workv%xid     ))   nullify(workv%xid     )
      if (associated(workv%yid     ))   nullify(workv%yid     )

      return
   end subroutine ed_nullify_work_vec
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   subroutine ed_dealloc_work_vec(workv)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type (work_vecs) :: workv
      !------------------------------------------------------------------------------------!


      if (associated(workv%glon    ))   deallocate(workv%glon    )
      if (associated(workv%glat    ))   deallocate(workv%glat    )
      if (associated(workv%work    ))   deallocate(workv%work    )
      if (associated(workv%landfrac))   deallocate(workv%landfrac)
      if (associated(workv%soilfrac))   deallocate(workv%soilfrac)
      if (associated(workv%ntext   ))   deallocate(workv%ntext   )
      if (associated(workv%xid     ))   deallocate(workv%xid     )
      if (associated(workv%yid     ))   deallocate(workv%yid     )

      return
   end subroutine ed_dealloc_work_vec
   !=======================================================================================!
   !=======================================================================================!

end module ed_work_vars
!==========================================================================================!
!==========================================================================================!
