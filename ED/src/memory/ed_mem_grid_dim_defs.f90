module ed_mem_grid_dim_defs
  implicit none

  public

  integer, pointer :: nmzp(:), nmxp(:),nmyp(:)
  integer          :: maxz,    maxx,   maxy

  contains

    subroutine define_grid_dim_pointer(proc_type, ngrids, maxgrds, &
         nnxp, nnyp, mmxp, mmyp)

      implicit none
      ! Arguments
      integer, intent(in) :: proc_type
      integer, intent(in) :: ngrids
      integer, intent(in) :: maxgrds
      integer, target, intent(in) :: nnxp(maxgrds)
      integer, target, intent(in) :: nnyp(maxgrds)
      integer, target, intent(in) :: mmxp(maxgrds)
      integer, target, intent(in) :: mmyp(maxgrds)
      ! Local Variables:
      integer             :: ng

      ! Defining pointers
      if (proc_type == 0 .or. proc_type == 1) then
         !  This is the call for either a single processor run or
         !    for the master process
         nmxp => nnxp
         nmyp => nnyp
      elseif (proc_type == 2) then
         !  This is the call for a initial compute node process
         nmxp => mmxp
         nmyp => mmyp

      endif

      ! Defining max grid dimensions
      maxx = nmxp(1)
      maxy = nmyp(1)

      do ng=2, ngrids
         if (nmxp(ng) > maxx) maxx = nmxp(ng)
         if (nmyp(ng) > maxy) maxy = nmyp(ng)
      enddo

    end subroutine define_grid_dim_pointer

end module ed_mem_grid_dim_defs
