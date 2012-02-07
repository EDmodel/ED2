!==========================================================================================!
!==========================================================================================!
!    Module mem_mnt_advec.f90                                                              !
!                                                                                          !
!    This module holds several variables that will help in the advection schemes.          !
!------------------------------------------------------------------------------------------!
module mem_mnt_advec
   implicit none

   !---------------------------------------------------------------------------------------!
   !     Structure with the auxiliary variables to be used by the advection scheme.        !
   !---------------------------------------------------------------------------------------!
   type advec_vars
      !----- scratch 3d wind variables to be used within the advection scheme. ------------!
      real,pointer,dimension  (:,:,:)  :: uavg
      real,pointer,dimension  (:,:,:)  :: vavg
      real,pointer,dimension  (:,:,:)  :: wavg
      !------------------------------------------------------------------------------------!


      !----- scratch 3d scalars to be used within the advection scheme. -------------------!
      real,pointer,dimension  (:,:,:)  :: scal_in
      real,pointer,dimension  (:,:,:)  :: scal_out
      !------------------------------------------------------------------------------------!



      !----- 3D version of density for the advection scheme. ------------------------------!
      real,pointer,dimension  (:,:,:)  :: denst
      real,pointer,dimension  (:,:,:)  :: densu
      real,pointer,dimension  (:,:,:)  :: densv
      real,pointer,dimension  (:,:,:)  :: densw
      !------------------------------------------------------------------------------------!


      !----- Densities as definied in the Walcek papers. ----------------------------------!
      real,pointer,dimension  (:,:,:)  :: den0_wal
      real,pointer,dimension  (:,:,:)  :: den1_wal
      real,pointer,dimension  (:,:,:)  :: den2_wal
      real,pointer,dimension  (:,:,:)  :: den3_wal
      !------------------------------------------------------------------------------------!



      !----- Grid spacing. ----------------------------------------------------------------!
      real,pointer,dimension  (:,:,:)    :: dxtw
      real,pointer,dimension  (:,:,:)    :: dytw
      real,pointer,dimension  (:,:,:)    :: dztw
      !------------------------------------------------------------------------------------!
   end type advec_vars
   type(advec_vars), dimension(:), allocatable   :: advec_g
   type(advec_vars), dimension(:), allocatable   :: advecm_g
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Variable to be read in from the namelist.                                         !
   !---------------------------------------------------------------------------------------!
   integer :: iadvec   ! 0 -- original advection scheme
                       ! 1 -- monotonic advection (Freitas et al, in press, JAMES)
   !---------------------------------------------------------------------------------------!


   !=======================================================================================!
   !=======================================================================================!


   contains



   !=======================================================================================!
   !=======================================================================================!
   subroutine alloc_advec(advec,n1,n2,n3)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(advec_vars), intent(inout) :: advec
      integer         , intent(in)    :: n1
      integer         , intent(in)    :: n2
      integer         , intent(in)    :: n3
      !------------------------------------------------------------------------------------!

      
      allocate(advec%uavg      (n1,n2,n3))
      allocate(advec%vavg      (n1,n2,n3))
      allocate(advec%wavg      (n1,n2,n3))
      allocate(advec%scal_in   (n1,n2,n3))
      allocate(advec%scal_out  (n1,n2,n3))
      allocate(advec%denst     (n1,n2,n3))
      allocate(advec%densu     (n1,n2,n3))
      allocate(advec%densv     (n1,n2,n3))
      allocate(advec%densw     (n1,n2,n3))
      allocate(advec%den0_wal  (n1,n2,n3))
      allocate(advec%den1_wal  (n1,n2,n3))
      allocate(advec%den2_wal  (n1,n2,n3))
      allocate(advec%den3_wal  (n1,n2,n3))
      allocate(advec%dxtw      (n1,n2,n3))
      allocate(advec%dytw      (n1,n2,n3))
      allocate(advec%dztw      (n1,n2,n3))

      return
   end subroutine alloc_advec
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   subroutine nullify_advec(advec)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(advec_vars), intent(inout) :: advec
      !------------------------------------------------------------------------------------!

      
      if (associated(advec%uavg      )) nullify(advec%uavg      )
      if (associated(advec%vavg      )) nullify(advec%vavg      )
      if (associated(advec%wavg      )) nullify(advec%wavg      )
      if (associated(advec%scal_in   )) nullify(advec%scal_in   )
      if (associated(advec%scal_out  )) nullify(advec%scal_out  )
      if (associated(advec%denst     )) nullify(advec%denst     )
      if (associated(advec%densu     )) nullify(advec%densu     )
      if (associated(advec%densv     )) nullify(advec%densv     )
      if (associated(advec%densw     )) nullify(advec%densw     )
      if (associated(advec%den0_wal  )) nullify(advec%den0_wal  )
      if (associated(advec%den1_wal  )) nullify(advec%den1_wal  )
      if (associated(advec%den2_wal  )) nullify(advec%den2_wal  )
      if (associated(advec%den3_wal  )) nullify(advec%den3_wal  )
      if (associated(advec%dxtw      )) nullify(advec%dxtw      )
      if (associated(advec%dytw      )) nullify(advec%dytw      )
      if (associated(advec%dztw      )) nullify(advec%dztw      )

      return
   end subroutine nullify_advec
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   subroutine zero_advec(advec)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(advec_vars), intent(inout) :: advec
      !------------------------------------------------------------------------------------!

      
      if (associated(advec%uavg      )) advec%uavg      = 0.0
      if (associated(advec%vavg      )) advec%vavg      = 0.0
      if (associated(advec%wavg      )) advec%wavg      = 0.0
      if (associated(advec%scal_in   )) advec%scal_in   = 0.0
      if (associated(advec%scal_out  )) advec%scal_out  = 0.0
      if (associated(advec%denst     )) advec%denst     = 0.0
      if (associated(advec%densu     )) advec%densu     = 0.0
      if (associated(advec%densv     )) advec%densv     = 0.0
      if (associated(advec%densw     )) advec%densw     = 0.0
      if (associated(advec%den0_wal  )) advec%den0_wal  = 0.0
      if (associated(advec%den1_wal  )) advec%den1_wal  = 0.0
      if (associated(advec%den2_wal  )) advec%den2_wal  = 0.0
      if (associated(advec%den3_wal  )) advec%den3_wal  = 0.0
      if (associated(advec%dxtw      )) advec%dxtw      = 0.0
      if (associated(advec%dytw      )) advec%dytw      = 0.0
      if (associated(advec%dztw      )) advec%dztw      = 0.0

      return
   end subroutine zero_advec
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   subroutine dealloc_advec(advec)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(advec_vars), intent(inout) :: advec
      !------------------------------------------------------------------------------------!

      
      if (associated(advec%uavg      )) deallocate(advec%uavg      )
      if (associated(advec%vavg      )) deallocate(advec%vavg      )
      if (associated(advec%wavg      )) deallocate(advec%wavg      )
      if (associated(advec%scal_in   )) deallocate(advec%scal_in   )
      if (associated(advec%scal_out  )) deallocate(advec%scal_out  )
      if (associated(advec%denst     )) deallocate(advec%denst     )
      if (associated(advec%densu     )) deallocate(advec%densu     )
      if (associated(advec%densv     )) deallocate(advec%densv     )
      if (associated(advec%densw     )) deallocate(advec%densw     )
      if (associated(advec%den0_wal  )) deallocate(advec%den0_wal  )
      if (associated(advec%den1_wal  )) deallocate(advec%den1_wal  )
      if (associated(advec%den2_wal  )) deallocate(advec%den2_wal  )
      if (associated(advec%den3_wal  )) deallocate(advec%den3_wal  )
      if (associated(advec%dxtw      )) deallocate(advec%dxtw      )
      if (associated(advec%dytw      )) deallocate(advec%dytw      )
      if (associated(advec%dztw      )) deallocate(advec%dztw      )

      return
   end subroutine dealloc_advec
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   subroutine filltab_advec(advec,advecm,imean,n1,n2,n3,ng)
      use var_tables

      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(advec_vars)   , intent(in) :: advec
      type(advec_vars)   , intent(in) :: advecm
      integer            , intent(in) :: imean
      integer            , intent(in) :: n1
      integer            , intent(in) :: n2
      integer            , intent(in) :: n3
      integer            , intent(in) :: ng
      !----- Local variables. -------------------------------------------------------------!
      integer                         :: npts
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !      3-D variables, exchanged at the main advectiom time step for.                                                                !
      !------------------------------------------------------------------------------------!
      npts = n1 * n2 * n3
      
      if (associated(advec%uavg      ))                                                    &
         call vtables2(advec%uavg,advecm%uavg,ng,npts,imean                                &
                      ,'UAVG :3:mpti:mpt3:advu')

      if (associated(advec%vavg      ))                                                    &
         call vtables2(advec%vavg,advecm%vavg,ng,npts,imean                                &
                      ,'VAVG :3:mpti:mpt3:advv')

      if (associated(advec%wavg      ))                                                    &
         call vtables2(advec%wavg,advecm%wavg,ng,npts,imean                                &
                      ,'WAVG :3:mpti:mpt3:advw')

      if (associated(advec%denst     ))                                                    &
         call vtables2(advec%denst,advecm%denst,ng,npts,imean                              &
                      ,'DENST :3:mpti:mpt3:advt')

      if (associated(advec%densu     ))                                                    &
         call vtables2(advec%densu,advecm%densu,ng,npts,imean                              &
                      ,'DENSU :3:mpti:mpt3:advu')

      if (associated(advec%densv     ))                                                    &
         call vtables2(advec%densv,advecm%densv,ng,npts,imean                              &
                      ,'DENSV :3:mpti:mpt3:advv')

      if (associated(advec%densw     ))                                                    &
         call vtables2(advec%densw,advecm%densw,ng,npts,imean                              &
                      ,'DENSW :3:mpti:mpt3:advw')

      if (associated(advec%den0_wal  ))                                                    &
         call vtables2(advec%den0_wal,advecm%den0_wal,ng,npts,imean                        &
                      ,'DEN0_WAL :3:mpti:mpt3:advt')

      if (associated(advec%den1_wal  ))                                                    &
         call vtables2(advec%den1_wal,advecm%den1_wal,ng,npts,imean                        &
                      ,'DEN1_WAL :3:mpti:mpt3:advt')

      if (associated(advec%den2_wal  ))                                                    &
         call vtables2(advec%den2_wal,advecm%den2_wal,ng,npts,imean                        &
                      ,'DEN2_WAL :3:mpti:mpt3:advt')

      if (associated(advec%den3_wal  ))                                                    &
         call vtables2(advec%den3_wal,advecm%den3_wal,ng,npts,imean                        &
                      ,'DEN3_WAL :3:mpti:mpt3:advt')
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !      Grid variables, save as 3-D variables for convenience.  They are exchanged    !
      ! only at the beginning.                                                             !
      !------------------------------------------------------------------------------------!
      if (associated(advec%dxtw      ))                                                    &
         call vtables2(advec%dxtw,advecm%dxtw,ng,npts,imean                                &
                      ,'DXTW :3:mpti:mpt3')

      if (associated(advec%dytw      ))                                                    &
         call vtables2(advec%dytw,advecm%dytw,ng,npts,imean                                &
                      ,'DYTW :3:mpti:mpt3')

      if (associated(advec%dztw      ))                                                    &
         call vtables2(advec%dztw,advecm%dztw,ng,npts,imean                                &
                      ,'DZTW :3:mpti:mpt3')
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Scalar temporary array.  We copy them at specific times within the time step. !
      !------------------------------------------------------------------------------------!
      if (associated(advec%scal_in   ))                                                    &
         call vtables2(advec%scal_in,advecm%scal_in,ng,npts,imean                          &
                      ,'SCAL_IN :3:mpti:mpt3')

      if (associated(advec%scal_out  ))                                                    &
         call vtables2(advec%scal_out,advecm%scal_out,ng,npts,imean                        &
                      ,'SCAL_OUT :3:mpti:mpt3')
      !------------------------------------------------------------------------------------!

      return
   end subroutine filltab_advec
   !=======================================================================================!
   !=======================================================================================!
end module mem_mnt_advec
!==========================================================================================!
!==========================================================================================!
