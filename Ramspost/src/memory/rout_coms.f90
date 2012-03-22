!==========================================================================================!
!==========================================================================================!
!     Module to allocate the output variables.                                             !
!------------------------------------------------------------------------------------------!
module rout_coms

   type rout_vars
      real   , pointer, dimension(:)          :: abuff
      real   , pointer, dimension(:)          :: bbuff
      real   , pointer, dimension(:,:)        :: r2
      real   , pointer, dimension(:,:,:)      :: r3
      real   , pointer, dimension(:,:,:,:)    :: r6
      real   , pointer, dimension(:,:,:)      :: r7
      real   , pointer, dimension(:,:,:,:)    :: r8
      real   , pointer, dimension(:,:,:)      :: r9
      real   , pointer, dimension(:,:,:)      :: r10
      integer, pointer, dimension(:,:)        :: iinf
      integer, pointer, dimension(:,:)        :: jinf
      real   , pointer, dimension(:,:,:)      :: rmi
      real   , pointer, dimension(:,:)        :: topo
      real   , pointer, dimension(:,:)        :: exner
      real   , pointer, dimension(:,:)        :: rlon
      real   , pointer, dimension(:,:)        :: rlat
      real   , pointer, dimension(:,:,:)      :: zplev
   end type rout_vars

   type(rout_vars), allocatable, dimension(:) :: rout
   type(rout_vars), allocatable, dimension(:) :: routgrads

   real           , parameter                 :: maxnormal =  1.e+06
   real           , parameter                 :: undefflg  = -1.e+34

   !=======================================================================================!
   !=======================================================================================!


   contains



   !=======================================================================================!
   !=======================================================================================!
   !      This routine allocates all buffers.                                              !
   !---------------------------------------------------------------------------------------!
   subroutine alloc_rout(this,nx,ny,nz,ngnd,npat,ncld,npl)
      implicit none
      !------ Arguments. ------------------------------------------------------------------!
      type(rout_vars), intent(inout) :: this
      integer        , intent(in)    :: nx
      integer        , intent(in)    :: ny
      integer        , intent(in)    :: nz
      integer        , intent(in)    :: ngnd
      integer        , intent(in)    :: npat
      integer        , intent(in)    :: ncld
      integer        , intent(in)    :: npl
      !------ Local variables. ------------------------------------------------------------!
      integer                        :: nbuff
      !------------------------------------------------------------------------------------!


      !----- Nullify all pointers. --------------------------------------------------------!
      call nullify_rout(this)
      !------------------------------------------------------------------------------------!


      !----- Find the maximum memory to allocate the generic buffers. ---------------------!
      nbuff = max(nx*ny*nz*ncld,nx*ny*nz*npat,nx*ny*ngnd*npat)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !      Then we can safely allocate them.                                             !
      !------------------------------------------------------------------------------------!
      allocate (this%abuff(nbuff)          )
      allocate (this%bbuff(nbuff)          )
      allocate (this%r2   (nx,ny)          )
      allocate (this%r3   (nx,ny,nz  )     )
      allocate (this%r6   (nx,ny,nz  ,ncld))
      allocate (this%r7   (nx,ny,npat)     )
      allocate (this%r8   (nx,ny,ngnd,npat))
      allocate (this%r9   (nx,ny,ncld)     )
      allocate (this%r10  (nx,ny,ngnd)     )
      allocate (this%iinf (nx,ny)          )
      allocate (this%jinf (nx,ny)          )
      allocate (this%rmi  (nx,ny,4)        )
      allocate (this%topo (nx,ny)          )
      allocate (this%exner(nx,ny)          )
      allocate (this%rlon (nx,ny)          )
      allocate (this%rlat (nx,ny)          )
      allocate (this%zplev(nx,ny,npl)      )
      !------------------------------------------------------------------------------------!
      
      call undef_rout (this,.true.)
      return
   end subroutine alloc_rout
   !=======================================================================================!
   !=======================================================================================!





   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine nullifies all pointers to ensure a safe allocation.               !
   !---------------------------------------------------------------------------------------!
   subroutine nullify_rout(this)
      implicit none
      !------ Arguments. ------------------------------------------------------------------!
      type(rout_vars), intent(inout) :: this
      !------------------------------------------------------------------------------------!


      !------ Nullify everything. ---------------------------------------------------------!
      nullify (this%abuff)
      nullify (this%bbuff)
      nullify (this%r2   )
      nullify (this%r3   )
      nullify (this%r6   )
      nullify (this%r7   )
      nullify (this%r8   )
      nullify (this%r9   )
      nullify (this%r10  )
      nullify (this%iinf )
      nullify (this%jinf )
      nullify (this%rmi  )
      nullify (this%topo )
      nullify (this%exner)
      nullify (this%rlon )
      nullify (this%rlat )
      nullify (this%zplev)
      !------------------------------------------------------------------------------------!

      return
   end subroutine nullify_rout
   !=======================================================================================!
   !=======================================================================================!





   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine deallocates all pointers.                                         !
   !---------------------------------------------------------------------------------------!
   subroutine undef_rout(this,everything)
      implicit none
      !------ Arguments. ------------------------------------------------------------------!
      type(rout_vars), intent(inout) :: this
      logical        , intent(in)    :: everything
      !------------------------------------------------------------------------------------!


      !------ Nullify everything. ---------------------------------------------------------!
      if (associated(this%abuff)) this%abuff = undefflg
      if (associated(this%bbuff)) this%bbuff = undefflg
      if (associated(this%r2   )) this%r2    = undefflg
      if (associated(this%r3   )) this%r3    = undefflg
      if (associated(this%r6   )) this%r6    = undefflg
      if (associated(this%r7   )) this%r7    = undefflg
      if (associated(this%r8   )) this%r8    = undefflg
      if (associated(this%r9   )) this%r9    = undefflg
      if (associated(this%r10  )) this%r10   = undefflg
      if (everything) then 
         if (associated(this%iinf )) this%iinf  = -1
         if (associated(this%jinf )) this%jinf  = -1
         if (associated(this%rmi  )) this%rmi   = undefflg
         if (associated(this%topo )) this%topo  = undefflg
         if (associated(this%exner)) this%exner = undefflg
         if (associated(this%rlon )) this%rlon  = undefflg
         if (associated(this%rlat )) this%rlat  = undefflg
         if (associated(this%zplev)) this%zplev = undefflg
      end if
      !------------------------------------------------------------------------------------!

      return
   end subroutine undef_rout
   !=======================================================================================!
   !=======================================================================================!





   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine deallocates all pointers.                                         !
   !---------------------------------------------------------------------------------------!
   subroutine dealloc_rout(this)
      implicit none
      !------ Arguments. ------------------------------------------------------------------!
      type(rout_vars), intent(inout) :: this
      !------------------------------------------------------------------------------------!


      !------ Nullify everything. ---------------------------------------------------------!
      if (associated(this%abuff)) deallocate (this%abuff)
      if (associated(this%bbuff)) deallocate (this%bbuff)
      if (associated(this%r2   )) deallocate (this%r2   )
      if (associated(this%r3   )) deallocate (this%r3   )
      if (associated(this%r6   )) deallocate (this%r6   )
      if (associated(this%r7   )) deallocate (this%r7   )
      if (associated(this%r8   )) deallocate (this%r8   )
      if (associated(this%r9   )) deallocate (this%r9   )
      if (associated(this%r10  )) deallocate (this%r10  )
      if (associated(this%iinf )) deallocate (this%iinf )
      if (associated(this%jinf )) deallocate (this%jinf )
      if (associated(this%rmi  )) deallocate (this%rmi  )
      if (associated(this%topo )) deallocate (this%topo )
      if (associated(this%exner)) deallocate (this%exner)
      if (associated(this%rlon )) deallocate (this%rlon )
      if (associated(this%rlat )) deallocate (this%rlat )
      if (associated(this%zplev)) deallocate (this%zplev)
      !------------------------------------------------------------------------------------!

      return
   end subroutine dealloc_rout
   !=======================================================================================!
   !=======================================================================================!
end module rout_coms
!==========================================================================================!
!==========================================================================================!
