!==========================================================================================!
!   Module ed_var_tables: this module contains all information about variables, includ-    !
! ing the kind of variable and whether they need to be passed in parallel runs, in which   !
! output they should be included etc.                                                      !
!   The following list is the number used to identify the type of variable.                !
!                                                                                          !
!     !----- POLYGON: --------------------------------------!                              !
!     !  10  : rank 1 : integer                             !                              !
!     !  11  : rank 1 : real                                !                              !
!     ! -11  : rank 2 : ndcycle                             !                              !
!     !  12  : rank 2 : s-layer                             !                              !
!     !  120 : rank 2 : s-layer, integer                    !                              !
!     ! -12  : rank 3 : ndcycle,s-layer                     !                              !
!     !  13  : rank 2 : w-layer                             !                              !
!     !  14  : rank 2 : pft                                 !                              !
!     !  146 : rank 3 : npft,ndbh                           !                              !
!     !  15  : rank 2 : disturbance                         !                              !
!     !  16  : rank 2 : dbh                                 !                              !
!     !  17  : rank 2 : age                                 !                              !
!     !  155 : rank 2 : max_lu_years                        !                              !
!     !  156 : rank 3 : max_lu_years, num_lu_transitions    !                              !
!     !-----------------------------------------------------!                              !
!                                                                                          !
!     !----- SITE: -----------------------------------------!                              !
!     !  20  : rank 1 : integer                             !                              !
!     !  21  : rank 1 : real                                !                              !
!     ! -21  : rank 2 : ndcycle                             !                              !
!     !  22  : rank 2 : s-layer                             !                              !
!     !  220 : rank 2 : s-layer, integer                    !                              !
!     !  23  : rank 2 : w-layer                             !                              !
!     !  24  : rank 2 : pft                                 !                              !
!     !  246 : rank 3 : pft, dbh                            !                              !
!     !  25  : rank 2 : disturbance                         !                              !
!     !  255 : rank 3 : disturbance,disturbance             !                              !
!     !  26  : rank 2 : dbh                                 !                              !
!     !  27  : rank 2 : age                                 !                              !
!     !  28  : rank 2 : mortality                           !                              !
!     !  29  : rank 2 : months (12)                         !                              !
!     !-----------------------------------------------------!                              !
!                                                                                          !
!     !----- PATCH: ----------------------------------------!                              !
!     !  30  : rank 1 : integer                             !                              !
!     !  31  : rank 1 : real                                !                              !
!     ! -31  : rank 2 : ndcycle                             !                              !
!     !  32  : rank 2 : s-layer                             !                              !
!     !  320 : rank 2 : s-layer, integer                    !                              !
!     !  33  : rank 2 : w-layer                             !                              !
!     !  34  : rank 2 : pft                                 !                              !
!     !  346 : rank 3 : pft,ff_hgt                          !                              !
!     !  35  : rank 2 : disturbance                         !                              !
!     !  36  : rank 2 : dbh                                 !                              !
!     !  37  : rank 2 : age                                 !                              !
!     !-----------------------------------------------------!                              !
!                                                                                          !
!     !----- COHORT: ---------------------------------------!                              !
!     !  40  : rank 1 : integer                             !                              !
!     !  41  : rank 1 : cohort (real)                       !                              !
!     ! -41  : rank 2 : ndcycle                             !                              !
!     !  44  : rank 2 : cohort, pft                         !                              !
!     !  46  : rank 2 : cohort, dbh                         !                              !
!     !  47  : rank 2 : cohort, age                         !                              !
!     !  49  : rank 2 : cohort, nmonths+1                   !                              !
!     !-----------------------------------------------------!                              !
!                                                                                          !
!     !----- OTHER: ----------------------------------------!                              !
!     ! 90  : rank 0 : integer scalar                       !                              !
!     ! 91  : rank 0 : real scalar                          !                              !
!     ! 92  : rank 1 : real, s-layer                        !                              !
!     ! 96  : rank 1 : real, height class                   !                              !
!     !-----------------------------------------------------!                              !
!                                                                                          !
!------------------------------------------------------------------------------------------!
module ed_var_tables
   use ed_max_dims, only : str_len
   !---------------------------------------------------------------------------------------!
   !    Define data type for main variable table                                           !
   !---------------------------------------------------------------------------------------!
   integer, parameter :: maxvars = 1500
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   type var_table
      logical             :: first
      integer             :: idim_type
      integer             :: nptrs
      integer             :: ihist
      integer             :: ianal
      integer             :: imean
      integer             :: ilite
      integer             :: impti
      integer             :: impt1
      integer             :: impt2
      integer             :: impt3
      integer             :: irecycle
      integer             :: iyear
      integer             :: iopti
      integer             :: imont
      integer             :: idcyc
      integer             :: idail
      integer             :: var_len_global
      character (len=64)  :: name
      character (len=2)   :: dtype
      character (len=64)  :: lname   ! Long name for description in file
      character (len=16)  :: units   ! Unit description of the data
      character (len=64)  :: dimlab
      !----- Multiple pointer defs (maxptrs) ----------------------------------------------!
      type(var_table_vector),pointer,dimension(:) :: vt_vector
   end type var_table
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   type var_table_vector
      real                   , dimension(:), pointer :: var_rp
      integer                , dimension(:), pointer :: var_ip
      character (len=str_len), dimension(:), pointer :: var_cp
      real(kind=8)           , dimension(:), pointer :: var_dp
      real                                 , pointer :: sca_rp
      integer                              , pointer :: sca_ip
      character (len=str_len)              , pointer :: sca_cp
      real(kind=8)                         , pointer :: sca_dp
      integer                                        :: globid
      integer                                        :: varlen
   end type var_table_vector
   !---------------------------------------------------------------------------------------!



   !----- Main variable table allocated to (maxvars,maxgrds) ------------------------------!
   type(var_table), dimension(:,:), allocatable :: vt_info
   !---------------------------------------------------------------------------------------!


   !----- Number of variables for each grid, allocated to "ngrids". -----------------------!
   integer, dimension(:), allocatable :: num_var
   !---------------------------------------------------------------------------------------!



   !=======================================================================================!
   !=======================================================================================!


   contains



   !=======================================================================================!
   !=======================================================================================!
   !   This sub-routine creates the pointers for output, when the variable is real.        !
   !                                                                                       !
   !     NPTS           - Size of this variable                                            !
   !     VAR            - The pointer of the current state variable                        !
   !     NV             - The variable type number                                         !
   !     IGR            - The number of the current grid                                   !
   !     INIT           - Initialize the vt_info?                                          !
   !     GLOB_ID        - The global index of the data                                     !
   !     VAR_LEN        - The length of the states current vector                          !
   !     VAR_LEN_GLOBAL - The length of the entire dataset's vector                        !
   !     MAX_PTRS       - The maximum possible number of pointers necessary for this       !
   !                      variable                                                         !
   !     TABSTR         - The string describing the variables usage                        !
   !---------------------------------------------------------------------------------------!
   recursive subroutine vtable_edio_r(npts,var,nv,igr,init,glob_id,var_len,var_len_global  &
                                     ,max_ptrs,tabstr)
      use ed_max_dims, only : str_len ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      integer                                , intent(in) :: npts
      real(kind=4)          , dimension(npts), target     :: var
      integer                                , intent(in) :: nv
      integer                                , intent(in) :: igr
      integer                                , intent(in) :: init
      integer                                , intent(in) :: glob_id
      integer                                , intent(in) :: var_len
      integer                                , intent(in) :: var_len_global
      integer                                , intent(in) :: max_ptrs
      character (len=*)                      , intent(in) :: tabstr
      !----- Local variables. -------------------------------------------------------------!
      integer                                             :: iptr
      character(len=str_len), dimension(10)               :: tokens
      character(len=8)                                    :: ctab
      integer                                             :: ntok
      integer                                             :: nt
      !----- Local constants. -------------------------------------------------------------!
      character (len=1)                      , parameter  :: toksep=':'
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Determine if this is the first time we view this variable.  If so, then fill  !
      ! some descriptors for the vtable and allocate some space for any pointers that may  !
      ! follow.                                                                            !
      !------------------------------------------------------------------------------------!
      if (init == 0) then

         !----- Count the number of variables. --------------------------------------------!
         num_var(igr) = num_var(igr) + 1
         call tokenize1(tabstr,tokens,ntok,toksep)

         vt_info(nv,igr)%name           = tokens(1)
         vt_info(nv,igr)%dtype          = 'R'  ! This is a real variable (vector)
         vt_info(nv,igr)%nptrs          = 0
         vt_info(nv,igr)%var_len_global = var_len_global

         nullify(vt_info(nv,igr)%vt_vector)
         allocate(vt_info(nv,igr)%vt_vector(max_ptrs))
         read(tokens(2),fmt=*) vt_info(nv,igr)%idim_type

         vt_info(nv,igr)%ihist    = 0
         vt_info(nv,igr)%ianal    = 0
         vt_info(nv,igr)%imean    = 0
         vt_info(nv,igr)%ilite    = 0
         vt_info(nv,igr)%impti    = 0
         vt_info(nv,igr)%impt1    = 0
         vt_info(nv,igr)%impt2    = 0
         vt_info(nv,igr)%impt3    = 0
         vt_info(nv,igr)%irecycle = 0
         vt_info(nv,igr)%imont    = 0
         vt_info(nv,igr)%idail    = 0
         vt_info(nv,igr)%iyear    = 0
         vt_info(nv,igr)%idcyc    = 0
         vt_info(nv,igr)%iopti    = 0
         
         do nt=3,ntok
            ctab=tokens(nt)

            select case (trim(ctab))
            case('hist') 
               vt_info(nv,igr)%ihist    = 1

            case('anal') 
               vt_info(nv,igr)%ianal    = 1

            case('lite') 
               vt_info(nv,igr)%ilite    = 1

            case('mpti') 
               vt_info(nv,igr)%impti    = 1

            case('mpt1') 
               vt_info(nv,igr)%impt1    = 1

            case('mpt2') 
               vt_info(nv,igr)%impt2    = 1

            case('mpt3') 
               vt_info(nv,igr)%impt3    = 1

            case('recycle') 
               vt_info(nv,igr)%irecycle = 1

            case('mont') 
               vt_info(nv,igr)%imont    = 1

            case('dail') 
               vt_info(nv,igr)%idail    = 1

            case('dcyc') 
               vt_info(nv,igr)%idcyc    = 1

            case('year') 
               vt_info(nv,igr)%iyear    = 1

            case('opti') 
               vt_info(nv,igr)%iopti    = 1

            case default
               print*, 'Illegal table specification for var:', tokens(1),ctab
               call fatal_error('Bad var table','vtable_edio_r','ed_var_tables.f90')
            end select
         end do
      else
         !---------------------------------------------------------------------------------!
         !     Make sure that vt_info is associated. If not, call the function with        !
         ! init = 0 then do this part.  Since I think this should never happen, I will     !
         ! also make a fuss to warn the user.                                              !
         !---------------------------------------------------------------------------------!
         if (.not.associated(vt_info(nv,igr)%vt_vector)) then
            write (unit=*,fmt='(a)') ' '
            write (unit=*,fmt='(a)') '----------------------------------------------'
            write (unit=*,fmt='(a)') ' WARNING! WARNING! WARNING! WARNING! WARNING! '
            write (unit=*,fmt='(a)') ' WARNING! WARNING! WARNING! WARNING! WARNING! '
            write (unit=*,fmt='(a)') ' WARNING! WARNING! WARNING! WARNING! WARNING! '
            write (unit=*,fmt='(a)') '----------------------------------------------'
            write (unit=*,fmt='(a)') ' - Subroutine vtable_edio_r (file ed_var_tables.f90)'
            write (unit=*,fmt='(a,1x,i4,1x,a,1x,i2,1x,a)')                                 &
                                     ' - Vt_vector for variable',nv,'of grid',igr          &
                                     ,'is not associated              !'
            write (unit=*,fmt='(a)') ' - I will allocate it now.'
            write (unit=*,fmt='(a,1x,i20,1x,a)') ' - MAX_PTRS=',max_ptrs,'...'
            write (unit=*,fmt='(a,1x,a,1x,a)') ' - Tabstr=',tabstr,'...'
            write (unit=*,fmt='(a)') '----------------------------------------------'
            write (unit=*,fmt='(a)') ' '
            call vtable_edio_r(npts,var,nv,igr,0,glob_id,var_len,var_len_global,max_ptrs   &
                              ,tabstr)
         end if
         
         vt_info(nv,igr)%nptrs                  =  vt_info(nv,igr)%nptrs + 1
         iptr                                   =  vt_info(nv,igr)%nptrs
         vt_info(nv,igr)%vt_vector(iptr)%globid =  glob_id

         vt_info(nv,igr)%vt_vector(iptr)%var_rp => var
         vt_info(nv,igr)%vt_vector(iptr)%varlen =  var_len
      end if
      return
   end subroutine vtable_edio_r
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This sub-routine creates the pointers for output, when the variable is double     !
   ! precision.                                                                            !
   !                                                                                       !
   !     NPTS           - Size of this variable                                            !
   !     VAR            - The pointer of the current state variable                        !
   !     NV             - The variable type number                                         !
   !     IGR            - The number of the current grid                                   !
   !     INIT           - Initialize the vt_info?                                          !
   !     GLOB_ID        - The global index of the data                                     !
   !     VAR_LEN        - The length of the states current vector                          !
   !     VAR_LEN_GLOBAL - The length of the entire dataset's vector                        !
   !     MAX_PTRS       - The maximum possible number of pointers necessary for this       !
   !                      variable                                                         !
   !     TABSTR         - The string describing the variables usage                        !
   !---------------------------------------------------------------------------------------!
   recursive subroutine vtable_edio_d(npts,var,nv,igr,init,glob_id,var_len,var_len_global  &
                                     ,max_ptrs,tabstr)
      use ed_max_dims, only : str_len ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      integer                                , intent(in) :: npts
      real(kind=8)          , dimension(npts), target     :: var
      integer                                , intent(in) :: nv
      integer                                , intent(in) :: igr
      integer                                , intent(in) :: init
      integer                                , intent(in) :: glob_id
      integer                                , intent(in) :: var_len
      integer                                , intent(in) :: var_len_global
      integer                                , intent(in) :: max_ptrs
      character (len=*)                      , intent(in) :: tabstr
      !----- Local variables. -------------------------------------------------------------!
      integer                                             :: iptr
      character(len=str_len), dimension(10)               :: tokens
      character(len=8)                                    :: ctab
      integer                                             :: ntok
      integer                                             :: nt
      !----- Local constants. -------------------------------------------------------------!
      character (len=1)                      , parameter  :: toksep=':'
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Determine if this is the first time we view this variable.  If so, then fill  !
      ! some descriptors for the vtable and allocate some space for any pointers that may  !
      ! follow.                                                                            !
      !------------------------------------------------------------------------------------!
      if (init == 0) then

         !----- Count the number of variables. --------------------------------------------!
         num_var(igr) = num_var(igr) + 1
         call tokenize1(tabstr,tokens,ntok,toksep)

         vt_info(nv,igr)%name           = tokens(1)
         vt_info(nv,igr)%dtype          = 'D'  ! Double precision variable (vector)
         vt_info(nv,igr)%nptrs          = 0
         vt_info(nv,igr)%var_len_global = var_len_global

         nullify(vt_info(nv,igr)%vt_vector)
         allocate(vt_info(nv,igr)%vt_vector(max_ptrs))
         read(tokens(2),fmt=*) vt_info(nv,igr)%idim_type

         vt_info(nv,igr)%ihist    = 0
         vt_info(nv,igr)%ianal    = 0
         vt_info(nv,igr)%imean    = 0
         vt_info(nv,igr)%ilite    = 0
         vt_info(nv,igr)%impti    = 0
         vt_info(nv,igr)%impt1    = 0
         vt_info(nv,igr)%impt2    = 0
         vt_info(nv,igr)%impt3    = 0
         vt_info(nv,igr)%irecycle = 0
         vt_info(nv,igr)%imont    = 0
         vt_info(nv,igr)%idail    = 0
         vt_info(nv,igr)%idcyc    = 0
         vt_info(nv,igr)%iyear    = 0
         vt_info(nv,igr)%iopti    = 0
         
         do nt=3,ntok
            ctab=tokens(nt)

            select case (trim(ctab))
            case('hist') 
               vt_info(nv,igr)%ihist    = 1

            case('anal') 
               vt_info(nv,igr)%ianal    = 1

            case('lite') 
               vt_info(nv,igr)%ilite    = 1

            case('mpti') 
               vt_info(nv,igr)%impti    = 1

            case('mpt1') 
               vt_info(nv,igr)%impt1    = 1

            case('mpt2') 
               vt_info(nv,igr)%impt2    = 1

            case('mpt3') 
               vt_info(nv,igr)%impt3    = 1

            case('recycle') 
               vt_info(nv,igr)%irecycle = 1

            case('mont') 
               vt_info(nv,igr)%imont    = 1

            case('dail') 
               vt_info(nv,igr)%idail    = 1

            case('dcyc') 
               vt_info(nv,igr)%idcyc    = 1

            case('year') 
               vt_info(nv,igr)%iyear    = 1

            case('opti') 
               vt_info(nv,igr)%iopti    = 1

            case default
               print*, 'Illegal table specification for var:', tokens(1),ctab
               call fatal_error('Bad var table','vtable_edio_r','ed_var_tables.f90')
            end select
         end do
      else
         !---------------------------------------------------------------------------------!
         !     Make sure that vt_info is associated. If not, call the function with        !
         ! init = 0 then do this part.  Since I think this should never happen, I will     !
         ! also make a fuss to warn the user.                                              !
         !---------------------------------------------------------------------------------!
         if (.not.associated(vt_info(nv,igr)%vt_vector)) then
            write (unit=*,fmt='(a)') ' '
            write (unit=*,fmt='(a)') '----------------------------------------------'
            write (unit=*,fmt='(a)') ' WARNING! WARNING! WARNING! WARNING! WARNING! '
            write (unit=*,fmt='(a)') ' WARNING! WARNING! WARNING! WARNING! WARNING! '
            write (unit=*,fmt='(a)') ' WARNING! WARNING! WARNING! WARNING! WARNING! '
            write (unit=*,fmt='(a)') '----------------------------------------------'
            write (unit=*,fmt='(a)') ' - Subroutine vtable_edio_r (file ed_var_tables.f90)'
            write (unit=*,fmt='(a,1x,i4,1x,a,1x,i2,1x,a)')                                 &
                                     ' - Vt_vector for variable',nv,'of grid',igr          &
                                     ,'is not associated              !'
            write (unit=*,fmt='(a)') ' - I will allocate it now.'
            write (unit=*,fmt='(a,1x,i20,1x,a)') ' - MAX_PTRS=',max_ptrs,'...'
            write (unit=*,fmt='(a,1x,a,1x,a)') ' - Tabstr=',tabstr,'...'
            write (unit=*,fmt='(a)') '----------------------------------------------'
            write (unit=*,fmt='(a)') ' '
            call vtable_edio_d(npts,var,nv,igr,0,glob_id,var_len,var_len_global,max_ptrs   &
                              ,tabstr)
         end if
         
         vt_info(nv,igr)%nptrs                  =  vt_info(nv,igr)%nptrs + 1
         iptr                                   =  vt_info(nv,igr)%nptrs
         vt_info(nv,igr)%vt_vector(iptr)%globid =  glob_id

         vt_info(nv,igr)%vt_vector(iptr)%var_dp => var
         vt_info(nv,igr)%vt_vector(iptr)%varlen =  var_len
      end if
      return
   end subroutine vtable_edio_d
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This sub-routine creates the pointers for output, when the variable is integer.   !
   !                                                                                       !
   !     NPTS           - Size of this variable                                            !
   !     VAR            - The pointer of the current state variable                        !
   !     NV             - The variable type number                                         !
   !     IGR            - The number of the current grid                                   !
   !     INIT           - Initialize the vt_info?                                          !
   !     GLOB_ID        - The global index of the data                                     !
   !     VAR_LEN        - The length of the states current vector                          !
   !     VAR_LEN_GLOBAL - The length of the entire dataset's vector                        !
   !     MAX_PTRS       - The maximum possible number of pointers necessary for this       !
   !                      variable                                                         !
   !     TABSTR         - The string describing the variables usage                        !
   !---------------------------------------------------------------------------------------!
   recursive subroutine vtable_edio_i(npts,var,nv,igr,init,glob_id,var_len,var_len_global  &
                                     ,max_ptrs,tabstr)
      use ed_max_dims, only : str_len ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      integer                                , intent(in) :: npts
      integer               , dimension(npts), target     :: var
      integer                                , intent(in) :: nv
      integer                                , intent(in) :: igr
      integer                                , intent(in) :: init
      integer                                , intent(in) :: glob_id
      integer                                , intent(in) :: var_len
      integer                                , intent(in) :: var_len_global
      integer                                , intent(in) :: max_ptrs
      character (len=*)                      , intent(in) :: tabstr
      !----- Local variables. -------------------------------------------------------------!
      integer                                             :: iptr
      character(len=str_len), dimension(10)               :: tokens
      character(len=8)                                    :: ctab
      integer                                             :: ntok
      integer                                             :: nt
      !----- Local constants. -------------------------------------------------------------!
      character (len=1)                      , parameter  :: toksep=':'
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Determine if this is the first time we view this variable.  If so, then fill  !
      ! some descriptors for the vtable and allocate some space for any pointers that may  !
      ! follow.                                                                            !
      !------------------------------------------------------------------------------------!
      if (init == 0) then

         !----- Count the number of variables. --------------------------------------------!
         num_var(igr) = num_var(igr) + 1
         call tokenize1(tabstr,tokens,ntok,toksep)

         vt_info(nv,igr)%name           = tokens(1)
         vt_info(nv,igr)%dtype          = 'I'  ! Integer variable (vector)
         vt_info(nv,igr)%nptrs          = 0
         vt_info(nv,igr)%var_len_global = var_len_global

         nullify(vt_info(nv,igr)%vt_vector)
         allocate(vt_info(nv,igr)%vt_vector(max_ptrs))
         read(tokens(2),fmt=*) vt_info(nv,igr)%idim_type

         vt_info(nv,igr)%ihist    = 0
         vt_info(nv,igr)%ianal    = 0
         vt_info(nv,igr)%imean    = 0
         vt_info(nv,igr)%ilite    = 0
         vt_info(nv,igr)%impti    = 0
         vt_info(nv,igr)%impt1    = 0
         vt_info(nv,igr)%impt2    = 0
         vt_info(nv,igr)%impt3    = 0
         vt_info(nv,igr)%irecycle = 0
         vt_info(nv,igr)%imont    = 0
         vt_info(nv,igr)%idail    = 0
         vt_info(nv,igr)%idcyc    = 0
         vt_info(nv,igr)%iyear    = 0
         vt_info(nv,igr)%iopti    = 0
         
         do nt=3,ntok
            ctab=tokens(nt)

            select case (trim(ctab))
            case('hist') 
               vt_info(nv,igr)%ihist    = 1

            case('anal') 
               vt_info(nv,igr)%ianal    = 1

            case('lite') 
               vt_info(nv,igr)%ilite    = 1

            case('mpti') 
               vt_info(nv,igr)%impti    = 1

            case('mpt1') 
               vt_info(nv,igr)%impt1    = 1

            case('mpt2') 
               vt_info(nv,igr)%impt2    = 1

            case('mpt3') 
               vt_info(nv,igr)%impt3    = 1

            case('recycle') 
               vt_info(nv,igr)%irecycle = 1

            case('mont') 
               vt_info(nv,igr)%imont    = 1

            case('dail') 
               vt_info(nv,igr)%idail    = 1

            case('dcyc') 
               vt_info(nv,igr)%idcyc    = 1

            case('year') 
               vt_info(nv,igr)%iyear    = 1

            case('opti') 
               vt_info(nv,igr)%iopti    = 1

            case default
               print*, 'Illegal table specification for var:', tokens(1),ctab
               call fatal_error('Bad var table','vtable_edio_r','ed_var_tables.f90')
            end select
         end do
      else
         !---------------------------------------------------------------------------------!
         !     Make sure that vt_info is associated. If not, call the function with        !
         ! init = 0 then do this part.  Since I think this should never happen, I will     !
         ! also make a fuss to warn the user.                                              !
         !---------------------------------------------------------------------------------!
         if (.not.associated(vt_info(nv,igr)%vt_vector)) then
            write (unit=*,fmt='(a)') ' '
            write (unit=*,fmt='(a)') '----------------------------------------------'
            write (unit=*,fmt='(a)') ' WARNING! WARNING! WARNING! WARNING! WARNING! '
            write (unit=*,fmt='(a)') ' WARNING! WARNING! WARNING! WARNING! WARNING! '
            write (unit=*,fmt='(a)') ' WARNING! WARNING! WARNING! WARNING! WARNING! '
            write (unit=*,fmt='(a)') '----------------------------------------------'
            write (unit=*,fmt='(a)') ' - Subroutine vtable_edio_r (file ed_var_tables.f90)'
            write (unit=*,fmt='(a,1x,i4,1x,a,1x,i2,1x,a)')                                 &
                                     ' - Vt_vector for variable',nv,'of grid',igr          &
                                     ,'is not associated              !'
            write (unit=*,fmt='(a)') ' - I will allocate it now.'
            write (unit=*,fmt='(a,1x,i20,1x,a)') ' - MAX_PTRS=',max_ptrs,'...'
            write (unit=*,fmt='(a,1x,a,1x,a)') ' - Tabstr=',tabstr,'...'
            write (unit=*,fmt='(a)') '----------------------------------------------'
            write (unit=*,fmt='(a)') ' '
            call vtable_edio_i(npts,var,nv,igr,0,glob_id,var_len,var_len_global,max_ptrs   &
                              ,tabstr)
         end if
         
         vt_info(nv,igr)%nptrs                  =  vt_info(nv,igr)%nptrs + 1
         iptr                                   =  vt_info(nv,igr)%nptrs
         vt_info(nv,igr)%vt_vector(iptr)%globid =  glob_id

         vt_info(nv,igr)%vt_vector(iptr)%var_ip => var
         vt_info(nv,igr)%vt_vector(iptr)%varlen =  var_len
      end if
      return
   end subroutine vtable_edio_i
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This sub-routine creates the pointers for output, when the variable is character. !
   !                                                                                       !
   !     NPTS           - Size of this variable                                            !
   !     VAR            - The pointer of the current state variable                        !
   !     NV             - The variable type number                                         !
   !     IGR            - The number of the current grid                                   !
   !     INIT           - Initialize the vt_info?                                          !
   !     GLOB_ID        - The global index of the data                                     !
   !     VAR_LEN        - The length of the states current vector                          !
   !     VAR_LEN_GLOBAL - The length of the entire dataset's vector                        !
   !     MAX_PTRS       - The maximum possible number of pointers necessary for this       !
   !                      variable                                                         !
   !     TABSTR         - The string describing the variables usage                        !
   !---------------------------------------------------------------------------------------!
   recursive subroutine vtable_edio_c(npts,var,nv,igr,init,glob_id,var_len,var_len_global  &
                                     ,max_ptrs,tabstr)
      use ed_max_dims, only : str_len ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      integer                                , intent(in) :: npts
      character(len=str_len), dimension(npts), target     :: var
      integer                                , intent(in) :: nv
      integer                                , intent(in) :: igr
      integer                                , intent(in) :: init
      integer                                , intent(in) :: glob_id
      integer                                , intent(in) :: var_len
      integer                                , intent(in) :: var_len_global
      integer                                , intent(in) :: max_ptrs
      character (len=*)                      , intent(in) :: tabstr
      !----- Local variables. -------------------------------------------------------------!
      integer                                             :: iptr
      character(len=str_len), dimension(10)               :: tokens
      character(len=8)                                    :: ctab
      integer                                             :: ntok
      integer                                             :: nt
      !----- Local constants. -------------------------------------------------------------!
      character (len=1)                      , parameter  :: toksep=':'
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Determine if this is the first time we view this variable.  If so, then fill  !
      ! some descriptors for the vtable and allocate some space for any pointers that may  !
      ! follow.                                                                            !
      !------------------------------------------------------------------------------------!
      if (init == 0) then

         !----- Count the number of variables. --------------------------------------------!
         num_var(igr) = num_var(igr) + 1
         call tokenize1(tabstr,tokens,ntok,toksep)

         vt_info(nv,igr)%name           = tokens(1)
         vt_info(nv,igr)%dtype          = 'C'  ! Character variable (vector)
         vt_info(nv,igr)%nptrs          = 0
         vt_info(nv,igr)%var_len_global = var_len_global

         nullify(vt_info(nv,igr)%vt_vector)
         allocate(vt_info(nv,igr)%vt_vector(max_ptrs))
         read(tokens(2),fmt=*) vt_info(nv,igr)%idim_type

         vt_info(nv,igr)%ihist    = 0
         vt_info(nv,igr)%ianal    = 0
         vt_info(nv,igr)%imean    = 0
         vt_info(nv,igr)%ilite    = 0
         vt_info(nv,igr)%impti    = 0
         vt_info(nv,igr)%impt1    = 0
         vt_info(nv,igr)%impt2    = 0
         vt_info(nv,igr)%impt3    = 0
         vt_info(nv,igr)%irecycle = 0
         vt_info(nv,igr)%imont    = 0
         vt_info(nv,igr)%idail    = 0
         vt_info(nv,igr)%idcyc    = 0
         vt_info(nv,igr)%iyear    = 0
         vt_info(nv,igr)%iopti    = 0
         
         do nt=3,ntok
            ctab=tokens(nt)

            select case (trim(ctab))
            case('hist') 
               vt_info(nv,igr)%ihist    = 1

            case('anal') 
               vt_info(nv,igr)%ianal    = 1

            case('lite') 
               vt_info(nv,igr)%ilite    = 1

            case('mpti') 
               vt_info(nv,igr)%impti    = 1

            case('mpt1') 
               vt_info(nv,igr)%impt1    = 1

            case('mpt2') 
               vt_info(nv,igr)%impt2    = 1

            case('mpt3') 
               vt_info(nv,igr)%impt3    = 1

            case('recycle') 
               vt_info(nv,igr)%irecycle = 1

            case('mont') 
               vt_info(nv,igr)%imont    = 1

            case('dail') 
               vt_info(nv,igr)%idail    = 1

            case('dcyc') 
               vt_info(nv,igr)%idcyc    = 1

            case('year') 
               vt_info(nv,igr)%iyear    = 1

            case('opti') 
               vt_info(nv,igr)%iopti    = 1

            case default
               print*, 'Illegal table specification for var:', tokens(1),ctab
               call fatal_error('Bad var table','vtable_edio_r','ed_var_tables.f90')
            end select
         end do
      else
         !---------------------------------------------------------------------------------!
         !     Make sure that vt_info is associated. If not, call the function with        !
         ! init = 0 then do this part.  Since I think this should never happen, I will     !
         ! also make a fuss to warn the user.                                              !
         !---------------------------------------------------------------------------------!
         if (.not.associated(vt_info(nv,igr)%vt_vector)) then
            write (unit=*,fmt='(a)') ' '
            write (unit=*,fmt='(a)') '----------------------------------------------'
            write (unit=*,fmt='(a)') ' WARNING! WARNING! WARNING! WARNING! WARNING! '
            write (unit=*,fmt='(a)') ' WARNING! WARNING! WARNING! WARNING! WARNING! '
            write (unit=*,fmt='(a)') ' WARNING! WARNING! WARNING! WARNING! WARNING! '
            write (unit=*,fmt='(a)') '----------------------------------------------'
            write (unit=*,fmt='(a)') ' - Subroutine vtable_edio_r (file ed_var_tables.f90)'
            write (unit=*,fmt='(a,1x,i4,1x,a,1x,i2,1x,a)')                                 &
                                     ' - Vt_vector for variable',nv,'of grid',igr          &
                                     ,'is not associated              !'
            write (unit=*,fmt='(a)') ' - I will allocate it now.'
            write (unit=*,fmt='(a,1x,i20,1x,a)') ' - MAX_PTRS=',max_ptrs,'...'
            write (unit=*,fmt='(a,1x,a,1x,a)') ' - Tabstr=',tabstr,'...'
            write (unit=*,fmt='(a)') '----------------------------------------------'
            write (unit=*,fmt='(a)') ' '
            call vtable_edio_c(npts,var,nv,igr,0,glob_id,var_len,var_len_global,max_ptrs   &
                              ,tabstr)
         end if
         
         vt_info(nv,igr)%nptrs                  =  vt_info(nv,igr)%nptrs + 1
         iptr                                   =  vt_info(nv,igr)%nptrs
         vt_info(nv,igr)%vt_vector(iptr)%globid =  glob_id

         vt_info(nv,igr)%vt_vector(iptr)%var_cp => var
         vt_info(nv,igr)%vt_vector(iptr)%varlen =  var_len
      end if
      return
   end subroutine vtable_edio_c
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !   This sub-routine creates the pointers for output, when the variable is real.        !
   !                                                                                       !
   !     VAR            - The pointer of the current state variable                        !
   !     NV             - The variable type number                                         !
   !     IGR            - The number of the current grid                                   !
   !     INIT           - Initialize the vt_info?                                          !
   !     GLOB_ID        - The global index of the data                                     !
   !     VAR_LEN        - The length of the states current vector                          !
   !     VAR_LEN_GLOBAL - The length of the entire dataset's vector                        !
   !     MAX_PTRS       - The maximum possible number of pointers necessary for this       !
   !                      variable                                                         !
   !     TABSTR         - The string describing the variables usage                        !
   !---------------------------------------------------------------------------------------!
   recursive subroutine vtable_edio_r_sca(var,nv,igr,init,glob_id,var_len,var_len_global   &
                                         ,max_ptrs,tabstr)
      use ed_max_dims, only : str_len ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real(kind=4)                           , target     :: var
      integer                                , intent(in) :: nv
      integer                                , intent(in) :: igr
      integer                                , intent(in) :: init
      integer                                , intent(in) :: glob_id
      integer                                , intent(in) :: var_len
      integer                                , intent(in) :: var_len_global
      integer                                , intent(in) :: max_ptrs
      character (len=*)                      , intent(in) :: tabstr
      !----- Local variables. -------------------------------------------------------------!
      integer                                             :: iptr
      character(len=str_len), dimension(10)               :: tokens
      character(len=8)                                    :: ctab
      integer                                             :: ntok
      integer                                             :: nt
      !----- Local constants. -------------------------------------------------------------!
      character (len=1)                      , parameter  :: toksep=':'
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Determine if this is the first time we view this variable.  If so, then fill  !
      ! some descriptors for the vtable and allocate some space for any pointers that may  !
      ! follow.                                                                            !
      !------------------------------------------------------------------------------------!
      if (init == 0) then

         !----- Count the number of variables. --------------------------------------------!
         num_var(igr) = num_var(igr) + 1
         call tokenize1(tabstr,tokens,ntok,toksep)

         vt_info(nv,igr)%name           = tokens(1)
         vt_info(nv,igr)%dtype          = 'r'  ! Real variable (scalar)
         vt_info(nv,igr)%nptrs          = 0
         vt_info(nv,igr)%var_len_global = var_len_global

         nullify(vt_info(nv,igr)%vt_vector)
         allocate(vt_info(nv,igr)%vt_vector(max_ptrs))
         read(tokens(2),fmt=*) vt_info(nv,igr)%idim_type

         vt_info(nv,igr)%ihist    = 0
         vt_info(nv,igr)%ianal    = 0
         vt_info(nv,igr)%imean    = 0
         vt_info(nv,igr)%ilite    = 0
         vt_info(nv,igr)%impti    = 0
         vt_info(nv,igr)%impt1    = 0
         vt_info(nv,igr)%impt2    = 0
         vt_info(nv,igr)%impt3    = 0
         vt_info(nv,igr)%irecycle = 0
         vt_info(nv,igr)%imont    = 0
         vt_info(nv,igr)%idail    = 0
         vt_info(nv,igr)%idcyc    = 0
         vt_info(nv,igr)%iyear    = 0
         vt_info(nv,igr)%iopti    = 0
         
         do nt=3,ntok
            ctab=tokens(nt)

            select case (trim(ctab))
            case('hist') 
               vt_info(nv,igr)%ihist    = 1

            case('anal') 
               vt_info(nv,igr)%ianal    = 1

            case('lite') 
               vt_info(nv,igr)%ilite    = 1

            case('mpti') 
               vt_info(nv,igr)%impti    = 1

            case('mpt1') 
               vt_info(nv,igr)%impt1    = 1

            case('mpt2') 
               vt_info(nv,igr)%impt2    = 1

            case('mpt3') 
               vt_info(nv,igr)%impt3    = 1

            case('recycle') 
               vt_info(nv,igr)%irecycle = 1

            case('mont') 
               vt_info(nv,igr)%imont    = 1

            case('dail') 
               vt_info(nv,igr)%idail    = 1

            case('dcyc') 
               vt_info(nv,igr)%idcyc    = 1

            case('year') 
               vt_info(nv,igr)%iyear    = 1

            case('opti') 
               vt_info(nv,igr)%iopti    = 1

            case default
               print*, 'Illegal table specification for var:', tokens(1),ctab
               call fatal_error('Bad var table','vtable_edio_r','ed_var_tables.f90')
            end select
         end do
      else
         !---------------------------------------------------------------------------------!
         !     Make sure that vt_info is associated. If not, call the function with        !
         ! init = 0 then do this part.  Since I think this should never happen, I will     !
         ! also make a fuss to warn the user.                                              !
         !---------------------------------------------------------------------------------!
         if (.not.associated(vt_info(nv,igr)%vt_vector)) then
            write (unit=*,fmt='(a)') ' '
            write (unit=*,fmt='(a)') '----------------------------------------------'
            write (unit=*,fmt='(a)') ' WARNING! WARNING! WARNING! WARNING! WARNING! '
            write (unit=*,fmt='(a)') ' WARNING! WARNING! WARNING! WARNING! WARNING! '
            write (unit=*,fmt='(a)') ' WARNING! WARNING! WARNING! WARNING! WARNING! '
            write (unit=*,fmt='(a)') '----------------------------------------------'
            write (unit=*,fmt='(a)') ' - Subroutine vtable_edio_r (file ed_var_tables.f90)'
            write (unit=*,fmt='(a,1x,i4,1x,a,1x,i2,1x,a)')                                 &
                                     ' - Vt_vector for variable',nv,'of grid',igr          &
                                     ,'is not associated              !'
            write (unit=*,fmt='(a)') ' - I will allocate it now.'
            write (unit=*,fmt='(a,1x,i20,1x,a)') ' - MAX_PTRS=',max_ptrs,'...'
            write (unit=*,fmt='(a,1x,a,1x,a)') ' - Tabstr=',tabstr,'...'
            write (unit=*,fmt='(a)') '----------------------------------------------'
            write (unit=*,fmt='(a)') ' '
            call vtable_edio_r_sca(var,nv,igr,0,glob_id,var_len,var_len_global,max_ptrs    &
                                  ,tabstr)
         end if
         
         vt_info(nv,igr)%nptrs                  =  vt_info(nv,igr)%nptrs + 1
         iptr                                   =  vt_info(nv,igr)%nptrs
         vt_info(nv,igr)%vt_vector(iptr)%globid =  glob_id

         vt_info(nv,igr)%vt_vector(iptr)%sca_rp => var
         vt_info(nv,igr)%vt_vector(iptr)%varlen =  var_len
      end if
      return
   end subroutine vtable_edio_r_sca
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This sub-routine creates the pointers for output, when the variable is double     !
   ! precision.                                                                            !
   !                                                                                       !
   !     VAR            - The pointer of the current state variable                        !
   !     NV             - The variable type number                                         !
   !     IGR            - The number of the current grid                                   !
   !     INIT           - Initialize the vt_info?                                          !
   !     GLOB_ID        - The global index of the data                                     !
   !     VAR_LEN        - The length of the states current vector                          !
   !     VAR_LEN_GLOBAL - The length of the entire dataset's vector                        !
   !     MAX_PTRS       - The maximum possible number of pointers necessary for this       !
   !                      variable                                                         !
   !     TABSTR         - The string describing the variables usage                        !
   !---------------------------------------------------------------------------------------!
   recursive subroutine vtable_edio_d_sca(var,nv,igr,init,glob_id,var_len,var_len_global   &
                                         ,max_ptrs,tabstr)
      use ed_max_dims, only : str_len ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real(kind=8)                           , target     :: var
      integer                                , intent(in) :: nv
      integer                                , intent(in) :: igr
      integer                                , intent(in) :: init
      integer                                , intent(in) :: glob_id
      integer                                , intent(in) :: var_len
      integer                                , intent(in) :: var_len_global
      integer                                , intent(in) :: max_ptrs
      character (len=*)                      , intent(in) :: tabstr
      !----- Local variables. -------------------------------------------------------------!
      integer                                             :: iptr
      character(len=str_len), dimension(10)               :: tokens
      character(len=8)                                    :: ctab
      integer                                             :: ntok
      integer                                             :: nt
      !----- Local constants. -------------------------------------------------------------!
      character (len=1)                      , parameter  :: toksep=':'
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Determine if this is the first time we view this variable.  If so, then fill  !
      ! some descriptors for the vtable and allocate some space for any pointers that may  !
      ! follow.                                                                            !
      !------------------------------------------------------------------------------------!
      if (init == 0) then

         !----- Count the number of variables. --------------------------------------------!
         num_var(igr) = num_var(igr) + 1
         call tokenize1(tabstr,tokens,ntok,toksep)

         vt_info(nv,igr)%name           = tokens(1)
         vt_info(nv,igr)%dtype          = 'd'  ! Double precision variable (scalar)
         vt_info(nv,igr)%nptrs          = 0
         vt_info(nv,igr)%var_len_global = var_len_global

         nullify(vt_info(nv,igr)%vt_vector)
         allocate(vt_info(nv,igr)%vt_vector(max_ptrs))
         read(tokens(2),fmt=*) vt_info(nv,igr)%idim_type

         vt_info(nv,igr)%ihist    = 0
         vt_info(nv,igr)%ianal    = 0
         vt_info(nv,igr)%imean    = 0
         vt_info(nv,igr)%ilite    = 0
         vt_info(nv,igr)%impti    = 0
         vt_info(nv,igr)%impt1    = 0
         vt_info(nv,igr)%impt2    = 0
         vt_info(nv,igr)%impt3    = 0
         vt_info(nv,igr)%irecycle = 0
         vt_info(nv,igr)%imont    = 0
         vt_info(nv,igr)%idail    = 0
         vt_info(nv,igr)%idcyc    = 0
         vt_info(nv,igr)%iyear    = 0
         vt_info(nv,igr)%iopti    = 0
         
         do nt=3,ntok
            ctab=tokens(nt)

            select case (trim(ctab))
            case('hist') 
               vt_info(nv,igr)%ihist    = 1

            case('anal') 
               vt_info(nv,igr)%ianal    = 1

            case('lite') 
               vt_info(nv,igr)%ilite    = 1

            case('mpti') 
               vt_info(nv,igr)%impti    = 1

            case('mpt1') 
               vt_info(nv,igr)%impt1    = 1

            case('mpt2') 
               vt_info(nv,igr)%impt2    = 1

            case('mpt3') 
               vt_info(nv,igr)%impt3    = 1

            case('recycle') 
               vt_info(nv,igr)%irecycle = 1

            case('mont') 
               vt_info(nv,igr)%imont    = 1

            case('dail') 
               vt_info(nv,igr)%idail    = 1

            case('dcyc') 
               vt_info(nv,igr)%idcyc    = 1

            case('year') 
               vt_info(nv,igr)%iyear    = 1

            case('opti') 
               vt_info(nv,igr)%iopti    = 1

            case default
               print*, 'Illegal table specification for var:', tokens(1),ctab
               call fatal_error('Bad var table','vtable_edio_r','ed_var_tables.f90')
            end select
         end do
      else
         !---------------------------------------------------------------------------------!
         !     Make sure that vt_info is associated. If not, call the function with        !
         ! init = 0 then do this part.  Since I think this should never happen, I will     !
         ! also make a fuss to warn the user.                                              !
         !---------------------------------------------------------------------------------!
         if (.not.associated(vt_info(nv,igr)%vt_vector)) then
            write (unit=*,fmt='(a)') ' '
            write (unit=*,fmt='(a)') '----------------------------------------------'
            write (unit=*,fmt='(a)') ' WARNING! WARNING! WARNING! WARNING! WARNING! '
            write (unit=*,fmt='(a)') ' WARNING! WARNING! WARNING! WARNING! WARNING! '
            write (unit=*,fmt='(a)') ' WARNING! WARNING! WARNING! WARNING! WARNING! '
            write (unit=*,fmt='(a)') '----------------------------------------------'
            write (unit=*,fmt='(a)') ' - Subroutine vtable_edio_r (file ed_var_tables.f90)'
            write (unit=*,fmt='(a,1x,i4,1x,a,1x,i2,1x,a)')                                 &
                                     ' - Vt_vector for variable',nv,'of grid',igr          &
                                     ,'is not associated              !'
            write (unit=*,fmt='(a)') ' - I will allocate it now.'
            write (unit=*,fmt='(a,1x,i20,1x,a)') ' - MAX_PTRS=',max_ptrs,'...'
            write (unit=*,fmt='(a,1x,a,1x,a)') ' - Tabstr=',tabstr,'...'
            write (unit=*,fmt='(a)') '----------------------------------------------'
            write (unit=*,fmt='(a)') ' '
            call vtable_edio_d_sca(var,nv,igr,0,glob_id,var_len,var_len_global,max_ptrs    &
                                  ,tabstr)
         end if
         
         vt_info(nv,igr)%nptrs                  =  vt_info(nv,igr)%nptrs + 1
         iptr                                   =  vt_info(nv,igr)%nptrs
         vt_info(nv,igr)%vt_vector(iptr)%globid =  glob_id

         vt_info(nv,igr)%vt_vector(iptr)%sca_dp => var
         vt_info(nv,igr)%vt_vector(iptr)%varlen =  var_len
      end if
      return
   end subroutine vtable_edio_d_sca
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This sub-routine creates the pointers for output, when the variable is integer.   !
   !                                                                                       !
   !     VAR            - The pointer of the current state variable                        !
   !     NV             - The variable type number                                         !
   !     IGR            - The number of the current grid                                   !
   !     INIT           - Initialize the vt_info?                                          !
   !     GLOB_ID        - The global index of the data                                     !
   !     VAR_LEN        - The length of the states current vector                          !
   !     VAR_LEN_GLOBAL - The length of the entire dataset's vector                        !
   !     MAX_PTRS       - The maximum possible number of pointers necessary for this       !
   !                      variable                                                         !
   !     TABSTR         - The string describing the variables usage                        !
   !---------------------------------------------------------------------------------------!
   recursive subroutine vtable_edio_i_sca(var,nv,igr,init,glob_id,var_len,var_len_global   &
                                         ,max_ptrs,tabstr)
      use ed_max_dims, only : str_len ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      integer                                , target     :: var
      integer                                , intent(in) :: nv
      integer                                , intent(in) :: igr
      integer                                , intent(in) :: init
      integer                                , intent(in) :: glob_id
      integer                                , intent(in) :: var_len
      integer                                , intent(in) :: var_len_global
      integer                                , intent(in) :: max_ptrs
      character (len=*)                      , intent(in) :: tabstr
      !----- Local variables. -------------------------------------------------------------!
      integer                                             :: iptr
      character(len=str_len), dimension(10)               :: tokens
      character(len=8)                                    :: ctab
      integer                                             :: ntok
      integer                                             :: nt
      !----- Local constants. -------------------------------------------------------------!
      character (len=1)                      , parameter  :: toksep=':'
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Determine if this is the first time we view this variable.  If so, then fill  !
      ! some descriptors for the vtable and allocate some space for any pointers that may  !
      ! follow.                                                                            !
      !------------------------------------------------------------------------------------!
      if (init == 0) then

         !----- Count the number of variables. --------------------------------------------!
         num_var(igr) = num_var(igr) + 1
         call tokenize1(tabstr,tokens,ntok,toksep)

         vt_info(nv,igr)%name           = tokens(1)
         vt_info(nv,igr)%dtype          = 'i'  ! Integer variable (scalar)
         vt_info(nv,igr)%nptrs          = 0
         vt_info(nv,igr)%var_len_global = var_len_global

         nullify(vt_info(nv,igr)%vt_vector)
         allocate(vt_info(nv,igr)%vt_vector(max_ptrs))
         read(tokens(2),fmt=*) vt_info(nv,igr)%idim_type

         vt_info(nv,igr)%ihist    = 0
         vt_info(nv,igr)%ianal    = 0
         vt_info(nv,igr)%imean    = 0
         vt_info(nv,igr)%ilite    = 0
         vt_info(nv,igr)%impti    = 0
         vt_info(nv,igr)%impt1    = 0
         vt_info(nv,igr)%impt2    = 0
         vt_info(nv,igr)%impt3    = 0
         vt_info(nv,igr)%irecycle = 0
         vt_info(nv,igr)%imont    = 0
         vt_info(nv,igr)%idail    = 0
         vt_info(nv,igr)%idcyc    = 0
         vt_info(nv,igr)%iyear    = 0
         vt_info(nv,igr)%iopti    = 0
         
         do nt=3,ntok
            ctab=tokens(nt)

            select case (trim(ctab))
            case('hist') 
               vt_info(nv,igr)%ihist    = 1

            case('anal') 
               vt_info(nv,igr)%ianal    = 1

            case('lite') 
               vt_info(nv,igr)%ilite    = 1

            case('mpti') 
               vt_info(nv,igr)%impti    = 1

            case('mpt1') 
               vt_info(nv,igr)%impt1    = 1

            case('mpt2') 
               vt_info(nv,igr)%impt2    = 1

            case('mpt3') 
               vt_info(nv,igr)%impt3    = 1

            case('recycle') 
               vt_info(nv,igr)%irecycle = 1

            case('mont') 
               vt_info(nv,igr)%imont    = 1

            case('dail') 
               vt_info(nv,igr)%idail    = 1

            case('dcyc') 
               vt_info(nv,igr)%idcyc    = 1

            case('year') 
               vt_info(nv,igr)%iyear    = 1

            case('opti') 
               vt_info(nv,igr)%iopti    = 1

            case default
               print*, 'Illegal table specification for var:', tokens(1),ctab
               call fatal_error('Bad var table','vtable_edio_r','ed_var_tables.f90')
            end select
         end do
      else
         !---------------------------------------------------------------------------------!
         !     Make sure that vt_info is associated. If not, call the function with        !
         ! init = 0 then do this part.  Since I think this should never happen, I will     !
         ! also make a fuss to warn the user.                                              !
         !---------------------------------------------------------------------------------!
         if (.not.associated(vt_info(nv,igr)%vt_vector)) then
            write (unit=*,fmt='(a)') ' '
            write (unit=*,fmt='(a)') '----------------------------------------------'
            write (unit=*,fmt='(a)') ' WARNING! WARNING! WARNING! WARNING! WARNING! '
            write (unit=*,fmt='(a)') ' WARNING! WARNING! WARNING! WARNING! WARNING! '
            write (unit=*,fmt='(a)') ' WARNING! WARNING! WARNING! WARNING! WARNING! '
            write (unit=*,fmt='(a)') '----------------------------------------------'
            write (unit=*,fmt='(a)') ' - Subroutine vtable_edio_r (file ed_var_tables.f90)'
            write (unit=*,fmt='(a,1x,i4,1x,a,1x,i2,1x,a)')                                 &
                                     ' - Vt_vector for variable',nv,'of grid',igr          &
                                     ,'is not associated              !'
            write (unit=*,fmt='(a)') ' - I will allocate it now.'
            write (unit=*,fmt='(a,1x,i20,1x,a)') ' - MAX_PTRS=',max_ptrs,'...'
            write (unit=*,fmt='(a,1x,a,1x,a)') ' - Tabstr=',tabstr,'...'
            write (unit=*,fmt='(a)') '----------------------------------------------'
            write (unit=*,fmt='(a)') ' '
            call vtable_edio_i_sca(var,nv,igr,0,glob_id,var_len,var_len_global,max_ptrs    &
                                  ,tabstr)
         end if
         
         vt_info(nv,igr)%nptrs                  =  vt_info(nv,igr)%nptrs + 1
         iptr                                   =  vt_info(nv,igr)%nptrs
         vt_info(nv,igr)%vt_vector(iptr)%globid =  glob_id

         vt_info(nv,igr)%vt_vector(iptr)%sca_ip => var
         vt_info(nv,igr)%vt_vector(iptr)%varlen =  var_len
      end if
      return
   end subroutine vtable_edio_i_sca
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This sub-routine creates the pointers for output, when the variable is character. !
   !                                                                                       !
   !     VAR            - The pointer of the current state variable                        !
   !     NV             - The variable type number                                         !
   !     IGR            - The number of the current grid                                   !
   !     INIT           - Initialize the vt_info?                                          !
   !     GLOB_ID        - The global index of the data                                     !
   !     VAR_LEN        - The length of the states current vector                          !
   !     VAR_LEN_GLOBAL - The length of the entire dataset's vector                        !
   !     MAX_PTRS       - The maximum possible number of pointers necessary for this       !
   !                      variable                                                         !
   !     TABSTR         - The string describing the variables usage                        !
   !---------------------------------------------------------------------------------------!
   recursive subroutine vtable_edio_c_sca(var,nv,igr,init,glob_id,var_len,var_len_global   &
                                         ,max_ptrs,tabstr)
      use ed_max_dims, only : str_len ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      character(len=str_len)                 , target     :: var
      integer                                , intent(in) :: nv
      integer                                , intent(in) :: igr
      integer                                , intent(in) :: init
      integer                                , intent(in) :: glob_id
      integer                                , intent(in) :: var_len
      integer                                , intent(in) :: var_len_global
      integer                                , intent(in) :: max_ptrs
      character (len=*)                      , intent(in) :: tabstr
      !----- Local variables. -------------------------------------------------------------!
      integer                                             :: iptr
      character(len=str_len), dimension(10)               :: tokens
      character(len=8)                                    :: ctab
      integer                                             :: ntok
      integer                                             :: nt
      !----- Local constants. -------------------------------------------------------------!
      character (len=1)                      , parameter  :: toksep=':'
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Determine if this is the first time we view this variable.  If so, then fill  !
      ! some descriptors for the vtable and allocate some space for any pointers that may  !
      ! follow.                                                                            !
      !------------------------------------------------------------------------------------!
      if (init == 0) then

         !----- Count the number of variables. --------------------------------------------!
         num_var(igr) = num_var(igr) + 1
         call tokenize1(tabstr,tokens,ntok,toksep)

         vt_info(nv,igr)%name           = tokens(1)
         vt_info(nv,igr)%dtype          = 'c'  ! Character variable (scalar)
         vt_info(nv,igr)%nptrs          = 0
         vt_info(nv,igr)%var_len_global = var_len_global

         nullify(vt_info(nv,igr)%vt_vector)
         allocate(vt_info(nv,igr)%vt_vector(max_ptrs))
         read(tokens(2),fmt=*) vt_info(nv,igr)%idim_type

         vt_info(nv,igr)%ihist    = 0
         vt_info(nv,igr)%ianal    = 0
         vt_info(nv,igr)%imean    = 0
         vt_info(nv,igr)%ilite    = 0
         vt_info(nv,igr)%impti    = 0
         vt_info(nv,igr)%impt1    = 0
         vt_info(nv,igr)%impt2    = 0
         vt_info(nv,igr)%impt3    = 0
         vt_info(nv,igr)%irecycle = 0
         vt_info(nv,igr)%imont    = 0
         vt_info(nv,igr)%idail    = 0
         vt_info(nv,igr)%idcyc    = 0
         vt_info(nv,igr)%iyear    = 0
         vt_info(nv,igr)%iopti    = 0
         
         do nt=3,ntok
            ctab=tokens(nt)

            select case (trim(ctab))
            case('hist') 
               vt_info(nv,igr)%ihist    = 1

            case('anal') 
               vt_info(nv,igr)%ianal    = 1

            case('lite') 
               vt_info(nv,igr)%ilite    = 1

            case('mpti') 
               vt_info(nv,igr)%impti    = 1

            case('mpt1') 
               vt_info(nv,igr)%impt1    = 1

            case('mpt2') 
               vt_info(nv,igr)%impt2    = 1

            case('mpt3') 
               vt_info(nv,igr)%impt3    = 1

            case('recycle') 
               vt_info(nv,igr)%irecycle = 1

            case('mont') 
               vt_info(nv,igr)%imont    = 1

            case('dail') 
               vt_info(nv,igr)%idail    = 1

            case('dcyc') 
               vt_info(nv,igr)%idcyc    = 1

            case('year') 
               vt_info(nv,igr)%iyear    = 1

            case('opti') 
               vt_info(nv,igr)%iopti    = 1

            case default
               print*, 'Illegal table specification for var:', tokens(1),ctab
               call fatal_error('Bad var table','vtable_edio_r','ed_var_tables.f90')
            end select
         end do
      else
         !---------------------------------------------------------------------------------!
         !     Make sure that vt_info is associated. If not, call the function with        !
         ! init = 0 then do this part.  Since I think this should never happen, I will     !
         ! also make a fuss to warn the user.                                              !
         !---------------------------------------------------------------------------------!
         if (.not.associated(vt_info(nv,igr)%vt_vector)) then
            write (unit=*,fmt='(a)') ' '
            write (unit=*,fmt='(a)') '----------------------------------------------'
            write (unit=*,fmt='(a)') ' WARNING! WARNING! WARNING! WARNING! WARNING! '
            write (unit=*,fmt='(a)') ' WARNING! WARNING! WARNING! WARNING! WARNING! '
            write (unit=*,fmt='(a)') ' WARNING! WARNING! WARNING! WARNING! WARNING! '
            write (unit=*,fmt='(a)') '----------------------------------------------'
            write (unit=*,fmt='(a)') ' - Subroutine vtable_edio_r (file ed_var_tables.f90)'
            write (unit=*,fmt='(a,1x,i4,1x,a,1x,i2,1x,a)')                                 &
                                     ' - Vt_vector for variable',nv,'of grid',igr          &
                                     ,'is not associated              !'
            write (unit=*,fmt='(a)') ' - I will allocate it now.'
            write (unit=*,fmt='(a,1x,i20,1x,a)') ' - MAX_PTRS=',max_ptrs,'...'
            write (unit=*,fmt='(a,1x,a,1x,a)') ' - Tabstr=',tabstr,'...'
            write (unit=*,fmt='(a)') '----------------------------------------------'
            write (unit=*,fmt='(a)') ' '
            call vtable_edio_c_sca(var,nv,igr,0,glob_id,var_len,var_len_global,max_ptrs    &
                                  ,tabstr)
         end if
         
         vt_info(nv,igr)%nptrs                  =  vt_info(nv,igr)%nptrs + 1
         iptr                                   =  vt_info(nv,igr)%nptrs
         vt_info(nv,igr)%vt_vector(iptr)%globid =  glob_id

         vt_info(nv,igr)%vt_vector(iptr)%sca_cp => var
         vt_info(nv,igr)%vt_vector(iptr)%varlen =  var_len
      end if
      return
   end subroutine vtable_edio_c_sca
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   subroutine metadata_edio(nv,igr,lname,units,dimstr)
      implicit none
     
      integer :: nv,igr
      character (len=*) :: lname
      character (len=*) :: units
      character (len=*) :: dimstr
     
      vt_info(nv,igr)%lname = trim(lname)
      vt_info(nv,igr)%units = trim(units)
      vt_info(nv,igr)%dimlab = trim(dimstr)
     
      return
   end subroutine metadata_edio
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   subroutine nullify_vt_vector_pointers(vt_vec)
      implicit none
      type(var_table_vector), target :: vt_vec
      
      if (associated(vt_vec%var_rp)) nullify(vt_vec%var_rp)
      if (associated(vt_vec%var_ip)) nullify(vt_vec%var_ip)
      if (associated(vt_vec%var_cp)) nullify(vt_vec%var_cp)
      if (associated(vt_vec%var_dp)) nullify(vt_vec%var_dp)
      return
   end subroutine nullify_vt_vector_pointers
   !=======================================================================================!
   !=======================================================================================!
end module ed_var_tables
!==========================================================================================!
!==========================================================================================!

