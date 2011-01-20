!==========================================================================================!
!   Module ed_var_tables: this module contains all information about variables, includ- !
! ing the kind of variable and whether they need to be passed in parallel runs, in which   !
! output they should be included etc.                                                      !
!   The following list is the number used to identify the type of variable.                !
!                                                                                          !
!     !----- POLYGON: --------------------------------------!                              !
!     ! 11  : rank 1 : real                                 !                              !
!     ! 12  : rank 2 : s-layer                              !                              !
!     ! 13  : rank 2 : w-layer                              !                              !
!     ! 14  : rank 2 : pft                                  !                              !
!     ! 146 : rank 3 : npft,ndbh                            !                              !
!     ! 15  : rank 2 : disturbance                          !                              !
!     ! 16  : rank 2 : dbh                                  !                              !
!     ! 17  : rank 2 : age                                  !                              !
!     ! 155 : rank 2 : max_lu_years                         !                              !
!     ! 156 : rank 3 : max_lu_years, num_lu_transitions     !                              !
!     !-----------------------------------------------------!                              !
!                                                                                          !
!     !----- SITE: -----------------------------------------!                              !
!     ! 20  : rank 1 : integer                              !                              !
!     ! 21  : rank 1 : real                                 !                              !
!     ! 22  : rank 2 : s-layer                              !                              !
!     ! 23  : rank 2 : w-layer                              !                              !
!     ! 24  : rank 2 : pft                                  !                              !
!     ! 246 : rank 3 : pft, dbh                             !                              !
!     ! 25  : rank 2 : disturbance                          !                              !
!     ! 255 : rank 3 : disturbance,disturbance              !                              !
!     ! 26  : rank 2 : dbh                                  !                              !
!     ! 27  : rank 2 : age                                  !                              !
!     ! 28  : rank 2 : months                               !                              !
!     !-----------------------------------------------------!                              !
!                                                                                          !
!     !----- SITE: -----------------------------------------!                              !
!     ! 30  : rank 1 : integer                              !                              !
!     ! 31  : rank 1 : real                                 !                              !
!     ! 32  : rank 2 : s-layer                              !                              !
!     ! 33  : rank 2 : w-layer                              !                              !
!     ! 34  : rank 2 : pft                                  !                              !
!     ! 346 : rank 3 : pft,ff_dbh                           !                              !
!     ! 35  : rank 2 : disturbance                          !                              !
!     ! 36  : rank 2 : dbh                                  !                              !
!     ! 37  : rank 2 : age                                  !                              !
!     !-----------------------------------------------------!                              !
!                                                                                          !
!     !----- COHORT: ---------------------------------------!                              !
!     ! 41  : rank 1 : cohort (real)                        !                              !
!     ! 44  : rank 2 : cohort, pft                          !                              !
!     ! 46  : rank 2 : cohort, dbh                          !                              !
!     ! 47  : rank 2 : cohort, age                          !                              !
!     ! 49  : rank 2 : cohort, nmonths+1                    !                              !
!     !-----------------------------------------------------!                              !
!                                                                                          !
!     !-----------------------------------------------------!                              !
!     ! 90 is a special flag for scalars                    !                              !
!     !-----------------------------------------------------!                              !
!                                                                                          !
!------------------------------------------------------------------------------------------!
module ed_var_tables
   use ed_max_dims, only : str_len
   !---------------------------------------------------------------------------------------!
   !    Define data type for main variable table                                           !
   !---------------------------------------------------------------------------------------!
   integer,parameter :: maxvars = 1500
  
   !---------------------------------------------------------------------------------------!
   type var_table
      integer :: idim_type
      integer :: nptrs
      integer :: ihist,ianal,imean,ilite,impti,impt1,impt2,impt3,irecycle,iyear,iopti
      character (len=64) :: name
      character (len=2) :: dtype
      integer :: imont,idail
      logical :: first
      integer :: var_len_global
      character (len=64) :: lname   ! Long name for description in file
      character (len=16)  :: units   ! Unit description of the data
      character (len=64)  :: dimlab
      !----- Multiple pointer defs (maxptrs) ----------------------------------------------!
      type(var_table_vector),pointer,dimension(:) :: vt_vector
   end type var_table
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   type var_table_vector
      real   , pointer            :: var_rp
      integer, pointer            :: var_ip
      character (len=str_len),pointer :: var_cp
      real(kind=8),pointer        :: var_dp
      integer                     :: globid
      integer                     :: varlen
   end type var_table_vector
   !---------------------------------------------------------------------------------------!

  
   !----- Main variable table allocated to (maxvars,maxgrds) ------------------------------!
   type(var_table), allocatable :: vt_info(:,:)


   !----- Number of variables for each grid, allocated to "ngrids". -----------------------!
   integer, allocatable :: num_var(:)

   contains
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   recursive subroutine vtable_edio_r( &
        var,      &    ! The pointer of the current state variable
        nv,       &    ! The variable type number
        igr,      &    ! The number of the current grid
        init,     &    ! Initialize the vt_info?
        glob_id,  &    ! The global index of the data
        var_len,  &    ! The length of the states current vector
        var_len_global, & ! THe length of the entire dataset's vector
        max_ptrs,  &    ! The maximum possible number of pointers
        ! necessary for this variable
        tabstr)        ! The string describing the variables usage
     
     implicit none
     
     real,target :: var
     
     integer :: init
     integer :: var_len,var_len_global,max_ptrs,glob_id,iptr,igr
     character (len=*) :: tabstr
     
     character (len=1), parameter ::toksep=':'
     character (len=128) ::tokens(10)
     character (len=8) :: ctab
     integer :: ntok,nt,nv
     
     ! ------------------------------------------------
     ! Determine if this is the first
     ! time we view this variable.  If so, then
     ! fill some descriptors for the vtable
     ! and allocate some space for any pointers
     ! that may follow
     ! ------------------------------------------------
     
     if (init == 0) then
        
        ! Count the number of variables
        num_var(igr) = num_var(igr) + 1
       
        call tokenize1(tabstr,tokens,ntok,toksep)
        
        vt_info(nv,igr)%name=tokens(1)

!       print*,num_var(igr),nv,trim(vt_info(nv,igr)%name)

        vt_info(nv,igr)%dtype='r'  ! This is a real variable

        vt_info(nv,igr)%nptrs = 0
        
        vt_info(nv,igr)%var_len_global = var_len_global

        nullify(vt_info(nv,igr)%vt_vector)
        allocate(vt_info(nv,igr)%vt_vector(max_ptrs))
        
        read(tokens(2),*) vt_info(nv,igr)%idim_type
        
        vt_info(nv,igr)%ihist=0
        vt_info(nv,igr)%ianal=0
        vt_info(nv,igr)%imean=0
        vt_info(nv,igr)%ilite=0
        vt_info(nv,igr)%impti=0
        vt_info(nv,igr)%impt1=0
        vt_info(nv,igr)%impt2=0
        vt_info(nv,igr)%impt3=0
        vt_info(nv,igr)%irecycle=0
        vt_info(nv,igr)%imont=0
        vt_info(nv,igr)%idail=0
        vt_info(nv,igr)%iyear=0
        vt_info(nv,igr)%iopti=0
        
        do nt=3,ntok
           ctab=tokens(nt)
           
           select case (trim(ctab))
           case('hist') 
              vt_info(nv,igr)%ihist=1
           case('anal') 
              vt_info(nv,igr)%ianal=1
           case('lite') 
              vt_info(nv,igr)%ilite=1
           case('mpti') 
              vt_info(nv,igr)%impti=1
           case('mpt1') 
              vt_info(nv,igr)%impt1=1
           case('mpt2') 
              vt_info(nv,igr)%impt2=1
           case('mpt3') 
              vt_info(nv,igr)%impt3=1
           case('recycle') 
              vt_info(nv,igr)%irecycle=1
           case('mont') 
              vt_info(nv,igr)%imont=1
           case('dail') 
              vt_info(nv,igr)%idail=1
           case('year') 
              vt_info(nv,igr)%iyear=1
           case('opti') 
              vt_info(nv,igr)%iopti=1
           case default
              print*, 'Illegal table specification for var:', tokens(1),ctab
              call fatal_error('Bad var table','vtable_edio_r','ed_var_tables.f90')
           end select
           
        enddo
        
        ! Set the first pass logical to false
         
     else
        !    Make sure that vt_info is associated. If not, call the function with init = 0 then 
        ! do this part. Since I think this should never happen, I will also make a fuss to warn 
        ! the user
        if (.not.associated(vt_info(nv,igr)%vt_vector)) then
          write (unit=*,fmt='(a)') ' '
          write (unit=*,fmt='(a)') '!-------------------------------------------------------------------------!'
          write (unit=*,fmt='(a)') '! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! !'
          write (unit=*,fmt='(a)') '! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! !'
          write (unit=*,fmt='(a)') '! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! !'
          write (unit=*,fmt='(a)') '!-------------------------------------------------------------------------!'
          write (unit=*,fmt='(a)') '! In subroutine vtable_edio_r (file ed_var_tables.f90)                 !'
          write (unit=*,fmt='(a,1x,i4,1x,a,1x,i2,1x,a)')  &
                                   '! Vt_vector for variable',nv,'of grid',igr,'is not associated                !'
          write (unit=*,fmt='(a)') '! I will allocate it now.                                                 !'
          write (unit=*,fmt='(a,1x,i20,1x,a)') '! MAX_PTRS=',max_ptrs,'...'
          write (unit=*,fmt='(a,1x,a,1x,a)') '! Tabstr=',tabstr,'...'
          write (unit=*,fmt='(a)') '!-------------------------------------------------------------------------!'
          write (unit=*,fmt='(a)') ' '
          call vtable_edio_r(var,nv,igr,0,glob_id,var_len,var_len_global,max_ptrs,tabstr)
        end if
        
        vt_info(nv,igr)%nptrs = vt_info(nv,igr)%nptrs + 1
        iptr = vt_info(nv,igr)%nptrs
        
        vt_info(nv,igr)%vt_vector(iptr)%globid = glob_id
        vt_info(nv,igr)%vt_vector(iptr)%var_rp   => var
        vt_info(nv,igr)%vt_vector(iptr)%varlen = var_len

     end if
     
     
     return
   end subroutine vtable_edio_r
   !=======================================================================================!
   !=======================================================================================!


 recursive subroutine vtable_edio_d( &
       var,      &    ! The pointer of the current state variable
       nv,       &    ! The variable type number
       igr,      &    ! The number of the current grid
       init,     &    ! Initialize the vt_info?
       glob_id,  &    ! The global index of the data
       var_len,  &    ! The length of the states current vector
       var_len_global, & ! THe length of the entire dataset's vector
       max_ptrs,  &    ! The maximum possible number of pointers
       ! necessary for this variable
       tabstr)        ! The string describing the variables usage
    
    implicit none
    
    real(kind=8),target :: var
    
    integer :: init
    integer :: var_len,var_len_global,max_ptrs,glob_id,iptr,igr
    character (len=*) :: tabstr
    
    character (len=1), parameter ::toksep=':'
    character (len=128) ::tokens(10)
    character (len=8) :: ctab
    integer :: ntok,nt,nv
    
    ! ------------------------------------------------
    ! Determine if this is the first
    ! time we view this variable.  If so, then
    ! fill some descriptors for the vtable
    ! and allocate some space for any pointers
    ! that may follow
    ! ------------------------------------------------
    
    if (init == 0) then
       
       ! Count the number of variables
       num_var(igr) = num_var(igr) + 1
      
       call tokenize1(tabstr,tokens,ntok,toksep)
       
       vt_info(nv,igr)%name=tokens(1)

!       print*,num_var(igr),nv,trim(vt_info(nv,igr)%name)

       vt_info(nv,igr)%dtype='d'  ! This is a double precision variable

       vt_info(nv,igr)%nptrs = 0
       
       vt_info(nv,igr)%var_len_global = var_len_global

       nullify(vt_info(nv,igr)%vt_vector)
       allocate(vt_info(nv,igr)%vt_vector(max_ptrs))
       
       read(tokens(2),*) vt_info(nv,igr)%idim_type
       
       vt_info(nv,igr)%ihist=0
       vt_info(nv,igr)%ianal=0
       vt_info(nv,igr)%imean=0
       vt_info(nv,igr)%ilite=0
       vt_info(nv,igr)%impti=0
       vt_info(nv,igr)%impt1=0
       vt_info(nv,igr)%impt2=0
       vt_info(nv,igr)%impt3=0
       vt_info(nv,igr)%irecycle=0
       vt_info(nv,igr)%imont=0
       vt_info(nv,igr)%idail=0
       vt_info(nv,igr)%iyear=0
       vt_info(nv,igr)%iopti=0
       
       do nt=3,ntok
          ctab=tokens(nt)
          
          select case (trim(ctab))
          case('hist') 
             vt_info(nv,igr)%ihist=1
          case('anal') 
             vt_info(nv,igr)%ianal=1
          case('lite') 
             vt_info(nv,igr)%ilite=1
          case('mpti') 
             vt_info(nv,igr)%impti=1
          case('mpt1') 
             vt_info(nv,igr)%impt1=1
          case('mpt2') 
             vt_info(nv,igr)%impt2=1
          case('mpt3') 
             vt_info(nv,igr)%impt3=1
          case('recycle') 
             vt_info(nv,igr)%irecycle=1
          case('mont') 
             vt_info(nv,igr)%imont=1
          case('dail') 
             vt_info(nv,igr)%idail=1
          case('year') 
             vt_info(nv,igr)%iyear=1
          case('opti') 
             vt_info(nv,igr)%iopti=1
          case default
             print*, 'Illegal table specification for var:', tokens(1),ctab
             call fatal_error('Bad var table','vtable_edio_r','ed_var_tables.f90')
          end select
          
       enddo
       
       ! Set the first pass logical to false
        
    else
       !    Make sure that vt_info is associated. If not, call the function with init = 0 then 
       ! do this part. Since I think this should never happen, I will also make a fuss to warn 
       ! the user
       if (.not.associated(vt_info(nv,igr)%vt_vector)) then
         write (unit=*,fmt='(a)') ' '
         write (unit=*,fmt='(a)') '!-------------------------------------------------------------------------!'
         write (unit=*,fmt='(a)') '! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! !'
         write (unit=*,fmt='(a)') '! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! !'
         write (unit=*,fmt='(a)') '! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! !'
         write (unit=*,fmt='(a)') '!-------------------------------------------------------------------------!'
         write (unit=*,fmt='(a)') '! In subroutine vtable_edio_r (file ed_var_tables.f90)                 !'
         write (unit=*,fmt='(a,1x,i4,1x,a,1x,i2,1x,a)')  &
                                  '! Vt_vector for variable',nv,'of grid',igr,'is not associated                !'
         write (unit=*,fmt='(a)') '! I will allocate it now.                                                 !'
         write (unit=*,fmt='(a,1x,i20,1x,a)') '! MAX_PTRS=',max_ptrs,'...'
         write (unit=*,fmt='(a,1x,a,1x,a)')   '! Tabstr=',tabstr,'...'
         write (unit=*,fmt='(a)') '!-------------------------------------------------------------------------!'
         write (unit=*,fmt='(a)') ' '
         call vtable_edio_d(var,nv,igr,0,glob_id,var_len,var_len_global,max_ptrs,tabstr)
       end if
       
       vt_info(nv,igr)%nptrs = vt_info(nv,igr)%nptrs + 1
       iptr = vt_info(nv,igr)%nptrs
       
       vt_info(nv,igr)%vt_vector(iptr)%globid = glob_id
       vt_info(nv,igr)%vt_vector(iptr)%var_dp   => var
       vt_info(nv,igr)%vt_vector(iptr)%varlen = var_len

    end if
    
    
    return
  end subroutine vtable_edio_d
!==============================================================================!
!==============================================================================!


!==============================================================================!
!==============================================================================!
  recursive subroutine vtable_edio_i( &
       var,      &    ! The pointer of the current state variable
       nv,       &    ! The variable type number
       igr,      &    ! The number of the current grid
       init,     &    ! Initialize the vt_info?
       glob_id,  &    ! The global index of the data
       var_len,  &    ! The length of the states current vector
       var_len_global, & ! THe length of the entire dataset's vector
       max_ptrs,  &    ! The maximum possible number of pointers
       ! necessary for this variable
       tabstr)        ! The string describing the variables usage
    
    implicit none
    
    integer,target :: var
    
    integer :: init
    integer :: var_len,var_len_global,max_ptrs,glob_id,iptr,igr
    character (len=*) :: tabstr
    
    character (len=1), parameter ::toksep=':'
    character (len=128) ::tokens(10)
    character (len=8) :: ctab
    
    integer :: ntok,nt,nv
    
    ! ------------------------------------------------
    ! Determine if this is the first
    ! time we view this variable.  If so, then
    ! fill some descriptors for the vtable
    ! and allocate some space for any pointers
    ! that may follow
    ! ------------------------------------------------
    if (init == 0) then
       
       ! Count the number of variables
       num_var(igr) = num_var(igr) + 1

       call tokenize1(tabstr,tokens,ntok,toksep)
       
       vt_info(nv,igr)%name=tokens(1)

!       print*,num_var(igr),nv,trim(vt_info(nv,igr)%name)

       vt_info(nv,igr)%dtype='i'  ! This is an integer variable

       vt_info(nv,igr)%nptrs = 0
       
       vt_info(nv,igr)%var_len_global = var_len_global

       nullify(vt_info(nv,igr)%vt_vector)
       allocate(vt_info(nv,igr)%vt_vector(max_ptrs))
       
       read(tokens(2),*) vt_info(nv,igr)%idim_type
       
       vt_info(nv,igr)%ihist=0
       vt_info(nv,igr)%ianal=0
       vt_info(nv,igr)%imean=0
       vt_info(nv,igr)%ilite=0
       vt_info(nv,igr)%impti=0
       vt_info(nv,igr)%impt1=0
       vt_info(nv,igr)%impt2=0
       vt_info(nv,igr)%impt3=0
       vt_info(nv,igr)%irecycle=0
       vt_info(nv,igr)%imont=0
       vt_info(nv,igr)%idail=0
       vt_info(nv,igr)%iyear=0
       vt_info(nv,igr)%iopti=0
       
       do nt=3,ntok
          ctab=tokens(nt)
          
          select case (trim(ctab))
          case('hist') 
             vt_info(nv,igr)%ihist=1
          case('anal') 
             vt_info(nv,igr)%ianal=1
          case('lite') 
             vt_info(nv,igr)%ilite=1
          case('mpti') 
             vt_info(nv,igr)%impti=1
          case('mpt1') 
             vt_info(nv,igr)%impt1=1
          case('mpt2') 
             vt_info(nv,igr)%impt2=1
          case('mpt3') 
             vt_info(nv,igr)%impt3=1
          case('recycle') 
             vt_info(nv,igr)%irecycle=1
          case('mont') 
             vt_info(nv,igr)%imont=1
          case('dail') 
             vt_info(nv,igr)%idail=1
          case('year') 
             vt_info(nv,igr)%iyear=1
          case('opti') 
             vt_info(nv,igr)%iopti=1
          case default
             print*, 'Illegal table specification for var:', tokens(1),ctab
             call fatal_error('Bad var table','vtable_edio_i','ed_var_tables.f90')
          end select
          
       enddo
       
       ! Set the first pass logical to false
       
       vt_info(nv,igr)%first = .false.
       
    else
       !    Make sure that vt_info is associated. If not, call the function with init = 0 then 
       ! do this part. Since I think this should never happen, I will also make a fuss to warn 
       ! the user
       if (.not.associated(vt_info(nv,igr)%vt_vector)) then
         write (unit=*,fmt='(a)') ' '
         write (unit=*,fmt='(a)') '!-------------------------------------------------------------------------!'
         write (unit=*,fmt='(a)') '! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! !'
         write (unit=*,fmt='(a)') '! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! !'
         write (unit=*,fmt='(a)') '! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! !'
         write (unit=*,fmt='(a)') '!-------------------------------------------------------------------------!'
         write (unit=*,fmt='(a)') '! In subroutine vtable_edio_i (file ed_var_tables.f90)                 !'
         write (unit=*,fmt='(a,1x,i4,1x,a,1x,i2,1x,a)')  &
                                  '! Vt_vector for variable',nv,'of grid',igr,'is not associated                !'
         write (unit=*,fmt='(a)') '! I will allocate it now.                                                 !'
         write (unit=*,fmt='(a)') '!-------------------------------------------------------------------------!'
         write (unit=*,fmt='(a)') ' '
         call vtable_edio_i(var,nv,igr,0,glob_id,var_len,var_len_global,max_ptrs,tabstr)
       end if
       vt_info(nv,igr)%nptrs = vt_info(nv,igr)%nptrs + 1
       iptr = vt_info(nv,igr)%nptrs
       vt_info(nv,igr)%vt_vector(iptr)%globid = glob_id
       vt_info(nv,igr)%vt_vector(iptr)%var_ip   => var
       vt_info(nv,igr)%vt_vector(iptr)%varlen = var_len

    endif

    
    return
  end subroutine vtable_edio_i

  ! =====================================================

  recursive subroutine vtable_edio_c( &
       var,      &    ! The pointer of the current state variable
       nv,       &    ! The variable type number
       igr,      &    ! The number of the current grid
       init,     &    ! Initialize the vt_info?
       glob_id,  &    ! The global index of the data
       var_len,  &    ! The length of the states current vector
       var_len_global, & ! THe length of the entire dataset's vector
       max_ptrs,  &    ! The maximum possible number of pointers
       ! necessary for this variable
       tabstr)        ! The string describing the variables usage
    use ed_max_dims, only : str_len
    implicit none
    
    character (len=str_len),target :: var
    
    integer :: init
    integer :: var_len,var_len_global,max_ptrs,glob_id,iptr,igr
    character (len=*) :: tabstr
    
    character (len=1), parameter ::toksep=':'
    character (len=128) ::tokens(10)
    character (len=8) :: ctab
    integer :: ntok,nt,nv
    
    ! ------------------------------------------------
    ! Determine if this is the first
    ! time we view this variable.  If so, then
    ! fill some descriptors for the vtable
    ! and allocate some space for any pointers
    ! that may follow
    ! ------------------------------------------------
    
    if (init == 0) then
       
       ! Count the number of variables
       num_var(igr) = num_var(igr) + 1
      
       call tokenize1(tabstr,tokens,ntok,toksep)
       
       vt_info(nv,igr)%name=tokens(1)

!       print*,num_var(igr),nv,trim(vt_info(nv,igr)%name)

       vt_info(nv,igr)%dtype='c'  ! This is a string variable

       vt_info(nv,igr)%nptrs = 0
       
       vt_info(nv,igr)%var_len_global = var_len_global

       nullify(vt_info(nv,igr)%vt_vector)
       allocate(vt_info(nv,igr)%vt_vector(max_ptrs))
       
       read(tokens(2),*) vt_info(nv,igr)%idim_type
       
       vt_info(nv,igr)%ihist=0
       vt_info(nv,igr)%ianal=0
       vt_info(nv,igr)%imean=0
       vt_info(nv,igr)%ilite=0
       vt_info(nv,igr)%impti=0
       vt_info(nv,igr)%impt1=0
       vt_info(nv,igr)%impt2=0
       vt_info(nv,igr)%impt3=0
       vt_info(nv,igr)%irecycle=0
       vt_info(nv,igr)%imont=0
       vt_info(nv,igr)%idail=0
       vt_info(nv,igr)%iyear=0
       vt_info(nv,igr)%iopti=0
       
       do nt=3,ntok
          ctab=tokens(nt)
          
          select case (trim(ctab))
          case('hist') 
             vt_info(nv,igr)%ihist=1
          case('anal') 
             vt_info(nv,igr)%ianal=1
          case('lite') 
             vt_info(nv,igr)%ilite=1
          case('mpti') 
             vt_info(nv,igr)%impti=1
          case('mpt1') 
             vt_info(nv,igr)%impt1=1
          case('mpt2') 
             vt_info(nv,igr)%impt2=1
          case('mpt3') 
             vt_info(nv,igr)%impt3=1
          case('recycle') 
             vt_info(nv,igr)%irecycle=1
          case('mont') 
             vt_info(nv,igr)%imont=1
          case('dail') 
             vt_info(nv,igr)%idail=1
          case('year') 
             vt_info(nv,igr)%iyear=1
          case('opti') 
             vt_info(nv,igr)%iopti=1
          case default
             print*, 'Illegal table specification for var:', tokens(1),ctab
             call fatal_error('Bad var table','vtable_edio_c','ed_var_tables.f90')
          end select
          
       enddo
       
       ! Set the first pass logical to false
        
    else
       !    Make sure that vt_info is associated. If not, call the function with init = 0 then 
       ! do this part. Since I think this should never happen, I will also make a fuss to warn 
       ! the user
       if (.not.associated(vt_info(nv,igr)%vt_vector)) then
         write (unit=*,fmt='(a)') ' '
         write (unit=*,fmt='(a)') '!-------------------------------------------------------------------------!'
         write (unit=*,fmt='(a)') '! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! !'
         write (unit=*,fmt='(a)') '! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! !'
         write (unit=*,fmt='(a)') '! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! !'
         write (unit=*,fmt='(a)') '!-------------------------------------------------------------------------!'
         write (unit=*,fmt='(a)') '! In subroutine vtable_edio_c (file ed_var_tables.f90)                 !'
         write (unit=*,fmt='(a,1x,i4,1x,a,1x,i2,1x,a)')  &
                                  '! Vt_vector for variable',nv,'of grid',igr,'is not associated                !'
         write (unit=*,fmt='(a)') '! I will allocate it now.                                                 !'
         write (unit=*,fmt='(a)') '!-------------------------------------------------------------------------!'
         write (unit=*,fmt='(a)') ' '
         call vtable_edio_c(var,nv,igr,0,glob_id,var_len,var_len_global,max_ptrs,tabstr)
       end if
       
       vt_info(nv,igr)%nptrs = vt_info(nv,igr)%nptrs + 1
       iptr = vt_info(nv,igr)%nptrs
       
       vt_info(nv,igr)%vt_vector(iptr)%globid = glob_id
       vt_info(nv,igr)%vt_vector(iptr)%var_cp   => var
       vt_info(nv,igr)%vt_vector(iptr)%varlen = var_len

    end if
    
    
    return
  end subroutine vtable_edio_c

  ! =====================================================
  
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

  ! =====================================================
  
  subroutine nullify_vt_vector_pointers(vt_vec)
     implicit none
     type(var_table_vector), target :: vt_vec
     
     if (associated(vt_vec%var_rp)) nullify(vt_vec%var_rp)
     if (associated(vt_vec%var_ip)) nullify(vt_vec%var_ip)
     if (associated(vt_vec%var_cp)) nullify(vt_vec%var_cp)
     if (associated(vt_vec%var_dp)) nullify(vt_vec%var_dp)
     return
  end subroutine nullify_vt_vector_pointers

end module ed_var_tables

