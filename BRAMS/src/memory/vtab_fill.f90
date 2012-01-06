!============================= Change Log =================================================!
! 5.0.0                                                                                    !
!                                                                                          !
!==========================================================================================!
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved                       !
!  Regional Atmospheric Modeling System - RAMS                                             !
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine will assign pointers to the variable table, and assign the values of !
! several flags to determine whether/when a variable should be written to files, and to    !
! be exchanged amongst nodes.                                                              !
!------------------------------------------------------------------------------------------!
subroutine vtables2(var,varm,ng,npts,imean,tabstr)
   use grid_dims , only : str_len ! ! intent(in)
   use var_tables
      
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                                 , intent(in) :: ng
   integer                                 , intent(in) :: npts
   integer                                 , intent(in) :: imean
   real                  , dimension (npts), target     :: var
   real                  , dimension (npts), target     :: varm
   character(len=*)                        , intent(in) :: tabstr
   !---- Local variables. -----------------------------------------------------------------!
   character(len=str_len)                               :: line
   character(len=1)                                     :: cdimen
   character(len=1)                                     :: ctype
   character(len=32)     , dimension(10)                :: tokens
   character(len=8)                                     :: cname
   character(len=8)                                     :: ctab
   integer                                              :: ntok
   integer                                              :: nt
   integer                                              :: nv
   character(len=1)                        , parameter  :: toksep = ':'
   !---------------------------------------------------------------------------------------!


   !----- Split the tokens and save them in tokens. ---------------------------------------!
   call tokenize1(tabstr,tokens,ntok,toksep)
   !---------------------------------------------------------------------------------------!



   !----- Add the variable to the list, and save nv as the variable "ID". -----------------!
   num_var(ng)=num_var(ng)+1
   nv=num_var(ng)
   !---------------------------------------------------------------------------------------!



   !----- Point the variables. ------------------------------------------------------------!
   vtab_r(nv,ng)%var_p => var
   vtab_r(nv,ng)%var_m => varm
   !---------------------------------------------------------------------------------------!


   !----- Token(1) must be variable name. -------------------------------------------------!
   vtab_r(nv,ng)%name = tokens(1)
   !---------------------------------------------------------------------------------------!


   !----- Variable size. ------------------------------------------------------------------!
   vtab_r(nv,ng)%npts = npts
   !---------------------------------------------------------------------------------------!


   !----- Token(2) is the dimension, read it so it becomes an integer. --------------------!
   read (tokens(2),fmt=*) vtab_r(nv,ng)%idim_type
   !---------------------------------------------------------------------------------------!



   !----- Initialise all flags as zeroes (or imean). --------------------------------------!
   vtab_r(nv,ng)%ihist    = 0
   vtab_r(nv,ng)%ianal    = 0
   vtab_r(nv,ng)%imean    = imean
   vtab_r(nv,ng)%ilite    = 0
   vtab_r(nv,ng)%impti    = 0
   vtab_r(nv,ng)%impt1    = 0
   vtab_r(nv,ng)%impt2    = 0
   vtab_r(nv,ng)%impt3    = 0
   vtab_r(nv,ng)%iadvt    = 0
   vtab_r(nv,ng)%iadvu    = 0
   vtab_r(nv,ng)%iadvv    = 0
   vtab_r(nv,ng)%iadvw    = 0
   vtab_r(nv,ng)%irecycle = 0
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !      Loop over the remaining tokens, and switch the flag to 1 for those actions that  !
   ! are listed in the token list.                                                         !
   !---------------------------------------------------------------------------------------!
   do nt=3,ntok
      ctab=tokens(nt)         
      
      select case(trim(ctab))
      case ('hist')
         vtab_r(nv,ng)%ihist    = 1

      case ('anal')
         vtab_r(nv,ng)%ianal    = 1

      case ('lite')
         vtab_r(nv,ng)%ilite    = 1

      case ('mpti')
         vtab_r(nv,ng)%impti    = 1

      case ('mpt1')
         vtab_r(nv,ng)%impt1    = 1

      case ('mpt2')
         vtab_r(nv,ng)%impt2    = 1

      case ('mpt3')
         vtab_r(nv,ng)%impt3    = 1

      case ('advt')
         vtab_r(nv,ng)%iadvt    = 1

      case ('advu')
         vtab_r(nv,ng)%iadvu    = 1

      case ('advv')
         vtab_r(nv,ng)%iadvv    = 1

      case ('advw')
         vtab_r(nv,ng)%iadvw    = 1

      case ('recycle')
         vtab_r(nv,ng)%irecycle = 1

      case default

         write(unit=*,fmt='(3(a,1x))')  'Illegal table specification for var:'             &
                                       ,tokens(1),ctab
         call abort_run('Bad settings.','vtables2','vtab_fill')
      end select
   end do
   !---------------------------------------------------------------------------------------!

   return
end subroutine vtables2
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine lite_varset(proc_type)

   use var_tables, only : nvgrids     & ! intent(in)
                        , vtab_r      & ! intent(inout)
                        , num_var     ! ! intent(in)
   use io_params , only : nlite_vars  & ! INTENT(IN)
                        , lite_vars   ! ! INTENT(IN)
   implicit none
  
   !----- Arguments. ----------------------------------------------------------------------!
   integer, intent(in) :: proc_type  
   !----- Local variables. ----------------------------------------------------------------!
   integer             :: nv
   integer             :: ng
   integer             :: nvl
   integer             :: ifound
   !---------------------------------------------------------------------------------------!

  
   !---------------------------------------------------------------------------------------!
   !      Loop over each variable input in namelist "LITE_VARS" and set lite flag in       !
   ! var_tables.                                                                           !
   !---------------------------------------------------------------------------------------!
   do ng = 1,nvgrids   
      vtab_r(1:num_var(ng),ng)%ilite = 0
   end do
   !---------------------------------------------------------------------------------------!
  

   !---------------------------------------------------------------------------------------!
   do nvl=1,nlite_vars
      ifound=0
      
      do ng=1,nvgrids
         do nv=1,num_var(ng)
            if (vtab_r(nv,ng)%name == lite_vars(nvl) ) then
               vtab_r(nv,ng)%ilite = 1
               ifound              = 1
            end if
         end do
      end do


      if (proc_type==0 .or. proc_type==1) then 
         !----- Output only in master process. --------------------------------------------!
         if (ifound == 0) then
            write(unit=*,fmt='(a)')  '!---------------------------------------------------'
            write(unit=*,fmt='(4a)') '! LITE_VARS variable does not exist in main '        &
                                        ,'variable table: -->',lite_vars(nvl),'<--'
            write(unit=*,fmt='(a)')  '!---------------------------------------------------'
         else
            write(unit=*,fmt='(a)')  '!---------------------------------------------------'
            write(unit=*,fmt='(3a)') '! LITE_VARS variable added -->', trim(lite_vars(nvl))
            write(unit=*,fmt='(a)')  '!---------------------------------------------------'
         end if
      end if
   end do

  return
end subroutine lite_varset
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine vtables_scalar(npts,varp,vart,ng,tabstr)
   use grid_dims , only : str_len ! ! intent(in)
   use var_tables
      
   implicit none
   !------ Arguments. ---------------------------------------------------------------------!
   integer                                 , intent(in) :: npts
   real                   , dimension(npts), target     :: varp
   real                   , dimension(npts), target     :: vart
   integer                                 , intent(in) :: ng
   character (len=*)                       , intent(in) :: tabstr
   !------ Local variables. ---------------------------------------------------------------!
   character (len=str_len)                              :: line
   character (len=32)     , dimension(10)               :: tokens
   character (len=16)                                   :: cname
   character (len=16)                                   :: ctab
   integer                                              :: ntok
   integer                                              :: nv
   integer                                              :: isnum
   integer                                              :: ns
   character (len=1)      , parameter                   :: toksep = ':'
   !---------------------------------------------------------------------------------------!


   !------ Split the tokens. --------------------------------------------------------------!
   call tokenize1(tabstr,tokens,ntok,toksep)
   cname = tokens(1)
   isnum = 0
   !---------------------------------------------------------------------------------------!


   !------ Fill in existing table slot or make new scalar slot. ---------------------------!
   num_scalar(ng)         = num_scalar(ng)+1
   nv                     = num_scalar(ng)
   scalar_tab(nv,ng)%name = cname
   !---------------------------------------------------------------------------------------!


   !------ Match the variable and the tendency. -------------------------------------------!
   scalar_tab(nv,ng)%var_p => varp
   scalar_tab(nv,ng)%var_t => vart
   !---------------------------------------------------------------------------------------!

   return
end subroutine vtables_scalar
!==========================================================================================!
!==========================================================================================!





!==========================================================================================!
!==========================================================================================!
subroutine vtables_scalar_new(varp,vart,ng,tabstr,elements)

   use var_tables
   use grid_dims , only : str_len      
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                                     , intent(in) :: elements
   real                   , dimension(elements), target     :: varp
   real                   , dimension(elements), target     :: vart
   integer                                     , intent(in) :: ng
   character (len=*)                           , intent(in) :: tabstr
   !----- Local variables. ----------------------------------------------------------------!
   character (len=str_len)                                  :: line
   character (len=32)     , dimension(10)                   :: tokens
   character (len=16)                                       :: cname
   character (len=16)                                       :: ctab
   character (len=1)                           , parameter  :: toksep = ':'
   !----- Local variables. ----------------------------------------------------------------!
   integer                                                  :: ntok
   integer                                                  :: nv
   integer                                                  :: isnum
   integer                                                  :: ns
   integer                                                  :: i
   !---------------------------------------------------------------------------------------!


   !----- Split the tokens. ---------------------------------------------------------------!
   call tokenize1(tabstr,tokens,ntok,toksep)
   !---------------------------------------------------------------------------------------!


   !----- Split the tokens. ---------------------------------------------------------------!
   cname = tokens(1)
   !---------------------------------------------------------------------------------------!


   !----- Check whether this scalar name is already in the table. -------------------------!
   isnum = 0
   !---------------------------------------------------------------------------------------!



   nv = num_scalar(ng)
   scalar_tab(nv,ng)%a_var_p => varp(:)
   scalar_tab(nv,ng)%a_var_t => vart(:)

   return
end subroutine vtables_scalar_new
!==========================================================================================!
!==========================================================================================!
