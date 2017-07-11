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
subroutine recycle()

   use mem_grid
   use mem_leaf
   use mem_scratch
   use var_tables
   use io_params
   use grid_dims  , only : str_len ! ! intent(in)
   use mem_aerad  , only : nwave   ! ! intent(in)
   use mem_cuparm , only : nclouds ! ! intent(in)

   implicit none
   !------ Local variables. ---------------------------------------------------------------!
   character(len=str_len)           :: flnm
   integer                          :: ng
   integer                          :: nvars
   integer                          :: lenf
   integer                          :: ierr
   integer                          :: nvert
   integer                          :: nenvi
   !------ External functions. ------------------------------------------------------------!
   integer               , external :: RAMS_getvar
   !---------------------------------------------------------------------------------------!

   flnm=pastfn(1:len_trim(pastfn)-9)

   write(unit=*,fmt='(2(a,1x))') 'Reading assimilation fields from analysis file '         &
                                ,trim(pastfn)

   call rams_read_header(trim(flnm))

   !----- Read the requested analysis file variables. -------------------------------------!
   gridloop: do ng=1,ngrids
      varloop: do nvars=1,num_var(ng)
      
         if (vtab_r(nvars,ng)%irecycle == 1) then

            write (unit=*,fmt='(2(a,1x),3(a,1x,i6,1x))')                                   &
                                      'Reading assimilation field:',vtab_r(nvars,ng)%name  &
                                     ,' for grid:',ng                                      &
                                     ,' dim:',vtab_r(nvars,ng)%idim_type                   &
                                     ,' npts:', vtab_r(nvars,ng)%npts

            ierr=RAMS_getvar(vtab_r(nvars,ng)%name,ng,scratch%scr1,scratch%scr2,flnm)

            !----- Find the vertical and environment dimensions. --------------------------!
            call ze_dims(ng,vtab_r(nvars,ng)%idim_type,.true.,nvert,nenvi)
            
            !----- Copy the result from scratch to the pointer. ---------------------------!
            call unarrange(nvert,nnxp(ng),nnyp(ng),nenvi,scratch%scr1                      &
                          ,vtab_r(nvars,ng)%var_p)
         end if
      end do varloop
   end do gridloop
   !---------------------------------------------------------------------------------------!

   return
end subroutine recycle
!==========================================================================================!
!==========================================================================================!
