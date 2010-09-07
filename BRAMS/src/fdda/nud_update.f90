!===================================== Change Log =========================================!
! 5.0.0                                                                                    !
!                                                                                          !
!==========================================================================================!
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved                       !
!  Regional Atmospheric Modeling System - RAMS                                             !
!==========================================================================================!
subroutine nud_update(iswap,nnud)
   use var_tables
   use an_header
   use mem_basic
   use mem_grid
   use mem_varinit
   use grid_struct
   use rconstants
   use mem_aerad  ,   only: nwave ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                                      , intent(in) :: iswap
   integer                                      , intent(in) :: nnud
   !----- Local variables. ----------------------------------------------------------------!
   character (len=256) :: hnameinh,prefix
   character (len=2) :: cng
   integer         , dimension(:)  , allocatable             :: nnxp1
   integer         , dimension(:)  , allocatable             :: nnyp1
   integer         , dimension(:)  , allocatable             :: nnzp1
   integer                                                   :: ngrids1
   integer                                                   :: ioutput1
   integer                                                   :: nzg1
   integer                                                   :: nzs1
   integer                                                   :: npatch1
   integer                                                   :: nclouds1
   integer                                                   :: iyr
   integer                                                   :: imn
   integer                                                   :: idy
   integer                                                   :: itm
   integer                                                   :: ie
   integer                                                   :: maxarr
   integer                                                   :: maxarr2
   integer                                                   :: ngr
   integer                                                   :: maxx1
   integer                                                   :: maxy1
   integer                                                   :: maxz1
   integer                                                   :: npts
   integer                                                   :: nptsh
   integer                                                   :: nv
   integer                                                   :: nvh
   integer                                                   :: i
   integer                                                   :: k
   integer                                                   :: nzpg1
   integer                                                   :: nc
   integer                                                   :: ierr
   integer                                                   :: ng
   integer                                                   :: ng_start
   real            , dimension(:)  , allocatable             :: platn1
   real            , dimension(:)  , allocatable             :: plonn1
   real            , dimension(:)  , allocatable             :: deltaxn1
   real            , dimension(:)  , allocatable             :: deltayn1
   real            , dimension(:)  , allocatable             :: scr
   real            , dimension(:,:), allocatable             :: xmn1
   real            , dimension(:,:), allocatable             :: xtn1
   real            , dimension(:,:), allocatable             :: ymn1
   real            , dimension(:,:), allocatable             :: ytn1
   real            , dimension(:,:), allocatable             :: zmn1
   real            , dimension(:,:), allocatable             :: ztn1
   real            , dimension(:,:), allocatable             :: topt1
   real(kind=8)                                              :: time1
   type(grid_def)  , dimension(:)  , allocatable             :: grdefh
   type(grid_def)  , dimension(:)  , allocatable             :: grdefn
   real                                                      :: ztop1
   !----- Locally saved variables. --------------------------------------------------------!
   type(head_table), dimension(:)  , allocatable, save       :: hr_table
   integer                                      , save       :: iunhd  = 11
   integer                                      , save       :: inhunt = 10
   !----- External functions. -------------------------------------------------------------!
   integer                                      , external   :: cio_i
   integer                                      , external   :: cio_f
   integer                                      , external   :: cio_i_sca
   integer                                      , external   :: cio_f_sca
   integer                                      , external   :: cio_f8_sca
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Put new fields into varinit future arrays.  If iswap == 1, swap future into past  !
   ! first.                                                                                !
   !---------------------------------------------------------------------------------------!
   if (iswap == 1) then
      if (nud_type == 1) then
         do ngr=1,ngrids
            npts = nnzp(ngr)*nnxp(ngr)*nnyp(ngr)
            call atob(npts,varinit_g(ngr)%varuf,varinit_g(ngr)%varup)
            call atob(npts,varinit_g(ngr)%varvf,varinit_g(ngr)%varvp)
            call atob(npts,varinit_g(ngr)%varpf,varinit_g(ngr)%varpp)
            call atob(npts,varinit_g(ngr)%vartf,varinit_g(ngr)%vartp)
            call atob(npts,varinit_g(ngr)%varrf,varinit_g(ngr)%varrp)
            if (co2_on) then
               call atob(npts,varinit_g(ngr)%varof,varinit_g(ngr)%varop)
            end if
         end do
      end if
   end if


   !----- Open the input history header file and read some of the info. -------------------!
   nc       = len_trim(fnames_nud(nnud))
   hnameinh = fnames_nud(nnud)(1:nc-9)//'.vfm'
   call rams_f_open(iunhd,fnames_nud(nnud),'FORMATTED','OLD','READ',0)

   ie=cio_i_sca(iunhd,1,'ngrids',ngrids1,1)
   ngridsh=ngrids1

   write(unit=*,fmt='(a,1x,i5)') 'ngrids1:',ngrids1
   !----- Allocate the number of nudging grids. -------------------------------------------!
   allocate (nnxp1(ngrids1),nnyp1(ngrids1),nnzp1(ngrids1))
   allocate (platn1(ngrids1),plonn1(ngrids1))

   ie=cio_i     (iunhd,1,'nnxp'   ,nnxp1   ,ngrids1)
   ie=cio_i     (iunhd,1,'nnyp'   ,nnyp1   ,ngrids1)
   ie=cio_i     (iunhd,1,'nnzp'   ,nnzp1   ,ngrids1)
   ie=cio_i_sca (iunhd,1,'npatch' ,npatch1 ,      1)
   ie=cio_i_sca (iunhd,1,'nclouds',nclouds1,      1)
   ie=cio_i_sca (iunhd,1,'nzg'    ,nzg1    ,      1)
   ie=cio_i_sca (iunhd,1,'nzs'    ,nzs1    ,      1)
   ie=cio_i_sca (iunhd,1,'ioutput',ioutput1,      1)
   ie=cio_f8_sca(iunhd,1,'time'   ,time1   ,      1)
   ie=cio_f_sca (iunhd,1,'ztop'   ,ztop1   ,      1)
   ie=cio_f     (iunhd,1,'platn'  ,platn1  ,ngrids1)
   ie=cio_f     (iunhd,1,'plonn'  ,plonn1  ,ngrids1)

   !---------------------------------------------------------------------------------------!
   !      Find maximum size of any array on history file. Allocate scratch array of this   !
   ! size.                                                                                 !
   !---------------------------------------------------------------------------------------!
   maxarr  = 0
   maxarr2 = 0
   maxx1   = maxval(nnxp1(1:ngrids1))
   maxy1   = maxval(nnyp1(1:ngrids1))
   maxz1   = maxval(nnzp1(1:ngrids1))
   do ngr=1,ngrids1
      maxarr  = max(maxarr,nnxp1(ngr)*nnyp1(ngr)*nnzp1(ngr)                                &
                          ,nnxp1(ngr)*nnyp1(ngr)*nzg1*npatch1                              &
                          ,nnxp1(ngr)*nnyp1(ngr)*nzs1*npatch1                              &
                          ,nnzp1(ngr)*nnxp1(ngr)*nnyp1(ngr)*nclouds1                       &
                          ,nnxp1(ngr)*nnyp1(ngr)*nwave )

      maxarr2 = max(maxarr2,nnxp1(ngr)*nnyp1(ngr))
   enddo

   allocate(xmn1(maxx1,ngrids1),xtn1(maxx1,ngrids1))
   allocate(ymn1(maxy1,ngrids1),ytn1(maxy1,ngrids1))
   allocate(zmn1(maxz1,ngrids1),ztn1(maxz1,ngrids1))

   do ngr=1,ngrids1
      write (cng,fmt='(i2.2)') ngr
      ie=cio_f(iunhd,1,'xmn'//cng,xmn1(:,ngr),nnxp1(ngr))
      ie=cio_f(iunhd,1,'xtn'//cng,xtn1(:,ngr),nnxp1(ngr))
      ie=cio_f(iunhd,1,'ymn'//cng,ymn1(:,ngr),nnyp1(ngr))
      ie=cio_f(iunhd,1,'ytn'//cng,ytn1(:,ngr),nnyp1(ngr))
      ie=cio_f(iunhd,1,'zmn'//cng,zmn1(:,ngr),nnzp1(ngr))
      ie=cio_f(iunhd,1,'ztn'//cng,ztn1(:,ngr),nnzp1(ngr))
   end do

   allocate (topt1(maxarr2,ngrids1))
   allocate (scr(maxarr))
   call rams_f_open(inhunt,hnameinh,'UNFORMATTED','OLD','READ',0)

   !----- Read variable header info. ------------------------------------------------------!
   rewind(unit=iunhd)

   read (unit=iunhd,fmt=*) nvbtab
   allocate (hr_table(nvbtab))
   do nv=1,nvbtab
      read(unit=iunhd,fmt=*)  hr_table(nv)%string,hr_table(nv)%npointer                    &
                             ,hr_table(nv)%idim_type,hr_table(nv)%ngrid                    &
                             ,hr_table(nv)%nvalues
   end do

   !---------------------------------------------------------------------------------------!
   !      Go through file and get all grids' TOPT's.  Don't know how else to get this      !
   ! before processing a field...                                                          !
   !---------------------------------------------------------------------------------------!
   vtableloop: do nvh=1,nvbtab
      ngr=hr_table(nvh)%ngrid
      nptsh=hr_table(nvh)%nvalues
      read(inhunt) (scr(i),i=1,nptsh)

      if (hr_table(nvh)%string == 'TOPTA') then
         call atob(nptsh, scr,topt1(:,ngr))
         if (ngr == ngrids1) exit vtableloop
      end if
   end do vtableloop


   !---------------------------------------------------------------------------------------!
   !      Set a flag array (for each grid on the history file) to determine:               !
   !   >= 1 = This grid is identical to a current grid;                                    !
   !     -1 = This grid is different.                                                      !
   !---------------------------------------------------------------------------------------!

   igrid_match(1:ngrids1)=0

   !----- Allocate grid structures. -------------------------------------------------------!
   allocate (grdefn(ngrids))
   allocate (grdefh(ngrids1))
   write(unit=*,fmt=*) '1============:', igrid_match(1:ngrids)

   do ngr=1,ngrids1
      call alloc_grid_def(grdefh(ngr),nnxp1(ngr),nnyp1(ngr),nnzp1(ngr))
      call  fill_grid_def(grdefh(ngr),nnxp1(ngr),nnyp1(ngr),nnzp1(ngr),nzg1,nzs1,npatch1   &
                         ,platn1(1),plonn1(1),xtn1(:,ngr),xmn1(:,ngr),ytn1(:,ngr)          &
                         ,ymn1(:,ngr),ztn1(:,ngr),zmn1(:,ngr),topt1(:,ngr))
   end do

   do ngr=1,ngrids
      call alloc_grid_def(grdefn(ngr),nnxp(ngr),nnyp(ngr),nnzp(ngr))
      call  fill_grid_def(grdefn(ngr),nnxp(ngr),nnyp(ngr),nnzp(ngr),nzg,nzs,npatch         &
                         ,platn(1),plonn(1),xtn(:,ngr),xmn(:,ngr),ytn(:,ngr),ymn(:,ngr)    &
                         ,ztn(:,ngr),zmn(:,ngr),grid_g(ngr)%topta)
   end do
   write(unit=*,fmt=*) '2============:',igrid_match(1:ngrids),ngrids,ngrids1

   !---------------------------------------------------------------------------------------!
   !      See if the history grids match any of the new grids...assuming 1:1 grid number   !
   ! correspondence for now.                                                               !
   !---------------------------------------------------------------------------------------!
   do ngr=1,min(ngrids,ngrids1)
      call compare_grid_def(grdefh(ngr),grdefn(ngr),'nud_update',ierr)
      if (ierr /= 0) then
         !----- No match... ---------------------------------------------------------------!
         igrid_match(ngr) = -1
      else
         !----- We have a match... --------------------------------------------------------!
         igrid_match(ngr) = ngr
      end if
   end do

   write(unit=*,fmt=*) '3============:', igrid_match(1:ngrids)

   !----- Finally, process the fields... --------------------------------------------------!


   if (nud_type == 1) then

      rewind (unit=inhunt)

      write(unit=*,fmt=*) 'matches:',igrid_match(1:ngrids1)

      !------------------------------------------------------------------------------------!
      !!!!!!!!!!!!!!!!!!!!!  NEED WIND ROTATION FOR THE GENERAL CASE!! !!!!!!!!!!!!!!!!!!!!!
      !------------------------------------------------------------------------------------!
      !      Loop through all variables - "normal" nudging.  This will interpolate from    !
      ! the coarsest grid first.  In case a point is on a nested grid, the value will be   !
      ! overwritten.                                                                       !
      !------------------------------------------------------------------------------------!

      read_loop: do nvh=1,nvbtab
         !----- Read a variable. ----------------------------------------------------------!
         nptsh=hr_table(nvh)%nvalues
         read(inhunt)(scr(i),i=1,nptsh)

         ngr=hr_table(nvh)%ngrid


         !---------------------------------------------------------------------------------!
         !     Okay, this gets complicated... there are a number of possibilities here,    !
         ! depending on if the new grids match the history file.                           !
         !     If the grids do not match, then we want to interpolate each new grid from   !
         ! every history file field, so as to get the highest resolution info on to each   !
         ! grid. This is the easy part                                                     !
         !     If some grids match, we will fill the fields directly into the correspond-  !
         ! ing grid.  But we still will need to interpolate non-matching grids.  So for    !
         ! each history field, we will check to see if the prior grid had matched.         !
         !     Here's some examples: we have a new 2 grid run using a history file that    !
         ! was only run with one grid, and grid 1 matches the new run. The grid 1 field be !
         ! can filled directly, therefore we don't need to interpolate the new grid 1. The !
         ! new grid 2, however, needs to be handled, so we interpolate it. The new grid    !
         ! loop needs to run through all grids.                                            !
         !     Let's now assume we have a new 3 grid run from a 2 grid history file and    !
         ! grids 1 and 2 match.  When we read a grid 1 history field, we will fill it      !
         ! directly.  Then we will interpolate grid 2 and 3 from the grid 1 field.  Then   !
         ! we read a grid 2 field.  We don't want to do anything with new grid 1: it had   !
         ! already been filled.  We want to start from the new grid 2, fill it, then       !
         ! interpolate grid 3.  So we start from the new grid 2 on the "grid_loop".        !
         !     I think the following logic will work, but IT ASSUMES THAT MATCHING GRIDS   !
         ! HAVE A GRID NUMBER CORRESPONDENCE WITH THE HISTORY FILE AND THE NEW RUN:        !
         ! history grid 1 matches new grid 1, etc.                                         !
         !---------------------------------------------------------------------------------!

         !----- Ngr is zero sometimes, that is the reason for the if (ngr >1)... ----------!
         ng_start=1
         if (ngr > 1) then
           if (igrid_match(ngr-1) > 0) ng_start = ngr
         end if

         grid_loop: do ng = ng_start, ngrids
            !----- Find which is the corresponding variable in the current run. -----------!
            var_loop: do nv = 1,num_var(ng)
               npts=vtab_r(nv,ng)%npts
               !---------------------------------------------------------------------------!
               !      See if this variable is active in the current run, but only interpo- !
               ! late if UP,VP,THP,RTP,PP,CO2P.  Then, for each of these variables, we     !
               ! must check whether we have a match, and if so, simply fill history field  !
               ! into array.  Otherwise, we call the interpolation procedure.              !
               !---------------------------------------------------------------------------!
               if (hr_table(nvh)%string == vtab_r(nv,ng)%name) then
                  select case (trim(vtab_r(nv,ng)%name))
                  case ('UP') !----- Zonal wind. ------------------------------------------!
                     if (igrid_match(ngr) == ng) then
                        write (unit=*,fmt='(a30,2i5,3x,a8,i8)')                            &
                              'NUD_UPDATE: Filling:       ',ngr,ng,vtab_r(nv,ngr)%name,npts
                        call atob(nptsh,scr,varinit_g(ng)%varuf)

                     else
                        write (unit=*,fmt='(a30,2i5,3x,a8,i8)')                            &
                              'NUD_UPDATE: Interpolating: ',ngr,ng,vtab_r(nv,ngr)%name,npts

                        call hi_interp(nnzp1(ngr),nnxp1(ngr),nnyp1(ngr),1,scr,xmn1(:,ngr)  &
                                      ,xtn1(:,ngr),ymn1(:,ngr),ytn1(:,ngr),zmn1(:,ngr)     &
                                      ,ztn1(:,ngr),platn1(ngr),plonn1(ngr),topt1(:,ngr)    &
                                      ,ztop1,nnzp(ng),nnxp(ng),nnyp(ng),1                  &
                                      ,varinit_g(ng)%varuf,ng,ngr,vtab_r(nv,ng)%name,3)
                     end if
                     cycle grid_loop

                 case ('VP') !----- Meridional wind. --------------------------------------!
                     if (igrid_match(ngr) == ng) then
                        write (unit=*,fmt='(a30,2i5,3x,a8,i8)')                            &
                              'NUD_UPDATE: Filling:       ',ngr,ng,vtab_r(nv,ngr)%name,npts
                        call atob(nptsh,scr(1),varinit_g(ng)%varvf)
                     else
                        write (unit=*,fmt='(a30,2i5,3x,a8,i8)')                            &
                              'NUD_UPDATE: Interpolating: ',ngr,ng,vtab_r(nv,ngr)%name,npts

                        call hi_interp(nnzp1(ngr),nnxp1(ngr),nnyp1(ngr),1,scr,xmn1(:,ngr)  &
                                      ,xtn1(:,ngr),ymn1(:,ngr),ytn1(:,ngr),zmn1(:,ngr)     &
                                      ,ztn1(:,ngr),platn1(ngr),plonn1(ngr),topt1(:,ngr)    &
                                      ,ztop1,nnzp(ng),nnxp(ng),nnyp(ng),1                  &
                                      ,varinit_g(ng)%varvf,ng,ngr,vtab_r(nv,ngr)%name,3)
                     end if
                     cycle grid_loop

                  case ('THP') !----- ("Ice-liquid") Potential temperature. ---------------!
                     if (igrid_match(ngr) == ng) then
                        write (unit=*,fmt='(a30,2i5,3x,a8,i8)')                            &
                              'NUD_UPDATE: Filling:       ',ngr,ng,vtab_r(nv,ngr)%name,npts
                        call atob(nptsh,scr(1),varinit_g(ng)%vartf)
                     else
                        write (unit=*,fmt='(a30,2i5,3x,a8,i8)')                            &
                              'NUD_UPDATE: Interpolating: ',ngr,ng,vtab_r(nv,ngr)%name,npts

                        call hi_interp(nnzp1(ngr),nnxp1(ngr),nnyp1(ngr),1,scr,xmn1(:,ngr)  &
                                      ,xtn1(:,ngr),ymn1(:,ngr),ytn1(:,ngr),zmn1(:,ngr)     &
                                      ,ztn1(:,ngr),platn1(ngr),plonn1(ngr),topt1(:,ngr)    &
                                      ,ztop1,nnzp(ng),nnxp(ng),nnyp(ng),1                  &
                                      ,varinit_g(ng)%vartf,ng,ngr,vtab_r(nv,ng)%name,3)
                     end if
                     cycle grid_loop

                  case ('RTP') !----- Total water mixing ratio. ------------------------------!
                     if (igrid_match(ngr) == ng) then
                        write (unit=*,fmt='(a30,2i5,3x,a8,i8)')                            &
                              'NUD_UPDATE: Filling:       ',ngr,ng,vtab_r(nv,ngr)%name,npts
                        call atob(nptsh,scr,varinit_g(ng)%varrf)

                     else
                        write (unit=*,fmt='(a30,2i5,3x,a8,i8)')                            &
                              'NUD_UPDATE: Interpolating: ',ngr,ng,vtab_r(nv,ngr)%name,npts
                        call hi_interp(nnzp1(ngr),nnxp1(ngr),nnyp1(ngr),1,scr,xmn1(:,ngr)  &
                                      ,xtn1(:,ngr),ymn1(:,ngr),ytn1(:,ngr),zmn1(:,ngr)     &
                                      ,ztn1(:,ngr),platn1(ngr),plonn1(ngr),topt1(:,ngr)    &
                                      ,ztop1,nnzp(ng),nnxp(ng),nnyp(ng),1                  &
                                      ,varinit_g(ng)%varrf,ng,ngr,vtab_r(nv,ng)%name,3)
                     end if
                     cycle grid_loop

                  case ('PP') !----- (Perturbation of) Exner function. -----------------------!
                     if ( igrid_match(ngr) == ng ) then
                        write (unit=*,fmt='(a30,2i5,3x,a8,i8)')                            &
                              'NUD_UPDATE: Filling:       ',ngr,ng,vtab_r(nv,ngr)%name,npts
                        call atob(nptsh,scr,varinit_g(ng)%varpf)

                     else
                        write (unit=*,fmt='(a30,2i5,3x,a8,i8)')                            &
                              'NUD_UPDATE: Interpolating: ',ngr,ng,vtab_r(nv,ngr)%name,npts
                        call hi_interp(nnzp1(ngr),nnxp1(ngr),nnyp1(ngr),1,scr,xmn1(:,ngr)  &
                                      ,xtn1(:,ngr),ymn1(:,ngr),ytn1(:,ngr),zmn1(:,ngr)     &
                                      ,ztn1(:,ngr),platn1(ngr),plonn1(ngr),topt1(:,ngr)    &
                                      ,ztop1,nnzp(ng),nnxp(ng),nnyp(ng),1                  &
                                      ,varinit_g(ng)%varpf,ng,ngr,vtab_r(nv,ngr)%name,3)
                     end if
                     cycle grid_loop

                  case ('CO2P') !----- CO2 mixing ratio. -------------------------------------!
                     if (igrid_match(ngr) == ng) then
                        write (unit=*,fmt='(a30,2i5,3x,a8,i8)')                            &
                              'NUD_UPDATE: Filling:       ',ngr,ng,vtab_r(nv,ngr)%name,npts
                        call atob(nptsh,scr,varinit_g(ng)%varof)

                     else
                        write (unit=*,fmt='(a30,2i5,3x,a8,i8)')                            &
                              'NUD_UPDATE: Interpolating: ',ngr,ng,vtab_r(nv,ngr)%name,npts
                        call hi_interp(nnzp1(ngr),nnxp1(ngr),nnyp1(ngr),1,scr,xmn1(:,ngr)  &
                                      ,xtn1(:,ngr),ymn1(:,ngr),ytn1(:,ngr),zmn1(:,ngr)     &
                                      ,ztn1(:,ngr),platn1(ngr),plonn1(ngr),topt1(:,ngr)    &
                                      ,ztop1,nnzp(ng),nnxp(ng),nnyp(ng),1                  &
                                      ,varinit_g(ng)%varof,ng,ngr,vtab_r(nv,ng)%name,3)
                     end if
                     cycle grid_loop
                  end select
               end if

            end do var_loop
         end do grid_loop
      end do read_loop
   end if

   !----- Close the input history file. ---------------------------------------------------!
   close(unit=inhunt,status='keep')
   close(unit=iunhd,status='keep')

   deallocate(scr,hr_table)
   deallocate(xmn1,xtn1,ymn1,ytn1,zmn1,ztn1)

   return
end subroutine nud_update
!==========================================================================================!
!==========================================================================================!
