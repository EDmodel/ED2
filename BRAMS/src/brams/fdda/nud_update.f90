!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################


subroutine nud_update(iswap,nnud)

use var_tables
use an_header
use mem_basic
use mem_grid
use mem_varinit
use grid_struct
use rconstants

! For use with CATT
use mem_aerad,   only: nwave      !intent(in)


implicit none

integer :: iswap,nnud

integer :: ngrids1,ioutput1,nzg1,nzs1,npatch1
real :: ztop1
real(kind=8) :: time1
integer, allocatable, dimension(:) :: nnxp1,nnyp1,nnzp1
real, allocatable, dimension(:) :: platn1,plonn1,deltaxn1,deltayn1
real, allocatable, dimension(:,:) :: xmn1,xtn1,ymn1,ytn1,zmn1,ztn1
real, allocatable, dimension(:,:) :: topt1


integer :: iyr,imn,idy,itm,ie,maxarr,maxarr2,ngr,maxx1,maxy1,maxz1
character (len=80) :: hnameinh,prefix
character (len=2) :: cng
integer, external :: cio_i,cio_f,cio_i_sca,cio_f_sca,cio_f8_sca
integer,save :: iunhd=11,inhunt=10

integer :: npts,nptsh,nv,nvh,i,k,nzpg1,nc,ierr,ng,ng_start
real, allocatable :: scr(:)
real :: t1,w1,cputime

type (head_table), allocatable,save :: hr_table(:)

type(grid_def), allocatable :: grdefh(:)
type(grid_def), allocatable :: grdefn(:)

! Put new fields into varinit future arrays. If iswap == 1,
!     swap future into past first

if (iswap == 1) then
   if (nud_type == 1) then
      do ngr=1,ngrids
         varinit_g(ngr)%varup(1:nnzp(ngr),1:nnxp(ngr),1:nnyp(ngr))=  &
            varinit_g(ngr)%varuf(1:nnzp(ngr),1:nnxp(ngr),1:nnyp(ngr))
         varinit_g(ngr)%varvp(1:nnzp(ngr),1:nnxp(ngr),1:nnyp(ngr))=  &
            varinit_g(ngr)%varvf(1:nnzp(ngr),1:nnxp(ngr),1:nnyp(ngr))
         varinit_g(ngr)%varpp(1:nnzp(ngr),1:nnxp(ngr),1:nnyp(ngr))=  &
            varinit_g(ngr)%varpf(1:nnzp(ngr),1:nnxp(ngr),1:nnyp(ngr))
         varinit_g(ngr)%vartp(1:nnzp(ngr),1:nnxp(ngr),1:nnyp(ngr))=  &
            varinit_g(ngr)%vartf(1:nnzp(ngr),1:nnxp(ngr),1:nnyp(ngr))
         varinit_g(ngr)%varrp(1:nnzp(ngr),1:nnxp(ngr),1:nnyp(ngr))=  &
            varinit_g(ngr)%varrf(1:nnzp(ngr),1:nnxp(ngr),1:nnyp(ngr))
      enddo
   endif
endif


! Open the input history header file and read some of the info.

nc=len_trim(fnames_nud(nnud))
hnameinh=fnames_nud(nnud)(1:nc-9)//'.vfm'

call rams_f_open(iunhd,fnames_nud(nnud),'FORMATTED','OLD','READ',0)

ie=cio_i_sca(iunhd,1,'ngrids',ngrids1,1)
ngridsh=ngrids1

print*,'ngrids1:',ngrids1
allocate (nnxp1(ngrids1),nnyp1(ngrids1),nnzp1(ngrids1))
allocate (platn1(ngrids1),plonn1(ngrids1))

ie=cio_i(iunhd,1,'nnxp',nnxp1,ngrids1)
ie=cio_i(iunhd,1,'nnyp',nnyp1,ngrids1)
ie=cio_i(iunhd,1,'nnzp',nnzp1,ngrids1)
ie=cio_i_sca(iunhd,1,'npatch',npatch1,1)
ie=cio_i_sca(iunhd,1,'nzg',nzg1,1)
ie=cio_i_sca(iunhd,1,'nzs',nzs1,1)
ie=cio_i_sca(iunhd,1,'ioutput',ioutput1,1)
ie=cio_f8_sca(iunhd,1,'time',time1,1)
ie=cio_f_sca(iunhd,1,'ztop',ztop1,1)
ie=cio_f(iunhd,1,'platn',platn1,ngrids1)
ie=cio_f(iunhd,1,'plonn',plonn1,ngrids1)

! Find maximum size of any array on history file. Allocate scratch array of
! this size.

maxarr=0
maxarr2=0
maxx1=maxval(nnxp1(1:ngrids1))
maxy1=maxval(nnyp1(1:ngrids1))
maxz1=maxval(nnzp1(1:ngrids1))
do ngr=1,ngrids1
   maxarr=max(maxarr,nnxp1(ngr)*nnyp1(ngr)*nnzp1(ngr)  &
         ,nnxp1(ngr)*nnyp1(ngr)*nzg1*npatch1 &
         ,nnxp1(ngr)*nnyp1(ngr)*nzs1*npatch1)

   !For use with CATT
   maxarr=max(maxarr,nnxp1(ngr)*nnyp1(ngr)*nwave)

   maxarr2=max(maxarr2,nnxp1(ngr)*nnyp1(ngr))
enddo

allocate(xmn1(maxx1,ngrids1),xtn1(maxx1,ngrids1))
allocate(ymn1(maxy1,ngrids1),ytn1(maxy1,ngrids1))
allocate(zmn1(maxz1,ngrids1),ztn1(maxz1,ngrids1))

do ngr=1,ngrids1
   write(cng,'(i2.2)') ngr
   ie=cio_f(iunhd,1,'xmn'//cng,xmn1(1,ngr),nnxp1(ngr))
   ie=cio_f(iunhd,1,'xtn'//cng,xtn1(1,ngr),nnxp1(ngr))
   ie=cio_f(iunhd,1,'ymn'//cng,ymn1(1,ngr),nnyp1(ngr))
   ie=cio_f(iunhd,1,'ytn'//cng,ytn1(1,ngr),nnyp1(ngr))
   ie=cio_f(iunhd,1,'zmn'//cng,zmn1(1,ngr),nnzp1(ngr))
   ie=cio_f(iunhd,1,'ztn'//cng,ztn1(1,ngr),nnzp1(ngr))
enddo

allocate (topt1(maxarr2,ngrids1))

allocate (scr(maxarr))

call rams_f_open(inhunt,hnameinh,'UNFORMATTED','OLD','READ',0)

!  Read variable header info

rewind(iunhd)

read(iunhd,*) nvbtab
allocate (hr_table(nvbtab))
do nv=1,nvbtab
   read(iunhd,*)  hr_table(nv)%string   &
                 ,hr_table(nv)%npointer  &
                 ,hr_table(nv)%idim_type  &
                 ,hr_table(nv)%ngrid  &
                 ,hr_table(nv)%nvalues
enddo

! Go through file and get all grids' TOPT's
!   Don't know how else to get this before processing a field...
do nvh=1,nvbtab
   ngr=hr_table(nvh)%ngrid
   nptsh=hr_table(nvh)%nvalues

   read(inhunt)(scr(i),i=1,nptsh)

   if (hr_table(nvh)%string == 'TOPTA') then
      call atob(nptsh, scr(1),topt1(1,ngr))
      if (ngr == ngrids1) exit
   endif
enddo


! Set a flag array (for each grid on the history file) to determine:
!   >= 1 = This grid is identical to a current grid
!     -1 = This grid is different.

igrid_match(1:ngrids1)=0

! Allocate grid structures
allocate (grdefn(ngrids))
allocate (grdefh(ngrids1))
print*,'1============:', igrid_match(1:ngrids)

do ngr=1,ngrids1
   call alloc_grid_def( grdefh(ngr),nnxp1(ngr),nnyp1(ngr),nnzp1(ngr) )
   call  fill_grid_def( grdefh(ngr),nnxp1(ngr),nnyp1(ngr),nnzp1(ngr) &
                       ,nzg1,nzs1,npatch1,platn1(1),plonn1(1)  &
                       ,xtn1(1,ngr),xmn1(1,ngr),ytn1(1,ngr),ymn1(1,ngr)  &
                       ,ztn1(1,ngr),zmn1(1,ngr),topt1(1,ngr) )
enddo
do ngr=1,ngrids
   call alloc_grid_def( grdefn(ngr),nnxp(ngr),nnyp(ngr),nnzp(ngr) )
   call  fill_grid_def( grdefn(ngr),nnxp(ngr),nnyp(ngr),nnzp(ngr) &
                       ,nzg,nzs,npatch,platn(1),plonn(1)  &
                       ,xtn(1,ngr),xmn(1,ngr),ytn(1,ngr),ymn(1,ngr)  &
                       ,ztn(1,ngr),zmn(1,ngr),grid_g(ngr)%topta )
enddo
print*,'2============:', igrid_match(1:ngrids),ngrids,ngrids1

! See if the history grids match any of the new grids...assuming 1:1 grid number
!   correspondence for now

do ngr=1,min(ngrids,ngrids1)
   call compare_grid_def(grdefh(ngr),grdefn(ngr),'nud_update',ierr)
   if (ierr /= 0) then
      ! No match...
      igrid_match(ngr)=-1
   else
      ! We have a match...
      igrid_match(ngr)=ngr
   endif
enddo

print*,'3============:', igrid_match(1:ngrids)

! Finally, process the fields...


if (nud_type == 1) then

rewind inhunt


print*,'matches:',igrid_match(1:ngrids1)

!!!!!!!!!!!!!!!!  Need wind rotation for the general case!!!!!!!!!!!!!!!!!!!!!!

! Loop through all variables - "normal" nudging

! This will interpolate from the coarsest grid first. If a point is on a nested
!   grid, the value will be overwritten.

read_loop: do nvh=1,nvbtab
   ! Read a variable
   nptsh=hr_table(nvh)%nvalues
   read(inhunt)(scr(i),i=1,nptsh)

   ngr=hr_table(nvh)%ngrid


   ! Okay, this gets complicated... there are a number of possibilities here,
   !   depending on if the new grids match the history file.

   !   If the grids do not match, then we want to interpolate each new grid
   !   from every history file field, so as to get the highest resolution info
   !   on to each grid. This is the easy part!

   !   If some grids match, we will fill the fields directly into the
   !   corresponding grid. But we still will need to interpolate non-matching
   !   grids. So for each history field, we will check to see if the prior grid
   !   matched.
   !   Here's some examples: we have a new 2 grid run using a history file that was
   !   only run with one grid, and grid 1 matches the new run. The grid 1 field
   !   can be filled directly, therefore we don't need to interpolate the new grid 1.
   !   But new grid 2 needs to be handled, so we interpolate it. The new grid loop
   !   needs to run through all grids.
   !   Let's now assume we have a new 3 grid run from a 2 grid history file and
   !   grids 1 and 2 match. When we read a grid 1 history field, we will fill it
   !   directly. Then we will interpolate grid 2 and 3 from the grid 1 field.
   !   Then we read a grid 2 field. We don't want to do anything with new grid 1;
   !   it was already filled. We want to start from the new grid 2, fill it, then
   !   interpolate grid 3. So we start from the new grid 2 on the "grid_loop".

   !   I think the following logic will work, but IT ASSUMES THAT MATCHING GRIDS
   !   HAVE A GRID NUMBER CORRESPONDENCE WITH THE HISTORY FILE AND THE NEW RUN:
   !   history grid 1 matches new grid 1, etc.

   ng_start=1
   ![MLO - Had to do this break to avoid check bounds error...
   if (ngr > 1) then
     if (igrid_match(ngr-1) > 0) ng_start = ngr
   end if
   grid_loop:   do ng = ng_start, ngrids


      ! Find which is the corresponding variable in the current run

      do nv = 1,num_var(ng)
         npts=vtab_r(nv,ng)%npts

         !  See if this variable is active in the current run,
         !      but only interpolate if UP,VP,THP,RTP,PP
         if(hr_table(nvh)%string == vtab_r(nv,ng)%name) then

            if (vtab_r(nv,ng)%name == 'UP') then
               ! If we have a match, simply fill history field into array
               if ( igrid_match(ngr) == ng ) then
                  print 33,'nud_update: filling: ',ngr &
                    ,ng,vtab_r(nv,ngr)%name,npts
                    33 format(a30,2i5,3x,a8,i8)
                  call atob(nptsh,scr(1),varinit_g(ng)%varuf(1,1,1))
               else
                  ! Otherwise, interpolate...
                  print 33,'nud_update: interpolating: ',ngr &
                       ,ng,vtab_r(nv,ng)%name,npts
                  t1=cputime(w1)
                  call hi_interp(nnzp1(ngr),nnxp1(ngr),nnyp1(ngr),1,scr(1)  &
                                   ,xmn1(1,ngr),xtn1(1,ngr)  &
                                   ,ymn1(1,ngr),ytn1(1,ngr)  &
                                   ,zmn1(1,ngr),ztn1(1,ngr)  &
                                   ,platn1(ngr),plonn1(ngr)  &
                                   ,topt1(1,ngr),ztop1  &
                                   ,nnzp(ng),nnxp(ng),nnyp(ng),1  &
                                   ,varinit_g(ng)%varuf(1,1,1)  &
                                   ,ng,ngr,vtab_r(nv,ng)%name,3)
               endif
               cycle grid_loop

            elseif (vtab_r(nv,ng)%name == 'VP') then
               if ( igrid_match(ngr) == ng ) then
                  print 33,'nud_update: filling: ',ngr &
                    ,ng,vtab_r(nv,ng)%name,npts
                  call atob(nptsh,scr(1),varinit_g(ng)%varvf(1,1,1))
               else
                  print 33,'nud_update: interpolating: ',ngr &
                       ,ng,vtab_r(nv,ng)%name,npts
                  call hi_interp(nnzp1(ngr),nnxp1(ngr),nnyp1(ngr),1,scr(1)  &
                                   ,xmn1(1,ngr),xtn1(1,ngr)  &
                                   ,ymn1(1,ngr),ytn1(1,ngr)  &
                                   ,zmn1(1,ngr),ztn1(1,ngr)  &
                                   ,platn1(ngr),plonn1(ngr)  &
                                   ,topt1(1,ngr),ztop1  &
                                   ,nnzp(ng),nnxp(ng),nnyp(ng),1  &
                                   ,varinit_g(ng)%varvf(1,1,1)  &
                                   ,ng,ngr,vtab_r(nv,ngr)%name,3)
               endif
               cycle grid_loop

            elseif (vtab_r(nv,ng)%name == 'THP') then
               if ( igrid_match(ngr) == ng ) then
                  print 33,'nud_update: filling: ',ngr &
                    ,ng,vtab_r(nv,ng)%name,npts
                  call atob(nptsh,scr(1),varinit_g(ng)%vartf(1,1,1))
               else
                  print 33,'nud_update: interpolating: ',ngr &
                       ,ng,vtab_r(nv,ng)%name,npts
                  call hi_interp(nnzp1(ngr),nnxp1(ngr),nnyp1(ngr),1,scr(1)  &
                                   ,xmn1(1,ngr),xtn1(1,ngr)  &
                                   ,ymn1(1,ngr),ytn1(1,ngr)  &
                                   ,zmn1(1,ngr),ztn1(1,ngr)  &
                                   ,platn1(ngr),plonn1(ngr)  &
                                   ,topt1(1,ngr),ztop1  &
                                   ,nnzp(ng),nnxp(ng),nnyp(ng),1  &
                                   ,varinit_g(ng)%vartf(1,1,1)  &
                                   ,ng,ngr,vtab_r(nv,ng)%name,3)
               endif
               cycle grid_loop

            elseif (vtab_r(nv,ng)%name == 'RTP') then
               if ( igrid_match(ngr) == ng ) then
                  print 33,'nud_update: filling: ',ngr &
                    ,ng,vtab_r(nv,ng)%name,npts
                  call atob(nptsh,scr(1),varinit_g(ng)%varrf(1,1,1))
               else
                  print 33,'nud_update: interpolating: ',ngr &
                       ,ng,vtab_r(nv,ng)%name,npts
                  call hi_interp(nnzp1(ngr),nnxp1(ngr),nnyp1(ngr),1,scr(1)  &
                                   ,xmn1(1,ngr),xtn1(1,ngr)  &
                                   ,ymn1(1,ngr),ytn1(1,ngr)  &
                                   ,zmn1(1,ngr),ztn1(1,ngr)  &
                                   ,platn1(ngr),plonn1(ngr)  &
                                   ,topt1(1,ngr),ztop1  &
                                   ,nnzp(ng),nnxp(ng),nnyp(ng),1  &
                                   ,varinit_g(ng)%varrf(1,1,1)  &
                                   ,ng,ngr,vtab_r(nv,ng)%name,3)
               endif
               cycle grid_loop

            elseif (vtab_r(nv,ng)%name == 'PP') then
               if ( igrid_match(ngr) == ng ) then
                  print 33,'nud_update: filling: ',ngr &
                    ,ng,vtab_r(nv,ng)%name,npts
                  call atob(nptsh,scr(1),varinit_g(ng)%varpf(1,1,1))
               else
                  print 33,'nud_update: interpolating: ',ngr &
                       ,ng,vtab_r(nv,ng)%name,npts
                  call hi_interp(nnzp1(ngr),nnxp1(ngr),nnyp1(ngr),1,scr(1)  &
                                   ,xmn1(1,ngr),xtn1(1,ngr)  &
                                   ,ymn1(1,ngr),ytn1(1,ngr)  &
                                   ,zmn1(1,ngr),ztn1(1,ngr)  &
                                   ,platn1(ngr),plonn1(ngr)  &
                                   ,topt1(1,ngr),ztop1  &
                                   ,nnzp(ng),nnxp(ng),nnyp(ng),1  &
                                   ,varinit_g(ng)%varpf(1,1,1)  &
                                   ,ng,ngr,vtab_r(nv,ngr)%name,3)
               endif
               cycle grid_loop

            endif

         endif
      enddo

   enddo grid_loop

enddo read_loop

endif

! Close the input history file

close(inhunt)
close(iunhd)

deallocate(scr,hr_table)
deallocate(xmn1,xtn1,ymn1,ytn1,zmn1,ztn1)

return
end

!-------------------

integer function check_real(xx, x, nx)

! Check two corresponding real arrays and see values are close enough

implicit none

integer :: nx,i
real :: xx(nx), x(nx), tol

tol = min( (maxval(xx)-minval(xx)),(maxval(x)-minval(x)) )*.0001

do i = 1, nx
   if (abs(xx(i)-x(i)) > tol) then
      check_real=i
      return
   endif
enddo

check_real = 0

return
end
