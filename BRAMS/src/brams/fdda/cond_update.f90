!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

subroutine cond_update(iswap,ncond)

use var_tables
use an_header
use mem_basic
use mem_grid
use mem_varinit
use grid_struct
use rconstants

implicit none

integer :: iswap,ncond

integer :: ngrids1,ioutput1,nzg1,nzs1,npatch1
real :: ztop1
integer, allocatable, dimension(:) :: nnxp1,nnyp1,nnzp1
real, allocatable, dimension(:) :: platn1,plonn1,deltaxn1,deltayn1
real, allocatable, dimension(:,:) :: xmn1,xtn1,ymn1,ytn1,zmn1,ztn1
real, allocatable, dimension(:,:) :: topt1
real(kind=8) :: time1

integer :: iyr,imn,idy,itm,ie,maxarr,maxarr2,ngr,maxx1,maxy1,maxz1
character (len=80) :: hnameinh,prefix
character (len=2) :: cng
integer, external :: cio_i,cio_f,cio_f8_sca,cio_i_sca,cio_f_sca
integer,save :: iunhd=11,inhunt=10

integer :: npts,nptsh,nv,nvh,i,k,nzpg1,nc,ierr,ng,ng_start
real, allocatable :: scr(:),scr2(:)
real :: t1,w1,cputime

type (head_table), allocatable,save :: hr_table(:)

type(grid_def), allocatable :: grdefh(:)
type(grid_def), allocatable :: grdefn(:)

! Put new fields into varinit future arrays. If iswap == 1, 
!     swap future into past first

if (iswap == 1) then
   do ngr=1,ngrids
      varinit_g(ngr)%varrph(1:nnzp(ngr),1:nnxp(ngr),1:nnyp(ngr))=  &
         varinit_g(ngr)%varrfh(1:nnzp(ngr),1:nnxp(ngr),1:nnyp(ngr))
      varinit_g(ngr)%varcph(1:nnzp(ngr),1:nnxp(ngr),1:nnyp(ngr))=  &
         varinit_g(ngr)%varcfh(1:nnzp(ngr),1:nnxp(ngr),1:nnyp(ngr))
   enddo
endif

! Initialize condensate array

do ngr=1,ngrids
   varinit_g(ngr)%varcfh(1:nnzp(ngr),1:nnxp(ngr),1:nnyp(ngr)) = 0.
enddo


! Open the input history header file and read some of the info.
print*,'ncond:',ncond,fnames_cond(ncond)
nc=len_trim(fnames_cond(ncond))
hnameinh=fnames_cond(ncond)(1:nc-9)//'.vfm'

call rams_f_open(iunhd,fnames_cond(ncond),'FORMATTED','OLD','READ',0)

ie=cio_i_sca(iunhd,1,'ngrids',ngrids1,1)
ngridsh=ngrids1


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

allocate (scr(maxarr),topt1(maxarr2,ngrids1))

! Find maximum size of any array in current run. Allocate scratch array of
! this size.

maxarr=0
do ngr=1,ngrids
   maxarr=max(maxarr,nnxp(ngr)*nnyp(ngr)*nnzp(ngr)  &
         ,nnxp(ngr)*nnyp(ngr)*nzg*npatch &
         ,nnxp(ngr)*nnyp(ngr)*nzs*npatch)
enddo

allocate (scr2(maxarr))

! Open history data file 

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
      call atob(nptsh, scr,topt1(:,ngr))
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

do ngr=1,ngrids1
   call alloc_grid_def( grdefh(ngr),nnxp1(ngr),nnyp1(ngr),nnzp1(ngr) )
   call  fill_grid_def( grdefh(ngr),nnxp1(ngr),nnyp1(ngr),nnzp1(ngr) &
                       ,nzg1,nzs1,npatch1,platn1(1),plonn1(1)  &
                       ,xtn1(:,ngr),xmn1(:,ngr),ytn1(:,ngr),ymn1(:,ngr)  &
                       ,ztn1(:,ngr),zmn1(:,ngr),topt1(:,ngr) )
enddo
do ngr=1,ngrids
   call alloc_grid_def( grdefn(ngr),nnxp(ngr),nnyp(ngr),nnzp(ngr) )
   call  fill_grid_def( grdefn(ngr),nnxp(ngr),nnyp(ngr),nnzp(ngr) &
                       ,nzg,nzs,npatch,platn(1),plonn(1)  &
                       ,xtn(:,ngr),xmn(:,ngr),ytn(:,ngr),ymn(:,ngr)  &
                       ,ztn(:,ngr),zmn(:,ngr),grid_g(ngr)%topta )
enddo

! See if the history grids match any of the new grids...assuming 1:1 grid number
!   correspondence for now

do ngr=1,min(ngrids,ngrids1)
   call compare_grid_def(grdefh(ngr),grdefn(ngr),'cond_update',ierr)
   if (ierr /= 0) then
      ! No match...
      igrid_match(ngr)=-1
   else
      ! We have a match...
      igrid_match(ngr)=ngr
   endif
enddo


! Finally, process the fields...


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


rewind inhunt

read_loop2: do nvh=1,nvbtab

   ! Read a variable
   nptsh=hr_table(nvh)%nvalues
   read(inhunt)(scr(i),i=1,nptsh)
   
   ngr=hr_table(nvh)%ngrid
   print*,'++++ read: ',trim(hr_table(nvh)%string), ' ',ngr
   
   ng_start=1
   if (ngr > 1 .and. igrid_match(ngr-1) > 0) ng_start = ngr

   grid_loop2:   do ng = ng_start, ngrids
      
      do nv = 1,num_var(ngr)
      
         npts=vtab_r(nv,ng)%npts
      
         !  See if this variable is active in the current run,
         !      but only interpolate if water-related
         if(hr_table(nvh)%string == vtab_r(nv,ng)%name) then

            if (vtab_r(nv,ng)%name == 'RTP') then
               if ( igrid_match(ngr) == ng ) then
                  print 33,'cond_update: filling: ',ngr &
                    ,ng,vtab_r(nv,ng)%name,npts
                  call atob(nptsh,scr,varinit_g(ng)%varrfh)

               else
                  print 33,'cond_update: interpolating: ',ngr &
                       ,ng,vtab_r(nv,ng)%name,npts
                  call hi_interp(nnzp1(ngr),nnxp1(ngr),nnyp1(ngr),1,scr  &
                                ,xmn1(:,ngr),xtn1(:,ngr)  &
                                ,ymn1(:,ngr),ytn1(:,ngr)  &
                                ,zmn1(:,ngr),ztn1(:,ngr)  &
                                ,platn1(ngr),plonn1(ngr)  &
                                ,topt1(:,ngr),ztop1  &
                                ,nnzp(ng),nnxp(ng),nnyp(ng),1  &
                                ,varinit_g(ng)%varrfh  &
                                ,ng,ngr,vtab_r(nv,ng)%name,3)
               endif
               cycle grid_loop2

            elseif (vtab_r(nv,ngr)%name == 'RCP' .or.  &
                    vtab_r(nv,ngr)%name == 'RRP' .or.  &
                    vtab_r(nv,ngr)%name == 'RPP' .or.  &
                    vtab_r(nv,ngr)%name == 'RAP' .or.  &
                    vtab_r(nv,ngr)%name == 'RSP' .or.  &
                    vtab_r(nv,ngr)%name == 'RGP' .or.  &
                    vtab_r(nv,ngr)%name == 'RHP'       ) then
               if ( igrid_match(ngr) == ng ) then
                  print 33,'cond_update: filling: ',ngr &
                    ,ng,vtab_r(nv,ng)%name,npts
                  call nud_cond_accum(nnzp(ng),nnxp(ng),nnyp(ng)  &
                                     ,varinit_g(ng)%varcfh(:,:,:),scr )
               else
                  print 33,'cond_update: interpolating: ',ngr &
                       ,ng,vtab_r(nv,ng)%name,npts
                  call hi_interp(nnzp1(ngr),nnxp1(ngr),nnyp1(ngr),1,scr  &
                                ,xmn1(:,ngr),xtn1(:,ngr)  &
                                ,ymn1(:,ngr),ytn1(:,ngr)  &
                                ,zmn1(:,ngr),ztn1(:,ngr)  &
                                ,platn1(ngr),plonn1(ngr)  &
                                ,topt1(:,ngr),ztop1  &
                                ,nnzp(ng),nnxp(ng),nnyp(ng),1  &
                                ,scr2  &
                                ,ng,ngr,vtab_r(nv,ng)%name,3)                         
                  call nud_cond_accum(nnzp(ng),nnxp(ng),nnyp(ng)  &
                                     ,varinit_g(ng)%varcfh(:,:,:),scr2 )
               endif
               cycle grid_loop2
               
            endif
         endif
      enddo
      
   enddo grid_loop2
print*,'+++++end read loop'
enddo read_loop2

33 format(a30,2i5,3x,a8,i8)

! Close the input history file

close(inhunt)
close(iunhd)

deallocate(scr,scr2,hr_table)
deallocate(xmn1,xtn1,ymn1,ytn1,zmn1,ztn1)

return
end

!-----------------------------------------------

subroutine nud_cond_accum(n1,n2,n3,a,b)

implicit none

integer :: n1,n2,n3
real, dimension(n1,n2,n3) :: a,b

a=a+b

return
end
