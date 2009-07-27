!===================================== Change Log =========================================!
! 5.0.0                                                                                    !
!                                                                                          !
! MLO - 09/29/08 Including the new Grell related variables.                                !
!                                                                                          !
!==========================================================================================!
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved                       !
!  Regional Atmospheric Modeling System - RAMS                                             !
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine will initialize the fields from history. This is called whenever the  !
! RAMSIN variable RUNTYPE is set to HISTORY. This will NOT be called when INITIAL is set   !
! to 3 (when subroutine inithis is called instead).                                        !
!------------------------------------------------------------------------------------------!
subroutine history_start(name_name)
   use grid_dims
   use var_tables
   use io_params
   use mem_grid
   use ref_sounding
   use mem_aerad, ONLY: nwave
   use mem_cuparm, only : nclouds,nnqparm
   use therm_lib , only : virtt,vapour_on
   use rconstants, only : day_sec,hr_sec,min_sec

   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   character(len=*)           , intent(in) :: name_name
   !----- Local variables -----------------------------------------------------------------!
   character (len=80)                      :: hnameinh,prefix
   character (len=2)                       :: cng
   integer                    , save       :: iunhd=11
   integer                                 :: ngrids1,ioutput1,nzg1,nzs1,npatch1,nclouds1
   integer                                 :: iyr,imn,idy,itm,ie,maxarr,ngr,nc
   integer, dimension(maxgrds)             :: nnxp1,nnyp1,nnzp1
   !----- External functions --------------------------------------------------------------!
   integer                    , external   :: cio_i,cio_f,cio_i_sca,cio_f_sca,cio_f8_sca
   !--------------------------------------------------------------------------------------
   !----- Open the input history header file and read some of the info. -------------------!
   nc=len_trim(hfilin)
   hnameinh=hfilin(1:nc-9)//'.vfm'

   call rams_f_open(iunhd,hfilin,'FORMATTED','OLD','READ',0)

   ie=cio_i_sca(iunhd,1,'ngrids',ngrids1,1)
   ngridsh=ngrids1
   ie=cio_i(iunhd,1,'nnxp',nnxp1,ngrids1)
   ie=cio_i(iunhd,1,'nnyp',nnyp1,ngrids1)
   ie=cio_i(iunhd,1,'nnzp',nnzp1,ngrids1)
   ie=cio_i_sca(iunhd,1,'npatch',npatch1,1)
   ie=cio_i_sca(iunhd,1,'nclouds',nclouds1,1)
   ie=cio_i_sca(iunhd,1,'nzg',nzg1,1)
   ie=cio_i_sca(iunhd,1,'nzs',nzs1,1)
   ie=cio_i_sca(iunhd,1,'ioutput',ioutput1,1)
   ie=cio_f8_sca(iunhd,1,'time',time,1)

   !----- Get the 1-d reference state -----------------------------------------------------!
   do ngr=1,ngridsh
      write(cng,fmt='(i2.2)') ngr
      ie=cio_f(iunhd,1,  'u01dn'//cng,  u01dn(:,ngr),nnzp(ngr))
      ie=cio_f(iunhd,1,  'v01dn'//cng,  v01dn(:,ngr),nnzp(ngr))
      ie=cio_f(iunhd,1, 'pi01dn'//cng, pi01dn(:,ngr),nnzp(ngr))
      ie=cio_f(iunhd,1, 'th01dn'//cng, th01dn(:,ngr),nnzp(ngr))
      ie=cio_f(iunhd,1, 'dn01dn'//cng, dn01dn(:,ngr),nnzp(ngr))
      ie=cio_f(iunhd,1, 'rt01dn'//cng, rt01dn(:,ngr),nnzp(ngr))
      ie=cio_f(iunhd,1,'co201dn'//cng,co201dn(:,ngr),nnzp(ngr))
   end do

   !----- Put these into regular arrays (for moving grids) --------------------------------!
   ie=cio_i(iunhd,1,'ninest',ninest,ngrids1)
   ie=cio_i(iunhd,1,'njnest',njnest,ngrids1)


   !----- Check time on file for time requested -------------------------------------------!
   if(dabs(time-timstr) > .1*dtlong)then
      write(unit=*,fmt=*) ' !!! History start error                     !!!'
      write(unit=*,fmt=*) ' !!! Requested time does not match file time !!!'
      write(unit=*,fmt=*) ' !!! Requested time, file time               !!!'
      write(unit=*,fmt=*) ' !!! TIMSTR,time,dtlong=',timstr,time,dtlong
      write(unit=*,fmt=*) ' !!! TIMSTR(m,h,d)=',time/min_sec,time/hr_sec,time/day_sec
      call abort_run ('HISTORY_START time error','history_start','rio.f90')
   endif

   !---------------------------------------------------------------------------------------!
   !     Find maximum size of any array on history file. Allocate scratch array of this    !
   ! size.                                                                                 !
   !---------------------------------------------------------------------------------------!
   maxarr=0
   do ngr=1,ngridsh
      maxarr=max(maxarr,nnxp1(ngr)*nnyp1(ngr)*nnzp1(ngr)                                   &
                ,nnxp1(ngr)*nnyp1(ngr)*nnzp1(ngr)*nclouds1                                 &
                ,nnxp1(ngr)*nnyp1(ngr)*nzg1*npatch1                                        &
                ,nnxp1(ngr)*nnyp1(ngr)*nzs1*npatch1                                        &
                ,nnxp1(ngr)*nnyp1(ngr)*nwave)
   end do

   !----- Read stuff here -----------------------------------------------------------------!
   call hist_read(maxarr,hnameinh,iunhd)

   write (unit=*,fmt=*) 'Back from read! '
   close (unit=iunhd)

   return
end subroutine history_start
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine reads the history file and load the data into the actual arrays.      !
!------------------------------------------------------------------------------------------!
subroutine hist_read(maxarr,hnamein,iunhd)

   use an_header
   use var_tables

   implicit none

   !----- External variable declaration ---------------------------------------------------!
   include 'interface.h'
   !----- Arguments -----------------------------------------------------------------------!
   integer                                     , intent(in) :: maxarr,iunhd
   character (len=*)                           , intent(in) :: hnamein
   !----- Local variables -----------------------------------------------------------------!
   integer                                                  :: ngr,npts,nptsh,nv,nvh,i
   character(len=10)                                        :: post
   real             , allocatable, dimension(:)             :: scr
   type (head_table), allocatable              , save       :: hr_table(:)
   !----- Local constants -----------------------------------------------------------------!
   integer                                     , parameter  :: inhunt=10
   !---------------------------------------------------------------------------------------!


   !----- Allocating scratch array --------------------------------------------------------!
   allocate (scr(maxarr))

   !----- Read variable header information ------------------------------------------------!
   rewind(unit= iunhd)
   read  (unit=iunhd,fmt=*) nvbtab
   allocate (hr_table(nvbtab))
   do nv=1,nvbtab
      read(unit=iunhd,fmt=*)  hr_table(nv)%string   ,hr_table(nv)%npointer                 &
                             ,hr_table(nv)%idim_type,hr_table(nv)%ngrid                    &
                             ,hr_table(nv)%nvalues
   end do


   !----- Open history data file ----------------------------------------------------------!
   call rams_f_open(inhunt,hnamein,'UNFORMATTED','OLD','READ',0)

   !----- Loop through all variables ------------------------------------------------------!
   do nvh=1,nvbtab
      !----- Read a variable --------------------------------------------------------------!
      nptsh = hr_table(nvh)%nvalues
      read(unit=inhunt) (scr(i),i=1,nptsh)

      !----- See if this variable is active in the current run ----------------------------!
      ngr = hr_table(nvh)%ngrid
      if(ngr > nvgrids) cycle

      do nv = 1,num_var(ngr)
         npts=vtab_r(nv,ngr)%npts
         if(hr_table(nvh)%string == vtab_r(nv,ngr)%name) then
            if(nptsh /= npts) then
               print*,'Grid point number mismatch on history field:',  &
                    vtab_r(nv,ngr)%name,npts,nptsh
               call abort_run ('History read number points error','hist_read','rio.f90')
            end if

            write (unit=*,fmt='(a25,2i5,3x,a18,i10)')                                      &
                  'History filling grid: ',ngr,nv,vtab_r(nv,ngr)%name,npts
            call atob(npts,scr(1),vtab_r(nv,ngr)%var_p)
            exit
         end if
      end do

   end do

   !----- Close the input history file then freeing some memory ---------------------------!
   close(unit=inhunt,status='keep')
   deallocate(scr,hr_table)

   return
end subroutine hist_read
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This routine writes the chosen variables on the history file.                        !
!------------------------------------------------------------------------------------------!
subroutine hiswrt(restart)

   use an_header
   use var_tables
   use mem_scratch
   use mem_grid
   use io_params

   implicit none

   !----- External variable declaration ---------------------------------------------------!
   include 'interface.h'
   !----- Arguments -----------------------------------------------------------------------!
   character (len=*)                           , intent(in) :: restart
   !----- Local variables -----------------------------------------------------------------!
   character(len=10)                                        :: c0, c1
   character(len=128)                          , save       :: hnameold,hnameoldh
   character(len=256)                                       :: hnamel,hnamelh
   logical                                                  :: hereitis
   integer                                                  :: nv,nwordh,ngr,nvcnt
   integer                                     , save       :: iohunt=10, ncall=0
   integer                                     , save       :: ncall_head=0,nvtoth=0
   type (head_table), dimension(:), allocatable, save       :: hw_table
   !----- Local constants -----------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
  
   if (ioutput == 0) return

   if (ncall_head == 0) then
      !-----  Find total number of fields to be written -----------------------------------!
      do ngr=1,ngrids
         do nv = 1,num_var(ngr)
            if (vtab_r(nv,ngr)%ihist == 1) nvtoth=nvtoth+1
         end do
      end do
      allocate (hw_table(nvtoth))
      ncall_head=1
   end if

   !----- Open a new output file. ---------------------------------------------------------!
   if(restart == 'yes') then
      call makefnam(hnamel ,hfilout,time,iyeara,imontha,idatea,itimea*100,'R'   ,'$','vfm')
      call makefnam(hnamelh,hfilout,time,iyeara,imontha,idatea,itimea*100,'R','head','txt')
   else
      call makefnam(hnamel ,hfilout,time,iyeara,imontha,idatea,itimea*100,'H'   ,'$','vfm')
      call makefnam(hnamelh,hfilout,time,iyeara,imontha,idatea,itimea*100,'H','head','txt')
   endif
   call rams_f_open(iohunt,hnamel,'UNFORMATTED','REPLACE','WRITE',iclobber)

   !----- Loop through each nest ----------------------------------------------------------!
   nwordh=0
   nvcnt=0
   do ngr=1,ngrids

      !----- Loop through the main variable table and write hist flagged variables. -------!
      do nv = 1,num_var(ngr)
         if ( vtab_r(nv,ngr)%ihist == 1) then
            call writebin(iohunt,vtab_r(nv,ngr)%var_p,vtab_r(nv,ngr)%npts)
            nwordh                    = nwordh+vtab_r(nv,ngr)%npts
            nvcnt                     = nvcnt+1
            hw_table(nvcnt)%string    = vtab_r(nv,ngr)%name
            hw_table(nvcnt)%npointer  = 0
            hw_table(nvcnt)%idim_type = vtab_r(nv,ngr)%idim_type
            hw_table(nvcnt)%ngrid     = ngr
            hw_table(nvcnt)%nvalues   = vtab_r(nv,ngr)%npts
         end if
      end do
   end do

   write(c0,"(f10.1)") time
   write(c1,"(i10)") nwordh
   write(*,"(/,a)") " === History write at Sim time "//trim(adjustl(c0))//" ==="
   write(*,"(a,/)") " === wrote "//trim(adjustl(c1))//" words to file "//&
                    &trim(adjustl(hnamel))//" ==="

   !----- Close history file --------------------------------------------------------------!
   close(iohunt)

   !----- Write the COMMON info out to the header file. -----------------------------------!
   call rams_f_open(iohunt,hnamelh,'FORMATTED','REPLACE','WRITE',iclobber)
   write(unit=iohunt,fmt='(i6)') nvcnt
   do nv=1,nvcnt
      write(unit=iohunt,fmt='(a16,1x,i12,i3,i3,1x,i9)')                                    &
            hw_table(nv)%string   ,hw_table(nv)%npointer                                   &
           ,hw_table(nv)%idim_type,hw_table(nv)%ngrid                                      &
           ,hw_table(nv)%nvalues
   end do
   call commio('HIST','WRITE',iohunt)
   close(unit=iohunt)

   !----- DO NOT remove the old history file if doing a restart or if IFLAG is set. -------!
   if (ihistdel == 1) then
      if (ncall == 0) then
         hnameold  = hnamel
         hnameoldh = hnamelh
      end if
      if (ncall == 1 .and. iflag == 0) then
         inquire (file=trim(hnameold),exist=hereitis)
         !----- This is to avoid system calls that cause infiniband to crash --------------!
         if (hereitis) then
            open  (unit=76,file=trim(hnameold))
            close (unit=76,status='delete')
         end if
         inquire (file=trim(hnameoldh),exist=hereitis)
         if (hereitis) then
            open  (unit=76,file=trim(hnameoldh))
            close (unit=76,status='delete')
         end if
         hnameold  = hnamel
         hnameoldh = hnamelh
      end if
      ncall = 1
   endif

   return
end subroutine hiswrt
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!   This subroutine simply adds two variables. The only reason we do it here instead of    !
! directly adding is because of the pointers, which may have different ranks.              !
!------------------------------------------------------------------------------------------!
subroutine rams_aprep_p (n1,a,b,c)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   integer              , intent(in) :: n1
   real   , dimension(*), intent(in)  :: a,b
   real   , dimension(*), intent(out) :: c
   !----- Local variables -----------------------------------------------------------------!
   integer :: i
   !---------------------------------------------------------------------------------------!

   do i=1,n1
      c(i)=a(i)+b(i)
   end do

   return
end subroutine rams_aprep_p
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine converts the turbulence coefficients (HKM to HKH). It depends on the  !
! turbulence closure used.                                                                 !
!------------------------------------------------------------------------------------------!
subroutine rams_aprep_hkh(n1,hkm,vkh,dn0,scr1,idiffk,xkhkm)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   integer              , intent(in)  :: n1,idiffk
   real                 , intent(in)  :: xkhkm
   real   , dimension(*), intent(in)  :: hkm,vkh,dn0
   real   , dimension(*), intent(out) :: scr1
   !----- Local variables -----------------------------------------------------------------!
   integer :: ind
   !---------------------------------------------------------------------------------------!
   
   select case (idiffk)
   case (1,2,3,7)
      do ind = 1,n1
         scr1(ind) = hkm(ind) * xkhkm / dn0(ind)
      end do
   case (4,5,6)
      do ind = 1,n1
         scr1(ind) = vkh(ind) / dn0(ind)
      end do
   end select

   return
end subroutine rams_aprep_hkh
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine converts the turbulence coefficients by removing the density use for  !
! weighting.                                                                               !
!------------------------------------------------------------------------------------------!
subroutine rams_aprep_vkh(n1,vkh,dn0,vt3dd)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   integer              , intent(in)  :: n1
   real   , dimension(*), intent(in)  :: vkh,dn0
   real   , dimension(*), intent(out) :: vt3dd
   !----- Local variables -----------------------------------------------------------------!
   integer :: ind
   !---------------------------------------------------------------------------------------!

   do ind = 1,n1
      vt3dd(ind) = vkh(ind) / dn0(ind)
   end do

   return
end subroutine rams_aprep_vkh
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine writes the chosen variables on the full analysis files.                                       !
!------------------------------------------------------------------------------------------!
subroutine anlwrt(restart,vtype)

   use an_header
   use var_tables
   use mem_scratch
   use mem_basic
   use mem_turb
   use mem_grid
   use io_params
   use mem_cuparm, only: nclouds
   use mem_aerad , only: nwave

   implicit none

   !----- External variable declaration ---------------------------------------------------!
   include 'interface.h'
   !----- Arguments -----------------------------------------------------------------------!
   character(len=*)                            , intent(in) :: restart,vtype
   !----- Local variables -----------------------------------------------------------------!
   character(len=1)                                         :: vnam
   character(len=2)                                         :: cgrid
   character(len=10)                                        :: c0
   character(len=16)                                        :: varn
   character(len=25)                                        :: subaname
   character(len=128)                                       :: anamel,anamelh
   integer                                                  :: ngr,nv,nvcnt,lenl,npointer
   integer                                                  :: n3d,indwrt,n2d
   integer                                     , save       :: ioaunt=10,ncall_head=0
   integer                                     , save       :: nvtota=0,nvtotl=0,nvtot
   logical                                                  :: exans
   logical                                                  :: found
   real(kind=8)                                             :: timeold
   type (head_table), dimension(:), allocatable, save       :: aw_table(:)
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !    To avoid rewriting analysis when it is a history restart. It will write if it      ! 
   ! doesn't find all the analysis (header and grid) there, though.                        !
   !---------------------------------------------------------------------------------------!
   found = .false.


   if (ioutput == 0) return

   if (ncall_head == 0) then
      !----- Find total number of fields to be written. -----------------------------------!
      do ngr=1,ngrids
         do nv = 1,num_var(ngr)
            if ( vtab_r(nv,ngr)%ianal == 1) nvtota=nvtota+1
            if ( vtab_r(nv,ngr)%ilite == 1) nvtotl=nvtotl+1
         end do
      end do
      nvtot=max(nvtota,nvtotl)
      allocate (aw_table(nvtot))
      ncall_head=1
   end if


   timeold=time
   if(vtype == 'MEAN'.or.vtype == 'BOTH') time=min(time,time-avgtim/2.)


   !----- Construct header file name. -----------------------------------------------------!
   select case (vtype)
   case ('INST')
      vnam='A'
   case ('LITE')
      vnam='L'
   case ('MEAN')
      vnam='M'
   case ('BOTH')
      vnam='B'
   end select
   call makefnam(anamelh,afilout,time,iyeara,imontha,idatea,itimea*100,vnam,'head','txt')


   !---------------------------------------------------------------------------------------!
   !     Here is a point that is called just at history start. It checks whether the       !
   ! analysis are already there. If they are, it will return without writing anything      !
   !---------------------------------------------------------------------------------------!
   if (restart == 'yes') then
      inquire (file=trim(anamelh),exist=found)
      if (found) then
         gridloop: do ngr=1,ngrids
            write(cgrid,fmt='(a1,i1)') 'g',ngr
            call makefnam(anamel,afilout,time,iyeara,imontha,idatea,itimea*100,vnam,cgrid  &
                         ,'vfm')
            inquire(file=trim(anamel),exist=found)
            if (.not.found) exit gridloop
         end do gridloop
      end if
      !------------------------------------------------------------------------------------!
      !    So if any of the files were missing, I continue, otherwise, I leave the sub-    !
      ! routine.                                                                           !
      !------------------------------------------------------------------------------------!
      if (found) return
   end if

   !----- Loop through each nest ----------------------------------------------------------!
   nvcnt=0

   do ngr=1,ngrids
      write(cgrid,'(a1,i1)') 'g',ngr
      call makefnam(anamel,afilout,time,iyeara,imontha,idatea,itimea*100,vnam,cgrid,'vfm')

      lenl = len_trim(anamel)

      !------------------------------------------------------------------------------------!
      !     Check whether the file already exists. In case there is such file already,     !
      ! decide between overwriting it or stopping the run.                                 !
      !------------------------------------------------------------------------------------!
      inquire(file=anamel,exist=exans)
      if(exans.and.iclobber.eq.0) then
         print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         print*,'!!!   trying to open file name :'
         print*,'!!!       ',anamel
         print*,'!!!   but it already exists. run is ended.'
         print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         call abort_run('Analysis file already existed and you said I can''t overwrite!'   &
                       ,'anlwrt','rio.f90')
      end if

      call rams_c_open(anamel(1:lenl)//char(0),'w'//char(0))
      npointer=0

      !------------------------------------------------------------------------------------!
      !     Loop through the main variable table and write those variables with the        !
      ! correct flag set.                                                                  !
      !------------------------------------------------------------------------------------!

      do nv = 1,num_var(ngr)

         !---------------------------------------------------------------------------------!
         !    Writing instantaneous analysis.                                              !
         !---------------------------------------------------------------------------------!
         if ((vtype == 'INST' .and. vtab_r(nv,ngr)%ianal == 1) .or.                        &
            (vtype == 'LITE' .and. vtab_r(nv,ngr)%ilite == 1)) then

            varn= vtab_r(nv,ngr)%name

            !------------------------------------------------------------------------------!
            !     First check whether the variable is a special case:                      !
            !------------------------------------------------------------------------------!
            !----- Exner perturbation, save full Exner function instead -------------------!
            if (varn == 'PP') then
              !----- Output total Exner function ------------------------------------------!
              call RAMS_aprep_p (nnxyzp(ngr),vtab_r(nv,ngr)%var_p,basic_g(ngr)%pi0         &
                                ,scratch%scr1)
              varn='PI'
              !----- Rearrange 3-d variables to (x,y,z) -----------------------------------!
              call rearrange(nnzp(ngr),nnxp(ngr),nnyp(ngr),scratch%scr1,scratch%scr2)

            !----- Removing density from the coefficient ----------------------------------!
            elseif(varn == 'HKM') then
               !----- Convert to HKM to HKH (note that VKH is HKH for Deardorff) ----------!
               call RAMS_aprep_hkh (nnxyzp(ngr),vtab_r(nv,ngr)%var_p,turb_g(ngr)%vkh       &
                                   ,basic_g(ngr)%dn0,scratch%scr1,idiffk(ngr),xkhkm(ngr))
               varn='HKH'
               !-----  Rearrange 3-d variables to (x,y,z) ---------------------------------!
               call rearrange(nnzp(ngr),nnxp(ngr),nnyp(ngr),scratch%scr1,scratch%scr2)

            !----- Removing density from the coefficient ----------------------------------!
            elseif(varn == 'VKH') then
               !----- Un-density weight VKH -----------------------------------------------!
               call RAMS_aprep_vkh (nnxyzp(ngr),vtab_r(nv,ngr)%var_p,basic_g(ngr)%dn0      &
                                   ,scratch%scr1)
               !----- Rearrange 3-d variables to (x,y,z) ----------------------------------!
               call rearrange(nnzp(ngr),nnxp(ngr),nnyp(ngr),scratch%scr1,scratch%scr2)

            !------------------------------------------------------------------------------!
            !    Now, the ordinary variables                                               !
            !------------------------------------------------------------------------------!
            !-----  Rearrange 3-d variables to (x,y,z) ------------------------------------!
            elseif(vtab_r(nv,ngr)%idim_type == 3) then
               call rearrange(nnzp(ngr),nnxp(ngr),nnyp(ngr),vtab_r(nv,ngr)%var_p           &
                             ,scratch%scr2)

            !----- Rearrange 4-d leaf%soil variables to (x,y,z,patch) ---------------------!
            elseif(vtab_r(nv,ngr)%idim_type == 4) then
               call rearrange_p(nnxp(ngr),nnyp(ngr),nzg,npatch,vtab_r(nv,ngr)%var_p        &
                               ,scratch%scr2)

            !----- Rearrange 4-d leaf%sfcwater variables to (x,y,z,patch) -----------------!
            elseif(vtab_r(nv,ngr)%idim_type == 5) then
               call rearrange_p(nnxp(ngr),nnyp(ngr),nzs,npatch,vtab_r(nv,ngr)%var_p        &
                               ,scratch%scr2)

            !----- Rearrange 4-d cuparm variables to (x,y,z,cloud) ------------------------!
            elseif (vtab_r(nv,ngr)%idim_type == 8) then
               call rearrange_p(nnxp(ngr),nnyp(ngr),nnzp(ngr),nclouds,vtab_r(nv,ngr)%var_p &
                               ,scratch%scr2)

            !------------------------------------------------------------------------------!
            !     For types 2, 6, 7, and 9 we don't need to change the order, but I need   !
            ! to copy to scr2, so I will use a dum function to do that.                    !
            !------------------------------------------------------------------------------!
            !----- Copy 2-d (x,y,1) -------------------------------------------------------!
            elseif (vtab_r(nv,ngr)%idim_type == 2) then
               call rearrange_dum(nnxp(ngr),nnyp(ngr),1,vtab_r(nv,ngr)%var_p               &
                                 ,scratch%scr2)

            !----- Copy 3-d (x,y,patch) ---------------------------------------------------!
            elseif (vtab_r(nv,ngr)%idim_type == 6) then
               call rearrange_dum(nnxp(ngr),nnyp(ngr),npatch,vtab_r(nv,ngr)%var_p          &
                                 ,scratch%scr2)

            !----- Copy 3-d (x,y,wave) ----------------------------------------------------!
            elseif (vtab_r(nv,ngr)%idim_type == 7) then
               call rearrange_dum(nnxp(ngr),nnyp(ngr),nwave,vtab_r(nv,ngr)%var_p           &
                                 ,scratch%scr2)

            !----- Copy 3-d (x,y,cloud) ---------------------------------------------------!
            elseif (vtab_r(nv,ngr)%idim_type == 9) then
               call rearrange_dum(nnxp(ngr),nnyp(ngr),nclouds,vtab_r(nv,ngr)%var_p         &
                                 ,scratch%scr2)
            end if
            
            !------------------------------------------------------------------------------!
            !    Adding the variable to the list and writing the variable to the data      !
            ! file.                                                                        !
            !------------------------------------------------------------------------------!
            nvcnt = nvcnt+1
            aw_table(nvcnt)%string    = varn
            aw_table(nvcnt)%npointer  = npointer
            aw_table(nvcnt)%idim_type = vtab_r(nv,ngr)%idim_type
            aw_table(nvcnt)%ngrid     = ngr
            aw_table(nvcnt)%nvalues   = vtab_r(nv,ngr)%npts
            !----- Writing ----------------------------------------------------------------!
            call vforecr(ioaunt,scratch%scr2 ,vtab_r(nv,ngr)%npts,18,scratch%scr1          &
                        ,scratch%scr1,'LIN',npointer)
            !------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !    Writing instantaneous analysis.                                              !
         !---------------------------------------------------------------------------------!
         elseif((vtype == 'MEAN' .and. vtab_r(nv,ngr)%ianal == 1) .or. &
                (vtype == 'BOTH' .and. vtab_r(nv,ngr)%ilite == 1)) then

            varn= vtab_r(nv,ngr)%name


            !------------------------------------------------------------------------------!
            !     First check whether the variable is a special case:                      !
            !------------------------------------------------------------------------------!
            !----- Exner perturbation, save full Exner function instead -------------------!
            if (varn == 'PP') then
              !----- Output total Exner function ------------------------------------------!
              call RAMS_aprep_p (nnxyzp(ngr),vtab_r(nv,ngr)%var_m,basic_g(ngr)%pi0         &
                                ,scratch%scr1)
              varn='PI'
              !----- Rearrange 3-d variables to (x,y,z) -----------------------------------!
              call rearrange(nnzp(ngr),nnxp(ngr),nnyp(ngr),scratch%scr1,scratch%scr2)

            !----- Removing density from the coefficient ----------------------------------!
            elseif(varn == 'HKM') then
               !----- Convert to HKM to HKH (note that VKH is HKH for Deardorff) ----------!
               call RAMS_aprep_hkh (nnxyzp(ngr),vtab_r(nv,ngr)%var_m                       &
                                   ,turb_g(ngr)%vkh,basic_g(ngr)%dn0,scratch%scr1          &
                                   ,idiffk(ngr),xkhkm(ngr))
               varn='HKH'
               !-----  Rearrange 3-d variables to (x,y,z) ---------------------------------!
               call rearrange(nnzp(ngr),nnxp(ngr),nnyp(ngr),scratch%scr1,scratch%scr2)

            !----- Removing density from the coefficient ----------------------------------!
            elseif(varn == 'VKH') then
               !----- Un-density weight VKH -----------------------------------------------!
               call RAMS_aprep_vkh (nnxyzp(ngr),vtab_r(nv,ngr)%var_m                       &
                                   ,basic_g(ngr)%dn0,scratch%scr1)
               !----- Rearrange 3-d variables to (x,y,z) ----------------------------------!
               call rearrange(nnzp(ngr),nnxp(ngr),nnyp(ngr),scratch%scr1,scratch%scr2)

            !------------------------------------------------------------------------------!
            !    Now, the ordinary variables                                               !
            !------------------------------------------------------------------------------!
            !-----  Rearrange 3-d variables to (x,y,z) ------------------------------------!
            elseif(vtab_r(nv,ngr)%idim_type == 3) then
               call rearrange(nnzp(ngr),nnxp(ngr),nnyp(ngr),vtab_r(nv,ngr)%var_m           &
                             ,scratch%scr2)

            !----- Rearrange 4-d leaf%soil variables to (x,y,z,patch) ---------------------!
            elseif(vtab_r(nv,ngr)%idim_type == 4) then
               call rearrange_p(nnxp(ngr),nnyp(ngr),nzg,npatch,vtab_r(nv,ngr)%var_m        &
                               ,scratch%scr2)

            !----- Rearrange 4-d leaf%sfcwater variables to (x,y,z,patch) -----------------!
            elseif(vtab_r(nv,ngr)%idim_type == 5) then
               call rearrange_p(nnxp(ngr),nnyp(ngr),nzs,npatch,vtab_r(nv,ngr)%var_m        &
                               ,scratch%scr2)

            !----- Rearrange 4-d cuparm variables to (x,y,z,cloud) ------------------------!
            elseif (vtab_r(nv,ngr)%idim_type == 8) then
               call rearrange_p(nnxp(ngr),nnyp(ngr),nnzp(ngr),nclouds,vtab_r(nv,ngr)%var_m &
                               ,scratch%scr2)

            !------------------------------------------------------------------------------!
            !     For types 2, 6, 7, and 9 we don't need to change the order, but I need   !
            ! to copy to scr2, so I will use a dum function to do that.                    !
            !------------------------------------------------------------------------------!
            !----- Copy 2-d (x,y,1) -------------------------------------------------------!
            elseif (vtab_r(nv,ngr)%idim_type == 2) then
               call rearrange_dum(nnxp(ngr),nnyp(ngr),1,vtab_r(nv,ngr)%var_m               &
                                 ,scratch%scr2)

            !----- Copy 3-d (x,y,patch) ---------------------------------------------------!
            elseif (vtab_r(nv,ngr)%idim_type == 6) then
               call rearrange_dum(nnxp(ngr),nnyp(ngr),npatch,vtab_r(nv,ngr)%var_m          &
                                 ,scratch%scr2)

            !----- Copy 3-d (x,y,wave) ----------------------------------------------------!
            elseif (vtab_r(nv,ngr)%idim_type == 7) then
               call rearrange_dum(nnxp(ngr),nnyp(ngr),nwave,vtab_r(nv,ngr)%var_m           &
                                 ,scratch%scr2)

            !----- Copy 3-d (x,y,cloud) ---------------------------------------------------!
            elseif (vtab_r(nv,ngr)%idim_type == 9) then
               call rearrange_dum(nnxp(ngr),nnyp(ngr),nclouds,vtab_r(nv,ngr)%var_m         &
                                 ,scratch%scr2)
            end if
            
            !------------------------------------------------------------------------------!
            !    Adding the variable to the list and writing the variable to the data      !
            ! file.                                                                        !
            !------------------------------------------------------------------------------!
            nvcnt=nvcnt+1
            aw_table(nvcnt)%string    = varn
            aw_table(nvcnt)%npointer  = npointer
            aw_table(nvcnt)%idim_type = vtab_r(nv,ngr)%idim_type
            aw_table(nvcnt)%ngrid     = ngr
            aw_table(nvcnt)%nvalues   = vtab_r(nv,ngr)%npts
            !----- Writing ----------------------------------------------------------------!
            call vforecr(ioaunt,scratch%scr2(1),vtab_r(nv,ngr)%npts,18,scratch%scr1        &
                        ,scratch%scr1,'LIN',npointer)
            !------------------------------------------------------------------------------!
         end if
      end do

      !----- Close the file and header ----------------------------------------------------!
      call rams_c_close()
      close(unit=ioaunt,status='keep')
   end do

   !----- Write the header information out to the file. -----------------------------------!
   call rams_f_open(ioaunt,anamelh,'FORMATTED','REPLACE','WRITE',iclobber)

   !----- Writing the header --------------------------------------------------------------!
   write(unit=ioaunt,fmt='(i6)') nvcnt
   do nv=1,nvcnt
      write(unit=ioaunt,fmt='(a16,1x,i12,i3,i3,1x,i9)')                                    &
            aw_table(nv)%string ,aw_table(nv)%npointer,aw_table(nv)%idim_type              &
           ,aw_table(nv)%ngrid,aw_table(nv)%nvalues
   end do
   call commio('ANAL','WRITE',ioaunt)
   close(unit=ioaunt)



   !----- Printing the banner -------------------------------------------------------------!
   select case (trim(vtype))
   case('LITE')
      subaname='  Analysis lite write'
   case('MEAN')
      subaname='  Averaged analysis write    '
   case('BOTH')
      subaname='  Averaged analysis lite write    '
   case default
      subaname='  Analysis write         '
   end select
   write(c0,"(f10.1)") time
   write(*,"(/,a)") " === "//trim(adjustl(subaname))//" at Sim time "//trim(adjustl(c0))//" ==="
   write(*,"(a,/)") " === wrote file "//&
        &trim(adjustl(anamel))//" ==="

   !----- Reset the time back to the original ---------------------------------------------!
   if(vtype == 'MEAN'.or.vtype == 'BOTH')  time=timeold

   return
end subroutine anlwrt
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine does nothing actually. It simply copies the data from the original    !
! place to the buffer for 3D variables.                                                    !
!------------------------------------------------------------------------------------------!
subroutine rearrange_dum(n2,n3,n4,mydata,buff)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   integer, intent(in)                       :: n2,n3,n4
   real,    intent(in) , dimension(n2,n3,n4) :: mydata
   real,    intent(out), dimension(n2,n3,n4) :: buff
   !----- Local variables -----------------------------------------------------------------!
   integer                                   :: i,j,k
   !---------------------------------------------------------------------------------------!
   do i=1,n2
      do j=1,n3
         do k=1,n4
           buff(i,j,k) = mydata(i,j,k)
         end do
      end do
   end do
   return
end subroutine rearrange_dum
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine does nothing actually. It simply copies the data from the original    !
! place to the buffer for 4D variables.                                                    !
!------------------------------------------------------------------------------------------!
subroutine rearrange_dum4(n2,n3,n4,n5,mydata,buff)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   integer, intent(in)                          :: n2,n3,n4,n5
   real,    intent(in) , dimension(n2,n3,n4,n5) :: mydata
   real,    intent(out), dimension(n2,n3,n4,n5) :: buff
   !----- Local variables -----------------------------------------------------------------!
   integer                                      :: i,j,k,l
   !---------------------------------------------------------------------------------------!
   do i=1,n2
      do j=1,n3
         do k=1,n4
            do l=1,n5
               buff(i,j,k,l) = mydata(i,j,k,l)
            end do
         end do
      end do
   end do
   return
end subroutine rearrange_dum4
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine moves the vertical dimension from the 1st to the 3rd rank.            !
!------------------------------------------------------------------------------------------!
subroutine rearrange_p(n2,n3,n4,n5,mydata,buff)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   integer                        , intent(in)  :: n2,n3,n4,n5
   real   , dimension(n4,n2,n3,n5), intent(in)  :: mydata
   real   , dimension(n2,n3,n4,n5), intent(out) :: buff
   !----- Local variables -----------------------------------------------------------------!
   integer                                      :: i,j,k,l
   !---------------------------------------------------------------------------------------!

   do l = 1,n5
      do k = 1,n4
         do j = 1,n3
            do i = 1,n2
               buff(i,j,k,l) = mydata(k,i,j,l)
            end do
         end do
      end do
   end do

   return
end subroutine rearrange_p
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine writes a variable on a binary file.                                   !
!------------------------------------------------------------------------------------------!
subroutine writebin(iun,var,npts)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   real   , dimension(*), intent(in) :: var
   !----- Local variables -----------------------------------------------------------------!
   integer              , intent(in) :: iun,npts
   integer                           :: i
   !---------------------------------------------------------------------------------------!

   write(iun) (var(i),i=1,npts)

   return
end subroutine writebin
!==========================================================================================!
!==========================================================================================!
