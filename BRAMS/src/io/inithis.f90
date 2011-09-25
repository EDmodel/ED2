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
!    This subroutine will initialize the fields with history files. This is not called in  !
! HISTORY run, it will be used only when the RAMSIN variable INITIAL is set to 3. If you   !
! are looking for the history runs, check the rio.f90 file, subroutine history_start.      !
!------------------------------------------------------------------------------------------!
subroutine inithis()
   use var_tables
   use an_header
   use mem_basic
   use mem_grid
   use rconstants
   use ref_sounding
   use io_params
   use mem_scratch
   use mem_aerad , only : nwave
   use grid_dims , only : maxgrds
   use mem_cuparm, only : nclouds,nnqparm
   use therm_lib , only : virtt,vapour_on

   implicit none

   !----- Local variables -----------------------------------------------------------------!
   character (len=256)                                  :: hnamel, hnamelh
   character (len=2)                                    :: cng
   integer                                              :: ngrids1,ioutput1,nzg1,nzs1
   integer                                              :: npatch1,nclouds1
   integer                                              :: iyr,imn,idy,itm,ie,ngrmin
   integer                                              :: maxarr,maxarr2
   integer                                              :: ngr,maxx1,maxy1,maxz1
   integer                                              :: npts,nptsh,nv,nvh,i,k,nzpg1,nc
   integer                                       , save :: iunhd=11,inhunt=10
   integer          , allocatable, dimension(:)         :: nnxp1,nnyp1,nnzp1
   integer          , allocatable, dimension(:)         :: nnqparm1
   real                                                 :: ztop1
   real             , allocatable, dimension(:)         :: scr,scr2,scr3
   real             , allocatable, dimension(:)         :: platn1,plonn1
   real             , allocatable, dimension(:)         :: u01dn1,v01dn1,rt01dn1,th01dn1
   real             , allocatable, dimension(:)         :: pi01dn1,dn01dn1,co201dn1
   real             , allocatable, dimension(:,:)       :: xmn1,xtn1,ymn1,ytn1,zmn1,ztn1
   real             , allocatable, dimension(:,:)       :: topt1,parea
   real(kind=8) :: time1
   type (head_table), allocatable, dimension(:)  , save :: hr_table
   !----- External functions --------------------------------------------------------------!
   integer, external :: cio_i,cio_f,cio_i_sca,cio_f_sca,cio_f8_sca
   !---------------------------------------------------------------------------------------!



   !----- Open the input history header file and read some of the info. -------------------!
   call makefnam(hnamel ,hfilin,0.d0,iyearh,imonthh,idateh,itimeh*100,'H'   ,'$','vfm')
   call makefnam(hnamelh,hfilin,0.d0,iyearh,imonthh,idateh,itimeh*100,'H','head','txt')

   call rams_f_open(iunhd,hnamelh,'FORMATTED','OLD','READ',0)

   ie=cio_i_sca(iunhd,1,'ngrids',ngrids1,1)

   allocate (nnxp1(ngrids1),nnyp1(ngrids1),nnzp1(ngrids1))
   allocate (platn1(ngrids1),plonn1(ngrids1),nnqparm1(ngrids1))

   ie=cio_i(iunhd,1,'nnxp',nnxp1,ngrids1)
   ie=cio_i(iunhd,1,'nnyp',nnyp1,ngrids1)
   ie=cio_i(iunhd,1,'nnzp',nnzp1,ngrids1)
   ie=cio_i_sca(iunhd,1,'npatch',npatch1,1)
   ie=cio_i_sca(iunhd,1,'nzg',nzg1,1)
   ie=cio_i_sca(iunhd,1,'nzs',nzs1,1)
   ie=cio_i_sca(iunhd,1,'nclouds',nclouds1,1)
   ie=cio_i_sca(iunhd,1,'ioutput',ioutput1,1)
   ie=cio_f8_sca(iunhd,1,'time',time1,1)
   ie=cio_f_sca(iunhd,1,'ztop',ztop1,1)
   ie=cio_f(iunhd,1,'platn',platn1,ngrids1)
   ie=cio_f(iunhd,1,'plonn',plonn1,ngrids1)
   ie=cio_i(iunhd,1,'nnqparm',nnqparm1,ngrids1)

   !----- Checking whether surface-related dimensions remain the same ---------------------!
   if (nzg /= nzg1 .or. nzs /= nzs1 .or. npatch /= npatch1) then
      write(unit=*,fmt='(a)')             'Surface dimension mismatch!'
      write(unit=*,fmt='(2(a,1x,i5,1x))') 'Your npatch:',npatch,'Hist npatch: ',npatch1
      write(unit=*,fmt='(2(a,1x,i5,1x))') 'Your nzg: ',nzg,'Hist nzg: ',nzg1
      write(unit=*,fmt='(2(a,1x,i5,1x))') 'Your nzs: ',nzs,'Hist nzs: ',nzs1
      call abort_run('LEAF parameters must be same for initial-history start'              &
                    ,'inithis','inithis.f90')
   end if

   !----- Checking whether cumulus scheme is consistent with history set up ---------------!
   ngrmin=min(ngrids,ngrids1)
   if (any(nnqparm(1:ngrmin) > 0) .and.                                                    &
       any(nnqparm(1:ngrmin) /= nnqparm1(1:ngrmin))) then
      write(unit=*,fmt='(a)') 'Cumulus scheme status has changed in at least 1 grid!'
      do ngr=1,ngrmin
         write(unit=*,fmt='(3(a,1x,i5,1x))') 'Grid        :',ngr                           &
                                            ,'Your nnqparm:',nnqparm(ngr)                  &
                                            ,'Hist nnqparm:',nnqparm1(ngr)
      end do
      call abort_run('NNQPARM must be the same for initial history start'                  &
                    ,'inithis','inithis.f90')
   elseif (any(nnqparm(1:ngrids1) > 0) .and. nclouds /= nclouds1) then
      write(unit=*,fmt='(a)')             'Cloud dimension mismatch!'
      write(unit=*,fmt='(2(a,1x,i5,1x))') 'Your nclouds:',nclouds,'Hist nclouds: ',nclouds1
      call abort_run('NCLOUDS must be same for initial-history start'                      &
                    ,'inithis','inithis.f90')
   end if

   !---------------------------------------------------------------------------------------!
   !     Find maximum size of any array on history file. Allocate scratch array of this    !
   ! size.                                                                                 !
   !---------------------------------------------------------------------------------------!
   maxarr=0
   maxarr2=0
   maxx1=0
   maxy1=0
   maxz1=0
   do ngr=1,ngrids1
      maxarr=max(maxarr,nnxp1(ngr)*nnyp1(ngr)*nnzp1(ngr)                                   &
                ,nnxp1(ngr)*nnyp1(ngr)*nnzp1(ngr)*nclouds1                                 &
                ,nnxp1(ngr)*nnyp1(ngr)*nzg1*npatch1                                        &
                ,nnxp1(ngr)*nnyp1(ngr)*nzs1*npatch1)
      maxarr2=max(maxarr2,nnxp1(ngr)*nnyp1(ngr))
      maxx1=max(maxx1,nnxp1(ngr))
      maxy1=max(maxy1,nnyp1(ngr))
      maxz1=max(maxz1,nnzp1(ngr))
   enddo

   !---- Allocating variables based on the maximum dimensions -----------------------------!
   allocate(xmn1(maxx1,ngrids1),xtn1(maxx1,ngrids1))
   allocate(ymn1(maxy1,ngrids1),ytn1(maxy1,ngrids1))
   allocate(zmn1(maxz1,ngrids1),ztn1(maxz1,ngrids1))




   do ngr=1,ngrids1
      write(cng,'(i2.2)') ngr
      ie=cio_f(iunhd,1,'xmn'//cng,xmn1(:,ngr),nnxp1(ngr))
      ie=cio_f(iunhd,1,'xtn'//cng,xtn1(:,ngr),nnxp1(ngr))
      ie=cio_f(iunhd,1,'ymn'//cng,ymn1(:,ngr),nnyp1(ngr))
      ie=cio_f(iunhd,1,'ytn'//cng,ytn1(:,ngr),nnyp1(ngr))
      ie=cio_f(iunhd,1,'zmn'//cng,zmn1(:,ngr),nnzp1(ngr))
      ie=cio_f(iunhd,1,'ztn'//cng,ztn1(:,ngr),nnzp1(ngr))
   enddo

   allocate (topt1(maxarr2,ngrids1))
   allocate (parea(maxarr,ngrids1))

   allocate (scr(maxarr))
   allocate (scr2(maxarr))
   allocate (scr3(maxarr))

   call rams_f_open(inhunt,hnamel,'UNFORMATTED','OLD','READ',0)

   !-----  Read variable header info ------------------------------------------------------!
   rewind(iunhd)

   read(unit=iunhd,fmt=*) nvbtab
   allocate (hr_table(nvbtab))
   do nv=1,nvbtab
      read(unit=iunhd,fmt=*)  hr_table(nv)%string,hr_table(nv)%npointer                    &
                             ,hr_table(nv)%idim_type,hr_table(nv)%ngrid                    &
                             ,hr_table(nv)%nvalues
   end do

   !---------------------------------------------------------------------------------------!
   !    Go through file and get all grids' TOPT's. Don't know how else to get this before  !
   ! processing a field...                                                                 !
   !---------------------------------------------------------------------------------------!
   do nvh=1,nvbtab 
      ngr=hr_table(nvh)%ngrid
      nptsh=hr_table(nvh)%nvalues

      read(unit=inhunt)(scr(i),i=1,nptsh)
      select case(trim(hr_table(nvh)%string))
      case ('TOPT')
         call atob(nptsh, scr,topt1(:,ngr))
      case ('PATCH_AREA')
         call atob(nptsh, scr,parea(:,ngr))
      end select
   end do

   rewind(unit=inhunt)


   !----- Need wind rotation for the general case -----------------------------------------!
   !----- Loop through all variables ------------------------------------------------------!
   do nvh=1,nvbtab
      !----- Read a variable --------------------------------------------------------------!
      nptsh=hr_table(nvh)%nvalues
      read(unit=inhunt)(scr(i),i=1,nptsh)

      !------------------------------------------------------------------------------------!
      !    See if this variable is active in the current run and interpolate to new grid   !
      ! structure. For the leaf variables, we will interpolate  everything except TOPxx,   !
      ! SOIL_TEXT, PATCH_AREA, and LEAF_CLASS.                                             !
      !------------------------------------------------------------------------------------!
      ngr=hr_table(nvh)%ngrid

      write (unit=*,fmt=*) 'INITHIS - Start: ',nvh,ngr,hr_table(nvh)%string                &
                                              ,hr_table(nvh)%idim_type,nptsh

      do nv = 1,num_var(1)
         npts=vtab_r(nv,1)%npts
         if(hr_table(nvh)%string == vtab_r(nv,1)%name) then

            write (unit=*,fmt=*) 'INITHIS - Interpolating: ',ngr,' ',nv,1                  &
                                 ,vtab_r(nv,1)%name,npts,vtab_r(nv,1)%idim_type

            !----- 2D variables (nxp,nyp) that can be interpolated ------------------------!
            if (vtab_r(nv,1)%idim_type == 2     .and. hr_table(nvh)%string /= 'TOPT' .and. &
                hr_table(nvh)%string /= 'TOPTA' .and. hr_table(nvh)%string /= 'TOPMA') then
               call hi_interp(1,nnxp1(ngr),nnyp1(ngr),1,scr,xmn1(:,ngr),xtn1(:,ngr)        &
                             ,ymn1(:,ngr),ytn1(:,ngr),zmn1(:,ngr),ztn1(:,ngr)              &
                             ,platn1(ngr),plonn1(ngr),topt1(:,ngr),ztop1,1,nnxp(1),nnyp(1) &
                             ,1,vtab_r(nv,1)%var_p,1,ngr,vtab_r(nv,1)%name,2)
            !----- 3D variables (nzp,nxp,nyp) ---------------------------------------------!
            elseif (vtab_r(nv,1)%idim_type == 3) then
               call hi_interp(nnzp1(ngr),nnxp1(ngr),nnyp1(ngr),1,scr,xmn1(:,ngr)           &
                             ,xtn1(:,ngr),ymn1(:,ngr),ytn1(:,ngr),zmn1(:,ngr),ztn1(:,ngr)  &
                             ,platn1(ngr),plonn1(ngr),topt1(:,ngr),ztop1,nnzp(1),nnxp(1)   &
                             ,nnyp(1),1,vtab_r(nv,1)%var_p,1,ngr,vtab_r(nv,1)%name,3)

            !----- 4D variables (nzg,nxp,nyp,npatch) that can be interpolated -------------!
            elseif (vtab_r(nv,1)%idim_type == 4 .and. hr_table(nvh)%string /= 'SOIL_TEXT') &
            then
               !----- First, interpolate patch 1. -----------------------------------------!
               call hi_interp(nzg1,nnxp1(ngr),nnyp1(ngr),1,scr,xmn1(:,ngr),xtn1(:,ngr)     &
                             ,ymn1(:,ngr),ytn1(:,ngr),zmn1(:,ngr),ztn1(:,ngr),platn1(ngr)  &
                             ,plonn1(ngr),topt1(:,ngr),ztop1,nzg,nnxp(1),nnyp(1),1         &
                             ,vtab_r(nv,1)%var_p,1,ngr,vtab_r(nv,1)%name,4)

               !---------------------------------------------------------------------------!
               !    Copy grid 1, patch 2 to scr3 - Irrelevant if this is the coarsest grid !
               ! on the history file, this will contain the land average if we have inter- !
               ! polated from the coarser grids. We will overwrite points as we inter-     !
               ! polate from finer grids.                                                  !
               !---------------------------------------------------------------------------!
               call patch_land_copy2(nzg,nnxp(1),nnyp(1),npatch,vtab_r(nv,1)%var_p,scr3)

               !----- Then average over history grid patches and interpolate --------------!
               call patch_land_average(nzg1,nnxp1(ngr),nnyp1(ngr),npatch1,parea(:,ngr)     &
                                      ,scr,scr2)
               call hi_interp(nzg1,nnxp1(ngr),nnyp1(ngr),1,scr2,xmn1(:,ngr),xtn1(:,ngr)    &
                             ,ymn1(:,ngr),ytn1(:,ngr),zmn1(:,ngr),ztn1(:,ngr),platn1(ngr)  &
                             ,plonn1(ngr),topt1(:,ngr),ztop1,nzg,nnxp(1),nnyp(1),1,scr3    &
                             ,1,ngr,vtab_r(nv,ngr)%name,4)
               call patch_land_unaverage(nzg,nnxp(ngr),nnyp(ngr),npatch,scr3               &
                                        ,vtab_r(nv,1)%var_p)

            !----- 4D variables (nzs,nxp,nyp,npatch) --------------------------------------!
            elseif (vtab_r(nv,1)%idim_type == 5 ) then
               call hi_interp(nzs1,nnxp1(ngr),nnyp1(ngr),1,scr,xmn1(:,ngr),xtn1(:,ngr)     &
                             ,ymn1(:,ngr),ytn1(:,ngr),zmn1(:,ngr),ztn1(:,ngr),platn1(ngr)  &
                             ,plonn1(ngr),topt1(:,ngr),ztop1,nzs,nnxp(1),nnyp(1),1         &
                             ,vtab_r(nv,1)%var_p,1,ngr,vtab_r(nv,1)%name,5)

               !----- Copy patch 2 to scr3 - This will contain the land average -----------!
               call patch_land_copy2(nzs,nnxp(1),nnyp(1),npatch,vtab_r(nv,1)%var_p,scr3)
               call patch_land_average(nzs1,nnxp1(ngr),nnyp1(ngr),npatch1,parea(:,ngr)     &
                                      ,scr,scr2)
               call hi_interp(nzs1,nnxp1(ngr),nnyp1(ngr),1,scr2,xmn1(1,ngr),xtn1(1,ngr)    &
                             ,ymn1(:,ngr),ytn1(:,ngr),zmn1(:,ngr),ztn1(:,ngr)              &
                             ,platn1(ngr),plonn1(ngr),topt1(:,ngr),ztop1,nzs,nnxp(1)       &
                             ,nnyp(1),1,scr3,1,ngr,vtab_r(nv,1)%name,5)
               call patch_land_unaverage(nzs,nnxp(1),nnyp(1),npatch,scr3                   &
                                        ,vtab_r(nv,1)%var_p)

            !----- 3D variables (nzp,nxp,npatch) that can be interpolated -----------------!
            elseif (vtab_r(nv,1)%idim_type == 6 .and. hr_table(nvh)%string /= 'LEAF_CLASS' &
                   .and. hr_table(nvh)%string /= 'PATCH_AREA'                              &
                   .and. hr_table(nvh)%string /= 'SOIL_COLOR' ) then

               call hi_interp(1,nnxp1(ngr),nnyp1(ngr),1,scr,xmn1(:,ngr),xtn1(:,ngr)        &
                             ,ymn1(:,ngr),ytn1(:,ngr),zmn1(:,ngr),ztn1(:,ngr)              &
                             ,platn1(ngr),plonn1(ngr),topt1(:,ngr),ztop1,1,nnxp(1),nnyp(1) &
                             ,1,vtab_r(nv,1)%var_p,1,ngr,vtab_r(nv,1)%name,6)

               !----- Copy patch 2 to scr3 - This will contain the land average -----------!
               call patch_land_copy2(1,nnxp(1),nnyp(1),npatch,vtab_r(nv,1)%var_p,scr3)

               call patch_land_average(1,nnxp1(ngr),nnyp1(ngr),npatch1,parea(1,ngr),scr    &
                                      ,scr2)
               call hi_interp(1,nnxp1(ngr),nnyp1(ngr),1,scr2,xmn1(:,ngr),xtn1(:,ngr)       &
                             ,ymn1(:,ngr),ytn1(:,ngr),zmn1(:,ngr),ztn1(:,ngr),platn1(ngr)  &
                             ,plonn1(ngr),topt1(:,ngr),ztop1,1,nnxp(1),nnyp(1),1,scr3      &
                             ,1,ngr,vtab_r(nv,1)%name,6)
               call patch_land_unaverage(1,nnxp(1),nnyp(1),npatch,scr3,vtab_r(nv,1)%var_p)

            !----- 3D variables (nxp,nyp,nwave) -------------------------------------------!
            elseif (vtab_r(nv,1)%idim_type == 7) then
               call hi_interp(1,nnxp1(ngr),nnyp1(ngr),nwave,scr,xmn1(:,ngr)                &
                             ,xtn1(:,ngr),ymn1(:,ngr),ytn1(:,ngr),zmn1(:,ngr),ztn1(:,ngr)  &
                             ,platn1(ngr),plonn1(ngr),topt1(:,ngr),ztop1,nnzp(1),nnxp(1)   &
                             ,nnyp(1),nwave,vtab_r(nv,1)%var_p,1,ngr,vtab_r(nv,1)%name,7)

            !----- 4D variables (nzp,nxp,nyp,nclouds) -------------------------------------!
            elseif (vtab_r(nv,1)%idim_type == 8) then
               call hi_interp(nnzp1(ngr),nnxp1(ngr),nnyp1(ngr),nclouds,scr,xmn1(:,ngr)     &
                             ,xtn1(:,ngr),ymn1(:,ngr),ytn1(:,ngr),zmn1(:,ngr),ztn1(:,ngr)  &
                             ,platn1(ngr),plonn1(ngr),topt1(:,ngr),ztop1,nnzp(1),nnxp(1)   &
                             ,nnyp(1),nclouds,vtab_r(nv,1)%var_p,1,ngr,vtab_r(nv,1)%name,8)

            !----- 3D variables (nzp,nxp,nclouds) that can be interpolated ----------------!
            elseif (      vtab_r(nv,1)%idim_type == 9                                      &
                    .and. hr_table(nvh)%string /= 'XIERR'                                  &
                    .and. hr_table(nvh)%string /= 'XJMIN'                                  &
                    .and. hr_table(nvh)%string /= 'XK22'                                   &
                    .and. hr_table(nvh)%string /= 'XKBCON'                                 &
                    .and. hr_table(nvh)%string /= 'XKTOP' ) then
               call hi_interp(1,nnxp1(ngr),nnyp1(ngr),nclouds,scr,xmn1(:,ngr)              &
                             ,xtn1(:,ngr),ymn1(:,ngr),ytn1(:,ngr),zmn1(:,ngr),ztn1(:,ngr)  &
                             ,platn1(ngr),plonn1(ngr),topt1(:,ngr),ztop1,1,nnxp(1),nnyp(1) &
                             ,1,vtab_r(nv,1)%var_p,nclouds,ngr,vtab_r(nv,1)%name,9)
            end if
            exit
         end if
      end do

   end do

   !----- Close the input history file and free some memory -------------------------------!
   close(unit=inhunt)
   deallocate(scr,scr2,scr3,hr_table)


   !-----  Prepare 1D reference sounding --------------------------------------------------!
   nzpg1=nnzp1(1)
   allocate(u01dn1(nzpg1), v01dn1(nzpg1),rt01dn1(nzpg1),th01dn1(nzpg1),pi01dn1(nzpg1)      &
           ,dn01dn1(nzpg1), co201dn1(nzpg1) )

   cng='01'
   ie=cio_f(iunhd,1,'u01dn'//cng  ,  u01dn1,nnzp1(1))
   ie=cio_f(iunhd,1,'v01dn'//cng  ,  v01dn1,nnzp1(1))
   ie=cio_f(iunhd,1,'pi01dn'//cng , pi01dn1,nnzp1(1))
   ie=cio_f(iunhd,1,'th01dn'//cng , th01dn1,nnzp1(1))
   ie=cio_f(iunhd,1,'dn01dn'//cng , dn01dn1,nnzp1(1))
   ie=cio_f(iunhd,1,'rt01dn'//cng , rt01dn1,nnzp1(1))
   ie=cio_f(iunhd,1,'co201dn'//cng,co201dn1,nnzp1(1))


   call htint(nzpg1,th01dn1,ztn1(1,1),nnzp(1),vctr1,ztn(1,1))
   call htint(nzpg1,u01dn1,ztn1(1,1) ,nnzp(1),u01dn(1,1),ztn(1,1))
   call htint(nzpg1,v01dn1,ztn1(1,1) ,nnzp(1),v01dn(1,1),ztn(1,1))

   !----- Assing vapour mixing ratio only if the user is running a "wet" run --------------!
   if (vapour_on) then
      call htint(nzpg1,rt01dn1,ztn1(1,1),nnzp(1),rt01dn(1,1),ztn(1,1))
   else
      rt01dn(1:nnzp(ngrid),1) = 0.
   end if

   !----- Assing CO2 mixing ratio only if the user is running a CO2 run -------------------!
   if (co2_on) then
      call htint(nzpg1,co201dn1,ztn1(1,1),nnzp(1),co201dn(1,1),ztn(1,1))
   else
      co201dn(1:nnzp(ngrid),1) = co2con(1)
   endif

   !----- Saving the virtual potential temperature ----------------------------------------!
   do k = 1,nnzp(ngrid)
      th01dn(k,1) = virtt(vctr1(k),rt01dn(k,1))
   end do

   !----- Lowest level: same as the one just above ----------------------------------------!
   u01dn(1,1)  = u01dn(2,1)
   v01dn(1,1)  = v01dn(2,1)
   rt01dn(1,1) = rt01dn(2,1)
   th01dn(1,1) = th01dn(2,1)
   pi01dn(1,1) = pi01dn1(1) + grav * (ztn1(1,1)-ztn(1,1))                                  &
               / (.5 * (th01dn(1,1) + virtt(th01dn1(1),rt01dn1(1)) ) )

   !----- Computing the ref. Exner function profile, based on hydrostatic equilibrium -----!
   do k = 2,nnzp(1)
      pi01dn(k,1) = pi01dn(k-1,1) - grav                                                   &
                                  / (dzmn(k-1,1)* .5 * (th01dn(k,1) + th01dn(k-1,1)))
   end do

   !----- Computing the ref. density profile, based on the perfect gas law ----------------!
   do k = 1,nnzp(1)
      vctr4(k) = (pi01dn(k,1) / cp) ** cpor * p00
      dn01dn(k,1) = cp * vctr4(k) / (rdry * th01dn(k,1) * pi01dn(k,1))
   end do

   close(unit=iunhd,status='keep')

   !----- Compute 3d reference state for grid 1 -------------------------------------------!
   call newgrid(1)
   call refs3d(nnzp(1),nnxp(1),nnyp(1) , basic_g(1)%pi0    ,basic_g(1)%dn0    &
                                       , basic_g(1)%dn0u   ,basic_g(1)%dn0v   &
                                       , basic_g(1)%th0    ,grid_g(1)%topt    &
                                       , grid_g(1)%rtgt    )

   return
end subroutine inithis
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!   This subroutine initialises the surface variables that are patch-dependent (aka LEAF   !
! arrays). It will also refill many variables, as the interpolated values may not be rele- !
! vant.                                                                                    !
!------------------------------------------------------------------------------------------!
subroutine sfcinit_hstart()

   use mem_leaf
   use mem_basic
   use mem_grid
   use leaf_coms
   use rconstants


   implicit none
   !----- Local variables -----------------------------------------------------------------!
   integer                                         :: i,j,ifm,ipat,k2,k,nveg
   !---------------------------------------------------------------------------------------!

   gridloop: do ifm=1,ngrids

      jloop: do j = 1,nnyp(ifm)
         iloop: do i = 1,nnxp(ifm)

            leaf_g(ifm)%patch_rough(i,j,1) = 0.001

            patchloop: do ipat = 2,npatch

               !---------------------------------------------------------------------------!
               !    Check whether the patch is also water.  If so, assign class 1 (inland  !
               ! water).                                                                   !
               !---------------------------------------------------------------------------!
               nveg = nint(leaf_g(ifm)%leaf_class(i,j,ipat))
               if(nveg == 0) then
                  leaf_g(ifm)%leaf_class(i,j,ipat) = 1.
                  nveg = nint(leaf_g(ifm)%leaf_class(i,j,ipat))
               end if

               !----- Assign roughnesses and a couple of other parameters. ----------------!
               leaf_g(ifm)%soil_rough  (i,j,ipat) = zrough
               leaf_g(ifm)%patch_rough (i,j,ipat) = max(zrough,grid_g(ifm)%topzo(i,j))
               leaf_g(ifm)%veg_rough   (i,j,ipat) = vh2vr * veg_ht(nveg)
               leaf_g(ifm)%veg_height  (i,j,ipat) = veg_ht(nveg)
               leaf_g(ifm)%veg_displace(i,j,ipat) = vh2dh * veg_ht(nveg)
               leaf_g(ifm)%stom_condct (i,j,ipat) = 1.e-6

               !----- Checking for temporary surface water/snow layers --------------------! 
               do k = 1,nzs
                  leaf_g(ifm)%sfcwater_nlev(i,j,ipat) = 0.
                  if (leaf_g(ifm)%sfcwater_mass(k,i,j,ipat) > 0.) then
                     leaf_g(ifm)%sfcwater_nlev(i,j,ipat) = real(k)
                  end if
               end do


               !----- Compute NDVI and LAI. -----------------------------------------------!
               if (ipat >= 2) then
                  call vegndvi(ifm                                                         &
                              ,leaf_g(ifm)%patch_area                  (i,j,ipat)          &
                              ,leaf_g(ifm)%leaf_class                  (i,j,ipat)          &
                              ,leaf_g(ifm)%veg_fracarea                (i,j,ipat)          &
                              ,leaf_g(ifm)%veg_lai                     (i,j,ipat)          &
                              ,leaf_g(ifm)%veg_tai                     (i,j,ipat)          &
                              ,leaf_g(ifm)%veg_rough                   (i,j,ipat)          &
                              ,leaf_g(ifm)%veg_height                  (i,j,ipat)          &
                              ,leaf_g(ifm)%veg_displace                (i,j,ipat)          &
                              ,leaf_g(ifm)%veg_albedo                  (i,j,ipat)          &
                              ,leaf_g(ifm)%veg_ndvip                   (i,j,ipat)          &
                              ,leaf_g(ifm)%veg_ndvic                   (i,j,ipat)          &
                              ,leaf_g(ifm)%veg_ndvif                   (i,j,ipat)          &
                              ,leaf_g(ifm)%psibar_10d                  (i,j,ipat) )
               end if

               !----- Find the surface saturation mixing ratio. ---------------------------!
               call leaf3_grndvap(leaf_g(ifm)%soil_energy           (nzg,i,j,ipat)         &
                                 ,leaf_g(ifm)%soil_water            (nzg,i,j,ipat)         &
                                 ,leaf_g(ifm)%soil_text             (nzg,i,j,ipat)         &
                                 ,leaf_g(ifm)%sfcwater_energy       (nzs,i,j,ipat)         &
                                 ,leaf_g(ifm)%sfcwater_nlev             (i,j,ipat)         &
                                 ,leaf_g(ifm)%can_rvap                  (i,j,ipat)         &
                                 ,leaf_g(ifm)%can_prss                  (i,j,ipat)         &
                                 ,leaf_g(ifm)%ground_rsat               (i,j,ipat)         &
                                 ,leaf_g(ifm)%ground_rvap               (i,j,ipat)         &
                                 ,leaf_g(ifm)%ground_temp               (i,j,ipat)         &
                                 ,leaf_g(ifm)%ground_fliq               (i,j,ipat)         )

            end do patchloop
         end do iloop
      end do jloop

   end do gridloop
   return
end subroutine sfcinit_hstart
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine performs the interpolation of history fields read at the initialisa-  !
! tion. This routine will interpolate the new run's coarse grid only.                      !
!------------------------------------------------------------------------------------------!
subroutine hi_interp(n1,n2,n3,n4,vn,xm1,xt1,ym1,yt1,zm1,zt1,plat1,plon1,topt1,ztop1        &
                    ,m1,m2,m3,m4,vm,ngm,ngr1,vname,idim)

   use mem_grid
   use mem_scratch

   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   character(len=*)               , intent(in)    :: vname
   integer                        , intent(in)    :: n1,n2,n3,n4,ngr1,idim
   integer                        , intent(in)    :: m1,m2,m3,m4,ngm
   real                           , intent(in)    :: plat1,plon1,ztop1
   real   , dimension(n1         ), intent(in)    :: zm1,zt1
   real   , dimension(   n2      ), intent(in)    :: xm1,xt1
   real   , dimension(      n3   ), intent(in)    :: ym1,yt1
   real   , dimension(   n2,n3   ), intent(in)    :: topt1
   real   , dimension(n1,n2,n3,n4), intent(in)    :: vn
   real   , dimension(m1,m2,m3,m4), intent(inout) :: vm
   !----- Local variables -----------------------------------------------------------------!
   integer                               :: i,j,k,np,ii,jj
   real                                  :: xxm,yym,fixxm,fiyym,topoh,rtgth
   real, allocatable, dimension(:,:,:,:) :: scr3
   !---------------------------------------------------------------------------------------!

   allocate(scr3(n2,n3,n1,n4))

   !---------------------------------------------------------------------------------------!
   !    We are going to swap indices due to excessive memory copies of non-contiguous      !
   ! memory.                                                                               !
   !---------------------------------------------------------------------------------------!
   do np=1,n4
      do j=1,n3
         do i=1,n2
            do k=1,n1
               scr3(i,j,k,np) = vn(k,i,j,np)
            end do
         end do
      end do
   end do

   !----- I'm going to interpolate to T points, then average to velocity points -----------!
   do j=1,m3
      do i=1,m2

         !----- Find real grid point x,y relative to history file -------------------------!
         call ll_xy(grid_g(ngm)%glat(i,j),grid_g(ngm)%glon(i,j),plat1,plon1,xxm,yym)

         !----- See if point is on this grid. ---------------------------------------------!
         if (xxm < xm1(1) .or. xxm > xm1(n2-1) .or. yym < ym1(1) .or. yym > ym1(n3-1) ) then 
            if (ngr1 == 1) then

               !---------------------------------------------------------------------------!
               !    We are on the input coarse grid and point is not on this grid. Stop    !
               ! immediately...                                                            !
               !---------------------------------------------------------------------------!
               print*,grid_g(ngm)%glat(i,j),grid_g(ngm)%glon(i,j) &
                    ,plat1,plon1,xxm,yym
               print*,xm1(1),xm1(n2-1),xxm
               print*,ym1(1),ym1(n3-1),yym
               call abort_run('INITHIS: grid point not on history file grids'              &
                             ,'hi_interp','inithis.f90')
            else
               !----- Otherwise, go to next point -----------------------------------------!
               cycle
            end if
         end if

         !---------------------------------------------------------------------------------!
         !    We are okay horizontally, now interpolate vertical column from field.        !
         !    Find x,y grid point locations on input field.                                !
         !    Assuming constant spacing and deal with stagger                              !
         !---------------------------------------------------------------------------------!
         if(trim(vname) == 'UP' .or. trim(vname) == 'UC') then
            fixxm=1.+(xxm-xm1(1))/(xm1(2)-xm1(1))
         else
            fixxm=1.+(xxm-xt1(1))/(xt1(2)-xt1(1))
         end if
         if(trim(vname) == 'VP' .or. trim(vname) == 'VC') then
            fiyym=1.+(yym-ym1(1))/(ym1(2)-ym1(1))
         else
            fiyym=1.+(yym-yt1(1))/(yt1(2)-yt1(1))
         end if

         do np=1,n4
            do k=1,n1
               call gdtost2(scr3(1:n2,1:n3,k,np),n2,n3,fixxm,fiyym,vctr1(k))    
            end do

            select case (idim)
            !----- 3D variable, interpolate this column vertically to actual grid ---------!
            case (3,8)
               
               call gdtost2(topt1(1:n2,1:n3),n2,n3,fixxm,fiyym,topoh)    
               rtgth=1.-topoh/ztop1
               do k=1,m1
                  !----- Actual grid level heights ----------------------------------------!
                  vctr2(k)= grid_g(ngm)%topt(i,j) + ztn(k,1) *grid_g(ngm)%rtgt(i,j)
               end do
               do k=1,n1
                  !----- History grid level heights ---------------------------------------!
                  vctr3(k)= topoh + zt1(k) *rtgth
               end do

               !----- Interpolate vertically ----------------------------------------------!
               call htint(n1,vctr1,vctr3,m1,vctr10,vctr2)
               vm(1:m1,i,j,np)=vctr10(1:m1)
            !----- 2d variable (or 2d+extra dimension), no vertical dimension -------------!
            case (2,6,7,9)
               vm(1,i,j,np)=vctr1(1)
            !---- 3d variable but no atmospheric level, just copying ----------------------!
            case (4,5)
               vm(1:m1,i,j,np)=vctr1(1:m1)
            end select

         end do

      end do
   end do

   !---- Special staggered variables, need averaging --------------------------------------!
   if     (vname == 'UP' .or. vname == 'UC') then
      call hi_avgu(m1,m2,m3,vm)
   elseif (vname == 'VP' .or. vname == 'VC') then
      call hi_avgv(m1,m2,m3,vm)
   elseif (vname == 'WP' .or. vname == 'WC') then
      call hi_avgw(m1,m2,m3,vm)
   endif

   deallocate(scr3)

   return
end subroutine hi_interp
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine performs the averaging for U-grid variables.                          !
!------------------------------------------------------------------------------------------!
subroutine hi_avgu(m1,m2,m3,u)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   integer, intent(in)                         :: m1,m2,m3
   real   , dimension(m1,m2,m3), intent(inout) :: u
   !----- Local variables -----------------------------------------------------------------!
   integer                                     :: i,j,k
   !---------------------------------------------------------------------------------------!

   do k=1,m1
      do j=1,m3
         do i=1,m2-1
            u(k,i,j)=(u(k,i,j)+u(k,i+1,j))*.5
         end do
         u(k,m2,j)=u(k,m2-1,j) !---- Boundary ---------------------------------------------!
      end do
   end do

   return
end subroutine hi_avgu
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine performs the averaging for V-grid variables.                          !
!------------------------------------------------------------------------------------------!
subroutine hi_avgv(m1,m2,m3,v)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   integer, intent(in)                         :: m1,m2,m3
   real   , dimension(m1,m2,m3), intent(inout) :: v
   !----- Local variables -----------------------------------------------------------------!
   integer                                     :: i,j,k
   !---------------------------------------------------------------------------------------!

   do k=1,m1
      do i=1,m2
         do j=1,m3-1
            v(k,i,j)=(v(k,i,j)+v(k,i,j+1))*.5
         end do
         v(k,i,m3)=v(k,i,m3-1) !---- Boundary ---------------------------------------------!
      end do
   end do

  return
end subroutine hi_avgv
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine performs the averaging for W-grid variables.                          !
!------------------------------------------------------------------------------------------!

subroutine hi_avgw(m1,m2,m3,w)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   integer, intent(in)                         :: m1,m2,m3
   real   , dimension(m1,m2,m3), intent(inout) :: w
   !----- Local variables -----------------------------------------------------------------!
   integer                                     :: i,j,k
   !---------------------------------------------------------------------------------------!

   do j=1,m3
      do i=1,m2
         do k=1,m1-1
            w(k,i,j)=(w(k,i,j)+w(k+1,i,j))*.5
         end do
         w(m1,i,j)=w(m1-1,i,j) !---- Boundary ---------------------------------------------!
      end do
   end do

  return
end subroutine hi_avgw
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine computes average of patch variable across all land patches. Recall    !
! water patch is first patch, all others are land.                                         !
!------------------------------------------------------------------------------------------!
subroutine patch_land_average(n1,n2,n3,n4,parea,pa1,paa)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   integer                     , intent(in)    :: n1,n2,n3,n4
   real, dimension(   n2,n3,n4), intent(in)    :: parea
   real, dimension(n1,n2,n3,n4), intent(in)    :: pa1
   real, dimension(n1,n2,n3   ), intent(inout) :: paa
   !----- Local variables -----------------------------------------------------------------!
   integer                                   :: i,j,k,np
   real                                      :: sw
   !---------------------------------------------------------------------------------------!

   write (unit=*,fmt=*) 'PLANDAV:',n1,n2,n3,n4

   do j=1,n3
      do i=1,n2
         do k=1,n1

            sw=0.
            paa(k,i,j)=0.
            do np=2,n4
               paa(k,i,j)= paa(k,i,j)+pa1(k,i,j,np)*parea(i,j,np)
               sw        = sw+parea(i,j,np)
            end do

            if(sw > 0.) then
               paa(k,i,j)=paa(k,i,j)/sw
            else
               paa(k,i,j)=pa1(k,i,j,2)
            end if

         end do
      end do
   end do

   return
end subroutine patch_land_average
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine assigns the average patch variable to all land patches. Recall water  !
! patch is first patch, all others are land.                                               !
!------------------------------------------------------------------------------------------!
subroutine patch_land_unaverage(n1,n2,n3,n4,paa,pa1)

   implicit none

   !----- Arguments -----------------------------------------------------------------------!
   integer                         , intent(in)    :: n1,n2,n3,n4
   real   , dimension (n1,n2,n3)   , intent(in)    :: paa
   real   , dimension (n1,n2,n3,n4), intent(inout) :: pa1
   !----- Local variables -----------------------------------------------------------------!
   integer                                       :: i,j,k,np
   real                                          :: sw
   !---------------------------------------------------------------------------------------!

   do j=1,n3
      do i=1,n2
         do k=1,n1
            do np=2,n4
               pa1(k,i,j,np)=paa(k,i,j)
            end do
         end do
      end do
   end do

   return
end subroutine patch_land_unaverage
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine copies the most common land patch value (patch #2). Recall water      !
! patch is first patch, all others are land.                                               !
!------------------------------------------------------------------------------------------!
subroutine patch_land_copy2(n1,n2,n3,n4,pa1,paa)

   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   integer                         , intent(in)  :: n1,n2,n3,n4
   real   , dimension (n1,n2,n3,n4), intent(in)  :: pa1
   real   , dimension (n1,n2,n3)   , intent(out) :: paa
   !----- Local variables -----------------------------------------------------------------!
   integer                                       :: i,j,k,np
   real                                          :: sw
   !---------------------------------------------------------------------------------------!

   do j=1,n3
      do i=1,n2
         do k=1,n1
            paa(k,i,j)=pa1(k,i,j,2)
         end do
      end do
   end do

   return
end subroutine patch_land_copy2
!==========================================================================================!
!==========================================================================================!
