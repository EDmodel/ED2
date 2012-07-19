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
!      This sub-routine initialises several parameters that define a simulation.           !
!------------------------------------------------------------------------------------------!
subroutine initlz (name_name)
   use mem_leaf
   use mem_grid
   use mem_scratch
   use mem_basic
   use mem_micro
   use var_tables
   use mem_varinit
   use mem_cuparm
   use mem_oda
   use io_params
   use micphys
   use therm_lib          , only : level
   use mem_turb           , only : turb_g                    & ! intent(in)
                                 , if_urban_canopy           ! ! intent(out)
   use catt_start         , only : catt                      ! ! intent(in)
   use mem_scalar         , only : scalar_g                  & ! intent(inout)
                                 , recycle_tracers           ! ! intent(in)
   use emission_source_map, only : read_emission_sources_map ! ! subroutine
   use teb_spm_start      , only : TEB_SPM                   ! ! intent(in)
   use mem_teb            , only : teb_g                     ! ! intent(inout)
   use mem_teb_common     , only : tebc_g                    ! ! intent(inout)
   use mem_mnt_advec      , only : iadvec                    ! ! intent(inout)
   use teb_vars_const     , only : iteb                      ! ! intent(in)
   use mem_gaspart        , only : gaspart_g                 ! ! intent(inout)
   use mem_emiss          , only : ichemi                    & ! intent(inout)
                                 , isource                   & ! intent(inout)
                                 , ichemi_in                 ! ! intent(inout)
   use mem_soil_moisture  , only : soil_moist                ! ! intent(in)
   use ref_sounding       , only : dn01dn                    ! ! intent(inout)

   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   character(len=*), intent(in) :: name_name
   !----- Local variables. ----------------------------------------------------------------!
   character(len=8)             :: rest
   integer                      :: ifm
   integer                      :: icm
   integer                      :: mynum
   integer                      :: ihm
   integer                      :: ngr
   integer                      :: nv
   integer                      :: ierr
   !---------------------------------------------------------------------------------------!



   !----- Set version number for common blocks. -------------------------------------------!
   iversion = 2
   !---------------------------------------------------------------------------------------!

   !----- Unit numbers for some I/O files. ------------------------------------------------!
   iopunt=6
   !---------------------------------------------------------------------------------------!




   select case (trim(runtype))
   case ('INITIAL')
      !------------------------------------------------------------------------------------!
      !       Initial start up.                                                            !
      !------------------------------------------------------------------------------------!

      time              = 0.d0
      ngbegun(1:ngrids) = 0


      !------------------------------------------------------------------------------------!
      !      Read surface, topo, sst, and ndvi files for all grids.  All the files were    !
      ! checked earlier, so they must be correct.                                          !
      !------------------------------------------------------------------------------------!
      do ifm = 1,ngrids
         call top_read(ifm)
      end do

      do ifm = 1,ngrids
         call sfc_read(ifm)
      end do

      !------------------------------------------------------------------------------------!
      !     Define grid topography, transform, latitude-longitude, and map factor arrays.  !
      !------------------------------------------------------------------------------------!
      call grid_setup(2)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Read in the SST files.                                                        !
      !------------------------------------------------------------------------------------!
      do ifm = 1,ngrids
         call sst_read(1,ifm,ierr)
      end do
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Read in the NDVI/LAI files.                                                   !
      !------------------------------------------------------------------------------------!
      do ifm = 1,ngrids
         call ndvi_read(1,ifm,ierr)
      end do
      !------------------------------------------------------------------------------------!



      !----- Initialize snowcover arrays. -------------------------------------------------!
      do ifm = 1,ngrids
         call snowinit(nnxp(ifm),nnyp(ifm),leaf_g(ifm)%snow_mass,leaf_g(ifm)%snow_depth)
      end do
      !------------------------------------------------------------------------------------!



      !----- If TEB is used, read in the local time zone (fuso) files. --------------------!
      if (TEB_SPM==1) then
         do ifm = 1,ngrids
            call fuso_read(ifm)
         end do
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !       Check which kind of initialisation we should use, based on the INITIAL       !
      ! value:                                                                             !
      !                                                                                    !
      ! 1 -- Horizontal homogenous ATM                                                     !
      ! 2 -- Varfile initialisation (non-BRAMS analysis or reanalysis, in RALPH2 format).  !
      ! 3 -- HFILIN initialisation (RAMS/BRAMS history file)                               !
      !------------------------------------------------------------------------------------!
      select case (initial)
      case (1,3)

         select case(initial)
         case (1)
            !------------------------------------------------------------------------------!
            !    If this is a horizontally homogeneous start up, call subroutine INITHH,   !
            ! which loops through all grids and initializes those for which nxtnest = 0.   !
            !------------------------------------------------------------------------------!
            write (unit=*,fmt='(a)') 'Horizontally-homogeneous-INITIAL start of grid- 1'
            call inithh()
         case (3)
            !------------------------------------------------------------------------------!
            !    If "history" initialization, call INITHIS.  This will define initial      !
            ! fields and reference state on grid 1 from history file.  Other grids will be !
            ! interpolated as in a INITIAL = 1 start.                                      !
            !------------------------------------------------------------------------------!
            write (unit=*,fmt='(a)') 'History-INITIAL start of grid- 1'
            call inithis()
         end select


         !---------------------------------------------------------------------------------!
         !     On all fine grids, initialize the surface layer characteristics, the 1-D    !
         ! reference state arrays, the 3-D reference state arrays, and the prognostic      !
         ! atmospheric fields by interpolation.                                            !
         !---------------------------------------------------------------------------------!
         call fmrefs1d(2,ngrids)
         do ifm = 2,ngrids
            icm = nxtnest(ifm)
            if (icm  >=  1) then
               call fmrefs3d(ifm,mynum)
               call prgintrp(nnzp(icm),nnxp(icm),nnyp(icm),nnzp(icm),nnxp(icm),nnyp(icm)   &
                            ,0,0,ifm,1,mynum)
               call fmdn0(ifm)

               write(unit=*,fmt='(a,1x,i6,a)') '  Initial interpolation of grid',ifm,'...'
            end if
         end do


      case (2)
         !---------------------------------------------------------------------------------!
         !     If this is a RALPH2 initialisation, call the driver that will do all the    !
         ! steps.                                                                          !
         !---------------------------------------------------------------------------------!
         call varf_read(0)
         !---------------------------------------------------------------------------------!
      end select

      !------------------------------------------------------------------------------------!
      !     Initialise past time level velocity and the perturbation of the Exner function !
      ! on all grids.                                                                      !
      !------------------------------------------------------------------------------------!
      do ifm=1,ngrids
         call newgrid(ifm)
         call fldinit(1)
         call negadj1(nzp,nxp,nyp,1,nxp,1,nyp)
         call thermo(nzp,nxp,nyp,1,nxp,1,nyp)

         if (level  ==  3) then
            call azero(nzp*nxp*nyp,scratch%vt3da)
            call azero(nzp*nxp*nyp,scratch%vt3dg)
            call azero(nzp*nxp*nyp,scratch%vt3dh)
            call azero(nzp*nxp*nyp,scratch%vt3dc)
            call azero(nzp*nxp*nyp,scratch%vt3di)
            !----- Use scratch variables to define cccnp and cifnp ------------------------!
            call initqin(nzp,nxp,nyp,scratch%vt3da,scratch%vt3dg,scratch%vt3dh             &
                        ,basic_g(ifm)%pi0,basic_g(ifm)%pp,basic_g(ifm)%theta               &
                        ,basic_g(ifm)%dn0,scratch%vt3dc,scratch%vt3di        )
            !----- Copying them to the micro arrays if they are allocated -----------------!
            if (irain  >= 1) call atob(nzp*nxp*nyp,scratch%vt3da,micro_g(ifm)%q2)
            if (igraup >= 1) call atob(nzp*nxp*nyp,scratch%vt3dg,micro_g(ifm)%q6)
            if (ihail  >= 1) call atob(nzp*nxp*nyp,scratch%vt3dh,micro_g(ifm)%q7)
            if (icloud == 7) call atob(nzp*nxp*nyp,scratch%vt3dc,micro_g(ifm)%cccnp)
            if (ipris  == 7) call atob(nzp*nxp*nyp,scratch%vt3di,micro_g(ifm)%cifnp)
         end if
      end do
      !------------------------------------------------------------------------------------!



      !------ If initializing some fields from previous runs... ---------------------------!


      !------------------------------------------------------------------------------------!
      !      Run recycle procedure for normal RAMS recycle or only to perform the          !
      ! assimilation of the tracers.                                                       !
      !------------------------------------------------------------------------------------!
      if ((ipastin == 1) .or. (catt==1 .and. recycle_tracers==1)) then
         call recycle()
      end if
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Fill land surface data for all grids that have no standard input files.        !
      !------------------------------------------------------------------------------------!
      select case (isfcl)
      case (1,2,5)
         call sfcdata
      case (3)
         !call sfcdata_sib_driver
      end select
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Initialize various LEAF variables.                                             !
      !------------------------------------------------------------------------------------!
      if (ipastin == 0) call geonest_nofile(1,ngrids)
      !------------------------------------------------------------------------------------!




      if (teb_spm == 1) then
         !---------------------------------------------------------------------------------!
         !      Initial values for Common use TEB vars.  For the time being, it is only    !
         ! being used to flush these four commonly used variables to zero.  These vari-    !
         ! ables will be used for other purposes later.                                    !
         ! Edmilson D. Freitas 07/07/2006                                                  !
         !---------------------------------------------------------------------------------!
         do ifm = 1,ngrids
            call tebc_init(ifm,nnxp(ifm),nnyp(ifm),npatch,leaf_g(ifm)%leaf_class           &
                          ,leaf_g(ifm)%g_urban,tebc_g(ifm)%emis_town,tebc_g(ifm)%alb_town  &
                          ,tebc_g(ifm)%ts_town)
         end do
         !---------------------------------------------------------------------------------!

         if (iteb==1) then
            do ifm = 1,ngrids
               call teb_init(ifm,nnzp(ifm),nnxp(ifm),nnyp(ifm),npatch                      &
                    ,leaf_g(ifm)%leaf_class        ,basic_g(ifm)%theta                     &
                    ,basic_g(ifm)%rv               ,basic_g(ifm)%pi0                       &
                    ,basic_g(ifm)%pp                                                       &
                    ,teb_g(ifm)%t_roof             ,teb_g(ifm)%t_road                      &
                    ,teb_g(ifm)%t_wall             ,teb_g(ifm)%ti_bld                      &
                    ,teb_g(ifm)%ti_road            ,teb_g(ifm)%t_canyon                    &
                    ,teb_g(ifm)%r_canyon           ,teb_g(ifm)%ts_roof                     &
                    ,teb_g(ifm)%ts_road            ,teb_g(ifm)%ts_wall                     &
                    ,teb_g(ifm)%h_traffic          ,teb_g(ifm)%le_traffic                  &
                    ,teb_g(ifm)%h_industry         ,teb_g(ifm)%le_industry                 &
                    ,teb_g(ifm)%ws_roof            ,teb_g(ifm)%ws_road                     &
                    ,tebc_g(ifm)%emis_town         ,tebc_g(ifm)%alb_town                   &
                    ,tebc_g(ifm)%ts_town           ,leaf_g(ifm)%g_urban                    )
            end do
         end if

         !----- Initialize gases and particulate matter. ----------------------------------!
         if (isource==1) then
            do ifm = 1,ngrids
               call init_conc1(1,ifm,nnzp(ifm),nnxp(ifm),nnyp(ifm),npatch                  &
                    ,leaf_g(ifm)%G_URBAN         ,gaspart_g(ifm)%pno                       &
                    ,gaspart_g(ifm)%pno2         ,gaspart_g(ifm)%ppm25                     &
                    ,gaspart_g(ifm)%pco          ,gaspart_g(ifm)%pvoc                      &
                    ,gaspart_g(ifm)%pso2         ,gaspart_g(ifm)%pso4                      &
                    ,gaspart_g(ifm)%paer         ,zt                                       )
               if (ichemi==1) then       !call more added scalars for chemistry
                  if (ichemi_in==1) then !read init.values from a previous run
                     call init_conc_prev(name_name)
                  else
                     call init_conc2(1,ifm,nnzp(ifm),nnxp(ifm),nnyp(ifm),npatch            &
                          ,leaf_g(ifm)%G_URBAN    ,gaspart_g(ifm)%po3                      &
                          ,gaspart_g(ifm)%prhco   ,gaspart_g(ifm)%pho2                     &
                          ,gaspart_g(ifm)%po3p    ,gaspart_g(ifm)%po1d                     &
                          ,gaspart_g(ifm)%pho     ,gaspart_g(ifm)%proo                     &
                          ,zt)
                  end if
               end if
            end do
         end if

      end if

      !----- Initialise turbulence factor akscal. -----------------------------------------!
      if (if_adap == 1) then
         do ifm=1,ngrids
            call akscal_init(nnxp(ifm),nnyp(ifm),ifm,grid_g(ifm)%topta,turb_g(ifm)%akscal)
         end do
      else
         do ifm=1,ngrids
            call akscal_init(nnxp(ifm),nnyp(ifm),ifm,grid_g(ifm)%topt,turb_g(ifm)%akscal)
         end do
      end if

      !----- If CATT is to be used, initialise tracer concentration and sources. ----------!
      if (catt==1) then
         if (recycle_tracers==0) then
            do ifm=1,ngrids
               !---------------------------------------------------------------------------!
               !     We shan't call the recycle procedure for tracers, instead we          !
               ! initialise them with zeroes.                                              !
               !---------------------------------------------------------------------------!
               scalar_g(1,ifm)%sclp(:,:,:) = 0.
               scalar_g(2,ifm)%sclp(:,:,:) = 0.
               scalar_g(3,ifm)%sclp(:,:,:) = 0.
               scalar_g(4,ifm)%sclp(:,:,:) = 0.
               !---------------------------------------------------------------------------!
            end do
         end if

         !----- Read emission maps for all grids. -----------------------------------------!
         call read_emission_sources_map()

      end if
      !----- Call the history start up driver for surface values. -------------------------!
      if (initial == 3) call sfcinit_hstart()

   case ('HISTORY')
      !------------------------------------------------------------------------------------!
      !       History initialisation.                                                      !
      !------------------------------------------------------------------------------------!
      call history_start(name_name)

      call grid_setup(1)

      !------------------------------------------------------------------------------------!
      !     Check surface, topo, sst, and ndvi files.  In case they are missing or don't   !
      ! match the current grid set up, re-create them.                                     !
      !------------------------------------------------------------------------------------!
      call make_sfcfiles()


      !----- Read surface and topo files for any new grid. --------------------------------!
      do ifm = ngridsh+1,ngrids
         call sfc_read(ifm)
      end do
      do ifm = ngridsh+1,ngrids
         call top_read(ifm)
      end do

      call grid_setup(2)


      !----- Read in SST and NDVI/LAI files for all grids. --------------------------------!
      do ifm = 1,ngrids
         call sst_read(1,ifm,ierr)
      end do
      do ifm = 1,ngrids
         call ndvi_read(1,ifm,ierr)
      end do

      !----- If TEB is used, read in the local time zone (fuso) files. --------------------!
      if (TEB_SPM==1) then
         do ifm = ngridsh+1, ngrids
            call fuso_read(ifm)
         end do
      end if
      !------------------------------------------------------------------------------------!


       !------ Initialise the 3-D reference states for density and the Exner function. ----!
       do ifm = 1,ngrids
         icm = nxtnest(ifm)
         if (icm  ==  0) then
            call newgrid(ifm)
            call refs3d (nzp,nxp,nyp,basic_g(ifm)%pi0,basic_g(ifm)%dn0,basic_g(ifm)%dn0u   &
                        ,basic_g(ifm)%dn0v,basic_g(ifm)%th0,grid_g(ifm)%topt               &
                        ,grid_g(ifm)%rtgt)
         end if
      end do
      !------------------------------------------------------------------------------------!

      do ifm = 1,min(ngrids,ngridsh)
         icm = nxtnest(ifm)
         if (icm  >  0) call fmrefs3d(ifm,0)
         call negadj1(nzp,nxp,nyp,1,nxp,1,nyp)
      enddo
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Fill land surface data for all grids that have no standard input files.        !
      !------------------------------------------------------------------------------------!
        select case (isfcl)
      case (1,2,5)
         call sfcdata
      case (3)
         !call sfcdata_sib_driver
      end select
      !------------------------------------------------------------------------------------!




      !----- Heterogenous soil moisture initialisation. -----------------------------------!
      select case (trim(soil_moist))
      case ('h','H','a','A')

         do ifm = 1,min(ngrids,ngridsh)
            call newgrid(ifm)
            call soil_moisture_init(nnzp(ifm),nnxp(ifm),nnyp(ifm),nzg,npatch,ifm           &
                                   ,leaf_g(ifm)%can_theta,leaf_g(ifm)%can_prss             &
                                   ,grid_g(ifm)%glat,grid_g(ifm)%glon                      &
                                   ,leaf_g(ifm)%soil_water,leaf_g(ifm)%soil_energy         &
                                   ,leaf_g(ifm)%soil_text,leaf_g(ifm)%psibar_10d           &
                                   ,leaf_g(ifm)%leaf_class)
         end do
      end select
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     If any grids are being added for this run, initialize their surface layer      !
      ! variables, 1-D reference state variables, and prognostic atmospheric and soil      !
      ! model fields.                                                                      !
      !------------------------------------------------------------------------------------!
      if (ngrids  >  ngridsh) then
         write (unit=*,fmt='(a)')        '------------------------------------------------'
         write (unit=*,fmt='(a)')        ' New grids will be added. '
         write (unit=*,fmt='(a)')        '------------------------------------------------'
         write (unit=*,fmt='(i6,1x,a)')  ngridsh,'grid(s) on history file.'
         write (unit=*,fmt='(i6,1x,a)')  ngrids ,'grids to be run.'
         write (unit=*,fmt='(a)')        '------------------------------------------------'
         call fmrefs1d(ngridsh+1,ngrids)
         do ifm = ngridsh+1,ngrids
            icm = nxtnest(ifm)
            if (icm  ==  0) then
               call abort_run('Attempted to add hemispheric grid on a history restart'     &
                             ,'initlz','rdint.f90')
            endif
            call fmrefs3d(ifm,0)
            call prgintrp(nnzp(icm),nnxp(icm),nnyp(icm),nnzp(icm),nnxp(icm),nnyp(icm),0,0  &
                         ,ifm,1,mynum)
            write(unit=*,fmt='(a,1x,i6,a)') 'History start interpolation of added grid '   &
                                           ,ngrid,'...'
            call fmdn0(ifm)
            call newgrid(ifm)
            call fldinit(0)
            call negadj1(nzp,nxp,nyp,1,nxp,1,nyp)
            call thermo(nzp,nxp,nyp,1,nxp,1,nyp)
            if (level  ==  3) then
               call azero(nzp*nxp*nyp,scratch%vt3da)
               call azero(nzp*nxp*nyp,scratch%vt3dg)
               call azero(nzp*nxp*nyp,scratch%vt3dh)
               call azero(nzp*nxp*nyp,scratch%vt3dc)
               call azero(nzp*nxp*nyp,scratch%vt3di)
               !----- Use scratch variables to define cccnp and cifnp ---------------------!
               call initqin(nzp,nxp,nyp,scratch%vt3da,scratch%vt3dg,scratch%vt3dh          &
                           ,basic_g(ifm)%pi0,basic_g(ifm)%pp,basic_g(ifm)%theta            &
                           ,basic_g(ifm)%dn0,scratch%vt3dc,scratch%vt3di        )
               !----- Copying them to the micro arrays if they are allocated --------------!
               if (irain  >= 1) call atob(nzp*nxp*nyp,scratch%vt3da,micro_g(ifm)%q2)
               if (igraup >= 1) call atob(nzp*nxp*nyp,scratch%vt3dg,micro_g(ifm)%q6)
               if (ihail  >= 1) call atob(nzp*nxp*nyp,scratch%vt3dh,micro_g(ifm)%q7)
               if (icloud == 7) call atob(nzp*nxp*nyp,scratch%vt3dc,micro_g(ifm)%cccnp)
               if (ipris  == 7) call atob(nzp*nxp*nyp,scratch%vt3di,micro_g(ifm)%cifnp)
            end if

            !----- Heterogenous Soil Moisture Initialisation. -----------------------------!
            select case (trim(soil_moist))
            case ('h','H','a','A')
               call soil_moisture_init(nnzp(ifm),nnxp(ifm),nnyp(ifm),nzg,npatch,ifm        &
                                      ,leaf_g(ifm)%can_theta    ,leaf_g(ifm)%can_prss      &
                                      ,grid_g(ifm)%glat         ,grid_g(ifm)%glon          &
                                      ,leaf_g(ifm)%soil_water   ,leaf_g(ifm)%soil_energy   &
                                      ,leaf_g(ifm)%soil_text    ,leaf_g(ifm)%psibar_10d    &
                                      ,leaf_g(ifm)%leaf_class   )
            end select
            !------------------------------------------------------------------------------!

            !----- Initialise turbulence factor akscal. -----------------------------------!
            if (if_adap == 1) then
               call akscal_init(nnxp(ifm),nnyp(ifm),ifm,grid_g(ifm)%topta                  &
                               ,turb_g(ifm)%akscal)
            else
               call akscal_init(nnxp(ifm),nnyp(ifm),ifm,grid_g(ifm)%topt                   &
                               ,turb_g(ifm)%akscal)
            end if
            !------------------------------------------------------------------------------!
         end do

         !----- Fill land surface data for all grids that have no standard input files. ---!
         call geonest_nofile(ngridsh+1,ngrids)

      elseif (ngrids  <  ngridsh) then
         write (unit=*,fmt='(a)')        '------------------------------------------------'
         write (unit=*,fmt='(a)')        ' Grids will be dismissed. '
         write (unit=*,fmt='(a)')        '------------------------------------------------'
         write (unit=*,fmt='(i6,1x,a)')  ngridsh,'grid(s) on history file.'
         write (unit=*,fmt='(i6,1x,a)')  ngrids ,'grids to be run.'
         write (unit=*,fmt='(a)')        '------------------------------------------------'
      end if

      !----- If CATT is turned on, read in emission maps for all grids. -------------------!
      if (CATT==1) then
         call read_emission_sources_map()
      end if

   case default
      call abort_run('Invalid runtype in initlz: '//trim(runtype)//'...'                   &
                    ,'initlz','rdint.f90')

   end select

   !---------------------------------------------------------------------------------------!
   !      For a global model domain, initialize interpolation table values for             !
   ! communication between hemispheric grids.  Subroutine hemintrp_cof will return         !
   ! immediately if nhemgrd2 is not greater than 1.                                        !
   !---------------------------------------------------------------------------------------!
   call newgrid(1)
   ihm = nnxyp(1)
   call hemintrp_cof (nnxp(1),nnyp(1),scratch%scr1,scratch%scr2,scratch%scr3,scratch%scr4  &
                     ,scratch%scr5,scratch%scr6)

   !---------------------------------------------------------------------------------------!
   !    Initialise a bunch of microphysics parameters. The new input is the density at the !
   ! top, to find the minimum concentration of hydrometeors to be considered.              !
   !---------------------------------------------------------------------------------------!
   call micro_master(dn01dn(nnzp(1),1))
   !---------------------------------------------------------------------------------------!



   !----- Fill latitude-longitude, map factor, and Coriolis arrays. -----------------------!
   do ifm = 1,ngrids
      call newgrid(ifm)
      call fcorio(nxp,nyp,basic_g(ifm)%fcoru,basic_g(ifm)%fcorv,grid_g(ifm)%glat)
   end do
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     If we are doing one-way nesting or varfile nudging, inventory, prepare history/   !
   ! varfile files and fill past/future nudging arrays for start of simulation.            !
   !---------------------------------------------------------------------------------------!
   select case (nud_type)
   case (1)
      call nud_read(1)
   case (2)
      call varf_read(1)
   end select
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Same thing if we are nudging hydrometeors.                                        !
   !---------------------------------------------------------------------------------------!
   if (nud_cond == 1) call cond_read(1)

   !----- Process and read observations for ODA - observational data assimilation. --------!
   if (if_oda == 1) call oda_read()

   !----- Read in cumulus heating fields. -------------------------------------------------!
   if (if_cuinv == 1) call cu_read(1)

   !----- Initialize urban canopy drag coefficients. --------------------------------------!
   if (if_urban_canopy == 1) call urb_drag_init()

   !----- Print locations of all grids. ---------------------------------------------------!
   call gridloc_prt()
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Save initial fields on history and analysis files.                                !
   !---------------------------------------------------------------------------------------!
   iinput = ioutput
   rest   = 'no'
   if (trim(runtype)  ==  'HISTORY') rest = 'yes'
   call hiswrt(rest)
   select case (ioutput)
   case (3)
      call anlhdf('INST')
   case default
      call anlwrt(rest,'INST')
   end select

   if (frqlite > 0.) then
      select case(ioutput)
      case (3)
         call anlhdf('LITE')
      case default
         call anlwrt(rest,'LITE')
      end select
   end if



   !---------------------------------------------------------------------------------------!
   !     Save initial fields on history and into the averaged arrays.                      !
   !---------------------------------------------------------------------------------------!
   if (avgtim /= 0.) then
      do ngr=1,ngrids
         do nv=1,num_var(ngr)
            if (vtab_r(nv,ngr)%imean == 1) then
               call atob(vtab_r(nv,ngr)%npts,vtab_r(nv,ngr)%var_p,vtab_r(nv,ngr)%var_m)
            end if
         end do
      end do
   end if
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Print the header information and initial fields.                                  !
   !---------------------------------------------------------------------------------------!
   ngrid=1
   call prtopt(6)

   if (initfld  ==  1) then
      do ifm = 1,ngrids
         call newgrid(ifm)
         call prtout()
      end do
   end if
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Check whether the settings make sense or not.                                     !
   !---------------------------------------------------------------------------------------!
   call opspec3()

   return
end subroutine initlz
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This sub-routine reads the namelist, and check whether some of the options make      !
! sense or not.                                                                            !
!------------------------------------------------------------------------------------------!
subroutine read_nl(filename)
   use io_params          , only : frqboth                 & ! intent(out)
                                 , afilout                 & ! intent(out)
                                 , avgtim                  & ! intent(out)
                                 , frqanl                  & ! intent(out)
                                 , frqhis                  & ! intent(out)
                                 , frqlite                 & ! intent(out)
                                 , frqmean                 & ! intent(out)
                                 , frqprt                  & ! intent(out)
                                 , hfilin                  & ! intent(out)
                                 , hfilout                 & ! intent(out)
                                 , iclobber                & ! intent(out)
                                 , ihistdel                & ! intent(out)
                                 , initfld                 & ! intent(out)
                                 , ioutput                 & ! intent(out)
                                 , ipastin                 & ! intent(out)
                                 , iplfld                  & ! intent(out)
                                 , isbval                  & ! intent(out)
                                 , isoilflg                & ! intent(out)
                                 , isoilfn                 & ! intent(out)
                                 , isstflg                 & ! intent(out)
                                 , isstfn                  & ! intent(out)
                                 , itopsflg                & ! intent(out)
                                 , itoptflg                & ! intent(out)
                                 , itoptfn                 & ! intent(out)
                                 , iupdndvi                & ! intent(out)
                                 , iupdsst                 & ! intent(out)
                                 , iuselai                 & ! intent(out)
                                 , ivegtflg                & ! intent(out)
                                 , ivegtfn                 & ! intent(out)
                                 , ixsctn                  & ! intent(out)
                                 , iz0flg                  & ! intent(out)
                                 , kwrite                  & ! intent(out)
                                 , lite_vars               & ! intent(out)
                                 , ndviflg                 & ! intent(out)
                                 , ndvifn                  & ! intent(out)
                                 , ndvifpfx                & ! intent(out)
                                 , nlite_vars              & ! intent(out)
                                 , nofilflg                & ! intent(out)
                                 , nplt                    & ! intent(out)
                                 , pastfn                  & ! intent(out)
                                 , sfcfiles                & ! intent(out)
                                 , sstfpfx                 & ! intent(out)
                                 , timstr                  & ! intent(out)
                                 , imonthh                 & ! intent(out)
                                 , iyearh                  & ! intent(out)
                                 , idateh                  & ! intent(out)
                                 , itimeh                  & ! intent(out)
                                 , topfiles                & ! intent(out)
                                 , toptenh                 & ! intent(out)
                                 , toptwvl                 & ! intent(out)
                                 , xlite                   & ! intent(out)
                                 , ylite                   & ! intent(out)
                                 , z0fact                  & ! intent(out)
                                 , z0max                   & ! intent(out)
                                 , zlite                   & ! intent(out)
                                 , ifusflg                 & ! intent(out)
                                 , ifusfn                  & ! intent(out)
                                 , fusfiles                ! ! intent(out)
   use isan_coms          , only : gobrad                  & ! intent(out)
                                 , gobsep                  & ! intent(out)
                                 , gridwt                  & ! intent(out)
                                 , guess1st                & ! intent(out)
                                 , hybbot                  & ! intent(out)
                                 , hybtop                  & ! intent(out)
                                 , i1st_flg                & ! intent(out)
                                 , iapr                    & ! intent(out)
                                 , iarawi                  & ! intent(out)
                                 , iasrfce                 & ! intent(out)
                                 , igridfl                 & ! intent(out)
                                 , iobswin                 & ! intent(out)
                                 , ioflgisz                & ! intent(out)
                                 , ioflgvar                & ! intent(out)
                                 , isan_inc                & ! intent(out)
                                 , isfc_flg                & ! intent(out)
                                 , iszstage                & ! intent(out)
                                 , iupa_flg                & ! intent(out)
                                 , ivrstage                & ! intent(out)
                                 , levth                   & ! intent(out)
                                 , maxsfc                  & ! intent(out)
                                 , maxsta                  & ! intent(out)
                                 , nfeedvar                & ! intent(out)
                                 , nigrids                 & ! intent(out)
                                 , nisn                    & ! intent(out)
                                 , notid                   & ! intent(out)
                                 , notsta                  & ! intent(out)
                                 , respon                  & ! intent(out)
                                 , sfcinf                  & ! intent(out)
                                 , sigzwt                  & ! intent(out)
                                 , stasep                  & ! intent(out)
                                 , swvlnth                 & ! intent(out)
                                 , topsigz                 & ! intent(out)
                                 , varpfx                  & ! intent(out)
                                 , wvlnth                  ! ! intent(out)
   use mem_cuparm         , only : confrq                  & ! intent(out)
                                 , cu_prefix               & ! intent(out)
                                 , cu_tel                  & ! intent(out)
                                 , cu_til                  & ! intent(out)
                                 , if_cuinv                & ! intent(out)
                                 , nnqparm                 & ! intent(out)
                                 , tcu_beg                 & ! intent(out)
                                 , tcu_end                 & ! intent(out)
                                 , tnudcu                  & ! intent(out)
                                 , wcldbs                  & ! intent(out)
                                 , wt_cu_grid              & ! intent(out)
                                 , nclouds                 & ! intent(out)
                                 , ndeepest                & ! intent(out)
                                 , nshallowest             & ! intent(out)
                                 , cptime                  ! ! intent(out)
   use grell_coms         , only : closure_type            & ! intent(out)
                                 , maxclouds               & ! intent(out)
                                 , iupmethod               & ! intent(out)
                                 , depth_min               & ! intent(out)
                                 , cap_maxs                & ! intent(out)
                                 , cld2prec                & ! intent(out)
                                 , maxens_lsf              & ! intent(out)
                                 , maxens_eff              & ! intent(out)
                                 , maxens_dyn              & ! intent(out)
                                 , maxens_cap              & ! intent(out)
                                 , iupmethod               & ! intent(out)
                                 , radius                  & ! intent(out)
                                 , zkbmax                  & ! intent(out)
                                 , max_heat                & ! intent(out)
                                 , zcutdown                & ! intent(out)
                                 , z_detr                  ! ! intent(out)
   use mem_grid ,           only : centlat                 & ! intent(out)
                                 , centlon                 & ! intent(out)
                                 , cphas                   & ! intent(out)
                                 , deltax                  & ! intent(out)
                                 , deltay                  & ! intent(out)
                                 , deltaz                  & ! intent(out)
                                 , distim                  & ! intent(out)
                                 , dtlong                  & ! intent(out)
                                 , dzmax                   & ! intent(out)
                                 , dzrat                   & ! intent(out)
                                 , expnme                  & ! intent(out)
                                 , gridu                   & ! intent(out)
                                 , gridv                   & ! intent(out)
                                 , ibnd                    & ! intent(out)
                                 , icorflg                 & ! intent(out)
                                 , idatea                  & ! intent(out)
                                 , idatez                  & ! intent(out)
                                 , ideltat                 & ! intent(out)
                                 , if_adap                 & ! intent(out)
                                 , ihtran                  & ! intent(out)
                                 , initial                 & ! intent(out)
                                 , imontha                 & ! intent(out)
                                 , imonthz                 & ! intent(out)
                                 , itimea                  & ! intent(out)
                                 , ihoura                  & ! intent(out)
                                 , itimez                  & ! intent(out)
                                 , iyeara                  & ! intent(out)
                                 , iyearz                  & ! intent(out)
                                 , jbnd                    & ! intent(out)
                                 , lsflg                   & ! intent(out)
                                 , nacoust                 & ! intent(out)
                                 , naddsc                  & ! intent(out)
                                 , nestz1                  & ! intent(out)
                                 , nestz2                  & ! intent(out)
                                 , nfpt                    & ! intent(out)
                                 , ngrids                  & ! intent(out)
                                 , ninest                  & ! intent(out)
                                 , njnest                  & ! intent(out)
                                 , nknest                  & ! intent(out)
                                 , nndtrat                 & ! intent(out)
                                 , nnstbot                 & ! intent(out)
                                 , nnsttop                 & ! intent(out)
                                 , nnxp                    & ! intent(out)
                                 , nnyp                    & ! intent(out)
                                 , nnzp                    & ! intent(out)
                                 , npatch                  & ! intent(out)
                                 , nstratx                 & ! intent(out)
                                 , nstraty                 & ! intent(out)
                                 , nstratz1                & ! intent(out)
                                 , nstratz2                & ! intent(out)
                                 , nxtnest                 & ! intent(out)
                                 , nzg                     & ! intent(out)
                                 , nzs                     & ! intent(out)
                                 , polelat                 & ! intent(out)
                                 , polelon                 & ! intent(out)
                                 , runtype                 & ! intent(out)
                                 , timeunit                & ! intent(out)
                                 , timmax                  & ! intent(out)
                                 , zz                      ! ! intent(out)
   use mem_leaf           , only : albedo                  & ! intent(out)
                                 , drtcon                  & ! intent(out)
                                 , dthcon                  & ! intent(out)
                                 , isfcl                   & ! intent(out)
                                 , dtleaf                  & ! intent(out)
                                 , istar                   & ! intent(out)
                                 , igrndvap                & ! intent(out)
                                 , nslcon                  & ! intent(out)
                                 , isoilcol                & ! intent(out)
                                 , nvegpat                 & ! intent(out)
                                 , nvgcon                  & ! intent(out)
                                 , pctlcon                 & ! intent(out)
                                 , seatmp                  & ! intent(out)
                                 , slmstr                  & ! intent(out)
                                 , slz                     & ! intent(out)
                                 , stgoff                  & ! intent(out)
                                 , zrough                  & ! intent(out)
                                 , isoilbc                 & ! intent(out)
                                 , sldrain                 & ! intent(out)
                                 , ipercol                 & ! intent(out)
                                 , runoff_time             ! ! intent(out)
   use leaf_coms          , only : ubmin                   & ! intent(out)
                                 , ugbmin                  & ! intent(out)
                                 , ustmin                  & ! intent(out)
                                 , gamm                    & ! intent(out)
                                 , gamh                    & ! intent(out)
                                 , tprandtl                & ! intent(out)
                                 , ribmax                  & ! intent(out)
                                 , leaf_maxwhc             & ! intent(out)
                                 , min_patch_area          ! ! intent(out)
   use mem_oda            , only : frqoda                  & ! intent(out)
                                 , if_oda                  & ! intent(out)
                                 , oda_sfc_tel             & ! intent(out)
                                 , oda_sfc_til             & ! intent(out)
                                 , oda_sfcprefix           & ! intent(out)
                                 , oda_upa_tel             & ! intent(out)
                                 , oda_upa_til             & ! intent(out)
                                 , oda_upaprefix           & ! intent(out)
                                 , roda_hgt                & ! intent(out)
                                 , roda_sfc0               & ! intent(out)
                                 , roda_sfce               & ! intent(out)
                                 , roda_upa0               & ! intent(out)
                                 , roda_upae               & ! intent(out)
                                 , roda_zfact              & ! intent(out)
                                 , tnudoda                 & ! intent(out)
                                 , todabeg                 & ! intent(out)
                                 , todaend                 & ! intent(out)
                                 , wt_oda_grid             & ! intent(out)
                                 , wt_oda_pi               & ! intent(out)
                                 , wt_oda_rt               & ! intent(out)
                                 , wt_oda_th               & ! intent(out)
                                 , wt_oda_uv               ! ! intent(out)
   use mem_radiate        , only : ilwrtyp                 & ! intent(out)
                                 , iswrtyp                 & ! intent(out)
                                 , icumfdbk                & ! intent(out)
                                 , lonrad                  & ! intent(out)
                                 , radfrq                  ! ! intent(out)
   use mem_soil_moisture  , only : soil_moist              & ! intent(out)
                                 , soil_moist_fail         & ! intent(out)
                                 , usdata_in               & ! intent(out)
                                 , usmodel_in              ! ! intent(out)
   use mem_turb           , only : akmin                   & ! intent(out)
                                 , akmax                   & ! intent(out)
                                 , hgtmin                  & ! intent(out)
                                 , hgtmax                  & ! intent(out)
                                 , csx                     & ! intent(out)
                                 , csz                     & ! intent(out)
                                 , idiffk                  & ! intent(out)
                                 , ibotflx                 & ! intent(out)
                                 , ibruvais                & ! intent(out)
                                 , if_urban_canopy         & ! intent(out)
                                 , ihorgrad                & ! intent(out)
                                 , xkhkm                   & ! intent(out)
                                 , zkhkm                   ! ! intent(out)
   use turb_coms          , only : nna                     & ! intent(out)
                                 , nnb                     & ! intent(out)
                                 , nnc                     ! ! intent(out)
   use mem_varinit        , only : cond_hfile              & ! intent(out)
                                 , nud_cond                & ! intent(out)
                                 , nud_hfile               & ! intent(out)
                                 , nud_type                & ! intent(out)
                                 , nudlat                  & ! intent(out)
                                 , t_nudge_rc              & ! intent(out)
                                 , tcond_beg               & ! intent(out)
                                 , tcond_end               & ! intent(out)
                                 , tnudcent                & ! intent(out)
                                 , tnudlat                 & ! intent(out)
                                 , tnudtop                 & ! intent(out)
                                 , varfpfx                 & ! intent(out)
                                 , vwait1                  & ! intent(out)
                                 , vwaittot                & ! intent(out)
                                 , wt_nudge_grid           & ! intent(out)
                                 , wt_nudge_pi             & ! intent(out)
                                 , wt_nudge_rt             & ! intent(out)
                                 , wt_nudge_th             & ! intent(out)
                                 , wt_nudge_uv             & ! intent(out)
                                 , wt_nudge_co2            & ! intent(out)
                                 , wt_nudgec_grid          & ! intent(out)
                                 , znudtop                 ! ! intent(out)
   use micphys            , only : aparm                   & ! intent(out)
                                 , coltabfn                & ! intent(out)
                                 , cparm                   & ! intent(out)
                                 , gnu                     & ! intent(out)
                                 , gparm                   & ! intent(out)
                                 , hparm                   & ! intent(out)
                                 , iaggr                   & ! intent(out)
                                 , icloud                  & ! intent(out)
                                 , igraup                  & ! intent(out)
                                 , ihail                   & ! intent(out)
                                 , ipris                   & ! intent(out)
                                 , irain                   & ! intent(out)
                                 , isnow                   & ! intent(out)
                                 , mkcoltab                & ! intent(out)
                                 , pparm                   & ! intent(out)
                                 , rparm                   & ! intent(out)
                                 , sparm                   ! ! intent(out)
   use node_mod           , only : load_bal                ! ! intent(out)
   use ref_sounding       , only : hs                      & ! intent(out)
                                 , ipsflg                  & ! intent(out)
                                 , irtsflg                 & ! intent(out)
                                 , itsflg                  & ! intent(out)
                                 , iusflg                  & ! intent(out)
                                 , ps                      & ! intent(out)
                                 , rts                     & ! intent(out)
                                 , ts                      & ! intent(out)
                                 , us                      & ! intent(out)
                                 , vs                      & ! intent(out)
                                 , co2s                    ! ! intent(out)
   use catt_start         , only : catt                    ! ! intent(out)
   use emission_source_map, only : firemapfn               & ! intent(out)
                                 , tracersfn               & ! intent(out)
                                 , plumerise               ! ! intent(out)
   use plume_utils        , only : prfrq                   ! ! intent(out)
   use mem_scalar         , only : recycle_tracers         ! ! intent(out)
   use teb_spm_start      , only : teb_spm                 ! ! intent(out)
   use mem_emiss          , only : ichemi                  & ! intent(out)
                                 , ichemi_in               & ! intent(out)
                                 , chemdata_in             & ! intent(out)
                                 , isource                 & ! intent(out)
                                 , weekdayin               & ! intent(out)
                                 , efsat                   & ! intent(out)
                                 , efsun                   & ! intent(out)
                                 , eindno                  & ! intent(out)
                                 , eindno2                 & ! intent(out)
                                 , eindpm                  & ! intent(out)
                                 , eindco                  & ! intent(out)
                                 , eindso2                 & ! intent(out)
                                 , eindvoc                 & ! intent(out)
                                 , eveino                  & ! intent(out)
                                 , eveino2                 & ! intent(out)
                                 , eveipm                  & ! intent(out)
                                 , eveico                  & ! intent(out)
                                 , eveiso2                 & ! intent(out)
                                 , eveivoc                 ! ! intent(out)
   use teb_vars_const     , only : rushh1                  & ! intent(out)
                                 , rushh2                  & ! intent(out)
                                 , daylight                & ! intent(out)
                                 , iteb                    & ! intent(out)
                                 , tminbld                 & ! intent(out)
                                 , nteb                    & ! intent(out)
                                 , hc_roof                 & ! intent(out)
                                 , tc_roof                 & ! intent(out)
                                 , d_roof                  & ! intent(out)
                                 , hc_road                 & ! intent(out)
                                 , d_road                  & ! intent(out)
                                 , tc_road                 & ! intent(out)
                                 , d_wall                  & ! intent(out)
                                 , tc_wall                 & ! intent(out)
                                 , hc_wall                 & ! intent(out)
                                 , nurbtype                & ! intent(out)
                                 , ileafcod                & ! intent(out)
                                 , z0_town                 & ! intent(out)
                                 , bld                     & ! intent(out)
                                 , bld_height              & ! intent(out)
                                 , bld_hl_ratio            & ! intent(out)
                                 , aroof                   & ! intent(out)
                                 , eroof                   & ! intent(out)
                                 , aroad                   & ! intent(out)
                                 , eroad                   & ! intent(out)
                                 , awall                   & ! intent(out)
                                 , ewall                   & ! intent(out)
                                 , htraf                   & ! intent(out)
                                 , hindu                   & ! intent(out)
                                 , pletraf                 & ! intent(out)
                                 , pleindu                 ! ! intent(out)
   use mem_mass           , only : iexev                   & ! intent(out)
                                 , imassflx                ! ! intent(out)
   use grid_dims          , only : maxsteb                 & ! intent(out)
                                 , maxubtp                 ! ! intent(out)
   use mem_basic          , only : co2_on                  & ! intent(out)
                                 , co2con                  & ! intent(out)
                                 , ico2                    ! ! intent(out)
   use domain_decomp      , only : domain_fname            ! ! intent(out)
   use therm_lib          , only : vapour_on               & ! intent(out)
                                 , cloud_on                & ! intent(out)
                                 , bulk_on                 & ! intent(out)
                                 , level                   ! ! intent(out)
   use therm_lib8         , only : vapour_on8 => vapour_on & ! intent(out)
                                 , cloud_on8  => cloud_on  & ! intent(out)
                                 , bulk_on8   => bulk_on   & ! intent(out)
                                 , level8     => level     ! ! intent(out)
   use rconstants         , only : vonk                    ! ! intent(in)
   use mem_mnt_advec      , only : iadvec                  ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   character(len=*) , intent(in)   :: filename           ! file name with namelists
   !----- Local variables. ----------------------------------------------------------------!
   character(len=10)               :: c0                 ! scratch
   character(len=2) , dimension(1) :: caux               ! Aux. variable
   integer                         :: i                  ! loop count
   integer                         :: iunit              ! io unit number
   logical                         :: op                 ! io unit number opened or not
   logical                         :: ex                 ! namelist file exists?
   integer                         :: err                ! return code on iostat
   !----- Local constants. ----------------------------------------------------------------!
   integer          , parameter    :: firstunit=20       ! lowest io unit number avail.
   integer          , parameter    :: lastunit=99        ! highest io unit number avail.
   character(len=*) , parameter    :: h="**(read_nl)**"  ! program unit name
   !---------------------------------------------------------------------------------------!

   !----- List of namelists. --------------------------------------------------------------!
   namelist /MODEL_GRIDS/         expnme,runtype,load_bal,imontha,idatea,iyeara,itimea     &
                                 ,imonthz,idatez,iyearz,itimez,ngrids,nnxp,nnyp,nnzp,nzg   &
                                 ,nzs,nxtnest,domain_fname,if_adap,ihtran,deltax,deltay    &
                                 ,deltaz,dzrat,dzmax,zz,dtlong,nacoust,ideltat,nstratx     &
                                 ,nstraty,nndtrat,nestz1,nstratz1,nestz2,nstratz2,polelat  &
                                 ,polelon,centlat,centlon,ninest,njnest,nknest,nnsttop     &
                                 ,nnstbot,gridu,gridv

   namelist /CATT_INFO/           catt,firemapfn, recycle_tracers,plumerise, prfrq

   namelist /TEB_SPM_INFO/        teb_spm,fusfiles,ifusflg,ifusfn,ichemi,ichemi_in         &
                                 ,chemdata_in,isource,weekdayin,rushh1,rushh2,daylight     &
                                 ,efsat,efsun,eindno,eindno2,eindpm,eindco,eindso2,eindvoc &
                                 ,eveino,eveino2,eveipm,eveico,eveiso2,eveivoc,iteb        &
                                 ,tminbld,nteb,hc_roof,tc_roof,d_roof,hc_road,tc_road      &
                                 ,d_road,hc_wall,tc_wall,d_wall,nurbtype,ileafcod,z0_town  &
                                 ,bld,bld_height,bld_hl_ratio,aroof,eroof,aroad,eroad      &
                                 ,awall,ewall,htraf,hindu,pletraf,pleindu

   namelist /MODEL_FILE_INFO/     initial,nud_type,varfpfx,vwait1,vwaittot,nud_hfile       &
                                 ,nudlat,tnudlat,tnudcent,tnudtop,znudtop,wt_nudge_grid    &
                                 ,wt_nudge_uv,wt_nudge_th,wt_nudge_pi,wt_nudge_rt          &
                                 ,wt_nudge_co2,nud_cond,cond_hfile,tcond_beg,tcond_end     &
                                 ,t_nudge_rc,wt_nudgec_grid,if_oda,oda_upaprefix           &
                                 ,oda_sfcprefix,frqoda,todabeg,todaend,tnudoda,wt_oda_grid &
                                 ,wt_oda_uv,wt_oda_th,wt_oda_pi,wt_oda_rt,roda_sfce        &
                                 ,roda_sfc0,roda_upae,roda_upa0,roda_hgt,roda_zfact        &
                                 ,oda_sfc_til,oda_sfc_tel,oda_upa_til,oda_upa_tel,if_cuinv &
                                 ,cu_prefix,tnudcu,wt_cu_grid,tcu_beg,tcu_end,cu_tel       &
                                 ,cu_til,imonthh,idateh,iyearh,itimeh,hfilin,ipastin       &
                                 ,pastfn,ioutput,hfilout,afilout,iclobber,ihistdel,frqhis  &
                                 ,frqanl,frqlite,xlite,ylite,zlite,nlite_vars,lite_vars    &
                                 ,avgtim,frqmean,frqboth,kwrite,frqprt,initfld,topfiles    &
                                 ,sfcfiles,sstfpfx,ndvifpfx,itoptflg,isstflg,ivegtflg      &
                                 ,isoilflg,ndviflg,nofilflg,iupdndvi,iupdsst,iuselai       &
                                 ,itoptfn,isstfn,ivegtfn,isoilfn,ndvifn,itopsflg,toptenh   &
                                 ,toptwvl,iz0flg,z0max,z0fact,mkcoltab,coltabfn

   namelist /CUPARM_OPTIONS/      nnqparm,nclouds,ndeepest,nshallowest,wcldbs,confrq       &
                                 ,cptime,iupmethod,radius,depth_min,cap_maxs,cld2prec      &
                                 ,zkbmax,zcutdown,z_detr,max_heat,closure_type,maxens_lsf  &
                                 ,maxens_eff,maxens_cap

   namelist /MODEL_OPTIONS/       naddsc,icorflg,iadvec,iexev,imassflx,ibnd,jbnd,cphas     &
                                 ,lsflg,nfpt,distim,iswrtyp,ilwrtyp,icumfdbk,radfrq,lonrad &
                                 ,npatch,nvegpat,min_patch_area,isfcl,dtleaf,istar         &
                                 ,igrndvap,ubmin,ugbmin,ustmin,gamm,gamh,tprandtl,ribmax   &
                                 ,leaf_maxwhc,ico2,co2con,nvgcon,pctlcon,nslcon,isoilcol   &
                                 ,drtcon,zrough,albedo,seatmp,dthcon,soil_moist            &
                                 ,soil_moist_fail,usdata_in,usmodel_in,slz,slmstr,stgoff   &
                                 ,isoilbc,sldrain,ipercol,runoff_time,if_urban_canopy      &
                                 ,idiffk,ibruvais,ibotflx,ihorgrad,csx,csz,xkhkm,zkhkm     &
                                 ,nna,nnb,nnc,akmin,akmax,hgtmin,hgtmax,level,icloud,irain &
                                 ,ipris,isnow,iaggr,igraup,ihail,cparm,rparm,pparm,sparm   &
                                 ,aparm,gparm,hparm,gnu

   namelist /MODEL_SOUND/         ipsflg,itsflg,irtsflg,iusflg,hs,ps,ts,rts,us,vs,co2s

   namelist /MODEL_PRINT/         nplt,iplfld,ixsctn,isbval

   namelist /ISAN_CONTROL/        iszstage,ivrstage,isan_inc,guess1st,i1st_flg,iupa_flg    &
                                 ,isfc_flg,iapr,iarawi,iasrfce,varpfx,ioflgisz,ioflgvar

   namelist /ISAN_ISENTROPIC/     nisn,levth,nigrids,topsigz,hybbot,hybtop,sfcinf,sigzwt   &
                                 ,nfeedvar,maxsta,maxsfc,notsta,notid,iobswin,stasep       &
                                 ,igridfl,gridwt,gobsep,gobrad,wvlnth,swvlnth,respon
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Defaults (just for arrays); they should be revised to be non-sense values so the  !
   ! user to set up the correct values.                                                    !
   !---------------------------------------------------------------------------------------!
   ihoura                 = 0
   nnxp                   = 0
   nnyp                   = 0
   nnzp                   = 0
   nxtnest                = 0
   zz                     = 0.0
   nstratx                = 0
   nstraty                = 0
   nndtrat                = 0
   nstratz1               = 0
   nstratz2               = 0
   centlat                = 0.0
   centlon                = 0.0
   ninest                 = 0
   njnest                 = 0
   nknest                 = 0
   nnsttop                = 0
   nnstbot                = 0
   gridu                  = 0.0
   gridv                  = 0.0
   wt_nudge_grid          = 0.0
   wt_nudgec_grid         = 0.0
   wt_oda_grid            = 0.0
   roda_sfce              = 0.0
   roda_sfc0              = 0.0
   roda_upae              = 0.0
   roda_upa0              = 0.0
   roda_hgt               = 0.0
   roda_zfact             = 0.0
   wt_cu_grid             = 0.0
   lite_vars              = " "
   itoptflg               = 0
   isstflg                = 0
   ivegtflg               = 0
   isoilflg               = 0
   ndviflg                = 0
   nofilflg               = 0
   itoptfn                = " "
   isstfn                 = " "
   ivegtfn                = " "
   isoilfn                = " "
   ndvifn                 = " "
   itopsflg               = 0
   toptenh                = 0.0
   toptwvl                = 0.0
   co2con                 = 0.0
   ico2                   = 0
   iz0flg                 = 0
   z0max                  = 0.0
   gnu                    = 0.0
   iplfld                 = " "
   ixsctn                 = 0
   isbval                 = 0
   us                     = 0.0
   vs                     = 0.0
   ts                     = 0.0
   ps                     = 0.0
   hs                     = 0.0
   rts                    = 0.0
   co2s                   = 0.0
   slz                    = 0.0
   slmstr                 = 0.0
   stgoff                 = 0.0
   idiffk                 = 0
   ibotflx                = 0
   ibruvais               = 0
   csx                    = 0.0
   csz                    = 0.0
   xkhkm                  = 0.0
   zkhkm                  = 0.0
   nna                    = 0.0
   nnb                    = 0.0
   nnc                    = 0.0
   akmin                  = 0.0
   akmax                  = 0.0
   hgtmin                 = 0.0
   hgtmax                 = 0.0
   levth                  = 0
   notid                  = " "
   gridwt                 = 0.0
   wvlnth                 = 0.0
   swvlnth                = 0.0
   respon                 = 0.0
   domain_fname           = ''
   prfrq                  = 3600. 
   recycle_tracers        = 0


   nnqparm      = 0
   ndeepest     = 0
   nshallowest  = 0
   confrq       = 0.
   cptime       = 0.
   radius       = 0.
   depth_min    = 0.
   cap_maxs     = 0.
   cld2prec     = 0.
   zkbmax       = 0.
   zcutdown     = 0.
   z_detr       = 0.
   max_heat     = 0.
   closure_type = ""
   maxens_lsf   = 0
   maxens_eff   = 0
   maxens_cap   = 0


   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !     TEB.  Initialise with default values as they may not be in the namelist.          !
   !---------------------------------------------------------------------------------------!
   fusfiles    = ''
   ifusflg     = 0
   ichemi      = 0      ! Photochemical module activation - 1=on, 0=off
   ichemi_in   = 0      ! Use initial values from previous run (1=yes,0=no)
   chemdata_in = ''
   isource     = 0      ! Emission module activation - 1=on, 0=off
   weekdayin   = 'SUN'  ! Initial weeakday of the simulation
   rushh1      = 7.81   ! Morning Rush Hour (Local Time in Hours)
   rushh2      = 16.0   ! Afternoon/Evening Rush Hour (Local Time)
   daylight    = 0.     ! Daylight saving time (horario de verao)
   !---------------------------------------------------------------------------------------!
   !      Emission factor (fraction of weekdays) for Saturdays and Sundays.  They are used !
   ! in the emission module and TEB. - EDF                                                 !
   !---------------------------------------------------------------------------------------!
   efsat      = 0.8
   efsun      = 0.5
   !----- Industrial emissions (kg/s/m2). -------------------------------------------------!
   eindno     = 2.6636227e-10
   eindno2    = 2.9595805e-11
   eindpm     = 4.3421278e-10
   eindco     = 8.1599860e-10
   eindso2    = 3.6149164e-10
   eindvoc    = 2.5367833e-10 
   !----- Veicular emissions (kg/day/m2). -------------------------------------------------!
   eveino     = 4.3196708e-04
   eveino2    = 6.8566209e-05
   eveipm     = 6.2648396e-06
   eveico     = 7.5433785e-03
   eveiso2    = 4.0730592e-05
   eveivoc    = 1.1892237e-03
   !----- Urban canopy parameterization using TEB (Masson, 2000). -------------------------!
   iteb       = 0     ! 1=on, 0=off
   tminbld    = 12.   ! Minimum internal building temperature (degrees Celsius)
   nteb       = 3     ! Number of roof,road and wall layers used in TEB, Max.3
   if (maxsteb == 3) then
      !------------------------------------------------------------------------------------!
      !     ROOF layers properties.                                                        !
      !------------------------------------------------------------------------------------!
      !----- Heat capacity. ---------------------------------------------------------------!
      hc_roof(1:3) = (/ 2110000., 280000., 290000. /)
      !----- Thermal conductivity. --------------------------------------------------------!
      tc_roof(1:3) = (/     0.41,    0.05,    0.03 /)
      !----- Depth. -----------------------------------------------------------------------!
      d_roof(1:3)  = (/     0.05,     0.4,    0.05 /)
      !------------------------------------------------------------------------------------!
      !     ROAD layers properties.                                                        !
      !------------------------------------------------------------------------------------!
      !----- Heat capacity. ---------------------------------------------------------------!
      ! heat capacity
      hc_road(1:3) = (/ 1240000., 1280000., 1280000. /)
      !----- Thermal conductivity. --------------------------------------------------------!
      tc_road      = 1.0103
      !----- Depth. -----------------------------------------------------------------------!
      d_road(1:3)  = (/     0.05,      0.1,      1.0 /)
      !------------------------------------------------------------------------------------!
      !     WALL layers properties.                                                        !
      !------------------------------------------------------------------------------------!
      !----- Heat capacity. ---------------------------------------------------------------!
      hc_wall = 1000000.
      !----- Thermal conductivity. --------------------------------------------------------!
      tc_wall = 0.81
      !----- Depth. -----------------------------------------------------------------------!
      d_wall(1:3)  = (/ 0.02, 0.125, 0.02 /)
   else
      !------------------------------------------------------------------------------------!
      !     ROOF layers properties.                                                        !
      !------------------------------------------------------------------------------------!
      hc_roof    = 0.    ! heat capacity
      tc_roof    = 0.    ! thermal conductivity
      d_roof     = 0.    ! depth
      !------------------------------------------------------------------------------------!
      !     ROAD layers properties.                                                        !
      !------------------------------------------------------------------------------------!
      hc_road    = 0.    ! heat capacity
      tc_road    = 0.    ! thermal conductivity 1.0103
      d_road     = 0.    ! depth
      !------------------------------------------------------------------------------------!
      !     WALL layers properties.                                                        !
      !------------------------------------------------------------------------------------!
      hc_wall    = 0.    ! heat capacity J/m3/K 10e6
      tc_wall    = 0.    ! thermal conductivity 0.81 W/m/K
      d_wall     = 0.    ! depth
   end if
   nurbtype   = 2         !Number of urban types (maximum of 3)
   if (maxubtp == 2) then
      !----- Leaf class code to identify each urban type. ---------------------------------!
      ileafcod(1:2)     = (/ 21, 19 /)
      !----- Urban type roughness length 5 and 1. -----------------------------------------!
      z0_town(1:2)      = (/ 3.0, 0.5 /)
      !----- Fraction occupied by buildings in the grid cell. -----------------------------!
      bld(1:2)          = (/ 0.5, 0.7 /)
      !----- Building Height. -------------------------------------------------------------!
      bld_height(1:2)   = (/ 50., 5.0 /)
      !----- Vertical/Horizontal rate 3 and 0.5. ------------------------------------------!
      bld_hl_ratio(1:2) = (/ 4.4, 2.4 /)
      !----- Roof albedo. -----------------------------------------------------------------!
      aroof = 0.15
      !----- Roof emissivity. -------------------------------------------------------------!
      eroof = 0.9
      !----- Road albedo. -----------------------------------------------------------------!
      aroad = 0.1
      !----- Road emissivity. -------------------------------------------------------------!
      eroad = 0.9
      !----- Wall albedo. -----------------------------------------------------------------!
      awall = 0.25
      !----- Wall emissivity. -------------------------------------------------------------!
      ewall = 0.85
      !----- Maximum value of sensible heat released by traffic (W/m2). -------------------!
      htraf(1:2)   = (/ 90.0, 60.0 /)
      !----- Maximum value of sensible heat released by industry (W/m2). ------------------!
      hindu(1:2)   = (/ 10.0, 14.0 /)
      !----- Maximum value of latent heat released by traffic (W/m2). ---------------------!
      pletraf(1:2) = (/ 10.0,  5.0 /)
      !----- Maximum value of latent heat released by industry (W/m2). --------------------!
      pleindu(1:2) = (/ 30.0, 50.0 /)
   else
      !----- Urban type properties. -------------------------------------------------------!
      z0_town      = 0.0   ! Urban type roughness length 5 e 1
      bld          = 0.0   ! Fraction occupied by buildings in the grid cell
      bld_height   = 0.0   ! Building Height
      bld_hl_ratio = 0.0   ! Vertical/Horizontal rate 3 e 0.5
      aroof        = 0.0   ! Roof albedo
      eroof        = 0.0   ! Roof emissivitiy
      aroad        = 0.0   ! Road albedo
      eroad        = 0.0   ! Road emissivity 90% masson
      awall        = 0.0   ! Wall albedo
      ewall        = 0.0   ! Wall emissivity
      htraf        = 0.0   ! Maximum value of sensible heat
      hindu        = 0.0   ! Maximum value of sensible heat
      pletraf      = 0.0   ! Maximum value of latent heat
      pleindu      = 0.0   ! Maximum value of latent heat
   end if
   !---------------------------------------------------------------------------------------!

   !----- select unused i/o unit. ---------------------------------------------------------!
   unitsearch: do iunit = firstunit, lastunit
      inquire(iunit,opened=op)
      if (.not. op) exit unitsearch
   end do unitsearch
   if (iunit > lastUnit) then
      call abort_run('All i/o units are already in use.','ReadNamelist','rdint.f90')
   end if
   !---------------------------------------------------------------------------------------!


   !----- If namelist file exists, open, read each section and close. ---------------------!
   inquire(file=trim(filename), exist=ex)
   if (.not. ex) then
      call abort_run('Namelist file '//trim(fileName)//' doesn''t exist.'                  &
                    ,'read_nl','rdint.f90')
   end if

   open(unit=iunit,file=trim(fileName), status='old', action='read',iostat=err)
   if (err /= 0) then
      write (unit=*,fmt='(a,1x,a)')   ' Namelist file: ',trim(filename)
      write (unit=*,fmt='(a,1x,i10)') ' Iostat       : ',err
      call abort_run(' Error reading the namelist...','read_nl','rdint.f90')
   end if
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !      Read in MODEL_GRIDS namelist.                                                    !
   !---------------------------------------------------------------------------------------!
   read (unit=iunit, iostat=err, nml=MODEL_GRIDS)
   if (err /= 0) then
      write (unit=*,fmt='(a)')        '--------------------------------------------------'
      write (unit=*,fmt='(a)')        ' ERROR reading the namelist! '
      write (unit=*,fmt='(a,1x,a)')   ' Namelist file   : ',trim(filename)
      write (unit=*,fmt='(a,1x,a)')   ' Namelist section: ','MODEL_GRIDS'
      write (unit=*,fmt='(a)')        ' Compare values read with file contents...'
      write (unit=*,fmt='(a)')        '--------------------------------------------------'
      write (unit=*,fmt='(a)')        ''
      write (unit=*,fmt=*) ' expnme       =',trim(expnme)
      write (unit=*,fmt=*) ' runtype      =',trim(runtype)
      write (unit=*,fmt=*) ' load_bal     =',load_bal
      write (unit=*,fmt=*) ' imontha      =',imontha
      write (unit=*,fmt=*) ' idatea       =',idatea
      write (unit=*,fmt=*) ' iyeara       =',iyeara
      write (unit=*,fmt=*) ' itimea       =',itimea
      write (unit=*,fmt=*) ' imonthz      =',imonthz
      write (unit=*,fmt=*) ' idatez       =',idatez
      write (unit=*,fmt=*) ' iyearz       =',iyearz
      write (unit=*,fmt=*) ' itimez       =',itimez
      write (unit=*,fmt=*) ' ngrids       =',ngrids
      write (unit=*,fmt=*) ' nnxp         =',nnxp
      write (unit=*,fmt=*) ' nnyp         =',nnyp
      write (unit=*,fmt=*) ' nnzp         =',nnzp
      write (unit=*,fmt=*) ' nzg          =',nzg
      write (unit=*,fmt=*) ' nzs          =',nzs
      write (unit=*,fmt=*) ' nxtnest      =',nxtnest
      write (unit=*,fmt=*) ' DOMAIN_FNAME =',DOMAIN_FNAME
      write (unit=*,fmt=*) ' if_adap      =',if_adap
      write (unit=*,fmt=*) ' ihtran       =',ihtran
      write (unit=*,fmt=*) ' deltax       =',deltax
      write (unit=*,fmt=*) ' deltay       =',deltay
      write (unit=*,fmt=*) ' deltaz       =',deltaz
      write (unit=*,fmt=*) ' dzrat        =',dzrat
      write (unit=*,fmt=*) ' dzmax        =',dzmax
      write (unit=*,fmt=*) ' zz           =',zz
      write (unit=*,fmt=*) ' dtlong       =',dtlong
      write (unit=*,fmt=*) ' nacoust      =',nacoust
      write (unit=*,fmt=*) ' ideltat      =',ideltat
      write (unit=*,fmt=*) ' nstratx      =',nstratx
      write (unit=*,fmt=*) ' nstraty      =',nstraty
      write (unit=*,fmt=*) ' nndtrat      =',nndtrat
      write (unit=*,fmt=*) ' nestz1       =',nestz1
      write (unit=*,fmt=*) ' nstratz1     =',nstratz1
      write (unit=*,fmt=*) ' nestz2       =',nestz2
      write (unit=*,fmt=*) ' nstratz2     =',nstratz2
      write (unit=*,fmt=*) ' polelat      =',polelat
      write (unit=*,fmt=*) ' polelon      =',polelon
      write (unit=*,fmt=*) ' centlat      =',centlat
      write (unit=*,fmt=*) ' centlon      =',centlon
      write (unit=*,fmt=*) ' ninest       =',ninest
      write (unit=*,fmt=*) ' njnest       =',njnest
      write (unit=*,fmt=*) ' nknest       =',nknest
      write (unit=*,fmt=*) ' nnsttop      =',nnsttop
      write (unit=*,fmt=*) ' nnstbot      =',nnstbot
      write (unit=*,fmt=*) ' gridu        =',gridu
      write (unit=*,fmt=*) ' gridv        =',gridv
      call abort_run('Error reading namelist...','read_nl','rdint.f90')
   end if
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !      Read in CATT_INFO namelist.                                                      !
   !---------------------------------------------------------------------------------------!
   read (unit=iunit, iostat=err, nml=CATT_INFO)
   if (err /= 0) then
      write (unit=*,fmt='(a)')        '--------------------------------------------------'
      write (unit=*,fmt='(a)')        ' ERROR reading the namelist! '
      write (unit=*,fmt='(a,1x,a)')   ' Namelist file   : ',trim(filename)
      write (unit=*,fmt='(a,1x,a)')   ' Namelist section: ','CATT_INFO'
      write (unit=*,fmt='(a)')        ' Compare values read with file contents...'
      write (unit=*,fmt='(a)')        '--------------------------------------------------'
      write (unit=*,fmt='(a)')        ''

      write (unit=*,fmt=*) ' catt            =', catt
      write (unit=*,fmt=*) ' firemapfn       =', firemapfn
      write (unit=*,fmt=*) ' recycle_tracers =', recycle_tracers
      write (unit=*,fmt=*) ' plumerise       =', plumerise
      write (unit=*,fmt=*) ' prfrq           =', prfrq
      call abort_run('Error reading namelist...','read_nl','rdint.f90')
   end if
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !      Read in TEB_SPM_INFO namelist.                                                   !
   !---------------------------------------------------------------------------------------!
   read (iunit, iostat=err, NML=TEB_SPM_INFO)
   if (err /= 0) then
      write (unit=*,fmt='(a)')        '--------------------------------------------------'
      write (unit=*,fmt='(a)')        ' ERROR reading the namelist! '
      write (unit=*,fmt='(a,1x,a)')   ' Namelist file   : ',trim(filename)
      write (unit=*,fmt='(a,1x,a)')   ' Namelist section: ','TEB_SPM_INFO'
      write (unit=*,fmt='(a)')        ' Compare values read with file contents...'
      write (unit=*,fmt='(a)')        '--------------------------------------------------'
      write (unit=*,fmt='(a)')        ''
      write (unit=*,fmt=*)  ' TEB_SPM      =',TEB_SPM
      write (unit=*,fmt=*)  ' ifusflg      =',ifusflg
      write (unit=*,fmt=*)  ' ifusfn       =',ifusfn
      write (unit=*,fmt=*)  ' fusfiles     =',trim(fusfiles)
      write (unit=*,fmt=*)  ' ICHEMI       =',ICHEMI
      write (unit=*,fmt=*)  ' ICHEMI_IN    =',ICHEMI_IN
      write (unit=*,fmt=*)  ' CHEMDATA_IN  =',CHEMDATA_IN
      write (unit=*,fmt=*)  ' ISOURCE      =',ISOURCE
      write (unit=*,fmt=*)  ' WEEKDAYIN    =',trim(WEEKDAYIN)
      write (unit=*,fmt=*)  ' RUSHH1       =',RUSHH1
      write (unit=*,fmt=*)  ' RUSHH2       =',RUSHH2
      write (unit=*,fmt=*)  ' DAYLIGHT     =',DAYLIGHT
      write (unit=*,fmt=*)  ' EFSAT        =',EFSAT
      write (unit=*,fmt=*)  ' EFSUN        =',EFSUN
      write (unit=*,fmt=*)  ' EINDNO       =',EINDNO
      write (unit=*,fmt=*)  ' EINDNO2      =',EINDNO2
      write (unit=*,fmt=*)  ' EINDPM       =',EINDPM
      write (unit=*,fmt=*)  ' EINDCO       =',EINDCO
      write (unit=*,fmt=*)  ' EINDSO2      =',EINDSO2
      write (unit=*,fmt=*)  ' EINDVOC      =',EINDVOC
      write (unit=*,fmt=*)  ' EVEINO       =',EVEINO
      write (unit=*,fmt=*)  ' EVEINO2      =',EVEINO2
      write (unit=*,fmt=*)  ' EVEIPM       =',EVEIPM
      write (unit=*,fmt=*)  ' EVEICO       =',EVEICO
      write (unit=*,fmt=*)  ' EVEISO2      =',EVEISO2
      write (unit=*,fmt=*)  ' EVEIVOC      =',EVEIVOC
      write (unit=*,fmt=*)  ' ITEB         =',ITEB
      write (unit=*,fmt=*)  ' TMINBLD      =',TMINBLD
      write (unit=*,fmt=*)  ' NTEB         =',NTEB
      write (unit=*,fmt=*)  ' HC_ROOF      =',HC_ROOF
      write (unit=*,fmt=*)  ' TC_ROOF      =',TC_ROOF
      write (unit=*,fmt=*)  ' D_ROOF       =',D_ROOF
      write (unit=*,fmt=*)  ' HC_ROAD      =',HC_ROAD
      write (unit=*,fmt=*)  ' TC_ROAD      =',TC_ROAD
      write (unit=*,fmt=*)  ' D_ROAD       =',D_ROAD
      write (unit=*,fmt=*)  ' HC_WALL      =',HC_WALL
      write (unit=*,fmt=*)  ' TC_WALL      =',TC_WALL
      write (unit=*,fmt=*)  ' D_WALL       =',D_WALL
      write (unit=*,fmt=*)  ' NURBTYPE     =',NURBTYPE
      write (unit=*,fmt=*)  ' ILEAFCOD     =',ILEAFCOD
      write (unit=*,fmt=*)  ' Z0_TOWN      =',Z0_TOWN
      write (unit=*,fmt=*)  ' BLD          =',BLD
      write (unit=*,fmt=*)  ' BLD_HEIGHT   =',BLD_HEIGHT
      write (unit=*,fmt=*)  ' BLD_HL_RATIO =',BLD_HL_RATIO
      write (unit=*,fmt=*)  ' AROOF        =',AROOF
      write (unit=*,fmt=*)  ' EROOF        =',EROOF
      write (unit=*,fmt=*)  ' AROAD        =',AROAD
      write (unit=*,fmt=*)  ' EROAD        =',EROAD
      write (unit=*,fmt=*)  ' AWALL        =',AWALL
      write (unit=*,fmt=*)  ' EWALL        =',EWALL
      write (unit=*,fmt=*)  ' HTRAF        =',HTRAF
      write (unit=*,fmt=*)  ' HINDU        =',HINDU
      write (unit=*,fmt=*)  ' PLETRAF      =',PLETRAF
      write (unit=*,fmt=*)  ' PLEINDU      =',PLEINDU
      call abort_run('Error reading namelist...','read_nl','rdint.f90')
   end if
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !      Read in MODEL_FILE_INFO namelist.                                                !
   !---------------------------------------------------------------------------------------!
   read (unit=iunit, iostat=err, nml=MODEL_FILE_INFO)
   if (err /= 0) then
      write (unit=*,fmt='(a)')        '--------------------------------------------------'
      write (unit=*,fmt='(a)')        ' ERROR reading the namelist! '
      write (unit=*,fmt='(a,1x,a)')   ' Namelist file   : ',trim(filename)
      write (unit=*,fmt='(a,1x,a)')   ' Namelist section: ','MODEL_FILE_INFO'
      write (unit=*,fmt='(a)')        ' Compare values read with file contents...'
      write (unit=*,fmt='(a)')        '--------------------------------------------------'
      write (unit=*,fmt='(a)')        ''
      write (unit=*,fmt=*) ' initial        =',initial
      write (unit=*,fmt=*) ' nud_type       =',nud_type
      write (unit=*,fmt=*) ' varfpfx        =',trim(varfpfx)
      write (unit=*,fmt=*) ' vwait1         =',vwait1
      write (unit=*,fmt=*) ' vwaittot       =',vwaittot
      write (unit=*,fmt=*) ' nud_hfile      =',trim(nud_hfile)
      write (unit=*,fmt=*) ' nudlat         =',nudlat
      write (unit=*,fmt=*) ' tnudlat        =',tnudlat
      write (unit=*,fmt=*) ' tnudcent       =',tnudcent
      write (unit=*,fmt=*) ' tnudtop        =',tnudtop
      write (unit=*,fmt=*) ' znudtop        =',znudtop
      write (unit=*,fmt=*) ' wt_nudge_grid  =',wt_nudge_grid
      write (unit=*,fmt=*) ' wt_nudge_uv    =',wt_nudge_uv
      write (unit=*,fmt=*) ' wt_nudge_th    =',wt_nudge_th
      write (unit=*,fmt=*) ' wt_nudge_pi    =',wt_nudge_pi
      write (unit=*,fmt=*) ' wt_nudge_rt    =',wt_nudge_rt
      write (unit=*,fmt=*) ' wt_nudge_co2   =',wt_nudge_co2
      write (unit=*,fmt=*) ' nud_cond       =',nud_cond
      write (unit=*,fmt=*) ' cond_hfile     =',trim(cond_hfile)
      write (unit=*,fmt=*) ' tcond_beg      =',tcond_beg
      write (unit=*,fmt=*) ' tcond_end      =',tcond_end
      write (unit=*,fmt=*) ' t_nudge_rc     =',t_nudge_rc
      write (unit=*,fmt=*) ' wt_nudgec_grid =',wt_nudgec_grid
      write (unit=*,fmt=*) ' if_oda         =',if_oda
      write (unit=*,fmt=*) ' oda_upaprefix  =',trim(oda_upaprefix)
      write (unit=*,fmt=*) ' oda_sfcprefix  =',trim(oda_sfcprefix)
      write (unit=*,fmt=*) ' frqoda         =',frqoda
      write (unit=*,fmt=*) ' todabeg        =',todabeg
      write (unit=*,fmt=*) ' todaend        =',todaend
      write (unit=*,fmt=*) ' tnudoda        =',tnudoda
      write (unit=*,fmt=*) ' wt_oda_grid    =',wt_oda_grid
      write (unit=*,fmt=*) ' wt_oda_uv      =',wt_oda_uv
      write (unit=*,fmt=*) ' wt_oda_th      =',wt_oda_th
      write (unit=*,fmt=*) ' wt_oda_pi      =',wt_oda_pi
      write (unit=*,fmt=*) ' wt_oda_rt      =',wt_oda_rt
      write (unit=*,fmt=*) ' roda_sfce      =',roda_sfce
      write (unit=*,fmt=*) ' roda_sfc0      =',roda_sfc0
      write (unit=*,fmt=*) ' roda_upae      =',roda_upae
      write (unit=*,fmt=*) ' roda_upa0      =',roda_upa0
      write (unit=*,fmt=*) ' roda_hgt       =',roda_hgt
      write (unit=*,fmt=*) ' roda_zfact     =',roda_zfact
      write (unit=*,fmt=*) ' oda_sfc_til    =',oda_sfc_til
      write (unit=*,fmt=*) ' oda_sfc_tel    =',oda_sfc_tel
      write (unit=*,fmt=*) ' oda_upa_til    =',oda_upa_til
      write (unit=*,fmt=*) ' oda_upa_tel    =',oda_upa_tel
      write (unit=*,fmt=*) ' if_cuinv       =',if_cuinv
      write (unit=*,fmt=*) ' cu_prefix      =',trim(cu_prefix)
      write (unit=*,fmt=*) ' tnudcu         =',tnudcu
      write (unit=*,fmt=*) ' wt_cu_grid     =',wt_cu_grid
      write (unit=*,fmt=*) ' tcu_beg        =',tcu_beg
      write (unit=*,fmt=*) ' tcu_end        =',tcu_end
      write (unit=*,fmt=*) ' cu_tel         =',cu_tel
      write (unit=*,fmt=*) ' cu_til         =',cu_til
      write (unit=*,fmt=*) ' imonth         =',imonthh
      write (unit=*,fmt=*) ' idateh         =',idateh
      write (unit=*,fmt=*) ' iyearh         =',iyearh
      write (unit=*,fmt=*) ' itimeh         =',itimeh
      write (unit=*,fmt=*) ' hfilin         =',trim(hfilin)
      write (unit=*,fmt=*) ' ipastin        =',ipastin
      write (unit=*,fmt=*) ' pastfn         =',trim(pastfn)
      write (unit=*,fmt=*) ' ioutput        =',ioutput
      write (unit=*,fmt=*) ' hfilout        =',trim(hfilout)
      write (unit=*,fmt=*) ' afilout        =',trim(afilout)
      write (unit=*,fmt=*) ' iclobber       =',iclobber
      write (unit=*,fmt=*) ' ihistdel       =',ihistdel
      write (unit=*,fmt=*) ' frqhis         =',frqhis
      write (unit=*,fmt=*) ' frqanl         =',frqanl
      write (unit=*,fmt=*) ' frqlite        =',frqlite
      write (unit=*,fmt=*) ' xlite          =',xlite
      write (unit=*,fmt=*) ' ylite          =',ylite
      write (unit=*,fmt=*) ' zlite          =',zlite
      write (unit=*,fmt=*) ' nlite_vars     =',nlite_vars
      write (unit=*,fmt=*) ' lite_vars      =',(trim(lite_vars(i))//';',i=1,size(lite_vars))
      write (unit=*,fmt=*) ' avgtim         =',avgtim
      write (unit=*,fmt=*) ' frqmean        =',frqmean
      write (unit=*,fmt=*) ' frqboth        =',frqboth
      write (unit=*,fmt=*) ' kwrite         =',kwrite
      write (unit=*,fmt=*) ' frqprt         =',frqprt
      write (unit=*,fmt=*) ' initfld        =',initfld
      write (unit=*,fmt=*) ' topfiles       =',trim(topfiles)
      write (unit=*,fmt=*) ' sfcfiles       =',trim(sfcfiles)
      write (unit=*,fmt=*) ' sstfpfx        =',trim(sstfpfx)
      write (unit=*,fmt=*) ' ndvifpfx       =',trim(ndvifpfx)
      write (unit=*,fmt=*) ' itoptflg       =',itoptflg
      write (unit=*,fmt=*) ' isstflg        =',isstflg
      write (unit=*,fmt=*) ' ivegtflg       =',ivegtflg
      write (unit=*,fmt=*) ' isoilflg       =',isoilflg
      write (unit=*,fmt=*) ' ndviflg        =',ndviflg
      write (unit=*,fmt=*) ' nofilflg       =',nofilflg
      write (unit=*,fmt=*) ' iupdndvi       =',iupdndvi
      write (unit=*,fmt=*) ' iupdsst        =',iupdsst
      write (unit=*,fmt=*) ' itoptfn        =',(trim(itoptfn(i))//';', i=1,size(itoptfn))
      write (unit=*,fmt=*) ' isstfn         =',(trim(isstfn(i))//';' , i=1,size(isstfn) )
      write (unit=*,fmt=*) ' ivegtfn        =',(trim(ivegtfn(i))//';', i=1,size(ivegtfn))
      write (unit=*,fmt=*) ' isoilfn        =',(trim(isoilfn(i))//';', i=1,size(isoilfn))
      write (unit=*,fmt=*) ' ndvifn         =',(trim(ndvifn(i))//';' , i=1,size(ndvifn) )
      write (unit=*,fmt=*) ' itopsflg       =',itopsflg
      write (unit=*,fmt=*) ' toptenh        =',toptenh
      write (unit=*,fmt=*) ' toptwvl        =',toptwvl
      write (unit=*,fmt=*) ' iz0flg         =',iz0flg
      write (unit=*,fmt=*) ' z0max          =',z0max
      write (unit=*,fmt=*) ' z0fact         =',z0fact
      write (unit=*,fmt=*) ' mkcoltab       =',mkcoltab
      write (unit=*,fmt=*) ' coltabfn       =',trim(coltabfn)
      call abort_run('Error reading namelist...','read_nl','rdint.f90')
   end if
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !      Read in CUPARM_OPTIONS namelist.                                                 !
   !---------------------------------------------------------------------------------------!
   read (unit=iunit, iostat=err, nml=CUPARM_OPTIONS)
   if (err /= 0) then
      write (unit=*,fmt='(a)')        '--------------------------------------------------'
      write (unit=*,fmt='(a)')        ' ERROR reading the namelist! '
      write (unit=*,fmt='(a,1x,a)')   ' Namelist file   : ',trim(filename)
      write (unit=*,fmt='(a,1x,a)')   ' Namelist section: ','CUPARM_OPTIONS'
      write (unit=*,fmt='(a)')        ' Compare values read with file contents...'
      write (unit=*,fmt='(a)')        '--------------------------------------------------'
      write (unit=*,fmt='(a)')        ''
      write (unit=*,fmt=*) ' nnqparm      =',nnqparm
      write (unit=*,fmt=*) ' nclouds      =',nclouds
      write (unit=*,fmt=*) ' ndeepest     =',ndeepest
      write (unit=*,fmt=*) ' nshallowest  =',nshallowest
      write (unit=*,fmt=*) ' wcldbs       =',wcldbs
      write (unit=*,fmt=*) ' confrq       =',confrq
      write (unit=*,fmt=*) ' cptime       =',cptime
      write (unit=*,fmt=*) ' iupmethod    =',iupmethod
      write (unit=*,fmt=*) ' radius       =',radius
      write (unit=*,fmt=*) ' depth_min    =',depth_min
      write (unit=*,fmt=*) ' cap_maxs     =',depth_min
      write (unit=*,fmt=*) ' cld2prec     =',cld2prec
      write (unit=*,fmt=*) ' zkbmax       =',zkbmax
      write (unit=*,fmt=*) ' zcutdown     =',zcutdown
      write (unit=*,fmt=*) ' z_detr       =',z_detr
      write (unit=*,fmt=*) ' max_heat     =',max_heat
      write (unit=*,fmt=*) ' closure_type =',closure_type
      write (unit=*,fmt=*) ' maxens_lsf   =',maxens_lsf
      write (unit=*,fmt=*) ' maxens_eff   =',maxens_eff
      write (unit=*,fmt=*) ' maxens_cap   =',maxens_cap
      call abort_run('Error reading namelist...','read_nl','rdint.f90')
   end if
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !      Read in MODEL_OPTIONS namelist.                                                  !
   !---------------------------------------------------------------------------------------!
   read (unit=iunit, iostat=err, nml=MODEL_OPTIONS)
   if (err /= 0) then
      write (unit=*,fmt='(a)')        '--------------------------------------------------'
      write (unit=*,fmt='(a)')        ' ERROR reading the namelist! '
      write (unit=*,fmt='(a,1x,a)')   ' Namelist file   : ',trim(filename)
      write (unit=*,fmt='(a,1x,a)')   ' Namelist section: ','MODEL_OPTIONS'
      write (unit=*,fmt='(a)')        ' Compare values read with file contents...'
      write (unit=*,fmt='(a)')        '--------------------------------------------------'
      write (unit=*,fmt='(a)')        ''
      write (unit=*,fmt=*) ' naddsc          =',naddsc
      write (unit=*,fmt=*) ' icorflg         =',icorflg
      write (unit=*,fmt=*) ' iadvec          =',iadvec
      write (unit=*,fmt=*) ' iexev           =',iexev
      write (unit=*,fmt=*) ' imassflx        =',imassflx
      write (unit=*,fmt=*) ' ibnd            =',ibnd
      write (unit=*,fmt=*) ' jbnd            =',jbnd
      write (unit=*,fmt=*) ' cphas           =',cphas
      write (unit=*,fmt=*) ' lsflg           =',lsflg
      write (unit=*,fmt=*) ' nfpt            =',nfpt
      write (unit=*,fmt=*) ' distim          =',distim
      write (unit=*,fmt=*) ' iswrtyp         =',iswrtyp
      write (unit=*,fmt=*) ' ilwrtyp         =',ilwrtyp
      write (unit=*,fmt=*) ' radfrq          =',radfrq
      write (unit=*,fmt=*) ' lonrad          =',lonrad
      write (unit=*,fmt=*) ' npatch          =',npatch
      write (unit=*,fmt=*) ' nvegpat         =',nvegpat
      write (unit=*,fmt=*) ' min_patch_area  =',min_patch_area
      write (unit=*,fmt=*) ' isfcl           =',isfcl
      write (unit=*,fmt=*) ' dtleaf          =',dtleaf
      write (unit=*,fmt=*) ' istar           =',istar
      write (unit=*,fmt=*) ' igrndvap        =',igrndvap
      write (unit=*,fmt=*) ' ubmin           =',ubmin
      write (unit=*,fmt=*) ' ugbmin          =',ugbmin
      write (unit=*,fmt=*) ' ustmin          =',ustmin
      write (unit=*,fmt=*) ' gamm            =',gamm
      write (unit=*,fmt=*) ' gamh            =',gamh
      write (unit=*,fmt=*) ' tprandtl        =',tprandtl
      write (unit=*,fmt=*) ' ribmax          =',ribmax
      write (unit=*,fmt=*) ' leaf_maxwhc     =',leaf_maxwhc
      write (unit=*,fmt=*) ' ico2            =',ico2
      write (unit=*,fmt=*) ' co2con          =',co2con
      write (unit=*,fmt=*) ' nvgcon          =',nvgcon
      write (unit=*,fmt=*) ' pctlcon         =',pctlcon
      write (unit=*,fmt=*) ' nslcon          =',nslcon
      write (unit=*,fmt=*) ' isoilcol        =',isoilcol
      write (unit=*,fmt=*) ' drtcon          =',drtcon
      write (unit=*,fmt=*) ' zrough          =',zrough
      write (unit=*,fmt=*) ' albedo          =',albedo
      write (unit=*,fmt=*) ' seatmp          =',seatmp
      write (unit=*,fmt=*) ' dthcon          =',dthcon
      write (unit=*,fmt=*) ' soil_moist      =',soil_moist
      write (unit=*,fmt=*) ' soil_moist_fail =',soil_moist_fail
      write (unit=*,fmt=*) ' usdata_in       =',trim(usdata_in)
      write (unit=*,fmt=*) ' usmodel_in      =',trim(usmodel_in)
      write (unit=*,fmt=*) ' slz             =',slz
      write (unit=*,fmt=*) ' slmstr          =',slmstr
      write (unit=*,fmt=*) ' stgoff          =',stgoff
      write (unit=*,fmt=*) ' isoilbc         =',isoilbc
      write (unit=*,fmt=*) ' sldrain         =',sldrain
      write (unit=*,fmt=*) ' ipercol         =',ipercol
      write (unit=*,fmt=*) ' runoff_time     =',runoff_time
      write (unit=*,fmt=*) ' if_urban_canopy =',if_urban_canopy
      write (unit=*,fmt=*) ' idiffk          =',idiffk
      write (unit=*,fmt=*) ' ibruvais        =',ibruvais
      write (unit=*,fmt=*) ' ibotflx         =',ibotflx
      write (unit=*,fmt=*) ' ihorgrad        =',ihorgrad
      write (unit=*,fmt=*) ' csx             =',csx
      write (unit=*,fmt=*) ' csz             =',csz
      write (unit=*,fmt=*) ' xkhkm           =',xkhkm
      write (unit=*,fmt=*) ' zkhkm           =',zkhkm
      write (unit=*,fmt=*) ' nna             =',nna
      write (unit=*,fmt=*) ' nnb             =',nnb
      write (unit=*,fmt=*) ' nnc             =',nnc
      write (unit=*,fmt=*) ' akmin           =',akmin
      write (unit=*,fmt=*) ' akmax           =',akmax
      write (unit=*,fmt=*) ' hgtmin          =',hgtmin
      write (unit=*,fmt=*) ' hgtmax          =',hgtmax
      write (unit=*,fmt=*) ' level           =',level
      write (unit=*,fmt=*) ' icloud          =',icloud
      write (unit=*,fmt=*) ' irain           =',irain
      write (unit=*,fmt=*) ' ipris           =',ipris
      write (unit=*,fmt=*) ' isnow           =',isnow
      write (unit=*,fmt=*) ' iaggr           =',iaggr
      write (unit=*,fmt=*) ' igraup          =',igraup
      write (unit=*,fmt=*) ' ihail           =',ihail
      write (unit=*,fmt=*) ' cparm           =',cparm
      write (unit=*,fmt=*) ' rparm           =',rparm
      write (unit=*,fmt=*) ' pparm           =',pparm
      write (unit=*,fmt=*) ' sparm           =',sparm
      write (unit=*,fmt=*) ' aparm           =',aparm
      write (unit=*,fmt=*) ' gparm           =',gparm
      write (unit=*,fmt=*) ' hparm           =',hparm
      write (unit=*,fmt=*) ' gnu             =',gnu
      call abort_run('Error reading namelist...','read_nl','rdint.f90')
   else
      !----- Switch closure type to lower case. -------------------------------------------!
      caux(1) = closure_type
      call tolower(caux,1)
      closure_type = caux(1)
   end if
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !      Read in ED2_INFO namelist.                                                       !
   !---------------------------------------------------------------------------------------!
   write (unit=*,fmt='(a)') 'Reading ED2 namelist information'
   call read_ednl(iunit,filename)
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !      Read in MODEL_OPTIONS namelist.                                                  !
   !---------------------------------------------------------------------------------------!
   read (unit=iunit, iostat=err, nml=MODEL_SOUND)
   if (err /= 0) then
      write (unit=*,fmt='(a)')        '--------------------------------------------------'
      write (unit=*,fmt='(a)')        ' ERROR reading the namelist! '
      write (unit=*,fmt='(a,1x,a)')   ' Namelist file   : ',trim(filename)
      write (unit=*,fmt='(a,1x,a)')   ' Namelist section: ','MODEL_SOUND'
      write (unit=*,fmt='(a)')        ' Compare values read with file contents...'
      write (unit=*,fmt='(a)')        '--------------------------------------------------'
      write (unit=*,fmt='(a)')        ''
      write (unit=*,fmt=*) ' ipsflg  =',ipsflg
      write (unit=*,fmt=*) ' itsflg  =',itsflg
      write (unit=*,fmt=*) ' irtsflg =',irtsflg
      write (unit=*,fmt=*) ' iusflg  =',iusflg
      write (unit=*,fmt=*) ' hs      =',hs
      write (unit=*,fmt=*) ' ps      =',ps
      write (unit=*,fmt=*) ' ts      =',ts
      write (unit=*,fmt=*) ' rts     =',rts
      write (unit=*,fmt=*) ' us      =',us
      write (unit=*,fmt=*) ' vs      =',vs
      write (unit=*,fmt=*) ' co2s    =',co2s
      call abort_run('Error reading namelist...','read_nl','rdint.f90')
   end if
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !      Read in MODEL_PRINT namelist.                                                    !
   !---------------------------------------------------------------------------------------!
   read (unit=iunit, iostat=err, nml=MODEL_PRINT)
   if (err /= 0) then
      write (unit=*,fmt='(a)')        '--------------------------------------------------'
      write (unit=*,fmt='(a)')        ' ERROR reading the namelist! '
      write (unit=*,fmt='(a,1x,a)')   ' Namelist file   : ',trim(filename)
      write (unit=*,fmt='(a,1x,a)')   ' Namelist section: ','MODEL_PRINT'
      write (unit=*,fmt='(a)')        ' Compare values read with file contents...'
      write (unit=*,fmt='(a)')        '--------------------------------------------------'
      write (unit=*,fmt='(a)')        ''
      write (unit=*,fmt=*) ' nplt   = ',nplt
      write (unit=*,fmt=*) ' iplfld = ',(trim(iplfld(i))//';', i=1,size(iplfld))
      write (unit=*,fmt=*) ' ixsctn = ',ixsctn
      write (unit=*,fmt=*) ' isbval = ',isbval
      call abort_run('Error reading namelist...','read_nl','rdint.f90')
   end if
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !      Read in ISAN_CONTROL namelist.                                                   !
   !---------------------------------------------------------------------------------------!
   read (unit=iunit, iostat=err, nml=ISAN_CONTROL)
   if (err /= 0) then
      write (unit=*,fmt='(a)')        '--------------------------------------------------'
      write (unit=*,fmt='(a)')        ' ERROR reading the namelist! '
      write (unit=*,fmt='(a,1x,a)')   ' Namelist file   : ',trim(filename)
      write (unit=*,fmt='(a,1x,a)')   ' Namelist section: ','ISAN_CONTROL'
      write (unit=*,fmt='(a)')        ' Compare values read with file contents...'
      write (unit=*,fmt='(a)')        '--------------------------------------------------'
      write (unit=*,fmt='(a)')        ''
      write (unit=*,fmt=*) ' iszstage  =',iszstage
      write (unit=*,fmt=*) ' ivrstage  =',ivrstage
      write (unit=*,fmt=*) ' isan_inc  =',isan_inc
      write (unit=*,fmt=*) ' guess1st  =',guess1st
      write (unit=*,fmt=*) ' i1st_flg  =',i1st_flg
      write (unit=*,fmt=*) ' iupa_flg  =',iupa_flg
      write (unit=*,fmt=*) ' isfc_flg  =',isfc_flg
      write (unit=*,fmt=*) ' iapr      =',trim(iapr)
      write (unit=*,fmt=*) ' iarawi    =',trim(iarawi)
      write (unit=*,fmt=*) ' iasrfce   =',trim(iasrfce)
      write (unit=*,fmt=*) ' varpfx    =',trim(varpfx)
      write (unit=*,fmt=*) ' ioflgisz  =',ioflgisz
      write (unit=*,fmt=*) ' ioflgvar  =',ioflgvar
      call abort_run('Error reading namelist...','read_nl','rdint.f90')
   end if
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !      Read in ISAN_ISENTROPIC namelist.                                                !
   !---------------------------------------------------------------------------------------!
   read (unit=iunit,iostat=err,nml=ISAN_ISENTROPIC)
   if (err /= 0) then
      write (unit=*,fmt='(a)')        '--------------------------------------------------'
      write (unit=*,fmt='(a)')        ' ERROR reading the namelist! '
      write (unit=*,fmt='(a,1x,a)')   ' Namelist file   : ',trim(filename)
      write (unit=*,fmt='(a,1x,a)')   ' Namelist section: ','ISAN_ISENTROPIC'
      write (unit=*,fmt='(a)')        ' Compare values read with file contents...'
      write (unit=*,fmt='(a)')        '--------------------------------------------------'
      write (unit=*,fmt='(a)')        ''
      write (unit=*,fmt=*) 'nisn     =',nisn
      write (unit=*,fmt=*) 'levth    =',levth
      write (unit=*,fmt=*) 'nigrids  =',nigrids
      write (unit=*,fmt=*) 'topsigz  =',topsigz
      write (unit=*,fmt=*) 'hybbot   =',hybbot
      write (unit=*,fmt=*) 'hybtop   =',hybtop
      write (unit=*,fmt=*) 'sfcinf   =',sfcinf
      write (unit=*,fmt=*) 'sigzwt   =',sigzwt
      write (unit=*,fmt=*) 'nfeedvar =',nfeedvar
      write (unit=*,fmt=*) 'maxsta   =',maxsta
      write (unit=*,fmt=*) 'maxsfc   =',maxsfc
      write (unit=*,fmt=*) 'notsta   =',notsta
      write (unit=*,fmt=*) 'notid    =',(trim(notid(i))//';',i=1,size(notid))
      write (unit=*,fmt=*) 'iobswin  =',iobswin
      write (unit=*,fmt=*) 'stasep   =',stasep
      write (unit=*,fmt=*) 'igridfl  =',igridfl
      write (unit=*,fmt=*) 'gridwt   =',gridwt
      write (unit=*,fmt=*) 'gobsep   =',gobsep
      write (unit=*,fmt=*) 'gobrad   =',gobrad
      write (unit=*,fmt=*) 'wvlnth   =',wvlnth
      write (unit=*,fmt=*) 'swvlnth  =',swvlnth
      write (unit=*,fmt=*) 'respon   =',respon
      call abort_run('Error reading namelist...','read_nl','rdint.f90')
   end if
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !      Close the namelist file.                                                         !
   !---------------------------------------------------------------------------------------!
   close(iunit,status='keep',iostat=err)
   if (err /= 0) then
      write (unit=*,fmt='(a,1x,a)')   ' Namelist file: ',trim(filename)
      write (unit=*,fmt='(a,1x,i10)') ' Iostat       : ',err
      call abort_run(' Error closing the namelist...','read_nl','rdint.f90')
   end if


   !------ Find the final and history time in seconds. ------------------------------------!
   call date_2_seconds (iyearz,imonthz,idatez,itimez*100,iyeara,imontha,idatea,itimea*100  &
                       ,timmax)
   call date_2_seconds (iyearh,imonthh,idateh,itimeh*100,iyeara,imontha,idatea,itimea*100  &
                       ,timstr)

   if (isfcl /= 5) then
      !---- Not an ED-BRAMS run, and isoilflg/ivegtflg are set to 3, switch them by 1. ----!
      where (isoilflg == 3) 
         isoilflg = 1
      end where
      where (ivegtflg == 3) 
         ivegtflg = 1
      end where
   end if
   !---- If someone accidentally defined itoptflg, isstflg or ndviflg to 3, make them 1. --!
   where (itoptflg == 3) itoptflg = 1
   where (isstflg  == 3) isstflg  = 1
   where (ndviflg  == 3) ndviflg  = 1

   !---------------------------------------------------------------------------------------!
   !     Save the moisture complexity level into logical variables. Note that vapour_on is !
   ! not level == 1, it will be true when level is 1, 2, and 3. (whenever vapour is on).   !
   ! Likewise, cloud_on will be true when level is either 2 or 3. Bulk microphysics will   !
   ! be true only when level >= 3 (levels = 4 and 5 exist too but rarely used).            !
   !---------------------------------------------------------------------------------------!
   vapour_on  = level >= 1
   cloud_on   = level >= 2
   bulk_on    = level >= 3
   vapour_on8 = vapour_on
   cloud_on8  = cloud_on
   bulk_on8   = bulk_on
   level8     = level
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    Save the CO2 complexity level into a logical variable.                             !
   !---------------------------------------------------------------------------------------!
   co2_on    = ico2 > 0
   !---------------------------------------------------------------------------------------!

   return
end subroutine read_nl
!==========================================================================================!
!==========================================================================================!

