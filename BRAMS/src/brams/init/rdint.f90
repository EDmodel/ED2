!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

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
  use mem_turb, only:if_urban_canopy
  use io_params
  use micphys
  use therm_lib, only : level
  use mem_turb, only : turb_g

  ! CATT
  use catt_start, only         : CATT                      ! intent(in)
  use mem_scalar, only         : &
       scalar_g, &                                         ! intent(inout)
       RECYCLE_TRACERS                                     ! intent(in)
  use emission_source_map, only: read_emission_sources_map ! Subroutine

  ! CARMA Radiation
  use mem_globrad, only: master_read_carma_data ! Subroutine

  ! TEB_SPM
  use teb_spm_start, only : TEB_SPM   ! intent(in)
  use mem_teb, only       : teb_g     ! intent(inout)
  use mem_teb_common, only: tebc_g    ! intent(inout)
  use teb_vars_const, only: iteb      ! intent(in)
  ! TEB_SPM: emission and chemistry
  use mem_gaspart, only   : gaspart_g ! intent(inout)
  use mem_emiss, only     : &
       ichemi,   &                    ! intent(inout)
       isource,  &                    ! intent(inout)
       ichemi_in                      ! intent(inout)

  ! Soil Moisture Init.
  use mem_soil_moisture, only : SOIL_MOIST ! INTENT(IN)

  use ref_sounding, only : dn01dn

  implicit none

  character(len=*), intent(in) :: name_name
 
  !---------------------------------------------------------------------
  !     *** This routine is the driver of all initialization packages
  !---------------------------------------------------------------------

  character(len=8) :: rest
  integer :: ifm,icm,mynum,ihm,ngr,nv,ierr

  !     Set version number for common blocks.
  iversion = 2

  !     Unit numbers for some I/O files.
  iopunt=6

  !CALL jnmbinit(level) ! Correctd by Craig.

  if (runtype(1:7) == 'INITIAL') then

     time=0.
     ngbegun(1:ngrids) = 0

     !----------------------------------------------------------------------
     !                 Initial startup
     !----------------------------------------------------------------------

     ! Read surface, topo, sst, and ndvi files for all grids. All the files
     !   were checked earlier, so they must be correct.,

     do ifm = 1,ngrids
        call top_read(ifm)
     enddo

     do ifm = 1,ngrids
        call sfc_read(ifm)
     enddo

     !     Define grid topography, transform, latitude-longitude,
     !        and map factor arrays.

     call grid_setup(2)

     ! read SST files

     do ifm = 1,ngrids
        call sst_read(1,ifm,ierr)
     enddo

     ! read NDVI files

     do ifm = 1,ngrids
        call ndvi_read(1,ifm,ierr)
     enddo

     ! Initialize snowcover arrays

     do ifm = 1,ngrids
        call snowinit(nnxp(ifm),nnyp(ifm)  &
             ,leaf_g(ifm)%snow_mass,leaf_g(ifm)%snow_depth)
     enddo

     ! TEB_SPM
     if (TEB_SPM==1) then
        ! read FUSO (Local Time) files
        do ifm = 1,ngrids
           call fuso_read(ifm)
        enddo
     endif

     ! The following things will be done for INITAIL = 1 or 3...

     ! Initial = 1, Horizontal homogenous ATM
     ! =2, Varfile initialization
     ! =3, HFILIN initialization

     if ( initial == 1 .or. initial == 3) then

        ! If horizontally homogeneous initialization,
        !    subroutine INITHH loops through all grids and initializes
        !    those for which nxtnest = 0.

        if(initial == 1) then
           print*,'Horizontally-homogeneous-INITIAL start of grid- 1'
           call inithh()
        endif

        !If "history" initialization, call INITHIS.
        !This will define initial fields and reference state on grid 1 from
        !history file. Other grids will be interpolated as in a INITIAL=1 start

        if (initial == 3) then
           print*,'History-INITIAL start of grid- 1'
           call inithis()
        endif


        !  On all fine grids, initialize the surface layer characteristics,
        !  the 1-D reference state arrays, the 3-D reference state arrays,
        !  and the prognostic atmospheric fields by interpolation.

        call fmrefs1d(2,ngrids)

        do ifm = 2,ngrids
           icm = nxtnest(ifm)
           if (icm  >=  1) then
              call fmrefs3d(ifm,mynum)

              call prgintrp(nnzp(icm),nnxp(icm),nnyp(icm)  &
                   ,nnzp(icm),nnxp(icm),nnyp(icm),0,0,ifm,1,mynum)

              call fmdn0(ifm)

              print*,'Initial interpolation of grid-',ifm
           endif
        enddo


     elseif(initial == 2) then

        ! If "variable initialization", do it all here

        call varf_read(0)

     endif

     !     Initialize past time level velocity and perturbation Exner function
     !     on all grids.

     do ifm=1,ngrids
        call newgrid(ifm)
        call fldinit(1)
        call negadj1(nzp,nxp,nyp,1,nxp,1,nyp)
        call thermo(nzp,nxp,nyp,1,nxp,1,nyp)

        if (level  ==  3) then
           call azero3(nzp*nxp*nyp,scratch%vt3da,scratch%vt3dg,scratch%vt3dh)
           call azero2(nzp*nxp*nyp,scratch%vt3dc,scratch%vt3di)
           !----- Using scratch variables to define cccnp and cifnp -----------------------!
           call initqin(nzp,nxp,nyp   &
                ,scratch%vt3da        &
                ,scratch%vt3dg        &
                ,scratch%vt3dh        &
                ,basic_g(ifm)%pi0     &
                ,basic_g(ifm)%pp      &
                ,basic_g(ifm)%theta   &
                ,basic_g(ifm)%dn0     &
                ,scratch%vt3dc        &
                ,scratch%vt3di        )
           !----- Copying them to the micro arrays if they are allocated ------------------!
           if (irain  >= 1) call atob(nzp*nxp*nyp,scratch%vt3da,micro_g(ifm)%rrp)
           if (igraup >= 1) call atob(nzp*nxp*nyp,scratch%vt3dg,micro_g(ifm)%rgp)
           if (ihail  >= 1) call atob(nzp*nxp*nyp,scratch%vt3dh,micro_g(ifm)%rhp)
           if (icloud == 7) call atob(nzp*nxp*nyp,scratch%vt3dc,micro_g(ifm)%cccnp)
           if (ipris  == 7) call atob(nzp*nxp*nyp,scratch%vt3di,micro_g(ifm)%cifnp)
        end if

     end do

     ! If initializing some fields from previous runs...

     ! Do recycle procedure for normal RAMS recycle or only to do
     ! tracers assimilation
     if ((ipastin == 1) .or. (CATT==1 .and. RECYCLE_TRACERS==1)) then
        call recycle()
     endif

     ! Fill land surface data for all grids that have no standard input files

     ! ALF - For use with SiB
     select case (isfcl)
     case (1,2,5)
        call sfcdata
     case (3)
        call sfcdata_sib_driver
     end select


     ! Initialize various LEAF variables.
     if (ipastin == 0) call geonest_nofile(1,ngrids)

     ! TEB
     if (TEB_SPM==1) then

        ! Initial values for Common use TEB vars 
        ! For now, it is only being used to zero out this four common 
        ! use variables.
        ! This variables will be used for other purposes later.
        ! Edmilson D. Freitas 07/07/2006

        do ifm = 1,ngrids
           call TEBC_INIT(ifm,nnxp(ifm),nnyp(ifm),npatch &
                ,leaf_g(ifm)%leaf_class                  &
                ,leaf_g(ifm)%G_URBAN                     &
                ,tebc_g(ifm)%EMIS_TOWN                   &
                ,tebc_g(ifm)%ALB_TOWN                    &
                ,tebc_g(ifm)%TS_TOWN                     )
        enddo

        if (iteb==1) then
           do ifm = 1,ngrids
              call TEB_INIT(ifm,nnzp(ifm),nnxp(ifm),nnyp(ifm),npatch          &
                   ,leaf_g(ifm)%leaf_class        ,basic_g(ifm)%theta         &
                   ,basic_g(ifm)%rv               ,basic_g(ifm)%pi0           &
                   ,basic_g(ifm)%pp                                           &
                   ,teb_g(ifm)%T_ROOF             ,teb_g(ifm)%T_ROAD          &
                   ,teb_g(ifm)%T_WALL             ,teb_g(ifm)%TI_BLD          &
                   ,teb_g(ifm)%TI_ROAD            ,teb_g(ifm)%T_CANYON        &
                   ,teb_g(ifm)%R_CANYON           ,teb_g(ifm)%TS_ROOF         &
                   ,teb_g(ifm)%TS_ROAD            ,teb_g(ifm)%TS_WALL         &
                   ,teb_g(ifm)%H_TRAFFIC          ,teb_g(ifm)%LE_TRAFFIC      &
                   ,teb_g(ifm)%H_INDUSTRY         ,teb_g(ifm)%LE_INDUSTRY     &
                   ,teb_g(ifm)%WS_ROOF            ,teb_g(ifm)%WS_ROAD         &
                   ,tebc_g(ifm)%EMIS_TOWN         ,tebc_g(ifm)%ALB_TOWN       &
                   ,tebc_g(ifm)%TS_TOWN           ,leaf_g(ifm)%G_URBAN        )
           enddo
        endif

        ! Initialize gases and particulate matter
        if (isource==1) then
           do ifm = 1,ngrids
              call init_conc1(1,ifm,nnzp(ifm),nnxp(ifm),nnyp(ifm),npatch     &
                   ,leaf_g(ifm)%G_URBAN         ,gaspart_g(ifm)%pno          &
                   ,gaspart_g(ifm)%pno2         ,gaspart_g(ifm)%ppm25        &
                   ,gaspart_g(ifm)%pco          ,gaspart_g(ifm)%pvoc         &
                   ,gaspart_g(ifm)%pso2         ,gaspart_g(ifm)%pso4         &
                   ,gaspart_g(ifm)%paer         ,zt                          )
              if (ichemi==1) then       !calling more added scalars for chemistry
                 if (ichemi_in==1) then !reading init.values from a previous run
                    call init_conc_prev(name_name)
                 else
                    call init_conc2(1,ifm,nnzp(ifm),nnxp(ifm),nnyp(ifm),npatch &
                         ,leaf_g(ifm)%G_URBAN    ,gaspart_g(ifm)%po3           &
                         ,gaspart_g(ifm)%prhco   ,gaspart_g(ifm)%pho2          &
                         ,gaspart_g(ifm)%po3p    ,gaspart_g(ifm)%po1d          &
                         ,gaspart_g(ifm)%pho     ,gaspart_g(ifm)%proo          &
                         ,zt)
                 end if
              end if
           end do
        end if

     endif

     ! Read Radiation Parameters if CARMA Radiation is selected
     call master_read_carma_data()

     ! Initialise turbulence factor akscal
     if (if_adap == 1) then
        do ifm=1,ngrids
           call akscal_init(nnxp(ifm),nnyp(ifm),ifm,grid_g(ifm)%topta,turb_g(ifm)%akscal)
        end do
     else
        do ifm=1,ngrids
           call akscal_init(nnxp(ifm),nnyp(ifm),ifm,grid_g(ifm)%topt,turb_g(ifm)%akscal)
        end do
     end if

     ! CATT
     if (CATT==1) then

        ! srf......................................................
        ! Initiating tracers concentration and sources
        if (RECYCLE_TRACERS==0) then
           do ifm=1,NGRIDS
              ! Do not perform recycle procedure for tracers
              ! Initiating with zeros
              ! Putting zero on specific scalar vectors
              scalar_g(1,ifm)%sclp(:,:,:) = 0.
              scalar_g(2,ifm)%sclp(:,:,:) = 0.
              scalar_g(3,ifm)%sclp(:,:,:) = 0.
              scalar_g(4,ifm)%sclp(:,:,:) = 0.
              !!scalar_g(5,ifm)%sclp(:,:,:) = 0.
           enddo
        endif

        ! Reading Emission Maps for all grids
        call read_emission_sources_map()

     endif
     !

     if (initial == 3) call sfcinit_hstart()

     !hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh

  elseif(runtype(1:7) == 'HISTORY') then

     !hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh

     !                  History file start
     !                  -----------------

     call history_start(name_name)

     call grid_setup(1)

     ! Check surface,topo,sst,ndvi files. Remake if necessary.

     call make_sfcfiles()

     ! Read surface and topo files for any added grids

     do ifm = ngridsh+1,ngrids
        call sfc_read(ifm)
     enddo

     do ifm = ngridsh+1,ngrids
        call top_read(ifm)
     enddo

     call grid_setup(2)

     ! Read in sst and ndvi files for all grids

     do ifm = 1,ngrids
        call sst_read(1,ifm,ierr)
     enddo

     do ifm = 1,ngrids
        call ndvi_read(1,ifm,ierr)
     enddo

     ! TEB
     if (TEB_SPM==1) then
        ! Read FUSO (Local Time) files for any added grids
        do ifm = ngridsh+1,ngrids
           call fuso_read(ifm)
        enddo
     endif

     do ifm = 1,ngrids
        icm = nxtnest(ifm)
        if (icm  ==  0) then
           call newgrid(ifm)
           call refs3d (nzp,nxp,nyp                     &
                ,basic_g(ifm)%pi0  ,basic_g(ifm)%dn0    &
                ,basic_g(ifm)%dn0u ,basic_g(ifm)%dn0v   &
                ,basic_g(ifm)%th0  ,grid_g(ifm)%topt    &
                ,grid_g(ifm)%rtgt  )
        endif
     enddo

     do ifm = 1,min(ngrids,ngridsh)
        icm = nxtnest(ifm)
        if (icm  >  0) call fmrefs3d(ifm,0)
        call negadj1(nzp,nxp,nyp,1,nxp,1,nyp)
     enddo

     ! ALF - For use with SiB
     select case (isfcl)
     case (1,2,5)
        call sfcdata
     case (3)
        call sfcdata_sib_driver
     end select

     ! Heterogenous Soil Moisture Init.
     if ((SOIL_MOIST == 'h').or.(SOIL_MOIST == 'H').or.  &
          (SOIL_MOIST == 'a').or.(SOIL_MOIST == 'A')) then

        do ifm = 1,min(ngrids,ngridsh)
           call newgrid(ifm)
           call soil_moisture_init(nnzp(ifm),nnxp(ifm),nnyp(ifm)              &
                ,nzg,nzs,npatch,ifm                                           &
                ,basic_g(ifm)%theta            ,basic_g(ifm)%pi0              &
                ,basic_g(ifm)%pp               ,leaf_g(ifm)%soil_water        &
                ,leaf_g(ifm)%soil_energy       ,leaf_g(ifm)%soil_text         &
                ,leaf_g(ifm)%sfcwater_mass     ,leaf_g(ifm)%sfcwater_energy   &
                ,leaf_g(ifm)%sfcwater_depth    ,grid_g(ifm)%glat              &
                ,grid_g(ifm)%glon              ,grid_g(ifm)%flpw              )
        enddo

     endif

     !     If any grids are being added for this run, initialize their
     !     surface layer variables, 1-D reference state variables, and
     !     prognostic atmospheric and soil model fields.

     if (ngrids  >  ngridsh) then
        print*,' +-------------------------------------'
        print*,'            !      New grids will be added.       '
        print*,'            !'
        print*,'            ! ',ngridsh,' grid(s) on history file.'
        print*,'            ! ',ngrids, ' grids to be run.        '
        print*,' +-------------------------------------'
        call fmrefs1d(ngridsh+1,ngrids)
        do ifm = ngridsh+1,ngrids
           icm = nxtnest(ifm)
           if (icm  ==  0) then
              call abort_run('Attempted to add hemispheric grid on a history restart' &
                            ,'initlz','rdint.f90')
           endif
           call fmrefs3d(ifm,0)
           call prgintrp(nnzp(icm),nnxp(icm),nnyp(icm)  &
                ,nnzp(icm),nnxp(icm),nnyp(icm),0,0,ifm,1,mynum)
           print*,'History start interpolation of added grid-',ngrid

           call fmdn0(ifm)
           call newgrid(ifm)
           call fldinit(0)
           call negadj1(nzp,nxp,nyp,1,nxp,1,nyp)
           call thermo(nzp,nxp,nyp,1,nxp,1,nyp)
           if (level  ==  3) then
              call initqin(nzp,nxp,nyp                         &
                   ,micro_g(ifm)%q2      ,micro_g(ifm)%q6      &
                   ,micro_g(ifm)%q7      ,basic_g(ifm)%pi0     &
                   ,basic_g(ifm)%pp      ,basic_g(ifm)%theta   &
                   ,basic_g(ifm)%dn0     ,micro_g(ifm)%cccnp   &
                   ,micro_g(ifm)%cifnp   )
           endif

           ! Heterogenous Soil Moisture Init.
           if ((SOIL_MOIST == 'h').or.(SOIL_MOIST == 'H').or.  &
              (SOIL_MOIST == 'a').or.(SOIL_MOIST == 'A')) then
              call soil_moisture_init(nnzp(ifm),nnxp(ifm),nnyp(ifm)            &
                   ,nzg,nzs,npatch,ifm                                         &
                   ,basic_g(ifm)%theta          ,basic_g(ifm)%pi0              &
                   ,basic_g(ifm)%pp             ,leaf_g(ifm)%soil_water        &
                   ,leaf_g(ifm)%soil_energy     ,leaf_g(ifm)%soil_text         &
                   ,leaf_g(ifm)%sfcwater_mass   ,leaf_g(ifm)%sfcwater_energy   &
                   ,leaf_g(ifm)%sfcwater_depth  ,grid_g(ifm)%glat              &
                   ,grid_g(ifm)%glon            ,grid_g(ifm)%flpw              )
           endif

          ! Initialise turbulence factor akscal
          if (if_adap == 1) then
             call akscal_init(nnxp(ifm),nnyp(ifm),ifm,grid_g(ifm)%topta,turb_g(ifm)%akscal)
          else
             call akscal_init(nnxp(ifm),nnyp(ifm),ifm,grid_g(ifm)%topt,turb_g(ifm)%akscal)
          end if

        enddo

        !Fill land surface data for all grids that have no standard input files
        call geonest_nofile(ngridsh+1,ngrids)

     elseif (ngrids  <  ngridsh) then
        print*,' +-------------------------------------'
        print*,'            !      Grids will be subtracted.       '
        print*,'            !'
        print*,'            ! ',NGRIDSH,' grid(s) on history file.'
        print*,'            ! ',NGRIDS, ' grids to be run.        '
        print*,' +-------------------------------------'
     endif

     ! Read Radiation Parameters if CARMA Radiation is selected
     call master_read_carma_data()
     !

     ! CATT
     if (CATT==1) then
        !srf......................................................
        ! Reading Emission Maps for all grids
        call read_emission_sources_map()
     endif

  else
     call abort_run('Wrong runtype in initlz: '//trim(runtype)//'...' &
                   ,'initlz','rdint.f90')

  endif

  ! For a global model domain, initialize interpolation table values for
  ! communication between hemispheric grids.  Subroutine hemintrp_cof will
  ! return immediately if nhemgrd2 is not greater than 1.

  call newgrid(1)
  ihm = nnxyp(1)
  call hemintrp_cof (nnxp(1),nnyp(1)                                                       &
                    ,scratch%scr1((1+0*ihm):(1*ihm))      ,scratch%scr1((1+1*ihm):(2*ihm)) &
                    ,scratch%scr1((1+2*ihm):(3*ihm))      ,scratch%scr1((1+3*ihm):(4*ihm)) &
                    ,scratch%scr1((1+4*ihm):(5*ihm))      ,scratch%scr1((1+5*ihm):(6*ihm)) )

  !----------------------------------------------------------------------------------------!
  !    Initialise a bunch of microphysics parameters. The new input is the density at the  !
  ! top, to find the minimum concentration of hydrometeors to be considered.               !
  !----------------------------------------------------------------------------------------!
  call micro_master(dn01dn(nnzp(1),1))
  !----------------------------------------------------------------------------------------!

  !       Fill latitude-longitude, map factor, and Coriolis arrays.

  do ifm = 1,ngrids
     call newgrid(ifm)
     call fcorio(nxp,nyp        &
          ,basic_g(ifm)%fcoru   &
          ,basic_g(ifm)%fcorv   &
          ,grid_g(ifm)%glat     )
  enddo


  !  If we are doing one-way nesting or varfile nudging, inventory,
  !     prepare history/varfile files
  !     and fill past/future nudging arrays for start of simulation

  if ( nud_type == 1 ) then
     call nud_read(1)
  elseif(nud_type == 2) then
     call varf_read(1)
  endif

  ! Do same if doing condensate nudging

  if ( nud_cond == 1 ) call cond_read(1)

  ! Process and read observations for ODA - observational data assimilation

  if (if_oda == 1) call oda_read()

  ! Read cumulus heating fields

  if (if_cuinv == 1) call cu_read(1)

  ! Initialize urban canopy drag coefficients

  if (if_urban_canopy == 1) call urb_drag_init()


  ! Print locations of all grids

  call gridloc_prt()

  !                  Save initial fields on history and analysis files
  !                  -------------------------------------------------
  iinput=ioutput
  rest  ='no'
  if (runtype  ==  'HISTORY') rest = 'yes'
  call hiswrt(rest)
  call anlwrt(rest,'INST')
  if(frqlite > 0.) call anlwrt(rest,'LITE')

  !                  Save initial fields into the averaged arrays
  !                  -------------------------------------------------
  if(avgtim /= 0.)then

     do ngr=1,ngrids

        do nv=1,num_var(ngr)
           if(vtab_r(nv,ngr)%imean == 1) then
              call atob(vtab_r(nv,ngr)%npts,vtab_r(nv,ngr)%var_p  &
                   ,vtab_r(nv,ngr)%var_m)
           endif
        enddo
     enddo
  endif


  !                  Print the header information and initial fields
  !                  -----------------------------------------------
  ngrid=1

  call prtopt(6)

!!$  call nameout()

  if (initfld  ==  1) then
     do ifm = 1,ngrids
        call newgrid(ifm)
        call prtout()
     enddo
  endif

  call opspec3()

  return
end subroutine initlz

!*****************************************************************************

subroutine read_nl(file)

  ! read_nl:
  !   reads namelist and checks consistency of a few options


  implicit none

  character(len=*), intent(in) :: file
  character(len=*), parameter :: h="**(read_nl)**"


  ! Deactivate namelist scanning and parsing
  ! Activate standard (fortran 90) namelist read
  ! Estabilishes default options and option consistency

  call ReadNamelist(file)

end subroutine read_nl



subroutine ReadNamelist(fileName)

  ! ReadNamelist:
  !    open, reads and close namelist file
  !    implements defaults for namelist variables
  !    check input options consistency

  use io_params, only: frqboth, &
       afilout,                 &
       avgtim,                  &
       frqanl,                  &
       frqhis,                  &
       frqlite,                 &
       frqmean,                 &
       frqprt,                  &
       hfilin,                  &
       hfilout,                 &
       iclobber,                &
       ihistdel,                &
       initfld,                 &
       ioutput,                 &
       ipastin,                 &
       iplfld,                  &
       isbval,                  &
       isoilflg,                &
       isoilfn,                 &
       isstflg,                 &
       isstfn,                  &
       itopsflg,                &
       itoptflg,                &
       itoptfn,                 &
       iupdndvi,                &
       iupdsst,                 &
       ivegtflg,                &
       ivegtfn,                 &
       ixsctn,                  &
       iz0flg,                  &
       kwrite,                  &
       lite_vars,               &
       ndviflg,                 &
       ndvifn,                  &
       ndvifpfx,                &
       nlite_vars,              &
       nofilflg,                &
       nplt,                    &
       pastfn,                  &
       sfcfiles,                &
       sstfpfx,                 &
       timstr,                  &
       imonthh,                 & !
       iyearh,                  & !
       idateh,                  & !
       itimeh,                  & !
       topfiles,                &
       toptenh,                 &
       toptwvl,                 &
       xlite,                   &
       ylite,                   &
       z0fact,                  &
       z0max,                   &
       zlite,                   &
       ! TEB
       ifusflg,                 &
       ifusfn,                  &
       fusfiles
  use isan_coms, only: gobrad, &
       gobsep, &
       gridwt, &
       guess1st, &
       hybbot, &
       hybtop, &
       i1st_flg, &
       iapr, &
       iarawi, &
       iasrfce, &
       igridfl, &
       iobswin, &
       ioflgisz, &
       ioflgvar, &
       isan_inc, &
       isfc_flg, &
       iszstage, &
       iupa_flg, &
       ivrstage, &
       levth, &
       maxsfc, &
       maxsta, &
       nfeedvar, &
       nigrids, &
       nisn, &
       notid, &
       notsta, &
       respon, &
       sfcinf, &
       sigzwt, &
       stasep, &
       swvlnth, &
       topsigz, &
       varpfx, &
       wvlnth
  use mem_cuparm, only: confrq, &
       cu_prefix, &
       cu_tel, &
       cu_til, &
       if_cuinv, &
       nnqparm, &
       tcu_beg, &
       tcu_end, &
       tnudcu, &
       wcldbs, &
       wt_cu_grid, &
       nclouds, &
       ndeepest, &
       nshallowest, &
       cptime
  use mem_globrad, only: raddatfn, & !CARMA
       rdatfnum                      !CARMA
  use grell_coms, only:  &
          closure_type,  & ! intent(out)
          maxclouds,     & ! intent(out)
          iupmethod,     & ! intent(out)
          depth_min,     & ! intent(out)
          cap_maxs,      & ! intent(out)
          maxens_lsf,    & ! intent(out)
          maxens_eff,    & ! intent(out)
          maxens_dyn,    & ! intent(out)
          maxens_cap,    & ! intent(out)
          iupmethod,     & ! intent(out)
          iupstrm,       & ! intent(out)
          radius,        & ! intent(out)
          zkbmax,        & ! intent(out)
          max_heat,      & ! intent(out)
          zcutdown,      & ! intent(out)
          z_detr         ! ! intent(out)
  use mem_grid, only: centlat, &
       centlon, &
       cphas, &
       deltax, &
       deltay, &
       deltaz, &
       distim, &
       dtlong, &
       dzmax, &
       dzrat, &
       expnme, &
       gridu, &
       gridv, &
       ibnd, &
       icorflg, &
       idatea, &
       idatez,  &
       ideltat, &
       if_adap, &
       ihtran, &
       initial, &
       imontha, &
       imonthz, &
       itimea,  &
       ihoura,  &
       itimez,  &
       iyeara,  &
       iyearz,  &
       jbnd, &
       lsflg, &
       nacoust, &
       naddsc, &
       nestz1, &
       nestz2, &
       nfpt, &
       ngrids, &
       ninest, &
       njnest, &
       nknest, &
       nndtrat, &
       nnstbot, &
       nnsttop, &
       nnxp, &
       nnyp, &
       nnzp, &
       npatch, &
       nstratx, &
       nstraty, &
       nstratz1, &
       nstratz2, &
       nxtnest, &
       nzg, &
       nzs, &
       polelat, &
       polelon, &
       runtype, &
       timeunit, &
       timmax, &
       zz
  use mem_leaf, only: albedo, &
       drtcon, &
       dthcon, &
       isfcl, &
       nslcon, &
       nvegpat, &
       nvgcon, &
       pctlcon, &
       seatmp, &
       slmstr, &
       slz, &
       stgoff, &
       zrough
  use mem_oda, only: frqoda, &
       if_oda, &
       oda_sfc_tel, &
       oda_sfc_til, &
       oda_sfcprefix, &
       oda_upa_tel, &
       oda_upa_til, &
       oda_upaprefix, &
       roda_hgt, &
       roda_sfc0, &
       roda_sfce, &
       roda_upa0, &
       roda_upae, &
       roda_zfact, &
       tnudoda, &
       todabeg, &
       todaend, &
       wt_oda_grid, &
       wt_oda_pi, &
       wt_oda_rt, &
       wt_oda_th, &
       wt_oda_uv
  use mem_radiate, only: ilwrtyp, &
       iswrtyp, &
       icumfdbk, & 
       lonrad, &
       radfrq
  use mem_soil_moisture, only: soil_moist, &
       soil_moist_fail, &
       usdata_in, &
       usmodel_in
  use mem_turb, only: akmin, &
       akmax, &
       hgtmin, &
       hgtmax, &
       csx, &
       csz, &
       idiffk, &
       ibotflx, &
       ibruvais, &
       if_urban_canopy, &
       ihorgrad, &
       xkhkm, &
       zkhkm
  use mem_varinit, only: cond_hfile, &
       nud_cond, &
       nud_hfile, &
       nud_type, &
       nudlat, &
       t_nudge_rc, &
       tcond_beg, &
       tcond_end, &
       tnudcent, &
       tnudlat, &
       tnudtop, &
       varfpfx, &
       vwait1, &
       vwaittot, &
       wt_nudge_grid, &
       wt_nudge_pi,  &
       wt_nudge_rt,  &
       wt_nudge_th,  &
       wt_nudge_uv,  &
       wt_nudge_co2, &
       wt_nudgec_grid, &
       znudtop
  use micphys, only: aparm, &
       coltabfn, &
       cparm, &
       gnu, &
       gparm, &
       hparm, &
       iaggr, &
       icloud, &
       igraup, &
       ihail, &
       ipris, &
       irain, &
       isnow, &
       mkcoltab, &
       pparm, &
       rparm, &
       sparm
  use node_mod, only: load_bal
  use ref_sounding, only: hs, &
       ipsflg, &
       irtsflg, &
       itsflg, &
       iusflg, &
       ps, &
       rts, &
       ts, &
       us, &
       vs

  ! CATT
  use catt_start, only: CATT
  use emission_source_map, only: firemapfn, &
       tracersfn,                           &
       plumerise
  
  use plume_utils, only: prfrq
  use mem_scalar, only: RECYCLE_TRACERS ! CATT

  !TEB_SPM
  use teb_spm_start, only: TEB_SPM

  use mem_emiss, only : &
       ichemi,          &
       ichemi_in,       &
       chemdata_in,     &
       isource,         &
       weekdayin,       &
       efsat,           &
       efsun,           &
       eindno,          &
       eindno2,         &
       eindpm,          &
       eindco,          &
       eindso2,         &
       eindvoc,         &
       eveino,          &
       eveino2,         &
       eveipm,          &
       eveico,          &
       eveiso2,         &
       eveivoc

  use teb_vars_const, only : &
       rushh1,               &
       rushh2,               &
       daylight,             &
       iteb,                 &
       tminbld,              &
       nteb,                 &
       hc_roof,              &
       tc_roof,              &
       d_roof,               &
       hc_road,              &
       d_road,               &
       tc_road,              &
       d_wall,               &
       tc_wall,              &
       hc_wall,              &
       nurbtype,             &
       ileafcod,             &
       z0_town,              &
       bld,                  &
       bld_height,           &
       bld_hl_ratio,         &
       aroof,                &
       eroof,                &
       aroad,                &
       eroad,                &
       awall,                &
       ewall,                &
       htraf,                &
       hindu,                &
       pletraf,              &
       pleindu
  
  use mem_mass, only : iexev,imassflx

  use grid_dims, only: &
       maxsteb,        &
       maxubtp

  ! Explicit domain decomposition
  use domain_decomp, only: domain_fname
  
  ! Logical tests for microphysics complexity 
  use therm_lib, only: vapour_on,cloud_on,bulk_on,level

  implicit none

  character(len=*), intent(in) :: fileName  ! file name with namelists

  integer :: i                        ! loop count
  integer :: iunit                    ! io unit number
  integer, parameter :: firstUnit=20  ! lowest io unit number available
  integer, parameter :: lastUnit=99   ! highest io unit number available
  logical :: op                       ! io unit number opened or not
  logical :: ex                       ! namelist file exists?
  integer :: err                      ! return code on iostat
  character(len=10) :: c0             ! scratch
  character(len=*), parameter :: h="**(ReadNamelist)**"  ! program unit name

  namelist /MODEL_GRIDS/                                               &
       expnme, runtype, load_bal, imontha, idatea,                     &
       iyeara, itimea, imonthz, idatez,                                &
       iyearz, itimez, ngrids, nnxp, nnyp, nnzp, nzg, nzs, nxtnest,    &
       domain_fname,                                                   &
       if_adap, ihtran, deltax, deltay, deltaz, dzrat, dzmax, zz,      &
       dtlong, nacoust, ideltat, nstratx, nstraty, nndtrat, nestz1,    &
       nstratz1, nestz2, nstratz2, polelat, polelon, centlat, centlon, &
       ninest, njnest, nknest, nnsttop, nnstbot, gridu, gridv

  namelist /CATT_INFO/                                                 &
       catt,                                                           &
       firemapfn, recycle_tracers,                                     &
       plumerise, prfrq

  namelist /TEB_SPM_INFO/                                              &
       teb_spm,                                                        &
       fusfiles, ifusflg, ifusfn,                                      &
       ichemi, ichemi_in, chemdata_in, isource, weekdayin, rushh1,     &
       rushh2, daylight, efsat, efsun, eindno, eindno2, eindpm,        &
       eindco, eindso2, eindvoc, eveino, eveino2, eveipm, eveico,      &
       eveiso2, eveivoc, iteb, tminbld, nteb, hc_roof, tc_roof,        &
       d_roof, hc_road, tc_road, d_road, hc_wall, tc_wall, d_wall,     &
       nurbtype, ileafcod, z0_town, bld, bld_height, bld_hl_ratio,     &
       aroof, eroof, aroad, eroad, awall, ewall, htraf, hindu,         &
       pletraf, pleindu

  namelist /MODEL_FILE_INFO/                                           &
       initial, nud_type, varfpfx, vwait1, vwaittot, nud_hfile, nudlat,&
       tnudlat, tnudcent, tnudtop, znudtop, wt_nudge_grid, wt_nudge_uv,&
       wt_nudge_th, wt_nudge_pi, wt_nudge_rt, wt_nudge_co2, nud_cond,  &
       cond_hfile,tcond_beg, tcond_end, t_nudge_rc, wt_nudgec_grid,    &
       if_oda,oda_upaprefix,oda_sfcprefix, frqoda, todabeg, todaend,   &
       tnudoda, wt_oda_grid, wt_oda_uv, wt_oda_th, wt_oda_pi,          &
       wt_oda_rt, roda_sfce, roda_sfc0, roda_upae,roda_upa0, roda_hgt, &
       roda_zfact, oda_sfc_til, oda_sfc_tel, oda_upa_til, oda_upa_tel, &
       if_cuinv, cu_prefix, tnudcu, wt_cu_grid, tcu_beg, tcu_end,      &
       cu_tel, cu_til, imonthh, idateh, iyearh, itimeh,                &
       hfilin, ipastin, pastfn, ioutput,                               &
       hfilout, afilout, iclobber, ihistdel, frqhis, frqanl, frqlite,  &
       xlite, ylite, zlite, nlite_vars, lite_vars, avgtim, frqmean,    &
       frqboth, kwrite, frqprt, initfld, topfiles, sfcfiles, sstfpfx,  &
       ndvifpfx, itoptflg, isstflg, ivegtflg, isoilflg, ndviflg,       &
       nofilflg, iupdndvi, iupdsst, itoptfn, isstfn, ivegtfn, isoilfn, &
       ndvifn, itopsflg, toptenh, toptwvl, iz0flg, z0max, z0fact,      &
       mkcoltab, coltabfn

  namelist /CUPARM_OPTIONS/ &
       nnqparm,nclouds,ndeepest,nshallowest,wcldbs,confrq,cptime,        &
       iupmethod,iupstrm,radius,depth_min,cap_maxs,zkbmax,zcutdown,      &
       z_detr,max_heat,closure_type,maxens_lsf,maxens_eff,maxens_cap

  namelist /MODEL_OPTIONS/ &
       naddsc, icorflg, iexev,imassflx, ibnd, jbnd, cphas, lsflg, nfpt,  &
       distim,iswrtyp, ilwrtyp,icumfdbk,                                 &
       raddatfn,radfrq, lonrad, npatch, nvegpat, isfcl,ico2,co2con       &
       nvgcon, pctlcon, nslcon, drtcon, zrough, albedo, seatmp, dthcon,  &
       soil_moist, soil_moist_fail, usdata_in, usmodel_in, slz, slmstr,  &
       stgoff, if_urban_canopy, idiffk, ibruvais, ibotflx, ihorgrad,     &
       csx, csz, xkhkm, zkhkm, akmin, akmax, hgtmin, hgtmax, level,      &
       icloud, irain, ipris, isnow, iaggr, igraup, ihail, cparm, rparm,  &
       pparm, sparm, aparm, gparm, hparm, gnu

  namelist /MODEL_SOUND/ &
       ipsflg, itsflg, irtsflg, iusflg, hs, ps, ts, rts, us, vs, co2s

  namelist /MODEL_PRINT/ &
       nplt, iplfld, ixsctn, isbval

  namelist /ISAN_CONTROL/ &
       iszstage, ivrstage, isan_inc, guess1st, i1st_flg, iupa_flg,       &
       isfc_flg, iapr, iarawi, iasrfce, varpfx, ioflgisz, ioflgvar

  namelist /ISAN_ISENTROPIC/ &
       nisn, levth, nigrids, topsigz, hybbot, hybtop, sfcinf, sigzwt,    &
       nfeedvar, maxsta, maxsfc, notsta, notid, iobswin, stasep, igridfl,&
       gridwt, gobsep, gobrad, wvlnth, swvlnth, respon


  ! defaults (just for arrays); should be revised to accomodate
  ! precise defaults

  ihoura = 0
  nnxp =0
  nnyp = 0
  nnzp = 0
  nxtnest = 0
  zz = 0.0
  nstratx = 0
  nstraty = 0
  nndtrat = 0
  nstratz1 = 0
  nstratz2 = 0
  centlat = 0.0
  centlon = 0.0
  ninest = 0
  njnest = 0
  nknest = 0
  nnsttop = 0
  nnstbot = 0
  gridu = 0.0
  gridv = 0.0
  wt_nudge_grid = 0.0
  wt_nudgec_grid = 0.0
  wt_oda_grid = 0.0
  roda_sfce = 0.0
  roda_sfc0 = 0.0
  roda_upae = 0.0
  roda_upa0 = 0.0
  roda_hgt  = 0.0
  roda_zfact = 0.0
  wt_cu_grid = 0.0
  lite_vars = " "
  itoptflg=0
  isstflg=0
  ivegtflg=0
  isoilflg=0
  ndviflg=0
  nofilflg=0
  itoptfn=" "
  isstfn=" "
  ivegtfn=" "
  isoilfn=" "
  ndvifn=" "
  itopsflg=0
  toptenh=0.0
  toptwvl=0.0
  iz0flg=0
  z0max=0.0
  gnu=0.0
  iplfld=" "
  ixsctn=0
  isbval=0
  us=0.0
  vs=0.0
  ts=0.0
  ps=0.0
  hs=0.0
  rts=0.0
  co2s=0.0
  co2_init=0.0
  slz=0.0
  slmstr=0.0
  stgoff=0.0
  idiffk=0
  ibotflx = 0
  ibruvais = 0
  csx=0.0
  csz=0.0
  xkhkm=0.0
  zkhkm=0.0
  akmin=0.0
  akmax=0.0
  hgtmin=0.0
  hgtmax=0.0
  levth=0
  notid=" "
  gridwt=0.0
  wvlnth=0.0
  swvlnth=0.0
  respon=0.0
  ! Initial Value
  domain_fname = ''
  ! CATT
  prfrq = 3600. ! Initial Value for PlumeRise Frequency - CATT
  recycle_tracers = 0
  ! TEB
  fusfiles    = ''
  ifusflg     = 0
  ichemi      = 0    !Photochemical module activation - 1=on, 0=off
  ichemi_in   = 0    !Use initial values from previous run (1=yes,0=no)
  chemdata_in = ''
  isource     = 0    !Emission module activation - 1=on, 0=off
  weekdayin = 'SUN'  !Initial weeakday of the simulation
  rushh1     = 7.81  !Morning Rush Hour (Local Time in Hours)
  rushh2     = 16.0  !Afternoon/Evening Rush Hour (Local Time)
  daylight   = 0.    !Daylight saving time (horario de verao)
  ! Emission factor (fraction of weekdays) for Saturdays and Sundays
  ! They are used in the emission module and TEB. - EDF
  efsat      = 0.8
  efsun      = 0.5
  ! Input GMT difference time variable (To define local time)
  !Industrial emissions (kg/s/m2)
  eindno     = 2.6636227e-10
  eindno2    = 2.9595805e-11
  eindpm     = 4.3421278e-10
  eindco     = 8.1599860e-10
  eindso2    = 3.6149164e-10
  eindvoc    = 2.5367833e-10 
  !Veicular emissions (kg/day/m2)
  eveino     = 4.3196708e-04
  eveino2    = 6.8566209e-05
  eveipm     = 6.2648396e-06
  eveico     = 7.5433785e-03
  eveiso2    = 4.0730592e-05
  eveivoc    = 1.1892237e-03
  !----- Urban canopy parameterization using TEB (Masson, 2000)-------------
  iteb       = 0     !1=on, 0=off
  tminbld    = 12.   !Minimum internal building temperature (degrees Celsius)
  nteb       = 3     !Number of roof,road and wall layers used in TEB, Max.3
  if (maxsteb==3) then
     ! ROOF layers properties
     ! heat capacity
     hc_roof(1) = 2110000.; hc_roof(2) = 280000.; hc_roof(3) = 290000.
     ! thermal conductivity
     tc_roof(1) = 0.41;     tc_roof(2) = 0.05;    tc_roof(3) = 0.03
     ! depth
     d_roof(1)  = 0.05;     d_roof(2)  = 0.4;     d_roof(3)  = 0.05
     ! ROAD layers properties
     ! heat capacity
     hc_road(1) = 1240000.; hc_road(2) = 1280000.; hc_road(3) = 1280000.
     ! thermal conductivity 1.01
     tc_road    = 1.0103
     ! depth
     d_road(1)  = 0.05;     d_road(2)  = 0.1;      d_road(3)  = 1.0
     ! WALL layers properties
     ! heat capacity J/m3/K 10e6
     hc_wall    = 1000000.
     ! thermal conductivity 0.81 W/m/K
     tc_wall    = 0.81
     ! depth
     d_wall(1)  = 0.02;     d_wall(2)  = 0.125;    d_wall(3)  = 0.02
  else
     ! roof layers properties
     hc_roof    = 0.    ! heat capacity
     tc_roof    = 0.    ! thermal conductivity
     d_roof     = 0.    ! depth
     ! road layers properties
     hc_road    = 0.    ! heat capacity
     tc_road    = 0.    ! thermal conductivity 1.0103
     d_road     = 0.    ! depth
     ! wall layers properties
     hc_wall    = 0.    ! heat capacity J/m3/K 10e6
     tc_wall    = 0.    ! thermal conductivity 0.81 W/m/K
     d_wall     = 0.    ! depth
  end if
  nurbtype   = 2         !Number of urban types (maximum of 3)
  if (maxubtp==2) then
     !Leaf class code to identify each urban type
     ileafcod(1) = 21;       ileafcod(2)     = 19
     !Urban type properties
     !Urban type roughness length 5 e 1
     z0_town(1) = 3.0;       z0_town(2)      = 0.5
     !Fraction occupied by buildings in the grid cell
     bld(1)    = 0.5;        bld(2)          = 0.7
     !Building Height
     bld_height(1) = 50.;    bld_height(2)   = 5.0
     !Vertical/Horizontal rate 3 e 0.5
     bld_hl_ratio(1) = 4.4;  bld_hl_ratio(2) = 2.4
     !Roof albedo
     aroof = 0.15
     !Roof emissivitiy
     eroof = 0.9
     !Road albedo
     aroad = 0.1
     !Road emissivity 90% masson
     eroad = 0.9
     !Wall albedo
     awall = 0.25
     !Wall emissivity
     ewall = 0.85
     !Maximum value of sensible heat
     htraf(1) = 90.0;        htraf(2) = 60.0
     !Maximum value of sensible heat
     hindu(1) = 10.0;        hindu(2) = 14.0
     !released by Industry (W/m2)
     !Maximum value of latent heat
     pletraf(1) = 10.0;      pletraf(2) = 5.0
     !released by Traffic (W/m2)
     !Maximum value of latent heat
     pleindu(1) = 30.0;      pleindu(1) = 50.0
  else
     !Urban type properties
     z0_town      = 0.0   !Urban type roughness length 5 e 1
     bld          = 0.0   !Fraction occupied by buildings in the grid cell
     bld_height   = 0.0   !Building Height
     bld_hl_ratio = 0.0   !Vertical/Horizontal rate 3 e 0.5
     aroof        = 0.0   !Roof albedo
     eroof        = 0.0   !Roof emissivitiy
     aroad        = 0.0   !Road albedo
     eroad        = 0.0   !Road emissivity 90% masson
     awall        = 0.0   !Wall albedo
     ewall        = 0.0   !Wall emissivity
     htraf        = 0.0   !Maximum value of sensible heat
     hindu        = 0.0   !Maximum value of sensible heat
     pletraf      = 0.0   !Maximum value of latent heat
     pleindu      = 0.0   !Maximum value of latent heat
  endif


  nnqparm=0
  ndeepest=0
  nshallowest=0
  confrq=0.
  cptime=0.
  radius=0.
  depth_min=0.
  cap_maxs=0.
  zkbmax=0.
  zcutdown=0.
  z_detr=0.
  max_heat=0.
  closure_type=""
  maxens_lsf=0
  maxens_eff=0
  maxens_cap=0

  ! select unused i/o unit

  do iunit = firstUnit, lastUnit
     inquire(iunit,opened=op)
     if (.not. op) exit
  end do

  if (iunit > lastUnit) then
     call abort_run('All i/o units are already in use.','ReadNamelist','rdint.f90')
  end if

  ! if namelist file exists, open, read each section and close

  inquire(file=trim(fileName), exist=ex)
  if (.not. ex) then
     call abort_run('Namelist file '//trim(fileName)//' does not exist.' &
                   ,'ReadNamelist','rdint.f90')
  end if

  open(iunit, file=trim(fileName), status="old", action="read",&
       iostat=err)
  if (err /= 0) then
     write(c0,"(i10)") err
     call abort_run(' Open namelist file '//trim(fileName)// &
                    ' returned iostat='//trim(adjustl(c0))//'...' &
                   ,'ReadNamelist','rdint.f90')
  end if

  read (iunit, iostat=err, NML=MODEL_GRIDS)
  if (err /= 0) then
     write(*,"(a)") "**(ERROR)** reading section MODEL_GRIDS "//&
          &"of namelist file "//trim(fileName)
     write(*,"(a)") " compare values read with file contents:"
     write(*,*) "expnme=",trim(expnme)
     write(*,*) "runtype=",trim(runtype)
     write(*,*) "load_bal=",load_bal
     write(*,*) "imontha=",imontha
     write(*,*) "idatea=",idatea
     write(*,*) "iyeara=",iyeara
     write(*,*) "itimea=",itimea
     write(*,*) "imonthz=",imonthz
     write(*,*) "idatez=",idatez
     write(*,*) "iyearz=",iyearz
     write(*,*) "itimez=",itimez
     write(*,*) "ngrids=",ngrids
     write(*,*) "nnxp=",nnxp
     write(*,*) "nnyp=",nnyp
     write(*,*) "nnzp=",nnzp
     write(*,*) "nzg=",nzg
     write(*,*) "nzs=",nzs
     write(*,*) "nxtnest=",nxtnest
     write(*,*) "DOMAIN_FNAME=", DOMAIN_FNAME
     write(*,*) "if_adap=",if_adap
     write(*,*) "ihtran=",ihtran
     write(*,*) "deltax=",deltax
     write(*,*) "deltay=",deltay
     write(*,*) "deltaz=",deltaz
     write(*,*) "dzrat=",dzrat
     write(*,*) "dzmax=",dzmax
     write(*,*) "zz=",zz
     write(*,*) "dtlong=",dtlong
     write(*,*) "nacoust=",nacoust
     write(*,*) "ideltat=",ideltat
     write(*,*) "nstratx=",nstratx
     write(*,*) "nstraty=",nstraty
     write(*,*) "nndtrat=",nndtrat
     write(*,*) "nestz1=",nestz1
     write(*,*) "nstratz1=",nstratz1
     write(*,*) "nestz2=",nestz2
     write(*,*) "nstratz2=",nstratz2
     write(*,*) "polelat=",polelat
     write(*,*) "polelon=",polelon
     write(*,*) "centlat=",centlat
     write(*,*) "centlon=",centlon
     write(*,*) "ninest=",ninest
     write(*,*) "njnest=",njnest
     write(*,*) "nknest=",nknest
     write(*,*) "nnsttop=",nnsttop
     write(*,*) "nnstbot=",nnstbot
     write(*,*) "gridu=",gridu
     write(*,*) "gridv=",gridv
     call abort_run('Error reading namelist, MODEL_GRIDS block.' &
                   ,'ReadNamelist','rdint.f90')
  end if

  read (iunit, iostat=err, NML=CATT_INFO)
  if (err /= 0) then
     write(*,"(a)") "**(ERROR)** reading section CATT_INFO "//&
          &"of namelist file "//trim(fileName)
     write(*,"(a)") " compare values read with file contents:"
     print *, "CATT=", CATT
     write(*,*) "FIREMAPFN=", FIREMAPFN
     !!write(*,*) "TRACERSFN=", TRACERSFN
     write(*,*) "RECYCLE_TRACERS=", RECYCLE_TRACERS
     write(*,*) "PLUMERISE=", PLUMERISE
     write(*,*) "PRFRQ=", PRFRQ
     call abort_run('Error reading namelist, CATT_INFO block.' &
                   ,'ReadNamelist','rdint.f90')
  end if

  read (iunit, iostat=err, NML=TEB_SPM_INFO)
  if (err /= 0) then
     write(*,"(a)") "**(ERROR)** reading section TEB_SPM_INFO "//&
          &"of namelist file "//trim(fileName)
     write(*,"(a)") " compare values read with file contents:"
     print *, "TEB_SPM=", TEB_SPM
     print *, "ifusflg=", ifusflg
     print *, "ifusfn=", ifusfn
     print *, "fusfiles=", trim(fusfiles)
     print *, "ICHEMI=", ICHEMI
     print *, "ICHEMI_IN=", ICHEMI_IN
     print *, "CHEMDATA_IN=", CHEMDATA_IN
     print *, "ISOURCE=", ISOURCE
     print *, "WEEKDAYIN=", trim(WEEKDAYIN)
     print *, "RUSHH1=", RUSHH1
     print *, "RUSHH2=", RUSHH2
     print *, "DAYLIGHT=", DAYLIGHT
     print *, "EFSAT=", EFSAT
     print *, "EFSUN=", EFSUN
     print *, "EINDNO=", EINDNO
     print *, "EINDNO2=", EINDNO2
     print *, "EINDPM=", EINDPM
     print *, "EINDCO=", EINDCO
     print *, "EINDSO2=", EINDSO2
     print *, "EINDVOC=", EINDVOC
     print *, "EVEINO=", EVEINO
     print *, "EVEINO2=", EVEINO2
     print *, "EVEIPM=", EVEIPM
     print *, "EVEICO=", EVEICO
     print *, "EVEISO2=", EVEISO2
     print *, "EVEIVOC=", EVEIVOC
     print *, "ITEB=", ITEB
     print *, "TMINBLD=", TMINBLD
     print *, "NTEB=", NTEB
     print *, "HC_ROOF=", HC_ROOF
     print *, "TC_ROOF=", TC_ROOF
     print *, "D_ROOF=", D_ROOF
     print *, "HC_ROAD=", HC_ROAD
     print *, "TC_ROAD=", TC_ROAD
     print *, "D_ROAD=", D_ROAD
     print *, "HC_WALL=", HC_WALL
     print *, "TC_WALL=", TC_WALL
     print *, "D_WALL=", D_WALL
     print *, "NURBTYPE=", NURBTYPE
     print *, "ILEAFCOD=", ILEAFCOD
     print *, "Z0_TOWN=", Z0_TOWN
     print *, "BLD=", BLD
     print *, "BLD_HEIGHT=", BLD_HEIGHT
     print *, "BLD_HL_RATIO=", BLD_HL_RATIO
     print *, "AROOF=", AROOF
     print *, "EROOF=", EROOF
     print *, "AROAD=", AROAD
     print *, "EROAD=", EROAD
     print *, "AWALL=", AWALL
     print *, "EWALL=", EWALL
     print *, "HTRAF=", HTRAF
     print *, "HINDU=", HINDU
     print *, "PLETRAF=", PLETRAF
     print *, "PLEINDU=", PLEINDU
     call abort_run('Error reading namelist, TEB_SPM_INFO block.' &
                   ,'ReadNamelist','rdint.f90')
  end if

  read (iunit, iostat=err, NML=MODEL_FILE_INFO)
  if (err /= 0) then
     write(*,"(a)") "**(ERROR)** reading section MODEL_FILE_INFO "//&
          &"of namelist file "//trim(fileName)
     write(*,"(a)") " compare values read with file contents:"
     write (*, "(a)") " namelist MODEL_FILE_INFO: "
     write (*,*) "initial=", initial
     write (*,*) "nud_type=", nud_type
     write (*,*) "varfpfx=", trim(varfpfx)
     write (*,*) "vwait1=", vwait1
     write (*,*) "vwaittot=", vwaittot
     write (*,*) "nud_hfile=", trim(nud_hfile)
     write (*,*) "nudlat=", nudlat
     write (*,*) "tnudlat=", tnudlat
     write (*,*) "tnudcent=", tnudcent
     write (*,*) "tnudtop=", tnudtop
     write (*,*) "znudtop=", znudtop
     write (*,*) "wt_nudge_grid=", wt_nudge_grid
     write (*,*) "wt_nudge_uv=", wt_nudge_uv
     write (*,*) "wt_nudge_th=", wt_nudge_th
     write (*,*) "wt_nudge_pi=", wt_nudge_pi
     write (*,*) "wt_nudge_rt=", wt_nudge_rt
     write (*,*) "wt_nudge_co2=", wt_nudge_co2
     write (*,*) "nud_cond=", nud_cond
     write (*,*) "cond_hfile=", trim(cond_hfile)
     write (*,*) "tcond_beg=", tcond_beg
     write (*,*) "tcond_end=", tcond_end
     write (*,*) "t_nudge_rc=", t_nudge_rc
     write (*,*) "wt_nudgec_grid=", wt_nudgec_grid
     write (*,*) "if_oda=", if_oda
     write (*,*) "oda_upaprefix=", trim(oda_upaprefix)
     write (*,*) "oda_sfcprefix=", trim(oda_sfcprefix)
     write (*,*) "frqoda=", frqoda
     write (*,*) "todabeg=", todabeg
     write (*,*) "todaend=", todaend
     write (*,*) "tnudoda=", tnudoda
     write (*,*) "wt_oda_grid=", wt_oda_grid
     write (*,*) "wt_oda_uv=", wt_oda_uv
     write (*,*) "wt_oda_th=", wt_oda_th
     write (*,*) "wt_oda_pi=", wt_oda_pi
     write (*,*) "wt_oda_rt=", wt_oda_rt
     write (*,*) "roda_sfce=", roda_sfce
     write (*,*) "roda_sfc0=", roda_sfc0
     write (*,*) "roda_upae=", roda_upae
     write (*,*) "roda_upa0=", roda_upa0
     write (*,*) "roda_hgt=", roda_hgt
     write (*,*) "roda_zfact=", roda_zfact
     write (*,*) "oda_sfc_til=", oda_sfc_til
     write (*,*) "oda_sfc_tel=", oda_sfc_tel
     write (*,*) "oda_upa_til=", oda_upa_til
     write (*,*) "oda_upa_tel=", oda_upa_tel
     write (*,*) "if_cuinv=", if_cuinv
     write (*,*) "cu_prefix=", trim(cu_prefix)
     write (*,*) "tnudcu=", tnudcu
     write (*,*) "wt_cu_grid=", wt_cu_grid
     write (*,*) "tcu_beg=", tcu_beg
     write (*,*) "tcu_end=", tcu_end
     write (*,*) "cu_tel=", cu_tel
     write (*,*) "cu_til=", cu_til
     write (*,*) "imonth=", imonthh          !
     write (*,*) "idateh=", idateh           !
     write (*,*) "iyearh=", iyearh           !
     write (*,*) "itimeh=", itimeh           !
     write (*,*) "hfilin=", trim(hfilin)
     write (*,*) "ipastin=", ipastin
     write (*,*) "pastfn=", trim(pastfn)
     write (*,*) "ioutput=", ioutput
     write (*,*) "hfilout=", trim(hfilout)
     write (*,*) "afilout=", trim(afilout)
     write (*,*) "iclobber=", iclobber
     write (*,*) "ihistdel=", ihistdel
     write (*,*) "frqhis=", frqhis
     write (*,*) "frqanl=", frqanl
     write (*,*) "frqlite=", frqlite
     write (*,*) "xlite=", xlite
     write (*,*) "ylite=", ylite
     write (*,*) "zlite=", zlite
     write (*,*) "nlite_vars=", nlite_vars
     write (*,*) "lite_vars=", (trim(lite_vars(i))//";", i=1,size(lite_vars))
     write (*,*) "avgtim=", avgtim
     write (*,*) "frqmean=", frqmean
     write (*,*) "frqboth=", frqboth
     write (*,*) "kwrite=", kwrite
     write (*,*) "frqprt=", frqprt
     write (*,*) "initfld=", initfld
     write (*,*) "topfiles=", trim(topfiles)
     write (*,*) "sfcfiles=", trim(sfcfiles)
     write (*,*) "sstfpfx=", trim(sstfpfx)
     write (*,*) "ndvifpfx=", trim(ndvifpfx)
     write (*,*) "itoptflg=", itoptflg
     write (*,*) "isstflg=", isstflg
     write (*,*) "ivegtflg=", ivegtflg
     write (*,*) "isoilflg=", isoilflg
     write (*,*) "ndviflg=", ndviflg
     write (*,*) "nofilflg=", nofilflg
     write (*,*) "iupdndvi=", iupdndvi
     write (*,*) "iupdsst=", iupdsst
     write (*,*) "itoptfn=", (trim(itoptfn(i))//";", i =1,size(itoptfn))
     write (*,*) "isstfn=", (trim(isstfn(i))//";", i=1,size(isstfn))
     write (*,*) "ivegtfn=", (trim(ivegtfn(i))//";", i = 1, size(ivegtfn))
     write (*,*) "isoilfn=", (trim(isoilfn(i))//";", i = 1, size(isoilfn))
     write (*,*) "ndvifn=", (trim(ndvifn(i))//";", i=1,size(ndvifn))
     write (*,*) "itopsflg=", itopsflg
     write (*,*) "toptenh=", toptenh
     write (*,*) "toptwvl=", toptwvl
     write (*,*) "iz0flg=", iz0flg
     write (*,*) "z0max=", z0max
     write (*,*) "z0fact=", z0fact
     write (*,*) "mkcoltab=", mkcoltab
     write (*,*) "coltabfn=", trim(coltabfn)
     call abort_run('Error reading namelist, MODEL_FILE_INFO block.' &
                   ,'ReadNamelist','rdint.f90')
  end if

  read (iunit, iostat=err, NML=CUPARM_OPTIONS)
  if (err /= 0) then
     write(*,"(a)") "**(ERROR)** reading section CUPARM_OPTIONS "//&
          &"of namelist file "//trim(fileName)
     write (*, *) "nnqparm=",nnqparm
     write (*, *) "nclouds=",nclouds
     write (*, *) "ndeepest=",ndeepest
     write (*, *) "nshallowest=",nshallowest
     write (*, *) "wcldbs=",wcldbs
     write (*, *) "confrq=",confrq
     write (*, *) "cptime=",cptime
     write (*, *) "iupmethod=",iupmethod
     write (*, *) "iupstrm=",iupstrm
     write (*, *) "radius=",radius
     write (*, *) "depth_min=",depth_min
     write (*, *) "cap_maxs=",depth_min
     write (*, *) "zkbmax=",zkbmax
     write (*, *) "zcutdown=",zcutdown
     write (*, *) "z_detr=",z_detr
     write (*, *) "max_heat=",max_heat
     write (*, *) "closure_type=",closure_type
     write (*, *) "maxens_lsf=",maxens_lsf
     write (*, *) "maxens_eff=",maxens_eff
     write (*, *) "maxens_cap=",maxens_cap

     call abort_run('Error reading namelist, CUPARM_OPTIONS block.' &
                   ,'ReadNamelist','rdint.f90')
  end if

  read (iunit, iostat=err, NML=MODEL_OPTIONS)
  if (err /= 0) then
     write(*,"(a)") "**(ERROR)** reading section MODEL_OPTIONS "//&
          &"of namelist file "//trim(fileName)
     write(*,"(a)") " compare values read with file contents:"
     write (*,*) "naddsc=",naddsc
     write (*, *) "icorflg=",icorflg
     ![MLO
     write (*, *) "iexev=",iexev
     write (*, *) "imassflx=",imassflx
     !MLO]
     write (*, *) "ibnd=",ibnd
     write (*, *) "jbnd=",jbnd
     write (*, *) "cphas=",cphas
     write (*, *) "lsflg=",lsflg
     write (*, *) "nfpt=",nfpt
     write (*, *) "distim=",distim
     write (*, *) "iswrtyp=",iswrtyp
     write (*, *) "ilwrtyp=",ilwrtyp
     write (*, *) "raddatfn=", RADDATFN
     write (*, *) "radfrq=",radfrq
     write (*, *) "lonrad=",lonrad
     write (*, *) "npatch=",npatch
     write (*, *) "nvegpat=",nvegpat
     write (*, *) "isfcl=",isfcl
     write (*, *) "ico2=",ico2
     write (*, *) "co2con=",co2con
     write (*, *) "nvgcon=",nvgcon
     write (*, *) "pctlcon=",pctlcon
     write (*, *) "nslcon=",nslcon
     write (*, *) "drtcon=",drtcon
     write (*, *) "zrough=",zrough
     write (*, *) "albedo=",albedo
     write (*, *) "seatmp=",seatmp
     write (*, *) "dthcon=",dthcon
     write (*, *) "soil_moist=",soil_moist
     write (*, *) "soil_moist_fail=",soil_moist_fail
     write (*, *) "usdata_in=",trim(usdata_in)
     write (*, *) "usmodel_in=",trim(usmodel_in)
     write (*, *) "slz=",slz
     write (*, *) "slmstr=",slmstr
     write (*, *) "stgoff=",stgoff
     write (*, *) "if_urban_canopy=",if_urban_canopy
     write (*, *) "idiffk=",idiffk
     write (*, *) "ibruvais=",ibruvais
     write (*, *) "ibotflx=",ibotflx
     write (*, *) "ihorgrad=",ihorgrad
     write (*, *) "csx=",csx
     write (*, *) "csz=",csz
     write (*, *) "xkhkm=",xkhkm
     write (*, *) "zkhkm=",zkhkm
     write (*, *) "akmin=",akmin
     write (*, *) "akmax=",akmax
     write (*, *) "hgtmin=",hgtmin
     write (*, *) "hgtmax=",hgtmax
     write (*, *) "level=",level
     write (*, *) "icloud=",icloud
     write (*, *) "irain=",irain
     write (*, *) "ipris=",ipris
     write (*, *) "isnow=",isnow
     write (*, *) "iaggr=",iaggr
     write (*, *) "igraup=",igraup
     write (*, *) "ihail=",ihail
     write (*, *) "cparm=",cparm
     write (*, *) "rparm=",rparm
     write (*, *) "pparm=",pparm
     write (*, *) "sparm=",sparm
     write (*, *) "aparm=",aparm
     write (*, *) "gparm=",gparm
     write (*, *) "hparm=",hparm
     write (*, *) "gnu=",gnu
     call abort_run('Error reading namelist, MODEL_OPTIONS block.' &
                   ,'ReadNamelist','rdint.f90')
  else
     !----- Switching closure type to lower case
     call tolower(closure_type,maxclouds)
  end if

  write (unit=*,fmt='(a)') 'Reading ED2 namelist information'
  call read_ednl(iunit)

  read (iunit, iostat=err, NML=MODEL_SOUND)
  if (err /= 0) then
     write(*,"(a)") "**(ERROR)** reading section MODEL_SOUND "//&
          &"of namelist file "//trim(fileName)
     write(*,"(a)") " compare values read with file contents:"
     write (*, *) "ipsflg=",ipsflg
     write (*, *) "itsflg=",itsflg
     write (*, *) "irtsflg=",irtsflg
     write (*, *) "iusflg=",iusflg
     write (*, *) "hs=",hs
     write (*, *) "ps=",ps
     write (*, *) "ts=",ts
     write (*, *) "rts=",rts
     write (*, *) "us=",us
     write (*, *) "vs=",vs
     write (*, *) "co2s=",co2s
     call abort_run('Error reading namelist, MODEL_SOUND block.' &
                   ,'ReadNamelist','rdint.f90')
  end if

  read (iunit, iostat=err, NML=MODEL_PRINT)
  if (err /= 0) then
     write(*,"(a)") "**(ERROR)** reading section MODEL_PRINT "//&
          &"of namelist file "//trim(fileName)
     write(*,"(a)") " compare values read with file contents:"
     write (*, *) "nplt=",nplt
     write (*, *) "iplfld=",(trim(iplfld(i))//";", i=1,size(iplfld))
     write (*, *) "ixsctn=",ixsctn
     write (*, *) "isbval=",isbval
     call abort_run('Error reading namelist, MODEL_PRINT block.' &
                   ,'ReadNamelist','rdint.f90')
  end if


  read (iunit, iostat=err, NML=ISAN_CONTROL)
  if (err /= 0) then
     write(*,"(a)") "**(ERROR)** reading section ISAN_CONTROL "//&
          &"of namelist file "//trim(fileName)
     write(*,"(a)") " compare values read with file contents:"
     write (*, *) "iszstage=",iszstage
     write (*, *) "ivrstage=",ivrstage
     write (*, *) "isan_inc=",isan_inc
     write (*, *) "guess1st=",guess1st
     write (*, *) "i1st_flg=",i1st_flg
     write (*, *) "iupa_flg=",iupa_flg
     write (*, *) "isfc_flg=",isfc_flg
     write (*, *) "iapr=",trim(iapr)
     write (*, *) "iarawi=",trim(iarawi)
     write (*, *) "iasrfce=",trim(iasrfce)
     write (*, *) "varpfx=",trim(varpfx)
     write (*, *) "ioflgisz=",ioflgisz
     write (*, *) "ioflgvar=",ioflgvar
     call abort_run('Error reading namelist, ISAN_CONTROL block.' &
                   ,'ReadNamelist','rdint.f90')
  end if

  read (iunit, iostat=err, NML=ISAN_ISENTROPIC)
  if (err /= 0) then
     write(*,"(a)") "**(ERROR)** reading section ISAN_ISENTROPIC "//&
          &"of namelist file "//trim(fileName)
     write(*,"(a)") " compare values read with file contents:"
     write (*, *) "nisn=",nisn
     write (*, *) "levth=",levth
     write (*, *) "nigrids=",nigrids
     write (*, *) "topsigz=",topsigz
     write (*, *) "hybbot=",hybbot
     write (*, *) "hybtop=",hybtop
     write (*, *) "sfcinf=",sfcinf
     write (*, *) "sigzwt=",sigzwt
     write (*, *) "nfeedvar=",nfeedvar
     write (*, *) "maxsta=",maxsta
     write (*, *) "maxsfc=",maxsfc
     write (*, *) "notsta=",notsta
     write (*, *) "notid=",(trim(notid(i))//";", i=1,size(notid))
     write (*, *) "iobswin=",iobswin
     write (*, *) "stasep=",stasep
     write (*, *) "igridfl=",igridfl
     write (*, *) "gridwt=",gridwt
     write (*, *) "gobsep=",gobsep
     write (*, *) "gobrad=",gobrad
     write (*, *) "wvlnth=",wvlnth
     write (*, *) "swvlnth=",swvlnth
     write (*, *) "respon=",respon
     call abort_run('Error reading namelist, ISAN_ISENTROPIC block.' &
                   ,'ReadNamelist','rdint.f90')
  end if

  close(iunit, iostat=err)
  if (err /= 0) then
     write(c0,"(i10)") err
     call abort_run('Closing file '//trim(fileName)//' returned iostat='//trim(adjustl(c0)) &
                   ,'ReadNamelist','rdint.f90')
  end if
  
  call date_2_seconds (iyearz,imonthz,idatez,itimez*100, &
  iyeara,imontha,idatea,itimea*100,timmax)
  
  call date_2_seconds (iyearh,imonthh,idateh,itimeh*100, &
  iyeara,imontha,idatea,itimea*100,timstr)

  !---- If this is a coupled run, make npatch = 2 and nvegpat= 1 --------------------------!
  if (isfcl == 5) then
     npatch = 2
     nvegpat = 1
  else !---- Not an ED-BRAMS run, and isoilflg/ivegtflg are set to 3, switch them by 1. ---!
     where (isoilflg == 3) isoilflg = 1
     where (ivegtflg == 3) ivegtflg = 1
  end if
  !---- If someone accidentally defined itoptflg, isstflg or ndviflg to 3, make them 1. ---!
  where (itoptflg == 3) itoptflg = 1
  where (isstflg  == 3) isstflg  = 1
  where (ndviflg  == 3) ndviflg  = 1

  !----------------------------------------------------------------------------------------!
  !    Saving the moisture complexity level into logical variables. Note that vapour_on is !
  ! not level == 1, it will be true when level is 1, 2, and 3. (whenever vapour is on).    !
  ! Likewise, cloud_on will be true when level is either 2 or 3. Bulk microphysics will be !
  ! true only when level >= 3 (levels = 4 and 5 exist too but rarely used).                !
  !----------------------------------------------------------------------------------------!
  vapour_on = level >= 1
  cloud_on  = level >= 2
  bulk_on   = level >= 3

  !----------------------------------------------------------------------------------------!
  !    Saving the CO2 complexity level into a logical variable.                            !
  !----------------------------------------------------------------------------------------!
  co2_on    = ico2 > 0

  return
end subroutine ReadNamelist
