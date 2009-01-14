
subroutine read_ednl(iunit)
  
  use max_dims , only: n_pft

  use soil_coms, only: ed_zrough => zrough, soil_database, &
       isoilstateinit, isoildepthflg, isoilbc, soilstate_db, soildepth_db,   &
       runoff_time,veg_database

  use met_driver_coms,only: ed_met_driver_db,imettype,metcyc1,metcycf,initial_co2 &
                           ,lapse_scheme

  use mem_sites, only: n_soi, soi_lat, soi_lon, n_ed_region, ed_reg_latmin,  &
       ed_reg_latmax, ed_reg_lonmin, ed_reg_lonmax, grid_res, grid_type, edres, &
       maxpatch, maxcohort
  
  use physiology_coms, only: istoma_scheme, n_plant_lim

  use phenology_coms, only: iphen_scheme,repro_scheme,iphenys1,iphenysf,iphenyf1,iphenyff,phenpath
 
  use decomp_coms, only: n_decomp_lim
  
  use disturb_coms, only: include_fire, ianth_disturb,   &
       treefall_disturbance_rate
  
  use pft_coms, only: include_these_pft,pft_1st_check
  
  use misc_coms, only:  ifoutput,idoutput,imoutput,iyoutput,isoutput, &
       frqfast, frqstate, outfast,outstate,unitfast,unitstate,        &
       ied_init_mode, current_time,ed_inputs_dir,                     &
       end_time, integration_scheme, ffilout,  dtlsm,                 &
       iprintpolys,printvars,npvars,pfmtstr,ipmax,ipmin,              &
       iedcnfgf,ffilout,sfilout,sfilin,event_file

  use grid_coms, only: timmax,time
  
  use ed_misc_coms,only: attach_metadata
  
  use optimiz_coms, only : ioptinpt


  !====== BRAMS modules ======!
  
  use mem_grid,only:expnme,runtype,itimez, idatez, imonthz, iyearz,  &
       itimea, idatea, imontha, iyeara,centlon,centlat,deltax,deltay,nnxp,nnyp,nstratx, &
       nstraty,polelat,polelon,ngrids,nzg, nzs,npatch
  
  use io_params,only: ioutput, iclobber, frqanl,frqhis,iyearh,idateh,itimeh,imonthh,isoilflg
  use mem_leaf,only:nslcon,leaf_zrough => zrough,slz, stgoff,slmstr,isfcl,nvegpat
  use mem_radiate, only: radfrq

  !============================!
  
  implicit none

  integer :: iunit                    ! io unit number
  integer :: err

  logical :: fexists,op
  namelist /ED2_INFO/  dtlsm,ifoutput,idoutput,imoutput,iyoutput,isoutput,   &
       attach_metadata,outfast,outstate,ffilout,sfilout,ied_init_mode,       &
       edres,sfilin,veg_database,soil_database,ed_inputs_dir,soilstate_db,   &
       soildepth_db,isoilstateinit,isoildepthflg,isoilbc,integration_scheme, &
       istoma_scheme,iphen_scheme,repro_scheme,lapse_scheme,n_plant_lim,     &
       n_decomp_lim,include_fire,ianth_disturb,include_these_pft,            &
       pft_1st_check,maxpatch,maxcohort,treefall_disturbance_rate,           &
       runoff_time,iprintpolys,npvars,printvars,pfmtstr,ipmin,ipmax,         &
       initial_co2,iphenys1,iphenysf,iphenyf1,iphenyff,iedcnfgf,event_file,  &
       phenpath

  read (iunit, iostat=err, NML=ED2_INFO)
  if (err /= 0) then
     write(*,"(a)") "**(ERROR)** reading section ED2_INFO "//&
          &"of namelist file "
     write(*,"(a)") " compare values read with file contents:"
     write(*,*) "dtlsm=",dtlsm
     write(*,*) "ifoutput=",ifoutput
     write(*,*) "idoutput=",idoutput
     write(*,*) "imoutput=",imoutput
     write(*,*) "iyoutput=",iyoutput
     write(*,*) "isoutput=",isoutput
     write(*,*) "attach_metadata=",attach_metadata
     write(*,*) "outfast=",outfast
     write(*,*) "outstate=",outstate
     write(*,*) "ffilout=",ffilout
     write(*,*) "sfilout=",sfilout
     write(*,*) "ied_init_mode=",ied_init_mode
     write(*,*) "edres=",edres
     write(*,*) "sfilin=",sfilin
     write(*,*) "veg_database=",veg_database
     write(*,*) "soil_database=",soil_database
     write(*,*) "ed_inputs_dir=",ed_inputs_dir
     write(*,*) "soilstate_db=",soilstate_db
     write(*,*) "soildepth_db=",soildepth_db
     write(*,*) "isoilstateinit=",isoilstateinit
     write(*,*) "isoildepthflg=",isoildepthflg
     write(*,*) "isoilbc=",isoilbc
     write(*,*) "integration_scheme=",integration_scheme
     write(*,*) "istoma_scheme=",istoma_scheme
     write(*,*) "iphen_scheme",iphen_scheme
     write(*,*) "repro_scheme",repro_scheme
     write(*,*) "lapse_scheme",lapse_scheme
     write(*,*) "n_plant_lim=",n_plant_lim
     write(*,*) "n_decomp_lim=",n_decomp_lim
     write(*,*) "include_fire=",include_fire
     write(*,*) "ianth_disturb=",ianth_disturb
     write(*,*) "include_these_pft=",include_these_pft
     write(*,*) "pft_1st_check=",pft_1st_check
     write(*,*) "maxpatch=",maxpatch
     write(*,*) "maxcohort=",maxcohort
     write(*,*) "treefall_disturbance_rate=",treefall_disturbance_rate
     write(*,*) "runoff_time=",runoff_time
     write(*,*) "iprintpolys=",iprintpolys
     write(*,*) "npvars=",npvars
     write(*,*) "printvars=",printvars
     write(*,*) "pfmtstr=",pfmtstr
     write(*,*) "ipmin=",ipmin
     write(*,*) "ipmax=",ipmax
     write(*,*) "initial_co2=",initial_co2
     write(*,*) "iphenys1=",iphenys1
     write(*,*) "iphenysf=",iphenysf
     write(*,*) "iphenyf1=",iphenyf1
     write(*,*) "iphenyff=",iphenyff
     write(*,*) "iedcnfgf=",iedcnfgf
     write(*,*) "event_file=",event_file
     write(*,*) "phenpath =",phenpath 
     call abort_run('Error reading namelist, ED2_INFO block.' &
          ,'read_ednl','ed_load_namelist.f90')
  end if

  call copy_in_bramsnl(expnme, runtype, itimez, idatez, imonthz, iyearz, &
       itimea, idatea, imontha, iyeara, radfrq, nnxp,nnyp,deltax,        &
       deltay,polelat,polelon,centlat,centlon,nstratx,nstraty,           &
       iclobber,nzg,nzs,isoilflg,nslcon,slz,slmstr,stgoff,leaf_zrough,ngrids)
  
  ! The following are standard namelist variables from ED2 in stand-alone
  ! for a coupled run, these can be fixed to any old value
  ! ---------------------------------------------------------------------
  
  
  n_ed_region   = ngrids   ! SAME THING
  grid_res      = 0        ! NOT USED - THIS IS FOR LAT-LON STYLE
  grid_type     = 1        ! This must be polar-stereo for now
  ed_reg_latmin = 0        ! NOT USED IN COUPLED
  ed_reg_latmax = 0        ! NOT USED IN COUPLED
  ed_reg_lonmin = 0        ! NOT USED IN COUPLED
  ed_reg_lonmax = 0        ! NOT USED IN COUPLED
  n_soi = 0                ! NO SOI IN COUPLED RUN
  soi_lat = 0              ! NO SOI IN COUPLED RUN
  soi_lon = 0              ! NO SOI IN COUPLED RUN
  ed_met_driver_db = ''    ! NOT USED IN COUPLED
  imettype = 1             ! NOT USED IN COUPLED
  metcyc1  = 0000          ! NOT USED IN COUPLED
  metcycf  = 0000          ! NOT USED IN COUPLED
  ioptinpt = ''            !NOT USED, DEPRICATED? It will be used once optimization is 
                           ! implemented in the new version.
  
  unitfast  = 0           ! Since BRAMS uses frqanl and frqhist in seconds, there is no
  unitstate = 0           ! reason to ask the user for units for outfast and outstate, the
                          ! special flags cover all possibilities.

  
  ! Force LEAF3's number of patches to be 2. This is necessary for
  ! filling the lower boundary condition arrays correctly
  ! --------------------------------------------------------------------
  if (isfcl==5) then
     npatch   = 2
     nvegpat  = 1
  endif


  ! These are ED2 variables that have the same function as certain BRAMS
  ! namelist variables under a different name
  ! --------------------------------------------------------------------
  
  frqfast  = frqanl
  frqstate = frqhis

  ! set current time to initial time here.  If this is a history run,
     ! reset current time in subroutine history_start.
     end_time%year = iyearz
     end_time%month = imonthz
     end_time%date = idatez
     end_time%time = int(itimez * 0.01) * 3600.0   &
          + (itimez * 0.01 - int(itimez*0.01))*100.0*60.0

     if (runtype == 'INITIAL') then
        
        ! The namelist variables in this section either must not be changed on a
        ! history restart or changing them would be irrelevant.  Thus, they are
        ! only copied to main model memory if this is not a history restart.
        
        ! set current time to initial time here.  If this is a history run,
        ! reset current time in subroutine history_start.
        current_time%year = iyeara
        current_time%month = imontha
        current_time%date = idatea
        current_time%time = int(itimea * 0.01) * 3600.0   &
             + (itimea * 0.01 - int(itimea*0.01))*100.0*60.0
        
        time = 0.0
     elseif (runtype == 'HISTORY') then
        
        ! The namelist variables in this section either must not be changed on a
        ! history restart or changing them would be irrelevant.  Thus, they are
        ! only copied to main model memory if this is not a history restart.
        ! If this is a history run,
        ! reset current time in subroutine history_start.
        
        current_time%year  = iyearh
        current_time%month = imonthh
        current_time%date  = idateh
        current_time%time  = int(itimeh * 0.01) * 3600.0   &
             + (itimeh * 0.01 - int(itimeh*0.01))*100.0*60.0
        
        ! Calculate the current time
        call date_2_seconds (iyearh,imonthh,idateh,itimeh*100, &
             iyeara,imontha,idatea,itimea*100,time)
        
     end if
     
     ! Sorting up the chosen PFTs
     where (include_these_pft < 1) include_these_pft=huge(1)
     call sort_up(include_these_pft,n_pft)
     
     !  Determine the length of simuation
     call date_2_seconds (iyearz,imonthz,idatez,itimez*100, &
          iyeara,imontha,idatea,itimea*100,timmax)

     return
   end subroutine read_ednl
   
!------------------------------------------------------------------------------------------!
  
subroutine copy_in_bramsnl(expnme_b, runtype_b, itimez_b, idatez_b, &
     imonthz_b, iyearz_b, itimea_b, idatea_b, imontha_b, iyeara_b,  &
     radfrq_b,nnxp_b,nnyp_b,deltax_b,deltay_b,polelat_b,polelon_b,  &
     centlat_b,centlon_b,nstratx_b,nstraty_b,iclobber_b,nzg_b,nzs_b,&
     isoilflg_b,nslcon_b,slz_b,slmstr_b,stgoff_b,zrough_b,ngrids_b)
  
  use misc_coms, only: expnme, runtype, itimez, idatez, imonthz, iyearz, &
       itimea, idatea, imontha, iyeara, iclobber,radfrq
  
  use grid_coms, only: centlon,centlat,deltax,deltay,nnxp,nnyp,nstratx, &
       nstraty,polelat,polelon,ngrids,timmax,time,nzg, nzs

  use soil_coms, only: isoilflg, nslcon, slmstr, zrough, slz,stgoff

  use grid_dims,  only: maxgrds,nzgmax

  implicit none
  
  character(len=64)         , intent(in) :: expnme_b   ! experiment name
  integer                   , intent(in) :: ngrids_b   ! how many grids
  character(len=16)         , intent(in) :: runtype_b
  integer                   , intent(in) :: iyeara_b
  integer                   , intent(in) :: imontha_b
  integer                   , intent(in) :: idatea_b
  integer                   , intent(in) :: itimea_b
  integer                   , intent(in) :: iyearz_b
  integer                   , intent(in) :: imonthz_b
  integer                   , intent(in) :: idatez_b
  integer                   , intent(in) :: itimez_b
  real                      , intent(in) :: radfrq_b
  integer                   , intent(in) :: iclobber_b
  integer,dimension(maxgrds), intent(in) :: nnxp_b
  integer,dimension(maxgrds), intent(in) :: nnyp_b
  real                      , intent(in) :: deltax_b
  real                      , intent(in) :: deltay_b
  real                      , intent(in) :: polelat_b
  real                      , intent(in) :: polelon_b
  real   ,dimension(maxgrds), intent(in) :: centlat_b ! grid center latitude (degrees)
  real   ,dimension(maxgrds), intent(in) :: centlon_b ! grid center longitude (degrees)
  integer,dimension(maxgrds), intent(in) :: nstratx_b ! nest ratio for next coarser grid
  integer,dimension(maxgrds), intent(in) :: nstraty_b ! nest ratio for next coarser grid
  integer                   , intent(in) :: nzg_b     ! soil layers
  integer                   , intent(in) :: nzs_b     ! snow layers
  integer,dimension(maxgrds), intent(in) :: isoilflg_b
  integer                   , intent(in) :: nslcon_b
  real                      , intent(in) :: zrough_b
  real, dimension(nzgmax)   , intent(in) :: slmstr_b
  real, dimension(nzgmax)   , intent(in) :: stgoff_b
  real, dimension(nzgmax)   , intent(in) :: slz_b

  expnme   = expnme_b
  ngrids   = ngrids_b
  runtype  = runtype_b
  itimez   = itimez_b
  idatez   = idatez_b
  imonthz  = imonthz_b
  iyearz   = iyearz_b
  itimea   = itimea_b
  idatea   = idatea_b
  imontha  = imontha_b
  iyeara   = iyeara_b

  radfrq   = radfrq_b

  iclobber = iclobber_b

  centlon  = centlon_b
  centlat  = centlat_b
  deltax   = deltax_b
  deltay   = deltay_b
  nnxp     = nnxp_b
  nnyp     = nnyp_b
  nstratx  = nstratx_b
  nstraty  = nstraty_b
  polelat  = polelat_b
  polelon  = polelon_b
  ngrids   = ngrids_b

  nzg      = nzg_b
  nzs      = nzs_b

  slz      = slz_b
  slmstr   = slmstr_b
  stgoff   = stgoff_b
  zrough   = zrough_b

  isoilflg = isoilflg_b


  
  return
end subroutine copy_in_bramsnl
