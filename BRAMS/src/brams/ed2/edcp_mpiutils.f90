subroutine masterput_ednl(mainnum)
  
  use max_dims,        only: str_len,n_pft,maxpvars,nzgmax,str_len_short,maxgrds,max_soi   &
                           ,max_ed_regions,n_pft
  use misc_coms,       only: expnme, runtype,itimea,iyeara,imontha,idatea ,itimez,iyearz   &
                           ,imonthz,idatez,dtlsm,radfrq,ifoutput,idoutput,imoutput         &
                           ,iclobber,frqfast,sfilin,ffilout,ied_init_mode,ed_inputs_dir    &
                           ,integration_scheme,end_time,current_time,sfilout,frqstate      &
                           ,isoutput,iprintpolys,printvars,npvars,pfmtstr,ipmin,ipmax      &
                           ,iedcnfgf,iyoutput,outfast,outstate,unitfast,unitstate
  use ed_misc_coms,only: attach_metadata
  use grid_coms,       only: nzg,nzs,ngrids,nnxp,nnyp,time,timmax
  use soil_coms,       only: isoilflg,nslcon,slz,slmstr,stgoff,veg_database,soil_database  &
                            ,soilstate_db,soildepth_db,isoilstateinit,isoildepthflg        &
                            ,isoilbc,runoff_time,zrough
  use met_driver_coms, only: initial_co2,lapse_scheme
  use mem_sites,       only: edres,maxpatch,maxcohort
  use physiology_coms, only: istoma_scheme,n_plant_lim
  use phenology_coms , only: iphen_scheme,repro_scheme,iphenys1,iphenysf,iphenyf1,iphenyff &
                             ,phenpath
  use decomp_coms,     only: n_decomp_lim
  use pft_coms,        only: include_these_pft,pft_1st_check
  use disturb_coms,    only: include_fire,ianth_disturb, treefall_disturbance_rate
  use optimiz_coms,    only: ioptinpt



  implicit none
  include 'mpif.h'
  integer :: ierr,mainnum
  integer :: npv
  
  call MPI_Bcast(time,1,MPI_DOUBLE_PRECISION,mainnum,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(timmax,1,MPI_DOUBLE_PRECISION,mainnum,MPI_COMM_WORLD,ierr)

  call MPI_Bcast(current_time%year,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(current_time%month,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(current_time%date,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(current_time%time,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(current_time%ifirst,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)

  call MPI_Bcast(end_time%year,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(end_time%month,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(end_time%date,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(end_time%time,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(end_time%ifirst,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)

  call MPI_Bcast(dtlsm,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(ifoutput,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(idoutput,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(imoutput,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(iyoutput,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(isoutput,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(attach_metadata,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(outfast,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(outstate,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(unitfast,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(unitstate,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(ffilout,str_len,MPI_CHARACTER,mainnum,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(sfilout,str_len,MPI_CHARACTER,mainnum,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(ied_init_mode,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(edres,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(sfilin,str_len,MPI_CHARACTER,mainnum,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(veg_database,str_len,MPI_CHARACTER,mainnum,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(soil_database,str_len,MPI_CHARACTER,mainnum,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(ed_inputs_dir,str_len,MPI_CHARACTER,mainnum,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(soilstate_db,str_len,MPI_CHARACTER,mainnum,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(soildepth_db,str_len,MPI_CHARACTER,mainnum,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(isoilstateinit,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(isoildepthflg,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(isoilbc,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(integration_scheme,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(istoma_scheme,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(iphen_scheme,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(repro_scheme,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(lapse_scheme,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(n_plant_lim,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(n_decomp_lim,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(include_fire,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(ianth_disturb,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(include_these_pft,n_pft,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(pft_1st_check,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(maxpatch,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(maxcohort,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(treefall_disturbance_rate,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(runoff_time,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(iprintpolys,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(npvars,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
  do npv=1,npvars
     call MPI_Bcast(printvars(npv),str_len_short,MPI_CHARACTER,mainnum,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(pfmtstr(npv),str_len_short,MPI_CHARACTER,mainnum,MPI_COMM_WORLD,ierr)
  end do
  call MPI_Bcast(ipmin,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(ipmax,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(initial_co2,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(iphenys1,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(iphenysf,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(iphenyf1,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(iphenyff,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(iedcnfgf,str_len,MPI_CHARACTER,mainnum,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(phenpath,str_len,MPI_CHARACTER,mainnum,MPI_COMM_WORLD,ierr)
  
  ! These are ED2 variables that have the same function as certain BRAMS
  ! namelist variables under a different name
  ! --------------------------------------------------------------------
  
  call MPI_Bcast(frqfast,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(frqstate,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
  
  call MPI_Bcast(expnme,str_len,MPI_CHARACTER,mainnum,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(ngrids,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(runtype,str_len,MPI_CHARACTER,mainnum,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(itimez,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(idatez,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(imonthz,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(iyearz,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(itimea,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(idatea,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(imontha,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(iyeara,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
  
  call MPI_Bcast(radfrq,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
  
  call MPI_Bcast(iclobber,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
  
  ! We should only need the grid sizes, because we expect to import
  ! the coordinates and sizes to allocate our polygons,
  ! The geometry is handled by BRAMS
  
  call MPI_Bcast(nnxp,maxgrds,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(nnyp,maxgrds,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
  
  call MPI_Bcast(nzg,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(nzs,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
  
  call MPI_Bcast(slz,nzgmax,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(slmstr,nzgmax,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(stgoff,nzgmax,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(zrough,nzgmax,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
  
  
  return
end subroutine masterput_ednl

!===================================================================
subroutine nodeget_ednl(master_num)
  
  use max_dims,        only: str_len,max_soi,max_ed_regions,nzgmax,n_pft,maxgrds,maxpvars &
                            ,str_len_short
  use misc_coms,       only: expnme, runtype,itimea,iyeara,imontha,idatea ,itimez,iyearz  &
                           ,imonthz,idatez,dtlsm,radfrq,ifoutput,idoutput,imoutput        &
                           ,iclobber,frqfast,sfilin,ffilout,ied_init_mode,ed_inputs_dir   &
                           ,integration_scheme,end_time,current_time,sfilout,frqstate     &
                           ,isoutput,iprintpolys,printvars,npvars,pfmtstr,ipmin,ipmax     &
                           ,iedcnfgf,iyoutput,outfast,outstate,unitfast,unitstate
  use ed_misc_coms,only: attach_metadata
  use grid_coms,       only: nzg,nzs,ngrids,nnxp,nnyp,time,timmax
  use soil_coms,       only: isoilflg,nslcon,slz,slmstr,stgoff,veg_database,soil_database &
                             ,soilstate_db,soildepth_db,isoilstateinit,isoildepthflg       &
                             ,isoilbc,runoff_time,zrough
  use met_driver_coms, only: initial_co2,lapse_scheme
  use mem_sites,       only: edres,maxpatch,maxcohort
  use physiology_coms, only: istoma_scheme,n_plant_lim
  use phenology_coms , only: iphen_scheme,repro_scheme,iphenys1,iphenysf,iphenyf1,iphenyff &
                            ,phenpath
  use decomp_coms,     only: n_decomp_lim
  use pft_coms,        only: include_these_pft,pft_1st_check
  use disturb_coms,    only: include_fire,ianth_disturb, treefall_disturbance_rate
  use optimiz_coms,    only: ioptinpt


  implicit none
  include 'mpif.h'
  integer :: ierr,master_num
  integer :: npv

  call MPI_Bcast(time,1,MPI_DOUBLE_PRECISION,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(timmax,1,MPI_DOUBLE_PRECISION,master_num,MPI_COMM_WORLD,ierr)  

  call MPI_Bcast(current_time%year,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(current_time%month,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(current_time%date,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(current_time%time,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(current_time%ifirst,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)

  call MPI_Bcast(end_time%year,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(end_time%month,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(end_time%date,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(end_time%time,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(end_time%ifirst,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)

  call MPI_Bcast(dtlsm,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(ifoutput,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(idoutput,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(imoutput,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(iyoutput,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(isoutput,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(attach_metadata,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(outfast,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(outstate,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(unitfast,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(unitstate,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(ffilout,str_len,MPI_CHARACTER,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(sfilout,str_len,MPI_CHARACTER,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(ied_init_mode,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(edres,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(sfilin,str_len,MPI_CHARACTER,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(veg_database,str_len,MPI_CHARACTER,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(soil_database,str_len,MPI_CHARACTER,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(ed_inputs_dir,str_len,MPI_CHARACTER,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(soilstate_db,str_len,MPI_CHARACTER,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(soildepth_db,str_len,MPI_CHARACTER,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(isoilstateinit,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(isoildepthflg,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(isoilbc,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(integration_scheme,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(istoma_scheme,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(iphen_scheme,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(repro_scheme,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(lapse_scheme,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(n_plant_lim,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(n_decomp_lim,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(include_fire,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(ianth_disturb,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(include_these_pft,n_pft,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(pft_1st_check,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(maxpatch,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(maxcohort,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(treefall_disturbance_rate,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(runoff_time,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(iprintpolys,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(npvars,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
  do npv=1,npvars
     call MPI_Bcast(printvars(npv),str_len_short,MPI_CHARACTER,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(pfmtstr(npv),str_len_short,MPI_CHARACTER,master_num,MPI_COMM_WORLD,ierr)
  end do
  call MPI_Bcast(ipmin,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(ipmax,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(initial_co2,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(iphenys1,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(iphenysf,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(iphenyf1,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(iphenyff,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(iedcnfgf,str_len,MPI_CHARACTER,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(phenpath,str_len,MPI_CHARACTER,master_num,MPI_COMM_WORLD,ierr)
  
  ! These are ED2 variables that have the same function as certain BRAMS
  ! namelist variables under a different name
  ! --------------------------------------------------------------------
  
  call MPI_Bcast(frqfast,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(frqstate,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
  
  call MPI_Bcast(expnme,str_len,MPI_CHARACTER,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(ngrids,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(runtype,str_len,MPI_CHARACTER,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(itimez,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(idatez,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(imonthz,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(iyearz,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(itimea,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(idatea,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(imontha,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(iyeara,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
  
  call MPI_Bcast(radfrq,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
  
  call MPI_Bcast(iclobber,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
  
  ! We should only need the grid sizes, because we expect to import
  ! the coordinates and sizes to allocate our polygons,
  ! The geometry is handled by BRAMS
  
  call MPI_Bcast(nnxp,maxgrds,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(nnyp,maxgrds,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
  
  call MPI_Bcast(nzg,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(nzs,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
  
  call MPI_Bcast(slz,nzgmax,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(slmstr,nzgmax,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(stgoff,nzgmax,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(zrough,nzgmax,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
  
  return
end subroutine nodeget_ednl
