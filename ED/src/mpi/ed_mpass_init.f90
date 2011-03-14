!==========================================================================================!
!==========================================================================================!
!    This routine gives basic processor ID info to the nodes.                              !
!------------------------------------------------------------------------------------------!
subroutine ed_masterput_processid(nproc,headnode_num,masterworks,par_run)

   use ed_para_coms, only: mainnum,nmachs,machsize,machnum
   use ed_node_coms, only: mynum,nnodetot,sendnum,recvnum,master_num,machs

   implicit none
   integer :: headnode_num,nproc
   include 'mpif.h'
   integer :: nm
   integer :: ierr
   integer :: par_run
   logical :: masterworks

   mainnum=headnode_num
   master_num=headnode_num
   nmachs=nproc

   if (masterworks) then
     mynum=nmachs+1
     nnodetot=machsize
     sendnum=1
     recvnum=nmachs
   else
     mynum=-1
     nnodetot=nmachs
     sendnum=-1
     recvnum=-1
   end if


   if (par_run == 0) return

   do nm=1,nmachs
      machnum(nm)=nm
      machs(nm)=nm
   enddo

   machs(machsize)=0  !Thats me!!

   do nm=1,nmachs
     call MPI_Send(mainnum,1,MPI_INTEGER,machnum(nm),311,MPI_COMM_WORLD,ierr)
     call MPI_Send(machnum(nm),1,MPI_INTEGER,machnum(nm),312,MPI_COMM_WORLD,ierr)
     call MPI_Send(nm,1,MPI_INTEGER,machnum(nm),313,MPI_COMM_WORLD,ierr)
     call MPI_Send(nmachs,1,MPI_INTEGER,machnum(nm),314,MPI_COMM_WORLD,ierr)
     call MPI_Send(machnum,nmachs,MPI_INTEGER,machnum(nm),315,MPI_COMM_WORLD,ierr)
     call MPI_Send(machsize,1,MPI_INTEGER,machnum(nm),316,MPI_COMM_WORLD,ierr)
   enddo



   return
end subroutine ed_masterput_processid
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!   This subroutine is responsible for sending all the namelist-related information to all !
! the nodes. This is done by the head node, because the head node has read this            !
! information to define the polygon distribution in regional runs (and in the future, the  !
! patches for the POI runs).                                                               !
!------------------------------------------------------------------------------------------!
subroutine ed_masterput_nl(par_run)
   use ed_para_coms,    only: mainnum
   use ed_max_dims,     only: str_len,max_poi,max_ed_regions,nzgmax,n_pft,maxgrds,maxpvars
   use ed_misc_coms,    only: expnme, runtype,itimea,iyeara,imontha,idatea ,itimez,iyearz  &
                             ,imonthz,idatez,dtlsm,radfrq,ifoutput,idoutput,imoutput       &
                             ,iqoutput,itoutput,iyoutput,iclobber,frqfast,sfilin,ffilout   &
                             ,ied_init_mode,thsums_database,integration_scheme,end_time    &
                             ,current_time,sfilout,frqstate,isoutput,iprintpolys,printvars &
                             ,pfmtstr,ipmin,ipmax,iedcnfgf,outfast,outstate,out_time_fast  &
                             ,out_time_state,nrec_fast,nrec_state,irec_fast,irec_state     &
                             ,unitfast,unitstate,event_file,itimeh,iyearh,imonthh,idateh   &
                             ,ndcycle

   use ed_misc_coms   , only: attach_metadata
   use canopy_air_coms, only: icanturb, i_blyr_condct, isfclyrm, ustmin, gamm, gamh        &
                             ,tprandtl,vkopr, ggfact, vh2vr, vh2dh
   use grid_coms,       only: nzg,nzs,ngrids,nnxp,nnyp,deltax,deltay,polelat,polelon       &
                             ,centlat,centlon,time,timmax,nstratx,nstraty
   use soil_coms,       only: isoilflg,nslcon,slxclay,slxsand,slz,slmstr,stgoff,veg_database,soil_database &
                             ,soilstate_db,soildepth_db,isoilstateinit,isoildepthflg       &
                             ,isoilbc,runoff_time,betapower,zrough,layer_index,nlon_lyr    &
                             ,nlat_lyr
   use met_driver_coms, only: ed_met_driver_db,imettype,ishuffle,metcyc1,metcycf           &
                             ,initial_co2, lapse_scheme
   use mem_polygons,    only: n_poi,n_ed_region,grid_type,grid_res,poi_lat,poi_lon         &
                             ,ed_reg_latmin,ed_reg_latmax,ed_reg_lonmin,ed_reg_lonmax      &
                             ,edres,maxpatch,maxcohort
   use physiology_coms, only: istoma_scheme, h2o_plant_lim, n_plant_lim, vmfact, mfact     &
                            , kfact, gamfact, lwfact, thioff, icomppt, quantum_efficiency_T
   use phenology_coms , only: iphen_scheme,iphenys1,iphenysf,iphenyf1,iphenyff,phenpath    &
                             ,repro_scheme, radint, radslp
   use decomp_coms,     only: n_decomp_lim
   use pft_coms,        only: include_these_pft,agri_stock,plantation_stock,pft_1st_check
   use disturb_coms,    only: include_fire,ianth_disturb, treefall_disturbance_rate        &
                             ,lu_database,plantation_file,lu_rescale_file
   use optimiz_coms,    only: ioptinpt
   use canopy_radiation_coms, only : crown_mod
   use rk4_coms,        only: rk4_tolerance, ibranch_thermo, ipercol

   implicit none
   include 'mpif.h'
   integer :: ierr
   integer :: par_run
   integer :: n
   if (par_run == 0 ) return

   !----- First, the namelist-derived type, before I forget... ----------------------------!
   call MPI_Bcast(ngrids,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
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

   !----- Now the namelist ----------------------------------------------------------------!
   call MPI_Bcast(expnme,str_len,MPI_CHARACTER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(runtype,str_len,MPI_CHARACTER,mainnum,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(imontha,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(idatea,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(iyeara,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(itimea,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(imonthh,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(idateh,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(iyearh,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(itimeh,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(imonthz,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(idatez,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(iyearz,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(itimez,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(dtlsm,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(radfrq,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(ifoutput,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(idoutput,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(imoutput,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(iqoutput,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(iyoutput,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(itoutput,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(isoutput,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ndcycle ,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(iclobber,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(frqfast,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(outfast,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(outstate,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(unitfast,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(unitstate,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
                      
   call MPI_Bcast(ffilout,str_len,MPI_CHARACTER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ied_init_mode,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(sfilout,str_len,MPI_CHARACTER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(frqstate,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(nzg ,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(nzs ,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(isoilflg,maxgrds,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(nslcon,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(slxclay,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(slxsand,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(slz ,nzgmax,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(stgoff,nzgmax,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(slmstr,nzgmax,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   
   do n=1, maxgrds
      call MPI_Bcast(sfilin         (n),str_len,MPI_CHARACTER,mainnum,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(veg_database   (n),str_len,MPI_CHARACTER,mainnum,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(soil_database  (n),str_len,MPI_CHARACTER,mainnum,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(lu_database    (n),str_len,MPI_CHARACTER,mainnum,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(plantation_file(n),str_len,MPI_CHARACTER,mainnum,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(lu_rescale_file(n),str_len,MPI_CHARACTER,mainnum,MPI_COMM_WORLD,ierr)
   end do

   call MPI_Bcast(thsums_database ,str_len,MPI_CHARACTER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(soilstate_db    ,str_len,MPI_CHARACTER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(soildepth_db    ,str_len,MPI_CHARACTER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ed_met_driver_db,str_len,MPI_CHARACTER,mainnum,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(isoilstateinit,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(isoildepthflg,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(isoilbc,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(n_poi,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(n_ed_region,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(grid_type,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   
   call MPI_Bcast(grid_res,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(poi_lat,max_poi,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(poi_lon,max_poi,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ed_reg_latmin,max_ed_regions,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ed_reg_latmax,max_ed_regions,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ed_reg_lonmin,max_ed_regions,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ed_reg_lonmax,max_ed_regions,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(nnxp,maxgrds,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(nnyp,maxgrds,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(nstratx,maxgrds,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(nstraty,maxgrds,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(deltax,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(deltay,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(polelat,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(polelon,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(centlat,maxgrds,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(centlon,maxgrds,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(integration_scheme,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(rk4_tolerance,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ibranch_thermo,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(istoma_scheme,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(iphen_scheme,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(repro_scheme,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(radint,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)   
   call MPI_Bcast(radslp,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
      
   call MPI_Bcast(lapse_scheme,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(crown_mod,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(h2o_plant_lim,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(vmfact,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(mfact,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(kfact,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(gamfact,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(lwfact,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(thioff,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(icomppt,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(quantum_efficiency_T,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(n_plant_lim,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(n_decomp_lim,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(include_fire,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ianth_disturb,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(include_these_pft,n_pft,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(agri_stock,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(plantation_stock,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(pft_1st_check,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(icanturb,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(i_blyr_condct,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(isfclyrm,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ipercol,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(treefall_disturbance_rate,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(runoff_time,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(betapower,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ustmin,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(gamm,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(gamh,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(tprandtl,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(vkopr,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(vh2vr,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(vh2dh,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ggfact,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(iprintpolys,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   do n=1,maxpvars
      call MPI_Bcast(printvars(n),str_len,MPI_CHARACTER,mainnum,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(pfmtstr(n),str_len,MPI_CHARACTER,mainnum,MPI_COMM_WORLD,ierr)
   end do
   call MPI_Bcast(ipmin,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ipmax,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   
   call MPI_Bcast(imettype,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ishuffle,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(metcyc1,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(metcycf,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(initial_co2,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(iphenys1,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(iphenysf,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(iphenyf1,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(iphenyff,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(iedcnfgf,str_len,MPI_CHARACTER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(phenpath,str_len,MPI_CHARACTER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(event_file,str_len,MPI_CHARACTER,mainnum,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(maxpatch,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(maxcohort,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)

  
   call MPI_Bcast(ioptinpt,str_len,MPI_CHARACTER,mainnum,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(zrough,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(edres,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(attach_metadata,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)

!------------------------------------------------------------------------------------------!
!   One last thing to send is the layer index based on the soil_depth. It is not really a  !
! namelist thing, but it is still a setup variable.                                        !
!------------------------------------------------------------------------------------------!
   call MPI_Barrier(MPI_COMM_WORLD,ierr) ! Just to wait until the matrix is allocated
   call MPI_Bcast(layer_index,nlat_lyr*nlon_lyr,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)

   return
end subroutine ed_masterput_nl
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine ed_masterput_met_header(par_run)
!------------------------------------------------------------------------------------------!
!    This subroutine sends the met driver information to the nodes, which can then read    !
! the hdf5 in parallel                                                                     !
!------------------------------------------------------------------------------------------!
   use ed_para_coms, only: mainnum
   use ed_max_dims, only: max_met_vars,str_len
   use met_driver_coms, only: nformats, met_names, met_nlon,   &
        met_nlat, met_dx, met_dy, met_xmin, met_ymin, met_nv,   &
        met_vars, met_frq, met_interp, ed_met_driver_db, no_ll,  &
        metname_len,metvars_len
   
   implicit none
   include 'mpif.h'
   integer, intent(in) :: par_run
   integer             :: ierr, nsize,f,v

   if (par_run == 0) return
   
   nsize=nformats*max_met_vars

!----- First I send the scalars -----------------------------------------------------------!
   call MPI_Bcast (nformats,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)

   
!------------------------------------------------------------------------------------------!
!   Here I need a MPI Barrier. The master has the variables already allocated, but I need  !
! the nodes with their structures already allocated before I proceed sending the info      !
!------------------------------------------------------------------------------------------!
   call MPI_Barrier(MPI_COMM_WORLD,ierr)

   do f=1,nformats
     call MPI_Bcast(met_names(f),metname_len,MPI_CHARACTER,mainnum,MPI_COMM_WORLD,ierr)
   end do
   
   call MPI_Bcast(met_nlon,nformats,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(met_nlat,nformats,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(met_dx,nformats,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(met_dy,nformats,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(met_xmin,nformats,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(met_ymin,nformats,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(met_nv,nformats,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(no_ll,nformats,MPI_LOGICAL,mainnum,MPI_COMM_WORLD,ierr)
  
   do f=1,nformats
      do v=1,max_met_vars
         call MPI_Bcast(met_vars(f,v),metvars_len,MPI_CHARACTER,mainnum,MPI_COMM_WORLD,ierr)
      end do
   end do
   
   call MPI_Bcast(met_frq,nsize,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(met_interp,nsize,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)

   return
end subroutine ed_masterput_met_header
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine ed_masterput_poly_dims(par_run,masterworks)
   use ed_state_vars , only : gdpy         & ! intent(in)
                            , py_off       ! ! intent(in)
   use grid_coms     , only : nnxp         & ! intent(in)
                            , nnyp         & ! intent(in)
                            , ngrids       ! ! intent(in)
   use ed_work_vars  , only : work_v       & ! intent(in)
                            , npolys_run   ! ! intent(in)
   use ed_para_coms  , only : mainnum      & ! intent(in)
                            , nmachs       & ! intent(in)
                            , loadmeth     ! ! intent(in)
   use mem_polygons  , only : n_ed_region  & ! intent(in)
                            , n_poi        ! ! intent(in)
   implicit none
   include 'mpif.h'
   !----- Local constants. ----------------------------------------------------------------!
   integer                     , parameter   :: nmethods = 3
   !----- Arguments. ----------------------------------------------------------------------!
   integer                     , intent(in)  :: par_run
   logical                     , intent(in)  :: masterworks
   !----- Local variables. ----------------------------------------------------------------!
   integer, dimension(:,:)     , allocatable :: machind
   integer, dimension(:)       , allocatable :: mpolys
   integer, dimension(:)       , allocatable :: moffset
   integer                                   :: ierr
   integer                                   :: npolys
   integer                                   :: ifm
   integer                                   :: ipy
   integer                                   :: imach
   integer                                   :: imeth
   integer                                   :: ibest
   integer                                   :: maxnmachs
   integer                                   :: ntotmachs
   real   , dimension(:,:)     , allocatable :: machload
   real   , dimension(:)       , allocatable :: thiswork
   real   , dimension(:)       , allocatable :: cumwork
   real   , dimension(nmethods)              :: maxload
   real                                      :: totalwork
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !    This is a logical flag to test wheter the master should also simulate polygons.    !
   ! This is true for offline runs, but false for coupled runs.                            !
   !---------------------------------------------------------------------------------------!
   if (masterworks) then
      ntotmachs=nmachs+1
   else
      ntotmachs=nmachs
   end if


   !---------------------------------------------------------------------------------------!
   !     We currently don't allow POI simulations to be run in parallel.                   !
   !---------------------------------------------------------------------------------------!
   if (par_run == 1 .and. n_poi > 0 ) then
      write (unit=*,fmt='(a)') '--------------------------------------------------------------------'
      write (unit=*,fmt='(a)') '                                                                    '
      write (unit=*,fmt='(a)') '  Dear ED user,                                                     '
      write (unit=*,fmt='(a)') '                                                                    '
      write (unit=*,fmt='(a)') '  Thank you for choosing ED, the Ecosystem Demography Model!        '
      write (unit=*,fmt='(a)') 'We know that choosing a model is more than looking for a            '
      write (unit=*,fmt='(a)') 'numeric solver, it is an investment in your academic career.        '
      write (unit=*,fmt='(a)') 'Therefore, it is our commitment to provide our valuable researchers ' 
      write (unit=*,fmt='(a)') 'with the very best in ecosystem dynamics modelling and we want to   '
      write (unit=*,fmt='(a)') 'reaffirm that your satisfaction is our number-one priority.         '
      write (unit=*,fmt='(a)') '  Unfortunately, the option of using parallelism under POI is not   '
      write (unit=*,fmt='(a)') 'available yet, and we would like to offer our most sincere apologies'
      write (unit=*,fmt='(a)') 'for any inconvenience that this may have caused. Our team is        '
      write (unit=*,fmt='(a)') 'currently working on implementing this capability, and we hope to   '
      write (unit=*,fmt='(a)') 'achieve this goal in a very near future, thus meeting your needs of '
      write (unit=*,fmt='(a)') 'a fast and reliable model. Meanwhile you may find useful to run the '
      write (unit=*,fmt='(a)') 'SOI runs in serial mode instead.                                    '
      write (unit=*,fmt='(a)') '                                                                    '
      write (unit=*,fmt='(a)') '--------------------------------------------------------------------' 
      call fatal_error('Parallel version of POI runs not available.'                       &
                      , 'ed_masterput_poly_dims','ed_mpass_init.f90')
   end if
  

   !---------------------------------------------------------------------------------------!
   !     Now we decide how many polygons we send for each node.  If this is a serial       !
   ! (single node) run, then all nodes go to the only node available.  Otherwise, we       !
   ! try to split the nodes evenly.  This doesn't necessarily mean giving each node the    !
   ! same number of polygons, because different environments may require very different    !
   ! number of operations.  So if we have information on how the model has been behaving   !
   ! (from a history start, for example), we use that information to split the work load   !
   ! more evenly.  Otherwise, if this is a near bare ground or an ED-1/ED-2.0 restart run, !
   ! we assume that all polygons are equally expensive.                                    !
   !---------------------------------------------------------------------------------------!
   if (par_run == 0) then
     do ifm=1,n_ed_region
        gdpy(1,ifm)   = npolys_run(ifm)
        py_off(1,ifm) = 0
     end do
     do ifm=n_ed_region+1,ngrids
        gdpy(1,ifm)   = 1
        py_off(1,ifm) = 0
     end do
   else
      !----- Parallel run, split the grids into the nodes. --------------------------------!
      gridloop: do ifm=1,ngrids
         npolys    = npolys_run(ifm)

         !---------------------------------------------------------------------------------!
         !     We will now try three different methods to put the nodes together.  The     !
         ! final choice will be the one with the lowest maximum load, because this is the  !
         ! one with the smallest load disparity.  But first we allocate some scratch       !
         ! arrays.                                                                         !
         !---------------------------------------------------------------------------------!
         allocate(thiswork(npolys),cumwork(npolys),machind(npolys,3))
         allocate(machload(ntotmachs,3),mpolys(ntotmachs),moffset(ntotmachs))
         !----- Copy the workload and also find the cumulative sum of workload. -----------!
         call atob(npolys,work_v(ifm)%work,thiswork)
         call atob(npolys,work_v(ifm)%work,cumwork)
         call cumsum(npolys,cumwork)
         totalwork = cumwork(npolys)

         !---------------------------------------------------------------------------------!
         !     Worksum is the sum of the relative workload. This will be always less than  !
         ! or equal to the number of polygons.  In case we have no idea on how the         !
         ! polygons behave (if they are slow or fast), then all polygons received workload !
         ! 1., and the total workload will be the same as the number of polygons.  If we   !
         ! do know the distribution from the history file, and the workload is very skew-  !
         ! ed (some polygons are extremely slow compared to the others), worksum will be   !
         ! significantly less than the number of polygons.                                 !
         !---------------------------------------------------------------------------------!
         maxnmachs = ceiling(totalwork)
         
         !---------------------------------------------------------------------------------!
         !     If nmachsmax is less than the number of machines, stop the run.  The MPI    !
         ! process contains way too many nodes, and this would be a waste of resources.    !
         ! At the run crash, give this information to the user...                          !
         !---------------------------------------------------------------------------------!
         if (maxnmachs < ntotmachs .and. loadmeth /= 2) then
            write (unit=*,fmt='(a)') '----------------------------------------------------'
            write (unit=*,fmt='(a)') '       The number of requested nodes exceeds the    '
            write (unit=*,fmt='(a)') ' maximum number of nodes in which you would still   '
            write (unit=*,fmt='(a)') ' take advantage of parallel processing.  This can   '
            write (unit=*,fmt='(a)') ' be because you have requested a number of nodes    '
            write (unit=*,fmt='(a)') ' that exceeds the number of polygons, or because    '
            write (unit=*,fmt='(a)') ' in the past some of your polygons were so much     '
            write (unit=*,fmt='(a)') ' slower than the others that having more nodes      '
            write (unit=*,fmt='(a)') ' wouldn''t help.  Reduce the number of nodes and    '
            write (unit=*,fmt='(a)') ' try again...                                       '
            write (unit=*,fmt='(a)') '----------------------------------------------------'
            write (unit=*,fmt='(a,1x,i6)') ' # of requested slave nodes :',ntotmachs
            write (unit=*,fmt='(a,1x,i6)') ' # of polygons              :',npolys 
            write (unit=*,fmt='(a,1x,i6)') ' Max. # of nodes needed     :',maxnmachs 
            write (unit=*,fmt='(a)') '----------------------------------------------------'
            call fatal_error('Requested number of nodes exceeds the maximum needed.'       &
                            ,'ed_masterput_poly_dims','ed_mpass_init.f90')
         end if


         do ipy = 1, npolys
            !------------------------------------------------------------------------------!
            ! Method 1.  Load will be evenly distributed accross the nodes.  The node      !
            !            index will be defined by the next integer of the relative         !
            !            position of the cumulative sum split into ntotmachs categories.   !
            !------------------------------------------------------------------------------!
            machind(ipy,1) = ceiling(cumwork(ipy)*real(ntotmachs)/totalwork)

            !------------------------------------------------------------------------------!
            ! Method 2.  We simply split the polygons evenly, without bothering about the  !
            !            workload.                                                         !
            !------------------------------------------------------------------------------!
            machind(ipy,2) = ceiling(real(ipy)*real(ntotmachs)/real(npolys))

            !------------------------------------------------------------------------------!
            ! Method 3.  Similar to 1, but here we round the index instead of finding the  !
            !            next integer.  Because of that we use the number of machines      !
            !            minus one.  This is theoretically less efficient, but sometimes   !
            !            it can be less just because a large workload polygon is next to a !
            !            small one such that it would cause an almost empty node next to   !
            !            very busy one.                                                    !
            !------------------------------------------------------------------------------!
            machind(ipy,3) = 1 + nint(cumwork(ipy)*real(ntotmachs-1)/totalwork)
         end do

         !---------------------------------------------------------------------------------!
         !    We now play it safe and ensure all indices are within bounds.                !
         !---------------------------------------------------------------------------------!
         where (machind < 1)
            machind = 1
         end where
         where (machind > ntotmachs)
            machind = ntotmachs
         end where

         !---------------------------------------------------------------------------------!
         !     We now loop over the machine indices and find the workload that each        !
         ! node will receive, as well as the maximum load that each method got.            !
         !---------------------------------------------------------------------------------!
         machload(:,:) = 0.0
         do imeth = 1, nmethods
            do imach = 1, ntotmachs
               machload(imach,imeth) = sum(thiswork,machind(:,imeth) == imach)
            end do
            maxload(imeth) = maxval(machload(:,imeth))
         end do

         !---------------------------------------------------------------------------------!
         !     The best method is the one that gives the lowest maximum load.  The average !
         ! load is actually the same for all the three methods, so the only difference is  !
         ! the range, that's why the smallest maximum is the preferred one.                !
         !---------------------------------------------------------------------------------!
         ibest = minloc(maxload,dim=1)

         if(loadmeth>0) ibest = loadmeth
         
         !----- Count how many polygons go to each node. ----------------------------------!
         do imach=1,ntotmachs
            mpolys(imach) = count(machind(:,ibest) == imach)
         end do

         !----- Now find the offset. ------------------------------------------------------!
         moffset(1) = 0
         do imach = 2, ntotmachs
            moffset(imach) = moffset(imach-1) + mpolys(imach-1)
         end do

         !---------------------------------------------------------------------------------!
         !     We now print a banner so the user can have some idea of the workload        !
         ! distribution.                                                                   !
         !---------------------------------------------------------------------------------!
         write(unit=*,fmt='(a)') '--------------------------------------------------------'
         write(unit=*,fmt='(a)') '                WORKLOAD PER POLYGON'
         write(unit=*,fmt='(a,1x,i6)')     ' - Maxnmachs   : ',maxnmachs
         write(unit=*,fmt='(a,1x,i6)')     ' - NPolys      : ',npolys
         write(unit=*,fmt='(a,1x,es13.6)') ' - Total Work  : ',totalwork
         write(unit=*,fmt='(a)') '--------------------------------------------------------'
         write(unit=*,fmt='(7(a,1x))') '   IPY','     WORKLOAD','       CUMSUM'            &
                                             ,'METH1','METH2','METH3'
         write(unit=*,fmt='(a)') '--------------------------------------------------------'
         do ipy=1,npolys
            write(unit=*,fmt='(i6,1x,2(es13.6,1x),3(i5,1x))')                              &
                                              ipy,thiswork(ipy),cumwork(ipy)               &
                                            , machind(ipy,1),machind(ipy,2),machind(ipy,3)
         end do
         write(unit=*,fmt='(a)') '--------------------------------------------------------'
         !---------------------------------------------------------------------------------!

         !---------------------------------------------------------------------------------!
         !     We now print a banner so the user can have some idea of the workload        !
         ! distribution.                                                                   !
         !---------------------------------------------------------------------------------!
         write(unit=*,fmt='(a)') '--------------------------------------------------------'
         write(unit=*,fmt='(a)') '                Node distribution'
         write(unit=*,fmt='(a,1x,i6)') ' - Best method :',ibest
         write(unit=*,fmt='(a)') '--------------------------------------------------------'
         write(unit=*,fmt='(6(a,1x))') ' IMACH','MPOLYS','OFFSET','    MACHLOAD1'          &
                                      ,'    MACHLOAD2','    MACHLOAD3'
         write(unit=*,fmt='(a)') '--------------------------------------------------------'
         do imach=1,ntotmachs
            write(unit=*,fmt='(3(i6,1x),3(es13.6,1x))') imach,mpolys(imach),moffset(imach) &
                                                       ,machload(imach,1)                  &
                                                       ,machload(imach,2)                  &
                                                       ,machload(imach,3)
         end do
         write(unit=*,fmt='(a)') '--------------------------------------------------------'
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !     Send the number of polygons and offset to the other nodes.                  !
         !---------------------------------------------------------------------------------!
         do imach=1,nmachs
            gdpy  (imach,ifm) = mpolys (imach)
            py_off(imach,ifm) = moffset(imach)
            call MPI_Bcast(mpolys (imach),1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
            call MPI_Bcast(moffset(imach),1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
         end do

         !----- Getting the message about this node itself. -------------------------------!
         gdpy  (ntotmachs,ifm) = mpolys(imach)
         py_off(ntotmachs,ifm) = moffset(imach)

         !---------------------------------------------------------------------------------!
         !     Free memory for re-allocation in the next grid (or to quit the subroutine). !
         !---------------------------------------------------------------------------------!
         deallocate(thiswork,cumwork,machind)
         deallocate(machload,mpolys,moffset)
      end do gridloop
   end if

   return
end subroutine ed_masterput_poly_dims
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine ed_masterput_worklist_info(par_run)

  use ed_max_dims, only: maxmach
  use grid_coms, only: ngrids
  use ed_work_vars , only : work_v                & ! intent(inout)
                          , work_e                & ! intent(inout)
                          , work_vecs             & ! structure
                          , ed_alloc_work_vec     & ! subroutine
                          , ed_nullify_work_vec   & ! subroutine
                          , ed_dealloc_work_vec   ! ! subroutine
  use ed_para_coms, only: nmachs,mainnum,machnum
  use soil_coms, only: isoilflg
  use ed_node_coms, only: mxp,myp,i0,j0
  use ed_state_vars, only: gdpy,py_off

  implicit none
  include 'mpif.h'
  integer, intent(in) :: par_run
  integer :: npoly
  integer :: offset
  integer :: nm
  integer :: ifm
  integer :: mpiid
  integer :: ierr

  type(work_vecs), allocatable, dimension(:)   :: sc_work

  real,            allocatable, dimension(:) :: rscratch
  integer,         allocatable, dimension(:) :: iscratch

  

  if (par_run == 1) then
     
     do nm=1,nmachs
        do ifm=1,ngrids
           
           npoly  = gdpy(nm,ifm)
           offset = py_off(nm,ifm)
           
           allocate(rscratch(npoly),iscratch(npoly))

           mpiid=maxmach*(ifm-1)+nm

           rscratch(1:npoly) = work_v(ifm)%glon(offset+1:offset+npoly)
           call MPI_Send(rscratch,npoly,MPI_REAL,machnum(nm),1300000+mpiid,MPI_COMM_WORLD,ierr)
  
           rscratch(1:npoly) = work_v(ifm)%glat(offset+1:offset+npoly)
           call MPI_Send(rscratch,npoly,MPI_REAL,machnum(nm),1400000+mpiid,MPI_COMM_WORLD,ierr)
           
           rscratch(1:npoly) = work_v(ifm)%landfrac(offset+1:offset+npoly)
           call MPI_Send(rscratch,npoly,MPI_REAL,machnum(nm),1500000+mpiid,MPI_COMM_WORLD,ierr)
           
           iscratch(1:npoly) = work_v(ifm)%ntext(offset+1:offset+npoly)
           call MPI_Send(iscratch,npoly,MPI_INTEGER,machnum(nm),1600000+mpiid,MPI_COMM_WORLD,ierr)

           iscratch(1:npoly) = work_v(ifm)%xid(offset+1:offset+npoly)
           call MPI_Send(iscratch,npoly,MPI_INTEGER,machnum(nm),1700000+mpiid,MPI_COMM_WORLD,ierr)
           
           iscratch(1:npoly) = work_v(ifm)%yid(offset+1:offset+npoly)
           call MPI_Send(iscratch,npoly,MPI_INTEGER,machnum(nm),1800000+mpiid,MPI_COMM_WORLD,ierr)



           deallocate(rscratch,iscratch)

        enddo
     enddo
     
  endif
  
  ! Deallocate each of the global vectors, and reallocate (and fill) as local vectors

  nm=nmachs+1
  allocate(sc_work(ngrids))

  do ifm=1,ngrids
     
     npoly  = gdpy(nm,ifm)
     offset = py_off(nm,ifm)
     call ed_nullify_work_vec(sc_work(ifm))
     call ed_alloc_work_vec(sc_work(ifm),npoly)

     sc_work(ifm)%glon(1:npoly)     = work_v(ifm)%glon(offset+1:offset+npoly)
     sc_work(ifm)%glat(1:npoly)     = work_v(ifm)%glat(offset+1:offset+npoly)
     sc_work(ifm)%landfrac(1:npoly) = work_v(ifm)%landfrac(offset+1:offset+npoly)
     sc_work(ifm)%ntext(1:npoly)    = work_v(ifm)%ntext(offset+1:offset+npoly)
     sc_work(ifm)%xid(1:npoly)      = work_v(ifm)%xid(offset+1:offset+npoly)
     sc_work(ifm)%yid(1:npoly)      = work_v(ifm)%yid(offset+1:offset+npoly)

     call ed_dealloc_work_vec(work_v(ifm))
     

  enddo

  deallocate(work_e,work_v)

  ! So here it will be 2 (node-style) if it is a parallel run, and 0 if it is a serial run
  call ed_mem_alloc(2*par_run) 

  allocate (work_v(ngrids))
  do ifm=1,ngrids

     npoly  = gdpy(nm,ifm)
     call ed_nullify_work_vec(work_v(ifm))
     call ed_alloc_work_vec(work_v(ifm),npoly)

     !----- Copying the scratch work structure back to the work structure. ----------------!
     work_v(ifm)%glon(1:npoly)     = sc_work(ifm)%glon(1:npoly)
     work_v(ifm)%glat(1:npoly)     = sc_work(ifm)%glat(1:npoly)
     work_v(ifm)%landfrac(1:npoly) = sc_work(ifm)%landfrac(1:npoly)
     work_v(ifm)%ntext(1:npoly)    = sc_work(ifm)%ntext(1:npoly)
     work_v(ifm)%xid(1:npoly)      = sc_work(ifm)%xid(1:npoly)
     work_v(ifm)%yid(1:npoly)      = sc_work(ifm)%yid(1:npoly)

     !----- Deallocating the scratch structure. -------------------------------------------!
     call ed_dealloc_work_vec(sc_work(ifm))

   end do
   
   deallocate(sc_work)


  
  return
end subroutine ed_masterput_worklist_info


!==========================================================================================!
!==========================================================================================!


subroutine ed_nodeget_processid(init)

  use ed_max_dims
  use ed_node_coms

  implicit none
  integer :: init

  include 'mpif.h'
  integer, dimension(MPI_STATUS_SIZE) :: status
  integer :: ierr

  if(init == 1) then
     
     call MPI_Recv(master_num,1,MPI_INTEGER,0,311,MPI_COMM_WORLD,status,ierr)
     call MPI_Recv(mchnum,1,MPI_INTEGER,0,312,MPI_COMM_WORLD,status,ierr)
     call MPI_Recv(mynum,1,MPI_INTEGER,0,313,MPI_COMM_WORLD,status,ierr)
     call MPI_Recv(nmachs,1,MPI_INTEGER,0,314,MPI_COMM_WORLD,status,ierr)
     call MPI_Recv(machs,nmachs,MPI_INTEGER,0,315,MPI_COMM_WORLD,status,ierr)
     call MPI_Recv(nnodetot,1,MPI_INTEGER,0,316,MPI_COMM_WORLD,status,ierr)
     
     recvnum = mynum-1
     sendnum = mynum+1
     if (mynum == nmachs) sendnum=0
  endif
  write(unit=*,fmt='(a,1x,i5,1x,a)') '---> Node',mynum,'got first message!'

  return
end subroutine ed_nodeget_processid
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine ed_nodeget_nl

  use ed_node_coms, only: master_num,mynum
!------------------------------------------------------------------------------------------!
!   This subroutine is responsible for getting all the namelist-related information in     !
! every node.                                                                              !
!------------------------------------------------------------------------------------------!
   use ed_max_dims,     only: str_len,max_poi,max_ed_regions,nzgmax,n_pft,maxgrds,maxpvars
   use ed_misc_coms,    only: expnme, runtype,itimea,iyeara,imontha,idatea ,itimez,iyearz  &
                             ,imonthz,idatez,dtlsm,radfrq,ifoutput,idoutput,imoutput       &
                             ,iqoutput, itoutput,iyoutput,iclobber,frqfast,sfilin,ffilout  &
                             ,ied_init_mode,thsums_database,integration_scheme,end_time    &
                             ,current_time,isoutput,sfilout,frqstate,iprintpolys,printvars &
                             ,pfmtstr,ipmin,ipmax,iedcnfgf,outfast,outstate,out_time_fast  &
                             ,out_time_state,nrec_fast,nrec_state,irec_fast,irec_state     &
                             ,unitfast,unitstate,event_file,itimeh,iyearh,imonthh,idateh   &
                             ,ndcycle

   use grid_coms,       only: nzg,nzs,ngrids,nnxp,nnyp,deltax,deltay,polelat,polelon       &
                             ,centlat,centlon,time,timmax,nstratx,nstraty
   use soil_coms,       only: isoilflg,nslcon,slxclay,slxsand,slz,slmstr,stgoff,veg_database,soil_database &
                             ,soilstate_db,soildepth_db,isoilstateinit,isoildepthflg       &
                             ,isoilbc,runoff_time,betapower,zrough,layer_index,nlon_lyr    &
                             ,nlat_lyr
   use met_driver_coms, only: ed_met_driver_db,imettype,ishuffle,metcyc1,metcycf           &
                             ,initial_co2,lapse_scheme
   use mem_polygons,    only: n_poi,n_ed_region,grid_type,grid_res,poi_lat,poi_lon         &
                             ,ed_reg_latmin,ed_reg_latmax,ed_reg_lonmin,ed_reg_lonmax      &
                             ,edres,maxpatch,maxcohort
   use physiology_coms, only: istoma_scheme, h2o_plant_lim, n_plant_lim, vmfact, mfact     &
                            , kfact, gamfact, lwfact, thioff, icomppt, quantum_efficiency_T
   use phenology_coms , only: iphen_scheme,iphenys1,iphenysf,iphenyf1,iphenyff,phenpath    &
                             ,repro_scheme, radint, radslp
   use decomp_coms,     only: n_decomp_lim
   use disturb_coms,    only: include_fire,ianth_disturb, treefall_disturbance_rate        &
                             ,lu_database,plantation_file,lu_rescale_file
   use optimiz_coms,    only: ioptinpt
   use ed_misc_coms,    only: attach_metadata
   use canopy_air_coms, only: icanturb, i_blyr_condct, isfclyrm, ustmin, gamm, gamh        &
                             ,tprandtl, vkopr, ggfact, vh2vr, vh2dh
   use pft_coms,        only: include_these_pft,agri_stock,plantation_stock,pft_1st_check
   use canopy_radiation_coms, only: crown_mod
   use rk4_coms,        only: rk4_tolerance, ibranch_thermo, ipercol

   implicit none
   include 'mpif.h'
   integer :: ierr
   integer :: n

!----- First, the namelist-derived type, before I forget... -------------------------------!
   call MPI_Bcast(ngrids,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
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

!----- Now the namelist -------------------------------------------------------------------!
   call MPI_Bcast(expnme,str_len,MPI_CHARACTER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(runtype,str_len,MPI_CHARACTER,master_num,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(imontha,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(idatea,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(iyeara,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(itimea,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(imonthh,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(idateh,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(iyearh,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(itimeh,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(imonthz,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(idatez,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(iyearz,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(itimez,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(dtlsm,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(radfrq,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(ifoutput,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(idoutput,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(imoutput,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(iqoutput,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(iyoutput,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(itoutput,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(isoutput,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ndcycle ,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(iclobber,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(frqfast,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(outfast,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(outstate,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(unitfast,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(unitstate,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
                      
   call MPI_Bcast(ffilout,str_len,MPI_CHARACTER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ied_init_mode,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(sfilout,str_len,MPI_CHARACTER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(frqstate,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(nzg ,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(nzs ,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(isoilflg,maxgrds,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(nslcon,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(slxclay,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(slxsand,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   
   call MPI_Bcast(slz ,nzgmax,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(stgoff,nzgmax,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(slmstr,nzgmax,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   
 
   
   do n=1, maxgrds
      call MPI_Bcast(sfilin         (n),str_len,MPI_CHARACTER,master_num                   &
                    ,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(veg_database   (n),str_len,MPI_CHARACTER,master_num                   &
                    ,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(soil_database  (n),str_len,MPI_CHARACTER,master_num                   &
                    ,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(lu_database    (n),str_len,MPI_CHARACTER,master_num                   &
                    ,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(plantation_file(n),str_len,MPI_CHARACTER,master_num                   &
                    ,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(lu_rescale_file(n),str_len,MPI_CHARACTER,master_num                   &
                    ,MPI_COMM_WORLD,ierr)
   end do

   call MPI_Bcast(thsums_database ,str_len,MPI_CHARACTER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(soilstate_db    ,str_len,MPI_CHARACTER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(soildepth_db    ,str_len,MPI_CHARACTER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ed_met_driver_db,str_len,MPI_CHARACTER,master_num,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(isoilstateinit,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(isoildepthflg,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(isoilbc,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(n_poi,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(n_ed_region,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(grid_type,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   
   call MPI_Bcast(grid_res,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(poi_lat,max_poi,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(poi_lon,max_poi,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ed_reg_latmin,max_ed_regions,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ed_reg_latmax,max_ed_regions,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ed_reg_lonmin,max_ed_regions,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ed_reg_lonmax,max_ed_regions,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(nnxp,maxgrds,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(nnyp,maxgrds,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(nstratx,maxgrds,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(nstraty,maxgrds,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(deltax,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(deltay,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(polelat,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(polelon,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(centlat,maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(centlon,maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(integration_scheme,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(rk4_tolerance,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ibranch_thermo,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(istoma_scheme,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(iphen_scheme,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(repro_scheme,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(radint,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(radslp,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)   
   
   call MPI_Bcast(lapse_scheme,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(crown_mod,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(h2o_plant_lim,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(vmfact,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(mfact,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(kfact,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(gamfact,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(lwfact,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(thioff,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(icomppt,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(quantum_efficiency_T,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(n_plant_lim,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(n_decomp_lim,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(include_fire,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ianth_disturb,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(include_these_pft,n_pft,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(agri_stock,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(plantation_stock,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(pft_1st_check,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(icanturb,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(i_blyr_condct,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(isfclyrm,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ipercol,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(treefall_disturbance_rate,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(runoff_time,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(betapower,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ustmin,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(gamm,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(gamh,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(tprandtl,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(vkopr,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(vh2vr,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(vh2dh,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ggfact,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(iprintpolys,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   do n=1,maxpvars
      call MPI_Bcast(printvars(n),str_len,MPI_CHARACTER,master_num,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(pfmtstr(n),str_len,MPI_CHARACTER,master_num,MPI_COMM_WORLD,ierr)
   end do
   call MPI_Bcast(ipmin,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ipmax,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   
   call MPI_Bcast(imettype,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ishuffle,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(metcyc1,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(metcycf,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(initial_co2,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(iphenys1,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(iphenysf,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(iphenyf1,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(iphenyff,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(iedcnfgf,str_len,MPI_CHARACTER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(phenpath,str_len,MPI_CHARACTER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(event_file,str_len,MPI_CHARACTER,master_num,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(maxpatch,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(maxcohort,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)

  
   call MPI_Bcast(ioptinpt,str_len,MPI_CHARACTER,master_num,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(zrough,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(edres,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   
   call MPI_Bcast(attach_metadata,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   
!------------------------------------------------------------------------------------------!
!     Receiving the layer index based on soil_depth. This is allocatable, so I first       !
! allocate, then let the master know that it is safe to send to me and I reveive the data. !
!------------------------------------------------------------------------------------------!
   if (allocated(layer_index)) deallocate(layer_index)
   allocate(layer_index(nlat_lyr,nlon_lyr))
   call MPI_Barrier(MPI_COMM_WORLD,ierr) ! Safe to receive the data.
   call MPI_Bcast(layer_index,nlat_lyr*nlon_lyr,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)

   return
end subroutine ed_nodeget_nl
!=========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine ed_nodeget_met_header()
!------------------------------------------------------------------------------------------!
!    This subroutine sends the met driver information to the nodes, which can then read    !
! the hdf5 in parallel                                                                     !
!------------------------------------------------------------------------------------------!
   use ed_node_coms, only: master_num,mynum
   use ed_max_dims, only: max_met_vars,str_len
   use met_driver_coms, only: nformats, met_names, met_nlon,   &
        met_nlat, met_dx, met_dy, met_xmin, met_ymin, met_nv,   &
        met_vars, met_frq, met_interp, ed_met_driver_db, no_ll,  &
        metname_len,metvars_len

   implicit none
   include 'mpif.h'
   integer             :: ierr, nsize,f,v

   

!----- First I get the scalars ------------------------------------------------------------!
   call MPI_Bcast (nformats,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)

   nsize=nformats*max_met_vars

!----- Allocate the vectors and matrices --------------------------------------------------!
   allocate(met_names(nformats))
   allocate(met_nlon(nformats))
   allocate(met_nlat(nformats))
   allocate(met_dx(nformats))
   allocate(met_dy(nformats))
   allocate(met_xmin(nformats))
   allocate(met_ymin(nformats))
   allocate(met_nv(nformats))
   allocate(met_vars(nformats, max_met_vars))
   allocate(met_frq(nformats, max_met_vars))
   allocate(met_interp(nformats, max_met_vars))
   allocate(no_ll(nformats))
       
!------------------------------------------------------------------------------------------!
!   Here I need a MPI Barrier. I don't want the master sending information before the      !
! variables are allocated in this node.                                                    !
!------------------------------------------------------------------------------------------!
   call MPI_Barrier(MPI_COMM_WORLD,ierr)

   do f=1,nformats
     call MPI_Bcast(met_names(f),metname_len,MPI_CHARACTER,master_num,MPI_COMM_WORLD,ierr)
   end do
   
   call MPI_Bcast(met_nlon,nformats,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(met_nlat,nformats,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(met_dx,nformats,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(met_dy,nformats,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(met_xmin,nformats,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(met_ymin,nformats,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(met_nv,nformats,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(no_ll, nformats, MPI_LOGICAL, master_num, MPI_COMM_WORLD, ierr)
   
   do f=1,nformats
      do v=1,max_met_vars
         call MPI_Bcast(met_vars(f,v),metvars_len,MPI_CHARACTER,master_num,MPI_COMM_WORLD,ierr)
      end do
   end do
   
   call MPI_Bcast(met_frq,nsize,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(met_interp,nsize,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)

   return
end subroutine ed_nodeget_met_header
!==========================================================================================!
!==========================================================================================!


!==========================================================================================!
!==========================================================================================!
subroutine ed_nodeget_poly_dims
   use ed_state_vars, only: gdpy,py_off
   use ed_node_coms, only: master_num,nmachs,mynum
   use grid_coms, only: ngrids
   implicit none
   include 'mpif.h'
   integer :: ierr
   integer :: ifm,nm
  
   do ifm=1,ngrids
      do nm=1,nmachs
         call MPI_Bcast(gdpy(nm,ifm),1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
         call MPI_Bcast(py_off(nm,ifm),1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
      end do
   end do
   return
end subroutine ed_nodeget_poly_dims
!==========================================================================================!
!==========================================================================================!


!==========================================================================================!
!==========================================================================================!
subroutine ed_nodeget_worklist_info

  use ed_max_dims, only: maxmach
  use grid_coms,  only: ngrids
  use ed_work_vars,  only: work_v              & ! intent(inout)
                         , ed_alloc_work_vec   & ! subroutine
                         , ed_nullify_work_vec ! ! subroutine
  use ed_node_coms,  only: mynum,nmachs,master_num
  use ed_state_vars, only: gdpy

  implicit none
  include 'mpif.h'
  
  integer,dimension(MPI_STATUS_SIZE) :: status
  integer :: ierr
  integer :: npolygons
  integer :: mpiid
  integer :: ifm
  
  ! Allocate the work vectors
  allocate(work_v(ngrids))

  do ifm=1,ngrids
     
     npolygons = gdpy(mynum,ifm)
     
     ! Allocate the work vectors
     call ed_nullify_work_vec(work_v(ifm))
     call ed_alloc_work_vec(work_v(ifm),npolygons)

     mpiid=maxmach*(ifm-1)+mynum

     call MPI_Recv(work_v(ifm)%glon(1:npolygons),npolygons, &
          MPI_REAL,0,1300000+mpiid,MPI_COMM_WORLD,status,ierr)
     
     call MPI_Recv(work_v(ifm)%glat(1:npolygons),npolygons, &
          MPI_REAL,0,1400000+mpiid,MPI_COMM_WORLD,status,ierr)
     
     call MPI_Recv(work_v(ifm)%landfrac(1:npolygons),npolygons, &
          MPI_REAL,0,1500000+mpiid,MPI_COMM_WORLD,status,ierr)
     
     call MPI_Recv(work_v(ifm)%ntext(1:npolygons),npolygons, &
          MPI_INTEGER,0,1600000+mpiid,MPI_COMM_WORLD,status,ierr)

     call MPI_Recv(work_v(ifm)%xid(1:npolygons),npolygons, &
          MPI_INTEGER,0,1700000+mpiid,MPI_COMM_WORLD,status,ierr)
     
     call MPI_Recv(work_v(ifm)%yid(1:npolygons),npolygons, &
          MPI_INTEGER,0,1800000+mpiid,MPI_COMM_WORLD,status,ierr)
  end do

  return
end subroutine ed_nodeget_worklist_info
!==========================================================================================!
!==========================================================================================!



!==========================================================================================!
!==========================================================================================!
subroutine mk_2_buff(a,b,n1,n2,m1,m2,i1,i2,j1,j2)
  implicit none
  integer :: n1,n2,m1,m2,i1,i2,j1,j2
  real :: a(n1,n2),b(m1,m2)
     
     b(1:m1,1:m2)=a(i1:i2,j1:j2)

  return
end subroutine mk_2_buff
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine mk_2p_buff(a,b,n1,n2,n3,m1,m2,m3,i1,i2,j1,j2)
  implicit none
  integer :: n1,n2,n3,m1,m2,m3,i1,i2,j1,j2
  real :: a(n1,n2,n3),b(m1,m2,m3)
     
     b(1:m1,1:m2,1:m3)=a(i1:i2,j1:j2,1:n3)

  return
end subroutine mk_2p_buff
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine mk_3_buff(a,b,n1,n2,n3,m1,m2,m3,i1,i2,j1,j2)
  implicit none
  integer :: n1,n2,n3,m1,m2,m3,i1,i2,j1,j2
  real :: a(n1,n2,n3),b(m1,m2,m3)
     
     b(1:m1,1:m2,1:m3)=a(1:n1,i1:i2,j1:j2)

  return
end subroutine mk_3_buff
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine mk_4_buff(a,b,n1,n2,n3,n4,m1,m2,m3,m4,i1,i2,j1,j2)
  implicit none
  integer :: n1,n2,n3,n4,m1,m2,m3,m4,i1,i2,j1,j2
  real :: a(n1,n2,n3,m4),b(m1,m2,m3,m4)
     
     b(1:m1,1:m2,1:m3,1:m4)=a(1:n1,i1:i2,j1:j2,1:n4)

  return
end subroutine mk_4_buff
!==========================================================================================!
!==========================================================================================!
