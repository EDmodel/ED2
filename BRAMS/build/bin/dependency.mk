# DO NOT DELETE THIS LINE - used by make depend
cyclic_mod.o: grid_dims.mod
rbnd.o: catt_start.mod mem_basic.mod mem_grid.mod mem_scratch.mod mem_tend.mod
rbnd.o: mem_turb.mod node_mod.mod ref_sounding.mod therm_lib.mod var_tables.mod
rbnd_adap.o: mem_grid.mod ref_sounding.mod
dry_dep.o: extras.mod leaf_coms.mod mem_basic.mod mem_grid.mod mem_leaf.mod
dry_dep.o: mem_micro.mod mem_scalar.mod mem_scratch.mod mem_turb.mod
dry_dep.o: rconstants.mod
emission_source_map.o: extras.mod mem_basic.mod mem_grid.mod
emission_source_map.o: mem_grid_dim_defs.mod mem_scalar.mod mem_scratch.mod
extra.o: var_tables.mod
plumerise_vector.o: extras.mod mem_basic.mod mem_grid.mod mem_scalar.mod
plumerise_vector.o: node_mod.mod rconstants.mod therm_lib.mod
coriolis.o: mem_basic.mod mem_grid.mod mem_scratch.mod mem_tend.mod
coriolis.o: rconstants.mod ref_sounding.mod
local_proc.o: io_params.mod mem_grid.mod node_mod.mod rconstants.mod
local_proc.o: ref_sounding.mod rpara.mod
mod_GhostBlock.o: mod_ghostblockpartition.mod
mod_advect_kit.o: mem_basic.mod mem_grid.mod mem_tend.mod mod_ghostblock.mod
mod_advect_kit.o: mod_ghostblockpartition.mod node_mod.mod var_tables.mod
model.o: advect_kit.mod catt_start.mod dtset.mod grid_dims.mod io_params.mod
model.o: mem_grid.mod mem_leaf.mod node_mod.mod rpara.mod
modsched.o: mem_basic.mod mem_grid.mod mem_scratch.mod
raco.o: mem_basic.mod mem_grid.mod mem_scratch.mod mem_tend.mod node_mod.mod
raco.o: rconstants.mod therm_lib.mod
raco_adap.o: mem_grid.mod mem_scratch.mod node_mod.mod rconstants.mod
radvc.o: mem_basic.mod mem_grid.mod mem_scratch.mod mem_tend.mod var_tables.mod
rams_master.o: catt_start.mod dtset.mod emission_source_map.mod grid_dims.mod
rams_master.o: io_params.mod mem_cuparm.mod mem_emiss.mod mem_grid.mod
rams_master.o: mem_leaf.mod mem_mass.mod mem_oda.mod mem_radiate.mod
rams_master.o: mem_varinit.mod node_mod.mod rpara.mod teb_spm_start.mod
ref_sounding.o: grid_dims.mod
rnode.o: advect_kit.mod catt_start.mod dtset.mod grid_dims.mod io_params.mod
rnode.o: mem_aerad.mod mem_cuparm.mod mem_grid.mod mem_leaf.mod mem_oda.mod
rnode.o: mem_radiate.mod node_mod.mod var_tables.mod
rthrm.o: mem_basic.mod mem_grid.mod mem_micro.mod mem_scratch.mod micphys.mod
rthrm.o: node_mod.mod rconstants.mod therm_lib.mod
rtimh.o: advect_kit.mod catt_start.mod emission_source_map.mod mem_all.mod
rtimh.o: mem_basic.mod mem_cuparm.mod mem_emiss.mod mem_grid.mod mem_leaf.mod
rtimh.o: mem_mass.mod mem_oda.mod mem_scalar.mod mem_turb.mod mem_varinit.mod
rtimh.o: node_mod.mod teb_spm_start.mod therm_lib.mod
rtimi.o: mem_basic.mod mem_grid.mod mem_scratch.mod mem_tend.mod node_mod.mod
rtimi.o: var_tables.mod
cu_read.o: mem_basic.mod mem_cuparm.mod mem_grid.mod
grell_coms.o: grid_dims.mod
grell_cupar_aux.o: grell_coms.mod mem_ensemble.mod mem_scratch_grell.mod
grell_cupar_aux.o: rconstants.mod therm_lib.mod
grell_cupar_downdraft.o: rconstants.mod therm_lib.mod
grell_cupar_driver.o: catt_start.mod extras.mod grell_coms.mod io_params.mod
grell_cupar_driver.o: mem_basic.mod mem_cuparm.mod mem_ensemble.mod mem_grid.mod
grell_cupar_driver.o: mem_mass.mod mem_micro.mod mem_scalar.mod mem_scratch.mod
grell_cupar_driver.o: mem_scratch_grell.mod mem_tend.mod mem_turb.mod
grell_cupar_driver.o: micphys.mod node_mod.mod therm_lib.mod
grell_cupar_dynamic.o: grell_coms.mod mem_ensemble.mod mem_scratch_grell.mod
grell_cupar_dynamic.o: rconstants.mod
grell_cupar_ensemble.o: rconstants.mod
grell_cupar_environment.o: grell_coms.mod rconstants.mod therm_lib.mod
grell_cupar_feedback.o: mem_ensemble.mod mem_scratch_grell.mod rconstants.mod
grell_cupar_static.o: mem_ensemble.mod mem_scratch_grell.mod rconstants.mod
grell_cupar_updraft.o: mem_cuparm.mod rconstants.mod therm_lib.mod
grell_extras_catt.o: grell_coms.mod mem_basic.mod mem_ensemble.mod mem_grid.mod
grell_extras_catt.o: mem_scalar.mod mem_scratch.mod mem_scratch_grell.mod
grell_extras_catt.o: mem_tconv.mod rconstants.mod
kuo_cupar_driver.o: conv_coms.mod mem_basic.mod mem_cuparm.mod mem_grid.mod
kuo_cupar_driver.o: mem_scratch.mod mem_tend.mod node_mod.mod rconstants.mod
kuo_cupar_driver.o: therm_lib.mod
mem_cuparm.o: grid_dims.mod var_tables.mod
rconv_driver.o: mem_basic.mod mem_cuparm.mod mem_grid.mod mem_scratch.mod
rconv_driver.o: mem_tend.mod mem_turb.mod node_mod.mod
shcu_vars_const.o: conv_coms.mod grid_dims.mod
souza_cupar_driver.o: conv_coms.mod mem_basic.mod mem_cuparm.mod mem_grid.mod
souza_cupar_driver.o: mem_micro.mod mem_scratch.mod mem_tend.mod mem_turb.mod
souza_cupar_driver.o: node_mod.mod shcu_vars_const.mod therm_lib.mod
edcp_driver.o: consts_coms.mod ed_misc_coms.mod ed_node_coms.mod
edcp_driver.o: ed_state_vars.mod ed_work_vars.mod grid_coms.mod io_params.mod
edcp_driver.o: mem_edcp.mod rk4_coms.mod soil_coms.mod
edcp_init.o: ed_max_dims.mod ed_node_coms.mod ed_para_coms.mod ed_state_vars.mod
edcp_init.o: ed_work_vars.mod grid_coms.mod mem_grid.mod mem_leaf.mod
edcp_init.o: node_mod.mod rpara.mod soil_coms.mod
edcp_load_namelist.o: canopy_air_coms.mod canopy_radiation_coms.mod
edcp_load_namelist.o: decomp_coms.mod disturb_coms.mod ed_max_dims.mod
edcp_load_namelist.o: ed_misc_coms.mod grid_coms.mod grid_dims.mod io_params.mod
edcp_load_namelist.o: leaf_coms.mod mem_edcp.mod mem_grid.mod mem_leaf.mod
edcp_load_namelist.o: mem_polygons.mod mem_radiate.mod met_driver_coms.mod
edcp_load_namelist.o: optimiz_coms.mod pft_coms.mod phenology_coms.mod
edcp_load_namelist.o: physiology_coms.mod rk4_coms.mod soil_coms.mod
edcp_met.o: ed_max_dims.mod ed_misc_coms.mod ed_node_coms.mod ed_state_vars.mod
edcp_met.o: leaf_coms.mod mem_basic.mod mem_cuparm.mod mem_edcp.mod mem_grid.mod
edcp_met.o: mem_leaf.mod mem_micro.mod mem_radiate.mod mem_turb.mod
edcp_met.o: met_driver_coms.mod micphys.mod node_mod.mod rconstants.mod
edcp_met.o: soil_coms.mod therm_lib.mod
edcp_met_init.o: ed_state_vars.mod ed_therm_lib.mod grid_coms.mod leaf_coms.mod
edcp_met_init.o: mem_grid.mod mem_leaf.mod mem_radiate.mod rconstants.mod
edcp_met_init.o: soil_coms.mod therm_lib.mod therm_lib8.mod
edcp_model.o: consts_coms.mod ed_max_dims.mod ed_misc_coms.mod ed_node_coms.mod
edcp_model.o: ed_state_vars.mod grid_coms.mod grid_dims.mod io_params.mod
edcp_model.o: mem_edcp.mod mem_grid.mod mem_polygons.mod rk4_coms.mod
edcp_model.o: rk4_driver.mod
edcp_mpiutils.o: canopy_air_coms.mod canopy_radiation_coms.mod decomp_coms.mod
edcp_mpiutils.o: disturb_coms.mod ed_max_dims.mod ed_misc_coms.mod grid_coms.mod
edcp_mpiutils.o: mem_edcp.mod mem_polygons.mod met_driver_coms.mod
edcp_mpiutils.o: optimiz_coms.mod pft_coms.mod phenology_coms.mod
edcp_mpiutils.o: physiology_coms.mod rk4_coms.mod soil_coms.mod
edcp_para_init.o: ed_node_coms.mod ed_work_vars.mod io_params.mod mem_grid.mod
edcp_para_init.o: mem_leaf.mod mem_polygons.mod node_mod.mod soil_coms.mod
edcp_water.o: canopy_air_coms.mod canopy_struct_dynamics.mod consts_coms.mod
edcp_water.o: io_params.mod leaf_coms.mod mem_basic.mod mem_cuparm.mod
edcp_water.o: mem_edcp.mod mem_grid.mod mem_leaf.mod mem_micro.mod
edcp_water.o: mem_radiate.mod met_driver_coms.mod node_mod.mod rk4_coms.mod
edcp_water.o: therm_lib.mod
mem_edcp.o: var_tables.mod
cond_read.o: mem_grid.mod mem_varinit.mod
cond_update.o: an_header.mod grid_struct.mod mem_basic.mod mem_grid.mod
cond_update.o: mem_varinit.mod rconstants.mod var_tables.mod
mem_oda.o: var_tables.mod
nud_analysis.o: mem_basic.mod mem_grid.mod mem_scratch.mod mem_tend.mod
nud_analysis.o: mem_varinit.mod node_mod.mod
nud_read.o: mem_grid.mod mem_varinit.mod
nud_update.o: an_header.mod grid_struct.mod mem_aerad.mod mem_basic.mod
nud_update.o: mem_grid.mod mem_varinit.mod rconstants.mod var_tables.mod
oda_krig.o: mem_oda.mod
oda_nudge.o: io_params.mod mem_basic.mod mem_grid.mod mem_oda.mod
oda_nudge.o: mem_scratch.mod mem_tend.mod node_mod.mod
oda_proc_obs.o: mem_grid.mod mem_oda.mod rconstants.mod therm_lib.mod
oda_read.o: mem_grid.mod mem_oda.mod
oda_sta_count.o: mem_grid.mod mem_oda.mod obs_input.mod
oda_sta_input.o: mem_grid.mod mem_oda.mod obs_input.mod
read_ralph.o: obs_input.mod rconstants.mod therm_lib.mod
varf_read.o: mem_grid.mod mem_varinit.mod
varf_update.o: mem_basic.mod mem_grid.mod mem_leaf.mod mem_scratch.mod
varf_update.o: mem_varinit.mod rconstants.mod ref_sounding.mod therm_lib.mod
adap_init.o: mem_leaf.mod
gridset.o: grid_dims.mod mem_grid.mod rconstants.mod
rams_grid.o: mem_grid.mod node_mod.mod rconstants.mod
rdint.o: catt_start.mod domain_decomp.mod emission_source_map.mod grell_coms.mod
rdint.o: grid_dims.mod io_params.mod isan_coms.mod leaf_coms.mod mem_basic.mod
rdint.o: mem_cuparm.mod mem_emiss.mod mem_gaspart.mod mem_grid.mod mem_leaf.mod
rdint.o: mem_mass.mod mem_micro.mod mem_oda.mod mem_radiate.mod mem_scalar.mod
rdint.o: mem_scratch.mod mem_soil_moisture.mod mem_teb.mod mem_teb_common.mod
rdint.o: mem_turb.mod mem_varinit.mod micphys.mod node_mod.mod plume_utils.mod
rdint.o: ref_sounding.mod teb_spm_start.mod teb_vars_const.mod therm_lib.mod
rdint.o: therm_lib8.mod turb_coms.mod var_tables.mod
rhhi.o: mem_basic.mod mem_grid.mod mem_scratch.mod rconstants.mod
rhhi.o: ref_sounding.mod therm_lib.mod
rinit.o: io_params.mod mem_basic.mod mem_grid.mod mem_micro.mod mem_scratch.mod
rinit.o: mem_turb.mod micphys.mod node_mod.mod rconstants.mod ref_sounding.mod
error_mess.o: node_mod.mod
inithis.o: an_header.mod grid_dims.mod io_params.mod leaf_coms.mod mem_aerad.mod
inithis.o: mem_basic.mod mem_cuparm.mod mem_grid.mod mem_leaf.mod
inithis.o: mem_scratch.mod rconstants.mod ref_sounding.mod therm_lib.mod
inithis.o: var_tables.mod
io_params.o: grid_dims.mod
opspec.o: catt_start.mod grell_coms.mod io_params.mod leaf_coms.mod
opspec.o: mem_basic.mod mem_cuparm.mod mem_emiss.mod mem_grid.mod mem_leaf.mod
opspec.o: mem_mass.mod mem_radiate.mod mem_turb.mod mem_varinit.mod micphys.mod
opspec.o: teb_spm_start.mod therm_lib.mod
rams_read_header.o: an_header.mod
ranlavg.o: io_params.mod mem_basic.mod mem_grid.mod mem_scratch.mod mem_turb.mod
ranlavg.o: node_mod.mod var_tables.mod
rcio.o: grell_coms.mod leaf_coms.mod mem_all.mod mem_mass.mod therm_lib.mod
rcio.o: turb_coms.mod
recycle.o: grid_dims.mod io_params.mod mem_aerad.mod mem_cuparm.mod mem_grid.mod
recycle.o: mem_leaf.mod mem_scratch.mod var_tables.mod
rhdf5.o: an_header.mod grid_dims.mod  io_params.mod mem_aerad.mod
rhdf5.o: mem_cuparm.mod mem_grid.mod var_tables.mod
rio.o: an_header.mod grid_dims.mod io_params.mod mem_aerad.mod mem_basic.mod
rio.o: mem_cuparm.mod mem_grid.mod mem_scratch.mod mem_turb.mod rconstants.mod
rio.o: ref_sounding.mod therm_lib.mod var_tables.mod
rname.o: catt_start.mod domain_decomp.mod emission_source_map.mod grell_coms.mod
rname.o: mem_all.mod mem_mass.mod mem_soil_moisture.mod plume_utils.mod
rname.o: teb_spm_start.mod therm_lib.mod turb_coms.mod
rprnt.o: io_params.mod leaf_coms.mod mem_all.mod mem_basic.mod mem_grid.mod
rprnt.o: mem_leaf.mod mem_scratch.mod mem_turb.mod rconstants.mod
rprnt.o: ref_sounding.mod therm_lib.mod var_tables.mod
aobj.o: isan_coms.mod rconstants.mod
asgen.o: io_params.mod isan_coms.mod mem_grid.mod
asnc.o: isan_coms.mod  rconstants.mod
asti.o: isan_coms.mod mem_grid.mod rconstants.mod therm_lib.mod
asti2.o: isan_coms.mod rconstants.mod therm_lib.mod
astp.o: isan_coms.mod rconstants.mod therm_lib.mod
avarf.o: isan_coms.mod mem_grid.mod rconstants.mod therm_lib.mod
file_inv.o: isan_coms.mod
first_rams.o: an_header.mod isan_coms.mod mem_grid.mod mem_scratch.mod
first_rams.o: rconstants.mod therm_lib.mod
isan_io.o: isan_coms.mod
refstate.o: rconstants.mod therm_lib.mod
v_interps.o: isan_coms.mod rconstants.mod therm_lib.mod
dateutils.o: rconstants.mod
getvar.o: an_header.mod
great_circle.o: rconstants.mod
hdf5_utils.o: hdf5_coms.mod
map_proj.o: rconstants.mod
numutils.o: rconstants.mod therm_lib.mod
polarst.o: rconstants.mod
therm_lib.o: rconstants.mod
therm_lib8.o: rconstants.mod therm_lib.mod
varutils.o: mem_aerad.mod mem_cuparm.mod mem_grid.mod node_mod.mod
mem_mass.o: grid_dims.mod var_tables.mod
rexev.o: mem_basic.mod mem_grid.mod mem_mass.mod mem_scratch.mod mem_tend.mod
rexev.o: rconstants.mod therm_lib.mod
rmass.o: mem_grid.mod mem_mass.mod mem_scratch.mod mem_scratch_grell.mod
rmass.o: mem_turb.mod
dealloc.o: catt_start.mod mem_aerad.mod mem_all.mod mem_ensemble.mod
dealloc.o: mem_gaspart.mod mem_globaer.mod mem_globrad.mod mem_mass.mod
dealloc.o: mem_opt.mod mem_scratch1_grell.mod mem_scratch_grell.mod mem_teb.mod
dealloc.o: mem_teb_common.mod mem_tend.mod teb_spm_start.mod
hdf5_coms.o: 
mem_all.o: io_params.mod mem_basic.mod mem_cuparm.mod mem_grid.mod mem_leaf.mod
mem_all.o: mem_micro.mod mem_nestb.mod mem_oda.mod mem_radiate.mod
mem_all.o: mem_scalar.mod mem_scratch.mod mem_scratch1.mod mem_tend.mod
mem_all.o: mem_turb.mod mem_varinit.mod micphys.mod ref_sounding.mod
mem_all.o: var_tables.mod
mem_basic.o: grid_dims.mod var_tables.mod
mem_grid.o: grid_dims.mod var_tables.mod
mem_scalar.o: var_tables.mod
mem_scratch.o: catt_start.mod grid_dims.mod mem_aerad.mod mem_radiate.mod
mem_scratch1_brams.o: var_tables.mod
mem_tend.o: mem_basic.mod mem_emiss.mod mem_gaspart.mod mem_micro.mod
mem_tend.o: mem_scalar.mod mem_turb.mod teb_spm_start.mod
mem_varinit.o: grid_dims.mod var_tables.mod
rams_mem_alloc.o: catt_start.mod extras.mod grell_coms.mod io_params.mod
rams_mem_alloc.o: leaf_coms.mod machine_arq.mod mem_aerad.mod mem_all.mod
rams_mem_alloc.o: mem_carma.mod mem_emiss.mod mem_ensemble.mod mem_gaspart.mod
rams_mem_alloc.o: mem_globaer.mod mem_globrad.mod mem_grell_param.mod
rams_mem_alloc.o: mem_grid_dim_defs.mod mem_mass.mod mem_opt.mod mem_scalar.mod
rams_mem_alloc.o: mem_scratch1_grell.mod mem_scratch2_grell.mod
rams_mem_alloc.o: mem_scratch2_grell_sh.mod mem_scratch3_grell.mod
rams_mem_alloc.o: mem_scratch3_grell_sh.mod mem_scratch_grell.mod mem_teb.mod
rams_mem_alloc.o: mem_teb_common.mod mem_turb_scalar.mod node_mod.mod
rams_mem_alloc.o: teb_spm_start.mod teb_vars_const.mod turb_coms.mod
vtab_fill.o: io_params.mod var_tables.mod
mem_micro.o: micphys.mod therm_lib.mod var_tables.mod
mic_coll.o: micphys.mod micro_coms.mod rconstants.mod therm_lib.mod
mic_driv.o: grid_dims.mod mem_basic.mod mem_grid.mod mem_micro.mod
mic_driv.o: mem_scratch.mod micphys.mod micro_coms.mod node_mod.mod
mic_driv.o: rconstants.mod therm_lib.mod
mic_gamma.o: therm_lib.mod
mic_init.o: mem_grid.mod mem_radiate.mod micphys.mod micro_coms.mod node_mod.mod
mic_init.o: rconstants.mod therm_lib.mod
mic_misc.o: mem_basic.mod mem_grid.mod mem_micro.mod mem_scratch.mod micphys.mod
mic_misc.o: micro_coms.mod rconstants.mod therm_lib.mod
mic_nuc.o: micphys.mod micro_coms.mod rconstants.mod therm_lib.mod
mic_tabs.o: micphys.mod micro_coms.mod rconstants.mod
mic_vap.o: micphys.mod micro_coms.mod rconstants.mod therm_lib.mod
micphys.o: grid_dims.mod
micro_coms.o: micphys.mod rconstants.mod
geodat.o: io_params.mod mem_grid.mod mem_leaf.mod rconstants.mod
geodat.o: teb_spm_start.mod
landuse_input.o: hdf5_utils.mod io_params.mod leaf_coms.mod mem_mksfc.mod
landuse_input.o: rconstants.mod
mem_mksfc.o: teb_spm_start.mod
mksfc_driver.o: io_params.mod mem_grid.mod mem_mksfc.mod teb_spm_start.mod
mksfc_fuso.o: io_params.mod mem_emiss.mod mem_gaspart.mod mem_grid.mod
mksfc_fuso.o: mem_mksfc.mod mem_teb.mod teb_vars_const.mod
mksfc_ndvi.o: io_params.mod mem_grid.mod mem_leaf.mod mem_mksfc.mod
mksfc_sfc.o: io_params.mod mem_grid.mod mem_leaf.mod mem_mksfc.mod
mksfc_sst.o: io_params.mod mem_grid.mod mem_leaf.mod mem_mksfc.mod
mksfc_top.o: io_params.mod mem_grid.mod mem_mksfc.mod
ndvi_read.o: io_params.mod mem_grid.mod mem_leaf.mod
nest_geosst.o: io_params.mod leaf_coms.mod mem_basic.mod mem_grid.mod
nest_geosst.o: mem_leaf.mod mem_mksfc.mod mem_scratch.mod mem_soil_moisture.mod
nest_init_aux.o: mem_basic.mod mem_grid.mod mem_leaf.mod mem_scratch.mod
sst_read.o: io_params.mod mem_grid.mod mem_leaf.mod
mpass_cyclic.o: cyclic_mod.mod grid_dims.mod mem_aerad.mod mem_basic.mod
mpass_cyclic.o: mem_cuparm.mod mem_grid.mod mem_scratch.mod node_mod.mod
mpass_cyclic.o: var_tables.mod
mpass_dtl.o: mem_grid.mod node_mod.mod rpara.mod
mpass_feed.o: grid_dims.mod mem_basic.mod mem_grid.mod mem_scratch1.mod
mpass_feed.o: node_mod.mod var_tables.mod
mpass_full.o: grid_dims.mod io_params.mod mem_aerad.mod mem_cuparm.mod
mpass_full.o: mem_grid.mod mem_scratch.mod mem_varinit.mod node_mod.mod
mpass_full.o: rpara.mod var_tables.mod
mpass_init.o: catt_start.mod cyclic_mod.mod emission_source_map.mod
mpass_init.o: grell_coms.mod grid_dims.mod leaf_coms.mod mem_all.mod
mpass_init.o: mem_cuparm.mod mem_emiss.mod mem_grid.mod mem_mass.mod micphys.mod
mpass_init.o: node_mod.mod plume_utils.mod ref_sounding.mod rpara.mod
mpass_init.o: teb_spm_start.mod teb_vars_const.mod therm_lib.mod turb_coms.mod
mpass_lbc.o: grid_dims.mod mem_aerad.mod mem_cuparm.mod mem_grid.mod
mpass_lbc.o: mem_scratch.mod node_mod.mod var_tables.mod
mpass_nest.o: grid_dims.mod mem_basic.mod mem_grid.mod mem_nestb.mod
mpass_nest.o: mem_scratch.mod node_mod.mod var_tables.mod
mpass_oda.o: grid_dims.mod mem_oda.mod node_mod.mod rpara.mod
mpass_st.o: grid_dims.mod mem_basic.mod mem_grid.mod mem_scratch.mod
mpass_st.o: node_mod.mod
node_mod.o: grid_dims.mod
par_decomp.o: cyclic_mod.mod domain_decomp.mod grid_dims.mod mem_grid.mod
par_decomp.o: rpara.mod
para_init.o: grid_dims.mod mem_aerad.mod mem_basic.mod mem_cuparm.mod
para_init.o: mem_grid.mod mem_scratch.mod node_mod.mod rpara.mod var_tables.mod
paral.o: grid_dims.mod mem_aerad.mod mem_cuparm.mod mem_grid.mod mem_scratch.mod
paral.o: node_mod.mod rpara.mod var_tables.mod
rnest_par.o: mem_grid.mod
rpara.o: grid_dims.mod
hemi2.o: grid_dims.mod mem_basic.mod mem_grid.mod var_tables.mod
mem_nestb.o: var_tables.mod
nest_drivers.o: mem_basic.mod mem_grid.mod mem_nestb.mod mem_scratch.mod
nest_drivers.o: mem_tend.mod node_mod.mod var_tables.mod
nest_feed.o: mem_grid.mod
nest_intrp.o: grid_dims.mod mem_basic.mod mem_grid.mod mem_nestb.mod
nest_intrp.o: mem_scratch.mod rconstants.mod ref_sounding.mod
nest_move.o: mem_basic.mod mem_grid.mod mem_leaf.mod mem_scratch.mod
nest_move.o: mem_tend.mod mem_turb.mod var_tables.mod
cup_dn.o: rconstants.mod
cup_env.o: rconstants.mod therm_lib.mod
cup_grell2.o: mem_grell_param.mod mem_scratch2_grell.mod mem_scratch3_grell.mod
cup_grell2.o: rconstants.mod
cup_grell2_shcu.o: mem_grell_param.mod mem_scratch2_grell_sh.mod
cup_grell2_shcu.o: mem_scratch3_grell_sh.mod rconstants.mod
cup_up.o: rconstants.mod
mem_grell_param2.o: grell_coms.mod
mem_scratch2_grell.o: mem_grell_param.mod node_mod.mod
mem_scratch2_grell_sh.o: mem_grell_param.mod node_mod.mod
mem_scratch3_grell.o: mem_grell_param.mod
mem_scratch3_grell_sh.o: mem_grell_param.mod
old_grell_cupar_driver.o: grell_coms.mod io_params.mod mem_basic.mod
old_grell_cupar_driver.o: mem_cuparm.mod mem_grell_param.mod mem_grid.mod
old_grell_cupar_driver.o: mem_leaf.mod mem_mass.mod mem_micro.mod mem_scalar.mod
old_grell_cupar_driver.o: mem_scratch.mod mem_scratch1_grell.mod
old_grell_cupar_driver.o: mem_scratch2_grell.mod mem_scratch2_grell_sh.mod
old_grell_cupar_driver.o: mem_scratch3_grell.mod mem_scratch3_grell_sh.mod
old_grell_cupar_driver.o: mem_tend.mod mem_turb.mod node_mod.mod rconstants.mod
old_grell_cupar_driver.o: therm_lib.mod
harr_coms.o: mem_harr.mod
harr_rad.o: harr_coms.mod mem_harr.mod rconstants.mod
harr_raddriv.o: harr_coms.mod mem_grid.mod mem_harr.mod mem_leaf.mod micphys.mod
harr_raddriv.o: rconstants.mod therm_lib.mod
harr_radinit.o: harr_coms.mod mem_cuparm.mod mem_grid.mod mem_harr.mod
harr_radinit.o: mem_radiate.mod micphys.mod
mem_aerad.o: mem_grid_dim_defs.mod
mem_carma.o: grid_dims.mod mem_aerad.mod mem_globrad.mod
mem_globaer.o: mem_aerad.mod
mem_globrad.o: mem_aerad.mod rconstants.mod
mem_mclat.o: rconstants.mod
mem_radiate.o: var_tables.mod
rad_carma.o: catt_start.mod grid_dims.mod mem_aerad.mod mem_carma.mod
rad_carma.o: mem_globaer.mod mem_globrad.mod mem_grid.mod mem_radiate.mod
rad_carma.o: node_mod.mod rconstants.mod
rad_ccmp.o: rconstants.mod
rad_driv.o: catt_start.mod mem_basic.mod mem_cuparm.mod mem_grid.mod
rad_driv.o: mem_harr.mod mem_leaf.mod mem_mclat.mod mem_micro.mod
rad_driv.o: mem_radiate.mod mem_scalar.mod mem_scratch.mod mem_teb_common.mod
rad_driv.o: mem_tend.mod micphys.mod rad_carma.mod rconstants.mod
rad_driv.o: teb_spm_start.mod therm_lib.mod
rad_mclat.o: harr_coms.mod mem_grid.mod mem_mclat.mod mem_radiate.mod
rad_mclat.o: rconstants.mod
mem_soil_moisture.o: leaf_coms.mod
soil_moisture_init.o: io_params.mod leaf_coms.mod mem_grid.mod mem_leaf.mod
soil_moisture_init.o: mem_soil_moisture.mod rconstants.mod
leaf3.o: io_params.mod leaf_coms.mod mem_basic.mod mem_cuparm.mod mem_grid.mod
leaf3.o: mem_leaf.mod mem_micro.mod mem_radiate.mod mem_scratch.mod mem_teb.mod
leaf3.o: mem_teb_common.mod mem_turb.mod node_mod.mod rconstants.mod
leaf3.o: teb_spm_start.mod therm_lib.mod
leaf3_can.o: catt_start.mod leaf_coms.mod rconstants.mod therm_lib.mod
leaf3_hyd.o: leaf_coms.mod mem_grid.mod mem_leaf.mod rconstants.mod
leaf3_hyd.o: therm_lib.mod
leaf3_init.o: io_params.mod leaf_coms.mod mem_grid.mod mem_leaf.mod
leaf3_init.o: rconstants.mod teb_spm_start.mod therm_lib.mod
leaf3_ocean.o: leaf_coms.mod rconstants.mod therm_lib.mod
leaf3_teb.o: mem_emiss.mod rconstants.mod teb_vars_const.mod therm_lib.mod
leaf3_temp.o: catt_start.mod leaf_coms.mod mem_all.mod mem_scratch.mod
leaf3_temp.o: mem_teb.mod mem_teb_common.mod node_mod.mod rconstants.mod
leaf3_temp.o: teb_spm_start.mod therm_lib.mod
leaf3_tw.o: catt_start.mod leaf_coms.mod mem_leaf.mod mem_scratch.mod
leaf3_tw.o: rconstants.mod therm_lib.mod
leaf3_utils.o: catt_start.mod io_params.mod leaf_coms.mod mem_grid.mod
leaf3_utils.o: mem_leaf.mod mem_radiate.mod mem_scratch.mod node_mod.mod
leaf3_utils.o: rconstants.mod teb_spm_start.mod therm_lib.mod
leaf_coms.o: grid_dims.mod mem_leaf.mod rconstants.mod therm_lib.mod
mem_leaf.o: grid_dims.mod io_params.mod var_tables.mod
ruser.o: catt_start.mod io_params.mod leaf_coms.mod mem_grid.mod mem_leaf.mod
ruser.o: rconstants.mod therm_lib.mod
urban.o: teb_vars_const.mod therm_lib.mod
urban_canopy.o: mem_basic.mod mem_grid.mod mem_scratch.mod mem_tend.mod
urban_canopy.o: mem_turb.mod node_mod.mod
gaspart.o: an_header.mod grid_dims.mod io_params.mod mem_basic.mod mem_emiss.mod
gaspart.o: mem_gaspart.mod mem_grid.mod mem_leaf.mod mem_tend.mod rconstants.mod
gaspart.o: ref_sounding.mod teb_vars_const.mod var_tables.mod
mem_gaspart.o: mem_emiss.mod var_tables.mod
mem_teb.o: var_tables.mod
mem_teb_common.o: var_tables.mod
mem_teb_vars_const.o: grid_dims.mod
ozone.o: mem_basic.mod mem_gaspart.mod mem_grid.mod mem_radiate.mod mem_tend.mod
ozone.o: ozone_const.mod rconstants.mod var_tables.mod
diffsclr.o: mem_grid.mod mem_scratch.mod
diffuse.o: ke_coms.mod mem_basic.mod mem_grid.mod mem_leaf.mod mem_mass.mod
diffuse.o: mem_micro.mod mem_opt.mod mem_scratch.mod mem_tend.mod mem_turb.mod
diffuse.o: node_mod.mod therm_lib.mod var_tables.mod
mem_turb.o: grid_dims.mod rconstants.mod var_tables.mod
mem_turb_scalar.o: grid_dims.mod var_tables.mod
rgrad.o: mem_grid.mod mem_scratch.mod
tkenn.o: leaf_coms.mod mem_grid.mod mem_scratch.mod rconstants.mod therm_lib.mod
tkenn.o: turb_coms.mod
turb_derivs.o: mem_grid.mod mem_scratch.mod mem_turb.mod rconstants.mod
turb_derivs.o: therm_lib.mod
turb_diff.o: catt_start.mod mem_cuparm.mod mem_grid.mod mem_opt.mod
turb_diff.o: mem_scratch.mod mem_turb.mod var_tables.mod
turb_diff_adap.o: mem_grid.mod mem_scratch.mod
turb_k.o: catt_start.mod ke_coms.mod mem_basic.mod mem_grid.mod mem_leaf.mod
turb_k.o: mem_mass.mod mem_micro.mod mem_scratch.mod mem_tend.mod mem_turb.mod
turb_k.o: mem_turb_scalar.mod node_mod.mod therm_lib.mod var_tables.mod
turb_ke.o: ke_coms.mod mem_grid.mod mem_scratch.mod mem_turb.mod rconstants.mod
turb_ke.o: turb_coms.mod
ed_1st.o: ed_misc_coms.mod ed_para_coms.mod ed_state_vars.mod
ed_driver.o: consts_coms.mod ed_misc_coms.mod ed_node_coms.mod ed_state_vars.mod
ed_driver.o: fuse_fiss_utils.mod grid_coms.mod soil_coms.mod
ed_met_driver.o: canopy_air_coms.mod consts_coms.mod ed_max_dims.mod
ed_met_driver.o: ed_misc_coms.mod ed_state_vars.mod grid_coms.mod hdf5_utils.mod
ed_met_driver.o: mem_polygons.mod met_driver_coms.mod therm_lib.mod
ed_model.o: consts_coms.mod disturb_coms.mod ed_misc_coms.mod ed_node_coms.mod
ed_model.o: ed_state_vars.mod grid_coms.mod mem_polygons.mod rk4_coms.mod
ed_model.o: rk4_driver.mod
canopy_struct_dynamics.o: allometry.mod canopy_air_coms.mod consts_coms.mod
canopy_struct_dynamics.o: ed_state_vars.mod met_driver_coms.mod pft_coms.mod
canopy_struct_dynamics.o: physiology_coms.mod rk4_coms.mod soil_coms.mod
disturbance.o: allometry.mod consts_coms.mod decomp_coms.mod disturb_coms.mod
disturbance.o: ed_max_dims.mod ed_misc_coms.mod ed_state_vars.mod
disturbance.o: ed_therm_lib.mod fuse_fiss_utils.mod grid_coms.mod
disturbance.o: mem_polygons.mod pft_coms.mod
euler_driver.o: canopy_air_coms.mod canopy_struct_dynamics.mod consts_coms.mod
euler_driver.o: ed_max_dims.mod ed_misc_coms.mod ed_state_vars.mod grid_coms.mod
euler_driver.o: hydrology_coms.mod met_driver_coms.mod rk4_coms.mod
euler_driver.o: rk4_driver.mod rk4_stepper.mod soil_coms.mod
events.o: allometry.mod consts_coms.mod decomp_coms.mod disturbance_utils.mod
events.o: ed_misc_coms.mod ed_state_vars.mod ed_therm_lib.mod
events.o: fuse_fiss_utils.mod grid_coms.mod pft_coms.mod therm_lib.mod
farq_leuning.o: c34constants.mod consts_coms.mod pft_coms.mod phenology_coms.mod
farq_leuning.o: physiology_coms.mod rk4_coms.mod therm_lib8.mod
fire.o: allometry.mod consts_coms.mod disturb_coms.mod ed_state_vars.mod
fire.o: grid_coms.mod pft_coms.mod soil_coms.mod
forestry.o: allometry.mod disturb_coms.mod disturbance_utils.mod ed_max_dims.mod
forestry.o: ed_state_vars.mod fuse_fiss_utils.mod grid_coms.mod
growth_balive.o: allometry.mod consts_coms.mod decomp_coms.mod ed_max_dims.mod
growth_balive.o: ed_misc_coms.mod ed_state_vars.mod ed_therm_lib.mod
growth_balive.o: grid_coms.mod mortality.mod pft_coms.mod physiology_coms.mod
heun_driver.o: canopy_air_coms.mod canopy_struct_dynamics.mod consts_coms.mod
heun_driver.o: ed_max_dims.mod ed_misc_coms.mod ed_state_vars.mod grid_coms.mod
heun_driver.o: hydrology_coms.mod met_driver_coms.mod rk4_coms.mod
heun_driver.o: rk4_driver.mod rk4_stepper.mod soil_coms.mod
lsm_hyd.o: consts_coms.mod ed_misc_coms.mod ed_node_coms.mod ed_state_vars.mod
lsm_hyd.o: grid_coms.mod hydrology_coms.mod hydrology_constants.mod pft_coms.mod
lsm_hyd.o: soil_coms.mod therm_lib.mod
mortality.o: consts_coms.mod disturb_coms.mod ed_max_dims.mod ed_misc_coms.mod
mortality.o: ed_state_vars.mod pft_coms.mod
phenology_aux.o: consts_coms.mod ed_max_dims.mod ed_state_vars.mod pft_coms.mod
phenology_aux.o: phenology_coms.mod
phenology_driv.o: allometry.mod consts_coms.mod decomp_coms.mod ed_max_dims.mod
phenology_driv.o: ed_misc_coms.mod ed_state_vars.mod ed_therm_lib.mod
phenology_driv.o: grid_coms.mod pft_coms.mod phenology_coms.mod soil_coms.mod
photosyn_driv.o: consts_coms.mod ed_max_dims.mod ed_misc_coms.mod
photosyn_driv.o: ed_state_vars.mod farq_leuning.mod met_driver_coms.mod
photosyn_driv.o: pft_coms.mod physiology_coms.mod soil_coms.mod
radiate_driver.o: allometry.mod canopy_radiation_coms.mod consts_coms.mod
radiate_driver.o: ed_max_dims.mod ed_misc_coms.mod ed_state_vars.mod
radiate_driver.o: grid_coms.mod soil_coms.mod
reproduction.o: allometry.mod consts_coms.mod decomp_coms.mod ed_max_dims.mod
reproduction.o: ed_state_vars.mod ed_therm_lib.mod fuse_fiss_utils.mod
reproduction.o: mem_polygons.mod pft_coms.mod phenology_coms.mod
rk4_derivs.o: canopy_struct_dynamics.mod consts_coms.mod ed_max_dims.mod
rk4_derivs.o: ed_misc_coms.mod ed_state_vars.mod grid_coms.mod pft_coms.mod
rk4_derivs.o: rk4_coms.mod soil_coms.mod therm_lib8.mod
rk4_driver.o: canopy_air_coms.mod canopy_struct_dynamics.mod consts_coms.mod
rk4_driver.o: ed_misc_coms.mod ed_state_vars.mod grid_coms.mod
rk4_driver.o: met_driver_coms.mod rk4_coms.mod soil_coms.mod therm_lib.mod
rk4_integ_utils.o: canopy_air_coms.mod consts_coms.mod ed_max_dims.mod
rk4_integ_utils.o: ed_misc_coms.mod ed_state_vars.mod grid_coms.mod
rk4_integ_utils.o: hydrology_coms.mod rk4_coms.mod rk4_stepper.mod soil_coms.mod
rk4_integ_utils.o: therm_lib8.mod
rk4_misc.o: allometry.mod canopy_air_coms.mod canopy_struct_dynamics.mod
rk4_misc.o: consts_coms.mod ed_max_dims.mod ed_misc_coms.mod ed_state_vars.mod
rk4_misc.o: ed_therm_lib.mod grid_coms.mod rk4_coms.mod soil_coms.mod
rk4_misc.o: therm_lib8.mod
rk4_stepper.o: consts_coms.mod ed_state_vars.mod grid_coms.mod rk4_coms.mod
rk4_stepper.o: soil_coms.mod therm_lib8.mod
soil_respiration.o: consts_coms.mod decomp_coms.mod ed_state_vars.mod
soil_respiration.o: grid_coms.mod pft_coms.mod soil_coms.mod
structural_growth.o: allometry.mod consts_coms.mod decomp_coms.mod
structural_growth.o: ed_max_dims.mod ed_misc_coms.mod ed_state_vars.mod
structural_growth.o: ed_therm_lib.mod pft_coms.mod
twostream_rad.o: canopy_radiation_coms.mod consts_coms.mod ed_max_dims.mod
twostream_rad.o: pft_coms.mod rk4_coms.mod
vegetation_dynamics.o: consts_coms.mod disturb_coms.mod disturbance_utils.mod
vegetation_dynamics.o: ed_misc_coms.mod ed_state_vars.mod fuse_fiss_utils.mod
vegetation_dynamics.o: grid_coms.mod growth_balive.mod mem_polygons.mod
ed_init.o: consts_coms.mod ed_misc_coms.mod ed_node_coms.mod ed_state_vars.mod
ed_init.o: ed_work_vars.mod grid_coms.mod phenology_coms.mod
ed_init.o: phenology_startup.mod rk4_coms.mod soil_coms.mod
ed_init_atm.o: canopy_struct_dynamics.mod consts_coms.mod ed_misc_coms.mod
ed_init_atm.o: ed_node_coms.mod ed_state_vars.mod ed_therm_lib.mod
ed_init_atm.o: fuse_fiss_utils.mod grid_coms.mod met_driver_coms.mod
ed_init_atm.o: pft_coms.mod soil_coms.mod therm_lib.mod
ed_nbg_init.o: allometry.mod consts_coms.mod ed_max_dims.mod ed_misc_coms.mod
ed_nbg_init.o: ed_state_vars.mod ed_therm_lib.mod fuse_fiss_utils.mod
ed_nbg_init.o: pft_coms.mod
ed_params.o: allometry.mod canopy_air_coms.mod canopy_radiation_coms.mod
ed_params.o: consts_coms.mod decomp_coms.mod disturb_coms.mod ed_max_dims.mod
ed_params.o: ed_misc_coms.mod fusion_fission_coms.mod grid_coms.mod
ed_params.o: hydrology_coms.mod met_driver_coms.mod pft_coms.mod
ed_params.o: phenology_coms.mod physiology_coms.mod rk4_coms.mod soil_coms.mod
ed_type_init.o: allometry.mod canopy_air_coms.mod consts_coms.mod
ed_type_init.o: ed_max_dims.mod ed_misc_coms.mod ed_state_vars.mod
ed_type_init.o: ed_therm_lib.mod grid_coms.mod pft_coms.mod phenology_coms.mod
ed_type_init.o: soil_coms.mod therm_lib.mod
init_hydro_sites.o: ed_max_dims.mod ed_misc_coms.mod ed_state_vars.mod
init_hydro_sites.o: grid_coms.mod mem_polygons.mod soil_coms.mod
landuse_init.o: consts_coms.mod disturb_coms.mod ed_max_dims.mod
landuse_init.o: ed_misc_coms.mod ed_state_vars.mod grid_coms.mod pft_coms.mod
phenology_startup.o: ed_max_dims.mod ed_misc_coms.mod ed_state_vars.mod
phenology_startup.o: grid_coms.mod phenology_coms.mod
average_utils.o: allometry.mod canopy_radiation_coms.mod consts_coms.mod
average_utils.o: ed_max_dims.mod ed_misc_coms.mod ed_state_vars.mod
average_utils.o: grid_coms.mod pft_coms.mod therm_lib.mod
ed_init_full_history.o: allometry.mod c34constants.mod consts_coms.mod
ed_init_full_history.o: ed_max_dims.mod ed_misc_coms.mod ed_node_coms.mod
ed_init_full_history.o: ed_state_vars.mod ed_therm_lib.mod
ed_init_full_history.o: fusion_fission_coms.mod grid_coms.mod 
ed_init_full_history.o: hdf5_coms.mod mem_polygons.mod pft_coms.mod
ed_init_full_history.o: soil_coms.mod therm_lib.mod
ed_load_namelist.o: canopy_air_coms.mod canopy_radiation_coms.mod
ed_load_namelist.o: consts_coms.mod decomp_coms.mod disturb_coms.mod
ed_load_namelist.o: ed_max_dims.mod ed_misc_coms.mod ed_para_coms.mod
ed_load_namelist.o: ename_coms.mod grid_coms.mod mem_polygons.mod
ed_load_namelist.o: met_driver_coms.mod optimiz_coms.mod pft_coms.mod
ed_load_namelist.o: phenology_coms.mod physiology_coms.mod rk4_coms.mod
ed_load_namelist.o: soil_coms.mod
ed_opspec.o: canopy_air_coms.mod canopy_radiation_coms.mod consts_coms.mod
ed_opspec.o: decomp_coms.mod disturb_coms.mod ed_max_dims.mod ed_misc_coms.mod
ed_opspec.o: ed_para_coms.mod grid_coms.mod mem_polygons.mod met_driver_coms.mod
ed_opspec.o: pft_coms.mod phenology_coms.mod physiology_coms.mod rk4_coms.mod
ed_opspec.o: soil_coms.mod
ed_print.o: ed_max_dims.mod ed_misc_coms.mod ed_node_coms.mod ed_state_vars.mod
ed_print.o: ed_var_tables.mod
ed_read_ed10_20_history.o: allometry.mod consts_coms.mod disturb_coms.mod
ed_read_ed10_20_history.o: ed_max_dims.mod ed_misc_coms.mod ed_state_vars.mod
ed_read_ed10_20_history.o: fuse_fiss_utils.mod grid_coms.mod mem_polygons.mod
ed_read_ed10_20_history.o: pft_coms.mod
ed_read_ed21_history.o: allometry.mod consts_coms.mod disturb_coms.mod
ed_read_ed21_history.o: ed_max_dims.mod ed_misc_coms.mod ed_state_vars.mod
ed_read_ed21_history.o: fuse_fiss_utils.mod grid_coms.mod  hdf5_coms.mod
ed_read_ed21_history.o: pft_coms.mod
ed_xml_config.o: canopy_radiation_coms.mod decomp_coms.mod disturb_coms.mod
ed_xml_config.o: ed_max_dims.mod ed_misc_coms.mod fusion_fission_coms.mod
ed_xml_config.o: grid_coms.mod hydrology_coms.mod met_driver_coms.mod
ed_xml_config.o: pft_coms.mod phenology_coms.mod physiology_coms.mod
ed_xml_config.o: rk4_coms.mod soil_coms.mod
edio.o: c34constants.mod consts_coms.mod ed_max_dims.mod ed_misc_coms.mod
edio.o: ed_node_coms.mod ed_state_vars.mod grid_coms.mod pft_coms.mod
edio.o: soil_coms.mod therm_lib.mod
h5_output.o: an_header.mod c34constants.mod ed_max_dims.mod ed_misc_coms.mod
h5_output.o: ed_node_coms.mod ed_state_vars.mod ed_var_tables.mod
h5_output.o: fusion_fission_coms.mod grid_coms.mod  hdf5_coms.mod
leaf_database.o: grid_coms.mod hdf5_utils.mod soil_coms.mod
canopy_air_coms.o: consts_coms.mod therm_lib.mod therm_lib8.mod
canopy_radiation_coms.o: ed_max_dims.mod
consts_coms.o: rconstants.mod
decomp_coms.o: ed_max_dims.mod
disturb_coms.o: ed_max_dims.mod
ed_max_dims.o: grid_dims.mod
ed_mem_alloc.o: ed_max_dims.mod ed_mem_grid_dim_defs.mod ed_misc_coms.mod
ed_mem_alloc.o: ed_node_coms.mod ed_state_vars.mod ed_work_vars.mod
ed_mem_alloc.o: grid_coms.mod mem_polygons.mod
ed_misc_coms.o: ed_max_dims.mod
ed_state_vars.o: c34constants.mod disturb_coms.mod ed_max_dims.mod
ed_state_vars.o: ed_misc_coms.mod ed_node_coms.mod ed_var_tables.mod
ed_state_vars.o: fusion_fission_coms.mod grid_coms.mod met_driver_coms.mod
ed_state_vars.o: phenology_coms.mod soil_coms.mod
ed_var_tables.o: ed_max_dims.mod
ed_work_vars.o: ed_max_dims.mod
ename_coms.o: ed_max_dims.mod
fusion_fission_coms.o: ed_max_dims.mod
grid_coms.o: ed_max_dims.mod
hdf5_coms.o: 
mem_polygons.o: ed_max_dims.mod
met_driver_coms.o: ed_max_dims.mod
optimiz_coms.o: ed_max_dims.mod
pft_coms.o: ed_max_dims.mod
phenology_coms.o: ed_max_dims.mod
physiology_coms.o: ed_max_dims.mod
rk4_coms.o: consts_coms.mod ed_max_dims.mod ed_misc_coms.mod grid_coms.mod
rk4_coms.o: soil_coms.mod therm_lib8.mod
soil_coms.o: ed_max_dims.mod grid_coms.mod leaf_coms.mod
ed_mpass_init.o: canopy_air_coms.mod canopy_radiation_coms.mod decomp_coms.mod
ed_mpass_init.o: disturb_coms.mod ed_max_dims.mod ed_misc_coms.mod
ed_mpass_init.o: ed_node_coms.mod ed_para_coms.mod ed_state_vars.mod
ed_mpass_init.o: ed_work_vars.mod grid_coms.mod mem_polygons.mod
ed_mpass_init.o: met_driver_coms.mod optimiz_coms.mod pft_coms.mod
ed_mpass_init.o: phenology_coms.mod physiology_coms.mod rk4_coms.mod
ed_mpass_init.o: soil_coms.mod
ed_node_coms.o: ed_max_dims.mod
ed_para_coms.o: ed_max_dims.mod
ed_para_init.o: ed_max_dims.mod ed_misc_coms.mod ed_node_coms.mod
ed_para_init.o: ed_para_coms.mod ed_work_vars.mod grid_coms.mod 
ed_para_init.o: hdf5_coms.mod mem_polygons.mod soil_coms.mod
allometry.o: consts_coms.mod grid_coms.mod pft_coms.mod rk4_coms.mod
allometry.o: soil_coms.mod
budget_utils.o: consts_coms.mod ed_max_dims.mod ed_misc_coms.mod
budget_utils.o: ed_state_vars.mod grid_coms.mod rk4_coms.mod soil_coms.mod
dateutils.o: consts_coms.mod
ed_filelist.o: ed_max_dims.mod
ed_grid.o: consts_coms.mod ed_max_dims.mod ed_node_coms.mod grid_coms.mod
ed_therm_lib.o: allometry.mod consts_coms.mod ed_state_vars.mod pft_coms.mod
ed_therm_lib.o: rk4_coms.mod soil_coms.mod therm_lib.mod therm_lib8.mod
fatal_error.o: ed_node_coms.mod
fuse_fiss_utils.o: allometry.mod consts_coms.mod decomp_coms.mod
fuse_fiss_utils.o: disturb_coms.mod ed_max_dims.mod ed_misc_coms.mod
fuse_fiss_utils.o: ed_node_coms.mod ed_state_vars.mod fusion_fission_coms.mod
fuse_fiss_utils.o: grid_coms.mod mem_polygons.mod pft_coms.mod soil_coms.mod
fuse_fiss_utils.o: therm_lib.mod
great_circle.o: consts_coms.mod
hdf5_utils.o: hdf5_coms.mod
invmondays.o: ed_misc_coms.mod
lapse.o: consts_coms.mod ed_state_vars.mod met_driver_coms.mod
numutils.o: consts_coms.mod therm_lib.mod
stable_cohorts.o: allometry.mod canopy_radiation_coms.mod ed_state_vars.mod
stable_cohorts.o: pft_coms.mod
therm_lib.o: consts_coms.mod
therm_lib8.o: consts_coms.mod therm_lib.mod
update_derived_props.o: allometry.mod canopy_air_coms.mod consts_coms.mod
update_derived_props.o: ed_misc_coms.mod ed_state_vars.mod ed_therm_lib.mod
update_derived_props.o: fuse_fiss_utils.mod grid_coms.mod soil_coms.mod
update_derived_props.o: therm_lib.mod
advect_kit.mod: mod_advect_kit.o
allometry.mod: allometry.o
an_header.mod: an_header.o
c34constants.mod: c34constants.o
canopy_air_coms.mod: canopy_air_coms.o
canopy_radiation_coms.mod: canopy_radiation_coms.o
canopy_struct_dynamics.mod: canopy_struct_dynamics.o
catt_start.mod: catt_start.o
consts_coms.mod: consts_coms.o
conv_coms.mod: conv_coms.o
cyclic_mod.mod: cyclic_mod.o
decomp_coms.mod: decomp_coms.o
disturb_coms.mod: disturb_coms.o
disturbance_utils.mod: disturbance.o
domain_decomp.mod: domain_decomp.o
dtset.mod: local_proc.o
ed_max_dims.mod: ed_max_dims.o
ed_mem_grid_dim_defs.mod: ed_mem_grid_dim_defs.o
ed_misc_coms.mod: ed_misc_coms.o
ed_node_coms.mod: ed_node_coms.o
ed_para_coms.mod: ed_para_coms.o
ed_state_vars.mod: ed_state_vars.o
ed_therm_lib.mod: ed_therm_lib.o
ed_var_tables.mod: ed_var_tables.o
ed_work_vars.mod: ed_work_vars.o
emission_source_map.mod: emission_source_map.o
ename_coms.mod: ename_coms.o
extras.mod: extra.o
farq_leuning.mod: farq_leuning.o
fuse_fiss_utils.mod: fuse_fiss_utils.o
fusion_fission_coms.mod: fusion_fission_coms.o
grell_coms.mod: grell_coms.o
grid_coms.mod: grid_coms.o
grid_dims.mod: grid_dims.o
grid_struct.mod: grid_struct.o
growth_balive.mod: growth_balive.o
harr_coms.mod: harr_coms.o
hdf5_coms.mod: hdf5_coms.o
hdf5_utils.mod: hdf5_utils.o
hydrology_coms.mod: hydrology_coms.o
hydrology_constants.mod: hydrology_constants.o
io_params.mod: io_params.o
isan_coms.mod: isan_coms.o
ke_coms.mod: ke_coms.o
leaf_coms.mod: leaf_coms.o
libxml2f90_interface_module.mod: libxml2f90.f90_pp.o
libxml2f90_module.mod: libxml2f90.f90_pp.o
libxml2f90_strings_module.mod: libxml2f90.f90_pp.o
ll_module.mod: libxml2f90.f90_pp.o
machine_arq.mod: machine_arq.o
mem_aerad.mod: mem_aerad.o
mem_all.mod: mem_all.o
mem_basic.mod: mem_basic.o
mem_carma.mod: mem_carma.o
mem_cuparm.mod: mem_cuparm.o
mem_edcp.mod: mem_edcp.o
mem_emiss.mod: mem_emiss.o
mem_ensemble.mod: mem_ensemble.o
mem_gaspart.mod: mem_gaspart.o
mem_globaer.mod: mem_globaer.o
mem_globrad.mod: mem_globrad.o
mem_grell_param.mod: mem_grell_param2.o
mem_grid.mod: mem_grid.o
mem_grid_dim_defs.mod: mem_grid_dim_defs.o
mem_harr.mod: mem_harr.o
mem_leaf.mod: mem_leaf.o
mem_mass.mod: mem_mass.o
mem_mclat.mod: mem_mclat.o
mem_micro.mod: mem_micro.o
mem_mksfc.mod: mem_mksfc.o
mem_nestb.mod: mem_nestb.o
mem_oda.mod: mem_oda.o
mem_opt.mod: mem_opt_scratch.o
mem_polygons.mod: mem_polygons.o
mem_radiate.mod: mem_radiate.o
mem_scalar.mod: mem_scalar.o
mem_scratch.mod: mem_scratch.o
mem_scratch1.mod: mem_scratch1_brams.o
mem_scratch1_grell.mod: mem_scratch1_grell.o
mem_scratch2_grell.mod: mem_scratch2_grell.o
mem_scratch2_grell_sh.mod: mem_scratch2_grell_sh.o
mem_scratch3_grell.mod: mem_scratch3_grell.o
mem_scratch3_grell_sh.mod: mem_scratch3_grell_sh.o
mem_scratch_grell.mod: mem_scratch_grell.o
mem_soil_moisture.mod: mem_soil_moisture.o
mem_tconv.mod: mem_tconv.o
mem_teb.mod: mem_teb.o
mem_teb_common.mod: mem_teb_common.o
mem_tend.mod: mem_tend.o
mem_turb.mod: mem_turb.o
mem_turb_scalar.mod: mem_turb_scalar.o
mem_varinit.mod: mem_varinit.o
met_driver_coms.mod: met_driver_coms.o
micphys.mod: micphys.o
micro_coms.mod: micro_coms.o
mod_ghostblock.mod: mod_GhostBlock.o
mod_ghostblockpartition.mod: mod_GhostBlockPartition.o
mortality.mod: mortality.o
node_mod.mod: node_mod.o
obs_input.mod: obs_input.o
optimiz_coms.mod: optimiz_coms.o
ozone_const.mod: mod_ozone.o
pft_coms.mod: pft_coms.o
phenology_coms.mod: phenology_coms.o
phenology_startup.mod: phenology_startup.o
physiology_coms.mod: physiology_coms.o
plume_utils.mod: plumerise_vector.o
rad_carma.mod: rad_carma.o
rconstants.mod: rconstants.o
ref_sounding.mod: ref_sounding.o
rk4_coms.mod: rk4_coms.o
rk4_driver.mod: rk4_driver.o
rk4_stepper.mod: rk4_stepper.o
rpara.mod: rpara.o
shcu_vars_const.mod: shcu_vars_const.o
soil_coms.mod: soil_coms.o
teb_spm_start.mod: teb_spm_start.o
teb_vars_const.mod: mem_teb_vars_const.o
therm_lib.mod: therm_lib.o
therm_lib8.mod: therm_lib8.o
turb_coms.mod: turb_coms.o
var_tables.mod: var_tables.o
