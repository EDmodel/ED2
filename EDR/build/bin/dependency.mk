# DO NOT DELETE THIS LINE - used by make depend
ed_1st.o: ed_mem_alloc.mod ed_misc_coms.mod ed_para_coms.mod ed_state_vars.mod
ed_driver.o: consts_coms.mod detailed_coms.mod ed_misc_coms.mod ed_node_coms.mod
ed_driver.o: ed_state_vars.mod fuse_fiss_utils.mod grid_coms.mod
ed_driver.o: phenology_aux.mod soil_coms.mod
ed_met_driver.o: canopy_air_coms.mod canopy_radiation_coms.mod consts_coms.mod
ed_met_driver.o: ed_max_dims.mod ed_misc_coms.mod ed_node_coms.mod
ed_met_driver.o: ed_state_vars.mod grid_coms.mod hdf5_utils.mod mem_polygons.mod
ed_met_driver.o: met_driver_coms.mod pft_coms.mod therm_lib.mod
ed_model.o: average_utils.mod consts_coms.mod ed_misc_coms.mod ed_node_coms.mod
ed_model.o: ed_state_vars.mod grid_coms.mod mem_polygons.mod rk4_coms.mod
ed_model.o: rk4_driver.mod
bdf2_solver.o: consts_coms.mod ed_misc_coms.mod ed_state_vars.mod
bdf2_solver.o: ed_therm_lib.mod grid_coms.mod rk4_coms.mod soil_coms.mod
bdf2_solver.o: therm_lib8.mod
canopy_struct_dynamics.o: allometry.mod canopy_air_coms.mod
canopy_struct_dynamics.o: canopy_layer_coms.mod consts_coms.mod ed_misc_coms.mod
canopy_struct_dynamics.o: ed_state_vars.mod grid_coms.mod met_driver_coms.mod
canopy_struct_dynamics.o: pft_coms.mod phenology_coms.mod physiology_coms.mod
canopy_struct_dynamics.o: rk4_coms.mod soil_coms.mod therm_lib.mod
disturbance.o: allometry.mod budget_utils.mod consts_coms.mod decomp_coms.mod
disturbance.o: disturb_coms.mod ed_max_dims.mod ed_misc_coms.mod
disturbance.o: ed_state_vars.mod ed_therm_lib.mod fuse_fiss_utils.mod
disturbance.o: grid_coms.mod mem_polygons.mod mortality.mod pft_coms.mod
disturbance.o: phenology_aux.mod phenology_coms.mod therm_lib.mod
euler_driver.o: budget_utils.mod canopy_air_coms.mod consts_coms.mod
euler_driver.o: ed_max_dims.mod ed_misc_coms.mod ed_state_vars.mod grid_coms.mod
euler_driver.o: hydrology_coms.mod met_driver_coms.mod rk4_coms.mod
euler_driver.o: rk4_driver.mod rk4_stepper.mod soil_coms.mod therm_lib.mod
euler_driver.o: therm_lib8.mod
events.o: allometry.mod budget_utils.mod consts_coms.mod decomp_coms.mod
events.o: disturbance_utils.mod ed_misc_coms.mod ed_state_vars.mod
events.o: ed_therm_lib.mod fuse_fiss_utils.mod grid_coms.mod pft_coms.mod
events.o: therm_lib.mod
farq_leuning.o: c34constants.mod consts_coms.mod pft_coms.mod phenology_coms.mod
farq_leuning.o: physiology_coms.mod rk4_coms.mod therm_lib8.mod
fire.o: consts_coms.mod disturb_coms.mod ed_misc_coms.mod ed_state_vars.mod
fire.o: grid_coms.mod soil_coms.mod
forestry.o: disturb_coms.mod ed_max_dims.mod ed_state_vars.mod
forestry.o: fuse_fiss_utils.mod grid_coms.mod
growth_balive.o: allometry.mod budget_utils.mod consts_coms.mod decomp_coms.mod
growth_balive.o: ed_max_dims.mod ed_misc_coms.mod ed_state_vars.mod
growth_balive.o: ed_therm_lib.mod fuse_fiss_utils.mod grid_coms.mod
growth_balive.o: mortality.mod pft_coms.mod phenology_coms.mod
growth_balive.o: physiology_coms.mod
heun_driver.o: budget_utils.mod canopy_air_coms.mod consts_coms.mod
heun_driver.o: ed_max_dims.mod ed_misc_coms.mod ed_state_vars.mod grid_coms.mod
heun_driver.o: hydrology_coms.mod met_driver_coms.mod rk4_coms.mod
heun_driver.o: rk4_driver.mod rk4_stepper.mod soil_coms.mod therm_lib.mod
heun_driver.o: therm_lib8.mod
hybrid_driver.o: budget_utils.mod consts_coms.mod ed_max_dims.mod
hybrid_driver.o: ed_misc_coms.mod ed_state_vars.mod grid_coms.mod
hybrid_driver.o: hydrology_coms.mod met_driver_coms.mod rk4_coms.mod
hybrid_driver.o: rk4_driver.mod rk4_stepper.mod soil_coms.mod therm_lib.mod
hybrid_driver.o: therm_lib8.mod
lsm_hyd.o: consts_coms.mod ed_misc_coms.mod ed_node_coms.mod ed_state_vars.mod
lsm_hyd.o: grid_coms.mod hydrology_coms.mod hydrology_constants.mod pft_coms.mod
lsm_hyd.o: soil_coms.mod therm_lib.mod
mortality.o: consts_coms.mod disturb_coms.mod ed_max_dims.mod ed_misc_coms.mod
mortality.o: ed_state_vars.mod pft_coms.mod
multiple_scatter.o: canopy_radiation_coms.mod consts_coms.mod ed_max_dims.mod
multiple_scatter.o: rk4_coms.mod
old_twostream_rad.o: canopy_radiation_coms.mod consts_coms.mod ed_max_dims.mod
old_twostream_rad.o: rk4_coms.mod
phenology_aux.o: allometry.mod consts_coms.mod ed_max_dims.mod ed_state_vars.mod
phenology_aux.o: ed_therm_lib.mod grid_coms.mod pft_coms.mod phenology_coms.mod
phenology_aux.o: soil_coms.mod therm_lib.mod
phenology_driv.o: allometry.mod consts_coms.mod decomp_coms.mod ed_max_dims.mod
phenology_driv.o: ed_misc_coms.mod ed_state_vars.mod ed_therm_lib.mod
phenology_driv.o: grid_coms.mod pft_coms.mod phenology_aux.mod
phenology_driv.o: phenology_coms.mod soil_coms.mod
photosyn_driv.o: allometry.mod consts_coms.mod ed_max_dims.mod ed_misc_coms.mod
photosyn_driv.o: ed_state_vars.mod farq_leuning.mod met_driver_coms.mod
photosyn_driv.o: pft_coms.mod phenology_coms.mod physiology_coms.mod
photosyn_driv.o: soil_coms.mod therm_lib.mod
radiate_driver.o: allometry.mod canopy_layer_coms.mod canopy_radiation_coms.mod
radiate_driver.o: consts_coms.mod ed_max_dims.mod ed_misc_coms.mod
radiate_driver.o: ed_state_vars.mod grid_coms.mod soil_coms.mod
reproduction.o: allometry.mod budget_utils.mod consts_coms.mod decomp_coms.mod
reproduction.o: ed_max_dims.mod ed_misc_coms.mod ed_state_vars.mod
reproduction.o: ed_therm_lib.mod fuse_fiss_utils.mod grid_coms.mod
reproduction.o: mem_polygons.mod pft_coms.mod phenology_aux.mod
reproduction.o: phenology_coms.mod therm_lib.mod
rk4_derivs.o: budget_utils.mod canopy_struct_dynamics.mod consts_coms.mod
rk4_derivs.o: ed_max_dims.mod ed_misc_coms.mod ed_state_vars.mod grid_coms.mod
rk4_derivs.o: pft_coms.mod physiology_coms.mod rk4_coms.mod soil_coms.mod
rk4_derivs.o: therm_lib8.mod
rk4_driver.o: allometry.mod budget_utils.mod canopy_air_coms.mod consts_coms.mod
rk4_driver.o: disturb_coms.mod ed_misc_coms.mod ed_state_vars.mod grid_coms.mod
rk4_driver.o: met_driver_coms.mod phenology_coms.mod rk4_coms.mod soil_coms.mod
rk4_driver.o: therm_lib.mod
rk4_integ_utils.o: c34constants.mod canopy_air_coms.mod canopy_layer_coms.mod
rk4_integ_utils.o: canopy_radiation_coms.mod consts_coms.mod ed_max_dims.mod
rk4_integ_utils.o: ed_misc_coms.mod ed_state_vars.mod grid_coms.mod
rk4_integ_utils.o: hydrology_coms.mod rk4_coms.mod rk4_stepper.mod soil_coms.mod
rk4_integ_utils.o: therm_lib8.mod
rk4_misc.o: canopy_air_coms.mod canopy_struct_dynamics.mod consts_coms.mod
rk4_misc.o: ed_max_dims.mod ed_misc_coms.mod ed_state_vars.mod ed_therm_lib.mod
rk4_misc.o: grid_coms.mod rk4_coms.mod soil_coms.mod therm_lib.mod
rk4_misc.o: therm_lib8.mod
rk4_stepper.o: ed_state_vars.mod grid_coms.mod rk4_coms.mod soil_coms.mod
soil_respiration.o: consts_coms.mod decomp_coms.mod ed_misc_coms.mod
soil_respiration.o: ed_state_vars.mod farq_leuning.mod pft_coms.mod
soil_respiration.o: physiology_coms.mod rk4_coms.mod soil_coms.mod therm_lib.mod
structural_growth.o: allometry.mod consts_coms.mod decomp_coms.mod
structural_growth.o: detailed_coms.mod ed_max_dims.mod ed_misc_coms.mod
structural_growth.o: ed_state_vars.mod ed_therm_lib.mod pft_coms.mod
structural_growth.o: physiology_coms.mod
twostream_rad.o: canopy_radiation_coms.mod consts_coms.mod ed_max_dims.mod
twostream_rad.o: rk4_coms.mod
vegetation_dynamics.o: average_utils.mod consts_coms.mod disturbance_utils.mod
vegetation_dynamics.o: ed_misc_coms.mod ed_state_vars.mod fuse_fiss_utils.mod
vegetation_dynamics.o: grid_coms.mod growth_balive.mod mem_polygons.mod
ed_bigleaf_init.o: allometry.mod consts_coms.mod ed_max_dims.mod
ed_bigleaf_init.o: ed_misc_coms.mod ed_node_coms.mod ed_state_vars.mod
ed_bigleaf_init.o: fuse_fiss_utils.mod pft_coms.mod
ed_init.o: consts_coms.mod ed_max_dims.mod ed_misc_coms.mod ed_node_coms.mod
ed_init.o: ed_state_vars.mod ed_work_vars.mod grid_coms.mod mem_polygons.mod
ed_init.o: phenology_coms.mod phenology_startup.mod rk4_coms.mod soil_coms.mod
ed_init_atm.o: budget_utils.mod canopy_layer_coms.mod canopy_struct_dynamics.mod
ed_init_atm.o: consts_coms.mod ed_misc_coms.mod ed_node_coms.mod
ed_init_atm.o: ed_state_vars.mod ed_therm_lib.mod fuse_fiss_utils.mod
ed_init_atm.o: grid_coms.mod met_driver_coms.mod pft_coms.mod soil_coms.mod
ed_init_atm.o: therm_lib.mod
ed_nbg_init.o: allometry.mod consts_coms.mod ed_max_dims.mod ed_misc_coms.mod
ed_nbg_init.o: ed_state_vars.mod fuse_fiss_utils.mod grid_coms.mod pft_coms.mod
ed_nbg_init.o: physiology_coms.mod
ed_params.o: allometry.mod canopy_air_coms.mod canopy_layer_coms.mod
ed_params.o: canopy_radiation_coms.mod consts_coms.mod decomp_coms.mod
ed_params.o: detailed_coms.mod disturb_coms.mod ed_max_dims.mod ed_misc_coms.mod
ed_params.o: ed_therm_lib.mod fusion_fission_coms.mod grid_coms.mod
ed_params.o: hydrology_coms.mod met_driver_coms.mod pft_coms.mod
ed_params.o: phenology_coms.mod physiology_coms.mod rk4_coms.mod soil_coms.mod
ed_type_init.o: allometry.mod canopy_air_coms.mod consts_coms.mod
ed_type_init.o: ed_max_dims.mod ed_misc_coms.mod ed_state_vars.mod
ed_type_init.o: ed_therm_lib.mod grid_coms.mod pft_coms.mod phenology_coms.mod
ed_type_init.o: soil_coms.mod therm_lib.mod
init_hydro_sites.o: ed_max_dims.mod ed_misc_coms.mod ed_state_vars.mod
init_hydro_sites.o: grid_coms.mod mem_polygons.mod soil_coms.mod
landuse_init.o: consts_coms.mod disturb_coms.mod ed_max_dims.mod
landuse_init.o: ed_misc_coms.mod ed_state_vars.mod grid_coms.mod
phenology_startup.o: ed_max_dims.mod ed_misc_coms.mod ed_state_vars.mod
phenology_startup.o: grid_coms.mod phenology_aux.mod phenology_coms.mod
average_utils.o: consts_coms.mod ed_max_dims.mod ed_misc_coms.mod
average_utils.o: ed_state_vars.mod grid_coms.mod met_driver_coms.mod
average_utils.o: physiology_coms.mod soil_coms.mod therm_lib.mod
ed_init_full_history.o: allometry.mod ed_max_dims.mod ed_misc_coms.mod
ed_init_full_history.o: ed_node_coms.mod ed_state_vars.mod
ed_init_full_history.o: fusion_fission_coms.mod grid_coms.mod 
ed_init_full_history.o: hdf5_coms.mod phenology_startup.mod soil_coms.mod
ed_load_namelist.o: canopy_air_coms.mod canopy_layer_coms.mod
ed_load_namelist.o: canopy_radiation_coms.mod consts_coms.mod decomp_coms.mod
ed_load_namelist.o: detailed_coms.mod disturb_coms.mod ed_max_dims.mod
ed_load_namelist.o: ed_misc_coms.mod ed_para_coms.mod ename_coms.mod
ed_load_namelist.o: grid_coms.mod mem_polygons.mod met_driver_coms.mod
ed_load_namelist.o: optimiz_coms.mod pft_coms.mod phenology_coms.mod
ed_load_namelist.o: physiology_coms.mod rk4_coms.mod soil_coms.mod
ed_opspec.o: canopy_air_coms.mod canopy_layer_coms.mod canopy_radiation_coms.mod
ed_opspec.o: consts_coms.mod decomp_coms.mod detailed_coms.mod disturb_coms.mod
ed_opspec.o: ed_max_dims.mod ed_misc_coms.mod ed_para_coms.mod grid_coms.mod
ed_opspec.o: mem_polygons.mod met_driver_coms.mod pft_coms.mod
ed_opspec.o: phenology_coms.mod physiology_coms.mod rk4_coms.mod soil_coms.mod
ed_print.o: ed_max_dims.mod ed_misc_coms.mod ed_node_coms.mod ed_state_vars.mod
ed_print.o: ed_var_tables.mod
ed_read_ed10_20_history.o: allometry.mod consts_coms.mod disturb_coms.mod
ed_read_ed10_20_history.o: ed_max_dims.mod ed_misc_coms.mod ed_state_vars.mod
ed_read_ed10_20_history.o: fuse_fiss_utils.mod grid_coms.mod mem_polygons.mod
ed_read_ed10_20_history.o: pft_coms.mod
ed_read_ed21_history.o: allometry.mod consts_coms.mod disturb_coms.mod
ed_read_ed21_history.o: ed_max_dims.mod ed_misc_coms.mod ed_state_vars.mod
ed_read_ed21_history.o: fuse_fiss_utils.mod grid_coms.mod  hdf5_coms.mod
ed_read_ed21_history.o: pft_coms.mod soil_coms.mod
ed_xml_config.o: canopy_radiation_coms.mod decomp_coms.mod disturb_coms.mod
ed_xml_config.o: ed_max_dims.mod ed_misc_coms.mod fusion_fission_coms.mod
ed_xml_config.o: grid_coms.mod hydrology_coms.mod met_driver_coms.mod
ed_xml_config.o: pft_coms.mod phenology_coms.mod physiology_coms.mod
ed_xml_config.o: rk4_coms.mod soil_coms.mod
edio.o: average_utils.mod ed_misc_coms.mod ed_node_coms.mod ed_state_vars.mod
edio.o: grid_coms.mod
h5_output.o: an_header.mod ed_max_dims.mod ed_misc_coms.mod ed_node_coms.mod
h5_output.o: ed_state_vars.mod ed_var_tables.mod fusion_fission_coms.mod
h5_output.o: grid_coms.mod  hdf5_coms.mod
leaf_database.o: grid_coms.mod hdf5_utils.mod soil_coms.mod
canopy_air_coms.o: consts_coms.mod therm_lib.mod therm_lib8.mod
canopy_radiation_coms.o: ed_max_dims.mod
consts_coms.o: 
decomp_coms.o: ed_max_dims.mod
disturb_coms.o: ed_max_dims.mod
ed_max_dims.o: 
ed_mem_alloc.o: ed_max_dims.mod ed_mem_grid_dim_defs.mod ed_node_coms.mod
ed_mem_alloc.o: ed_state_vars.mod ed_work_vars.mod grid_coms.mod
ed_mem_alloc.o: mem_polygons.mod
ed_misc_coms.o: ed_max_dims.mod
ed_state_vars.o: disturb_coms.mod ed_max_dims.mod ed_misc_coms.mod
ed_state_vars.o: ed_node_coms.mod ed_var_tables.mod fusion_fission_coms.mod
ed_state_vars.o: grid_coms.mod met_driver_coms.mod phenology_coms.mod
ed_state_vars.o: soil_coms.mod
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
soil_coms.o: consts_coms.mod ed_max_dims.mod grid_coms.mod 
ed_mpass_init.o: canopy_air_coms.mod canopy_layer_coms.mod
ed_mpass_init.o: canopy_radiation_coms.mod decomp_coms.mod detailed_coms.mod
ed_mpass_init.o: disturb_coms.mod ed_max_dims.mod ed_mem_alloc.mod
ed_mpass_init.o: ed_misc_coms.mod ed_node_coms.mod ed_para_coms.mod
ed_mpass_init.o: ed_state_vars.mod ed_work_vars.mod grid_coms.mod
ed_mpass_init.o: mem_polygons.mod met_driver_coms.mod optimiz_coms.mod
ed_mpass_init.o: pft_coms.mod phenology_coms.mod physiology_coms.mod
ed_mpass_init.o: rk4_coms.mod soil_coms.mod
ed_node_coms.o: ed_max_dims.mod
ed_para_coms.o: ed_max_dims.mod
ed_para_init.o: ed_max_dims.mod ed_misc_coms.mod ed_node_coms.mod
ed_para_init.o: ed_para_coms.mod ed_work_vars.mod grid_coms.mod 
ed_para_init.o: hdf5_coms.mod mem_polygons.mod soil_coms.mod
allometry.o: consts_coms.mod ed_misc_coms.mod grid_coms.mod pft_coms.mod
allometry.o: rk4_coms.mod soil_coms.mod
budget_utils.o: consts_coms.mod ed_max_dims.mod ed_misc_coms.mod
budget_utils.o: ed_state_vars.mod grid_coms.mod rk4_coms.mod soil_coms.mod
budget_utils.o: therm_lib.mod
dateutils.o: consts_coms.mod
ed_filelist.o: ed_max_dims.mod
ed_grid.o: consts_coms.mod ed_max_dims.mod ed_node_coms.mod grid_coms.mod
ed_therm_lib.o: canopy_air_coms.mod consts_coms.mod ed_max_dims.mod
ed_therm_lib.o: ed_misc_coms.mod ed_state_vars.mod grid_coms.mod pft_coms.mod
ed_therm_lib.o: rk4_coms.mod soil_coms.mod therm_lib.mod therm_lib8.mod
fatal_error.o: ed_node_coms.mod
fuse_fiss_utils.o: allometry.mod budget_utils.mod canopy_layer_coms.mod
fuse_fiss_utils.o: consts_coms.mod decomp_coms.mod disturb_coms.mod
fuse_fiss_utils.o: ed_max_dims.mod ed_misc_coms.mod ed_node_coms.mod
fuse_fiss_utils.o: ed_state_vars.mod fusion_fission_coms.mod grid_coms.mod
fuse_fiss_utils.o: mem_polygons.mod pft_coms.mod rk4_coms.mod soil_coms.mod
fuse_fiss_utils.o: therm_lib.mod
great_circle.o: consts_coms.mod
hdf5_utils.o: hdf5_coms.mod
invmondays.o: ed_misc_coms.mod
lapse.o: consts_coms.mod ed_misc_coms.mod ed_state_vars.mod met_driver_coms.mod
numutils.o: consts_coms.mod therm_lib.mod
radiate_utils.o: canopy_radiation_coms.mod consts_coms.mod ed_max_dims.mod
radiate_utils.o: ed_misc_coms.mod ed_state_vars.mod met_driver_coms.mod
stable_cohorts.o: ed_max_dims.mod ed_state_vars.mod pft_coms.mod
stable_cohorts.o: phenology_coms.mod
therm_lib.o: consts_coms.mod
therm_lib8.o: consts_coms.mod therm_lib.mod
update_derived_props.o: allometry.mod canopy_air_coms.mod consts_coms.mod
update_derived_props.o: decomp_coms.mod ed_max_dims.mod ed_misc_coms.mod
update_derived_props.o: ed_state_vars.mod ed_therm_lib.mod fuse_fiss_utils.mod
update_derived_props.o: grid_coms.mod pft_coms.mod soil_coms.mod therm_lib.mod
utils_c.o: ../../src/include/utils_sub_names.h
utils_c.o:
allometry.mod: allometry.o
an_header.mod: an_header.o
average_utils.mod: average_utils.o
budget_utils.mod: budget_utils.o
c34constants.mod: c34constants.o
canopy_air_coms.mod: canopy_air_coms.o
canopy_layer_coms.mod: canopy_layer_coms.o
canopy_radiation_coms.mod: canopy_radiation_coms.o
canopy_struct_dynamics.mod: canopy_struct_dynamics.o
consts_coms.mod: consts_coms.o
decomp_coms.mod: decomp_coms.o
detailed_coms.mod: detailed_coms.o
disturb_coms.mod: disturb_coms.o
disturbance_utils.mod: disturbance.o
ed_max_dims.mod: ed_max_dims.o
ed_mem_alloc.mod: ed_mem_alloc.o
ed_mem_grid_dim_defs.mod: ed_mem_grid_dim_defs.o
ed_misc_coms.mod: ed_misc_coms.o
ed_node_coms.mod: ed_node_coms.o
ed_para_coms.mod: ed_para_coms.o
ed_state_vars.mod: ed_state_vars.o
ed_therm_lib.mod: ed_therm_lib.o
ed_var_tables.mod: ed_var_tables.o
ed_work_vars.mod: ed_work_vars.o
ename_coms.mod: ename_coms.o
farq_leuning.mod: farq_leuning.o
fuse_fiss_utils.mod: fuse_fiss_utils.o
fusion_fission_coms.mod: fusion_fission_coms.o
grid_coms.mod: grid_coms.o
growth_balive.mod: growth_balive.o
hdf5_coms.mod: hdf5_coms.o
hdf5_utils.mod: hdf5_utils.o
hydrology_coms.mod: hydrology_coms.o
hydrology_constants.mod: hydrology_constants.o
libxml2f90_interface_module.mod: libxml2f90.f90_pp.o
libxml2f90_module.mod: libxml2f90.f90_pp.o
libxml2f90_strings_module.mod: libxml2f90.f90_pp.o
ll_module.mod: libxml2f90.f90_pp.o
mem_polygons.mod: mem_polygons.o
met_driver_coms.mod: met_driver_coms.o
mortality.mod: mortality.o
optimiz_coms.mod: optimiz_coms.o
pft_coms.mod: pft_coms.o
phenology_aux.mod: phenology_aux.o
phenology_coms.mod: phenology_coms.o
phenology_startup.mod: phenology_startup.o
physiology_coms.mod: physiology_coms.o
rk4_coms.mod: rk4_coms.o
rk4_driver.mod: rk4_driver.o
rk4_stepper.mod: rk4_stepper.o
soil_coms.mod: soil_coms.o
therm_lib.mod: therm_lib.o
therm_lib8.mod: therm_lib8.o
