#Makefile objects_edofl.mk

# Define main source.

MAIN = $(ED_DRIVER)/edmain.f90
MAINOBJ = edmain.o


# Define objects.

OBJ_MODEL =                        \
	allometry.o                \
	an_header.o                \
	average_utils.o            \
	budget_utils.o             \
	canopy_air_coms.o          \
	canopy_layer_coms.o        \
	canopy_radiation_coms.o    \
	canopy_struct_dynamics.o   \
	c34constants.o             \
	charutils.o                \
	consts_coms.o              \
	dateutils.o                \
	decomp_coms.o              \
	disturbance.o              \
	disturb_coms.o             \
	edio.o                     \
	ed_1st.o                   \
	ed_driver.o                \
	ed_filelist.o              \
	ed_grid.o                  \
	ed_init.o                  \
	ed_init_atm.o              \
	ed_init_full_history.o     \
	ed_load_namelist.o         \
	ed_max_dims.o              \
	ed_mem_alloc.o             \
	ed_mem_grid_dim_defs.o     \
	ed_met_driver.o            \
	ed_misc_coms.o             \
	ed_model.o                 \
	ed_mpass_init.o            \
	ed_nbg_init.o              \
	ed_node_coms.o             \
	ed_opspec.o                \
	ed_params.o                \
	ed_para_coms.o             \
	ed_para_init.o             \
	ed_print.o                 \
	ed_read_ed10_20_history.o  \
	ed_read_ed21_history.o     \
	ed_state_vars.o            \
	ed_therm_lib.o             \
	ed_type_init.o             \
	ed_var_tables.o            \
	ed_work_vars.o             \
	ed_xml_config.o            \
	ename_coms.o               \
	euler_driver.o             \
	events.o                   \
	farq_leuning.o             \
	fatal_error.o              \
	fire.o                     \
	forestry.o                 \
	fusion_fission_coms.o      \
	fuse_fiss_utils.o          \
	great_circle.o             \
	grid_coms.o                \
	growth_balive.o            \
	h5_output.o                \
	hdf5_coms.o                \
	hdf5_utils.o               \
	heun_driver.o              \
	hydrology_coms.o           \
	hydrology_constants.o      \
	init_hydro_sites.o         \
	invmondays.o               \
	landuse_init.o             \
	lapse.o                    \
	leaf_database.o            \
	libxml2f90.f90_pp.o        \
	lsm_hyd.o                  \
	mem_polygons.o             \
	met_driver_coms.o          \
	mortality.o                \
	numutils.o                 \
	optimiz_coms.o             \
	phenology_aux.o            \
	phenology_coms.o           \
	phenology_driv.o           \
	phenology_startup.o        \
	photosyn_driv.o            \
	physiology_coms.o          \
	pft_coms.o                 \
	radiate_driver.o           \
	radiate_utils.o            \
	reproduction.o             \
	rk4_coms.o                 \
	rk4_derivs.o               \
	rk4_driver.o               \
	rk4_integ_utils.o          \
	rk4_misc.o                 \
	rk4_stepper.o              \
	rsys.o                     \
	soil_coms.o                \
	soil_respiration.o         \
	stable_cohorts.o           \
	structural_growth.o        \
	therm_lib.o                \
	therm_lib8.o               \
	twostream_rad.o            \
	update_derived_props.o     \
	utils_c.o                  \
	utils_f.o                  \
	vegetation_dynamics.o
