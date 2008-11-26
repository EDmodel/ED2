allometry.o : $(ED_DYNAMICS)/allometry.f90 pft_coms.o grid_coms.o soil_coms.o 
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

an_header.o: $(ED_IO)/an_header.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

average_utils.o : $(ED_IO)/average_utils.f90 grid_coms.o misc_coms.o ed_state_vars.o       \
	ed_misc_coms.o max_dims.o canopy_radiation_coms.o consts_coms.o 
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

c34constants.o : $(ED_MEMORY)/c34constants.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

canopy_air_coms.o : $(ED_MEMORY)/canopy_air_coms.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

canopy_radiation_coms.o : $(ED_MEMORY)/canopy_radiation_coms.f90 max_dims.o
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

canopy_update_euler.o : $(ED_DYNAMICS)/canopy_update_euler.f90 ed_state_vars.o grid_coms.o \
	canopy_radiation_coms.o consts_coms.o max_dims.o grid_coms.o therm_lib.o
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

charutils.o: $(ED_UTILS)/charutils.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

consts_coms.o : $(ED_MEMORY)/consts_coms.F90
	cp -f $< $(<F:.F90=.F90)
	$(FPP_COMMAND) $(<F:.F90=.F90)
	rm -f $(<F:.F90=.F90)

dateutils.o: $(ED_UTILS)/dateutils.f90 consts_coms.o
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

decomposition.o : $(ED_DYNAMICS)/decomposition.f90 ed_state_vars.o soil_coms.o grid_coms.o \
	pft_coms.o decomp_coms.o
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

decomp_coms.o : $(ED_MEMORY)/decomp_coms.f90 max_dims.o
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

disturbance.o : $(ED_DYNAMICS)/disturbance.f90 ed_state_vars.o fuse_fiss_utils.o           \
	misc_coms.o disturb_coms.o max_dims.o pft_coms.o consts_coms.o grid_coms.o         \
	decomp_coms.o canopy_air_coms.o
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

disturb_coms.o : $(ED_MEMORY)/disturb_coms.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

edio.o : $(ED_IO)/edio.f90 ed_state_vars.o grid_coms.o ed_misc_coms.o ed_node_coms.o       \
	misc_coms.o canopy_radiation_coms.o consts_coms.o var_tables_array.o 
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

ed_1st.o : $(ED_DRIVER)/ed_1st.f90 ed_para_coms.o misc_coms.o ed_state_vars.o mem_sites.o
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

ed_bare_restart.o : $(ED_INIT)/ed_bare_restart.f90 ed_state_vars.o max_dims.o pft_coms.o   \
	consts_coms.o canopy_air_coms.o
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

ed_driver.o : $(ED_DRIVER)/ed_driver.f90 grid_coms.o ed_state_vars.o misc_coms.o           \
	soil_coms.o ed_node_coms.o consts_coms.o
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

ed_filelist.o : $(ED_UTILS)/ed_filelist.F90 max_dims.o
	cp -f $< $(<F:.F90=.F90)
	$(FPP_COMMAND) $(<F:.F90=.F90)
	rm -f $(<F:.F90=.F90) 

ed_grid.o : $(ED_UTILS)/ed_grid.f90 grid_coms.o ed_node_coms.o max_dims.o consts_coms.o
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90) 

ed_history_io.o : $(ED_IO)/ed_history_io.f90 max_dims.o pft_coms.o misc_coms.o mem_sites.o \
	consts_coms.o ed_misc_coms.o ed_state_vars.o grid_coms.o max_dims.o soil_coms.o    \
	ed_node_coms.o hdf5_coms.o fusion_fission_coms.o 
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(HDF5_INCS) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

ed_init.o : $(ED_INIT)/ed_init.f90 grid_coms.o ed_work_vars.o soil_coms.o ed_node_coms.o   \
	ed_state_vars.o phenology_coms.o phenology_init.o misc_coms.o 
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

ed_init_atm.o : $(ED_INIT)/ed_init_atm.f90 misc_coms.o ed_state_vars.o soil_coms.o         \
	consts_coms.o grid_coms.o fuse_fiss_utils.o ed_node_coms.o pft_coms.o therm_lib.o  \
	canopy_air_coms.o
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

ed_load_namelist.o : $(ED_IO)/ed_load_namelist.f90 ename_coms.o max_dims.o soil_coms.o     \
	met_driver_coms.o mem_sites.o physiology_coms.o phenology_coms.o decomp_coms.o     \
	disturb_coms.o pft_coms.o misc_coms.o grid_coms.o ed_misc_coms.o optimiz_coms.o
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

ed_mem_alloc.o : $(ED_MEMORY)/ed_mem_alloc.f90 max_dims.o ed_mem_grid_dim_defs.o           \
	ed_state_vars.o grid_coms.o mem_sites.o ed_work_vars.o misc_coms.o ed_node_coms.o  \
	soil_coms.o 
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90) 

ed_mem_grid_dim_defs.o : $(ED_MEMORY)/ed_mem_grid_dim_defs.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

ed_met_driver.o : $(ED_DRIVER)/ed_met_driver.f90 max_dims.o met_driver_coms.o              \
	ed_state_vars.o grid_coms.o hdf5_utils.o mem_sites.o misc_coms.o consts_coms.o     \
	therm_lib.o canopy_radiation_coms.o
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(HDF5_INCS) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

ed_misc_coms.o : $(ED_MEMORY)/ed_misc_coms.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

ed_model.o : $(ED_DRIVER)/ed_model.f90 misc_coms.o grid_coms.o ed_state_vars.o             \
	rk4_driver.o ed_node_coms.o disturb_coms.o mem_sites.o consts_coms.o disturbance.o \
	fuse_fiss_utils.o growth_balive.o ed_misc_coms.o consts_coms.o events.o
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

ed_mpass_init.o : $(ED_MPI)/ed_mpass_init.f90 ed_para_coms.o ed_node_coms.o max_dims.o     \
	misc_coms.o ed_misc_coms.o grid_coms.o soil_coms.o met_driver_coms.o mem_sites.o   \
	physiology_coms.o phenology_coms.o decomp_coms.o pft_coms.o disturb_coms.o         \
	optimiz_coms.o met_driver_coms.o ed_state_vars.o ed_work_vars.o
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

ed_node_coms.o : $(ED_MPI)/ed_node_coms.f90 max_dims.o 
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90) 

ed_opspec.o : $(ED_IO)/ed_opspec.f90 max_dims.o grid_coms.o mem_sites.o soil_coms.o        \
	ed_para_coms.o misc_coms.o consts_coms.o physiology_coms.o decomp_coms.o           \
	disturb_coms.o phenology_coms.o pft_coms.o
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

ed_params.o : $(ED_INIT)/ed_params.f90 max_dims.o pft_coms.o disturb_coms.o ed_misc_coms.o \
	met_driver_coms.o canopy_radiation_coms.o decomp_coms.o hydrology_coms.o           \
	misc_coms.o soil_coms.o phenology_coms.o fusion_fission_coms.o 
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

ed_para_coms.o : $(ED_MPI)/ed_para_coms.f90 max_dims.o 
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90) 

ed_para_init.o : $(ED_MPI)/ed_para_init.f90 grid_coms.o ed_node_coms.o ed_para_coms.o      \
	mem_sites.o ed_work_vars.o soil_coms.o ed_state_vars.o
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

ed_state_vars.o : $(ED_MEMORY)/ed_state_vars.f90 grid_coms.o max_dims.o c34constants.o     \
	disturb_coms.o met_driver_coms.o fusion_fission_coms.o phenology_coms.o            \
	misc_coms.o var_tables_array.o ed_node_coms.o soil_coms.o
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

ed_type_init.o : $(ED_INIT)/ed_type_init.f90 ed_state_vars.o max_dims.o grid_coms.o        \
	soil_coms.o therm_lib.o
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

ed_work_vars.o : $(ED_MEMORY)/ed_work_vars.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

ed_xml_config.o : $(ED_IO)/ed_xml_config.f90 pft_coms.o hydrology_coms.o met_driver_coms.o \
	canopy_radiation_coms.o disturb_coms.o phenology_coms.o physiology_coms.o          \
	decomp_coms.o fusion_fission_coms.o ed_misc_coms.o soil_coms.o misc_coms.o
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

ename_coms.o : $(ED_MEMORY)/ename_coms.f90 max_dims.o
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

euler_driver.o : $(ED_DYNAMICS)/euler_driver.f90 ed_state_vars.o misc_coms.o soil_coms.o   \
	consts_coms.o grid_coms.o max_dims.o therm_lib.o
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

events.o : $(ED_DYNAMICS)/events.f90 misc_coms.o ed_state_vars.o grid_coms.o \
	pft_coms.o ed_init_atm.o disturbance.o
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

farq_leuning.o : $(ED_DYNAMICS)/farq_leuning.f90 c34constants.o pft_coms.o                 \
	physiology_coms.o therm_lib.o consts_coms.o
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

fatal_error.o : $(ED_UTILS)/fatal_error.f90 ed_node_coms.o misc_coms.o
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

fire.o : $(ED_DYNAMICS)/fire.f90 ed_state_vars.o pft_coms.o grid_coms.o soil_coms.o        \
	disturb_coms.o
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

forestry.o : $(ED_DYNAMICS)/forestry.f90 ed_state_vars.o disturb_coms.o disturbance.o      \
	fuse_fiss_utils.o max_dims.o grid_coms.o
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

fusion_fission_coms.o : $(ED_MEMORY)/fusion_fission_coms.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

fuse_fiss_utils.o : $(ED_UTILS)/fuse_fiss_utils.f90 ed_state_vars.o pft_coms.o             \
	decomp_coms.o fusion_fission_coms.o disturb_coms.o max_dims.o mem_sites.o          \
	soil_coms.o grid_coms.o consts_coms.o
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

great_circle.o : $(ED_UTILS)/great_circle.f90 consts_coms.o
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90) 

grid_coms.o : $(ED_MEMORY)/grid_coms.f90 max_dims.o
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90) 

growth_balive.o : $(ED_DYNAMICS)/growth_balive.f90 ed_state_vars.o pft_coms.o              \
	decomp_coms.o consts_coms.o physiology_coms.o grid_coms.o
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

h5_output.o : $(ED_IO)/h5_output.F90 an_header.o var_tables_array.o mem_sites.o            \
	misc_coms.o ed_misc_coms.o grid_coms.o hdf5_coms.o ed_node_coms.o max_dims.o       \
	ed_state_vars.o fusion_fission_coms.o
	cp -f $< $(<F:.F90=.F90)
	$(FPP_COMMAND) $(HDF5_INCS) $(<F:.F90=.F90)
	rm -f $(<F:.F90=.F90)

hdf5_coms.o : $(ED_MEMORY)/hdf5_coms.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(HDF5_INCS) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

hdf5_utils.o : $(ED_UTILS)/hdf5_utils.f90 hdf5_coms.o
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(HDF5_INCS) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

hydrology_coms.o: $(ED_MEMORY)/hydrology_coms.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

hydrology_constants.o: $(ED_MEMORY)/hydrology_constants.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

init_hydro_sites.o : $(ED_INIT)/init_hydro_sites.f90 soil_coms.o grid_coms.o misc_coms.o   \
	mem_sites.o ed_state_vars.o max_dims.o
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90) 
	rm -f $(<F:.f90=.f90)

invmondays.o : $(ED_UTILS)/invmondays.f90 misc_coms.o
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

landuse_init.o : $(ED_INIT)/landuse_init.f90 ed_state_vars.o consts_coms.o disturb_coms.o  \
	misc_coms.o grid_coms.o
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

leaf_database.o : $(ED_IO)/leaf_database.f90 hdf5_utils.o soil_coms.o grid_coms.o
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

libxml2f90.f90_pp.o : $(ED_UTILS)/libxml2f90.f90_pp.f90 
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

lsm_hyd.o : $(ED_DYNAMICS)/lsm_hyd.f90 hydrology_constants.o grid_coms.o hydrology_coms.o  \
	ed_state_vars.o ed_node_coms.o soil_coms.o misc_coms.o therm_lib.o consts_coms.o
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

max_dims.o : $(ED_MEMORY)/max_dims.F90
	cp -f $< $(<F:.F90=.F90)
	$(FPP_COMMAND) $(<F:.F90=.F90)
	rm -f $(<F:.F90=.F90)

mem_sites.o : $(ED_MEMORY)/mem_sites.f90 max_dims.o
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

met_driver_coms.o : $(ED_MEMORY)/met_driver_coms.f90 max_dims.o
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

misc_coms.o : $(ED_MEMORY)/misc_coms.f90 max_dims.o 
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

mortality.o : $(ED_DYNAMICS)/mortality.f90 ed_state_vars.o pft_coms.o disturb_coms.o       \
	consts_coms.o
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

numutils.o: $(ED_UTILS)/numutils.f90 therm_lib.o consts_coms.o
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

optimiz_coms.o : $(ED_MEMORY)/optimiz_coms.f90 max_dims.o
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90) 

phenology_aux.o : $(ED_DYNAMICS)/phenology_aux.f90 phenology_coms.o max_dims.o pft_coms.o  \
	ed_state_vars.o consts_coms.o
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

phenology_coms.o : $(ED_MEMORY)/phenology_coms.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

phenology_driv.o : $(ED_DYNAMICS)/phenology_driv.f90 \
	soil_coms.o pft_coms.o decomp_coms.o phenology_coms.o \
	ed_state_vars.o grid_coms.o canopy_air_coms.o
	cp -f $< $(<F:.f90=.f90) 
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

phenology_init.o : $(ED_INIT)/phenology_init.f90 grid_coms.o \
	phenology_coms.o misc_coms.o ed_state_vars.o
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

photosyn_driv.o : $(ED_DYNAMICS)/photosyn_driv.f90 ed_state_vars.o phenology_coms.o        \
	misc_coms.o grid_coms.o pft_coms.o decomp_coms.o soil_coms.o 
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

physiology_coms.o : $(ED_MEMORY)/physiology_coms.f90 max_dims.o
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

pft_coms.o : $(ED_MEMORY)/pft_coms.f90 max_dims.o
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

radiate_driver.o : $(ED_DYNAMICS)/radiate_driver.f90 misc_coms.o ed_state_vars.o           \
	canopy_radiation_coms.o consts_coms.o grid_coms.o soil_coms.o max_dims.o
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

reproduction.o : $(ED_DYNAMICS)/reproduction.f90 ed_state_vars.o misc_coms.o pft_coms.o    \
	decomp_coms.o fusion_fission_coms.o max_dims.o fuse_fiss_utils.o \
	phenology_coms.o consts_coms.o canopy_air_coms.o
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

rk4_derivs.o : $(ED_DYNAMICS)/rk4_derivs.F90 ed_state_vars.o consts_coms.o grid_coms.o     \
	max_dims.o consts_coms.o grid_coms.o soil_coms.o misc_coms.o                       \
	canopy_radiation_coms.o therm_lib.o pft_coms.o ed_misc_coms.o canopy_air_coms.o
	cp -f $< $(<F:.F90=.F90)
	$(FPP_COMMAND) $(<F:.F90=.F90)
	rm -f $(<F:.F90=.F90)

rk4_driver.o : $(ED_DYNAMICS)/rk4_driver.f90 ed_state_vars.o grid_coms.o max_dims.o        \
	misc_coms.o consts_coms.o soil_coms.o canopy_radiation_coms.o ed_misc_coms.o       \
	canopy_air_coms.o therm_lib.o
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

rk4_integ_utils.o : $(ED_DYNAMICS)/rk4_integ_utils.f90 ed_state_vars.o rk4_stepper.o       \
	hydrology_coms.o grid_coms.o soil_coms.o consts_coms.o grid_coms.o ed_misc_coms.o  \
	canopy_radiation_coms.o therm_lib.o max_dims.o misc_coms.o
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

rk4_stepper.o : $(ED_DYNAMICS)/rk4_stepper.F90 ed_state_vars.o grid_coms.o soil_coms.o     \
	canopy_radiation_coms.o consts_coms.o canopy_air_coms.o therm_lib.o
	cp -f $< $(<F:.F90=.F90)
	$(FPP_COMMAND) $(<F:.F90=.F90)
	rm -f $(<F:.F90=.F90)

rsys.o: $(ED_UTILS)/rsys.F90
	cp -f $< $(<F:.F90=.F90)
	$(FPP_COMMAND) $(<F:.F90=.F90)
	rm -f $(<F:.F90=.F90)

soil_coms.o : $(ED_MEMORY)/soil_coms.F90 max_dims.o grid_coms.o
	cp -f $< $(<F:.F90=.F90)
	$(FPP_COMMAND) $(<F:.F90=.F90)
	rm -f $(<F:.F90=.F90)

structural_growth.o : $(ED_DYNAMICS)/structural_growth.f90 ed_state_vars.o pft_coms.o      \
	decomp_coms.o max_dims.o consts_coms.o
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

therm_lib.o: $(ED_UTILS)/therm_lib.f90 consts_coms.o
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

twostream_rad.o : $(ED_DYNAMICS)/twostream_rad.f90 pft_coms.o canopy_radiation_coms.o      \
	consts_coms.o
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

utils_c.o: $(ED_UTILS)/utils_c.c
	cp -f $< $(<F:.c=.c)
	$(CXX_COMMAND) $<
	rm -f $(<F:.c=.c)

utils_f.o: $(ED_UTILS)/utils_f.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

var_tables_array.o : $(ED_MEMORY)/var_tables_array.f90 
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)
