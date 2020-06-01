allometry.o : $(ED_UTILS)/allometry.f90
	/bin/rm -f $(<F:.f90=.f90)
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)

average_utils.o : $(ED_IO)/average_utils.f90
	/bin/rm -f $(<F:.f90=.f90)
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_LOWO_COMMAND) $(<F:.f90=.f90)

bdf2_solver.o : $(ED_DYNAMICS)/bdf2_solver.f90
	/bin/rm -f $(<F:.f90=.f90)
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)

budget_utils.o : $(ED_UTILS)/budget_utils.f90
	/bin/rm -f $(<F:.f90=.f90)
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)

canopy_air_coms.o : $(ED_MEMORY)/canopy_air_coms.f90
	/bin/rm -f $(<F:.f90=.f90)
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)

canopy_layer_coms.o : $(ED_MEMORY)/canopy_layer_coms.f90
	/bin/rm -f $(<F:.f90=.f90)
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)

canopy_radiation_coms.o : $(ED_MEMORY)/canopy_radiation_coms.f90
	/bin/rm -f $(<F:.f90=.f90)
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)

canopy_struct_dynamics.o : $(ED_DYNAMICS)/canopy_struct_dynamics.f90
	/bin/rm -f $(<F:.f90=.f90)
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)

c34constants.o : $(ED_MEMORY)/c34constants.f90
	/bin/rm -f $(<F:.f90=.f90)
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)

charutils.o: $(ED_UTILS)/charutils.f90
	/bin/rm -f $(<F:.f90=.f90)
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)

consts_coms.o : $(ED_MEMORY)/consts_coms.F90
	/bin/rm -f $(<F:.F90=.F90)
	/bin/cp -f $< $(<F:.F90=.F90)
	$(FPP_COMMAND) $(<F:.F90=.F90)

dateutils.o: $(ED_UTILS)/dateutils.f90
	/bin/rm -f $(<F:.f90=.f90)
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)

decomp_coms.o : $(ED_MEMORY)/decomp_coms.f90
	/bin/rm -f $(<F:.f90=.f90)
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)

detailed_coms.o : $(ED_MEMORY)/detailed_coms.f90
	/bin/rm -f $(<F:.f90=.f90)
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)

disturbance.o : $(ED_DYNAMICS)/disturbance.f90
	/bin/rm -f $(<F:.f90=.f90)
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)

disturb_coms.o : $(ED_MEMORY)/disturb_coms.f90
	/bin/rm -f $(<F:.f90=.f90)
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)

edio.o : $(ED_IO)/edio.f90
	/bin/rm -f $(<F:.f90=.f90)
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)

ed_1st.o : $(ED_DRIVER)/ed_1st.F90
	/bin/rm -f $(<F:.F90=.F90)
	/bin/cp -f $< $(<F:.F90=.F90)
	$(FPP_COMMAND) $(<F:.F90=.F90)

ed_bigleaf_init.o : $(ED_INIT)/ed_bigleaf_init.f90
	/bin/rm -f $(<F:.f90=.f90)
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)

ed_cn_utils.o : $(ED_UTILS)/ed_cn_utils.f90
	/bin/rm -f $(<F:.f90=.f90)
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)

ed_driver.o : $(ED_DRIVER)/ed_driver.F90
	/bin/rm -f $(<F:.F90=.F90)
	/bin/cp -f $< $(<F:.F90=.F90)
	$(FPP_COMMAND) $(<F:.F90=.F90)

ed_filelist.o : $(ED_UTILS)/ed_filelist.F90
	/bin/rm -f $(<F:.F90=.F90) 
	/bin/cp -f $< $(<F:.F90=.F90)
	$(FPP_COMMAND) $(<F:.F90=.F90)

ed_grid.o : $(ED_UTILS)/ed_grid.f90
	/bin/rm -f $(<F:.f90=.f90) 
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)

ed_init.o : $(ED_INIT)/ed_init.F90
	/bin/rm -f $(<F:.F90=.F90)
	/bin/cp -f $< $(<F:.F90=.F90)
	$(FPP_COMMAND) $(<F:.F90=.F90)

ed_init_atm.o : $(ED_INIT)/ed_init_atm.F90
	/bin/rm -f $(<F:.F90=.F90) 
	/bin/cp -f $< $(<F:.F90=.F90)
	$(FPP_COMMAND) $(<F:.F90=.F90)

ed_init_history.o : $(ED_IO)/ed_init_history.f90
	/bin/rm -f $(<F:.f90=.f90)
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_LOWO_COMMAND) $(<F:.f90=.f90)

ed_load_namelist.o : $(ED_IO)/ed_load_namelist.f90
	/bin/rm -f $(<F:.f90=.f90)
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)

ed_max_dims.o : $(ED_MEMORY)/ed_max_dims.F90
	/bin/rm -f $(<F:.F90=.F90)
	/bin/cp -f $< $(<F:.F90=.F90)
	$(FPP_COMMAND) $(<F:.F90=.F90)

ed_mem_alloc.o : $(ED_MEMORY)/ed_mem_alloc.f90
	/bin/rm -f $(<F:.f90=.f90) 
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)

ed_met_driver.o : $(ED_DRIVER)/ed_met_driver.f90
	/bin/rm -f $(<F:.f90=.f90)
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(HDF5_INCS) $(<F:.f90=.f90)

ed_misc_coms.o : $(ED_MEMORY)/ed_misc_coms.f90
	/bin/rm -f $(<F:.f90=.f90)
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)

ed_model.o : $(ED_DRIVER)/ed_model.F90
	/bin/rm -f $(<F:.F90=.F90)
	/bin/cp -f $< $(<F:.F90=.F90)
	$(FPP_COMMAND) $(<F:.F90=.F90)

ed_mpass_init.o : $(ED_MPI)/ed_mpass_init.F90
	/bin/rm -f $(<F:.F90=.F90)
	/bin/cp -f $< $(<F:.F90=.F90)
	$(FPP_COMMAND) $(<F:.F90=.F90)

ed_nbg_init.o : $(ED_INIT)/ed_nbg_init.f90
	/bin/rm -f $(<F:.f90=.f90)
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)

ed_node_coms.o : $(ED_MPI)/ed_node_coms.f90
	/bin/rm -f $(<F:.f90=.f90) 
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)

ed_opspec.o : $(ED_IO)/ed_opspec.F90
	/bin/rm -f $(<F:.F90=.F90)
	/bin/cp -f $< $(<F:.F90=.F90)
	$(FPP_COMMAND) $(<F:.F90=.F90)

ed_params.o : $(ED_INIT)/ed_params.f90
	/bin/rm -f $(<F:.f90=.f90)
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)

ed_para_coms.o : $(ED_MPI)/ed_para_coms.f90
	/bin/rm -f $(<F:.f90=.f90) 
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)

ed_para_init.o : $(ED_MPI)/ed_para_init.f90
	/bin/rm -f $(<F:.f90=.f90)
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)

ed_print.o : $(ED_IO)/ed_print.F90
	/bin/rm -f $(<F:.F90=.F90)
	/bin/cp -f $< $(<F:.F90=.F90)
	$(FPP_COMMAND) $(<F:.F90=.F90)

ed_read_ed10_20_history.o : $(ED_IO)/ed_read_ed10_20_history.f90
	/bin/rm -f $(<F:.f90=.f90)
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)

ed_read_ed21_history.o : $(ED_IO)/ed_read_ed21_history.f90
	/bin/rm -f $(<F:.f90=.f90)
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)

ed_state_vars.o : $(ED_MEMORY)/ed_state_vars.F90
	/bin/rm -f $(<F:.F90=.F90)
	/bin/cp -f $< $(<F:.F90=.F90)
	$(FPP_LOWO_COMMAND) $(<F:.F90=.F90)

ed_therm_lib.o : $(ED_UTILS)/ed_therm_lib.f90
	/bin/rm -f $(<F:.f90=.f90)
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)

ed_type_init.o : $(ED_INIT)/ed_type_init.f90
	/bin/rm -f $(<F:.f90=.f90)
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_LOWO_COMMAND) $(<F:.f90=.f90)

ed_var_tables.o : $(ED_MEMORY)/ed_var_tables.f90 
	/bin/rm -f $(<F:.f90=.f90)
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)

ed_work_vars.o : $(ED_MEMORY)/ed_work_vars.f90
	/bin/rm -f $(<F:.f90=.f90)
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)

ed_xml_config.o : $(ED_IO)/ed_xml_config.f90
	/bin/rm -f $(<F:.f90=.f90)
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)

ename_coms.o : $(ED_MEMORY)/ename_coms.f90
	/bin/rm -f $(<F:.f90=.f90)
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)

euler_driver.o : $(ED_DYNAMICS)/euler_driver.f90
	/bin/rm -f $(<F:.f90=.f90)
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)

events.o : $(ED_DYNAMICS)/events.f90
	/bin/rm -f $(<F:.f90=.f90)
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)

farq_katul.o : $(ED_DYNAMICS)/farq_katul.f90
	/bin/rm -f $(<F:.f90=.f90)
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)

farq_leuning.o : $(ED_DYNAMICS)/farq_leuning.f90
	/bin/rm -f $(<F:.f90=.f90)
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)

fatal_error.o : $(ED_UTILS)/fatal_error.F90
	/bin/rm -f $(<F:.F90=.F90)
	/bin/cp -f $< $(<F:.F90=.F90)
	$(FPP_COMMAND) $(<F:.F90=.F90)

fire.o : $(ED_DYNAMICS)/fire.f90
	/bin/rm -f $(<F:.f90=.f90)
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)

forestry.o : $(ED_DYNAMICS)/forestry.f90
	/bin/rm -f $(<F:.f90=.f90)
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)

fusion_fission_coms.o : $(ED_MEMORY)/fusion_fission_coms.f90
	/bin/rm -f $(<F:.f90=.f90)
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)

fuse_fiss_utils.o : $(ED_UTILS)/fuse_fiss_utils.f90
	/bin/rm -f $(<F:.f90=.f90)
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)

great_circle.o : $(ED_UTILS)/great_circle.f90
	/bin/rm -f $(<F:.f90=.f90) 
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)

grid_coms.o : $(ED_MEMORY)/grid_coms.f90
	/bin/rm -f $(<F:.f90=.f90) 
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)

growth_balive.o : $(ED_DYNAMICS)/growth_balive.f90
	/bin/rm -f $(<F:.f90=.f90)
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)

h5_output.o : $(ED_IO)/h5_output.F90
	/bin/rm -f $(<F:.F90=.F90)
	/bin/cp -f $< $(<F:.F90=.F90)
	$(FPP_COMMAND) $(HDF5_INCS) $(<F:.F90=.F90)

hdf5_coms.o : $(ED_MEMORY)/hdf5_coms.f90
	/bin/rm -f $(<F:.f90=.f90)
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)

hdf5_utils.o : $(ED_UTILS)/hdf5_utils.f90
	/bin/rm -f $(<F:.f90=.f90)
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)

heun_driver.o: $(ED_DYNAMICS)/heun_driver.f90
	/bin/rm -f $(<F:.f90=.f90)
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)

hrzshade_utils.o: $(ED_UTILS)/hrzshade_utils.f90
	/bin/rm -f $(<F:.f90=.f90)
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)

hybrid_driver.o : $(ED_DYNAMICS)/hybrid_driver.f90
	/bin/rm -f $(<F:.f90=.f90)
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)

hydrology_coms.o: $(ED_MEMORY)/hydrology_coms.f90
	/bin/rm -f $(<F:.f90=.f90)
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)

hydrology_constants.o: $(ED_MEMORY)/hydrology_constants.f90
	/bin/rm -f $(<F:.f90=.f90)
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)

invmondays.o : $(ED_UTILS)/invmondays.f90
	/bin/rm -f $(<F:.f90=.f90)
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)

landuse_init.o : $(ED_INIT)/landuse_init.f90
	/bin/rm -f $(<F:.f90=.f90)
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)

lapse.o : $(ED_UTILS)/lapse.f90
	/bin/rm -f $(<F:.f90=.f90)
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)

leaf_database.o : $(ED_IO)/leaf_database.f90
	/bin/rm -f $(<F:.f90=.f90)
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)

libxml2f90.f90_pp.o : $(ED_UTILS)/libxml2f90.f90_pp.f90
	/bin/rm -f $(<F:.f90=.f90)
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)

lsm_hyd.o : $(ED_DYNAMICS)/lsm_hyd.f90
	/bin/rm -f $(<F:.f90=.f90)
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)

mem_polygons.o : $(ED_MEMORY)/mem_polygons.f90
	/bin/rm -f $(<F:.f90=.f90)
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)

met_driver_coms.o : $(ED_MEMORY)/met_driver_coms.f90
	/bin/rm -f $(<F:.f90=.f90)
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)

mortality.o : $(ED_DYNAMICS)/mortality.f90
	/bin/rm -f $(<F:.f90=.f90)
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)

multiple_scatter.o : $(ED_DYNAMICS)/multiple_scatter.f90
	/bin/rm -f $(<F:.f90=.f90)
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)

numutils.o: $(ED_UTILS)/numutils.f90
	/bin/rm -f $(<F:.f90=.f90)
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)

old_twostream_rad.o : $(ED_DYNAMICS)/old_twostream_rad.f90
	/bin/rm -f $(<F:.f90=.f90)
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)

pft_coms.o : $(ED_MEMORY)/pft_coms.f90
	/bin/rm -f $(<F:.f90=.f90)
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)

phenology_aux.o : $(ED_DYNAMICS)/phenology_aux.f90
	/bin/rm -f $(<F:.f90=.f90)
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)

phenology_coms.o : $(ED_MEMORY)/phenology_coms.f90
	/bin/rm -f $(<F:.f90=.f90)
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)

phenology_driv.o : $(ED_DYNAMICS)/phenology_driv.f90
	/bin/rm -f $(<F:.f90=.f90)
	/bin/cp -f $< $(<F:.f90=.f90) 
	$(F90_COMMAND) $(<F:.f90=.f90)

phenology_startup.o : $(ED_INIT)/phenology_startup.f90
	/bin/rm -f $(<F:.f90=.f90)
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)

photosyn_driv.o : $(ED_DYNAMICS)/photosyn_driv.f90
	/bin/rm -f $(<F:.f90=.f90)
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)

physiology_coms.o : $(ED_MEMORY)/physiology_coms.f90
	/bin/rm -f $(<F:.f90=.f90)
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)

plant_hydro.o : $(ED_DYNAMICS)/plant_hydro.f90
	/bin/rm -f $(<F:.f90=.f90)
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)

radiate_driver.o : $(ED_DYNAMICS)/radiate_driver.f90
	/bin/rm -f $(<F:.f90=.f90)
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)

radiate_utils.o : $(ED_UTILS)/radiate_utils.f90
	/bin/rm -f $(<F:.f90=.f90)
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)

random_utils.o : $(ED_UTILS)/random_utils.F90
	/bin/rm -f $(<F:.F90=.F90)
	/bin/cp -f $< $(<F:.F90=.F90)
	$(FPP_COMMAND) $(<F:.F90=.F90)

reproduction.o : $(ED_DYNAMICS)/reproduction.f90
	/bin/rm -f $(<F:.f90=.f90)
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)

rk4_coms.o : $(ED_MEMORY)/rk4_coms.f90
	/bin/rm -f $(<F:.f90=.f90)
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)

rk4_copy_patch.o : $(ED_DYNAMICS)/rk4_copy_patch.f90
	/bin/rm -f $(<F:.f90=.f90)
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)

rk4_derivs.o : $(ED_DYNAMICS)/rk4_derivs.f90
	/bin/rm -f $(<F:.f90=.f90)
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)

rk4_driver.o : $(ED_DYNAMICS)/rk4_driver.F90
	/bin/rm -f $(<F:.F90=.F90)
	/bin/cp -f $< $(<F:.F90=.F90)
	$(FPP_COMMAND) $(<F:.F90=.F90)

rk4_integ_utils.o : $(ED_DYNAMICS)/rk4_integ_utils.f90
	/bin/rm -f $(<F:.f90=.f90)
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)

rk4_misc.o : $(ED_DYNAMICS)/rk4_misc.f90
	/bin/rm -f $(<F:.f90=.f90)
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)

rsys.o: $(ED_UTILS)/rsys.F90
	/bin/rm -f $(<F:.F90=.F90)
	/bin/cp -f $< $(<F:.F90=.F90)
	$(FPP_COMMAND) $(<F:.F90=.F90)

soil_coms.o : $(ED_MEMORY)/soil_coms.F90
	/bin/rm -f $(<F:.F90=.F90)
	/bin/cp -f $< $(<F:.F90=.F90)
	$(FPP_COMMAND) $(<F:.F90=.F90)

soil_respiration.o : $(ED_DYNAMICS)/soil_respiration.f90
	/bin/rm -f $(<F:.f90=.f90)
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)

stable_cohorts.o : $(ED_UTILS)/stable_cohorts.f90
	/bin/rm -f $(<F:.f90=.f90)
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)

stem_resp_driv.o : $(ED_DYNAMICS)/stem_resp_driv.f90
	/bin/rm -f $(<F:.f90=.f90)
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)

structural_growth.o : $(ED_DYNAMICS)/structural_growth.f90
	/bin/rm -f $(<F:.f90=.f90)
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)

therm_lib.o: $(ED_UTILS)/therm_lib.f90
	/bin/rm -f $(<F:.f90=.f90)
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)

therm_lib8.o: $(ED_UTILS)/therm_lib8.f90
	/bin/rm -f $(<F:.f90=.f90)
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)

twostream_rad.o : $(ED_DYNAMICS)/twostream_rad.f90
	/bin/rm -f $(<F:.f90=.f90)
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)

update_derived_utils.o : $(ED_UTILS)/update_derived_utils.f90
	/bin/rm -f $(<F:.f90=.f90)
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)

utils_c.o: $(ED_UTILS)/utils_c.c
	/bin/rm -f $(<F:.c=.c)
	/bin/cp -f $< $(<F:.c=.c)
	$(CXX_COMMAND) $(<F:.c=.c)

utils_f.o: $(ED_UTILS)/utils_f.f90
	/bin/rm -f $(<F:.f90=.f90)
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)

vegetation_dynamics.o : $(ED_DYNAMICS)/vegetation_dynamics.f90
	/bin/rm -f $(<F:.f90=.f90)
	/bin/cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)


include dependency.mk
