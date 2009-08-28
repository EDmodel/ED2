an_header.o : $(RAPP_MODULES)/an_header.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

charutils.o : $(RAPP_UTILS)/charutils.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

dateutils.o : $(RAPP_UTILS)/dateutils.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

dealloc_driver.o : $(RAPP_DRIVER)/dealloc_driver.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

ed_metd_header.o : $(RAPP_IO)/ed_metd_header.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

fatal_error.o : $(RAPP_UTILS)/fatal_error.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

great_circle.o : $(RAPP_UTILS)/great_circle.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

hdf5_coms.o : $(RAPP_MODULES)/hdf5_coms.F90
	cp -f $< $(<F:.F90=.F90)
	$(FPP_COMMAND) $(<F:.F90=.F90)
	rm -f $(<F:.F90=.F90)

hdf5_utils.o : $(RAPP_UTILS)/hdf5_utils.F90
	cp -f $< $(<F:.F90=.F90)
	$(FPP_COMMAND) $(<F:.F90=.F90)
	rm -f $(<F:.F90=.F90)

interp_driver.o : $(RAPP_INTERP)/interp_driver.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

load_namelist.o : $(RAPP_IO)/load_namelist.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

mod_grid.o : $(RAPP_MODULES)/mod_grid.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

mod_interp.o : $(RAPP_MODULES)/mod_interp.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

mod_ioopts.o : $(RAPP_MODULES)/mod_ioopts.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

mod_maxdims.o : $(RAPP_MODULES)/mod_maxdims.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

mod_model.o : $(RAPP_MODULES)/mod_model.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

mod_namelist.o : $(RAPP_MODULES)/mod_namelist.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

mod_ncdf_globio.o: $(RAPP_MODULES)/mod_ncdf_globio.F90
	cp -f $< $(<F:.F90=.F90)
	$(FPP_COMMAND) $(<F:.F90=.F90)
	rm -f $(<F:.F90=.F90)

mod_ncep.o : $(RAPP_MODULES)/mod_ncep.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

mod_netcdf.o : $(RAPP_MODULES)/mod_netcdf.F90
	cp -f $< $(<F:.F90=.F90)
	$(FPP_COMMAND) $(<F:.F90=.F90)
	rm -f $(<F:.F90=.F90)

mod_time.o : $(RAPP_MODULES)/mod_time.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

ncepcio.o : $(RAPP_NCEP)/ncepcio.F90
	cp -f $< $(<F:.F90=.F90)
	$(FPP_COMMAND) $(<F:.F90=.F90)
	rm -f $(<F:.F90=.F90)

ncep_alloc.o : $(RAPP_NCEP)/ncep_alloc.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

ncep_coordinates.o: $(RAPP_NCEP)/ncep_coordinates.f90 
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

ncep_fill_infotable.o: $(RAPP_NCEP)/ncep_fill_infotable.F90
	cp -f $< $(<F:.F90=.F90)
	$(FPP_COMMAND) $(<F:.F90=.F90)
	rm -f $(<F:.F90=.F90)

ncep_loadvars.o: $(RAPP_NCEP)/ncep_loadvars.F90 
	cp -f $< $(<F:.F90=.F90)
	$(FPP_COMMAND) $(<F:.F90=.F90)
	rm -f $(<F:.F90=.F90)

ncep_output.o : $(RAPP_NCEP)/ncep_output.F90 
	cp -f $< $(<F:.F90=.F90)
	$(FPP_COMMAND) $(<F:.F90=.F90)
	rm -f $(<F:.F90=.F90)

numutils.o : $(RAPP_UTILS)/numutils.f90 
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

rain_downscale.o : $(RAPP_INTERP)/rain_downscale.f90 
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

rapp_coord.o : $(RAPP_DRIVER)/rapp_coord.f90 
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

rapp_driver.o : $(RAPP_DRIVER)/rapp_driver.f90 
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

rapp_opspec.o : $(RAPP_IO)/rapp_opspec.f90 
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

rconstants.o : $(RAPP_MODULES)/rconstants.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

rsys.o: $(RAPP_UTILS)/rsys.F90
	cp -f $< $(<F:.F90=.F90)
	$(FPP_COMMAND) $(<F:.F90=.F90)
	rm -f $(<F:.F90=.F90)

space_interp.o : $(RAPP_INTERP)/space_interp.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

therm_lib.o: $(RAPP_UTILS)/therm_lib.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

time_interp.o : $(RAPP_INTERP)/time_interp.f90 
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

utils_f.o: $(RAPP_UTILS)/utils_f.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

include dependency.mk
