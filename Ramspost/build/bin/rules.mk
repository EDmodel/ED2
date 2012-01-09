abort_run.o: $(RPOST_UTILS)/abort_run.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

an_header.o: $(RPOST_MEMORY)/an_header.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

brams_data.o: $(RPOST_MEMORY)/brams_data.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

charutils.o: $(RPOST_UTILS)/charutils.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

comp_lib.o: $(RPOST_LIB)/comp_lib.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

dateutils.o: $(RPOST_UTILS)/dateutils.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

dted.o: $(RPOST_LIB)/dted.c 
	cp -f $< $(<F:.c=.c)
	$(CXX_COMMAND) $< $(F:.c=.c)
	rm -f $(<F:.c=.c)

eenviron.o: $(RPOST_LIB)/eenviron.c 
	cp -f $< $(<F:.c=.c)
	$(CXX_COMMAND) $< $(F:.c=.c)
	rm -f $(<F:.c=.c)

error_mess.o: $(RPOST_UTILS)/error_mess.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

interp_lib.o: $(RPOST_LIB)/interp_lib.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

leaf_coms.o: $(RPOST_LIB)/leaf_coms.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

micro_coms.o: $(RPOST_MEMORY)/micro_coms.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

misc_coms.o: $(RPOST_MEMORY)/misc_coms.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

numutils.o: $(RPOST_UTILS)/numutils.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

parlib.o: $(RPOST_LIB)/parlib.c 
	cp -f $< $(<F:.c=.c)
	$(CXX_COMMAND) $< $(F:.c=.c)
	rm -f $(<F:.c=.c)

polarst.o: $(RPOST_UTILS)/polarst.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

rams_read_header.o: $(RPOST_LIB)/rams_read_header.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

rfvec.o: $(RPOST_LIB)/rfvec.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

rcio.o: $(RPOST_DRIVER)/rcio.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

rconstants.o: $(RPOST_MEMORY)/rconstants.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

rgrad.o: $(RPOST_DRIVER)/rgrad.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

rnamel.o: $(RPOST_LIB)/rnamel.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

rnumr.o: $(RPOST_LIB)/rnumr.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

rpost_coms.o: $(RPOST_MEMORY)/rpost_coms.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

rpost_dims.o: $(RPOST_MEMORY)/rpost_dims.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

rpost_filelist.o: $(RPOST_UTILS)/rpost_filelist.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

rpost_misc.o: $(RPOST_DRIVER)/rpost_misc.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

rsys.o: $(RPOST_UTILS)/rsys.F90
	cp -f $< $(<F:.F90=.F90)
	$(FPP_COMMAND) $(<F:.F90=.F90)
	rm -f $(<F:.F90=.F90)

rutil.o: $(RPOST_UTILS)/rutil.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

scratch_coms.o: $(RPOST_MEMORY)/scratch_coms.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

soil_coms.o: $(RPOST_MEMORY)/soil_coms.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

somevars.o: $(RPOST_MEMORY)/somevars.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

therm_lib.o: $(RPOST_LIB)/therm_lib.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

tmpname.o: $(RPOST_LIB)/tmpname.c 
	cp -f $< $(<F:.c=.c)
	$(CXX_COMMAND) $< $(F:.c=.c)
	rm -f $(<F:.c=.c)

utils_c.o: $(RPOST_UTILS)/utils_c.c 
	cp -f $< $(<F:.c=.c)
	$(CXX_COMMAND) $< $(F:.c=.c)
	rm -f $(<F:.c=.c)

utils_f.o: $(RPOST_UTILS)/utils_f.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

variables.o: $(RPOST_DRIVER)/variables.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

vformat_brams3.3.o: $(RPOST_LIB)/vformat_brams3.3.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

include dependency.mk
