# DO NOT DELETE THIS LINE - used by make depend
rcio.o: leaf_coms.mod micro_coms.mod rconstants.mod rpost_coms.mod somevars.mod
rcio.o: therm_lib.mod
rpost_main.o: brams_data.mod leaf_coms.mod misc_coms.mod rpost_coms.mod
rpost_main.o: rpost_dims.mod
rpost_misc.o: misc_coms.mod rpost_dims.mod
variables.o: an_header.mod brams_data.mod leaf_coms.mod micro_coms.mod
variables.o: misc_coms.mod rconstants.mod rpost_coms.mod rpost_dims.mod
variables.o: scratch_coms.mod somevars.mod
comp_lib.o: leaf_coms.mod rconstants.mod rpost_coms.mod soil_coms.mod
comp_lib.o: somevars.mod therm_lib.mod
dted.o: /n/Moorcroft_Lab/Users/mlongo/EDBRAMS/Ramspost/src/include/utils_sub_names.h
dted.o:
eenviron.o: /n/Moorcroft_Lab/Users/mlongo/EDBRAMS/Ramspost/src/include/utils_sub_names.h
eenviron.o:
interp_lib.o: misc_coms.mod
leaf_coms.o: rconstants.mod therm_lib.mod
parlib.o: /n/Moorcroft_Lab/Users/mlongo/EDBRAMS/Ramspost/src/include/utils_sub_names.h
parlib.o:
rams_read_header.o: an_header.mod
rnamel.o: misc_coms.mod
therm_lib.o: rconstants.mod
tmpname.o: /n/Moorcroft_Lab/Users/mlongo/EDBRAMS/Ramspost/src/include/utils_sub_names.h
tmpname.o:
vformat_brams3.3.o: misc_coms.mod
brams_data.o: rpost_dims.mod
micro_coms.o: rpost_dims.mod
misc_coms.o: rpost_dims.mod
rpost_coms.o: rpost_dims.mod
soil_coms.o: rpost_dims.mod
dateutils.o: rconstants.mod
numutils.o: rconstants.mod therm_lib.mod
polarst.o: rconstants.mod
utils_c.o: /n/Moorcroft_Lab/Users/mlongo/EDBRAMS/Ramspost/src/include/utils_sub_names.h
utils_c.o:
an_header.mod: an_header.o
brams_data.mod: brams_data.o
leaf_coms.mod: leaf_coms.o
micro_coms.mod: micro_coms.o
misc_coms.mod: misc_coms.o
rconstants.mod: rconstants.o
rpost_coms.mod: rpost_coms.o
rpost_dims.mod: rpost_dims.o
scratch_coms.mod: scratch_coms.o
soil_coms.mod: soil_coms.o
somevars.mod: somevars.o
therm_lib.mod: therm_lib.o
