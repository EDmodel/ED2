
# Define file/object containing the main program
MAIN    = $(RAPP_DRIVER)/rapp_main.f90
MAINOBJ = rapp_main.o

# Define objects:
OBJECTS =                       \
	an_header.o             \
	charutils.o             \
	dateutils.o             \
	fatal_error.o           \
	load_namelist.o         \
	mod_grid.o              \
	mod_ioopts.o            \
	mod_maxdims.o           \
	mod_model.o             \
	mod_namelist.o          \
	mod_ncdf_globio.o       \
	mod_ncep.o              \
	mod_netcdf.o            \
	mod_time.o              \
	numutils.o              \
	ncepcio.o               \
	ncep_coordinates.o      \
	ncep_fill_infotable.o   \
	rapp_driver.o           \
	rapp_opspec.o           \
	rconstants.o            \
	rsys.o                  \
	therm_lib.o             \
	utils_f.o
