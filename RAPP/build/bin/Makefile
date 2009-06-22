#------------------------------------------------------------------------------------------#
# RAPP Makefile                                                                            #
# Compilation instructions to RAPP.                                                        #
#                                                                                          #
# WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING WARNING! WARNING!!#
# WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING WARNING! WARNING!!#
# WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING WARNING! WARNING!!#
#                                                                                          #
#    Unless you are developing the code, you should not change anything in here. You can   #
# set up the compilation to your needs by changing include.mk only!                        #
#                                                                                          #
# WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING WARNING! WARNING!!#
# WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING WARNING! WARNING!!#
# WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING WARNING! WARNING!!#
#------------------------------------------------------------------------------------------#

include include.mk

#----- Configuring paths based on the provided root folder --------------------------------#

RAPP_VERSION = 1.0

RAPP_BASE    = $(RAPP_ROOT)/build/

RAPP_DRIVER  = $(RAPP_ROOT)/src/driver
RAPP_HDF5    = $(RAPP_ROOT)/src/hdf5
RAPP_INTERP  = $(RAPP_ROOT)/src/interp
RAPP_IO      = $(RAPP_ROOT)/src/io
RAPP_MODULES = $(RAPP_ROOT)/src/modules
RAPP_NCEP    = $(RAPP_ROOT)/src/ncep
RAPP_UTILS   = $(RAPP_ROOT)/src/utils

#----- Define the executable name ---------------------------------------------------------#
EXE      = $(RAPP_BASE)/rapp_$(RAPP_VERSION)
LIBMODEL = $(RAPP_BASE)/rapp_$(RAPP_VERSION).a


#----- Compiler commands. -----------------------------------------------------------------#
INCLUDES    = $(RAPP_INCS) $(HDF5_INCS) $(NCDF_INCS)
ALL_LIBS    = $(LIBS) $(HDF5_LIBS) $(NCDF_LIBS)
F90_COMMAND = $(F_COMP) $(F_OPTS) $(INCLUDES) -c
FPP_COMMAND = $(F_COMP) -DUSE_NCDF=$(USE_NCDF) -D$(CMACH) -DUSE_HDF5=$(USE_HDF5)           \
              $(F_OPTS) $(INCLUDES) -c


#----- Including the object list ----------------------------------------------------------#
include objects.mk

#----- Defining targets -------------------------------------------------------------------#
all: $(EXE)

$(EXE): $(LIBMODEL) $(MAINOBJ)
	@echo ""
	$(LOADER) -o $(EXE) $(MAINOBJ) $(LOADER_OPTS) $(LIBMODEL) $(ALL_LIBS)
	@echo ""
	@echo Finished building === $(EXE)
	@echo ""

$(MAINOBJ): $(MAIN)
	@echo ""
	cp -f $< $(<F:.f90=.f90)
	$(LOADER) $(LOADER_OPTS) -c $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90) 

$(LIBMODEL): $(OBJECTS)
	$(ARCHIVE) $(LIBMODEL) $(OBJECTS)

FORCE: 

clean: 
	@echo ""
	rm -fv $(LIBMODEL) $(EXE) *.o *.mod *.F90 *.f90 *.stb *.d
	@echo ""

#------ Defining the dependencies ---------------------------------------------------------#
include depend.mk


