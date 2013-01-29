###################################### Change Log ##########################################
#   4.3.0.1                                                                                #
#                                                                                          #
############################################################################################

#----- Define path and compilation --------------------------------------------------------#
include paths.mk
include include.mk.$(OPT)

#----- Compiler commands. -----------------------------------------------------------------#
INCLUDES  = -I$(UTILS_INCS) $(PAR_INCS) $(HDF5_INCS) $(MPI_INCS)
F90_COMMAND = $(F_COMP) -c $(F_OPTS) $(INCLUDES) $(PAR_DEFS)
FPP_COMMAND = $(F_COMP) -c  -DUSE_INTERF=$(USE_INTERF) -DUSENC=$(USENC) -D$(CMACH)         \
               -DUSE_HDF4=$(USE_HDF4) -DUSE_HDF5=$(USE_HDF5)                               \
               -DUSE_MPIWTIME=$(USE_MPIWTIME) -DCOUPLED $(F_OPTS) $(INCLUDES) $(PAR_DEFS)  \
               -DUSE_COLLECTIVE_IO=$(USE_COLLECTIVE_MPIO)
CXX_COMMAND = $(C_COMP) -c $(C_OPTS) -D$(CMACH) $(INCLUDES) $(PAR_DEFS)

#----- Define archive and executable names. -----------------------------------------------#
EXE=$(BASE)/edbrams-$(OPT)
LIBMODEL=$(BASE)/edbrams-$(OPT).a

include objects.mk

#----- Define targets. --------------------------------------------------------------------#
all: 
	make gendep
	make $(EXE)
	make $(EXE)
	make $(EXE)
	make $(EXE)
	make $(EXE)

gendep:
	@echo ""
	`pwd`/generate_deps.sh $(BRAMS_ROOT) $(ED_ROOT)
	@echo Finished dependencies

$(EXE): $(LIBMODEL) $(MAINOBJ)
	@echo ""
	$(LOADER) -o $(EXE) rammain.o $(LOADER_OPTS) $(LIBMODEL) \
	$(HDF4_LIBS) $(HDF5_LIBS) $(PAR_LIBS) $(NC_LIBS) $(LIBS)
	@echo ""
	@echo Finished building === $(EXE)
	@echo ""

$(MAINOBJ):  $(MAIN) 
	@echo ""
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90) 

$(LIBMODEL): $(OBJ_MODEL)
	$(ARCHIVE) $(LIBMODEL) $(OBJ_MODEL)


FORCE: 

install:
	@echo ""
	ln -fs `pwd`/$(EXE) ../run/$(BASE)
	ln -fs `pwd`/$(EXE) ../test/$(BASE)
	@echo ""

clean:
	@echo ""
	#(cd ./utils ; $(MAKE) -f Make.utils.$(OPT) OPT=$(OPT) clean)
	rm -f $(LIBMODEL) $(EXE) *.o *.mod *.f *.f90 *.c *.F90 *.stb *.d dependency.mk 
	rm -f ../$(EXE) ../$(LIBMODEL)
	touch dependency.mk
	@echo ""

#----- Define rules -----------------------------------------------------------------------#
include rules.mk

