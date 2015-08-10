###################################### Change Log ##########################################
#   4.3.0.1                                                                                #
#                                                                                          #
############################################################################################

#----- Define path and compilation --------------------------------------------------------#
include paths.mk
include include.mk.$(OPT)


#------------------------------------------------------------------------------------------#
#      Double check that the "LOWO" flags have been set.  In case they have not, clone the #
# standard options.  LOWO stands for LOWer Optimisation, and these flags are used for a    #
# subroutines that are taking several hours to compile with ifort-13 (ed_state_vars.f90    #
# and a few others).                                                                       #
#------------------------------------------------------------------------------------------#
ifeq ($(F_LOWO_OPTS),)
   F_LOWO_OPTS = $(F_OPTS)
endif
#------------------------------------------------------------------------------------------#



#----- Compiler commands. -----------------------------------------------------------------#
INCLUDES         = $(PAR_INCS) -I$(ED_INCS) $(HDF5_INCS) $(MPI_INCS) 
F90_COMMAND      = $(F_COMP) -c $(F_OPTS) $(INCLUDES) $(PAR_DEFS)
F90_LOWO_COMMAND = $(F_COMP) -c $(F_LOWO_OPTS) $(INCLUDES) $(PAR_DEFS)
FPP_COMMAND      = $(F_COMP) -c -DUSE_INTERF=$(USE_INTERF) -DUSENC=$(USENC) -D$(CMACH)     \
                   -DUSE_HDF5=$(USE_HDF5) -DUSE_COLLECTIVE_MPIO=$(USE_COLLECTIVE_MPIO)     \
                   -DUSE_MPIWTIME=$(USE_MPIWTIME) $(F_OPTS) $(INCLUDES) $(PAR_DEFS)
FPP_LOWO_COMMAND = $(F_COMP) -c -DUSE_INTERF=$(USE_INTERF) -DUSENC=$(USENC) -D$(CMACH)     \
                   -DUSE_HDF5=$(USE_HDF5) -DUSE_COLLECTIVE_MPIO=$(USE_COLLECTIVE_MPIO)     \
                   -DUSE_MPIWTIME=$(USE_MPIWTIME) $(F_LOWO_OPTS) $(INCLUDES) $(PAR_DEFS)
CXX_COMMAND      = $(C_COMP) -c $(C_OPTS) -D$(CMACH) $(HDF5_INCS) $(INCLUDES) $(PAR_DEFS)
#------------------------------------------------------------------------------------------#




#----- Define archive and executable names. -----------------------------------------------#
EXE=$(BASE)/ed_$(ED_VERSION)-$(OPT)
LIBMODEL=$(BASE)/ed_$(ED_VERSION)-$(OPT).a
#------------------------------------------------------------------------------------------#


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
	./generate_deps.sh $(ED_ROOT)
	@echo === Finished dependencies ===

$(EXE): $(LIBMODEL) $(MAINOBJ) 
	@echo ""
	$(LOADER) -o $(EXE) edmain.o $(LOADER_OPTS) $(LIBMODEL) $(HDF5_LIBS) $(PAR_LIBS)   \
	$(NC_LIBS) $(LIBS) $(LOADER_OPTS)
	@echo ""
	@echo Finished building === $(EXE)
	@echo ""

$(MAINOBJ):  $(MAIN) 
	@echo ""
	/bin/rm -f $(<F:.F90=.F90)
	/bin/cp -f $< $(<F:.F90=.F90)
	$(FPP_COMMAND) $(<F:.F90=.F90)

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
	rm -f $(LIBMODEL) $(EXE) *.o *.mod *.F90 *.f90 *.stb *.d dependency.mk 
	rm -f ../$(EXE) ../$(LIBMODEL)
	touch dependency.mk
	@echo ""

#----- Define rules -----------------------------------------------------------------------#
include rules.mk
