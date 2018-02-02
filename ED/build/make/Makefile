###################################### Change Log ##########################################
#   4.3.0.1                                                                                #
#                                                                                          #
############################################################################################

#----- Define path and compilation --------------------------------------------------------#
include paths.mk
include include.mk
#------------------------------------------------------------------------------------------#


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

ifeq ($(OPT),"dbg")
	INCLUDES = ""
endif

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
EXE=$(BASE)/ed_$(ED_VERSION)-$(OPT)$(GIT_TAG)
LIBMODEL=$(BASE)/ed_$(ED_VERSION)-$(OPT)$(GIT_TAG).a
#------------------------------------------------------------------------------------------#

include objects.mk

#----- Define targets. --------------------------------------------------------------------#

#------------------------------------------------------------------------------------------#
#     MLO.  The six "make $(EXE)" below may look like a typo, but it is actually a rather  #
#           ugly solution for a problem with make. When someone changes a fortran file     #
#           that is a module, all files that use the module must be recompiled too. Some   #
#           files are updated in the first pass, but not all of them are recompiled, and   #
#           nested modules (a module used by another module) may require multiple          #
#           compilations so the final executable is consistent.  Failing to do so may      #
#           cause the executable to be inconsistent, with dangerous results.  If anyone    #
#           knows a better way to ensure all dependencies are recompiled, feel free to     #
#           implement, but for the time being please keep all the "make $(EXE)" calls.     #
#------------------------------------------------------------------------------------------#
all:
	make gendep
	make $(EXE)
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
	$(LOADER) -o $(EXE) edmain.o $(LOADER_OPTS) $(INCLUDES) $(LIBMODEL) $(HDF5_LIBS)   \
	$(PAR_LIBS) $(NC_LIBS) $(LIBS) $(LOADER_OPTS)
	@echo ""
	@echo Finished building === $(EXE)
	@echo ""

$(MAINOBJ):  $(MAIN) 
	@echo ""
	cp -f $< $(<F:.F90=.F90)
	$(FPP_COMMAND) $(<F:.F90=.F90)
	rm -f $(<F:.F90=.F90)

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
	rm -f $(LIBMODEL) $(EXE)
	rm -f ./$(BIN)/*.mk ./$(BIN)/*.sh ./$(BIN)/*.c ./$(BIN)/*.f90 ./$(BIN)/*.F90
	rm -f ./$(BIN)/*.mod ./$(BIN)/*.o ./$(BIN)/dependency.mk~
	rm -f ../$(EXE) ../$(LIBMODEL)
	rm -f include.* Makefile *.mk *.sh configure *.pl
	@echo ""

#----- Define rules -----------------------------------------------------------------------#
include rules.mk
