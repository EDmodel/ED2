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



#------------------------------------------------------------------------------------------#
#     Make variables CMACH and FC_TYPE are correctly set.                                  #
#------------------------------------------------------------------------------------------#
FATAL_ERROR=0
ifeq ($(CMACH),)
   CMACH_ERROR=1
   FATAL_ERROR=1
else ifeq ($(CMACH),$(filter $(CMACH), LINUX MACOS DOCKER CINTEG))
   CMACH_ERROR=0

   ifeq ($(FC_TYPE),$(filter $(FC_TYPE), GNU INTEL PGI))
      FCTYPE_ERROR=0
   else
      FCTYPE_ERROR=1
      FATAL_ERROR=1
   endif
else
   CMACH_ERROR=1
   FATAL_ERROR=1
endif
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#     Gracefully fail if CMACH and FC_TYPE are not correctly set. This will allow users to #
# make sure they selected a supported version.                                             #
#------------------------------------------------------------------------------------------#
ifeq ($(CMACH_ERROR),1)
   @echo ""
   @echo "---~---"
   @echo "   FATAL ERROR!"
   @echo "---~---"
   @echo "   Variable CMACH is either not defined or has an invalid value in your"
   @echo " include.mk.[PLATFORM]. Please set to one of the following values"
   @echo ""
   @echo " CMACH=LINUX  (For Linux machines)"
   @echo " CMACH=MACOS  (For MacOS machines)"
   @echo ""
   @echo "   It is possible to set CMACH to DOCKER and CINTEG, but these are special and"
   @echo " specific cases and should not be used unless you are working on these specific"
   @echo " cases."
   @echo "---~---"
   @echo ""
endif
ifeq ($(FCTYPE_ERROR),1)
      @echo ""
      @echo "---~---"
      @echo "   FATAL ERROR!"
      @echo "---~---"
      @echo "   Variable FC_TYPE is either not defined or has an invalid value in your"
      @echo " include.mk.[PLATFORM]. Please set to one of the following values"
      @echo ""
      @echo " FC_TYPE=GNU   # To compile ED2 with gfortran/gcc"
      @echo " FC_TYPE=INTEL # To compile ED2 with ifx/icx (or ifort/icc in older systems)"
      @echo " FC_TYPE=PGI   # To compile ED2 with pgfortran/pgcc"
      @echo ""
      @echo "    New types may be added in the future (feel free to contribute with other"
      @echo " compilers)."
      @echo "---~---"
      @echo ""
endif
ifeq ($(FATAL_ERROR),1)
   $(error CMACH and/or FC_TYPE are not properly set in your include.mk.[PLATFORM] file.)"
endif
#------------------------------------------------------------------------------------------#



#----- Set parallel flag based on the settings. -------------------------------------------#
ifeq ($(F_COMP),$(filter $(F_COMP), mpif90 mpifortran))
   PAR_DEFS=-DRAMS_MPI
else ifdef $(PAR_LIBS)
   PAR_DEFS=-DRAMS_MPI
else
   PAR_DEFS=
endif
#------------------------------------------------------------------------------------------#



#----- Compiler commands. -----------------------------------------------------------------#
INCLUDES         = $(PAR_INCS) -I$(ED_INCS) $(HDF5_INCS) $(MPI_INCS) $(LAPACK_INCS)
F90_COMMAND      = $(F_COMP) -c $(F_OPTS) $(INCLUDES) $(PAR_DEFS)
F90_LOWO_COMMAND = $(F_COMP) -c $(F_LOWO_OPTS) $(INCLUDES) $(PAR_DEFS)
FPP_COMMAND      = $(F_COMP) -c -D$(CMACH) -D$(FC_TYPE)                                    \
                   -DUSE_COLLECTIVE_MPIO=$(USE_COLLECTIVE_MPIO)                            \
                   -DUSE_MPIWTIME=$(USE_MPIWTIME) $(F_OPTS) $(INCLUDES) $(PAR_DEFS)
FPP_LOWO_COMMAND = $(F_COMP) -c -D$(CMACH) -D$(FC_TYPE)                                    \
                   -DUSE_COLLECTIVE_MPIO=$(USE_COLLECTIVE_MPIO)                            \
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
	$(LOADER) -o $(EXE) edmain.o $(LOADER_OPTS) $(INCLUDES) $(LIBMODEL) $(HDF5_LIBS)        \
	$(LAPACK_LIBS) $(PAR_LIBS) $(LIBS) $(LOADER_OPTS)
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
