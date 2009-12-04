#------------------------------------------------------------------------------------------#
# Ramspost Makefile                                                                        #
# Compilation instructions to Ramspost.                                                    #
#                                                                                          #
# WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING WARNING! WARNING!!#
# WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING WARNING! WARNING!!#
# WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING WARNING! WARNING!!#
#                                                                                          #
#    Unless you are developing the code, you should not change anything in here. You can   #
# set up the compilation to your needs by changing include_ramspost.mk only!               #
#                                                                                          #
# WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING WARNING! WARNING!!#
# WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING WARNING! WARNING!!#
# WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING WARNING! WARNING!!#
#------------------------------------------------------------------------------------------#
include include_ramspost.mk

#----- Configuring paths based on the provided root folder --------------------------------#
RPOST_VERSION = 6.2

RPOST_BASE    = $(RPOST_ROOT)/build/

RPOST_DRIVER  = $(RPOST_ROOT)/src/driver
RPOST_INCLUDE = $(RPOST_ROOT)/src/include
RPOST_MEMORY  = $(RPOST_ROOT)/src/memory
RPOST_LIB     = $(RPOST_ROOT)/src/lib
RPOST_UTILS   = $(RPOST_ROOT)/src/utils


#----- Define the executable name ---------------------------------------------------------#
EXE      = $(RPOST_BASE)/ramspost_$(RPOST_VERSION)
LIBMODEL = $(EXE).a

INCLUDES    = -I$(RPOST_INCLUDE) $(MPI_INCS)
ALL_LIBS    = $(LIBS) $(MPI_LIBS)
F90_COMMAND = $(F_COMP) $(F_OPTS) $(INCLUDES) -c
CXX_COMMAND = $(C_COMP) -c $(C_OPTS) -D$(CMACH) $(INCLUDES) $(PAR_DEFS)


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
	`pwd`/generate_deps.sh $(RPOST_ROOT)
	@echo Finished dependencies

$(EXE): $(LIBMODEL) $(MAINOBJ)
	@echo ""
	$(LOADER) -o $(EXE) $(MAINOBJ) $(LOADER_OPTS) $(INCLUDES) $(LIBMODEL) $(ALL_LIBS)
	@echo ""
	@echo Finished building === $(EXE)
	@echo ""

$(MAINOBJ):  $(MAIN) 
	@echo ""
	cp -f $< $(<F:.f90=.f90)
	$(LOADER) $(LOADER_OPTS) $(INCLUDES) -c  $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90) 

$(LIBMODEL): $(OBJECTS)
	$(ARCHIVE) $(LIBMODEL) $(OBJECTS)


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
	rm -f $(EXE) $(LIBMODEL)
	touch dependency.mk
	@echo ""

#----- Define rules -----------------------------------------------------------------------#
include rules.mk

