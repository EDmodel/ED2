# Include the objects.mk file first to define $(SRC_MAIN) and $(OBJ_MODEL)
include ED/build/shell/sources.mk

SRC_MAIN = ED/src/driver/edmain.F90

AUTOMAKE_OPTIONS = subdir-objects

# Define the main program and static library variables
MAIN = ed2
STATICLIB = libed2.a

# Define the programs and libraries
bin_PROGRAMS = $(MAIN)
lib_LIBRARIES = $(STATICLIB)

# Compiler and linker flags
AM_CFLAGS = @CFLAGS@
AM_FCFLAGS = @CFLAGS@
AM_LDFLAGS = @LIBS@

# Define source files for the main program
ed2_SOURCES = $(SRC_MAIN)
ed2_LDADD = $(STATICLIB)

# Define source files for the static library
libed2_a_SOURCES = $(SOURCES)

# Additional Fortran clean files
CLEANFILES = *.mod

# Additional files for distribution
EXTRA_DIST = ED/src/include EDR EDTS R-utils

# Include the dependency makefile
include ED/build/shell/dependency.mk
