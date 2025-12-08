#==========================================================================================#
#==========================================================================================#
#    Makefile include.mk.gfortran                                                          #
#                                                                                          #
#    Compilation controls for GNU-based Fortran/C.                                         #
#------------------------------------------------------------------------------------------#



#----- Define make (gnu make works best). -------------------------------------------------#
MAKE=/usr/bin/make
#------------------------------------------------------------------------------------------#



#----- Main path for compilation. ---------------------------------------------------------#
BASE=$(ED_ROOT)/build/
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#    HDF 5 libraries.                                                                      #
#                                                                                          #
#    Since ED-2.1, this is no longer optional.  You must provide HDF5 libraries, which     #
# have been compiled with the same compiler defined in F_COMP and C_COMP (see below).      #
#------------------------------------------------------------------------------------------#
HDF5_PATH=/path/to/hdf5-compiled-with-gcc-gfortran
HDF5_INCS=-I$(HDF5_PATH)/include
HDF5_LIBS= -lz -lm -L$(HDF5_PATH)/lib -lhdf5 -lhdf5_fortran -lhdf5_hl
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     If you have a version of hdf5 compiled in parallel, then you may benefit from        #
# collective I/O, then use this flag = 1.  Otherwise, set it to zero.                      #
#------------------------------------------------------------------------------------------#
USE_COLLECTIVE_MPIO=0
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Linear Algebra Package (LAPACK) libraries.                                           #
#                                                                                          #
#     Lapack is a well-established package for solving linear systems in Fortran. This is  #
# more efficient than the former built-in solution, and thus it became the new default.    #
#                                                                                          #
#     For those compiling the code with Intel compilers (ifort/icc or ifx/icx), leave      #
# these empty and compile the code with the Math Kernel Library option instead (-mkl or    #
# -qmkl depending on the ifort/icc version). Otherwise, provide the path to Lapack         #
# installation.                                                                            #
#------------------------------------------------------------------------------------------#
LAPACK_PATH=/usr/local/opt/lapack
LAPACK_INCS=-I$(LAPACK_PATH)/include
LAPACK_LIBS=-L$(LAPACK_PATH)/lib -llapack -lblas
#------------------------------------------------------------------------------------------#



#################################### MACHINE SETTINGS ######################################
#------------------------------------------------------------------------------------------#
#     CMACH. This tells the code which type of machine is being used. This is useful for   #
# setting up a few system-dependent pre-compilation instructions.  Current options include #
# the following:                                                                           #
#                                                                                          #
# LINUX  -- Linux or Unix systems.                                                         #
# MACOS  -- MacOS system. This applies more limited memory requests so the code can        #
#           compile and run on personal computers.                                         #
#                                                                                          #
#     We currently do not support compilation in Windows. If you know how to set up ED2 to #
# run on Windows machines, please submit a pull request.                                   #
#------------------------------------------------------------------------------------------#
CMACH=LINUX
#------------------------------------------------------------------------------------------#





#################################### COMPILER SETTINGS #####################################
#------------------------------------------------------------------------------------------#
#   FC_TYPE -- Specify from which family of compilers the code should be built. Current    #
#              options include:                                                            #
#              GNU       -- gfortran and gcc.                                              #
#              GFORTRAN  -- same as GNU.                                                   #
#              INTEL     -- ifx and icx (or ifort and icc in older systems)                #
#              IFORT     -- same as INTEL.                                                 #
#              PGI       -- pgf90 and pgcc.                                                #
#              PGF90     -- same as PGI.                                                   #
#   F_COMP  -- Fortran compiler. If empty, the code builder will automatically select the  #
#              default compiler names based on FC_TYPE. If you intend to compile ED2 with  #
#              openmpi, then set this to mpif90.                                           #
#   C_COMP  -- C compiler. If empty, the code builder will automatically select the        #
#              default compiler names based on FC_TYPE. If you intend to compile ED2 with  #
#              openmpi, then set this to mpicc.                                            #
#   LOADER  -- Fortran compiler to build the executable. By default, this is the same as   #
#              F_COMP.                                                                     #
#   LIBS    -- Libraries that should be linked to the compilation. Most of the time, this  #
#              should be left empty.                                                       #
#   MOD_EXT -- Suffix for external modules. This should be typically set to mod.           #
#------------------------------------------------------------------------------------------#
FC_TYPE=GNU
F_COMP=gfortran
C_COMP=gcc
LOADER=$(F_COMP)
LIBS=
MOD_EXT=mod
#------------------------------------------------------------------------------------------#





##################################### COMPILER OPTIONS #####################################
#------------------------------------------------------------------------------------------#
# A/B/C/D. Debugging, strictest compilation flags, lowest performance.                     #
# E.       Running, most relaxed compilation flags, highest performance.                   #
#------------------------------------------------------------------------------------------#
ifeq ($(KIND_COMP),)
   KIND_COMP=E
endif
ifeq ($(KIND_COMP),$(filter $(KIND_COMP), A B C D))
   F_OPTS=-g -ffree-line-length-none -fno-whole-file -O0 -fopenmp -ffpe-trap=invalid,zero,overflow -fbounds-check #-std=f2003
   C_OPTS=-fopenmp -g -O0 -ffpe-trap=invalid,zero,overflow -fbounds-check
   LOADER_OPTS=${F_OPTS}
else ifeq ($(KIND_COMP),E)
   F_OPTS=-ffree-line-length-none -fno-whole-file -O2 -fopenmp #-ffpe-trap=invalid,zero,overflow -fbounds-check  #-O2
   C_OPTS=-O2 -fopenmp #-ffpe-trap=invalid,zero,overflow -fbounds-check #-O2
   LOADER_OPTS=${F_OPTS}
else
   $(error Option KIND_COMP provided ($(KIND_COMP)) is invalid)."
endif
#------------------------------------------------------------------------------------------#
############################################################################################




#------------------------------------------------------------------------------------------#
#    MPI configuration. Most users should leave all of the settings below empty.           #
# The only exception is when users want to compile ED2 with OpenMPI/MPICH and, for         #
# whichever reason, they cannot set F_COMP=mpif90 and C_COMP=mpicc.                        #
#                                                                                          #
#   In case you truly need to set these variables, the commented-out example below shows   #
# how to set these variables:                                                              #
#                                                                                          #
# MPI_PATH=/path/to/openmpi                                                                #
# PAR_INCS=-I$(MPI_PATH)/include                                                           #
# PAR_LIBS=-L$(MPI_PATH)/lib -lmpi (-lmpi_f90 or -lmpi_usempif08 may be needed too).       #
#------------------------------------------------------------------------------------------#
MPI_PATH=
PAR_INCS=
PAR_LIBS=
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Archive options. This is typically "ar rs", so do not change it unless you know your #
# system has a different standard.                                                         #
#------------------------------------------------------------------------------------------#
ARCHIVE=ar rs
#------------------------------------------------------------------------------------------#
