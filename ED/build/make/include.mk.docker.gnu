#==========================================================================================#
#==========================================================================================#
#     Makefile include.mk.docker.gnu                                                       #
#                                                                                          #
#    Compilation controls for GNU-based Fortran/C for building the Docker container.       #
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
HDF5_INCS=-I/usr/include/hdf5/openmpi
HDF5_LIBS= -L/usr/lib/$(shell uname -m)-linux-gnu/hdf5/openmpi -lhdf5_fortran -lhdf5_hl -lhdf5 -lz -lm
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Linear Algebra Package (LAPACK) libraries.                                           #
#                                                                                          #
#     Lapack is a well-established package for solving linear systems in Fortran. This is  #
# more efficient than the former built-in solution, and thus it became the new default.    #
#                                                                                          #
#     For those compiling the code with Intel compilers (ifx/icx, or ifort/icc in older    #
# systems), leave these empty and compile the code with the Math Kernel Library option     #
# instead (-mkl or -qmkl depending on the ifx/ifort version). Otherwise, provide the path  #
# to Lapack installation.                                                                  #
#------------------------------------------------------------------------------------------#
LAPACK_PATH=/usr
LAPACK_INCS=-I$(LAPACK_PATH)/include
LAPACK_LIBS=-L$(LAPACK_PATH)/lib -llapack -lblas
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#      If you have a version of hdf5 compiled in parallel, then you may benefit from       #
# collective I/O, then use this flag = 1.  Otherwise, set it to zero.                      #
#------------------------------------------------------------------------------------------#
USE_COLLECTIVE_MPIO=0
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
# DOCKER -- Specific setting for building the containerised version of ED2. Do not use     #
#           this for regular builds.                                                       #
# CINTEG -- Specific setting for building the continuous integration version of ED2. Do    #
#           not use this for regular builds.                                               #
#                                                                                          #
#     We currently do not support compilation in Windows. If you know how to set up ED2 to #
# run on Windows machines, please submit a pull request!.                                  #
#------------------------------------------------------------------------------------------#
CMACH=DOCKER
#------------------------------------------------------------------------------------------#


#################################### COMPILER SETTINGS #####################################
#------------------------------------------------------------------------------------------#
#   FC_TYPE -- Specify from which family of compilers the code should be built. Current    #
#              options include:                                                            #
#              GNU       -- gfortran and gcc.                                              #
#              INTEL     -- ifx and icx (or ifort and icc in older systems)                #
#              PGI       -- pgfortran and pgcc.                                            #
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
F_COMP=mpif90.openmpi
C_COMP=mpicc.openmpi
LOADER=$(F_COMP)
LIBS=
MOD_EXT=mod
#------------------------------------------------------------------------------------------#





##################################### COMPILER OPTIONS #####################################
#------------------------------------------------------------------------------------------#
# A/B/C/D. Debugging, strictest compilation flags, lowest performance.                     #
# E.       Running, most relaxed compilation flags, highest performance.                   #
#------------------------------------------------------------------------------------------#
ifeq ($(KIND_COMP),$(filter $(KIND_COMP), A B C D))
   F_OPTS= -O0 -ffree-line-length-none -g -fimplicit-none -Wall -finit-real=snan           \
           -finit-integer=-2147483648 -ffpe-trap=invalid,zero,overflow,underflow           \
           -fcheck=all -frecursive -fsignaling-nans -Werror -fopenmp -static
   C_OPTS= -O0 -DLITTLE  -g -static
   LOADER_OPTS= -O0 -ffree-line-length-none -g -fimplicit-none -Wall -finit-real=snan      \
                -finit-integer=-2147483648 -ffpe-trap=invalid,zero,overflow,underflow      \
                -fcheck=all -frecursive -fsignaling-nans -Werror -fopenmp
   #---------------------------------------------------------------------------------------#
else ifeq ($(KIND_COMP),E)
   F_OPTS= -O3 -ffree-line-length-none -frecursive -fopenmp -static
   C_OPTS= -O0 -DLITTLE  -g  -static
   LOADER_OPTS= -O3 -ffree-line-length-none -frecursive  -fopenmp
   #---------------------------------------------------------------------------------------#
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
