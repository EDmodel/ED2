#==========================================================================================#
#==========================================================================================#
#     Make file include.mk.docker.gnu                                                      #
#------------------------------------------------------------------------------------------#

#----- Define make (gnu make works best). -------------------------------------------------#
MAKE=/usr/bin/make
#------------------------------------------------------------------------------------------#

#----- Libraries. -------------------------------------------------------------------------#
BASE=$(ED_ROOT)/build/
#------------------------------------------------------------------------------------------#

#------------------------------------------------------------------------------------------#
#     HDF5 libraries                                                                       #
#                                                                                          #
#     Since ED-2.1, this is no longer optional for real simulations.  You must have HDF5   #
# libraries compiled with the same compiler you set for F_COMP and C_COMP.  You may still  #
# be able to compile without HDF5 but it will not run.                                     #
#------------------------------------------------------------------------------------------#
HDF5_INCS=-I/usr/include/hdf5/openmpi
HDF5_LIBS= -L/usr/lib/$(shell uname -m)-linux-gnu/hdf5/openmpi -lhdf5_fortran -lhdf5_hl -lhdf5 -lz -lm

#------------------------------------------------------------------------------------------#

#------------------------------------------------------------------------------------------#
#      If you have a version of hdf5 compiled in parallel, then you may benefit from       #
# collective I/O, then use this flag = 1.  Otherwise, set it to zero.                      #
#------------------------------------------------------------------------------------------#
USE_COLLECTIVE_MPIO=0
#------------------------------------------------------------------------------------------#

#------------------------------------------------------------------------------------------#
#      This should be 1 unless you are running with -gen-interfaces.  Interfaces usually   #
# make the compilation to crash when the -gen-interfaces option are on, so this flag       #
# bypass all interfaces in the code.                                                       #
#------------------------------------------------------------------------------------------#
USE_INTERF=1
#------------------------------------------------------------------------------------------#

#################################### COMPILER SETTINGS #####################################
CMACH=DOCKER_GNU
FC_TYPE=GNU
#F_COMP=gfortran
F_COMP=mpif90.openmpi
#C_COMP=gcc
C_COMP=mpicc.openmpi
#LOADER=gfortran
LOADER=mpif90.openmpi
#C_LOADER=gcc
C_LOADER=mpicc.openmpi
LIBS=
MOD_EXT=mod
############################################################################################





##################################### COMPILER OPTIONS #####################################
#------------------------------------------------------------------------------------------#
# A. Pickiest   - Use this whenever you change arguments on functions and subroutines.     #
#                 This will perform the same tests as B but it will also check whether all #
#                 arguments match between subroutine declaration and subroutine calls.     #
#                 WARNING: In order to really check all interfaces you must compile with   #
#                          this option twice:                                              #
#                 1. Compile (./install.sh A)                                              #
#                 2. Prepare second compilation(./2ndcomp.sh)                              #
#                 3. Compile one more time (./install.sh B)                                #
#                 If the compilation fails either at step 3, then your code has interface  #
#                 problems. If it successfully compiles, then the code is fine for         #
#                 interfaces.                                                              #
# E. Fast - This is all about performance, use only when you are sure that the model has   #
#           no code problem, and you want results asap. This will not check for any        #
#           problems, which means that this is an option suitable for end users, not de-   #
#           velopers.                                                                      #
#------------------------------------------------------------------------------------------#
ifeq ($(KIND_COMP),)
   KIND_COMP=E
endif
#------------------------------------------------------------------------------------------#
ifeq ($(KIND_COMP),A)
   F_OPTS= -O0 -ffree-line-length-none -g -fimplicit-none -Wall -finit-real=snan           \
           -finit-integer=-2147483648 -ffpe-trap=invalid,zero,overflow,underflow           \
           -fcheck=all -frecursive -fsignaling-nans -Werror -fopenmp -static
   C_OPTS= -O0 -DLITTLE  -g -static
   LOADER_OPTS= -O0 -ffree-line-length-none -g -fimplicit-none -Wall -finit-real=snan      \
                -finit-integer=-2147483648 -ffpe-trap=invalid,zero,overflow,underflow      \
                -fcheck=all -frecursive -fsignaling-nans -Werror -fopenmp
   #---------------------------------------------------------------------------------------#
endif
ifeq ($(KIND_COMP),E)
   F_OPTS= -O3 -ffree-line-length-none -frecursive -fopenmp -static
   C_OPTS= -O0 -DLITTLE  -g  -static
   LOADER_OPTS= -O3 -ffree-line-length-none -frecursive  -fopenmp
   #---------------------------------------------------------------------------------------#
endif
#------------------------------------------------------------------------------------------#
############################################################################################

#------------------------------------------------------------------------------------------#
#     If using mpicc and mpif90 as compilers (recommended), leave MPI_PATH, PAR_INCS, and  #
# PAR_LIBS blank, otherwise provide the includes and libraries for mpi.  Either way, don't #
# change PAR_DEFS.                                                                         #
#------------------------------------------------------------------------------------------#
MPI_PATH=
PAR_INCS=
PAR_LIBS=
PAR_DEFS=
#------------------------------------------------------------------------------------------#

#------------------------------------------------------------------------------------------#
#      Archive options.                                                                    #
#------------------------------------------------------------------------------------------#
#ARCHIVE=libtool -c -static -stack_size 0x1000000 -o
ARCHIVE=ar rs
#------------------------------------------------------------------------------------------#
