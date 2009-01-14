###Makefile Make_model

## paths for Harvard
##HDF5_DIR=/n/Moorcroft_Lab/Users/dmm2/hdf5/hdf5-1.6.5/hdf5
#HDF5_DIR=/n/Moorcroft_Lab/Lab/apps/i91
#F9X=/n/Moorcroft_Lab/Lab/apps/intel9.1/bin/ifort
#CX=/n/Moorcroft_Lab/Lab/apps/intel9.1/bin/icc

## paths for UofI
HDF5_DIR=/home/mdietze/hdf5-1.8.1/hdf5
F9X=/usr/bin/gfortran
CX=/usr/bin/gcc
LOADER_OPTS= -g -Wl,-rpath -Wl,/home/mdietze/hdf5-1.8.1/hdf5/lib

all: convert_ascii.f90 hdf5_coms.f90 hdf5_utils_f.f90
	$(F9X) -c hdf5_coms.f90 -I$(HDF5_DIR)/include
	$(F9X) -c hdf5_utils_f.f90 -I$(HDF5_DIR)/include
	$(F9X) -c convert_ascii.f90
	$(F9X) -o convert_ascii $(LOADER_OPTS) convert_ascii.o hdf5_coms.o hdf5_utils_f.o -L$(HDF5_DIR)/lib -lhdf5 -lz -lm -lhdf5_fortran


