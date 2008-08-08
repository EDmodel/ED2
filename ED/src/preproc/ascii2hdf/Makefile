#Makefile Make_model
##HDF5_DIR=/n/Moorcroft_Lab/Users/dmm2/hdf5/hdf5-1.6.5/hdf5
HDF5_DIR=/n/Moorcroft_Lab/Lab/apps/i91
F9X=/n/Moorcroft_Lab/Lab/apps/intel9.1/bin/ifort
CX=/n/Moorcroft_Lab/Lab/apps/intel9.1/bin/icc

all: convert_ascii.f90 hdf5_coms.f90 hdf5_utils_f.f90
	$(F9X) -c convert_ascii.f90
	$(F9X) -c hdf5_coms.f90 -I$(HDF5_DIR)/lib
	$(F9X) -c hdf5_utils_f.f90 -I$(HDF5_DIR)/lib
	$(F9X) -o convert_ascii convert_ascii.o hdf5_coms.o hdf5_utils_f.o -L$(HDF5_DIR)/lib -lhdf5 -lz -lm -lhdf5_fortran

