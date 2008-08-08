#Makefile Make_model
NETCDF_DIR=/n/Moorcroft_Lab/Lab/apps
HDF5_DIR=/n/Moorcroft_Lab/Lab/apps/i91
##HDF5_DIR=/n/Moorcroft_Lab/Users/dmm2/hdf5/hdf5-1.6.5/hdf5
ZLIB_DIR=/n/Moorcroft_Lab/Users/mlongo/util/hdflib/hdf4_64/zlib
F9X=/n/Moorcroft_Lab/Lab/apps/intel9.1/bin/ifort
CX=/n/Moorcroft_Lab/Lab/apps/intel9.1/bin/icc
F_OPTS=-g -traceback
C_OPTS=-g -traceback


FORCE:
	$(CX)  $(C_OPTS) -c read_ncep.c -I$(NETCDF_DIR)/include
	$(CX)  $(C_OPTS) -o read_ncep read_ncep.o -L$(NETCDF_DIR)/lib -lnetcdf
	$(F9X) $(F_OPTS) -c hdf5_coms.f90 -I$(HDF5_DIR)/lib
	$(F9X) $(F_OPTS) -c hdf5_utils_f.f90 -I$(HDF5_DIR)/lib
	$(F9X) $(F_OPTS) -c convert_radiation_gen.f90
	$(F9X) $(F_OPTS) -o convert_radiation_gen hdf5_coms.o hdf5_utils_f.o convert_radiation_gen.o -L$(HDF5_DIR)/lib -lhdf5  -lm -lhdf5_fortran -L$(ZLIB_DIR)/lib -lz


test:
	$(F9X) $(F_OPTS) -c test_coms.f90
	$(F9X) $(F_OPTS) -c test.f90
	$(F9X) $(F_OPTS) -o test test.o -L$(NETCDF_DIR)/lib -lnetcdf

clean:
	rm -fv *.o *.mod convert_radiation_gen read_ncep
