#Makefile Make_model
NETCDF_INCS=-I/odyssey/apps/netcdf/4.0/netcdf-4.0/include \
            -I/odyssey/apps/netcdf/4.0/netcdf-4.0/lib
NETCDF_LIBS=-lnetcdf

HDF5_INCS=-I/n/sw/hdf5/hdf5-1.8.1_intel-10.0.015/include \
          -I/n/sw/hdf5/hdf5-1.8.1_intel-10.0.015/lib 
HDF5_LIBS=-lhdf5  -lm -lhdf5_fortran -lhdf5_hl -lz

ZLIB_DIR=/n/Moorcroft_Lab/Users/mlongo/util/hdflib/hdf4_64/zlib

F9X=ifort
F_OPTS=-g -traceback

CX=icc
C_OPTS=-g -traceback


FORCE:
	$(CX)  $(C_OPTS) -c read_ncep.c $(NETCDF_INCS) $(HDF5_INCS)
	$(CX)  $(C_OPTS) -o read_ncep read_ncep.o $(NETCDF_LIBS) $(HDF5_LIBS)
	$(F9X) $(F_OPTS) -c hdf5_coms.f90 $(NETCDF_INCS) $(HDF5_INCS)
	$(F9X) $(F_OPTS) -c hdf5_utils_f.f90 $(NETCDF_INCS) $(HDF5_INCS)
	$(F9X) $(F_OPTS) -c consts_coms.f90 
	$(F9X) $(F_OPTS) -c therm_lib.f90
	$(F9X) $(F_OPTS) -c convert_radiation_gen.f90
	$(F9X) $(F_OPTS) -o convert_radiation_gen hdf5_coms.o hdf5_utils_f.o consts_coms.o therm_lib.o convert_radiation_gen.o $(NETCDF_LIBS) $(HDF5_LIBS)


test:
	$(F9X) $(F_OPTS) -c test_coms.f90
	$(F9X) $(F_OPTS) -c test.f90
	$(F9X) $(F_OPTS) -o test test.o -L$(NETCDF_DIR)/lib -lnetcdf

clean:
	rm -fv *.o *.mod convert_radiation_gen read_ncep
