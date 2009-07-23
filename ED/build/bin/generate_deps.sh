#!/bin/sh
./sfmakedepend.pl -I ../../src/include -f Make.depend.1 ../../src/*/*.f90 ../../src/*/*.F90 ../../src/*/*.c
sed s/hdf5.mod// Make.depend.1 > Make.depend.2
sed s/leaf_coms.mod// Make.depend.2 > Make.depend.3
sed s/grid_dims.mod// Make.depend.3 > Make.depend.4
sed s/rconstants.mod// Make.depend.4 > Make.depend
rm Make.depend.1
rm Make.depend.2
rm Make.depend.3
rm Make.depend.4
