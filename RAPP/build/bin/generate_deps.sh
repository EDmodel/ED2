#!/bin/sh
rapproot=${1}
includes="-I ${rapproot}/src/include"
rappsrc="${rapproot}/src/*/*.f90 ${rapproot}/src/*/*.F90 ${rapproot}/src/*/*.c"
rm -f dependency.mk
./sfmakedepend.pl ${includes} -f dependency.mk ${rappsrc}
cat dependency.mk   | sed s@\ hdf5.mod@\ @g        > dependency.alt
cat dependency.alt  | sed s@\ netcdf.mod@\ @g      > dependency.mk
/bin/rm -f dependency.alt dependency.mk.old*
