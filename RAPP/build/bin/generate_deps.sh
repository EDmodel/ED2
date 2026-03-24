#!/bin/sh
rapproot=${1}
includes="-I ${rapproot}/src/include"
rappsrc="${rapproot}/src/*/*.f90 ${rapproot}/src/*/*.F90 ${rapproot}/src/*/*.c"
/bin/rm -f dependency.mk
./sfmakedepend.pl ${includes} -f dependency.mk ${rappsrc} 2> /dev/null
sed -i~ s@hdf5.mod@@g dependency.mk
sed -i~ s@ifport.mod@@g dependency.mk
sed -i~ 's@[^_]netcdf.mod@ @g' dependency.mk
/bin/rm -f dependency.alt dependency.mk.old*

