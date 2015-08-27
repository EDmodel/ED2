#!/bin/sh
bramsroot=${1}
edroot=${2}

includes="-I ${bramsroot}/src/utils/include"
bramssrc="${bramsroot}/src/*/*.f90 ${bramsroot}/src/*/*.F90 ${bramsroot}/src/*/*.c"
edsrc="${edroot}/src/*/*.f90 ${edroot}/src/*/*.F90 ${edroot}/src/*/*.c"
rm -f dependency.mk
./sfmakedepend.pl ${includes} -f dependency.mk ${bramssrc} ${edsrc}
sed -i~ s@hdf5.mod@@g   dependency.mk
sed -i~ s@netcdf.mod@@g dependency.mk
/bin/rm -f dependency.mk.old*
