#!/bin/sh
edroot=${1}
includes="-I ${edroot}/src/include"
edsrc="${edroot}/src/*/*.f90 ${edroot}/src/*/*.F90 ${edroot}/src/*/*.c"
rm -f dependency.mk
./sfmakedepend.pl ${includes} -f dependency.mk ${bramssrc} ${edsrc}
sed -i s@hdf5.mod@@g dependency.mk
sed -i s@leaf.coms.mod@@g dependency.mk
sed -i s@grid_dims.mod@@g dependency.mk
sed -i s@rconstants.mod@@g dependency.mk
/bin/rm -f dependency.alt dependency.mk.old*
