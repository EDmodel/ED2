#!/bin/sh
rpostroot=${1}
includes="-I ${rpostroot}/src/include"
rpostsrc="${rpostroot}/src/*/*.f90 ${rpostroot}/src/*/*.F90 ${rpostroot}/src/*/*.c"
rm -f dependency.mk
./sfmakedepend.pl ${includes} -f dependency.mk ${rpostsrc}
sed -i s@hdf5.mod@@g dependency.mk
