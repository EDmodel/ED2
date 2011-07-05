#!/bin/bash

#----- Load totalview. --------------------------------------------------------------------#
module load devel/totalview-8.8.0-1   1> /dev/null  2> /dev/null

#----- Get the executable name. -----------------------------------------------------------#
if [ 'x'${1} == 'x' ]
then
   /bin/ls --color=auto
   echo -n ' Enter the debugger executable: '
   read debugexe
else
   debugexe=${1}
fi

#----- Find the number of cores. ----------------------------------------------------------#
nmach=`wc -l < ${LSB_DJOB_HOSTFILE}`


#---- Run the parallel debugger. ----------------------------------------------------------#
mpirun -machinefile ${LSB_DJOB_HOSTFILE} --nooversubscribe -np ${nmach} --debug ${debugexe}
