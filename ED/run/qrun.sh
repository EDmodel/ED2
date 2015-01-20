#!/bin/sh

timestamp=`date +%y%m%d.%H%M%S`
logfile="/usr2/postdoc/rykelly/edoutputs/log/log.${timestamp}.txt"

qsub -v GFORTRAN_UNBUFFERED_ALL='y' -V -o ${logfile} -j y -N ED2 ~/ED2/ED/run/run.sh $1
#qsub -l h_rt=72:00:00 -v GFORTRAN_UNBUFFERED_ALL='y' -V -o ${logfile} -j y -N ED2 ./run.sh $1
