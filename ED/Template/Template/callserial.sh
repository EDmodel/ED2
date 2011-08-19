#!/bin/sh

#------------------------------------- CUSTOMIZE HERE -------------------------------------#
root='thisroot'
here=${root}/'thispoly'                        # Directory to start the run
exe=${here}'/myexec'                               # Executable
logfile=${here}'/serial_out.out'               # Log file of the executable run
errfile=${here}'/serial_out.err'               # Executable error file
currloc=`pwd`                                  # Current location
#------------------------------------------------------------------------------------------#


#----- Erasing old logfiles and joblogs ---------------------------------------------------#
if [ -s ${logfile} ]
then
  rm -fv ${logfile}
fi

if [ -s ${errfile} ]
then
  rm -fv ${errfile}
fi

#----- Submitting the jobs to the nodes ---------------------------------------------------#
cd ${here}
${exe} -f ${here}/ED2IN 1> ${logfile} 2> ${errfile}
cd ${currloc}
