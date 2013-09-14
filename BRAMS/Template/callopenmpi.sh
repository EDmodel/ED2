#!/bin/bash
. ${HOME}/.bashrc i11

#------------------------------- CHANGE YOUR SETTINGS HERE. -------------------------------#
here='thispath'                                        # Folder to start the run
brams=${here}'/edbrams'                                # Executable
logfile=${here}'/out.out'                              # Log file
errfile=${here}'/out.err'                              # MPI error file
mpirun='mpirun'
#------------------------------------------------------------------------------------------#


#----- Find the number of cores. ----------------------------------------------------------#
nmach=`wc -l < ${LSB_DJOB_HOSTFILE}`

#----- Change the output filename to change for different processors ----------------------#
if [ ${nmach} -lt 10 ]
then
  logfile=`dirname ${logfile}`'/00'${nmach}'_'`basename ${logfile}` 
  errfile=`dirname ${errfile}`'/00'${nmach}'_'`basename ${errfile}` 
elif [ ${nmach} -lt 100 ]
then
  logfile=`dirname ${logfile}`'/0'${nmach}'_'`basename ${logfile}` 
  errfile=`dirname ${errfile}`'/0'${nmach}'_'`basename ${errfile}` 
else
  logfile=`dirname ${logfile}`'/'${nmach}'_'`basename ${logfile}` 
  errfile=`dirname ${errfile}`'/'${nmach}'_'`basename ${errfile}` 
fi

#----- Erase old logfiles and joblogs -----------------------------------------------------#
if [ -s ${logfile} ]
then
  rm -fv ${logfile}
fi
if [ -s ${errfile} ]
then
  rm -fv ${errfile}
fi

#----- Submit--- the jobs to the nodes ----------------------------------------------------#
#${mpirun} -np ${nmach} --nooversubscribe -machinefile ${LSB_DJOB_HOSTFILE}  -mca btl_openib_ib_timeout 20 -mca btl openib,sm ${brams} 1> ${logfile} 2> ${errfile}
${mpirun} -np ${nmach} --nooversubscribe -machinefile ${LSB_DJOB_HOSTFILE}   ${brams} 1> ${logfile} 2> ${errfile}
#${mpirun} -np ${nmach} --nooversubscribe -machinefile ${LSB_DJOB_HOSTFILE}  -mca btl_openib_ib_timeout 20 -mca btl self,tcp ${brams} 1> ${logfile} 2> ${errfile}

