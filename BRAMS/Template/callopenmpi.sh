#!/bin/sh
#BSUB -u mlongo@fas.harvard.edu
#BSUB -J edbrams
#BSUB -o hello_lsf.out
#BSUB -e hello_lsf.err
#BSUB -n 1024
#BSUB -R "span[ptile=8]"

#module add hpc/openmpi-intel
#module add hpc/hdf5-intel

#------------------------------- CHANGE YOUR SETTINGS HERE. -------------------------------#
here=`pwd`                                             # Folder to start the run
brams=${here}'/edbrams'                            # Executable
logfile=${here}'/out.out'                              # Log file
errfile=${here}'/out.err'                              # MPI error file
mpirun='mpirun'
#------------------------------------------------------------------------------------------#



#----- Checking available spots -----------------------------------------------------------#
host=/odyssey/home/mlongo/.lsbatch/hosts.${LSB_JOBID}


/bin/cp ${LSB_DJOB_HOSTFILE} ${host}
let NPROC=0
let NNODE=0
let NCORE=0
for machine in `cat ${LSB_DJOB_HOSTFILE}`
do
   let NPROC=${NPROC}+1
   let NCORE=${NCORE}+1
   if [ ${NCORE} -eq 8 ]
   then
      let NCORE=0
      let NNODE=${NNODE}+1
   fi
done

#----- Changing the output filename to change for different processors --------------------#
if [ ${NPROC} -lt 10 ]
then
  logfile=`dirname ${logfile}`'/00'${NPROC}'_'`basename ${logfile}` 
  errfile=`dirname ${errfile}`'/00'${NPROC}'_'`basename ${errfile}` 
elif [ ${NPROC} -lt 100 ]
then
  logfile=`dirname ${logfile}`'/0'${NPROC}'_'`basename ${logfile}` 
  errfile=`dirname ${errfile}`'/0'${NPROC}'_'`basename ${errfile}` 
else
  logfile=`dirname ${logfile}`'/'${NPROC}'_'`basename ${logfile}` 
  errfile=`dirname ${errfile}`'/'${NPROC}'_'`basename ${errfile}` 
fi
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
#${mpirun} -np ${NPROC} --nooversubscribe -machinefile ${host}  -mca btl_openib_ib_timeout 20 -mca btl openib,sm ${brams} 1> ${logfile} 2> ${errfile}
${mpirun} -np ${NPROC} --nooversubscribe -machinefile ${host}   ${brams} 1> ${logfile} 2> ${errfile}
#${mpirun} -np ${NPROC} --nooversubscribe -machinefile ${host}  -mca btl_openib_ib_timeout 20 -mca btl self,tcp ${brams} 1> ${logfile} 2> ${errfile}


rm -f ${host}
