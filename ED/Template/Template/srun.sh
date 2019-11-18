#!/bin/bash
#--------------------------------- Change settings here -----------------------------------#
here="pathhere/thispoly"                  # Folder to start the run
queue="thisqueue"                         # Queue name
joblog="${here}/serial_lsf.out"           # Name of the job output
jobname="thisdesc-thispoly"               # Job name
callserial="${here}/callserial.sh"        # Name of executable
initrc="myinitrc"                         # Script to load before doing anything
thisnum=myorder                           # No longer used
sbatch=$(which sbatch)                    # SLURM command to submit job.
memory=thismemory                         # Requested memory (per CPU)
runtime=thistime                          # Requested time
#------------------------------------------------------------------------------------------#


#----- Source script. ---------------------------------------------------------------------#
. ${initrc}
#------------------------------------------------------------------------------------------#


#----- Erase old logfiles and joblogs -----------------------------------------------------#
if [ -s ${joblog} ]
then
  rm -fv ${joblog}
fi
#------------------------------------------------------------------------------------------#



#----- Submit the job, we use a shell script so we can track the run on the fly -----------#
${sbatch} -p ${queue} --mem-per-cpu=${memory} -t ${runtime} -o ${joblog} -J ${jobname}     \
          -n 1 --wrap="${callserial} zzzzzzzz"
#------------------------------------------------------------------------------------------#
