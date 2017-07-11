#!/bin/bash
. ${HOME}/.bashrc

#------------------------------------- CUSTOMISE HERE -------------------------------------#
here='myoutpath/rpost'                                     # Path to start the run
queue='thisqueue'                                   # Queue name
options=''                                          # Options, or leave it blank...
bsub='/lsf/7.0/linux2.6-glibc2.3-x86_64/bin/bsub'   # bsub, command to submit the job
joblog=${here}'/serial_lsf.out'                     # Name of the job output
callserial=${here}'/1eachtime-sigma.sh'             # Name of executable
jobname='RP-thissim'                               # Job name
#------------------------------------------------------------------------------------------#


#----- Erasing old logfiles and joblogs ---------------------------------------------------#
if [ -s ${joblog} ]
then
  rm -fv ${joblog}
fi

#------------------------------------------------------------------------------------------#
#    Check whether there is job name or not...                                             #
#------------------------------------------------------------------------------------------#
if [ 'x'${jobname} != 'x' ]
then
   jobcmd="-J ${jobname}"
else
   jobcmd=""
fi
#------------------------------------------------------------------------------------------#




#----- Submitting the job, we use a shell script so we can track the run on the fly -------#
if [ 'x'${options} == 'x' ]
then
   ${bsub} ${jobcmd} -q ${queue} -o ${joblog} -n 1 ${callserial}
else
   ${bsub} ${jobcmd} -q ${queue} -R ${options} -o ${joblog} -n 1 ${callserial}
fi
