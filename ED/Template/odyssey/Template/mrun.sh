#!/bin/bash
#-----Change your settings here  ----------------------------------------------------------#
here=`pwd` # Folder to start the run
queue='thisqueue'                                 # Queue name
#------------------------------------------------------------------------------------------#

if [ 'z'${1} = 'z' ]
then
  echo 'Which m-file do you want to run (needs to be in current dir include .m)?'
  read filenameIN
else
  filenameIN=${1}
fi

if [ 'z'${2} = 'z' ]
then
  echo 'What run name do you want to use?'
  read jobname
else
  jobname=${2}
fi


joberrname=${jobname}'_err.out'

module load math/matlab-R2009a 1>/dev/null 2>/dev/null

#----- Submitting the job, we use a shell script so we can track the run on the fly -------#
bsub -q ${queue}    -o ${joberrname} -J ${jobname} "matlab -nodesktop -nojvm -nosplash < ${filenameIN} > ${jobname}.out"

echo " - Running ${filenameIN} under job name ${jobname}"
