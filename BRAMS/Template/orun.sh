#!/bin/bash


#----  Change settings here ---------------------------------------------------------------#
here='thispath'                                       # Directory to start the run
queue='thisqueue'                           # Queue name
options=''                                       # Options, or leave it blank...
bsub='bsub'   # bsub, command to submit the job
joblog=${here}'/lsf.out'                         # Name of the job output
callopenmpi=${here}'/callopenmpi.sh'             # Script that calls MPI to run BRAMS
namelist='RAMSIN'                                # Namelist
span='12'                                        # Force span per node
nodespread='no'                                 # Force to spread accross as few nodes
                                                 # as possible
jobname='EB-thisjob'                      # Job name
numcore=60                                       # Number of cores
#------------------------------------------------------------------------------------------#





#----- Get the number of nodes from prompt if given... ------------------------------------#
if [ 'x'${1} != 'x' ]
then
 numcore=${1}
fi
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#    Check whether to spread cores to as few nodes as possible or go for first come first  #
# serve...                                                                                 #
#------------------------------------------------------------------------------------------#
if [ ${nodespread} == 'yes' ]
then
   nspcmd="-R span[ptile=${span}]"
else
   nspcmd=""
fi
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#    Check whether there is any option or not...                                           #
#------------------------------------------------------------------------------------------#
if [ 'x'${options} != 'x' ]
then
   optcmd="-R ${options}"
else
   optcmd=""
fi
#------------------------------------------------------------------------------------------#





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




#----- Erase old logfiles and joblogs -----------------------------------------------------#
if [ -s ${joblog} ]
then
  rm -fv ${joblog}
fi
#------------------------------------------------------------------------------------------#




#----- Check whether this is a MAKEVFILE/MAKESFC run or an actual run ---------------------#
makevfile=`grep "'MAKEVFILE'," ${namelist} | wc -l`
makesfc=`grep "'MAKESFC'," ${namelist}   | wc -l`
preponly=`expr ${makevfile} + ${makesfc}`
if [ ${preponly} -eq 1 ]
then
  numcore=1
  echo '+--------------------------------------------------------------------------------+'
  echo '+ MAKEVFILE or MAKESFC run, I will use only one core!                            +'
  echo '+--------------------------------------------------------------------------------+'
elif [ ${preponly} -eq 0 ]
then
  echo '+--------------------------------------------------------------------------------+'
  echo '+ INITIAL or HISTORY run, I will use '${numcore}' cores!                         +'
  echo '+--------------------------------------------------------------------------------+'
else
  echo '+--------------------------------------------------------------------------------+'
  echo '+ Something strange happened here... Check your RUNTYPE variable.                +'
  echo '+ Do not comment out lines with RUNTYPE, I feel confused when they exist, sorry! +'
  echo '+--------------------------------------------------------------------------------+'
  exit
fi
#------------------------------------------------------------------------------------------#




#----- Change the output filename to change for different processors. ---------------------#
if [ ${numcore} -lt 10 ]
then
  joblog=`dirname ${joblog}`'/00'${numcore}'_'`basename ${joblog}` 
elif [ ${numcore} -lt 100 ]
then
  joblog=`dirname ${joblog}`'/0'${numcore}'_'`basename ${joblog}` 
else
  joblog=`dirname ${joblog}`'/'${numcore}'_'`basename ${joblog}` 
fi
#------------------------------------------------------------------------------------------#





#----- This part will submit the MPI submitter script to LSF. -----------------------------#
command="${bsub} -q ${queue} ${optcmd} ${nspcmd} -o ${joblog}"
command="${command} ${jobcmd} -n ${numcore} ${callopenmpi}"
${command}
#------------------------------------------------------------------------------------------#
