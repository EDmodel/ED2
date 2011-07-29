#------------------------------------- CUSTOMIZE HERE -------------------------------------#
here='pathhere/thispoly' # Folder to start the run
queue='thisqueue'                                   # Queue name
options=''                                          # Options, or leave it blank...
bsub='/lsf/7.0/linux2.6-glibc2.3-x86_64/bin/bsub'   # bsub, command to submit the job
joblog=${here}'/serial_lsf.out'                     # Name of the job output
jobname='thisdesc-thispoly'
callserial=${here}'/callserial.sh'                        # Name of executable
#------------------------------------------------------------------------------------------#


#----- Erasing old logfiles and joblogs ---------------------------------------------------#
if [ -s ${joblog} ]
then
  rm -fv ${joblog}
fi

#----- Submitting the job, we use a shell script so we can track the run on the fly -------#
if [ 'x'${options} == 'x' ]
then
   ${bsub} -q ${queue} -o ${joblog} -J ${jobname} -n 1 ${callserial} zzzzzzzz
else
   ${bsub} -q ${queue} -R ${options} -o ${joblog} -J ${jobname} -n 1 ${callserial} zzzzzzzz
fi
