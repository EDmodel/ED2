#--------------------------------- Change settings here -----------------------------------#
here='pathhere/thispoly'                            # Folder to start the run
queue='thisqueue'                                   # Queue name
options=''                                          # Options, or leave it blank...
bsub='/lsf/7.0/linux2.6-glibc2.3-x86_64/bin/bsub'   # bsub, command to submit the job
joblog=${here}'/serial_lsf.out'                     # Name of the job output
jobname='thisdesc-thispoly'                         # Job name
callserial=${here}'/callserial.sh'                  # Name of executable
thisnum=myorder
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#    Special queues, make sure they use the least amount of nodes.                         #
#------------------------------------------------------------------------------------------#
if [ ${queue} == "wofsy" ]
then
   if [ ${thisnum} -le 12 ]
   then
      options='hname=wofsy011'
   elif [ ${thisnum} -le 24 ]
   then
      options='hname=wofsy012'
   elif [ ${thisnum} -le 36 ]
   then
      options='hname=wofsy013'
   fi
elif [ ${queue} == "camd" ]
then
   if [ ${thisnum} -le 48 ]
   then
      options='hname=camd06'
   elif [ ${thisnum} -le 96 ]
   then
      options='hname=camd08'
   fi
fi
#------------------------------------------------------------------------------------------#


#----- Erase old logfiles and joblogs -----------------------------------------------------#
if [ -s ${joblog} ]
then
  rm -fv ${joblog}
fi
#------------------------------------------------------------------------------------------#



#----- Submit the job, we use a shell script so we can track the run on the fly -----------#
if [ 'x'${options} == 'x' ]
then
   ${bsub} -q ${queue} -o ${joblog} -J ${jobname} -n 1 ${callserial} zzzzzzzz
else
   ${bsub} -q ${queue} -R ${options} -o ${joblog} -J ${jobname} -n 1 ${callserial} zzzzzzzz
fi
#------------------------------------------------------------------------------------------#
