#$ -S /bin/bash
#$ -q thisqueue
#$ -o pathhere/thispoly
#$ -N thisdesc-thispoly
#$ -j y
#$ -r n
#---------------------------------- Change settings here ----------------------------------#
root="thisroot"                           # Main directory
moi="myname"                              # User's account
here="${root}/thispoly"                   # Directory to start the run
exe="${here}/myexec"                      # Executable
initrc="myinitrc"                         # Script to load before doing anything
logfile="${here}/serial_out.out"          # Log file of the executable run
errfile="${here}/serial_out.err"          # Executable error file
currloc=$(pwd)                            # Current location
mddir="met_driver"                        # Path with met drivers
scendir="myscenmain"                      # Path with scenarios
datasrc="mypackdata"                      # Path where the original data is
datadest="/scratch/${moi}"                # Scratch area where the data will be copied
scenario="myscenario"                     # Which scenario to use
naptime=zzzzzzzz                          # Nap time
#------------------------------------------------------------------------------------------#


#----- Source script. ---------------------------------------------------------------------#
. ${initrc} -i
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Check whether there is rsync in this machine or not.                                 #
#------------------------------------------------------------------------------------------#
rsync_here=$(which rsync 2> /dev/null)
if [ "x${rsync_here}" == "x" ]
then
   cp_here=$(which cp | grep -v alias | sed s@"\t"@""@g)
   rsync="${cp_here} -ruva"
else
   rsync="${rsync_here} -Pruvaz"
fi
#------------------------------------------------------------------------------------------#


#----- Erase old logfiles and joblogs -----------------------------------------------------#
if [ -s ${logfile} ]
then
  rm -fv ${logfile}
fi

if [ -s ${errfile} ]
then
  rm -fv ${errfile}
fi
#------------------------------------------------------------------------------------------#


#----- Find the waiting time before the first check.  If none was given, make up one... ---#
if [ "x${naptime}" == "x999999" ]
then
  zzz=0
  copy="n"

elif [ "x${naptime}" == "x" ]
then
   zzz=$(date +%S)
   if [ ${zzz} -lt 10 ]
   then
      zzz=$(echo ${zzz} | awk '{print substr($1,2,1)}')
   fi
   let zzz=${zzz}%15
   let zzz=${zzz}+2
   copy="y"
else
   zzz=${naptime}
   copy="y"
fi
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Check size based on the met driver.                                                  #
#------------------------------------------------------------------------------------------#
if [ ${scenario} == "sheffield" ]
then
   datasize=39000000
else
   datasize=300000
fi
#------------------------------------------------------------------------------------------#


#----- These are the names of the lock and unlock files. ----------------------------------#
dontcopy=${datadest}/dontcopy_${scenario}.txt
diskfull=${datadest}/fulldisk_${scenario}.txt
unlocked=${datadest}/unlocked_${scenario}.txt
#------------------------------------------------------------------------------------------#


#----- These variables are the scenario paths (global and node). --------------------------#
source_thisscen="${datasrc}/${mddir}/${scendir}/${scenario}"
node_driverroot="${datadest}/${mddir}"
node_scenroot="${node_driverroot}/${scendir}"
node_thisscen="${node_scenroot}/${scenario}"
#------------------------------------------------------------------------------------------#


#----- Make the code more "personal", so it is easy to spot where the runs are running. ---#
thismach=$(hostname -s)
echo "------------------------------------------------------" 1>> ${logfile} 2>> ${errfile}
echo " Script callserial.sh starts on node ${thismach}..."    1>> ${logfile} 2>> ${errfile}
echo "------------------------------------------------------" 1>> ${logfile} 2>> ${errfile}
echo " "                                                      1>> ${logfile} 2>> ${errfile}
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#    Decide whether to copy files or not.                                                  #
#------------------------------------------------------------------------------------------#
if [ ${copy} == "n" ]
then
   success="true"
else
   success="false"

   #----- Take a break before checking, so we avoid several jobs going simultaneously. ----#
   blah="Waiting ${zzz} seconds before checking for the first time..."
   echo ${blah} 1>> ${logfile} 2>> ${errfile}
   sleep ${zzz}
   #---------------------------------------------------------------------------------------#

   if [ ! -s ${dontcopy} ] && [ ! -s ${diskfull} ]
   then
      #----- Don't remove the directory if it is there... ---------------------------------#
      if [ ! -s ${datadest}          ]
      then
         mkdir ${datadest}
      fi
      if [ ! -s ${node_driverroot}   ]
      then
         mkdir ${node_driverroot}
      fi
      if [ ! -s ${node_scenroot} ]
      then
         mkdir ${node_scenroot}
      fi
      #------------------------------------------------------------------------------------#

      echo "Copying files to node ${thismach}... " > ${dontcopy}

      blah=" + Files are not all there, cleaning path before we copy again."
      blah="${blah}  Please hold!"
      echo ${blah} 1>> ${logfile} 2>> ${errfile}

      if [ -s ${node_thisscen} ]
      then
         blah="  - Deleting old stuff from the meterological forcing driver..."
         echo ${blah} 1>> ${logfile} 2>> ${errfile}
         /bin/rm -fr ${node_thisscen} 1>> ${logfile} 2>> ${errfile}
      fi

      #----- First thing we do is to check whether this disk is full. ---------------------#
      blah="  - Checking whether there is enough disk space..."
      echo ${blah} 1>> ${logfile} 2>> ${errfile}
      ans=$(df /scratch)
      nwords=$(echo ${ans} | wc -w)
      let avail=${nwords}-2
      space=$(echo ${ans} | awk '{print $'${avail}'}')

      if [ ${space} -gt ${datasize} ]
      then

         #----- Copy the meteorological forcing. ------------------------------------------#
         blah="  - Copying the meterological forcing driver..."
         echo ${blah} 1>> ${logfile} 2>> ${errfile}
         ${rsync} ${source_thisscen} ${node_scenroot} 1>> ${logfile} 2>> ${errfile}


         #----- Copy finished.  Create a file to unlock this node. ------------------------#
         echo "All the data needed are here!" > ${unlocked}

         blah=" + The files were successfully copied to ${thismach}."
         blah="${blah}  The run will start soon!"
         echo ${blah} 1>> ${logfile} 2>> ${errfile}
         
         success="true"
      else
         blah="  - Sorry, but there is not enough disk space here..."
         echo ${blah} 1>> ${logfile} 2>> ${errfile}

         echo "The available disk size is not enough to run here..." > ${diskfull}
         /bin/rm -f ${dontcopy}
      fi

   elif [ -s ${diskfull} ]
   then
      last=$(date +%s -r ${diskfull})
      now=$(date +%s)
      let howlong=${now}-${last}

      if [ ${howlong} -gt 43200 ]
      then

         echo "Copying files to node ${thismach}... " > ${dontcopy}

         blah=" + Files are not all here, cleaning the directory before we copy again."
         blah="${blah}  Please hold!"
         echo ${blah} 1>> ${logfile} 2>> ${errfile}

         if [ -s ${node_thisscen} ]
         then
            blah="  - Deleting old stuff from the meterological forcing driver..."
            echo ${blah} 1>> ${logfile} 2>> ${errfile}
            /bin/rm -fr ${node_thisscen} 1>> ${logfile} 2>> ${errfile}
         fi

         #----- First thing we do is to check whether this disk is full. ------------------------#
         blah="  - Checking whether there is enough disk space..."
         echo ${blah} 1>> ${logfile} 2>> ${errfile}
         ans=$(df /scratch)
         nwords=$(echo ${ans} | wc -w)
         let avail=${nwords}-2
         space=$(echo ${ans} | awk '{print $'${avail}'}')

         if [ ${space} -gt ${datasize} ]
         then

            #----- Copy the meteorological forcing. ---------------------------------------#
            blah="  - Copying the meterological forcing driver..."
            echo ${blah} 1>> ${logfile} 2>> ${errfile}
            ${rsync} ${source_thisscen} ${node_scenroot} 1>> ${logfile} 2>> ${errfile}


            #----- Copy finished.  Create a file to unlock this node. ---------------------#
            echo "All the data needed are here!" > ${unlocked}

            blah=" + The files were successfully copied to ${thismach}."
            blah="${blah}  The run will start soon!"
            echo ${blah} 1>> ${logfile} 2>> ${errfile}
            
            success="true"
         else
            blah="  - Sorry, but there is still not enough disk space here..."
            echo ${blah} 1>> ${logfile} 2>> ${errfile}

            echo "The available disk size is not enough to run here..." > ${diskfull}
            /bin/rm -f ${dontcopy}
         fi
      fi

   elif [ ! -s ${unlocked} ]
   then

      pp=-1
      while [ ! -s ${unlocked} ] && [ ! -s ${diskfull} ]
      do
         let pp=${pp}+1
         
         #----- Variety is the key when it comes to automated messages... --------------------#
         let isfour=${pp}%4
         let iseight=${pp}%8
         let istwelve=${pp}%12
         let issixty=${pp}%60

         #----- Print some messages to entertain the bored user. -----------------------------#
         if   [ ${pp} -eq 0 ]
         then 
            blah=" ~ Another job is copying the files to here (${thismach})."
            blah="${blah}  Please hold!"
         elif [ ${issixty} -eq 0 ]
         then
            blah=" + Patience is a virtue, isn't it?"
         elif [ ${istwelve} -eq 0 ]
         then
            blah=" + Your simulation is very important to us.  Please hold!"
         elif [ ${iseight} -eq 0 ]
         then
            blah=" + Thank you for choosing ${thismach}.  Please hold!"
         elif [ ${isfour} -eq 0 ]
         then
            blah=" + All our representatives are busy, but we will be with you shortly."
            blah="${blah}  Please hold!"
         else
            blah=" - Input data are not ready yet... (${pp} minutes waiting so far)"
         fi

         echo ${blah} 1>> ${logfile} 2>> ${errfile}

         sleep 60
      done

      if [ -s ${unlocked} ]
      then
         blah=" + Thank you for waiting.  The files are here and the run will start soon!"
         echo ${blah} 1>> ${logfile} 2>> ${errfile}
         success='true'
      else
         blah=" + Sorry... There was not enough disk space here."
         echo ${blah} 1>> ${logfile} 2>> ${errfile}
      fi
      #------------------------------------------------------------------------------------#
   else
      success="true"
   fi
   #---------------------------------------------------------------------------------------#
fi
#------------------------------------------------------------------------------------------#



#----- Start the run. ---------------------------------------------------------------------#
if [ ${success} == "true" ]
then
   blah="Waiting ${zzz} seconds before starting the run..."
   echo ${blah} 1>> ${logfile} 2>> ${errfile}
   sleep ${zzz}

   #----- The time we were waiting. -------------------------------------------------------#
   cd ${here}
   ${exe} -f ${here}/ED2IN 1>> ${logfile} 2>> ${errfile}
   cd ${currloc}
else
   blah="Execution will not happen this time... Good luck next time!!!"
   echo ${blah} 1>> ${logfile} 2>> ${errfile}
fi
#------------------------------------------------------------------------------------------#
