#!/bin/bash

#---------------------------------- Change settings here ----------------------------------#
root="thisroot"                  # Main directory
moi="myname"                     # User's account
here="${root}/thispoly"          # Directory to start the run
exec="myexec"                    # Executable
logfile="${here}/serial_out.out" # Log file of the executable run
errfile="${here}/serial_out.err" # Executable error file
currloc=$(pwd)                   # Current location
mddir="met_driver"               # Name of met driver root directory
scendir="myscenmain"             # Root directory of met driver type
datasrc="mypackdata"             # Source
datadest="/scratch/${moi}"       # Destination
scenario="myscenario"            # Actual scenario
n_cpt=mycpus                     # Number of CPUs per task (used only if not automatic)
copy=mycopy                      # Copy to scratch (true / false)
#----- Initialisation scripts. ------------------------------------------------------------#
optsrc="myoptsrc"                # Option for .bashrc (for special submission settings)
                                 #   In case none is needed, leave it blank ("").
#------------------------------------------------------------------------------------------#


#----- Source script. ---------------------------------------------------------------------#
case "x${HOME}" in
   x) HOME=$(echo ~) ;;
esac
. ${HOME}/.bashrc ${optsrc}
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
if ${copy}
then
   let zzz=${RANDOM}%60+1
else
   let zzz=${RANDOM}%10+1
fi
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Check size based on the met driver.                                                  #
#------------------------------------------------------------------------------------------#
case "${scenario}" in
sheffield)   datasize=39000000   ;;
WFDEI*)      datasize=27500000   ;;
ERAINT*)     datasize=21000000   ;;
MERRA2*)     datasize=57000000   ;;
*)           datasize=300000     ;;
esac
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
echo " Script callserial.sh starts on node ${thismach}."      1>> ${logfile} 2>> ${errfile}
echo " "                                                      1>> ${logfile} 2>> ${errfile}
echo " Settings:"                                             1>> ${logfile} 2>> ${errfile}
echo " Job:             ${SLURM_JOB_NAME} (${SLURM_JOB_ID})"  1>> ${logfile} 2>> ${errfile}
echo " Queue:           ${SLURM_JOB_PARTITION}"               1>> ${logfile} 2>> ${errfile}
echo " Node list:       ${SLURM_NODELIST}"                    1>> ${logfile} 2>> ${errfile}
echo " Number of nodes: ${SLURM_NNODES}"                      1>> ${logfile} 2>> ${errfile}
echo " Number of tasks: ${SLURM_NTASKS}"                      1>> ${logfile} 2>> ${errfile}
echo " Memory per CPU:  ${SLURM_MEM_PER_CPU}"                 1>> ${logfile} 2>> ${errfile}
echo " CPUs per task:   ${SLURM_CPUS_PER_TASK}"               1>> ${logfile} 2>> ${errfile}
echo " ED-2 output:     ${logfile}"                           1>> ${logfile} 2>> ${errfile}
echo " ED-2 error:      ${errfile}"                           1>> ${logfile} 2>> ${errfile}
echo " "                                                      1>> ${logfile} 2>> ${errfile}
echo " Copy files to scratch: ${copy}"                        1>> ${logfile} 2>> ${errfile}
echo "------------------------------------------------------" 1>> ${logfile} 2>> ${errfile}
echo " "                                                      1>> ${logfile} 2>> ${errfile}
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#    Decide whether to copy files or not.                                                  #
#------------------------------------------------------------------------------------------#
if ${copy}
then
   success=false

   #----- Take a break before checking, so we avoid several jobs going simultaneously. ----#
   blah="Wait ${zzz} seconds before checking for the first time."
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

      echo "Copy files to node ${thismach}. " > ${dontcopy}

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
      blah="  - Check whether there is enough disk space."
      echo ${blah} 1>> ${logfile} 2>> ${errfile}
      ans=$(df /scratch)
      nwords=$(echo ${ans} | wc -w)
      let avail=${nwords}-2
      space=$(echo ${ans} | awk '{print $'${avail}'}')

      if [ ${space} -gt ${datasize} ]
      then

         #----- Copy the meteorological forcing. ------------------------------------------#
         blah="  - Copy the meterological forcing driver."
         echo ${blah} 1>> ${logfile} 2>> ${errfile}
         ${rsync} ${source_thisscen} ${node_scenroot} 1>> ${logfile} 2>> ${errfile}


         #----- Copy finished.  Create a file to unlock this node. ------------------------#
         echo "All the data needed are here!" > ${unlocked}

         blah=" + The files have been successfully copied to ${thismach}."
         blah="${blah}  The run will start soon!"
         echo ${blah} 1>> ${logfile} 2>> ${errfile}
         
         success=true
      else
         blah="  - Sorry, but there is not enough disk space here."
         echo ${blah} 1>> ${logfile} 2>> ${errfile}

         echo "The available disk size is not enough to run here." > ${diskfull}
         /bin/rm -f ${dontcopy}
      fi

   elif [ -s ${diskfull} ]
   then
      last=$(date +%s -r ${diskfull})
      now=$(date +%s)
      let howlong=${now}-${last}

      if [ ${howlong} -gt 43200 ]
      then

         echo "Copy files to node ${thismach}." > ${dontcopy}

         blah=" + Files are not all here, cleaning the directory before we copy again."
         blah="${blah}  Please hold!"
         echo ${blah} 1>> ${logfile} 2>> ${errfile}

         if [ -s ${node_thisscen} ]
         then
            blah="  - Delete old stuff from the meterological forcing driver."
            echo ${blah} 1>> ${logfile} 2>> ${errfile}
            /bin/rm -fr ${node_thisscen} 1>> ${logfile} 2>> ${errfile}
         fi

         #----- First thing we do is to check whether this disk is full. ------------------------#
         blah="  - Check whether there is enough disk space."
         echo ${blah} 1>> ${logfile} 2>> ${errfile}
         ans=$(df /scratch)
         nwords=$(echo ${ans} | wc -w)
         let avail=${nwords}-2
         space=$(echo ${ans} | awk '{print $'${avail}'}')

         if [ ${space} -gt ${datasize} ]
         then

            #----- Copy the meteorological forcing. ---------------------------------------#
            blah="  - Copy the meterological forcing driver."
            echo ${blah} 1>> ${logfile} 2>> ${errfile}
            ${rsync} ${source_thisscen} ${node_scenroot} 1>> ${logfile} 2>> ${errfile}


            #----- Copy finished.  Create a file to unlock this node. ---------------------#
            echo "All the data needed are here!" > ${unlocked}

            blah=" + The files have been successfully copied to ${thismach}."
            blah="${blah}  The run will start soon!"
            echo ${blah} 1>> ${logfile} 2>> ${errfile}
            
            success=true
         else
            blah="  - Sorry, but there is still not enough disk space here."
            echo ${blah} 1>> ${logfile} 2>> ${errfile}

            echo "The available disk size is not enough to run here." > ${diskfull}
            /bin/rm -f ${dontcopy}
         fi
      fi

   elif [ ! -s ${unlocked} ]
   then

      pp=-1
      while [ ! -s ${unlocked} ] && [ ! -s ${diskfull} ]
      do
         let pp=${pp}+1

         #----- Print some messages to entertain the bored user. -----------------------------#
         if [ ${pp} -lt 2 ]
         then
            blah=" - Input data are not ready yet. (${pp} minute waiting so far)"
         else
            blah=" - Input data are not ready yet. (${pp} minutes waiting so far)"
         fi

         echo ${blah} 1>> ${logfile} 2>> ${errfile}

         sleep 60
      done

      if [ -s ${unlocked} ]
      then
         blah=" + Thank you for waiting.  The files are here and the run will start soon!"
         echo ${blah} 1>> ${logfile} 2>> ${errfile}
         success=true
      else
         blah=" + Bad news, there wasn't enough disk space in this node."
         echo ${blah} 1>> ${logfile} 2>> ${errfile}
      fi
      #------------------------------------------------------------------------------------#
   else
      success=true
   fi
   #---------------------------------------------------------------------------------------#
else
   success=true
fi
#------------------------------------------------------------------------------------------#



#----- Start the run. ---------------------------------------------------------------------#
if ${success}
then
   blah="Wait ${zzz} seconds before starting the run."
   echo ${blah} 1>> ${logfile} 2>> ${errfile}
   sleep ${zzz}

   #---- Get number of threads. -----------------------------------------------------------#
   if [ "${SLURM_CPUS_PER_TASK}" == "" ] && [ "${OMP_NUM_THREADS}" == "" ]
   then
      export OMP_NUM_THREADS=${n_cpt}
   elif [ "${SLURM_CPUS_PER_TASK}" != "" ]
   then
      export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
   fi
   #---------------------------------------------------------------------------------------#



   #----- Start the execution. ------------------------------------------------------------#
   cd ${here}
   blah=" Run ED-2 using ${OMP_NUM_THREADS} threads."
   echo ${blah} 1>> ${logfile} 2>> ${errfile}
   comm="${exec} -f ${here}/ED2IN"
   blah=" Command to be called: \"${comm}\""
   echo ${blah} 1>> ${logfile} 2>> ${errfile}
   ${comm} 1>> ${logfile} 2>> ${errfile}
   #---------------------------------------------------------------------------------------#


   blah=" Model simulation has ended."
   echo ${blah} 1>> ${logfile} 2>> ${errfile}
   cd ${currloc}
else
   blah="Execution will not happen this time. Better luck next time!"
   echo ${blah} 1>> ${logfile} 2>> ${errfile}
fi
#------------------------------------------------------------------------------------------#
