#!/bin/bash
#==========================================================================================#
#==========================================================================================#
#  Simple script to check simulation status.                                               #
#------------------------------------------------------------------------------------------#
here=$(pwd)

#----- Find the number of arguments. ------------------------------------------------------#
nargs=$#
args=$@
#------------------------------------------------------------------------------------------#


#----- Initialise variables. --------------------------------------------------------------#
SUITE=""
PLATFORM=""
OUTFORM="JobName%200,State%12"
VERBOSE=false
#------------------------------------------------------------------------------------------#



#----- Argument parsing. ------------------------------------------------------------------#
while [[ ${#} > 0 ]]
do
   key="${1}"
   case ${key} in
   -p|--platform)
      PLATFORM="${2}"
      shift 2
      ;;
   -s|--suite)
      SUITE="${2}"
      shift 2
      ;;
   -o|--outform)
      OUTFORM="${2}"
      shift 2
      ;;
   -v|--verbose)
      VERBOSE=true
      shift 1
      ;;
   *)
      echo "Unrecognised option ${key}"
      echo "Usage: ./check_runs.sh -p <platform> -s <suite> [-o <outform>]"
      exit 1
      ;;
   esac
done
#------------------------------------------------------------------------------------------#


#---- Ask user in case they didn't provide the platform. ----------------------------------#
if [ "${PLATFORM}" == "" ]
then
   echo -n  "Which platform are you using? (e.g. odyssey, sdumont)   "
   read PLATFORM
fi
#------------------------------------------------------------------------------------------#



#---- Ask the user in case they didn't provide the test suite path. -----------------------#
if [ "${SUITE}" == "" ]
then
   echo ""
   echo ""
   /bin/ls ${here}
   echo ""
   echo -n  "Which test suite to check? (possibly listed above)   "
   read SUITE
fi
#------------------------------------------------------------------------------------------#

#----- Delete trailing slashes from SUITE. ------------------------------------------------#
SUITE=$(echo ${SUITE} | sed s@"/$"@""@g)
#------------------------------------------------------------------------------------------#

#----- List sites. ------------------------------------------------------------------------#
sites=$(ls ${SUITE}/ED2IN-???-MAIN | sed s@"${SUITE}/"@""@g | cut -c 7-9)
sites=$(echo ${sites} | tr '[:upper:]' '[:lower:]')
nsites=$(ls ${SUITE}/ED2IN-???-MAIN | wc -l)
#------------------------------------------------------------------------------------------#

#---- List tests. -------------------------------------------------------------------------#
sims="dbug test main"
nsims=$(echo ${sims} | wc -w)
#------------------------------------------------------------------------------------------#

let nsuite=${nsites}*${nsims}

#----- Go through all sites. --------------------------------------------------------------#
ffout=
for site in ${sites}
do

   echo " Site ${site}:"

   for sim in ${sims}
   do

      #------------------------------------------------------------------------------------#
      #    Format count.                                                                   #
      #------------------------------------------------------------------------------------#
      if   [ ${nsuite} -ge 10   ] && [ ${nsuite} -lt 100   ]
      then
         ffout=$(printf '%2.2i' ${ff})
      elif [ ${nsuite} -ge 100  ] && [ ${nsuite} -lt 1000  ]
      then
         ffout=$(printf '%2.2i' ${ff})
      elif [ ${nsuite} -ge 100  ] && [ ${nsuite} -lt 10000 ]
      then
         ffout=$(printf '%2.2i' ${ff})
      else
         ffout=${ff}
      fi
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Set some variables to check whether the simulation is running.                 #
      #------------------------------------------------------------------------------------#
      stdout="${here}/${SUITE}/${sim}_${site}.out"
      stderr="${here}/${SUITE}/${sim}_${site}.err"
      jobname="${SUITE}_${site}_${sim}"
      #------------------------------------------------------------------------------------#


      if [ -s ${stdout} ]
      then
         case ${PLATFORM} in
         odyssey|sdumont)
            stask="stask --noheader -u ${USER} -j ${jobname} -t ${jobname}"
            running=$(${stask} -o "${OUTFORM}" | grep "RUNNING"   | wc -l)
            pending=$(${stask} -o "${OUTFORM}" | grep "PENDING"   | wc -l)
            suspend=$(${stask} -o "${OUTFORM}" | grep "SUSPENDED" | wc -l)
            ;;
         sun-lncc)
            running=$(qcheck   -s R -n | grep ${jobname} 2> /dev/null | wc -l)
            pending=$(qcheck   -s P -n | grep ${jobname} 2> /dev/null | wc -l)
            suspend=$(qcheck -s S -n | grep ${jobname} 2> /dev/null | wc -l)
            ;;
         # lbl)
         #    DCHECK=$(qstat -u ${USER} | grep -is $tag | grep -is dbug | wc -l)
         #    if [ $DCHECK -gt 0 ]; then
         #        DRES="RUNN"
         #    else
         #        lastfile=`ls ${folder}/dbug_${tag}*out | awk '{ f=$NF };END{ print f }'`
         #        if [ -f $lastfile ]; then
         #            DCHECK=$(grep -is 'ED execution ends' $lastfile | wc -l)
         #            if [ $DCHECK -gt 0 ];then
         #                DRES="COMP"
         #            else
         #                DRES="FAIL"
         #            fi
         #        else
         #            DRES="----"
         #        fi
         #    fi
         #    ;;
         *)
            echo "Platform is not recognised.  Check script and add your commands."
            exit 1
            ;;
         esac
         #---------------------------------------------------------------------------------#


         #---- Check the current simulation time. -----------------------------------------#
         simline=$(grep "Simulating: "   ${stdout} | tail -1)
         runtime=$(echo ${simline} | awk '{print $3}')
         #---------------------------------------------------------------------------------#


         #----- Check for segmentation violations. ----------------------------------------#
         if [ -s ${stderr} ]
         then
            segv1=$(grep -i "sigsegv"            ${stderr} | wc -l)
            segv2=$(grep -i "segmentation fault" ${stderr} | wc -l)
            let sigsegv=${segv1}+${segv2}
         else
            sigsegv=0
         fi
         #---------------------------------------------------------------------------------#



         #----- Check whether met files are missing... (bad start) ------------------------#
         metbs1=$(grep "Cannot open met driver input file" ${stdout} | wc -l)
         metbs2=$(grep "Specify ED_MET_DRIVER_DB properly" ${stdout} | wc -l)
         let metmiss=${metbs1}+${metbs2}
         #---------------------------------------------------------------------------------#



         #----- Check whether met files are missing... (bad start) ------------------------#
         metbs1=$(grep "Cannot open met driver input file" ${stdout} | wc -l)
         metbs2=$(grep "Specify ED_MET_DRIVER_DB properly" ${stdout} | wc -l)
         let metmiss=${metbs1}+${metbs2}
         #---------------------------------------------------------------------------------#



         #----- Check for other possible outcomes. ----------------------------------------#
         stopped=$(grep "FATAL ERROR"                       ${stdout} | wc -l)
         crashed=$(grep "IFLAG1 problem."                   ${stdout} | wc -l)
         bad_met=$(grep "Meteorological forcing has issues" ${stdout} | wc -l)
         the_end=$(grep "ED-2.2 execution ends"             ${stdout} | wc -l)
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Plot a message so the user knows what is going on.                          #
         #---------------------------------------------------------------------------------#
         if [ ${sigsegv} -gt 0 ]
         then
            echo  "    ${sim} HAD SEGMENTATION VIOLATION."
         elif [ ${crashed} -gt 0 ]
         then 
            echo  "    ${sim} HAS CRASHED (RK4 PROBLEM)."
         elif [ ${bad_met} -gt 0 ]
         then 
            echo  "    ${sim} HAS BAD MET DRIVERS."
         elif [ ${metmiss} -gt 0 ]
         then 
            echo  "    ${sim} DID NOT FIND MET DRIVERS."
         elif [ ${stopped} -gt 0 ]
         then
            echo  "    ${sim} STOPPED (UNKNOWN REASON)!"
         elif [ ${the_end} -gt 0 ]
         then
            echo  "    ${sim} has finished."
         elif [ ${suspend} -gt 0 ] 
         then
            echo  "    ${sim} is SUSPENDED (${runtime})."
         elif [ ${running} -gt 0 ] 
         then
            echo  "    ${sim} is running (${runtime})."
         elif ${VERBOSE}
         then
            echo  "    ${sim} status is UNKNOWN (Last time ${runtime})."
            echo  "       STASK   -- ${stask}"
            echo  "       OUTFORM -- ${OUTFORM}"
            echo  "       RUNNING -- ${running}"
            echo  "       PENDING -- ${pending}"
            echo  "       SUSPEND -- ${suspend}"
            echo  "       SIGSEGV -- ${sigsegv}"
            echo  "       METMISS -- ${metmiss}"
            echo  "       STOPPED -- ${stopped}"
            echo  "       CRASHED -- ${crashed}"
            echo  "       THE_END -- ${the_end}"
         else
            echo  "    ${sim} status is UNKNOWN (Last time ${runtime})."
         fi
         #---------------------------------------------------------------------------------#
   else
      echo "    ${sim} is pending ."
   fi
   #---------------------------------------------------------------------------------------#
   done
   #---------------------------------------------------------------------------------------#
done
#------------------------------------------------------------------------------------------#
