#!/bin/bash
. ${HOME}/.bashrc


#------ Simulation settings. --------------------------------------------------------------#
main=""
running_msg="is running"
the_end_msg="has finished"
unknown_msg="status is unknown"
sigsegv_msg="HAD SEGMENTATION VIOLATION"
crashed_msg="RK4 PROBLEM"
metmiss_msg="DID NOT FIND MET DRIVERS"
stopped_msg="UNKNOWN REASON"
fstline_msg="Number of polygons:"
wait_minutes=30
frqpost=4
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#       Define some unique and temporary output files.                                     #
#------------------------------------------------------------------------------------------#
check_out="/tmp/check_run_${$}.out"
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#       First check that the main path has been set.  If now, don't even start.            #
#------------------------------------------------------------------------------------------#
if [ "x${main}" == "x" ]
then
   echo " + You must set variable \"main\" in your script before running!"
   exit 99
fi
#------------------------------------------------------------------------------------------#


#------ Move to the current directory. ----------------------------------------------------#
cd ${main}
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#       Main loop: the script stays inside until all simulations have finished.            #
#------------------------------------------------------------------------------------------#
n_ongoing=9999
it=0
while [ ${n_ongoing} -gt 0 ]
do
   #------ Update iteration counter. ------------------------------------------------------#
   let it=${it}+1
   #---------------------------------------------------------------------------------------#


   #------ Run check_run.sh. --------------------------------------------------------------#
   ./check_run.sh 1> ${check_out} 2>&1
   #---------------------------------------------------------------------------------------#

   #------ Find the number of simulations with different statuses. ------------------------#
   header=$(cat ${check_out} | grep "${fstline_msg}")
   n_polygon=$(echo "${header}" | awk '{print $4}' | sed s@"\.\.\."@@g)
   n_running=$(cat ${check_out} | grep "${running_msg}" | wc -l)
   n_the_end=$(cat ${check_out} | grep "${the_end_msg}" | wc -l)
   n_unknown=$(cat ${check_out} | grep "${unknown_msg}" | wc -l)
   n_sigsegv=$(cat ${check_out} | grep "${sigsegv_msg}" | wc -l)
   n_crashed=$(cat ${check_out} | grep "${crashed_msg}" | wc -l)
   n_metmiss=$(cat ${check_out} | grep "${metmiss_msg}" | wc -l)
   n_stopped=$(cat ${check_out} | grep "${stopped_msg}" | wc -l)
   let n_ongoing=${n_polygon}-${n_the_end}
   let n_problem=${n_unknown}+${n_sigsegv}+${n_crashed}+${n_metmiss}+${n_stopped}
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #       Print report.                                                                   #
   #---------------------------------------------------------------------------------------#
   echo "==================================================================="
   echo "==================================================================="
   echo " Iteration: ${it}."
   echo " Current time: $(date +'%a %d-%b-%Y %H:%M:%S %Z')"
   echo " Total number of simulations: ${n_polygon}"
   echo " "
   echo " Individual polygon check: "
   cat ${check_out} | tail -${n_polygon}
   echo " "
   echo " Status count:"
   echo " Running:                       ${n_running}"
   echo " Finished:                      ${n_the_end}"
   echo " Unknown:                       ${n_unknown}"
   echo " Segmentation violation:        ${n_sigsegv}"
   echo " Crashed (RK4 problem):         ${n_crashed}"
   echo " Missing meteorological driver: ${n_metmiss}"
   echo " Stopped for some other reason: ${n_stopped}"
   echo "-------------------------------------------------------------------"
   echo ""
   #---------------------------------------------------------------------------------------#

   rm -f ${check_out}
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      In case any simulation is flagged as 'problemati
   #---------------------------------------------------------------------------------------#
   if [ ${n_problem} -gt 0 ]
   then
      echo "    + Some polygons had problems, re-submit simulations."
      ./spawn_poly.sh
   fi
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Run the post processing.                                                         #
   #---------------------------------------------------------------------------------------#
   let xpost=${it}%${frqpost}
   if [ -s "${main}/epost.sh" ] && [ ${xpost} -eq 0 -o ${n_ongoing} -eq 0 ]
   then
      echo "    + Run post-processing."
      ./epost.sh
   fi
   #---------------------------------------------------------------------------------------#

   #---------------------------------------------------------------------------------------#
   #       In case any simulation is still on the works, take a break before checking      #
   #   again.                                                                              #
   #---------------------------------------------------------------------------------------#
   if [ ${n_ongoing} -gt 0 ]
   then
      echo "     + Take a ${wait_minutes}-minute break before checking again."

      echo ""
      echo "==================================================================="
      echo "==================================================================="
      echo ""
      sleep "${wait_minutes}m"
   else
      echo "     + All simulations have finished."
      echo ""
      echo "==================================================================="
      echo "==================================================================="
      echo ""
   fi
   #---------------------------------------------------------------------------------------#
done
#------------------------------------------------------------------------------------------#



#------ Compress files to save disk space. ------------------------------------------------#
echo " + Compress files."
./last_histo.sh
#------------------------------------------------------------------------------------------#
