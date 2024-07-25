#!/bin/bash
. ~/.bashrc
#----- Main path, usually set by $(pwd) so you don't need to change it. -------------------#
here=""
#----- Computer name, to decide which queues are available. -------------------------------#
ordinateur=$(hostname -s)
#----- User name, usually set by $(whoami) so you don't need to change it. ----------------#
myself=$(whoami)
#----- Description of this simulation, used to create unique job names. -------------------#
desc=$(basename ${here})
#----- File containing the list of jobs and their settings: -------------------------------#
joborder="${here}/joborder.txt"         # ! File with the job instructions
#----- POV-Ray include files. -------------------------------------------------------------#
pov_incs="${HOME}/Util/Modules/povray/3.6/share/povray-3.6/include"
#----- How should the post-processing be handled? -----------------------------------------#
submit=false     # true -- Try to submit the script to the queue
                 # false -- Prepare script for batch, but don't dispatch job.
#----- Time frame for post-processing. ----------------------------------------------------#
useperiod="a"    # Which bounds should I use? (Ignored by plot_eval_ed.r)
                 # "a" -- All period
                 # "t" -- One eddy flux tower met cycle
                 # "u" -- User defined period, defined by the variables below.
                 # "f" -- Force the tower cycle.  You may need to edit the script, though
                 # "b" -- Force one biometry cycle.
yusera=1972      # First year to use
yuserz=2011      # Last year to use
#----- Yearly comparison . ----------------------------------------------------------------#
seasonmona=1
#----- Census comparison. -----------------------------------------------------------------#
varcycle="TRUE"  # Find the average mortality for various cycles (TRUE/FALSE).
#----- Hourly comparison. -----------------------------------------------------------------#
usedistrib="edf" # Which distribution to plot on top of histograms:
                 #   norm -- Normal distribution
                 #   sn   -- Skewed normal distribution      (requires package sn)
                 #   edf  -- Empirical distribution function (function density)
#----- Output format. ---------------------------------------------------------------------#
outform="c(\"pdf\")"           # x11 - On screen (deprecated on shell scripts)
                               # png - Portable Network Graphics
                               # tif - 
                               # eps - Encapsulated Post Script
                               # pdf - Portable Document Format
#----- DBH classes. -----------------------------------------------------------------------#
idbhtype=4                     # Type of DBH class
                               # 1 -- Every 10 cm until 100cm; > 100cm
                               # 2 -- 0-10; 10-20; 20-35; 35-55; 55-80; > 80 (cm)
                               # 3 -- 0-10; 10-35; 35-55; > 55 (cm)
                               # 4 -- 0-10; 10-30; 30-50; 50-80; > 80 (cm)
#----- Default background colour. ---------------------------------------------------------#
background=0                   # 0 -- White
                               # 1 -- Pitch black
                               # 2 -- Dark grey
#----- Select integration interval for some photosynthesis-related variables. -------------#
iint_photo=1                   # 0 -- 24h
                               # 1 -- daytime only
#----- Trim the year comparison for tower years only? -------------------------------------#
efttrim="TRUE"
#----- Correction factor for respiration. -------------------------------------------------#
correct_gs=1.0                 # Correction factor for growth and storage respiration
#----- Use only old-growth patches for census comparison? (plot_census.r only). -----------#
oldgrowth="FALSE"
#----- Path with R scripts that are useful. -----------------------------------------------#
rscpath="${HOME}/EDBRAMS/R-utils"
#----- bashrc (usually ${HOME}/.bashrc). --------------------------------------------------#
initrc="${HOME}/.bashrc"
#----- Initialisation scripts. ------------------------------------------------------------#
optsrc="-n"                   # Option for .bashrc (for special submission settings)
                              #   In case none is needed, leave it blank ("").
#----- Settings for this group of polygons. -----------------------------------------------#
global_queue=""               # Queue
reservation=""                # Reservation
overcommit=""                 # Ignore maximum task limit for partition? (true/false)
partial=false                 # Partial submission (false will ignore polya and npartial
                              #    and send all polygons.
skip_end=true                 # Skip processing in case the R script has already loaded
                              #    all files through the end of the simulation
                              #    (true/false).  This is only used for "monthly"-based
                              #    scripts (those marked with (*) in the table below).
polya=501                     # First polygon to submit
npartial=100                  # Maximum number of polygons to include in this bundle
                              #    (actual number will be adjusted for total number of 
                              #     polygons if needed be).
dttask=2                      # Time to wait between task submission
sim_memory=0                  # Memory per simulation.  If zero, then it will be 
                              #    automatically determined by the maximum number of tasks
                              #    per node.
runtime="00:00:00"            # Simulation time in hours.  If zero, then it will be
                              #    automatically determined by the maximum walltime for 
                              #    the queue.
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#     Which scripts to run.                                                                #
#                                                                                          #
#   - read_monthly.r  - (*) This reads the monthly mean files (results can then be used    #
#                       for plot_monthly.r, plot_yearly.r, and others, but it doesn't plot #
#                       anything.)                                                         #
#   - yearly_ascii.r  - (*) This creates three ascii (csv) files with annual averages of   #
#                       various variables.  It doesn't have all possible variables as it   #
#                       is intended to simplify the output for learning purposes.          #
#   - monthly_ascii.r - (*) This creates three ascii (csv) files with annual averages of   #
#                       various variables.  It doesn't have all possible variables as it   #
#                       is intended to simplify the output for learning purposes.          #
#   - plot_monthly.r  - (*) This creates several plots based on the monthly mean output.   #
#   - plot_yearly.r   - (*) This creates plots with year time series.                      #
#   - plot_ycomp.r    - (*) This creates yearly comparisons based on the monthly mean      #
#                       output.                                                            #
#   - plot_povray.r   - (*) This creates yearly plots of the polygon using POV-Ray.        #
#   - plot_rk4.r      - This creates plots from the detailed output for Runge-Kutta.       #
#                       (patch-level only).                                                #
#   - plot_photo.r    - This creates plots from the detailed output for Farquhar-Leuning.  #
#   - plot_rk4pc.r    - This creates plots from the detailed output for Runge-Kutta.       #
#                       (patch- and cohort-level).                                         #
#   - plot_budget.r   - This creates plots from the detailed budget for Runge-Kutta.       #
#                       (patch-level only).                                                #
#   - plot_eval_ed.r  - This creates plots comparing model with eddy flux observations.    #
#   - plot_census.r   - This creates plots comparing model with biometric data.            #
#   - whichrun.r      - This checks the run status.                                        #
#                                                                                          #
#   The following scripts should work too, but I haven't tested them.                      #
#   - plot_daily.r    - This creates plots from the daily mean output.                     #
#   - plot_fast.r     - This creates plots from the analysis files.                        #
#   - patchprops.r    - This creates simple plots showing the patch structure.             #
#   - reject_ed.r     - This tracks the number of steps that were rejected, and what       #
#                       caused the step to be rejected.                                    #
#------------------------------------------------------------------------------------------#
rscript=""
#rscript="yearly_ascii.r"
#rscript="plot_monthly.r"
#rscript="plot_census.r" 
#rscript="plot_ycomp.r"
#rscript="plot_eval_ed.r"
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#    Tell whether to plot pseudo-drought or not.                                           #
#------------------------------------------------------------------------------------------#
droughtmark="FALSE"         # Should I plot a rectangle to show the drought?
                            #     capital letters only: TRUE means yes, FALSE means no
droughtyeara=1605           # Year that the first drought instance happens (even if it is 
                            #     just the last bit)
droughtyearz=1609           # Year that the last drought instance happens (even if it 
                            #     partial)
monthsdrought="c(12,1,2,3)" # List of months that get drought, if it starts late in the
                            #     year, put the last month first.
#------------------------------------------------------------------------------------------#


#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#     Unless you are going to modify the way jobs are submitted, you don't need to change  #
# anything beyond this point.                                                              #
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#


#------------------------------------------------------------------------------------------#
#       First check that the main path and e-mail have been set.  If not, don't run.       #
#------------------------------------------------------------------------------------------#
if [[ "x${here}" == "x" ]] || [[ "x${global_queue}" == "x" ]] || [[ "x${rscript}" == "x" ]]
then
   echo " You must set some variables before running the script:"
   echo " Check variables \"here\", \"global_queue\" and \"rscript\"!"
   exit 99
fi
#------------------------------------------------------------------------------------------#


#----- Load settings. ---------------------------------------------------------------------#
if [[ -s ${initrc} ]]
then
   . ${initrc}
fi
#------------------------------------------------------------------------------------------#


#----- Find out which platform we are using. ----------------------------------------------#
host=$(hostname -s)
case ${host} in
rclogin*|holy*|moorcroft*|rcnx*)
   cluster="CANNON"
   ;;
sdumont*)
   cluster="SDUMONT"
   ;;
*)
   echo -n "Failed guessing cluster from node name.  Please type the name:   "
   read cluster
   ;;
esac
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Configurations depend on the global_queue and cluster.                               #
#------------------------------------------------------------------------------------------#
case "${cluster}" in
SDUMONT)
   #----- Set parameters according to partitions. -----------------------------------------#
   case ${global_queue} in
   cpu_long|nvidia_long)
      n_nodes_max=10
      n_cpt=1
      n_tpn=24
      runtime_max="31-00:00:00"
      node_memory=64000
      ;;
   cpu|nvidia|phi)
      n_nodes_max=50
      n_cpt=1
      n_tpn=24
      runtime_max="2-00:00:00"
      node_memory=64000
      ;;
   cpu_dev)
      n_nodes_max=20
      n_cpt=1
      n_tpn=24
      runtime_max="02:00:00"
      node_memory=64000
      ;;
   nvidia_dev|phi_dev)
      n_nodes_max=2
      n_cpt=1
      n_tpn=24
      runtime_max="02:00:00"
      node_memory=64000
      ;;
   cpu_scal|nvidia_scal)
      n_nodes_max=128
      n_cpt=1
      n_tpn=24
      runtime_max="18:00:00"
      node_memory=64000
      ;;
   *)
      echo "Global queue ${global_queue} is not recognised!"
      exit
      ;;
   esac
   #---------------------------------------------------------------------------------------#
   ;;
CANNON)
   #----- Set parameters according to partitions. -----------------------------------------#
   case ${global_queue} in
   "commons")
      n_nodes_max=5
      n_cpt=1
      n_tpn=32
      runtime_max="14-00:00:00"
      node_memory=126820
      ;;
   "huce_cascade")
      n_nodes_max=30
      n_cpt=1
      n_tpn=48
      runtime_max="14-00:00:00"
      node_memory=183240
      ;;
   "huce_intel")
      n_nodes_max=150
      n_cpt=1
      n_tpn=24
      runtime_max="14-00:00:00"
      node_memory=122090
      ;;
   "shared,huce_intel"|"huce_intel,shared")
      n_nodes_max=150
      n_cpt=1
      n_tpn=24
      runtime_max="7-00:00:00"
      node_memory=122090
      ;;
   "shared")
      n_nodes_max=220
      n_cpt=1
      n_tpn=48
      runtime_max="7-00:00:00"
      node_memory=183240
      ;;
   "test")
      n_nodes_max=9
      n_cpt=1
      n_tpn=48
      runtime_max="8:00:00"
      node_memory=183240
      ;;
   "unrestricted")
      n_nodes_max=4
      n_cpt=1
      n_tpn=48
      runtime_max="infinite"
      node_memory=183240
      ;;
   *)
      echo "Global queue ${global_queue} is not recognised!"
      exit
      ;;
   esac
   #---------------------------------------------------------------------------------------#
   ;;
esac
#------------------------------------------------------------------------------------------#


#----- Set the number of tasks. -----------------------------------------------------------#
let n_tasks_max=${n_nodes_max}*${n_tpn}
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#    Set time.                                                                             #
#------------------------------------------------------------------------------------------#
runtime=$(echo     ${runtime}     | tr '[:upper:]' '[:lower:]')
runtime_max=$(echo ${runtime_max} | tr '[:upper:]' '[:lower:]')
case "${runtime}" in
infinite)
   #----- Infinite runtime.  Make sure the queue supports this type of submission. --------#
   case "${runtime_max}" in
   infinite)
      echo "" > /dev/null
      ;;
   *)
      echo " Requested partition:       ${global_queue}"
      echo " Maximum runtime permitted: ${runtime_max}"
      echo " Requested runtime:         ${runtime}"
      echo " Partition ${global_queue} does not support infinite time."
      exit 91
      ;;
   esac
   #---------------------------------------------------------------------------------------#
   ;;
*)
   #----- Find out the format provided. ---------------------------------------------------#
   case "${runtime}" in
   *-*:*)
      #----- dd-hh:mm:ss. -----------------------------------------------------------------#
      ndays=$(echo ${runtime} | sed s@"-.*"@@g)
      nhours=$(echo ${runtime} | sed s@"^.*-"@@g | sed s@":.*"@@g)
      #------------------------------------------------------------------------------------#
      ;;
   *:*)
      #----- hh:mm:ss. --------------------------------------------------------------------#
      ndays=0
      nhours=$(echo ${runtime} | sed s@":.*"@@g)
      #------------------------------------------------------------------------------------#
      ;;
   *)
      #----- Hours. -----------------------------------------------------------------------#
      let ndays="10#${runtime}"/24
      let nhours="10#${runtime}"%24
      #------------------------------------------------------------------------------------#
      ;;
   esac
   #---------------------------------------------------------------------------------------#


   #----- Find the walltime in hours, and the runtime in nice format. ---------------------#
   let wall="10#${nhours}"+24*"10#${ndays}"
   let ndays="10#${wall}"/24
   let nhours="10#${wall}"%24
   if [[ ${ndays} -gt 0 ]]
   then
      fmtday=$(printf '%2.2i' ${ndays})
      fmthr=$(printf '%2.2i' ${nhours})
      runtime="${fmtday}-${fmthr}:00:00"
   else
      fmthr=$(printf '%2.2i' ${nhours})
      runtime="${fmthr}:00:00"
   fi
   #---------------------------------------------------------------------------------------#



   #----- Find the maximum number of hours allowed in the partition. ----------------------#
   case "${runtime_max}" in
   infinite)
      let ndays_max="10#${ndays}"+1
      let nhours_max="10#${nhours}"
      ;;
   *-*:*)
      #----- dd-hh:mm:ss. -----------------------------------------------------------------#
      ndays_max=$(echo ${runtime_max} | sed s@"-.*"@@g)
      nhours_max=$(echo ${runtime_max} | sed s@"^.*-"@@g | sed s@":.*"@@g)
      #------------------------------------------------------------------------------------#
      ;;
   *:*)
      #----- hh:mm:ss. --------------------------------------------------------------------#
      ndays_max=0
      nhours_max=$(echo ${runtime_max} | sed s@":.*"@@g)
      #------------------------------------------------------------------------------------#
      ;;
   *)
      ndays_max=0
      nhours_max=$(echo ${runtime_max} | sed s@":.*"@@g)
      ;;
   esac
   let wall_max="10#${nhours_max}"+24*"10#${ndays_max}"
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Check requested walltime and the availability.                                    #
   #---------------------------------------------------------------------------------------#
   if [[ ${wall} -eq 0 ]]
   then
      case "${runtime_max}" in
      infinite) runtime="infinite"     ;;
      *)        runtime=${runtime_max} ;;
      esac
   elif [[ ${wall} -gt ${wall_max} ]]
   then
      echo " Requested partition:       ${global_queue}"
      echo " Maximum runtime permitted: ${runtime_max}"
      echo " Requested runtime:         ${runtime}"
      echo " - Requested time exceeds limits."
      exit 92
   fi
   #---------------------------------------------------------------------------------------#
   ;;
esac
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#    Use the general path.                                                                 #
#------------------------------------------------------------------------------------------#
case ${myself} in
mlongo|marcosl)
   rscpath="${HOME}/Util/Rsc"
   ;;
marcos.longo)
   rscpath="${SCRATCH}/Util/Rsc"
   rsync -Prutv ${R_SCRP}/* ${rscpath}
   ;;
esac
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#      Define the job name, and the names of the output files.                             #
#------------------------------------------------------------------------------------------#
case ${rscript} in
read_monthly.r)
   epostkey="rmon"
   ;;
yearly_ascii.r)
   epostkey="yasc"
   ;;
monthly_ascii.r)
   epostkey="masc"
   ;;
r10_monthly.r)
   epostkey="rm10"
   ;;
plot_monthly.r)
   epostkey="pmon"
   ;;
plot_yearly.r)
   epostkey="pyrs"
   ;;
plot_ycomp.r)
   epostkey="pycp"
   ;;
plot_census.r)
   epostkey="pcen"
   ;;
plot_povray.r)
   epostkey="ppov"
   ;;
plot_eval_ed.r)
   epostkey="peed"
   ;;
plot_budget.r)
   epostkey="pbdg"
   ;;
plot_rk4.r)
   epostkey="prk4"
   ;;
plot_rk4pc.r)
   epostkey="prpc"
   ;;
plot_photo.r)
   epostkey="ppht"
   ;;
reject_ed.r)
   epostkey="prej"
   ;;
patchprops.r)
  epostkey="ppro"
  ;;
whichrun.r)
  epostkey="pwhr"
  ;;
plot_daily.r)
   epostkey="pday"
   ;;
plot_fast.r)
   epostkey="pfst"
   ;;
*)
   #---------------------------------------------------------------------------------------#
   #     If the script is here, then it could not find the script... And this should never #
   # happen, so interrupt the script!                                                      #
   #---------------------------------------------------------------------------------------#
   echo " Script ${rscript} is not recognised by epost.sh!"
   exit 1
   #---------------------------------------------------------------------------------------#
   ;;
esac
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#   Make sure memory does not exceed maximum amount that can be requested.                 #
#------------------------------------------------------------------------------------------#
if [[ ${sim_memory} -eq 0 ]]
then
   let sim_memory=${node_memory}/${n_tpn}
   let node_memory=${n_tpn}*${sim_memory}
elif [[ ${sim_memory} -gt ${node_memory} ]]
then 
   echo "Simulation memory ${sim_memory} cannot exceed node memory ${node_memory}!"
   exit 99
else
   #------ Set memory and number of CPUs per task. ----------------------------------------#
   let n_tpn_try=${node_memory}/${sim_memory}
   if [[ ${n_tpn_try} -le ${n_tpn} ]]
   then
      n_tpn=${n_tpn_try}
      let sim_memory=${node_memory}/${n_tpn}
   else
      let node_memory=${n_tpn}*${sim_memory}
   fi
   #---------------------------------------------------------------------------------------#
fi
#------------------------------------------------------------------------------------------#




#----- Determine the number of polygons to run. -------------------------------------------#
let npolys=$(wc -l ${joborder} | awk '{print $1 }')-3
if [[ ${npolys} -lt 100 ]]
then
   ndig=2
elif [[ ${npolys} -lt 1000 ]]
then
   ndig=3
elif [[ ${npolys} -lt 10000 ]]
then
   ndig=4
else
   ndig=5
fi
#------------------------------------------------------------------------------------------#



#---- Partial or complete. ----------------------------------------------------------------#
rprefix=$(basename ${rscript} .r)
if ${partial}
then
   let ff=${polya}-1
   let polyz=${ff}+${npartial}
   if [[ ${polyz} -gt ${npolys} ]]
   then
      polyz=${npolys}
   fi
   partlabel="$(printf '%3.3i' ${polya})-$(printf '%3.3i' ${polyz})"
   sbatch="${here}/sub_${rprefix}_${partlabel}.sh"
   obatch="${here}/out_${rprefix}_${partlabel}.log"
   ebatch="${here}/err_${rprefix}_${partlabel}.log"
   epostjob="${epostkey}-${desc}_${partlabel}"
else
   ff=0
   polya=1
   polyz=${npolys}
   sbatch="${here}/sub_${rprefix}.sh"
   obatch="${here}/out_${rprefix}.log"
   ebatch="${here}/err_${rprefix}.log"
   epostjob="${epostkey}-${desc}"
fi
let ntasks=1+${polyz}-${polya}
#------------------------------------------------------------------------------------------#




#----- Summary for this submission preparation.  Then give 5 seconds for user to cancel. --#
echo "------------------------------------------------"
echo "  Submission summary: "
echo ""
echo "  Memory per cpu:      ${sim_memory}"
echo "  Tasks per node:      ${n_tpn}"
echo "  Partition:           ${global_queue}"
echo "  Run time:            ${runtime}"
echo "  First polygon:       ${polya}"
echo "  Last polygon:        ${polyz}"
echo "  Total polygon count: ${npolys}"
echo " "
echo " Partial submission:   ${partial}"
echo " Automatic submission: ${submit}"
echo " "
echo " R script:             ${rscript}"
echo " Submission script:    $(basename ${sbatch})"
echo "------------------------------------------------"
echo ""
echo ""
echo -n " Waiting five seconds before proceeding... "
sleep 5
echo "Done!"
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#    Initialise executable.  The main executable will have different settings depending    #
# on whether this is to be a multi-task single job, or multiple single-task jobs.          #
#------------------------------------------------------------------------------------------#
case "${cluster}" in
SDUMONT)
   #----- Initialise script to submit a single multi-task job. ----------------------------#
   rm -fr ${sbatch}
   touch ${sbatch}
   chmod u+x ${sbatch}
   echo "#!/bin/bash" >> ${sbatch}
   echo "#SBATCH --ntasks=${ntasks}              # Number of tasks"            >> ${sbatch}
   echo "#SBATCH --cpus-per-task=1               # Number of CPUs per task"    >> ${sbatch}
   echo "#SBATCH --partition=${global_queue}     # Queue that will run job"    >> ${sbatch}
   if [[ "${reservation}" != "" ]]
   then
      echo "#SBATCH --reservation=${reservation}    # Reserved nodes"          >> ${sbatch}
   fi
   echo "#SBATCH --job-name=${epostjob}          # Job name"                   >> ${sbatch}
   echo "#SBATCH --mem-per-cpu=${sim_memory}     # Memory per CPU"             >> ${sbatch}
   echo "#SBATCH --time=${runtime}               # Time for job"               >> ${sbatch}
   echo "#SBATCH --output=${obatch}              # Standard output path"       >> ${sbatch}
   echo "#SBATCH --error=${ebatch}               # Standard error path"        >> ${sbatch}
   echo ""                                                                     >> ${sbatch}
   echo "#--- Get plenty of memory."                                           >> ${sbatch}
   echo "ulimit -s unlimited"                                                  >> ${sbatch}
   echo ""                                                                     >> ${sbatch}
   echo "#--- Initial settings."                                               >> ${sbatch}
   echo "here=\"${here}\"                            # Main path"              >> ${sbatch}
   echo "rscript=\"${rscript}\"                      # R Script"               >> ${sbatch}
   echo "rstdout=\"${epostout}\"                     # Standard output"        >> ${sbatch}
   echo ""                                                                     >> ${sbatch}
   echo "#--- Print information about this job."                               >> ${sbatch}
   echo "echo \"\""                                                            >> ${sbatch}
   echo "echo \"\""                                                            >> ${sbatch}
   echo "echo \"----- Summary of current job ------------------------------\"" >> ${sbatch}
   echo "echo \" CPUs per task:   \${SLURM_CPUS_PER_TASK}\""                   >> ${sbatch}
   echo "echo \" Job:             \${SLURM_JOB_NAME} (\${SLURM_JOB_ID})\""     >> ${sbatch}
   echo "echo \" Queue:           \${SLURM_JOB_PARTITION}\""                   >> ${sbatch}
   echo "echo \" Number of nodes: \${SLURM_NNODES}\""                          >> ${sbatch}
   echo "echo \" Number of tasks: \${SLURM_NTASKS}\""                          >> ${sbatch}
   echo "echo \" Memory per CPU:  \${SLURM_MEM_PER_CPU}\""                     >> ${sbatch}
   echo "echo \" Memory per node: \${SLURM_MEM_PER_NODE}\""                    >> ${sbatch}
   echo "echo \" Node list:       \${SLURM_JOB_NODELIST}\""                    >> ${sbatch}
   echo "echo \" Time limit:      \${SLURM_TIMELIMIT}\""                       >> ${sbatch}
   echo "echo \" Std. Output:     \${SLURM_STDOUTMODE}\""                      >> ${sbatch}
   echo "echo \" Std. Error:      \${SLURM_STDERRMODE}\""                      >> ${sbatch}
   echo "echo \"-----------------------------------------------------------\"" >> ${sbatch}
   echo "echo \"\""                                                            >> ${sbatch}
   echo "echo \"\""                                                            >> ${sbatch}
   echo "echo \"\""                                                            >> ${sbatch}
   echo "echo \"\""                                                            >> ${sbatch}
   echo ""                                                                     >> ${sbatch}
   echo ""                                                                     >> ${sbatch}
   echo "#--- Define home in case home is not set"                             >> ${sbatch}
   echo "if [[ \"x\${HOME}\" == \"x\" ]]"                                      >> ${sbatch}
   echo "then"                                                                 >> ${sbatch}
   echo "   export HOME=\$(echo ~)"                                            >> ${sbatch}
   echo "fi"                                                                   >> ${sbatch}
   echo ""                                                                     >> ${sbatch}
   echo "#--- Load modules and settings."                                      >> ${sbatch}
   echo ". \${HOME}/.bashrc ${optsrc}"                                         >> ${sbatch}
   echo ""                                                                     >> ${sbatch}
   echo "#----- Task list."                                                    >> ${sbatch}
   #---------------------------------------------------------------------------------------#
   ;;
CANNON)
   #----- Initialise script to submit multiple single-task jobs. --------------------------#
   rm -f ${sbatch}
   touch ${sbatch}
   chmod u+x ${sbatch}
   echo "#!/bin/bash"                                                          >> ${sbatch}
   echo ""                                                                     >> ${sbatch}
   echo "#--- Initial settings."                                               >> ${sbatch}
   echo "here=\"${here}\"                            # Main path"              >> ${sbatch}
   echo "exec=\"${exec_full}\"                       # Executable"             >> ${sbatch}
   echo ""                                                                     >> ${sbatch}
   echo "echo \"----- Global settings for this array of simulations -------\"" >> ${sbatch}
   echo "echo \" Main path:       \${here}\""                                  >> ${sbatch}
   echo "echo \" Executable:      \${exec}\""                                  >> ${sbatch}
   echo "echo \"-----------------------------------------------------------\"" >> ${sbatch}
   echo "echo \"\""                                                            >> ${sbatch}
   echo "echo \"\""                                                            >> ${sbatch}
   echo ""                                                                     >> ${sbatch}
   echo ""                                                                     >> ${sbatch}
   echo "#--- Define home in case home is not set"                             >> ${sbatch}
   echo "if [[ \"x\${HOME}\" == \"x\" ]]"                                      >> ${sbatch}
   echo "then"                                                                 >> ${sbatch}
   echo "   export HOME=\$(echo ~)"                                            >> ${sbatch}
   echo "fi"                                                                   >> ${sbatch}
   echo ""                                                                     >> ${sbatch}
   echo "#--- Load modules and settings."                                      >> ${sbatch}
   echo ". \${HOME}/.bashrc ${optsrc}"                                         >> ${sbatch}
   echo ""                                                                     >> ${sbatch}
   echo "#--- Get plenty of memory."                                           >> ${sbatch}
   echo "ulimit -s unlimited"                                                  >> ${sbatch}
   echo "ulimit -u unlimited"                                                  >> ${sbatch}
   echo ""                                                                     >> ${sbatch}
   echo "#----- Task list."                                                    >> ${sbatch}
   #---------------------------------------------------------------------------------------#
   ;;
esac
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#      Loop over all polygons.                                                             #
#------------------------------------------------------------------------------------------#
n_submit=0
while [[ ${ff} -lt ${polyz} ]]
do
   let ff=${ff}+1
   let line=${ff}+3


   #---------------------------------------------------------------------------------------#
   #    Format count.                                                                      #
   #---------------------------------------------------------------------------------------#
   if   [[ ${npolys} -ge 10   ]] && [[ ${npolys} -lt 100   ]]
   then
      ffout=$(printf '%2.2i' ${ff})
   elif [[ ${npolys} -ge 100  ]] && [[ ${npolys} -lt 1000  ]]
   then
      ffout=$(printf '%3.3i' ${ff})
   elif [[ ${npolys} -ge 100  ]] && [[ ${npolys} -lt 10000 ]]
   then
      ffout=$(printf '%4.4i' ${ff})
   else
      ffout=${ff}
   fi
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #      Read the ffth line of the polygon list.  There must be smarter ways of doing     #
   # this, but this works.  Here we obtain the polygon name, and its longitude and         #
   # latitude.                                                                             #
   #---------------------------------------------------------------------------------------#
   oi=$(head -${line} ${joborder} | tail -1)
   polyname=$(echo ${oi}      | awk '{print $1  }')
   polyiata=$(echo ${oi}      | awk '{print $2  }')
   polylon=$(echo ${oi}       | awk '{print $3  }')
   polylat=$(echo ${oi}       | awk '{print $4  }')
   yeara=$(echo ${oi}         | awk '{print $5  }')
   montha=$(echo ${oi}        | awk '{print $6  }')
   datea=$(echo ${oi}         | awk '{print $7  }')
   timea=$(echo ${oi}         | awk '{print $8  }')
   yearz=$(echo ${oi}         | awk '{print $9  }')
   monthz=$(echo ${oi}        | awk '{print $10 }')
   datez=$(echo ${oi}         | awk '{print $11 }')
   timez=$(echo ${oi}         | awk '{print $12 }')
   initmode=$(echo ${oi}      | awk '{print $13 }')
   iscenario=$(echo ${oi}     | awk '{print $14 }')
   isizepft=$(echo ${oi}      | awk '{print $15 }')
   iage=$(echo ${oi}          | awk '{print $16 }')
   imaxcohort=$(echo ${oi}    | awk '{print $17 }')
   polyisoil=$(echo ${oi}     | awk '{print $18 }')
   polyntext=$(echo ${oi}     | awk '{print $19 }')
   polysand=$(echo ${oi}      | awk '{print $20 }')
   polyclay=$(echo ${oi}      | awk '{print $21 }')
   polyslsoc=$(echo ${oi}     | awk '{print $22 }')
   polyslph=$(echo ${oi}      | awk '{print $23 }')
   polyslcec=$(echo ${oi}     | awk '{print $24 }')
   polysldbd=$(echo ${oi}     | awk '{print $25 }')
   polydepth=$(echo ${oi}     | awk '{print $26 }')
   polyslhydro=$(echo ${oi}   | awk '{print $27 }')
   polysoilbc=$(echo ${oi}    | awk '{print $28 }')
   polysldrain=$(echo ${oi}   | awk '{print $29 }')
   polycol=$(echo ${oi}       | awk '{print $30 }')
   slzres=$(echo ${oi}        | awk '{print $31 }')
   queue=$(echo ${oi}         | awk '{print $32 }')
   metdriver=$(echo ${oi}     | awk '{print $33 }')
   dtlsm=$(echo ${oi}         | awk '{print $34 }')
   monyrstep=$(echo ${oi}     | awk '{print $35 }')
   iphysiol=$(echo ${oi}      | awk '{print $36 }')
   vmfactc3=$(echo ${oi}      | awk '{print $37 }')
   vmfactc4=$(echo ${oi}      | awk '{print $38 }')
   mphototrc3=$(echo ${oi}    | awk '{print $39 }')
   mphototec3=$(echo ${oi}    | awk '{print $40 }')
   mphotoc4=$(echo ${oi}      | awk '{print $41 }')
   bphotoblc3=$(echo ${oi}    | awk '{print $42 }')
   bphotonlc3=$(echo ${oi}    | awk '{print $43 }')
   bphotoc4=$(echo ${oi}      | awk '{print $44 }')
   kwgrass=$(echo ${oi}       | awk '{print $45 }')
   kwtree=$(echo ${oi}        | awk '{print $46 }')
   gammac3=$(echo ${oi}       | awk '{print $47 }')
   gammac4=$(echo ${oi}       | awk '{print $48 }')
   d0grass=$(echo ${oi}       | awk '{print $49 }')
   d0tree=$(echo ${oi}        | awk '{print $50 }')
   alphac3=$(echo ${oi}       | awk '{print $51 }')
   alphac4=$(echo ${oi}       | awk '{print $52 }')
   klowco2=$(echo ${oi}       | awk '{print $53 }')
   decomp=$(echo ${oi}        | awk '{print $54 }')
   rrffact=$(echo ${oi}       | awk '{print $55 }')
   growthresp=$(echo ${oi}    | awk '{print $56 }')
   lwidthgrass=$(echo ${oi}   | awk '{print $57 }')
   lwidthbltree=$(echo ${oi}  | awk '{print $58 }')
   lwidthnltree=$(echo ${oi}  | awk '{print $59 }')
   q10c3=$(echo ${oi}         | awk '{print $60 }')
   q10c4=$(echo ${oi}         | awk '{print $61 }')
   h2olimit=$(echo ${oi}      | awk '{print $62 }')
   imortscheme=$(echo ${oi}   | awk '{print $63 }')
   ddmortconst=$(echo ${oi}   | awk '{print $64 }')
   cbrscheme=$(echo ${oi}     | awk '{print $65 }')
   isfclyrm=$(echo ${oi}      | awk '{print $66 }')
   icanturb=$(echo ${oi}      | awk '{print $67 }')
   ubmin=$(echo ${oi}         | awk '{print $68 }')
   ugbmin=$(echo ${oi}        | awk '{print $69 }')
   ustmin=$(echo ${oi}        | awk '{print $70 }')
   gamm=$(echo ${oi}          | awk '{print $71 }')
   gamh=$(echo ${oi}          | awk '{print $72 }')
   tprandtl=$(echo ${oi}      | awk '{print $73 }')
   ribmax=$(echo ${oi}        | awk '{print $74 }')
   atmco2=$(echo ${oi}        | awk '{print $75 }')
   thcrit=$(echo ${oi}        | awk '{print $76 }')
   smfire=$(echo ${oi}        | awk '{print $77 }')
   ifire=$(echo ${oi}         | awk '{print $78 }')
   fireparm=$(echo ${oi}      | awk '{print $79 }')
   ipercol=$(echo ${oi}       | awk '{print $80 }')
   runoff=$(echo ${oi}        | awk '{print $81 }')
   imetrad=$(echo ${oi}       | awk '{print $82 }')
   ibranch=$(echo ${oi}       | awk '{print $83 }')
   icanrad=$(echo ${oi}       | awk '{print $84 }')
   ihrzrad=$(echo ${oi}       | awk '{print $85 }')
   crown=$(echo   ${oi}       | awk '{print $86 }')
   ltransvis=$(echo ${oi}     | awk '{print $87 }')
   lreflectvis=$(echo ${oi}   | awk '{print $88 }')
   ltransnir=$(echo ${oi}     | awk '{print $89 }')
   lreflectnir=$(echo ${oi}   | awk '{print $90 }')
   orienttree=$(echo ${oi}    | awk '{print $91 }')
   orientgrass=$(echo ${oi}   | awk '{print $92 }')
   clumptree=$(echo ${oi}     | awk '{print $93 }')
   clumpgrass=$(echo ${oi}    | awk '{print $94 }')
   igoutput=$(echo ${oi}      | awk '{print $95 }')
   ivegtdyn=$(echo ${oi}      | awk '{print $96 }')
   ihydro=$(echo ${oi}        | awk '{print $97 }')
   istemresp=$(echo ${oi}     | awk '{print $98 }')
   istomata=$(echo ${oi}      | awk '{print $99 }')
   iplastic=$(echo ${oi}      | awk '{print $100}')
   icarbonmort=$(echo ${oi}   | awk '{print $101}')
   ihydromort=$(echo ${oi}    | awk '{print $102}')
   igndvap=$(echo ${oi}       | awk '{print $103}')
   iphen=$(echo ${oi}         | awk '{print $104}')
   iallom=$(echo ${oi}        | awk '{print $105}')
   ieconomics=$(echo ${oi}    | awk '{print $106}')
   igrass=$(echo ${oi}        | awk '{print $107}')
   ibigleaf=$(echo ${oi}      | awk '{print $108}')
   integscheme=$(echo ${oi}   | awk '{print $109}')
   nsubeuler=$(echo ${oi}     | awk '{print $110}')
   irepro=$(echo ${oi}        | awk '{print $111}')
   treefall=$(echo ${oi}      | awk '{print $112}')
   ianthdisturb=$(echo ${oi}  | awk '{print $113}')
   ianthdataset=$(echo ${oi}  | awk '{print $114}')
   slscale=$(echo ${oi}       | awk '{print $115}')
   slyrfirst=$(echo ${oi}     | awk '{print $116}')
   slnyrs=$(echo ${oi}        | awk '{print $117}')
   bioharv=$(echo ${oi}       | awk '{print $118}')
   skidarea=$(echo ${oi}      | awk '{print $119}')
   skiddbhthresh=$(echo ${oi} | awk '{print $120}')
   skidsmall=$(echo ${oi}     | awk '{print $121}')
   skidlarge=$(echo ${oi}     | awk '{print $122}')
   fellingsmall=$(echo ${oi}  | awk '{print $123}')
   #---------------------------------------------------------------------------------------#



   #------ Last month and year for monthly-based scripts. ---------------------------------#
   if [[ ${monthz} -eq 1 ]]
   then
      rm_monthz=12
      let rm_yearz=${yearz}-1
   else
      let rm_monthz=${monthz}-1
      rm_yearz=${yearz}
   fi
   #------ Update the time. ---------------------------------------------------------------#
   let rm_whenz=12*${rm_yearz}+${rm_monthz}
   #---------------------------------------------------------------------------------------#


   #----- Find time and minute. -----------------------------------------------------------#
   houra=$(echo ${timea}  | awk '{print substr($1,1,2)}')
   minua=$(echo ${timea}  | awk '{print substr($1,3,2)}')
   hourz=$(echo ${timez}  | awk '{print substr($1,1,2)}')
   minuz=$(echo ${timez}  | awk '{print substr($1,3,2)}')
   #---------------------------------------------------------------------------------------#


   #----- Retrieve some information from ED2IN. -------------------------------------------#
   ED2IN="${here}/${polyname}/ED2IN"
   iphysiol=$(grep -i NL%IPHYSIOL          ${ED2IN} | awk '{print $3}')
   iallom=$(grep   -i NL%IALLOM            ${ED2IN} | awk '{print $3}')
   islhydro=$(grep -i NL%SOIL_HYDRO_SCHEME ${ED2IN} | awk '{print $3}')
   metcyca=$(grep  -i NL%METCYC1           ${ED2IN} | awk '{print $3}')
   metcycz=$(grep  -i NL%METCYCF           ${ED2IN} | awk '{print $3}')
   klight=$(grep   -i NL%DDMORT_CONST      ${ED2IN} | awk '{print $3}')
   #---------------------------------------------------------------------------------------#


   #---- Find the forest inventory cycle. -------------------------------------------------#
   case ${polyiata} in
   gyf|s67)
      biocyca=2004
      biocycz=2009
      subcens=1
      ;;
   s67)
      biocyca=2001
      biocycz=2011
      subcens=1
      ;;
   *)
      biocyca=${metcyca}
      biocycz=${metcycz}
      subcens=0
      ;;
   esac
   #---------------------------------------------------------------------------------------#



   #---- The eddy flux tower cycles. ------------------------------------------------------#
   case ${polyiata} in
   gyf)
      eftyeara=2004
      eftyearz=2012
      ;;
   cax)
      eftyeara=1999
      eftyearz=2003
      ;;
   l[0-5][0-3])
      eftyeara=2006
      eftyearz=2016
      ;;
   m34)
      eftyeara=1999
      eftyearz=2006
      ;;
   s67)
      eftyeara=2001
      eftyearz=2010
      ;;
   s77)
      eftyeara=2001
      eftyearz=2005
      ;;
   s83)
      eftyeara=2000
      eftyearz=2003
      ;;
   pnz)
      eftyeara=2004
      eftyearz=2004
      ;;
   ban)
      eftyeara=2004
      eftyearz=2006
      ;;
   rja)
      eftyeara=1999
      eftyearz=2002
      ;;
   fns)
      eftyeara=1999
      eftyearz=2002
      ;;
   pdg)
      eftyeara=2001
      eftyearz=2003
      ;;
   bsb)
      eftyeara=2006
      eftyearz=2011
      ;;
   tb0|tbx)
      eftyeara=2014
      eftyearz=2018
      ;;
   hvd)
      eftyeara=1992
      eftyearz=2003
      ;;
   *)
      eftyeara=${metcyca}
      eftyearz=${metcycz}
      ;;
   esac
   #---------------------------------------------------------------------------------------#



   #---- The eddy flux tower cycles. ------------------------------------------------------#
   case ${polyiata} in
   gyf)
      bioyeara=2004
      bioyearz=2013
      ;;
   s67)
      bioyeara=1999
      bioyearz=2011
      ;;
   *)
      bioyeara=${eftcyca}
      bioyearz=${eftcycz}
      ;;
   esac
   #---------------------------------------------------------------------------------------#



   #---- Cheat and force the met cycle to be the tower cycle. -----------------------------#
   if [[ ${useperiod} == "f" ]]
   then
      metcyca=${eftyeara}
      metcycz=${eftyearz}
   elif [[ ${useperiod} == "b" ]]
   then
      metcyca=${bioyeara}
      metcycz=${bioyearz}
   fi
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Switch years in case this is a specific drought run.                              #
   #---------------------------------------------------------------------------------------#
   if [[ ${droughtmark} == "TRUE" ]]
   then 
      let yeara=${droughtyeara}-1
      let yearz=${droughtyearz}+1
   fi
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Set up the time and output variables according to the script.                     #
   #---------------------------------------------------------------------------------------#
   case ${rscript} in
   read_monthly.r|yearly_ascii.r|monthly_ascii.r|plot_monthly.r|plot_yearly.r|plot_ycomp.r|plot_census.r|plot_povray.r|r10_monthly.r)
      #------------------------------------------------------------------------------------#
      #     Scripts that are based on monthly means.  The set up is the same, the only     #
      # difference is in the output names.                                                 #
      #------------------------------------------------------------------------------------#
      #------ Check which period to use. --------------------------------------------------#
      if [[ ${useperiod} == "t" ]]
      then
         #------ One meteorological cycle.  Check the type of meteorological driver. ------#
         case ${metdriver} in
         ERA5*|ERAINT*|MERRA2*|PGMF3*|Sheffield|WFDE5*|WFDEI*)
            thisyeara=${metcyca}
            thisyearz=${metcycz}
            ;;
         Tanguro_*)
            thisyeara=2014
            thisyearz=2018
            ;;
         *)
            thisyeara=${metcyca}
            thisyearz=${metcycz}
            for i in ${shiftiata}
            do
               if [[ "x${i}" == "x${polyiata}" ]]
               then
                  echo "     -> Shifting met cycle"
                  let metcycle=${metcycz}-${metcyca}+1
                  let deltayr=${shiftcycle}*${metcycle}
                  let thisyeara=${metcyca}+${deltayr}
                  let thisyearz=${metcycz}+${deltayr}
               fi # end [[ ${i} == ${iata} ]]
            done #end for i in ${shiftiata}
            ;;
         esac #  ${metdriver} in
         #---------------------------------------------------------------------------------#

      elif [[ ${useperiod} == "u" ]]
      then
         #----- The user said which period to use. ----------------------------------------#
         thisyeara=${yusera}
         thisyearz=${yuserz}
         #---------------------------------------------------------------------------------#

      elif [[ ${useperiod} == "f" ]]
      then
         #----- The user said to use the eddy flux period. --------------------------------#
         thisyeara=${eftyeara}
         thisyearz=${eftyearz}
         #---------------------------------------------------------------------------------#

      elif [[ ${useperiod} == "b" ]]
      then
         #----- The user said to use the eddy flux period. --------------------------------#
         thisyeara=${bioyeara}
         thisyearz=${bioyearz}
         #---------------------------------------------------------------------------------#

      else
         #----- Grab all years that the simulation is supposed to run. --------------------#
         thisyeara=${yeara}
         thisyearz=${yearz}
         #---------------------------------------------------------------------------------#
      fi # end [[ ${useperiod} == "t" ]]
      #------------------------------------------------------------------------------------#



      #----- Set up months and days. ------------------------------------------------------#
      thismontha=${montha}
      thismonthz=${monthz}
      thisdatea=${datea}
      #------------------------------------------------------------------------------------#



      #----- Check whether or not to submit the task. -------------------------------------#
      if [[ ${rscript} == "plot_census.r" ]] && [[ ${subcens} -eq 0 ]]
      then
         #---- No need to submit the job if plot_census.r and place doesn't have census. --#
         submit_now=false
         #---------------------------------------------------------------------------------#
      elif ${skip_end}
      then
         status="${here}/${polyname}/rdata_month/status_${polyname}.txt"
         if [[ -s ${status} ]]
         then
            #----- Retrieve current status of the post-processing. ------------------------#
            st_yearz=$(cat ${status}  | awk '{print $1}')
            st_monthz=$(cat ${status} | awk '{print $2}')
            let st_whenz=12*${st_yearz}+${st_monthz}
            #------------------------------------------------------------------------------#

            #------------------------------------------------------------------------------#
            #    Compare the processed time with the last time needed for processing.      #
            #------------------------------------------------------------------------------#
            if [[ ${st_whenz} -ge ${rm_whenz} ]]
            then
               #----- Skip submission because it has reached the end. ---------------------#
               submit_now=false
               #---------------------------------------------------------------------------#
            else
               #----- Run script as it has not reached the end yet. -----------------------#
               submit_now=true
               #---------------------------------------------------------------------------#
            fi
            #------------------------------------------------------------------------------#

         else
            #----- File not find, run the script. -----------------------------------------#
            submit_now=true
            #------------------------------------------------------------------------------#
         fi
         #---------------------------------------------------------------------------------#
      else
         #----- Submit job. ---------------------------------------------------------------#
         submit_now=true
         #---------------------------------------------------------------------------------#
      fi
      #------------------------------------------------------------------------------------#
      ;;
   plot_eval_ed.r)
      #------------------------------------------------------------------------------------#
      #     Cheat by changing metcyca and metcycz in case the meteorological driver is     #
      # Petrolina (output variables exist only for 2004, so we don't need to process       #
      # all years).                                                                        #
      #------------------------------------------------------------------------------------#
      if [[ ${metdriver} == "Petrolina" ]]
      then 
         thismetcyca=2004
         thismetcycz=2004
      else
         thismetcyca=${metcyca}
         thismetcycz=${metcycz}
      fi
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     The period should be equivalent to one meteorological driver period, so we     #
      # compare apples to apples.  The ED2 years don't need to match as long as we pick    #
      # one cycle.                                                                         #
      #------------------------------------------------------------------------------------#
      thisyeara=${thismetcyca}
      thisyearz=${thismetcycz}
      for i in ${shiftiata}
      do
         if [[ "x${i}" == "x${polyiata}" ]]
         then
            #----- Always use the true met driver to find the cycle shift. ----------------#
            echo "     -> Shifting met cycle"
            let metcycle=${metcycz}-${metcyca}+1
            let deltayr=${shiftcycle}*${metcycle}
            let thisyeara=${thismetcyca}+${deltayr}
            let thisyearz=${thismetcycz}+${deltayr}
            #------------------------------------------------------------------------------#
         fi # end [[ ${i} == ${iata} ]]
      done #end for i in ${shiftiata}
      #------------------------------------------------------------------------------------#



      #----- Set up months and days. ------------------------------------------------------#
      thismontha=1
      thismonthz=12
      thisdatea=${datea}
      #------------------------------------------------------------------------------------#



      #----- Assume this should be submitted. ---------------------------------------------#
      submit_now=true
      #------------------------------------------------------------------------------------#
      ;;

   plot_budget.r|plot_rk4.r|plot_rk4pc.r|plot_photo.r|reject_ed.r)
      #------------------------------------------------------------------------------------#
      #     Scripts with very high frequency output (dtlsm or shorter).  The first day     #
      # usually has initialisation problems (for example incoming longwave may be zero     #
      # at the first time step), so we normally skip the first day.                        #
      #------------------------------------------------------------------------------------#
      #----- Check whether to use the user choice of year or the default. -----------------#
      if [[ ${useperiod} == "u" ]]
      then
         thisyeara=${yusera}
         thisyearz=${yuserz}
      else
         thisyeara=${yeara}
         thisyearz=${yearz}
      fi
      #------------------------------------------------------------------------------------#



      #----- Set up months and days. ------------------------------------------------------#
      thismontha=${montha}
      thismonthz=${monthz}
      let thisdatea=${datea}+1
      #------------------------------------------------------------------------------------#



      #----- Assume this should be submitted. ---------------------------------------------#
      submit_now=true
      #------------------------------------------------------------------------------------#
      ;;


   whichrun.r|patchprops.r)
      #------------------------------------------------------------------------------------#
      #     Script with time-independent patch properties.  No need to skip anything.      #
      #------------------------------------------------------------------------------------#
      #----- Check whether to use the user choice of year or the default. -----------------#
      if [[ ${useperiod} == "u" ]]
      then
         thisyeara=${yusera}
         thisyearz=${yuserz}
      else
         thisyeara=${yeara}
         thisyearz=${yearz}
      fi
      #------------------------------------------------------------------------------------#



      #----- Set up months and days. ------------------------------------------------------#
      thismontha=${montha}
      thismonthz=${monthz}
      thisdatea=${datea}
      #------------------------------------------------------------------------------------#



      #----- Assume this should be submitted. ---------------------------------------------#
      submit_now=true
      #------------------------------------------------------------------------------------#
      ;;
   plot_daily.r)
      #------------------------------------------------------------------------------------#
      #     Script with daily means.  No need to skip anything.                            #
      #------------------------------------------------------------------------------------#
      #----- Check whether to use the user choice of year or the default. -----------------#
      if [[ ${useperiod} == "u" ]]
      then
         thisyeara=${yusera}
         thisyearz=${yuserz}
      else
         thisyeara=${yeara}
         thisyearz=${yearz}
      fi
      #------------------------------------------------------------------------------------#



      #----- Set up months and days. ------------------------------------------------------#
      thismontha=${montha}
      thismonthz=${monthz}
      thisdatea=${datea}
      #------------------------------------------------------------------------------------#



      #----- Assume this should be submitted. ---------------------------------------------#
      submit_now=true
      #------------------------------------------------------------------------------------#
      ;;

   plot_fast.r)
      #------------------------------------------------------------------------------------#
      #     Script with short-term averages (usually hourly).  No need to skip any-        #
      # thing.                                                                             #
      #------------------------------------------------------------------------------------#
      if [[ ${useperiod} == "u" ]]
      then
         thisyeara=${yusera}
         thisyearz=${yuserz}
      else
         thisyeara=${yeara}
         thisyearz=${yearz}
      fi
      #------------------------------------------------------------------------------------#



      #----- Set up months and days. ------------------------------------------------------#
      thismontha=${montha}
      thismonthz=${monthz}
      thisdatea=${datea}
      #------------------------------------------------------------------------------------#



      #----- Assume this should be submitted. ---------------------------------------------#
      submit_now=true
      #------------------------------------------------------------------------------------#
      ;;
   esac
   #---------------------------------------------------------------------------------------#



   #----- Copy the R script from the Template folder to the local path. -------------------#
   cp -f ${here}/Template/${rscript} ${here}/${polyname}
   scriptnow="${here}/${polyname}/${rscript}"
   #---------------------------------------------------------------------------------------#



   #----- Switch the keywords by the current settings. ------------------------------------#
   sed -i s@thispoly@${polyname}@g             ${scriptnow}
   sed -i s@thisoutroot@${here}@g              ${scriptnow}
   sed -i s@thispath@${here}@g                 ${scriptnow}
   sed -i s@thatpath@${here}@g                 ${scriptnow}
   sed -i s@thisrscpath@${rscpath}@g           ${scriptnow}
   sed -i s@thispovincs@${pov_incs}@g          ${scriptnow}
   sed -i s@thisyeara@${thisyeara}@g           ${scriptnow}
   sed -i s@thismontha@${thismontha}@g         ${scriptnow}
   sed -i s@thisdatea@${thisdatea}@g           ${scriptnow}
   sed -i s@thishoura@${houra}@g               ${scriptnow}
   sed -i s@thisminua@${minua}@g               ${scriptnow}
   sed -i s@thisyearz@${thisyearz}@g           ${scriptnow}
   sed -i s@thismonthz@${thismonthz}@g         ${scriptnow}
   sed -i s@thisdatez@${datez}@g               ${scriptnow}
   sed -i s@thishourz@${hourz}@g               ${scriptnow}
   sed -i s@thisminuz@${minuz}@g               ${scriptnow}
   sed -i s@thisseasonmona@${seasonmona}@g     ${scriptnow}
   sed -i s@myphysiol@${iphysiol}@g            ${scriptnow}
   sed -i s@myallom@${iallom}@g                ${scriptnow}
   sed -i s@myslhydro@${islhydro}@g            ${scriptnow}
   sed -i s@mydroughtmark@${droughtmark}@g     ${scriptnow}
   sed -i s@mydroughtyeara@${droughtyeara}@g   ${scriptnow}
   sed -i s@mydroughtyearz@${droughtyearz}@g   ${scriptnow}
   sed -i s@mymonthsdrought@${monthsdrought}@g ${scriptnow}
   sed -i s@myvarcycle@${varcycle}@g           ${scriptnow}
   sed -i s@thisoutform@${outform}@g           ${scriptnow}
   sed -i s@mydistrib@${usedistrib}@g          ${scriptnow}
   sed -i s@mymetcyca@${metcyca}@g             ${scriptnow}
   sed -i s@mymetcycz@${metcycz}@g             ${scriptnow}
   sed -i s@mybiocyca@${biocyca}@g             ${scriptnow}
   sed -i s@mybiocycz@${biocycz}@g             ${scriptnow}
   sed -i s@myidbhtype@${idbhtype}@g           ${scriptnow}
   sed -i s@mybackground@${background}@g       ${scriptnow}
   sed -i s@mycorrection@${correct_gs}@g       ${scriptnow}
   sed -i s@myiintphoto@${iint_photo}@g        ${scriptnow}
   sed -i s@myklight@${klight}@g               ${scriptnow}
   sed -i s@myefttrim@${efttrim}@g             ${scriptnow}
   sed -i s@myoldgrowth@${oldgrowth}@g         ${scriptnow}
   sed -i s@myeftyeara@${eftyeara}@g           ${scriptnow}
   sed -i s@myeftyearz@${eftyearz}@g           ${scriptnow}
   #---------------------------------------------------------------------------------------#


   #----- Initialise script. --------------------------------------------------------------#
   epostsh="${here}/${polyname}/exec_$(basename ${rscript} .r).sh"
   complete="${here}/${polyname}/eval_load_complete.txt"
   epoststo="${here}/${polyname}/${epostkey}_epost.sto"
   epostste="${here}/${polyname}/${epostkey}_epost.ste"
   epostout="${here}/${polyname}/${epostkey}_epost.out"
   epostexe="R CMD BATCH --no-save --no-restore ${rscript} ${epostout}"
   rm -fr ${epostsh}
   touch ${epostsh}
   chmod u+x ${epostsh}
   #---------------------------------------------------------------------------------------#


   #----- The script header must check which type of submission to use. -------------------#
   case "${cluster}" in
   SDUMONT)
      #----- Task name with the prefix. ---------------------------------------------------#
      eposttask="${polyname}"
      #------------------------------------------------------------------------------------#


      #----- Header doesn't need to have SBATCH instructions. -----------------------------#
      echo "#!/bin/bash"                                                     >> ${epostsh}
      echo "main=\"${here}/${polyname}\""                                    >> ${epostsh}
      echo "complete=\"\${main}/eval_load_complete.txt\""                    >> ${epostsh}
      echo "yeara=${thisyeara}"                                              >> ${epostsh}
      echo "yearz=${thisyearz}"                                              >> ${epostsh}
      echo ""                                                                >> ${epostsh}
      echo ""                                                                >> ${epostsh}
      #------------------------------------------------------------------------------------#

      ;;
   CANNON)
      #----- Task name with the prefix. ---------------------------------------------------#
      eposttask="${epostkey}-${desc}-${polyname}"
      #------------------------------------------------------------------------------------#


      #----- Header must have SBATCH instructions. ----------------------------------------#
      echo "#!/bin/bash"                                                     >> ${epostsh}
      echo "#SBATCH --ntasks=1                   # Number of tasks"          >> ${epostsh}
      echo "#SBATCH --cpus-per-task=1            # Number of CPUs per task"  >> ${epostsh}
      echo "#SBATCH --partition=${global_queue}  # Queue that will run job"  >> ${epostsh}
      if [[ "${reservation}" != "" ]]
      then
         echo "#SBATCH --reservation=${reservation} # Reserved nodes"        >> ${epostsh}
      fi
      echo "#SBATCH --job-name=${eposttask}      # Task name"                >> ${epostsh}
      echo "#SBATCH --mem-per-cpu=${sim_memory}  # Memory per CPU"           >> ${epostsh}
      echo "#SBATCH --time=${runtime}            # Time for job"             >> ${epostsh}
      echo "#SBATCH --output=${epoststo}         # Standard output path"     >> ${epostsh}
      echo "#SBATCH --error=${epostste}          # Standard error path"      >> ${epostsh}
      echo "#SBATCH --chdir=${here}/${polyname}  # Main directory"           >> ${epostsh}
      echo ""                                                                >> ${epostsh}
      echo "#--- Get plenty of memory."                                      >> ${epostsh}
      echo "ulimit -s unlimited"                                             >> ${epostsh}
      echo ""                                                                >> ${epostsh}
      echo "#--- Initial settings."                                          >> ${epostsh}
      echo "here=\"${here}\"                     # Main path"                >> ${epostsh}
      echo "rscript=\"${rscript}\"               # R Script"                 >> ${epostsh}
      echo "rstdout=\"${epostout}\"              # Standard output"          >> ${epostsh}
      echo ""                                                                >> ${epostsh}
      echo "#--- Print information about this job."                          >> ${epostsh}
      echo "echo \"\""                                                       >> ${epostsh}
      echo "echo \"\""                                                       >> ${epostsh}
      echo "echo \"----- Summary of current job -------------------------\"" >> ${epostsh}
      echo "echo \" CPUs per task:   \${SLURM_CPUS_PER_TASK}\""              >> ${epostsh}
      echo "echo \" Job name:        \${SLURM_JOB_NAME}\""                   >> ${epostsh}
      echo "echo \" Job ID:          \${SLURM_JOB_ID}\""                     >> ${epostsh}
      echo "echo \" Queue:           \${SLURM_JOB_PARTITION}\""              >> ${epostsh}
      echo "echo \" Number of nodes: \${SLURM_NNODES}\""                     >> ${epostsh}
      echo "echo \" Number of tasks: \${SLURM_NTASKS}\""                     >> ${epostsh}
      echo "echo \" Memory per CPU:  \${SLURM_MEM_PER_CPU}\""                >> ${epostsh}
      echo "echo \" Memory per node: \${SLURM_MEM_PER_NODE}\""               >> ${epostsh}
      echo "echo \" Node list:       \${SLURM_JOB_NODELIST}\""               >> ${epostsh}
      echo "echo \" Time limit:      \${SLURM_TIMELIMIT}\""                  >> ${epostsh}
      echo "echo \" Std. Output:     \${SLURM_STDOUTMODE}\""                 >> ${epostsh}
      echo "echo \" Std. Error:      \${SLURM_STDERRMODE}\""                 >> ${epostsh}
      echo "echo \"------------------------------------------------------\"" >> ${epostsh}
      echo "echo \"\""                                                       >> ${epostsh}
      echo "echo \"\""                                                       >> ${epostsh}
      echo "echo \"\""                                                       >> ${epostsh}
      echo "echo \"\""                                                       >> ${epostsh}
      echo ""                                                                >> ${epostsh}
      echo ""                                                                >> ${epostsh}
      echo "#--- Define home in case home is not set"                        >> ${epostsh}
      echo "if [[ \"x\${HOME}\" == \"x\" ]]"                                 >> ${epostsh}
      echo "then"                                                            >> ${epostsh}
      echo "   export HOME=\$(echo ~)"                                       >> ${epostsh}
      echo "fi"                                                              >> ${epostsh}
      echo ""                                                                >> ${epostsh}
      echo "#--- Load modules and settings."                                 >> ${epostsh}
      echo ". \${HOME}/.bashrc ${optsrc}"                                    >> ${epostsh}
      echo ""                                                                >> ${epostsh}
      echo "main=\"${here}/${polyname}\""                                    >> ${epostsh}
      echo "complete=\"\${main}/eval_load_complete.txt\""                    >> ${epostsh}
      echo "yeara=${thisyeara}"                                              >> ${epostsh}
      echo "yearz=${thisyearz}"                                              >> ${epostsh}
      echo ""                                                                >> ${epostsh}
      echo ""                                                                >> ${epostsh}
      #------------------------------------------------------------------------------------#

      ;;
   esac
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      plot_eval_ed won't run all at once due to the sheer number of HDF5 files.        #
   # Run it several times until it is complete.                                            #
   #---------------------------------------------------------------------------------------#
   case ${rscript} in
   plot_eval_ed.r)
      #----- Create script that will run R until all files have been read. ----------------#
      echo "let itmax=\${yearz}-\${yeara}+2"                             >> ${epostsh}
      echo ""                                                            >> ${epostsh}
      echo "/bin/rm -fr \${complete}"                                    >> ${epostsh}
      echo "it=0"                                                        >> ${epostsh}
      echo "while [[ ! -s \${complete} ]] && [[ \${it} -lt \${itmax} ]]" >> ${epostsh}
      echo "do"                                                          >> ${epostsh}
      echo "   let it=\${it}+1"                                          >> ${epostsh}
      echo "   sleep 3"                                                  >> ${epostsh}
      echo "   ${epostexe}"                                              >> ${epostsh}
      echo "done"                                                        >> ${epostsh}
      #------------------------------------------------------------------------------------#
      ;;
   *)
      #----- Create script that will run R until all files have been read. ----------------#
      echo "${epostexe}"                                                 >> ${epostsh}
      echo ""                                                            >> ${epostsh}
      #------------------------------------------------------------------------------------#
      ;;
   esac
   #---------------------------------------------------------------------------------------#




   #----- In case of single task jobs, add commands to complete the job. ------------------#
   case "${cluster}" in
   CANNON)
      #----- Wait for the processing to finish before killing it, and write runtime info. -#
      echo ""                                                               >> ${epostsh}
      echo ""                                                               >> ${epostsh}
      echo "#----- Make sure that jobs complete before terminating script"  >> ${epostsh}
      echo "wait"                                                           >> ${epostsh}
      echo ""                                                               >> ${epostsh}
      echo "#----- Report efficiency of this job"                           >> ${epostsh}
      echo "seff \${SLURM_JOBID}"                                           >> ${epostsh}
      echo ""                                                               >> ${epostsh}
      #------------------------------------------------------------------------------------#
      ;;
   esac
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Decide whether or not to submit task.                                             #
   #---------------------------------------------------------------------------------------#
   if ${submit_now}
   then
      #----- Submit script. ---------------------------------------------------------------#
      echo "${ffout} - Copying script ${rscript} to polygon: ${polyname}..."
      #------------------------------------------------------------------------------------#


      #----- Update the list of scripts to be included in the batch. ----------------------#
      let n_submit=${n_submit}+1
      #------------------------------------------------------------------------------------#




      #----- Decide how to submit the job. ------------------------------------------------#
      case "${cluster}" in
      SDUMONT)
         #----- Append task to main job submission script. --------------------------------#
         srun="srun --nodes=1 --ntasks=1"
         srun="${srun} --cpus-per-task=\${SLURM_CPUS_PER_TASK}"
         srun="${srun} --mem-per-cpu=\${SLURM_MEM_PER_CPU}"
         srun="${srun} --job-name=${eposttask}"
         srun="${srun} --chdir=\${here}/${polyname}"
         srun="${srun} --output=\${here}/${polyname}/${epoststo}"
         srun="${srun} --error=\${here}/${polyname}/${epostste}"
         echo "${srun} ${epostsh} &" >> ${sbatch}
         #---------------------------------------------------------------------------------#
         ;;
      CANNON)

         #----- Add task submission command to the main script. ---------------------------#
         echo "echo \" + Submit post-processing for: ${polyname}.\""    >> ${sbatch}
         echo "sbatch ${epostsh}" >> ${sbatch}
         #---------------------------------------------------------------------------------#
         ;;
      esac
      #------------------------------------------------------------------------------------#
   else
      #----- Skip submission. -------------------------------------------------------------#
      echo "${ffout} - Skipping submission of ${rscript} for polygon: ${polyname}..."
      #------------------------------------------------------------------------------------#
   fi
   #---------------------------------------------------------------------------------------#
done
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#      Make sure job list doesn't request too many nodes.                                  #
#------------------------------------------------------------------------------------------#
if ! ${overcommit} && [[ ${n_submit} -gt ${n_tasks_max} ]]
then
   echo " Number of jobs to submit: ${n_submit}"
   echo " Maximum number of tasks in queue ${global_queue}: ${n_tasks_max}"
   echo " Reduce the number of simulations or try another queue..."
   exit 99
elif ${submit}
then
   echo " Submitting jobs."
   #------- How we submit depends on the cluster. -----------------------------------------#
   case "${cluster}" in
   SDUMONT)
      #------------------------------------------------------------------------------------#
      #     Single multi-task job.  Submit the main script using sbatch.                   #
      #------------------------------------------------------------------------------------#
      sbatch ${sbatch}
      #------------------------------------------------------------------------------------#
      ;;
   CANNON)
      #------------------------------------------------------------------------------------#
      #     Multiple single-task jobs.  Run the script, the sbatch commands are in there.  #
      #------------------------------------------------------------------------------------#
      ${sbatch}
      #------------------------------------------------------------------------------------#
      ;;
   esac
   #---------------------------------------------------------------------------------------#
else
   #------- Provide instructions to the user for a later submission. ----------------------#
   case "${cluster}" in
   SDUMONT)
      #----- Instruct the user to use sbatch. ---------------------------------------------#
      echo "-------------------------------------------------------------------------"
      echo " To submit the simulations, you must use sbatch. "
      echo " Copy and paste the following command in the terminal. "
      echo " "
      echo "       sbatch ${sbatch}"
      echo " "
      #------------------------------------------------------------------------------------#
      ;;
   CANNON)
      #----- Instruct the user NOT to use sbatch. -----------------------------------------#
      echo "-------------------------------------------------------------------------"
      echo " WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! "
      echo " WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! "
      echo " WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! "
      echo " WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! "
      echo "-------------------------------------------------------------------------"
      echo " Do NOT use sbatch ${sbatch}."
      echo " Instead, call the following script directly in your terminal."
      echo " "
      echo "       ${sbatch}"
      echo " "
      #------------------------------------------------------------------------------------#
      ;;
   esac
   #---------------------------------------------------------------------------------------#
fi
#------------------------------------------------------------------------------------------#
