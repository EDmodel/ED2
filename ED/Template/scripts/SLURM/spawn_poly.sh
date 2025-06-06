#!/bin/bash

#==========================================================================================#
#==========================================================================================#
#    Main settings:                                                                        #
#------------------------------------------------------------------------------------------#
#----- Main path, usually set by $(pwd) so you don't need to change it. -------------------#
here=$(pwd)
#----- User name, usually set by $(whoami) so you don't need to change it. ----------------#
moi=$(whoami)
#----- Description of this simulation, used to create unique job names. -------------------#
desc=$(basename ${here})
#----- Source path for ED code and executable. --------------------------------------------#
srcpath="${HOME}/EDBRAMS/ED"
#----- Original and scratch main data paths. ----------------------------------------------#
ordinateur=$(hostname -s)
d_path="${HOME}/data"
#----- Path where biomass initialisation files are: ---------------------------------------#
bioinit="${d_path}/ed2_data/site_bio_data"
alsinit="${d_path}/ed2_data/lidar_spline_bio_data"
intinit="${d_path}/ed2_data/lidar_intensity_bio_data"
lutinit="${d_path}/ed2_data/lidar_lookup_bio_data"
ebainit="${d_path}/ed2_data/lidar_eba_bio_data"
biotype=0      # 0 -- "default" setting (isizepft controls default/nounder)
               # 1 -- isizepft controls number of PFTs, whereas iage controls patches.
               # 2 -- airborne lidar initialisation using return counts ("default"). 
               # 3 -- airborne lidar initialisation using intensity counts.
               # 4 -- airborne lidar/inventory hybrid initialisation ("lookup table"). 
               # 5 -- airborne lidar (EBA Transects).
               # For lidar initialisation (2-4), isizepft is the disturbance history key.
#----- Path and file prefix for init_mode = 5. --------------------------------------------#
restart="${d_path}/ed2_data/restarts_XXX"
#----- File containing the list of jobs and their settings: -------------------------------#
joborder="${here}/joborder.txt"
#----- This is the header with the Sheffield data. ----------------------------------------#
shefhead='SHEF_NCEP_DRIVER_DS314'
#----- Path with drivers for each scenario. -----------------------------------------------#
metmaindef="${d_path}/ed2_data"
packdatasrc="${d_path}/to_scratch"
#----- Path with land use scenarios. ------------------------------------------------------#
lumain="${d_path}/ed2_data/lcluc_scenarios"
#----- Path with other input data bases (soil texture, DGD, land mask, etc). --------------#
inpmain="${d_path}/ed2_data"
#----- Should the met driver be copied to local scratch disks? ----------------------------#
copy2scratch=false
#------------------------------------------------------------------------------------------#

#------------------------------------------------------------------------------------------#
#      History run variables.                                                              #
#------------------------------------------------------------------------------------------#
#----- Force history run (0 = no, 1 = yes). -----------------------------------------------#
forcehisto=0
#----- Path with the history file to be used. ---------------------------------------------#
fullygrown="${HOME}/Simulations/Debug/D001_Debug/xyz_settings/histo/xyz_settings"
#----- Time that we shall use. ------------------------------------------------------------#
yearh="1510"  # Year
monthh="07"   # Month
dateh="01"    # Day
timeh="0000"  # Hour
#----- Default tolerance. -----------------------------------------------------------------#
toldef="0.01"
#----- Initialisation scripts. ------------------------------------------------------------#
initrc="${HOME}/.bashrc"          # Initialisation script for most nodes
#----- Initialisation scripts. ------------------------------------------------------------#
optsrc="-n"                   # Option for .bashrc (for special submission settings)
                              #   In case none is needed, leave it blank ("").
#----- Submit job automatically? (It may become false if something prevents submission). --#
submit=false      # This checks whether to submit runs or not.
on_the_fly=false  # In case submit=true and this is the CANNON cluster, on_the_fly allows
                  #   directories are jobs to be submitted as the directories are ready.
#----- Executable names. ------------------------------------------------------------------#
execname="ed_2.2-opt"            # Normal executable, for most queues
checkexec=true                   # Check executable
#----- Settings for this group of polygons. -----------------------------------------------#
global_queue="shared,huce_intel" # Queue
partial=false                    # Partial submission (false will ignore polya and npartial
                                 #    and send all polygons.
polya=21                         # First polygon to submit
npartial=300                     # Maximum number of polygons to include in this bundle
                                 #    (actual number will be adjusted for total number of 
                                 #     polygons if needed be).
dttask=2                         # Time to wait between task submission
init_only=false                  # Run model initialisation only?
                                 #    This can be useful when the initialisation step
                                 #    demands much more memory than the time steps (e.g.,
                                 #    when using lidar data to initialisate the model).
if ${init_only}
then
   sim_memory=64000              # Memory per simulation. Zero uses queue's default
   n_cpt=1                       # Number of cpus per task (Zero uses queue's maximum)
   runtime="01-00:00:00"         # Requested runtime.  Zero uses the queue's maximum.
else
   sim_memory=3072               # Memory per simulation. Zero uses queue's default
   n_cpt=12                      # Number of cpus per task (Zero uses queue's maximum)
   runtime="07-00:00:00"         # Requested runtime.  Zero uses the queue's maximum.
fi
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


#----- squeue settings to check whether a directory exists. -------------------------------#
sqout="%.200j %.8T"
squeue="squeue --noheader -u ${moi}"
#------------------------------------------------------------------------------------------#


#----- Find out which platform we are using. ----------------------------------------------#
host=$(hostname -s)
case "${host}" in
rclogin*|holy*|moorcroft*|rcnx*)
   #----- Cluster is CANNON (formerly known as ODYSSEY). ----------------------------------#
   cluster="CANNON"
   #---------------------------------------------------------------------------------------#
   ;;
sdumont*)
   #----- Cluster is SDUMONT. Ironically perhaps, we cannot submit jobs on the fly. -------#
   cluster="SDUMONT"
   on_the_fly=false
   #---------------------------------------------------------------------------------------#
   ;;
*)
   echo -n "Failed guessing cluster from node name.  Please type the name:   "
   read cluster
   cluster=$(echo ${cluster} | tr '[:lower:]' '[:upper:]')
   case "${cluster}" in
   CANNON|CANNO|CANN|CAN|CA|C)
      #----- Set cluster to CANNON. -------------------------------------------------------#
      cluster="CANNON"
      #------------------------------------------------------------------------------------#
      ;;
   ODYSSEY|ODYSSE|ODYSS|ODYS|ODY|OD|O)
      #----- Set cluster to CANNON. -------------------------------------------------------#
      cluster="CANNON"
      echo " - Odyssey is now known as Cannon.  Using Cannon settings."
      #------------------------------------------------------------------------------------#
      ;;
   SDUMONT|SDUMON|SDUMO|SDUM|SDU|SD|S)
      #----- Cluster is SDUMONT. Ironically perhaps, we cannot submit jobs on the fly. ----#
      cluster="SDUMONT"
      on_the_fly=false
      #------------------------------------------------------------------------------------#
      ;;
   *)
      #----- Unknown cluster. -------------------------------------------------------------#
      echo " Cluster ${cluster} is not recognised.  Pick either CANNON or SDUMONT."
      echo "    (or edit the script to add the new cluster)."
      exit 92
      #------------------------------------------------------------------------------------#
      ;;
   esac
   ;;
esac
#------------------------------------------------------------------------------------------#

#------------------------------------------------------------------------------------------#
#     Do not let on_the_fly if submit is false.                                            #
#------------------------------------------------------------------------------------------#
if ! ${submit}
then
   on_the_fly=false
fi
#------------------------------------------------------------------------------------------#


#----- Set the main path for the site, pseudo past and Sheffield met drivers. -------------#
if ${copy2scratch}
then
   metmain="/scratch/$(whoami)"
else
   metmain=${metmaindef}
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


#------------------------------------------------------------------------------------------#
#   Change parameters for first and last polygons.                                         #
#------------------------------------------------------------------------------------------#
while [[ ${#} > 0 ]]
do
   key="${1}"
   case ${key} in
   -a)
      polya=${2}
      shift 2
      ;;
   -n)
      npartial=${2}
      shift 2
      ;;
   -p)
      partial=true
      shift 1
      ;;
   -f)
      partial=false
      shift 1
      ;;
   *)
      echo "Nothing" > /dev/null
      shift 1
      ;;
   esac
done

#------------------------------------------------------------------------------------------#
#   Check whether the executable is copied.  If not, let the user know and stop the        #
# script.                                                                                  #
#------------------------------------------------------------------------------------------#
exec_full="${here}/executable/${execname}"
exec_src="${srcpath}/build/${execname}"
if [[ ! -s ${exec_full} ]]
then
   echo "Executable file : ${exec_full} is not in the executable directory"
   echo "Copy the executable to the file before running this script!"
   exit 99
elif ${checkexec}
then
   #----- Check whether the executable is up to date. -------------------------------------#
   idiff=$(diff ${exec_full} ${exec_src} 2> /dev/null | wc -l)
   if [[ ${idiff} -gt 0 ]]
   then
      #----- Executable is not up to date, warn user and offer to update it. --------------#
      echo "Executable ${exec_full} is not the same as the most recently compiled file."
      echo "Most recent: ${exec_src}."
      echo "Do you want to update the executable? [y/N]"
      read update
      update=$(echo ${update} | tr '[:upper:]' '[:lower:]')
      case "${update}" in
         y*) rsync - Putv ${exec_src} ${exec_full} ;;
      esac
      #------------------------------------------------------------------------------------#
   fi
   #---------------------------------------------------------------------------------------#
fi
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
      n_cpt_max=12
      n_cpn=24
      runtime_max="31-00:00:00"
      node_memory=64000
      ;;
   cpu|nvidia|phi)
      n_nodes_max=50
      n_cpt_max=12
      n_cpn=24
      runtime_max="2-00:00:00"
      node_memory=64000
      ;;
   cpu_dev)
      n_nodes_max=20
      n_cpt_max=12
      n_cpn=24
      runtime_max="02:00:00"
      node_memory=64000
      ;;
   nvidia_dev|phi_dev)
      n_nodes_max=2
      n_cpt_max=12
      n_cpn=24
      runtime_max="02:00:00"
      node_memory=64000
      ;;
   cpu_scal|nvidia_scal)
      n_nodes_max=128
      n_cpt_max=12
      n_cpn=24
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
      n_cpt_max=18
      n_cpn=32
      runtime_max="14-00:00:00"
      node_memory=244660
      ;;
   "huce_cascade")
      n_nodes_max=30
      n_cpt_max=48
      n_cpn=24
      runtime_max="14-00:00:00"
      node_memory=183240
      ;;
   "huce_intel")
      n_nodes_max=150
      n_cpt_max=12
      n_cpn=24
      runtime_max="14-00:00:00"
      node_memory=122090
      ;;
   "serial_requeue")
      n_nodes_max=900
      n_cpt_max=24
      n_cpn=48
      runtime_max="7-00:00:00"
      node_memory=30362
      ;;
   "shared,huce_intel"|"huce_intel,shared")
      n_nodes_max=150
      n_cpt_max=12
      n_cpn=24
      runtime_max="7-00:00:00"
      node_memory=122090
      ;;
   "shared")
      n_nodes_max=220
      n_cpt_max=24
      n_cpn=48
      runtime_max="7-00:00:00"
      node_memory=183240
      ;;
   "test")
      n_nodes_max=9
      n_cpt_max=24
      n_cpn=48
      runtime_max="08:00:00"
      node_memory=183240
      ;;
   "unrestricted")
      n_nodes_max=4
      n_cpt_max=24
      n_cpn=48
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


#----- Check that the resources requested are reasonable. ---------------------------------#
if [[ ${n_cpt} -gt ${n_cpt_max} ]]
then
   echo " Too many CPUs per task requested:"
   echo " Queue                   = ${global_queue}"
   echo " Maximum CPUs per task   = ${n_cpt_max}"
   echo " Requested CPUs per task = ${n_cpt}"
   exit 99
else
   if [[ ${n_cpt} -eq 0 ]]
   then
      n_cpt=${n_cpt_max}
   fi
   let n_tasks_max=${n_nodes_max}*${n_cpn}
   let n_tasks_max=${n_tasks_max}/${n_cpt}
fi
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
#   Make sure memory does not exceed maximum amount that can be requested.                 #
#------------------------------------------------------------------------------------------#
if [[ ${sim_memory} -eq 0 ]]
then
   let sim_memory=${node_memory}/${n_cpn}
   let node_memory=${n_cpn}*${sim_memory}
elif [[ ${sim_memory} -gt ${node_memory} ]]
then 
   echo "Simulation memory ${sim_memory} cannot exceed node memory ${node_memory}!"
   exit 99
fi
#------------------------------------------------------------------------------------------#


#----- Limit the memory per CPU based on the simulation memory. ---------------------------#
let cpu_memory=${sim_memory}/${n_cpt}
#------------------------------------------------------------------------------------------#



#---- Partial or complete. ----------------------------------------------------------------#

if ${partial}
then
   let ff=${polya}-1
   let polyz=${ff}+${npartial}
   if [[ ${polyz} -gt ${npolys} ]]
   then
      polyz=${npolys}
   fi
   pfmt="%${ndig}.${ndig}i"
   partlabel=$(printf "${pfmt}" ${polya})-$(printf "${pfmt}" ${polyz})
   sbatch="${here}/sub_batch_${partlabel}.sh"
   obatch="${here}/out_batch_${partlabel}.log"
   ebatch="${here}/err_batch_${partlabel}.log"
   jobname="${desc}-sims_${partlabel}"
else
   ff=0
   polya=1
   polyz=${npolys}
   sbatch="${here}/sub_batch.sh"
   obatch="${here}/out_batch.log"
   ebatch="${here}/err_batch.log"
   jobname="${desc}-sims"
fi
let ntasks=1+${polyz}-${polya}
#------------------------------------------------------------------------------------------#


#----- Summary for this submission preparation.  Then give 5 seconds for user to cancel. --#
echo "------------------------------------------------"
echo "  Submission summary: "
echo ""
echo "  Memory per task:     ${sim_memory}"
echo "  Memory per cpu:      ${cpu_memory}"
echo "  CPUs per node:       ${n_cpn}"
echo "  CPUs per task:       ${n_cpt}"
echo "  Partition:           ${global_queue}"
echo "  Run time:            ${runtime}"
echo "  First polygon:       ${polya}"
echo "  Last polygon:        ${polyz}"
echo "  Potl. task count:    ${ntasks}"
echo "  Total polygon count: ${npolys}"
echo " "
echo " Partial submission:   ${partial}"
echo " Automatic submission: ${submit}"
echo " Script:               $(basename ${sbatch})"
echo "------------------------------------------------"
echo ""
echo -n " Waiting five seconds before proceeding... "
sleep 5
echo "Done!"
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#    Initialise executable.  This depends on the cluster because in some of them it is     #
# better to submit a single multi-task job, whereas in others it makes more sense to       #
# submit individual jobs.                                                                  #
#------------------------------------------------------------------------------------------#
case "${cluster}" in
SDUMONT)
   #---------------------------------------------------------------------------------------#
   #   Check whether there is already a job submitted that looks like this one.            #
   #---------------------------------------------------------------------------------------#
   queued=$(${squeue} -o "${outform}" | grep ${jobname} | wc -l)
   if [[ ${queued} -gt 0 ]]
   then
      echo "There is already a job called \"${jobname}\" running."
      echo "New submissions must have different names: be creative!"
      exit 99
   fi
   #---------------------------------------------------------------------------------------#



   #----- Initialise script to submit a single multi-task job. ----------------------------#
   rm -f ${sbatch}
   touch ${sbatch}
   chmod u+x ${sbatch}
   echo "#!/bin/bash" >> ${sbatch}
   echo "#SBATCH --ntasks=myntasks               # Number of tasks"            >> ${sbatch}
   echo "#SBATCH --cpus-per-task=${n_cpt}        # Number of CPUs per task"    >> ${sbatch}
   echo "#SBATCH --cores-per-socket=${n_cpt}     # Min. # of cores per socket" >> ${sbatch}
   echo "#SBATCH --partition=${global_queue}     # Queue that will run job"    >> ${sbatch}
   echo "#SBATCH --job-name=${jobname}           # Job name"                   >> ${sbatch}
   echo "#SBATCH --mem-per-cpu=${cpu_memory}     # Memory per CPU"             >> ${sbatch}
   echo "#SBATCH --time=${runtime}               # Time for job"               >> ${sbatch}
   echo "#SBATCH --output=${obatch}              # Standard output path"       >> ${sbatch}
   echo "#SBATCH --error=${ebatch}               # Standard error path"        >> ${sbatch}
   echo ""                                                                     >> ${sbatch}
   echo "#--- Initial settings."                                               >> ${sbatch}
   echo "here=\"${here}\"                            # Main path"              >> ${sbatch}
   echo "exec=\"${exec_full}\"                       # Executable"             >> ${sbatch}
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
   echo "#--- Set OpenMP parameters"                                           >> ${sbatch}
   echo "if [ \"\${SLURM_CPUS_PER_TASK}\" == \"\" ]"                           >> ${sbatch}
   echo "then"                                                                 >> ${sbatch}
   echo "   export OMP_NUM_THREADS=1"                                          >> ${sbatch}
   echo "else"                                                                 >> ${sbatch}
   echo "   export OMP_NUM_THREADS=\${SLURM_CPUS_PER_TASK}"                    >> ${sbatch}
   echo "fi"                                                                   >> ${sbatch}
   echo ""                                                                     >> ${sbatch}
   echo "#----- Task list."                                                    >> ${sbatch}
   #---------------------------------------------------------------------------------------#
   ;;
CANNON)
   #----- Initialise script to submit multiple single-task jobs. --------------------------#
   rm -f ${sbatch}
   touch ${sbatch}
   chmod u+x ${sbatch}
   echo "#!/bin/bash" >> ${sbatch}
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
   echo "#----- Job list."                                                     >> ${sbatch}
   #---------------------------------------------------------------------------------------#
   ;;
esac
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Loop over all polygons.                                                              #
#------------------------------------------------------------------------------------------#
n_submit=0
while [[ ${ff} -lt ${polyz} ]]
do
   let ff=${ff}+1
   let line=${ff}+3


   #---------------------------------------------------------------------------------------#
   #    Format count.                                                                      #
   #---------------------------------------------------------------------------------------#
   if   [[ ${npolys} -lt 100   ]]
   then
      ffout=$(printf '%2.2i' ${ff})
   elif [[ ${npolys} -lt 1000  ]]
   then
      ffout=$(printf '%3.3i' ${ff})
   elif [[ ${npolys} -lt 10000 ]]
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



   #---------------------------------------------------------------------------------------#
   #      In case this is an "init_only" run, impose final time to be the initial time.    #
   #---------------------------------------------------------------------------------------#
   if ${init_only}
   then
      yearz=${yeara}
      monthz=${montha}
      datez=${datea}
      timez=${timea}
   fi
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #    Task name.  This depends on the type of submission (just to avoid clutter).        #
   #---------------------------------------------------------------------------------------#
   case "${cluster}" in
   SDUMONT) 
      #----- Task name is bound to the global jobname, no need to add prefix. -------------#
      taskname="${polyname}"
      #------------------------------------------------------------------------------------#
      ;;
   CANNON )
      #----- Add prefix to the task name, which will name this job. -----------------------#
      taskname="${desc}-${polyname}"
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #   Check whether there is already a job submitted that looks like this one.         #
      #------------------------------------------------------------------------------------#
      queued=$(${squeue} -o "${outform}" | grep ${taskname} | wc -l)
      if [[ ${queued} -gt 0 ]]
      then
         echo "There is already a job called \"${taskname}\" running."
         echo "New submissions must have different names: be creative!"
         exit 99
      fi
      #------------------------------------------------------------------------------------#
      ;;
   esac
   #---------------------------------------------------------------------------------------#





   #---------------------------------------------------------------------------------------#
   #     Decide units for restart files.                                                   #
   #---------------------------------------------------------------------------------------#
   let elapseda=12*${yeara}+${montha}
   let elapsedz=12*${yearz}+${monthz}
   let nmonths=${elapsedz}-${elapseda}
   if [[ ${nmonths} -ge 240 ]]
   then
      iunitstate=3
   elif [[ ${nmonths} -ge 3 ]]
   then
      iunitstate=2
   else
      iunitstate=1
   fi
   #---------------------------------------------------------------------------------------#


   #----- Check whether the directories exist or not, and stop the script if they do. -----#
   if [[ -s ${here}/${polyname} ]]
   then
      echo -n "${ffout} ${polyname}: updating files..."
      
      #----- Save the last tolerance in case we are going to make it more strict. ---------#
      oldtol=$(grep NL%RK4_TOLERANCE ${here}/${polyname}/ED2IN | awk '{print $3}')
      rm -f ${here}/${polyname}/ED2IN 
      cp ${here}/Template/ED2IN         ${here}/${polyname}/ED2IN
   else
      echo -n "${ffout} ${polyname}: creating directory..."
      #----- Copy the Template directory to a unique polygon directory. -------------------#
      oldtol=""
      cp -r ${here}/Template ${here}/${polyname}
   fi
   #---------------------------------------------------------------------------------------#






   #---------------------------------------------------------------------------------------#
   #   Make sure that we have a reasonable tolerance.                                      #
   #---------------------------------------------------------------------------------------#
   if [[ "x${oldtol}" == "xmytoler" ]] || [[ "x${oldtol}" == "x" ]]
   then
      toler=${toldef}
   else
      toler=${oldtol}
   fi
   #---------------------------------------------------------------------------------------#






   #---------------------------------------------------------------------------------------#
   #      Here we check whether the an old directory exists or not.  If it doesn't, then   #
   # we simply copy the template, and assume initial run.  Otherwise, we must find out     #
   # where the simulation was when it stopped.                                             #
   #---------------------------------------------------------------------------------------#
   if [[ -s  ${here}/${polyname} ]]
   then

      #------------------------------------------------------------------------------------#
      #      This step is necessary because we may have killed the run while it was        #
      # writing, and as a result, the file may be corrupt.                                 #
      #------------------------------------------------------------------------------------#
      nhdf5=$(ls -1 ${here}/${polyname}/histo/*.h5 2> /dev/null | wc -l)
      if [[ ${nhdf5} -gt 0 ]]
      then
         h5fine=0

         while [[ ${h5fine} -eq 0 ]]
         do
            lasthdf5=$(ls -1 ${here}/${polyname}/histo/*.h5 | tail -1)
            h5dump -H ${lasthdf5} 1> /dev/null 2> ${here}/badfile.txt

            if [[ -s ${here}/badfile.txt ]]
            then
               /bin/rm -fv ${lasthdf5}
               nhdf5=$(ls -1 ${here}/${polyname}/histo/*.h5 2> /dev/null | wc -l)
               if [[ ${nhdf5} -eq 0 ]]
               then
                  h5fine=1
               fi
            else
               h5fine=1
            fi

            /bin/rm -f ${here}/badfile.txt
         done
         #---------------------------------------------------------------------------------#
      fi
      #------------------------------------------------------------------------------------#
   fi
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Run the small R script to check whether the simulation was running or not, and   #
   # whether there was any cohort left by the time the runs were stopped.                  #
   #---------------------------------------------------------------------------------------#
   /bin/rm -f ${here}/${polyname}/statusrun.txt
   /bin/rm -f ${here}/${polyname}/whichrun.r
   /bin/cp -f ${here}/Template/whichrun.r ${here}/${polyname}/whichrun.r
   whichrun="${here}/${polyname}/whichrun.r"
   outwhich="${here}/${polyname}/outwhichrun.txt"
   sed -i~ s@thispoly@${polyname}@g           ${whichrun}
   sed -i~ s@thisqueue@${queue}@g             ${whichrun}
   sed -i~ s@pathhere@${here}@g               ${whichrun}
   sed -i~ s@paththere@${here}@g              ${whichrun}
   sed -i~ s@thisyeara@${yeara}@g             ${whichrun}
   sed -i~ s@thismontha@${montha}@g           ${whichrun}
   sed -i~ s@thisdatea@${datea}@g             ${whichrun}
   sed -i~ s@thistimea@${timea}@g             ${whichrun}
   sed -i~ s@thischecksteady@FALSE@g          ${whichrun}
   sed -i~ s@thismetcyc1@${metcyc1}@g         ${whichrun}
   sed -i~ s@thismetcycf@${metcycf}@g         ${whichrun}
   sed -i~ s@thisnyearmin@10000@g             ${whichrun}
   sed -i~ s@thisststcrit@0.0@g               ${whichrun}
   R CMD BATCH --no-save --no-restore ${whichrun} ${outwhich}
   while [[ ! -s ${here}/${polyname}/statusrun.txt ]]
   do
      sleep 0.2
   done
   year=$(cat  ${here}/${polyname}/statusrun.txt | awk '{print $2}')
   month=$(cat ${here}/${polyname}/statusrun.txt | awk '{print $3}')
   date=$(cat  ${here}/${polyname}/statusrun.txt | awk '{print $4}')
   time=$(cat  ${here}/${polyname}/statusrun.txt | awk '{print $5}')
   runt=$(cat  ${here}/${polyname}/statusrun.txt | awk '{print $6}')
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Run the small R script to generate specific observation times.                   #
   #---------------------------------------------------------------------------------------#
   /bin/rm -f ${here}/${polyname}/gen_obstimes.r
   /bin/cp -f ${here}/Template/gen_obstimes.r ${here}/${polyname}/gen_obstimes.r
   obstimes="${here}/${polyname}/gen_obstimes.r"
   outtimes="${here}/${polyname}/out_obstimes.txt"
   sed -i~ s@thispoly@${polyname}@g           ${obstimes}
   sed -i~ s@thisqueue@${queue}@g             ${obstimes}
   sed -i~ s@pathhere@${here}@g               ${obstimes}
   sed -i~ s@paththere@${here}@g              ${obstimes}
   sed -i~ s@thisyeara@${yeara}@g             ${obstimes}
   sed -i~ s@thismontha@${montha}@g           ${obstimes}
   sed -i~ s@thisdatea@${datea}@g             ${obstimes}
   sed -i~ s@thistimea@${timea}@g             ${obstimes}
   sed -i~ s@thisyearz@${yearz}@g             ${obstimes}
   sed -i~ s@thismonthz@${monthz}@g           ${obstimes}
   sed -i~ s@thisdatez@${datez}@g             ${obstimes}
   sed -i~ s@thistimez@${timez}@g             ${obstimes}
   sed -i~ s@thislon@${polylon}@g             ${obstimes}
   sed -i~ s@thislat@${polylat}@g             ${obstimes}
   R CMD BATCH --no-save --no-restore ${obstimes} ${outtimes}
   while [[ ! -s ${here}/${polyname}/${polyname}_obstimes.txt ]]
   do
      sleep 0.2
   done
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #    To ensure simulations can be requeued, we no longer set RUNTYPE to INITIAL or      #
   # HISTORY (except when we should force history).  Instead, we select RESTORE and let    #
   # the model decide between initial or history.                                          #
   #---------------------------------------------------------------------------------------#
   if [[ "${runt}" == "INITIAL" ]] || [[ "${runt}" == "HISTORY" ]]
   then
      runt="RESTORE"
   fi
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Determine which PFTs to use based on the "iata" code, isizepft, and biotype.      #
   #---------------------------------------------------------------------------------------#
   case ${biotype} in
   2)
      #------------------------------------------------------------------------------------#
      #     Airborne lidar.  isizepft is not controlling isizepft, but disturbance         #
      # history. Initialise PFTs using the default settings.                               #
      #------------------------------------------------------------------------------------#
      case ${polyiata} in
      tzi|zmh|nqn)
         pfts="1,7,8,9,10,11,15,16"
         pasture=16
         crop=16
         plantation=15
         logging="7,8,9,10,11,15"
         probharv="1.0,1.0,1.0,1.0,1.0,1.0"
         dbhharv="0.0,0.0,0.0,0.0,0.0,0.0"
         seedharv=0.75
         storharv=0.00
         leafharv=0.00
         ;;
      hvd|wch|tqh)
         pfts="6,8,9,10,11,16"
         pasture=16
         crop=16
         plantation=11
         logging="6,8,9,10,11"
         probharv="1.0,1.0,1.0,1.0,1.0"
         dbhharv="0.0,0.0,0.0,0.0,0.0"
         seedharv=0.75
         storharv=0.00
         leafharv=0.00
         ;;
      l[0-5][0-3])
         pfts="6,8,9,10,11"
         pasture=16
         crop=16
         plantation=6
         logging="6,8,9,10,11"
         probharv="1.0,1.0,1.0,1.0,1.0"
         dbhharv="0.0,0.0,0.0,0.0,0.0"
         seedharv=0.75
         storharv=0.00
         leafharv=0.00
         ;;
      asu|cnf|bnu|cwb|erm|iqq|ipv|mgf|rao|sla|zpe|kna|sfn)
         pfts="1,2,3,4,15,16"
         pasture=1
         crop=16
         plantation=15
         logging="2,3,4,15"
         probharv="1.0,1.0,1.0,1.0,1.0,1.0"
         dbhharv="50.,50.,50.,50.,50.,50."
         seedharv=0.75
         storharv=0.00
         leafharv=0.00
         ;;
      fns*)
         pfts="1,16"
         pasture=1
         crop=1
         plantation=13
         logging="13"
         probharv="1.0"
         dbhharv="50."
         seedharv=0.75
         storharv=0.00
         leafharv=0.00
         ;;
      s77*)
         pfts="1,16"
         pasture=1
         crop=16
         plantation=13
         logging="13"
         probharv="1.0"
         dbhharv="50."
         seedharv=0.75
         storharv=0.00
         leafharv=0.00
         ;;
      *)
         pfts="1,2,3,4,16"
         pasture=1
         crop=1
         plantation=13
         logging="2,3,4"
         probharv="1.0,1.0,1.0,1.0,1.0"
         dbhharv="50.,50.,50.,50.,50."
         seedharv=0.75
         storharv=0.00
         leafharv=0.00
         ;;
      esac
      ;;
      #------------------------------------------------------------------------------------#
   *)
      case ${polyiata} in
      tzi|zmh|nqn)
         pfts="1,7,8,9,10,11,15,16"
         pasture=16
         crop=16
         plantation=15
         logging="7,8,9,10,11,15"
         probharv="1.0,1.0,1.0,1.0,1.0,1.0"
         dbhharv="0.0,0.0,0.0,0.0,0.0,0.0"
         seedharv=0.75
         storharv=0.00
         leafharv=0.00
         ;;
      hvd|wch|tqh)
         pfts="6,8,9,10,11,16"
         pasture=16
         crop=16
         plantation=10
         logging="6,8,9,10,11"
         probharv="1.0,1.0,1.0,1.0,1.0"
         dbhharv="0.0,0.0,0.0,0.0,0.0"
         seedharv=0.75
         storharv=0.00
         leafharv=0.00
         ;;
      l[0-5][0-3])
         pfts="6,8,9,10,11"
         pasture=16
         crop=16
         plantation=6
         logging="6,8,9,10,11"
         probharv="1.0,1.0,1.0,1.0,1.0"
         dbhharv="0.0,0.0,0.0,0.0,0.0"
         seedharv=0.75
         storharv=0.00
         leafharv=0.00
         ;;
      asu|cnf|bnu|cwb|erm|iqq|ipv|mgf|rao|sla|zpe|kna|sfn)
         pfts="1,2,3,4,15,16"
         pasture=1
         crop=16
         plantation=15
         logging="2,3,4,15"
         probharv="1.0,1.0,1.0,1.0,1.0,1.0"
         dbhharv="50.,50.,50.,50.,50.,50."
         seedharv=0.75
         storharv=0.00
         leafharv=0.00
         ;;
      fns*)
         pfts="1,16"
         pasture=1
         crop=1
         plantation=13
         logging="13"
         probharv="1.0"
         dbhharv="50."
         seedharv=0.75
         storharv=0.00
         leafharv=0.00
         ;;
      s77*)
         pfts="1,16"
         pasture=1
         crop=16
         plantation=13
         logging="13"
         probharv="1.0"
         dbhharv="50."
         seedharv=0.75
         storharv=0.00
         leafharv=0.00
         ;;
      *)
         pfts="1,2,3,4,16"
         pasture=1
         crop=1
         plantation=13
         logging="2,3,4"
         probharv="1.0,1.0,1.0,1.0,1.0"
         dbhharv="50.,50.,50.,50.,50."
         seedharv=0.75
         storharv=0.00
         leafharv=0.00
         ;;
      esac
      ;;
      #------------------------------------------------------------------------------------#
   esac
   #---------------------------------------------------------------------------------------#





   #---------------------------------------------------------------------------------------#
   #     Determine the census structure.                                                   #
   #---------------------------------------------------------------------------------------#
   let yodd=${yeara}%2
   case ${polyiata} in
   gyf)
      dtcensus=24
      let yr1stcensus=${yeara}+${yodd}
      mon1stcensus=3
      minrecruitdbh=10
      ;;
   s67|dcm)
      dtcensus=24
      let yr1stcensus=${yeara}+1-${yodd}
      mon1stcensus=7
      minrecruitdbh=10
      ;;
   *)
      dtcensus=1
      yr1stcensus=${yeara}
      mon1stcensus=1
      minrecruitdbh=10
      ;;
   esac
   #---------------------------------------------------------------------------------------#





   #---------------------------------------------------------------------------------------#
   #     Determine the scenario paths.                                                     #
   #---------------------------------------------------------------------------------------#
   case ${iscenario} in
   default)
      case ${metdriver} in
      ERA5_CHIRPS)
         #----- ERA5 (CHIRPS precipitation, Brazilian Amazon only). -----------------------#
         scentype="ERA5"
         iscenario="ERA5_SOUTHAM_CHIRPS"
         ;;
      ERA5_CHIRPS)
         #----- ERA-Interim (CHIRPS precipitation). ---------------------------------------#
         scentype="ERA_Interim"
         iscenario="ERAINT_SOUTHAM_CHIRPS"
         ;;
      ERAINT_CHIRPS)
         #----- ERA-Interim (CHIRPS precipitation). ---------------------------------------#
         scentype="ERA_Interim"
         iscenario="ERAINT_SOUTHAM_CHIRPS"
         ;;
      ERAINT_MSWEP2)
         #----- ERA-Interim (MSWEP2 precipitation). ---------------------------------------#
         scentype="ERA_Interim"
         iscenario="ERAINT_SOUTHAM_MSWEP2"
         ;;
      ERAINT_NATIVE)
         #----- ERA-Interim (native precipitation). ---------------------------------------#
         scentype="ERA_Interim"
         iscenario="ERAINT_SOUTHAM_NATIVE"
         ;;
      MERRA2_CHIRPS)
         #----- MERRA2 (CHIRPS precipitation). --------------------------------------------#
         scentype="MERRA2"
         iscenario="MERRA2_SOUTHAM_CHIRPS"
         ;;
      MERRA2_MSWEP2)
         #----- MERRA2 (MSWEP2 precipitation). --------------------------------------------#
         scentype="MERRA2"
         iscenario="MERRA2_SOUTHAM_MSWEP2"
         ;;
      MERRA2_NATIVE)
         #----- MERRA-2 (native precipitation). -------------------------------------------#
         scentype="MERRA2"
         iscenario="MERRA2_SOUTHAM_NATIVE"
         ;;
      PGMF3_CHIRPS)
         #----- PGMF-3 (CHIRPS precipitation). --------------------------------------------#
         scentype="PGMF3"
         iscenario="PGMF3_SOUTHAM_CHIRPS"
         ;;
      PGMF3_MSWEP2)
         #----- PGMF-3 (CHIRPS precipitation). --------------------------------------------#
         scentype="PGMF3"
         iscenario="PGMF3_SOUTHAM_MSWEP2"
         ;;
      PGMF3_NATIVE)
         #----- PGMF-3 (native precipitation). --------------------------------------------#
         scentype="PGMF3"
         iscenario="PGMF3_SOUTHAM_NATIVE"
         ;;
      Sheffield)
         #----- Sheffield. ----------------------------------------------------------------#
         scentype="sheffield"
         iscenario="sheffield"
         ;;
      WFDE5_CHIRPS)
         #----- WFDEI (CHIRPS Precipitation). ---------------------------------------------#
         scentype="WFDE5"
         iscenario="WFDE5_SOUTHAM_CHIRPS"
         ;;
      WFDEI_CHIRPS)
         #----- WFDEI (CHIRPS Precipitation). ---------------------------------------------#
         scentype="WFDEI"
         iscenario="WFDEI_SOUTHAM_CHIRPS"
         ;;
      WFDEI_CRUP)
         #----- WFDEI (CRU Precipitation). ------------------------------------------------#
         scentype="WFDEI"
         iscenario="WFDEI_SOUTHAM_CRUP"
         ;;
      WFDEI_GPCC)
         #----- WFDEI (GPCC Precipitation). -----------------------------------------------#
         scentype="WFDEI"
         iscenario="WFDEI_SOUTHAM_GPCC"
         ;;
      WFDEI_MSWEP2)
         #----- WFDEI (MSWEP2 Precipitation). ---------------------------------------------#
         scentype="WFDEI"
         iscenario="WFDEI_SOUTHAM_MSWEP2"
         ;;
      *)
         #----- Tower data. ---------------------------------------------------------------#
         scentype="wmo+eft"
         iscenario="eft"
         ;;
      esac
      ;;
   eft|wmo|shr)
      #----- Tower data, keep scenario as is. ---------------------------------------------#
      scentype="wmo+eft"
      ;;
   *)
      #----- Rainfall scenario, keep scenario as is. --------------------------------------#
      scentype="realisation_scen_driver"
      ;;
   esac
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      Choose the scenario to use.  Iscenario follows the following convention:         #
   # "default"    -- No scenario.  Use the tower/Sheffield data.                           #
   # "wmo"        -- No scenario.  Use the WMO-based data.                                 #
   # "rRRR_tTTT   -- Use scenarios, with rRRRR controlling the rainfall, and tTTTT         #
   #                 controlling temperature.                                              #
   #                                                                                       #
   # rRRR         -- Rainfall scenarios, where rRRRR means:                                #
   #                 r+000: INMET-based time series, no resampling.                        #
   #                 r+010: Re-sampling with substitution, but equal chances for all years #
   #                 r-XXX: Shift the location (similar to mean) of the distribution by    #
   #                        -X.XX units of scale (similar to standard deviation), so the   #
   #                        time series becomes drier.                                     #
   #                 r+XXX: Similar to above, but make the time series wetter.             #
   #                                                                                       #
   # tTTT         -- Temperature scenarios, where tTTTT means:                             #
   #                 t+000: No change in temperature                                       #
   #                 t-YYY: Change temperature by -Y.YY Kelvin.  Keep relative humidity    #
   #                        the same and correct specific humidity.                        #
   #                 t+YYY: Change temperature by +Y.YY Kelvin.  Keep relative humidity    #
   #                        the same and correct by +Y.YY Kelvin.                          #
   #                 r+XXX: Similar to above, but make the time series wetter.             #
   #---------------------------------------------------------------------------------------#
   #----- Find out which scenario to use. -------------------------------------------------#
   fullscen="${metmain}/met_driver/${scentype}/${iscenario}"
   #---------------------------------------------------------------------------------------#
   #     Determine which meteorological data set to use.  Default is the Sheffield/NCEP    #
   # dataset, otherwise the site-level tower data is used.                                 #
   #---------------------------------------------------------------------------------------#
   case ${metdriver} in
   Bananal)
      metdriverdb="${fullscen}/Bananal/Bananal_HEADER"
      metcyc1=2004
      metcycf=2006
      ;;
   Brasilia)
      metdriverdb="${fullscen}/Brasilia/Brasilia_HEADER"
      metcyc1=2006
      metcycf=2012
      ;;
   Caxiuana)
      metdriverdb="${fullscen}/Caxiuana/Caxiuana_HEADER"
      metcyc1=1999
      metcycf=2003
      ;;
   ERA5_CHIRPS)
      metdriverdb="${fullscen}/${iscenario}_HEADER"
      metcyc1=1981
      metcycf=2019
      ;;
   ERAINT_CHIRPS)
      metdriverdb="${fullscen}/${iscenario}_HEADER"
      metcyc1=1981
      metcycf=2017
      ;;
   ERAINT_NATIVE)
      metdriverdb="${fullscen}/${iscenario}_HEADER"
      metcyc1=1979
      metcycf=2017
      ;;
   ERAINT_MSWEP2)
      metdriverdb="${fullscen}/${iscenario}_HEADER"
      metcyc1=1979
      metcycf=2016
      ;;
   Fazenda_Nossa_Senhora)
      metdriverdb="${fullscen}/Fazenda_Nossa_Senhora/Fazenda_Nossa_Senhora_HEADER"
      metcyc1=1999
      metcycf=2002
      ;;
   Harvard)
      metdriverdb="${fullscen}/Harvard/Harvard_HEADER"
      metcyc1=1992
      metcycf=2003
      ;;
   Laegern)
      metdriverdb="${fullscen}/Laegern/Laegern_HEADER"
      metcyc1=2006
      metcycf=2016
      ;;
   Manaus_Km34)
      metdriverdb="${fullscen}/Manaus_Km34/Manaus_Km34_HEADER"
      metcyc1=1999
      metcycf=2006
      ;;
   MERRA2_CHIRPS)
      metdriverdb="${fullscen}/${iscenario}_HEADER"
      metcyc1=1981
      metcycf=2017
      ;;
   MERRA2_MSWEP2)
      metdriverdb="${fullscen}/${iscenario}_HEADER"
      metcyc1=1980
      metcycf=2016
      ;;
   MERRA2_NATIVE)
      metdriverdb="${fullscen}/${iscenario}_HEADER"
      metcyc1=1980
      metcycf=2017
      ;;
   Natal)
      metdriverdb="${fullscen}/Natal/Natal_HEADER"
      metcyc1=2009
      metcycf=2012
      ;;
   Paracou)
      metdriverdb="${fullscen}/Paracou/Paracou_HEADER"
      metcyc1=2004
      metcycf=2014
      ;;
   Pe-de-Gigante)
      metdriverdb="${fullscen}/Pe-de-Gigante/Pe-de-Gigante_HEADER"
      metcyc1=2001
      metcycf=2003
      ;;
   Petrolina)
      metdriverdb="${fullscen}/Petrolina/Petrolina_HEADER"
      metcyc1=2004
      metcycf=2012
      ;;
   PGMF3_CHIRPS)
      metdriverdb="${fullscen}/${iscenario}_HEADER"
      metcyc1=1981
      metcycf=2016
      ;;
   PGMF3_MSWEP2)
      metdriverdb="${fullscen}/${iscenario}_HEADER"
      metcyc1=1979
      metcycf=2016
      ;;
   PGMF3_NATIVE)
      metdriverdb="${fullscen}/${iscenario}_HEADER"
      metcyc1=1979
      metcycf=2016
      ;;
   Rebio_Jaru)
      metdriverdb="${fullscen}/Rebio_Jaru/Rebio_Jaru_HEADER"
      metcyc1=1999
      metcycf=2002
      ;;
   Santarem_Km67)
      metdriverdb="${fullscen}/Santarem_Km67/Santarem_Km67_HEADER"
      metcyc1=2001
      metcycf=2011
      ;;
   Santarem_Km77)
      metdriverdb="${fullscen}/Santarem_Km77/Santarem_Km77_HEADER"
      metcyc1=2001
      metcycf=2005
      ;;
   Santarem_Km83)
      metdriverdb="${fullscen}/Santarem_Km83/Santarem_Km83_HEADER"
      metcyc1=2000
      metcycf=2003
      ;;
   Sheffield)
      metdriverdb=${fullscen}/${shefhead}
      metcyc1=1969
      metcycf=2008
      ;;
   Tanguro_Burn)
      metdriverdb="${fullscen}/Tanguro_Burn/Tanguro_Burn_HEADER"
      metcyc1=2008
      metcycf=2018
      ;;
   Tanguro_Ctrl)
      metdriverdb="${fullscen}/Tanguro_Ctrl/Tanguro_Ctrl_HEADER"
      metcyc1=2008
      metcycf=2018
      ;;
   Tonzi)
      metdriverdb="${fullscen}/Tonzi/Tonzi_HEADER"
      metcyc1=2000
      metcycf=2010
      ;;
   WFDE5_CHIRPS)
      metdriverdb="${fullscen}/${iscenario}_HEADER"
      metcyc1=1981
      metcycf=2018
      ;;
   WFDEI_CHIRPS)
      metdriverdb="${fullscen}/${iscenario}_HEADER"
      metcyc1=1981
      metcycf=2016
      ;;
   WFDEI_CRUP)
      metdriverdb="${fullscen}/${iscenario}_HEADER"
      metcyc1=1979
      metcycf=2016
      ;;
   WFDEI_GPCC)
      metdriverdb="${fullscen}/${iscenario}_HEADER"
      metcyc1=1979
      metcycf=2016
      ;;
   WFDEI_MSWEP2)
      metdriverdb="${fullscen}/${iscenario}_HEADER"
      metcyc1=1979
      metcycf=2016
      ;;
   *)
      echo "Met driver: ${metdriver}"
      echo "Sorry, this met driver is not valid for regular runs"
      exit 85
      ;;
   esac
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Correct years so it is not tower-based or Sheffield.                              #
   #---------------------------------------------------------------------------------------#
   case ${iscenario} in
   default|eft|shr|ERA5*|ERAINT*|MERRA2*|PGMF3*|Sheffield|WFDE5*|WFDEI*)
      echo "Nothing" > /dev/null
      ;;
   *)
      metcyc1=1972
      metcycf=2012
      imetavg=1
      ;;
   esac
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Set the land use data set.                                                        #
   #---------------------------------------------------------------------------------------#
   case ${ianthdataset} in
   glu-331)
      case ${polyiata} in
      tzi|zmh|nqn|hvd|wch|tqh|l[0-5][0-3])
         ludatabase="${lumain}/glu/glu-"
         ;;
      *)
         ludatabase="${lumain}/glu-3.3.1/glu-3.3.1-"
         ;;
      esac
      ;;
   glu-sa1)
      case ${polyiata} in
      tzi|zmh|nqn|hvd|wch|tqh|l[0-5][0-3])
         ludatabase="${lumain}/glu/glu-"
         ;;
      *)
         ludatabase="${lumain}/glu-3.3.1+sa1.bau/glu-3.3.1+sa1.bau-"
         ;;
      esac
      ;;
   glu-sag)
      case ${polyiata} in
      tzi|zmh|nqn|hvd|wch|tqh|l[0-5][0-3])
         ludatabase="${lumain}/glu/glu-"
         ;;
      *)
         ludatabase="${lumain}/glu-3.3.1+sa1.gov/glu-3.3.1+sa1.gov-"
         ;;
      esac
      ;;
   glu-sa2)
      case ${polyiata} in
      tzi|zmh|nqn|hvd|wch|tqh|l[0-5][0-3])
         ludatabase="${lumain}/glu/glu-"
         ;;
      *)
         ludatabase="${lumain}/glu-3.3.1+sa2.bau/glu-3.3.1+sa2.bau-"
         ;;
      esac
      ;;
   sa2-ril)
      ludatabase="${lumain}/SimAmazonia2/ril/sa2_ril_"
      ;;
   sa2-cvl)
      ludatabase="${lumain}/SimAmazonia2/cvl/sa2_cvl_"
      ;;
   lurcp26)
      ludatabase="${lumain}/luh-1.1+rcp26_image/luh-1.1+rcp26_image-"
      ;;
   lurcp45)
      ludatabase="${lumain}/luh-1.1+rcp45_minicam/luh-1.1+rcp45_minicam-"
      ;;
   lurcp60)
      ludatabase="${lumain}/luh-1.1+rcp60_aim/luh-1.1+rcp60_aim-"
      ;;
   lurcp85)
      ludatabase="${lumain}/luh-1.1+rcp85_message/luh-1.1+rcp85_message-"
      ;;
   *)
      #------------------------------------------------------------------------------------#
      #     Stop the script if anthropogenic dataset is invalid and this is a simulation   #
      # with anthropogenic disturbance.                                                    #
      #------------------------------------------------------------------------------------#
      if [[ ${ianthdisturb} -eq 1 ]]
      then
         echo " Polygon:       ${polyname}"
         echo " IATA:          ${polyiata}"
         echo " IANTH_DATASET: ${iage}"
         echo 'Invalid anthropogenic disturbance data set!'
         exit 53
      fi
      #------------------------------------------------------------------------------------#
      ;;
   esac
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Define whether we use the met cycle to define the first and last year, or the     #
   # default year.                                                                         #
   #---------------------------------------------------------------------------------------#
   if [[ ${yeara} -eq 0 ]]
   then
      thisyeara=${metcyc1}
   else
      thisyeara=${yeara}
   fi
   if [[ ${yearz} -eq 0 ]]
   then
      thisyearz=${metcycf}
   else
      thisyearz=${yearz}
   fi
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Determine which soil profile to use.  We have eight categories (A-H), and the     #
   # soil resolution.                                                                      #
   #---------------------------------------------------------------------------------------#
   polyslz1=""
   polyslz2=""
   polyslz3=""
   polyslz4=""
   polyslz5=""
   polyslz6=""
   polyslz7=""
   polyslm1=""
   polyslm2=""
   polyslm3=""
   polyslm4=""
   polyslm5=""
   polyslm6=""
   polyslm7=""
   polyslt1=""
   polyslt2=""
   polyslt3=""
   polyslt4=""
   polyslt5=""
   polyslt6=""
   polyslt7=""
   case ${slzres} in
   0)
      case ${polydepth} in
      A)
         polynzg=16
         polyslz1="-1.250,-1.154,-1.059,-0.966,-0.875,-0.785,-0.697,-0.612,-0.529,-0.448,"
         polyslz2="-0.370,-0.295,-0.224,-0.156,-0.095,-0.040"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslm2=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslt2=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000"
         ;;
      B)
         polynzg=16
         polyslz1="-2.000,-1.826,-1.657,-1.492,-1.333,-1.179,-1.030,-0.888,-0.752,-0.623,"
         polyslz2="-0.501,-0.388,-0.283,-0.188,-0.106,-0.040"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslm2=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslt2=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000"
         ;;
      C)
         polynzg=16
         polyslz1="-3.000,-2.713,-2.437,-2.171,-1.917,-1.674,-1.443,-1.225,-1.019,-0.828,"
         polyslz2="-0.651,-0.490,-0.346,-0.221,-0.118,-0.040"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslm2=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslt2=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000"
         ;;
      D)
         polynzg=16
         polyslz1="-4.000,-3.593,-3.204,-2.833,-2.481,-2.147,-1.832,-1.538,-1.265,-1.013,"
         polyslz2="-0.784,-0.579,-0.400,-0.248,-0.126,-0.040"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslm2=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslt2=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000"
         ;;
      E)
         polynzg=16
         polyslz1="-5.000,-4.468,-3.963,-3.483,-3.030,-2.604,-2.205,-1.836,-1.495,-1.185,"
         polyslz2="-0.906,-0.660,-0.447,-0.271,-0.134,-0.040"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslm2=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslt2=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000"
         ;;
      F)
         polynzg=16
         polyslz1="-6.000,-5.339,-4.714,-4.123,-3.567,-3.048,-2.566,-2.121,-1.714,-1.347,"
         polyslz2="-1.019,-0.733,-0.490,-0.291,-0.140,-0.040"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslm2=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslt2=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000"
         ;;
      G)
         polynzg=16
         polyslz1="-7.000,-6.207,-5.458,-4.755,-4.096,-3.483,-2.917,-2.397,-1.925,-1.501,"
         polyslz2="-1.126,-0.802,-0.529,-0.310,-0.145,-0.040"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslm2=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslt2=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000"
         ;;
      H)
         polynzg=16
         polyslz1="-8.000,-7.072,-6.198,-5.380,-4.617,-3.910,-3.259,-2.664,-2.127,-1.648,"
         polyslz2="-1.228,-0.866,-0.566,-0.326,-0.150,-0.040"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslm2=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslt2=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000"
         ;;
      I)
         polynzg=16
         polyslz1="-9.000,-7.934,-6.934,-5.999,-5.131,-4.329,-3.593,-2.925,-2.324,-1.790,"
         polyslz2="-1.325,-0.928,-0.600,-0.342,-0.155,-0.040"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslm2=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslt2=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000"
         ;;
      J)
         polynzg=16
         polyslz1="-10.50,-9.223,-8.029,-6.919,-5.891,-4.946,-4.084,-3.305,-2.609,-1.995,"
         polyslz2="-1.464,-1.015,-0.648,-0.364,-0.161,-0.040"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslm2=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslt2=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000"
         ;;
      K)
         polynzg=16
         polyslz1="-12.00,-10.51,-9.118,-7.828,-6.640,-5.552,-4.563,-3.674,-2.883,-2.191"
         polyslz2="-1.595,-1.096,-0.693,-0.383,-0.166,-0.040"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslm2=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslt2=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000"
         ;;
      P)
         polynzg=7
         polyslz1="-0.580,-0.420,-0.300,-0.200,-0.120,-0.060,-0.020"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000"
         ;;
      Q)
         polynzg=8
         polyslz1="-0.830,-0.630,-0.450,-0.300,-0.200,-0.120,-0.060,-0.020"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000"
         ;;
      U)
         polynzg=6
         polyslz1="-0.400,-0.300,-0.200,-0.120,-0.060,-0.020"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000"
         ;;
      V)
         polynzg=6
         polyslz1="-0.420,-0.300,-0.200,-0.120,-0.060,-0.020"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000"
         ;;
      W)
         polynzg=7
         polyslz1="-0.510,-0.400,-0.300,-0.200,-0.120,-0.060,-0.020"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000"
         ;;
      X)
         polynzg=8
         polyslz1="-0.670,-0.520,-0.400,-0.300,-0.200,-0.120,-0.060,-0.020"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000"
         ;;
      Y)
         polynzg=8
         polyslz1="-0.700,-0.550,-0.400,-0.300,-0.200,-0.120,-0.060,-0.020"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000"
         ;;
      Z)
         polynzg=8
         polyslz1="-0.750,-0.600,-0.450,-0.300,-0.200,-0.120,-0.060,-0.020"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000"
         ;;
      *)
         polynzg=16
         polyslz1="-8.000,-7.072,-6.198,-5.380,-4.617,-3.910,-3.259,-2.664,-2.127,-1.648,"
         polyslz2="-1.228,-0.866,-0.566,-0.326,-0.150,-0.040"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslm2=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslt2=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000"
         ;;
      esac
      ;;
   1)
      case ${polydepth} in
      A)
         polynzg=16
         polyslz1="-1.250,-1.135,-1.024,-0.917,-0.814,-0.715,-0.620,-0.530,-0.445,-0.364,"
         polyslz2="-0.289,-0.221,-0.158,-0.103,-0.056,-0.020"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslm2=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslt2=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000"
         ;;
      B)
         polynzg=16
         polyslz1="-2.000,-1.797,-1.602,-1.417,-1.240,-1.073,-0.916,-0.769,-0.632,-0.507,"
         polyslz2="-0.392,-0.290,-0.200,-0.124,-0.063,-0.020"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslm2=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslt2=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000"
         ;;
      C)
         polynzg=16
         polyslz1="-3.000,-2.670,-2.357,-2.061,-1.784,-1.524,-1.283,-1.061,-0.857,-0.673,"
         polyslz2="-0.510,-0.367,-0.245,-0.146,-0.070,-0.020"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslm2=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslt2=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000"
         ;;
      D)
         polynzg=16
         polyslz1="-4.000,-3.536,-3.099,-2.690,-2.308,-1.955,-1.629,-1.332,-1.064,-0.824,"
         polyslz2="-0.614,-0.433,-0.283,-0.163,-0.075,-0.020"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslm2=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslt2=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000"
         ;;
      E)
         polynzg=16
         polyslz1="-5.000,-4.397,-3.833,-3.307,-2.819,-2.371,-1.961,-1.590,-1.257,-0.964,"
         polyslz2="-0.709,-0.493,-0.316,-0.178,-0.080,-0.020"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslm2=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslt2=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000"
         ;;
      F)
         polynzg=16
         polyslz1="-6.000,-5.254,-4.559,-3.914,-3.320,-2.776,-2.282,-1.837,-1.442,-1.095,"
         polyslz2="-0.798,-0.548,-0.346,-0.192,-0.083,-0.020"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslm2=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslt2=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000"
         ;;
      G)
         polynzg=16
         polyslz1="-7.000,-6.108,-5.279,-4.514,-3.812,-3.172,-2.593,-2.076,-1.618,-1.221,"
         polyslz2="-0.881,-0.600,-0.374,-0.204,-0.087,-0.020"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslm2=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslt2=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000"
         ;;
      H)
         polynzg=16
         polyslz1="-8.000,-6.959,-5.995,-5.108,-4.296,-3.560,-2.897,-2.307,-1.789,-1.340,"
         polyslz2="-0.961,-0.648,-0.400,-0.215,-0.089,-0.020"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslm2=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslt2=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000"
         ;;
      I)
         polynzg=16
         polyslz1="-9.000,-7.807,-6.706,-5.696,-4.775,-3.942,-3.195,-2.533,-1,954,-1.456,"
         polyslz2="-1.037,-0.694,-0.424,-0.225,-0.092,-0.020"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslm2=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslt2=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000"
         ;;
      J)
         polynzg=16
         polyslz1="-10.50,-9.076,-7.766,-6.569,-5.482,-4.504,-3.631,-2.862,-2.194,-1.622,"
         polyslz2="-1.145,-0.759,-0.458,-0.239,-0.096,-0.020"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslm2=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslt2=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000"
         ;;
      K)
         polynzg=16
         polyslz1="-12.00,-10.34,-8.818,-7.432,-6.179,-5.055,-4.057,-3.182,-2.425,-1.782"
         polyslz2="-1.248,-0.820,-0.490,-0.252,-0.099,-0.020"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslm2=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslt2=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000"
         ;;
      P)
         polynzg=7
         polyslz1="-0.580,-0.420,-0.300,-0.200,-0.120,-0.060,-0.020"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000"
         ;;
      Q)
         polynzg=8
         polyslz1="-0.830,-0.630,-0.450,-0.300,-0.200,-0.120,-0.060,-0.020"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000"
         ;;
      U)
         polynzg=6
         polyslz1="-0.400,-0.300,-0.200,-0.120,-0.060,-0.020"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000"
         ;;
      V)
         polynzg=6
         polyslz1="-0.420,-0.300,-0.200,-0.120,-0.060,-0.020"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000"
         ;;
      W)
         polynzg=7
         polyslz1="-0.510,-0.400,-0.300,-0.200,-0.120,-0.060,-0.020"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000"
         ;;
      X)
         polynzg=8
         polyslz1="-0.670,-0.520,-0.400,-0.300,-0.200,-0.120,-0.060,-0.020"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000"
         ;;
      Y)
         polynzg=8
         polyslz1="-0.700,-0.550,-0.400,-0.300,-0.200,-0.120,-0.060,-0.020"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000"
         ;;
      Z)
         polynzg=8
         polyslz1="-0.750,-0.600,-0.450,-0.300,-0.200,-0.120,-0.060,-0.020"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000"
         ;;
      *)
         polynzg=16
         polyslz1="-8.000,-6.959,-5.995,-5.108,-4.296,-3.560,-2.897,-2.307,-1.789,-1.340,"
         polyslz2="-0.961,-0.648,-0.400,-0.215,-0.089,-0.020"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslm2=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslt2=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000"
         ;;
      esac
      ;;
   2)
      case ${polydepth} in
      A)
         polynzg="30"
         polyslz1="-1.329,-1.246,-1.168,-1.093,-1.022,-0.955,-0.890,-0.829,-0.770,-0.714,"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslz2="-0.661,-0.611,-0.563,-0.517,-0.473,-0.432,-0.392,-0.354,-0.318,-0.284,"
         polyslm2=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslt2=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslz3="-0.252,-0.221,-0.191,-0.163,-0.136,-0.111,-0.086,-0.063,-0.041,-0.020"
         polyslm3=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000"
         polyslt3=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000"
         ;;
      B)
         polynzg="37"
         polyslz1="-2.033,-1.917,-1.806,-1.701,-1.601,-1.506,-1.415,-1.329,-1.246,-1.168,"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslz2="-1.093,-1.022,-0.955,-0.890,-0.829,-0.770,-0.714,-0.661,-0.611,-0.563,"
         polyslm2=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslt2=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslz3="-0.517,-0.473,-0.432,-0.392,-0.354,-0.318,-0.284,-0.252,-0.221,-0.191,"
         polyslm3=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslt3=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslz4="-0.163,-0.136,-0.111,-0.086,-0.063,-0.041,-0.020"
         polyslm4=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000"
         polyslt4=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000"
         ;;
      C)
         polynzg="44"
         polyslz1="-3.023,-2.860,-2.705,-2.557,-2.416,-2.282,-2.154,-2.033,-1.917,-1.806,"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslz2="-1.701,-1.601,-1.506,-1.415,-1.329,-1.246,-1.168,-1.093,-1.022,-0.955,"
         polyslm2=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslt2=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslz3="-0.890,-0.829,-0.770,-0.714,-0.661,-0.611,-0.563,-0.517,-0.473,-0.432,"
         polyslm3=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslt3=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslz4="-0.392,-0.354,-0.318,-0.284,-0.252,-0.221,-0.191,-0.163,-0.136,-0.111,"
         polyslm4=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslt4=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslz5="-0.086,-0.063,-0.041,-0.020"
         polyslm5=" 1.000, 1.000, 1.000, 1.000"
         polyslt5=" 0.000, 0.000, 0.000, 0.000"
         ;;
      D)
         polynzg="50"
         polyslz1="-4.187,-3.969,-3.761,-3.562,-3.374,-3.194,-3.023,-2.860,-2.705,-2.557,"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslz2="-2.416,-2.282,-2.154,-2.033,-1.917,-1.806,-1.701,-1.601,-1.506,-1.415,"
         polyslm2=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslt2=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslz3="-1.329,-1.246,-1.168,-1.093,-1.022,-0.955,-0.890,-0.829,-0.770,-0.714,"
         polyslm3=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslt3=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslz4="-0.661,-0.611,-0.563,-0.517,-0.473,-0.432,-0.392,-0.354,-0.318,-0.284,"
         polyslm4=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslt4=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslz5="-0.252,-0.221,-0.191,-0.163,-0.136,-0.111,-0.086,-0.063,-0.041,-0.020"
         polyslm5=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000"
         polyslt5=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000"
         ;;
      E)
         polynzg="52"
         polyslz1="-4.657,-4.416,-4.187,-3.969,-3.761,-3.562,-3.374,-3.194,-3.023,-2.860,"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslz2="-2.705,-2.557,-2.416,-2.282,-2.154,-2.033,-1.917,-1.806,-1.701,-1.601,"
         polyslm2=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslt2=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslz3="-1.506,-1.415,-1.329,-1.246,-1.168,-1.093,-1.022,-0.955,-0.890,-0.829,"
         polyslm3=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslt3=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslz4="-0.770,-0.714,-0.661,-0.611,-0.563,-0.517,-0.473,-0.432,-0.392,-0.354,"
         polyslm4=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslt4=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslz5="-0.318,-0.284,-0.252,-0.221,-0.191,-0.163,-0.136,-0.111,-0.086,-0.063,"
         polyslm5=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslt5=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslz6="-0.041,-0.020"
         polyslm6=" 1.000, 1.000"
         polyslt6=" 0.000, 0.000"
         ;;
      F)
         polynzg="58"
         polyslz1="-6.157,-5.907,-5.657,-5.407,-5.157,-4.907,-4.657,-4.416,-4.187,-3.969,"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslz2="-3.761,-3.562,-3.374,-3.194,-3.023,-2.860,-2.705,-2.557,-2.416,-2.282,"
         polyslm2=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslt2=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslz3="-2.154,-2.033,-1.917,-1.806,-1.701,-1.601,-1.506,-1.415,-1.329,-1.246,"
         polyslm3=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslt3=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslz4="-1.168,-1.093,-1.022,-0.955,-0.890,-0.829,-0.770,-0.714,-0.661,-0.611,"
         polyslm4=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslt4=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslz5="-0.563,-0.517,-0.473,-0.432,-0.392,-0.354,-0.318,-0.284,-0.252,-0.221,"
         polyslm5=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslt5=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslz6="-0.191,-0.163,-0.136,-0.111,-0.086,-0.063,-0.041,-0.020"
         polyslm6=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000"
         polyslt6=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000"
         ;;
      G)
         polynzg="62"
         polyslz1="-7.157,-6.907,-6.657,-6.407,-6.157,-5.907,-5.657,-5.407,-5.157,-4.907,"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslz2="-4.657,-4.416,-4.187,-3.969,-3.761,-3.562,-3.374,-3.194,-3.023,-2.860,"
         polyslm2=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslt2=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslz3="-2.705,-2.557,-2.416,-2.282,-2.154,-2.033,-1.917,-1.806,-1.701,-1.601,"
         polyslm3=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslt3=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslz4="-1.506,-1.415,-1.329,-1.246,-1.168,-1.093,-1.022,-0.955,-0.890,-0.829,"
         polyslm4=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslt4=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslz5="-0.770,-0.714,-0.661,-0.611,-0.563,-0.517,-0.473,-0.432,-0.392,-0.354,"
         polyslm5=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslt5=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslz6="-0.318,-0.284,-0.252,-0.221,-0.191,-0.163,-0.136,-0.111,-0.086,-0.063,"
         polyslm6=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslt6=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslz7="-0.041,-0.020"
         polyslm7=" 1.000, 1.000"
         polyslt7=" 0.000, 0.000"
         ;;
      H)
         polynzg="66"
         polyslz1="-8.157,-7.907,-7.657,-7.407,-7.157,-6.907,-6.657,-6.407,-6.157,-5.907,"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslz2="-5.657,-5.407,-5.157,-4.907,-4.657,-4.416,-4.187,-3.969,-3.761,-3.562,"
         polyslm2=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslt2=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslz3="-3.374,-3.194,-3.023,-2.860,-2.705,-2.557,-2.416,-2.282,-2.154,-2.033,"
         polyslm3=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslt3=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslz4="-1.917,-1.806,-1.701,-1.601,-1.506,-1.415,-1.329,-1.246,-1.168,-1.093,"
         polyslm4=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslt4=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslz5="-1.022,-0.955,-0.890,-0.829,-0.770,-0.714,-0.661,-0.611,-0.563,-0.517,"
         polyslm5=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslt5=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslz6="-0.473,-0.432,-0.392,-0.354,-0.318,-0.284,-0.252,-0.221,-0.191,-0.163,"
         polyslm6=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslt6=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslz7="-0.136,-0.111,-0.086,-0.063,-0.041,-0.020"
         polyslm7=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000"
         polyslt7=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000"
         ;;
      *)
         polynzg="66"
         polyslz1="-8.157,-7.907,-7.657,-7.407,-7.157,-6.907,-6.657,-6.407,-6.157,-5.907,"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslz2="-5.657,-5.407,-5.157,-4.907,-4.657,-4.416,-4.187,-3.969,-3.761,-3.562,"
         polyslm2=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslt2=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslz3="-3.374,-3.194,-3.023,-2.860,-2.705,-2.557,-2.416,-2.282,-2.154,-2.033,"
         polyslm3=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslt3=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslz4="-1.917,-1.806,-1.701,-1.601,-1.506,-1.415,-1.329,-1.246,-1.168,-1.093,"
         polyslm4=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslt4=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslz5="-1.022,-0.955,-0.890,-0.829,-0.770,-0.714,-0.661,-0.611,-0.563,-0.517,"
         polyslm5=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslt5=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslz6="-0.473,-0.432,-0.392,-0.354,-0.318,-0.284,-0.252,-0.221,-0.191,-0.163,"
         polyslm6=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslt6=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslz7="-0.136,-0.111,-0.086,-0.063,-0.041,-0.020"
         polyslm7=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000"
         polyslt7=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000"
         ;;
      esac
      ;;
   esac
   #---------------------------------------------------------------------------------------#



   #----- Check whether to use SFILIN as restart or history. ------------------------------#
   if [[ ${runt} == "RESTORE" ]] && [[ ${forcehisto} -eq 1 ]]
   then
      runt="HISTORY"
      year=${yearh}
      month=${monthh}
      date=${dateh}
      time=${timeh}
      thissfilin=${fullygrown}
   elif [[ ${runt} == "RESTORE" ]]
   then
      case ${initmode} in
      5|7)
         #---------------------------------------------------------------------------------#
         #     Initialise ED2 with hdf5 files, which should be in directory restart.       #
         #---------------------------------------------------------------------------------#
         if [[ ! -s ${restart} ]]
         then
            echo " Directory restart does not exist!"
            echo " Change the variable restart at the beginning of the script"
            exit 44
         else
            runt="RESTORE"
            thissfilin=${restart}
         fi
         #---------------------------------------------------------------------------------#
         ;;

      1|2|3|6|8)
         #---------------------------------------------------------------------------------#
         #    Find the biometric files.  This has been checked in spawn_poly.sh so they    #
         # are correct.  Add a dummy name in case this is not supposed to be a biomass     #
         # initialisation run.                                                             #
         #---------------------------------------------------------------------------------#
         case ${biotype} in
         0)
            #----- isizepft controls everything, and iage is ignored. ---------------------#
            case ${isizepft} in
            0)
               #----- Frankeinstein's under storey. ---------------------------------------#
               thissfilin="${bioinit}/${polyiata}_default."
               ;;
            1)
               #----- No under storey. ----------------------------------------------------#
               thissfilin="${bioinit}/${polyiata}_nounder."
               ;;
            2)
               #----- ALS initialisation. -------------------------------------------------#
               thissfilin="${bioinit}/${polyiata}_alsinit."
               ;;
            *)
               #----- Invalid option. Stop the script. ------------------------------------#
               echo " Polygon:  ${polyname}"
               echo " IATA:     ${polyiata}"
               echo " ISIZEPFT: ${isizepft}"
               echo " INITMODE: ${initmode}"
               echo "This IATA cannot be initialised with these ISIZEPFT and INITMODE!"
               exit 57
               ;;
            esac
            #------------------------------------------------------------------------------#
            ;;

         1)
            #------------------------------------------------------------------------------#
            #    'isizepft' controls how many PFTs to use.                                 #
            #------------------------------------------------------------------------------#
            case ${isizepft} in 
            0|5)
               pftname="pft05"
               ;;
            2)
               pftname="pft02"
               ;;
            esac
            #------------------------------------------------------------------------------#

            #------------------------------------------------------------------------------#
            #     'iage' controls how many patches to use.                                 #
            #------------------------------------------------------------------------------#
            case ${iage} in
            1)
               agename="age01"
               ;;
            *)
               agename="age00"
               ;;
            esac
            #------------------------------------------------------------------------------#




            #------------------------------------------------------------------------------#
            #      Check whether the site has the PFT and age structure.                   #
            #------------------------------------------------------------------------------#
            case ${polyiata} in
            hvd|s77|fns|cau|and|par|tap|dcm)
               thissfilin="${bioinit}/${polyiata}_default."
               ;;
            cax|s67|s83|m34|gyf|pdg|rja|pnz|ban)
               thissfilin="${bioinit}/${polyiata}_${pftname}+${agename}."
               ;;
            *)
               echo " Polygon:  ${polyname}"
               echo " IATA:     ${polyiata}"
               echo " IAGE:     ${iage}"
               echo " ISIZEPFT: ${isizepft}"
               echo " INITMODE: ${initmode}"
               echo "This IATA cannot be initiealised with these settings!"
               exit 59
               ;;
            esac
            #------------------------------------------------------------------------------#
            ;;
         2)
            #------------------------------------------------------------------------------#
            #     ALS initialisation. ISIZEPFT has disturbance history information.        #
            #------------------------------------------------------------------------------#
            case ${polyiata} in
            l[0-5][0-3])
               thissfilin="${alsinit}/${polyiata}."
               ;;
            *)
               thissfilin="${alsinit}/${polyiata}_${isizepft}."
               ;;
            esac
            #------------------------------------------------------------------------------#
            ;;
         3)
            #------------------------------------------------------------------------------#
            #     ALS initialisation using intensity. ISIZEPFT has disturbance history     #
            # information.                                                                 #
            #------------------------------------------------------------------------------#
            thissfilin="${intinit}/${polyiata}_${isizepft}."
            #------------------------------------------------------------------------------#
            ;;
         4)
            #------------------------------------------------------------------------------#
            #     ALS initialisation using the lookup table. ISIZEPFT has disturbance      #
            # history information.                                                         #
            #------------------------------------------------------------------------------#
            thissfilin="${lutinit}/${polyiata}_${isizepft}."
            #------------------------------------------------------------------------------#
            ;;
         5)
            #----- isizepft controls actual or majestic initialisation. -------------------#
            case ${isizepft} in
            0)
               #----- Actual sampling. ----------------------------------------------------#
               thissfilin="${ebainit}/eba_actual_default."
               #---------------------------------------------------------------------------#
               ;;
            1)
               #----- Majestic sampling. --------------------------------------------------#
               thissfilin="${ebainit}/eba_majestic_default."
               #---------------------------------------------------------------------------#
               ;;
            2)
               #----- Actual sampling + land use (pasture/cropland/plantation). -----------#
               thissfilin="${ebainit}/eba_luactual_default."
               #---------------------------------------------------------------------------#
               ;;
            3)
               #----- Majestic sampling + land use (pasture/cropland/plantation). ---------#
               thissfilin="${ebainit}/eba_lumajestic_default."
               #---------------------------------------------------------------------------#
               ;;
            *)
               #----- Invalid option. Stop the script. ------------------------------------#
               echo " Polygon:  ${polyname}"
               echo " IATA:     ${polyiata}"
               echo " ISIZEPFT: ${isizepft}"
               echo " INITMODE: ${initmode}"
               echo "This IATA cannot be initiealised with these settings!"
               exit 57
               ;;
            esac
            #------------------------------------------------------------------------------#
            ;;
         esac
         #---------------------------------------------------------------------------------#
         ;;
      *)
         #----- No initial file needed; set the history file. -----------------------------#
         thissfilin=${here}/${polyname}/histo/${polyname}
         #---------------------------------------------------------------------------------#
         ;;
         #---------------------------------------------------------------------------------#
      esac
      #------------------------------------------------------------------------------------#
   else
      #----- Set the history file. --------------------------------------------------------#
      thissfilin=${here}/${polyname}/histo/${polyname}
      #------------------------------------------------------------------------------------#
   fi
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Set gfilout prefix according to ihrzrad.                                          #
   #---------------------------------------------------------------------------------------#
   case ${ihrzrad} in
   1) 
      gpref="gap"
      ;;
   2)
      gpref="pix"
      ;;
   *)
      gpref="dum"
      ;;
   esac
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Replace the flags in ED2IN.                                                       #
   #---------------------------------------------------------------------------------------#
   ED2IN="${here}/${polyname}/ED2IN"
   sed -i~ s@paththere@${here}@g                 ${ED2IN}
   sed -i~ s@myinpmain@${inpmain}@g              ${ED2IN}
   sed -i~ s@myyeara@${thisyeara}@g              ${ED2IN}
   sed -i~ s@mymontha@${montha}@g                ${ED2IN}
   sed -i~ s@mydatea@${datea}@g                  ${ED2IN}
   sed -i~ s@mytimea@${timea}@g                  ${ED2IN}
   sed -i~ s@myyearz@${thisyearz}@g              ${ED2IN}
   sed -i~ s@mymonthz@${monthz}@g                ${ED2IN}
   sed -i~ s@mydatez@${datez}@g                  ${ED2IN}
   sed -i~ s@mytimez@${timez}@g                  ${ED2IN}
   sed -i~ s@mydtlsm@${dtlsm}@g                  ${ED2IN}
   sed -i~ s@mymonyrstep@${monyrstep}@g          ${ED2IN}
   sed -i~ s@thispoly@${polyname}@g              ${ED2IN}
   sed -i~ s@plonflag@${polylon}@g               ${ED2IN}
   sed -i~ s@platflag@${polylat}@g               ${ED2IN}
   sed -i~ s@timehhhh@${time}@g                  ${ED2IN}
   sed -i~ s@datehhhh@${date}@g                  ${ED2IN}
   sed -i~ s@monthhhh@${month}@g                 ${ED2IN}
   sed -i~ s@yearhhhh@${year}@g                  ${ED2IN}
   sed -i~ s@myunitstate@${iunitstate}@g         ${ED2IN}
   sed -i~ s@myinitmode@${initmode}@g            ${ED2IN}
   sed -i~ s@mysfilin@${thissfilin}@g            ${ED2IN}
   sed -i~ s@mytrees@${pfts}@g                   ${ED2IN}
   sed -i~ s@mycrop@${crop}@g                    ${ED2IN}
   sed -i~ s@mypasture@${pasture}@g              ${ED2IN}
   sed -i~ s@myplantation@${plantation}@g        ${ED2IN}
   sed -i~ s@myiphen@${iphen}@g                  ${ED2IN}
   sed -i~ s@myallom@${iallom}@g                 ${ED2IN}
   sed -i~ s@myeconomics@${ieconomics}@g         ${ED2IN}
   sed -i~ s@mygrass@${igrass}@g                 ${ED2IN}
   sed -i~ s@myisoilflg@${polyisoil}@g           ${ED2IN}
   sed -i~ s@mynslcon@${polyntext}@g             ${ED2IN}
   sed -i~ s@myslxsand@${polysand}@g             ${ED2IN}
   sed -i~ s@myslxclay@${polyclay}@g             ${ED2IN}
   sed -i~ s@myslsoc@${polyslsoc}@g              ${ED2IN}
   sed -i~ s@myslph@${polyslph}@g                ${ED2IN}
   sed -i~ s@myslcec@${polyslcec}@g              ${ED2IN}
   sed -i~ s@mysldbd@${polysldbd}@g              ${ED2IN}
   sed -i~ s@myslhydro@${polyslhydro}@g          ${ED2IN}
   sed -i~ s@mysoilbc@${polysoilbc}@g            ${ED2IN}
   sed -i~ s@mysldrain@${polysldrain}@g          ${ED2IN}
   sed -i~ s@mysoilcol@${polycol}@g              ${ED2IN}
   sed -i~ s@mynzg@${polynzg}@g                  ${ED2IN}
   sed -i~ s@mymetdriverdb@${metdriverdb}@g      ${ED2IN}
   sed -i~ s@mymetcyc1@${metcyc1}@g              ${ED2IN}
   sed -i~ s@mymetcycf@${metcycf}@g              ${ED2IN}
   sed -i~ s@mytoler@${toler}@g                  ${ED2IN}
   sed -i~ s@RUNFLAG@${runt}@g                   ${ED2IN}
   sed -i~ s@myiphysiol@${iphysiol}@g            ${ED2IN}
   sed -i~ s@myvmfactc3@${vmfactc3}@g            ${ED2IN}
   sed -i~ s@myvmfactc4@${vmfactc4}@g            ${ED2IN}
   sed -i~ s@mymphototrc3@${mphototrc3}@g        ${ED2IN}
   sed -i~ s@mymphototec3@${mphototec3}@g        ${ED2IN}
   sed -i~ s@mymphotoc4@${mphotoc4}@g            ${ED2IN}
   sed -i~ s@mybphotoblc3@${bphotoblc3}@g        ${ED2IN}
   sed -i~ s@mybphotonlc3@${bphotonlc3}@g        ${ED2IN}
   sed -i~ s@mybphotoc4@${bphotoc4}@g            ${ED2IN}
   sed -i~ s@mykwgrass@${kwgrass}@g              ${ED2IN}
   sed -i~ s@mykwtree@${kwtree}@g                ${ED2IN}
   sed -i~ s@mygammac3@${gammac3}@g              ${ED2IN}
   sed -i~ s@mygammac4@${gammac4}@g              ${ED2IN}
   sed -i~ s@myd0grass@${d0grass}@g              ${ED2IN}
   sed -i~ s@myd0tree@${d0tree}@g                ${ED2IN}
   sed -i~ s@myalphac3@${alphac3}@g              ${ED2IN}
   sed -i~ s@myalphac4@${alphac4}@g              ${ED2IN}
   sed -i~ s@myklowco2@${klowco2}@g              ${ED2IN}
   sed -i~ s@mydecomp@${decomp}@g                ${ED2IN}
   sed -i~ s@myrrffact@${rrffact}@g              ${ED2IN}
   sed -i~ s@mygrowthresp@${growthresp}@g        ${ED2IN}
   sed -i~ s@mylwidthgrass@${lwidthgrass}@g      ${ED2IN}
   sed -i~ s@mylwidthbltree@${lwidthbltree}@g    ${ED2IN}
   sed -i~ s@mylwidthnltree@${lwidthnltree}@g    ${ED2IN}
   sed -i~ s@myq10c3@${q10c3}@g                  ${ED2IN}
   sed -i~ s@myq10c4@${q10c4}@g                  ${ED2IN}
   sed -i~ s@myh2olimit@${h2olimit}@g            ${ED2IN}
   sed -i~ s@mymortscheme@${imortscheme}@g       ${ED2IN}
   sed -i~ s@myddmortconst@${ddmortconst}@g      ${ED2IN}
   sed -i~ s@mycbrscheme@${cbrscheme}@g          ${ED2IN}
   sed -i~ s@mysfclyrm@${isfclyrm}@g             ${ED2IN}
   sed -i~ s@myicanturb@${icanturb}@g            ${ED2IN}
   sed -i~ s@myatmco2@${atmco2}@g                ${ED2IN}
   sed -i~ s@mythcrit@${thcrit}@g                ${ED2IN}
   sed -i~ s@mysmfire@${smfire}@g                ${ED2IN}
   sed -i~ s@myfire@${ifire}@g                   ${ED2IN}
   sed -i~ s@myfuel@${fireparm}@g                ${ED2IN}
   sed -i~ s@mymetavg@${imetavg}@g               ${ED2IN}
   sed -i~ s@mypercol@${ipercol}@g               ${ED2IN}
   sed -i~ s@myrunoff@${runoff}@g                ${ED2IN}
   sed -i~ s@mymetrad@${imetrad}@g               ${ED2IN}
   sed -i~ s@mybranch@${ibranch}@g               ${ED2IN}
   sed -i~ s@mycanrad@${icanrad}@g               ${ED2IN}
   sed -i~ s@myhrzrad@${ihrzrad}@g               ${ED2IN}
   sed -i~ s@mycrown@${crown}@g                  ${ED2IN}
   sed -i~ s@myltransvis@${ltransvis}@g          ${ED2IN}
   sed -i~ s@myltransnir@${ltransnir}@g          ${ED2IN}
   sed -i~ s@mylreflectvis@${lreflectvis}@g      ${ED2IN}
   sed -i~ s@mylreflectnir@${lreflectnir}@g      ${ED2IN}
   sed -i~ s@myorienttree@${orienttree}@g        ${ED2IN}
   sed -i~ s@myorientgrass@${orientgrass}@g      ${ED2IN}
   sed -i~ s@myclumptree@${clumptree}@g          ${ED2IN}
   sed -i~ s@myclumpgrass@${clumpgrass}@g        ${ED2IN}
   sed -i~ s@myigoutput@${igoutput}@g            ${ED2IN}
   sed -i~ s@mygpref@${gpref}@g                  ${ED2IN}
   sed -i~ s@myvegtdyn@${ivegtdyn}@g             ${ED2IN}
   sed -i~ s@myhydroscheme@${ihydro}@g           ${ED2IN}
   sed -i~ s@mystemresp@${istemresp}@g           ${ED2IN}
   sed -i~ s@mystomata@${istomata}@g             ${ED2IN}
   sed -i~ s@myplastic@${iplastic}@g             ${ED2IN}
   sed -i~ s@mycarbonmort@${icarbonmort}@g       ${ED2IN}
   sed -i~ s@myhydromort@${ihydromort}@g         ${ED2IN}
   sed -i~ s@mybigleaf@${ibigleaf}@g             ${ED2IN}
   sed -i~ s@myintegscheme@${integscheme}@g      ${ED2IN}
   sed -i~ s@mynsubeuler@${nsubeuler}@g          ${ED2IN}
   sed -i~ s@myrepro@${irepro}@g                 ${ED2IN}
   sed -i~ s@myubmin@${ubmin}@g                  ${ED2IN}
   sed -i~ s@myugbmin@${ugbmin}@g                ${ED2IN}
   sed -i~ s@myustmin@${ustmin}@g                ${ED2IN}
   sed -i~ s@mygamm@${gamm}@g                    ${ED2IN}
   sed -i~ s@mygamh@${gamh}@g                    ${ED2IN}
   sed -i~ s@mytprandtl@${tprandtl}@g            ${ED2IN}
   sed -i~ s@myribmax@${ribmax}@g                ${ED2IN}
   sed -i~ s@mygndvap@${igndvap}@g               ${ED2IN}
   sed -i~ s@mydtcensus@${dtcensus}@g            ${ED2IN}
   sed -i~ s@myyr1stcensus@${yr1stcensus}@g      ${ED2IN}
   sed -i~ s@mymon1stcensus@${mon1stcensus}@g    ${ED2IN}
   sed -i~ s@myminrecruitdbh@${minrecruitdbh}@g  ${ED2IN}
   sed -i~ s@mytreefall@${treefall}@g            ${ED2IN}
   sed -i~ s@mymaxpatch@${iage}@g                ${ED2IN}
   sed -i~ s@mymaxcohort@${imaxcohort}@g         ${ED2IN}
   sed -i~ s@myanthdisturb@${ianthdisturb}@g     ${ED2IN}
   sed -i~ s@myludatabase@${ludatabase}@g        ${ED2IN}
   sed -i~ s@myslscale@${slscale}@g              ${ED2IN}
   sed -i~ s@myslyrfirst@${slyrfirst}@g          ${ED2IN}
   sed -i~ s@myslnyrs@${slnyrs}@g                ${ED2IN}
   sed -i~ s@mylogging@${logging}@g              ${ED2IN}
   sed -i~ s@myprobharv@${probharv}@g            ${ED2IN}
   sed -i~ s@mydbhharv@${dbhharv}@g              ${ED2IN}
   sed -i~ s@mybioharv@${bioharv}@g              ${ED2IN}
   sed -i~ s@myskidarea@${skidarea}@g            ${ED2IN}
   sed -i~ s@myskiddbhthresh@${skiddbhthresh}@g  ${ED2IN}
   sed -i~ s@myskidsmall@${skidsmall}@g          ${ED2IN}
   sed -i~ s@myskidlarge@${skidlarge}@g          ${ED2IN}
   sed -i~ s@myfellingsmall@${fellingsmall}@g    ${ED2IN}
   sed -i~ s@myseedharv@${seedharv}@g            ${ED2IN}
   sed -i~ s@mystorharv@${storharv}@g            ${ED2IN}
   sed -i~ s@myleafharv@${leafharv}@g            ${ED2IN}
   #---------------------------------------------------------------------------------------#

   #------ Soil variables. ----------------------------------------------------------------#
   sed -i~ s@myslz1@"${polyslz1}"@g           ${ED2IN}
   sed -i~ s@myslz2@"${polyslz2}"@g           ${ED2IN}
   sed -i~ s@myslz3@"${polyslz3}"@g           ${ED2IN}
   sed -i~ s@myslz4@"${polyslz4}"@g           ${ED2IN}
   sed -i~ s@myslz5@"${polyslz5}"@g           ${ED2IN}
   sed -i~ s@myslz6@"${polyslz6}"@g           ${ED2IN}
   sed -i~ s@myslz7@"${polyslz7}"@g           ${ED2IN}
   sed -i~ s@myslmstr1@"${polyslm1}"@g        ${ED2IN}
   sed -i~ s@myslmstr2@"${polyslm2}"@g        ${ED2IN}
   sed -i~ s@myslmstr3@"${polyslm3}"@g        ${ED2IN}
   sed -i~ s@myslmstr4@"${polyslm4}"@g        ${ED2IN}
   sed -i~ s@myslmstr5@"${polyslm5}"@g        ${ED2IN}
   sed -i~ s@myslmstr6@"${polyslm6}"@g        ${ED2IN}
   sed -i~ s@myslmstr7@"${polyslm7}"@g        ${ED2IN}
   sed -i~ s@mystgoff1@"${polyslt1}"@g        ${ED2IN}
   sed -i~ s@mystgoff2@"${polyslt2}"@g        ${ED2IN}
   sed -i~ s@mystgoff3@"${polyslt3}"@g        ${ED2IN}
   sed -i~ s@mystgoff4@"${polyslt4}"@g        ${ED2IN}
   sed -i~ s@mystgoff5@"${polyslt5}"@g        ${ED2IN}
   sed -i~ s@mystgoff6@"${polyslt6}"@g        ${ED2IN}
   sed -i~ s@mystgoff7@"${polyslt7}"@g        ${ED2IN}
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #   In case this is a multithreaded run, copy executables to each directory.            #
   #---------------------------------------------------------------------------------------#
   case ${n_cpt} in
   1)
      exec_sub="${here}/${polyname}/${execname}"
      cp ${exec_full} ${exec_sub}
      ;;
   *)
      exec_sub=${exec_full}
      ;;
   esac
   #---------------------------------------------------------------------------------------#



   #----- Script names. -------------------------------------------------------------------#
   tempserial="${here}/Template/callserial.sh"
   callserial="${here}/${polyname}/callserial.sh"
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #    Set script callserial.sh.  Depending on the settings (single multi-task jobs or    #
   # multiple single-task jobs), we must update the header.                                #
   #---------------------------------------------------------------------------------------#
   case "${cluster}" in
   SDUMONT)
      #------------------------------------------------------------------------------------#
      #       Single multi task job, callserial.sh doesn't need header updates.            #
      #------------------------------------------------------------------------------------#
      /bin/rm -f ${callserial}
      /bin/cp -f ${tempserial} ${callserial}
      #------------------------------------------------------------------------------------#
      ;;
   CANNON)
      #------------------------------------------------------------------------------------#
      #       Update header with sbatch instructions for multiple single task jobs.        #
      #------------------------------------------------------------------------------------#

      #----- Output files. ----------------------------------------------------------------#
      obatch="${here}/${polyname}/serial_slm.out"
      ebatch="${here}/${polyname}/serial_slm.err"
      #------------------------------------------------------------------------------------#

      #----- Initialise the callserial.sh file. -------------------------------------------#
      rm -f     ${callserial}
      touch     ${callserial}
      chmod u+x ${callserial}
      echo "#!/bin/bash" >> ${callserial}
      echo "#SBATCH --ntasks=1                  # Number of tasks"         >> ${callserial}
      echo "#SBATCH --cpus-per-task=${n_cpt}    # Number of CPUs per task" >> ${callserial}
      echo "#SBATCH --cores-per-socket=${n_cpt} # Min. # of cores/socket"  >> ${callserial}
      echo "#SBATCH --partition=${global_queue} # Queue that will run job" >> ${callserial}
      echo "#SBATCH --job-name=${taskname}      # Job name"                >> ${callserial}
      echo "#SBATCH --mem-per-cpu=${cpu_memory} # Memory per CPU"          >> ${callserial}
      echo "#SBATCH --time=${runtime}           # Time for job"            >> ${callserial}
      echo "#SBATCH --output=${obatch}          # Standard output path"    >> ${callserial}
      echo "#SBATCH --error=${ebatch}           # Standard error path"     >> ${callserial}
      echo "#SBATCH --chdir=${here}/${polyname} # Main directory"          >> ${callserial}
      echo ""                                                              >> ${callserial}
      echo "#--- Initial settings."                                        >> ${callserial}
      echo "here=\"${here}\"                    # Main path"               >> ${callserial}
      echo "exec=\"${exec_full}\"               # Executable"              >> ${callserial}
      echo ""                                                              >> ${callserial}
      echo "#--- Print information about this job."                        >> ${callserial}
      echo "echo \"\""                                                     >> ${callserial}
      echo "echo \"\""                                                     >> ${callserial}
      echo "echo \"----- Summary of current job -----------------------\"" >> ${callserial}
      echo "echo \" CPUs per task:   \${SLURM_CPUS_PER_TASK}\""            >> ${callserial}
      echo "echo \" Job name:        \${SLURM_JOB_NAME}\""                 >> ${callserial}
      echo "echo \" Job ID:          \${SLURM_JOB_ID}\""                   >> ${callserial}
      echo "echo \" Partition:       \${SLURM_JOB_PARTITION}\""            >> ${callserial}
      echo "echo \" Number of nodes: \${SLURM_NNODES}\""                   >> ${callserial}
      echo "echo \" Number of tasks: \${SLURM_NTASKS}\""                   >> ${callserial}
      echo "echo \" Memory per CPU:  \${SLURM_MEM_PER_CPU}\""              >> ${callserial}
      echo "echo \" Memory per node: \${SLURM_MEM_PER_NODE}\""             >> ${callserial}
      echo "echo \" Node list:       \${SLURM_JOB_NODELIST}\""             >> ${callserial}
      echo "echo \" Time limit:      \${SLURM_TIMELIMIT}\""                >> ${callserial}
      echo "echo \" Std. Output:     \${SLURM_STDOUTMODE}\""               >> ${callserial}
      echo "echo \" Std. Error:      \${SLURM_STDERRMODE}\""               >> ${callserial}
      echo "echo \"----------------------------------------------------\"" >> ${callserial}
      echo "echo \"\""                                                     >> ${callserial}
      echo "echo \"\""                                                     >> ${callserial}
      echo ""                                                              >> ${callserial}
      echo "echo \"----- Global settings for this simulation. ---\""       >> ${callserial}
      echo "echo \" Main path:       \${here}\""                           >> ${callserial}
      echo "echo \" Executable:      \${exec}\""                           >> ${callserial}
      echo "echo \"----------------------------------------------------\"" >> ${callserial}
      echo "echo \"\""                                                     >> ${callserial}
      echo "echo \"\""                                                     >> ${callserial}
      echo ""                                                              >> ${callserial}
      echo ""                                                              >> ${callserial}
      echo "#--- Define home in case home is not set"                      >> ${callserial}
      echo "if [[ \"x\${HOME}\" == \"x\" ]]"                               >> ${callserial}
      echo "then"                                                          >> ${callserial}
      echo "   export HOME=\$(echo ~)"                                     >> ${callserial}
      echo "fi"                                                            >> ${callserial}
      echo ""                                                              >> ${callserial}
      echo "#--- Load modules and settings."                               >> ${callserial}
      echo ". \${HOME}/.bashrc ${optsrc}"                                  >> ${callserial}
      echo ""                                                              >> ${callserial}
      echo "#--- Get plenty of memory."                                    >> ${callserial}
      echo "ulimit -s unlimited"                                           >> ${callserial}
      echo "ulimit -u unlimited"                                           >> ${callserial}
      echo ""                                                              >> ${callserial}
      echo "#--- Set OpenMP parameters"                                    >> ${callserial}
      echo "if [[ \"\${SLURM_CPUS_PER_TASK}\" == \"\" ]]"                  >> ${callserial}
      echo "then"                                                          >> ${callserial}
      echo "   export OMP_NUM_THREADS=1"                                   >> ${callserial}
      echo "else"                                                          >> ${callserial}
      echo "   export OMP_NUM_THREADS=\${SLURM_CPUS_PER_TASK}"             >> ${callserial}
      echo "fi"                                                            >> ${callserial}
      echo ""                                                              >> ${callserial}
      #------------------------------------------------------------------------------------#



      #----- Append the template callserial.sh, skipping the first line. ------------------#
      tail -n +2 ${tempserial} >> ${callserial}
      #------------------------------------------------------------------------------------#


      #----- Make sure the script waits until all tasks are completed... ------------------#
      echo ""                                                              >> ${callserial}
      echo ""                                                              >> ${callserial}
      echo "#----- Ensure that jobs complete before terminating script."   >> ${callserial}
      echo "wait"                                                          >> ${callserial}
      echo ""                                                              >> ${callserial}
      echo "#----- Report efficiency of this job."                         >> ${callserial}
      echo "seff \${SLURM_JOBID}"                                          >> ${callserial}
      echo ""                                                              >> ${callserial}
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      ;;
   esac
   #---------------------------------------------------------------------------------------#


   #----- Substitute placeholders with specific data. -------------------------------------#
   sed -i s@thisroot@${here}@g          ${callserial}
   sed -i s@thispoly@${polyname}@g      ${callserial}
   sed -i s@myname@${moi}@g             ${callserial}
   sed -i s@myexec@${exec_sub}@g        ${callserial}
   sed -i s@mypackdata@${packdatasrc}@g ${callserial}
   sed -i s@myscenario@${iscenario}@g   ${callserial}
   sed -i s@myscenmain@${scentype}@g    ${callserial}
   sed -i s@mycopy@${copy2scratch}@g    ${callserial}
   sed -i s@mycpus@${n_cpt}@g           ${callserial}
   sed -i s@myoptsrc@${optsrc}@g        ${callserial}
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     We will not even consider the files that have gone extinct.                       #
   #---------------------------------------------------------------------------------------#
   case "${runt}" in
   "THE_END")
      echo "Polygon has reached the end.  No need to re-submit it."
      ;;
   "STSTATE")
      echo "Polygon has reached steady state.  No need to re-submit it."
      ;;
   "EXTINCT")
      echo "Polygon population has gone extinct.  No need to re-submit it."
      ;;
   "CRASHED"|"HYDFAIL"|"METMISS"|"SIGSEGV"|"BAD_MET"|"STOPPED")
      #----- Decide whether to skip this job submission or every job submission. ----------#
      if ${on_the_fly}
      then
         echo "Polygon has serious errors.  Script will not submit this job."
      else
         echo "Polygon has serious errors.  Script will not submit any job this time."
         submit=false
      fi
      #------------------------------------------------------------------------------------#
      ;;

   "RESTORE"|"HISTORY")

      #------------------------------------------------------------------------------------#
      #      Update job count.                                                             #
      #------------------------------------------------------------------------------------#
      echo "  Polygon scheduled for submission."
      let n_submit=${n_submit}+1
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #    Set command to the main script. Again, this depends on the type of job          #
      # submission.                                                                        #
      #------------------------------------------------------------------------------------#
      case "${cluster}" in
      SDUMONT)
         #---------------------------------------------------------------------------------#
         #     Append task in the single multi-task job script.                            #
         #---------------------------------------------------------------------------------#
         srun="srun --nodes=1 --ntasks=1 --cpu_bind=cores"
         srun="${srun} --cpus-per-task=\${SLURM_CPUS_PER_TASK}"
         srun="${srun} --mem-per-cpu=\${SLURM_MEM_PER_CPU}"
         srun="${srun} --job-name=${polyname}"
         srun="${srun} --chdir=\${here}/${polyname}"
         srun="${srun} --output=\${here}/${polyname}/serial_slm.out"
         srun="${srun} --error=\${here}/${polyname}/serial_slm.err"
         echo "${srun} \${here}/${polyname}/callserial.sh &" >> ${sbatch}
         echo "sleep ${dttask}" >> ${sbatch}
         #---------------------------------------------------------------------------------#
         ;;
      CANNON)
         #---------------------------------------------------------------------------------#
         #      Make sure we don't exceed the maximum number of jobs.                      #
         #---------------------------------------------------------------------------------#
         if [[ ${n_submit} -gt ${n_tasks_max} ]] && ${on_the_fly}
         then
            echo " Number of jobs to submit: ${n_submit}"
            echo " Maximum number of tasks in queue ${global_queue}: ${n_tasks_max}"
            echo " Reduce the number of simulations or try another queue..."
            echo " Appending remaining jobs to $(basename ${sbatch}) for later submission."
            on_the_fly=false
         fi
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #       If on the fly, submit the job now.  Otherwise, append job to sbatch.      #
         #---------------------------------------------------------------------------------#
         if ${on_the_fly}
         then
            #----- Submit job on the fly. -------------------------------------------------#
            sbatch ${callserial}
            #------------------------------------------------------------------------------#
         else
            #------------------------------------------------------------------------------#
            #     Append submission to the script to submit the multiple single-task jobs. #
            #------------------------------------------------------------------------------#
            echo "echo \" + ${ffout}/${ntasks}. Submit job: ${polyname}.\""    >> ${sbatch}
            echo "sbatch ${callserial}"                                        >> ${sbatch}
            echo "sleep  ${dttask}"                                            >> ${sbatch}
            #------------------------------------------------------------------------------#
         fi
         #---------------------------------------------------------------------------------#

         ;;
      esac
      #------------------------------------------------------------------------------------#
      ;;
   *)
      #------------------------------------------------------------------------------------#
      #    Unknown state, either skip this job submission, or every job submission.        #
      #------------------------------------------------------------------------------------#
      if ${on_the_fly}
      then
         echo "Unknown polygon state (${runt})! Script will not submit this job."
      else
         echo "Unknown polygon state (${runt})! Script will not submit any job this time."
         submit=false
      fi
      #------------------------------------------------------------------------------------#
      ;;
   esac
   #---------------------------------------------------------------------------------------#
done
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#     Update the main submission script to reflect the correct number of tasks (done only  #
# when running a single multi-task job).                                                   #
#------------------------------------------------------------------------------------------#
case "${cluster}" in
SDUMONT)
   #----- Update the number of tasks in batch script. -------------------------------------#
   sed -i~ s@myntasks@${n_submit}@g ${sbatch}
   #---------------------------------------------------------------------------------------#
   ;;
esac
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#      Make sure job list doesn't request too many nodes.                                  #
#------------------------------------------------------------------------------------------#
if ! ${on_the_fly}
then
   if [[ ${n_submit} -gt ${n_tasks_max} ]]
   then
      echo " Number of jobs to submit: ${n_submit}"
      echo " Maximum number of tasks in queue ${global_queue}: ${n_tasks_max}"
      echo " Reduce the number of simulations or try another queue..."
      exit 99
   elif ${submit}
   then
      echo " Submitting jobs."
      #------- How we submit depends on the cluster. --------------------------------------#
      case "${cluster}" in
      SDUMONT)
         #---------------------------------------------------------------------------------#
         #     Single multi-task job.  Submit the main script using sbatch.                #
         #---------------------------------------------------------------------------------#
         sbatch ${sbatch}
         #---------------------------------------------------------------------------------#
         ;;
      CANNON)
         #---------------------------------------------------------------------------------#
         #     Multiple single-task jobs.  Run the script, the sbatch commands are in      #
         # there.                                                                          #
         #---------------------------------------------------------------------------------#
         ${sbatch}
         #---------------------------------------------------------------------------------#
         ;;
      esac
      #------------------------------------------------------------------------------------#
   else
      #------- Provide instructions to the user for a later submission. -------------------#
      case "${cluster}" in
      SDUMONT)
         #----- Instruct the user to use sbatch. ------------------------------------------#
         echo "-------------------------------------------------------------------------"
         echo " To submit the simulations, you must use sbatch. "
         echo " Copy and paste the following command in the terminal. "
         echo " "
         echo "       sbatch ${sbatch}"
         echo " "
         #---------------------------------------------------------------------------------#
         ;;
      CANNON)
         #----- Instruct the user NOT to use sbatch. --------------------------------------#
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
         #---------------------------------------------------------------------------------#
         ;;
      esac
      #------------------------------------------------------------------------------------#
   fi
   #---------------------------------------------------------------------------------------#
fi
#------------------------------------------------------------------------------------------#
