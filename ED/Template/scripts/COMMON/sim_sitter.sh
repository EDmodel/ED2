#!/bin/bash



#------------------------------------------------------------------------------------------#
#     Which scripts to run.                                                                #
#                                                                                          #
#   - read_monthly.r - This reads the monthly mean files (results can then be used for     #
#                      plot_monthly.r, plot_yearly.r, and others, but it doesn't plot any- #
#                      thing.)                                                             #
#   - yearly_ascii.r - This creates three ascii (csv) files with annual averages of        #
#                      various variables.  It doesn't have all possible variables as it is #
#                      intended to simplify the output for learning purposes.              #
#   - plot_monthly.r - This creates several plots based on the monthly mean output.        #
#   - plot_yearly.r  - This creates plots with year time series.                           #
#   - plot_ycomp.r   - This creates yearly comparisons based on the monthly mean output.   #
#   - plot_povray.r  - This creates yearly plots of the polygon using POV-Ray.             #
#   - plot_rk4.r     - This creates plots from the detailed output for Runge-Kutta.        #
#                      (patch-level only).                                                 #
#   - plot_photo.r   - This creates plots from the detailed output for Farquhar-Leuning.   #
#   - plot_rk4pc.r   - This creates plots from the detailed output for Runge-Kutta.        #
#                      (patch- and cohort-level).                                          #
#   - plot_budget.r  - This creates plots from the detailed budget for Runge-Kutta.        #
#                      (patch-level only).                                                 #
#   - plot_eval_ed.r - This creates plots comparing model with eddy flux observations.     #
#   - plot_census.r  - This creates plots comparing model with biometric data.             #
#   - whichrun.r     - This checks the run status.                                         #
#                                                                                          #
#   The following scripts should work too, but I haven't tested them.                      #
#   - plot_daily.r   - This creates plots from the daily mean output.                      #
#   - plot_fast.r    - This creates plots from the analysis files.                         #
#   - patchprops.r   - This creates simple plots showing the patch structure.              #
#   - reject_ed.r    - This tracks the number of steps that were rejected, and what caused #
#                      the step to be rejected.                                            #
#------------------------------------------------------------------------------------------#


#----- Main settings.  Do look at run_sitter.sh and epost.sh for additional settings. -----#
here=$(pwd)                        # Current path.
there="/path/to/permanent/storage" # Permanent storage path (leave empty for no transfer)
email="myself\@myserver.com"       # Your e-mail. Put a backslash before the @
sbatch=$(which sbatch)             # SLURM command to submit job.
sitter_queue="myqueue"             # Queue to run run_sitter.sh
sitter_runtime="mytime"            # Run time request for run_sitter.sh
sitter_memory=0                    # Requested memory (Mb) for run_sitter.sh
lhisto_queue="myqueue"             # Queue to run last_histo.sh
lhisto_runtime="mytime"            # Run time request for last_histo.sh
lhisto_memory=0                    # Requested memory (Mb) for last_histo.sh
transfer_full=false                # Flag to decide between full and partial transfer
epost_queue="myqueue"              # Queue to run epost.sh
epost_runtime="mytime"             # Run time request for epost.sh
epost_memory=0                     # Requested memory (Mb) for epost.sh
epost_reserve=""                   # Reservation flag (leave blank unless you have 
                                   #    reservation privileges).
rscript="nothing.r"                # Which script to run with epost.sh.  See options above. 
                                   #    Multiple scripts are allowed, put spaces between
                                   #    them. (e.g. rscript="plot_monthly.r plot_fast.r")
overcommit=false                   # Ignore maximum number of tasks when submitting
                                   #    epost tasks?
frqemail=43200                     # How often to send emails on simulation status?
delay1st_min=20                    # Time (in minutes) to wait before the first check
                                   #    (needed in case the run_sitter script is submitted
                                   #    before the simulation starts.  In case the
                                   #    simulation is already running, it is fine to set
                                   #    this to zero).
wait_minutes=60                    # Waiting time before checking run again (in minutes)
frqpost=3                          # How often to run post-processing and file management
                                   #    This number is in iterations.  Zero means never.
frqtouch=0                         # How often to touch executable, scripts, and Template?
                                   #    This number is in iterations.  Zero means never.
checkhourly="n"                    # Check hourly files.
checkstatus="y"                    # Check status before compressing
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
#                  THINK TWICE BEFORE CHANGING ANYTHING BEYOND THIS POINT.                 #
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


#----- Make sure e-mail and queue are pre-defined. ----------------------------------------#
if [[ "${email}"           == "myself\@myserver.com"       ]] ||
   [[ "${there}"           == "/path/to/permanent/storage" ]] ||
   [[ "${rscript}"         == "nothing.r"                  ]] ||
   [[ "${sitter_queue}"    == "myqueue"                    ]] ||
   [[ "${sitter_runtime}"  == "mytime"                     ]] ||
   [[  ${sitter_memory}    -eq 0                           ]] ||
   [[ "${lhisto_queue}"    ==  "myqueue"                   ]] ||
   [[ "${lhisto_runtime}"  ==  "mytime"                    ]] ||
   [[  ${lhisto_memory}    -eq 0                           ]] ||
   [[ "${epost_queue}"     ==  "myqueue"                   ]] ||
   [[ "${epost_runtime}"   ==  "mytime"                    ]] ||
   [[  ${epost_memory}     -eq 0                           ]]
then
   echo "---------------------------------------------------------------------------------"
   echo "    The following variables must be set.  In case any of them have dummy values, "
   echo " check your script sim_sitter.sh."
   echo " "
   echo " email          = ${email}"
   echo " there          = ${there}"
   echo " rscript        = ${rscript}"
   echo " sitter_queue   = ${sitter_queue}"
   echo " epost_queue    = ${epost_queue}"
   echo " lhisto_queue   = ${lhisto_queue}"
   echo " sitter_queue   = ${sitter_queue}"
   echo " sitter_runtime = ${sitter_runtime}"
   echo " sitter_memory  = ${sitter_memory}"
   echo " lhisto_queue   = ${lhisto_queue}"
   echo " lhisto_runtime = ${lhisto_runtime}"
   echo " lhisto_memory  = ${lhisto_memory}"
   echo " epost_queue    = ${epost_queue}"
   echo " epost_runtime  = ${epost_runtime}"
   echo " epost_memory   = ${epost_memory}"
   echo " "
   echo "---------------------------------------------------------------------------------"
   exit 99
fi
#------------------------------------------------------------------------------------------#



#----- Find out which platform we are using. ----------------------------------------------#
if [ "${1}" == "" ]
then
   #------ No platform provided.  Try to guess, and if failed, then prompts the user. -----#
   host=$(hostname -s)
   case ${host} in
      rclogin*|holy*|moorcroft*|rcnx*|sdumont*)
         #----- Use SLURM scripts. --------------------------------------------------------#
         platform="SLURM"
         #---------------------------------------------------------------------------------#
         ;;
      au*|ha*|sun-master|cmm*)
         #----- Use PBS scripts. ----------------------------------------------------------#
         platform="PBS"
         #---------------------------------------------------------------------------------#
         ;;
      *)
         #----- Host name is not one of the known ones. -----------------------------------#
         echo -n "Failed guessing platform from node name.  Please type the name:   "
         read platform
         #---------------------------------------------------------------------------------#
         ;;
   esac
   #---------------------------------------------------------------------------------------#
else
   #------ Platform is provided as argument. ----------------------------------------------#
   platform=${1}
   #---------------------------------------------------------------------------------------#
fi
#------------------------------------------------------------------------------------------#




#------ List of scripts. ------------------------------------------------------------------#
run_sitter="${here}/run_sitter.sh"
epost="${here}/epost.sh"
transfer="${here}/transfer.sh"
orig_histo="${here}/scripts/COMMON/last_histo.sh"
last_histo="${here}/last_histo.sh"
#------------------------------------------------------------------------------------------#



#----- Set job prefix. --------------------------------------------------------------------#
desc=$(basename ${here})
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#    In case this is a SLURM cluster, append header to last_histo.sh.                      #
#------------------------------------------------------------------------------------------#
case "${platform}" in
SLURM)
   #----- Reset last_histo. ---------------------------------------------------------------#
   echo " + Reset $(basename ${last_histo})"
   /bin/rm -f ${last_histo}
   touch ${last_histo}
   chmod u+x ${last_histo}
   #---------------------------------------------------------------------------------------#


   #----- Set some job instructions. ------------------------------------------------------#
   lhisto_task="${desc}-last_histo.sh"
   lhisto_sto="${here}/out_last_histo.out"
   lhisto_ste="${here}/out_last_histo.err"
   #---------------------------------------------------------------------------------------#



   #----- Add header. ---------------------------------------------------------------------#
   echo "#!/bin/bash"                                                     >> ${last_histo}
   echo "#SBATCH --ntasks=1                      # Number of tasks"       >> ${last_histo}
   echo "#SBATCH --cpus-per-task=1               # CPUs per task"         >> ${last_histo}
   echo "#SBATCH --partition=${lhisto_queue}     # Job partition"         >> ${last_histo}
   echo "#SBATCH --job-name=${lhisto_task}       # Task name"             >> ${last_histo}
   echo "#SBATCH --mem-per-cpu=${lhisto_memory}  # Memory per CPU"        >> ${last_histo}
   echo "#SBATCH --time=${lhisto_runtime}        # Time for job"          >> ${last_histo}
   echo "#SBATCH --output=${lhisto_sto}          # Standard output path"  >> ${last_histo}
   echo "#SBATCH --error=${lhisto_ste}           # Standard error path"   >> ${last_histo}
   echo "#SBATCH --chdir=${here}                 # Main directory"        >> ${last_histo}
   echo ""                                                                >> ${last_histo}
   echo "#--- Initial settings."                                          >> ${last_histo}
   echo "here=\"${here}\"                     # Main path"                >> ${last_histo}
   echo ""                                                                >> ${last_histo}
   echo "#--- Print information about this job."                          >> ${last_histo}
   echo "echo \"\""                                                       >> ${last_histo}
   echo "echo \"\""                                                       >> ${last_histo}
   echo "echo \"----- Summary of current job -------------------------\"" >> ${last_histo}
   echo "echo \" CPUs per task:   \${SLURM_CPUS_PER_TASK}\""              >> ${last_histo}
   echo "echo \" Job name:        \${SLURM_JOB_NAME}\""                   >> ${last_histo}
   echo "echo \" Job ID:          \${SLURM_JOB_ID}\""                     >> ${last_histo}
   echo "echo \" Queue:           \${SLURM_JOB_PARTITION}\""              >> ${last_histo}
   echo "echo \" Number of nodes: \${SLURM_NNODES}\""                     >> ${last_histo}
   echo "echo \" Number of tasks: \${SLURM_NTASKS}\""                     >> ${last_histo}
   echo "echo \" Memory per CPU:  \${SLURM_MEM_PER_CPU}\""                >> ${last_histo}
   echo "echo \" Memory per node: \${SLURM_MEM_PER_NODE}\""               >> ${last_histo}
   echo "echo \" Node list:       \${SLURM_JOB_NODELIST}\""               >> ${last_histo}
   echo "echo \" Time limit:      \${SLURM_TIMELIMIT}\""                  >> ${last_histo}
   echo "echo \" Std. Output:     \${SLURM_STDOUTMODE}\""                 >> ${last_histo}
   echo "echo \" Std. Error:      \${SLURM_STDERRMODE}\""                 >> ${last_histo}
   echo "echo \"------------------------------------------------------\"" >> ${last_histo}
   echo "echo \"\""                                                       >> ${last_histo}
   echo "echo \"\""                                                       >> ${last_histo}
   echo "echo \"\""                                                       >> ${last_histo}
   echo "echo \"\""                                                       >> ${last_histo}
   echo ""                                                                >> ${last_histo}
   echo ""                                                                >> ${last_histo}
   echo "#--- Define home in case home is not set"                        >> ${last_histo}
   echo "if [[ \"x\${HOME}\" == \"x\" ]]"                                 >> ${last_histo}
   echo "then"                                                            >> ${last_histo}
   echo "   export HOME=\$(echo ~)"                                       >> ${last_histo}
   echo "fi"                                                              >> ${last_histo}
   echo ""                                                                >> ${last_histo}
   echo "#--- Load modules and settings."                                 >> ${last_histo}
   echo ". \${HOME}/.bashrc"                                              >> ${last_histo}
   #---------------------------------------------------------------------------------------#



   #----- Append the template last_histo.sh, skipping the first two lines. ----------------#
   tail -n +3 ${orig_histo} >> ${last_histo}
   #---------------------------------------------------------------------------------------#

   ;;
esac
#------------------------------------------------------------------------------------------#






#------------------------------------------------------------------------------------------#
#     Make substitutions.                                                                  #
#     IMPORTANT: make sure "there" substitutions precede "here", otherwise the script will #
#                not work as intended.                                                     #
#------------------------------------------------------------------------------------------#
sed -i~ s@"here=\"\""@"here=\"${here}\""@g                            ${run_sitter}
sed -i~ s@"recipient=\"\""@"recipient=\"${email}\""@g                 ${run_sitter}
sed -i~ s@"frqemail=\"\""@"frqemail=${frqemail}"@g                    ${run_sitter}
sed -i~ s@"delay1st_min=\"\""@"delay1st_min=${delay1st_min}"@g        ${run_sitter}
sed -i~ s@"wait_minutes=\"\""@"wait_minutes=${wait_minutes}"@g        ${run_sitter}
sed -i~ s@"frqpost=\"\""@"frqpost=${frqpost}"@g                       ${run_sitter}
sed -i~ s@"frqtouch=\"\""@"frqtouch=${frqtouch}"@g                    ${run_sitter}
sed -i~ s@"here=\"\""@"here=\"${here}\""@g                            ${epost}
sed -i~ s@"global_queue=\"\""@"global_queue=\"${epost_queue}\""@g     ${epost}
sed -i~ s@"rscript=\"\""@"rscript=\"${rscript}\""@g                   ${epost}
sed -i~ s@"submit=false"@"submit=true"@g                              ${epost}
sed -i~ s@"reservation=\"\""@"reservation=\"${epost_reserve}\""@g     ${epost}
sed -i~ s@"sim_memory=0"@"sim_memory=${epost_memory}"@g               ${epost}
sed -i~ s@"runtime=\"00:00:00\""@"runtime=\"${epost_runtime}\""@g     ${epost}
sed -i~ s@"overcommit=\"\""@"overcommit=${overcommit}"@g              ${epost}
sed -i~ s@"there=\"\""@"there=\"${there}\""@g                         ${transfer}
sed -i~ s@"here=\"\""@"here=\"${here}\""@g                            ${transfer}
sed -i~ s@"full_transfer=boolean"@"full_transfer=${transfer_full}"@g  ${transfer}
sed -i~ s@"here=\"\""@"here=\"${here}\""@g                            ${last_histo}
sed -i~ s@"checkhourly=\"\""@"checkhourly=\"${checkhourly}\""@g       ${last_histo}
sed -i~ s@"checkstatus=\"\""@"checkstatus=\"${checkstatus}\""@g       ${last_histo}
#------------------------------------------------------------------------------------------#


#----- Substitute paths in sit_utils R scripts. -------------------------------------------#
sed -i~ s@"mypath"@"${here}"@g      ${here}/sit_utils/plot.region.r
sed -i~ s@"mypath"@"${here}"@g      ${here}/sit_utils/plot.status.r
#------------------------------------------------------------------------------------------#


#----- Job preferences. -------------------------------------------------------------------#
sitter_joblog="${here}/out_sitter.out"
sitter_jobpref=$(basename ${here})
sitter_jobname="${desc}-sitter"
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#    Check whether to submit the job.                                                      #
#------------------------------------------------------------------------------------------#
if [[ -s "${here}/run_sitter.lock" ]] || [[ -s "${here}/transfer.lock" ]]
then
   echo " Lock file found for run_sitter.sh and/or transfer.lock."
   echo " The scripts may be already running."
   echo " Submit the script anyway [y|N]?"
   read goahead
else
   goahead="y"
fi
goahead="$(echo ${goahead} | tr '[:upper:]' '[:lower:]')"
#------------------------------------------------------------------------------------------#


#----- Submit run_sitter.sh in batch mode. ------------------------------------------------#
comm="${sbatch} -p ${sitter_queue} --mem-per-cpu=${sitter_memory} -t ${sitter_runtime}"
comm="${comm} -o ${sitter_joblog} -J ${sitter_jobname} -n 1 "
comm="${comm} --wrap=\"${here}/run_sitter.sh\""
case ${goahead} in
y|yes)
   /bin/rm -f ${here}/run_sitter.lock
   /bin/rm -f ${here}/transfer.lock
   /bin/rm -f ${here}/sit_utils/mycheck.txt
   /bin/rm -f ${here}/sit_utils/lastcheck.txt
   echo ${comm}
   ${comm}
   ;;
*)
   bann="Script run_sitter was not submitted.  In case you want to submit manually"
   bann="${bann} this is the command:"
   echo ${bann}
   echo " "
   echo ${comm}
   ;;
esac
#------------------------------------------------------------------------------------------#
