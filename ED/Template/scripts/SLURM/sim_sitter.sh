#!/bin/bash

#----- Main settings.  Do look at run_sitter.sh and epost.sh for additional settings. -----#
here=$(pwd)                        # Current path.
there="/path/to/permanent/storage" # Permanent storage path (leave empty for no transfer)
email="myself\@myserver.com"       # Your e-mail. Put a backslash before the @
queue="myqueue"                    # Queue to run run_sitter.sh and epost.sh
runtime="infinite"                 # Run time request
memory=2048                        # Requested memory (Mb)
sbatch=$(which sbatch)             # SLURM command to submit job.
rscript="read_monthly.r"           # Which script to run with epost.sh.  See options below. 
                                   #    Multiple scripts are allowed, put spaces between
                                   #    them. (e.g. rscript="plot_monthly.r plot_fast.r")
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


#----- Make sure e-mail and queue are pre-defined. ----------------------------------------#
if [[ "x${email}" == "x" ]] || [[ "x${queue}" == "x" ]] || [[ "x${rscript}" == "x" ]]
then
   echo "---------------------------------------------------------------------------------"
   echo "    The following variables must be set.  In case any of them are empty, check"
   echo " your script sim_sitter.sh."
   echo " "
   echo " email   = ${email}"
   echo " queue   = ${queue}"
   echo " rscript = ${rscript}"
   echo " "
   echo "---------------------------------------------------------------------------------"
   exit 99
fi
#------------------------------------------------------------------------------------------#



#----- Make substitutions. ----------------------------------------------------------------#
sed -i~ s@"here=\"\""@"here=\"${here}\""@g                       ${here}/run_sitter.sh
sed -i~ s@"recipient=\"\""@"recipient=\"${email}\""@g            ${here}/run_sitter.sh
sed -i~ s@"frqemail=\"\""@"frqemail=${frqemail}"@g               ${here}/run_sitter.sh
sed -i~ s@"delay1st_min=\"\""@"delay1st_min=${delay1st_min}"@g   ${here}/run_sitter.sh
sed -i~ s@"wait_minutes=\"\""@"wait_minutes=${wait_minutes}"@g   ${here}/run_sitter.sh
sed -i~ s@"frqpost=\"\""@"frqpost=${frqpost}"@g                  ${here}/run_sitter.sh
sed -i~ s@"frqtouch=\"\""@"frqtouch=${frqtouch}"@g               ${here}/run_sitter.sh
sed -i~ s@"here=\"\""@"here=\"${here}\""@g                       ${here}/epost.sh
sed -i~ s@"global_queue=\"\""@"global_queue=\"${queue}\""@g      ${here}/epost.sh
sed -i~ s@"rscript=\"\""@"rscript=\"${rscript}\""@g              ${here}/epost.sh
sed -i~ s@"submit=false"@"submit=true"@g                         ${here}/epost.sh
sed -i~ s@"here=\"\""@"here=\"${here}\""@g                       ${here}/transfer.sh
sed -i~ s@"there=\"\""@"there=\"${there}\""@g                    ${here}/transfer.sh
sed -i~ s@"here=\"\""@"here=\"${here}\""@g                       ${here}/last_histo.sh
sed -i~ s@"there=\"\""@"there=\"${there}\""@g                    ${here}/last_histo.sh
sed -i~ s@"checkhourly=\"\""@"checkhourly=\"${checkhourly}\""@g  ${here}/last_histo.sh
sed -i~ s@"checkstatus=\"\""@"checkstatus=\"${checkstatus}\""@g  ${here}/last_histo.sh
#------------------------------------------------------------------------------------------#


#----- Job preferences. -------------------------------------------------------------------#
joblog="${here}/out_sitter.out"
jobpref=$(basename ${here})
jobname="${jobpref}-control"
jobmain="${jobpref}-sims"
#------------------------------------------------------------------------------------------#


#------ Check ID of main job name, so it includes a dependency. ---------------------------#
mainid=$(sjobs | grep -i ${jobmain} | grep -v CANCELLED | grep -v COMPLETED | grep -v COMPLETING | grep -v FAILED | awk '{print $1}')
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#    Check whether to submit the job.                                                      #
#------------------------------------------------------------------------------------------#
if [[ "x${mainid}" == "x" ]]
then
   echo " Job name ${jobmain} was not found, so the run_sitter.sh job will not be queued."
   goahead="n"
elif [[ -s "${here}/run_sitter.lock" ]] || [[ -s "${here}/transfer.lock" ]]
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
comm="${sbatch} -p ${queue} --mem-per-cpu=${memory} -t ${runtime} -o ${joblog}"
comm="${comm} -J ${jobname} -n 1 --dependency=after:${mainid}"
comm="${comm} --wrap=\"${here}/run_sitter.sh\""
case ${goahead} in
y|yes)
   /bin/rm -f ${here}/run_sitter.lock
   /bin/rm -f ${here}/transfer.lock
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
