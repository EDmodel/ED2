#!/bin/bash
. ~/.bashrc

#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#    MAIN SETTINGS:                                                                        #
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#


#----- Main path, you must actually type the path in case you add to chron. ---------------#
here=""
#------------------------------------------------------------------------------------------#


#----- Path with some utilities for run_sitter.sh (this script). --------------------------#
situtils="${here}/sit_utils"
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#   desc     -- Unique job description for these simulations (it really must be unique).   #
#               Normally it is the path name.                                              #
#   runtitle -- Full name of this simulation, this is used only in the e-mail subject.     #
#------------------------------------------------------------------------------------------#
desc=$(basename ${here})
jobname="${desc}-sims*"
runtitle="Simulation group: ${desc}"
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#     Useful files to check the simulation progress:                                       #
# joborder  -- File with the job instructions and model settings                           #
# lastcheck -- File with the report from last time                                         #
# outcheck  -- File with the current report                                                #
# situation -- Situation of the run                                                        #
#------------------------------------------------------------------------------------------#
joborder="${here}/joborder.txt"
lastcheck="${situtils}/lastcheck.txt"
outcheck="${situtils}/mycheck.txt"
check_out="/tmp/check_run_${$}.out"
#------------------------------------------------------------------------------------------#


#------ User name. ------------------------------------------------------------------------#
moi=$(whoami)  # User name
#------------------------------------------------------------------------------------------#


#------ Calculator. -----------------------------------------------------------------------#
ccc="${HOME}/util/calc.sh"  # Calculator
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#      Variables that tell whether to check for steady state and deserts.                  #
#------------------------------------------------------------------------------------------#
checksteady="FALSE"   #   Check whether to stop simulations? (TRUE of FALSE, R style)
nyearmin=160          #   Minimum number of years before we check
ststcrit=0.01         #   Maximum change allowed between two cycles
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#    E-mail options (normally only the first three options may require changes.            #
#    emailevery   -- Should I e-mail everytime (true) or every frqemail (false)?           #
#    recipient    -- To which e-mail should I send the update?                             #
#    mailprog     -- Which e-mail program to use? mutt works best, I installed in my util  #
#                    directory.  Give a try, if it doesn't work you may need to install    #
#                    locally on your directory.                                            #
#    frqemail     -- minimum time interval between emails (in seconds).                    #
#    plotstatus   -- Plot the current status of the run (you will need to create some R    #
#                    script for each simulation). 0: no; 1: yes                            #
#    Rscript_plot -- Script that you want to run to generate some plots.                   #
#    R_figlist    -- List with figure names (you must list all files you want to append)   #
#    emailbody    -- File that will contain the e-mail.                                    #
#    headfile     -- File with header                                                      #
#    tailfile     -- File with "bye"                                                       #
#    recefile     -- File with the recent activity                                         #
#    statfile     -- File with the current status                                          #
#    queuefile    -- File with queue status                                                #
#    pendfile     -- File with pending status                                              #
#    email1day    -- Reminder so the script knows whether an e-mail has been sent or not.  #
#------------------------------------------------------------------------------------------#
emailevery=false
recipient=""
mailprog=$(which mutt)
frqemail=""
plotstatus=true
Rscript_plot="${situtils}/plot.region.r"
R_figlist="${situtils}/status_region.png
           ${situtils}/agb_region.png
           ${situtils}/bsa_region.png
           ${situtils}/lai_region.png
           ${situtils}/scb_region.png
           ${situtils}/npa_region.png"
emailbody="${situtils}/email.txt"
headfile="${situtils}/head.txt"
tailfile="${situtils}/tail.txt"
recefile="${situtils}/rece.txt"
statfile="${situtils}/stat.txt"
queuefile="${situtils}/queue.txt"
tablefile="${situtils}/table_queue.txt"
pendfile="${here}/pending.txt"
emailday="${here}/emailday.txt"
#------------------------------------------------------------------------------------------#


#------ Simulation settings. --------------------------------------------------------------#
delay1st_min=""
wait_minutes=""
frqpost=""
frqtouch=""
#------------------------------------------------------------------------------------------#


#----- Format for checking simulations. ---------------------------------------------------#
outform="JobId%20,JobName%200,State%12"
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#     Optional feature to reset stalled polygons.  This feature is particularly useful     #
# when using serial_requeue, as some jobs land in faulty nodes.  This feature is currently #
# available for CANNON only.                                                               #
#------------------------------------------------------------------------------------------#
stall_reset=true
nstall=2
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


#------------------------------------------------------------------------------------------#
#       First check that the main path and e-mail have been set.  If not, don't run.       #
#------------------------------------------------------------------------------------------#
if [[ "x${here}"         == "x" ]] || [[ "x${recipient}"    == "x" ]] ||
   [[ "x${frqemail}"     == "x" ]] || [[ "x${delay1st_min}" == "x" ]] ||
   [[ "x${wait_minutes}" == "x" ]] || [[ "x${frqpost}"      == "x" ]] ||
   [[ "x${frqtouch}"     == "x" ]]
then
   echo "-------------------=========---------------------------------------------------"
   echo "    The following variables must be set.  In case any of them are empty, check "
   echo " your script run_sitter.sh."
   echo " "
   echo " here         = ${here}"
   echo " recipient    = ${recipient}"
   echo " frqemail     = ${frqemail}"
   echo " delay1st_min = ${delay1st_min}"
   echo " wait_minutes = ${wait_minutes}"
   echo " frqpost      = ${frqpost}"
   echo " frqtouch     = ${frqtouch}"
   echo "-------------------------------------------------------------------------------"
   exit 99
fi
#------------------------------------------------------------------------------------------#




#----- Check whether run_sitter.sh is still running or not.  If it is, exit. --------------#
if [[ -s ${here}/run_sitter.lock ]]
then
   exit
else
   echo "I am going to check your run." > ${here}/run_sitter.lock
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
   stall_reset=false
   ;;
*)
   echo -n "Failed guessing cluster from node name.  Please type the name:   "
   read cluster
   ;;
esac
#------------------------------------------------------------------------------------------#



#------ Move to the current directory. ----------------------------------------------------#
cd ${here}
#------------------------------------------------------------------------------------------#


#------ Make sure transfer is set up at the right path. -----------------------------------#
transfer="${here}/transfer.sh"
translock="${here}/transfer.lock"
hereline=$(cat ${here}/transfer.sh | grep "^here=")
sed -i~ s@"${hereline}"@"here=\"${here}\""@g ${transfer}
#------------------------------------------------------------------------------------------#


#----- Determine the number of polygons to run. -------------------------------------------#
let n_polygon=$(wc -l ${joborder} | awk '{print $1 }')-3
if [[ ${n_polygon} -lt 100 ]]
then
   ndig=2
elif [[ ${n_polygon} -lt 1000 ]]
then
   ndig=3
elif [[ ${n_polygon} -lt 10000 ]]
then
   ndig=4
else
   ndig=5
fi
pfmt="%${ndig}.${ndig}i"
echo "Number of polygons: ${n_polygon}..."
#------------------------------------------------------------------------------------------#



#----- Set job prefix. --------------------------------------------------------------------#
desc=$(basename ${here})
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#      Check whether to wait a bit before the first check.                                 #
#------------------------------------------------------------------------------------------#
if [[ ${delay1st_min} -gt 0 ]]
then
   let nsl=${delay1st_min}+1
   while [[ ${nsl} -gt 1 ]]
   do
      let nsl=${nsl}-1
      case ${nsl} in
      1)
         echo "       - ${nsl} minute before starting."
         ;;
      *)
         echo "       - ${nsl} minutes before starting."
         ;;
      esac
      sleep 1m
   done
fi
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#       Main loop: the script stays inside until all simulations have finished.            #
#------------------------------------------------------------------------------------------#
n_ongoing=9999
iter=0
while [[ ${n_ongoing} -gt 0 ]]
do
   #------ Update iteration counter. ------------------------------------------------------#
   let iter=${iter}+1
   #---------------------------------------------------------------------------------------#

   #---- Reset count of polygons. ---------------------------------------------------------#
   n_running=0
   n_stalled=0
   n_suspend=0
   n_the_end=0
   n_unknown=0
   n_sigsegv=0
   n_crashed=0
   n_hydfail=0
   n_metmiss=0
   n_stopped=0
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Check whether it's time to update simulation setting files.                       #
   #---------------------------------------------------------------------------------------#
   if [[ ${frqtouch} -gt 0 ]]
   then
      let xtouch=${iter}%${frqtouch}
      if [[ ${xtouch} -eq 0 ]]
      then
         echo " Touch main files."
         touch ${here}/*.sh
         touch ${here}/*.r
         touch ${here}/*.txt
         touch ${here}/Template/*
         touch ${here}/executable/*
         touch ${here}/sit_utils/*
      fi
   else
      xtouch=999
   fi
   #---------------------------------------------------------------------------------------#


   #----- Transfer files. -----------------------------------------------------------------#
   if [[ ! -s ${translock} ]]
   then
      echo " Backup files."
      nice nohup ${transfer} 1> ${here}/out_transfer.out 2>&1 &
   fi
   #---------------------------------------------------------------------------------------#


   #----- Move current check to the last check. -------------------------------------------#
   if [[ -s "${outcheck}" ]]
   then
      nlnout=$(cat ${outcheck} | wc -l)
      if [[ ${nlnout} -eq ${n_polygon} ]]
      then
         rm -f ${lastcheck}
         mv ${outcheck} ${lastcheck}
         touch ${outcheck}
      else
         rm -f ${outcheck}
         touch ${outcheck}
      fi
   fi
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Go through all polygons.                                                          #
   #---------------------------------------------------------------------------------------#
   ff=0
   while [[ ${ff} -lt ${n_polygon} ]]
   do

      #------------------------------------------------------------------------------------#
      #    Format count.                                                                   #
      #------------------------------------------------------------------------------------#
      let ff=${ff}+1
      let line=${ff}+3
      ffout=$(printf "${pfmt}" ${ff})
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #      Read the ffth line of the polygon list.  There must be smarter ways of doing  #
      # this, but this works.  Here we obtain the polygon name, and its longitude and      #
      # latitude.                                                                          #
      #------------------------------------------------------------------------------------#
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
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #    Task name.  This depends on the type of submission (just to avoid clutter).     #
      #------------------------------------------------------------------------------------#
      case "${cluster}" in
      SDUMONT) 
         #----- Task name is bound to the global jobname, no need to add prefix. ----------#
         taskname="${polyname}"
         #---------------------------------------------------------------------------------#
         ;;
      CANNON )
         #----- Add prefix to the task name, which will name this job. --------------------#
         taskname="${desc}-${polyname}"
         #---------------------------------------------------------------------------------#
         ;;
      esac
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Determine the scenario paths.                                                  #
      #------------------------------------------------------------------------------------#
      case ${iscenario} in
      default)
         case ${metdriver} in
         ERA5_CHIRPS)
            #----- ERA5 (CHIRPS precipitation, Brazilian Amazon only). --------------------#
            scentype="ERA5"
            iscenario="ERA5_SOUTHAM_CHIRPS"
            ;;
         ERAINT_CHIRPS)
            #----- ERA-Interim (CHIRPS precipitation). ------------------------------------#
            scentype="ERA_Interim"
            iscenario="ERAINT_SOUTHAM_CHIRPS"
            ;;
         ERAINT_MSWEP2)
            #----- ERA-Interim (MSWEP2 precipitation). ------------------------------------#
            scentype="ERA_Interim"
            iscenario="ERAINT_SOUTHAM_MSWEP2"
            ;;
         ERAINT_NATIVE)
            #----- ERA-Interim (native precipitation). ------------------------------------#
            scentype="ERA_Interim"
            iscenario="ERAINT_SOUTHAM_NATIVE"
            ;;
         MERRA2_CHIRPS)
            #----- MERRA2 (CHIRPS precipitation). -----------------------------------------#
            scentype="MERRA2"
            iscenario="MERRA2_SOUTHAM_CHIRPS"
            ;;
         MERRA2_MSWEP2)
            #----- MERRA2 (MSWEP2 precipitation). -----------------------------------------#
            scentype="MERRA2"
            iscenario="MERRA2_SOUTHAM_MSWEP2"
            ;;
         MERRA2_NATIVE)
            #----- MERRA-2 (native precipitation). ----------------------------------------#
            scentype="MERRA2"
            iscenario="MERRA2_SOUTHAM_NATIVE"
            ;;
         PGMF3_CHIRPS)
            #----- PGMF-3 (CHIRPS precipitation). -----------------------------------------#
            scentype="PGMF3"
            iscenario="PGMF3_SOUTHAM_CHIRPS"
            ;;
         PGMF3_MSWEP2)
            #----- PGMF-3 (CHIRPS precipitation). -----------------------------------------#
            scentype="PGMF3"
            iscenario="PGMF3_SOUTHAM_MSWEP2"
            ;;
         PGMF3_NATIVE)
            #----- PGMF-3 (native precipitation). -----------------------------------------#
            scentype="PGMF3"
            iscenario="PGMF3_SOUTHAM_NATIVE"
            ;;
         Sheffield)
            #----- Sheffield. -------------------------------------------------------------#
            scentype="sheffield"
            iscenario="sheffield"
            ;;
         WFDE5_CHIRPS)
            #----- WFDEI (CHIRPS Precipitation). ------------------------------------------#
            scentype="WFDE5"
            iscenario="WFDE5_SOUTHAM_CHIRPS"
            ;;
         WFDEI_CHIRPS)
            #----- WFDEI (CHIRPS Precipitation). ------------------------------------------#
            scentype="WFDEI"
            iscenario="WFDEI_SOUTHAM_CHIRPS"
            ;;
         WFDEI_CRUP)
            #----- WFDEI (CRU Precipitation). ---------------------------------------------#
            scentype="WFDEI"
            iscenario="WFDEI_SOUTHAM_CRUP"
            ;;
         WFDEI_GPCC)
            #----- WFDEI (GPCC Precipitation). --------------------------------------------#
            scentype="WFDEI"
            iscenario="WFDEI_SOUTHAM_GPCC"
            ;;
         WFDEI_MSWEP2)
            #----- WFDEI (MSWEP2 Precipitation). ------------------------------------------#
            scentype="WFDEI"
            iscenario="WFDEI_SOUTHAM_MSWEP2"
            ;;
         *)
            #----- Tower data. ------------------------------------------------------------#
            scentype="wmo+eft"
            iscenario="eft"
            ;;
         esac
         ;;
      eft|wmo|shr)
         #----- Tower data, keep scenario as is. ------------------------------------------#
         scentype="wmo+eft"
         ;;
      *)
         #----- Rainfall scenario, keep scenario as is. -----------------------------------#
         scentype="realisation_scen_driver"
         ;;
      esac
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #      Choose the scenario to use.  Iscenario follows the following convention:      #
      # "default"    -- No scenario.  Use the tower/Sheffield data.                        #
      # "wmo"        -- No scenario.  Use the WMO-based data.                              #
      # "rRRR_tTTT   -- Use scenarios, with rRRRR controlling the rainfall, and tTTTT      #
      #                 controlling temperature.                                           #
      #                                                                                    #
      # rRRR         -- Rainfall scenarios, where rRRRR means:                             #
      #                 r+000: INMET-based time series, no resampling.                     #
      #                 r+010: Re-sampling with substitution, but equal chances for all    #
      #                        years                                                       #
      #                 r-XXX: Shift the location (similar to mean) of the distribution by #
      #                        -X.XX units of scale (similar to standard deviation), so    #
      #                        the time series becomes drier.                              #
      #                 r+XXX: Similar to above, but make the time series wetter.          #
      #                                                                                    #
      # tTTT         -- Temperature scenarios, where tTTTT means:                          #
      #                 t+000: No change in temperature                                    #
      #                 t-YYY: Change temperature by -Y.YY Kelvin.  Keep relative humidity #
      #                        the same and correct specific humidity.                     #
      #                 t+YYY: Change temperature by +Y.YY Kelvin.  Keep relative humidity #
      #                        the same and correct by +Y.YY Kelvin.                       #
      #                 r+XXX: Similar to above, but make the time series wetter.          #
      #------------------------------------------------------------------------------------#
      #----- Find out which scenario to use. ----------------------------------------------#
      fullscen="${metmain}/met_driver/${scentype}/${iscenario}"
      #------------------------------------------------------------------------------------#
      #     Determine which meteorological data set to use.  Default is the Sheffield/NCEP #
      # dataset, otherwise the site-level tower data is used.                              #
      #------------------------------------------------------------------------------------#
      case ${metdriver} in
      Bananal)
         metdriverdb="${fullscen}/Bananal/Bananal_HEADER"
         metcyc1=2004
         metcycf=2006
         imetavg=1
         ;;
      Brasilia)
         metdriverdb="${fullscen}/Brasilia/Brasilia_HEADER"
         metcyc1=2006
         metcycf=2012
         imetavg=1
         ;;
      Caxiuana)
         metdriverdb="${fullscen}/Caxiuana/Caxiuana_HEADER"
         metcyc1=1999
         metcycf=2003
         imetavg=1
         ;;
      ERA5_CHIRPS)
         metdriverdb="${fullscen}/${iscenario}_HEADER"
         metcyc1=1981
         metcycf=2019
         imetavg=1
         ;;
      ERAINT_CHIRPS)
         metdriverdb="${fullscen}/${iscenario}_HEADER"
         metcyc1=1981
         metcycf=2017
         imetavg=2
         ;;
      ERAINT_NATIVE)
         metdriverdb="${fullscen}/${iscenario}_HEADER"
         metcyc1=1979
         metcycf=2017
         imetavg=2
         ;;
      ERAINT_MSWEP2)
         metdriverdb="${fullscen}/${iscenario}_HEADER"
         metcyc1=1979
         metcycf=2016
         imetavg=2
         ;;
      Fazenda_Nossa_Senhora)
         metdriverdb="${fullscen}/Fazenda_Nossa_Senhora/Fazenda_Nossa_Senhora_HEADER"
         metcyc1=1999
         metcycf=2002
         imetavg=1
         ;;
      Harvard)
         metdriverdb="${fullscen}/Harvard/Harvard_HEADER"
         metcyc1=1992
         metcycf=2003
         imetavg=3
         ;;
      Laegern)
         metdriverdb="${fullscen}/Laegern/Laegern_HEADER"
         metcyc1=2006
         metcycf=2016
         imetavg=1
         ;;
      Manaus_Km34)
         metdriverdb="${fullscen}/Manaus_Km34/Manaus_Km34_HEADER"
         metcyc1=1999
         metcycf=2006
         imetavg=1
         ;;
      MERRA2_CHIRPS)
         metdriverdb="${fullscen}/${iscenario}_HEADER"
         metcyc1=1981
         metcycf=2017
         imetavg=3
         ;;
      MERRA2_MSWEP2)
         metdriverdb="${fullscen}/${iscenario}_HEADER"
         metcyc1=1980
         metcycf=2016
         imetavg=3
         ;;
      MERRA2_NATIVE)
         metdriverdb="${fullscen}/${iscenario}_HEADER"
         metcyc1=1980
         metcycf=2017
         imetavg=3
         ;;
      Natal)
         metdriverdb="${fullscen}/Natal/Natal_HEADER"
         metcyc1=2009
         metcycf=2012
         imetavg=1
         ;;
      Paracou)
         metdriverdb="${fullscen}/Paracou/Paracou_HEADER"
         metcyc1=2004
         metcycf=2014
         imetavg=1
         ;;
      Pe-de-Gigante)
         metdriverdb="${fullscen}/Pe-de-Gigante/Pe-de-Gigante_HEADER"
         metcyc1=2001
         metcycf=2003
         imetavg=1
         ;;
      Petrolina)
         metdriverdb="${fullscen}/Petrolina/Petrolina_HEADER"
         metcyc1=2004
         metcycf=2012
         imetavg=1
         ;;
      PGMF3_CHIRPS)
         metdriverdb="${fullscen}/${iscenario}_HEADER"
         metcyc1=1981
         metcycf=2016
         imetavg=3
         ;;
      PGMF3_MSWEP2)
         metdriverdb="${fullscen}/${iscenario}_HEADER"
         metcyc1=1979
         metcycf=2016
         imetavg=3
         ;;
      PGMF3_NATIVE)
         metdriverdb="${fullscen}/${iscenario}_HEADER"
         metcyc1=1979
         metcycf=2016
         imetavg=3
         ;;
      Rebio_Jaru)
         metdriverdb="${fullscen}/Rebio_Jaru/Rebio_Jaru_HEADER"
         metcyc1=1999
         metcycf=2002
         imetavg=1
         ;;
      Santarem_Km67)
         metdriverdb="${fullscen}/Santarem_Km67/Santarem_Km67_HEADER"
         metcyc1=2001
         metcycf=2011
         imetavg=1
         ;;
      Santarem_Km77)
         metdriverdb="${fullscen}/Santarem_Km77/Santarem_Km77_HEADER"
         metcyc1=2001
         metcycf=2005
         imetavg=1
         ;;
      Santarem_Km83)
         metdriverdb="${fullscen}/Santarem_Km83/Santarem_Km83_HEADER"
         metcyc1=2000
         metcycf=2003
         imetavg=1
         ;;
      Sheffield)
         metdriverdb=${fullscen}/${shefhead}
         metcyc1=1969
         metcycf=2008
         imetavg=2
         ;;
      Tanguro_Burn)
         metdriverdb="${fullscen}/Tanguro_Burn/Tanguro_Burn_HEADER"
         metcyc1=2008
         metcycf=2018
         imetavg=1
         ;;
      Tanguro_Ctrl)
         metdriverdb="${fullscen}/Tanguro_Ctrl/Tanguro_Ctrl_HEADER"
         metcyc1=2008
         metcycf=2018
         imetavg=1
         ;;
      Tonzi)
         metdriverdb="${fullscen}/Tonzi/Tonzi_HEADER"
         metcyc1=2000
         metcycf=2010
         imetavg=1
         ;;
      WFDE5_CHIRPS)
         metdriverdb="${fullscen}/${iscenario}_HEADER"
         metcyc1=1981
         metcycf=2018
         imetavg=2
         ;;
      WFDEI_CHIRPS)
         metdriverdb="${fullscen}/${iscenario}_HEADER"
         metcyc1=1981
         metcycf=2016
         imetavg=1
         ;;
      WFDEI_CRUP)
         metdriverdb="${fullscen}/${iscenario}_HEADER"
         metcyc1=1979
         metcycf=2016
         imetavg=1
         ;;
      WFDEI_GPCC)
         metdriverdb="${fullscen}/${iscenario}_HEADER"
         metcyc1=1979
         metcycf=2016
         imetavg=1
         ;;
      WFDEI_MSWEP2)
         metdriverdb="${fullscen}/${iscenario}_HEADER"
         metcyc1=1979
         metcycf=2016
         imetavg=1
         ;;
      *)
         echo "Met driver: ${metdriver}"
         echo "Sorry, this met driver is not valid for regular runs"
         exit 85
         ;;
      esac
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Correct years so it is not tower-based or Sheffield.                           #
      #------------------------------------------------------------------------------------#
      case ${iscenario} in
      default|eft|shr|ERA5*|ERAINT*|MERRA2*|PGMF3*|WFDE5*|WFDEI*|Sheffield)
         echo "Nothing" > /dev/null
         ;;
      *)
         metcyc1=1972
         metcycf=2012
         imetavg=1
         ;;
      esac
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Set some variables to check whether the simulation is running.                 #
      #------------------------------------------------------------------------------------#
      stdout="${here}/${polyname}/serial_out.out"
      stderr="${here}/${polyname}/serial_out.err"
      lsfout="${here}/${polyname}/serial_lsf.out"
      #------------------------------------------------------------------------------------#



      #----- Update the R script that controls the run. -----------------------------------#
      /bin/rm -f ${here}/${polyname}/whichrun.r
      /bin/cp -f ${here}/Template/whichrun.r    ${here}/${polyname}/whichrun.r
      sed -i s@thispoly@${polyname}@g           ${here}/${polyname}/whichrun.r
      sed -i s@thisqueue@${queue}@g             ${here}/${polyname}/whichrun.r
      sed -i s@pathhere@${here}@g               ${here}/${polyname}/whichrun.r
      sed -i s@paththere@${here}@g              ${here}/${polyname}/whichrun.r
      sed -i s@thisyeara@${yeara}@g             ${here}/${polyname}/whichrun.r
      sed -i s@thismontha@${montha}@g           ${here}/${polyname}/whichrun.r
      sed -i s@thisdatea@${datea}@g             ${here}/${polyname}/whichrun.r
      sed -i s@thistimea@${timea}@g             ${here}/${polyname}/whichrun.r
      sed -i s@thischecksteady@${checksteady}@g ${here}/${polyname}/whichrun.r
      sed -i s@thismetcyc1@${metcyc1}@g         ${here}/${polyname}/whichrun.r
      sed -i s@thismetcycf@${metcycf}@g         ${here}/${polyname}/whichrun.r
      sed -i s@thisnyearmin@${nyearmin}@g       ${here}/${polyname}/whichrun.r
      sed -i s@thisststcrit@${ststcrit}@g       ${here}/${polyname}/whichrun.r
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #    Retrieve the information on statusrun.txt because it may tell whether the run   #
      # was in steady state or all PFTs went extinct, and we don't need to check them      #
      # again.                                                                             #
      #------------------------------------------------------------------------------------#
      statrun=${here}/${polyname}/statusrun.txt
      if [[ -s ${statrun} ]]
      then
         #----- Obtain previous status. ---------------------------------------------------#
         yearh_old=$(cat ${statrun}  | awk '{print  $2}')
         monthh_old=$(cat ${statrun} | awk '{print  $3}')
         dateh_old=$(cat ${statrun}  | awk '{print  $4}')
         timeh_old=$(cat ${statrun}  | awk '{print  $5}')
         runt_old=$(cat ${statrun}   | awk '{print  $6}')
         agb_old=$(cat ${statrun}    | awk '{print  $7}')
         bsa_old=$(cat ${statrun}    | awk '{print  $8}')
         lai_old=$(cat ${statrun}    | awk '{print  $9}')
         scb_old=$(cat ${statrun}    | awk '{print $10}')
         npa_old=$(cat ${statrun}    | awk '{print $11}')
         #---------------------------------------------------------------------------------#



         #----- In case the simulation never started, force it to be INITIAL. -------------#
         if [[ ${yearh_old} -eq ${yeara} ]] && [[ ${monthh_old} -eq ${montha} ]]
         then
            yearh_old=${yeara}
            monthh_old=${montha}
            dateh_old=${datea}
            timeh_old=${timea}
            runt_old="INITIAL"
            agb_old="NA"
            bsa_old="NA"
            lai_old="NA"
            scb_old="NA"
            npa_old="NA"
         fi
         #---------------------------------------------------------------------------------#
      else
         yearh_old=${yeara}
         monthh_old=${montha}
         dateh_old=${datea}
         timeh_old=${timea}
         runt_old="INITIAL"
         agb_old="NA"
         bsa_old="NA"
         lai_old="NA"
         scb_old="NA"
         npa_old="NA"
      fi
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Initially assume stall=0.  This will be checked and updated only if the        #
      # simulation is running.  Also retrieve the previous stall count if there is one,    #
      # in case we must update the stall check.                                            #
      #------------------------------------------------------------------------------------#
      stall=0
      #----- Load previous stall. ---------------------------------------------------------#
      if [[ -s ${lastcheck} ]]
      then
         stall_old=$(cat ${lastcheck} | grep ${polyname} | awk '{print $8}')
      else
         stall_old=0
      fi
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #    We only check those polygons that may be running, so the check doesn't take too #
      # long.  Once a polygon has reached a steady state / gone extinct / finished, then   #
      # its  #
      # status should remain the same.                                                        #
      #------------------------------------------------------------------------------------#
      #----- Call R to check status. ------------------------------------------------------#
      whichrun=${here}/${polyname}/whichrun.r
      outwhich=${here}/${polyname}/outwhich.txt
      R CMD BATCH --no-save --no-restore ${whichrun} ${outwhich}
      if [[ -s ${statrun} ]]
      then
         yearh=$(cat ${statrun}  | awk '{print  $2}')
         monthh=$(cat ${statrun} | awk '{print  $3}')
         dateh=$(cat ${statrun}  | awk '{print  $4}')
         timeh=$(cat ${statrun}  | awk '{print  $5}')
         runt=$(cat ${statrun}   | awk '{print  $6}')
         agb=$(cat ${statrun}    | awk '{print  $7}')
         bsa=$(cat ${statrun}    | awk '{print  $8}')
         lai=$(cat ${statrun}    | awk '{print  $9}')
         scb=$(cat ${statrun}    | awk '{print $10}')
         npa=$(cat ${statrun}    | awk '{print $11}')
      else
         yearh=${yeara}
         monthh=${montha}
         dateh=${datea}
         timeh=${timea}
         runt="INITIAL"
         agb="NA"
         bsa="NA"
         lai="NA"
         scb="NA"
         npa="NA"
      fi
      #------------------------------------------------------------------------------------#




      #----- Update files that may be used to check status or to re-submit. ---------------#
      if [[ ${xtouch} -eq 0 ]]
      then
         touch ${here}/${polyname}/*.r
         touch ${here}/${polyname}/ED2IN
         touch ${here}/${polyname}/*.sh
         touch ${here}/${polyname}/SHEF_NCEP_DRIVER_DS314
         touch ${here}/${polyname}/statusrun.txt
      fi
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Check whether the simulation is still running, and if not, why it isn't.       #
      #------------------------------------------------------------------------------------#
      if [[ -s ${stdout} ]]
      then
         #----- Set task inquire command. -------------------------------------------------#
         case "${cluster}" in
         SDUMONT)
            #----- Multi-task job, search for specific task. ------------------------------#
            stask="stask --noheader -u ${moi} -t ${taskname} -j ${jobname}"
            #------------------------------------------------------------------------------#
            ;;
         CANNON)
            #----- Single-task jobs, search the task using the job name. ------------------#
            stask="stask --noheader -u ${moi} -j ${taskname}"
            #------------------------------------------------------------------------------#
            ;;
         esac
         #---------------------------------------------------------------------------------#


         #----- Check whether the simulation is running, and when in model time it is. ----#
         running=$(${stask}   -o "${outform}" | grep "RUNNING"   | wc -l)
         pending=$(${stask}   -o "${outform}" | grep "PENDING"   | wc -l)
         suspended=$(${stask} -o "${outform}" | grep "SUSPENDED" | wc -l)
         simline=$(grep "Simulating: "   ${stdout} | tail -1)
         runtime=$(echo ${simline} | awk '{print $3}')
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Check whether the simulations are stalled.                                 #
         #---------------------------------------------------------------------------------#
         if [[ ${running} -gt 0 ]]
         then
            #----- Check whether the runs are stalled. ------------------------------------#
            if [[ ${yearh} == ${yearh_old} ]] && [[ ${monthh} == ${monthh_old} ]] &&
               [[ ${dateh} == ${dateh_old} ]] && [[ ${timeh}  == ${timeh_old}  ]]
            then
               let stall=${stall_old}+1
            else
               stall=0
            fi
            #------------------------------------------------------------------------------#




            #------------------------------------------------------------------------------#
            #     Check whether to kill and resubmit the job.                              #
            #------------------------------------------------------------------------------#
            if ${stall_reset} && [[ ${stall} -ge ${nstall} ]]
            then
               #------- Set stalled flag to one so it is clear the run is not running. ----#
               stalled=true
               #---------------------------------------------------------------------------#
            else
               #------- Set stalled flag to zero. It is either running or something else. -#
               stalled=false
               #---------------------------------------------------------------------------#
            fi
            #------------------------------------------------------------------------------#



         else
            #----- Not in history, simulations are unlikely to be stalled. ----------------#
            stall=0
            stalled=false
            #------------------------------------------------------------------------------#
         fi
         #---------------------------------------------------------------------------------#




         #----- Check for segmentation violations. ----------------------------------------#
         if [[ -s ${stderr} ]]
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



         #----- Check for other possible outcomes. ----------------------------------------#
         stopped=$(grep "FATAL ERROR"                       ${stdout} | wc -l)
         crashed=$(grep "IFLAG1 problem."                   ${stdout} | wc -l)
         hydfail=$(grep "Plant Hydrodynamics is off-track." ${stdout} | wc -l)
         the_end=$(grep "ED-2.2 execution ends"             ${stdout} | wc -l)
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #     Plot a message so the user knows what is going on.                          #
         #---------------------------------------------------------------------------------#
         if [[ ${pending} -gt 0 ]]
         then
            let n_pending=${n_pending}+1
            echo -e "${ffout}: ${polyname} is pending..."
         elif [[ ${suspended} -gt 0 ]]
         then
            let n_suspend=${n_suspend}+1
            echo -e "${ffout}: ${polyname} is SUSPENDED."
         elif ${stalled}
         then
            let n_stalled=${n_stalled}+1
            echo -e "${ffout}: ${polyname} is STALLED (${stall} - ${runtime})."

            #------------------------------------------------------------------------------#
            #      Kill jobs, it must be stalled.                                          #
            #------------------------------------------------------------------------------#
            jobid="enter_loop"
            kcnt=0
            while [[ "${jobid}" != "" ]] && [[ ${kcnt} -lt 10 ]]
            do
               #----- Ensure to eliminate potential duplicates. ---------------------------#
               let kcnt=${kcnt}+1
               jobln=$(${stask} -o "${outform}" | grep -v "CANCEL" | grep -v "COMPLET"     \
                      | head -1)
               jobid=$(echo ${jobln} | awk '{print $1}')
               if [[ "${jobid}" != "" ]]
               then
                  echo "   - Cancel Job ${jobid} (${taskname})."
                  scancel ${jobid}
               fi
               #---------------------------------------------------------------------------#
            done
            if [[ ${kcnt} -ge 10 ]]
            then
               #------ The script is bogus. -----------------------------------------------#
               echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
               echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
               echo "  JOB NUMBER IS LIKELY INCORRECT IN STASK. STOPPING THE SCRIPT.    "
               echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
               echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
               echo "   Output of stask command:"
               ${stask} -o "${outform}"
               echo "   Resulting JobId:"
               jobln=$(${stask} -o "${outform}" | grep -v "CANCEL" | grep -v "COMPLET"     \
                      | head -1)
               jobid=$(echo ${jobln} | awk '{print $1}')
               echo " Job Line=${jobln}"
               echo " JobId=${jobid}"
               echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
               exit 99
               #---------------------------------------------------------------------------#
            fi
            #------------------------------------------------------------------------------#


            #----- Re-submit job. ---------------------------------------------------------#
            echo "   - Submit New job (${taskname})."
            callserial="${here}/${polyname}/callserial.sh"
            sbatch  ${callserial}
            stall=0
            #------------------------------------------------------------------------------#

         elif [[ ${running} -gt 0 ]] && [[ ${sigsegv} -eq 0 ]]
         then
            let n_running=${n_running}+1
            echo -e "${ffout}: ${polyname} is running (${runtime})."
         elif [[ ${sigsegv} -gt 0 ]]
         then
            let n_sigsegv=${n_sigsegv}+1
            echo -e "${ffout}: ${polyname} HAD SEGMENTATION VIOLATION."
         elif [[ ${crashed} -gt 0 ]]
         then 
            let n_crashed=${n_crashed}+1
            echo -e "${ffout}: ${polyname} HAS CRASHED (RK4 PROBLEM)."
         elif [[ ${hydfail} -gt 0 ]]
         then 
            let n_hydfail=${n_hydfail}+1
            echo -e "${ffout}: ${polyname} HAS CRASHED (HYDRODYNAMICS)."
         elif [[ ${metmiss} -gt 0 ]]
         then 
            let n_metmiss=${n_metmiss}+1
            echo -e "${ffout}: ${polyname} DID NOT FIND MET DRIVERS."
         elif [[ ${stopped} -gt 0 ]]
         then
            let n_stopped=${n_stopped}+1
            echo -e "${ffout}: ${polyname} STOPPED (UNKNOWN REASON)."
         elif [[ ${the_end} -gt 0 ]]
         then
            let n_the_end=${n_the_end}+1
            echo -e "${ffout}: ${polyname} has finished."
         else
            let n_unknown=${n_unknown}+1
            echo -e "${ffout}: ${polyname} status is UNKNOWN."

            #------------------------------------------------------------------------------#
            #      Kill jobs, it must be stalled.                                          #
            #------------------------------------------------------------------------------#
            jobid="enter_loop"
            kcnt=0
            while [[ "${jobid}" != "" ]] && [[ ${kcnt} -lt 10 ]]
            do
               #----- Ensure to eliminate potential duplicates. ---------------------------#
               let kcnt=${kcnt}+1
               jobln=$(${stask} -o "${outform}" | grep -v "CANCEL" | grep -v "COMPLET"     \
                      | head -1)
               jobid=$(echo ${jobln} | awk '{print $1}')
               if [[ "${jobid}" != "" ]]
               then
                  echo "   - Cancel Job ${jobid} (${taskname})."
                  scancel ${jobid}
               fi
               #---------------------------------------------------------------------------#
            done
            if [[ ${kcnt} -ge 10 ]]
            then
               #------ The script is bogus. -----------------------------------------------#
               echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
               echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
               echo "  JOB NUMBER IS LIKELY INCORRECT IN STASK. STOPPING THE SCRIPT.    "
               echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
               echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
               echo "   Output of stask command:"
               ${stask} -o "${outform}"
               echo "   Resulting JobId:"
               jobln=$(${stask} -o "${outform}" | grep -v "CANCEL" | grep -v "COMPLET"     \
                      | head -1)
               jobid=$(echo ${jobln} | awk '{print $1}')
               echo " JobLine=${jobid}"
               echo " JobId=${jobid}"
               echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
               exit 99
               #---------------------------------------------------------------------------#
            fi
            #------------------------------------------------------------------------------#


            #----- Re-submit job. ---------------------------------------------------------#
            echo "   - Submit New job (${taskname})."
            callserial="${here}/${polyname}/callserial.sh"
            sbatch  ${callserial}
            #------------------------------------------------------------------------------#
         fi
         #---------------------------------------------------------------------------------#
      else
         let n_pending=${n_pending}+1
         echo -e "${ffout}: ${polyname} is pending."
      fi
      #------------------------------------------------------------------------------------#



      #----- Write polygon check into a single table. -------------------------------------#
      output="${polyname} ${polylon} ${polylat} ${yearh} ${monthh} ${dateh} ${timeh}"
      output="${output} ${stall} ${runt} ${agb} ${bsa} ${lai} ${scb} ${npa}"
      echo ${output} >> ${outcheck}
      #------------------------------------------------------------------------------------#
   done
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      Run the post processing.                                                         #
   #---------------------------------------------------------------------------------------#
   if [[ ${frqpost} -gt 0 ]]
   then
      let xpost=${iter}%${frqpost}
      
      if [[ ${xpost} -eq 0 ]] || [[ ${n_ongoing} -eq 0 ]]
      then
         echo " Run post-processing."
         ${here}/epost.sh


         echo " Compress/clean files."
         sbatch ${here}/last_histo.sh
      fi
   fi
   #---------------------------------------------------------------------------------------#



   #----- Run R to make the status check. -------------------------------------------------#
   if ${plotstatus}
   then
      echo "     + Run the status check."
      Rscript_out="$(dirname ${Rscript_plot})/out_$(basename ${Rscript_plot} .r).txt"
      R CMD BATCH --no-save --no-restore ${Rscript_plot} ${Rscript_out}
   fi
   #---------------------------------------------------------------------------------------#



   #----- Current simulation status. ------------------------------------------------------#
   n_initial=$(grep INITIAL ${outcheck} | wc -l)
   n_history=$(grep HISTORY ${outcheck} | wc -l)
   n_metmiss=$(grep METMISS ${outcheck} | wc -l)
   n_bad_met=$(grep BAD_MET ${outcheck} | wc -l)
   n_crashed=$(grep CRASHED ${outcheck} | wc -l)
   n_hydfail=$(grep HYDFAIL ${outcheck} | wc -l)
   n_stopped=$(grep STOPPED ${outcheck} | wc -l)
   n_extinct=$(grep EXTINCT ${outcheck} | wc -l)
   n_ststate=$(grep STSTATE ${outcheck} | wc -l)
   n_the_end=$(grep THE_END ${outcheck} | wc -l)
   let n_ongoing=${n_polygon}-${n_the_end}-${n_extinct}-${n_ststate}
   /bin/rm -f ${statfile}
   touch ${statfile}
   echo "------- Simulation status. --------------------------------------" >> ${statfile}
   echo " Number of polygons that have never started        : ${n_initial}" >> ${statfile}
   echo " Number of polygons that have partially run        : ${n_running}" >> ${statfile}
   echo " Number of polygons that have stalled              : ${n_stalled}" >> ${statfile}
   echo " Number of polygons that haven't found met drivers : ${n_metmiss}" >> ${statfile}
   echo " Number of polygons that have bad met drivers      : ${n_bad_met}" >> ${statfile}
   echo " Number of polygons that have crashed (RK4)        : ${n_crashed}" >> ${statfile}
   echo " Number of polygons that have crashed (Hydrodyn.)  : ${n_hydfail}" >> ${statfile}
   echo " Number of polygons that have mysteriously stopped : ${n_stopped}" >> ${statfile}
   echo " Number of polygons that became desert             : ${n_extinct}" >> ${statfile}
   echo " Number of polygons that have reached steady state : ${n_ststate}" >> ${statfile}
   echo " Number of polygons that have reached the end      : ${n_the_end}" >> ${statfile}
   echo "-----------------------------------------------------------------" >> ${statfile}
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Build the e-mail body.                                                            #
   #---------------------------------------------------------------------------------------#
   echo "     + Build e-mail..."
   /bin/rm -f ${emailbody}
   touch ${emailbody}
   cat ${headfile}  >> ${emailbody}
   echo " "         >> ${emailbody}
   cat ${statfile}  >> ${emailbody}
   echo " "         >> ${emailbody}
   cat ${tailfile}  >> ${emailbody}
   echo " "         >> ${emailbody}
   #----- Check whether to append some plots. ---------------------------------------------#
   if ${plotstatus}
   then
      attach=""
      for fichier in ${R_figlist}
      do
         if [[ -s ${fichier} ]]
         then
            attach="${attach} -a ${fichier}"
         fi
      done
   else
      attach=""
   fi
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Create the subject for the e-mail.                                                #
   #---------------------------------------------------------------------------------------#
   when=$(date +'%d %B %Y - %R %Z')
   elapsed=$(date +'%s')
   subject="${runtitle} run status as of ${when}"
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Check whether to send the e-mail or not.                                          #
   #---------------------------------------------------------------------------------------#
   if ${emailevery}
   then
     #----- Always send e-mail. -----------------------------------------------------------#
     sendemail=true
     #-------------------------------------------------------------------------------------#
   elif [[ -s ${emailday} ]]
   then
      #----- E-mail has been sent.  Check whether it is time to send again. ---------------#
      elast=$(cat ${emailday})
      let ediff=${elapsed}-${elast}
      if [[ ${ediff} -ge ${frqemail} ]]
      then
         #----- Time to send another e-mail. ----------------------------------------------#
         sendemail=true
         #---------------------------------------------------------------------------------#
      else
         #----- E-mail has been sent recently.  Don't send another one now. ---------------#
         sendemail=false
         #---------------------------------------------------------------------------------#
      fi
      #------------------------------------------------------------------------------------#
   else
      #----- E-mail has not been sent yet.  Send it. --------------------------------------#
      sendemail=true
      #------------------------------------------------------------------------------------#
   fi
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Send/skip the e-mail.                                                            #
   #---------------------------------------------------------------------------------------#
   if ${sendemail}
   then
      echo ${elapsed} > ${emailday}
      echo "     + Send e-mail."
      ${mailprog} -s "${subject}" ${attach} ${recipient} < ${emailbody}
   else
      echo "     + Skip e-mail."
   fi
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #       In case any simulation is still on the works, take a break before checking      #
   #   again.                                                                              #
   #---------------------------------------------------------------------------------------#
   if [[ ${n_ongoing} -gt 0 ]]
   then
      echo "     + Take a break before checking again."
      let nsl=${wait_minutes}+1
      while [[ ${nsl} -gt 1 ]]
      do
         let nsl=${nsl}-1
         case ${nsl} in
         1)
            echo "       - ${nsl} minute before next iteration."
            ;;
         *)
            echo "       - ${nsl} minutes before next iteration."
            ;;
         esac
         sleep 1m
      done
   else
      echo "     + All simulations have finished."
   fi
   echo "==================================================================="
   echo "==================================================================="
   echo ""
   #---------------------------------------------------------------------------------------#
done
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Run backup one last time before we end the execution.  Mind that the script may be   #
# transferring the files.                                                                  #
#------------------------------------------------------------------------------------------#
nwait=0
while [[ -s ${translock} ]]
do
   let nwait=${nwait}+1
   echo " + Attempt ${nwait}, script transfer is running.  Wait 10 minutes and try again."
   sleep 10m
done
#---- This time we don't run in background.  ----------------------------------------------#
echo " Backup files one last time."
${transfer}
#------------------------------------------------------------------------------------------#



#----- Clean-up stuff. --------------------------------------------------------------------#
echo " Unlock run_sitter.sh."
/bin/rm -f ${here}/run_sitter.lock
echo "==== run_sitter.sh execution ends. ===="
#------------------------------------------------------------------------------------------#
