#!/bin/bash
. ${HOME}/.bashrc

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
here="xxxxxxxxxxxxxxxxxxxxx"
if [ ${here} == "xxxxxxxxxxxxxxxxxxxxx" ]
then
   echo " You must set up variable here before using run_sitter.sh!!!"
   exit 99
fi
#------------------------------------------------------------------------------------------#


#----- Where the output is. ---------------------------------------------------------------#
there=${here}
#------------------------------------------------------------------------------------------#


#----- Path with some utilities for run_sitter.sh (this script). --------------------------#
situtils="${here}/sit_utils"
#------------------------------------------------------------------------------------------#


#----- Path where biomass initialisation files are: ---------------------------------------#
bioinit="/n/home00/mlongo/data/ed2_data/site_bio_data"
biotype=0      # 0 -- "default" setting (isizepft controls default/nounder)
               # 1 -- isizepft controls number of PFTs, whereas iage controls patches.
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#   desc     -- Unique job description for these simulations (it really must be unique).   #
#               Normally it is the path name.                                              #
#   runtitle -- Full name of this simulation, this is used only in the e-mail subject.     #
#------------------------------------------------------------------------------------------#
desc=$(basename ${here})
runtitle="ED-2 simulation"
#------------------------------------------------------------------------------------------#


#----- This is the header with the Sheffield data. ----------------------------------------#
shefhead="SHEF_NCEP_DRIVER_DS314"
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#      Variables that tell whether to check for steady state and deserts.                  #
#------------------------------------------------------------------------------------------#
checksteady="FALSE"   #   Check whether to stop simulations? (TRUE of FALSE, R style)
nyearmin=160          #   Minimum number of years before we check
ststcrit=0.01         #   Maximum change allowed between two cycles
#------------------------------------------------------------------------------------------#



#----- User name, usually set by $(whoami) so you don't need to change it. ----------------#
moi=$(whoami)
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#    Met driver location settings.                                                         #
#    copy2scratch -- Should the met driver be copied to local scratch disks? ("y"|"n")     #
#    metmaindef   -- Source path if we are not going to copy.                              #
#    packdatasrc  -- Source path from where to copy to scratch.                            #
#------------------------------------------------------------------------------------------#
copy2scratch="y"
metmaindef="/n/home00/mlongo/data/ed2_data"
packdatasrc="/n/home00/mlongo/data/2scratch"
#------------------------------------------------------------------------------------------#



#----- Path with land use scenarios. ------------------------------------------------------#
lumain="/n/gstore/Labs/moorcroft_lab_protected/mlongo/scenarios"
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#    Main actions.                                                                         #
#    printevery -- when looping through jobs, print something on screen every printevery   #
#                  polygons                                                                #
#    bigloop    -- Go through the big loop (0 - no, 1 - yes)                               #
#    fixqueue   -- 0: perform the full queue management)                                   #
#                  1: simply fix the queues on joborder.                                   #
#------------------------------------------------------------------------------------------#
printevery=1
bigloop=1
fixqueue=0
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#  waitpot -- Should polygons wait until some previous run is done? ("n" will start them   #
#             right away).                                                                 #
#  potveg  -- Path where to look for initialisation files in case waitpot = "y".           #
#------------------------------------------------------------------------------------------#
waitpot="n"
potveg="/n/home00/mlongo/data/ed2_data/restarts_sci_006/potveg"
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
# copyrestart -- Should I copy the latest history file somewhere else ("y"/"N").           #
# restart     -- Path to where to copy restart files in case copyrestart = "y".            #
#------------------------------------------------------------------------------------------#
copyrestart="n"
restart="/n/home00/mlongo/data/ed2_data/restarts_XXX"
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
situation="${situtils}/situation.txt"
#------------------------------------------------------------------------------------------#


#------ Calculator. -----------------------------------------------------------------------#
ccc="${HOME}/util/calc.sh"  # Calculator
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#    Maximum number of jobs in the queue that I am allowed to run in each queue.           #
#------------------------------------------------------------------------------------------#
#------ Queue: moorcroft_6100. ------------------------------------------------------------#
m61max=60   # moorcroft_6100
amdmax=60   # moorcroft_amd
uremax=60   # unrestricted
#------------------------------------------------------------------------------------------#



#----- Requested memory and time. ---------------------------------------------------------#
memory=2048                 # Requested memory per cpu
runtime_own="30-00:00:00"   # Runtime for lab owned nodes
runtime_gen="7-00:00:00"    # Runtime for general nodes
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#    E-mail options (normally only the first three options may require changes.            #
#    email1day    -- Should I e-mail once a day (1) or every time I run (0)?               #
#    recipient    -- To which e-mail should I send the update?                             #
#    mailprog     -- Which e-mail program to use? mutt works best, I installed in my util  #
#                    directory.  Give a try, if it doesn't work you may need to install    #
#                    locally on your directory.                                            #
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
email1day=1
recipient="xxxxxx@xxxx.com"
mailprog="/n/home00/mlongo/util/mutt"
plotstatus=1
Rscript_plot="${situtils}/plot.status.r"
R_figlist="${situtils}/stt_t+000.png
           ${situtils}/agb_t+000.png
           ${situtils}/bsa_t+000.png
           ${situtils}/lai_t+000.png
           ${situtils}/scb_t+000.png
           ${situtils}/npa_t+000.png"
emailbody="${situtils}/email.txt"
headfile="${situtils}/head.txt"
tailfile="${situtils}/tail.txt"
recefile="${situtils}/rece.txt"
statfile="${situtils}/stat.txt"
queuefile="${situtils}/queue.txt"
tablefile="${situtils}/table_queue.txt"
pendfile="${here}/pending.txt"
emailday="${here}/emailday.txt"
if [ ${recipient} == "xxxxxx@xxxx.com" ]
then
   echo "You must specify the e-mail (variable recipient)"
   stop 92
fi
#------------------------------------------------------------------------------------------#



#----- Executable name. -------------------------------------------------------------------#
execname="ed_2.1-opt"
initrc="${HOME}/.bashrc"
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
#   Check whether the executable is copied.  If not, let the user know and stop the        #
# script.                                                                                  #
#------------------------------------------------------------------------------------------#
if [ ! -s ${here}/executable/${execname} ]
then
   echo "Executable file : ${execname} is not in the executable directory"
   echo "Copy the executable to the file before running this script!"
   exit 99
fi
#------------------------------------------------------------------------------------------#





#----- Check whether run_sitter.sh is still running or not.  If it is, exit. --------------#
if [ -s ${here}/run_sitter.lock ]
then
   exit
else
   echo "I am going to check your run." > ${here}/run_sitter.lock
fi
#------------------------------------------------------------------------------------------#


#----- Set the main path for the site, pseudo past and Sheffield met drivers. -------------#
if [ ${copy2scratch} == "y" -o ${copy2scratch} == "Y" ]
then
   metmain="/scratch/mlongo"
else
   metmain=${metmaindef}
fi
#------------------------------------------------------------------------------------------#



#----- Determine the number of polygons to check. -----------------------------------------#
let npolys=$(wc -l ${joborder} | awk '{print $1 }')-3
echo "Number of polygons: ${npolys}..."
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#     Reset most counters.                                                                 #
#------------------------------------------------------------------------------------------#
ff=0           # Polygon counter.
nstart=0       # Number of new polygons that initialised
nmetmiss=0     # Number of new polygons that crashed due to missing met driver
nbad_met=0     # Number of new polygons that crashed due to bad met driver
nsigsegv=0     # Number of new polygons that crashed due to segmentation violation
nstopped=0     # Number of new polygons that crashed due to unknown reason
ncrashed=0     # Number of new polygons that crashed due to numeric instability
newsuspend=0   # Number of new polygons that have been suspended by priority people
newunknown=0   # Number of new polygons whose status is unknown
newststate=0   # Number of new polygons at steady state
newextinct=0   # Number of new extinct polygons
nresubmit=0    # Number of moorcroft2c polygons that were submitted again.
deathrow=""    # List of polygons to be killed (they went extinct or reached steady state)
#------------------------------------------------------------------------------------------#



if [ ${bigloop} -ne 0 ]
then
   #---------------------------------------------------------------------------------------#
   #     Eliminate the latest check, and move the most recent check to the latest check,   #
   # and start the new check.                                                              #
   #---------------------------------------------------------------------------------------#
   rm -f ${lastcheck}
   mv ${outcheck} ${lastcheck}
   touch ${outcheck}
   rm -f ${situation}
   touch ${situation}
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Loop over all polygons.                                                           #
   #---------------------------------------------------------------------------------------#
   while [ ${ff} -lt ${npolys} ]
   do
      let ff=${ff}+1
      let line=${ff}+3


      #------------------------------------------------------------------------------------#
      #    Format count.                                                                   #
      #------------------------------------------------------------------------------------#
      if   [ ${npolys} -ge 10   ] && [ ${npolys} -lt 100   ]
      then
         ffout=$(printf '%2.2i' ${ff})
      elif [ ${npolys} -ge 100  ] && [ ${npolys} -lt 1000  ]
      then
         ffout=$(printf '%2.2i' ${ff})
      elif [ ${npolys} -ge 100  ] && [ ${npolys} -lt 10000 ]
      then
         ffout=$(printf '%2.2i' ${ff})
      else
         ffout=${ff}
      fi
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #      Read the polyth line of the polygon list.  There must be smarter ways of do-  #
      # ing this, but this works.  Here we obtain the polygon name, and its longitude and  #
      # latitude.                                                                          #
      #------------------------------------------------------------------------------------#
      oi=$(head -${line} ${lonlat} | tail -1)
      polyname=$(echo ${oi}     | awk '{print $1 }')
      polyiata=$(echo ${oi}     | awk '{print $2 }')
      polylon=$(echo ${oi}      | awk '{print $3 }')
      polylat=$(echo ${oi}      | awk '{print $4 }')
      yeara=$(echo ${oi}        | awk '{print $5 }')
      montha=$(echo ${oi}       | awk '{print $6 }')
      datea=$(echo ${oi}        | awk '{print $7 }')
      timea=$(echo ${oi}        | awk '{print $8 }')
      yearz=$(echo ${oi}        | awk '{print $9 }')
      monthz=$(echo ${oi}       | awk '{print $10}')
      datez=$(echo ${oi}        | awk '{print $11}')
      timez=$(echo ${oi}        | awk '{print $12}')
      initmode=$(echo ${oi}     | awk '{print $13}')
      iscenario=$(echo ${oi}    | awk '{print $14}')
      isizepft=$(echo ${oi}     | awk '{print $15}')
      iage=$(echo ${oi}         | awk '{print $16}')
      polyisoil=$(echo ${oi}    | awk '{print $17}')
      polyntext=$(echo ${oi}    | awk '{print $18}')
      polysand=$(echo ${oi}     | awk '{print $19}')
      polyclay=$(echo ${oi}     | awk '{print $20}')
      polydepth=$(echo ${oi}    | awk '{print $21}')
      polysoilbc=$(echo ${oi}   | awk '{print $22}')
      polysldrain=$(echo ${oi}  | awk '{print $23}')
      polycol=$(echo ${oi}      | awk '{print $24}')
      slzres=$(echo ${oi}       | awk '{print $25}')
      queue=$(echo ${oi}        | awk '{print $26}')
      metdriver=$(echo ${oi}    | awk '{print $27}')
      dtlsm=$(echo ${oi}        | awk '{print $28}')
      vmfactc3=$(echo ${oi}     | awk '{print $29}')
      vmfactc4=$(echo ${oi}     | awk '{print $30}')
      mphototrc3=$(echo ${oi}   | awk '{print $31}')
      mphototec3=$(echo ${oi}   | awk '{print $32}')
      mphotoc4=$(echo ${oi}     | awk '{print $33}')
      bphotoblc3=$(echo ${oi}   | awk '{print $34}')
      bphotonlc3=$(echo ${oi}   | awk '{print $35}')
      bphotoc4=$(echo ${oi}     | awk '{print $36}')
      kwgrass=$(echo ${oi}      | awk '{print $37}')
      kwtree=$(echo ${oi}       | awk '{print $38}')
      gammac3=$(echo ${oi}      | awk '{print $39}')
      gammac4=$(echo ${oi}      | awk '{print $40}')
      d0grass=$(echo ${oi}      | awk '{print $41}')
      d0tree=$(echo ${oi}       | awk '{print $42}')
      alphac3=$(echo ${oi}      | awk '{print $43}')
      alphac4=$(echo ${oi}      | awk '{print $44}')
      klowco2=$(echo ${oi}      | awk '{print $45}')
      decomp=$(echo ${oi}       | awk '{print $46}')
      rrffact=$(echo ${oi}      | awk '{print $47}')
      growthresp=$(echo ${oi}   | awk '{print $48}')
      lwidthgrass=$(echo ${oi}  | awk '{print $49}')
      lwidthbltree=$(echo ${oi} | awk '{print $50}')
      lwidthnltree=$(echo ${oi} | awk '{print $51}')
      q10c3=$(echo ${oi}        | awk '{print $52}')
      q10c4=$(echo ${oi}        | awk '{print $53}')
      h2olimit=$(echo ${oi}     | awk '{print $54}')
      imortscheme=$(echo ${oi}  | awk '{print $55}')
      ddmortconst=$(echo ${oi}  | awk '{print $56}')
      cbrscheme=$(echo ${oi}    | awk '{print $57}')
      isfclyrm=$(echo ${oi}     | awk '{print $58}')
      icanturb=$(echo ${oi}     | awk '{print $59}')
      ubmin=$(echo ${oi}        | awk '{print $60}')
      ugbmin=$(echo ${oi}       | awk '{print $61}')
      ustmin=$(echo ${oi}       | awk '{print $62}')
      gamm=$(echo ${oi}         | awk '{print $63}')
      gamh=$(echo ${oi}         | awk '{print $64}')
      tprandtl=$(echo ${oi}     | awk '{print $65}')
      ribmax=$(echo ${oi}       | awk '{print $66}')
      atmco2=$(echo ${oi}       | awk '{print $67}')
      thcrit=$(echo ${oi}       | awk '{print $68}')
      smfire=$(echo ${oi}       | awk '{print $69}')
      ifire=$(echo ${oi}        | awk '{print $70}')
      fireparm=$(echo ${oi}     | awk '{print $71}')
      ipercol=$(echo ${oi}      | awk '{print $72}')
      runoff=$(echo ${oi}       | awk '{print $73}')
      imetrad=$(echo ${oi}      | awk '{print $74}')
      ibranch=$(echo ${oi}      | awk '{print $75}')
      icanrad=$(echo ${oi}      | awk '{print $76}')
      crown=$(echo   ${oi}      | awk '{print $77}')
      ltransvis=$(echo ${oi}    | awk '{print $78}')
      lreflectvis=$(echo ${oi}  | awk '{print $79}')
      ltransnir=$(echo ${oi}    | awk '{print $80}')
      lreflectnir=$(echo ${oi}  | awk '{print $81}')
      orienttree=$(echo ${oi}   | awk '{print $82}')
      orientgrass=$(echo ${oi}  | awk '{print $83}')
      clumptree=$(echo ${oi}    | awk '{print $84}')
      clumpgrass=$(echo ${oi}   | awk '{print $85}')
      ivegtdyn=$(echo ${oi}     | awk '{print $86}')
      igndvap=$(echo ${oi}      | awk '{print $87}')
      iphen=$(echo ${oi}        | awk '{print $88}')
      iallom=$(echo ${oi}       | awk '{print $89}')
      ibigleaf=$(echo ${oi}     | awk '{print $90}')
      irepro=$(echo ${oi}       | awk '{print $91}')
      treefall=$(echo ${oi}     | awk '{print $92}')
      ianthdisturb=$(echo ${oi} | awk '{print $93}')
      ianthdataset=$(echo ${oi} | awk '{print $94}')
      #------------------------------------------------------------------------------------#



      #----- Check whether the directories exist or not. ----------------------------------#
      if [ ! -s ${here}/${polyname} ]
      then
         echo "Order: ${ff}, creating directory for polygon ${polyname}..."
         #----- Copy the Template directory to a unique polygon directory. -------------------#
         oldtol=""
         cp -r ${here}/Template ${here}/${polyname}
      fi
      #---------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Define whether we use the met cycle to define the first and last year, or the  #
      # default year.                                                                      #
      #------------------------------------------------------------------------------------#
      if [ ${yeara} -eq 0 ]
      then
         thisyeara=${metcyc1}
      else
         thisyeara=${yeara}
      fi
      if [ ${yearz} -eq 0 ]
      then
         thisyearz=${metcycf}
      else
         thisyearz=${yearz}
      fi
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Set the land use data set.                                                     #
      #------------------------------------------------------------------------------------#
      case ${ianthdataset} in
      glu-331)
         case ${polyiata} in
         tzi|zmh|nqn|hvd|wch|tqh)
            ludatabase="${lumain}/glu/outglu/glu-"
            ;;
         *)
            ludatabase="${lumain}/glu-3.3.1/one/glu-3.3.1-"
            ;;
         esac
         ;;
      glu-sa1)
         case ${polyiata} in
         tzi|zmh|nqn|hvd|wch|tqh)
            ludatabase="${lumain}/glu/outglu/glu-"
            ;;
         *)
            ludatabase="${lumain}/glu-3.3.1+sa1.bau/one/glu-3.3.1+sa1.bau-"
            ;;
         esac
         ;;
      glu-sag)
         case ${polyiata} in
         tzi|zmh|nqn|hvd|wch|tqh)
            ludatabase="${lumain}/glu/outglu/glu-"
            ;;
         *)
            ludatabase="${lumain}/glu-3.3.1+sa1.gov/one/glu-3.3.1+sa1.gov-"
            ;;
         esac
         ;;
      glu-sa2)
         case ${polyiata} in
         tzi|zmh|nqn|hvd|wch|tqh)
            ludatabase="${lumain}/glu/outglu/glu-"
            ;;
         *)
            ludatabase="${lumain}/glu-3.3.1+sa2.bau/one/glu-3.3.1+sa2.bau-"
            ;;
         esac
         ;;
      lurcp26)
         ludatabase="${lumain}/luha-v1/luh-1.1+rcp26_image/half/luh-1.1+rcp26_image-"
         ;;
      lurcp45)
         ludatabase="${lumain}/luha-v1/luh-1.1+rcp45_minicam/half/luh-1.1+rcp45_minicam-"
         ;;
      lurcp60)
         ludatabase="${lumain}/luha-v1/luh-1.1+rcp60_aim/half/luh-1.1+rcp60_aim-"
         ;;
      lurcp85)
         ludatabase="${lumain}/luha-v1/luh-1.1+rcp85_message/half/luh-1.1+rcp85_message-"
         ;;
      *)
         #---------------------------------------------------------------------------------#
         #     Stop the script if anthropogenic dataset is invalid and this is a           #
         # simulation with anthropogenic disturbance.                                      #
         #---------------------------------------------------------------------------------#
         if [ ${ianthdisturb} -eq 1 ]
         then
            echo " Polygon:       ${polyname}"
            echo " IATA:          ${polyiata}"
            echo " IANTH_DATASET: ${iage}"
            echo "Invalid anthropogenic disturbance data set!"
            exit 53
         fi
         #---------------------------------------------------------------------------------#
         ;;
      esac
      #------------------------------------------------------------------------------------#




      #------------------------------------------------------------------------------------#
      #     Determine which soil profile to use.  We have eight categories (A-H), and the  #
      # soil resolution.                                                                   #
      #------------------------------------------------------------------------------------#
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
            polyslz1="-4.500,-3.967,-3.467,-3.000,-2.565,-2.164,-1.797,-1.462,-1.162,-0.895,"
            polyslz2="-0.662,-0.464,-0.300,-0.171,-0.077,-0.020"
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
      1)
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
            ;;
         esac
         ;;
      esac
      #------------------------------------------------------------------------------------#





      #------------------------------------------------------------------------------------#
      #     Determine which PFTs to use based on the "iata" code and isizepft.             #
      #------------------------------------------------------------------------------------#
      case ${isizepft} in
      0|1|5)
         case ${polyiata} in
         tzi|zmh|nqn)
            pfts="6,7,9,10,11,16,17"
            crop=16
            plantation=17
            ;;
         hvd|wch|tqh)
            pfts="6,8,9,10,11,16,17"
            crop=16
            plantation=17
            ;;
         asu|cnf|bnu|cwb|erm|iqq|ipv|mgf|rao|sla|zpe|kna|sfn)
            pfts="1,2,3,4,16,17"
            crop=16
            plantation=17
            ;;
         fns*)
            pfts="1,16"
            crop=1
            plantation=17
            ;;
         s77*)
            pfts="1,16"
            crop=16
            plantation=17
            ;;
         *)
            pfts="1,2,3,4,16"
            crop=1
            plantation=3
            ;;
         esac
         ;;
      2)
         case ${polyiata} in
         tzi|zmh|nqn|hvd|wch|tqh)
            pfts="10,16"
            crop=16
            plantation=17
            ;;
         fns*|s77*)
            pfts="1"
            crop=1
            plantation=17
            ;;
         *)
            pfts="1,3"
            crop=1
            plantation=3
            ;;
         esac
         ;;
      esac
      #------------------------------------------------------------------------------------#





      #------------------------------------------------------------------------------------#
      #     Determine the census structure.                                                #
      #------------------------------------------------------------------------------------#
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
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #    Find the biometric files.  This has been checked in spawn_poly.sh so they are   #
      # correct.  Add a dummy name in case this is not supposed to be a biomass            #
      # initialisation run.                                                                #
      #------------------------------------------------------------------------------------#
      case ${biotype} in
      0)
         #----- isizepft controls everything, and iage is ignored. ------------------------#
         case ${isizepft} in
         0)
            #----- Frankeinstein's under storey. ------------------------------------------#
            thisbiomin="${bioinit}/${polyiata}_default."
            ;;
         1)
            #----- No under storey. -------------------------------------------------------#
            thisbiomin="${bioinit}/${polyiata}_nounder."
            ;;
         2)
            #----- Same as default, but with only one grass and one tree. -----------------#
            thisbiomin="${bioinit}/${polyiata}_twopft."
            ;;
         *)
            #----- Invalid option. Stop the script. ---------------------------------------#
            echo ' Polygon:  '${polyname}
            echo ' IATA:     '${polyiata}
            echo ' ISIZEPFT: '${isizepft}
            echo 'This IATA cannot be used by biomass initialisation with this ISIZEPFT!'
            exit 57
            ;;
         esac
         #---------------------------------------------------------------------------------#
         ;;

      1)
         #---------------------------------------------------------------------------------#
         #    'isizepft' controls how many PFTs to use.                                    #
         #---------------------------------------------------------------------------------#
         case ${isizepft} in 
         0|5)
            pftname="pft05"
            ;;
         2)
            pftname="pft02"
            ;;
         esac
         #---------------------------------------------------------------------------------#

         #---------------------------------------------------------------------------------#
         #     'iage' controls how many patches to use.                                    #
         #---------------------------------------------------------------------------------#
         case ${iage} in
         1)
            agename="age01"
            ;;
         *)
            agename="age00"
            ;;
         esac
         #---------------------------------------------------------------------------------#




         #---------------------------------------------------------------------------------#
         #      Check whether the site has the PFT and age structure.                      #
         #---------------------------------------------------------------------------------#
         case ${polyiata} in
         hvd|s77|fns|cau|and|par|tap|dcm)
            thisbiomin="${bioinit}/${polyiata}_default."
            ;;
         cax|s67|s83|m34|gyf|pdg|rja|pnz|ban)
            thisbiomin="${bioinit}/${polyiata}_${pftname}+${agename}."
            ;;
         *)
            echo " Polygon:  ${polyname}"
            echo " IATA:     ${polyiata}"
            echo " IAGE:     ${iage}"
            echo " ISIZEPFT: ${isizepft}"
            echo "This IATA cannot be used by biomass initialisation with this ISIZEPFT!"
            exit 59
            ;;
         esac
         #---------------------------------------------------------------------------------#
         ;;
      *)
         echo " Invalid biotype:  ${biotype}"
         exit 58
         ;;
      esac
      #------------------------------------------------------------------------------------#





      #------------------------------------------------------------------------------------#
      #     Determine the scenario paths.                                                  #
      #------------------------------------------------------------------------------------#
      case ${iscenario} in
      default)
         case ${metdriver} in
         Sheffield)
            #----- Sheffield. -------------------------------------------------------------#
            scentype="sheffield"
            iscenario="sheffield"
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
      Manaus_Km34)
         metdriverdb="${fullscen}/Manaus_Km34/Manaus_Km34_HEADER"
         metcyc1=1999
         metcycf=2006
         imetavg=1
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
         metcycf=2012
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
      Tonzi)
         metdriverdb="${fullscen}/Tonzi/Tonzi_HEADER"
         metcyc1=2000
         metcycf=2010
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
      if [ ${iscenario} != "default"   ] && [ ${iscenario} != "eft"       ] && 
         [ ${iscenario} != "shr"       ] && [ ${iscenario} != "sheffield" ]
      then
         metcyc1=1972
         metcycf=2012
         imetavg=1
      fi
      #------------------------------------------------------------------------------------#



      #----- Define the job name. ---------------------------------------------------------#
      jobname="${desc}-${polyname}"
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
      #    Retrie ve the information on statusrun.txt because it may tell whether the      #
      # run was in steady state or all PFTs went extinct.  These will never start          #
      # again.                                                                             #
      #------------------------------------------------------------------------------------#
      statrun=${here}/${polyname}/statusrun.txt
      if [ -s ${statrun} ]
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




      #------------------------------------------------------------------------------------#
      #    We only check those polygons that may be running, so the check doesn't take     #
      # too long.  Once a polygon has reached a steady state / gone extinct / finished,    #
      # then its status should remain the same.                                            #
      #------------------------------------------------------------------------------------#
      if [ ${runt} != "EXTINCT" ] && [ ${runt} != "STSTATE" ] && 
         [ ${runt} != "THE_END" ]
      then
         #----- Call R to check status. ---------------------------------------------------#
         whichrun=${here}/${polyname}/whichrun.r
         outwhich=${here}/${polyname}/outwhich.txt
         R CMD BATCH --no-save --no-restore ${whichrun} ${outwhich}
         year=$(cat ${statrun}  | awk '{print  $2}')
         month=$(cat ${statrun} | awk '{print  $3}')
         day=$(cat ${statrun}   | awk '{print  $4}')
         hhmm=$(cat ${statrun}  | awk '{print  $5}')
         runt=$(cat ${statrun}  | awk '{print  $6}')
         agb=$(cat ${statrun}   | awk '{print  $7}')
         bsa=$(cat ${statrun}   | awk '{print  $8}')
         lai=$(cat ${statrun}   | awk '{print  $9}')
         scb=$(cat ${statrun}   | awk '{print $10}')
         npa=$(cat ${statrun}   | awk '{print $11}')
      fi
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #      Variables with the simulation output.                                         #
      #------------------------------------------------------------------------------------#
      serial_out="${here}/${polyname}/serial_out.out"
      serial_err="${here}/${polyname}/serial_out.err"
      serial_lsf="${here}/${polyname}/serial_lsf.out"
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #      Check whether the job has been suspended because a priority job by a          #
      # priority user have priority over yours, second-class user.  If so, kill the job    #
      # and re-submit, it may be faster.                                                   #
      #------------------------------------------------------------------------------------#
      jobstat=$(squeue -h -u ${moi} -n ${jobname} -o "%T")
      if [ "x${jobstat}" == "x" ]
      then
         jobstat="NOT_FOUND"
      fi
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #    In case the polygon has crashed, we may need to re-submit...                    #
      #------------------------------------------------------------------------------------#
      case ${jobstat} in
      SUSPENDED)
         #----- Kill the job to resubmit it. ----------------------------------------------#
         echo "${ffout}: ${polyname} HAS BEEN SUSPENDED... <========="
         deathrow="${deathrow} ${jobname}"
         let newsuspend=${newsuspend}+1
         ;;
         #---------------------------------------------------------------------------------#
      UNKNOWN)
         #----- Kill the job to resubmit it. ----------------------------------------------#
         echo "${ffout}: ${polyname} status is UNKNOWN..."
         deathrow="${deathrow} ${jobname}"
         let newunknown=${newunknown}+1
         ;;
         #---------------------------------------------------------------------------------#
      *)
         #---------------------------------------------------------------------------------#
         #     Either the job is still running or it has stopped.  Use runt to decide what #
         # to do.                                                                          #
         #---------------------------------------------------------------------------------#
         if [ ${runt} == "CRASHED" ] || [ ${runt} == "STOPPED" ] ||
            [ ${runt} == "METMISS" ] || [ ${runt} == "SIGSEGV" ] ||
            [ ${runt} == "BAD_MET" ]
         then
            #------------------------------------------------------------------------------#
            #     If serial_out exists, then it must be because the simulation is not      #
            # running.  We must re-submit.                                                 #
            #------------------------------------------------------------------------------#
            if [ -s ${serial_out} ]
            then
               #---------------------------------------------------------------------------#
               #     Find out the tolerance used in the unfortunate run, then submit the   #
               # job with a tolerance 10 times more strict.                                #
               #---------------------------------------------------------------------------#
               ED2IN="${here}/${polyname}/ED2IN"
               #----- Save the tolerance before we delete the file. -----------------------#
               toler=$(grep NL%RK4_TOLERANCE ${ED2IN} | awk '{print $3}')
               if [ "x${toler}" == "xmytoler" -o "x${toler}" == "x" ]
               then
                  toler=0.01
               fi # [ "x${toler}" == "xmytoler" -o "x${toler}" == "x" ]
               #---------------------------------------------------------------------------#


               if [ ${runt} == "CRASHED" ]
               then
                  echo "${ffout}: ${polyname} HAS CRASHED (RK4 PROBLEM)... <==========="
                  let ncrashed=${ncrashed}+1
                  #----- Make the tolerance 10 times smaller. -----------------------------#
                  toler=$(${ccc} ${toler}/10)
                  echo "      - New tolerance = ${toler}"
                  /bin/mv ${serial_out} ${here}/${polyname}/crashed_out.out
                  /bin/mv ${serial_err} ${here}/${polyname}/crashed_out.err
               elif [ ${runt} == "SIGSEGV" ]
               then
                  let nsigsegv=${nsigsegv}+1
                  echo "${ffout}: ${polyname} HAD SEGMENTATION VIOLATION... <==========="
                  /bin/mv ${serial_out} ${here}/${polyname}/sigsegv_out.out
                  /bin/mv ${serial_err} ${here}/${polyname}/sigsegv_out.err
               elif [ ${runt} == "METMISS" ]
               then
                  let nmetmiss=${nmetmiss}+1
                  echo "${ffout}: ${polyname} DID NOT FIND MET DRIVERS... <==========="
                  /bin/mv ${serial_out} ${here}/${polyname}/metmiss_out.out
                  /bin/mv ${serial_err} ${here}/${polyname}/metmiss_out.err
               elif [ ${runt} == "BAD_MET" ]
               then
                  let nbad_met=${nbad_met}+1
                  echo "${ffout}: ${polyname} HAS BAD MET DRIVERS... <==========="
                  /bin/mv ${serial_out} ${here}/${polyname}/bad_met_out.out
                  /bin/mv ${serial_err} ${here}/${polyname}/bad_met_out.err
                  #----- Delete files as they are the met driver is bad. ------------------#
                  /bin/rm -fvr ${here}/${polyname}/analy/*
                  /bin/rm -fvr ${here}/${polyname}/histo/*
               elif [ ${runt} == "STOPPED" ]
               then
                  let nstopped=${nstopped}+1
                  echo "${ffout}: ${polyname} STOPPED (UNKNOWN REASON)... <==========="
                  /bin/mv ${serial_out} ${here}/${polyname}/stopped_out.out
                  /bin/mv ${serial_err} ${here}/${polyname}/stopped_out.err
               fi
               #---------------------------------------------------------------------------#

               if [ ${toler} != ".0000100" ] && [ ${toler} != ".0000010" ] && 
                  [ ${toler} != ".0000001" ] && [ ${toler} != "0" ]        && 
                  [ ${runt} != "BAD_MET" ]
               then

                  #----- Re-generate ED2IN. -----------------------------------------------#
                  rm -f ${ED2IN} 
                  rm -f ${here}/${polyname}/serial_lsf.out

                  #----- Copy the Template to the right directory. ------------------------#
                  cp ${here}/Template/ED2IN ${ED2IN}


                  #----- If last history is the first year, switch it to initial. ---------#
                  if [ ${year} -eq ${yeara} ] || [ ${runt} == "INITIAL" ]
                  then
                     runflag="INITIAL"
                  else
                     runflag="HISTORY"
                  fi
                  #------------------------------------------------------------------------#


                  #----- Check whether to use SFILIN as restart or history. ---------------#
                  if [ ${runflag} == "INITIAL" ] && [ ${initmode} -eq 6 ]
                  then
                     thissfilin=${thisbiomin}
                  elif [ ${runflag} == "INITIAL" ] && [ ${initmode} -eq 5 ]
                  then
                     thissfilin=${restart}
                  else
                     thissfilin=${there}/${polyname}/histo/${polyname}
                  fi
                  #------------------------------------------------------------------------#


                  #----- Update the polygon information on ED2IN. -------------------------#
                  sed -i s@paththere@${there}@g                ${ED2IN}
                  sed -i s@myyeara@${thisyeara}@g              ${ED2IN}
                  sed -i s@mymontha@${montha}@g                ${ED2IN}
                  sed -i s@mydatea@${datea}@g                  ${ED2IN}
                  sed -i s@mytimea@${timea}@g                  ${ED2IN}
                  sed -i s@myyearz@${thisyearz}@g              ${ED2IN}
                  sed -i s@mymonthz@${monthz}@g                ${ED2IN}
                  sed -i s@mydatez@${datez}@g                  ${ED2IN}
                  sed -i s@mytimez@${timez}@g                  ${ED2IN}
                  sed -i s@mydtlsm@${dtlsm}@g                  ${ED2IN}
                  sed -i s@thispoly@${polyname}@g              ${ED2IN}
                  sed -i s@plonflag@${polylon}@g               ${ED2IN}
                  sed -i s@platflag@${polylat}@g               ${ED2IN}
                  sed -i s@timehhhh@${timeh}@g                 ${ED2IN}
                  sed -i s@datehhhh@${dateh}@g                 ${ED2IN}
                  sed -i s@monthhhh@${monthh}@g                ${ED2IN}
                  sed -i s@yearhhhh@${yearh}@g                 ${ED2IN}
                  sed -i s@myinitmode@${initmode}@g            ${ED2IN}
                  sed -i s@mysfilin@${thissfilin}@g            ${ED2IN}
                  sed -i s@mytrees@${pfts}@g                   ${ED2IN}
                  sed -i s@mycrop@${crop}@g                    ${ED2IN}
                  sed -i s@myplantation@${plantation}@g        ${ED2IN}
                  sed -i s@myiphen@${iphen}@g                  ${ED2IN}
                  sed -i s@myallom@${iallom}@g                 ${ED2IN}
                  sed -i s@myisoilflg@${polyisoil}@g           ${ED2IN}
                  sed -i s@mynslcon@${polyntext}@g             ${ED2IN}
                  sed -i s@myslxsand@${polysand}@g             ${ED2IN}
                  sed -i s@myslxclay@${polyclay}@g             ${ED2IN}
                  sed -i s@mysoilbc@${polysoilbc}@g            ${ED2IN}
                  sed -i s@mysldrain@${polysldrain}@g          ${ED2IN}
                  sed -i s@mysoilcol@${polycol}@g              ${ED2IN}
                  sed -i s@mynzg@${polynzg}@g                  ${ED2IN}
                  sed -i s@mymetdriverdb@${metdriverdb}@g      ${ED2IN}
                  sed -i s@mymetcyc1@${metcyc1}@g              ${ED2IN}
                  sed -i s@mymetcycf@${metcycf}@g              ${ED2IN}
                  sed -i s@mytoler@${toler}@g                  ${ED2IN}
                  sed -i s@RUNFLAG@${runflag}@g                ${ED2IN}
                  sed -i s@myvmfactc3@${vmfactc3}@g            ${ED2IN}
                  sed -i s@myvmfactc4@${vmfactc4}@g            ${ED2IN}
                  sed -i s@mymphototrc3@${mphototrc3}@g        ${ED2IN}
                  sed -i s@mymphototec3@${mphototec3}@g        ${ED2IN}
                  sed -i s@mymphotoc4@${mphotoc4}@g            ${ED2IN}
                  sed -i s@mybphotoblc3@${bphotoblc3}@g        ${ED2IN}
                  sed -i s@mybphotonlc3@${bphotonlc3}@g        ${ED2IN}
                  sed -i s@mybphotoc4@${bphotoc4}@g            ${ED2IN}
                  sed -i s@mykwgrass@${kwgrass}@g              ${ED2IN}
                  sed -i s@mykwtree@${kwtree}@g                ${ED2IN}
                  sed -i s@mygammac3@${gammac3}@g              ${ED2IN}
                  sed -i s@mygammac4@${gammac4}@g              ${ED2IN}
                  sed -i s@myd0grass@${d0grass}@g              ${ED2IN}
                  sed -i s@myd0tree@${d0tree}@g                ${ED2IN}
                  sed -i s@myalphac3@${alphac3}@g              ${ED2IN}
                  sed -i s@myalphac4@${alphac4}@g              ${ED2IN}
                  sed -i s@myklowco2@${klowco2}@g              ${ED2IN}
                  sed -i s@mydecomp@${decomp}@g                ${ED2IN}
                  sed -i s@myrrffact@${rrffact}@g              ${ED2IN}
                  sed -i s@mygrowthresp@${growthresp}@g        ${ED2IN}
                  sed -i s@mylwidthgrass@${lwidthgrass}@g      ${ED2IN}
                  sed -i s@mylwidthbltree@${lwidthbltree}@g    ${ED2IN}
                  sed -i s@mylwidthnltree@${lwidthnltree}@g    ${ED2IN}
                  sed -i s@myq10c3@${q10c3}@g                  ${ED2IN}
                  sed -i s@myq10c4@${q10c4}@g                  ${ED2IN}
                  sed -i s@myh2olimit@${h2olimit}@g            ${ED2IN}
                  sed -i s@mymortscheme@${imortscheme}@g       ${ED2IN}
                  sed -i s@myddmortconst@${ddmortconst}@g      ${ED2IN}
                  sed -i s@mycbrscheme@${cbrscheme}@g          ${ED2IN}
                  sed -i s@mysfclyrm@${isfclyrm}@g             ${ED2IN}
                  sed -i s@myicanturb@${icanturb}@g            ${ED2IN}
                  sed -i s@myatmco2@${atmco2}@g                ${ED2IN}
                  sed -i s@mythcrit@${thcrit}@g                ${ED2IN}
                  sed -i s@mysmfire@${smfire}@g                ${ED2IN}
                  sed -i s@myfire@${ifire}@g                   ${ED2IN}
                  sed -i s@myfuel@${fireparm}@g                ${ED2IN}
                  sed -i s@mymetavg@${imetavg}@g               ${ED2IN}
                  sed -i s@mypercol@${ipercol}@g               ${ED2IN}
                  sed -i s@myrunoff@${runoff}@g                ${ED2IN}
                  sed -i s@mymetrad@${imetrad}@g               ${ED2IN}
                  sed -i s@mybranch@${ibranch}@g               ${ED2IN}
                  sed -i s@mycanrad@${icanrad}@g               ${ED2IN}
                  sed -i s@mycrown@${crown}@g                  ${ED2IN}
                  sed -i s@myltransvis@${ltransvis}@g          ${ED2IN}
                  sed -i s@myltransnir@${ltransnir}@g          ${ED2IN}
                  sed -i s@mylreflectvis@${lreflectvis}@g      ${ED2IN}
                  sed -i s@mylreflectnir@${lreflectnir}@g      ${ED2IN}
                  sed -i s@myorienttree@${orienttree}@g        ${ED2IN}
                  sed -i s@myorientgrass@${orientgrass}@g      ${ED2IN}
                  sed -i s@myclumptree@${clumptree}@g          ${ED2IN}
                  sed -i s@myclumpgrass@${clumpgrass}@g        ${ED2IN}
                  sed -i s@myvegtdyn@${ivegtdyn}@g             ${ED2IN}
                  sed -i s@mybigleaf@${ibigleaf}@g             ${ED2IN}
                  sed -i s@myrepro@${irepro}@g                 ${ED2IN}
                  sed -i s@myubmin@${ubmin}@g                  ${ED2IN}
                  sed -i s@myugbmin@${ugbmin}@g                ${ED2IN}
                  sed -i s@myustmin@${ustmin}@g                ${ED2IN}
                  sed -i s@mygamm@${gamm}@g                    ${ED2IN}
                  sed -i s@mygamh@${gamh}@g                    ${ED2IN}
                  sed -i s@mytprandtl@${tprandtl}@g            ${ED2IN}
                  sed -i s@myribmax@${ribmax}@g                ${ED2IN}
                  sed -i s@mygndvap@${igndvap}@g               ${ED2IN}
                  sed -i s@mydtcensus@${dtcensus}@g            ${ED2IN}
                  sed -i s@myyr1stcensus@${yr1stcensus}@g      ${ED2IN}
                  sed -i s@mymon1stcensus@${mon1stcensus}@g    ${ED2IN}
                  sed -i s@myminrecruitdbh@${minrecruitdbh}@g  ${ED2IN}
                  sed -i s@mytreefall@${treefall}@g            ${ED2IN}
                  sed -i s@mymaxpatch@${iage}@g                ${ED2IN}
                  sed -i s@myanthdisturb@${ianthdisturb}@g     ${ED2IN}
                  sed -i s@myludatabase@${ludatabase}@g        ${ED2IN}
                  #------ Soil variables. -------------------------------------------------#
                  sed -i s@myslz1@"${polyslz1}"@g           ${ED2IN}
                  sed -i s@myslz2@"${polyslz2}"@g           ${ED2IN}
                  sed -i s@myslz3@"${polyslz3}"@g           ${ED2IN}
                  sed -i s@myslz4@"${polyslz4}"@g           ${ED2IN}
                  sed -i s@myslz5@"${polyslz5}"@g           ${ED2IN}
                  sed -i s@myslz6@"${polyslz6}"@g           ${ED2IN}
                  sed -i s@myslz7@"${polyslz7}"@g           ${ED2IN}
                  sed -i s@myslmstr1@"${polyslm1}"@g        ${ED2IN}
                  sed -i s@myslmstr2@"${polyslm2}"@g        ${ED2IN}
                  sed -i s@myslmstr3@"${polyslm3}"@g        ${ED2IN}
                  sed -i s@myslmstr4@"${polyslm4}"@g        ${ED2IN}
                  sed -i s@myslmstr5@"${polyslm5}"@g        ${ED2IN}
                  sed -i s@myslmstr6@"${polyslm6}"@g        ${ED2IN}
                  sed -i s@myslmstr7@"${polyslm7}"@g        ${ED2IN}
                  sed -i s@mystgoff1@"${polyslt1}"@g        ${ED2IN}
                  sed -i s@mystgoff2@"${polyslt2}"@g        ${ED2IN}
                  sed -i s@mystgoff3@"${polyslt3}"@g        ${ED2IN}
                  sed -i s@mystgoff4@"${polyslt4}"@g        ${ED2IN}
                  sed -i s@mystgoff5@"${polyslt5}"@g        ${ED2IN}
                  sed -i s@mystgoff6@"${polyslt6}"@g        ${ED2IN}
                  sed -i s@mystgoff7@"${polyslt7}"@g        ${ED2IN}
                  #------------------------------------------------------------------------#






                  #------------------------------------------------------------------------#
                  #      A crashed polygon will inevitably go to the general queue.        #
                  # Re-write the srun.sh file.                                             #
                  #------------------------------------------------------------------------#
                  srun="${here}/${polyname}/srun.sh"
                  callserial="${here}/${polyname}/callserial.sh"
                  #------------------------------------------------------------------------#



                  #----- Get the waiting time used before. --------------------------------#
                  thiswait=$(grep "Sleep time is" ${polyname}/srun.sh | awk '{print $5}')
                  #------------------------------------------------------------------------#


                  #----- Copy the data set again. -----------------------------------------#
                  /bin/rm -f ${srun}
                  /bin/rm -f ${callserial}
                  /bin/cp ${here}/Template/srun.sh       ${srun}
                  /bin/cp ${here}/Template/callserial.sh ${callserial}
                  #------------------------------------------------------------------------#


                  #----- The new queue is by default the general. -------------------------#
                  newqueue="general"
                  runtime=${runtime_gen}
                  #------------------------------------------------------------------------#



                  #----- Change some settings in srun.sh. ---------------------------------#
                  sed -i s@pathhere@${here}@g      ${srun}
                  sed -i s@paththere@${here}@g     ${srun}
                  sed -i s@thispoly@${polyname}@g  ${srun}
                  sed -i s@thisdesc@${desc}@g      ${srun}
                  sed -i s@thisqueue@${newqueue}@g ${srun}
                  sed -i s@zzzzzzzz@${wtime}@g     ${srun}
                  sed -i s@myorder@${ff}@g         ${srun}
                  sed -i s@myinitrc@${initrc}@g    ${srun}
                  sed -i s@thismemory@${memory}@g  ${srun}
                  sed -i s@thistime@${runtime}@g   ${srun}
                  #------------------------------------------------------------------------#




                  #----- Change the callserial.sh file. -----------------------------------#
                  sed -i s@thisroot@${here}@g          ${callserial}
                  sed -i s@thispoly@${polyname}@g      ${callserial}
                  sed -i s@myexec@${execname}@g        ${callserial}
                  sed -i s@myname@${moi}@g             ${callserial}
                  sed -i s@mypackdata@${packdatasrc}@g ${callserial}
                  sed -i s@myscenario@${iscenario}@g   ${callserial}
                  sed -i s@myscenmain@${scentype}@g    ${callserial}
                  sed -i s@myinitrc@${initrc}@g        ${callserial}
                  #------------------------------------------------------------------------#



                  #----- Re-submit the job to moorcroft2c queue. --------------------------#
                  ${srun}
                  #------------------------------------------------------------------------#



                  #----- Remove information but keep the crashing sign. -------------------#
                  year=${yeara}
                  month=${montha}
                  day=${datea}
                  hhmm=${timea}
                  agb="NA"
                  bsa="NA"
                  lai="NA"
                  scb="NA"
                  npa="NA"
                  #------------------------------------------------------------------------#

               else
                  #----- Give up on this node. --------------------------------------------#
                  echo "${ff}  :-{ Tolerance is tiny and it still crashes.  I gave up..."
                  #----- Delete files as they are the simulation can't run. ---------------#
                  /bin/rm -fvr ${here}/${polyname}/analy/*
                  /bin/rm -fvr ${here}/${polyname}/histo/*
                  #------------------------------------------------------------------------#
               fi # [ ${toler} != "0.0000010","0.0000001","00" ]
               #---------------------------------------------------------------------------#
            else
               #---------------------------------------------------------------------------#
               #      Something bad happened, but it has already been re-submitted.        #
               #---------------------------------------------------------------------------#
               if [ ${runt} == "CRASHED" ]
               then
                  echo "${ffout}: ${polyname} HAS CRASHED (RK4 PROBLEM)... <==========="
                  let ncrashed=${ncrashed}+1
               elif [ ${runt} == "SIGSEGV" ]
               then
                  let nsigsegv=${nsigsegv}+1
                  echo "${ffout}: ${polyname} HAD SEGMENTATION VIOLATION... <==========="
               elif [ ${runt} == "METMISS" ]
               then
                  let nmetmiss=${nmetmiss}+1
                  echo "${ffout}: ${polyname} DID NOT FIND MET DRIVERS... <==========="
               elif [ ${runt} == "BAD_MET" ]
               then
                  let nbad_met=${nbad_met}+1
                  echo "${ffout}: ${polyname} HAS BAD MET DRIVERS... <==========="
               elif [ ${runt} == "STOPPED" ]
               then
                  let nstopped=${nstopped}+1
                  echo "${ffout}: ${polyname} STOPPED (UNKNOWN REASON)... <==========="
               fi
               #---------------------------------------------------------------------------#
            fi # [ -s ${serial_out} ]
            #------------------------------------------------------------------------------#

         elif [ ${runt} == "THE_END" ]
         then 
            echo "${ffout}: ${polyname} has finished o/\o"

         elif [ ${runt} == "EXTINCT" -o ${runt} == "STSTATE" ] 
         then
            running=$(squeue -h -u ${moi} -n ${jobname} -o "%T")
            if [ "x${running}" == "xRUNNING" ]
            then
               deathrow="${deathrow} ${jobname}"
               if [ ${runt} == "EXTINCT" ]
               then
                  let newextinct=${newextinct}+1
                  echo "${ffout}: ${polyname} has become a desert..."
               elif [ ${runt} == "STSTATE" ]
               then
                  let newststate=${newststate}+1
                  echo "${ffout}: ${polyname} has reached steady state..."
               fi
            fi
         elif [ -s ${serial_out} ]
         then
            #----- Check whether the simulation is running, and when in model time it is. -#
            stdout="${here}/${polyname}/serial_out.out"
            simline=$(grep "Simulating: "   ${stdout} | tail -1)
            runtime=$(echo ${simline} | awk '{print $3}')
            simline=$(grep "Simulating: "   ${stdout} | tail -1)
            runtime=$(echo ${simline} | awk '{print $3}')
            echo "${ffout}: ${polyname} is running (${runtime})..."
            #------------------------------------------------------------------------------#
         elif [ ${runt} == "HISTORY" ] 
         then
            echo "${ffout}: ${polyname} is pending again ..."
         elif [ ${runt} == "INITIAL" ]
         then
            echo "${ffout}: ${polyname} is pending (never started)..."
         else
            echo "${ffout}: ${polyname} status is unclear..."
         fi # [ ${runt} == "CRASHED" ]...
         #---------------------------------------------------------------------------------#
         ;;
      esac
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #    Last, we check whether the job is still running or not.  In case it is not      #
      # but it should, update the ED2IN with the latest history and re-submit the job.     #
      #------------------------------------------------------------------------------------#
      if [ ${runt} != "THE_END" ] && [ ${runt} != "STSTATE" ] && 
         [ ${runt} != "EXTINCT" ] && [ ${runt} != "BAD_MET" ]
      then
         running=$(squeue -h -u ${moi} -n ${jobname} -o "%T")
         if [ "x${running}" == "x" ]
         then
            running="NOT_FOUND"
         fi

         if [ "x${running}" != "xRUNNING" ]
         then
            echo "${polyname} is missing.  Re-submitting..." >> ${situation}

            let nresubmit=${nresubmit}+1

            #------------------------------------------------------------------------------#
            #     Job was kicked out, submit it again.                                     #
            #------------------------------------------------------------------------------#
            ED2IN="${here}/${polyname}/ED2IN"
            #----- Save the tolerance before we delete the file. --------------------------#
            toler=$(grep NL%RK4_TOLERANCE ${ED2IN} | awk '{print $3}')
            #------------------------------------------------------------------------------#
            #   This should never happen, but...                                           #
            #------------------------------------------------------------------------------#
            if [ "x${oldtol}" == "xmytoler" -o "x${oldtol}" == "x" ]
            then
               toler=0.01
            else
               toler=${oldtol}
            fi # [ "x${oldtol}" == "xmytoler" -o "x${oldtol}" == "x" ]
            #------------------------------------------------------------------------------#


            #----- Re-generate ED2IN. -----------------------------------------------------#
            rm -f ${ED2IN} 
            rm -f ${here}/${polyname}/serial_out.out
            rm -f ${here}/${polyname}/serial_lsf.out
            rm -f ${here}/${polyname}/serial_out.err
            #------------------------------------------------------------------------------#



            #----- Copy the Template to the right directory. ------------------------------#
            cp ${here}/Template/ED2IN ${ED2IN}
            #------------------------------------------------------------------------------#




            #----- If last history is the first year, switch it to initial. ---------------#
            if [ ${year} -eq ${yeara} ] || [ ${runt} == "INITIAL" ]
            then
               runflag="INITIAL"
            else
               runflag="HISTORY"
            fi
            #------------------------------------------------------------------------------#




            #----- Check whether to use SFILIN as restart or history. ---------------------#
            if [ ${runflag} == "INITIAL" ] && [ ${initmode} -eq 6 ]
            then
               thissfilin=${thisbiomin}
            elif [ ${runflag} == "INITIAL" ] && [ ${initmode} -eq 5 ]
            then
               thissfilin=${restart}
            else
               thissfilin=${there}/${polyname}/histo/${polyname}
            fi
            #------------------------------------------------------------------------------#




            #----- Update the polygon information on ED2IN. -------------------------------#
            sed -i s@paththere@${there}@g                ${ED2IN}
            sed -i s@RUNFLAG@${runflag}@g                ${ED2IN}
            sed -i s@myyeara@${thisyeara}@g              ${ED2IN}
            sed -i s@mymontha@${montha}@g                ${ED2IN}
            sed -i s@mydatea@${datea}@g                  ${ED2IN}
            sed -i s@mytimea@${timea}@g                  ${ED2IN}
            sed -i s@myyearz@${thisyearz}@g              ${ED2IN}
            sed -i s@mymonthz@${monthz}@g                ${ED2IN}
            sed -i s@mydatez@${datez}@g                  ${ED2IN}
            sed -i s@mytimez@${timez}@g                  ${ED2IN}
            sed -i s@mydtlsm@${dtlsm}@g                  ${ED2IN}
            sed -i s@thispoly@${polyname}@g              ${ED2IN}
            sed -i s@plonflag@${polylon}@g               ${ED2IN}
            sed -i s@platflag@${polylat}@g               ${ED2IN}
            sed -i s@timehhhh@${timeh}@g                 ${ED2IN}
            sed -i s@datehhhh@${dateh}@g                 ${ED2IN}
            sed -i s@monthhhh@${monthh}@g                ${ED2IN}
            sed -i s@yearhhhh@${yearh}@g                 ${ED2IN}
            sed -i s@myinitmode@${initmode}@g            ${ED2IN}
            sed -i s@mysfilin@${thissfilin}@g            ${ED2IN}
            sed -i s@mytrees@${pfts}@g                   ${ED2IN}
            sed -i s@mycrop@${crop}@g                    ${ED2IN}
            sed -i s@myplantation@${plantation}@g        ${ED2IN}
            sed -i s@myiphen@${iphen}@g                  ${ED2IN}
            sed -i s@myallom@${iallom}@g                 ${ED2IN}
            sed -i s@myisoilflg@${polyisoil}@g           ${ED2IN}
            sed -i s@mynslcon@${polyntext}@g             ${ED2IN}
            sed -i s@myslxsand@${polysand}@g             ${ED2IN}
            sed -i s@myslxclay@${polyclay}@g             ${ED2IN}
            sed -i s@mysoilbc@${polysoilbc}@g            ${ED2IN}
            sed -i s@mysldrain@${polysldrain}@g          ${ED2IN}
            sed -i s@mysoilcol@${polycol}@g              ${ED2IN}
            sed -i s@mynzg@${polynzg}@g                  ${ED2IN}
            sed -i s@mymetdriverdb@${metdriverdb}@g      ${ED2IN}
            sed -i s@mymetcyc1@${metcyc1}@g              ${ED2IN}
            sed -i s@mymetcycf@${metcycf}@g              ${ED2IN}
            sed -i s@mytoler@${toler}@g                  ${ED2IN}
            sed -i s@myvmfactc3@${vmfactc3}@g            ${ED2IN}
            sed -i s@myvmfactc4@${vmfactc4}@g            ${ED2IN}
            sed -i s@mymphototrc3@${mphototrc3}@g        ${ED2IN}
            sed -i s@mymphototec3@${mphototec3}@g        ${ED2IN}
            sed -i s@mymphotoc4@${mphotoc4}@g            ${ED2IN}
            sed -i s@mybphotoblc3@${bphotoblc3}@g        ${ED2IN}
            sed -i s@mybphotonlc3@${bphotonlc3}@g        ${ED2IN}
            sed -i s@mybphotoc4@${bphotoc4}@g            ${ED2IN}
            sed -i s@mykwgrass@${kwgrass}@g              ${ED2IN}
            sed -i s@mykwtree@${kwtree}@g                ${ED2IN}
            sed -i s@mygammac3@${gammac3}@g              ${ED2IN}
            sed -i s@mygammac4@${gammac4}@g              ${ED2IN}
            sed -i s@myd0grass@${d0grass}@g              ${ED2IN}
            sed -i s@myd0tree@${d0tree}@g                ${ED2IN}
            sed -i s@myalphac3@${alphac3}@g              ${ED2IN}
            sed -i s@myalphac4@${alphac4}@g              ${ED2IN}
            sed -i s@myklowco2@${klowco2}@g              ${ED2IN}
            sed -i s@mydecomp@${decomp}@g                ${ED2IN}
            sed -i s@myrrffact@${rrffact}@g              ${ED2IN}
            sed -i s@mygrowthresp@${growthresp}@g        ${ED2IN}
            sed -i s@mylwidthgrass@${lwidthgrass}@g      ${ED2IN}
            sed -i s@mylwidthbltree@${lwidthbltree}@g    ${ED2IN}
            sed -i s@mylwidthnltree@${lwidthnltree}@g    ${ED2IN}
            sed -i s@myq10c3@${q10c3}@g                  ${ED2IN}
            sed -i s@myq10c4@${q10c4}@g                  ${ED2IN}
            sed -i s@myh2olimit@${h2olimit}@g            ${ED2IN}
            sed -i s@mymortscheme@${imortscheme}@g       ${ED2IN}
            sed -i s@myddmortconst@${ddmortconst}@g      ${ED2IN}
            sed -i s@mycbrscheme@${cbrscheme}@g          ${ED2IN}
            sed -i s@mysfclyrm@${isfclyrm}@g             ${ED2IN}
            sed -i s@myicanturb@${icanturb}@g            ${ED2IN}
            sed -i s@myatmco2@${atmco2}@g                ${ED2IN}
            sed -i s@mythcrit@${thcrit}@g                ${ED2IN}
            sed -i s@mysmfire@${smfire}@g                ${ED2IN}
            sed -i s@myfire@${ifire}@g                   ${ED2IN}
            sed -i s@myfuel@${fireparm}@g                ${ED2IN}
            sed -i s@mymetavg@${imetavg}@g               ${ED2IN}
            sed -i s@mypercol@${ipercol}@g               ${ED2IN}
            sed -i s@myrunoff@${runoff}@g                ${ED2IN}
            sed -i s@mymetrad@${imetrad}@g               ${ED2IN}
            sed -i s@mybranch@${ibranch}@g               ${ED2IN}
            sed -i s@mycanrad@${icanrad}@g               ${ED2IN}
            sed -i s@mycrown@${crown}@g                  ${ED2IN}
            sed -i s@myltransvis@${ltransvis}@g          ${ED2IN}
            sed -i s@myltransnir@${ltransnir}@g          ${ED2IN}
            sed -i s@mylreflectvis@${lreflectvis}@g      ${ED2IN}
            sed -i s@mylreflectnir@${lreflectnir}@g      ${ED2IN}
            sed -i s@myorienttree@${orienttree}@g        ${ED2IN}
            sed -i s@myorientgrass@${orientgrass}@g      ${ED2IN}
            sed -i s@myclumptree@${clumptree}@g          ${ED2IN}
            sed -i s@myclumpgrass@${clumpgrass}@g        ${ED2IN}
            sed -i s@myvegtdyn@${ivegtdyn}@g             ${ED2IN}
            sed -i s@mybigleaf@${ibigleaf}@g             ${ED2IN}
            sed -i s@myrepro@${irepro}@g                 ${ED2IN}
            sed -i s@myubmin@${ubmin}@g                  ${ED2IN}
            sed -i s@myugbmin@${ugbmin}@g                ${ED2IN}
            sed -i s@myustmin@${ustmin}@g                ${ED2IN}
            sed -i s@mygamm@${gamm}@g                    ${ED2IN}
            sed -i s@mygamh@${gamh}@g                    ${ED2IN}
            sed -i s@mytprandtl@${tprandtl}@g            ${ED2IN}
            sed -i s@myribmax@${ribmax}@g                ${ED2IN}
            sed -i s@mygndvap@${igndvap}@g               ${ED2IN}
            sed -i s@mydtcensus@${dtcensus}@g            ${ED2IN}
            sed -i s@myyr1stcensus@${yr1stcensus}@g      ${ED2IN}
            sed -i s@mymon1stcensus@${mon1stcensus}@g    ${ED2IN}
            sed -i s@myminrecruitdbh@${minrecruitdbh}@g  ${ED2IN}
            sed -i s@mytreefall@${treefall}@g            ${ED2IN}
            sed -i s@mymaxpatch@${iage}@g                ${ED2IN}
            sed -i s@myanthdisturb@${ianthdisturb}@g     ${ED2IN}
            sed -i s@myludatabase@${ludatabase}@g        ${ED2IN}
            #------ Soil variables. -------------------------------------------------------#
            sed -i s@myslz1@"${polyslz1}"@g           ${ED2IN}
            sed -i s@myslz2@"${polyslz2}"@g           ${ED2IN}
            sed -i s@myslz3@"${polyslz3}"@g           ${ED2IN}
            sed -i s@myslz4@"${polyslz4}"@g           ${ED2IN}
            sed -i s@myslz5@"${polyslz5}"@g           ${ED2IN}
            sed -i s@myslz6@"${polyslz6}"@g           ${ED2IN}
            sed -i s@myslz7@"${polyslz7}"@g           ${ED2IN}
            sed -i s@myslmstr1@"${polyslm1}"@g        ${ED2IN}
            sed -i s@myslmstr2@"${polyslm2}"@g        ${ED2IN}
            sed -i s@myslmstr3@"${polyslm3}"@g        ${ED2IN}
            sed -i s@myslmstr4@"${polyslm4}"@g        ${ED2IN}
            sed -i s@myslmstr5@"${polyslm5}"@g        ${ED2IN}
            sed -i s@myslmstr6@"${polyslm6}"@g        ${ED2IN}
            sed -i s@myslmstr7@"${polyslm7}"@g        ${ED2IN}
            sed -i s@mystgoff1@"${polyslt1}"@g        ${ED2IN}
            sed -i s@mystgoff2@"${polyslt2}"@g        ${ED2IN}
            sed -i s@mystgoff3@"${polyslt3}"@g        ${ED2IN}
            sed -i s@mystgoff4@"${polyslt4}"@g        ${ED2IN}
            sed -i s@mystgoff5@"${polyslt5}"@g        ${ED2IN}
            sed -i s@mystgoff6@"${polyslt6}"@g        ${ED2IN}
            sed -i s@mystgoff7@"${polyslt7}"@g        ${ED2IN}
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #      A halted polygon will inevitably go to the moorcroft2c queue.           #
            # Re-write the srun.sh file.                                                   #
            #------------------------------------------------------------------------------#
            srun="${here}/${polyname}/srun.sh"
            callserial="${here}/${polyname}/callserial.sh"

            #----- Generate random waiting time. ------------------------------------------#
            let wtime=${ff}%8
            let wtime=${wtime}*20
            nudge=$(date +%S)
            if [ ${nudge} -lt 10 ]
            then 
               nudge=$(echo ${nudge} | awk '{print substr($1,2,1)}')
            fi
            let nudge=${nudge}%15
            let wtime=${wtime}+${nudge}
            let wtime=${wtime}+2
            #------------------------------------------------------------------------------#



            #----- Refresh srun and callserial. -------------------------------------------#
            /bin/rm -f ${srun}
            /bin/rm -f ${callserial}
            cp ${here}/Template/srun.sh       ${srun}
            cp ${here}/Template/callserial.sh ${callserial}
            #------------------------------------------------------------------------------#



            #----- The new queue is by default the general. -------------------------------#
            newqueue="general"
            runtime=${runtime_gen}
            #------------------------------------------------------------------------------#



            #----- Change some settings in srun.sh. ---------------------------------------#
            sed -i s@pathhere@${here}@g      ${srun}
            sed -i s@paththere@${here}@g     ${srun}
            sed -i s@thispoly@${polyname}@g  ${srun}
            sed -i s@thisqueue@${newqueue}@g ${srun}
            sed -i s@thisdesc@${desc}@g      ${srun}
            sed -i s@zzzzzzzz@${wtime}@g     ${srun}
            sed -i s@myorder@${ff}@g         ${srun}
            sed -i s@myinitrc@${initrc}@g    ${srun}
            sed -i s@thismemory@${memory}@g  ${srun}
            sed -i s@thistime@${runtime}@g   ${srun}
            #------------------------------------------------------------------------------#




            #----- Change the callserial.sh file. -----------------------------------------#
            sed -i s@thisroot@${here}@g          ${callserial}
            sed -i s@thispoly@${polyname}@g      ${callserial}
            sed -i s@myexec@${execname}@g        ${callserial}
            sed -i s@myname@${moi}@g             ${callserial}
            sed -i s@mypackdata@${packdatasrc}@g ${callserial}
            sed -i s@myscenario@${iscenario}@g   ${callserial}
            sed -i s@myscenmain@${scentype}@g    ${callserial}
            sed -i s@myinitrc@${initrc}@g        ${callserial}
            #------------------------------------------------------------------------------#


            #----- Re-submit the job to moorcroft2c queue. --------------------------------#
            ${srun}

            #----- Switch the status to pending. ------------------------------------------#
            year=${yeara}
            month=${montha}
            day=${datea}
            hhmm=${timea}
            runt="INITIAL"
            agb="NA"
            bsa="NA"
            lai="NA"
            scb="NA"
            npa="NA"
         else
            echo "${polyname} is running/pending..." >> ${situation}
         fi # [ "x${running}" != "xRUNNING" ]
      else
         echo "${polyname} status is ${runt}.  No need to resubmit." >> ${situation}
      fi # [ ${runt} == "HISTORY" -a ${queue} == "moorcroft2c" ]
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Test whether the job can be submitted.                                         #
      #------------------------------------------------------------------------------------#
      if [ ${runt} == "INITIAL" ] && [ -s ${here}/${polyname} ]
      then
         running=$(squeue -h -u ${moi} -n ${jobname} -o "%T")
         if [ "x${running}" == "x" ]
         then
            running="NOT_FOUND"
         fi
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Check whether to start a new simulation or wait.                            #
         #---------------------------------------------------------------------------------#
         if [ "x${waitpot}" == "xy" ] || [ "x${waitpot}" == "xy" ]
         then
            readytostart=$(/bin/ls -1 ${potveg}/${oldname}*h5 2> /dev/null | wc -l)
         else
            readytostart=1
         fi
         if [ "x${running}" != "xRUNNING" ] && [ ${readytostart} -eq 1 ]
         then
            let nstart=${nstart}+1
            echo " >>> Initial submission of polygon ${polyname}!!!!"

            #----- Re-submit the job to moorcroft2c queue. --------------------------------#
            "${here}/${polyname}/srun.sh"

         elif [ ${readytostart} -gt 1 ]
         then
            echo "  Ouch!  Something strange with the restart for polygon ${polyname}..."
         fi
         #---------------------------------------------------------------------------------#
      fi
      #---------------------------------------------------------------------------------------#



      #----- Write polygon check into a single table. -------------------------------------#
      output="${polyname} ${polylon} ${polylat} ${year} ${month} ${day} ${hhmm}"
      output="${output} ${runt} ${agb} ${bsa} ${lai} ${scb} ${npa}"
      echo ${output} >> ${outcheck}
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #    Update joborder.txt file.                                                       #
      #------------------------------------------------------------------------------------#
      if [ ${runt} != "THE_END" ] && [ ${runt} != "STSTATE" ] && [ ${runt} != "EXTINCT" ]
      then
         polyIATA=$(echo ${polyiata} | tr '[:lower:]' '[:upper:]')
         polyqueue=$(squeue -h -u ${moi} -n ${jobname} -o "%P")

         if [ "x${polyqueue}" != "x" ]
         then
            #----- "Normal" job, with regular name. ---------------------------------------#
            newqueue=${polyqueue}
            success="TRUE"
         else
            #------------------------------------------------------------------------------#
            #    Job was interrupted.  We resubmit the job in moorcroft2c.                 #
            #------------------------------------------------------------------------------#
            success="FALSE"
         fi
         if [ ${success} == "TRUE" ] 
         then
            oldline=${oi}
            newline=$(echo ${oldline} | sed s@${queue}@${newqueue}@g)
            sed -i s@"${oldline}"@"${newline}"@g ${joborder}
         fi
      elif [ ! -s ${here}/${polyname} ]
      then
         newqueue="general"
         oldline=${oi}
         newline=$(echo ${oldline} | sed s@${queue}@${newqueue}@g)
         sed -i s@"${oldline}"@"${newline}"@g ${joborder}
      fi
      
      #------------------------------------------------------------------------------------#
      #    Check whether the steady-state restart file is copied to the destination        #
      # directory.                                                                         #
      #------------------------------------------------------------------------------------#
      if [ "x${copyrestart}" == "xy" ] || [ "x${copyrestart}" == "xY" ] &&
         [ -s ${here}/${polyname} ]
      then
         if [ ${runt} == "STSTATE" -o ${runt} == "THE_END" -o ${runt} == "EXTINCT" ]
         then
            nrestart=$(/bin/ls -1 ${restart}/${polyname}-S-*h5 2> /dev/null | wc -l)
            #----- Check whether there is any file in the restart directory ---------------#
            if [ ${nrestart} -eq 0 ]
            then
               lastht=$(/bin/ls -1 ${here}/${polyname}/histo/*-S-*h5 2> /dev/null | tail -1)
               baseht=$(basename ${lastht})
               destht="${restart}/${baseht}"
               echo " - Copying file ${baseht} to the restart directory..."
               echo " - LASTHISTO = ${lastht}"
               echo " - DESTHISTO = ${destht}"
               /bin/cp -uv ${lastht} ${destht}
            fi
            #------------------------------------------------------------------------------#
         fi
         #---------------------------------------------------------------------------------#
      fi
      #------------------------------------------------------------------------------------#
   done
   #---------------------------------------------------------------------------------------#
fi # end if
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#     Check whether to fix the queues in joborder.                                         #
#------------------------------------------------------------------------------------------#
if [ ${fixqueue} -eq 1 ]
then
   while [ ${ff} -lt ${npolys} ]
   do
      let ff=${ff}+1
      let line=${ff}+3
      let rem=${ff}%${printevery}
      
      if [ ${rem} -eq 0 -o ${ff} -eq ${npolys} ]
      then
         echo " - Fixed queues for ${ff} polygons so far..."
      fi

      #------------------------------------------------------------------------------------#
      #      Read the lineth line of the polygon list.  There must be smarter ways of do-  #
      # ing this, but this works.  Here we obtain the polygon name, and its longitude and  #
      # latitude.                                                                          #
      #------------------------------------------------------------------------------------#
      oi=$(head -${line} ${lonlat} | tail -1)
      polyname=$(echo ${oi}     | awk '{print $1 }')
      polyiata=$(echo ${oi}     | awk '{print $2 }')
      polylon=$(echo ${oi}      | awk '{print $3 }')
      polylat=$(echo ${oi}      | awk '{print $4 }')
      yeara=$(echo ${oi}        | awk '{print $5 }')
      montha=$(echo ${oi}       | awk '{print $6 }')
      datea=$(echo ${oi}        | awk '{print $7 }')
      timea=$(echo ${oi}        | awk '{print $8 }')
      yearz=$(echo ${oi}        | awk '{print $9 }')
      monthz=$(echo ${oi}       | awk '{print $10}')
      datez=$(echo ${oi}        | awk '{print $11}')
      timez=$(echo ${oi}        | awk '{print $12}')
      initmode=$(echo ${oi}     | awk '{print $13}')
      iscenario=$(echo ${oi}    | awk '{print $14}')
      isizepft=$(echo ${oi}     | awk '{print $15}')
      iage=$(echo ${oi}         | awk '{print $16}')
      polyisoil=$(echo ${oi}    | awk '{print $17}')
      polyntext=$(echo ${oi}    | awk '{print $18}')
      polysand=$(echo ${oi}     | awk '{print $19}')
      polyclay=$(echo ${oi}     | awk '{print $20}')
      polydepth=$(echo ${oi}    | awk '{print $21}')
      polysoilbc=$(echo ${oi}   | awk '{print $22}')
      polysldrain=$(echo ${oi}  | awk '{print $23}')
      polycol=$(echo ${oi}      | awk '{print $24}')
      slzres=$(echo ${oi}       | awk '{print $25}')
      queue=$(echo ${oi}        | awk '{print $26}')
      metdriver=$(echo ${oi}    | awk '{print $27}')
      dtlsm=$(echo ${oi}        | awk '{print $28}')
      vmfactc3=$(echo ${oi}     | awk '{print $29}')
      vmfactc4=$(echo ${oi}     | awk '{print $30}')
      mphototrc3=$(echo ${oi}   | awk '{print $31}')
      mphototec3=$(echo ${oi}   | awk '{print $32}')
      mphotoc4=$(echo ${oi}     | awk '{print $33}')
      bphotoblc3=$(echo ${oi}   | awk '{print $34}')
      bphotonlc3=$(echo ${oi}   | awk '{print $35}')
      bphotoc4=$(echo ${oi}     | awk '{print $36}')
      kwgrass=$(echo ${oi}      | awk '{print $37}')
      kwtree=$(echo ${oi}       | awk '{print $38}')
      gammac3=$(echo ${oi}      | awk '{print $39}')
      gammac4=$(echo ${oi}      | awk '{print $40}')
      d0grass=$(echo ${oi}      | awk '{print $41}')
      d0tree=$(echo ${oi}       | awk '{print $42}')
      alphac3=$(echo ${oi}      | awk '{print $43}')
      alphac4=$(echo ${oi}      | awk '{print $44}')
      klowco2=$(echo ${oi}      | awk '{print $45}')
      decomp=$(echo ${oi}       | awk '{print $46}')
      rrffact=$(echo ${oi}      | awk '{print $47}')
      growthresp=$(echo ${oi}   | awk '{print $48}')
      lwidthgrass=$(echo ${oi}  | awk '{print $49}')
      lwidthbltree=$(echo ${oi} | awk '{print $50}')
      lwidthnltree=$(echo ${oi} | awk '{print $51}')
      q10c3=$(echo ${oi}        | awk '{print $52}')
      q10c4=$(echo ${oi}        | awk '{print $53}')
      h2olimit=$(echo ${oi}     | awk '{print $54}')
      imortscheme=$(echo ${oi}  | awk '{print $55}')
      ddmortconst=$(echo ${oi}  | awk '{print $56}')
      cbrscheme=$(echo ${oi}    | awk '{print $57}')
      isfclyrm=$(echo ${oi}     | awk '{print $58}')
      icanturb=$(echo ${oi}     | awk '{print $59}')
      ubmin=$(echo ${oi}        | awk '{print $60}')
      ugbmin=$(echo ${oi}       | awk '{print $61}')
      ustmin=$(echo ${oi}       | awk '{print $62}')
      gamm=$(echo ${oi}         | awk '{print $63}')
      gamh=$(echo ${oi}         | awk '{print $64}')
      tprandtl=$(echo ${oi}     | awk '{print $65}')
      ribmax=$(echo ${oi}       | awk '{print $66}')
      atmco2=$(echo ${oi}       | awk '{print $67}')
      thcrit=$(echo ${oi}       | awk '{print $68}')
      smfire=$(echo ${oi}       | awk '{print $69}')
      ifire=$(echo ${oi}        | awk '{print $70}')
      fireparm=$(echo ${oi}     | awk '{print $71}')
      ipercol=$(echo ${oi}      | awk '{print $72}')
      runoff=$(echo ${oi}       | awk '{print $73}')
      imetrad=$(echo ${oi}      | awk '{print $74}')
      ibranch=$(echo ${oi}      | awk '{print $75}')
      icanrad=$(echo ${oi}      | awk '{print $76}')
      crown=$(echo   ${oi}      | awk '{print $77}')
      ltransvis=$(echo ${oi}    | awk '{print $78}')
      lreflectvis=$(echo ${oi}  | awk '{print $79}')
      ltransnir=$(echo ${oi}    | awk '{print $80}')
      lreflectnir=$(echo ${oi}  | awk '{print $81}')
      orienttree=$(echo ${oi}   | awk '{print $82}')
      orientgrass=$(echo ${oi}  | awk '{print $83}')
      clumptree=$(echo ${oi}    | awk '{print $84}')
      clumpgrass=$(echo ${oi}   | awk '{print $85}')
      ivegtdyn=$(echo ${oi}     | awk '{print $86}')
      igndvap=$(echo ${oi}      | awk '{print $87}')
      iphen=$(echo ${oi}        | awk '{print $88}')
      iallom=$(echo ${oi}       | awk '{print $89}')
      ibigleaf=$(echo ${oi}     | awk '{print $90}')
      irepro=$(echo ${oi}       | awk '{print $91}')
      treefall=$(echo ${oi}     | awk '{print $92}')
      ianthdisturb=$(echo ${oi} | awk '{print $93}')
      ianthdataset=$(echo ${oi} | awk '{print $94}')
      #------------------------------------------------------------------------------------#


      jobname="${desc}-${polyname}"
      jobqueue=$(squeue -h -u ${moi} -n ${jobname} -o "%P")

      #----- Update the queue in joborder.txt. --------------------------------------------#
      if [ "x${jobqueue}" != "x" ]
      then
         newqueue=${jobqueue}
      else
         newqueue=${queue}
      fi
      oldline=${oi}
      newline=$(echo ${oldline} | sed s@${queue}@${newqueue}@g)
      sed -i s@"${oldline}"@"${newline}"@g ${joborder}
      #------------------------------------------------------------------------------------#
   done
   #---------------------------------------------------------------------------------------#
   /bin/rm -f ${here}/run_sitter.lock
   exit
   #---------------------------------------------------------------------------------------#
fi
#------------------------------------------------------------------------------------------#



#----- Run R to make the status check. ----------------------------------------------------#
if [ ${plotstatus} -eq 1 ]
then
   echo "Running the status check..."
   Rscript_out="$(dirname ${Rscript_plot})/$(basename ${Rscript_plot} .r).txt"
   R CMD BATCH --no-save --no-restore ${Rscript_plot} ${Rscript_out}
fi
#------------------------------------------------------------------------------------------#



#----- Kill all jobs that can be killed. --------------------------------------------------#
echo "Killing polygons that have been suspended, went extinct, or reached steady state..."
if [ "x${deathrow}" != "x" ]
then
   for killme in ${deathrow}
   do
      scancel -u ${moi} -n ${killme}
   done # killme in ${deathrow}
fi # [ "x${deathrow}" != "x" ]
#------------------------------------------------------------------------------------------#



#----- Give a three minute break. ---------------------------------------------------------#
echo "Taking a nap before counting jobs..."
sleep 180

#----- Check how many jobs are left on our preferred queues. ------------------------------#
echo "Counting the number of polygons on queues"
#----- Moorcroft_6100. --------------------------------------------------------------------#
m61run=$(sinfo -h -p moorcroft_6100 -o "%C" | sed s@"/"@" "@g | awk '{print $1}')
m61idl=$(sinfo -h -p moorcroft_6100 -o "%C" | sed s@"/"@" "@g | awk '{print $2}')
m61oth=$(sinfo -h -p moorcroft_6100 -o "%C" | sed s@"/"@" "@g | awk '{print $3}')
m61tot=$(sinfo -h -p moorcroft_6100 -o "%C" | sed s@"/"@" "@g | awk '{print $4}')
m61moi=$(scount   -p moorcroft_6100 -u ${moi} )
#----- Moorcroft_amd. ---------------------------------------------------------------------#
amdrun=$(sinfo -h -p moorcroft_amd -o "%C" | sed s@"/"@" "@g | awk '{print $1}')
amdidl=$(sinfo -h -p moorcroft_amd -o "%C" | sed s@"/"@" "@g | awk '{print $2}')
amdoth=$(sinfo -h -p moorcroft_amd -o "%C" | sed s@"/"@" "@g | awk '{print $3}')
amdtot=$(sinfo -h -p moorcroft_amd -o "%C" | sed s@"/"@" "@g | awk '{print $4}')
amdmoi=$(scount   -p moorcroft_amd -u ${moi} )
#----- Moorcroft_amd. ---------------------------------------------------------------------#
urerun=$(sinfo -h -p unrestricted -o "%C" | sed s@"/"@" "@g | awk '{print $1}')
ureidl=$(sinfo -h -p unrestricted -o "%C" | sed s@"/"@" "@g | awk '{print $2}')
ureoth=$(sinfo -h -p unrestricted -o "%C" | sed s@"/"@" "@g | awk '{print $3}')
uretot=$(sinfo -h -p unrestricted -o "%C" | sed s@"/"@" "@g | awk '{print $4}')
uremoi=$(scount   -p unrestricted -u ${moi} )
#----- Moorcroft_amd. ---------------------------------------------------------------------#
urerun=$(sinfo -h -p unrestricted -o "%C" | sed s@"/"@" "@g | awk '{print $1}')
ureidl=$(sinfo -h -p unrestricted -o "%C" | sed s@"/"@" "@g | awk '{print $2}')
ureoth=$(sinfo -h -p unrestricted -o "%C" | sed s@"/"@" "@g | awk '{print $3}')
uretot=$(sinfo -h -p unrestricted -o "%C" | sed s@"/"@" "@g | awk '{print $4}')
uremoi=$(scount   -p unrestricted -u ${moi} )
#----- Create a file with the pending jobs. -----------------------------------------------#
echo "Counting the number of polygons in queue on general..."
genpen=$(scount -g "${desc}-" -p general -t PENDING)
genrun=$(scount -g "${desc}-" -p general -t RUNNING)
/bin/rm -f ${pendfile}
squeue -h -t PENDING -u ${moi} -p general -o "%.j %.N %.P" | grep ${desc} > ${pendfile}
#------------------------------------------------------------------------------------------#



#----- Find the available space on available queues. --------------------------------------#
echo "Looking for empty cores on the queues..."
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#     Check the amount of room left for each queue.  Make sure that we don't exceed the    #
# maximum availability or the maximum use allowed.                                         #
#------------------------------------------------------------------------------------------#
/bin/rm -f ${tablefile}
touch ${tablefile}
#------ moorcroft_6100. -------------------------------------------------------------------#
let m61pot=${m61max}-${m61moi}
if [ ${m61pot} -lt 0 ] || [ ${m61idl} -lt 0 ]
then
   m61avl=0
elif [ ${m61pot} -lt ${m61idl} ]
then
   m61avl=${m61pot}
else
   m61avl=${m61idl}
fi
echo " -----------------------------------------------------------" >> ${tablefile}
echo "   moorcroft_6100"                                            >> ${tablefile}
echo " "                                                            >> ${tablefile}
echo "   TOTAL     = ${m61tot}"                                     >> ${tablefile}
echo "   RUNNING   = ${m61run}"                                     >> ${tablefile}
echo "   PENDING   = ${m61pen}"                                     >> ${tablefile}
echo "   User max  = ${m61max}"                                     >> ${tablefile}
echo "   MY RUNS   = ${m61moi}"                                     >> ${tablefile}
echo "   POTENTIAL = ${m61pot}"                                     >> ${tablefile}
echo "   AVAILABLE = ${m61avl}"                                     >> ${tablefile}
echo " -----------------------------------------------------------" >> ${tablefile}
echo " "                                                            >> ${tablefile}
#------ moorcroft_amd. --------------------------------------------------------------------#
let amdpot=${amdmax}-${amdmoi}
if [ ${amdpot} -lt 0 ] || [ ${amdidl} -lt 0 ]
then
   amdavl=0
elif [ ${amdpot} -lt ${amdidl} ]
then
   amdavl=${amdpot}
else
   amdavl=${amdidl}
fi
echo " -----------------------------------------------------------" >> ${tablefile}
echo "   moorcroft_amd"                                             >> ${tablefile}
echo " "                                                            >> ${tablefile}
echo "   TOTAL     = ${amdtot}"                                     >> ${tablefile}
echo "   RUNNING   = ${amdrun}"                                     >> ${tablefile}
echo "   PENDING   = ${amdpen}"                                     >> ${tablefile}
echo "   User max  = ${amdmax}"                                     >> ${tablefile}
echo "   MY RUNS   = ${amdmoi}"                                     >> ${tablefile}
echo "   POTENTIAL = ${amdpot}"                                     >> ${tablefile}
echo "   AVAILABLE = ${amdavl}"                                     >> ${tablefile}
echo " -----------------------------------------------------------" >> ${tablefile}
echo " "                                                            >> ${tablefile}
#------ unrestricted. ---------------------------------------------------------------------#
let urepot=${uremax}-${uremoi}
if [ ${urepot} -lt 0 ] || [ ${ureidl} -lt 0 ]
then
   ureavl=0
elif [ ${urepot} -lt ${ureidl} ]
then
   ureavl=${urepot}
else
   ureavl=${ureidl}
fi
echo " -----------------------------------------------------------" >> ${tablefile}
echo "   unrestricted"                                              >> ${tablefile}
echo " "                                                            >> ${tablefile}
echo "   TOTAL     = ${uretot}"                                     >> ${tablefile}
echo "   RUNNING   = ${urerun}"                                     >> ${tablefile}
echo "   PENDING   = ${urepen}"                                     >> ${tablefile}
echo "   User max  = ${uremax}"                                     >> ${tablefile}
echo "   MY RUNS   = ${uremoi}"                                     >> ${tablefile}
echo "   POTENTIAL = ${urepot}"                                     >> ${tablefile}
echo "   AVAILABLE = ${ureavl}"                                     >> ${tablefile}
echo " -----------------------------------------------------------" >> ${tablefile}
echo " "                                                            >> ${tablefile}
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#      Count the total free room.                                                          #
#------------------------------------------------------------------------------------------#
let totavl=${m61avl}+${amdavl}+${ureavl}
let m61pamd=${m61avl}+${amdavl}
#------------------------------------------------------------------------------------------#



#----- Find out how many cores can be moved. ----------------------------------------------#
if [ ${genpen} -le ${totavl} ]
then
   nfill=${genpen}
else
   nfill=${totavl}
fi
#------------------------------------------------------------------------------------------#


#----- This loop is where the job's queue switching happens. ------------------------------#
echo "Switching queues for ${nfill} jobs..."
ff=0
ll=0
while [ ${ff} -lt ${nfill} ]
do
   let ll=${ll}+1
   let line=${ll}+3

   #---------------------------------------------------------------------------------------#
   #      Read the lineth line of the polygon list.  There must be smarter ways of doing   #
   # this, but this works.  Here we obtain the polygon name, and its longitude and         #
   # latitude.                                                                             #
   #---------------------------------------------------------------------------------------#
   oi=$(head -${line} ${lonlat} | tail -1)
   polyname=$(echo ${oi}     | awk '{print $1 }')
   polyiata=$(echo ${oi}     | awk '{print $2 }')
   polylon=$(echo ${oi}      | awk '{print $3 }')
   polylat=$(echo ${oi}      | awk '{print $4 }')
   yeara=$(echo ${oi}        | awk '{print $5 }')
   montha=$(echo ${oi}       | awk '{print $6 }')
   datea=$(echo ${oi}        | awk '{print $7 }')
   timea=$(echo ${oi}        | awk '{print $8 }')
   yearz=$(echo ${oi}        | awk '{print $9 }')
   monthz=$(echo ${oi}       | awk '{print $10}')
   datez=$(echo ${oi}        | awk '{print $11}')
   timez=$(echo ${oi}        | awk '{print $12}')
   initmode=$(echo ${oi}     | awk '{print $13}')
   iscenario=$(echo ${oi}    | awk '{print $14}')
   isizepft=$(echo ${oi}     | awk '{print $15}')
   iage=$(echo ${oi}         | awk '{print $16}')
   polyisoil=$(echo ${oi}    | awk '{print $17}')
   polyntext=$(echo ${oi}    | awk '{print $18}')
   polysand=$(echo ${oi}     | awk '{print $19}')
   polyclay=$(echo ${oi}     | awk '{print $20}')
   polydepth=$(echo ${oi}    | awk '{print $21}')
   polysoilbc=$(echo ${oi}   | awk '{print $22}')
   polysldrain=$(echo ${oi}  | awk '{print $23}')
   polycol=$(echo ${oi}      | awk '{print $24}')
   slzres=$(echo ${oi}       | awk '{print $25}')
   queue=$(echo ${oi}        | awk '{print $26}')
   metdriver=$(echo ${oi}    | awk '{print $27}')
   dtlsm=$(echo ${oi}        | awk '{print $28}')
   vmfactc3=$(echo ${oi}     | awk '{print $29}')
   vmfactc4=$(echo ${oi}     | awk '{print $30}')
   mphototrc3=$(echo ${oi}   | awk '{print $31}')
   mphototec3=$(echo ${oi}   | awk '{print $32}')
   mphotoc4=$(echo ${oi}     | awk '{print $33}')
   bphotoblc3=$(echo ${oi}   | awk '{print $34}')
   bphotonlc3=$(echo ${oi}   | awk '{print $35}')
   bphotoc4=$(echo ${oi}     | awk '{print $36}')
   kwgrass=$(echo ${oi}      | awk '{print $37}')
   kwtree=$(echo ${oi}       | awk '{print $38}')
   gammac3=$(echo ${oi}      | awk '{print $39}')
   gammac4=$(echo ${oi}      | awk '{print $40}')
   d0grass=$(echo ${oi}      | awk '{print $41}')
   d0tree=$(echo ${oi}       | awk '{print $42}')
   alphac3=$(echo ${oi}      | awk '{print $43}')
   alphac4=$(echo ${oi}      | awk '{print $44}')
   klowco2=$(echo ${oi}      | awk '{print $45}')
   decomp=$(echo ${oi}       | awk '{print $46}')
   rrffact=$(echo ${oi}      | awk '{print $47}')
   growthresp=$(echo ${oi}   | awk '{print $48}')
   lwidthgrass=$(echo ${oi}  | awk '{print $49}')
   lwidthbltree=$(echo ${oi} | awk '{print $50}')
   lwidthnltree=$(echo ${oi} | awk '{print $51}')
   q10c3=$(echo ${oi}        | awk '{print $52}')
   q10c4=$(echo ${oi}        | awk '{print $53}')
   h2olimit=$(echo ${oi}     | awk '{print $54}')
   imortscheme=$(echo ${oi}  | awk '{print $55}')
   ddmortconst=$(echo ${oi}  | awk '{print $56}')
   cbrscheme=$(echo ${oi}    | awk '{print $57}')
   isfclyrm=$(echo ${oi}     | awk '{print $58}')
   icanturb=$(echo ${oi}     | awk '{print $59}')
   ubmin=$(echo ${oi}        | awk '{print $60}')
   ugbmin=$(echo ${oi}       | awk '{print $61}')
   ustmin=$(echo ${oi}       | awk '{print $62}')
   gamm=$(echo ${oi}         | awk '{print $63}')
   gamh=$(echo ${oi}         | awk '{print $64}')
   tprandtl=$(echo ${oi}     | awk '{print $65}')
   ribmax=$(echo ${oi}       | awk '{print $66}')
   atmco2=$(echo ${oi}       | awk '{print $67}')
   thcrit=$(echo ${oi}       | awk '{print $68}')
   smfire=$(echo ${oi}       | awk '{print $69}')
   ifire=$(echo ${oi}        | awk '{print $70}')
   fireparm=$(echo ${oi}     | awk '{print $71}')
   ipercol=$(echo ${oi}      | awk '{print $72}')
   runoff=$(echo ${oi}       | awk '{print $73}')
   imetrad=$(echo ${oi}      | awk '{print $74}')
   ibranch=$(echo ${oi}      | awk '{print $75}')
   icanrad=$(echo ${oi}      | awk '{print $76}')
   crown=$(echo   ${oi}      | awk '{print $77}')
   ltransvis=$(echo ${oi}    | awk '{print $78}')
   lreflectvis=$(echo ${oi}  | awk '{print $79}')
   ltransnir=$(echo ${oi}    | awk '{print $80}')
   lreflectnir=$(echo ${oi}  | awk '{print $81}')
   orienttree=$(echo ${oi}   | awk '{print $82}')
   orientgrass=$(echo ${oi}  | awk '{print $83}')
   clumptree=$(echo ${oi}    | awk '{print $84}')
   clumpgrass=$(echo ${oi}   | awk '{print $85}')
   ivegtdyn=$(echo ${oi}     | awk '{print $86}')
   igndvap=$(echo ${oi}      | awk '{print $87}')
   iphen=$(echo ${oi}        | awk '{print $88}')
   iallom=$(echo ${oi}       | awk '{print $89}')
   ibigleaf=$(echo ${oi}     | awk '{print $90}')
   irepro=$(echo ${oi}       | awk '{print $91}')
   treefall=$(echo ${oi}     | awk '{print $92}')
   ianthdisturb=$(echo ${oi} | awk '{print $93}')
   ianthdataset=$(echo ${oi} | awk '{print $94}')
   #---------------------------------------------------------------------------------------#



   #----- Find the job name (or the alternative job name). --------------------------------#
   jobname="${desc}-${polyname}"
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #       If the job is pending, switch the queue.                                        #
   #---------------------------------------------------------------------------------------#
   jobpending=$(grep ${jobname} ${pendfile} | wc -l)
   if [ ${jobpending} -gt 0 ]
   then
      let ff=${ff}+1

      #----- Decide the new queue. --------------------------------------------------------#
      if [ ${ff} -le ${m61avl} ]
      then
         newqueue="moorcroft_6100"
      elif [ ${ff} -le ${m61pamd} ]
      then
         newqueue="moorcroft_amd"
      else
         newqueue="unrestricted"
      fi # [ ${ff} -le ${m2croom} ]
      #------------------------------------------------------------------------------------#



      #----- Moving to a queue without time restrictions, use longer runtime. -------------#
      runtime=${runtime_own}
      #------------------------------------------------------------------------------------#


      #----- Kill the job. ----------------------------------------------------------------#
      scancel -u ${moi} -n ${jobname}
      #------------------------------------------------------------------------------------#


      #----- Update the queue in joborder.txt. --------------------------------------------#
      oldline=${oi}
      newline=$(echo ${oldline} | sed s@${queue}@${newqueue}@g)
      sed -i s@"${oldline}"@"${newline}"@g ${joborder}
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #      Re-write srun.sh and callserial.sh to update queue and runtime information.   #
      #------------------------------------------------------------------------------------#
      srun="${here}/${polyname}/srun.sh"
      callserial="${here}/${polyname}/callserial.sh"
      #------------------------------------------------------------------------------------#


      #----- Set random waiting time. -----------------------------------------------------#
      let wtime=${ff}%8
      let wtime=${wtime}*20
      nudge=$(date +%S)
      if [ ${nudge} -lt 10 ]
      then 
         nudge=$(echo ${nudge} | awk '{print substr($1,2,1)}')
      fi
      let nudge=${nudge}%15
      let wtime=${wtime}+${nudge}
      let wtime=${wtime}+2
      #------------------------------------------------------------------------------------#




      #----- Reset srun.sh and callserial.sh. ---------------------------------------------#
      /bin/rm -f ${srun}
      /bin/rm -f ${callserial}
      cp ${here}/Template/srun.sh       ${srun}
      cp ${here}/Template/callserial.sh ${callserial}
      #------------------------------------------------------------------------------------#



      #----- Change some settings in srun.sh. ---------------------------------------------#
      sed -i s@pathhere@${here}@g      ${srun}
      sed -i s@paththere@${here}@g     ${srun}
      sed -i s@thispoly@${polyname}@g  ${srun}
      sed -i s@thisqueue@${newqueue}@g ${srun}
      sed -i s@thisdesc@${desc}@g      ${srun}
      sed -i s@zzzzzzzz@${wtime}@g     ${srun}
      sed -i s@myorder@${ff}@g         ${srun}
      sed -i s@myinitrc@${initrc}@g    ${srun}
      sed -i s@thismemory@${memory}@g  ${srun}
      sed -i s@thistime@${runtime}@g   ${srun}
      #------------------------------------------------------------------------------------#




      #----- Change the callserial.sh file. -----------------------------------------------#
      sed -i s@thisroot@${here}@g          ${callserial}
      sed -i s@thispoly@${polyname}@g      ${callserial}
      sed -i s@myexec@${execname}@g        ${callserial}
      sed -i s@myname@${moi}@g             ${callserial}
      sed -i s@mypackdata@${packdatasrc}@g ${callserial}
      sed -i s@myscenario@${iscenario}@g   ${callserial}
      sed -i s@myscenmain@${scentype}@g    ${callserial}
      sed -i s@myinitrc@${initrc}@g        ${callserial}
      #------------------------------------------------------------------------------------#


      #----- Re-submit the job to new queue. ----------------------------------------------#
      ${srun}
      #------------------------------------------------------------------------------------#
   fi # [ ${ispending} -gt 0 ]
   #---------------------------------------------------------------------------------------#
done # [ ${ff} -lt ${nfill} ]
#------------------------------------------------------------------------------------------#


#----- Get rid of the temporary file. -----------------------------------------------------#
/bin/rm -f ${pendfile}
#------------------------------------------------------------------------------------------#



#----- Sleep 5 more minutes. --------------------------------------------------------------#
echo "Taking another nap before checking the queues again..."
sleep 300
#------------------------------------------------------------------------------------------#



#----- Quick check the status of the queues. ----------------------------------------------#
/bin/rm -f ${recefile}
touch ${recefile}
echo "------- Polygon recent activity status. ----------------------------" >> ${recefile}
echo "Number of polygons that were re-submitted:    ${nresubmit}"           >> ${recefile}
echo "Number of polygons that switched queues:      ${nfill}"               >> ${recefile}
echo "Number of polygons that have been suspended:  ${newsuspend}"          >> ${recefile}
echo "Number of polygons that were unknown:         ${newunknown}"          >> ${recefile}
echo "Number of polygons that went extinct:         ${newextinct}"          >> ${recefile}
echo "Number of polygons that reached steady state: ${newststate}"          >> ${recefile}
echo "--------------------------------------------------------------------" >> ${recefile}
#------------------------------------------------------------------------------------------#




#----- Check the queue status. ------------------------------------------------------------#
echo "Counting the jobs..."
m61run=$(scount  -u ${moi} -g "${desc}-" -p moorcroft_6100 -t RUNNING)
amdrun=$(scount  -u ${moi} -g "${desc}-" -p moorcroft_amd  -t RUNNING)
urerun=$(scount  -u ${moi} -g "${desc}-" -p unrestricted   -t RUNNING)
genrun=$(scount  -u ${moi} -g "${desc}-" -p general        -t RUNNING)
m61pen=$(scount  -u ${moi} -g "${desc}-" -p moorcroft_6100 -t PENDING)
amdpen=$(scount  -u ${moi} -g "${desc}-" -p moorcroft_amd  -t PENDING)
urepen=$(scount  -u ${moi} -g "${desc}-" -p unrestricted   -t PENDING)
genpen=$(scount  -u ${moi} -g "${desc}-" -p general        -t PENDING)
m61rtot=$(scount -u ${moi}               -p moorcroft_6100 -t RUNNING)
amdrtot=$(scount -u ${moi}               -p moorcroft_amd  -t RUNNING)
urertot=$(scount -u ${moi}               -p unrestricted   -t RUNNING)
genrtot=$(scount -u ${moi}               -p general        -t RUNNING)
m61ptot=$(scount -u ${moi}               -p moorcroft_6100 -t PENDING)
amdptot=$(scount -u ${moi}               -p moorcroft_amd  -t PENDING)
ureptot=$(scount -u ${moi}               -p unrestricted   -t PENDING)
genptot=$(scount -u ${moi}               -p general        -t PENDING)
let m61tot=${m61rtot}+${m61ptot}
let amdtot=${amdrtot}+${amdptot}
let uretot=${urertot}+${ureptot}
let gentot=${genrtot}+${genptot}
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Make the table with the queue statuses.                                              #
#------------------------------------------------------------------------------------------#
/bin/rm -f ${queuefile}
touch ${queuefile}
echo "------- Queue status. --------------------------------------------" >> ${queuefile}
echo "Moorcroft_6100         RUN=${m61run}  PEN=${m61pen}  TOT=${m61tot}" >> ${queuefile}
echo "Moorcroft_AMD          RUN=${amdrun}  PEN=${amdpen}  TOT=${amdtot}" >> ${queuefile}
echo "Unrestricted           RUN=${urerun}  PEN=${urepen}  TOT=${uretot}" >> ${queuefile}
echo "General                RUN=${genrun}  PEN=${genpen}  TOT=${gentot}" >> ${queuefile}
echo "------------------------------------------------------------------" >> ${queuefile}
#------------------------------------------------------------------------------------------#



#----- Current simulation status. ---------------------------------------------------------#
n_initial=$(grep INITIAL ${outcheck} | wc -l)
n_history=$(grep HISTORY ${outcheck} | wc -l)
n_metmiss=$(grep METMISS ${outcheck} | wc -l)
n_bad_met=$(grep BAD_MET ${outcheck} | wc -l)
n_crashed=$(grep CRASHED ${outcheck} | wc -l)
n_stopped=$(grep STOPPED ${outcheck} | wc -l)
n_extinct=$(grep EXTINCT ${outcheck} | wc -l)
n_ststate=$(grep STSTATE ${outcheck} | wc -l)
n_the_end=$(grep THE_END ${outcheck} | wc -l)
/bin/rm -f ${statfile}
touch ${statfile}
echo "------- Simulation status. -----------------------------------------" >> ${statfile}
echo " Number of polygons that have never started           : ${n_initial}" >> ${statfile}
echo " Number of polygons that have partially run           : ${n_history}" >> ${statfile}
echo " Number of polygons that haven't found met drivers    : ${n_metmiss}" >> ${statfile}
echo " Number of polygons that have bad met drivers         : ${n_bad_met}" >> ${statfile}
echo " Number of polygons that have crashed                 : ${n_crashed}" >> ${statfile}
echo " Number of polygons that have mysteriously stopped    : ${n_stopped}" >> ${statfile}
echo " Number of polygons that became desert                : ${n_extinct}" >> ${statfile}
echo " Number of polygons that have reached steady state    : ${n_ststate}" >> ${statfile}
echo " Number of polygons that have reached the end         : ${n_the_end}" >> ${statfile}
echo "--------------------------------------------------------------------" >> ${statfile}
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Build the e-mail body.                                                               #
#------------------------------------------------------------------------------------------#
echo "Building the e-mail..."
/bin/rm -f ${emailbody}
touch ${emailbody}
cat ${headfile}  >> ${emailbody}
echo " "         >> ${emailbody}
cat ${statfile}  >> ${emailbody}
echo " "         >> ${emailbody}
cat ${recefile}  >> ${emailbody}
echo " "         >> ${emailbody}
cat ${queuefile} >> ${emailbody}
echo " "         >> ${emailbody}
cat ${tailfile}  >> ${emailbody}
echo " "         >> ${emailbody}
#----- Check whether to append some plots. ------------------------------------------------#
if [ ${plotstatus} -eq 1 ]
then
   attach=""
   for fichier in ${R_figlist}
   do
      if [ -s ${fichier} ]
      then
         attach="${attach} -a ${fichier}"
      fi
   done
else
   attach=""
fi
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#     Create the subject for the e-mail.                                                   #
#------------------------------------------------------------------------------------------#
when=$(date +'%d %B %Y - %R %Z')
today=$(date +'%Y-%m-%d')
subject="${runtitle} run status as of ${when}"
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Check whether to send the e-mail or not.                                             #
#------------------------------------------------------------------------------------------#
if [ ${email1day} -eq 0 ]
then
  #----- Always send e-mail. --------------------------------------------------------------#
  sendemail="y"
  #----------------------------------------------------------------------------------------#
elif [ -s ${emailday} ]
then
   #----- E-mail has been sent.  Check whether it is time to send again. ------------------#
   lastemail=$(cat ${emailday})
   if [ ${today} != ${lastemail} ]
   then
      #----- Time to send another e-mail. -------------------------------------------------#
      sendemail="y"
      #------------------------------------------------------------------------------------#
   else
      #----- E-mail has been sent recently.  Don't send another one now. ------------------#
      sendemail="n"
      #------------------------------------------------------------------------------------#
   fi
   #---------------------------------------------------------------------------------------#
else
   #----- E-mail has not been sent yet.  Send it. -----------------------------------------#
   sendemail="y"
   #---------------------------------------------------------------------------------------#
fi
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#      Send/skip the e-mail.                                                               #
#------------------------------------------------------------------------------------------#
if [ ${sendemail} == "y" ]
then
   echo ${today} > ${emailday}
   echo "Sending e-mail..."
   ${mailprog} -s "${subject}" ${attach} ${recipient} < ${emailbody}
else
   echo "Skipping e-mail..."
fi
#------------------------------------------------------------------------------------------#

#----- Clean-up stuff. --------------------------------------------------------------------#
echo "Deleting some temporary files..."
/bin/rm -f ${queuefile} ${recefile}
/bin/rm -f ${here}/run_sitter.lock
echo "==== run_sitter.sh execution ends. ===="
