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
here="/n/regal/moorcroft_lab/mlongo/pve+sas"
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
runtitle="ED-2 potential vegetation simulation"
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


#------ Calculator. -----------------------------------------------------------------------#
ccc="${HOME}/util/calc.sh"  # Calculator
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#      Variables that tell whether to check for steady state and deserts.                  #
#------------------------------------------------------------------------------------------#
checksteady="TRUE"    #   Check whether to stop simulations? (TRUE of FALSE, R style)
nyearmin=160          #   Minimum number of years before we check
ststcrit=0.01         #   Maximum change allowed between two cycles
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
recipient="mdplongo@gmail.com"
mailprog="${HOME}/util/mutt"
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
wait_minutes=120
frqpost=0
frqtouch=5
#------------------------------------------------------------------------------------------#


#----- Format for checking simulations. ---------------------------------------------------#
outform="JobName%200,State%12"
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
if [ "x${here}" == "x" ] || [ "x${recipient}" == "x" ]
then
   echo " + Set variables \"here\" and \"recipient\" in your script before running!"
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



#------ Move to the current directory. ----------------------------------------------------#
cd ${main}
#------------------------------------------------------------------------------------------#


#----- Determine the number of polygons to run. -------------------------------------------#
let n_polygon=$(wc -l ${joborder} | awk '{print $1 }')-3
echo "Number of polygons: ${n_polygon}..."
#------------------------------------------------------------------------------------------#

desc=$(basename ${here})
jobname="${desc}-sims"


#------------------------------------------------------------------------------------------#
#       Main loop: the script stays inside until all simulations have finished.            #
#------------------------------------------------------------------------------------------#
n_ongoing=9999
iter=0
while [ ${n_ongoing} -gt 0 ]
do
   #------ Update iteration counter. ------------------------------------------------------#
   let iter=${iter}+1
   #---------------------------------------------------------------------------------------#

   #---- Reset count of polygons. ---------------------------------------------------------#
   n_running=0
   n_suspend=0
   n_the_end=0
   n_unknown=0
   n_sigsegv=0
   n_crashed=0
   n_metmiss=0
   n_stopped=0
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Check whether it's time to touch files.                                           #
   #---------------------------------------------------------------------------------------#
   if [ ${frqtouch} -gt 0 ]
   then
      let xtouch=${iter}%${frqtouch}
      if [ ${xtouch} -eq 0 ]
      then
         echo " Touch files."
         touch ${here}/*
         touch ${here}/Template/*
         touch ${here}/executable/*
         touch ${here}/sit_utils/*
      fi
   else
      xtouch=999
   fi
   #---------------------------------------------------------------------------------------#


   rm -f ${lastcheck}
   mv ${outcheck} ${lastcheck}
   touch ${outcheck}
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Go through all polygons.                                                          #
   #---------------------------------------------------------------------------------------#
   ff=0
   while [ ${ff} -lt ${n_polygon} ]
   do
      let ff=${ff}+1
      let line=${ff}+3

      #------------------------------------------------------------------------------------#
      #    Format count.                                                                   #
      #------------------------------------------------------------------------------------#
      if   [ ${n_polygon} -ge 10   ] && [ ${n_polygon} -lt 100   ]
      then
         ffout=$(printf '%2.2i' ${ff})
      elif [ ${n_polygon} -ge 100  ] && [ ${n_polygon} -lt 1000  ]
      then
         ffout=$(printf '%2.2i' ${ff})
      elif [ ${n_polygon} -ge 100  ] && [ ${n_polygon} -lt 10000 ]
      then
         ffout=$(printf '%2.2i' ${ff})
      else
         ffout=${ff}
      fi
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #      Read the ffth line of the polygon list.  There must be smarter ways of doing  #
      # this, but this works.  Here we obtain the polygon name, and its longitude and      #
      # latitude.                                                                          #
      #------------------------------------------------------------------------------------#
      oi=$(head -${line} ${joborder} | tail -1)
      polyname=$(echo ${oi}     | awk '{print $1  }')
      polyiata=$(echo ${oi}     | awk '{print $2  }')
      polylon=$(echo ${oi}      | awk '{print $3  }')
      polylat=$(echo ${oi}      | awk '{print $4  }')
      yeara=$(echo ${oi}        | awk '{print $5  }')
      montha=$(echo ${oi}       | awk '{print $6  }')
      datea=$(echo ${oi}        | awk '{print $7  }')
      timea=$(echo ${oi}        | awk '{print $8  }')
      yearz=$(echo ${oi}        | awk '{print $9  }')
      monthz=$(echo ${oi}       | awk '{print $10 }')
      datez=$(echo ${oi}        | awk '{print $11 }')
      timez=$(echo ${oi}        | awk '{print $12 }')
      initmode=$(echo ${oi}     | awk '{print $13 }')
      iscenario=$(echo ${oi}    | awk '{print $14 }')
      isizepft=$(echo ${oi}     | awk '{print $15 }')
      iage=$(echo ${oi}         | awk '{print $16 }')
      imaxcohort=$(echo ${oi}   | awk '{print $17 }')
      polyisoil=$(echo ${oi}    | awk '{print $18 }')
      polyntext=$(echo ${oi}    | awk '{print $19 }')
      polysand=$(echo ${oi}     | awk '{print $20 }')
      polyclay=$(echo ${oi}     | awk '{print $21 }')
      polydepth=$(echo ${oi}    | awk '{print $22 }')
      polysoilbc=$(echo ${oi}   | awk '{print $23 }')
      polysldrain=$(echo ${oi}  | awk '{print $24 }')
      polycol=$(echo ${oi}      | awk '{print $25 }')
      slzres=$(echo ${oi}       | awk '{print $26 }')
      queue=$(echo ${oi}        | awk '{print $27 }')
      metdriver=$(echo ${oi}    | awk '{print $28 }')
      dtlsm=$(echo ${oi}        | awk '{print $29 }')
      monyrstep=$(echo ${oi}    | awk '{print $30 }')
      vmfactc3=$(echo ${oi}     | awk '{print $31 }')
      vmfactc4=$(echo ${oi}     | awk '{print $32 }')
      mphototrc3=$(echo ${oi}   | awk '{print $33 }')
      mphototec3=$(echo ${oi}   | awk '{print $34 }')
      mphotoc4=$(echo ${oi}     | awk '{print $35 }')
      bphotoblc3=$(echo ${oi}   | awk '{print $36 }')
      bphotonlc3=$(echo ${oi}   | awk '{print $37 }')
      bphotoc4=$(echo ${oi}     | awk '{print $38 }')
      kwgrass=$(echo ${oi}      | awk '{print $39 }')
      kwtree=$(echo ${oi}       | awk '{print $40 }')
      gammac3=$(echo ${oi}      | awk '{print $41 }')
      gammac4=$(echo ${oi}      | awk '{print $42 }')
      d0grass=$(echo ${oi}      | awk '{print $43 }')
      d0tree=$(echo ${oi}       | awk '{print $44 }')
      alphac3=$(echo ${oi}      | awk '{print $45 }')
      alphac4=$(echo ${oi}      | awk '{print $46 }')
      klowco2=$(echo ${oi}      | awk '{print $47 }')
      decomp=$(echo ${oi}       | awk '{print $48 }')
      rrffact=$(echo ${oi}      | awk '{print $49 }')
      growthresp=$(echo ${oi}   | awk '{print $50 }')
      lwidthgrass=$(echo ${oi}  | awk '{print $51 }')
      lwidthbltree=$(echo ${oi} | awk '{print $52 }')
      lwidthnltree=$(echo ${oi} | awk '{print $53 }')
      q10c3=$(echo ${oi}        | awk '{print $54 }')
      q10c4=$(echo ${oi}        | awk '{print $55 }')
      h2olimit=$(echo ${oi}     | awk '{print $56 }')
      imortscheme=$(echo ${oi}  | awk '{print $57 }')
      ddmortconst=$(echo ${oi}  | awk '{print $58 }')
      cbrscheme=$(echo ${oi}    | awk '{print $59 }')
      isfclyrm=$(echo ${oi}     | awk '{print $60 }')
      icanturb=$(echo ${oi}     | awk '{print $61 }')
      ubmin=$(echo ${oi}        | awk '{print $62 }')
      ugbmin=$(echo ${oi}       | awk '{print $63 }')
      ustmin=$(echo ${oi}       | awk '{print $64 }')
      gamm=$(echo ${oi}         | awk '{print $65 }')
      gamh=$(echo ${oi}         | awk '{print $66 }')
      tprandtl=$(echo ${oi}     | awk '{print $67 }')
      ribmax=$(echo ${oi}       | awk '{print $68 }')
      atmco2=$(echo ${oi}       | awk '{print $69 }')
      thcrit=$(echo ${oi}       | awk '{print $70 }')
      smfire=$(echo ${oi}       | awk '{print $71 }')
      ifire=$(echo ${oi}        | awk '{print $72 }')
      fireparm=$(echo ${oi}     | awk '{print $73 }')
      ipercol=$(echo ${oi}      | awk '{print $74 }')
      runoff=$(echo ${oi}       | awk '{print $75 }')
      imetrad=$(echo ${oi}      | awk '{print $76 }')
      ibranch=$(echo ${oi}      | awk '{print $77 }')
      icanrad=$(echo ${oi}      | awk '{print $78 }')
      ihrzrad=$(echo ${oi}      | awk '{print $79 }')
      crown=$(echo   ${oi}      | awk '{print $80 }')
      ltransvis=$(echo ${oi}    | awk '{print $81 }')
      lreflectvis=$(echo ${oi}  | awk '{print $82 }')
      ltransnir=$(echo ${oi}    | awk '{print $83 }')
      lreflectnir=$(echo ${oi}  | awk '{print $84 }')
      orienttree=$(echo ${oi}   | awk '{print $85 }')
      orientgrass=$(echo ${oi}  | awk '{print $86 }')
      clumptree=$(echo ${oi}    | awk '{print $87 }')
      clumpgrass=$(echo ${oi}   | awk '{print $88 }')
      igoutput=$(echo ${oi}     | awk '{print $89 }')
      ivegtdyn=$(echo ${oi}     | awk '{print $90 }')
      igndvap=$(echo ${oi}      | awk '{print $91 }')
      iphen=$(echo ${oi}        | awk '{print $92 }')
      iallom=$(echo ${oi}       | awk '{print $93 }')
      ibigleaf=$(echo ${oi}     | awk '{print $94 }')
      integscheme=$(echo ${oi}  | awk '{print $95 }')
      nsubeuler=$(echo ${oi}    | awk '{print $96 }')
      irepro=$(echo ${oi}       | awk '{print $97 }')
      treefall=$(echo ${oi}     | awk '{print $98 }')
      ianthdisturb=$(echo ${oi} | awk '{print $99 }')
      ianthdataset=$(echo ${oi} | awk '{print $100}')
      slscale=$(echo ${oi}      | awk '{print $101}')
      slyrfirst=$(echo ${oi}    | awk '{print $102}')
      slnyrs=$(echo ${oi}       | awk '{print $103}')
      bioharv=$(echo ${oi}      | awk '{print $104}')
      skidarea=$(echo ${oi}     | awk '{print $105}')
      skidsmall=$(echo ${oi}    | awk '{print $106}')
      skidlarge=$(echo ${oi}    | awk '{print $107}')
      fellingsmall=$(echo ${oi} | awk '{print $108}')
      #------------------------------------------------------------------------------------#



      #----- Check whether the directories exist or not. ----------------------------------#
      if [ ${xtouch} -eq 0 ]
      then
         touch ${here}/${polyname}/analy/*
         touch ${here}/${polyname}/histo/*
         touch ${here}/${polyname}/*
      fi
      #---------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Set some variables to check whether the simulation is running.                 #
      #------------------------------------------------------------------------------------#
      stdout="${here}/${polyname}/serial_out.out"
      stderr="${here}/${polyname}/serial_out.err"
      lsfout="${here}/${polyname}/serial_lsf.out"
      #------------------------------------------------------------------------------------#



      #----- Update the R script that controls the run. --------------------------------------#
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
      #---------------------------------------------------------------------------------------#



      #---------------------------------------------------------------------------------------#
      #    Retrieve the information on statusrun.txt because it may tell whether the run was  #
      # in steady state or all PFTs went extinct, and we don't need to check them again.      #
      #---------------------------------------------------------------------------------------#
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
      #---------------------------------------------------------------------------------------#




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
      #------------------------------------------------------------------------------------#



      #----- Write polygon check into a single table. -------------------------------------#
      output="${polyname} ${polylon} ${polylat} ${yearh} ${monthh} ${dateh} ${timeh}"
      output="${output} ${runt} ${agb} ${bsa} ${lai} ${scb} ${npa}"
      echo ${output} >> ${outcheck}
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Check whether the simulation is still running, and if not, why it isn't.       #
      #------------------------------------------------------------------------------------#
      if [ -s ${stdout} ]
      then
         #----- Check whether the simulation is running, and when in model time it is. ----#
         stask="stask --noheader -u $(whoami) -t ${polyname} -j ${jobname}"
         running=$(${stask}   -o "${outform}" | grep "RUNNING"   | wc -l)
         pending=$(${stask}   -o "${outform}" | grep "PENDING"   | wc -l)
         suspended=$(${stask} -o "${outform}" | grep "SUSPENDED" | wc -l)
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



         #----- Check for other possible outcomes. ----------------------------------------#
         stopped=$(grep "FATAL ERROR"           ${stdout} | wc -l)
         crashed=$(grep "IFLAG1 problem."       ${stdout} | wc -l)
         the_end=$(grep "ED-2.2 execution ends" ${stdout} | wc -l)
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #     Plot a message so the user knows what is going on.                          #
         #---------------------------------------------------------------------------------#
         if [ ${pending} -gt 0 ]
         then
            let n_pending=${n_pending}+1
            echo -e "${ffout}: ${polyname} is pending..."
         elif [ ${suspended} -gt 0 ]
         then
            let n_suspend=${n_suspend}+1
            echo -e "${ffout}: ${polyname} is suspended!!!"
         elif [ ${running} -gt 0 ] || [ -s ${skipper} ] && [ ${sigsegv} -eq 0 ]
         then
            let n_running=${n_running}+1
            echo -e "${ffout}: ${polyname} is running (${runtime})..."
         elif [ ${sigsegv} -gt 0 ]
         then
            let n_sigsegv=${n_sigsegv}+1
            echo -e "${ffout}: ${polyname} HAD SEGMENTATION VIOLATION... <==========="
         elif [ ${crashed} -gt 0 ]
         then 
            let n_crashed=${n_crashed}+1
            echo -e "${ffout}: ${polyname} HAS CRASHED (RK4 PROBLEM)... <==========="
         elif [ ${metmiss} -gt 0 ]
         then 
            let n_metmiss=${n_metmiss}+1
            echo -e "${ffout}: ${polyname} DID NOT FIND MET DRIVERS... <==========="
         elif [ ${stopped} -gt 0 ]
         then
            let n_stopped=${n_stopped}+1
            echo -e "${ffout}: ${polyname} STOPPED (UNKNOWN REASON)... <==========="
         elif [ ${the_end} -gt 0 ]
         then
            let n_the_end=${n_the_end}+1
            echo -e "${ffout}: ${polyname} has finished o/\o..."
         else
            let n_unknown=${n_unknown}+1
            echo -e "${ffout}: ${polyname} status is unknown..."
         fi
         #---------------------------------------------------------------------------------#
      else
         let n_pending=${n_pending}+1
         echo -e "${ffout}: ${polyname} is pending ..."
      fi
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Update the history time for ED2IN so it doesn't start from the beginning in    #
      # case the job is resubmitted automatically.                                         #
      #------------------------------------------------------------------------------------#
      ED2IN="${here}/${polyname}/ED2IN"
      runtype_new="   NL%RUNTYPE = 'HISTORY'"
      sfilin_new="   NL%SFILIN = '${here}/${polyname}/histo/${polyname}'"
      itimeh_new="   NL%ITIMEH = ${timeh}"
      idateh_new="   NL%IDATEH = ${dateh}"
      imonthh_new="   NL%IMONTHH = ${monthh}"
      iyearh_new="   NL%IYEARH = ${yearh}"
      runtype_old=$(grep  -i "NL%RUNTYPE" ${ED2IN} | grep -v "\!")
      sfilin_old=$(grep  -i "NL%SFILIN"   ${ED2IN} | grep -v "\!")
      itimeh_old=$(grep  -i "NL%ITIMEH"   ${ED2IN} )
      idateh_old=$(grep  -i "NL%IDATEH"   ${ED2IN} )
      imonthh_old=$(grep -i "NL%IMONTHH"  ${ED2IN} )
      iyearh_old=$(grep -i  "NL%IYEARH"   ${ED2IN} )
      sed -i~ s@"${runtype_old}"@"${runtype_new}"@g ${ED2IN}
      sed -i~ s@"${sfilin_old}"@"${sfilin_new}"@g   ${ED2IN}
      sed -i~ s@"${itimeh_old}"@"${itimeh_new}"@g   ${ED2IN}
      sed -i~ s@"${idateh_old}"@"${idateh_new}"@g   ${ED2IN}
      sed -i~ s@"${imonthh_old}"@"${imonthh_new}"@g ${ED2IN}
      sed -i~ s@"${iyearh_old}"@"${iyearh_new}"@g   ${ED2IN}
      #------------------------------------------------------------------------------------#
   done
   #---------------------------------------------------------------------------------------#



   #----- Run R to make the status check. -------------------------------------------------#
   if ${plotstatus}
   then
      echo "     + Run the status check."
      Rscript_out="$(dirname ${Rscript_plot})/$(basename ${Rscript_plot} .r).txt"
      R CMD BATCH --no-save --no-restore ${Rscript_plot} ${Rscript_out}
   fi
   #---------------------------------------------------------------------------------------#



   #----- Current simulation status. ------------------------------------------------------#
   n_initial=$(grep INITIAL ${outcheck} | wc -l)
   n_bad_met=$(grep BAD_MET ${outcheck} | wc -l)
   n_extinct=$(grep EXTINCT ${outcheck} | wc -l)
   n_ststate=$(grep STSTATE ${outcheck} | wc -l)
   let n_ongoing=${n_polygon}-${n_the_end}-${n_extinct}-${n_ststate}
   /bin/rm -f ${statfile}
   touch ${statfile}
   echo "------- Simulation status. --------------------------------------" >> ${statfile}
   echo " Number of polygons that have never started        : ${n_initial}" >> ${statfile}
   echo " Number of polygons that have partially run        : ${n_running}" >> ${statfile}
   echo " Number of polygons that haven't found met drivers : ${n_metmiss}" >> ${statfile}
   echo " Number of polygons that have bad met drivers      : ${n_bad_met}" >> ${statfile}
   echo " Number of polygons that have crashed              : ${n_crashed}" >> ${statfile}
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
   cat ${recefile}  >> ${emailbody}
   echo " "         >> ${emailbody}
   cat ${queuefile} >> ${emailbody}
   echo " "         >> ${emailbody}
   cat ${tailfile}  >> ${emailbody}
   echo " "         >> ${emailbody}
   #----- Check whether to append some plots. ---------------------------------------------#
   if ${plotstatus}
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
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Create the subject for the e-mail.                                                #
   #---------------------------------------------------------------------------------------#
   when=$(date +'%d %B %Y - %R %Z')
   today=$(date +'%Y-%m-%d')
   subject="${runtitle} run status as of ${when}"
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Check whether to send the e-mail or not.                                          #
   #---------------------------------------------------------------------------------------#
   if [ ${email1day} -eq 0 ]
   then
     #----- Always send e-mail. -----------------------------------------------------------#
     sendemail=true
     #-------------------------------------------------------------------------------------#
   elif [ -s ${emailday} ]
   then
      #----- E-mail has been sent.  Check whether it is time to send again. ---------------#
      lastemail=$(cat ${emailday})
      if [ ${today} != ${lastemail} ]
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
      echo ${today} > ${emailday}
      echo "     + Send e-mail."
      ${mailprog} -s "${subject}" ${attach} ${recipient} < ${emailbody}
   else
      echo "     + Skip e-mail."
   fi
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      Run the post processing.                                                         #
   #---------------------------------------------------------------------------------------#
   if [ ${frqpost} -gt 0 ]
   then
      let xpost=${iter}%${frqpost}
      if [ -s "${here}/epost.sh" ] && [ ${xpost} -eq 0 -o ${n_ongoing} -eq 0 ]
      then
         echo " Run post-processing."
         ${here}/epost.sh


         echo " Compress/clean files."
         ${here}/last_histo.sh
      fi
   fi
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #       In case any simulation is still on the works, take a break before checking      #
   #   again.                                                                              #
   #---------------------------------------------------------------------------------------#
   if [ ${n_ongoing} -gt 0 ]
   then
      echo "     + Take a break before checking again."
      let nsl=${wait_minutes}+1
      while [ ${nsl} -gt 1 ]
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


   #----- Clean-up stuff. -----------------------------------------------------------------#
   echo "     + Delete some temporary files..."
   /bin/rm -f ${queuefile} ${recefile}
   #---------------------------------------------------------------------------------------#
done
#------------------------------------------------------------------------------------------#

#----- Clean-up stuff. --------------------------------------------------------------------#
echo " Unlock run_sitter.sh."
/bin/rm -f ${here}/run_sitter.lock
echo "==== run_sitter.sh execution ends. ===="
#------------------------------------------------------------------------------------------#
