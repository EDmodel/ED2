#!/bin/bash
here=$(pwd)
joborder="${here}/joborder.txt"
desc=$(basename ${here})
jobname="${desc}-sims*"
moi=$(whoami)
outform="JobName%200,State%12"
#----- Determine the number of polygons to run. -------------------------------------------#
let npolys=$(wc -l ${joborder} | awk '{print $1 }')-3
echo "Number of polygons: ${npolys}..."




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



polya=1
polyz=${npolys}

#----- Argument parsing. ------------------------------------------------------------------#
while [[ ${#} > 0 ]]
do
key="${1}"
   case ${key} in
   -a)
      polya="${2}"
      shift
      ;;
   -z)
      polyz="${2}"
      shift
      ;;
   --polya=*)
      polya=$(echo ${key} | sed s@"-polya=="@""@g)
      ;;
   --polyz=*)
      polyz=$(echo ${key} | sed s@"-polyz=="@""@g)
      ;;
   *)
      echo "Unknown key-value argument pair."
      exit 2
      ;;
   esac

   shift # past argument or value
done
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#     Loop over all polygons.                                                              #
#------------------------------------------------------------------------------------------#
let ff=${polya}-1
while [ ${ff} -lt ${polyz} ]
do
   let ff=${ff}+1
   let line=${ff}+3

   #---------------------------------------------------------------------------------------#
   #    Format count.                                                                      #
   #---------------------------------------------------------------------------------------#
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
      ;;
   esac
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Set some variables to check whether the simulation is running.                    #
   #---------------------------------------------------------------------------------------#
   stdout="${here}/${polyname}/serial_out.out"
   stderr="${here}/${polyname}/serial_out.err"
   skipper="${here}/${polyname}/skipper.txt"
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Check whether the simulation is still running, and if not, why it isn't.          #
   #---------------------------------------------------------------------------------------#
   if [ -s ${stdout} ]
   then

      #----- Set task inquire command. ----------------------------------------------------#
      case "${cluster}" in
      SDUMONT)
         #----- Multi-task job, search for specific task. ---------------------------------#
         stask="stask --noheader -u ${moi} -t ${taskname} -j ${jobname}"
         #---------------------------------------------------------------------------------#
         ;;
      CANNON)
         #----- Single-task jobs, search the task using the job name. ---------------------#
         stask="stask --noheader -u ${moi} -j ${taskname}"
         #---------------------------------------------------------------------------------#
         ;;
      esac
      #------------------------------------------------------------------------------------#


      #----- Check whether the simulation is running, and when in model time it is. -------#
      running=$(${stask}   -o "${outform}" | grep "RUNNING"   | wc -l)
      pending=$(${stask}   -o "${outform}" | grep "PENDING"   | wc -l)
      suspended=$(${stask} -o "${outform}" | grep "SUSPENDED" | wc -l)
      simline=$(grep "Simulating: "   ${stdout} | tail -1)
      runtime=$(echo ${simline} | awk '{print $3}')
      #------------------------------------------------------------------------------------#



      #----- Check for segmentation violations. -------------------------------------------#
      if [ -s ${stderr} ]
      then
         segv1=$(grep -i "sigsegv"            ${stderr} | wc -l)
         segv2=$(grep -i "segmentation fault" ${stderr} | wc -l)
         let sigsegv=${segv1}+${segv2}
      else
         sigsegv=0
      fi
      #------------------------------------------------------------------------------------#



      #----- Check whether met files are missing... (bad start) ---------------------------#
      metbs1=$(grep "Cannot open met driver input file" ${stdout} | wc -l)
      metbs2=$(grep "Specify ED_MET_DRIVER_DB properly" ${stdout} | wc -l)
      let metmiss=${metbs1}+${metbs2}
      #------------------------------------------------------------------------------------#



      #----- Check for other possible outcomes. -------------------------------------------#
      stopped=$(grep "FATAL ERROR"                       ${stdout} | wc -l)
      crashed=$(grep "IFLAG1 problem."                   ${stdout} | wc -l)
      hydfail=$(grep "Plant Hydrodynamics is off-track." ${stdout} | wc -l)
      the_end=$(grep "ED-2.2 execution ends"             ${stdout} | wc -l)
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Plot a message so the user knows what is going on.                             #
      #------------------------------------------------------------------------------------#
      if [ ${pending} -gt 0 ]
      then
         echo -e "${ffout}: ${polyname} is pending..."
      elif [ ${suspended} -gt 0 ]
      then
         echo -e "${ffout}: ${polyname} is suspended!!!"
      elif [ ${running} -gt 0 ] || [ -s ${skipper} ] && [ ${sigsegv} -eq 0 ]
      then
         echo -e "${ffout}: ${polyname} is running (${runtime})..."
      elif [ ${sigsegv} -gt 0 ]
      then
         echo -e "${ffout}: ${polyname} HAD SEGMENTATION VIOLATION... <==========="
      elif [ ${crashed} -gt 0 ]
      then 
         echo -e "${ffout}: ${polyname} HAS CRASHED (RK4 PROBLEM)... <==========="
      elif [ ${hydfail} -gt 0 ]
      then 
         echo -e "${ffout}: ${polyname} HAS CRASHED (HYDRODYNAMICS)... <==========="
      elif [ ${metmiss} -gt 0 ]
      then 
         echo -e "${ffout}: ${polyname} DID NOT FIND MET DRIVERS... <==========="
      elif [ ${stopped} -gt 0 ]
      then
         echo -e "${ffout}: ${polyname} STOPPED (UNKNOWN REASON)... <==========="
      elif [ ${the_end} -gt 0 ]
      then
         echo -e "${ffout}: ${polyname} has finished o/\o..."
      else
         echo -e "${ffout}: ${polyname} status is unknown..."
      fi
      #------------------------------------------------------------------------------------#
   else
      echo -e "${ffout}: ${polyname} is pending ..."
   fi
   #---------------------------------------------------------------------------------------#
done
#------------------------------------------------------------------------------------------#

