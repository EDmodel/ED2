#!/bin/bash

#==========================================================================================#
#==========================================================================================#
#     Main settings:                                                                       #
#------------------------------------------------------------------------------------------#
#----- Main path, usually set by $(pwd) so you don't need to change it. -------------------#
here=$(pwd)
#----- User name, usually set by $(whoami) so you don't need to change it. ----------------#
moi=$(whoami)
#----- Description of this simulation, used to create unique job names. -------------------#
desc=$(basename ${here})
#----- File with job schedule, normally ${here}/joborder.txt. -----------------------------#
joborder="${here}/joborder.txt"
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#     Sometimes you may want to print more than one column at a time.  In this case        #
# provide the number of columns.                                                           #
#------------------------------------------------------------------------------------------#
if [ "x${1}" == "x" ]
then
   ncol=1
else
   ncol=${1}
fi
#------------------------------------------------------------------------------------------#



#----- Determine the number of polygons to run. -------------------------------------------#
let npolys=$(wc -l ${joborder} | awk '{print $1 }')-3
echo "Number of polygons: ${npolys}..."
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#     Loop over all polygons.                                                              #
#------------------------------------------------------------------------------------------#
ff=0
while [ ${ff} -lt ${npolys} ]
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



   #----- Make two columns. ---------------------------------------------------------------#
   let col=${ff}%${ncol}
   if [ ${ncol} -eq 1 ]
   then
      opt=""
      off=""
   elif [ ${col} -eq 0 ]
   then
      opt=""
      off=".\t"
   else
      opt="-n"
      off=""
   fi

   #---------------------------------------------------------------------------------------#
   #      Read the ffth line of the polygon list.  There must be smarter ways of doing     #
   # this, but this works.  Here we obtain the polygon name, and its longitude and         #
   # latitude.                                                                             #
   #---------------------------------------------------------------------------------------#
   oi=$(head -${line} ${joborder} | tail -1)
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
   ihrzrad=$(echo ${oi}      | awk '{print $77}')
   crown=$(echo   ${oi}      | awk '{print $78}')
   ltransvis=$(echo ${oi}    | awk '{print $79}')
   lreflectvis=$(echo ${oi}  | awk '{print $80}')
   ltransnir=$(echo ${oi}    | awk '{print $81}')
   lreflectnir=$(echo ${oi}  | awk '{print $82}')
   orienttree=$(echo ${oi}   | awk '{print $83}')
   orientgrass=$(echo ${oi}  | awk '{print $84}')
   clumptree=$(echo ${oi}    | awk '{print $85}')
   clumpgrass=$(echo ${oi}   | awk '{print $86}')
   ivegtdyn=$(echo ${oi}     | awk '{print $87}')
   igndvap=$(echo ${oi}      | awk '{print $88}')
   iphen=$(echo ${oi}        | awk '{print $89}')
   iallom=$(echo ${oi}       | awk '{print $90}')
   ibigleaf=$(echo ${oi}     | awk '{print $91}')
   irepro=$(echo ${oi}       | awk '{print $92}')
   treefall=$(echo ${oi}     | awk '{print $93}')
   ianthdisturb=$(echo ${oi} | awk '{print $94}')
   ianthdataset=$(echo ${oi} | awk '{print $95}')
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Set some variables to check whether the simulation is running.                    #
   #---------------------------------------------------------------------------------------#
   jobname="${desc}-${polyname}"
   stdout="${here}/${polyname}/serial_out.out"
   stderr="${here}/${polyname}/serial_out.err"
   lsfout="${here}/${polyname}/serial_lsf.out"
   skipper="${here}/${polyname}/skipper.txt"
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Check whether the simulation is still running, and if not, why it isn't.          #
   #---------------------------------------------------------------------------------------#
   if [ -s ${stdout} ]
   then
      #----- Check whether the simulation is running, and when in model time it is. -------#
      running=$(qstat -j ${jobname} -s r 2> /dev/null | wc -l)
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
      stopped=$(grep "FATAL ERROR"           ${stdout} | wc -l)
      crashed=$(grep "IFLAG1 problem."       ${stdout} | wc -l)
      the_end=$(grep "ED-2.2 execution ends" ${stdout} | wc -l)
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Plot a message so the user knows what is going on.                             #
      #------------------------------------------------------------------------------------#
      if [ ${running} -gt 0 ] || [ -s ${skipper} ] && [ ${sigsegv} -eq 0 ]
      then
         echo -e ${opt} "${off} ${ffout}: ${polyname} is running (${runtime})..."
      elif [ ${sigsegv} -gt 0 ]
      then
         echo -e ${opt} "${off} ${ffout}: ${polyname} HAD SEGMENTATION VIOLATION... <===="
      elif [ ${crashed} -gt 0 ]
      then 
         echo -e ${opt} "${off} ${ffout}: ${polyname} HAS CRASHED (RK4 PROBLEM)... <====="
      elif [ ${metmiss} -gt 0 ]
      then 
         echo -e ${opt} "${off} ${ffout}: ${polyname} DID NOT FIND MET DRIVERS... <======"
      elif [ ${stopped} -gt 0 ]
      then
         echo -e ${opt} "${off} ${ffout}: ${polyname} STOPPED (UNKNOWN REASON)... <======"
      elif [ ${the_end} -gt 0 ]
      then
         echo -e ${opt} "${off} ${ffout}: ${polyname} has finished o/\o..."
      else
         echo -e ${opt} "${off} ${ffout}: ${polyname} status is unknown..."
      fi
      #------------------------------------------------------------------------------------#
   else
      echo -e ${opt} "${off} ${ffout}: ${polyname} is pending ..."
   fi
   #---------------------------------------------------------------------------------------#
done
#------------------------------------------------------------------------------------------#

