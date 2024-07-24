#!/bin/sh
here=$(pwd)
moi=$(whoami)
diskthere="/n/moorcroftfs2"
joborder="${here}/joborder.txt"

desc=$(basename ${here})

#----- Executable name. -------------------------------------------------------------------#
execname="ed_2.2-opt"
execsrc="${HOME}/EDBRAMS/ED/build"
#------------------------------------------------------------------------------------------#


#----- Find out which platform we are using. ----------------------------------------------#
host=$(hostname -s)
case ${host} in
rclogin*|holy*|moorcroft*|rcnx*|sdumont*)
   platform="SLURM"
   ;;
au*|ha*|sun-master|cmm*)
   platform="PBS"
   ;;
*)
   echo -n "Failed guessing platform from node name.  Please type the name:   "
   read platform
   ;;
esac
#------------------------------------------------------------------------------------------#


#----- Find the output path (both local and remote paths will be cleaned). ----------------#
basehere=$(basename ${here})
dirhere=$(dirname ${here})
while [ ${basehere} != ${moi} ]
do
   basehere=$(basename ${dirhere})
   dirhere=$(dirname ${dirhere})
done
diskhere=${dirhere}
echo "-------------------------------------------------------------------------------"
echo " - Simulation control on disk: ${diskhere}"
echo " - Output on disk:             ${diskthere}"
echo "-------------------------------------------------------------------------------"
there=$(echo ${here} | sed s@${diskhere}@${diskthere}@g)
#------------------------------------------------------------------------------------------#


#----- Determine the number of polygons to run. -------------------------------------------#
let npolys=$(wc -l ${joborder} | awk '{print $1 }')-3
#------------------------------------------------------------------------------------------#


#----- Check that the user is aware that it will remove everything... ---------------------#
if [ "x${1}" == "x-d" ]
then
   echo "Are you sure you want to stop all jobs, and remove all files and folders? [y/N]"
else
   echo "Are you sure you want to stop all jobs, and remove all files? [y/N]"
fi
read proceed

if [ "x${proceed}" != "xy" ] && [ "x${proceed}" != "xY" ]
then
   exit
fi

echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
echo " "
echo "     Look, this will really stop ALL your jobs and delete all files!!!"
echo "     Are you sure? [y/N]"
echo " "
echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
read proceed

if [ "x${proceed}" != "xy" ] && [ "x${proceed}" != "xY" ]
then
   exit
fi

echo "Okay then, but if you regret later don't say that I did not warn you..."
echo "I'm giving you a few seconds to kill this script in case you change your mind..."
delfun=11
while [ ${delfun} -gt 1 ]
do
   let delfun=${delfun}-1
   echo "  - Job stopping will begin in ${delfun} seconds..."
   sleep 1
done
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


   #------ Define job name. ---------------------------------------------------------------#
   jobname="${desc}-${polyname}"
   #---------------------------------------------------------------------------------------#


   #------- Delete jobs. ------------------------------------------------------------------#
   case "${platform}" in
   SLURM)
      scancel -n ${jobname}
      ;;
   PBS)
      jobid=$(qjobs -j ${jobname} -n | awk '{print $1}')
      if [[ "${jobid}" != "" ]]
      then
         qdel ${jobid}
      fi
      ;;
   esac
   #---------------------------------------------------------------------------------------#
done
#------------------------------------------------------------------------------------------#


delfun=16
while [ ${delfun} -gt 1 ]
do
   let delfun=${delfun}-1
   echo "  - Files will be deleted in ${delfun} seconds..."
   sleep 1
done


#------------------------------------------------------------------------------------------#
#     Loop over all polygons.                                                              #
#------------------------------------------------------------------------------------------#
ff=0
while [ ${ff} -lt ${npolys} ]
do
   let ff=${ff}+1
   let line=${ff}+3
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



   if [ "x${1}" == "x-d" ]
   then
      find "${here}/${polyname}" -print -delete
      find "${there}/${polyname}" -print -delete
   else
      /bin/cp "${here}/Template/purge.sh" "${here}/${polyname}/purge.sh"
      /bin/cp "${here}/Template/purge.sh" "${there}/${polyname}/purge.sh"
      cd "${here}/${polyname}"
      ./purge.sh
      cd "${there}/${polyname}"
      ./purge.sh
   fi
done
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#      Replace the executable.                                                             #
#------------------------------------------------------------------------------------------#
cd ${here}
if [ -s ${here}/executable/${execname} ]
then
   rm -frv ${here}/executable/${execname}
   cp -fv ${execsrc}/${execname} ${here}/executable
fi
#------------------------------------------------------------------------------------------#
