#!/bin/sh
here=$(pwd)
moi=$(whoami)
diskthere="/n/moorcroftfs2"
lonlat="${here}/joborder.txt"

desc=$(basename ${here})

#----- Executable name. -------------------------------------------------------------------#
execname="ed_2.1-opt"
execsrc="${HOME}/EDBRAMS/ED/build"
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
let npolys=$(wc -l ${lonlat} | awk '{print $1 }')-3
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

   scancel -n ${desc}-${polyname} -p ${queue}
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



   if [ "x${1}" == "x-d" ]
   then
      rm -frv "${here}/${polyname}"
      rm -frv "${there}/${polyname}"
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
