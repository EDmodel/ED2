#!/bin/sh
here=`pwd`
moi=`whoami`
diskthere='/n/moorcroftfs2'
lonlat=${here}'/joborder.txt'

#----- Find the output path (both local and remote paths will be cleaned). ----------------#
basehere=`basename ${here}`
dirhere=`dirname ${here}`
while [ ${basehere} != ${moi} ]
do
   basehere=`basename ${dirhere}`
   dirhere=`dirname ${dirhere}`
done
diskhere=${dirhere}
echo '-------------------------------------------------------------------------------'
echo ' - Simulation control on disk: '${diskhere}
echo ' - Output on disk:             '${diskthere}
echo '-------------------------------------------------------------------------------'
there=`echo ${here} | sed s@${diskhere}@${diskthere}@g`
#------------------------------------------------------------------------------------------#




#----- Determine the number of polygons to run. -------------------------------------------#
let npolys=`wc -l ${lonlat} | awk '{print $1 }'`-3
#------------------------------------------------------------------------------------------#




#----- Check that the user is aware that it will remove everything... ---------------------#
if [ 'x'${1} == 'x-d' ]
then
   echo 'Are you sure that you want to remove all files and directories? [y/N]'
else
   echo 'Are you sure that you want to remove all files? [y/N]'
fi
read proceed
if [ ${proceed} != 'y' -a ${proceed} != 'Y' ]
then
   exit
fi
#------------------------------------------------------------------------------------------#



#----- Check that the user is aware that it will remove everything... ---------------------#
echo ' '
if [ 'x'${1} == 'x-d' ]
then
   echo '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   echo '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   echo ' '
   echo '     Look, this will REALLY delete all '${npolys}' output directories and files...'
   echo ' '
   echo '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   echo '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
else
   echo '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   echo '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   echo ' '
   echo '     Look, this will REALLY delete all '${npolys}' output files...'
echo ' '
   echo '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   echo '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
fi

echo 'This is PERMANENT, once they are gone, adieu, no chance to recover them!'
echo 'Is that what you really want? [y/N]'
read proceed

echo ' '

if [ ${proceed} != 'y' -a ${proceed} != 'Y' ]
then 
   exit
fi

echo 'Okay then, but if you regret later do not say that I did not warn you...'
echo 'I am giving you a few seconds to kill this script in case you change your mind...'
delfun=16
while [ ${delfun} -gt 1 ]
do
   delfun=`expr ${delfun} - 1`
   echo '  - Deletion will begin in '${delfun}' seconds...'
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
   oi=`head -${line} ${lonlat} | tail -1`
   polyname=`echo ${oi}     | awk '{print $1 }'`
   polyiata=`echo ${oi}     | awk '{print $2 }'`
   polylon=`echo ${oi}      | awk '{print $3 }'`
   polylat=`echo ${oi}      | awk '{print $4 }'`
   yeara=`echo ${oi}        | awk '{print $5 }'`
   montha=`echo ${oi}       | awk '{print $6 }'`
   datea=`echo ${oi}        | awk '{print $7 }'`
   timea=`echo ${oi}        | awk '{print $8 }'`
   yearz=`echo ${oi}        | awk '{print $9 }'`
   monthz=`echo ${oi}       | awk '{print $10}'`
   datez=`echo ${oi}        | awk '{print $11}'`
   timez=`echo ${oi}        | awk '{print $12}'`
   initmode=`echo ${oi}     | awk '{print $13}'`
   iscenario=`echo ${oi}    | awk '{print $14}'`
   isizepft=`echo ${oi}     | awk '{print $15}'`
   polyisoil=`echo ${oi}    | awk '{print $16}'`
   polyntext=`echo ${oi}    | awk '{print $17}'`
   polysand=`echo ${oi}     | awk '{print $18}'`
   polyclay=`echo ${oi}     | awk '{print $19}'`
   polydepth=`echo ${oi}    | awk '{print $20}'`
   polysoilbc=`echo ${oi}   | awk '{print $21}'`
   polysldrain=`echo ${oi}  | awk '{print $22}'`
   polycol=`echo ${oi}      | awk '{print $23}'`
   slzres=`echo ${oi}       | awk '{print $24}'`
   queue=`echo ${oi}        | awk '{print $25}'`
   metdriver=`echo ${oi}    | awk '{print $26}'`
   dtlsm=`echo ${oi}        | awk '{print $27}'`
   vmfactc3=`echo ${oi}     | awk '{print $28}'`
   vmfactc4=`echo ${oi}     | awk '{print $29}'`
   mphototrc3=`echo ${oi}   | awk '{print $30}'`
   mphototec3=`echo ${oi}   | awk '{print $31}'`
   mphotoc4=`echo ${oi}     | awk '{print $32}'`
   bphotoblc3=`echo ${oi}   | awk '{print $33}'`
   bphotonlc3=`echo ${oi}   | awk '{print $34}'`
   bphotoc4=`echo ${oi}     | awk '{print $35}'`
   kwgrass=`echo ${oi}      | awk '{print $36}'`
   kwtree=`echo ${oi}       | awk '{print $37}'`
   gammac3=`echo ${oi}      | awk '{print $38}'`
   gammac4=`echo ${oi}      | awk '{print $39}'`
   d0grass=`echo ${oi}      | awk '{print $40}'`
   d0tree=`echo ${oi}       | awk '{print $41}'`
   alphac3=`echo ${oi}      | awk '{print $42}'`
   alphac4=`echo ${oi}      | awk '{print $43}'`
   klowco2=`echo ${oi}      | awk '{print $44}'`
   decomp=`echo ${oi}       | awk '{print $45}'`
   rrffact=`echo ${oi}      | awk '{print $46}'`
   growthresp=`echo ${oi}   | awk '{print $47}'`
   lwidthgrass=`echo ${oi}  | awk '{print $48}'`
   lwidthbltree=`echo ${oi} | awk '{print $49}'`
   lwidthnltree=`echo ${oi} | awk '{print $50}'`
   q10c3=`echo ${oi}        | awk '{print $51}'`
   q10c4=`echo ${oi}        | awk '{print $52}'`
   h2olimit=`echo ${oi}     | awk '{print $53}'`
   imortscheme=`echo ${oi}  | awk '{print $54}'`
   ddmortconst=`echo ${oi}  | awk '{print $55}'`
   isfclyrm=`echo ${oi}     | awk '{print $56}'`
   icanturb=`echo ${oi}     | awk '{print $57}'`
   ubmin=`echo ${oi}        | awk '{print $58}'`
   ugbmin=`echo ${oi}       | awk '{print $59}'`
   ustmin=`echo ${oi}       | awk '{print $60}'`
   gamm=`echo ${oi}         | awk '{print $61}'`
   gamh=`echo ${oi}         | awk '{print $62}'`
   tprandtl=`echo ${oi}     | awk '{print $63}'`
   ribmax=`echo ${oi}       | awk '{print $64}'`
   atmco2=`echo ${oi}       | awk '{print $65}'`
   thcrit=`echo ${oi}       | awk '{print $66}'`
   smfire=`echo ${oi}       | awk '{print $67}'`
   ifire=`echo ${oi}        | awk '{print $68}'`
   fireparm=`echo ${oi}     | awk '{print $69}'`
   ipercol=`echo ${oi}      | awk '{print $70}'`
   runoff=`echo ${oi}       | awk '{print $71}'`
   imetrad=`echo ${oi}      | awk '{print $72}'`
   ibranch=`echo ${oi}      | awk '{print $73}'`
   icanrad=`echo ${oi}      | awk '{print $74}'`
   crown=`echo   ${oi}      | awk '{print $75}'`
   ltransvis=`echo ${oi}    | awk '{print $76}'`
   lreflectvis=`echo ${oi}  | awk '{print $77}'`
   ltransnir=`echo ${oi}    | awk '{print $78}'`
   lreflectnir=`echo ${oi}  | awk '{print $79}'`
   orienttree=`echo ${oi}   | awk '{print $80}'`
   orientgrass=`echo ${oi}  | awk '{print $81}'`
   clumptree=`echo ${oi}    | awk '{print $82}'`
   clumpgrass=`echo ${oi}   | awk '{print $83}'`
   ivegtdyn=`echo ${oi}     | awk '{print $84}'`
   igndvap=`echo ${oi}      | awk '{print $85}'`
   iphen=`echo ${oi}        | awk '{print $86}'`
   iallom=`echo ${oi}       | awk '{print $87}'`
   ibigleaf=`echo ${oi}     | awk '{print $88}'`
   irepro=`echo ${oi}       | awk '{print $89}'`
   treefall=`echo ${oi}     | awk '{print $90}'`
   #---------------------------------------------------------------------------------------#



   if [ 'x'${1} == 'x-d' ]
   then
      rm -frv ${here}'/'${polyname} 
      rm -frv ${there}'/'${polyname} 
   else
      cd ${here}'/'${polyname}
      ./purge.sh
      cd ${there}'/'${polyname}
      ./purge.sh
   fi
done
#------------------------------------------------------------------------------------------#
