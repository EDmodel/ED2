#!/bin/sh
here=`pwd`
lonlat=${here}'/joborder.txt'
desc=`basename ${here}`

#----- Determine the number of polygons to run. -------------------------------------------#
let npolys=`wc -l ${lonlat} | awk '{print $1 }'`-3

echo 'Are you sure that you want to stop all jobs? [y/N]'
read proceed

if [ ${proceed} != 'y' -a ${proceed} != 'Y' ]
then
   exit
fi

echo '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
echo '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
echo ' '
echo '     Look, this will really stop ALL your jobs... Are you sure? [y/N]'
echo ' '
echo '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
echo '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
read proceed

if [ ${proceed} != 'y' -a ${proceed} != 'Y' ]
then
   exit
fi

echo 'Okay then, but if you regret later do not say that I did not warn you...'
echo 'I am giving you a few seconds to kill this script in case you change your mind...'
delfun=11
while [ ${delfun} -gt 1 ]
do
   delfun=`expr ${delfun} - 1`
   echo '  - Job stopping will begin in '${delfun}' seconds...'
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
   iage=`echo ${oi}         | awk '{print $16}'`
   polyisoil=`echo ${oi}    | awk '{print $17}'`
   polyntext=`echo ${oi}    | awk '{print $18}'`
   polysand=`echo ${oi}     | awk '{print $19}'`
   polyclay=`echo ${oi}     | awk '{print $20}'`
   polydepth=`echo ${oi}    | awk '{print $21}'`
   polysoilbc=`echo ${oi}   | awk '{print $22}'`
   polysldrain=`echo ${oi}  | awk '{print $23}'`
   polycol=`echo ${oi}      | awk '{print $24}'`
   slzres=`echo ${oi}       | awk '{print $25}'`
   queue=`echo ${oi}        | awk '{print $26}'`
   metdriver=`echo ${oi}    | awk '{print $27}'`
   dtlsm=`echo ${oi}        | awk '{print $28}'`
   vmfactc3=`echo ${oi}     | awk '{print $29}'`
   vmfactc4=`echo ${oi}     | awk '{print $30}'`
   mphototrc3=`echo ${oi}   | awk '{print $31}'`
   mphototec3=`echo ${oi}   | awk '{print $32}'`
   mphotoc4=`echo ${oi}     | awk '{print $33}'`
   bphotoblc3=`echo ${oi}   | awk '{print $34}'`
   bphotonlc3=`echo ${oi}   | awk '{print $35}'`
   bphotoc4=`echo ${oi}     | awk '{print $36}'`
   kwgrass=`echo ${oi}      | awk '{print $37}'`
   kwtree=`echo ${oi}       | awk '{print $38}'`
   gammac3=`echo ${oi}      | awk '{print $39}'`
   gammac4=`echo ${oi}      | awk '{print $40}'`
   d0grass=`echo ${oi}      | awk '{print $41}'`
   d0tree=`echo ${oi}       | awk '{print $42}'`
   alphac3=`echo ${oi}      | awk '{print $43}'`
   alphac4=`echo ${oi}      | awk '{print $44}'`
   klowco2=`echo ${oi}      | awk '{print $45}'`
   decomp=`echo ${oi}       | awk '{print $46}'`
   rrffact=`echo ${oi}      | awk '{print $47}'`
   growthresp=`echo ${oi}   | awk '{print $48}'`
   lwidthgrass=`echo ${oi}  | awk '{print $49}'`
   lwidthbltree=`echo ${oi} | awk '{print $50}'`
   lwidthnltree=`echo ${oi} | awk '{print $51}'`
   q10c3=`echo ${oi}        | awk '{print $52}'`
   q10c4=`echo ${oi}        | awk '{print $53}'`
   h2olimit=`echo ${oi}     | awk '{print $54}'`
   imortscheme=`echo ${oi}  | awk '{print $55}'`
   ddmortconst=`echo ${oi}  | awk '{print $56}'`
   isfclyrm=`echo ${oi}     | awk '{print $57}'`
   icanturb=`echo ${oi}     | awk '{print $58}'`
   ubmin=`echo ${oi}        | awk '{print $59}'`
   ugbmin=`echo ${oi}       | awk '{print $60}'`
   ustmin=`echo ${oi}       | awk '{print $61}'`
   gamm=`echo ${oi}         | awk '{print $62}'`
   gamh=`echo ${oi}         | awk '{print $63}'`
   tprandtl=`echo ${oi}     | awk '{print $64}'`
   ribmax=`echo ${oi}       | awk '{print $65}'`
   atmco2=`echo ${oi}       | awk '{print $66}'`
   thcrit=`echo ${oi}       | awk '{print $67}'`
   smfire=`echo ${oi}       | awk '{print $68}'`
   ifire=`echo ${oi}        | awk '{print $69}'`
   fireparm=`echo ${oi}     | awk '{print $70}'`
   ipercol=`echo ${oi}      | awk '{print $71}'`
   runoff=`echo ${oi}       | awk '{print $72}'`
   imetrad=`echo ${oi}      | awk '{print $73}'`
   ibranch=`echo ${oi}      | awk '{print $74}'`
   icanrad=`echo ${oi}      | awk '{print $75}'`
   crown=`echo   ${oi}      | awk '{print $76}'`
   ltransvis=`echo ${oi}    | awk '{print $77}'`
   lreflectvis=`echo ${oi}  | awk '{print $78}'`
   ltransnir=`echo ${oi}    | awk '{print $79}'`
   lreflectnir=`echo ${oi}  | awk '{print $80}'`
   orienttree=`echo ${oi}   | awk '{print $81}'`
   orientgrass=`echo ${oi}  | awk '{print $82}'`
   clumptree=`echo ${oi}    | awk '{print $83}'`
   clumpgrass=`echo ${oi}   | awk '{print $84}'`
   ivegtdyn=`echo ${oi}     | awk '{print $85}'`
   igndvap=`echo ${oi}      | awk '{print $86}'`
   iphen=`echo ${oi}        | awk '{print $87}'`
   iallom=`echo ${oi}       | awk '{print $88}'`
   ibigleaf=`echo ${oi}     | awk '{print $89}'`
   irepro=`echo ${oi}       | awk '{print $90}'`
   treefall=`echo ${oi}     | awk '{print $91}'`
   ianthdisturb=`echo ${oi} | awk '{print $92}'`
   ianthdataset=`echo ${oi} | awk '{print $93}'`
   #---------------------------------------------------------------------------------------#

   bkill -J ${desc}-${polyname} -q ${queue}
done
#------------------------------------------------------------------------------------------#



