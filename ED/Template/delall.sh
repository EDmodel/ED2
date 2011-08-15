#!/bin/sh
mainpath=`pwd`
lonlat=${mainpath}'/joborder.txt'
let npolys=`wc -l ${lonlat} | awk '{print $1 }'`-3

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
   polyname=`echo ${oi}    | awk '{print $1 }'`
   polyiata=`echo ${oi}    | awk '{print $2 }'`
   polylon=`echo ${oi}     | awk '{print $3 }'`
   polylat=`echo ${oi}     | awk '{print $4 }'`
   yeara=`echo ${oi}       | awk '{print $5 }'`
   montha=`echo ${oi}      | awk '{print $6 }'`
   datea=`echo ${oi}       | awk '{print $7 }'`
   timea=`echo ${oi}       | awk '{print $8 }'`
   yearz=`echo ${oi}       | awk '{print $9 }'`
   monthz=`echo ${oi}      | awk '{print $10}'`
   datez=`echo ${oi}       | awk '{print $11}'`
   timez=`echo ${oi}       | awk '{print $12}'`
   polyisoil=`echo ${oi}   | awk '{print $13}'`
   polyntext=`echo ${oi}   | awk '{print $14}'`
   polysand=`echo ${oi}    | awk '{print $15}'`
   polyclay=`echo ${oi}    | awk '{print $16}'`
   polydepth=`echo ${oi}   | awk '{print $17}'`
   slzres=`echo ${oi}      | awk '{print $18}'`
   queue=`echo ${oi}       | awk '{print $19}'`
   metdriver=`echo ${oi}   | awk '{print $20}'`
   dtlsm=`echo ${oi}       | awk '{print $21}'`
   vmfact=`echo ${oi}      | awk '{print $22}'`
   mfact=`echo ${oi}       | awk '{print $23}'`
   kfact=`echo ${oi}       | awk '{print $24}'`
   gamfact=`echo ${oi}     | awk '{print $25}'`
   d0fact=`echo ${oi}      | awk '{print $26}'`
   alphafact=`echo ${oi}   | awk '{print $27}'`
   rrffact=`echo ${oi}     | awk '{print $28}'`
   growthresp=`echo ${oi}  | awk '{print $29}'`
   isfclyrm=`echo ${oi}    | awk '{print $30}'`
   icanturb=`echo ${oi}    | awk '{print $31}'`
   ubmin=`echo ${oi}       | awk '{print $32}'`
   ugbmin=`echo ${oi}      | awk '{print $33}'`
   ustmin=`echo ${oi}      | awk '{print $34}'`
   gamm=`echo ${oi}        | awk '{print $35}'`
   gamh=`echo ${oi}        | awk '{print $36}'`
   tprandtl=`echo ${oi}    | awk '{print $37}'`
   ribmax=`echo ${oi}      | awk '{print $38}'`
   atmco2=`echo ${oi}      | awk '{print $39}'`
   thcrit=`echo ${oi}      | awk '{print $40}'`
   smfire=`echo ${oi}      | awk '{print $41}'`
   isoilbc=`echo ${oi}     | awk '{print $42}'`
   imetrad=`echo ${oi}     | awk '{print $43}'`
   ibranch=`echo ${oi}     | awk '{print $44}'`
   icanrad=`echo ${oi}     | awk '{print $45}'`
   ltransvis=`echo ${oi}   | awk '{print $46}'`
   lreflectvis=`echo ${oi} | awk '{print $47}'`
   ltransnir=`echo ${oi}   | awk '{print $48}'`
   lreflectnir=`echo ${oi} | awk '{print $49}'`
   orienttree=`echo ${oi}  | awk '{print $50}'`
   orientgrass=`echo ${oi} | awk '{print $51}'`
   clumptree=`echo ${oi}   | awk '{print $52}'`
   clumpgrass=`echo ${oi}  | awk '{print $53}'`
   ivegtdyn=`echo ${oi}    | awk '{print $54}'`
   igndvap=`echo ${oi}     | awk '{print $55}'`
   #---------------------------------------------------------------------------------------#



   therepath=`echo ${mainpath} | sed s@/n/Moorcroft_Lab/Users@/n/moorcroft_scratch@g`
   if [ 'x'${1} == 'x-d' ]
   then
      rm -frv ${mainpath}'/'${polyname} 
      rm -frv ${therepath}'/'${polyname} 
   else
      cd ${mainpath}'/'${polyname}
      ./purge.sh
   fi
done
#------------------------------------------------------------------------------------------#
