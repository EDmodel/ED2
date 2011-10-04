#!/bin/sh

if [ 'x'${1} == 'x' ]
then
   ncol=1
else
   ncol=${1}
fi
echo ${ncol}

here=`pwd`
lonlat=${here}'/joborder.txt'

#----- Determine the number of polygons to run. -------------------------------------------#
let npolys=`wc -l ${lonlat} | awk '{print $1 }'`-3
echo 'Number of polygons: '${npolys}'...'


#------------------------------------------------------------------------------------------#
#     Loop over all polygons.                                                              #
#------------------------------------------------------------------------------------------#
ff=0
while [ ${ff} -lt ${npolys} ]
do
   let ff=${ff}+1
   let line=${ff}+3

   #----- Make two columns. ---------------------------------------------------------------#
   let col=${ff}%${ncol}
   if [ ${ncol} -eq 1 ]
   then
      opt=''
      off=''
   elif [ ${col} -eq 0 ]
   then
      opt=''
      off='.\t'
   else
      opt='-n'
      off=''
   fi

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
   polycol=`echo ${oi}     | awk '{print $18}'`
   slzres=`echo ${oi}      | awk '{print $19}'`
   queue=`echo ${oi}       | awk '{print $20}'`
   metdriver=`echo ${oi}   | awk '{print $21}'`
   dtlsm=`echo ${oi}       | awk '{print $22}'`
   vmfactc3=`echo ${oi}    | awk '{print $23}'`
   vmfactc4=`echo ${oi}    | awk '{print $24}'`
   mphotoc3=`echo ${oi}    | awk '{print $25}'`
   mphotoc4=`echo ${oi}    | awk '{print $26}'`
   kwgrass=`echo ${oi}     | awk '{print $27}'`
   kwtree=`echo ${oi}      | awk '{print $28}'`
   gammac3=`echo ${oi}     | awk '{print $29}'`
   gammac4=`echo ${oi}     | awk '{print $30}'`
   d0grass=`echo ${oi}     | awk '{print $31}'`
   d0tree=`echo ${oi}      | awk '{print $32}'`
   d0decay=`echo ${oi}     | awk '{print $33}'`
   alphac3=`echo ${oi}     | awk '{print $34}'`
   alphac4=`echo ${oi}     | awk '{print $35}'`
   klowco2=`echo ${oi}     | awk '{print $36}'`
   rrffact=`echo ${oi}     | awk '{print $37}'`
   growthresp=`echo ${oi}  | awk '{print $38}'`
   h2olimit=`echo ${oi}    | awk '{print $39}'`
   isfclyrm=`echo ${oi}    | awk '{print $40}'`
   icanturb=`echo ${oi}    | awk '{print $41}'`
   ubmin=`echo ${oi}       | awk '{print $42}'`
   ugbmin=`echo ${oi}      | awk '{print $43}'`
   ustmin=`echo ${oi}      | awk '{print $44}'`
   gamm=`echo ${oi}        | awk '{print $45}'`
   gamh=`echo ${oi}        | awk '{print $46}'`
   tprandtl=`echo ${oi}    | awk '{print $47}'`
   ribmax=`echo ${oi}      | awk '{print $48}'`
   atmco2=`echo ${oi}      | awk '{print $49}'`
   thcrit=`echo ${oi}      | awk '{print $50}'`
   smfire=`echo ${oi}      | awk '{print $51}'`
   isoilbc=`echo ${oi}     | awk '{print $52}'`
   imetrad=`echo ${oi}     | awk '{print $53}'`
   ibranch=`echo ${oi}     | awk '{print $54}'`
   icanrad=`echo ${oi}     | awk '{print $55}'`
   crown=`echo   ${oi}     | awk '{print $56}'`
   ltransvis=`echo ${oi}   | awk '{print $57}'`
   lreflectvis=`echo ${oi} | awk '{print $58}'`
   ltransnir=`echo ${oi}   | awk '{print $59}'`
   lreflectnir=`echo ${oi} | awk '{print $60}'`
   orienttree=`echo ${oi}  | awk '{print $61}'`
   orientgrass=`echo ${oi} | awk '{print $62}'`
   clumptree=`echo ${oi}   | awk '{print $63}'`
   clumpgrass=`echo ${oi}  | awk '{print $64}'`
   ivegtdyn=`echo ${oi}    | awk '{print $65}'`
   igndvap=`echo ${oi}     | awk '{print $66}'`
   iphen=`echo ${oi}       | awk '{print $67}'`
   #---------------------------------------------------------------------------------------#



   if [ -s ${here}/${polyname}/serial_out.out ]
   then
      fatal=`grep "FATAL ERROR" ${here}/${polyname}/serial_out.out | wc -l`
      simulating=`grep "Simulating: " ${here}/${polyname}/serial_out.out | tail -1`
      if [ ${fatal} -gt 0 ]
      then 
         echo -e ${opt} ${off}':-( '${polyname}' HAS CRASHED  ... <======================'
      elif [ -s ${here}/${polyname}/serial_lsf.out ]
      then
         echo -e ${opt} ${off}':-D '${polyname}' has finished  ...'
      else
         echo -e ${opt} ${off}':-) '${polyname}' is running. '${simulating}'...'
      fi
   else
      echo -e ${opt} ${off}':-| '${polyname}' is pending ...'
   fi
done

