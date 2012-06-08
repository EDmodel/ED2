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
   polyisoil=`echo ${oi}    | awk '{print $13}'`
   polyntext=`echo ${oi}    | awk '{print $14}'`
   polysand=`echo ${oi}     | awk '{print $15}'`
   polyclay=`echo ${oi}     | awk '{print $16}'`
   polydepth=`echo ${oi}    | awk '{print $17}'`
   polycol=`echo ${oi}      | awk '{print $18}'`
   slzres=`echo ${oi}       | awk '{print $19}'`
   queue=`echo ${oi}        | awk '{print $20}'`
   metdriver=`echo ${oi}    | awk '{print $21}'`
   dtlsm=`echo ${oi}        | awk '{print $22}'`
   vmfactc3=`echo ${oi}     | awk '{print $23}'`
   vmfactc4=`echo ${oi}     | awk '{print $24}'`
   mphototrc3=`echo ${oi}   | awk '{print $25}'`
   mphototec3=`echo ${oi}   | awk '{print $26}'`
   mphotoc4=`echo ${oi}     | awk '{print $27}'`
   bphotoblc3=`echo ${oi}   | awk '{print $28}'`
   bphotonlc3=`echo ${oi}   | awk '{print $29}'`
   bphotoc4=`echo ${oi}     | awk '{print $30}'`
   kwgrass=`echo ${oi}      | awk '{print $31}'`
   kwtree=`echo ${oi}       | awk '{print $32}'`
   gammac3=`echo ${oi}      | awk '{print $33}'`
   gammac4=`echo ${oi}      | awk '{print $34}'`
   d0grass=`echo ${oi}      | awk '{print $35}'`
   d0tree=`echo ${oi}       | awk '{print $36}'`
   alphac3=`echo ${oi}      | awk '{print $37}'`
   alphac4=`echo ${oi}      | awk '{print $38}'`
   klowco2=`echo ${oi}      | awk '{print $39}'`
   rrffact=`echo ${oi}      | awk '{print $40}'`
   growthresp=`echo ${oi}   | awk '{print $41}'`
   lwidthgrass=`echo ${oi}  | awk '{print $42}'`
   lwidthbltree=`echo ${oi} | awk '{print $43}'`
   lwidthnltree=`echo ${oi} | awk '{print $44}'`
   q10c3=`echo ${oi}        | awk '{print $45}'`
   q10c4=`echo ${oi}        | awk '{print $46}'`
   h2olimit=`echo ${oi}     | awk '{print $47}'`
   ddmort=`echo ${oi}       | awk '{print $48}'`
   isfclyrm=`echo ${oi}     | awk '{print $49}'`
   icanturb=`echo ${oi}     | awk '{print $50}'`
   ubmin=`echo ${oi}        | awk '{print $51}'`
   ugbmin=`echo ${oi}       | awk '{print $52}'`
   ustmin=`echo ${oi}       | awk '{print $53}'`
   gamm=`echo ${oi}         | awk '{print $54}'`
   gamh=`echo ${oi}         | awk '{print $55}'`
   tprandtl=`echo ${oi}     | awk '{print $56}'`
   ribmax=`echo ${oi}       | awk '{print $57}'`
   atmco2=`echo ${oi}       | awk '{print $58}'`
   thcrit=`echo ${oi}       | awk '{print $59}'`
   smfire=`echo ${oi}       | awk '{print $60}'`
   ifire=`echo ${oi}        | awk '{print $61}'`
   fireparm=`echo ${oi}     | awk '{print $62}'`
   ipercol=`echo ${oi}      | awk '{print $63}'`
   isoilbc=`echo ${oi}      | awk '{print $64}'`
   runoff=`echo ${oi}       | awk '{print $65}'`
   imetrad=`echo ${oi}      | awk '{print $66}'`
   ibranch=`echo ${oi}      | awk '{print $67}'`
   icanrad=`echo ${oi}      | awk '{print $68}'`
   crown=`echo   ${oi}      | awk '{print $69}'`
   ltransvis=`echo ${oi}    | awk '{print $70}'`
   lreflectvis=`echo ${oi}  | awk '{print $71}'`
   ltransnir=`echo ${oi}    | awk '{print $72}'`
   lreflectnir=`echo ${oi}  | awk '{print $73}'`
   orienttree=`echo ${oi}   | awk '{print $74}'`
   orientgrass=`echo ${oi}  | awk '{print $75}'`
   clumptree=`echo ${oi}    | awk '{print $76}'`
   clumpgrass=`echo ${oi}   | awk '{print $77}'`
   ivegtdyn=`echo ${oi}     | awk '{print $78}'`
   igndvap=`echo ${oi}      | awk '{print $79}'`
   iphen=`echo ${oi}        | awk '{print $80}'`
   iallom=`echo ${oi}       | awk '{print $81}'`
   ibigleaf=`echo ${oi}     | awk '{print $82}'`
   irepro=`echo ${oi}       | awk '{print $83}'`
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

