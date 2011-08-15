#!/bin/sh
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

   #----- This is just to make sure every node is going to be out of phase. ---------------#
   sleep 1
   #---------------------------------------------------------------------------------------#



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


   #---------------------------------------------------------------------------------------#
   #    Get the correct status of this simulation.                                         #
   #---------------------------------------------------------------------------------------#
   if [ -s ${here}/${polyname}/statusrun.txt ]
   then
      year=`cat ${here}/${polyname}/statusrun.txt  | awk '{print $2}'`
      month=`cat ${here}/${polyname}/statusrun.txt | awk '{print $3}'`
      date=`cat ${here}/${polyname}/statusrun.txt  | awk '{print $4}'`
      time=`cat ${here}/${polyname}/statusrun.txt  | awk '{print $5}'`
      runt=`cat ${here}/${polyname}/statusrun.txt  | awk '{print $6}'`
   else
      year=${yeara}
      month='01'
      date='01'
      time='0000'
      runt='INITIAL'
   fi
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     We will not even consider the files that have gone extinct.                       #
   #---------------------------------------------------------------------------------------#
   if [ ${runt} == 'INITIAL' -o ${runt} == 'HISTORY' ]
   then
      #----- Check whether I should submit from this path or not. -------------------------#
      if [ ! -s ${here}/${polyname}/skipper.txt ]
      then
         blah='Order: '${ff}' Polygon '${polyname}' - Regular job submitted. '
         ${here}/${polyname}/srun.sh 1> /dev/null 2> /dev/null
      elif [ -s ${here}/${polyname}/unparun.sh ]
      then
         blah='Order: '${ff}' Polygon '${polyname}' - Unrestricted_parallel job submitted. '
         ${here}/${polyname}/unparun.sh 1> /dev/null 2> /dev/null
      elif [ ${queue} == 'unrestricted_parallel' ]
      then
         blah='Order: '${ff}' Polygon '${polyname}' - Job will be submitted with unrestricted_parallel'
      else
         blah='Order: '${ff}' Polygon '${polyname}' - Job will not be submitted this time.'
      fi
   elif [ ${runt} == 'THE_END' ]
   then
      blah='Order: '${ff}' Polygon '${polyname}' - This polygon has already finished.'
   elif [ ${runt} == 'STSTATE' ]
   then
      blah='Order: '${ff}' Polygon '${polyname}' - This polygon has already reached steady state.'
   elif [ ${runt} == 'EXTINCT' ]
   then
      blah='Order: '${ff}' Polygon '${polyname}' - This polygon has gone extinct.'
   else
      blah='Order: '${ff}' Polygon '${polyname}' - No idea of what is going on.'
   fi
   #---------------------------------------------------------------------------------------#

   echo ${blah}
done
#------------------------------------------------------------------------------------------#

