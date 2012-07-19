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
   polysoilbc=`echo ${oi}   | awk '{print $18}'`
   polysldrain=`echo ${oi}  | awk '{print $19}'`
   polycol=`echo ${oi}      | awk '{print $20}'`
   slzres=`echo ${oi}       | awk '{print $21}'`
   queue=`echo ${oi}        | awk '{print $22}'`
   metdriver=`echo ${oi}    | awk '{print $23}'`
   dtlsm=`echo ${oi}        | awk '{print $24}'`
   vmfactc3=`echo ${oi}     | awk '{print $25}'`
   vmfactc4=`echo ${oi}     | awk '{print $26}'`
   mphototrc3=`echo ${oi}   | awk '{print $27}'`
   mphototec3=`echo ${oi}   | awk '{print $28}'`
   mphotoc4=`echo ${oi}     | awk '{print $29}'`
   bphotoblc3=`echo ${oi}   | awk '{print $30}'`
   bphotonlc3=`echo ${oi}   | awk '{print $31}'`
   bphotoc4=`echo ${oi}     | awk '{print $32}'`
   kwgrass=`echo ${oi}      | awk '{print $33}'`
   kwtree=`echo ${oi}       | awk '{print $34}'`
   gammac3=`echo ${oi}      | awk '{print $35}'`
   gammac4=`echo ${oi}      | awk '{print $36}'`
   d0grass=`echo ${oi}      | awk '{print $37}'`
   d0tree=`echo ${oi}       | awk '{print $38}'`
   alphac3=`echo ${oi}      | awk '{print $39}'`
   alphac4=`echo ${oi}      | awk '{print $40}'`
   klowco2=`echo ${oi}      | awk '{print $41}'`
   rrffact=`echo ${oi}      | awk '{print $42}'`
   growthresp=`echo ${oi}   | awk '{print $43}'`
   lwidthgrass=`echo ${oi}  | awk '{print $44}'`
   lwidthbltree=`echo ${oi} | awk '{print $45}'`
   lwidthnltree=`echo ${oi} | awk '{print $46}'`
   q10c3=`echo ${oi}        | awk '{print $47}'`
   q10c4=`echo ${oi}        | awk '{print $48}'`
   h2olimit=`echo ${oi}     | awk '{print $49}'`
   imortscheme=`echo ${oi}  | awk '{print $50}'`
   ddmortconst=`echo ${oi}  | awk '{print $51}'`
   isfclyrm=`echo ${oi}     | awk '{print $52}'`
   icanturb=`echo ${oi}     | awk '{print $53}'`
   ubmin=`echo ${oi}        | awk '{print $54}'`
   ugbmin=`echo ${oi}       | awk '{print $55}'`
   ustmin=`echo ${oi}       | awk '{print $56}'`
   gamm=`echo ${oi}         | awk '{print $57}'`
   gamh=`echo ${oi}         | awk '{print $58}'`
   tprandtl=`echo ${oi}     | awk '{print $59}'`
   ribmax=`echo ${oi}       | awk '{print $60}'`
   atmco2=`echo ${oi}       | awk '{print $61}'`
   thcrit=`echo ${oi}       | awk '{print $62}'`
   smfire=`echo ${oi}       | awk '{print $63}'`
   ifire=`echo ${oi}        | awk '{print $64}'`
   fireparm=`echo ${oi}     | awk '{print $65}'`
   ipercol=`echo ${oi}      | awk '{print $66}'`
   runoff=`echo ${oi}       | awk '{print $67}'`
   imetrad=`echo ${oi}      | awk '{print $68}'`
   ibranch=`echo ${oi}      | awk '{print $69}'`
   icanrad=`echo ${oi}      | awk '{print $70}'`
   crown=`echo   ${oi}      | awk '{print $71}'`
   ltransvis=`echo ${oi}    | awk '{print $72}'`
   lreflectvis=`echo ${oi}  | awk '{print $73}'`
   ltransnir=`echo ${oi}    | awk '{print $74}'`
   lreflectnir=`echo ${oi}  | awk '{print $75}'`
   orienttree=`echo ${oi}   | awk '{print $76}'`
   orientgrass=`echo ${oi}  | awk '{print $77}'`
   clumptree=`echo ${oi}    | awk '{print $78}'`
   clumpgrass=`echo ${oi}   | awk '{print $79}'`
   ivegtdyn=`echo ${oi}     | awk '{print $80}'`
   igndvap=`echo ${oi}      | awk '{print $81}'`
   iphen=`echo ${oi}        | awk '{print $82}'`
   iallom=`echo ${oi}       | awk '{print $83}'`
   ibigleaf=`echo ${oi}     | awk '{print $84}'`
   irepro=`echo ${oi}       | awk '{print $85}'`
   treefall=`echo ${oi}     | awk '{print $86}'`
   #---------------------------------------------------------------------------------------#


   ncname=${polyname}'.out'
   ncerror=${polyname}'_err.out'
    
   if [ -s ${here}/${polyname}/${ncname}  ]
   then
      if [ -s ${here}/${polyname}/${ncerror} ]
      then
         ncDONE=`grep "Successfully completed" ${here}/${polyname}/${ncerror} | wc -l`
         if [ ${ncDONE} -gt 0 ]
         then
            echo ':-D ncfile for '${polyname}' has finished ...'
         else
            echo ':-( looks like ncfile '${polyname}' has crashed ... <============='
         fi
      else
         WAITING=`grep "Waiting for file..." ${here}/${polyname}/${ncname} | tail -n1`
         wcheck=`grep "Waiting for file..." ${here}/${polyname}/${ncname} | tail -n1 | wc -l`
         simulating=`grep "file /" ${here}/${polyname}/${ncname} | tail -n1`
         scheck=`grep "file /" ${here}/${polyname}/${ncname}  | tail -n1 | wc -l`
         if [ ${wcheck} -gt 0 ]
         then 
            if [ ${scheck} -gt 0 ]
            then
                echo ':-) '${polyname}' has loaded '${simulating:(-47):23}' and is waiting. '${WAITING: -13}'  ...'
            else
                echo ':-| '${polyname}' has waited '${WAITING: -6}'  ...'
            fi
         elif [ ${scheck} -gt 0 ]
         then
            echo ':-) '${polyname}' is running. '${simulating: -47}'...'
         else
            echo ':-| '${polyname}' is pending ...'
         fi
      fi
   fi

done

