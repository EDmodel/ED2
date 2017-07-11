 #!/bin/sh

here=`pwd`
there=${here}
lonlat=${here}'/joborder.txt'

#----- Determine the number of polygons to run. -------------------------------------------#
let npolys=`wc -l ${lonlat} | awk '{print $1 }'`-3
echo 'Number of polygons: '${npolys}'...'

#------------------------------------------------------------------------------------------#
#     Loop over all polygons.                                                              #
#------------------------------------------------------------------------------------------#
ff=0
if [ -s ${here}/nc_files ]
then
    echo 'nc_files already exists, make sure you are not going to overwrite nc files'
    echo 'I am giving you a few seconds to kill this script in case you change your mind...'
    repfun=7
    while [ ${repfun} -gt 1 ]
    do
       let repfun=${repfun}-1
       echo '  - Script will start in '${repfun}' seconds...'
       sleep 1
    done
else
    mkdir ${here}/nc_files
fi

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



   #----- Find time and minute. -----------------------------------------------------------#
   houra=`echo ${timea}  | awk '{print substr($1,1,2)}'`
   minua=`echo ${timea}  | awk '{print substr($1,3,2)}'`
   hourz=`echo ${timez}  | awk '{print substr($1,1,2)}'`
   minuz=`echo ${timez}  | awk '{print substr($1,3,2)}'`
   #---------------------------------------------------------------------------------------#




   echo 'Order: '${ff}', creating gen_nc file for polygon '${polyname}'...'
   #----- Copy the Template directory to a unique polygon directory. ----------------------#
   cp -f ${here}/Template/gen_netcdf.m ${here}/${polyname}/gen_netcdf.m
   cp -f ${here}/Template/mrun.sh ${here}/${polyname}/mrun.sh
   #---------------------------------------------------------------------------------------#


   #----- Use default year for beginning of the time period. ------------------------------#
   thisyeara=${yeara}
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Find the last day with full diurnal cycle.                                        #
   #---------------------------------------------------------------------------------------#
   if [ ${hourz} -lt 23 ]
   then
      #------------------------------------------------------------------------------------#
      #     Last day is not complete, use previous day instead.                            #
      #------------------------------------------------------------------------------------#
      let thisdatez=${datez}-1
      if [ ${thisdatez} -eq 0 ]
      then
         #----- Last time is on the first day of the month, go back one month. ------------#
         let thismonthz=${monthz}-1
         if [ ${thismonthz} -eq 0 ]
         then
            #----- Last time is in January, go back to December on the previous year. -----#
            thismonthz=12
            let thisyearz=${yearz}-1
         else
            #----- Last time is not January, keep the same year. --------------------------#
            thisyearz=${yearz}
         fi

         #----- Find the last date with full diurnal cycle. -------------------------------#
         case thismonth in
         1|3|5|7|8|10|12)
            thisdatez=31
            ;;
         4|6|9|11)
            thisdatez=30
            ;;
         2)
            #------ Previous month is February, find out whether it's a leap year. --------#
            let rem400=${thisyear}%400
            let rem100=${thisyear}%100
            let rem4=${thisyear}%4
            
            if [ ${rem400} -eq 0 ] || [ ${rem4} -eq 0 -a ${rem100} -ne 0 ]
            then
               thisdatez=29
            else
               thisdatez=28
            fi
            ;;
         esac
      else
         thismonthz=${monthz}
         thisyearz=${yearz}
      fi
      #------------------------------------------------------------------------------------#
   else
      #------------------------------------------------------------------------------------#
      #     Last day is complete, use the last day as is.                                  #
      #------------------------------------------------------------------------------------#
      thisdatez=${datez}
      thismonthz=${monthz}
      thisyearz=${yearz}
   fi
   #---------------------------------------------------------------------------------------#



   #----- Pad zeroes to the left of date and month if needed be. --------------------------#
   usethismonthz=`printf "%02d" ${thismonthz}`
   usethisdatez=`printf "%02d" ${thisdatez}`

   #---------------------------------------------------------------------------------------#
   #     Change the gen_netcdf file.                                                       #
   #---------------------------------------------------------------------------------------#
   genNC=${here}'/'${polyname}'/gen_netcdf.m'
   sed -i s@paththere@${there}@g         ${genNC}
   sed -i s@myyeara@${thisyeara}@g       ${genNC}
   sed -i s@mymontha@${montha}@g         ${genNC}
   sed -i s@myyearz@${thisyearz}@g       ${genNC}
   sed -i s@mydatez@${usethisdatez}@g    ${genNC}
   sed -i s@mymonthz@${usethismonthz}@g  ${genNC}
   sed -i s@thispoly@${polyname}@g       ${genNC}
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #    Prepare the script to run the netCDF
   #---------------------------------------------------------------------------------------#
   mrun=${here}'/'${polyname}'/mrun.sh'
   sed -i s@thisqueue@${queue}@g ${mrun}
   #---------------------------------------------------------------------------------------#
done
