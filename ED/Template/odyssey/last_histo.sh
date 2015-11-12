#!/bin/bash
. ${HOME}/.bashrc

#==========================================================================================#
#==========================================================================================#
#     This script keeps only the latest few history files (in case it needs to resume the  #
# simulation), and compresses analysis files.                                              #
#------------------------------------------------------------------------------------------#
here="/x/xxxxxxxxxxxx/xxxxxx/xxxxxxx/xxxxxxxx" # Main path
diskthere=""                                   # Disk where the output files are
joborder="${here}/joborder.txt"                # File with the job instructions
bzip2="/bin/gzip -9"                           # Program to compress files (with options)
checkhourly="y"                                # Check hourly files.
checkstatus="y"                                # Check status before compressing
#------ Calculator. -----------------------------------------------------------------------#
ccc="${HOME}/util/calc.sh"  # Calculator
#------ # of history files to keep (in odyssey's world, always keep more than one). -------#
retain=3
#------------------------------------------------------------------------------------------#




#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#  No need to change anything beyond this point unless you are developing the code.        #
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#



#------------------------------------------------------------------------------------------#
#     Set the maximum number of days for each month.  day[0] is a dummy variable.          #
#------------------------------------------------------------------------------------------#
daymax=(0 31 28 31 30 31 30 31 31 30 31 30 31)
#------------------------------------------------------------------------------------------#



#----- This script is normally ran with crontab.  Make sure the main path is set. ---------#
if [ ${here} == "/x/xxxxxxxxxxxx/xxxxxx/xxxxxxx/xxxxxxxx" ]
then
   echo " here = ${here} "
   echo " Set up variable here before running email_reset.sh"
   exit 92
fi
#------------------------------------------------------------------------------------------#



#-----Make sure last_histo.sh isn't running. ----------------------------------------------#
if [ -s ${here}/last_histo.lock ]
then
   exit
else
   echo "I'm going to clean or compress files. Lots of them!" > ${here}/last_histo.lock
fi
#------------------------------------------------------------------------------------------#



#----- Make yes/no choices case insensitive. ----------------------------------------------#
checkhourly=$(echo ${checkhourly} | awk '{print substr($1,1,1)}' | tr '[:upper:]' '[:lower:]')
checkstatus=$(echo ${checkstatus} | awk '{print substr($1,1,1)}' | tr '[:upper:]' '[:lower:]')
#------------------------------------------------------------------------------------------#


#----- Find the disk here to create the "there" path. -------------------------------------#
moi=$(whoami)
namehere=$(basename ${here})
diskhere=$(dirname ${here})
while [ ${namehere} != ${moi} ]
do
   namehere=$(basename ${diskhere})
   diskhere=$(dirname ${diskhere})
done
if [ "x${diskthere}" == "x" ]
then
   there=${here}
else
   there=$(echo ${here} | sed s@${diskhere}@${diskthere}@g)
fi
#------------------------------------------------------------------------------------------#




#----- Determine the number of polygons to run. -------------------------------------------#
let npolys=$(wc -l ${joborder} | awk '{print $1 }')-3
echo "Number of polygons: ${npolys}..."
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#      Loop over all polygons.                                                             #
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
   #---------------------------------------------------------------------------------------


   #----- Find time and minute. -----------------------------------------------------------#
   houra=$(echo ${timea}  | awk '{print substr($1,1,2)}')
   minua=$(echo ${timea}  | awk '{print substr($1,3,2)}')
   hourz=$(echo ${timez}  | awk '{print substr($1,1,2)}')
   minuz=$(echo ${timez}  | awk '{print substr($1,3,2)}')
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #    If the directory exists, delete all history files but the last.                    #
   #---------------------------------------------------------------------------------------#
   if [ -s ${here}/${polyname} ]
   then
      echo "${ff} - Deleting history files for polygon ${polyname}:"

      #---- First we delete all -Z- files. ------------------------------------------------#
      nzed=$(/bin/ls -1 ${there}/${polyname}/histo/${polyname}-Z-*h5 2> /dev/null | wc -l)
      if [ ${nzed} -gt 0 ]
      then
         zeds=$(/bin/ls -1 ${there}/${polyname}/histo/${polyname}-Z-*h5 2> /dev/null)
         for zed in ${zeds}
         do
            echo -n "    - Deleting: $(basename ${zed})..."
            /bin/nice /bin/rm -f ${zed}
            echo " Gone!"
         done
      fi
      #------------------------------------------------------------------------------------#



      #---- Now we delete all -S- files except the last ${retain} ones. -------------------#
      ness=$(/bin/ls -1 ${there}/${polyname}/histo/${polyname}-S-*h5 2> /dev/null | wc -l)
      if [ ${ness} -gt ${retain} ]
      then
         let head=${ness}-${retain}
         esses=$(/bin/ls -1 ${there}/${polyname}/histo/${polyname}-S-*h5 | head -${head})
         for ess in ${esses}
         do
            echo -n "    - Deleting: $(basename ${ess})..."
            /bin/nice /bin/rm -f ${ess}
            echo " Gone!"
         done
      fi
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #      Decide whether to check the status before compressing files.                  #
      #------------------------------------------------------------------------------------#
      analy="${there}/${polyname}/analy"
      status="${there}/${polyname}/rdata_month/status_${polyname}.txt"

      #----- Check the status before compressing?. ----------------------------------------#
      if [ "x${checkstatus}" == "xy" ]
      then
         #----- Check status of read_monthly.sh. ------------------------------------------#
         if [ -s ${status} ]
         then
            #----- Compress files until the last processed file. --------------------------#
            yeare=$(cat ${status} | awk '{print $1}')
            monthe=$(cat ${status} | awk '{print $2}')
            compress="y"
            #------------------------------------------------------------------------------#
         else
            #----- Do not compress any file. ----------------------------------------------#
            compress="n"
            #------------------------------------------------------------------------------#
         fi
         #---------------------------------------------------------------------------------#
      else
         #----- Ignore checks, just go on and compress everything. ------------------------#
         yeare=${yearz}
         monthe=${monthz}
         compress="y"
         #---------------------------------------------------------------------------------#
      fi
      #------------------------------------------------------------------------------------#




      #------------------------------------------------------------------------------------#
      #     If compress is "y", then go on and compress files.                             #
      #------------------------------------------------------------------------------------#
      if [ "x${compress}" == "xy" ]
      then
         let year=${yeara}-1


         #---------------------------------------------------------------------------------#
         #     Loop over the years.                                                        #
         #---------------------------------------------------------------------------------#
         while [ ${year} -lt ${yeare} ]
         do
            #----- Update year and year string. -------------------------------------------#
            let year=${year}+1
            yyyy=$(printf "%4.4i" ${year})
            #------------------------------------------------------------------------------#

            #----- Update daymax for February (check for leap year). ----------------------#
            let leap400=${year}%400
            let leap100=${year}%100
            let leap004=${year}%4
            if [ ${leap400} -eq 0 ] || [ ${leap100} -ne 0 -a ${leap004} -eq 0 ]
            then
               daymax[2]=29
            else
               daymax[2]=28
            fi
            #------------------------------------------------------------------------------#




            #------ Check last month. -----------------------------------------------------#
            if [ ${year} -eq ${yeare} ]
            then
               monthl=${monthe}
            else
               monthl=12
            fi
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #     Loop over the months.                                                    #
            #------------------------------------------------------------------------------#
            month=0
            while [ ${month} -lt ${monthl} ]
            do
               #----- Update month and month string. --------------------------------------#
               let month=${month}+1
               mm=$(printf "%2.2i" ${month})
               datel=${daymax[${month}]}
               #---------------------------------------------------------------------------#


               #----- Compress monthly mean file. -----------------------------------------#
               qfile=${analy}/${polyname}-Q-${yyyy}-${mm}-00-000000-g01.h5
               if [ -s ${qfile} ]
               then
                  echo -n "   - Compress file: $(basename ${qfile})..."
                  ${bzip2} ${qfile} 2> /dev/null
                  echo "Zipped!"
               fi
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #      Check hourly files.                                                  #
               #---------------------------------------------------------------------------#
               if [ "x${checkhourly}" == "xy" ] || [ "x${checkhourly}" == "xY" ]
               then

                  echo "   - Check hourly files: ${mm}/${yyyy}"


                  #------------------------------------------------------------------------#
                  #     Loop over days.                                                    #
                  #------------------------------------------------------------------------#
                  date=0
                  while [ ${date} -lt ${datel} ]
                  do
                     let date=${date}+1
                     dd=$(printf "%2.2i" ${date})
                     hourl=23
                     hour=-1
                     #---------------------------------------------------------------------#
                     #     Loop over hours.                                                #
                     #---------------------------------------------------------------------#
                     while [ ${hour} -lt ${hourl} ]
                     do
                        let hour=${hour}+1
                        hh=$(printf "%2.2i" ${hour})

                        ifile=${analy}/${polyname}-I-${yyyy}-${mm}-${dd}-${hh}0000-g01.h5
                        if [ -s ${ifile} ]
                        then
                           echo -n "     * Compress file: $(basename ${ifile})..."
                           ${bzip2} ${ifile} 2> /dev/null
                           echo "Zipped!"
                        fi
                        #------------------------------------------------------------------#
                     done
                     #---------------------------------------------------------------------#
                  done
                  #------------------------------------------------------------------------#
               fi
               #---------------------------------------------------------------------------#
            done
            #------------------------------------------------------------------------------#
         done
         #---------------------------------------------------------------------------------#
      fi
      #------------------------------------------------------------------------------------#
   else
      echo "${ff} - Directory doesn't exist"
   fi
   #---------------------------------------------------------------------------------------#
done
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#      Keep this as the very last command.  Release the directory, so last_histo.sh can    #
# be called again.                                                                         #
#------------------------------------------------------------------------------------------#
/bin/rm -f ${here}/last_histo.lock
#------------------------------------------------------------------------------------------#
