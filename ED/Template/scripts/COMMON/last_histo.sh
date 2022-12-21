#!/bin/bash
. ~/.bashrc

#==========================================================================================#
#==========================================================================================#
#     This script keeps only the latest few history files (in case it needs to resume the  #
# simulation), and compresses analysis files.                                              #
#------------------------------------------------------------------------------------------#
here=""                                        # Main path
diskthere=""                                   # Disk where the output files are
joborder="${here}/joborder.txt"                # File with the job instructions
bzip2="/bin/gzip -9"                           # Program to compress files (with options)
checkhourly=""                                 # Check hourly files.
checkstatus=""                                 # Check status before compressing
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



#------------------------------------------------------------------------------------------#
#     Look for variables that may nThis script is normally ran with crontab.  Make sure the main path is set. ---------#
#------------------------------------------------------------------------------------------#
bad=0
if [[ "x${here}"        == "x" ]] || [[ "x${checkhourly}" == "x" ]] || 
   [[ "x${checkstatus}" == "x" ]]
then
   echo "------------------------------------------------------------"
   echo " here        = ${here}"
   echo " checkhourly = ${checkhourly}"
   echo " checkstatus = ${checkstatus}"
   echo " None of the variables above should be empty."
   echo " Please check last_histo.sh settings."
   echo "------------------------------------------------------------"
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
if [[ "x${diskthere}" == "x" ]] || [[ "x${diskthere}" == "x${here}" ]]
then
   there=${here}
else
   there=$(echo ${here} | sed s@${diskhere}@${diskthere}@g)
fi
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Find the extention of the compressed file.                                           #
#------------------------------------------------------------------------------------------#
case ${bzip2} in
*gzip*) 
   ezip="gz"
   ;;
*bzip2*)
   ezip="bz2"
   ;;
*compress*)
   ezip="Z"
   ;;
*)
   echo " Compressing tool (${bzip2}) is not recognised."
   exit 1
   ;;
esac
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
   polyname=$(echo ${oi}      | awk '{print $1  }')
   polyiata=$(echo ${oi}      | awk '{print $2  }')
   polylon=$(echo ${oi}       | awk '{print $3  }')
   polylat=$(echo ${oi}       | awk '{print $4  }')
   yeara=$(echo ${oi}         | awk '{print $5  }')
   montha=$(echo ${oi}        | awk '{print $6  }')
   datea=$(echo ${oi}         | awk '{print $7  }')
   timea=$(echo ${oi}         | awk '{print $8  }')
   yearz=$(echo ${oi}         | awk '{print $9  }')
   monthz=$(echo ${oi}        | awk '{print $10 }')
   datez=$(echo ${oi}         | awk '{print $11 }')
   timez=$(echo ${oi}         | awk '{print $12 }')
   initmode=$(echo ${oi}      | awk '{print $13 }')
   iscenario=$(echo ${oi}     | awk '{print $14 }')
   isizepft=$(echo ${oi}      | awk '{print $15 }')
   iage=$(echo ${oi}          | awk '{print $16 }')
   imaxcohort=$(echo ${oi}    | awk '{print $17 }')
   polyisoil=$(echo ${oi}     | awk '{print $18 }')
   polyntext=$(echo ${oi}     | awk '{print $19 }')
   polysand=$(echo ${oi}      | awk '{print $20 }')
   polyclay=$(echo ${oi}      | awk '{print $21 }')
   polyslsoc=$(echo ${oi}     | awk '{print $22 }')
   polyslph=$(echo ${oi}      | awk '{print $23 }')
   polyslcec=$(echo ${oi}     | awk '{print $24 }')
   polysldbd=$(echo ${oi}     | awk '{print $25 }')
   polydepth=$(echo ${oi}     | awk '{print $26 }')
   polyslhydro=$(echo ${oi}   | awk '{print $27 }')
   polysoilbc=$(echo ${oi}    | awk '{print $28 }')
   polysldrain=$(echo ${oi}   | awk '{print $29 }')
   polycol=$(echo ${oi}       | awk '{print $30 }')
   slzres=$(echo ${oi}        | awk '{print $31 }')
   queue=$(echo ${oi}         | awk '{print $32 }')
   metdriver=$(echo ${oi}     | awk '{print $33 }')
   dtlsm=$(echo ${oi}         | awk '{print $34 }')
   monyrstep=$(echo ${oi}     | awk '{print $35 }')
   iphysiol=$(echo ${oi}      | awk '{print $36 }')
   vmfactc3=$(echo ${oi}      | awk '{print $37 }')
   vmfactc4=$(echo ${oi}      | awk '{print $38 }')
   mphototrc3=$(echo ${oi}    | awk '{print $39 }')
   mphototec3=$(echo ${oi}    | awk '{print $40 }')
   mphotoc4=$(echo ${oi}      | awk '{print $41 }')
   bphotoblc3=$(echo ${oi}    | awk '{print $42 }')
   bphotonlc3=$(echo ${oi}    | awk '{print $43 }')
   bphotoc4=$(echo ${oi}      | awk '{print $44 }')
   kwgrass=$(echo ${oi}       | awk '{print $45 }')
   kwtree=$(echo ${oi}        | awk '{print $46 }')
   gammac3=$(echo ${oi}       | awk '{print $47 }')
   gammac4=$(echo ${oi}       | awk '{print $48 }')
   d0grass=$(echo ${oi}       | awk '{print $49 }')
   d0tree=$(echo ${oi}        | awk '{print $50 }')
   alphac3=$(echo ${oi}       | awk '{print $51 }')
   alphac4=$(echo ${oi}       | awk '{print $52 }')
   klowco2=$(echo ${oi}       | awk '{print $53 }')
   decomp=$(echo ${oi}        | awk '{print $54 }')
   rrffact=$(echo ${oi}       | awk '{print $55 }')
   growthresp=$(echo ${oi}    | awk '{print $56 }')
   lwidthgrass=$(echo ${oi}   | awk '{print $57 }')
   lwidthbltree=$(echo ${oi}  | awk '{print $58 }')
   lwidthnltree=$(echo ${oi}  | awk '{print $59 }')
   q10c3=$(echo ${oi}         | awk '{print $60 }')
   q10c4=$(echo ${oi}         | awk '{print $61 }')
   h2olimit=$(echo ${oi}      | awk '{print $62 }')
   imortscheme=$(echo ${oi}   | awk '{print $63 }')
   ddmortconst=$(echo ${oi}   | awk '{print $64 }')
   cbrscheme=$(echo ${oi}     | awk '{print $65 }')
   isfclyrm=$(echo ${oi}      | awk '{print $66 }')
   icanturb=$(echo ${oi}      | awk '{print $67 }')
   ubmin=$(echo ${oi}         | awk '{print $68 }')
   ugbmin=$(echo ${oi}        | awk '{print $69 }')
   ustmin=$(echo ${oi}        | awk '{print $70 }')
   gamm=$(echo ${oi}          | awk '{print $71 }')
   gamh=$(echo ${oi}          | awk '{print $72 }')
   tprandtl=$(echo ${oi}      | awk '{print $73 }')
   ribmax=$(echo ${oi}        | awk '{print $74 }')
   atmco2=$(echo ${oi}        | awk '{print $75 }')
   thcrit=$(echo ${oi}        | awk '{print $76 }')
   smfire=$(echo ${oi}        | awk '{print $77 }')
   ifire=$(echo ${oi}         | awk '{print $78 }')
   fireparm=$(echo ${oi}      | awk '{print $79 }')
   ipercol=$(echo ${oi}       | awk '{print $80 }')
   runoff=$(echo ${oi}        | awk '{print $81 }')
   imetrad=$(echo ${oi}       | awk '{print $82 }')
   ibranch=$(echo ${oi}       | awk '{print $83 }')
   icanrad=$(echo ${oi}       | awk '{print $84 }')
   ihrzrad=$(echo ${oi}       | awk '{print $85 }')
   crown=$(echo   ${oi}       | awk '{print $86 }')
   ltransvis=$(echo ${oi}     | awk '{print $87 }')
   lreflectvis=$(echo ${oi}   | awk '{print $88 }')
   ltransnir=$(echo ${oi}     | awk '{print $89 }')
   lreflectnir=$(echo ${oi}   | awk '{print $90 }')
   orienttree=$(echo ${oi}    | awk '{print $91 }')
   orientgrass=$(echo ${oi}   | awk '{print $92 }')
   clumptree=$(echo ${oi}     | awk '{print $93 }')
   clumpgrass=$(echo ${oi}    | awk '{print $94 }')
   igoutput=$(echo ${oi}      | awk '{print $95 }')
   ivegtdyn=$(echo ${oi}      | awk '{print $96 }')
   ihydro=$(echo ${oi}        | awk '{print $97 }')
   istemresp=$(echo ${oi}     | awk '{print $98 }')
   istomata=$(echo ${oi}      | awk '{print $99 }')
   iplastic=$(echo ${oi}      | awk '{print $100}')
   icarbonmort=$(echo ${oi}   | awk '{print $101}')
   ihydromort=$(echo ${oi}    | awk '{print $102}')
   igndvap=$(echo ${oi}       | awk '{print $103}')
   iphen=$(echo ${oi}         | awk '{print $104}')
   iallom=$(echo ${oi}        | awk '{print $105}')
   ieconomics=$(echo ${oi}    | awk '{print $106}')
   igrass=$(echo ${oi}        | awk '{print $107}')
   ibigleaf=$(echo ${oi}      | awk '{print $108}')
   integscheme=$(echo ${oi}   | awk '{print $109}')
   nsubeuler=$(echo ${oi}     | awk '{print $110}')
   irepro=$(echo ${oi}        | awk '{print $111}')
   treefall=$(echo ${oi}      | awk '{print $112}')
   ianthdisturb=$(echo ${oi}  | awk '{print $113}')
   ianthdataset=$(echo ${oi}  | awk '{print $114}')
   slscale=$(echo ${oi}       | awk '{print $115}')
   slyrfirst=$(echo ${oi}     | awk '{print $116}')
   slnyrs=$(echo ${oi}        | awk '{print $117}')
   bioharv=$(echo ${oi}       | awk '{print $118}')
   skidarea=$(echo ${oi}      | awk '{print $119}')
   skiddbhthresh=$(echo ${oi} | awk '{print $120}')
   skidsmall=$(echo ${oi}     | awk '{print $121}')
   skidlarge=$(echo ${oi}     | awk '{print $122}')
   fellingsmall=$(echo ${oi}  | awk '{print $123}')
   #---------------------------------------------------------------------------------------


   #----- Find time and minute. -----------------------------------------------------------#
   houra=$(echo ${timea}  | awk '{print substr($1,1,2)}')
   minua=$(echo ${timea}  | awk '{print substr($1,3,2)}')
   hourz=$(echo ${timez}  | awk '{print substr($1,1,2)}')
   minuz=$(echo ${timez}  | awk '{print substr($1,3,2)}')
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Set gfilout prefix according to ihrzrad.                                          #
   #---------------------------------------------------------------------------------------#
   case ${ihrzrad} in
   1) 
      gpref="gap"
      ;;
   2)
      gpref="pix"
      ;;
   *)
      gpref="dum"
      ;;
   esac
   #---------------------------------------------------------------------------------------#


   #----- Define the main path for this polygon. ------------------------------------------#
   pypath="${there}/${polyname}"
   sfilout="${pypath}/histo/${polyname}"
   ffilout="${pypath}/analy/${polyname}"
   gfilout="${pypath}/shade/${gpref}"
   #---------------------------------------------------------------------------------------#

   #---------------------------------------------------------------------------------------#
   #    If the directory exists, delete all history files but the last.                    #
   #---------------------------------------------------------------------------------------#
   if [ -s ${pypath} ]
   then
      echo "${ff} - Delete history files for polygon ${polyname}:"

      #---- First we delete all -Z- files. ------------------------------------------------#
      nzed=$(/bin/ls -1 ${sfilout}-Z-*h5 2> /dev/null | wc -l)
      if [ ${nzed} -gt 0 ]
      then
         zeds=$(/bin/ls -1 ${sfilout}-Z-*h5 2> /dev/null)
         for zed in ${zeds}
         do
            echo -n "    - Delete: $(basename ${zed})..."
            /bin/nice /bin/rm -f ${zed}
            echo " Gone!"
         done
      fi
      #------------------------------------------------------------------------------------#



      #---- Now we delete all -S- files except the last ${retain} ones. -------------------#
      ness=$(/bin/ls -1 ${sfilout}-S-*h5 2> /dev/null | wc -l)
      if [ ${ness} -gt ${retain} ]
      then
         let head=${ness}-${retain}
         esses=$(/bin/ls -1 ${sfilout}-S-*h5 | head -${head})
         for ess in ${esses}
         do
            echo -n "    - Delete: $(basename ${ess})..."
            /bin/nice /bin/rm -f ${ess}
            echo " Gone!"
         done
      fi
      #------------------------------------------------------------------------------------#



      #---- Now we compress all raster files except for the last one. ---------------------#
      nrst=$(/bin/ls -1 ${gfilout}_raster_isi001_????-01.txt 2> /dev/null | wc -l)
      if [ ${nrst} -gt 1 ]
      then
         let head=${nrst}-1
         rasters=$(/bin/ls -1 ${gfilout}_raster_isi001_????-01.txt | head -${head})
         for rst in ${rasters}
         do
            zipped="${rst}.${ezip}"
            if [[ -s ${rst} ]] && [[ -f ${zipped} ]]
            then
               echo " - Remove previously compressed file ($(basename ${zipped}))."
               /bin/rm -f ${zipped}
            fi
            echo -n "    - Compress: $(basename ${rst})..."
            ${bzip2} ${rst}
            echo " Compressed!"
         done
      fi
      #------------------------------------------------------------------------------------#



      #---- Now we compress all ptable files except for the last one. ---------------------#
      nptb=$(/bin/ls -1 ${gfilout}_ptable_isi001_????-01.txt 2> /dev/null | wc -l)
      if [ ${nptb} -gt 1 ]
      then
         let head=${nptb}-1
         ptables=$(/bin/ls -1 ${gfilout}_ptable_isi001_????-01.txt | head -${head})
         for ptb in ${ptables}
         do
            zipped="${ptb}.${ezip}"
            if [[ -s ${ptb} ]] && [[ -f ${zipped} ]]
            then
               echo " - Remove previously compressed file ($(basename ${zipped}))."
               /bin/rm -f ${zipped}
            fi
            echo -n "    - Compress: $(basename ${ptb})..."
            ${bzip2} ${ptb}
            echo " Compressed!"
         done
      fi
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #      Decide whether to check the status before compressing files.                  #
      #------------------------------------------------------------------------------------#
      status="${pypath}/rdata_month/status_${polyname}.txt"

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
               qfile=${ffilout}-Q-${yyyy}-${mm}-00-000000-g01.h5
               if [ -s ${qfile} ]
               then
                  zipped="${qfile}.${ezip}"
                  if [[ -f ${zipped} ]]
                  then
                     echo " - Remove previously compressed file ($(basename ${zipped}))."
                     /bin/rm -f ${zipped}
                  fi
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

                        ifile=${ffilout}-I-${yyyy}-${mm}-${dd}-${hh}0000-g01.h5
                        if [ -s ${ifile} ]
                        then
                           zipped="${ifile}.${ezip}"
                           basezip=$(basename ${zipped})
                           if [[ -f ${zipped} ]]
                           then
                              echo " - Remove previously compressed file (${basezip})."
                              /bin/rm -f ${zipped}
                           fi
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
