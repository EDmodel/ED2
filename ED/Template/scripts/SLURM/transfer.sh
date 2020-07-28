#!/bin/bash

#==========================================================================================#
#==========================================================================================#
#    Main settings:                                                                        #
#------------------------------------------------------------------------------------------#
#----- Main path, usually set by $(pwd) so you don't need to change it. -------------------#
here=""
#----- Location to create the copy of these files. ----------------------------------------#
there=""
#----- File containing the list of jobs and their settings: -------------------------------#
joborder="${here}/joborder.txt"
#----- Command to be used for rsync. ------------------------------------------------------#
frsync="rsync -Putq  --links --copy-unsafe-links"
rrsync="rsync -Prutq --links --copy-unsafe-links"
#------------------------------------------------------------------------------------------#





#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#     Unless you are going to modify the scripts, you don't need to change anything beyond #
# this point.                                                                              #
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#


#------------------------------------------------------------------------------------------#
#       First check that the main path and e-mail have been set.  If not, don't run.       #
#------------------------------------------------------------------------------------------#
if [ "x${here}" == "x" ] || [ "x${there}" == "x" ]
then
   echo " You must set some variables before running the script:"
   echo " Check variables \"here\" and \"there\"!"
   exit 99
fi
#------------------------------------------------------------------------------------------#




#----- Check whether run_sitter.sh is still running or not.  If it is, exit. --------------#
if [[ "${here}" == "${there}" ]]
then
   echo " Source and destination paths are the same. Transfer is not needed."
   exit 0
elif [[ -s ${here}/transfer.lock ]]
then
   echo " Script transfer is running. Skip transfer for the time being."
   exit
else
   echo "I am going to back up your run." > ${here}/transfer.lock
fi
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#   Create the output path in case it isn't there.                                         #
#------------------------------------------------------------------------------------------#
if [ ! -s ${there} ]
then
   echo "Create backup path: ${there}."
   mkdir -p ${there}
fi
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     First, copy the files at the basal directory.                                        #
#------------------------------------------------------------------------------------------#
echo " + Copy files from the main directory."
${frsync} ${here}/* ${there}
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Then copy the fixed directories.                                                     #
#------------------------------------------------------------------------------------------#
echo " + Copy executable."
${rrsync} ${here}/executable ${there}
echo " + Copy sit_utils."
${rrsync} ${here}/sit_utils  ${there}
echo " + Copy Template."
${rrsync} ${here}/Template   ${there}
#------------------------------------------------------------------------------------------#




#----- Determine the number of polygons to run. -------------------------------------------#
let npolys=$(wc -l ${joborder} | awk '{print $1 }')-3
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
   #    Format count.                                                                      #
   #---------------------------------------------------------------------------------------#
   if   [ ${npolys} -ge 10   ] && [ ${npolys} -lt 100   ]
   then
      ffout=$(printf '%2.2i' ${ff})
   elif [ ${npolys} -ge 100  ] && [ ${npolys} -lt 1000  ]
   then
      ffout=$(printf '%3.3i' ${ff})
   elif [ ${npolys} -ge 100  ] && [ ${npolys} -lt 10000 ]
   then
      ffout=$(printf '%4.4i' ${ff})
   else
      ffout=${ff}
   fi
   ffout="${ffout} of ${npolys}"
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      Read the ffth line of the polygon list.  There must be smarter ways of doing     #
   # this, but this works.  Here we obtain the polygon name, and its longitude and         #
   # latitude.                                                                             #
   #---------------------------------------------------------------------------------------#
   oi=$(head -${line} ${joborder} | tail -1)
   polyname=$(echo ${oi}     | awk '{print $1  }')
   polyiata=$(echo ${oi}     | awk '{print $2  }')
   polylon=$(echo ${oi}      | awk '{print $3  }')
   polylat=$(echo ${oi}      | awk '{print $4  }')
   yeara=$(echo ${oi}        | awk '{print $5  }')
   montha=$(echo ${oi}       | awk '{print $6  }')
   datea=$(echo ${oi}        | awk '{print $7  }')
   timea=$(echo ${oi}        | awk '{print $8  }')
   yearz=$(echo ${oi}        | awk '{print $9  }')
   monthz=$(echo ${oi}       | awk '{print $10 }')
   datez=$(echo ${oi}        | awk '{print $11 }')
   timez=$(echo ${oi}        | awk '{print $12 }')
   initmode=$(echo ${oi}     | awk '{print $13 }')
   iscenario=$(echo ${oi}    | awk '{print $14 }')
   isizepft=$(echo ${oi}     | awk '{print $15 }')
   iage=$(echo ${oi}         | awk '{print $16 }')
   imaxcohort=$(echo ${oi}   | awk '{print $17 }')
   polyisoil=$(echo ${oi}    | awk '{print $18 }')
   polyntext=$(echo ${oi}    | awk '{print $19 }')
   polysand=$(echo ${oi}     | awk '{print $20 }')
   polyclay=$(echo ${oi}     | awk '{print $21 }')
   polydepth=$(echo ${oi}    | awk '{print $22 }')
   polysoilbc=$(echo ${oi}   | awk '{print $23 }')
   polysldrain=$(echo ${oi}  | awk '{print $24 }')
   polycol=$(echo ${oi}      | awk '{print $25 }')
   slzres=$(echo ${oi}       | awk '{print $26 }')
   queue=$(echo ${oi}        | awk '{print $27 }')
   metdriver=$(echo ${oi}    | awk '{print $28 }')
   dtlsm=$(echo ${oi}        | awk '{print $29 }')
   monyrstep=$(echo ${oi}    | awk '{print $30 }')
   iphysiol=$(echo ${oi}     | awk '{print $31 }')
   vmfactc3=$(echo ${oi}     | awk '{print $32 }')
   vmfactc4=$(echo ${oi}     | awk '{print $33 }')
   mphototrc3=$(echo ${oi}   | awk '{print $34 }')
   mphototec3=$(echo ${oi}   | awk '{print $35 }')
   mphotoc4=$(echo ${oi}     | awk '{print $36 }')
   bphotoblc3=$(echo ${oi}   | awk '{print $37 }')
   bphotonlc3=$(echo ${oi}   | awk '{print $38 }')
   bphotoc4=$(echo ${oi}     | awk '{print $39 }')
   kwgrass=$(echo ${oi}      | awk '{print $40 }')
   kwtree=$(echo ${oi}       | awk '{print $41 }')
   gammac3=$(echo ${oi}      | awk '{print $42 }')
   gammac4=$(echo ${oi}      | awk '{print $43 }')
   d0grass=$(echo ${oi}      | awk '{print $44 }')
   d0tree=$(echo ${oi}       | awk '{print $45 }')
   alphac3=$(echo ${oi}      | awk '{print $46 }')
   alphac4=$(echo ${oi}      | awk '{print $47 }')
   klowco2=$(echo ${oi}      | awk '{print $48 }')
   decomp=$(echo ${oi}       | awk '{print $49 }')
   rrffact=$(echo ${oi}      | awk '{print $50 }')
   growthresp=$(echo ${oi}   | awk '{print $51 }')
   lwidthgrass=$(echo ${oi}  | awk '{print $52 }')
   lwidthbltree=$(echo ${oi} | awk '{print $53 }')
   lwidthnltree=$(echo ${oi} | awk '{print $54 }')
   q10c3=$(echo ${oi}        | awk '{print $55 }')
   q10c4=$(echo ${oi}        | awk '{print $56 }')
   h2olimit=$(echo ${oi}     | awk '{print $57 }')
   imortscheme=$(echo ${oi}  | awk '{print $58 }')
   ddmortconst=$(echo ${oi}  | awk '{print $59 }')
   cbrscheme=$(echo ${oi}    | awk '{print $60 }')
   isfclyrm=$(echo ${oi}     | awk '{print $61 }')
   icanturb=$(echo ${oi}     | awk '{print $62 }')
   ubmin=$(echo ${oi}        | awk '{print $63 }')
   ugbmin=$(echo ${oi}       | awk '{print $64 }')
   ustmin=$(echo ${oi}       | awk '{print $65 }')
   gamm=$(echo ${oi}         | awk '{print $66 }')
   gamh=$(echo ${oi}         | awk '{print $67 }')
   tprandtl=$(echo ${oi}     | awk '{print $68 }')
   ribmax=$(echo ${oi}       | awk '{print $69 }')
   atmco2=$(echo ${oi}       | awk '{print $70 }')
   thcrit=$(echo ${oi}       | awk '{print $71 }')
   smfire=$(echo ${oi}       | awk '{print $72 }')
   ifire=$(echo ${oi}        | awk '{print $73 }')
   fireparm=$(echo ${oi}     | awk '{print $74 }')
   ipercol=$(echo ${oi}      | awk '{print $75 }')
   runoff=$(echo ${oi}       | awk '{print $76 }')
   imetrad=$(echo ${oi}      | awk '{print $77 }')
   ibranch=$(echo ${oi}      | awk '{print $78 }')
   icanrad=$(echo ${oi}      | awk '{print $79 }')
   ihrzrad=$(echo ${oi}      | awk '{print $80 }')
   crown=$(echo   ${oi}      | awk '{print $81 }')
   ltransvis=$(echo ${oi}    | awk '{print $82 }')
   lreflectvis=$(echo ${oi}  | awk '{print $83 }')
   ltransnir=$(echo ${oi}    | awk '{print $84 }')
   lreflectnir=$(echo ${oi}  | awk '{print $85 }')
   orienttree=$(echo ${oi}   | awk '{print $86 }')
   orientgrass=$(echo ${oi}  | awk '{print $87 }')
   clumptree=$(echo ${oi}    | awk '{print $88 }')
   clumpgrass=$(echo ${oi}   | awk '{print $89 }')
   igoutput=$(echo ${oi}     | awk '{print $90 }')
   ivegtdyn=$(echo ${oi}     | awk '{print $91 }')
   ihydro=$(echo ${oi}       | awk '{print $92 }')
   istemresp=$(echo ${oi}    | awk '{print $93 }')
   istomata=$(echo ${oi}     | awk '{print $94 }')
   iplastic=$(echo ${oi}     | awk '{print $95 }')
   icarbonmort=$(echo ${oi}  | awk '{print $96 }')
   ihydromort=$(echo ${oi}   | awk '{print $97 }')
   igndvap=$(echo ${oi}      | awk '{print $98 }')
   iphen=$(echo ${oi}        | awk '{print $99 }')
   iallom=$(echo ${oi}       | awk '{print $100}')
   ieconomics=$(echo ${oi}   | awk '{print $101}')
   igrass=$(echo ${oi}       | awk '{print $102}')
   ibigleaf=$(echo ${oi}     | awk '{print $103}')
   integscheme=$(echo ${oi}  | awk '{print $104}')
   nsubeuler=$(echo ${oi}    | awk '{print $105}')
   irepro=$(echo ${oi}       | awk '{print $106}')
   treefall=$(echo ${oi}     | awk '{print $107}')
   ianthdisturb=$(echo ${oi} | awk '{print $108}')
   ianthdataset=$(echo ${oi} | awk '{print $109}')
   slscale=$(echo ${oi}      | awk '{print $110}')
   slyrfirst=$(echo ${oi}    | awk '{print $111}')
   slnyrs=$(echo ${oi}       | awk '{print $112}')
   bioharv=$(echo ${oi}      | awk '{print $113}')
   skidarea=$(echo ${oi}     | awk '{print $114}')
   skidsmall=$(echo ${oi}    | awk '{print $115}')
   skidlarge=$(echo ${oi}    | awk '{print $116}')
   fellingsmall=$(echo ${oi} | awk '{print $117}')
   #---------------------------------------------------------------------------------------#


   #----- Check whether the directories exist or not, and stop the script if they do. -----#
   if [ -s ${here}/${polyname} ]
   then
      echo -n " + Copy ${polyname} (${ffout})."
      
      #----- Sync this directory. ---------------------------------------------------------#
      ${rrsync} ${here}/${polyname} ${there}
      #------------------------------------------------------------------------------------#

      echo "Done!"
   else
      echo " + Skip ${polyname} (${ffout})."
   fi
   #---------------------------------------------------------------------------------------#
done
#------------------------------------------------------------------------------------------#


#----- Clean-up stuff. --------------------------------------------------------------------#
echo " Unlock transfer.sh."
/bin/rm -f ${here}/transfer.lock
echo "==== transfer.sh execution ends. ===="
#------------------------------------------------------------------------------------------#
