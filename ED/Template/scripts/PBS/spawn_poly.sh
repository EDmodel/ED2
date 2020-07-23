#!/bin/bash

#==========================================================================================#
#==========================================================================================#
#    Main settings:                                                                        #
#------------------------------------------------------------------------------------------#
#----- Main path, usually set by $(pwd) so you don't need to change it. -------------------#
here=$(pwd)
#----- User name, usually set by $(whoami) so you don't need to change it. ----------------#
moi=$(whoami)
#----- Description of this simulation, used to create unique job names. -------------------#
desc=$(basename ${here})
#----- Select main file system path. ------------------------------------------------------#
ordinateur=$(hostname -s)
case ${ordinateur} in
  sun-master|cmm*)  export fs0="/prj/prjidfca/${moi}"                       ;;
  ha*)              export fs0="/halo_nobackup/jpl-gewc/${moi}"             ;;
  au*)              export fs0="/aurora_nobackup/jpl-gewc/${moi}"           ;;
  *)    echo " Invalid computer ${ordinateur}.  Check script header."; exit ;;
esac
#----- Path where biomass initialisation files are: ---------------------------------------#
bioinit="${fs0}/Data/ed2_data/site_bio_data"
alsinit="${fs0}/Data/ed2_data/lidar_spline_bio_data"
intinit="${fs0}/Data/ed2_data/lidar_intensity_bio_data"
lutinit="${fs0}/Data/ed2_data/lidar_lookup_bio_data"
biotype=0      # 0 -- "default" setting (isizepft controls default/nounder)
               # 1 -- isizepft controls number of PFTs, whereas iage controls patches.
               # 2 -- airborne lidar initialisation using return counts ("default"). 
               # 3 -- airborne lidar initialisation using intensity counts.
               # 4 -- airborne lidar/inventory hybrid initialisation ("lookup table"). 
               # For lidar initialisation (2-4), isizepft is the disturbance history key.
#----- Path and file prefix for init_mode = 5. --------------------------------------------#
restart="${fs0}/Data/ed2_data/restarts_XXX"
#----- File containing the list of jobs and their settings: -------------------------------#
joborder="${here}/joborder.txt"
#----- This is the header with the Sheffield data. ----------------------------------------#
shefhead='SHEF_NCEP_DRIVER_DS314'
#----- Path with drivers for each scenario. -----------------------------------------------#
metmaindef="${fs0}/Data/ed2_data"
packdatasrc="${fs0}/Data/2scratch"
#----- Path with land use scenarios. ------------------------------------------------------#
lumain="${fs0}/Data/lu_scenarios"
#----- Path with other input data bases (soil texture, DGD, land mask, etc). --------------#
inpmain="${fs0}/Data/ed2_data"
#----- If submit is "n", we create paths but skip submission. -----------------------------#
submit="n"
#----- Maximum number of attempts before giving up. ---------------------------------------#
nsubtry_max=5
#------------------------------------------------------------------------------------------#

#------------------------------------------------------------------------------------------#
#      History run variables.                                                              #
#------------------------------------------------------------------------------------------#
#----- Force history run (0 = no, 1 = yes). -----------------------------------------------#
forcehisto=0
#----- Path with the history file to be used. ---------------------------------------------#
fullygrown="${fs0}/Simulations/Debug/D001_Debug/xyz_settings/histo/xyz_settings"
#----- Time that we shall use. ------------------------------------------------------------#
yearh="1510"  # Year
monthh="07"   # Month
dateh="01"    # Day
timeh="0000"  # Hour
#----- Default tolerance. -----------------------------------------------------------------#
toldef="0.01"
#----- Executable names. ------------------------------------------------------------------#
execname="ed_2.2-opt"             # Normal executable, for most queues
#----- Initialisation scripts. ------------------------------------------------------------#
initrc="${HOME}/.bashrc"          # Initialisation script for most nodes
#----- Settings for this group of polygons. -----------------------------------------------#
global_queue="verylongq"      # Queue
n_cpt=8                       # Number of cpus per task (it will be limited by maximum)
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
#==========================================================================================#
#==========================================================================================#
#     Unless you are going to modify the way jobs are submitted, you don't need to change  #
# anything beyond this point.                                                              #
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


#----- Dummy, currently only odyssey needs copying data to scratch. -----------------------#
metmain=${metmaindef}
copy2scratch=false
#------------------------------------------------------------------------------------------#




#----- Determine the number of polygons to run. -------------------------------------------#
let npolys=$(wc -l ${joborder} | awk '{print $1 }')-3
if [ ${npolys} -lt 100 ]
then
   ndig=2
elif [ ${npolys} -lt 1000 ]
then
   ndig=3
elif [ ${npolys} -lt 10000 ]
then
   ndig=4
else
   ndig=5
fi
pfmt="%${ndig}.${ndig}i"
echo "Number of polygons: ${npolys}..."
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#   Check whether the executable is copied.  If not, let the user know and stop the        #
# script.                                                                                  #
#------------------------------------------------------------------------------------------#
exec_full="${here}/executable/${execname}"
if [ ! -s ${exec_full} ]
then
   echo "Executable file : ${exec_full} is not in the executable directory"
   echo "Copy the executable to the file before running this script!"
   exit 99
fi
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Set general limits.  Currently only the maximum time is restricted.                  #
#------------------------------------------------------------------------------------------#
case ${ordinateur} in
sun-master|cmm*)
   #----- SunHPC (LNCC). ------------------------------------------------------------------#
   n_nodes_max=32
   n_cpn=16
   n_cpt_max=8
   node_memory=65536
   case ${global_queue} in
      linuxq)    runtime="168:00:00"                              ;;
      *)         echo " Queue ${global_queue} is invalid!"; exit  ;;
   esac
   #---------------------------------------------------------------------------------------#
   ;;
au*|ha*)
   #----- JPL (Aurora and Halo). ----------------------------------------------------------#
   n_nodes_max=32
   n_cpn=16
   n_cpt_max=8
   node_memory=65536
   case ${global_queue} in
      verylongq) runtime="192:00:00"                              ;;
      longq)     runtime="48:00:00"                               ;;
      mediumq)   runtime="12:00:00"                               ;;
      shortq)    runtime="03:00:00"                               ;;
      debugq)    runtime="01:00:00"                               ;;
      *)         echo " Queue ${global_queue} is invalid!"; exit  ;;
   esac
   #---------------------------------------------------------------------------------------#
   ;;
*)
   #----- Computer is not listed.  Crash. -------------------------------------------------#
   echo " Invalid computer ${ordinateur}.  Check queue settings in the script."
   exit 39
   #---------------------------------------------------------------------------------------#
   ;;
esac
if [ ${n_cpt} -gt ${n_cpt_max} ]
then
   echo " Too many CPUs per task requested:"
   echo " Queue                   = ${global_queue}"
   echo " Maximum CPUs per task   = ${n_cpt_max}"
   echo " Requested CPUs per task = ${n_cpt}"
   exit 99
else
   if [[ ${n_cpt} -eq 0 ]]
   then
      n_cpt=${n_cpt_max}
   fi
   let memory_max=${node_memory}*${n_cpt_max}/${n_cpn}
   let n_tasks_max=${n_nodes_max}*${n_cpn}
   let n_tasks_max=${n_tasks_max}/${n_cpt}
fi
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#   here/there is not an option at sunhpc.                                                 #
#------------------------------------------------------------------------------------------#
basehere=$(basename ${here})
dirhere=$(dirname ${here})
while [ ${basehere} != ${moi} ]
do
   basehere=$(basename ${dirhere})
   dirhere=$(dirname ${dirhere})
done
diskhere=${dirhere}
echo "-------------------------------------------------------------------------------"
echo " - Path for simulation control and output: ${diskhere}"
echo "-------------------------------------------------------------------------------"
there=${here}
#------------------------------------------------------------------------------------------#


#----- Save queue status to a temporary file. ---------------------------------------------#
jobstat="/tmp/spjs.$$"
qjobs > ${jobstat}
#------------------------------------------------------------------------------------------#


#----- Summary for this submission preparation.  Then give 5 seconds for user to cancel. --#
echo "------------------------------------------------"
echo "  Submission summary: "
echo ""
echo "  Memory per cpu:      ${sim_memory}"
echo "  CPUs per node:       ${n_cpn}"
echo "  CPUs per task:       ${n_cpt}"
echo "  Queue:               ${global_queue}"
echo "  Job Name:            ${jobname}"
echo "  Total polygon count: ${npolys}"
echo " "
echo "------------------------------------------------"
echo ""
sleep 2
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#     Loop over all polygons.                                                              #
#------------------------------------------------------------------------------------------#
ff=0
mc2=0
while [ ${ff} -lt ${npolys} ]
do
   let ff=${ff}+1
   let line=${ff}+3
   ffout=$(printf ${pfmt} ${ff})


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



   #---------------------------------------------------------------------------------------#
   #     Decide units for restart files.                                                   #
   #---------------------------------------------------------------------------------------#
   let elapseda=12*${yeara}+${montha}
   let elapsedz=12*${yearz}+${monthz}
   let nmonths=${elapsedz}-${elapseda}
   if [ ${nmonths} -ge 240 ]
   then
      iunitstate=3
   elif [ ${nmonths} -ge 3 ]
   then
      iunitstate=2
   else
      iunitstate=1
   fi
   #---------------------------------------------------------------------------------------#


   #----- Check whether the directories exist or not, and stop the script if they do. -----#
   if [ -s ${here}/${polyname} ]
   then
      echo -n "${ffout} ${polyname}: updating files..."
      
      #----- Save the last tolerance in case we are going to make it more strict. ---------#
      oldtol=$(grep NL%RK4_TOLERANCE ${here}/${polyname}/ED2IN | awk '{print $3}')
      rm -f ${here}/${polyname}/ED2IN 
      rm -f ${here}/${polyname}/callserial.sh
      rm -f ${here}/${polyname}/callunpa.sh 
      rm -f ${here}/${polyname}/skipper.txt
      cp ${here}/Template/ED2IN         ${here}/${polyname}/ED2IN
      cp ${here}/Template/callserial.sh ${here}/${polyname}/callserial.sh
   else
      echo -n "${ffout} ${polyname}: creating directory..."
      #----- Copy the Template directory to a unique polygon directory. -------------------#
      oldtol=""
      cp -r ${here}/Template ${here}/${polyname}
   fi
   #---------------------------------------------------------------------------------------#






   #---------------------------------------------------------------------------------------#
   #   Make sure that we have a reasonable tolerance.                                      #
   #---------------------------------------------------------------------------------------#
   if [ "x${oldtol}" == "xmytoler" -o "x${oldtol}" == "x" ]
   then
      toler=${toldef}
   else
      toler=${oldtol}
   fi
   #---------------------------------------------------------------------------------------#






   #---------------------------------------------------------------------------------------#
   #      Here we check whether the an old directory exists or not.  If it doesn't, then   #
   # we simply copy the template, and assume initial run.  Otherwise, we must find out     #
   # where the simulation was when it stopped.                                             #
   #---------------------------------------------------------------------------------------#
   if [ -s  ${here}/${polyname} ]
   then

      #------------------------------------------------------------------------------------#
      #      This step is necessary because we may have killed the run while it was        #
      # writing, and as a result, the file may be corrupt.                                 #
      #------------------------------------------------------------------------------------#
      nhdf5=$(ls -1 ${here}/${polyname}/histo/* 2> /dev/null | wc -l)
      if [ ${nhdf5} -gt 0 ]
      then
         h5fine=0

         while [ ${h5fine} -eq 0 ]
         do
            lasthdf5=$(ls -1 ${here}/${polyname}/histo/* | tail -1)
            h5dump -H ${lasthdf5} 1> /dev/null 2> ${here}/badfile.txt

            if [ -s ${here}/badfile.txt ]
            then
               /bin/rm -fv ${lasthdf5}
               nhdf5=$(ls -1 ${here}/${polyname}/histo/* 2> /dev/null | wc -l)
               if [ ${nhdf5} -eq 0 ]
               then
                  h5fine=1
               fi
            else
               h5fine=1
            fi

            /bin/rm -f ${here}/badfile.txt
         done
         #---------------------------------------------------------------------------------#
      fi
      #------------------------------------------------------------------------------------#
   fi
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Run the small R script to check whether the simulation was running or not, and   #
   # whether there was any cohort left by the time the runs were stopped.                  #
   #---------------------------------------------------------------------------------------#
   /bin/rm -f ${here}/${polyname}/statusrun.txt
   /bin/rm -f ${here}/${polyname}/whichrun.r
   /bin/cp -f ${here}/Template/whichrun.r ${here}/${polyname}/whichrun.r
   whichrun="${here}/${polyname}/whichrun.r"
   outwhich="${here}/${polyname}/outwhichrun.txt"
   sed -i s@thispoly@${polyname}@g           ${whichrun}
   sed -i s@thisqueue@${queue}@g             ${whichrun}
   sed -i s@pathhere@${here}@g               ${whichrun}
   sed -i s@paththere@${here}@g              ${whichrun}
   sed -i s@thisyeara@${yeara}@g             ${whichrun}
   sed -i s@thismontha@${montha}@g           ${whichrun}
   sed -i s@thisdatea@${datea}@g             ${whichrun}
   sed -i s@thistimea@${timea}@g             ${whichrun}
   sed -i s@thischecksteady@FALSE@g          ${whichrun}
   sed -i s@thismetcyc1@${metcyc1}@g         ${whichrun}
   sed -i s@thismetcycf@${metcycf}@g         ${whichrun}
   sed -i s@thisnyearmin@10000@g             ${whichrun}
   sed -i s@thisststcrit@0.0@g               ${whichrun}
   R CMD BATCH --no-save --no-restore ${whichrun} ${outwhich}
   while [ ! -s ${here}/${polyname}/statusrun.txt ]
   do
      sleep 0.5
   done
   year=$(cat ${here}/${polyname}/statusrun.txt  | awk '{print $2}')
   month=$(cat ${here}/${polyname}/statusrun.txt | awk '{print $3}')
   date=$(cat ${here}/${polyname}/statusrun.txt  | awk '{print $4}')
   time=$(cat ${here}/${polyname}/statusrun.txt  | awk '{print $5}')
   runt=$(cat ${here}/${polyname}/statusrun.txt  | awk '{print $6}')
   #---------------------------------------------------------------------------------------#





   #---------------------------------------------------------------------------------------#
   #     Determine which PFTs to use based on the "iata" code, isizepft, and biotype.      #
   #---------------------------------------------------------------------------------------#
   case ${biotype} in
   2)
      #------------------------------------------------------------------------------------#
      #     Airborne lidar.  isizepft is not controlling isizepft, but disturbance         #
      # history. Initialise PFTs using the default settings.                               #
      #------------------------------------------------------------------------------------#
      case ${polyiata} in
      tzi|zmh|nqn)
         pfts="1,7,8,9,10,11,15,16"
         pasture=16
         crop=16
         plantation=15
         logging="7,8,9,10,11,15"
         probharv="1.0,1.0,1.0,1.0,1.0,1.0"
         dbhharv="0.0,0.0,0.0,0.0,0.0,0.0"
         seedharv=0.75
         storharv=0.00
         leafharv=0.00
         ;;
      hvd|wch|tqh)
         pfts="6,8,9,10,11,16"
         pasture=16
         crop=16
         plantation=11
         logging="6,8,9,10,11"
         probharv="1.0,1.0,1.0,1.0,1.0"
         dbhharv="0.0,0.0,0.0,0.0,0.0"
         seedharv=0.75
         storharv=0.00
         leafharv=0.00
         ;;
      l[0-5][0-3])
         pfts="6,8,9,10,11"
         pasture=16
         crop=16
         plantation=6
         logging="6,8,9,10,11"
         probharv="1.0,1.0,1.0,1.0,1.0"
         dbhharv="0.0,0.0,0.0,0.0,0.0"
         seedharv=0.75
         storharv=0.00
         leafharv=0.00
         ;;
      asu|cnf|bnu|cwb|erm|iqq|ipv|mgf|rao|sla|zpe|kna|sfn)
         pfts="1,2,3,4,15,16"
         pasture=1
         crop=16
         plantation=15
         logging="2,3,4,15"
         probharv="1.0,1.0,1.0,1.0,1.0,1.0"
         dbhharv="50.,50.,50.,50.,50.,50."
         seedharv=0.75
         storharv=0.00
         leafharv=0.00
         ;;
      fns*)
         pfts="1,16"
         pasture=1
         crop=1
         plantation=13
         logging="13"
         probharv="1.0"
         dbhharv="50."
         seedharv=0.75
         storharv=0.00
         leafharv=0.00
         ;;
      s77*)
         pfts="1,16"
         pasture=1
         crop=16
         plantation=13
         logging="13"
         probharv="1.0"
         dbhharv="50."
         seedharv=0.75
         storharv=0.00
         leafharv=0.00
         ;;
      *)
         pfts="1,2,3,4,16"
         pasture=1
         crop=1
         plantation=13
         logging="2,3,4"
         probharv="1.0,1.0,1.0,1.0,1.0"
         dbhharv="50.,50.,50.,50.,50."
         seedharv=0.75
         storharv=0.00
         leafharv=0.00
         ;;
      esac
      ;;
      #------------------------------------------------------------------------------------#
   *)
      case ${polyiata} in
      tzi|zmh|nqn)
         pfts="1,7,8,9,10,11,15,16"
         pasture=16
         crop=16
         plantation=15
         logging="7,8,9,10,11,15"
         probharv="1.0,1.0,1.0,1.0,1.0,1.0"
         dbhharv="0.0,0.0,0.0,0.0,0.0,0.0"
         seedharv=0.75
         storharv=0.00
         leafharv=0.00
         ;;
      hvd|wch|tqh)
         pfts="6,8,9,10,11,16"
         pasture=16
         crop=16
         plantation=10
         logging="6,8,9,10,11"
         probharv="1.0,1.0,1.0,1.0,1.0"
         dbhharv="0.0,0.0,0.0,0.0,0.0"
         seedharv=0.75
         storharv=0.00
         leafharv=0.00
         ;;
      l[0-5][0-3])
         pfts="6,8,9,10,11"
         pasture=16
         crop=16
         plantation=6
         logging="6,8,9,10,11"
         probharv="1.0,1.0,1.0,1.0,1.0"
         dbhharv="0.0,0.0,0.0,0.0,0.0"
         seedharv=0.75
         storharv=0.00
         leafharv=0.00
         ;;
      asu|cnf|bnu|cwb|erm|iqq|ipv|mgf|rao|sla|zpe|kna|sfn)
         pfts="1,2,3,4,15,16"
         pasture=1
         crop=16
         plantation=15
         logging="2,3,4,15"
         probharv="1.0,1.0,1.0,1.0,1.0,1.0"
         dbhharv="50.,50.,50.,50.,50.,50."
         seedharv=0.75
         storharv=0.00
         leafharv=0.00
         ;;
      fns*)
         pfts="1,16"
         pasture=1
         crop=1
         plantation=13
         logging="13"
         probharv="1.0"
         dbhharv="50."
         seedharv=0.75
         storharv=0.00
         leafharv=0.00
         ;;
      s77*)
         pfts="1,16"
         pasture=1
         crop=16
         plantation=13
         logging="13"
         probharv="1.0"
         dbhharv="50."
         seedharv=0.75
         storharv=0.00
         leafharv=0.00
         ;;
      *)
         pfts="1,2,3,4,16"
         pasture=1
         crop=1
         plantation=13
         logging="2,3,4"
         probharv="1.0,1.0,1.0,1.0,1.0"
         dbhharv="50.,50.,50.,50.,50."
         seedharv=0.75
         storharv=0.00
         leafharv=0.00
         ;;
      esac
      ;;
      #------------------------------------------------------------------------------------#
   esac
   #---------------------------------------------------------------------------------------#





   #---------------------------------------------------------------------------------------#
   #     Determine the census structure.                                                   #
   #---------------------------------------------------------------------------------------#
   let yodd=${yeara}%2
   case ${polyiata} in
   gyf)
      dtcensus=24
      let yr1stcensus=${yeara}+${yodd}
      mon1stcensus=3
      minrecruitdbh=10
      ;;
   s67|dcm)
      dtcensus=24
      let yr1stcensus=${yeara}+1-${yodd}
      mon1stcensus=7
      minrecruitdbh=10
      ;;
   *)
      dtcensus=1
      yr1stcensus=${yeara}
      mon1stcensus=1
      minrecruitdbh=10
      ;;
   esac
   #---------------------------------------------------------------------------------------#





   #---------------------------------------------------------------------------------------#
   #     Determine the scenario paths.                                                     #
   #---------------------------------------------------------------------------------------#
   case ${iscenario} in
   default)
      case ${metdriver} in
      ERAINT_CHIRPS)
         #----- ERA-Interim (CHIRPS precipitation). ---------------------------------------#
         scentype="ERA_Interim"
         iscenario="ERAINT_SOUTHAM_CHIRPS"
         ;;
      ERAINT_MSWEP2)
         #----- ERA-Interim (MSWEP2 precipitation). ---------------------------------------#
         scentype="ERA_Interim"
         iscenario="ERAINT_SOUTHAM_MSWEP2"
         ;;
      ERAINT_NATIVE)
         #----- ERA-Interim (native precipitation). ---------------------------------------#
         scentype="ERA_Interim"
         iscenario="ERAINT_SOUTHAM_NATIVE"
         ;;
      MERRA2_CHIRPS)
         #----- MERRA2 (CHIRPS precipitation). --------------------------------------------#
         scentype="MERRA2"
         iscenario="MERRA2_SOUTHAM_CHIRPS"
         ;;
      MERRA2_MSWEP2)
         #----- MERRA2 (MSWEP2 precipitation). --------------------------------------------#
         scentype="MERRA2"
         iscenario="MERRA2_SOUTHAM_MSWEP2"
         ;;
      MERRA2_NATIVE)
         #----- MERRA-2 (native precipitation). -------------------------------------------#
         scentype="MERRA2"
         iscenario="MERRA2_SOUTHAM_NATIVE"
         ;;
      PGMF3_CHIRPS)
         #----- PGMF-3 (CHIRPS precipitation). --------------------------------------------#
         scentype="PGMF3"
         iscenario="PGMF3_SOUTHAM_CHIRPS"
         ;;
      PGMF3_MSWEP2)
         #----- PGMF-3 (CHIRPS precipitation). --------------------------------------------#
         scentype="PGMF3"
         iscenario="PGMF3_SOUTHAM_MSWEP2"
         ;;
      PGMF3_NATIVE)
         #----- PGMF-3 (native precipitation). --------------------------------------------#
         scentype="PGMF3"
         iscenario="PGMF3_SOUTHAM_NATIVE"
         ;;
      Sheffield)
         #----- Sheffield. ----------------------------------------------------------------#
         scentype="sheffield"
         iscenario="sheffield"
         ;;
      WFDEI_CHIRPS)
         #----- WFDEI (CHIRPS Precipitation). ---------------------------------------------#
         scentype="WFDEI"
         iscenario="WFDEI_SOUTHAM_CHIRPS"
         ;;
      WFDEI_CRUP)
         #----- WFDEI (CRU Precipitation). ------------------------------------------------#
         scentype="WFDEI"
         iscenario="WFDEI_SOUTHAM_CRUP"
         ;;
      WFDEI_GPCC)
         #----- WFDEI (GPCC Precipitation). -----------------------------------------------#
         scentype="WFDEI"
         iscenario="WFDEI_SOUTHAM_GPCC"
         ;;
      WFDEI_MSWEP2)
         #----- WFDEI (MSWEP2 Precipitation). ---------------------------------------------#
         scentype="WFDEI"
         iscenario="WFDEI_SOUTHAM_MSWEP2"
         ;;
      *)
         #----- Tower data. ---------------------------------------------------------------#
         scentype="wmo+eft"
         iscenario="eft"
         ;;
      esac
      ;;
   eft|wmo|shr)
      #----- Tower data, keep scenario as is. ---------------------------------------------#
      scentype="wmo+eft"
      ;;
   *)
      #----- Rainfall scenario, keep scenario as is. --------------------------------------#
      scentype="realisation_scen_driver"
      ;;
   esac
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      Choose the scenario to use.  Iscenario follows the following convention:         #
   # "default"    -- No scenario.  Use the tower/Sheffield data.                           #
   # "wmo"        -- No scenario.  Use the WMO-based data.                                 #
   # "rRRR_tTTT   -- Use scenarios, with rRRRR controlling the rainfall, and tTTTT         #
   #                 controlling temperature.                                              #
   #                                                                                       #
   # rRRR         -- Rainfall scenarios, where rRRRR means:                                #
   #                 r+000: INMET-based time series, no resampling.                        #
   #                 r+010: Re-sampling with substitution, but equal chances for all years #
   #                 r-XXX: Shift the location (similar to mean) of the distribution by    #
   #                        -X.XX units of scale (similar to standard deviation), so the   #
   #                        time series becomes drier.                                     #
   #                 r+XXX: Similar to above, but make the time series wetter.             #
   #                                                                                       #
   # tTTT         -- Temperature scenarios, where tTTTT means:                             #
   #                 t+000: No change in temperature                                       #
   #                 t-YYY: Change temperature by -Y.YY Kelvin.  Keep relative humidity    #
   #                        the same and correct specific humidity.                        #
   #                 t+YYY: Change temperature by +Y.YY Kelvin.  Keep relative humidity    #
   #                        the same and correct by +Y.YY Kelvin.                          #
   #                 r+XXX: Similar to above, but make the time series wetter.             #
   #---------------------------------------------------------------------------------------#
   #----- Find out which scenario to use. -------------------------------------------------#
   fullscen="${metmain}/met_driver/${scentype}/${iscenario}"
   #---------------------------------------------------------------------------------------#
   #     Determine which meteorological data set to use.  Default is the Sheffield/NCEP    #
   # dataset, otherwise the site-level tower data is used.                                 #
   #---------------------------------------------------------------------------------------#
   case ${metdriver} in
   Bananal)
      metdriverdb="${fullscen}/Bananal/Bananal_HEADER"
      metcyc1=2004
      metcycf=2006
      imetavg=1
      ;;
   Brasilia)
      metdriverdb="${fullscen}/Brasilia/Brasilia_HEADER"
      metcyc1=2006
      metcycf=2012
      imetavg=1
      ;;
   Caxiuana)
      metdriverdb="${fullscen}/Caxiuana/Caxiuana_HEADER"
      metcyc1=1999
      metcycf=2003
      imetavg=1
      ;;
   ERAINT_CHIRPS)
      metdriverdb="${fullscen}/${iscenario}_HEADER"
      metcyc1=1981
      metcycf=2017
      imetavg=2
      ;;
   ERAINT_NATIVE)
      metdriverdb="${fullscen}/${iscenario}_HEADER"
      metcyc1=1979
      metcycf=2017
      imetavg=2
      ;;
   ERAINT_MSWEP2)
      metdriverdb="${fullscen}/${iscenario}_HEADER"
      metcyc1=1979
      metcycf=2016
      imetavg=2
      ;;
   Fazenda_Nossa_Senhora)
      metdriverdb="${fullscen}/Fazenda_Nossa_Senhora/Fazenda_Nossa_Senhora_HEADER"
      metcyc1=1999
      metcycf=2002
      imetavg=1
      ;;
   Harvard)
      metdriverdb="${fullscen}/Harvard/Harvard_HEADER"
      metcyc1=1992
      metcycf=2003
      imetavg=3
      ;;
   Laegern)
      metdriverdb="${fullscen}/Laegern/Laegern_HEADER"
      metcyc1=2006
      metcycf=2016
      imetavg=1
      ;;
   Manaus_Km34)
      metdriverdb="${fullscen}/Manaus_Km34/Manaus_Km34_HEADER"
      metcyc1=1999
      metcycf=2006
      imetavg=1
      ;;
   MERRA2_CHIRPS)
      metdriverdb="${fullscen}/${iscenario}_HEADER"
      metcyc1=1981
      metcycf=2017
      imetavg=3
      ;;
   MERRA2_MSWEP2)
      metdriverdb="${fullscen}/${iscenario}_HEADER"
      metcyc1=1980
      metcycf=2016
      imetavg=3
      ;;
   MERRA2_NATIVE)
      metdriverdb="${fullscen}/${iscenario}_HEADER"
      metcyc1=1980
      metcycf=2017
      imetavg=3
      ;;
   Natal)
      metdriverdb="${fullscen}/Natal/Natal_HEADER"
      metcyc1=2009
      metcycf=2012
      imetavg=1
      ;;
   Paracou)
      metdriverdb="${fullscen}/Paracou/Paracou_HEADER"
      metcyc1=2004
      metcycf=2014
      imetavg=1
      ;;
   Pe-de-Gigante)
      metdriverdb="${fullscen}/Pe-de-Gigante/Pe-de-Gigante_HEADER"
      metcyc1=2001
      metcycf=2003
      imetavg=1
      ;;
   Petrolina)
      metdriverdb="${fullscen}/Petrolina/Petrolina_HEADER"
      metcyc1=2004
      metcycf=2012
      imetavg=1
      ;;
   PGMF3_CHIRPS)
      metdriverdb="${fullscen}/${iscenario}_HEADER"
      metcyc1=1981
      metcycf=2016
      imetavg=3
      ;;
   PGMF3_MSWEP2)
      metdriverdb="${fullscen}/${iscenario}_HEADER"
      metcyc1=1979
      metcycf=2016
      imetavg=3
      ;;
   PGMF3_NATIVE)
      metdriverdb="${fullscen}/${iscenario}_HEADER"
      metcyc1=1979
      metcycf=2016
      imetavg=3
      ;;
   Rebio_Jaru)
      metdriverdb="${fullscen}/Rebio_Jaru/Rebio_Jaru_HEADER"
      metcyc1=1999
      metcycf=2002
      imetavg=1
      ;;
   Santarem_Km67)
      metdriverdb="${fullscen}/Santarem_Km67/Santarem_Km67_HEADER"
      metcyc1=2001
      metcycf=2011
      imetavg=1
      ;;
   Santarem_Km77)
      metdriverdb="${fullscen}/Santarem_Km77/Santarem_Km77_HEADER"
      metcyc1=2001
      metcycf=2005
      imetavg=1
      ;;
   Santarem_Km83)
      metdriverdb="${fullscen}/Santarem_Km83/Santarem_Km83_HEADER"
      metcyc1=2000
      metcycf=2003
      imetavg=1
      ;;
   Sheffield)
      metdriverdb=${fullscen}/${shefhead}
      metcyc1=1969
      metcycf=2008
      imetavg=2
      ;;
   Tanguro_Burn)
      metdriverdb="${fullscen}/Tanguro_Burn/Tanguro_Burn_HEADER"
      metcyc1=2008
      metcycf=2018
      imetavg=1
      ;;
   Tanguro_Ctrl)
      metdriverdb="${fullscen}/Tanguro_Ctrl/Tanguro_Ctrl_HEADER"
      metcyc1=2008
      metcycf=2018
      imetavg=1
      ;;
   Tonzi)
      metdriverdb="${fullscen}/Tonzi/Tonzi_HEADER"
      metcyc1=2000
      metcycf=2010
      imetavg=1
      ;;
   WFDEI_CHIRPS)
      metdriverdb="${fullscen}/${iscenario}_HEADER"
      metcyc1=1981
      metcycf=2016
      imetavg=1
      ;;
   WFDEI_CRUP)
      metdriverdb="${fullscen}/${iscenario}_HEADER"
      metcyc1=1979
      metcycf=2016
      imetavg=1
      ;;
   WFDEI_GPCC)
      metdriverdb="${fullscen}/${iscenario}_HEADER"
      metcyc1=1979
      metcycf=2016
      imetavg=1
      ;;
   WFDEI_MSWEP2)
      metdriverdb="${fullscen}/${iscenario}_HEADER"
      metcyc1=1979
      metcycf=2016
      imetavg=1
      ;;
   *)
      echo "Met driver: ${metdriver}"
      echo "Sorry, this met driver is not valid for regular runs"
      exit 85
      ;;
   esac
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Correct years so it is not tower-based or Sheffield.                              #
   #---------------------------------------------------------------------------------------#
   case ${iscenario} in
   default|eft|shr|sheffield|WFDEI*|ERAINT*|MERRA2*|PGMF3*)
      echo "Nothing" > /dev/null
      ;;
   *)
      metcyc1=1972
      metcycf=2012
      imetavg=1
      ;;
   esac
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Set the land use data set.                                                        #
   #---------------------------------------------------------------------------------------#
   case ${ianthdataset} in
   glu-331)
      case ${polyiata} in
      tzi|zmh|nqn|hvd|wch|tqh|l[0-5][0-3])
         ludatabase="${lumain}/glu/glu-"
         ;;
      *)
         ludatabase="${lumain}/glu-3.3.1/glu-3.3.1-"
         ;;
      esac
      ;;
   glu-sa1)
      case ${polyiata} in
      tzi|zmh|nqn|hvd|wch|tqh|l[0-5][0-3])
         ludatabase="${lumain}/glu/glu-"
         ;;
      *)
         ludatabase="${lumain}/glu-3.3.1+sa1.bau/glu-3.3.1+sa1.bau-"
         ;;
      esac
      ;;
   glu-sag)
      case ${polyiata} in
      tzi|zmh|nqn|hvd|wch|tqh|l[0-5][0-3])
         ludatabase="${lumain}/glu/glu-"
         ;;
      *)
         ludatabase="${lumain}/glu-3.3.1+sa1.gov/glu-3.3.1+sa1.gov-"
         ;;
      esac
      ;;
   glu-sa2)
      case ${polyiata} in
      tzi|zmh|nqn|hvd|wch|tqh|l[0-5][0-3])
         ludatabase="${lumain}/glu/glu-"
         ;;
      *)
         ludatabase="${lumain}/glu-3.3.1+sa2.bau/glu-3.3.1+sa2.bau-"
         ;;
      esac
      ;;
   lurcp26)
      ludatabase="${lumain}/luh-1.1+rcp26_image/luh-1.1+rcp26_image-"
      ;;
   lurcp45)
      ludatabase="${lumain}/luh-1.1+rcp45_minicam/luh-1.1+rcp45_minicam-"
      ;;
   lurcp60)
      ludatabase="${lumain}/luh-1.1+rcp60_aim/luh-1.1+rcp60_aim-"
      ;;
   lurcp85)
      ludatabase="${lumain}/luh-1.1+rcp85_message/luh-1.1+rcp85_message-"
      ;;
   *)
      #------------------------------------------------------------------------------------#
      #     Stop the script if anthropogenic dataset is invalid and this is a simulation   #
      # with anthropogenic disturbance.                                                    #
      #------------------------------------------------------------------------------------#
      if [ ${ianthdisturb} -eq 1 ]
      then
         echo " Polygon:       ${polyname}"
         echo " IATA:          ${polyiata}"
         echo " IANTH_DATASET: ${iage}"
         echo 'Invalid anthropogenic disturbance data set!'
         exit 53
      fi
      #------------------------------------------------------------------------------------#
      ;;
   esac
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Define whether we use the met cycle to define the first and last year, or the     #
   # default year.                                                                         #
   #---------------------------------------------------------------------------------------#
   if [ ${yeara} -eq 0 ]
   then
      thisyeara=${metcyc1}
   else
      thisyeara=${yeara}
   fi
   if [ ${yearz} -eq 0 ]
   then
      thisyearz=${metcycf}
   else
      thisyearz=${yearz}
   fi
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Determine which soil profile to use.  We have eight categories (A-H), and the     #
   # soil resolution.                                                                      #
   #---------------------------------------------------------------------------------------#
   polyslz1=""
   polyslz2=""
   polyslz3=""
   polyslz4=""
   polyslz5=""
   polyslz6=""
   polyslz7=""
   polyslm1=""
   polyslm2=""
   polyslm3=""
   polyslm4=""
   polyslm5=""
   polyslm6=""
   polyslm7=""
   polyslt1=""
   polyslt2=""
   polyslt3=""
   polyslt4=""
   polyslt5=""
   polyslt6=""
   polyslt7=""
   case ${slzres} in
   0)
      case ${polydepth} in
      A)
         polynzg=16
         polyslz1="-1.250,-1.154,-1.059,-0.966,-0.875,-0.785,-0.697,-0.612,-0.529,-0.448,"
         polyslz2="-0.370,-0.295,-0.224,-0.156,-0.095,-0.040"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslm2=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslt2=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000"
         ;;
      B)
         polynzg=16
         polyslz1="-2.000,-1.826,-1.657,-1.492,-1.333,-1.179,-1.030,-0.888,-0.752,-0.623,"
         polyslz2="-0.501,-0.388,-0.283,-0.188,-0.106,-0.040"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslm2=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslt2=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000"
         ;;
      C)
         polynzg=16
         polyslz1="-3.000,-2.713,-2.437,-2.171,-1.917,-1.674,-1.443,-1.225,-1.019,-0.828,"
         polyslz2="-0.651,-0.490,-0.346,-0.221,-0.118,-0.040"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslm2=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslt2=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000"
         ;;
      D)
         polynzg=16
         polyslz1="-4.000,-3.593,-3.204,-2.833,-2.481,-2.147,-1.832,-1.538,-1.265,-1.013,"
         polyslz2="-0.784,-0.579,-0.400,-0.248,-0.126,-0.040"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslm2=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslt2=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000"
         ;;
      E)
         polynzg=16
         polyslz1="-4.500,-4.032,-3.584,-3.159,-2.757,-2.377,-2.021,-1.689,-1.382,-1.101,"
         polyslz2="-0.846,-0.620,-0.424,-0.260,-0.130,-0.040"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslm2=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslt2=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000"
         ;;
      F)
         polynzg=16
         polyslz1="-6.000,-5.339,-4.714,-4.123,-3.567,-3.048,-2.566,-2.121,-1.714,-1.347,"
         polyslz2="-1.019,-0.733,-0.490,-0.291,-0.140,-0.040"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslm2=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslt2=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000"
         ;;
      G)
         polynzg=16
         polyslz1="-7.000,-6.207,-5.458,-4.755,-4.096,-3.483,-2.917,-2.397,-1.925,-1.501,"
         polyslz2="-1.126,-0.802,-0.529,-0.310,-0.145,-0.040"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslm2=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslt2=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000"
         ;;
      H)
         polynzg=16
         polyslz1="-8.000,-7.072,-6.198,-5.380,-4.617,-3.910,-3.259,-2.664,-2.127,-1.648,"
         polyslz2="-1.228,-0.866,-0.566,-0.326,-0.150,-0.040"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslm2=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslt2=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000"
         ;;
      I)
         polynzg=16
         polyslz1="-10.50,-9.223,-8.029,-6.919,-5.891,-4.946,-4.084,-3.305,-2.609,-1.995,"
         polyslz2="-1.464,-1.015,-0.648,-0.364,-0.161,-0.040"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslm2=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslt2=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000"
         ;;
      P)
         polynzg=7
         polyslz1="-0.580,-0.420,-0.300,-0.200,-0.120,-0.060,-0.020"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000"
         ;;
      Q)
         polynzg=8
         polyslz1="-0.830,-0.630,-0.450,-0.300,-0.200,-0.120,-0.060,-0.020"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000"
         ;;
      U)
         polynzg=6
         polyslz1="-0.400,-0.300,-0.200,-0.120,-0.060,-0.020"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000"
         ;;
      V)
         polynzg=6
         polyslz1="-0.420,-0.300,-0.200,-0.120,-0.060,-0.020"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000"
         ;;
      W)
         polynzg=7
         polyslz1="-0.510,-0.400,-0.300,-0.200,-0.120,-0.060,-0.020"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000"
         ;;
      X)
         polynzg=8
         polyslz1="-0.670,-0.520,-0.400,-0.300,-0.200,-0.120,-0.060,-0.020"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000"
         ;;
      Y)
         polynzg=8
         polyslz1="-0.700,-0.550,-0.400,-0.300,-0.200,-0.120,-0.060,-0.020"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000"
         ;;
      Z)
         polynzg=8
         polyslz1="-0.750,-0.600,-0.450,-0.300,-0.200,-0.120,-0.060,-0.020"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000"
         ;;
      *)
         polynzg=16
         polyslz1="-8.000,-7.072,-6.198,-5.380,-4.617,-3.910,-3.259,-2.664,-2.127,-1.648,"
         polyslz2="-1.228,-0.866,-0.566,-0.326,-0.150,-0.040"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslm2=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslt2=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000"
         ;;
      esac
      ;;
   1)
      case ${polydepth} in
      A)
         polynzg=16
         polyslz1="-1.250,-1.135,-1.024,-0.917,-0.814,-0.715,-0.620,-0.530,-0.445,-0.364,"
         polyslz2="-0.289,-0.221,-0.158,-0.103,-0.056,-0.020"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslm2=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslt2=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000"
         ;;
      B)
         polynzg=16
         polyslz1="-2.000,-1.797,-1.602,-1.417,-1.240,-1.073,-0.916,-0.769,-0.632,-0.507,"
         polyslz2="-0.392,-0.290,-0.200,-0.124,-0.063,-0.020"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslm2=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslt2=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000"
         ;;
      C)
         polynzg=16
         polyslz1="-3.000,-2.670,-2.357,-2.061,-1.784,-1.524,-1.283,-1.061,-0.857,-0.673,"
         polyslz2="-0.510,-0.367,-0.245,-0.146,-0.070,-0.020"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslm2=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslt2=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000"
         ;;
      D)
         polynzg=16
         polyslz1="-4.000,-3.536,-3.099,-2.690,-2.308,-1.955,-1.629,-1.332,-1.064,-0.824,"
         polyslz2="-0.614,-0.433,-0.283,-0.163,-0.075,-0.020"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslm2=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslt2=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000"
         ;;
      E)
         polynzg=16
         polyslz1="-4.500,-3.967,-3.467,-3.000,-2.565,-2.164,-1.797,-1.462,-1.162,-0.895,"
         polyslz2="-0.662,-0.464,-0.300,-0.171,-0.077,-0.020"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslm2=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslt2=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000"
         ;;
      F)
         polynzg=16
         polyslz1="-6.000,-5.254,-4.559,-3.914,-3.320,-2.776,-2.282,-1.837,-1.442,-1.095,"
         polyslz2="-0.798,-0.548,-0.346,-0.192,-0.083,-0.020"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslm2=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslt2=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000"
         ;;
      G)
         polynzg=16
         polyslz1="-7.000,-6.108,-5.279,-4.514,-3.812,-3.172,-2.593,-2.076,-1.618,-1.221,"
         polyslz2="-0.881,-0.600,-0.374,-0.204,-0.087,-0.020"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslm2=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslt2=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000"
         ;;
      H)
         polynzg=16
         polyslz1="-8.000,-6.959,-5.995,-5.108,-4.296,-3.560,-2.897,-2.307,-1.789,-1.340,"
         polyslz2="-0.961,-0.648,-0.400,-0.215,-0.089,-0.020"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslm2=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslt2=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000"
         ;;
      I)
         polynzg=16
         polyslz1="-10.50,-9.076,-7.766,-6.569,-5.482,-4.504,-3.631,-2.862,-2.194,-1.622,"
         polyslz2="-1.145,-0.759,-0.458,-0.239,-0.096,-0.020"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslm2=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslt2=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000"
         ;;
      P)
         polynzg=7
         polyslz1="-0.580,-0.420,-0.300,-0.200,-0.120,-0.060,-0.020"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000"
         ;;
      Q)
         polynzg=8
         polyslz1="-0.830,-0.630,-0.450,-0.300,-0.200,-0.120,-0.060,-0.020"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000"
         ;;
      U)
         polynzg=6
         polyslz1="-0.400,-0.300,-0.200,-0.120,-0.060,-0.020"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000"
         ;;
      V)
         polynzg=6
         polyslz1="-0.420,-0.300,-0.200,-0.120,-0.060,-0.020"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000"
         ;;
      W)
         polynzg=7
         polyslz1="-0.510,-0.400,-0.300,-0.200,-0.120,-0.060,-0.020"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000"
         ;;
      X)
         polynzg=8
         polyslz1="-0.670,-0.520,-0.400,-0.300,-0.200,-0.120,-0.060,-0.020"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000"
         ;;
      Y)
         polynzg=8
         polyslz1="-0.700,-0.550,-0.400,-0.300,-0.200,-0.120,-0.060,-0.020"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000"
         ;;
      Z)
         polynzg=8
         polyslz1="-0.750,-0.600,-0.450,-0.300,-0.200,-0.120,-0.060,-0.020"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000"
         ;;
      *)
         polynzg=16
         polyslz1="-8.000,-6.959,-5.995,-5.108,-4.296,-3.560,-2.897,-2.307,-1.789,-1.340,"
         polyslz2="-0.961,-0.648,-0.400,-0.215,-0.089,-0.020"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslm2=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslt2=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000"
         ;;
      esac
      ;;
   2)
      case ${polydepth} in
      A)
         polynzg="30"
         polyslz1="-1.329,-1.246,-1.168,-1.093,-1.022,-0.955,-0.890,-0.829,-0.770,-0.714,"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslz2="-0.661,-0.611,-0.563,-0.517,-0.473,-0.432,-0.392,-0.354,-0.318,-0.284,"
         polyslm2=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslt2=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslz3="-0.252,-0.221,-0.191,-0.163,-0.136,-0.111,-0.086,-0.063,-0.041,-0.020"
         polyslm3=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000"
         polyslt3=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000"
         ;;
      B)
         polynzg="37"
         polyslz1="-2.033,-1.917,-1.806,-1.701,-1.601,-1.506,-1.415,-1.329,-1.246,-1.168,"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslz2="-1.093,-1.022,-0.955,-0.890,-0.829,-0.770,-0.714,-0.661,-0.611,-0.563,"
         polyslm2=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslt2=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslz3="-0.517,-0.473,-0.432,-0.392,-0.354,-0.318,-0.284,-0.252,-0.221,-0.191,"
         polyslm3=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslt3=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslz4="-0.163,-0.136,-0.111,-0.086,-0.063,-0.041,-0.020"
         polyslm4=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000"
         polyslt4=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000"
         ;;
      C)
         polynzg="44"
         polyslz1="-3.023,-2.860,-2.705,-2.557,-2.416,-2.282,-2.154,-2.033,-1.917,-1.806,"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslz2="-1.701,-1.601,-1.506,-1.415,-1.329,-1.246,-1.168,-1.093,-1.022,-0.955,"
         polyslm2=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslt2=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslz3="-0.890,-0.829,-0.770,-0.714,-0.661,-0.611,-0.563,-0.517,-0.473,-0.432,"
         polyslm3=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslt3=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslz4="-0.392,-0.354,-0.318,-0.284,-0.252,-0.221,-0.191,-0.163,-0.136,-0.111,"
         polyslm4=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslt4=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslz5="-0.086,-0.063,-0.041,-0.020"
         polyslm5=" 1.000, 1.000, 1.000, 1.000"
         polyslt5=" 0.000, 0.000, 0.000, 0.000"
         ;;
      D)
         polynzg="50"
         polyslz1="-4.187,-3.969,-3.761,-3.562,-3.374,-3.194,-3.023,-2.860,-2.705,-2.557,"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslz2="-2.416,-2.282,-2.154,-2.033,-1.917,-1.806,-1.701,-1.601,-1.506,-1.415,"
         polyslm2=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslt2=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslz3="-1.329,-1.246,-1.168,-1.093,-1.022,-0.955,-0.890,-0.829,-0.770,-0.714,"
         polyslm3=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslt3=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslz4="-0.661,-0.611,-0.563,-0.517,-0.473,-0.432,-0.392,-0.354,-0.318,-0.284,"
         polyslm4=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslt4=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslz5="-0.252,-0.221,-0.191,-0.163,-0.136,-0.111,-0.086,-0.063,-0.041,-0.020"
         polyslm5=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000"
         polyslt5=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000"
         ;;
      E)
         polynzg="52"
         polyslz1="-4.657,-4.416,-4.187,-3.969,-3.761,-3.562,-3.374,-3.194,-3.023,-2.860,"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslz2="-2.705,-2.557,-2.416,-2.282,-2.154,-2.033,-1.917,-1.806,-1.701,-1.601,"
         polyslm2=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslt2=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslz3="-1.506,-1.415,-1.329,-1.246,-1.168,-1.093,-1.022,-0.955,-0.890,-0.829,"
         polyslm3=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslt3=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslz4="-0.770,-0.714,-0.661,-0.611,-0.563,-0.517,-0.473,-0.432,-0.392,-0.354,"
         polyslm4=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslt4=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslz5="-0.318,-0.284,-0.252,-0.221,-0.191,-0.163,-0.136,-0.111,-0.086,-0.063,"
         polyslm5=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslt5=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslz6="-0.041,-0.020"
         polyslm6=" 1.000, 1.000"
         polyslt6=" 0.000, 0.000"
         ;;
      F)
         polynzg="58"
         polyslz1="-6.157,-5.907,-5.657,-5.407,-5.157,-4.907,-4.657,-4.416,-4.187,-3.969,"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslz2="-3.761,-3.562,-3.374,-3.194,-3.023,-2.860,-2.705,-2.557,-2.416,-2.282,"
         polyslm2=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslt2=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslz3="-2.154,-2.033,-1.917,-1.806,-1.701,-1.601,-1.506,-1.415,-1.329,-1.246,"
         polyslm3=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslt3=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslz4="-1.168,-1.093,-1.022,-0.955,-0.890,-0.829,-0.770,-0.714,-0.661,-0.611,"
         polyslm4=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslt4=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslz5="-0.563,-0.517,-0.473,-0.432,-0.392,-0.354,-0.318,-0.284,-0.252,-0.221,"
         polyslm5=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslt5=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslz6="-0.191,-0.163,-0.136,-0.111,-0.086,-0.063,-0.041,-0.020"
         polyslm6=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000"
         polyslt6=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000"
         ;;
      G)
         polynzg="62"
         polyslz1="-7.157,-6.907,-6.657,-6.407,-6.157,-5.907,-5.657,-5.407,-5.157,-4.907,"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslz2="-4.657,-4.416,-4.187,-3.969,-3.761,-3.562,-3.374,-3.194,-3.023,-2.860,"
         polyslm2=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslt2=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslz3="-2.705,-2.557,-2.416,-2.282,-2.154,-2.033,-1.917,-1.806,-1.701,-1.601,"
         polyslm3=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslt3=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslz4="-1.506,-1.415,-1.329,-1.246,-1.168,-1.093,-1.022,-0.955,-0.890,-0.829,"
         polyslm4=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslt4=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslz5="-0.770,-0.714,-0.661,-0.611,-0.563,-0.517,-0.473,-0.432,-0.392,-0.354,"
         polyslm5=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslt5=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslz6="-0.318,-0.284,-0.252,-0.221,-0.191,-0.163,-0.136,-0.111,-0.086,-0.063,"
         polyslm6=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslt6=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslz7="-0.041,-0.020"
         polyslm7=" 1.000, 1.000"
         polyslt7=" 0.000, 0.000"
         ;;
      H)
         polynzg="66"
         polyslz1="-8.157,-7.907,-7.657,-7.407,-7.157,-6.907,-6.657,-6.407,-6.157,-5.907,"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslz2="-5.657,-5.407,-5.157,-4.907,-4.657,-4.416,-4.187,-3.969,-3.761,-3.562,"
         polyslm2=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslt2=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslz3="-3.374,-3.194,-3.023,-2.860,-2.705,-2.557,-2.416,-2.282,-2.154,-2.033,"
         polyslm3=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslt3=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslz4="-1.917,-1.806,-1.701,-1.601,-1.506,-1.415,-1.329,-1.246,-1.168,-1.093,"
         polyslm4=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslt4=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslz5="-1.022,-0.955,-0.890,-0.829,-0.770,-0.714,-0.661,-0.611,-0.563,-0.517,"
         polyslm5=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslt5=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslz6="-0.473,-0.432,-0.392,-0.354,-0.318,-0.284,-0.252,-0.221,-0.191,-0.163,"
         polyslm6=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslt6=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslz7="-0.136,-0.111,-0.086,-0.063,-0.041,-0.020"
         polyslm7=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000"
         polyslt7=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000"
         ;;
      *)
         polynzg="66"
         polyslz1="-8.157,-7.907,-7.657,-7.407,-7.157,-6.907,-6.657,-6.407,-6.157,-5.907,"
         polyslm1=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslt1=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslz2="-5.657,-5.407,-5.157,-4.907,-4.657,-4.416,-4.187,-3.969,-3.761,-3.562,"
         polyslm2=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslt2=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslz3="-3.374,-3.194,-3.023,-2.860,-2.705,-2.557,-2.416,-2.282,-2.154,-2.033,"
         polyslm3=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslt3=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslz4="-1.917,-1.806,-1.701,-1.601,-1.506,-1.415,-1.329,-1.246,-1.168,-1.093,"
         polyslm4=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslt4=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslz5="-1.022,-0.955,-0.890,-0.829,-0.770,-0.714,-0.661,-0.611,-0.563,-0.517,"
         polyslm5=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslt5=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslz6="-0.473,-0.432,-0.392,-0.354,-0.318,-0.284,-0.252,-0.221,-0.191,-0.163,"
         polyslm6=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,"
         polyslt6=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,"
         polyslz7="-0.136,-0.111,-0.086,-0.063,-0.041,-0.020"
         polyslm7=" 1.000, 1.000, 1.000, 1.000, 1.000, 1.000"
         polyslt7=" 0.000, 0.000, 0.000, 0.000, 0.000, 0.000"
         ;;
      esac
      ;;
   esac
   #---------------------------------------------------------------------------------------#



   #----- Check whether to use SFILIN as restart or history. ------------------------------#
   if [ ${runt} == "INITIAL" ] && [ ${forcehisto} -eq 1 ]
   then
      runt="HISTORY"
      year=${yearh}
      month=${monthh}
      date=${dateh}
      time=${timeh}
      thissfilin=${fullygrown}
   elif [ ${runt} == "INITIAL" ] && [ ${initmode} -eq 5 ]
   then
      if [ ! -s ${restart} ]
      then
         echo " Directory restart does not exist!"
         echo " Change the variable restart at the beginning of the script"
         exit 44
      else
         runt="INITIAL"
         thissfilin=${restart}
      fi
   elif [ ${runt} == "INITIAL" ] && [ ${initmode} -eq 6 ]
   then
      thissfilin=${fullygrown}



      #------------------------------------------------------------------------------------#
      #    Find the biometric files.  This has been checked in spawn_poly.sh so they are   #
      # correct.  Add a dummy name in case this is not supposed to be a biomass            #
      # initialisation run.                                                                #
      #------------------------------------------------------------------------------------#
      case ${biotype} in
      0)
         #----- isizepft controls everything, and iage is ignored. ------------------------#
         case ${isizepft} in
         0)
            #----- Frankeinstein's under storey. ------------------------------------------#
            thissfilin="${bioinit}/${polyiata}_default."
            ;;
         1)
            #----- No under storey. -------------------------------------------------------#
            thissfilin="${bioinit}/${polyiata}_nounder."
            ;;
         2)
            #----- ALS initialisation. ----------------------------------------------------#
            thissfilin="${bioinit}/${polyiata}_alsinit."
            ;;
         *)
            #----- Invalid option. Stop the script. ---------------------------------------#
            echo " Polygon:  ${polyname}"
            echo " IATA:     ${polyiata}"
            echo " ISIZEPFT: ${isizepft}"
            echo "This IATA cannot be used by biomass initialisation with this ISIZEPFT!"
            exit 57
            ;;
         esac
         #---------------------------------------------------------------------------------#
         ;;

      1)
         #---------------------------------------------------------------------------------#
         #    'isizepft' controls how many PFTs to use.                                    #
         #---------------------------------------------------------------------------------#
         case ${isizepft} in 
         0|5)
            pftname="pft05"
            ;;
         2)
            pftname="pft02"
            ;;
         esac
         #---------------------------------------------------------------------------------#

         #---------------------------------------------------------------------------------#
         #     'iage' controls how many patches to use.                                    #
         #---------------------------------------------------------------------------------#
         case ${iage} in
         1)
            agename="age01"
            ;;
         *)
            agename="age00"
            ;;
         esac
         #---------------------------------------------------------------------------------#




         #---------------------------------------------------------------------------------#
         #      Check whether the site has the PFT and age structure.                      #
         #---------------------------------------------------------------------------------#
         case ${polyiata} in
         hvd|s77|fns|cau|and|par|tap|dcm)
            thissfilin="${bioinit}/${polyiata}_default."
            ;;
         cax|s67|s83|m34|gyf|pdg|rja|pnz|ban)
            thissfilin="${bioinit}/${polyiata}_${pftname}+${agename}."
            ;;
         *)
            echo " Polygon:  ${polyname}"
            echo " IATA:     ${polyiata}"
            echo " IAGE:     ${iage}"
            echo " ISIZEPFT: ${isizepft}"
            echo "This IATA cannot be used by biomass initialisation with this ISIZEPFT!"
            exit 59
            ;;
         esac
         #---------------------------------------------------------------------------------#
         ;;
      2)
         #---------------------------------------------------------------------------------#
         #     ALS initialisation. ISIZEPFT has disturbance history information.           #
         #---------------------------------------------------------------------------------#
         case ${polyiata} in
         l[0-5][0-3])
            thissfilin="${alsinit}/${polyiata}."
            ;;
         *)
            thissfilin="${alsinit}/${polyiata}_${isizepft}."
            ;;
         esac
         #---------------------------------------------------------------------------------#
         ;;
      esac
      #------------------------------------------------------------------------------------#
   else
      thissfilin=${here}/${polyname}/histo/${polyname}
   fi
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



   #---------------------------------------------------------------------------------------#
   #     Replace the flags in ED2IN.                                                       #
   #---------------------------------------------------------------------------------------#
   ED2IN="${here}/${polyname}/ED2IN"
   sed -i~ s@paththere@${here}@g                 ${ED2IN}
   sed -i~ s@myinpmain@${inpmain}@g              ${ED2IN}
   sed -i~ s@myyeara@${thisyeara}@g              ${ED2IN}
   sed -i~ s@mymontha@${montha}@g                ${ED2IN}
   sed -i~ s@mydatea@${datea}@g                  ${ED2IN}
   sed -i~ s@mytimea@${timea}@g                  ${ED2IN}
   sed -i~ s@myyearz@${thisyearz}@g              ${ED2IN}
   sed -i~ s@mymonthz@${monthz}@g                ${ED2IN}
   sed -i~ s@mydatez@${datez}@g                  ${ED2IN}
   sed -i~ s@mytimez@${timez}@g                  ${ED2IN}
   sed -i~ s@mydtlsm@${dtlsm}@g                  ${ED2IN}
   sed -i~ s@mymonyrstep@${monyrstep}@g          ${ED2IN}
   sed -i~ s@thispoly@${polyname}@g              ${ED2IN}
   sed -i~ s@plonflag@${polylon}@g               ${ED2IN}
   sed -i~ s@platflag@${polylat}@g               ${ED2IN}
   sed -i~ s@timehhhh@${time}@g                  ${ED2IN}
   sed -i~ s@datehhhh@${date}@g                  ${ED2IN}
   sed -i~ s@monthhhh@${month}@g                 ${ED2IN}
   sed -i~ s@yearhhhh@${year}@g                  ${ED2IN}
   sed -i~ s@myunitstate@${iunitstate}@g         ${ED2IN}
   sed -i~ s@myinitmode@${initmode}@g            ${ED2IN}
   sed -i~ s@mysfilin@${thissfilin}@g            ${ED2IN}
   sed -i~ s@mytrees@${pfts}@g                   ${ED2IN}
   sed -i~ s@mycrop@${crop}@g                    ${ED2IN}
   sed -i~ s@mypasture@${pasture}@g              ${ED2IN}
   sed -i~ s@myplantation@${plantation}@g        ${ED2IN}
   sed -i~ s@myiphen@${iphen}@g                  ${ED2IN}
   sed -i~ s@myallom@${iallom}@g                 ${ED2IN}
   sed -i~ s@myeconomics@${ieconomics}@g         ${ED2IN}
   sed -i~ s@mygrass@${igrass}@g                 ${ED2IN}
   sed -i~ s@myisoilflg@${polyisoil}@g           ${ED2IN}
   sed -i~ s@mynslcon@${polyntext}@g             ${ED2IN}
   sed -i~ s@myslxsand@${polysand}@g             ${ED2IN}
   sed -i~ s@myslxclay@${polyclay}@g             ${ED2IN}
   sed -i~ s@mysoilbc@${polysoilbc}@g            ${ED2IN}
   sed -i~ s@mysldrain@${polysldrain}@g          ${ED2IN}
   sed -i~ s@mysoilcol@${polycol}@g              ${ED2IN}
   sed -i~ s@mynzg@${polynzg}@g                  ${ED2IN}
   sed -i~ s@mymetdriverdb@${metdriverdb}@g      ${ED2IN}
   sed -i~ s@mymetcyc1@${metcyc1}@g              ${ED2IN}
   sed -i~ s@mymetcycf@${metcycf}@g              ${ED2IN}
   sed -i~ s@mytoler@${toler}@g                  ${ED2IN}
   sed -i~ s@RUNFLAG@${runt}@g                   ${ED2IN}
   sed -i~ s@myiphysiol@${iphysiol}@g            ${ED2IN}
   sed -i~ s@myvmfactc3@${vmfactc3}@g            ${ED2IN}
   sed -i~ s@myvmfactc4@${vmfactc4}@g            ${ED2IN}
   sed -i~ s@mymphototrc3@${mphototrc3}@g        ${ED2IN}
   sed -i~ s@mymphototec3@${mphototec3}@g        ${ED2IN}
   sed -i~ s@mymphotoc4@${mphotoc4}@g            ${ED2IN}
   sed -i~ s@mybphotoblc3@${bphotoblc3}@g        ${ED2IN}
   sed -i~ s@mybphotonlc3@${bphotonlc3}@g        ${ED2IN}
   sed -i~ s@mybphotoc4@${bphotoc4}@g            ${ED2IN}
   sed -i~ s@mykwgrass@${kwgrass}@g              ${ED2IN}
   sed -i~ s@mykwtree@${kwtree}@g                ${ED2IN}
   sed -i~ s@mygammac3@${gammac3}@g              ${ED2IN}
   sed -i~ s@mygammac4@${gammac4}@g              ${ED2IN}
   sed -i~ s@myd0grass@${d0grass}@g              ${ED2IN}
   sed -i~ s@myd0tree@${d0tree}@g                ${ED2IN}
   sed -i~ s@myalphac3@${alphac3}@g              ${ED2IN}
   sed -i~ s@myalphac4@${alphac4}@g              ${ED2IN}
   sed -i~ s@myklowco2@${klowco2}@g              ${ED2IN}
   sed -i~ s@mydecomp@${decomp}@g                ${ED2IN}
   sed -i~ s@myrrffact@${rrffact}@g              ${ED2IN}
   sed -i~ s@mygrowthresp@${growthresp}@g        ${ED2IN}
   sed -i~ s@mylwidthgrass@${lwidthgrass}@g      ${ED2IN}
   sed -i~ s@mylwidthbltree@${lwidthbltree}@g    ${ED2IN}
   sed -i~ s@mylwidthnltree@${lwidthnltree}@g    ${ED2IN}
   sed -i~ s@myq10c3@${q10c3}@g                  ${ED2IN}
   sed -i~ s@myq10c4@${q10c4}@g                  ${ED2IN}
   sed -i~ s@myh2olimit@${h2olimit}@g            ${ED2IN}
   sed -i~ s@mymortscheme@${imortscheme}@g       ${ED2IN}
   sed -i~ s@myddmortconst@${ddmortconst}@g      ${ED2IN}
   sed -i~ s@mycbrscheme@${cbrscheme}@g          ${ED2IN}
   sed -i~ s@mysfclyrm@${isfclyrm}@g             ${ED2IN}
   sed -i~ s@myicanturb@${icanturb}@g            ${ED2IN}
   sed -i~ s@myatmco2@${atmco2}@g                ${ED2IN}
   sed -i~ s@mythcrit@${thcrit}@g                ${ED2IN}
   sed -i~ s@mysmfire@${smfire}@g                ${ED2IN}
   sed -i~ s@myfire@${ifire}@g                   ${ED2IN}
   sed -i~ s@myfuel@${fireparm}@g                ${ED2IN}
   sed -i~ s@mymetavg@${imetavg}@g               ${ED2IN}
   sed -i~ s@mypercol@${ipercol}@g               ${ED2IN}
   sed -i~ s@myrunoff@${runoff}@g                ${ED2IN}
   sed -i~ s@mymetrad@${imetrad}@g               ${ED2IN}
   sed -i~ s@mybranch@${ibranch}@g               ${ED2IN}
   sed -i~ s@mycanrad@${icanrad}@g               ${ED2IN}
   sed -i~ s@myhrzrad@${ihrzrad}@g               ${ED2IN}
   sed -i~ s@mycrown@${crown}@g                  ${ED2IN}
   sed -i~ s@myltransvis@${ltransvis}@g          ${ED2IN}
   sed -i~ s@myltransnir@${ltransnir}@g          ${ED2IN}
   sed -i~ s@mylreflectvis@${lreflectvis}@g      ${ED2IN}
   sed -i~ s@mylreflectnir@${lreflectnir}@g      ${ED2IN}
   sed -i~ s@myorienttree@${orienttree}@g        ${ED2IN}
   sed -i~ s@myorientgrass@${orientgrass}@g      ${ED2IN}
   sed -i~ s@myclumptree@${clumptree}@g          ${ED2IN}
   sed -i~ s@myclumpgrass@${clumpgrass}@g        ${ED2IN}
   sed -i~ s@myigoutput@${igoutput}@g            ${ED2IN}
   sed -i~ s@mygpref@${gpref}@g                  ${ED2IN}
   sed -i~ s@myvegtdyn@${ivegtdyn}@g             ${ED2IN}
   sed -i~ s@myhydroscheme@${ihydro}@g           ${ED2IN}
   sed -i~ s@mystemresp@${istemresp}@g           ${ED2IN}
   sed -i~ s@mystomata@${istomata}@g             ${ED2IN}
   sed -i~ s@myplastic@${iplastic}@g             ${ED2IN}
   sed -i~ s@mycarbonmort@${icarbonmort}@g       ${ED2IN}
   sed -i~ s@myhydromort@${ihydromort}@g         ${ED2IN}
   sed -i~ s@mybigleaf@${ibigleaf}@g             ${ED2IN}
   sed -i~ s@myintegscheme@${integscheme}@g      ${ED2IN}
   sed -i~ s@mynsubeuler@${nsubeuler}@g          ${ED2IN}
   sed -i~ s@myrepro@${irepro}@g                 ${ED2IN}
   sed -i~ s@myubmin@${ubmin}@g                  ${ED2IN}
   sed -i~ s@myugbmin@${ugbmin}@g                ${ED2IN}
   sed -i~ s@myustmin@${ustmin}@g                ${ED2IN}
   sed -i~ s@mygamm@${gamm}@g                    ${ED2IN}
   sed -i~ s@mygamh@${gamh}@g                    ${ED2IN}
   sed -i~ s@mytprandtl@${tprandtl}@g            ${ED2IN}
   sed -i~ s@myribmax@${ribmax}@g                ${ED2IN}
   sed -i~ s@mygndvap@${igndvap}@g               ${ED2IN}
   sed -i~ s@mydtcensus@${dtcensus}@g            ${ED2IN}
   sed -i~ s@myyr1stcensus@${yr1stcensus}@g      ${ED2IN}
   sed -i~ s@mymon1stcensus@${mon1stcensus}@g    ${ED2IN}
   sed -i~ s@myminrecruitdbh@${minrecruitdbh}@g  ${ED2IN}
   sed -i~ s@mytreefall@${treefall}@g            ${ED2IN}
   sed -i~ s@mymaxpatch@${iage}@g                ${ED2IN}
   sed -i~ s@mymaxcohort@${imaxcohort}@g         ${ED2IN}
   sed -i~ s@myanthdisturb@${ianthdisturb}@g     ${ED2IN}
   sed -i~ s@myludatabase@${ludatabase}@g        ${ED2IN}
   sed -i~ s@myslscale@${slscale}@g              ${ED2IN}
   sed -i~ s@myslyrfirst@${slyrfirst}@g          ${ED2IN}
   sed -i~ s@myslnyrs@${slnyrs}@g                ${ED2IN}
   sed -i~ s@mylogging@${logging}@g              ${ED2IN}
   sed -i~ s@myprobharv@${probharv}@g            ${ED2IN}
   sed -i~ s@mydbhharv@${dbhharv}@g              ${ED2IN}
   sed -i~ s@mybioharv@${bioharv}@g              ${ED2IN}
   sed -i~ s@myskidarea@${skidarea}@g            ${ED2IN}
   sed -i~ s@myskidsmall@${skidsmall}@g          ${ED2IN}
   sed -i~ s@myskidlarge@${skidlarge}@g          ${ED2IN}
   sed -i~ s@myfellingsmall@${fellingsmall}@g    ${ED2IN}
   sed -i~ s@myseedharv@${seedharv}@g            ${ED2IN}
   sed -i~ s@mystorharv@${storharv}@g            ${ED2IN}
   sed -i~ s@myleafharv@${leafharv}@g            ${ED2IN}
   #---------------------------------------------------------------------------------------#

   #------ Soil variables. ----------------------------------------------------------------#
   sed -i~ s@myslz1@"${polyslz1}"@g           ${ED2IN}
   sed -i~ s@myslz2@"${polyslz2}"@g           ${ED2IN}
   sed -i~ s@myslz3@"${polyslz3}"@g           ${ED2IN}
   sed -i~ s@myslz4@"${polyslz4}"@g           ${ED2IN}
   sed -i~ s@myslz5@"${polyslz5}"@g           ${ED2IN}
   sed -i~ s@myslz6@"${polyslz6}"@g           ${ED2IN}
   sed -i~ s@myslz7@"${polyslz7}"@g           ${ED2IN}
   sed -i~ s@myslmstr1@"${polyslm1}"@g        ${ED2IN}
   sed -i~ s@myslmstr2@"${polyslm2}"@g        ${ED2IN}
   sed -i~ s@myslmstr3@"${polyslm3}"@g        ${ED2IN}
   sed -i~ s@myslmstr4@"${polyslm4}"@g        ${ED2IN}
   sed -i~ s@myslmstr5@"${polyslm5}"@g        ${ED2IN}
   sed -i~ s@myslmstr6@"${polyslm6}"@g        ${ED2IN}
   sed -i~ s@myslmstr7@"${polyslm7}"@g        ${ED2IN}
   sed -i~ s@mystgoff1@"${polyslt1}"@g        ${ED2IN}
   sed -i~ s@mystgoff2@"${polyslt2}"@g        ${ED2IN}
   sed -i~ s@mystgoff3@"${polyslt3}"@g        ${ED2IN}
   sed -i~ s@mystgoff4@"${polyslt4}"@g        ${ED2IN}
   sed -i~ s@mystgoff5@"${polyslt5}"@g        ${ED2IN}
   sed -i~ s@mystgoff6@"${polyslt6}"@g        ${ED2IN}
   sed -i~ s@mystgoff7@"${polyslt7}"@g        ${ED2IN}
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #   In case this is a multithreaded run, copy executables to each directory.            #
   #---------------------------------------------------------------------------------------#
   case ${n_cpt} in
   1)
      exec_sub="${here}/${polyname}/${execname}"
      cp ${exec_full} ${exec_sub}
      ;;
   *)
      exec_sub=${exec_full}
      ;;
   esac

   #----- Change the callserial.sh file. --------------------------------------------------#
   callserial="${here}/${polyname}/callserial.sh"
   rm -f ${callserial}
   cp -f ${here}/Template/callserial.sh ${callserial}
   sed -i~ s@thisroot@${here}@g          ${callserial}
   sed -i~ s@thispoly@${polyname}@g      ${callserial}
   sed -i~ s@myname@${moi}@g             ${callserial}
   sed -i~ s@myexec@${exec_sub}@g        ${callserial}
   sed -i~ s@mypackdata@${packdatasrc}@g ${callserial}
   sed -i~ s@myscenario@${iscenario}@g   ${callserial}
   sed -i~ s@myscenmain@${scentype}@g    ${callserial}
   sed -i~ s@mycopy@${copy2scratch}@g    ${callserial}
   sed -i~ s@mycpus@${n_cpt}@g           ${callserial}
   sed -i~ s@myoptsrc@${optsrc}@g        ${callserial}
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     We will not even consider the files that have gone extinct.                       #
   #---------------------------------------------------------------------------------------#
   jobname="${desc}-${polyname}"
   running=$(cat ${jobstat} | grep ${jobname} 2> /dev/null | wc -l)
   case ${runt} in
   "THE_END")
      echo "Polygon has reached the end.  No need to re-submit it."
      submit_now=false
      ;;
   "STSTATE")
      echo "Polygon has reached steady state.  No need to re-submit it."
      submit_now=false
      ;;
   "EXTINCT")
      echo "Polygon population has gone extinct.  No need to re-submit it."
      submit_now=false
      ;;
   "CRASHED"|"METMISS"|"SIGSEGV"|"BAD_MET"|"STOPPED")
      echo "Polygon has serious errors.  Script will not submit the job this time."
      submit_now=false
      ;;
   "INITIAL"|"HISTORY")
      #------------------------------------------------------------------------------------#
      #     Check whether the job is still running
      #------------------------------------------------------------------------------------#
      case "${running}" in
      ""|"0")
         case "${submit}" in
         n|N)
            echo "Polygon will be prepared and ready to be submitted."
            submit_now=false
            ;;
         *)
            echo "Submit polygon:"
            submit_now=true
            ;;
         esac
         ;;
      *)
         echo "Polygon is running.  Do not submit this time."
         submit_now=false
         ;;
      esac
      #------------------------------------------------------------------------------------#
      ;;
   esac
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Submit the job using the specific comments.                                       #
   #---------------------------------------------------------------------------------------#
   if ${submit_now}
   then
      #----- Job options. -----------------------------------------------------------------#
      options="walltime=${runtime},mem=${memory_max}mb,ncpus=${n_cpt}"
      pbsout="${here}/${polyname}/serial_pbs.out"
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Submit, then check whether it went through.  If not, keep trying until         #
      # it works (or give up after nsubtry_max attempts).                                  #
      #------------------------------------------------------------------------------------#
      nfail=1
      attempt=0
      while [[ ${nfail} -gt 0 ]] && [[ ${attempt} -lt ${nsubtry_max} ]]
      do
         let attempt=${attempt}+1

         #----- Submit job. ---------------------------------------------------------------#
         echo -n "  + Attempt number: ${attempt}..."
         qsub -q ${global_queue} -o ${pbsout} -j oe -N ${jobname} -l ${options}            \
                 ${callserial} 1> /dev/null 2> /dev/null
         #---------------------------------------------------------------------------------#


         #----- Wait a bit, then check whether the submission went through. ---------------#
         sleep 3
         nfail=$(qclean | wc -l)
         #---------------------------------------------------------------------------------#



         #------ Check result. ------------------------------------------------------------#
         if [[ ${nfail} -gt 0 ]] && [[ ${attempt} -eq ${nsubtry_max} ]]
         then
            echo "  Failed.  Giving up, check for errors in your script."
         elif [ ${nfail} -eq 0 ]
         then
            echo "  Success."
         else
            echo "  Failed."
         fi
         #---------------------------------------------------------------------------------#
      done
      #------------------------------------------------------------------------------------#
   fi
   #---------------------------------------------------------------------------------------#
done
#------------------------------------------------------------------------------------------#

#----- Delete temporary job status. -------------------------------------------------------#
/bin/rm -f ${jobstat}
#------------------------------------------------------------------------------------------#
