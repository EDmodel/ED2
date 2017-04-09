#!/bin/bash
. ${HOME}/.bashrc
here="xxxxxxxxxxxxxxxxxxxxx"          # ! Main path
myself=$(whoami)                      # ! You
diskthere=""                          # ! Disk where the output files are
thisqueue="qqqqqqqqqq"                # ! Queue where jobs should be submitted
joborder="${here}/joborder.txt"       # ! File with the job instructions
#----- Outroot is the main output directory. ----------------------------------------------#
outroot="xxxxxxxxxxxxxxxxxxxx"
submit="y"       # y = Submit the script; n = Copy the script
#----- Plot only one meteorological cycle. ------------------------------------------------#
useperiod="a"    # Which bounds should I use? (Ignored by plot_eval_ed.r)
                 # "a" -- All period
                 # "t" -- One eddy flux tower met cycle
                 # "u" -- User defined period, defined by the variables below.
                 # "f" -- Force the tower cycle.  You may need to edit the script, though
yusera=1972      # First year to use
yuserz=2011      # Last year to use
#----- Yearly comparison . ----------------------------------------------------------------#
seasonmona=1
#----- Census comparison. -----------------------------------------------------------------#
varcycle="FALSE" # Find the average mortality for various cycles (TRUE/FALSE).
#----- Hourly comparison. -----------------------------------------------------------------#
usedistrib="edf" # Which distribution to plot on top of histograms:
                 #   norm -- Normal distribution
                 #   sn   -- Skewed normal distribution      (requires package sn)
                 #   edf  -- Empirical distribution function (function density)
#----- Output format. ---------------------------------------------------------------------#
outform="c(\"pdf\")"           # x11 - On screen (deprecated on shell scripts)
                               # png - Portable Network Graphics
                               # eps - Encapsulated Post Script
                               # pdf - Portable Document Format
#----- DBH classes. -----------------------------------------------------------------------#
idbhtype=3                     # Type of DBH class
                               # 1 -- Every 10 cm until 100cm; > 100cm
                               # 2 -- 0-10; 10-20; 20-35; 35-50; 50-70; > 70 (cm)
                               # 3 -- 0-10; 10-35; 35-55; > 55 (cm)
#----- Force to run again from scratch. ---------------------------------------------------#
irerun=0                       # Options for re-running:
                               # 0 -- never; updates only.
                               # 1 -- re-run only those that have not finished yet
                               # 2 -- re-run everything, including the finished ones.
#----- Default background colour. ---------------------------------------------------------#
background=0                   # 0 -- White
                               # 1 -- Pitch black
                               # 2 -- Dark grey
#----- Trim the year comparison for tower years only? -------------------------------------#
efttrim="FALSE"
#----- Correction factor for respiration. -------------------------------------------------#
correct_gs=1.0                 # Correction factor for growth and storage respiration
#----- Simple = 1 means that the output is going to be simple. ----------------------------#
simple=0                       # 0 -- default
                               # 1 -- simplified output
#----- Path with R scripts that are useful. -----------------------------------------------#
rscpath="${HOME}/EDBRAMS/R-utils"
#----- bashrc (usually ${HOME}/.bashrc). --------------------------------------------------#
initrc="${HOME}/.bashrc"
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#    Tell whether to plot pseudo-drought or not.                                           #
#------------------------------------------------------------------------------------------#
droughtmark="FALSE"         # Should I plot a rectangle to show the drought?
                            #     capital letters only: TRUE means yes, FALSE means no
droughtyeara=1605           # Year that the first drought instance happens (even if it is 
                            #     just the last bit)
droughtyearz=1609           # Year that the last drought instance happens (even if it 
                            #     partial)
monthsdrought="c(12,1,2,3)" # List of months that get drought, if it starts late in the
                            #     year, put the last month first.
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#    Use the general path.                                                                 #
#------------------------------------------------------------------------------------------#
if [ ${myself} == "marcosl" ]
then
   rscpath="/prj/prjidfca/marcosl/Util/Rsc"
fi
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#     Make sure the paths are set.                                                         #
#------------------------------------------------------------------------------------------#
if [ "y${here}"      == "yxxxxxxxxxxxxxxxxxxxxx" ] || [ "y${here}"      == "y" ]
   [ "y${outroot}"   == "yxxxxxxxxxxxxxxxxxxxxx" ] || [ "y${outroot}"   == "y" ]
   [ "y${thisqueue}" == "yqqqqqqqqqq" ]            || [ "y${thisqueue}" == "y" ]
then
   echo " here    = ${here}"
   echo " outroot = ${outroot}"
   echo " queue   = ${queue}"
   echo " Set up variables here, outroot, and queue before using read_monthly.sh!!!"
   exit 99
fi
#------------------------------------------------------------------------------------------#




#----- Check whether run_sitter.sh is still running or not.  If it is, exit. --------------#
lock="${here}/read_monthly.lock"
if [ -s ${lock} ]
then
   exit
else
   echo "I am going to submit post-processors. Lots of them!" > ${lock}
fi
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#    Make sure that the directory there exists, if not, create all parent directories      #
# needed.                                                                                  #
#------------------------------------------------------------------------------------------#
while [ ! -s ${outroot} ]
do
   namecheck=$(basename ${outroot})
   dircheck=$(dirname ${outroot})
   while [ ! -s ${dircheck} ] && [ ${namecheck} != "/" ]
   do
      namecheck=$(basename ${dircheck})
      dircheck=$(dirname ${dircheck})
   done

   if [ ${namecheck} == "/" ]
   then
      echo "Invalid disk for variable outroot:"
      echo " DISK = ${diskhere}"
      exit 58
   elif [ ${namecheck} == "xxxxxxxx" ] || [ ${namecheck} == "xxx_XXX" ] ||
        [ ${namecheck} == "XXXXXXXXXXX" ]
   then
      echo " - Found this directory in your path: ${namecheck} ..."
      echo " - Outroot given: ${outroot} ..."
      echo " - It looks like you forgot to set up your outroot path, check it!"
      exit 92
   else
      echo "Making directory: ${dircheck}/${namecheck}"
      mkdir ${dircheck}/${namecheck}
   fi
done
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
#      Set the correct script (full or simple).                                            #
#------------------------------------------------------------------------------------------#
case ${simple} in
0)
   read_monthly="read_monthly.r"
   rmon="rmon"
   rdata_path="rdata_month"
   ;;
1)
   read_monthly="read_simple.r"
   rmon="rsim"
   rdata_path="rdata_simple"
   ;;
esac
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
   #    Format count.                                                                      #
   #---------------------------------------------------------------------------------------#
   if   [ ${npolys} -ge 10   ] && [ ${npolys} -lt 100   ]
   then
      ffout=$(printf '%2.2i' ${ff})
   elif [ ${npolys} -ge 100  ] && [ ${npolys} -lt 1000  ]
   then
      ffout=$(printf '%2.2i' ${ff})
   elif [ ${npolys} -ge 100  ] && [ ${npolys} -lt 10000 ]
   then
      ffout=$(printf '%2.2i' ${ff})
   else
      ffout=${ff}
   fi
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
   vmfactc3=$(echo ${oi}     | awk '{print $31 }')
   vmfactc4=$(echo ${oi}     | awk '{print $32 }')
   mphototrc3=$(echo ${oi}   | awk '{print $33 }')
   mphototec3=$(echo ${oi}   | awk '{print $34 }')
   mphotoc4=$(echo ${oi}     | awk '{print $35 }')
   bphotoblc3=$(echo ${oi}   | awk '{print $36 }')
   bphotonlc3=$(echo ${oi}   | awk '{print $37 }')
   bphotoc4=$(echo ${oi}     | awk '{print $38 }')
   kwgrass=$(echo ${oi}      | awk '{print $39 }')
   kwtree=$(echo ${oi}       | awk '{print $40 }')
   gammac3=$(echo ${oi}      | awk '{print $41 }')
   gammac4=$(echo ${oi}      | awk '{print $42 }')
   d0grass=$(echo ${oi}      | awk '{print $43 }')
   d0tree=$(echo ${oi}       | awk '{print $44 }')
   alphac3=$(echo ${oi}      | awk '{print $45 }')
   alphac4=$(echo ${oi}      | awk '{print $46 }')
   klowco2=$(echo ${oi}      | awk '{print $47 }')
   decomp=$(echo ${oi}       | awk '{print $48 }')
   rrffact=$(echo ${oi}      | awk '{print $49 }')
   growthresp=$(echo ${oi}   | awk '{print $50 }')
   lwidthgrass=$(echo ${oi}  | awk '{print $51 }')
   lwidthbltree=$(echo ${oi} | awk '{print $52 }')
   lwidthnltree=$(echo ${oi} | awk '{print $53 }')
   q10c3=$(echo ${oi}        | awk '{print $54 }')
   q10c4=$(echo ${oi}        | awk '{print $55 }')
   h2olimit=$(echo ${oi}     | awk '{print $56 }')
   imortscheme=$(echo ${oi}  | awk '{print $57 }')
   ddmortconst=$(echo ${oi}  | awk '{print $58 }')
   cbrscheme=$(echo ${oi}    | awk '{print $59 }')
   isfclyrm=$(echo ${oi}     | awk '{print $60 }')
   icanturb=$(echo ${oi}     | awk '{print $61 }')
   ubmin=$(echo ${oi}        | awk '{print $62 }')
   ugbmin=$(echo ${oi}       | awk '{print $63 }')
   ustmin=$(echo ${oi}       | awk '{print $64 }')
   gamm=$(echo ${oi}         | awk '{print $65 }')
   gamh=$(echo ${oi}         | awk '{print $66 }')
   tprandtl=$(echo ${oi}     | awk '{print $67 }')
   ribmax=$(echo ${oi}       | awk '{print $68 }')
   atmco2=$(echo ${oi}       | awk '{print $69 }')
   thcrit=$(echo ${oi}       | awk '{print $70 }')
   smfire=$(echo ${oi}       | awk '{print $71 }')
   ifire=$(echo ${oi}        | awk '{print $72 }')
   fireparm=$(echo ${oi}     | awk '{print $73 }')
   ipercol=$(echo ${oi}      | awk '{print $74 }')
   runoff=$(echo ${oi}       | awk '{print $75 }')
   imetrad=$(echo ${oi}      | awk '{print $76 }')
   ibranch=$(echo ${oi}      | awk '{print $77 }')
   icanrad=$(echo ${oi}      | awk '{print $78 }')
   ihrzrad=$(echo ${oi}      | awk '{print $79 }')
   crown=$(echo   ${oi}      | awk '{print $80 }')
   ltransvis=$(echo ${oi}    | awk '{print $81 }')
   lreflectvis=$(echo ${oi}  | awk '{print $82 }')
   ltransnir=$(echo ${oi}    | awk '{print $83 }')
   lreflectnir=$(echo ${oi}  | awk '{print $84 }')
   orienttree=$(echo ${oi}   | awk '{print $85 }')
   orientgrass=$(echo ${oi}  | awk '{print $86 }')
   clumptree=$(echo ${oi}    | awk '{print $87 }')
   clumpgrass=$(echo ${oi}   | awk '{print $88 }')
   igoutput=$(echo ${oi}     | awk '{print $89 }')
   ivegtdyn=$(echo ${oi}     | awk '{print $90 }')
   igndvap=$(echo ${oi}      | awk '{print $91 }')
   iphen=$(echo ${oi}        | awk '{print $92 }')
   iallom=$(echo ${oi}       | awk '{print $93 }')
   ibigleaf=$(echo ${oi}     | awk '{print $94 }')
   irepro=$(echo ${oi}       | awk '{print $95 }')
   treefall=$(echo ${oi}     | awk '{print $96 }')
   ianthdisturb=$(echo ${oi} | awk '{print $97 }')
   ianthdataset=$(echo ${oi} | awk '{print $98 }')
   slscale=$(echo ${oi}      | awk '{print $99 }')
   slnyrs=$(echo ${oi}       | awk '{print $100}')
   bioharv=$(echo ${oi}      | awk '{print $101}')
   skidarea=$(echo ${oi}     | awk '{print $102}')
   #---------------------------------------------------------------------------------------#


   #----- Find time and minute. -----------------------------------------------------------#
   houra=$(echo ${timea}  | awk '{print substr($1,1,2)}')
   minua=$(echo ${timea}  | awk '{print substr($1,3,2)}')
   hourz=$(echo ${timez}  | awk '{print substr($1,1,2)}')
   minuz=$(echo ${timez}  | awk '{print substr($1,3,2)}')
   #---------------------------------------------------------------------------------------#


   #----- Find the last output time. ------------------------------------------------------#
   let monthf=${monthz}-1
   if [ ${monthf} == 0 ]
   then
      monthf=12
      let yearf=${yearz}-1
   else
      yearf=${yearz}
   fi
   #---------------------------------------------------------------------------------------#


   #----- Retrieve some information from ED2IN. -------------------------------------------#
   iphysiol=$(grep -i NL%IPHYSIOL     ${here}/${polyname}/ED2IN | awk '{print $3}')
   iallom=$(grep   -i NL%IALLOM       ${here}/${polyname}/ED2IN | awk '{print $3}')
   metcyca=$(grep  -i NL%METCYC1      ${here}/${polyname}/ED2IN | awk '{print $3}')
   metcycz=$(grep  -i NL%METCYCF      ${here}/${polyname}/ED2IN | awk '{print $3}')
   klight=$(grep   -i NL%DDMORT_CONST ${here}/${polyname}/ED2IN | awk '{print $3}')
   #---------------------------------------------------------------------------------------#


   #---- Find the forest inventory cycle. -------------------------------------------------#
   case ${polyiata} in
   gyf|s67)
      biocyca=2004
      biocycz=2009
      subcens=1
      ;;
   s67)
      biocyca=2001
      biocycz=2011
      subcens=1
      ;;
   *)
      biocyca=${metcyca}
      biocycz=${metcycz}
      subcens=0
      ;;
   esac
   #---------------------------------------------------------------------------------------#



   #---- The eddy flux tower cycles. ------------------------------------------------------#
   case ${polyiata} in
   gyf)
      eftyeara=2004
      eftyearz=2012
      ;;
   cax)
      eftyeara=1999
      eftyearz=2003
      ;;
   m34)
      eftyeara=1999
      eftyearz=2006
      ;;
   s67)
      eftyeara=2001
      eftyearz=2010
      ;;
   s77)
      eftyeara=2001
      eftyearz=2005
      ;;
   s83)
      eftyeara=2000
      eftyearz=2003
      ;;
   pnz)
      eftyeara=2004
      eftyearz=2004
      ;;
   ban)
      eftyeara=2004
      eftyearz=2006
      ;;
   rja)
      eftyeara=1999
      eftyearz=2002
      ;;
   fns)
      eftyeara=1999
      eftyearz=2002
      ;;
   bsb)
      eftyeara=2006
      eftyearz=2011
      ;;
   pdg)
      eftyeara=2001
      eftyearz=2003
      ;;
   hvd)
      eftyeara=1992
      eftyearz=2003
      ;;
   *)
      eftyeara=${metcyca}
      eftyearz=${metcycz}
      ;;
   esac
   #---------------------------------------------------------------------------------------#



   #---- Cheat and force the met cycle to be the tower cycle. -----------------------------#
   if [ ${useperiod} == "f" ]
   then
      metcyca=${eftyeara}
      metcycz=${eftyearz}
   fi
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Switch years in case this is a specific drought run.                              #
   #---------------------------------------------------------------------------------------#
   if [ ${droughtmark} == "TRUE" ]
   then 
      let yeara=${droughtyeara}-1
      let yearz=${droughtyearz}+1
   fi
   #---------------------------------------------------------------------------------------#


   if [ -s ${here}/${polyname} ]
   then


      #------ Check which period to use. --------------------------------------------------#
      if [ ${useperiod} == "t" ]
      then
         #------ One meteorological cycle.  Check the type of meteorological driver. ------#
         if [ ${metdriver} != "Sheffield" ]
         then
            thisyeara=${metcyca}
            thisyearz=${metcycz}
            for i in ${shiftiata}
            do
               if [ "x${i}" == "x${polyiata}" ]
               then
                  echo "     -> Shifting met cycle"
                  let metcycle=${metcycz}-${metcyca}+1
                  let deltayr=${shiftcycle}*${metcycle}
                  let thisyeara=${metcyca}+${deltayr}
                  let thisyearz=${metcycz}+${deltayr}
               fi # end [ ${i} == ${iata} ]
            done #end for i in ${shiftiata}
         else
            thisyeara=${metcyca}
            thisyearz=${metcycz}
         fi # end [ ${metdriver} != "Sheffield" ]
         #---------------------------------------------------------------------------------#

      elif [ ${useperiod} == "u" ]
      then
         #----- The user said which period to use. ----------------------------------------#
         thisyeara=${yusera}
         thisyearz=${yuserz}
         #---------------------------------------------------------------------------------#

      elif [ ${useperiod} == "f" ]
      then
         #----- The user said to use the eddy flux period. --------------------------------#
         thisyeara=${eftyeara}
         thisyearz=${eftyearz}
         #---------------------------------------------------------------------------------#

      else
         #----- Grab all years that the simulation is supposed to run. --------------------#
         thisyeara=${yeara}
         thisyearz=${yearz}
         #---------------------------------------------------------------------------------#
      fi # end [ ${useperiod} == "t" ]
      #------------------------------------------------------------------------------------#



      #----- Set up months and days. ------------------------------------------------------#
      thismontha=${montha}
      thismonthz=${monthz}
      thisdatea=${datea}
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #      Define the job name, and the names of the output files.                       #
      #------------------------------------------------------------------------------------#
      epostout="${rmon}_epost.out"
      epostsh="${rmon}_epost.sh"
      epostlsf="${rmon}_epost.lsf"
      epostjob="eb-${rmon}-${polyname}"
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Check the status of the run.                                                   #
      #------------------------------------------------------------------------------------#
      /bin/rm -f ${here}/${polyname}/statusrun.txt
      /bin/rm -f ${here}/${polyname}/whichrun.r
      /bin/cp -f ${here}/Template/whichrun.r ${here}/${polyname}/whichrun.r
      whichrun="${here}/${polyname}/whichrun.r"
      outwhich="${here}/${polyname}/outwhichrun.txt"
      sed -i s@thispoly@${polyname}@g           ${whichrun}
      sed -i s@thisqueue@${queue}@g             ${whichrun}
      sed -i s@pathhere@${here}@g               ${whichrun}
      sed -i s@paththere@${there}@g             ${whichrun}
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
      if [ ${runt} != "INITIAL" ]
      then
         runt="HISTORY"
      fi
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #      We submit only the jobs that haven't finished.  If the job has just finished, #
      # we submit once again, but save a file to remember that this polygon is loaded.     #
      #------------------------------------------------------------------------------------#
      status="${here}/${polyname}/${rdata_path}/status_${polyname}.txt"
      rdata="${here}/${polyname}/${rdata_path}/${polyname}.RData"
      if [ -s ${status} ]
      then
         yearl=$(cat ${status} | awk '{print $1}')
         monthl=$(cat ${status} | awk '{print $2}')
         if [ ${yearl} -eq ${yearf} ] && [ ${monthl} -eq ${monthf} ]
         then
            cestfini="y"
         else
            cestfini="n"
         fi
         lasttime=${monthl}/${yearl}
      else
         lasttime="Never"
         cestfini="n"
      fi
      #------------------------------------------------------------------------------------#




      #----- Print banner. ----------------------------------------------------------------#
      echo " ${fflab} - ${polyname}"
      echo "   - ED-2.2 Status:             ${runt}"
      echo "   - Last time processed:       ${lasttime}.  Finished: ${cestfini}"
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #       Decide whether to change the status to force running again.                  #
      #------------------------------------------------------------------------------------#
      case ${irerun} in
      0)
         byeprev="n"
         ;;
      1)
         byeprev="y"
         ;;
      2)
         cestfini="n"
         byeprev="y"
         ;;
      esac
      #------------------------------------------------------------------------------------#


      if [ ${runt} == "INITIAL" ]
      then
         submitnow="n"
         echo "${ff} - ${polyname} : polygon hasn't started yet"

      elif [ ${runt} == "THE_END" ] && [ ${cestfini} == "y" ]
      then
         #----- Job has ended and all files have been processed. --------------------------#
         submitnow="n"
         #---------------------------------------------------------------------------------#

         echo "${ff} - ${polyname} : polygon is already loaded or queued for the last time"
      else
         #---------------------------------------------------------------------------------#
         #      Job is still running or it has started again...  Remove the blocker and    #
         # re-submit if the post-processor is not queued.                                  #
         #---------------------------------------------------------------------------------#
         #----- Check that the script is not in the queue. --------------------------------#
         inqueue=$(bjobs -w -q ${thisqueue} -J ${epostjob} 2> /dev/null | wc -l)
         if [ ${inqueue} -eq 0 ]
         then
            submitnow="y"
            echo "${ff} - ${polyname}: submit post-processor script."
         else
            submitnow="n"
            echo "${ff} - ${polyname}: post-processor job has already been queued."
         fi
         #---------------------------------------------------------------------------------#
      fi
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Find out whether the job is on the queue.  In case it is not, re-submit.       #
      #------------------------------------------------------------------------------------#
      if [ "x${submitnow}" == "xy" ]
      then
         #----- Check whether to delete the previous post-processing or not. --------------#
         if [ ${byeprev} == "y" ]
         then
            echo "     * Delete previous post-processing..."
            /bin/rm -fr ${status} ${rdata}
         elif [ -s ${rdata} ]
         then
            echo "     * Continuing previous post-processing..."
         else
            echo "     * Starting new post-processing..."
         fi
         #---------------------------------------------------------------------------------#

         #----- Copy the R script from the Template folder to the local path. -------------#
         cp -f ${here}/Template/${read_monthly} ${here}/${polyname}
         scriptnow=${here}/${polyname}/${read_monthly}
         #---------------------------------------------------------------------------------#



         #----- Switch the keywords by the current settings. ------------------------------#
         sed -i s@thispoly@${polyname}@g             ${scriptnow}
         sed -i s@thisoutroot@${outroot}@g           ${scriptnow}
         sed -i s@thispath@${here}@g                 ${scriptnow}
         sed -i s@thatpath@${there}@g                ${scriptnow}
         sed -i s@thisrscpath@${rscpath}@g           ${scriptnow}
         sed -i s@thisyeara@${thisyeara}@g           ${scriptnow}
         sed -i s@thismontha@${thismontha}@g         ${scriptnow}
         sed -i s@thisdatea@${thisdatea}@g           ${scriptnow}
         sed -i s@thishoura@${houra}@g               ${scriptnow}
         sed -i s@thisminua@${minua}@g               ${scriptnow}
         sed -i s@thisyearz@${thisyearz}@g           ${scriptnow}
         sed -i s@thismonthz@${thismonthz}@g         ${scriptnow}
         sed -i s@thisdatez@${datez}@g               ${scriptnow}
         sed -i s@thishourz@${hourz}@g               ${scriptnow}
         sed -i s@thisminuz@${minuz}@g               ${scriptnow}
         sed -i s@thisseasonmona@${seasonmona}@g     ${scriptnow}
         sed -i s@myphysiol@${iphysiol}@g            ${scriptnow}
         sed -i s@myallom@${iallom}@g                ${scriptnow}
         sed -i s@mydroughtmark@${droughtmark}@g     ${scriptnow}
         sed -i s@mydroughtyeara@${droughtyeara}@g   ${scriptnow}
         sed -i s@mydroughtyearz@${droughtyearz}@g   ${scriptnow}
         sed -i s@mymonthsdrought@${monthsdrought}@g ${scriptnow}
         sed -i s@myvarcycle@${varcycle}@g           ${scriptnow}
         sed -i s@thisoutform@${outform}@g           ${scriptnow}
         sed -i s@mydistrib@${usedistrib}@g          ${scriptnow}
         sed -i s@mymetcyca@${metcyca}@g             ${scriptnow}
         sed -i s@mymetcycz@${metcycz}@g             ${scriptnow}
         sed -i s@mybiocyca@${biocyca}@g             ${scriptnow}
         sed -i s@mybiocycz@${biocycz}@g             ${scriptnow}
         sed -i s@myidbhtype@${idbhtype}@g           ${scriptnow}
         sed -i s@mybackground@${background}@g       ${scriptnow}
         sed -i s@mycorrection@${correct_gs}@g       ${scriptnow}
         sed -i s@myklight@${klight}@g               ${scriptnow}
         sed -i s@myefttrim@${efttrim}@g             ${scriptnow}
         sed -i s@myeftyeara@${eftyeara}@g           ${scriptnow}
         sed -i s@myeftyearz@${eftyearz}@g           ${scriptnow}
         #---------------------------------------------------------------------------------#



         #----- Run R to get the plots. ---------------------------------------------------#
         rbin="R CMD BATCH --no-save --no-restore"
         comm="${rbin} ${scriptnow} ${here}/${polyname}/${epostout}"
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Create shell script with job instructions.                                 #
         #---------------------------------------------------------------------------------#
         sbatchout="${here}/${polyname}/${epostlsf}"
         epostnow="${here}/${polyname}/${epostsh}"
         complete="${here}/${polyname}/eval_load_complete.txt"
         rm -f ${epostnow}
         echo "#$ -S /bin/bash"                > ${epostnow}
         echo "#$ -q ${thisqueue}"            >> ${epostnow}
         echo "#$ -o ${sbatchout}"            >> ${epostnow}
         echo "#$ -N ${epostjob}"             >> ${epostnow}
         echo "#$ -j y"                       >> ${epostnow}
         echo "#$ -r n"                       >> ${epostnow}
         echo " "                             >> ${epostnow}
         echo ". ${initrc}"                   >> ${epostnow}
         echo ${comm}                         >> ${epostnow}
         chmod +x ${epostnow}
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Submit the job.                                                             #
         #---------------------------------------------------------------------------------#
         if [ "x${submit}" == "xy" ] || [ "x${submit}" == "xY" ]
         then
            qsub ${epostnow} 1> /dev/null 2> /dev/null
            sleep 3
         fi
         #---------------------------------------------------------------------------------#

      fi
      #------------------------------------------------------------------------------------#
   else
      echo "${fflab} - ${polyname}: directory not found."
   fi
   #---------------------------------------------------------------------------------------#
done
#------------------------------------------------------------------------------------------#

/bin/rm -f ${lock}

