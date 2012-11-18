#!/bin/bash
. ${HOME}/.bashrc
here='/xxxxxxxx/xxxxxxxx/xxx_XXX/XXXXXXXXXXX' # ! Main path
diskthere='/n/moorcroftfs2'           # ! Disk where the output files are
thisqueue='moorcroft2c'                 # ! Queue where jobs should be submitted
lonlat=${here}'/joborder.txt'         # ! File with the job instructions
#----- Outroot is the main output directory. ----------------------------------------------#
outroot='/xxxxxxxx/xxxxxxxx/xxx_XXX/XXXXXXXXXXX/figures'
submit='y'       # y = Submit the script; n = Copy the script
#----- Plot only one meteorological cycle. ------------------------------------------------#
useperiod='t'    # Which bounds should I use? (Ignored by plot_eval_ed.r)
                 # 'a' -- All period
                 # 't' -- One eddy flux tower met cycle
                 # 'u' -- User defined period, defined by the variables below.
yusera=1972      # First year to use
yuserz=2011      # Last year to use
#----- Check whether to use openlava or typical job submission. ---------------------------#
openlava='n'
#----- Yearly comparison . ----------------------------------------------------------------#
seasonmona=1
#----- Census comparison. -----------------------------------------------------------------#
varcycle='TRUE'  # Find the average mortality for various cycles (TRUE/FALSE).
#----- Hourly comparison. -----------------------------------------------------------------#
usedistrib='edf' # Which distribution to plot on top of histograms:
                 #   norm -- Normal distribution
                 #   sn   -- Skewed normal distribution      (requires package sn)
                 #   edf  -- Empirical distribution function (function density)
#----- Output format. ---------------------------------------------------------------------#
outform='c("eps","png","pdf")' # x11 - On screen (deprecated on shell scripts)
                               # png - Portable Network Graphics
                               # eps - Encapsulated Post Script
                               # pdf - Portable Document Format
#----- DBH classes. -----------------------------------------------------------------------#
idbhtype=2                     # Type of DBH class
                               # 1 -- Every 10 cm until 100cm; > 100cm
                               # 2 -- 0-10; 10-20; 20-35; 35-50; 50-70; > 70 (cm)
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
#    If this is an openlava run, load the openlava stuff.                                  #
#------------------------------------------------------------------------------------------#
if [ 'x'${openlava} == 'xy' ] || [ 'x'${openlava} == 'xY' ]
then
   . /opt/openlava-2.0/etc/openlava-client.sh
fi
#------------------------------------------------------------------------------------------#




#----- Check whether run_sitter.sh is still running or not.  If it is, exit. --------------#
if [ -s ${here}/read_monthly.lock ]
then
   exit
else
   echo 'I am going to submit post-processors. Lots of them!' > ${here}/read_monthly.lock
fi
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#    Make sure that the directory there exists, if not, create all parent directories      #
# needed.                                                                                  #
#------------------------------------------------------------------------------------------#
while [ ! -s ${outroot} ]
do
   namecheck=`basename ${outroot}`
   dircheck=`dirname ${outroot}`
   while [ ! -s ${dircheck} ] && [ ${namecheck} != '/' ]
   do
      namecheck=`basename ${dircheck}`
      dircheck=`dirname ${dircheck}`
   done

   if [ ${namecheck} == '/' ]
   then
      echo 'Invalid disk for variable outroot:'
      echo ' DISK ='${diskhere}
      exit 58
   elif [ ${namecheck} == 'xxxxxxxx' ] || [ ${namecheck} == 'xxx_XXX' ] ||
        [ ${namecheck} == 'XXXXXXXXXXX' ]
   then
      echo " - Found this directory in your path: ${namecheck} ..."
      echo " - Outroot given: ${outroot} ..."
      echo " - It looks like you forgot to set up your outroot path, check it!"
      exit 92
   else
      echo 'Making directory: '${dircheck}/${namecheck}
      mkdir ${dircheck}/${namecheck}
   fi
done
#------------------------------------------------------------------------------------------#


#----- Find the disk here to create the "there" path. -------------------------------------#
moi=`whoami`
namehere=`basename ${here}`
diskhere=`dirname ${here}`
while [ ${namehere} != ${moi} ]
do
   namehere=`basename ${diskhere}`
   diskhere=`dirname ${diskhere}`
done
if [ 'x'${diskthere} == 'x' ]
then
   there=${here}
else
   there=`echo ${here} | sed s@${diskhere}@${diskthere}@g`
fi
#------------------------------------------------------------------------------------------#


#----- Determine the number of polygons to run. -------------------------------------------#
let npolys=`wc -l ${lonlat} | awk '{print $1 }'`-3
echo 'Number of polygons: '${npolys}'...'
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
   initmode=`echo ${oi}     | awk '{print $13}'`
   iscenario=`echo ${oi}    | awk '{print $14}'`
   isizepft=`echo ${oi}     | awk '{print $15}'`
   polyisoil=`echo ${oi}    | awk '{print $16}'`
   polyntext=`echo ${oi}    | awk '{print $17}'`
   polysand=`echo ${oi}     | awk '{print $18}'`
   polyclay=`echo ${oi}     | awk '{print $19}'`
   polydepth=`echo ${oi}    | awk '{print $20}'`
   polysoilbc=`echo ${oi}   | awk '{print $21}'`
   polysldrain=`echo ${oi}  | awk '{print $22}'`
   polycol=`echo ${oi}      | awk '{print $23}'`
   slzres=`echo ${oi}       | awk '{print $24}'`
   queue=`echo ${oi}        | awk '{print $25}'`
   metdriver=`echo ${oi}    | awk '{print $26}'`
   dtlsm=`echo ${oi}        | awk '{print $27}'`
   vmfactc3=`echo ${oi}     | awk '{print $28}'`
   vmfactc4=`echo ${oi}     | awk '{print $29}'`
   mphototrc3=`echo ${oi}   | awk '{print $30}'`
   mphototec3=`echo ${oi}   | awk '{print $31}'`
   mphotoc4=`echo ${oi}     | awk '{print $32}'`
   bphotoblc3=`echo ${oi}   | awk '{print $33}'`
   bphotonlc3=`echo ${oi}   | awk '{print $34}'`
   bphotoc4=`echo ${oi}     | awk '{print $35}'`
   kwgrass=`echo ${oi}      | awk '{print $36}'`
   kwtree=`echo ${oi}       | awk '{print $37}'`
   gammac3=`echo ${oi}      | awk '{print $38}'`
   gammac4=`echo ${oi}      | awk '{print $39}'`
   d0grass=`echo ${oi}      | awk '{print $40}'`
   d0tree=`echo ${oi}       | awk '{print $41}'`
   alphac3=`echo ${oi}      | awk '{print $42}'`
   alphac4=`echo ${oi}      | awk '{print $43}'`
   klowco2=`echo ${oi}      | awk '{print $44}'`
   decomp=`echo ${oi}       | awk '{print $45}'`
   rrffact=`echo ${oi}      | awk '{print $46}'`
   growthresp=`echo ${oi}   | awk '{print $47}'`
   lwidthgrass=`echo ${oi}  | awk '{print $48}'`
   lwidthbltree=`echo ${oi} | awk '{print $49}'`
   lwidthnltree=`echo ${oi} | awk '{print $50}'`
   q10c3=`echo ${oi}        | awk '{print $51}'`
   q10c4=`echo ${oi}        | awk '{print $52}'`
   h2olimit=`echo ${oi}     | awk '{print $53}'`
   imortscheme=`echo ${oi}  | awk '{print $54}'`
   ddmortconst=`echo ${oi}  | awk '{print $55}'`
   isfclyrm=`echo ${oi}     | awk '{print $56}'`
   icanturb=`echo ${oi}     | awk '{print $57}'`
   ubmin=`echo ${oi}        | awk '{print $58}'`
   ugbmin=`echo ${oi}       | awk '{print $59}'`
   ustmin=`echo ${oi}       | awk '{print $60}'`
   gamm=`echo ${oi}         | awk '{print $61}'`
   gamh=`echo ${oi}         | awk '{print $62}'`
   tprandtl=`echo ${oi}     | awk '{print $63}'`
   ribmax=`echo ${oi}       | awk '{print $64}'`
   atmco2=`echo ${oi}       | awk '{print $65}'`
   thcrit=`echo ${oi}       | awk '{print $66}'`
   smfire=`echo ${oi}       | awk '{print $67}'`
   ifire=`echo ${oi}        | awk '{print $68}'`
   fireparm=`echo ${oi}     | awk '{print $69}'`
   ipercol=`echo ${oi}      | awk '{print $70}'`
   runoff=`echo ${oi}       | awk '{print $71}'`
   imetrad=`echo ${oi}      | awk '{print $72}'`
   ibranch=`echo ${oi}      | awk '{print $73}'`
   icanrad=`echo ${oi}      | awk '{print $74}'`
   crown=`echo   ${oi}      | awk '{print $75}'`
   ltransvis=`echo ${oi}    | awk '{print $76}'`
   lreflectvis=`echo ${oi}  | awk '{print $77}'`
   ltransnir=`echo ${oi}    | awk '{print $78}'`
   lreflectnir=`echo ${oi}  | awk '{print $79}'`
   orienttree=`echo ${oi}   | awk '{print $80}'`
   orientgrass=`echo ${oi}  | awk '{print $81}'`
   clumptree=`echo ${oi}    | awk '{print $82}'`
   clumpgrass=`echo ${oi}   | awk '{print $83}'`
   ivegtdyn=`echo ${oi}     | awk '{print $84}'`
   igndvap=`echo ${oi}      | awk '{print $85}'`
   iphen=`echo ${oi}        | awk '{print $86}'`
   iallom=`echo ${oi}       | awk '{print $87}'`
   ibigleaf=`echo ${oi}     | awk '{print $88}'`
   irepro=`echo ${oi}       | awk '{print $89}'`
   treefall=`echo ${oi}     | awk '{print $90}'`
   #---------------------------------------------------------------------------------------#


   #----- Find time and minute. -----------------------------------------------------------#
   houra=`echo ${timea}  | awk '{print substr($1,1,2)}'`
   minua=`echo ${timea}  | awk '{print substr($1,3,2)}'`
   hourz=`echo ${timez}  | awk '{print substr($1,1,2)}'`
   minuz=`echo ${timez}  | awk '{print substr($1,3,2)}'`
   #---------------------------------------------------------------------------------------#


   #----- Retrieve some information from ED2IN. -------------------------------------------#
   iphysiol=`grep -i NL%IPHYSIOL ${here}/${polyname}/ED2IN | awk '{print $3}'`
   iallom=`grep -i NL%IALLOM ${here}/${polyname}/ED2IN | awk '{print $3}'`
   metcyca=`grep -i NL%METCYC1 ${here}/${polyname}/ED2IN | awk '{print $3}'`
   metcycz=`grep -i NL%METCYCF ${here}/${polyname}/ED2IN | awk '{print $3}'`
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



   #---------------------------------------------------------------------------------------#
   #     Switch years in case this is a specific drought run.                              #
   #---------------------------------------------------------------------------------------#
   if [ ${droughtmark} == "TRUE" ]
   then 
      let yeara=${droughtyeara}-1
      let yearz=${droughtyearz}+1
   fi
   #---------------------------------------------------------------------------------------#


   #------ Check which period to use. -----------------------------------------------------#
   if [ ${useperiod} == 't' ]
   then
      #------ One meteorological cycle.  Check the type of meteorological driver. ---------#
      if [ ${metdriver} != 'Sheffield' ]
      then
         thisyeara=${metcyca}
         thisyearz=${metcycz}
         for i in ${shiftiata}
         do
            if [ 'x'${i} == 'x'${polyiata} ]
            then
               echo '     -> Shifting met cycle'
               let metcycle=${metcycz}-${metcyca}+1
               let deltayr=${shiftcycle}*${metcycle}
               let thisyeara=${metcyca}+${deltayr}
               let thisyearz=${metcycz}+${deltayr}
            fi # end [ ${i} == ${iata} ]
         done #end for i in ${shiftiata}
      else
         thisyeara=${metcyca}
         thisyearz=${metcycz}
      fi # end [ ${metdriver} != 'Sheffield' ]
      #------------------------------------------------------------------------------------#

   elif [ ${useperiod} == 'u' ]
   then
      #----- The user said which period to use. -------------------------------------------#
      thisyeara=${yusera}
      thisyearz=${yuserz}
      #------------------------------------------------------------------------------------#
   else
      #----- Grab all years that the simulation is supposed to run. -----------------------#
      thisyeara=${yeara}
      thisyearz=${yearz}
      #------------------------------------------------------------------------------------#
   fi # end [ ${useperiod} == 't' ]
   #---------------------------------------------------------------------------------------#



   #----- Set up months and days. ---------------------------------------------------------#
   thismontha=${montha}
   thismonthz=${monthz}
   thisdatea=${datea}
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Define the job name, and the names of the output files.                          #
   #---------------------------------------------------------------------------------------#
   epostout='rmon_epost.out'
   epostsh='rmon_epost.sh'
   epostlsf='rmon_epost.lsf'
   epostjob='eb-rmon-'${polyname}
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Check the status of the run.                                                      #
   #---------------------------------------------------------------------------------------#
   statrun=${here}/${polyname}/statusrun.txt
   if [ -s ${statrun} ]
   then
      runt=`cat ${statrun} | awk '{print $6}'`
   else
      runt='INITIAL'
   fi
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      We submit only the jobs that haven't finished.  If the job has just finished, we #
   # submit once again, but save a file to remember that this polygon is loaded.           #
   #---------------------------------------------------------------------------------------#
   fullload="${here}/${polyname}/qfiles_loaded.txt"
   if [ ${runt} == "INITIAL" ]
   then
      submitnow="n"
      echo "${ff} - ${polyname} : polygon hasn't started yet"

   elif [ ${runt} == "THE_END" ] && [ -s ${fullload} ]
   then
      #----- Job has ended and all files have been processed. -----------------------------#
      submitnow="n"
      #------------------------------------------------------------------------------------#

      echo "${ff} - ${polyname} : polygon is already loaded or queued for the last time"

   elif [ ${runt} == "THE_END" ]
   then
      #------------------------------------------------------------------------------------#
      #      Job has ended but loading is not complete.  Run for one last time.            #
      #------------------------------------------------------------------------------------#
      #----- Check that the script is not in the queue. -----------------------------------#
      inqueue=`bjobs -w -q ${thisqueue} -J ${epostjob} 2> /dev/null | wc -l`
      if [ ${inqueue} -eq 0 ]
      then
         #----- Save the time to the file that will block future submissions. -------------#
         when=`date +'%d %B %Y - %R %Z'`
         echo "Last submission on ${when}" > ${fullload}
         #---------------------------------------------------------------------------------#

         submitnow="y"
         echo "${ff} - ${polyname}: run has finished! Submit script for the last time."
      else
         submitnow="n"
         echo "${ff} - ${polyname}: post-processor job has already been queued."
      fi
      #------------------------------------------------------------------------------------#
   else
      #------------------------------------------------------------------------------------#
      #      Job is still running or it has started again...  Remove the blocker and       #
      # re-submit if the post-processor is not queued.                                     #
      #------------------------------------------------------------------------------------#
      #----- Delete the blocker. ----------------------------------------------------------#
      /bin/rm -f ${fullload}
      #----- Check that the script is not in the queue. -----------------------------------#
      inqueue=`bjobs -w -q ${thisqueue} -J ${epostjob} 2> /dev/null | wc -l`
      if [ ${inqueue} -eq 0 ]
      then
         submitnow="y"
         echo "${ff} - ${polyname}: submit post-processor script."
      else
         submitnow="n"
         echo "${ff} - ${polyname}: post-processor job has already been queued."
      fi
      #------------------------------------------------------------------------------------#
   fi
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Find out whether the job is on the queue.  In case it is not, re-submit.          #
   #---------------------------------------------------------------------------------------#
   if [ "x${submitnow}" == "xy" ]
   then

      #----- Copy the R script from the Template folder to the local path. ----------------#
      cp -f ${here}/Template/read_monthly.r ${here}/${polyname}
      scriptnow=${here}/${polyname}/read_monthly.r
      #------------------------------------------------------------------------------------#



      #----- Switch the keywords by the current settings. ---------------------------------#
      sed -i s@thispoly@${polyname}@g             ${scriptnow}
      sed -i s@thisoutroot@${outroot}@g           ${scriptnow}
      sed -i s@thispath@${here}@g                 ${scriptnow}
      sed -i s@thatpath@${there}@g                ${scriptnow}
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
      #------------------------------------------------------------------------------------#



      #----- Run R to get the plots. ------------------------------------------------------#
      rbin="R CMD BATCH --no-save --no-restore"
      comm="${rbin} ${scriptnow} ${here}/${polyname}/${epostout}"
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #      plot_eval_ed won't run all at once due to the sheer number of HDF5 files.     #
      # Run it several times until it is complete.                                         #
      #------------------------------------------------------------------------------------#
      echo '#!/bin/bash' > ${here}/${polyname}/${epostsh}
      echo ${comm} >> ${here}/${polyname}/${epostsh}
      chmod +x ${here}/${polyname}/${epostsh}
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Submit the job according to the style (LSF or openlava).                       #
      #------------------------------------------------------------------------------------#
      if [ 'x'${submit} == 'xy' ] || [ 'x'${submit} == 'xY' ]
      then
         #------ Check whether to use openlava or LSF. ------------------------------------#
         if [ 'x'${openlava} == 'xy' ] || [ 'x'${openlava} == 'xY' ]
         then
            bsub="iobsub -J ${epostjob} -o ${here}/${polyname}/${epostlsf}"
            bsub="${bsub} ${here}/${polyname}/${epostsh} 1> /dev/null 2> /dev/null"
         else
            bsub="bsub -q ${thisqueue} -J ${epostjob} -o ${here}/${polyname}/${epostlsf}"
            bsub="${bsub} ${here}/${polyname}/${epostsh} 1> /dev/null 2> /dev/null"
         fi
         #---------------------------------------------------------------------------------#
         ${bsub}
      fi
      #------------------------------------------------------------------------------------#
   fi
   #---------------------------------------------------------------------------------------#
done
#------------------------------------------------------------------------------------------#

/bin/rm -f ${here}/read_monthly.lock
