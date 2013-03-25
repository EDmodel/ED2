#!/bin/bash
. ${HOME}/.bashrc
here='xxxxxxxxxxxxxxxxxxxxx'                  # ! Main path
myself=`whoami`                               # ! You
diskthere=''                                  # ! Disk where the output files are
thisqueue='qqqqqqqqqq'                        # ! Queue where jobs should be submitted
lonlat=${here}'/joborder.txt'                 # ! File with the job instructions
#----- Outroot is the main output directory. ----------------------------------------------#
outroot='xxxxxxxxxxxxxxxxxxxx'
submit='y'       # y = Submit the script; n = Copy the script
#----- Plot only one meteorological cycle. ------------------------------------------------#
useperiod='a'    # Which bounds should I use? (Ignored by plot_eval_ed.r)
                 # 'a' -- All period
                 # 't' -- One eddy flux tower met cycle
                 # 'u' -- User defined period, defined by the variables below.
                 # 'f' -- Force the tower cycle.  You may need to edit the script, though
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
idbhtype=3                     # Type of DBH class
                               # 1 -- Every 10 cm until 100cm; > 100cm
                               # 2 -- 0-10; 10-20; 20-35; 35-50; 50-70; > 70 (cm)
                               # 3 -- 0-10; 10-35; 35-55; > 55 (cm)
#----- Default background colour. ---------------------------------------------------------#
background=0                   # 0 -- White
                               # 1 -- Pitch black
                               # 2 -- Dark grey
#----- Trim the year comparison for tower years only? -------------------------------------#
efttrim="FALSE"
#----- Path with R scripts that are useful. -----------------------------------------------#
rscpath="${HOME}/EDBRAMS/R-utils"
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
if [ ${myself} == "mlongo" ]
then
   rscpath="/n/home00/mlongo/util/Rsc"
fi
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#    If this is an openlava run, load the openlava stuff.                                  #
#------------------------------------------------------------------------------------------#
if [ 'x'${openlava} == 'xy' ] || [ 'x'${openlava} == 'xY' ]
then
   . /opt/openlava-2.0/etc/openlava-client.sh
fi
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#     Make sure the paths are set.                                                         #
#------------------------------------------------------------------------------------------#
if [ ${here} == 'xxxxxxxxxxxxxxxxxxxxx' ] || [ ${outroot} == 'xxxxxxxxxxxxxxxxxxxxx' ] ||
   [ ${thisqueue} =='qqqqqqqqqq' ]
then
   echo " here    = ${here}"
   echo " outroot = ${outroot}"
   echo " queue   = ${queue}"
   echo " Set up variables here, outroot, and queue before using read_monthly.sh!!!"
   exit 99
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
   iage=`echo ${oi}         | awk '{print $16}'`
   polyisoil=`echo ${oi}    | awk '{print $17}'`
   polyntext=`echo ${oi}    | awk '{print $18}'`
   polysand=`echo ${oi}     | awk '{print $19}'`
   polyclay=`echo ${oi}     | awk '{print $20}'`
   polydepth=`echo ${oi}    | awk '{print $21}'`
   polysoilbc=`echo ${oi}   | awk '{print $22}'`
   polysldrain=`echo ${oi}  | awk '{print $23}'`
   polycol=`echo ${oi}      | awk '{print $24}'`
   slzres=`echo ${oi}       | awk '{print $25}'`
   queue=`echo ${oi}        | awk '{print $26}'`
   metdriver=`echo ${oi}    | awk '{print $27}'`
   dtlsm=`echo ${oi}        | awk '{print $28}'`
   vmfactc3=`echo ${oi}     | awk '{print $29}'`
   vmfactc4=`echo ${oi}     | awk '{print $30}'`
   mphototrc3=`echo ${oi}   | awk '{print $31}'`
   mphototec3=`echo ${oi}   | awk '{print $32}'`
   mphotoc4=`echo ${oi}     | awk '{print $33}'`
   bphotoblc3=`echo ${oi}   | awk '{print $34}'`
   bphotonlc3=`echo ${oi}   | awk '{print $35}'`
   bphotoc4=`echo ${oi}     | awk '{print $36}'`
   kwgrass=`echo ${oi}      | awk '{print $37}'`
   kwtree=`echo ${oi}       | awk '{print $38}'`
   gammac3=`echo ${oi}      | awk '{print $39}'`
   gammac4=`echo ${oi}      | awk '{print $40}'`
   d0grass=`echo ${oi}      | awk '{print $41}'`
   d0tree=`echo ${oi}       | awk '{print $42}'`
   alphac3=`echo ${oi}      | awk '{print $43}'`
   alphac4=`echo ${oi}      | awk '{print $44}'`
   klowco2=`echo ${oi}      | awk '{print $45}'`
   decomp=`echo ${oi}       | awk '{print $46}'`
   rrffact=`echo ${oi}      | awk '{print $47}'`
   growthresp=`echo ${oi}   | awk '{print $48}'`
   lwidthgrass=`echo ${oi}  | awk '{print $49}'`
   lwidthbltree=`echo ${oi} | awk '{print $50}'`
   lwidthnltree=`echo ${oi} | awk '{print $51}'`
   q10c3=`echo ${oi}        | awk '{print $52}'`
   q10c4=`echo ${oi}        | awk '{print $53}'`
   h2olimit=`echo ${oi}     | awk '{print $54}'`
   imortscheme=`echo ${oi}  | awk '{print $55}'`
   ddmortconst=`echo ${oi}  | awk '{print $56}'`
   isfclyrm=`echo ${oi}     | awk '{print $57}'`
   icanturb=`echo ${oi}     | awk '{print $58}'`
   ubmin=`echo ${oi}        | awk '{print $59}'`
   ugbmin=`echo ${oi}       | awk '{print $60}'`
   ustmin=`echo ${oi}       | awk '{print $61}'`
   gamm=`echo ${oi}         | awk '{print $62}'`
   gamh=`echo ${oi}         | awk '{print $63}'`
   tprandtl=`echo ${oi}     | awk '{print $64}'`
   ribmax=`echo ${oi}       | awk '{print $65}'`
   atmco2=`echo ${oi}       | awk '{print $66}'`
   thcrit=`echo ${oi}       | awk '{print $67}'`
   smfire=`echo ${oi}       | awk '{print $68}'`
   ifire=`echo ${oi}        | awk '{print $69}'`
   fireparm=`echo ${oi}     | awk '{print $70}'`
   ipercol=`echo ${oi}      | awk '{print $71}'`
   runoff=`echo ${oi}       | awk '{print $72}'`
   imetrad=`echo ${oi}      | awk '{print $73}'`
   ibranch=`echo ${oi}      | awk '{print $74}'`
   icanrad=`echo ${oi}      | awk '{print $75}'`
   crown=`echo   ${oi}      | awk '{print $76}'`
   ltransvis=`echo ${oi}    | awk '{print $77}'`
   lreflectvis=`echo ${oi}  | awk '{print $78}'`
   ltransnir=`echo ${oi}    | awk '{print $79}'`
   lreflectnir=`echo ${oi}  | awk '{print $80}'`
   orienttree=`echo ${oi}   | awk '{print $81}'`
   orientgrass=`echo ${oi}  | awk '{print $82}'`
   clumptree=`echo ${oi}    | awk '{print $83}'`
   clumpgrass=`echo ${oi}   | awk '{print $84}'`
   ivegtdyn=`echo ${oi}     | awk '{print $85}'`
   igndvap=`echo ${oi}      | awk '{print $86}'`
   iphen=`echo ${oi}        | awk '{print $87}'`
   iallom=`echo ${oi}       | awk '{print $88}'`
   ibigleaf=`echo ${oi}     | awk '{print $89}'`
   irepro=`echo ${oi}       | awk '{print $90}'`
   treefall=`echo ${oi}     | awk '{print $91}'`
   #---------------------------------------------------------------------------------------#


   #----- Find time and minute. -----------------------------------------------------------#
   houra=`echo ${timea}  | awk '{print substr($1,1,2)}'`
   minua=`echo ${timea}  | awk '{print substr($1,3,2)}'`
   hourz=`echo ${timez}  | awk '{print substr($1,1,2)}'`
   minuz=`echo ${timez}  | awk '{print substr($1,3,2)}'`
   #---------------------------------------------------------------------------------------#


   #----- Retrieve some information from ED2IN. -------------------------------------------#
   iphysiol=`grep -i NL%IPHYSIOL     ${here}/${polyname}/ED2IN | awk '{print $3}'`
   iallom=`grep   -i NL%IALLOM       ${here}/${polyname}/ED2IN | awk '{print $3}'`
   metcyca=`grep  -i NL%METCYC1      ${here}/${polyname}/ED2IN | awk '{print $3}'`
   metcycz=`grep  -i NL%METCYCF      ${here}/${polyname}/ED2IN | awk '{print $3}'`
   klight=`grep   -i NL%DDMORT_CONST ${here}/${polyname}/ED2IN | awk '{print $3}'`
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
      eftyearz=2009
      ;;
   cax)
      eftyeara=1999
      eftyearz=2003
      ;;
   m34)
      eftyeara=1999
      eftyearz=2005
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
      if [ ${useperiod} == 't' ]
      then
         #------ One meteorological cycle.  Check the type of meteorological driver. ------#
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
         #---------------------------------------------------------------------------------#

      elif [ ${useperiod} == 'u' ]
      then
         #----- The user said which period to use. ----------------------------------------#
         thisyeara=${yusera}
         thisyearz=${yuserz}
         #---------------------------------------------------------------------------------#

      elif [ ${useperiod} == 'f' ]
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
      fi # end [ ${useperiod} == 't' ]
      #------------------------------------------------------------------------------------#



      #----- Set up months and days. ------------------------------------------------------#
      thismontha=${montha}
      thismonthz=${monthz}
      thisdatea=${datea}
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #      Define the job name, and the names of the output files.                       #
      #------------------------------------------------------------------------------------#
      epostout='rmon_epost.out'
      epostsh='rmon_epost.sh'
      epostlsf='rmon_epost.lsf'
      epostjob='eb-rmon-'${polyname}
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Check the status of the run.                                                   #
      #------------------------------------------------------------------------------------#
      statrun=${here}/${polyname}/statusrun.txt
      if [ -s ${statrun} ]
      then
         runt=`cat ${statrun} | awk '{print $6}'`
      else
         runt='INITIAL'
      fi
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #      We submit only the jobs that haven't finished.  If the job has just finished, #
      # we submit once again, but save a file to remember that this polygon is loaded.     #
      #------------------------------------------------------------------------------------#
      status="${here}/${polyname}/rdata_month/status_${polyname}.txt"
      if [ -s ${status} ]
      then
         yearl=`cat ${status} | awk '{print $1}'`
         monthl=`cat ${status} | awk '{print $2}'`
         if [ ${yearl} -eq ${yearz} ] && [ ${monthl} -eq ${monthz} ]
         then
            cestfini="y"
         else
            cestfini="n"
         fi
      else
         cestfini="n"
      fi
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
         inqueue=`bjobs -w -q ${thisqueue} -J ${epostjob} 2> /dev/null | wc -l`
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

         #----- Copy the R script from the Template folder to the local path. -------------#
         cp -f ${here}/Template/read_monthly.r ${here}/${polyname}
         scriptnow=${here}/${polyname}/read_monthly.r
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
         #      plot_eval_ed won't run all at once due to the sheer number of HDF5 files.  #
         # Run it several times until it is complete.                                      #
         #---------------------------------------------------------------------------------#
         echo '#!/bin/bash' > ${here}/${polyname}/${epostsh}
         echo ${comm} >> ${here}/${polyname}/${epostsh}
         chmod +x ${here}/${polyname}/${epostsh}
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Submit the job according to the style (LSF or openlava).                    #
         #---------------------------------------------------------------------------------#
         if [ 'x'${submit} == 'xy' ] || [ 'x'${submit} == 'xY' ]
         then
            #------ Check whether to use openlava or LSF. ---------------------------------#
            if [ 'x'${openlava} == 'xy' ] || [ 'x'${openlava} == 'xY' ]
            then
               bsub="iobsub -J ${epostjob} -o ${here}/${polyname}/${epostlsf}"
               bsub="${bsub} ${here}/${polyname}/${epostsh} 1> /dev/null 2> /dev/null"
            else
               bsub="bsub -q ${thisqueue} -J ${epostjob} -o ${here}/${polyname}/${epostlsf}"
               bsub="${bsub} ${here}/${polyname}/${epostsh} 1> /dev/null 2> /dev/null"
            fi
            #------------------------------------------------------------------------------#
            ${bsub}
         fi
         #---------------------------------------------------------------------------------#
      fi
      #------------------------------------------------------------------------------------#
   else
      echo "${ff} - ${polyname}: directory not found."
   fi
   #---------------------------------------------------------------------------------------#
done
#------------------------------------------------------------------------------------------#

/bin/rm -f ${here}/read_monthly.lock
