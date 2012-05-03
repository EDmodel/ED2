#!/bin/bash
. ${HOME}/.bashrc
here=`pwd`                            # ! Main path
diskthere='/n/scratch2/moorcroft_lab' # ! Disk where the output files are
thisqueue='wofsy'                     # ! Queue where jobs should be submitted
lonlat=${here}'/joborder.txt'         # ! File with the job instructions
#----- Outroot is the main output directory. ----------------------------------------------#
outroot='/n/moorcroftfs2/mlongo/diary/XXXXXXXXXXX/figures/xxx_XXX/xxxxxxxx'
submit='y'       # y = Submit the script; n = Copy the script
#----- Plot only one meteorological cycle. ------------------------------------------------#
onemetcycle='n'  # Plot only one met cycle only (ignored by plot_eval_ed.r/plot_census.r)
shiftiata=''     # Places that we must shift the cycle
shiftcycle=-1    # In case your met driver doesn't match the model simulation
#----- Check whether to use openlava or typical job submission. ---------------------------#
openlava='n'
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
#    List all the R scripts you want to run.                                               #
#   - plot_monthly.r - This creates several plots based on the monthly mean output.        #
#   - plot_rk4.r     - This creates plots from the detailed output for Runge-Kutta.        #
#                      (patch-level only).                                                 #
#   - plot_photo.r   - This creates plots from the detailed output for Farquhar-Leuning.   #
#   - plot_rk4pc.r   - This creates plots from the detailed output for Runge-Kutta.        #
#                      (patch- and cohort-level).                                          #
#   - plot_budget.r  - This creates plots from the detailed budget for Runge-Kutta.        #
#                      (patch-level only).                                                 #
#   - plot_eval_ed.r - This creates plots comparing model with eddy flux observations.     #
#   - plot_census.r  - This creates plots comparing model with biometric data.             #
#                                                                                          #
#   The following scripts should work too, but I haven't tested them.                      #
#   - plot_daily.r   - This creates plots from the daily mean output.                      #
#   - plot_fast.r    - This creates plots from the analysis files.                         #
#   - patchprops.r   - This creates simple plots showing the patch structure.              #
#   - reject_ed.r    - This tracks the number of steps that were rejected, and what caused #
#                      the step to be rejected.                                            #
#------------------------------------------------------------------------------------------#
rscripts="plot_monthly.r"
#rscripts="patchprops.r plot_photo.r"
#rscripts="patchprops.r"

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
   polycol=`echo ${oi}      | awk '{print $18}'`
   slzres=`echo ${oi}       | awk '{print $19}'`
   queue=`echo ${oi}        | awk '{print $20}'`
   metdriver=`echo ${oi}    | awk '{print $21}'`
   dtlsm=`echo ${oi}        | awk '{print $22}'`
   vmfactc3=`echo ${oi}     | awk '{print $23}'`
   vmfactc4=`echo ${oi}     | awk '{print $24}'`
   mphototrc3=`echo ${oi}   | awk '{print $25}'`
   mphototec3=`echo ${oi}   | awk '{print $26}'`
   mphotoc4=`echo ${oi}     | awk '{print $27}'`
   bphotoblc3=`echo ${oi}   | awk '{print $28}'`
   bphotonlc3=`echo ${oi}   | awk '{print $29}'`
   bphotoc4=`echo ${oi}     | awk '{print $30}'`
   kwgrass=`echo ${oi}      | awk '{print $31}'`
   kwtree=`echo ${oi}       | awk '{print $32}'`
   gammac3=`echo ${oi}      | awk '{print $33}'`
   gammac4=`echo ${oi}      | awk '{print $34}'`
   d0grass=`echo ${oi}      | awk '{print $35}'`
   d0tree=`echo ${oi}       | awk '{print $36}'`
   alphac3=`echo ${oi}      | awk '{print $37}'`
   alphac4=`echo ${oi}      | awk '{print $38}'`
   klowco2=`echo ${oi}      | awk '{print $39}'`
   rrffact=`echo ${oi}      | awk '{print $40}'`
   growthresp=`echo ${oi}   | awk '{print $41}'`
   lwidthgrass=`echo ${oi}  | awk '{print $42}'`
   lwidthbltree=`echo ${oi} | awk '{print $43}'`
   lwidthnltree=`echo ${oi} | awk '{print $44}'`
   q10c3=`echo ${oi}        | awk '{print $45}'`
   q10c4=`echo ${oi}        | awk '{print $46}'`
   h2olimit=`echo ${oi}     | awk '{print $47}'`
   isfclyrm=`echo ${oi}     | awk '{print $48}'`
   icanturb=`echo ${oi}     | awk '{print $49}'`
   ubmin=`echo ${oi}        | awk '{print $50}'`
   ugbmin=`echo ${oi}       | awk '{print $51}'`
   ustmin=`echo ${oi}       | awk '{print $52}'`
   gamm=`echo ${oi}         | awk '{print $53}'`
   gamh=`echo ${oi}         | awk '{print $54}'`
   tprandtl=`echo ${oi}     | awk '{print $55}'`
   ribmax=`echo ${oi}       | awk '{print $56}'`
   atmco2=`echo ${oi}       | awk '{print $57}'`
   thcrit=`echo ${oi}       | awk '{print $58}'`
   smfire=`echo ${oi}       | awk '{print $59}'`
   ifire=`echo ${oi}        | awk '{print $60}'`
   fireparm=`echo ${oi}     | awk '{print $61}'`
   ipercol=`echo ${oi}      | awk '{print $62}'`
   isoilbc=`echo ${oi}      | awk '{print $63}'`
   runoff=`echo ${oi}       | awk '{print $64}'`
   imetrad=`echo ${oi}      | awk '{print $65}'`
   ibranch=`echo ${oi}      | awk '{print $66}'`
   icanrad=`echo ${oi}      | awk '{print $67}'`
   crown=`echo   ${oi}      | awk '{print $68}'`
   ltransvis=`echo ${oi}    | awk '{print $69}'`
   lreflectvis=`echo ${oi}  | awk '{print $70}'`
   ltransnir=`echo ${oi}    | awk '{print $71}'`
   lreflectnir=`echo ${oi}  | awk '{print $72}'`
   orienttree=`echo ${oi}   | awk '{print $73}'`
   orientgrass=`echo ${oi}  | awk '{print $74}'`
   clumptree=`echo ${oi}    | awk '{print $75}'`
   clumpgrass=`echo ${oi}   | awk '{print $76}'`
   ivegtdyn=`echo ${oi}     | awk '{print $77}'`
   igndvap=`echo ${oi}      | awk '{print $78}'`
   iphen=`echo ${oi}        | awk '{print $79}'`
   iallom=`echo ${oi}       | awk '{print $80}'`
   ibigleaf=`echo ${oi}     | awk '{print $81}'`
   irepro=`echo ${oi}       | awk '{print $82}'`
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

   if [ ${droughtmark} == "TRUE" ]
   then 
      let yeara=${droughtyeara}-1
      let yearz=${droughtyearz}+1
   fi

   for script in ${rscripts}
   do
      if [ 'x'${submit} == 'xy' ] || [ 'x'${submit} == 'xY' ]
      then
         echo "Submitting script ${script} for polygon: ${polyname}..."
      else
         echo "Copying script ${script} to polygon: ${polyname}..."
      fi

      case ${script} in
      plot_monthly.r)
         if [ ${onemetcycle} == 'y' ]
         then
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
         else 
            let thisyeara=${yeara}+0
            thisyearz=${yearz}
         fi # end [ ${onemetcycle} == 'y' ]
         thismontha=${montha}
         thismonthz=${monthz}
         thisdatea=${datea}
         epostout='pmon_epost.out'
         epostsh='pmon_epost.sh'
         epostlsf='pmon_epost.lsf'
         epostjob='eb-pmon-'${polyname}
         ;;
      plot_eval_ed.r)
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
         thismontha=1
         thismonthz=12
         thisdatea=${datea}
         epostout='peed_epost.out'
         epostsh='peed_epost.sh'
         epostlsf='peed_epost.lsf'
         epostjob='eb-peed-'${polyname}
         ;;
      plot_census.r)
         thismontha=${montha}
         thismontha=${monthz}
         thisyeara=${yeara}
         thisyearz=${yearz}

         thisdatea=${datea}
         epostout='pcen_epost.out'
         epostsh='pcen_epost.sh'
         epostlsf='pcen_epost.lsf'
         epostjob='eb-pcen-'${polyname}
         ;;
      plot_budget.r)
         thisyeara=${yeara}
         thisyearz=${yearz}
         thismontha=${montha}
         thismonthz=${monthz}
         let thisdatea=${datea}+1
         epostout='pbdg_epost.out'
         epostsh='pbdg_epost.sh'
         epostlsf='pbdg_epost.lsf'
         epostjob='eb-pbdg-'${polyname}
         ;;
      plot_rk4.r)
         thisyeara=${yeara}
         thisyearz=${yearz}
         thismontha=${montha}
         thismonthz=${monthz}
         let thisdatea=${datea}+1
         epostout='prk4_epost.out'
         epostsh='prk4_epost.sh'
         epostlsf='prk4_epost.lsf'
         epostjob='eb-prk4-'${polyname}
         ;;
      plot_rk4pc.r)
         thisyeara=${yeara}
         thisyearz=${yearz}
         thismontha=${montha}
         thismonthz=${monthz}
         let thisdatea=${datea}+1
         epostout='prpc_epost.out'
         epostsh='prpc_epost.sh'
         epostlsf='prpc_epost.lsf'
         epostjob='eb-prpc-'${polyname}
         ;;
      plot_photo.r)
         thisyeara=${yeara}
         thisyearz=${yearz}
         thismontha=${montha}
         thismonthz=${monthz}
         let thisdatea=${datea}+1
         epostout='ppht_epost.out'
         epostsh='ppht_epost.sh'
         epostlsf='ppht_epost.lsf'
         epostjob='eb-ppht-'${polyname}
         ;;
      patchprops.r)
         thisyeara=${yeara}
         thisyearz=${yearz}
         thismontha=${montha}
         thismonthz=${monthz}
         thisdatea=${datea}
         epostout='ppro_epost.out'
         epostsh='ppro_epost.sh'
         epostlsf='ppro_epost.lsf'
         epostjob='eb-ppro-'${polyname}
         ;;
      plot_daily.r)
         thisyeara=${yeara}
         thisyearz=${yearz}
         thismontha=${montha}
         thismonthz=${monthz}
         thisdatea=${datea}
         epostout='pday_epost.out'
         epostsh='pday_epost.sh'
         epostlsf='pday_epost.lsf'
         epostjob='eb-pday-'${polyname}
         ;;
      plot_fast.r)
         thisyeara=${yeara}
         thisyearz=${yearz}
         thismontha=${montha}
         thismonthz=${monthz}
         thisdatea=${datea}
         epostout='pfst_epost.out'
         epostsh='pfst_epost.sh'
         epostlsf='pfst_epost.lsf'
         epostjob='eb-pfst-'${polyname}
         ;;
      reject_ed.r)
         thisyeara=${yeara}
         thisyearz=${yearz}
         thismontha=${montha}
         thismonthz=${monthz}
         thisdatea=${datea}
         epostout='prej_epost.out'
         epostsh='prej_epost.sh'
         epostlsf='prej_epost.lsf'
         epostjob='eb-prej-'${polyname}
         ;;
      *)
         thisyeara=${yeara}
         thisyearz=${yearz}
         thismontha=${montha}
         thismonthz=${monthz}
         thisdatea=${datea}
         epostout='pidn_epost.out'
         epostsh='pidn_epost.sh'
         epostlsf='pidn_epost.lsf'
         epostjob='eb-pidn-'${polyname}
         ;;
      esac


      #----- Copy the R script from the Template folder to the local path. ----------------#
      cp -f ${here}/Template/${script} ${here}/${polyname}

      #----- Switch the keywords by the current settings. ---------------------------------#
      sed -i s@thispoly@${polyname}@g             ${here}/${polyname}/${script}
      sed -i s@thisoutroot@${outroot}@g           ${here}/${polyname}/${script}
      sed -i s@thispath@${here}@g                 ${here}/${polyname}/${script}
      sed -i s@thatpath@${there}@g                ${here}/${polyname}/${script}
      sed -i s@thisyeara@${thisyeara}@g           ${here}/${polyname}/${script}
      sed -i s@thismontha@${thismontha}@g         ${here}/${polyname}/${script}
      sed -i s@thisdatea@${thisdatea}@g           ${here}/${polyname}/${script}
      sed -i s@thishoura@${houra}@g               ${here}/${polyname}/${script}
      sed -i s@thisminua@${minua}@g               ${here}/${polyname}/${script}
      sed -i s@thisyearz@${thisyearz}@g           ${here}/${polyname}/${script}
      sed -i s@thismonthz@${thismonthz}@g         ${here}/${polyname}/${script}
      sed -i s@thisdatez@${datez}@g               ${here}/${polyname}/${script}
      sed -i s@thishourz@${hourz}@g               ${here}/${polyname}/${script}
      sed -i s@thisminuz@${minuz}@g               ${here}/${polyname}/${script}
      sed -i s@myphysiol@${iphysiol}@g            ${here}/${polyname}/${script}
      sed -i s@myallom@${iallom}@g                ${here}/${polyname}/${script}
      sed -i s@mydroughtmark@${droughtmark}@g     ${here}/${polyname}/${script}
      sed -i s@mydroughtyeara@${droughtyeara}@g   ${here}/${polyname}/${script}
      sed -i s@mydroughtyearz@${droughtyearz}@g   ${here}/${polyname}/${script}
      sed -i s@mymonthsdrought@${monthsdrought}@g ${here}/${polyname}/${script}

      #----- Run R to get the plots. ------------------------------------------------------#
      comm="R CMD BATCH ${here}/${polyname}/${script} ${here}/${polyname}/${epostout}"

      #------------------------------------------------------------------------------------#
      #      plot_eval_ed won't run all at once due to the sheer number of HDF5 files.     #
      # Run it several times until it is complete.                                         #
      #------------------------------------------------------------------------------------#
      case ${script} in
      plot_eval_ed.r)
         complete="${here}/${polyname}/eval_load_complete.txt"
         echo '#!/bin/bash'                >  ${here}/${polyname}/${epostsh}
         echo "/bin/rm -fr ${complete}"    >> ${here}/${polyname}/${epostsh}
         echo "while [ ! -s ${complete} ]" >> ${here}/${polyname}/${epostsh}
         echo "do"                         >> ${here}/${polyname}/${epostsh}
         echo "   sleep 3"                 >> ${here}/${polyname}/${epostsh}
         echo "   ${comm}"                 >> ${here}/${polyname}/${epostsh}
         echo "done"                       >> ${here}/${polyname}/${epostsh}
         chmod +x ${here}/${polyname}/${epostsh}
         ;;
      *)
         echo '#!/bin/bash' > ${here}/${polyname}/${epostsh}
         echo ${comm} >> ${here}/${polyname}/${epostsh}
         chmod +x ${here}/${polyname}/${epostsh}
         ;;
      esac
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Submit the job according to the style (LSF or openlava).                       #
      #------------------------------------------------------------------------------------#
      if [ 'x'${submit} == 'xy' ] || [ 'x'${submit} == 'xY' ]
      then
         if [ 'x'${openlava} == 'xy' ] || [ 'x'${openlava} == 'xY' ]
         then
            bsub="iobsub -J ${epostjob} -o ${here}/${polyname}/${epostlsf}"
            bsub="${bsub} ${here}/${polyname}/${epostsh} 1> /dev/null 2> /dev/null"
         else
            bsub="bsub -q ${thisqueue} -J ${epostjob} -o ${polyname}/${epostlsf}"
            bsub="${bsub} ${here}/${polyname}/${epostsh} 1> /dev/null 2> /dev/null"
         fi
         ${bsub}
      fi
      #------------------------------------------------------------------------------------#
   done
   #---------------------------------------------------------------------------------------#
done
#------------------------------------------------------------------------------------------#
