#!/bin/bash
here=`pwd`                            # ! Main path
diskthere='/n/scratch2/moorcroft_lab' # ! Disk where the output files are
thisqueue='moorcroft'                 # ! Queue where jobs should be submitted
lonlat=${here}'/joborder.txt'         # ! File with the job instructions
#----- Outroot is the main output directory. ----------------------------------------------#
outroot='/n/moorcroftfs1/mlongo/diary/XXXXXXXXXXXXX/figures/Template'
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
   isoilbc=`echo ${oi}      | awk '{print $60}'`
   imetrad=`echo ${oi}      | awk '{print $61}'`
   ibranch=`echo ${oi}      | awk '{print $62}'`
   icanrad=`echo ${oi}      | awk '{print $63}'`
   crown=`echo   ${oi}      | awk '{print $64}'`
   ltransvis=`echo ${oi}    | awk '{print $65}'`
   lreflectvis=`echo ${oi}  | awk '{print $66}'`
   ltransnir=`echo ${oi}    | awk '{print $67}'`
   lreflectnir=`echo ${oi}  | awk '{print $68}'`
   orienttree=`echo ${oi}   | awk '{print $69}'`
   orientgrass=`echo ${oi}  | awk '{print $70}'`
   clumptree=`echo ${oi}    | awk '{print $71}'`
   clumpgrass=`echo ${oi}   | awk '{print $72}'`
   ivegtdyn=`echo ${oi}     | awk '{print $73}'`
   igndvap=`echo ${oi}      | awk '{print $74}'`
   iphen=`echo ${oi}        | awk '{print $75}'`
   iallom=`echo ${oi}       | awk '{print $76}'`
   #---------------------------------------------------------------------------------------#



   #----- Find time and minute. -----------------------------------------------------------#
   houra=`echo ${timea}  | awk '{print substr($1,1,2)}'`
   minua=`echo ${timea}  | awk '{print substr($1,3,2)}'`
   hourz=`echo ${timez}  | awk '{print substr($1,1,2)}'`
   minuz=`echo ${timez}  | awk '{print substr($1,3,2)}'`

   #----- Retrieve some information from ED2IN. -------------------------------------------#
   iphysiol=`grep -i NL%IPHYSIOL ${here}/${polyname}/ED2IN | awk '{print $3}'`
   iallom=`grep -i NL%IALLOM ${here}/${polyname}/ED2IN | awk '{print $3}'`
   
   if [ ${droughtmark} == "TRUE" ]
   then 
      let yeara=${droughtyeara}-1
      let yearz=${droughtyearz}+1
   fi

   for script in ${rscripts}
   do
      echo "Submitting script ${script} for polygon: ${polyname}..."

      case ${script} in
      plot_monthly.r)
         let thisyeara=${yeara}+0
         thisdatea=${datea}
         epostout='pmon_epost.out'
         epostlsf='pmon_epost.lsf'
         epostjob='eb-pmon-'${polyiata}
         ;;
      plot_budget.r)
         thisyeara=${yeara}
         let thisdatea=${datea}+1
         epostout='pbdg_epost.out'
         epostlsf='pbdg_epost.lsf'
         epostjob='eb-prk4-'${polyiata}
         ;;
      plot_rk4.r)
         thisyeara=${yeara}
         let thisdatea=${datea}+1
         epostout='prk4_epost.out'
         epostlsf='prk4_epost.lsf'
         epostjob='eb-prk4-'${polyiata}
         ;;
      plot_rk4pc.r)
         thisyeara=${yeara}
         let thisdatea=${datea}+1
         epostout='prpc_epost.out'
         epostlsf='prpc_epost.lsf'
         epostjob='eb-prpc-'${polyiata}
         ;;
      plot_photo.r)
         thisyeara=${yeara}
         let thisdatea=${datea}+1
         epostout='ppht_epost.out'
         epostlsf='ppht_epost.lsf'
         epostjob='eb-ppht-'${polyiata}
         ;;
      patchprops.r)
         thisyeara=${yeara}
         thisdatea=${datea}
         epostout='ppro_epost.out'
         epostlsf='ppro_epost.lsf'
         epostjob='eb-ppro-'${polyiata}
         ;;
      plot_daily.r)
         thisyeara=${yeara}
         thisdatea=${datea}
         epostout='pday_epost.out'
         epostlsf='pday_epost.lsf'
         epostjob='eb-pday-'${polyiata}
         ;;
      plot_fast.r)
         thisyeara=${yeara}
         thisdatea=${datea}
         epostout='pfst_epost.out'
         epostlsf='pfst_epost.lsf'
         epostjob='eb-pfst-'${polyiata}
         ;;
      reject_ed.r)
         thisyeara=${yeara}
         thisdatea=${datea}
         epostout='prej_epost.out'
         epostlsf='prej_epost.lsf'
         epostjob='eb-prej-'${polyiata}
         ;;
      *)
         thisyeara=${yeara}
         thisdatea=${datea}
         epostout='pidn_epost.out'
         epostlsf='pidn_epost.lsf'
         epostjob='eb-pidn-'${polyiata}
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
      sed -i s@thismontha@${montha}@g             ${here}/${polyname}/${script}
      sed -i s@thisdatea@${thisdatea}@g           ${here}/${polyname}/${script}
      sed -i s@thishoura@${houra}@g               ${here}/${polyname}/${script}
      sed -i s@thisminua@${minua}@g               ${here}/${polyname}/${script}
      sed -i s@thisyearz@${yearz}@g               ${here}/${polyname}/${script}
      sed -i s@thismonthz@${monthz}@g             ${here}/${polyname}/${script}
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
      bsub -q ${thisqueue} -J ${epostjob} -o ${polyname}/${epostlsf} "${comm}" 1> /dev/null 2> /dev/null
   done
done
#------------------------------------------------------------------------------------------#
