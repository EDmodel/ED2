#!/bin/bash
here=`pwd`                            # ! Main path
diskthere='/n/scratch2/moorcroft_lab' # ! Disk where the output files are
thisqueue='moorcroft'                 # ! Queue where jobs should be submitted
lonlat=${here}'/joborder.txt'         # ! File with the job instructions
#----- Outroot is the main output directory. ----------------------------------------------#
outroot='/n/moorcroftfs1/mlongo/diary/XXXXXXXXXXXXXX/figures/Template'
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
# 
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
   mphotoc3=`echo ${oi}     | awk '{print $25}'`
   mphotoc4=`echo ${oi}     | awk '{print $26}'`
   kwgrass=`echo ${oi}      | awk '{print $27}'`
   kwtree=`echo ${oi}       | awk '{print $28}'`
   gammac3=`echo ${oi}      | awk '{print $29}'`
   gammac4=`echo ${oi}      | awk '{print $30}'`
   d0grass=`echo ${oi}      | awk '{print $31}'`
   d0tree=`echo ${oi}       | awk '{print $32}'`
   alphac3=`echo ${oi}      | awk '{print $33}'`
   alphac4=`echo ${oi}      | awk '{print $34}'`
   klowco2=`echo ${oi}      | awk '{print $35}'`
   rrffact=`echo ${oi}      | awk '{print $36}'`
   growthresp=`echo ${oi}   | awk '{print $37}'`
   lwidthgrass=`echo ${oi}  | awk '{print $38}'`
   lwidthbltree=`echo ${oi} | awk '{print $39}'`
   lwidthnltree=`echo ${oi} | awk '{print $40}'`
   h2olimit=`echo ${oi}     | awk '{print $41}'`
   isfclyrm=`echo ${oi}     | awk '{print $42}'`
   icanturb=`echo ${oi}     | awk '{print $43}'`
   ubmin=`echo ${oi}        | awk '{print $44}'`
   ugbmin=`echo ${oi}       | awk '{print $45}'`
   ustmin=`echo ${oi}       | awk '{print $46}'`
   gamm=`echo ${oi}         | awk '{print $47}'`
   gamh=`echo ${oi}         | awk '{print $48}'`
   tprandtl=`echo ${oi}     | awk '{print $49}'`
   ribmax=`echo ${oi}       | awk '{print $50}'`
   atmco2=`echo ${oi}       | awk '{print $51}'`
   thcrit=`echo ${oi}       | awk '{print $52}'`
   smfire=`echo ${oi}       | awk '{print $53}'`
   isoilbc=`echo ${oi}      | awk '{print $54}'`
   imetrad=`echo ${oi}      | awk '{print $55}'`
   ibranch=`echo ${oi}      | awk '{print $56}'`
   icanrad=`echo ${oi}      | awk '{print $57}'`
   crown=`echo   ${oi}      | awk '{print $58}'`
   ltransvis=`echo ${oi}    | awk '{print $59}'`
   lreflectvis=`echo ${oi}  | awk '{print $60}'`
   ltransnir=`echo ${oi}    | awk '{print $61}'`
   lreflectnir=`echo ${oi}  | awk '{print $62}'`
   orienttree=`echo ${oi}   | awk '{print $63}'`
   orientgrass=`echo ${oi}  | awk '{print $64}'`
   clumptree=`echo ${oi}    | awk '{print $65}'`
   clumpgrass=`echo ${oi}   | awk '{print $66}'`
   ivegtdyn=`echo ${oi}     | awk '{print $67}'`
   igndvap=`echo ${oi}      | awk '{print $68}'`
   iphen=`echo ${oi}        | awk '{print $69}'`
   iallom=`echo ${oi}       | awk '{print $70}'`
   #---------------------------------------------------------------------------------------#



   #----- Find time and minute. -----------------------------------------------------------#
   houra=`echo ${timea}  | awk '{print substr($1,1,2)}'`
   minua=`echo ${timea}  | awk '{print substr($1,3,2)}'`
   hourz=`echo ${timez}  | awk '{print substr($1,1,2)}'`
   minuz=`echo ${timez}  | awk '{print substr($1,3,2)}'`

   #----- Retrieve some information from ED2IN. -------------------------------------------#
   iphysiol=`grep -i NL%IPHYSIOL ${here}/${polyname}/ED2IN | awk '{print $3}'`
   iallom=`grep -i NL%IALLOM ${here}/${polyname}/ED2IN | awk '{print $3}'`

   for script in ${rscripts}
   do
      echo "Submitting script ${script} for polygon: ${polyname}..."

      case ${script} in
      plot_monthly.r)
         let thisyeara=${yeara}+1
         thisdatea=${datea}
         epostout='pmon_epost.out'
         epostlsf='pmon_epost.lsf'
         epostjob='eb-pmon-'${polyiata}
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
      sed -i s@thispoly@${polyname}@g    ${here}/${polyname}/${script}
      sed -i s@thisoutroot@${outroot}@g  ${here}/${polyname}/${script}
      sed -i s@thispath@${here}@g        ${here}/${polyname}/${script}
      sed -i s@thatpath@${there}@g       ${here}/${polyname}/${script}
      sed -i s@thisyeara@${thisyeara}@g  ${here}/${polyname}/${script}
      sed -i s@thismontha@${montha}@g    ${here}/${polyname}/${script}
      sed -i s@thisdatea@${thisdatea}@g  ${here}/${polyname}/${script}
      sed -i s@thishoura@${houra}@g      ${here}/${polyname}/${script}
      sed -i s@thisminua@${minua}@g      ${here}/${polyname}/${script}
      sed -i s@thisyearz@${yearz}@g      ${here}/${polyname}/${script}
      sed -i s@thismonthz@${monthz}@g    ${here}/${polyname}/${script}
      sed -i s@thisdatez@${datez}@g      ${here}/${polyname}/${script}
      sed -i s@thishourz@${hourz}@g      ${here}/${polyname}/${script}
      sed -i s@thisminuz@${minuz}@g      ${here}/${polyname}/${script}
      sed -i s@myphysiol@${iphysiol}@g   ${here}/${polyname}/${script}
      sed -i s@myallom@${iallom}@g       ${here}/${polyname}/${script}

      #----- Run R to get the plots. ------------------------------------------------------#
      comm="R CMD BATCH ${here}/${polyname}/${script} ${here}/${polyname}/${epostout}"
      bsub -q ${thisqueue} -J ${epostjob} -o ${polyname}/${epostlsf} "${comm}" 1> /dev/null 2> /dev/null
   done
done
#------------------------------------------------------------------------------------------#
