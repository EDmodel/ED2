#!/bin/bash
. ${HOME}/.bashrc
here=`pwd`                            # ! Main path
diskthere='/n/moorcroftfs2'           # ! Disk where the output files are
thisqueue='moorcroft_6100b'           # ! Queue where jobs should be submitted
lonlat=${here}'/joborder.txt'         # ! File with the job instructions
#----- Outroot is the main output directory. ----------------------------------------------#
outroot='/n/moorcroftfs2/mlongo/diary/xxxxxxxx/figures/xxx_XXX/XXXXXXXXXXX'
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
#     Which scripts to run.                                                                #
#                                                                                          #
#   - read_monthly.r - This reads the monthly mean files (results can then be used for     #
#                      plot_monthly.r, plot_yearly.r, and others, but it doesn't plot any- #
#                      thing.)                                                             #
#   - plot_monthly.r - This creates several plots based on the monthly mean output.        #
#   - plot_yearly.r  - This creates plots with year time series.                           #
#   - plot_ycomp.r   - This creates yearly comparisons based on the monthly mean output.   #
#   - plot_rk4.r     - This creates plots from the detailed output for Runge-Kutta.        #
#                      (patch-level only).                                                 #
#   - plot_photo.r   - This creates plots from the detailed output for Farquhar-Leuning.   #
#   - plot_rk4pc.r   - This creates plots from the detailed output for Runge-Kutta.        #
#                      (patch- and cohort-level).                                          #
#   - plot_budget.r  - This creates plots from the detailed budget for Runge-Kutta.        #
#                      (patch-level only).                                                 #
#   - plot_eval_ed.r - This creates plots comparing model with eddy flux observations.     #
#   - plot_census.r  - This creates plots comparing model with biometric data.             #
#   - whichrun.r     - This checks the run status.                                         #
#                                                                                          #
#   The following scripts should work too, but I haven't tested them.                      #
#   - plot_daily.r   - This creates plots from the daily mean output.                      #
#   - plot_fast.r    - This creates plots from the analysis files.                         #
#   - patchprops.r   - This creates simple plots showing the patch structure.              #
#   - reject_ed.r    - This tracks the number of steps that were rejected, and what caused #
#                      the step to be rejected.                                            #
#------------------------------------------------------------------------------------------#
rscripts="plot_yearly.r"
#rscripts="plot_monthly.r"
#rscripts="plot_census.r" 
#rscripts="plot_ycomp.r"
#rscripts="plot_eval_ed.r"
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
#      Loop over all polygons.                                                             #
#------------------------------------------------------------------------------------------#
ff=0
while [ ${ff} -lt ${npolys} ]
do
   let ff=${ff}+1
   let line=${ff}+3
   
   fflab="${ff}/${npolys}"

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


   #---------------------------------------------------------------------------------------#
   #     Loop over all scripts.                                                            #
   #---------------------------------------------------------------------------------------#
   for script in ${rscripts}
   do
      #----- Print a banner. --------------------------------------------------------------#
      if [ ${script} == 'plot_census.r' ] && [ ${subcens} -eq 0 ]
      then
         echo "${fflab} - Skipping submission of ${script} for polygon: ${polyname}..."
      elif [ 'x'${submit} == 'xy' ] || [ 'x'${submit} == 'xY' ]
      then
         echo "${fflab} - Submitting script ${script} for polygon: ${polyname}..."
      else
         echo "${fflab} - Copying script ${script} to polygon: ${polyname}..."
      fi
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Set up the time and output variables according to the script.                  #
      #------------------------------------------------------------------------------------#
      case ${script} in
      read_monthly.r|plot_monthly.r|plot_yearly.r|plot_ycomp.r|plot_census.r)
         #---------------------------------------------------------------------------------#
         #     Scripts that are based on monthly means.  The set up is the same, the only  #
         # difference is in the output names.                                              #
         #---------------------------------------------------------------------------------#
         #------ Check which period to use. -----------------------------------------------#
         if [ ${useperiod} == 't' ]
         then
            #------ One meteorological cycle.  Check the type of meteorological driver. ---#
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
            #------------------------------------------------------------------------------#

         elif [ ${useperiod} == 'u' ]
         then
            #----- The user said which period to use. -------------------------------------#
            thisyeara=${yusera}
            thisyearz=${yuserz}
            #------------------------------------------------------------------------------#
         else
            #----- Grab all years that the simulation is supposed to run. -----------------#
            thisyeara=${yeara}
            thisyearz=${yearz}
            #------------------------------------------------------------------------------#
         fi # end [ ${useperiod} == 't' ]
         #---------------------------------------------------------------------------------#



         #----- Set up months and days. ---------------------------------------------------#
         thismontha=${montha}
         thismonthz=${monthz}
         thisdatea=${datea}
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #      Define the job name, and the names of the output files.                    #
         #---------------------------------------------------------------------------------#
         case ${script} in
         read_monthly.r)
            epostout='rmon_epost.out'
            epostsh='rmon_epost.sh'
            epostlsf='rmon_epost.lsf'
            epostjob='eb-rmon-'${polyname}
            ;;
         plot_monthly.r)
            epostout='pmon_epost.out'
            epostsh='pmon_epost.sh'
            epostlsf='pmon_epost.lsf'
            epostjob='eb-pmon-'${polyname}
            ;;
         plot_yearly.r)
            epostout='pyrs_epost.out'
            epostsh='pyrs_epost.sh'
            epostlsf='pyrs_epost.lsf'
            epostjob='eb-pyrs-'${polyname}
            ;;
         plot_ycomp.r)
            epostout='pycp_epost.out'
            epostsh='pycp_epost.sh'
            epostlsf='pycp_epost.lsf'
            epostjob='eb-pycp-'${polyname}
            ;;
         plot_census.r)
            epostout='pcen_epost.out'
            epostsh='pcen_epost.sh'
            epostlsf='pcen_epost.lsf'
            epostjob='eb-pcen-'${polyname}
            ;;
         esac
         #---------------------------------------------------------------------------------#
         ;;


      plot_eval_ed.r)
         #---------------------------------------------------------------------------------#
         #     Cheat by changing metcyca and metcycz in case the meteorological driver is  #
         # Petrolina (output variables exist only for 2004, so we don't need to process    #
         # all years).                                                                     #
         #---------------------------------------------------------------------------------#
         if [ ${metdriver} == "Petrolina" ]
         then 
            thismetcyca=2004
            thismetcycz=2004
         else
            thismetcyca=${metcyca}
            thismetcycz=${metcycz}
         fi
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #     The period should be equivalent to one meteorological driver period, so we  #
         # compare apples to apples.  The ED2 years don't need to match as long as we pick #
         # one cycle.                                                                      #
         #---------------------------------------------------------------------------------#
         thisyeara=${thismetcyca}
         thisyearz=${thismetcycz}
         for i in ${shiftiata}
         do
            if [ 'x'${i} == 'x'${polyiata} ]
            then
               #----- Always use the true met driver to find the cycle shift. -------------#
               echo '     -> Shifting met cycle'
               let metcycle=${metcycz}-${metcyca}+1
               let deltayr=${shiftcycle}*${metcycle}
               let thisyeara=${thismetcyca}+${deltayr}
               let thisyearz=${thismetcycz}+${deltayr}
               #---------------------------------------------------------------------------#
            fi # end [ ${i} == ${iata} ]
         done #end for i in ${shiftiata}
         #---------------------------------------------------------------------------------#



         #----- Set up months and days. ---------------------------------------------------#
         thismontha=1
         thismonthz=12
         thisdatea=${datea}
         #---------------------------------------------------------------------------------#


         #----- Define the job name, and the names of the output files. -------------------#
         epostout='peed_epost.out'
         epostsh='peed_epost.sh'
         epostlsf='peed_epost.lsf'
         epostjob='eb-peed-'${polyname}
         #---------------------------------------------------------------------------------#

         ;;

      plot_budget.r|plot_rk4.r|plot_rk4pc.r|plot_photo.r|reject_ed.r)
         #---------------------------------------------------------------------------------#
         #     Scripts with very high frequency output (dtlsm or shorter).  The first day  #
         # usually has initialisation problems (for example incoming longwave may be zero  #
         # at the first time step), so we normally skip the first day.                     #
         #---------------------------------------------------------------------------------#
         #----- Check whether to use the user choice of year or the default. --------------#
         if [ ${useperiod} == 'u' ]
         then
            thisyeara=${yusera}
            thisyearz=${yuserz}
         else
            thisyeara=${yeara}
            thisyearz=${yearz}
         fi
         #---------------------------------------------------------------------------------#



         #----- Set up months and days. ---------------------------------------------------#
         thismontha=${montha}
         thismonthz=${monthz}
         let thisdatea=${datea}+1
         #---------------------------------------------------------------------------------#



         #----- Define the job name, and the names of the output files. -------------------#
         case ${script} in 
         plot_budget)
            epostout='pbdg_epost.out'
            epostsh='pbdg_epost.sh'
            epostlsf='pbdg_epost.lsf'
            epostjob='eb-pbdg-'${polyname}
            ;;
         plot_rk4.r)
            epostout='prk4_epost.out'
            epostsh='prk4_epost.sh'
            epostlsf='prk4_epost.lsf'
            epostjob='eb-prk4-'${polyname}
            ;;
         plot_rk4pc.r)
            epostout='prpc_epost.out'
            epostsh='prpc_epost.sh'
            epostlsf='prpc_epost.lsf'
            epostjob='eb-prpc-'${polyname}
            ;;
         plot_photo.r)
            epostout='ppht_epost.out'
            epostsh='ppht_epost.sh'
            epostlsf='ppht_epost.lsf'
            epostjob='eb-ppht-'${polyname}
            ;;
         reject_ed.r)
            epostout='prej_epost.out'
            epostsh='prej_epost.sh'
            epostlsf='prej_epost.lsf'
            epostjob='eb-prej-'${polyname}
            ;;
         esac
         #---------------------------------------------------------------------------------#
         ;;


      whichrun.r|patchprops.r)
         #---------------------------------------------------------------------------------#
         #     Script with time-independent patch properties.  No need to skip anything.   #
         #---------------------------------------------------------------------------------#
         #----- Check whether to use the user choice of year or the default. --------------#
         if [ ${useperiod} == 'u' ]
         then
            thisyeara=${yusera}
            thisyearz=${yuserz}
         else
            thisyeara=${yeara}
            thisyearz=${yearz}
         fi
         #---------------------------------------------------------------------------------#



         #----- Set up months and days. ---------------------------------------------------#
         thismontha=${montha}
         thismonthz=${monthz}
         thisdatea=${datea}
         #---------------------------------------------------------------------------------#



         #----- Define the job name, and the names of the output files. -------------------#
         case ${script} in
         patchprops.r)
           epostout='ppro_epost.out'
           epostsh='ppro_epost.sh'
           epostlsf='ppro_epost.lsf'
           epostjob='eb-ppro-'${polyname}
           ;;
         whichrun.r)
           epostout='pwhr_epost.out'
           epostsh='pwhr_epost.sh'
           epostlsf='pwhr_epost.lsf'
           epostjob='eb-pwhr-'${polyname}
           ;;
         esac
         #---------------------------------------------------------------------------------#
         ;;
      plot_daily.r)
         #---------------------------------------------------------------------------------#
         #     Script with daily means.  No need to skip anything.                         #
         #---------------------------------------------------------------------------------#
         #----- Check whether to use the user choice of year or the default. --------------#
         if [ ${useperiod} == 'u' ]
         then
            thisyeara=${yusera}
            thisyearz=${yuserz}
         else
            thisyeara=${yeara}
            thisyearz=${yearz}
         fi
         #---------------------------------------------------------------------------------#



         #----- Set up months and days. ---------------------------------------------------#
         thismontha=${montha}
         thismonthz=${monthz}
         thisdatea=${datea}
         #---------------------------------------------------------------------------------#



         #----- Define the job name, and the names of the output files. -------------------#
         epostout='pday_epost.out'
         epostsh='pday_epost.sh'
         epostlsf='pday_epost.lsf'
         epostjob='eb-pday-'${polyname}
         #---------------------------------------------------------------------------------#
         ;;

      plot_fast.r)
         #---------------------------------------------------------------------------------#
         #     Script with short-term averages (usually hourly).  No need to skip any-     #
         # thing.                                                                          #
         #---------------------------------------------------------------------------------#
         if [ ${useperiod} == 'u' ]
         then
            thisyeara=${yusera}
            thisyearz=${yuserz}
         else
            thisyeara=${yeara}
            thisyearz=${yearz}
         fi
         #---------------------------------------------------------------------------------#



         #----- Set up months and days. ---------------------------------------------------#
         thismontha=${montha}
         thismonthz=${monthz}
         thisdatea=${datea}
         #---------------------------------------------------------------------------------#



         #----- Define the job name, and the names of the output files. -------------------#
         epostout='pfst_epost.out'
         epostsh='pfst_epost.sh'
         epostlsf='pfst_epost.lsf'
         epostjob='eb-pfst-'${polyname}
         #---------------------------------------------------------------------------------#

         ;;
      *)
         #---------------------------------------------------------------------------------#
         #     If the script is here, then it could not find the script... And this should #
         # never happen, crash!                                                            #
         #---------------------------------------------------------------------------------#
         echo " Script ${script} is not recognised by epost.sh!"
         exit 193
         #---------------------------------------------------------------------------------#
         ;;
      esac
      #------------------------------------------------------------------------------------#



      #----- Copy the R script from the Template folder to the local path. ----------------#
      cp -f ${here}/Template/${script} ${here}/${polyname}
      #------------------------------------------------------------------------------------#



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
      sed -i s@thisseasonmona@${seasonmona}@g     ${here}/${polyname}/${script}
      sed -i s@myphysiol@${iphysiol}@g            ${here}/${polyname}/${script}
      sed -i s@myallom@${iallom}@g                ${here}/${polyname}/${script}
      sed -i s@mydroughtmark@${droughtmark}@g     ${here}/${polyname}/${script}
      sed -i s@mydroughtyeara@${droughtyeara}@g   ${here}/${polyname}/${script}
      sed -i s@mydroughtyearz@${droughtyearz}@g   ${here}/${polyname}/${script}
      sed -i s@mymonthsdrought@${monthsdrought}@g ${here}/${polyname}/${script}
      sed -i s@myvarcycle@${varcycle}@g           ${here}/${polyname}/${script}
      sed -i s@thisoutform@${outform}@g           ${here}/${polyname}/${script}
      sed -i s@mydistrib@${usedistrib}@g          ${here}/${polyname}/${script}
      sed -i s@mymetcyca@${metcyca}@g             ${here}/${polyname}/${script}
      sed -i s@mymetcycz@${metcycz}@g             ${here}/${polyname}/${script}
      sed -i s@mybiocyca@${biocyca}@g             ${here}/${polyname}/${script}
      sed -i s@mybiocycz@${biocycz}@g             ${here}/${polyname}/${script}
      sed -i s@myidbhtype@${idbhtype}@g           ${here}/${polyname}/${script}
      #------------------------------------------------------------------------------------#



      #----- Run R to get the plots. ------------------------------------------------------#
      rbin="R CMD BATCH --no-save --no-restore"
      comm="${rbin} ${here}/${polyname}/${script} ${here}/${polyname}/${epostout}"
      #------------------------------------------------------------------------------------#



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



      #----- Make sure this is not the census script for a site we don't have census. -----#
      if [ ${script} == "plot_census.r" ] && [ ${subcens} -eq 0 ]
      then
         submitnow='n'
      else
         submitnow=${submit}
      fi
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Submit the job according to the style (LSF or openlava).                       #
      #------------------------------------------------------------------------------------#
      if [ 'x'${submitnow} == 'xy' ] || [ 'x'${submitnow} == 'xY' ]
      then
         #------ Check whether to use openlava or LSF. ------------------------------------#
         if [ 'x'${openlava} == 'xy' ] || [ 'x'${openlava} == 'xY' ]
         then
            bsub="iobsub -J ${epostjob} -o ${here}/${polyname}/${epostlsf}"
            bsub="${bsub} ${here}/${polyname}/${epostsh} 1> /dev/null 2> /dev/null"
         else
            bsub="bsub -q ${thisqueue} -J ${epostjob} -o ${polyname}/${epostlsf}"
            bsub="${bsub} ${here}/${polyname}/${epostsh} 1> /dev/null 2> /dev/null"
         fi
         #---------------------------------------------------------------------------------#
         ${bsub}
      fi
      #------------------------------------------------------------------------------------#
   done
   #---------------------------------------------------------------------------------------#
done
#------------------------------------------------------------------------------------------#
