#!/bin/bash
. ${HOME}/.bashrc
here="/nowhere"                       # ! Main path
myself=$(whoami)                      # ! You
diskthere=""                          # ! Disk where the output files are
thisqueue="anyqueue"                  # ! Queue where jobs should be submitted
runtime="7-00:00:00"                  # ! Run time request
memory=2048                           # ! Requested memory (Mb)
sbatch=$(which sbatch)                # ! SLURM command to submit job.
joborder="${here}/joborder.txt"         # ! File with the job instructions
#----- Outroot is the main output directory. ----------------------------------------------#
outroot="/nowhere"
submit="y"       # y = Submit the script; n = Copy the script
#----- Plot only one meteorological cycle. ------------------------------------------------#
useperiod="a"    # Which bounds should I use? (Ignored by plot_eval_ed.r)
                 # "a" -- All period
                 # "t" -- One eddy flux tower met cycle
                 # "u" -- User defined period, defined by the variables below.
                 # "f" -- Force the tower cycle.  You may need to edit the script, though
                 # "b" -- Force one biometry cycle.
yusera=1972      # First year to use
yuserz=2011      # Last year to use
#----- Yearly comparison . ----------------------------------------------------------------#
seasonmona=1
#----- Census comparison. -----------------------------------------------------------------#
varcycle="TRUE"  # Find the average mortality for various cycles (TRUE/FALSE).
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
idbhtype=4                     # Type of DBH class
                               # 1 -- Every 10 cm until 100cm; > 100cm
                               # 2 -- 0-10; 10-20; 20-35; 35-50; 50-70; > 70 (cm)
                               # 3 -- 0-10; 10-35; 35-55; > 55 (cm)
                               # 4 -- 0-10; 10-30; 30-50; 50-80; > 80 (cm)
#----- Default background colour. ---------------------------------------------------------#
background=0                   # 0 -- White
                               # 1 -- Pitch black
                               # 2 -- Dark grey
#----- Select integration interval for some photosynthesis-related variables. -------------#
iint_photo=1                   # 0 -- 24h
                               # 1 -- daytime only
#----- Trim the year comparison for tower years only? -------------------------------------#
efttrim="TRUE"
#----- Correction factor for respiration. -------------------------------------------------#
correct_gs=1.0                 # Correction factor for growth and storage respiration
#----- Use only old-growth patches for census comparison? (plot_census.r only). -----------#
oldgrowth="FALSE"
#----- Path with R scripts that are useful. -----------------------------------------------#
rscpath="${HOME}/EDBRAMS/R-utils"
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#     Which scripts to run.                                                                #
#                                                                                          #
#   - read_monthly.r - This reads the monthly mean files (results can then be used for     #
#                      plot_monthly.r, plot_yearly.r, and others, but it doesn't plot any- #
#                      thing.)                                                             #
#   - yearly_ascii.r - This creates three ascii (csv) files with annual averages of        #
#                      various variables.  It doesn't have all possible variables as it is #
#                      intended to simplify the output for learning purposes.              #
#   - plot_monthly.r - This creates several plots based on the monthly mean output.        #
#   - plot_yearly.r  - This creates plots with year time series.                           #
#   - plot_ycomp.r   - This creates yearly comparisons based on the monthly mean output.   #
#   - plot_povray.r  - This creates yearly plots of the polygon using POV-Ray.             #
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
rscripts="nothing"
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
#     Make sure the main path is set.                                                      #
#------------------------------------------------------------------------------------------#
if [ "x${here}" == "x/nowhere" ]
then
   echo "You must set variable \"here\"!"
   exit 91
fi
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#     Make sure the queue is set.                                                          #
#------------------------------------------------------------------------------------------#
if [ "x${queue}" == "xanyqueue" ]
then
   echo "You must set variable \"queue\"!"
   exit 91
fi
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#     Make sure the R script is set.                                                       #
#------------------------------------------------------------------------------------------#
if [ "x${rscript}" == "xrscript" ]
then
   echo "You must set variable \"rscript\"!"
   exit 91
fi
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#    Make sure that the directory there exists, if not, create all parent directories      #
# needed.                                                                                  #
#------------------------------------------------------------------------------------------#
if [ "x${outroot}" == "x" ]
then
   outroot=${here}
elif [ "x${outroot}" == "x/nowhere" ]
then
   echo "You must set variable \"outroot\"!"
   exit 91
elif [ ! -s ${outroot} ]
then
   mkdir -p ${outroot}
fi
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
   oi=$(head -${line} ${joborder} | tail -1)
   polyname=$(echo ${oi}     | awk '{print $1 }')
   polyiata=$(echo ${oi}     | awk '{print $2 }')
   polylon=$(echo ${oi}      | awk '{print $3 }')
   polylat=$(echo ${oi}      | awk '{print $4 }')
   yeara=$(echo ${oi}        | awk '{print $5 }')
   montha=$(echo ${oi}       | awk '{print $6 }')
   datea=$(echo ${oi}        | awk '{print $7 }')
   timea=$(echo ${oi}        | awk '{print $8 }')
   yearz=$(echo ${oi}        | awk '{print $9 }')
   monthz=$(echo ${oi}       | awk '{print $10}')
   datez=$(echo ${oi}        | awk '{print $11}')
   timez=$(echo ${oi}        | awk '{print $12}')
   initmode=$(echo ${oi}     | awk '{print $13}')
   iscenario=$(echo ${oi}    | awk '{print $14}')
   isizepft=$(echo ${oi}     | awk '{print $15}')
   iage=$(echo ${oi}         | awk '{print $16}')
   imaxcohort=$(echo ${oi}   | awk '{print $17}')
   polyisoil=$(echo ${oi}    | awk '{print $18}')
   polyntext=$(echo ${oi}    | awk '{print $19}')
   polysand=$(echo ${oi}     | awk '{print $20}')
   polyclay=$(echo ${oi}     | awk '{print $21}')
   polydepth=$(echo ${oi}    | awk '{print $22}')
   polysoilbc=$(echo ${oi}   | awk '{print $23}')
   polysldrain=$(echo ${oi}  | awk '{print $24}')
   polycol=$(echo ${oi}      | awk '{print $25}')
   slzres=$(echo ${oi}       | awk '{print $26}')
   queue=$(echo ${oi}        | awk '{print $27}')
   metdriver=$(echo ${oi}    | awk '{print $28}')
   dtlsm=$(echo ${oi}        | awk '{print $29}')
   vmfactc3=$(echo ${oi}     | awk '{print $30}')
   vmfactc4=$(echo ${oi}     | awk '{print $31}')
   mphototrc3=$(echo ${oi}   | awk '{print $32}')
   mphototec3=$(echo ${oi}   | awk '{print $33}')
   mphotoc4=$(echo ${oi}     | awk '{print $34}')
   bphotoblc3=$(echo ${oi}   | awk '{print $35}')
   bphotonlc3=$(echo ${oi}   | awk '{print $36}')
   bphotoc4=$(echo ${oi}     | awk '{print $37}')
   kwgrass=$(echo ${oi}      | awk '{print $38}')
   kwtree=$(echo ${oi}       | awk '{print $39}')
   gammac3=$(echo ${oi}      | awk '{print $40}')
   gammac4=$(echo ${oi}      | awk '{print $41}')
   d0grass=$(echo ${oi}      | awk '{print $42}')
   d0tree=$(echo ${oi}       | awk '{print $43}')
   alphac3=$(echo ${oi}      | awk '{print $44}')
   alphac4=$(echo ${oi}      | awk '{print $45}')
   klowco2=$(echo ${oi}      | awk '{print $46}')
   decomp=$(echo ${oi}       | awk '{print $47}')
   rrffact=$(echo ${oi}      | awk '{print $48}')
   growthresp=$(echo ${oi}   | awk '{print $49}')
   lwidthgrass=$(echo ${oi}  | awk '{print $50}')
   lwidthbltree=$(echo ${oi} | awk '{print $51}')
   lwidthnltree=$(echo ${oi} | awk '{print $52}')
   q10c3=$(echo ${oi}        | awk '{print $53}')
   q10c4=$(echo ${oi}        | awk '{print $54}')
   h2olimit=$(echo ${oi}     | awk '{print $55}')
   imortscheme=$(echo ${oi}  | awk '{print $56}')
   ddmortconst=$(echo ${oi}  | awk '{print $57}')
   cbrscheme=$(echo ${oi}    | awk '{print $58}')
   isfclyrm=$(echo ${oi}     | awk '{print $59}')
   icanturb=$(echo ${oi}     | awk '{print $60}')
   ubmin=$(echo ${oi}        | awk '{print $61}')
   ugbmin=$(echo ${oi}       | awk '{print $62}')
   ustmin=$(echo ${oi}       | awk '{print $63}')
   gamm=$(echo ${oi}         | awk '{print $64}')
   gamh=$(echo ${oi}         | awk '{print $65}')
   tprandtl=$(echo ${oi}     | awk '{print $66}')
   ribmax=$(echo ${oi}       | awk '{print $67}')
   atmco2=$(echo ${oi}       | awk '{print $68}')
   thcrit=$(echo ${oi}       | awk '{print $69}')
   smfire=$(echo ${oi}       | awk '{print $70}')
   ifire=$(echo ${oi}        | awk '{print $71}')
   fireparm=$(echo ${oi}     | awk '{print $72}')
   ipercol=$(echo ${oi}      | awk '{print $73}')
   runoff=$(echo ${oi}       | awk '{print $74}')
   imetrad=$(echo ${oi}      | awk '{print $75}')
   ibranch=$(echo ${oi}      | awk '{print $76}')
   icanrad=$(echo ${oi}      | awk '{print $77}')
   ihrzrad=$(echo ${oi}      | awk '{print $78}')
   crown=$(echo   ${oi}      | awk '{print $79}')
   ltransvis=$(echo ${oi}    | awk '{print $80}')
   lreflectvis=$(echo ${oi}  | awk '{print $81}')
   ltransnir=$(echo ${oi}    | awk '{print $82}')
   lreflectnir=$(echo ${oi}  | awk '{print $83}')
   orienttree=$(echo ${oi}   | awk '{print $84}')
   orientgrass=$(echo ${oi}  | awk '{print $85}')
   clumptree=$(echo ${oi}    | awk '{print $86}')
   clumpgrass=$(echo ${oi}   | awk '{print $87}')
   igoutput=$(echo ${oi}     | awk '{print $88}')
   ivegtdyn=$(echo ${oi}     | awk '{print $89}')
   igndvap=$(echo ${oi}      | awk '{print $90}')
   iphen=$(echo ${oi}        | awk '{print $91}')
   iallom=$(echo ${oi}       | awk '{print $92}')
   ibigleaf=$(echo ${oi}     | awk '{print $93}')
   irepro=$(echo ${oi}       | awk '{print $94}')
   treefall=$(echo ${oi}     | awk '{print $95}')
   ianthdisturb=$(echo ${oi} | awk '{print $96}')
   ianthdataset=$(echo ${oi} | awk '{print $97}')
   #---------------------------------------------------------------------------------------#


   #----- Find time and minute. -----------------------------------------------------------#
   houra=$(echo ${timea}  | awk '{print substr($1,1,2)}')
   minua=$(echo ${timea}  | awk '{print substr($1,3,2)}')
   hourz=$(echo ${timez}  | awk '{print substr($1,1,2)}')
   minuz=$(echo ${timez}  | awk '{print substr($1,3,2)}')
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
   pdg)
      eftyeara=2001
      eftyearz=2003
      ;;
   bsb)
      eftyeara=2006
      eftyearz=2011
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



   #---- The eddy flux tower cycles. ------------------------------------------------------#
   case ${polyiata} in
   gyf)
      bioyeara=2004
      bioyearz=2010
      ;;
   s67)
      bioyeara=1999
      bioyearz=2011
      ;;
   *)
      bioyeara=${eftcyca}
      bioyearz=${eftcycz}
      ;;
   esac
   #---------------------------------------------------------------------------------------#



   #---- Cheat and force the met cycle to be the tower cycle. -----------------------------#
   if [ ${useperiod} == "f" ]
   then
      metcyca=${eftyeara}
      metcycz=${eftyearz}
   elif [ ${useperiod} == "b" ]
   then
      metcyca=${bioyeara}
      metcycz=${bioyearz}
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


   #---------------------------------------------------------------------------------------#
   #     Loop over all scripts.                                                            #
   #---------------------------------------------------------------------------------------#
   for script in ${rscripts}
   do
      #----- Print a banner. --------------------------------------------------------------#
      if [ ${script} == "plot_census.r" ] && [ ${subcens} -eq 0 ]
      then
         echo "${fflab} - Skipping submission of ${script} for polygon: ${polyname}..."
      elif [ "x${submit}" == "xy" ] || [ "x${submit}" == "xY" ]
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
      read_monthly.r|yearly_ascii.r|plot_monthly.r|plot_yearly.r|plot_ycomp.r|plot_census.r|plot_povray.r|r10_monthly.r)
         #---------------------------------------------------------------------------------#
         #     Scripts that are based on monthly means.  The set up is the same, the only  #
         # difference is in the output names.                                              #
         #---------------------------------------------------------------------------------#
         #------ Check which period to use. -----------------------------------------------#
         if [ ${useperiod} == "t" ]
         then
            #------ One meteorological cycle.  Check the type of meteorological driver. ---#
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
            #------------------------------------------------------------------------------#

         elif [ ${useperiod} == "u" ]
         then
            #----- The user said which period to use. -------------------------------------#
            thisyeara=${yusera}
            thisyearz=${yuserz}
            #------------------------------------------------------------------------------#

         elif [ ${useperiod} == "f" ]
         then
            #----- The user said to use the eddy flux period. -----------------------------#
            thisyeara=${eftyeara}
            thisyearz=${eftyearz}
            #------------------------------------------------------------------------------#

         elif [ ${useperiod} == "b" ]
         then
            #----- The user said to use the eddy flux period. -----------------------------#
            thisyeara=${bioyeara}
            thisyearz=${bioyearz}
            #------------------------------------------------------------------------------#

         else
            #----- Grab all years that the simulation is supposed to run. -----------------#
            thisyeara=${yeara}
            thisyearz=${yearz}
            #------------------------------------------------------------------------------#
         fi # end [ ${useperiod} == "t" ]
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
            epostout="rmon_epost.out"
            epostsh="rmon_epost.sh"
            epostlsf="rmon_epost.lsf"
            epostjob="eb-rmon-${polyname}"
            ;;
         yearly_ascii.r)
            epostout="yasc_epost.out"
            epostsh="yasc_epost.sh"
            epostlsf="yasc_epost.lsf"
            epostjob="eb-yasc-${polyname}"
            ;;
         r10_monthly.r)
            epostout="rm10_epost.out"
            epostsh="rm10_epost.sh"
            epostlsf="rm10_epost.lsf"
            epostjob="eb-rm10-${polyname}"
            ;;
         plot_monthly.r)
            epostout="pmon_epost.out"
            epostsh="pmon_epost.sh"
            epostlsf="pmon_epost.lsf"
            epostjob="eb-pmon-${polyname}"
            ;;
         plot_yearly.r)
            epostout="pyrs_epost.out"
            epostsh="pyrs_epost.sh"
            epostlsf="pyrs_epost.lsf"
            epostjob="eb-pyrs-${polyname}"
            ;;
         plot_ycomp.r)
            epostout="pycp_epost.out"
            epostsh="pycp_epost.sh"
            epostlsf="pycp_epost.lsf"
            epostjob="eb-pycp-${polyname}"
            ;;
         plot_census.r)
            epostout="pcen_epost.out"
            epostsh="pcen_epost.sh"
            epostlsf="pcen_epost.lsf"
            epostjob="eb-pcen-${polyname}"
            ;;
         plot_povray.r)
            epostout="ppov_epost.out"
            epostsh="ppov_epost.sh"
            epostlsf="ppov_epost.lsf"
            epostjob="eb-ppov-${polyname}"
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
            if [ "x${i}" == "x${polyiata}" ]
            then
               #----- Always use the true met driver to find the cycle shift. -------------#
               echo "     -> Shifting met cycle"
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
         epostout="peed_epost.out"
         epostsh="peed_epost.sh"
         epostlsf="peed_epost.lsf"
         epostjob="eb-peed-${polyname}"
         #---------------------------------------------------------------------------------#

         ;;

      plot_budget.r|plot_rk4.r|plot_rk4pc.r|plot_photo.r|reject_ed.r)
         #---------------------------------------------------------------------------------#
         #     Scripts with very high frequency output (dtlsm or shorter).  The first day  #
         # usually has initialisation problems (for example incoming longwave may be zero  #
         # at the first time step), so we normally skip the first day.                     #
         #---------------------------------------------------------------------------------#
         #----- Check whether to use the user choice of year or the default. --------------#
         if [ ${useperiod} == "u" ]
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
         plot_budget.r)
            epostout="pbdg_epost.out"
            epostsh="pbdg_epost.sh"
            epostlsf="pbdg_epost.lsf"
            epostjob="eb-pbdg-${polyname}"
            ;;
         plot_rk4.r)
            epostout="prk4_epost.out"
            epostsh="prk4_epost.sh"
            epostlsf="prk4_epost.lsf"
            epostjob="eb-prk4-${polyname}"
            ;;
         plot_rk4pc.r)
            epostout="prpc_epost.out"
            epostsh="prpc_epost.sh"
            epostlsf="prpc_epost.lsf"
            epostjob="eb-prpc-${polyname}"
            ;;
         plot_photo.r)
            epostout="ppht_epost.out"
            epostsh="ppht_epost.sh"
            epostlsf="ppht_epost.lsf"
            epostjob="eb-ppht-${polyname}"
            ;;
         reject_ed.r)
            epostout="prej_epost.out"
            epostsh="prej_epost.sh"
            epostlsf="prej_epost.lsf"
            epostjob="eb-prej-${polyname}"
            ;;
         esac
         #---------------------------------------------------------------------------------#
         ;;


      whichrun.r|patchprops.r)
         #---------------------------------------------------------------------------------#
         #     Script with time-independent patch properties.  No need to skip anything.   #
         #---------------------------------------------------------------------------------#
         #----- Check whether to use the user choice of year or the default. --------------#
         if [ ${useperiod} == "u" ]
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
           epostout="ppro_epost.out"
           epostsh="ppro_epost.sh"
           epostlsf="ppro_epost.lsf"
           epostjob="eb-ppro-${polyname}"
           ;;
         whichrun.r)
           epostout="pwhr_epost.out"
           epostsh="pwhr_epost.sh"
           epostlsf="pwhr_epost.lsf"
           epostjob="eb-pwhr-${polyname}"
           ;;
         esac
         #---------------------------------------------------------------------------------#
         ;;
      plot_daily.r)
         #---------------------------------------------------------------------------------#
         #     Script with daily means.  No need to skip anything.                         #
         #---------------------------------------------------------------------------------#
         #----- Check whether to use the user choice of year or the default. --------------#
         if [ ${useperiod} == "u" ]
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
         epostout="pday_epost.out"
         epostsh="pday_epost.sh"
         epostlsf="pday_epost.lsf"
         epostjob="eb-pday-${polyname}"
         #---------------------------------------------------------------------------------#
         ;;

      plot_fast.r)
         #---------------------------------------------------------------------------------#
         #     Script with short-term averages (usually hourly).  No need to skip any-     #
         # thing.                                                                          #
         #---------------------------------------------------------------------------------#
         if [ ${useperiod} == "u" ]
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
         epostout="pfst_epost.out"
         epostsh="pfst_epost.sh"
         epostlsf="pfst_epost.lsf"
         epostjob="eb-pfst-${polyname}"
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
      scriptnow="${here}/${polyname}/${script}"
      #------------------------------------------------------------------------------------#



      #----- Switch the keywords by the current settings. ---------------------------------#
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
      sed -i s@myiintphoto@${iint_photo}@g        ${scriptnow}
      sed -i s@myklight@${klight}@g               ${scriptnow}
      sed -i s@myefttrim@${efttrim}@g             ${scriptnow}
      sed -i s@myoldgrowth@${oldgrowth}@g         ${scriptnow}
      sed -i s@myeftyeara@${eftyeara}@g           ${scriptnow}
      sed -i s@myeftyearz@${eftyearz}@g           ${scriptnow}
      #------------------------------------------------------------------------------------#



      #----- Run R to get the plots. ------------------------------------------------------#
      rbin="R CMD BATCH --no-save --no-restore"
      comm="${rbin} ${scriptnow} ${here}/${polyname}/${epostout}"
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #      plot_eval_ed won't run all at once due to the sheer number of HDF5 files.     #
      # Run it several times until it is complete.                                         #
      #------------------------------------------------------------------------------------#
      case ${script} in
      plot_eval_ed.r)
         complete="${here}/${polyname}/eval_load_complete.txt"
         echo "#!/bin/bash"                >  ${here}/${polyname}/${epostsh}
         echo "/bin/rm -fr ${complete}"    >> ${here}/${polyname}/${epostsh}
         echo "while [ ! -s ${complete} ]" >> ${here}/${polyname}/${epostsh}
         echo "do"                         >> ${here}/${polyname}/${epostsh}
         echo "   sleep 3"                 >> ${here}/${polyname}/${epostsh}
         echo "   ${comm}"                 >> ${here}/${polyname}/${epostsh}
         echo "done"                       >> ${here}/${polyname}/${epostsh}
         chmod +x ${here}/${polyname}/${epostsh}
         ;;
      *)
         echo "#!/bin/bash" > ${here}/${polyname}/${epostsh}
         echo ${comm} >> ${here}/${polyname}/${epostsh}
         chmod +x ${here}/${polyname}/${epostsh}
         ;;
      esac
      #------------------------------------------------------------------------------------#



      #----- Make sure this is not the census script for a site we don't have census. -----#
      if [ ${script} == "plot_census.r" ] && [ ${subcens} -eq 0 ]
      then
         submitnow="n"
      else
         submitnow=${submit}
      fi
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Submit the job.                                                                #
      #------------------------------------------------------------------------------------#
      if [ "x${submitnow}" == "xy" ] || [ "x${submitnow}" == "xY" ]
      then
         sbatchout="${here}/${polyname}/${epostlsf}"
         epostnow="${here}/${polyname}/${epostsh}"
         sbatch -p ${thisqueue} --mem-per-cpu=${memory} -t ${runtime} -J ${epostjob}       \
             -o ${sbatchout} -n 1 --wrap="${epostnow}"
      fi
      #------------------------------------------------------------------------------------#
   done
   #---------------------------------------------------------------------------------------#
done
#------------------------------------------------------------------------------------------#
