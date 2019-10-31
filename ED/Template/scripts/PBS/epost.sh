#!/bin/bash
. ~/.bashrc
here=$(pwd)                           # ! Main path
myself=$(whoami)                      # ! You
joborder="${here}/joborder.txt"       # ! File with the job instructions
#----- POV-Ray include files. -------------------------------------------------------------#
pov_incs="${HOME}/Util/Modules/povray/3.6/share/povray-3.6/include"
#----- Outroot is the main output directory. ----------------------------------------------#
submit="y"       # y -- Submit the script to the queue
                 # l -- Run the script locally
                 # n -- Copy the script
#----- Maximum number of attempts before giving up. ---------------------------------------#
nsubtry_max=5
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
                               # 2 -- 0-10; 10-20; 20-35; 35-55; 55-80; > 80 (cm)
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
#----- bashrc (usually ${HOME}/.bashrc). --------------------------------------------------#
initrc="${HOME}/.bashrc"
#----- Settings for this group of polygons. -----------------------------------------------#
global_queue=""               # Queue
sim_memory=0                  # Memory per simulation.  If zero, then it will be 
                              #    automatically determined by the maximum number of tasks
                              #    per node.
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
rscripts="plot_yearly.r"
#rscripts="yearly_ascii.r"
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


#------------------------------------------------------------------------------------------#
#       First check that the main path and e-mail have been set.  If not, don't run.       #
#------------------------------------------------------------------------------------------#
if [ "x${here}" == "x" ] || [ "x${global_queue}" == "x" ] || [ "x${rscript}" == "x" ]
then
   echo " You must set some variables before running the script:"
   echo " Check variables \"here\", \"global_queue\" and \"rscript\"!"
   exit 99
fi
#------------------------------------------------------------------------------------------#


#----- Load settings. ---------------------------------------------------------------------#
if [ -s ${initrc} ]
then
   . ${initrc}
fi
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#    Use the general path.                                                                 #
#------------------------------------------------------------------------------------------#
case ${myself} in
mlongo)
   rscpath="${HOME}/util/Rsc"
   ;;
marcos.longo)
   rscpath="${SCRATCH}/Util/Rsc"
   rlibs="${SCRATCH}/Util/Rlibs"
   rsync -Prutv ${R_SCRP}/* ${rscpath}
   ;;
marcosl)
   rscpath="${HOME}/Util/Rsc"
   ;;
esac
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Set general limits.  Currently only the maximum time is restricted.                  #
#------------------------------------------------------------------------------------------#
case ${ordinateur} in
sun-master|cmm*)
   #----- SunHPC (LNCC). ------------------------------------------------------------------#
   n_nodes_max=32
   n_cpt=1
   n_tpn=16
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
   n_cpt=1
   n_tpn=16
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
   let n_tasks_max=${n_nodes_max}*${n_tpn}
fi
;;
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#    Use the general path.                                                                 #
#------------------------------------------------------------------------------------------#
case ${myself} in
mlongo|marcosl)
   rscpath="${HOME}/Util/Rsc"
   ;;
marcos.longo)
   rscpath="${SCRATCH}/Util/Rsc"
   rsync -Prutv ${R_SCRP}/* ${rscpath}
   ;;
esac
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#      Define the job name, and the names of the output files.                             #
#------------------------------------------------------------------------------------------#
case ${rscript} in
read_monthly.r)
   epostkey="rmon"
   ;;
yearly_ascii.r)
   epostkey="yasc"
   ;;
r10_monthly.r)
   epostkey="rm10"
   ;;
plot_monthly.r)
   epostkey="pmon"
   ;;
plot_yearly.r)
   epostkey="pyrs"
   ;;
plot_ycomp.r)
   epostkey="pycp"
   ;;
plot_census.r)
   epostkey="pcen"
   ;;
plot_povray.r)
   epostkey="ppov"
   ;;
plot_eval_ed.r)
   epostkey="peed"
   ;;
plot_budget.r)
   epostkey="pbdg"
   ;;
plot_rk4.r)
   epostkey="prk4"
   ;;
plot_rk4pc.r)
   epostkey="prpc"
   ;;
plot_photo.r)
   epostkey="ppht"
   ;;
reject_ed.r)
   epostkey="prej"
   ;;
patchprops.r)
  epostkey="ppro"
  ;;
whichrun.r)
  epostkey="pwhr"
  ;;
plot_daily.r)
   epostkey="pday"
   ;;
plot_fast.r)
   epostkey="pfst"
   ;;
*)
   #---------------------------------------------------------------------------------------#
   #     If the script is here, then it could not find the script... And this should never #
   # happen, so interrupt the script!                                                      #
   #---------------------------------------------------------------------------------------#
   echo " Script ${script} is not recognised by epost.sh!"
   exit 1
   #---------------------------------------------------------------------------------------#
   ;;
esac
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#   Make sure memory does not exceed maximum amount that can be requested.                 #
#------------------------------------------------------------------------------------------#
if [ ${sim_memory} -eq 0 ]
then
   let sim_memory=${node_memory}/${n_tpn}
   let node_memory=${n_tpn}*${sim_memory}
elif [ ${sim_memory} -gt ${node_memory} ]
then 
   echo "Simulation memory ${sim_memory} cannot exceed node memory ${node_memory}!"
   exit 99
else
   #------ Set memory and number of CPUs per task. ----------------------------------------#
   let n_tpn_try=${node_memory}/${sim_memory}
   if [ ${n_tpn_try} -le ${n_tpn} ]
   then
      n_tpn=${n_tpn_try}
      let sim_memory=${node_memory}/${n_tpn}
   else
      let node_memory=${n_tpn}*${sim_memory}
   fi
   #---------------------------------------------------------------------------------------#
fi
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
#      Loop over all polygons.                                                             #
#------------------------------------------------------------------------------------------#
ff=0
while [ ${ff} -lt ${npolys} ]
do

   #---------------------------------------------------------------------------------------#
   #    Format count.                                                                      #
   #---------------------------------------------------------------------------------------#
   let ff=${ff}+1
   let line=${ff}+3
   ffout=$(printf ${pfmt} ${ff})
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
   istomata=$(echo ${oi}     | awk '{print $93 }')
   iplastic=$(echo ${oi}     | awk '{print $94 }')
   igndvap=$(echo ${oi}      | awk '{print $95 }')
   iphen=$(echo ${oi}        | awk '{print $96 }')
   iallom=$(echo ${oi}       | awk '{print $97 }')
   ieconomics=$(echo ${oi}   | awk '{print $98 }')
   igrass=$(echo ${oi}       | awk '{print $99 }')
   ibigleaf=$(echo ${oi}     | awk '{print $100}')
   integscheme=$(echo ${oi}  | awk '{print $101}')
   nsubeuler=$(echo ${oi}    | awk '{print $102}')
   irepro=$(echo ${oi}       | awk '{print $103}')
   treefall=$(echo ${oi}     | awk '{print $104}')
   ianthdisturb=$(echo ${oi} | awk '{print $105}')
   ianthdataset=$(echo ${oi} | awk '{print $106}')
   slscale=$(echo ${oi}      | awk '{print $107}')
   slyrfirst=$(echo ${oi}    | awk '{print $108}')
   slnyrs=$(echo ${oi}       | awk '{print $109}')
   bioharv=$(echo ${oi}      | awk '{print $110}')
   skidarea=$(echo ${oi}     | awk '{print $111}')
   skidsmall=$(echo ${oi}    | awk '{print $112}')
   skidlarge=$(echo ${oi}    | awk '{print $113}')
   fellingsmall=$(echo ${oi} | awk '{print $114}')
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
   tb0|tbx)
      eftyeara=2014
      eftyearz=2017
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



   #----- Print a banner. -----------------------------------------------------------------#
   if [ ${rscript} == "plot_census.r" ] && [ ${subcens} -eq 0 ]
   then
      echo "${ffout} - Skipping submission of ${rscript} for polygon: ${polyname}..."
   else
      echo "${ffout} - Copying script ${rscript} to polygon: ${polyname}..."
   fi
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Set up the time and output variables according to the script.                     #
   #---------------------------------------------------------------------------------------#
   case ${rscript} in
   read_monthly.r|yearly_ascii.r|plot_monthly.r|plot_yearly.r|plot_ycomp.r|plot_census.r|plot_povray.r|r10_monthly.r)
      #------------------------------------------------------------------------------------#
      #     Scripts that are based on monthly means.  The set up is the same, the only     #
      # difference is in the output names.                                                 #
      #------------------------------------------------------------------------------------#
      #------ Check which period to use. --------------------------------------------------#
      if [ ${useperiod} == "t" ]
      then
         #------ One meteorological cycle.  Check the type of meteorological driver. ------#
         case ${metdriver} in
         Sheffield|WFDEI*|ERAINT*|MERRA2*|PGMF3*)
            thisyeara=${metcyca}
            thisyearz=${metcycz}
            ;;
         *)
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
            ;;
         esac #  ${metdriver} in
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

      elif [ ${useperiod} == "b" ]
      then
         #----- The user said to use the eddy flux period. --------------------------------#
         thisyeara=${bioyeara}
         thisyearz=${bioyearz}
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
      ;;
   plot_eval_ed.r)
      #------------------------------------------------------------------------------------#
      #     Cheat by changing metcyca and metcycz in case the meteorological driver is     #
      # Petrolina (output variables exist only for 2004, so we don't need to process       #
      # all years).                                                                        #
      #------------------------------------------------------------------------------------#
      if [ ${metdriver} == "Petrolina" ]
      then 
         thismetcyca=2004
         thismetcycz=2004
      else
         thismetcyca=${metcyca}
         thismetcycz=${metcycz}
      fi
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     The period should be equivalent to one meteorological driver period, so we     #
      # compare apples to apples.  The ED2 years don't need to match as long as we pick    #
      # one cycle.                                                                         #
      #------------------------------------------------------------------------------------#
      thisyeara=${thismetcyca}
      thisyearz=${thismetcycz}
      for i in ${shiftiata}
      do
         if [ "x${i}" == "x${polyiata}" ]
         then
            #----- Always use the true met driver to find the cycle shift. ----------------#
            echo "     -> Shifting met cycle"
            let metcycle=${metcycz}-${metcyca}+1
            let deltayr=${shiftcycle}*${metcycle}
            let thisyeara=${thismetcyca}+${deltayr}
            let thisyearz=${thismetcycz}+${deltayr}
            #------------------------------------------------------------------------------#
         fi # end [ ${i} == ${iata} ]
      done #end for i in ${shiftiata}
      #------------------------------------------------------------------------------------#



      #----- Set up months and days. ------------------------------------------------------#
      thismontha=1
      thismonthz=12
      thisdatea=${datea}
      #------------------------------------------------------------------------------------#
      ;;

   plot_budget.r|plot_rk4.r|plot_rk4pc.r|plot_photo.r|reject_ed.r)
      #------------------------------------------------------------------------------------#
      #     Scripts with very high frequency output (dtlsm or shorter).  The first day     #
      # usually has initialisation problems (for example incoming longwave may be zero     #
      # at the first time step), so we normally skip the first day.                        #
      #------------------------------------------------------------------------------------#
      #----- Check whether to use the user choice of year or the default. -----------------#
      if [ ${useperiod} == "u" ]
      then
         thisyeara=${yusera}
         thisyearz=${yuserz}
      else
         thisyeara=${yeara}
         thisyearz=${yearz}
      fi
      #------------------------------------------------------------------------------------#



      #----- Set up months and days. ------------------------------------------------------#
      thismontha=${montha}
      thismonthz=${monthz}
      let thisdatea=${datea}+1
      #------------------------------------------------------------------------------------#
      ;;


   whichrun.r|patchprops.r)
      #------------------------------------------------------------------------------------#
      #     Script with time-independent patch properties.  No need to skip anything.      #
      #------------------------------------------------------------------------------------#
      #----- Check whether to use the user choice of year or the default. -----------------#
      if [ ${useperiod} == "u" ]
      then
         thisyeara=${yusera}
         thisyearz=${yuserz}
      else
         thisyeara=${yeara}
         thisyearz=${yearz}
      fi
      #------------------------------------------------------------------------------------#



      #----- Set up months and days. ------------------------------------------------------#
      thismontha=${montha}
      thismonthz=${monthz}
      thisdatea=${datea}
      #------------------------------------------------------------------------------------#
      ;;
   plot_daily.r)
      #------------------------------------------------------------------------------------#
      #     Script with daily means.  No need to skip anything.                            #
      #------------------------------------------------------------------------------------#
      #----- Check whether to use the user choice of year or the default. -----------------#
      if [ ${useperiod} == "u" ]
      then
         thisyeara=${yusera}
         thisyearz=${yuserz}
      else
         thisyeara=${yeara}
         thisyearz=${yearz}
      fi
      #------------------------------------------------------------------------------------#



      #----- Set up months and days. ------------------------------------------------------#
      thismontha=${montha}
      thismonthz=${monthz}
      thisdatea=${datea}
      #------------------------------------------------------------------------------------#
      ;;

   plot_fast.r)
      #------------------------------------------------------------------------------------#
      #     Script with short-term averages (usually hourly).  No need to skip any-        #
      # thing.                                                                             #
      #------------------------------------------------------------------------------------#
      if [ ${useperiod} == "u" ]
      then
         thisyeara=${yusera}
         thisyearz=${yuserz}
      else
         thisyeara=${yeara}
         thisyearz=${yearz}
      fi
      #------------------------------------------------------------------------------------#



      #----- Set up months and days. ------------------------------------------------------#
      thismontha=${montha}
      thismonthz=${monthz}
      thisdatea=${datea}
      #------------------------------------------------------------------------------------#
      ;;
   esac
   #---------------------------------------------------------------------------------------#



   #----- Copy the R script from the Template folder to the local path. -------------------#
   cp -f ${here}/Template/${rscript} ${here}/${polyname}
   scriptnow="${here}/${polyname}/${rscript}"
   #---------------------------------------------------------------------------------------#



   #----- Switch the keywords by the current settings. ------------------------------------#
   sed -i s@thispoly@${polyname}@g             ${scriptnow}
   sed -i s@thisoutroot@${here}@g              ${scriptnow}
   sed -i s@thispath@${here}@g                 ${scriptnow}
   sed -i s@thatpath@${here}@g                 ${scriptnow}
   sed -i s@thisrscpath@${rscpath}@g           ${scriptnow}
   sed -i s@thispovincs@${pov_incs}@g          ${scriptnow}
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
   #---------------------------------------------------------------------------------------#




   #----- Make sure this is not the census script for a site we don't have census. --------#
   if [ ${rscript} != "plot_census.r" ] || [ ${subcens} -ne 0 ]
   then
      #----- Set script- and site-specific variables. -------------------------------------#
      epostjob="${epostkey}-${desc}-${polyname}"
      epostnow="${here}/${polyname}/${epostkey}_epost.sh"
      epostout="${here}/${polyname}/${epostkey}_epost.out"
      epostrsc="${here}/${polyname}/${rscript}"
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #      plot_eval_ed won't run all at once due to the sheer number of HDF5 files.     #
      # Run it several times until it is complete.                                         #
      #------------------------------------------------------------------------------------#
      epostexe="R CMD BATCH --no-save --no-restore ${epostrsc} ${epostout}"
      complete="${here}/${polyname}/eval_load_complete.txt"
      options="walltime=${runtime},mem=${sim_memory}mb,ncpus=${n_cpt}"
      rm -f ${epostnow}
      echo "#!/bin/bash"                    > ${epostnow}
      echo "#PBS -S /bin/bash"             >> ${epostnow}
      echo "#PBS -q ${global_queue}"       >> ${epostnow}
      echo "#PBS -o ${epostout}"           >> ${epostnow}
      echo "#PBS -N ${epostjob}"           >> ${epostnow}
      echo "#PBS -j oe"                    >> ${epostnow}
      echo "#PBS -r n"                     >> ${epostnow}
      echo "#PBS -l ${options}"            >> ${epostnow}
      echo " "                             >> ${epostnow}
      echo ". ${initrc}"                   >> ${epostnow}


      case ${rscript} in
      plot_eval_ed.r)
         echo "/bin/rm -fr ${complete}"    >> ${epostnow}
         echo "while [ ! -s ${complete} ]" >> ${epostnow}
         echo "do"                         >> ${epostnow}
         echo "   sleep 3"                 >> ${epostnow}
         echo "   ${epostexe}"             >> ${epostnow}
         echo "done"                       >> ${epostnow}
         ;;
      *)
         echo ${epostexe}                  >> ${epostnow}
         ;;
      esac
      chmod +x ${epostnow}
      #------------------------------------------------------------------------------------#



      #----- Decide whether to submit job or not. -----------------------------------------#
      case "${submit}" in
      n|N) 
         echo " + Files are ready for submission"
         ;;
      *)
         #---------------------------------------------------------------------------------#
         #     Submit, then check whether it went through.  If not, keep trying until      #
         # it works (or give up after nsubtry_max attempts).                               #
         #---------------------------------------------------------------------------------#
         nfail=1
         attempt=0
         while [[ ${nfail} -gt 0 ]] && [[ ${attempt} -lt ${nsubtry_max} ]]
         do
            let attempt=${attempt}+1

            #----- Submit job. ------------------------------------------------------------#
            echo -n "  + Attempt number: ${attempt}..."
            qsub -q ${global_queue} -l ${options} ${callserial} 1> /dev/null 2> /dev/null
            #------------------------------------------------------------------------------#


            #----- Wait a bit, then check whether the submission went through. ------------#
            sleep 3
            nfail=$(qclean | wc -l)
            #------------------------------------------------------------------------------#



            #------ Check result. ---------------------------------------------------------#
            if [[ ${nfail} -gt 0 ]] && [[ ${attempt} -eq ${nsubtry_max} ]]
            then
               echo "  Failed.  Giving up, check for errors in your script."
            elif [ ${nfail} -eq 0 ]
            then
               echo "  Success."
            else
               echo "  Failed."
            fi
            #------------------------------------------------------------------------------#
         done
         #---------------------------------------------------------------------------------#

         ;;
      esac
      #------------------------------------------------------------------------------------#
   fi
   #---------------------------------------------------------------------------------------#

done
#------------------------------------------------------------------------------------------#
