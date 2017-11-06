#!/bin/bash
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#     User control variables                                                               #
#------------------------------------------------------------------------------------------#

#------------------------------------------------------------------------------------------#
# VERSION -- This is the unique identifier tag for this test environment.  It should       #
#            indicate the revision number branched from, the users initials and a test     #
#            number.  For git versions, use the first 7 characters instance if you         #
#            branched from 3e31dd3, your initials are xyz, your local branch is b13aac2    #
#            and this is your 1st attempt to run the test suite, you may use:              #
#------------------------------------------------------------------------------------------#
VERSION="3e31dd3-xyz-b13aac2-v1"
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#      These are the executables.  They must exist.  Duh.                                  #
#                                                                                          #
#  EXE_PATH      -- Path where the three executables are located.                          #
#  MAIN_EXE_PATH -- MainLine, the code you downloaded from https://github.com/EDmodel/ED2. #
#                   The out-of-the-box version, no changes in the source code.             #
#  TEST_EXE_PATH -- The version you are willing to push to the MainLine, compiled the same #
#                   way most normal people would compile.                                  #
#  DBUG_EXE_PATH -- The same source code as TEST_EXE, but with every possible debugging    #
#                   flag you can think of.  If your F_OPTS doesn't have at least a dozen   #
#                   options or have no idea what FPE stands for, then you are not          #
#                   debugging it properly.                                                 #
#------------------------------------------------------------------------------------------#
EXE_PATH="$(pwd)/executables"
MAIN_EXE_PATH="${EXE_PATH}/ed_2.1-main"
TEST_EXE_PATH="${EXE_PATH}/ed_2.1-test"
DBUG_EXE_PATH="${EXE_PATH}/ed_2.1-dbug"
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
# TESTTYPE -- Which tests are you going to run?  Options are:                              #
#             "rapid"  -- 2-3 years of simulation                                          #
#             "medium" -- ~ 125 years of simulation.   The long tests will take a long     #
#                         time, beware.                                                    #
#             "long"  -- ~ 300 years of simulation.   The long tests will take a really    #
#                        really long time, beware.                                         #
#------------------------------------------------------------------------------------------#
TESTTYPE="rapid"
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
# RUNTYPE -- Should the test suite start from scratch or continue from where it stopped?   #
#            Options are:                                                                  #
#            INITIAL -- start from the beginning.  Most often the case                     #
#            HISTORY -- continue from where it stopped.  Only use this option in case      #
#                       you want to resume runs that were stopped by reasons other than    #
#                       the model crashing (like power outage, queue terminated job, full  #
#                       disk, etc.)                                                        #
#------------------------------------------------------------------------------------------#
RUNTYPE="HISTORY"
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#  DATAPATH -- Location with all data sets.                                                #
#------------------------------------------------------------------------------------------#
DATAPATH="/n/moorcroftfs5/mlongo/edts_datasets"
DATAPATH="/n/regal/moorcroft_lab/mlongo/data/edts_datasets"
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#     The following flags are switches to decide which sites to run.  You should           #
# ultimately have results for all of them, but this is helpful if some sites give you      #
# problems after the first try.  Y for Yes and N for No (case insensitive).                #
#------------------------------------------------------------------------------------------#
USE_M34="y"     # Manaus K34 POI
USE_S67="y"     # Santarem km 67 POI
USE_HAR="y"     # Harvard Forest POI
USE_PDG="y"     # Pe-de-Gigante POI
USE_TON="y"     # Tonzi POI
USE_CAX="y"     # Caxiuana POI
USE_TNF="y"     # Tapajos National Forest POI
USE_ATA="y"     # Atacama Desert POI
USE_PET="y"     # Petrolina POI
USE_GYF="y"     # Paracou POI
USE_S83="y"     # Santarem km 83 (logging) POI
USE_PRG="y"     # Paragominas (ALS init) POI
USE_TL2="y"     # Toolik (boreal) POI
USE_HIP="y"     # Petrolina High Frequency Detailed Short POI
USE_HIM="y"     # Manaus High Frequency Detailed Short POI
USE_HIH="y"     # Manaus High Frequency Detailed Short POI (Hybrid)
USE_RJG="y"     # Gridded 12x12 simulation centred on Reserva Jaru
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#     The following flags are the queue names for each simulation.                         #
#   Main, test, and dbug should all use the same queue.                                    #
#------------------------------------------------------------------------------------------#
Q_M34="moorcroft_amd"      # Manaus K34 POI
Q_S67="moorcroft_6100"     # Santarem km 67 POI
Q_HAR="moorcroft_amd"      # Harvard Forest POI
Q_PDG="moorcroft_6100"     # Pe-de-Gigante POI
Q_TON="moorcroft_6100"     # Tonzi POI
Q_CAX="moorcroft_6100"     # Caxiuana POI
Q_TNF="moorcroft_amd"      # Tapajos National Forest POI
Q_ATA="moorcroft_6100"     # Atacama Desert POI
Q_PET="moorcroft_6100"     # Petrolina POI
Q_GYF="moorcroft_6100"     # Paracou POI
Q_S83="moorcroft_amd"      # Santarem Km 83 (logging) POI
Q_PRG="moorcroft_amd"      # Paragominas (ALS init) POI
Q_TL2="moorcroft_amd"      # Toolik (boreal) POI
Q_HIP="moorcroft_6100"     # Petrolina High Frequency Detailed Short POI
Q_HIM="moorcroft_6100"     # Manaus High Frequency Detailed Short POI
Q_HIH="moorcroft_6100"     # Manaus High Frequency Detailed Short POI (Hybrid)
Q_RJG="moorcroft_6100"     # Gridded 12x12 simulation centred on Reserva Jaru
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#     The following flags are the memory per cpu for each simulation.  Check the queue     #
# memory and number of CPUs to define these variables.  Also remember that simulation PRG  #
# requires at least 5550 Mb of memory because of the large initialisation.                 #
#------------------------------------------------------------------------------------------#
MEM_M34="1845"     # Manaus K34 POI
MEM_S67="1845"     # Santarem km 67 POI
MEM_HAR="1845"     # Harvard Forest POI
MEM_PDG="1845"     # Pe-de-Gigante POI
MEM_TON="1845"     # Tonzi POI
MEM_CAX="1845"     # Caxiuana POI
MEM_TNF="1845"     # Tapajos National Forest POI
MEM_ATA="1845"     # Atacama Desert POI
MEM_PET="1845"     # Petrolina POI
MEM_GYF="1845"     # Paracou POI
MEM_S83="1845"     # Santarem Km 83 (logging) POI
MEM_PRG="5500"     # Paragominas (ALS init) POI
MEM_TL2="1845"     # Toolik (Boreal) POI
MEM_HIP="1845"     # Petrolina High Frequency Detailed Short POI
MEM_HIM="1845"     # Manaus High Frequency Detailed Short POI
MEM_HIH="1845"     # Manaus High Frequency Detailed Short POI (Hybrid)
MEM_RJG="1845"     # Gridded 12x12 simulation centred on Reserva Jaru
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#     How many cores for runs.  For POI and High-frequency, this sets the number of        #
# threads for each run (using by some parallel do loops).  This will only work if the code #
# is compiled with -openmp (ifort) or -fopenmp (gfortran), otherwise it will be ignored by #
# the code, although it will request the resources.  For gridded runs, this is the number  #
# of processors to split the domain. This requires the code to be compiled with mpif90 and #
# option PAR_DEFS=-DRAMS_MPI set in the include.mk file.                                   #
#                                                                                          #
# High-frequency runs by default require 1 CPU, as they have only one patch.               #
#------------------------------------------------------------------------------------------#
CPU_M34="12"   # Manaus K34 POI (MAXPATCH=20; RAPID_INIT_MODE=5)
CPU_S67="12"   # Santarem km 67 POI (MAXPATCH=20; RAPID_INIT_MODE=0)
CPU_HAR="12"   # Harvard Forest POI (MAXPATCH=20; RAPID_INIT_MODE=6)
CPU_PDG="12"   # Pe-de-Gigante POI (MAXPATCH=20; RAPID_INIT_MODE=0)
CPU_TON="12"   # Tonzi POI (MAXPATCH=8; RAPID_INIT_MODE=5)
CPU_CAX="12"   # Caxiuana POI (MAXPATCH=20; RAPID_INIT_MODE=0)
CPU_TNF="12"   # Tapajos National Forest POI (MAXPATCH=20; RAPID_INIT_MODE=5)
CPU_ATA="12"   # Atacama Desert POI (MAXPATCH=20; RAPID_INIT_MODE=0)
CPU_PET="12"   # Petrolina POI (MAXPATCH=20; RAPID_INIT_MODE=6)
CPU_GYF="12"   # Paracou POI (MAXPATCH=20; RAPID_INIT_MODE=0)
CPU_S83="12"   # Santarem Km 83 (logging) POI (MAXPATCH=20; RAPID_INIT_MODE=5)
CPU_PRG="24"   # Paragominas (ALS init) POI (MAXPATCH=24; RAPID_INIT_MODE=6)
CPU_TL2="12"   # Toolik (Boreal) POI (MAXPATCH=20; RAPID_INIT_MODE=0)
CPU_RJG="27"   # Gridded 12x12 simulation centred on Reserva Jaru (NPTS=144)
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#     The following variables are the time required by each simulation, hh:mm:ss.          #
#------------------------------------------------------------------------------------------#
POI_TIME="Infinite"
GRID_TIME="Infinite"
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#  TEST_DESCRIPTION -- A short but comprehensive explanation of the tests.  Explain what   #
#                      the commits had involved.                                           #
#------------------------------------------------------------------------------------------#
TEST_DESCRIPTION="1. Changed minimum height for reproduction based on BCI measurements (tropical trees only).  2. Changes in minimum reproduction size, so seed_rain works.  3. Changed compilation instructions so it can compile most files with -O3 but uses -O2 for files that would otherwise take days to compile with ifort 13.  4.  New dist_type categories, which now distinguishes tree fall, logging, fires, abandonment and forest plantations.  5.  Changed forestry.f90, this file computes the disturbance rates but lets disturbance.f90 to apply them do cpoly%disturbance_rates is now correctly updated.  6.   Several minor bug fixes."
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
# VERSION_BRANCHED_FROM -- The version you branched from, in case you didn't guess.        #
#                          We know that may had indicated which version you branched from, #
#                          but indicate it here too.                                       #
#------------------------------------------------------------------------------------------#
VERSION_BRANCHED_FROM="3e31dd3"
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
# TESTER_NAME    -- Who is running this test?                                              #
# COMMITTER_NAME -- Who were the developers that actually made the changes to the code     #
#                   that is being tested?                                                  #
#------------------------------------------------------------------------------------------#
TESTER_NAME="John Harvard"
COMMITTER_NAME="Drew Faust"
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
# SUBMIT_JOBS    -- Should the script submit the jobs?                                     #
#                   y or Y -- Generate script then run the script.                         #
#                   n or N -- Generate script but don't submit the jobs.                   #
#------------------------------------------------------------------------------------------#
SUBMIT_JOBS="n"
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
#     Unless you are developing the script, you don't need to change anything beyond this  #
# point.                                                                                   #
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




echo ""
echo ""
echo "==========================================================================="
echo "                  Starting the EDM Dev Test Suite  (EDTS)                  "
echo "==========================================================================="
echo ""


#---- Make test type lower case. ----------------------------------------------------------#
TESTTYPE=$(echo ${TESTTYPE} | tr '[:upper:]' '[:lower:]')
#------------------------------------------------------------------------------------------#



#---- Define some runtime variables for POI. ----------------------------------------------#
declare -a USE_SITE=( ${USE_M34} ${USE_S67} ${USE_HAR} ${USE_PDG} ${USE_TON} ${USE_CAX} \
                      ${USE_TNF} ${USE_ATA} ${USE_PET} ${USE_GYF} ${USE_S83} ${USE_TL2} \
                      ${USE_PRG} )
declare -a SITEID=(m34 s67 har pdg ton cax tnf ata pet gyf s83 prg tl2)
declare -a SITEPFX=(M34 S67 HAR PDG TON CAX TNF ATA PET GYF S83 PRG TL2)
declare -a SITEQ=( ${Q_M34} ${Q_S67} ${Q_HAR} ${Q_PDG} ${Q_TON} ${Q_CAX} \
                   ${Q_TNF} ${Q_ATA} ${Q_PET} ${Q_GYF} ${Q_S83} ${Q_PRG} \
                   ${Q_TL2} )
declare -a SITEMEM=( ${MEM_M34} ${MEM_S67} ${MEM_HAR} ${MEM_PDG} ${MEM_TON} ${MEM_CAX} \
                     ${MEM_TNF} ${MEM_ATA} ${MEM_PET} ${MEM_GYF} ${MEM_S83} ${MEM_PRG} \
                     ${MEM_TL2} )
declare -a SITECPU=( ${CPU_M34} ${CPU_S67} ${CPU_HAR} ${CPU_PDG} ${CPU_TON} ${CPU_CAX} \
                     ${CPU_TNF} ${CPU_ATA} ${CPU_PET} ${CPU_GYF} ${CPU_S83} ${CPU_PRG} \
                     ${CPU_TL2} )
#------------------------------------------------------------------------------------------#



#----- POI debug time. --------------------------------------------------------------------#
declare -a D_IYEARAS=(1500 1500 2007 1500 2000 2000 2002 1500 2005 2007 2000 2011 2008)
declare -a D_IYEARZS=(1501 1501 2008 1501 2001 2001 2002 1501 2006 2008 2001 2012 2009)
declare -a D_INITMDS=(5    0    6    0    5    0    5    0    6    0    5    6    0   )
declare -a D_RUNTYPS=( INITIAL INITIAL INITIAL INITIAL INITIAL INITIAL \
                       INITIAL INITIAL INITIAL INITIAL INITIAL INITIAL INITIAL )
#------------------------------------------------------------------------------------------#



#----- High-frequency detailed runs (HI-PET, HI-M34, and HI-M34/Hybrid). ------------------#
declare -a USE_HIFR=(${USE_HIP} ${USE_HIM} ${USE_HIH})
declare -a HIFRID=(hip him hih)
declare -a HIFRPFX=(HIP HIM HIH)
declare -a HIFRQ=(${Q_HIP} ${Q_HIM} ${Q_HIH})
declare -a HIFRMEM=(${MEM_HIP} ${MEM_HIM} ${MEM_HIH})
declare -a HIFRCPU=(1 1 1)
declare -a IDATEAH=(21 01 01)
declare -a IDATEZH=(28 08 08)
declare -a INITMDH=(6  5  5)
declare -a RUNTYPH=(INITIAL INITIAL INITIAL)
#------------------------------------------------------------------------------------------#


#----- Gridded runs. ----------------------------------------------------------------------#
declare -a USE_GRID=(${USE_RJG})
declare -a GRIDID=(rjg)
declare -a GRIDPFX=(RJG)
declare -a GRIDQ=(${Q_RJG})
declare -a GRIDMEM=(${MEM_RJG})
declare -a GRIDCPU=(${CPU_RJG})
declare -a D_IYEARAG=(2008)
declare -a D_IYEARZG=(2008)
declare -a D_INITMDG=(5)
declare -a D_RUNTYPG=(INITIAL)
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#      Additional settings, which depend on whether this is rapid or long tests.           #
#------------------------------------------------------------------------------------------#
case ${TESTTYPE} in
rapid)

   echo " - Performing rapid tests (2 years for POI, 1 year for grid)"

   #----- POI tests will run for two years. -----------------------------------------------#
   declare -a IYEARAS=(1500 1500 2007 1500 2000 2000 2002 1500 2005 2007 2000 2011 2008)
   declare -a IYEARZS=(1502 1502 2009 1502 2002 2002 2004 1502 2007 2009 2002 2013 2010)
   declare -a INITMDS=(5    0    6    0    5    0    5    0    6    0    5    6    0   )
   declare -a RUNTYPS=(INITIAL INITIAL INITIAL INITIAL INITIAL INITIAL INITIAL \
                       INITIAL INITIAL INITIAL INITIAL INITIAL INITIAL )
   #---------------------------------------------------------------------------------------#



   #----- Gridded tests will run for only 1 year. -----------------------------------------#
   declare -a IYEARAG=(2008)
   declare -a IYEARZG=(2008)
   declare -a INITMDG=(5)
   declare -a RUNTYPG=(INITIAL)
   #---------------------------------------------------------------------------------------#


   #----- Monthly output for state files. -------------------------------------------------#
   UNITSTATE=2
   #---------------------------------------------------------------------------------------#
   ;;

medium)

   echo " - Performing intermediate tests (125 years for POI, 15 years for grid)"

   #----- POI tests will run for two years. -----------------------------------------------#
   declare -a IYEARAS=(1975 1975 1975 1975 1975 1975 1975 1975 1975 1975 1975 1975 1975)
   declare -a IYEARZS=(2100 2100 2100 2100 2100 2100 2100 2100 2100 2100 2100 2100 2100)
   declare -a INITMDS=(5    0    6    0    5    0    5    0    6    0    5    6    0   )
   declare -a RUNTYPS=(INITIAL INITIAL INITIAL INITIAL INITIAL INITIAL INITIAL \
                       INITIAL INITIAL INITIAL INITIAL INITIAL INITIAL )
   #---------------------------------------------------------------------------------------#



   #----- Gridded tests will run for only 1 year. -----------------------------------------#
   declare -a IYEARAG=(1995)
   declare -a IYEARZG=(2010)
   declare -a INITMDG=(5)
   declare -a RUNTYPG=(INITIAL)
   #---------------------------------------------------------------------------------------#


   #----- Monthly output for state files. -------------------------------------------------#
   UNITSTATE=3
   #---------------------------------------------------------------------------------------#
   ;;

long)

   echo " - Performing long tests (300 years for POI, 35 year for grid)"

   #----- POI tests will run for 300 years. -----------------------------------------------#
   declare -a IYEARAS=(1800 1800 1800 1800 1800 1800 1800 1800 1800 1800 1800 1800 1800)
   declare -a IYEARZS=(2100 2100 2100 2100 2100 2100 2100 2100 2100 2100 2100 2100 2100)
   declare -a INITMDS=(0    0    0    0    0    0    0    0    0    0    0    0    0   )
   declare -a RUNTYPS=(INITIAL INITIAL INITIAL INITIAL INITIAL INITIAL INITIAL \
                       INITIAL INITIAL INITIAL INITIAL INITIAL INITIAL)
   #---------------------------------------------------------------------------------------#



   #----- Gridded tests will run for 35 years. --------------------------------------------#
   declare -a IYEARAG=(1500)
   declare -a IYEARZG=(1535)
   declare -a INITMDG=(5)
   declare -a RUNTYPG=(INITIAL)
   #---------------------------------------------------------------------------------------#


   #----- Yearly output for state files. --------------------------------------------------#
   UNITSTATE=3
   #---------------------------------------------------------------------------------------#
   ;;
*)
   #----- I don't know what to do... ------------------------------------------------------#
   echo " !Invalid TESTTYPE: ${TESTTYPE}.  Please specify either 'rapid', 'medium' or 'long'."
   exit 99
   #---------------------------------------------------------------------------------------#
   ;;
esac
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#     Loop over all sites, make sure they are properly set.                                #
#------------------------------------------------------------------------------------------#
for i in ${!SITEID[@]}
do
   #---- Let the user know what is going to happen. ---------------------------------------#
   case ${USE_SITE[i]} in
   y|Y)
      echo " - Processing site: ${SITEID[i]}"
      ;;
   n|N)
      echo " - Skipping site: ${SITEID[i]}"
      ;;
   *)
      echo " Improper use specifier: ${SITEID[i]}! Stopping..."
      exit 99
      ;;
   esac
   #---------------------------------------------------------------------------------------#
done
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#     Loop over all grids, make sure they are properly set.                                #
#------------------------------------------------------------------------------------------#
for i in ${!GRIDID[@]}
do
   #---- Let the user know what is going to happen. ---------------------------------------#
   case ${USE_GRID[i]} in
   y|Y)
      echo " - Processing grid: ${GRIDID[i]}"
      ;;
   n|N)
      echo " - Skipping grid: ${GRIDID[i]}"
      ;;
   *)
      echo " Improper use specifier: ${GRIDID[i]}! Stopping..."
      exit 99
      ;;
   esac
   #---------------------------------------------------------------------------------------#
done
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#     Loop over all high frequency sites, make sure they are properly set.                 #
#------------------------------------------------------------------------------------------#
for i in ${!HIFRID[@]}
do
   #---- Let the user know what is going to happen. ---------------------------------------#
   case ${USE_HIFR[i]} in
   y|Y)
      echo " - Processing high-frequency site: ${HIFRID[i]}"
      ;;
   n|N)
      echo " - Skipping high-frequency site: ${HIFRID[i]}"
      ;;
   *)
      echo " Improper use specifier: ${HIFRID[i]}! Stopping..."
      exit 99
      ;;
   esac
   #---------------------------------------------------------------------------------------#
done
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#     Update version string to contain the test type, and re-name the executables so they  #
# have the same names but without the path.                                                #
#------------------------------------------------------------------------------------------#
HERE=$(pwd)
VERSION="${VERSION}_${TESTTYPE}"
MAIN_EXE=$(basename ${MAIN_EXE_PATH})
TEST_EXE=$(basename ${TEST_EXE_PATH})
DBUG_EXE=$(basename ${DBUG_EXE_PATH})
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#     Create the working folder.                                                           #
#------------------------------------------------------------------------------------------#
RUNTYPE=$(echo ${RUNTYPE} | tr '[:lower:]' '[:upper:]')
case "${RUNTYPE}" in
INITIAL)
   if [ -s ${HERE}/${VERSION} ]
   then
      echo ""
      echo "     Delete previous data and generate pristine directory: ${VERSION}? [y|n]"
      echo " In case you answer y ALL previous data will be lost!"
      read flusher
      flusher=$(echo ${flusher} | tr '[:upper:]' '[:lower:]')
      case ${flusher} in
      y|Y)
         #----- Be nice and give a chance to the user to change their minds. --------------#
         echo "Fine, but if you regret later don't say I didn't warn you..."
         echo "I'm giving you a few seconds to kill this script in case you aren't sure."
         delfun=16
         while [ ${delfun} -gt 1 ]
         do
            let delfun=${delfun}-1
            echo "  - Deletion will begin in '${delfun}' seconds..."
            sleep 1
         done
         #---------------------------------------------------------------------------------#


         #----- Purge directory then recreate it. -----------------------------------------#
         find ${HERE}/${VERSION} -print -delete
         #---------------------------------------------------------------------------------#

         ;;
      n|N)
         #---------------------------------------------------------------------------------#
         #     Overwrite stuff.                                                            #
         #---------------------------------------------------------------------------------#
         echo " Keeping the current directory.  Files may be overwritten, beware..."
         ;;
         #---------------------------------------------------------------------------------#
      *)
         #----- User doesn't want to continue, stop script. -------------------------------#
         echo " Invalid option. Hint: 'y' is short of yes, and 'n' is short of no."
         echo " 'Maybe' and 'whatever' have not been implemented yet."
         exit 0
         #---------------------------------------------------------------------------------#
         ;;
      esac
      #------------------------------------------------------------------------------------#
   fi
   #---------------------------------------------------------------------------------------#
   ;; 
esac


#----- Start a new directory, create links and copy xml files. ----------------------------#
mkdir -p ${HERE}/${VERSION}
LNK_MAIN_EXE="${HERE}/${VERSION}/${MAIN_EXE}"
LNK_TEST_EXE="${HERE}/${VERSION}/${TEST_EXE}"
LNK_DBUG_EXE="${HERE}/${VERSION}/${DBUG_EXE}"
ln -sfv ${DATAPATH}            ${HERE}/${VERSION}/edts_datasets
ln -sfv ${MAIN_EXE_PATH}       ${LNK_MAIN_EXE}
ln -sfv ${TEST_EXE_PATH}       ${LNK_TEST_EXE}
ln -sfv ${DBUG_EXE_PATH}       ${LNK_DBUG_EXE}
cp -v  ${HERE}/Templates/*xml ${HERE}/${VERSION}/
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#     Write out the text info to the report directory.                                     #
#------------------------------------------------------------------------------------------#
TEST_TEXT="${HERE}/${VERSION}/test_text.xml"
/bin/rm -fv ${TEST_TEXT}
touch ${TEST_TEXT}
echo "<?xml version=\"1.0\" ?>"                                    >> ${TEST_TEXT}
echo "<description> "                                              >> ${TEST_TEXT}
echo "<branch_version> ${VERSION_BRANCHED_FROM} </branch_version>" >> ${TEST_TEXT}
echo "<tester_name> ${TESTER_NAME} </tester_name>"                 >> ${TEST_TEXT}
echo "<committer_name> ${COMMITTER_NAME} </committer_name>"        >> ${TEST_TEXT}
echo "<test_description> ${TEST_DESCRIPTION} </test_description>"  >> ${TEST_TEXT}
echo "</description>"                                              >> ${TEST_TEXT}
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#     Loop over the different POI cases.                                                   #
#------------------------------------------------------------------------------------------#
ntasks=0
nnodes=1
for i in ${!SITEID[@]}
do
   case ${USE_SITE[i]} in
   y|Y)
      let ntasks=${ntasks}+3
      ;;
   esac
done
for i in ${!HIFRID[@]}
do
   case ${USE_HIFR[i]} in
   y|Y)
      let ntasks=${ntasks}+3
      ;;
   esac
done

if [ ${ntasks} -gt 0 ]
then
   echo " - Total number of POI tasks = ${ntasks}"
   tasksnow=$(awk  'BEGIN { rounded = sprintf("%.0f", '${ntasks}/24' ); print rounded }')
   let nnodes=${nnodes}+${tasksnow}
   echo " - Total number of nodes = "${nnodes}
   mppwidth=$(awk  'BEGIN { rounded = sprintf("%.0f", '${nnodes}*24' ); print rounded }')
   echo " - mppwidth = ${mppwidth}"
else
   echo "- No poi tasks..."
fi
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#     Reset the script with job submission instructions.                                   #
#------------------------------------------------------------------------------------------#
subbatch="${HERE}/${VERSION}/submit_batch.sh"
/bin/rm -f ${subbatch}
touch ${subbatch}
echo "#!/bin/bash"       >> ${subbatch}
echo ". ${HOME}/.bashrc" >> ${subbatch}
#------------------------------------------------------------------------------------------#


#----- Go through all sites, generate command for job queue submission. -------------------#
tcount=0
for i in ${!SITEID[@]}
do
   case ${USE_SITE[i]} in
   y|Y)
      #----- Template and working namelists. ----------------------------------------------#
      TEMPLATEMAIN="${HERE}/Templates/ED2IN-${SITEPFX[i]}-MAIN"
      TEMPLATETEST="${HERE}/Templates/ED2IN-${SITEPFX[i]}-TEST"
      TEMPLATEDBUG="${HERE}/Templates/ED2IN-${SITEPFX[i]}-DBUG"
      FILEMAIN="${HERE}/${VERSION}/ED2IN-${SITEPFX[i]}-MAIN"
      FILETEST="${HERE}/${VERSION}/ED2IN-${SITEPFX[i]}-TEST"
      FILEDBUG="${HERE}/${VERSION}/ED2IN-${SITEPFX[i]}-DBUG"
      #------------------------------------------------------------------------------------#

      #----- Entertain user: --------------------------------------------------------------#
      echo ""
      echo ""
      echo " + Processing ${SITEID[i]}"
      echo "   - Path:${HERE}/${VERSION}"
      echo "   - Control Files:"
      echo "     * Main:  $(basename ${FILEMAIN})"
      echo "     * Test:  $(basename ${FILETEST})"
      echo "     * Debug: $(basename ${FILEDBUG})"
      echo ""
      #------------------------------------------------------------------------------------#


      #----- Copy the Template to the working directory. ----------------------------------#
      /bin/cp -f ${TEMPLATEMAIN} ${FILEMAIN}
      /bin/cp -f ${TEMPLATETEST} ${FILETEST}
      /bin/cp -f ${TEMPLATEDBUG} ${FILEDBUG}
      #------------------------------------------------------------------------------------#


      #----- Output folders. --------------------------------------------------------------#
      F_MAIN_PATH="${HERE}/${VERSION}/F_main_${SITEID[i]}"
      S_MAIN_PATH="${HERE}/${VERSION}/S_main_${SITEID[i]}"
      F_TEST_PATH="${HERE}/${VERSION}/F_test_${SITEID[i]}"
      S_TEST_PATH="${HERE}/${VERSION}/S_test_${SITEID[i]}"
      F_DBUG_PATH="${HERE}/${VERSION}/F_dbug_${SITEID[i]}"
      S_DBUG_PATH="${HERE}/${VERSION}/S_dbug_${SITEID[i]}"
      #------------------------------------------------------------------------------------#


      #----- Output prefixes. -------------------------------------------------------------#
      F_MAIN_PREF="${F_MAIN_PATH}/main_${SITEID[i]}"
      S_MAIN_PREF="${S_MAIN_PATH}/main_${SITEID[i]}"
      F_TEST_PREF="${F_TEST_PATH}/test_${SITEID[i]}"
      S_TEST_PREF="${S_TEST_PATH}/test_${SITEID[i]}"
      F_DBUG_PREF="${F_DBUG_PATH}/dbug_${SITEID[i]}"
      S_DBUG_PREF="${S_DBUG_PATH}/dbug_${SITEID[i]}"
      #------------------------------------------------------------------------------------#


      #----- Reset and flush the output folders. ------------------------------------------#
      mkdir -p ${F_MAIN_PATH}
      mkdir -p ${F_TEST_PATH}
      mkdir -p ${F_DBUG_PATH}
      mkdir -p ${S_MAIN_PATH}
      mkdir -p ${S_TEST_PATH}
      mkdir -p ${S_DBUG_PATH}
      #------------------------------------------------------------------------------------#



      #----- Decide whether to continue the run or start from beginning. ------------------#
      case ${RUNTYPE} in
      INITIAL)
         #----- Use default initialisation. -----------------------------------------------#
         sed -i '/NL%RUNTYPE/c\   NL%RUNTYPE = '\'${RUNTYPS[i]}\'   ${FILEMAIN}
         sed -i '/NL%RUNTYPE/c\   NL%RUNTYPE = '\'${RUNTYPS[i]}\'   ${FILETEST}
         sed -i '/NL%RUNTYPE/c\   NL%RUNTYPE = '\'${D_RUNTYPS[i]}\' ${FILEDBUG}
         #---------------------------------------------------------------------------------#

         #----- Initial.  All simulations must be submitted. ------------------------------#
         SUBMIT_MAIN="y"
         SUBMIT_TEST="y"
         SUBMIT_DBUG="y"
         #---------------------------------------------------------------------------------#

         ;;
      HISTORY)
         #----- Look for history files (MAIN) then decide how to set up the run. ----------#
         LHIST_MAIN=$(/bin/ls -1 ${S_MAIN_PREF}-S-*h5 2> /dev/null | xargs -n 1 basename   \
                      2> /dev/null | tail -1)
         if [ "${LHIST_MAIN}" != "" ]
         then
            IYEARH_MAIN=$(echo  ${LHIST_MAIN} | cut -c 12-15)
            IMONTHH_MAIN=$(echo ${LHIST_MAIN} | cut -c 17-18)
            IDATEH_MAIN=$(echo  ${LHIST_MAIN} | cut -c 20-21)
            ITIMEH_MAIN=$(echo  ${LHIST_MAIN} | cut -c 23-26)
            sed -i '/NL%RUNTYPE/c\   NL%RUNTYPE = '\''HISTORY'\'      ${FILEMAIN}
            sed -i '/NL%IYEARH/c\   NL%IYEARH = '${IYEARH_MAIN}   ${FILEMAIN}
            sed -i '/NL%IMONTHH/c\   NL%IMONTHH = '${IMONTHH_MAIN} ${FILEMAIN}
            sed -i '/NL%IDATEH/c\   NL%IDATEH = '${IDATEH_MAIN}   ${FILEMAIN}
            sed -i '/NL%ITIMEH/c\   NL%ITIMEH = '${ITIMEH_MAIN}   ${FILEMAIN}
            sed -i '/NL%SFILIN/c\   NL%SFILIN = '\'${S_MAIN_PREF}\'   ${FILEMAIN}
         else
            sed -i '/NL%RUNTYPE/c\   NL%RUNTYPE = '\'${RUNTYPS[i]}\'   ${FILEMAIN}
         fi
         #---------------------------------------------------------------------------------#



         #----- Look for history files (TEST) then decide how to set up the run. ----------#
         LHIST_TEST=$(/bin/ls -1 ${S_TEST_PREF}-S-*h5 2> /dev/null | xargs -n 1 basename   \
                      2> /dev/null | tail -1)
         
         if [ "${LHIST_TEST}" != "" ]
         then
            IYEARH_TEST=$(echo  ${LHIST_TEST} | cut -c 12-15)
            IMONTHH_TEST=$(echo ${LHIST_TEST} | cut -c 17-18)
            IDATEH_TEST=$(echo  ${LHIST_TEST} | cut -c 20-21)
            ITIMEH_TEST=$(echo  ${LHIST_TEST} | cut -c 23-26)
            sed -i '/NL%RUNTYPE/c\   NL%RUNTYPE = '\''HISTORY'\'      ${FILETEST}
            sed -i '/NL%IYEARH/c\   NL%IYEARH = '${IYEARH_TEST}   ${FILETEST}
            sed -i '/NL%IMONTHH/c\   NL%IMONTHH = '${IMONTHH_TEST} ${FILETEST}
            sed -i '/NL%IDATEH/c\   NL%IDATEH = '${IDATEH_TEST}   ${FILETEST}
            sed -i '/NL%ITIMEH/c\   NL%ITIMEH = '${ITIMEH_TEST}   ${FILETEST}
            sed -i '/NL%SFILIN/c\   NL%SFILIN = '\'${S_TEST_PREF}\'   ${FILETEST}
         else
            sed -i '/NL%RUNTYPE/c\   NL%RUNTYPE = '\'${RUNTYPS[i]}\'   ${FILETEST}
         fi
         #---------------------------------------------------------------------------------#



         #----- Look for history files (DBUG) then decide how to set up the run. ----------#
         LHIST_DBUG=$(/bin/ls -1 ${S_DBUG_PREF}-S-*h5 2> /dev/null | xargs -n 1 basename   \
                      2> /dev/null | tail -1)
         if [ "${LHIST_DBUG}" != "" ]
         then
            IYEARH_DBUG=$(echo  ${LHIST_DBUG} | cut -c 12-15)
            IMONTHH_DBUG=$(echo ${LHIST_DBUG} | cut -c 17-18)
            IDATEH_DBUG=$(echo  ${LHIST_DBUG} | cut -c 20-21)
            ITIMEH_DBUG=$(echo  ${LHIST_DBUG} | cut -c 23-26)
            sed -i '/NL%RUNTYPE/c\   NL%RUNTYPE = '\''HISTORY'\'         ${FILEDBUG}
            sed -i '/NL%IYEARH/c\   NL%IYEARH = '${IYEARH_DBUG}      ${FILEDBUG}
            sed -i '/NL%IMONTHH/c\   NL%IMONTHH = '${IMONTHH_DBUG}    ${FILEDBUG}
            sed -i '/NL%IDATEH/c\   NL%IDATEH = '${IDATEH_DBUG}      ${FILEDBUG}
            sed -i '/NL%ITIMEH/c\   NL%ITIMEH = '${ITIMEH_DBUG}      ${FILEDBUG}
            sed -i '/NL%SFILIN/c\   NL%SFILIN = '\'${S_DBUG_PREF}\'      ${FILEDBUG}
         else
            sed -i '/NL%RUNTYPE/c\   NL%RUNTYPE = '\'${D_RUNTYPS[i]}\'   ${FILEDBUG}
         fi
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #     Find out whether to submit the run or not.  Default is yes, unless it has   #
         # reached the end.                                                                #
         #---------------------------------------------------------------------------------#
         #----- Main. ---------------------------------------------------------------------#
         STDOUT_MAIN="${HERE}/${VERSION}/main_${SITEID[i]}.out"
         THE_END_MAIN=$(grep "ED-2.2 execution ends" ${STDOUT_MAIN} 2> /dev/null | wc -l)
         if [ ${THE_END_MAIN} -gt 0 ]
         then 
            SUBMIT_MAIN="n"
         else
            SUBMIT_MAIN="y"
         fi
         #----- Test. ---------------------------------------------------------------------#
         STDOUT_TEST="${HERE}/${VERSION}/test_${SITEID[i]}.out"
         THE_END_TEST=$(grep "ED-2.2 execution ends" ${STDOUT_TEST} 2> /dev/null | wc -l)
         if [ ${THE_END_TEST} -gt 0 ]
         then 
            SUBMIT_TEST="n"
         else
            SUBMIT_TEST="y"
         fi
         #----- DBUG. ---------------------------------------------------------------------#
         STDOUT_DBUG="${HERE}/${VERSION}/dbug_${SITEID[i]}.out"
         THE_END_DBUG=$(grep "ED-2.2 execution ends" ${STDOUT_DBUG} 2> /dev/null | wc -l)
         if [ ${THE_END_DBUG} -gt 0 ]
         then 
            SUBMIT_DBUG="n"
         else
            SUBMIT_DBUG="y"
         fi
         #---------------------------------------------------------------------------------#
         ;;
      esac
      #------------------------------------------------------------------------------------#


      #----- Modify the ED2IN IED_INIT_MODE. ----------------------------------------------#
      sed -i '/NL%IED_INIT_MODE/c\   NL%IED_INIT_MODE = '${INITMDS[i]}   ${FILEMAIN}
      sed -i '/NL%IED_INIT_MODE/c\   NL%IED_INIT_MODE = '${INITMDS[i]}   ${FILETEST}
      sed -i '/NL%IED_INIT_MODE/c\   NL%IED_INIT_MODE = '${D_INITMDS[i]} ${FILEDBUG}
      #------------------------------------------------------------------------------------#


      #----- Modify the ED2IN start years. ------------------------------------------------#
      sed -i '/NL%IYEARA/c\   NL%IYEARA   = '${IYEARAS[i]}   ${FILEMAIN}
      sed -i '/NL%IYEARA/c\   NL%IYEARA   = '${IYEARAS[i]}   ${FILETEST}
      sed -i '/NL%IYEARA/c\   NL%IYEARA   = '${D_IYEARAS[i]} ${FILEDBUG}
      #------------------------------------------------------------------------------------#


      #----- Modify the ED2IN end years. --------------------------------------------------#
      sed -i '/NL%IYEARZ/c\   NL%IYEARZ   = '${IYEARZS[i]}   ${FILEMAIN}
      sed -i '/NL%IYEARZ/c\   NL%IYEARZ   = '${IYEARZS[i]}   ${FILETEST}
      sed -i '/NL%IYEARZ/c\   NL%IYEARZ   = '${D_IYEARZS[i]} ${FILEDBUG}
      #------------------------------------------------------------------------------------#


      #----- Modify the ED2IN files to point to the desired output directories. -----------#
      sed -i '/NL%FFILOUT/c\   NL%FFILOUT = '\'${F_MAIN_PREF}\' ${FILEMAIN}
      sed -i '/NL%SFILOUT/c\   NL%SFILOUT = '\'${S_MAIN_PREF}\' ${FILEMAIN}
      sed -i '/NL%FFILOUT/c\   NL%FFILOUT = '\'${F_TEST_PREF}\' ${FILETEST}
      sed -i '/NL%SFILOUT/c\   NL%SFILOUT = '\'${S_TEST_PREF}\' ${FILETEST}
      sed -i '/NL%FFILOUT/c\   NL%FFILOUT = '\'${F_DBUG_PREF}\' ${FILEDBUG}
      sed -i '/NL%SFILOUT/c\   NL%SFILOUT = '\'${S_DBUG_PREF}\' ${FILEDBUG}
      #------------------------------------------------------------------------------------#


      #----- Modify the state file output frequency. --------------------------------------#
      sed -i '/NL%UNITSTATE/c\   NL%UNITSTATE  = '${UNITSTATE} ${FILEMAIN}
      sed -i '/NL%UNITSTATE/c\   NL%UNITSTATE  = '${UNITSTATE} ${FILETEST}
      sed -i '/NL%UNITSTATE/c\   NL%UNITSTATE  = '${UNITSTATE} ${FILEDBUG}
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #   Make sure medium and long tests won't write fast or daily files.                 #
      #------------------------------------------------------------------------------------#
      case ${TESTTYPE} in
      medium|long)
         sed -i '/NL%IFOUTPUT/c\   NL%IFOUTPUT  = 0' ${FILETEST}
         sed -i '/NL%IDOUTPUT/c\   NL%IDOUTPUT  = 0' ${FILETEST}
         sed -i '/NL%IFOUTPUT/c\   NL%IFOUTPUT  = 0' ${FILEMAIN}
         sed -i '/NL%IDOUTPUT/c\   NL%IDOUTPUT  = 0' ${FILEMAIN}
         ;;
      esac
      #------------------------------------------------------------------------------------#

      #------------------------------------------------------------------------------------#
      #    Prepare job (MAIN)
      #------------------------------------------------------------------------------------#
      case ${SUBMIT_MAIN} in
      y|Y)
         jobout=${HERE}/${VERSION}/main_${SITEID[i]}.out
         joberr=${HERE}/${VERSION}/main_${SITEID[i]}.err
         jobname=${VERSION}_${SITEID[i]}_main
         jobopts="-t ${POI_TIME} --mem-per-cpu=${SITEMEM[i]} --cpus-per-task=${SITECPU[i]}"
         jobopts="${jobopts} -p ${SITEQ[i]} -n 1"
         jobwrap=". ${HOME}/.bashrc; cd ${HERE}/${VERSION}"
         jobwrap="${jobwrap}; export OMP_NUM_THREADS=${SITECPU[i]}"
         jobwrap="${jobwrap}; mpirun -np 1 ${LNK_MAIN_EXE} -f ${FILEMAIN}"
         jobwrap="\"(${jobwrap})\""
         jobcomm="sbatch -o ${jobout} -e ${joberr} -J ${jobname}"
         jobcomm="${jobcomm} ${jobopts} --wrap=${jobwrap}"
         echo ${jobcomm} >> ${subbatch}
         ;;
      *)
         echo "     * Main has already reached the end and won't be submitted."
         ;;
      esac
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #    Prepare job (TEST)
      #------------------------------------------------------------------------------------#
      case ${SUBMIT_TEST} in
      y|Y)
         jobout=${HERE}/${VERSION}/test_${SITEID[i]}.out
         joberr=${HERE}/${VERSION}/test_${SITEID[i]}.err
         jobname=${VERSION}_${SITEID[i]}_test
         jobopts="-t ${POI_TIME} --mem-per-cpu=${SITEMEM[i]} --cpus-per-task=${SITECPU[i]}"
         jobopts="${jobopts} -p ${SITEQ[i]} -n 1"
         jobwrap=". ${HOME}/.bashrc; cd ${HERE}/${VERSION}"
         jobwrap="${jobwrap}; export OMP_NUM_THREADS=${SITECPU[i]}"
         jobwrap="${jobwrap}; mpirun -np 1 ${LNK_TEST_EXE} -f ${FILETEST}"
         jobwrap="\"(${jobwrap})\""
         jobcomm="sbatch -o ${jobout} -e ${joberr} -J ${jobname}"
         jobcomm="${jobcomm} ${jobopts} --wrap=${jobwrap}"
         echo ${jobcomm} >> ${subbatch}
         ;;
      *)
         echo "     * Test has already reached the end and won't be submitted."
         ;;
      esac
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #    Prepare job (DBUG)
      #------------------------------------------------------------------------------------#
      case ${SUBMIT_DBUG} in
      y|Y)
         jobout=${HERE}/${VERSION}/dbug_${SITEID[i]}.out
         joberr=${HERE}/${VERSION}/dbug_${SITEID[i]}.err
         jobname=${VERSION}_${SITEID[i]}_dbug
         jobopts="-t ${POI_TIME} --mem-per-cpu=${SITEMEM[i]} --cpus-per-task=${SITECPU[i]}"
         jobopts="${jobopts} -p ${SITEQ[i]} -n 1"
         jobwrap=". ${HOME}/.bashrc; cd ${HERE}/${VERSION}"
         jobwrap="${jobwrap}; export OMP_NUM_THREADS=${SITECPU[i]}"
         jobwrap="${jobwrap}; mpirun -np 1 ${LNK_DBUG_EXE} -f ${FILEDBUG}"
         jobwrap="\"(${jobwrap})\""
         jobcomm="sbatch -o ${jobout} -e ${joberr} -J ${jobname}"
         jobcomm="${jobcomm} ${jobopts} --wrap=${jobwrap}"
         echo ${jobcomm} >> ${subbatch}
         ;;
      *)
         echo "     * Dbug has already reached the end and won't be submitted."
         ;;
      esac
      #------------------------------------------------------------------------------------#


      #----- Update counter. --------------------------------------------------------------
      let tcount=${tcount}+3
      #------------------------------------------------------------------------------------#

      ;;
   esac
   #---------------------------------------------------------------------------------------#
done
#------------------------------------------------------------------------------------------#



#----- Go through all high-frequency sites, generate command for job queue submission. ----#
for i in ${!HIFRID[@]}
do
   case ${USE_HIFR[i]} in
   y|Y)
      #----- Template and working namelists. ----------------------------------------------#
      TEMPLATEMAIN="${HERE}/Templates/ED2IN-${HIFRPFX[i]}-MAIN"
      TEMPLATETEST="${HERE}/Templates/ED2IN-${HIFRPFX[i]}-TEST"
      TEMPLATEDBUG="${HERE}/Templates/ED2IN-${HIFRPFX[i]}-DBUG"
      FILEMAIN="${HERE}/${VERSION}/ED2IN-${HIFRPFX[i]}-MAIN"
      FILETEST="${HERE}/${VERSION}/ED2IN-${HIFRPFX[i]}-TEST"
      FILEDBUG="${HERE}/${VERSION}/ED2IN-${HIFRPFX[i]}-DBUG"
      #------------------------------------------------------------------------------------#

      #----- Entertain user: --------------------------------------------------------------#
      echo ""
      echo ""
      echo " + Processing ${HIFRID[i]}"
      echo "   - Path:${HERE}/${VERSION}"
      echo "   - Control Files:"
      echo "     * Main:  $(basename ${FILEMAIN})"
      echo "     * Test:  $(basename ${FILETEST})"
      echo "     * Debug: $(basename ${FILEDBUG})"
      echo ""
      #------------------------------------------------------------------------------------#


      #----- Copy the Template to the working directory. ----------------------------------#
      /bin/cp -f ${TEMPLATEMAIN} ${FILEMAIN}
      /bin/cp -f ${TEMPLATETEST} ${FILETEST}
      /bin/cp -f ${TEMPLATEDBUG} ${FILEDBUG}
      #------------------------------------------------------------------------------------#


      #----- Output folders. --------------------------------------------------------------#
      F_MAIN_PATH="${HERE}/${VERSION}/F_main_${HIFRID[i]}"
      S_MAIN_PATH="${HERE}/${VERSION}/S_main_${HIFRID[i]}"
      F_TEST_PATH="${HERE}/${VERSION}/F_test_${HIFRID[i]}"
      S_TEST_PATH="${HERE}/${VERSION}/S_test_${HIFRID[i]}"
      F_DBUG_PATH="${HERE}/${VERSION}/F_dbug_${HIFRID[i]}"
      S_DBUG_PATH="${HERE}/${VERSION}/S_dbug_${HIFRID[i]}"
      #------------------------------------------------------------------------------------#


      #----- Output prefixes. -------------------------------------------------------------#
      F_MAIN_PREF="${F_MAIN_PATH}/main_${HIFRID[i]}"
      S_MAIN_PREF="${S_MAIN_PATH}/main_${HIFRID[i]}"
      F_TEST_PREF="${F_TEST_PATH}/test_${HIFRID[i]}"
      S_TEST_PREF="${S_TEST_PATH}/test_${HIFRID[i]}"
      F_DBUG_PREF="${F_DBUG_PATH}/dbug_${HIFRID[i]}"
      S_DBUG_PREF="${S_DBUG_PATH}/dbug_${HIFRID[i]}"
      #------------------------------------------------------------------------------------#


      #----- Reset and flush the output folders. ------------------------------------------#
      mkdir -p ${F_MAIN_PATH}
      mkdir -p ${F_TEST_PATH}
      mkdir -p ${F_DBUG_PATH}
      mkdir -p ${S_MAIN_PATH}
      mkdir -p ${S_TEST_PATH}
      mkdir -p ${S_DBUG_PATH}
      #------------------------------------------------------------------------------------#



      #----- Decide whether to continue the run or start from beginning. ------------------#
      case ${RUNTYPE} in
      INITIAL)
         #----- Use default initialisation. -----------------------------------------------#
         sed -i '/NL%RUNTYPE/c\   NL%RUNTYPE = '\'${RUNTYPH[i]}\'   ${FILEMAIN}
         sed -i '/NL%RUNTYPE/c\   NL%RUNTYPE = '\'${RUNTYPH[i]}\'   ${FILETEST}
         sed -i '/NL%RUNTYPE/c\   NL%RUNTYPE = '\'${D_RUNTYPH[i]}\' ${FILEDBUG}
         #---------------------------------------------------------------------------------#

         #----- Initial.  All simulations must be submitted. ------------------------------#
         SUBMIT_MAIN="y"
         SUBMIT_TEST="y"
         SUBMIT_DBUG="y"
         #---------------------------------------------------------------------------------#

         ;;
      HISTORY)
         #----- Look for history files (MAIN) then decide how to set up the run. ----------#
         LHIST_MAIN=$(/bin/ls -1 ${S_MAIN_PREF}-S-*h5 2> /dev/null | xargs -n 1 basename   \
                      2> /dev/null | tail -1)
         if [ "${LHIST_MAIN}" != "" ]
         then
            IYEARH_MAIN=$(echo  ${LHIST_MAIN} | cut -c 12-15)
            IMONTHH_MAIN=$(echo ${LHIST_MAIN} | cut -c 17-18)
            IDATEH_MAIN=$(echo  ${LHIST_MAIN} | cut -c 20-21)
            ITIMEH_MAIN=$(echo  ${LHIST_MAIN} | cut -c 23-26)
            sed -i '/NL%RUNTYPE/c\   NL%RUNTYPE = '\''HISTORY'\'      ${FILEMAIN}
            sed -i '/NL%IYEARH/c\   NL%IYEARH = '${IYEARH_MAIN}   ${FILEMAIN}
            sed -i '/NL%IMONTHH/c\   NL%IMONTHH = '${IMONTHH_MAIN} ${FILEMAIN}
            sed -i '/NL%IDATEH/c\   NL%IDATEH = '${IDATEH_MAIN}   ${FILEMAIN}
            sed -i '/NL%ITIMEH/c\   NL%ITIMEH = '${ITIMEH_MAIN}   ${FILEMAIN}
            sed -i '/NL%SFILIN/c\   NL%SFILIN = '\'${S_MAIN_PREF}\'   ${FILEMAIN}
         else
            sed -i '/NL%RUNTYPE/c\   NL%RUNTYPE = '\'${RUNTYPH[i]}\'   ${FILEMAIN}
         fi
         #---------------------------------------------------------------------------------#



         #----- Look for history files (TEST) then decide how to set up the run. ----------#
         LHIST_TEST=$(/bin/ls -1 ${S_TEST_PREF}-S-*h5 2> /dev/null | xargs -n 1 basename   \
                      2> /dev/null | tail -1)
         if [ "${LHIST_TEST}" != "" ]
         then
            IYEARH_TEST=$(echo  ${LHIST_TEST} | cut -c 12-15)
            IMONTHH_TEST=$(echo ${LHIST_TEST} | cut -c 17-18)
            IDATEH_TEST=$(echo  ${LHIST_TEST} | cut -c 20-21)
            ITIMEH_TEST=$(echo  ${LHIST_TEST} | cut -c 23-26)
            sed -i '/NL%RUNTYPE/c\   NL%RUNTYPE = '\''HISTORY'\'      ${FILETEST}
            sed -i '/NL%IYEARH/c\   NL%IYEARH = '${IYEARH_TEST}   ${FILETEST}
            sed -i '/NL%IMONTHH/c\   NL%IMONTHH = '${IMONTHH_TEST} ${FILETEST}
            sed -i '/NL%IDATEH/c\   NL%IDATEH = '${IDATEH_TEST}   ${FILETEST}
            sed -i '/NL%ITIMEH/c\   NL%ITIMEH = '${ITIMEH_TEST}   ${FILETEST}
            sed -i '/NL%SFILIN/c\   NL%SFILIN = '\'${S_TEST_PREF}\'   ${FILETEST}
         else
            sed -i '/NL%RUNTYPE/c\   NL%RUNTYPE = '\'${RUNTYPH[i]}\'   ${FILETEST}
         fi
         #---------------------------------------------------------------------------------#



         #----- Look for history files (DBUG) then decide how to set up the run. ----------#
         LHIST_DBUG=$(/bin/ls -1 ${S_DBUG_PREF}-S-*h5 2> /dev/null | xargs -n 1 basename   \
                      2> /dev/null | tail -1)
         if [ "${LHIST_DBUG}" != "" ]
         then
            IYEARH_DBUG=$(echo  ${LHIST_DBUG} | cut -c 12-15)
            IMONTHH_DBUG=$(echo ${LHIST_DBUG} | cut -c 17-18)
            IDATEH_DBUG=$(echo  ${LHIST_DBUG} | cut -c 20-21)
            ITIMEH_DBUG=$(echo  ${LHIST_DBUG} | cut -c 23-26)
            sed -i '/NL%RUNTYPE/c\   NL%RUNTYPE = '\''HISTORY'\'         ${FILEDBUG}
            sed -i '/NL%IYEARH/c\   NL%IYEARH = '${IYEARH_DBUG}      ${FILEDBUG}
            sed -i '/NL%IMONTHH/c\   NL%IMONTHH = '${IMONTHH_DBUG}    ${FILEDBUG}
            sed -i '/NL%IDATEH/c\   NL%IDATEH = '${IDATEH_DBUG}      ${FILEDBUG}
            sed -i '/NL%ITIMEH/c\   NL%ITIMEH = '${ITIMEH_DBUG}      ${FILEDBUG}
            sed -i '/NL%SFILIN/c\   NL%SFILIN = '\'${S_DBUG_PREF}\'      ${FILEDBUG}
         else
            sed -i '/NL%RUNTYPE/c\   NL%RUNTYPE = '\'${D_RUNTYPH[i]}\'   ${FILEDBUG}
         fi
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Find out whether to submit the run or not.  Default is yes, unless it has   #
         # reached the end.                                                                #
         #---------------------------------------------------------------------------------#
         #----- Main. ---------------------------------------------------------------------#
         STDOUT_MAIN="${HERE}/${VERSION}/main_${HIFRID[i]}.out"
         THE_END_MAIN=$(grep "ED-2.2 execution ends" ${STDOUT_MAIN} 2> /dev/null | wc -l)
         if [ ${THE_END_MAIN} -gt 0 ]
         then 
            SUBMIT_MAIN="n"
         else
            SUBMIT_MAIN="y"
         fi
         #----- Test. ---------------------------------------------------------------------#
         STDOUT_TEST="${HERE}/${VERSION}/test_${HIFRID[i]}.out"
         THE_END_TEST=$(grep "ED-2.2 execution ends" ${STDOUT_TEST} 2> /dev/null | wc -l)
         if [ ${THE_END_TEST} -gt 0 ]
         then 
            SUBMIT_TEST="n"
         else
            SUBMIT_TEST="y"
         fi
         #----- DBUG. ---------------------------------------------------------------------#
         STDOUT_DBUG="${HERE}/${VERSION}/dbug_${HIFRID[i]}.out"
         THE_END_DBUG=$(grep "ED-2.2 execution ends" ${STDOUT_DBUG} 2> /dev/null | wc -l)
         if [ ${THE_END_DBUG} -gt 0 ]
         then 
            SUBMIT_DBUG="n"
         else
            SUBMIT_DBUG="y"
         fi
         #---------------------------------------------------------------------------------#
         ;;
      esac
      #------------------------------------------------------------------------------------#


      #----- Modify the ED2IN IED_INIT_MODE. ----------------------------------------------#
      sed -i '/NL%IED_INIT_MODE/c\   NL%IED_INIT_MODE = '${INITMDH[i]} ${FILEMAIN}
      sed -i '/NL%IED_INIT_MODE/c\   NL%IED_INIT_MODE = '${INITMDH[i]} ${FILETEST}
      sed -i '/NL%IED_INIT_MODE/c\   NL%IED_INIT_MODE = '${INITMDH[i]} ${FILEDBUG}
      #------------------------------------------------------------------------------------#


      #----- Modify the ED2IN start dates. ------------------------------------------------#
      sed -i '/NL%IDATEA/c\   NL%IDATEA   = '${IDATEAH[i]} ${FILEMAIN}
      sed -i '/NL%IDATEA/c\   NL%IDATEA   = '${IDATEAH[i]} ${FILETEST}
      sed -i '/NL%IDATEA/c\   NL%IDATEA   = '${IDATEAH[i]} ${FILEDBUG}
      #------------------------------------------------------------------------------------#


      #----- Modify the ED2IN end DATEs. --------------------------------------------------#
      sed -i '/NL%IDATEZ/c\   NL%IDATEZ   = '${IDATEZH[i]} ${FILEMAIN}
      sed -i '/NL%IDATEZ/c\   NL%IDATEZ   = '${IDATEZH[i]} ${FILETEST}
      sed -i '/NL%IDATEZ/c\   NL%IDATEZ   = '${IDATEZH[i]} ${FILEDBUG}
      #------------------------------------------------------------------------------------#


      #----- Modify the ED2IN files to point to the desired output directories. -----------#
      sed -i '/NL%FFILOUT/c\   NL%FFILOUT = '\'${F_MAIN_PREF}\' ${FILEMAIN}
      sed -i '/NL%SFILOUT/c\   NL%SFILOUT = '\'${S_MAIN_PREF}\' ${FILEMAIN}
      sed -i '/NL%FFILOUT/c\   NL%FFILOUT = '\'${F_TEST_PREF}\' ${FILETEST}
      sed -i '/NL%SFILOUT/c\   NL%SFILOUT = '\'${S_TEST_PREF}\' ${FILETEST}
      sed -i '/NL%FFILOUT/c\   NL%FFILOUT = '\'${F_DBUG_PREF}\' ${FILEDBUG}
      sed -i '/NL%SFILOUT/c\   NL%SFILOUT = '\'${S_DBUG_PREF}\' ${FILEDBUG}
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #    Prepare job (MAIN)
      #------------------------------------------------------------------------------------#
      case ${SUBMIT_MAIN} in
      y|Y)
         jobout=${HERE}/${VERSION}/main_${HIFRID[i]}.out
         joberr=${HERE}/${VERSION}/main_${HIFRID[i]}.err
         jobname=${VERSION}_${HIFRID[i]}_main
         jobopts="-t ${POI_TIME} --mem-per-cpu=${HIFRMEM[i]} --cpus-per-task=${HIFRCPU[i]}"
         jobopts="${jobopts} -p ${HIFRQ[i]} -n 1"
         jobwrap=". ${HOME}/.bashrc; cd ${HERE}/${VERSION}"
         jobwrap="${jobwrap}; export OMP_NUM_THREADS=${HIFRCPU[i]}"
         jobwrap="${jobwrap}; mpirun -np 1 ${LNK_MAIN_EXE} -f ${FILEMAIN}"
         jobwrap="\"(${jobwrap})\""
         jobcomm="sbatch -o ${jobout} -e ${joberr} -J ${jobname}"
         jobcomm="${jobcomm} ${jobopts} --wrap=${jobwrap}"
         echo ${jobcomm} >> ${subbatch}
         ;;
      *)
         echo "     * Main has already reached the end and won't be submitted."
         ;;
      esac
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #    Prepare job (TEST)
      #------------------------------------------------------------------------------------#
      case ${SUBMIT_TEST} in
      y|Y)
         jobout=${HERE}/${VERSION}/test_${HIFRID[i]}.out
         joberr=${HERE}/${VERSION}/test_${HIFRID[i]}.err
         jobname=${VERSION}_${HIFRID[i]}_test
         jobopts="-t ${POI_TIME} --mem-per-cpu=${HIFRMEM[i]} --cpus-per-task=${HIFRCPU[i]}"
         jobopts="${jobopts} -p ${HIFRQ[i]} -n 1"
         jobwrap=". ${HOME}/.bashrc; cd ${HERE}/${VERSION}"
         jobwrap="${jobwrap}; export OMP_NUM_THREADS=${HIFRCPU[i]}"
         jobwrap="${jobwrap}; mpirun -np 1 ${LNK_TEST_EXE} -f ${FILETEST}"
         jobwrap="\"(${jobwrap})\""
         jobcomm="sbatch -o ${jobout} -e ${joberr} -J ${jobname}"
         jobcomm="${jobcomm} ${jobopts} --wrap=${jobwrap}"
         echo ${jobcomm} >> ${subbatch}
         ;;
      *)
         echo "     * Test has already reached the end and won't be submitted."
         ;;
      esac
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #    Prepare job (DBUG)
      #------------------------------------------------------------------------------------#
      case ${SUBMIT_DBUG} in
      y|Y)
         jobout=${HERE}/${VERSION}/dbug_${HIFRID[i]}.out
         joberr=${HERE}/${VERSION}/dbug_${HIFRID[i]}.err
         jobname=${VERSION}_${HIFRID[i]}_dbug
         jobopts="-t ${POI_TIME} --mem-per-cpu=${HIFRMEM[i]} --cpus-per-task=${HIFRCPU[i]}"
         jobopts="${jobopts} -p ${HIFRQ[i]} -n 1"
         jobwrap=". ${HOME}/.bashrc; cd ${HERE}/${VERSION}"
         jobwrap="${jobwrap}; export OMP_NUM_THREADS=${HIFRCPU[i]}"
         jobwrap="${jobwrap}; mpirun -np 1 ${LNK_DBUG_EXE} -f ${FILEDBUG}"
         jobwrap="\"(${jobwrap})\""
         jobcomm="sbatch -o ${jobout} -e ${joberr} -J ${jobname}"
         jobcomm="${jobcomm} ${jobopts} --wrap=${jobwrap}"
         echo ${jobcomm} >> ${subbatch}
         ;;
      *)
         echo "     * Dbug has already reached the end and won't be submitted."
         ;;
      esac
      #------------------------------------------------------------------------------------#

      ;;
   esac
   #---------------------------------------------------------------------------------------#
done
#------------------------------------------------------------------------------------------#



#----- Go through all grids, generate command for job queue submission. -------------------#
for i in ${!GRIDID[@]}
do
   case ${USE_GRID[i]} in
   y|Y)
      #----- Template and working namelists. ----------------------------------------------#
      TEMPLATEMAIN="${HERE}/Templates/ED2IN-${GRIDPFX[i]}-MAIN"
      TEMPLATETEST="${HERE}/Templates/ED2IN-${GRIDPFX[i]}-TEST"
      TEMPLATEDBUG="${HERE}/Templates/ED2IN-${GRIDPFX[i]}-DBUG"
      FILEMAIN="${HERE}/${VERSION}/ED2IN-${GRIDPFX[i]}-MAIN"
      FILETEST="${HERE}/${VERSION}/ED2IN-${GRIDPFX[i]}-TEST"
      FILEDBUG="${HERE}/${VERSION}/ED2IN-${GRIDPFX[i]}-DBUG"
      #------------------------------------------------------------------------------------#

      #----- Entertain user: --------------------------------------------------------------#
      echo ""
      echo ""
      echo " + Processing ${GRIDID[i]}"
      echo "   - Path:${HERE}/${VERSION}"
      echo "   - Control Files:"
      echo "     * Main:  $(basename ${FILEMAIN})"
      echo "     * Test:  $(basename ${FILETEST})"
      echo "     * Debug: $(basename ${FILEDBUG})"
      echo ""
      #------------------------------------------------------------------------------------#


      #----- Copy the Template to the working directory. ----------------------------------#
      /bin/cp -f ${TEMPLATEMAIN} ${FILEMAIN}
      /bin/cp -f ${TEMPLATETEST} ${FILETEST}
      /bin/cp -f ${TEMPLATEDBUG} ${FILEDBUG}
      #------------------------------------------------------------------------------------#


      #----- Output folders. --------------------------------------------------------------#
      F_MAIN_PATH="${HERE}/${VERSION}/F_main_${GRIDID[i]}"
      S_MAIN_PATH="${HERE}/${VERSION}/S_main_${GRIDID[i]}"
      F_TEST_PATH="${HERE}/${VERSION}/F_test_${GRIDID[i]}"
      S_TEST_PATH="${HERE}/${VERSION}/S_test_${GRIDID[i]}"
      F_DBUG_PATH="${HERE}/${VERSION}/F_dbug_${GRIDID[i]}"
      S_DBUG_PATH="${HERE}/${VERSION}/S_dbug_${GRIDID[i]}"
      #------------------------------------------------------------------------------------#


      #----- Output prefixes. -------------------------------------------------------------#
      F_MAIN_PREF="${F_MAIN_PATH}/main_${GRIDID[i]}"
      S_MAIN_PREF="${S_MAIN_PATH}/main_${GRIDID[i]}"
      F_TEST_PREF="${F_TEST_PATH}/test_${GRIDID[i]}"
      S_TEST_PREF="${S_TEST_PATH}/test_${GRIDID[i]}"
      F_DBUG_PREF="${F_DBUG_PATH}/dbug_${GRIDID[i]}"
      S_DBUG_PREF="${S_DBUG_PATH}/dbug_${GRIDID[i]}"
      #------------------------------------------------------------------------------------#


      #----- Reset and flush the output folders. ------------------------------------------#
      mkdir -p ${F_MAIN_PATH}
      mkdir -p ${F_TEST_PATH}
      mkdir -p ${F_DBUG_PATH}
      mkdir -p ${S_MAIN_PATH}
      mkdir -p ${S_TEST_PATH}
      mkdir -p ${S_DBUG_PATH}
      #------------------------------------------------------------------------------------#



      #----- Decide whether to continue the run or start from beginning. ------------------#
      case ${RUNTYPE} in
      INITIAL)
         #----- Use default initialisation. -----------------------------------------------#
         sed -i '/NL%RUNTYPE/c\   NL%RUNTYPE = '\'${RUNTYPG[i]}\'   ${FILEMAIN}
         sed -i '/NL%RUNTYPE/c\   NL%RUNTYPE = '\'${RUNTYPG[i]}\'   ${FILETEST}
         sed -i '/NL%RUNTYPE/c\   NL%RUNTYPE = '\'${D_RUNTYPG[i]}\' ${FILEDBUG}
         #---------------------------------------------------------------------------------#

         #----- Initial.  All simulations must be submitted. ------------------------------#
         SUBMIT_MAIN="y"
         SUBMIT_TEST="y"
         SUBMIT_DBUG="y"
         #---------------------------------------------------------------------------------#

         ;;
      HISTORY)
         #----- Look for history files (MAIN) then decide how to set up the run. ----------#
         LHIST_MAIN=$(/bin/ls -1 ${S_MAIN_PREF}-S-*h5 2> /dev/null | xargs -n 1 basename   \
                      2> /dev/null | tail -1)
         if [ "${LHIST_MAIN}" != "" ]
         then
            IYEARH_MAIN=$(echo  ${LHIST_MAIN} | cut -c 12-15)
            IMONTHH_MAIN=$(echo ${LHIST_MAIN} | cut -c 17-18)
            IDATEH_MAIN=$(echo  ${LHIST_MAIN} | cut -c 20-21)
            ITIMEH_MAIN=$(echo  ${LHIST_MAIN} | cut -c 23-26)
            sed -i '/NL%RUNTYPE/c\   NL%RUNTYPE = '\''HISTORY'\'      ${FILEMAIN}
            sed -i '/NL%IYEARH/c\   NL%IYEARH = '${IYEARH_MAIN}   ${FILEMAIN}
            sed -i '/NL%IMONTHH/c\   NL%IMONTHH = '${IMONTHH_MAIN} ${FILEMAIN}
            sed -i '/NL%IDATEH/c\   NL%IDATEH = '${IDATEH_MAIN}   ${FILEMAIN}
            sed -i '/NL%ITIMEH/c\   NL%ITIMEH = '${ITIMEH_MAIN}   ${FILEMAIN}
            sed -i '/NL%SFILIN/c\   NL%SFILIN = '\'${S_MAIN_PREF}\'   ${FILEMAIN}
         else
            sed -i '/NL%RUNTYPE/c\   NL%RUNTYPE = '\'${RUNTYPG[i]}\'   ${FILEMAIN}
         fi
         #---------------------------------------------------------------------------------#



         #----- Look for history files (TEST) then decide how to set up the run. ----------#
         LHIST_TEST=$(/bin/ls -1 ${S_TEST_PREF}-S-*h5 2> /dev/null | xargs -n 1 basename   \
                      2> /dev/null | tail -1)
         if [ "${LHIST_TEST}" != "" ]
         then
            IYEARH_TEST=$(echo  ${LHIST_TEST} | cut -c 12-15)
            IMONTHH_TEST=$(echo ${LHIST_TEST} | cut -c 17-18)
            IDATEH_TEST=$(echo  ${LHIST_TEST} | cut -c 20-21)
            ITIMEH_TEST=$(echo  ${LHIST_TEST} | cut -c 23-26)
            sed -i '/NL%RUNTYPE/c\   NL%RUNTYPE = '\''HISTORY'\'      ${FILETEST}
            sed -i '/NL%IYEARH/c\   NL%IYEARH = '${IYEARH_TEST}   ${FILETEST}
            sed -i '/NL%IMONTHH/c\   NL%IMONTHH = '${IMONTHH_TEST} ${FILETEST}
            sed -i '/NL%IDATEH/c\   NL%IDATEH = '${IDATEH_TEST}   ${FILETEST}
            sed -i '/NL%ITIMEH/c\   NL%ITIMEH = '${ITIMEH_TEST}   ${FILETEST}
            sed -i '/NL%SFILIN/c\   NL%SFILIN = '\'${S_TEST_PREF}\'   ${FILETEST}
         else
            sed -i '/NL%RUNTYPE/c\   NL%RUNTYPE = '\'${RUNTYPG[i]}\'   ${FILETEST}
         fi
         #---------------------------------------------------------------------------------#



         #----- Look for history files (DBUG) then decide how to set up the run. ----------#
         LHIST_DBUG=$(/bin/ls -1 ${S_DBUG_PREF}-S-*h5 2> /dev/null | xargs -n 1 basename   \
                      2> /dev/null | tail -1)
         if [ "${LHIST_DBUG}" != "" ]
         then
            IYEARH_DBUG=$(echo  ${LHIST_DBUG} | cut -c 12-15)
            IMONTHH_DBUG=$(echo ${LHIST_DBUG} | cut -c 17-18)
            IDATEH_DBUG=$(echo  ${LHIST_DBUG} | cut -c 20-21)
            ITIMEH_DBUG=$(echo  ${LHIST_DBUG} | cut -c 23-26)
            sed -i '/NL%RUNTYPE/c\   NL%RUNTYPE = '\''HISTORY'\'         ${FILEDBUG}
            sed -i '/NL%IYEARH/c\   NL%IYEARH = '${IYEARH_DBUG}      ${FILEDBUG}
            sed -i '/NL%IMONTHH/c\   NL%IMONTHH = '${IMONTHH_DBUG}    ${FILEDBUG}
            sed -i '/NL%IDATEH/c\   NL%IDATEH = '${IDATEH_DBUG}      ${FILEDBUG}
            sed -i '/NL%ITIMEH/c\   NL%ITIMEH = '${ITIMEH_DBUG}      ${FILEDBUG}
            sed -i '/NL%SFILIN/c\   NL%SFILIN = '\'${S_DBUG_PREF}\'      ${FILEDBUG}
         else
            sed -i '/NL%RUNTYPE/c\   NL%RUNTYPE = '\'${D_RUNTYPG[i]}\'   ${FILEDBUG}
         fi
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Find out whether to submit the run or not.  Default is yes, unless it has   #
         # reached the end.                                                                #
         #---------------------------------------------------------------------------------#
         #----- Main. ---------------------------------------------------------------------#
         STDOUT_MAIN="${HERE}/${VERSION}/main_${GRIDID[i]}.out"
         THE_END_MAIN=$(grep "ED-2.2 execution ends" ${STDOUT_MAIN} 2> /dev/null | wc -l)
         if [ ${THE_END_MAIN} -gt 0 ]
         then 
            SUBMIT_MAIN="n"
         else
            SUBMIT_MAIN="y"
         fi
         #----- Test. ---------------------------------------------------------------------#
         STDOUT_TEST="${HERE}/${VERSION}/test_${GRIDID[i]}.out"
         THE_END_TEST=$(grep "ED-2.2 execution ends" ${STDOUT_TEST} 2> /dev/null | wc -l)
         if [ ${THE_END_TEST} -gt 0 ]
         then 
            SUBMIT_TEST="n"
         else
            SUBMIT_TEST="y"
         fi
         #----- DBUG. ---------------------------------------------------------------------#
         STDOUT_DBUG="${HERE}/${VERSION}/dbug_${GRIDID[i]}.out"
         THE_END_DBUG=$(grep "ED-2.2 execution ends" ${STDOUT_DBUG} 2> /dev/null | wc -l)
         if [ ${THE_END_DBUG} -gt 0 ]
         then 
            SUBMIT_DBUG="n"
         else
            SUBMIT_DBUG="y"
         fi
         #---------------------------------------------------------------------------------#
         ;;
      esac
      #------------------------------------------------------------------------------------#



      #----- Modify the ED2IN IED_INIT_MODE. ----------------------------------------------#
      sed -i '/NL%IED_INIT_MODE/c\   NL%IED_INIT_MODE = '${INITMDG[i]}   ${FILEMAIN}
      sed -i '/NL%IED_INIT_MODE/c\   NL%IED_INIT_MODE = '${INITMDG[i]}   ${FILETEST}
      sed -i '/NL%IED_INIT_MODE/c\   NL%IED_INIT_MODE = '${D_INITMDG[i]} ${FILEDBUG}
      #------------------------------------------------------------------------------------#


      #----- Modify the ED2IN start years. ------------------------------------------------#
      sed -i '/NL%IYEARA/c\   NL%IYEARA   = '${IYEARAG[i]}   ${FILEMAIN}
      sed -i '/NL%IYEARA/c\   NL%IYEARA   = '${IYEARAG[i]}   ${FILETEST}
      sed -i '/NL%IYEARA/c\   NL%IYEARA   = '${D_IYEARAG[i]} ${FILEDBUG}
      #------------------------------------------------------------------------------------#


      #----- Modify the ED2IN end years. --------------------------------------------------#
      sed -i '/NL%IYEARZ/c\   NL%IYEARZ   = '${IYEARZG[i]}   ${FILEMAIN}
      sed -i '/NL%IYEARZ/c\   NL%IYEARZ   = '${IYEARZG[i]}   ${FILETEST}
      sed -i '/NL%IYEARZ/c\   NL%IYEARZ   = '${D_IYEARZG[i]} ${FILEDBUG}
      #------------------------------------------------------------------------------------#


      #----- Modify the ED2IN files to point to the desired output directories. -----------#
      sed -i '/NL%FFILOUT/c\   NL%FFILOUT = '\'${F_MAIN_PREF}\' ${FILEMAIN}
      sed -i '/NL%SFILOUT/c\   NL%SFILOUT = '\'${S_MAIN_PREF}\' ${FILEMAIN}
      sed -i '/NL%FFILOUT/c\   NL%FFILOUT = '\'${F_TEST_PREF}\' ${FILETEST}
      sed -i '/NL%SFILOUT/c\   NL%SFILOUT = '\'${S_TEST_PREF}\' ${FILETEST}
      sed -i '/NL%FFILOUT/c\   NL%FFILOUT = '\'${F_DBUG_PREF}\' ${FILEDBUG}
      sed -i '/NL%SFILOUT/c\   NL%SFILOUT = '\'${S_DBUG_PREF}\' ${FILEDBUG}
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #   Make sure medium and long tests won't write fast or daily files.                 #
      #------------------------------------------------------------------------------------#
      case ${TESTTYPE} in
      medium|long)
         sed -i '/NL%IFOUTPUT/c\   NL%IFOUTPUT  = 0' ${FILETEST}
         sed -i '/NL%IDOUTPUT/c\   NL%IDOUTPUT  = 0' ${FILETEST}
         sed -i '/NL%IFOUTPUT/c\   NL%IFOUTPUT  = 0' ${FILEMAIN}
         sed -i '/NL%IDOUTPUT/c\   NL%IDOUTPUT  = 0' ${FILEMAIN}
         ;;
      esac
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #    Prepare job (MAIN)
      #------------------------------------------------------------------------------------#
      case ${SUBMIT_MAIN} in
      y|Y)
         jobout=${HERE}/${VERSION}/main_${GRIDID[i]}.out
         joberr=${HERE}/${VERSION}/main_${GRIDID[i]}.err
         jobname=${VERSION}_${GRIDID[i]}_main
         jobopts="-t ${POI_TIME} --mem-per-cpu=${GRIDMEM[i]} --cpus-per-task=1"
         jobopts="${jobopts} -p ${GRIDQ[i]} -n ${GRIDCPU[i]}"
         jobwrap=". ${HOME}/.bashrc; cd ${HERE}/${VERSION}"
         jobwrap="${jobwrap}; export OMP_NUM_THREADS=1"
         jobwrap="${jobwrap}; mpirun -np ${GRIDCPU[i]} ${LNK_MAIN_EXE} -f ${FILEMAIN}"
         jobwrap="\"(${jobwrap})\""
         jobcomm="sbatch -o ${jobout} -e ${joberr} -J ${jobname}"
         jobcomm="${jobcomm} ${jobopts} --wrap=${jobwrap}"
         echo ${jobcomm} >> ${subbatch}
         ;;
      *)
         echo "     * Main has already reached the end and won't be submitted."
         ;;
      esac
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #    Prepare job (TEST)
      #------------------------------------------------------------------------------------#
      case ${SUBMIT_TEST} in
      y|Y)
         jobout=${HERE}/${VERSION}/test_${GRIDID[i]}.out
         joberr=${HERE}/${VERSION}/test_${GRIDID[i]}.err
         jobname=${VERSION}_${GRIDID[i]}_test
         jobopts="-t ${POI_TIME} --mem-per-cpu=${GRIDMEM[i]} --cpus-per-task=1"
         jobopts="${jobopts} -p ${GRIDQ[i]} -n ${GRIDCPU[i]}"
         jobwrap=". ${HOME}/.bashrc; cd ${HERE}/${VERSION}"
         jobwrap="${jobwrap}; export OMP_NUM_THREADS=1"
         jobwrap="${jobwrap}; mpirun -np ${GRIDCPU[i]} ${LNK_TEST_EXE} -f ${FILETEST}"
         jobwrap="\"(${jobwrap})\""
         jobcomm="sbatch -o ${jobout} -e ${joberr} -J ${jobname}"
         jobcomm="${jobcomm} ${jobopts} --wrap=${jobwrap}"
         echo ${jobcomm} >> ${subbatch}
         ;;
      *)
         echo "     * Test has already reached the end and won't be submitted."
         ;;
      esac
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #    Prepare job (DBUG)
      #------------------------------------------------------------------------------------#
      case ${SUBMIT_DBUG} in
      y|Y)
         jobout=${HERE}/${VERSION}/dbug_${GRIDID[i]}.out
         joberr=${HERE}/${VERSION}/dbug_${GRIDID[i]}.err
         jobname=${VERSION}_${GRIDID[i]}_dbug
         jobopts="-t ${POI_TIME} --mem-per-cpu=${GRIDMEM[i]} --cpus-per-task=1"
         jobopts="${jobopts} -p ${GRIDQ[i]} -n ${GRIDCPU[i]}"
         jobwrap=". ${HOME}/.bashrc; cd ${HERE}/${VERSION}"
         jobwrap="${jobwrap}; export OMP_NUM_THREADS=1"
         jobwrap="${jobwrap}; mpirun -np ${GRIDCPU[i]} ${LNK_DBUG_EXE} -f ${FILEDBUG}"
         jobwrap="\"(${jobwrap})\""
         jobcomm="sbatch -o ${jobout} -e ${joberr} -J ${jobname}"
         jobcomm="${jobcomm} ${jobopts} --wrap=${jobwrap}"
         echo ${jobcomm} >> ${subbatch}
         ;;
      *)
         echo "     * Dbug has already reached the end and won't be submitted."
         ;;
      esac
      #------------------------------------------------------------------------------------#

      ;;
   esac
   #---------------------------------------------------------------------------------------#
done
#------------------------------------------------------------------------------------------#



#----- Check that the model is running. ---------------------------------------------------#
echo "sleep 10"                                                              >> ${subbatch}
echo "echo \"\""                                                             >> ${subbatch}
echo "echo \"\""                                                             >> ${subbatch}
echo "echo \"============================================================\"" >> ${subbatch}
echo "echo \" Job list - SLURM \""                                           >> ${subbatch}
echo "squeue -u ${USER}"                                                     >> ${subbatch}
echo "echo \"============================================================\"" >> ${subbatch}
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Change permissions of subbatch.                                                      #
#------------------------------------------------------------------------------------------#
chmod u+x ${subbatch}
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#      Check whether to submit jobs or not.                                                #
#------------------------------------------------------------------------------------------#
case ${SUBMIT_JOBS} in
y|Y)
   #----- Submit the jobs. ----------------------------------------------------------------#
   echo " + Submitting jobs."
   (cd ${HERE}/${VERSION}; ${subbatch})
   #---------------------------------------------------------------------------------------#
   ;;
n|N)
   #----- Tell the user where the script is located. --------------------------------------#
   echo " + Script with EDTS jobs available at ${subbatch}."
   echo " + The jobs haven't been submitted, you must run the script."
   #---------------------------------------------------------------------------------------#
   ;;
*)
esac
#------------------------------------------------------------------------------------------#
