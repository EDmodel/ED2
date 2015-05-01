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
#             "rapid" -- 2-3 years of simulation                                           #
#             "long"  -- ~ 300 years of simulation.   The long tests will take a really    #
#                        really long time, beware.                                         #
#------------------------------------------------------------------------------------------#
TESTTYPE="rapid"
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#  DATAPATH -- Location with all data sets.                                                #
#------------------------------------------------------------------------------------------#
DATAPATH="/n/gstore/Labs/moorcroft_lab_protected/ed_mainline/test_suite/datasets"
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
USE_HIP="y"     # Petrolina High Frequency Detailed Short POI
USE_HIM="y"     # Manaus High Frequency Detailed Short POI
USE_RJG="y"     # Gridded 12x12 simulation centred on Reserva Jaru
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#     The following flags are the queue names for each simulation.                         #
#   Main, test, and dbug should all use the same queue.                                    #
#------------------------------------------------------------------------------------------#
Q_M34="general"     # Manaus K34 POI
Q_S67="general"     # Santarem km 67 POI
Q_HAR="general"     # Harvard Forest POI
Q_PDG="general"     # Pe-de-Gigante POI
Q_TON="general"     # Tonzi POI
Q_CAX="general"     # Caxiuana POI
Q_TNF="general"     # Tapajos National Forest POI
Q_ATA="general"     # Atacama Desert POI
Q_PET="general"     # Petrolina POI
Q_HIP="general"     # Petrolina High Frequency Detailed Short POI
Q_HIM="general"     # Manaus High Frequency Detailed Short POI
Q_RJG="general"     # Gridded 12x12 simulation centred on Reserva Jaru
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
# NPROC -- How many cores do you want to use for the gridded simulations?                  #
#------------------------------------------------------------------------------------------#
NPROC=27
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#     The following variables are the time required by each simulation, hh:mm:ss.          #
#------------------------------------------------------------------------------------------#
POI_FAST_TIME="168:00:00"
POI_LONG_TIME="168:00:00"
GRID_FAST_TIME="168:00:00"
GRID_LONG_TIME="312:00:00"
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#   The amount of memory per cpu needed (POI and gridded runs).                            #
#------------------------------------------------------------------------------------------#
POI_MEMORY=8192
GRID_MEMORY=4096
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#  TEST_DESCRIPTION -- A short but comprehensive explanation of the tests.  Explain what   #
#                      the commits had involved.                                           #
#------------------------------------------------------------------------------------------#
TEST_DESCRIPTION="1. Changed minimum height for reproduction based on BCI measurements (tropical trees only).  2. Changes in minimum reproduction size, so seed_rain works.  3. Changed compilation instructions so it can compile most files with -O3 but uses -O2 for files that would otherwise take days to compile with ifort 13.  4.  New dist_type categories, which now distinguishes tree fall, logging, fires, abandonment and forest plantations.  5.  Changed forestry.f90, this file computes the disturbance rates but lets disturbance.f90 to apply them do cpoly%disturbance_rates is now correctly updated.   6.  Several minor bug fixes."
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
declare -a USE_SITE=( ${USE_M34} ${USE_S67} ${USE_HAR} ${USE_PDG} ${USE_TON} \
                      ${USE_CAX} ${USE_TNF} ${USE_ATA} ${USE_PET} )
declare -a SITEID=(m34 s67 har pdg ton cax tnf ata pet)
declare -a SITEPFX=(M34 S67 HAR PDG TON CAX TNF ATA PET)
declare -a SITEQ=( ${Q_M34} ${Q_S67} ${Q_HAR} ${Q_PDG} ${Q_TON} \
                   ${Q_CAX} ${Q_TNF} ${Q_ATA} ${Q_PET})
#------------------------------------------------------------------------------------------#



#----- POI debug time. --------------------------------------------------------------------#
declare -a D_IYEARAS=(1500 1500 2007 1500 2000 2000 2002 1500 2005)
declare -a D_IYEARZS=(1501 1501 2008 1501 2001 2001 2002 1501 2006)
declare -a D_INITMDS=(5    0    6    0    5    0    5    0    6   )
declare -a D_RUNTYPS=( INITIAL INITIAL INITIAL INITIAL INITIAL \
                       INITIAL INITIAL INITIAL INITIAL )
#------------------------------------------------------------------------------------------#



#----- High-frequency detailed runs (HI-PET and HI-M34). ----------------------------------#
declare -a USE_HIFR=(${USE_HIP} ${USE_HIM})
declare -a HIFRID=(hip him)
declare -a HIFRPFX=(HIP HIM)
declare -a HIFRQ=(${Q_HIP} ${Q_HIM})
declare -a IDATEAH=(21 01)
declare -a IDATEZH=(28 08)
declare -a INITMDH=(6  5)
declare -a RUNTYPH=(INITIAL INITIAL)
#------------------------------------------------------------------------------------------#


#----- Gridded runs. ----------------------------------------------------------------------#
declare -a USE_GRID=($USE_RJG)
declare -a GRIDID=(rjg)
declare -a GRIDPFX=(RJG)
declare -a GRIDQ=($Q_RJG)
declare -a GRIDPROC=($NPROC)
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
   declare -a IYEARAS=(1500 1500 2007 1500 2000 2000 2002 1500 2005)
   declare -a IYEARZS=(1502 1502 2009 1502 2002 2002 2004 1502 2007)
   declare -a INITMDS=(5    0    6    0    6    0    5    0    6   )
   declare -a RUNTYPS=(INITIAL INITIAL INITIAL INITIAL INITIAL \
                       INITIAL INITIAL INITIAL INITIAL)
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

   #----- Select fast time requests. ------------------------------------------------------#
   POI_TIME=${POI_FAST_TIME}
   GRID_TIME=${GRID_FAST_TIME}
   #---------------------------------------------------------------------------------------#
   ;;

long)

   echo " - Performing long tests (300 years for POI, 35 year for grid)"

   #----- POI tests will run for 300 years. -----------------------------------------------#
   declare -a IYEARAS=(1500 1500 1500 1500 1500 1500 1500 1500 1500)
   declare -a IYEARZS=(1800 1800 1800 1800 1800 1800 1800 1800 1800)
   declare -a INITMDS=(0    0    0    0    0    0    0    0    0   )
   declare -a RUNTYPS=( \
       INITIAL INITIAL INITIAL \
       INITIAL INITIAL INITIAL \
       INITIAL INITIAL INITIAL)
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


   #----- Select fast time requests. ------------------------------------------------------#
   POI_TIME=${POI_LONG_TIME}
   GRID_TIME=${GRID_LONG_TIME}
   #---------------------------------------------------------------------------------------#
   ;;
*)
   #----- I don't know what to do... ------------------------------------------------------#
   echo " !Invalid TESTTYPE: ${TESTTYPE}.  Please specify either 'rapid' or 'long'."
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
if [ -s ${HERE}/${VERSION} ]
then
   echo ""
   echo "     Delete previous data and generate pristine directory: ${VERSION}? [y|n]"
   echo " In case you answer y ALL previous data will be lost!"
   read flusher
   flusher=$(echo ${flusher} | tr '[:upper:]' '[:lower:]')
   case ${flusher} in
   y|Y)
      #----- Be nice and give a chance to the user to change their minds. -----------------#
      echo "Fine, but if you regret later don't say I didn't warn you..."
      echo "I am giving you a few seconds to kill this script in case you aren't sure..."
      delfun=16
      while [ ${delfun} -gt 1 ]
      do
         let delfun=${delfun}-1
         echo "  - Deletion will begin in '${delfun}' seconds..."
         sleep 1
      done
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Purge directory then recreate it.  Use find to erase it because it's faster.   #
      #------------------------------------------------------------------------------------#
      find ${HERE}/${VERSION} -print -delete
      #------------------------------------------------------------------------------------#

      ;;
   n|N)
      #------------------------------------------------------------------------------------#
      #     Overwrite stuff.                                                               #
      #------------------------------------------------------------------------------------#
      echo " Keeping the current directory.  Files may be overwritten, beware..."
      ;;
      #------------------------------------------------------------------------------------#
   *)
      #----- User doesn't want to continue, stop script. ----------------------------------#
      echo " Invalid option. Hint: 'y' is short of yes, and 'n' is short of no."
      echo " 'Maybe' and 'whatever' have not been implemented yet."
      exit 0
      #------------------------------------------------------------------------------------#
      ;;
   esac
   #---------------------------------------------------------------------------------------#
fi
#------------------------------------------------------------------------------------------#


#----- Start a new directory, create links and copy xml files. ----------------------------#
mkdir -p ${HERE}/${VERSION}
LNK_MAIN_EXE="${HERE}/${VERSION}/${MAIN_EXE}"
LNK_TEST_EXE="${HERE}/${VERSION}/${TEST_EXE}"
LNK_DBUG_EXE="${HERE}/${VERSION}/${DBUG_EXE}"
ln -sv ${DATAPATH}            ${HERE}/${VERSION}/edts_datasets
ln -sv ${MAIN_EXE_PATH}       ${LNK_MAIN_EXE}
ln -sv ${TEST_EXE_PATH}       ${LNK_TEST_EXE}
ln -sv ${DBUG_EXE_PATH}       ${LNK_DBUG_EXE}
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



      #----- Modify the ED2IN RUNTYPE. ----------------------------------------------------#
      sed -i '/NL%RUNTYPE/c\   NL%RUNTYPE = '\'${RUNTYPS[i]}\'   ${FILEMAIN}
      sed -i '/NL%RUNTYPE/c\   NL%RUNTYPE = '\'${RUNTYPS[i]}\'   ${FILETEST}
      sed -i '/NL%RUNTYPE/c\   NL%RUNTYPE = '\'${D_RUNTYPS[i]}\' ${FILEDBUG}
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
      sed -i '/NL%FFILOUT/c\   NL%FFILOUT = '\'${HERE}'/'${VERSION}'/F_main_'${SITEID[i]}'/main_'${SITEID[i]}''\' ${FILEMAIN}
      sed -i '/NL%SFILOUT/c\   NL%SFILOUT = '\'${HERE}'/'${VERSION}'/S_main_'${SITEID[i]}'/main_'${SITEID[i]}''\' ${FILEMAIN}
      sed -i '/NL%FFILOUT/c\   NL%FFILOUT = '\'${HERE}'/'${VERSION}'/F_test_'${SITEID[i]}'/test_'${SITEID[i]}''\' ${FILETEST}
      sed -i '/NL%SFILOUT/c\   NL%SFILOUT = '\'${HERE}'/'${VERSION}'/S_test_'${SITEID[i]}'/test_'${SITEID[i]}''\' ${FILETEST}
      sed -i '/NL%FFILOUT/c\   NL%FFILOUT = '\'${HERE}'/'${VERSION}'/F_dbug_'${SITEID[i]}'/dbug_'${SITEID[i]}''\' ${FILEDBUG}
      sed -i '/NL%SFILOUT/c\   NL%SFILOUT = '\'${HERE}'/'${VERSION}'/S_dbug_'${SITEID[i]}'/dbug_'${SITEID[i]}''\' ${FILEDBUG}
      #------------------------------------------------------------------------------------#



      #----- Modify the state file output frequency. --------------------------------------#
      sed -i '/NL%UNITSTATE/c\   NL%UNITSTATE  = '${UNITSTATE} ${FILEMAIN}
      sed -i '/NL%UNITSTATE/c\   NL%UNITSTATE  = '${UNITSTATE} ${FILETEST}
      sed -i '/NL%UNITSTATE/c\   NL%UNITSTATE  = '${UNITSTATE} ${FILEDBUG}
      #------------------------------------------------------------------------------------#


      #----- Reset and flush the output folders. ------------------------------------------#
      mkdir -p ${HERE}/${VERSION}/F_main_${SITEID[i]}
      mkdir -p ${HERE}/${VERSION}/F_test_${SITEID[i]}
      mkdir -p ${HERE}/${VERSION}/F_dbug_${SITEID[i]}
      mkdir -p ${HERE}/${VERSION}/S_main_${SITEID[i]}
      mkdir -p ${HERE}/${VERSION}/S_test_${SITEID[i]}
      mkdir -p ${HERE}/${VERSION}/S_dbug_${SITEID[i]}
      #------------------------------------------------------------------------------------#


      #----- Update count. ----------------------------------------------------------------#
      let maincount=${tcount}+0
      let testcount=${tcount}+1
      let dbugcount=${tcount}+2
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #    Prepare job (MAIN)
      #------------------------------------------------------------------------------------#
      jobout=${HERE}/${VERSION}/main_${SITEID[i]}.out
      joberr=${HERE}/${VERSION}/main_${SITEID[i]}.err
      jobname=${VERSION}_${SITEID[i]}_main
      jobopts="-t ${POI_TIME} --mem-per-cpu=${POI_MEMORY} -p ${SITEQ[i]} -n 1"
      jobwrap="\"(. ${HOME}/.bashrc; cd ${HERE}/${VERSION}; ${LNK_MAIN_EXE} -f ${FILEMAIN})\""
      jobcomm="sbatch -o ${jobout} -e ${joberr} -J ${jobname}"
      jobcomm="${jobcomm} ${jobopts} --wrap=${jobwrap}"
      echo ${jobcomm} >> ${subbatch}
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #    Prepare job (TEST)
      #------------------------------------------------------------------------------------#
      jobout=${HERE}/${VERSION}/test_${SITEID[i]}.out
      joberr=${HERE}/${VERSION}/test_${SITEID[i]}.err
      jobname=${VERSION}_${SITEID[i]}_test
      jobopts="-t ${POI_TIME} --mem-per-cpu=${POI_MEMORY} -p ${SITEQ[i]} -n 1"
      jobwrap="\"(. ${HOME}/.bashrc; cd ${HERE}/${VERSION}; ${LNK_TEST_EXE} -f ${FILETEST})\""
      jobcomm="sbatch -o ${jobout} -e ${joberr} -J ${jobname}"
      jobcomm="${jobcomm} ${jobopts} --wrap=${jobwrap}"
      echo ${jobcomm} >> ${subbatch}
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #    Prepare job (DBUG)
      #------------------------------------------------------------------------------------#
      jobout=${HERE}/${VERSION}/dbug_${SITEID[i]}.out
      joberr=${HERE}/${VERSION}/dbug_${SITEID[i]}.out
      jobname=${VERSION}_${SITEID[i]}_dbug
      jobopts="-t ${POI_TIME} --mem-per-cpu=${POI_MEMORY} -p ${SITEQ[i]} -n 1"
      jobwrap="\"(. ${HOME}/.bashrc; cd ${HERE}/${VERSION}; ${LNK_DBUG_EXE} -f ${FILEDBUG})\""
      jobcomm="sbatch -o ${jobout} -e ${joberr} -J ${jobname}"
      jobcomm="${jobcomm} ${jobopts} --wrap=${jobwrap}"
      echo ${jobcomm} >> ${subbatch}
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



      #----- Modify the ED2IN RUNTYPE. ----------------------------------------------------#
      sed -i '/NL%RUNTYPE/c\   NL%RUNTYPE = '\'${RUNTYPH[i]}\' ${FILEMAIN}
      sed -i '/NL%RUNTYPE/c\   NL%RUNTYPE = '\'${RUNTYPH[i]}\' ${FILETEST}
      sed -i '/NL%RUNTYPE/c\   NL%RUNTYPE = '\'${RUNTYPH[i]}\' ${FILEDBUG}
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
      sed -i '/NL%FFILOUT/c\   NL%FFILOUT = '\'${HERE}'/'${VERSION}'/F_main_'${HIFRID[i]}'/main_'${HIFRID[i]}''\' ${FILEMAIN}
      sed -i '/NL%SFILOUT/c\   NL%SFILOUT = '\'${HERE}'/'${VERSION}'/S_main_'${HIFRID[i]}'/main_'${HIFRID[i]}''\' ${FILEMAIN}
      sed -i '/NL%FFILOUT/c\   NL%FFILOUT = '\'${HERE}'/'${VERSION}'/F_test_'${HIFRID[i]}'/test_'${HIFRID[i]}''\' ${FILETEST}
      sed -i '/NL%SFILOUT/c\   NL%SFILOUT = '\'${HERE}'/'${VERSION}'/S_test_'${HIFRID[i]}'/test_'${HIFRID[i]}''\' ${FILETEST}
      sed -i '/NL%FFILOUT/c\   NL%FFILOUT = '\'${HERE}'/'${VERSION}'/F_dbug_'${HIFRID[i]}'/dbug_'${HIFRID[i]}''\' ${FILEDBUG}
      sed -i '/NL%SFILOUT/c\   NL%SFILOUT = '\'${HERE}'/'${VERSION}'/S_dbug_'${HIFRID[i]}'/dbug_'${HIFRID[i]}''\' ${FILEDBUG}
      #------------------------------------------------------------------------------------#


      #----- Reset and flush the output folders. ------------------------------------------#
      mkdir -p ${HERE}/${VERSION}/F_main_${HIFRID[i]}
      mkdir -p ${HERE}/${VERSION}/F_test_${HIFRID[i]}
      mkdir -p ${HERE}/${VERSION}/F_dbug_${HIFRID[i]}
      mkdir -p ${HERE}/${VERSION}/S_main_${HIFRID[i]}
      mkdir -p ${HERE}/${VERSION}/S_test_${HIFRID[i]}
      mkdir -p ${HERE}/${VERSION}/S_dbug_${HIFRID[i]}
      #------------------------------------------------------------------------------------#


      #----- Update count. ----------------------------------------------------------------#
      let maincount=${tcount}+0
      let testcount=${tcount}+1
      let dbugcount=${tcount}+2
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #    Prepare job (MAIN)
      #------------------------------------------------------------------------------------#
      jobout=${HERE}/${VERSION}/main_${HIFRID[i]}.out
      joberr=${HERE}/${VERSION}/main_${HIFRID[i]}.out
      jobname=${VERSION}_${HIFRID[i]}_main
      jobopts="-t ${POI_TIME} --mem-per-cpu=${POI_MEMORY} -p ${HIFRQ[i]} -n 1"
      jobwrap="\"(. ${HOME}/.bashrc; cd ${HERE}/${VERSION}; ${LNK_MAIN_EXE} -f ${FILEMAIN})\""
      jobcomm="sbatch -o ${jobout} -e ${joberr} -J ${jobname}"
      jobcomm="${jobcomm} ${jobopts} --wrap=${jobwrap}"
      echo ${jobcomm} >> ${subbatch}
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #    Prepare job (TEST)
      #------------------------------------------------------------------------------------#
      jobout=${HERE}/${VERSION}/test_${HIFRID[i]}.out
      joberr=${HERE}/${VERSION}/test_${HIFRID[i]}.out
      jobname=${VERSION}_${HIFRID[i]}_test
      jobopts="-t ${POI_TIME} --mem-per-cpu=${POI_MEMORY} -p ${HIFRQ[i]} -n 1"
      jobwrap="\"(. ${HOME}/.bashrc; cd ${HERE}/${VERSION}; ${LNK_TEST_EXE} -f ${FILETEST})\""
      jobcomm="sbatch -o ${jobout} -e ${joberr} -J ${jobname}"
      jobcomm="${jobcomm} ${jobopts} --wrap=${jobwrap}"
      echo ${jobcomm} >> ${subbatch}
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #    Prepare job (DBUG)
      #------------------------------------------------------------------------------------#
      jobout=${HERE}/${VERSION}/dbug_${HIFRID[i]}.out
      joberr=${HERE}/${VERSION}/dbug_${HIFRID[i]}.out
      jobname=${VERSION}_${HIFRID[i]}_dbug
      jobopts="-t ${POI_TIME} --mem-per-cpu=${POI_MEMORY} -p ${HIFRQ[i]} -n 1"
      jobwrap="\"(. ${HOME}/.bashrc; cd ${HERE}/${VERSION}; ${LNK_DBUG_EXE} -f ${FILEDBUG})\""
      jobcomm="sbatch -o ${jobout} -e ${joberr} -J ${jobname}"
      jobcomm="${jobcomm} ${jobopts} --wrap=${jobwrap}"
      echo ${jobcomm} >> ${subbatch}
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



      #----- Modify the ED2IN RUNTYPE. ----------------------------------------------------#
      sed -i '/NL%RUNTYPE/c\   NL%RUNTYPE = '\'${RUNTYPG[i]}\'   ${FILEMAIN}
      sed -i '/NL%RUNTYPE/c\   NL%RUNTYPE = '\'${RUNTYPG[i]}\'   ${FILETEST}
      sed -i '/NL%RUNTYPE/c\   NL%RUNTYPE = '\'${D_RUNTYPG[i]}\' ${FILEDBUG}
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
      sed -i '/NL%FFILOUT/c\   NL%FFILOUT = '\'${HERE}'/'${VERSION}'/F_main_'${GRIDID[i]}'/main_'${GRIDID[i]}''\' ${FILEMAIN}
      sed -i '/NL%SFILOUT/c\   NL%SFILOUT = '\'${HERE}'/'${VERSION}'/S_main_'${GRIDID[i]}'/main_'${GRIDID[i]}''\' ${FILEMAIN}
      sed -i '/NL%FFILOUT/c\   NL%FFILOUT = '\'${HERE}'/'${VERSION}'/F_test_'${GRIDID[i]}'/test_'${GRIDID[i]}''\' ${FILETEST}
      sed -i '/NL%SFILOUT/c\   NL%SFILOUT = '\'${HERE}'/'${VERSION}'/S_test_'${GRIDID[i]}'/test_'${GRIDID[i]}''\' ${FILETEST}
      sed -i '/NL%FFILOUT/c\   NL%FFILOUT = '\'${HERE}'/'${VERSION}'/F_dbug_'${GRIDID[i]}'/dbug_'${GRIDID[i]}''\' ${FILEDBUG}
      sed -i '/NL%SFILOUT/c\   NL%SFILOUT = '\'${HERE}'/'${VERSION}'/S_dbug_'${GRIDID[i]}'/dbug_'${GRIDID[i]}''\' ${FILEDBUG}
      #------------------------------------------------------------------------------------#


      #----- Reset and flush the output folders. ------------------------------------------#
      mkdir -p ${HERE}/${VERSION}/F_main_${GRIDID[i]}
      mkdir -p ${HERE}/${VERSION}/F_test_${GRIDID[i]}
      mkdir -p ${HERE}/${VERSION}/F_dbug_${GRIDID[i]}
      mkdir -p ${HERE}/${VERSION}/S_main_${GRIDID[i]}
      mkdir -p ${HERE}/${VERSION}/S_test_${GRIDID[i]}
      mkdir -p ${HERE}/${VERSION}/S_dbug_${GRIDID[i]}
      #------------------------------------------------------------------------------------#


      #----- Update count. ----------------------------------------------------------------#
      let maincount=${tcount}+0
      let testcount=${tcount}+1
      let dbugcount=${tcount}+2
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #    Prepare job (MAIN)
      #------------------------------------------------------------------------------------#
      jobout=${HERE}/${VERSION}/main_${GRIDID[i]}.out
      joberr=${HERE}/${VERSION}/main_${GRIDID[i]}.out
      jobname=${VERSION}_${GRIDID[i]}_main
      jobopts="-t ${GRID_TIME} --mem-per-cpu=${GRID_MEMORY} -p ${GRIDQ[i]} -n ${GRIDPROC[i]}"
      jobwrap="\"(. ${HOME}/.bashrc; cd ${HERE}/${VERSION}; mpirun -np ${GRIDPROC[i]} ${LNK_MAIN_EXE} -f ${FILEMAIN})\""
      jobcomm="sbatch -o ${jobout} -e ${joberr} -J ${jobname}"
      jobcomm="${jobcomm} ${jobopts} --wrap=${jobwrap}"
      echo ${jobcomm} >> ${subbatch}
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #    Prepare job (TEST)
      #------------------------------------------------------------------------------------#
      jobout=${HERE}/${VERSION}/test_${GRIDID[i]}.out
      joberr=${HERE}/${VERSION}/test_${GRIDID[i]}.out
      jobname=${VERSION}_${GRIDID[i]}_test
      jobopts="-t ${GRID_TIME} --mem-per-cpu=${GRID_MEMORY} -p ${GRIDQ[i]} -n ${GRIDPROC[i]}"
      jobwrap="\"(. ${HOME}/.bashrc; cd ${HERE}/${VERSION}; mpirun -np ${GRIDPROC[i]} ${LNK_TEST_EXE} -f ${FILETEST})\""
      jobcomm="sbatch -o ${jobout} -e ${joberr} -J ${jobname}"
      jobcomm="${jobcomm} ${jobopts} --wrap=${jobwrap}"
      echo ${jobcomm} >> ${subbatch}
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #    Prepare job (DBUG)
      #------------------------------------------------------------------------------------#
      jobout=${HERE}/${VERSION}/dbug_${GRIDID[i]}.out
      joberr=${HERE}/${VERSION}/dbug_${GRIDID[i]}.out
      jobname=${VERSION}_${GRIDID[i]}_dbug
      jobopts="-t ${GRID_TIME} --mem-per-cpu=${GRID_MEMORY} -p ${GRIDQ[i]} -n ${GRIDPROC[i]}"
      jobwrap="\"(. ${HOME}/.bashrc; cd ${HERE}/${VERSION}; mpirun -np ${GRIDPROC[i]} ${LNK_DBUG_EXE} -f ${FILEDBUG})\""
      jobcomm="sbatch -o ${jobout} -e ${joberr} -J ${jobname}"
      jobcomm="${jobcomm} ${jobopts} --wrap=${jobwrap}"
      echo ${jobcomm} >> ${subbatch}
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
