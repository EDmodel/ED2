#!/bin/bash


#===============================================================================
#                    User control variables
#===============================================================================

# This is the unique identifier tag for this test environment
# It should indicate the revision number branched from, the users initials
# and a test number.  For instance if you branched from r81, your initials
# are rk and you tried 4 iterations of this test before your commits
# were verified, you would use:

VERSION="r85v2ghub"

# FILE PATHS TO YOUR THREE EXECUTABLES

TEST_EXE_PATH=/global/u2/r/rgknox/Models/ghub_r85v2_rgknox/ED2/ED/build/ed_2.1-opt
DBUG_EXE_PATH=/global/u2/r/rgknox/Models/ghub_r85v2_rgknox/ED2/ED/build/ed_2.1-dbg
MAIN_EXE_PATH=/global/u2/r/rgknox/Models/ghub_r84_master/ED/build/ed_2.1-opt

# Provide the path where the test-suit driver archive is stored

DATAPATH='/global/scratch2/sd/rgknox/edts_data'

# Decide on a "rapid" (2 year) or "long" (300 year)
# THE long tests will take a really really long time, beware

TESTTYPE="rapid"

# Decide on which sites to use.  You should ultimately
# have results for all of them, but this is helpfull if some
# sites give you problems after the first try
# Y for Yes and N for No

USE_M34=N     #Manaus km 34 SOI
USE_S67=N     #Santarem km 67 SOI
USE_HAR=N     #Harvard forest SOI
USE_PDG=N     #Pe de Gigante SOI 
USE_TON=N     #Tonzi SOI
USE_CAX=N     #Caxuana SOI
USE_TNF=N     #Tapajos National Forest SOI
USE_ATA=N     #Atacama Desert SOI
USE_PET=N     #Petrolina SOI
USE_HIP=N     #Petrolina High Frequency Detailed Short SOI
USE_HIM=N     #Manaus High Frequency Detailed Short SOI
USE_RJG=Y     #Gridded 12x12 simulation centerd on Reserva Jaru

# How many cores do you want to use for the gridded simulations
# Currently there are 3 (RJG-MAIN RJG-TEST and RJG-DBUG)

NNODES=3
PPN=8


# WHICH QUEUES WOULD USE LIKE TO USE?
# MAIN-TEST-DBUG MUST ALL USE THE SAME QUEUE

Q_M34=regular     #Manaus km 34 SOI
Q_S67=regular     #Santarem km 67 SOI
Q_HAR=regular     #Harvard forest SOI
Q_PDG=regular     #Pe de Gigante SOI
Q_TON=regular     #Tonzi SOI
Q_CAX=regular     #Caxuana SOI
Q_TNF=regular     #Tapajos National Forest SOI
Q_ATA=regular     #Atacama Desert SOI
Q_PET=regular     #Petrolina SOI
Q_HIP=regular     #Petrolina High Frequency Detailed Short SOI
Q_HIM=regular     #Manaus High Frequency Detailed Short SOI
Q_RJG=regular     #Gridded 12x12 simulation centerd on Reserva Jaru

# Give an explanation of the tests.  Explain what the commits had involved.

TEST_DESCRIPTION="This is a test of the final commits from Odyssey, which pretty much amounts
to Marcos changes with a couple of minor tweaks from Ryan."

# The identifier may had indicated which version you branched from, but indicate it here
# also
    
VERSION_BRANCHED_FROM='r83'

# Who is running this test?
TESTER_NAME='Ryan Knox'

# Who was the developer(s) that actually made the changes to the code that is being tested?
COMMITTER_NAME='Marcos Longo'



#===============================================================================


echo ""
echo ""
echo "========================================================================="
echo "                 Starting the EDM Dev Test Suit  (EDTS)                  "
echo "========================================================================="
echo ""

# DEFINE SOME RUNTIME VARIABLES FOR SOI's

declare -a USE_SITE=( \
    $USE_M34 $USE_S67 $USE_HAR \
    $USE_PDG $USE_TON $USE_CAX \
    $USE_TNF $USE_ATA $USE_PET)
declare -a SITEID=(m34 s67 har pdg ton cax tnf ata pet)
declare -a SITEPFX=(M34 S67 HAR PDG TON CAX TNF ATA PET)
declare -a SITEQ=( \
    $Q_M34 $Q_S67 $Q_HAR \
    $Q_PDG $Q_TON $Q_CAX \
    $Q_TNF $Q_ATA $Q_PET)

# SOI DEBUG 
declare -a D_IYEARAS=(1500 1500 2007 1500 2000 2000 2002 1500 2005)
declare -a D_IYEARZS=(1501 1501 2008 1501 2001 2001 2002 1501 2006)
declare -a D_INITMDS=(5    0    6    0    6    0    5    0    6   )
declare -a D_RUNTYPS=( \
    INITIAL INITIAL INITIAL \
    INITIAL INITIAL INITIAL \
    INITIAL INITIAL INITIAL)

# HI FREQUENCY DETAILED RUNS (HI-PET and HI-M34)
declare -a USE_HIFR=($USE_HIP $USE_HIM)
declare -a HIFRID=(hip him)
declare -a HIFRPFX=(HIP HIM)
declare -a HIFRQ=($Q_HIP $Q_HIM)
declare -a IDATEAH=(21 01)
declare -a IDATEZH=(28 08)
declare -a INITMDH=(6  5)
declare -a RUNTYPH=(INITIAL INITIAL)


let NPROC=NNODES*PPN

# GRID RUNS
declare -a USE_GRID=($USE_RJG)
declare -a GRIDID=(rjg)
declare -a GRIDPFX=(RJG)
declare -a GRIDQ=($Q_RJG)
declare -a GRIDPROC=($NPROC)
declare -a D_IYEARAG=(2008)
declare -a D_IYEARZG=(2008)
declare -a D_INITMDG=(5)
declare -a D_RUNTYPG=(INITIAL)


# SOME SPECIFICATION DEPEND ON RAPID OR LONG TYPE SIMULATIONS
if [ $TESTTYPE == "rapid" ]; then

    echo "PERFORMING 2 YEAR RAPID TESTS"

    declare -a IYEARAS=(1500 1500 2007 1500 2000 2000 2002 1500 2005)
    declare -a IYEARZS=(1502 1502 2009 1502 2002 2002 2004 1502 2007)
    declare -a INITMDS=(5    0    6    0    6    0    5    0    6   )
    declare -a RUNTYPS=( \
        INITIAL INITIAL INITIAL \
        INITIAL INITIAL INITIAL \
        INITIAL INITIAL INITIAL)

    # The gridded runs will only be 1 year
    declare -a IYEARAG=(2008)
    declare -a IYEARZG=(2008)
    declare -a INITMDG=(5)
    declare -a RUNTYPG=(INITIAL)

    UNITSTATE=2

elif [ $TESTTYPE == "long" ]; then

    echo "PERFORMING 300 YEAR LONG TESTS"

    declare -a IYEARAS=(1500 1500 1500 1500 1500 1500 1500 1500 1500)
    declare -a IYEARZS=(1800 1800 1800 1800 1800 1800 1800 1800 1800)
    declare -a INITMDS=(0    0    0    0    0    0    0    0    0   )
    declare -a RUNTYPS=( \
        INITIAL INITIAL INITIAL \
        INITIAL INITIAL INITIAL \
        INITIAL INITIAL INITIAL)

    declare -a IYEARAG=(1500)
    declare -a IYEARZG=(1535)
    declare -a INITMDG=(5)
    declare -a RUNTYPG=(INITIAL)

    UNITSTATE=3
else
    echo "PLEASE SPECIFY rapid OR long for TESTTYPE"
    exit
fi


for i in ${!SITEID[@]}
do
    if [ ${USE_SITE[i]} == "Y" ]; then
	echo "PROCESSING SITE: "${SITEID[i]}
    elif [ ${USE_SITE[i]} == "N" ]; then
	echo "SKIPPING SITE: "${SITEID[i]}
    else
	echo "IMPROPER USE SPECIFIER: "${SITEID[i]}
	echo "STOPPING"
	exit
    fi
done

echo ""
for i in ${!GRIDID[@]}
do
  if [ ${USE_GRID[i]} == "Y" ]; then
      echo "PROCESSING GRID: "${GRIDID[i]}
  elif [ ${USE_GRID[i]} == "N" ]; then
      echo "SKIPPING GRID: "${GRIDID[i]}
  else
      echo "IMPROPER USE SPECIFIER: "${GRIDID[i]}
      echo "STOPPING"
      exit
  fi
done

echo ""
for i in ${!HIFRID[@]}
do
    if [ ${USE_HIFR[i]} == "Y" ]; then
	echo "PROCESSING HI-FREQ: "${HIFRID[i]}
    elif [ ${USE_HIFR[i]} == "N" ]; then
	echo "SKIPPING HI-FREQ: "${HIFRID[i]}
    else
	echo "IMPROPER USE SPECIFIER: "${HIFRID[i]}
	echo "STOPPING"
	exit
    fi
done


#===============================================================================
#===============================================================================


# Update version string to contain the test type
VERSION=${VERSION}${TESTTYPE}


MAIN_EXE="./ed_2.1-main"
TEST_EXE="./ed_2.1-test"
DBUG_EXE="./ed_2.1-dbug"



# Create the working folder

echo ""
echo ""
echo "============================================="
echo "  GENERATE PRISTINE DIRECTORY: ${VERSION} ???"
echo ""
echo " ALL DATA WILL BE LOST!"
echo " ANSWER [Y/N]"
echo "============================================="


read flusher

if [ ${flusher} == "Y" ]; then
    echo "FLUSHING"
    sleep 2
    rm -rf ${VERSION}
    mkdir -p ${VERSION}

    # Generate the symbolic link to edts_datasets

    ln -s ${DATAPATH} ${VERSION}/edts_datasets

    ln -s $MAIN_EXE_PATH ${VERSION}/${MAIN_EXE}
    ln -s $TEST_EXE_PATH ${VERSION}/${TEST_EXE}
    ln -s $DBUG_EXE_PATH ${VERSION}/${DBUG_EXE}

    cp Templates/*xml ${VERSION}/

else
    if [ ${flusher} == "N" ]; then
	echo "WILL NOT FLUSH"
	sleep 2
    else
	echo "NEITHER Y OR N, EXITING"
	exit
    fi
fi


# Write out the text info to the report directory

echo '<?xml version="1.0" ?>' > ${VERSION}/test_text.xml
echo '<description> '  >> ${VERSION}/test_text.xml
echo '<branch_version> '$VERSION_BRANCHED_FROM' </branch_version>' >> ${VERSION}/test_text.xml
echo '<tester_name> '$TESTER_NAME' </tester_name>' >> ${VERSION}/test_text.xml
echo '<committer_name> '$COMMITTER_NAME' </committer_name>' >> ${VERSION}/test_text.xml
echo '<test_description> '$TEST_DESCRIPTION' </test_description>' >> ${VERSION}/test_text.xml
echo '</description>' >> ${VERSION}/test_text.xml


# Loop over the different SOI CASES
# ==============================================================================

ntasks=0
nnodes=1

for i in ${!SITEID[@]}
do
  if [ ${USE_SITE[i]} == "Y" ]; then
      let ntasks=ntasks+3
  fi
done
for i in ${!HIFRID[@]}
do

  if [ ${USE_HIFR[i]} == "Y" ]; then
      let ntasks=ntasks+3
fi
done

if [ ${ntasks} > 0 ]; then
    echo "TOTAL NUMBER OF SOI TASKS = "${ntasks}
    let nnodes=nnodes+`awk  'BEGIN { rounded = sprintf("%.0f", '${ntasks}/24' ); print rounded }'`
    echo "TOTAL NUMBER OF NODES = "${nnodes}
    let mppwidth=`awk  'BEGIN { rounded = sprintf("%.0f", '${nnodes}*24' ); print rounded }'`
    echo "MPPWIDTH = "${mppwidth}

else
    echo "NO SOI TASKS"
fi




echo '#!/bin/sh' > ${VERSION}/submit_batch.sh
chmod +x ${VERSION}/submit_batch.sh


tcount=0
for i in ${!SITEID[@]}
do

  if [ ${USE_SITE[i]} == "Y" ]; then

  TEMPLATEMAIN=Templates/ED2IN-${SITEPFX[i]}-MAIN
  TEMPLATETEST=Templates/ED2IN-${SITEPFX[i]}-TEST
  TEMPLATEDBUG=Templates/ED2IN-${SITEPFX[i]}-DBUG

  FILEMAIN=${VERSION}/ED2IN-${SITEPFX[i]}-MAIN
  FILETEST=${VERSION}/ED2IN-${SITEPFX[i]}-TEST
  FILEDBUG=${VERSION}/ED2IN-${SITEPFX[i]}-DBUG

  echo ""
  echo ""
  echo "Processing "${SITEID[i]}
  echo ""
  echo "Control Files: "$FILEMAIN", "$FILETEST", "$FILEDBUG
  echo "IO/Err Files: "${VERSION}/${SITEID[i]}"*"
  echo ""

  cp $TEMPLATEMAIN $FILEMAIN
  cp $TEMPLATETEST $FILETEST
  cp $TEMPLATEDBUG $FILEDBUG

# Modify the ED2IN RUNTYPE
# ========================

  sed -i '/NL%RUNTYPE/c\   NL%RUNTYPE = '\'${RUNTYPS[i]}\' $FILEMAIN
  sed -i '/NL%RUNTYPE/c\   NL%RUNTYPE = '\'${RUNTYPS[i]}\' $FILETEST
  sed -i '/NL%RUNTYPE/c\   NL%RUNTYPE = '\'${D_RUNTYPS[i]}\' $FILEDBUG

# Modify the ED2IN IED_INIT_MODE
# ==============================

  sed -i '/NL%IED_INIT_MODE/c\   NL%IED_INIT_MODE = '${INITMDS[i]} $FILEMAIN
  sed -i '/NL%IED_INIT_MODE/c\   NL%IED_INIT_MODE = '${INITMDS[i]} $FILETEST
  sed -i '/NL%IED_INIT_MODE/c\   NL%IED_INIT_MODE = '${D_INITMDS[i]} $FILEDBUG

# Modify the ED2IN start years
# ============================

  sed -i '/NL%IYEARA/c\   NL%IYEARA   = '${IYEARAS[i]} $FILEMAIN
  sed -i '/NL%IYEARA/c\   NL%IYEARA   = '${IYEARAS[i]} $FILETEST
  sed -i '/NL%IYEARA/c\   NL%IYEARA   = '${D_IYEARAS[i]} $FILEDBUG

# Modify the ED2IN end years
# ==========================

  sed -i '/NL%IYEARZ/c\   NL%IYEARZ   = '${IYEARZS[i]} $FILEMAIN
  sed -i '/NL%IYEARZ/c\   NL%IYEARZ   = '${IYEARZS[i]} $FILETEST
  sed -i '/NL%IYEARZ/c\   NL%IYEARZ   = '${D_IYEARZS[i]} $FILEDBUG

# Modify the ED2IN files to point to the desired output directories   
# =================================================================   
  
  sed -i '/NL%FFILOUT/c\   NL%FFILOUT = '\''F_main_'${SITEID[i]}'/main_'${SITEID[i]}''\' $FILEMAIN
  sed -i '/NL%SFILOUT/c\   NL%SFILOUT = '\''S_main_'${SITEID[i]}'/main_'${SITEID[i]}''\' $FILEMAIN
  sed -i '/NL%FFILOUT/c\   NL%FFILOUT = '\''F_test_'${SITEID[i]}'/test_'${SITEID[i]}''\' $FILETEST
  sed -i '/NL%SFILOUT/c\   NL%SFILOUT = '\''S_test_'${SITEID[i]}'/test_'${SITEID[i]}''\' $FILETEST
  sed -i '/NL%FFILOUT/c\   NL%FFILOUT = '\''F_dbug_'${SITEID[i]}'/dbug_'${SITEID[i]}''\' $FILEDBUG
  sed -i '/NL%SFILOUT/c\   NL%SFILOUT = '\''S_dbug_'${SITEID[i]}'/dbug_'${SITEID[i]}''\' $FILEDBUG

# Modify the state file output frequency
# ======================================

  sed -i '/NL%UNITSTATE/c\   NL%UNITSTATE  = '$UNITSTATE $FILEMAIN
  sed -i '/NL%UNITSTATE/c\   NL%UNITSTATE  = '$UNITSTATE $FILETEST
  sed -i '/NL%UNITSTATE/c\   NL%UNITSTATE  = '$UNITSTATE $FILEDBUG


# Reset and flush the output folders
# ====================================

mkdir -p ${VERSION}/F_main_${SITEID[i]}
mkdir -p ${VERSION}/F_test_${SITEID[i]}
mkdir -p ${VERSION}/F_dbug_${SITEID[i]}
mkdir -p ${VERSION}/S_main_${SITEID[i]}
mkdir -p ${VERSION}/S_test_${SITEID[i]}
mkdir -p ${VERSION}/S_dbug_${SITEID[i]}

# Execute
# =======

  let maincount=tcount+0
  let testcount=tcount+1
  let dbugcount=tcount+2
  
  # CREATE THE PBS FILE

  # FOR MAIN

  PBSFILE=${VERSION}/batch_main_${SITEID[i]}
  echo "#PBS -q serial"                                         > $PBSFILE
  echo "#PBS -l walltime=12:00:00"                             >> $PBSFILE
  echo "#PBS -l pvmem=3GB"                                     >> $PBSFILE
  echo "#PBS -N main_${SITEID[i]}"                             >> $PBSFILE
  echo '#PBS -e main_'${SITEID[i]}.'$PBS_JOBID.err' >> $PBSFILE
  echo '#PBS -o main_'${SITEID[i]}.'$PBS_JOBID.out' >> $PBSFILE
  echo "#PBS -V"                                               >> $PBSFILE
  echo ""                                                      >> $PBSFILE
  echo 'cd $PBS_O_WORKDIR'                                     >> $PBSFILE
  echo ""                                                      >> $PBSFILE
  echo "${MAIN_EXE} -f ED2IN-${SITEPFX[i]}-MAIN"                            >> $PBSFILE
  echo "" >> ${VERSION}/submit_batch.sh
  echo "qsub batch_main_${SITEID[i]}" >> ${VERSION}/submit_batch.sh
  
  # FOR TEST

  PBSFILE=${VERSION}/batch_test_${SITEID[i]}
  echo "#PBS -q serial"                                  > $PBSFILE
  echo "#PBS -l walltime=12:00:00"                      >> $PBSFILE
  echo "#PBS -l pvmem=3GB"                              >> $PBSFILE
  echo "#PBS -N test_${SITEID[i]}"                      >> $PBSFILE
  echo '#PBS -e test_'${SITEID[i]}.'$PBS_JOBID.err' >> $PBSFILE
  echo '#PBS -o test_'${SITEID[i]}.'$PBS_JOBID.out' >> $PBSFILE
  echo "#PBS -V"                                        >> $PBSFILE
  echo ""                                               >> $PBSFILE
  echo 'cd $PBS_O_WORKDIR'                              >> $PBSFILE
  echo ""                                               >> $PBSFILE
  echo "${TEST_EXE} -f ED2IN-${SITEPFX[i]}-TEST"                     >> $PBSFILE
  echo "" >> ${VERSION}/submit_batch.sh
  echo "qsub batch_test_${SITEID[i]}" >> ${VERSION}/submit_batch.sh

  # FOR DBUG

  PBSFILE=${VERSION}/batch_dbug_${SITEID[i]}
  echo "#PBS -q serial"                                  > $PBSFILE
  echo "#PBS -l walltime=12:00:00"                      >> $PBSFILE
  echo "#PBS -l pvmem=3GB"                              >> $PBSFILE
  echo "#PBS -N dbug_${SITEID[i]}"                      >> $PBSFILE
  echo '#PBS -e dbug_'${SITEID[i]}'.$PBS_JOBID.err' >> $PBSFILE
  echo '#PBS -o dbug_'${SITEID[i]}'.$PBS_JOBID.out' >> $PBSFILE
  echo "#PBS -V"                                        >> $PBSFILE
  echo ""                                               >> $PBSFILE
  echo 'cd $PBS_O_WORKDIR'                              >> $PBSFILE
  echo ""                                               >> $PBSFILE
  echo "${DBUG_EXE} -f ED2IN-${SITEPFX[i]}-DBUG"                     >> $PBSFILE
  echo "" >> ${VERSION}/submit_batch.sh
  echo "qsub batch_dbug_${SITEID[i]}" >> ${VERSION}/submit_batch.sh


  let tcount=tcount+3
  
fi

done


# Loop over the different HI-FREQUENCY CASES
# ==============================================================================

for i in ${!HIFRID[@]}
do

  if [ ${USE_HIFR[i]} == "Y" ]; then

  TEMPLATEMAIN=Templates/ED2IN-${HIFRPFX[i]}-MAIN
  TEMPLATETEST=Templates/ED2IN-${HIFRPFX[i]}-TEST
  TEMPLATEDBUG=Templates/ED2IN-${HIFRPFX[i]}-DBUG

  FILEMAIN=${VERSION}/ED2IN-${HIFRPFX[i]}-MAIN
  FILETEST=${VERSION}/ED2IN-${HIFRPFX[i]}-TEST
  FILEDBUG=${VERSION}/ED2IN-${HIFRPFX[i]}-DBUG

  cp $TEMPLATEMAIN $FILEMAIN
  cp $TEMPLATETEST $FILETEST
  cp $TEMPLATEDBUG $FILEDBUG

  echo ""
  echo ""
  echo "Processing "${HIFRID[i]}
  echo ""
  echo "Control Files: "$FILEMAIN", "$FILETEST
  echo ""

# Modify the ED2IN RUNTYPE
# ========================
  sed -i '/NL%RUNTYPE/c\   NL%RUNTYPE = '\'${RUNTYPH[i]}\' $FILEMAIN
  sed -i '/NL%RUNTYPE/c\   NL%RUNTYPE = '\'${RUNTYPH[i]}\' $FILETEST
  sed -i '/NL%RUNTYPE/c\   NL%RUNTYPE = '\'${RUNTYPH[i]}\' $FILEDBUG

# Modify the ED2IN IED_INIT_MODE
# ==============================

  sed -i '/NL%IED_INIT_MODE/c\   NL%IED_INIT_MODE = '${INITMDH[i]} $FILEMAIN
  sed -i '/NL%IED_INIT_MODE/c\   NL%IED_INIT_MODE = '${INITMDH[i]} $FILETEST
  sed -i '/NL%IED_INIT_MODE/c\   NL%IED_INIT_MODE = '${INITMDH[i]} $FILEDBUG


# Modify the ED2IN start dates
# ============================

  sed -i '/NL%IDATEA/c\   NL%IDATEA   = '${IDATEAH[i]} $FILEMAIN
  sed -i '/NL%IDATEA/c\   NL%IDATEA   = '${IDATEAH[i]} $FILETEST
  sed -i '/NL%IDATEA/c\   NL%IDATEA   = '${IDATEAH[i]} $FILEDBUG

# Modify the ED2IN end date
# ==========================

  sed -i '/NL%IDATEZ/c\   NL%IDATEZ   = '${IDATEZH[i]} $FILEMAIN
  sed -i '/NL%IDATEZ/c\   NL%IDATEZ   = '${IDATEZH[i]} $FILETEST
  sed -i '/NL%IDATEZ/c\   NL%IDATEZ   = '${IDATEZH[i]} $FILEDBUG

# Modify the ED2IN files to point to the desired output directories   
# =================================================================   
  
  sed -i '/NL%FFILOUT/c\   NL%FFILOUT = '\''F_main_'${HIFRID[i]}'/main_'${HIFRID[i]}''\' $FILEMAIN
  sed -i '/NL%SFILOUT/c\   NL%SFILOUT = '\''S_main_'${HIFRID[i]}'/main_'${HIFRID[i]}''\' $FILEMAIN
  sed -i '/NL%FFILOUT/c\   NL%FFILOUT = '\''F_test_'${HIFRID[i]}'/test_'${HIFRID[i]}''\' $FILETEST
  sed -i '/NL%SFILOUT/c\   NL%SFILOUT = '\''S_test_'${HIFRID[i]}'/test_'${HIFRID[i]}''\' $FILETEST
  sed -i '/NL%FFILOUT/c\   NL%FFILOUT = '\''F_dbug_'${HIFRID[i]}'/dbug_'${HIFRID[i]}''\' $FILEDBUG
  sed -i '/NL%SFILOUT/c\   NL%SFILOUT = '\''S_dbug_'${HIFRID[i]}'/dbug_'${HIFRID[i]}''\' $FILEDBUG


# Reset and flush the output folders
# ====================================

mkdir -p ${VERSION}/F_main_${HIFRID[i]}
mkdir -p ${VERSION}/F_test_${HIFRID[i]}
mkdir -p ${VERSION}/F_dbug_${HIFRID[i]}
mkdir -p ${VERSION}/S_main_${HIFRID[i]}
mkdir -p ${VERSION}/S_test_${HIFRID[i]}
mkdir -p ${VERSION}/S_dbug_${HIFRID[i]}


  let maincount=tcount+0
  let testcount=tcount+1
  let dbugcount=tcount+2

    # CREATE THE PBS FILE

  # FOR MAIN

  PBSFILE=${VERSION}/batch_main_${HIFRID[i]}
  echo "#PBS -q serial"                                         > $PBSFILE
  echo "#PBS -l walltime=12:00:00"                             >> $PBSFILE
  echo "#PBS -l pvmem=3GB"                                     >> $PBSFILE
  echo "#PBS -N main_${HIFRID[i]}"                             >> $PBSFILE
  echo '#PBS -e main_'${HIFRID[i]}.'$PBS_JOBID.err' >> $PBSFILE
  echo '#PBS -o main_'${HIFRID[i]}.'$PBS_JOBID.out' >> $PBSFILE
  echo "#PBS -V"                                               >> $PBSFILE
  echo ""                                                      >> $PBSFILE
  echo 'cd $PBS_O_WORKDIR'                                     >> $PBSFILE
  echo ""                                                      >> $PBSFILE
  echo "${MAIN_EXE} -f ED2IN-${HIFRPFX[i]}-MAIN"                            >> $PBSFILE
  echo "" >> ${VERSION}/submit_batch.sh
  echo "qsub batch_main_${HIFRID[i]}" >> ${VERSION}/submit_batch.sh

  # FOR TEST

  PBSFILE=${VERSION}/batch_test_${HIFRID[i]}
  echo "#PBS -q serial"                                  > $PBSFILE
  echo "#PBS -l walltime=12:00:00"                      >> $PBSFILE
  echo "#PBS -l pvmem=3GB"                              >> $PBSFILE
  echo "#PBS -N test_${HIFRID[i]}"                      >> $PBSFILE
  echo '#PBS -e test_'${HIFRID[i]}.'$PBS_JOBID.err' >> $PBSFILE
  echo '#PBS -o test_'${HIFRID[i]}.'$PBS_JOBID.out' >> $PBSFILE
  echo "#PBS -V"                                        >> $PBSFILE
  echo ""                                               >> $PBSFILE
  echo 'cd $PBS_O_WORKDIR'                              >> $PBSFILE
  echo ""                                               >> $PBSFILE
  echo "${TEST_EXE} -f ED2IN-${HIFRPFX[i]}-TEST"                     >> $PBSFILE
  echo "" >> ${VERSION}/submit_batch.sh
  echo "qsub batch_test_${HIFRID[i]}" >> ${VERSION}/submit_batch.sh

  # FOR DBUG

  PBSFILE=${VERSION}/batch_dbug_${HIFRID[i]}
  echo "#PBS -q serial"                                  > $PBSFILE
  echo "#PBS -l walltime=12:00:00"                      >> $PBSFILE
  echo "#PBS -l pvmem=3GB"                              >> $PBSFILE
  echo "#PBS -N dbug_${SITEID[i]}"                      >> $PBSFILE
  echo '#PBS -e dbug_'${HIFRID[i]}'.$PBS_JOBID.err' >> $PBSFILE
  echo '#PBS -o dbug_'${HIFRID[i]}'.$PBS_JOBID.out' >> $PBSFILE
  echo "#PBS -V"                                        >> $PBSFILE
  echo ""                                               >> $PBSFILE
  echo 'cd $PBS_O_WORKDIR'                              >> $PBSFILE
  echo ""                                               >> $PBSFILE
  echo "${DBUG_EXE} -f ED2IN-${HIFRPFX[i]}-DBUG"                     >> $PBSFILE
  echo "" >> ${VERSION}/submit_batch.sh
  echo "qsub batch_dbug_${HIFRID[i]}" >> ${VERSION}/submit_batch.sh

  let tcount=tcount+3


elif [ ${USE_HIFR[i]}=="N" ]; then
      echo "Skipping"
else
      echo "YOUR USE_HIFR VECTOR HAS AN ERROR"
      exit
fi

done

# Loop over the different GRIDDED CASES
# ==============================================================================

for i in ${!GRIDID[@]}
do

  if [ ${USE_GRID[i]} == "Y" ]; then

  TEMPLATEMAIN=Templates/ED2IN-${GRIDPFX[i]}-MAIN
  TEMPLATETEST=Templates/ED2IN-${GRIDPFX[i]}-TEST
  TEMPLATEDBUG=Templates/ED2IN-${GRIDPFX[i]}-DBUG

  FILEMAIN=${VERSION}/ED2IN-${GRIDPFX[i]}-MAIN
  FILETEST=${VERSION}/ED2IN-${GRIDPFX[i]}-TEST
  FILEDBUG=${VERSION}/ED2IN-${GRIDPFX[i]}-DBUG

  cp $TEMPLATEMAIN $FILEMAIN
  cp $TEMPLATETEST $FILETEST
  cp $TEMPLATEDBUG $FILEDBUG

  echo ""
  echo ""
  echo "Processing "${GRIDID[i]}
  echo ""
  echo "Control Files: "$FILEMAIN", "$FILETEST", "$FILEDBUG
  echo ""

# Modify the ED2IN RUNTYPE
# ========================

  sed -i '/NL%RUNTYPE/c\   NL%RUNTYPE = '\'${RUNTYPG[i]}\' $FILEMAIN
  sed -i '/NL%RUNTYPE/c\   NL%RUNTYPE = '\'${RUNTYPG[i]}\' $FILETEST
  sed -i '/NL%RUNTYPE/c\   NL%RUNTYPE = '\'${D_RUNTYPG[i]}\' $FILEDBUG

# Modify the ED2IN IED_INIT_MODE
# ==============================

  sed -i '/NL%IED_INIT_MODE/c\   NL%IED_INIT_MODE = '${INITMDG[i]} $FILEMAIN
  sed -i '/NL%IED_INIT_MODE/c\   NL%IED_INIT_MODE = '${INITMDG[i]} $FILETEST
  sed -i '/NL%IED_INIT_MODE/c\   NL%IED_INIT_MODE = '${D_INITMDG[i]} $FILEDBUG

# Modify the ED2IN start years
# ============================

  sed -i '/NL%IYEARA/c\   NL%IYEARA   = '${IYEARAG[i]} $FILEMAIN
  sed -i '/NL%IYEARA/c\   NL%IYEARA   = '${IYEARAG[i]} $FILETEST
  sed -i '/NL%IYEARA/c\   NL%IYEARA   = '${D_IYEARAG[i]} $FILEDBUG

# Modify the ED2IN end years
# ==========================

  sed -i '/NL%IYEARZ/c\   NL%IYEARZ   = '${IYEARZG[i]} $FILEMAIN
  sed -i '/NL%IYEARZ/c\   NL%IYEARZ   = '${IYEARZG[i]} $FILETEST
  sed -i '/NL%IYEARZ/c\   NL%IYEARZ   = '${D_IYEARZG[i]} $FILEDBUG

# Modify the ED2IN files to point to the desired output directories
# =================================================================

  sed -i '/NL%FFILOUT/c\   NL%FFILOUT = '\''F_main_'${GRIDID[i]}'/main_'${GRIDID[i]}''\' $FILEMAIN
  sed -i '/NL%SFILOUT/c\   NL%SFILOUT = '\''S_main_'${GRIDID[i]}'/main_'${GRIDID[i]}''\' $FILEMAIN
  sed -i '/NL%FFILOUT/c\   NL%FFILOUT = '\''F_test_'${GRIDID[i]}'/test_'${GRIDID[i]}''\' $FILETEST
  sed -i '/NL%SFILOUT/c\   NL%SFILOUT = '\''S_test_'${GRIDID[i]}'/test_'${GRIDID[i]}''\' $FILETEST
  sed -i '/NL%FFILOUT/c\   NL%FFILOUT = '\''F_dbug_'${GRIDID[i]}'/dbug_'${GRIDID[i]}''\' $FILEDBUG
  sed -i '/NL%SFILOUT/c\   NL%SFILOUT = '\''S_dbug_'${GRIDID[i]}'/dbug_'${GRIDID[i]}''\' $FILEDBUG

# Reset and flush the output folders
# ====================================

mkdir -p ${VERSION}/F_main_${GRIDID[i]}
mkdir -p ${VERSION}/F_test_${GRIDID[i]}
mkdir -p ${VERSION}/F_dbug_${GRIDID[i]}
mkdir -p ${VERSION}/S_main_${GRIDID[i]}
mkdir -p ${VERSION}/S_test_${GRIDID[i]}
mkdir -p ${VERSION}/S_dbug_${GRIDID[i]}

PBSFILE=${VERSION}/batch_main_${GRIDID[i]}

echo ""
echo "#PBS -q regular"                                > $PBSFILE
echo "#PBS -l nodes=${NNODES}:ppn=${PPN}"             >> $PBSFILE
echo "#PBS -l walltime=12:00:00"                     >> $PBSFILE
echo "#PBS -N main_${GRIDID[i]}"                      >> $PBSFILE
echo '#PBS -e main_'${GRIDID[i]}'.$PBS_JOBID.err' >> $PBSFILE
echo '#PBS -o main_'${GRIDID[i]}'.$PBS_JOBID.out' >> $PBSFILE
echo ""                                              >> $PBSFILE
echo 'cd $PBS_O_WORKDIR'                             >> $PBSFILE
echo ""
echo "mpirun -n ${NPROC} ${MAIN_EXE} -f ED2IN-${GRIDPFX[i]}-MAIN" >> $PBSFILE
echo "" >> ${VERSION}/submit_batch.sh
echo "qsub batch_main_${GRIDID[i]}" >> ${VERSION}/submit_batch.sh

PBSFILE=${VERSION}/batch_test_${GRIDID[i]}
echo ""
echo "#PBS -q regular"                                > $PBSFILE
echo "#PBS -l nodes=${NNODES}:ppn=${PPN}"             >> $PBSFILE
echo "#PBS -l walltime=12:00:00"                     >> $PBSFILE
echo "#PBS -N test_${GRIDID[i]}"                  >> $PBSFILE
echo '#PBS -e test_'${GRIDID[i]}'.$PBS_JOBID.err' >> $PBSFILE
echo '#PBS -o test_'${GRIDID[i]}'.$PBS_JOBID.out' >> $PBSFILE
echo ""                                              >> $PBSFILE
echo 'cd $PBS_O_WORKDIR'                             >> $PBSFILE
echo ""
echo "mpirun -n ${NPROC} ${TEST_EXE} -f ED2IN-${GRIDPFX[i]}-TEST" >> $PBSFILE
echo "" >> ${VERSION}/submit_batch.sh
echo "qsub batch_test_${GRIDID[i]}" >> ${VERSION}/submit_batch.sh

PBSFILE=${VERSION}/batch_dbug_${GRIDID[i]}
echo ""
echo "#PBS -q regular"                                > $PBSFILE
echo "#PBS -l nodes=${NNODES}:ppn=${PPN}"             >> $PBSFILE
echo "#PBS -l walltime=12:00:00"                     >> $PBSFILE
echo "#PBS -N dbug_${GRIDID[i]}"                     >> $PBSFILE
echo '#PBS -e dbug_'${GRIDID[i]}'.$PBS_JOBID.err' >> $PBSFILE
echo '#PBS -o dbug_'${GRIDID[i]}'.$PBS_JOBID.out' >> $PBSFILE
echo ""                                              >> $PBSFILE
echo 'cd $PBS_O_WORKDIR'                             >> $PBSFILE
echo ""
echo "mpirun -n ${NPROC} ${DBUG_EXE} -f ED2IN-${GRIDPFX[i]}-DBUG" >> $PBSFILE
echo "" >> ${VERSION}/submit_batch.sh
echo "qsub batch_dbug_${GRIDID[i]}" >> ${VERSION}/submit_batch.sh

fi

done
