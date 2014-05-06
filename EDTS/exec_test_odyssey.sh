#!/bin/bash


#===============================================================================
#                    User control variables
#===============================================================================

# This is the unique identifier tag for this test environment
# It should indicate the revision number branched from, the users initials
# and a test number.  For instance if you branched from r81, your initials
# are rk and you tried 4 iterations of this test before your commits
# were verified, you would use:

VERSION="83ml327v3"

# These are you edm executables.  They must exist.  Duh.

MAIN_EXE="./ed_2.1-main"
TEST_EXE="./ed_2.1-test"
DBUG_EXE="./ed_2.1-dbug"

# Decide on a "rapid" (2 year) or "long" (300 year)
# THE long tests will take a really really long time, beware

TESTTYPE="rapid"

# Decide on which sites to use.  You should ultimately
# have results for all of them, but this is helpfull if some
# sites give you problems after the first try
# Y for Yes and N for No

USE_M34=Y     #Manaus km 34 POI
USE_S67=Y     #Santarem km 67 POI
USE_HAR=Y     #Harvard forest SOI
USE_PDG=Y     #Pe de Gigante SOI 
USE_TON=Y     #Tonzi SOI
USE_CAX=Y     #Caxuana SOI
USE_TNF=Y     #Tapajos National Forest SOI
USE_ATA=Y     #Atacama Desert SOI
USE_PET=Y     #Petrolina SOI
USE_HIP=Y     #Petrolina High Frequency Detailed Short SOI
USE_HIM=Y     #Manaus High Frequency Detailed Short SOI
USE_RJG=Y     #Gridded 12x12 simulation centerd on Reserva Jaru

# How many cores do you want to use for the gridded simulations
# Currently there are 3 (RJG-MAIN RJG-TEST and RJG-DBUG)

NPROC=48

# HOW MANY MINUTES DO THE QUEUES NEED

SOI_FAST_TIME=720
SOI_LONG_TIME=14400

GRID_FAST_TIME=1080
GRID_LONG_TIME=14400


# WHICH QUEUES WOULD USE LIKE TO USE?
# MAIN-TEST-DBUG MUST ALL USE THE SAME QUEUE

Q_M34=general     #Manaus km 34 SOI
Q_S67=general     #Santarem km 67 SOI
Q_HAR=general     #Harvard forest SOI
Q_PDG=general     #Pe de Gigante SOI
Q_TON=general     #Tonzi SOI
Q_CAX=general     #Caxuana SOI
Q_TNF=general     #Tapajos National Forest SOI
Q_ATA=general     #Atacama Desert SOI
Q_PET=general     #Petrolina SOI
Q_HIP=general     #Petrolina High Frequency Detailed Short SOI
Q_HIM=general     #Manaus High Frequency Detailed Short SOI
Q_RJG=general     #Gridded 12x12 simulation centerd on Reserva Jaru

# Give an explanation of the tests.  Explain what the commits had involved.

TEST_DESCRIPTION="Committers Fixes: Fixed bug in fmean_ carbon values and some _py accounting variables.  revno: 327: 1 Minor bug fixes in the heat capacity of soils and vegetation. 2 Minor bug fix: mmean_broot wasn't being updated. 3 Updates of many R scripts. 4 Included many utility R scripts.   revno: 326  1. Various changes in the post-processing shell and R scripts. 2. Added POV-Ray representations of the plant community.  3. Added the option of turning off the site-, patch- and cohort-level means (which make the output considerably larger).  revno: 325:  1.  Updated R scripts. 2.  Added a new flag for MAXPATCH: if it is -1 or 1, it will force the fusion down to one patch for each land use type. revno: 324: Minor bug fixes on BRAMS, so it also works with the new variables.  revno: 323:  1. Added R scripts  2. Changed ed_state_vars.f90 so it only writes fmean, dmean, mmean, qmean variables to history if it makes sense to do so. revno: 322: 1.  Intermediate commit, to avoid losing stuff.  New output standards. 2. Minor bug fixes in fusion, particularly with mortality.  Fast means are now part of the patch/cohort dynamics because they are written to the output after the dynamics."

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

    SOI_TIME=${SOI_FAST_TIME}
    GRID_TIME=${GRID_FAST_TIME}

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

    SOI_TIME=${SOI_LONG_TIME}
    GRID_TIME=${GRID_LONG_TIME}

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

echo ""
echo ""
echo "CREATING FOLDERS:"
echo "S_"${VERSION}
echo "F_"${VERSION}
echo ${VERSION}"_report"
echo ""


# Create the output directories

mkdir -p S_${VERSION}
mkdir -p F_${VERSION}

REPORTDIR=${VERSION}_report
mkdir -p ${REPORTDIR}


# Write out the text info to the report directory

echo '<?xml version="1.0" ?>' > ${REPORTDIR}/test_text.xml
echo '<description> '  >> ${REPORTDIR}/test_text.xml
echo '<branch_version> '$VERSION_BRANCHED_FROM' </branch_version>' >> ${REPORTDIR}/test_text.xml
echo '<tester_name> '$TESTER_NAME' </tester_name>' >> ${REPORTDIR}/test_text.xml
echo '<committer_name> '$COMMITTER_NAME' </committer_name>' >> ${REPORTDIR}/test_text.xml
echo '<test_description> '$TEST_DESCRIPTION' </test_description>' >> ${REPORTDIR}/test_text.xml
echo '</description>' >> ${REPORTDIR}/test_text.xml


# Loop over the different SOI CASES
# ==============================================================================

for i in ${!SITEID[@]}
do

  if [ ${USE_SITE[i]} == "Y" ]; then

  FILEMAIN=ED2IN-${SITEPFX[i]}-MAIN
  FILETEST=ED2IN-${SITEPFX[i]}-TEST
  FILEDBUG=ED2IN-${SITEPFX[i]}-DBUG

  echo ""
  echo ""
  echo "Processing "${SITEID[i]}
  echo ""
  echo "Control Files: "$FILEMAIN", "$FILETEST", "$FILEDBUG
  echo "IO/Err Files: "${REPORTDIR}/${SITEID[i]}"*"
  echo ""

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
  
  sed -i '/NL%FFILOUT/c\   NL%FFILOUT = '\''F_'${VERSION}'/main_'${SITEID[i]}''\' $FILEMAIN
  sed -i '/NL%SFILOUT/c\   NL%SFILOUT = '\''S_'${VERSION}'/main_'${SITEID[i]}''\' $FILEMAIN
  sed -i '/NL%FFILOUT/c\   NL%FFILOUT = '\''F_'${VERSION}'/test_'${SITEID[i]}''\' $FILETEST
  sed -i '/NL%SFILOUT/c\   NL%SFILOUT = '\''S_'${VERSION}'/test_'${SITEID[i]}''\' $FILETEST
  sed -i '/NL%FFILOUT/c\   NL%FFILOUT = '\''F_'${VERSION}'/dbug_'${SITEID[i]}''\' $FILEDBUG
  sed -i '/NL%SFILOUT/c\   NL%SFILOUT = '\''S_'${VERSION}'/dbug_'${SITEID[i]}''\' $FILEDBUG

# Modify the state file output frequency

  sed -i '/NL%UNITSTATE/c\   NL%UNITSTATE  = '$UNITSTATE $FILEMAIN
  sed -i '/NL%UNITSTATE/c\   NL%UNITSTATE  = '$UNITSTATE $FILETEST
  sed -i '/NL%UNITSTATE/c\   NL%UNITSTATE  = '$UNITSTATE $FILEDBUG

# Remove all the relevant output files
# ====================================

rm -f F_${VERSION}/main_${SITEID[i]}*
rm -f F_${VERSION}/test_${SITEID[i]}*
rm -f F_${VERSION}/dbug_${SITEID[i]}*

rm -f S_${VERSION}/main_${SITEID[i]}*
rm -f S_${VERSION}/test_${SITEID[i]}*
rm -f S_${VERSION}/dbug_${SITEID[i]}*

# Remove std io and err files
# ===========================

  rm -f ${REPORTDIR}/${SITEID[i]}main_out
  rm -f ${REPORTDIR}/${SITEID[i]}main_err
  rm -f ${REPORTDIR}/${SITEID[i]}test_out
  rm -f ${REPORTDIR}/${SITEID[i]}test_err
  rm -f ${REPORTDIR}/${SITEID[i]}dbug_out
  rm -f ${REPORTDIR}/${SITEID[i]}dbug_err
  
# Remove output files
  
  rm -f F_${VERSION}/main_${SITEID[i]}*
  rm -f S_${VERSION}/main_${SITEID[i]}*
  rm -f F_${VERSION}/test_${SITEID[i]}*
  rm -f S_${VERSION}/test_${SITEID[i]}*
  rm -f F_${VERSION}/dbug_${SITEID[i]}*
  rm -f S_${VERSION}/dbug_${SITEID[i]}*

# Execute
# =======

sbatch -o ${REPORTDIR}/${SITEID[i]}main_out \
 -e ${REPORTDIR}/${SITEID[i]}main_err \
 -J ${VERSION}${SITEID[i]}main -t ${SOI_TIME} --mem-per-cpu=2048 \
 -p ${SITEQ[i]} -n 1 \
 --wrap="mpirun -np 1 ${MAIN_EXE} -f ${FILEMAIN}"
  
sbatch -o ${REPORTDIR}/${SITEID[i]}test_out \
 -e ${REPORTDIR}/${SITEID[i]}test_err \
 -J ${VERSION}${SITEID[i]}test -t ${SOI_TIME} --mem-per-cpu=2048 \
 -p ${SITEQ[i]} -n 1 \
 --wrap="mpirun -np 1 ${TEST_EXE} -f ${FILETEST}"
  
sbatch -o ${REPORTDIR}/${SITEID[i]}dbug_out \
 -e ${REPORTDIR}/${SITEID[i]}dbug_err \
 -J ${VERSION}${SITEID[i]}dbug -t ${SOI_TIME} --mem-per-cpu=2048 \
 -p ${SITEQ[i]} -n 1 \
 --wrap="mpirun -np 1 ${DBUG_EXE} -f ${FILEDBUG}"
  
elif [ ${USE_SITE[i]}=="N" ]; then
      BSUBSTR=""
else
      echo "YOUR USE_SITE VECTOR HAS AN ERROR"
      exit
fi

done


# Loop over the different GRIDDED CASES
# ==============================================================================

for i in ${!GRIDID[@]}
do

  if [ ${USE_GRID[i]} == "Y" ]; then

  FILEMAIN=ED2IN-${GRIDPFX[i]}-MAIN
  FILETEST=ED2IN-${GRIDPFX[i]}-TEST
  FILEDBUG=ED2IN-${GRIDPFX[i]}-DBUG

  echo ""
  echo ""
  echo "Processing "${GRIDID[i]}
  echo ""
  echo "Control Files: "$FILEMAIN", "$FILETEST", "$FILEDBUG
  echo "IO/Err Files: "${REPORTDIR}/${GRIDID[i]}"*"
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
  
  sed -i '/NL%FFILOUT/c\   NL%FFILOUT = '\''F_'${VERSION}'/main_'${GRIDID[i]}''\' $FILEMAIN
  sed -i '/NL%SFILOUT/c\   NL%SFILOUT = '\''S_'${VERSION}'/main_'${GRIDID[i]}''\' $FILEMAIN
  sed -i '/NL%FFILOUT/c\   NL%FFILOUT = '\''F_'${VERSION}'/test_'${GRIDID[i]}''\' $FILETEST
  sed -i '/NL%SFILOUT/c\   NL%SFILOUT = '\''S_'${VERSION}'/test_'${GRIDID[i]}''\' $FILETEST
  sed -i '/NL%FFILOUT/c\   NL%FFILOUT = '\''F_'${VERSION}'/dbug_'${GRIDID[i]}''\' $FILEDBUG
  sed -i '/NL%SFILOUT/c\   NL%SFILOUT = '\''S_'${VERSION}'/dbug_'${GRIDID[i]}''\' $FILEDBUG

# Remove std io and err files
# ===========================

  rm -f ${REPORTDIR}/${GRIDID[i]}main_out
  rm -f ${REPORTDIR}/${GRIDID[i]}main_err
  rm -f ${REPORTDIR}/${GRIDID[i]}test_out
  rm -f ${REPORTDIR}/${GRIDID[i]}test_err
  rm -f ${REPORTDIR}/${GRIDID[i]}dbug_out
  rm -f ${REPORTDIR}/${GRIDID[i]}dbug_err

# Remove output files

  rm -f F_${VERSION}/main_${GRIDID[i]}*
  rm -f S_${VERSION}/main_${GRIDID[i]}*
  rm -f F_${VERSION}/test_${GRIDID[i]}*
  rm -f S_${VERSION}/test_${GRIDID[i]}*
  rm -f F_${VERSION}/dbug_${GRIDID[i]}*
  rm -f S_${VERSION}/dbug_${GRIDID[i]}*

sbatch -o ${REPORTDIR}/${GRIDID[i]}main_out \
 -e ${REPORTDIR}/${GRIDID[i]}main_err \
 -J ${VERSION}${GRIDID[i]}main -t ${GRID_TIME} --mem-per-cpu=1024\
 -p ${GRIDQ[i]} -n ${GRIDPROC[i]} \
 --wrap="mpirun -np ${GRIDPROC[i]} ${MAIN_EXE} -f ${FILEMAIN}"
  
sbatch -o ${REPORTDIR}/${GRIDID[i]}test_out \
 -e ${REPORTDIR}/${GRIDID[i]}test_err \
 -J ${VERSION}${GRIDID[i]}test -t ${GRID_TIME} --mem-per-cpu=1024\
 -p ${GRIDQ[i]} -n ${GRIDPROC[i]} \
 --wrap="mpirun -np ${GRIDPROC[i]} ${TEST_EXE} -f ${FILETEST}"
  
sbatch -o ${REPORTDIR}/${GRIDID[i]}dbug_out \
 -e ${REPORTDIR}/${GRIDID[i]}dbug_err \
 -J ${VERSION}${GRIDID[i]}dbug -t ${GRID_TIME} --mem-per-cpu=1024\
 -p ${GRIDQ[i]} -n ${GRIDPROC[i]} \
 --wrap="mpirun -np ${GRIDPROC[i]} ${DBUG_EXE} -f ${FILEDBUG}"

  
elif [ ${USE_GRID[i]}=="N" ]; then
      BSUBSTR=""
else
      echo "YOUR USE_GRID VECTOR HAS AN ERROR"
      exit
fi

done


# Loop over the different HI-FREQUENCY CASES
# ==============================================================================

for i in ${!HIFRID[@]}
do

  if [ ${USE_HIFR[i]} == "Y" ]; then

  FILEMAIN=ED2IN-${HIFRPFX[i]}-MAIN
  FILETEST=ED2IN-${HIFRPFX[i]}-TEST
  FILEDBUG=ED2IN-${HIFRPFX[i]}-DBUG

  echo ""
  echo ""
  echo "Processing "${HIFRID[i]}
  echo ""
  echo "Control Files: "$FILEMAIN", "$FILETEST
  echo "IO/Err Files: "${REPORTDIR}/${HIFRID[i]}"*"
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
  
  sed -i '/NL%FFILOUT/c\   NL%FFILOUT = '\''F_'${VERSION}'/main_'${HIFRID[i]}''\' $FILEMAIN
  sed -i '/NL%SFILOUT/c\   NL%SFILOUT = '\''S_'${VERSION}'/main_'${HIFRID[i]}''\' $FILEMAIN
  sed -i '/NL%FFILOUT/c\   NL%FFILOUT = '\''F_'${VERSION}'/test_'${HIFRID[i]}''\' $FILETEST
  sed -i '/NL%SFILOUT/c\   NL%SFILOUT = '\''S_'${VERSION}'/test_'${HIFRID[i]}''\' $FILETEST
  sed -i '/NL%FFILOUT/c\   NL%FFILOUT = '\''F_'${VERSION}'/dbug_'${HIFRID[i]}''\' $FILEDBUG
  sed -i '/NL%SFILOUT/c\   NL%SFILOUT = '\''S_'${VERSION}'/dbug_'${HIFRID[i]}''\' $FILEDBUG


# Remove std io and err files
# ===========================

  rm -f ${REPORTDIR}/${HIFRID[i]}main_out
  rm -f ${REPORTDIR}/${HIFRID[i]}main_err
  rm -f ${REPORTDIR}/${HIFRID[i]}test_out
  rm -f ${REPORTDIR}/${HIFRID[i]}test_err
  rm -f ${REPORTDIR}/${HIFRID[i]}dbug_out
  rm -f ${REPORTDIR}/${HIFRID[i]}dbug_err

# Remove output files
  
  rm -f F_${VERSION}/main_${HIFRID[i]}*
  rm -f S_${VERSION}/main_${HIFRID[i]}*
  rm -f F_${VERSION}/test_${HIFRID[i]}*
  rm -f S_${VERSION}/test_${HIFRID[i]}*
  rm -f F_${VERSION}/dbug_${HIFRID[i]}*
  rm -f S_${VERSION}/dbug_${HIFRID[i]}*

sbatch -o ${REPORTDIR}/${HIFRID[i]}main_out \
 -e ${REPORTDIR}/${HIFRID[i]}main_err \
 -J ${VERSION}${HIFRID[i]}main -t ${SOI_TIME} --mem-per-cpu=2048 \
 -p ${HIFRQ[i]} -n 1 \
 --wrap="mpirun -np 1 ${MAIN_EXE} -f ${FILEMAIN}"
  
sbatch -o ${REPORTDIR}/${HIFRID[i]}test_out \
 -e ${REPORTDIR}/${HIFRID[i]}test_err \
 -J ${VERSION}${HIFRID[i]}test -t ${SOI_TIME} --mem-per-cpu=2048 \
 -p ${HIFRQ[i]} -n 1 \
 --wrap="mpirun -np 1 ${TEST_EXE} -f ${FILETEST}"

sbatch -o ${REPORTDIR}/${HIFRID[i]}dbug_out \
 -e ${REPORTDIR}/${HIFRID[i]}dbug_err \
 -J ${VERSION}${HIFRID[i]}dbug -t ${SOI_TIME} --mem-per-cpu=2048 \
 -p ${HIFRQ[i]} -n 1 \
 --wrap="mpirun -np 1 ${DBUG_EXE} -f ${FILEDBUG}"
  
elif [ ${USE_HIFR[i]}=="N" ]; then
      BSUBSTR=""
else
      echo "YOUR USE_HIFR VECTOR HAS AN ERROR"
      exit
fi

done

echo ""
echo ""
echo "===================================================================="
echo "EVERYONE LOVES SLURRRRRRRRRRM"
squeue -u ${USER}
echo "===================================================================="
