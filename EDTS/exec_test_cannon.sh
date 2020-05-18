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
MAIN_EXE_PATH="${EXE_PATH}/ed_2.2-main"
TEST_EXE_PATH="${EXE_PATH}/ed_2.2-test"
DBUG_EXE_PATH="${EXE_PATH}/ed_2.2-dbug"
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
#  RC_OPTS -- Optional argument to pass to .bashrc.  This only matters if you based your   #
#             .bashrc on MLO's file.  It allows the job to load different modules on       #
#             odyssey.  Otherwise, you can leave it blank, then the jobs will load the     #
#             default configuration.                                                       #
#------------------------------------------------------------------------------------------#
RC_OPTS="-n"
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
USE_HBG="y"     # Harvard (Bare ground) POI
USE_HIP="y"     # Petrolina High Frequency Detailed Short POI
USE_HIM="y"     # Manaus High Frequency Detailed Short POI
USE_HIG="y"     # Paracou High Frequency Detailed Short POI
USE_HIH="y"     # Harvard High Frequency Detailed Short POI
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
Q_HBG="moorcroft_amd"      # Harvard (bare ground) POI
Q_HIP="moorcroft_6100"     # Petrolina High Frequency Detailed Short POI
Q_HIM="moorcroft_6100"     # Manaus High Frequency Detailed Short POI
Q_HIG="moorcroft_6100"     # Paracou High Frequency Detailed Short POI
Q_HIH="moorcroft_6100"     # Harvard High Frequency Detailed Short POI
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
MEM_HBG="1845"     # Harvard (Bare ground) POI
MEM_HIP="1845"     # Petrolina High Frequency Detailed Short POI
MEM_HIM="1845"     # Manaus High Frequency Detailed Short POI
MEM_HIG="1845"     # Paracou High Frequency Detailed Short POI
MEM_HIH="1845"     # Harvard High Frequency Detailed Short POI
MEM_RJG="2000"     # Gridded 12x12 simulation centred on Reserva Jaru
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
CPU_HBG="1"    # Harvard (Bare Ground) POI (MAXPATCH=20; RAPID_INIT_MODE=-1)
CPU_RJG="16"   # Gridded 12x12 simulation centred on Reserva Jaru (NPTS=144)
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#     The following variables are the time required by each simulation, hh:mm:ss.          #
#------------------------------------------------------------------------------------------#
POI_TIME="Infinite"
GRID_TIME="Infinite"
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#     Number of nodes for the gridded simulations.  Fewer nodes means that the CPUs will   #
# be more clustered.  In case you don't mind how many nodes will be used, set it to 0.     #
#------------------------------------------------------------------------------------------#
GRID_NODE="1"
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



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#    List of functions to handle specific namelist requirements.                           #
#------------------------------------------------------------------------------------------#
   #---------------------------------------------------------------------------------------#
   #     This function extracts the namelist settings from the list.  This allows using    #
   # different settings for the mainline and the test simulation (e.g., when a namelist    #
   # setting doesn't exist in one or the other version).                                   #
   #---------------------------------------------------------------------------------------#
   function nml_setting(){
      usc=$(echo "${1}" | grep -i "__" | wc -l | awk '{print $1}')
      case ${usc} in
      0)
         nl_main="${1}"
         nl_test="${1}"
         ;;
      *) 
         nl_main=$(echo ${1} | sed s@"__"@" "@g | awk '{print $1}')
         nl_test=$(echo ${1} | sed s@"__"@" "@g | awk '{print $2}')
         ;;
      esac
   }
   #---------------------------------------------------------------------------------------#





   #---------------------------------------------------------------------------------------#
   #      This function sets soil variables based on soil depth flag.                      #
   #---------------------------------------------------------------------------------------#
   function soil_setting(){
      case "${1}" in
      A)
         nzg=7
         slz="-0.623,-0.501,-0.388,-0.283,-0.188,-0.106,-0.040"
         slm="  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00"
         slt="  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00"
         ;;
      B)
         nzg=10
         slz="-2.021,-1.689,-1.382,-1.101,-0.846,-0.620,-0.424,-0.260,-0.130,-0.040"
         slm="  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00"
         slt="  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00"
         ;;
      C)
         nzg=9
         slz="-3.305,-2.609,-1.995,-1.464,-1.015,-0.648,-0.364,-0.161,-0.040"
         slm="  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00"
         slt="  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00"
         ;;
      D)
         nzg=10
         slz="-4.084,-3.305,-2.609,-1.995,-1.464,-1.015,-0.648,-0.364,-0.161,-0.040"
         slm="  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00"
         slt="  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00"
         ;;
      E)
         nzg=11
         slz="-4.946,-4.084,-3.305,-2.609,-1.995,-1.464"
         slz="${slz},-1.015,-0.648,-0.364,-0.161,-0.040"
         slm="  1.00,  1.00,  1.00,  1.00,  1.00,  1.00"
         slm="${slm},  1.00,  1.00,  1.00,  1.00,  1.00"
         slt="  0.00,  0.00,  0.00,  0.00,  0.00,  0.00"
         slt="${slt},  0.00,  0.00,  0.00,  0.00,  0.00"
         ;;
      F)
         nzg=12
         slz="-5.891,-4.946,-4.084,-3.305,-2.609,-1.995,-1.464"
         slz="${slz},-1.015,-0.648,-0.364,-0.161,-0.040"
         slm="  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00"
         slm="${slm},  1.00,  1.00,  1.00,  1.00,  1.00"
         slt="  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00"
         slt="${slt},  0.00,  0.00,  0.00,  0.00,  0.00"
         ;;
      G)
         nzg=13
         slz="-6.919,-5.891,-4.946,-4.084,-3.305,-2.609,-1.995"
         slz="${slz},-1.464,-1.015,-0.648,-0.364,-0.161,-0.040"
         slm="  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00"
         slm="${slm},  1.00,  1.00,  1.00,  1.00,  1.00,  1.00"
         slt="  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00"
         slt="${slt},  0.00,  0.00,  0.00,  0.00,  0.00,  0.00"
         ;;
      H)
         nzg=14
         slz="-8.029,-6.919,-5.891,-4.946,-4.084,-3.305,-2.609,-1.995"
         slz="${slz},-1.464,-1.015,-0.648,-0.364,-0.161,-0.040"
         slm="  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00"
         slm="${slm},  1.00,  1.00,  1.00,  1.00,  1.00,  1.00"
         slt="  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00"
         slt="${slt},  0.00,  0.00,  0.00,  0.00,  0.00,  0.00"
         ;;
      I)
         nzg=16
         slz="-10.50,-9.223,-8.029,-6.919,-5.891,-4.946,-4.084,-3.305,-2.609"
         slz="${slz},-1.995,-1.464,-1.015,-0.648,-0.364,-0.161,-0.040"
         slm="  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00"
         slm="${slm},  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00"
         slt="  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00"
         slt="${slt},  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00"
         ;;
      *)
         nzg=12
         slz="-5.891,-4.946,-4.084,-3.305,-2.609,-1.995,-1.464"
         slz="${slz},-1.015,-0.648,-0.364,-0.161,-0.040"
         slm="  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00"
         slm="${slm},  1.00,  1.00,  1.00,  1.00,  1.00"
         slt="  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00"
         slt="${slt},  0.00,  0.00,  0.00,  0.00,  0.00"
         ;;
      esac
   }
   #---------------------------------------------------------------------------------------#





   #---------------------------------------------------------------------------------------#
   #      This function sets PFT variables based on the PFT flag.                          #
   #---------------------------------------------------------------------------------------#
   function pft_setting(){
      case "${1}" in
      A)
         inc_pft="1,2,3,4"
         pasture="1"
         agri="1"
         plantation="3"
         sl_pft="2,3,4"
         sl_probharv="1.,1.,1."
         sl_mindbh="50.,50.,50."
         ;;
      B)
         inc_pft="1,2,3,4,16"
         pasture="1"
         agri="16"
         plantation="3"
         sl_pft="2,3,4"
         sl_probharv="1.,1.,1."
         sl_mindbh="50.,50.,50."
         ;;
      C)
         inc_pft="5,6,8,9,10,11"
         pasture="5"
         agri="5"
         plantation="6"
         sl_pft="6,8,9,10,11"
         sl_probharv="1.,1.,1.,1.,1."
         sl_mindbh="0.,0.,0.,0.,0."
         ;;
      D)
         inc_pft="5,8,11"
         pasture="5"
         agri="5"
         plantation="8"
         sl_pft="8,11"
         sl_probharv="1.,1."
         sl_mindbh="0.,0."
         ;;
      E)
         inc_pft="5,6,7,8,9,10,11"
         pasture="5"
         agri="5"
         plantation="7"
         sl_pft="6,7,8,9,10,11"
         sl_probharv="1.,1.,1.,1.,1.,1."
         sl_mindbh="0.,0.,0.,0.,0.,0."
         ;;
      esac
   }
   #---------------------------------------------------------------------------------------#







   #---------------------------------------------------------------------------------------#
   #      This function sets the initial conditions file depending on the initialisation   #
   # type.  Sometimes this is just a dummy variable because the model does not use it      #
   # (e.g., when running near-bare ground or bare ground initial conditions.               #
   #---------------------------------------------------------------------------------------#
   function sfilin_setting()
   {
      case "${2}" in
      -1|0)
         #----- Dummy, not needed. --------------------------------------------------------#
         init_file="edts_datasets/inits/${1}/${1}."
         #---------------------------------------------------------------------------------#
         ;;
      5)
         #----- You may need to edit this block if you add new sites. ---------------------#
         case "${1}" in
         ton)     init_file="edts_datasets/inits/kz_tonzi/test_ton-S-1912" ;;
         s83)     init_file="edts_datasets/inits/ts83_hist/ts83_log"       ;;
         gyf|hig) init_file="edts_datasets/inits/gyf/tgyf"                 ;;
         *)       init_file="edts_datasets/inits/sci_006_potveg_links/PVE" ;;
         esac
         #---------------------------------------------------------------------------------#
         ;;
      6)
         #----- You may need to edit this block if you add new sites. ---------------------#
         case "${1}" in
         har|hih) init_file="edts_datasets/inits/har/ems." ;;
         pet|hip) init_file="edts_datasets/inits/pet/pnz_default." ;;
         prg)     init_file="edts_datasets/inits/prg_als/prg_als." ;;
         *)       echo " - Site ${1} not configured for IED_INIT_MODE=6"; exit 1 ;;
         esac
         #---------------------------------------------------------------------------------#
         ;;
      esac
   }
   #---------------------------------------------------------------------------------------#







   #---------------------------------------------------------------------------------------#
   #      This function sets the meteorological drivers given the site.                    #
   #---------------------------------------------------------------------------------------#
   function metdriv_setting()
   {
      case "${1}" in
      cax)
         metdriv="edts_datasets/met/caxtfe/CAXTFE_HEADER"
         metcyca=2000
         metcycz=2008
         imetavg=1
         imetrad=0
         ;;
      gyf|hig)
         metdriv="edts_datasets/met/wfdei/frguiana/WFDEI_FRGUIANA_NOLAND_GPCC_HEADER"
         metcyca=1979
         metcycz=2013
         imetavg=1
         imetrad=2
         ;;
      har|hih|hbg)
         metdriv="edts_datasets/met/har/HARVARD_MET_93_09_CO2"
         metcyca=1993
         metcycz=2009
         imetavg=2
         imetrad=0
         ;;
      pet|hip)
         metdriv="edts_datasets/met/pet/Petrolina_HEADER"
         metcyca=2005
         metcycz=2011
         imetavg=1
         imetrad=2
         ;;
      prg)
         metdriv="edts_datasets/met/wfdei/paragominas/WFDEI_PARAGOMINAS_NOLAND_GPCC_HEADER"
         metcyca=1979
         metcycz=2013
         imetavg=1
         imetrad=2
         ;;
      s67|tnf)
         metdriv="edts_datasets/met/s67_nlevine/Santarem_KM67_HEADER"
         metcyca=2002
         metcycz=2004
         imetavg=2
         imetrad=0
         ;;
      s83)
         metdriv="edts_datasets/met/wfdei/tapajos/WFDEI_TAPAJOS_NOLAND_GPCC_HEADER"
         metcyca=1979
         metcycz=2013
         imetavg=1
         imetrad=2
         ;;
      tl2)
         metdriv="edts_datasets/met/toolik2/ED_MET_DRIVER_HEADER-toolik.obs.rad"
         metcyca=1997
         metcycz=2010
         imetavg=3
         imetrad=0
         ;;
      ton)
         metdriv="edts_datasets/met/tonzi_kz/tonzi_driver_co2"
         metcyca=2000
         metcycz=2012
         imetavg=0
         imetrad=0
         ;;
      *)
         metdriv="edts_datasets/met/mod_ds314/SHEF_NCEP_DRIVER_DS314"
         metcyca=1969
         metcycz=2008
         imetavg=2
         imetrad=2
         ;;
      esac
   }
   #---------------------------------------------------------------------------------------#







   #---------------------------------------------------------------------------------------#
   #      This function sets the land use files given the site.  Currently this is a dummy #
   # function, but could adjust land use data sets in the future.                          #
   #---------------------------------------------------------------------------------------#
   function lu_setting()
   {
      case "${1}" in
         *) lu_database="edts_datasets/land_use/glu+sa2/blend-" ;;
      esac
   }
   #---------------------------------------------------------------------------------------#







   #---------------------------------------------------------------------------------------#
   #      This function sets the phenology files given the site.                           #
   #---------------------------------------------------------------------------------------#
   function phenol_setting()
   {
      case "${1}" in
         har|hih|hbg)
            phen_path="edts_datasets/phenology/HarOneCycle/phenology"
            iphenysa=1992
            iphenysz=2009
            iphenyfa=1992
            iphenyfz=2009
            ;;
         *) 
            phen_path="edts_datasets/phenology/${1}/phenology"
            iphenysa=1992
            iphenysz=2003
            iphenyfa=1992
            iphenyfz=2003
            ;;
      esac
   }
   #---------------------------------------------------------------------------------------#







   #---------------------------------------------------------------------------------------#
   #      This function sets the XML files given the site.                                 #
   #---------------------------------------------------------------------------------------#
   function xml_setting()
   {
      case "${1}" in
         tl2) edcnfgf="tl2_test.xml" ;;
         *)   edcnfgf="config.xml"   ;;
      esac
   }
   #---------------------------------------------------------------------------------------#







   #---------------------------------------------------------------------------------------#
   #      Temporary function to set water conductance depending on the water limitation    #
   # scheme and latitude.                                                                  #
   #---------------------------------------------------------------------------------------#
   function kwd0_setting()
   {
      case "${1}" in
      har|ton|tl2|hbg|hih)
         case "${2}" in
            5) kw_grass=12.  ; kw_tree=12.  ;;
            *) kw_grass=300. ; kw_tree=300. ;;
         esac
         d0_grass=0.010
         d0_tree=0.010
         ;;
      *)
         case "${2}" in
            5) kw_grass=25.  ; kw_tree=20.  ;;
            *) kw_grass=900. ; kw_tree=600. ;;
         esac
         d0_grass=0.016
         d0_tree=0.016
         ;;
      esac
   }
   #---------------------------------------------------------------------------------------#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#



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
                      ${USE_PRG} ${USE_HBG})
declare -a SITEID=(m34 s67 har pdg ton cax tnf ata pet gyf s83 prg tl2 hbg)
declare -a SITEPFX=(M34 S67 HAR PDG TON CAX TNF ATA PET GYF S83 PRG TL2 HBG)
declare -a SITEQ=( ${Q_M34} ${Q_S67} ${Q_HAR} ${Q_PDG} ${Q_TON} ${Q_CAX} \
                   ${Q_TNF} ${Q_ATA} ${Q_PET} ${Q_GYF} ${Q_S83} ${Q_PRG} \
                   ${Q_TL2} ${Q_HBG})
declare -a SITEMEM=( ${MEM_M34} ${MEM_S67} ${MEM_HAR} ${MEM_PDG} ${MEM_TON} ${MEM_CAX} \
                     ${MEM_TNF} ${MEM_ATA} ${MEM_PET} ${MEM_GYF} ${MEM_S83} ${MEM_PRG} \
                     ${MEM_TL2} ${MEM_HBG})
declare -a SITECPU=( ${CPU_M34} ${CPU_S67} ${CPU_HAR} ${CPU_PDG} ${CPU_TON} ${CPU_CAX} \
                     ${CPU_TNF} ${CPU_ATA} ${CPU_PET} ${CPU_GYF} ${CPU_S83} ${CPU_PRG} \
                     ${CPU_TL2} ${CPU_HBG} )
#---- Initial/final month and day (years are defined below, depending on the suite. -------#
declare -a IMONTHAS=(01 01 05 01 03 01 01 01 01 01 01 01 05 01)
declare -a IMONTHZS=(01 01 05 01 03 01 01 01 01 01 01 01 05 01)
declare -a IDATEAS=(01 01 01 01 01 01 01 01 01 01 01 01 01 01)
declare -a IDATEZS=(01 01 01 01 01 01 01 01 01 01 01 01 01 01)
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     POI debug year and settings. We ought to run at least two years to ensure at least   #
# one patch dynamics step is tested.                                                       #
#------------------------------------------------------------------------------------------#
declare -a D_IYEARAS=(2010 2010 2010 2010 2010 2010 2010 2010 2010 2010 2010 2010 2010 2010)
declare -a D_IYEARZS=(2012 2012 2012 2012 2012 2012 2012 2012 2012 2012 2012 2012 2012 2012)
declare -a D_INITMDS=(5    0    6    0    5    0    5    0    6    0    5    6    0    -1  )
declare -a D_RUNTYPS=( INITIAL INITIAL INITIAL INITIAL INITIAL INITIAL INITIAL \
                       INITIAL INITIAL INITIAL INITIAL INITIAL INITIAL INITIAL )
#------------------------------------------------------------------------------------------#




#----- High-frequency detailed runs. ------------------------------------------------------#
declare -a USE_HIFR=(${USE_HIP} ${USE_HIM} ${USE_HIG} ${USE_HIH})
declare -a HIFRID=(hip him hig hih)
declare -a HIFRPFX=(HIP HIM HIG HIH)
declare -a HIFRQ=(${Q_HIP} ${Q_HIM} ${Q_HIG} ${Q_HIH})
declare -a HIFRMEM=(${MEM_HIP} ${MEM_HIM} ${MEM_HIG} ${MEM_HIH})
declare -a HIFRCPU=(1 1 1 1)
declare -a INITMDH=(6  5  5   6)
declare -a RUNTYPH=(INITIAL INITIAL INITIAL INITIAL)
#------------------------------------------------------------------------------------------#


#----- Gridded runs. ----------------------------------------------------------------------#
declare -a USE_GRID=(${USE_RJG})
declare -a GRIDID=(rjg)
declare -a GRIDPFX=(RJG)
declare -a GRIDQ=(${Q_RJG})
declare -a GRIDMEM=(${MEM_RJG})
declare -a GRIDCPU=(${CPU_RJG})
declare -a IMONTHAG=(01)
declare -a IMONTHZG=(01)
declare -a IDATEAG=(01)
declare -a IDATEZG=(01)
declare -a D_IYEARAG=(2010)
declare -a D_IYEARZG=(2011)
declare -a D_INITMDG=(5)
declare -a D_RUNTYPG=(INITIAL)
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#      Variable settings (site runs).   In case the test and mainline must be set          #
# differently (e.g., new feature in ED2IN), use | to separate the two entries.  The order  #
# should be always main__test and debugging option is always set exactly the same as test. #
#------------------------------------------------------------------------------------------#
declare -a PLONS=(-60.209 -54.959 -72.170 -47.650 -130.97 -51.460 -54.959 -67.478 -40.370 -52.912 -54.971 -46.837 -149.60 -72.170)
declare -a PLATS=( -2.609  -2.857  42.540 -21.619  38.432  -1.720  -2.857 -20.509  -9.165   5.282  -3.018  -2.552  68.630  42.540)
declare -a ISOILFLGS=(1 2 2 2 2 2 2 1 2 2 2 2 2 2)
declare -a NSLCONS=(6 6 6 6 6 6 6 6 6 6 6 6 6 6)
declare -a SLXCLAYS=( -1.00 0.590 0.060 0.460 0.448 0.150 0.590 -1.00 0.052 0.345 0.590 -1.00 0.060 0.060)
declare -a SLXSANDS=( -1.00 0.390 0.660 0.410 0.255 0.780 0.390 -1.00 0.821 0.562 0.390 -1.00 0.100 0.660)
declare -a ISOILCOLS=(2 2 2 21 21 2 2 21 14 14 14 14 2 9)
declare -a SDEPTHS=(H I E E B H I B C G I H D A)
declare -a ISOILBCS=(1 1 0 1 0 1 1 2 2 1 1 1 2 0)
declare -a INTEGRS=(3 1 1 1 1 1 1 1 1 1 1 1 1 1)
declare -a IBRANCHS=(1 0 0 0 0 0 0 0 1 1 1 1 2 0)
declare -a IPHYSIOLS=(2 2 0 2 2 2 2 2 2 2__3 2__3 2__3 0__1 2)
declare -a IALLOMS=(2 2 0 2 0 2 2 2 2 3 3 3 0 0)
declare -a IECONS=(0 0 0 0 0 0 0 0 0 1 1 1 0 0)
declare -a IGRASSS=(1 0 0 0 0 1 0 0 0 1 1 1 1 0)
declare -a IREPROS=(2 2 2 2 2 2 2 2 2 2__3 2__3 2__3 2 0)
declare -a DECOMPS=(0 0 0 0 0 0 0 0 2 2__5 2__5 2__5 0 2)
declare -a H2OLIMS=(2 1 1 1 1 1 1 1 2 2__5 2 2__5 1 2)
declare -a STGROWS=(0 0 0 0 0 0 0 0 0 1 1 1 0 0)
declare -a PLASTICS=(0 0 0 0 0 0 0 0 0 2 2 2 0 0)
declare -a IANTHS=(1 0 0 0 0 0 0 0 0 0 1__2 0 0 0)
declare -a IFIRES=(2 0 0 2 0 0 0 2 0 0 0 2 0 0)
declare -a SMFIRES=(-1.40 -1.40 0.068 0.068 -1.40 -1.40 0.068 -1.40 -1.40 -1.40 -1.40 -1.40 0.068 -1.40)
declare -a PHENOLS=(2 2 1 2 2 2 2 2 2 -1 -1 -1 -1 1)
declare -a THCRITS=(-1.20 -1.20 0.09 -1.20 -1.20 0.09 0.09 -1.20 -1.20 -1.20 -1.20 -1.20 0.09 0.09)
declare -a ICANTURBS=(3 2 2 2 2 2 2 3 2 0 0 0 2 2)
declare -a ISFCLYRMS=(3 3 3 3 3 3 3 3 3 4 4 4 3 4)
declare -a IPERCOLS=(0 1 1 1 1 0 0 1 0 0 0 0 1 0)
declare -a PFTFLAGS=(B B C B E A A A B B B B D C)
declare -a TREEFALLS=(0.014 0.014 0.0 0.014 0.0 0.014 0.014 0.014 0.014 0.014 0.014 0.014 0.0 0.0)
declare -a IFUSIONS=(0 0 0 0 0 0 0 0 0 1 1 1 0 0)
declare -a IVEGTDYNS=(1 1 1 1 1 1 1 1 1 1 1 1 1 1)
declare -a RK4TOLERS=(0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01)
declare -a MAXSITES=(3 1 1 1 1 1 1 1 1 1 1 1 1 1)
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#      Variable settings (high-frequency runs).   In case the test and mainline must be    #
# set differently (e.g., new feature in ED2IN), use | to separate the two entries.  The    #
# order should be main__test and debugging option is always set exactly the same as test.  #
#------------------------------------------------------------------------------------------#
declare -a PLONH=(-40.370 -60.209 -52.912 -72.170)
declare -a PLATH=( -9.165  -2.609   5.282  42.540)
declare -a ISOILFLGH=(2 2 2 2)
declare -a NSLCONH=(6 6 6 6)
declare -a SLXCLAYH=(0.052 0.680 0.345 0.060)
declare -a SLXSANDH=(0.821 0.200 0.562 0.660)
declare -a ISOILCOLH=(14 2 14 2)
declare -a SDEPTHH=(C H G E)
declare -a ISOILBCH=(2 1 1 0)
declare -a INTEGRH=(1 3 1 1)
declare -a IBRANCHH=(1 1 1 0)
declare -a IPHYSIOLH=(2 2 2__3)
declare -a IALLOMH=(2 2 3 0)
declare -a IECONH=(0 0 1 0)
declare -a IGRASSH=(0 1 1 0)
declare -a IREPROH=(2 2 2__3 2)
declare -a DECOMPH=(2 2 2__5 0)
declare -a H2OLIMH=(2 2 2__5 1)
declare -a STGROWH=(0 0 1 0)
declare -a PLASTICH=(0 0 2 0)
declare -a IANTHH=(0 0 0 0)
declare -a IFIREH=(0 0 0 0)
declare -a SMFIREH=(-1.40 -1.40 -1.40 0.068)
declare -a PHENOLH=(2 2 -1 1)
declare -a THCRITH=(-1.20 -1.20 -1.20 0.09)
declare -a ICANTURBH=(2 3 0 2)
declare -a ISFCLYRMH=(3 3 4 3)
declare -a IPERCOLH=(0 0 0 1)
declare -a PFTFLAGH=(B B B C)
declare -a TREEFALLH=(0.014 0.014 0.014 0.0)
declare -a IFUSIONH=(0 0 1 0)
declare -a IVEGTDYNH=(1 1 1 1)
declare -a RK4TOLERH=(0.01 0.01 0.01 0.01)
declare -a MAXSITEH=(1 1 1 1)
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#      Variable settings (gridded runs).   In case the test and mainline must be set       #
# differently (e.g., new feature in ED2IN), use | to separate the two entries.  The order  #
# should be main__test and debugging option is always set exactly the same as test.        #
#------------------------------------------------------------------------------------------#
declare -a GRIDRESG=(1.0)
declare -a REGWLONG=(-66.0)
declare -a REGELONG=(-49.0)
declare -a REGSLATG=(-12.0)
declare -a REGNLATG=(  1.0)
declare -a ISOILFLGG=(1)
declare -a ISOILCOLG=(14)
declare -a NSLCONG=(6)
declare -a SLXCLAYG=(-1.00)
declare -a SLXSANDG=(-1.00)
declare -a SDEPTHG=(F)
declare -a ISOILBCG=(2)
declare -a INTEGRG=(1)
declare -a IBRANCHG=(1)
declare -a IPHYSIOLG=(2)
declare -a IALLOMG=(2)
declare -a IECONG=(0)
declare -a IGRASSG=(1)
declare -a IREPROG=(2)
declare -a DECOMPG=(2)
declare -a H2OLIMG=(2)
declare -a STGROWG=(0)
declare -a PLASTICG=(0 )
declare -a IANTHG=(1)
declare -a IFIREG=(2)
declare -a SMFIREG=(-1.40)
declare -a PHENOLG=(2)
declare -a THCRITG=(-1.20)
declare -a ICANTURBG=(3)
declare -a ISFCLYRMG=(3)
declare -a IPERCOLG=(0)
declare -a PFTFLAGG=(B)
declare -a TREEFALLG=(0.014)
declare -a IFUSIONG=(0)
declare -a IVEGTDYNG=(1)
declare -a RK4TOLERG=(0.01)
declare -a MAXSITEG=(1)
#------------------------------------------------------------------------------------------#





#------------------------------------------------------------------------------------------#
#      Additional settings, which depend on whether this is rapid or long tests.           #
#------------------------------------------------------------------------------------------#
case ${TESTTYPE} in
rapid)

   echo " - Performing rapid tests (2 years for POI, 1 year for grid)"

   #----- POI tests will run for two years. -----------------------------------------------#
   declare -a IYEARAS=(2010 2010 2010 2010 2010 2010 2010 2010 2010 2010 2010 2010 2010 2010)
   declare -a IYEARZS=(2012 2012 2012 2012 2012 2012 2012 2012 2012 2012 2012 2012 2012 2012)
   declare -a INITMDS=(5    0    6    0    5    0    5    0    6    5    5    6    0    -1  )
   declare -a RUNTYPS=(INITIAL INITIAL INITIAL INITIAL INITIAL INITIAL INITIAL \
                       INITIAL INITIAL INITIAL INITIAL INITIAL INITIAL INITIAL )
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      High-frequency test will be run for a few days, and save the budget, the Runge-  #
   # Kutta output and the photosynthesis output.                                           #
   #---------------------------------------------------------------------------------------#
   declare -a IYEARAH=(2010 2010 2010 2010)
   declare -a IYEARZH=(2010 2010 2010 2010)
   declare -a IMONTHAH=(10 01 05 02)
   declare -a IMONTHZH=(10 01 05 02)
   declare -a IDATEAH=(21 01 01 01)
   declare -a IDATEZH=(28 08 08 08)
   declare -a D_IYEARAH=(2010 2010 2010 2010)
   declare -a D_IYEARZH=(2010 2010 2010 2010)
   declare -a IDETAILEDH=(7 7 7 7)
   #---------------------------------------------------------------------------------------#


   #----- Gridded tests will run for only 1 year. -----------------------------------------#
   declare -a IYEARAG=(2010)
   declare -a IYEARZG=(2011)
   declare -a INITMDG=(5)
   declare -a RUNTYPG=(INITIAL)
   #---------------------------------------------------------------------------------------#


   #----- Monthly output for state files. -------------------------------------------------#
   UNITSTATE=2
   #---------------------------------------------------------------------------------------#
   ;;

medium)

   echo " - Performing intermediate tests (125 years for POI, 15 years for grid)"

   #----- POI tests will run for 125 years. -----------------------------------------------#
   declare -a IYEARAS=(1975 1975 1975 1975 1975 1975 1975 1975 1975 1975 1975 1975 1975 1975)
   declare -a IYEARZS=(2100 2100 2100 2100 2100 2100 2100 2100 2100 2100 2100 2100 2100 2100)
   declare -a INITMDS=(5    0    6    0    5    0    5    0    6    5    5    6    0    -1  )
   declare -a RUNTYPS=(INITIAL INITIAL INITIAL INITIAL INITIAL INITIAL INITIAL \
                       INITIAL INITIAL INITIAL INITIAL INITIAL INITIAL INITIAL )
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      High-frequency test will be run for 10 years, but we only save the budget        #
   # output.                                                                               #
   #---------------------------------------------------------------------------------------#
   declare -a IYEARAH=(2005 2005 2005 2005)
   declare -a IYEARZH=(2015 2015 2015 2015)
   declare -a IMONTHAH=(01 01 01 01)
   declare -a IMONTHZH=(01 01 01 01)
   declare -a IDATEAH=(01 01 01 01)
   declare -a IDATEZH=(01 01 01 01)
   declare -a D_IYEARAH=(2005 2005 2005 2005)
   declare -a D_IYEARZH=(2006 2006 2006 2006)
   declare -a IDETAILEDH=(1 1 1 1)
   #---------------------------------------------------------------------------------------#



   #----- Gridded tests will run for 15 years. --------------------------------------------#
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
   declare -a IYEARAS=(1800 1800 1800 1800 1800 1800 1800 1800 1800 1800 1800 1800 1800 1800)
   declare -a IYEARZS=(2100 2100 2100 2100 2100 2100 2100 2100 2100 2100 2100 2100 2100 2100)
   declare -a INITMDS=(0    0    0    0    0    0    0    0    0    0    0    0    0    -1  )
   declare -a RUNTYPS=(INITIAL INITIAL INITIAL INITIAL INITIAL INITIAL INITIAL \
                       INITIAL INITIAL INITIAL INITIAL INITIAL INITIAL INITIAL )
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      High-frequency test will be run for 50 years, but we only save the budget        #
   # output.                                                                               #
   #---------------------------------------------------------------------------------------#
   declare -a IYEARAH=(2000 2000 2000 2000)
   declare -a IYEARZH=(2050 2050 2050 2050)
   declare -a IMONTHAH=(01 01 01 01)
   declare -a IMONTHZH=(01 01 01 01)
   declare -a IDATEAH=(01 01 01 01)
   declare -a IDATEZH=(01 01 01 01)
   declare -a D_IYEARAH=(2000 2000 2000 2000)
   declare -a D_IYEARZH=(2001 2001 2001 2001)
   declare -a IDETAILEDH=(1 1 1 1)
   #---------------------------------------------------------------------------------------#



   #----- Gridded tests will run for 35 years. --------------------------------------------#
   declare -a IYEARAG=(1500)
   declare -a IYEARZG=(1550)
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
echo "#!/bin/bash"                  >> ${subbatch}
echo ". ${HOME}/.bashrc ${RC_OPTS}" >> ${subbatch}
#------------------------------------------------------------------------------------------#


#----- Go through all sites, generate command for job queue submission. -------------------#
tcount=0
for i in ${!SITEID[@]}
do
   case ${USE_SITE[i]} in
   y|Y)
      #----- Template and working namelists. ----------------------------------------------#
      TEMPLATEMAIN="${HERE}/Templates/ED2IN-MAIN"
      TEMPLATETEST="${HERE}/Templates/ED2IN-TEST"
      TEMPLATEDBUG="${TEMPLATETEST}"
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



      #------------------------------------------------------------------------------------#
      #     Extract most variable settings.                                                #
      #------------------------------------------------------------------------------------#
      nml_setting ${PLONS[i]}    ; PLON_MAIN=${nl_main}    ; PLON_TEST=${nl_test}
      nml_setting ${PLATS[i]}    ; PLAT_MAIN=${nl_main}    ; PLAT_TEST=${nl_test}
      nml_setting ${ISOILFLGS[i]}; ISOILFLG_MAIN=${nl_main}; ISOILFLG_TEST=${nl_test}
      nml_setting ${NSLCONS[i]}  ; NSLCON_MAIN=${nl_main}  ; NSLCON_TEST=${nl_test}
      nml_setting ${SLXCLAYS[i]} ; SLXCLAY_MAIN=${nl_main} ; SLXCLAY_TEST=${nl_test}
      nml_setting ${SLXSANDS[i]} ; SLXSAND_MAIN=${nl_main} ; SLXSAND_TEST=${nl_test}
      nml_setting ${ISOILCOLS[i]}; ISOILCOL_MAIN=${nl_main}; ISOILCOL_TEST=${nl_test}
      nml_setting ${SDEPTHS[i]}  ; SDEPTH_MAIN=${nl_main}  ; SDEPTH_TEST=${nl_test}
      nml_setting ${ISOILBCS[i]} ; ISOILBC_MAIN=${nl_main} ; ISOILBC_TEST=${nl_test}
      nml_setting ${INTEGRS[i]}  ; INTEGR_MAIN=${nl_main}  ; INTEGR_TEST=${nl_test}
      nml_setting ${IBRANCHS[i]} ; IBRANCH_MAIN=${nl_main} ; IBRANCH_TEST=${nl_test}
      nml_setting ${IALLOMS[i]}  ; IALLOM_MAIN=${nl_main}  ; IALLOM_TEST=${nl_test}
      nml_setting ${IECONS[i]}   ; IECON_MAIN=${nl_main}   ; IECON_TEST=${nl_test}
      nml_setting ${IGRASSS[i]}  ; IGRASS_MAIN=${nl_main}  ; IGRASS_TEST=${nl_test}
      nml_setting ${IREPROS[i]}  ; IREPRO_MAIN=${nl_main}  ; IREPRO_TEST=${nl_test}
      nml_setting ${DECOMPS[i]}  ; DECOMP_MAIN=${nl_main}  ; DECOMP_TEST=${nl_test}
      nml_setting ${H2OLIMS[i]}  ; H2OLIM_MAIN=${nl_main}  ; H2OLIM_TEST=${nl_test}
      nml_setting ${STGROWS[i]}  ; STGROW_MAIN=${nl_main}  ; STGROW_TEST=${nl_test}
      nml_setting ${PLASTICS[i]} ; PLASTIC_MAIN=${nl_main} ; PLASTIC_TEST=${nl_test}
      nml_setting ${IANTHS[i]}   ; IANTH_MAIN=${nl_main}   ; IANTH_TEST=${nl_test}
      nml_setting ${IFIRES[i]}   ; IFIRE_MAIN=${nl_main}   ; IFIRE_TEST=${nl_test}
      nml_setting ${SMFIRES[i]}  ; SMFIRE_MAIN=${nl_main}  ; SMFIRE_TEST=${nl_test}
      nml_setting ${PHENOLS[i]}  ; PHENOL_MAIN=${nl_main}  ; PHENOL_TEST=${nl_test}
      nml_setting ${THCRITS[i]}  ; THCRIT_MAIN=${nl_main}  ; THCRIT_TEST=${nl_test}
      nml_setting ${ICANTURBS[i]}; ICANTURB_MAIN=${nl_main}; ICANTURB_TEST=${nl_test}
      nml_setting ${ISFCLYRMS[i]}; ISFCLYRM_MAIN=${nl_main}; ISFCLYRM_TEST=${nl_test}
      nml_setting ${IPERCOLS[i]} ; IPERCOL_MAIN=${nl_main} ; IPERCOL_TEST=${nl_test}
      nml_setting ${PFTFLAGS[i]} ; PFTFLAG_MAIN=${nl_main} ; PFTFLAG_TEST=${nl_test}
      nml_setting ${TREEFALLS[i]}; TFALL_MAIN=${nl_main}   ; TFALL_TEST=${nl_test}
      nml_setting ${IFUSIONS[i]} ; IFUSION_MAIN=${nl_main} ; IFUSION_TEST=${nl_test}
      nml_setting ${IVEGTDYNS[i]}; IVEGTDYN_MAIN=${nl_main}; IVEGTDYN_TEST=${nl_test}
      nml_setting ${RK4TOLERS[i]}; RK4TOLER_MAIN=${nl_main}; RK4TOLER_TEST=${nl_test}
      nml_setting ${MAXSITES[i]} ; MAXSITE_MAIN=${nl_main} ; MAXSITE_TEST=${nl_test}
      #----- Derived parameters. ----------------------------------------------------------#
      kwd0_setting ${H2OLIM_MAIN}
      KWGRASS_MAIN=${kw_grass}; KWTREE_MAIN=${kw_tree}
      D0GRASS_MAIN=${d0_grass}; D0TREE_MAIN=${d0_tree}
      kwd0_setting ${H2OLIM_TEST}
      KWGRASS_TEST=${kw_grass}; KWTREE_MAIN=${kw_tree}
      D0GRASS_TEST=${d0_grass}; D0TREE_MAIN=${d0_tree}
      #----- Retrieve PFT settings. -------------------------------------------------------#
      pft_setting ${PFTFLAG_MAIN}
      INCPFT_MAIN=${inc_pft}  ; PASTURE_MAIN=${pasture}       ; AGRI_MAIN=${agri}
      PLANT_MAIN=${plantation}
      SLPFT_MAIN=${sl_pft}    ; SLPROBHARV_MAIN=${sl_probharv}; SLMINDBH_MAIN=${sl_mindbh}
      pft_setting ${PFTFLAG_TEST}
      INCPFT_TEST=${inc_pft}; PASTURE_TEST=${pasture}       ; AGRI_TEST=${agri}
      PLANT_TEST=${plantation}
      SLPFT_TEST=${sl_pft}  ; SLPROBHARV_TEST=${sl_probharv}; SLMINDBH_TEST=${sl_mindbh}
      #----- Obtain soil layers. ----------------------------------------------------------#
      soil_setting ${SDEPTH_MAIN}
      NZG_MAIN=${nzg}; SLZ_MAIN=${slz}; SLMSTR_MAIN=${slm}; STGOFF_MAIN=${slt}
      soil_setting ${SDEPTH_TEST}
      NZG_TEST=${nzg}; SLZ_TEST=${slz}; SLMSTR_TEST=${slm}; STGOFF_TEST=${slt}
      #----- Define inital file. ----------------------------------------------------------#
      sfilin_setting ${SITEID[i]} ${INITMDS[i]}; SFILIN=${init_file}
      #----- Define meteorological driver settings. ---------------------------------------#
      metdriv_setting ${SITEID[i]}
      METDRIV=${metdriv}; METCYCA=${metcyca}; METCYCZ=${metcycz}
      IMETAVG=${imetavg}; IMETRAD=${imetrad}
      #----- Define land use driver settings. ---------------------------------------------#
      lu_setting ${SITEID[i]}; LU_DATABASE=${lu_database}
      #----- Define prescribed-phenology settings. ----------------------------------------#
      phenol_setting ${SITEID[i]}
      PHEN_PATH=${phen_path}
      IPHENYSA=${iphenysa}
      IPHENYSZ=${iphenysz}
      IPHENYFA=${iphenyfa}
      IPHENYFZ=${iphenyfz}
      #----- Define parameter file settings. ----------------------------------------------#
      xml_setting ${SITEID[i]}
      EDCNFGF=${edcnfgf}
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Extract most variable settings.                                                #
      #------------------------------------------------------------------------------------#
      #----- Number of ED regions. --------------------------------------------------------#
      sed -i '/NL%N_ED_REGION/c\   NL%N_ED_REGION = 0'                          ${FILEMAIN}
      sed -i '/NL%N_ED_REGION/c\   NL%N_ED_REGION = 0'                          ${FILETEST}
      sed -i '/NL%N_ED_REGION/c\   NL%N_ED_REGION = 0'                          ${FILEDBUG}
      #----- Number of ED Polygons. -------------------------------------------------------#
      sed -i '/NL%N_POI/c\   NL%N_POI = 1'                                      ${FILEMAIN}
      sed -i '/NL%N_POI/c\   NL%N_POI = 1'                                      ${FILETEST}
      sed -i '/NL%N_POI/c\   NL%N_POI = 1'                                      ${FILEDBUG}
      #----- Longitude. -------------------------------------------------------------------#
      sed -i '/NL%POI_LON/c\   NL%POI_LON = '"${PLON_MAIN}"                     ${FILEMAIN}
      sed -i '/NL%POI_LON/c\   NL%POI_LON = '"${PLON_TEST}"                     ${FILETEST}
      sed -i '/NL%POI_LON/c\   NL%POI_LON = '"${PLON_TEST}"                     ${FILEDBUG}
      #----- Latitude. --------------------------------------------------------------------#
      sed -i '/NL%POI_LAT/c\   NL%POI_LAT = '"${PLAT_MAIN}"                     ${FILEMAIN}
      sed -i '/NL%POI_LAT/c\   NL%POI_LAT = '"${PLAT_TEST}"                     ${FILETEST}
      sed -i '/NL%POI_LAT/c\   NL%POI_LAT = '"${PLAT_TEST}"                     ${FILEDBUG}
      #----- Soil texture flag. -----------------------------------------------------------#
      sed -i '/NL%ISOILFLG/c\   NL%ISOILFLG = '"${ISOILFLG_MAIN}"               ${FILEMAIN}
      sed -i '/NL%ISOILFLG/c\   NL%ISOILFLG = '"${ISOILFLG_TEST}"               ${FILETEST}
      sed -i '/NL%ISOILFLG/c\   NL%ISOILFLG = '"${ISOILFLG_TEST}"               ${FILEDBUG}
      #----- Default soil type. -----------------------------------------------------------#
      sed -i '/NL%NSLCON/c\   NL%NSLCON = '"${NSLCON_MAIN}"                     ${FILEMAIN}
      sed -i '/NL%NSLCON/c\   NL%NSLCON = '"${NSLCON_TEST}"                     ${FILETEST}
      sed -i '/NL%NSLCON/c\   NL%NSLCON = '"${NSLCON_TEST}"                     ${FILEDBUG}
      #----- Clay fraction. ---------------------------------------------------------------#
      sed -i '/NL%SLXCLAY/c\   NL%SLXCLAY = '"${SLXCLAY_MAIN}"                  ${FILEMAIN}
      sed -i '/NL%SLXCLAY/c\   NL%SLXCLAY = '"${SLXCLAY_TEST}"                  ${FILETEST}
      sed -i '/NL%SLXCLAY/c\   NL%SLXCLAY = '"${SLXCLAY_TEST}"                  ${FILEDBUG}
      #----- Sand fraction. ---------------------------------------------------------------#
      sed -i '/NL%SLXSAND/c\   NL%SLXSAND = '"${SLXSAND_MAIN}"                  ${FILEMAIN}
      sed -i '/NL%SLXSAND/c\   NL%SLXSAND = '"${SLXSAND_TEST}"                  ${FILETEST}
      sed -i '/NL%SLXSAND/c\   NL%SLXSAND = '"${SLXSAND_TEST}"                  ${FILEDBUG}
      #----- Default soil colour. ---------------------------------------------------------#
      sed -i '/NL%ISOILCOL/c\   NL%ISOILCOL = '"${ISOILCOL_MAIN}"               ${FILEMAIN}
      sed -i '/NL%ISOILCOL/c\   NL%ISOILCOL = '"${ISOILCOL_TEST}"               ${FILETEST}
      sed -i '/NL%ISOILCOL/c\   NL%ISOILCOL = '"${ISOILCOL_TEST}"               ${FILEDBUG}
      #----- Number of soil layers. -------------------------------------------------------#
      sed -i '/NL%NZG/c\   NL%NZG = '"${NZG_MAIN}"                              ${FILEMAIN}
      sed -i '/NL%NZG/c\   NL%NZG = '"${NZG_TEST}"                              ${FILETEST}
      sed -i '/NL%NZG/c\   NL%NZG = '"${NZG_TEST}"                              ${FILEDBUG}
      #----- Soil layers. -----------------------------------------------------------------#
      sed -i '/NL%SLZ/c\   NL%SLZ = '"${SLZ_MAIN}"                              ${FILEMAIN}
      sed -i '/NL%SLZ/c\   NL%SLZ = '"${SLZ_TEST}"                              ${FILETEST}
      sed -i '/NL%SLZ/c\   NL%SLZ = '"${SLZ_TEST}"                              ${FILEDBUG}
      #----- Soil moisture. ---------------------------------------------------------------#
      sed -i '/NL%SLMSTR/c\   NL%SLMSTR = '"${SLMSTR_MAIN}"                     ${FILEMAIN}
      sed -i '/NL%SLMSTR/c\   NL%SLMSTR = '"${SLMSTR_TEST}"                     ${FILETEST}
      sed -i '/NL%SLMSTR/c\   NL%SLMSTR = '"${SLMSTR_TEST}"                     ${FILEDBUG}
      #----- Soil temperature offset. -----------------------------------------------------#
      sed -i '/NL%STGOFF/c\   NL%STGOFF = '"${STGOFF_MAIN}"                     ${FILEMAIN}
      sed -i '/NL%STGOFF/c\   NL%STGOFF = '"${STGOFF_TEST}"                     ${FILETEST}
      sed -i '/NL%STGOFF/c\   NL%STGOFF = '"${STGOFF_TEST}"                     ${FILEDBUG}
      #----- Soil boundary condition. -----------------------------------------------------#
      sed -i '/NL%ISOILBC/c\   NL%ISOILBC = '"${ISOILBC_MAIN}"                  ${FILEMAIN}
      sed -i '/NL%ISOILBC/c\   NL%ISOILBC = '"${ISOILBC_TEST}"                  ${FILETEST}
      sed -i '/NL%ISOILBC/c\   NL%ISOILBC = '"${ISOILBC_TEST}"                  ${FILEDBUG}
      #----- Integration scheme. ----------------------------------------------------------#
      sed -i '/NL%INTEGRATION_SCHEME/c\   NL%INTEGRATION_SCHEME = '"${INTEGR_MAIN}"        \
                                                                                ${FILEMAIN}
      sed -i '/NL%INTEGRATION_SCHEME/c\   NL%INTEGRATION_SCHEME = '"${INTEGR_TEST}"        \
                                                                                ${FILETEST}
      sed -i '/NL%INTEGRATION_SCHEME/c\   NL%INTEGRATION_SCHEME = '"${INTEGR_TEST}"        \
                                                                                ${FILEDBUG}
      #----- Branch themodynamics. --------------------------------------------------------#
      sed -i '/NL%IBRANCH_THERMO/c\   NL%IBRANCH_THERMO = '"${IBRANCH_MAIN}"    ${FILEMAIN}
      sed -i '/NL%IBRANCH_THERMO/c\   NL%IBRANCH_THERMO = '"${IBRANCH_TEST}"    ${FILETEST}
      sed -i '/NL%IBRANCH_THERMO/c\   NL%IBRANCH_THERMO = '"${IBRANCH_TEST}"    ${FILEDBUG}
      #----- Allometry. -------------------------------------------------------------------#
      sed -i '/NL%IALLOM/c\   NL%IALLOM = '"${IALLOM_MAIN}"                     ${FILEMAIN}
      sed -i '/NL%IALLOM/c\   NL%IALLOM = '"${IALLOM_TEST}"                     ${FILETEST}
      sed -i '/NL%IALLOM/c\   NL%IALLOM = '"${IALLOM_TEST}"                     ${FILEDBUG}
      #----- Economic. --------------------------------------------------------------------#
      sed -i '/NL%ECONOMICS_SCHEME/c\   NL%ECONOMICS_SCHEME = '"${IECON_MAIN}"  ${FILEMAIN}
      sed -i '/NL%ECONOMICS_SCHEME/c\   NL%ECONOMICS_SCHEME = '"${IECON_TEST}"  ${FILETEST}
      sed -i '/NL%ECONOMICS_SCHEME/c\   NL%ECONOMICS_SCHEME = '"${IECON_TEST}"  ${FILEDBUG}
      #----- Grass. -----------------------------------------------------------------------#
      sed -i '/NL%IGRASS/c\   NL%IGRASS = '"${IGRASS_MAIN}"                     ${FILEMAIN}
      sed -i '/NL%IGRASS/c\   NL%IGRASS = '"${IGRASS_TEST}"                     ${FILETEST}
      sed -i '/NL%IGRASS/c\   NL%IGRASS = '"${IGRASS_TEST}"                     ${FILEDBUG}
      #----- Reproduction. ----------------------------------------------------------------#
      sed -i '/NL%REPRO_SCHEME/c\   NL%REPRO_SCHEME = '"${IREPRO_MAIN}"         ${FILEMAIN}
      sed -i '/NL%REPRO_SCHEME/c\   NL%REPRO_SCHEME = '"${IREPRO_TEST}"         ${FILETEST}
      sed -i '/NL%REPRO_SCHEME/c\   NL%REPRO_SCHEME = '"${IREPRO_TEST}"         ${FILEDBUG}
      #----- Decomposition. ---------------------------------------------------------------#
      sed -i '/NL%DECOMP_SCHEME/c\   NL%DECOMP_SCHEME = '"${DECOMP_MAIN}"       ${FILEMAIN}
      sed -i '/NL%DECOMP_SCHEME/c\   NL%DECOMP_SCHEME = '"${DECOMP_TEST}"       ${FILETEST}
      sed -i '/NL%DECOMP_SCHEME/c\   NL%DECOMP_SCHEME = '"${DECOMP_TEST}"       ${FILEDBUG}
      #----- Water limitation. ------------------------------------------------------------#
      sed -i '/NL%H2O_PLANT_LIM/c\   NL%H2O_PLANT_LIM = '"${H2OLIM_MAIN}"       ${FILEMAIN}
      sed -i '/NL%H2O_PLANT_LIM/c\   NL%H2O_PLANT_LIM = '"${H2OLIM_TEST}"       ${FILETEST}
      sed -i '/NL%H2O_PLANT_LIM/c\   NL%H2O_PLANT_LIM = '"${H2OLIM_TEST}"       ${FILEDBUG}
      #----- Structural growth. -----------------------------------------------------------#
      sed -i '/NL%ISTRUCT_GROWTH_SCHEME/c\   NL%ISTRUCT_GROWTH_SCHEME = '"${STGROW_MAIN}"  \
                                                                                ${FILEMAIN}
      sed -i '/NL%ISTRUCT_GROWTH_SCHEME/c\   NL%ISTRUCT_GROWTH_SCHEME = '"${STGROW_TEST}"  \
                                                                                ${FILETEST}
      sed -i '/NL%ISTRUCT_GROWTH_SCHEME/c\   NL%ISTRUCT_GROWTH_SCHEME = '"${STGROW_TEST}"  \
                                                                                ${FILEDBUG}
      #----- Trait plasticity. ------------------------------------------------------------#
      sed -i                                                                               \
        '/NL%TRAIT_PLASTICITY_SCHEME/c\   NL%TRAIT_PLASTICITY_SCHEME = '"${PLASTIC_MAIN}"  \
                                                                                ${FILEMAIN}
      sed -i                                                                               \
        '/NL%TRAIT_PLASTICITY_SCHEME/c\   NL%TRAIT_PLASTICITY_SCHEME = '"${PLASTIC_TEST}"  \
                                                                                ${FILETEST}
      sed -i                                                                               \
        '/NL%TRAIT_PLASTICITY_SCHEME/c\   NL%TRAIT_PLASTICITY_SCHEME = '"${PLASTIC_TEST}"  \
                                                                                ${FILEDBUG}
      #----- Anthropogenic disturbance. ---------------------------------------------------#
      sed -i '/NL%INCLUDE_FIRE/c\   NL%INCLUDE_FIRE = '"${IFIRE_MAIN}"          ${FILEMAIN}
      sed -i '/NL%INCLUDE_FIRE/c\   NL%INCLUDE_FIRE = '"${IFIRE_TEST}"          ${FILETEST}
      sed -i '/NL%INCLUDE_FIRE/c\   NL%INCLUDE_FIRE = '"${IFIRE_TEST}"          ${FILEDBUG}
      #----- Fire disturbance. ------------------------------------------------------------#
      sed -i '/NL%INCLUDE_FIRE/c\   NL%INCLUDE_FIRE = '"${IFIRE_MAIN}"          ${FILEMAIN}
      sed -i '/NL%INCLUDE_FIRE/c\   NL%INCLUDE_FIRE = '"${IFIRE_TEST}"          ${FILETEST}
      sed -i '/NL%INCLUDE_FIRE/c\   NL%INCLUDE_FIRE = '"${IFIRE_TEST}"          ${FILEDBUG}
      #----- Soil moisture threshold for fire. --------------------------------------------#
      sed -i '/NL%SM_FIRE/c\   NL%SM_FIRE = '"${SMFIRE_MAIN}"                   ${FILEMAIN}
      sed -i '/NL%SM_FIRE/c\   NL%SM_FIRE = '"${SMFIRE_TEST}"                   ${FILETEST}
      sed -i '/NL%SM_FIRE/c\   NL%SM_FIRE = '"${SMFIRE_TEST}"                   ${FILEDBUG}
      #----- Phenology. -------------------------------------------------------------------#
      sed -i '/NL%IPHEN_SCHEME/c\   NL%IPHEN_SCHEME = '"${PHENOL_MAIN}"         ${FILEMAIN}
      sed -i '/NL%IPHEN_SCHEME/c\   NL%IPHEN_SCHEME = '"${PHENOL_TEST}"         ${FILETEST}
      sed -i '/NL%IPHEN_SCHEME/c\   NL%IPHEN_SCHEME = '"${PHENOL_TEST}"         ${FILEDBUG}
      #----- Soil moisture threshold for drought deciduous. -------------------------------#
      sed -i '/NL%THETACRIT/c\   NL%THETACRIT = '"${THCRIT_MAIN}"               ${FILEMAIN}
      sed -i '/NL%THETACRIT/c\   NL%THETACRIT = '"${THCRIT_TEST}"               ${FILETEST}
      sed -i '/NL%THETACRIT/c\   NL%THETACRIT = '"${THCRIT_TEST}"               ${FILEDBUG}
      #----- Canopy drag parametrisation. -------------------------------------------------#
      sed -i '/NL%ICANTURB/c\   NL%ICANTURB = '"${ICANTURB_MAIN}"               ${FILEMAIN}
      sed -i '/NL%ICANTURB/c\   NL%ICANTURB = '"${ICANTURB_TEST}"               ${FILETEST}
      sed -i '/NL%ICANTURB/c\   NL%ICANTURB = '"${ICANTURB_TEST}"               ${FILEDBUG}
      #----- Canopy turbulence parametrisation. -------------------------------------------#
      sed -i '/NL%ISFCLYRM/c\   NL%ISFCLYRM = '"${ISFCLYRM_MAIN}"               ${FILEMAIN}
      sed -i '/NL%ISFCLYRM/c\   NL%ISFCLYRM = '"${ISFCLYRM_TEST}"               ${FILETEST}
      sed -i '/NL%ISFCLYRM/c\   NL%ISFCLYRM = '"${ISFCLYRM_TEST}"               ${FILEDBUG}
      #----- Surface water percolation parametrisation. -----------------------------------#
      sed -i '/NL%IPERCOL/c\   NL%IPERCOL = '"${IPERCOL_MAIN}"                  ${FILEMAIN}
      sed -i '/NL%IPERCOL/c\   NL%IPERCOL = '"${IPERCOL_TEST}"                  ${FILETEST}
      sed -i '/NL%IPERCOL/c\   NL%IPERCOL = '"${IPERCOL_TEST}"                  ${FILEDBUG}
      #----- Tree fall disturbance rate. --------------------------------------------------#
      sed -i                                                                               \
       '/NL%TREEFALL_DISTURBANCE_RATE/c\   NL%TREEFALL_DISTURBANCE_RATE = '"${TFALL_MAIN}" \
                                                                                ${FILEMAIN}
      sed -i                                                                               \
       '/NL%TREEFALL_DISTURBANCE_RATE/c\   NL%TREEFALL_DISTURBANCE_RATE = '"${TFALL_TEST}" \
                                                                                ${FILETEST}
      sed -i                                                                               \
       '/NL%TREEFALL_DISTURBANCE_RATE/c\   NL%TREEFALL_DISTURBANCE_RATE = '"${TFALL_TEST}" \
                                                                                ${FILEDBUG}
      #----- Patch/cohort fusion method. --------------------------------------------------#
      sed -i '/NL%IFUSION/c\   NL%IFUSION = '"${IFUSION_MAIN}"                  ${FILEMAIN}
      sed -i '/NL%IFUSION/c\   NL%IFUSION = '"${IFUSION_TEST}"                  ${FILETEST}
      sed -i '/NL%IFUSION/c\   NL%IFUSION = '"${IFUSION_TEST}"                  ${FILEDBUG}
      #----- Vegetation dynamics flag. ----------------------------------------------------#
      sed -i '/NL%IVEGT_DYNAMICS/c\   NL%IVEGT_DYNAMICS = '"${IVEGTDYN_MAIN}"   ${FILEMAIN}
      sed -i '/NL%IVEGT_DYNAMICS/c\   NL%IVEGT_DYNAMICS = '"${IVEGTDYN_TEST}"   ${FILETEST}
      sed -i '/NL%IVEGT_DYNAMICS/c\   NL%IVEGT_DYNAMICS = '"${IVEGTDYN_TEST}"   ${FILEDBUG}
      #----- Error tolerance for Runge-Kutta integrator. ----------------------------------#
      sed -i '/NL%RK4_TOLERANCE/c\   NL%RK4_TOLERANCE = '"${RK4TOLER_MAIN}"     ${FILEMAIN}
      sed -i '/NL%RK4_TOLERANCE/c\   NL%RK4_TOLERANCE = '"${RK4TOLER_TEST}"     ${FILETEST}
      sed -i '/NL%RK4_TOLERANCE/c\   NL%RK4_TOLERANCE = '"${RK4TOLER_TEST}"     ${FILEDBUG}
      #----- Maximum number of sites. -----------------------------------------------------#
      sed -i '/NL%MAXSITE/c\   NL%MAXSITE = '"${MAXSITE_MAIN}"                  ${FILEMAIN}
      sed -i '/NL%MAXSITE/c\   NL%MAXSITE = '"${MAXSITE_TEST}"                  ${FILETEST}
      sed -i '/NL%MAXSITE/c\   NL%MAXSITE = '"${MAXSITE_TEST}"                  ${FILEDBUG}
      #----- Root conductance (grass). ----------------------------------------------------#
      sed -i '/NL%KW_GRASS/c\   NL%KW_GRASS = '"${KWGRASS_MAIN}"                ${FILEMAIN}
      sed -i '/NL%KW_GRASS/c\   NL%KW_GRASS = '"${KWGRASS_TEST}"                ${FILETEST}
      sed -i '/NL%KW_GRASS/c\   NL%KW_GRASS = '"${KWGRASS_TEST}"                ${FILEDBUG}
      #----- Root conductance (grass). ----------------------------------------------------#
      sed -i '/NL%KW_TREE/c\   NL%KW_TREE = '"${KWGRASS_MAIN}"                  ${FILEMAIN}
      sed -i '/NL%KW_TREE/c\   NL%KW_TREE = '"${KWGRASS_TEST}"                  ${FILETEST}
      sed -i '/NL%KW_TREE/c\   NL%KW_TREE = '"${KWGRASS_TEST}"                  ${FILEDBUG}
      #----- VPD threshold (grass). -------------------------------------------------------#
      sed -i '/NL%D0_GRASS/c\   NL%D0_GRASS = '"${D0GRASS_MAIN}"                ${FILEMAIN}
      sed -i '/NL%D0_GRASS/c\   NL%D0_GRASS = '"${D0GRASS_TEST}"                ${FILETEST}
      sed -i '/NL%D0_GRASS/c\   NL%D0_GRASS = '"${D0GRASS_TEST}"                ${FILEDBUG}
      #----- VPD threshold (grass). -------------------------------------------------------#
      sed -i '/NL%D0_TREE/c\   NL%D0_TREE = '"${D0GRASS_MAIN}"                  ${FILEMAIN}
      sed -i '/NL%D0_TREE/c\   NL%D0_TREE = '"${D0GRASS_TEST}"                  ${FILETEST}
      sed -i '/NL%D0_TREE/c\   NL%D0_TREE = '"${D0GRASS_TEST}"                  ${FILEDBUG}
      #----- Initial condition file (it may be overwritten if this is HISTORY run). -------#
      sed -i '/NL%SFILIN/c\   NL%SFILIN = '\'"${SFILIN}"\'                      ${FILEMAIN}
      sed -i '/NL%SFILIN/c\   NL%SFILIN = '\'"${SFILIN}"\'                      ${FILETEST}
      sed -i '/NL%SFILIN/c\   NL%SFILIN = '\'"${SFILIN}"\'                      ${FILEDBUG}
      #----- Meteorological driver. -------------------------------------------------------#
      sed -i '/NL%ED_MET_DRIVER_DB/c\   NL%ED_MET_DRIVER_DB = '\'"${METDRIV}"\' ${FILEMAIN}
      sed -i '/NL%ED_MET_DRIVER_DB/c\   NL%ED_MET_DRIVER_DB = '\'"${METDRIV}"\' ${FILETEST}
      sed -i '/NL%ED_MET_DRIVER_DB/c\   NL%ED_MET_DRIVER_DB = '\'"${METDRIV}"\' ${FILEDBUG}
      #----- First year of meteorological driver. -----------------------------------------#
      sed -i '/NL%METCYC1/c\   NL%METCYC1 = '"${METCYCA}"                       ${FILEMAIN}
      sed -i '/NL%METCYC1/c\   NL%METCYC1 = '"${METCYCA}"                       ${FILETEST}
      sed -i '/NL%METCYC1/c\   NL%METCYC1 = '"${METCYCA}"                       ${FILEDBUG}
      #----- Last year of meteorological driver. ------------------------------------------#
      sed -i '/NL%METCYCF/c\   NL%METCYCF = '"${METCYCZ}"                       ${FILEMAIN}
      sed -i '/NL%METCYCF/c\   NL%METCYCF = '"${METCYCZ}"                       ${FILETEST}
      sed -i '/NL%METCYCF/c\   NL%METCYCF = '"${METCYCZ}"                       ${FILEDBUG}
      #----- Time-averaging flag (meteorological driver). ---------------------------------#
      sed -i '/NL%IMETAVG/c\   NL%IMETAVG = '"${IMETAVG}"                       ${FILEMAIN}
      sed -i '/NL%IMETAVG/c\   NL%IMETAVG = '"${IMETAVG}"                       ${FILETEST}
      sed -i '/NL%IMETAVG/c\   NL%IMETAVG = '"${IMETAVG}"                       ${FILEDBUG}
      #----- Radiation decomposition scheme. ----------------------------------------------#
      sed -i '/NL%IMETRAD/c\   NL%IMETRAD = '"${IMETRAD}"                       ${FILEMAIN}
      sed -i '/NL%IMETRAD/c\   NL%IMETRAD = '"${IMETRAD}"                       ${FILETEST}
      sed -i '/NL%IMETRAD/c\   NL%IMETRAD = '"${IMETRAD}"                       ${FILEDBUG}
      #----- PFTs to include in this simulation. ------------------------------------------#
      sed -i '/NL%INCLUDE_THESE_PFT/c\   NL%INCLUDE_THESE_PFT = '"${INCPFT_MAIN}"          \
                                                                                ${FILEMAIN}
      sed -i '/NL%INCLUDE_THESE_PFT/c\   NL%INCLUDE_THESE_PFT = '"${INCPFT_TEST}"          \
                                                                                ${FILETEST}
      sed -i '/NL%INCLUDE_THESE_PFT/c\   NL%INCLUDE_THESE_PFT = '"${INCPFT_TEST}"          \
                                                                                ${FILEDBUG}
      #----- PFT for pastures. ------------------------------------------------------------#
      sed -i '/NL%PASTURE_STOCK/c\   NL%PASTURE_STOCK = '"${PASTURE_MAIN}"      ${FILEMAIN}
      sed -i '/NL%PASTURE_STOCK/c\   NL%PASTURE_STOCK = '"${PASTURE_TEST}"      ${FILETEST}
      sed -i '/NL%PASTURE_STOCK/c\   NL%PASTURE_STOCK = '"${PASTURE_TEST}"      ${FILEDBUG}
      #----- PFT for croplands. -----------------------------------------------------------#
      sed -i '/NL%AGRI_STOCK/c\   NL%AGRI_STOCK = '"${AGRI_MAIN}"               ${FILEMAIN}
      sed -i '/NL%AGRI_STOCK/c\   NL%AGRI_STOCK = '"${AGRI_TEST}"               ${FILETEST}
      sed -i '/NL%AGRI_STOCK/c\   NL%AGRI_STOCK = '"${AGRI_TEST}"               ${FILEDBUG}
      #----- PFT for forest plantation. ---------------------------------------------------#
      sed -i '/NL%PLANTATION_STOCK/c\   NL%PLANTATION_STOCK = '"${PLANT_MAIN}"  ${FILEMAIN}
      sed -i '/NL%PLANTATION_STOCK/c\   NL%PLANTATION_STOCK = '"${PLANT_TEST}"  ${FILETEST}
      sed -i '/NL%PLANTATION_STOCK/c\   NL%PLANTATION_STOCK = '"${PLANT_TEST}"  ${FILEDBUG}
      #----- PFT for selective logging. ---------------------------------------------------#
      sed -i '/NL%SL_PFT/c\   NL%SL_PFT = '"${SLPFT_MAIN}"                      ${FILEMAIN}
      sed -i '/NL%SL_PFT/c\   NL%SL_PFT = '"${SLPFT_TEST}"                      ${FILETEST}
      sed -i '/NL%SL_PFT/c\   NL%SL_PFT = '"${SLPFT_TEST}"                      ${FILEDBUG}
      #----- Harvesting probability for selective logging. --------------------------------#
      sed -i '/NL%SL_PROB_HARVEST/c\   NL%SL_PROB_HARVEST = '"${SLPROBHARV_MAIN}"          \
                                                                                ${FILEMAIN}
      sed -i '/NL%SL_PROB_HARVEST/c\   NL%SL_PROB_HARVEST = '"${SLPROBHARV_TEST}"          \
                                                                                ${FILETEST}
      sed -i '/NL%SL_PROB_HARVEST/c\   NL%SL_PROB_HARVEST = '"${SLPROBHARV_TEST}"          \
                                                                                ${FILEDBUG}
      #----- Minimum DBH for selective logging. -------------------------------------------#
      sed -i '/NL%SL_MINDBH_HARVEST/c\   NL%SL_MINDBH_HARVEST = '"${SLMINDBH_MAIN}"        \
                                                                                ${FILEMAIN}
      sed -i '/NL%SL_MINDBH_HARVEST/c\   NL%SL_MINDBH_HARVEST = '"${SLMINDBH_TEST}"        \
                                                                                ${FILETEST}
      sed -i '/NL%SL_MINDBH_HARVEST/c\   NL%SL_MINDBH_HARVEST = '"${SLMINDBH_TEST}"        \
                                                                                ${FILEDBUG}
      #----- Land use data base. ----------------------------------------------------------#
      sed -i '/NL%LU_DATABASE/c\   NL%LU_DATABASE = '\'"${LU_DATABASE}"\'       ${FILEMAIN}
      sed -i '/NL%LU_DATABASE/c\   NL%LU_DATABASE = '\'"${LU_DATABASE}"\'       ${FILETEST}
      sed -i '/NL%LU_DATABASE/c\   NL%LU_DATABASE = '\'"${LU_DATABASE}"\'       ${FILEDBUG}
      #----- Phenology data base. ---------------------------------------------------------#
      sed -i '/NL%PHENPATH/c\   NL%PHENPATH = '\'"${PHEN_PATH}"\'               ${FILEMAIN}
      sed -i '/NL%PHENPATH/c\   NL%PHENPATH = '\'"${PHEN_PATH}"\'               ${FILETEST}
      sed -i '/NL%PHENPATH/c\   NL%PHENPATH = '\'"${PHEN_PATH}"\'               ${FILEDBUG}
      #----- First year of the phenology (Spring). ----------------------------------------#
      sed -i '/NL%IPHENYS1/c\   NL%IPHENYS1 = '"${IPHENYSA}"                    ${FILEMAIN}
      sed -i '/NL%IPHENYS1/c\   NL%IPHENYS1 = '"${IPHENYSA}"                    ${FILETEST}
      sed -i '/NL%IPHENYS1/c\   NL%IPHENYS1 = '"${IPHENYSA}"                    ${FILEDBUG}
      #----- Last year of the phenology (Spring). -----------------------------------------#
      sed -i '/NL%IPHENYSF/c\   NL%IPHENYSF = '"${IPHENYSZ}"                    ${FILEMAIN}
      sed -i '/NL%IPHENYSF/c\   NL%IPHENYSF = '"${IPHENYSZ}"                    ${FILETEST}
      sed -i '/NL%IPHENYSF/c\   NL%IPHENYSF = '"${IPHENYSZ}"                    ${FILEDBUG}
      #----- First year of the phenology (Autumn). ----------------------------------------#
      sed -i '/NL%IPHENYF1/c\   NL%IPHENYF1 = '"${IPHENYFA}"                    ${FILEMAIN}
      sed -i '/NL%IPHENYF1/c\   NL%IPHENYF1 = '"${IPHENYFA}"                    ${FILETEST}
      sed -i '/NL%IPHENYF1/c\   NL%IPHENYF1 = '"${IPHENYFA}"                    ${FILEDBUG}
      #----- Last year of the phenology (Autumn). -----------------------------------------#
      sed -i '/NL%IPHENYFF/c\   NL%IPHENYFF = '"${IPHENYFZ}"                    ${FILEMAIN}
      sed -i '/NL%IPHENYFF/c\   NL%IPHENYFF = '"${IPHENYFZ}"                    ${FILETEST}
      sed -i '/NL%IPHENYFF/c\   NL%IPHENYFF = '"${IPHENYFZ}"                    ${FILEDBUG}
      #----- XML file. --------------------------------------------------------------------#
      sed -i '/NL%IEDCNFGF/c\   NL%IEDCNFGF = '\'"${EDCNFGF}"\'                 ${FILEMAIN}
      sed -i '/NL%IEDCNFGF/c\   NL%IEDCNFGF = '\'"${EDCNFGF}"\'                 ${FILETEST}
      sed -i '/NL%IEDCNFGF/c\   NL%IEDCNFGF = '\'"${EDCNFGF}"\'                 ${FILEDBUG}
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #      Detailed output settings.                                                     #
      #------------------------------------------------------------------------------------#
      #----- Detailed output list. --------------------------------------------------------#
      sed -i '/NL%IDETAILED/c\   NL%IDETAILED = 0'                          ${FILEMAIN}
      sed -i '/NL%IDETAILED/c\   NL%IDETAILED = 0'                          ${FILETEST}
      sed -i '/NL%IDETAILED/c\   NL%IDETAILED = 0'                          ${FILEDBUG}
      #----- Patches to keep. -------------------------------------------------------------#
      sed -i '/NL%PATCH_KEEP/c\   NL%PATCH_KEEP = 0'                        ${FILEMAIN}
      sed -i '/NL%PATCH_KEEP/c\   NL%PATCH_KEEP = 0'                        ${FILETEST}
      sed -i '/NL%PATCH_KEEP/c\   NL%PATCH_KEEP = 0'                        ${FILEDBUG}
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


      #----- Modify the ED2IN start times. ------------------------------------------------#
      sed -i '/NL%IYEARA/c\   NL%IYEARA   = '${IYEARAS[i]}    ${FILEMAIN}
      sed -i '/NL%IYEARA/c\   NL%IYEARA   = '${IYEARAS[i]}    ${FILETEST}
      sed -i '/NL%IYEARA/c\   NL%IYEARA   = '${D_IYEARAS[i]}  ${FILEDBUG}
      sed -i '/NL%IMONTHA/c\   NL%IMONTHA   = '${IMONTHAS[i]} ${FILEMAIN}
      sed -i '/NL%IMONTHA/c\   NL%IMONTHA   = '${IMONTHAS[i]} ${FILETEST}
      sed -i '/NL%IMONTHA/c\   NL%IMONTHA   = '${IMONTHAS[i]} ${FILEDBUG}
      sed -i '/NL%IDATEA/c\   NL%IDATEA   = '${IDATEAS[i]}    ${FILEMAIN}
      sed -i '/NL%IDATEA/c\   NL%IDATEA   = '${IDATEAS[i]}    ${FILETEST}
      sed -i '/NL%IDATEA/c\   NL%IDATEA   = '${IDATEAS[i]}    ${FILEDBUG}
      sed -i '/NL%ITIMEA/c\   NL%ITIMEA   = 0000'             ${FILEMAIN}
      sed -i '/NL%ITIMEA/c\   NL%ITIMEA   = 0000'             ${FILETEST}
      sed -i '/NL%ITIMEA/c\   NL%ITIMEA   = 0000'             ${FILEDBUG}
      #------------------------------------------------------------------------------------#


      #----- Modify the ED2IN end years. --------------------------------------------------#
      sed -i '/NL%IYEARZ/c\   NL%IYEARZ   = '${IYEARZS[i]}    ${FILEMAIN}
      sed -i '/NL%IYEARZ/c\   NL%IYEARZ   = '${IYEARZS[i]}    ${FILETEST}
      sed -i '/NL%IYEARZ/c\   NL%IYEARZ   = '${D_IYEARZS[i]}  ${FILEDBUG}
      sed -i '/NL%IMONTHZ/c\   NL%IMONTHZ   = '${IMONTHZS[i]} ${FILEMAIN}
      sed -i '/NL%IMONTHZ/c\   NL%IMONTHZ   = '${IMONTHZS[i]} ${FILETEST}
      sed -i '/NL%IMONTHZ/c\   NL%IMONTHZ   = '${IMONTHZS[i]} ${FILEDBUG}
      sed -i '/NL%IDATEZ/c\   NL%IDATEZ   = '${IDATEZS[i]}    ${FILEMAIN}
      sed -i '/NL%IDATEZ/c\   NL%IDATEZ   = '${IDATEZS[i]}    ${FILETEST}
      sed -i '/NL%IDATEZ/c\   NL%IDATEZ   = '${IDATEZS[i]}    ${FILEDBUG}
      sed -i '/NL%ITIMEZ/c\   NL%ITIMEZ   = 0000'             ${FILEMAIN}
      sed -i '/NL%ITIMEZ/c\   NL%ITIMEZ   = 0000'             ${FILETEST}
      sed -i '/NL%ITIMEZ/c\   NL%ITIMEZ   = 0000'             ${FILEDBUG}
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
         jobwrap=". ${HOME}/.bashrc ${RC_OPTS}; cd ${HERE}/${VERSION}"
         jobwrap="${jobwrap}; export OMP_NUM_THREADS=${SITECPU[i]}"
         jobwrap="${jobwrap}; ${LNK_MAIN_EXE} -f ${FILEMAIN}"
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
         jobwrap=". ${HOME}/.bashrc ${RC_OPTS}; cd ${HERE}/${VERSION}"
         jobwrap="${jobwrap}; export OMP_NUM_THREADS=${SITECPU[i]}"
         jobwrap="${jobwrap}; ${LNK_TEST_EXE} -f ${FILETEST}"
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
         jobwrap=". ${HOME}/.bashrc ${RC_OPTS}; cd ${HERE}/${VERSION}"
         jobwrap="${jobwrap}; export OMP_NUM_THREADS=${SITECPU[i]}"
         jobwrap="${jobwrap}; ${LNK_DBUG_EXE} -f ${FILEDBUG}"
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
      TEMPLATEMAIN="${HERE}/Templates/ED2IN-MAIN"
      TEMPLATETEST="${HERE}/Templates/ED2IN-TEST"
      TEMPLATEDBUG="${TEMPLATETEST}"
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



      #------------------------------------------------------------------------------------#
      #     Extract most variable settings.                                                #
      #------------------------------------------------------------------------------------#
      nml_setting ${PLONH[i]}    ; PLON_MAIN=${nl_main}    ; PLON_TEST=${nl_test}
      nml_setting ${PLATH[i]}    ; PLAT_MAIN=${nl_main}    ; PLAT_TEST=${nl_test}
      nml_setting ${ISOILFLGH[i]}; ISOILFLG_MAIN=${nl_main}; ISOILFLG_TEST=${nl_test}
      nml_setting ${NSLCONH[i]}  ; NSLCON_MAIN=${nl_main}  ; NSLCON_TEST=${nl_test}
      nml_setting ${SLXCLAYH[i]} ; SLXCLAY_MAIN=${nl_main} ; SLXCLAY_TEST=${nl_test}
      nml_setting ${SLXSANDH[i]} ; SLXSAND_MAIN=${nl_main} ; SLXSAND_TEST=${nl_test}
      nml_setting ${ISOILCOLH[i]}; ISOILCOL_MAIN=${nl_main}; ISOILCOL_TEST=${nl_test}
      nml_setting ${SDEPTHH[i]}  ; SDEPTH_MAIN=${nl_main}  ; SDEPTH_TEST=${nl_test}
      nml_setting ${ISOILBCH[i]} ; ISOILBC_MAIN=${nl_main} ; ISOILBC_TEST=${nl_test}
      nml_setting ${INTEGRH[i]}  ; INTEGR_MAIN=${nl_main}  ; INTEGR_TEST=${nl_test}
      nml_setting ${IBRANCHH[i]} ; IBRANCH_MAIN=${nl_main} ; IBRANCH_TEST=${nl_test}
      nml_setting ${IALLOMH[i]}  ; IALLOM_MAIN=${nl_main}  ; IALLOM_TEST=${nl_test}
      nml_setting ${IECONH[i]}   ; IECON_MAIN=${nl_main}   ; IECON_TEST=${nl_test}
      nml_setting ${IGRASSH[i]}  ; IGRASS_MAIN=${nl_main}  ; IGRASS_TEST=${nl_test}
      nml_setting ${IREPROH[i]}  ; IREPRO_MAIN=${nl_main}  ; IREPRO_TEST=${nl_test}
      nml_setting ${DECOMPH[i]}  ; DECOMP_MAIN=${nl_main}  ; DECOMP_TEST=${nl_test}
      nml_setting ${H2OLIMH[i]}  ; H2OLIM_MAIN=${nl_main}  ; H2OLIM_TEST=${nl_test}
      nml_setting ${STGROWH[i]}  ; STGROW_MAIN=${nl_main}  ; STGROW_TEST=${nl_test}
      nml_setting ${PLASTICH[i]} ; PLASTIC_MAIN=${nl_main} ; PLASTIC_TEST=${nl_test}
      nml_setting ${IANTHH[i]}   ; IANTH_MAIN=${nl_main}   ; IANTH_TEST=${nl_test}
      nml_setting ${IFIREH[i]}   ; IFIRE_MAIN=${nl_main}   ; IFIRE_TEST=${nl_test}
      nml_setting ${SMFIREH[i]}  ; SMFIRE_MAIN=${nl_main}  ; SMFIRE_TEST=${nl_test}
      nml_setting ${PHENOLH[i]}  ; PHENOL_MAIN=${nl_main}  ; PHENOL_TEST=${nl_test}
      nml_setting ${THCRITH[i]}  ; THCRIT_MAIN=${nl_main}  ; THCRIT_TEST=${nl_test}
      nml_setting ${ICANTURBH[i]}; ICANTURB_MAIN=${nl_main}; ICANTURB_TEST=${nl_test}
      nml_setting ${ISFCLYRMH[i]}; ISFCLYRM_MAIN=${nl_main}; ISFCLYRM_TEST=${nl_test}
      nml_setting ${IPERCOLH[i]} ; IPERCOL_MAIN=${nl_main} ; IPERCOL_TEST=${nl_test}
      nml_setting ${PFTFLAGH[i]} ; PFTFLAG_MAIN=${nl_main} ; PFTFLAG_TEST=${nl_test}
      nml_setting ${TREEFALLH[i]}; TFALL_MAIN=${nl_main}   ; TFALL_TEST=${nl_test}
      nml_setting ${IFUSIONH[i]} ; IFUSION_MAIN=${nl_main} ; IFUSION_TEST=${nl_test}
      nml_setting ${IVEGTDYNH[i]}; IVEGTDYN_MAIN=${nl_main}; IVEGTDYN_TEST=${nl_test}
      nml_setting ${RK4TOLERH[i]}; RK4TOLER_MAIN=${nl_main}; RK4TOLER_TEST=${nl_test}
      nml_setting ${MAXSITEH[i]} ; MAXSITE_MAIN=${nl_main} ; MAXSITE_TEST=${nl_test}
      #----- Derived parameters. ----------------------------------------------------------#
      kwd0_setting ${H2OLIM_MAIN}
      KWGRASS_MAIN=${kw_grass}; KWTREE_MAIN=${kw_tree}
      D0GRASS_MAIN=${d0_grass}; D0TREE_MAIN=${d0_tree}
      kwd0_setting ${H2OLIM_TEST}
      KWGRASS_TEST=${kw_grass}; KWTREE_MAIN=${kw_tree}
      D0GRASS_TEST=${d0_grass}; D0TREE_MAIN=${d0_tree}
      #----- Retrieve PFT settings. -------------------------------------------------------#
      pft_setting ${PFTFLAG_MAIN}
      INCPFT_MAIN=${inc_pft}  ; PASTURE_MAIN=${pasture}       ; AGRI_MAIN=${agri}
      PLANT_MAIN=${plantation}
      SLPFT_MAIN=${sl_pft}    ; SLPROBHARV_MAIN=${sl_probharv}; SLMINDBH_MAIN=${sl_mindbh}
      pft_setting ${PFTFLAG_TEST}
      INCPFT_TEST=${inc_pft}; PASTURE_TEST=${pasture}       ; AGRI_TEST=${agri}
      PLANT_TEST=${plantation}
      SLPFT_TEST=${sl_pft}  ; SLPROBHARV_TEST=${sl_probharv}; SLMINDBH_TEST=${sl_mindbh}
      #----- Obtain soil layers. ----------------------------------------------------------#
      soil_setting ${SDEPTH_MAIN}
      NZG_MAIN=${nzg}; SLZ_MAIN=${slz}; SLMSTR_MAIN=${slm}; STGOFF_MAIN=${slt}
      soil_setting ${SDEPTH_TEST}
      NZG_TEST=${nzg}; SLZ_TEST=${slz}; SLMSTR_TEST=${slm}; STGOFF_TEST=${slt}
      #----- Define inital file. ----------------------------------------------------------#
      sfilin_setting ${HIFRID[i]} ${INITMDH[i]}; SFILIN=${init_file}
      #----- Define meteorological driver settings. ---------------------------------------#
      metdriv_setting ${HIFRID[i]}
      METDRIV=${metdriv}; METCYCA=${metcyca}; METCYCZ=${metcycz}
      IMETAVG=${imetavg}; IMETRAD=${imetrad}
      #----- Define land use driver settings. ---------------------------------------------#
      lu_setting ${HIFRID[i]}; LU_DATABASE=${lu_database}
      #----- Define prescribed-phenology settings. ----------------------------------------#
      phenol_setting ${HIFRID[i]}
      PHEN_PATH=${phen_path}
      IPHENYSA=${iphenysa}
      IPHENYSZ=${iphenysz}
      IPHENYFA=${iphenyfa}
      IPHENYFZ=${iphenyfz}
      #----- Define parameter file settings. ----------------------------------------------#
      xml_setting ${HIFRID[i]}
      EDCNFGF=${edcnfgf}
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Extract most variable settings.                                                #
      #------------------------------------------------------------------------------------#
      #----- Number of ED regions. --------------------------------------------------------#
      sed -i '/NL%N_ED_REGION/c\   NL%N_ED_REGION = 0'                          ${FILEMAIN}
      sed -i '/NL%N_ED_REGION/c\   NL%N_ED_REGION = 0'                          ${FILETEST}
      sed -i '/NL%N_ED_REGION/c\   NL%N_ED_REGION = 0'                          ${FILEDBUG}
      #----- Number of ED Polygons. -------------------------------------------------------#
      sed -i '/NL%N_POI/c\   NL%N_POI = 1'                                      ${FILEMAIN}
      sed -i '/NL%N_POI/c\   NL%N_POI = 1'                                      ${FILETEST}
      sed -i '/NL%N_POI/c\   NL%N_POI = 1'                                      ${FILEDBUG}
      #----- Longitude. -------------------------------------------------------------------#
      sed -i '/NL%POI_LON/c\   NL%POI_LON = '"${PLON_MAIN}"                     ${FILEMAIN}
      sed -i '/NL%POI_LON/c\   NL%POI_LON = '"${PLON_TEST}"                     ${FILETEST}
      sed -i '/NL%POI_LON/c\   NL%POI_LON = '"${PLON_TEST}"                     ${FILEDBUG}
      #----- Latitude. --------------------------------------------------------------------#
      sed -i '/NL%POI_LAT/c\   NL%POI_LAT = '"${PLAT_MAIN}"                     ${FILEMAIN}
      sed -i '/NL%POI_LAT/c\   NL%POI_LAT = '"${PLAT_TEST}"                     ${FILETEST}
      sed -i '/NL%POI_LAT/c\   NL%POI_LAT = '"${PLAT_TEST}"                     ${FILEDBUG}
      #----- Soil texture flag. -----------------------------------------------------------#
      sed -i '/NL%ISOILFLG/c\   NL%ISOILFLG = '"${ISOILFLG_MAIN}"               ${FILEMAIN}
      sed -i '/NL%ISOILFLG/c\   NL%ISOILFLG = '"${ISOILFLG_TEST}"               ${FILETEST}
      sed -i '/NL%ISOILFLG/c\   NL%ISOILFLG = '"${ISOILFLG_TEST}"               ${FILEDBUG}
      #----- Default soil type. -----------------------------------------------------------#
      sed -i '/NL%NSLCON/c\   NL%NSLCON = '"${NSLCON_MAIN}"                     ${FILEMAIN}
      sed -i '/NL%NSLCON/c\   NL%NSLCON = '"${NSLCON_TEST}"                     ${FILETEST}
      sed -i '/NL%NSLCON/c\   NL%NSLCON = '"${NSLCON_TEST}"                     ${FILEDBUG}
      #----- Clay fraction. ---------------------------------------------------------------#
      sed -i '/NL%SLXCLAY/c\   NL%SLXCLAY = '"${SLXCLAY_MAIN}"                  ${FILEMAIN}
      sed -i '/NL%SLXCLAY/c\   NL%SLXCLAY = '"${SLXCLAY_TEST}"                  ${FILETEST}
      sed -i '/NL%SLXCLAY/c\   NL%SLXCLAY = '"${SLXCLAY_TEST}"                  ${FILEDBUG}
      #----- Sand fraction. ---------------------------------------------------------------#
      sed -i '/NL%SLXSAND/c\   NL%SLXSAND = '"${SLXSAND_MAIN}"                  ${FILEMAIN}
      sed -i '/NL%SLXSAND/c\   NL%SLXSAND = '"${SLXSAND_TEST}"                  ${FILETEST}
      sed -i '/NL%SLXSAND/c\   NL%SLXSAND = '"${SLXSAND_TEST}"                  ${FILEDBUG}
      #----- Default soil colour. ---------------------------------------------------------#
      sed -i '/NL%ISOILCOL/c\   NL%ISOILCOL = '"${ISOILCOL_MAIN}"               ${FILEMAIN}
      sed -i '/NL%ISOILCOL/c\   NL%ISOILCOL = '"${ISOILCOL_TEST}"               ${FILETEST}
      sed -i '/NL%ISOILCOL/c\   NL%ISOILCOL = '"${ISOILCOL_TEST}"               ${FILEDBUG}
      #----- Number of soil layers. -------------------------------------------------------#
      sed -i '/NL%NZG/c\   NL%NZG = '"${NZG_MAIN}"                              ${FILEMAIN}
      sed -i '/NL%NZG/c\   NL%NZG = '"${NZG_TEST}"                              ${FILETEST}
      sed -i '/NL%NZG/c\   NL%NZG = '"${NZG_TEST}"                              ${FILEDBUG}
      #----- Soil layers. -----------------------------------------------------------------#
      sed -i '/NL%SLZ/c\   NL%SLZ = '"${SLZ_MAIN}"                              ${FILEMAIN}
      sed -i '/NL%SLZ/c\   NL%SLZ = '"${SLZ_TEST}"                              ${FILETEST}
      sed -i '/NL%SLZ/c\   NL%SLZ = '"${SLZ_TEST}"                              ${FILEDBUG}
      #----- Soil moisture. ---------------------------------------------------------------#
      sed -i '/NL%SLMSTR/c\   NL%SLMSTR = '"${SLMSTR_MAIN}"                     ${FILEMAIN}
      sed -i '/NL%SLMSTR/c\   NL%SLMSTR = '"${SLMSTR_TEST}"                     ${FILETEST}
      sed -i '/NL%SLMSTR/c\   NL%SLMSTR = '"${SLMSTR_TEST}"                     ${FILEDBUG}
      #----- Soil temperature offset. -----------------------------------------------------#
      sed -i '/NL%STGOFF/c\   NL%STGOFF = '"${STGOFF_MAIN}"                     ${FILEMAIN}
      sed -i '/NL%STGOFF/c\   NL%STGOFF = '"${STGOFF_TEST}"                     ${FILETEST}
      sed -i '/NL%STGOFF/c\   NL%STGOFF = '"${STGOFF_TEST}"                     ${FILEDBUG}
      #----- Soil boundary condition. -----------------------------------------------------#
      sed -i '/NL%ISOILBC/c\   NL%ISOILBC = '"${ISOILBC_MAIN}"                  ${FILEMAIN}
      sed -i '/NL%ISOILBC/c\   NL%ISOILBC = '"${ISOILBC_TEST}"                  ${FILETEST}
      sed -i '/NL%ISOILBC/c\   NL%ISOILBC = '"${ISOILBC_TEST}"                  ${FILEDBUG}
      #----- Integration scheme. ----------------------------------------------------------#
      sed -i '/NL%INTEGRATION_SCHEME/c\   NL%INTEGRATION_SCHEME = '"${INTEGR_MAIN}"        \
                                                                                ${FILEMAIN}
      sed -i '/NL%INTEGRATION_SCHEME/c\   NL%INTEGRATION_SCHEME = '"${INTEGR_TEST}"        \
                                                                                ${FILETEST}
      sed -i '/NL%INTEGRATION_SCHEME/c\   NL%INTEGRATION_SCHEME = '"${INTEGR_TEST}"        \
                                                                                ${FILEDBUG}
      #----- Branch themodynamics. --------------------------------------------------------#
      sed -i '/NL%IBRANCH_THERMO/c\   NL%IBRANCH_THERMO = '"${IBRANCH_MAIN}"    ${FILEMAIN}
      sed -i '/NL%IBRANCH_THERMO/c\   NL%IBRANCH_THERMO = '"${IBRANCH_TEST}"    ${FILETEST}
      sed -i '/NL%IBRANCH_THERMO/c\   NL%IBRANCH_THERMO = '"${IBRANCH_TEST}"    ${FILEDBUG}
      #----- Allometry. -------------------------------------------------------------------#
      sed -i '/NL%IALLOM/c\   NL%IALLOM = '"${IALLOM_MAIN}"                     ${FILEMAIN}
      sed -i '/NL%IALLOM/c\   NL%IALLOM = '"${IALLOM_TEST}"                     ${FILETEST}
      sed -i '/NL%IALLOM/c\   NL%IALLOM = '"${IALLOM_TEST}"                     ${FILEDBUG}
      #----- Economic. --------------------------------------------------------------------#
      sed -i '/NL%ECONOMICS_SCHEME/c\   NL%ECONOMICS_SCHEME = '"${IECON_MAIN}"  ${FILEMAIN}
      sed -i '/NL%ECONOMICS_SCHEME/c\   NL%ECONOMICS_SCHEME = '"${IECON_TEST}"  ${FILETEST}
      sed -i '/NL%ECONOMICS_SCHEME/c\   NL%ECONOMICS_SCHEME = '"${IECON_TEST}"  ${FILEDBUG}
      #----- Grass. -----------------------------------------------------------------------#
      sed -i '/NL%IGRASS/c\   NL%IGRASS = '"${IGRASS_MAIN}"                     ${FILEMAIN}
      sed -i '/NL%IGRASS/c\   NL%IGRASS = '"${IGRASS_TEST}"                     ${FILETEST}
      sed -i '/NL%IGRASS/c\   NL%IGRASS = '"${IGRASS_TEST}"                     ${FILEDBUG}
      #----- Reproduction. ----------------------------------------------------------------#
      sed -i '/NL%REPRO_SCHEME/c\   NL%REPRO_SCHEME = '"${IREPRO_MAIN}"         ${FILEMAIN}
      sed -i '/NL%REPRO_SCHEME/c\   NL%REPRO_SCHEME = '"${IREPRO_TEST}"         ${FILETEST}
      sed -i '/NL%REPRO_SCHEME/c\   NL%REPRO_SCHEME = '"${IREPRO_TEST}"         ${FILEDBUG}
      #----- Decomposition. ---------------------------------------------------------------#
      sed -i '/NL%DECOMP_SCHEME/c\   NL%DECOMP_SCHEME = '"${DECOMP_MAIN}"       ${FILEMAIN}
      sed -i '/NL%DECOMP_SCHEME/c\   NL%DECOMP_SCHEME = '"${DECOMP_TEST}"       ${FILETEST}
      sed -i '/NL%DECOMP_SCHEME/c\   NL%DECOMP_SCHEME = '"${DECOMP_TEST}"       ${FILEDBUG}
      #----- Water limitation. ------------------------------------------------------------#
      sed -i '/NL%H2O_PLANT_LIM/c\   NL%H2O_PLANT_LIM = '"${H2OLIM_MAIN}"       ${FILEMAIN}
      sed -i '/NL%H2O_PLANT_LIM/c\   NL%H2O_PLANT_LIM = '"${H2OLIM_TEST}"       ${FILETEST}
      sed -i '/NL%H2O_PLANT_LIM/c\   NL%H2O_PLANT_LIM = '"${H2OLIM_TEST}"       ${FILEDBUG}
      #----- Structural growth. -----------------------------------------------------------#
      sed -i '/NL%ISTRUCT_GROWTH_SCHEME/c\   NL%ISTRUCT_GROWTH_SCHEME = '"${STGROW_MAIN}"  \
                                                                                ${FILEMAIN}
      sed -i '/NL%ISTRUCT_GROWTH_SCHEME/c\   NL%ISTRUCT_GROWTH_SCHEME = '"${STGROW_TEST}"  \
                                                                                ${FILETEST}
      sed -i '/NL%ISTRUCT_GROWTH_SCHEME/c\   NL%ISTRUCT_GROWTH_SCHEME = '"${STGROW_TEST}"  \
                                                                                ${FILEDBUG}
      #----- Trait plasticity. ------------------------------------------------------------#
      sed -i                                                                               \
        '/NL%TRAIT_PLASTICITY_SCHEME/c\   NL%TRAIT_PLASTICITY_SCHEME = '"${PLASTIC_MAIN}"  \
                                                                                ${FILEMAIN}
      sed -i                                                                               \
        '/NL%TRAIT_PLASTICITY_SCHEME/c\   NL%TRAIT_PLASTICITY_SCHEME = '"${PLASTIC_TEST}"  \
                                                                                ${FILETEST}
      sed -i                                                                               \
        '/NL%TRAIT_PLASTICITY_SCHEME/c\   NL%TRAIT_PLASTICITY_SCHEME = '"${PLASTIC_TEST}"  \
                                                                                ${FILEDBUG}
      #----- Anthropogenic disturbance. ---------------------------------------------------#
      sed -i '/NL%INCLUDE_FIRE/c\   NL%INCLUDE_FIRE = '"${IFIRE_MAIN}"          ${FILEMAIN}
      sed -i '/NL%INCLUDE_FIRE/c\   NL%INCLUDE_FIRE = '"${IFIRE_TEST}"          ${FILETEST}
      sed -i '/NL%INCLUDE_FIRE/c\   NL%INCLUDE_FIRE = '"${IFIRE_TEST}"          ${FILEDBUG}
      #----- Fire disturbance. ------------------------------------------------------------#
      sed -i '/NL%INCLUDE_FIRE/c\   NL%INCLUDE_FIRE = '"${IFIRE_MAIN}"          ${FILEMAIN}
      sed -i '/NL%INCLUDE_FIRE/c\   NL%INCLUDE_FIRE = '"${IFIRE_TEST}"          ${FILETEST}
      sed -i '/NL%INCLUDE_FIRE/c\   NL%INCLUDE_FIRE = '"${IFIRE_TEST}"          ${FILEDBUG}
      #----- Soil moisture threshold for fire. --------------------------------------------#
      sed -i '/NL%SM_FIRE/c\   NL%SM_FIRE = '"${SMFIRE_MAIN}"                   ${FILEMAIN}
      sed -i '/NL%SM_FIRE/c\   NL%SM_FIRE = '"${SMFIRE_TEST}"                   ${FILETEST}
      sed -i '/NL%SM_FIRE/c\   NL%SM_FIRE = '"${SMFIRE_TEST}"                   ${FILEDBUG}
      #----- Phenology. -------------------------------------------------------------------#
      sed -i '/NL%IPHEN_SCHEME/c\   NL%IPHEN_SCHEME = '"${PHENOL_MAIN}"         ${FILEMAIN}
      sed -i '/NL%IPHEN_SCHEME/c\   NL%IPHEN_SCHEME = '"${PHENOL_TEST}"         ${FILETEST}
      sed -i '/NL%IPHEN_SCHEME/c\   NL%IPHEN_SCHEME = '"${PHENOL_TEST}"         ${FILEDBUG}
      #----- Soil moisture threshold for drought deciduous. -------------------------------#
      sed -i '/NL%THETACRIT/c\   NL%THETACRIT = '"${THCRIT_MAIN}"               ${FILEMAIN}
      sed -i '/NL%THETACRIT/c\   NL%THETACRIT = '"${THCRIT_TEST}"               ${FILETEST}
      sed -i '/NL%THETACRIT/c\   NL%THETACRIT = '"${THCRIT_TEST}"               ${FILEDBUG}
      #----- Canopy drag parametrisation. -------------------------------------------------#
      sed -i '/NL%ICANTURB/c\   NL%ICANTURB = '"${ICANTURB_MAIN}"               ${FILEMAIN}
      sed -i '/NL%ICANTURB/c\   NL%ICANTURB = '"${ICANTURB_TEST}"               ${FILETEST}
      sed -i '/NL%ICANTURB/c\   NL%ICANTURB = '"${ICANTURB_TEST}"               ${FILEDBUG}
      #----- Canopy turbulence parametrisation. -------------------------------------------#
      sed -i '/NL%ISFCLYRM/c\   NL%ISFCLYRM = '"${ISFCLYRM_MAIN}"               ${FILEMAIN}
      sed -i '/NL%ISFCLYRM/c\   NL%ISFCLYRM = '"${ISFCLYRM_TEST}"               ${FILETEST}
      sed -i '/NL%ISFCLYRM/c\   NL%ISFCLYRM = '"${ISFCLYRM_TEST}"               ${FILEDBUG}
      #----- Surface water percolation parametrisation. -----------------------------------#
      sed -i '/NL%IPERCOL/c\   NL%IPERCOL = '"${IPERCOL_MAIN}"                  ${FILEMAIN}
      sed -i '/NL%IPERCOL/c\   NL%IPERCOL = '"${IPERCOL_TEST}"                  ${FILETEST}
      sed -i '/NL%IPERCOL/c\   NL%IPERCOL = '"${IPERCOL_TEST}"                  ${FILEDBUG}
      #----- Tree fall disturbance rate. --------------------------------------------------#
      sed -i                                                                               \
       '/NL%TREEFALL_DISTURBANCE_RATE/c\   NL%TREEFALL_DISTURBANCE_RATE = '"${TFALL_MAIN}" \
                                                                                ${FILEMAIN}
      sed -i                                                                               \
       '/NL%TREEFALL_DISTURBANCE_RATE/c\   NL%TREEFALL_DISTURBANCE_RATE = '"${TFALL_TEST}" \
                                                                                ${FILETEST}
      sed -i                                                                               \
       '/NL%TREEFALL_DISTURBANCE_RATE/c\   NL%TREEFALL_DISTURBANCE_RATE = '"${TFALL_TEST}" \
                                                                                ${FILEDBUG}
      #----- Patch/cohort fusion method. --------------------------------------------------#
      sed -i '/NL%IFUSION/c\   NL%IFUSION = '"${IFUSION_MAIN}"                  ${FILEMAIN}
      sed -i '/NL%IFUSION/c\   NL%IFUSION = '"${IFUSION_TEST}"                  ${FILETEST}
      sed -i '/NL%IFUSION/c\   NL%IFUSION = '"${IFUSION_TEST}"                  ${FILEDBUG}
      #----- Vegetation dynamics flag. ----------------------------------------------------#
      sed -i '/NL%IVEGT_DYNAMICS/c\   NL%IVEGT_DYNAMICS = '"${IVEGTDYN_MAIN}"   ${FILEMAIN}
      sed -i '/NL%IVEGT_DYNAMICS/c\   NL%IVEGT_DYNAMICS = '"${IVEGTDYN_TEST}"   ${FILETEST}
      sed -i '/NL%IVEGT_DYNAMICS/c\   NL%IVEGT_DYNAMICS = '"${IVEGTDYN_TEST}"   ${FILEDBUG}
      #----- Error tolerance for Runge-Kutta integrator. ----------------------------------#
      sed -i '/NL%RK4_TOLERANCE/c\   NL%RK4_TOLERANCE = '"${RK4TOLER_MAIN}"     ${FILEMAIN}
      sed -i '/NL%RK4_TOLERANCE/c\   NL%RK4_TOLERANCE = '"${RK4TOLER_TEST}"     ${FILETEST}
      sed -i '/NL%RK4_TOLERANCE/c\   NL%RK4_TOLERANCE = '"${RK4TOLER_TEST}"     ${FILEDBUG}
      #----- Maximum number of sites. -----------------------------------------------------#
      sed -i '/NL%MAXSITE/c\   NL%MAXSITE = '"${MAXSITE_MAIN}"                  ${FILEMAIN}
      sed -i '/NL%MAXSITE/c\   NL%MAXSITE = '"${MAXSITE_TEST}"                  ${FILETEST}
      sed -i '/NL%MAXSITE/c\   NL%MAXSITE = '"${MAXSITE_TEST}"                  ${FILEDBUG}
      #----- Root conductance (grass). ----------------------------------------------------#
      sed -i '/NL%KW_GRASS/c\   NL%KW_GRASS = '"${KWGRASS_MAIN}"                ${FILEMAIN}
      sed -i '/NL%KW_GRASS/c\   NL%KW_GRASS = '"${KWGRASS_TEST}"                ${FILETEST}
      sed -i '/NL%KW_GRASS/c\   NL%KW_GRASS = '"${KWGRASS_TEST}"                ${FILEDBUG}
      #----- Root conductance (grass). ----------------------------------------------------#
      sed -i '/NL%KW_TREE/c\   NL%KW_TREE = '"${KWGRASS_MAIN}"                  ${FILEMAIN}
      sed -i '/NL%KW_TREE/c\   NL%KW_TREE = '"${KWGRASS_TEST}"                  ${FILETEST}
      sed -i '/NL%KW_TREE/c\   NL%KW_TREE = '"${KWGRASS_TEST}"                  ${FILEDBUG}
      #----- VPD threshold (grass). -------------------------------------------------------#
      sed -i '/NL%D0_GRASS/c\   NL%D0_GRASS = '"${D0GRASS_MAIN}"                ${FILEMAIN}
      sed -i '/NL%D0_GRASS/c\   NL%D0_GRASS = '"${D0GRASS_TEST}"                ${FILETEST}
      sed -i '/NL%D0_GRASS/c\   NL%D0_GRASS = '"${D0GRASS_TEST}"                ${FILEDBUG}
      #----- VPD threshold (grass). -------------------------------------------------------#
      sed -i '/NL%D0_TREE/c\   NL%D0_TREE = '"${D0GRASS_MAIN}"                  ${FILEMAIN}
      sed -i '/NL%D0_TREE/c\   NL%D0_TREE = '"${D0GRASS_TEST}"                  ${FILETEST}
      sed -i '/NL%D0_TREE/c\   NL%D0_TREE = '"${D0GRASS_TEST}"                  ${FILEDBUG}
      #----- Initial condition file (it may be overwritten if this is HISTORY run). -------#
      sed -i '/NL%SFILIN/c\   NL%SFILIN = '\'"${SFILIN}"\'                      ${FILEMAIN}
      sed -i '/NL%SFILIN/c\   NL%SFILIN = '\'"${SFILIN}"\'                      ${FILETEST}
      sed -i '/NL%SFILIN/c\   NL%SFILIN = '\'"${SFILIN}"\'                      ${FILEDBUG}
      #----- Meteorological driver. -------------------------------------------------------#
      sed -i '/NL%ED_MET_DRIVER_DB/c\   NL%ED_MET_DRIVER_DB = '\'"${METDRIV}"\' ${FILEMAIN}
      sed -i '/NL%ED_MET_DRIVER_DB/c\   NL%ED_MET_DRIVER_DB = '\'"${METDRIV}"\' ${FILETEST}
      sed -i '/NL%ED_MET_DRIVER_DB/c\   NL%ED_MET_DRIVER_DB = '\'"${METDRIV}"\' ${FILEDBUG}
      #----- First year of meteorological driver. -----------------------------------------#
      sed -i '/NL%METCYC1/c\   NL%METCYC1 = '"${METCYCA}"                       ${FILEMAIN}
      sed -i '/NL%METCYC1/c\   NL%METCYC1 = '"${METCYCA}"                       ${FILETEST}
      sed -i '/NL%METCYC1/c\   NL%METCYC1 = '"${METCYCA}"                       ${FILEDBUG}
      #----- Last year of meteorological driver. ------------------------------------------#
      sed -i '/NL%METCYCF/c\   NL%METCYCF = '"${METCYCZ}"                       ${FILEMAIN}
      sed -i '/NL%METCYCF/c\   NL%METCYCF = '"${METCYCZ}"                       ${FILETEST}
      sed -i '/NL%METCYCF/c\   NL%METCYCF = '"${METCYCZ}"                       ${FILEDBUG}
      #----- Time-averaging flag (meteorological driver). ---------------------------------#
      sed -i '/NL%IMETAVG/c\   NL%IMETAVG = '"${IMETAVG}"                       ${FILEMAIN}
      sed -i '/NL%IMETAVG/c\   NL%IMETAVG = '"${IMETAVG}"                       ${FILETEST}
      sed -i '/NL%IMETAVG/c\   NL%IMETAVG = '"${IMETAVG}"                       ${FILEDBUG}
      #----- Radiation decomposition scheme. ----------------------------------------------#
      sed -i '/NL%IMETRAD/c\   NL%IMETRAD = '"${IMETRAD}"                       ${FILEMAIN}
      sed -i '/NL%IMETRAD/c\   NL%IMETRAD = '"${IMETRAD}"                       ${FILETEST}
      sed -i '/NL%IMETRAD/c\   NL%IMETRAD = '"${IMETRAD}"                       ${FILEDBUG}
      #----- PFTs to include in this simulation. ------------------------------------------#
      sed -i '/NL%INCLUDE_THESE_PFT/c\   NL%INCLUDE_THESE_PFT = '"${INCPFT_MAIN}"          \
                                                                                ${FILEMAIN}
      sed -i '/NL%INCLUDE_THESE_PFT/c\   NL%INCLUDE_THESE_PFT = '"${INCPFT_TEST}"          \
                                                                                ${FILETEST}
      sed -i '/NL%INCLUDE_THESE_PFT/c\   NL%INCLUDE_THESE_PFT = '"${INCPFT_TEST}"          \
                                                                                ${FILEDBUG}
      #----- PFT for pastures. ------------------------------------------------------------#
      sed -i '/NL%PASTURE_STOCK/c\   NL%PASTURE_STOCK = '"${PASTURE_MAIN}"      ${FILEMAIN}
      sed -i '/NL%PASTURE_STOCK/c\   NL%PASTURE_STOCK = '"${PASTURE_TEST}"      ${FILETEST}
      sed -i '/NL%PASTURE_STOCK/c\   NL%PASTURE_STOCK = '"${PASTURE_TEST}"      ${FILEDBUG}
      #----- PFT for croplands. -----------------------------------------------------------#
      sed -i '/NL%AGRI_STOCK/c\   NL%AGRI_STOCK = '"${AGRI_MAIN}"               ${FILEMAIN}
      sed -i '/NL%AGRI_STOCK/c\   NL%AGRI_STOCK = '"${AGRI_TEST}"               ${FILETEST}
      sed -i '/NL%AGRI_STOCK/c\   NL%AGRI_STOCK = '"${AGRI_TEST}"               ${FILEDBUG}
      #----- PFT for forest plantation. ---------------------------------------------------#
      sed -i '/NL%PLANTATION_STOCK/c\   NL%PLANTATION_STOCK = '"${PLANT_MAIN}"  ${FILEMAIN}
      sed -i '/NL%PLANTATION_STOCK/c\   NL%PLANTATION_STOCK = '"${PLANT_TEST}"  ${FILETEST}
      sed -i '/NL%PLANTATION_STOCK/c\   NL%PLANTATION_STOCK = '"${PLANT_TEST}"  ${FILEDBUG}
      #----- PFT for selective logging. ---------------------------------------------------#
      sed -i '/NL%SL_PFT/c\   NL%SL_PFT = '"${SLPFT_MAIN}"                      ${FILEMAIN}
      sed -i '/NL%SL_PFT/c\   NL%SL_PFT = '"${SLPFT_TEST}"                      ${FILETEST}
      sed -i '/NL%SL_PFT/c\   NL%SL_PFT = '"${SLPFT_TEST}"                      ${FILEDBUG}
      #----- Harvesting probability for selective logging. --------------------------------#
      sed -i '/NL%SL_PROB_HARVEST/c\   NL%SL_PROB_HARVEST = '"${SLPROBHARV_MAIN}"          \
                                                                                ${FILEMAIN}
      sed -i '/NL%SL_PROB_HARVEST/c\   NL%SL_PROB_HARVEST = '"${SLPROBHARV_TEST}"          \
                                                                                ${FILETEST}
      sed -i '/NL%SL_PROB_HARVEST/c\   NL%SL_PROB_HARVEST = '"${SLPROBHARV_TEST}"          \
                                                                                ${FILEDBUG}
      #----- Minimum DBH for selective logging. -------------------------------------------#
      sed -i '/NL%SL_MINDBH_HARVEST/c\   NL%SL_MINDBH_HARVEST = '"${SLMINDBH_MAIN}"        \
                                                                                ${FILEMAIN}
      sed -i '/NL%SL_MINDBH_HARVEST/c\   NL%SL_MINDBH_HARVEST = '"${SLMINDBH_TEST}"        \
                                                                                ${FILETEST}
      sed -i '/NL%SL_MINDBH_HARVEST/c\   NL%SL_MINDBH_HARVEST = '"${SLMINDBH_TEST}"        \
                                                                                ${FILEDBUG}
      #----- Land use data base. ----------------------------------------------------------#
      sed -i '/NL%LU_DATABASE/c\   NL%LU_DATABASE = '\'"${LU_DATABASE}"\'       ${FILEMAIN}
      sed -i '/NL%LU_DATABASE/c\   NL%LU_DATABASE = '\'"${LU_DATABASE}"\'       ${FILETEST}
      sed -i '/NL%LU_DATABASE/c\   NL%LU_DATABASE = '\'"${LU_DATABASE}"\'       ${FILEDBUG}
      #----- Phenology data base. ---------------------------------------------------------#
      sed -i '/NL%PHENPATH/c\   NL%PHENPATH = '\'"${PHEN_PATH}"\'               ${FILEMAIN}
      sed -i '/NL%PHENPATH/c\   NL%PHENPATH = '\'"${PHEN_PATH}"\'               ${FILETEST}
      sed -i '/NL%PHENPATH/c\   NL%PHENPATH = '\'"${PHEN_PATH}"\'               ${FILEDBUG}
      #----- First year of the phenology (Spring). ----------------------------------------#
      sed -i '/NL%IPHENYS1/c\   NL%IPHENYS1 = '"${IPHENYSA}"                    ${FILEMAIN}
      sed -i '/NL%IPHENYS1/c\   NL%IPHENYS1 = '"${IPHENYSA}"                    ${FILETEST}
      sed -i '/NL%IPHENYS1/c\   NL%IPHENYS1 = '"${IPHENYSA}"                    ${FILEDBUG}
      #----- Last year of the phenology (Spring). -----------------------------------------#
      sed -i '/NL%IPHENYSF/c\   NL%IPHENYSF = '"${IPHENYSZ}"                    ${FILEMAIN}
      sed -i '/NL%IPHENYSF/c\   NL%IPHENYSF = '"${IPHENYSZ}"                    ${FILETEST}
      sed -i '/NL%IPHENYSF/c\   NL%IPHENYSF = '"${IPHENYSZ}"                    ${FILEDBUG}
      #----- First year of the phenology (Autumn). ----------------------------------------#
      sed -i '/NL%IPHENYF1/c\   NL%IPHENYF1 = '"${IPHENYFA}"                    ${FILEMAIN}
      sed -i '/NL%IPHENYF1/c\   NL%IPHENYF1 = '"${IPHENYFA}"                    ${FILETEST}
      sed -i '/NL%IPHENYF1/c\   NL%IPHENYF1 = '"${IPHENYFA}"                    ${FILEDBUG}
      #----- Last year of the phenology (Autumn). -----------------------------------------#
      sed -i '/NL%IPHENYFF/c\   NL%IPHENYFF = '"${IPHENYFZ}"                    ${FILEMAIN}
      sed -i '/NL%IPHENYFF/c\   NL%IPHENYFF = '"${IPHENYFZ}"                    ${FILETEST}
      sed -i '/NL%IPHENYFF/c\   NL%IPHENYFF = '"${IPHENYFZ}"                    ${FILEDBUG}
      #----- XML file. --------------------------------------------------------------------#
      sed -i '/NL%IEDCNFGF/c\   NL%IEDCNFGF = '\'"${EDCNFGF}"\'                 ${FILEMAIN}
      sed -i '/NL%IEDCNFGF/c\   NL%IEDCNFGF = '\'"${EDCNFGF}"\'                 ${FILETEST}
      sed -i '/NL%IEDCNFGF/c\   NL%IEDCNFGF = '\'"${EDCNFGF}"\'                 ${FILEDBUG}
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #      Detailed output settings.                                                     #
      #------------------------------------------------------------------------------------#
      #----- Detailed output list. --------------------------------------------------------#
      sed -i '/NL%IDETAILED/c\   NL%IDETAILED = '"${IDETAILEDH[i]}"            ${FILEMAIN}
      sed -i '/NL%IDETAILED/c\   NL%IDETAILED = '"${IDETAILEDH[i]}"            ${FILETEST}
      sed -i '/NL%IDETAILED/c\   NL%IDETAILED = '"${IDETAILEDH[i]}"            ${FILEDBUG}
      #----- Patches to keep. -------------------------------------------------------------#
      sed -i '/NL%PATCH_KEEP/c\   NL%PATCH_KEEP = -1'                          ${FILEMAIN}
      sed -i '/NL%PATCH_KEEP/c\   NL%PATCH_KEEP = -1'                          ${FILETEST}
      sed -i '/NL%PATCH_KEEP/c\   NL%PATCH_KEEP = -1'                          ${FILEDBUG}
      #----- Maximum number of patches. ---------------------------------------------------#
      sed -i '/NL%MAXPATCH/c\   NL%MAXPATCH = 1'                               ${FILEMAIN}
      sed -i '/NL%MAXPATCH/c\   NL%MAXPATCH = 1'                               ${FILETEST}
      sed -i '/NL%MAXPATCH/c\   NL%MAXPATCH = 1'                               ${FILEDBUG}
      #------------------------------------------------------------------------------------#


      #----- Decide whether to continue the run or start from beginning. ------------------#
      case ${RUNTYPE} in
      INITIAL)
         #----- Use default initialisation. -----------------------------------------------#
         sed -i '/NL%RUNTYPE/c\   NL%RUNTYPE = '\'${RUNTYPH[i]}\'   ${FILEMAIN}
         sed -i '/NL%RUNTYPE/c\   NL%RUNTYPE = '\'${RUNTYPH[i]}\'   ${FILETEST}
         sed -i '/NL%RUNTYPE/c\   NL%RUNTYPE = '\'${RUNTYPH[i]}\' ${FILEDBUG}
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
            sed -i '/NL%RUNTYPE/c\   NL%RUNTYPE = '\'${RUNTYPH[i]}\'   ${FILEDBUG}
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


      #----- Modify the ED2IN start times. ------------------------------------------------#
      sed -i '/NL%IYEARA/c\   NL%IYEARA   = '${IYEARAH[i]}    ${FILEMAIN}
      sed -i '/NL%IYEARA/c\   NL%IYEARA   = '${IYEARAH[i]}    ${FILETEST}
      sed -i '/NL%IYEARA/c\   NL%IYEARA   = '${D_IYEARAH[i]}  ${FILEDBUG}
      sed -i '/NL%IMONTHA/c\   NL%IMONTHA   = '${IMONTHAH[i]} ${FILEMAIN}
      sed -i '/NL%IMONTHA/c\   NL%IMONTHA   = '${IMONTHAH[i]} ${FILETEST}
      sed -i '/NL%IMONTHA/c\   NL%IMONTHA   = '${IMONTHAH[i]} ${FILEDBUG}
      sed -i '/NL%IDATEA/c\   NL%IDATEA   = '${IDATEAH[i]}    ${FILEMAIN}
      sed -i '/NL%IDATEA/c\   NL%IDATEA   = '${IDATEAH[i]}    ${FILETEST}
      sed -i '/NL%IDATEA/c\   NL%IDATEA   = '${IDATEAH[i]}    ${FILEDBUG}
      sed -i '/NL%ITIMEA/c\   NL%ITIMEA   = 0000'             ${FILEMAIN}
      sed -i '/NL%ITIMEA/c\   NL%ITIMEA   = 0000'             ${FILETEST}
      sed -i '/NL%ITIMEA/c\   NL%ITIMEA   = 0000'             ${FILEDBUG}
      #------------------------------------------------------------------------------------#


      #----- Modify the ED2IN end years. --------------------------------------------------#
      sed -i '/NL%IYEARZ/c\   NL%IYEARZ   = '${IYEARZH[i]}    ${FILEMAIN}
      sed -i '/NL%IYEARZ/c\   NL%IYEARZ   = '${IYEARZH[i]}    ${FILETEST}
      sed -i '/NL%IYEARZ/c\   NL%IYEARZ   = '${D_IYEARZH[i]}  ${FILEDBUG}
      sed -i '/NL%IMONTHZ/c\   NL%IMONTHZ   = '${IMONTHZH[i]} ${FILEMAIN}
      sed -i '/NL%IMONTHZ/c\   NL%IMONTHZ   = '${IMONTHZH[i]} ${FILETEST}
      sed -i '/NL%IMONTHZ/c\   NL%IMONTHZ   = '${IMONTHZH[i]} ${FILEDBUG}
      sed -i '/NL%IDATEZ/c\   NL%IDATEZ   = '${IDATEZH[i]}    ${FILEMAIN}
      sed -i '/NL%IDATEZ/c\   NL%IDATEZ   = '${IDATEZH[i]}    ${FILETEST}
      sed -i '/NL%IDATEZ/c\   NL%IDATEZ   = '${IDATEZH[i]}    ${FILEDBUG}
      sed -i '/NL%ITIMEZ/c\   NL%ITIMEZ   = 0000'             ${FILEMAIN}
      sed -i '/NL%ITIMEZ/c\   NL%ITIMEZ   = 0000'             ${FILETEST}
      sed -i '/NL%ITIMEZ/c\   NL%ITIMEZ   = 0000'             ${FILEDBUG}
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
         jobwrap=". ${HOME}/.bashrc ${RC_OPTS}; cd ${HERE}/${VERSION}"
         jobwrap="${jobwrap}; export OMP_NUM_THREADS=${HIFRCPU[i]}"
         jobwrap="${jobwrap}; ${LNK_MAIN_EXE} -f ${FILEMAIN}"
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
         jobwrap=". ${HOME}/.bashrc ${RC_OPTS}; cd ${HERE}/${VERSION}"
         jobwrap="${jobwrap}; export OMP_NUM_THREADS=${HIFRCPU[i]}"
         jobwrap="${jobwrap}; ${LNK_TEST_EXE} -f ${FILETEST}"
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
         jobwrap=". ${HOME}/.bashrc ${RC_OPTS}; cd ${HERE}/${VERSION}"
         jobwrap="${jobwrap}; export OMP_NUM_THREADS=${HIFRCPU[i]}"
         jobwrap="${jobwrap}; ${LNK_DBUG_EXE} -f ${FILEDBUG}"
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
      TEMPLATEMAIN="${HERE}/Templates/ED2IN-MAIN"
      TEMPLATETEST="${HERE}/Templates/ED2IN-TEST"
      TEMPLATEDBUG="${TEMPLATETEST}"
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



      #------------------------------------------------------------------------------------#
      #     Extract most variable settings.                                                #
      #------------------------------------------------------------------------------------#
      nml_setting ${GRIDRESG[i]} ; GRIDRES_MAIN=${nl_main} ; GRIDRES_TEST=${nl_test}
      nml_setting ${REGWLONG[i]} ; REGWLON_MAIN=${nl_main} ; REGWLON_TEST=${nl_test}
      nml_setting ${REGELONG[i]} ; REGELON_MAIN=${nl_main} ; REGELON_TEST=${nl_test}
      nml_setting ${REGSLATG[i]} ; REGSLAT_MAIN=${nl_main} ; REGSLAT_TEST=${nl_test}
      nml_setting ${REGNLATG[i]} ; REGNLAT_MAIN=${nl_main} ; REGNLAT_TEST=${nl_test}
      nml_setting ${ISOILFLGG[i]}; ISOILFLG_MAIN=${nl_main}; ISOILFLG_TEST=${nl_test}
      nml_setting ${NSLCONG[i]}  ; NSLCON_MAIN=${nl_main}  ; NSLCON_TEST=${nl_test}
      nml_setting ${SLXCLAYG[i]} ; SLXCLAY_MAIN=${nl_main} ; SLXCLAY_TEST=${nl_test}
      nml_setting ${SLXSANDG[i]} ; SLXSAND_MAIN=${nl_main} ; SLXSAND_TEST=${nl_test}
      nml_setting ${ISOILCOLG[i]}; ISOILCOL_MAIN=${nl_main}; ISOILCOL_TEST=${nl_test}
      nml_setting ${SDEPTHG[i]}  ; SDEPTH_MAIN=${nl_main}  ; SDEPTH_TEST=${nl_test}
      nml_setting ${ISOILBCG[i]} ; ISOILBC_MAIN=${nl_main} ; ISOILBC_TEST=${nl_test}
      nml_setting ${INTEGRG[i]}  ; INTEGR_MAIN=${nl_main}  ; INTEGR_TEST=${nl_test}
      nml_setting ${IBRANCHG[i]} ; IBRANCH_MAIN=${nl_main} ; IBRANCH_TEST=${nl_test}
      nml_setting ${IALLOMG[i]}  ; IALLOM_MAIN=${nl_main}  ; IALLOM_TEST=${nl_test}
      nml_setting ${IECONG[i]}   ; IECON_MAIN=${nl_main}   ; IECON_TEST=${nl_test}
      nml_setting ${IGRASSG[i]}  ; IGRASS_MAIN=${nl_main}  ; IGRASS_TEST=${nl_test}
      nml_setting ${IREPROG[i]}  ; IREPRO_MAIN=${nl_main}  ; IREPRO_TEST=${nl_test}
      nml_setting ${DECOMPG[i]}  ; DECOMP_MAIN=${nl_main}  ; DECOMP_TEST=${nl_test}
      nml_setting ${H2OLIMG[i]}  ; H2OLIM_MAIN=${nl_main}  ; H2OLIM_TEST=${nl_test}
      nml_setting ${STGROWG[i]}  ; STGROW_MAIN=${nl_main}  ; STGROW_TEST=${nl_test}
      nml_setting ${PLASTICG[i]} ; PLASTIC_MAIN=${nl_main} ; PLASTIC_TEST=${nl_test}
      nml_setting ${IANTHG[i]}   ; IANTH_MAIN=${nl_main}   ; IANTH_TEST=${nl_test}
      nml_setting ${IFIREG[i]}   ; IFIRE_MAIN=${nl_main}   ; IFIRE_TEST=${nl_test}
      nml_setting ${SMFIREG[i]}  ; SMFIRE_MAIN=${nl_main}  ; SMFIRE_TEST=${nl_test}
      nml_setting ${PHENOLG[i]}  ; PHENOL_MAIN=${nl_main}  ; PHENOL_TEST=${nl_test}
      nml_setting ${THCRITG[i]}  ; THCRIT_MAIN=${nl_main}  ; THCRIT_TEST=${nl_test}
      nml_setting ${ICANTURBG[i]}; ICANTURB_MAIN=${nl_main}; ICANTURB_TEST=${nl_test}
      nml_setting ${ISFCLYRMG[i]}; ISFCLYRM_MAIN=${nl_main}; ISFCLYRM_TEST=${nl_test}
      nml_setting ${IPERCOLG[i]} ; IPERCOL_MAIN=${nl_main} ; IPERCOL_TEST=${nl_test}
      nml_setting ${PFTFLAGG[i]} ; PFTFLAG_MAIN=${nl_main} ; PFTFLAG_TEST=${nl_test}
      nml_setting ${TREEFALLG[i]}; TFALL_MAIN=${nl_main}   ; TFALL_TEST=${nl_test}
      nml_setting ${IFUSIONG[i]} ; IFUSION_MAIN=${nl_main} ; IFUSION_TEST=${nl_test}
      nml_setting ${IVEGTDYNG[i]}; IVEGTDYN_MAIN=${nl_main}; IVEGTDYN_TEST=${nl_test}
      nml_setting ${RK4TOLERG[i]}; RK4TOLER_MAIN=${nl_main}; RK4TOLER_TEST=${nl_test}
      nml_setting ${MAXSITEG[i]} ; MAXSITE_MAIN=${nl_main} ; MAXSITE_TEST=${nl_test}
      #----- Derived parameters. ----------------------------------------------------------#
      kwd0_setting ${H2OLIM_MAIN}
      KWGRASS_MAIN=${kw_grass}; KWTREE_MAIN=${kw_tree}
      D0GRASS_MAIN=${d0_grass}; D0TREE_MAIN=${d0_tree}
      kwd0_setting ${H2OLIM_TEST}
      KWGRASS_TEST=${kw_grass}; KWTREE_MAIN=${kw_tree}
      D0GRASS_TEST=${d0_grass}; D0TREE_MAIN=${d0_tree}
      #----- Retrieve PFT settings. -------------------------------------------------------#
      pft_setting ${PFTFLAG_MAIN}
      INCPFT_MAIN=${inc_pft}  ; PASTURE_MAIN=${pasture}       ; AGRI_MAIN=${agri}
      PLANT_MAIN=${plantation}
      SLPFT_MAIN=${sl_pft}    ; SLPROBHARV_MAIN=${sl_probharv}; SLMINDBH_MAIN=${sl_mindbh}
      pft_setting ${PFTFLAG_TEST}
      INCPFT_TEST=${inc_pft}; PASTURE_TEST=${pasture}       ; AGRI_TEST=${agri}
      PLANT_TEST=${plantation}
      SLPFT_TEST=${sl_pft}  ; SLPROBHARV_TEST=${sl_probharv}; SLMINDBH_TEST=${sl_mindbh}
      #----- Obtain soil layers. ----------------------------------------------------------#
      soil_setting ${SDEPTH_MAIN}
      NZG_MAIN=${nzg}; SLZ_MAIN=${slz}; SLMSTR_MAIN=${slm}; STGOFF_MAIN=${slt}
      soil_setting ${SDEPTH_TEST}
      NZG_TEST=${nzg}; SLZ_TEST=${slz}; SLMSTR_TEST=${slm}; STGOFF_TEST=${slt}
      #----- Define inital file. ----------------------------------------------------------#
      sfilin_setting ${GRIDID[i]} ${INITMDG[i]}; SFILIN=${init_file}
      #----- Define meteorological driver settings. ---------------------------------------#
      metdriv_setting ${GRIDID[i]}
      METDRIV=${metdriv}; METCYCA=${metcyca}; METCYCZ=${metcycz}
      IMETAVG=${imetavg}; IMETRAD=${imetrad}
      #----- Define land use driver settings. ---------------------------------------------#
      lu_setting ${GRIDID[i]}; LU_DATABASE=${lu_database}
      #----- Define prescribed-phenology settings. ----------------------------------------#
      phenol_setting ${GRIDID[i]}
      PHEN_PATH=${phen_path}
      IPHENYSA=${iphenysa}
      IPHENYSZ=${iphenysz}
      IPHENYFA=${iphenyfa}
      IPHENYFZ=${iphenyfz}
      #----- Define parameter file settings. ----------------------------------------------#
      xml_setting ${GRIDID[i]}
      EDCNFGF=${edcnfgf}
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Extract most variable settings.                                                #
      #------------------------------------------------------------------------------------#
      #----- Number of ED regions. --------------------------------------------------------#
      sed -i '/NL%N_ED_REGION/c\   NL%N_ED_REGION = 1'                          ${FILEMAIN}
      sed -i '/NL%N_ED_REGION/c\   NL%N_ED_REGION = 1'                          ${FILETEST}
      sed -i '/NL%N_ED_REGION/c\   NL%N_ED_REGION = 1'                          ${FILEDBUG}
      #----- Grid type (for now we only test regular lon/lat grid). -----------------------#
      sed -i '/NL%GRID_TYPE/c\   NL%GRID_TYPE = 0'                              ${FILEMAIN}
      sed -i '/NL%GRID_TYPE/c\   NL%GRID_TYPE = 0'                              ${FILETEST}
      sed -i '/NL%GRID_TYPE/c\   NL%GRID_TYPE = 0'                              ${FILEDBUG}
      #----- Number of ED Polygons. -------------------------------------------------------#
      sed -i '/NL%N_POI/c\   NL%N_POI = 0'                                      ${FILEMAIN}
      sed -i '/NL%N_POI/c\   NL%N_POI = 0'                                      ${FILETEST}
      sed -i '/NL%N_POI/c\   NL%N_POI = 0'                                      ${FILEDBUG}
      #----- Grid resolution. -------------------------------------------------------------#
      sed -i '/NL%GRID_RES/c\   NL%GRID_RES = '"${GRIDRES_MAIN}"                ${FILEMAIN}
      sed -i '/NL%GRID_RES/c\   NL%GRID_RES = '"${GRIDRES_TEST}"                ${FILETEST}
      sed -i '/NL%GRID_RES/c\   NL%GRID_RES = '"${GRIDRES_TEST}"                ${FILEDBUG}
      #----- Westernmost longitude. -------------------------------------------------------#
      sed -i '/NL%ED_REG_LONMIN/c\   NL%ED_REG_LONMIN = '"${REGWLON_MAIN}"      ${FILEMAIN}
      sed -i '/NL%ED_REG_LONMIN/c\   NL%ED_REG_LONMIN = '"${REGWLON_TEST}"      ${FILETEST}
      sed -i '/NL%ED_REG_LONMIN/c\   NL%ED_REG_LONMIN = '"${REGWLON_TEST}"      ${FILEDBUG}
      #----- Easternmost longitude. -------------------------------------------------------#
      sed -i '/NL%ED_REG_LONMAX/c\   NL%ED_REG_LONMAX = '"${REGELON_MAIN}"      ${FILEMAIN}
      sed -i '/NL%ED_REG_LONMAX/c\   NL%ED_REG_LONMAX = '"${REGELON_TEST}"      ${FILETEST}
      sed -i '/NL%ED_REG_LONMAX/c\   NL%ED_REG_LONMAX = '"${REGELON_TEST}"      ${FILEDBUG}
      #----- Southernmost latitude. -------------------------------------------------------#
      sed -i '/NL%ED_REG_LATMIN/c\   NL%ED_REG_LATMIN = '"${REGSLAT_MAIN}"      ${FILEMAIN}
      sed -i '/NL%ED_REG_LATMIN/c\   NL%ED_REG_LATMIN = '"${REGSLAT_TEST}"      ${FILETEST}
      sed -i '/NL%ED_REG_LATMIN/c\   NL%ED_REG_LATMIN = '"${REGSLAT_TEST}"      ${FILEDBUG}
      #----- Northernmost latitude. -------------------------------------------------------#
      sed -i '/NL%ED_REG_LATMAX/c\   NL%ED_REG_LATMAX = '"${REGNLAT_MAIN}"      ${FILEMAIN}
      sed -i '/NL%ED_REG_LATMAX/c\   NL%ED_REG_LATMAX = '"${REGNLAT_TEST}"      ${FILETEST}
      sed -i '/NL%ED_REG_LATMAX/c\   NL%ED_REG_LATMAX = '"${REGNLAT_TEST}"      ${FILEDBUG}
      #----- Longitude. -------------------------------------------------------------------#
      sed -i '/NL%POI_LON/c\   NL%POI_LON = '"${PLON_MAIN}"                     ${FILEMAIN}
      sed -i '/NL%POI_LON/c\   NL%POI_LON = '"${PLON_TEST}"                     ${FILETEST}
      sed -i '/NL%POI_LON/c\   NL%POI_LON = '"${PLON_TEST}"                     ${FILEDBUG}
      #----- Latitude. --------------------------------------------------------------------#
      sed -i '/NL%POI_LAT/c\   NL%POI_LAT = '"${PLAT_MAIN}"                     ${FILEMAIN}
      sed -i '/NL%POI_LAT/c\   NL%POI_LAT = '"${PLAT_TEST}"                     ${FILETEST}
      sed -i '/NL%POI_LAT/c\   NL%POI_LAT = '"${PLAT_TEST}"                     ${FILEDBUG}
      #----- Soil texture flag. -----------------------------------------------------------#
      sed -i '/NL%ISOILFLG/c\   NL%ISOILFLG = '"${ISOILFLG_MAIN}"               ${FILEMAIN}
      sed -i '/NL%ISOILFLG/c\   NL%ISOILFLG = '"${ISOILFLG_TEST}"               ${FILETEST}
      sed -i '/NL%ISOILFLG/c\   NL%ISOILFLG = '"${ISOILFLG_TEST}"               ${FILEDBUG}
      #----- Default soil type. -----------------------------------------------------------#
      sed -i '/NL%NSLCON/c\   NL%NSLCON = '"${NSLCON_MAIN}"                     ${FILEMAIN}
      sed -i '/NL%NSLCON/c\   NL%NSLCON = '"${NSLCON_TEST}"                     ${FILETEST}
      sed -i '/NL%NSLCON/c\   NL%NSLCON = '"${NSLCON_TEST}"                     ${FILEDBUG}
      #----- Clay fraction. ---------------------------------------------------------------#
      sed -i '/NL%SLXCLAY/c\   NL%SLXCLAY = '"${SLXCLAY_MAIN}"                  ${FILEMAIN}
      sed -i '/NL%SLXCLAY/c\   NL%SLXCLAY = '"${SLXCLAY_TEST}"                  ${FILETEST}
      sed -i '/NL%SLXCLAY/c\   NL%SLXCLAY = '"${SLXCLAY_TEST}"                  ${FILEDBUG}
      #----- Sand fraction. ---------------------------------------------------------------#
      sed -i '/NL%SLXSAND/c\   NL%SLXSAND = '"${SLXSAND_MAIN}"                  ${FILEMAIN}
      sed -i '/NL%SLXSAND/c\   NL%SLXSAND = '"${SLXSAND_TEST}"                  ${FILETEST}
      sed -i '/NL%SLXSAND/c\   NL%SLXSAND = '"${SLXSAND_TEST}"                  ${FILEDBUG}
      #----- Default soil colour. ---------------------------------------------------------#
      sed -i '/NL%ISOILCOL/c\   NL%ISOILCOL = '"${ISOILCOL_MAIN}"               ${FILEMAIN}
      sed -i '/NL%ISOILCOL/c\   NL%ISOILCOL = '"${ISOILCOL_TEST}"               ${FILETEST}
      sed -i '/NL%ISOILCOL/c\   NL%ISOILCOL = '"${ISOILCOL_TEST}"               ${FILEDBUG}
      #----- Number of soil layers. -------------------------------------------------------#
      sed -i '/NL%NZG/c\   NL%NZG = '"${NZG_MAIN}"                              ${FILEMAIN}
      sed -i '/NL%NZG/c\   NL%NZG = '"${NZG_TEST}"                              ${FILETEST}
      sed -i '/NL%NZG/c\   NL%NZG = '"${NZG_TEST}"                              ${FILEDBUG}
      #----- Soil layers. -----------------------------------------------------------------#
      sed -i '/NL%SLZ/c\   NL%SLZ = '"${SLZ_MAIN}"                              ${FILEMAIN}
      sed -i '/NL%SLZ/c\   NL%SLZ = '"${SLZ_TEST}"                              ${FILETEST}
      sed -i '/NL%SLZ/c\   NL%SLZ = '"${SLZ_TEST}"                              ${FILEDBUG}
      #----- Soil moisture. ---------------------------------------------------------------#
      sed -i '/NL%SLMSTR/c\   NL%SLMSTR = '"${SLMSTR_MAIN}"                     ${FILEMAIN}
      sed -i '/NL%SLMSTR/c\   NL%SLMSTR = '"${SLMSTR_TEST}"                     ${FILETEST}
      sed -i '/NL%SLMSTR/c\   NL%SLMSTR = '"${SLMSTR_TEST}"                     ${FILEDBUG}
      #----- Soil temperature offset. -----------------------------------------------------#
      sed -i '/NL%STGOFF/c\   NL%STGOFF = '"${STGOFF_MAIN}"                     ${FILEMAIN}
      sed -i '/NL%STGOFF/c\   NL%STGOFF = '"${STGOFF_TEST}"                     ${FILETEST}
      sed -i '/NL%STGOFF/c\   NL%STGOFF = '"${STGOFF_TEST}"                     ${FILEDBUG}
      #----- Soil boundary condition. -----------------------------------------------------#
      sed -i '/NL%ISOILBC/c\   NL%ISOILBC = '"${ISOILBC_MAIN}"                  ${FILEMAIN}
      sed -i '/NL%ISOILBC/c\   NL%ISOILBC = '"${ISOILBC_TEST}"                  ${FILETEST}
      sed -i '/NL%ISOILBC/c\   NL%ISOILBC = '"${ISOILBC_TEST}"                  ${FILEDBUG}
      #----- Integration scheme. ----------------------------------------------------------#
      sed -i '/NL%INTEGRATION_SCHEME/c\   NL%INTEGRATION_SCHEME = '"${INTEGR_MAIN}"        \
                                                                                ${FILEMAIN}
      sed -i '/NL%INTEGRATION_SCHEME/c\   NL%INTEGRATION_SCHEME = '"${INTEGR_TEST}"        \
                                                                                ${FILETEST}
      sed -i '/NL%INTEGRATION_SCHEME/c\   NL%INTEGRATION_SCHEME = '"${INTEGR_TEST}"        \
                                                                                ${FILEDBUG}
      #----- Branch themodynamics. --------------------------------------------------------#
      sed -i '/NL%IBRANCH_THERMO/c\   NL%IBRANCH_THERMO = '"${IBRANCH_MAIN}"    ${FILEMAIN}
      sed -i '/NL%IBRANCH_THERMO/c\   NL%IBRANCH_THERMO = '"${IBRANCH_TEST}"    ${FILETEST}
      sed -i '/NL%IBRANCH_THERMO/c\   NL%IBRANCH_THERMO = '"${IBRANCH_TEST}"    ${FILEDBUG}
      #----- Allometry. -------------------------------------------------------------------#
      sed -i '/NL%IALLOM/c\   NL%IALLOM = '"${IALLOM_MAIN}"                     ${FILEMAIN}
      sed -i '/NL%IALLOM/c\   NL%IALLOM = '"${IALLOM_TEST}"                     ${FILETEST}
      sed -i '/NL%IALLOM/c\   NL%IALLOM = '"${IALLOM_TEST}"                     ${FILEDBUG}
      #----- Economic. --------------------------------------------------------------------#
      sed -i '/NL%ECONOMICS_SCHEME/c\   NL%ECONOMICS_SCHEME = '"${IECON_MAIN}"  ${FILEMAIN}
      sed -i '/NL%ECONOMICS_SCHEME/c\   NL%ECONOMICS_SCHEME = '"${IECON_TEST}"  ${FILETEST}
      sed -i '/NL%ECONOMICS_SCHEME/c\   NL%ECONOMICS_SCHEME = '"${IECON_TEST}"  ${FILEDBUG}
      #----- Grass. -----------------------------------------------------------------------#
      sed -i '/NL%IGRASS/c\   NL%IGRASS = '"${IGRASS_MAIN}"                     ${FILEMAIN}
      sed -i '/NL%IGRASS/c\   NL%IGRASS = '"${IGRASS_TEST}"                     ${FILETEST}
      sed -i '/NL%IGRASS/c\   NL%IGRASS = '"${IGRASS_TEST}"                     ${FILEDBUG}
      #----- Reproduction. ----------------------------------------------------------------#
      sed -i '/NL%REPRO_SCHEME/c\   NL%REPRO_SCHEME = '"${IREPRO_MAIN}"         ${FILEMAIN}
      sed -i '/NL%REPRO_SCHEME/c\   NL%REPRO_SCHEME = '"${IREPRO_TEST}"         ${FILETEST}
      sed -i '/NL%REPRO_SCHEME/c\   NL%REPRO_SCHEME = '"${IREPRO_TEST}"         ${FILEDBUG}
      #----- Decomposition. ---------------------------------------------------------------#
      sed -i '/NL%DECOMP_SCHEME/c\   NL%DECOMP_SCHEME = '"${DECOMP_MAIN}"       ${FILEMAIN}
      sed -i '/NL%DECOMP_SCHEME/c\   NL%DECOMP_SCHEME = '"${DECOMP_TEST}"       ${FILETEST}
      sed -i '/NL%DECOMP_SCHEME/c\   NL%DECOMP_SCHEME = '"${DECOMP_TEST}"       ${FILEDBUG}
      #----- Water limitation. ------------------------------------------------------------#
      sed -i '/NL%H2O_PLANT_LIM/c\   NL%H2O_PLANT_LIM = '"${H2OLIM_MAIN}"       ${FILEMAIN}
      sed -i '/NL%H2O_PLANT_LIM/c\   NL%H2O_PLANT_LIM = '"${H2OLIM_TEST}"       ${FILETEST}
      sed -i '/NL%H2O_PLANT_LIM/c\   NL%H2O_PLANT_LIM = '"${H2OLIM_TEST}"       ${FILEDBUG}
      #----- Structural growth. -----------------------------------------------------------#
      sed -i '/NL%ISTRUCT_GROWTH_SCHEME/c\   NL%ISTRUCT_GROWTH_SCHEME = '"${STGROW_MAIN}"  \
                                                                                ${FILEMAIN}
      sed -i '/NL%ISTRUCT_GROWTH_SCHEME/c\   NL%ISTRUCT_GROWTH_SCHEME = '"${STGROW_TEST}"  \
                                                                                ${FILETEST}
      sed -i '/NL%ISTRUCT_GROWTH_SCHEME/c\   NL%ISTRUCT_GROWTH_SCHEME = '"${STGROW_TEST}"  \
                                                                                ${FILEDBUG}
      #----- Trait plasticity. ------------------------------------------------------------#
      sed -i                                                                               \
        '/NL%TRAIT_PLASTICITY_SCHEME/c\   NL%TRAIT_PLASTICITY_SCHEME = '"${PLASTIC_MAIN}"  \
                                                                                ${FILEMAIN}
      sed -i                                                                               \
        '/NL%TRAIT_PLASTICITY_SCHEME/c\   NL%TRAIT_PLASTICITY_SCHEME = '"${PLASTIC_TEST}"  \
                                                                                ${FILETEST}
      sed -i                                                                               \
        '/NL%TRAIT_PLASTICITY_SCHEME/c\   NL%TRAIT_PLASTICITY_SCHEME = '"${PLASTIC_TEST}"  \
                                                                                ${FILEDBUG}
      #----- Anthropogenic disturbance. ---------------------------------------------------#
      sed -i '/NL%INCLUDE_FIRE/c\   NL%INCLUDE_FIRE = '"${IFIRE_MAIN}"          ${FILEMAIN}
      sed -i '/NL%INCLUDE_FIRE/c\   NL%INCLUDE_FIRE = '"${IFIRE_TEST}"          ${FILETEST}
      sed -i '/NL%INCLUDE_FIRE/c\   NL%INCLUDE_FIRE = '"${IFIRE_TEST}"          ${FILEDBUG}
      #----- Fire disturbance. ------------------------------------------------------------#
      sed -i '/NL%INCLUDE_FIRE/c\   NL%INCLUDE_FIRE = '"${IFIRE_MAIN}"          ${FILEMAIN}
      sed -i '/NL%INCLUDE_FIRE/c\   NL%INCLUDE_FIRE = '"${IFIRE_TEST}"          ${FILETEST}
      sed -i '/NL%INCLUDE_FIRE/c\   NL%INCLUDE_FIRE = '"${IFIRE_TEST}"          ${FILEDBUG}
      #----- Soil moisture threshold for fire. --------------------------------------------#
      sed -i '/NL%SM_FIRE/c\   NL%SM_FIRE = '"${SMFIRE_MAIN}"                   ${FILEMAIN}
      sed -i '/NL%SM_FIRE/c\   NL%SM_FIRE = '"${SMFIRE_TEST}"                   ${FILETEST}
      sed -i '/NL%SM_FIRE/c\   NL%SM_FIRE = '"${SMFIRE_TEST}"                   ${FILEDBUG}
      #----- Phenology. -------------------------------------------------------------------#
      sed -i '/NL%IPHEN_SCHEME/c\   NL%IPHEN_SCHEME = '"${PHENOL_MAIN}"         ${FILEMAIN}
      sed -i '/NL%IPHEN_SCHEME/c\   NL%IPHEN_SCHEME = '"${PHENOL_TEST}"         ${FILETEST}
      sed -i '/NL%IPHEN_SCHEME/c\   NL%IPHEN_SCHEME = '"${PHENOL_TEST}"         ${FILEDBUG}
      #----- Soil moisture threshold for drought deciduous. -------------------------------#
      sed -i '/NL%THETACRIT/c\   NL%THETACRIT = '"${THCRIT_MAIN}"               ${FILEMAIN}
      sed -i '/NL%THETACRIT/c\   NL%THETACRIT = '"${THCRIT_TEST}"               ${FILETEST}
      sed -i '/NL%THETACRIT/c\   NL%THETACRIT = '"${THCRIT_TEST}"               ${FILEDBUG}
      #----- Canopy drag parametrisation. -------------------------------------------------#
      sed -i '/NL%ICANTURB/c\   NL%ICANTURB = '"${ICANTURB_MAIN}"               ${FILEMAIN}
      sed -i '/NL%ICANTURB/c\   NL%ICANTURB = '"${ICANTURB_TEST}"               ${FILETEST}
      sed -i '/NL%ICANTURB/c\   NL%ICANTURB = '"${ICANTURB_TEST}"               ${FILEDBUG}
      #----- Canopy turbulence parametrisation. -------------------------------------------#
      sed -i '/NL%ISFCLYRM/c\   NL%ISFCLYRM = '"${ISFCLYRM_MAIN}"               ${FILEMAIN}
      sed -i '/NL%ISFCLYRM/c\   NL%ISFCLYRM = '"${ISFCLYRM_TEST}"               ${FILETEST}
      sed -i '/NL%ISFCLYRM/c\   NL%ISFCLYRM = '"${ISFCLYRM_TEST}"               ${FILEDBUG}
      #----- Surface water percolation parametrisation. -----------------------------------#
      sed -i '/NL%IPERCOL/c\   NL%IPERCOL = '"${IPERCOL_MAIN}"                  ${FILEMAIN}
      sed -i '/NL%IPERCOL/c\   NL%IPERCOL = '"${IPERCOL_TEST}"                  ${FILETEST}
      sed -i '/NL%IPERCOL/c\   NL%IPERCOL = '"${IPERCOL_TEST}"                  ${FILEDBUG}
      #----- Tree fall disturbance rate. --------------------------------------------------#
      sed -i                                                                               \
       '/NL%TREEFALL_DISTURBANCE_RATE/c\   NL%TREEFALL_DISTURBANCE_RATE = '"${TFALL_MAIN}" \
                                                                                ${FILEMAIN}
      sed -i                                                                               \
       '/NL%TREEFALL_DISTURBANCE_RATE/c\   NL%TREEFALL_DISTURBANCE_RATE = '"${TFALL_TEST}" \
                                                                                ${FILETEST}
      sed -i                                                                               \
       '/NL%TREEFALL_DISTURBANCE_RATE/c\   NL%TREEFALL_DISTURBANCE_RATE = '"${TFALL_TEST}" \
                                                                                ${FILEDBUG}
      #----- Patch/cohort fusion method. --------------------------------------------------#
      sed -i '/NL%IFUSION/c\   NL%IFUSION = '"${IFUSION_MAIN}"                  ${FILEMAIN}
      sed -i '/NL%IFUSION/c\   NL%IFUSION = '"${IFUSION_TEST}"                  ${FILETEST}
      sed -i '/NL%IFUSION/c\   NL%IFUSION = '"${IFUSION_TEST}"                  ${FILEDBUG}
      #----- Vegetation dynamics flag. ----------------------------------------------------#
      sed -i '/NL%IVEGT_DYNAMICS/c\   NL%IVEGT_DYNAMICS = '"${IVEGTDYN_MAIN}"   ${FILEMAIN}
      sed -i '/NL%IVEGT_DYNAMICS/c\   NL%IVEGT_DYNAMICS = '"${IVEGTDYN_TEST}"   ${FILETEST}
      sed -i '/NL%IVEGT_DYNAMICS/c\   NL%IVEGT_DYNAMICS = '"${IVEGTDYN_TEST}"   ${FILEDBUG}
      #----- Error tolerance for Runge-Kutta integrator. ----------------------------------#
      sed -i '/NL%RK4_TOLERANCE/c\   NL%RK4_TOLERANCE = '"${RK4TOLER_MAIN}"     ${FILEMAIN}
      sed -i '/NL%RK4_TOLERANCE/c\   NL%RK4_TOLERANCE = '"${RK4TOLER_TEST}"     ${FILETEST}
      sed -i '/NL%RK4_TOLERANCE/c\   NL%RK4_TOLERANCE = '"${RK4TOLER_TEST}"     ${FILEDBUG}
      #----- Maximum number of sites. -----------------------------------------------------#
      sed -i '/NL%MAXSITE/c\   NL%MAXSITE = '"${MAXSITE_MAIN}"                  ${FILEMAIN}
      sed -i '/NL%MAXSITE/c\   NL%MAXSITE = '"${MAXSITE_TEST}"                  ${FILETEST}
      sed -i '/NL%MAXSITE/c\   NL%MAXSITE = '"${MAXSITE_TEST}"                  ${FILEDBUG}
      #----- Root conductance (grass). ----------------------------------------------------#
      sed -i '/NL%KW_GRASS/c\   NL%KW_GRASS = '"${KWGRASS_MAIN}"                ${FILEMAIN}
      sed -i '/NL%KW_GRASS/c\   NL%KW_GRASS = '"${KWGRASS_TEST}"                ${FILETEST}
      sed -i '/NL%KW_GRASS/c\   NL%KW_GRASS = '"${KWGRASS_TEST}"                ${FILEDBUG}
      #----- Root conductance (grass). ----------------------------------------------------#
      sed -i '/NL%KW_TREE/c\   NL%KW_TREE = '"${KWGRASS_MAIN}"                  ${FILEMAIN}
      sed -i '/NL%KW_TREE/c\   NL%KW_TREE = '"${KWGRASS_TEST}"                  ${FILETEST}
      sed -i '/NL%KW_TREE/c\   NL%KW_TREE = '"${KWGRASS_TEST}"                  ${FILEDBUG}
      #----- VPD threshold (grass). -------------------------------------------------------#
      sed -i '/NL%D0_GRASS/c\   NL%D0_GRASS = '"${D0GRASS_MAIN}"                ${FILEMAIN}
      sed -i '/NL%D0_GRASS/c\   NL%D0_GRASS = '"${D0GRASS_TEST}"                ${FILETEST}
      sed -i '/NL%D0_GRASS/c\   NL%D0_GRASS = '"${D0GRASS_TEST}"                ${FILEDBUG}
      #----- VPD threshold (grass). -------------------------------------------------------#
      sed -i '/NL%D0_TREE/c\   NL%D0_TREE = '"${D0GRASS_MAIN}"                  ${FILEMAIN}
      sed -i '/NL%D0_TREE/c\   NL%D0_TREE = '"${D0GRASS_TEST}"                  ${FILETEST}
      sed -i '/NL%D0_TREE/c\   NL%D0_TREE = '"${D0GRASS_TEST}"                  ${FILEDBUG}
      #----- Initial condition file (it may be overwritten if this is HISTORY run). -------#
      sed -i '/NL%SFILIN/c\   NL%SFILIN = '\'"${SFILIN}"\'                      ${FILEMAIN}
      sed -i '/NL%SFILIN/c\   NL%SFILIN = '\'"${SFILIN}"\'                      ${FILETEST}
      sed -i '/NL%SFILIN/c\   NL%SFILIN = '\'"${SFILIN}"\'                      ${FILEDBUG}
      #----- Meteorological driver. -------------------------------------------------------#
      sed -i '/NL%ED_MET_DRIVER_DB/c\   NL%ED_MET_DRIVER_DB = '\'"${METDRIV}"\' ${FILEMAIN}
      sed -i '/NL%ED_MET_DRIVER_DB/c\   NL%ED_MET_DRIVER_DB = '\'"${METDRIV}"\' ${FILETEST}
      sed -i '/NL%ED_MET_DRIVER_DB/c\   NL%ED_MET_DRIVER_DB = '\'"${METDRIV}"\' ${FILEDBUG}
      #----- First year of meteorological driver. -----------------------------------------#
      sed -i '/NL%METCYC1/c\   NL%METCYC1 = '"${METCYCA}"                       ${FILEMAIN}
      sed -i '/NL%METCYC1/c\   NL%METCYC1 = '"${METCYCA}"                       ${FILETEST}
      sed -i '/NL%METCYC1/c\   NL%METCYC1 = '"${METCYCA}"                       ${FILEDBUG}
      #----- Last year of meteorological driver. ------------------------------------------#
      sed -i '/NL%METCYCF/c\   NL%METCYCF = '"${METCYCZ}"                       ${FILEMAIN}
      sed -i '/NL%METCYCF/c\   NL%METCYCF = '"${METCYCZ}"                       ${FILETEST}
      sed -i '/NL%METCYCF/c\   NL%METCYCF = '"${METCYCZ}"                       ${FILEDBUG}
      #----- Time-averaging flag (meteorological driver). ---------------------------------#
      sed -i '/NL%IMETAVG/c\   NL%IMETAVG = '"${IMETAVG}"                       ${FILEMAIN}
      sed -i '/NL%IMETAVG/c\   NL%IMETAVG = '"${IMETAVG}"                       ${FILETEST}
      sed -i '/NL%IMETAVG/c\   NL%IMETAVG = '"${IMETAVG}"                       ${FILEDBUG}
      #----- Radiation decomposition scheme. ----------------------------------------------#
      sed -i '/NL%IMETRAD/c\   NL%IMETRAD = '"${IMETRAD}"                       ${FILEMAIN}
      sed -i '/NL%IMETRAD/c\   NL%IMETRAD = '"${IMETRAD}"                       ${FILETEST}
      sed -i '/NL%IMETRAD/c\   NL%IMETRAD = '"${IMETRAD}"                       ${FILEDBUG}
      #----- PFTs to include in this simulation. ------------------------------------------#
      sed -i '/NL%INCLUDE_THESE_PFT/c\   NL%INCLUDE_THESE_PFT = '"${INCPFT_MAIN}"          \
                                                                                ${FILEMAIN}
      sed -i '/NL%INCLUDE_THESE_PFT/c\   NL%INCLUDE_THESE_PFT = '"${INCPFT_TEST}"          \
                                                                                ${FILETEST}
      sed -i '/NL%INCLUDE_THESE_PFT/c\   NL%INCLUDE_THESE_PFT = '"${INCPFT_TEST}"          \
                                                                                ${FILEDBUG}
      #----- PFT for pastures. ------------------------------------------------------------#
      sed -i '/NL%PASTURE_STOCK/c\   NL%PASTURE_STOCK = '"${PASTURE_MAIN}"      ${FILEMAIN}
      sed -i '/NL%PASTURE_STOCK/c\   NL%PASTURE_STOCK = '"${PASTURE_TEST}"      ${FILETEST}
      sed -i '/NL%PASTURE_STOCK/c\   NL%PASTURE_STOCK = '"${PASTURE_TEST}"      ${FILEDBUG}
      #----- PFT for croplands. -----------------------------------------------------------#
      sed -i '/NL%AGRI_STOCK/c\   NL%AGRI_STOCK = '"${AGRI_MAIN}"               ${FILEMAIN}
      sed -i '/NL%AGRI_STOCK/c\   NL%AGRI_STOCK = '"${AGRI_TEST}"               ${FILETEST}
      sed -i '/NL%AGRI_STOCK/c\   NL%AGRI_STOCK = '"${AGRI_TEST}"               ${FILEDBUG}
      #----- PFT for forest plantation. ---------------------------------------------------#
      sed -i '/NL%PLANTATION_STOCK/c\   NL%PLANTATION_STOCK = '"${PLANT_MAIN}"  ${FILEMAIN}
      sed -i '/NL%PLANTATION_STOCK/c\   NL%PLANTATION_STOCK = '"${PLANT_TEST}"  ${FILETEST}
      sed -i '/NL%PLANTATION_STOCK/c\   NL%PLANTATION_STOCK = '"${PLANT_TEST}"  ${FILEDBUG}
      #----- PFT for selective logging. ---------------------------------------------------#
      sed -i '/NL%SL_PFT/c\   NL%SL_PFT = '"${SLPFT_MAIN}"                      ${FILEMAIN}
      sed -i '/NL%SL_PFT/c\   NL%SL_PFT = '"${SLPFT_TEST}"                      ${FILETEST}
      sed -i '/NL%SL_PFT/c\   NL%SL_PFT = '"${SLPFT_TEST}"                      ${FILEDBUG}
      #----- Harvesting probability for selective logging. --------------------------------#
      sed -i '/NL%SL_PROB_HARVEST/c\   NL%SL_PROB_HARVEST = '"${SLPROBHARV_MAIN}"          \
                                                                                ${FILEMAIN}
      sed -i '/NL%SL_PROB_HARVEST/c\   NL%SL_PROB_HARVEST = '"${SLPROBHARV_TEST}"          \
                                                                                ${FILETEST}
      sed -i '/NL%SL_PROB_HARVEST/c\   NL%SL_PROB_HARVEST = '"${SLPROBHARV_TEST}"          \
                                                                                ${FILEDBUG}
      #----- Minimum DBH for selective logging. -------------------------------------------#
      sed -i '/NL%SL_MINDBH_HARVEST/c\   NL%SL_MINDBH_HARVEST = '"${SLMINDBH_MAIN}"        \
                                                                                ${FILEMAIN}
      sed -i '/NL%SL_MINDBH_HARVEST/c\   NL%SL_MINDBH_HARVEST = '"${SLMINDBH_TEST}"        \
                                                                                ${FILETEST}
      sed -i '/NL%SL_MINDBH_HARVEST/c\   NL%SL_MINDBH_HARVEST = '"${SLMINDBH_TEST}"        \
                                                                                ${FILEDBUG}
      #----- Land use data base. ----------------------------------------------------------#
      sed -i '/NL%LU_DATABASE/c\   NL%LU_DATABASE = '\'"${LU_DATABASE}"\'       ${FILEMAIN}
      sed -i '/NL%LU_DATABASE/c\   NL%LU_DATABASE = '\'"${LU_DATABASE}"\'       ${FILETEST}
      sed -i '/NL%LU_DATABASE/c\   NL%LU_DATABASE = '\'"${LU_DATABASE}"\'       ${FILEDBUG}
      #----- Phenology data base. ---------------------------------------------------------#
      sed -i '/NL%PHENPATH/c\   NL%PHENPATH = '\'"${PHEN_PATH}"\'               ${FILEMAIN}
      sed -i '/NL%PHENPATH/c\   NL%PHENPATH = '\'"${PHEN_PATH}"\'               ${FILETEST}
      sed -i '/NL%PHENPATH/c\   NL%PHENPATH = '\'"${PHEN_PATH}"\'               ${FILEDBUG}
      #----- First year of the phenology (Spring). ----------------------------------------#
      sed -i '/NL%IPHENYS1/c\   NL%IPHENYS1 = '"${IPHENYSA}"                    ${FILEMAIN}
      sed -i '/NL%IPHENYS1/c\   NL%IPHENYS1 = '"${IPHENYSA}"                    ${FILETEST}
      sed -i '/NL%IPHENYS1/c\   NL%IPHENYS1 = '"${IPHENYSA}"                    ${FILEDBUG}
      #----- Last year of the phenology (Spring). -----------------------------------------#
      sed -i '/NL%IPHENYSF/c\   NL%IPHENYSF = '"${IPHENYSZ}"                    ${FILEMAIN}
      sed -i '/NL%IPHENYSF/c\   NL%IPHENYSF = '"${IPHENYSZ}"                    ${FILETEST}
      sed -i '/NL%IPHENYSF/c\   NL%IPHENYSF = '"${IPHENYSZ}"                    ${FILEDBUG}
      #----- First year of the phenology (Autumn). ----------------------------------------#
      sed -i '/NL%IPHENYF1/c\   NL%IPHENYF1 = '"${IPHENYFA}"                    ${FILEMAIN}
      sed -i '/NL%IPHENYF1/c\   NL%IPHENYF1 = '"${IPHENYFA}"                    ${FILETEST}
      sed -i '/NL%IPHENYF1/c\   NL%IPHENYF1 = '"${IPHENYFA}"                    ${FILEDBUG}
      #----- Last year of the phenology (Autumn). -----------------------------------------#
      sed -i '/NL%IPHENYFF/c\   NL%IPHENYFF = '"${IPHENYFZ}"                    ${FILEMAIN}
      sed -i '/NL%IPHENYFF/c\   NL%IPHENYFF = '"${IPHENYFZ}"                    ${FILETEST}
      sed -i '/NL%IPHENYFF/c\   NL%IPHENYFF = '"${IPHENYFZ}"                    ${FILEDBUG}
      #----- XML file. --------------------------------------------------------------------#
      sed -i '/NL%IEDCNFGF/c\   NL%IEDCNFGF = '\'"${EDCNFGF}"\'                 ${FILEMAIN}
      sed -i '/NL%IEDCNFGF/c\   NL%IEDCNFGF = '\'"${EDCNFGF}"\'                 ${FILETEST}
      sed -i '/NL%IEDCNFGF/c\   NL%IEDCNFGF = '\'"${EDCNFGF}"\'                 ${FILEDBUG}
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #      Detailed output settings.                                                     #
      #------------------------------------------------------------------------------------#
      #----- Detailed output list. --------------------------------------------------------#
      sed -i '/NL%IDETAILED/c\   NL%IDETAILED = 0'                          ${FILEMAIN}
      sed -i '/NL%IDETAILED/c\   NL%IDETAILED = 0'                          ${FILETEST}
      sed -i '/NL%IDETAILED/c\   NL%IDETAILED = 0'                          ${FILEDBUG}
      #----- Patches to keep. -------------------------------------------------------------#
      sed -i '/NL%PATCH_KEEP/c\   NL%PATCH_KEEP = 0'                        ${FILEMAIN}
      sed -i '/NL%PATCH_KEEP/c\   NL%PATCH_KEEP = 0'                        ${FILETEST}
      sed -i '/NL%PATCH_KEEP/c\   NL%PATCH_KEEP = 0'                        ${FILEDBUG}
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


      #----- Modify the ED2IN start times. ------------------------------------------------#
      sed -i '/NL%IYEARA/c\   NL%IYEARA   = '${IYEARAG[i]}    ${FILEMAIN}
      sed -i '/NL%IYEARA/c\   NL%IYEARA   = '${IYEARAG[i]}    ${FILETEST}
      sed -i '/NL%IYEARA/c\   NL%IYEARA   = '${D_IYEARAG[i]}  ${FILEDBUG}
      sed -i '/NL%IMONTHA/c\   NL%IMONTHA   = '${IMONTHAG[i]} ${FILEMAIN}
      sed -i '/NL%IMONTHA/c\   NL%IMONTHA   = '${IMONTHAG[i]} ${FILETEST}
      sed -i '/NL%IMONTHA/c\   NL%IMONTHA   = '${IMONTHAG[i]} ${FILEDBUG}
      sed -i '/NL%IDATEA/c\   NL%IDATEA   = '${IDATEAG[i]}    ${FILEMAIN}
      sed -i '/NL%IDATEA/c\   NL%IDATEA   = '${IDATEAG[i]}    ${FILETEST}
      sed -i '/NL%IDATEA/c\   NL%IDATEA   = '${IDATEAG[i]}    ${FILEDBUG}
      sed -i '/NL%ITIMEA/c\   NL%ITIMEA   = 0000'             ${FILEMAIN}
      sed -i '/NL%ITIMEA/c\   NL%ITIMEA   = 0000'             ${FILETEST}
      sed -i '/NL%ITIMEA/c\   NL%ITIMEA   = 0000'             ${FILEDBUG}
      #------------------------------------------------------------------------------------#


      #----- Modify the ED2IN end years. --------------------------------------------------#
      sed -i '/NL%IYEARZ/c\   NL%IYEARZ   = '${IYEARZG[i]}    ${FILEMAIN}
      sed -i '/NL%IYEARZ/c\   NL%IYEARZ   = '${IYEARZG[i]}    ${FILETEST}
      sed -i '/NL%IYEARZ/c\   NL%IYEARZ   = '${D_IYEARZG[i]}  ${FILEDBUG}
      sed -i '/NL%IMONTHZ/c\   NL%IMONTHZ   = '${IMONTHZG[i]} ${FILEMAIN}
      sed -i '/NL%IMONTHZ/c\   NL%IMONTHZ   = '${IMONTHZG[i]} ${FILETEST}
      sed -i '/NL%IMONTHZ/c\   NL%IMONTHZ   = '${IMONTHZG[i]} ${FILEDBUG}
      sed -i '/NL%IDATEZ/c\   NL%IDATEZ   = '${IDATEZG[i]}    ${FILEMAIN}
      sed -i '/NL%IDATEZ/c\   NL%IDATEZ   = '${IDATEZG[i]}    ${FILETEST}
      sed -i '/NL%IDATEZ/c\   NL%IDATEZ   = '${IDATEZG[i]}    ${FILEDBUG}
      sed -i '/NL%ITIMEZ/c\   NL%ITIMEZ   = 0000'             ${FILEMAIN}
      sed -i '/NL%ITIMEZ/c\   NL%ITIMEZ   = 0000'             ${FILETEST}
      sed -i '/NL%ITIMEZ/c\   NL%ITIMEZ   = 0000'             ${FILEDBUG}
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
      #    Define the MPI command based on the MPI Type.                                   #
      #------------------------------------------------------------------------------------#
      MPI_COMM="mpirun -np ${GRIDCPU[i]}"
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #    Define the node settings.                                                       #
      #------------------------------------------------------------------------------------#
      case "${GRID_NODE}" in
         0)  NNODES_REQ=""                ;;
         *)  NNODES_REQ="-N ${GRID_NODE}" ;;
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
         jobopts="${jobopts} -p ${GRIDQ[i]} -n ${GRIDCPU[i]} ${NNODES_REQ}"
         jobwrap=". ${HOME}/.bashrc ${RC_OPTS}; cd ${HERE}/${VERSION}"
         jobwrap="${jobwrap}; export OMP_NUM_THREADS=1"
         jobwrap="${jobwrap}; ${MPI_COMM} ${LNK_MAIN_EXE} -f ${FILEMAIN}"
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
         jobopts="${jobopts} -p ${GRIDQ[i]} -n ${GRIDCPU[i]} ${NNODES_REQ}"
         jobwrap=". ${HOME}/.bashrc ${RC_OPTS}; cd ${HERE}/${VERSION}"
         jobwrap="${jobwrap}; export OMP_NUM_THREADS=1"
         jobwrap="${jobwrap}; ${MPI_COMM} ${LNK_TEST_EXE} -f ${FILETEST}"
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
         jobopts="${jobopts} -p ${GRIDQ[i]} -n ${GRIDCPU[i]} ${NNODES_REQ}"
         jobwrap=". ${HOME}/.bashrc ${RC_OPTS}; cd ${HERE}/${VERSION}"
         jobwrap="${jobwrap}; export OMP_NUM_THREADS=1"
         jobwrap="${jobwrap}; ${MPI_COMM} ${LNK_DBUG_EXE} -f ${FILEDBUG}"
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
