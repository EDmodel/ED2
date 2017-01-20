#!/bin/bash
################ COMPILER OPTIONS (INTEL, ODYSSEY, SUN-HPC, SDUMONT ONLY) ##################
#------------------------------------------------------------------------------------------#
# A/B. Pickiest.  Use this whenever you change arguments on functions and subroutines, and #
#      perform the three steps to make sure the code has consistent interfaces (i.e.       #
#      arguments match between subroutine declaration and subroutine calls).  A and B are  #
#      the same now (stick with one), but you must provide another argument, -s <value>,   #
#      which will tell which step of the compilation is to be run.                         #
#      Step 1. Compile and generate interfaces                                             #
#              (e.g. ./install.sh -k A -p odyssey -s 1)                                    #
#      Step 2. Remove original files, but keep modules (only run this step after you       #
#              create the executable using step 1).                                        #
#              (e.g. ./install.sh -k A -p odyssey -s 2)                                    #
#      Step 3. Compile again with the interface verification.                              #
#              (e.g. ./install.sh -k A -p odyssey -s 3)                                    #
#                                                                                          #
# C. Pickiest with no interface - This will compile fast but the run will be slow due to   #
#    the -O0 option. However, by setting -O0 you will take full advantage of the intel     #
#    debugger.                                                                             #
#    Ideally, you should compile your code with this option whenever you make any changes. #
#    Note, however, that if you change arguments you should first try A.                   #
#                                                                                          #
# D. Fast check - This will check pretty much the same as C, but it will not set up        #
#    anything for the debugger. Use this only if you really don't want to deal with idb or #
#    if you have a good idea of which problem you are dealing with.                        #
#                                                                                          #
# E. Fast - This is all about performance, use only when you are sure that the model has   #
#           no code problem, and you want results fast. This will not check for any        #
#           problems, which means that this is an OPTION SUITABLE for end users.           #
#           DEVELOPERS, THIS IS NOT THE OPTION YOU SHOULD USE TO CHECK WHETHER YOUR        #
#           NEW CODE IS WORKING OR NOT.                                                    #
#------------------------------------------------------------------------------------------#

#----- Define the number of arguments. ----------------------------------------------------#
nargs=$#
args=$@
#------------------------------------------------------------------------------------------#

# Initialize vars
CLEAN=""
KIND=""
PLATFORM=""
OPT=""
USE_GIT=true
STEP=0

# Argument parsing
while [[ $# > 0 ]]
do
key="$1"
   case $key in
   -p|--platform)
      PLATFORM="$2"
      shift # past argument
      ;;
   -k|--kind)
      KIND="$2"
      shift # past argument
      ;;
   -c|--clean)
      CLEAN="clean"
      ;;
   -g|--gitoff)
      USE_GIT=false
      ;;
   -s|--step)
      STEP="$2"
      shift # past argument
      ;;
   *)
      echo "Unknown key-value argument pair."
      exit 2
      ;;
   esac

   shift # past argument or value
done

if [ "x${PLATFORM}" == "x" ]
then
   echo "No platform specified, defaulting to gfortran."
   PLATFORM="gfortran"
fi

if [ "x${KIND}" == "x" ]
then  
   echo "No optimization level specified, defaulting to E."
   KIND="E"
fi


# Check that the step is properly set 
case ${KIND} in
['A','B']*)
   case ${STEP} in
   0)
      echo "You must provide step (option -s or --step) when use \"-k A\" or \"-k B\""
      exit 1
      ;;
   1)
      LKIND="A"
      echo "Step ${STEP}, generate interfaces."
      ;;
   2)
      LKIND="A"
      echo "Step ${STEP}, prepare for interface check."
      ;;
   3)
      LKIND="B"
      echo "Step ${STEP}, check interfaces."
      ;;
   *)
      echo "Invalid step! It must be one of the following options:"
      echo "-s 1 (--step 1): Generate interfaces"
      echo "-s 2 (--step 2): Prepare for interface check"
      echo "-s 3 (--step 3): Check interfaces"
      exit 1
      ;;
   esac
   ;;
*)
   LKIND=${KIND}
   ;;
esac

      
# Set opt and bin
case ${KIND} in
['A','B','C','D']*)
   OPT='dbg'
   ;;
['E']*)
   OPT='opt'
   ;;
*)
   # Default to opt
   echo "Compiler optimization not recognized as opt or dbg."
   exit 1
   ;;
esac

# Tag executables with a git version and branch name if possible.
GIT_EXIST=$(git rev-parse --is-inside-work-tree)
if [ "x${GIT_EXIST}" == "xtrue" ] && ${USE_GIT}
then
   GIT_TAG=$(git branch -v | awk '/\*/ {print "-" $2 "-" $3}')
   GIT_TAG=$(echo ${GIT_TAG} | tr -d '()/[]')
   echo "Git found, it will be used to tag things."
   echo "To disable revision tagging, use --gitoff or -g."
else
   GIT_TAG=""
fi

BIN=bin-${OPT}-${KIND}${GIT_TAG}

# Move to the binary directory
if [ ! -d "${BIN}" ]; then
   mkdir ${BIN}
fi
cd ${BIN}

case ${STEP} in
2)
   ./2ndcomp.sh
   ;;
*)
   # Link to makefiles, includes, and shell scripts
   ln -sf ../make/*.mk ./
   ln -sf ../make/Makefile ./
   ln -sf ../make/include.mk.${OPT}.${PLATFORM} ./include.mk
   ln -sf ../shell/* ./
   touch dependency.mk

   #----- Launch the compiler. ------------------------------------------------------------#
   make OPT=${OPT} KIND_COMP=${LKIND} ${CLEAN} GIT_TAG=${GIT_TAG}
   make_exit_code=$?
   #---------------------------------------------------------------------------------------#

   if [ ${make_exit_code} != 0 ]
   then
      exit 1
   else
      echo "Installation Complete."
      exit 0
   fi
   ;;
esac
