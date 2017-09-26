#!/bin/bash
################ COMPILER OPTIONS (INTEL, ODYSSEY, and EMBRAPA ONLY) #######################
#------------------------------------------------------------------------------------------#
# A/B. Pickiest - Use this whenever you change arguments on functions and subroutines.     #
#                 This will perform the same tests as B but it will also check whether all #
#                 arguments match between subroutine declaration and subroutine calls.     #
#                 WARNING: In order to really check all interfaces you must compile with   #
#                          this option twice:                                              #
#                 1. Compile (./install.sh A)                                              #
#                 2. Prepare second compilation(./2ndcomp.sh)                              #
#                 3. Compile one more time (./install.sh B)                                #
#                 If the compilation fails either at step 3, then your code has interface  #
#                 problems. If it successfully compiles, then the code is fine for         #
#                 interfaces.                                                              #
# C. Pickiest with no interface - This will compile fast but the run will be slow due to   #
#    the -O0 option. However, by setting -O0 you will take full advantage of the intel     #
#    debugger.                                                                             #
#    Ideally, you should compile your code with this option whenever you make any changes. #
#    Note, however, that if you change arguments you should first try A.                   #
# D. Fast check - This will check pretty much the same as C, but it will not set up        #
#    anything for the debugger. Use this only if you really don't want to deal with idb or #
#    if you have a good idea of which problem you are dealing with.                        #
# E. Fast - This is all about performance, use only when you are sure that the model has   #
#           no code problem, and you want results asap. This will not check for any        #
#           problems, which means that this is an option suitable for end users, not de-   #
#           velopers.                                                                      #
#                                                                                          #
#                                                                                          #
# usage example:                                                                           #
#                                                                                          #
# install.sh --kind A --platform intel                                                     #
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
   *)
      echo "Unknown key-value argument pair."
      exit 2
      ;;
   esac

   shift # past argument or value
done

if [ "${PLATFORM}" == "" ]
then
   echo "No platform specified, defaulting to intel."
   PLATFORM="intel"
fi

if [ "${KIND}" == "" ]
then  
   echo "No optimization level specified, defaulting to E."
   KIND="E"
fi

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
GIT_EXIST=`git rev-parse --is-inside-work-tree`
if [ ${GIT_EXIST} == "true" -a ${USE_GIT} ]
then
   GIT_TAG=`git branch -v | awk '/\*/ {print "-" $2 "-" $3}'`
   GIT_TAG=`echo ${GIT_TAG} | tr -d '()/[]'`
   echo "Git found, it will be used to tag things."
   echo "To disable revision tagging, use --gitoff or -g."
else
   GIT_TAG=''
fi

BIN=bin-${OPT}-${KIND}${GIT_TAG}

# Move to the binary directory
if [ ! -d "$BIN" ]; then
   mkdir ${BIN}
fi
cd ${BIN}


# Link to makefiles, includes, and shell scripts
ln -sf ../make/*.mk ./
ln -sf ../make/Makefile ./
ln -sf ../make/include.mk.${PLATFORM} ./include.mk
ln -sf ../shell/* ./
touch dependency.mk

#----- Launch the compiler. ---------------------------------------------------------------#
make OPT=${OPT} KIND_COMP=${KIND} ${CLEAN} GIT_TAG=${GIT_TAG}
make_exit_code=$?
#------------------------------------------------------------------------------------------#

if [ ${make_exit_code} != 0 ]
then
   exit 1
else
   echo "Installation Complete."
   exit 0
fi
