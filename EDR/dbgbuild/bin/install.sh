#!/bin/bash
######################### COMPILER OPTIONS (INTEL AND ODYSSEY ONLY) ########################
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
#------------------------------------------------------------------------------------------#

#----- Define the number of arguments. ----------------------------------------------------#
nargs=$#
args=$@
#------------------------------------------------------------------------------------------#


#----- Default parameters. ----------------------------------------------------------------#
kind="B"
clean=""
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Check whether the user gave options.                                                 #
#------------------------------------------------------------------------------------------#
if [ ${nargs} -gt 0 ]
then
   for arg in ${args}
   do
      #----- Check whether which argument this is. ----------------------------------------#
      if [ "x${arg}" == "xclean" ]
      then
         clean="clean"
      else
         kind=`echo ${arg} | tr '[:lower:]' '[:upper:]'`
      fi
      #------------------------------------------------------------------------------------#
   done
   #---------------------------------------------------------------------------------------#
fi
#------------------------------------------------------------------------------------------#



#----- Launch the compiler. ---------------------------------------------------------------#
make OPT=opt KIND_COMP=${kind} ${clean}
#------------------------------------------------------------------------------------------#
