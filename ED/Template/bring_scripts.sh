#!/bin/bash

#----- Save script name. ------------------------------------------------------------------#
this_script=$(basename ${0})
#------------------------------------------------------------------------------------------#


#----- Find out which platform we are using. ----------------------------------------------#
if [ "${1}" == "" ]
then
   #------ No platform provided.  Try to guess, and if failed, then prompts the user. -----#
   host=$(hostname -s)
   case ${host} in
      rclogin*|holy*|moorcroft*|rcnx*) platform="odyssey" ;;
      aurora|halo)                     platform="jpl"     ;;
      sdumont*)                        platform="sdumont" ;;
      sun-master|cmm*)                 platform="sun-lncc" ;;
      *)
         echo -n "Failed guessing platform from node name.  Please type the name:   "
         read platform
         ;;
   esac

else
   platform=${1}
fi
#------------------------------------------------------------------------------------------#



#----- Bring scripts from directory to the main directory, and delete the others. ---------#
sh_here=$(/bin/ls -1 *.sh | grep -v ${this_script})
if [ "${sh_here}" != "" ]
then
   echo -n "Shell scripts were found in this directory.  Replace them? [y|N]  "
   read reset
   echo ""
   reset=$(echo ${reset} | tr '[:upper:]' '[:lower:]')
   echo ${reset}
   
   case ${reset} in
   y|yes) 
      #----- Delete shell scripts. --------------------------------------------------------#
      for shell in ${sh_here}
      do
         echo " Delete current file ${shell}."
         find ${shell} -delete
      done
      #------------------------------------------------------------------------------------#

      ;;
   *)
      echo " "
      echo "Cannot copy scripts here."
      echo "Either delete old shell scripts or let me overwrite."
      exit 0
      ;;
   esac
fi
#------------------------------------------------------------------------------------------#


#----- Copy all scripts from the script pool. ---------------------------------------------#
/bin/cp -v ./scripts/${platform}/*.sh .
#------------------------------------------------------------------------------------------#
