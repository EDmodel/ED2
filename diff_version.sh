#!/bin/bash

#------------------------------------------------------------------------------------------#
#       This is how it works:  you give two directories with EDBRAMS.  It will go through  #
# all files and check whether the files exists on both locations.                          #
# 1.  If they do and are exactly the same, nothing happens;                                #
# 2.  If they do and they are different, the editor will open three files, two of them are #
#     the source code at the locations, and a third with the result of "diff -ib".  Look   #
#     at the file, when you are done, close all the three files so it continues.           #
# 3.  If the file exists in only one location, a message will appear on screen.            #
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#      Settings.                                                                           #
#------------------------------------------------------------------------------------------#
#------ Extensions. -----------------------------------------------------------------------#
exts=".c .f90 .F90 .r .txt"
#------ Sub-directories to be compared. ---------------------------------------------------#
subdirs="ED BRAMS Ramspost R-utils"
#------ Editor to use (I've only tested with nedit, use others at your own risk). ---------#
editor="nedit"
#------ Two paths with EDBRAMS (full path). -----------------------------------------------#
here="/n/home00/mlongo/EDBRAMS"
there="/n/home00/mlongo/MainLine/mpaiao-EDBRAMS"
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#      Remove any comparison file that may exist.                                          #
#------------------------------------------------------------------------------------------#
/bin/rm -f notthesame.txt
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#      Loop over the sub-directories.                                                      #
#------------------------------------------------------------------------------------------#
for subdir in ${subdirs}
do
   #---------------------------------------------------------------------------------------#
   #      Loop over the extensions.                                                        #
   #---------------------------------------------------------------------------------------#
   for ext in ${exts}
   do
     if [ ${subdir} == "R-utils" ] && [ ${ext} == ".txt" ]
     then
        srchere="${here}/${subdir}/samap"
        srcthere="${there}/${subdir}/samap"
     elif [ ${subdir} == "R-utils" ]
     then
        srchere="${here}/${subdir}"
        srcthere="${there}/${subdir}"
     else
        srchere="${here}/${subdir}/src"
        srcthere="${there}/${subdir}/src"
     fi
     #-------------------------------------------------------------------------------------#


     #-------------------------------------------------------------------------------------#
     #    Look for files that were changed or were removed in the remote test.             #
     #-------------------------------------------------------------------------------------#
     lookuptable=$(find ${srchere} -name "*${ext}")
     for filehere in ${lookuptable}
     do
       file=$(basename ${filehere})
       newpath=$(dirname ${filehere} | sed s@${here}@${there}@g)
       filethere="${newpath}/${file}"
       if [ -s ${filethere} ] 
       then
         ldif=$(diff -ib ${filethere} ${filehere} | wc -l)
         if [ ${ldif} -gt 0 ]
         then
            woroot=$(echo ${filehere} | sed s@"${srchere}/"@""@g)
            echo "${woroot} has changed..."
            diff -ib ${filehere} ${filethere} > notthesame.txt
            ${editor} ${filehere} ${filethere} notthesame.txt 1> /dev/null 2> /dev/null
         fi
       else
           woroot=$(echo ${filehere} | sed s@"${srchere}/"@""@g)
           echo "${woroot} is exclusive to ${here} version..."
       fi #if [ -s ${filethere} ]
     done #for filehere in ${lookuptable}
     #-------------------------------------------------------------------------------------#



     #-------------------------------------------------------------------------------------#
     #    Look for files that were changed or were removed in the remote test.             #
     #-------------------------------------------------------------------------------------#
     lookuptable=$(find ${srcthere} -name "*${ext}")
     for filethere in ${lookuptable}
     do
       file=$(basename ${filethere})
       newpath=$(dirname ${filethere} | sed s@${there}@${here}@g)
       filehere="${newpath}/${file}"
       if [ ! -s ${filehere} ] 
       then
           woroot=$(echo ${filethere} | sed s@"${srcthere}/"@""@g)
           echo "${woroot} is exclusive to ${there} version..."
       fi #if [ -s ${filethere} ]
     done #for filehere in ${lookuptable}
     #-------------------------------------------------------------------------------------#
   done #for ext in ${exts}
   #---------------------------------------------------------------------------------------#
done #for subdir in ${subdirs}
#------------------------------------------------------------------------------------------#
