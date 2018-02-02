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
our="${HOME}/EDBRAMS"
their="${HOME}/MainLine/pr_EDBRAMS"
ournew=true
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
        srcour="${our}/${subdir}/samap"
        srctheir="${their}/${subdir}/samap"
     elif [ ${subdir} == "R-utils" ]
     then
        srcour="${our}/${subdir}"
        srctheir="${their}/${subdir}"
     else
        srcour="${our}/${subdir}/src"
        srctheir="${their}/${subdir}/src"
     fi
     #-------------------------------------------------------------------------------------#


     #-------------------------------------------------------------------------------------#
     #    Look for files that were changed or were removed in the remote test.             #
     #-------------------------------------------------------------------------------------#
     lookuptable=$(find ${srcour} -name "*${ext}")
     for fileour in ${lookuptable}
     do
       file=$(basename ${fileour})
       newpath=$(dirname ${fileour} | sed s@${our}@${their}@g)
       case ${ext} in
       ".f90")
          alt_file=$(echo ${file} | sed s@"\\.f90"@".F90"@g)
          if [ -s "${newpath}/${alt_file}" ]
          then
             filetheir="${newpath}/${alt_file}"
             alt=".F90"
          else
             filetheir="${newpath}/${file}"
             alt=""
          fi
          ;;
       ".F90")
          alt_file=$(echo ${file} | sed s@"\\.F90"@".f90"@g)
          if [ -s "${newpath}/${alt_file}" ]
          then
             filetheir="${newpath}/${alt_file}"
             alt=".f90"
          else
             filetheir="${newpath}/${file}"
             alt=""
          fi
          ;;
       *)
          filetheir="${newpath}/${file}"
          alt=""
          ;;
       esac
       filetheir="${newpath}/${file}"
       if [ -s ${filetheir} ] 
       then
         ldif=$(diff -ibB <(grep -vE "^\s*!" ${fileour}) <(grep -vE "^\s*!" ${filetheir}) | wc -l)
         if [ ${ldif} -gt 0 ]
         then
            woroot=$(echo ${fileour} | sed s@"${srcour}/"@""@g)
            case "${alt}" in
            .f90|.F90)
               echo "${woroot} has changed.  New extension is ${alt}."
               ;;
            *)
               echo "${woroot} has changed."
               ;;
            esac
            if ${ournew}
            then
               diff -uibB <(grep -vE "^\s*!" ${filetheir}) <(grep -vE "^\s*!" ${fileour})  \
                  > notthesame.txt
            else
               diff -uibB <(grep -vE "^\s*!" ${fileour}) <(grep -vE "^\s*!" ${filetheir})  \
                  > notthesame.txt
            fi
            ${editor} ${fileour} ${filetheir} notthesame.txt 1> /dev/null 2> /dev/null
         fi
       else
           woroot=$(echo ${fileour} | sed s@"${srcour}/"@""@g)
           echo "${woroot} is exclusive to ${our} version."
       fi #if [ -s ${filetheir} ]
     done #for fileour in ${lookuptable}
     #-------------------------------------------------------------------------------------#



     #-------------------------------------------------------------------------------------#
     #    Look for files that were changed or were removed in the remote test.             #
     #-------------------------------------------------------------------------------------#
     lookuptable=$(find ${srctheir} -name "*${ext}")
     for filetheir in ${lookuptable}
     do
       file=$(basename ${filetheir})
       newpath=$(dirname ${filetheir} | sed s@${their}@${our}@g)
       fileour="${newpath}/${file}"
       if [ ! -s ${fileour} ] 
       then
           woroot=$(echo ${filetheir} | sed s@"${srctheir}/"@""@g)
           echo "${woroot} is exclusive to ${their} version."
       fi #if [ -s ${filetheir} ]
     done #for fileour in ${lookuptable}
     #-------------------------------------------------------------------------------------#
   done #for ext in ${exts}
   #---------------------------------------------------------------------------------------#
done #for subdir in ${subdirs}
#------------------------------------------------------------------------------------------#
