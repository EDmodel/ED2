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
#------ System. ---------------------------------------------------------------------------#
myos=$(uname -s)
#------ How detailed should the check be (partial - c, f90, and F90; total - everything. --#
comptype="partial"
#------ Sub-directories to be compared. ---------------------------------------------------#
subdirs="ED BRAMS Ramspost R-utils"
#------ Editor to use (I've only tested with nedit, use others at your own risk). ---------#
editor="nedit"
#------ Two paths with EDBRAMS (full path). -----------------------------------------------#
ours="${HOME}/Models/EDBRAMS"
theirs="${HOME}/Models/MainLine-ED2"
ournew=true
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#      Extensions.  Decide how much to check.                                              #
#------------------------------------------------------------------------------------------#
case "${comptype}" in
partial) exts=".c .f90 .F90" ;;
*)       exts=".c .f90 .F90 .r .txt" ;;
esac
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#      Remove any comparison file that may exist.                                          #
#------------------------------------------------------------------------------------------#
/bin/rm -f notthesame.*
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
        srcours="${ours}/${subdir}/samap"
        srctheirs="${theirs}/${subdir}/samap"
     elif [ ${subdir} == "R-utils" ]
     then
        srcours="${ours}/${subdir}"
        srctheirs="${theirs}/${subdir}"
     else
        srcours="${ours}/${subdir}/src"
        srctheirs="${theirs}/${subdir}/src"
     fi
     #-------------------------------------------------------------------------------------#


     #-------------------------------------------------------------------------------------#
     #    Look for files that were changed or were removed in the remote test.             #
     #-------------------------------------------------------------------------------------#
     lookuptable=$(find ${srcours} -name "*${ext}")
     for fileours in ${lookuptable}
     do
       file=$(basename ${fileours})
       newpath=$(dirname ${fileours} | sed s@${ours}@${theirs}@g)
       case "${myos}" in
       Darwin|Cygwin)
          #---- Case-insensitive file system. ---------------------------------------------#
          filetheirs="${newpath}/${file}"
          alt=""
          ;;
          #--------------------------------------------------------------------------------#
       *)
          #----- Probably Linux or Unix, assume case-sensitive file system. ---------------#
          case "${ext}" in
          ".f90")
             alt_file=$(echo ${file} | sed s@"\\.f90"@".F90"@g)
             if [ -s "${newpath}/${alt_file}" ]
             then
                filetheirs="${newpath}/${alt_file}"
                alt=".F90"
             else
                filetheirs="${newpath}/${file}"
                alt=""
             fi
             ;;
          ".F90")
             alt_file=$(echo ${file} | sed s@"\\.F90"@".f90"@g)
             if [ -s "${newpath}/${alt_file}" ]
             then
                filetheirs="${newpath}/${alt_file}"
                alt=".f90"
             else
                filetheirs="${newpath}/${file}"
                alt=""
             fi
             ;;
          *)
             filetheirs="${newpath}/${file}"
             alt=""
             ;;
          esac
          ;;
          #--------------------------------------------------------------------------------#
       esac
       filetheirs="${newpath}/${file}"
       if [ -s ${filetheirs} ] 
       then
         ldif=$(diff -ibB <(grep -vE "^\s*!" ${fileours}) <(grep -vE "^\s*!" ${filetheirs}) | wc -l)
         if [ ${ldif} -gt 0 ]
         then
            woroot=$(echo ${fileours} | sed s@"${srcours}/"@""@g)
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
               diff -uibB <(grep -vE "^\s*!" ${filetheirs}) <(grep -vE "^\s*!" ${fileours})\
                  > notthesame${ext}
            else
               diff -uibB <(grep -vE "^\s*!" ${fileours}) <(grep -vE "^\s*!" ${filetheirs})\
                  > notthesame${ext}
            fi
            ${editor} ${fileours} ${filetheirs} notthesame${ext} 1> /dev/null 2> /dev/null
            /bin/rm -f notthesame${ext}
         fi
       else
           woroot=$(echo ${fileours} | sed s@"${srcours}/"@""@g)
           echo "${woroot} is exclusive to ${ours} version."
       fi #if [ -s ${filetheirs} ]
     done #for fileours in ${lookuptable}
     #-------------------------------------------------------------------------------------#



     #-------------------------------------------------------------------------------------#
     #    Look for files that were changed or were removed in the remote test.             #
     #-------------------------------------------------------------------------------------#
     lookuptable=$(find ${srctheirs} -name "*${ext}")
     for filetheirs in ${lookuptable}
     do
       file=$(basename ${filetheirs})
       newpath=$(dirname ${filetheirs} | sed s@${theirs}@${ours}@g)
       fileours="${newpath}/${file}"
       if [ ! -s ${fileours} ] 
       then
           woroot=$(echo ${filetheirs} | sed s@"${srctheirs}/"@""@g)
           echo "${woroot} is exclusive to ${theirs} version."
       fi #if [ -s ${filetheirs} ]
     done #for fileours in ${lookuptable}
     #-------------------------------------------------------------------------------------#
   done #for ext in ${exts}
   #---------------------------------------------------------------------------------------#
done #for subdir in ${subdirs}
#------------------------------------------------------------------------------------------#
