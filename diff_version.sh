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
ours="${HOME}/Downloads/MLN-EDBRAMS"
theirs="${HOME}/Downloads/MLO-EDBRAMS"
ournew=false
#------ Working path when implementing a partial merge (otherwise, leave it blank). -------#
work="${HOME}/Models/EDBRAMS"
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#      Extensions.  Decide how much to check.                                              #
#------------------------------------------------------------------------------------------#
case "${comptype}" in
partial) exts=".c .F90 .f90" ;;
*)       exts=".c .F90 .f90 .r .txt" ;;
esac
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#      Remove any comparison file that may exist.                                          #
#------------------------------------------------------------------------------------------#
/bin/rm -f notthesame.* \~notthesame.*
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#    Add settings for partial merge (if sought).                                           #
#------------------------------------------------------------------------------------------#
if [[ "${work}" == "" ]]
then
   #---- Disable work version (preferred method, much less messy). ------------------------#
   load_work=false
   work=${theirs}
   #---------------------------------------------------------------------------------------#
else
   #---- Enable work version (more prone to complications, be careful). -------------------#
   load_work=true
   #---------------------------------------------------------------------------------------#
fi
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
     if [[ ${subdir} == "R-utils" ]] && [[ ${ext} == ".txt" ]]
     then
        srcours="${ours}/${subdir}/samap"
        srctheirs="${theirs}/${subdir}/samap"
     elif [[ ${subdir} == "R-utils" ]]
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
       theirpath=$(dirname ${fileours} | sed s@${ours}@${theirs}@g)
       workpath=$(dirname ${fileours} | sed s@${ours}@${work}@g)
       case "${myos}" in
       Darwin|Cygwin)
          #---- Case-insensitive file system. ---------------------------------------------#
          filetheirs="${theirpath}/${file}"
          filework="${workpath}/${file}"
          theirs_alt=""
          ;;
          #--------------------------------------------------------------------------------#
       *)
          #----- Probably Linux or Unix, assume case-sensitive file system. ---------------#
          case "${ext}" in
          ".f90")
             alt_file=$(echo ${file} | sed s@"\\.f90"@".F90"@g)
             if [[ -s "${theirpath}/${alt_file}" ]]
             then
                filetheirs="${theirpath}/${alt_file}"
                theirs_alt=".F90"
             else
                filetheirs="${theirpath}/${file}"
                theirs_alt=""
             fi
             ;;
          ".F90")
             alt_file=$(echo ${file} | sed s@"\\.F90"@".f90"@g)
             if [[ -s "${theirpath}/${alt_file}" ]]
             then
                filetheirs="${theirpath}/${alt_file}"
                theirs_alt=".f90"
             else
                filetheirs="${theirpath}/${file}"
                theirs_alt=""
             fi
             ;;
          *)
             filetheirs="${theirpath}/${file}"
             theirs_alt=""
             ;;
          esac
          #--------------------------------------------------------------------------------#


          #----- Probably Linux or Unix, assume case-sensitive file system. ---------------#
          case "${ext}" in
          ".f90")
             alt_file=$(echo ${file} | sed s@"\\.f90"@".F90"@g)
             if [[ -s "${workpath}/${alt_file}" ]]
             then
                filework="${workpath}/${alt_file}"
                work_alt=".F90"
             else
                filework="${workpath}/${file}"
                work_alt=""
             fi
             ;;
          ".F90")
             alt_file=$(echo ${file} | sed s@"\\.F90"@".f90"@g)
             if [[ -s "${workpath}/${alt_file}" ]]
             then
                filework="${workpath}/${alt_file}"
                work_alt=".f90"
             else
                filework="${workpath}/${file}"
                work_alt=""
             fi
             ;;
          *)
             filework="${workpath}/${file}"
             work_alt=""
             ;;
          esac
          #--------------------------------------------------------------------------------#

          ;;
       esac


       #-----------------------------------------------------------------------------------#
       #    If we are doing partial merge, we must check the existence three-way,          #
       # otherwise we compare just the "ours" and "theirs".                                #
       #-----------------------------------------------------------------------------------#
       if [[ -s ${filetheirs} ]]
       then
          #---- Check for extension change. -----------------------------------------------#
          woroot=$(echo ${fileours} | sed s@"${srcours}/"@""@g)
          case "${theirs_alt}" in
          .f90|.F90)
             if ${ournew}
             then
                ext_mess="Old extension was ${theirs_alt}."
             else
                ext_mess="New extension is ${theirs_alt}."
             fi
             ;;
          *)
             ext_mess=""
             ;;
          esac
          #--------------------------------------------------------------------------------#



          #--------------------------------------------------------------------------------#
          #     Check for changes (ignoring spaces and commented lines).                   #
          #--------------------------------------------------------------------------------#
          ldif=$(diff -ibB <(grep -vE "^\s*!" ${fileours}) <(grep -vE "^\s*!" ${filetheirs}) | wc -l)
          if [[ ${ldif} -gt 0 ]]
          then
             change=true
             if [[ "${ext_mess}" == "" ]]
             then
                echo "File ${woroot}.  Code has changed."
             else
                echo "File ${woroot}.  Code has changed.  ${ext_mess}"
             fi
          elif [[ "${ext_mess}" == "" ]]
          then
             change=false
          else
             change=false
             echo "File ${woroot}.  ${ext_mess}"
          fi
          #--------------------------------------------------------------------------------#


          #--------------------------------------------------------------------------------#
          #    If files have changes, report differences.                                  #
          #--------------------------------------------------------------------------------#
          if ${change}
          then
             #----- Write differences to a temporary file. --------------------------------#
             if ${ournew}
             then
                diff -uibB <(grep -vE "^\s*!" ${filetheirs})                               \
                           <(grep -vE "^\s*!" ${fileours}) > notthesame${ext}
             else
                diff -uibB <(grep -vE "^\s*!" ${fileours})                                 \
                           <(grep -vE "^\s*!" ${filetheirs}) > notthesame${ext}
             fi
             #-----------------------------------------------------------------------------#


             #----- Open editor with all files needed. ------------------------------------#
             if ${load_work} && [[ -s ${filework} ]]
             then
                ${editor} ${filework} ${filetheirs} notthesame${ext}                       \
                   1> /dev/null 2> /dev/null
                /bin/rm -f notthesame${ext}
             elif ${load_work}
             then
                ${editor} ${filetheirs} notthesame${ext} 1> /dev/null 2> /dev/null
                echo " ---> Working file is missing!"
             else
                ${editor} ${fileours} ${filetheirs} notthesame${ext}                       \
                   1> /dev/null 2> /dev/null
             fi
             #-----------------------------------------------------------------------------#
          fi
       else
          woroot=$(echo ${fileours} | sed s@"${srcours}/"@""@g)
          echo "${woroot} is exclusive to ${ours} version."
       fi #if [[ -s ${filetheirs} ]]
     done #for fileours in ${lookuptable}
     #-------------------------------------------------------------------------------------#



     #-------------------------------------------------------------------------------------#
     #    Look for files that were changed or were removed in the remote test.             #
     #-------------------------------------------------------------------------------------#
     lookuptable=$(find ${srctheirs} -name "*${ext}")
     for filetheirs in ${lookuptable}
     do
       file=$(basename ${filetheirs})
       ourpath=$(dirname ${filetheirs} | sed s@${theirs}@${ours}@g)
       fileours="${ourpath}/${file}"
       if [[ ! -s ${fileours} ]] 
       then
           woroot=$(echo ${filetheirs} | sed s@"${srctheirs}/"@""@g)
           echo "${woroot} is exclusive to ${theirs} version."
       fi #if [[ -s ${filetheirs} ]]
     done #for fileours in ${lookuptable}
     #-------------------------------------------------------------------------------------#
   done #for ext in ${exts}
   #---------------------------------------------------------------------------------------#
done #for subdir in ${subdirs}
#------------------------------------------------------------------------------------------#
