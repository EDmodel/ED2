#!/bin/bash
. ${HOME}/.bashrc
here=`pwd`
thisqueue='moorcroft2b'               # ! Queue where jobs should be submitted
#----- Check whether to use openlava or typical job submission. ---------------------------#
openlava='y'
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#    Grab arguments and also set some defaults.                                            #
#------------------------------------------------------------------------------------------#
nargs=$#               # Number of arguments
args=$@                # List of arguments
ibg="0"                # Background is white by default
rdata="RData_scenario"
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Check whether this machine supports openlava or not (in case this is a openlava.     #
#------------------------------------------------------------------------------------------#
thismach=`hostname -s`
if [ ${thismach} == "rclogin01" ] && [ ${openlava} == "y" -o ${openlava} == "Y" ]
then
   . /opt/openlava-2.0/etc/openlava-client.sh
elif [ ${thismach} != "rclogin01" ] && [ ${openlava} == "y" -o ${openlava} == "Y" ]
then
   echo " This machine (${hostname}) doesn't support openlava, sorry!"
   exit 39
fi
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#     Loop over all arguments, and assign soil texture and temperature based on the        #
# options.                                                                                 #
#------------------------------------------------------------------------------------------#
istext=16
itemp=0
irain=0
if [ ${nargs} -gt 0 ]
then
   nextstext="n"
   nexttemp="n"
   nextrain="n"
   nextbg="n"
   n=0
   for arg in ${args}
   do
      let n=${n}+1

      if [ ${nextstext} == "y" ]
      then
         istext=${arg}
         nextstext="n"
      elif [ ${nexttemp} == "y" ]
      then
         itemp=${arg}
         nexttemp="n"
      elif [ ${nextrain} == "y" ]
      then
         irain=${arg}
         nextrain="n"
      elif [ ${nextbg} == "y" ]
      then
         ibg=${arg}
         nextbg="n"
      elif [ "x${arg}" == "x-s" ]
      then
         nextstext="y"
      elif [ "x${arg}" == "x-r" ]
      then
         nextrain="y"
      elif [ "x${arg}" == "x-t" ]
      then
         nexttemp="y"
      elif [ "x${arg}" == "x-b" ]
      then
         nextbg="y"
      else
         echo "Not sure what you mean..."
         echo "Argument ${n}  (${arg}) doesn't make sense!"
         exit 91
      fi
   done
fi
#------------------------------------------------------------------------------------------#


bad=0


#------------------------------------------------------------------------------------------#
#     Check that soil texture makes sense.                                                 #
#------------------------------------------------------------------------------------------#
case ${istext} in
2|stext02)
   stext="stext02"
   ;;
6|stext06)
   stext="stext06"
   ;;
8|stext08)
   stext="stext08"
   ;;
11|stext11)
   stext="stext11"
   ;;
16|stext16)
   stext="stext16"
   ;;
*)
   echo "Invalid 'stext' argument (${istext})"
   let bad=${bad}+1
   ;;
esac
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#     Check that temperature makes sense.                                                  #
#------------------------------------------------------------------------------------------#
case ${itemp} in
0|"t+000")
   temp="t+000"
   idxtemp="1"
   ;;
1|"t+100")
   temp="t+100"
   idxtemp="2"
   ;;
2|"t+200")
   temp="t+200"
   idxtemp="3"
   ;;
3|"t+300")
   temp="t+300"
   idxtemp="4"
   ;;
*)
   echo "Invalid 'temp' argument (${itemp})"
   let bad=${bad}+1
   ;;
esac
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#     Check that rainfall makes sense.                                                     #
#------------------------------------------------------------------------------------------#
case ${irain} in
0|00|"r+000")
   rain="r+000"
   ;;
2|02|"r-020")
   rain="r-020"
   ;;
4|04|"r-040")
   rain="r-040"
   ;;
6|06|"r-060")
   rain="r-060"
   ;;
8|08|"r-080")
   rain="r-080"
   ;;
10|"r-100")
   rain="r-100"
   ;;
12|"r-120")
   rain="r-120"
   ;;
14|"r-140")
   rain="r-140"
   ;;
16|"r-160")
   rain="r-160"
   ;;
*)
   echo "Invalid 'rain' argument (${irain})"
   let bad=${bad}+1
   ;;
esac
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#     Check that ibackground makes sense.                                                  #
#------------------------------------------------------------------------------------------#

case ${ibg} in
0)
   bg="ibg00"
   ;;
1)
   bg="ibg01"
   ;;
2)
   bg="ibg02"
   ;;
*)
   echo "Invalid 'bg' argument (${ibg})"
   let bad=${bad}+1
   ;;
esac
#------------------------------------------------------------------------------------------#





#------------------------------------------------------------------------------------------#
#     Stop in case anything is wrongly set.                                                #
#------------------------------------------------------------------------------------------#
if [ ${bad} -gt 0 ]
then 
   exit 99
fi
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#      Create a unique job name.                                                           #
#------------------------------------------------------------------------------------------#
job="cscen_${stext}_${rain}_${temp}_${bg}"
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#      Define some R-related settings.                                                     #
#------------------------------------------------------------------------------------------#
R_Comm="R CMD BATCH --no-save --no-restore"
R_Out="${here}/.${job}_out.out"
R_Orig="${here}/compare_scenarios.r"
R_Submit="${here}/.${job}.r"
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#     Switch settings to the current configuration.                                        #
#------------------------------------------------------------------------------------------#
/bin/cp ${R_Orig} ${R_Submit}
stext_ostr=`grep stext.default ${R_Submit} | head -1`
temp_ostr=`grep use.global ${R_Submit} | head -1`
rain_ostr=`grep drain.default ${R_Submit} | head -1`
bg_ostr=`grep ibackground ${R_Submit} | head -1`
rdata_ostr=`grep rdata.path ${R_Submit} | head -1`

stext_nstr="stext.default = \"${stext}\"                   # Default soil texture"
temp_nstr="use.global       = ${idxtemp}   # Which global to use (TRUE means all of them)"
rain_nstr="drain.default = \"${rain}\"                   # Default rainfall scenario"
bg_nstr="ibackground   = ${ibg}                           # Target background colour:"
rdata_nstr="rdata.path       = file.path(here,\"${rdata}\") # Path with R object."
sed -i s@"${stext_ostr}"@"${stext_nstr}"@g ${R_Submit}
sed -i s@"${rain_ostr}"@"${rain_nstr}"@g   ${R_Submit}
sed -i s@"${temp_ostr}"@"${temp_nstr}"@g   ${R_Submit}
sed -i s@"${bg_ostr}"@"${bg_nstr}"@g       ${R_Submit}
sed -i s@"${rdata_ostr}"@"${rdata_nstr}"@g ${R_Submit}
#------------------------------------------------------------------------------------------#





#------------------------------------------------------------------------------------------#
#      Create the script to submit to LSF.                                                 #
#------------------------------------------------------------------------------------------#
cspost="${here}/.${job}.sh"
csout="${here}/.${job}_lsf.out"
cscomm="${R_Comm} ${R_Submit} ${R_Out}"
status="${here}/${rdata}/status_stext_sim_${temp}.txt"
echo "#!/bin/bash"              >  ${cspost}
echo "/bin/rm -fr ${status}"    >> ${cspost}
echo "while [ ! -s ${status} ]" >> ${cspost}
echo "do"                       >> ${cspost}
echo "   sleep 3"               >> ${cspost}
echo "   ${cscomm}"             >> ${cspost}
echo "done"                     >> ${cspost}
chmod +x ${cspost}
#------------------------------------------------------------------------------------------#




#------ Check whether to use openlava or LSF. ---------------------------------------------#
if [ 'x'${openlava} == 'xy' ] || [ 'x'${openlava} == 'xY' ]
then
   bcall="iobsub -J ${job} -o ${csout}"
else
   bcall="bsub -q ${thisqueue} -J ${job} -o ${csout}"
fi
bsub="${bcall} ${cspost} 1> /dev/null 2> /dev/null"
${bsub}
#------------------------------------------------------------------------------------------#
