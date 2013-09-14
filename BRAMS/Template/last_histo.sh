#!/bin/bash
. ${HOME}/.bashrc
here="thispath"             # ! Main path
there="myoutpath"           # ! Disk where the output files are
bzip2="/bin/gzip -9"        # ! Program to compress files (with options)
#------ Calculator. -----------------------------------------------------------------------#
ccc="${HOME}/util/calc.sh"  # Calculator
#------ # of history files to keep (in odyssey's world, always keep more than one). -------#
retain=3
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#     Set the maximum number of days for each month.  day[0] is a dummy variable.          #
#------------------------------------------------------------------------------------------#
daymax=(0 31 28 31 30 31 30 31 31 30 31 30 31)
#------------------------------------------------------------------------------------------#


#----- Check whether run_sitter.sh is still running or not.  If it is, exit. --------------#
if [ -s ${here}/last_histo.lock ]
then
   exit
else
   echo 'I am going to clean history files. Lots of them!' > ${here}/last_histo.lock
fi
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#       Set up RAMSIN file.                                                                #
#------------------------------------------------------------------------------------------#
RAMSIN="${here}/RAMSIN"
polyname=`basename ${here}`
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#     Find out the first and last time.                                                    #
#------------------------------------------------------------------------------------------#
montha=`cat ${RAMSIN} | grep IMONTHA | grep "=" | sed s@,@@g | awk '{print $3}'`
datea=`cat  ${RAMSIN} | grep IDATEA  | grep "=" | sed s@,@@g | awk '{print $3}'`
yeara=`cat  ${RAMSIN} | grep IYEARA  | grep "=" | sed s@,@@g | awk '{print $3}'`
timea=`cat  ${RAMSIN} | grep ITIMEA  | grep "=" | sed s@,@@g | awk '{print $3}'`
monthz=`cat ${RAMSIN} | grep IMONTHZ | grep "=" | sed s@,@@g | awk '{print $3}'`
datez=`cat  ${RAMSIN} | grep IDATEZ  | grep "=" | sed s@,@@g | awk '{print $3}'`
yearz=`cat  ${RAMSIN} | grep IYEARZ  | grep "=" | sed s@,@@g | awk '{print $3}'`
timez=`cat  ${RAMSIN} | grep ITIMEZ  | grep "=" | sed s@,@@g | awk '{print $3}'`
#----- Find time and minute. --------------------------------------------------------------#
houra=`echo ${timea}  | awk '{print substr($1,1,2)}'`
minua=`echo ${timea}  | awk '{print substr($1,3,2)}'`
hourz=`echo ${timez}  | awk '{print substr($1,1,2)}'`
minuz=`echo ${timez}  | awk '{print substr($1,3,2)}'`
#----- Find the return period of analysis files. ------------------------------------------#
dtime=`cat  ${RAMSIN} | grep FRQANL  | grep "=" | sed s@,@@g | awk '{print $3}'`
dtime=`echo ${dtime} | sed s@"\\."@@g`
let dtime=${dtime}/3600
#------------------------------------------------------------------------------------------#



echo " - Deleting history files for simulation ${polyname}:"

#---- First we delete all -Z- and -R- files. ----------------------------------------------#
nzed=`/bin/ls -1    ${there}/ecort/${polyname}-Z-*h5       2> /dev/null | wc -l`
narhead=`/bin/ls -1 ${there}/histo/${polyname}-R-*head.txt 2> /dev/null | wc -l`
narvfm=`/bin/ls -1  ${there}/histo/${polyname}-R-*.vfm     2> /dev/null | wc -l`
if [ ${nzed} -gt 0 ]
then
   zeds=`/bin/ls -1 ${there}/ecort/${polyname}-Z-*h5 2> /dev/null`
   for zed in ${zeds}
   do
      echo -n "    - Deleting: `basename ${zed}`..."
      /bin/nice /bin/rm -f ${zed}
      echo " Gone!"
   done
fi
if [ ${narhead} -gt 0 ]
then
   arheads=`/bin/ls -1 ${there}/histo/${polyname}-R-*head.txt 2> /dev/null`
   for arhead in ${arheads}
   do
      echo -n "    - Deleting: `basename ${arhead}`..."
      /bin/nice /bin/rm -f ${arhead}
      echo " Gone!"
   done
fi
if [ ${narvfm} -gt 0 ]
then
   arvfms=`/bin/ls -1 ${there}/histo/${polyname}-R-*.vfm 2> /dev/null`
   for arvfm in ${arvfms}
   do
      echo -n "    - Deleting: `basename ${arvfm}`..."
      /bin/nice /bin/rm -f ${arvfm}
      echo " Gone!"
   done
fi
#------------------------------------------------------------------------------------------#




#---- Now we delete all -S- and -H- files except the last one. ----------------------------#
nessh5=`/bin/ls     -1 ${there}/ecort/${polyname}-S-*h5        2> /dev/null | wc -l`
nesscmp=`/bin/ls    -1 ${there}/ecort/${polyname}-S-*cmp       2> /dev/null | wc -l`
naitchhead=`/bin/ls -1 ${there}/histo/${polyname}-H-*head.txt  2> /dev/null | wc -l`
naitchvfm=`/bin/ls  -1 ${there}/histo/${polyname}-H-*.vfm      2> /dev/null | wc -l`
if [ ${nessh5} -gt ${retain} ]
then
   let head=${nessh5}-${retain}
   esses=`/bin/ls -1 ${there}/ecort/${polyname}-S-*h5 | head -${head}`
   for ess in ${esses}
   do
      echo -n "    - Deleting: `basename ${ess}`..."
      /bin/nice /bin/rm -f ${ess}
      echo " Gone!"
   done
fi
if [ ${nesscmp} -gt ${retain} ]
then
   let head=${nesscmp}-${retain}
   esses=`/bin/ls -1 ${there}/ecort/${polyname}-S-*cmp | head -${head}`
   for ess in ${esses}
   do
      echo -n "    - Deleting: `basename ${ess}`..."
      /bin/nice /bin/rm -f ${ess}
      echo " Gone!"
   done
fi
if [ ${naitchhead} -gt ${retain} ]
then
   let head=${naitchhead}-${retain}
   aitches=`/bin/ls -1 ${there}/histo/${polyname}-H-*head.txt | head -${head}`
   for aitch in ${aitches}
   do
      echo -n "    - Deleting: `basename ${aitch}`..."
      /bin/nice /bin/rm -f ${aitch}
      echo " Gone!"
   done
fi
if [ ${naitchvfm} -gt ${retain} ]
then
   let head=${naitchvfm}-${retain}
   aitches=`/bin/ls -1 ${there}/histo/${polyname}-H-*.vfm | head -${head}`
   for aitch in ${aitches}
   do
      echo -n "    - Deleting: `basename ${aitch}`..."
      /bin/nice /bin/rm -f ${aitch}
      echo " Gone!"
   done
fi
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#      Check whether we can zip files.                                                     #
#------------------------------------------------------------------------------------------#
analy="${there}/analy"
binpref="${there}/rpost/binary/${polyname}"
nbin=`ls ${binpref}*.gra 2> /dev/null | wc -l`
if [ ${nbin} -gt 1 ]
then
   lastbin=`ls ${binpref}*.gra | tail -1 | head -1`
   lastbin=`echo ${lastbin} | sed s@${binpref}@@g | sed s@_g1.gra@@g`
   yeare=`echo  ${lastbin} | awk '{print substr($1,1,4)}'`
   monthe=`echo ${lastbin} | awk '{print substr($1,6,2)}'`
   datee=`echo  ${lastbin} | awk '{print substr($1,9,2)}'`
   houre=`echo  ${lastbin} | awk '{print substr($1,12,2)}'`
   minue=`echo  ${lastbin} | awk '{print substr($1,14,2)}'`

   let year=${yeara}-1


   #---------------------------------------------------------------------------------------#
   #     Loop over the years.                                                              #
   #---------------------------------------------------------------------------------------#
   while [ ${year} -lt ${yeare} ]
   do
      let year=${year}+1
      yyyy=`printf "%2.2i" ${year}`

      #----- Update daymax for February depending on whether it is a leap year or not. ----#
      let leap400=${year}%400
      let leap100=${year}%100
      let leap004=${year}%4
      if [ ${leap400} -eq 0 ] || [ ${leap100} -ne 0 -a ${leap004} -eq 0 ]
      then
         daymax[2]=29
      else
         daymax[2]=28
      fi
      #------------------------------------------------------------------------------------#

      if [ ${year} -eq ${yeare} ]
      then
         monthl=${monthe}
      else
         monthl=12
      fi
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Loop over the months.                                                          #
      #------------------------------------------------------------------------------------#
      month=0
      while [ ${month} -lt ${monthl} ]
      do
         let month=${month}+1
         mm=`printf "%2.2i" ${month}`


         if [ ${year} -eq ${yeare} ] && [ ${month} -eq ${monthe} ]
         then
            datel=${datee}
         else
            datel=${daymax[${month}]}
         fi
         #---------------------------------------------------------------------------------#


         date=0
         while [ ${date} -lt ${datel} ]
         do
            let date=${date}+1
            dd=`printf "%2.2i" ${date}`

            if [ ${year} -eq ${yeare} ] && [ ${month} -eq ${monthe} ] && 
               [ ${date} -eq ${datee} ]
            then
               hourl=${houre}
            else
               let hourl=24-${dtime}
            fi

            hour=-${dtime}
            while [ ${hour} -lt ${hourl} ]
            do
               let hour=${hour}+${dtime}
               hh=`printf "%2.2i" ${hour}`

               ahead=${analy}/${polyname}-A-${yyyy}-${mm}-${dd}-${hh}0000-head.txt
               avfm=${analy}/${polyname}-A-${yyyy}-${mm}-${dd}-${hh}0000-g1.vfm

               if [ -s ${ahead} ]
               then
                  echo -n "   - Squeezing file: `basename ${ahead}`..."
                  ${bzip2} ${ahead} 2> /dev/null
                  echo "Zipped!"
               fi

               if [ -s ${avfm} ]
               then
                  echo -n "   - Squeezing file: `basename ${avfm}`..."
                  ${bzip2} ${avfm} 2> /dev/null
                  echo "Zipped!"
               fi
               #---------------------------------------------------------------------------#
            done
            #------------------------------------------------------------------------------#
         done
         #---------------------------------------------------------------------------------#
      done
      #------------------------------------------------------------------------------------#
   done
   #---------------------------------------------------------------------------------------#
fi
#------------------------------------------------------------------------------------------#

/bin/rm -f ${here}/last_histo.lock
