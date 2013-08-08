#!/bin/bash

#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#   ADJUST SOME VARIABLES HERE!                                                            #
#------------------------------------------------------------------------------------------#
here=`pwd`                             # Local disk
moi=`whoami`                           # User name
diskthere='/n/moorcroftfs2'            # Output directory
queue='moorcroft_6100b'                # Queue to be used
whena='01-01-2008 00:00'               # Initial time for simulation
whenz='01-01-2009 00:00'               # Final time for simulation
isfcl=5                                # 1 = LEAF-3 run, 5 = ED-2.2 run

#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#

#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#   NO NEED TO CHANGE ANYTHING BEYOND THIS POINT UNLESS YOU'RE DEVELOPING THE SCRIPT!      #
#------------------------------------------------------------------------------------------#

#----- Define the name of this simulation. ------------------------------------------------#
thissim=`basename ${here}`
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Check the directories.                                                               #
#------------------------------------------------------------------------------------------#
basehere=`basename ${here}`
dirhere=`dirname ${here}`
while [ ${basehere} != ${moi} ]
do
   basehere=`basename ${dirhere}`
   dirhere=`dirname ${dirhere}`
done
diskhere=${dirhere}
echo '-------------------------------------------------------------------------------'
echo ' - Simulation control on disk: '${diskhere}
echo ' - Output on disk:             '${diskthere}
echo '-------------------------------------------------------------------------------'
there=`echo ${here} | sed s@${diskhere}@${diskthere}@g`
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#    Make sure that the directory there exists, if not, create all parent directories      #
# needed.                                                                                  #
#------------------------------------------------------------------------------------------#
while [ ! -s ${there} ]
do
   namecheck=`basename ${there}`
   dircheck=`dirname ${there}`
   while [ ! -s ${dircheck} ] && [ ${namecheck} != '/' ]
   do
      namecheck=`basename ${dircheck}`
      dircheck=`dirname ${dircheck}`
   done
   
   if [ ${namecheck} == '/' ]
   then
      echo 'Invalid disk for variable there:'
      echo ' DISK ='${diskhere}
      exit 58
   else
      echo 'Making directory: '${dircheck}/${namecheck}
      mkdir ${dircheck}/${namecheck}
   fi
done
#------------------------------------------------------------------------------------------#



#----- Extract the time from the whena and whenz variables. -------------------------------#
montha=`echo ${whena} | awk '{print substr($1,1,2)}'`
datea=`echo ${whena}  | awk '{print substr($1,4,2)}'`
yeara=`echo ${whena}  | awk '{print substr($1,7,4)}'`
houra=`echo ${whena}  | awk '{print substr($2,1,2)}'`
minua=`echo ${whena}  | awk '{print substr($2,4,2)}'`
timea=${houra}${minua}
monthz=`echo ${whenz} | awk '{print substr($1,1,2)}'`
datez=`echo ${whenz}  | awk '{print substr($1,4,2)}'`
yearz=`echo ${whenz}  | awk '{print substr($1,7,4)}'`
hourz=`echo ${whenz}  | awk '{print substr($2,1,2)}'`
minuz=`echo ${whenz}  | awk '{print substr($2,4,2)}'`
timez=${hourz}${minuz}
#------------------------------------------------------------------------------------------#




#------ Decide which soil moisture data set to use. ---------------------------------------#
if [ ${yeara} -le 2004 ]
then
   smds='GPCP'
else
   smds='GPNR'
fi
#------------------------------------------------------------------------------------------#



#----- Decide the number of patches based on isfcl. ---------------------------------------#
if [ ${isfcl} -eq 1 ]
then
   npatch=5
   nvegpat=4
elif [ ${isfcl} -eq 5 ]
then
   npatch=2
   nvegpat=1
else
   echo ' Invalid ISFCL -> '${isfcl}'...'
   exit 1
fi
#------------------------------------------------------------------------------------------#


echo 'I am going to set up the '${thissim}' simulation...'
echo ' - Current directory: '${here}
echo ' - Output main directory: '${there}
echo ' - Initial time: '${montha}'/'${datea}'/'${yeara}' '${timea}' UTC'
echo ' - Final time:   '${monthz}'/'${datez}'/'${yearz}' '${timez}' UTC'
echo ' - Soil moisture: '${smds}


#----- Define RAMSIN. ---------------------------------------------------------------------#
RAMSIN=${here}/RAMSIN
#------------------------------------------------------------------------------------------#



#----- Replace some flags in RAMSIN. ------------------------------------------------------#
sed -i s@myoutpath@${there}@g     ${RAMSIN}
sed -i s@mysimul@${thissim}@g     ${RAMSIN}
sed -i s@myyeara@${yeara}@g       ${RAMSIN}
sed -i s@mymontha@${montha}@g     ${RAMSIN}
sed -i s@mydatea@${datea}@g       ${RAMSIN}
sed -i s@mytimea@${timea}@g       ${RAMSIN}
sed -i s@myyearz@${yearz}@g       ${RAMSIN}
sed -i s@mymonthz@${monthz}@g     ${RAMSIN}
sed -i s@mydatez@${datez}@g       ${RAMSIN}
sed -i s@mytimez@${timez}@g       ${RAMSIN}
sed -i s@myisfcl@${isfcl}@g       ${RAMSIN}
sed -i s@mysmds@${smds}@g         ${RAMSIN}
sed -i s@mynpatch@${npatch}@g     ${RAMSIN}
sed -i s@mynvegpat@${nvegpat}@g   ${RAMSIN}
#------------------------------------------------------------------------------------------#


#----- Replace some flags in orun.sh. -----------------------------------------------------#
sed -i s@thispath@${here}@g    ${here}/orun.sh
sed -i s@thisqueue@${queue}@g  ${here}/orun.sh
sed -i s@thisjob@${thissim}@g ${here}/orun.sh
#------------------------------------------------------------------------------------------#


#----- Replace some flags in callopenmpi.sh. ----------------------------------------------#
sed -i s@thispath@${here}@g    ${here}/callopenmpi.sh
sed -i s@thisqueue@${queue}@g  ${here}/callopenmpi.sh
sed -i s@thisjob@${thissim}@g ${here}/callopenmpi.sh
#------------------------------------------------------------------------------------------#


#----- Replace some flags in purge.sh -----------------------------------------------------#
sed -i s@myoutpath@${there}@g  ${here}/purge.sh
#------------------------------------------------------------------------------------------#


#----- Replace some flags in ramspost.inp -------------------------------------------------#
sed -i s@myoutpath@${there}@g ${here}/tothere/rpost/ramspost.inp
sed -i s@mysimul@${thissim}@g ${here}/tothere/rpost/ramspost.inp
#------------------------------------------------------------------------------------------#


#----- Replace some flags in srun.sh ------------------------------------------------------#
sed -i s@myoutpath@${there}@g  ${here}/tothere/rpost/srun.sh
sed -i s@thisqueue@${queue}@g  ${here}/tothere/rpost/srun.sh
sed -i s@thisjob@${thissim}@g  ${here}/tothere/rpost/srun.sh
#------------------------------------------------------------------------------------------#


#----- Replace some flags in 1eachtime-sigma.sh -------------------------------------------#
sed -i s@myoutpath@${there}@g ${here}/tothere/rpost/1eachtime-sigma.sh
#------------------------------------------------------------------------------------------#


#----- Replace some flags in last_histo.sh ------------------------------------------------#
sed -i s@thispath@${here}@g   ${here}/last_histo.sh
sed -i s@myoutpath@${there}@g ${here}/last_histo.sh
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Check whether the output directory is already "there" (the disk with large storage   #
# capacity).  This is the 
#------------------------------------------------------------------------------------------#
if [ ! -s ${here}/tothere ]
then 
   echo 'It seems you already moved tothere to '${there}'...'
elif [ ! -s ${there} ]
then
   mv tothere ${there}
elif [ ${here} != ${there} ]
then
   echo ' There is already a directory called '${there}'...'
   echo ' Do you want to delete it? [y/N]'
   read  proceed
   if [ ${proceed} != 'y' -a ${proceed} != 'Y' ]
   then
      exit
   fi

   echo '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   echo '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   echo ' '
   echo '     Look, this will REALLY delete all '${npolys}' output directories and files...'
   echo ' '
   echo '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   echo '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'

   echo 'This is PERMANENT, once they are gone, adieu, no chance to recover them!'
   echo 'Is that what you really want? [y/N]'
   read proceed

   echo ' '

   if [ ${proceed} != 'y' -a ${proceed} != 'Y' ]
   then 
      exit
   fi

   echo 'Okay then, but if you regret later do not say that I did not warn you...'
   echo 'I am giving you a few seconds to kill this script in case you change your mind...'
   delfun=16
   while [ ${delfun} -gt 1 ]
   do
      delfun=`expr ${delfun} - 1`
      echo '  - Deletion will begin in '${delfun}' seconds...'
      sleep 1
   done
   rm -frv ${there}

   mv tothere ${there}
else
   mv tothere/* .
fi
#------------------------------------------------------------------------------------------#


