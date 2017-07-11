#!/bin/sh
here=`pwd`
datum=myoutpath

#----- Check that the user is aware that it will remove everything... ---------------------#
if [ 'x'${1} == 'x-a' ] || [ 'x'${1} == 'xall' ]
then
   echo '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   echo ' '
   echo 'This script will remove files and sub-directories from the following directory:'
   echo ' '
   echo ${datum}
   echo ' '
   echo 'Are you sure that you want to remove all files and directories? [y/N]'
   echo ' '
   echo '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
else
   echo '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   echo ' '
   echo 'This script will remove files from the following directory:'
   echo ' '
   echo ${datum}
   echo ' '
   echo 'Are you sure that you want to remove all files? [y/N]'
   echo ' '
   echo '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
fi
read proceed
if [ ${proceed} != 'y' -a ${proceed} != 'Y' ]
then
   exit
fi
#----- Check that the user is aware that it will remove everything... ---------------------#
echo ' '
if [ 'x'${1} == 'x-d' ]
then
   echo '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   echo '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   echo '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   echo '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   echo '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   echo '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   echo '!!! '
   echo '!!!     Look, this will REALLY delete all output directories and files from'
   echo '!!! '
   echo '!!!     '${datum}
   echo '!!! '
   echo '!!!     This is PERMANENT, once they are gone, adieu, no chance to recover them!'
   echo '!!!     Is that what you really want? [y/N]'
   echo '!!! '
   echo '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   echo '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   echo '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   echo '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   echo '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   echo '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
else
   echo '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   echo '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   echo '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   echo '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   echo '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   echo '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   echo '!!! '
   echo '!!!     Look, this will REALLY delete all '${npolys}' output files from:'
   echo '!!! '
   echo '!!!     '${datum}
   echo '!!! '
   echo '!!!     This is PERMANENT, once they are gone, adieu, no chance to recover them!'
   echo '!!!     Is that what you really want? [y/N]'
   echo '!!! '
   echo '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   echo '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   echo '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   echo '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   echo '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   echo '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
fi

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


if [ 'x'${1} == 'x-a' ] || [ 'x'${1} == 'xall' ]
then
   rm -frv ${datum}/analy
   rm -frv ${datum}/histo
   rm -frv ${datum}/isean
   rm -frv ${datum}/surfa
   rm -frv ${datum}/micro
   rm -frv ${datum}/ecort
   rm -frv ${datum}/ecoss
   mkdir   ${datum}/analy
   mkdir   ${datum}/histo
   mkdir   ${datum}/isean
   mkdir   ${datum}/surfa
   mkdir   ${datum}/micro
   mkdir   ${datum}/ecort
   mkdir   ${datum}/ecoss
   rm -fv core.* ???_out.out ???_out.err ???_lsf.out fort.*
else
   rm -frv ${datum}/analy
   rm -frv ${datum}/histo
   rm -frv ${datum}/micro
   rm -frv ${datum}/ecoss
   rm -frv ${datum}/ecort
   mkdir   ${datum}/analy
   mkdir   ${datum}/histo
   mkdir   ${datum}/micro
   mkdir   ${datum}/ecoss
   mkdir   ${datum}/ecort
   rm -fv core.* ???_out.out ???_out.err ???_lsf.out fort.*
fi
