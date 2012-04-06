#!/bin/sh
#------------------------------------------------------------------------------------------#
# 1eachtime.sh                                                                             #
# Developed by Marcos Longo - Lab. MASTER/IAG/USP - May 19, 2004                           #
#                                                                                          #
#   This script runs Ramspost for each file separatedly, then creates a template that read #
# all of them.  It deals with compressed files, but those files are uncompressed in a      #
# scratch folder, so it saves some time.  It also skips files that were already generated. #
#------------------------------------------------------------------------------------------#

#------------------------------------------------------------------------------------------#
#    CHANGE LOG                                                                            #
#------------------------------------------------------------------------------------------#
ramspost='myoutpath/rpost/ramspost_6.2'      # Name of executable file
tmpfolder='myoutpath/rpost/.temp'            # Name of a timestrrary folder
nice=''                                      # Command to "nice" the job.  Put nothing
                                             #    if you don't want to be nice
runoutput='myoutpath/rpost/ramspost.out'     # Name of a renewable output file
compression='none'                           # Kind of compression:
                                             #    (Z, bz2, zip, gz, or none)   
title='EDBRAMS-1.4'                          # Title to appear in the header 
                                             #    (no practical relevance)
deleteintctl='y'                             # Delete intermediate ctl [y/N]
                                             #    (a template will be provided)
outshell='y'                                 # Create an output file for shell
shellout='myoutpath/rpost/serial_out.out'    # File for 1eachtime-sigma.sh output

#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#-------------------- YOU SHOULD NOT CHANGE ANYTHING BEYOND THIS POINT --------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#

#------------------------------------------------------------------------------------------#
#      Define the kind of compression was applied to the files. The default is no          #
# compression.                                                                             #
#------------------------------------------------------------------------------------------#
case ${compression} in
  'zip' ) ending='.zip' ;unzip='unzip'      ;;
  'ZIP' ) ending='.ZIP' ;unzip='unzip'      ;;
  'bz2' ) ending='.bz2' ;unzip='bunzip2'    ;;
  'BZ2' ) ending='.BZ2' ;unzip='bunzip2'    ;;
  'gz'  ) ending='.gz'  ;unzip='gunzip'     ;;
  'GZ'  ) ending='.GZ'  ;unzip='gunzip'     ;;
  'Z'   ) ending='.Z'   ;unzip='uncompress' ;;
  'z'   ) ending='.z'   ;unzip='gunzip'     ;;
  *     ) ending=''     ;unzip='touch'      ;; # I know that this is horrible, but hey,
                                               #   it works ;-)
esac
#------------------------------------------------------------------------------------------#



#----- Clean the temporary directory. -----------------------------------------------------#
if [ ! -s ${tmpfolder} ]
then
   mkdir ${tmpfolder}
else
   rm -f ${tmpfolder}/*
fi
#------------------------------------------------------------------------------------------#


#----- Reset the shell output file in case we want one. -----------------------------------#
rm -f ${shellout}
if [ ${outshell} == 'y' -o ${outshell} == 'Y' ]
then
   touch ${shellout}
fi
#------------------------------------------------------------------------------------------#


#----- Check whether ramspost.inp exists or not. ------------------------------------------#
if [ -s ramspost.inp-backup ]
then
   rm -f ramspost.inp
   cp ramspost.inp-backup ramspost.inp
elif [ -s ramspost.inp ]
then
  cp -f ramspost.inp ramspost.inp-backup
else
  if [ ${outshell} == 'y' -o ${outshell} == 'Y' ]
  then
     echo 'There should be a file called ramspost.inp here. Exitting...' >> ${shellout}
  else
     echo 'There should be a file called ramspost.inp here. Exitting...'
  fi
  exit
fi
#------------------------------------------------------------------------------------------#



#----- Determine the analysis prefix from the list. ---------------------------------------#
fprefix=`grep -i FPREFIX ramspost.inp | grep -vi "\-\-"`
fprefix=`echo ${fprefix} | sed s/" "/""/g |sed s/"'"/""/g`
ext=`echo ${fprefix} |wc -c`
p=0

#----- Find the comma position. -----------------------------------------------------------#
while [ ${p} -lt ${ext} ]  # Finding the comma position
do
  p=`expr ${p} + 1`
  comma=`echo "${fprefix}" | awk '{print substr($1,'${p}',1)}'`
  if [ y${comma} == 'y,' ]
  then 
     comma=${p}
     p=${ext}
  fi
done
ext=`expr ${comma} - 9` # 9 is the position after the comma
fprefix=`echo ${fprefix} | awk '{print substr($1,9,'${ext}')}'`
#------------------------------------------------------------------------------------------#

#------------------------------------------------------------------------------------------#
# Determine the output file prefix from the namelist.                                      #
#------------------------------------------------------------------------------------------#
gprefix=`grep -i GPREFIX ramspost.inp | grep -vi "\-\-"`
gprefix=`echo ${gprefix} | sed s/" "/""/g |sed s/"'"/""/g`
ext=`echo ${gprefix} |wc -c`
p=0
#----- Find the comma position. -----------------------------------------------------------#
while [ ${p} -lt ${ext} ]                           
do
  p=`expr ${p} + 1`
  comma=`echo "${gprefix}" | awk '{print substr($1,'${p}',1)}'`
  if [ y${comma} == 'y,' ]
  then 
     comma=${p}
     p=${ext}
  fi
done
ext=`expr ${comma} - 9` #----- 9 corresponds to the position after the comma
gprefix=`echo ${gprefix} | awk '{print substr($1,9,'${ext}')}'`
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#     Run for each time, provided that it hasn't been run yet.                             #
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
donecount=0
#----- The list mylist is built in a way that doesn't crash when the list is too long. ----#
directory=`dirname ${fprefix}`
myprefix=`basename ${fprefix}`
mylist=`ls -1 ${directory}/ |grep ${myprefix} | grep head.txt${ending}`
#----- Loop for each file. ----------------------------------------------------------------#
for analysis in ${mylist}
do
  #----- Find the corresponding vfm file. -------------------------------------------------#
  analysis=`basename ${analysis}`
  vfm=`basename ${analysis} head.txt${ending}`'g*.vfm'${ending}
  #----- Rearrange the the prefix of both analysis and output for a given time... ---------#
  anapref=${tmpfolder}'/'`basename ${analysis} \-head.txt${ending}`
  ext=`echo ${anapref} | wc -c`;p0=`expr ${ext} - 17`
  timestr=`echo ${anapref} | awk '{print substr($1,'${p0}',15)}'`
  outpref=${gprefix}${timestr}
  #----- Here we check whether the files already exist. -----------------------------------#
  if [ ! -s ${outpref}'_g1.gra' ]
  then
    #--------------------------------------------------------------------------------------#
    #     In case they don't, we copy the header and data files into a temporary folder    #
    # and uncompress them...                                                               #
    #--------------------------------------------------------------------------------------#
    cp ${directory}'/'${analysis} ${directory}'/'${vfm} ${tmpfolder}
    ${nice} ${unzip} ${tmpfolder}'/'`basename ${analysis}` ${tmpfolder}'/'`basename ${vfm}`
    donecount=`expr ${donecount} + 1`
    if [ ${outshell} == 'y' -o ${outshell} == 'Y' ]
    then
       echo 'Running ramspost for time '${timestr}'...' >> ${shellout}
    else
       echo 'Running ramspost for time '${timestr}'...'
    fi

    #--------------------------------------------------------------------------------------#
    #     Switch the fprefix by the timestrrary folder, with the added restriction to the  #
    # time.                                                                                #
    #--------------------------------------------------------------------------------------#
    rm -f ramspost.inp
    cp ramspost.inp-backup ramspost.inp
    
    sed -i s@${fprefix}@${anapref}@g ramspost.inp
    sed -i s@${gprefix}@${outpref}@g ramspost.inp

    ${nice} ${ramspost} > ${runoutput}

    #----- Clean the temporary folder, so I save disk space. ------------------------------#
    rm -f ${tmpfolder}/*
  else
    #--------------------------------------------------------------------------------------#
    #      If the files exist, I skip this time                                            #
    #--------------------------------------------------------------------------------------#
    skipmess='I will not run for time '${timestr}' because these files already exist...'
    if [ ${outshell} == 'y' -o ${outshell} == 'Y' ]
    then
       echo ${skipmess} >> ${shellout}
    else
       echo ${skipmess} >> ${shellout}
    fi
  fi
done
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#      Determine the number of grids.                                                      #
#------------------------------------------------------------------------------------------#
gmax=`ls -1 ${fprefix}${timestr}* |wc -l`
gmax=`expr ${gmax} - 1`
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#

#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#      Change the output folder to avoid possible "./".                                    #
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
ext=`echo ${gprefix} | wc -c`
startpnt=`echo ${gprefix} | awk '{print substr($1,1,2)}'`

rad=${gprefix}
if [ 'z'${startpnt} == 'z./' ] 
then
  ext=`expr ${ext} - 2`
  gprefix=`echo ${gprefix} |awk '{print substr($1,3,'${ext}')}'`
fi
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#      Find the time interval in sec.                                                      #
#------------------------------------------------------------------------------------------#
txtpath=`dirname ${fprefix}`
txtbase=`basename ${fprefix}`
allheads=`ls -1 ${txtpath} | grep ${txtbase} | grep txt${ending}`
first=${txtpath}/`echo ${allheads}| awk '{print $1}'`
second=${txtpath}/`echo ${allheads}| awk '{print $2}'`
if [ ${second} == ${txtpath}/ ]
then
   #----- Only one file exists, assign 1hr for time step. ---------------------------------#
   deltat='1hr'
else
   #---------------------------------------------------------------------------------------#
   #      Uncompress into a scratch file.  Though it may look inefficient, it's just a txt #
   # file...                                                                               #
   #---------------------------------------------------------------------------------------#
   cp ${first} './deleteme.txt'${ending}
   ${unzip} './deleteme.txt'${ending}
   first='./deleteme.txt'
   lmax=`cat ${first} | wc -l`
   l=0
   #----- Get the time of the first analysis ----------------------------------------------#
   while [ ${l} -lt ${lmax} ]
   do 
     l=`expr ${l} + 1`
     thisline=`head -${l} ${first} | tail -1`
     whichone=`echo ${thisline} | awk '{print $1}'`
     if [ z$whichone = 'z__time' ]
     then 
       l=`expr ${l} + 2`
       thisline=`head -${l} ${first} | tail -1`
       first=`echo ${thisline} | awk '{print $1}'`
       l=${lmax}
     fi
   done
   rm -f ./deleteme.txt
   #----- Do the same thing for the second file... ----------------------------------------#
   cp ${second} './deleteme.txt'${ending}
   #---------------------------------------------------------------------------------------#
   #      Uncompress into a scratch file.  Though it may look inefficient, it's just a txt #
   # file...                                                                               #
   #---------------------------------------------------------------------------------------#
   ${unzip} './deleteme.txt'${ending}
   second='./deleteme.txt'
   #----- Determine how many lines are in the header file ---------------------------------#
   lmax=`cat ${second} | wc -l`
   l=0
   #----- Get the time of the second analysis. --------------------------------------------#
   while [ ${l} -lt ${lmax} ];
   do 
     l=`expr ${l} + 1`
     thisline=`head -${l} ${second} | tail -1`
     whichone=`echo ${thisline} | awk '{print $1}'`
     if [ z$whichone = 'z__time' ]
     then
       l=`expr ${l} + 2`
       thisline=`head -${l} ${second} | tail -1`
       second=`echo ${thisline} | awk '{print $1}'`
       l=${lmax}
     fi
   done
   rm -f ./deleteme.txt

   #---------------------------------------------------------------------------------------#
   #     Find the interval between the first and the second time.  We must switch E+ by    #
   # 10^, so bc will work.  The time interval must be given in minutes or hours in the CTL #
   # file, as GrADS doesn't accept seconds as units.                                       #
   #---------------------------------------------------------------------------------------#
   first=`echo ${first} | sed s@'E+'@'*(10^'@g`')'
   second=`echo ${second} | sed s@'E+'@'*(10^'@g`')'
   deltat=`echo ${second}' - '${first}| bc`
   #----- Delta-t is the time interval in seconds. ----------------------------------------#
   deltat=`echo ${deltat} \/ 60 |bc`
   #----- Decide whether to show the interval in minutes or hours. ------------------------#      
   deltathour=`echo ${deltat} \/ 60 |bc`
   wouldbeinteger=`echo ${deltathour} \* 60 |bc`
   if [ ${deltat} -eq ${wouldbeinteger} ]
   then
      deltat=${deltathour}'hr'
   else 
      deltat=${deltat}'mn'
   fi
   #------------------------------------------------------------------------------------------#
   #------------------------------------------------------------------------------------------#
fi
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#     In case binary files were created, we generate a new template for each grid.         #
#------------------------------------------------------------------------------------------#
if [ ${donecount} -gt 0 ]
then
  g=0
  while [ ${g} -lt ${gmax} ]
  do
    g=`expr ${g} + 1`
    ctl=${rad}'_g'${g}'.ctl'
    ctl=`basename ${ctl}`

    binpath=`dirname ${gprefix}`
    basepath=`basename ${gprefix}`

    #--------------------------------------------------------------------------------------#
    #      Take the first file as a "canvas"...  To make sure we can find the ctl even     #
    # when there are too many files, we list the directory, and grep the files.            #
    #--------------------------------------------------------------------------------------#
    ref=`ls -1 ${binpath} | grep ${basepath} | grep ${g}'.ctl'`
    ref=${binpath}/`echo ${ref} | awk '{print $1}'`
    lmax=`cat ${ref} | wc -l`

    #--------------------------------------------------------------------------------------#
    #     Go through the file line by line and check whether the line must be modified for #
    # a good template ctl. If not, we just paste the line as is into the template.         #
    #--------------------------------------------------------------------------------------#
    l=0
    while [ ${l} -lt ${lmax} ]
    do
      l=`expr ${l} + 1`
      thisline=`head -${l} ${ref} | tail -1`
      whichone=`echo ${thisline} |awk '{print $1}'`

      if [ z${whichone} == 'zdset' ]
      then         
        #----- Filename line. Switch it and also add the "template" line ------------------#
        if [ 'z'`echo ${gprefix} | awk '{print substr($1,1,1)}'` == 'z/' ]
        then 
           echo 'dset '${gprefix}'%y4-%m2-%d2-%h2%n2_g'${g}'.gra' > ${ctl}
        else 
           echo 'dset ^'${gprefix}'%y4-%m2-%d2-%h2%n2_g'${g}'.gra' > ${ctl}
        fi
        echo 'options template' >> ${ctl}

      elif [ z${whichone} == 'ztitle' ] 
      then
        #------ Not an important change, it just changes the title of this dataset. -------#
        echo 'title '${title} >> ${ctl}
      elif [ z${whichone} == 'ztdef' ]
      then
        #----- Time line, switch it by the total number of times. -------------------------#
        amax=`ls -1 ${gprefix}*${g}.gra |wc -l`
        timestra=`echo ${thisline} | awk '{print $4}'`
        echo 'tdef '${amax}' linear '${timestra}' '${deltat} >> ${ctl}
      else
        #----- No need to change, simply copy it to the template. -------------------------#
        echo ${thisline} >> ${ctl}
      fi
    done  

    #--------------------------------------------------------------------------------------#
    #     Display the command line in the screen so we know when the shell script is done, #
    # and also allow lazy people to save time by copying the command with the mouse ;-).   #
    #--------------------------------------------------------------------------------------#
    if [ ${outshell} == 'y' -o ${outshell} == 'Y' ]
    then
       echo "grads -cl 'open ${ctl}'" >> ${shellout}
    else
       echo "grads -cl 'open ${ctl}'"
    fi
  done 
fi

#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#      Final steps, cleaning up the CTL files that are not required for the next run...    #
# We always leave the first CTL file so the next time the script is called we have the     #
# best file to start when we re-create the template.                                       #
#------------------------------------------------------------------------------------------#
if [ ${deleteintclt}='y' -o ${deleteintctl}='Y' ]
then 
  g=0
  while [ ${g} -lt ${gmax} ]
  do
     g=`expr ${g} + 1`

     binpath=`dirname ${gprefix}`
     basepath=`basename ${gprefix}`

     howmany=`ls ${binpath} | grep ${basepath} | grep ${g}'.ctl' 2> /dev/null | wc -l`
     if [ ${howmany} -ge 1 ]
     then
        files=`ls ${binpath} | grep ${basepath} | grep ${g}'.ctl'`
        n=1
        while [ ${n} -lt ${howmany} ]
        do
          let n=${n}+1
          file=${binpath}/`echo ${files} | awk '{print $'${n}'}'`
          rm -f ${file}
        done #file in ${files}
     fi # [ ${howmany} -ge 1 ]
  done #while [ ${g} -lt ${gmax} ]
fi #if [ ${deleteintclt}='y' -o ${deleteintctl}='Y' ]
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#      Return ramspost.inp to the original file                                            #
#------------------------------------------------------------------------------------------#
rm -f ramspost.inp
mv -f ramspost.inp-backup ramspost.inp 
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
