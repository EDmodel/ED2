#!/bin/sh
##########################################################################################################
# 1eachtime.sh                                                                                           #
# Developed by Marcos Longo - Lab. MASTER/IAG/USP - May 19, 2004                                         #
#                                                                                                        #
#   This script runs Ramspost for each file separatedly, then creates a template that read all of them.  #
# It deals with compressed files, but those files are uncompressed in a scratch folder, so it saves some #
# time. It also skips files that were already generated.                                                 #
##########################################################################################################

##########################################################################################################
#######                                         CHANGE LOG                                         #######
##########################################################################################################
ramspost='./ramspost_6.2'                                 # Name of executable file
tmpfolder=${HOME}'/.temp'                                 # Name of a temporary folder
nice=''                    # Command to "nice" the job. Put nothing if you don't want to be nice
runoutput='ramspost.out'                                  # Name of a renewable output file
compression='none'                                        # Kind of compression:(Z, bz2, zip, gz, or none)   
title='My EDBRAMS simulation'                             # Title to appear in the header 
                                                          #   (no practical relevance)
deleteintctl='y'                                          # Delete intermediate ctl [y/N]
                                                          # (a template will be provided)

##########################################################################################################
##########################################################################################################
############################ YOU SHOULD NOT CHANGE ANYTHING BEYOND THIS POINT ############################
##########################################################################################################
##########################################################################################################

##########################################################################################################
# Defining the kind of compression was applied to the files. The default is no compression               #
##########################################################################################################
case ${compression} in
  'zip' ) ending='.zip' ;unzip='unzip'      ;;
  'ZIP' ) ending='.ZIP' ;unzip='unzip'      ;;
  'bz2' ) ending='.bz2' ;unzip='bunzip2'    ;;
  'BZ2' ) ending='.BZ2' ;unzip='bunzip2'    ;;
  'gz'  ) ending='.gz'  ;unzip='gunzip'     ;;
  'GZ'  ) ending='.GZ'  ;unzip='gunzip'     ;;
  'Z'   ) ending='.Z'   ;unzip='uncompress' ;;
  'z'   ) ending='.z'   ;unzip='gunzip'     ;;
  *     ) ending=''     ;unzip='touch'      ;; # I know that this is dirty, but it works ;-)
esac

##########################################################################################################
# Checking whether the folder ${tmpfolder} is empty or not.                                              #
##########################################################################################################
#if ! test -s ${tmpfolder}
#then
#  mkdir ${tmpfolder}
#else
#  echo -n 'The folder '${tmpfolder}' is not empty. Can I remove all its contents and proceed [y/N]? '
#  read response
#  if [ ${response} = 'y' -o ${response} = 'Y' ]
#  then
    rm -f ${tmpfolder}/*
#  else
#    echo 'The script will not be executed...'
#    exit
#  fi
#fi

##########################################################################################################
# Checking whether ramspost.inp exists or not.                                                           #
##########################################################################################################
if test -s ramspost.inp
then
  cp -f ramspost.inp ramspost.inp-salvaguarda
else
  echo 'There should be a file called ramspost.inp here. Exitting...'
  exit
fi

##########################################################################################################
# Determining the analysis prefix from the list                                                          #
##########################################################################################################
fprefix=`grep -i FPREFIX ramspost.inp`; fprefix=`echo ${fprefix} | sed s/" "/""/g |sed s/"'"/""/g`
ext=`echo ${fprefix} |wc -c`
p=0
###### Finding the comma position ########################################################################
while [ ${p} -lt ${ext} ]  # Finding the comma position
do
  p=`expr ${p} + 1`
  virgula=`echo "${fprefix}" | awk '{print substr($1,'${p}',1)}'`
  if [ y${virgula} = 'y,' ]; then virgula=${p}; p=${ext}; fi
done
ext=`expr ${virgula} - 9` # 9 is the position after the comma
fprefix=`echo ${fprefix} | awk '{print substr($1,9,'${ext}')}'`
##########################################################################################################

##########################################################################################################
# Determining the output file prefix from the namelist                                                   #
##########################################################################################################
gprefix=`grep -i GPREFIX ramspost.inp`; gprefix=`echo ${gprefix} | sed s/" "/""/g |sed s/"'"/""/g`
ext=`echo ${gprefix} |wc -c`
p=0
###### Finding the comma position ########################################################################
while [ ${p} -lt ${ext} ]                           
do
  p=`expr ${p} + 1`
  virgula=`echo "${gprefix}" | awk '{print substr($1,'${p}',1)}'`
  if [ y${virgula} = 'y,' ]; then virgula=${p}; p=${ext}; fi
done
ext=`expr ${virgula} - 9` # 9 corresponds to the position after the comma
gprefix=`echo ${gprefix} | awk '{print substr($1,9,'${ext}')}'`
##########################################################################################################


##########################################################################################################
# Running for each time, provided that it hasn't been run yet.                                           #
##########################################################################################################
quantosfiz=0
###### The list list is built in a way that doesn't crash when the list is too long ######################
diretorio=`dirname ${fprefix}`
radical=`basename ${fprefix}`
lista=`ls -1 ${diretorio}/ |grep ${radical} | grep head.txt${ending}`
###### Loop for each file ################################################################################
for analise in ${lista}
do
###### Finding the corresponding vfm file ################################################################
  analise=`basename ${analise}`
  vfm=`basename ${analise} head.txt${ending}`'g*.vfm'${ending}
###### Rearranging the the prefix of both analysis and output for a given time... ########################
  anapref=${tmpfolder}'/'`basename ${analise} \-head.txt${ending}`
  ext=`echo ${anapref} | wc -c`;p0=`expr ${ext} - 17`
  tempo=`echo ${anapref} | awk '{print substr($1,'${p0}',15)}'`
  saipref=${gprefix}${tempo}
###### Here I check whether the files already exist ######################################################
  if [ ! -s ${saipref}'_g1.gra' ]
  then
###### If they don't, I copy the header and data files into a temporary folder and uncompress them... ####
    cp ${diretorio}'/'${analise} ${diretorio}'/'${vfm} ${tmpfolder}
    ${nice} ${unzip} ${tmpfolder}'/'`basename ${analise}` ${tmpfolder}'/'`basename ${vfm}`
    quantosfiz=`expr ${quantosfiz} + 1`
    echo 'Running ramspost for time '${tempo}'...'
###### Switch the fprefix by the temporary folder, with the added restriction to the time ################
    cat ramspost.inp-salvaguarda | sed s@${fprefix}@${anapref}@g | sed s@${gprefix}@${saipref}@g > ramspost.inp
    ${nice} ${ramspost} > ${runoutput}
###### Cleaning the temporary folder, so I save disk space ###############################################
    rm -f ${tmpfolder}/*
  else
###### If the files exist, I skip this time ##############################################################
    echo 'I will not run for time '${tempo}' because those files already exist...'
  fi
done


##########################################################################################################
# Determining the number of grids                                                                        #
##########################################################################################################
gmax=`ls -1 ${fprefix}${tempo}* |wc -l`
gmax=`expr ${gmax} - 1`

##########################################################################################################
# Changing the output folder to avoid possible "./"                                                      #
##########################################################################################################
ext=`echo ${gprefix} | wc -c`
inicio=`echo ${gprefix} | awk '{print substr($1,1,2)}'`

rad=${gprefix}
if [ 'z'${inicio} = 'z./' ] 
then
  ext=`expr ${ext} - 2`
  gprefix=`echo ${gprefix} |awk '{print substr($1,3,'${ext}')}'`
fi
##########################################################################################################
# Here I find the time interval in sec.                                                                  #
##########################################################################################################
###### Time of the first file ############################################################################
primeiro=`ls -1 ${fprefix}*txt${ending} | head -1`
cp ${primeiro} './deleteme.txt'${ending}
###### Here I uncompress into a scratch file. Though it might sound inefficient, it's just a txt file... #
${unzip} './deleteme.txt'${ending}
primeiro='./deleteme.txt'
lmax=`cat ${primeiro} | wc -l`
l=0
###### Getting the time of the first analysis ############################################################
while [ ${l} -lt ${lmax} ]
do 
  l=`expr ${l} + 1`
  linha=`head -${l} ${primeiro} | tail -1`
  quale=`echo ${linha} | awk '{print $1}'`
  if [ z$quale = 'z__time' ]
  then 
    l=`expr ${l} + 2`
    linha=`head -${l} ${primeiro} | tail -1`
    primeiro=`echo ${linha} | awk '{print $1}'`
    l=${lmax}
  fi
done
rm ./deleteme.txt
###### Here I do the same thing for the second file... ###################################################
segundo=`ls -1 ${fprefix}*txt${ending} | head -2 | tail -1`
cp ${segundo} './deleteme.txt'${ending}
###### Here I uncompress into a scratch file. Though it might sound inefficient, it's just a txt file... #
${unzip} './deleteme.txt'${ending}
segundo='./deleteme.txt'
###### Finding out how many lines there are in the header file ###########################################
lmax=`cat ${segundo} | wc -l`
l=0
###### Getting the time of the second analysis ###########################################################
while [ ${l} -lt ${lmax} ];
do 
  l=`expr ${l} + 1`
  linha=`head -${l} ${segundo} | tail -1`
  quale=`echo ${linha} | awk '{print $1}'`
  if [ z$quale = 'z__time' ]
  then
    l=`expr ${l} + 2`
    linha=`head -${l} ${segundo} | tail -1`
    segundo=`echo ${linha} | awk '{print $1}'`
    l=${lmax}
  fi
done
rm -f ./deleteme.txt
##########################################################################################################
# Here I find the interval between the first and the second time...                                      #
###### Getting the time of the second analysis ###########################################################
###### Getting the time of the second analysis ###########################################################
primeiro=`echo ${primeiro} | sed s@'E+'@'*(10^'@g`')'          # Neither sh nor bc understand E+
segundo=`echo ${segundo} | sed s@'E+'@'*(10^'@g`')'
deltat=`echo ${segundo}' - '${primeiro}| bc`
deltat=`echo ${deltat} \/ 60 |bc`                              # Here is the interval in seconds.
###### Here I decide if I show the interval in minutes or seconds (GrADS doesn't deal with seconds...)         
deltathor=`echo ${deltat} \/ 60 |bc`
seriainteiro=`echo ${deltathor} \* 60 |bc`
if [ ${deltat} -eq ${seriainteiro} ]
  then deltat=${deltathor}'hr'
  else deltat=${deltat}'mn'
fi

##########################################################################################################
# Here I will generate a template for each grid, provided that some binary was created...                #
##########################################################################################################
if [ ${quantosfiz} -gt 0 ]
then
  g=0
  while [ ${g} -lt ${gmax} ]
  do
    g=`expr ${g} + 1`
    ctl=${rad}'_g'${g}'.ctl'
    ctl=`basename ${ctl}`
###### Taking the first file as a "canvas"... ############################################################
    ref=`ls -1 ${gprefix}????-??-??-????_?${g}'.ctl' |head -1`
    lmax=`cat ${ref} | wc -l`
###### Checking line by line the need to change. If it is unecessary, I just paste into the template #####
    l=0
    while [ ${l} -lt ${lmax} ]
    do
      l=`expr ${l} + 1`
      linha=`head -${l} ${ref} | tail -1`
      quale=`echo ${linha} |awk '{print $1}'`
###### Filename, I also use the opportunity to add the "template" line ################################### 
      if [ z${quale} = 'zdset' ]                                                             
      then         
        if [ 'z'`echo ${gprefix} | awk '{print substr($1,1,1)}'` = 'z/' ]
          then echo 'dset '${gprefix}'%y4-%m2-%d2-%h2%n2_g'${g}'.gra' > ${ctl}
          else echo 'dset ^'${gprefix}'%y4-%m2-%d2-%h2%n2_g'${g}'.gra' > ${ctl}
        fi
        echo 'options template' >> ${ctl}
###### Useless change, but you customize the title in the ctl... #########################################
      elif [ z${quale} = 'ztitle' ] 
      then
        echo 'title '${title} >> ${ctl}
      elif [ z${quale} = 'ztdef' ] #Altero o delta-T para o template
      then
        amax=`ls -1 ${gprefix}*${g}.gra |wc -l`
        tempoa=`echo ${linha} | awk '{print $4}'`
        echo 'tdef '${amax}' linear '${tempoa}' '${deltat} >> ${ctl}
      else
        echo ${linha} >> ${ctl}
      fi
    done  
###### Displaying the command line in the screen so lazy people can copy with the mouse ;-) ##############
    echo "gradsc -cl 'open ${ctl}'"
  done 
fi

##########################################################################################################
# Final touches...                                                                                       #
##########################################################################################################
###### Removing the intermediate ctl if you asked me for. I will always leave the first for template #####
if [ ${deleteintclt}='y' -o ${deleteintctl}='Y' ] 
then 
  howmany=`ls -1 ${gprefix}*????-??-??-????_g?.ctl 2> /dev/null | wc -l`
  if [ ${howmany} -ge 1 ]
  then
     howmany=`expr ${howmany} - 1`
     files=`ls -1 ${gprefix}*????-??-??-????_g?.ctl | tail -${howmany}`
     for file in ${files}
     do
       rm -f ${file}
     done
  fi
fi
###### Returning ramspost.inp to the original file #######################################################
mv -f ramspost.inp-salvaguarda ramspost.inp 
##########################################################################################################
