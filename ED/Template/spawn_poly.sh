#!/bin/sh

here=`pwd`
there=`echo ${here} | sed s@/n/Moorcroft_Lab/Users@/n/moorcroft_scratch@g`
desc=`basename ${here}`
sitemet='/n/moorcroft_scratch/nlevine/data/ed2_data/site_met_driver'
hvdmet='/n/home11/aantonarakis/EDrelease65/run/'
bioinit='/n/moorcroft_scratch/nlevine/data/ed2_data/site_bio_data'
sheffield='/n/moorcroft_scratch/kzhang/Sheffield_Met_Data/southam_shef_driver/SHEF_NCEP_DRIVER_DS314'
lonlat=${here}'/joborder.txt'

#----- History run variables. -------------------------------------------------------------#
forcehisto=0     # Impose history start (0 = no, 1 = yes).  The following variables will 
                 #    be used only whent forcehisto is 1.
fullygrown='/n/moorcroftfs1/mlongo/EDBRAMS/debug/dbg_033/pdg_crash/histo/pedegigante'
yearh='1510'
monthh='07'
dateh='01'
timeh='0000'

toldef='0.001'
initmode=6

callunpa=${here}'/callunpa.sh'
unparun=${here}'/unparun.sh'
hostlist=${here}'/hostlist.txt'

#----- Executable name. -------------------------------------------------------------------#
execname='ed_2.1-opt'


#----- Determine the number of polygons to run. -------------------------------------------#
let npolys=`wc -l ${lonlat} | awk '{print $1 }'`-3
echo 'Number of polygons: '${npolys}'...'

#----- Start a new callunpa.sh. -----------------------------------------------------------#
rm -f ${callunpa}
echo '#!/bin/sh' > ${callunpa}

#----- Start a new host list. -------------------------------------------------------------#
rm -f ${hostlist}
touch ${hostlist}
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#   Check whether the executable is copied.  If not, let the user know and stop the        #
# script.                                                                                  #
#------------------------------------------------------------------------------------------#
if [ ! -s ${here}/executable/${execname} ]
then
   echo 'Executable file : '${execname}' is not in the executable directory'
   echo 'Copy the executable to the file before running this script!'
   exit 99
fi





#------------------------------------------------------------------------------------------#
#     Loop over all polygons.                                                              #
#------------------------------------------------------------------------------------------#
ff=0
unpa=0
mc2=0
lastunpa='none'
while [ ${ff} -lt ${npolys} ]
do
   let ff=${ff}+1
   let line=${ff}+3

   #---------------------------------------------------------------------------------------#
   #      Read the ffth line of the polygon list.  There must be smarter ways of doing     #
   # this, but this works.  Here we obtain the polygon name, and its longitude and         #
   # latitude.                                                                             #
   #---------------------------------------------------------------------------------------#
   oi=`head -${line} ${lonlat} | tail -1`
   polyname=`echo ${oi}  | awk '{print $1 }'`
   polyiata=`echo ${oi}  | awk '{print $2 }'`
   polylon=`echo ${oi}   | awk '{print $3 }'`
   polylat=`echo ${oi}   | awk '{print $4 }'`
   yeara=`echo ${oi}     | awk '{print $5 }'`
   montha=`echo ${oi}    | awk '{print $6 }'`
   datea=`echo ${oi}     | awk '{print $7 }'`
   timea=`echo ${oi}     | awk '{print $8 }'`
   yearz=`echo ${oi}     | awk '{print $9 }'`
   monthz=`echo ${oi}    | awk '{print $10}'`
   datez=`echo ${oi}     | awk '{print $11}'`
   timez=`echo ${oi}     | awk '{print $12}'`
   polyisoil=`echo ${oi} | awk '{print $13}'`
   polyntext=`echo ${oi} | awk '{print $14}'`
   polysand=`echo ${oi}  | awk '{print $15}'`
   polyclay=`echo ${oi}  | awk '{print $16}'`
   polydepth=`echo ${oi} | awk '{print $17}'`
   queue=`echo ${oi}     | awk '{print $18}'`
   metdriver=`echo ${oi} | awk '{print $19}'`
   dtlsm=`echo ${oi}     | awk '{print $20}'`
   vmfact=`echo ${oi}    | awk '{print $21}'`
   mfact=`echo ${oi}     | awk '{print $22}'`
   kfact=`echo ${oi}     | awk '{print $23}'`
   gamfact=`echo ${oi}   | awk '{print $24}'`
   lwfact=`echo ${oi}    | awk '{print $25}'`
   betaflag=`echo ${oi}  | awk '{print $26}'`
   thioff=`echo ${oi}    | awk '{print $27}'`
   ustmin=`echo ${oi}    | awk '{print $28}'`
   ggfact=`echo ${oi}    | awk '{print $29}'`
   wlimit=`echo ${oi}    | awk '{print $30}'`
   blyrcnd=`echo ${oi}   | awk '{print $31}'`
   iallom=`echo ${oi}    | awk '{print $32}'`
   icanturb=`echo ${oi}  | awk '{print $33}'`
   isfclyrm=`echo ${oi}  | awk '{print $34}'`
   gamm=`echo ${oi}      | awk '{print $35}'`
   gamh=`echo ${oi}      | awk '{print $36}'`
   tprandtl=`echo ${oi}  | awk '{print $37}'`
   vh2vr=`echo ${oi}     | awk '{print $38}'`
   vh2dh=`echo ${oi}     | awk '{print $39}'`
   ribmax=`echo ${oi}    | awk '{print $40}'`
   maxwhc=`echo ${oi}    | awk '{print $41}'`
   runoff=`echo ${oi}    | awk '{print $42}'`
   atmco2=`echo ${oi}    | awk '{print $43}'`
   thcrit=`echo ${oi}    | awk '{print $44}'`
   smfire=`echo ${oi}    | awk '{print $45}'`
   agefall=`echo ${oi}   | awk '{print $46}'`
   grndvap=`echo ${oi}   | awk '{print $47}'`
   crownmod=`echo ${oi}  | awk '{print $48}'`
   quantum=`echo ${oi}   | awk '{print $49}'`
   #---------------------------------------------------------------------------------------#



   #----- Find the actual beta power that goes to the namelist. ---------------------------#
   bpower=`calc.sh ${betaflag}/2`
   #---------------------------------------------------------------------------------------#



   #----- Check whether the directories exist or not, and stop the script if they do. -----#
   if [ -s ${here}/${polyname} ]
   then
      echo 'Order: '${ff}', updating a couple of files on '${here}/${polyname}'...'
      
      #----- Save the last tolerance in case we are going to make it more strict. ---------#
      oldtol=`grep RK4_TOLERANCE ${here}/${polyname}/ED2IN | awk '{print $3}'`
      rm -f ${here}/${polyname}/ED2IN 
      rm -f ${here}/${polyname}/callserial.sh
      rm -f ${here}/${polyname}/callunpa.sh 
      rm -f ${here}/${polyname}/skipper.txt
      cp ${here}/Template/ED2IN         ${here}/${polyname}/ED2IN
      cp ${here}/Template/callserial.sh ${here}/${polyname}/callserial.sh
   else
      echo 'Order: '${ff}', creating directory for polygon '${polyname}'...'
      #----- Copy the Template directory to a unique polygon directory. -------------------#
      oldtol=''
      cp -r ${here}/Template ${here}/${polyname}
   fi
   #---------------------------------------------------------------------------------------#






   #---------------------------------------------------------------------------------------#
   #   Make sure that we have a reasonable tolerance.                                      #
   #---------------------------------------------------------------------------------------#
   if [ 'x'${oldtol} == 'xmytoler' -o 'x'${oldtol} == 'x' ]
   then
      toler=${toldef}
   else
      toler=${oldtol}
   fi
   #---------------------------------------------------------------------------------------#






   #---------------------------------------------------------------------------------------#
   #      Here we check whether the an old directory exists or not.  If it doesn't, then   #
   # we simply copy the template, and assume initial run.  Otherwise, we must find out     #
   # where the simulation was when it stopped.                                             #
   #---------------------------------------------------------------------------------------#
   if [ -s  ${here}/${polyname} ]
   then

      #------------------------------------------------------------------------------------#
      #      This step is necessary because we may have killed the run while it was        #
      # writing, and as a result, the file may be corrupt.                                 #
      #------------------------------------------------------------------------------------#
      nhdf5=`ls -1 ${here}/${polyname}/histo/* 2> /dev/null | wc -l`
      if [ ${nhdf5} -gt 0 ]
      then

         h5fine=0

         while [ $h5fine -eq 0 ]
         do
            lasthdf5=`ls -1 ${here}/${polyname}/histo/* | tail -1`
            h5dump -H ${lasthdf5} 1> /dev/null 2> ${here}/badfile.txt

            if [ -s ${here}/badfile.txt ]
            then
               /bin/rm -fv ${lasthdf5}
               nhdf5=`ls -1 ${here}/${polyname}/histo/* 2> /dev/null | wc -l`
               if [ ${nhdf5} -eq 0 ]
               then
                  hdf5fine=1
               fi
            else
               h5fine=1
            fi

            /bin/rm -f ${here}/badfile.txt
         done
      fi
      #------------------------------------------------------------------------------------#






      #------------------------------------------------------------------------------------#
      #      Run the small R script to check whether the simulation was running or not,    #
      # and whether there was any cohort left by the time the runs were stopped.           #
      #------------------------------------------------------------------------------------#
      sed -i s@thispoly@${polyname}@g ${here}/${polyname}/whichrun.r
      sed -i s@thisqueue@${queue}@g   ${here}/${polyname}/whichrun.r
      sed -i s@pathhere@${here}@g     ${here}/${polyname}/whichrun.r
      sed -i s@paththere@${there}@g   ${here}/${polyname}/whichrun.r
      R CMD BATCH ${here}/${polyname}/whichrun.r ${here}/${polyname}/outwhichrun.txt
      while [ ! -s ${here}/${polyname}/statusrun.txt ]
      do
         sleep 0.5
      done
      year=`cat ${here}/${polyname}/statusrun.txt  | awk '{print $2}'`
      month=`cat ${here}/${polyname}/statusrun.txt | awk '{print $3}'`
      date=`cat ${here}/${polyname}/statusrun.txt  | awk '{print $4}'`
      time=`cat ${here}/${polyname}/statusrun.txt  | awk '{print $5}'`
      runt=`cat ${here}/${polyname}/statusrun.txt  | awk '{print $6}'`
      #------------------------------------------------------------------------------------#

   else
      sed -i s@thispoly@${polyname}@g ${here}/${polyname}/whichrun.r
      sed -i s@thisqueue@${queue}@g   ${here}/${polyname}/whichrun.r
      sed -i s@pathhere@${here}@g     ${here}/${polyname}/whichrun.r
      sed -i s@paththere@${there}@g   ${here}/${polyname}/whichrun.r
      year=${yeara}
      month=${montha}
      date=${datea}
      time=${timea}
      runt='INITIAL'
   fi
   #---------------------------------------------------------------------------------------#












   #---------------------------------------------------------------------------------------#
   #     Determine which PFTs to use based on the "iata" code.                             #
   #---------------------------------------------------------------------------------------#
   case ${polyiata} in
   wch|zmh|nqn)
      pfts='5,6,7,8,9,10,11,17'
      crop=5
      plantation=17
      ;;
   hvd)
      pfts='5,6,8,9,10,11'
      crop=5
      plantation=8
      ;;
   aei|asu|cnf|bnu|cwb|erm|igp|iqq|ipv|mgf|rao|sla|zpe|kna)
      pfts='1,2,3,4,16,17'
      crop=16
      plantation=17
      ;;
   fns)
      pfts='1,16'
      crop=1
      plantation=17
      ;;
   s77)
      pfts='1,16'
      crop=16
      plantation=17
      ;;
   *)
      pfts='1,2,3,4,16'
      crop=1
      plantation=3
      ;;
   esac
   #---------------------------------------------------------------------------------------#






   #---------------------------------------------------------------------------------------#
   #     Determine which meteorological data set to use.  Default is the Sheffield/NCEP    #
   # dataset, otherwise the site-level tower data is used.                                 #
   #---------------------------------------------------------------------------------------#
   case ${metdriver} in
   Bananal_Island)
      metdriverdb=${sitemet}'/Bananal_Island/Bananal_HEADER'
      metcyc1=2004
      metcycf=2006
      imetavg=1
      iphen=2
      ;;
   Caxiuana)
      metdriverdb=${sitemet}'/Caxiuana/Caxiuana06_HEADER'
      metcyc1=1999
      metcycf=2003
      imetavg=1
      iphen=2
      ;;
   Harvard)
      metdriverdb=${hvdmet}'/HARVARD_MET_93_09'
      metcyc1=1993
      metcycf=2008
      imetavg=1
      iphen=1
      ;;
   Manaus_KM34)
      metdriverdb=${sitemet}'/Manaus_KM34/Manaus_KM34_HEADER'
      metcyc1=2002
      metcycf=2005
      imetavg=1
      iphen=2
      ;;
   Reserva_Jaru)
      metdriverdb=${sitemet}'/Reserva_Jaru/Reserva_Jaru_HEADER'
      metcyc1=2000
      metcycf=2002
      imetavg=1
      iphen=2
      ;;
   Reserva_Pe-de-Gigante)
      metdriverdb=${sitemet}'/Reserva_Pe_de_Gigante/Reserva_Pe-de-Gigante_HEADER'
      metcyc1=2001
      metcycf=2003
      imetavg=1
      iphen=2
      ;;
   Santarem_KM67)
      metdriverdb=${sitemet}'/Santarem_KM67/Santarem_KM67_HEADER'
      metcyc1=2002
      metcycf=2004
      imetavg=1
      iphen=2
      ;;
   Santarem_KM77)
      metdriverdb=${sitemet}'/Fazenda_Nossa_Senhora/Santarem_KM77_HEADER'
      metcyc1=2001
      metcycf=2005
      imetavg=1
      iphen=2
      ;;
   Fazenda_NS)
      metdriverdb=${sitemet}'/Fazenda_Nossa_Senhora/Fazenda_Nossa_Senhora_HEADER'
      metcyc1=1999
      metcycf=2001
      imetavg=1
      iphen=2
      ;;
   *)
      metdriverdb=${sheffield}
      metcyc1=1969
      metcycf=2008
      imetavg=0
      iphen=2
      ;;
   esac
   #---------------------------------------------------------------------------------------#






   #---------------------------------------------------------------------------------------#
   #     Determine which soil profile to use.  We have eight categories (A-H)              #
   #---------------------------------------------------------------------------------------#
   case ${polydepth} in
   A)
      polynzg=10
      polyslz='-1.25,-1.00,-0.80,-0.60,-0.45,-0.30,-0.20,-0.12,-0.06,-0.02'
      polyslm=' 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60'
      polyslt=' 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00'
      ;;
   B)
      polynzg=11
      polyslz='-2.00,-1.50,-1.00,-0.80,-0.60,-0.45,-0.30,-0.20,-0.12,-0.06,-0.02'
      polyslm=' 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60'
      polyslt=' 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00'
      ;;
   C)
      polynzg=13
      polyslz='-3.00,-2.25,-1.75,-1.35,-1.00,-0.80,-0.60,-0.45,-0.30,-0.20,-0.12,-0.06,-0.02'
      polyslm=' 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60'
      polyslt=' 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00'
      ;;
   D)
      polynzg=14
      polyslz='-4.00,-3.25,-2.50,-2.00,-1.35,-1.00,-0.80,-0.60,-0.45,-0.30,-0.20,-0.12,-0.06,-0.02'
      polyslm=' 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60'
      polyslt=' 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00'
      ;;
   E)
      polynzg=14
      polyslz='-4.50,-3.50,-2.75,-2.00,-1.35,-1.00,-0.80,-0.60,-0.45,-0.30,-0.20,-0.12,-0.06,-0.02'
      polyslm=' 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60'
      polyslt=' 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00'
      ;;
   F)
      polynzg=15
      polyslz='-6.00,-4.75,-3.75,-2.75,-2.00,-1.35,-1.00,-0.80,-0.60,-0.45,-0.30,-0.20,-0.12,-0.06,-0.02'
      polyslm=' 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60'
      polyslt=' 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00'
      ;;
   G)
      polynzg=15
      polyslz='-7.00,-5.50,-4.50,-3.50,-2.80,-2.20,-1.70,-1.20,-0.80,-0.50,-0.30,-0.20,-0.12,-0.06,-0.02'
      polyslm=' 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60'
      polyslt=' 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00'
      ;;
   H)
      polynzg=16
      polyslz='-8.00,-6.50,-5.50,-4.50,-3.50,-2.80,-2.20,-1.70,-1.20,-0.80,-0.50,-0.30,-0.20,-0.12,-0.06,-0.02'
      polyslm=' 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60'
      polyslt=' 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00'
      ;;
   AM1)
      polynzg=9
      polyslz='-1.25,-1.00,-0.80,-0.60,-0.45,-0.30,-0.20,-0.12,-0.06'
      polyslm=' 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60'
      polyslt=' 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00'
      ;;
   BM1)
      polynzg=10
      polyslz='-2.00,-1.50,-1.00,-0.80,-0.60,-0.45,-0.30,-0.20,-0.12,-0.06'
      polyslm=' 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60'
      polyslt=' 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00'
      ;;
   CM1)
      polynzg=12
      polyslz='-3.00,-2.25,-1.75,-1.35,-1.00,-0.80,-0.60,-0.45,-0.30,-0.20,-0.12,-0.06'
      polyslm=' 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60'
      polyslt=' 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00'
      ;;
   DM1)
      polynzg=13
      polyslz='-4.00,-3.25,-2.50,-2.00,-1.35,-1.00,-0.80,-0.60,-0.45,-0.30,-0.20,-0.12,-0.06'
      polyslm=' 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60'
      polyslt=' 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00'
      ;;
   EM1)
      polynzg=13
      polyslz='-4.50,-3.50,-2.75,-2.00,-1.35,-1.00,-0.80,-0.60,-0.45,-0.30,-0.20,-0.12,-0.06'
      polyslm=' 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60'
      polyslt=' 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00'
      ;;
   FM1)
      polynzg=14
      polyslz='-6.00,-4.75,-3.75,-2.75,-2.00,-1.35,-1.00,-0.80,-0.60,-0.45,-0.30,-0.20,-0.12,-0.06'
      polyslm=' 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60'
      polyslt=' 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00'
      ;;
   GM1)
      polynzg=14
      polyslz='-7.00,-5.50,-4.50,-3.50,-2.80,-2.20,-1.70,-1.20,-0.80,-0.50,-0.30,-0.20,-0.12,-0.06'
      polyslm=' 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60'
      polyslt=' 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00'
      ;;
   HM1)
      polynzg=15
      polyslz='-8.00,-6.50,-5.50,-4.50,-3.50,-2.80,-2.20,-1.70,-1.20,-0.80,-0.50,-0.30,-0.20,-0.12,-0.06'
      polyslm=' 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60'
      polyslt=' 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00'
      ;;
   *)
      polynzg=16
      polyslz='-8.00,-6.50,-5.50,-4.50,-3.50,-2.80,-2.20,-1.70,-1.20,-0.80,-0.50,-0.30,-0.20,-0.12,-0.06,-0.02'
      polyslm=' 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60'
      polyslt=' 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00'
      ;;
   esac
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Define whether we use the met cycle to define the first and last year, or the     #
   # default year.                                                                         #
   #---------------------------------------------------------------------------------------#
   if [ ${yeara} -eq 0 ]
   then
      thisyeara=${metcyc1}
   else
      thisyeara=${yeara}
   fi
   if [ ${yearz} -eq 0 ]
   then
      thisyearz=${metcycf}
   else
      thisyearz=${yearz}
   fi



   #---------------------------------------------------------------------------------------#
   #     Change the ED2IN file.                                                            #
   #---------------------------------------------------------------------------------------#
   if [ ${runt}  == 'CRASHED' ]
   then
      sed -i s@CRASHED@HISTORY@g ${here}/${polyname}/statusrun.txt
      runt='HISTORY'
      toler=`calc.sh ${toler}/10`
   fi
   #---------------------------------------------------------------------------------------#

   #----- Check whether to use SFILIN as restart or history. ------------------------------#
   if [ ${runt} == 'INITIAL' ] && [ ${forcehisto} -eq 1 ]
   then
      runt='HISTORY'
      year=${yearh}
      month=${monthh}
      date=${dateh}
      time=${timeh}
      thissfilin=${fullygrown}
   elif [ ${initmode} -eq 6 ]
   then
      thissfilin=${fullygrown}
      case ${polyiata} in
      s67)
         thissfilin=${bioinit}'/km67_ustein.'
         ;;
      m34)
         thissfilin=${bioinit}'/pBDFFP1_1983_ustein.'
         ;;
      pdg)
         thissfilin=${bioinit}'/pdg_grass.'
         ;;
      fns)
         thissfilin=${bioinit}'/fns.'
         ;;
      *)
         echo ' Polygon: '${polyname}
         echo ' IATA: '${polyiata}
         echo 'This IATA cannot be used by biomass initialisation!'
         exit 59
         ;;
      esac
   else
      thissfilin=${here}/${polyname}/histo/${polyname}
   fi
   #---------------------------------------------------------------------------------------#

   ED2IN=${here}'/'${polyname}'/ED2IN'
   sed -i s@paththere@${there}@g           ${ED2IN}
   sed -i s@myyeara@${thisyeara}@g         ${ED2IN}
   sed -i s@mymontha@${montha}@g           ${ED2IN}
   sed -i s@mydatea@${datea}@g             ${ED2IN}
   sed -i s@mytimea@${timea}@g             ${ED2IN}
   sed -i s@myyearz@${thisyearz}@g         ${ED2IN}
   sed -i s@mymonthz@${monthz}@g           ${ED2IN}
   sed -i s@mydatez@${datez}@g             ${ED2IN}
   sed -i s@mytimez@${timez}@g             ${ED2IN}
   sed -i s@mydtlsm@${dtlsm}@g             ${ED2IN}
   sed -i s@thispoly@${polyname}@g         ${ED2IN}
   sed -i s@plonflag@${polylon}@g          ${ED2IN}
   sed -i s@platflag@${polylat}@g          ${ED2IN}
   sed -i s@timehhhh@${time}@g             ${ED2IN}
   sed -i s@datehhhh@${date}@g             ${ED2IN}
   sed -i s@monthhhh@${month}@g            ${ED2IN}
   sed -i s@yearhhhh@${year}@g             ${ED2IN}
   sed -i s@myinitmode@${initmode}@g       ${ED2IN}
   sed -i s@mysfilin@${thissfilin}@g       ${ED2IN}
   sed -i s@mytrees@${pfts}@g              ${ED2IN}
   sed -i s@mycrop@${crop}@g               ${ED2IN}
   sed -i s@myplantation@${plantation}@g   ${ED2IN}
   sed -i s@myiphen@${iphen}@g             ${ED2IN}
   sed -i s@myisoilflg@${polyisoil}@g      ${ED2IN}
   sed -i s@mynslcon@${polyntext}@g        ${ED2IN}
   sed -i s@myslxsand@${polysand}@g        ${ED2IN}
   sed -i s@myslxclay@${polyclay}@g        ${ED2IN}
   sed -i s@mynzg@${polynzg}@g             ${ED2IN}
   sed -i s@myslz@"${polyslz}"@g           ${ED2IN}
   sed -i s@myslmstr@"${polyslm}"@g        ${ED2IN}
   sed -i s@mystgoff@"${polyslt}"@g        ${ED2IN}
   sed -i s@mymetdriverdb@${metdriverdb}@g ${ED2IN}
   sed -i s@mymetcyc1@${metcyc1}@g         ${ED2IN}
   sed -i s@mymetcycf@${metcycf}@g         ${ED2IN}
   sed -i s@mytoler@${toler}@g             ${ED2IN}
   sed -i s@RUNFLAG@${runt}@g              ${ED2IN}
   sed -i s@myvmfact@${vmfact}@g           ${ED2IN}
   sed -i s@mymfact@${mfact}@g             ${ED2IN}
   sed -i s@mykfact@${kfact}@g             ${ED2IN}
   sed -i s@mywlimit@${wlimit}@g           ${ED2IN}
   sed -i s@mygamfact@${gamfact}@g         ${ED2IN}
   sed -i s@mylwfact@${lwfact}@g           ${ED2IN}
   sed -i s@mythioff@${thioff}@g           ${ED2IN}
   sed -i s@myallom@${iallom}@g            ${ED2IN}
   sed -i s@myblyrcnd@${blyrcnd}@g         ${ED2IN}
   sed -i s@mybpower@${bpower}@g           ${ED2IN}
   sed -i s@myicanturb@${icanturb}@g       ${ED2IN}
   sed -i s@myisfclyrm@${isfclyrm}@g       ${ED2IN}
   sed -i s@myustmin@${ustmin}@g           ${ED2IN}
   sed -i s@myggfact@${ggfact}@g           ${ED2IN}
   sed -i s@mygamm@${gamm}@g               ${ED2IN}
   sed -i s@mygamh@${gamh}@g               ${ED2IN}
   sed -i s@mytprandtl@${tprandtl}@g       ${ED2IN}
   sed -i s@myvh2vr@${vh2vr}@g             ${ED2IN}
   sed -i s@myvh2dh@${vh2dh}@g             ${ED2IN}
   sed -i s@myribmax@${ribmax}@g           ${ED2IN}
   sed -i s@mymaxwhc@${maxwhc}@g           ${ED2IN}
   sed -i s@myrunoff@${runoff}@g           ${ED2IN}
   sed -i s@myatmco2@${atmco2}@g           ${ED2IN}
   sed -i s@mythcrit@${thcrit}@g           ${ED2IN}
   sed -i s@mysmfire@${smfire}@g           ${ED2IN}
   sed -i s@myagefall@${agefall}@g         ${ED2IN}
   sed -i s@mygrndvap@${grndvap}@g         ${ED2IN}
   sed -i s@mycrownmod@${crownmod}@g       ${ED2IN}
   sed -i s@mymetavg@${imetavg}@g          ${ED2IN}
   sed -i s@myquantum@${quantum}@g         ${ED2IN}

   #----- Change the srun.sh file. --------------------------------------------------------#
   srun=${here}'/'${polyname}'/srun.sh'
   sed -i s@pathhere@${here}@g     ${srun}
   sed -i s@thispoly@${polyname}@g ${srun}
   sed -i s@thisdesc@${desc}@g     ${srun}

   #----- Change the callserial.sh file. --------------------------------------------------#
   callserial=${here}'/'${polyname}'/callserial.sh'
   sed -i s@thisroot@${here}@g ${callserial}
   sed -i s@thispoly@${polyname}@g ${callserial}

   if [ ${queue} == 'GC3' ]
   then
      queue='long_serial'
   fi

   if [ ${queue} == 'unrestricted_parallel' ]
   then
      sed -i s@thisqueue@GC3@g ${srun}
      echo 'Skip this node for now.' > ${here}/${polyname}/skipper.txt
   else
      sed -i s@thisqueue@${queue}@g ${srun}
   fi

   #----- Make the shell script for the unrestricted_parallel. ----------------------------#
   if [ ${queue} == 'unrestricted_parallel' ] && [ ${runt} == 'HISTORY' -o ${runt} == 'INITIAL' ]
   then
      touch ${callunpa}

      #----- Add command to the node. -----------------------------------------------------#
      let unpa=${unpa}+1
      let submit=${unpa}%8 # % is mod operator
      let wtime=${submit}*15
      let wtime=${wtime}+2
      mycomm="${here}/${polyname}/callserial.sh ${wtime} &"
      echo ${mycomm} >> ${callunpa}
      
      lastunpa=${polyname}

      #----- Finish the shell script and put it inside the node. --------------------------#
      if [ ${submit} -eq 0 ]
      then
         #----- Move callunpa to the polygon. ---------------------------------------------#
         echo 'wait' >> ${callunpa}
         chmod u+x ${callunpa}
         mv ${callunpa} ${here}/${polyname} 

         #----- Start a new callunpa. -----------------------------------------------------#
         echo '#!/bin/sh' > ${callunpa}

         #----- Create the shell script that will call bsub. ------------------------------#
         echo '#!/bin/sh' > ${unparun}
         bsub='bsub -q unrestricted_parallel'
         bsub=${bsub}' -J '${polyname}
         bsub=${bsub}' -o '${here}'/'${polyname}'/serial_lsf.out -n '${unpa}
         bsub=${bsub}' < '${here}'/'${polyname}'/'`basename ${callunpa}`
         echo ${bsub} >> ${unparun}
         chmod u+x ${unparun}
         mv ${unparun} ${here}/${polyname}

         #----- Reset unpa for next job. --------------------------------------------------#
         lastunpa='none'
         unpa=0
      fi
   fi
done

#----- Make sure that all jobs were submitted to unrestricted serial. ---------------------#
if [ ${lastunpa} == 'none' ]
then
   /bin/rm -f ${unparun}
else
   #----- Move callunpa to the polygon. ---------------------------------------------------#
   echo 'wait' >> ${callunpa}
   chmod u+x ${callunpa}
   mv ${callunpa} ${here}/${lastunpa} 

   #----- Create the shell script that will call bsub. ------------------------------------#
   echo '#!/bin/sh' > ${unparun}
   bsub='bsub -q unrestricted_parallel'
   bsub=${bsub}' -J '${polyname}
   bsub=${bsub}' -o '${here}'/'${lastunpa}'/serial_lsf.out -n '${unpa}
   bsub=${bsub}' < '${here}'/'${lastunpa}'/'`basename ${callunpa}`
   echo ${bsub} >> ${unparun}
   chmod u+x ${unparun}
   mv ${unparun} ${here}/${lastunpa}

   #----- Reset unpa for next job. --------------------------------------------------------#
   lastunpa='none'
   unpa=0
fi
