#!/bin/sh

here=`pwd`
there=`echo ${here} | sed s@/n/Moorcroft_Lab/Users@/n/moorcroft_scratch@g`
desc=`basename ${here}`
sitemet='/n/moorcroft_data/mlongo/data/ed2_data/site_met_driver'
hvdmet='/n/home11/aantonarakis/EDrelease65/run/'
bioinit='/n/moorcroft_data/mlongo/data/ed2_data/site_bio_data'
sheffield='SHEF_NCEP_DRIVER_DS314'
lonlat=${here}'/joborder.txt'

#----- History run variables. -------------------------------------------------------------#
forcehisto=0     # Impose history start (0 = no, 1 = yes).  The following variables will 
                 #    be used only whent forcehisto is 1.
fullygrown='/n/moorcroftfs1/mlongo/EDBRAMS/debug/dbg_033/pdg_crash/histo/pedegigante'
yearh='1510'
monthh='07'
dateh='01'
timeh='0000'

toldef='0.01'
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
   polyname=`echo ${oi}    | awk '{print $1 }'`
   polyiata=`echo ${oi}    | awk '{print $2 }'`
   polylon=`echo ${oi}     | awk '{print $3 }'`
   polylat=`echo ${oi}     | awk '{print $4 }'`
   yeara=`echo ${oi}       | awk '{print $5 }'`
   montha=`echo ${oi}      | awk '{print $6 }'`
   datea=`echo ${oi}       | awk '{print $7 }'`
   timea=`echo ${oi}       | awk '{print $8 }'`
   yearz=`echo ${oi}       | awk '{print $9 }'`
   monthz=`echo ${oi}      | awk '{print $10}'`
   datez=`echo ${oi}       | awk '{print $11}'`
   timez=`echo ${oi}       | awk '{print $12}'`
   polyisoil=`echo ${oi}   | awk '{print $13}'`
   polyntext=`echo ${oi}   | awk '{print $14}'`
   polysand=`echo ${oi}    | awk '{print $15}'`
   polyclay=`echo ${oi}    | awk '{print $16}'`
   polydepth=`echo ${oi}   | awk '{print $17}'`
   polycol=`echo ${oi}     | awk '{print $18}'`
   slzres=`echo ${oi}      | awk '{print $19}'`
   queue=`echo ${oi}       | awk '{print $20}'`
   metdriver=`echo ${oi}   | awk '{print $21}'`
   dtlsm=`echo ${oi}       | awk '{print $22}'`
   vmfactc3=`echo ${oi}    | awk '{print $23}'`
   vmfactc4=`echo ${oi}    | awk '{print $24}'`
   mphotoc3=`echo ${oi}    | awk '{print $25}'`
   mphotoc4=`echo ${oi}    | awk '{print $26}'`
   kwgrass=`echo ${oi}     | awk '{print $27}'`
   kwtree=`echo ${oi}      | awk '{print $28}'`
   gammac3=`echo ${oi}     | awk '{print $29}'`
   gammac4=`echo ${oi}     | awk '{print $30}'`
   d0grass=`echo ${oi}     | awk '{print $31}'`
   d0tree=`echo ${oi}      | awk '{print $32}'`
   d0decay=`echo ${oi}     | awk '{print $33}'`
   alphac3=`echo ${oi}     | awk '{print $34}'`
   alphac4=`echo ${oi}     | awk '{print $35}'`
   klowco2=`echo ${oi}     | awk '{print $36}'`
   rrffact=`echo ${oi}     | awk '{print $37}'`
   growthresp=`echo ${oi}  | awk '{print $38}'`
   h2olimit=`echo ${oi}    | awk '{print $39}'`
   isfclyrm=`echo ${oi}    | awk '{print $40}'`
   icanturb=`echo ${oi}    | awk '{print $41}'`
   ubmin=`echo ${oi}       | awk '{print $42}'`
   ugbmin=`echo ${oi}      | awk '{print $43}'`
   ustmin=`echo ${oi}      | awk '{print $44}'`
   gamm=`echo ${oi}        | awk '{print $45}'`
   gamh=`echo ${oi}        | awk '{print $46}'`
   tprandtl=`echo ${oi}    | awk '{print $47}'`
   ribmax=`echo ${oi}      | awk '{print $48}'`
   atmco2=`echo ${oi}      | awk '{print $49}'`
   thcrit=`echo ${oi}      | awk '{print $50}'`
   smfire=`echo ${oi}      | awk '{print $51}'`
   isoilbc=`echo ${oi}     | awk '{print $52}'`
   imetrad=`echo ${oi}     | awk '{print $53}'`
   ibranch=`echo ${oi}     | awk '{print $54}'`
   icanrad=`echo ${oi}     | awk '{print $55}'`
   crown=`echo   ${oi}     | awk '{print $56}'`
   ltransvis=`echo ${oi}   | awk '{print $57}'`
   lreflectvis=`echo ${oi} | awk '{print $58}'`
   ltransnir=`echo ${oi}   | awk '{print $59}'`
   lreflectnir=`echo ${oi} | awk '{print $60}'`
   orienttree=`echo ${oi}  | awk '{print $61}'`
   orientgrass=`echo ${oi} | awk '{print $62}'`
   clumptree=`echo ${oi}   | awk '{print $63}'`
   clumpgrass=`echo ${oi}  | awk '{print $64}'`
   ivegtdyn=`echo ${oi}    | awk '{print $65}'`
   igndvap=`echo ${oi}     | awk '{print $66}'`
   iphen=`echo ${oi}       | awk '{print $67}'`
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
      ;;
   Caxiuana)
      metdriverdb=${sitemet}'/Caxiuana/Caxiuana06_HEADER'
      metcyc1=1999
      metcycf=2003
      imetavg=1
      ;;
   Harvard)
      metdriverdb=${hvdmet}'/HARVARD_MET_93_09'
      metcyc1=1993
      metcycf=2008
      imetavg=1
      ;;
   Manaus_KM34)
      metdriverdb=${sitemet}'/Manaus_KM34/Manaus_KM34_HEADER'
      metcyc1=2002
      metcycf=2005
      imetavg=1
      ;;
   Reserva_Jaru)
      metdriverdb=${sitemet}'/Reserva_Jaru/Reserva_Jaru_HEADER'
      metcyc1=2000
      metcycf=2002
      imetavg=1
      ;;
   Reserva_Pe-de-Gigante)
      metdriverdb=${sitemet}'/Reserva_Pe_de_Gigante/Reserva_Pe-de-Gigante_HEADER'
      metcyc1=2001
      metcycf=2003
      imetavg=1
      ;;
   Santarem_KM67)
      metdriverdb=${sitemet}'/Santarem_KM67/Santarem_KM67_HEADER'
      metcyc1=2002
      metcycf=2004
      imetavg=1
      ;;
   Santarem_KM77)
      metdriverdb=${sitemet}'/Santarem_KM77/Santarem_KM77_HEADER'
      metcyc1=2001
      metcycf=2005
      imetavg=1
      ;;
   Santarem_KM83)
      metdriverdb=${sitemet}'/Santarem_KM83/Santarem_KM83_HEADER'
      metcyc1=2001
      metcycf=2003
      imetavg=1
      ;;
   Fazenda_NS)
      metdriverdb=${sitemet}'/Fazenda_Nossa_Senhora/Fazenda_Nossa_Senhora_HEADER'
      metcyc1=1999
      metcycf=2001
      imetavg=1
      ;;
   Guyaflux)
      metdriverdb=${sitemet}'/Guyaflux/Guyaflux_HEADER'
      metcyc1=2007
      metcycf=2009
      imetavg=1
      ;;
   *)
      metdriverdb=${here}/${polyname}/${sheffield}
      metcyc1=1969
      metcycf=2008
      imetavg=2
      ;;
   esac
   #---------------------------------------------------------------------------------------#






   #---------------------------------------------------------------------------------------#
   #     Determine which soil profile to use.  We have eight categories (A-H), and the     #
   # soil resolution.                                                                      #
   #---------------------------------------------------------------------------------------#
   polyslz1=''
   polyslz2=''
   polyslz3=''
   polyslz4=''
   polyslz5=''
   polyslz6=''
   polyslz7=''
   polyslm1=''
   polyslm2=''
   polyslm3=''
   polyslm4=''
   polyslm5=''
   polyslm6=''
   polyslm7=''
   polyslt1=''
   polyslt2=''
   polyslt3=''
   polyslt4=''
   polyslt5=''
   polyslt6=''
   polyslt7=''
   case ${slzres} in
   0)
      case ${polydepth} in
      A)
         polynzg=16
         polyslz1='-1.250,-1.135,-1.024,-0.917,-0.814,-0.715,-0.620,-0.530,-0.445,-0.364,'
         polyslz2='-0.289,-0.221,-0.158,-0.103,-0.056,-0.020'
         polyslm1=' 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,'
         polyslm2=' 1.000, 1.000, 1.000, 1.000, 1.000, 1.000'
         polyslt1=' 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,'
         polyslt2=' 0.000, 0.000, 0.000, 0.000, 0.000, 0.000'
         ;;
      B)
         polynzg=16
         polyslz1='-2.000,-1.797,-1.602,-1.417,-1.240,-1.073,-0.916,-0.769,-0.632,-0.507,'
         polyslz2='-0.392,-0.290,-0.200,-0.124,-0.063,-0.020'
         polyslm1=' 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,'
         polyslm2=' 1.000, 1.000, 1.000, 1.000, 1.000, 1.000'
         polyslt1=' 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,'
         polyslt2=' 0.000, 0.000, 0.000, 0.000, 0.000, 0.000'
         ;;
      C)
         polynzg=16
         polyslz1='-3.000,-2.670,-2.357,-2.061,-1.784,-1.524,-1.283,-1.061,-0.857,-0.673,'
         polyslz2='-0.510,-0.367,-0.245,-0.146,-0.070,-0.020'
         polyslm1=' 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,'
         polyslm2=' 1.000, 1.000, 1.000, 1.000, 1.000, 1.000'
         polyslt1=' 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,'
         polyslt2=' 0.000, 0.000, 0.000, 0.000, 0.000, 0.000'
         ;;
      D)
         polynzg=16
         polyslz1='-4.000,-3.536,-3.099,-2.690,-2.308,-1.955,-1.629,-1.332,-1.064,-0.824,'
         polyslz2='-0.614,-0.433,-0.283,-0.163,-0.075,-0.020'
         polyslm1=' 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,'
         polyslm2=' 1.000, 1.000, 1.000, 1.000, 1.000, 1.000'
         polyslt1=' 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,'
         polyslt2=' 0.000, 0.000, 0.000, 0.000, 0.000, 0.000'
         ;;
      E)
         polynzg=16
         polyslz1='-4.500,-3.967,-3.467,-3.000,-2.565,-2.164,-1.797,-1.462,-1.162,-0.895,'
         polyslz2='-0.662,-0.464,-0.300,-0.171,-0.077,-0.020'
         polyslm1=' 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,'
         polyslm2=' 1.000, 1.000, 1.000, 1.000, 1.000, 1.000'
         polyslt1=' 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,'
         polyslt2=' 0.000, 0.000, 0.000, 0.000, 0.000, 0.000'
         ;;
      F)
         polynzg=16
         polyslz1='-6.000,-5.254,-4.559,-3.914,-3.320,-2.776,-2.282,-1.837,-1.442,-1.095,'
         polyslz2='-0.798,-0.548,-0.346,-0.192,-0.083,-0.020'
         polyslm1=' 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,'
         polyslm2=' 1.000, 1.000, 1.000, 1.000, 1.000, 1.000'
         polyslt1=' 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,'
         polyslt2=' 0.000, 0.000, 0.000, 0.000, 0.000, 0.000'
         ;;
      G)
         polynzg=16
         polyslz1='-7.000,-6.108,-5.279,-4.514,-3.812,-3.172,-2.593,-2.076,-1.618,-1.221,'
         polyslz2='-0.881,-0.600,-0.374,-0.204,-0.087,-0.020'
         polyslm1=' 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,'
         polyslm2=' 1.000, 1.000, 1.000, 1.000, 1.000, 1.000'
         polyslt1=' 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,'
         polyslt2=' 0.000, 0.000, 0.000, 0.000, 0.000, 0.000'
         ;;
      H)
         polynzg=16
         polyslz1='-8.000,-6.959,-5.995,-5.108,-4.296,-3.560,-2.897,-2.307,-1.789,-1.340,'
         polyslz2='-0.961,-0.648,-0.400,-0.215,-0.089,-0.020'
         polyslm1=' 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,'
         polyslm2=' 1.000, 1.000, 1.000, 1.000, 1.000, 1.000'
         polyslt1=' 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,'
         polyslt2=' 0.000, 0.000, 0.000, 0.000, 0.000, 0.000'
         ;;
      *)
         polynzg=16
         polyslz1='-8.000,-6.959,-5.995,-5.108,-4.296,-3.560,-2.897,-2.307,-1.789,-1.340,'
         polyslz2='-0.961,-0.648,-0.400,-0.215,-0.089,-0.020'
         polyslm1=' 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,'
         polyslm2=' 1.000, 1.000, 1.000, 1.000, 1.000, 1.000'
         polyslt1=' 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,'
         polyslt2=' 0.000, 0.000, 0.000, 0.000, 0.000, 0.000'
         ;;
      esac
      ;;
   1)
      case ${polydepth} in
      A)
         polynzg='30'
         polyslz1='-1.329,-1.246,-1.168,-1.093,-1.022,-0.955,-0.890,-0.829,-0.770,-0.714,'
         polyslm1=' 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,'
         polyslt1=' 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,'
         polyslz2='-0.661,-0.611,-0.563,-0.517,-0.473,-0.432,-0.392,-0.354,-0.318,-0.284,'
         polyslm2=' 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,'
         polyslt2=' 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,'
         polyslz3='-0.252,-0.221,-0.191,-0.163,-0.136,-0.111,-0.086,-0.063,-0.041,-0.020'
         polyslm3=' 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000'
         polyslt3=' 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000'
         ;;
      B)
         polynzg='37'
         polyslz1='-2.033,-1.917,-1.806,-1.701,-1.601,-1.506,-1.415,-1.329,-1.246,-1.168,'
         polyslm1=' 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,'
         polyslt1=' 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,'
         polyslz2='-1.093,-1.022,-0.955,-0.890,-0.829,-0.770,-0.714,-0.661,-0.611,-0.563,'
         polyslm2=' 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,'
         polyslt2=' 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,'
         polyslz3='-0.517,-0.473,-0.432,-0.392,-0.354,-0.318,-0.284,-0.252,-0.221,-0.191,'
         polyslm3=' 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,'
         polyslt3=' 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,'
         polyslz4='-0.163,-0.136,-0.111,-0.086,-0.063,-0.041,-0.020'
         polyslm4=' 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000'
         polyslt4=' 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000'
         ;;
      C)
         polynzg='44'
         polyslz1='-3.023,-2.860,-2.705,-2.557,-2.416,-2.282,-2.154,-2.033,-1.917,-1.806,'
         polyslm1=' 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,'
         polyslt1=' 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,'
         polyslz2='-1.701,-1.601,-1.506,-1.415,-1.329,-1.246,-1.168,-1.093,-1.022,-0.955,'
         polyslm2=' 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,'
         polyslt2=' 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,'
         polyslz3='-0.890,-0.829,-0.770,-0.714,-0.661,-0.611,-0.563,-0.517,-0.473,-0.432,'
         polyslm3=' 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,'
         polyslt3=' 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,'
         polyslz4='-0.392,-0.354,-0.318,-0.284,-0.252,-0.221,-0.191,-0.163,-0.136,-0.111,'
         polyslm4=' 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,'
         polyslt4=' 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,'
         polyslz5='-0.086,-0.063,-0.041,-0.020'
         polyslm5=' 1.000, 1.000, 1.000, 1.000'
         polyslt5=' 0.000, 0.000, 0.000, 0.000'
         ;;
      D)
         polynzg='50'
         polyslz1='-4.187,-3.969,-3.761,-3.562,-3.374,-3.194,-3.023,-2.860,-2.705,-2.557,'
         polyslm1=' 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,'
         polyslt1=' 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,'
         polyslz2='-2.416,-2.282,-2.154,-2.033,-1.917,-1.806,-1.701,-1.601,-1.506,-1.415,'
         polyslm2=' 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,'
         polyslt2=' 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,'
         polyslz3='-1.329,-1.246,-1.168,-1.093,-1.022,-0.955,-0.890,-0.829,-0.770,-0.714,'
         polyslm3=' 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,'
         polyslt3=' 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,'
         polyslz4='-0.661,-0.611,-0.563,-0.517,-0.473,-0.432,-0.392,-0.354,-0.318,-0.284,'
         polyslm4=' 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,'
         polyslt4=' 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,'
         polyslz5='-0.252,-0.221,-0.191,-0.163,-0.136,-0.111,-0.086,-0.063,-0.041,-0.020'
         polyslm5=' 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000'
         polyslt5=' 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000'
         ;;
      E)
         polynzg='52'
         polyslz1='-4.657,-4.416,-4.187,-3.969,-3.761,-3.562,-3.374,-3.194,-3.023,-2.860,'
         polyslm1=' 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,'
         polyslt1=' 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,'
         polyslz2='-2.705,-2.557,-2.416,-2.282,-2.154,-2.033,-1.917,-1.806,-1.701,-1.601,'
         polyslm2=' 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,'
         polyslt2=' 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,'
         polyslz3='-1.506,-1.415,-1.329,-1.246,-1.168,-1.093,-1.022,-0.955,-0.890,-0.829,'
         polyslm3=' 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,'
         polyslt3=' 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,'
         polyslz4='-0.770,-0.714,-0.661,-0.611,-0.563,-0.517,-0.473,-0.432,-0.392,-0.354,'
         polyslm4=' 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,'
         polyslt4=' 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,'
         polyslz5='-0.318,-0.284,-0.252,-0.221,-0.191,-0.163,-0.136,-0.111,-0.086,-0.063,'
         polyslm5=' 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,'
         polyslt5=' 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,'
         polyslz6='-0.041,-0.020'
         polyslm6=' 1.000, 1.000'
         polyslt6=' 0.000, 0.000'
         ;;
      F)
         polynzg='58'
         polyslz1='-6.157,-5.907,-5.657,-5.407,-5.157,-4.907,-4.657,-4.416,-4.187,-3.969,'
         polyslm1=' 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,'
         polyslt1=' 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,'
         polyslz2='-3.761,-3.562,-3.374,-3.194,-3.023,-2.860,-2.705,-2.557,-2.416,-2.282,'
         polyslm2=' 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,'
         polyslt2=' 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,'
         polyslz3='-2.154,-2.033,-1.917,-1.806,-1.701,-1.601,-1.506,-1.415,-1.329,-1.246,'
         polyslm3=' 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,'
         polyslt3=' 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,'
         polyslz4='-1.168,-1.093,-1.022,-0.955,-0.890,-0.829,-0.770,-0.714,-0.661,-0.611,'
         polyslm4=' 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,'
         polyslt4=' 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,'
         polyslz5='-0.563,-0.517,-0.473,-0.432,-0.392,-0.354,-0.318,-0.284,-0.252,-0.221,'
         polyslm5=' 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,'
         polyslt5=' 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,'
         polyslz6='-0.191,-0.163,-0.136,-0.111,-0.086,-0.063,-0.041,-0.020'
         polyslm6=' 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000'
         polyslt6=' 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000'
         ;;
      G)
         polynzg='62'
         polyslz1='-7.157,-6.907,-6.657,-6.407,-6.157,-5.907,-5.657,-5.407,-5.157,-4.907,'
         polyslm1=' 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,'
         polyslt1=' 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,'
         polyslz2='-4.657,-4.416,-4.187,-3.969,-3.761,-3.562,-3.374,-3.194,-3.023,-2.860,'
         polyslm2=' 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,'
         polyslt2=' 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,'
         polyslz3='-2.705,-2.557,-2.416,-2.282,-2.154,-2.033,-1.917,-1.806,-1.701,-1.601,'
         polyslm3=' 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,'
         polyslt3=' 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,'
         polyslz4='-1.506,-1.415,-1.329,-1.246,-1.168,-1.093,-1.022,-0.955,-0.890,-0.829,'
         polyslm4=' 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,'
         polyslt4=' 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,'
         polyslz5='-0.770,-0.714,-0.661,-0.611,-0.563,-0.517,-0.473,-0.432,-0.392,-0.354,'
         polyslm5=' 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,'
         polyslt5=' 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,'
         polyslz6='-0.318,-0.284,-0.252,-0.221,-0.191,-0.163,-0.136,-0.111,-0.086,-0.063,'
         polyslm6=' 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,'
         polyslt6=' 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,'
         polyslz7='-0.041,-0.020'
         polyslm7=' 1.000, 1.000'
         polyslt7=' 0.000, 0.000'
         ;;
      H)
         polynzg='66'
         polyslz1='-8.157,-7.907,-7.657,-7.407,-7.157,-6.907,-6.657,-6.407,-6.157,-5.907,'
         polyslm1=' 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,'
         polyslt1=' 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,'
         polyslz2='-5.657,-5.407,-5.157,-4.907,-4.657,-4.416,-4.187,-3.969,-3.761,-3.562,'
         polyslm2=' 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,'
         polyslt2=' 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,'
         polyslz3='-3.374,-3.194,-3.023,-2.860,-2.705,-2.557,-2.416,-2.282,-2.154,-2.033,'
         polyslm3=' 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,'
         polyslt3=' 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,'
         polyslz4='-1.917,-1.806,-1.701,-1.601,-1.506,-1.415,-1.329,-1.246,-1.168,-1.093,'
         polyslm4=' 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,'
         polyslt4=' 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,'
         polyslz5='-1.022,-0.955,-0.890,-0.829,-0.770,-0.714,-0.661,-0.611,-0.563,-0.517,'
         polyslm5=' 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,'
         polyslt5=' 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,'
         polyslz6='-0.473,-0.432,-0.392,-0.354,-0.318,-0.284,-0.252,-0.221,-0.191,-0.163,'
         polyslm6=' 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,'
         polyslt6=' 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,'
         polyslz7='-0.136,-0.111,-0.086,-0.063,-0.041,-0.020'
         polyslm7=' 1.000, 1.000, 1.000, 1.000, 1.000, 1.000'
         polyslt7=' 0.000, 0.000, 0.000, 0.000, 0.000, 0.000'
         ;;
      *)
         ;;
      esac
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
         thissfilin=${bioinit}'/km67_ustein_newallom.'
         ;;
      m34)
         thissfilin=${bioinit}'/k34_ustein_newallom.'
         ;;
      bdf)
         thissfilin=${bioinit}'/pBDFFP1_1983_ustein.'
         ;;
      pdg)
         thissfilin=${bioinit}'/pdg_ustein_newallom.'
         ;;
      fns)
         thissfilin=${bioinit}'/fns.'
         ;;
      s77)
         thissfilin=${bioinit}'/km77.'
         ;;
      cax)
         thissfilin=${bioinit}'/cax06_ustein_newallom.'
         ;;
      rja)
         thissfilin=${bioinit}'/rja_ustein_newallom.'
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
   sed -i s@paththere@${there}@g             ${ED2IN}
   sed -i s@myyeara@${thisyeara}@g           ${ED2IN}
   sed -i s@mymontha@${montha}@g             ${ED2IN}
   sed -i s@mydatea@${datea}@g               ${ED2IN}
   sed -i s@mytimea@${timea}@g               ${ED2IN}
   sed -i s@myyearz@${thisyearz}@g           ${ED2IN}
   sed -i s@mymonthz@${monthz}@g             ${ED2IN}
   sed -i s@mydatez@${datez}@g               ${ED2IN}
   sed -i s@mytimez@${timez}@g               ${ED2IN}
   sed -i s@mydtlsm@${dtlsm}@g               ${ED2IN}
   sed -i s@thispoly@${polyname}@g           ${ED2IN}
   sed -i s@plonflag@${polylon}@g            ${ED2IN}
   sed -i s@platflag@${polylat}@g            ${ED2IN}
   sed -i s@timehhhh@${time}@g               ${ED2IN}
   sed -i s@datehhhh@${date}@g               ${ED2IN}
   sed -i s@monthhhh@${month}@g              ${ED2IN}
   sed -i s@yearhhhh@${year}@g               ${ED2IN}
   sed -i s@myinitmode@${initmode}@g         ${ED2IN}
   sed -i s@mysfilin@${thissfilin}@g         ${ED2IN}
   sed -i s@mytrees@${pfts}@g                ${ED2IN}
   sed -i s@mycrop@${crop}@g                 ${ED2IN}
   sed -i s@myplantation@${plantation}@g     ${ED2IN}
   sed -i s@myiphen@${iphen}@g               ${ED2IN}
   sed -i s@myisoilflg@${polyisoil}@g        ${ED2IN}
   sed -i s@mynslcon@${polyntext}@g          ${ED2IN}
   sed -i s@myslxsand@${polysand}@g          ${ED2IN}
   sed -i s@myslxclay@${polyclay}@g          ${ED2IN}
   sed -i s@mysoilcol@${polycol}@g           ${ED2IN}
   sed -i s@mynzg@${polynzg}@g               ${ED2IN}
   sed -i s@mymetdriverdb@${metdriverdb}@g   ${ED2IN}
   sed -i s@mymetcyc1@${metcyc1}@g           ${ED2IN}
   sed -i s@mymetcycf@${metcycf}@g           ${ED2IN}
   sed -i s@mytoler@${toler}@g               ${ED2IN}
   sed -i s@RUNFLAG@${runt}@g                ${ED2IN}
   sed -i s@myvmfactc3@${vmfactc3}@g         ${ED2IN}
   sed -i s@myvmfactc4@${vmfactc4}@g         ${ED2IN}
   sed -i s@mymphotoc3@${mphotoc3}@g         ${ED2IN}
   sed -i s@mymphotoc4@${mphotoc4}@g         ${ED2IN}
   sed -i s@mykwgrass@${kwgrass}@g           ${ED2IN}
   sed -i s@mykwtree@${kwtree}@g             ${ED2IN}
   sed -i s@mygammac3@${gammac3}@g           ${ED2IN}
   sed -i s@mygammac4@${gammac4}@g           ${ED2IN}
   sed -i s@myd0grass@${d0grass}@g           ${ED2IN}
   sed -i s@myd0tree@${d0tree}@g             ${ED2IN}
   sed -i s@myd0decay@${d0decay}@g           ${ED2IN}
   sed -i s@myalphac3@${alphac3}@g           ${ED2IN}
   sed -i s@myalphac4@${alphac4}@g           ${ED2IN}
   sed -i s@myklowco2@${klowco2}@g           ${ED2IN}
   sed -i s@myrrffact@${rrffact}@g           ${ED2IN}
   sed -i s@mygrowthresp@${growthresp}@g     ${ED2IN}
   sed -i s@myh2olimit@${h2olimit}@g         ${ED2IN}
   sed -i s@mysfclyrm@${isfclyrm}@g          ${ED2IN}
   sed -i s@myicanturb@${icanturb}@g         ${ED2IN}
   sed -i s@myatmco2@${atmco2}@g             ${ED2IN}
   sed -i s@mythcrit@${thcrit}@g             ${ED2IN}
   sed -i s@mysmfire@${smfire}@g             ${ED2IN}
   sed -i s@mymetavg@${imetavg}@g            ${ED2IN}
   sed -i s@mysoilbc@${isoilbc}@g            ${ED2IN}
   sed -i s@mymetrad@${imetrad}@g            ${ED2IN}
   sed -i s@mybranch@${ibranch}@g            ${ED2IN}
   sed -i s@mycanrad@${icanrad}@g            ${ED2IN}
   sed -i s@mycrown@${crown}@g               ${ED2IN}
   sed -i s@myltransvis@${ltransvis}@g       ${ED2IN}
   sed -i s@myltransnir@${ltransnir}@g       ${ED2IN}
   sed -i s@mylreflectvis@${lreflectvis}@g   ${ED2IN}
   sed -i s@mylreflectnir@${lreflectnir}@g   ${ED2IN}
   sed -i s@myorienttree@${orienttree}@g     ${ED2IN}
   sed -i s@myorientgrass@${orientgrass}@g   ${ED2IN}
   sed -i s@myclumptree@${clumptree}@g       ${ED2IN}
   sed -i s@myclumpgrass@${clumpgrass}@g     ${ED2IN}
   sed -i s@myvegtdyn@${ivegtdyn}@g          ${ED2IN}
   sed -i s@myubmin@${ubmin}@g               ${ED2IN}
   sed -i s@myugbmin@${ugbmin}@g             ${ED2IN}
   sed -i s@myustmin@${ustmin}@g             ${ED2IN}
   sed -i s@mygamm@${gamm}@g                 ${ED2IN}
   sed -i s@mygamh@${gamh}@g                 ${ED2IN}
   sed -i s@mytprandtl@${tprandtl}@g         ${ED2IN}
   sed -i s@myribmax@${ribmax}@g             ${ED2IN}
   sed -i s@mygndvap@${igndvap}@g            ${ED2IN}

   #------ Soil variables. ----------------------------------------------------------------#
   sed -i s@myslz1@"${polyslz1}"@g           ${ED2IN}
   sed -i s@myslz2@"${polyslz2}"@g           ${ED2IN}
   sed -i s@myslz3@"${polyslz3}"@g           ${ED2IN}
   sed -i s@myslz4@"${polyslz4}"@g           ${ED2IN}
   sed -i s@myslz5@"${polyslz5}"@g           ${ED2IN}
   sed -i s@myslz6@"${polyslz6}"@g           ${ED2IN}
   sed -i s@myslz7@"${polyslz7}"@g           ${ED2IN}
   sed -i s@myslmstr1@"${polyslm1}"@g        ${ED2IN}
   sed -i s@myslmstr2@"${polyslm2}"@g        ${ED2IN}
   sed -i s@myslmstr3@"${polyslm3}"@g        ${ED2IN}
   sed -i s@myslmstr4@"${polyslm4}"@g        ${ED2IN}
   sed -i s@myslmstr5@"${polyslm5}"@g        ${ED2IN}
   sed -i s@myslmstr6@"${polyslm6}"@g        ${ED2IN}
   sed -i s@myslmstr7@"${polyslm7}"@g        ${ED2IN}
   sed -i s@mystgoff1@"${polyslt1}"@g        ${ED2IN}
   sed -i s@mystgoff2@"${polyslt2}"@g        ${ED2IN}
   sed -i s@mystgoff3@"${polyslt3}"@g        ${ED2IN}
   sed -i s@mystgoff4@"${polyslt4}"@g        ${ED2IN}
   sed -i s@mystgoff5@"${polyslt5}"@g        ${ED2IN}
   sed -i s@mystgoff6@"${polyslt6}"@g        ${ED2IN}
   sed -i s@mystgoff7@"${polyslt7}"@g        ${ED2IN}


   #----- Change the srun.sh file. --------------------------------------------------------#
   srun=${here}'/'${polyname}'/srun.sh'
   sed -i s@pathhere@${here}@g     ${srun}
   sed -i s@thispoly@${polyname}@g ${srun}
   sed -i s@thisdesc@${desc}@g     ${srun}

   #----- Change the callserial.sh file. --------------------------------------------------#
   callserial=${here}'/'${polyname}'/callserial.sh'
   sed -i s@thisroot@${here}@g ${callserial}
   sed -i s@thispoly@${polyname}@g ${callserial}
   sed -i s@myexec@${execname}@g ${callserial}

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
