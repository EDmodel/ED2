#!/bin/bash
source $HOME/.bashrc

#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#    MAIN SETTINGS:                                                                        #
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#


#----- Main path, you must actually type the path in case you add to chron. ---------------#
here='xxxxxxxxxxxxxxxxxxxxx'
if [ ${here} == 'xxxxxxxxxxxxxxxxxxxxx' ]
then
   echo ' You must set up variable here before using run_sitter.sh!!!'
   exit 99
fi
#------------------------------------------------------------------------------------------#


#----- Path with some utilities for run_sitter.sh (this script). --------------------------#
situtils="${here}/sit_utils"
#------------------------------------------------------------------------------------------#


#----- Where the output is. ---------------------------------------------------------------#
there=${here}
#------------------------------------------------------------------------------------------#

#----- User name, usually set by `whoami` so you don't need to change it. -----------------#
moi=`whoami`
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#    Main actions.                                                                         #
#    printevery -- when looping through jobs, print something on screen every printevery   #
#                  polygons                                                                #
#    bigloop    -- Go through the big loop (0 - no, 1 - yes)                               #
#    fixqueue   -- 0: perform the full queue management)                                   #
#                  1: simply fix the queues on joborder.                                   #
#------------------------------------------------------------------------------------------#
printevery=1
bigloop=1
fixqueue=0
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#  waitpot -- Should polygons wait until some previous run is done? ('n' will start them   #
#             right away).                                                                 #
#  potveg  -- Path where to look for initialisation files in case waitpot = 'y'.           #
#------------------------------------------------------------------------------------------#
waitpot='n'
potveg='/n/moorcroft_data/data/ed2_data/restarts_sci_006/potveg'
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
# copyrestart -- Should I copy the latest history file somewhere else ('y'/'N').           #
# restart     -- Path to where to copy restart files in case copyrestart = 'y'.            #
#------------------------------------------------------------------------------------------#
copyrestart='n'
restart='/n/moorcroft_data/mlongo/ed2_data/restarts_XXX'
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#     Useful files to check the simulation progress:                                       #
# joborder  -- File with the job instructions and model settings                           #
# lastcheck -- File with the report from last time                                         #
# outcheck  -- File with the current report                                                #
# situation -- Situation of the run                                                        #
#------------------------------------------------------------------------------------------#
joborder="${here}/joborder.txt"
lastcheck="${situtils}/lastcheck.txt"
outcheck="${situtils}/mycheck.txt"
situation="${situtils}/situation.txt"
#------------------------------------------------------------------------------------------#


#------ Calculator. -----------------------------------------------------------------------#
ccc='/n/Moorcroft_Lab/Users/mlongo/util/calc.sh'  # Calculator
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#    Maximum load in each queue.                                                           #
#    - full: maximum number of jobs in the queue                                           #
#    - umax: maximum number of jobs in the queue that I am allowed to run                  #
#    - pmax: maximum number of pending jobs that we may leave                              #
#------------------------------------------------------------------------------------------#
#------ Queue: moorcroft2a. ---------------------------------------------------------------#
m2afull=96
m2aumax=80
#------ Queue: moorcroft2b. ---------------------------------------------------------------#
m2bfull=88
m2bumax=88
#------ Queue: moorcroft2c. ---------------------------------------------------------------#
m2cfull=64
m2cumax=64
#------ Queue: moorcroft_6100a. -----------------------------------------------------------#
a61full=48
a61umax=36
#------ Queue: moorcroft_6100b. -----------------------------------------------------------#
b61full=480
b61umax=380
#------ Queue: wofsy. ---------------------------------------------------------------------#
wsyfull=96
wsyumax=24
#------ Queue: unrestricted_serial. -------------------------------------------------------#
usefull=80
useumax=80
#------ Queue: unrestricted_parallel. -----------------------------------------------------#
upapmax=0
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#    E-mail options (normally only the first three options may require changes.            #
#    email1day  -- Should I e-mail once a day (1) or every time I run (0)?                 #
#    recipient  -- To which e-mail should I send the update?                               #
#    plotstatus -- Plot the current status of the run (you will need to edit the R script  #
#                  for each simulation). 0: no; 1: yes                                     #
#    emailbody  -- File that will contain the e-mail.                                      #
#    headfile   -- File with header                                                        #
#    tailfile   -- File with "bye"                                                         #
#    recefile   -- File with the recent activity                                           #
#    statfile   -- File with the current status                                            #
#    queuefile  -- File with queue status                                                  #
#    pendfile   -- File with pending status                                                #
#    email1day  -- Reminder so the script knows whether an e-mail has been sent or not.    #
#------------------------------------------------------------------------------------------#
email1day=1
recipient='xxxxxx@xxxx.com'
plotstatus='n'
emailbody="${situtils}/email.txt"
headfile="${situtils}/head.txt"
tailfile="${situtils}/tail.txt"
recefile="${situtils}/rece.txt"
statfile="${situtils}/stat.txt"
queuefile="${situtils}/queue.txt"
pendfile=${here}'/pending.txt'
emailday=${here}'/emailday.txt'
if [ ${recipient} == "xxxxxx@xxxx.com" ]
then
   echo 'You must specify the e-mail (variable recipient)'
   stop 92
fi
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#      Names for the special unrestricted_parallel scripts.                                #
#------------------------------------------------------------------------------------------#
callunpa=${here}'/callunpa.sh'
unparun=${here}'/unparun.sh'
#------------------------------------------------------------------------------------------#



#----- Unique job description for these simulations (it really must be unique). -----------#
desc=`basename ${here}`
#------------------------------------------------------------------------------------------#





#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#                  THINK TWICE BEFORE CHANGING ANYTHING BEYOND THIS POINT.                 #
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#





#----- Check whether run_sitter.sh is still running or not.  If it is, exit. --------------#
if [ -s ${here}/run_sitter.lock ]
then
   exit
else
   echo 'I am going to check your run.' > ${here}/run_sitter.lock
fi
#------------------------------------------------------------------------------------------#



#----- Determine the number of polygons to check. -----------------------------------------#
let npolys=`wc -l ${joborder} | awk '{print $1 }'`-3
echo 'Number of polygons: '${npolys}'...'
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#     Reset most counters.                                                                 #
#------------------------------------------------------------------------------------------#
ff=0          # Polygon counter.
nstart=0      # Number of new polygons that initialised
nmetmiss=0    # Number of new polygons that crashed due to missing met driver
nstopped=0    # Number of new polygons that crashed due to unknown reason
ncrashed=0    # Number of new polygons that crashed due to numeric instability
newststate=0  # Number of new polygons at steady state
newextinct=0  # Number of new extinct polygons
nresubmit=0   # Number of long_serial polygons that were submitted again.
deathrow=""   # List of polygons to be killed (they went extinct or reached steady state)
#------------------------------------------------------------------------------------------#



if [ ${bigloop} -ne 0 ]
then
   #---------------------------------------------------------------------------------------#
   #     Eliminate the latest check, and move the most recent check to the latest check,   #
   # and start the new check.                                                              #
   #---------------------------------------------------------------------------------------#
   rm -f ${lastcheck}
   mv ${outcheck} ${lastcheck}
   touch ${outcheck}
   rm -f ${situation}
   touch ${situation}
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Loop over all polygons.                                                           #
   #---------------------------------------------------------------------------------------#
   while [ ${ff} -lt ${npolys} ]
   do
      let ff=${ff}+1
      let line=${ff}+3
      let rem=${ff}%${printevery}

      #------------------------------------------------------------------------------------#
      #     Plot a banner to entertain the user.                                           #
      #------------------------------------------------------------------------------------#
      if [ ${rem} -eq 0 -o ${ff} -eq ${npolys} ]
      then
         echo ' - Checked '${ff}' polygons so far...'
      fi
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #      Read the polyth line of the polygon list.  There must be smarter ways of do-  #
      # ing this, but this works.  Here we obtain the polygon name, and its longitude and  #
      # latitude.                                                                          #
      #------------------------------------------------------------------------------------#
      oi=`head -${line} ${lonlat} | tail -1`
      polyname=`echo ${oi}     | awk '{print $1 }'`
      polyiata=`echo ${oi}     | awk '{print $2 }'`
      polylon=`echo ${oi}      | awk '{print $3 }'`
      polylat=`echo ${oi}      | awk '{print $4 }'`
      yeara=`echo ${oi}        | awk '{print $5 }'`
      montha=`echo ${oi}       | awk '{print $6 }'`
      datea=`echo ${oi}        | awk '{print $7 }'`
      timea=`echo ${oi}        | awk '{print $8 }'`
      yearz=`echo ${oi}        | awk '{print $9 }'`
      monthz=`echo ${oi}       | awk '{print $10}'`
      datez=`echo ${oi}        | awk '{print $11}'`
      timez=`echo ${oi}        | awk '{print $12}'`
      initmode=`echo ${oi}     | awk '{print $13}'`
      iscenario=`echo ${oi}    | awk '{print $14}'`
      isizepft=`echo ${oi}     | awk '{print $15}'`
      polyisoil=`echo ${oi}    | awk '{print $16}'`
      polyntext=`echo ${oi}    | awk '{print $17}'`
      polysand=`echo ${oi}     | awk '{print $18}'`
      polyclay=`echo ${oi}     | awk '{print $19}'`
      polydepth=`echo ${oi}    | awk '{print $20}'`
      polysoilbc=`echo ${oi}   | awk '{print $21}'`
      polysldrain=`echo ${oi}  | awk '{print $22}'`
      polycol=`echo ${oi}      | awk '{print $23}'`
      slzres=`echo ${oi}       | awk '{print $24}'`
      queue=`echo ${oi}        | awk '{print $25}'`
      metdriver=`echo ${oi}    | awk '{print $26}'`
      dtlsm=`echo ${oi}        | awk '{print $27}'`
      vmfactc3=`echo ${oi}     | awk '{print $28}'`
      vmfactc4=`echo ${oi}     | awk '{print $29}'`
      mphototrc3=`echo ${oi}   | awk '{print $30}'`
      mphototec3=`echo ${oi}   | awk '{print $31}'`
      mphotoc4=`echo ${oi}     | awk '{print $32}'`
      bphotoblc3=`echo ${oi}   | awk '{print $33}'`
      bphotonlc3=`echo ${oi}   | awk '{print $34}'`
      bphotoc4=`echo ${oi}     | awk '{print $35}'`
      kwgrass=`echo ${oi}      | awk '{print $36}'`
      kwtree=`echo ${oi}       | awk '{print $37}'`
      gammac3=`echo ${oi}      | awk '{print $38}'`
      gammac4=`echo ${oi}      | awk '{print $39}'`
      d0grass=`echo ${oi}      | awk '{print $40}'`
      d0tree=`echo ${oi}       | awk '{print $41}'`
      alphac3=`echo ${oi}      | awk '{print $42}'`
      alphac4=`echo ${oi}      | awk '{print $43}'`
      klowco2=`echo ${oi}      | awk '{print $44}'`
      decomp=`echo ${oi}       | awk '{print $45}'`
      rrffact=`echo ${oi}      | awk '{print $46}'`
      growthresp=`echo ${oi}   | awk '{print $47}'`
      lwidthgrass=`echo ${oi}  | awk '{print $48}'`
      lwidthbltree=`echo ${oi} | awk '{print $49}'`
      lwidthnltree=`echo ${oi} | awk '{print $50}'`
      q10c3=`echo ${oi}        | awk '{print $51}'`
      q10c4=`echo ${oi}        | awk '{print $52}'`
      h2olimit=`echo ${oi}     | awk '{print $53}'`
      imortscheme=`echo ${oi}  | awk '{print $54}'`
      ddmortconst=`echo ${oi}  | awk '{print $55}'`
      isfclyrm=`echo ${oi}     | awk '{print $56}'`
      icanturb=`echo ${oi}     | awk '{print $57}'`
      ubmin=`echo ${oi}        | awk '{print $58}'`
      ugbmin=`echo ${oi}       | awk '{print $59}'`
      ustmin=`echo ${oi}       | awk '{print $60}'`
      gamm=`echo ${oi}         | awk '{print $61}'`
      gamh=`echo ${oi}         | awk '{print $62}'`
      tprandtl=`echo ${oi}     | awk '{print $63}'`
      ribmax=`echo ${oi}       | awk '{print $64}'`
      atmco2=`echo ${oi}       | awk '{print $65}'`
      thcrit=`echo ${oi}       | awk '{print $66}'`
      smfire=`echo ${oi}       | awk '{print $67}'`
      ifire=`echo ${oi}        | awk '{print $68}'`
      fireparm=`echo ${oi}     | awk '{print $69}'`
      ipercol=`echo ${oi}      | awk '{print $70}'`
      runoff=`echo ${oi}       | awk '{print $71}'`
      imetrad=`echo ${oi}      | awk '{print $72}'`
      ibranch=`echo ${oi}      | awk '{print $73}'`
      icanrad=`echo ${oi}      | awk '{print $74}'`
      crown=`echo   ${oi}      | awk '{print $75}'`
      ltransvis=`echo ${oi}    | awk '{print $76}'`
      lreflectvis=`echo ${oi}  | awk '{print $77}'`
      ltransnir=`echo ${oi}    | awk '{print $78}'`
      lreflectnir=`echo ${oi}  | awk '{print $79}'`
      orienttree=`echo ${oi}   | awk '{print $80}'`
      orientgrass=`echo ${oi}  | awk '{print $81}'`
      clumptree=`echo ${oi}    | awk '{print $82}'`
      clumpgrass=`echo ${oi}   | awk '{print $83}'`
      ivegtdyn=`echo ${oi}     | awk '{print $84}'`
      igndvap=`echo ${oi}      | awk '{print $85}'`
      iphen=`echo ${oi}        | awk '{print $86}'`
      iallom=`echo ${oi}       | awk '{print $87}'`
      ibigleaf=`echo ${oi}     | awk '{print $88}'`
      irepro=`echo ${oi}       | awk '{print $89}'`
      treefall=`echo ${oi}     | awk '{print $90}'`
      #------------------------------------------------------------------------------------#


      #----- Define the job name. ---------------------------------------------------------#
      jobname="${desc}-${polyname}"
      #------------------------------------------------------------------------------------#


      #----- Update the R script that controls the run. -----------------------------------#
      /bin/rm -f ${here}/${polyname}/whichrun.r
      /bin/cp -f ${here}/Template/whichrun.r ${here}/${polyname}/whichrun.r
      sed -i s@thispoly@${polyname}@g        ${here}/${polyname}/whichrun.r
      sed -i s@thisqueue@${queue}@g          ${here}/${polyname}/whichrun.r
      sed -i s@pathhere@${here}@g            ${here}/${polyname}/whichrun.r
      sed -i s@paththere@${here}@g           ${here}/${polyname}/whichrun.r
      sed -i s@thisfirst@${yr1sthisto}@g     ${here}/${polyname}/whichrun.r
      sed -i s@thisdthisto@${dthisto}@g      ${here}/${polyname}/whichrun.r
      sed -i s@thismetyeara@${metyeara}@g    ${here}/${polyname}/whichrun.r
      sed -i s@thismetcycle@${metcycle}@g    ${here}/${polyname}/whichrun.r
      sed -i s@thisststcrit@${ststcrit}@g    ${here}/${polyname}/whichrun.r
      sed -i s@thisncycle@${ncycle}@g        ${here}/${polyname}/whichrun.r
      sed -i s@thisyeara@${yeara}@g          ${here}/${polyname}/whichrun.r
      sed -i s@thismontha@${montha}@g        ${here}/${polyname}/whichrun.r
      sed -i s@thisdatea@${datea}@g          ${here}/${polyname}/whichrun.r
      sed -i s@thistimea@${timea}@g          ${here}/${polyname}/whichrun.r
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #    Retrie ve the information on statusrun.txt because it may tell whether the run  #
      # was in steady state or all PFTs went extinct.  These will never start again.       #
      #------------------------------------------------------------------------------------#
      statrun=${here}/${polyname}/statusrun.txt
      if [ -s ${statrun} ]
      then
         yearh=`cat ${statrun}  | awk '{print $2}'`
         monthh=`cat ${statrun} | awk '{print $3}'`
         dateh=`cat ${statrun}  | awk '{print $4}'`
         timeh=`cat ${statrun}  | awk '{print $5}'`
         runt=`cat ${statrun}   | awk '{print $6}'`
         agb=`cat ${statrun}    | awk '{print $7}'`
         bsa=`cat ${statrun}    | awk '{print $8}'`
         lai=`cat ${statrun}    | awk '{print $9}'`
      else
         yearh=${yeara}
         monthh=${montha}
         dateh=${datea}
         timeh=${timea}
         runt='INITIAL'
         agb='NA'
         bsa='NA'
         lai='NA'
      fi
      #------------------------------------------------------------------------------------#




      #------------------------------------------------------------------------------------#
      #    We only check those polygons that may be running, so the check doesn't take too #
      # long.  Once a polygon has reached a steady state / gone extinct / crashed (oh      #
      # well, that happens unfortunately), then its status should remain the same.         #
      #------------------------------------------------------------------------------------#
      if [ ${runt} != 'EXTINCT' ] && [ ${runt} != 'STSTATE' ] && [ ${runt} != 'THE_END' ]
      then
         #----- Call R to check status. ---------------------------------------------------#
         R CMD BATCH ${here}/${polyname}/whichrun.r ${here}/${polyname}/outwhich.txt
         year=`cat ${statrun}  | awk '{print $2}'`
         month=`cat ${statrun} | awk '{print $3}'`
         day=`cat ${statrun}   | awk '{print $4}'`
         hhmm=`cat ${statrun}  | awk '{print $5}'`
         runt=`cat ${statrun}  | awk '{print $6}'`
         agb=`cat ${statrun}   | awk '{print $7}'`
         bsa=`cat ${statrun}   | awk '{print $8}'`
         lai=`cat ${statrun}   | awk '{print $9}'`
      fi
      #------------------------------------------------------------------------------------#




      #------------------------------------------------------------------------------------#
      #    In case the polygon is extinct or it has crashed, we must do some stuff...      #
      #------------------------------------------------------------------------------------#
      if [ ${runt} == 'CRASHED' ] || [ ${runt} == 'STOPPED' ] || [ ${runt} == 'METMISS' ]
      then

         #---------------------------------------------------------------------------------#
         #     Find out the tolerance used in the unfortunate run, then submit the job     #
         # with a tolerance 10 times more strict.                                          #
         #---------------------------------------------------------------------------------#
         ED2IN=${here}'/'${polyname}'/ED2IN'
         #----- Save the tolerance before we delete the file. -----------------------------#
         toler=`grep NL%RK4_TOLERANCE ${ED2IN} | awk '{print $3}'`
         if [ 'x'${toler} == 'xmytoler' -o 'x'${toler} == 'x' ]
         then
            toler=0.01
         fi # [ 'x'${oldtol} == 'xmytoler' -o 'x'${oldtol} == 'x' ]
         #---------------------------------------------------------------------------------#


         if [ ${runt} == 'CRASHED' ]
         then
            echo '  :-( Polygon '${polyname}' has crashed !!!'
            let ncrashed=${ncrashed}+1
            #----- Make the tolerance 10 times smaller. -----------------------------------#
            toler=`${ccc} ${toler}/10`
            echo '      - New tolerance = '${toler}
            /bin/mv ${here}/${polyname}/serial_out.out ${here}/${polyname}/crashed_out.out
            /bin/mv ${here}/${polyname}/serial_out.err ${here}/${polyname}/crashed_out.err
         elif [ ${runt} == 'METMISS' ]
         then
            let nmetmiss=${nmetmiss}+1
            echo '  :-{ Polygon '${polyname}' stopped due to missing met driver !!!'
            /bin/mv ${here}/${polyname}/serial_out.out ${here}/${polyname}/metmiss_out.out
            /bin/mv ${here}/${polyname}/serial_out.err ${here}/${polyname}/metmiss_out.err
         elif [ ${runt} == 'STOPPED' ]
         then
            let nstopped=${nstopped}+1
            echo '  o_0 Polygon '${polyname}' mysteriously stopped!!!'
            /bin/mv ${here}/${polyname}/serial_out.out ${here}/${polyname}/stopped_out.out
            /bin/mv ${here}/${polyname}/serial_out.err ${here}/${polyname}/stopped_out.err
         fi
         #---------------------------------------------------------------------------------#

         if [ ${toler} != '0.0000010' -a ${toler} != '0.0000001' -a ${toler} != '00' ]
         then

            #----- Re-generate ED2IN. -----------------------------------------------------#
            rm -f ${ED2IN} 
            rm -f ${here}/${polyname}/serial_lsf.out

            #----- Copy the Template to the right directory. ------------------------------#
            cp ${here}/Template/ED2IN ${ED2IN}


            #------------------------------------------------------------------------------#
            #     Determine which PFTs to use.                                             #
            #------------------------------------------------------------------------------#
            lonint=${polylon/.*}
            latint=${polylat/.*}
            pfts='1,2,3,4,16'
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #    Find the biometric files.  This has been checked in spawn_poly.sh so they #
            # are correct.  Add a dummy name in case this is not supposed to be a biomass  #
            # initialisation run.                                                          #
            #------------------------------------------------------------------------------#
            case ${isizepft} in
            0)
               #----- Frankeinstein's under storey. ---------------------------------------#
               thisbiomin="${bioinit}/${polyiata}_default."
               ;;
            1)
               #----- No under storey. ----------------------------------------------------#
               thisbiomin="${bioinit}/${polyiata}_nounder."
               ;;
            2)
               #----- Same as default, but with only one grass and one tree. --------------#
               thisbiomin="${bioinit}/${polyiata}_twopft."
               ;;
            *)
               thisbiomin="${bioinit}/${polyiata}_nothing."
               ;;
            esac
            #------------------------------------------------------------------------------#




            #----- Check whether to use SFILIN as restart or history. ---------------------#
            if [ ${runt} == 'INITIAL' ] && [ ${initmode} -eq 6 ]
            then
               thissfilin=${thisbiomin}
            else
               thissfilin=${there}/${polyname}/histo/${polyname}
            fi
            #---------------------------------------------------------------------------------------#


            #----- Update the polygon information on ED2IN. -------------------------------#
            sed -i s@paththere@${there}@g                ${ED2IN}
            sed -i s@myyeara@${thisyeara}@g              ${ED2IN}
            sed -i s@mymontha@${montha}@g                ${ED2IN}
            sed -i s@mydatea@${datea}@g                  ${ED2IN}
            sed -i s@mytimea@${timea}@g                  ${ED2IN}
            sed -i s@myyearz@${thisyearz}@g              ${ED2IN}
            sed -i s@mymonthz@${monthz}@g                ${ED2IN}
            sed -i s@mydatez@${datez}@g                  ${ED2IN}
            sed -i s@mytimez@${timez}@g                  ${ED2IN}
            sed -i s@mydtlsm@${dtlsm}@g                  ${ED2IN}
            sed -i s@thispoly@${polyname}@g              ${ED2IN}
            sed -i s@plonflag@${polylon}@g               ${ED2IN}
            sed -i s@platflag@${polylat}@g               ${ED2IN}
            sed -i s@timehhhh@${time}@g                  ${ED2IN}
            sed -i s@datehhhh@${date}@g                  ${ED2IN}
            sed -i s@monthhhh@${month}@g                 ${ED2IN}
            sed -i s@yearhhhh@${year}@g                  ${ED2IN}
            sed -i s@myinitmode@${initmode}@g            ${ED2IN}
            sed -i s@mysfilin@${thissfilin}@g            ${ED2IN}
            sed -i s@mytrees@${pfts}@g                   ${ED2IN}
            sed -i s@mycrop@${crop}@g                    ${ED2IN}
            sed -i s@myplantation@${plantation}@g        ${ED2IN}
            sed -i s@myiphen@${iphen}@g                  ${ED2IN}
            sed -i s@myallom@${iallom}@g                 ${ED2IN}
            sed -i s@myisoilflg@${polyisoil}@g           ${ED2IN}
            sed -i s@mynslcon@${polyntext}@g             ${ED2IN}
            sed -i s@myslxsand@${polysand}@g             ${ED2IN}
            sed -i s@myslxclay@${polyclay}@g             ${ED2IN}
            sed -i s@mysoilbc@${polysoilbc}@g            ${ED2IN}
            sed -i s@mysldrain@${polysldrain}@g          ${ED2IN}
            sed -i s@mysoilcol@${polycol}@g              ${ED2IN}
            sed -i s@mynzg@${polynzg}@g                  ${ED2IN}
            sed -i s@mymetdriverdb@${metdriverdb}@g      ${ED2IN}
            sed -i s@mymetcyc1@${metcyc1}@g              ${ED2IN}
            sed -i s@mymetcycf@${metcycf}@g              ${ED2IN}
            sed -i s@mytoler@${toler}@g                  ${ED2IN}
            sed -i s@RUNFLAG@HISTORY@g                   ${ED2IN}
            sed -i s@myvmfactc3@${vmfactc3}@g            ${ED2IN}
            sed -i s@myvmfactc4@${vmfactc4}@g            ${ED2IN}
            sed -i s@mymphototrc3@${mphototrc3}@g        ${ED2IN}
            sed -i s@mymphototec3@${mphototec3}@g        ${ED2IN}
            sed -i s@mymphotoc4@${mphotoc4}@g            ${ED2IN}
            sed -i s@mybphotoblc3@${bphotoblc3}@g        ${ED2IN}
            sed -i s@mybphotonlc3@${bphotonlc3}@g        ${ED2IN}
            sed -i s@mybphotoc4@${bphotoc4}@g            ${ED2IN}
            sed -i s@mykwgrass@${kwgrass}@g              ${ED2IN}
            sed -i s@mykwtree@${kwtree}@g                ${ED2IN}
            sed -i s@mygammac3@${gammac3}@g              ${ED2IN}
            sed -i s@mygammac4@${gammac4}@g              ${ED2IN}
            sed -i s@myd0grass@${d0grass}@g              ${ED2IN}
            sed -i s@myd0tree@${d0tree}@g                ${ED2IN}
            sed -i s@myalphac3@${alphac3}@g              ${ED2IN}
            sed -i s@myalphac4@${alphac4}@g              ${ED2IN}
            sed -i s@myklowco2@${klowco2}@g              ${ED2IN}
            sed -i s@mydecomp@${decomp}@g                ${ED2IN}
            sed -i s@myrrffact@${rrffact}@g              ${ED2IN}
            sed -i s@mygrowthresp@${growthresp}@g        ${ED2IN}
            sed -i s@mylwidthgrass@${lwidthgrass}@g      ${ED2IN}
            sed -i s@mylwidthbltree@${lwidthbltree}@g    ${ED2IN}
            sed -i s@mylwidthnltree@${lwidthnltree}@g    ${ED2IN}
            sed -i s@myq10c3@${q10c3}@g                  ${ED2IN}
            sed -i s@myq10c4@${q10c4}@g                  ${ED2IN}
            sed -i s@myh2olimit@${h2olimit}@g            ${ED2IN}
            sed -i s@mymortscheme@${imortscheme}@g       ${ED2IN}
            sed -i s@myddmortconst@${ddmortconst}@g      ${ED2IN}
            sed -i s@mysfclyrm@${isfclyrm}@g             ${ED2IN}
            sed -i s@myicanturb@${icanturb}@g            ${ED2IN}
            sed -i s@myatmco2@${atmco2}@g                ${ED2IN}
            sed -i s@mythcrit@${thcrit}@g                ${ED2IN}
            sed -i s@mysmfire@${smfire}@g                ${ED2IN}
            sed -i s@myfire@${ifire}@g                   ${ED2IN}
            sed -i s@myfuel@${fireparm}@g                ${ED2IN}
            sed -i s@mymetavg@${imetavg}@g               ${ED2IN}
            sed -i s@mypercol@${ipercol}@g               ${ED2IN}
            sed -i s@myrunoff@${runoff}@g                ${ED2IN}
            sed -i s@mymetrad@${imetrad}@g               ${ED2IN}
            sed -i s@mybranch@${ibranch}@g               ${ED2IN}
            sed -i s@mycanrad@${icanrad}@g               ${ED2IN}
            sed -i s@mycrown@${crown}@g                  ${ED2IN}
            sed -i s@myltransvis@${ltransvis}@g          ${ED2IN}
            sed -i s@myltransnir@${ltransnir}@g          ${ED2IN}
            sed -i s@mylreflectvis@${lreflectvis}@g      ${ED2IN}
            sed -i s@mylreflectnir@${lreflectnir}@g      ${ED2IN}
            sed -i s@myorienttree@${orienttree}@g        ${ED2IN}
            sed -i s@myorientgrass@${orientgrass}@g      ${ED2IN}
            sed -i s@myclumptree@${clumptree}@g          ${ED2IN}
            sed -i s@myclumpgrass@${clumpgrass}@g        ${ED2IN}
            sed -i s@myvegtdyn@${ivegtdyn}@g             ${ED2IN}
            sed -i s@mybigleaf@${ibigleaf}@g             ${ED2IN}
            sed -i s@myrepro@${irepro}@g                 ${ED2IN}
            sed -i s@myubmin@${ubmin}@g                  ${ED2IN}
            sed -i s@myugbmin@${ugbmin}@g                ${ED2IN}
            sed -i s@myustmin@${ustmin}@g                ${ED2IN}
            sed -i s@mygamm@${gamm}@g                    ${ED2IN}
            sed -i s@mygamh@${gamh}@g                    ${ED2IN}
            sed -i s@mytprandtl@${tprandtl}@g            ${ED2IN}
            sed -i s@myribmax@${ribmax}@g                ${ED2IN}
            sed -i s@mygndvap@${igndvap}@g               ${ED2IN}
            sed -i s@mydtcensus@${dtcensus}@g            ${ED2IN}
            sed -i s@myyr1stcensus@${yr1stcensus}@g      ${ED2IN}
            sed -i s@mymon1stcensus@${mon1stcensus}@g    ${ED2IN}
            sed -i s@myminrecruitdbh@${minrecruitdbh}@g  ${ED2IN}
            sed -i s@mytreefall@${treefall}@g            ${ED2IN}
            #------ Soil variables. -------------------------------------------------------#
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
            #------------------------------------------------------------------------------#






            #------------------------------------------------------------------------------#
            #      A crashed polygon will inevitably go to the long_serial queue.  Check   #
            # whether this is a head of a unrestricted parallel run.  If it is, then call  #
            # this something else, otherwise use the polygon name.  In any case, re-write  #
            # the srun.sh file.                                                            #
            #------------------------------------------------------------------------------#
            srun=${here}'/'${polyname}'/srun.sh'

            unpahead=`bjobs -w -J ${jobname} 2> /dev/null | wc -l`
            if [ ${unpahead} -eq 0 ]
            then
               thisname=${jobname}
            else
               thisname="app-${thisname}"
            fi
            #------------------------------------------------------------------------------#



            #----- Get the waiting time used before. --------------------------------------#
            thiswait=`grep "Sleep time is" ${polyname}/srun.sh | awk '{print $5}'`
            #------------------------------------------------------------------------------#


            #----- Copy the data set again. -----------------------------------------------#
            /bin/rm -f ${srun}
            /bin/rm -f ${here}/${polyname}/skipper.txt
            /bin/cp ${here}/Template/srun.sh ${here}/${polyname}/srun.sh
            #------------------------------------------------------------------------------#


            #----- The new queue is by default the long_serial. ---------------------------#
            newqueue='long_serial'
            #------------------------------------------------------------------------------#



            #----- Change some settings in srun.sh. ---------------------------------------#
            sed -i s@"jobname='thisjob'"@"jobname='${thisname}'"@g ${srun}
            sed -i s@pathhere@${here}@g      ${srun}
            sed -i s@paththere@${here}@g     ${srun}
            sed -i s@thispoly@${polyname}@g  ${srun}
            sed -i s@thisdesc@${desc}@g      ${srun}
            sed -i s@thisqueue@${newqueue}@g ${srun}
            sed -i s@zzzzzzzz@${wtime}@g     ${srun}
            sed -i s@myorder@${poly}@g       ${srun}
            #------------------------------------------------------------------------------#

         else
            #----- Give up on this node. --------------------------------------------------#
            echo '   :-/ Tolerance is tiny and it still does not work.  Giving up...'
            runt='THE_END'
            #------------------------------------------------------------------------------#
         fi # [ ${toler} != '0.0000010' -a ${toler} != '0.0000001' -a ${toler} != '00' ]
         #---------------------------------------------------------------------------------#

      elif [ ${runt} == 'THE_END' ]
      then 
         echo '  :-) Polygon '${polyname}' finished !!!'

      elif [ ${runt} == 'EXTINCT' -o ${runt} == 'STSTATE' ] && 
           [ ${queue} != 'unrestricted_parallel' ]
      then
         running=`bjobs -w -J ${jobname} 2> /dev/null | grep RUN | wc -l`
         if [ ${running} -gt 0 ]
         then
            deathrow="${deathrow} ${jobname}"
            if [ ${runt} == 'EXTINCT' ]
            then
               let newextinct=${newextinct}+1
            elif [ ${runt} == 'STSTATE' ]
            then
               let newststate=${newststate}+1
            fi
         fi
      fi # [ ${runt} == 'CRASHED' ]
      #------------------------------------------------------------------------------------#




      #------------------------------------------------------------------------------------#
      #    Last, we check whether the job is still running or not.  In case it is not but  #
      # it should, update the ED2IN with the latest history and re-submit the job.         #
      #------------------------------------------------------------------------------------#
      if [ ${runt} != 'THE_END' ] && [ ${runt} != 'STSTATE' ] && [ ${runt} != 'EXTINCT' ]
      then
         if [ ${queue} == 'unrestricted_parallel' ]
         then
            nlist=`grep ${polyname} ${here}/*/callunpa.sh | wc -l`
            if [ ${nlist} -eq 0 ]
            then
               running=0
            else
               j=0
               running=0
               while [ ${j} -lt ${nlist} ] && [ ${running} -eq 0 ]
               do
                  let j=${j}+1
                  thisline=`grep ${polyname} ${here}/*/callunpa.sh | head -${j} | tail -1`
                  parent=`echo ${thisline} | awk '{print substr($1,1,8)}'`
                  parjob="${desc}-${parent}"
                  running=`bjobs -w -q ${queue} -J ${parjob} 2> /dev/null | tail -1 | wc -l`
               done
            fi
         else
            running=`bjobs -w -J ${jobname} 2> /dev/null | tail -1 | wc -l`
         fi

         if [ ${running} -eq 0 ]
         then
            echo "${polyname} is missing.  Re-submitting..." >> ${situation}

            let nresubmit=${nresubmit}+1

            #------------------------------------------------------------------------------#
            #     Job was kicked out, submit it again.                                     #
            #------------------------------------------------------------------------------#
            ED2IN=${here}'/'${polyname}'/ED2IN'
            #----- Save the tolerance before we delete the file. --------------------------#
            toler=`grep NL%RK4_TOLERANCE ${ED2IN} | awk '{print $3}'`
            #------------------------------------------------------------------------------#
            #   This should never happen, but...                                           #
            #------------------------------------------------------------------------------#
            if [ 'x'${oldtol} == 'xmytoler' -o 'x'${oldtol} == 'x' ]
            then
               toler=0.01
            else
               toler=${oldtol}
            fi # [ 'x'${oldtol} == 'xmytoler' -o 'x'${oldtol} == 'x' ]
            #------------------------------------------------------------------------------#


            #----- Re-generate ED2IN. -----------------------------------------------------#
            rm -f ${ED2IN} 
            rm -f ${here}/${polyname}/serial_out.out
            rm -f ${here}/${polyname}/serial_lsf.out
            rm -f ${here}/${polyname}/serial_out.err

            #----- Copy the Template to the right directory. ------------------------------#
            cp ${here}/Template/ED2IN ${ED2IN}

            #------------------------------------------------------------------------------#
            #     Determine which PFTs to use.                                             #
            #------------------------------------------------------------------------------#
            lonint=${polylon/.*}
            latint=${polylat/.*}
            pfts='1,2,3,4,16'
            #------------------------------------------------------------------------------#




            #----- Check whether to use SFILIN as restart or history. ---------------------#
            if [ ${runt} == 'INITIAL' ] && [ ${initmode} -eq 6 ]
            then
               thissfilin=${thisbiomin}
            else
               thissfilin=${there}/${polyname}/histo/${polyname}
            fi
            #---------------------------------------------------------------------------------------#

            #----- Update the polygon information on ED2IN. -------------------------------#
            sed -i s@paththere@${there}@g                ${ED2IN}
            sed -i s@myyeara@${thisyeara}@g              ${ED2IN}
            sed -i s@mymontha@${montha}@g                ${ED2IN}
            sed -i s@mydatea@${datea}@g                  ${ED2IN}
            sed -i s@mytimea@${timea}@g                  ${ED2IN}
            sed -i s@myyearz@${thisyearz}@g              ${ED2IN}
            sed -i s@mymonthz@${monthz}@g                ${ED2IN}
            sed -i s@mydatez@${datez}@g                  ${ED2IN}
            sed -i s@mytimez@${timez}@g                  ${ED2IN}
            sed -i s@mydtlsm@${dtlsm}@g                  ${ED2IN}
            sed -i s@thispoly@${polyname}@g              ${ED2IN}
            sed -i s@plonflag@${polylon}@g               ${ED2IN}
            sed -i s@platflag@${polylat}@g               ${ED2IN}
            sed -i s@timehhhh@${time}@g                  ${ED2IN}
            sed -i s@datehhhh@${date}@g                  ${ED2IN}
            sed -i s@monthhhh@${month}@g                 ${ED2IN}
            sed -i s@yearhhhh@${year}@g                  ${ED2IN}
            sed -i s@myinitmode@${initmode}@g            ${ED2IN}
            sed -i s@mysfilin@${thissfilin}@g            ${ED2IN}
            sed -i s@mytrees@${pfts}@g                   ${ED2IN}
            sed -i s@mycrop@${crop}@g                    ${ED2IN}
            sed -i s@myplantation@${plantation}@g        ${ED2IN}
            sed -i s@myiphen@${iphen}@g                  ${ED2IN}
            sed -i s@myallom@${iallom}@g                 ${ED2IN}
            sed -i s@myisoilflg@${polyisoil}@g           ${ED2IN}
            sed -i s@mynslcon@${polyntext}@g             ${ED2IN}
            sed -i s@myslxsand@${polysand}@g             ${ED2IN}
            sed -i s@myslxclay@${polyclay}@g             ${ED2IN}
            sed -i s@mysoilbc@${polysoilbc}@g            ${ED2IN}
            sed -i s@mysldrain@${polysldrain}@g          ${ED2IN}
            sed -i s@mysoilcol@${polycol}@g              ${ED2IN}
            sed -i s@mynzg@${polynzg}@g                  ${ED2IN}
            sed -i s@mymetdriverdb@${metdriverdb}@g      ${ED2IN}
            sed -i s@mymetcyc1@${metcyc1}@g              ${ED2IN}
            sed -i s@mymetcycf@${metcycf}@g              ${ED2IN}
            sed -i s@mytoler@${toler}@g                  ${ED2IN}
            sed -i s@myvmfactc3@${vmfactc3}@g            ${ED2IN}
            sed -i s@myvmfactc4@${vmfactc4}@g            ${ED2IN}
            sed -i s@mymphototrc3@${mphototrc3}@g        ${ED2IN}
            sed -i s@mymphototec3@${mphototec3}@g        ${ED2IN}
            sed -i s@mymphotoc4@${mphotoc4}@g            ${ED2IN}
            sed -i s@mybphotoblc3@${bphotoblc3}@g        ${ED2IN}
            sed -i s@mybphotonlc3@${bphotonlc3}@g        ${ED2IN}
            sed -i s@mybphotoc4@${bphotoc4}@g            ${ED2IN}
            sed -i s@mykwgrass@${kwgrass}@g              ${ED2IN}
            sed -i s@mykwtree@${kwtree}@g                ${ED2IN}
            sed -i s@mygammac3@${gammac3}@g              ${ED2IN}
            sed -i s@mygammac4@${gammac4}@g              ${ED2IN}
            sed -i s@myd0grass@${d0grass}@g              ${ED2IN}
            sed -i s@myd0tree@${d0tree}@g                ${ED2IN}
            sed -i s@myalphac3@${alphac3}@g              ${ED2IN}
            sed -i s@myalphac4@${alphac4}@g              ${ED2IN}
            sed -i s@myklowco2@${klowco2}@g              ${ED2IN}
            sed -i s@mydecomp@${decomp}@g                ${ED2IN}
            sed -i s@myrrffact@${rrffact}@g              ${ED2IN}
            sed -i s@mygrowthresp@${growthresp}@g        ${ED2IN}
            sed -i s@mylwidthgrass@${lwidthgrass}@g      ${ED2IN}
            sed -i s@mylwidthbltree@${lwidthbltree}@g    ${ED2IN}
            sed -i s@mylwidthnltree@${lwidthnltree}@g    ${ED2IN}
            sed -i s@myq10c3@${q10c3}@g                  ${ED2IN}
            sed -i s@myq10c4@${q10c4}@g                  ${ED2IN}
            sed -i s@myh2olimit@${h2olimit}@g            ${ED2IN}
            sed -i s@mymortscheme@${imortscheme}@g       ${ED2IN}
            sed -i s@myddmortconst@${ddmortconst}@g      ${ED2IN}
            sed -i s@mysfclyrm@${isfclyrm}@g             ${ED2IN}
            sed -i s@myicanturb@${icanturb}@g            ${ED2IN}
            sed -i s@myatmco2@${atmco2}@g                ${ED2IN}
            sed -i s@mythcrit@${thcrit}@g                ${ED2IN}
            sed -i s@mysmfire@${smfire}@g                ${ED2IN}
            sed -i s@myfire@${ifire}@g                   ${ED2IN}
            sed -i s@myfuel@${fireparm}@g                ${ED2IN}
            sed -i s@mymetavg@${imetavg}@g               ${ED2IN}
            sed -i s@mypercol@${ipercol}@g               ${ED2IN}
            sed -i s@myrunoff@${runoff}@g                ${ED2IN}
            sed -i s@mymetrad@${imetrad}@g               ${ED2IN}
            sed -i s@mybranch@${ibranch}@g               ${ED2IN}
            sed -i s@mycanrad@${icanrad}@g               ${ED2IN}
            sed -i s@mycrown@${crown}@g                  ${ED2IN}
            sed -i s@myltransvis@${ltransvis}@g          ${ED2IN}
            sed -i s@myltransnir@${ltransnir}@g          ${ED2IN}
            sed -i s@mylreflectvis@${lreflectvis}@g      ${ED2IN}
            sed -i s@mylreflectnir@${lreflectnir}@g      ${ED2IN}
            sed -i s@myorienttree@${orienttree}@g        ${ED2IN}
            sed -i s@myorientgrass@${orientgrass}@g      ${ED2IN}
            sed -i s@myclumptree@${clumptree}@g          ${ED2IN}
            sed -i s@myclumpgrass@${clumpgrass}@g        ${ED2IN}
            sed -i s@myvegtdyn@${ivegtdyn}@g             ${ED2IN}
            sed -i s@mybigleaf@${ibigleaf}@g             ${ED2IN}
            sed -i s@myrepro@${irepro}@g                 ${ED2IN}
            sed -i s@myubmin@${ubmin}@g                  ${ED2IN}
            sed -i s@myugbmin@${ugbmin}@g                ${ED2IN}
            sed -i s@myustmin@${ustmin}@g                ${ED2IN}
            sed -i s@mygamm@${gamm}@g                    ${ED2IN}
            sed -i s@mygamh@${gamh}@g                    ${ED2IN}
            sed -i s@mytprandtl@${tprandtl}@g            ${ED2IN}
            sed -i s@myribmax@${ribmax}@g                ${ED2IN}
            sed -i s@mygndvap@${igndvap}@g               ${ED2IN}
            sed -i s@mydtcensus@${dtcensus}@g            ${ED2IN}
            sed -i s@myyr1stcensus@${yr1stcensus}@g      ${ED2IN}
            sed -i s@mymon1stcensus@${mon1stcensus}@g    ${ED2IN}
            sed -i s@myminrecruitdbh@${minrecruitdbh}@g  ${ED2IN}
            sed -i s@mytreefall@${treefall}@g            ${ED2IN}
            #------ Soil variables. -------------------------------------------------------#
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
            #----- Change run type. -------------------------------------------------------#
            if [ ${runt} == 'INITIAL' ]
            then
               sed -i s@RUNFLAG@${runt}@g      ${ED2IN}
           else
               sed -i s@RUNFLAG@HISTORY@g      ${ED2IN}
            fi
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #      A halted polygon will inevitably go to the long_serial queue.  Re-write #
            # the srun.sh file.                                                            #
            #------------------------------------------------------------------------------#
            srun=${here}'/'${polyname}'/srun.sh'

            #----- Get the waiting time used before. --------------------------------------#
            thiswait=`grep 'Sleep time is' ${polyname}/srun.sh | awk '{print $5}'`
            if [ 'x'${thiswait} == 'x' ]
            then
               let wtime=${poly}%8
               let wtime=${wtime}*20
               nudge=`date +%S`
               if [ ${nudge} -lt 10 ]
               then 
                  nudge=`echo ${nudge} | awk '{print substr($1,2,1)}'`
               fi
               let nudge=${nudge}%15
               let wtime=${wtime}+${nudge}
               let wtime=${wtime}+2
            fi
            /bin/rm -f ${srun}
            /bin/rm -f ${here}/${polyname}/skipper.txt

            cp ${here}/Template/srun.sh ${here}/${polyname}/srun.sh

            #----- The new queue is by default the long_serial. ---------------------------#
            newqueue='long_serial'

            #----- Change some settings in srun.sh. ---------------------------------------#
            sed -i s@pathhere@${here}@g      ${srun}
            sed -i s@paththere@${here}@g     ${srun}
            sed -i s@thispoly@${polyname}@g  ${srun}
            sed -i s@thisqueue@${newqueue}@g ${srun}
            sed -i s@thisdesc@${desc}@g      ${srun}
            sed -i s@zzzzzzzz@${wtime}@g     ${srun}
            sed -i s@myorder@${poly}@g       ${srun}

            #----- Re-submit the job to long_serial queue. --------------------------------#
            ${here}'/'${polyname}'/srun.sh'

            #----- Switch the status to pending. ------------------------------------------#
            year=${yeara}
            month=${montha}
            day=${datea}
            hhmm=${houra}
            runt='INITIAL'
            agb='NA'
            bsa='NA'
            lai='NA'
         else
            echo "${polyname} is running/pending..." >> ${situation}
         fi # [ ${running} -eq 0 ]
      else
         echo "${polyname} status is ${runt}.  No need to resubmit." >> ${situation}
      fi # [ ${runt} == 'HISTORY' -a ${queue} == 'long_serial' ]
      #------------------------------------------------------------------------------------#




      #------------------------------------------------------------------------------------#
      #     Test whether the job can be submitted.                                         #
      #------------------------------------------------------------------------------------#
      if [ ${runt} == 'INITIAL' ]
      then
         if [ ${queue} == 'unrestricted_parallel' ]
         then
            nlist=`grep ${polyname} ${here}/*/callunpa.sh | wc -l`
            if [ ${nlist} -eq 0 ]
            then
               running=0
            else
               j=0
               running=0
               while [ ${j} -lt ${nlist} ] && [ ${running} -eq 0 ]
               do
                  let j=${j}+1
                  thisline=`grep ${polyname} ${here}/*/callunpa.sh | head -${j} | tail -1`
                  parent=`echo ${thisline} | awk '{print substr($1,1,8)}'`
                  parjob="${desc}-${parent}"
                  running=`bjobs -w -q ${queue} -J ${parjob} 2> /dev/null | tail -1 | wc -l`
               done
            fi
         else
            running=`bjobs -w -J ${jobname} 2> /dev/null | tail -1 | wc -l`
         fi
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Check whether to start a new simulation or wait.                            #
         #---------------------------------------------------------------------------------#
         if [ "x${waitpot}" == "xy" ] || [ "x${waitpot}" == "xy" ]
         then
            readytostart=`/bin/ls -1 ${potveg}/${oldname}*h5 2> /dev/null | wc -l`
         else
            readytostart=1
         fi
         if [ ${running} -eq 0 ] && [ ${readytostart} -eq 1 ]
         then
            let nstart=${nstart}+1
            echo " >>> Initial submission of polygon ${polyname}!!!!"

            #----- Re-submit the job to long_serial queue. --------------------------------#
            ${here}'/'${polyname}'/srun.sh'

         elif [ ${readytostart} -gt 1 ]
         then
            echo ":-S Something strange with the restart for polygon ${polyname}..."
         fi
         #---------------------------------------------------------------------------------#
      fi
      #------------------------------------------------------------------------------------#



      #----- Add the speed test for the polygons that ran at least two cycles. ------------#
      let dyear=${year}-${yeara}+1
      if [ ${dyear} -gt 2 ]
      then
         premier=${here}/${polyname}/histo/${polyname}-S-${yeara}-01-01-000000-g01.h5
         dernier=${here}/${polyname}/histo/${polyname}-S-${year}-01-01-000000-g01.h5

         #----- Find time it took so far. -------------------------------------------------#
         whenh1st=`stat -c %Z ${premier} | awk '{print $1}'`
         whenhlast=`stat -c %Z ${dernier} | awk '{print $1}'`
         let dwhen=${whenhlast}-${whenh1st}
         dwhen=`${ccc} ${dwhen}/86400.`

         #----- Find the number of years. -------------------------------------------------#
         let yearfull=${yearz}-${yeara}
         let yearsofar=${year}-${yeara}

         #----- Estimate the number of years until the last year. -------------------------#
         whenend=`${ccc} ${dwhen}*${yearfull}/${yearsofar}`
      else
         whenend='NA'
      fi
      #---------------------------------------------------------------------------------------#



      #----- Write polygon check into a single table. -------------------------------------#
      output="${polyname} ${polylon} ${polylat} ${year} ${month} ${day} ${hhmm}"
      output="${output} ${runt} ${agb} ${bsa} ${lai} ${whenend}"
      echo ${output} >> ${outcheck}
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #    Update joborder.txt file.                                                       #
      #------------------------------------------------------------------------------------#
      if [ ${runt} != 'THE_END' ] && [ ${runt} != 'STSTATE' ] && [ ${runt} != 'EXTINCT' ]
      then
         polyIATA=`echo ${polyiata} | tr '[:lower:]' '[:upper:]'`
         altname=`echo ${jobname} | sed s@${polyiata}@${polyIATA}@g`
         polyqueue=`bjobs -w -J ${jobname}  2> /dev/null | tail -1 | awk '{print $4}'`
         altequeue=`bjobs -w -J ${altname}  2> /dev/null | tail -1 | awk '{print $4}'`

         if [ -s ${here}/${polyname}/skipper.txt ]
         then
            #----- No job with this name, it is an unrestricted parallel job. -------------#
            newqueue='unrestricted_parallel'
            success='TRUE'
         elif [ 'x'${polyqueue} != 'x' ] && [ 'x'${altequeue} == 'x' ]
         then
            #----- "Normal" job, with regular name. ---------------------------------------#
            newqueue=${polyqueue}
            success='TRUE'
         elif [ 'x'${altequeue} != 'x' ]
         then
            #------------------------------------------------------------------------------#
            #    This will happen only when the job used to be the "head" of an un-        #
            # restricted parallel job, but the head crashed.                               #
            #------------------------------------------------------------------------------#
            newqueue=${altequeue}
            success='TRUE'
         else
            #------------------------------------------------------------------------------#
            #    Job was interrupted.  We resubmit the job in long_serial.                 #
            #------------------------------------------------------------------------------#
            success='FALSE'
         fi
         if [ ${success} == 'TRUE' ] 
         then
            oldline=${oi}
            newline=`echo ${oldline} | sed s@${queue}@${newqueue}@g`
            sed -i s@"${oldline}"@"${newline}"@g ${joborder}
         fi
      fi
      
      #------------------------------------------------------------------------------------#
      #    Check whether the steady-state restart file is copied to the destination        #
      # directory.                                                                         #
      #------------------------------------------------------------------------------------#
      if [ "x${copyrestart}" == "xy" ] || [ "x${copyrestart}" == "xY" ]
      then
         if [ ${runt} == 'STSTATE' -o ${runt} == 'THE_END' -o ${runt} == 'EXTINCT' ]
         then
            nrestart=`/bin/ls -1 ${restart}/${polyname}-S-*h5 2> /dev/null | wc -l`
            #----- Check whether there is any file in the restart directory ---------------#
            if [ ${nrestart} -eq 0 ]
            then
               lastht=`/bin/ls -1 ${here}/${polyname}/histo/*-S-*h5 2> /dev/null | tail -1`
               baseht=`basename ${lastht}`
               destht=${restart}'/'${baseht}
               echo ' - Copying file '${baseht}' to the restart directory...'
               echo ' - LASTHISTO = '${lastht}
               echo ' - DESTHISTO = '${destht}
               /bin/cp -uv ${lastht} ${destht}
            fi
            #------------------------------------------------------------------------------#
         fi
         #---------------------------------------------------------------------------------#
      fi
      #------------------------------------------------------------------------------------#
   done
   #---------------------------------------------------------------------------------------#
fi # end if
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#     Check whether to fix the queues in joborder.                                         #
#------------------------------------------------------------------------------------------#
if [ ${fixqueue} -eq 1 ]
then
   while [ ${ff} -lt ${npolys} ]
   do
      let ff=${ff}+1
      let line=${ff}+3
      let rem=${ff}%${printevery}
      
      if [ ${rem} -eq 0 -o ${ff} -eq ${npolys} ]
      then
         echo ' - Fixed queues for '${ff}' polygons so far...'
      fi

      #------------------------------------------------------------------------------------#
      #      Read the polyth line of the polygon list.  There must be smarter ways of do-  #
      # ing this, but this works.  Here we obtain the polygon name, and its longitude and  #
      # latitude.                                                                          #
      #------------------------------------------------------------------------------------#
      oi=`head -${line} ${lonlat} | tail -1`
      polyname=`echo ${oi}     | awk '{print $1 }'`
      polyiata=`echo ${oi}     | awk '{print $2 }'`
      polylon=`echo ${oi}      | awk '{print $3 }'`
      polylat=`echo ${oi}      | awk '{print $4 }'`
      yeara=`echo ${oi}        | awk '{print $5 }'`
      montha=`echo ${oi}       | awk '{print $6 }'`
      datea=`echo ${oi}        | awk '{print $7 }'`
      timea=`echo ${oi}        | awk '{print $8 }'`
      yearz=`echo ${oi}        | awk '{print $9 }'`
      monthz=`echo ${oi}       | awk '{print $10}'`
      datez=`echo ${oi}        | awk '{print $11}'`
      timez=`echo ${oi}        | awk '{print $12}'`
      initmode=`echo ${oi}     | awk '{print $13}'`
      iscenario=`echo ${oi}    | awk '{print $14}'`
      isizepft=`echo ${oi}     | awk '{print $15}'`
      polyisoil=`echo ${oi}    | awk '{print $16}'`
      polyntext=`echo ${oi}    | awk '{print $17}'`
      polysand=`echo ${oi}     | awk '{print $18}'`
      polyclay=`echo ${oi}     | awk '{print $19}'`
      polydepth=`echo ${oi}    | awk '{print $20}'`
      polysoilbc=`echo ${oi}   | awk '{print $21}'`
      polysldrain=`echo ${oi}  | awk '{print $22}'`
      polycol=`echo ${oi}      | awk '{print $23}'`
      slzres=`echo ${oi}       | awk '{print $24}'`
      queue=`echo ${oi}        | awk '{print $25}'`
      metdriver=`echo ${oi}    | awk '{print $26}'`
      dtlsm=`echo ${oi}        | awk '{print $27}'`
      vmfactc3=`echo ${oi}     | awk '{print $28}'`
      vmfactc4=`echo ${oi}     | awk '{print $29}'`
      mphototrc3=`echo ${oi}   | awk '{print $30}'`
      mphototec3=`echo ${oi}   | awk '{print $31}'`
      mphotoc4=`echo ${oi}     | awk '{print $32}'`
      bphotoblc3=`echo ${oi}   | awk '{print $33}'`
      bphotonlc3=`echo ${oi}   | awk '{print $34}'`
      bphotoc4=`echo ${oi}     | awk '{print $35}'`
      kwgrass=`echo ${oi}      | awk '{print $36}'`
      kwtree=`echo ${oi}       | awk '{print $37}'`
      gammac3=`echo ${oi}      | awk '{print $38}'`
      gammac4=`echo ${oi}      | awk '{print $39}'`
      d0grass=`echo ${oi}      | awk '{print $40}'`
      d0tree=`echo ${oi}       | awk '{print $41}'`
      alphac3=`echo ${oi}      | awk '{print $42}'`
      alphac4=`echo ${oi}      | awk '{print $43}'`
      klowco2=`echo ${oi}      | awk '{print $44}'`
      decomp=`echo ${oi}       | awk '{print $45}'`
      rrffact=`echo ${oi}      | awk '{print $46}'`
      growthresp=`echo ${oi}   | awk '{print $47}'`
      lwidthgrass=`echo ${oi}  | awk '{print $48}'`
      lwidthbltree=`echo ${oi} | awk '{print $49}'`
      lwidthnltree=`echo ${oi} | awk '{print $50}'`
      q10c3=`echo ${oi}        | awk '{print $51}'`
      q10c4=`echo ${oi}        | awk '{print $52}'`
      h2olimit=`echo ${oi}     | awk '{print $53}'`
      imortscheme=`echo ${oi}  | awk '{print $54}'`
      ddmortconst=`echo ${oi}  | awk '{print $55}'`
      isfclyrm=`echo ${oi}     | awk '{print $56}'`
      icanturb=`echo ${oi}     | awk '{print $57}'`
      ubmin=`echo ${oi}        | awk '{print $58}'`
      ugbmin=`echo ${oi}       | awk '{print $59}'`
      ustmin=`echo ${oi}       | awk '{print $60}'`
      gamm=`echo ${oi}         | awk '{print $61}'`
      gamh=`echo ${oi}         | awk '{print $62}'`
      tprandtl=`echo ${oi}     | awk '{print $63}'`
      ribmax=`echo ${oi}       | awk '{print $64}'`
      atmco2=`echo ${oi}       | awk '{print $65}'`
      thcrit=`echo ${oi}       | awk '{print $66}'`
      smfire=`echo ${oi}       | awk '{print $67}'`
      ifire=`echo ${oi}        | awk '{print $68}'`
      fireparm=`echo ${oi}     | awk '{print $69}'`
      ipercol=`echo ${oi}      | awk '{print $70}'`
      runoff=`echo ${oi}       | awk '{print $71}'`
      imetrad=`echo ${oi}      | awk '{print $72}'`
      ibranch=`echo ${oi}      | awk '{print $73}'`
      icanrad=`echo ${oi}      | awk '{print $74}'`
      crown=`echo   ${oi}      | awk '{print $75}'`
      ltransvis=`echo ${oi}    | awk '{print $76}'`
      lreflectvis=`echo ${oi}  | awk '{print $77}'`
      ltransnir=`echo ${oi}    | awk '{print $78}'`
      lreflectnir=`echo ${oi}  | awk '{print $79}'`
      orienttree=`echo ${oi}   | awk '{print $80}'`
      orientgrass=`echo ${oi}  | awk '{print $81}'`
      clumptree=`echo ${oi}    | awk '{print $82}'`
      clumpgrass=`echo ${oi}   | awk '{print $83}'`
      ivegtdyn=`echo ${oi}     | awk '{print $84}'`
      igndvap=`echo ${oi}      | awk '{print $85}'`
      iphen=`echo ${oi}        | awk '{print $86}'`
      iallom=`echo ${oi}       | awk '{print $87}'`
      ibigleaf=`echo ${oi}     | awk '{print $88}'`
      irepro=`echo ${oi}       | awk '{print $89}'`
      treefall=`echo ${oi}     | awk '{print $90}'`
      #------------------------------------------------------------------------------------#


      polyIATA=`echo ${polyiata} | tr '[:lower:]' '[:upper:]'`
      jobname="${desc}-${polyname}"
      altname=`echo ${jobname} | sed s@${polyiata}@${polyIATA}@g`

      jobqueue=`bjobs -w -J ${jobname} 2> /dev/null | tail -1 | awk '{print $4}'`
      altqueue=`bjobs -w -J ${altname} 2> /dev/null | tail -1 | awk '{print $4}'`


      #----- Update the queue in joborder.txt. --------------------------------------------#
      if [ 'x'${jobqueue} != 'x' ]
      then
         newqueue=${jobqueue}
      elif [ 'x'${altqueue} != 'x' ]
      then
         newqueue=${altqueue}
      elif [ -s ${here}/${polyname}/skipper.txt ]
      then
         newqueue='unrestricted_parallel'
      fi
      oldline=${oi}
      newline=`echo ${oldline} | sed s@${queue}@${newqueue}@g`
      sed -i s@"${oldline}"@"${newline}"@g ${joborder}
      #------------------------------------------------------------------------------------#
   done
   #---------------------------------------------------------------------------------------#
   /bin/rm -f ${here}/run_sitter.lock
   exit
   #---------------------------------------------------------------------------------------#
fi
#------------------------------------------------------------------------------------------#



#----- Run R to make the status check. ----------------------------------------------------#
if [ ${plotstatus} -eq 1 ]
then
   echo 'Running the status check...'
   R CMD BATCH "${situtils}/status.r" "${situtils}/rout.rout"
fi
#------------------------------------------------------------------------------------------#



#----- Kill all jobs that can be killed. --------------------------------------------------#
echo 'Killing polygons that went extinct or reached steady state...'
if [ "x${deathrow}" != "x" ]
then
   for killme in ${deathrow}
   do
      bkill -J ${killme}
   done # killme in ${deathrow}
fi # [ "x${deathrow}" != "x" ]
#------------------------------------------------------------------------------------------#



#----- Give a three minute break. ---------------------------------------------------------#
echo 'Taking a nap before counting jobs...'
sleep 180

#----- Check how many jobs are left on our preferred queues. ------------------------------#
echo 'Counting the number of polygons on queues'
#----- Moorcroft_6100a. -------------------------------------------------------------------#
a61run=`bqueues moorcroft_6100a               | tail -1 | awk '{print $10}'`
a61pen=`bqueues moorcroft_6100a               | tail -1 | awk '{print  $9}'`
a61moi=`bjobs -q moorcroft_6100a 2> /dev/null | wc -l`
#----- Moorcroft_6100b. -------------------------------------------------------------------#
b61run=`bqueues moorcroft_6100b  | tail -1    | awk '{print $10}'`
b61pen=`bqueues moorcroft_6100b  | tail -1    | awk '{print  $9}'`
b61moi=`bjobs -q moorcroft_6100b 2> /dev/null | wc -l`
#----- Moorcroft2a queue.  We must count the "mortal" and the priority queues. ------------#
mmorun=`bqueues moorcroft2a            | tail -1    | awk '{print $10}'`
mmopen=`bqueues moorcroft2a            | tail -1    | awk '{print  $9}'`
mprrun=`bqueues moorcroft2a_priority   | tail -1    | awk '{print $10}'`
mprpen=`bqueues moorcroft2a_priority   | tail -1    | awk '{print  $9}'`
m2amoi=`bjobs -q moorcroft2a           2> /dev/null | wc -l`
let m2arun=${mmorun}+${mprrun}
let m2apen=${mmopen}+${mprpen}
#----- Moorcroft2b. -----------------------------------------------------------------------#
m2brun=`bqueues moorcroft2b      | tail -1    | awk '{print $10}'`
m2bpen=`bqueues moorcroft2b      | tail -1    | awk '{print  $9}'`
m2bmoi=`bjobs -q moorcroft2b     2> /dev/null | wc -l`
#----- Moorcroft2c. -----------------------------------------------------------------------#
m2crun=`bqueues moorcroft2c      | tail -1    | awk '{print $10}'`
m2cpen=`bqueues moorcroft2c      | tail -1    | awk '{print  $9}'`
m2cmoi=`bjobs -q moorcroft2c     2> /dev/null | wc -l`
#----- Wofsy queue.  We must count the "mortal" and the priority queues. ------------------#
wmorun=`bqueues wofsy            | tail -1    | awk '{print $10}'`
wmopen=`bqueues wofsy            | tail -1    | awk '{print  $9}'`
wprrun=`bqueues wofsy_priority   | tail -1    | awk '{print $10}'`
wprpen=`bqueues wofsy_priority   | tail -1    | awk '{print  $9}'`
wsymoi=`bjobs -q wofsy           2> /dev/null | wc -l`
let wsyrun=${wmorun}+${wprrun}
let wsypen=${wmopen}+${wprpen}
#----- Unrestricted serial.  Count only the jobs by the user. -----------------------------#
userun=`bjobs -q unrestricted_serial          2> /dev/null | grep RUN  | wc -l`
usepen=`bjobs -q unrestricted_serial          2> /dev/null | grep PEND | wc -l`


#----- Create a file with the pending jobs. -----------------------------------------------#
echo 'Counting the number of polygons in queue on long_serial...'
lsepen=`bjobs -w -J ${desc}-* -q long_serial 2> /dev/null | grep PEND  | wc -l`
lserun=`bjobs -w -J ${desc}-* -q long_serial 2> /dev/null | grep RUN   | wc -l`
/bin/rm -f ${pendfile}
bjobs -w -J ${desc}-* -q long_serial 2> /dev/null | grep PEND > ${pendfile}
#------------------------------------------------------------------------------------------#



#----- Find the available space on available queues. --------------------------------------#
echo 'Looking for empty cores on the queues...'
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#     Check the amount of room left for each queue.  Make sure that we don't exceed the    #
# maximum availability or the maximum use allowed.                                         #
#------------------------------------------------------------------------------------------#
#------ Moorcroft2a. ----------------------------------------------------------------------#
let m2atot=${m2arun}+${m2apen}
let m2aroom=${m2afull}-${m2atot}
let m2apot=${m2aumax}-${m2amoi}
if [ ${m2apot} -lt 0 ]
then
   m2aroom=0
elif [ ${m2apot} -lt ${m2aroom} ]
then
   m2aroom=${m2apot}
fi
echo " -----------------------------------------------------------"
echo "   Moorcroft2a"
echo " "
echo "   TOTAL     = ${m2atot} "
echo "   RUN       = ${m2arun} "
echo "   PEND      = ${m2apen} "
echo "   User max  = ${m2aumax}"
echo "   My runs   = ${m2amoi} "
echo "   Room      = ${m2aroom}"
echo "   Potential = ${m2apot} "
echo "   User max  = ${m2aumax}"
echo " -----------------------------------------------------------"
echo " "
#------ Moorcroft2b. ----------------------------------------------------------------------#
let m2btot=${m2brun}+${m2bpen}
let m2broom=${m2bfull}-${m2btot}
let m2bpot=${m2bumax}-${m2bmoi}
if [ ${m2bpot} -lt 0 ]
then
   m2broom=0
elif [ ${m2bpot} -lt ${m2broom} ]
then
   m2broom=${m2bpot}
fi
echo " -----------------------------------------------------------"
echo "   Moorcroft2b"
echo " "
echo "   TOTAL     = ${m2btot} "
echo "   RUN       = ${m2brun} "
echo "   PEND      = ${m2bpen} "
echo "   User max  = ${m2bumax}"
echo "   My runs   = ${m2bmoi} "
echo "   Room      = ${m2broom}"
echo "   Potential = ${m2bpot} "
echo "   User max  = ${m2bumax}"
echo " -----------------------------------------------------------"
echo " "
#------ Moorcroft2c. ----------------------------------------------------------------------#
let m2ctot=${m2crun}+${m2cpen}
let m2croom=${m2cfull}-${m2ctot}
let m2cpot=${m2cumax}-${m2cmoi}
if [ ${m2cpot} -lt 0 ]
then
   m2croom=0
elif [ ${m2cpot} -lt ${m2croom} ]
then
   m2croom=${m2cpot}
fi
echo " -----------------------------------------------------------"
echo "   Moorcroft2c"
echo " "
echo "   TOTAL     = ${m2ctot} "
echo "   RUN       = ${m2crun} "
echo "   PEND      = ${m2cpen} "
echo "   User max  = ${m2cumax}"
echo "   My runs   = ${m2cmoi} "
echo "   Room      = ${m2croom}"
echo "   Potential = ${m2cpot} "
echo "   User max  = ${m2cumax}"
echo " -----------------------------------------------------------"
echo " "
#------ Moorcroft_6100a. ------------------------------------------------------------------#
let a61tot=${a61run}+${a61pen}
let a61room=${a61full}-${a61tot}
let a61pot=${a61umax}-${a61moi}
if [ ${a61pot} -lt 0 ]
then
   a61room=0
elif [ ${a61pot} -lt ${a61room} ]
then
   a61room=${a61pot}
fi
echo " -----------------------------------------------------------"
echo "   Moorcroft_6100a"
echo " "
echo "   TOTAL     = ${a61tot} "
echo "   RUN       = ${a61run} "
echo "   PEND      = ${a61pen} "
echo "   User max  = ${a61umax}"
echo "   My runs   = ${a61moi} "
echo "   Room      = ${a61room}"
echo "   Potential = ${a61pot} "
echo "   User max  = ${a61umax}"
echo " -----------------------------------------------------------"
echo " "
#------ Moorcroft_6100b. ------------------------------------------------------------------#
let b61tot=${b61run}+${b61pen}
let b61room=${b61full}-${b61tot}
let b61pot=${b61umax}-${b61moi}
if [ ${b61pot} -lt 0 ]
then
   b61room=0
elif [ ${b61pot} -lt ${b61room} ]
then
   b61room=${b61pot}
fi
echo " -----------------------------------------------------------"
echo "   Moorcroft_6100b"
echo " "
echo "   TOTAL     = ${b61tot} "
echo "   RUN       = ${b61run} "
echo "   PEND      = ${b61pen} "
echo "   User max  = ${b61umax}"
echo "   My runs   = ${b61moi} "
echo "   Room      = ${b61room}"
echo "   Potential = ${b61pot} "
echo "   User max  = ${b61umax}"
echo " -----------------------------------------------------------"
echo " "
#------ Wofsy. ----------------------------------------------------------------------------#
let wsytot=${wsyrun}+${wsypen}
let wsyroom=${wsyfull}-${wsytot}
let wsypot=${wsyumax}-${wsymoi}
if [ ${wsypot} -lt 0 ]
then
   wsyroom=0
elif [ ${wsypot} -lt ${wsyroom} ]
then
   wsyroom=${wsypot}
fi
echo " -----------------------------------------------------------"
echo "   Wofsy"
echo " "
echo "   TOTAL     = ${wsytot} "
echo "   RUN       = ${wsyrun} "
echo "   PEND      = ${wsypen} "
echo "   User max  = ${wsyumax}"
echo "   My runs   = ${wsymoi} "
echo "   Room      = ${wsyroom}"
echo "   Potential = ${wsypot} "
echo "   User max  = ${wsyumax}"
echo " -----------------------------------------------------------"
echo " "
#------ Unrestricted_serial. --------------------------------------------------------------#
let usetot=${userun}+${usepen}
let useroom=${usefull}-${usetot}
let usepot=${useumax}-${usemoi}
if [ ${usepot} -lt 0 ]
then
   useroom=0
elif [ ${usepot} -lt ${useroom} ]
then
   useroom=${usepot}
fi
echo " -----------------------------------------------------------"
echo "   Unrestricted_serial"
echo " "
echo "   TOTAL     = ${usetot} "
echo "   RUN       = ${userun} "
echo "   PEND      = ${usepen} "
echo "   User max  = ${useumax}"
echo "   My runs   = ${usemoi} "
echo "   Room      = ${useroom}"
echo "   Potential = ${usepot} "
echo "   User max  = ${useumax}"
echo " -----------------------------------------------------------"
echo " "
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#      Count the total free room.                                                          #
#------------------------------------------------------------------------------------------#
let freeroom=${useroom}+${b61room}+${a61room}+${wsyroom}+${m2croom}+${m2broom}+${m2aroom}
let b61pa61=${b61room}+${a61room}
let b61pa61pwsy=${b61pa61}+${wsyroom}
let b61pa61pwsypm2c=${b61pa61pwsy}+${m2croom}
let b61pa61pwsypm2cpm2b=${b61pa61pwsypm2c}+${m2broom}
let b61pa61pwsypm2cpm2bpm2a=${b61pa61pwsypm2cpm2b}+${m2aroom}

#----- Find out how many cores can be moved. ----------------------------------------------#
if [ ${lsepen} -le ${freeroom} ]
then
   nfill=${lsepen}
else
   nfill=${freeroom}
fi

#----- This loop is where the job's queue switching happens. ------------------------------#
echo 'Switching queues for '${nfill}' jobs...'
ff=0
ll=0
while [ ${ff} -lt ${nfill} ]
do
   let ll=${ll}+1

   oi=`head -${ll} ${joborder} | tail -1`
   polyname=`echo ${oi} | awk '{print $1}'`
   thisjob="${desc}-${polyname}"

   #----- Check whether the job is pending or not.  If it is, then switch queues. ---------#
   ispending=`grep ${thisjob} ${pendfile} | wc -l`
   
   if [ ${ispending} -gt 0 ]
   then
      let ff=${ff}+1

      #----- Decide the new queue. --------------------------------------------------------#
      if [ ${ff} -le ${b61room} ]
      then
         newqueue='moorcroft_6100b'
      elif [ ${ff} -le ${b61pa61} ]
      then
         newqueue='moorcroft_6100a'
      elif [ ${ff} -le ${b61pa61pwsy} ]
      then
         newqueue='wofsy'
      elif [ ${ff} -le ${b61pa61pwsypm2c} ]
      then
         newqueue='moorcroft2c'
      elif [ ${ff} -le ${b61pa61pwsypm2cpm2b} ]
      then
         newqueue='moorcroft2b'
      elif [ ${ff} -le ${b61pa61pwsypm2cpm2bpm2a} ]
      then
         newqueue='moorcroft2a'
      else
         newqueue='unrestricted_serial'
      fi # [ ${ff} -le ${m2croom} ]

      bswitch -J ${thisjob} ${newqueue}

      #----- Update the queue in joborder.txt. --------------------------------------------#
      line=`grep ${polyname} ${joborder}`
      thislon=`echo ${line}  | awk '{print $2}'`
      thislat=`echo ${line}  | awk '{print $3}'`
      oldqueue=`echo ${line} | awk '{print $4}'`

      oldline=${polyname}' '${thislon}' '${thislat}' '${oldqueue}
      newline=${polyname}' '${thislon}' '${thislat}' '${newqueue}
      sed -i s@"${line}"@"${newline}"@g ${joborder}
   fi # [ ${ispending} -gt 0 ]
done # [ ${ff} -lt ${nfill} ]

#------------------------------------------------------------------------------------------#
#     If unrestricted_parallel has some room, we also create new jobs there.               #
#------------------------------------------------------------------------------------------#
uparun=`bjobs -q unrestricted_parallel 2> /dev/null | grep RUN  | wc -l`
upapen=`bjobs -q unrestricted_parallel 2> /dev/null | grep PEND | wc -l`
let uparoom=${upapmax}-${upapen}
if [ ${uparoom} -lt 0 ]
then
   uparoom=0
fi
let lsepen=${lsepen}-${nfill}
let lsenode=${lsepen}/8
if [ ${lsenode} -le ${uparoom} ]
then
   nodefill=${lsenode}
else
   nodefill=${uparoom}
fi

echo 'Creating '${nodefill}' new jobs to the unrestricted_parallel queue...'
nnff=0
while [ ${nnff} -lt ${nodefill} ]
do
   let ll=${ll}+1

   oi=`head -${ll} ${joborder} | tail -1`
   polyname=`echo ${oi} | awk '{print $1}'`
   thisjob="${desc}-${polyname}"

   #----- Check whether the job is pending or not.  If it is, then switch queues. ---------#
   ispending1=`grep ${thisjob} ${pendfile} | wc -l`
   
   if [ ${ispending1} -gt 0 ]
   then
      let nnff=${nnff}+1

      #----- Start a new callunpa. --------------------------------------------------------#
      /bin/rm -f ${callunpa}
      echo '#!/bin/sh' > ${callunpa}

      unpa=0
      while [ ${unpa} -lt 8 ]
      do            

         let unpa=${unpa}+1
         let remain=${unpa}%8 # % is mod operator
         let wtime=${remain}*15
         let wtime=${wtime}+2

         #---------------------------------------------------------------------------------#
         #     This part is done only when unpa is greater than 1 (because when unpa is 1  #
         # the pending is already defined).                                                #
         #---------------------------------------------------------------------------------#
         if [ ${unpa} -ne 1 ]
         then
            ispending2=0
            while [ ${ispending2} -eq 0 ]
            do
               let ll=${ll}+1

               oi=`head -${ll} ${joborder} | tail -1`
               polyname=`echo ${oi} | awk '{print $1}'`
               thisjob="${desc}-${polyname}"

               #---------------------------------------------------------------------------#
               #     Check whether the job is pending or not.  If it is, then switch       #
               # queues.                                                                   #
               #---------------------------------------------------------------------------#
               ispending2=`grep ${thisjob} ${pendfile} | wc -l`
            done
         fi

         #----- Put the skipping flag to the job. -----------------------------------------#
         echo 'Skip this node for now.' > ${here}/${polyname}/skipper.txt

         #----- Kill the job because it will go to the unrestricted_parallel. -------------#
         bkill -J ${thisjob}

         #----- Add the command to call this job from unrestricted_parallel driver. -------#
         mycomm="${here}/${polyname}/callserial.sh ${wtime} &"
         echo ${mycomm} >> ${callunpa}

         #----- Update the queue in joborder.txt. --------------------------------------------#
         line=`grep ${polyname} ${joborder}`
         thislon=`echo ${line} | awk '{print $2}'`
         thislat=`echo ${line} | awk '{print $3}'`
         oldqueue=`echo ${line} | awk '{print $4}'`
         newqueue='unrestricted_parallel'

         oldline=${polyname}' '${thislon}' '${thislat}' '${oldqueue}
         newline=${polyname}' '${thislon}' '${thislat}' '${newqueue}
         sed -i s@"${line}"@"${newline}"@g ${joborder}
      done
   
      #----- Finalise the script that controls submits the jobs, then submit the job... ------#
      echo 'wait' >> ${callunpa}
      chmod u+x ${callunpa}
      mv ${callunpa} ${here}/${polyname}


      #----- Create the shell script that will call bsub. ------------------------------------#
      /bin/rm -f ${unparun}
      echo '#!/bin/sh' > ${unparun}
      bsub="bsub -q unrestricted_parallel -R 'hname!=hero3114'"
      bsub="${bsub} -J ${thisjob}"
      bsub="${bsub} -o ${here}/${polyname}/serial_lsf.out -n 8"
      bsub="${bsub} < ${here}/${polyname}/"`basename ${callunpa}`
      echo ${bsub} >> ${unparun}
      chmod u+x ${unparun}
      mv ${unparun} ${here}/${polyname}

      #----- Submit the jobs to the queue. ---------------------------------------------------#
      ${here}/${polyname}/unparun.sh
   fi # [ ${ispending1} -gt 0 ]
done
let nptounpa=${nodefill}*8


#----- Get rid of the temporary file. -----------------------------------------------------#
/bin/rm -f ${pendfile}

#----- Sleep 5 more minutes. --------------------------------------------------------------#
echo 'Taking another nap before checking the queues again...'
sleep 300

#----- Quick check the status of the queues. ----------------------------------------------#
/bin/rm -f ${recefile}
touch ${recefile}
echo '------- Polygon recent activity status. ----------------------------' >> ${recefile}
echo 'Number of polygons that were re-submitted:    '${nresubmit}  >> ${recefile}
echo 'Number of polygons that switched queues:      '${nfill}      >> ${recefile}
echo 'Number of polygons to unrestricted_parallel:  '${nptounpa}   >> ${recefile}
echo 'Number of polygons that went extinct:         '${newextinct} >> ${recefile}
echo 'Number of polygons that reached steady state: '${newststate} >> ${recefile}
echo '--------------------------------------------------------------------' >> ${recefile}

#----- Check the queue status. ------------------------------------------------------------#
echo 'Counting the jobs...'
a61run=`bjobs -J ${desc}-* -w -q moorcroft_6100a       2> /dev/null | grep RUN   | wc -l`
b61run=`bjobs -J ${desc}-* -w -q moorcroft_6100b       2> /dev/null | grep RUN   | wc -l`
m2arun=`bjobs -J ${desc}-* -w -q moorcroft2a           2> /dev/null | grep RUN   | wc -l`
m2brun=`bjobs -J ${desc}-* -w -q moorcroft2b           2> /dev/null | grep RUN   | wc -l`
m2crun=`bjobs -J ${desc}-* -w -q moorcroft2c           2> /dev/null | grep RUN   | wc -l`
wsyrun=`bjobs -J ${desc}-* -w -q wofsy                 2> /dev/null | grep RUN   | wc -l`
lserun=`bjobs -J ${desc}-* -w -q long_serial           2> /dev/null | grep RUN   | wc -l`
userun=`bjobs -J ${desc}-* -w -q unrestricted_serial   2> /dev/null | grep RUN   | wc -l`
uparun=`bjobs -J ${desc}-* -w -q unrestricted_parallel 2> /dev/null | grep RUN   | wc -l`
a61pen=`bjobs -J ${desc}-* -w -q moorcroft_6100a       2> /dev/null | grep PEND  | wc -l`
b61pen=`bjobs -J ${desc}-* -w -q moorcroft_6100b       2> /dev/null | grep PEND  | wc -l`
m2apen=`bjobs -J ${desc}-* -w -q moorcroft2a           2> /dev/null | grep PEND  | wc -l`
m2bpen=`bjobs -J ${desc}-* -w -q moorcroft2b           2> /dev/null | grep PEND  | wc -l`
m2cpen=`bjobs -J ${desc}-* -w -q moorcroft2c           2> /dev/null | grep PEND  | wc -l`
wsypen=`bjobs -J ${desc}-* -w -q wofsy                 2> /dev/null | grep PEND  | wc -l`
lsepen=`bjobs -J ${desc}-* -w -q long_serial           2> /dev/null | grep PEND  | wc -l`
usepen=`bjobs -J ${desc}-* -w -q unrestricted_serial   2> /dev/null | grep PEND  | wc -l`
upapen=`bjobs -J ${desc}-* -w -q unrestricted_parallel 2> /dev/null | grep PEND  | wc -l`
lsertot=`bjobs -q long_serial                          2> /dev/null | grep RUN  | wc -l`
m2artot=`bjobs -q moorcroft2a                          2> /dev/null | grep RUN  | wc -l`
m2brtot=`bjobs -q moorcroft2b                          2> /dev/null | grep RUN  | wc -l`
m2crtot=`bjobs -q moorcroft2c                          2> /dev/null | grep RUN  | wc -l`
wsyrtot=`bjobs -q wofsy                                2> /dev/null | grep RUN  | wc -l`
a61rtot=`bjobs -q moorcroft_6100a                      2> /dev/null | grep RUN  | wc -l`
b61rtot=`bjobs -q moorcroft_6100b                      2> /dev/null | grep RUN  | wc -l`
usertot=`bjobs -q unrestricted_serial                  2> /dev/null | grep RUN  | wc -l`
upartot=`bjobs -q unrestricted_parallel                2> /dev/null | grep RUN  | wc -l`
lseptot=`bjobs -q long_serial                          2> /dev/null | grep PEND | wc -l`
m2aptot=`bjobs -q moorcroft2a                          2> /dev/null | grep PEND | wc -l`
m2bptot=`bjobs -q moorcroft2b                          2> /dev/null | grep PEND | wc -l`
m2cptot=`bjobs -q moorcroft2c                          2> /dev/null | grep PEND | wc -l`
wsyptot=`bjobs -q wofsy                                2> /dev/null | grep PEND | wc -l`
a61ptot=`bjobs -q moorcroft_6100a                      2> /dev/null | grep PEND | wc -l`
b61ptot=`bjobs -q moorcroft_6100b                      2> /dev/null | grep PEND | wc -l`
useptot=`bjobs -q unrestricted_serial                  2> /dev/null | grep PEND | wc -l`
upaptot=`bjobs -q unrestricted_parallel                2> /dev/null | grep PEND | wc -l`
let lsertot=${lsertot}+${lseptot}
let m2artot=${m2artot}+${m2aptot}
let m2brtot=${m2brtot}+${m2bptot}
let m2crtot=${m2crtot}+${m2cptot}
let wsyrtot=${wsyrtot}+${wsyptot}
let a61rtot=${a61rtot}+${a61ptot}
let b61rtot=${b61rtot}+${b61ptot}
let usertot=${usertot}+${useptot}
let upartot=${upartot}+${upaptot}

/bin/rm -f ${queuefile}
touch ${queuefile}
echo "------- Queue status. --------------------------------------------" >> ${queuefile}
echo "Moorcroft_6100a        RUN=${a61run}  PEN=${a61pen}  TOT=${a61tot}" >> ${queuefile}
echo "Moorcroft_6100b        RUN=${b61run}  PEN=${b61pen}  TOT=${b61tot}" >> ${queuefile}
echo "Moorcroft2a            RUN=${m2arun}  PEN=${m2apen}  TOT=${m2atot}" >> ${queuefile}
echo "Moorcroft2b            RUN=${m2brun}  PEN=${m2bpen}  TOT=${m2btot}" >> ${queuefile}
echo "Moorcroft2c            RUN=${m2crun}  PEN=${m2cpen}  TOT=${m2ctot}" >> ${queuefile}
echo "Wofsy                  RUN=${wsyrun}  PEN=${wsypen}  TOT=${wsytot}" >> ${queuefile}
echo "Long_serial            RUN=${lserun}  PEN=${lsepen}  TOT=${lsetot}" >> ${queuefile}
echo "Unrestricted_serial    RUN=${userun}  PEN=${usepen}  TOT=${usetot}" >> ${queuefile}
echo "Unrestricted_parallel  RUN=${uparun}  PEN=${upapen}  TOT=${upatot}" >> ${queuefile}
echo "------------------------------------------------------------------" >> ${queuefile}

#----- Current simulation status. ---------------------------------------------------------#
n_initial=`grep INITIAL ${outcheck} | wc -l`
n_history=`grep HISTORY ${outcheck} | wc -l`
n_metmiss=`grep METMISS ${outcheck} | wc -l`
n_crashed=`grep CRASHED ${outcheck} | wc -l`
n_stopped=`grep STOPPED ${outcheck} | wc -l`
n_extinct=`grep EXTINCT ${outcheck} | wc -l`
n_ststate=`grep STSTATE ${outcheck} | wc -l`
n_the_end=`grep THE_END ${outcheck} | wc -l`
/bin/rm -f ${statfile}
touch ${statfile}
echo "------- Simulation status. -----------------------------------------" >> ${statfile}
echo " Number of polygons that have never started           : ${n_initial}" >> ${statfile}
echo " Number of polygons that have partially run           : ${n_history}" >> ${statfile}
echo " Number of polygons that haven't found met drivers    : ${n_metmiss}" >> ${statfile}
echo " Number of polygons that have crashed                 : ${n_crashed}" >> ${statfile}
echo " Number of polygons that have mysteriously stopped    : ${n_stopped}" >> ${statfile}
echo " Number of polygons that became desert                : ${n_extinct}" >> ${statfile}
echo " Number of polygons that have reached steady state    : ${n_ststate}" >> ${statfile}
echo " Number of polygons that have reached the end         : ${n_the_end}" >> ${statfile}
echo "--------------------------------------------------------------------" >> ${statfile}

#----- Convert the PNG files to an e-mail friendly format. --------------------------------#
if [ ${plotstatus} -eq 1 ]
then
   echo 'Converting the PNG files to e-mail format...'
   uuencode ${situtils}/polystatus.png polystat.png > ${situtils}/polystat.txt
   uuencode ${situtils}/polyagb.png    polyagb.png  > ${situtils}/polyagb.txt
   uuencode ${situtils}/polybsa.png    polybsa.png  > ${situtils}/polybsa.txt
   uuencode ${situtils}/polylai.png    polylai.png  > ${situtils}/polylai.txt
   uuencode ${situtils}/polyvel.png    polyvel.png  > ${situtils}/polyvel.txt
fi

#------ Build the e-mail body. ------------------------------------------------------------#
echo 'Building the e-mail...'
/bin/rm -f ${emailbody}
touch ${emailbody}
cat ${headfile}  >> ${emailbody}
echo ' '         >> ${emailbody}
cat ${statfile}  >> ${emailbody}
echo ' '         >> ${emailbody}
cat ${recefile}  >> ${emailbody}
echo ' '         >> ${emailbody}
cat ${queuefile} >> ${emailbody}
echo ' '         >> ${emailbody}
cat ${tailfile}  >> ${emailbody}
echo ' '         >> ${emailbody}

if [ ${plotstatus} -eq 1 ]
then
   cat ${situtils}/polystat.txt >> ${emailbody}
   echo ' '         >> ${emailbody}
   cat ${situtils}/polyagb.txt  >> ${emailbody}
   echo ' '         >> ${emailbody}
   cat ${situtils}/polybsa.txt  >> ${emailbody}
   echo ' '         >> ${emailbody}
   cat ${situtils}/polylai.txt  >> ${emailbody}
   echo ' '         >> ${emailbody}
   cat ${situtils}/polyvel.txt  >> ${emailbody}
   echo ' '         >> ${emailbody}
fi

when=`date +'%B %d, %Y - %R %Z'`
today=`date +'%Y-%m-%d'`
subject="${desc} run status as of ${when}"

if [ ${email1day} -eq 0 ]
then
  echo ${today} > ${emailday}
  echo 'Sending e-mail...'
  mail -s "${subject}" ${recipient} < ${emailbody}
elif [ -s ${emailday} ]
then
   lastemail=`cat ${emailday}`
   if [ ${today} != ${lastemail} ]
   then
      echo ${today} > ${emailday}
      echo 'Sending e-mail...'
      mail -s "${subject}" ${recipient} < ${emailbody}
   else
      echo 'Skipping e-mail...'
   fi
else
  echo ${today} > ${emailday}
  echo 'Sending e-mail...'
  mail -s "${subject}" ${recipient} < ${emailbody}
fi

#----- Clean-up stuff. --------------------------------------------------------------------#
echo 'Deleting some temporary files...'
/bin/rm -f ${emailbody} ${queuefile} ${recefile}
/bin/rm -f ${situtils}/polystat.txt
/bin/rm -f ${situtils}/polyagb.txt
/bin/rm -f ${situtils}/polybsa.txt
/bin/rm -f ${situtils}/polylai.txt
/bin/rm -f ${here}/run_sitter.lock
echo '==== run_sitter.sh execution ends. ===='
