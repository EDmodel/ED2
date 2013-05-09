#==========================================================================================#
#==========================================================================================#
#     This function reads the ED2 monthly mean files that contain mean diurnal cycle.      #
#   Inputs:                                                                                #
#   - datum   -- The monthly structure that will contain the data.  It must be initialised #
#                by create.monthly, otherwise it won't work.                               #
#   - ntimes  -- Total number of times (including previously loaded times).                #
#   - tresume -- The first time to read (in case data have been partially loaded.          #
#------------------------------------------------------------------------------------------#
read.q.files <<- function(datum,ntimes,tresume=1,sasmonth=5){


   #----- Copy some dimensions to scalars. ------------------------------------------------#
   nzg        = datum$nzg
   nzs        = datum$nzs
   ndcycle    = datum$ndcycle
   isoilflg   = datum$isoilflg
   slz        = datum$slz
   slxsand    = datum$slxsand
   slxclay    = datum$slxclay
   ntext      = datum$ntext
   soil.prop  = datum$soil.prop
   dslz       = datum$dslz
   soil.depth = datum$soil.depth
   soil.dry   = datum$soil.dry
   soil.poro  = datum$soil.poro
   ka         = datum$ka
   kz         = datum$kz
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Copy the variables to scratch lists, we will copy them back once we are done.     #
   #---------------------------------------------------------------------------------------#
   emean  = datum$emean
   emsqu  = datum$emsqu
   szpft  = datum$szpft
   lu     = datum$lu
   qmean  = datum$qmean
   qmsqu  = datum$qmsqu
   patch  = datum$patch
   cohort = datum$cohort
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Loop over all times that haven't been read yet.                                   #
   #---------------------------------------------------------------------------------------#
   for (m in tresume:ntimes){

      #----- Print a banner to entertain the bored user staring at the screen. ------------#
      if (m == tresume | datum$month[m] == 1){
         cat("    - Reading data from year ",datum$year[m],"...","\n")
      }#end if
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Number of days in a month.                                                     #
      #------------------------------------------------------------------------------------#
      mondays   = daymax(datum$when[m])
      thismonth = datum$month[m]
      lastmonth = 1 + (thismonth - 2) %% 12
      thisyear  = datum$year [m]
      #------------------------------------------------------------------------------------#




      #----- Read data and close connection immediately after. ----------------------------#
      h5file       = datum$input[m]
      h5file.bz2   = paste(datum$input[m],"bz2",sep=".")
      h5file.gz    = paste(datum$input[m],"gz" ,sep=".")
      if (file.exists(h5file)){
         mymont    = hdf5load(file=h5file,load=FALSE,verbosity=0,tidy=TRUE)

      }else if(file.exists(h5file.bz2)){
         temp.file = file.path(tempdir(),basename(h5file))
         dummy     = bunzip2(filename=h5file.bz2,destname=temp.file,remove=FALSE)
         mymont    = hdf5load(file=temp.file,load=FALSE,verbosity=0,tidy=TRUE)
         dummy     = file.remove(temp.file)

      }else if(file.exists(h5file.gz)){
         temp.file = file.path(tempdir(),basename(h5file))
         dummy     = gunzip(filename=h5file.gz,destname=temp.file,remove=FALSE)
         mymont    = hdf5load(file=temp.file,load=FALSE,verbosity=0,tidy=TRUE)
         dummy     = file.remove(temp.file)
      }else{
         cat (" - File      : ",basename(h5file)    ,"\n")
         cat (" - File (bz2): ",basename(h5file.bz2),"\n")
         stop(" Neither the expanded nor the compressed files were found!")

      }#end if
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #      Create some additional radiation variables.                                   #
      #------------------------------------------------------------------------------------#
      #----- Direct radiation. ------------------------------------------------------------#
      mymont$MMEAN.ATM.RSHORT.BEAM.PY = ( mymont$MMEAN.ATM.RSHORT.PY 
                                        - mymont$MMEAN.ATM.RSHORT.DIFF.PY )
      mymont$MMEAN.ATM.PAR.BEAM.PY    = ( mymont$MMEAN.ATM.PAR.PY 
                                        - mymont$MMEAN.ATM.PAR.DIFF.PY    )
      mymont$MMEAN.ATM.RSHORT.BEAM.SI = ( mymont$MMEAN.ATM.RSHORT.SI 
                                        - mymont$MMEAN.ATM.RSHORT.DIFF.SI )
      mymont$MMEAN.ATM.PAR.BEAM.SI    = ( mymont$MMEAN.ATM.PAR.SI 
                                        - mymont$MMEAN.ATM.PAR.DIFF.SI    )
      mymont$QMEAN.ATM.RSHORT.BEAM.PY = ( mymont$QMEAN.ATM.RSHORT.PY 
                                        - mymont$QMEAN.ATM.RSHORT.DIFF.PY )
      mymont$QMEAN.ATM.PAR.BEAM.PY    = ( mymont$QMEAN.ATM.PAR.PY 
                                        - mymont$QMEAN.ATM.PAR.DIFF.PY    )
      mymont$QMEAN.ATM.RSHORT.BEAM.SI = ( mymont$QMEAN.ATM.RSHORT.SI 
                                        - mymont$QMEAN.ATM.RSHORT.DIFF.SI )
      mymont$QMEAN.ATM.PAR.BEAM.SI    = ( mymont$QMEAN.ATM.PAR.SI 
                                        - mymont$QMEAN.ATM.PAR.DIFF.SI    )
      #----- Near infrared. ---------------------------------------------------------------#
      mymont$MMEAN.ATM.NIR.PY         = ( mymont$MMEAN.ATM.RSHORT.PY
                                        - mymont$MMEAN.ATM.PAR.PY          )
      mymont$MMEAN.ATM.NIR.DIFF.PY    = ( mymont$MMEAN.ATM.RSHORT.DIFF.PY
                                        - mymont$MMEAN.ATM.PAR.DIFF.PY     )
      mymont$MMEAN.ATM.NIR.BEAM.PY    = ( mymont$MMEAN.ATM.RSHORT.BEAM.PY
                                        - mymont$MMEAN.ATM.PAR.BEAM.PY     )
      mymont$QMEAN.ATM.NIR.PY         = ( mymont$QMEAN.ATM.RSHORT.PY
                                        - mymont$QMEAN.ATM.PAR.PY          )
      mymont$QMEAN.ATM.NIR.DIFF.PY    = ( mymont$QMEAN.ATM.RSHORT.DIFF.PY
                                        - mymont$QMEAN.ATM.PAR.DIFF.PY     )
      mymont$QMEAN.ATM.NIR.BEAM.PY    = ( mymont$QMEAN.ATM.RSHORT.BEAM.PY
                                        - mymont$QMEAN.ATM.PAR.BEAM.PY     )
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Find the mean latent heat of vaporisation.  Because we assume it to be a       #
      # linear function of temperature, the mean can be found a posteriori.  The mean      #
      # fluxes won't be exact though, because the covariance part is missing.              #
      #------------------------------------------------------------------------------------#
      mmean.can.alvli.py  = alvli(mymont$MMEAN.CAN.TEMP.PY)
      qmean.can.alvli.py  = alvli(mymont$QMEAN.CAN.TEMP.PY)
      mmean.can.alvli.pa  = alvli(mymont$MMEAN.CAN.TEMP.PA)
      qmean.can.alvli.pa  = alvli(mymont$QMEAN.CAN.TEMP.PA)
      mmean.can.alvli.py2 = mmean.can.alvli.py * mmean.can.alvli.py
      qmean.can.alvli.py2 = qmean.can.alvli.py * qmean.can.alvli.py
      mmean.can.alvli.pa2 = mmean.can.alvli.pa * mmean.can.alvli.pa
      qmean.can.alvli.pa2 = qmean.can.alvli.pa * qmean.can.alvli.pa
      #------------------------------------------------------------------------------------#



      #----- Load the total number of patches and cohorts. --------------------------------#
      emean$npat.global[m] = mymont$NPATCHES.GLOBAL
      emean$ncoh.global[m] = mymont$NCOHORTS.GLOBAL
      #------------------------------------------------------------------------------------#


      #----- Load the simple variables. ---------------------------------------------------#
      emean$fast.soil.c     [m] =   mymont$MMEAN.FAST.SOIL.C.PY
      emean$slow.soil.c     [m] =   mymont$MMEAN.SLOW.SOIL.C.PY
      emean$struct.soil.c   [m] =   mymont$MMEAN.STRUCT.SOIL.C.PY
      emean$het.resp        [m] =   mymont$MMEAN.RH.PY
      emean$cwd.resp        [m] =   mymont$MMEAN.CWD.RH.PY
      emean$gpp             [m] =   mymont$MMEAN.GPP.PY
      emean$npp             [m] =   mymont$MMEAN.NPP.PY
      emean$plant.resp      [m] =   mymont$MMEAN.PLRESP.PY
      emean$leaf.resp       [m] =   mymont$MMEAN.LEAF.RESP.PY
      emean$root.resp       [m] =   mymont$MMEAN.ROOT.RESP.PY
      emean$growth.resp     [m] =   mymont$MMEAN.GROWTH.RESP.PY
      emean$reco            [m] =   mymont$MMEAN.PLRESP.PY    + mymont$MMEAN.RH.PY
      emean$nep             [m] =   mymont$MMEAN.NEP.PY
      emean$nee             [m] =   mymont$MMEAN.CARBON.ST.PY - mymont$MMEAN.CARBON.AC.PY
      emean$ustar           [m] =   mymont$MMEAN.USTAR.PY
      emean$cflxca          [m] = - mymont$MMEAN.CARBON.AC.PY
      emean$cflxst          [m] =   mymont$MMEAN.CARBON.ST.PY
      emean$ustar           [m] =   mymont$MMEAN.USTAR.PY
      emean$atm.vels        [m] =   mymont$MMEAN.ATM.VELS.PY
      emean$atm.prss        [m] =   mymont$MMEAN.ATM.PRSS.PY   * 0.01
      emean$atm.temp        [m] =   mymont$MMEAN.ATM.TEMP.PY   - t00
      emean$atm.shv         [m] =   mymont$MMEAN.ATM.SHV.PY    * kg2g
      emean$atm.co2         [m] =   mymont$MMEAN.ATM.CO2.PY
      emean$atm.vpd         [m] =   mymont$MMEAN.ATM.VPDEF.PY  * 0.01
      emean$can.prss        [m] =   mymont$MMEAN.CAN.PRSS.PY   * 0.01
      emean$can.temp        [m] =   mymont$MMEAN.CAN.TEMP.PY   - t00
      emean$can.shv         [m] =   mymont$MMEAN.CAN.SHV.PY    * kg2g
      emean$can.co2         [m] =   mymont$MMEAN.CAN.CO2.PY
      emean$can.vpd         [m] =   mymont$MMEAN.CAN.VPDEF.PY  * 0.01
      emean$gnd.temp        [m] =   mymont$MMEAN.GND.TEMP.PY   - t00
      emean$gnd.shv         [m] =   mymont$MMEAN.GND.SHV.PY    * kg2g
      emean$leaf.temp       [m] =   mymont$MMEAN.LEAF.TEMP.PY  - t00
      emean$leaf.vpd        [m] =   mymont$MMEAN.LEAF.VPDEF.PY * 0.01
      emean$wood.temp       [m] =   mymont$MMEAN.WOOD.TEMP.PY  - t00
      emean$hflxca          [m] = - mymont$MMEAN.SENSIBLE.AC.PY
      emean$qwflxca         [m] = - mymont$MMEAN.VAPOR.AC.PY   * mmean.can.alvli.py
      emean$hflxgc          [m] =   mymont$MMEAN.SENSIBLE.GC.PY
      emean$hflxlc          [m] =   mymont$MMEAN.SENSIBLE.LC.PY
      emean$hflxwc          [m] =   mymont$MMEAN.SENSIBLE.WC.PY
      emean$wflxca          [m] = - mymont$MMEAN.VAPOR.AC.PY   * day.sec
      emean$wflxgc          [m] =   mymont$MMEAN.VAPOR.GC.PY   * day.sec
      emean$wflxlc          [m] =   mymont$MMEAN.VAPOR.LC.PY   * day.sec
      emean$wflxwc          [m] =   mymont$MMEAN.VAPOR.WC.PY   * day.sec
      emean$runoff          [m] = ( mymont$MMEAN.RUNOFF.PY
                                  + mymont$MMEAN.DRAINAGE.PY       ) * mondays * day.sec
      emean$intercepted     [m] = ( mymont$MMEAN.INTERCEPTED.AL.PY
                                  + mymont$MMEAN.INTERCEPTED.AW.PY ) * mondays * day.sec
      emean$wshed           [m] = ( mymont$MMEAN.WSHED.LG.PY
                                  + mymont$MMEAN.WSHED.WG.PY       ) * mondays * day.sec
      emean$evap            [m] = ( mymont$MMEAN.VAPOR.GC.PY
                                  + mymont$MMEAN.VAPOR.LC.PY
                                  + mymont$MMEAN.VAPOR.WC.PY ) * day.sec
      emean$transp          [m] =   mymont$MMEAN.TRANSP.PY     * day.sec

      #----- Save ET and the past ET [mm/day]. --------------------------------------------#
      emean$et              [m] = emean$evap[m] + emean$transp[m]
      emean$last.1yr.et     [m] = mean(emean$et  [max(m-11,1):m])
      emean$last.2yr.et     [m] = mean(emean$et  [max(m-23,1):m])
      emean$last.3yr.et     [m] = mean(emean$et  [max(m-35,1):m])
      #------------------------------------------------------------------------------------#


      #----- Save rainfall and the past rainfall. -----------------------------------------#
      emean$rain            [m] =   mymont$MMEAN.PCPG.PY * mondays * day.sec
      emean$last.1yr.rain   [m] = sum(emean$rain[max(m-11,1):m])
      emean$last.2yr.rain   [m] = sum(emean$rain[max(m-23,1):m]) * 0.5
      emean$last.3yr.rain   [m] = sum(emean$rain[max(m-35,1):m]) * onethird
      #------------------------------------------------------------------------------------#

      if ( m == 1){
         emean$ra.rain    [m] = emean$rain[m]
         emean$nmon.lt.090[m] = as.numeric(emean$ra.rain[m] <  90  )
         emean$nmon.lt.100[m] = as.numeric(emean$ra.rain[m] < 100  )
         emean$nmon.lt.110[m] = as.numeric(emean$ra.rain[m] < 110  )
         emean$nmon.lt.120[m] = as.numeric(emean$ra.rain[m] < 120  )
      }else{
         emean$ra.rain    [m] = ( emean$rain[m] + 11. * emean$ra.rain[m-1] ) / 12
         emean$nmon.lt.090[m] = ( as.numeric(emean$ra.rain[m] < 90  )
                                * ( emean$nmon.lt.090[m-1] + 1 ) )
         emean$nmon.lt.100[m] = ( as.numeric(emean$ra.rain[m] < 100 )
                                * ( emean$nmon.lt.100[m-1] + 1 ) )
         emean$nmon.lt.110[m] = ( as.numeric(emean$ra.rain[m] < 110 )
                                * ( emean$nmon.lt.110[m-1] + 1 ) )
         emean$nmon.lt.120[m] = ( as.numeric(emean$ra.rain[m] < 120 )
                                * ( emean$nmon.lt.120[m-1] + 1 ) )
      }#end if


      emean$fs.open         [m] =   mymont$MMEAN.FS.OPEN.PY
      emean$rshort          [m] =   mymont$MMEAN.ATM.RSHORT.PY
      emean$rshort.beam     [m] = ( mymont$MMEAN.ATM.RSHORT.PY
                                  - mymont$MMEAN.ATM.RSHORT.DIFF.PY )
      emean$rshort.diff     [m] =   mymont$MMEAN.ATM.RSHORT.DIFF.PY
      emean$rshortup        [m] =   mymont$MMEAN.RSHORTUP.PY
      emean$rlong           [m] =   mymont$MMEAN.ATM.RLONG.PY
      emean$rshort.gnd      [m] =   mymont$MMEAN.RSHORT.GND.PY
      emean$rlong.gnd       [m] =   mymont$MMEAN.RLONG.GND.PY
      emean$rlongup         [m] =   mymont$MMEAN.RLONGUP.PY
      emean$par.tot         [m] =   mymont$MMEAN.ATM.PAR.PY        * Watts.2.Ein * 1.e6
      emean$par.beam        [m] = ( mymont$MMEAN.ATM.PAR.PY
                                  - mymont$MMEAN.ATM.PAR.DIFF.PY ) * Watts.2.Ein * 1.e6
      emean$par.diff        [m] =   mymont$MMEAN.ATM.PAR.DIFF.PY   * Watts.2.Ein * 1.e6
      emean$par.gnd         [m] =   mymont$MMEAN.PAR.GND.PY        * Watts.2.Ein * 1.e6
      emean$parup           [m] =   mymont$MMEAN.PARUP.PY          * Watts.2.Ein * 1.e6
      emean$rnet            [m] =   mymont$MMEAN.RNET.PY
      emean$albedo          [m] =   mymont$MMEAN.ALBEDO.PY
      if (all(c("MMEAN.ALBEDO.PAR.PY","MMEAN.ALBEDO.NIR.PY") %in% names(mymont))){
         emean$albedo.par   [m] =   mymont$MMEAN.ALBEDO.PAR.PY
         emean$albedo.nir   [m] =   mymont$MMEAN.ALBEDO.NIR.PY
      }else{
         emean$albedo.par   [m] = ifelse( mymont$MMEAN.ATM.PAR.PY > 0.5
                                        , mymont$MMEAN.PARUP.PY / mymont$MMEAN.ATM.PAR.PY
                                        , mymont$MMEAN.ALBEDO.PY
                                        )#end ifelse
         emean$albedo.nir   [m] = ifelse( mymont$MMEAN.ATM.NIR.PY > 0.5
                                        , mymont$MMEAN.NIRUP.PY / mymont$MMEAN.ATM.NIR.PY
                                        , mymont$MMEAN.ALBEDO.PY
                                        )#end ifelse
      }#end if
      emean$rlong.albedo    [m] =   mymont$MMEAN.RLONG.ALBEDO.PY
      emean$leaf.gbw        [m] =   mymont$MMEAN.LEAF.GBW.PY * day.sec
      emean$leaf.gsw        [m] =   mymont$MMEAN.LEAF.GSW.PY * day.sec
      emean$wood.gbw        [m] =   mymont$MMEAN.WOOD.GBW.PY * day.sec
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     The following variables must be aggregated because the polygon-level is split  #
      # by PFT and DBH class.                                                              #
      #------------------------------------------------------------------------------------#
      emean$mco             [m] = apply( X      = ( mymont$MMEAN.LEAF.MAINTENANCE.PY
                                                  + mymont$MMEAN.ROOT.MAINTENANCE.PY
                                                  ) * yr.day
                                       , MARGIN = 1
                                       , FUN    = sum
                                       )#end if
      emean$ldrop           [m] = apply( X      = mymont$MMEAN.LEAF.DROP.PY * yr.day
                                       , MARGIN = 1
                                       , FUN    = sum
                                       )#end if
      #------------------------------------------------------------------------------------#




      #------ Read in soil properties. ----------------------------------------------------#
      emean$soil.temp  [m,] =   mymont$MMEAN.SOIL.TEMP.PY - t00
      emean$soil.water [m,] =   mymont$MMEAN.SOIL.WATER.PY
      emean$soil.mstpot[m,] = - mymont$MMEAN.SOIL.MSTPOT.PY
      #------------------------------------------------------------------------------------#



      #----- Find averaged soil properties. -----------------------------------------------#
      swater.now     = rev(cumsum(rev(mymont$MMEAN.SOIL.WATER.PY * wdns * dslz)))
      smoist.avg     = swater.now / (wdns * soil.depth)
      emean$paw  [m] = 100. * ( ( swater.now[ka] - soil.dry [ka] )
                              / ( soil.poro [ka] - soil.dry [ka] ) )
      emean$smpot[m] = ( - smoist2mpot(smoist=smoist.avg[ka],mysoil=soil.prop)
                       * 0.001 * grav )
      #------------------------------------------------------------------------------------#



      #----- Read workload, and retrieve only the current month. --------------------------#
      emean$workload  [m] = mymont$WORKLOAD[thismonth]
      emean$specwork  [m] = mymont$WORKLOAD[thismonth] / sum(mymont$SIPA.N,na.rm=TRUE)
      #------------------------------------------------------------------------------------#





      #------------------------------------------------------------------------------------#
      #     Retrieve the sum of squares (that will be used to find standard deviation.     #
      #------------------------------------------------------------------------------------#
      emsqu$gpp       [m] =   mymont$MMSQU.GPP.PY
      emsqu$plant.resp[m] =   mymont$MMSQU.PLRESP.PY
      emsqu$het.resp  [m] =   mymont$MMSQU.RH.PY
      emsqu$cwd.resp  [m] =   mymont$MMSQU.CWD.RH.PY
      emsqu$cflxca    [m] =   mymont$MMSQU.CARBON.AC.PY
      emsqu$cflxst    [m] =   mymont$MMSQU.CARBON.ST.PY
      emsqu$hflxca    [m] =   mymont$MMSQU.SENSIBLE.AC.PY
      emsqu$hflxlc    [m] =   mymont$MMSQU.SENSIBLE.LC.PY
      emsqu$hflxwc    [m] =   mymont$MMSQU.SENSIBLE.WC.PY
      emsqu$hflxgc    [m] =   mymont$MMSQU.SENSIBLE.GC.PY
      emsqu$wflxca    [m] =   mymont$MMSQU.VAPOR.AC.PY  * day.sec2
      emsqu$qwflxca   [m] =   mymont$MMSQU.VAPOR.AC.PY  * mmean.can.alvli.py2
      emsqu$wflxlc    [m] =   mymont$MMSQU.VAPOR.LC.PY  * day.sec2
      emsqu$wflxwc    [m] =   mymont$MMSQU.VAPOR.WC.PY  * day.sec2
      emsqu$wflxgc    [m] =   mymont$MMSQU.VAPOR.GC.PY  * day.sec2
      emsqu$evap      [m] = ( mymont$MMSQU.VAPOR.LC.PY
                            + mymont$MMSQU.VAPOR.WC.PY
                            + mymont$MMSQU.VAPOR.GC.PY
                            + 2. * ( mymont$MMEAN.VAPOR.LC.PY
                                   * mymont$MMEAN.VAPOR.WC.PY
                                   + mymont$MMEAN.VAPOR.LC.PY
                                   * mymont$MMEAN.VAPOR.GC.PY
                                   + mymont$MMEAN.VAPOR.WC.PY
                                   * mymont$MMEAN.VAPOR.GC.PY ) ) * day.sec2
      emsqu$transp    [m] =   mymont$MMSQU.TRANSP.PY    * day.sec2
      emsqu$ustar     [m] =   mymont$MMSQU.USTAR.PY
      emsqu$albedo    [m] =   mymont$MMSQU.ALBEDO.PY
      emsqu$rshortup  [m] =   mymont$MMSQU.RSHORTUP.PY
      emsqu$rlongup   [m] =   mymont$MMSQU.RLONGUP.PY
      emsqu$parup     [m] =   mymont$MMSQU.PARUP.PY * Watts.2.Ein^2 * 1.e12
      emsqu$rnet      [m] =   mymont$MMSQU.RNET.PY
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #       Read the mean diurnal cycle and the mean sum of the squares.                 #
      #------------------------------------------------------------------------------------#
      qmean$gpp         [m,] =   mymont$QMEAN.GPP.PY
      qmean$plant.resp  [m,] =   mymont$QMEAN.PLRESP.PY
      qmean$npp         [m,] =   mymont$QMEAN.NPP.PY
      qmean$leaf.resp   [m,] =   mymont$QMEAN.LEAF.RESP.PY
      qmean$root.resp   [m,] =   mymont$QMEAN.ROOT.RESP.PY
      qmean$het.resp    [m,] =   mymont$QMEAN.RH.PY
      qmean$cwd.resp    [m,] =   mymont$QMEAN.CWD.RH.PY
      qmean$nep         [m,] =   mymont$QMEAN.NEP.PY
      qmean$nee         [m,] = ( mymont$QMEAN.CARBON.ST.PY
                               - mymont$QMEAN.CARBON.AC.PY )
      qmean$reco        [m,] =   mymont$QMEAN.PLRESP.PY + mymont$QMEAN.RH.PY
      qmean$cflxca      [m,] = - mymont$QMEAN.CARBON.AC.PY
      qmean$cflxst      [m,] = - mymont$QMEAN.CARBON.ST.PY
      qmean$hflxca      [m,] = - mymont$QMEAN.SENSIBLE.AC.PY
      qmean$hflxlc      [m,] =   mymont$QMEAN.SENSIBLE.LC.PY
      qmean$hflxwc      [m,] =   mymont$QMEAN.SENSIBLE.WC.PY
      qmean$hflxgc      [m,] =   mymont$QMEAN.SENSIBLE.GC.PY
      qmean$wflxca      [m,] = - mymont$QMEAN.VAPOR.AC.PY     * day.sec
      qmean$qwflxca     [m,] = - mymont$QMEAN.VAPOR.AC.PY     * qmean.can.alvli.py
      qmean$wflxlc      [m,] =   mymont$QMEAN.VAPOR.LC.PY     * day.sec
      qmean$wflxwc      [m,] =   mymont$QMEAN.VAPOR.WC.PY     * day.sec
      qmean$wflxgc      [m,] =   mymont$QMEAN.VAPOR.GC.PY     * day.sec
      qmean$runoff      [m,] = ( mymont$QMEAN.RUNOFF.PY
                               + mymont$QMEAN.DRAINAGE.PY )   * day.sec
      qmean$intercepted [m,] = ( mymont$QMEAN.INTERCEPTED.AL.PY
                               + mymont$QMEAN.INTERCEPTED.AW.PY ) * day.sec
      qmean$wshed       [m,] = ( mymont$QMEAN.WSHED.LG.PY
                               + mymont$QMEAN.WSHED.WG.PY       ) * day.sec
      qmean$evap        [m,] = ( mymont$QMEAN.VAPOR.GC.PY
                               + mymont$QMEAN.VAPOR.WC.PY
                               + mymont$QMEAN.VAPOR.LC.PY )   * day.sec
      qmean$transp      [m,] =   mymont$QMEAN.TRANSP.PY       * day.sec
      qmean$atm.temp    [m,] =   mymont$QMEAN.ATM.TEMP.PY     - t00
      qmean$can.temp    [m,] =   mymont$QMEAN.CAN.TEMP.PY     - t00
      qmean$leaf.temp   [m,] =   mymont$QMEAN.LEAF.TEMP.PY    - t00
      qmean$wood.temp   [m,] =   mymont$QMEAN.WOOD.TEMP.PY    - t00
      qmean$gnd.temp    [m,] =   mymont$QMEAN.GND.TEMP.PY     - t00
      qmean$atm.shv     [m,] =   mymont$QMEAN.ATM.SHV.PY      * kg2g
      qmean$can.shv     [m,] =   mymont$QMEAN.CAN.SHV.PY      * kg2g
      qmean$gnd.shv     [m,] =   mymont$QMEAN.GND.SHV.PY      * kg2g
      qmean$atm.vpd     [m,] =   mymont$QMEAN.ATM.VPDEF.PY    * 0.01
      qmean$can.vpd     [m,] =   mymont$QMEAN.CAN.VPDEF.PY    * 0.01
      qmean$leaf.vpd    [m,] =   mymont$QMEAN.LEAF.VPDEF.PY   * 0.01
      qmean$atm.co2     [m,] =   mymont$QMEAN.ATM.CO2.PY
      qmean$can.co2     [m,] =   mymont$QMEAN.CAN.CO2.PY
      qmean$atm.vels    [m,] =   mymont$QMEAN.ATM.VELS.PY
      qmean$ustar       [m,] =   mymont$QMEAN.USTAR.PY
      qmean$atm.prss    [m,] =   mymont$QMEAN.ATM.PRSS.PY     * 0.01
      qmean$can.prss    [m,] =   mymont$QMEAN.CAN.PRSS.PY     * 0.01
      qmean$fs.open     [m,] =   mymont$QMEAN.FS.OPEN.PY
      qmean$rain        [m,] =   mymont$QMEAN.PCPG.PY         * day.sec
      qmean$rshort      [m,] =   mymont$QMEAN.ATM.RSHORT.PY
      qmean$rshort.beam [m,] = ( mymont$QMEAN.ATM.RSHORT.PY 
                               - mymont$QMEAN.ATM.RSHORT.DIFF.PY )
      qmean$rshort.diff [m,] =   mymont$QMEAN.ATM.RSHORT.DIFF.PY
      qmean$rshort.gnd  [m,] =   mymont$QMEAN.RSHORT.GND.PY
      qmean$rshortup    [m,] =   mymont$QMEAN.RSHORTUP.PY
      qmean$rlong       [m,] =   mymont$QMEAN.ATM.RLONG.PY
      qmean$rlong.gnd   [m,] =   mymont$QMEAN.RLONG.GND.PY
      qmean$rlongup     [m,] =   mymont$QMEAN.RLONGUP.PY
      qmean$par.tot     [m,] =   mymont$QMEAN.ATM.PAR.PY        * Watts.2.Ein * 1.e6
      qmean$par.beam    [m,] = ( mymont$QMEAN.ATM.PAR.PY
                               - mymont$QMEAN.ATM.PAR.DIFF.PY ) * Watts.2.Ein * 1.e6
      qmean$par.diff    [m,] =   mymont$QMEAN.ATM.PAR.DIFF.PY   * Watts.2.Ein * 1.e6
      qmean$par.gnd     [m,] =   mymont$QMEAN.PAR.GND.PY        * Watts.2.Ein * 1.e6
      qmean$parup       [m,] =   mymont$QMEAN.PARUP.PY          * Watts.2.Ein * 1.e6
      qmean$rnet        [m,] =   mymont$QMEAN.RNET.PY
      qmean$albedo      [m,] =   mymont$QMEAN.ALBEDO.PY
      if (all(c("QMEAN.ALBEDO.PAR.PY","QMEAN.ALBEDO.NIR.PY") %in% names(mymont))){
         qmean$albedo.par   [m,] =   mymont$QMEAN.ALBEDO.PAR.PY
         qmean$albedo.nir   [m,] =   mymont$QMEAN.ALBEDO.NIR.PY
      }else{
         qmean$albedo.par   [m,] = ifelse( mymont$QMEAN.ATM.PAR.PY > 0.5
                                         , mymont$QMEAN.PARUP.PY / mymont$QMEAN.ATM.PAR.PY
                                         , mymont$QMEAN.ALBEDO.PY
                                         )#end ifelse
         qmean$albedo.nir   [m,] = ifelse( mymont$QMEAN.ATM.NIR.PY > 0.5
                                         , mymont$QMEAN.NIRUP.PY / mymont$QMEAN.ATM.NIR.PY
                                         , mymont$QMEAN.ALBEDO.PY
                                         )#end ifelse
      }#end if
      qmean$rlong.albedo[m,] =   mymont$QMEAN.RLONG.ALBEDO.PY
      qmean$leaf.gbw    [m,] =   mymont$QMEAN.LEAF.GBW.PY       * day.sec
      qmean$leaf.gsw    [m,] =   mymont$QMEAN.LEAF.GSW.PY       * day.sec
      qmean$wood.gbw    [m,] =   mymont$QMEAN.WOOD.GBW.PY       * day.sec
      #------------------------------------------------------------------------------------#



      #------ Read the mean sum of squares for diel. --------------------------------------#
      qmsqu$gpp         [m,] =   mymont$QMSQU.GPP.PY
      qmsqu$npp         [m,] =   mymont$QMSQU.NPP.PY
      qmsqu$plant.resp  [m,] =   mymont$QMSQU.PLRESP.PY
      qmsqu$het.resp    [m,] =   mymont$QMSQU.RH.PY
      qmsqu$cwd.resp    [m,] =   mymont$QMSQU.CWD.RH.PY
      qmsqu$nep         [m,] =   mymont$QMSQU.NEP.PY
      qmsqu$cflxca      [m,] =   mymont$QMSQU.CARBON.AC.PY
      qmsqu$cflxst      [m,] =   mymont$QMSQU.CARBON.ST.PY
      qmsqu$hflxca      [m,] =   mymont$QMSQU.SENSIBLE.AC.PY
      qmsqu$hflxlc      [m,] =   mymont$QMSQU.SENSIBLE.LC.PY
      qmsqu$hflxwc      [m,] =   mymont$QMSQU.SENSIBLE.WC.PY
      qmsqu$hflxgc      [m,] =   mymont$QMSQU.SENSIBLE.GC.PY
      qmsqu$wflxca      [m,] =   mymont$QMSQU.VAPOR.AC.PY  * day.sec2
      qmsqu$qwflxca     [m,] =   mymont$QMSQU.VAPOR.AC.PY  * qmean.can.alvli.py2
      qmsqu$wflxlc      [m,] =   mymont$QMSQU.VAPOR.WC.PY  * day.sec2
      qmsqu$wflxwc      [m,] =   mymont$QMSQU.VAPOR.LC.PY  * day.sec2
      qmsqu$wflxgc      [m,] =   mymont$QMSQU.VAPOR.GC.PY  * day.sec2
      qmsqu$transp      [m,] =   mymont$QMSQU.TRANSP.PY    * day.sec2
      qmsqu$ustar       [m,] =   mymont$QMSQU.USTAR.PY
      qmsqu$albedo      [m,] =   mymont$QMSQU.ALBEDO.PY
      qmsqu$rshortup    [m,] =   mymont$QMSQU.RSHORTUP.PY
      qmsqu$rlongup     [m,] =   mymont$QMSQU.RLONGUP.PY
      qmsqu$parup       [m,] =   mymont$QMSQU.PARUP.PY     * Watts.2.Ein^2 * 1e12
      #------------------------------------------------------------------------------------#


      #---- Read in the site-level area. --------------------------------------------------#
      areasi     = mymont$AREA.SI
      npatches   = mymont$SIPA.N
      #------------------------------------------------------------------------------------#


      #----- Read a few patch-level variables. --------------------------------------------#
      areapa      = mymont$AREA * rep(areasi,times=npatches)
      areapa      = areapa / sum(areapa)
      ipa         = sequence(mymont$NPATCHES.GLOBAL)
      lupa        = mymont$DIST.TYPE
      agepa       = mymont$AGE
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Get the water deficit, and the estimate using Malhi's ET (100mm/month).        #
      #------------------------------------------------------------------------------------#
      if (m == 1){
         emean$water.deficit [m] = max( 0., emean$wflxca[m] * mondays - emean$rain[m])
         emean$malhi.deficit [m] = max( 0., et.malhi - emean$rain [m])
      }else{
         emean$water.deficit [m] = max( 0., emean$water.deficit [m-1] 
                                          + emean$wflxca[m] * mondays - emean$rain[m])
         emean$malhi.deficit [m] = max( 0., emean$malhi.deficit [m-1] 
                                          + et.malhi - emean$rain [m] )
      }#end if
      #------------------------------------------------------------------------------------#


      #----- Keep the maximum water deficit of the past period. ---------------------------#
      emean$last.1yr.mwd [m] = max(emean$water.deficit[max(m-11,1):m])
      emean$last.2yr.mwd [m] = max(emean$water.deficit[max(m-23,1):m])
      emean$last.3yr.mwd [m] = max(emean$water.deficit[max(m-35,1):m])
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #    If this is a biomass initialisation, or a run with anthropogenic                #
      # disturbance, we must jitter the age so we can distinguish the patches.             #
      #------------------------------------------------------------------------------------#
      sameage        = duplicated(agepa)
      agepa[sameage] = jitter(x=agepa[sameage],amount=0.4)
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Get the total number of cohorts.                                               #
      #------------------------------------------------------------------------------------#
      ncohorts    = mymont$PACO.N
      ipaconow    = rep(sequence(mymont$NPATCHES.GLOBAL),times=mymont$PACO.N)
      icoconow    = unlist(sapply(X = mymont$PACO.N, FUN = sequence))
      idx         = match(unique(ipaconow),sequence(mymont$NPATCHES.GLOBAL))
      #------------------------------------------------------------------------------------#



      #----- Disturbance rates. -----------------------------------------------------------#
      lu$dist  [m,,] = apply ( X      = mymont$DISTURBANCE.RATES
                             , MARGIN = c(2,3)
                             , FUN    = weighted.mean
                             , w      = areasi
                             )#end apply
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #      Build the cohort-level lists if this is the right month.                      #
      #------------------------------------------------------------------------------------#
      plab = paste( "y",sprintf("%4.4i",thisyear )
                  , "m",sprintf("%2.2i",thismonth),sep="")
      #----- Bind the current patches. ----------------------------------------------------#
      patch$ipa          [[plab]] =   ipa
      patch$age          [[plab]] =   agepa
      patch$area         [[plab]] =   areapa
      patch$lu           [[plab]] =   lupa
      patch$nep          [[plab]] =   mymont$MMEAN.NEP.PA
      patch$het.resp     [[plab]] =   mymont$MMEAN.RH.PA
      patch$can.temp     [[plab]] =   mymont$MMEAN.CAN.TEMP.PA    - t00
      patch$gnd.temp     [[plab]] =   mymont$MMEAN.GND.TEMP.PA    - t00
      patch$can.shv      [[plab]] =   mymont$MMEAN.CAN.SHV.PA     * 1000.
      patch$gnd.shv      [[plab]] =   mymont$MMEAN.GND.SHV.PA     * 1000.
      patch$can.vpd      [[plab]] =   mymont$MMEAN.CAN.VPDEF.PA   * 0.01
      patch$can.co2      [[plab]] =   mymont$MMEAN.CAN.CO2.PA
      patch$can.prss     [[plab]] =   mymont$MMEAN.CAN.PRSS.PA    * 0.01
      patch$cflxca       [[plab]] = - mymont$MMEAN.CARBON.AC.PA
      patch$cflxst       [[plab]] =   mymont$MMEAN.CARBON.ST.PA
      patch$nee          [[plab]] = ( mymont$MMEAN.CARBON.ST.PA
                                    - mymont$MMEAN.CARBON.AC.PA )
      patch$hflxca       [[plab]] = - mymont$MMEAN.SENSIBLE.AC.PA
      patch$hflxgc       [[plab]] =   mymont$MMEAN.SENSIBLE.GC.PA
      patch$qwflxca      [[plab]] = - mymont$MMEAN.VAPOR.AC.PA    * mmean.can.alvli.pa
      patch$wflxca       [[plab]] = - mymont$MMEAN.VAPOR.AC.PA
      patch$wflxgc       [[plab]] =   mymont$MMEAN.VAPOR.GC.PA
      patch$ustar        [[plab]] =   mymont$MMEAN.USTAR.PA
      patch$albedo       [[plab]] =   mymont$MMEAN.ALBEDO.PA
      patch$rshortup     [[plab]] =   mymont$MMEAN.RSHORTUP.PA
      patch$rlongup      [[plab]] =   mymont$MMEAN.RLONGUP.PA
      patch$parup        [[plab]] =   mymont$MMEAN.PARUP.PA       * Watts.2.Ein * 1e6
      patch$rnet         [[plab]] =   mymont$MMEAN.RNET.PA
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Find the patch-level properties that are derived from cohort-level.            #
      #------------------------------------------------------------------------------------#
      patch$lai       [[plab]] = rep(0.,times=mymont$NPATCHES.GLOBAL)
      patch$wai       [[plab]] = rep(0.,times=mymont$NPATCHES.GLOBAL)
      patch$leaf.temp [[plab]] = mymont$MMEAN.CAN.TEMP.PA  - t00
      patch$leaf.vpd  [[plab]] = mymont$MMEAN.CAN.VPDEF.PA * 0.01
      patch$wood.temp [[plab]] = mymont$MMEAN.CAN.TEMP.PA  - t00
      patch$gpp       [[plab]] = rep(0.,times=mymont$NPATCHES.GLOBAL)
      patch$npp       [[plab]] = rep(0.,times=mymont$NPATCHES.GLOBAL)
      patch$plant.resp[[plab]] = rep(0.,times=mymont$NPATCHES.GLOBAL)
      patch$leaf.temp [[plab]] = mymont$MMEAN.CAN.TEMP.PA - t00
      patch$hflxlc    [[plab]] = rep(0.,times=mymont$NPATCHES.GLOBAL)
      patch$hflxwc    [[plab]] = rep(0.,times=mymont$NPATCHES.GLOBAL)
      patch$wflxlc    [[plab]] = rep(0.,times=mymont$NPATCHES.GLOBAL)
      patch$wflxwc    [[plab]] = rep(0.,times=mymont$NPATCHES.GLOBAL)
      patch$transp    [[plab]] = rep(0.,times=mymont$NPATCHES.GLOBAL)


      if (any(ncohorts >0)){
         #----- Find some auxiliary patch-level properties. -------------------------------#
         lai.pa         = tapply( X     = mymont$MMEAN.LAI.CO
                                , INDEX = ipaconow
                                , FUN   = sum
                                )#end tapply
         wai.pa         = tapply( X     = mymont$WAI.CO
                                , INDEX = ipaconow
                                , FUN   = sum
                                )#end tapply
         leaf.energy.pa = tapply( X     = mymont$MMEAN.LEAF.ENERGY.CO
                                , INDEX = ipaconow
                                , FUN   = sum
                                )#end tapply
         leaf.water.pa  = tapply( X     = mymont$MMEAN.LEAF.WATER.CO
                                , INDEX = ipaconow
                                , FUN   = sum
                                )#end tapply
         leaf.hcap.pa   = tapply( X     = mymont$MMEAN.LEAF.HCAP.CO
                                , INDEX = ipaconow
                                , FUN   = sum
                                )#end tapply
         wood.energy.pa = tapply( X     = mymont$MMEAN.WOOD.ENERGY.CO
                                , INDEX = ipaconow
                                , FUN   = sum
                                )#end tapply
         wood.water.pa  = tapply( X     = mymont$MMEAN.WOOD.WATER.CO
                                , INDEX = ipaconow
                                , FUN   = sum
                                )#end tapply
         wood.hcap.pa   = tapply( X     = mymont$MMEAN.WOOD.HCAP.CO
                                , INDEX = ipaconow
                                , FUN   = sum
                                )#end tapply
         #---------------------------------------------------------------------------------#


         #----- Find the temperature and liquid fraction of leaf and wood. ----------------#
         leaf.empty               = leaf.hcap.pa == 0
         wood.empty               = wood.hcap.pa == 0
         leaf.temp.pa             = uextcm2tl( uext    = leaf.energy.pa
                                             , wmass   = leaf.water.pa
                                             , dryhcap = leaf.hcap.pa   )$temp - t00
         wood.temp.pa             = uextcm2tl( uext    = wood.energy.pa
                                             , wmass   = wood.water.pa
                                             , dryhcap = wood.hcap.pa   )$temp - t00
         leaf.temp.pa[leaf.empty] = NA
         wood.temp.pa[wood.empty] = NA
         #---------------------------------------------------------------------------------#





         #----- Find the variables that must be rendered extensive. -----------------------#
         gpp.pa        = tapply( X     = mymont$MMEAN.GPP.CO * mymont$NPLANT
                               , INDEX = ipaconow
                               , FUN   = sum
                               )#end tapply
         npp.pa        = tapply( X     = mymont$MMEAN.NPP.CO * mymont$NPLANT
                               , INDEX = ipaconow
                               , FUN   = sum
                               )#end tapply
         plant.resp.pa = tapply( X     = mymont$MMEAN.PLRESP.CO * mymont$NPLANT
                               , INDEX = ipaconow
                               , FUN   = sum
                               )#end tapply
         #---------------------------------------------------------------------------------#





         #----- Add the variables that are already extensive. -----------------------------#
         hflxlc.pa = tapply( X     = mymont$MMEAN.SENSIBLE.LC.CO
                           , INDEX = ipaconow
                           , FUN   = sum
                           )#end tapply
         hflxwc.pa = tapply( X     = mymont$MMEAN.SENSIBLE.WC.CO
                           , INDEX = ipaconow
                           , FUN   = sum
                           )#end tapply
         wflxlc.pa = tapply( X     = mymont$MMEAN.VAPOR.LC.CO  * day.sec
                           , INDEX = ipaconow
                           , FUN   = sum
                           )#end tapply
         wflxwc.pa = tapply( X     = mymont$MMEAN.VAPOR.WC.CO  * day.sec
                           , INDEX = ipaconow
                           , FUN   = sum
                           )#end tapply
         transp.pa = tapply( X     = mymont$MMEAN.TRANSP.CO    * day.sec
                           , INDEX = ipaconow
                           , FUN   = sum
                           )#end tapply
         #---------------------------------------------------------------------------------#





         #---------------------------------------------------------------------------------#
         #      Vapour pressure deficit is found using weighted averages.                  #
         #---------------------------------------------------------------------------------#
         leaf.vpd.pa = mapply( FUN      = weighted.mean
                             , x        = split(mymont$MMEAN.LEAF.VPDEF.CO,ipaconow)
                             , w        = split(mymont$MMEAN.LAI.CO       ,ipaconow)
                             , SIMPLIFY = TRUE
                             )#end mapply
         leaf.vpd.pa[leaf.empty] = NA
         leaf.vpd.pa             = leaf.vpd.pa * 0.01
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #     Copy the data back to the patch.                                            #
         #---------------------------------------------------------------------------------#
         patch$lai       [[plab]][idx              ] = lai.pa
         patch$wai       [[plab]][idx              ] = wai.pa
         patch$leaf.temp [[plab]][idx[! leaf.empty]] = leaf.temp.pa[! leaf.empty]
         patch$leaf.vpd  [[plab]][idx[! leaf.empty]] = leaf.vpd.pa [! leaf.empty]
         patch$wood.temp [[plab]][idx[! wood.empty]] = wood.temp.pa[! wood.empty]
         patch$gpp       [[plab]][idx              ] = gpp.pa
         patch$npp       [[plab]][idx              ] = npp.pa
         patch$plant.resp[[plab]][idx              ] = npp.pa
         patch$hflxlc    [[plab]][idx              ] = hflxlc.pa
         patch$hflxwc    [[plab]][idx              ] = hflxwc.pa
         patch$wflxlc    [[plab]][idx              ] = wflxlc.pa
         patch$wflxwc    [[plab]][idx              ] = wflxwc.pa
         patch$transp    [[plab]][idx              ] = transp.pa
         #---------------------------------------------------------------------------------#
      }#end if
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Ecosystem respiration, which is a combination of plant respiration (cohort-    #
      # -based) and heterotrophic respiration (patch-based).                               #
      #------------------------------------------------------------------------------------#
      patch$reco[[plab]] = patch$plant.resp[[plab]] + patch$het.resp[[plab]]
      #------------------------------------------------------------------------------------#




      #------------------------------------------------------------------------------------#
      #     Read the cohort-level variables.  Because empty patchs do exist (deserts),     #
      # we must check whether there is any cohort to be read.  If not, assign NA to        #
      # all variables.                                                                     #
      #------------------------------------------------------------------------------------#
      if (any (ncohorts > 0)){




         areaconow  = rep(areapa,times=ncohorts)

         #----- Define the land use classes. ----------------------------------------------#
         luconow    = rep(lupa,times=ncohorts)

         #----- Define the DBH classes. ---------------------------------------------------#
         dbhconow        = mymont$DBH
         dbhcut          = cut(dbhconow,breaks=breakdbh)
         dbhlevs         = levels(dbhcut)
         dbhfac          = match(dbhcut,dbhlevs)
         #---------------------------------------------------------------------------------#



         #----- Define the previous DBH class (for recruitment). --------------------------#
         dbhconow.lastmon = mymont$DBH * exp(-pmax(0,mymont$DLNDBH.DT/12))
         dbhconow.1ago    = mymont$DBH * exp(-pmax(0,mymont$DLNDBH.DT))
         dbhcut.1ago      = cut(dbhconow.1ago,breaks=breakdbh)
         dbhlevs.1ago     = levels(dbhcut.1ago)
         dbhfac.1ago      = match(dbhcut.1ago,dbhlevs.1ago)
         #---------------------------------------------------------------------------------#


         #----- Define the age classes. ---------------------------------------------------#
         ageconow          = rep(x=agepa,times=ncohorts)
         #---------------------------------------------------------------------------------#



         #----- Read the cohort level variables. ------------------------------------------#
         pftconow          = mymont$PFT
         nplantconow       = mymont$NPLANT
         heightconow       = mymont$HITE
         baconow           = mymont$BA.CO
         agbconow          = mymont$AGB.CO
         laiconow          = mymont$MMEAN.LAI.CO
         waiconow          = mymont$WAI.CO
         taiconow          = laiconow + waiconow
         gppconow          = mymont$MMEAN.GPP.CO
         leaf.respconow    = mymont$MMEAN.LEAF.RESP.CO
         root.respconow    = mymont$MMEAN.ROOT.RESP.CO
         growth.respconow  = mymont$MMEAN.GROWTH.RESP.CO
         plant.respconow   = ( mymont$MMEAN.LEAF.RESP.CO   + mymont$MMEAN.ROOT.RESP.CO
                             + mymont$MMEAN.GROWTH.RESP.CO + mymont$MMEAN.STORAGE.RESP.CO
                             + mymont$MMEAN.VLEAF.RESP.CO  )
         nppconow          = gppconow-plant.respconow
         
         if (kludgecbal){
            cbaconow          = mymont$MMEAN.CB - (mondays - 1) * mymont$MMEAN.BSTORAGE.CO
            #---------------------------------------------------------------------------------#
            #     Temporary fix to correct the carbon balance.                                #
            #---------------------------------------------------------------------------------#
            if (sum(ncohorts) == 1){
               cbalightconow     = ( mymont$CB.LIGHTMAX[thismonth] 
                                   - (mondays - 1) * mymont$MMEAN.BSTORAGE.CO )
               cbamoistconow     = ( mymont$CB.MOISTMAX[thismonth]
                                   - (mondays - 1) * mymont$MMEAN.BSTORAGE.CO )
            }else{
               cbalightconow     = ( mymont$CB.LIGHTMAX[,thismonth] 
                                   - (mondays - 1) * mymont$MMEAN.BSTORAGE.CO )
               cbamoistconow     = ( mymont$CB.MOISTMAX[,thismonth]
                                   - (mondays - 1) * mymont$MMEAN.BSTORAGE.CO )
            }#end if
            #---------------------------------------------------------------------------------#
         }else{
            cbaconow             = mymont$MMEAN.CB
            #---------------------------------------------------------------------------------#
            #     Temporary fix to correct the carbon balance.                                #
            #---------------------------------------------------------------------------------#
            if (sum(ncohorts) == 1){
               cbalightconow     = mymont$CB.LIGHTMAX[thismonth]
               cbamoistconow     = mymont$CB.MOISTMAX[thismonth]
            }else{
               cbalightconow     = mymont$CB.LIGHTMAX[,thismonth]
               cbamoistconow     = mymont$CB.MOISTMAX[,thismonth]
            }#end if
            #---------------------------------------------------------------------------------#
         }#end if

         cbamaxconow       = klight * cbalightconow + (1. - klight) * cbamoistconow
         cbarelconow       = mymont$CBR.BAR
         mcostconow        = ( mymont$MMEAN.LEAF.MAINTENANCE.CO
                             + mymont$MMEAN.ROOT.MAINTENANCE.CO ) * yr.day
         ldropconow        = mymont$MMEAN.LEAF.DROP.CO * yr.day
         fs.openconow      = mymont$MMEAN.FS.OPEN.CO
         lightconow        = mymont$MMEAN.LIGHT.LEVEL.CO
         light.beamconow   = mymont$MMEAN.LIGHT.LEVEL.BEAM.CO
         light.diffconow   = mymont$MMEAN.LIGHT.LEVEL.DIFF.CO

         bdeadconow        = mymont$BDEAD
         bleafconow        = mymont$MMEAN.BLEAF.CO
         brootconow        = mymont$MMEAN.BROOT.CO
         bsapwoodconow     = mymont$BSAPWOODA+mymont$BSAPWOODB
         baliveconow       = bleafconow + brootconow + bsapwoodconow
         bstorageconow     = mymont$MMEAN.BSTORAGE.CO
         bseedsconow       = mymont$BSEEDS.CO
         bgbconow          = ( brootconow + mymont$BSAPWOODB
                             + (1. - pft$agf.bs[pftconow]) * bdeadconow )
         biomassconow      = baliveconow + bstorageconow + bseedsconow + bdeadconow
         #---------------------------------------------------------------------------------#



         #----- Allocation and productivity relative to the total biomass. ----------------#
         f.gppconow        =  100. * gppconow        / biomassconow
         f.plant.respconow =  100. * plant.respconow / biomassconow
         f.nppconow        =  100. * nppconow        / biomassconow
         f.cbaconow        =  100. * cbaconow        / biomassconow
         f.bstorageconow   = 1000. * bstorageconow   / biomassconow
         f.bleafconow      = 1000. * bleafconow      / biomassconow
         f.brootconow      = 1000. * brootconow      / biomassconow
         f.bseedsconow     = 1000. * bseedsconow     / biomassconow
         #---------------------------------------------------------------------------------#



         #----- Energy and water fluxes: convert them to plant area. ----------------------#
         hflxlcconow       = mymont$MMEAN.SENSIBLE.LC.CO           / nplantconow
         wflxlcconow       = mymont$MMEAN.VAPOR.LC.CO    * day.sec / nplantconow
         transpconow       = mymont$MMEAN.TRANSP.CO      * day.sec / nplantconow
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #     Find the leaf interstitial space and boundary layer specific humidities to  #
         # convert conductance to kgW/m2/day.                                              #
         #---------------------------------------------------------------------------------#
         lpsiconow      = ifelse(laiconow > 1.e-10,mymont$MMEAN.TRANSP.CO/laiconow,0.)
         wpsiconow      = ifelse(waiconow > 1.e-10,mymont$MMEAN.WFLXWC.CO/waiconow,0.)
         can.shv.conow  = rep(mymont$MMEAN.CAN.SHV.PA,times=ncohorts)
         #---- Net conductance, combining stomatal and boundary layer. --------------------#
         leaf.gnw.conow = ifelse( mymont$MMEAN.LEAF.GBW.CO+mymont$MMEAN.LEAF.GSW.CO>1.e-10
                                ,   mymont$MMEAN.LEAF.GBW.CO * mymont$MMEAN.LEAF.GSW.CO
                                / ( mymont$MMEAN.LEAF.GBW.CO + mymont$MMEAN.LEAF.GSW.CO )
                                , pmin(mymont$MMEAN.LEAF.GBW.CO,mymont$MMEAN.LEAF.GSW.CO)
                                )#end ifelse
         lbl.shv.conow  = can.shv.conow + lpsiconow / pmax(mymont$MMEAN.LEAF.GBW.CO,1.e-10)
         wbl.shv.conow  = can.shv.conow + wpsiconow / pmax(mymont$MMEAN.WOOD.GBW.CO,1.e-10)
         lis.shv.conow  = can.shv.conow + lpsiconow / pmax(leaf.gnw.conow,1.e-10)
         #---------------------------------------------------------------------------------#


         #----- Find the conductances in kgW/m2/day. --------------------------------------#
         leaf.gbwconow  = ( mymont$MMEAN.LEAF.GBW.CO * day.sec * ep
                          * (1 + epim1 * can.shv.conow) * (1 + epim1 * lbl.shv.conow) )
         leaf.gswconow  = ( mymont$MMEAN.LEAF.GSW.CO * day.sec * ep
                          * (1 + epim1 * lbl.shv.conow) * (1 + epim1 * lis.shv.conow) )
         wood.gbwconow  = ( mymont$MMEAN.WOOD.GBW.CO * day.sec * ep
                          * (1 + epim1 * can.shv.conow) * (1 + epim1 * wbl.shv.conow) )
         #---------------------------------------------------------------------------------#


         #----- Find the net radiation for leaves (in m2/leaf!). --------------------------#
         leaf.parconow    = ifelse( laiconow > 1.e-10
                                  , mymont$MMEAN.PAR.L.CO / laiconow *  Watts.2.Ein
                                  , 0.
                                  )#end ifelse
         leaf.rshortconow = ifelse( laiconow > 1.e-10
                                  , mymont$MMEAN.RSHORT.L.CO / laiconow
                                  , 0.
                                  )#end ifelse
         leaf.rlongconow  = ifelse( laiconow > 1.e-10
                                  , mymont$MMEAN.RLONG.L.CO  / laiconow
                                  , 0.
                                  )#end ifelse
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Leaf/wood thermal properties.                                               #
         #---------------------------------------------------------------------------------#
         leaf.tempconow    = mymont$MMEAN.LEAF.TEMP.CO  - t00
         wood.tempconow    = mymont$MMEAN.WOOD.TEMP.CO  - t00
         leaf.vpdconow     = mymont$MMEAN.LEAF.VPDEF.CO * 0.01
         #---------------------------------------------------------------------------------#


         #------ Find the demand and supply by m2gnd. -------------------------------------#
         demandconow       = mymont$MMEAN.PSI.OPEN.CO     * laiconow * day.sec
         supplyconow       = mymont$MMEAN.WATER.SUPPLY.CO * day.sec
         #---------------------------------------------------------------------------------#


         #------ Find the demographic rates. ----------------------------------------------#
         mortconow       = pmax(0,rowSums(mymont$MMEAN.MORT.RATE.CO))
         ncbmortconow    = pmax(0,mymont$MMEAN.MORT.RATE.CO[,2])
         dimortconow     = pmax(0,mortconow - ncbmortconow)
         recruitconow    = mymont$RECRUIT.DBH
         growthconow     = 100. * pmax(0,mymont$DLNDBH.DT)
         agb.growthconow = 100. * pmax(0,mymont$DLNAGB.DT)
         bsa.growthconow = 100. * pmax(0,mymont$DLNBA.DT )
         #---------------------------------------------------------------------------------#


         #------ Find ANPP, BNPP and rainfall use efficiencies. ---------------------------#
         oldconow        = dbhconow.lastmon >= recruit.dbh.min
         agb.deltaconow  = ( mortconow - 0.01 * agb.growthconow ) / 12.
         bsa.deltaconow  = ( mortconow - 0.01 * bsa.growthconow ) / 12.
         
         anppconow       = 12. * agbconow * ( 1. - oldconow * exp(agb.deltaconow) )
         bnppconow       = 12. * baconow  * ( 1. - oldconow * exp(bsa.deltaconow) )
         arueconow       = 1000. * anppconow / max(12. * emean$rain[m],0.01)
         brueconow       = 100.  * bnppconow / max(12. * emean$rain[m],0.01)
         #---------------------------------------------------------------------------------#



         #----- Find the water use efficiency. --------------------------------------------#
         wueconow = ifelse(transpconow > 1.e-10, 1000. * anppconow/transpconow/yr.day, 0.)
         #---------------------------------------------------------------------------------#




      }else{
         #----- Make everything NA. -------------------------------------------------------#
         ipaconow          = NA
         icoconow          = NA
         areaconow         = NA
         luconow           = NA
         dbhconow          = NA
         dbhcut            = NA
         dbhlevs           = NA
         dbhfac            = NA
         dbhconow.1ago     = NA
         dbhcut.1ago       = NA
         dbhlevs.1ago      = NA
         dbhfac.1ago       = NA
         ageconow          = NA
         pftconow          = NA
         nplantconow       = NA
         heightconow       = NA
         baconow           = NA
         agbconow          = NA
         bgbconow          = NA
         biomassconow      = NA
         laiconow          = NA
         waiconow          = NA
         taiconow          = NA
         gppconow          = NA
         leaf.respconow    = NA
         root.respconow    = NA
         growth.respconow  = NA
         plant.respconow   = NA
         nppconow          = NA
         cbaconow          = NA
         cbamaxconow       = NA
         cbalightconow     = NA
         cbamoistconow     = NA
         cbarelconow       = NA
         mcostconow        = NA
         ldropconow        = NA
         fs.openconow      = NA
         lightconow        = NA
         light.beamconow   = NA
         light.diffconow   = NA
         baliveconow       = NA
         bdeadconow        = NA
         bleafconow        = NA
         brootconow        = NA
         bsapwoodconow     = NA
         bstorageconow     = NA
         bseedsconow       = NA
         hflxlcconow       = NA
         wflxlcconow       = NA
         transpconow       = NA
         wueconow          = NA
         leaf.tempconow    = NA
         wood.tempconow    = NA
         leaf.vpdconow     = NA
         demandconow       = NA
         supplyconow       = NA
         mortconow         = NA
         ncbmortconow      = NA
         dimortconow       = NA
         recruitconow      = NA
         growthconow       = NA
         agb.growthconow   = NA
         bsa.growthconow   = NA
         leaf.gbwconow     = NA
         leaf.gswconow     = NA
         wood.gbwconow     = NA
         f.gppconow        = NA
         f.plant.respconow = NA
         f.nppconow        = NA
         f.cbaconow        = NA
         f.bstorageconow   = NA
         f.bleafconow      = NA
         f.brootconow      = NA
         f.bseedsconow     = NA
         leaf.parconow     = NA
         leaf.rshortconow  = NA
         leaf.rlongconow   = NA
         anppconow         = NA
         bnppconow         = NA
         arueconow         = NA
         brueconow         = NA
      }#end if
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     The following two variables are used to scale "intensive" properties           #
      # (whatever/plant) to "extensive" (whatever/m2).  Sometimes they may be used to      #
      # build weighted averages.                                                           #
      #------------------------------------------------------------------------------------#
      w.nplant  = nplantconow  * areaconow
      w.lai     = laiconow     * areaconow
      w.wai     = waiconow     * areaconow
      w.tai     = taiconow     * areaconow
      w.biomass = biomassconow * nplantconow
      #------------------------------------------------------------------------------------#




      #------------------------------------------------------------------------------------#
      #     Build the LU arrays.                                                           #
      #------------------------------------------------------------------------------------#
      for (l in sequence(nlu+1)){
         selpa    = lupa    == l | l == (nlu+1)
         if (all(is.na(luconow))){
            sel      = rep(FALSE,times=length(luconow))
         }else{
            sel      = luconow == l | l == (nlu+1)
         }#end if
         
         if (any(sel)){
            lu$lai    [m,l] = sum( laiconow[sel] * areaconow    [sel] )
            lu$ba     [m,l] = sum( w.nplant[sel] * baconow      [sel] )
            lu$agb    [m,l] = sum( w.nplant[sel] * agbconow     [sel] )
            lu$bgb    [m,l] = sum( w.nplant[sel] * bgbconow     [sel] )
            lu$biomass[m,l] = sum( w.nplant[sel] * biomassconow [sel] )
            lu$gpp    [m,l] = sum( w.nplant[sel] * gppconow     [sel] )
            lu$npp    [m,l] = sum( w.nplant[sel] * nppconow     [sel] )
         }#end if
         lu$area      [m,l] = lu$area [m,l] + sum(areapa[selpa])
      }#end for
      #------------------------------------------------------------------------------------#




      #------------------------------------------------------------------------------------#
      #     Build the size (DBH) structure arrays.                                         #
      #------------------------------------------------------------------------------------#
      for (d in sequence(ndbh+1)){
         #----- Decide which DBH to use. --------------------------------------------------#
         if (all(is.na(dbhfac))){
            sel.dbh       = rep(FALSE,times=length(dbhfac     ))
            sel.dbh.1ago  = rep(FALSE,times=length(dbhfac.1ago))

            #----- Define the minimum DBH. ------------------------------------------------#
            dbhminconow   = rep(Inf,times=length(pftconow))
            #------------------------------------------------------------------------------#
         }else{
            sel.dbh       = dbhfac      == d | d == (ndbh+1)
            sel.dbh.1ago  = dbhfac.1ago == d | d == (ndbh+1)

            #----- Define the minimum DBH. ------------------------------------------------#
            dbhminconow   = pft$dbh.min[pftconow] * (d == 1) + census.dbh.min * (d != 1)
            #------------------------------------------------------------------------------#
         }#end if
         #---------------------------------------------------------------------------------#


         #----- Decide which PFT to use. --------------------------------------------------#
         for (p in sequence(npft+1)){
            sel.pft   = pftconow == p | p == (npft+1)
            sel       = sel.pft & sel.dbh
            if (any(sel)){
               #----- Extensive properties. -----------------------------------------------#
               szpft$lai        [m,d,p] = sum( laiconow         [sel] 
                                             * areaconow        [sel] )
               szpft$wai        [m,d,p] = sum( waiconow         [sel]  
                                             * areaconow        [sel] )
               szpft$tai        [m,d,p] = sum( taiconow         [sel]  
                                             * areaconow        [sel] )
               szpft$nplant     [m,d,p] = sum( nplantconow      [sel]  
                                             * areaconow        [sel] )
               szpft$demand     [m,d,p] = sum( demandconow      [sel]  
                                             * areaconow        [sel] )
               szpft$supply     [m,d,p] = sum( supplyconow      [sel]  
                                             * areaconow        [sel] )
               #----- Intensive properties, use nplant to make them extensive. ------------#
               szpft$agb        [m,d,p] = sum( w.nplant         [sel]  
                                             * agbconow         [sel] )
               szpft$bgb        [m,d,p] = sum( w.nplant         [sel]  
                                             * bgbconow         [sel] )
               szpft$biomass    [m,d,p] = sum( w.nplant         [sel]  
                                             * biomassconow     [sel] )
               szpft$ba         [m,d,p] = sum( w.nplant         [sel]  
                                             * baconow          [sel] )
               szpft$gpp        [m,d,p] = sum( w.nplant         [sel]  
                                             * gppconow         [sel] )
               szpft$npp        [m,d,p] = sum( w.nplant         [sel]  
                                             * nppconow         [sel] )
               szpft$mco        [m,d,p] = sum( w.nplant         [sel]  
                                             * mcostconow       [sel] )
               szpft$cba        [m,d,p] = sum( w.nplant         [sel]  
                                             * cbaconow         [sel] )
               szpft$cbamax     [m,d,p] = sum( w.nplant         [sel]  
                                             * cbamaxconow      [sel] )
               szpft$cbalight   [m,d,p] = sum( w.nplant         [sel]  
                                             * cbalightconow    [sel] )
               szpft$cbamoist   [m,d,p] = sum( w.nplant         [sel]  
                                             * cbamoistconow    [sel] )
               szpft$ldrop      [m,d,p] = sum( w.nplant         [sel]  
                                             * ldropconow       [sel] )
               szpft$leaf.resp  [m,d,p] = sum( w.nplant         [sel] 
                                             * leaf.respconow   [sel] )
               szpft$root.resp  [m,d,p] = sum( w.nplant         [sel] 
                                             * root.respconow   [sel] )
               szpft$growth.resp[m,d,p] = sum( w.nplant         [sel] 
                                             * growth.respconow [sel] )
               szpft$plant.resp [m,d,p] = sum( w.nplant         [sel] 
                                             * plant.respconow  [sel] )
               szpft$bdead      [m,d,p] = sum( w.nplant         [sel] 
                                             * bdeadconow       [sel] )
               szpft$bleaf      [m,d,p] = sum( w.nplant         [sel]
                                             * bleafconow       [sel] )
               szpft$broot      [m,d,p] = sum( w.nplant         [sel]
                                             * brootconow       [sel] )
               szpft$bsapwood   [m,d,p] = sum( w.nplant         [sel]
                                             * bsapwoodconow    [sel] )
               szpft$bstorage   [m,d,p] = sum( w.nplant         [sel] 
                                             * bstorageconow    [sel] )
               szpft$bseeds     [m,d,p] = sum( w.nplant         [sel] 
                                             * bseedsconow      [sel] )
               szpft$hflxlc     [m,d,p] = sum( w.nplant         [sel]
                                             * hflxlcconow      [sel] )
               szpft$wflxlc     [m,d,p] = sum( w.nplant         [sel]
                                             * wflxlcconow      [sel] )
               szpft$transp     [m,d,p] = sum( w.nplant         [sel]
                                             * transpconow      [sel] )
               szpft$anpp       [m,d,p] = sum( w.nplant         [sel]
                                             * anppconow        [sel] )
               szpft$bnpp       [m,d,p] = sum( w.nplant         [sel]
                                             * bnppconow        [sel] )
               szpft$arue       [m,d,p] = sum( w.nplant         [sel]
                                             * arueconow        [sel] )
               szpft$brue       [m,d,p] = sum( w.nplant         [sel]
                                             * brueconow        [sel] )
               #---------------------------------------------------------------------------#



               #----- Leaf/wood intensive properties , weighted means using LAI/WAI. ------#
               szpft$fs.open    [m,d,p] = weighted.mean( x = fs.openconow    [sel] 
                                                       , w = laiconow        [sel]
                                                       )#end weighted.mean
               szpft$leaf.par   [m,d,p] = weighted.mean( x = leaf.parconow   [sel] 
                                                       , w = laiconow        [sel]
                                                       )#end weighted.mean
               szpft$leaf.rshort[m,d,p] = weighted.mean( x = leaf.rshortconow[sel] 
                                                       , w = laiconow        [sel]
                                                       )#end weighted.mean
               szpft$leaf.rlong [m,d,p] = weighted.mean( x = leaf.rlongconow [sel] 
                                                       , w = laiconow        [sel]
                                                       )#end weighted.mean
               szpft$leaf.temp  [m,d,p] = weighted.mean( x = leaf.tempconow  [sel]
                                                       , w = laiconow        [sel]
                                                       )#end weighted.mean
               szpft$wood.temp  [m,d,p] = weighted.mean( x = wood.tempconow  [sel]
                                                       , w = waiconow        [sel]
                                                       )#end weighted.mean
               szpft$leaf.vpd   [m,d,p] = weighted.mean( x = leaf.vpdconow   [sel]
                                                       , w = laiconow        [sel]
                                                       )#end weighted.mean
               #---------------------------------------------------------------------------#


               #----- Individual-based properties, use weighted mean by Nplant. -----------#
               szpft$cbarel       [m,d,p] = weighted.mean( x = cbarelconow     [sel]
                                                         , w = w.nplant        [sel])
               szpft$i.gpp        [m,d,p] = weighted.mean( x = gppconow        [sel]
                                                         , w = w.nplant        [sel])
               szpft$i.npp        [m,d,p] = weighted.mean( x = nppconow        [sel]
                                                         , w = w.nplant        [sel])
               szpft$i.plant.resp [m,d,p] = weighted.mean( x = plant.respconow [sel]
                                                         , w = w.nplant        [sel])
               szpft$i.mco        [m,d,p] = weighted.mean( x = mcostconow      [sel]
                                                         , w = w.nplant        [sel])
               szpft$i.cba        [m,d,p] = weighted.mean( x = cbaconow        [sel]
                                                         , w = w.nplant        [sel])
               szpft$i.cbamax     [m,d,p] = weighted.mean( x = cbamaxconow     [sel]
                                                         , w = w.nplant        [sel])
               szpft$i.cbalight   [m,d,p] = weighted.mean( x = cbalightconow   [sel]
                                                         , w = w.nplant        [sel])
               szpft$i.cbamoist   [m,d,p] = weighted.mean( x = cbamoistconow   [sel]
                                                         , w = w.nplant        [sel])
               szpft$i.transp     [m,d,p] = weighted.mean( x = transpconow     [sel]
                                                         , w = w.nplant        [sel])
               szpft$i.wflxlc     [m,d,p] = weighted.mean( x = wflxlcconow     [sel]
                                                         , w = w.nplant        [sel])
               szpft$i.hflxlc     [m,d,p] = weighted.mean( x = hflxlcconow     [sel]
                                                         , w = w.nplant        [sel])
               szpft$leaf.gbw     [m,d,p] = weighted.mean( x = leaf.gbwconow   [sel]
                                                         , w = w.lai           [sel])
               szpft$leaf.gsw     [m,d,p] = weighted.mean( x = leaf.gswconow   [sel]
                                                         , w = w.lai           [sel])
               szpft$wood.gbw     [m,d,p] = weighted.mean( x = wood.gbwconow   [sel]
                                                         , w = w.wai           [sel])
               #---------------------------------------------------------------------------#


               #---------------------------------------------------------------------------#
               # WUE: use the bulk value because this is closer to what would be           #
               #      measured.                                                            #
               #---------------------------------------------------------------------------#
               szpft$wue          [m,d,p] = ifelse( szpft$transp[m,d,p] > 1.e-10
                                                  , 1000. * szpft$anpp  [m,d,p]
                                                          / szpft$transp[m,d,p]
                                                          / yr.day
                                                  , 0.
                                                  )#end ifelse
               #---------------------------------------------------------------------------#



               #----- Fractional biomass, use biomass as the weighting. -------------------#
               szpft$f.gpp        [m,d,p] = weighted.mean( x = f.gppconow        [sel]
                                                         , w = w.biomass         [sel]) 
               szpft$f.plant.resp [m,d,p] = weighted.mean( x = f.plant.respconow [sel]
                                                         , w = w.biomass         [sel]) 
               szpft$f.npp        [m,d,p] = weighted.mean( x = f.nppconow        [sel]
                                                         , w = w.biomass         [sel]) 
               szpft$f.cba        [m,d,p] = weighted.mean( x = f.cbaconow        [sel]
                                                         , w = w.biomass         [sel]) 
               szpft$f.bstorage   [m,d,p] = weighted.mean( x = f.bstorageconow   [sel]
                                                         , w = w.biomass         [sel]) 
               szpft$f.bleaf      [m,d,p] = weighted.mean( x = f.bleafconow      [sel]
                                                         , w = w.biomass         [sel]) 
               szpft$f.broot      [m,d,p] = weighted.mean( x = f.brootconow      [sel]
                                                         , w = w.biomass         [sel]) 
               szpft$f.bseeds     [m,d,p] = weighted.mean( x = f.bseedsconow     [sel]
                                                         , w = w.biomass         [sel]) 
               #---------------------------------------------------------------------------#
            }#end if
            #------------------------------------------------------------------------------#




            #------------------------------------------------------------------------------#
            #    For mortality and growth, we keep deleting the tiny guys because they     #
            # skew the rates quite significantly.                                          #
            #------------------------------------------------------------------------------#
            sel = sel.pft & sel.dbh & dbhconow >= dbhminconow
            if (any(sel)){
               #----- Growth rates are weighted by population. ----------------------------#
               szpft$growth     [m,d,p] = weighted.mean( x = growthconow    [sel]
                                                       , w = w.nplant       [sel]
                                                       )#end weighted.mean
               szpft$agb.growth [m,d,p] = weighted.mean( x = agb.growthconow[sel]
                                                       , w = w.nplant       [sel]
                                                       )#end weighted.mean
               szpft$bsa.growth [m,d,p] = weighted.mean( x = bsa.growthconow[sel]
                                                       , w = w.nplant       [sel]
                                                       )#end weighted.mean
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #      Find the total number of plants and previous population if the only  #
               # mortality was the mortality we test.                                      #
               #---------------------------------------------------------------------------#
               survivor             = sum( w.nplant[sel]                          )
               previous             = sum( w.nplant[sel] * exp(mortconow   [sel]) )
               ncb.previous         = sum( w.nplant[sel] * exp(ncbmortconow[sel]) )
               di.previous          = sum( w.nplant[sel] * exp(dimortconow [sel]) )
               szpft$mort   [m,d,p] = log( previous     / survivor )
               szpft$ncbmort[m,d,p] = log( ncb.previous / survivor )
               szpft$dimort [m,d,p] = log( di.previous  / survivor )
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #      Find the total AGB and previous AGB if the only mortality was the    #
               # mortality we test.                                                        #
               #---------------------------------------------------------------------------#
               survivor                 = sum( w.nplant[sel] * agbconow[sel])
               previous                 = sum( w.nplant[sel] * agbconow[sel]
                                             * exp(mortconow           [sel] ) )
               ncb.previous             = sum( w.nplant[sel] * agbconow[sel]
                                             * exp(ncbmortconow        [sel] ) )
               di.previous              = sum( w.nplant[sel] * agbconow[sel]
                                             * exp(dimortconow         [sel] ) )
               szpft$agb.mort   [m,d,p] = log( previous     / survivor )
               szpft$agb.ncbmort[m,d,p] = log( ncb.previous / survivor )
               szpft$agb.dimort [m,d,p] = log( di.previous  / survivor )
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #      Find the total basal area and previous basal area if the only        #
               # mortality was the mortality we test.                                      #
               #---------------------------------------------------------------------------#
               survivor                 = sum( w.nplant[sel] * baconow[sel])
               previous                 = sum( w.nplant[sel] * baconow[sel]
                                             * exp(mortconow          [sel] ) )
               ncb.previous             = sum( w.nplant[sel] * baconow[sel]
                                             * exp(ncbmortconow       [sel] ) )
               di.previous              = sum( w.nplant[sel] * baconow[sel]
                                             * exp(dimortconow        [sel] ) )
               szpft$bsa.mort   [m,d,p] = log( previous     / survivor )
               szpft$bsa.ncbmort[m,d,p] = log( ncb.previous / survivor )
               szpft$bsa.dimort [m,d,p] = log( di.previous  / survivor )
               #---------------------------------------------------------------------------#
            }#end if
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #    Recruitment: we must determine whether the plant grew into the new        #
            # category or not.                                                             #
            #------------------------------------------------------------------------------#
            sel.pop = sel.pft & sel.dbh      & dbhconow      >= dbhminconow
            sel.est = sel.pop & sel.dbh.1ago & dbhconow.1ago >= dbhminconow
            if (any(sel.pop) & any(sel.est)){
               #----- Recruitment rate in terms of individuals. ---------------------------#
               population             = sum(w.nplant[sel.pop])
               established            = sum(w.nplant[sel.est])
               szpft$recr     [m,d,p] = log(population / established)
               #---------------------------------------------------------------------------#


               #----- Recruitment rate in terms of above-ground biomass. ------------------#
               population             = sum(w.nplant[sel.pop] * agbconow[sel.pop])
               established            = sum(w.nplant[sel.est] * agbconow[sel.est])
               szpft$agb.recr [m,d,p] = log(population / established)
               #---------------------------------------------------------------------------#


               #----- Recruitment rate in terms of basal area. ----------------------------#
               population             = sum(w.nplant[sel.pop] * baconow [sel.pop])
               established            = sum(w.nplant[sel.est] * baconow [sel.est])
               szpft$bsa.recr [m,d,p] = log(population / established)
               #---------------------------------------------------------------------------#
            }#end if
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Census variables.  They have additional restrictions.                    #
            #------------------------------------------------------------------------------#
            sel.tall = heightconow >= census.height.min
            sel.fat  = dbhconow    >= census.dbh.min
            #----- "Census" LAI, WAI, and TAI discard small trees. ------------------------#
            sel = sel.pft & sel.dbh & sel.tall
            if (any(sel)){
               szpft$census.lai[m,d,p] = sum(laiconow[sel] * areaconow[sel])
               szpft$census.wai[m,d,p] = sum(waiconow[sel] * areaconow[sel])
               szpft$census.tai[m,d,p] = sum(taiconow[sel] * areaconow[sel])
            }#end if
            #----- "Census" AGB and BA discard skinny trees. ------------------------------#
            sel = sel.pft & sel.dbh & sel.fat
            if (any(sel)){
               szpft$census.agb[m,d,p] = sum(agbconow[sel] * w.nplant[sel])
               szpft$census.ba [m,d,p] = sum(baconow [sel] * w.nplant[sel])
            }#end if
            #------------------------------------------------------------------------------#
         }#end for PFT
         #---------------------------------------------------------------------------------#
      }#end for DBH
      #------------------------------------------------------------------------------------#




      #------------------------------------------------------------------------------------#
      #       Find the above-ground NPP and rain use efficiency.                           #
      #------------------------------------------------------------------------------------#
      szpft$last.1yr.anpp  [m,,] = apply( X      = szpft$anpp  [max(1,m-12):m,,,drop=FALSE]
                                        , MARGIN = c(2,3)
                                        , FUN    = mean
                                        )
      szpft$last.2yr.anpp  [m,,] = apply( X      = szpft$anpp  [max(1,m-24):m,,,drop=FALSE]
                                        , MARGIN = c(2,3)
                                        , FUN    = mean
                                        )
      szpft$last.3yr.anpp  [m,,] = apply( X      = szpft$anpp  [max(1,m-36):m,,,drop=FALSE]
                                        , MARGIN = c(2,3)
                                        , FUN    = mean
                                        )
      szpft$last.1yr.bnpp  [m,,] = apply( X      = szpft$bnpp  [max(1,m-12):m,,,drop=FALSE]
                                        , MARGIN = c(2,3)
                                        , FUN    = mean
                                        )
      szpft$last.2yr.bnpp  [m,,] = apply( X      = szpft$bnpp  [max(1,m-24):m,,,drop=FALSE]
                                        , MARGIN = c(2,3)
                                        , FUN    = mean
                                        )
      szpft$last.3yr.bnpp  [m,,] = apply( X      = szpft$bnpp  [max(1,m-36):m,,,drop=FALSE]
                                        , MARGIN = c(2,3)
                                        , FUN    = mean
                                        )
      szpft$last.1yr.transp[m,,] = apply( X      = szpft$transp[max(1,m-12):m,,,drop=FALSE]
                                      , MARGIN = c(2,3)
                                      , FUN    = mean
                                      )
      szpft$last.2yr.transp[m,,] = apply( X      = szpft$transp[max(1,m-24):m,,,drop=FALSE]
                                      , MARGIN = c(2,3)
                                      , FUN    = mean
                                      )
      szpft$last.3yr.transp[m,,] = apply( X      = szpft$transp[max(1,m-36):m,,,drop=FALSE]
                                      , MARGIN = c(2,3)
                                      , FUN    = mean
                                      )
      szpft$last.1yr.arue  [m,,] = 1000. * ( szpft$last.1yr.anpp[m,,] 
                                           / max(emean$last.1yr.rain[m],0.01) )
      szpft$last.2yr.arue  [m,,] = 1000. * ( szpft$last.2yr.anpp[m,,] 
                                           / max(emean$last.2yr.rain[m],0.01) )
      szpft$last.3yr.arue  [m,,] = 1000. * ( szpft$last.3yr.anpp[m,,]
                                           / max(emean$last.3yr.rain[m],0.01) )
      szpft$last.1yr.brue  [m,,] = 100.  * ( szpft$last.1yr.bnpp[m,,] 
                                           / max(emean$last.1yr.rain[m],0.01) )
      szpft$last.2yr.brue  [m,,] = 100.  * ( szpft$last.2yr.bnpp[m,,] 
                                           / max(emean$last.2yr.rain[m],0.01) )
      szpft$last.3yr.brue  [m,,] = 100.  * ( szpft$last.3yr.bnpp[m,,] 
                                           / max(emean$last.3yr.rain[m],0.01) )
      #------------------------------------------------------------------------------------#



      #----- Find the water use efficiency. -----------------------------------------------#
      szpft$last.1yr.wue[m,,] = ifelse( szpft$last.1yr.transp [m,,] > 1.e-10
                                      , 1000 * szpft$last.1yr.anpp   [m,,]
                                             / szpft$last.1yr.transp [m,,] / yr.day
                                      , 0
                                      )#end ifelse
      szpft$last.2yr.wue[m,,] = ifelse( szpft$last.2yr.transp [m,,] > 1.e-10
                                      , 1000 * szpft$last.2yr.anpp   [m,,]
                                             / szpft$last.2yr.transp [m,,] / yr.day
                                      , 0
                                      )#end ifelse
      szpft$last.3yr.wue[m,,] = ifelse( szpft$last.3yr.transp [m,,] > 1.e-10
                                      , 1000 * szpft$last.3yr.anpp   [m,,]
                                             / szpft$last.3yr.transp [m,,] / yr.day
                                      , 0
                                      )#end ifelse
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #       Build the derived variables.                                                 #
      #------------------------------------------------------------------------------------#
      emean$npp             [m] = szpft$npp            [m,ndbh+1,npft+1]
      emean$mco             [m] = szpft$mco            [m,ndbh+1,npft+1]
      emean$cba             [m] = szpft$cba            [m,ndbh+1,npft+1]
      emean$cbamax          [m] = szpft$cbamax         [m,ndbh+1,npft+1]
      emean$cbalight        [m] = szpft$cbalight       [m,ndbh+1,npft+1]
      emean$cbamoist        [m] = szpft$cbamoist       [m,ndbh+1,npft+1]
      emean$cbarel          [m] = szpft$cbarel         [m,ndbh+1,npft+1]
      emean$nplant          [m] = szpft$nplant         [m,ndbh+1,npft+1]
      emean$lai             [m] = szpft$lai            [m,ndbh+1,npft+1]
      emean$wai             [m] = szpft$wai            [m,ndbh+1,npft+1]
      emean$tai             [m] = szpft$tai            [m,ndbh+1,npft+1]
      emean$agb             [m] = szpft$agb            [m,ndbh+1,npft+1]
      emean$bgb             [m] = szpft$bgb            [m,ndbh+1,npft+1]
      emean$biomass         [m] = szpft$biomass        [m,ndbh+1,npft+1]
      emean$ldrop           [m] = szpft$ldrop          [m,ndbh+1,npft+1]
      emean$demand          [m] = szpft$demand         [m,ndbh+1,npft+1]
      emean$supply          [m] = szpft$supply         [m,ndbh+1,npft+1]
      emean$i.gpp           [m] = szpft$i.gpp          [m,ndbh+1,npft+1]
      emean$i.npp           [m] = szpft$i.npp          [m,ndbh+1,npft+1]
      emean$i.plresp        [m] = szpft$i.plresp       [m,ndbh+1,npft+1]
      emean$i.mco           [m] = szpft$i.mco          [m,ndbh+1,npft+1]
      emean$i.cba           [m] = szpft$i.cba          [m,ndbh+1,npft+1]
      emean$i.cbamax        [m] = szpft$i.cbamax       [m,ndbh+1,npft+1]
      emean$i.cbalight      [m] = szpft$i.cbalight     [m,ndbh+1,npft+1]
      emean$i.cbamoist      [m] = szpft$i.cbamoist     [m,ndbh+1,npft+1]
      emean$i.transp        [m] = szpft$i.transp       [m,ndbh+1,npft+1]
      emean$i.wflxlc        [m] = szpft$i.wflxlc       [m,ndbh+1,npft+1]
      emean$i.hflxlc        [m] = szpft$i.hflxlc       [m,ndbh+1,npft+1]
      emean$f.gpp           [m] = szpft$f.gpp          [m,ndbh+1,npft+1]
      emean$f.plant.resp    [m] = szpft$f.plant.resp   [m,ndbh+1,npft+1]
      emean$f.npp           [m] = szpft$f.npp          [m,ndbh+1,npft+1]
      emean$f.cba           [m] = szpft$f.cba          [m,ndbh+1,npft+1]
      emean$f.bstorage      [m] = szpft$f.bstorage     [m,ndbh+1,npft+1]
      emean$f.bleaf         [m] = szpft$f.bleaf        [m,ndbh+1,npft+1]
      emean$f.broot         [m] = szpft$f.broot        [m,ndbh+1,npft+1]
      emean$f.bseeds        [m] = szpft$f.bseeds       [m,ndbh+1,npft+1]
      emean$leaf.par        [m] = szpft$leaf.par       [m,ndbh+1,npft+1]
      emean$leaf.rshort     [m] = szpft$leaf.rshort    [m,ndbh+1,npft+1]
      emean$leaf.rlong      [m] = szpft$leaf.rlong     [m,ndbh+1,npft+1]
      emean$transp          [m] = szpft$transp         [m,ndbh+1,npft+1]
      emean$last.1yr.transp [m] = szpft$last.1yr.transp[m,ndbh+1,npft+1]
      emean$last.2yr.transp [m] = szpft$last.2yr.transp[m,ndbh+1,npft+1]
      emean$last.3yr.transp [m] = szpft$last.3yr.transp[m,ndbh+1,npft+1]
      emean$wue             [m] = szpft$wue            [m,ndbh+1,npft+1]
      emean$last.1yr.wue    [m] = szpft$last.1yr.wue   [m,ndbh+1,npft+1]
      emean$last.2yr.wue    [m] = szpft$last.2yr.wue   [m,ndbh+1,npft+1]
      emean$last.3yr.wue    [m] = szpft$last.3yr.wue   [m,ndbh+1,npft+1]
      emean$anpp            [m] = szpft$anpp           [m,ndbh+1,npft+1]
      emean$last.1yr.anpp   [m] = szpft$last.1yr.anpp  [m,ndbh+1,npft+1]
      emean$last.2yr.anpp   [m] = szpft$last.2yr.anpp  [m,ndbh+1,npft+1]
      emean$last.3yr.anpp   [m] = szpft$last.3yr.anpp  [m,ndbh+1,npft+1]
      emean$bnpp            [m] = szpft$bnpp           [m,ndbh+1,npft+1]
      emean$last.1yr.bnpp   [m] = szpft$last.1yr.bnpp  [m,ndbh+1,npft+1]
      emean$last.2yr.bnpp   [m] = szpft$last.2yr.bnpp  [m,ndbh+1,npft+1]
      emean$last.3yr.bnpp   [m] = szpft$last.3yr.bnpp  [m,ndbh+1,npft+1]
      emean$arue            [m] = szpft$arue           [m,ndbh+1,npft+1]
      emean$last.1yr.arue   [m] = szpft$last.1yr.arue  [m,ndbh+1,npft+1]
      emean$last.2yr.arue   [m] = szpft$last.2yr.arue  [m,ndbh+1,npft+1]
      emean$last.3yr.arue   [m] = szpft$last.3yr.arue  [m,ndbh+1,npft+1]
      emean$brue            [m] = szpft$brue           [m,ndbh+1,npft+1]
      emean$last.1yr.brue   [m] = szpft$last.1yr.brue  [m,ndbh+1,npft+1]
      emean$last.2yr.brue   [m] = szpft$last.2yr.brue  [m,ndbh+1,npft+1]
      emean$last.3yr.brue   [m] = szpft$last.3yr.brue  [m,ndbh+1,npft+1]
      #------------------------------------------------------------------------------------#




      #------------------------------------------------------------------------------------#
      #     Define evapotranspiration use efficiency.                                      #
      #------------------------------------------------------------------------------------#
      emean$etue         [m] =   emean$anpp[m] / max(emean$et[m]*yr.day,0.01)
      emean$last.1yr.etue[m] = ( 1000. * emean$last.1yr.anpp[m] 
                               / max(emean$last.1yr.et[m]*yr.day,0.01) )
      emean$last.2yr.etue[m] = ( 1000. * emean$last.2yr.anpp[m] 
                               / max(emean$last.2yr.et[m]*yr.day,0.01) )
      emean$last.3yr.etue[m] = ( 1000. * emean$last.3yr.anpp[m] 
                               / max(emean$last.3yr.et[m]*yr.day,0.01) )
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #      Build the cohort-level lists if this is the right month.                      #
      #------------------------------------------------------------------------------------#
      if (thismonth %in% sasmonth){
         clab = paste( "y",sprintf("%4.4i",thisyear )
                     , "m",sprintf("%2.2i",thismonth),sep="")
         #----- Binding the current cohorts. ----------------------------------------------#
         cohort$ipa          [[clab]] = ipaconow
         cohort$ico          [[clab]] = icoconow
         cohort$area         [[clab]] = areaconow
         cohort$lu           [[clab]] = luconow
         cohort$dbh          [[clab]] = dbhconow
         cohort$age          [[clab]] = ageconow
         cohort$pft          [[clab]] = pftconow
         cohort$nplant       [[clab]] = nplantconow * areaconow
         cohort$height       [[clab]] = heightconow
         cohort$ba           [[clab]] = nplantconow * baconow * areaconow
         cohort$agb          [[clab]] = agbconow
         cohort$bgb          [[clab]] = bgbconow
         cohort$biomass      [[clab]] = biomassconow
         cohort$lai          [[clab]] = laiconow
         cohort$wai          [[clab]] = waiconow
         cohort$tai          [[clab]] = taiconow
         cohort$gpp          [[clab]] = gppconow
         cohort$leaf.resp    [[clab]] = leaf.respconow
         cohort$root.resp    [[clab]] = root.respconow
         cohort$plant.resp   [[clab]] = plant.respconow
         cohort$npp          [[clab]] = nppconow
         cohort$cba          [[clab]] = cbaconow
         cohort$cbamax       [[clab]] = cbamaxconow
         cohort$cbalight     [[clab]] = cbalightconow
         cohort$cbamoist     [[clab]] = cbamoistconow
         cohort$cbarel       [[clab]] = cbarelconow
         cohort$mcost        [[clab]] = mcostconow
         cohort$ldrop        [[clab]] = ldropconow
         cohort$fs.open      [[clab]] = fs.openconow
         cohort$light        [[clab]] = lightconow
         cohort$light.beam   [[clab]] = light.beamconow
         cohort$light.diff   [[clab]] = light.diffconow
         cohort$balive       [[clab]] = baliveconow
         cohort$bdead        [[clab]] = bdeadconow
         cohort$bleaf        [[clab]] = bleafconow
         cohort$broot        [[clab]] = brootconow
         cohort$bsapwood     [[clab]] = bsapwoodconow
         cohort$bstorage     [[clab]] = bstorageconow
         cohort$bseeds       [[clab]] = bseedsconow
         cohort$hflxlc       [[clab]] = hflxlcconow
         cohort$wflxlc       [[clab]] = wflxlcconow
         cohort$transp       [[clab]] = transpconow
         cohort$wue          [[clab]] = wueconow
         cohort$demand       [[clab]] = demandconow
         cohort$supply       [[clab]] = supplyconow
         cohort$mort         [[clab]] = 100. * (1.0 - exp(-mortconow      ))
         cohort$ncbmort      [[clab]] = 100. * (1.0 - exp(-ncbmortconow   ))
         cohort$dimort       [[clab]] = 100. * (1.0 - exp(-dimortconow    ))
         cohort$recruit      [[clab]] = 100. * (exp(recruitconow) - 1.0)
         cohort$growth       [[clab]] = growthconow
         cohort$agb.growth   [[clab]] = agb.growthconow
         cohort$bsa.growth   [[clab]] = bsa.growthconow
         cohort$f.gpp        [[clab]] = f.gppconow
         cohort$f.plant.resp [[clab]] = f.plant.respconow
         cohort$f.npp        [[clab]] = f.nppconow
         cohort$f.cba        [[clab]] = f.cbaconow
         cohort$f.bstorage   [[clab]] = f.bstorageconow
         cohort$f.bleaf      [[clab]] = f.bleafconow
         cohort$f.broot      [[clab]] = f.brootconow
         cohort$f.bseeds     [[clab]] = f.bseedsconow
         cohort$leaf.par     [[clab]] = leaf.parconow
         cohort$leaf.rshort  [[clab]] = leaf.rshortconow
         cohort$leaf.rlong   [[clab]] = leaf.rlongconow
         cohort$anpp         [[clab]] = anppconow
         cohort$bnpp         [[clab]] = bnppconow
         cohort$arue         [[clab]] = arueconow
         cohort$brue         [[clab]] = brueconow
      } #end if month=sasmonth
      #------------------------------------------------------------------------------------#
   }# end for (m in tresume,ntimes)
   #---------------------------------------------------------------------------------------#







   #---------------------------------------------------------------------------------------#
   #     Copy the variables back to datum.                                                 #
   #---------------------------------------------------------------------------------------#
   datum$emean  = emean
   datum$emsqu  = emsqu
   datum$szpft  = szpft
   datum$lu     = lu
   datum$qmean  = qmean
   datum$qmsqu  = qmsqu
   datum$patch  = patch
   datum$cohort = cohort
   #---------------------------------------------------------------------------------------#

   return(datum)
   #---------------------------------------------------------------------------------------#
}#end function read.q.files
#==========================================================================================#
#==========================================================================================#
