 #==========================================================================================#
#==========================================================================================#
#     This function reads the ED2 monthly mean files that contain mean diurnal cycle.      #
#   Inputs:                                                                                #
#   - datum   -- The monthly structure that will contain the data.  It must be initialised #
#                by create.monthly, otherwise it won't work.                               #
#   - ntimes  -- Total number of times (including previously loaded times).                #
#   - tresume -- The first time to read (in case data have been partially loaded.          #
#------------------------------------------------------------------------------------------#
read.q.simple <<- function(datum,ntimes,tresume=1,sasmonth=5){


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
   qmean  = datum$qmean
   qmsqu  = datum$qmsqu
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
      #     Set daytime flag.                                                              #
      #------------------------------------------------------------------------------------#
      phap        = mymont$QMEAN.ATM.RSHORT.PY > phap.min | iint.photo == 0
      polar.night = ! any(phap,na.rm=TRUE)
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
      emean$leaf.water      [m] =   mymont$MMEAN.LEAF.WATER.PY
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
      emean$et              [m] = emean$evap[m] + emean$transp[m]
      emean$rain            [m] = mymont$MMEAN.PCPG.PY * mondays * day.sec

      emean$sm.stress       [m] =   1. - mymont$MMEAN.FS.OPEN.PY
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



      #----- Read cohort-derived properties. ----------------------------------------------#
      emean$agb       [m] = sum(mymont$AGB.PY       )
      emean$lai       [m] = sum(mymont$MMEAN.LAI.PY )
      emean$wai       [m] = sum(mymont$WAI.PY       )
      emean$basarea   [m] = sum(mymont$BASAL.AREA.PY)
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
      qmean$leaf.water  [m,] =   mymont$QMEAN.LEAF.WATER.PY
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
      qmean$sm.stress   [m,] =   1. - mymont$QMEAN.FS.OPEN.PY
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



      #------------------------------------------------------------------------------------#
      #      Convert leaf water to kg/m2leaf.                                              #
      #------------------------------------------------------------------------------------#
      emean$leaf.water  [m ] = emean$leaf.water[m ] / pmax(emean$lai[m],0.01)
      qmean$leaf.water  [m,] = qmean$leaf.water[m,] / pmax(emean$lai[m],0.01)
      #------------------------------------------------------------------------------------#




      #----- Find drought length using rainfall running average. --------------------------#
      if ( m == 1){
         emean$nmon.lt.090[m] = as.numeric(emean$rain[m] <  90  )
         emean$nmon.lt.100[m] = as.numeric(emean$rain[m] < 100  )
         emean$nmon.lt.110[m] = as.numeric(emean$rain[m] < 110  )
         emean$nmon.lt.120[m] = as.numeric(emean$rain[m] < 120  )
         emean$nmon.wdef  [m] = as.numeric(emean$water.deficit[m] >  10  )
         emean$nmon.mdef  [m] = as.numeric(emean$malhi.deficit[m] >  10  )
      }else{
         ma                   = max(1,m-11)
         ra.rain              = mean(emean$rain[ma:m])
         wdef                 = emean$water.deficit[m]
         mdef                 = emean$malhi.deficit[m]
         emean$nmon.lt.090[m] = as.numeric(ra.rain <  90) * ( emean$nmon.lt.090[m-1] + 1 )
         emean$nmon.lt.100[m] = as.numeric(ra.rain < 100) * ( emean$nmon.lt.100[m-1] + 1 )
         emean$nmon.lt.110[m] = as.numeric(ra.rain < 110) * ( emean$nmon.lt.110[m-1] + 1 )
         emean$nmon.lt.120[m] = as.numeric(ra.rain < 120) * ( emean$nmon.lt.120[m-1] + 1 )
         emean$nmon.wdef  [m] = as.numeric(wdef    >  10) * ( emean$nmon.wdef  [m-1] + 1 )
         emean$nmon.mdef  [m] = as.numeric(mdef    >  10) * ( emean$nmon.mdef  [m-1] + 1 )
      }#end if
      #------------------------------------------------------------------------------------#
   }# end for (m in tresume,ntimes)
   #---------------------------------------------------------------------------------------#







   #---------------------------------------------------------------------------------------#
   #     Copy the variables back to datum.                                                 #
   #---------------------------------------------------------------------------------------#
   datum$emean  = emean
   datum$emsqu  = emsqu
   datum$qmean  = qmean
   datum$qmsqu  = qmsqu
   #---------------------------------------------------------------------------------------#

   return(datum)
   #---------------------------------------------------------------------------------------#
}#end function read.q.files
#==========================================================================================#
#==========================================================================================#
