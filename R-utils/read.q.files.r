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
   qpatch = datum$qpatch
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
      #----- Make a growth respiration variable. ------------------------------------------#
      if (! "MMEAN.GROWTH.RESP.PY" %in% names(mymont)){
         mymont$MMEAN.GROWTH.RESP.PY  = ( mymont$MMEAN.LEAF.GROWTH.RESP.PY
                                        + mymont$MMEAN.ROOT.GROWTH.RESP.PY
                                        + mymont$MMEAN.SAPA.GROWTH.RESP.PY
                                        + mymont$MMEAN.SAPB.GROWTH.RESP.PY )
         mymont$QMEAN.GROWTH.RESP.PY  = ( mymont$QMEAN.LEAF.GROWTH.RESP.PY
                                        + mymont$QMEAN.ROOT.GROWTH.RESP.PY
                                        + mymont$QMEAN.SAPA.GROWTH.RESP.PY
                                        + mymont$QMEAN.SAPB.GROWTH.RESP.PY )
         #----- Apply the same correction term to cohort-level variables. -----------------#
         if (mymont$NCOHORTS.GLOBAL > 0){
         mymont$MMEAN.GROWTH.RESP.CO  = ( mymont$MMEAN.LEAF.GROWTH.RESP.CO
                                        + mymont$MMEAN.ROOT.GROWTH.RESP.CO
                                        + mymont$MMEAN.SAPA.GROWTH.RESP.CO
                                        + mymont$MMEAN.SAPB.GROWTH.RESP.CO )
         mymont$QMEAN.GROWTH.RESP.CO  = ( mymont$QMEAN.LEAF.GROWTH.RESP.CO
                                        + mymont$QMEAN.ROOT.GROWTH.RESP.CO
                                        + mymont$QMEAN.SAPA.GROWTH.RESP.CO
                                        + mymont$QMEAN.SAPB.GROWTH.RESP.CO )
         }#end if
         #---------------------------------------------------------------------------------#
      }#end if (! "MMEAN.GROWTH.RESP.PY" %in% names(mymont))
      #----- Make a growth respiration variable. ------------------------------------------#
      if (! "MMEAN.STORAGE.RESP.PY" %in% names(mymont)){
         mymont$MMEAN.STORAGE.RESP.PY  = ( mymont$MMEAN.LEAF.STORAGE.RESP.PY
                                         + mymont$MMEAN.ROOT.STORAGE.RESP.PY
                                         + mymont$MMEAN.SAPA.STORAGE.RESP.PY
                                         + mymont$MMEAN.SAPB.STORAGE.RESP.PY )
         mymont$QMEAN.STORAGE.RESP.PY  = ( mymont$QMEAN.LEAF.STORAGE.RESP.PY
                                         + mymont$QMEAN.ROOT.STORAGE.RESP.PY
                                         + mymont$QMEAN.SAPA.STORAGE.RESP.PY
                                         + mymont$QMEAN.SAPB.STORAGE.RESP.PY )
         #----- Apply the same correction term to cohort-level variables. -----------------#
         if (mymont$NCOHORTS.GLOBAL > 0){
         mymont$MMEAN.STORAGE.RESP.CO  = ( mymont$MMEAN.LEAF.STORAGE.RESP.CO
                                         + mymont$MMEAN.ROOT.STORAGE.RESP.CO
                                         + mymont$MMEAN.SAPA.STORAGE.RESP.CO
                                         + mymont$MMEAN.SAPB.STORAGE.RESP.CO )
         mymont$QMEAN.STORAGE.RESP.CO  = ( mymont$QMEAN.LEAF.STORAGE.RESP.CO
                                         + mymont$QMEAN.ROOT.STORAGE.RESP.CO
                                         + mymont$QMEAN.SAPA.STORAGE.RESP.CO
                                         + mymont$QMEAN.SAPB.STORAGE.RESP.CO )
         }#end if
         #---------------------------------------------------------------------------------#
      }#end if (! "MMEAN.STORAGE.RESP.PY" %in% names(mymont))
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #      Fix units for some variables.                                                 #
      #------------------------------------------------------------------------------------#
      if (corr.growth.storage != 1.){
         cat("     - Correcting growth and storage respiration...","\n")
         #----- Correct the rates. --------------------------------------------------------#
         multfac = corr.growth.storage
         mymont$MMEAN.GROWTH.RESP.PY  =   mymont$MMEAN.GROWTH.RESP.PY  * multfac
         mymont$MMEAN.STORAGE.RESP.PY =   mymont$MMEAN.STORAGE.RESP.PY * multfac
         mymont$QMEAN.GROWTH.RESP.PY  =   mymont$QMEAN.GROWTH.RESP.PY  * multfac
         mymont$QMEAN.STORAGE.RESP.PY =   mymont$QMEAN.STORAGE.RESP.PY * multfac
         if (mymont$NCOHORTS.GLOBAL > 0){
            mymont$MMEAN.GROWTH.RESP.CO  = mymont$MMEAN.GROWTH.RESP.CO  * multfac
            mymont$MMEAN.STORAGE.RESP.CO = mymont$MMEAN.STORAGE.RESP.CO * multfac
            mymont$QMEAN.GROWTH.RESP.CO  = mymont$QMEAN.GROWTH.RESP.CO  * multfac
            mymont$QMEAN.STORAGE.RESP.CO = mymont$QMEAN.STORAGE.RESP.CO * multfac
         }#end if
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #      Correct Plant respiration.                                                 #
         #---------------------------------------------------------------------------------#
         mymont$MMEAN.PLRESP.PY       = ( mymont$MMEAN.LEAF.RESP.PY
                                        + mymont$MMEAN.ROOT.RESP.PY
                                        + mymont$MMEAN.GROWTH.RESP.PY
                                        + mymont$MMEAN.STORAGE.RESP.PY
                                        )#end 
         mymont$QMEAN.PLRESP.PY       = ( mymont$QMEAN.LEAF.RESP.PY
                                        + mymont$QMEAN.ROOT.RESP.PY
                                        + mymont$QMEAN.GROWTH.RESP.PY
                                        + mymont$QMEAN.STORAGE.RESP.PY
                                        )#end 
         mymont$MMSQU.PLRESP.PY       = mymont$MMSQU.PLRESP.PY + NA
         mymont$QMSQU.PLRESP.PY       = mymont$QMSQU.PLRESP.PY + NA
         if (mymont$NCOHORTS.GLOBAL > 0){
            mymont$MMEAN.PLRESP.CO       = ( mymont$MMEAN.LEAF.RESP.CO
                                           + mymont$MMEAN.ROOT.RESP.CO
                                           + mymont$MMEAN.GROWTH.RESP.CO
                                           + mymont$MMEAN.STORAGE.RESP.CO
                                           )#end 
            mymont$QMEAN.PLRESP.CO       = ( mymont$QMEAN.LEAF.RESP.CO
                                           + mymont$QMEAN.ROOT.RESP.CO
                                           + mymont$QMEAN.GROWTH.RESP.CO
                                           + mymont$QMEAN.STORAGE.RESP.CO
                                           )#end 
            mymont$MMSQU.PLRESP.CO       = mymont$MMSQU.PLRESP.CO + NA
            mymont$QMSQU.PLRESP.CO       = mymont$QMSQU.PLRESP.CO + NA
         }#end if
         #---------------------------------------------------------------------------------#




         #---------------------------------------------------------------------------------#
         #      Correct NPP.                                                               #
         #---------------------------------------------------------------------------------#
         mymont$MMEAN.NPP.PY    = mymont$MMEAN.GPP.PY - mymont$MMEAN.PLRESP.PY
         mymont$QMEAN.NPP.PY    = mymont$QMEAN.GPP.PY - mymont$QMEAN.PLRESP.PY
         mymont$MMSQU.NPP.PY    = mymont$MMSQU.NPP.PY + NA
         mymont$QMSQU.NPP.PY    = mymont$QMSQU.NPP.PY + NA
         if (mymont$NCOHORTS.GLOBAL > 0){
            mymont$MMEAN.NPP.CO = mymont$MMEAN.GPP.CO - mymont$MMEAN.PLRESP.CO
            mymont$QMEAN.NPP.CO = mymont$QMEAN.GPP.CO - mymont$QMEAN.PLRESP.CO
            mymont$MMSQU.NPP.CO = mymont$MMSQU.NPP.CO + NA
            mymont$QMSQU.NPP.CO = mymont$QMSQU.NPP.CO + NA
         }#end if
         #---------------------------------------------------------------------------------#




         #---------------------------------------------------------------------------------#
         #      Correct NEP.                                                               #
         #---------------------------------------------------------------------------------#
         mymont$MMEAN.NEP.PY    = mymont$MMEAN.NPP.PY - mymont$MMEAN.RH.PY
         mymont$QMEAN.NEP.PY    = mymont$QMEAN.NPP.PY - mymont$QMEAN.RH.PY
         mymont$MMSQU.NEP.PY    = mymont$MMSQU.NEP.PY + NA
         mymont$QMSQU.NEP.PY    = mymont$QMSQU.NEP.PY + NA

         mymont$MMEAN.NEP.PA    = -mymont$MMEAN.RH.PA
         mymont$QMEAN.NEP.PA    = -mymont$QMEAN.RH.PA
         mymont$MMSQU.NEP.PA    = mymont$MMSQU.NEP.PA + NA
         mymont$QMSQU.NEP.PA    = mymont$QMSQU.NEP.PA + NA
         if (mymont$NCOHORTS.GLOBAL > 0){
            ipaconow    = rep(sequence(mymont$NPATCHES.GLOBAL),times=mymont$PACO.N)
            idx         = match(unique(ipaconow),sequence(mymont$NPATCHES.GLOBAL))
            mymont$MMEAN.NEP.PA[idx]  = ( tapply( X     = mymont$MMEAN.NPP.CO
                                                        * mymont$NPLANT
                                                , INDEX = ipaconow
                                                , FUN   = sum
                                                )#end tapply
                                        - mymont$MMEAN.RH.PA[idx]
                                        )#end
            mymont$QMEAN.NEP.PA[idx,] = ( qapply( X     = mymont$QMEAN.NPP.CO
                                                        * mymont$NPLANT
                                                , INDEX = ipaconow
                                                , DIM   = 1
                                                , FUN   = sum
                                                )#end tapply
                                        - mymont$QMEAN.RH.PA[idx,]
                                        )#end
         }#end if
         #---------------------------------------------------------------------------------#
      }#end if
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
      emean$nep             [m] =   mymont$MMEAN.NEP.PY
      emean$nee             [m] =   mymont$MMEAN.CARBON.ST.PY - mymont$MMEAN.CARBON.AC.PY
      emean$plant.resp      [m] =   mymont$MMEAN.PLRESP.PY
      emean$growth.resp     [m] =   mymont$MMEAN.GROWTH.RESP.PY
      emean$storage.resp    [m] =   mymont$MMEAN.STORAGE.RESP.PY
      emean$assim.light     [m] =   mymont$MMEAN.A.LIGHT.PY
      emean$assim.rubp      [m] =   mymont$MMEAN.A.RUBP.PY
      emean$assim.co2       [m] =   mymont$MMEAN.A.CO2.PY
      emean$assim.ratio     [m] = ( mymont$MMEAN.A.LIGHT.PY
                                  / max( 1e-6,min( mymont$MMEAN.A.RUBP.PY
                                                 , mymont$MMEAN.A.CO2.PY  )))
      #----- (Leaf, root, stem, and soil respiration are corrected for growth+storage). ---#
      emean$reco            [m] =   mymont$MMEAN.PLRESP.PY    + mymont$MMEAN.RH.PY
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
      emean$par.leaf        [m] =   mymont$MMEAN.PAR.L.PY          * Watts.2.Ein * 1.e6
      emean$par.leaf.beam   [m] =   mymont$MMEAN.PAR.L.BEAM.PY     * Watts.2.Ein * 1.e6
      emean$par.leaf.diff   [m] =   mymont$MMEAN.PAR.L.DIFF.PY     * Watts.2.Ein * 1.e6
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
      emean$soil.temp     [m,] =   mymont$MMEAN.SOIL.TEMP.PY - t00
      emean$soil.water    [m,] =   mymont$MMEAN.SOIL.WATER.PY
      emean$soil.mstpot   [m,] = - mymont$MMEAN.SOIL.MSTPOT.PY * grav * wdnsi
      emean$soil.extracted[m,] = - mymont$MMEAN.TRANSLOSS.PY * day.sec / dslz
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
      qmean$gpp          [m,] =   mymont$QMEAN.GPP.PY
      qmean$npp          [m,] =   mymont$QMEAN.NPP.PY
      qmean$het.resp     [m,] =   mymont$QMEAN.RH.PY
      qmean$cwd.resp     [m,] =   mymont$QMEAN.CWD.RH.PY
      qmean$assim.light  [m,] =   mymont$QMEAN.A.LIGHT.PY
      qmean$assim.rubp   [m,] =   mymont$QMEAN.A.RUBP.PY
      qmean$assim.co2    [m,] =   mymont$QMEAN.A.CO2.PY
      qmean$assim.ratio  [m,] = ( mymont$QMEAN.A.LIGHT.PY
                                / pmax(1e-6, pmin( mymont$QMEAN.A.RUBP.PY
                                                 , mymont$QMEAN.A.CO2.PY  )))
      qmean$nee          [m,] = ( mymont$QMEAN.CARBON.ST.PY
                                - mymont$QMEAN.CARBON.AC.PY )
      qmean$reco         [m,] =   mymont$QMEAN.PLRESP.PY + mymont$QMEAN.RH.PY
      qmean$cflxca       [m,] = - mymont$QMEAN.CARBON.AC.PY
      qmean$cflxst       [m,] = - mymont$QMEAN.CARBON.ST.PY
      qmean$hflxca       [m,] = - mymont$QMEAN.SENSIBLE.AC.PY
      qmean$hflxlc       [m,] =   mymont$QMEAN.SENSIBLE.LC.PY
      qmean$hflxwc       [m,] =   mymont$QMEAN.SENSIBLE.WC.PY
      qmean$hflxgc       [m,] =   mymont$QMEAN.SENSIBLE.GC.PY
      qmean$wflxca       [m,] = - mymont$QMEAN.VAPOR.AC.PY     * day.sec
      qmean$qwflxca      [m,] = - mymont$QMEAN.VAPOR.AC.PY     * qmean.can.alvli.py
      qmean$wflxlc       [m,] =   mymont$QMEAN.VAPOR.LC.PY     * day.sec
      qmean$wflxwc       [m,] =   mymont$QMEAN.VAPOR.WC.PY     * day.sec
      qmean$wflxgc       [m,] =   mymont$QMEAN.VAPOR.GC.PY     * day.sec
      qmean$runoff       [m,] = ( mymont$QMEAN.RUNOFF.PY
                                + mymont$QMEAN.DRAINAGE.PY )   * day.sec
      qmean$intercepted  [m,] = ( mymont$QMEAN.INTERCEPTED.AL.PY
                                + mymont$QMEAN.INTERCEPTED.AW.PY ) * day.sec
      qmean$wshed        [m,] = ( mymont$QMEAN.WSHED.LG.PY
                                + mymont$QMEAN.WSHED.WG.PY       ) * day.sec
      qmean$evap         [m,] = ( mymont$QMEAN.VAPOR.GC.PY
                                + mymont$QMEAN.VAPOR.WC.PY
                                + mymont$QMEAN.VAPOR.LC.PY )   * day.sec
      qmean$transp       [m,] =   mymont$QMEAN.TRANSP.PY       * day.sec
      qmean$atm.temp     [m,] =   mymont$QMEAN.ATM.TEMP.PY     - t00
      qmean$can.temp     [m,] =   mymont$QMEAN.CAN.TEMP.PY     - t00
      qmean$leaf.temp    [m,] =   mymont$QMEAN.LEAF.TEMP.PY    - t00
      qmean$leaf.water   [m,] =   mymont$QMEAN.LEAF.WATER.PY
      qmean$wood.temp    [m,] =   mymont$QMEAN.WOOD.TEMP.PY    - t00
      qmean$gnd.temp     [m,] =   mymont$QMEAN.GND.TEMP.PY     - t00
      qmean$atm.shv      [m,] =   mymont$QMEAN.ATM.SHV.PY      * kg2g
      qmean$can.shv      [m,] =   mymont$QMEAN.CAN.SHV.PY      * kg2g
      qmean$gnd.shv      [m,] =   mymont$QMEAN.GND.SHV.PY      * kg2g
      qmean$atm.vpd      [m,] =   mymont$QMEAN.ATM.VPDEF.PY    * 0.01
      qmean$can.vpd      [m,] =   mymont$QMEAN.CAN.VPDEF.PY    * 0.01
      qmean$leaf.vpd     [m,] =   mymont$QMEAN.LEAF.VPDEF.PY   * 0.01
      qmean$atm.co2      [m,] =   mymont$QMEAN.ATM.CO2.PY
      qmean$can.co2      [m,] =   mymont$QMEAN.CAN.CO2.PY
      qmean$atm.vels     [m,] =   mymont$QMEAN.ATM.VELS.PY
      qmean$ustar        [m,] =   mymont$QMEAN.USTAR.PY
      qmean$atm.prss     [m,] =   mymont$QMEAN.ATM.PRSS.PY     * 0.01
      qmean$can.prss     [m,] =   mymont$QMEAN.CAN.PRSS.PY     * 0.01
      qmean$sm.stress    [m,] =   1. - mymont$QMEAN.FS.OPEN.PY
      qmean$rain         [m,] =   mymont$QMEAN.PCPG.PY         * day.sec
      qmean$rshort       [m,] =   mymont$QMEAN.ATM.RSHORT.PY
      qmean$rshort.beam  [m,] = ( mymont$QMEAN.ATM.RSHORT.PY 
                                - mymont$QMEAN.ATM.RSHORT.DIFF.PY )
      qmean$rshort.diff  [m,] =   mymont$QMEAN.ATM.RSHORT.DIFF.PY
      qmean$rshort.gnd   [m,] =   mymont$QMEAN.RSHORT.GND.PY
      qmean$rshortup     [m,] =   mymont$QMEAN.RSHORTUP.PY
      qmean$rlong        [m,] =   mymont$QMEAN.ATM.RLONG.PY
      qmean$rlong.gnd    [m,] =   mymont$QMEAN.RLONG.GND.PY
      qmean$rlongup      [m,] =   mymont$QMEAN.RLONGUP.PY
      qmean$par.tot      [m,] =   mymont$QMEAN.ATM.PAR.PY        * Watts.2.Ein * 1.e6
      qmean$par.beam     [m,] = ( mymont$QMEAN.ATM.PAR.PY
                                - mymont$QMEAN.ATM.PAR.DIFF.PY ) * Watts.2.Ein * 1.e6
      qmean$par.diff     [m,] =   mymont$QMEAN.ATM.PAR.DIFF.PY   * Watts.2.Ein * 1.e6
      qmean$par.gnd      [m,] =   mymont$QMEAN.PAR.GND.PY        * Watts.2.Ein * 1.e6
      qmean$parup        [m,] =   mymont$QMEAN.PARUP.PY          * Watts.2.Ein * 1.e6
      qmean$par.leaf     [m,] =   mymont$QMEAN.PAR.L.PY          * Watts.2.Ein * 1.e6
      qmean$par.leaf.beam[m,] =   mymont$QMEAN.PAR.L.BEAM.PY     * Watts.2.Ein * 1.e6
      qmean$par.leaf.diff[m,] =   mymont$QMEAN.PAR.L.DIFF.PY     * Watts.2.Ein * 1.e6
      qmean$rnet         [m,] =   mymont$QMEAN.RNET.PY
      qmean$albedo       [m,] =   mymont$QMEAN.ALBEDO.PY
      if (all(c("QMEAN.ALBEDO.PAR.PY","QMEAN.ALBEDO.NIR.PY") %in% names(mymont))){
         qmean$albedo.par[m,] =   mymont$QMEAN.ALBEDO.PAR.PY
         qmean$albedo.nir[m,] =   mymont$QMEAN.ALBEDO.NIR.PY
      }else{
         qmean$albedo.par[m,] = ifelse( mymont$QMEAN.ATM.PAR.PY > 0.5
                                      , mymont$QMEAN.PARUP.PY / mymont$QMEAN.ATM.PAR.PY
                                      , mymont$QMEAN.ALBEDO.PY
                                      )#end ifelse
         qmean$albedo.nir[m,] = ifelse( mymont$QMEAN.ATM.NIR.PY > 0.5
                                      , mymont$QMEAN.NIRUP.PY / mymont$QMEAN.ATM.NIR.PY
                                      , mymont$QMEAN.ALBEDO.PY
                                      )#end ifelse
      }#end if
      qmean$rlong.albedo [m,] =   mymont$QMEAN.RLONG.ALBEDO.PY
      qmean$leaf.gbw     [m,] =   mymont$QMEAN.LEAF.GBW.PY       * day.sec
      qmean$leaf.gsw     [m,] =   mymont$QMEAN.LEAF.GSW.PY       * day.sec
      qmean$wood.gbw     [m,] =   mymont$QMEAN.WOOD.GBW.PY       * day.sec
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
      q.ipaconow  = matrix( data  = ipaconow
                          , nrow  = mymont$NCOHORTS.GLOBAL
                          , ncol  = mymont$NDCYC
                          , byrow = FALSE
                          )#end matrix
      icoconow    = unlist(sapply(X = mymont$PACO.N, FUN = sequence))
      idx         = match(unique(ipaconow),sequence(mymont$NPATCHES.GLOBAL))
      #------------------------------------------------------------------------------------#



      #----- Disturbance rates. -----------------------------------------------------------#
      lu$dist  [m,,] = apply ( X      = mymont$DISTURBANCE.RATES
                             , MARGIN = c(2,3)
                             , FUN    = weighted.mean
                             , w      = areasi
                             )#end apply
#                             browser()
      #------------------------------------------------------------------------------------#


      #====================================================================================#
      #====================================================================================#
      #====================================================================================#
      #====================================================================================#
      #    Cohort-level variables.                                                         #
      #------------------------------------------------------------------------------------#

      #------------------------------------------------------------------------------------#
      #     Read the cohort-level variables.  Because empty patchs do exist (deserts),     #
      # we must check whether there is any cohort to be read.  If not, assign NA to        #
      # all variables.                                                                     #
      #------------------------------------------------------------------------------------#
      one.cohort = sum(ncohorts) == 1
      one.patch  = sum(npatches) == 1
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
         dbhconow.lmon = mymont$DBH * exp(-pmax(0,mymont$DLNDBH.DT/12))
         dbhconow.1ago = mymont$DBH * exp(-pmax(0,mymont$DLNDBH.DT))
         dbhcut.1ago   = cut   (dbhconow.1ago,breaks=breakdbh)
         dbhlevs.1ago  = levels(dbhcut.1ago)
         dbhfac.1ago   = match (dbhcut.1ago,dbhlevs.1ago)
         dbhcut.lmon   = cut   (dbhconow.lmon,breaks=breakdbh)
         dbhlevs.lmon  = levels(dbhcut.lmon)
         dbhfac.lmon   = match (dbhcut.lmon,dbhlevs.lmon)
         #---------------------------------------------------------------------------------#


         #----- Define the age classes. ---------------------------------------------------#
         ageconow          = rep(x=agepa,times=ncohorts)
         #---------------------------------------------------------------------------------#



         #----- Read the cohort level variables. ------------------------------------------#
         pftconow          = mymont$PFT
         nplantconow       = mymont$NPLANT
         heightconow       = mymont$HITE
         wood.densconow    = pft$rho[pftconow]
         baconow           = mymont$BA.CO
         agbconow          = mymont$AGB.CO
         laiconow          = mymont$MMEAN.LAI.CO
         waiconow          = mymont$WAI.CO
         caiconow          = pmin(1.,nplantconow * dbh2ca(dbh=dbhconow,ipft=pftconow))
         taiconow          = laiconow + waiconow
         #------ Auxiliary variables for mean diurnal cycle. ------------------------------#
         q.pftconow        = matrix( data  = pftconow
                                   , nrow  = mymont$NCOHORTS.GLOBAL
                                   , ncol  = mymont$NDCYC
                                   , byrow = FALSE
                                   )#end matrix
         q.nplantconow     = matrix( data  = nplantconow
                                   , nrow  = mymont$NCOHORTS.GLOBAL
                                   , ncol  = mymont$NDCYC
                                   , byrow = FALSE
                                   )#end matrix
         q.laiconow        = matrix( data  = laiconow
                                   , nrow  = mymont$NCOHORTS.GLOBAL
                                   , ncol  = mymont$NDCYC
                                   , byrow = FALSE
                                   )#end matrix
         q.waiconow        = matrix( data  = waiconow
                                   , nrow  = mymont$NCOHORTS.GLOBAL
                                   , ncol  = mymont$NDCYC
                                   , byrow = FALSE
                                   )#end matrix
         q.taiconow        = matrix( data  = taiconow
                                   , nrow  = mymont$NCOHORTS.GLOBAL
                                   , ncol  = mymont$NDCYC
                                   , byrow = FALSE
                                   )#end matrix
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #     Find biomass of all tissues.                                                #
         #---------------------------------------------------------------------------------#
         bdeadconow        = mymont$BDEAD
         bleafconow        = mymont$MMEAN.BLEAF.CO
         bsapwoodconow     = mymont$BSAPWOODA+mymont$BSAPWOODB
         if (all(mymont$MMEAN.BROOT.CO == 0)){
            bfrootconow    = ( dbh2bl(dbh=dbhconow.lmon,ipft=pftconow)
                             * pft$qroot[pftconow] )
         }else{
            bfrootconow    = mymont$MMEAN.BROOT.CO
         }#end if
         bcrootconow       = mymont$BSAPWOODB + (1. - pft$agf.bs[pftconow]) * bdeadconow
         bstemconow        = mymont$BSAPWOODA +       pft$agf.bs[pftconow]  * bdeadconow
         brootconow        = bfrootconow + bcrootconow
         baliveconow       = bleafconow + bfrootconow + bsapwoodconow
         bstorageconow     = mymont$MMEAN.BSTORAGE.CO
         bseedsconow       = mymont$BSEEDS.CO
         biomassconow      = baliveconow + bstorageconow + bseedsconow + bdeadconow
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #      Find total NPP then total autotrophic, growth and storage respiration.     #
         # The latter two will be distributed amongst tissues.                             #
         #---------------------------------------------------------------------------------#
         #----- Monthly means. ------------------------------------------------------------#
         gppconow          = mymont$MMEAN.GPP.CO
         nppconow           = mymont$MMEAN.NPP.CO
         plant.respconow    = mymont$MMEAN.PLRESP.CO
         growth.respconow   = mymont$MMEAN.GROWTH.RESP.CO
         storage.respconow  = mymont$MMEAN.STORAGE.RESP.CO
         gs.respconow       = growth.respconow + storage.respconow
         #----- Mean diurnal cycle. -------------------------------------------------------#
         q.gppconow           = mymont$QMEAN.GPP.CO
         q.nppconow           = mymont$QMEAN.NPP.CO
         q.plant.respconow    = mymont$QMEAN.PLRESP.CO
         q.growth.respconow   = mymont$QMEAN.GROWTH.RESP.CO
         q.storage.respconow  = mymont$QMEAN.STORAGE.RESP.CO
         q.gs.respconow       = q.growth.respconow + q.storage.respconow
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Find the fractions that go to each pool.                                   #
         #---------------------------------------------------------------------------------#
         fg.leaf    = bleafconow  / ( baliveconow + bdeadconow )
         fg.stem    = bstemconow  / ( baliveconow + bdeadconow )
         fg.froot   = bfrootconow / ( baliveconow + bdeadconow )
         fg.croot   = bcrootconow / ( baliveconow + bdeadconow )
         fs.stem    = bstemconow  / ( bcrootconow + bstemconow )
         fs.croot   = bcrootconow / ( bcrootconow + bstemconow )
         q.fg.leaf  = matrix(data=fg.leaf ,nrow=mymont$NCOHORTS.GLOBAL,ncol=mymont$NDCYC)
         q.fg.stem  = matrix(data=fg.stem ,nrow=mymont$NCOHORTS.GLOBAL,ncol=mymont$NDCYC)
         q.fg.froot = matrix(data=fg.froot,nrow=mymont$NCOHORTS.GLOBAL,ncol=mymont$NDCYC)
         q.fg.croot = matrix(data=fg.croot,nrow=mymont$NCOHORTS.GLOBAL,ncol=mymont$NDCYC)
         q.fs.stem  = matrix(data=fs.stem ,nrow=mymont$NCOHORTS.GLOBAL,ncol=mymont$NDCYC)
         q.fs.croot = matrix(data=fs.croot,nrow=mymont$NCOHORTS.GLOBAL,ncol=mymont$NDCYC)
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Attribute respiration to the different pools.  Assuming that non-          #
         # structural carbon respiration 
         #---------------------------------------------------------------------------------#
         if ("MMEAN.LEAF.GROWTH.CO" %in% names(mymont)){
            #----- Mean monthly cycle. ----------------------------------------------------#
            leaf.respconow  = ( mymont$MMEAN.LEAF.RESP.CO
                              + mymont$MMEAN.LEAF.GROWTH.RESP.CO
                              + mymont$MMEAN.LEAF.STORAGE.RESP.CO
                              )#end leaf.respconow
            stem.respconow  = ( mymont$MMEAN.SAPA.GROWTH.RESP.CO
                              + mymont$MMEAN.SAPA.STORAGE.RESP.CO
                              )#end leaf.respconow
            froot.respconow = ( mymont$MMEAN.ROOT.RESP.CO
                              + mymont$MMEAN.ROOT.GROWTH.RESP.CO
                              + mymont$MMEAN.ROOT.STORAGE.RESP.CO
                              )#end leaf.respconow
            croot.respconow = ( mymont$MMEAN.SAPB.GROWTH.RESP.CO
                              + mymont$MMEAN.SAPB.STORAGE.RESP.CO
                              )#end leaf.respconow
            root.respconow  = froot.respconow + croot.respconow
            #----- Mean diurnal cycle. ----------------------------------------------------#
            q.leaf.respconow  = ( mymont$QMEAN.LEAF.RESP.CO
                                + mymont$QMEAN.LEAF.GROWTH.RESP.CO
                                + mymont$QMEAN.LEAF.STORAGE.RESP.CO
                                )#end q.leaf.respconow
            q.stem.respconow  = ( mymont$QMEAN.SAPA.GROWTH.RESP.CO
                                + mymont$QMEAN.SAPA.STORAGE.RESP.CO
                                )#end q.leaf.respconow
            q.froot.respconow = ( mymont$QMEAN.ROOT.RESP.CO
                                + mymont$QMEAN.ROOT.GROWTH.RESP.CO
                                + mymont$QMEAN.ROOT.STORAGE.RESP.CO
                                )#end q.leaf.respconow
            q.croot.respconow = ( mymont$QMEAN.SAPB.GROWTH.RESP.CO
                                + mymont$QMEAN.SAPB.STORAGE.RESP.CO
                                )#end q.leaf.respconow
            q.root.respconow  = q.froot.respconow + q.croot.respconow
         }else{
            #----- Mean monthly cycle. ----------------------------------------------------#
            leaf.respconow  = mymont$MMEAN.LEAF.RESP.CO    + fg.leaf  * growth.respconow
            stem.respconow  = fs.stem  * storage.respconow + fg.stem  * growth.respconow
            froot.respconow = mymont$MMEAN.ROOT.RESP.CO    + fg.froot * growth.respconow
            croot.respconow = fs.croot * storage.respconow + fg.croot * growth.respconow
            root.respconow  = froot.respconow + croot.respconow
            #----- Mean diurnal cycle. ----------------------------------------------------#
            q.leaf.respconow  = ( mymont$QMEAN.LEAF.RESP.CO   
                                + q.fg.leaf  * q.growth.respconow
                                )#end q.leaf.respconow
            q.stem.respconow  = ( q.fs.stem  * q.storage.respconow
                                + q.fg.stem  * q.growth.respconow
                                )#end q.stem.respconow
            q.froot.respconow = ( mymont$QMEAN.ROOT.RESP.CO
                                + q.fg.froot * q.growth.respconow
                                )#end q.froot.respconow
            q.croot.respconow = ( q.fs.croot * q.storage.respconow
                                + q.fg.croot * q.growth.respconow
                                )#end q.croot.respconow
            q.root.respconow  = q.froot.respconow + q.croot.respconow
         }#end if ("MMEAN.LEAF.GROWTH.CO" %in% names(mymont))
         #---------------------------------------------------------------------------------#



         #----- Flags to tell whether leaves and branchwood were resolvable. --------------#
         leaf.okconow   = mymont$MMEAN.LEAF.HCAP.CO %>=% pft$veg.hcap.min[pftconow  ]
         wood.okconow   = mymont$MMEAN.WOOD.HCAP.CO %>=% pft$veg.hcap.min[pftconow  ]
         q.leaf.okconow = mymont$QMEAN.LEAF.HCAP.CO %>=% pft$veg.hcap.min[q.pftconow]
         q.wood.okconow = mymont$QMEAN.WOOD.HCAP.CO %>=% pft$veg.hcap.min[q.pftconow]
         #---------------------------------------------------------------------------------#



         if (kludgecbal){
            cbaconow          = mymont$MMEAN.CB - (mondays - 1) * mymont$MMEAN.BSTORAGE.CO
            #------------------------------------------------------------------------------#
            #     Temporary fix to correct the carbon balance.                             #
            #------------------------------------------------------------------------------#
            if (one.cohort){
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
            #------------------------------------------------------------------------------#
         }else{
            cbaconow             = mymont$MMEAN.CB
            #------------------------------------------------------------------------------#
            #     Temporary fix to correct the carbon balance.                             #
            #------------------------------------------------------------------------------#
            if (one.cohort){
               cbalightconow     = mymont$CB.LIGHTMAX[thismonth]
               cbamoistconow     = mymont$CB.MOISTMAX[thismonth]
            }else{
               cbalightconow     = mymont$CB.LIGHTMAX[,thismonth]
               cbamoistconow     = mymont$CB.MOISTMAX[,thismonth]
            }#end if
            #------------------------------------------------------------------------------#
         }#end if
         #---------------------------------------------------------------------------------#

         cbamaxconow       = klight * cbalightconow + (1. - klight) * cbamoistconow
         cbarelconow       = mymont$CBR.BAR
         mcostconow        = ( mymont$MMEAN.LEAF.MAINTENANCE.CO
                             + mymont$MMEAN.ROOT.MAINTENANCE.CO ) * yr.day
         ldropconow        = mymont$MMEAN.LEAF.DROP.CO * yr.day


         sm.stressconow      = 1. - mymont$MMEAN.FS.OPEN.CO
         lightconow          = mymont$MMEAN.LIGHT.LEVEL.CO
         light.beamconow     = mymont$MMEAN.LIGHT.LEVEL.BEAM.CO
         light.diffconow     = mymont$MMEAN.LIGHT.LEVEL.DIFF.CO
         q.sm.stressconow    = 1. - mymont$QMEAN.FS.OPEN.CO
         q.lightconow        = mymont$QMEAN.LIGHT.LEVEL.CO
         q.light.beamconow   = mymont$QMEAN.LIGHT.LEVEL.BEAM.CO
         q.light.diffconow   = mymont$QMEAN.LIGHT.LEVEL.DIFF.CO
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #      Solve the change in storage .                                              #
         #---------------------------------------------------------------------------------#
         dcbadtconow       = nppconow - mcostconow - ldropconow
         #---------------------------------------------------------------------------------#


         #----- Allocation and productivity relative to the total living biomass. ---------#
         f.gppconow        =  100. * gppconow        / pmax(baliveconow,0.01)
         f.plant.respconow =  100. * plant.respconow / pmax(baliveconow,0.01)
         f.nppconow        =  100. * nppconow        / pmax(baliveconow,0.01)
         f.mcoconow        =  100. * mcostconow      / pmax(baliveconow,0.01)
         f.dcbadtconow     =  100. * dcbadtconow     / pmax(baliveconow,0.01)
         f.cbaconow        =         cbaconow        / pmax(baliveconow,0.01)
         f.bstorageconow   =         bstorageconow   / pmax(baliveconow,0.01)
         f.bleafconow      =         bleafconow      / pmax(baliveconow,0.01)
         f.bstemconow      =         bstemconow      / pmax(baliveconow,0.01)
         f.brootconow      =         brootconow      / pmax(baliveconow,0.01)
         f.bseedsconow     =         bseedsconow     / pmax(baliveconow,0.01)
         #---------------------------------------------------------------------------------#



         #----- Energy and water fluxes: convert them to plant area. ----------------------#
         hflxlcconow       = mymont$MMEAN.SENSIBLE.LC.CO           / nplantconow
         wflxlcconow       = mymont$MMEAN.VAPOR.LC.CO    * day.sec / nplantconow
         transpconow       = mymont$MMEAN.TRANSP.CO      * day.sec / nplantconow
         i.hflxlcconow     = ifelse( leaf.okconow
                                   , mymont$MMEAN.SENSIBLE.LC.CO           / laiconow
                                   , NA
                                   )#end ifelse
         i.wflxlcconow     = ifelse( leaf.okconow
                                   , mymont$MMEAN.VAPOR.LC.CO    * day.sec / laiconow
                                   , NA
                                   )#end ifelse
         i.transpconow     = ifelse( leaf.okconow
                                   , mymont$MMEAN.TRANSP.CO      * day.sec / laiconow
                                   , NA
                                   )#end ifelse
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Find the leaf interstitial space and boundary layer specific humidities to  #
         # convert conductance to kgW/m2/day.                                              #
         #---------------------------------------------------------------------------------#
         lpsiconow      = ifelse( test = leaf.okconow
                                , yes  = mymont$MMEAN.TRANSP.CO / laiconow
                                , no   = NA
                                )#end ifelse
         wpsiconow      = ifelse( test = wood.okconow
                                , yes  = mymont$MMEAN.VAPOR.WC / waiconow
                                , no   = NA
                                )#end ifelse
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
         #------ Mean diurnal cycle. ------------------------------------------------------#
         q.lpsiconow      = ifelse( test = q.leaf.okconow
                                  , yes  = mymont$QMEAN.TRANSP.CO / q.laiconow
                                  , no   = NA
                                  )#end ifelse
         q.wpsiconow      = ifelse( test = q.wood.okconow
                                  , yes  = mymont$QMEAN.VAPOR.WC.CO / q.waiconow
                                  , no   = NA
                                  )#end ifelse
         q.can.shv.conow  = matrix( data = rep( x     = mymont$QMEAN.CAN.SHV.PA
                                              , times = rep( x     = ncohorts
                                                           , times = mymont$NDCYC
                                                           )#end rep
                                              )#end rep
                                  , nrow  = mymont$NCOHORTS.GLOBAL
                                  , ncol  = mymont$NDCYC
                                  )#end matrix
         #---- Net conductance, combining stomatal and boundary layer. --------------------#
         q.leaf.gnw.conow = ifelse( test = mymont$QMEAN.LEAF.GBW.CO
                                         + mymont$QMEAN.LEAF.GSW.CO > 1.e-10
                                  , yes  = mymont$QMEAN.LEAF.GBW.CO 
                                         * mymont$QMEAN.LEAF.GSW.CO
                                         / ( mymont$QMEAN.LEAF.GBW.CO
                                           + mymont$QMEAN.LEAF.GSW.CO )
                                  , no   = pmin( mymont$QMEAN.LEAF.GBW.CO
                                               , mymont$QMEAN.LEAF.GSW.CO
                                               )#end pmin
                                )#end ifelse
         q.lbl.shv.conow  = ( q.can.shv.conow 
                            + q.lpsiconow / pmax(mymont$QMEAN.LEAF.GBW.CO,1.e-10)
                            )#end q.lbl.shv.conow
         q.wbl.shv.conow  = ( q.can.shv.conow 
                            + q.wpsiconow / pmax(mymont$QMEAN.WOOD.GBW.CO,1.e-10)
                            )#end q.wbl.shv.conow
         q.lis.shv.conow  = ( q.can.shv.conow 
                            + q.lpsiconow / pmax(q.leaf.gnw.conow,1.e-10)
                            )#end q.lis.shv.conow
         #---------------------------------------------------------------------------------#


         #----- Find the conductances in kgW/m2/day. --------------------------------------#
         leaf.gbwconow    = ( mymont$MMEAN.LEAF.GBW.CO * day.sec * ep
                            * (1 + epim1 * can.shv.conow) * (1 + epim1 * lbl.shv.conow)
                            )#end leaf.gbwconow
         leaf.gswconow    = ( mymont$MMEAN.LEAF.GSW.CO * day.sec * ep
                            * (1 + epim1 * lbl.shv.conow) * (1 + epim1 * lis.shv.conow)
                            )#end leaf.gswconow
         wood.gbwconow    = ( mymont$MMEAN.WOOD.GBW.CO * day.sec * ep
                            * (1 + epim1 * can.shv.conow) * (1 + epim1 * wbl.shv.conow)
                            )#end wood.gbwconow
         q.leaf.gbwconow  = ( mymont$QMEAN.LEAF.GBW.CO * day.sec * ep
                            * (1 + epim1 * q.can.shv.conow) * (1 + epim1 * q.lbl.shv.conow)
                            )#end q.leaf.gbwconow
         q.leaf.gswconow  = ( mymont$QMEAN.LEAF.GSW.CO * day.sec * ep
                            * (1 + epim1 * q.lbl.shv.conow) * (1 + epim1 * q.lis.shv.conow)
                            )#end q.leaf.gswconow
         q.wood.gbwconow  = ( mymont$QMEAN.WOOD.GBW.CO * day.sec * ep
                            * (1 + epim1 * q.can.shv.conow) * (1 + epim1 * q.wbl.shv.conow)
                            )#end q.wood.gbwconow
         #---------------------------------------------------------------------------------#


         #----- Find the net radiation for leaves (in m2/leaf!). --------------------------#
         par.mult           = Watts.2.Ein * 1.e6
         leaf.parconow      = ifelse( leaf.okconow
                                    , mymont$MMEAN.PAR.L.CO / laiconow *  par.mult
                                    , NA
                                    )#end ifelse
         leaf.par.beamconow = ifelse( leaf.okconow
                                    , mymont$MMEAN.PAR.L.BEAM.CO / laiconow *  par.mult
                                    , NA
                                    )#end ifelse
         leaf.par.diffconow = ifelse( leaf.okconow
                                    , mymont$MMEAN.PAR.L.DIFF.CO / laiconow *  par.mult
                                    , NA
                                    )#end ifelse
         leaf.rshortconow   = ifelse( leaf.okconow
                                    , mymont$MMEAN.RSHORT.L.CO / laiconow
                                    , NA
                                    )#end ifelse
         leaf.rlongconow    = ifelse( leaf.okconow
                                    , mymont$MMEAN.RLONG.L.CO  / laiconow
                                    , NA
                                    )#end ifelse
         leaf.gppconow      = ifelse( leaf.okconow
                                    , mymont$MMEAN.GPP.CO * nplantconow / laiconow
                                    , NA
                                    )#end ifelse
         #----- Mean diurnal cycle. -------------------------------------------------------#
         q.leaf.parconow      = ifelse( q.leaf.okconow
                                      , mymont$QMEAN.PAR.L.CO / q.laiconow *  par.mult
                                      , NA
                                      )#end ifelse
         q.leaf.par.beamconow = ifelse( q.leaf.okconow
                                      , mymont$QMEAN.PAR.L.BEAM.CO / q.laiconow *  par.mult
                                      , NA
                                      )#end ifelse
         q.leaf.par.diffconow = ifelse( q.leaf.okconow
                                      , mymont$QMEAN.PAR.L.DIFF.CO / q.laiconow *  par.mult
                                      , NA
                                      )#end ifelse
         q.leaf.rshortconow   = ifelse( q.leaf.okconow
                                      , mymont$QMEAN.RSHORT.L.CO / q.laiconow
                                      , NA
                                      )#end ifelse
         q.leaf.rlongconow    = ifelse( q.leaf.okconow
                                      , mymont$QMEAN.RLONG.L.CO  / q.laiconow
                                      , NA
                                      )#end ifelse
         q.leaf.gppconow      = ifelse( q.leaf.okconow
                                      , mymont$QMEAN.GPP.CO * q.nplantconow / q.laiconow
                                      , NA
                                      )#end ifelse
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Leaf/wood thermal properties.                                               #
         #---------------------------------------------------------------------------------#
         leaf.waterconow   = ifelse( test = leaf.okconow
                                   , yes  = mymont$MMEAN.LEAF.WATER.CO / laiconow
                                   , no   = NA
                                   )#end ifelse
         leaf.tempconow    = ifelse( test = leaf.okconow
                                   , yes  = mymont$MMEAN.LEAF.TEMP.CO  - t00
                                   , no   = NA
                                   )#end ifelse
         wood.tempconow    = ifelse( test = wood.okconow
                                   , yes  = mymont$MMEAN.WOOD.TEMP.CO  - t00
                                   , no   = NA
                                   )#end ifelse
         leaf.vpdconow     = ifelse( test = leaf.okconow
                                   , yes  = mymont$MMEAN.LEAF.VPDEF.CO  * 0.01
                                   , no   = NA
                                   )#end ifelse
         #----- Mean diurnal cycle. -------------------------------------------------------#
         q.leaf.waterconow   = ifelse( test = q.leaf.okconow
                                     , yes  = mymont$QMEAN.LEAF.WATER.CO / q.laiconow
                                     , no   = NA
                                     )#end ifelse
         q.leaf.tempconow    = ifelse( test = q.leaf.okconow
                                     , yes  = mymont$QMEAN.LEAF.TEMP.CO  - t00
                                     , no   = NA
                                     )#end ifelse
         q.wood.tempconow    = ifelse( test = q.wood.okconow
                                     , yes  = mymont$QMEAN.WOOD.TEMP.CO  - t00
                                     , no   = NA
                                     )#end ifelse
         q.leaf.vpdconow     = ifelse( test = q.leaf.okconow
                                     , yes  = mymont$QMEAN.LEAF.VPDEF.CO  * 0.01
                                     , no   = NA
                                     )#end ifelse
         #---------------------------------------------------------------------------------#



         #----- Assimilation rates by limitation. -----------------------------------------#
         assim.lightconow  = mymont$MMEAN.A.LIGHT.CO
         assim.rubpconow   = mymont$MMEAN.A.RUBP.CO
         assim.co2conow    = mymont$MMEAN.A.CO2.CO
         assim.ratioconow  = with(mymont, MMEAN.A.LIGHT.CO
                                        / pmax(1e-6,pmin(MMEAN.A.RUBP.CO,MMEAN.A.CO2.CO)))
         #----- Mean diurnal cycle. -------------------------------------------------------#
         q.assim.lightconow  = mymont$QMEAN.A.LIGHT.CO
         q.assim.rubpconow   = mymont$QMEAN.A.RUBP.CO
         q.assim.co2conow    = mymont$QMEAN.A.CO2.CO
         q.assim.ratioconow  = with(mymont, QMEAN.A.LIGHT.CO
                                          / pmax(1e-6,pmin(QMEAN.A.RUBP.CO,QMEAN.A.CO2.CO)))
         #---------------------------------------------------------------------------------#



         #------ Find the demand and supply by m2gnd. -------------------------------------#
         demandconow       = mymont$MMEAN.PSI.OPEN.CO     * laiconow * day.sec
         supplyconow       = mymont$MMEAN.WATER.SUPPLY.CO * day.sec
         q.demandconow     = mymont$QMEAN.PSI.OPEN.CO     * q.laiconow * day.sec
         q.supplyconow     = mymont$QMEAN.WATER.SUPPLY.CO * day.sec
         #---------------------------------------------------------------------------------#


         #------ Find the demographic rates. ----------------------------------------------#
         if (one.cohort){
            mortconow    = sum(mymont$MMEAN.MORT.RATE.CO)
            mortconow    = max(0,mortconow)
         }else{
            mortconow    = try(rowSums(mymont$MMEAN.MORT.RATE.CO))
            if ("try-error" %in% is(mortconow)) browser()
            mortconow    = pmax(0,mortconow)
         }#end if
         ncbmortconow    = pmax(0,mymont$MMEAN.MORT.RATE.CO[,2])
         dimortconow     = pmax(0,mortconow - ncbmortconow)
         recruitconow    = mymont$RECRUIT.DBH
         growthconow     = pmax(0,mymont$DLNDBH.DT)
         agb.growthconow = pmax(0,mymont$DLNAGB.DT)
         bsa.growthconow = pmax(0,mymont$DLNBA.DT )
         #---------------------------------------------------------------------------------#


         #------ Find the AGB and basal area of the previous month. -----------------------#
         agbcolmon        = agbconow * exp(-agb.growthconow/12.)
         bacolmon         = baconow  * exp(-bsa.growthconow/12.)
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #     Use efficiencies.  Do not calculate if transpiration is insignificant.      #
         #---------------------------------------------------------------------------------#
         tfineconow        = transpconow*yr.day >= 1.0
         rfineconow        = 12*emean$rain[m]   >= 1.0
         etconow           = transpconow + wflxlcconow
         #------ Find Rainfall use Efficiencies. ------------------------------------------#
         rueconow  = 1000. * nppconow    / ifelse( rfineconow, 12. * emean$rain[m]  , NA )
         wueconow  = 1000. * nppconow    / ifelse( tfineconow, transpconow * yr.day , NA )
         etueconow = 1000. * nppconow    / ifelse( tfineconow, etconow     * yr.day , NA )
         cueconow  =         nppconow    / ifelse( tfineconow, gppconow             , NA )
         ecueconow =         dcbadtconow / ifelse( tfineconow, gppconow             , NA )
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #   Canopy height and open canopy fraction.                                       #
         #---------------------------------------------------------------------------------#
         opencanconow            = unlist( tapply( X     = 1 - caiconow
                                                 , INDEX = ipaconow
                                                 , FUN   = cumprod
                                                 )#end tapply
                                         )#end unlist
         names(opencanconow)     = NULL
         zeroconow               = is.finite(opencanconow) & opencanconow <= tiny.num
         opencanconow[zeroconow] = 0.
         #---------------------------------------------------------------------------------#



         #----- Find some averages for photoperiod. ---------------------------------------#
         if (polar.night){
            phap.lparconow     = NA + pftconow
            phap.ltempconow    = NA + pftconow
            phap.lwaterconow   = NA + pftconow
            phap.lvpdconow     = NA + pftconow
            phap.fs.openconow  = NA + pftconow
            phap.lpsiconow     = NA + pftconow
            phap.leaf.gbaconow = NA + pftconow
            phap.leaf.gsaconow = NA + pftconow
            phap.can.shv.conow = NA + pftconow
         }else if (one.cohort){
            phap.lparconow     = mean(mymont$QMEAN.PAR.L.CO     [phap])
            phap.ltempconow    = mean(mymont$QMEAN.LEAF.TEMP.CO [phap])
            phap.lwaterconow   = mean(mymont$QMEAN.LEAF.WATER.CO[phap])
            phap.lvpdconow     = mean(mymont$QMEAN.LEAF.VPDEF.CO[phap])
            phap.fs.openconow  = mean(mymont$QMEAN.FS.OPEN.CO   [phap])
            phap.lpsiconow     = mean(mymont$QMEAN.TRANSP.CO    [phap])
            phap.leaf.gbaconow = mean(mymont$QMEAN.LEAF.GBW.CO  [phap])
            phap.leaf.gsaconow = mean(mymont$QMEAN.LEAF.GSW.CO  [phap])
            phap.can.shv.conow = rep(mean(mymont$QMEAN.CAN.SHV.PA[phap]),times=ncohorts)
         }else{
            phap.lparconow     = rowMeans(mymont$QMEAN.PAR.L.CO     [,phap])
            phap.ltempconow    = rowMeans(mymont$QMEAN.LEAF.TEMP.CO [,phap])
            phap.lwaterconow   = rowMeans(mymont$QMEAN.LEAF.WATER.CO[,phap])
            phap.lvpdconow     = rowMeans(mymont$QMEAN.LEAF.VPDEF.CO[,phap])
            phap.fs.openconow  = rowMeans(mymont$QMEAN.FS.OPEN.CO   [,phap])
            phap.lpsiconow     = rowMeans(mymont$QMEAN.TRANSP.CO    [,phap])
            phap.leaf.gbaconow = rowMeans(mymont$QMEAN.LEAF.GBW.CO  [,phap])
            phap.leaf.gsaconow = rowMeans(mymont$QMEAN.LEAF.GSW.CO  [,phap])
            if (one.patch){
               phap.can.shv.conow = rep( x     = mean(mymont$QMEAN.CAN.SHV.PA[phap])
                                       , times = ncohorts
                                       )#end rep
            }else{
               phap.can.shv.conow = rep( x     = rowMeans(mymont$QMEAN.CAN.SHV.PA[,phap])
                                       , times = ncohorts
                                       )#end rep
            }#end if
            #------------------------------------------------------------------------------#
         }#end if
         phap.lparconow      = ifelse( leaf.okconow
                                     , phap.lparconow  / laiconow * Watts.2.Ein * 1.e6
                                     , NA
                                     )#end ifelse
         phap.ltempconow     = ifelse( leaf.okconow, phap.ltempconow - t00      , NA )
         phap.lwaterconow    = ifelse( leaf.okconow, phap.lwaterconow / laiconow, NA )
         phap.lvpdconow      = ifelse( leaf.okconow, phap.lvpdconow * 0.01      , NA )
         phap.smsconow       = ifelse( leaf.okconow, 1 - phap.fs.openconow      , NA )
         phap.lpsiconow      = ifelse( leaf.okconow, phap.lpsiconow / laiconow  , NA )
         phap.leaf.gbaconow  = ifelse( leaf.okconow, phap.leaf.gbaconow         , NA )
         phap.leaf.gsaconow  = ifelse( leaf.okconow, phap.leaf.gsaconow         , NA )
         phap.can.shv.conow  = ifelse( leaf.okconow, phap.can.shv.conow         , NA )
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Find the leaf interstitial space and boundary layer specific humidities to  #
         # convert conductance to kgW/m2/day.                                              #
         #---------------------------------------------------------------------------------#
         #---- Net conductance, combining stomatal and boundary layer. --------------------#
         fine.cond          = is.finite(phap.leaf.gbaconow) & is.finite(phap.leaf.gsaconow)
         enough.cond        = ( fine.cond
                              & ( phap.leaf.gbaconow + phap.leaf.gsaconow ) > 1.e-10 )
         phap.leaf.gnaconow = ifelse( fine.cond
                                    , ifelse( enough.cond
                                            , ( phap.leaf.gbaconow * phap.leaf.gsaconow )
                                            / ( phap.leaf.gbaconow + phap.leaf.gsaconow )
                                            , pmin(phap.leaf.gbaconow,phap.leaf.gsaconow)
                                            )#end ifelse
                                    , NA
                                    )#end ifelse
         phap.lbl.shv.conow = ( phap.can.shv.conow
                              + phap.lpsiconow / pmax(phap.leaf.gbaconow,1.e-10) )
         phap.lis.shv.conow = ( phap.can.shv.conow
                              + phap.lpsiconow / pmax(phap.leaf.gnaconow,1.e-10) )
         #---------------------------------------------------------------------------------#


         #----- Find the conductances in kgW/m2/day. --------------------------------------#
         phap.lgbwconow  = ( phap.leaf.gbaconow * day.sec * ep
                           * (1 + epim1 * phap.can.shv.conow) 
                           * (1 + epim1 * phap.lbl.shv.conow)
                           )#end phap.lgbwconow
         phap.lgswconow  = ( phap.leaf.gsaconow * day.sec * ep
                           * (1 + epim1 * phap.lbl.shv.conow) 
                           * (1 + epim1 * phap.lis.shv.conow)
                           )#end phap.lgswconow
         #---------------------------------------------------------------------------------#




      }else{
         #----- Make everything NA. -------------------------------------------------------#
         ipaconow            = NA
         icoconow            = NA
         areaconow           = NA
         luconow             = NA
         dbhconow            = NA
         dbhcut              = NA
         dbhlevs             = NA
         dbhfac              = NA
         dbhconow.1ago       = NA
         dbhcut.1ago         = NA
         dbhlevs.1ago        = NA
         dbhfac.1ago         = NA
         dbhconow.lmon       = NA
         dbhcut.lmon         = NA
         dbhlevs.lmon        = NA
         dbhfac.lmon         = NA
         ageconow            = NA
         pftconow            = NA
         nplantconow         = NA
         heightconow         = NA
         wood.densconow      = NA
         baconow             = NA
         agbconow            = NA
         biomassconow        = NA
         laiconow            = NA
         waiconow            = NA
         taiconow            = NA
         gppconow            = NA
         leaf.respconow      = NA
         stem.respconow      = NA
         root.respconow      = NA
         froot.respconow     = NA
         croot.respconow     = NA
         growth.respconow    = NA
         storage.respconow   = NA
         plant.respconow     = NA
         assim.lightconow    = NA
         assim.rubpconow     = NA
         assim.co2conow      = NA
         assim.ratioconow    = NA
         nppconow            = NA
         cbaconow            = NA
         cbamaxconow         = NA
         cbalightconow       = NA
         cbamoistconow       = NA
         cbarelconow         = NA
         mcostconow          = NA
         ldropconow          = NA
         dcbadtconow         = NA
         sm.stressconow      = NA
         lightconow          = NA
         light.beamconow     = NA
         light.diffconow     = NA
         baliveconow         = NA
         bdeadconow          = NA
         bleafconow          = NA
         bsapwoodconow       = NA
         bfrootconow         = NA
         bcrootconow         = NA
         brootconow          = NA
         bstemconow          = NA
         bstorageconow       = NA
         bseedsconow         = NA
         hflxlcconow         = NA
         wflxlcconow         = NA
         transpconow         = NA
         i.hflxlcconow       = NA
         i.wflxlcconow       = NA
         i.transpconow       = NA
         wueconow            = NA
         cueconow            = NA
         ecueconow           = NA
         etueconow           = NA
         leaf.tempconow      = NA
         leaf.waterconow     = NA
         wood.tempconow      = NA
         leaf.vpdconow       = NA
         demandconow         = NA
         supplyconow         = NA
         mortconow           = NA
         ncbmortconow        = NA
         dimortconow         = NA
         recruitconow        = NA
         growthconow         = NA
         agb.growthconow     = NA
         bsa.growthconow     = NA
         leaf.gbwconow       = NA
         leaf.gswconow       = NA
         wood.gbwconow       = NA
         f.gppconow          = NA
         f.plant.respconow   = NA
         f.nppconow          = NA
         f.mcoconow          = NA
         f.dcbadtconow       = NA
         f.cbaconow          = NA
         f.bstorageconow     = NA
         f.bleafconow        = NA
         f.bstemconow        = NA
         f.brootconow        = NA
         f.bseedsconow       = NA
         leaf.parconow       = NA
         leaf.par.beamconow  = NA
         leaf.par.diffconow  = NA
         leaf.rshortconow    = NA
         leaf.rlongconow     = NA
         leaf.gppconow       = NA
         rueconow            = NA
         opencanconow        = NA
         useconow            = NA
         phap.lparconow      = NA
         phap.ltempconow     = NA
         phap.lwaterconow    = NA
         phap.lvpdconow      = NA
         phap.smsconow       = NA
         phap.lgbwconow      = NA
         phap.lgswconow      = NA
      }#end if
      #------------------------------------------------------------------------------------#
      #====================================================================================#
      #====================================================================================#
      #====================================================================================#
      #====================================================================================#






      #====================================================================================#
      #====================================================================================#
      #====================================================================================#
      #====================================================================================#
      #     Patch-level variables.                                                         #
      #------------------------------------------------------------------------------------#
      plab = paste0("y",sprintf("%4.4i",thisyear),"m",sprintf("%2.2i",thismonth))
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
      patch$wflxca       [[plab]] = - mymont$MMEAN.VAPOR.AC.PA    * day.sec
      patch$wflxgc       [[plab]] =   mymont$MMEAN.VAPOR.GC.PA    * day.sec
      patch$ustar        [[plab]] =   mymont$MMEAN.USTAR.PA
      patch$albedo       [[plab]] =   mymont$MMEAN.ALBEDO.PA
      patch$rshortup     [[plab]] =   mymont$MMEAN.RSHORTUP.PA
      patch$rlongup      [[plab]] =   mymont$MMEAN.RLONGUP.PA
      patch$parup        [[plab]] =   mymont$MMEAN.PARUP.PA       * Watts.2.Ein * 1e6
      patch$rshort.gnd   [[plab]] =   mymont$MMEAN.RSHORT.GND.PA
      patch$par.gnd      [[plab]] =   mymont$MMEAN.PAR.GND.PA     * Watts.2.Ein * 1e6
      patch$rnet         [[plab]] =   mymont$MMEAN.RNET.PA
      #----- Bind the current mean diurnal cycle patch. -----------------------------------#
      qpatch$nep          [[plab]] =   mymont$QMEAN.NEP.PA
      qpatch$het.resp     [[plab]] =   mymont$QMEAN.RH.PA
      qpatch$can.temp     [[plab]] =   mymont$QMEAN.CAN.TEMP.PA    - t00
      qpatch$gnd.temp     [[plab]] =   mymont$QMEAN.GND.TEMP.PA    - t00
      qpatch$can.shv      [[plab]] =   mymont$QMEAN.CAN.SHV.PA     * 1000.
      qpatch$gnd.shv      [[plab]] =   mymont$QMEAN.GND.SHV.PA     * 1000.
      qpatch$can.vpd      [[plab]] =   mymont$QMEAN.CAN.VPDEF.PA   * 0.01
      qpatch$can.co2      [[plab]] =   mymont$QMEAN.CAN.CO2.PA
      qpatch$can.prss     [[plab]] =   mymont$QMEAN.CAN.PRSS.PA    * 0.01
      qpatch$cflxca       [[plab]] = - mymont$QMEAN.CARBON.AC.PA
      qpatch$cflxst       [[plab]] =   mymont$QMEAN.CARBON.ST.PA
      qpatch$nee          [[plab]] = ( mymont$QMEAN.CARBON.ST.PA
                                     - mymont$QMEAN.CARBON.AC.PA )
      qpatch$hflxca       [[plab]] = - mymont$QMEAN.SENSIBLE.AC.PA
      qpatch$hflxgc       [[plab]] =   mymont$QMEAN.SENSIBLE.GC.PA
      qpatch$qwflxca      [[plab]] = - mymont$QMEAN.VAPOR.AC.PA    * qmean.can.alvli.pa
      qpatch$wflxca       [[plab]] = - mymont$QMEAN.VAPOR.AC.PA    * day.sec
      qpatch$wflxgc       [[plab]] =   mymont$QMEAN.VAPOR.GC.PA    * day.sec
      qpatch$ustar        [[plab]] =   mymont$QMEAN.USTAR.PA
      qpatch$albedo       [[plab]] =   mymont$QMEAN.ALBEDO.PA
      qpatch$rshortup     [[plab]] =   mymont$QMEAN.RSHORTUP.PA
      qpatch$rlongup      [[plab]] =   mymont$QMEAN.RLONGUP.PA
      qpatch$parup        [[plab]] =   mymont$QMEAN.PARUP.PA       * Watts.2.Ein * 1e6
      qpatch$rshort.gnd   [[plab]] =   mymont$QMEAN.RSHORT.GND.PA
      qpatch$par.gnd      [[plab]] =   mymont$QMEAN.PAR.GND.PA     * Watts.2.Ein * 1e6
      qpatch$rnet         [[plab]] =   mymont$QMEAN.RNET.PA
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Initialise patch-level properties that are derived from cohort-level.          #
      #------------------------------------------------------------------------------------#
      patch$lai          [[plab]] = rep(0.,times=mymont$NPATCHES.GLOBAL)
      patch$wai          [[plab]] = rep(0.,times=mymont$NPATCHES.GLOBAL)
      patch$agb          [[plab]] = rep(0.,times=mymont$NPATCHES.GLOBAL)
      patch$ba           [[plab]] = rep(0.,times=mymont$NPATCHES.GLOBAL)
      patch$wood.dens    [[plab]] = rep(NA,times=mymont$NPATCHES.GLOBAL)
      patch$can.depth    [[plab]] = rep(0.,times=mymont$NPATCHES.GLOBAL)
      patch$can.area     [[plab]] = rep(0.,times=mymont$NPATCHES.GLOBAL)
      patch$sm.stress    [[plab]] = rep(NA,times=mymont$NPATCHES.GLOBAL)
      patch$leaf.temp    [[plab]] = mymont$MMEAN.CAN.TEMP.PA  - t00
      patch$leaf.water   [[plab]] = rep(0.,times=mymont$NPATCHES.GLOBAL)
      patch$leaf.vpd     [[plab]] = mymont$MMEAN.CAN.VPDEF.PA * 0.01
      patch$wood.temp    [[plab]] = mymont$MMEAN.CAN.TEMP.PA  - t00
      patch$par.leaf     [[plab]] = rep(0.,times=mymont$NPATCHES.GLOBAL)
      patch$par.leaf.beam[[plab]] = rep(0.,times=mymont$NPATCHES.GLOBAL)
      patch$par.leaf.diff[[plab]] = rep(0.,times=mymont$NPATCHES.GLOBAL)
      patch$phap.lpar    [[plab]] = rep(NA,times=mymont$NPATCHES.GLOBAL)
      patch$phap.ltemp   [[plab]] = rep(NA,times=mymont$NPATCHES.GLOBAL)
      patch$phap.lwater  [[plab]] = rep(NA,times=mymont$NPATCHES.GLOBAL)
      patch$phap.lvpd    [[plab]] = rep(NA,times=mymont$NPATCHES.GLOBAL)
      patch$phap.sms     [[plab]] = rep(NA,times=mymont$NPATCHES.GLOBAL)
      patch$phap.lgbw    [[plab]] = rep(NA,times=mymont$NPATCHES.GLOBAL)
      patch$phap.lgsw    [[plab]] = rep(NA,times=mymont$NPATCHES.GLOBAL)
      patch$leaf.gpp     [[plab]] = rep(NA,times=mymont$NPATCHES.GLOBAL)
      patch$leaf.gsw     [[plab]] = rep(NA,times=mymont$NPATCHES.GLOBAL)
      patch$leaf.par     [[plab]] = rep(NA,times=mymont$NPATCHES.GLOBAL)
      patch$leaf.par.beam[[plab]] = rep(NA,times=mymont$NPATCHES.GLOBAL)
      patch$leaf.par.diff[[plab]] = rep(NA,times=mymont$NPATCHES.GLOBAL)
      patch$assim.light  [[plab]] = rep(NA,times=mymont$NPATCHES.GLOBAL)
      patch$assim.rubp   [[plab]] = rep(NA,times=mymont$NPATCHES.GLOBAL)
      patch$assim.co2    [[plab]] = rep(NA,times=mymont$NPATCHES.GLOBAL)
      patch$gpp          [[plab]] = rep(0.,times=mymont$NPATCHES.GLOBAL)
      patch$npp          [[plab]] = rep(0.,times=mymont$NPATCHES.GLOBAL)
      patch$cba          [[plab]] = rep(0.,times=mymont$NPATCHES.GLOBAL)
      patch$plant.resp   [[plab]] = rep(0.,times=mymont$NPATCHES.GLOBAL)
      patch$hflxlc       [[plab]] = rep(0.,times=mymont$NPATCHES.GLOBAL)
      patch$hflxwc       [[plab]] = rep(0.,times=mymont$NPATCHES.GLOBAL)
      patch$wflxlc       [[plab]] = rep(0.,times=mymont$NPATCHES.GLOBAL)
      patch$wflxwc       [[plab]] = rep(0.,times=mymont$NPATCHES.GLOBAL)
      patch$transp       [[plab]] = rep(0.,times=mymont$NPATCHES.GLOBAL)
      patch$soil.resp    [[plab]] = mymont$MMEAN.RH.PA
      patch$fast.soil.c  [[plab]] = mymont$MMEAN.FAST.SOIL.C
      patch$slow.soil.c  [[plab]] = mymont$MMEAN.SLOW.SOIL.C
      patch$struct.soil.c[[plab]] = mymont$MMEAN.STRUCT.SOIL.C
      #------ Mean diurnal cycle. ---------------------------------------------------------#
      zero.qpatch = matrix(data=0., nrow=mymont$NPATCHES.GLOBAL,ncol=mymont$NDCYC)
      na.qpatch   = matrix(data=NA, nrow=mymont$NPATCHES.GLOBAL,ncol=mymont$NDCYC)
      qpatch$sm.stress    [[plab]] = na.qpatch
      qpatch$leaf.temp    [[plab]] = mymont$QMEAN.CAN.TEMP.PA  - t00
      qpatch$leaf.water   [[plab]] = zero.qpatch
      qpatch$leaf.vpd     [[plab]] = mymont$QMEAN.CAN.VPDEF.PA * 0.01
      qpatch$wood.temp    [[plab]] = mymont$QMEAN.CAN.TEMP.PA  - t00
      qpatch$par.leaf     [[plab]] = zero.qpatch
      qpatch$par.leaf.beam[[plab]] = zero.qpatch
      qpatch$par.leaf.diff[[plab]] = zero.qpatch
      qpatch$leaf.gpp     [[plab]] = na.qpatch
      qpatch$leaf.gsw     [[plab]] = na.qpatch
      qpatch$leaf.par     [[plab]] = na.qpatch
      qpatch$leaf.par.beam[[plab]] = na.qpatch
      qpatch$leaf.par.diff[[plab]] = na.qpatch
      qpatch$assim.light  [[plab]] = na.qpatch
      qpatch$assim.rubp   [[plab]] = na.qpatch
      qpatch$assim.co2    [[plab]] = na.qpatch
      qpatch$gpp          [[plab]] = zero.qpatch
      qpatch$npp          [[plab]] = zero.qpatch
      qpatch$plant.resp   [[plab]] = zero.qpatch
      qpatch$hflxlc       [[plab]] = zero.qpatch
      qpatch$hflxwc       [[plab]] = zero.qpatch
      qpatch$wflxlc       [[plab]] = zero.qpatch
      qpatch$wflxwc       [[plab]] = zero.qpatch
      qpatch$transp       [[plab]] = zero.qpatch
      qpatch$soil.resp    [[plab]] = mymont$QMEAN.RH.PA
      #------------------------------------------------------------------------------------#


      if (any(ncohorts >0)){
         #---------------------------------------------------------------------------------#
         #     Monthly means -- no mean diurnal cycle.                                     #
         #---------------------------------------------------------------------------------#

         #----- Find some auxiliary patch-level properties. -------------------------------#
         lai.pa           = tapply( X     = mymont$MMEAN.LAI.CO
                                  , INDEX = ipaconow
                                  , FUN   = sum
                                  )#end tapply
         wai.pa           = tapply( X     = mymont$WAI.CO
                                  , INDEX = ipaconow
                                  , FUN   = sum
                                  )#end tapply
         leaf.energy.pa   = tapply( X     = mymont$MMEAN.LEAF.ENERGY.CO
                                  , INDEX = ipaconow
                                  , FUN   = sum
                                  )#end tapply
         leaf.water.pa    = tapply( X     = mymont$MMEAN.LEAF.WATER.CO
                                  , INDEX = ipaconow
                                  , FUN   = sum
                                  )#end tapply
         leaf.hcap.pa     = tapply( X     = mymont$MMEAN.LEAF.HCAP.CO
                                  , INDEX = ipaconow
                                  , FUN   = sum
                                  )#end tapply
         par.leaf.pa      = tapply( X     = mymont$MMEAN.PAR.L.CO
                                  , INDEX = ipaconow
                                  , FUN   = sum
                                  )#end tapply
         par.leaf.beam.pa = tapply( X     = mymont$MMEAN.PAR.L.BEAM.CO
                                  , INDEX = ipaconow
                                  , FUN   = sum
                                  )#end tapply
         par.leaf.diff.pa = tapply( X     = mymont$MMEAN.PAR.L.DIFF.CO
                                  , INDEX = ipaconow
                                  , FUN   = sum
                                  )#end tapply
         wood.energy.pa   = tapply( X     = mymont$MMEAN.WOOD.ENERGY.CO
                                  , INDEX = ipaconow
                                  , FUN   = sum
                                  )#end tapply
         wood.water.pa    = tapply( X     = mymont$MMEAN.WOOD.WATER.CO
                                  , INDEX = ipaconow
                                  , FUN   = sum
                                  )#end tapply
         wood.hcap.pa     = tapply( X     = mymont$MMEAN.WOOD.HCAP.CO
                                  , INDEX = ipaconow
                                  , FUN   = sum
                                  )#end tapply
         #----- Make root respiration extensive. ------------------------------------------#
         root.resp.pa     = tapply( X     = root.respconow*nplantconow
                                  , INDEX = ipaconow
                                  , FUN   = sum
                                  )#end tapply
         #---------------------------------------------------------------------------------#




         #---------------------------------------------------------------------------------#
         #       Find the polygon-average depth and area.                                  #
         #---------------------------------------------------------------------------------#
         useconow     = as.numeric(opencanconow > tiny.num)
         xconow       = heightconow  * nplantconow * baconow * useconow
         wconow       = nplantconow  * baconow * useconow
         oconow       = opencanconow * useconow
         can.depth.pa = (    tapply(X=xconow,INDEX=ipaconow,FUN=sum,na.rm=TRUE)
                        /    tapply(X=wconow,INDEX=ipaconow,FUN=sum,na.rm=TRUE)
                        )#end can.depth.pa
         can.area.pa  = 1. - tapply(X=oconow,INDEX=ipaconow,FUN=min,na.rm=TRUE)
         #---------------------------------------------------------------------------------#



         #----- Find the temperature and liquid fraction of leaf and wood. ----------------#
         leaf.empty                = leaf.hcap.pa == 0
         wood.empty                = wood.hcap.pa == 0
         leaf.temp.pa              = uextcm2tl( uext    = leaf.energy.pa
                                              , wmass   = leaf.water.pa
                                              , dryhcap = leaf.hcap.pa   )$temp - t00
         wood.temp.pa              = uextcm2tl( uext    = wood.energy.pa
                                              , wmass   = wood.water.pa
                                              , dryhcap = wood.hcap.pa   )$temp - t00
         leaf.water.pa             = leaf.water.pa / lai.pa
         leaf.temp.pa [leaf.empty] = NA
         leaf.water.pa[leaf.empty] = NA
         wood.temp.pa [wood.empty] = NA
         #---------------------------------------------------------------------------------#





         #----- Find the variables that must be rendered extensive. -----------------------#
         agb.pa        = tapply(X=agbconow       *nplantconow,INDEX=ipaconow,FUN=sum)
         ba.pa         = tapply(X=baconow        *nplantconow,INDEX=ipaconow,FUN=sum)
         gpp.pa        = tapply(X=gppconow       *nplantconow,INDEX=ipaconow,FUN=sum)
         npp.pa        = tapply(X=nppconow       *nplantconow,INDEX=ipaconow,FUN=sum)
         cba.pa        = tapply(X=cbaconow       *nplantconow,INDEX=ipaconow,FUN=sum)
         plant.resp.pa = tapply(X=plant.respconow*nplantconow,INDEX=ipaconow,FUN=sum)
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
         #      Wood density is found using weighted averages (basal areas are the         #
         # weights).                                                                       #
         #---------------------------------------------------------------------------------#
         wood.dens.pa = mapply( FUN      = weighted.mean
                              , x        = split(wood.densconow,ipaconow)
                              , w        = split(baconow       ,ipaconow)
                              , SIMPLIFY = TRUE
                              )#end mapply
         #---------------------------------------------------------------------------------#





         #---------------------------------------------------------------------------------#
         #      Aggregate variables that must be weighted by LAI.                          #
         #---------------------------------------------------------------------------------#
         leaf.gpp.pa      = mapply( FUN      = weighted.mean
                                  , x        = split(leaf.gppconow          ,ipaconow)
                                  , w        = split(laiconow               ,ipaconow)
                                  , na.rm    = TRUE
                                  , SIMPLIFY = TRUE
                                  )#end mapply
         leaf.par.pa      = mapply( FUN      = weighted.mean
                                  , x        = split(leaf.parconow          ,ipaconow)
                                  , w        = split(laiconow               ,ipaconow)
                                  , na.rm    = TRUE
                                  , SIMPLIFY = TRUE
                                  )#end mapply
         leaf.par.beam.pa = mapply( FUN      = weighted.mean
                                  , x        = split(leaf.par.beamconow     ,ipaconow)
                                  , w        = split(laiconow               ,ipaconow)
                                  , na.rm    = TRUE
                                  , SIMPLIFY = TRUE
                                  )#end mapply
         leaf.par.diff.pa = mapply( FUN      = weighted.mean
                                  , x        = split(leaf.par.diffconow     ,ipaconow)
                                  , w        = split(laiconow               ,ipaconow)
                                  , na.rm    = TRUE
                                  , SIMPLIFY = TRUE
                                  )#end mapply
         leaf.gsw.pa      = mapply( FUN      = weighted.mean
                                  , x        = split(leaf.gswconow          ,ipaconow)
                                  , w        = split(laiconow               ,ipaconow)
                                  , na.rm    = TRUE
                                  , SIMPLIFY = TRUE
                                  )#end mapply
         leaf.vpd.pa      = mapply( FUN      = weighted.mean
                                  , x        = split(leaf.vpdconow          ,ipaconow)
                                  , w        = split(laiconow               ,ipaconow)
                                  , na.rm    = TRUE
                                  , SIMPLIFY = TRUE
                                  )#end mapply
         sm.stress.pa     = mapply( FUN      = weighted.mean
                                  , x        = split(sm.stressconow         ,ipaconow)
                                  , w        = split(laiconow               ,ipaconow)
                                  , na.rm    = TRUE
                                  , SIMPLIFY = TRUE
                                  )#end mapply
         assim.light.pa   = mapply( FUN      = weighted.mean
                                  , x        = split(assim.lightconow       ,ipaconow)
                                  , w        = split(laiconow               ,ipaconow)
                                  , na.rm    = TRUE
                                  , SIMPLIFY = TRUE
                                  )#end mapply
         assim.rubp.pa    = mapply( FUN      = weighted.mean
                                  , x        = split(assim.rubpconow        ,ipaconow)
                                  , w        = split(laiconow               ,ipaconow)
                                  , na.rm    = TRUE
                                  , SIMPLIFY = TRUE
                                  )#end mapply
         assim.co2.pa     = mapply( FUN      = weighted.mean
                                  , x        = split(assim.co2conow         ,ipaconow)
                                  , w        = split(laiconow               ,ipaconow)
                                  , na.rm    = TRUE
                                  , SIMPLIFY = TRUE
                                  )#end mapply
         phap.lpar.pa     = mapply( FUN      = weighted.mean
                                  , x        = split(phap.lparconow         ,ipaconow)
                                  , w        = split(laiconow               ,ipaconow)
                                  , na.rm    = TRUE
                                  , SIMPLIFY = TRUE
                                  )#end mapply
         phap.ltemp.pa    = mapply( FUN      = weighted.mean
                                  , x        = split(phap.ltempconow        ,ipaconow)
                                  , w        = split(laiconow               ,ipaconow)
                                  , na.rm    = TRUE
                                  , SIMPLIFY = TRUE
                                  )#end mapply
         phap.lwater.pa   = mapply( FUN      = weighted.mean
                                  , x        = split(phap.lwaterconow       ,ipaconow)
                                  , w        = split(laiconow               ,ipaconow)
                                  , na.rm    = TRUE
                                  , SIMPLIFY = TRUE
                                  )#end mapply
         phap.lvpd.pa     = mapply( FUN      = weighted.mean
                                  , x        = split(phap.lvpdconow         ,ipaconow)
                                  , w        = split(laiconow               ,ipaconow)
                                  , na.rm    = TRUE
                                  , SIMPLIFY = TRUE
                                  )#end mapply
         phap.sms.pa      = mapply( FUN      = weighted.mean
                                  , x        = split(phap.smsconow          ,ipaconow)
                                  , w        = split(laiconow               ,ipaconow)
                                  , na.rm    = TRUE
                                  , SIMPLIFY = TRUE
                                  )#end mapply
         phap.lgbw.pa     = mapply( FUN      = weighted.mean
                                  , x        = split(phap.lgbwconow         ,ipaconow)
                                  , w        = split(laiconow               ,ipaconow)
                                  , na.rm    = TRUE
                                  , SIMPLIFY = TRUE
                                  )#end mapply
         phap.lgsw.pa     = mapply( FUN      = weighted.mean
                                  , x        = split(phap.lgswconow         ,ipaconow)
                                  , w        = split(laiconow               ,ipaconow)
                                  , na.rm    = TRUE
                                  , SIMPLIFY = TRUE
                                  )#end mapply
         #----- Discard data from empty cohorts.  -----------------------------------------#
         leaf.gpp.pa     [leaf.empty] = NA
         leaf.par.pa     [leaf.empty] = NA
         leaf.par.beam.pa[leaf.empty] = NA
         leaf.par.diff.pa[leaf.empty] = NA
         leaf.gsw.pa     [leaf.empty] = NA
         leaf.vpd.pa     [leaf.empty] = NA
         sm.stress.pa    [leaf.empty] = NA
         assim.light.pa  [leaf.empty] = NA
         assim.rubp.pa   [leaf.empty] = NA
         assim.co2.pa    [leaf.empty] = NA
         phap.lpar.pa    [leaf.empty] = NA
         phap.ltemp.pa   [leaf.empty] = NA
         phap.lwater.pa  [leaf.empty] = NA
         phap.lvpd.pa    [leaf.empty] = NA
         phap.sms.pa     [leaf.empty] = NA
         phap.lgbw.pa    [leaf.empty] = NA
         phap.lgsw.pa    [leaf.empty] = NA
         #---------------------------------------------------------------------------------#




         #---------------------------------------------------------------------------------#
         #     Copy the data back to the patch.                                            #
         #---------------------------------------------------------------------------------#
         idx.leaf                              = idx[! leaf.empty]
         idx.wood                              = idx[! wood.empty]
         patch$lai          [[plab]][idx     ] = lai.pa
         patch$wai          [[plab]][idx     ] = wai.pa
         patch$agb          [[plab]][idx     ] = agb.pa
         patch$ba           [[plab]][idx     ] = ba.pa
         patch$can.depth    [[plab]][idx     ] = can.depth.pa
         patch$can.area     [[plab]][idx     ] = can.area.pa
         patch$wood.dens    [[plab]][idx     ] = wood.dens.pa
         patch$par.leaf     [[plab]][idx     ] = par.leaf.pa
         patch$par.leaf.beam[[plab]][idx     ] = par.leaf.beam.pa
         patch$par.leaf.diff[[plab]][idx     ] = par.leaf.diff.pa
         patch$phap.lpar    [[plab]][idx.leaf] = phap.lpar.pa    [! leaf.empty]
         patch$phap.ltemp   [[plab]][idx.leaf] = phap.ltemp.pa   [! leaf.empty]
         patch$phap.lwater  [[plab]][idx.leaf] = phap.lwater.pa  [! leaf.empty]
         patch$phap.lvpd    [[plab]][idx.leaf] = phap.lvpd.pa    [! leaf.empty]
         patch$phap.sms     [[plab]][idx.leaf] = phap.sms.pa     [! leaf.empty]
         patch$phap.lgbw    [[plab]][idx.leaf] = phap.lgbw.pa    [! leaf.empty]
         patch$phap.lgsw    [[plab]][idx.leaf] = phap.lgsw.pa    [! leaf.empty]
         patch$sm.stress    [[plab]][idx.leaf] = sm.stress.pa    [! leaf.empty]
         patch$leaf.temp    [[plab]][idx.leaf] = leaf.temp.pa    [! leaf.empty]
         patch$leaf.water   [[plab]][idx.leaf] = leaf.water.pa   [! leaf.empty]
         patch$leaf.vpd     [[plab]][idx.leaf] = leaf.vpd.pa     [! leaf.empty]
         patch$leaf.gpp     [[plab]][idx.leaf] = leaf.gpp.pa     [! leaf.empty]
         patch$leaf.gsw     [[plab]][idx.leaf] = leaf.gsw.pa     [! leaf.empty]
         patch$leaf.par     [[plab]][idx.leaf] = leaf.par.pa     [! leaf.empty]
         patch$leaf.par.beam[[plab]][idx.leaf] = leaf.par.beam.pa[! leaf.empty]
         patch$leaf.par.diff[[plab]][idx.leaf] = leaf.par.diff.pa[! leaf.empty]
         patch$assim.light  [[plab]][idx.leaf] = assim.light.pa  [! leaf.empty]
         patch$assim.rubp   [[plab]][idx.leaf] = assim.rubp.pa   [! leaf.empty]
         patch$assim.co2    [[plab]][idx.leaf] = assim.co2.pa    [! leaf.empty]
         patch$wood.temp    [[plab]][idx.wood] = wood.temp.pa    [! wood.empty]
         patch$gpp          [[plab]][idx     ] = gpp.pa
         patch$npp          [[plab]][idx     ] = npp.pa
         patch$cba          [[plab]][idx     ] = cba.pa
         patch$plant.resp   [[plab]][idx     ] = plant.resp.pa
         patch$hflxlc       [[plab]][idx     ] = hflxlc.pa
         patch$hflxwc       [[plab]][idx     ] = hflxwc.pa
         patch$wflxlc       [[plab]][idx     ] = wflxlc.pa
         patch$wflxwc       [[plab]][idx     ] = wflxwc.pa
         patch$transp       [[plab]][idx     ] = transp.pa
         #------ Soil respiration mixes cohort (root) and patch (hetetrophic). ------------#
         patch$soil.resp [[plab]][idx] = patch$soil.resp [[plab]][idx] + root.resp.pa
         #---------------------------------------------------------------------------------#




         #---------------------------------------------------------------------------------#
         #---------------------------------------------------------------------------------#
         #---------------------------------------------------------------------------------#
         #    Mean diurnal cycle.                                                          #
         #---------------------------------------------------------------------------------#
         #---------------------------------------------------------------------------------#
         #---------------------------------------------------------------------------------#



         #----- Find some auxiliary patch-level properties. -------------------------------#
         q.lai.pa           = qapply( X     = q.laiconow
                                    , DIM   = 1
                                    , INDEX = ipaconow
                                    , FUN   = sum
                                    )#end tapply
         q.wai.pa           = qapply( X     = q.waiconow
                                    , DIM   = 1
                                    , INDEX = ipaconow
                                    , FUN   = sum
                                    )#end tapply
         q.leaf.energy.pa   = qapply( X     = mymont$QMEAN.LEAF.ENERGY.CO
                                    , DIM   = 1
                                    , INDEX = ipaconow
                                    , FUN   = sum
                                    )#end tapply
         q.leaf.water.pa    = qapply( X     = mymont$QMEAN.LEAF.WATER.CO
                                    , DIM   = 1
                                    , INDEX = ipaconow
                                    , FUN   = sum
                                    )#end tapply
         q.leaf.hcap.pa     = qapply( X     = mymont$QMEAN.LEAF.HCAP.CO
                                    , DIM   = 1
                                    , INDEX = ipaconow
                                    , FUN   = sum
                                    )#end tapply
         q.par.leaf.pa      = qapply( X     = mymont$QMEAN.PAR.L.CO
                                    , DIM   = 1
                                    , INDEX = ipaconow
                                    , FUN   = sum
                                    )#end tapply
         q.par.leaf.beam.pa = qapply( X     = mymont$QMEAN.PAR.L.BEAM.CO
                                    , DIM   = 1
                                    , INDEX = ipaconow
                                    , FUN   = sum
                                    )#end tapply
         q.par.leaf.diff.pa = qapply( X     = mymont$QMEAN.PAR.L.DIFF.CO
                                    , DIM   = 1
                                    , INDEX = ipaconow
                                    , FUN   = sum
                                    )#end tapply
         q.wood.energy.pa   = qapply( X     = mymont$QMEAN.WOOD.ENERGY.CO
                                    , DIM   = 1
                                    , INDEX = ipaconow
                                    , FUN   = sum
                                    )#end tapply
         q.wood.water.pa    = qapply( X     = mymont$QMEAN.WOOD.WATER.CO
                                    , DIM   = 1
                                    , INDEX = ipaconow
                                    , FUN   = sum
                                    )#end tapply
         q.wood.hcap.pa     = qapply( X     = mymont$QMEAN.WOOD.HCAP.CO
                                    , DIM   = 1
                                    , INDEX = ipaconow
                                    , FUN   = sum
                                    )#end tapply
         #----- Make root respiration extensive. ------------------------------------------#
         q.root.resp.pa     = qapply( X     = q.root.respconow * q.nplantconow
                                    , DIM   = 1
                                    , INDEX = ipaconow
                                    , FUN   = sum
                                    )#end tapply
         #---------------------------------------------------------------------------------#



         #----- Find the temperature and liquid fraction of leaf and wood. ----------------#
         q.leaf.empty    = q.leaf.hcap.pa == 0
         q.wood.empty    = q.wood.hcap.pa == 0
         q.leaf.temp.pa  = uextcm2tl( uext    = q.leaf.energy.pa
                                    , wmass   = q.leaf.water.pa
                                    , dryhcap = q.leaf.hcap.pa   )$temp - t00
         q.wood.temp.pa  = uextcm2tl( uext    = q.wood.energy.pa
                                    , wmass   = q.wood.water.pa
                                    , dryhcap = q.wood.hcap.pa   )$temp - t00
         q.leaf.water.pa = q.leaf.water.pa / q.lai.pa
         q.leaf.temp.pa  = ifelse( test = q.leaf.hcap.pa == 0
                                 , yes  = NA
                                 , no   = q.leaf.temp.pa
                                 )#end ifelse
         q.leaf.water.pa = ifelse( test = q.leaf.hcap.pa == 0
                                 , yes  = NA
                                 , no   = q.leaf.water.pa
                                 )#end ifelse
         q.wood.temp.pa  = ifelse( test = q.wood.hcap.pa == 0
                                 , yes  = NA
                                 , no   = q.wood.temp.pa
                                 )#end ifelse
         #---------------------------------------------------------------------------------#





         #----- Find the variables that must be rendered extensive. -----------------------#
         q.gpp.pa        = qapply( X     = q.gppconow * q.nplantconow
                                 , DIM   = 1
                                 , INDEX = ipaconow
                                 , FUN   = sum
                                 )#end qapply
         q.npp.pa        = qapply( X     = q.nppconow * q.nplantconow
                                 , DIM   = 1
                                 , INDEX = ipaconow
                                 , FUN   = sum
                                 )#end qapply
         q.plant.resp.pa = qapply( X     = q.plant.respconow * q.nplantconow
                                 , DIM   = 1
                                 , INDEX = ipaconow
                                 , FUN   = sum
                                 )#end qapply
         #---------------------------------------------------------------------------------#





         #----- Add the variables that are already extensive. -----------------------------#
         q.hflxlc.pa = qapply( X     = mymont$QMEAN.SENSIBLE.LC.CO
                             , DIM   = 1
                             , INDEX = ipaconow
                             , FUN   = sum
                             )#end tapply
         q.hflxwc.pa = qapply( X     = mymont$QMEAN.SENSIBLE.WC.CO
                             , DIM   = 1
                             , INDEX = ipaconow
                             , FUN   = sum
                             )#end tapply
         q.wflxlc.pa = qapply( X     = mymont$QMEAN.VAPOR.LC.CO  * day.sec
                             , DIM   = 1
                             , INDEX = ipaconow
                             , FUN   = sum
                             )#end tapply
         q.wflxwc.pa = qapply( X     = mymont$QMEAN.VAPOR.WC.CO  * day.sec
                             , DIM   = 1
                             , INDEX = ipaconow
                             , FUN   = sum
                             )#end tapply
         q.transp.pa = qapply( X     = mymont$QMEAN.TRANSP.CO    * day.sec
                             , DIM   = 1
                             , INDEX = ipaconow
                             , FUN   = sum
                             )#end tapply
         #---------------------------------------------------------------------------------#




         #---------------------------------------------------------------------------------#
         #      Aggregate variables that must be weighted by LAI.                          #
         #---------------------------------------------------------------------------------#
         q.leaf.gpp.pa      = t( mapply( FUN      = weighted.colMeans
                                       , X        = msplit(q.leaf.gppconow     ,ipaconow)
                                       , wr       = msplit(laiconow            ,ipaconow)
                                       , na.rm    = TRUE
                                       , SIMPLIFY = TRUE
                                       )#end mapply
                               )#end t
         q.leaf.par.pa      = t( mapply( FUN      = weighted.colMeans
                                       , X        = msplit(q.leaf.parconow     ,ipaconow)
                                       , w        = msplit(laiconow            ,ipaconow)
                                       , na.rm    = TRUE
                                       , SIMPLIFY = TRUE
                                       )#end mapply
                               )#end t
         q.leaf.par.beam.pa = t( mapply( FUN      = weighted.colMeans
                                       , X        = msplit(q.leaf.par.beamconow,ipaconow)
                                       , w        = msplit(laiconow            ,ipaconow)
                                       , na.rm    = TRUE
                                       , SIMPLIFY = TRUE
                                       )#end mapply
                               )#end t
         q.leaf.par.diff.pa = t( mapply( FUN      = weighted.colMeans
                                       , X        = msplit(q.leaf.par.diffconow,ipaconow)
                                       , w        = msplit(laiconow            ,ipaconow)
                                       , na.rm    = TRUE
                                       , SIMPLIFY = TRUE
                                       )#end mapply
                               )#end t
         q.leaf.gsw.pa      = t( mapply( FUN      = weighted.colMeans
                                       , X        = msplit(q.leaf.gswconow     ,ipaconow)
                                       , w        = msplit(laiconow            ,ipaconow)
                                       , na.rm    = TRUE
                                       , SIMPLIFY = TRUE
                                       )#end mapply
                               )#end t
         q.leaf.vpd.pa      = t( mapply( FUN      = weighted.colMeans
                                       , X        = msplit(q.leaf.vpdconow     ,ipaconow)
                                       , w        = msplit(laiconow            ,ipaconow)
                                       , na.rm    = TRUE
                                       , SIMPLIFY = TRUE
                                       )#end mapply
                               )#end t
         q.sm.stress.pa     = t( mapply( FUN      = weighted.colMeans
                                       , X        = msplit(q.sm.stressconow    ,ipaconow)
                                       , w        = msplit(laiconow            ,ipaconow)
                                       , na.rm    = TRUE
                                       , SIMPLIFY = TRUE
                                       )#end mapply
                               )#end t
         q.assim.light.pa   = t( mapply( FUN      = weighted.colMeans
                                       , X        = msplit(q.assim.lightconow  ,ipaconow)
                                       , w        = msplit(laiconow            ,ipaconow)
                                       , na.rm    = TRUE
                                       , SIMPLIFY = TRUE
                                       )#end mapply
                               )#end t
         q.assim.rubp.pa    = t( mapply( FUN      = weighted.colMeans
                                       , X        = msplit(q.assim.rubpconow   ,ipaconow)
                                       , w        = msplit(laiconow            ,ipaconow)
                                       , na.rm    = TRUE
                                       , SIMPLIFY = TRUE
                                       )#end mapply
                               )#end t
         q.assim.co2.pa     = t( mapply( FUN      = weighted.colMeans
                                       , X        = msplit(q.assim.co2conow    ,ipaconow)
                                       , w        = msplit(laiconow            ,ipaconow)
                                       , na.rm    = TRUE
                                       , SIMPLIFY = TRUE
                                       )#end mapply
                               )#end t
         #----- Discard data from empty cohorts.  -----------------------------------------#
         q.leaf.gpp.pa     [q.leaf.empty] = NA
         q.leaf.par.pa     [q.leaf.empty] = NA
         q.leaf.par.beam.pa[q.leaf.empty] = NA
         q.leaf.par.diff.pa[q.leaf.empty] = NA
         q.leaf.gsw.pa     [q.leaf.empty] = NA
         q.leaf.vpd.pa     [q.leaf.empty] = NA
         q.sm.stress.pa    [q.leaf.empty] = NA
         q.assim.light.pa  [q.leaf.empty] = NA
         q.assim.rubp.pa   [q.leaf.empty] = NA
         q.assim.co2.pa    [q.leaf.empty] = NA
         #---------------------------------------------------------------------------------#




         #---------------------------------------------------------------------------------#
         #     Copy the data back to the patch.                                            #
         #---------------------------------------------------------------------------------#
         qpatch$par.leaf     [[plab]][idx     ,] = q.par.leaf.pa
         qpatch$par.leaf.beam[[plab]][idx     ,] = q.par.leaf.beam.pa
         qpatch$par.leaf.diff[[plab]][idx     ,] = q.par.leaf.diff.pa
         qpatch$sm.stress    [[plab]][idx.leaf,] = q.sm.stress.pa    [! q.leaf.empty]
         qpatch$leaf.temp    [[plab]][idx.leaf,] = q.leaf.temp.pa    [! q.leaf.empty]
         qpatch$leaf.water   [[plab]][idx.leaf,] = q.leaf.water.pa   [! q.leaf.empty]
         qpatch$leaf.vpd     [[plab]][idx.leaf,] = q.leaf.vpd.pa     [! q.leaf.empty]
         qpatch$leaf.gpp     [[plab]][idx.leaf,] = q.leaf.gpp.pa     [! q.leaf.empty]
         qpatch$leaf.gsw     [[plab]][idx.leaf,] = q.leaf.gsw.pa     [! q.leaf.empty]
         qpatch$leaf.par     [[plab]][idx.leaf,] = q.leaf.par.pa     [! q.leaf.empty]
         qpatch$leaf.par.beam[[plab]][idx.leaf,] = q.leaf.par.beam.pa[! q.leaf.empty]
         qpatch$leaf.par.diff[[plab]][idx.leaf,] = q.leaf.par.diff.pa[! q.leaf.empty]
         qpatch$assim.light  [[plab]][idx.leaf,] = q.assim.light.pa  [! q.leaf.empty]
         qpatch$assim.rubp   [[plab]][idx.leaf,] = q.assim.rubp.pa   [! q.leaf.empty]
         qpatch$assim.co2    [[plab]][idx.leaf,] = q.assim.co2.pa    [! q.leaf.empty]
         qpatch$wood.temp    [[plab]][idx.wood,] = q.wood.temp.pa    [! q.wood.empty]
         qpatch$gpp          [[plab]][idx     ,] = q.gpp.pa
         qpatch$npp          [[plab]][idx     ,] = q.npp.pa
         qpatch$plant.resp   [[plab]][idx     ,] = q.plant.resp.pa
         qpatch$hflxlc       [[plab]][idx     ,] = q.hflxlc.pa
         qpatch$hflxwc       [[plab]][idx     ,] = q.hflxwc.pa
         qpatch$wflxlc       [[plab]][idx     ,] = q.wflxlc.pa
         qpatch$wflxwc       [[plab]][idx     ,] = q.wflxwc.pa
         qpatch$transp       [[plab]][idx     ,] = q.transp.pa
         #------ Soil respiration mixes cohort (root) and patch (hetetrophic). ------------#
         qpatch$soil.resp[[plab]][idx,] = qpatch$soil.resp[[plab]][idx,] + q.root.resp.pa
         #---------------------------------------------------------------------------------#
      }#end if (any(ncohorts) > 0)
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Ecosystem respiration, which is a combination of plant respiration (cohort-    #
      # -based) and heterotrophic respiration (patch-based).                               #
      #------------------------------------------------------------------------------------#
      patch$reco [[plab]] = patch$plant.resp [[plab]] + patch$het.resp [[plab]]
      qpatch$reco[[plab]] = qpatch$plant.resp[[plab]] + qpatch$het.resp[[plab]]
      #------------------------------------------------------------------------------------#

      #====================================================================================#
      #====================================================================================#
      #====================================================================================#
      #====================================================================================#










      #------------------------------------------------------------------------------------#
      #     The following two variables are used to scale "intensive" properties           #
      # (whatever/plant) to "extensive" (whatever/m2).  Sometimes they may be used to      #
      # build weighted averages.                                                           #
      #------------------------------------------------------------------------------------#
      w.nplant  = nplantconow  * areaconow
      w.lai     = laiconow     * areaconow
      w.wai     = waiconow     * areaconow
      w.tai     = taiconow     * areaconow
      w.biomass = biomassconow * w.nplant
      w.balive  = baliveconow  * w.nplant
      w.basarea = baconow      * w.nplant
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
            arealu.i          = 1. / sum(areapa[selpa])
            lu$lai      [m,l] = sum( w.lai   [sel]                      )
            lu$ba       [m,l] = sum( w.nplant[sel] * baconow      [sel] )
            lu$agb      [m,l] = sum( w.nplant[sel] * agbconow     [sel] )
            lu$biomass  [m,l] = sum( w.nplant[sel] * biomassconow [sel] )
            lu$gpp      [m,l] = sum( w.nplant[sel] * gppconow     [sel] )
            lu$npp      [m,l] = sum( w.nplant[sel] * nppconow     [sel] )
            lu$f.lai    [m,l] = sum( w.lai   [sel]                      ) * arealu.i
            lu$f.ba     [m,l] = sum( w.nplant[sel] * baconow      [sel] ) * arealu.i
            lu$f.agb    [m,l] = sum( w.nplant[sel] * agbconow     [sel] ) * arealu.i
            lu$f.biomass[m,l] = sum( w.nplant[sel] * biomassconow [sel] ) * arealu.i
            lu$f.gpp    [m,l] = sum( w.nplant[sel] * gppconow     [sel] ) * arealu.i
            lu$f.npp    [m,l] = sum( w.nplant[sel] * nppconow     [sel] ) * arealu.i
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
            sel.dbh.lmon  = rep(FALSE,times=length(dbhfac.lmon))

            #----- Define the minimum DBH. ------------------------------------------------#
            dbhminconow   = rep(Inf,times=length(pftconow))
            #------------------------------------------------------------------------------#
         }else{
            sel.dbh       = dbhfac      == d | d == (ndbh+1)
            sel.dbh.1ago  = dbhfac.1ago == d | d == (ndbh+1)
            sel.dbh.lmon  = dbhfac.lmon == d | d == (ndbh+1)

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
               szpft$lai          [m,d,p] = sum( laiconow          [sel]
                                               * areaconow         [sel]
                                               , na.rm = TRUE
                                               )#end sum
               szpft$wai          [m,d,p] = sum( waiconow          [sel]
                                               * areaconow         [sel]
                                               , na.rm = TRUE
                                               )#end sum
               szpft$tai          [m,d,p] = sum( taiconow          [sel]
                                               * areaconow         [sel]
                                               , na.rm = TRUE
                                               )#end sum
               szpft$nplant       [m,d,p] = sum( nplantconow       [sel]
                                               * areaconow         [sel]
                                               , na.rm = TRUE
                                               )#end sum
               szpft$demand       [m,d,p] = sum( demandconow       [sel]
                                               * areaconow         [sel]
                                               , na.rm = TRUE
                                               )#end sum
               szpft$supply       [m,d,p] = sum( supplyconow       [sel]
                                               * areaconow         [sel]
                                               , na.rm = TRUE
                                               )#end sum
               #----- Intensive properties, use nplant to make them extensive. ------------#
               szpft$agb          [m,d,p] = sum( w.nplant          [sel]
                                               * agbconow          [sel]
                                               , na.rm = TRUE
                                               )#end sum
               szpft$biomass      [m,d,p] = sum( w.nplant          [sel]
                                               * biomassconow      [sel]
                                               , na.rm = TRUE
                                               )#end sum
               szpft$ba           [m,d,p] = sum( w.nplant          [sel]
                                               * baconow           [sel]
                                               , na.rm = TRUE
                                               )#end sum
               szpft$gpp          [m,d,p] = sum( w.nplant          [sel]
                                               * gppconow          [sel]
                                               , na.rm = TRUE
                                               )#end sum
               szpft$npp          [m,d,p] = sum( w.nplant          [sel]
                                               * nppconow          [sel]
                                               , na.rm = TRUE
                                               )#end sum
               szpft$mco          [m,d,p] = sum( w.nplant          [sel]
                                               * mcostconow        [sel]
                                               , na.rm = TRUE
                                               )#end sum
               szpft$dcbadt       [m,d,p] = sum( w.nplant          [sel]
                                               * dcbadtconow       [sel]
                                               , na.rm = TRUE
                                               )#end sum
               szpft$cba          [m,d,p] = sum( w.nplant          [sel]
                                               * cbaconow          [sel]
                                               , na.rm = TRUE
                                               )#end sum
               szpft$cbamax       [m,d,p] = sum( w.nplant          [sel]
                                               * cbamaxconow       [sel]
                                               , na.rm = TRUE
                                               )#end sum
               szpft$cbalight     [m,d,p] = sum( w.nplant          [sel]
                                               * cbalightconow     [sel]
                                               , na.rm = TRUE
                                               )#end sum
               szpft$cbamoist     [m,d,p] = sum( w.nplant          [sel]
                                               * cbamoistconow     [sel]
                                               , na.rm = TRUE
                                               )#end sum
               szpft$ldrop        [m,d,p] = sum( w.nplant          [sel]
                                               * ldropconow        [sel]
                                               , na.rm = TRUE
                                               )#end sum
               szpft$leaf.resp    [m,d,p] = sum( w.nplant          [sel]
                                               * leaf.respconow    [sel]
                                               , na.rm = TRUE
                                               )#end sum
               szpft$stem.resp    [m,d,p] = sum( w.nplant          [sel]
                                               * stem.respconow    [sel]
                                               , na.rm = TRUE
                                               )#end sum
               szpft$root.resp    [m,d,p] = sum( w.nplant          [sel]
                                               * root.respconow    [sel]
                                               , na.rm = TRUE
                                               )#end sum
               szpft$froot.resp   [m,d,p] = sum( w.nplant          [sel]
                                               * froot.respconow   [sel]
                                               , na.rm = TRUE
                                               )#end sum
               szpft$croot.resp   [m,d,p] = sum( w.nplant          [sel]
                                               * croot.respconow   [sel]
                                               , na.rm = TRUE
                                               )#end sum
               szpft$growth.resp  [m,d,p] = sum( w.nplant          [sel]
                                               * growth.respconow  [sel]
                                               , na.rm = TRUE
                                               )#end sum
               szpft$storage.resp [m,d,p] = sum( w.nplant          [sel]
                                               * storage.respconow [sel]
                                               , na.rm = TRUE
                                               )#end sum
               szpft$plant.resp   [m,d,p] = sum( w.nplant          [sel]
                                               * plant.respconow   [sel]
                                               , na.rm = TRUE
                                               )#end sum
               szpft$bdead        [m,d,p] = sum( w.nplant          [sel]
                                               * bdeadconow        [sel]
                                               , na.rm = TRUE
                                               )#end sum
               szpft$balive       [m,d,p] = sum( w.nplant          [sel]
                                               * baliveconow       [sel]
                                               , na.rm = TRUE
                                               )#end sum
               szpft$bleaf        [m,d,p] = sum( w.nplant          [sel]
                                               * bleafconow        [sel]
                                               , na.rm = TRUE
                                               )#end sum
               szpft$bstem        [m,d,p] = sum( w.nplant          [sel]
                                               * bstemconow        [sel]
                                               , na.rm = TRUE
                                               )#end sum
               szpft$broot        [m,d,p] = sum( w.nplant          [sel]
                                               * brootconow        [sel]
                                               , na.rm = TRUE
                                               )#end sum
               szpft$bfroot       [m,d,p] = sum( w.nplant          [sel]
                                               * bfrootconow       [sel]
                                               , na.rm = TRUE
                                               )#end sum
               szpft$bcroot       [m,d,p] = sum( w.nplant          [sel]
                                               * bcrootconow       [sel]
                                               , na.rm = TRUE
                                               )#end sum
               szpft$bsapwood     [m,d,p] = sum( w.nplant          [sel]
                                               * bsapwoodconow     [sel]
                                               , na.rm = TRUE
                                               )#end sum
               szpft$bstorage     [m,d,p] = sum( w.nplant          [sel]
                                               * bstorageconow     [sel]
                                               , na.rm = TRUE
                                               )#end sum
               szpft$bseeds       [m,d,p] = sum( w.nplant          [sel]
                                               * bseedsconow       [sel]
                                               , na.rm = TRUE
                                               )#end sum
               szpft$hflxlc       [m,d,p] = sum( w.nplant          [sel]
                                               * hflxlcconow       [sel]
                                               , na.rm = TRUE
                                               )#end sum
               szpft$wflxlc       [m,d,p] = sum( w.nplant          [sel]
                                               * wflxlcconow       [sel]
                                               , na.rm = TRUE
                                               )#end sum
               szpft$transp       [m,d,p] = sum( w.nplant          [sel]
                                               * transpconow       [sel]
                                               , na.rm = TRUE
                                               )#end sum
               #----- Extensive properties that must be scaled by LAI, not nplant. --------#
               szpft$par.leaf     [m,d,p] = sum( w.lai             [sel]
                                               * leaf.parconow     [sel]
                                               , na.rm = TRUE
                                               )#end sum
               szpft$par.leaf.beam[m,d,p] = sum( w.lai             [sel]
                                               * leaf.par.beamconow[sel]
                                               , na.rm = TRUE
                                               )#end sum
               szpft$par.leaf.diff[m,d,p] = sum( w.lai             [sel]
                                               * leaf.par.diffconow[sel]
                                               , na.rm = TRUE
                                               )#end sum
               #---------------------------------------------------------------------------#



               #----- Leaf/wood intensive properties , weighted means using LAI/WAI. ------#
               szpft$sm.stress    [m,d,p] = weighted.mean( x     = sm.stressconow    [sel]
                                                         , w     = w.lai             [sel]
                                                         , na.rm = TRUE
                                                         )#end weighted.mean
               szpft$phap.sms     [m,d,p] = weighted.mean( x     = phap.smsconow     [sel]
                                                         , w     = w.lai             [sel]
                                                         , na.rm = TRUE
                                                         )#end weighted.mean
               szpft$leaf.par     [m,d,p] = weighted.mean( x     = leaf.parconow     [sel]
                                                         , w     = w.lai             [sel]
                                                         , na.rm = TRUE
                                                         )#end weighted.mean
               szpft$leaf.par.beam[m,d,p] = weighted.mean( x     = leaf.par.beamconow[sel]
                                                         , w     = w.lai             [sel]
                                                         , na.rm = TRUE
                                                         )#end weighted.mean
               szpft$leaf.par.diff[m,d,p] = weighted.mean( x     = leaf.par.diffconow[sel]
                                                         , w     = w.lai             [sel]
                                                         , na.rm = TRUE
                                                         )#end weighted.mean
               szpft$phap.lpar    [m,d,p] = weighted.mean( x     = phap.lparconow    [sel]
                                                         , w     = w.lai             [sel]
                                                         , na.rm = TRUE
                                                         )#end weighted.mean
               szpft$leaf.rshort  [m,d,p] = weighted.mean( x     = leaf.rshortconow  [sel]
                                                         , w     = w.lai             [sel]
                                                         , na.rm = TRUE
                                                         )#end weighted.mean
               szpft$leaf.rlong   [m,d,p] = weighted.mean( x     = leaf.rlongconow   [sel]
                                                         , w     = w.lai             [sel]
                                                         , na.rm = TRUE
                                                         )#end weighted.mean
               szpft$leaf.gpp     [m,d,p] = weighted.mean( x     = leaf.gppconow     [sel]
                                                         , w     = w.lai             [sel]
                                                         , na.rm = TRUE
                                                         )#end weighted.mean
               szpft$leaf.temp    [m,d,p] = weighted.mean( x     = leaf.tempconow    [sel]
                                                         , w     = w.lai             [sel]
                                                         , na.rm = TRUE
                                                         )#end weighted.mean
               szpft$phap.ltemp   [m,d,p] = weighted.mean( x     = phap.ltempconow   [sel]
                                                         , w     = w.lai             [sel]
                                                         , na.rm = TRUE
                                                         )#end weighted.mean
               szpft$leaf.water   [m,d,p] = weighted.mean( x     = leaf.waterconow   [sel]
                                                         , w     = w.lai             [sel]
                                                         , na.rm = TRUE
                                                         )#end weighted.mean
               szpft$phap.lwater  [m,d,p] = weighted.mean( x     = phap.lwaterconow  [sel]
                                                         , w     = w.lai             [sel]
                                                         , na.rm = TRUE
                                                         )#end weighted.mean
               szpft$wood.temp    [m,d,p] = weighted.mean( x     = wood.tempconow    [sel]
                                                         , w     = w.lai             [sel]
                                                         , na.rm = TRUE
                                                         )#end weighted.mean
               szpft$leaf.vpd     [m,d,p] = weighted.mean( x     = leaf.vpdconow     [sel]
                                                         , w     = w.lai             [sel]
                                                         , na.rm = TRUE
                                                         )#end weighted.mean
               szpft$phap.lvpd    [m,d,p] = weighted.mean( x     = phap.lvpdconow    [sel]
                                                         , w     = w.lai             [sel]
                                                         , na.rm = TRUE
                                                         )#end weighted.mean
               szpft$i.transp     [m,d,p] = weighted.mean( x     = i.transpconow     [sel]
                                                         , w     = w.lai             [sel]
                                                         , na.rm = TRUE
                                                         )#end weighted.mean
               szpft$i.wflxlc     [m,d,p] = weighted.mean( x     = i.wflxlcconow     [sel]
                                                         , w     = w.lai             [sel]
                                                         , na.rm = TRUE
                                                         )#end weighted.mean
               szpft$i.hflxlc     [m,d,p] = weighted.mean( x     = i.hflxlcconow     [sel]
                                                         , w     = w.lai             [sel]
                                                         , na.rm = TRUE
                                                         )#end weighted.mean
               szpft$leaf.gbw     [m,d,p] = weighted.mean( x     = leaf.gbwconow     [sel]
                                                         , w     = w.lai             [sel]
                                                         , na.rm = TRUE
                                                         )#end weighted.mean
               szpft$phap.lgbw    [m,d,p] = weighted.mean( x     = phap.lgbwconow    [sel]
                                                         , w     = w.lai             [sel]
                                                         , na.rm = TRUE
                                                         )#end weighted.mean
               szpft$leaf.gsw     [m,d,p] = weighted.mean( x     = leaf.gswconow     [sel]
                                                         , w     = w.lai             [sel]
                                                         , na.rm = TRUE
                                                         )#end weighted.mean
               szpft$phap.lgsw    [m,d,p] = weighted.mean( x     = phap.lgswconow    [sel]
                                                         , w     = w.lai             [sel]
                                                         , na.rm = TRUE
                                                         )#end weighted.mean
               szpft$assim.light  [m,d,p] = weighted.mean( x     = assim.lightconow  [sel]
                                                         , w     = w.lai             [sel]
                                                         , na.rm = TRUE
                                                         )#end weighted.mean
               szpft$assim.rubp   [m,d,p] = weighted.mean( x     = assim.rubpconow   [sel]
                                                         , w     = w.lai             [sel]
                                                         , na.rm = TRUE
                                                         )#end weighted.mean
               szpft$assim.co2    [m,d,p] = weighted.mean( x     = assim.co2conow    [sel]
                                                         , w     = w.lai             [sel]
                                                         , na.rm = TRUE
                                                         )#end weighted.mean
               szpft$wood.gbw     [m,d,p] = weighted.mean( x     = wood.gbwconow     [sel]
                                                         , w     = w.wai             [sel]
                                                         , na.rm = TRUE
                                                         )#end weighted.mean
               #---------------------------------------------------------------------------#


               #----- Individual-based properties, use weighted mean by Nplant. -----------#
               szpft$cbarel       [m,d,p] = weighted.mean( x     = cbarelconow      [sel]
                                                         , w     = w.nplant         [sel]
                                                         , na.rm = TRUE
                                                         )#end weighted.mean
               szpft$i.gpp        [m,d,p] = weighted.mean( x     = gppconow         [sel]
                                                         , w     = w.nplant         [sel]
                                                         , na.rm = TRUE
                                                         )#end weighted.mean
               szpft$i.npp        [m,d,p] = weighted.mean( x     = nppconow         [sel]
                                                         , w     = w.nplant         [sel]
                                                         , na.rm = TRUE
                                                         )#end weighted.mean
               szpft$i.plant.resp [m,d,p] = weighted.mean( x     = plant.respconow  [sel]
                                                         , w     = w.nplant         [sel]
                                                         , na.rm = TRUE
                                                         )#end weighted.mean
               szpft$i.mco        [m,d,p] = weighted.mean( x     = mcostconow       [sel]
                                                         , w     = w.nplant         [sel]
                                                         , na.rm = TRUE
                                                         )#end weighted.mean
               szpft$i.cba        [m,d,p] = weighted.mean( x     = cbaconow         [sel]
                                                         , w     = w.nplant         [sel]
                                                         , na.rm = TRUE
                                                         )#end weighted.mean
               szpft$i.cbamax     [m,d,p] = weighted.mean( x     = cbamaxconow      [sel]
                                                         , w     = w.nplant         [sel]
                                                         , na.rm = TRUE
                                                         )#end weighted.mean
               szpft$i.cbalight   [m,d,p] = weighted.mean( x     = cbalightconow    [sel]
                                                         , w     = w.nplant         [sel]
                                                         , na.rm = TRUE
                                                         )#end weighted.mean
               szpft$i.cbamoist   [m,d,p] = weighted.mean( x     = cbamoistconow    [sel]
                                                         , w     = w.nplant         [sel]
                                                         , na.rm = TRUE
                                                         )#end weighted.mean
               #---------------------------------------------------------------------------#


               #----- Wood density: averaged by basal area. -------------------------------#
               szpft$wood.dens    [m,d,p] = weighted.mean( x     = wood.densconow  [sel]
                                                         , w     = w.basarea       [sel]
                                                         , na.rm = TRUE
                                                         )#end weighted.mean
               #---------------------------------------------------------------------------#
            }#end if
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            # Use efficiency: use the bulk value to reduce biases towards large numbers.   #
            #                 Also, if none of the trees exist, or if the denominator      #
            #                 would be too close to zero, we don't calculate: instead we   #
            #                 make them NA.  Also, we use the past year data to reduce     #
            #                 noise.                                                       #
            #------------------------------------------------------------------------------#
            last.12              = seq (from=max(m-11,1),to=m,by=1)
            #----- Last year average, use na.rm=FALSE to skip extinctions. ----------------#
            last.1yr.transp      = mean(szpft$transp[last.12,d,p], na.rm=FALSE) * yr.day
            last.1yr.rain        = sum (emean$rain  [last.12]    , na.rm=FALSE)
            last.1yr.gpp         = mean(szpft$gpp   [last.12,d,p], na.rm=FALSE)
            last.1yr.npp         = mean(szpft$npp   [last.12,d,p], na.rm=FALSE)
            last.1yr.dcbadt      = mean(szpft$dcbadt[last.12,d,p], na.rm=FALSE)
            if ( d == ndbh+1 & p == npft+1){
                last.1yr.et      = mean(emean$et[last.12], na.rm=FALSE) * yr.day
            }else{
                last.1yr.et      = mean( szpft$transp[last.12,d,p] 
                                       + szpft$wflxlc[last.12,d,p], na.rm=FALSE ) * yr.day
            }#end if
            #----- Make sure that plants were doing something. ----------------------------#
            tfine                = is.finite(last.1yr.transp) & last.1yr.transp >= 1.0
            rfine                = is.finite(last.1yr.rain)   & last.1yr.rain   >= 1.0
            #----- Find use efficiencies only if things are fine. -------------------------#
            if (tfine){
               szpft$wue [m,d,p] = 1000. * last.1yr.npp    / last.1yr.transp
               szpft$etue[m,d,p] = 1000. * last.1yr.npp    / last.1yr.et
               szpft$cue [m,d,p] =         last.1yr.npp    / last.1yr.gpp
               szpft$ecue[m,d,p] =         last.1yr.dcbadt / last.1yr.gpp
            }#end if
            if (rfine){
               szpft$rue [m,d,p] = 1000. * last.1yr.npp / last.1yr.rain
            }#end if
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            # Fractional biomass: use the total variable and divide by the total biomass   #
            #                     so it gives a full community fraction.  Like in the UE   #
            #                     case, force values to be NA if no cohort matches this    #
            #                     class, or if it didn't have much living biomass.         #
            #------------------------------------------------------------------------------#
            balive.szpft              = szpft$balive[m,d,p]
            balive.szpft              = ifelse(balive.szpft > 1.e-7,balive.szpft,NA)
            szpft$f.gpp       [m,d,p] = 100. * szpft$gpp       [m,d,p] / balive.szpft
            szpft$f.plant.resp[m,d,p] = 100. * szpft$plant.resp[m,d,p] / balive.szpft
            szpft$f.npp       [m,d,p] = 100. * szpft$npp       [m,d,p] / balive.szpft
            szpft$f.mco       [m,d,p] = 100. * szpft$mco       [m,d,p] / balive.szpft
            szpft$f.dcbadt    [m,d,p] = 100. * szpft$dcbadt    [m,d,p] / balive.szpft
            szpft$f.cba       [m,d,p] =        szpft$cba       [m,d,p] / balive.szpft
            szpft$f.bstorage  [m,d,p] =        szpft$bstorage  [m,d,p] / balive.szpft
            szpft$f.bleaf     [m,d,p] =        szpft$bleaf     [m,d,p] / balive.szpft
            szpft$f.bstem     [m,d,p] =        szpft$bstem     [m,d,p] / balive.szpft
            szpft$f.broot     [m,d,p] =        szpft$broot     [m,d,p] / balive.szpft
            szpft$f.bseeds    [m,d,p] =        szpft$bseeds    [m,d,p] / balive.szpft
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #     Assimilation ratio: use the averaged values.                             #
            #------------------------------------------------------------------------------#
            szpft$assim.ratio [m,d,p] = ( szpft$assim.light[m,d,p]
                                        / max( 1e-6, min( szpft$assim.rubp[m,d,p]
                                                        , szpft$assim.co2 [m,d,p] ) ) )
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #    For mortality and growth, we keep deleting the tiny guys because they     #
            # skew the rates quite significantly.                                          #
            #------------------------------------------------------------------------------#
            sel = sel.pft & sel.dbh & dbhconow >= dbhminconow
            if (any(sel)){
               #----- Growth rates are weighted by population. ----------------------------#
               dbh.growth = - 100. * log( weighted.mean( x = exp(-growthconow    [sel])
                                                       , w = w.nplant            [sel]
                                                           * dbhconow            [sel]
                                                       )#end weighted.mean
                                        )#end log
               agb.growth = - 100. * log( weighted.mean( x = exp(-agb.growthconow[sel])
                                                       , w = w.nplant            [sel]
                                                           * agbconow            [sel]
                                                       )#end weighted.mean
                                        )#end log
               bsa.growth = - 100. * log( weighted.mean( x = exp(-bsa.growthconow[sel])
                                                       , w = w.nplant            [sel]
                                                           * baconow             [sel]
                                                       )#end weighted.mean
                                        )#end log
               acc.growth = sum( w.nplant[sel]
                               * agbconow[sel] * (1.-exp(-agb.growthconow[sel]))
                               )#end sum
               szpft$growth     [m,d,p] = dbh.growth
               szpft$agb.growth [m,d,p] = agb.growth
               szpft$acc.growth [m,d,p] = acc.growth
               szpft$bsa.growth [m,d,p] = bsa.growth
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
               survivor                 = sum( w.nplant[sel] * agbcolmon[sel])
               previous                 = sum( w.nplant[sel] * agbcolmon[sel]
                                             * exp(mortconow            [sel] ) )
               ncb.previous             = sum( w.nplant[sel] * agbcolmon[sel]
                                             * exp(ncbmortconow         [sel] ) )
               di.previous              = sum( w.nplant[sel] * agbcolmon[sel]
                                             * exp(dimortconow          [sel] ) )
               szpft$agb.mort   [m,d,p] = log( previous     / survivor )
               szpft$agb.ncbmort[m,d,p] = log( ncb.previous / survivor )
               szpft$agb.dimort [m,d,p] = log( di.previous  / survivor )



               #---------------------------------------------------------------------------#
               #      Find the total AGB and previous AGB if the only mortality was the    #
               # mortality we test.                                                        #
               #---------------------------------------------------------------------------#
               survivor                 = sum( w.nplant[sel] * agbcolmon[sel])
               previous                 = sum( w.nplant[sel] * agbcolmon[sel]
                                             * exp(mortconow            [sel] / 12.) )
               ncb.previous             = sum( w.nplant[sel] * agbcolmon[sel]
                                             * exp(ncbmortconow         [sel] / 12. ) )
               di.previous              = sum( w.nplant[sel] * agbcolmon[sel]
                                             * exp(dimortconow          [sel] / 12. ) )
               szpft$acc.mort   [m,d,p] = 12. * (previous     - survivor)
               szpft$acc.ncbmort[m,d,p] = 12. * (ncb.previous - survivor)
               szpft$acc.dimort [m,d,p] = 12. * (di.previous  - survivor)
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #      Find the total basal area and previous basal area if the only        #
               # mortality was the mortality we test.                                      #
               #---------------------------------------------------------------------------#
               survivor                 = sum( w.nplant[sel] * bacolmon[sel])
               previous                 = sum( w.nplant[sel] * bacolmon[sel]
                                             * exp(mortconow           [sel] ) )
               ncb.previous             = sum( w.nplant[sel] * bacolmon[sel]
                                             * exp(ncbmortconow        [sel] ) )
               di.previous              = sum( w.nplant[sel] * bacolmon[sel]
                                             * exp(dimortconow         [sel] ) )
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
            sel.elm = sel.pop & sel.dbh.lmon & dbhconow.lmon >= dbhminconow
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


               #----- Recruitment rate in terms of above-ground biomass. ------------------#
               population             = sum(w.nplant[sel.pop] * agbconow[sel.pop])
               established            = sum(w.nplant[sel.elm] * agbconow[sel.elm])
               szpft$acc.recr [m,d,p] = 12. * (population - established)
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



            #------------------------------------------------------------------------------#
            #     Change of basic properties.                                              #
            #------------------------------------------------------------------------------#
            if (m == 1){
               szpft$change     [m,d,p] = 0.
               szpft$agb.change [m,d,p] = 0.
               szpft$acc.change [m,d,p] = 0.
               szpft$bsa.change [m,d,p] = 0.
            }else{
               szpft$change     [m,d,p] = ( 12. * log( szpft$nplant[m  ,d,p]
                                                     / szpft$nplant[m-1,d,p] ) )
               szpft$agb.change [m,d,p] = ( 12. * log( szpft$agb   [m  ,d,p]
                                                     / szpft$agb   [m-1,d,p] ) )
               szpft$acc.change [m,d,p] = ( 12. *    ( szpft$agb   [m  ,d,p]
                                                     - szpft$agb   [m-1,d,p] ) )
               szpft$bsa.change [m,d,p] = ( 12. * log( szpft$ba    [m  ,d,p]
                                                     / szpft$ba    [m-1,d,p] ) )
            }#end if
            #------------------------------------------------------------------------------#
         }#end for PFT
         #---------------------------------------------------------------------------------#
      }#end for DBH
      #------------------------------------------------------------------------------------#




      #------------------------------------------------------------------------------------#
      #       Find the polygon-average depth and area.                                     #
      #------------------------------------------------------------------------------------#
      if (any(ncohorts > 0)){
         useconow      = as.numeric(opencanconow > tiny.num)
         xconow        = heightconow  * nplantconow * baconow * useconow
         wconow        = nplantconow  * baconow * useconow
         oconow        = opencanconow * useconow
         can.depth.idx = (    tapply(X = xconow, INDEX = ipaconow, FUN = sum, na.rm = TRUE)
                         /    tapply(X = wconow, INDEX = ipaconow, FUN = sum, na.rm = TRUE)
                         )#end can.depth.pa
         can.area.idx  = 1. - tapply(X = oconow, INDEX = ipaconow, FUN = min, na.rm = TRUE)
         
         can.depth.pa      = rep(NA,times=sum(npatches))
         can.area.pa       = rep(NA,times=sum(npatches))

         idx               = as.numeric(names(can.depth.idx))
         can.depth.pa[idx] = can.depth.idx
         can.area.pa [idx] = can.area.idx

         emean$can.depth [m] = sum(can.depth.pa * areapa)
         emean$can.area  [m] = sum(can.area.pa  * areapa)
      }else{
         emean$can.depth [m] = 0.
         emean$can.area  [m] = 0.
      }#end if
      #------------------------------------------------------------------------------------#




      #------------------------------------------------------------------------------------#
      #       Build the derived variables.                                                 #
      #------------------------------------------------------------------------------------#
      emean$leaf.resp       [m] = szpft$leaf.resp      [m,ndbh+1,npft+1]
      emean$stem.resp       [m] = szpft$stem.resp      [m,ndbh+1,npft+1]
      emean$root.resp       [m] = szpft$root.resp      [m,ndbh+1,npft+1]
      emean$froot.resp      [m] = szpft$froot.resp     [m,ndbh+1,npft+1]
      emean$croot.resp      [m] = szpft$croot.resp     [m,ndbh+1,npft+1]
      emean$soil.resp       [m] = emean$root.resp [m] + emean$het.resp[m]
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
      emean$f.mco           [m] = szpft$f.mco          [m,ndbh+1,npft+1]
      emean$f.cba           [m] = szpft$f.cba          [m,ndbh+1,npft+1]
      emean$f.bstorage      [m] = szpft$f.bstorage     [m,ndbh+1,npft+1]
      emean$f.bleaf         [m] = szpft$f.bleaf        [m,ndbh+1,npft+1]
      emean$f.bstem         [m] = szpft$f.bstem        [m,ndbh+1,npft+1]
      emean$f.broot         [m] = szpft$f.broot        [m,ndbh+1,npft+1]
      emean$f.bseeds        [m] = szpft$f.bseeds       [m,ndbh+1,npft+1]
      emean$f.dcbadt        [m] = szpft$f.dcbadt       [m,ndbh+1,npft+1]
      emean$leaf.par        [m] = szpft$leaf.par       [m,ndbh+1,npft+1]
      emean$leaf.par.beam   [m] = szpft$leaf.par.beam  [m,ndbh+1,npft+1]
      emean$leaf.par.diff   [m] = szpft$leaf.par.diff  [m,ndbh+1,npft+1]
      emean$leaf.rshort     [m] = szpft$leaf.rshort    [m,ndbh+1,npft+1]
      emean$leaf.rlong      [m] = szpft$leaf.rlong     [m,ndbh+1,npft+1]
      emean$leaf.gpp        [m] = szpft$leaf.gpp       [m,ndbh+1,npft+1]
      emean$transp          [m] = szpft$transp         [m,ndbh+1,npft+1]
      emean$wue             [m] = szpft$wue            [m,ndbh+1,npft+1]
      emean$npp             [m] = szpft$npp            [m,ndbh+1,npft+1]
      emean$dcbadt          [m] = szpft$dcbadt         [m,ndbh+1,npft+1]
      emean$rue             [m] = szpft$rue            [m,ndbh+1,npft+1]
      emean$etue            [m] = szpft$etue           [m,ndbh+1,npft+1]
      emean$cue             [m] = szpft$cue            [m,ndbh+1,npft+1]
      emean$ecue            [m] = szpft$ecue           [m,ndbh+1,npft+1]
      emean$agb.growth      [m] = szpft$agb.growth     [m,ndbh+1,npft+1]
      emean$agb.mort        [m] = szpft$agb.mort       [m,ndbh+1,npft+1]
      emean$agb.dimort      [m] = szpft$agb.dimort     [m,ndbh+1,npft+1]
      emean$agb.ncbmort     [m] = szpft$agb.ncbmort    [m,ndbh+1,npft+1]
      emean$agb.change      [m] = szpft$agb.change     [m,ndbh+1,npft+1]
      emean$acc.change      [m] = szpft$acc.change     [m,ndbh+1,npft+1]
      emean$acc.growth      [m] = szpft$acc.growth     [m,ndbh+1,npft+1]
      emean$acc.mort        [m] = szpft$acc.mort       [m,ndbh+1,npft+1]
      emean$acc.ncbmort     [m] = szpft$acc.ncbmort    [m,ndbh+1,npft+1]
      emean$acc.dimort      [m] = szpft$acc.dimort     [m,ndbh+1,npft+1]
      emean$acc.recr        [m] = szpft$acc.recr       [m,ndbh+1,npft+1]
      emean$wood.dens       [m] = szpft$wood.dens      [m,ndbh+1,npft+1]
      emean$phap.lpar       [m] = szpft$phap.lpar      [m,ndbh+1,npft+1]
      emean$phap.ltemp      [m] = szpft$phap.ltemp     [m,ndbh+1,npft+1]
      emean$phap.lvpd       [m] = szpft$phap.lvpd      [m,ndbh+1,npft+1]
      emean$phap.lwater     [m] = szpft$phap.lwater    [m,ndbh+1,npft+1]
      emean$phap.lgsw       [m] = szpft$phap.lgsw      [m,ndbh+1,npft+1]
      emean$phap.lgbw       [m] = szpft$phap.lgbw      [m,ndbh+1,npft+1]
      emean$phap.sms        [m] = szpft$phap.sms       [m,ndbh+1,npft+1]
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #      Convert leaf water to kg/m2leaf.                                              #
      #------------------------------------------------------------------------------------#
      emean$leaf.water  [m ] = emean$leaf.water[m ] / pmax(emean$lai[m],0.01)
      qmean$leaf.water  [m,] = qmean$leaf.water[m,] / pmax(emean$lai[m],0.01)
      #------------------------------------------------------------------------------------#





      #------------------------------------------------------------------------------------#
      #      Find the last-123 variables that are directly averaged/summed/maximised.      #
      #------------------------------------------------------------------------------------#
      last.12 = seq(from=max(m-11,1),to=m,by=1)
      last.24 = seq(from=max(m-23,1),to=m,by=1)
      last.36 = seq(from=max(m-35,1),to=m,by=1)
      #----- Gross primary productivity. --------------------------------------------------#
      emean$last.1yr.gpp     [m] = mean(emean$gpp           [last.12],na.rm=TRUE)
      emean$last.2yr.gpp     [m] = mean(emean$gpp           [last.24],na.rm=TRUE)
      emean$last.3yr.gpp     [m] = mean(emean$gpp           [last.36],na.rm=TRUE)
      #----- Plant respiration. -----------------------------------------------------------#
      emean$last.1yr.plresp  [m] = mean(emean$plant.resp    [last.12],na.rm=TRUE)
      emean$last.2yr.plresp  [m] = mean(emean$plant.resp    [last.24],na.rm=TRUE)
      emean$last.3yr.plresp  [m] = mean(emean$plant.resp    [last.36],na.rm=TRUE)
      #----- Carbon balance. --------------------------------------------------------------#
      emean$last.1yr.cba     [m] = mean(emean$cba           [last.12],na.rm=TRUE)
      emean$last.2yr.cba     [m] = mean(emean$cba           [last.24],na.rm=TRUE)
      emean$last.3yr.cba     [m] = mean(emean$cba           [last.36],na.rm=TRUE)
      #----- Evaporation. -----------------------------------------------------------------#
      emean$last.1yr.evap    [m] = mean(emean$evap          [last.12],na.rm=TRUE)
      emean$last.2yr.evap    [m] = mean(emean$evap          [last.24],na.rm=TRUE)
      emean$last.3yr.evap    [m] = mean(emean$evap          [last.36],na.rm=TRUE)
      #----- Net primary production. ------------------------------------------------------#
      emean$last.1yr.npp     [m] = mean(emean$npp           [last.12],na.rm=TRUE)
      emean$last.2yr.npp     [m] = mean(emean$npp           [last.24],na.rm=TRUE)
      emean$last.3yr.npp     [m] = mean(emean$npp           [last.36],na.rm=TRUE)
      #----- Change in carbon balance. ----------------------------------------------------#
      emean$last.1yr.dcbadt  [m] = mean(emean$dcbadt        [last.12],na.rm=TRUE)
      emean$last.2yr.dcbadt  [m] = mean(emean$dcbadt        [last.24],na.rm=TRUE)
      emean$last.3yr.dcbadt  [m] = mean(emean$dcbadt        [last.36],na.rm=TRUE)
      #----- Evapotranspiration. ----------------------------------------------------------#
      emean$last.1yr.et      [m] = mean(emean$et            [last.12],na.rm=TRUE)
      emean$last.2yr.et      [m] = mean(emean$et            [last.24],na.rm=TRUE)
      emean$last.3yr.et      [m] = mean(emean$et            [last.36],na.rm=TRUE)
      #----- Transpiration. ---------------------------------------------------------------#
      emean$last.1yr.transp  [m] = mean(emean$transp        [last.12],na.rm=TRUE)
      emean$last.2yr.transp  [m] = mean(emean$transp        [last.24],na.rm=TRUE)
      emean$last.3yr.transp  [m] = mean(emean$transp        [last.36],na.rm=TRUE)
      #----- Rainfall. --------------------------------------------------------------------#
      emean$last.1yr.rain    [m] = mean(emean$rain          [last.12],na.rm=TRUE) * 12.
      emean$last.2yr.rain    [m] = mean(emean$rain          [last.24],na.rm=TRUE) * 12.
      emean$last.3yr.rain    [m] = mean(emean$rain          [last.36],na.rm=TRUE) * 12.
      #----- Shortwave radiation. ---------------------------------------------------------#
      emean$last.1yr.rshort  [m] = mean(emean$rshort        [last.12],na.rm=TRUE)
      emean$last.2yr.rshort  [m] = mean(emean$rshort        [last.24],na.rm=TRUE)
      emean$last.3yr.rshort  [m] = mean(emean$rshort        [last.36],na.rm=TRUE)
      #----- Soil matric potential. -------------------------------------------------------#
      emean$last.1yr.smpot   [m] = mean(emean$smpot         [last.12],na.rm=TRUE)
      emean$last.2yr.smpot   [m] = mean(emean$smpot         [last.24],na.rm=TRUE)
      emean$last.3yr.smpot   [m] = mean(emean$smpot         [last.36],na.rm=TRUE)
      #----- Maximum water deficit of the past period. ------------------------------------#
      emean$last.1yr.mwd     [m] = max (emean$water.deficit [last.12],na.rm=TRUE)
      emean$last.2yr.mwd     [m] = max (emean$water.deficit [last.24],na.rm=TRUE)
      emean$last.3yr.mwd     [m] = max (emean$water.deficit [last.36],na.rm=TRUE)
      #----- Growth rate. -----------------------------------------------------------------#
      emean$last.1yr.growth  [m] = mean(emean$agb.growth    [last.12],na.rm=TRUE)
      emean$last.2yr.growth  [m] = mean(emean$agb.growth    [last.24],na.rm=TRUE)
      emean$last.3yr.growth  [m] = mean(emean$agb.growth    [last.36],na.rm=TRUE)
      #----- Mortality rate. --------------------------------------------------------------#
      emean$last.1yr.mort    [m] = mean(emean$agb.mort      [last.12],na.rm=TRUE)
      emean$last.2yr.mort    [m] = mean(emean$agb.mort      [last.24],na.rm=TRUE)
      emean$last.3yr.mort    [m] = mean(emean$agb.mort      [last.36],na.rm=TRUE)
      #----- Density-independent mortality rate. ------------------------------------------#
      emean$last.1yr.dimort  [m] = mean(emean$agb.dimort    [last.12],na.rm=TRUE)
      emean$last.2yr.dimort  [m] = mean(emean$agb.dimort    [last.24],na.rm=TRUE)
      emean$last.3yr.dimort  [m] = mean(emean$agb.dimort    [last.36],na.rm=TRUE)
      #----- Density-dependent mortality rate. --------------------------------------------#
      emean$last.1yr.ncbmort [m] = mean(emean$agb.ncbmort   [last.12],na.rm=TRUE)
      emean$last.2yr.ncbmort [m] = mean(emean$agb.ncbmort   [last.24],na.rm=TRUE)
      emean$last.3yr.ncbmort [m] = mean(emean$agb.ncbmort   [last.36],na.rm=TRUE)
      #----- AGB change. ------------------------------------------------------------------#
      emean$last.1yr.change  [m] = mean(emean$agb.change    [last.12],na.rm=TRUE)
      emean$last.2yr.change  [m] = mean(emean$agb.change    [last.24],na.rm=TRUE)
      emean$last.3yr.change  [m] = mean(emean$agb.change    [last.36],na.rm=TRUE)
      #----- The following variables depend on whether to use PhAP or 24 hours. -----------#
      if (iint.photo == 0){
         #----- Leaf absorbed PAR. --------------------------------------------------------#
         emean$last.1yr.lpar    [m] = mean(emean$leaf.par      [last.12],na.rm=TRUE)
         emean$last.2yr.lpar    [m] = mean(emean$leaf.par      [last.24],na.rm=TRUE)
         emean$last.3yr.lpar    [m] = mean(emean$leaf.par      [last.36],na.rm=TRUE)
         #----- Leaf absorbed PAR. --------------------------------------------------------#
         emean$last.1yr.lgpp    [m] = mean(emean$leaf.gpp      [last.12],na.rm=TRUE)
         emean$last.2yr.lgpp    [m] = mean(emean$leaf.gpp      [last.24],na.rm=TRUE)
         emean$last.3yr.lgpp    [m] = mean(emean$leaf.gpp      [last.36],na.rm=TRUE)
         #----- Leaf temperature. ---------------------------------------------------------#
         emean$last.1yr.ltemp   [m] = mean(emean$leaf.temp     [last.12],na.rm=TRUE)
         emean$last.2yr.ltemp   [m] = mean(emean$leaf.temp     [last.24],na.rm=TRUE)
         emean$last.3yr.ltemp   [m] = mean(emean$leaf.temp     [last.36],na.rm=TRUE)
         #----- Leaf vapour pressure deficit. ---------------------------------------------#
         emean$last.1yr.lvpd    [m] = mean(emean$leaf.vpd      [last.12],na.rm=TRUE)
         emean$last.2yr.lvpd    [m] = mean(emean$leaf.vpd      [last.24],na.rm=TRUE)
         emean$last.3yr.lvpd    [m] = mean(emean$leaf.vpd      [last.36],na.rm=TRUE)
         #----- Leaf water. ---------------------------------------------------------------#
         emean$last.1yr.lwater  [m] = mean(emean$leaf.water    [last.12],na.rm=TRUE)
         emean$last.2yr.lwater  [m] = mean(emean$leaf.water    [last.24],na.rm=TRUE)
         emean$last.3yr.lwater  [m] = mean(emean$leaf.water    [last.36],na.rm=TRUE)
         #----- Stomatal conductance. -----------------------------------------------------#
         emean$last.1yr.lgsw    [m] = mean(emean$leaf.gsw      [last.12],na.rm=TRUE)
         emean$last.2yr.lgsw    [m] = mean(emean$leaf.gsw      [last.24],na.rm=TRUE)
         emean$last.3yr.lgsw    [m] = mean(emean$leaf.gsw      [last.36],na.rm=TRUE)
         #----- Soil moisture stress. -----------------------------------------------------#
         emean$last.1yr.sms     [m] = mean(emean$sm.stress     [last.12],na.rm=TRUE)
         emean$last.2yr.sms     [m] = mean(emean$sm.stress     [last.24],na.rm=TRUE)
         emean$last.3yr.sms     [m] = mean(emean$sm.stress     [last.36],na.rm=TRUE)
         #---------------------------------------------------------------------------------#
      }else{
         #----- Leaf absorbed PAR. --------------------------------------------------------#
         emean$last.1yr.lpar    [m] = mean(emean$phap.lpar     [last.12],na.rm=TRUE)
         emean$last.2yr.lpar    [m] = mean(emean$phap.lpar     [last.24],na.rm=TRUE)
         emean$last.3yr.lpar    [m] = mean(emean$phap.lpar     [last.36],na.rm=TRUE)
         #----- Leaf temperature. ---------------------------------------------------------#
         emean$last.1yr.ltemp   [m] = mean(emean$phap.ltemp    [last.12],na.rm=TRUE)
         emean$last.2yr.ltemp   [m] = mean(emean$phap.ltemp    [last.24],na.rm=TRUE)
         emean$last.3yr.ltemp   [m] = mean(emean$phap.ltemp    [last.36],na.rm=TRUE)
         #----- Leaf vapour pressure deficit. ---------------------------------------------#
         emean$last.1yr.lvpd    [m] = mean(emean$phap.lvpd     [last.12],na.rm=TRUE)
         emean$last.2yr.lvpd    [m] = mean(emean$phap.lvpd     [last.24],na.rm=TRUE)
         emean$last.3yr.lvpd    [m] = mean(emean$phap.lvpd     [last.36],na.rm=TRUE)
         #----- Leaf water. ---------------------------------------------------------------#
         emean$last.1yr.lwater  [m] = mean(emean$phap.lwater   [last.12],na.rm=TRUE)
         emean$last.2yr.lwater  [m] = mean(emean$phap.lwater   [last.24],na.rm=TRUE)
         emean$last.3yr.lwater  [m] = mean(emean$phap.lwater   [last.36],na.rm=TRUE)
         #----- Stomatal conductance. -----------------------------------------------------#
         emean$last.1yr.lgsw    [m] = mean(emean$phap.lgsw     [last.12],na.rm=TRUE)
         emean$last.2yr.lgsw    [m] = mean(emean$phap.lgsw     [last.24],na.rm=TRUE)
         emean$last.3yr.lgsw    [m] = mean(emean$phap.lgsw     [last.36],na.rm=TRUE)
         #----- Soil moisture stress. -----------------------------------------------------#
         emean$last.1yr.sms     [m] = mean(emean$phap.sms      [last.12],na.rm=TRUE)
         emean$last.2yr.sms     [m] = mean(emean$phap.sms      [last.24],na.rm=TRUE)
         emean$last.3yr.sms     [m] = mean(emean$phap.sms      [last.36],na.rm=TRUE)
         #---------------------------------------------------------------------------------#
      }#end if
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Find derived rates using the full period and the community averages, to avoid  #
      # biases towards large numbers.                                                      #
      #------------------------------------------------------------------------------------#
      #----- Carbon use efficiency. -------------------------------------------------------#
      last.1yr.gpp           = ifelse( emean$last.1yr.transp[m] * yr.day >= 1.0
                                     , emean$last.1yr.gpp   [m]
                                     , NA
                                     )#end ifelse
      last.2yr.gpp           = ifelse( emean$last.2yr.transp[m] * yr.day >= 1.0
                                     , emean$last.2yr.gpp   [m]
                                     , NA
                                  )#end ifelse
      last.3yr.gpp           = ifelse( emean$last.3yr.transp[m] * yr.day >= 1.0
                                     , emean$last.3yr.gpp   [m]
                                     , NA
                                     )#end ifelse
      emean$last.1yr.cue [m] =         emean$last.1yr.npp[m] / last.1yr.gpp
      emean$last.2yr.cue [m] =         emean$last.2yr.npp[m] / last.2yr.gpp
      emean$last.3yr.cue [m] =         emean$last.3yr.npp[m] / last.3yr.gpp
      #----- Effective CUE. ---------------------------------------------------------------#
      last.1yr.gpp           = ifelse( emean$last.1yr.transp[m] * yr.day >= 1.0
                                     , emean$last.1yr.gpp   [m]
                                     , NA
                                     )#end ifelse
      last.2yr.gpp           = ifelse( emean$last.2yr.transp[m] * yr.day >= 1.0
                                     , emean$last.2yr.gpp   [m]
                                     , NA
                                  )#end ifelse
      last.3yr.gpp           = ifelse( emean$last.3yr.transp[m] * yr.day >= 1.0
                                     , emean$last.3yr.gpp   [m]
                                     , NA
                                     )#end ifelse
      emean$last.1yr.ecue[m] =         emean$last.1yr.dcbadt[m] / last.1yr.gpp
      emean$last.2yr.ecue[m] =         emean$last.2yr.dcbadt[m] / last.2yr.gpp
      emean$last.3yr.ecue[m] =         emean$last.3yr.dcbadt[m] / last.3yr.gpp
      #----- Water use efficiency. --------------------------------------------------------#
      last.1yr.transp        = ifelse( emean$last.1yr.transp[m] * yr.day >= 1.0
                                     , emean$last.1yr.transp[m] * yr.day
                                     , NA
                                     )#end ifelse
      last.2yr.transp        = ifelse( emean$last.2yr.transp[m] * yr.day >= 1.0
                                     , emean$last.2yr.transp[m] * yr.day
                                     , NA
                                  )#end ifelse
      last.3yr.transp        = ifelse( emean$last.3yr.transp[m] * yr.day >= 1.0
                                     , emean$last.3yr.transp[m] * yr.day
                                     , NA
                                     )#end ifelse
      emean$last.1yr.wue [m] = 1000. * emean$last.1yr.npp[m] / last.1yr.transp
      emean$last.2yr.wue [m] = 1000. * emean$last.2yr.npp[m] / last.2yr.transp
      emean$last.3yr.wue [m] = 1000. * emean$last.3yr.npp[m] / last.3yr.transp
      #----- Bulk water use efficiency. ---------------------------------------------------#
      last.1yr.et            = ifelse( emean$last.1yr.transp[m] * yr.day >= 1.0
                                     , emean$last.1yr.et    [m] * yr.day
                                     , NA
                                     )#end ifelse
      last.2yr.et            = ifelse( emean$last.2yr.transp[m] * yr.day >= 1.0
                                     , emean$last.2yr.et    [m] * yr.day
                                     , NA
                                  )#end ifelse
      last.3yr.et            = ifelse( emean$last.3yr.transp[m] * yr.day >= 1.0
                                     , emean$last.3yr.et    [m] * yr.day
                                     , NA
                                     )#end ifelse
      emean$last.1yr.etue[m] = 1000. * emean$last.1yr.npp[m] / last.1yr.et
      emean$last.2yr.etue[m] = 1000. * emean$last.2yr.npp[m] / last.2yr.et
      emean$last.3yr.etue[m] = 1000. * emean$last.3yr.npp[m] / last.3yr.et
      #----- Rain use efficiency. ---------------------------------------------------------#
      last.1yr.rain          = ifelse( emean$last.1yr.rain  [m] >= 1.0
                                     , emean$last.1yr.rain  [m]
                                     , NA
                                     )#end ifelse
      last.2yr.rain          = ifelse( emean$last.2yr.rain  [m] >= 1.0
                                     , emean$last.2yr.rain  [m]
                                     , NA
                                  )#end ifelse
      last.3yr.rain          = ifelse( emean$last.3yr.rain  [m] >= 1.0
                                     , emean$last.3yr.rain  [m]
                                     , NA
                                     )#end ifelse
      emean$last.1yr.rue [m] = 1000. * emean$last.1yr.npp[m] / last.1yr.rain
      emean$last.2yr.rue [m] = 1000. * emean$last.2yr.npp[m] / last.2yr.rain
      emean$last.3yr.rue [m] = 1000. * emean$last.3yr.npp[m] / last.3yr.rain
      #------------------------------------------------------------------------------------#




      #----- Find drought length using rainfall running average. --------------------------#
      if ( m == 1){
         ra.rain              = emean$last.1yr.rain[m] / 12
         emean$nmon.lt.090[m] = as.numeric(emean$last.1yr.rain[m] <  90  )
         emean$nmon.lt.100[m] = as.numeric(emean$last.1yr.rain[m] < 100  )
         emean$nmon.lt.110[m] = as.numeric(emean$last.1yr.rain[m] < 110  )
         emean$nmon.lt.120[m] = as.numeric(emean$last.1yr.rain[m] < 120  )
         emean$nmon.wdef  [m] = as.numeric(emean$water.deficit[m] >  10  )
         emean$nmon.mdef  [m] = as.numeric(emean$malhi.deficit[m] >  10  )
      }else{
         ra.rain              = emean$last.1yr.rain[m] / 12
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
         cohort$biomass      [[clab]] = biomassconow
         cohort$lai          [[clab]] = laiconow
         cohort$wai          [[clab]] = waiconow
         cohort$tai          [[clab]] = taiconow
         cohort$gpp          [[clab]] = gppconow
         cohort$leaf.resp    [[clab]] = leaf.respconow
         cohort$stem.resp    [[clab]] = stem.respconow
         cohort$root.resp    [[clab]] = root.respconow
         cohort$froot.resp   [[clab]] = froot.respconow
         cohort$croot.resp   [[clab]] = croot.respconow
         cohort$plant.resp   [[clab]] = plant.respconow
         cohort$assim.light  [[clab]] = assim.lightconow
         cohort$assim.rubp   [[clab]] = assim.rubpconow
         cohort$assim.co2    [[clab]] = assim.co2conow
         cohort$assim.ratio  [[clab]] = assim.ratioconow
         cohort$npp          [[clab]] = nppconow
         cohort$cba          [[clab]] = cbaconow
         cohort$cbamax       [[clab]] = cbamaxconow
         cohort$cbalight     [[clab]] = cbalightconow
         cohort$cbamoist     [[clab]] = cbamoistconow
         cohort$cbarel       [[clab]] = cbarelconow
         cohort$mcost        [[clab]] = mcostconow
         cohort$ldrop        [[clab]] = ldropconow
         cohort$dcbadt       [[clab]] = dcbadtconow
         cohort$sm.stress    [[clab]] = sm.stressconow
         cohort$light        [[clab]] = lightconow
         cohort$light.beam   [[clab]] = light.beamconow
         cohort$light.diff   [[clab]] = light.diffconow
         cohort$balive       [[clab]] = baliveconow
         cohort$bdead        [[clab]] = bdeadconow
         cohort$bleaf        [[clab]] = bleafconow
         cohort$bstem        [[clab]] = bstemconow
         cohort$broot        [[clab]] = brootconow
         cohort$bfroot       [[clab]] = bfrootconow
         cohort$bcroot       [[clab]] = bcrootconow
         cohort$bsapwood     [[clab]] = bsapwoodconow
         cohort$bstorage     [[clab]] = bstorageconow
         cohort$bseeds       [[clab]] = bseedsconow
         cohort$hflxlc       [[clab]] = hflxlcconow
         cohort$wflxlc       [[clab]] = wflxlcconow
         cohort$transp       [[clab]] = transpconow
         cohort$wue          [[clab]] = wueconow
         cohort$cue          [[clab]] = cueconow
         cohort$ecue         [[clab]] = ecueconow
         cohort$etue         [[clab]] = etueconow
         cohort$demand       [[clab]] = demandconow
         cohort$supply       [[clab]] = supplyconow
         cohort$mort         [[clab]] = 100. * (1.0 - exp(-mortconow      ))
         cohort$ncbmort      [[clab]] = 100. * (1.0 - exp(-ncbmortconow   ))
         cohort$dimort       [[clab]] = 100. * (1.0 - exp(-dimortconow    ))
         cohort$recruit      [[clab]] = recruitconow
         cohort$growth       [[clab]] = 100. * growthconow
         cohort$agb.growth   [[clab]] = 100. * agb.growthconow
         cohort$bsa.growth   [[clab]] = 100. * bsa.growthconow
         cohort$f.gpp        [[clab]] = f.gppconow
         cohort$f.plant.resp [[clab]] = f.plant.respconow
         cohort$f.npp        [[clab]] = f.nppconow
         cohort$f.mco        [[clab]] = f.mcoconow
         cohort$f.cba        [[clab]] = f.cbaconow
         cohort$f.bstorage   [[clab]] = f.bstorageconow
         cohort$f.bleaf      [[clab]] = f.bleafconow
         cohort$f.bstem      [[clab]] = f.bstemconow
         cohort$f.broot      [[clab]] = f.brootconow
         cohort$f.bseeds     [[clab]] = f.bseedsconow
         cohort$f.dcbadt     [[clab]] = f.dcbadtconow
         cohort$leaf.par     [[clab]] = leaf.parconow
         cohort$leaf.par.beam[[clab]] = leaf.par.beamconow
         cohort$leaf.par.diff[[clab]] = leaf.par.diffconow
         cohort$leaf.rshort  [[clab]] = leaf.rshortconow
         cohort$leaf.rlong   [[clab]] = leaf.rlongconow
         cohort$leaf.gpp     [[clab]] = leaf.gppconow
         cohort$rue          [[clab]] = rueconow
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
   datum$qpatch = qpatch
   datum$cohort = cohort
   #---------------------------------------------------------------------------------------#

   return(datum)
   #---------------------------------------------------------------------------------------#
}#end function read.q.files
#==========================================================================================#
#==========================================================================================#
