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
        cat(h5file,"\n")
        mymont    = H5Fopen(h5file)

      }else if(file.exists(h5file.bz2)){
         temp.file = file.path(tempdir(),basename(h5file))
         dummy     = bunzip2(filename=h5file.bz2,destname=temp.file,remove=FALSE)
         mymont    = H5Fopen(temp.file)
         dummy     = file.remove(temp.file)

      }else if(file.exists(h5file.gz)){
         temp.file = file.path(tempdir(),basename(h5file))
         dummy     = gunzip(filename=h5file.gz,destname=temp.file,remove=FALSE)
         mymont    = H5Fopen(temp.file)
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
      mymont$MMEAN_ATM_RSHORT_BEAM_PY = ( mymont$MMEAN_ATM_RSHORT_PY 
                                        - mymont$MMEAN_ATM_RSHORT_DIFF_PY )
      mymont$MMEAN_ATM_PAR_BEAM_PY    = ( mymont$MMEAN_ATM_PAR_PY 
                                        - mymont$MMEAN_ATM_PAR_DIFF_PY    )
      mymont$MMEAN_ATM_RSHORT_BEAM_SI = ( mymont$MMEAN_ATM_RSHORT_SI 
                                        - mymont$MMEAN_ATM_RSHORT_DIFF_SI )
      mymont$MMEAN_ATM_PAR_BEAM_SI    = ( mymont$MMEAN_ATM_PAR_SI 
                                        - mymont$MMEAN_ATM_PAR_DIFF_SI    )
      mymont$QMEAN_ATM_RSHORT_BEAM_PY = ( mymont$QMEAN_ATM_RSHORT_PY 
                                        - mymont$QMEAN_ATM_RSHORT_DIFF_PY )
      mymont$QMEAN_ATM_PAR_BEAM_PY    = ( mymont$QMEAN_ATM_PAR_PY 
                                        - mymont$QMEAN_ATM_PAR_DIFF_PY    )
      mymont$QMEAN_ATM_RSHORT_BEAM_SI = ( mymont$QMEAN_ATM_RSHORT_SI 
                                        - mymont$QMEAN_ATM_RSHORT_DIFF_SI )
      mymont$QMEAN_ATM_PAR_BEAM_SI    = ( mymont$QMEAN_ATM_PAR_SI 
                                        - mymont$QMEAN_ATM_PAR_DIFF_SI    )
      #----- Near infrared. ---------------------------------------------------------------#
      mymont$MMEAN_ATM_NIR_PY         = ( mymont$MMEAN_ATM_RSHORT_PY
                                        - mymont$MMEAN_ATM_PAR_PY          )
      mymont$MMEAN_ATM_NIR_DIFF_PY    = ( mymont$MMEAN_ATM_RSHORT_DIFF_PY
                                        - mymont$MMEAN_ATM_PAR_DIFF_PY     )
      mymont$MMEAN_ATM_NIR_BEAM_PY    = ( mymont$MMEAN_ATM_RSHORT_BEAM_PY
                                        - mymont$MMEAN_ATM_PAR_BEAM_PY     )
      mymont$QMEAN_ATM_NIR_PY         = ( mymont$QMEAN_ATM_RSHORT_PY
                                        - mymont$QMEAN_ATM_PAR_PY          )
      mymont$QMEAN_ATM_NIR_DIFF_PY    = ( mymont$QMEAN_ATM_RSHORT_DIFF_PY
                                        - mymont$QMEAN_ATM_PAR_DIFF_PY     )
      mymont$QMEAN_ATM_NIR_BEAM_PY    = ( mymont$QMEAN_ATM_RSHORT_BEAM_PY
                                        - mymont$QMEAN_ATM_PAR_BEAM_PY     )
      #------------------------------------------------------------------------------------#

      #----- Make a growth respiration variable. ------------------------------------------#
      if (! "MMEAN_GROWTH_RESP_PY" %in% names(mymont)){
         mymont$MMEAN_GROWTH_RESP_PY  = ( mymont$MMEAN_LEAF_GROWTH_RESP_PY
                                        + mymont$MMEAN_ROOT_GROWTH_RESP_PY
                                        + mymont$MMEAN_SAPA_GROWTH_RESP_PY
                                        + mymont$MMEAN_SAPB_GROWTH_RESP_PY )
         mymont$QMEAN_GROWTH_RESP_PY  = ( mymont$QMEAN_LEAF_GROWTH_RESP_PY
                                        + mymont$QMEAN_ROOT_GROWTH_RESP_PY
                                        + mymont$QMEAN_SAPA_GROWTH_RESP_PY
                                        + mymont$QMEAN_SAPB_GROWTH_RESP_PY )
      }#end if (! "MMEAN_GROWTH_RESP_PY" %in% names(mymont))
      if (! "MMEAN_VLEAF_RESP_PY" %in% names(mymont)){
         mymont$MMEAN_VLEAF_RESP_PY = 0. * mymont$MMEAN_LEAF_RESP_PY
         mymont$MMEAN_VLEAF_RESP_CO = 0. * mymont$MMEAN_LEAF_RESP_CO
         mymont$QMEAN_VLEAF_RESP_PY = 0. * mymont$QMEAN_LEAF_RESP_PY
         mymont$QMEAN_VLEAF_RESP_CO = 0. * mymont$QMEAN_LEAF_RESP_CO
      }#end if
      #------------------------------------------------------------------------------------#

      #fabio Add a variable for storage respiration
      mymont$MMEAN_STORAGE_RESP_PY = 0

      #------------------------------------------------------------------------------------#
      #      Fix units for some variables.                                                 #
      #------------------------------------------------------------------------------------#
      if (corr.growth.storage != 1.){
         cat("     - Correcting growth and storage respiration...","\n")
         #----- Correct the rates. --------------------------------------------------------#
         multfac = corr.growth.storage
         mymont$MMEAN_GROWTH_RESP_PY  =   mymont$MMEAN_GROWTH_RESP_PY  * multfac
         #mymont$MMEAN_STORAGE_RESP_PY =   mymont$MMEAN_STORAGE_RESP_PY * multfac
         mymont$MMEAN_VLEAF_RESP_PY   =   mymont$MMEAN_VLEAF_RESP_PY   * multfac
         mymont$QMEAN_GROWTH_RESP_PY  =   mymont$QMEAN_GROWTH_RESP_PY  * multfac
         mymont$QMEAN_STORAGE_RESP_PY = 0
         #mymont$QMEAN_STORAGE_RESP_PY =   mymont$QMEAN_STORAGE_RESP_PY * multfac
         mymont$QMEAN_VLEAF_RESP_PY   =   mymont$QMEAN_VLEAF_RESP_PY   * multfac
         if (mymont$NCOHORTS_GLOBAL > 0){
            mymont$MMEAN_GROWTH_RESP_CO  = mymont$MMEAN_GROWTH_RESP_CO  * multfac
            mymont$MMEAN_STORAGE_RESP_CO = 0
            mymont$MMEAN_STORAGE_RESP_CO = mymont$MMEAN_STORAGE_RESP_CO * multfac
            mymont$MMEAN_VLEAF_RESP_CO   = mymont$MMEAN_VLEAF_RESP_CO   * multfac
            mymont$QMEAN_GROWTH_RESP_CO  = mymont$QMEAN_GROWTH_RESP_CO  * multfac
            mymont$QMEAN_STORAGE_RESP_CO = mymont$QMEAN_STORAGE_RESP_CO * multfac
            mymont$QMEAN_VLEAF_RESP_CO   = mymont$QMEAN_VLEAF_RESP_CO   * multfac
         }#end if
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #      Correct Plant respiration_                                                 #
         #---------------------------------------------------------------------------------#
         mymont$MMEAN_PLRESP_PY       = ( mymont$MMEAN_LEAF_RESP_PY
                                        + mymont$MMEAN_ROOT_RESP_PY
                                        + mymont$MMEAN_GROWTH_RESP_PY
                                        + mymont$MMEAN_STORAGE_RESP_PY
                                        + mymont$MMEAN_VLEAF_RESP_PY
                                        )#end 
         mymont$QMEAN_PLRESP_PY       = ( mymont$QMEAN_LEAF_RESP_PY
                                        + mymont$QMEAN_ROOT_RESP_PY
                                        + mymont$QMEAN_GROWTH_RESP_PY
                                        + mymont$QMEAN_STORAGE_RESP_PY
                                        + mymont$QMEAN_VLEAF_RESP_PY
                                        )#end 
         mymont$MMSQU_PLRESP_PY       = mymont$MMSQU_PLRESP_PY + NA
         mymont$QMSQU_PLRESP_PY       = mymont$QMSQU_PLRESP_PY + NA
         if (mymont$NCOHORTS_GLOBAL > 0){
            mymont$MMEAN_PLRESP_CO       = ( mymont$MMEAN_LEAF_RESP_CO
                                           + mymont$MMEAN_ROOT_RESP_CO
                                           + mymont$MMEAN_GROWTH_RESP_CO
                                           + mymont$MMEAN_STORAGE_RESP_CO
                                           + mymont$MMEAN_VLEAF_RESP_CO
                                           )#end 
            mymont$QMEAN_PLRESP_CO       = ( mymont$QMEAN_LEAF_RESP_CO
                                           + mymont$QMEAN_ROOT_RESP_CO
                                           + mymont$QMEAN_GROWTH_RESP_CO
                                           + mymont$QMEAN_STORAGE_RESP_CO
                                           + mymont$QMEAN_VLEAF_RESP_CO
                                           )#end 
            mymont$MMSQU_PLRESP_CO       = mymont$MMSQU_PLRESP_CO + NA
            mymont$QMSQU_PLRESP_CO       = mymont$QMSQU_PLRESP_CO + NA
         }#end if
         #---------------------------------------------------------------------------------#




         #---------------------------------------------------------------------------------#
         #      Correct NPP_                                                               #
         #---------------------------------------------------------------------------------#
         mymont$MMEAN_NPP_PY    = mymont$MMEAN_GPP_PY - mymont$MMEAN_PLRESP_PY
         mymont$QMEAN_NPP_PY    = mymont$QMEAN_GPP_PY - mymont$QMEAN_PLRESP_PY
         mymont$MMSQU_NPP_PY    = mymont$MMSQU_NPP_PY + NA
         mymont$QMSQU_NPP_PY    = mymont$QMSQU_NPP_PY + NA
         if (mymont$NCOHORTS_GLOBAL > 0){
            mymont$MMEAN_NPP_CO = mymont$MMEAN_GPP_CO - mymont$MMEAN_PLRESP_CO
            mymont$QMEAN_NPP_CO = mymont$QMEAN_GPP_CO - mymont$QMEAN_PLRESP_CO
            mymont$MMSQU_NPP_CO = mymont$MMSQU_NPP_CO + NA
            mymont$QMSQU_NPP_CO = mymont$QMSQU_NPP_CO + NA
         }#end if
         #---------------------------------------------------------------------------------#




         #---------------------------------------------------------------------------------#
         #      Correct NEP_                                                               #
         #---------------------------------------------------------------------------------#
         mymont$MMEAN_NEP_PY    = mymont$MMEAN_NPP_PY - mymont$MMEAN_RH_PY
         mymont$QMEAN_NEP_PY    = mymont$QMEAN_NPP_PY - mymont$QMEAN_RH_PY
         mymont$MMSQU_NEP_PY    = mymont$MMSQU_NEP_PY + NA
         mymont$QMSQU_NEP_PY    = mymont$QMSQU_NEP_PY + NA

         mymont$MMEAN_NEP_PA    = -mymont$MMEAN_RH_PA
         mymont$QMEAN_NEP_PA    = -mymont$QMEAN_RH_PA
         mymont$MMSQU_NEP_PA    = mymont$MMSQU_NEP_PA + NA
         mymont$QMSQU_NEP_PA    = mymont$QMSQU_NEP_PA + NA
         if (mymont$NCOHORTS_GLOBAL > 0){
            ipaconow    = rep(sequence(mymont$NPATCHES_GLOBAL),times=mymont$PACO_N)
            idx         = match(unique(ipaconow),sequence(mymont$NPATCHES_GLOBAL))
            mymont$MMEAN_NEP_PA[idx]  = ( tapply( X     = mymont$MMEAN_NPP_CO
                                                        * mymont$NPLANT
                                                , INDEX = ipaconow
                                                , FUN   = sum
                                                )#end tapply
                                        - mymont$MMEAN_RH_PA[idx]
                                        )#end
            mymont$QMEAN_NEP_PA[idx,] = ( qapply( X     = mymont$QMEAN_NPP_CO
                                                        * mymont$NPLANT
                                                , INDEX = ipaconow
                                                , DIM   = 1
                                                , FUN   = sum
                                                )#end tapply
                                        - mymont$QMEAN_RH_PA[idx,]
                                        )#end
         }#end if
         #---------------------------------------------------------------------------------#
      }#end if
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Set daytime flag.                                                              #
      #------------------------------------------------------------------------------------#
      phap        = mymont$QMEAN_ATM_RSHORT_PY > phap.min | iint.photo == 0
      polar.night = ! any(phap,na.rm=TRUE)
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Find the mean latent heat of vaporisation.  Because we assume it to be a       #
      # linear function of temperature, the mean can be found a posteriori.  The mean      #
      # fluxes won't be exact though, because the covariance part is missing.              #
      #------------------------------------------------------------------------------------#
      mmean.can.alvli.py  = alvli(mymont$MMEAN_CAN_TEMP_PY)
      qmean.can.alvli.py  = alvli(mymont$QMEAN_CAN_TEMP_PY)
      mmean.can.alvli.pa  = alvli(mymont$MMEAN_CAN_TEMP_PA)
      qmean.can.alvli.pa  = alvli(mymont$QMEAN_CAN_TEMP_PA)
      mmean.can.alvli.py2 = mmean.can.alvli.py * mmean.can.alvli.py
      qmean.can.alvli.py2 = qmean.can.alvli.py * qmean.can.alvli.py
      mmean.can.alvli.pa2 = mmean.can.alvli.pa * mmean.can.alvli.pa
      qmean.can.alvli.pa2 = qmean.can.alvli.pa * qmean.can.alvli.pa
      #------------------------------------------------------------------------------------#



      #----- Load the total number of patches and cohorts. --------------------------------#
      emean$npat.global[m] = mymont$NPATCHES_GLOBAL
      emean$ncoh.global[m] = mymont$NCOHORTS_GLOBAL
      #------------------------------------------------------------------------------------#


      #----- Load the simple variables_ ---------------------------------------------------#
      emean$fast.soil.c     [m] =   mymont$MMEAN_FAST_SOIL_C_PY
      emean$slow.soil.c     [m] =   mymont$MMEAN_SLOW_SOIL_C_PY
      emean$struct.soil.c   [m] =   mymont$MMEAN_STRUCT_SOIL_C_PY
      emean$het.resp        [m] =   mymont$MMEAN_RH_PY
      emean$cwd.resp        [m] =   mymont$MMEAN_CWD_RH_PY
      emean$gpp             [m] =   mymont$MMEAN_GPP_PY
      emean$npp             [m] =   mymont$MMEAN_NPP_PY
      emean$nep             [m] =   mymont$MMEAN_NEP_PY
      emean$nee             [m] =   mymont$MMEAN_CARBON_ST_PY - mymont$MMEAN_CARBON_AC_PY
      emean$plant.resp      [m] =   mymont$MMEAN_PLRESP_PY
      emean$growth.resp     [m] =   mymont$MMEAN_GROWTH_RESP_PY
      emean$storage.resp    [m] =   mymont$MMEAN_STORAGE_RESP_PY
      emean$assim.light     [m] =   mymont$MMEAN_A_LIGHT_PY
      emean$assim.rubp      [m] =   mymont$MMEAN_A_RUBP_PY
      emean$assim.co2       [m] =   mymont$MMEAN_A_CO2_PY
      emean$assim.ratio     [m] = ( mymont$MMEAN_A_LIGHT_PY
                                  / max( 1e-6,min( mymont$MMEAN_A_RUBP_PY
                                                 , mymont$MMEAN_A_CO2_PY  )))
      #----- (Leaf, root, stem, and soil respiration are corrected for growth+storage)_ ---#
      emean$reco            [m] =   mymont$MMEAN_PLRESP_PY    + mymont$MMEAN_RH_PY
      emean$ustar           [m] =   mymont$MMEAN_USTAR_PY
      emean$cflxca          [m] = - mymont$MMEAN_CARBON_AC_PY
      emean$cflxst          [m] =   mymont$MMEAN_CARBON_ST_PY
      emean$ustar           [m] =   mymont$MMEAN_USTAR_PY
      emean$atm.vels        [m] =   mymont$MMEAN_ATM_VELS_PY
      emean$atm.prss        [m] =   mymont$MMEAN_ATM_PRSS_PY   * 0.01
      emean$atm.temp        [m] =   mymont$MMEAN_ATM_TEMP_PY   - t00
      emean$atm.shv         [m] =   mymont$MMEAN_ATM_SHV_PY    * kg2g
      emean$atm.co2         [m] =   mymont$MMEAN_ATM_CO2_PY
      emean$atm.vpd         [m] =   mymont$MMEAN_ATM_VPDEF_PY  * 0.01
      emean$can.prss        [m] =   mymont$MMEAN_CAN_PRSS_PY   * 0.01
      emean$can.temp        [m] =   mymont$MMEAN_CAN_TEMP_PY   - t00
      emean$can.shv         [m] =   mymont$MMEAN_CAN_SHV_PY    * kg2g
      emean$can.co2         [m] =   mymont$MMEAN_CAN_CO2_PY
      emean$can.vpd         [m] =   mymont$MMEAN_CAN_VPDEF_PY  * 0.01
      emean$gnd.temp        [m] =   mymont$MMEAN_GND_TEMP_PY   - t00
      emean$gnd.shv         [m] =   mymont$MMEAN_GND_SHV_PY    * kg2g
      emean$leaf.temp       [m] =   mymont$MMEAN_LEAF_TEMP_PY  - t00
      emean$leaf.water      [m] =   mymont$MMEAN_LEAF_WATER_PY
      emean$leaf.vpd        [m] =   mymont$MMEAN_LEAF_VPDEF_PY * 0.01
      emean$wood.temp       [m] =   mymont$MMEAN_WOOD_TEMP_PY  - t00
      emean$hflxca          [m] = - mymont$MMEAN_SENSIBLE_AC_PY
      emean$qwflxca         [m] = - mymont$MMEAN_VAPOR_AC_PY   * mmean.can.alvli.py
      emean$hflxgc          [m] =   mymont$MMEAN_SENSIBLE_GC_PY
      emean$hflxlc          [m] =   mymont$MMEAN_SENSIBLE_LC_PY
      emean$hflxwc          [m] =   mymont$MMEAN_SENSIBLE_WC_PY
      emean$wflxca          [m] = - mymont$MMEAN_VAPOR_AC_PY   * day.sec
      emean$wflxgc          [m] =   mymont$MMEAN_VAPOR_GC_PY   * day.sec
      emean$wflxlc          [m] =   mymont$MMEAN_VAPOR_LC_PY   * day.sec
      emean$wflxwc          [m] =   mymont$MMEAN_VAPOR_WC_PY   * day.sec
      emean$runoff          [m] = ( mymont$MMEAN_RUNOFF_PY
                                  + mymont$MMEAN_DRAINAGE_PY       ) * mondays * day.sec
      emean$intercepted     [m] = ( mymont$MMEAN_INTERCEPTED_AL_PY
                                  + mymont$MMEAN_INTERCEPTED_AW_PY ) * mondays * day.sec
      emean$wshed           [m] = ( mymont$MMEAN_WSHED_LG_PY
                                  + mymont$MMEAN_WSHED_WG_PY       ) * mondays * day.sec
      emean$evap            [m] = ( mymont$MMEAN_VAPOR_GC_PY
                                  + mymont$MMEAN_VAPOR_LC_PY
                                  + mymont$MMEAN_VAPOR_WC_PY ) * day.sec
      emean$transp          [m] =   mymont$MMEAN_TRANSP_PY     * day.sec
      emean$et              [m] = emean$evap[m] + emean$transp[m]
      emean$rain            [m] = mymont$MMEAN_PCPG_PY * mondays * day.sec

      emean$sm.stress       [m] =   1. - mymont$MMEAN_FS_OPEN_PY
      emean$rshort          [m] =   mymont$MMEAN_ATM_RSHORT_PY
      emean$rshort.beam     [m] = ( mymont$MMEAN_ATM_RSHORT_PY
                                  - mymont$MMEAN_ATM_RSHORT_DIFF_PY )
      emean$rshort.diff     [m] =   mymont$MMEAN_ATM_RSHORT_DIFF_PY
      emean$rshortup        [m] =   mymont$MMEAN_RSHORTUP_PY
      emean$rlong           [m] =   mymont$MMEAN_ATM_RLONG_PY
      emean$rshort.gnd      [m] =   mymont$MMEAN_RSHORT_GND_PY
      emean$rlong.gnd       [m] =   mymont$MMEAN_RLONG_GND_PY
      emean$rlongup         [m] =   mymont$MMEAN_RLONGUP_PY
      emean$par.tot         [m] =   mymont$MMEAN_ATM_PAR_PY        * Watts.2.Ein * 1.e6
      emean$par.beam        [m] = ( mymont$MMEAN_ATM_PAR_PY
                                  - mymont$MMEAN_ATM_PAR_DIFF_PY ) * Watts.2.Ein * 1.e6
      emean$par.diff        [m] =   mymont$MMEAN_ATM_PAR_DIFF_PY   * Watts.2.Ein * 1.e6
      emean$par.gnd         [m] =   mymont$MMEAN_PAR_GND_PY        * Watts.2.Ein * 1.e6
      emean$parup           [m] =   mymont$MMEAN_PARUP_PY          * Watts.2.Ein * 1.e6
      emean$rnet            [m] =   mymont$MMEAN_RNET_PY
      emean$albedo          [m] =   mymont$MMEAN_ALBEDO_PY
      if (all(c("MMEAN_ALBEDO_PAR_PY","MMEAN_ALBEDO_NIR_PY") %in% names(mymont))){
         emean$albedo.par   [m] =   mymont$MMEAN_ALBEDO_PAR_PY
         emean$albedo.nir   [m] =   mymont$MMEAN_ALBEDO_NIR_PY
      }else{
         emean$albedo.par   [m] = ifelse( mymont$MMEAN_ATM_PAR_PY > 0.5
                                        , mymont$MMEAN_PARUP_PY / mymont$MMEAN_ATM_PAR_PY
                                        , mymont$MMEAN_ALBEDO_PY
                                        )#end ifelse
         emean$albedo.nir   [m] = ifelse( mymont$MMEAN_ATM_NIR_PY > 0.5
                                        , mymont$MMEAN_NIRUP_PY / mymont$MMEAN_ATM_NIR_PY
                                        , mymont$MMEAN_ALBEDO_PY
                                        )#end ifelse
      }#end if
      emean$rlong.albedo    [m] =   mymont$MMEAN_RLONG_ALBEDO_PY
      emean$leaf.gbw        [m] =   mymont$MMEAN_LEAF_GBW_PY * day.sec
      emean$leaf.gsw        [m] =   mymont$MMEAN_LEAF_GSW_PY * day.sec
      emean$wood.gbw        [m] =   mymont$MMEAN_WOOD_GBW_PY * day.sec
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     The following variables must be aggregated because the polygon-level is split  #
      # by PFT and DBH class.                                                              #
      #------------------------------------------------------------------------------------#
      ## fabio warnings number of items to replace is not a multiple of replacement
      emean$mco             [m] = apply( X      = ( mymont$MMEAN_LEAF_MAINTENANCE_PY
                                                  + mymont$MMEAN_ROOT_MAINTENANCE_PY
                                                  ) * yr.day
                                       , MARGIN = 1
                                       , FUN    = sum
                                       )#end if
      emean$ldrop           [m] = apply( X      = mymont$MMEAN_LEAF_DROP_PY * yr.day
                                       , MARGIN = 1
                                       , FUN    = sum
                                       )#end if
      #------------------------------------------------------------------------------------#




      #------ Read in soil properties. ----------------------------------------------------#
      emean$soil.temp     [m,] =   mymont$MMEAN_SOIL_TEMP_PY - t00
      emean$soil.water    [m,] =   mymont$MMEAN_SOIL_WATER_PY
      emean$soil.mstpot   [m,] = - mymont$MMEAN_SOIL_MSTPOT_PY * grav * wdnsi
      emean$soil.extracted[m,] = - mymont$MMEAN_TRANSLOSS_PY * day.sec / dslz
      #------------------------------------------------------------------------------------#



      #----- Find averaged soil properties. -----------------------------------------------#
      swater.now     = rev(cumsum(rev(mymont$MMEAN_SOIL_WATER_PY * wdns * dslz)))
      smoist.avg     = swater.now / (wdns * soil.depth)
      emean$paw  [m] = 100. * ( ( swater.now[ka] - soil.dry [ka] )
                              / ( soil.poro [ka] - soil.dry [ka] ) )
      emean$smpot[m] = ( - smoist2mpot(smoist=smoist.avg[ka],mysoil=soil.prop)
                       * 0.001 * grav )
      #------------------------------------------------------------------------------------#



      #----- Read workload, and retrieve only the current month. --------------------------#
      emean$workload  [m] = mymont$WORKLOAD[thismonth]
      emean$specwork  [m] = mymont$WORKLOAD[thismonth] / sum(mymont$SIPA_N,na.rm=TRUE)
      #------------------------------------------------------------------------------------#





      #------------------------------------------------------------------------------------#
      #     Retrieve the sum of squares (that will be used to find standard deviation.     #
      #------------------------------------------------------------------------------------#
      emsqu$gpp       [m] =   mymont$MMSQU_GPP_PY
      emsqu$plant.resp[m] =   mymont$MMSQU_PLRESP_PY
      emsqu$het.resp  [m] =   mymont$MMSQU_RH_PY
      ##fabio MMSQU WD doesn't exist in h5 file 
      #emsqu$cwd.resp  [m] =   mymont$MMSQU_WD_RH_PY
      emsqu$cflxca    [m] =   mymont$MMSQU_CARBON_AC_PY
      emsqu$cflxst    [m] =   mymont$MMSQU_CARBON_ST_PY
      emsqu$hflxca    [m] =   mymont$MMSQU_SENSIBLE_AC_PY
      emsqu$hflxlc    [m] =   mymont$MMSQU_SENSIBLE_LC_PY
      emsqu$hflxwc    [m] =   mymont$MMSQU_SENSIBLE_WC_PY
      emsqu$hflxgc    [m] =   mymont$MMSQU_SENSIBLE_GC_PY
      emsqu$wflxca    [m] =   mymont$MMSQU_VAPOR_AC_PY  * day.sec2
      emsqu$qwflxca   [m] =   mymont$MMSQU_VAPOR_AC_PY  * mmean.can.alvli.py2
      emsqu$wflxlc    [m] =   mymont$MMSQU_VAPOR_LC_PY  * day.sec2
      emsqu$wflxwc    [m] =   mymont$MMSQU_VAPOR_WC_PY  * day.sec2
      emsqu$wflxgc    [m] =   mymont$MMSQU_VAPOR_GC_PY  * day.sec2
      emsqu$transp    [m] =   mymont$MMSQU_TRANSP_PY    * day.sec2
      emsqu$ustar     [m] =   mymont$MMSQU_USTAR_PY
      emsqu$albedo    [m] =   mymont$MMSQU_ALBEDO_PY
      emsqu$rshortup  [m] =   mymont$MMSQU_RSHORTUP_PY
      emsqu$rlongup   [m] =   mymont$MMSQU_RLONGUP_PY
      emsqu$parup     [m] =   mymont$MMSQU_PARUP_PY * Watts.2.Ein^2 * 1.e12
      emsqu$rnet      [m] =   mymont$MMSQU_RNET_PY
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #       Read the mean diurnal cycle and the mean sum of the squares.                 #
      #------------------------------------------------------------------------------------#
      qmean$gpp         [m,] =   mymont$QMEAN_GPP_PY
      qmean$npp         [m,] =   mymont$QMEAN_NPP_PY
      qmean$het.resp    [m,] =   mymont$QMEAN_RH_PY
      qmean$cwd.resp    [m,] =   mymont$QMEAN_CWD_RH_PY
      qmean$assim.light [m,] =   mymont$QMEAN_A_LIGHT_PY
      qmean$assim.rubp  [m,] =   mymont$QMEAN_A_RUBP_PY
      qmean$assim.co2   [m,] =   mymont$QMEAN_A_CO2_PY
      qmean$assim.ratio [m,] = ( mymont$QMEAN_A_LIGHT_PY
                               / pmax(1e-6, pmin( mymont$QMEAN_A_RUBP_PY
                                                , mymont$QMEAN_A_CO2_PY  )))
      qmean$nee         [m,] = ( mymont$QMEAN_CARBON_ST_PY
                               - mymont$QMEAN_CARBON_AC_PY )
      qmean$reco        [m,] =   mymont$QMEAN_PLRESP_PY + mymont$QMEAN_RH_PY
      qmean$cflxca      [m,] = - mymont$QMEAN_CARBON_AC_PY
      qmean$cflxst      [m,] = - mymont$QMEAN_CARBON_ST_PY
      qmean$hflxca      [m,] = - mymont$QMEAN_SENSIBLE_AC_PY
      qmean$hflxlc      [m,] =   mymont$QMEAN_SENSIBLE_LC_PY
      qmean$hflxwc      [m,] =   mymont$QMEAN_SENSIBLE_WC_PY
      qmean$hflxgc      [m,] =   mymont$QMEAN_SENSIBLE_GC_PY
      qmean$wflxca      [m,] = - mymont$QMEAN_VAPOR_AC_PY     * day.sec
      qmean$qwflxca     [m,] = - mymont$QMEAN_VAPOR_AC_PY     * qmean.can.alvli.py
      qmean$wflxlc      [m,] =   mymont$QMEAN_VAPOR_LC_PY     * day.sec
      qmean$wflxwc      [m,] =   mymont$QMEAN_VAPOR_WC_PY     * day.sec
      qmean$wflxgc      [m,] =   mymont$QMEAN_VAPOR_GC_PY     * day.sec
      qmean$runoff      [m,] = ( mymont$QMEAN_RUNOFF_PY
                               + mymont$QMEAN_DRAINAGE_PY )   * day.sec
      qmean$intercepted [m,] = ( mymont$QMEAN_INTERCEPTED_AL_PY
                               + mymont$QMEAN_INTERCEPTED_AW_PY ) * day.sec
      qmean$wshed       [m,] = ( mymont$QMEAN_WSHED_LG_PY
                               + mymont$QMEAN_WSHED_WG_PY       ) * day.sec
      qmean$evap        [m,] = ( mymont$QMEAN_VAPOR_GC_PY
                               + mymont$QMEAN_VAPOR_WC_PY
                               + mymont$QMEAN_VAPOR_LC_PY )   * day.sec
      qmean$transp      [m,] =   mymont$QMEAN_TRANSP_PY       * day.sec
      qmean$atm.temp    [m,] =   mymont$QMEAN_ATM_TEMP_PY     - t00
      qmean$can.temp    [m,] =   mymont$QMEAN_CAN_TEMP_PY     - t00
      qmean$leaf.temp   [m,] =   mymont$QMEAN_LEAF_TEMP_PY    - t00
      qmean$leaf.water  [m,] =   mymont$QMEAN_LEAF_WATER_PY
      qmean$wood.temp   [m,] =   mymont$QMEAN_WOOD_TEMP_PY    - t00
      qmean$gnd.temp    [m,] =   mymont$QMEAN_GND_TEMP_PY     - t00
      qmean$atm.shv     [m,] =   mymont$QMEAN_ATM_SHV_PY      * kg2g
      qmean$can.shv     [m,] =   mymont$QMEAN_CAN_SHV_PY      * kg2g
      qmean$gnd.shv     [m,] =   mymont$QMEAN_GND_SHV_PY      * kg2g
      qmean$atm.vpd     [m,] =   mymont$QMEAN_ATM_VPDEF_PY    * 0.01
      qmean$can.vpd     [m,] =   mymont$QMEAN_CAN_VPDEF_PY    * 0.01
      qmean$leaf.vpd    [m,] =   mymont$QMEAN_LEAF_VPDEF_PY   * 0.01
      qmean$atm.co2     [m,] =   mymont$QMEAN_ATM_CO2_PY
      qmean$can.co2     [m,] =   mymont$QMEAN_CAN_CO2_PY
      qmean$atm.vels    [m,] =   mymont$QMEAN_ATM_VELS_PY
      qmean$ustar       [m,] =   mymont$QMEAN_USTAR_PY
      qmean$atm.prss    [m,] =   mymont$QMEAN_ATM_PRSS_PY     * 0.01
      qmean$can.prss    [m,] =   mymont$QMEAN_CAN_PRSS_PY     * 0.01
      qmean$sm.stress   [m,] =   1. - mymont$QMEAN_FS_OPEN_PY
      qmean$rain        [m,] =   mymont$QMEAN_PCPG_PY         * day.sec
      qmean$rshort      [m,] =   mymont$QMEAN_ATM_RSHORT_PY
      qmean$rshort.beam [m,] = ( mymont$QMEAN_ATM_RSHORT_PY 
                               - mymont$QMEAN_ATM_RSHORT_DIFF_PY )
      qmean$rshort.diff [m,] =   mymont$QMEAN_ATM_RSHORT_DIFF_PY
      qmean$rshort.gnd  [m,] =   mymont$QMEAN_RSHORT_GND_PY
      qmean$rshortup    [m,] =   mymont$QMEAN_RSHORTUP_PY
      qmean$rlong       [m,] =   mymont$QMEAN_ATM_RLONG_PY
      qmean$rlong.gnd   [m,] =   mymont$QMEAN_RLONG_GND_PY
      qmean$rlongup     [m,] =   mymont$QMEAN_RLONGUP_PY
      qmean$par.tot     [m,] =   mymont$QMEAN_ATM_PAR_PY        * Watts.2.Ein * 1.e6
      qmean$par.beam    [m,] = ( mymont$QMEAN_ATM_PAR_PY
                               - mymont$QMEAN_ATM_PAR_DIFF_PY ) * Watts.2.Ein * 1.e6
      qmean$par.diff    [m,] =   mymont$QMEAN_ATM_PAR_DIFF_PY   * Watts.2.Ein * 1.e6
      qmean$par.gnd     [m,] =   mymont$QMEAN_PAR_GND_PY        * Watts.2.Ein * 1.e6
      qmean$parup       [m,] =   mymont$QMEAN_PARUP_PY          * Watts.2.Ein * 1.e6
      qmean$rnet        [m,] =   mymont$QMEAN_RNET_PY
      qmean$albedo      [m,] =   mymont$QMEAN_ALBEDO_PY
      if (all(c("QMEAN_ALBEDO_PAR_PY","QMEAN_ALBEDO_NIR_PY") %in% names(mymont))){
         qmean$albedo.par   [m,] =   mymont$QMEAN_ALBEDO_PAR_PY
         qmean$albedo.nir   [m,] =   mymont$QMEAN_ALBEDO_NIR_PY
      }else{
         qmean$albedo.par   [m,] = ifelse( mymont$QMEAN_ATM_PAR_PY > 0.5
                                         , mymont$QMEAN_PARUP_PY / mymont$QMEAN_ATM_PAR_PY
                                         , mymont$QMEAN_ALBEDO_PY
                                         )#end ifelse
         qmean$albedo.nir   [m,] = ifelse( mymont$QMEAN_ATM_NIR_PY > 0.5
                                         , mymont$QMEAN_NIRUP_PY / mymont$QMEAN_ATM_NIR_PY
                                         , mymont$QMEAN_ALBEDO_PY
                                         )#end ifelse
      }#end if
      qmean$rlong.albedo[m,] =   mymont$QMEAN_RLONG_ALBEDO_PY
      qmean$leaf.gbw    [m,] =   mymont$QMEAN_LEAF_GBW_PY       * day.sec
      qmean$leaf.gsw    [m,] =   mymont$QMEAN_LEAF_GSW_PY       * day.sec
      qmean$wood.gbw    [m,] =   mymont$QMEAN_WOOD_GBW_PY       * day.sec
      #------------------------------------------------------------------------------------#



      #------ Read the mean sum of squares for diel. --------------------------------------#
      qmsqu$gpp         [m,] =   mymont$QMSQU_GPP_PY
      qmsqu$npp         [m,] =   mymont$QMSQU_NPP_PY
      qmsqu$plant.resp  [m,] =   mymont$QMSQU_PLRESP_PY
      qmsqu$het.resp    [m,] =   mymont$QMSQU_RH_PY
      qmsqu$cwd.resp    [m,] =   mymont$QMSQU_CWD_RH_PY
      qmsqu$nep         [m,] =   mymont$QMSQU_NEP_PY
      qmsqu$cflxca      [m,] =   mymont$QMSQU_CARBON_AC_PY
      qmsqu$cflxst      [m,] =   mymont$QMSQU_CARBON_ST_PY
      qmsqu$hflxca      [m,] =   mymont$QMSQU_SENSIBLE_AC_PY
      qmsqu$hflxlc      [m,] =   mymont$QMSQU_SENSIBLE_LC_PY
      qmsqu$hflxwc      [m,] =   mymont$QMSQU_SENSIBLE_WC_PY
      qmsqu$hflxgc      [m,] =   mymont$QMSQU_SENSIBLE_GC_PY
      qmsqu$wflxca      [m,] =   mymont$QMSQU_VAPOR_AC_PY  * day.sec2
      qmsqu$qwflxca     [m,] =   mymont$QMSQU_VAPOR_AC_PY  * qmean.can.alvli.py2
      qmsqu$wflxlc      [m,] =   mymont$QMSQU_VAPOR_WC_PY  * day.sec2
      qmsqu$wflxwc      [m,] =   mymont$QMSQU_VAPOR_LC_PY  * day.sec2
      qmsqu$wflxgc      [m,] =   mymont$QMSQU_VAPOR_GC_PY  * day.sec2
      qmsqu$transp      [m,] =   mymont$QMSQU_TRANSP_PY    * day.sec2
      qmsqu$ustar       [m,] =   mymont$QMSQU_USTAR_PY
      qmsqu$albedo      [m,] =   mymont$QMSQU_ALBEDO_PY
      qmsqu$rshortup    [m,] =   mymont$QMSQU_RSHORTUP_PY
      qmsqu$rlongup     [m,] =   mymont$QMSQU_RLONGUP_PY
      qmsqu$parup       [m,] =   mymont$QMSQU_PARUP_PY     * Watts.2.Ein^2 * 1e12
      #------------------------------------------------------------------------------------#


      #---- Read in the site-level area. --------------------------------------------------#
      areasi     = mymont$AREA_SI
      npatches   = mymont$SIPA_N
      #------------------------------------------------------------------------------------#


      #----- Read a few patch-level variables. --------------------------------------------#
      areapa      = mymont$AREA * rep(areasi,times=npatches)
      areapa      = areapa / sum(areapa)
      ipa         = sequence(mymont$NPATCHES_GLOBAL)
      lupa        = mymont$DIST_TYPE
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
      ncohorts    = mymont$PACO_N
      ipaconow    = rep(sequence(mymont$NPATCHES_GLOBAL),times=mymont$PACO_N)
      icoconow    = unlist(sapply(X = mymont$PACO_N, FUN = sequence))
      idx         = match(unique(ipaconow),sequence(mymont$NPATCHES_GLOBAL))
      #------------------------------------------------------------------------------------#



      #----- Disturbance rates. -----------------------------------------------------------#
      ## fabio, problem with different length X and w
      lu$dist  [m,,] = 0.0125 #apply ( X      = mymont$DISTURBANCE_RATES
                       #      , MARGIN = c(2,3)
                        #     , FUN    = weighted.mean
                         #    , w      = areasi
                          #   )#end apply
      #         browser()
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
      patch$nep          [[plab]] =   mymont$MMEAN_NEP_PA
      patch$het.resp     [[plab]] =   mymont$MMEAN_RH_PA
      patch$can.temp     [[plab]] =   mymont$MMEAN_CAN_TEMP_PA    - t00
      patch$gnd.temp     [[plab]] =   mymont$MMEAN_GND_TEMP_PA    - t00
      patch$can.shv      [[plab]] =   mymont$MMEAN_CAN_SHV_PA     * 1000.
      patch$gnd.shv      [[plab]] =   mymont$MMEAN_GND_SHV_PA     * 1000.
      patch$can.vpd      [[plab]] =   mymont$MMEAN_CAN_VPDEF_PA   * 0.01
      patch$can.co2      [[plab]] =   mymont$MMEAN_CAN_CO2_PA
      patch$can.prss     [[plab]] =   mymont$MMEAN_CAN_PRSS_PA    * 0.01
      patch$cflxca       [[plab]] = - mymont$MMEAN_CARBON_AC_PA
      patch$cflxst       [[plab]] =   mymont$MMEAN_CARBON_ST_PA
      patch$nee          [[plab]] = ( mymont$MMEAN_CARBON_ST_PA
                                    - mymont$MMEAN_CARBON_AC_PA )
      patch$hflxca       [[plab]] = - mymont$MMEAN_SENSIBLE_AC_PA
      patch$hflxgc       [[plab]] =   mymont$MMEAN_SENSIBLE_GC_PA
      patch$qwflxca      [[plab]] = - mymont$MMEAN_VAPOR_AC_PA    * mmean.can.alvli.pa
      patch$wflxca       [[plab]] = - mymont$MMEAN_VAPOR_AC_PA    * day.sec
      patch$wflxgc       [[plab]] =   mymont$MMEAN_VAPOR_GC_PA    * day.sec
      patch$ustar        [[plab]] =   mymont$MMEAN_USTAR_PA
      patch$albedo       [[plab]] =   mymont$MMEAN_ALBEDO_PA
      patch$rshortup     [[plab]] =   mymont$MMEAN_RSHORTUP_PA
      patch$rlongup      [[plab]] =   mymont$MMEAN_RLONGUP_PA
      patch$parup        [[plab]] =   mymont$MMEAN_PARUP_PA       * Watts.2.Ein * 1e6
      patch$rnet         [[plab]] =   mymont$MMEAN_RNET_PA
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Find the patch-level properties that are derived from cohort-level.            #
      #------------------------------------------------------------------------------------#
      patch$lai          [[plab]] = rep(0.,times=mymont$NPATCHES_GLOBAL)
      patch$wai          [[plab]] = rep(0.,times=mymont$NPATCHES_GLOBAL)
      patch$agb          [[plab]] = rep(0.,times=mymont$NPATCHES_GLOBAL)
      patch$ba           [[plab]] = rep(0.,times=mymont$NPATCHES_GLOBAL)
      patch$wood.dens    [[plab]] = rep(NA,times=mymont$NPATCHES_GLOBAL)
      patch$can.depth    [[plab]] = rep(0.,times=mymont$NPATCHES_GLOBAL)
      patch$can.area     [[plab]] = rep(0.,times=mymont$NPATCHES_GLOBAL)
      patch$leaf.temp    [[plab]] = mymont$MMEAN_CAN_TEMP_PA  - t00
      patch$leaf.vpd     [[plab]] = mymont$MMEAN_CAN_VPDEF_PA * 0.01
      patch$wood.temp    [[plab]] = mymont$MMEAN_CAN_TEMP_PA  - t00
      patch$gpp          [[plab]] = rep(0.,times=mymont$NPATCHES_GLOBAL)
      patch$npp          [[plab]] = rep(0.,times=mymont$NPATCHES_GLOBAL)
      patch$plant.resp   [[plab]] = rep(0.,times=mymont$NPATCHES_GLOBAL)
      patch$leaf.temp    [[plab]] = mymont$MMEAN_CAN_TEMP_PA - t00
      patch$hflxlc       [[plab]] = rep(0.,times=mymont$NPATCHES_GLOBAL)
      patch$hflxwc       [[plab]] = rep(0.,times=mymont$NPATCHES_GLOBAL)
      patch$wflxlc       [[plab]] = rep(0.,times=mymont$NPATCHES_GLOBAL)
      patch$wflxwc       [[plab]] = rep(0.,times=mymont$NPATCHES_GLOBAL)
      patch$transp       [[plab]] = rep(0.,times=mymont$NPATCHES_GLOBAL)
      patch$soil.resp    [[plab]] = mymont$MMEAN_RH_PA
      patch$fast.soil.c  [[plab]] = mymont$MMEAN_FAST_SOIL_C
      patch$slow.soil.c  [[plab]] = mymont$MMEAN_SLOW_SOIL_C
      patch$struct.soil.c[[plab]] = mymont$MMEAN_STRUCT_SOIL_C


      if (any(ncohorts >0)){
         #----- Find some auxiliary patch-level properties. -------------------------------#
         lai.pa         = tapply( X     = mymont$MMEAN_LAI_CO
                                , INDEX = ipaconow
                                , FUN   = sum
                                )#end tapply
         wai.pa         = tapply( X     = mymont$WAI_CO
                                , INDEX = ipaconow
                                , FUN   = sum
                                )#end tapply
         leaf.energy.pa = tapply( X     = mymont$MMEAN_LEAF_ENERGY_CO
                                , INDEX = ipaconow
                                , FUN   = sum
                                )#end tapply
         leaf.water.pa  = tapply( X     = mymont$MMEAN_LEAF_WATER_CO
                                , INDEX = ipaconow
                                , FUN   = sum
                                )#end tapply
         leaf.hcap.pa   = tapply( X     = mymont$MMEAN_LEAF_HCAP_CO
                                , INDEX = ipaconow
                                , FUN   = sum
                                )#end tapply
         wood.energy.pa = tapply( X     = mymont$MMEAN_WOOD_ENERGY_CO
                                , INDEX = ipaconow
                                , FUN   = sum
                                )#end tapply
         wood.water.pa  = tapply( X     = mymont$MMEAN_WOOD_WATER_CO
                                , INDEX = ipaconow
                                , FUN   = sum
                                )#end tapply
         wood.hcap.pa   = tapply( X     = mymont$MMEAN_WOOD_HCAP_CO
                                , INDEX = ipaconow
                                , FUN   = sum
                                )#end tapply
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Load some cohort-level structures that we will use multiple times.          #
         #---------------------------------------------------------------------------------#
         dbhconow          = mymont$DBH
         dbhconow.lastmon  = mymont$DBH * exp(-pmax(0,mymont$DLNDBH_DT/12))
         pftconow          = mymont$PFT
         nplantconow       = mymont$NPLANT
         heightconow       = mymont$HITE
         wood.densconow    = pft$rho[pftconow]
         baconow           = mymont$BA_CO
         agbconow          = mymont$AGB_CO
         laiconow          = mymont$MMEAN_LAI_CO
         waiconow          = mymont$WAI_CO
         caiconow          = pmin(1.,nplantconow * dbh2ca(dbh=dbhconow,ipft=pftconow))
         taiconow          = laiconow + waiconow
         gppconow          = mymont$MMEAN_GPP_CO
         nppconow          = mymont$MMEAN_NPP_CO
         plrespconow       = mymont$MMEAN_PLRESP_CO
         assim.lightconow  = mymont$MMEAN_A_LIGHT_CO
         assim.rubpconow   = mymont$MMEAN_A_RUBP_CO
         assim.co2conow    = mymont$MMEAN_A_CO2_CO
         #fabio used newly defined variables instead of mymont
         assim.ratioconow  = with(mymont, assim.lightconow
                                        / pmax(1e-6,pmin(assim.rubpconow,assim.co2conow)))


         #---------------------------------------------------------------------------------#
         #     Find soil respiration.  We must aggregate the components, some of which are #
         # cohort-level, whilst others are already patch level.                            #
         #---------------------------------------------------------------------------------#
            #----- Find biomass of some tissues. ------------------------------------------#
            bdeadconow        = mymont$BDEAD
            bleafconow        = mymont$MMEAN_BLEAF_CO
            bsapwoodconow     = mymont$BSAPWOODA+mymont$BSAPWOODB
            if (all(mymont$MMEAN_BROOT_CO == 0)){
               bfrootconow    = ( dbh2bl(dbh=dbhconow.lastmon,ipft=pftconow)
                                * pft$qroot[pftconow] )
            }else{
               bfrootconow    = mymont$MMEAN_BROOT_CO
            }#end if
            bcrootconow       = mymont$BSAPWOODB + (1. - pft$agf.bs[pftconow]) * bdeadconow
            bstemconow        = mymont$BSAPWOODA +       pft$agf.bs[pftconow]  * bdeadconow
            brootconow        = bfrootconow + bcrootconow
            baliveconow       = bleafconow + bfrootconow + bsapwoodconow
            bstorageconow     = mymont$MMEAN_BSTORAGE_CO
            bseedsconow       = mymont$BSEEDS_CO
            biomassconow      = baliveconow + bstorageconow + bseedsconow + bdeadconow
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #    Get storage and growth respiration, which will be distributed amongst     #
            # tissues.                                                                     #
            #------------------------------------------------------------------------------#
            #growth.respconow  = mymont$MMEAN.GROWTH.RESP.CO
            #storage.respconow = mymont$MMEAN_STORAGE_RESP_CO
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #      Find the fractions that go to each pool.                                #
            #------------------------------------------------------------------------------#
            fg.froot = bfrootconow / (baliveconow + bdeadconow )
            fg.croot = bcrootconow / (baliveconow + bdeadconow )
            fs.croot = bcrootconow / (bcrootconow + bstemconow )
            #------------------------------------------------------------------------------#



             #-----------------------------------------------------------------------------#
             #      Attribute respiration to the different pools.  Assuming that non-      #
             # structural carbon respiration occurs in the woody components.               #
             #-----------------------------------------------------------------------------#
             #froot.respconow = mymont$MMEAN.ROOT.RESP.CO    + fg.froot * growth.respconow
             #croot.respconow = fs.croot * storage.respconow + fg.croot * growth.respconow
             #root.respconow  = froot.respconow + croot.respconow
             #-----------------------------------------------------------------------------#





             #----- Make root respiration extensive. --------------------------------------#
             #root.resp.pa = tapply(X=root.respconow*nplantconow,INDEX=ipaconow,FUN=sum)
             #-----------------------------------------------------------------------------#
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
         agb.pa        = tapply( X= agbconow    * nplantconow, INDEX = ipaconow, FUN = sum)
         ba.pa         = tapply( X= baconow     * nplantconow, INDEX = ipaconow, FUN = sum)
         gpp.pa        = tapply( X= gppconow    * nplantconow, INDEX = ipaconow, FUN = sum)
         npp.pa        = tapply( X= nppconow    * nplantconow, INDEX = ipaconow, FUN = sum)
         plant.resp.pa = tapply( X= plrespconow * nplantconow, INDEX = ipaconow, FUN = sum)
         #---------------------------------------------------------------------------------#





         #----- Add the variables that are already extensive. -----------------------------#
         hflxlc.pa = tapply( X     = mymont$MMEAN_SENSIBLE_LC_CO
                           , INDEX = ipaconow
                           , FUN   = sum
                           )#end tapply
         hflxwc.pa = tapply( X     = mymont$MMEAN_SENSIBLE_WC_CO
                           , INDEX = ipaconow
                           , FUN   = sum
                           )#end tapply
         wflxlc.pa = tapply( X     = mymont$MMEAN_VAPOR_LC_CO  * day.sec
                           , INDEX = ipaconow
                           , FUN   = sum
                           )#end tapply
         wflxwc.pa = tapply( X     = mymont$MMEAN_VAPOR_WC_CO  * day.sec
                           , INDEX = ipaconow
                           , FUN   = sum
                           )#end tapply
         transp.pa = tapply( X     = mymont$MMEAN_TRANSP_CO    * day.sec
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
         #      Vapour pressure deficit is found using weighted averages (LAIs are the     #
         # weights).                                                                       #
         #---------------------------------------------------------------------------------#
         leaf.vpd.pa = mapply( FUN      = weighted.mean
                             , x        = split(mymont$MMEAN_LEAF_VPDEF_CO,ipaconow)
                             , w        = split(laiconow                  ,ipaconow)
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
         patch$agb       [[plab]][idx              ] = agb.pa
         patch$ba        [[plab]][idx              ] = ba.pa
         patch$can.depth [[plab]][idx              ] = can.depth.pa
         patch$can.area  [[plab]][idx              ] = can.area.pa
         patch$wood.dens [[plab]][idx              ] = wood.dens.pa
         patch$leaf.temp [[plab]][idx[! leaf.empty]] = leaf.temp.pa [! leaf.empty]
         patch$leaf.water[[plab]][idx[! leaf.empty]] = leaf.water.pa[! leaf.empty]
         patch$leaf.vpd  [[plab]][idx[! leaf.empty]] = leaf.vpd.pa  [! leaf.empty]
         patch$wood.temp [[plab]][idx[! wood.empty]] = wood.temp.pa [! wood.empty]
         patch$gpp       [[plab]][idx              ] = gpp.pa
         patch$npp       [[plab]][idx              ] = npp.pa
         patch$plant.resp[[plab]][idx              ] = plant.resp.pa
         patch$hflxlc    [[plab]][idx              ] = hflxlc.pa
         patch$hflxwc    [[plab]][idx              ] = hflxwc.pa
         patch$wflxlc    [[plab]][idx              ] = wflxlc.pa
         patch$wflxwc    [[plab]][idx              ] = wflxwc.pa
         patch$transp    [[plab]][idx              ] = transp.pa
         #------ Soil respiration mixes cohort (root) and patch (hetetrophic). ------------#
         #patch$soil.resp [[plab]][idx] = patch$soil.resp [[plab]][idx] + root.resp.pa
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
      one.cohort = sum(ncohorts) == 1
      one.patch  = sum(npatches) == 1
      #fabio have to comment part of this section
      if (any (ncohorts > 0)){
      #if (any (0 > 0)){

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
         dbhconow.lastmon = mymont$DBH * exp(-pmax(0,mymont$DLNDBH_DT/12))
         dbhconow.1ago    = mymont$DBH * exp(-pmax(0,mymont$DLNDBH_DT))
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
         wood.densconow    = pft$rho[pftconow]
         baconow           = mymont$BA_CO
         agbconow          = mymont$AGB_CO
         laiconow          = mymont$MMEAN_LAI_CO
         waiconow          = mymont$WAI_CO
         caiconow          = pmin(1.,nplantconow * dbh2ca(dbh=dbhconow,ipft=pftconow))
         taiconow          = laiconow + waiconow
         gppconow          = mymont$MMEAN_GPP_CO



         #---------------------------------------------------------------------------------#
         #     Find biomass of all tissues.                                                #
         #---------------------------------------------------------------------------------#
         bdeadconow        = mymont$BDEAD
         bleafconow        = mymont$MMEAN_BLEAF_CO
         bsapwoodconow     = mymont$BSAPWOODA+mymont$BSAPWOODB
         if (all(mymont$MMEAN_BROOT_CO == 0)){
            bfrootconow    = ( dbh2bl(dbh=dbhconow.lastmon,ipft=pftconow)
                             * pft$qroot[pftconow] )
         }else{
            bfrootconow    = mymont$MMEAN_BROOT_CO
         }#end if
         bcrootconow       = mymont$BSAPWOODB + (1. - pft$agf.bs[pftconow]) * bdeadconow
         bstemconow        = mymont$BSAPWOODA +       pft$agf.bs[pftconow]  * bdeadconow
         brootconow        = bfrootconow + bcrootconow
         baliveconow       = bleafconow + bfrootconow + bsapwoodconow
         bstorageconow     = mymont$MMEAN_BSTORAGE_CO
         bseedsconow       = mymont$BSEEDS_CO
         biomassconow      = baliveconow + bstorageconow + bseedsconow + bdeadconow
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #      Find total NPP then total autotrophic, growth and storage respiration.     #
         # The latter two will be distributed amongst tissues.                             #
         #---------------------------------------------------------------------------------#
         nppconow          = mymont$MMEAN_NPP_CO
         plant.respconow   = mymont$MMEAN_PLRESP_CO
         #growth.respconow  = mymont$MMEAN.GROWTH.RESP.CO
         #storage.respconow = mymont$MMEAN_STORAGE_RESP_CO
         #gs.respconow      = growth.respconow + storage.respconow
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Find the fractions that go to each pool.                                   #
         #---------------------------------------------------------------------------------#
         fg.leaf  = bleafconow  / ( baliveconow + bdeadconow )
         fg.stem  = bstemconow  / ( baliveconow + bdeadconow )
         fg.froot = bfrootconow / ( baliveconow + bdeadconow )
         fg.croot = bcrootconow / ( baliveconow + bdeadconow )
         fs.stem  = bstemconow  / ( bcrootconow + bstemconow )
         fs.croot = bcrootconow / ( bcrootconow + bstemconow )
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Attribute respiration to the different pools.  Assuming that non-          #
         # structural carbon respiration 
         #---------------------------------------------------------------------------------#
         #leaf.respconow  = mymont$MMEAN.LEAF.RESP.CO    + fg.leaf  * growth.respconow
         #stem.respconow  = fs.stem  * storage.respconow + fg.stem  * growth.respconow
         #froot.respconow = mymont$MMEAN.ROOT.RESP.CO    + fg.froot * growth.respconow
         #croot.respconow = fs.croot * storage.respconow + fg.croot * growth.respconow
         #root.respconow  = froot.respconow + croot.respconow
         #---------------------------------------------------------------------------------#



         #----- Flags to tell whether leaves and branchwood were resolvable. --------------#
         leaf.okconow      = ( is.finite(mymont$MMEAN_LEAF_HCAP_CO)
                             & mymont$MMEAN_LEAF_HCAP_CO >= pft$veg.hcap.min[pftconow] )
         wood.okconow      = ( is.finite(mymont$MMEAN_WOOD_HCAP_CO)
                             & mymont$MMEAN_WOOD_HCAP_CO >= pft$veg.hcap.min[pftconow] )
         #---------------------------------------------------------------------------------#


         ##fabio mmean_cb doesnt exist i added a _CO, also get a subscript out of bounds
         ## for the thismonth operation
         if (kludgecbal){
            cbaconow          = mymont$MMEAN_CB_CO - (mondays - 1) * mymont$MMEAN_BSTORAGE_CO
            #------------------------------------------------------------------------------#
            #     Temporary fix to correct the carbon balance_                             #
            #------------------------------------------------------------------------------#
            # if (one.cohort){
            #    cbalightconow     = ( mymont$CB_LIGHTMAX[thismonth] 
            #                        - (mondays - 1) * mymont$MMEAN_BSTORAGE_CO )
            #    cbamoistconow     = ( mymont$CB_MOISTMAX[thismonth]
            #                        - (mondays - 1) * mymont$MMEAN_BSTORAGE_CO )
            # }else{
            #    cbalightconow     = ( mymont$CB_LIGHTMAX[,thismonth] 
            #                        - (mondays - 1) * mymont$MMEAN_BSTORAGE_CO )
            #    cbamoistconow     = ( mymont$CB_MOISTMAX[,thismonth]
            #                        - (mondays - 1) * mymont$MMEAN_BSTORAGE_CO )
            # }#end if
            #------------------------------------------------------------------------------#
         }else{
            cbaconow             = mymont$MMEAN_CB_CO
            #------------------------------------------------------------------------------#
            #     Temporary fix to correct the carbon balance_                             #
            #------------------------------------------------------------------------------#
            # if (one.cohort){
            #    cbalightconow     = mymont$CB_LIGHTMAX[thismonth]
            #    cbamoistconow     = mymont$CB_MOISTMAX[thismonth]
            # }else{
            #    cbalightconow     = mymont$CB_LIGHTMAX[,thismonth]
            #    cbamoistconow     = mymont$CB_MOISTMAX[,thismonth]
            # }#end if
            #------------------------------------------------------------------------------#
         }#end if

         #cbamaxconow       = klight * cbalightconow + (1. - klight) * cbamoistconow
         cbarelconow       = mymont$CBR_BAR
         mcostconow        = ( mymont$MMEAN_LEAF_MAINTENANCE_CO
                             + mymont$MMEAN_ROOT_MAINTENANCE_CO ) * yr.day
         ldropconow        = mymont$MMEAN_LEAF_DROP_CO * yr.day
         sm.stressconow    = 1. - mymont$MMEAN_FS_OPEN_CO
         lightconow        = mymont$MMEAN_LIGHT_LEVEL_CO
         light.beamconow   = mymont$MMEAN_LIGHT_LEVEL_BEAM_CO
         light.diffconow   = mymont$MMEAN_LIGHT_LEVEL_DIFF_CO


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
         hflxlcconow       = mymont$MMEAN_SENSIBLE_LC_CO           / nplantconow
         wflxlcconow       = mymont$MMEAN_VAPOR_LC_CO    * day.sec / nplantconow
         transpconow       = mymont$MMEAN_TRANSP_CO      * day.sec / nplantconow
         i.hflxlcconow     = ifelse( leaf.okconow
                                   , mymont$MMEAN_SENSIBLE_LC_CO           / laiconow
                                   , NA
                                   )#end ifelse
         i.wflxlcconow     = ifelse( leaf.okconow
                                   , mymont$MMEAN_VAPOR_LC_CO    * day.sec / laiconow
                                   , NA
                                   )#end ifelse
         i.transpconow     = ifelse( leaf.okconow
                                   , mymont$MMEAN_TRANSP_CO      * day.sec / laiconow
                                   , NA
                                   )#end ifelse
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Find the leaf interstitial space and boundary layer specific humidities to  #
         # convert conductance to kgW/m2/day.                                              #
         #---------------------------------------------------------------------------------#
         lpsiconow      = ifelse( leaf.okconow
                                , mymont$MMEAN_TRANSP_CO/laiconow
                                , NA
                                )#end ifelse
         ##fabio MEAN_WFLXWC_CO doesn't exist in h5
         # wpsiconow      = ifelse( wood.okconow
         #                        , mymont$MMEAN_WFLXWC_CO/waiconow
         #                        , NA
         #                        )#end ifelse
         can.shv.conow  = rep(mymont$MMEAN_CAN_SHV_PA,times=ncohorts)
         #---- Net conductance, combining stomatal and boundary layer. --------------------#
         leaf.gnw.conow = ifelse( mymont$MMEAN_LEAF_GBW_CO+mymont$MMEAN_LEAF_GSW_CO>1.e-10
                                ,   mymont$MMEAN_LEAF_GBW_CO * mymont$MMEAN_LEAF_GSW_CO
                                / ( mymont$MMEAN_LEAF_GBW_CO + mymont$MMEAN_LEAF_GSW_CO )
                                , pmin(mymont$MMEAN_LEAF_GBW_CO,mymont$MMEAN_LEAF_GSW_CO)
                                )#end ifelse
         lbl.shv.conow  = can.shv.conow + lpsiconow / pmax(mymont$MMEAN_LEAF_GBW_CO,1.e-10)
         ##fabio see above
         #wbl.shv.conow  = can.shv.conow + wpsiconow / pmax(mymont$MMEAN_WOOD_GBW_CO,1.e-10)
         lis.shv.conow  = can.shv.conow + lpsiconow / pmax(leaf.gnw.conow,1.e-10)
         #---------------------------------------------------------------------------------#


         #----- Find the conductances in kgW/m2/day. --------------------------------------#
         leaf.gbwconow  = ( mymont$MMEAN_LEAF_GBW_CO * day.sec * ep
                          * (1 + epim1 * can.shv.conow) * (1 + epim1 * lbl.shv.conow) )
         leaf.gswconow  = ( mymont$MMEAN_LEAF_GSW_CO * day.sec * ep
                          * (1 + epim1 * lbl.shv.conow) * (1 + epim1 * lis.shv.conow) )
         ##fabio see above
         # wood.gbwconow  = ( mymont$MMEAN_WOOD_GBW_CO * day.sec * ep
         #                  * (1 + epim1 * can.shv.conow) * (1 + epim1 * wbl.shv.conow) )
         #---------------------------------------------------------------------------------#


         #----- Find the net radiation for leaves (in m2/leaf!). --------------------------#
         leaf.parconow    = ifelse( leaf.okconow
                                  , mymont$MMEAN_PAR_L_CO / laiconow *  Watts.2.Ein
                                  , NA
                                  )#end ifelse
         leaf.rshortconow = ifelse( leaf.okconow
                                  , mymont$MMEAN_RSHORT_L_CO / laiconow
                                  , NA
                                  )#end ifelse
         leaf.rlongconow  = ifelse( leaf.okconow
                                  , mymont$MMEAN_RLONG_L_CO  / laiconow
                                  , NA
                                  )#end ifelse
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Leaf/wood thermal properties.                                               #
         #---------------------------------------------------------------------------------#
         leaf.waterconow   = ifelse( leaf.okconow
                                   , mymont$MMEAN_LEAF_WATER_CO / laiconow
                                   , NA
                                   )#end ifelse
         leaf.tempconow    = ifelse( leaf.okconow
                                   , mymont$MMEAN_LEAF_TEMP_CO  - t00
                                   , NA
                                   )#end ifelse
         wood.tempconow    = ifelse( wood.okconow
                                   , mymont$MMEAN_WOOD_TEMP_CO  - t00
                                   , NA
                                   )#end ifelse
         leaf.vpdconow     = ifelse( leaf.okconow
                                   , mymont$MMEAN_LEAF_VPDEF_CO  * 0.01
                                   , NA
                                   )#end ifelse
         #---------------------------------------------------------------------------------#




         #------ Find the demand and supply by m2gnd. -------------------------------------#
         demandconow       = mymont$MMEAN_PSI_OPEN_CO     * laiconow * day.sec
         supplyconow       = mymont$MMEAN_WATER_SUPPLY_CO * day.sec
         #---------------------------------------------------------------------------------#


         #------ Find the demographic rates. ----------------------------------------------#
         if (one.cohort){
            mortconow    = sum(mymont$MMEAN_MORT_RATE_CO)
            mortconow    = max(0,mortconow)
         }else{
            mortconow    = try(rowSums(mymont$MMEAN_MORT_RATE_CO))
            if ("try-error" %in% is(mortconow)) browser()
            mortconow    = pmax(0,mortconow)
         }#end if
         ncbmortconow    = pmax(0,mymont$MMEAN_MORT_RATE_CO[,2])
         dimortconow     = pmax(0,mortconow - ncbmortconow)
         recruitconow    = mymont$RECRUIT_DBH
         growthconow     = pmax(0,mymont$DLNDBH_DT)
         agb.growthconow = pmax(0,mymont$DLNAGB_DT)
         bsa.growthconow = pmax(0,mymont$DLNBA_DT )
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


         #fabio last else problem logical subscript too long, added drop = FALSE
         # and changed the subscript but it didn't work because later it didn't match
         # the length of the product in phap.lparconow  / laiconow * Watts.2.Ein * 1.e6
         #----- Find some averages for photoperiod. ---------------------------------------#
         # if (polar.night){
         #    phap.lparconow     = NA + pftconow
         #    phap.ltempconow    = NA + pftconow
         #    phap.lwaterconow   = NA + pftconow
         #    phap.lvpdconow     = NA + pftconow
         #    phap.fs.openconow  = NA + pftconow
         #    phap.lpsiconow     = NA + pftconow
         #    phap.leaf.gbaconow = NA + pftconow
         #    phap.leaf.gsaconow = NA + pftconow
         #    phap.can.shv.conow = NA + pftconow
         # }else if (one.cohort){
         #    phap.lparconow     = mean(mymont$QMEAN_PAR_L_CO     [phap])
         #    phap.ltempconow    = mean(mymont$QMEAN_LEAF_TEMP_CO [phap])
         #    phap.lwaterconow   = mean(mymont$QMEAN_LEAF_WATER_CO[phap])
         #    phap.lvpdconow     = mean(mymont$QMEAN_LEAF_VPDEF_CO[phap])
         #    phap.fs.openconow  = mean(mymont$QMEAN_FS_OPEN_CO   [phap])
         #    phap.lpsiconow     = mean(mymont$QMEAN_TRANSP_CO    [phap])
         #    phap.leaf.gbaconow = mean(mymont$QMEAN_LEAF_GBW_CO  [phap])
         #    phap.leaf.gsaconow = mean(mymont$QMEAN_LEAF_GSW_CO  [phap])
         #    phap.can.shv.conow = rep(mean(mymont$QMEAN_CAN_SHV_PA[phap]),times=ncohorts)
         # }else{
         #    phap.lparconow     = rowMeans(mymont$QMEAN_PAR_L_CO     [phap, ,drop=FALSE])
         #    phap.ltempconow    = rowMeans(mymont$QMEAN_LEAF_TEMP_CO [phap, ,drop=FALSE])
         #    phap.lwaterconow   = rowMeans(mymont$QMEAN_LEAF_WATER_CO[phap, ,drop=FALSE])
         #    phap.lvpdconow     = rowMeans(mymont$QMEAN_LEAF_VPDEF_CO[phap, ,drop=FALSE])
         #    phap.fs.openconow  = rowMeans(mymont$QMEAN_FS_OPEN_CO   [phap, ,drop=FALSE])
         #    phap.lpsiconow     = rowMeans(mymont$QMEAN_TRANSP_CO    [phap, ,drop=FALSE])
         #    phap.leaf.gbaconow = rowMeans(mymont$QMEAN_LEAF_GBW_CO  [phap, ,drop=FALSE])
         #    phap.leaf.gsaconow = rowMeans(mymont$QMEAN_LEAF_GSW_CO  [phap, ,drop=FALSE])
         #    if (one.patch){
         #       phap.can.shv.conow = rep( x     = mean(mymont$QMEAN_CAN_SHV_PA[phap])
         #                               , times = ncohorts
         #                               )#end rep
         #    }else{
         #       phap.can.shv.conow = rep( x     = rowMeans(mymont$QMEAN_CAN_SHV_PA[,phap])
         #                               , times = ncohorts
         #                               )#end rep
         #    }#end if
         #    #------------------------------------------------------------------------------#
         # }#end if
         # phap.lparconow      = ifelse( leaf.okconow
         #                             , phap.lparconow  / laiconow * Watts.2.Ein * 1.e6
         #                             , NA
         #                             )#end ifelse
         # phap.ltempconow     = ifelse( leaf.okconow, phap.ltempconow - t00      , NA )
         # phap.lwaterconow    = ifelse( leaf.okconow, phap.lwaterconow / laiconow, NA )
         # phap.lvpdconow      = ifelse( leaf.okconow, phap.lvpdconow * 0.01      , NA )
         # phap.smsconow       = ifelse( leaf.okconow, 1 - phap.fs.openconow      , NA )
         # phap.lpsiconow      = ifelse( leaf.okconow, phap.lpsiconow / laiconow  , NA )
         # phap.leaf.gbaconow  = ifelse( leaf.okconow, phap.leaf.gbaconow         , NA )
         # phap.leaf.gsaconow  = ifelse( leaf.okconow, phap.leaf.gsaconow         , NA )
         # phap.can.shv.conow  = ifelse( leaf.okconow, phap.can.shv.conow         , NA )
         # #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Find the leaf interstitial space and boundary layer specific humidities to  #
         # convert conductance to kgW/m2/day.                                              #
         #---------------------------------------------------------------------------------#
         #---- Net conductance, combining stomatal and boundary layer. --------------------#
         # fine.cond          = is.finite(phap.leaf.gbaconow) & is.finite(phap.leaf.gsaconow)
         # enough.cond        = ( fine.cond
         #                      & ( phap.leaf.gbaconow + phap.leaf.gsaconow ) > 1.e-10 )
         # phap.leaf.gnaconow = ifelse( fine.cond
         #                            , ifelse( enough.cond
         #                                    , ( phap.leaf.gbaconow * phap.leaf.gsaconow )
         #                                    / ( phap.leaf.gbaconow + phap.leaf.gsaconow )
         #                                    , pmin(phap.leaf.gbaconow,phap.leaf.gsaconow)
         #                                    )#end ifelse
         #                            , NA
         #                            )#end ifelse
         # phap.lbl.shv.conow = ( phap.can.shv.conow
         #                      + phap.lpsiconow / pmax(phap.leaf.gbaconow,1.e-10) )
         # phap.lis.shv.conow = ( phap.can.shv.conow
         #                      + phap.lpsiconow / pmax(phap.leaf.gnaconow,1.e-10) )
         # #---------------------------------------------------------------------------------#


         # #----- Find the conductances in kgW/m2/day. --------------------------------------#
         # phap.lgbwconow  = ( phap.leaf.gbaconow * day.sec * ep
         #                   * (1 + epim1 * phap.can.shv.conow) 
         #                   * (1 + epim1 * phap.lbl.shv.conow)
         #                   )#end phap.lgbwconow
         # phap.lgswconow  = ( phap.leaf.gsaconow * day.sec * ep
         #                   * (1 + epim1 * phap.lbl.shv.conow) 
         #                   * (1 + epim1 * phap.lis.shv.conow)
         #                   )#end phap.lgswconow
         # #---------------------------------------------------------------------------------#




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
         leaf.rshortconow    = NA
         leaf.rlongconow     = NA
         rueconow            = NA
         opencanconow        = NA
         useconow            = NA
         phap.lparconow      = NA
         phap.ltempconow     = NA
         phap.lwaterconow    = NA
         phap.lvpdconow      = NA
         phap.sm.stressconow = NA
         phap.leaf.gbwconow  = NA
         phap.leaf.gswconow  = NA
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

         ##fabio had to comment several things in here because they either didnt
         ## have a variable because it was commented before or becuse the length
         ## of the two operators in the sum function didn't match
         ## the problem is that sel for some reason doesn't subset the same number
         ## of objects, in one variable subsets two elements in another four
         #----- Decide which PFT to use. --------------------------------------------------#
         for (p in sequence(npft+1)){
            sel.pft   = pftconow == p | p == (npft+1)
            sel       = sel.pft & sel.dbh
            if (any(sel)){
               #----- Extensive properties. -----------------------------------------------#
               szpft$lai         [m,d,p] = sum( laiconow          [sel]
                                              * areaconow         [sel]
                                              , na.rm = TRUE
                                              )#end if
               szpft$wai         [m,d,p] = sum( waiconow          [sel]
                                              * areaconow         [sel]
                                              , na.rm = TRUE
                                              )#end if
               szpft$tai         [m,d,p] = sum( taiconow          [sel]
                                              * areaconow         [sel]
                                              , na.rm = TRUE
                                              )#end if
               szpft$nplant      [m,d,p] = sum( nplantconow       [sel]
                                              * areaconow         [sel]
                                              , na.rm = TRUE
                                              )#end if
               szpft$demand      [m,d,p] = sum( demandconow       [sel]
                                              * areaconow         [sel]
                                              , na.rm = TRUE
                                              )#end if
               szpft$supply      [m,d,p] = sum( supplyconow       [sel]
                                              * areaconow         [sel]
                                              , na.rm = TRUE
                                              )#end if
               #----- Intensive properties, use nplant to make them extensive. ------------#
               szpft$agb         [m,d,p] = sum( w.nplant          [sel]
                                              * agbconow          [sel]
                                              , na.rm = TRUE
                                              )#end if
               szpft$biomass     [m,d,p] = sum( w.nplant          [sel]
                                              * biomassconow      [sel]
                                              , na.rm = TRUE
                                              )#end if
               szpft$ba          [m,d,p] = sum( w.nplant          [sel]
                                              * baconow           [sel]
                                              , na.rm = TRUE
                                              )#end if
               szpft$gpp         [m,d,p] = sum( w.nplant          [sel]
                                              * gppconow          [sel]
                                              , na.rm = TRUE
                                              )#end if
               szpft$npp         [m,d,p] = sum( w.nplant          [sel]
                                              * nppconow          [sel]
                                              , na.rm = TRUE
                                              )#end if
               szpft$mco         [m,d,p] = sum( w.nplant          [sel]
                                              * mcostconow        [sel]
                                              , na.rm = TRUE
                                              )#end if
               szpft$dcbadt      [m,d,p] = sum( w.nplant          [sel]
                                              * dcbadtconow       [sel]
                                              , na.rm = TRUE
                                              )#end if
               # szpft$cba         [m,d,p] = sum( w.nplant          [sel]
               #                                * cbaconow          [sel]
               #                                , na.rm = TRUE
               #                                )#end if
               # szpft$cbamax      [m,d,p] = sum( w.nplant          [sel]
               #                                * cbamaxconow       [sel]
               #                                , na.rm = TRUE
               #                                )#end if
               # szpft$cbalight    [m,d,p] = sum( w.nplant          [sel]
               #                                * cbalightconow     [sel]
               #                                , na.rm = TRUE
               #                                )#end if
               # szpft$cbamoist    [m,d,p] = sum( w.nplant          [sel]
               #                                * cbamoistconow     [sel]
               #                                , na.rm = TRUE
               #                                )#end if
               # szpft$ldrop       [m,d,p] = sum( w.nplant          [sel]
               #                                * ldropconow        [sel]
               #                                , na.rm = TRUE
               #                                )#end if
               # szpft$leaf.resp   [m,d,p] = sum( w.nplant          [sel]
               #                               * leaf.respconow    [sel]
               #                               , na.rm = TRUE
               #                               )#end if
               # szpft$stem.resp   [m,d,p] = sum( w.nplant          [sel]
               #                                * stem.respconow    [sel]
               #                                , na.rm = TRUE
               #                                )#end if
               #  # szpft$root.resp   [m,d,p] = sum( w.nplant          [sel]
                #                               * root.respconow    [sel]
                #                               , na.rm = TRUE
                #                               )#end if
                # szpft$froot.resp  [m,d,p] = sum( w.nplant          [sel]
                #                               * froot.respconow   [sel]
                #                               , na.rm = TRUE
                #                               )#end if
                # szpft$croot.resp  [m,d,p] = sum( w.nplant          [sel]
                #                               * croot.respconow   [sel]
                #                               , na.rm = TRUE
                #                               )#end if
               # szpft$growth.resp [m,d,p] = sum( w.nplant          [sel]
               #                                * growth.respconow  [sel]
               #                                , na.rm = TRUE
               #                                )#end if
               # szpft$storage.resp[m,d,p] = sum( w.nplant          [sel]
               #                                * storage.respconow [sel]
               #                                , na.rm = TRUE
               #                                )#end if
               # szpft$plant.resp  [m,d,p] = sum( w.nplant          [sel]
               #                                * plant.respconow   [sel]
               #                                , na.rm = TRUE
               #                                )#end if
               # szpft$bdead       [m,d,p] = sum( w.nplant          [sel]
               #                                * bdeadconow        [sel]
               #                                , na.rm = TRUE
               #                                )#end if
               # szpft$balive      [m,d,p] = sum( w.nplant          [sel]
               #                                * baliveconow       [sel]
               #                                , na.rm = TRUE
               #                                )#end if
               # szpft$bleaf       [m,d,p] = sum( w.nplant          [sel]
               #                                * bleafconow        [sel]
               #                                , na.rm = TRUE
               #                                )#end if
               # szpft$bstem       [m,d,p] = sum( w.nplant          [sel]
               #                                * bstemconow        [sel]
               #                                , na.rm = TRUE
               #                                )#end if
               # szpft$broot       [m,d,p] = sum( w.nplant          [sel]
               #                                * brootconow        [sel]
               #                                , na.rm = TRUE
               #                                )#end if
               # szpft$bfroot      [m,d,p] = sum( w.nplant          [sel]
               #                                * bfrootconow       [sel]
               #                                , na.rm = TRUE
               #                                )#end if
               # szpft$bcroot      [m,d,p] = sum( w.nplant          [sel]
               #                                * bcrootconow       [sel]
               #                                , na.rm = TRUE
               #                                )#end if
               # szpft$bsapwood    [m,d,p] = sum( w.nplant          [sel]
               #                                * bsapwoodconow     [sel]
               #                                , na.rm = TRUE
               #                                )#end if
               # szpft$bstorage    [m,d,p] = sum( w.nplant          [sel]
               #                                * bstorageconow     [sel]
               #                                , na.rm = TRUE
               #                                )#end if
               # szpft$bseeds      [m,d,p] = sum( w.nplant          [sel]
               #                                * bseedsconow       [sel]
               #                                , na.rm = TRUE
               #                                )#end if
               # szpft$hflxlc      [m,d,p] = sum( w.nplant          [sel]
               #                                * hflxlcconow       [sel]
               #                                , na.rm = TRUE
               #                                )#end if
               # szpft$wflxlc      [m,d,p] = sum( w.nplant          [sel]
               #                                * wflxlcconow       [sel]
               #                                , na.rm = TRUE
               #                                )#end if
               # szpft$transp      [m,d,p] = sum( w.nplant          [sel]
               #                                * transpconow       [sel]
               #                                , na.rm = TRUE
               #                                )#end if
               #---------------------------------------------------------------------------#



               #----- Leaf/wood intensive properties , weighted means using LAI/WAI. ------#
               szpft$sm.stress  [m,d,p] = weighted.mean( x     = sm.stressconow  [sel]
                                                       , w     = w.lai           [sel]
                                                       , na.rm = TRUE
                                                       )#end weighted.mean
               # szpft$phap.sms   [m,d,p] = weighted.mean( x     = phap.smsconow   [sel]
               #                                         , w     = w.lai           [sel]
               #                                         , na.rm = TRUE
               #                                         )#end weighted.mean
               szpft$leaf.par   [m,d,p] = weighted.mean( x     = leaf.parconow   [sel] 
                                                       , w     = w.lai           [sel]
                                                       , na.rm = TRUE
                                                       )#end weighted.mean
               # szpft$phap.lpar  [m,d,p] = weighted.mean( x     = phap.lparconow  [sel] 
               #                                         , w     = w.lai           [sel]
               #                                         , na.rm = TRUE
               #                                         )#end weighted.mean
               szpft$leaf.rshort[m,d,p] = weighted.mean( x     = leaf.rshortconow[sel] 
                                                       , w     = w.lai           [sel]
                                                       , na.rm = TRUE
                                                       )#end weighted.mean
               szpft$leaf.rlong [m,d,p] = weighted.mean( x     = leaf.rlongconow [sel] 
                                                       , w     = w.lai           [sel]
                                                       , na.rm = TRUE
                                                       )#end weighted.mean
               szpft$leaf.temp  [m,d,p] = weighted.mean( x     = leaf.tempconow  [sel]
                                                       , w     = w.lai           [sel]
                                                       , na.rm = TRUE
                                                       )#end weighted.mean
               # szpft$phap.ltemp [m,d,p] = weighted.mean( x     = phap.ltempconow [sel]
               #                                         , w     = w.lai           [sel]
               #                                         , na.rm = TRUE
               #                                         )#end weighted.mean
               szpft$leaf.water [m,d,p] = weighted.mean( x     = leaf.waterconow [sel]
                                                       , w     = w.lai           [sel]
                                                       , na.rm = TRUE
                                                       )#end weighted.mean
               # szpft$phap.lwater[m,d,p] = weighted.mean( x     = phap.lwaterconow[sel]
               #                                         , w     = w.lai           [sel]
               #                                         , na.rm = TRUE
               #                                         )#end weighted.mean
               szpft$wood.temp  [m,d,p] = weighted.mean( x     = wood.tempconow  [sel]
                                                       , w     = w.lai           [sel]
                                                       , na.rm = TRUE
                                                       )#end weighted.mean
               szpft$leaf.vpd   [m,d,p] = weighted.mean( x     = leaf.vpdconow   [sel]
                                                       , w     = w.lai           [sel]
                                                       , na.rm = TRUE
                                                       )#end weighted.mean
               # szpft$phap.lvpd  [m,d,p] = weighted.mean( x     = phap.lvpdconow  [sel]
               #                                         , w     = w.lai           [sel]
               #                                         , na.rm = TRUE
               #                                         )#end weighted.mean
               szpft$i.transp   [m,d,p] = weighted.mean( x     = i.transpconow   [sel]
                                                       , w     = w.lai           [sel]
                                                       , na.rm = TRUE
                                                       )#end weighted.mean
               szpft$i.wflxlc   [m,d,p] = weighted.mean( x     = i.wflxlcconow   [sel]
                                                       , w     = w.lai           [sel]
                                                       , na.rm = TRUE
                                                       )#end weighted.mean
               szpft$i.hflxlc   [m,d,p] = weighted.mean( x     = i.hflxlcconow   [sel]
                                                       , w     = w.lai           [sel]
                                                       , na.rm = TRUE
                                                       )#end weighted.mean
               szpft$leaf.gbw   [m,d,p] = weighted.mean( x     = leaf.gbwconow   [sel]
                                                       , w     = w.lai           [sel]
                                                       , na.rm = TRUE
                                                       )#end weighted.mean
               # szpft$phap.lgbw  [m,d,p] = weighted.mean( x     = phap.lgbwconow  [sel]
               #                                         , w     = w.lai           [sel]
               #                                         , na.rm = TRUE
               #                                         )#end weighted.mean
               szpft$leaf.gsw   [m,d,p] = weighted.mean( x     = leaf.gswconow   [sel]
                                                       , w     = w.lai           [sel]
                                                       , na.rm = TRUE
                                                       )#end weighted.mean
               # szpft$phap.lgsw  [m,d,p] = weighted.mean( x     = phap.lgswconow  [sel]
               #                                         , w     = w.lai           [sel]
               #                                         , na.rm = TRUE
               #                                         )#end weighted.mean
               szpft$assim.light[m,d,p] = weighted.mean( x     = assim.lightconow[sel]
                                                       , w     = w.lai           [sel]
                                                       , na.rm = TRUE
                                                       )#end weighted.mean
               szpft$assim.rubp [m,d,p] = weighted.mean( x     = assim.rubpconow [sel]
                                                       , w     = w.lai           [sel]
                                                       , na.rm = TRUE
                                                       )#end weighted.mean
               szpft$assim.co2  [m,d,p] = weighted.mean( x     = assim.co2conow  [sel]
                                                       , w     = w.lai           [sel]
                                                       , na.rm = TRUE
                                                       )#end weighted.mean
               # szpft$wood.gbw   [m,d,p] = weighted.mean( x     = wood.gbwconow   [sel]
               #                                         , w     = w.wai           [sel]
               #                                         , na.rm = TRUE
               #                                         )#end weighted.mean
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
               # szpft$i.cbamax     [m,d,p] = weighted.mean( x     = cbamaxconow      [sel]
               #                                           , w     = w.nplant         [sel]
               #                                           , na.rm = TRUE
               #                                           )#end weighted.mean
               # szpft$i.cbalight   [m,d,p] = weighted.mean( x     = cbalightconow    [sel]
               #                                           , w     = w.nplant         [sel]
               #                                           , na.rm = TRUE
               #                                           )#end weighted.mean
               # szpft$i.cbamoist   [m,d,p] = weighted.mean( x     = cbamoistconow    [sel]
               #                                           , w     = w.nplant         [sel]
               #                                           , na.rm = TRUE
               #                                           )#end weighted.mean
               # #---------------------------------------------------------------------------#


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
               szpft$growth     [m,d,p] = dbh.growth
               szpft$agb.growth [m,d,p] = agb.growth
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



            #------------------------------------------------------------------------------#
            #     Change of basic properties.                                              #
            #------------------------------------------------------------------------------#
            if (m == 1){
               szpft$change     [m,d,p] = 0.
               szpft$agb.change [m,d,p] = 0.
               szpft$bsa.change [m,d,p] = 0.
            }else{
               szpft$change     [m,d,p] = ( 12. * log( szpft$nplant[m  ,d,p]
                                                     / szpft$nplant[m-1,d,p] ) )
               szpft$agb.change [m,d,p] = ( 12. * log( szpft$agb   [m  ,d,p]
                                                     / szpft$agb   [m-1,d,p] ) )
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
      emean$leaf.rshort     [m] = szpft$leaf.rshort    [m,ndbh+1,npft+1]
      emean$leaf.rlong      [m] = szpft$leaf.rlong     [m,ndbh+1,npft+1]
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
         #cohort$leaf.resp    [[clab]] = leaf.respconow
         #cohort$stem.resp    [[clab]] = stem.respconow
         # cohort$root.resp    [[clab]] = root.respconow
         # cohort$froot.resp   [[clab]] = froot.respconow
         # cohort$croot.resp   [[clab]] = croot.respconow
         cohort$plant.resp   [[clab]] = plant.respconow
         cohort$assim.light  [[clab]] = assim.lightconow
         cohort$assim.rubp   [[clab]] = assim.rubpconow
         cohort$assim.co2    [[clab]] = assim.co2conow
         cohort$assim.ratio  [[clab]] = assim.ratioconow
         cohort$npp          [[clab]] = nppconow
         cohort$cba          [[clab]] = cbaconow
         #cohort$cbamax       [[clab]] = cbamaxconow
          # cohort$cbalight     [[clab]] = cbalightconow
         #cohort$cbamoist     [[clab]] = cbamoistconow
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
         cohort$leaf.rshort  [[clab]] = leaf.rshortconow
         cohort$leaf.rlong   [[clab]] = leaf.rlongconow
         cohort$rue          [[clab]] = rueconow
      } #end if month=sasmonth
      #------------------------------------------------------------------------------------#
    H5close()
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
