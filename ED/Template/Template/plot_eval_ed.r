#==========================================================================================#
#==========================================================================================#
#     Leave these commands at the beginning.  They will refresh the session.               #
#------------------------------------------------------------------------------------------#
rm(list=ls())
graphics.off()
options(warn=0)
gc()
#==========================================================================================#
#==========================================================================================#



#==========================================================================================#
#==========================================================================================#
#      Here is the user defined variable section.                                          #
#------------------------------------------------------------------------------------------#


#----- Paths. -----------------------------------------------------------------------------#
here           = "thispath"     # Current directory.
there          = "thatpath"     # Directory where analyses/history are
srcdir         = "thisrscpath"  # Source  directory.
outroot        = "thisoutroot"  # Directory for figures
#------------------------------------------------------------------------------------------#


#----- Time options. ----------------------------------------------------------------------#
monthbeg       = thismontha   # First month to use
yearbeg        = thisyeara    # First year to consider
yearend        = thisyearz    # Maximum year to consider
reload.data    = TRUE         # Should I reload partially loaded data?
sasmonth.short = c(2,5,8,11)  # Months for SAS plots (short runs)
sasmonth.long  = 5            # Months for SAS plots (long runs)
nyears.long    = 25           # Runs longer than this are considered long runs.
#------------------------------------------------------------------------------------------#



#----- Name of the simulations. -----------------------------------------------------------#
myplaces       = c("thispoly")
#------------------------------------------------------------------------------------------#



#----- Plot options. ----------------------------------------------------------------------#
outform        = thisoutform    # Formats for output file.  Supported formats are:
                                #   - "X11" - for printing on screen
                                #   - "eps" - for postscript printing
                                #   - "png" - for PNG printing
                                #   - "pdf" - for PDF printing

depth          = 96             # PNG resolution, in pixels per inch
paper          = "letter"       # Paper size, to define the plot shape
ptsz           = 14             # Font size.
lwidth         = 2.5            # Line width
plotgrid       = TRUE           # Should I plot the grid in the background? 

legwhere       = "topleft"      # Where should I place the legend?
inset          = 0.01           # inset distance between legend and edge of plot region.
scalleg        = 0.20
cex.main       = 0.8            # Scale coefficient for the title
ibackground    = mybackground   # Background settings (check load_everything.r)
#------------------------------------------------------------------------------------------#




#------ Settings for the model assessment. ------------------------------------------------#
use.distrib    = "mydistrib"    # Which distribution to use.  Valid options are:
                                #   - "norm" -- normal distribution
                                #   - "sn"   -- skewed normal distribution
                                #   - "edf"  -- empirical distribution function
hourblock.len  = 3              # Length of the time blocks, in hours
                                #
n.quant        = 1024           # Number of quantiles to produce the density function.
                                #    We strongly advise to choose a number that is a power
                                #    of two, especially when using EDF (otherwise 
                                #    distributions will be interpolated).
#------------------------------------------------------------------------------------------#



#------ Reload the data? ------------------------------------------------------------------#
reload.data    = TRUE           # Reload data?
#------------------------------------------------------------------------------------------#




#------ Miscellaneous settings. -----------------------------------------------------------#
slz.min             = -5.0         # The deepest depth that trees access water.
idbh.type           = myidbhtype   # Type of DBH class
                                   # 1 -- Every 10 cm until 100cm; > 100cm
                                   # 2 -- 0-10; 10-20; 20-35; 35-50; 50-70; > 70 (cm)
                                   # 3 -- 0-10; 10-35; 35-55; > 55 (cm)
corr.growth.storage = mycorrection # Correction factor to be applied to growth and
                                   #   storage respiration
#------------------------------------------------------------------------------------------#





#------------------------------------------------------------------------------------------#
#     Eddy flux comparisons.                                                               #
#------------------------------------------------------------------------------------------#
n            = 0
compvar      = list()
n            = n + 1
compvar[[n]] = list( vnam       = "hflxca"
                   , desc       = "Sensible heat flux"
                   , unit       = "wom2"
                   , col.obser  = "#404040"
                   , col.model  = "#F87856"
                   , leg.corner = "topleft"
                   , sunvar     = FALSE
                   )#end list
n            = n + 1
compvar[[n]] = list( vnam       = "wflxca"
                   , desc       = "Water vapour flux"
                   , unit       = "kgwom2oday"
                   , col.obser  = "#404040"
                   , col.model  = "#1BA2F7"
                   , leg.corner = "topleft"
                   , sunvar     = FALSE
                   )#end list
n            = n + 1
compvar[[n]] = list( vnam       = "cflxca"
                   , desc       = "Carbon dioxide flux"
                   , unit       = "umolcom2os"
                   , col.obser  = "#404040"
                   , col.model  = "#2BD2DB"
                   , leg.corner = "bottomright"
                   , sunvar     = FALSE
                   )#end list
n            = n + 1
compvar[[n]] = list( vnam       = "cflxst"
                   , desc       = "Carbon dioxide storage"
                   , unit       = "umolcom2os"
                   , col.obser  = "#404040"
                   , col.model  = "#0E6E81"
                   , leg.corner = "topleft"
                   , sunvar     = FALSE
                   )#end list
n            = n + 1
compvar[[n]] = list( vnam       = "gpp"
                   , desc       = "Gross primary productivity"
                   , unit       = "kgcom2oyr"
                   , col.obser  = "#404040"
                   , col.model  = "#0E6E81"
                   , leg.corner = "topleft"
                   , sunvar     = TRUE
                   )#end list
n            = n + 1
compvar[[n]] = list( vnam       = "reco"
                   , desc       = "Ecosystem respiration"
                   , unit       = "kgcom2oyr"
                   , col.obser  = "#404040"
                   , col.model  = "#F87856"
                   , leg.corner = "topleft"
                   , sunvar     = FALSE
                   )#end list
n            = n + 1
compvar[[n]] = list( vnam       = "nep"
                   , desc       = "Net ecosystem productivity"
                   , unit       = "kgcom2oyr"
                   , col.obser  = "#404040"
                   , col.model  = "#2BD2DB"
                   , leg.corner = "topleft"
                   , sunvar     = FALSE
                   )#end list
n            = n + 1
compvar[[n]] = list( vnam       = "nee"
                   , desc       = "Net ecosystem exchange"
                   , unit       = "umolcom2os"
                   , col.obser  = "#404040"
                   , col.model  = "#2BD2DB"
                   , leg.corner = "bottomright"
                   , sunvar     = FALSE
                   )#end list
n            = n + 1
compvar[[n]] = list( vnam       = "ustar"
                   , desc       = "Friction velocity"
                   , unit       = "mos"
                   , col.obser  = "#404040"
                   , col.model  = "#811F9E"
                   , leg.corner = "topleft"
                   , sunvar     = FALSE
                   )#end list
n            = n + 1
compvar[[n]] = list( vnam       = "rlongup"
                   , desc       = "Outgoing longwave radiation"
                   , unit       = "wom2"
                   , col.obser  = "#404040"
                   , col.model  = "#880D32"
                   , leg.corner = "topleft"
                   , sunvar     = FALSE
                   )#end list
n            = n + 1
compvar[[n]] = list( vnam       = "rnet"
                   , desc       = "Net radiation"
                   , unit       = "wom2"
                   , col.obser  = "#404040"
                   , col.model  = "#F87856"
                   , leg.corner = "topleft"
                   , sunvar     = FALSE
                   )#end list
n            = n + 1
compvar[[n]] = list( vnam       = "albedo"
                   , desc       = "Albedo"
                   , unit       = "empty"
                   , col.obser  = "#404040"
                   , col.model  = "#F87856"
                   , leg.corner = "topleft"
                   , sunvar     = TRUE
                   )#end list
n            = n + 1
compvar[[n]] = list( vnam       = "parup"
                   , desc       = "Outgoing PAR"
                   , unit       = "umolom2os"
                   , col.obser  = "#404040"
                   , col.model  = "#2BD2DB"
                   , leg.corner = "topleft"
                   , sunvar     = TRUE
                   )#end list
n            = n + 1
compvar[[n]] = list( vnam       = "rshortup"
                   , desc       = "Outgoing SW"
                   , unit       = "wom2"
                   , col.obser  = "#404040"
                   , col.model  = "#F87856"
                   , leg.corner = "topleft"
                   , sunvar     = TRUE
                   )#end list
n            = n + 1
compvar[[n]] = list( vnam       = "can.tmp"
                   , desc       = "CAS temperature"
                   , unit       = "degC"
                   , col.obser  = "#404040"
                   , col.model  = "#F87856"
                   , leg.corner = "topleft"
                   , sunvar     = FALSE
                   )#end list
n            = n + 1
compvar[[n]] = list( vnam       = "can.shv"
                   , desc       = "CAS specific humidity"
                   , unit       = "gwokg"
                   , col.obser  = "#404040"
                   , col.model  = "#1BA2F7"
                   , leg.corner = "topleft"
                   , sunvar     = FALSE
                   )#end list
n            = n + 1
compvar[[n]] = list( vnam       = "can.co2"
                   , desc       = "CAS CO2"
                   , unit       = "umolcomol"
                   , col.obser  = "#404040"
                   , col.model  = "#2BD2DB"
                   , leg.corner = "topleft"
                   , sunvar     = FALSE
                   )#end list
#------------------------------------------------------------------------------------------#



#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#      NO NEED TO CHANGE ANYTHING BEYOND THIS POINT UNLESS YOU ARE DEVELOPING THE CODE...  #
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#



#----- Load some packages and scripts. ----------------------------------------------------#
source(file.path(srcdir,"load.everything.r"))
#------------------------------------------------------------------------------------------#


#----- Set how many formats we must output. -----------------------------------------------#
outform = tolower(outform)
nout = length(outform)
#------------------------------------------------------------------------------------------#


#----- Set how many variables we will compare. --------------------------------------------#
compvar  = list.2.data.frame(compvar)
ncompvar = length(compvar)
#------------------------------------------------------------------------------------------#



#----- Load observations. -----------------------------------------------------------------#
obser.file = file.path(srcdir,"LBA_MIP.nogapfill.RData")
load(file=obser.file)
#------------------------------------------------------------------------------------------#



#----- Load census data. ------------------------------------------------------------------#
census.file = file.path(srcdir,"LBA_MIP.census_summ.RData")
load(file=census.file)
#------------------------------------------------------------------------------------------#



#----- Define plot window size ------------------------------------------------------------#
size = plotsize(proje=FALSE,paper=paper)
#------------------------------------------------------------------------------------------#



#---- Create the main output directory in case there is none. -----------------------------#
if (! file.exists(outroot)) dir.create(outroot)
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#     Big place loop starts here...                                                        #
#------------------------------------------------------------------------------------------#
for (place in myplaces){

   #----- Retrieve default information about this place and set up some variables. --------#
   thispoi   = locations(where=place,here=there,yearbeg=yearbeg,yearend=yearend
                        ,monthbeg=monthbeg)
   inpref    = thispoi$pathin
   outmain   = file.path(outroot,place)
   outpref   = file.path(outmain,"eval_ed")
   lieu      = thispoi$lieu
   lon       = thispoi$lon
   lat       = thispoi$lat
   iata      = thispoi$iata
   suffix    = thispoi$iata
   yeara     = thispoi$yeara
   yearz     = thispoi$yearz
   meszz     = thispoi$monz

   eft.yeara = poilist$yeara[match(iata,poilist$iata)]+5
   eft.yearz = poilist$yearz[match(iata,poilist$iata)]-1

   #----- Find the observations for this particular site. ---------------------------------#
   if (iata %in% c("mao","bdf")){
      obs.name = "obs.m34"
   }else if(iata %in% "stm"){
      obs.name = "obs.s67"
   }else if(iata %in% "rao"){
      obs.name = "obs.pdg"
   }else if(iata %in% "jpr"){
      obs.name = "obs.fns"
   }else if(iata %in% "btr"){
      obs.name = "obs.s77"
   }else{
      obs.name = paste0("obs.",iata)
   }#end if
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     We only run this part of the code if there are observations to compare with the   #
   # model.                                                                                #
   #---------------------------------------------------------------------------------------#
   if (obs.name %in% ls()){


      #----- Print a banner to entretain the user. ----------------------------------------#
      cat0(" + Evaluate model for ",lieu,".")

      #----- Get the observed variables. --------------------------------------------------#
      obser  = get(obs.name)
      ntimes = length(obser$when)
      keep   = obser$year %in% seq(from=eft.yeara,to=eft.yearz)
      if ("nzg" %in% names(obser)){
         obs.nzg = obser$nzg
      }else{
         obs.nzg = Inf
      }#end if
      if ("ust.filter" %in% names(obser)){
         nustar  = length(obser$ust.filter)
      }else{
         nustar  = Inf
      }#end if
      for (vname in names(obser)){
         if (length(obser[[vname]]) == ntimes){
            obser[[vname]] = obser[[vname]][keep]
         }else if (length(obser[[vname]]) %in%  (ntimes*c(obs.nzg,nustar))){
            obser[[vname]] = obser[[vname]][keep,]
         }else if ( ! ( length(obser[[vname]]) %in% c(0,1,obs.nzg,nustar) ) ){
            stay = numyears(obser[[vname]]) %in% seq(from=eft.yeara,to=eft.yearz)
            obser[[vname]] = obser[[vname]][stay]
         }#end if
      }#end for
      ntimes = length(obser$when)
      #------------------------------------------------------------------------------------#



      #----- Find some additional time info. ----------------------------------------------#
      obser        = alltimes(datin=obser,lon=lon,lat=lat,ed21=TRUE,zeronight=FALSE
                             ,meanval=TRUE,imetavg=1,nmean=12,na.rm=TRUE)
      obser$hr.idx = period.day(obser$when,dtblock=hourblock.len)
      obser$yr.idx = season(obser$when,add.year=FALSE)
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #      Make the RData file name, then we check whether we must read the files again  #
      # or use the stored RData.                                                           #
      #------------------------------------------------------------------------------------#
      path.data  = file.path(here,place,"rdata_hour")
      if (! file.exists(path.data)) dir.create(path.data)
      ed22.rdata = file.path(path.data,paste0(place,".RData"))
      if (reload.data && file.exists(ed22.rdata)){
         #----- Load the modelled dataset. ------------------------------------------------#
         cat0("   - Load previous session.")
         load(ed22.rdata)
         if ((! "eddy.complete" %in% ls()) && "complete" %in% ls()) eddy.complete = complete
         if ((! "eddy.tresume"  %in% ls()) && "tresume"  %in% ls()) eddy.tresume  = tresume
      }else{
         cat0("   - Start new session.")
         eddy.tresume    = 1
         eddy.complete   = FALSE
      }#end if
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     In case the met driver has become longer, we may need to expand the model      #
      # structure and continue reading.                                                    #
      #------------------------------------------------------------------------------------#
      if (eddy.tresume > 1){
         eddy.extend   = length(model$when) < ntimes
         eddy.complete = eddy.complete && ! eddy.extend
      }else{
         eddy.extend   = FALSE
      }#end if
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     In case we must extend the list...                                             #
      #------------------------------------------------------------------------------------#
      if (eddy.extend){
         cat0("   - Observed met driver became longer, extend model structure.")
         partial           = model
         npartial          = length(partial$when)
         rm(model)
         model             = list()
         model$when        = obser$when
         model$hr.idx      = period.day(model$when,dtblock=hourblock.len)
         model$yr.idx      = season(model$when,add.year=FALSE)
         model$nzg         = partial$nzg
         model$slz         = partial$slz
         model$soil        = partial$soil
         na.pad.vec        = rep(x=NA,times=ntimes-npartial)
         na.pad.mat        = matrix(data=NA,nrow=ntimes-npartial,ncol=model$nzg)
         model$atm.tmp     = c    (partial$atm.tmp     ,na.pad.vec)
         model$atm.shv     = c    (partial$atm.shv     ,na.pad.vec)
         model$atm.prss    = c    (partial$atm.prss    ,na.pad.vec)
         model$rain        = c    (partial$rain        ,na.pad.vec)
         model$atm.co2     = c    (partial$atm.co2     ,na.pad.vec)
         model$atm.vels    = c    (partial$atm.vels    ,na.pad.vec)
         model$atm.vpdef   = c    (partial$atm.vpdef   ,na.pad.vec)
         model$can.tmp     = c    (partial$can.tmp     ,na.pad.vec)
         model$can.shv     = c    (partial$can.shv     ,na.pad.vec)
         model$can.co2     = c    (partial$can.co2     ,na.pad.vec)
         model$rshort      = c    (partial$rshort      ,na.pad.vec)
         model$rlong       = c    (partial$rlong       ,na.pad.vec)
         model$par         = c    (partial$par         ,na.pad.vec)
         model$hflxca      = c    (partial$hflxca      ,na.pad.vec)
         model$wflxca      = c    (partial$wflxca      ,na.pad.vec)
         model$cflxca      = c    (partial$cflxca      ,na.pad.vec)
         model$cflxst      = c    (partial$cflxst      ,na.pad.vec)
         model$gpp         = c    (partial$gpp         ,na.pad.vec)
         model$reco        = c    (partial$reco        ,na.pad.vec)
         model$nep         = c    (partial$nep         ,na.pad.vec)
         model$nee         = c    (partial$nee         ,na.pad.vec)
         model$ustar       = c    (partial$ustar       ,na.pad.vec)
         model$rnet        = c    (partial$rnet        ,na.pad.vec)
         model$rlongup     = c    (partial$rlongup     ,na.pad.vec)
         model$albedo      = c    (partial$albedo      ,na.pad.vec)
         model$parup       = c    (partial$parup       ,na.pad.vec)
         model$rshortup    = c    (partial$rshortup    ,na.pad.vec)
         model$soil.tmp    = rbind(partial$soil.tmp    ,na.pad.mat)
         model$soil.water  = rbind(partial$soil.water  ,na.pad.mat)
         model$soil.matpot = rbind(partial$soil.matpot ,na.pad.mat)

         #---------------------------------------------------------------------------------#
         #     Sanity check.                                                               #
         #---------------------------------------------------------------------------------#
         if (length(model$atm.tmp) != length(obser$atm.temp)){
            cat0(" Mismatch in lengths!!!")
            cat0(" Length (obser):   ",length(obser$atm.temp),".")
            cat0(" Length (model):   ",length(model$atm.tmp ),".")
            cat0(" Length (partial): ",length(partial$when)  ,".")
            cat0(" Length (na.pad):  ",length(na.pad)        ,".")
            cat0(" TRESUME:          ",eddy.tresume          ,".")
            cat0(" NTIMES:           ",ntimes                ,".")
            stop("Structures should match, check your code.")
         }else{
            rm(partial,na.pad)
         }#end if
         #---------------------------------------------------------------------------------#
      }#end if
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #      Read data.                                                                    #
      #------------------------------------------------------------------------------------#
      if (! eddy.complete){

         #---------------------------------------------------------------------------------#
         #      Load data for all times.                                                   #
         #---------------------------------------------------------------------------------#
         last.cday   = "00"
         last.cmonth = sprintf("%4.4i",numyears (obser$when[eddy.tresume]))
         last.cyear  = sprintf("%4.4i",numyears (obser$when[eddy.tresume]))
         cat0("   - Read in files.")

         loop.times = seq(from=eddy.tresume,to=ntimes,by=1)
         for (tt in loop.times){
            cyear  = sprintf("%4.4i",numyears (obser$when[tt]))
            cmonth = sprintf("%2.2i",nummonths(obser$when[tt]))
            cday   = sprintf("%2.2i",numdays  (obser$when[tt]))
            chour  = sprintf("%2.2i",hours    (obser$when[tt]))
            cminu  = sprintf("%2.2i",minutes  (obser$when[tt]))
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #      Read the data in case the file exists.  In case it doesn't, save what   #
            # has been already read and quit.                                              #
            #------------------------------------------------------------------------------#
            h5file       = paste0(inpref,"-I-",cyear,"-",cmonth,"-",cday,"-",chour,cminu
                                 ,"00-g01.h5")
            h5file.bz2   = paste0(h5file,".bz2")
            h5file.gz    = paste0(h5file,".gz" )
            if (file.exists(h5file)){
               myinst    = hdf5load(file=h5file,load=FALSE,verbosity=0,tidy=TRUE)

            }else if(file.exists(h5file.bz2)){
               temp.file = file.path(tempdir(),basename(h5file))
               dummy     = bunzip2(filename=h5file.bz2,destname=temp.file,remove=FALSE)
               myinst    = hdf5load(file=temp.file,load=FALSE,verbosity=0,tidy=TRUE)
               dummy     = file.remove(temp.file)

            }else if(file.exists(h5file.gz)){
               temp.file = file.path(tempdir(),basename(h5file))
               dummy     = gunzip(filename=h5file.gz,destname=temp.file,remove=FALSE)
               myinst    = hdf5load(file=temp.file,load=FALSE,verbosity=0,tidy=TRUE)
               dummy     = file.remove(temp.file)
            }else{
               eddy.complete = FALSE
               eddy.tresume  = tt
               cat0("   - Simulation not ready yet.  Save partial ED-2.2 data to "
                   ,basename(ed22.rdata),".")
               save( list = c("model","eddy.complete","eddy.tresume")
                   , file = ed22.rdata)
               cat0("Quit.")
               q("no")
            }#end if
            #------------------------------------------------------------------------------#



            #----- Initialise the model structure. ----------------------------------------#
            if (tt == 1){
               model        = list()
               model$when   = obser$when
               model$hr.idx = period.day(model$when,dtblock=hourblock.len)
               model$yr.idx = season(model$when,add.year=FALSE)
               model$nzg    = myinst$NZG
               model$slz    = myinst$SLZ


               #---------------------------------------------------------------------------#
               #      Save the soil properties into a data.frame.                          #
               #---------------------------------------------------------------------------#
               ntext    = myinst$NTEXT.SOIL
               oargs    = list( isoilflg = myinst$ISOILFLG
                              , slxsand  = myinst$SLXSAND
                              , slxclay  = myinst$SLXCLAY
                              )#end list
               model$soil = data.frame(lapply(X  =apply(X     =mapply(FUN     =soil.params
                                                                     ,ntext   =ntext
                                                                     ,MoreArgs=oargs
                                                                     )#end mapply
                                                       ,MARGIN=1
                                                       ,FUN   =list
                                                       )#end apply
                                             ,FUN=unlist
                                             )#end lapply
                                      ,stringsAsFactors = FALSE
                                      )#end data.frame
               #---------------------------------------------------------------------------#

               empty.vec         = rep(NA,times=ntimes)
               empty.mat         = matrix(data=NA,nrow=ntimes,ncol=model$nzg)
               model$atm.tmp     = empty.vec
               model$atm.shv     = empty.vec
               model$atm.prss    = empty.vec
               model$rain        = empty.vec
               model$atm.co2     = empty.vec
               model$atm.vels    = empty.vec
               model$atm.vpdef   = empty.vec
               model$can.tmp     = empty.vec
               model$can.shv     = empty.vec
               model$can.co2     = empty.vec
               model$rshort      = empty.vec
               model$rlong       = empty.vec
               model$par         = empty.vec
               model$hflxca      = empty.vec
               model$wflxca      = empty.vec
               model$cflxca      = empty.vec
               model$cflxst      = empty.vec
               model$gpp         = empty.vec
               model$reco        = empty.vec
               model$nep         = empty.vec
               model$nee         = empty.vec
               model$ustar       = empty.vec
               model$rnet        = empty.vec
               model$rlongup     = empty.vec
               model$albedo      = empty.vec
               model$parup       = empty.vec
               model$rshortup    = empty.vec
               model$soil.tmp    = empty.mat
               model$soil.water  = empty.mat
               model$soil.matpot = empty.mat
            }#end if
            #------------------------------------------------------------------------------#







            if (last.cday != cday) cat0("     * ",basename(h5file),".")
            last.cday   = cday
            #------------------------------------------------------------------------------#
            model$atm.tmp    [tt ] =   myinst$FMEAN.ATM.TEMP.PY       - t00
            model$atm.shv    [tt ] =   myinst$FMEAN.ATM.SHV.PY        * 1000.
            model$atm.prss   [tt ] =   myinst$FMEAN.ATM.PRSS.PY       * 0.01
            model$rain       [tt ] =   myinst$FMEAN.PCPG.PY           * hr.sec
            model$atm.co2    [tt ] =   myinst$FMEAN.ATM.CO2.PY
            model$atm.vels   [tt ] =   myinst$FMEAN.ATM.VELS.PY
            model$atm.vpdef  [tt ] =   myinst$FMEAN.ATM.VPDEF.PY      * 0.01
            model$can.tmp    [tt ] =   myinst$FMEAN.CAN.TEMP.PY       - t00
            model$can.shv    [tt ] =   myinst$FMEAN.CAN.SHV.PY        * 1000.
            model$can.co2    [tt ] =   myinst$FMEAN.CAN.CO2.PY
            model$rshort     [tt ] =   myinst$FMEAN.ATM.RSHORT.PY
            model$rlong      [tt ] =   myinst$FMEAN.ATM.RLONG.PY
            model$par        [tt ] =   myinst$FMEAN.ATM.PAR.PY        * Watts.2.Ein * 1.e6
            model$hflxca     [tt ] = - myinst$FMEAN.SENSIBLE.AC.PY
            model$wflxca     [tt ] = - myinst$FMEAN.VAPOR.AC.PY       * day.sec
            model$cflxca     [tt ] = - myinst$FMEAN.CARBON.AC.PY 
            model$cflxst     [tt ] = + myinst$FMEAN.CARBON.ST.PY
            model$gpp        [tt ] =   myinst$FMEAN.GPP.PY
            model$reco       [tt ] = ( myinst$FMEAN.PLRESP.PY
                                     + myinst$FMEAN.RH.PY           )
            model$nep        [tt ] =  model$gpp[tt] - model$reco[tt]
            model$nee        [tt ] = ( myinst$FMEAN.CARBON.ST.PY 
                                     - myinst$FMEAN.CARBON.AC.PY )
            model$ustar      [tt ] =   myinst$FMEAN.USTAR.PY
            model$rlongup    [tt ] =   myinst$FMEAN.RLONGUP.PY
            model$albedo     [tt ] =   myinst$FMEAN.ALBEDO.PY
            model$rnet       [tt ] =   myinst$FMEAN.RNET.PY
            model$parup      [tt ] =   myinst$FMEAN.PARUP.PY          * Watts.2.Ein * 1.e6
            model$rshortup   [tt ] =   myinst$FMEAN.RSHORTUP.PY
            model$soil.tmp   [tt,] =   myinst$FMEAN.SOIL.TEMP.PY      - t00
            model$soil.water [tt,] =   myinst$FMEAN.SOIL.WATER.PY
            model$soil.matpot[tt,] =   myinst$FMEAN.SOIL.MSTPOT.PY

            if (tt == ntimes){
               eddy.complete = TRUE
               eddy.tresume  = ntimes+1
               cat0("   - Save full ED-2.2 data to ",basename(ed22.rdata),".")
               save( list = c("model","eddy.complete","eddy.tresume")
                   , file = ed22.rdata)
            }else if (last.cyear != cyear){
               eddy.complete = FALSE
               eddy.tresume  = tt+1
               cat("   - Save partial ED-2.2 data to ",basename(ed22.rdata),".")
               save( list = c("model","eddy.complete","eddy.tresume")
                   , file = ed22.rdata)
               cat0("Quit.")
               q("no")
            }#end if

            last.cyear = cyear
         }#end for
         #---------------------------------------------------------------------------------#
      }#end if
      #------------------------------------------------------------------------------------#


      #----- Create the directories in case they don't exist. -----------------------------#
      if (! file.exists(outmain)) dir.create(outmain)
      if (! file.exists(outpref)) dir.create(outpref)
      outboxmain = file.path(outpref,"boxplot")
      outpdfmain = file.path(outpref,"pdfplot")
      outqqpmain = file.path(outpref,"qqplot" )
      outlight   = file.path(outpref,"light"  )
      if (! file.exists(outboxmain)) dir.create(outboxmain)
      if (! file.exists(outpdfmain)) dir.create(outpdfmain)
      if (! file.exists(outqqpmain)) dir.create(outqqpmain)
      if (! file.exists(outlight  )) dir.create(outlight  )
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #      Find out how many diel and seasonal blocks to create.                         #
      #------------------------------------------------------------------------------------#
      mydiel    = sort(unique(obser$hr.idx))
      ndiel     = length(mydiel)
      hr.end    = seq(from=hourblock.len,to=day.hr,by=hourblock.len) %% day.hr
      hr.beg    = (hr.end - hourblock.len + 1) %% day.hr
      diel.list = paste0(sprintf("%2.2i",hr.beg),"-",sprintf("%2.2i",hr.end))
      dl.name   = diel.list[mydiel]
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #       Compare the light response curve.                                            #
      #------------------------------------------------------------------------------------#
      cat0("   - Compare the light response curve.")
      for (s in sequence(nseasons)){
         #---------------------------------------------------------------------------------#
         #     Select the data points for which we have observations, then fit the light   #
         # response curve to both observed and modelled GPP.                               #
         #---------------------------------------------------------------------------------#
         s.sel        = obser$yr.idx == s | s == nseasons
         sel          = s.sel & is.finite(obser$gpp) & is.finite(obser$par) & obser$highsun
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #      Make sure there are at least a few points.                                 #
         #---------------------------------------------------------------------------------#
         if (sum(sel) > 10){


            #------------------------------------------------------------------------------#
            #      Observations.                                                           #
            #------------------------------------------------------------------------------#
            #----- Select and sort the data. ----------------------------------------------#
            obser.data  = data.frame( par = obser$par[sel]
                                    , gpp = obser$gpp[sel]  * kgCyr.2.umols)
            o = order(obser.data$par)
            obser.data  = obser.data[o,]
            #----- Give more weight to light-limited side. --------------------------------#
            obser.wgts  = 1.0 / (( obser.data$par * 0.002 + 3.0) )
            #----- Fit the observations. --------------------------------------------------#
            obser.fit   = try( expr   = nls( formula = gpp ~ a1 + a2 * par / (a3 + par)
                                           , data    = obser.data
                                           , weights = obser.wgts
                                           , start   = list(a1 = 1., a2 = -40., a3 = 500.)
                                           )#end nls
                             , silent = TRUE
                             )#end try
            #----- Find the light-response curve. -----------------------------------------#
            if ("try-error" %in% is(obser.fit)){
               obser.pred  = rep(NA,times=length(obser.data$par))
            }else{
               obser.pred  = predict(obser.fit,data=obser.data)
            }#end if
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #      Model results.                                                          #
            #------------------------------------------------------------------------------#
            #----- Select and sort the data. ----------------------------------------------#
            model.data  = data.frame( par = model$par[sel]
                                    , gpp = model$gpp[sel]  * kgCyr.2.umols)
            o = order(model.data$par)
            model.data  = model.data[o,]
            #----- Give more weight to light-limited side. --------------------------------#
            model.wgts  = 1.0 / (( model.data$par * 0.002 + 3.0) )
            #----- Fit the modelled results. ----------------------------------------------#
            model.fit   = try( expr   = nls( formula = gpp ~ a1 + a2 * par / (a3 + par)
                                           , data    = model.data
                                           , weights = model.wgts
                                           , start   = list(a1 = 1., a2 = -40., a3 = 500.)
                                           )#end nls
                             , silent = TRUE
                             )#end try
            #----- Find the light-response curve. -----------------------------------------#
            if ("try-error" %in% is(model.fit)){
               model.pred  = rep(NA,times=length(model.data$par))
            }else{
               model.pred  = predict(model.fit,data=model.data)
            }#end if
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #     Find the X and Y limits.                                                 #
            #------------------------------------------------------------------------------#
            xlimit  = pretty.xylim(u=c(obser.data$par,model.data$par))
            xat     = pretty(xlimit)
            xlabels = sprintf("%g",xat)
            ylimit  = pretty.xylim(u=c(obser.data$gpp,model.data$gpp))
            yat     = pretty(ylimit)
            ylabels = sprintf("%g",yat)
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Set up the title and axes labels.                                        #
            #------------------------------------------------------------------------------#
            letitre = paste0(lieu,"\n Polygon-level light response curve - ",season.full[s])
            lex     = desc.unit( desc = "Photosynthetically Active Radiation"
                               , unit = untab$umolom2os
                               )#end desc.unit
            ley     = desc.unit( desc = "Gross Primary Productivity"
                               , unit = untab$umolcom2os
                               )#end desc.unit
            #---------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #     Plot the predicted curves and scatter plot.                              #
            #------------------------------------------------------------------------------#
            for (o in sequence(nout)){
               #----- Make the file name. -------------------------------------------------#
               fichier = paste0("gpp_light","-",season.label[s],"-",suffix,".",outform[o])
               fichier = file.path(outlight,fichier)
               dummy   = open.plot( fichier = fichier
                                  , outform = outform[o]
                                  , size    = size
                                  , ptsz    = ptsz
                                  , depth   = depth
                                  )#end open.plot
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     Open an empty plotting area.                                          #
               #---------------------------------------------------------------------------#
               par(par.user)
               par(mar=c(4.1,4.6,2.1,1.1))
               plot.new()
               plot.window(xlim=xlimit,ylim=ylimit)
               axis(side=1,las=1,at=xat,labels=xlabels)
               axis(side=2,las=1,at=yat,labels=ylabels)
               title(main=letitre,cex.main=cex.main)
               title(xlab=lex,ylab=ley)
               #----- Add the observations. -----------------------------------------------#
               points(x=obser.data$par,y=obser.data$gpp,pch=16,col=alpha("#404040",0.40))
               points(x=model.data$par,y=model.data$gpp,pch=16,col=alpha("#0E6E81",0.40))
               #----- Add the density functions. ------------------------------------------#
               lines(x=obser.data$par,y=obser.pred,lwd=3.0,col="#404040")
               lines(x=model.data$par,y=model.pred,lwd=3.0,col="#0E6E81")
               #----- Add the legend. -----------------------------------------------------#
               legend( x      = "topleft"
                     , bg     = "transparent"
                     , inset  = 0.01
                     , legend = c("Tower","Fit - Tower","ED-2.2","Fit - ED-2.2")
                     , pch    = c(16,NA,16,NA)
                     , lwd    = 3
                     , lty    = c("blank","solid","blank","solid")
                     , col    = c(alpha("#404040",c(0.40,1.0)),alpha("#0E6E81",c(0.40,1.00)))
                     , cex    = 1.0
                     , bty    = "n"
                     )#end legend
               box()
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     Close the plotting device.                                            #
               #---------------------------------------------------------------------------#
               dummy = close.plot(outform=outform[o])
               #---------------------------------------------------------------------------#
            }#end for
            #------------------------------------------------------------------------------#
         }#end if
         #---------------------------------------------------------------------------------#
      }#end for
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Loop over all variables, and hours blocks, and all seasons to determine the    #
      # distribution of the data.                                                          #
      #------------------------------------------------------------------------------------#
      cat0("   - Compare season and diurnal distributions.")
      dist.comp = list()
      for (cc in sequence(ncompvar)){
         #----- Handy aliases. ------------------------------------------------------------#
         this.vnam  = compvar$vnam      [cc]
         this.desc  = compvar$desc      [cc]
         this.unit  = compvar$unit      [cc]
         col.obser  = compvar$col.obser [cc]
         col.model  = compvar$col.model [cc]
         leg.corner = compvar$leg.corner[cc]
         sunvar     = compvar$sunvar    [cc]
         col.obser  = alpha(col.obser,c(0.40,1.00))
         col.model  = alpha(col.model,c(0.40,1.00))
         this.meas  = paste0("measured.",this.vnam)
         #---------------------------------------------------------------------------------#

         #---------------------------------------------------------------------------------#
         #     Grab the variables that we are going to use now.                            #
         #---------------------------------------------------------------------------------#
         cat0("     * ",this.desc,".")
         this.obser    = obser[[this.vnam]]
         this.measured = obser[[this.meas]]
         this.model    = model[[this.vnam]]
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Initialise the list to keep comparisons between model and observations.     #
         #---------------------------------------------------------------------------------#
         comp = list()
         mat              = matrix( data     = NA
                                  , nrow     = ndiel+4
                                  , ncol     = nseasons
                                  , dimnames = list( c(dl.name,"night","rise.set","day"
                                                              ,"all.hrs")
                                                   , season.label
                                                   )#end list
                                  )#end matrix
         arr              = array ( data     = NA
                                  , dim      = c(ndiel+4,nseasons,4)
                                  , dimnames = list( c(dl.name,"night","rise.set","day"
                                                              ,"all.hrs")
                                                   , season.label
                                                   , c("mean","variance"
                                                      ,"skewness","kurtosis")
                                                   )
                                  )#end array
         comp$n           = mat
         comp$residuals   = rep(NA,times=length(this.obser))
         comp$bias        = mat
         comp$rmse        = mat
         comp$r.squared   = mat
         comp$fvue        = mat
         comp$sw.stat     = mat
         comp$sw.p.value  = mat
         comp$ks.stat     = mat
         comp$ks.p.value  = mat
         comp$lsq.lnlike  = mat
         comp$sn.lnlike   = mat
         comp$norm.lnlike = mat
         comp$obs.moment  = arr
         comp$mod.moment  = arr
         comp$res.moment  = arr
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #      Check whether the variable directories exist.  If not, create them.        #
         #---------------------------------------------------------------------------------#
         outboxvar = file.path(outboxmain,this.vnam)
         outpdfvar = file.path(outpdfmain,this.vnam)
         outqqpvar = file.path(outqqpmain,this.vnam)
         if (! file.exists(outboxvar)) dir.create(outboxvar)
         if (! file.exists(outpdfvar)) dir.create(outpdfvar)
         if (! file.exists(outqqpvar)) dir.create(outqqpvar)
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Season block.                                                               #
         #---------------------------------------------------------------------------------#
         for (s in sequence(nseasons)){
            #------------------------------------------------------------------------------#
            #     Select the season.                                                       #
            #------------------------------------------------------------------------------#
            cat0("       # Season: ",season.full[s],".")
            sel.season = obser$yr.idx == s | s == nseasons
            #------------------------------------------------------------------------------#


            bp.list = list()
            #----- These lists will contain the data for the box plot. --------------------#
            ylimit.bp = NULL
            #------------------------------------------------------------------------------#
            #     Diel block.                                                              #
            #------------------------------------------------------------------------------#
            for (d in sequence(ndiel+4)){
               #---------------------------------------------------------------------------#
               #     Select the period of the day to plot.                                 #
               #---------------------------------------------------------------------------#
               if (d <= ndiel){
                  cat0("         ~ Hour: ",dl.name[d]," UTC.")
                  sel.diel   = obser$hr.idx == mydiel   [d]
                  diel.label = paste0("hr_"    ,dl.name[d])
                  diel.desc  = paste0("Hours: ",dl.name[d]," UTC")
               }else if (d == ndiel+1){
                  cat0("         ~ Night time.")
                  sel.diel   = obser$nighttime
                  diel.label = "night"
                  diel.desc  = "Nighttime"
               }else if (d == ndiel+2){
                  cat0("         ~ Sunrise/sunset time.")
                  sel.diel   = obser$riseset
                  diel.label = "riseset"
                  diel.desc  = "Sunrise/sunset time"
               }else if (d == ndiel+3){
                  cat0("         ~ Day time.")
                  sel.diel   = obser$highsun
                  diel.label = "day"
                  diel.desc  = "Daytime"
               }else if (d == ndiel+4){
                  cat0("         ~ 24 hours.")
                  sel.diel   = rep(TRUE,times=length(obser$hr.idx))
                  diel.label = "allhrs"
                  diel.desc  = "All hours"
               }#end if
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     Correct sel in case this is a "Sun" variable.  We only check Sun      #
               # variables during daytime.                                                 #
               #---------------------------------------------------------------------------#
               if (sunvar){
                  sel.sun = obser$daytime
               }else{
                  sel.sun = rep(TRUE,times=length(obser$daytime))
               }#end if
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     Number of observations that we use to build the comparison metrics.   #
               #---------------------------------------------------------------------------#
               sel   = sel.season & sel.diel & sel.sun & this.measured
               n.sel = sel
               #---------------------------------------------------------------------------#


               #---------------------------------------------------------------------------#
               #     Append the data to the box plot lists, except if this is day time     #
               # or night time plot.                                                       #
               #---------------------------------------------------------------------------#
               if (d <= ndiel){
                  o.bp.name = paste0("Tower ",dl.name[d])
                  m.bp.name = paste0("ED22 " ,dl.name[d])
                  bp.list[[o.bp.name]] = this.obser[sel]
                  bp.list[[m.bp.name]] = this.model[sel]
               }#end if
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     Find the statistics of the observed quantities.                       #
               #---------------------------------------------------------------------------#
               if (any(sel)){


                  #----- Find and plot the distribution function for this hour. -----------#
                  sd.obser = sd(this.obser[sel],na.rm=TRUE)
                  if (sd.obser %>% 1.0e-6){
                     #----- Find the residuals. -------------------------------------------#
                     this.resid          = this.obser - this.model
                     comp$residuals[sel] = this.resid[sel]
                     #---------------------------------------------------------------------#


                     #---------------------------------------------------------------------#
                     #      Find multiple statistics that may be used for finding the      #
                     # support function.                                                   #
                     #---------------------------------------------------------------------#
                     o.location = sn.location(this.obser[sel],na.rm=TRUE)
                     o.scale    = sn.scale   (this.obser[sel],na.rm=TRUE)
                     o.shape    = sn.shape   (this.obser[sel],na.rm=TRUE)
                     o.mean     = mean       (this.obser[sel],na.rm=TRUE)
                     o.vari     = var        (this.obser[sel],na.rm=TRUE)
                     o.sdev     = sd         (this.obser[sel],na.rm=TRUE)
                     o.skew     = skew       (this.obser[sel],na.rm=TRUE)
                     o.kurt     = kurt       (this.obser[sel],na.rm=TRUE)
                     m.location = sn.location(this.model[sel],na.rm=TRUE)
                     m.scale    = sn.scale   (this.model[sel],na.rm=TRUE)
                     m.shape    = sn.shape   (this.model[sel],na.rm=TRUE)
                     m.mean     = mean       (this.model[sel],na.rm=TRUE)
                     m.vari     = var        (this.model[sel],na.rm=TRUE)
                     m.sdev     = sd         (this.model[sel],na.rm=TRUE)
                     m.skew     = skew       (this.model[sel],na.rm=TRUE)
                     m.kurt     = kurt       (this.model[sel],na.rm=TRUE)
                     r.location = sn.location(this.resid[sel],na.rm=TRUE)
                     r.scale    = sn.scale   (this.resid[sel],na.rm=TRUE)
                     r.shape    = sn.shape   (this.resid[sel],na.rm=TRUE)
                     r.mean     = mean       (this.resid[sel],na.rm=TRUE)
                     r.vari     = var        (this.resid[sel],na.rm=TRUE)
                     r.sdev     = sd         (this.resid[sel],na.rm=TRUE)
                     r.skew     = skew       (this.resid[sel],na.rm=TRUE)
                     r.kurt     = kurt       (this.resid[sel],na.rm=TRUE)
                     #---------------------------------------------------------------------#



                     #---------------------------------------------------------------------#
                     #     Standardise the quantiles for the density function.  Always     #
                     # choose a power of two for the total length, because it reduces the  #
                     # amount of interpolation if we go with the empirical density         #
                     # function.                                                           #
                     #---------------------------------------------------------------------#
                     qlimit  = pretty.xylim(u=c(this.obser[sel],this.model[sel]))
                     quant   = seq(from=qlimit[1],to=qlimit[2],length.out=n.quant)
                     qat     = pretty(qlimit)
                     qlabels = sprintf("%g",qat)
                     #---------------------------------------------------------------------#



                     #---------------------------------------------------------------------#
                     #      Find the distribution curves for model and observations.       #
                     #---------------------------------------------------------------------#
                     if ( use.distrib %in% "sn"){
                        dfunc.obser = dsn(x=quant,dp=c(o.location,o.scale,o.shape))
                        dfunc.model = dsn(x=quant,dp=c(m.location,m.scale,m.shape))
                     }else if (use.distrib %in% "norm"){
                        dfunc.obser = dnorm(x=quant,mean=o.mean,sd=o.sdev)
                        dfunc.model = dnorm(x=quant,mean=m.mean,sd=m.sdev)
                     }else if (use.distrib %in% "edf"){
                        dfunc.obser = density( x    = this.obser[sel], n  = n.quant
                                             , from = qlimit[1]      , to = qlimit[2])$y
                        dfunc.model = density( x    = this.model[sel], n  = n.quant
                                             , from = qlimit[1]      , to = qlimit[2])$y
                     }#end if
                     #---------------------------------------------------------------------#



                     #---------------------------------------------------------------------#
                     #      Run a Kolgomorov-Smirnov test comparing the two distributions. #
                     #---------------------------------------------------------------------#
                     this.ks               = ks.test(x=this.obser[sel],y=this.model[sel])
                     comp$n          [d,s] = sum(sel)
                     comp$ks.stat    [d,s] = this.ks$statistic
                     comp$ks.p.value [d,s] = this.ks$p.value
                     comp$norm.lnlike[d,s] = sum( x     = dnorm( x         = this.obser[sel]
                                                               , mean      = m.mean
                                                               , sd        = m.sdev
                                                               , log       = TRUE
                                                               )#end dnorm
                                                , na.rm = TRUE
                                                )#end sum
                     comp$sn.lnlike  [d,s] = sum( x     =  dsn  ( x        = this.obser[sel]
                                                                , location = m.location
                                                                , scale    = m.scale
                                                                , shape    = m.shape
                                                                , log      = TRUE
                                                                )#end dnorm
                                                , na.rm = TRUE
                                                )#end sum
                     #---------------------------------------------------------------------#



                     #---------------------------------------------------------------------#
                     #      Find the mean bias, the root mean square error, the coeffi-    #
                     # cient of determination, and the fraction of variance unexplained    #
                     # for this simulation.                                                #
                     #---------------------------------------------------------------------#
                     goodness = test.goodness ( x.mod        = this.model[sel]
                                              , x.obs        = this.obser[sel]
                                              , n.parameters = NULL
                                              )#end test.goodness
                     comp$bias      [d,s ] = goodness$bias
                     comp$rmse      [d,s ] = goodness$rmse
                     comp$lsq.lnlike[d,s ] = goodness$lsq.lnlike
                     comp$r.squared [d,s ] = goodness$r.squared
                     comp$fvue      [d,s ] = goodness$fvue
                     comp$sw.stat   [d,s ] = goodness$sw.statistic
                     comp$sw.pvalue [d,s ] = goodness$sw.pvalue
                     #---------------------------------------------------------------------#


                     #---------------------------------------------------------------------#
                     #      Find the first four moments of the distribution for            #
                     # observations, model, and residuals.                                 #
                     #---------------------------------------------------------------------#
                     comp$obs.moment[d,s,] = c(o.mean,o.vari,o.skew,o.kurt)
                     comp$mod.moment[d,s,] = c(m.mean,m.vari,m.skew,m.kurt)
                     comp$res.moment[d,s,] = c(r.mean,r.vari,r.skew,r.kurt)
                     #---------------------------------------------------------------------#



                     #---------------------------------------------------------------------#
                     #      Find the range of observations/modelled variables.             #
                     #---------------------------------------------------------------------#
                     xbreaks    = pretty(qlimit,n=20)
                     freq.obser = hist(this.obser[sel],breaks=xbreaks,plot=FALSE)$density
                     freq.model = hist(this.model[sel],breaks=xbreaks,plot=FALSE)$density
                     yrange     = c(dfunc.obser,dfunc.model,freq.obser,freq.model)
                     ylimit     = pretty.xylim(u=yrange,fracexp=scalleg,is.log=FALSE)
                     yat        = pretty(ylimit)
                     ylabels    = sprintf("%g",yat)
                     #---------------------------------------------------------------------#



                     #---------------------------------------------------------------------#
                     #     Create edges for the histogram rectangles (we use rect instead  #
                     # of hist for histograms because we cannot control the line width in  #
                     # hist).                                                              #
                     #---------------------------------------------------------------------#
                     xleft   = xbreaks[-length(xbreaks)]
                     xright  = xbreaks[-1]
                     ybottom = 0. * xleft
                     #---------------------------------------------------------------------#




                     #---------------------------------------------------------------------#
                     #     Set up the title and axes labels.                               #
                     #---------------------------------------------------------------------#
                     letitre = paste(lieu,"\n",this.desc," - ",season.full[s]
                                    ," - ",diel.desc,sep="")
                     lex     = desc.unit(desc=this.desc,unit=untab[[this.unit]])
                     ley     = desc.unit(desc="Density function",unit=untab$empty)
                     #---------------------------------------------------------------------#



                     #---------------------------------------------------------------------#
                     #     Plot the histograms and density curves.                         #
                     #---------------------------------------------------------------------#
                     for (o in sequence(nout)){
                        #----- Make the file name. ----------------------------------------#
                        fichier = paste0("histcomp_",this.vnam,"-",season.label[s],"-"
                                        ,diel.label,"-",suffix,".",outform[o])
                        fichier = file.path(outpdfvar,fichier)
                        dummy   = open.plot( fichier = fichier
                                           , outform = outform[o]
                                           , size    = size
                                           , ptsz    = ptsz
                                           , depth   = depth
                                           )#end open.plot
                        #------------------------------------------------------------------#



                        #------------------------------------------------------------------#
                        #     Open an empty plotting area.                                 #
                        #------------------------------------------------------------------#
                        par(par.user)
                        par(mar=c(4.1,4.6,2.1,1.1))
                        plot.new()
                        plot.window(xlim=qlimit,ylim=ylimit)
                        axis(side=1,las=1,at=qat,labels=qlabels)
                        axis(side=2,las=1,at=yat,labels=ylabels)
                        title(main=letitre,cex.main=cex.main)
                        title(xlab=xlab,ylab=ylab)
                        #------------------------------------------------------------------#



                        #------------------------------------------------------------------#
                        #      Plot the histograms using "rect.                            #
                        #------------------------------------------------------------------#
                        #----- Observations. ----------------------------------------------#
                        rect(xleft,ybottom,xright,freq.obser,density=12,angle=-45
                            ,col=col.obser[1],border=col.obser[1],lwd=1.6)
                        #----- Model. -----------------------------------------------------#
                        rect(xleft,ybottom,xright,freq.model,density=12,angle=+45
                            ,col=col.model[1],border=col.model[1],lwd=1.6)
                        #----- Add the density functions. ---------------------------------#
                        lines(x=quant,y=dfunc.obser,lwd=3.0,col=col.obser[2])
                        lines(x=quant,y=dfunc.model,lwd=3.0,col=col.model[2])
                        #----- Add the legend. --------------------------------------------#
                        legend( x       = "topleft"
                              , inset   = 0.01
                              , legend  = c("Tower","ED-2.2")
                              , fill    = c(col.obser[2],col.model[2])
                              , border  = c(col.obser[2],col.model[2])
                              , angle   = c(-45,45)
                              , density = 30
                              , lwd     = 2.0
                              , col     = c(col.obser[2],col.model[2])
                              , bg      = "transparent"
                              , bty     = "n"
                              , cex     = 1.0
                              )#end legend
                        box()
                        #------------------------------------------------------------------#



                        #------------------------------------------------------------------#
                        #     Close the plotting device.                                   #
                        #------------------------------------------------------------------#
                        dummy = close.plot(outform=outform[o])
                        #------------------------------------------------------------------#
                     }#end for (o in sequence(nout))
                     #---------------------------------------------------------------------#




                     #---------------------------------------------------------------------#
                     #     Organise the quantiles for plotting.                            #
                     #---------------------------------------------------------------------#
                     qq       = qqplot(x=this.obser[sel],y=this.model[sel],plot.it=FALSE)
                     xylimit  = pretty.xylim(c(qq$x,qq$y))
                     xyat     = pretty(xylimit)
                     xylabels = sprintf("%g",xyat)
                     #---------------------------------------------------------------------#



                     #---------------------------------------------------------------------#
                     #     Set up the title and axes labels.                               #
                     #---------------------------------------------------------------------#
                     letitre = paste0(lieu,"\n","QQ Plot for ",this.desc
                                     ," - ",season.full[s]," - ",diel.desc)
                     lex     = desc.unit(desc="Tower" ,unit=untab[[this.unit]])
                     ley     = desc.unit(desc="ED-2.2",unit=untab[[this.unit]])
                     #---------------------------------------------------------------------#



                     #---------------------------------------------------------------------#
                     #     Plot the QQ-Plot of the distributions.                          #
                     #---------------------------------------------------------------------#
                     for (o in sequence(nout)){
                        #----- Make the file name. ----------------------------------------#
                        fichier = paste0("qqplot_",this.vnam,"-",season.label[s]
                                       ,"-",diel.label,"-",suffix,".",outform[o])
                        fichier = file.path(outqqpvar,fichier)
                        dummy   = open.plot( fichier = fichier
                                           , outform = outform[o]
                                           , size    = size
                                           , ptsz    = ptsz
                                           , depth   = depth
                                           )#end open.plot
                        #------------------------------------------------------------------#



                        #------------------------------------------------------------------#
                        #     Open an empty plotting area.                                 #
                        #------------------------------------------------------------------#
                        par(par.user)
                        par(mar=c(4.1,4.6,2.1,1.1))
                        plot.new()
                        plot.window(xlim=xylimit,ylim=xylimit)
                        axis(side=1,las=1,at=xyat,labels=xylabels)
                        axis(side=2,las=1,at=xyat,labels=xylabels)
                        lines (x=qq$y,y=qq$y,type="l",col=col.obser[2],lwd=3.0)
                        points(x=qq$x,y=qq$y,type="p",pch=16,cex=0.8,col=col.model[1])
                        box()
                        #------------------------------------------------------------------#



                        #------------------------------------------------------------------#
                        #     Close the plotting device.                                   #
                        #------------------------------------------------------------------#
                        dummy = close.plot(outform=outform[o])
                        #------------------------------------------------------------------#
                     }#end for (o in sequenece(nout))
                     #---------------------------------------------------------------------#
                  }#end if (sd(this.obser[sel]) >= 1.e-6)
                  #------------------------------------------------------------------------#
               }#end if (any(sel))
               #---------------------------------------------------------------------------#
            }#end for (d in sequence(ndiel+4))
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Make the box plot comparing observations with model.                     #
            #------------------------------------------------------------------------------#
            if (length(bp.list) > 0){
               #----- Set axis annotation. ------------------------------------------------#
               xlimit  = c(0,2*ndiel)
               xat     = seq(from=1.5,to=2*ndiel-0.5,by=2)
               xgrid   = seq(from=0.5,to=2*ndiel+0.5,by=2)
               xlabels = dl.name
               ylimit  = pretty.xylim(u=bp.list,fracexp=scalleg,is.log=FALSE)
               yat     = pretty(ylimit)
               ylabels = sprintf("%g",yat)
               #---------------------------------------------------------------------------#


               #----- Set colours for background. -----------------------------------------#
               bpcolour = rep(c(col.obser[1],col.model[1]),times=ndiel)
               #---------------------------------------------------------------------------#



               #----- Set up the title and axes labels. -----------------------------------#
               letitre = paste0(lieu,"\n",this.desc," - ",season.full[s])
               lex     = desc.unit( desc = paste0(hourblock.len,"-hour period")
                                  , unit = untab$gmt
                                  )#end desc.unit
               ley     = desc.unit(desc=this.desc,unit=untab[[this.unit]])
               #------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #    Plot box plots.                                                        #
               #---------------------------------------------------------------------------#
               for (o in sequence(nout)){
                  #----- Make the file name. ----------------------------------------------#
                  fichier = paste0("bpcomp_",this.vnam,"-",season.label[s]
                                  ,"-",suffix,".",outform[o])
                  fichier = file.path(outboxvar,fichier)
                  dummy   = open.plot( fichier = fichier
                                     , outform = outform[o]
                                     , size    = size
                                     , ptsz    = ptsz
                                     , depth   = depth
                                     )#end open.plot
                  #------------------------------------------------------------------------#




                  #----- Plot the box plot. -----------------------------------------------#
                  par(par.user)
                  par(mar=c(4.1,4.6,2.1,1.1))
                  plot.new()
                  plot.window(xlim=xlimit,ylim=ylimit)
                  axis(side=1,at=xat,labels=dl.name)
                  axis(side=2,at=yat,labels=ylabels)
                  title(main=letitre,cex.main=cex.main)
                  title(xlab=lex,ylab=ley)
                  abline(v=xgrid,col=grid.colour,lty="solid")
                  boxplot(x=bp.list,col=bpcolour,notch=TRUE,add=TRUE,show.names=FALSE)
                  legend( x      = leg.corner
                        , inset  = 0.01
                        , legend = c("Tower","ED-2.2")
                        , bg     = "transparent"
                        , fill   = c(col.obser[1],col.model[1])
                        )#end legend
                  box()
                  #------------------------------------------------------------------------#



                  #------------------------------------------------------------------------#
                  #     Close the plot device.                                             #
                  #------------------------------------------------------------------------#
                  dummy = close.plot(outform=outform[o])
                  #------------------------------------------------------------------------#
               }#end for outform
            #------------------------------------------------------------------------------#
            }#end if (length(bp.list) > 0)
            #------------------------------------------------------------------------------#
         }#end for
         #---------------------------------------------------------------------------------#


         #----- Save the comparison list for this variable. -------------------------------#
         dist.comp[[this.vnam]] = comp
         #---------------------------------------------------------------------------------#
      }#end for (cc in sequence(ncompvar))
      #------------------------------------------------------------------------------------#
   }#end if (obs.name %in% ls())
   #---------------------------------------------------------------------------------------#


   #----- Write output. -------------------------------------------------------------------#
   dum = write( x    = "Finished"
              , file = file.path(here,place,"eval_load_complete.txt")
              )#end write

   stat.rdata = file.path(path.data,paste0("comp-",place,".RData"))
   cat0(" + Save statistics on model comparison to ",basename(stat.rdata),".") 
   dum        = save(dist.comp,file=stat.rdata)
   #---------------------------------------------------------------------------------------#
}#end for places
#------------------------------------------------------------------------------------------#

#q("no")
