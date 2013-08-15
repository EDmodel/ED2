#==========================================================================================#
#==========================================================================================#
#     Leave these commands at the beginning.  They will refresh the session.               #
#------------------------------------------------------------------------------------------#
rm(list=ls())
graphics.off()
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
slz.min        = -5.0           # The deepest depth that trees access water.
idbh.type      = myidbhtype     # Type of DBH class
                                # 1 -- Every 10 cm until 100cm; > 100cm
                                # 2 -- 0-10; 10-20; 20-35; 35-50; 50-70; > 70 (cm)
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





#------------------------------------------------------------------------------------------#
#     Eddy flux comparisons.                                                               #
#------------------------------------------------------------------------------------------#
compvar       = list()
compvar[[ 1]] = list( vnam       = "hflxca"
                    , desc       = "Sensible heat flux"
                    , unit       = untab$wom2
                    , col.obser  = c(grey.bg,grey.fg)
                    , col.model  = c("orange1","chocolate4")
                    , leg.corner = "topleft"
                    , sunvar     = FALSE
                    )#end list
compvar[[ 2]] = list( vnam       = "wflxca"
                    , desc       = "Water vapour flux"
                    , unit       = untab$kgwom2oday
                    , col.obser  = c(grey.bg,grey.fg)
                    , col.model  = c("deepskyblue","royalblue4")
                    , leg.corner = "topleft"
                    , sunvar     = FALSE
                    )#end list
compvar[[ 3]] = list( vnam       = "cflxca"
                    , desc       = "Carbon dioxide flux"
                    , unit       = untab$umolcom2os
                    , col.obser  = c(grey.bg,grey.fg)
                    , col.model  = c("chartreuse2","darkgreen")
                    , leg.corner = "bottomright"
                    , sunvar     = FALSE
                    )#end list
compvar[[ 4]] = list( vnam       = "cflxst"
                    , desc       = "Carbon dioxide storage"
                    , unit       = untab$umolcom2os
                    , col.obser  = c(grey.bg,grey.fg)
                    , col.model  = c("lightgoldenrod3","darkorange1")
                    , leg.corner = "topleft"
                    , sunvar     = FALSE
                    )#end list
compvar[[ 5]] = list( vnam       = "gpp"
                    , desc       = "Gross primary productivity"
                    , unit       = untab$kgcom2oyr
                    , col.obser  = c(grey.bg,grey.fg)
                    , col.model  = c("green3","darkgreen")
                    , leg.corner = "topleft"
                    , sunvar     = TRUE
                    )#end list
compvar[[ 6]] = list( vnam       = "reco"
                    , desc       = "Ecosystem respiration"
                    , unit       = untab$kgcom2oyr
                    , col.obser  = c(grey.bg,grey.fg)
                    , col.model  = c("yellow3","peru")
                    , leg.corner = "topleft"
                    , sunvar     = FALSE
                    )#end list
compvar[[ 7]] = list( vnam       = "nep"
                    , desc       = "Net ecosystem productivity"
                    , unit       = untab$kgcom2oyr
                    , col.obser  = c(grey.bg,grey.fg)
                    , col.model  = c("olivedrab2","darkolivegreen4")
                    , leg.corner = "topleft"
                    , sunvar     = FALSE
                    )#end list
compvar[[ 8]] = list( vnam       = "nee"
                    , desc       = "Net ecosystem exchange"
                    , unit       = untab$umolcom2os
                    , col.obser  = c(grey.bg,grey.fg)
                    , col.model  = c("chartreuse","chartreuse4")
                    , leg.corner = "bottomright"
                    , sunvar     = FALSE
                    )#end list
compvar[[ 9]] = list( vnam       = "ustar"
                    , desc       = "Friction velocity"
                    , unit       = untab$mos
                    , col.obser  = c(grey.bg,grey.fg)
                    , col.model  = c("mediumpurple1","purple4")
                    , leg.corner = "topleft"
                    , sunvar     = FALSE
                    )#end list
compvar[[10]] = list( vnam       = "rlongup"
                    , desc       = "Outgoing longwave radiation"
                    , unit       = untab$wom2
                    , col.obser  = c(grey.bg,grey.fg)
                    , col.model  = c("gold","orangered")
                    , leg.corner = "topleft"
                    , sunvar     = FALSE
                    )#end list
compvar[[11]] = list( vnam       = "rnet"
                    , desc       = "Net radiation"
                    , unit       = untab$wom2
                    , col.obser  = c(grey.bg,grey.fg)
                    , col.model  = c("gold","orangered")
                    , leg.corner = "topleft"
                    , sunvar     = FALSE
                    )#end list
compvar[[12]] = list( vnam       = "albedo"
                    , desc       = "Albedo"
                    , unit       = untab$empty
                    , col.obser  = c(grey.bg,grey.fg)
                    , col.model  = c("orange1","chocolate4")
                    , leg.corner = "topleft"
                    , sunvar     = TRUE
                    )#end list
compvar[[13]] = list( vnam       = "parup"
                    , desc       = "Outgoing PAR"
                    , unit       = untab$umolom2os
                    , col.obser  = c(grey.bg,grey.fg)
                    , col.model  = c("chartreuse","chartreuse4")
                    , leg.corner = "topleft"
                    , sunvar     = TRUE
                    )#end list
compvar[[14]] = list( vnam       = "rshortup"
                    , desc       = "Outgoing SW"
                    , unit       = untab$wom2
                    , col.obser  = c(grey.bg,grey.fg)
                    , col.model  = c("deepskyblue","royalblue3")
                    , leg.corner = "topleft"
                    , sunvar     = TRUE
                    )#end list
#------------------------------------------------------------------------------------------#


#----- Set how many formats we must output. -----------------------------------------------#
outform = tolower(outform)
nout = length(outform)
#------------------------------------------------------------------------------------------#


#----- Set how many variables we will compare. --------------------------------------------#
ncompvar = length(compvar)
#------------------------------------------------------------------------------------------#


#----- Avoid unecessary and extremely annoying beeps. -------------------------------------#
options(locatorBell=FALSE)
#------------------------------------------------------------------------------------------#



#----- Load observations. -----------------------------------------------------------------#
obser.file = paste(srcdir,"LBA_MIP.nogapfill.RData",sep="/")
load(file=obser.file)
#------------------------------------------------------------------------------------------#



#----- Load census data. ------------------------------------------------------------------#
census.file = paste(srcdir,"LBA_MIP.census_summ.RData",sep="/")
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
   thispoi = locations(where=place,here=there,yearbeg=yearbeg,yearend=yearend
                      ,monthbeg=monthbeg)
   inpref  = thispoi$pathin
   outmain = paste(outroot,place,sep="/")
   outpref = paste(outmain,"eval_ed",sep="/")
   lieu    = thispoi$lieu
   lon     = thispoi$lon
   lat     = thispoi$lat
   iata    = thispoi$iata
   suffix  = thispoi$iata
   yeara   = thispoi$yeara
   yearz   = thispoi$yearz
   meszz   = thispoi$monz


   #----- Find the observations for this particular site. ---------------------------------#
   if (iata == "mao" | iata == "bdf"){
      obs.name = "obs.m34"
   }else if(iata == "stm"){
      obs.name = "obs.s67"
   }else if(iata == "rao"){
      obs.name = "obs.pdg"
   }else if(iata == "jpr"){
      obs.name = "obs.fns"
   }else if(iata == "btr"){
      obs.name = "obs.s77"
   }else{
      obs.name = paste("obs.",iata,sep="")
   }#end if
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     We only run this part of the code if there are observations to compare with the   #
   # model.                                                                                #
   #---------------------------------------------------------------------------------------#
   if (obs.name %in% ls()){


      #----- Print a banner to entretain the user. ----------------------------------------#
      cat(" + Evaluating model for ",lieu,"...","\n")

      #----- Get the observed variables. --------------------------------------------------#
      obser        = get(obs.name)
      ntimes       = length(obser$when)
      obser        = alltimes(datin=obser,lon=lon,lat=lat,ed21=TRUE,zeronight=FALSE
                             ,meanval=TRUE,imetavg=1,nmean=12,na.rm=TRUE)
      obser$hr.idx = period.day(obser$when,dtblock=hourblock.len)
      obser$yr.idx = season(obser$when,add.year=FALSE)
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #      Make the RData file name, then we check whether we must read the files again  #
      # or use the stored RData.                                                           #
      #------------------------------------------------------------------------------------#
      path.data  = paste(here,place,"rdata_hour",sep="/")
      if (! file.exists(path.data)) dir.create(path.data)
      ed22.rdata = paste(path.data,paste(place,"RData",sep="."),sep="/")
      if (reload.data && file.exists(ed22.rdata)){
         #----- Load the modelled dataset. ------------------------------------------------#
         cat("   - Loading previous session...","\n")
         load(ed22.rdata)
         if ((! "eddy.complete" %in% ls()) && "complete" %in% ls()) eddy.complete = complete
         if ((! "eddy.tresume"  %in% ls()) && "tresume"  %in% ls()) eddy.tresume  = tresume
      }else{
         cat("   - Starting new session...","\n")
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
         cat("   - Observed met driver became longer, extend model structure...","\n")
         partial        = model
         npartial       = length(partial$when)
         model          = list()
         model$when     = obser$when
         model$hr.idx   = period.day(model$when,dtblock=hourblock.len)
         model$yr.idx   = season(model$when,add.year=FALSE)
         na.pad         = rep(NA,times=ntimes-npartial)
         model$atm.tmp  = c(partial$atm.tmp ,na.pad)
         model$atm.shv  = c(partial$atm.shv ,na.pad)
         model$atm.prss = c(partial$atm.prss,na.pad)
         model$rain     = c(partial$rain    ,na.pad)
         model$atm.co2  = c(partial$atm.co2 ,na.pad)
         model$atm.vels = c(partial$atm.vels,na.pad)
         model$rshort   = c(partial$rshort  ,na.pad)
         model$rlong    = c(partial$rlong   ,na.pad)
         model$par      = c(partial$par     ,na.pad)
         model$hflxca   = c(partial$hflxca  ,na.pad)
         model$wflxca   = c(partial$wflxca  ,na.pad)
         model$cflxca   = c(partial$cflxca  ,na.pad)
         model$cflxst   = c(partial$cflxst  ,na.pad)
         model$gpp      = c(partial$gpp     ,na.pad)
         model$reco     = c(partial$reco    ,na.pad)
         model$nep      = c(partial$nep     ,na.pad)
         model$nee      = c(partial$nee     ,na.pad)
         model$ustar    = c(partial$ustar   ,na.pad)
         model$rnet     = c(partial$rnet    ,na.pad)
         model$rlongup  = c(partial$rlongup ,na.pad)
         model$albedo   = c(partial$albedo  ,na.pad)
         model$parup    = c(partial$parup   ,na.pad)
         model$rshortup = c(partial$rshortup,na.pad)

         #---------------------------------------------------------------------------------#
         #     Sanity check.                                                               #
         #---------------------------------------------------------------------------------#
         if (length(model$atm.tmp) != length(obser$atm.temp)){
            cat(" Mismatch in lengths!!!","\n")
            cat(" Length (obser):   ",length(obser$atm.temp),"\n")
            cat(" Length (model):   ",length(model$atm.tmp) ,"\n")
            cat(" Length (partial): ",length(partial$when)  ,"\n")
            cat(" Length (na.pad):  ",length(na.pad)        ,"\n")
            cat(" TRESUME:          ",eddy.tresume          ,"\n")
            cat(" NTIMES:           ",ntimes                ,"\n")
            stop("Structures should match, check your code...")
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
         #----- Initialise the model structure. -------------------------------------------#
         if (eddy.tresume == 1){
            model          = list()
            model$when     = obser$when
            model$hr.idx   = period.day(model$when,dtblock=hourblock.len)
            model$yr.idx   = season(model$when,add.year=FALSE)
            empty          = rep(NA,times=ntimes)
            model$atm.tmp  = empty
            model$atm.shv  = empty
            model$atm.prss = empty
            model$rain     = empty
            model$atm.co2  = empty
            model$atm.vels = empty
            model$rshort   = empty
            model$rlong    = empty
            model$par      = empty
            model$hflxca   = empty
            model$wflxca   = empty
            model$cflxca   = empty
            model$cflxst   = empty
            model$gpp      = empty
            model$reco     = empty
            model$nep      = empty
            model$nee      = empty
            model$ustar    = empty
            model$rnet     = empty
            model$rlongup  = empty
            model$albedo   = empty
            model$parup    = empty
            model$rshortup = empty
         }#end if
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Load data for all times.                                                   #
         #---------------------------------------------------------------------------------#
         last.cday   = "00"
         last.cmonth = sprintf("%4.4i",numyears (model$when[eddy.tresume]))
         last.cyear  = sprintf("%4.4i",numyears (model$when[eddy.tresume]))
         cat("   - Reading in files...","\n")
         for (tt in eddy.tresume:ntimes){
            cyear  = sprintf("%4.4i",numyears (model$when[tt]))
            cmonth = sprintf("%2.2i",nummonths(model$when[tt]))
            cday   = sprintf("%2.2i",numdays  (model$when[tt]))
            chour  = sprintf("%2.2i",hours    (model$when[tt]))
            cminu  = sprintf("%2.2i",minutes  (model$when[tt]))
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #      Read the data in case the file exists.  In case it doesn't, save what   #
            # has been already read and quit.                                              #
            #------------------------------------------------------------------------------#
            h5file       = paste(inpref,"-I-",cyear,"-",cmonth,"-",cday,"-",chour,cminu
                                ,"00-g01.h5",sep="")
            h5file.bz2   = paste(h5file,"bz2",sep=".")
            h5file.gz    = paste(h5file,"gz" ,sep=".")
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
               cat("   - Simulation not ready yet.  Saving partial ED-2.2 data to "
                  ,basename(ed22.rdata),"...","\n")
               save( list = c("model","eddy.complete","eddy.tresume")
                   , file = ed22.rdata)
               cat("Quitting","\n")
               q("no")
            }#end if
            if (last.cday != cday) cat("     * ",basename(myinst),"...","\n")
            last.cday   = cday
            #------------------------------------------------------------------------------#
            model$atm.tmp  [tt] =   myinst$FMEAN.ATM.TEMP.PY       - t00
            model$atm.shv  [tt] =   myinst$FMEAN.ATM.SHV.PY        * 1000.
            model$atm.prss [tt] =   myinst$FMEAN.ATM.PRSS.PY       * 0.01
            model$rain     [tt] =   myinst$FMEAN.PCPG.PY           * hr.sec
            model$atm.co2  [tt] =   myinst$FMEAN.ATM.CO2.PY
            model$atm.vels [tt] =   myinst$FMEAN.ATM.VELS.PY
            model$rshort   [tt] =   myinst$FMEAN.ATM.RSHORT.PY
            model$rlong    [tt] =   myinst$FMEAN.ATM.RLONG.PY
            model$par      [tt] =   myinst$FMEAN.ATM.PAR.PY        * Watts.2.Ein * 1.e6
            model$hflxca   [tt] = - myinst$FMEAN.SENSIBLE.AC.PY
            model$wflxca   [tt] = - myinst$FMEAN.VAPOR.AC.PY       * day.sec
            model$cflxca   [tt] = - myinst$FMEAN.CARBON.AC.PY 
            model$cflxst   [tt] = + myinst$FMEAN.CARBON.ST.PY
            model$gpp      [tt] =   myinst$FMEAN.GPP.PY
            model$reco     [tt] = ( myinst$FMEAN.PLRESP.PY
                                  + myinst$FMEAN.RH.PY           )
            model$nep      [tt] =   myinst$FMEAN.NEP.PY
            model$nee      [tt] = ( myinst$FMEAN.CARBON.ST.PY 
                                  - myinst$FMEAN.CARBON.AC.PY )
            model$ustar    [tt] =   myinst$FMEAN.USTAR.PY
            model$rlongup  [tt] =   myinst$FMEAN.RLONGUP.PY
            model$albedo   [tt] =   myinst$FMEAN.ALBEDO.PY
            model$rnet     [tt] =   myinst$FMEAN.RNET.PY
            model$parup    [tt] =   myinst$FMEAN.PARUP.PY          * Watts.2.Ein * 1.e6
            model$rshortup [tt] =   myinst$FMEAN.RSHORTUP.PY

            if (tt == ntimes){
               eddy.complete = TRUE
               eddy.tresume  = ntimes+1
               cat("   - Saving full ED-2.2 data to ",basename(ed22.rdata),"...","\n")
               save( list = c("model","eddy.complete","eddy.tresume")
                   , file = ed22.rdata)
            }else if (last.cyear != cyear){
               eddy.complete = FALSE
               eddy.tresume  = tt+1
               cat("   - Saving partial ED-2.2 data to ",basename(ed22.rdata),"...","\n")
               save( list = c("model","eddy.complete","eddy.tresume")
                   , file = ed22.rdata)
               cat("Quitting","\n")
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
      outboxmain = paste(outpref,"boxplot",sep="/")
      outpdfmain = paste(outpref,"pdfplot",sep="/")
      outqqpmain = paste(outpref,"qqplot" ,sep="/")
      outlight   = paste(outpref,"light"  ,sep="/")
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
      diel.list = paste(sprintf("%2.2i",hr.beg),sprintf("%2.2i",hr.end),sep="-")
      dl.name   = diel.list[mydiel]
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #       Compare the light response curve.                                            #
      #------------------------------------------------------------------------------------#
      cat("   - Comparing the light response curve...","\n")
      for (s in 1:nseasons){
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
            xlimit = range(c(obser.data$par,model.data$par),na.rm=TRUE)
            ylimit = range(c(obser.data$gpp,model.data$gpp),na.rm=TRUE)
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #     Plot the predicted curves and scatter plot.                              #
            #------------------------------------------------------------------------------#
            for (o in 1:nout){
               #----- Make the file name. -------------------------------------------------#
               fichier = paste(outlight,"/gpp_light","-",season.label[s],"-",suffix,"."
                              ,outform[o],sep="")
               if (outform[o] == "x11"){
                  X11(width=size$width,height=size$height,pointsize=ptsz)
               }else if(outform[o] == "png"){
                  png(filename=fichier,width=size$width*depth
                     ,height=size$height*depth,pointsize=ptsz,res=depth)
               }else if(outform[o] == "eps"){
                  postscript(file=fichier,width=size$width,height=size$height
                            ,pointsize=ptsz,paper=size$paper)
               }else if(outform[o] == "pdf"){
                  pdf(file=fichier,onefile=FALSE
                     ,width=size$width,height=size$height,pointsize=ptsz,paper=size$paper)
               }#end if
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     Set up the title and axes labels.                                     #
               #---------------------------------------------------------------------------#
               letitre = paste(lieu,"\n Polygon-level light response curve - "
                                   ,season.full[s],sep="")
               lex     = paste("Photosynthetically Active Radiation [umol/m2/s]")
               ley     = paste("Gross Primary Productivity [umol/m2/s]")
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     Open an empty plotting area.                                          #
               #---------------------------------------------------------------------------#
               par(par.user)
               plot(x=obser.data$par,y=obser.data$gpp,type="n",main=letitre
                   ,xlab=lex,ylab=ley,cex.main=cex.main,ylim=ylimit)
               grid(col=grid.colour,lty="solid")
               #----- Add the observations. -----------------------------------------------#
               points(x=obser.data$par,y=obser.data$gpp,pch=16,col=grey.bg)
               points(x=model.data$par,y=model.data$gpp,pch=16,col="chartreuse")
               #----- Add the density functions. ------------------------------------------#
               lines(x=obser.data$par,y=obser.pred,lwd=3.0,col=grey.fg)
               lines(x=model.data$par,y=model.pred,lwd=3.0,col="chartreuse4")
               #----- Add the legend. -----------------------------------------------------#
               legend( x      = "topleft"
                     , bg     = background
                     , inset  = 0.01
                     , legend = c("Observation","Fit - Observation","Model","Fit - Model")
                     , pch    = c(16,NA,16,NA)
                     , lwd    = c(NA, 3,NA, 3)
                     , col    = c(grey.bg,grey.fg,"chartreuse","chartreuse4")
                     , cex    = 1.0)
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     Close the plotting device.                                            #
               #---------------------------------------------------------------------------#
               if (outform[o] == "x11"){
                  locator(n=1)
                  dev.off()
               }else{
                  dev.off()
               }#end if
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
      cat("   - Comparing season and diurnal distributions...","\n")
      dist.comp = list()
      for (cc in 1:ncompvar){

         this.comp  = compvar[[cc]]
         this.vnam  = this.comp$vnam
         this.meas  = paste("measured",this.vnam,sep=".")
         this.desc  = this.comp$desc
         this.unit  = this.comp$unit
         col.obser  = this.comp$col.obser
         col.model  = this.comp$col.model
         leg.corner = this.comp$leg.corner
         sunvar     = this.comp$sunvar

         #---------------------------------------------------------------------------------#
         #     Grab the variables that we are going to use now.                            #
         #---------------------------------------------------------------------------------#
         cat("     * ",this.desc,"...","\n")
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
         outboxvar = paste(outboxmain,this.vnam,sep="/")
         outpdfvar = paste(outpdfmain,this.vnam,sep="/")
         outqqpvar = paste(outqqpmain,this.vnam,sep="/")
         if (! file.exists(outboxvar)) dir.create(outboxvar)
         if (! file.exists(outpdfvar)) dir.create(outpdfvar)
         if (! file.exists(outqqpvar)) dir.create(outqqpvar)
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Season block.                                                               #
         #---------------------------------------------------------------------------------#
         for (s in 1:nseasons){
            #------------------------------------------------------------------------------#
            #     Select the season.                                                       #
            #------------------------------------------------------------------------------#
            cat("       # Season: ",season.full[s],"...","\n")
            sel.season = obser$yr.idx == s | s == nseasons
            #------------------------------------------------------------------------------#


            bp.list = list()
            #----- These lists will contain the data for the box plot. --------------------#
            ylimit.bp = NULL
            #------------------------------------------------------------------------------#
            #     Diel block.                                                              #
            #------------------------------------------------------------------------------#
            for (d in 1:(ndiel+4)){
               #---------------------------------------------------------------------------#
               #     Select the period of the day to plot.                                 #
               #---------------------------------------------------------------------------#
               if (d <= ndiel){
                  cat("         ~ Hour: ",dl.name[d]," UTC...","\n")
                  sel.diel   = obser$hr.idx == mydiel   [d]
                  diel.label = paste("hr",dl.name[d],sep="_")
                  diel.desc  = paste("Hours: ",dl.name[d]," UTC",sep="")
               }else if (d == ndiel+1){
                  cat("         ~ Night time...","\n")
                  sel.diel   = obser$nighttime
                  diel.label = paste("night")
                  diel.desc  = paste("Nighttime",sep="")
               }else if (d == ndiel+2){
                  cat("         ~ Sunrise/sunset time...","\n")
                  sel.diel   = obser$riseset
                  diel.label = paste("riseset")
                  diel.desc  = paste("Sunrise/sunset time",sep="")
               }else if (d == ndiel+3){
                  cat("         ~ Day time...","\n")
                  sel.diel   = obser$highsun
                  diel.label = paste("day")
                  diel.desc  = paste("Daytime",sep="")
               }else if (d == ndiel+4){
                  cat("         ~ 24 hours...","\n")
                  sel.diel   = rep(TRUE,times=length(obser$hr.idx))
                  diel.label = paste("allhrs")
                  diel.desc  = paste("All hours",sep="")
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
               #     Find the statistics of the observed quantities.                       #
               #---------------------------------------------------------------------------#
               if (any(sel)){
                  #------------------------------------------------------------------------#
                  #     Append the data to the box plot lists, except if this is day time  #
                  # or night time plot.                                                    #
                  #------------------------------------------------------------------------#
                  if (d <= ndiel){
                     o.bp.name = paste("Obs.",dl.name[d],sep=" ")
                     m.bp.name = paste("ED22",dl.name[d],sep=" ")
                     bp.list[[o.bp.name]] = this.obser[sel]
                     bp.list[[m.bp.name]] = this.model[sel]
                  }#end if
                  #------------------------------------------------------------------------#



                  #----- Find and plot the distribution function for this hour. -----------#
                  sd.obser = sd(this.obser[sel],na.rm=TRUE)
                  if (is.finite(sd.obser) && sd.obser > 1.0e-6){
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
                     qlimit     = range(c(this.obser[sel],this.model[sel]))
                     quant      = seq(from=qlimit[1],to=qlimit[2],length.out=n.quant)
                     #---------------------------------------------------------------------#



                     #---------------------------------------------------------------------#
                     #      Find the distribution curves for model and observations.       #
                     #---------------------------------------------------------------------#
                     if ( use.distrib == "sn"){
                        dfunc.obser = dsn(x=quant,dp=c(o.location,o.scale,o.shape))
                        dfunc.model = dsn(x=quant,dp=c(m.location,m.scale,m.shape))
                     }else if (use.distrib == "norm"){
                        dfunc.obser = dnorm(x=quant,mean=o.mean,sd=o.sdev)
                        dfunc.model = dnorm(x=quant,mean=m.mean,sd=m.sdev)
                     }else if (use.distrib == "edf"){
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
                     xbreaks     = pretty(qlimit,n=20)
                     freq.obser  = hist(this.obser[sel],breaks=xbreaks,plot=FALSE)$density
                     freq.model  = hist(this.model[sel],breaks=xbreaks,plot=FALSE)$density
                     yrange      = c(dfunc.obser,dfunc.model,freq.obser,freq.model)
                     ylimit      = pretty.xylim(u=yrange,fracexp=scalleg,is.log=FALSE)
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
                     #     Plot the histograms and density curves.                         #
                     #---------------------------------------------------------------------#
                     for (o in 1:nout){
                        #----- Make the file name. ----------------------------------------#
                        fichier = paste(outpdfvar,"/histcomp_",this.vnam,"-"
                                       ,season.label[s],"-",diel.label,"-",suffix,"."
                                       ,outform[o],sep="")
                        if (outform[o] == "x11"){
                           X11(width=size$width,height=size$height,pointsize=ptsz)
                        }else if(outform[o] == "png"){
                           png(filename=fichier,width=size$width*depth
                              ,height=size$height*depth,pointsize=ptsz,res=depth)
                        }else if(outform[o] == "eps"){
                           postscript(file=fichier,width=size$width,height=size$height
                                     ,pointsize=ptsz,paper=size$paper)
                        }else if(outform[o] == "pdf"){
                           pdf(file=fichier,onefile=FALSE
                              ,width=size$width,height=size$height,pointsize=ptsz
                              ,paper=size$paper)
                        }#end if
                        #------------------------------------------------------------------#



                        #------------------------------------------------------------------#
                        #     Set up the title and axes labels.                            #
                        #------------------------------------------------------------------#
                        letitre = paste(lieu,"\n",this.desc," - ",season.full[s]
                                       ," - ",diel.desc,sep="")
                        lex     = desc.unit(desc=this.desc,unit=this.unit)
                        ley     = "Density function [ ]"
                        #------------------------------------------------------------------#



                        #------------------------------------------------------------------#
                        #     Open an empty plotting area.                                 #
                        #------------------------------------------------------------------#
                        par(par.user)
                        plot(x=quant,y=dfunc.obser,type="n",main=letitre,xlab=lex,ylab=ley
                            ,cex.main=cex.main,ylim=ylimit)
                        grid(col=grid.colour,lty="solid")
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
                        legend(x="topleft",inset=0.01,legend=c("Observation","Model")
                              ,fill  =c(col.obser[2],col.model[2])
                              ,border=c(col.obser[2],col.model[2])
                              ,angle =c(-45,45),density=30
                              ,lwd=2.0,col=c(col.obser[2],col.model[2]),bg=background
                              ,cex=1.0)
                        #------------------------------------------------------------------#



                        #------------------------------------------------------------------#
                        #     Close the plotting device.                                   #
                        #------------------------------------------------------------------#
                        if (outform[o] == "x11"){
                           locator(n=1)
                           dev.off()
                        }else{
                           dev.off()
                        }#end if
                        #------------------------------------------------------------------#
                     }#end for (o in 1:nout)
                     #---------------------------------------------------------------------#




                     #---------------------------------------------------------------------#
                     #     Organise the quantiles for plotting.                            #
                     #---------------------------------------------------------------------#
                     qq     = qqplot(x=this.obser[sel],y=this.model[sel],plot.it=FALSE)
                     xlimit = range(qq$x)
                     ylimit = range(qq$y)
                     #---------------------------------------------------------------------#



                     #---------------------------------------------------------------------#
                     #     Plot the QQ-Plot of the distributions.                          #
                     #---------------------------------------------------------------------#
                     for (o in 1:nout){
                        #----- Make the file name. ----------------------------------------#
                        fichier = paste(outqqpvar,"/qqplot_",this.vnam,"-",season.label[s]
                                       ,"-",diel.label,"-",suffix,".",outform[o],sep="")
                        if (outform[o] == "x11"){
                           X11(width=size$width,height=size$height,pointsize=ptsz)
                        }else if(outform[o] == "png"){
                           png(filename=fichier,width=size$width*depth
                              ,height=size$height*depth,pointsize=ptsz,res=depth)
                        }else if(outform[o] == "eps"){
                           postscript(file=fichier,width=size$width,height=size$height
                                     ,pointsize=ptsz,paper=size$paper)
                        }else if(outform[o] == "pdf"){
                           pdf(file=fichier,onefile=FALSE
                              ,width=size$width,height=size$height,pointsize=ptsz
                              ,paper=size$paper)
                        }#end if
                        #------------------------------------------------------------------#



                        #------------------------------------------------------------------#
                        #     Set up the title and axes labels.                            #
                        #------------------------------------------------------------------#
                        letitre = paste(lieu,"\n","QQ Plot for ",this.desc
                                       ," - ",season.full[s]," - ",diel.desc,sep="")
                        lex     = desc.unit(desc="Observed",unit=this.unit)
                        ley     = desc.unit(desc="Model"   ,unit=this.unit)
                        #------------------------------------------------------------------#



                        #------------------------------------------------------------------#
                        #     Open an empty plotting area.                                 #
                        #------------------------------------------------------------------#
                        par(par.user)
                        plot(x=xlimit,y=ylimit,type="n",main=letitre,xlab=lex,ylab=ley
                            ,cex.main=cex.main)
                        grid(col=grid.colour,lty="solid")
                        lines(x=qq$y,y=qq$y,type="l",col=col.obser[2],lwd=3.0)
                        points(x=qq$x,y=qq$y,type="p",pch=16,cex=0.8,col=col.model[1])
                        #------------------------------------------------------------------#



                        #------------------------------------------------------------------#
                        #     Close the plotting device.                                   #
                        #------------------------------------------------------------------#
                        if (outform[o] == "x11"){
                           locator(n=1)
                           dev.off()
                        }else{
                           dev.off()
                        }#end if
                        #------------------------------------------------------------------#
                     }#end for (o in 1:nout)
                     #---------------------------------------------------------------------#



                  }#end if (sd(this.obser[sel]) >= 1.e-6)
                  #------------------------------------------------------------------------#
               }#end if (any(sel))
               #---------------------------------------------------------------------------#

            }#end for (d in 1:(ndiel+4))
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Make the box plot comparing observations with model.                     #
            #------------------------------------------------------------------------------#
            if (length(bp.list) > 0){
               for (o in 1:nout){
                  #----- Make the file name. ----------------------------------------------#
                  fichier = paste(outboxvar,"/bpcomp_",this.vnam,"-",season.label[s]
                                           ,"-",suffix,".",outform[o],sep="")
                  if (outform[o] == "x11"){
                     X11(width=size$width,height=size$height,pointsize=ptsz)
                  }else if(outform[o] == "png"){
                     png(filename=fichier,width=size$width*depth,height=size$height*depth
                        ,pointsize=ptsz,res=depth)
                  }else if(outform[o] == "eps"){
                     postscript(file=fichier,width=size$width,height=size$height
                               ,pointsize=ptsz,paper=size$paper)
                  }else if(outform[o] == "pdf"){
                     pdf(file=fichier,onefile=FALSE
                        ,width=size$width,height=size$height,pointsize=ptsz
                        ,paper=size$paper)
                  }#end if
                  #------------------------------------------------------------------------#


                  #----- Set up some plot defaults. ---------------------------------------#
                  xlimit   = c(0,2*ndiel)
                  ylimit   = pretty.xylim(u=bp.list,fracexp=scalleg,is.log=FALSE)
                  bpcolour = rep(c(col.obser[1],col.model[1]),times=ndiel)
                  xat      = seq(from=1.5,to=2*ndiel-0.5,by=2)
                  xgrid    = seq(from=0.5,to=2*ndiel+0.5,by=2)
                  #------------------------------------------------------------------------#



                  #----- Set up the title and axes labels. --------------------------------#
                  letitre = paste(lieu,"\n",this.desc," - ",season.full[s],sep="")
                  lex     = desc.unit( desc = paste(hourblock.len,"-hour period",sep="")
                                     , unit = untab$gmt)
                  ley     = desc.unit(desc=this.desc,unit=this.unit)
                  #------------------------------------------------------------------------#




                  #----- Plot the box plot. -----------------------------------------------#
                  par(par.user)
                  plot(x=xlimit,y=ylimit,type="n",main=letitre,xlab=lex,ylab=ley
                      ,cex.main=cex.main,xaxt="n")
                  axis(side=1,at=xat,labels=dl.name)
                  abline(h=axTicks(side=2),v=xgrid,col=grid.colour,lty="solid")
                  boxplot(x=bp.list,col=bpcolour,notch=TRUE,add=TRUE,show.names=FALSE)
                  legend(x=leg.corner,inset=0.01,legend=c("Observation","Model")
                        ,bg=background,fill=c(col.obser[1],col.model[1]))
                  #------------------------------------------------------------------------#

                  #------------------------------------------------------------------------#
                  #     Close the plot device.                                             #
                  #------------------------------------------------------------------------#
                  if (outform[o] == "x11"){
                     locator(n=1)
                     dev.off()
                  }else{
                     dev.off()
                  }#end if
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
      }#end for (c in 1:ncompvar)
      #------------------------------------------------------------------------------------#
   }#end if (obs.name %in% ls())
   #---------------------------------------------------------------------------------------#

   dum = write(x = "Finished",file=paste(here,place,"eval_load_complete.txt",sep="/"))

   stat.rdata = paste(path.data,paste("comp-",place,".RData",sep=""),sep="/")
   cat(" + Saving statistics on model comparison to ",basename(stat.rdata),"...","\n") 
   dum        = save(dist.comp,file=stat.rdata)
}#end for places
#------------------------------------------------------------------------------------------#

#q("no")
