rm(list=ls())

#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#      Here is the user defined variable section.                                          #
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
here           = "thispath"    # Current directory.
there          = "thatpath"    # Directory where analyses/history are 
srcdir         = "/n/moorcroft_data/mlongo/util/Rsc"      # Source  directory.
outroot        = "thisoutroot"
monthbeg       = thismontha
yearbeg        = thisyeara         # First year to consider
yearend        = thisyearz         # Maximum year to consider
myplaces       = c("thispoly")
sasmonth.short = c(2,5,8,11)
sasmonth.long  = 5
nyears.long    = 25
outform        = thisoutform          # Formats for output file.  Supported formats are:
                                #   - "X11" - for printing on screen
                                #   - "eps" - for postscript printing
                                #   - "png" - for PNG printing
                                #   - "pdf" - for PDF printing

byeold         = TRUE           # Remove old files of the given format?

depth          = 96             # PNG resolution, in pixels per inch
paper          = "letter"       # Paper size, to define the plot shape
ptsz           = 14             # Font size.
lwidth         = 2.5            # Line width
plotgrid       = TRUE           # Should I plot the grid in the background? 

legwhere       = "topleft"      # Where should I place the legend?
inset          = 0.01           # inset distance between legend and edge of plot region.
legbg          = "white"        # Legend background colour.
scalleg        = 0.20
cex.main       = 0.8            # Scale coefficient for the title

hourblock.len  = 3              # Length of the time blocks, in hours
#------------------------------------------------------------------------------------------#

reload.data    = TRUE           # Reload data?
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#      NO NEED TO CHANGE ANYTHING BEYOND THIS POINT UNLESS YOU ARE DEVELOPING THE CODE...  #
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Eddy flux comparisons.                                                               #
#------------------------------------------------------------------------------------------#
compvar       = list()
compvar[[ 1]] = list( vnam       = "hflxca"
                    , desc       = "Sensible heat flux"
                    , unit       = "[W/m2]"
                    , col.obser  = c("gray42","gray21")
                    , col.model  = c("orange1","chocolate4")
                    , leg.corner = "topleft"
                    )#end list
compvar[[ 2]] = list( vnam       = "wflxca"
                    , desc       = "Water vapour flux"
                    , unit       = "[kg/m2/day]"
                    , col.obser  = c("gray42","gray21")
                    , col.model  = c("deepskyblue","royalblue4")
                    , leg.corner = "topleft"
                    )#end list
compvar[[ 3]] = list( vnam       = "cflxca"
                    , desc       = "Carbon dioxide flux"
                    , unit       = "[umol/m2/s]"
                    , col.obser  = c("gray42","gray21")
                    , col.model  = c("chartreuse2","darkgreen")
                    , leg.corner = "bottomright"
                    )#end list
compvar[[ 4]] = list( vnam       = "cflxst"
                    , desc       = "Carbon dioxide storage"
                    , unit       = "[umol/m2/s]"
                    , col.obser  = c("gray42","gray21")
                    , col.model  = c("lightgoldenrod3","darkorange1")
                    , leg.corner = "topleft"
                    )#end list
compvar[[ 5]] = list( vnam       = "gpp"
                    , desc       = "Gross primary productivity"
                    , unit       = "[kgC/m2/yr]"
                    , col.obser  = c("gray42","gray21")
                    , col.model  = c("green3","darkgreen")
                    , leg.corner = "topleft"
                    )#end list
compvar[[ 6]] = list( vnam       = "reco"
                    , desc       = "Ecosystem respiration"
                    , unit       = "[kgC/m2/yr]"
                    , col.obser  = c("gray42","gray21")
                    , col.model  = c("yellow3","peru")
                    , leg.corner = "topleft"
                    )#end list
compvar[[ 7]] = list( vnam       = "nep"
                    , desc       = "Net ecosystem productivity"
                    , unit       = "[kgC/m2/yr]"
                    , col.obser  = c("gray42","gray21")
                    , col.model  = c("olivedrab2","darkolivegreen4")
                    , leg.corner = "topleft"
                    )#end list
compvar[[ 8]] = list( vnam       = "nee"
                    , desc       = "Net ecosystem exchange"
                    , unit       = "[umol/m2/s]"
                    , col.obser  = c("gray42","gray21")
                    , col.model  = c("chartreuse","chartreuse4")
                    , leg.corner = "bottomright"
                    )#end list
compvar[[ 9]] = list( vnam       = "ustar"
                    , desc       = "Friction velocity"
                    , unit       = "[m/s]"
                    , col.obser  = c("gray42","gray21")
                    , col.model  = c("mediumpurple1","purple4")
                    , leg.corner = "topleft"
                    )#end list
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Plot-level comparisons.                                                              #
#------------------------------------------------------------------------------------------#
plotvar       = list()
plotvar[[ 1]] = list( vnam.ed    = "recr"
                    , vnam.obs   = "recr"
                    , desc       = "Recruitment rate"
                    , unit       = "[%pop/yr]"
                    , col.obser  = c("gray42","gray21")
                    , col.model  = "chartreuse4"
                    , leg.corner = "topleft"
                    , plog       = TRUE
                    )#end list
plotvar[[ 2]] = list( vnam.ed    = "mort.plot"
                    , vnam.obs   = "mort.plot"
                    , desc       = "Total mortality rate"
                    , unit       = "[%pop/yr]"
                    , col.obser  = c("gray42","gray21")
                    , col.model  = "purple4"
                    , leg.corner = "topleft"
                    , plog       = TRUE
                    )#end list
plotvar[[ 3]] = list( vnam.ed    = "ddmort.plot"
                    , vnam.obs   = "mort.plot"
                    , desc       = "Density dependent mortality rate"
                    , unit       = "[%pop/yr]"
                    , col.obser  = c("gray42","gray21")
                    , col.model  = "mediumpurple1"
                    , leg.corner = "topleft"
                    , plog       = TRUE
                    )#end list
plotvar[[ 4]] = list( vnam.ed    = "dimort.plot"
                    , vnam.obs   = "mort.plot"
                    , desc       = "Density independent mortality rate"
                    , unit       = "[%pop/yr]"
                    , col.obser  = c("gray42","gray21")
                    , col.model  = "mediumpurple1"
                    , leg.corner = "topleft"
                    , plog       = TRUE
                    )#end list
plotvar[[ 5]] = list( vnam.ed    = "growth.plot"
                    , vnam.obs   = "growth.plot"
                    , desc       = "Growth rate"
                    , unit       = "[%DBH/yr]"
                    , col.obser  = c("gray42","gray21")
                    , col.model  = "royalblue4"
                    , leg.corner = "topleft"
                    , plog       = TRUE
                    )#end list
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Plot-level comparisons.                                                              #
#------------------------------------------------------------------------------------------#
sizevar       = list()
sizevar[[ 1]] = list( vnam.ed    = "mort.size"
                    , vnam.obs   = "mort.size"
                    , desc       = "Total mortality rate"
                    , unit       = "[%pop/yr]"
                    , col.obser  = c("gray42","gray21")
                    , col.model  = "purple4"
                    , leg.corner = "topleft"
                    , plog       = TRUE
                    )#end list
sizevar[[ 2]] = list( vnam.ed    = "ddmort.size"
                    , vnam.obs   = "mort.size"
                    , desc       = "Density dependent mortality rate"
                    , unit       = "[%pop/yr]"
                    , col.obser  = c("gray42","gray21")
                    , col.model  = "mediumpurple1"
                    , leg.corner = "topleft"
                    , plog       = TRUE
                    )#end list
sizevar[[ 3]] = list( vnam.ed    = "dimort.size"
                    , vnam.obs   = "mort.size"
                    , desc       = "Density independent mortality rate"
                    , unit       = "[%pop/yr]"
                    , col.obser  = c("gray42","gray21")
                    , col.model  = "mediumpurple1"
                    , leg.corner = "topleft"
                    , plog       = TRUE
                    )#end list
sizevar[[ 4]] = list( vnam.ed    = "growth.size"
                    , vnam.obs   = "growth.size"
                    , desc       = "Growth rate"
                    , unit       = "[%DBH/yr]"
                    , col.obser  = c("gray42","gray21")
                    , col.model  = "royalblue4"
                    , leg.corner = "topleft"
                    , plog       = TRUE
                    )#end list
#------------------------------------------------------------------------------------------#


#----- Load some packages. ----------------------------------------------------------------#
library(hdf5)
library(chron)
library(scatterplot3d)
library(lattice)
library(maps)
library(mapdata)
library(akima)
library(Hmisc)
library(sn)
#------------------------------------------------------------------------------------------#


#----- In case there is some graphic still opened. ----------------------------------------#
graphics.off()
#------------------------------------------------------------------------------------------#


#----- Set how many formats we must output. -----------------------------------------------#
outform = tolower(outform)
nout = length(outform)
#------------------------------------------------------------------------------------------#


#----- Set how many variables we will compare. --------------------------------------------#
ncompvar = length(compvar)
nplotvar = length(plotvar)
nsizevar = length(sizevar)
#------------------------------------------------------------------------------------------#


#----- Avoid unecessary and extremely annoying beeps. -------------------------------------#
options(locatorBell=FALSE)
#------------------------------------------------------------------------------------------#


#----- Load some files with functions. ----------------------------------------------------#
source(paste(srcdir,"atlas.r"           ,sep="/"))
source(paste(srcdir,"charutils.r"       ,sep="/"))
source(paste(srcdir,"census.r"          ,sep="/"))
source(paste(srcdir,"cloudy.r"          ,sep="/"))
source(paste(srcdir,"epolygon.r"       ,sep="/"))
source(paste(srcdir,"error.bar.r"       ,sep="/"))
source(paste(srcdir,"globdims.r"        ,sep="/"))
source(paste(srcdir,"locations.r"       ,sep="/"))
source(paste(srcdir,"muitas.r"          ,sep="/"))
source(paste(srcdir,"plotsize.r"        ,sep="/"))
source(paste(srcdir,"pretty.log.r"      ,sep="/"))
source(paste(srcdir,"pretty.time.r"     ,sep="/"))
source(paste(srcdir,"qapply.r"          ,sep="/"))
source(paste(srcdir,"rconstants.r"      ,sep="/"))
source(paste(srcdir,"skewnorm.stats.r"  ,sep="/"))
source(paste(srcdir,"soilutils.r"       ,sep="/"))
source(paste(srcdir,"sombreado.r"       ,sep="/"))
source(paste(srcdir,"southammap.r"      ,sep="/"))
source(paste(srcdir,"thermlib.r"        ,sep="/"))
source(paste(srcdir,"timeutils.r"       ,sep="/"))
source(paste(srcdir,"zen.r"             ,sep="/"))
#----- These should be called after the others. -------------------------------------------#
source(paste(srcdir,"pft.coms.r"        ,sep="/"))
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
      print(paste(" + Evaluating model for ",lieu,"...",sep=""))

      #----- Get the observed variables. --------------------------------------------------#
      obser        = get(obs.name)
      ntimes       = length(obser$when)
      obser        = alltimes(datin=obser,lon=lon,lat=lat,ed21=TRUE,zeronight=FALSE
                             ,meanval=TRUE,imetavg=1,nmean=12,na.rm=TRUE)
      obser$hr.idx = period.day(obser$when,dtblock=hourblock.len)
      obser$yr.idx = season(obser$when,add.year=FALSE)
      #----- Convert the NEE back to umol/m2/s, normal units for eddy flux. ---------------#
      obser$nee    = obser$nee  * kgCyr.2.umols
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
         print(paste("   - Loading previous session...",sep=""))
         load(ed22.rdata)
         if ((! "eddy.complete" %in% ls()) && "complete" %in% ls()) eddy.complete = complete
         if ((! "eddy.tresume"  %in% ls()) && "tresume"  %in% ls()) eddy.tresume  = tresume
      }else{
         print(paste("   - Starting new session...",sep=""))
         eddy.tresume    = 1
         eddy.complete   = FALSE
      }#end if

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
         }#end if
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Load data for all times.                                                   #
         #---------------------------------------------------------------------------------#
         last.cday   = "00"
         last.cmonth = sprintf("%4.4i",numyears (model$when[eddy.tresume]))
         last.cyear  = sprintf("%4.4i",numyears (model$when[eddy.tresume]))
         print(paste("   - Reading in files...",sep=""))
         for (tt in eddy.tresume:ntimes){
            cyear  = sprintf("%4.4i",numyears (model$when[tt]))
            cmonth = sprintf("%2.2i",nummonths(model$when[tt]))
            cday   = sprintf("%2.2i",numdays  (model$when[tt]))
            chour  = sprintf("%2.2i",hours    (model$when[tt]))
            cminu  = sprintf("%2.2i",minutes  (model$when[tt]))
            myfile = paste(inpref,"-I-",cyear,"-",cmonth,"-",cday,"-",chour,cminu
                           ,"00-g01.h5",sep="")
            if (last.cday != cday) print(paste("     * ",basename(myfile),"...",sep=""))
            last.cday   = cday
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #      Read the data in case the file exists.  In case it doesn't, save what   #
            # has been already read and quit.                                              #
            #------------------------------------------------------------------------------#
            if (! file.exists(myfile)){
               eddy.complete = FALSE
               eddy.tresume  = tt
               print(paste("   - Simulation not ready yet.  Saving partial ED-2.2 data to "
                    ,basename(ed22.rdata),"...",sep=""))
               save( list = c("model","eddy.complete","eddy.tresume")
                   , file = ed22.rdata)
               print("Quitting")
               q("no")
            }else{
               myinst = hdf5load(file=myfile,load=FALSE,verbosity=0,tidy=TRUE)
            }#end if
            #------------------------------------------------------------------------------#
            model$atm.tmp  [tt] = myinst$AVG.ATM.TMP  - t00
            model$atm.shv  [tt] = myinst$AVG.ATM.SHV  * 1000.
            model$atm.prss [tt] = myinst$AVG.ATM.PRSS * 0.01
            model$rain     [tt] = myinst$AVG.PCPG     * hr.sec
            model$atm.co2  [tt] = myinst$AVG.ATM.CO2
            model$atm.vels [tt] = myinst$AVG.VELS
            model$rshort   [tt] = myinst$AVG.RSHORT
            model$rlong    [tt] = myinst$AVG.RLONG
            model$par      [tt] = ( ( myinst$AVG.PAR.BEAM + myinst$AVG.PAR.DIFF ) 
                                  * Watts.2.Ein * 1.e6)
            model$hflxca   [tt] = - myinst$AVG.SENSIBLE.AC
            model$wflxca   [tt] = - myinst$AVG.VAPOR.AC    * day.sec
            model$cflxca   [tt] = - myinst$AVG.CARBON.AC 
            model$cflxst   [tt] = + myinst$AVG.CARBON.ST
            model$gpp      [tt] = myinst$AVG.GPP * umols.2.kgCyr
            model$reco     [tt] = ( ( myinst$AVG.PLANT.RESP + myinst$AVG.HTROPH.RESP )
                                  * umols.2.kgCyr )
            model$nep      [tt] = model$gpp[tt] - model$reco[tt]
            model$nee      [tt] = - model$nep[tt] * kgCyr.2.umols
            model$ustar    [tt] = myinst$AVG.USTAR

            if (tt == ntimes){
               eddy.complete = TRUE
               eddy.tresume  = ntimes+1
               print(paste("   - Saving full ED-2.2 data to ",basename(ed22.rdata)
                          ,"...",sep=""))
               save( list = c("model","eddy.complete","eddy.tresume")
                   , file = ed22.rdata)
            }else if (last.cyear != cyear){
               eddy.complete = FALSE
               eddy.tresume  = tt+1
               print(paste("   - Saving partial ED-2.2 data to ",basename(ed22.rdata)
                          ,"...",sep=""))
               save( list = c("model","eddy.complete","eddy.tresume")
                   , file = ed22.rdata)
               print("Quitting")
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
      outlight   = paste(outpref,"light"  ,sep="/")
      if (! file.exists(outboxmain)) dir.create(outboxmain)
      if (! file.exists(outpdfmain)) dir.create(outpdfmain)
      if (! file.exists(outlight  )) dir.create(outlight)
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #      Find out how many diel and seasonal blocks to create.                         #
      #------------------------------------------------------------------------------------#
      myseasons = sort(unique(obser$yr.idx))
      nseasons  = length(myseasons)
      ss.name   = paste(sprintf("%2.2i",myseasons),season.list[myseasons],sep="-")
      ss.title  = season.full[myseasons] 
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
      print(paste("   - Comparing the light response curve...",sep=""))
      for (s in 1:nseasons){
         #---------------------------------------------------------------------------------#
         #     Select the data points for which we have observations, then fit the light   #
         # response curve to both observed and modelled GPP.                               #
         #---------------------------------------------------------------------------------#
         sel         = ( obser$yr.idx == myseasons[s] & is.finite(obser$gpp)
                       & (is.finite(obser$rshort) & obser$daytime) )

         #---------------------------------------------------------------------------------#
         #      Observations.                                                              #
         #---------------------------------------------------------------------------------#
         #----- Select and sort the data. -------------------------------------------------#
         obser.data  = data.frame( par = obser$par[sel]
                                 , gpp = obser$gpp[sel]  * kgCyr.2.umols)
         o = order(obser.data$par)
         obser.data  = obser.data[o,]
         #----- Give more weight to light-limited side. -----------------------------------#
         obser.wgts  = 1.0 / (( obser.data$par * 0.002 + 3.0) )
         #----- Fit the observations. -----------------------------------------------------#
         obser.fit   = nls( formula = gpp ~ a1 + a2 * par / (a3 + par)
                          , data    = obser.data
                          , weights = obser.wgts
                          , start   = list(a1 = 1., a2 = -40., a3 = 500.) )
         #----- Find the light-response curve. --------------------------------------------#
         obser.pred  = predict(obser.fit)
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #      Model results.                                                             #
         #-----name----------------------------------------------------------------------------#
         #----- Select and sort the data. -------------------------------------------------#
         model.data  = data.frame( par = model$par[sel]
                                 , gpp = model$gpp[sel]  * kgCyr.2.umols)
         o = order(model.data$par)
         model.data  = model.data[o,]
         #----- Give more weight to light-limited side. -----------------------------------#
         model.wgts  = 1.0 / (( model.data$par * 0.002 + 3.0) )
         #----- Fit the modelled results. -------------------------------------------------#
         model.fit   = nls( formula = gpp ~ a1 + a2 * par / (a3 + par)
                          , data    = model.data
                          , weights = model.wgts
                          , start   = list(a1 = 1., a2 = -40., a3 = 500.) )
         #----- Find the light-response curve. --------------------------------------------#
         model.pred  = predict(model.fit)
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #     Find the X and Y limits.                                                    #
         #---------------------------------------------------------------------------------#
         xlimit = range(c(obser.data$par,model.data$par))
         ylimit = range(c(obser.data$gpp,model.data$gpp))
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #     Plot the predicted curves and scatter plot.                                 #
         #---------------------------------------------------------------------------------#
         for (o in 1:nout){
            #----- Make the file name. ----------------------------------------------------#
            fichier = paste(outlight,"/gpp_light","-",ss.name[s],"-",suffix,"."
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
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Set up the title and axes labels.                                        #
            #------------------------------------------------------------------------------#
            letitre = paste(lieu,"\n Polygon-level light response curve - "
                                ,ss.title[s],sep="")
            lex     = paste("Photosynthetically Active Radiation [umol/m2/s]")
            ley     = paste("Gross Primary Productivity [umol/m2/s]")
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Open an empty plotting area.                                             #
            #------------------------------------------------------------------------------#
            plot(x=obser.data$par,y=obser.data$gpp,type="n",main=letitre,xlab=lex,ylab=ley
                ,cex.main=cex.main,ylim=ylimit)
            grid(col="gray76",lty="solid")
            #----- Add the observations. --------------------------------------------------#
            points(x=obser.data$par,y=obser.data$gpp,pch=16,col="gray50")
            points(x=model.data$par,y=model.data$gpp,pch=16,col="chartreuse")
            #----- Add the density functions. ---------------------------------------------#
            lines(x=obser.data$par,y=obser.pred,lwd=3.0,col="gray21")
            lines(x=model.data$par,y=model.pred,lwd=3.0,col="chartreuse4")
            #----- Add the legend. --------------------------------------------------------#
            legend( x      = "topleft"
                  , inset  = 0.01
                  , legend = c("Observation","Fit - Observation","Model","Fit - Model")
                  , pch    = c(16,NA,16,NA)
                  , lwd    = c(NA, 3,NA, 3)
                  , col    = c("gray50","gray21","chartreuse","chartreuse4")
                  , bg     = "white"
                  , cex    = 1.0)
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Close the plotting device.                                               #
            #------------------------------------------------------------------------------#
            if (outform[o] == "x11"){
               locator(n=1)
               dev.off()
            }else{
               dev.off()
            }#end if
            #------------------------------------------------------------------------------#
         }#end for
         #---------------------------------------------------------------------------------#
      }#end for
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Loop over all variables, and hours blocks, and all seasons to determine the    #
      # distribution of the data.                                                          #
      #------------------------------------------------------------------------------------#
      print(paste("   - Comparing season and diurnal distributions...",sep=""))
      for (cc in 1:ncompvar){

         this.comp  = compvar[[cc]]
         this.vnam  = this.comp$vnam
         this.desc  = this.comp$desc
         this.unit  = this.comp$unit
         col.obser  = this.comp$col.obser
         col.model  = this.comp$col.model
         leg.corner = this.comp$leg.corner
         
         print(paste("     * ",this.desc,"...",sep=""))
         this.obser = obser[[this.vnam]]
         this.model = model[[this.vnam]]


         #---------------------------------------------------------------------------------#
         #      Check whether the variable directories exist.  If not, create them.        #
         #---------------------------------------------------------------------------------#
         outboxvar = paste(outboxmain,this.vnam,sep="/")
         outpdfvar = paste(outpdfmain,this.vnam,sep="/")
         if (! file.exists(outboxvar)) dir.create(outboxvar)
         if (! file.exists(outpdfvar)) dir.create(outpdfvar)
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Season block.                                                               #
         #---------------------------------------------------------------------------------#
         for (s in 1:nseasons){
            print(paste("       # Season: ",ss.name[s],"...",sep=""))
            bp.list = list()
            #----- These lists will contain the data for the box plot. --------------------#
            ylimit.bp = NULL
            #------------------------------------------------------------------------------#
            #     Diel block.                                                              #
            #------------------------------------------------------------------------------#
            for (d in 1:(ndiel+3)){
               #---------------------------------------------------------------------------#
               #     Select the period of the day to plot.                                 #
               #---------------------------------------------------------------------------#
               if (d <= ndiel){
                  print(paste("         ~ Hour: ",dl.name[d]," UTC...",sep=""))
                  sel        = ( is.finite(this.obser) & obser$yr.idx == myseasons[s]
                                                       & obser$hr.idx == mydiel   [d] )
                  diel.label = paste("hr",dl.name[d],sep="_")
                  diel.desc  = paste("Hours: ",dl.name[d]," UTC",sep="")
               }else if (d == ndiel+1){
                  print(paste("         ~ Night time...",sep=""))
                  sel        = ( is.finite(this.obser) & obser$yr.idx == myseasons[s]
                               & (is.finite(obser$rshort) & obser$nighttime ) )
                  diel.label = paste("night")
                  diel.desc  = paste("Nighttime",sep="")
               }else if (d == ndiel+2){
                  print(paste("         ~ Sunrise/sunset time...",sep=""))
                  sel        = ( is.finite(this.obser) & obser$yr.idx == myseasons[s]
                               & (is.finite(obser$rshort) & obser$riseset ) )
                  diel.label = paste("riseset")
                  diel.desc  = paste("Sunrise/sunset time",sep="")
               }else if (d == ndiel+3){
                  print(paste("         ~ Day time...",sep=""))
                  sel        = ( is.finite(this.obser) & obser$yr.idx == myseasons[s]
                               & (is.finite(obser$rshort) & obser$highsun ) )
                  diel.label = paste("day")
                  diel.desc  = paste("Daytime",sep="")
               }#end if
               #---------------------------------------------------------------------------#


               #---------------------------------------------------------------------------#
               #     Find the skewed normal distribution of the observed quantities.       #
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
                  if (sd(this.obser[sel]) >= 1.e-6){
                     #---------------------------------------------------------------------#
                     #      Skewed-Gaussian statistics.                                    #
                     #---------------------------------------------------------------------#
                     o.location  = sn.location(this.obser[sel],na.rm=TRUE)
                     o.scale     = sn.scale   (this.obser[sel],na.rm=TRUE)
                     o.shape     = sn.shape   (this.obser[sel],na.rm=TRUE)
                     m.location  = sn.location(this.model[sel],na.rm=TRUE)
                     m.scale     = sn.scale   (this.model[sel],na.rm=TRUE)
                     m.shape     = sn.shape   (this.model[sel],na.rm=TRUE)
                     #---------------------------------------------------------------------#


                     #---------------------------------------------------------------------#
                     #      Find the range of observations/modelled variables.             #
                     #---------------------------------------------------------------------#
                     qlimit      = range(c(this.obser[sel],this.model[sel]))
                     quant       = seq(from=qlimit[1],to=qlimit[2],length.out=1000)
                     dfunc.obser = dsn(x=quant,dp=c(o.location,o.scale,o.shape))
                     dfunc.model = dsn(x=quant,dp=c(m.location,m.scale,m.shape))
                     xbreaks     = pretty(qlimit,n=20)
                     freq.obser  = hist(this.obser[sel],breaks=xbreaks,plot=FALSE)$density
                     freq.model  = hist(this.model[sel],breaks=xbreaks,plot=FALSE)$density
                     ylimit      = range(c(dfunc.obser,dfunc.model,freq.obser,freq.model))
                     if ( any(! is.finite(ylimit)) 
                        || (ylimit[1] == ylimit[2] && ylimit[1] == 0)){
                        ylimit = c(-1,1)
                     }else if (ylimit[1] == ylimit[2] ){
                        ylimit[1] = ylimit[1] * ( 1. - sign(ylimit[1]) * ylnudge)
                        ylimit[2] = ylimit[2] * ( 1. + sign(ylimit[2]) * ylnudge)
                     }else{
                        ylimit[2] = ylimit[2] + scalleg * (ylimit[2] - ylimit[1])
                     }#end if
                     #---------------------------------------------------------------------#


                     #---------------------------------------------------------------------#
                     #     Plot the predicted curves and histograms.                       #
                     #---------------------------------------------------------------------#
                     for (o in 1:nout){
                        #----- Make the file name. ----------------------------------------#
                        fichier = paste(outpdfvar,"/histcomp_",this.vnam,"-",ss.name[s]    
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
                        letitre = paste(lieu,"\n",this.desc," - ",ss.title[s]
                                       ," - ",diel.desc,sep="")
                        lex     = paste(this.desc,this.unit,sep=" ")
                        ley     = "Density function [ ]"
                        #------------------------------------------------------------------#



                        #------------------------------------------------------------------#
                        #     Open an empty plotting area.                                 #
                        #------------------------------------------------------------------#
                        plot(x=quant,y=dfunc.obser,type="n",main=letitre,xlab=lex,ylab=ley
                            ,cex.main=cex.main,ylim=ylimit)
                        grid(col="gray76",lty="solid")
                        #----- Add the histograms. ----------------------------------------#
                        hist(this.obser[sel],breaks=xbreaks,freq=FALSE,col=col.obser[1]
                            ,angle=-45,density=12,add=TRUE)
                        hist(this.model[sel],breaks=xbreaks,freq=FALSE,col=col.model[1]
                            ,angle=45,density=12,add=TRUE)
                        #----- Add the density functions. ---------------------------------#
                        lines(x=quant,y=dfunc.obser,lwd=3.0,col=col.obser[2])
                        lines(x=quant,y=dfunc.model,lwd=3.0,col=col.model[2])
                        #----- Add the legend. --------------------------------------------#
                        legend(x="topleft",inset=0.01,legend=c("Observation","Model")
                              ,fill  =c(col.obser[2],col.model[2])
                              ,border=c(col.obser[2],col.model[2])
                              ,angle=c(-45,45),density=30
                              ,lwd=2.0,col=c(col.obser[2],col.model[2]),bg="white"
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
                     }#end for
                     #---------------------------------------------------------------------#
                  }#end if
                  #------------------------------------------------------------------------#
               }#end if
               #---------------------------------------------------------------------------#

            }#end for
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Make the box plot comparing observations with model.                     #
            #------------------------------------------------------------------------------#
            if (length(bp.list) > 0){
               for (o in 1:nout){
                  #----- Make the file name. ----------------------------------------------#
                  fichier = paste(outboxvar,"/bpcomp_",this.vnam,"-",ss.name[s],"-",suffix
                                           ,".",outform[o],sep="")
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
                  ylimit   = range(bp.list)
                  if ( any(! is.finite(ylimit)) 
                     || (ylimit[1] == ylimit[2] && ylimit[1] == 0)){
                     ylimit = c(-1,1)
                  }else if (ylimit[1] == ylimit[2] ){
                     ylimit[1] = ylimit[1] * ( 1. - sign(ylimit[1]) * ylnudge)
                     ylimit[2] = ylimit[2] * ( 1. + sign(ylimit[2]) * ylnudge)
                  }else{
                     ylimit[2] = ylimit[2] + scalleg * (ylimit[2] - ylimit[1])
                  }#end if
                  bpcolour = rep(c(col.obser[1],col.model[1]),times=ndiel)
                  xat      = seq(from=1.5,to=2*ndiel-0.5,by=2)
                  xgrid    = seq(from=0.5,to=2*ndiel+0.5,by=2)
                  #------------------------------------------------------------------------#



                  #----- Set up the title and axes labels. --------------------------------#
                  letitre = paste(lieu,"\n",this.desc," - ",ss.title[s],sep="")
                  lex     = paste(hourblock.len,"-hour period [UTC]",sep="")
                  ley     = paste(this.desc,this.unit,sep=" ")
                  #------------------------------------------------------------------------#




                  #----- Plot the box plot. -----------------------------------------------#
                  plot(x=xlimit,y=ylimit,type="n",main=letitre,xlab=lex,ylab=ley
                      ,cex.main=cex.main,xaxt="n")
                  axis(side=1,at=xat,labels=dl.name)
                  abline(h=axTicks(side=2),v=xgrid,col="gray66",lty="solid")
                  boxplot(x=bp.list,col=bpcolour,notch=TRUE,add=TRUE,show.names=FALSE)
                  legend(x=leg.corner,inset=0.01,legend=c("Observation","Model"),bg="white"
                        ,fill=c(col.obser[1],col.model[1]))
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
      }#end for
      #------------------------------------------------------------------------------------#
   }#end if (obs.name %in% ls())
   #---------------------------------------------------------------------------------------#

   dum = write(x = "Finished",file=paste(here,place,"eval_load_complete.txt",sep="/"))
}#end for places
#------------------------------------------------------------------------------------------#

#q("no")
