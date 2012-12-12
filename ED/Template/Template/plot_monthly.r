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
here           = "thispath"                          # Current directory.
there          = "thatpath"                          # Directory where analyses/history are 
srcdir         = "/n/moorcroft_data/mlongo/util/Rsc" # Source  directory.
outroot        = "thisoutroot"                       # Directory for figures
#------------------------------------------------------------------------------------------#


#----- Time options. ----------------------------------------------------------------------#
monthbeg       = thismontha   # First month to use
yearbeg        = thisyeara    # First year to consider
yearend        = thisyearz    # Maximum year to consider
reload.data    = TRUE         # Should I reload partially loaded data?
sasmonth.short = c(2,5,8,11)  # Months for SAS plots (short runs)
sasmonth.long  = 5            # Months for SAS plots (long runs)
nyears.long    = 25           # Runs longer than this are considered long runs.
n.density      = 256          # Number of density points
#------------------------------------------------------------------------------------------#



#----- Name of the simulations. -----------------------------------------------------------#
myplaces       = c("thispoly")
#------------------------------------------------------------------------------------------#



#----- Plot options. ----------------------------------------------------------------------#
outform        = thisoutform            # Formats for output file.  Supported formats are:
                                        #   - "X11" - for printing on screen
                                        #   - "eps" - for postscript printing
                                        #   - "png" - for PNG printing
                                        #   - "pdf" - for PDF printing
depth          = 96                     # PNG resolution, in pixels per inch
paper          = "letter"               # Paper size, to define the plot shape
ptsz           = 14                     # Font size.
lwidth         = 2.5                    # Line width
plotgrid       = TRUE                   # Should I plot the grid in the background? 
sasfixlimits   = FALSE                  # Use a fixed scale for size and age-structure
                                        #    plots? (FALSE will set a suitable scale for
                                        #    each plot)
fcgrid         = TRUE                   # Include a grid on the filled contour plots?
ncolsfc        = 80                     # Target number of colours for filled contour.
hovgrid        = TRUE                   # Include a grid on the Hovmoller plots?
legwhere       = "topleft"              # Where should I place the legend?
inset          = 0.01                   # Inset between legend and edge of plot region.
scalleg        = 0.40                   # Expand y limits by this relative amount to fit
                                        #    the legend
cex.main       = 0.8                    # Scale coefficient for the title
theta          = 315.                   # Azimuth for perspective projection
phi            = 30.                    # Vertical angle for perspective projection
ltheta         = -210.                  # Azimuth angle for light
shade          = 0.125                  # Shade intensity
expz           = 0.5                    # Expansion factor for Z axis
cexmin         = 0.5                    # Minimum "head" size of the lollipop
cexmax         = 3.0                    # Maximum "head" size of the lollipop
ylnudge         = 0.05                  # Nudging factor for ylimit
ptype          = "l"                    # Type of plot
ptyped         = "p"                    # Type of plot
ptypeb         = "o"                    # Type of plot
drought.mark   = mydroughtmark          # Put a background to highlight droughts?
drought.yeara  = mydroughtyeara         # First year that has drought
drought.yearz  = mydroughtyearz         # Last year that has drought
months.drought = mymonthsdrought        # Months with drought
ibackground    = mybackground           # Background settings (check load_everything.r)
#------------------------------------------------------------------------------------------#


#------ Miscellaneous settings. -----------------------------------------------------------#
slz.min        = -5.0         # The deepest depth that trees access water.
idbh.type      = myidbhtype   # Type of DBH class
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



#----- Loading some packages and scripts. -------------------------------------------------#
source(file.path(srcdir,"load.everything.r"))
#------------------------------------------------------------------------------------------#



#----- Set how many formats we must output. -----------------------------------------------#
outform = tolower(outform)
nout    = length (outform)
#------------------------------------------------------------------------------------------#


#----- Avoid unecessary and extremely annoying beeps. -------------------------------------#
options(locatorBell=FALSE)
#------------------------------------------------------------------------------------------#



#----- Load observations. -----------------------------------------------------------------#
obsrfile = paste(srcdir,"LBA_MIP.v8.RData",sep="/")
load(file=obsrfile)
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
   outpref = paste(outmain,"monthly",sep="/")
   lieu    = thispoi$lieu
   iata    = thispoi$iata
   suffix  = thispoi$iata
   yeara   = thispoi$yeara
   yearz   = thispoi$yearz
   meszz   = thispoi$monz
   #---------------------------------------------------------------------------------------#



   #----- Create the directories in case they don't exist. --------------------------------#
   if (! file.exists(outmain)) dir.create(outmain)
   if (! file.exists(outpref)) dir.create(outpref)
   #---------------------------------------------------------------------------------------#



   #----- Decide how frequently the cohort-level variables should be saved. ---------------#
   if (yearend - yearbeg + 1 <= nyears.long){
      sasmonth   = sasmonth.short
      emean.line = TRUE
   }else{
      sasmonth   = sasmonth.long
      emean.line = FALSE
   }#end if
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Find the total number of months that can be loaded this time.                     #
   #---------------------------------------------------------------------------------------#
   ntimes     = (yearz-yeara-1)*12+meszz+(12-monthbeg+1)
   #---------------------------------------------------------------------------------------#



   #----- Print a banner to entretain the user. -------------------------------------------#
   cat(" + Post-processing output from ",lieu,"...","\n")
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      Make the RData file name, then we check whether we must read the files again     #
   # or use the stored RData.                                                              #
   #---------------------------------------------------------------------------------------#
   path.data  = paste(here,place,"rdata_month",sep="/")
   if (! file.exists(path.data)) dir.create(path.data)
   ed22.rdata = paste(path.data,paste(place,"RData",sep="."),sep="/")
   if (reload.data && file.exists(ed22.rdata)){
      #----- Load the modelled dataset. ---------------------------------------------------#
      cat("   - Loading previous session...","\n")
      load(ed22.rdata)
      tresume = datum$ntimes + 1
      datum   = update.monthly( new.ntimes = ntimes 
                              , old.datum  = datum
                              , montha     = monthbeg
                              , yeara      = yeara
                              , inpref     = inpref
                              , slz.min    = slz.min
                              )#end update.monthly
   }else{
      cat("   - Starting new session...","\n")
      tresume    = 1
      datum      = create.monthly( ntimes  = ntimes
                                 , montha  = monthbeg
                                 , yeara   = yeara
                                 , inpref  = inpref
                                 , slz.min = slz.min
                                 )#end create.monthly
   }#end if
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Make a list with the time span of each drought so we can plot rectangles showing  #
   # the drought.                                                                          #
   #---------------------------------------------------------------------------------------#
   drought = list()
   year    = drought.yeara
   ndrought = length(months.drought)
   n        = 0
   overyear = months.drought[1] > months.drought[ndrought]
   for (year in seq(from=drought.yeara,to=drought.yearz-as.integer(overyear),by=1)){
      n             = n + 1
      
      #----- Define the beginning and the end of the drought. -----------------------------#
      month.start   = months.drought[1]
      month.end     = 1 + (months.drought[ndrought] %% 12)
      year.start    = year
      year.end      = year + as.integer(month.end == 1) + 1

      drought.whena = chron(dates=paste(month.start,1,year.start,sep="/"))
      drought.whenz = chron(dates=paste(month.end  ,1,year.end  ,sep="/"))
      drought[[n]]  = c(drought.whena,drought.whenz)
   }#end for
   #----- ndrought becomes the number of blocks with drought. -----------------------------#
   ndrought = length(drought)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Check whether we have anything to update.                                         #
   #---------------------------------------------------------------------------------------#
   complete = tresume > ntimes
   #---------------------------------------------------------------------------------------#



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
   #     Loop over all times in case there is anything new to be read.                     #
   #---------------------------------------------------------------------------------------#
   if (! complete){

      #------------------------------------------------------------------------------------#
      #     This function will read the files.                                             #
      #------------------------------------------------------------------------------------#
      datum = read.q.files(datum=datum,ntimes=ntimes,tresume=tresume,sasmonth=sasmonth)
      #------------------------------------------------------------------------------------#

      #------ Save the data to the R object. ----------------------------------------------#
      cat(" + Saving data to ",basename(ed22.rdata),"...","\n")
      save(datum,file=ed22.rdata)
      #------------------------------------------------------------------------------------#
   }#end if (! complete)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      Define a suitable scale for those time series that uses datum$tomonth...             #
   #---------------------------------------------------------------------------------------#
   whenplot6 = pretty.time(datum$tomonth,n=6)
   whenplot8 = pretty.time(datum$tomonth,n=8)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      Define a suitable scale for diurnal cycle...                                     #
   #---------------------------------------------------------------------------------------#
   thisday = seq(from=0,to=ndcycle,by=1) * 24 / ndcycle
   uplot = list()
   uplot$levels = c(0,4,8,12,16,20,24)
   uplot$n      = 7
   uplot$scale  = "hours"
   uplot$padj   = rep(0,times=uplot$n)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Define a suitable scale for soil profile layers...                               #
   #---------------------------------------------------------------------------------------#
   znice  = -pretty.log(-slz,n=8)
   znice  = sort(c(znice,slz[1],slz[nzg]))
   sel    = znice >= slz[1] & znice <= slz[nzg]
   znice  = znice[sel]
   zat    = -log(-znice)
   nznice = length(znice)
   znice  = sprintf("%.2f",znice)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Define a suitable scale for monthly means...                                     #
   #---------------------------------------------------------------------------------------#
   montmont  = seq(from=1,to=12,by=1)
   mplot        = list()
   mplot$levels = montmont
   mplot$labels = capwords(mon2mmm(montmont))
   mplot$n      = 12
   mplot$scale  = "months"
   mplot$padj   = rep(0,times=mplot$n)
   #---------------------------------------------------------------------------------------#




   #----- Make some shorter versions of some variables. -----------------------------------#
   mfac   = datum$month
   emean  = datum$emean
   emsqu  = datum$emsqu
   qmean  = datum$qmean
   qmsqu  = datum$qmsqu
   szpft  = datum$szpft
   lu     = datum$lu
   patch  = datum$patch
   cohort = datum$cohort
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Here we find the monthly means for month, then compute the standard deviation.   #
   #---------------------------------------------------------------------------------------#
   cat ("    - Finding the monthly mean...","\n")
   cat ("      * Aggregating the monthly mean...","\n")
   mmean               = list()
   mmean$fast.soil.c   = tapply(X=emean$fast.soil.c  ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$slow.soil.c   = tapply(X=emean$slow.soil.c  ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$struct.soil.c = tapply(X=emean$struct.soil.c,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$het.resp      = tapply(X=emean$het.resp     ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$cwd.resp      = tapply(X=emean$cwd.resp     ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$gpp           = tapply(X=emean$gpp          ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$npp           = tapply(X=emean$npp          ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$plant.resp    = tapply(X=emean$plant.resp   ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$leaf.resp     = tapply(X=emean$leaf.resp    ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$root.resp     = tapply(X=emean$root.resp    ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$growth.resp   = tapply(X=emean$growth.resp  ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$reco          = tapply(X=emean$reco         ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$mco           = tapply(X=emean$mco          ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$cba           = tapply(X=emean$cba          ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$cbalight      = tapply(X=emean$cbalight     ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$cbamoist      = tapply(X=emean$cbamoist     ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$cbarel        = tapply(X=emean$cbarel       ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$ldrop         = tapply(X=emean$ldrop        ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$nep           = tapply(X=emean$nep          ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$nee           = tapply(X=emean$nee          ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$cflxca        = tapply(X=emean$cflxca       ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$cflxst        = tapply(X=emean$cflxst       ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$ustar         = tapply(X=emean$ustar        ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$atm.vels      = tapply(X=emean$atm.vels     ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$atm.prss      = tapply(X=emean$atm.prss     ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$atm.temp      = tapply(X=emean$atm.temp     ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$atm.shv       = tapply(X=emean$atm.shv      ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$atm.vpd       = tapply(X=emean$atm.vpd      ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$atm.co2       = tapply(X=emean$atm.co2      ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$can.prss      = tapply(X=emean$can.prss     ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$can.temp      = tapply(X=emean$can.temp     ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$can.co2       = tapply(X=emean$can.co2      ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$can.shv       = tapply(X=emean$can.shv      ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$can.vpd       = tapply(X=emean$can.vpd      ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$gnd.temp      = tapply(X=emean$gnd.temp     ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$gnd.shv       = tapply(X=emean$gnd.shv      ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$leaf.temp     = tapply(X=emean$leaf.temp    ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$leaf.vpd      = tapply(X=emean$leaf.vpd     ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$wood.temp     = tapply(X=emean$wood.temp    ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$hflxca        = tapply(X=emean$hflxca       ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$qwflxca       = tapply(X=emean$qwflxca      ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$hflxgc        = tapply(X=emean$hflxgc       ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$hflxlc        = tapply(X=emean$hflxlc       ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$hflxwc        = tapply(X=emean$hflxwc       ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$wflxca        = tapply(X=emean$wflxca       ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$wflxgc        = tapply(X=emean$wflxgc       ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$wflxlc        = tapply(X=emean$wflxlc       ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$wflxwc        = tapply(X=emean$wflxwc       ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$evap          = tapply(X=emean$evap         ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$transp        = tapply(X=emean$transp       ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$et            = tapply(X=emean$et           ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$wue           = tapply(X=emean$wue          ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$rain          = tapply(X=emean$rain         ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$fs.open       = tapply(X=emean$fs.open      ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$rshort        = tapply(X=emean$rshort       ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$rshort.beam   = tapply(X=emean$rshort.beam  ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$rshort.diff   = tapply(X=emean$rshort.diff  ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$rshortup      = tapply(X=emean$rshortup     ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$rshort.gnd    = tapply(X=emean$rshort.gnd   ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$rlong         = tapply(X=emean$rlong        ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$rlong.gnd     = tapply(X=emean$rlong.gnd    ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$rlongup       = tapply(X=emean$rlongup      ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$par.tot       = tapply(X=emean$par.tot      ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$par.beam      = tapply(X=emean$par.beam     ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$par.diff      = tapply(X=emean$par.diff     ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$par.gnd       = tapply(X=emean$par.gnd      ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$parup         = tapply(X=emean$parup        ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$rnet          = tapply(X=emean$rnet         ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$albedo        = tapply(X=emean$albedo       ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$albedo.beam   = tapply(X=emean$albedo.beam  ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$albedo.diff   = tapply(X=emean$albedo.diff  ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$rlong.albedo  = tapply(X=emean$rlong.albedo ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$nplant        = tapply(X=emean$nplant       ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$agb           = tapply(X=emean$agb          ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$lai           = tapply(X=emean$lai          ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$wai           = tapply(X=emean$wai          ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$tai           = tapply(X=emean$tai          ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$area          = tapply(X=emean$area         ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$workload      = tapply(X=emean$workload     ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$specwork      = tapply(X=emean$specwork     ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$demand        = tapply(X=emean$demand       ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$supply        = tapply(X=emean$supply       ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$paw           = tapply(X=emean$paw          ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$smpot         = tapply(X=emean$smpot        ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$leaf.gsw      = tapply(X=emean$leaf.gsw     ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$leaf.gbw      = tapply(X=emean$leaf.gbw     ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$wood.gbw      = tapply(X=emean$wood.gbw     ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$npat.global   = tapply(X=emean$npat.global  ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$ncoh.global   = tapply(X=emean$ncoh.global  ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$water.deficit = tapply(X=emean$water.deficit,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$i.gpp         = tapply(X=emean$i.gpp        ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$i.npp         = tapply(X=emean$i.npp        ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$i.plant.resp  = tapply(X=emean$i.plant.resp ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$i.mco         = tapply(X=emean$i.mco        ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$i.cba         = tapply(X=emean$i.cba        ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$i.cbalight    = tapply(X=emean$i.cbalight   ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$i.cbamoist    = tapply(X=emean$i.cbamoist   ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$i.transp      = tapply(X=emean$i.transp     ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$i.wflxlc      = tapply(X=emean$i.wflxlc     ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmean$i.hflxlc      = tapply(X=emean$i.hflxlc     ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   #----- Soil variables. -----------------------------------------------------------------#
   mmean$soil.temp     = qapply(X=emean$soil.temp    ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   mmean$soil.water    = qapply(X=emean$soil.water   ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   mmean$soil.mstpot   = qapply(X=emean$soil.mstpot  ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   #----- Find the mean sum of squares. ---------------------------------------------------#
   cat ("      * Aggregating the monthly mean sum of squares...","\n")
   mmsqu               = list()
   mmsqu$gpp           = tapply(X=emsqu$gpp          , INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmsqu$plant.resp    = tapply(X=emsqu$plant.resp   , INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmsqu$het.resp      = tapply(X=emsqu$het.resp     , INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmsqu$cwd.resp      = tapply(X=emsqu$cwd.resp     , INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmsqu$cflxca        = tapply(X=emsqu$cflxca       , INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmsqu$cflxst        = tapply(X=emsqu$cflxst       , INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmsqu$hflxca        = tapply(X=emsqu$hflxca       , INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmsqu$hflxlc        = tapply(X=emsqu$hflxlc       , INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmsqu$hflxwc        = tapply(X=emsqu$hflxwc       , INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmsqu$hflxgc        = tapply(X=emsqu$hflxgc       , INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmsqu$wflxca        = tapply(X=emsqu$wflxca       , INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmsqu$qwflxca       = tapply(X=emsqu$qwflxca      , INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmsqu$wflxlc        = tapply(X=emsqu$wflxlc       , INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmsqu$wflxwc        = tapply(X=emsqu$wflxwc       , INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmsqu$wflxgc        = tapply(X=emsqu$wflxgc       , INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmsqu$evap          = tapply(X=emsqu$evap         , INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmsqu$transp        = tapply(X=emsqu$transp       , INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmsqu$ustar         = tapply(X=emsqu$ustar        , INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmsqu$albedo        = tapply(X=emsqu$albedo       , INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmsqu$rshortup      = tapply(X=emsqu$rshortup     , INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmsqu$rlongup       = tapply(X=emsqu$rlongup      , INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmsqu$parup         = tapply(X=emsqu$parup        , INDEX=mfac,FUN=mean,na.rm=TRUE)
   mmsqu$rnet          = tapply(X=emsqu$rnet         , INDEX=mfac,FUN=mean,na.rm=TRUE)
   #---------------------------------------------------------------------------------------#
   #   Here we convert the sum of squares into standard deviation. The standard devi-      #
   # ation can be written in two different ways, and we will use the latter because it     #
   # doesn't require previous knowledge of the mean.                                       #
   #              ____________________          _____________________________________      #
   #             / SUM_i[X_i - Xm]^2           /  / SUM_i[X_i^2]         \      1          #
   # sigma = \  /  ------------------   =  \  /  |  ------------  - Xm^2 | ---------       #
   #          \/       N - 1                \/    \      N               /   1 - 1/N       #
   #                                                                                       #
   # srnonm1 is the square root of 1 / (1 - 1/N)                                           #
   #     Find the standard deviation.                                                      #
   #---------------------------------------------------------------------------------------#
   cat ("      * Finding the standard deviation...","\n")
   srnorm1 = sqrt(1./(1. - 1. / datum$montable))
   srnorm1[!is.finite(srnorm1)] = 0.

   msdev = list()
   msdev$gpp           = sqrt( mmsqu$gpp        - mmean$gpp        ^ 2 ) * srnorm1 
   msdev$plant.resp    = sqrt( mmsqu$plant.resp - mmean$plant.resp ^ 2 ) * srnorm1 
   msdev$het.resp      = sqrt( mmsqu$het.resp   - mmean$het.resp   ^ 2 ) * srnorm1 
   msdev$cwd.resp      = sqrt( mmsqu$cwd.resp   - mmean$cwd.resp   ^ 2 ) * srnorm1 
   msdev$cflxca        = sqrt( mmsqu$cflxca     - mmean$cflxca     ^ 2 ) * srnorm1 
   msdev$cflxst        = sqrt( mmsqu$cflxst     - mmean$cflxst     ^ 2 ) * srnorm1 
   msdev$hflxca        = sqrt( mmsqu$hflxca     - mmean$hflxca     ^ 2 ) * srnorm1 
   msdev$hflxlc        = sqrt( mmsqu$hflxlc     - mmean$hflxlc     ^ 2 ) * srnorm1 
   msdev$hflxwc        = sqrt( mmsqu$hflxwc     - mmean$hflxwc     ^ 2 ) * srnorm1 
   msdev$hflxgc        = sqrt( mmsqu$hflxgc     - mmean$hflxgc     ^ 2 ) * srnorm1 
   msdev$wflxca        = sqrt( mmsqu$wflxca     - mmean$wflxca     ^ 2 ) * srnorm1 
   msdev$qwflxca       = sqrt( mmsqu$qwflxca    - mmean$qwflxca    ^ 2 ) * srnorm1 
   msdev$wflxlc        = sqrt( mmsqu$wflxlc     - mmean$wflxlc     ^ 2 ) * srnorm1 
   msdev$wflxwc        = sqrt( mmsqu$wflxwc     - mmean$wflxwc     ^ 2 ) * srnorm1 
   msdev$wflxgc        = sqrt( mmsqu$wflxgc     - mmean$wflxgc     ^ 2 ) * srnorm1 
   msdev$evap          = sqrt( mmsqu$evap       - mmean$evap       ^ 2 ) * srnorm1 
   msdev$transp        = sqrt( mmsqu$transp     - mmean$transp     ^ 2 ) * srnorm1 
   msdev$ustar         = sqrt( mmsqu$ustar      - mmean$ustar      ^ 2 ) * srnorm1 
   msdev$albedo        = sqrt( mmsqu$albedo     - mmean$albedo     ^ 2 ) * srnorm1 
   msdev$rshortup      = sqrt( mmsqu$rshortup   - mmean$rshortup   ^ 2 ) * srnorm1 
   msdev$rlongup       = sqrt( mmsqu$rlongup    - mmean$rlongup    ^ 2 ) * srnorm1 
   msdev$parup         = sqrt( mmsqu$parup      - mmean$parup      ^ 2 ) * srnorm1 
   msdev$rnet          = sqrt( mmsqu$rnet       - mmean$rnet       ^ 2 ) * srnorm1 
   #---------------------------------------------------------------------------------------#
   #     Set standard deviations that became NaN to 0.  This usually happens when we run   #
   # the post-processing script when the simulation hasn't run for more than 2 years.  We  #
   # can't find the standard deviation because the number of degrees of freedom becomes 0. #
   #---------------------------------------------------------------------------------------#
   msdev$gpp           [! is.finite(msdev$gpp       )] = 0.
   msdev$plant.resp    [! is.finite(msdev$plant.resp)] = 0.
   msdev$leaf.resp     [! is.finite(msdev$leaf.resp )] = 0.
   msdev$root.resp     [! is.finite(msdev$root.resp )] = 0.
   msdev$het.resp      [! is.finite(msdev$het.resp  )] = 0.
   msdev$cwd.resp      [! is.finite(msdev$cwd.resp  )] = 0.
   msdev$cflxca        [! is.finite(msdev$cflxca    )] = 0.
   msdev$cflxst        [! is.finite(msdev$cflxst    )] = 0.
   msdev$hflxca        [! is.finite(msdev$hflxca    )] = 0.
   msdev$hflxlc        [! is.finite(msdev$hflxlc    )] = 0.
   msdev$hflxwc        [! is.finite(msdev$hflxwc    )] = 0.
   msdev$hflxgc        [! is.finite(msdev$hflxgc    )] = 0.
   msdev$wflxca        [! is.finite(msdev$wflxca    )] = 0.
   msdev$qwflxca       [! is.finite(msdev$qwflxca   )] = 0.
   msdev$wflxlc        [! is.finite(msdev$wflxlc    )] = 0.
   msdev$wflxwc        [! is.finite(msdev$wflxwc    )] = 0.
   msdev$wflxgc        [! is.finite(msdev$wflxgc    )] = 0.
   msdev$transp        [! is.finite(msdev$transp    )] = 0.
   msdev$ustar         [! is.finite(msdev$ustar     )] = 0.
   msdev$albedo        [! is.finite(msdev$albedo    )] = 0.
   msdev$rshortup      [! is.finite(msdev$rshortup  )] = 0.
   msdev$rlongup       [! is.finite(msdev$rlongup   )] = 0.
   msdev$parup         [! is.finite(msdev$parup     )] = 0.
   msdev$rnet          [! is.finite(msdev$rnet      )] = 0.
   #---------------------------------------------------------------------------------------#
   #     Estimate the standard deviation of NEE, REco, and evaporation.                    #
   #---------------------------------------------------------------------------------------#
   msdev$nee  = sqrt(msdev$cflxca^2     + msdev$cflxst^2                        )
   msdev$reco = sqrt(msdev$plant.resp^2 + msdev$het.resp^2                      )
   msdev$evap = sqrt(msdev$wflxgc^2     + msdev$wflxlc^2    + msdev$wflxwc^2    )
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Here we find the Mean diurnal cycle for each month, then compute the standard    #
   # deviation.                                                                            #
   #---------------------------------------------------------------------------------------#
   cat ("    - Aggregating the monthly mean of the diurnal cycle...","\n")
   umean              = list()
   umean$gpp          = qapply(X=qmean$gpp          ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   umean$npp          = qapply(X=qmean$npp          ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   umean$plant.resp   = qapply(X=qmean$plant.resp   ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   umean$leaf.resp    = qapply(X=qmean$leaf.resp    ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   umean$root.resp    = qapply(X=qmean$root.resp    ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   umean$growth.resp  = qapply(X=qmean$growth.resp  ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   umean$het.resp     = qapply(X=qmean$het.resp     ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   umean$cwd.resp     = qapply(X=qmean$cwd.resp     ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   umean$nep          = qapply(X=qmean$nep          ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   umean$nee          = qapply(X=qmean$nee          ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   umean$reco         = qapply(X=qmean$reco         ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   umean$cflxca       = qapply(X=qmean$cflxca       ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   umean$cflxst       = qapply(X=qmean$cflxst       ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   umean$hflxca       = qapply(X=qmean$hflxca       ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   umean$hflxlc       = qapply(X=qmean$hflxlc       ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   umean$hflxwc       = qapply(X=qmean$hflxwc       ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   umean$hflxgc       = qapply(X=qmean$hflxgc       ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   umean$qwflxca      = qapply(X=qmean$qwflxca      ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   umean$wflxca       = qapply(X=qmean$wflxca       ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   umean$wflxlc       = qapply(X=qmean$wflxlc       ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   umean$wflxwc       = qapply(X=qmean$wflxwc       ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   umean$wflxgc       = qapply(X=qmean$wflxgc       ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   umean$evap         = qapply(X=qmean$evap         ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   umean$transp       = qapply(X=qmean$transp       ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   umean$atm.temp     = qapply(X=qmean$atm.temp     ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   umean$can.temp     = qapply(X=qmean$can.temp     ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   umean$leaf.temp    = qapply(X=qmean$leaf.temp    ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   umean$wood.temp    = qapply(X=qmean$wood.temp    ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   umean$gnd.temp     = qapply(X=qmean$gnd.temp     ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   umean$atm.shv      = qapply(X=qmean$atm.shv      ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   umean$can.shv      = qapply(X=qmean$can.shv      ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   umean$gnd.shv      = qapply(X=qmean$gnd.shv      ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   umean$atm.vpd      = qapply(X=qmean$atm.vpd      ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   umean$can.vpd      = qapply(X=qmean$can.vpd      ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   umean$leaf.vpd     = qapply(X=qmean$leaf.vpd     ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   umean$atm.co2      = qapply(X=qmean$atm.co2      ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   umean$can.co2      = qapply(X=qmean$can.co2      ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   umean$atm.prss     = qapply(X=qmean$atm.prss     ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   umean$can.prss     = qapply(X=qmean$can.prss     ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   umean$atm.vels     = qapply(X=qmean$atm.vels     ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   umean$ustar        = qapply(X=qmean$ustar        ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   umean$fs.open      = qapply(X=qmean$fs.open      ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   umean$rain         = qapply(X=qmean$rain         ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   umean$rshort       = qapply(X=qmean$rshort       ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   umean$rshort.beam  = qapply(X=qmean$rshort.beam  ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   umean$rshort.diff  = qapply(X=qmean$rshort.diff  ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   umean$rshort.gnd   = qapply(X=qmean$rshort.gnd   ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   umean$rshortup     = qapply(X=qmean$rshortup     ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   umean$rlong        = qapply(X=qmean$rlong        ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   umean$rlong.gnd    = qapply(X=qmean$rlong.gnd    ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   umean$rlongup      = qapply(X=qmean$rlongup      ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   umean$par.tot      = qapply(X=qmean$par.tot      ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   umean$par.beam     = qapply(X=qmean$par.beam     ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   umean$par.diff     = qapply(X=qmean$par.diff     ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   umean$par.gnd      = qapply(X=qmean$par.gnd      ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   umean$parup        = qapply(X=qmean$parup        ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   umean$rnet         = qapply(X=qmean$rnet         ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   umean$albedo       = qapply(X=qmean$albedo       ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   umean$albedo.beam  = qapply(X=qmean$albedo.beam  ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   umean$albedo.diff  = qapply(X=qmean$albedo.diff  ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   umean$rlong.albedo = qapply(X=qmean$rlong.albedo ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   umean$leaf.gsw     = qapply(X=qmean$leaf.gsw     ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   umean$leaf.gbw     = qapply(X=qmean$leaf.gbw     ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   umean$wood.gbw     = qapply(X=qmean$wood.gbw     ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   #----- Find the mean sum of squares. ---------------------------------------------------#
   cat ("    - Aggregating the monthly mean sum of squares...","\n")
   umsqu              = list()
   umsqu$gpp          = qapply(X=qmsqu$gpp          ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   umsqu$npp          = qapply(X=qmsqu$npp          ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   umsqu$plant.resp   = qapply(X=qmsqu$plant.resp   ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   umsqu$het.resp     = qapply(X=qmsqu$het.resp     ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   umsqu$cwd.resp     = qapply(X=qmsqu$cwd.resp     ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   umsqu$nep          = qapply(X=qmsqu$nep          ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   umsqu$cflxca       = qapply(X=qmsqu$cflxca       ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   umsqu$cflxst       = qapply(X=qmsqu$cflxst       ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   umsqu$hflxca       = qapply(X=qmsqu$hflxca       ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   umsqu$hflxlc       = qapply(X=qmsqu$hflxlc       ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   umsqu$hflxwc       = qapply(X=qmsqu$hflxwc       ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   umsqu$hflxgc       = qapply(X=qmsqu$hflxgc       ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   umsqu$qwflxca      = qapply(X=qmsqu$qwflxca      ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   umsqu$wflxca       = qapply(X=qmsqu$wflxca       ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   umsqu$wflxlc       = qapply(X=qmsqu$wflxlc       ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   umsqu$wflxwc       = qapply(X=qmsqu$wflxwc       ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   umsqu$wflxgc       = qapply(X=qmsqu$wflxgc       ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   umsqu$transp       = qapply(X=qmsqu$transp       ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   umsqu$ustar        = qapply(X=qmsqu$ustar        ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   umsqu$albedo       = qapply(X=qmsqu$albedo       ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   umsqu$rshortup     = qapply(X=qmsqu$rshortup     ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   umsqu$rlongup      = qapply(X=qmsqu$rlongup      ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   umsqu$parup        = qapply(X=qmsqu$parup        ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   umsqu$rnet         = qapply(X=qmsqu$rnet         ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   #---------------------------------------------------------------------------------------#
   #   Here we convert the sum of squares into standard deviation. The standard devi-      #
   # ation can be written in two different ways, and we will use the latter because it     #
   # doesn't require previous knowledge of the mean.                                       #
   #              ____________________          _____________________________________      #
   #             / SUM_i[X_i - Xm]^2           /  / SUM_i[X_i^2]         \      1          #
   # sigma = \  /  ------------------   =  \  /  |  ------------  - Xm^2 | ---------       #
   #          \/       N - 1                \/    \      N               /   1 - 1/N       #
   #                                                                                       #
   # srnonm1 is the square root of 1 / (1 - 1/N)                                           #
   #     Find the standard deviation.                                                      #
   #---------------------------------------------------------------------------------------#
   cat ("    - Finding the standard deviation...","\n")
   srnorm1 = sqrt(1./(1. - 1. / datum$moncnt))
   srnorm1[!is.finite(srnorm1)] = 0.
   usdev              = list()
   usdev$gpp          = sqrt(umsqu$gpp           - umean$gpp           ^ 2) * srnorm1
   usdev$npp          = sqrt(umsqu$npp           - umean$npp           ^ 2) * srnorm1
   usdev$plant.resp   = sqrt(umsqu$plant.resp    - umean$plant.resp    ^ 2) * srnorm1
   usdev$het.resp     = sqrt(umsqu$het.resp      - umean$het.resp      ^ 2) * srnorm1
   usdev$cwd.resp     = sqrt(umsqu$cwd.resp      - umean$cwd.resp      ^ 2) * srnorm1
   usdev$nep          = sqrt(umsqu$nep           - umean$nep           ^ 2) * srnorm1
   usdev$cflxca       = sqrt(umsqu$cflxca        - umean$cflxca        ^ 2) * srnorm1
   usdev$cflxst       = sqrt(umsqu$cflxst        - umean$cflxst        ^ 2) * srnorm1
   usdev$hflxca       = sqrt(umsqu$hflxca        - umean$hflxca        ^ 2) * srnorm1
   usdev$hflxlc       = sqrt(umsqu$hflxlc        - umean$hflxlc        ^ 2) * srnorm1
   usdev$hflxwc       = sqrt(umsqu$hflxwc        - umean$hflxwc        ^ 2) * srnorm1
   usdev$hflxgc       = sqrt(umsqu$hflxgc        - umean$hflxgc        ^ 2) * srnorm1
   usdev$qwflxca      = sqrt(umsqu$qwflxca       - umean$qwflxca       ^ 2) * srnorm1
   usdev$wflxca       = sqrt(umsqu$wflxca        - umean$wflxca        ^ 2) * srnorm1
   usdev$wflxlc       = sqrt(umsqu$wflxlc        - umean$wflxlc        ^ 2) * srnorm1
   usdev$wflxwc       = sqrt(umsqu$wflxwc        - umean$wflxwc        ^ 2) * srnorm1
   usdev$wflxgc       = sqrt(umsqu$wflxgc        - umean$wflxgc        ^ 2) * srnorm1
   usdev$transp       = sqrt(umsqu$transp        - umean$transp        ^ 2) * srnorm1
   usdev$ustar        = sqrt(umsqu$ustar         - umean$ustar         ^ 2) * srnorm1
   usdev$albedo       = sqrt(umsqu$albedo        - umean$albedo        ^ 2) * srnorm1
   usdev$rshortup     = sqrt(umsqu$rshortup      - umean$rshortup      ^ 2) * srnorm1
   usdev$rlongup      = sqrt(umsqu$rlongup       - umean$rlongup       ^ 2) * srnorm1
   usdev$parup        = sqrt(umsqu$parup         - umean$parup         ^ 2) * srnorm1
   usdev$rnet         = sqrt(umsqu$rnet          - umean$rnet          ^ 2) * srnorm1
   #---------------------------------------------------------------------------------------#
   #     Set standard deviations that became NaN to 0.  This usually happens when we run   #
   # the post-processing script when the simulation hasn't run for more than 2 years.  We  #
   # can't find the standard deviation because the number of degrees of freedom becomes 0. #
   #---------------------------------------------------------------------------------------#
   usdev$gpp          [! is.finite(usdev$gpp          )] = 0.0
   usdev$npp          [! is.finite(usdev$npp          )] = 0.0
   usdev$plant.resp   [! is.finite(usdev$plant.resp   )] = 0.0
   usdev$het.resp     [! is.finite(usdev$het.resp     )] = 0.0
   usdev$cwd.resp     [! is.finite(usdev$cwd.resp     )] = 0.0
   usdev$nep          [! is.finite(usdev$nep          )] = 0.0
   usdev$cflxca       [! is.finite(usdev$cflxca       )] = 0.0
   usdev$cflxst       [! is.finite(usdev$cflxst       )] = 0.0
   usdev$hflxca       [! is.finite(usdev$hflxca       )] = 0.0
   usdev$hflxlc       [! is.finite(usdev$hflxlc       )] = 0.0
   usdev$hflxwc       [! is.finite(usdev$hflxwc       )] = 0.0
   usdev$hflxgc       [! is.finite(usdev$hflxgc       )] = 0.0
   usdev$qwflxca      [! is.finite(usdev$qwflxca      )] = 0.0
   usdev$wflxca       [! is.finite(usdev$wflxca       )] = 0.0
   usdev$wflxlc       [! is.finite(usdev$wflxlc       )] = 0.0
   usdev$wflxwc       [! is.finite(usdev$wflxwc       )] = 0.0
   usdev$wflxgc       [! is.finite(usdev$wflxgc       )] = 0.0
   usdev$transp       [! is.finite(usdev$transp       )] = 0.0
   usdev$ustar        [! is.finite(usdev$ustar        )] = 0.0
   usdev$albedo       [! is.finite(usdev$albedo       )] = 0.0
   usdev$rshortup     [! is.finite(usdev$rshortup     )] = 0.0
   usdev$rlongup      [! is.finite(usdev$rlongup      )] = 0.0
   usdev$parup        [! is.finite(usdev$parup        )] = 0.0
   usdev$rnet         [! is.finite(usdev$rnet         )] = 0.0
   #---------------------------------------------------------------------------------------#
   #      Estimate NPP and NEE standard deviation.                                         #
   #---------------------------------------------------------------------------------------#
   usdev$nee  = sqrt(usdev$cflxca^2     + usdev$cflxst^2                        )
   usdev$reco = sqrt(usdev$plant.resp^2 + usdev$het.resp^2                      )
   usdev$evap = sqrt(usdev$wflxgc^2     + usdev$wflxlc^2    + usdev$wflxwc^2    )
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Remove all elements of the DBH/PFT class that do not have a single valid cohort   #
   # at any given time.                                                                    #
   #---------------------------------------------------------------------------------------#
   empty = szpft$nplant == 0
   szpft$agb               [empty] = NA
   szpft$lai               [empty] = NA
   szpft$wai               [empty] = NA
   szpft$tai               [empty] = NA
   szpft$ba                [empty] = NA
   szpft$gpp               [empty] = NA
   szpft$npp               [empty] = NA
   szpft$plant.resp        [empty] = NA
   szpft$mco               [empty] = NA
   szpft$cba               [empty] = NA
   szpft$cbalight          [empty] = NA
   szpft$cbamoist          [empty] = NA
   szpft$cbarel            [empty] = NA
   szpft$ldrop             [empty] = NA
   szpft$fs.open           [empty] = NA
   szpft$leaf.gbw          [empty] = NA
   szpft$leaf.gsw          [empty] = NA
   szpft$wood.gbw          [empty] = NA
   szpft$demand            [empty] = NA
   szpft$supply            [empty] = NA
   szpft$nplant            [empty] = NA
   szpft$mort              [empty] = NA
   szpft$dimort            [empty] = NA
   szpft$ncbmort           [empty] = NA
   szpft$growth            [empty] = NA
   szpft$bdead             [empty] = NA
   szpft$bleaf             [empty] = NA
   szpft$broot             [empty] = NA
   szpft$bsapwood          [empty] = NA
   szpft$bstorage          [empty] = NA
   szpft$bseeds            [empty] = NA
   szpft$hflxlc            [empty] = NA
   szpft$wflxlc            [empty] = NA
   szpft$transp            [empty] = NA
   szpft$wue               [empty] = NA
   szpft$i.gpp             [empty] = NA
   szpft$i.npp             [empty] = NA
   szpft$i.plant.resp      [empty] = NA
   szpft$i.mco             [empty] = NA
   szpft$i.cba             [empty] = NA
   szpft$i.cbalight        [empty] = NA
   szpft$i.cbamoist        [empty] = NA
   szpft$i.transp          [empty] = NA
   szpft$i.wflxlc          [empty] = NA
   szpft$i.hflxlc          [empty] = NA
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Replace the mortality and recruitment exponential rates by the "interests" rates. #
   #---------------------------------------------------------------------------------------#
   szpft$mort      = 100. * (1.0 - exp(- szpft$mort     )      )
   szpft$dimort    = 100. * (1.0 - exp(- szpft$dimort   )      )
   szpft$ncbmort   = 100. * (1.0 - exp(- szpft$ncbmort  )      )
   szpft$recrpft   = 100. * (      exp(  szpft$recr     ) - 1.0)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      Find the patch density function for all patch-level data.                        #
   #---------------------------------------------------------------------------------------#
   cat ("    - Finding the distribution function of patch properties...","\n")
   patchpdf = list()
   for (pp in 1:nplotpatch){
      this        = plotpatch[[pp]]
      vname       = this$vnam
      col.scheme  = get(this$col.scheme)(n=ncolsfc)
      emean.area  = patch$area
      emean.vname = patch[[vname]]
      mmean.area  = tapply(X=emean.area ,INDEX=mfac,FUN=unlist)
      mmean.vname = tapply(X=emean.vname,INDEX=mfac,FUN=unlist)

      #----- Find the range for which we find the density function. -----------------------#
      low.vname   = min(unlist(emean.vname),na.rm=TRUE)
      high.vname  = max(unlist(emean.vname),na.rm=TRUE)
      #------------------------------------------------------------------------------------#



      #----- Find the density function for each time. -------------------------------------#
      edfun.now   = mapply( FUN      = density.safe
                          , x        = emean.vname
                          , weights  = emean.area
                          , MoreArgs = list(n=n.density,from=low.vname,to=high.vname)
                          )#end mapply
      mdfun.now   = mapply( FUN      = density.safe
                          , x        = mmean.vname
                          , weights  = mmean.area
                          , MoreArgs = list(n=n.density,from=low.vname,to=high.vname)
                          )#end mapply
      #------------------------------------------------------------------------------------#




      #----- Save the density function. ---------------------------------------------------#
      edfun        = list()
      edfun$x      = chron(datum$when)
      edfun$y      = seq(from=low.vname,to=high.vname,length.out=n.density)
      edfun$z      = t(sapply(X=edfun.now["y",],FUN=cbind))
      #------------------------------------------------------------------------------------#




      #----- Save the density function. ---------------------------------------------------#
      mdfun        = list()
      mdfun$x      = montmont
      mdfun$y      = seq(from=low.vname,to=high.vname,length.out=n.density)
      mdfun$z      = t(sapply(X=mdfun.now["y",],FUN=cbind))
      #------------------------------------------------------------------------------------#



      #----- Remove tiny values (even with log scale values can be very hard to see. ------#
      bye         = is.finite(edfun$z) & edfun$z < 1.e-6 * max(unlist(edfun$z),na.rm=TRUE)
      edfun$z[bye] = NA
      #------------------------------------------------------------------------------------#


      #----- Remove tiny values (even with log scale values can be very hard to see. ------#
      bye         = is.finite(mdfun$z) & mdfun$z < 1.e-6 * max(unlist(mdfun$z),na.rm=TRUE)
      mdfun$z[bye] = NA
      #------------------------------------------------------------------------------------#
      patchpdf[[vname]] = list(edensity=edfun,mdensity=mdfun)
   }#end for 
   #---------------------------------------------------------------------------------------#




   #----- Find which PFTs, land uses and transitions we need to consider ------------------#
   pftave  = apply( X      = szpft$agb[,ndbh+1,]
                  , MARGIN = 2
                  , FUN    = mean
                  , na.rm  = TRUE
                  )#end apply
   luave   = apply( X      = lu$agb 
                  , MARGIN = 2
                  , FUN    = mean
                  , na.rm  = TRUE
                  )#end apply
   distave = apply(X=lu$dist,MARGIN=c(2,3),FUN=mean)
   selpft  = is.finite(pftave ) & pftave  > 0.
   sellu   = is.finite(luave  ) & luave   > 0.
   seldist = is.finite(distave) & distave > 0.
   n.selpft  = sum(selpft )
   n.sellu   = sum(sellu  )
   n.seldist = sum(seldist)
   #---------------------------------------------------------------------------------------#






   #=======================================================================================#
   #=======================================================================================#
   #=======================================================================================#
   #      Plotting section begins here...                                                  #
   #---------------------------------------------------------------------------------------#
   cat ("    - Plotting figures...","\n")
   #---------------------------------------------------------------------------------------#






   #---------------------------------------------------------------------------------------#
   #      Time series by PFT.                                                              #
   #---------------------------------------------------------------------------------------#
   for (v in 1:ntspftdbh){
      thistspft   = tspftdbh[[v]]
      vnam        = thistspft$vnam
      description = thistspft$desc
      unit        = thistspft$e.unit
      plog        = thistspft$plog
      plotit      = thistspft$pft

      #----- Check whether the user wants to have this variable plotted. ------------------#
      if (plotit && any(selpft)){

         #---------------------------------------------------------------------------------#
         #    Check whether the time series directory exists.  If not, create it.          #
         #---------------------------------------------------------------------------------#
         outdir = paste(outpref,"tspft",sep="/")
         if (! file.exists(outdir)) dir.create(outdir)
         cat("      + ",description," time series for all PFTs...","\n")

         #----- Load variable -------------------------------------------------------------#
         thisvar = szpft[[vnam]][,ndbh+1,]
         if (plog){
            #----- Eliminate non-positive values in case it is a log plot. ----------------#
            thisvar[thisvar <= 0] = NA
         }#end if
         #---------------------------------------------------------------------------------#



         #----- Loop over output formats. -------------------------------------------------#
         for (o in 1:nout){
            fichier = paste(outdir,"/",vnam,"-",suffix,".",outform[o],sep="")
            if(outform[o] == "x11"){
               X11(width=size$width,height=size$height,pointsize=ptsz)
            }else if(outform[o] == "png"){
               png(filename=fichier,width=size$width*depth,height=size$height*depth
                  ,pointsize=ptsz,res=depth)
            }else if(outform[o] == "eps"){
               postscript(file=fichier,width=size$width,height=size$height
                         ,pointsize=ptsz,paper=size$paper)
            }else if(outform[o] == "pdf"){
               pdf(file=fichier,onefile=FALSE
                  ,width=size$width,height=size$height,pointsize=ptsz,paper=size$paper)
            }#end if


            #------------------------------------------------------------------------------#
            #     Find the limit, make some room for the legend, and in case the field is  #
            # a constant, nudge the limits so the plot command will not complain.          #
            #------------------------------------------------------------------------------#
            ylimit = pretty.xylim(u=thisvar[,selpft],fracexp=scalleg,is.log=plog)
            if (plog){
               xylog    = "y"
               ydrought = c( exp(ylimit[1] * sqrt(ylimit[1]/ylimit[2]))
                           , exp(ylimit[2] * sqrt(ylimit[2]/ylimit[1]))
                           )#end c
            }else{
               xylog    = ""
               ydrought = c(ylimit[1] - 0.5 * diff(ylimit), ylimit[2] + 0.5 * diff(ylimit))
            }#end if
            #------------------------------------------------------------------------------#



            letitre = paste(description,lieu,sep=" - ")
            cols    = pft$colour[selpft]
            legs    = pft$name  [selpft]
            par(par.user)
            plot(x=datum$tomonth,y=thisvar[,1],type="n",main=letitre,ylim=ylimit
                ,xlab="Time",xaxt="n",ylab=unit,cex.main=0.7,log=xylog)
            axis(side=1,at=whenplot8$levels,labels=whenplot8$labels,padj=whenplot8$padj)

            if (drought.mark){
               for (n in 1:ndrought){
                  rect(xleft  = drought[[n]][1],ybottom = ydrought[1]
                      ,xright = drought[[n]][2],ytop    = ydrought[2]
                      ,col    = grid.colour,border=NA)
               }#end for
            }#end if
            
            if (plotgrid){ 
               abline(v=whenplot8$levels,h=axTicks(side=2),col=grid.colour,lty="solid")
            }#end if
            for (n in 1:(npft+1)){
               if (selpft[n]){
                  lines(datum$tomonth,thisvar[,n],type="l",col=pft$colour[n],lwd=lwidth)
               }#end if
            }#end for
            legend( x      = legwhere
                  , inset  = inset
                  , bg     = background
                  , legend = legs
                  , col    = cols
                  , lwd    = lwidth
                  , ncol   = min(pretty.box(n.selpft)$ncol,3)
                  , title  = expression(bold("Plant Functional Type"))
                  )#end legend

            if (outform[o] == "x11"){
               locator(n=1)
               dev.off()
            }else{
               dev.off()
            }#end if
            dummy = clean.tmp()
         } #end for outform
      }#end if (tseragbpft)
   } #end for tseries
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #      Time series by DBH, by PFT.                                                      #
   #---------------------------------------------------------------------------------------#
   #----- Find the PFTs to plot. ----------------------------------------------------------#
   pftuse  = which(apply(X=szpft$nplant,MARGIN=3,FUN=sum,na.rm=TRUE) > 0.)
   pftuse  = pftuse[pftuse != (npft+1)]
   for (v in 1:ntspftdbh){
      thistspftdbh   = tspftdbh[[v]]
      vnam           = thistspftdbh$vnam
      description    = thistspftdbh$desc
      unit           = thistspftdbh$e.unit
      plog           = thistspftdbh$plog
      plotit         = thistspftdbh$pftdbh
      
      #----- Load variable ----------------------------------------------------------------#
      thisvar = szpft[[vnam]]
      if (plog){
         xylog="y"
         thisvar[thisvar <= 0] = NA
      }else{
         xylog=""
      }#end if
      #----- Check whether the user wants to have this variable plotted. ------------------#
      if (plotit && length(pftuse) > 0 && any(is.finite(thisvar))){

         #---------------------------------------------------------------------------------#
         #    Check whether the time series directory exists.  If not, create it.          #
         #---------------------------------------------------------------------------------#
         outdir = paste(outpref,"tsdbh",sep="/")
         if (! file.exists(outdir)) dir.create(outdir)
         outvar = paste(outdir,vnam,sep="/")
         if (! file.exists(outvar)) dir.create(outvar)
         #---------------------------------------------------------------------------------#

         cat("      + ",description," time series for DBH class...","\n")


         #---------------------------------------------------------------------------------#
         #     Find the limit, make some room for the legend, and in case the field is a   #
         # constant, nudge the limits so the plot command will not complain.               #
         #---------------------------------------------------------------------------------#
         ylimit = pretty.xylim(u=thisvar[,,pftuse],fracexp=scalleg,is.log=plog)
         if (plog){
            xylog    = "y"
            ydrought = c( exp(ylimit[1] * sqrt(ylimit[1]/ylimit[2]))
                        , exp(ylimit[2] * sqrt(ylimit[2]/ylimit[1]))
                        )#end c
         }else{
            xylog    = ""
            ydrought = c(ylimit[1] - 0.5 * diff(ylimit), ylimit[2] + 0.5 * diff(ylimit))
         }#end if
         #---------------------------------------------------------------------------------#

         for (p in pftuse){

            cpp = substring(100+p,2,3)
            pftlab = paste("pft-",cpp,sep="")

            cat("        - ",pft$name[p],"\n")


            #----- Loop over output formats. ----------------------------------------------#
            for (o in 1:nout){
               fichier = paste(outvar,"/",vnam,"-",pftlab,"-",suffix,".",outform[o],sep="")
               if(outform[o] == "x11"){
                  X11(width=size$width,height=size$height,pointsize=ptsz)
               }else if(outform[o] == "png"){
                  png(filename=fichier,width=size$width*depth,height=size$height*depth
                     ,pointsize=ptsz,res=depth)
               }else if(outform[o] == "eps"){
                  postscript(file=fichier,width=size$width,height=size$height
                            ,pointsize=ptsz,paper=size$paper)
               }else if(outform[o] == "pdf"){
                  pdf(file=fichier,onefile=FALSE
                     ,width=size$width,height=size$height,pointsize=ptsz,paper=size$paper)
               }#end if

               letitre = paste(description,pft$name[p],lieu,sep=" - ")
               par(par.user)
               plot(x=datum$tomonth,y=thisvar[,1,p],type="n",main=letitre,ylim=ylimit
                   ,xlab="Time",xaxt="n",ylab=unit,cex.main=0.7,log=xylog)
               axis(side=1,at=whenplot8$levels,labels=whenplot8$labels,padj=whenplot8$padj)
               if (drought.mark){
                  for (n in 1:ndrought){
                     rect(xleft  = drought[[n]][1],ybottom = ydrought[1]
                         ,xright = drought[[n]][2],ytop    = ydrought[2]
                         ,col    = grid.colour,border=NA)
                  }#end for
               }#end if
               if (plotgrid){ 
                  abline(v=whenplot8$levels,h=axTicks(side=2),col=grid.colour,lty="solid")
               }#end if
               for (d in seq(from=1,to=ndbh+1,by=1)){
                  lines(datum$tomonth,thisvar[,d,p],type="l",col=dbhcols[d],lwd=lwidth)
               }#end for
               legend( x      = legwhere
                     , inset  = inset
                     , bg     = background
                     , legend = dbhnames
                     , col    = dbhcols
                     , ncol   = min(pretty.box(ndbh+1)$ncol,3)
                     , title  = expression(bold("DBH class"))
                     , lwd    = lwidth
                     )#end legend

               if (outform[o] == "x11"){
                  locator(n=1)
                  dev.off()
               }else{
                  dev.off()
               }#end if
               dummy = clean.tmp()
            } #end for outform
         }#end for (p in pftuse)
      }#end if (tseragbpft)
   } #end for tseries
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #   Plot the comparison between observations and model.                                 #
   #---------------------------------------------------------------------------------------#
   cat("    + Comparisons of time series (model vs. observations)...","\n")
   for (cc in 1:ncompmodel){

      #----- Retrieve variable information from the list. ---------------------------------#
      compnow      = compmodel[[cc]]
      vname        = compnow$vnam  
      description  = compnow$desc  
      unit         = compnow$unit  
      lcolours     = compnow$colour
      llwd         = compnow$lwd
      llwd         = compnow$lwd
      ltype        = compnow$type
      plog         = compnow$plog
      legpos       = compnow$legpos
      plotit       = compnow$emean

      #----- Check whether there are observations for this particular site. ---------------#
      if (iata == "mao" | iata == "bdf"){
         obsnow = "obs.m34"
      }else if(iata == "stm"){
         obsnow = "obs.s67"
      }else if(iata == "rao"){
         obsnow = "obs.pdg"
      }else if(iata == "jpr"){
         obsnow = "obs.fns"
      }else if(iata == "btr"){
         obsnow = "obs.s77"
      }else{
         obsnow = paste("obs.",iata,sep="")
      }#end if

      #------------------------------------------------------------------------------------#
      #     Last check to see if we should plot it or not.                                 #
      #------------------------------------------------------------------------------------#
      plotit       = plotit && obsnow %in% ls()
      if (plotit){
         thisobs = get(obsnow)
         obswhen = thisobs$tomonth
         sel     = datum$tomonth >= min(obswhen) & datum$tomonth <= max(obswhen)
         plotit  = any(sel)
      }#end if
      #------------------------------------------------------------------------------------#





      #------------------------------------------------------------------------------------#
      #     Enter here only if there is any overlap of time between observations and       #
      # model.                                                                             #
      #------------------------------------------------------------------------------------#
      if (plotit){
         #---------------------------------------------------------------------------------#
         #    Copy the observations to a scratch variable.                                 #
         #---------------------------------------------------------------------------------#
         mnvar   = paste("emean",vname,sep=".")
         obsmean = thisobs[[mnvar]]
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #    Check whether the time series directory exists.  If not, create it.          #
         #---------------------------------------------------------------------------------#
         outdir   = paste(outpref,"compemean",sep="/")
         if (! file.exists(outdir)) dir.create(outdir)
         cat("      - ",description," comparison...","\n")
         #---------------------------------------------------------------------------------#



         #----- Define the number of layers. ----------------------------------------------#
         thiswhen  = datum$tomonth [sel]
         thismean  = emean[[vname]][sel]
         #---------------------------------------------------------------------------------# 



         #----- Find the plot range. ------------------------------------------------------#
         ylimit = pretty.xylim(u=c(thismean,obsmean),fracexp=scalleg,is.log=FALSE)
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Find the nice scale for time.                                               #
         #---------------------------------------------------------------------------------#
         whenplote = pretty.time(obswhen,n=8)
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Check if the directory exists.  If not, create it.                          #
         #---------------------------------------------------------------------------------#

         #----- Loop over formats. --------------------------------------------------------#
         for (o in 1:nout){
            fichier = paste(outdir,"/",vname,".",outform[o],sep="")
            if(outform[o] == "x11"){
               X11(width=size$width,height=size$height,pointsize=ptsz)
            }else if(outform[o] == "png"){
               png(filename=fichier,width=size$width*depth,height=size$height*depth
                  ,pointsize=ptsz,res=depth)
            }else if(outform[o] == "eps"){
               postscript(file=fichier,width=size$width,height=size$height
                         ,pointsize=ptsz,paper=size$paper)
            }else if(outform[o] == "pdf"){
               pdf(file=fichier,onefile=FALSE
                  ,width=size$width,height=size$height,pointsize=ptsz,paper=size$paper)
            }#end if

            #----- Load variable ----------------------------------------------------------#
            letitre = paste(description," - ",lieu,"\n","Monthly mean",sep="")
            par(par.user)
            plot(x=thiswhen,y=thismean,type="n",main=letitre,xlab="Time"
                ,ylim=ylimit,ylab=paste("[",unit,"]",sep=""),log=plog,xaxt="n"
                ,cex.main=cex.main)
            axis(side=1,at=whenplote$levels,labels=whenplote$labels,padj=whenplote$padj)
            if (plotgrid){
               abline(v=whenplote$levels,h=axTicks(side=2),col=grid.colour,lty="solid")
            }#end if
            points(x=thiswhen,y=thismean,col=lcolours[1],lwd=llwd[1],type=ltype
                  ,pch=16,cex=1.0)
            points(x=obswhen,y=obsmean ,col=lcolours[2],lwd=llwd[2],type=ltype
                  ,pch=16,cex=1.0)
            legend(x=legpos,inset=inset,legend=c("Model","Observation")
                  ,col=lcolours,lwd=llwd,cex=1.0,pch=16)
            if (outform[o] == "x11"){
               locator(n=1)
               dev.off()
            }else{
               dev.off()
            }#end if
            dummy = clean.tmp()
         } #end for outform
      }#end if plotit
   }#end for ncompare
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #   Plot the comparison between observations and model.                                 #
   #---------------------------------------------------------------------------------------#
   cat("    + Comparisons of monthly means (model vs. observations)...","\n")
   for (cc in 1:ncompmodel){

      #----- Retrieve variable information from the list. ---------------------------------#
      compnow      = compmodel[[cc]]
      vname        = compnow$vnam  
      description  = compnow$desc  
      unit         = compnow$unit  
      plotsd       = compnow$plotsd
      lcolours     = compnow$colour
      errcolours   = compnow$errcol
      angle        = compnow$angle
      dens         = compnow$dens
      llwd         = compnow$lwd
      shwd         = compnow$shwd
      llwd         = compnow$lwd
      ltype        = compnow$type
      plog         = compnow$plog
      legpos       = compnow$legpos
      plotit       = compnow$mmean

      #----- Check whether there are observations for this particular site. ---------------#
      if (iata == "mao" | iata == "bdf"){
         obsnow = "obs.m34"
      }else if(iata == "stm"){
         obsnow = "obs.s67"
      }else if(iata == "rao"){
         obsnow = "obs.pdg"
      }else if(iata == "jpr"){
         obsnow = "obs.fns"
      }else if(iata == "btr"){
         obsnow = "obs.s77"
      }else{
         obsnow = paste("obs.",iata,sep="")
      }#end if

      plotit       = plotit && obsnow %in% ls()





      if (plotit){
         #---------------------------------------------------------------------------------#
         #    Copy the observations to a scratch variable.                                 #
         #---------------------------------------------------------------------------------#
         thisobs = get(obsnow)
         mnvar   = paste("mmean",vname,sep=".")
         sdvar   = paste("msdev",vname,sep=".")
         obsmean = thisobs[[mnvar]]
         obssdev = thisobs[[sdvar]]
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Some observations do not have enough measurements to make a full year.  If #
         # this is the case, then we must split the observations into smaller intervals so #
         # the polygon works.  In case no observation is available, make the vectors NULL  #
         # so we will not plot observations at all.                                        #
         #---------------------------------------------------------------------------------#
         if (all(is.na(obsmean+obssdev))){
            obs.x     = NULL
            obs.ylow  = NULL
            obs.yhigh = NULL
         }else{
            #------ Find the periods with continous data. ---------------------------------#
            ok        = is.finite(obsmean+obssdev)
            obs.x     = montmont[ok]
            obs.ylow  = obsmean [ok] - obssdev[ok]
            obs.yhigh = obsmean [ok] + obssdev[ok]
            #------------------------------------------------------------------------------#
         }#end if
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #    Check whether the time series directory exists.  If not, create it.          #
         #---------------------------------------------------------------------------------#
         outdir   = paste(outpref,"compmmean",sep="/")
         if (! file.exists(outdir)) dir.create(outdir)
         cat("      - ",description," comparison...","\n")
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------# 
         #    Define the number of layers.  Some variables have no standard deviation in   #
         # the model, so Make them 0 if this is the case.                                  #
         #---------------------------------------------------------------------------------# 
         thismean = mmean[[vname]]
         thissdev = msdev[[vname]]
         if (length(msdev[[vname]]) == 0){
            thissdev = 0. * thismean
         }else{
            thissdev = msdev[[vname]]
         }#end if
         mod.x     = montmont
         mod.ylow  = thismean - thissdev
         mod.yhigh = thismean + thissdev 
         #---------------------------------------------------------------------------------# 



         #----- Find the plot range. ------------------------------------------------------#
         if (plotsd){
            ylimit    = c(mod.ylow,mod.yhigh,obs.ylow,obs.yhigh)
         }else{
            ylimit    = c(thismean,obsmean)
         }#end if
         ylimit = pretty.xylim(u=ylimit,fracexp=scalleg,is.log=FALSE)
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Check if the directory exists.  If not, create it.                          #
         #---------------------------------------------------------------------------------#

         #----- Loop over formats. --------------------------------------------------------#
         for (o in 1:nout){
            fichier = paste(outdir,"/",vname,".",outform[o],sep="")
            if(outform[o] == "x11"){
               X11(width=size$width,height=size$height,pointsize=ptsz)
            }else if(outform[o] == "png"){
               png(filename=fichier,width=size$width*depth,height=size$height*depth
                  ,pointsize=ptsz,res=depth)
            }else if(outform[o] == "eps"){
               postscript(file=fichier,width=size$width,height=size$height
                         ,pointsize=ptsz,paper=size$paper)
            }else if(outform[o] == "pdf"){
               pdf(file=fichier,onefile=FALSE
                  ,width=size$width,height=size$height,pointsize=ptsz,paper=size$paper)
            }#end if

            #----- Load variable ----------------------------------------------------------#
            letitre = paste(description," - ",lieu,"\n","Monthly mean",sep="")
            par(par.user)
            plot(x=montmont,y=thismean,type="n",main=letitre,xlab="Time"
                ,ylim=ylimit,ylab=paste("[",unit,"]",sep=""),log=plog,xaxt="n"
                ,cex.main=cex.main)
            axis(side=1,at=mplot$levels,labels=mplot$labels,padj=mplot$padj)
            if (plotgrid){ 
               abline(v=mplot$levels,h=axTicks(side=2),col=grid.colour,lty="solid")
            }#end if
            if (plotsd){
               if (is.null(obs.x)){
                  mod.x.poly = c(mod.x,rev(mod.x))
                  mod.y.poly = c(mod.ylow,rev(mod.yhigh))
                  mod.keep   = is.finite(mod.y.poly)
                  err.x      = mod.x.poly[mod.keep]
                  err.y      = mod.y.poly[mod.keep]
                  polygon(x=err.x,y=err.y,col=errcolours[1],angle=angle[1],density=dens[1]
                         ,lty="solid",lwd=shwd[1])
               }else{
                  mod.x.poly = c(mod.x,rev(mod.x))
                  mod.y.poly = c(mod.ylow,rev(mod.yhigh))
                  mod.keep   = is.finite(mod.y.poly)
                  obs.x.poly = c(obs.x,rev(obs.x))
                  obs.y.poly = c(obs.ylow,rev(obs.yhigh))
                  obs.keep   = is.finite(obs.y.poly)

                  err.x = c(mod.x.poly[mod.keep],NA,obs.x.poly[obs.keep])
                  err.y = c(mod.y.poly[mod.keep],NA,obs.y.poly[obs.keep])
                  polygon(x=err.x,y=err.y,col=errcolours,angle=angle,density=dens
                         ,lty="solid",lwd=shwd)
               }#end if
            }#end if
            points(x=montmont,y=thismean,col=lcolours[1],lwd=llwd[1],type=ltype
                  ,pch=16,cex=1.0)
            points(x=montmont,y=obsmean ,col=lcolours[2],lwd=llwd[2],type=ltype
                  ,pch=16,cex=1.0)
            if (plotsd){
               legend(x=legpos,inset=inset,legend=c("Model","Observation")
                     ,fill=errcolours,angle=angle,density=dens,lwd=llwd,col=lcolours
                     ,bg=background,title="Shaded areas = 1 SD",cex=1.0,pch=16)
            }else{
               legend(x=legpos,inset=inset,legend=c("Model","Observation")
                     ,col=lcolours,lwd=llwd,cex=1.0,pch=16)
            }#end if
            if (outform[o] == "x11"){
               locator(n=1)
               dev.off()
            }else{
               dev.off()
            }#end if
            dummy = clean.tmp()
         } #end for outform
      }#end if plotit
   }#end for ncompare
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #   Plot the comparison between observations and model.                                 #
   #---------------------------------------------------------------------------------------#
   cat("      * Comparisons of mean diurnal cycle (model vs. observations)...","\n")
   for (cc in 1:ncompmodel){

      #----- Retrieve variable information from the list. ---------------------------------#
      compnow      = compmodel[[cc]]
      vname        = compnow$vnam  
      description  = compnow$desc  
      unit         = compnow$unit  
      plotsd       = compnow$plotsd
      lcolours     = compnow$colour
      errcolours   = compnow$errcol
      angle        = compnow$angle
      dens         = compnow$dens
      llwd         = compnow$lwd
      shwd         = compnow$shwd
      llwd         = compnow$lwd
      ltype        = compnow$type
      plog         = compnow$plog
      legpos       = compnow$legpos
      plotit       = compnow$qmean

      #----- Check whether there are observations for this particular site. ---------------#
      if (iata == "mao" | iata == "bdf"){
         obsnow = "obs.m34"
      }else if(iata == "stm"){
         obsnow = "obs.s67"
      }else if(iata == "rao"){
         obsnow = "obs.pdg"
      }else if(iata == "jpr"){
         obsnow = "obs.fns"
      }else if(iata == "btr"){
         obsnow = "obs.s77"
      }else{
         obsnow = paste("obs.",iata,sep="")
      }#end if
      plotit       = plotit && obsnow %in% ls()

      if (plotit){
         #---------------------------------------------------------------------------------#
         #    Copy the observations to a scratch variable.                                 #
         #---------------------------------------------------------------------------------#
         thisobs = get(obsnow)
         mnvar   = paste("qmean",vname,sep=".")
         sdvar   = paste("qsdev",vname,sep=".")
         obsmean = thisobs[[mnvar]]
         obssdev = thisobs[[sdvar]]
         #----- Append 1st hour after the last. -------------------------------------------#
         obsmean = cbind(obsmean,obsmean[,1])
         obssdev = cbind(obssdev,obssdev[,1])
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Some observations do not have enough measurements to make a full year.  If #
         # this is the case, then we must split the observations into smaller intervals so #
         # the polygon works.  In case no observation is available, make the vectors NULL  #
         # so we will not plot observations at all.                                        #
         #---------------------------------------------------------------------------------#
         obs.ylow  = obsmean - obssdev
         obs.yhigh = obsmean + obssdev
         #---------------------------------------------------------------------------------#




         #---------------------------------------------------------------------------------#
         #    Check whether the time series directory exists.  If not, create it.          #
         #---------------------------------------------------------------------------------#
         outdir   = paste(outpref,"compdcyc",sep="/")
         if (! file.exists(outdir)) dir.create(outdir)
         outtheme = paste(outdir,vname,sep="/")
         if (! file.exists(outtheme)) dir.create(outtheme)
         cat("      + ",description," comparison...","\n")
         #---------------------------------------------------------------------------------#




         #---------------------------------------------------------------------------------# 
         #    Define the number of layers.  Some variables have no standard deviation in   #
         # the model, so Make them 0 if this is the case.  We also append the last hour    #
         # before the first one so 00 UTC appears in the left.                             #
         #---------------------------------------------------------------------------------# 
         thismean  = umean[[vname]]
         thismean  = cbind(thismean[,ndcycle],thismean)
         if (length(usdev[[vname]]) == 0){
            thissdev = 0. * thismean
         }else{
            thissdev = usdev[[vname]]
            thissdev  = cbind(thissdev[,ndcycle],thissdev)
         }#end if
         mod.ylow  = thismean - thissdev
         mod.yhigh = thismean + thissdev 
         #---------------------------------------------------------------------------------# 


         #----- Find the plot range. ------------------------------------------------------#
         if (plotsd){
            ylimit    = c(mod.ylow,mod.yhigh,obs.ylow,obs.yhigh)
         }else{
            ylimit    = c(thismean,obsmean)
         }#end if
         ylimit = pretty.xylim(u=ylimit,fracexp=scalleg,is.log=FALSE)
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #      Loop over all months.                                                      #
         #---------------------------------------------------------------------------------#
         for (pmon in 1:12){
            cmon    = substring(100+pmon,2,3)
            namemon = mlist[pmon]

            #------------------------------------------------------------------------------#
            #     Check if the directory exists.  If not, create it.                       #
            #------------------------------------------------------------------------------#
            cat("        > ",description," time series - ",namemon,"...","\n")

            #----- Loop over formats. -----------------------------------------------------#
            for (o in 1:nout){
               fichier = paste(outtheme,"/",vname,"-",cmon,".",outform[o]
                              ,sep="")
               if(outform[o] == "x11"){
                  X11(width=size$width,height=size$height,pointsize=ptsz)
               }else if(outform[o] == "png"){
                  png(filename=fichier,width=size$width*depth,height=size$height*depth
                     ,pointsize=ptsz,res=depth)
               }else if(outform[o] == "eps"){
                  postscript(file=fichier,width=size$width,height=size$height
                            ,pointsize=ptsz,paper=size$paper)
               }else if(outform[o] == "pdf"){
                  pdf(file=fichier,onefile=FALSE
                     ,width=size$width,height=size$height,pointsize=ptsz,paper=size$paper)
               }#end if

               #----- Load variable -------------------------------------------------------#
               letitre = paste(description," - ",lieu,"\n"
                              ,"Mean diurnal cycle - ",namemon,sep="")
               par(par.user)
               plot(x=thisday,y=thismean[pmon,],type="n",main=letitre,xlab="Time"
                   ,ylim=ylimit,ylab=paste("[",unit,"]",sep=""),log=plog,xaxt="n"
                   ,cex.main=cex.main)
               axis(side=1,at=uplot$levels,labels=uplot$labels,padj=uplot$padj)
               if (plotgrid){ 
                  abline(v=uplot$levels,h=axTicks(side=2),col=grid.colour,lty="solid")
               }#end if
               if (plotsd){
                  mod.x.now     = thisday
                  mod.ylow.now  = mod.ylow [pmon,]
                  mod.yhigh.now = mod.yhigh[pmon,]
                  #------ Find the periods with continous data. ---------------------------#
                  ok        = is.finite(obs.ylow[pmon,]) & is.finite(obs.yhigh[pmon,])
                  if (any(ok)){
                     obs.x.now     = thisday       [ok]
                     obs.ylow.now  = obs.ylow [pmon,ok]
                     obs.yhigh.now = obs.yhigh[pmon,ok]
                  }else{
                     obs.x.now     = NULL
                     obs.ylow.now  = NULL
                     obs.yhigh.now = NULL
                  }#end if
                  #------------------------------------------------------------------------#

                  if (is.null(obs.x.now)){
                     mod.x.poly = c(mod.x.now,rev(mod.x.now))
                     mod.y.poly = c(mod.ylow.now,rev(mod.yhigh.now))
                     mod.keep   = is.finite(mod.y.poly)
                     err.x      = mod.x.poly[mod.keep]
                     err.y      = mod.y.poly[mod.keep]
                     polygon(x=err.x,y=err.y,col=errcolours[1],angle=angle[1]
                            ,density=dens[1],lty="solid",lwd=shwd[1])
                  }else{
                     mod.x.poly = c(mod.x.now,rev(mod.x.now))
                     mod.y.poly = c(mod.ylow.now,rev(mod.yhigh.now))
                     mod.keep   = is.finite(mod.y.poly)
                     obs.x.poly = c(obs.x.now,rev(obs.x.now))
                     obs.y.poly = c(obs.ylow.now,rev(obs.yhigh.now))
                     obs.keep   = is.finite(obs.y.poly)
                     err.x = c(mod.x.poly[mod.keep],NA,obs.x.poly[obs.keep])
                     err.y = c(mod.y.poly[mod.keep],NA,obs.y.poly[obs.keep])
                     polygon(x=err.x,y=err.y,col=errcolours,angle=angle,density=dens
                            ,lty="solid",lwd=shwd)
                  }#end if
               }#end if
               points(x=thisday,y=thismean[pmon,],col=lcolours[1]
                     ,lwd=llwd[1],type=ltype,pch=16,cex=1.0)
               points(x=thisday,y=obsmean[pmon,],col=lcolours[2]
                     ,lwd=llwd[2],type=ltype,pch=16,cex=1.0)
               if (plotsd){
                  legend(x=legpos,inset=inset,legend=c("Model","Observation")
                        ,fill=errcolours,angle=angle,density=dens,lwd=llwd,col=lcolours
                        ,bg=background,title="Shaded areas = 1 SD",cex=1.0,pch=16)
               }else{
                  legend(x=legpos,inset=inset,legend=c("Model","Observation")
                        ,col=lcolours,lwd=llwd,cex=1.0,pch=16)
               }#end if
               if (outform[o] == "x11"){
                  locator(n=1)
                  dev.off()
               }else{
                  dev.off()
               }#end if
               dummy = clean.tmp()
            } #end for outform
         }#end for pmon
      }#end if plotit
   }#end for ncompare
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #      Time series by LU.                                                               #
   #---------------------------------------------------------------------------------------#
   for (v in 1:ntslu){
      thistslu    = tslu[[v]]
      vnam        = thistslu$vnam
      description = thistslu$desc
      unit        = thistslu$unit
      plog        = thistslu$plog
      plotit      = thistslu$plt

      #----- Check whether the user wants to have this variable plotted. ------------------#
      if (plotit && any(sellu)){

         #---------------------------------------------------------------------------------#
         #    Check whether the time series directory exists.  If not, create it.          #
         #---------------------------------------------------------------------------------#
         outdir = paste(outpref,"tslu",sep="/")
         if (! file.exists(outdir)) dir.create(outdir)
         cat("      + ",description," time series for all LUs...","\n")



         #----- Load variable -------------------------------------------------------------#
         thisvar = lu[[vnam]]
         if (plog){
            #----- Eliminate non-positive values in case it is a log plot. ----------------#
            thisvar[thisvar <= 0] = NA
         }#end if
         #---------------------------------------------------------------------------------#

         #----- Loop over output formats. -------------------------------------------------#
         for (o in 1:nout){
            fichier = paste(outdir,"/",vnam,"-",suffix,".",outform[o],sep="")
            if(outform[o] == "x11"){
               X11(width=size$width,height=size$height,pointsize=ptsz)
            }else if(outform[o] == "png"){
               png(filename=fichier,width=size$width*depth,height=size$height*depth
                  ,pointsize=ptsz,res=depth)
            }else if(outform[o] == "eps"){
               postscript(file=fichier,width=size$width,height=size$height
                         ,pointsize=ptsz,paper=size$paper)
            }else if(outform[o] == "pdf"){
               pdf(file=fichier,onefile=FALSE
                  ,width=size$width,height=size$height,pointsize=ptsz,paper=size$paper)
            }#end if


            #------------------------------------------------------------------------------#
            #     Find the limit, make some room for the legend, and in case the field is  #
            # a constant, nudge the limits so the plot command will not complain.          #
            #------------------------------------------------------------------------------#
            ylimit = pretty.xylim(thisvar[,sellu],fracexp=scalleg,is.log=plog)
            if (plog){
               xylog    = "y"
               ydrought = c( exp(ylimit[1] * sqrt(ylimit[1]/ylimit[2]))
                           , exp(ylimit[2] * sqrt(ylimit[2]/ylimit[1]))
                           )#end c
            }else{
               xylog    = ""
               ydrought = c(ylimit[1] - 0.5 * diff(ylimit), ylimit[2] + 0.5 * diff(ylimit))
            }#end if
            #------------------------------------------------------------------------------#

            letitre = paste(description,lieu,sep=" - ")
            cols    = lucols[sellu]
            legs    = lunames[sellu]
            par(par.user)
            plot(datum$tomonth,thisvar[,1],type="n",main=letitre,ylim=ylimit
                ,xlab="Time",ylab=unit,xaxt="n",cex.main=0.7)
            axis(side=1,at=whenplot8$levels,labels=whenplot8$labels,padj=whenplot8$padj)

            if (drought.mark){
               for (n in 1:ndrought){
                  rect(xleft  = drought[[n]][1],ybottom = ydrought[1]
                      ,xright = drought[[n]][2],ytop    = ydrought[2]
                      ,col    = grid.colour,border=NA)
               }#end for
            }#end if
            if (plotgrid){ 
               abline(v=whenplot8$levels,h=axTicks(side=2),col=grid.colour,lty="solid")
            }#end if
            for (n in 1:(nlu+1)){
               if (sellu[n]){
                  lines(datum$tomonth,thisvar[,n],type="l",col=lucols[n],lwd=lwidth)
               }#end if
            }#end for
            legend( x      = legwhere
                  , inset  = inset
                  , bg     = background
                  , legend = legs
                  , col    = cols
                  , lwd    = lwidth
                  , ncol   = min(3,pretty.box(n.sellu)$ncol)
                  , title  = expression(bold("Land use type"))
                  )#end legend

            if (outform[o] == "x11"){
               locator(n=1)
               dev.off()
            }else{
               dev.off()
            }#end if
            dummy = clean.tmp()
         } #end for outform
      }#end if (tseragbpft)
   } #end for tseries
   #---------------------------------------------------------------------------------------#





   #---------------------------------------------------------------------------------------#
   #   Plot disturbance rate by disturbance transition.                                    #
   #---------------------------------------------------------------------------------------#
   if (tserdist && any(seldist)){
      cat("      + Disturbance rate time series for all disturbances...","\n")
      for (o in 1:nout){
         fichier = paste(outpref,"/disturb-",suffix,".",outform[o],sep="")
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
               ,width=size$width,height=size$height,pointsize=ptsz,paper=size$paper)
         }#end if

         #---------------------------------------------------------------------------------#
         #     Find the limit, make some room for the legend, and in case the field is a   #
         #  constant, nudge the limits so the plot command will not complain.              #
         #---------------------------------------------------------------------------------#
         ylimit  = NULL
         for (jlu in 1:nlu){
            for (ilu in 1:nlu){
               if (seldist[ilu,jlu]){
                  ylimit = c(ylimit,lu$dist[,ilu,jlu])
               }#end if
            }#end for
         }#end for
         ylimit   = pretty.xylim(u=ylimit,fracexp=scalleg,is.log=FALSE)
         ydrought = c(ylimit[1] - 0.5 * diff(ylimit), ylimit[2] + 0.5 * diff(ylimit))
         #---------------------------------------------------------------------------------#

         letitre = paste("Disturbance rates",lieu,sep=" - ")
         cols    = NULL
         legs    = NULL
         par(par.user)
         plot(datum$tomonth,lu$dist[,1,1],type="n",main=letitre,ylim=ylimit
             ,xlab="Time",ylab="[1/yr]",xaxt="n",cex.main=0.7)
            axis(side=1,at=whenplot8$levels,labels=whenplot8$labels,padj=whenplot8$padj)
            if (drought.mark){
               for (n in 1:ndrought){
                  rect(xleft  = drought[[n]][1],ybottom = ydrought[1]
                      ,xright = drought[[n]][2],ytop    = ydrought[2]
                      ,col    = grid.colour,border=NA)
               }#end for
            }#end if
            if (plotgrid){ 
               abline(v=whenplot8$levels,h=axTicks(side=2),col=grid.colour,lty="solid")
            }#end if
         n = 0
         for (jlu in 1:nlu){
            for (ilu in 1:nlu){
               n = n + 1
               if (seldist[ilu,jlu]){
                  cols = c(cols,distcols[n])
                  legs = c(legs,distnames[n])
                  lines(datum$tomonth,lu$dist[,ilu,jlu],type="l"
                       ,col=distcols[n],lwd=lwidth)
               }#end if
            }#end for
         }#end for
         legend(x      = legwhere
               ,inset  = inset
               ,bg     = background
               ,legend = legs
               ,col    = cols
               ,lwd    = lwidth
               , ncol  = min(3,pretty.box(n)$ncol)
               , title = expression(bold("Transition"))
               )#end legend

         if (outform[o] == "x11"){
            locator(n=1)
            dev.off()
         }else{
            dev.off()
         }#end if
         dummy = clean.tmp()
      } #end for outform
   }#end if
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #   Plot the time series diagrams showing months and years.                             #
   #---------------------------------------------------------------------------------------#
   cat("      * Plot some time series with groups of variables...","\n")
   for (hh in 1:ntheme){

      #----- Retrieve variable information from the list. ---------------------------------#
      themenow     = theme[[hh]]
      vnames       = themenow$vnam  
      description  = themenow$desc  
      lcolours     = themenow$colour
      llwd         = themenow$lwd
      if (emean.line){
         ltype        = "l"
      }else{
         ltype        = themenow$type 
      }#end if
      plog         = themenow$plog
      prefix       = themenow$prefix
      group        = themenow$title
      unit         = themenow$unit  
      legpos       = themenow$legpos
      plotit       = themenow$emean
   
      if (plotit){

         #---------------------------------------------------------------------------------#
         #    Check whether the time series directory exists.  If not, create it.          #
         #---------------------------------------------------------------------------------#
         outdir = paste(outpref,"theme_emean",sep="/")
         if (! file.exists(outdir)) dir.create(outdir)
         cat("      + ",group," time series for several variables...","\n")


         #----- Define the number of layers. ----------------------------------------------#
         nlayers   = length(vnames)
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Find the limit, make some room for the legend, and in case the field is a   #
         # constant, nudge the limits so the plot command will not complain.               #
         #---------------------------------------------------------------------------------#
         ylimit    = NULL
         for (l in 1:nlayers) ylimit  = c(ylimit,emean[[vnames[l]]])
         ylimit = pretty.xylim(u=ylimit,fracexp=scalleg,is.log=plog)
         if (plog){
            xylog    = "y"
            ydrought = c( exp(ylimit[1] * sqrt(ylimit[1]/ylimit[2]))
                        , exp(ylimit[2] * sqrt(ylimit[2]/ylimit[1]))
                        )#end c
         }else{
            xylog    = ""
            ydrought = c(ylimit[1] - 0.5 * diff(ylimit), ylimit[2] + 0.5 * diff(ylimit))
         }#end if
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Check if the directory exists.  If not, create it.                          #
         #---------------------------------------------------------------------------------#

         #----- Loop over formats. --------------------------------------------------------#
         for (o in 1:nout){
            fichier = paste(outdir,"/",prefix,"-",suffix,".",outform[o],sep="")
            if(outform[o] == "x11"){
               X11(width=size$width,height=size$height,pointsize=ptsz)
            }else if(outform[o] == "png"){
               png(filename=fichier,width=size$width*depth,height=size$height*depth
                  ,pointsize=ptsz,res=depth)
            }else if(outform[o] == "eps"){
               postscript(file=fichier,width=size$width,height=size$height
                         ,pointsize=ptsz,paper=size$paper)
            }else if(outform[o] == "pdf"){
               pdf(file=fichier,onefile=FALSE
                  ,width=size$width,height=size$height,pointsize=ptsz,paper=size$paper)
            }#end if

            #----- Load variable ----------------------------------------------------------#
            thisvar = emean[[vnames[1]]]

            letitre = paste(" Time series: ",group," \n",lieu,sep="")

            par(par.user)
            plot(x=datum$tomonth,y=thisvar,type="n",main=letitre,xlab="Time"
                ,ylim=ylimit,ylab=paste("[",unit,"]",sep=""),log=xylog,xaxt="n"
                ,cex.main=cex.main)
            axis(side=1,at=whenplot8$levels,labels=whenplot8$labels,padj=whenplot8$padj)
            if (drought.mark){
               for (n in 1:ndrought){
                  rect(xleft  = drought[[n]][1],ybottom = ydrought[1]
                      ,xright = drought[[n]][2],ytop    = ydrought[2]
                      ,col    = grid.colour,border=NA)
               }#end for
            }#end if
            if (plotgrid){ 
               abline(v=whenplot8$levels,h=axTicks(side=2),col=grid.colour,lty="solid")
            }#end if
            for (l in 1:nlayers){
               thisvar = emean[[vnames[l]]]
               points(x=datum$tomonth,y=thisvar,col=lcolours[l]
                     ,lwd=llwd[l],type=ltype,pch=16)
            }#end for
            legend( x      = legpos
                  , inset  = inset
                  , legend = description
                  , col    = lcolours
                  , lwd    = llwd
                  , ncol   = min(3,pretty.box(nlayers)$ncol)
                  )#end legend
            if (outform[o] == "x11"){
               locator(n=1)
               dev.off()
            }else{
               dev.off()
            }#end if
            dummy = clean.tmp()
         } #end for outform
      }#end if plotit
   }#end for ntser
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #   Plot the time series diagrams showing months and years.                             #
   #---------------------------------------------------------------------------------------#
   cat("      * Plot some monthly means of groups of variables ...","\n")
   for (hh in 1:ntheme){

      #----- Retrieve variable information from the list. ---------------------------------#
      themenow     = theme[[hh]]
      vnames       = themenow$vnam  
      description  = themenow$desc  
      lcolours     = themenow$colour
      llwd         = themenow$lwd
      ltype        = themenow$type
      plog         = themenow$plog
      prefix       = themenow$prefix
      group        = themenow$title
      unit         = themenow$unit  
      legpos       = themenow$legpos
      plotit       = themenow$mmean
   
      if (plotit){

         #---------------------------------------------------------------------------------#
         #    Check whether the time series directory exists.  If not, create it.          #
         #---------------------------------------------------------------------------------#
         outdir = paste(outpref,"theme_mmean",sep="/")
         if (! file.exists(outdir)) dir.create(outdir)
         cat("      + ",group," time series for several variables...","\n")


         #----- Define the number of layers. ----------------------------------------------#
         nlayers   = length(vnames)
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Find the limit, make some room for the legend, and in case the field is a   #
         # constant, nudge the limits so the plot command will not complain.               #
         #---------------------------------------------------------------------------------#
         ylimit    = NULL
         for (l in 1:nlayers) ylimit  = c(ylimit,mmean[[vnames[l]]])
         ylimit = pretty.xylim(u=ylimit,fracexp=scalleg,is.log=plog)
         if (plog){
            xylog    = "y"
            ydrought = c( exp(ylimit[1] * sqrt(ylimit[1]/ylimit[2]))
                        , exp(ylimit[2] * sqrt(ylimit[2]/ylimit[1]))
                        )#end c
         }else{
            xylog    = ""
            ydrought = c(ylimit[1] - 0.5 * diff(ylimit), ylimit[2] + 0.5 * diff(ylimit))
         }#end if
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Check if the directory exists.  If not, create it.                          #
         #---------------------------------------------------------------------------------#

         #----- Loop over formats. --------------------------------------------------------#
         for (o in 1:nout){
            fichier = paste(outdir,"/",prefix,"-",suffix,".",outform[o],sep="")
            if(outform[o] == "x11"){
               X11(width=size$width,height=size$height,pointsize=ptsz)
            }else if(outform[o] == "png"){
               png(filename=fichier,width=size$width*depth,height=size$height*depth
                  ,pointsize=ptsz,res=depth)
            }else if(outform[o] == "eps"){
               postscript(file=fichier,width=size$width,height=size$height
                         ,pointsize=ptsz,paper=size$paper)
            }else if(outform[o] == "pdf"){
               pdf(file=fichier,onefile=FALSE
                  ,width=size$width,height=size$height,pointsize=ptsz,paper=size$paper)
            }#end if

            #----- Load variable ----------------------------------------------------------#
            thisvar = mmean[[vnames[1]]]

            letitre = paste(" Time series: ",group," \n",lieu,sep="")

            par(par.user)
            plot(x=montmont,y=thisvar,type="n",main=letitre,xlab="Month"
                ,ylim=ylimit,ylab=paste("[",unit,"]",sep=""),log=xylog,xaxt="n"
                ,cex.main=cex.main)
            axis(side=1,at=mplot$levels,labels=mplot$labels,padj=mplot$padj)
            if (drought.mark){
               for (n in 1:ndrought){
                  rect(xleft  = drought[[n]][1],ybottom = ydrought[1]
                      ,xright = drought[[n]][2],ytop    = ydrought[2]
                      ,col    = grid.colour,border=NA)
               }#end for
            }#end if
            if (plotgrid){ 
               abline(v=mplot$levels,h=axTicks(side=2),col=grid.colour,lty="solid")
            }#end if
            for (l in 1:nlayers){
               thisvar = mmean[[vnames[l]]]
               points(x=montmont,y=thisvar,col=lcolours[l],lwd=llwd[l],type=ltype,pch=16)
            }#end for
            legend( x      = legpos
                  , inset  = inset
                  , legend = description
                  , col    = lcolours
                  , lwd    = llwd
                  , ncol   = min(3,pretty.box(nlayers)$ncol)
                  )#end legend
            if (outform[o] == "x11"){
               locator(n=1)
               dev.off()
            }else{
               dev.off()
            }#end if
            dummy = clean.tmp()
         } #end for outform
      }#end if plotit
   }#end for ntser
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #   Plot the climatology of the mean diurnal cycle.                                     #
   #---------------------------------------------------------------------------------------#
   cat("      * Plot the mean diel of groups of variables...","\n")
   for (hh in 1:ntheme){

      #----- Retrieve variable information from the list. ---------------------------------#
      themenow     = theme[[hh]]
      vnames       = themenow$vnam  
      description  = themenow$desc  
      lcolours     = themenow$colour
      llwd         = themenow$lwd
      ltype        = themenow$type
      plog         = themenow$plog
      prefix       = themenow$prefix
      group        = themenow$title
      unit         = themenow$unit  
      legpos       = themenow$legpos
      plotit       = themenow$qmean 
      if (plog){ 
         xylog = "y"
      }else{
         xylog = ""
      }#end if


      if (plotit){

         #---------------------------------------------------------------------------------#
         #    Check whether the time series directory exists.  If not, create it.          #
         #---------------------------------------------------------------------------------#
         outdir   = paste(outpref,"theme_qmean",sep="/")
         if (! file.exists(outdir)) dir.create(outdir)
         outtheme = paste(outdir,prefix,sep="/")
         if (! file.exists(outtheme)) dir.create(outtheme)
         cat("      + ",group," diurnal cycle for several variables...","\n")


         #----- Define the number of layers. ----------------------------------------------#
         nlayers   = length(vnames)
         ylimit    = NULL
         for (l in 1:nlayers) ylimit = c(ylimit,umean[[vnames[l]]])
         ylimit = pretty.xylim(u=ylimit,fracexp=scalleg,is.log=FALSE)
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #      Loop over all months.                                                      #
         #---------------------------------------------------------------------------------#
         for (pmon in 1:12){
            cmon    = substring(100+pmon,2,3)
            namemon = mlist[pmon]

            #------------------------------------------------------------------------------#
            #     Check if the directory exists.  If not, create it.                       #
            #------------------------------------------------------------------------------#

            #----- Loop over formats. -----------------------------------------------------#
            for (o in 1:nout){
               fichier = paste(outtheme,"/",prefix,"-",cmon,"-",suffix,".",outform[o]
                              ,sep="")
               if(outform[o] == "x11"){
                  X11(width=size$width,height=size$height,pointsize=ptsz)
               }else if(outform[o] == "png"){
                  png(filename=fichier,width=size$width*depth,height=size$height*depth
                     ,pointsize=ptsz,res=depth)
               }else if(outform[o] == "eps"){
                  postscript(file=fichier,width=size$width,height=size$height
                            ,pointsize=ptsz,paper=size$paper)
               }else if(outform[o] == "pdf"){
                  pdf(file=fichier,onefile=FALSE
                     ,width=size$width,height=size$height,pointsize=ptsz,paper=size$paper)
               }#end if

               #----- Load variable -------------------------------------------------------#
               thisvar = umean[[vnames[1]]]
               thisvar = cbind(thisvar[,ndcycle],thisvar)

               letitre = paste(group," - ",lieu,"\n"
                              ,"Mean diurnal cycle - ",namemon,sep="")

               par(par.user)
               plot(x=thisday,y=thisvar[pmon,],type="n",main=letitre,xlab="Time"
                   ,ylim=ylimit,ylab=paste("[",unit,"]",sep=""),log=xylog,xaxt="n"
                   ,cex.main=cex.main)
               axis(side=1,at=uplot$levels,labels=uplot$labels,padj=uplot$padj)
               if (plotgrid){ 
                  abline(v=uplot$levels,h=axTicks(side=2),col=grid.colour,lty="solid")
               }#end if
               for (l in 1:nlayers){
                  thisvar = umean[[vnames[l]]]
                  thisvar = cbind(thisvar[,ndcycle],thisvar)
                  points(x=thisday,y=thisvar[pmon,],col=lcolours[l]
                        ,lwd=llwd[l],type=ltype,pch=16)
               }#end for
               legend( x      = legpos
                     , inset  = inset
                     , legend = description
                     , col    = lcolours
                     , lwd    = llwd
                     , ncol   = min(3,pretty.box(nlayers)$ncol)
                     )#end legend
               if (outform[o] == "x11"){
                  locator(n=1)
                  dev.off()
               }else{
                  dev.off()
               }#end if
               dummy = clean.tmp()
            } #end for outform
         }#end for pmon
      }#end if plotit
   }#end for ntser
   #---------------------------------------------------------------------------------------#





   #---------------------------------------------------------------------------------------#
   #   Plot the climatology of the soil properties.                                        #
   #---------------------------------------------------------------------------------------#
   for (v in 1:nsoilplot){

      #----- Retrieve variable information from the list. ---------------------------------#
      thissoil    = soilplot[[v]]
      vnam        = thissoil$vnam
      description = thissoil$desc
      unit        = thissoil$unit
      vcscheme    = thissoil$csch
      pnlog       = thissoil$pnlog
      plotit      = thissoil$mmean

      if (plotit){

         #---------------------------------------------------------------------------------#
         #     Check if the directory exists.  If not, create it.                          #
         #---------------------------------------------------------------------------------#
         outdir  =  paste(outpref,"soil_mmean",sep="/")
         if (! file.exists(outdir)) dir.create(outdir)
         cat("      + Climatology profile of ",description,"...","\n")

         #----- Find the number of rows and columns, and the axes. ------------------------#
         monaxis  = sort(unique(datum$month))
         soilaxis = slz
         nmon     = length(monaxis)
         nsoil    = nzg

         #----- Save the meaningful months and years. -------------------------------------#
         monat   = 1:12
         monlab  = c("J","F","M","A","M","J","J","A","S","O","N","D")

         #----- Convert the vector data into an array. ------------------------------------#
         vararr  = mmean[[vnam]]

         #----- Copy Decembers ans Januaries to make the edges buffered. ------------------#
         january  = vararr[1,]
         january  = c(january,january[nzg],january[nzg])

         december = vararr[12,]
         december = c(december[1],december[1],december)

         #----- Bind first and last year to the array, to make the edges buffered. ---------#
         varbuff  = cbind(vararr[,1],vararr,vararr[,nzg])
         varbuff  = rbind(december,varbuff,january)

         #----------------------------------------------------------------------------------#
         #   Expand the month and year axes.  Make the -------------------------------------------#
         monaxis  = c(min(monaxis)-1,monaxis,max(monaxis)+1)
         soilaxis = -log(-1.0 * c( slz[1]*(slz[1]/slz[2])
                                 , soilaxis
                                 , slz[nzg]*(slz[nzg]/slz[nzg-1]) ))

         if (pnlog){
            vrange  = range(varbuff,na.rm=TRUE)
            vlevels = pretty.log(x=vrange,n=ncolsfc)
            vnlev   = length(vlevels)
         }else{
            vrange  = range(varbuff,na.rm=TRUE)
            vlevels = pretty(x=vrange,n=ncolsfc)
            vnlev   = length(vlevels)
         }#end if

         #----- Loop over formats. --------------------------------------------------------#
         for (o in 1:nout){
            fichier = paste(outdir,"/",vnam,"-",suffix,".",outform[o],sep="")
            if(outform[o] == "x11"){
               X11(width=size$width,height=size$height,pointsize=ptsz)
            }else if(outform[o] == "png"){
               png(filename=fichier,width=size$width*depth,height=size$height*depth
                  ,pointsize=ptsz,res=depth)
            }else if(outform[o] == "eps"){
               postscript(file=fichier,width=size$width,height=size$height
                         ,pointsize=ptsz,paper=size$paper)
            }else if(outform[o] == "pdf"){
               pdf(file=fichier,onefile=FALSE
                  ,width=size$width,height=size$height,pointsize=ptsz,paper=size$paper)
            }#end if

            letitre = paste(description," - ",lieu,sep="")
            par(par.user)
            sombreado(x=monaxis,y=soilaxis,z=varbuff,levels=vlevels,nlevels=vnlev
                     ,colour.palette=get(vcscheme)
                     ,plot.title=title(main=letitre,xlab="Month",ylab="Soil depth [m]"
                                      ,cex.main=0.7)
                     ,key.title=title(main=unit,cex.main=0.8)
                     ,key.log=pnlog
                     ,plot.axes={axis(side=1,at=monat,labels=monlab)
                                 axis(side=2,at=zat,labels=znice)
                                 if (fcgrid){
                                    abline(h=zat,v=monat,col=grid.colour,lty="dotted")
                                 }#end if fcgrid
                                }#end plot.axes
                     )

            if (outform[o] == "x11"){
               locator(n=1)
               dev.off()
            }else{
               dev.off()
            }#end if
            dummy = clean.tmp()
         } #end for outform
      }#end if plotit
   }#end for (v in 1:nsoilplot)
   #---------------------------------------------------------------------------------------#





   #---------------------------------------------------------------------------------------#
   #   Plot the climatology of the soil properties.                                        #
   #---------------------------------------------------------------------------------------#
   for (sts in 1:nsoilplot){

      #----- Retrieve variable information from the list. ---------------------------------#
      thissoil    = soilplot[[sts]]
      vnam        = thissoil$vnam
      description = thissoil$desc
      unit        = thissoil$unit
      vcscheme    = thissoil$csch
      pnlog       = thissoil$pnlog
      plotit      = thissoil$emean

      if (plotit){

         #---------------------------------------------------------------------------------#
         #     Check if the directory exists.  If not, create it.                          #
         #---------------------------------------------------------------------------------#
         outdir  =  paste(outpref,"soil_emean",sep="/")
         if (! file.exists(outdir)) dir.create(outdir)
         cat("      + Time series profile of ",description,"...","\n")

         #----- Find the number of rows and columns, and the axes. ------------------------#
         timeaxis  = datum$tomonth
         soilaxis  = slz
         nmon      = length(timeaxis)
         nsoil     = nzg

         #----- Convert the vector data into an array. ------------------------------------#
         vararr  = emean[[vnam]]

         #----- Copy Decembers ans Januaries to make the edges buffered. ------------------#
         first    = vararr[1,]
         first    = c(first,first[nzg],first[nzg])

         last     = vararr[ntimes,]
         last     = c(last[1],last[1],last)

         #----- Bind first and last year to the array, to make the edges buffered. --------#
         varbuff  = cbind(vararr[,1],vararr,vararr[,nzg])
         varbuff  = rbind(first,varbuff,last)

         #---------------------------------------------------------------------------------#
         #      Expand the month and year axes.  Make the first and last time equal time   #
         # steps.                                                                          #
         #---------------------------------------------------------------------------------#
         dwhen    = as.numeric(datum$tomonth[2]-datum$tomonth[1])
         whenaxis = c(chron(as.numeric(datum$tomonth[1]-dwhen))
                     ,timeaxis
                     ,chron(as.numeric(datum$tomonth[ntimes]+dwhen)))
         soilaxis = -log(-1.0 * c( slz[1]*(slz[1]/slz[2])
                                 , soilaxis
                                 , slz[nzg]*(slz[nzg]/slz[nzg-1]) ))

         if (pnlog){
            vrange  = range(varbuff,na.rm=TRUE)
            vlevels = pretty.log(x=vrange,n=ncolsfc)
            vnlev   = length(vlevels)
         }else{
            vrange  = range(varbuff,na.rm=TRUE)
            vlevels = pretty(x=vrange,n=ncolsfc)
            vnlev   = length(vlevels)
         }#end if

         #----- Loop over formats. --------------------------------------------------------#
         for (o in 1:nout){
            fichier = paste(outdir,"/",vnam,"-",suffix,".",outform[o],sep="")
            if(outform[o] == "x11"){
               X11(width=size$width,height=size$height,pointsize=ptsz)
            }else if(outform[o] == "png"){
               png(filename=fichier,width=size$width*depth,height=size$height*depth
                  ,pointsize=ptsz,res=depth)
            }else if(outform[o] == "eps"){
               postscript(file=fichier,width=size$width,height=size$height
                         ,pointsize=ptsz,paper=size$paper)
            }else if(outform[o] == "pdf"){
               pdf(file=fichier,onefile=FALSE
                  ,width=size$width,height=size$height,pointsize=ptsz,paper=size$paper)
            }#end if

            letitre = paste(description," - ",lieu,sep="")
            par(par.user)
            sombreado(x=whenaxis,y=soilaxis,z=varbuff,levels=vlevels,nlevels=vnlev
                     ,colour.palette=get(vcscheme)
                     ,plot.title=title(main=letitre,xlab="Month",ylab="Soil depth [m]"
                                      ,cex.main=0.7)
                     ,key.title=title(main=unit,cex.main=0.8)
                     ,key.log=pnlog
                     ,plot.axes={axis(side=1,at=whenplot6$levels
                                     ,labels=whenplot6$labels,padj=whenplot6$padj)
                                 axis(side=2,at=zat,labels=znice)
                                 if (fcgrid){
                                    abline(h=zat,v=whenplot6$levels,col=grid.colour
                                          ,lty="dotted")
                                 }#end if fcgrid
                                }#end plot.axes
                     )#end sombreado

            if (outform[o] == "x11"){
               locator(n=1)
               dev.off()
            }else{
               dev.off()
            }#end if
            dummy = clean.tmp()
         } #end for outform
      }#end if plotit
   }#end for nhov
   #---------------------------------------------------------------------------------------#





   #---------------------------------------------------------------------------------------#
   #   Plot a filled contour plot showing months and years.                                #
   #---------------------------------------------------------------------------------------#
   for (v in 1:nsqueeze){

      #----- Retrieve variable information from the list. ---------------------------------#
      thisfillc   = squeeze[[v]]
      vnam        = thisfillc$vnam
      description = thisfillc$desc
      unit        = thisfillc$unit
      vcscheme    = thisfillc$col.scheme
      plotit      = thisfillc$fco.mmean

      #------------------------------------------------------------------------------------#
      #     Find the first and the last full years.  These will be the actual first and    #
      # last year only if the years are complete, otherwise the first and the last year    #
      # will be taken out.                                                                 #
      #------------------------------------------------------------------------------------#
      if (monthbeg == 1){
         yearaa = yeara
      }else{
         yearaa = yeara + 1
      }# end if
      if (meszz == 12){
         yearzz = yearz
      }else{
         yearzz = yearz - 1
      }#end if
      sel      = datum$year >= yearaa & datum$year <= yearzz
      twoyears = sum(sel) >= 24

      if (plotit && twoyears){

         #---------------------------------------------------------------------------------#
         #     Check if the directory exists.  If not, create it.                          #
         #---------------------------------------------------------------------------------#
         outdir  =  paste(outpref,"fillc_mmean",sep="/")
         if (! file.exists(outdir)) dir.create(outdir)
         cat("      + ",description," time series in filled contour...","\n")

         #----- Load this variable into "thisvar". ----------------------------------------#
         thisvar = emean[[vnam]]

         #----- Find the number of rows and columns, and the axes. ------------------------#
         monaxis = sort(unique(datum$month[sel]))
         yraxis  = sort(unique(datum$year[sel]))
         nmon    = length(monaxis)
         nyear   = length(yraxis)

         #----- Save the meaningful months and years. -------------------------------------#
         monat   = 1:12
         monlab  = c("J","F","M","A","M","J","J","A","S","O","N","D")
         yrat    = pretty(yraxis)

         #----- Convert the vector data into an array. ------------------------------------#
         vararr  = array(thisvar[sel],c(nmon,nyear))

         #----- Copy Decembers ans Januaries to make the edges buffered. ------------------#
         january  = vararr[1,]
         january  = c(january,january[nyear],january[nyear])

         december = vararr[12,]
         december = c(december[1],december[1],december)

         #----- Bind first and last year to the array, to make the edges buffered. --------#
         varbuff  = cbind(vararr[,1],vararr,vararr[,nyear])
         varbuff  = rbind(december,varbuff,january)

         #----- Expand the month and year axes. -------------------------------------------#
         monaxis = c(min(monaxis)-1,monaxis,max(monaxis)+1)
         yraxis  = c(min(yraxis)-1,yraxis,max(yraxis)+1)

         vrange  = range(varbuff,na.rm=TRUE)
         vlevels = pretty(x=vrange,n=ncolsfc)
         vnlev   = length(vlevels)

         #----- Loop over formats. --------------------------------------------------------#
         for (o in 1:nout){
            fichier = paste(outdir,"/",vnam,"-",suffix,".",outform[o],sep="")
            if(outform[o] == "x11"){
               X11(width=size$width,height=size$height,pointsize=ptsz)
            }else if(outform[o] == "png"){
               png(filename=fichier,width=size$width*depth,height=size$height*depth
                  ,pointsize=ptsz,res=depth)
            }else if(outform[o] == "eps"){
               postscript(file=fichier,width=size$width,height=size$height
                         ,pointsize=ptsz,paper=size$paper)
            }else if(outform[o] == "pdf"){
               pdf(file=fichier,onefile=FALSE
                  ,width=size$width,height=size$height,pointsize=ptsz,paper=size$paper)
            }#end if

            letitre = paste(description," - ",lieu,sep="")
            par(par.user)
            sombreado(x=monaxis,y=yraxis,z=varbuff,levels=vlevels,nlevels=vnlev
                     ,colour.palette=get(vcscheme)
                     ,plot.title=title(main=letitre,xlab="Month",ylab="Year",cex.main=0.7)
                     ,key.title=title(main=unit,cex.main=0.8)
                     ,plot.axes={axis(side=1,at=monat,labels=monlab)
                                 axis(side=2,at=yrat)
                                 if (fcgrid){
                                    for (yl in yrat){
                                       abline(h=yl,col=grid.colour,lty="dotted")
                                    } #end for yl
                                    for (ml in monat){
                                       abline(v=ml,col=grid.colour,lty="dotted")
                                    } #end for ml
                                 }#end if fcgrid
                                }#end plot.axes
                     )

            if (outform[o] == "x11"){
               locator(n=1)
               dev.off()
            }else{
               dev.off()
            }#end if
            dummy = clean.tmp()
         } #end for outform
      }#end if plotit
   }#end for nhov
   #---------------------------------------------------------------------------------------#






   #---------------------------------------------------------------------------------------#
   #   Plot the filled contour diagrams showing time of day and time.                      #
   #---------------------------------------------------------------------------------------#
   for (v in 1:nsqueeze){

      #----- Retrieve variable information from the list. ---------------------------------#
      thisfillc   = squeeze[[v]]
      vnam        = thisfillc$vnam
      description = thisfillc$desc
      unit        = thisfillc$unit
      vcscheme    = thisfillc$col.scheme
      plotit      = thisfillc$fco.qmean

      if (plotit){

         #---------------------------------------------------------------------------------#
         #     Check if the directory exists.  If not, create it.                          #
         #---------------------------------------------------------------------------------#
         outdir  =  paste(outpref,"fillc_qmean",sep="/")
         if (! file.exists(outdir)) dir.create(outdir)
         cat("      + ",description," time series of diurnal cycle...","\n")

         #----- Load this variable into "thisvar". ----------------------------------------#
         vararr   = qmean[[vnam]]

         #----- Copy Decembers ans Januaries to make the edges buffered. ------------------#
         firsthr  = vararr[,1]
         firsthr  = c(firsthr,firsthr[ntimes],firsthr[ntimes])

         lasthr   = vararr[,ndcycle]
         lasthr   = c(lasthr[1],lasthr[1],lasthr)

         #----- Bind first and last year to the array, to make the edges buffered. --------#
         varbuff  = rbind(vararr[1,],vararr,vararr[ntimes,])
         varbuff  = cbind(lasthr,varbuff,firsthr)

         #----- Expand the month and year axes. -------------------------------------------#
         hraxis    = seq(from=0,to=ndcycle+1,by=1) * 24 / ndcycle
         dwhen     = datum$tomonth[2]-datum$tomonth[1]
         whenaxis  = c(datum$tomonth[1]-dwhen,datum$tomonth,datum$tomonth[ntimes]+dwhen)
         huplot    = pretty.time(whenaxis,n=8)

         vrange  = range(varbuff,na.rm=TRUE)
         vlevels = pretty(x=vrange,n=ncolsfc)
         vnlev   = length(vlevels)

         #----- Loop over formats. --------------------------------------------------------#
         for (o in 1:nout){
            fichier = paste(outdir,"/",vnam,"-",suffix,".",outform[o],sep="")
            if(outform[o] == "x11"){
               X11(width=size$width,height=size$height,pointsize=ptsz)
            }else if(outform[o] == "png"){
               png(filename=fichier,width=size$width*depth,height=size$height*depth
                  ,pointsize=ptsz,res=depth)
            }else if(outform[o] == "eps"){
               postscript(file=fichier,width=size$width,height=size$height
                         ,pointsize=ptsz,paper=size$paper)
            }else if(outform[o] == "pdf"){
               pdf(file=fichier,onefile=FALSE
                  ,width=size$width,height=size$height,pointsize=ptsz,paper=size$paper)
            }#end if

            letitre = paste("Mean diurnal cycle \n ",description," - ",lieu,sep="")
            par(par.user)
            sombreado(x=whenaxis,y=hraxis,z=varbuff,levels=vlevels,nlevels=vnlev
                     ,colour.palette=get(vcscheme)
                     ,plot.title=title(main=letitre,ylab="Time of day (GMT)"
                                      ,xlab="Time",cex.main=0.7)
                     ,key.title=title(main=unit,cex.main=0.8)
                     ,plot.axes={axis(side=1,at=huplot$level,labels=huplot$labels)
                                 axis(side=2,at=uplot$levels,labels=uplot$labels)
                                 if (fcgrid){
                                    abline(v=huplot$levels,h=uplot$levels
                                          ,col=grid.colour,lty="dotted")
                                 }#end if fcgrid
                                }#end plot.axes
                     )

            if (outform[o] == "x11"){
               locator(n=1)
               dev.off()
            }else{
               dev.off()
            }#end if
            dummy = clean.tmp()
         } #end for outform
      }#end if plotit
   }#end for nhov
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #   Plot the monthly boxplots.                                                          #
   #---------------------------------------------------------------------------------------#
   for (v in 1:nsqueeze){

      #----- Retrieve variable information from the list. ---------------------------------#
      thisbplot   = squeeze[[v]]
      vnam        = thisbplot$vnam
      description = thisbplot$desc
      unit        = thisbplot$unit
      plotit      = thisbplot$box.plot

      if (plotit){
         #---------------------------------------------------------------------------------#
         #     Check if the directory exists.  If not, create it.                          #
         #---------------------------------------------------------------------------------#
         outdir  =  paste(outpref,"boxplot",sep="/")
         if (! file.exists(outdir)) dir.create(outdir)
         cat("      + ",description," box plot...","\n")

         #----- Load this variable into "thisvar". ----------------------------------------#
         thisvar = emean[[vnam]]

         for (o in 1:nout){
            fichier = paste(outdir,"/",vnam,"-",suffix,".",outform[o],sep="")
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
                  ,width=size$width,height=size$height,pointsize=ptsz,paper=size$paper)
            }#end if


            ylimit  = pretty.xylim(u=thisvar,fracexp=0.0,is.log=FALSE)
            letitre = paste(description,lieu,sep=" - ")
            par(par.user)
            plot(mmonth,thisvar,main=letitre,ylim=ylimit,cex.main=0.7
                ,xlab="Time",ylab=paste("[",unit,"]",sep=""))

            if (outform[o] == "x11"){
               locator(n=1)
               dev.off()
            }else{
               dev.off()
            }#end if
            dummy = clean.tmp()
         } #end for outform
      }#end if
   }#end for nbox
   #---------------------------------------------------------------------------------------#





   #---------------------------------------------------------------------------------------#
   #   Plot the PDF of patch-level properties as a function of time.                       #
   #---------------------------------------------------------------------------------------#
   cat ("      * Time series of PDF of properties by patch...","\n")
   for (v in 1:nplotpatch){

      #----- Retrieve variable information from the list. ---------------------------------#
      thispatch   = plotpatch[[v]]
      vnam        = thispatch$vnam
      description = thispatch$desc
      unit        = thispatch$unit
      vcscheme    = thispatch$col.scheme
      plog        = thispatch$plog
      plotit      = thispatch$emean

      if (plotit){

         #---------------------------------------------------------------------------------#
         #     Check if the directory exists.  If not, create it.                          #
         #---------------------------------------------------------------------------------#
         outdir  =  paste(outpref,"patch_emean",sep="/")
         if (! file.exists(outdir)) dir.create(outdir)
         cat("      + PDF plot of ",description,"...","\n")

         this = patchpdf[[vnam]]$edensity


         #----- Loop over formats. --------------------------------------------------------#
         for (o in 1:nout){
            fichier = paste(outdir,"/",vnam,"-",suffix,".",outform[o],sep="")
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
                  ,width=size$width,height=size$height,pointsize=ptsz,paper=size$paper)
            }#end if

            letitre = paste("Density function of ",description," \ ",lieu,sep="")
            lex     = "Time"
            ley     = paste(description," [",unit,"]",sep="")


            #------------------------------------------------------------------------------#
            #     Plot the PDF distribution.                                               #
            #------------------------------------------------------------------------------#
            par(par.user)
            sombreado( x              = this$x
                     , y              = this$y
                     , z              = this$z
                     , nlevels        = ncolsfc
                     , colour.palette = get(vcscheme)
                     , plot.title     = title(main=letitre,xlab=lex,ylab=ley,cex.main=0.7)
                     , key.title      = title(main="Density",cex.main=0.8)
                     , key.log        = plog 
                     , plot.axes      = {  axis( side   = 1
                                               , at     = whenplot8$levels
                                               , labels = whenplot8$labels
                                               , padj   = whenplot8$padj
                                               )#end axis
                                           axis( side   = 2
                                               , at     = pretty(this$y)
                                               , labels = NULL
                                               )#end axis
                                           if (fcgrid){
                                              abline( v   = whenplot8$levels
                                                    , h   = pretty(this$y)
                                                    , col = grid.colour
                                                    , lty = "dotted"
                                                    )#end abline
                                           }#end if fcgrid
                                        }#end plot.axes
                     )#end sombreado
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #     Close the device.                                                        #
            #------------------------------------------------------------------------------#
            if (outform[o] == "x11"){
               locator(n=1)
               dev.off()
            }else{
               dev.off()
            }#end if
            #------------------------------------------------------------------------------#
            dummy = clean.tmp()
            #------------------------------------------------------------------------------#
         } #end for outform
         #---------------------------------------------------------------------------------#
      }#end if plotit
      #------------------------------------------------------------------------------------#
   }#end for (v in 1:npatchplot)
   #---------------------------------------------------------------------------------------#





   #---------------------------------------------------------------------------------------#
   #   Plot the PDF of patch-level properties as a function of time.                       #
   #---------------------------------------------------------------------------------------#
   cat ("      * Monthly PDF of properties by patch...","\n")
   for (v in 1:nplotpatch){

      #----- Retrieve variable information from the list. ---------------------------------#
      thispatch   = plotpatch[[v]]
      vnam        = thispatch$vnam
      description = thispatch$desc
      unit        = thispatch$unit
      vcscheme    = thispatch$col.scheme
      plog        = thispatch$plog
      plotit      = thispatch$mmean

      if (plotit){

         #---------------------------------------------------------------------------------#
         #     Check if the directory exists.  If not, create it.                          #
         #---------------------------------------------------------------------------------#
         outdir  =  paste(outpref,"patch_mmean",sep="/")
         if (! file.exists(outdir)) dir.create(outdir)
         cat("      + PDF plot of ",description,"...","\n")

         this = patchpdf[[vnam]]$mdensity


         #----- Find the month tick marks. ------------------------------------------------#
         monat  = 1:12
         monlab = c("J","F","M","A","M","J","J","A","S","O","N","D")
         #---------------------------------------------------------------------------------#


         #----- Loop over formats. --------------------------------------------------------#
         for (o in 1:nout){
            fichier = paste(outdir,"/",vnam,"-",suffix,".",outform[o],sep="")
            if(outform[o] == "x11"){
               X11(width=size$width,height=size$height,pointsize=ptsz)
            }else if(outform[o] == "png"){
               png(filename=fichier,width=size$width*depth,height=size$height*depth
                  ,pointsize=ptsz,res=depth)
            }else if(outform[o] == "eps"){
               postscript(file=fichier,width=size$width,height=size$height
                         ,pointsize=ptsz,paper=size$paper)
            }else if(outform[o] == "pdf"){
               pdf(file=fichier,onefile=FALSE
                  ,width=size$width,height=size$height,pointsize=ptsz,paper=size$paper)
            }#end if

            letitre = paste("Density function of ",description," \ ",lieu,sep="")
            lex     = "Months"
            ley     = paste(description," [",unit,"]",sep="")


            #------------------------------------------------------------------------------#
            #     Plot the PDF distribution.                                               #
            #------------------------------------------------------------------------------#
            par(par.user)
            sombreado( x              = this$x
                     , y              = this$y
                     , z              = this$z
                     , nlevels        = ncolsfc
                     , colour.palette = get(vcscheme)
                     , plot.title     = title(main=letitre,xlab=lex,ylab=ley,cex.main=0.7)
                     , key.title      = title(main="Density",cex.main=0.8)
                     , key.log        = plog 
                     , plot.axes      = {  axis( side   = 1
                                               , at     = monat
                                               , labels = monlab
                                               )#end axis
                                           axis( side   = 2
                                               , at     = pretty(this$y)
                                               , labels = NULL
                                               )#end axis
                                           if (fcgrid){
                                              abline( v   = monat
                                                    , h   = pretty(this$y)
                                                    , col = grid.colour
                                                    , lty = "dotted"
                                                    )#end abline
                                           }#end if fcgrid
                                        }#end plot.axes
                     )#end sombreado
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #     Close the device.                                                        #
            #------------------------------------------------------------------------------#
            if (outform[o] == "x11"){
               locator(n=1)
               dev.off()
            }else{
               dev.off()
            }#end if
            #------------------------------------------------------------------------------#
            dummy = clean.tmp()
            #------------------------------------------------------------------------------#
         } #end for outform
         #---------------------------------------------------------------------------------#
      }#end if plotit
      #------------------------------------------------------------------------------------#
   }#end for (v in 1:npatchplot)
   #---------------------------------------------------------------------------------------#






   #---------------------------------------------------------------------------------------#
   #      Bar plot by DBH class.                                                           #
   #---------------------------------------------------------------------------------------#
   cat("    + Bar plot by DBH classes...","\n")
   monbplot    = which(nummonths(datum$tomonth) %in% sasmonth)
   nmonbplot   = length(monbplot)
   pftuse      = which(apply(X=szpft$nplant,MARGIN=3,FUN=sum,na.rm=TRUE) > 0.)
   pftuse      = pftuse[pftuse != (npft+1)]
   npftuse     = length(pftuse)
   pftname.use = pft$name  [pftuse]
   pftcol.use  = pft$colour[pftuse]
   for (v in 1:ntspftdbh){
      #----- Load settings for this variable.----------------------------------------------#
      thisbar     = tspftdbh[[v]]
      vnam        = thisbar$vnam
      description = thisbar$desc
      unit        = thisbar$unit
      stacked     = thisbar$stack
      plotit      = thisbar$bar.plot
      plog        = thisbar$plog
      if (plog){
         xylog   = "y"
         stacked = FALSE
      }else{
         xylog   = ""
      }#end if
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #      Check whether to plot this 
      #------------------------------------------------------------------------------------#
      if (plotit){
         cat("      - ",description,"...","\n")


         #---------------------------------------------------------------------------------#
         #     Retrieve the variable, and keep only the part that is usable.               #
         #---------------------------------------------------------------------------------#
         thisvnam                  = szpft[[vnam]][monbplot,,]
         thisvnam                  = thisvnam [,,pftuse]
         thisvnam                  = thisvnam [,-(ndbh+1),]
         
         thisvnam[is.na(thisvnam)] = 0.
         thiswhen                  = datum$tomonth[monbplot]
         #---------------------------------------------------------------------------------#
       

         #---------------------------------------------------------------------------------#
         #      Find the limits for the plots.  We use the same axis so it is easier to    #
         # compare different times.                                                        #
         #---------------------------------------------------------------------------------#
         if (stacked){
            ylimit   = c(0,max(apply(X=thisvnam,MARGIN=c(1,2),FUN=sum,na.rm=TRUE)))
         }else{
            ylimit   = range(x=thisvnam,na.rm=TRUE)
         }#end if
         ylimit = pretty.xylim(u=ylimit,fracexp=scalleg,is.log=plog)
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Check if the directory exists.  If not, create it.                          #
         #---------------------------------------------------------------------------------#
         barplotdir = paste(outpref,"barplot_dbh",sep="/")
         if (! file.exists(barplotdir)) dir.create(barplotdir)
         outdir = paste(barplotdir,vnam,sep="/")
         if (! file.exists(outdir)) dir.create(outdir)
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Loop over all possible months.                                             #
         #---------------------------------------------------------------------------------#
         for (m in 1:nmonbplot){

            #----- Find which year we are plotting. ---------------------------------------#
            cmonth    = sprintf("%2.2i",(nummonths(thiswhen[m])))
            cyear     = sprintf("%4.4i",(numyears(thiswhen[m])))
            mm        = as.numeric(cmonth)
            yy        = as.numeric(cyear)
            whentitle = paste(mon2mmm(mm,cap1=TRUE),cyear,sep="-")
            #------------------------------------------------------------------------------#


            #----- Loop over output formats. ----------------------------------------------#
            for (o in 1:nout){
               #------ Open the plot. -----------------------------------------------------#
               fichier = paste(outdir,"/",vnam,"-",cyear,"-",cmonth,"-",suffix
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
                     ,width=size$width,height=size$height,pointsize=ptsz,paper=size$paper)
               }#end if
               #---------------------------------------------------------------------------#


               #------ Set up the title and axis labels. ----------------------------------#
               letitre = paste(lieu,"\n",description," - Time : ",whentitle,sep="")
               lexlab  = "DBH Classes"
               leylab  = paste(description," [",unit,"]",sep="")
               #---------------------------------------------------------------------------#


               #----- Plot all monthly means together. ------------------------------------#
               par(par.user)
               barplot(height=t(thisvnam[m,,]),names.arg=dbhnames[1:ndbh],width=1.0
                      ,main=letitre,xlab=lexlab,ylab=leylab,ylim=ylimit,legend.text=FALSE
                      ,beside=(! stacked),col=pftcol.use,log=xylog
                      ,border=grey.fg,xpd=FALSE,cex.main=cex.main)
               if (plotgrid & (! stacked)){
                  xgrid=0.5+(1:ndbh)*(1+npftuse)
                  abline(v=xgrid,col=grid.colour,lty="solid")
               }#end if
               box()
               legend( x      = "topleft"
                     , inset  = inset
                     , legend = pftname.use
                     , fill   = pftcol.use
                     , ncol   = min(3,pretty.box(n.selpft)$ncol)
                     , title  = expression(bold("Plant functional type"))
                     , cex    = 1.0
                     , bg     = background
                     )#end legend
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     Close the device.                                                     #
               #---------------------------------------------------------------------------#
               if (outform[o] == "x11"){
                  locator(n=1)
                  dev.off()
               }else{
                  dev.off()
               }#end if
               dummy = clean.tmp()
               #---------------------------------------------------------------------------#
            } #end for outform
            #------------------------------------------------------------------------------#
         }#end for
         #---------------------------------------------------------------------------------#
      }#end if
      #------------------------------------------------------------------------------------#
   }#end for
   #---------------------------------------------------------------------------------------#





   #---------------------------------------------------------------------------------------#
   #    Plot the 3-D size and age structure of light level.                                #
   #---------------------------------------------------------------------------------------#
   for (v in 1:ntspftdbh){
      #----- Retrieve variable information from the list. ---------------------------------#
      thissas     = tspftdbh[[v]]
      vnam        = thissas$vnam
      description = thissas$desc
      unit        = thissas$i.unit
      plotit      = thissas$sas

      #----- If this variable is to be plotted, then go through this if block. ------------#
      if (plotit){

         cat("      + ",description," size and age structure plot...","\n")

         #---------------------------------------------------------------------------------#
         #     Check if the directory exists.  If not, create it.                          #
         #---------------------------------------------------------------------------------#
         sasdir = paste(outpref,"sas",sep="/")
         if (! file.exists(sasdir)) dir.create(sasdir)
         outdir = paste(sasdir,vnam,sep="/")
         if (! file.exists(outdir)) dir.create(outdir)

         #----- Load this list into "thislist". -------------------------------------------#
         varco =  cohort[[vnam]]


         for (ww in names(cohort$age)){

            #----- Find which year we are plotting. ---------------------------------------#
            cmonth   = substring(ww,7,8)
            thisyear = substring(ww,2,5)
            mm       = as.numeric(cmonth)
            yy       = as.numeric(thisyear)

            #----- Retrieve variable list, age, DBH, and PFT for this year. ---------------#
            ageww   = cohort$age   [[ww]]
            dbhww   = cohort$dbh   [[ww]]
            pftww   = cohort$pft   [[ww]]
            varww   = varco        [[ww]]
            popww   = cohort$nplant[[ww]] * cohort$area[[ww]]

            #------------------------------------------------------------------------------#
            #     We only plot the SAS figures when the polygon is not an absolute desert. #
            #------------------------------------------------------------------------------#
            if (any (! is.na(varww))){
               #---------------------------------------------------------------------------#
               #      Find the range.  If the user wants the range to be fixed, then use   #
               # the global range, otherwise, simply use the range for this year.          #
               #---------------------------------------------------------------------------#
               if (sasfixlimits){
                  xlimit  = range(unlist(cohort$age)                 , na.rm=TRUE)
                  ylimit  = range(unlist(cohort$dbh)                 , na.rm=TRUE)
                  zlimit  = range(unlist(varco)                      , na.rm=TRUE)
                  popmin  = min  (unlist(cohort$nplant * cohort$area), na.rm=TRUE)
                  popmax  = max  (unlist(cohort$nplant * cohort$area), na.rm=TRUE)
               }else{
                  xlimit  = range(ageww  ,na.rm=TRUE)
                  ylimit  = range(dbhww  ,na.rm=TRUE)
                  zlimit  = range(varww  ,na.rm=TRUE)
                  popmin  = min  (popww  ,na.rm=TRUE)
                  popmax  = max  (popww  ,na.rm=TRUE)
               }#end if

               #----- Define the scale-dependent population size. -------------------------#
               cexww = cexmin + (cexmax - cexmin) * log(popww/popmin) / log(popmax/popmin)

               #----- Define the floor location. ------------------------------------------#
               if (zlimit[1] == zlimit[2]){
                  if (zlimit[1] == 0){
                     zlimit = c(-1.,1.)
                  }else{
                     zlimit = sort(c(0.9,1.1)*zlimit[1])
                  }#end if
               }#end if
               if ((zlimit[1] > 0) != (zlimit[2] > 0)){
                  floor3d = 0.
               }else if (zlimit[1] > 0){
                  floor3d = zlimit[1]
               }else{
                  floor3d = zlimit[2]
               }#end if

               #----- Define the grid information for the 3-D plot. -----------------------#
               ageaxis   = pretty(xlimit,n=20)
               dbhaxis   = pretty(ylimit,n=20)
               xlimit    = range(ageaxis)
               ylimit    = range(dbhaxis)
               flooraxis = matrix(floor3d,nrow=length(ageaxis),ncol=length(dbhaxis))

               #----- Expand the lines to make the lollipops. -----------------------------#
               ncohnow  = length(varww)
               ageww    = rep(ageww,each=3)
               dbhww    = rep(dbhww,each=3)
               pftww    = rep(pftww,each=3)
               varww    = as.vector(rbind(rep(floor3d,times=ncohnow)
                                         ,varco[[ww]]
                                         ,rep(NA,times=ncohnow)))
               pchww    = rep(c(NA,16,NA),times=ncohnow)
               cexww    = rep(cexww,each=3)
               colww    = pft$colour[pftww]

               pftin   = sort(unique(cohort$pft[[ww]]))
               colleg  = pft$colour[pftin]
               pftleg  = pft$name  [pftin]


               #----- Loop over output formats. -------------------------------------------#
               for (o in 1:nout){
                  fichier = paste(outdir,"/",vnam,"-",thisyear,"-",cmonth,"-",suffix
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

                  stcol   = pft$colour[pftww]
                  letitre = paste(description," - ",lieu,
                                  "\n Time :",mlist[mm],"/",thisyear,sep=" ")
                  lezlab  = paste(description," [",unit,"]",sep="")

                  #----- First plot: the box. ---------------------------------------------#
                  par(par.user)
                  pout = persp(x=ageaxis,y=dbhaxis,z=flooraxis,xlim=xlimit,ylim=ylimit
                              ,zlim=zlimit,theta=theta,phi=phi,col=gcol,expand=expz
                              ,ticktype="detailed",border=NA,xlab="Gap age [yr]"
                              ,ylab="DBH [cm]",zlab=lezlab,shade=shade,ltheta=ltheta
                              ,main=letitre,cex.main=0.7)
                  #----- Second plot, the actual data (aka my lollipop trees). ------------#
                  lines (trans3d(x=ageww,y=dbhww,z=varww,pmat=pout),type="l"
                        ,col=grey.fg,lwd=2)
                  points(trans3d(x=ageww,y=dbhww,z=varww,pmat=pout),type="p"
                        ,pch=pchww,col=colww,cex=cexww)
                  legend( x      = "bottomright"
                        , inset  = inset
                        , legend = pftleg
                        , fill   = colleg
                        , ncol   = 1
                        , title  = expression(bold("PFT"))
                        , bg     = background
                        , cex    = 0.9
                        )#end legend


                  if (outform[o] == "x11"){
                     locator(n=1)
                     dev.off()
                  }else{
                     dev.off()
                  }#end if
                  dummy = clean.tmp()
               } #end for outform
            }#end if is.na(varww)
         }#end for nameco
      } #end if
   }#end for npsas
   #---------------------------------------------------------------------------------------#
}#end for places
#==========================================================================================#
#==========================================================================================#
