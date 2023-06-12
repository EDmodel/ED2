
#----- Paths. -----------------------------------------------------------------------------#
here           = "thispath"    # Current directory.
there          = "thatpath"    # Directory where analyses/history are 
srcdir         = "thisrscpath" # Source  directory.
outroot        = "thisoutroot" # Directory for figures
#------------------------------------------------------------------------------------------#


#----- Time options. ----------------------------------------------------------------------#
monthbeg       = thismontha   # First month to use
yearbeg        = thisyeara    # First year to consider
yearend        = thisyearz    # Maximum year to consider
reload.data    = TRUE         # Should I reload partially loaded data?
sasmonth.short = c(2,5,8,11)  # Months for SAS plots (short runs)
sasmonth.long  = 5            # Months for SAS plots (long runs)
nyears.long    = 15           # Runs longer than this are considered long runs.
#------------------------------------------------------------------------------------------#



#----- Name of the simulations. -----------------------------------------------------------#
myplaces       = c("thispoly")
#------------------------------------------------------------------------------------------#



#----- Plot options. ----------------------------------------------------------------------#
outform        = thisoutform            # Formats for output file.  Supported formats are:
                                        #   - "X11"    - for printing on screen
                                        #   - "quartz" - for printing on Mac OS screen
                                        #   - "eps"    - for postscript printing
                                        #   - "png"    - for PNG printing
                                        #   - "tif"    - for TIFF printing
                                        #   - "pdf"    - for PDF printing
depth          = 96                     # PNG resolution, in pixels per inch
paper          = "letter"               # Paper size, to define the plot shape
ptsz           = 16                     # Font size.
lwidth         = 2.5                    # Line width
plotgrid       = TRUE                   # Should I plot the grid in the background? 
sasfixlimits   = FALSE                  # Use a fixed scale for size and age-structure
                                        #    plots? (FALSE will set a suitable scale for
                                        #    each plot)
fcgrid         = TRUE                   # Include a grid on the filled contour plots?
ncolshov       = 200                    # Target number of colours for Hovmoller diagrams.
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
ylnudge        = 0.05                   # Nudging factor for ylimit
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
slz.min             = -5.0         # The deepest depth that trees access water.
idbh.type           = myidbhtype   # Type of DBH class
                                   # 1 -- Every 10 cm until 100cm; > 100cm
                                   # 2 -- 0-10; 10-20; 20-35; 35-50; 50-70; > 70 (cm)
                                   # 3 -- 0-10; 10-35; 35-55; > 55 (cm)
klight              = myklight     # Weighting factor for maximum carbon balance
corr.growth.storage = mycorrection # Correction factor to be applied to growth and
                                   #   storage respiration
iallom              = myallom      # Allometry to use
isoil.hydro         = myslhydro    # Soil hydrology method
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
nout    = length (outform)
#------------------------------------------------------------------------------------------#


#----- Avoid unecessary and extremely annoying beeps. -------------------------------------#
options(locatorBell=FALSE)
#------------------------------------------------------------------------------------------#


#----- Load observations. -----------------------------------------------------------------#
obsrfile = file.path(srcdir,"LBA_MIP.v9.RData")
load(file=obsrfile)

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
   outmain = file.path(outroot,place)
   outcsv  = file.path(outmain,"csv_month")
   lieu    = thispoi$lieu
   iata    = thispoi$iata
   suffix  = thispoi$iata
   yeara   = thispoi$yeara
   yearz   = thispoi$yearz
   meszz   = thispoi$monz

   #---------------------------------------------------------------------------------------#
   #     Make sure we only deal with full years.                                           #
   #---------------------------------------------------------------------------------------#
   if (monthbeg >  1) yeara = yeara + 1
   if (meszz    < 12) yearz = yearz - 1
   monthbeg = 1
   meszz    = 12
   if (yeara > yearz){
      cat(" - Yeara:  ",yeara,"\n")
      cat(" - Yearz:  ",yearz,"\n")
      cat(" - Prefix: ",inpref,"\n")
      cat(" - Invalid years, will not process data...","\n")
      q("no")
   }#end if
   #---------------------------------------------------------------------------------------#


   #----- Create the directories in case they don't exist. --------------------------------#
   if (! file.exists(outmain)) dir.create(outmain)
   if (! file.exists(outcsv )) dir.create(outcsv )
   #---------------------------------------------------------------------------------------#



   #----- Decide how frequently the cohort-level variables should be saved. ---------------#
   if ((yearend - yearbeg + 1) <= nyears.long){
      sasmonth   = sasmonth.short
      plot.ycomp = TRUE
   }else{
      sasmonth   = sasmonth.long
      plot.ycomp = FALSE
   }#end if
   #---------------------------------------------------------------------------------------#




   #----- Print a banner to entretain the user. -------------------------------------------#
   cat(" + Post-processing output from ",lieu,"...","\n")


   #---------------------------------------------------------------------------------------#
   #     Flush all variables that will hold the data.                                      #
   #---------------------------------------------------------------------------------------#
   ntimes      = (yearz-yeara-1)*12+meszz+(12-monthbeg+1)
   nyears      =  yearz-yeara+1
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      Make the RData file name, then we check whether we must read the files again     #
   # or use the stored RData.  Notice that the path is the same for plot_ycomp.r and       #
   # plot_monthly, so you don't need to read in the data twice.                            #
   #---------------------------------------------------------------------------------------#
   path.data  = file.path(here,place,"rdata_month")
   if (! file.exists(path.data)) dir.create(path.data)
   ed22.rdata  = file.path(path.data,paste(place,"RData",sep="."))
   ed22.status = file.path(path.data,paste("status_",place,".txt",sep=""))
   if (reload.data && file.exists(ed22.rdata)){
      #----- Load the modelled dataset. ---------------------------------------------------#
      cat("   - Loading previous session...","\n")
      load(ed22.rdata)
      tresume = datum$ntimes + 1
      if (ntimes > datum$ntimes){
         datum   = update.monthly( new.ntimes = ntimes
                                 , old.datum  = datum
                                 , montha     = monthbeg
                                 , yeara      = yeara
                                 , inpref     = inpref
                                 , slz.min    = slz.min
                                 )#end update.monthly
      }#end if
      #------------------------------------------------------------------------------------#
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
   
   
   #----- Update status file with latest data converted into R. ---------------------------#
   latest = paste(datum$year[ntimes],datum$month[ntimes],sep=" ")
   dummy  = write(x=latest,file=ed22.status,append=FALSE)
   #---------------------------------------------------------------------------------------#




   #----- Make some shorter versions of some variables. -----------------------------------#
   mfac   = datum$month
   yfac   = datum$year
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
   #     Remove all elements of the DBH/PFT class that do not have a single valid cohort   #
   # at any given time.                                                                    #
   #---------------------------------------------------------------------------------------#
   empty = is.na(szpft$nplant) | szpft$nplant == 0
   for (vname in names(szpft)) szpft[[vname]][empty] = NA
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Convert mortality and recruitment so it is scaled between 0 and 100%.             #
   #---------------------------------------------------------------------------------------#
   struct    = c("szpft","emean","mmean","ymean")
   struct    = struct[struct %in% ls()]
   nstruct   = length(struct)
   mort.list = c(    "mort",    "dimort",    "ncbmort",    "hydmort","fire.lethal"
                 ,"agb.mort","agb.dimort","agb.ncbmort","agb.hydmort"
                 ,"bsa.mort","bsa.dimort","bsa.ncbmort","bsa.hydmort"
                 )#end c
   recr.list = c(    "recr","agb.recr","bsa.recr")
   for (s in sequence(nstruct)){
      #----- Copy structure to a temporary variable. --------------------------------------#
      stnow = struct[s]
      xmean = get(struct[s])
      #------------------------------------------------------------------------------------#


      #----- Select mortality and recruitment variables to update. ------------------------#
      mort.check  = mort.list[mort.list %in% names(xmean)]
      nmort.check = length(mort.check)
      recr.check  = recr.list[recr.list %in% names(xmean)]
      nrecr.check = length(recr.check)
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Loop through mortality variables, make them "interest rates".                  #
      #------------------------------------------------------------------------------------#
      for (m in sequence(nmort.check)){
         mort.now = mort.check[m]
         xmean[[mort.now]] = 100. * ( 1.0 - exp( - xmean[[mort.now]]) )
      }#end for (m in sequence(nmort.check))
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Loop through recruitment variables, make them "interest rates".                #
      #------------------------------------------------------------------------------------#
      for (r in sequence(nrecr.check)){
         recr.now = recr.check[r]
         xmean[[recr.now]] = 100. * ( exp( + xmean[[recr.now]] ) - 1.0)
      }#end for (m in sequence(nmort.check))
      #------------------------------------------------------------------------------------#


      #------ Update structure. -----------------------------------------------------------#
      dummy = assign(x=stnow,value=xmean)
      #------------------------------------------------------------------------------------#
   }#end for (s in seq_along(struct))
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




   #---------------------------------------------------------------------------------------#
   #     Create data frame with yearly averages.                                           #
   #---------------------------------------------------------------------------------------#
   cat0(" + Create data frame with monthly averages.")
   mglob  = data.table( year          = datum$year
                      , month         = datum$month
                      , numdays       = daymax(datum$month,datum$year)
                      , atm.prss      = emean$atm.prss      * 100.
                      , atm.temp      = emean$atm.temp      + t00
                      , atm.shv       = emean$atm.shv       * 0.001
                      , atm.vels      = emean$atm.vels
                      , precip        = emean$rain/daymax(datum$month,datum$year)/day.sec
                      , rshort.in     = emean$rshort
                      , rlong.in      = emean$rlong
                      , agb           = emean$agb
                      , lai           = emean$lai
                      , wai           = emean$wai
                      , rshort.out    = emean$rshort
                      , rlong.out     = emean$rlong
                      , cas.height    = emean$can.depth
                      , cas.prss      = emean$can.prss      * 100.
                      , cas.temp      = emean$can.temp      + t00
                      , cas.shv       = emean$can.shv       * 0.001
                      , gnd.temp      = emean$gnd.temp      + t00
                      , leaf.temp     = emean$leaf.temp     + t00
                      , wood.temp     = emean$wood.temp     + t00
                      , soil.temp     = emean$soil.temp.top + t00
                      , soil.water    = emean$soil.water.top
                      , soil.rmoist   = emean$soil.wetness.top
                      , sfcw.temp     = emean$sfcw.temp     + t00
                      , sfcw.fliq     = emean$sfcw.fliq
                      , sfcw.mass     = emean$sfcw.mass
                      , sfcw.depth    = emean$sfcw.depth
                      , sfcw.cover    = emean$sfcw.cover
                      , runoff        = emean$runoff / day.sec
                      , leaf.water    = emean$leaf.water
                      , wood.water    = 0. * emean$leaf.water
                      , hflxca        = emean$hflxca
                      , wflxca        = emean$wflxca / day.sec
                      , hflxgc        = emean$hflxgc
                      , wflxgc        = emean$wflxgc / day.sec
                      , hflxlc        = emean$hflxlc
                      , wflxlc        = emean$wflxlc / day.sec
                      , hflxwc        = emean$hflxwc
                      , wflxwc        = emean$wflxwc / day.sec
                      , transp        = emean$transp / day.sec
                      )#end data.frame
   #---------------------------------------------------------------------------------------#

   #---------------------------------------------------------------------------------------#
   #    Write csv files with the output.                                                   #
   #---------------------------------------------------------------------------------------#
   cat0(" + Write output files.")
   csv.emean = file.path(outcsv,paste0(place,"_emean.csv" ))
   dummy = write.table(x=mglob ,file=csv.emean,quote=FALSE,sep=",",row.names=FALSE)
   #---------------------------------------------------------------------------------------#

}#end for (place in myplaces)
#------------------------------------------------------------------------------------------#
