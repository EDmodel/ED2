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
   outcsv  = file.path(outmain,"csv_year")
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
   #     Consolidate the yearly means for the long-term dynamics (the PFT and DBH/PFT      #
   # stuff).                                                                               #
   #---------------------------------------------------------------------------------------#
   cat0("    - Find the annual statistics for multi-dimensional variables...")
   cat0("      * Aggregate the annual mean of PFT-DBH variables...")
   for (vname in names(szpft)){
      szpft[[vname]] = qapply(X=szpft[[vname]],INDEX=yfac,DIM=1,FUN=mean,na.rm=TRUE)
   }#end for
   #----- LU arrays.   The "+1" column contains the total. --------------------------------#
   cat0("      * Aggregate the annual mean of LU variables...")
   for (vname in names(lu)){
      lu   [[vname]] = qapply(X=lu   [[vname]],INDEX=yfac,DIM=1,FUN=mean,na.rm=TRUE)
   }#end for
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Here we find the monthly means for month, then compute the standard deviation.   #
   #---------------------------------------------------------------------------------------#
   cat0("    - Find the monthly and annual means...")
   cat0("      * Aggregating the monthly mean and standard deviation...")
   mmean = list()
   msdev = list()
   ymean = list()
   ysdev = list()
   for (vname in names(emean)){
      if (vname %in% c("soil.temp","soil.water","soil.mstpot","soil.extracted")){
         mmean[[vname]] = qapply(X=emean[[vname]], INDEX=mfac, DIM=1, FUN=mean, na.rm=TRUE)
         msdev[[vname]] = qapply(X=emean[[vname]], INDEX=mfac, DIM=1, FUN=sd  , na.rm=TRUE)
         ymean[[vname]] = qapply(X=emean[[vname]], INDEX=yfac, DIM=1, FUN=mean, na.rm=TRUE)
         ysdev[[vname]] = qapply(X=emean[[vname]], INDEX=yfac, DIM=1, FUN=sd  , na.rm=TRUE)
      }else if (vname %in% c("rain","runoff","intercepted","wshed")){
         mmean[[vname]] = tapply(X=emean[[vname]], INDEX=mfac, FUN=mean, na.rm=TRUE)
         msdev[[vname]] = tapply(X=emean[[vname]], INDEX=mfac, FUN=sd  , na.rm=TRUE)
         ymean[[vname]] = tapply(X=emean[[vname]], INDEX=yfac, FUN=sum , na.rm=TRUE)
         ysdev[[vname]] = tapply(X=emean[[vname]], INDEX=yfac, FUN=sd  , na.rm=TRUE)
      }else{
         mmean[[vname]] = tapply(X=emean[[vname]], INDEX=mfac, FUN=mean, na.rm=TRUE)
         msdev[[vname]] = tapply(X=emean[[vname]], INDEX=mfac, FUN=sd  , na.rm=TRUE)
         ymean[[vname]] = tapply(X=emean[[vname]], INDEX=yfac, FUN=mean, na.rm=TRUE)
         ysdev[[vname]] = tapply(X=emean[[vname]], INDEX=yfac, FUN=sd  , na.rm=TRUE)
      }#end if
      #------------------------------------------------------------------------------------#


      #----- Fix the bad data. ------------------------------------------------------------#
      bad.mmean = ! is.finite(mmean[[vname]])
      bad.msdev = ! is.finite(msdev[[vname]])
      bad.ymean = ! is.finite(ymean[[vname]])
      bad.ysdev = ! is.finite(ysdev[[vname]])
      mmean[[vname]][bad.mmean] = NA
      msdev[[vname]][bad.msdev] = 0.
      ymean[[vname]][bad.ymean] = NA
      ysdev[[vname]][bad.ysdev] = 0.
      #------------------------------------------------------------------------------------#
   }#end for
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #      Here we find the Mean diurnal cycle for each month, then compute the standard    #
   # deviation.                                                                            #
   #---------------------------------------------------------------------------------------#
   cat0("    - Aggregate the annual mean and std. dev. of the diurnal cycle...")
   umean = list()
   usdev = list()
   for (vname in names(qmean)){
      umean[[vname]] = qapply(qmean[[vname]],INDEX=yfac,DIM=1,FUN=mean,na.rm=TRUE)
      usdev[[vname]] = qapply(qmean[[vname]],INDEX=yfac,DIM=1,FUN=sd  ,na.rm=TRUE)
      bad.umean      = ! is.finite(umean[[vname]])
      bad.usdev      = ! is.finite(usdev[[vname]])
      umean[[vname]][bad.umean] = NA
      usdev[[vname]][bad.usdev] = 0.
   }#end for
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
   szpft$mort          = 100. * (1.0 - exp(- szpft$mort         )      )
   szpft$dimort        = 100. * (1.0 - exp(- szpft$dimort       )      )
   szpft$ncbmort       = 100. * (1.0 - exp(- szpft$ncbmort      )      )
   szpft$recrpft       = 100. * (      exp(  szpft$recr         ) - 1.0)
   szpft$agb.mort      = 100. * (1.0 - exp(- szpft$agb.mort     )      )
   szpft$agb.dimort    = 100. * (1.0 - exp(- szpft$agb.dimort   )      )
   szpft$agb.ncbmort   = 100. * (1.0 - exp(- szpft$agb.ncbmort  )      )
   szpft$agb.recrpft   = 100. * (      exp(  szpft$agb.recr     ) - 1.0)
   szpft$bsa.mort      = 100. * (1.0 - exp(- szpft$bsa.mort     )      )
   szpft$bsa.dimort    = 100. * (1.0 - exp(- szpft$bsa.dimort   )      )
   szpft$bsa.ncbmort   = 100. * (1.0 - exp(- szpft$bsa.ncbmort  )      )
   szpft$bsa.recrpft   = 100. * (      exp(  szpft$bsa.recr     ) - 1.0)
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
   cat0(" + Create data frame with yearly averages.")
   ymean$year = as.numeric(names(ymean$agb))
   yglob  = data.frame( year          = ymean$year
                      , atm.temp      = ymean$atm.temp      + t00
                      , atm.vpd       = ymean$atm.vpd       * 100.
                      , atm.co2       = ymean$atm.co2
                      , rshort.in     = ymean$rshort
                      , rlong.in      = ymean$rlong
                      , par.in        = ymean$par.tot       * 1e-6
                      , rain          = ymean$rain
                      , runoff        = ymean$runoff
                      , agb           = ymean$agb
                      , lai           = ymean$lai
                      , nplant        = ymean$nplant
                      , cba           = ymean$cba
                      , fast.soil.c   = ymean$fast.soil.c
                      , struct.soil.c = ymean$struct.soil.c
                      , slow.soil.c   = ymean$slow.soil.c
                      , sensible.flux = ymean$hflxca
                      , vapour.flux   = ymean$wflxca        / day.sec
                      , latent.flux   = ymean$qwflxca
                      , co2.flux      = ymean$cflxca        * 1.0e-6
                      , co2.storage   = ymean$cflxst        * 1.0e-6
                      , nee           = ymean$nee           * 1.0e-6
                      , rshort.out    = ymean$rshortup
                      , rlong.out     = ymean$rlongup
                      , par.out       = ymean$parup         * 1.0e-6
                      , can.temp      = ymean$can.temp      + t00
                      , can.vpd       = ymean$can.vpd       * 100.
                      , leaf.temp     = ymean$leaf.temp     + t00
                      , leaf.vpd      = ymean$leaf.vpd      * 100.
                      , leaf.gsw      = ymean$leaf.gsw      / day.sec
                      , leaf.par      = ymean$par.leaf      * 1.0e-6
                      , evap          = ymean$evap          / day.sec
                      , transp        = ymean$transp        / day.sec
                      , gpp           = ymean$gpp
                      , auto.resp     = ymean$plant.resp
                      , het.resp      = ymean$het.resp
                      , mco           = ymean$mco
                      , growth        = ymean$acc.growth
                      , di.mortality  = ymean$acc.dimort
                      , ncb.mortality = ymean$acc.ncbmort
                      , recruitment   = ymean$acc.recr
                      , smpot         = ymean$smpot         * (-1. * 1.e6)
                      )#end data.frame
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Create data frame with yearly averages of size and age.                           #
   #---------------------------------------------------------------------------------------#
   cat0(" + Create data frame with yearly averages by PFT.")
   szpft$year       = array(data=ymean$year,dim=dim(szpft$agb),dimnames=dimnames(szpft))
   szpft.idx        = data.frame(arrayInd(seq_along(szpft$year),.dim=dim(szpft$year)))
   names(szpft.idx) = c("year","idbh","ipft")
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #    List variables by PFT.                                                             #
   #---------------------------------------------------------------------------------------#
   sel   = with(szpft.idx, idbh %in% (ndbh+1) & ipft %in% which(selpft & ipft <= npft))
   ypft  = data.frame( year             = c(szpft$year)[sel]
                     , ipft             = pft$key[c(szpft.idx$ipft     )[sel]]
                     , agb              = c(szpft$agb          )[sel]
                     , ba               = c(szpft$ba           )[sel]
                     , lai              = c(szpft$lai          )[sel]
                     , nplant           = c(szpft$nplant       )[sel]
                     , cba              = c(szpft$cba          )[sel]
                     , leaf.temp        = c(szpft$leaf.temp    )[sel] + t00
                     , leaf.vpd         = c(szpft$leaf.vpd     )[sel] * 100.
                     , leaf.gsw         = c(szpft$leaf.gsw     )[sel] / day.sec
                     , leaf.par         = c(szpft$par.leaf     )[sel] * 1.0e-6
                     , transp           = c(szpft$transp       )[sel] / day.sec
                     , gpp              = c(szpft$gpp          )[sel]
                     , auto.resp        = c(szpft$plant.resp   )[sel]
                     , mco              = c(szpft$mco          )[sel]
                     , l.transp         = c(szpft$i.transp     )[sel] / day.sec
                     , i.gpp            = c(szpft$i.gpp        )[sel]
                     , i.auto.resp      = c(szpft$i.plant.resp )[sel]
                     , i.mco            = c(szpft$i.mco        )[sel]
                     , i.cba            = c(szpft$cba          )[sel]
                     , growth           = c(szpft$acc.growth   )[sel]
                     , di.mortality     = c(szpft$acc.dimort   )[sel]
                     , ncb.mortality    = c(szpft$acc.ncbmort  )[sel]
                     , recruitment      = c(szpft$acc.recr     )[sel]
                     , cue              = c(szpft$cue          )[sel]
                     , wue              = c(szpft$wue          )[sel]
                     , stringsAsFactors = FALSE
                     )#end data.frame
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #    List variables by DBH class.                                                       #
   #---------------------------------------------------------------------------------------#
   cat0(" + Create data frame with yearly averages by DBH class.")
   sel   = with( szpft.idx, idbh <= ndbh & ipft %in% (npft+1))
   ydbh  = data.frame( year             = c(szpft$year)[sel]
                     , idbh             = dbhkeys[c(szpft.idx$idbh     )[sel]]
                     , agb              = c(szpft$agb          )[sel]
                     , ba               = c(szpft$ba           )[sel]
                     , lai              = c(szpft$lai          )[sel]
                     , nplant           = c(szpft$nplant       )[sel]
                     , cba              = c(szpft$cba          )[sel]
                     , leaf.temp        = c(szpft$leaf.temp    )[sel] + t00
                     , leaf.vpd         = c(szpft$leaf.vpd     )[sel] * 100.
                     , leaf.gsw         = c(szpft$leaf.gsw     )[sel] / day.sec
                     , leaf.par         = c(szpft$par.leaf     )[sel] * 1.0e-6
                     , transp           = c(szpft$transp       )[sel] / day.sec
                     , gpp              = c(szpft$gpp          )[sel]
                     , auto.resp        = c(szpft$plant.resp   )[sel]
                     , mco              = c(szpft$mco          )[sel]
                     , l.transp         = c(szpft$i.transp     )[sel] / day.sec
                     , i.gpp            = c(szpft$i.gpp        )[sel]
                     , i.auto.resp      = c(szpft$i.plant.resp )[sel]
                     , i.mco            = c(szpft$i.mco        )[sel]
                     , i.cba            = c(szpft$cba          )[sel]
                     , growth           = c(szpft$acc.growth   )[sel]
                     , di.mortality     = c(szpft$acc.dimort   )[sel]
                     , ncb.mortality    = c(szpft$acc.ncbmort  )[sel]
                     , recruitment      = c(szpft$acc.recr     )[sel]
                     , cue              = c(szpft$cue          )[sel]
                     , wue              = c(szpft$wue          )[sel]
                     , stringsAsFactors = FALSE
                     )#end data.frame
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #    List of variables to include in patch.                                             #
   #---------------------------------------------------------------------------------------#
   cat0(" + Create data frame with yearly averages by patch.")
   pvars = list( list(vin="agb"          ,vout="agb"          ,add=0. ,mult=1.        )
               , list(vin="lai"          ,vout="lai"          ,add=0. ,mult=1.        )
               , list(vin="ba"           ,vout="ba"           ,add=0. ,mult=1.        )
               , list(vin="fast.soil.c"  ,vout="fast.soil.c"  ,add=0. ,mult=1.        )
               , list(vin="struct.soil.c",vout="struct.soil.c",add=0. ,mult=1.        )
               , list(vin="slow.soil.c"  ,vout="slow.soil.c"  ,add=0. ,mult=1.        )
               , list(vin="can.temp"     ,vout="can.temp"     ,add=t00,mult=1.        )
               , list(vin="can.vpd"      ,vout="can.vpd"      ,add=0. ,mult=100.      )
               , list(vin="nee"          ,vout="nee"          ,add=0. ,mult=1.e-6     )
               , list(vin="hflxca"       ,vout="sensible.flux",add=0. ,mult=1.        )
               , list(vin="wflxca"       ,vout="vapour.flux"  ,add=0. ,mult=1./day.sec)
               , list(vin="qwflxca"      ,vout="latent.flux"  ,add=0. ,mult=1.        )
               , list(vin="cflxca"       ,vout="co2.flux"     ,add=0. ,mult=1.e-6     )
               , list(vin="cflxst"       ,vout="co2.storage"  ,add=0. ,mult=1.e-6     )
               , list(vin="rshortup"     ,vout="rshort.out"   ,add=0. ,mult=1.        )
               , list(vin="rlongup"      ,vout="rlong.out"    ,add=0. ,mult=1.        )
               , list(vin="leaf.temp"    ,vout="leaf.temp"    ,add=t00,mult=1.        )
               , list(vin="leaf.vpd"     ,vout="leaf.vpd"     ,add=0. ,mult=100.      )
               , list(vin="leaf.gsw"     ,vout="leaf.gsw"     ,add=0. ,mult=1./day.sec)
               , list(vin="sm.stress"    ,vout="sm.stress"    ,add=0. ,mult=1.        )
               , list(vin="transp"       ,vout="transp"       ,add=0. ,mult=1./day.sec)
               , list(vin="gpp"          ,vout="gpp"          ,add=0. ,mult=1.        )
               , list(vin="plant.resp"   ,vout="auto.resp"    ,add=0. ,mult=1.        )
               , list(vin="het.resp"     ,vout="het.resp"     ,add=0. ,mult=1.        )
               )#end list
   pvars  = list.2.data.frame(pvars)
   npvars = nrow(pvars)
   #---------------------------------------------------------------------------------------#



   #----- Patch year. ---------------------------------------------------------------------#
   epatch.year  = as.numeric(substring(names(patch$area),2,5))
   ypatch.year  = sort(unique(epatch.year))
   nypatch.year = length(ypatch.year)
   #---------------------------------------------------------------------------------------#

   ypatch = data.frame(matrix(nrow=3*nypatch.year,ncol=2+npvars))
   names(ypatch) = c("year", "ipat",pvars$vout)
   ypatch$year   = rep(ypatch.year,each=3)
   ypatch$ipat   = rep(c("lwr","mid","upr"),times=nypatch.year)
   for (ip in sequence(npvars)){
      #----- Alias for variables. ---------------------------------------------------------#
      vin   = pvars$vin [ip]
      vout  = pvars$vout[ip]
      vadd  = pvars$add [ip]
      vmult = pvars$mult[ip]
      #------------------------------------------------------------------------------------#



      #----- Find lower, median, and upper bounds. ----------------------------------------#
      qlwr  = mapply( FUN = weighted.quantile
                    , x   = patch[[vin]]
                    , w   = patch$area
                    , MoreArgs = list(qu=0.025,na.rm=TRUE)
                    )#end mapply
      qmid  = mapply( FUN = weighted.quantile
                    , x   = patch[[vin]]
                    , w   = patch$area
                    , MoreArgs = list(qu=0.500,na.rm=TRUE)
                    )#end mapply
      qupr  = mapply( FUN = weighted.quantile
                    , x   = patch[[vin]]
                    , w   = patch$area
                    , MoreArgs = list(qu=0.975,na.rm=TRUE)
                    )#end mapply
      qlwr  = tapply(X=qlwr,INDEX=epatch.year,FUN=mean)
      qmid  = tapply(X=qmid,INDEX=epatch.year,FUN=mean)
      qupr  = tapply(X=qupr,INDEX=epatch.year,FUN=mean)
      #------------------------------------------------------------------------------------#


      #----- Append patch. ----------------------------------------------------------------#
      ypatch[[vout]] = vadd + c(t(cbind(qlwr,qmid,qupr))) * vmult
      #------------------------------------------------------------------------------------#
   }#end for
   #---------------------------------------------------------------------------------------#

   #---------------------------------------------------------------------------------------#
   #    Write csv files with the output.                                                   #
   #---------------------------------------------------------------------------------------#
   cat0(" + Write output files.")
   csv.ymean = file.path(outcsv,paste0(place,"_ymean.csv" ))
   csv.ympft = file.path(outcsv,paste0(place,"_ympft.csv" ))
   csv.ymdbh = file.path(outcsv,paste0(place,"_ymdbh.csv" ))
   csv.ympat = file.path(outcsv,paste0(place,"_ympat.csv" ))

   dummy = write.table(x=ymean ,file=csv.ymean,quote=FALSE,sep=",",row.names=FALSE)
   dummy = write.table(x=ypft  ,file=csv.ympft,quote=FALSE,sep=",",row.names=FALSE)
   dummy = write.table(x=ydbh  ,file=csv.ymdbh,quote=FALSE,sep=",",row.names=FALSE)
   dummy = write.table(x=yglob ,file=csv.ympat,quote=FALSE,sep=",",row.names=FALSE)
   #---------------------------------------------------------------------------------------#

}#end for (place in myplaces)
#------------------------------------------------------------------------------------------#
