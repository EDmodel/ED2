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
nyears.long    = 15           # Runs longer than this are considered long runs.
n.density      = 512          # Number of density points
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
paper          = "square"               # Paper size, to define the plot shape
ptsz           = 17                     # Font size.
lwidth         = 2.5                    # Line width
plotgrid       = FALSE                  # Should I plot the grid in the background? 
sasfixlimits   = FALSE                  # Use a fixed scale for size and age-structure
                                        #    plots? (FALSE will set a suitable scale for
                                        #    each plot)
fcgrid         = FALSE                  # Include a grid on the filled contour plots?
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
f.leg          = 1/6                    # Fraction of device for legend
emean.line     = TRUE                   # Use lines instead of points and lines for theme
                                        #    plots? 
                                        #    TRUE  - Lines only
                                        #    FALSE - Use type defined in pmonthly_varlist.r
                                        #    NA    - Let the script decide
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
obsrfile = file.path(srcdir,"LBA_MIP.v9.RData")
load(file=obsrfile)
#------------------------------------------------------------------------------------------#



#----- Define plot window size ------------------------------------------------------------#
f.ext  = f.leg / (1. - f.leg)
sqsize = plotsize(proje=FALSE,paper=paper)
exsize = plotsize(proje=FALSE,paper=paper,extendfc="lon",extfactor=f.ext)
eysize = plotsize(proje=FALSE,paper=paper,extendfc="lat",extfactor=f.ext)
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
   outpref = file.path(outmain,"monthly")
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
   if ((yearend - yearbeg + 1) <= nyears.long){
      sasmonth   = sasmonth.short
      if (is.na(emean.line)) emean.line = FALSE
   }else{
      sasmonth   = sasmonth.long
      if (is.na(emean.line)) emean.line = TRUE
   }#end if
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Find the total number of months that can be loaded this time.                     #
   #---------------------------------------------------------------------------------------#
   ntimes     = (yearz-yeara-1)*12+meszz+(12-monthbeg+1)
   #---------------------------------------------------------------------------------------#



   #----- Print a banner to entretain the user. -------------------------------------------#
   cat0(" + Post-process output from ",lieu,".")
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      Make the RData file name, then we check whether we must read the files again     #
   # or use the stored RData.                                                              #
   #---------------------------------------------------------------------------------------#
   path.data  = file.path(here,place,"2019_rdata_month")
   if (! file.exists(path.data)) dir.create(path.data)
   ed22.rdata  = file.path(path.data,paste0(place,".RData"))
   ed22.status = file.path(path.data,paste0("status_",place,".txt"))
   if (reload.data && file.exists(ed22.rdata)){
      #----- Load the modelled dataset. ---------------------------------------------------#
      cat0("   - Load previous session.")
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
      cat0("   - Start new session.")
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
      cat0(" + Save data to ",basename(ed22.rdata),".")
      save(datum,file=ed22.rdata)
      #------------------------------------------------------------------------------------#
   }#end if (! complete)
   #---------------------------------------------------------------------------------------#


   #----- Update status file with latest data converted into R. ---------------------------#
   latest = paste(datum$year[ntimes],datum$month[ntimes],sep=" ")
   dummy  = write(x=latest,file=ed22.status,append=FALSE)
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
   #     Find the mean and standard deviation.  For variables for which we did not track   #
   # the mean sum of squares, the standard deviation is just the standard deviation of the #
   # means, otherwise we convert the mean sum of squares:                                  #
   #              ____________________          _____________________________________      #
   #             / SUM_i[X_i - Xm]^2           /  / SUM_i[X_i^2]         \      1          #
   # sigma = \  /  ------------------   =  \  /  |  ------------  - Xm^2 | ---------       #
   #          \/       N - 1                \/    \      N               /   1 - 1/N       #
   #                                                                                       #
   # srnonm1 is the square root of 1 / (1 - 1/N)                                           #
   #     Find the standard deviation.                                                      #
   #---------------------------------------------------------------------------------------#
   cat0("    - Find the monthly means.")

   srnorm1 = sqrt(1./(1. - 1. / datum$montable))
   srnorm1[!is.finite(srnorm1)] = 0.

   mmean = list()
   msdev = list()
   for (vname in names(emean)){
      if (vname %in% names(emsqu)){
         has.emsqu = any(is.finite(emsqu[[vname]]))
      }else{
         has.emsqu = FALSE
      }#end if

      #------------------------------------------------------------------------------------#
      #     Soil variables are multi-dimensional.  Use qapply.  Otherwise, check whether   #
      # the mean sum of squares is available or not.                                       #
      #------------------------------------------------------------------------------------#
      if (vname %in% c("soil.temp","soil.water","soil.mstpot","soil.extracted")){
         mmean[[vname]] = qapply(X=emean[[vname]], INDEX=mfac, DIM=1, FUN=mean, na.rm=TRUE)
         msdev[[vname]] = qapply(X=emean[[vname]], INDEX=mfac, DIM=1, FUN=sd  , na.rm=TRUE)
      }else if (has.emsqu){
         mmean[[vname]] = tapply(X=emean[[vname]], INDEX=mfac, FUN=mean, na.rm=TRUE)
         mmsqu          = tapply(X=emsqu[[vname]], INDEX=mfac, FUN=mean, na.rm=TRUE)
         msdev[[vname]] = sqrt  ( mmsqu - mmean[[vname]]^ 2 ) * srnorm1
      }else{
         mmean[[vname]] = tapply(X=emean[[vname]], INDEX=mfac, FUN=mean, na.rm=TRUE)
         msdev[[vname]] = tapply(X=emean[[vname]], INDEX=mfac, FUN=sd  , na.rm=TRUE)
      }#end if
      #------------------------------------------------------------------------------------#


      #----- Fix the bad data. ------------------------------------------------------------#
      bad.mmean = ! is.finite(mmean[[vname]])
      bad.msdev = ! is.finite(msdev[[vname]])
      mmean[[vname]][bad.mmean] = NA
      msdev[[vname]][bad.msdev] = 0.
      #------------------------------------------------------------------------------------#
   }#end for
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #      Here we find the Mean diurnal cycle for each month, then compute the standard    #
   # deviation, similar to the monthly mean.                                               #
   #---------------------------------------------------------------------------------------#
   cat0("    - Aggregate the monthly mean of the diurnal cycle.")
   umean              = list()
   usdev              = list()
   for (vname in names(qmean)){

      #------------------------------------------------------------------------------------#
      #     Soil variables are multi-dimensional.  Use qapply.  Otherwise, check whether   #
      # the mean sum of squares is available or not.                                       #
      #------------------------------------------------------------------------------------#
      if (vname %in% names(qmsqu)){
         umean[[vname]] = qapply(X=qmean[[vname]], INDEX=mfac, DIM=1,FUN=mean, na.rm=TRUE)
         umsqu          = qapply(X=qmsqu[[vname]], INDEX=mfac, DIM=1,FUN=mean, na.rm=TRUE)
         usdev[[vname]] = sqrt  ( umsqu - umean[[vname]]^ 2 ) * srnorm1
      }else{
         umean[[vname]] = qapply(X=qmean[[vname]], INDEX=mfac, DIM=1,FUN=mean, na.rm=TRUE)
         usdev[[vname]] = qapply(X=qmean[[vname]], INDEX=mfac, DIM=1,FUN=sd  , na.rm=TRUE)
      }#end if
      #------------------------------------------------------------------------------------#


      #----- Fix the bad data. ------------------------------------------------------------#
      bad.umean = ! is.finite(umean[[vname]])
      bad.usdev = ! is.finite(usdev[[vname]])
      umean[[vname]][bad.umean] = NA
      usdev[[vname]][bad.usdev] = 0.
      #------------------------------------------------------------------------------------#
   }#end for
   #---------------------------------------------------------------------------------------#





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
   empty = is.na(szpft$nplant) | szpft$nplant == 0
   for (vname in names(szpft)) szpft[[vname]][empty] = NA
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Replace the mortality and recruitment exponential rates by the "interests" rates. #
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



   #---------------------------------------------------------------------------------------#
   #      Find the patch density function for all patch-level data.                        #
   #---------------------------------------------------------------------------------------#
   cat0("    - Find the distribution function of patch properties.")
   patchpdf = list()
   for (pp in sequence(nplotpatch)){
      #----- Get settings for the current variable. ---------------------------------------#
      this        = plotpatch[[pp]]
      vname       = this$vnam
      col.scheme  = get(this$col.scheme)(n=ncolsfc)
      #------------------------------------------------------------------------------------#



      #----- Retrieve patch values, and find monthly means. -------------------------------#
      emean.area  = patch$area
      emean.vname = patch[[vname]]
      #------------------------------------------------------------------------------------#



      #----- Find the range for which we find the density function. -----------------------#
      espan.vname = pretty.xylim(u=c(unlist(emean.vname)),fracexp=c(-0.04,0.04))
      elwr.vname  = espan.vname[1]
      eupr.vname  = espan.vname[2]
      #------------------------------------------------------------------------------------#



      #----- Find the density function for each time. -------------------------------------#
      edfun.now   = mapply( FUN      = density.safe
                          , x        = emean.vname
                          , weights  = emean.area
                          , MoreArgs = list(n=n.density,from=elwr.vname,to=eupr.vname)
                          )#end mapply
      #------------------------------------------------------------------------------------#




      #----- Save the density function. ---------------------------------------------------#
      edfun        = list()
      edfun$x      = chron(datum$when)
      edfun$y      = seq(from=elwr.vname,to=eupr.vname,length.out=n.density)
      edfun$z      = t(sapply(X=edfun.now["y",],FUN=cbind))
      #------------------------------------------------------------------------------------#




      #----- Save the density function. ---------------------------------------------------#
      mdfun        = list()
      mdfun$x      = sort(unique(nummonths(edfun$x)))
      mdfun$y      = edfun$y
      mdfun$z      = qapply(X=edfun$z,DIM=1,INDEX=nummonths(edfun$x),FUN=mean,na.rm=TRUE)
      #------------------------------------------------------------------------------------#



      #----- Remove tiny values (even with log scale values can be very hard to see. ------#
      bye          = ! ( edfun$z %>=% 1.e-4 * max(unlist(edfun$z),na.rm=TRUE) )
      edfun$z[bye] = NA_real_
      #------------------------------------------------------------------------------------#


      #----- Remove tiny values (even with log scale values can be very hard to see. ------#
      bye          = ! ( mdfun$z %>=% 1.e-4 * max(unlist(mdfun$z),na.rm=TRUE) )
      mdfun$z[bye] = NA_real_
      #------------------------------------------------------------------------------------#
      patchpdf[[vname]] = list(edensity=edfun,mdensity=mdfun)
   }#end for (pp in sequence(nplotpatch))
   #---------------------------------------------------------------------------------------#




   #----- Find which PFTs, land uses and transitions we need to consider ------------------#
   pftave    = apply( X      = szpft$agb[,ndbh+1,]
                    , MARGIN = 2
                    , FUN    = mean
                    , na.rm  = TRUE
                    )#end apply
   luave     = apply( X      = lu$agb 
                    , MARGIN = 2
                    , FUN    = mean
                    , na.rm  = TRUE
                    )#end apply
   distave   = apply(X=lu$dist,MARGIN=c(2,3),FUN=mean)
   selpft    = pftave  %>% 0.
   sellu     = luave   %>% 0.
   seldist   = distave %>% 0.
   n.selpft  = sum(selpft )
   n.sellu   = sum(sellu  )
   n.seldist = sum(seldist)
   #---------------------------------------------------------------------------------------#





   #=======================================================================================#
   #=======================================================================================#
   #=======================================================================================#
   #      Plotting section begins here...                                                  #
   #---------------------------------------------------------------------------------------#
   cat0("    - Plot figures.")
   #---------------------------------------------------------------------------------------#






   #---------------------------------------------------------------------------------------#
   #      Time series by PFT.                                                              #
   #---------------------------------------------------------------------------------------#
   for (v in sequence(ntspftdbh)){
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
         outdir = file.path(outpref,"tspft")
         if (! file.exists(outdir)) dir.create(outdir)
         cat0("      + ",description," time series for all PFTs.")

         #----- Load variable -------------------------------------------------------------#
         if (vnam %in% names(szpft)){
            thisvar = szpft[[vnam]][,ndbh+1,]
            if (plog){
               #----- Eliminate non-positive values in case it is a log plot. -------------#
               badlog          = ! (thisvar %>% 0)
               thisvar[badlog] = NA
            }#end if
         }else{
            thisvar = matrix(NA,ncol=npft+1,nrow=ntimes)
         }#end if
         #---------------------------------------------------------------------------------#



         #----- Loop over output formats. -------------------------------------------------#
         for (o in sequence(nout)){
            #----- Open file. -------------------------------------------------------------#
            fichier = file.path(outdir,paste0(vnam,"-",suffix,".",outform[o]))
            dummy   = open.plot( fichier = fichier
                               , outform = outform[o]
                               , size    = eysize
                               , ptsz    = ptsz
                               , depth   = depth
                               )#end open.plot
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #     Find the limit, make some room for the legend, and in case the field is  #
            # a constant, nudge the limits so the plot command will not complain.          #
            #------------------------------------------------------------------------------#
            xlimit = pretty.xylim(u=as.numeric(datum$tomonth),fracexp=0.0,is.log=FALSE)
            ylimit = pretty.xylim(u=thisvar[,selpft]         ,fracexp=0.0,is.log=plog )
            if (plog){
               xylog    = "y"
               ydrought = c( exp(sqrt(ylimit[1]^3/ylimit[2]))
                           , exp(sqrt(ylimit[2]^3/ylimit[1]))
                           )#end c
            }else{
               xylog    = ""
               ydrought = c(ylimit[1] - 0.5 * diff(ylimit), ylimit[2] + 0.5 * diff(ylimit))
            }#end if
            #------------------------------------------------------------------------------#


            #----- Plot settings. ---------------------------------------------------------#
            letitre       = paste(description,lieu,sep=" - ")
            ley           = desc.unit(desc=description,unit=unit)
            cols          = pft$colour[selpft]
            legs          = pft$name  [selpft]
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #     Split the plot into two windows.                                         #
            #------------------------------------------------------------------------------#
            par(par.user)
            layout(mat=rbind(2,1),heights=c(1-f.leg,f.leg))
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #      First plot: legend.                                                     #
            #------------------------------------------------------------------------------#
            par(mar=c(0.1,4.6,0.1,2.1))
            plot.new()
            plot.window(xlim=c(0,1),ylim=c(0,1))
            legend( x      = "bottom"
                  , inset  = 0.0
                  , legend = legs
                  , col    = cols
                  , lwd    = lwidth
                  , ncol   = min(pretty.box(n.selpft)$ncol,3)
                  , title  = expression(bold("Plant Functional Type"))
                  , xpd    = TRUE
                  , bty    = "n"
                  )#end legend
            #------------------------------------------------------------------------------#




            #------------------------------------------------------------------------------#
            #      Main plot.                                                              #
            #------------------------------------------------------------------------------#
            par(mar=c(4.1,4.6,4.1,2.1))
            plot.new()
            plot.window(xlim=xlimit,ylim=ylimit,log=xylog)
            axis(side=1,at=whenplot8$levels,labels=whenplot8$labels,padj=whenplot8$padj)
            axis(side=2,las=1)
            box()
            title(main=letitre,xlab="Year",ylab=ley,cex.main=0.7)
            if (drought.mark){
               for (n in sequence(ndrought)){
                  rect(xleft  = drought[[n]][1],ybottom = ydrought[1]
                      ,xright = drought[[n]][2],ytop    = ydrought[2]
                      ,col    = grid.colour,border=NA)
               }#end for
            }#end if
            #----- Plot grid. -------------------------------------------------------------#
            if (plotgrid){ 
               abline(v=whenplot8$levels,h=axTicks(side=2),col=grid.colour,lty="solid")
            }#end if
            #----- Plot lines. ------------------------------------------------------------#
            for (n in 1:(npft+1)){
               if (selpft[n]){
                  lines(datum$tomonth,thisvar[,n],type="l",col=pft$colour[n],lwd=lwidth)
               }#end if
            }#end for
            #------------------------------------------------------------------------------#


            #----- Close the device. ------------------------------------------------------#
            dummy = close.plot(outform=outform[o])
            #------------------------------------------------------------------------------#
         } #end for outform
      }#end if (tseragbpft)
   } #end for tseries
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #      Time series by DBH, by PFT.                                                      #
   #---------------------------------------------------------------------------------------#
   #----- Find the PFTs to plot. ----------------------------------------------------------#
   pftmat  = apply( X      = szpft$nplant[,-(ndbh+1),-(npft+1),drop=FALSE]
                  , MARGIN = c(1,3)
                  , FUN    = function(x) sum(is.finite(x))
                  )#end apply
   pftuse  = which(apply(X=pftmat,MARGIN=2,FUN=function(x) max(x) > 1))
   for (v in sequence(ntspftdbh)){
      thistspftdbh   = tspftdbh[[v]]
      vnam           = thistspftdbh$vnam
      description    = thistspftdbh$desc
      unit           = thistspftdbh$e.unit
      plog           = thistspftdbh$plog
      plotit         = thistspftdbh$pftdbh
      
      #----- Load variable ----------------------------------------------------------------#
      if (vnam %in% names(szpft)){
         thisvar = szpft[[vnam]]
         if (plog){
            xylog           = "y"
            badlog          = ! (thisvar %>% 0)
            thisvar[badlog] = NA
         }else{
            xylog           = ""
         }#end if
      }else{
         thisvar = array(NA,dim=c(ntimes,ndbh+1,npft+1))
      }#end if
      #----- Check whether the user wants to have this variable plotted. ------------------#
      if (plotit && length(pftuse) > 0 && any(is.finite(thisvar))){
         #------ Append "all" to the PFT classes to plot. ---------------------------------#
         pftplt = c(pftuse,npft+1)
         #---------------------------------------------------------------------------------#

         #---------------------------------------------------------------------------------#
         #    Check whether the time series directory exists.  If not, create it.          #
         #---------------------------------------------------------------------------------#
         outdir = file.path(outpref,"tsdbh")
         if (! file.exists(outdir)) dir.create(outdir)
         outvar = file.path(outdir,vnam)
         if (! file.exists(outvar)) dir.create(outvar)
         #---------------------------------------------------------------------------------#

         cat0("      + ",description," time series for DBH class.")


         #---------------------------------------------------------------------------------#
         #     Find the limit, make some room for the legend, and in case the field is a   #
         # constant, nudge the limits so the plot command will not complain.               #
         #---------------------------------------------------------------------------------#
         xlimit     = pretty.xylim(u=datum$tomonth    ,fracexp=0.0,is.log=FALSE)
         ylimit.pft = pretty.xylim(u=thisvar[,,pftuse],fracexp=0.0,is.log=plog)
         ylimit.all = pretty.xylim(u=thisvar[,,npft+1],fracexp=0.0,is.log=plog)
         #---------------------------------------------------------------------------------#




         #---------------------------------------------------------------------------------#
         #       Loop over plant functional types.                                         #
         #---------------------------------------------------------------------------------#
         for (p in pftplt){
            if (p == (npft + 1)){
               pftlab = "pft-00"
               ylimit = ylimit.all
            }else{
               pftlab = paste0("pft-",sprintf("%2.2i",p))
               ylimit = ylimit.pft
            }#end if
            cat0("        - ",pft$name[p],".")


            #------------------------------------------------------------------------------#
            #     Set y-limits for drought rectangle.                                      #
            #------------------------------------------------------------------------------#
            if (plog){
               xylog    = "y"
               ydrought = c( exp(sqrt(ylimit[1]^3/ylimit[2]))
                           , exp(sqrt(ylimit[2]^3/ylimit[1]))
                           )#end c
            }else{
               xylog    = ""
               ydrought = c(ylimit[1] - 0.5 * diff(ylimit), ylimit[2] + 0.5 * diff(ylimit))
            }#end if
            #------------------------------------------------------------------------------#


            #----- Loop over output formats. ----------------------------------------------#
            for (o in sequence(nout)){
               #----- Open file. ----------------------------------------------------------#
               fichier = file.path(outvar,paste0(vnam,"-",pftlab,"-",suffix,".",outform[o]))
               dummy   = open.plot( fichier = fichier
                                 , outform = outform[o]
                                 , size    = eysize
                                 , ptsz    = ptsz
                                 , depth   = depth
                                 )#end open.plot
               #---------------------------------------------------------------------------#



               #-----  Plot annotation. ---------------------------------------------------#
               letitre = paste(description,pft$name[p],lieu,sep=" - ")
               ley     = desc.unit(desc=description,unit=unit)
               #---------------------------------------------------------------------------#


               #---------------------------------------------------------------------------#
               #     Split the plot into two windows.                                      #
               #---------------------------------------------------------------------------#
               par(par.user)
               layout(mat=rbind(2,1),heights=c(1.-f.leg,f.leg))
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #      First plot: legend.                                                  #
               #---------------------------------------------------------------------------#
               par(mar=c(0.1,4.6,0.1,2.1))
               plot.new()
               plot.window(xlim=c(0,1),ylim=c(0,1))
               legend( x      = "bottom"
                     , inset  = 0.0
                     , bg     = background
                     , legend = dbhnames
                     , col    = dbhcols
                     , ncol   = min(pretty.box(ndbh+1)$ncol,3)
                     , title  = expression(bold("DBH class"))
                     , lwd    = lwidth
                     , bty    = "n"
                     , xpd    = TRUE
                     )#end legend
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #      Main plot.                                                           #
               #---------------------------------------------------------------------------#
               par(mar=c(4.1,4.6,4.1,2.1))
               plot.new()
               plot.window(xlim=xlimit,ylim=ylimit,log=xylog)
               axis(side=1,at=whenplot8$levels,labels=whenplot8$labels,padj=whenplot8$padj)
               axis(side=2,las=1)
               box()
               title(main=letitre,xlab="Year",ylab=ley,cex.main=0.7)
               if (drought.mark){
                  for (n in sequence(ndrought)){
                     rect(xleft  = drought[[n]][1],ybottom = ydrought[1]
                         ,xright = drought[[n]][2],ytop    = ydrought[2]
                         ,col    = grid.colour,border=NA)
                  }#end for
               }#end if
               #----- Plot grid. ----------------------------------------------------------#
               if (plotgrid){ 
                  abline(v=whenplot8$levels,h=axTicks(side=2),col=grid.colour,lty="solid")
               }#end if
               #----- Plot lines. ---------------------------------------------------------#
               for (d in seq(from=1,to=ndbh+1,by=1)){
                  lines(datum$tomonth,thisvar[,d,p],type="l",col=dbhcols[d],lwd=lwidth)
               }#end for
               #---------------------------------------------------------------------------#


               #----- Close the device. ---------------------------------------------------#
               dummy = close.plot(outform=outform[o])
               #---------------------------------------------------------------------------#
            }#end for outform
            #------------------------------------------------------------------------------#
         }#end for (p in pftuse)
         #---------------------------------------------------------------------------------#
      }#end if (tseragbpft)
      #------------------------------------------------------------------------------------#
   } #end for tseries
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #   Plot the comparison between observations and model.                                 #
   #---------------------------------------------------------------------------------------#
   cat0("    + Comparisons of time series (model vs. observations).")
   for (cc in sequence(ncompmodel)){

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
      if (iata %in% c("mao","bdf")){
         obsnow = "obs.m34"
      }else if(iata %in% c("stm")){
         obsnow = "obs.s67"
      }else if(iata %in% c("rao")){
         obsnow = "obs.pdg"
      }else if(iata %in% c("jpr")){
         obsnow = "obs.fns"
      }else if(iata %in% c("btr")){
         obsnow = "obs.s77"
      }else{
         obsnow = paste0("obs.",iata)
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
         outdir   = file.path(outpref,"compemean")
         if (! file.exists(outdir)) dir.create(outdir)
         cat0("      - ",description," comparison.")
         #---------------------------------------------------------------------------------#



         #----- Define the number of layers. ----------------------------------------------#
         thiswhen  = datum$tomonth [sel]
         thismean  = emean[[vname]][sel]
         #---------------------------------------------------------------------------------# 



         #----- Find the plot range. ------------------------------------------------------#
         xlimit = range(thiswhen)
         ylimit = pretty.xylim(u=c(thismean,obsmean),fracexp=0.0,is.log=FALSE)
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Find the nice scale for time.                                               #
         #---------------------------------------------------------------------------------#
         whenplote = pretty.time(obswhen,n=8)
         #---------------------------------------------------------------------------------#



         #----- Plot annotation. ----------------------------------------------------------#
         letitre = paste0(description," - ",lieu,"\n","Monthly mean")
         ley     = desc.unit(desc=description,unit=unit)
         #---------------------------------------------------------------------------------#


         #----- Loop over formats. --------------------------------------------------------#
         for (o in sequence(nout)){
            fichier = file.path(outdir,paste0(vname,".",outform[o]))
            dummy   = open.plot( fichier = fichier
                               , outform = outform[o]
                               , size    = eysize
                               , ptsz    = ptsz
                               , depth   = depth
                               )#end open.plot

            #----- Split window into two. -------------------------------------------------#
            par(par.user)
            layout(mat=rbind(2,1),heights=c(1.-f.leg,f.leg))
            #------------------------------------------------------------------------------#


            #----- First plot: the legend. ------------------------------------------------#
            par(mar=c(0.1,4.6,0.1,2.1))
            plot.new()
            plot.window(xlim=c(0,1),ylim=c(0,1))
            legend(x="bottom",inset=0.0,legend=c("Model","Observation")
                  ,col=lcolours,lwd=llwd,cex=cex.ptsz,pch=16,xpd=TRUE,bty="n")
            #------------------------------------------------------------------------------#



            #----- Second panel: the actual plot. -----------------------------------------#
            par(mar=c(4.1,4.6,4.1,2.1))
            plot.new()
            plot.window(xlim=xlimit,ylim=ylimit,log=plog)
            if (plotgrid){
               abline(v=whenplote$levels,h=axTicks(side=2),col=grid.colour,lty="solid")
            }#end if
            points(x=thiswhen,y=thismean,col=lcolours[1],lwd=llwd[1],type=ltype
                  ,pch=16,cex=1.0)
            points(x=obswhen,y=obsmean ,col=lcolours[2],lwd=llwd[2],type=ltype
                  ,pch=16,cex=1.0)
            box()
            axis(side=1,at=whenplote$levels,labels=whenplote$labels,padj=whenplote$padj)
            axis(side=2,las=1)
            title(main=letitre,xlab="Time",ylab=ley,cex.main=0.8)
            #------------------------------------------------------------------------------#



            #----- Close plot. ------------------------------------------------------------#
            dummy = close.plot(outform=outform[o])
            #------------------------------------------------------------------------------#
         }#end for outform
         #---------------------------------------------------------------------------------#
      }#end if plotit
      #------------------------------------------------------------------------------------#
   }#end for ncompare
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #   Plot the comparison between observations and model.                                 #
   #---------------------------------------------------------------------------------------#
   cat0("    + Comparisons of monthly means (model vs. observations).")
   for (cc in sequence(ncompmodel)){

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
      if (iata %in% c("mao","bdf")){
         obsnow = "obs.m34"
      }else if(iata %in% c("stm")){
         obsnow = "obs.s67"
      }else if(iata %in% c("rao")){
         obsnow = "obs.pdg"
      }else if(iata %in% c("jpr")){
         obsnow = "obs.fns"
      }else if(iata %in% c("btr")){
         obsnow = "obs.s77"
      }else{
         obsnow = paste0("obs.",iata)
      }#end if
      #------------------------------------------------------------------------------------#

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
         outdir   = file.path(outpref,"compmmean")
         if (! file.exists(outdir)) dir.create(outdir)
         cat0("      - ",description," comparison.")
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
         xlimit    = range(montmont)
         if (plotsd){
            ylimit = c(mod.ylow,mod.yhigh,obs.ylow,obs.yhigh)
         }else{
            ylimit = c(thismean,obsmean)
         }#end if
         ylimit    = pretty.xylim(u=ylimit,fracexp=0.0,is.log=FALSE)
         #---------------------------------------------------------------------------------#



         #------  Plot annotation. --------------------------------------------------------#
         letitre = paste0(description," - ",lieu,"\n","Monthly mean")
         ley     = desc.unit(desc=description,unit=unit)
         #---------------------------------------------------------------------------------#




         #----- Loop over formats. --------------------------------------------------------#
         for (o in sequence(nout)){
            fichier = file.path(outdir,paste0(vname,".",outform[o]))
            dummy   = open.plot( fichier = fichier
                               , outform = outform[o]
                               , size    = eysize
                               , ptsz    = ptsz
                               , depth   = depth
                               )#end open.plot


            #----- Split window into two. -------------------------------------------------#
            par(par.user)
            layout(mat=rbind(2,1),heights=c(1.-f.leg,f.leg))
            #------------------------------------------------------------------------------#


            #----- First plot: the legend. ------------------------------------------------#
            par(mar=c(0.1,4.6,0.1,2.1))
            plot.new()
            plot.window(xlim=c(0,1),ylim=c(0,1))
            if (plotsd){
               legend(x="bottom",inset=0.0,legend=c("Model","Observation")
                     ,fill=errcolours,angle=angle,density=dens,lwd=llwd,col=lcolours
                     ,bg=background,title="Shaded areas = 1 SD",cex=cex.ptsz
                     ,xpd=TRUE,pch=16,bty="n")
            }else{
               legend(x="bottom",inset=0.0,legend=c("Model","Observation")
                     ,col=lcolours,lwd=llwd,cex=cex.ptsz,xpd=TRUE,pch=16,bty="n")
            }#end if
            #------------------------------------------------------------------------------#



            #----- Second panel: the actual plot. -----------------------------------------#
            par(mar=c(4.1,4.6,4.1,2.1))
            plot.new()
            plot.window(xlim=xlimit,ylim=ylimit,log=plog)
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
            box()
            title(main=letitre,xlab="Time",ylab=ley,cex.main=0.7)
            axis(side=1,at=mplot$levels,labels=mplot$labels,padj=mplot$padj)
            axis(side=2,las=1)


            dummy = close.plot(outform=outform[o])
         } #end for outform
      }#end if plotit
   }#end for ncompare
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #   Plot the comparison between observations and model.                                 #
   #---------------------------------------------------------------------------------------#
   cat0("      * Comparisons of mean diurnal cycle (model vs. observations).")
   for (cc in sequence(ncompmodel)){

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
      if (iata %in% c("mao","bdf")){
         obsnow = "obs.m34"
      }else if(iata %in% c("stm")){
         obsnow = "obs.s67"
      }else if(iata %in% c("rao")){
         obsnow = "obs.pdg"
      }else if(iata %in% c("jpr")){
         obsnow = "obs.fns"
      }else if(iata %in% c("btr")){
         obsnow = "obs.s77"
      }else{
         obsnow = paste0("obs.",iata)
      }#end if
      plotit       = plotit && obsnow %in% ls()
      #------------------------------------------------------------------------------------#

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
         outdir   = file.path(outpref,"compdcyc")
         if (! file.exists(outdir)) dir.create(outdir)
         outtheme = file.path(outdir,vname)
         if (! file.exists(outtheme)) dir.create(outtheme)
         cat0("      + ",description," comparison.")
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
         xlimit    = range(thisday)
         if (plotsd){
            ylimit = c(mod.ylow,mod.yhigh,obs.ylow,obs.yhigh)
         }else{
            ylimit = c(thismean,obsmean)
         }#end if
         ylimit    = pretty.xylim(u=ylimit,fracexp=0.0,is.log=FALSE)
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #      Loop over all months.                                                      #
         #---------------------------------------------------------------------------------#
         for (pmon in sequence(12)){
            cmon    = sprintf("%2.2i",pmon)
            namemon = mlist[pmon]

            #------------------------------------------------------------------------------#
            #     Check if the directory exists.  If not, create it.                       #
            #------------------------------------------------------------------------------#
            cat0("        > ",description," time series - ",namemon,".")


            #------ Plot annotation. ------------------------------------------------------#
            letitre = paste0(description," - ",lieu,"\n","Mean diurnal cycle - ",namemon)
            ley     = desc.unit(desc=description,unit=unit)
            #------------------------------------------------------------------------------#


            #----- Loop over formats. -----------------------------------------------------#
            for (o in sequence(nout)){
               fichier = file.path(outtheme,paste0(vname,"-",cmon,".",outform[o]))
               dummy   = open.plot( fichier = fichier
                                  , outform = outform[o]
                                  , size    = eysize
                                  , ptsz    = ptsz
                                  , depth   = depth
                                  )#end open.plot



               #----- Split window into two. ----------------------------------------------#
               par(par.user)
               layout(mat=rbind(2,1),heights=c(1.-f.leg,f.leg))
               #---------------------------------------------------------------------------#


               #----- First plot: the legend. ---------------------------------------------#
               par(mar=c(0.1,4.6,0.1,2.1))
               plot.new()
               plot.window(xlim=c(0,1),ylim=c(0,1))
               if (plotsd){
                  legend(x="bottom",inset=0.0,legend=c("Model","Observation")
                        ,fill=errcolours,angle=angle,density=dens,lwd=llwd,col=lcolours
                        ,bg=background,title="Shaded areas = 1 SD",cex=cex.ptsz,xpd=TRUE
                        ,pch=16,bty="n")
               }else{
                  legend(x="bottom",inset=0.0,legend=c("Model","Observation")
                        ,col=lcolours,lwd=llwd,cex=cex.ptsz,xpd=TRUE,pch=16,bty="n")
               }#end if
               #---------------------------------------------------------------------------#



               #----- Second panel: the actual plot. --------------------------------------#
               par(mar=c(4.1,4.6,4.1,2.1))
               plot.new()
               plot.window(xlim=xlimit,ylim=ylimit,log=plog)
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
               box()
               title(main=letitre,xlab="Time",ylab=ley,cex.main=0.7)
               axis(side=1,at=uplot$levels,labels=uplot$labels,padj=uplot$padj)
               axis(side=2,las=1)
               #---------------------------------------------------------------------------#



               dummy = close.plot(outform=outform[o])
               #---------------------------------------------------------------------------#
            } #end for outform
         }#end for pmon
      }#end if plotit
   }#end for ncompare
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #      Time series by LU.                                                               #
   #---------------------------------------------------------------------------------------#
   for (v in sequence(ntslu)){
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
         outdir = file.path(outpref,"tslu")
         if (! file.exists(outdir)) dir.create(outdir)
         cat0("      + ",description," time series for all LUs.")



         #----- Load variable -------------------------------------------------------------#
         if (vnam %in% names(lu)){
            thisvar = lu[[vnam]]
            if (plog){
               xylog           = "y"
               badlog          = ! (thisvar %>% 0)
               thisvar[badlog] = NA
            }else{
               xylog           = ""
            }#end if
         }else{
            thisvar = array(NA,dim=c(ntimes,nlu+1))
         }#end if
         #---------------------------------------------------------------------------------#



         #----- Loop over output formats. -------------------------------------------------#
         for (o in sequence(nout)){
            fichier = file.path(outdir,paste0(vnam,"-",suffix,".",outform[o]))
            dummy   = open.plot( fichier = fichier
                               , outform = outform[o]
                               , size    = eysize
                               , ptsz    = ptsz
                               , depth   = depth
                               )#end open.plot


            #------------------------------------------------------------------------------#
            #     Find the limit, make some room for the legend, and in case the field is  #
            # a constant, nudge the limits so the plot command will not complain.          #
            #------------------------------------------------------------------------------#
            xlimit = pretty.xylim(u=as.numeric(datum$tomonth),fracexp=0.0,is.log=FALSE)
            ylimit = pretty.xylim(u=thisvar[,sellu]          ,fracexp=0.0,is.log=plog )
            if (plog){
               xylog    = "y"
               ydrought = c( exp(sqrt(ylimit[1]^3/ylimit[2]))
                           , exp(sqrt(ylimit[2]^3/ylimit[1]))
                           )#end c
            }else{
               xylog    = ""
               ydrought = c(ylimit[1] - 0.5 * diff(ylimit), ylimit[2] + 0.5 * diff(ylimit))
            }#end if
            #------------------------------------------------------------------------------#



            #----- Plot settings. ---------------------------------------------------------#
            letitre = paste(description,lieu,sep=" - ")
            ley     = desc.unit(desc=description,unit=unit)
            cols    = lucols[sellu]
            legs    = lunames[sellu]
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #     Split the plot into two windows.                                         #
            #------------------------------------------------------------------------------#
            par(par.user)
            layout(mat=rbind(2,1),heights=c(1.-f.leg,f.leg))
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #      First plot: legend.                                                     #
            #------------------------------------------------------------------------------#
            par(mar=c(0.1,4.6,0.1,2.1))
            plot.new()
            plot.window(xlim=c(0,1),ylim=c(0,1))
            legend( x      = "bottom"
                  , inset  = 0.0
                  , legend = legs
                  , col    = cols
                  , lwd    = lwidth
                  , ncol   = min(3,pretty.box(n.sellu)$ncol)
                  , title  = expression(bold("Land use type"))
                  , xpd    = TRUE
                  , bty    = "n"
                  )#end legend
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #      Main plot.                                                              #
            #------------------------------------------------------------------------------#
            par(mar=c(4.1,4.6,4.1,2.1))
            plot.new()
            plot.window(xlim=xlimit,ylim=ylimit,log=xylog)
            axis(side=1,at=whenplot8$levels,labels=whenplot8$labels,padj=whenplot8$padj)
            axis(side=2,las=1)
            box()
            title(main=letitre,xlab="Year",ylab=ley,cex.main=0.7,log=xylog)
            if (drought.mark){
               for (n in sequence(ndrought)){
                  rect(xleft  = drought[[n]][1],ybottom = ydrought[1]
                      ,xright = drought[[n]][2],ytop    = ydrought[2]
                      ,col    = grid.colour,border=NA)
               }#end for
            }#end if
            #----- Plot grid. -------------------------------------------------------------#
            if (plotgrid){ 
               abline(v=whenplot8$levels,h=axTicks(side=2),col=grid.colour,lty="solid")
            }#end if
            #----- Plot lines. ------------------------------------------------------------#
            for (n in sequence(nlu+1)){
               if (sellu[n]){
                  lines(datum$tomonth,thisvar[,n],type="l",col=lucols[n],lwd=lwidth)
               }#end if
            }#end for
            #------------------------------------------------------------------------------#


            #----- Close the device. ------------------------------------------------------#
            dummy = close.plot(outform=outform[o])
            #------------------------------------------------------------------------------#
         }#end for outform
         #---------------------------------------------------------------------------------#
      }#end if (tseragbpft)
      #------------------------------------------------------------------------------------#
   }#end for tseries
   #---------------------------------------------------------------------------------------#





   #---------------------------------------------------------------------------------------#
   #   Plot disturbance rate by disturbance transition.                                    #
   #---------------------------------------------------------------------------------------#
   if (tserdist && any(seldist)){
      cat0("      + Disturbance rate time series for all disturbances.")
      for (o in sequence(nout)){
         fichier = file.path(outpref,paste0("disturb-",suffix,".",outform[o]))
         dummy   = open.plot( fichier = fichier
                            , outform = outform[o]
                            , size    = eysize
                            , ptsz    = ptsz
                            , depth   = depth
                            )#end open.plot

         #---------------------------------------------------------------------------------#
         #     Find the limit, make some room for the legend, and in case the field is a   #
         #  constant, nudge the limits so the plot command will not complain.              #
         #---------------------------------------------------------------------------------#
         xlimit   = pretty.xylim(u=as.numeric(datum$tomonth),fracexp=0.0,is.log=FALSE)
         ylimit   = NULL
         mylucols = NULL
         mylulegs = NULL
         n       = 0
         for (jlu in sequence(nlu)){
            for (ilu in sequence(nlu)){
               n = n + 1
               if (seldist[ilu,jlu]){
                  ylimit   = c(ylimit,lu$dist[,ilu,jlu])
                  mylucols = c(mylucols,distcols [n])
                  mylulegs = c(mylulegs,distnames[n])
               }#end if
            }#end for
         }#end for
         ylimit   = pretty.xylim(u=ylimit,fracexp=0.0,is.log=FALSE)
         ydrought = c(ylimit[1] - 0.5 * diff(ylimit), ylimit[2] + 0.5 * diff(ylimit))
         #---------------------------------------------------------------------------------#



         #----- Plot settings. ------------------------------------------------------------#
         letitre = paste("Disturbance rates",lieu,sep=" - ")
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #     Split the plot into two windows.                                            #
         #---------------------------------------------------------------------------------#
         par(par.user)
         layout(mat=rbind(2,1),heights=c(1.-f.leg,f.leg))
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      First plot: legend.                                                        #
         #---------------------------------------------------------------------------------#
         par(mar=c(0.1,4.6,0.1,2.1))
         plot.new()
         plot.window(xlim=c(0,1),ylim=c(0,1))
         legend( x      = "bottom"
               , inset  = 0.0
               , bg     = background
               , legend = mylulegs
               , col    = mylucols
               , lwd    = lwidth
               , ncol   = min(3,pretty.box(n)$ncol)
               , title  = expression(bold("Transition"))
               , xpd    = TRUE
               , bty    = "n"
               )#end legend
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Main plot.                                                                 #
         #---------------------------------------------------------------------------------#
         par(mar=c(4.1,4.6,4.1,2.1))
         plot.new()
         plot.window(xlim=xlimit,ylim=ylimit,log=xylog)
         axis(side=1,at=whenplot8$levels,labels=whenplot8$labels,padj=whenplot8$padj)
         axis(side=2,las=1)
         box()
         title( main     = letitre
              , xlab     = "Year"
              , ylab     = desc.unit(desc="Disturbance rate",unit=untab$oneoyr)
              , cex.main = 0.7
              )#end title
         if (drought.mark){
            for (n in sequence(ndrought)){
               rect(xleft  = drought[[n]][1],ybottom = ydrought[1]
                   ,xright = drought[[n]][2],ytop    = ydrought[2]
                   ,col    = grid.colour,border=NA)
            }#end for
         }#end if
         #----- Plot grid. ----------------------------------------------------------------#
         if (plotgrid){ 
            abline(v=whenplot8$levels,h=axTicks(side=2),col=grid.colour,lty="solid")
         }#end if
         #----- Plot lines. ---------------------------------------------------------------#
         n = 0
         for (jlu in sequence(nlu)){
            for (ilu in sequence(nlu)){
               n = n + 1
               if (seldist[ilu,jlu]){
                  lines(datum$tomonth,lu$dist[,ilu,jlu],type="l"
                       ,col=distcols[n],lwd=lwidth)
               }#end if
            }#end for
         }#end for
         #---------------------------------------------------------------------------------#


         #----- Close the device. ---------------------------------------------------------#
         dummy = close.plot(outform=outform[o])
         #---------------------------------------------------------------------------------#
      }#end for outform
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #   Plot the time series diagrams showing months and years.                             #
   #---------------------------------------------------------------------------------------#
   cat0("      * Plot some time series with groups of variables.")
   for (hh in sequence(ntheme)){

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
      ylimit.fix   = themenow$emean.lim
   
      if (plotit){

         #---------------------------------------------------------------------------------#
         #    Check whether the time series directory exists.  If not, create it.          #
         #---------------------------------------------------------------------------------#
         outdir = file.path(outpref,"theme_emean")
         if (! file.exists(outdir)) dir.create(outdir)
         cat0("      + ",group," time series for several variables.")


         #----- Define the number of layers. ----------------------------------------------#
         nlayers   = length(vnames)
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Find the limit, make some room for the legend, and in case the field is a   #
         # constant, nudge the limits so the plot command will not complain.               #
         #---------------------------------------------------------------------------------#
         xlimit   = pretty.xylim(u=as.numeric(datum$tomonth),fracexp=0.0,is.log=FALSE)
         if (any(! is.finite(ylimit.fix))){
            ylimit    = NULL
            for (l in 1:nlayers) ylimit  = c(ylimit,emean[[vnames[l]]])
            ylimit = pretty.xylim(u=ylimit,fracexp=0.0,is.log=plog)
         }else{
            ylimit = ylimit.fix
         }#end if
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
         for (o in sequence(nout)){
            #------ Open file. ------------------------------------------------------------#
            fichier = file.path(outdir,paste0(prefix,"-",suffix,".",outform[o]))
            dummy   = open.plot( fichier = fichier
                               , outform = outform[o]
                               , size    = eysize
                               , ptsz    = ptsz
                               , depth   = depth
                               )#end open.plot
            #------------------------------------------------------------------------------#



            #----- Plot settings. ---------------------------------------------------------#
            letitre = paste0(" Time series: ",group,"\n",lieu)
            ley     = desc.unit(desc=group,unit=unit)
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #     Split the plot into two windows.                                         #
            #------------------------------------------------------------------------------#
            par(par.user)
            layout(mat=rbind(2,1),heights=c(1.-f.leg,f.leg))
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #      First plot: legend.                                                     #
            #------------------------------------------------------------------------------#
            par(mar=c(0.1,4.6,0.1,2.1))
            plot.new()
            plot.window(xlim=c(0,1),ylim=c(0,1))
            legend( x      = "bottom"
                  , inset  = 0.0
                  , legend = description
                  , col    = lcolours
                  , lwd    = llwd
                  , ncol   = min(3,pretty.box(nlayers)$ncol)
                  , xpd    = TRUE
                  , bty    = "n"
                  )#end legend
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #      Main plot.                                                              #
            #------------------------------------------------------------------------------#
            par(mar=c(4.1,4.6,4.1,2.1))
            plot.new()
            plot.window(xlim=xlimit,ylim=ylimit,log=xylog)
            axis(side=1,at=whenplot8$levels,labels=whenplot8$labels,padj=whenplot8$padj)
            axis(side=2,las=1)
            box()
            title(main=letitre,xlab="Year",ylab=ley,cex.main=0.7)
            if (drought.mark){
               for (n in sequence(ndrought)){
                  rect(xleft  = drought[[n]][1],ybottom = ydrought[1]
                      ,xright = drought[[n]][2],ytop    = ydrought[2]
                      ,col    = grid.colour,border=NA)
               }#end for
            }#end if
            #----- Plot grid. -------------------------------------------------------------#
            if (plotgrid){ 
               abline(v=whenplot8$levels,h=axTicks(side=2),col=grid.colour,lty="solid")
            }#end if
            #----- Plot lines. ------------------------------------------------------------#
            for (l in sequence(nlayers)){
               thisvar = emean[[vnames[l]]]
               points(x=datum$tomonth,y=thisvar,col=lcolours[l],lwd=llwd[l],type=ltype
                     ,pch=16,cex=0.8)
            }#end for
            #------------------------------------------------------------------------------#


            #----- Close the device. ------------------------------------------------------#
            dummy = close.plot(outform=outform[o])
            #------------------------------------------------------------------------------#
         } #end for outform
         #---------------------------------------------------------------------------------#
      }#end if plotit
      #------------------------------------------------------------------------------------#
   }#end for ntser
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #   Plot the time series diagrams showing months and years.                             #
   #---------------------------------------------------------------------------------------#
   cat0("      * Plot some monthly means of groups of variables.")
   for (hh in sequence(ntheme)){

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
      ylimit.fix   = themenow$mmean.lim
   
      if (plotit){

         #---------------------------------------------------------------------------------#
         #    Check whether the time series directory exists.  If not, create it.          #
         #---------------------------------------------------------------------------------#
         outdir = file.path(outpref,"theme_mmean")
         if (! file.exists(outdir)) dir.create(outdir)
         cat0("      + ",group," time series for several variables.")


         #----- Define the number of layers. ----------------------------------------------#
         nlayers   = length(vnames)
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Find the limit, make some room for the legend, and in case the field is a   #
         # constant, nudge the limits so the plot command will not complain.               #
         #---------------------------------------------------------------------------------#
         xlimit    = pretty.xylim(u=montmont,fracexp=0.0,is.log=plog)
         if (any (! is.finite(ylimit.fix))){
            ylimit    = NULL
            for (l in 1:nlayers) ylimit  = c(ylimit,mmean[[vnames[l]]])
            ylimit = pretty.xylim(u=ylimit,fracexp=0.0,is.log=plog)
         }else{
            ylimit = ylimit.fix
         }#end if
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
         for (o in sequence(nout)){
            #------ Open file. ------------------------------------------------------------#
            fichier = file.path(outdir,paste0(prefix,"-",suffix,".",outform[o]))
            dummy   = open.plot( fichier = fichier
                               , outform = outform[o]
                               , size    = eysize
                               , ptsz    = ptsz
                               , depth   = depth
                               )#end open.plot
            #------------------------------------------------------------------------------#



            #----- Plot settings. ---------------------------------------------------------#
            letitre = paste0(" Time series: ",group,"\n",lieu)
            ley     = desc.unit(desc=group,unit=unit)
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #     Split the plot into two windows.                                         #
            #------------------------------------------------------------------------------#
            par(par.user)
            layout(mat=rbind(2,1),heights=c(1.-f.leg,f.leg))
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #      First plot: legend.                                                     #
            #------------------------------------------------------------------------------#
            par(mar=c(0.1,4.6,0.1,2.1))
            plot.new()
            plot.window(xlim=c(0,1),ylim=c(0,1))
            legend( x      = "bottom"
                  , inset  = 0.0
                  , legend = description
                  , col    = lcolours
                  , lwd    = llwd
                  , pch    = 16
                  , ncol   = min(3,pretty.box(nlayers)$ncol)
                  , cex    = 0.9*cex.ptsz
                  , xpd    = TRUE
                  , bty    = "n"
                  )#end legend
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #      Main plot.                                                              #
            #------------------------------------------------------------------------------#
            par(mar=c(4.1,4.6,4.1,2.1))
            plot.new()
            plot.window(xlim=xlimit,ylim=ylimit,log=xylog)
            axis(side=1,at=mplot$levels,labels=substring(mplot$labels,1,1),padj=mplot$padj)
            axis(side=2,las=1)
            box()
            title(main=letitre,xlab="Year",ylab=ley,cex.main=0.7)
            if (drought.mark){
               for (n in sequence(ndrought)){
                  rect(xleft  = drought[[n]][1],ybottom = ydrought[1]
                      ,xright = drought[[n]][2],ytop    = ydrought[2]
                      ,col    = grid.colour,border=NA)
               }#end for
            }#end if
            #----- Plot grid. -------------------------------------------------------------#
            if (plotgrid){ 
               abline(v=mplot$levels,h=axTicks(side=2),col=grid.colour,lty="solid")
            }#end if
            #----- Plot lines. ------------------------------------------------------------#
            for (l in sequence(nlayers)){
               thisvar = mmean[[vnames[l]]]
               points(x=montmont,y=thisvar,col=lcolours[l],lwd=llwd[l],type=ltype
                     ,pch=16,cex=0.8)
            }#end for
            #------------------------------------------------------------------------------#


            #----- Close the device. ------------------------------------------------------#
            dummy = close.plot(outform=outform[o])
            #------------------------------------------------------------------------------#
         } #end for outform
         #---------------------------------------------------------------------------------#
      }#end if plotit
      #------------------------------------------------------------------------------------#
   }#end for ntser
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #   Plot the climatology of the mean diurnal cycle.                                     #
   #---------------------------------------------------------------------------------------#
   cat0("      * Plot the mean diel of groups of variables.")
   for (hh in sequence(ntheme)){

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
         outdir   = file.path(outpref,"theme_qmean")
         if (! file.exists(outdir)) dir.create(outdir)
         outtheme = file.path(outdir,prefix)
         if (! file.exists(outtheme)) dir.create(outtheme)
         cat0("      + ",group," diurnal cycle for several variables.")


         #----- Define the number of layers. ----------------------------------------------#
         nlayers   = length(vnames)
         xlimit    = range(thisday)
         ylimit    = NULL
         for (l in sequence(nlayers)) ylimit = c(ylimit,umean[[vnames[l]]])
         ylimit = pretty.xylim(u=ylimit,fracexp=0.0,is.log=FALSE)
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #      Loop over all months.                                                      #
         #---------------------------------------------------------------------------------#
         for (pmon in sequence(12)){
            cmon    = sprintf("%2.2i",pmon)
            namemon = mlist[pmon]

            #------------------------------------------------------------------------------#
            #     Check if the directory exists.  If not, create it.                       #
            #------------------------------------------------------------------------------#

            #----- Loop over formats. -----------------------------------------------------#
            for (o in sequence(nout)){
               #------ Open file. ---------------------------------------------------------#
               fichier = file.path( outtheme
                                  , paste0(prefix,"-",cmon,"-",suffix,".",outform[o])
                                  )#end file.path
               dummy   = open.plot( fichier = fichier
                                  , outform = outform[o]
                                  , size    = eysize
                                  , ptsz    = ptsz
                                  , depth   = depth
                                  )#end open.plot
               #---------------------------------------------------------------------------#



               #----- Plot settings. ------------------------------------------------------#
               letitre = paste0(group," - ",lieu,"\n","Mean diurnal cycle - ",namemon)
               ley     = desc.unit(desc=group,unit=unit)
               #---------------------------------------------------------------------------#


               #---------------------------------------------------------------------------#
               #     Split the plot into two windows.                                      #
               #---------------------------------------------------------------------------#
               par(par.user)
               layout(mat=rbind(2,1),heights=c(1.-f.leg,f.leg))
               #------------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #      First plot: legend.                                                  #
               #---------------------------------------------------------------------------#
               par(mar=c(0.1,4.6,0.1,2.1))
               plot.new()
               plot.window(xlim=c(0,1),ylim=c(0,1))
               legend( x      = "bottom"
                     , inset  = 0.0
                     , legend = description
                     , col    = lcolours
                     , lwd    = llwd
                     , ncol   = min(3,pretty.box(nlayers)$ncol)
                     , xpd    = TRUE
                     , bty    = "n"
                     )#end legend
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #      Main plot.                                                           #
               #---------------------------------------------------------------------------#
               par(mar=c(4.1,4.6,4.1,2.1))
               plot.new()
               plot.window(xlim=xlimit,ylim=ylimit,log=xylog)
               axis(side=1,at=uplot$levels,labels=uplot$labels,padj=uplot$padj)
               axis(side=2,las=1)
               box()
               title(main=letitre,xlab="Year",ylab=ley,cex.main=0.7)
               if (drought.mark){
                  for (n in sequence(ndrought)){
                     rect(xleft  = drought[[n]][1],ybottom = ydrought[1]
                         ,xright = drought[[n]][2],ytop    = ydrought[2]
                         ,col    = grid.colour,border=NA)
                  }#end for
               }#end if
               #----- Plot grid. ----------------------------------------------------------#
               if (plotgrid){ 
                  abline(v=uplot$levels,h=axTicks(side=2),col=grid.colour,lty="solid")
               }#end if
               #----- Plot lines. ---------------------------------------------------------#
               for (l in sequence(nlayers)){
                  thisvar = umean[[vnames[l]]]
                  thisvar = cbind(thisvar[,ndcycle],thisvar)
                  points(x=thisday,y=thisvar[pmon,],col=lcolours[l]
                        ,lwd=llwd[l],type=ltype,pch=16)
               }#end for
               #---------------------------------------------------------------------------#


               #----- Close the device. ---------------------------------------------------#
               dummy = close.plot(outform=outform[o])
               #---------------------------------------------------------------------------#
            } #end for outform
            #------------------------------------------------------------------------------#
         }#end for pmon
         #---------------------------------------------------------------------------------#
      }#end if plotit
      #------------------------------------------------------------------------------------#
   }#end for ntser
   #---------------------------------------------------------------------------------------#





   #---------------------------------------------------------------------------------------#
   #   Plot the climatology of the soil properties.                                        #
   #---------------------------------------------------------------------------------------#
   for (v in sequence(nsoilplot)){

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
         outdir  =  file.path(outpref,"soil_mmean")
         if (! file.exists(outdir)) dir.create(outdir)
         cat0("      + Climatology profile of ",description,".")

         #----- Find the number of rows and columns, and the axes. ------------------------#
         monaxis  = sort(unique(datum$month))
         soilaxis = slz
         nmon     = length(monaxis)
         nsoil    = nzg

         #----- Save the meaningful months and years. -------------------------------------#
         monat   = sequence(12)
         monlab  = substring(month.abb,1,1)

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
            vlevels = pretty.log(x=vrange,n=ncolsfc,forcelog=TRUE)
            vnlev   = length(vlevels)
         }else{
            vrange  = range(varbuff,na.rm=TRUE)
            vlevels = pretty(x=vrange,n=ncolsfc)
            vnlev   = length(vlevels)
         }#end if

         #----- Loop over formats. --------------------------------------------------------#
         for (o in sequence(nout)){
            fichier = file.path(outdir,paste0(vnam,"-",suffix,".",outform[o]))
            dummy   = open.plot( fichier = fichier
                               , outform = outform[o]
                               , size    = exsize
                               , ptsz    = ptsz
                               , depth   = depth
                               )#end open.plot

            letitre = paste0(description," - ",lieu)
            ley     = desc.unit(desc="Soil depth",unit=untab$m)
            lacle   = desc.unit(desc=NULL,unit=unit)
            par(par.user)
            sombreado(x=monaxis,y=soilaxis,z=varbuff,levels=vlevels,nlevels=vnlev
                     ,colour.palette=get(vcscheme)
                     ,plot.title=title(main=letitre,xlab="Month",ylab=ley,cex.main=0.7)
                     ,key.title=title(main=lacle,cex.main=0.8)
                     ,key.log=pnlog
                     ,useRaster=TRUE
                     ,plot.axes={axis(side=1,at=monat,labels=monlab)
                                 axis(side=2,at=zat,labels=znice)
                                 if (fcgrid){
                                    abline(h=zat,v=monat,col=grid.colour,lty="dotted")
                                 }#end if fcgrid
                                }#end plot.axes
                     )

            dummy = close.plot(outform=outform[o])
         } #end for outform
      }#end if plotit
   }#end for (v in 1:nsoilplot)
   #---------------------------------------------------------------------------------------#





   #---------------------------------------------------------------------------------------#
   #   Plot the climatology of the soil properties.                                        #
   #---------------------------------------------------------------------------------------#
   for (sts in sequence(nsoilplot)){

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
         outdir  =  file.path(outpref,"soil_emean")
         if (! file.exists(outdir)) dir.create(outdir)
         cat0("      + Time series profile of ",description,".")

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
            vlevels = pretty.log(x=vrange,n=ncolsfc,forcelog=TRUE)
            vnlev   = length(vlevels)
         }else{
            vrange  = range(varbuff,na.rm=TRUE)
            vlevels = pretty(x=vrange,n=ncolsfc)
            vnlev   = length(vlevels)
         }#end if

         #----- Loop over formats. --------------------------------------------------------#
         for (o in sequence(nout)){
            fichier = file.path(outdir,paste0(vnam,"-",suffix,".",outform[o]))
            dummy   = open.plot( fichier = fichier
                               , outform = outform[o]
                               , size    = exsize
                               , ptsz    = ptsz
                               , depth   = depth
                               )#end open.plot

            letitre = paste0(description," - ",lieu)
            ley     = desc.unit(desc="Soil depth",unit=untab$m)
            lacle   = desc.unit(desc=NULL,unit=unit)
            par(par.user)
            sombreado(x=whenaxis,y=soilaxis,z=varbuff,levels=vlevels,nlevels=vnlev
                     ,colour.palette=get(vcscheme)
                     ,plot.title=title(main=letitre,xlab="Month",ylab=ley,cex.main=0.7)
                     ,key.title=title(main=lacle,cex.main=0.8)
                     ,key.log=pnlog
                     ,useRaster=TRUE
                     ,plot.axes={axis(side=1,at=whenplot6$levels
                                     ,labels=whenplot6$labels,padj=whenplot6$padj)
                                 axis(side=2,at=zat,labels=znice)
                                 if (fcgrid){
                                    abline(h=zat,v=whenplot6$levels,col=grid.colour
                                          ,lty="dotted")
                                 }#end if fcgrid
                                }#end plot.axes
                     )#end sombreado

            dummy = close.plot(outform=outform[o])
         } #end for outform
      }#end if plotit
   }#end for nhov
   #---------------------------------------------------------------------------------------#





   #---------------------------------------------------------------------------------------#
   #   Plot a filled contour plot showing months and years.                                #
   #---------------------------------------------------------------------------------------#
   for (v in sequence(nsqueeze)){

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
         outdir  =  file.path(outpref,"fillc_mmean")
         if (! file.exists(outdir)) dir.create(outdir)
         cat0("      + ",description," time series in filled contour.")

         #----- Load this variable into "thisvar". ----------------------------------------#
         thisvar = emean[[vnam]]

         #----- Find the number of rows and columns, and the axes. ------------------------#
         monaxis = sort(unique(datum$month[sel]))
         yraxis  = sort(unique(datum$year[sel]))
         nmon    = length(monaxis)
         nyear   = length(yraxis)

         #----- Save the meaningful months and years. -------------------------------------#
         monat   = sequence(12)
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
         for (o in sequence(nout)){
            fichier = file.path(outdir,paste0(vnam,"-",suffix,".",outform[o]))
            dummy   = open.plot( fichier = fichier
                               , outform = outform[o]
                               , size    = exsize
                               , ptsz    = ptsz
                               , depth   = depth
                               )#end open.plot

            letitre = paste0(description," - ",lieu)
            lacle   = desc.unit(desc=NULL,unit=unit)
            par(par.user)
            sombreado(x=monaxis,y=yraxis,z=varbuff,levels=vlevels,nlevels=vnlev
                     ,colour.palette=get(vcscheme)
                     ,plot.title=title(main=letitre,xlab="Month",ylab="Year",cex.main=0.7)
                     ,key.title=title(main=lacle,cex.main=0.8)
                     ,useRaster=TRUE
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

            dummy = close.plot(outform=outform[o])
         } #end for outform
      }#end if plotit
   }#end for nhov
   #---------------------------------------------------------------------------------------#






   #---------------------------------------------------------------------------------------#
   #   Plot the filled contour diagrams showing time of day and time.                      #
   #---------------------------------------------------------------------------------------#
   for (v in sequence(nsqueeze)){

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
         outdir  =  file.path(outpref,"fillc_qmean")
         if (! file.exists(outdir)) dir.create(outdir)
         cat0("      + ",description," time series of diurnal cycle.")

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
         for (o in sequence(nout)){
            fichier = file.path(outdir,paste0(vnam,"-",suffix,".",outform[o]))
            dummy   = open.plot( fichier = fichier
                               , outform = outform[o]
                               , size    = exsize
                               , ptsz    = ptsz
                               , depth   = depth
                               )#end open.plot

            letitre = paste0("Mean diurnal cycle \n ",description," - ",lieu)
            ley     = desc.unit(desc="Time of day",unit=untab$gmt)
            lacle   = desc.unit(desc=NULL         ,unit=unit)
            par(par.user)
            sombreado(x=whenaxis,y=hraxis,z=varbuff,levels=vlevels,nlevels=vnlev
                     ,colour.palette=get(vcscheme)
                     ,plot.title=title(main=letitre,ylab=ley,xlab="Time",cex.main=0.7)
                     ,key.title=title(main=lacle,cex.main=0.8)
                     ,useRaster=TRUE
                     ,plot.axes={axis(side=1,at=huplot$level,labels=huplot$labels)
                                 axis(side=2,at=uplot$levels,labels=uplot$labels)
                                 if (fcgrid){
                                    abline(v=huplot$levels,h=uplot$levels
                                          ,col=grid.colour,lty="dotted")
                                 }#end if fcgrid
                                }#end plot.axes
                     )

            dummy = close.plot(outform=outform[o])
         } #end for outform
      }#end if plotit
   }#end for nhov
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #   Plot the monthly boxplots.                                                          #
   #---------------------------------------------------------------------------------------#
   for (v in sequence(nsqueeze)){

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
         outdir  =  file.path(outpref,"boxplot")
         if (! file.exists(outdir)) dir.create(outdir)
         cat0("      + ",description," box plot.")

         #----- Load this variable into "thisvar". ----------------------------------------#
         thisvar = emean[[vnam]]

         for (o in sequence(nout)){
            fichier = file.path(outdir,paste0(vnam,"-",suffix,".",outform[o]))
            dummy   = open.plot( fichier = fichier
                               , outform = outform[o]
                               , size    = sqsize
                               , ptsz    = ptsz
                               , depth   = depth
                               )#end open.plot

            ylimit  = pretty.xylim(u=thisvar,fracexp=0.0,is.log=FALSE)
            letitre = paste(description,lieu,sep=" - ")
            ley     = desc.unit(desc=description,unit=unit)
            par(par.user)
            plot(mmonth,thisvar,main=letitre,ylim=ylimit,cex.main=0.7
                ,xlab="Time",ylab=ley,las=1)

            dummy = close.plot(outform=outform[o])
         } #end for outform
      }#end if
   }#end for nbox
   #---------------------------------------------------------------------------------------#





   #---------------------------------------------------------------------------------------#
   #   Plot the PDF of patch-level properties as a function of time.                       #
   #---------------------------------------------------------------------------------------#
   cat0("      * Time series of PDF of properties by patch.")
   for (v in sequence(nplotpatch)){

      #----- Retrieve variable information from the list. ---------------------------------#
      thispatch   = plotpatch[[v]]
      vnam        = thispatch$vnam
      description = thispatch$desc
      unit        = thispatch$unit
      vcscheme    = thispatch$col.scheme
      plog        = thispatch$plog
      plotit      = thispatch$emean

      this        = patchpdf[[vnam]]$edensity
      plotit      = ( plotit && any(is.finite(this$x),na.rm=TRUE)
                             && any(is.finite(this$y),na.rm=TRUE) 
                             && any(is.finite(this$z),na.rm=TRUE) )


      #----- Select polygon average. ------------------------------------------------------#
      if (vnam %in% names(emean)){
         thisemean  = emean[[vnam]]
      }else{
         thisemean  = rep(NA_real_,times=datum$when)
      }#end if (vnam %in% names(emean))
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Find levels, and expand PDF scale in case it is a constant.                    #
      #------------------------------------------------------------------------------------#
      if (plog){
         vlevs = sort(unique(pretty.log(this$z,n=ncolsfc,forcelog=TRUE)))
      }else{
         vlevs = sort(unique(pretty(this$z,n=ncolsfc)))
      }#end if
      if (length(vlevs) == 1) vlevs = pretty.xylim(u=vlevs,fracexp=0.0,is.log=plog)
      #------------------------------------------------------------------------------------#

      if (plotit){

         #---------------------------------------------------------------------------------#
         #     Check if the directory exists.  If not, create it.                          #
         #---------------------------------------------------------------------------------#
         outdir  =  file.path(outpref,"patch_emean")
         if (! file.exists(outdir)) dir.create(outdir)
         cat0("      + PDF plot of ",description,".")


         #----- Limits for y axis. --------------------------------------------------------#
         ylimit  = pretty.xylim(u=this$y,fracexp=c(-0.04,0.04))
         yat     = pretty(ylimit)
         ylabels = sprintf("%g",yat)
         #---------------------------------------------------------------------------------#


         #----- Set annotation. -----------------------------------------------------------#
         zlimit  = pretty.xylim(u=vlevs,is.log=plog)
         if (plog){
            zat  = pretty.log(zlimit)
            zlwr = zlimit[1] * (zlimit[2]/zlimit[1]) ^ sqrt(.Machine$double.eps)
            zupr = zlimit[2] * (zlimit[1]/zlimit[2]) ^ sqrt(.Machine$double.eps)
         }else{
            zlwr = zlimit[1] + sqrt(.Machine$double.eps) * diff(zlimit)
            zupr = zlimit[2] - sqrt(.Machine$double.eps) * diff(zlimit)
            zat  = pretty(zlimit)
         }#end if (plog)
         zlabels = sprintf("%g",zat)
         this$z  = pmax(zlwr,pmin(zupr,this$z)) + 0. * this$z
         #---------------------------------------------------------------------------------#


         #----- Set annotation. -----------------------------------------------------------#
         letitre = paste0("Density function of ",description," \ ",lieu)
         lex     = "Time"
         ley     = desc.unit(desc=description,unit=unit)
         #---------------------------------------------------------------------------------#



         #----- Commands to plot after. ---------------------------------------------------#
         plot.after = list( abline = list( v   = whenplot8$levels
                                         , h   = yat
                                         , col = if(fcgrid){grid.colour}else{"transparent"}
                                         , lty = if(fcgrid){"dotted"   }else{"blank"      }
                                         )#end abline
                          , lines  = list( x   = datum$when
                                         , y   = thisemean
                                         , col = "grey16"
                                         , lwd = 2.0
                                         )#end lines
                          )#end list 
         #---------------------------------------------------------------------------------#



         #----- Loop over formats. --------------------------------------------------------#
         for (o in sequence(nout)){
            fichier = file.path(outdir,paste0(vnam,"-",suffix,".",outform[o]))
            dummy   = open.plot( fichier = fichier
                               , outform = outform[o]
                               , size    = exsize
                               , ptsz    = ptsz
                               , depth   = depth
                               )#end open.plot


            #------------------------------------------------------------------------------#
            #     Plot the PDF distribution.                                               #
            #------------------------------------------------------------------------------#
            par(par.user)
            gridded.plot( x                = as.numeric(this$x)
                        , y                = this$y
                        , z                = this$z
                        , ylim             = ylimit
                        , levels           = vlevs
                        , colour.palette   = get(vcscheme)
                        , main.title       = list( main     = letitre
                                                 , xlab     = lex
                                                 , ylab     = ley
                                                 , cex.main = 0.7
                                                 )#end list
                        , key.title        = list( main     = "Density"
                                                 , cex.main = 0.8
                                                 )#end list
                        , key.log          = plog 
                        , useRaster        = FALSE
                        , na.col           = "black"
                        , x.axis.options   = list( side   = 1
                                                 , las    = 1
                                                 , at     = whenplot8$levels
                                                 , labels = whenplot8$labels
                                                 , padj   = whenplot8$padj
                                                 )#end list
                        , y.axis.options   = list( side   = 2
                                                 , las    = 1
                                                 , at     = yat
                                                 , labels = ylabels
                                                 )#end list
                        , key.axis.options = list( side   = 4
                                                 , las    = 1
                                                 , at     = zat
                                                 , labels = zlabels
                                                 )#end list
                        , plot.after       = plot.after
                        , mar              = c(4.1,4.6,2.1,0.6)
                        , mar.key          = c(4.1,0.6,2.1,4.1)
                        )#end gridded.plot
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #     Close the device.                                                        #
            #------------------------------------------------------------------------------#
            dummy = close.plot(outform=outform[o])
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
   cat0("      * Monthly PDF of properties by patch.")
   for (v in sequence(nplotpatch)){

      #----- Retrieve variable information from the list. ---------------------------------#
      thispatch   = plotpatch[[v]]
      vnam        = thispatch$vnam
      description = thispatch$desc
      unit        = thispatch$unit
      vcscheme    = thispatch$col.scheme
      plog        = thispatch$plog
      plotit      = thispatch$mmean

      this        = patchpdf[[vnam]]$mdensity
      plotit      = ( plotit && any(is.finite(this$x),na.rm=TRUE)
                             && any(is.finite(this$y),na.rm=TRUE) 
                             && any(is.finite(this$z),na.rm=TRUE) )


      #----- Select polygon average. ------------------------------------------------------#
      if (vnam %in% names(mmean)){
         thismmean  = mmean[[vnam]]
      }else{
         thismmean  = rep(NA_real_,times=this$x)
      }#end if (vnam %in% names(emean))
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Find levels, and expand PDF scale in case it is a constant.                    #
      #------------------------------------------------------------------------------------#
      if (plog){
         vlevs = sort(unique(pretty.log(this$z,n=ncolsfc,forcelog=TRUE)))
      }else{
         vlevs = sort(unique(pretty(this$z,n=ncolsfc)))
      }#end if
      if (length(vlevs) == 1) vlevs = pretty.xylim(u=vlevs,fracexp=0.0,is.log=plog)
      #------------------------------------------------------------------------------------#

      if (plotit){

         #---------------------------------------------------------------------------------#
         #     Check if the directory exists.  If not, create it.                          #
         #---------------------------------------------------------------------------------#
         outdir  =  file.path(outpref,"patch_mmean")
         if (! file.exists(outdir)) dir.create(outdir)
         cat0("      + PDF plot of ",description,".")


         #----- Find the month tick marks. ------------------------------------------------#
         monlim = c(0.5,12.5)
         monat  = sequence(12)
         monlab = c("J","F","M","A","M","J","J","A","S","O","N","D")
         #---------------------------------------------------------------------------------#


         #----- Limits for y axis. --------------------------------------------------------#
         ylimit  = pretty.xylim(u=this$y,fracexp=c(-0.04,0.04))
         yat     = pretty(ylimit)
         ylabels = sprintf("%g",yat)
         #---------------------------------------------------------------------------------#


         #----- Set annotation. -----------------------------------------------------------#
         zlimit  = pretty.xylim(u=vlevs,is.log=plog)
         if (plog){
            zat  = pretty.log(zlimit)
            zlwr = zlimit[1] * (zlimit[2]/zlimit[1]) ^ sqrt(.Machine$double.eps)
            zupr = zlimit[2] * (zlimit[1]/zlimit[2]) ^ sqrt(.Machine$double.eps)
         }else{
            zlwr = zlimit[1] + sqrt(.Machine$double.eps) * diff(zlimit)
            zupr = zlimit[2] - sqrt(.Machine$double.eps) * diff(zlimit)
            zat  = pretty(zlimit)
         }#end if (plog)
         zlabels = sprintf("%g",zat)
         this$z  = pmax(zlwr,pmin(zupr,this$z)) + 0. * this$z
         #---------------------------------------------------------------------------------#


         #----- Set annotation. -----------------------------------------------------------#
         letitre = paste0("Density function of ",description," \ ",lieu)
         lex     = "Month"
         ley     = desc.unit(desc=description,unit=unit)
         #---------------------------------------------------------------------------------#



         #----- Commands to plot after. ---------------------------------------------------#
         plot.after = list( abline = list( v    = monat
                                         , h    = yat
                                         , col  = if(fcgrid){grid.colour}else{"transparent"}
                                         , lty  = if(fcgrid){"dotted"   }else{"blank"      }
                                         )#end abline
                          , lines  = list( x    = this$x
                                         , y    = thismmean
                                         , col  = "grey16"
                                         , lwd  = 2.0
                                         , type = "o"
                                         , pch  = 16
                                         )#end lines
                          )#end list 
         #---------------------------------------------------------------------------------#


         #----- Loop over formats. --------------------------------------------------------#
         for (o in sequence(nout)){
            fichier = file.path(outdir,paste0(vnam,"-",suffix,".",outform[o]))
            dummy   = open.plot( fichier = fichier
                               , outform = outform[o]
                               , size    = exsize
                               , ptsz    = ptsz
                               , depth   = depth
                               )#end open.plot


            #------------------------------------------------------------------------------#
            #     Plot the PDF distribution.                                               #
            #------------------------------------------------------------------------------#
            par(par.user)
            gridded.plot( x                = this$x
                        , y                = this$y
                        , z                = this$z
                        , xlim             = monlim
                        , ylim             = ylimit
                        , levels           = vlevs
                        , colour.palette   = get(vcscheme)
                        , main.title       = list( main     = letitre
                                                 , xlab     = lex
                                                 , ylab     = ley
                                                 , cex.main = 0.7
                                                 )#end list
                        , key.title        = list(main="Density",cex.main=0.8)
                        , key.log          = plog 
                        , useRaster        = TRUE
                        , na.col           = "black"
                        , x.axis.options   = list( side   = 1
                                                 , las    = 1
                                                 , at     = monat
                                                 , labels = monlab
                                                 )#end list
                        , y.axis.options   = list( side   = 2
                                                 , las    = 1
                                                 , at     = yat
                                                 , labels = ylabels
                                                 )#end list
                        , key.axis.options = list( side   = 4
                                                 , las    = 1
                                                 , at     = zat
                                                 , labels = zlabels
                                                 )#end list
                        , plot.after       = plot.after
                        , mar              = c(4.1,4.6,2.1,0.6)
                        , mar.key          = c(4.1,0.6,2.1,4.1)
                        )#end gridded.plot
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #     Close the device.                                                        #
            #------------------------------------------------------------------------------#
            dummy = close.plot(outform=outform[o])
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
   cat0("    + Bar plot by DBH classes.")
   monbplot    = which(nummonths(datum$tomonth) %in% sasmonth)
   nmonbplot   = length(monbplot)
   pftuse      = which(apply(X=szpft$nplant,MARGIN=3,FUN=sum,na.rm=TRUE) > 0)
   pftuse      = pftuse[pftuse != (npft+1)]
   npftuse     = length(pftuse)
   pftname.use = pft$name  [pftuse]
   pftcol.use  = pft$colour[pftuse]
   for (v in sequence(ntspftdbh)){
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
         cat0("      - ",description,".")


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
         ylimit = pretty.xylim(u=ylimit,fracexp=0.0,is.log=plog)
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Check if the directory exists.  If not, create it.                          #
         #---------------------------------------------------------------------------------#
         barplotdir = file.path(outpref,"barplot_dbh")
         if (! file.exists(barplotdir)) dir.create(barplotdir)
         outdir = file.path(barplotdir,vnam)
         if (! file.exists(outdir)) dir.create(outdir)
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Loop over all possible months.                                             #
         #---------------------------------------------------------------------------------#
         for (m in sequence(nmonbplot)){

            #----- Find which year we are plotting. ---------------------------------------#
            cmonth    = sprintf("%2.2i",(nummonths(thiswhen[m])))
            cyear     = sprintf("%4.4i",(numyears(thiswhen[m])))
            mm        = as.numeric(cmonth)
            yy        = as.numeric(cyear)
            whentitle = paste(mon2mmm(mm,cap1=TRUE),cyear,sep="-")
            #------------------------------------------------------------------------------#


            #----- Loop over output formats. ----------------------------------------------#
            for (o in sequence(nout)){
               #------ Open the plot. -----------------------------------------------------#
               fichier = file.path(outdir
                                  ,paste0(vnam,"-",cyear,"-",cmonth,"-",suffix
                                         ,".",outform[o])
                                  )#end file.path
               dummy   = open.plot( fichier = fichier
                                  , outform = outform[o]
                                  , size    = eysize
                                  , ptsz    = ptsz
                                  , depth   = depth
                                  )#end open.plot
               #---------------------------------------------------------------------------#


               #------ Set up the title and axis labels. ----------------------------------#
               letitre = paste0(lieu,"\n",description," - Time : ",whentitle)
               lexlab  = "DBH Classes"
               leylab  = desc.unit(desc=description,unit=unit)
               #---------------------------------------------------------------------------#


               #------ Split window. ------------------------------------------------------#
               par(par.user)
               layout(mat=rbind(2,1),heights=c(1.-f.leg,f.leg))
               #---------------------------------------------------------------------------#



               #------ Legend. ------------------------------------------------------------#
               par(mar=c(0.1,4.6,0.1,2.1))
               plot.new()
               plot.window(xlim=c(0,1),ylim=c(0,1))
               legend( x      = "bottom"
                     , inset  = 0.0
                     , legend = pftname.use
                     , fill   = pftcol.use
                     , ncol   = min(3,pretty.box(n.selpft)$ncol)
                     , title  = expression(bold("Plant functional type"))
                     , cex    = cex.ptsz
                     , bg     = background
                     , xpd    = TRUE
                     , bty    = "n"
                     )#end legend
               #---------------------------------------------------------------------------#


               #----- Plot all monthly means together. ------------------------------------#
               par(mar=c(4.1,4.6,4.1,2.1))
               barplot(height=t(thisvnam[m,,]),names.arg=dbhnames[1:ndbh],width=1.0
                      ,main=letitre,xlab=lexlab,ylab=leylab,ylim=ylimit,legend.text=FALSE
                      ,beside=(! stacked),col=pftcol.use,log=xylog
                      ,border=grey.fg,xpd=FALSE,cex.main=0.7,las=1)
               if (plotgrid & (! stacked)){
                  xgrid=0.5+(1:ndbh)*(1+npftuse)
                  abline(v=xgrid,col=grid.colour,lty="solid")
               }#end if
               box()
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     Close the device.                                                     #
               #---------------------------------------------------------------------------#
               dummy = close.plot(outform=outform[o])
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
   #    Plot the 3-D size and age structure of various variables.                          #
   #---------------------------------------------------------------------------------------#
   for (v in sequence(ntspftdbh)){
      #----- Retrieve variable information from the list. ---------------------------------#
      thissas     = tspftdbh[[v]]
      vnam        = thissas$vnam
      description = thissas$desc
      unit        = thissas$i.unit
      plotit      = thissas$sas
      plog        = thissas$plog

      #----- If this variable is to be plotted, then go through this if block. ------------#
      if (plotit){

         cat0("      + Size and age structure plot: ",description,".")

         #---------------------------------------------------------------------------------#
         #     Check if the directory exists.  If not, create it.                          #
         #---------------------------------------------------------------------------------#
         sasdir = file.path(outpref,"sas")
         if (! file.exists(sasdir)) dir.create(sasdir)
         outdir = file.path(sasdir,vnam)
         if (! file.exists(outdir)) dir.create(outdir)
         #---------------------------------------------------------------------------------#


         #----- Load this list into "thislist". -------------------------------------------#
         varco =  cohort[[vnam]]
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #      Loop over all times.                                                       #
         #---------------------------------------------------------------------------------#
         for (ww in names(cohort$age)){

            #----- Find which year we are plotting. ---------------------------------------#
            cmonth   = substring(ww,7,8)
            thisyear = substring(ww,2,5)
            mm       = as.numeric(cmonth)
            yy       = as.numeric(thisyear)

            #----- Retrieve variable list, age, DBH, and PFT for this year. ---------------#
            ageww   = cohort$age   [[ww]]
            if (any(ageww <= 0,na.rm=TRUE)){
               minww = min(ageww,na.rm=TRUE)
               ageww = ageww - minww + 0.01
            }#end if
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
                  xlimit  = pretty.xylim(u=unlist(cohort$age),fracexp=0.0,is.log=TRUE )
                  ylimit  = pretty.xylim(u=unlist(cohort$dbh),fracexp=0.0,is.log=FALSE)
                  zlimit  = pretty.xylim(u=unlist(varco)     ,fracexp=0.0,is.log=plog )
                  popmin  = min  (unlist(cohort$nplant * cohort$area), na.rm=TRUE)
                  popmax  = max  (unlist(cohort$nplant * cohort$area), na.rm=TRUE)
               }else{
                  xlimit  = pretty.xylim(u=ageww             ,fracexp=0.0,is.log=TRUE )
                  ylimit  = pretty.xylim(u=dbhww             ,fracexp=0.0,is.log=FALSE)
                  zlimit  = pretty.xylim(u=varww             ,fracexp=0.0,is.log=plog )
                  popmin  = min  (popww  ,na.rm=TRUE)
                  popmax  = max  (popww  ,na.rm=TRUE)
               }#end if
               #---------------------------------------------------------------------------#


               #----- Define the scale-dependent population size. -------------------------#
               cexww = cexmin + (cexmax - cexmin) * log(popww/popmin) / log(popmax/popmin)
               #---------------------------------------------------------------------------#



               #----- Define the floor location. ------------------------------------------#
               if ((zlimit[1] > 0) != (zlimit[2] > 0)){
                  floor3d = 0.
               }else if (zlimit[1] > 0){
                  floor3d = zlimit[1]
               }else{
                  floor3d = zlimit[2]
               }#end if
               #---------------------------------------------------------------------------#



               #----- Define the grid information for the 3-D plot. -----------------------#
               xlabels = pretty.log(xlimit,n=5)
               ylabels = pretty(ylimit,n=5)
               zlabels = if(plog){pretty.log(zlimit,n=5)}else{pretty(zlimit,n=5)}
               xat     = log(xlabels)
               yat     = ylabels
               zat     = if(plog){log(zlabels)}else{zlabels}
               xlimit  = range(x=xat)
               ylimit  = range(x=yat)
               zlimit  = range(x=zat)
               xfloor  = seq(from=xlimit[1],to=xlimit[2],length.out=16)
               yfloor  = seq(from=ylimit[1],to=ylimit[2],length.out=16)
               zfloor  = matrix(floor3d,nrow=length(xfloor),ncol=length(yfloor))
               #---------------------------------------------------------------------------#



               #----- Expand the lines to make the lollipops. -----------------------------#
               ncohnow  = length(varww)
               ageww    = rep(ageww,each=3)
               dbhww    = rep(dbhww,each=3)
               pftww    = rep(pftww,each=3)
               varww    = as.vector(rbind(rep(floor3d,times=ncohnow)
                                         ,varco[[ww]]
                                         ,rep(NA,times=ncohnow)))
               xww      = log(ageww)
               yww      = dbhww
               zww      = if(plog){log(varww)}else{varww}
               pchww    = rep(c(NA,16,NA),times=ncohnow)
               cexww    = rep(cexww,each=3)
               colww    = pft$colour[pftww]

               pftin   = sort(unique(cohort$pft[[ww]]))
               colleg  = pft$colour[pftin]
               pftleg  = pft$name  [pftin]
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #   Plot annotation.                                                        #
               #---------------------------------------------------------------------------#
               letitre = paste0(description," - ",lieu,
                               "\n Time :",month.abb[mm],"/",thisyear)
               lexlab  = desc.unit(desc="Gap age",unit=untab$yr)
               leylab  = desc.unit(desc="DBH",unit=untab$cm)
               lezlab  = desc.unit(desc=description,unit=unit)
               #---------------------------------------------------------------------------#


               #----- Loop over output formats. -------------------------------------------#
               for (o in sequence(nout)){
                  #----- Open file. -------------------------------------------------------#
                  fichier = file.path(outdir
                                     ,paste0(vnam,"-",thisyear,"-",cmonth,"-",suffix
                                            ,".",outform[o])
                                     )#end file.path
                  dummy   = open.plot( fichier = fichier
                                     , outform = outform[o]
                                     , size    = eysize
                                     , ptsz    = ptsz
                                     , depth   = depth
                                     )#end open.plot
                  #------------------------------------------------------------------------#


                  #----- Split the domain into 2. -----------------------------------------#
                  par(par.user)
                  layout(mat=rbind(2,1),heights=c(1.-f.leg,f.leg))
                  #------------------------------------------------------------------------#


                  #------------------------------------------------------------------------#
                  #     Plot legend.                                                       #
                  #------------------------------------------------------------------------#
                  par(mar=c(0.1,0.1,0.1,0.1))
                  plot.new()
                  plot.window(xlim=c(0,1),ylim=c(0,1))
                  legend( x      = "center"
                        , inset  = 0.0
                        , legend = pftleg
                        , fill   = colleg
                        , ncol   = min(4,pretty.box(length(pftleg))$ncol)
                        , title  = expression(bold("Plant functional type"))
                        , cex    = cex.ptsz
                        , xpd    = TRUE
                        , bty    = "n"
                        )#end legend
                  #------------------------------------------------------------------------#


                  #------------------------------------------------------------------------#
                  #     Plot the 3-D plot.                                                 #
                  #------------------------------------------------------------------------#
                  par(mar=c(1.1,1.1,4.1,1.1))
                  pout = perspx( x         = xfloor
                               , y         = yfloor
                               , z         = zfloor
                               , xlim      = xlimit
                               , ylim      = ylimit
                               , zlim      = zlimit
                               , theta     = theta
                               , phi       = phi
                               , col       = gcol
                               , expand    = expz
                               , ticktype  = "detailed"
                               , border    = NA
                               , shade     = shade
                               , ltheta    = ltheta
                               , main      = letitre
                               , cex.main  = 0.8*cex.ptsz
                               , axes      = FALSE
                               )#end perspx
                  #----- Add axes. --------------------------------------------------------#
                  paxis3d(edge="X--",pmat=pout,at=xat,cex=0.9*cex.ptsz,labels=xlabels)
                  paxis3d(edge="Y--",pmat=pout,at=yat,cex=0.9*cex.ptsz,labels=ylabels)
                  paxis3d(edge="Z-+",pmat=pout,at=zat,cex=0.9*cex.ptsz,labels=zlabels)
                  mtext3d(edge="X--",pmat=pout,labels=lexlab,cex=cex.ptsz,srt=theta+90)
                  mtext3d(edge="Y--",pmat=pout,labels=leylab,cex=cex.ptsz,srt=theta)
                  mtext3d(edge="Z-+",pmat=pout,labels=lezlab,cex=cex.ptsz,srt=-75)
                  #------------------------------------------------------------------------#


                  #----- Add the cohorts. -------------------------------------------------#
                  lines (trans3d(x=xww,y=yww,z=zww,pmat=pout),type="l",col=grey.fg,lwd=2)
                  points(trans3d(x=xww,y=yww,z=zww,pmat=pout),type="p",pch=pchww
                        ,col=colww,cex=cexww)
                  #------------------------------------------------------------------------#



                  #----- Close the device. ------------------------------------------------#
                  dummy = close.plot(outform=outform[o])
                  #------------------------------------------------------------------------#
               }#end for outform
               #---------------------------------------------------------------------------#
            }#end if is.na(varww)
            #------------------------------------------------------------------------------#
         }#end for nameco
         #---------------------------------------------------------------------------------#
      }#end if
      #------------------------------------------------------------------------------------#
   }#end for npsas
   #---------------------------------------------------------------------------------------#
}#end for places
#==========================================================================================#
#==========================================================================================#
