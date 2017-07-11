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
obsrfile = file.path(srcdir,"LBA_MIP.v8.RData")
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
   outpref = file.path(outmain,"yearly")
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
   if (! file.exists(outpref)) dir.create(outpref)
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
   cat ("    - Finding the annual statistics for multi-dimensional variables...","\n")
   cat ("      * Aggregating the annual mean of PFT-DBH variables...","\n")
   for (vname in names(szpft)){
      szpft[[vname]] = qapply(X=szpft[[vname]],INDEX=yfac,DIM=1,FUN=mean,na.rm=TRUE)
   }#end for
   #----- LU arrays.   The "+1" column contains the total. --------------------------------#
   cat ("      * Aggregating the annual mean of LU variables...","\n")
   for (vname in names(lu)){
      lu   [[vname]] = qapply(X=lu   [[vname]],INDEX=yfac,DIM=1,FUN=mean,na.rm=TRUE)
   }#end for
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Here we find the monthly means for month, then compute the standard deviation.   #
   #---------------------------------------------------------------------------------------#
   cat ("    - Finding the monthly and annual means...","\n")
   cat ("      * Aggregating the monthly mean and standard deviation...","\n")
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
   cat ("    - Aggregating the annual mean and std. dev. of the diurnal cycle...","\n")
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
   mplot  = list()
   mplot$levels = montmont
   mplot$labels = capwords(mon2mmm(montmont))
   mplot$n      = 12
   mplot$scale  = "months"
   mplot$padj   = rep(0,times=mplot$n)
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
         cat("      +",description,"time series for all PFTs...","\n")

         #----- Load variable -------------------------------------------------------------#
         if (vnam %in% names(szpft)){
            thisvar = szpft[[vnam]][,ndbh+1,]
            if (plog){
               #----- Eliminate non-positive values in case it is a log plot. -------------#
               badlog          = is.finite(thisvar) & thisvar <= 0
               thisvar[badlog] = NA
            }#end if
         }else{
            thisvar = matrix(NA,ncol=npft+1,nrow=nyears)
         }#end if
         #---------------------------------------------------------------------------------#



         #----- Loop over output formats. -------------------------------------------------#
         for (o in sequence(nout)){
            fichier = file.path(outdir,paste0(vnam,"-",suffix,".",outform[o]))
            if(outform[o] %in% "x11"){
               X11(width=size$width,height=size$height,pointsize=ptsz)
            }else if(outform[o] %in% "quartz"){
               quartz(width=size$width,height=size$height,pointsize=ptsz)
            }else if(outform[o] %in% "png"){
               png(filename=fichier,width=size$width*depth,height=size$height*depth
                  ,pointsize=ptsz,res=depth,bg="transparent")
            }else if(outform[o] %in% "tif"){
               tiff(filename=fichier,width=size$width*depth,height=size$height*depth
                   ,pointsize=ptsz,res=depth,bg="transparent",compression="lzw")
            }else if(outform[o] %in% "eps"){
               postscript(file=fichier,width=size$width,height=size$height
                         ,pointsize=ptsz,paper=size$paper)
            }else if(outform[o] %in% "pdf"){
               pdf(file=fichier,onefile=FALSE
                  ,width=size$width,height=size$height,pointsize=ptsz,paper=size$paper)
            }#end if


            #------------------------------------------------------------------------------#
            #     Find the limit, make some room for the legend, and in case the field is  #
            # a constant, nudge the limits so the plot command will not complain.          #
            #------------------------------------------------------------------------------#
            xlimit = pretty.xylim(u = datum$toyear    ,fracexp=0.0,is.log=FALSE)
            ylimit = pretty.xylim(u = thisvar[,selpft],fracexp=0.0,is.log=plog )
            if (plog){
               xylog    = "y"
               ydrought = c( exp(sqrt(ylimit[1]^3/ylimit[2]))
                           , exp(sqrt(ylimit[2]^3/ylimit[1]))
                           )#end c
            }else{
               xylog    = ""
               ydrought = c( ylimit[1] - 0.5 * diff(ylimit),ylimit[2] + 0.5 * diff(ylimit) )
            }#end if
            #------------------------------------------------------------------------------#


            #----- Plot settings. ---------------------------------------------------------#
            letitre = paste(description,lieu,sep=" - ")
            ley     = desc.unit(desc=description,unit=unit)
            cols    = pft$colour[selpft]
            legs    = pft$name  [selpft]
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #     Split the plot into two windows.                                         #
            #------------------------------------------------------------------------------#
            par(par.user)
            layout(mat=rbind(2,1),heights=c(5,1))
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
            axis(side=1)
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
               abline(v=axTicks(side=1),h=axTicks(side=2),col=grid.colour,lty="solid")
            }#end if
            #----- Plot lines. ------------------------------------------------------------#
            for (n in sequence(npft+1)){
               if (selpft[n]){
                  lines(datum$toyear,thisvar[,n],type="l",col=pft$colour[n],lwd=lwidth)
               }#end if
            }#end for
            #------------------------------------------------------------------------------#


            #----- Close the device. ------------------------------------------------------#
            if (outform[o] %in% c("x11","quartz")){
               locator(n=1)
               dev.off()
            }else{
               dev.off()
            }#end if
            dummy=clean.tmp()
            #------------------------------------------------------------------------------#
         } #end for outform
      }#end if (tseragbpft)
   } #end for tseries
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #      Time series by DBH, by PFT.                                                      #
   #---------------------------------------------------------------------------------------#
   #----- Find the PFTs to plot. ----------------------------------------------------------#
   pftuse  = which(apply(X=is.na(szpft$nplant),MARGIN=3,FUN=sum,na.rm=TRUE) == 0.)
   pftuse  = pftuse[pftuse != (npft+1)]
   for (v in sequence(ntspftdbh)){
      thistspftdbh   = tspftdbh[[v]]
      vnam        = thistspftdbh$vnam
      description = thistspftdbh$desc
      unit        = thistspftdbh$e.unit
      plog        = thistspftdbh$plog
      plotit      = thistspftdbh$pftdbh
      
      #----- Load variable ----------------------------------------------------------------#
      if (vnam %in% names(szpft)){
         thisvar = szpft[[vnam]]
         if (plog){
            xylog="y"
            badlog = is.finite(thisvar) & thisvar <= 0
            thisvar[badlog] = NA
         }else{
            xylog=""
         }#end if
      }else{
         thisvar = array(NA,dim=c(nyears,ndbh+1,npft+1))
      }#end if
      #----- Check whether the user wants to have this variable plotted. ------------------#
      if (plotit && length(pftuse) > 0 && any(is.finite(thisvar))){

         #---------------------------------------------------------------------------------#
         #    Check whether the time series directory exists.  If not, create it.          #
         #---------------------------------------------------------------------------------#
         outdir = file.path(outpref,"tsdbh")
         if (! file.exists(outdir)) dir.create(outdir)
         outvar = file.path(outdir,vnam)
         if (! file.exists(outvar)) dir.create(outvar)
         #---------------------------------------------------------------------------------#

         cat("      +",description,"time series for DBH class...","\n")


         #---------------------------------------------------------------------------------#
         #     Find the limit, make some room for the legend, and in case the field is a   #
         # constant, nudge the limits so the plot command will not complain.               #
         #---------------------------------------------------------------------------------#
         xlimit = pretty.xylim(u=datum$toyear     ,fracexp=0.0,is.log=FALSE)
         ylimit = pretty.xylim(u=thisvar[,,pftuse],fracexp=0.0,is.log=plog)
         if (plog){
            xylog    = "y"
            ydrought = c( exp(sqrt(ylimit[1]^3/ylimit[2]))
                        , exp(sqrt(ylimit[2]^3/ylimit[1]))
                        )#end c
         }else{
            xylog    = ""
            ydrought = c( ylimit[1] - 0.5 * diff(ylimit),ylimit[2] + 0.5 * diff(ylimit) )
         }#end if
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #       Loop over plant functional types.                                         #
         #---------------------------------------------------------------------------------#
         for (p in pftuse){
            pftlab = paste0("pft-",sprintf("%2.2i",p))
            cat("        - ",pft$name[p],"\n")


            #----- Loop over output formats. ----------------------------------------------#
            for (o in sequence(nout)){
               #----- Open file. ----------------------------------------------------------#
               fichier = file.path( outvar
                                  , paste0(vnam,"-",pftlab,"-",suffix,".",outform[o])
                                  )#end file.path
               if (outform[o] %in% "x11"){
                  X11(width=size$width,height=size$height,pointsize=ptsz)
               }else if (outform[o] %in% "quartz"){
                  quartz(width=size$width,height=size$height,pointsize=ptsz)
               }else if (outform[o] %in% "png"){
                  png(filename=fichier,width=size$width*depth,height=size$height*depth
                     ,pointsize=ptsz,res=depth,bg="transparent")
               }else if (outform[o] %in% "tif"){
                  tiff(filename=fichier,width=size$width*depth,height=size$height*depth
                      ,pointsize=ptsz,res=depth,bg="transparent",compression="lzw")
               }else if (outform[o] %in% "eps"){
                  postscript(file=fichier,width=size$width,height=size$height
                            ,pointsize=ptsz,paper=size$paper)
               }else if (outform[o] %in% "pdf"){
                  pdf(file=fichier,onefile=FALSE
                     ,width=size$width,height=size$height,pointsize=ptsz,paper=size$paper)
               }#end if
               #---------------------------------------------------------------------------#



               #-----  Plot annotation. ---------------------------------------------------#
               letitre = paste(description,pft$name[p],lieu,sep=" - ")
               ley     = desc.unit(desc=description,unit=unit)
               #---------------------------------------------------------------------------#


               #---------------------------------------------------------------------------#
               #     Split the plot into two windows.                                      #
               #---------------------------------------------------------------------------#
               par(par.user)
               layout(mat=rbind(2,1),heights=c(5,1))
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
               axis(side=1)
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
               #----- Plot grid. ----------------------------------------------------------#
               if (plotgrid){ 
                  abline(v=axTicks(side=1),h=axTicks(side=2),col=grid.colour,lty="solid")
               }#end if
               #----- Plot lines. ---------------------------------------------------------#
               for (d in seq(from=1,to=ndbh+1,by=1)){
                  lines(datum$toyear,thisvar[,d,p],type="l",col=dbhcols[d],lwd=lwidth)
               }#end for
               #---------------------------------------------------------------------------#


               #----- Close the device. ---------------------------------------------------#
               if (outform[o] %in% c("x11","quartz")){
                  locator(n=1)
                  dev.off()
               }else{
                  dev.off()
               }#end if
               dummy=clean.tmp()
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
   cat("    + Year-by-year comparisons of monthly means...","\n")
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

      plotit       = ( plotit && vname %in% names(emean) && vname %in% names(mmean)
                              && plot.ycomp )

      if (plotit){
         #---------------------------------------------------------------------------------#
         #    Copy the observations to a scratch variable.                                 #
         #---------------------------------------------------------------------------------#
         thisvar     = emean [[vname]]
         thismean    = mmean [[vname]]
         if (length(msdev[[vname]]) == 0){
            thissdev = 0. * thismean
         }else{
            thissdev = msdev[[vname]]
         }#end if
         mod.x       = montmont
         mod.ylow    = thismean - thissdev
         mod.yhigh   = thismean + thissdev
         mod.x.poly  = c(mod.x,rev(mod.x))
         mod.y.poly  = c(mod.ylow,rev(mod.yhigh))
         mod.keep    = is.finite(mod.y.poly)
         mod.x.poly  = mod.x.poly[mod.keep]
         mod.y.poly  = mod.y.poly[mod.keep]
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #    Check whether the time series directory exists.  If not, create it.          #
         #---------------------------------------------------------------------------------#
         outdir   = file.path(outpref,"ycomp")
         outvar   = file.path(outdir,vname)
         if (! file.exists(outdir)) dir.create(outdir)
         if (! file.exists(outvar)) dir.create(outvar)
         cat("      - ",description,"comparison...","\n")
         #---------------------------------------------------------------------------------#



         #----- Find the plot range. ------------------------------------------------------#
         if (plotsd){
            ylimit    = range(c(mod.ylow,mod.yhigh,thisvar),na.rm=TRUE)
         }else{
            ylimit    = range(thisvar,na.rm=TRUE)
         }#end if
         #----- Expand the upper range in so the legend doesn't hide things. --------------#
         ylimit = pretty.xylim(u=ylimit,fracexp=0.0,is.log=FALSE)
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Loop over all years, and make one plot per year.                            #
         #---------------------------------------------------------------------------------#
         for (y in sequence(nyears)){
            #----- Retrieve the year and the variable for this year. ----------------------#
            year.now = datum$toyear[y]
            cyear    = sprintf("%4.4i",year.now)
            var.year = thisvar[yfac == year.now]
            #------------------------------------------------------------------------------#



            #----- Load variable ----------------------------------------------------------#
            letitre = paste0(description," - ",lieu,"\n","Monthly mean - ",cyear)
            ley     = desc.unit(desc=description,unit=unit)
            #------------------------------------------------------------------------------#


            #----- Loop over formats. -----------------------------------------------------#
            for (o in sequence(nout)){
               fichier = file.path(outvar,paste0(vname,"-",cyear,".",outform[o]))
               if (outform[o] %in% "x11"){
                  X11(width=size$width,height=size$height,pointsize=ptsz)
               }else if (outform[o] %in% "quartz"){
                  quartz(width=size$width,height=size$height,pointsize=ptsz)
               }else if (outform[o] %in% "png"){
                  png(filename=fichier,width=size$width*depth,height=size$height*depth
                     ,pointsize=ptsz,res=depth,bg="transparent")
               }else if (outform[o] %in% "tif"){
                  tiff(filename=fichier,width=size$width*depth,height=size$height*depth
                      ,pointsize=ptsz,res=depth,bg="transparent",compression="lzw")
               }else if (outform[o] %in% "eps"){
                  postscript(file=fichier,width=size$width,height=size$height
                            ,pointsize=ptsz,paper=size$paper)
               }else if (outform[o] %in% "pdf"){
                  pdf(file=fichier,onefile=FALSE
                     ,width=size$width,height=size$height,pointsize=ptsz,paper=size$paper)
               }#end if



               #----- Split plot into two windows. ----------------------------------------#
               par(par.user)
               layout(mat=rbind(2,1),heights=c(5,1))
               #---------------------------------------------------------------------------#



               #------ First plot: the legend. --------------------------------------------#
               par(mar=c(0.1,4.6,0.1,2.1))
               plot.new()
               plot.window(xlim=c(0,1),ylim=c(0,1))
               if (plotsd){
                  legend( x       = "bottom"
                        , inset   = 0.0
                        , legend  = c(cyear,paste0("Mean: ",yeara,"-",yearz))
                        , fill    = errcolours
                        , angle   = angle
                        , density = dens
                        , lwd     = llwd
                        , col     = lcolours
                        , bg      = background
                        , title   = expression(bold("Shaded areas = 1 SD"))
                        , cex     = cex.ptsz
                        , pch     = 16
                        , xpd     = TRUE
                        , bty     = "n"
                        )#end legend
               }else{
                  legend( x      = "bottom"
                        , inset  = 0.0
                        , legend = c(cyear,paste0("Mean: ",yeara,"-",yearz))
                        , col    = lcolours
                        , lwd    = llwd
                        , cex    = cex.ptsz
                        , pch    = 16
                        , xpd    = TRUE
                        , bty    = "n"
                        )#end legend
               }#end if
               #---------------------------------------------------------------------------#




               #------ Second plot: the comparison. ---------------------------------------#
               par(mar=c(4.1,4.6,4.1,2.1))
               plot.new()
               plot.window(xlim=range(montmont),ylim=ylimit,log=plog)
               if (plotgrid){ 
                  abline(v=mplot$levels,h=axTicks(side=2),col=grid.colour,lty="solid")
               }#end if
               if (plotsd){
                  polygon(x=mod.x.poly,y=mod.y.poly,col=errcolours[2],angle=angle[2]
                         ,density=dens[1],lty="solid",lwd=shwd[1])
               }#end if
               points(x=montmont,y=var.year,col=lcolours[1],lwd=llwd[1],type=ltype
                     ,pch=16,cex=1.0)
               points(x=montmont,y=thismean,col=lcolours[2],lwd=llwd[2],type=ltype
                     ,pch=16,cex=1.0)
               axis(side=1,at=mplot$levels,labels=mplot$labels,padj=mplot$padj)
               axis(side=2,las=1)
               title(main=letitre,xlab="Time",ylab=ley,cex.main=cex.main)
               box()
               #---------------------------------------------------------------------------#




               #----- Close plotting window. ----------------------------------------------#
               if (outform[o] %in% c("x11","quartz")){
                  locator(n=1)
                  dev.off()
               }else{
                  dev.off()
               }#end if
               dummy=clean.tmp()
               #---------------------------------------------------------------------------#
            }#end for outform
            #------------------------------------------------------------------------------#
         }#end for years
         #---------------------------------------------------------------------------------#
      }#end if plotit
      #------------------------------------------------------------------------------------#
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
         cat("      +",description,"time series for all LUs...","\n")



         #----- Load variable -------------------------------------------------------------#
         thisvar = lu[[vnam]]
         if (plog){
            #----- Eliminate non-positive values in case it is a log plot. ----------------#
            thisvar[thisvar <= 0] = NA
         }#end if
         #---------------------------------------------------------------------------------#

         #----- Loop over output formats. -------------------------------------------------#
         for (o in sequence(nout)){
            fichier = file.path(outdir,paste0(vnam,"-",suffix,".",outform[o]))
            if (outform[o] %in% "x11"){
               X11(width=size$width,height=size$height,pointsize=ptsz)
            }else if (outform[o] %in% "quartz"){
               quartz(width=size$width,height=size$height,pointsize=ptsz)
            }else if(outform[o] %in% "png"){
               png(filename=fichier,width=size$width*depth,height=size$height*depth
                  ,pointsize=ptsz,res=depth,bg="transparent")
            }else if(outform[o] %in% "tif"){
               tiff(filename=fichier,width=size$width*depth,height=size$height*depth
                   ,pointsize=ptsz,res=depth,bg="transparent",compression="lzw")
            }else if(outform[o] %in% "eps"){
               postscript(file=fichier,width=size$width,height=size$height
                         ,pointsize=ptsz,paper=size$paper)
            }else if(outform[o] %in% "pdf"){
               pdf(file=fichier,onefile=FALSE
                  ,width=size$width,height=size$height,pointsize=ptsz,paper=size$paper)
            }#end if


            #------------------------------------------------------------------------------#
            #     Find the limit, make some room for the legend, and in case the field is  #
            # a constant, nudge the limits so the plot command will not complain.          #
            #------------------------------------------------------------------------------#
            xlimit = pretty.xylim(u = datum$toyear ,fracexp=0.0,is.log=FALSE)
            ylimit = pretty.xylim(u=thisvar[,sellu],fracexp=0.0,is.log=plog)
            if (plog){
               xylog    = "y"
               ydrought = c( exp(sqrt(ylimit[1]^3/ylimit[2]))
                           , exp(sqrt(ylimit[2]^3/ylimit[1]))
                           )#end c
            }else{
               xylog    = ""
               ydrought = c( ylimit[1] - 0.5 * diff(ylimit),ylimit[2] + 0.5 * diff(ylimit) )
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
            layout(mat=rbind(2,1),heights=c(5,1))
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
            axis(side=1)
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
               abline(v=axTicks(side=1),h=axTicks(side=2),col=grid.colour,lty="solid")
            }#end if
            #----- Plot lines. ------------------------------------------------------------#
            for (n in sequence(nlu+1)){
               if (sellu[n]){
                  lines(datum$toyear,thisvar[,n],type="l",col=lucols[n],lwd=lwidth)
               }#end if
            }#end for
            #------------------------------------------------------------------------------#


            #----- Close the device. ------------------------------------------------------#
            if (outform[o] %in% c("x11","quartz")){
               locator(n=1)
               dev.off()
            }else{
               dev.off()
            }#end if
            dummy=clean.tmp()
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
      cat("      + Disturbance rate time series for all disturbances...","\n")
      for (o in sequence(nout)){
         fichier = file.path(outpref,paste0("disturb-",suffix,".",outform[o]))
         if (outform[o] %in% "x11"){
            X11(width=size$width,height=size$height,pointsize=ptsz)
         }else if (outform[o] %in% "quartz"){
            quartz(width=size$width,height=size$height,pointsize=ptsz)
         }else if(outform[o] %in% "png"){
            png(filename=fichier,width=size$width*depth,height=size$height*depth
               ,pointsize=ptsz,res=depth,bg="transparent")
         }else if(outform[o] %in% "tif"){
            tiff(filename=fichier,width=size$width*depth,height=size$height*depth
                ,pointsize=ptsz,res=depth,bg="transparent",compression="lzw")
         }else if(outform[o] %in% "eps"){
            postscript(file=fichier,width=size$width,height=size$height
                      ,pointsize=ptsz,paper=size$paper)
         }else if(outform[o] %in% "pdf"){
            pdf(file=fichier,onefile=FALSE
               ,width=size$width,height=size$height,pointsize=ptsz,paper=size$paper)
         }#end if

         #---------------------------------------------------------------------------------#
         #     Find the limit, make some room for the legend, and in case the field is a   #
         #  constant, nudge the limits so the plot command will not complain.              #
         #---------------------------------------------------------------------------------#
         xlimit   = pretty.xylim(u=datum$toyear,fracexp=0.0,is.log=FALSE)
         ylimit   = NULL
         n        = 0
         mylucols = NULL
         mylulegs = NULL
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
         layout(mat=rbind(2,1),heights=c(5,1))
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
         axis(side=1)
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
            abline(v=axTicks(side=1),h=axTicks(side=2),col=grid.colour,lty="solid")
         }#end if
         #----- Plot lines. ---------------------------------------------------------------#
         n = 0
         for (jlu in sequence(nlu)){
            for (ilu in sequence(nlu)){
               n = n + 1
               if (seldist[ilu,jlu]){
                  lines(datum$toyear,lu$dist[,ilu,jlu],type="l"
                       ,col=distcols[n],lwd=lwidth)
               }#end if
            }#end for
         }#end for
         #---------------------------------------------------------------------------------#


         #----- Close the device. ---------------------------------------------------------#
         if (outform[o] %in% c("x11","quartz")){
            locator(n=1)
            dev.off()
         }else{
            dev.off()
         }#end if
         dummy=clean.tmp()
         #---------------------------------------------------------------------------------#
      } #end for outform
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #   Plot the time series diagrams showing annual means.                                 #
   #---------------------------------------------------------------------------------------#
   cat("      * Plot time series of groups of variables...","\n")
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
      plotit       = themenow$ymean
      ylimit.fix   = themenow$ymean.lim

      if (plotit){

         #---------------------------------------------------------------------------------#
         #    Check whether the time series directory exists.  If not, create it.          #
         #---------------------------------------------------------------------------------#
         outdir = file.path(outpref,"theme_ymean")
         if (! file.exists(outdir)) dir.create(outdir)
         cat("      +",group,"time series...","\n")


         #----- Define the number of layers. ----------------------------------------------#
         nlayers   = length(vnames)
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Find the limit, make some room for the legend, and in case the field is a   #
         # constant, nudge the limits so the plot command will not complain.               #
         #---------------------------------------------------------------------------------#
         xlimit   = pretty.xylim(u=datum$toyear,fracexp=0.0,is.log=FALSE)
         if (any(! is.finite(ylimit.fix))){
            ylimit    = NULL
            for (l in sequence(nlayers)){
               thisvar = ymean[[vnames[l]]]
               ylimit  = range(c(ylimit,thisvar),na.rm=TRUE)
            }#end for
            ylimit = pretty.xylim(u=ylimit,fracexp=0.0,is.log=plog)
         }else{
            ylimit = ylimit.fix
         }#end if
         if (plog) {
            xylog    = "y"
            ydrought = c( exp(sqrt(ylimit[1]^3/ylimit[2]))
                        , exp(sqrt(ylimit[2]^3/ylimit[1]))
                        )#end c
         }else{
            xylog    = ""
            ydrought = c(ylimit[1]-0.5*diff(ylimit),ylimit[2]+0.5*diff(ylimit))
         }#end if
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Check if the directory exists.  If not, create it.                          #
         #---------------------------------------------------------------------------------#

         #----- Loop over formats. --------------------------------------------------------#
         for (o in sequence(nout)){
            #------ Open file. ------------------------------------------------------------#
            fichier = file.path(outdir,paste0(prefix,"-",suffix,".",outform[o]))
            if (outform[o] %in% "x11"){
               X11(width=size$width,height=size$height,pointsize=ptsz)
            }else if (outform[o] %in% "quartz"){
               quartz(width=size$width,height=size$height,pointsize=ptsz)
            }else if (outform[o] %in% "png"){
               png(filename=fichier,width=size$width*depth,height=size$height*depth
                  ,pointsize=ptsz,res=depth,bg="transparent")
            }else if (outform[o] %in% "tiff"){
               tiff(filename=fichier,width=size$width*depth,height=size$height*depth
                  ,pointsize=ptsz,res=depth,bg="transparent",compression="lzw")
            }else if (outform[o] %in% "eps"){
               postscript(file=fichier,width=size$width,height=size$height
                         ,pointsize=ptsz,paper=size$paper)
            }else if (outform[o] %in% "pdf"){
               pdf(file=fichier,onefile=FALSE
                  ,width=size$width,height=size$height,pointsize=ptsz,paper=size$paper)
            }#end if
            #------------------------------------------------------------------------------#



            #----- Plot settings. ---------------------------------------------------------#
            letitre = paste0(" Time series: ",group,"\n",lieu)
            ley     = desc.unit(desc=group,unit=unit)
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #     Split the plot into two windows.                                         #
            #------------------------------------------------------------------------------#
            par(par.user)
            layout(mat=rbind(2,1),heights=c(5,1))
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
            axis(side=1)
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
               abline(v=axTicks(side=1),h=axTicks(side=2),col=grid.colour,lty="solid")
            }#end if
            #----- Plot lines. ------------------------------------------------------------#
            for (l in sequence(nlayers)){
               thisvar = ymean[[vnames[l]]]
               points(x=datum$toyear,y=thisvar,col=lcolours[l],lwd=llwd[l],type=ltype
                     ,pch=16,cex=0.8)
            }#end for
            #------------------------------------------------------------------------------#


            #----- Close the device. ------------------------------------------------------#
            if (outform[o] %in% c("x11","quartz")){
               locator(n=1)
               dev.off()
            }else{
               dev.off()
            }#end if
            dummy=clean.tmp()
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
   cat("      * Plot mean diel for groups of variables...","\n")
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
      ylimit.fix   = themenow$qmean.lim
      plotit       = themenow$qmean && plot.ycomp
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
         cat("      +",group," diurnal cycle...","\n")


         #----- Define the number of layers. ----------------------------------------------#
         nlayers   = length(vnames)
         xlimit    = range(thisday)
         if (any(! is.finite(ylimit.fix))){
            ylimit    = NULL
            for (l in sequence(nlayers)) ylimit  = c(ylimit,umean[[vnames[l]]])
            ylimit = pretty.xylim(u=ylimit,fracexp=0.0,is.log=plog)
         }else{
            ylimit = ylimit.fix
         }#end if
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #      Loop over all months.                                                      #
         #---------------------------------------------------------------------------------#
         yplot = as.numeric(dimnames(umean[[vnames[1]]])[[1]])
         for (yy in seq_along(yplot)){
            cyear    = sprintf("%4.4i",yplot[yy])

            #------------------------------------------------------------------------------#
            #     Check if the directory exists.  If not, create it.                       #
            #------------------------------------------------------------------------------#

            #----- Loop over formats. -----------------------------------------------------#
            for (o in sequence(nout)){
               #------ Open file. ---------------------------------------------------------#
               fichier = file.path( outtheme
                                  , paste0(prefix,"-",cyear,"-",suffix,".",outform[o])
                                  )#end file.path
               if (outform[o] %in% "x11"){
                  X11(width=size$width,height=size$height,pointsize=ptsz)
               }else if (outform[o] %in% "quartz"){
                  quartz(width=size$width,height=size$height,pointsize=ptsz)
               }else if (outform[o] %in% "png"){
                  png(filename=fichier,width=size$width*depth,height=size$height*depth
                     ,pointsize=ptsz,res=depth,bg="transparent")
               }else if (outform[o] %in% "tif"){
                  tiff(filename=fichier,width=size$width*depth,height=size$height*depth
                      ,pointsize=ptsz,res=depth,bg="transparent",compression="lzw")
               }else if (outform[o] %in% "eps"){
                  postscript(file=fichier,width=size$width,height=size$height
                            ,pointsize=ptsz,paper=size$paper)
               }else if (outform[o] %in% "pdf"){
                  pdf(file=fichier,onefile=FALSE
                     ,width=size$width,height=size$height,pointsize=ptsz,paper=size$paper)
               }#end if
               #---------------------------------------------------------------------------#



               #----- Plot settings. ------------------------------------------------------#
               letitre = paste0(group," - ",lieu,"\n","Mean diurnal cycle - ",cyear)
               ley     = desc.unit(desc=group,unit=unit)
               #---------------------------------------------------------------------------#


               #---------------------------------------------------------------------------#
               #     Split the plot into two windows.                                      #
               #---------------------------------------------------------------------------#
               par(par.user)
               layout(mat=rbind(2,1),heights=c(5,1))
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
               title(main=letitre,xlab="Year",ylab=ley,cex.main=0.7,log=xylog)
               #----- Plot grid. ----------------------------------------------------------#
               if (plotgrid){ 
                  abline(v=uplot$levels,h=axTicks(side=2),col=grid.colour,lty="solid")
               }#end if
               #----- Plot lines. ---------------------------------------------------------#
               for (l in sequence(nlayers)){
                  thisvar = umean[[vnames[l]]]
                  thisvar = cbind(thisvar[,ndcycle],thisvar)
                  points(x=thisday,y=thisvar[yy,],col=lcolours[l]
                        ,lwd=llwd[l],type=ltype,pch=16)
               }#end for
               #---------------------------------------------------------------------------#


               #----- Close the device. ---------------------------------------------------#
               if (outform[o] %in% c("x11","quartz")){
                  locator(n=1)
                  dev.off()
               }else{
                  dev.off()
               }#end if
               dummy=clean.tmp()
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
      thisclim    = soilplot[[v]]
      vnam        = thisclim$vnam
      description = thisclim$desc
      unit        = thisclim$unit
      vcscheme    = thisclim$csch
      pnlog       = thisclim$pnlog
      plotit      = thisclim$ymean

      if (plotit){

         #---------------------------------------------------------------------------------#
         #     Check if the directory exists.  If not, create it.                          #
         #---------------------------------------------------------------------------------#
         outdir  =  file.path(outpref,"soil_ymean")
         if (! file.exists(outdir)) dir.create(outdir)
         cat("      + Climatology profile of ",description,"...","\n")

         #----- Find the number of rows and columns, and the axes. ------------------------#
         monaxis  = sort(unique(datum$year))
         soilaxis = slz
         nmon     = length(monaxis)
         nsoil    = nzg

         #----- Convert the vector data into an array. ------------------------------------#
         vararr  = ymean[[vnam]]

         #----- Copy the first and the last year to make the edges buffered. --------------#
         first    = vararr[1,]
         first    = c(first,first[nzg],first[nzg])

         last     = vararr[nyears,]
         last     = c(last[1],last[1],last)
         #---------------------------------------------------------------------------------#



         #----- Bind first and last year to the array, to make the edges buffered. --------#
         varbuff  = cbind(vararr[,1],vararr,vararr[,nzg])
         varbuff  = rbind(last,varbuff,first)
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #   Expand the month and year axes.                                               #
         #---------------------------------------------------------------------------------#
         yearaxis = c(min(datum$toyear)-1,datum$toyear,max(datum$toyear)+1)
         soilaxis = -log(-1.0 * c( slz[1]*(slz[1]/slz[2])
                                 , soilaxis
                                 , slz[nzg]*(slz[nzg]/slz[nzg-1]) ))

         if (pnlog){
            vrange  = range(varbuff,na.rm=TRUE)
            vlevels = pretty.log(x=vrange,n=ncolshov)
            vnlev   = length(vlevels)
         }else{
            vrange  = range(varbuff,na.rm=TRUE)
            vlevels = pretty(x=vrange,n=ncolshov)
            vnlev   = length(vlevels)
         }#end if

         #----- Loop over formats. --------------------------------------------------------#
         for (o in sequence(nout)){
            fichier = file.path(outdir,paste0(vnam,"-",suffix,".",outform[o]))
            if (outform[o] %in% "x11"){
               X11(width=size$width,height=size$height,pointsize=ptsz)
            }else if (outform[o] %in% "quartz"){
               quartz(width=size$width,height=size$height,pointsize=ptsz)
            }else if (outform[o] %in% "png"){
               png(filename=fichier,width=size$width*depth,height=size$height*depth
                  ,pointsize=ptsz,res=depth,bg="transparent")
            }else if (outform[o] %in% "tif"){
               tiff(filename=fichier,width=size$width*depth,height=size$height*depth
                   ,pointsize=ptsz,res=depth,bg="transparent",compression="lzw")
            }else if (outform[o] %in% "eps"){
               postscript(file=fichier,width=size$width,height=size$height
                         ,pointsize=ptsz,paper=size$paper)
            }else if (outform[o] %in% "pdf"){
               pdf(file=fichier,onefile=FALSE
                  ,width=size$width,height=size$height,pointsize=ptsz,paper=size$paper)
            }#end if

            letitre = paste0(description,"\n",lieu)
            ley     = desc.unit(desc="Soil depth",unit=untab$m)
            lacle   = desc.unit(desc=NULL,unit=unit)
            par(par.user)
            sombreado(x=yearaxis,y=soilaxis,z=varbuff,levels=vlevels,nlevels=vnlev
                     ,color.palette=get(vcscheme)
                     ,plot.title=title(main=letitre,xlab="Month",ylab=ley,cex.main=0.7)
                     ,key.title=title(main=lacle,cex.main=0.8)
                     ,key.log=pnlog
                     ,plot.axes={axis(side=1)
                                 axis(side=2,at=zat,labels=znice)
                                 if (hovgrid){
                                    abline(h=zat,v=axTicks(1),col=grid.colour,lty="dotted")
                                 }#end if hovgrid
                                }#end plot.axes
                     )

            if (outform[o] %in% c("x11","quartz")){
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
   #      Bar plot by DBH class.                                                           #
   #---------------------------------------------------------------------------------------#
   cat("    + Bar plot by DBH classes...","\n")
   pftuse      = which(apply(X=szpft$nplant,MARGIN=3,FUN=sum,na.rm=TRUE) > 0.)
   pftuse      = pftuse[pftuse != (npft+1)]
   npftuse     = length(pftuse)
   pftname.use = pft$name  [pftuse]
   pftcol.use  = pft$colour[pftuse]
   for (v in sequence(ntspftdbh)){
      #----- Load settings for this variable.----------------------------------------------#
      thisbar     = tspftdbh[[v]]
      vnam        = thisbar$vnam
      description = thisbar$desc
      unit        = thisbar$e.unit
      stacked     = thisbar$stack
      plotit      = thisbar$bar.plot && plot.ycomp 
      plog        = thisbar$plog
      if (plog){
         stacked = FALSE
         xylog   = "y"
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
         thisvnam                  = szpft[[vnam]]
         thisvnam                  = thisvnam [,,pftuse]
         thisvnam                  = thisvnam [,-(ndbh+1),]
         
         thisvnam[is.na(thisvnam)] = 0.
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
         for (y in sequence(nyears)){

            #----- Find which year we are plotting. ---------------------------------------#
            cyear     = sprintf("%4.4i",datum$toyear[y])
            yy        = as.numeric(cyear)
            #------------------------------------------------------------------------------#


            #----- Loop over output formats. ----------------------------------------------#
            for (o in sequence(nout)){
               #------ Open the plot. -----------------------------------------------------#
               fichier = file.path( outdir
                                  , paste0(vnam,"-",cyear,"-",suffix,".",outform[o])
                                  )#end file.path
               if (outform[o] %in% "x11"){
                  X11(width=size$width,height=size$height,pointsize=ptsz)
               }else if (outform[o] %in% "quartz"){
                  quartz(width=size$width,height=size$height,pointsize=ptsz)
               }else if (outform[o] %in% "png"){
                  png(filename=fichier,width=size$width*depth,height=size$height*depth
                     ,pointsize=ptsz,res=depth,bg="transparent")
               }else if (outform[o] %in% "tif"){
                  tiff(filename=fichier,width=size$width*depth,height=size$height*depth
                      ,pointsize=ptsz,res=depth,bg="transparent",compression="lzw")
               }else if (outform[o] %in% "eps"){
                  postscript(file=fichier,width=size$width,height=size$height
                            ,pointsize=ptsz,paper=size$paper)
               }else if (outform[o] %in% "pdf"){
                  pdf(file=fichier,onefile=FALSE
                     ,width=size$width,height=size$height,pointsize=ptsz,paper=size$paper)
               }#end if
               #---------------------------------------------------------------------------#


               #----- Split plotting area into two. ---------------------------------------#
               par(par.user)
               layout(mat=rbind(2,1),heights=c(5,1))
               #---------------------------------------------------------------------------#



               #----- Plot the legend. ----------------------------------------------------#
               par(mar=c(0.1,4.1,0.1,2.1))
               plot.new()
               plot.window(xlim=c(0,1),ylim=c(0,1))
               legend( x      = "topleft"
                     , inset  = 0.0
                     , legend = pftname.use
                     , fill   = pftcol.use
                     , ncol   = min(3,pretty.box(n.selpft)$ncol)
                     , title  = expression(bold("Plant functional type"))
                     , cex    = 1.0
                     , bg     = background
                     )#end legend
               #---------------------------------------------------------------------------#


               #------ Set up the title and axis labels. ----------------------------------#
               letitre = paste0(lieu,"\n",description," - Year : ",cyear)
               lexlab  = "DBH Classes"
               leylab  = desc.unit(desc=description,unit=unit)
               #---------------------------------------------------------------------------#


               #----- Plot all monthly means together. ------------------------------------#
               par(mar=c(4.1,4.1,4.1,2.1))
               plot.new()
               barplot(height=t(thisvnam[y,,]),names.arg=dbhnames[sequence(ndbh)],width=1.0
                      ,main=letitre,xlab=lexlab,ylab=leylab,ylim=ylimit,legend.text=FALSE
                      ,beside=(! stacked),col=pftcol.use,log=xylog
                      ,border=grid.colour,xpd=FALSE,cex.main=cex.main,las=1)
               if (plotgrid & (! stacked)){
                  xgrid=0.5+sequence(ndbh)*(1+npftuse)
                  abline(v=xgrid,col=grid.colour,lty="solid")
               }#end if
               box()
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     Close the device.                                                     #
               #---------------------------------------------------------------------------#
               if (outform[o] %in% c("x11","quartz")){
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

         cat("      + Size and age structure plot: ",description,"...","\n")

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
            #------------------------------------------------------------------------------#




            #----- Retrieve variable list, age, DBH, and PFT for this year. ---------------#
            ageww   = cohort$age   [[ww]]
            if (any(ageww <= 0,na.rm=TRUE)){
               minww = min(ageww,na.rm=TRUE)
               ageww = ageww - minww + 0.01
            }#end if
            dbhww    = cohort$dbh   [[ww]]
            pftww    = cohort$pft   [[ww]]
            varww    = varco        [[ww]]
            popww    = cohort$nplant[[ww]] * cohort$area[[ww]]
            #------------------------------------------------------------------------------#




            #------------------------------------------------------------------------------#
            #     We only plot the SAS figures when the polygon is not an absolute desert. #
            #------------------------------------------------------------------------------#
            if (any (! is.na(varww))){
               #---------------------------------------------------------------------------#
               #      Find the range.  If the user wants the range to be fixed, then use   #
               # the global range, otherwise, simply use the range for this year.          #
               #---------------------------------------------------------------------------#
               if (sasfixlimits){
                  xlimit = pretty.xylim(u=unlist(cohort$age),fracexp=0.0,is.log=TRUE )
                  ylimit = pretty.xylim(u=unlist(cohort$dbh),fracexp=0.0,is.log=FALSE)
                  zlimit = pretty.xylim(u=unlist(varco)     ,fracexp=0.0,is.log=plog )
                  popmin = min(unlist(cohort$nplant * cohort$area), na.rm=TRUE)
                  popmax = max(unlist(cohort$nplant * cohort$area), na.rm=TRUE)
               }else{
                  xlimit = pretty.xylim(u=ageww             ,fracexp=0.0,is.log=TRUE )
                  ylimit = pretty.xylim(u=dbhww             ,fracexp=0.0,is.log=FALSE)
                  zlimit = pretty.xylim(u=varww             ,fracexp=0.0,is.log=plog )
                  popmin = min(popww  ,na.rm=TRUE)
                  popmax = max(popww  ,na.rm=TRUE)
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
               ncohnow = length(varww)
               ageww   = rep(ageww,each=3)
               dbhww   = rep(dbhww,each=3)
               pftww   = rep(pftww,each=3)
               varww   = as.vector( rbind( rep(floor3d,times=ncohnow)
                                         , varco[[ww]]
                                         , rep(NA,times=ncohnow)
                                         )#end rbind
                                  )#end as.vector
               xww     = log(ageww)
               yww     = dbhww
               zww     = if(plog){log(varww)}else{varww}
               pchww   = rep(c(NA,16,NA),times=ncohnow)
               cexww   = rep(cexww,each=3)
               colww   = pft$colour[pftww]

               pftin   = sort(unique(cohort$pft[[ww]]))
               colleg  = pft$colour[pftin]
               pftleg  = pft$name  [pftin]
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #   Plot annotation.                                                        #
               #---------------------------------------------------------------------------#
               letitre = paste(description," - ",lieu
                              ,"\n Time :",mlist[mm],"/",thisyear,sep=" ")
               lexlab  = desc.unit(desc="Gap age",unit=untab$yr)
               leylab  = desc.unit(desc="DBH",unit=untab$cm)
               lezlab  = desc.unit(desc=description,unit=unit)
               #---------------------------------------------------------------------------#


               #----- Loop over output formats. -------------------------------------------#
               for (o in sequence(nout)){
                  #----- Open file. -------------------------------------------------------#
                  fichier = file.path( outdir
                                     , paste0( vnam,"-",thisyear,"-",cmonth,"-",suffix
                                             ,".",outform[o]
                                             )#end paste0
                                     )#end file.path
                  if (outform[o] %in% "x11"){
                     X11(width=size$width,height=size$height,pointsize=ptsz)
                  }else if (outform[o] %in% "quartz"){
                     quartz(width=size$width,height=size$height,pointsize=ptsz)
                  }else if(outform[o] %in% "png"){
                     png(filename=fichier,width=size$width*depth,height=size$height*depth
                        ,pointsize=ptsz,res=depth,bg="transparent")
                  }else if(outform[o] %in% "tif"){
                     tiff(filename=fichier,width=size$width*depth,height=size$height*depth
                         ,pointsize=ptsz,res=depth,bg="transparent",compression="lzw")
                  }else if(outform[o] %in% "eps"){
                     postscript(file=fichier,width=size$width,height=size$height
                               ,pointsize=ptsz,paper=size$paper)
                  }else if(outform[o] %in% "pdf"){
                     pdf(file=fichier,onefile=FALSE,width=size$width,height=size$height
                        ,pointsize=ptsz,paper=size$paper)
                  }#end if
                  #------------------------------------------------------------------------#


                  #----- Split the domain into 2. -----------------------------------------#
                  par(par.user)
                  layout(mat=rbind(2,1),heights=c(5,1))
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
                  if (outform[o] %in% c("x11","quartz")){
                     locator(n=1)
                     dev.off()
                  }else{
                     dev.off()
                  }#end if
                  dummy = clean.tmp()
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
#q("no")
