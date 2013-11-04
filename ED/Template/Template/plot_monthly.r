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
ncolsfc        = 200                    # Target number of colours for filled contour.
fcgrid         = TRUE                   # Include a grid on the filled contour plots?
ncolshov       = 200                    # Target number of colours for Hovmoller diagrams.
hovgrid        = TRUE                   # Include a grid on the Hovmoller plots?
legwhere       = "topleft"              # Where should I place the legend?
inset          = 0.01                   # Inset between legend and edge of plot region.
legbg          = "white"                # Legend background colour.
scalleg        = 0.40                   # Expand y limits by this relative amount to fit
                                        #    the legend
cex.main       = 0.8                    # Scale coefficient for the title
theta          = 315.                   # Azimuth for perspective projection
phi            = 30.                    # Vertical angle for perspective projection
ltheta         = -210.                  # Azimuth angle for light
shade          = 0.125                  # Shade intensity
expz           = 0.5                    # Expansion factor for Z axis
gcol           = c("lightblue","white") # Colours for the fifties style floor
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
      sasmonth = sasmonth.short
   }else{
      sasmonth = sasmonth.long
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




   #----- Make some shorter versions of some variables. -----------------------------------#
   mfac     = datum$month
   dcycmean = datum$dcycmean
   dcycmsqu = datum$dcycmsqu
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Here we find the monthly means for month, then compute the standard deviation.   #
   #---------------------------------------------------------------------------------------#
   cat ("    - Finding the monthly mean...","\n")
   cat ("      * Aggregating the monthly mean...","\n")
   mont12mn             = list()
   mont12mn$gpp         = tapply(X=datum$gpp          ,INDEX=mfac      ,FUN=mean,na.rm=T)
   mont12mn$npp         = tapply(X=datum$npp          ,INDEX=mfac      ,FUN=mean,na.rm=T)
   mont12mn$nep         = tapply(X=datum$nep          ,INDEX=mfac      ,FUN=mean,na.rm=T)
   mont12mn$plresp      = tapply(X=datum$plresp       ,INDEX=mfac      ,FUN=mean,na.rm=T)
   mont12mn$leaf.resp   = tapply(X=datum$leaf.resp    ,INDEX=mfac      ,FUN=mean,na.rm=T)
   mont12mn$root.resp   = tapply(X=datum$root.resp    ,INDEX=mfac      ,FUN=mean,na.rm=T)
   mont12mn$growth.resp = tapply(X=datum$growth.resp  ,INDEX=mfac      ,FUN=mean,na.rm=T)
   mont12mn$hetresp     = tapply(X=datum$hetresp      ,INDEX=mfac      ,FUN=mean,na.rm=T)
   mont12mn$cwdresp     = tapply(X=datum$cwdresp      ,INDEX=mfac      ,FUN=mean,na.rm=T)
   mont12mn$cflxca      = tapply(X=datum$cflxca       ,INDEX=mfac      ,FUN=mean,na.rm=T)
   mont12mn$cflxst      = tapply(X=datum$cflxst       ,INDEX=mfac      ,FUN=mean,na.rm=T)
   mont12mn$nee         = tapply(X=datum$nee          ,INDEX=mfac      ,FUN=mean,na.rm=T)
   mont12mn$reco        = tapply(X=datum$reco         ,INDEX=mfac      ,FUN=mean,na.rm=T)
   mont12mn$hflxca      = tapply(X=datum$hflxca       ,INDEX=mfac      ,FUN=mean,na.rm=T)
   mont12mn$hflxlc      = tapply(X=datum$hflxlc       ,INDEX=mfac      ,FUN=mean,na.rm=T)
   mont12mn$hflxwc      = tapply(X=datum$hflxwc       ,INDEX=mfac      ,FUN=mean,na.rm=T)
   mont12mn$hflxgc      = tapply(X=datum$hflxgc       ,INDEX=mfac      ,FUN=mean,na.rm=T)
   mont12mn$wflxca      = tapply(X=datum$wflxca       ,INDEX=mfac      ,FUN=mean,na.rm=T)
   mont12mn$qwflxca     = tapply(X=datum$qwflxca      ,INDEX=mfac      ,FUN=mean,na.rm=T)
   mont12mn$wflxlc      = tapply(X=datum$wflxlc       ,INDEX=mfac      ,FUN=mean,na.rm=T)
   mont12mn$wflxwc      = tapply(X=datum$wflxwc       ,INDEX=mfac      ,FUN=mean,na.rm=T)
   mont12mn$wflxgc      = tapply(X=datum$wflxgc       ,INDEX=mfac      ,FUN=mean,na.rm=T)
   mont12mn$evap        = tapply(X=datum$evap         ,INDEX=mfac      ,FUN=mean,na.rm=T)
   mont12mn$transp      = tapply(X=datum$transp       ,INDEX=mfac      ,FUN=mean,na.rm=T)
   mont12mn$rain        = tapply(X=datum$rain         ,INDEX=mfac      ,FUN=mean,na.rm=T)
   mont12mn$atm.temp    = tapply(X=datum$atm.temp     ,INDEX=mfac      ,FUN=mean,na.rm=T)
   mont12mn$rshort      = tapply(X=datum$rshort       ,INDEX=mfac      ,FUN=mean,na.rm=T)
   mont12mn$rshortup    = tapply(X=datum$rshortup     ,INDEX=mfac      ,FUN=mean,na.rm=T)
   mont12mn$rlong       = tapply(X=datum$rlong        ,INDEX=mfac      ,FUN=mean,na.rm=T)
   mont12mn$rlongup     = tapply(X=datum$rlongup      ,INDEX=mfac      ,FUN=mean,na.rm=T)
   # mont12mn$par.tot     = tapply(X=datum$par.tot      ,INDEX=mfac      ,FUN=mean,na.rm=T)
   mont12mn$parup       = tapply(X=datum$parup        ,INDEX=mfac      ,FUN=mean,na.rm=T)
   mont12mn$rnet        = tapply(X=datum$rnet         ,INDEX=mfac      ,FUN=mean,na.rm=T)
   mont12mn$albedo      = tapply(X=datum$albedo       ,INDEX=mfac      ,FUN=mean,na.rm=T)
   mont12mn$atm.shv     = tapply(X=datum$atm.shv      ,INDEX=mfac      ,FUN=mean,na.rm=T)
   mont12mn$atm.co2     = tapply(X=datum$atm.co2      ,INDEX=mfac      ,FUN=mean,na.rm=T)
   mont12mn$atm.prss    = tapply(X=datum$atm.prss     ,INDEX=mfac      ,FUN=mean,na.rm=T)
   mont12mn$atm.vels    = tapply(X=datum$atm.vels     ,INDEX=mfac      ,FUN=mean,na.rm=T)
   mont12mn$ustar       = tapply(X=datum$ustar        ,INDEX=mfac      ,FUN=mean,na.rm=T)
   mont12mn$soil.temp   = qapply(X=datum$soil.temp    ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   mont12mn$soil.water  = qapply(X=datum$soil.water   ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   mont12mn$soil.mstpot = qapply(X=datum$soil.mstpot  ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   #----- Find the mean sum of squares. ---------------------------------------------------#
   cat ("      * Aggregating the monthly mean sum of squares...","\n")
   mont12sq           = list()
   mont12sq$gpp       = tapply(X=datum$mmsqu.gpp       ,INDEX=mfac     ,FUN=mean,na.rm=T)
   mont12sq$plresp    = tapply(X=datum$mmsqu.plresp    ,INDEX=mfac     ,FUN=mean,na.rm=T)
   mont12sq$leaf.resp = tapply(X=datum$mmsqu.leaf.resp ,INDEX=mfac     ,FUN=mean,na.rm=T)
   mont12sq$root.resp = tapply(X=datum$mmsqu.root.resp ,INDEX=mfac     ,FUN=mean,na.rm=T)
   mont12sq$hetresp   = tapply(X=datum$mmsqu.hetresp   ,INDEX=mfac     ,FUN=mean,na.rm=T)
   mont12sq$cwdresp   = tapply(X=datum$mmsqu.cwdresp   ,INDEX=mfac     ,FUN=mean,na.rm=T)
   mont12sq$cflxca    = tapply(X=datum$mmsqu.cflxca    ,INDEX=mfac     ,FUN=mean,na.rm=T)
   mont12sq$cflxst    = tapply(X=datum$mmsqu.cflxst    ,INDEX=mfac     ,FUN=mean,na.rm=T)
   mont12sq$hflxca    = tapply(X=datum$mmsqu.hflxca    ,INDEX=mfac     ,FUN=mean,na.rm=T)
   mont12sq$hflxlc    = tapply(X=datum$mmsqu.hflxlc    ,INDEX=mfac     ,FUN=mean,na.rm=T)
   mont12sq$hflxwc    = tapply(X=datum$mmsqu.hflxwc    ,INDEX=mfac     ,FUN=mean,na.rm=T)
   mont12sq$hflxgc    = tapply(X=datum$mmsqu.hflxgc    ,INDEX=mfac     ,FUN=mean,na.rm=T)
   mont12sq$wflxca    = tapply(X=datum$mmsqu.wflxca    ,INDEX=mfac     ,FUN=mean,na.rm=T)
   mont12sq$qwflxca   = tapply(X=datum$mmsqu.qwflxca   ,INDEX=mfac     ,FUN=mean,na.rm=T)
   mont12sq$wflxlc    = tapply(X=datum$mmsqu.wflxlc    ,INDEX=mfac     ,FUN=mean,na.rm=T)
   mont12sq$wflxwc    = tapply(X=datum$mmsqu.wflxwc    ,INDEX=mfac     ,FUN=mean,na.rm=T)
   mont12sq$wflxgc    = tapply(X=datum$mmsqu.wflxgc    ,INDEX=mfac     ,FUN=mean,na.rm=T)
   mont12sq$evap      = tapply(X=datum$mmsqu.evap      ,INDEX=mfac     ,FUN=mean,na.rm=T)
   mont12sq$transp    = tapply(X=datum$mmsqu.transp    ,INDEX=mfac     ,FUN=mean,na.rm=T)
   mont12sq$ustar     = tapply(X=datum$mmsqu.ustar     ,INDEX=mfac     ,FUN=mean,na.rm=T)
   mont12sq$albedo    = tapply(X=datum$mmsqu.albedo    ,INDEX=mfac     ,FUN=mean,na.rm=T)
   mont12sq$rshortup  = tapply(X=datum$mmsqu.rshortup  ,INDEX=mfac     ,FUN=mean,na.rm=T)
   mont12sq$rlongup   = tapply(X=datum$mmsqu.rlongup   ,INDEX=mfac     ,FUN=mean,na.rm=T)
   mont12sq$parup     = tapply(X=datum$mmsqu.parup     ,INDEX=mfac     ,FUN=mean,na.rm=T)
   mont12sq$rnet      = tapply(X=datum$mmsqu.rnet      ,INDEX=mfac     ,FUN=mean,na.rm=T)
   #---------------------------------------------------------------------------------------#
   #   Here we convert the sum of squares into standard deviation. The standard devi-      #
   # ation can be written in two different ways, and we will use the latter because it     #
   # doesn't require previous knowledge of the mean.                                       #
   #              __________________          ____________________________________         #
   #             / SUM_i[X_i - Xm]²          /  / SUM_i[X_i²]        \      1              #
   # sigma = \  /  ----------------   =  \  /  |  -----------  - Xm²  | ---------          #
   #          \/       N - 1              \/    \      N             /   1 - 1/N           #
   #                                                                                       #
   # srnonm1 is the square root of 1 / (1 - 1/N)                                           #
   #     Find the standard deviation.                                                      #
   #---------------------------------------------------------------------------------------#
   cat ("      * Finding the standard deviation...","\n")
   srnorm1 = sqrt(1./(1. - 1. / datum$montable))
   srnorm1[!is.finite(srnorm1)] = 0.
   mont12sd            = list()
   mont12sd$gpp        = sqrt(mont12sq$gpp        - mont12mn$gpp^2        ) * srnorm1
   mont12sd$plresp     = sqrt(mont12sq$plresp     - mont12mn$plresp^2     ) * srnorm1
   mont12sd$leaf.resp  = sqrt(mont12sq$leaf.resp  - mont12mn$leaf.resp^2  ) * srnorm1
   mont12sd$root.resp  = sqrt(mont12sq$root.resp  - mont12mn$root.resp^2  ) * srnorm1
   mont12sd$hetresp    = sqrt(mont12sq$hetresp    - mont12mn$hetresp^2    ) * srnorm1
   mont12sd$cwdresp    = sqrt(mont12sq$cwdresp    - mont12mn$cwdresp^2    ) * srnorm1
   mont12sd$cflxca     = sqrt(mont12sq$cflxca     - mont12mn$cflxca^2     ) * srnorm1
   mont12sd$cflxst     = sqrt(mont12sq$cflxst     - mont12mn$cflxst^2     ) * srnorm1
   mont12sd$hflxca     = sqrt(mont12sq$hflxca     - mont12mn$hflxca^2     ) * srnorm1
   mont12sd$hflxlc     = sqrt(mont12sq$hflxlc     - mont12mn$hflxlc^2     ) * srnorm1
   mont12sd$hflxwc     = sqrt(mont12sq$hflxwc     - mont12mn$hflxwc^2     ) * srnorm1
   mont12sd$hflxgc     = sqrt(mont12sq$hflxgc     - mont12mn$hflxgc^2     ) * srnorm1
   mont12sd$wflxca     = sqrt(mont12sq$wflxca     - mont12mn$wflxca^2     ) * srnorm1
   mont12sd$qwflxca    = sqrt(mont12sq$qwflxca    - mont12mn$qwflxca^2    ) * srnorm1
   mont12sd$wflxlc     = sqrt(mont12sq$wflxlc     - mont12mn$wflxlc^2     ) * srnorm1
   mont12sd$wflxwc     = sqrt(mont12sq$wflxwc     - mont12mn$wflxwc^2     ) * srnorm1
   mont12sd$wflxgc     = sqrt(mont12sq$wflxgc     - mont12mn$wflxgc^2     ) * srnorm1
   mont12sd$evap       = sqrt(mont12sq$evap       - mont12mn$evap^2       ) * srnorm1
   mont12sd$transp     = sqrt(mont12sq$transp     - mont12mn$transp^2     ) * srnorm1
   mont12sd$ustar      = sqrt(mont12sq$ustar      - mont12mn$ustar^2      ) * srnorm1
   mont12sd$albedo     = sqrt(mont12sq$albedo     - mont12mn$albedo^2     ) * srnorm1
   mont12sd$rshortup   = sqrt(mont12sq$rshortup   - mont12mn$rshortup^2   ) * srnorm1
   mont12sd$rlongup    = sqrt(mont12sq$rlongup    - mont12mn$rlongup^2    ) * srnorm1
   mont12sd$parup      = sqrt(mont12sq$parup      - mont12mn$parup^2      ) * srnorm1
   mont12sd$rnet       = sqrt(mont12sq$rnet       - mont12mn$rnet^2       ) * srnorm1
   #---------------------------------------------------------------------------------------#
   #     Set standard deviations that became NaN to 0.  This usually happens when we run   #
   # the post-processing script when the simulation hasn't run for more than 2 years.  We  #
   # can't find the standard deviation because the number of degrees of freedom becomes 0. #
   #---------------------------------------------------------------------------------------#
   mont12sd$gpp        [!is.finite(mont12sd$gpp        )] = 0.
   mont12sd$plresp     [!is.finite(mont12sd$plresp     )] = 0.
   mont12sd$leaf.resp  [!is.finite(mont12sd$leaf.resp  )] = 0.
   mont12sd$root.resp  [!is.finite(mont12sd$root.resp  )] = 0.
   mont12sd$hetresp    [!is.finite(mont12sd$hetresp    )] = 0.
   mont12sd$cwdresp    [!is.finite(mont12sd$cwdresp    )] = 0.
   mont12sd$cflxca     [!is.finite(mont12sd$cflxca     )] = 0.
   mont12sd$cflxst     [!is.finite(mont12sd$cflxst     )] = 0.
   mont12sd$hflxca     [!is.finite(mont12sd$hflxca     )] = 0.
   mont12sd$hflxlc     [!is.finite(mont12sd$hflxlc     )] = 0.
   mont12sd$hflxlc     [!is.finite(mont12sd$hflxwc     )] = 0.
   mont12sd$hflxgc     [!is.finite(mont12sd$hflxgc     )] = 0.
   mont12sd$wflxca     [!is.finite(mont12sd$wflxca     )] = 0.
   mont12sd$qwflxca    [!is.finite(mont12sd$qwflxca    )] = 0.
   mont12sd$wflxlc     [!is.finite(mont12sd$wflxlc     )] = 0.
   mont12sd$wflxwc     [!is.finite(mont12sd$wflxwc     )] = 0.
   mont12sd$wflxgc     [!is.finite(mont12sd$wflxgc     )] = 0.
   mont12sd$evap       [!is.finite(mont12sd$evap       )] = 0.
   mont12sd$transp     [!is.finite(mont12sd$transp     )] = 0.
   mont12sd$ustar      [!is.finite(mont12sd$ustar      )] = 0.
   mont12sd$albedo     [!is.finite(mont12sd$albedo     )] = 0.
   mont12sd$rshortup   [!is.finite(mont12sd$rshortup   )] = 0.
   mont12sd$rlongup    [!is.finite(mont12sd$rlongup    )] = 0.
   mont12sd$parup      [!is.finite(mont12sd$parup      )] = 0.
   mont12sd$rnet       [!is.finite(mont12sd$rnet       )] = 0.
   #---------------------------------------------------------------------------------------#
   #     Estimate the standard deviation of NPP, NEP, NEE, and REco.                       #
   #---------------------------------------------------------------------------------------#
   mont12sd$npp  = sqrt(mont12sd$gpp^2    + mont12sd$plresp^2                       )
   mont12sd$nep  = sqrt(mont12sd$gpp^2    + mont12sd$plresp^2  + mont12sd$hetresp^2 )
   mont12sd$nee  = sqrt(mont12sd$cflxca^2 + mont12sd$cflxst^2                       )
   mont12sd$reco = sqrt(mont12sd$plresp^2 + mont12sd$hetresp^2                      )
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Here we find the Mean diurnal cycle for each month, then compute the standard    #
   # deviation.                                                                            #
   #---------------------------------------------------------------------------------------#
   cat ("    - Aggregating the monthly mean of the diurnal cycle...","\n")
   dcyc12mn             =list()
   dcyc12mn$gpp         =qapply(X=dcycmean$gpp         ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$npp         =qapply(X=dcycmean$npp         ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$plresp      =qapply(X=dcycmean$plresp      ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$leaf.resp   =qapply(X=dcycmean$leaf.resp   ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$root.resp   =qapply(X=dcycmean$root.resp   ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$hetresp     =qapply(X=dcycmean$hetresp     ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$cwdresp     =qapply(X=dcycmean$cwdresp     ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$nep         =qapply(X=dcycmean$nep         ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$nee         =qapply(X=dcycmean$nee         ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$reco        =qapply(X=dcycmean$reco        ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$cflxca      =qapply(X=dcycmean$cflxca      ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$cflxst      =qapply(X=dcycmean$cflxst      ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$hflxca      =qapply(X=dcycmean$hflxca      ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$hflxlc      =qapply(X=dcycmean$hflxlc      ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$hflxwc      =qapply(X=dcycmean$hflxwc      ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$hflxgc      =qapply(X=dcycmean$hflxgc      ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$wflxca      =qapply(X=dcycmean$wflxca      ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$qwflxca     =qapply(X=dcycmean$qwflxca     ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$wflxlc      =qapply(X=dcycmean$wflxlc      ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$wflxwc      =qapply(X=dcycmean$wflxwc      ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$wflxgc      =qapply(X=dcycmean$wflxgc      ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$evap        =qapply(X=dcycmean$evap        ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$transp      =qapply(X=dcycmean$transp      ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$atm.temp    =qapply(X=dcycmean$atm.temp    ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$can.temp    =qapply(X=dcycmean$can.temp    ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$leaf.temp   =qapply(X=dcycmean$leaf.temp   ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$wood.temp   =qapply(X=dcycmean$wood.temp   ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$gnd.temp    =qapply(X=dcycmean$gnd.temp    ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$atm.shv     =qapply(X=dcycmean$atm.shv     ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$can.shv     =qapply(X=dcycmean$can.shv     ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$gnd.shv     =qapply(X=dcycmean$gnd.shv     ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$atm.vpd     =qapply(X=dcycmean$atm.vpd     ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$can.vpd     =qapply(X=dcycmean$can.vpd     ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$leaf.vpd    =qapply(X=dcycmean$leaf.vpd    ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$atm.co2     =qapply(X=dcycmean$atm.co2     ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$can.co2     =qapply(X=dcycmean$can.co2     ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$atm.prss    =qapply(X=dcycmean$atm.prss    ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$can.prss    =qapply(X=dcycmean$can.prss    ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$atm.vels    =qapply(X=dcycmean$atm.vels    ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$ustar       =qapply(X=dcycmean$ustar       ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$fs.open     =qapply(X=dcycmean$fs.open     ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$rain        =qapply(X=dcycmean$rain        ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$rshort      =qapply(X=dcycmean$rshort      ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$rshort.beam =qapply(X=dcycmean$rshort.beam ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$rshort.diff =qapply(X=dcycmean$rshort.diff ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$rshort.gnd  =qapply(X=dcycmean$rshort.gnd  ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$rshortup    =qapply(X=dcycmean$rshortup    ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$rlong       =qapply(X=dcycmean$rlong       ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$rlong.gnd   =qapply(X=dcycmean$rlong.gnd   ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$rlongup     =qapply(X=dcycmean$rlongup     ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   # dcyc12mn$par.tot     =qapply(X=dcycmean$par.tot     ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   # dcyc12mn$par.beam    =qapply(X=dcycmean$par.beam    ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   # dcyc12mn$par.diff    =qapply(X=dcycmean$par.diff    ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   # dcyc12mn$par.gnd     =qapply(X=dcycmean$par.gnd     ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$parup       =qapply(X=dcycmean$parup       ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$rnet        =qapply(X=dcycmean$rnet        ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$albedo      =qapply(X=dcycmean$albedo      ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$albedo.beam =qapply(X=dcycmean$albedo.beam ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$albedo.diff =qapply(X=dcycmean$albedo.diff ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$rlong.albedo=qapply(X=dcycmean$rlong.albedo,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)

   #----- Find the mean sum of squares. ---------------------------------------------------#
   cat ("    - Aggregating the monthly mean sum of squares...","\n")
   dcyc12sq            = list()
   dcyc12sq$gpp        = qapply(X=dcycmsqu$gpp       ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12sq$plresp     = qapply(X=dcycmsqu$plresp    ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12sq$leaf.resp  = qapply(X=dcycmsqu$leaf.resp ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12sq$root.resp  = qapply(X=dcycmsqu$root.resp ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12sq$hetresp    = qapply(X=dcycmsqu$hetresp   ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12sq$cwdresp    = qapply(X=dcycmsqu$cwdresp   ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12sq$nep        = qapply(X=dcycmsqu$nep       ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12sq$cflxca     = qapply(X=dcycmsqu$cflxca    ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12sq$cflxst     = qapply(X=dcycmsqu$cflxst    ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12sq$hflxca     = qapply(X=dcycmsqu$hflxca    ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12sq$hflxlc     = qapply(X=dcycmsqu$hflxlc    ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12sq$hflxwc     = qapply(X=dcycmsqu$hflxwc    ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12sq$hflxgc     = qapply(X=dcycmsqu$hflxgc    ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12sq$wflxca     = qapply(X=dcycmsqu$wflxca    ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12sq$qwflxca    = qapply(X=dcycmsqu$qwflxca   ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12sq$wflxlc     = qapply(X=dcycmsqu$wflxlc    ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12sq$wflxwc     = qapply(X=dcycmsqu$wflxwc    ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12sq$wflxgc     = qapply(X=dcycmsqu$wflxgc    ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12sq$transp     = qapply(X=dcycmsqu$transp    ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12sq$ustar      = qapply(X=dcycmsqu$ustar     ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12sq$albedo     = qapply(X=dcycmsqu$albedo    ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12sq$rshortup   = qapply(X=dcycmsqu$rshortup  ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12sq$rlongup    = qapply(X=dcycmsqu$rlongup   ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12sq$parup      = qapply(X=dcycmsqu$parup     ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12sq$rnet       = qapply(X=dcycmsqu$rnet      ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)

   #---------------------------------------------------------------------------------------#
   #   Here we convert the sum of squares into standard deviation. The standard devi-      #
   # ation can be written in two different ways, and we will use the latter because it     #
   # doesn't require previous knowledge of the mean.                                       #
   #              __________________          ____________________________________         #
   #             / SUM_i[X_i - Xm]²          /  / SUM_i[X_i²]        \      1              #
   # sigma = \  /  ----------------   =  \  /  |  -----------  - Xm²  | ---------          #
   #          \/       N - 1              \/    \      N             /   1 - 1/N           #
   #                                                                                       #
   # srnonm1 is the square root of 1 / (1 - 1/N)                                           #
   #     Find the standard deviation.                                                      #
   #---------------------------------------------------------------------------------------#
   cat ("    - Finding the standard deviation...","\n")
   srnorm1 = sqrt(1./(1. - 1. / datum$moncnt))
   srnorm1[!is.finite(srnorm1)] = 0.
   dcyc12sd            = list()
   dcyc12sd$gpp        = sqrt(dcyc12sq$gpp       - dcyc12mn$gpp^2           )*srnorm1
   dcyc12sd$plresp     = sqrt(dcyc12sq$plresp    - dcyc12mn$plresp^2        )*srnorm1
   dcyc12sd$leaf.resp  = sqrt(dcyc12sq$leaf.resp - dcyc12mn$leaf.resp^2     )*srnorm1
   dcyc12sd$root.resp  = sqrt(dcyc12sq$root.resp - dcyc12mn$root.resp^2     )*srnorm1
   dcyc12sd$hetresp    = sqrt(dcyc12sq$hetresp   - dcyc12mn$hetresp^2       )*srnorm1
   dcyc12sd$cwdresp    = sqrt(dcyc12sq$cwdresp   - dcyc12mn$cwdresp^2       )*srnorm1
   dcyc12sd$nep        = sqrt(dcyc12sq$nep       - dcyc12mn$nep^2           )*srnorm1
   dcyc12sd$cflxca     = sqrt(dcyc12sq$cflxca    - dcyc12mn$cflxca^2        )*srnorm1
   dcyc12sd$cflxst     = sqrt(dcyc12sq$cflxst    - dcyc12mn$cflxst^2        )*srnorm1
   dcyc12sd$hflxca     = sqrt(dcyc12sq$hflxca    - dcyc12mn$hflxca^2        )*srnorm1
   dcyc12sd$hflxlc     = sqrt(dcyc12sq$hflxlc    - dcyc12mn$hflxlc^2        )*srnorm1
   dcyc12sd$hflxwc     = sqrt(dcyc12sq$hflxwc    - dcyc12mn$hflxwc^2        )*srnorm1
   dcyc12sd$hflxgc     = sqrt(dcyc12sq$hflxgc    - dcyc12mn$hflxgc^2        )*srnorm1
   dcyc12sd$wflxca     = sqrt(dcyc12sq$wflxca    - dcyc12mn$wflxca^2        )*srnorm1
   dcyc12sd$qwflxca    = sqrt(dcyc12sq$qwflxca   - dcyc12mn$qwflxca^2       )*srnorm1
   dcyc12sd$wflxlc     = sqrt(dcyc12sq$wflxlc    - dcyc12mn$wflxlc^2        )*srnorm1
   dcyc12sd$wflxwc     = sqrt(dcyc12sq$wflxwc    - dcyc12mn$wflxwc^2        )*srnorm1
   dcyc12sd$wflxgc     = sqrt(dcyc12sq$wflxgc    - dcyc12mn$wflxgc^2        )*srnorm1
   dcyc12sd$transp     = sqrt(dcyc12sq$transp    - dcyc12mn$transp^2        )*srnorm1
   dcyc12sd$ustar      = sqrt(dcyc12sq$ustar     - dcyc12mn$ustar^2         )*srnorm1
   dcyc12sd$albedo     = sqrt(dcyc12sq$albedo    - dcyc12mn$albedo^2        )*srnorm1
   dcyc12sd$rshortup   = sqrt(dcyc12sq$rshortup  - dcyc12mn$rshortup^2      )*srnorm1
   dcyc12sd$rlongup    = sqrt(dcyc12sq$rlongup   - dcyc12mn$rlongup^2       )*srnorm1
   dcyc12sd$parup      = sqrt(dcyc12sq$parup     - dcyc12mn$parup^2         )*srnorm1
   dcyc12sd$rnet       = sqrt(dcyc12sq$rnet      - dcyc12mn$rnet^2          )*srnorm1
   #---------------------------------------------------------------------------------------#
   #     Set standard deviations that became NaN to 0.  This usually happens when we run   #
   # the post-processing script when the simulation hasn't run for more than 2 years.  We  #
   # can't find the standard deviation because the number of degrees of freedom becomes 0. #
   #---------------------------------------------------------------------------------------#
   dcyc12sd$gpp        [! is.finite(dcyc12sd$gpp       )] = 0.
   dcyc12sd$plresp     [! is.finite(dcyc12sd$plresp    )] = 0.
   dcyc12sd$leaf.resp  [! is.finite(dcyc12sd$leaf.resp )] = 0.
   dcyc12sd$root.resp  [! is.finite(dcyc12sd$root.resp )] = 0.
   dcyc12sd$hetresp    [! is.finite(dcyc12sd$hetresp   )] = 0.
   dcyc12sd$cwdresp    [! is.finite(dcyc12sd$cwdresp   )] = 0.
   dcyc12sd$nep        [! is.finite(dcyc12sd$nep       )] = 0.
   dcyc12sd$cflxca     [! is.finite(dcyc12sd$cflxca    )] = 0.
   dcyc12sd$cflxst     [! is.finite(dcyc12sd$cflxst    )] = 0.
   dcyc12sd$hflxca     [! is.finite(dcyc12sd$hflxca    )] = 0.
   dcyc12sd$hflxlc     [! is.finite(dcyc12sd$hflxlc    )] = 0.
   dcyc12sd$hflxwc     [! is.finite(dcyc12sd$hflxwc    )] = 0.
   dcyc12sd$hflxgc     [! is.finite(dcyc12sd$hflxgc    )] = 0.
   dcyc12sd$wflxca     [! is.finite(dcyc12sd$wflxca    )] = 0.
   dcyc12sd$qwflxca    [! is.finite(dcyc12sd$qwflxca   )] = 0.
   dcyc12sd$wflxlc     [! is.finite(dcyc12sd$wflxlc    )] = 0.
   dcyc12sd$wflxwc     [! is.finite(dcyc12sd$wflxwc    )] = 0.
   dcyc12sd$wflxgc     [! is.finite(dcyc12sd$wflxgc    )] = 0.
   dcyc12sd$transp     [! is.finite(dcyc12sd$transp    )] = 0.
   dcyc12sd$ustar      [! is.finite(dcyc12sd$ustar     )] = 0.
   dcyc12sd$albedo     [! is.finite(dcyc12sd$albedo    )] = 0.
   dcyc12sd$rshortup   [! is.finite(dcyc12sd$rshortup  )] = 0.
   dcyc12sd$rlongup    [! is.finite(dcyc12sd$rlongup   )] = 0.
   dcyc12sd$parup      [! is.finite(dcyc12sd$parup     )] = 0.
   dcyc12sd$rnet       [! is.finite(dcyc12sd$rnet      )] = 0.
   #---------------------------------------------------------------------------------------#
   #      Estimate NPP and NEE standard deviation.                                         #
   #---------------------------------------------------------------------------------------#
   dcyc12sd$npp  = sqrt(dcyc12sd$gpp^2    + dcyc12sd$plresp^2 )
   dcyc12sd$nee  = sqrt(dcyc12sd$cflxca^2 + dcyc12sd$cflxst^2 )
   dcyc12sd$reco = sqrt(dcyc12sd$plresp^2 + dcyc12sd$hetresp^2)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Remove all elements of the DBH/PFT class that do not have a single valid cohort   #
   # at any given time.                                                                    #
   #---------------------------------------------------------------------------------------#
   empty = datum$nplantpftdbh == 0
   datum$agbpftdbh          [empty] = NA
   datum$basareapftdbh      [empty] = NA
   datum$laipftdbh          [empty] = NA
   datum$waipftdbh          [empty] = NA
   datum$taipftdbh          [empty] = NA
   datum$gpppftdbh          [empty] = NA
   datum$npppftdbh          [empty] = NA
   datum$mcopftdbh          [empty] = NA
   datum$cbapftdbh          [empty] = NA
   datum$cbalightpftdbh     [empty] = NA
   datum$cbamoistpftdbh     [empty] = NA
   datum$cbal12lightpftdbh  [empty] = NA
   datum$cbal12moistpftdbh  [empty] = NA
   datum$cbarelpftdbh       [empty] = NA
   datum$ldroppftdbh        [empty] = NA
   datum$fsopftdbh          [empty] = NA
   datum$demandpftdbh       [empty] = NA
   datum$supplypftdbh       [empty] = NA
   datum$mortpftdbh         [empty] = NA
   datum$agemortpftdbh      [empty] = NA
   datum$ncbmortpftdbh      [empty] = NA
   datum$tfallmortpftdbh    [empty] = NA
   datum$coldmortpftdbh     [empty] = NA
   datum$distmortpftdbh     [empty] = NA
   datum$growthpftdbh       [empty] = NA
   datum$plresppftdbh       [empty] = NA
   datum$bstorepftdbh       [empty] = NA
   datum$hflxlcpftdbh       [empty] = NA
   datum$wflxlcpftdbh       [empty] = NA
   datum$transppftdbh       [empty] = NA
   datum$i.gpppftdbh        [empty] = NA
   datum$i.npppftdbh        [empty] = NA
   datum$i.plresppftdbh     [empty] = NA
   datum$i.mcopftdbh        [empty] = NA
   datum$i.cbapftdbh        [empty] = NA
   datum$i.cbalightpftdbh   [empty] = NA
   datum$i.cbamoistpftdbh   [empty] = NA
   datum$i.cbal12lightpftdbh[empty] = NA
   datum$i.cbam12lightpftdbh[empty] = NA
   datum$i.transppftdbh     [empty] = NA
   datum$i.wflxlcpftdbh     [empty] = NA
   datum$i.hflxlcpftdbh     [empty] = NA
   datum$nplantpftdbh       [empty] = NA
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Replace the mortality and recruitment exponential rates by the "interests" rates. #
   #---------------------------------------------------------------------------------------#
   datum$mortpftdbh      = 100. * (1.0 - exp(- datum$mortpftdbh     ))
   datum$agemortpftdbh   = 100. * (1.0 - exp(- datum$agemortpftdbh  ))
   datum$ncbmortpftdbh   = 100. * (1.0 - exp(- datum$ncbmortpftdbh  ))
   datum$tfallmortpftdbh = 100. * (1.0 - exp(- datum$tfallmortpftdbh))
   datum$coldmortpftdbh  = 100. * (1.0 - exp(- datum$coldmortpftdbh ))
   datum$distmortpftdbh  = 100. * (1.0 - exp(- datum$distmortpftdbh ))
   datum$mortpft         = 100. * (1.0 - exp(- datum$mortpft        ))
   datum$agemortpft      = 100. * (1.0 - exp(- datum$agemortpft     ))
   datum$ncbmortpft      = 100. * (1.0 - exp(- datum$ncbmortpft     ))
   datum$tfallmortpft    = 100. * (1.0 - exp(- datum$tfallmortpft   ))
   datum$coldmortpft     = 100. * (1.0 - exp(- datum$coldmortpft    ))
   datum$distmortpft     = 100. * (1.0 - exp(- datum$distmortpft    ))
   datum$recrpft         = 100. * (exp(  datum$recrpft        ) - 1.0)
   #---------------------------------------------------------------------------------------#


   #----- Find which PFTs, land uses and transitions we need to consider ------------------#
   pftave  = colMeans(datum$agbpft,na.rm=TRUE)
   luave   = colMeans(datum$agblu ,na.rm=TRUE)
   distave = matrix(NA,nrow=3,ncol=3)
   for (jlu in 1:nlu){
      for (ilu in 1:nlu){
          distave[ilu,jlu] = mean(datum$dist[,ilu,jlu])
      }#end for
   }#end for
   selpft  = pftave  > 0.
   sellu   = luave   > 0.
   seldist = distave > 0.
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
   dcycplot = list()
   dcycplot$levels = c(0,4,8,12,16,20,24)
   dcycplot$n      = 7
   dcycplot$scale  = "hours"
   dcycplot$padj   = rep(0,times=dcycplot$n)
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
   montplot  = list()
   montplot$levels = montmont
   montplot$labels = capwords(mon2mmm(montmont))
   montplot$n      = 12
   montplot$scale  = "months"
   montplot$padj   = rep(0,times=dcycplot$n)
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
   for (v in 1:ntspft){
      thistspft   = tspft[[v]]
      vnam        = thistspft$vnam
      description = thistspft$desc
      unit        = thistspft$unit
      plog        = thistspft$plog
      plotit      = thistspft$plt

      #----- Check whether the user wants to have this variable plotted. ------------------#
      if (plotit && any(selpft)){

         #---------------------------------------------------------------------------------#
         #    Check whether the time series directory exists.  If not, create it.          #
         #---------------------------------------------------------------------------------#
         outdir = paste(outpref,"tspft",sep="/")
         if (! file.exists(outdir)) dir.create(outdir)
         cat("      + ",description," time series for all PFTs...","\n")

         #----- Load variable -------------------------------------------------------------#
         thisvar = datum[[vnam]]
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
            plot(x=datum$tomonth,y=thisvar[,1],type="n",main=letitre,ylim=ylimit
                ,xlab="Time",xaxt="n",ylab=unit,cex.main=0.7,log=xylog)
            axis(side=1,at=whenplot8$levels,labels=whenplot8$labels,padj=whenplot8$padj)

            if (drought.mark){
               for (n in 1:ndrought){
                  rect(xleft  = drought[[n]][1],ybottom = ydrought[1]
                      ,xright = drought[[n]][2],ytop    = ydrought[2]
                      ,col    = "gray84",border=NA)
               }#end for
            }#end if
            
            if (plotgrid){ 
               abline(v=whenplot8$levels,h=axTicks(side=2),col="gray52",lty="solid")
            }#end if
            for (n in 1:(npft+1)){
               if (selpft[n]){
                  lines(datum$tomonth,thisvar[,n],type="l",col=pft$colour[n],lwd=lwidth)
               }#end if
            }#end for
            legend(x=legwhere,inset=inset,bg=legbg,legend=legs,col=cols,lwd=lwidth)

            if (outform[o] == "x11"){
               locator(n=1)
               dev.off()
            }else{
               dev.off()
            }#end if
         } #end for outform
      }#end if (tseragbpft)
   } #end for tseries
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #      Time series by DBH, by PFT.                                                      #
   #---------------------------------------------------------------------------------------#
   #----- Find the PFTs to plot. ----------------------------------------------------------#
   pftuse  = which(apply(X=datum$agbpftdbh,MARGIN=3,FUN=sum,na.rm=TRUE) > 0.)
   pftuse  = pftuse[pftuse != (npft+1)]
   for (v in 1:ntspftdbh){
      thistspftdbh   = tspftdbh[[v]]
      vnam           = thistspftdbh$vnam
      description    = thistspftdbh$desc
      unit           = thistspftdbh$unit
      plog           = thistspftdbh$plog
      plotit         = thistspftdbh$plt
      
      #----- Load variable ----------------------------------------------------------------#
      thisvar = datum[[vnam]]
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
               plot(x=datum$tomonth,y=thisvar[,1,p],type="n",main=letitre,ylim=ylimit
                   ,xlab="Time",xaxt="n",ylab=unit,cex.main=0.7,log=xylog)
               axis(side=1,at=whenplot8$levels,labels=whenplot8$labels,padj=whenplot8$padj)
               if (drought.mark){
                  for (n in 1:ndrought){
                     rect(xleft  = drought[[n]][1],ybottom = ydrought[1]
                         ,xright = drought[[n]][2],ytop    = ydrought[2]
                         ,col    = "gray84",border=NA)
                  }#end for
               }#end if
               if (plotgrid){ 
                  abline(v=whenplot8$levels,h=axTicks(side=2),col="gray52",lty="solid")
               }#end if
               for (d in seq(from=1,to=ndbh+1,by=1)){
                  lines(datum$tomonth,thisvar[,d,p],type="l",col=dbhcols[d],lwd=lwidth)
               }#end for
               legend(x=legwhere,inset=inset,bg=legbg,legend=dbhnames,col=dbhcols
                     ,ncol=2,title="DBH class",lwd=lwidth,cex=0.8)

               if (outform[o] == "x11"){
                  locator(n=1)
                  dev.off()
               }else{
                  dev.off()
               }#end if
            } #end for outform
         }#end for (p in pftuse)
      }#end if (tseragbpft)
   } #end for tseries
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #   Plot the comparison between observations and model.                                 #
   #---------------------------------------------------------------------------------------#
   cat("    + Comparisons of time series (model vs. observations)...","\n")
   for (cc in 1:ncompemean){

      #----- Retrieve variable information from the list. ---------------------------------#
      compnow      = compemean[[cc]]
      vname        = compnow$vnam  
      description  = compnow$desc  
      unit         = compnow$unit  
      lcolours     = compnow$colour
      llwd         = compnow$lwd
      llwd         = compnow$lwd
      ltype        = compnow$type
      plog         = compnow$plog
      legpos       = compnow$legpos
      plotit       = compnow$plt

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
         thismean  = datum[[vname]][sel]
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
            plot(x=thiswhen,y=thismean,type="n",main=letitre,xlab="Time"
                ,ylim=ylimit,ylab=paste("[",unit,"]",sep=""),log=plog,xaxt="n"
                ,cex.main=cex.main)
            axis(side=1,at=whenplote$levels,labels=whenplote$labels,padj=whenplote$padj)
            if (plotgrid){
               abline(v=whenplote$levels,h=axTicks(side=2),col="gray52",lty="solid")
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
         } #end for outform
      }#end if plotit
   }#end for ncompare
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #   Plot the comparison between observations and model.                                 #
   #---------------------------------------------------------------------------------------#
   cat("    + Comparisons of monthly means (model vs. observations)...","\n")
   for (cc in 1:ncompmmean){

      #----- Retrieve variable information from the list. ---------------------------------#
      compnow      = compmmean[[cc]]
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
      plotit       = compnow$plt

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
         thismean  = mont12mn[[vname]]
         thissdev = mont12sd[[vname]]
         if (length(mont12sd[[vname]]) == 0){
            thissdev = 0. * thismean
         }else{
            thissdev = mont12sd[[vname]]
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
            plot(x=montmont,y=thismean,type="n",main=letitre,xlab="Time"
                ,ylim=ylimit,ylab=paste("[",unit,"]",sep=""),log=plog,xaxt="n"
                ,cex.main=cex.main)
            axis(side=1,at=montplot$levels,labels=montplot$labels,padj=montplot$padj)
            if (plotgrid){ 
               abline(v=montplot$levels,h=axTicks(side=2),col="gray52",lty="solid")
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
                     ,bg="white",title="Shaded areas = 1 SD",cex=1.0,pch=16)
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
         } #end for outform
      }#end if plotit
   }#end for ncompare
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #   Plot the comparison between observations and model.                                 #
   #---------------------------------------------------------------------------------------#
   cat("      * Comparisons of mean diurnal cycle (model vs. observations)...","\n")
   for (cc in 1:ncompdcyc){

      #----- Retrieve variable information from the list. ---------------------------------#
      compnow      = compdcyc[[cc]]
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
      plotit       = compnow$plt

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
         thismean  = dcyc12mn[[vname]]
         thismean  = cbind(thismean[,ndcycle],thismean)
         if (length(dcyc12sd[[vname]]) == 0){
            thissdev = 0. * thismean
         }else{
            thissdev = dcyc12sd[[vname]]
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
               plot(x=thisday,y=thismean[pmon,],type="n",main=letitre,xlab="Time"
                   ,ylim=ylimit,ylab=paste("[",unit,"]",sep=""),log=plog,xaxt="n"
                   ,cex.main=cex.main)
               axis(side=1,at=dcycplot$levels,labels=dcycplot$labels,padj=dcycplot$padj)
               if (plotgrid){ 
                  abline(v=dcycplot$levels,h=axTicks(side=2),col="gray52",lty="solid")
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
                        ,bg="white",title="Shaded areas = 1 SD",cex=1.0,pch=16)
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
         thisvar = datum[[vnam]]
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
            ylimit = pretty.xylim(thisvar[,selpft],fracexp=scalleg,is.log=plog)
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
            plot(datum$tomonth,thisvar[,1],type="n",main=letitre,ylim=ylimit
                ,xlab="Time",ylab=unit,xaxt="n",cex.main=0.7)
            axis(side=1,at=whenplot8$levels,labels=whenplot8$labels,padj=whenplot8$padj)

            if (drought.mark){
               for (n in 1:ndrought){
                  rect(xleft  = drought[[n]][1],ybottom = ydrought[1]
                      ,xright = drought[[n]][2],ytop    = ydrought[2]
                      ,col    = "gray84",border=NA)
               }#end for
            }#end if
            if (plotgrid){ 
               abline(v=whenplot8$levels,h=axTicks(side=2),col="gray52",lty="solid")
            }#end if
            for (n in 1:(nlu+1)){
               if (sellu[n]){
                  lines(datum$tomonth,thisvar[,n],type="l",col=lucols[n],lwd=lwidth)
               }#end if
            }#end for
            legend(x=legwhere,inset=inset,bg=legbg,legend=legs,col=cols,lwd=lwidth)

            if (outform[o] == "x11"){
               locator(n=1)
               dev.off()
            }else{
               dev.off()
            }#end if
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
                  ylimit = c(ylimit,datum$dist[,ilu,jlu])
               }#end if
            }#end for
         }#end for
         ylimit   = pretty.xylim(u=ylimit,fracexp=scalleg,is.log=FALSE)
         ydrought = c(ylimit[1] - 0.5 * diff(ylimit), ylimit[2] + 0.5 * diff(ylimit))
         #---------------------------------------------------------------------------------#

         letitre = paste("Disturbance rates",lieu,sep=" - ")
         cols    = NULL
         legs    = NULL
         plot(datum$tomonth,datum$dist[,1,1],type="n",main=letitre,ylim=ylimit
             ,xlab="Time",ylab="[1/yr]",xaxt="n",cex.main=0.7)
            axis(side=1,at=whenplot8$levels,labels=whenplot8$labels,padj=whenplot8$padj)
            if (drought.mark){
               for (n in 1:ndrought){
                  rect(xleft  = drought[[n]][1],ybottom = ydrought[1]
                      ,xright = drought[[n]][2],ytop    = ydrought[2]
                      ,col    = "gray84",border=NA)
               }#end for
            }#end if
            if (plotgrid){ 
               abline(v=whenplot8$levels,h=axTicks(side=2),col="gray52",lty="solid")
            }#end if
         n = 0
         for (jlu in 1:nlu){
            for (ilu in 1:nlu){
               n = n + 1
               if (seldist[ilu,jlu]){
                  cols = c(cols,distcols[n])
                  legs = c(legs,distnames[n])
                  lines(datum$tomonth,datum$dist[,ilu,jlu],type="l"
                       ,col=distcols[n],lwd=lwidth)
               }#end if
            }#end for
         }#end for
         legend(x=legwhere,inset=inset,bg=legbg,legend=legs,col=cols,lwd=lwidth)

         if (outform[o] == "x11"){
            locator(n=1)
            dev.off()
         }else{
            dev.off()
         }#end if
      } #end for outform
   }#end if
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #   Plot the time series diagrams showing months and years.                             #
   #---------------------------------------------------------------------------------------#
   cat("      * Plot some time series figures...","\n")
   for (hh in 1:ntser){

      #----- Retrieve variable information from the list. ---------------------------------#
      tsernow      = tser[[hh]]
      vnames       = tsernow$vnam  
      description  = tsernow$desc  
      lcolours     = tsernow$colour
      llwd         = tsernow$lwd
      ltype        = tsernow$type
      plog         = tsernow$plog
      prefix       = tsernow$prefix
      theme        = tsernow$theme 
      unit         = tsernow$unit  
      legpos       = tsernow$legpos
      plotit       = tsernow$plt   
   
      if (plotit){

         #---------------------------------------------------------------------------------#
         #    Check whether the time series directory exists.  If not, create it.          #
         #---------------------------------------------------------------------------------#
         outdir = paste(outpref,"tseries",sep="/")
         if (! file.exists(outdir)) dir.create(outdir)
         cat("      + ",theme," time series for several variables...","\n")


         #----- Define the number of layers. ----------------------------------------------#
         nlayers   = length(vnames)
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Find the limit, make some room for the legend, and in case the field is a   #
         # constant, nudge the limits so the plot command will not complain.               #
         #---------------------------------------------------------------------------------#
         ylimit    = NULL
         for (l in 1:nlayers) ylimit  = c(ylimit,datum[[vnames[l]]])
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
         cat("        > ",theme," time series ...","\n")

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
            thisvar = datum[[vnames[1]]]

            letitre = paste(theme," - ",lieu," \n"," Time series: ",theme,sep="")

            plot(x=datum$tomonth,y=thisvar,type="n",main=letitre,xlab="Time"
                ,ylim=ylimit,ylab=paste("[",unit,"]",sep=""),log=xylog,xaxt="n"
                ,cex.main=cex.main)
            axis(side=1,at=whenplot8$levels,labels=whenplot8$labels,padj=whenplot8$padj)
            if (drought.mark){
               for (n in 1:ndrought){
                  rect(xleft  = drought[[n]][1],ybottom = ydrought[1]
                      ,xright = drought[[n]][2],ytop    = ydrought[2]
                      ,col    = "gray84",border=NA)
               }#end for
            }#end if
            if (plotgrid){ 
               abline(v=whenplot8$levels,h=axTicks(side=2),col="gray52",lty="solid")
            }#end if
            for (l in 1:nlayers){
               thisvar = datum[[vnames[l]]]
               points(x=datum$tomonth,y=thisvar,col=lcolours[l]
                     ,lwd=llwd[l],type=ltype,pch=16,cex=0.8)
            }#end for
            legend(x=legpos,inset=inset,legend=description,col=lcolours,lwd=llwd,cex=0.8)
            if (outform[o] == "x11"){
               locator(n=1)
               dev.off()
            }else{
               dev.off()
            }#end if
         } #end for outform
      }#end if plotit
   }#end for ntser
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #   Plot the climatology of the mean diurnal cycle.                                     #
   #---------------------------------------------------------------------------------------#
   cat("      * Plot some climatology of diurnal cycle...","\n")
   for (hh in 1:nclim){

      #----- Retrieve variable information from the list. ---------------------------------#
      climnow      = clim[[hh]]
      vnames       = climnow$vnam  
      description  = climnow$desc  
      lcolours     = climnow$colour
      llwd         = climnow$lwd
      ltype        = climnow$type
      plog         = climnow$plog
      prefix       = climnow$prefix
      theme        = climnow$theme 
      unit         = climnow$unit  
      legpos       = climnow$legpos
      plotit       = climnow$plt   
   
      if (plotit){

         #---------------------------------------------------------------------------------#
         #    Check whether the time series directory exists.  If not, create it.          #
         #---------------------------------------------------------------------------------#
         outdir   = paste(outpref,"climdcyc",sep="/")
         if (! file.exists(outdir)) dir.create(outdir)
         outtheme = paste(outdir,prefix,sep="/")
         if (! file.exists(outtheme)) dir.create(outtheme)
         cat("      + ",description," diurnal cycle for several variables...","\n")


         #----- Define the number of layers. ----------------------------------------------#
         nlayers   = length(vnames)
         ylimit    = NULL
         for (l in 1:nlayers) ylimit = c(ylimit,dcyc12mn[[vnames[l]]])
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
            cat("        > ",theme," time series - ",namemon,"...","\n")

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
               thisvar = dcyc12mn[[vnames[1]]]
               thisvar = cbind(thisvar[,ndcycle],thisvar)

               letitre = paste(theme," - ",lieu,"\n"
                              ,"Mean diurnal cycle - ",namemon,sep="")

               plot(x=thisday,y=thisvar[pmon,],type="n",main=letitre,xlab="Time"
                   ,ylim=ylimit,ylab=paste("[",unit,"]",sep=""),log=plog,xaxt="n"
                   ,cex.main=cex.main)
               axis(side=1,at=dcycplot$levels,labels=dcycplot$labels,padj=dcycplot$padj)
               if (plotgrid){ 
                  abline(v=dcycplot$levels,h=axTicks(side=2),col="gray52",lty="solid")
               }#end if
               for (l in 1:nlayers){
                  thisvar = dcyc12mn[[vnames[l]]]
                  thisvar = cbind(thisvar[,ndcycle],thisvar)
                  points(x=thisday,y=thisvar[pmon,],col=lcolours[l]
                        ,lwd=llwd[l],type=ltype,pch=16,cex=0.8)
               }#end for
               legend(x=legpos,inset=inset,legend=description,col=lcolours,lwd=llwd)
               if (outform[o] == "x11"){
                  locator(n=1)
                  dev.off()
               }else{
                  dev.off()
               }#end if
            } #end for outform
         }#end for pmon
      }#end if plotit
   }#end for ntser
   #---------------------------------------------------------------------------------------#





   #---------------------------------------------------------------------------------------#
   #   Plot the climatology of the soil properties.                                        #
   #---------------------------------------------------------------------------------------#
   for (v in 1:nsoilclim){

      #----- Retrieve variable information from the list. ---------------------------------#
      thisclim    = soilclim[[v]]
      vnam        = thisclim$vnam
      description = thisclim$desc
      unit        = thisclim$unit
      vcscheme    = thisclim$csch
      pnlog       = thisclim$pnlog
      plotit      = thisclim$plt

      if (plotit){

         #---------------------------------------------------------------------------------#
         #     Check if the directory exists.  If not, create it.                          #
         #---------------------------------------------------------------------------------#
         outdir  =  paste(outpref,"soilclim",sep="/")
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
         vararr  = mont12mn[[vnam]]

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
            vlevels = pretty.log(x=vrange,n=ncolshov)
            vnlev   = length(vlevels)
         }else{
            vrange  = range(varbuff,na.rm=TRUE)
            vlevels = pretty(x=vrange,n=ncolshov)
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
            sombreado(x=monaxis,y=soilaxis,z=varbuff,levels=vlevels,nlevels=vnlev
                     ,color.palette=get(vcscheme)
                     ,plot.title=title(main=letitre,xlab="Month",ylab="Soil depth [m]"
                                      ,cex.main=0.7)
                     ,key.title=title(main=unit,cex.main=0.8)
                     ,key.log=pnlog
                     ,plot.axes={axis(side=1,at=monat,labels=monlab)
                                 axis(side=2,at=zat,labels=znice)
                                 if (hovgrid){
                                    abline(h=zat,v=monat,col="gray52",lty="dotted")
                                 }#end if hovgrid
                                }#end plot.axes
                     )

            if (outform[o] == "x11"){
               locator(n=1)
               dev.off()
            }else{
               dev.off()
            }#end if
         } #end for outform
      }#end if plotit
   }#end for nhov
   #---------------------------------------------------------------------------------------#





   #---------------------------------------------------------------------------------------#
   #   Plot the climatology of the soil properties.                                        #
   #---------------------------------------------------------------------------------------#
   for (sts in 1:nsoilts){

      #----- Retrieve variable information from the list. ---------------------------------#
      thissts    = soilts[[sts]]
      vnam        = thissts$vnam
      description = thissts$desc
      unit        = thissts$unit
      vcscheme    = thissts$csch
      pnlog       = thissts$pnlog
      plotit      = thissts$plt

      if (plotit){

         #---------------------------------------------------------------------------------#
         #     Check if the directory exists.  If not, create it.                          #
         #---------------------------------------------------------------------------------#
         outdir  =  paste(outpref,"soilts",sep="/")
         if (! file.exists(outdir)) dir.create(outdir)
         cat("      + Time series profile of ",description,"...","\n")

         #----- Find the number of rows and columns, and the axes. ------------------------#
         timeaxis  = datum$tomonth
         soilaxis  = slz
         nmon      = length(timeaxis)
         nsoil     = nzg

         #----- Convert the vector data into an array. ------------------------------------#
         vararr  = datum[[vnam]]

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
            vlevels = pretty.log(x=vrange,n=ncolshov)
            vnlev   = length(vlevels)
         }else{
            vrange  = range(varbuff,na.rm=TRUE)
            vlevels = pretty(x=vrange,n=ncolshov)
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
            sombreado(x=whenaxis,y=soilaxis,z=varbuff,levels=vlevels,nlevels=vnlev
                     ,color.palette=get(vcscheme)
                     ,plot.title=title(main=letitre,xlab="Month",ylab="Soil depth [m]"
                                      ,cex.main=0.7)
                     ,key.title=title(main=unit,cex.main=0.8)
                     ,key.log=pnlog
                     ,plot.axes={axis(side=1,at=whenplot6$levels
                                     ,labels=whenplot6$labels,padj=whenplot6$padj)
                                 axis(side=2,at=zat,labels=znice)
                                 if (hovgrid){
                                    abline(h=zat,v=whenplot6$levels,col="gray52"
                                          ,lty="dotted")
                                 }#end if hovgrid
                                }#end plot.axes
                     )#end sombreado

            if (outform[o] == "x11"){
               locator(n=1)
               dev.off()
            }else{
               dev.off()
            }#end if
         } #end for outform
      }#end if plotit
   }#end for nhov
   #---------------------------------------------------------------------------------------#





   #---------------------------------------------------------------------------------------#
   #   Plot the Hovmoller looking diagrams showing months and years.                       #
   #---------------------------------------------------------------------------------------#
   for (v in 1:nhov){

      #----- Retrieve variable information from the list. ---------------------------------#
      thishovdi   = hovdi[[v]]
      vnam        = thishovdi$vnam
      description = thishovdi$desc
      unit        = thishovdi$unit
      vcscheme    = thishovdi$csch
      plotit      = thishovdi$plt

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
         outdir  =  paste(outpref,"hovmoller",sep="/")
         if (! file.exists(outdir)) dir.create(outdir)
         cat("      + ",description," Hovmoller time series ...","\n")

         #----- Load this variable into "thisvar". ----------------------------------------#
         thisvar = datum[[vnam]]

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
         vlevels = pretty(x=vrange,n=ncolshov)
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
            sombreado(x=monaxis,y=yraxis,z=varbuff,levels=vlevels,nlevels=vnlev
                     ,color.palette=get(vcscheme)
                     ,plot.title=title(main=letitre,xlab="Month",ylab="Year",cex.main=0.7)
                     ,key.title=title(main=unit,cex.main=0.8)
                     ,plot.axes={axis(side=1,at=monat,labels=monlab)
                                 axis(side=2,at=yrat)
                                 if (hovgrid){
                                    for (yl in yrat){
                                       abline(h=yl,col="gray52",lty="dotted")
                                    } #end for yl
                                    for (ml in monat){
                                       abline(v=ml,col="gray52",lty="dotted")
                                    } #end for ml
                                 }#end if hovgrid
                                }#end plot.axes
                     )

            if (outform[o] == "x11"){
               locator(n=1)
               dev.off()
            }else{
               dev.off()
            }#end if
         } #end for outform
      }#end if plotit
   }#end for nhov
   #---------------------------------------------------------------------------------------#






   #---------------------------------------------------------------------------------------#
   #   Plot the Hovmoller looking diagrams showing time of day and time.                   #
   #---------------------------------------------------------------------------------------#

   for (v in 1:nhdcyc){

      #----- Retrieve variable information from the list. ---------------------------------#
      thishdcyc   = hdcyc[[v]]
      vnam        = thishdcyc$vnam
      description = thishdcyc$desc
      unit        = thishdcyc$unit
      vcscheme    = thishdcyc$csch
      plotit      = thishdcyc$plt

      if (plotit){

         #---------------------------------------------------------------------------------#
         #     Check if the directory exists.  If not, create it.                          #
         #---------------------------------------------------------------------------------#
         outdir  =  paste(outpref,"hovdcycle",sep="/")
         if (! file.exists(outdir)) dir.create(outdir)
         cat("      + ",description," time series of diurnal cycle...","\n")

         #----- Load this variable into "thisvar". ----------------------------------------#
         vararr   = dcycmean[[vnam]]

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
         hdcycplot = pretty.time(whenaxis,n=8)

         vrange  = range(varbuff,na.rm=TRUE)
         vlevels = pretty(x=vrange,n=ncolshov)
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
            sombreado(x=whenaxis,y=hraxis,z=varbuff,levels=vlevels,nlevels=vnlev
                     ,color.palette=get(vcscheme)
                     ,plot.title=title(main=letitre,ylab="Time of day (GMT)"
                                      ,xlab="Time",cex.main=0.7)
                     ,key.title=title(main=unit,cex.main=0.8)
                     ,plot.axes={axis(side=1,at=hdcycplot$level,labels=hdcycplot$labels)
                                 axis(side=2,at=dcycplot$levels,labels=dcycplot$labels)
                                 if (hovgrid){
                                    abline(v=hdcycplot$levels,h=dcycplot$levels
                                          ,col="gray52",lty="dotted")
                                 }#end if hovgrid
                                }#end plot.axes
                     )

            if (outform[o] == "x11"){
               locator(n=1)
               dev.off()
            }else{
               dev.off()
            }#end if
         } #end for outform
      }#end if plotit
   }#end for nhov
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #   Plot the monthly boxplots.                                                          #
   #---------------------------------------------------------------------------------------#
   for (v in 1:nbox){

      #----- Retrieve variable information from the list. ---------------------------------#
      thisbplot   = pm.bplot[[v]]
      vnam        = thisbplot$vnam
      description = thisbplot$desc
      unit        = thisbplot$unit
      plotit      = thisbplot$plt

      if (plotit){
         #---------------------------------------------------------------------------------#
         #     Check if the directory exists.  If not, create it.                          #
         #---------------------------------------------------------------------------------#
         outdir  =  paste(outpref,"boxplot",sep="/")
         if (! file.exists(outdir)) dir.create(outdir)
         cat("      + ",description," box plot...","\n")

         #----- Load this variable into "thisvar". ----------------------------------------#
         thisvar = datum[[vnam]]

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
            plot(mmonth,thisvar,main=letitre,ylim=ylimit,cex.main=0.7
                ,xlab="Time",ylab=paste("[",unit,"]",sep=""))

            if (outform[o] == "x11"){
               locator(n=1)
               dev.off()
            }else{
               dev.off()
            }#end if
         } #end for outform
      }#end if
   }#end for nbox
   #---------------------------------------------------------------------------------------#






   #---------------------------------------------------------------------------------------#
   #      Bar plot by DBH class.                                                           #
   #---------------------------------------------------------------------------------------#
   cat("    + Bar plot by DBH classes...","\n")
   monbplot    = which(nummonths(datum$tomonth) %in% sasmonth)
   nmonbplot   = length(monbplot)
   pftuse      = which(apply(X=datum$nplantpftdbh,MARGIN=3,FUN=sum,na.rm=TRUE) > 0.)
   pftuse      = pftuse[pftuse != (npft+1)]
   npftuse     = length(pftuse)
   pftname.use = pft$name  [pftuse]
   pftcol.use  = pft$colour[pftuse]
   for (v in 1:nbarplotdbh){
      #----- Load settings for this variable.----------------------------------------------#
      thisbar     = barplotdbh[[v]]
      vnam        = thisbar$vnam
      description = thisbar$desc
      unit        = thisbar$unit
      stacked     = thisbar$stack
      plotit      = thisbar$plt
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
         thisvnam                  = datum[[vnam]][monbplot,,]
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
         barplotdir = paste(outpref,"barplotdbh",sep="/")
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
               barplot(height=t(thisvnam[m,,]),names.arg=dbhnames[1:ndbh],width=1.0
                      ,main=letitre,xlab=lexlab,ylab=leylab,ylim=ylimit,legend.text=FALSE
                      ,beside=(! stacked),col=pftcol.use,log=xylog
                      ,border="gray23",xpd=FALSE,cex.main=cex.main)
               if (plotgrid & (! stacked)){
                  xgrid=0.5+(1:ndbh)*(1+npftuse)
                  abline(v=xgrid,col="gray46",lty="solid")
               }#end if
               box()
               legend(x="topleft",inset=inset,legend=pftname.use,fill=pftcol.use
                     ,ncol=1,title="PFT",cex=1.0,bg="white")
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
   for (v in 1:npsas){
      #----- Retrieve variable information from the list. ---------------------------------#
      thissas     = sasplot[[v]]
      vnam        = thissas$vnam
      description = thissas$desc
      unit        = thissas$unit
      plotit      = thissas$plt

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
         varco =  datum[[vnam]]


         for (ww in names(datum$ageco)){

            #----- Find which year we are plotting. ---------------------------------------#
            cmonth   = substring(ww,7,8)
            thisyear = substring(ww,2,5)
            mm       = as.numeric(cmonth)
            yy       = as.numeric(thisyear)

            #----- Retrieve variable list, age, DBH, and PFT for this year. ---------------#
            ageww   = datum$ageco[[ww]]
            dbhww   = datum$dbhco[[ww]]
            pftww   = datum$pftco[[ww]]
            varww   = datum$varco[[ww]]
            popww   = datum$nplantco[[ww]] * datum$areaco[[ww]]

            #------------------------------------------------------------------------------#
            #     We only plot the SAS figures when the polygon is not an absolute desert. #
            #------------------------------------------------------------------------------#
            if (any (! is.na(varww))){
               #---------------------------------------------------------------------------#
               #      Find the range.  If the user wants the range to be fixed, then use   #
               # the global range, otherwise, simply use the range for this year.          #
               #---------------------------------------------------------------------------#
               if (sasfixlimits){
                  xlimit  = range(datum$ageco                  , na.rm=TRUE)
                  ylimit  = range(datum$dbhco                  , na.rm=TRUE)
                  zlimit  = range(datum$varco                  , na.rm=TRUE)
                  popmin  = min  (datum$nplantco * datum$areaco, na.rm=TRUE)
                  popmax  = max  (datum$nplantco * datum$areaco, na.rm=TRUE)
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

               pftin   = sort(unique(pftco[[ww]]))
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
                  pout = persp(x=ageaxis,y=dbhaxis,z=flooraxis,xlim=xlimit,ylim=ylimit
                              ,zlim=zlimit,theta=theta,phi=phi,col=gcol,expand=expz
                              ,ticktype="detailed",border=NA,xlab="Gap age [yr]"
                              ,ylab="DBH [cm]",zlab=lezlab,shade=shade,ltheta=ltheta
                              ,main=letitre,cex.main=0.7)
                  #----- Second plot, the actual data (aka my lollipop trees). ------------#
                  lines (trans3d(x=ageww,y=dbhww,z=varww,pmat=pout),type="l"
                        ,col="gray29",lwd=2)
                  points(trans3d(x=ageww,y=dbhww,z=varww,pmat=pout),type="p"
                        ,pch=pchww,col=colww,cex=cexww)
                  legend(x="bottomright",inset=inset,legend=pftleg,fill=colleg
                        ,ncol=1,bg="white",cex=0.9)


                  if (outform[o] == "x11"){
                     locator(n=1)
                     dev.off()
                  }else{
                     dev.off()
                  }#end if
               } #end for outform
            }#end if is.na(varww)
         }#end for nameco
      } #end if
   }#end for npsas
   #---------------------------------------------------------------------------------------#
}#end for places

#q("no")
