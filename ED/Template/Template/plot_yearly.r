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
fcgrid         = TRUE                   # Include a grid on the filled contour plots?
ncolshov       = 200                    # Target number of colours for Hovmoller diagrams.
hovgrid        = TRUE                   # Include a grid on the Hovmoller plots?
legwhere       = "topleft"              # Where should I place the legend?
inset          = 0.01                   # Inset between legend and edge of plot region.
legbg          = "white"                # Legend background colour.
scalleg        = 0.40                   # Expand y limits by this relative amount to fit
                                        #    the legend
cex.main       = 0.8                    # Scale coefficient for the title
ylnudge        = 0.05                  # Nudging factor for ylimit
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
obsrfile = paste(srcdir,"LBA_MIP.v8.RData",sep="/")
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
   outmain = paste(outroot,place,sep="/")
   outpref = paste(outmain,"yearly",sep="/")
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
   if (yearend - yearbeg + 1 <= nyears.long){
      sasmonth = sasmonth.short
   }else{
      sasmonth = sasmonth.long
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
   yfac     = datum$year
   dcycmean = datum$dcycmean
   dcycmsqu = datum$dcycmsqu
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Consolidate the yearly means for the long-term dynamics (the PFT and DBH/PFT      #
   # stuff).                                                                               #
   #---------------------------------------------------------------------------------------#
   cat ("    - Finding the annual statistics for multi-dimensional variables...","\n")
   cat ("      * Aggregating the annual mean of PFT-DBH variables...","\n")
   datum$agbpftdbh       = qapply(X=datum$agbpftdbh      ,DIM=1,INDEX=yfac,FUN=mean,na.rm=T)
   datum$basareapftdbh   = qapply(X=datum$basareapftdbh  ,DIM=1,INDEX=yfac,FUN=mean,na.rm=T)
   datum$laipftdbh       = qapply(X=datum$laipftdbh      ,DIM=1,INDEX=yfac,FUN=mean,na.rm=T)
   datum$waipftdbh       = qapply(X=datum$waipftdbh      ,DIM=1,INDEX=yfac,FUN=mean,na.rm=T)
   datum$taipftdbh       = qapply(X=datum$taipftdbh      ,DIM=1,INDEX=yfac,FUN=mean,na.rm=T)
   datum$gpppftdbh       = qapply(X=datum$gpppftdbh      ,DIM=1,INDEX=yfac,FUN=mean,na.rm=T)
   datum$npppftdbh       = qapply(X=datum$npppftdbh      ,DIM=1,INDEX=yfac,FUN=mean,na.rm=T)
   datum$mcopftdbh       = qapply(X=datum$mcopftdbh      ,DIM=1,INDEX=yfac,FUN=mean,na.rm=T)
   datum$cbapftdbh       = qapply(X=datum$cbapftdbh      ,DIM=1,INDEX=yfac,FUN=mean,na.rm=T)
   datum$cbalightpftdbh  = qapply(X=datum$cbalightpftdbh ,DIM=1,INDEX=yfac,FUN=mean,na.rm=T)
   datum$cbamoistpftdbh  = qapply(X=datum$cbamoistpftdbh ,DIM=1,INDEX=yfac,FUN=mean,na.rm=T)
   datum$cbarelpftdbh    = qapply(X=datum$cbarelpftdbh   ,DIM=1,INDEX=yfac,FUN=mean,na.rm=T)
   datum$ldroppftdbh     = qapply(X=datum$ldroppftdbh    ,DIM=1,INDEX=yfac,FUN=mean,na.rm=T)
   datum$fsopftdbh       = qapply(X=datum$fsopftdbh      ,DIM=1,INDEX=yfac,FUN=mean,na.rm=T)
   datum$demandpftdbh    = qapply(X=datum$demandpftdbh   ,DIM=1,INDEX=yfac,FUN=mean,na.rm=T)
   datum$supplypftdbh    = qapply(X=datum$supplypftdbh   ,DIM=1,INDEX=yfac,FUN=mean,na.rm=T)
   datum$nplantpftdbh    = qapply(X=datum$nplantpftdbh   ,DIM=1,INDEX=yfac,FUN=mean,na.rm=T)
   datum$mortpftdbh      = qapply(X=datum$mortpftdbh     ,DIM=1,INDEX=yfac,FUN=mean,na.rm=T)
   datum$agemortpftdbh   = qapply(X=datum$agemortpftdbh  ,DIM=1,INDEX=yfac,FUN=mean,na.rm=T)
   datum$ncbmortpftdbh   = qapply(X=datum$ncbmortpftdbh  ,DIM=1,INDEX=yfac,FUN=mean,na.rm=T)
   datum$tfallmortpftdbh = qapply(X=datum$tfallmortpftdbh,DIM=1,INDEX=yfac,FUN=mean,na.rm=T)
   datum$coldmortpftdbh  = qapply(X=datum$coldmortpftdbh ,DIM=1,INDEX=yfac,FUN=mean,na.rm=T)
   datum$distmortpftdbh  = qapply(X=datum$distmortpftdbh ,DIM=1,INDEX=yfac,FUN=mean,na.rm=T)
   datum$growthpftdbh    = qapply(X=datum$growthpftdbh   ,DIM=1,INDEX=yfac,FUN=mean,na.rm=T)
   #----- PFT arrays.   The "+1" column contains the total. -------------------------------#
   cat ("      * Aggregating the annual mean of PFT variables...","\n")
   datum$agbpft          = qapply(X=datum$agbpft         ,DIM=1,INDEX=yfac,FUN=mean,na.rm=T)
   datum$bseedspft       = qapply(X=datum$bseedspft      ,DIM=1,INDEX=yfac,FUN=mean,na.rm=T)
   datum$nplantpft       = qapply(X=datum$nplantpft      ,DIM=1,INDEX=yfac,FUN=mean,na.rm=T)
   datum$laipft          = qapply(X=datum$laipft         ,DIM=1,INDEX=yfac,FUN=mean,na.rm=T)
   datum$waipft          = qapply(X=datum$waipft         ,DIM=1,INDEX=yfac,FUN=mean,na.rm=T)
   datum$taipft          = qapply(X=datum$taipft         ,DIM=1,INDEX=yfac,FUN=mean,na.rm=T)
   datum$gpppft          = qapply(X=datum$gpppft         ,DIM=1,INDEX=yfac,FUN=mean,na.rm=T)
   datum$npppft          = qapply(X=datum$npppft         ,DIM=1,INDEX=yfac,FUN=mean,na.rm=T)
   datum$mcopft          = qapply(X=datum$mcopft         ,DIM=1,INDEX=yfac,FUN=mean,na.rm=T)
   datum$cbapft          = qapply(X=datum$cbapft         ,DIM=1,INDEX=yfac,FUN=mean,na.rm=T)
   datum$cbalightpft     = qapply(X=datum$cbalightpft    ,DIM=1,INDEX=yfac,FUN=mean,na.rm=T)
   datum$cbamoistpft     = qapply(X=datum$cbamoistpft    ,DIM=1,INDEX=yfac,FUN=mean,na.rm=T)
   datum$cbarelpft       = qapply(X=datum$cbarelpft      ,DIM=1,INDEX=yfac,FUN=mean,na.rm=T)
   datum$ldroppft        = qapply(X=datum$ldroppft       ,DIM=1,INDEX=yfac,FUN=mean,na.rm=T)
   datum$fsopft          = qapply(X=datum$fsopft         ,DIM=1,INDEX=yfac,FUN=mean,na.rm=T)
   datum$balivepft       = qapply(X=datum$balivepft      ,DIM=1,INDEX=yfac,FUN=mean,na.rm=T)
   datum$bdeadpft        = qapply(X=datum$bdeadpft       ,DIM=1,INDEX=yfac,FUN=mean,na.rm=T)
   datum$bleafpft        = qapply(X=datum$bleafpft       ,DIM=1,INDEX=yfac,FUN=mean,na.rm=T)
   datum$brootpft        = qapply(X=datum$brootpft       ,DIM=1,INDEX=yfac,FUN=mean,na.rm=T)
   datum$bswoodpft       = qapply(X=datum$bswoodpft      ,DIM=1,INDEX=yfac,FUN=mean,na.rm=T)
   datum$bstorepft       = qapply(X=datum$bstorepft      ,DIM=1,INDEX=yfac,FUN=mean,na.rm=T)
   datum$basareapft      = qapply(X=datum$basareapft     ,DIM=1,INDEX=yfac,FUN=mean,na.rm=T)
   datum$leafresppft     = qapply(X=datum$leafresppft    ,DIM=1,INDEX=yfac,FUN=mean,na.rm=T)
   datum$rootresppft     = qapply(X=datum$rootresppft    ,DIM=1,INDEX=yfac,FUN=mean,na.rm=T)
   datum$growthresppft   = qapply(X=datum$growthresppft  ,DIM=1,INDEX=yfac,FUN=mean,na.rm=T)
   datum$mortpft         = qapply(X=datum$mortpft        ,DIM=1,INDEX=yfac,FUN=mean,na.rm=T)
   datum$agemortpft      = qapply(X=datum$agemortpft     ,DIM=1,INDEX=yfac,FUN=mean,na.rm=T)
   datum$ncbmortpft      = qapply(X=datum$ncbmortpft     ,DIM=1,INDEX=yfac,FUN=mean,na.rm=T)
   datum$tfallmortpft    = qapply(X=datum$tfallmortpft   ,DIM=1,INDEX=yfac,FUN=mean,na.rm=T)
   datum$coldmortpft     = qapply(X=datum$coldmortpft    ,DIM=1,INDEX=yfac,FUN=mean,na.rm=T)
   datum$distmortpft     = qapply(X=datum$distmortpft    ,DIM=1,INDEX=yfac,FUN=mean,na.rm=T)
   datum$recrpft         = qapply(X=datum$recrpft        ,DIM=1,INDEX=yfac,FUN=mean,na.rm=T)
   datum$growthpft       = qapply(X=datum$growthpft      ,DIM=1,INDEX=yfac,FUN=mean,na.rm=T)
   datum$censusnplantpft = qapply(X=datum$censusnplantpft,DIM=1,INDEX=yfac,FUN=mean,na.rm=T)
   datum$censuslaipft    = qapply(X=datum$censuslaipft   ,DIM=1,INDEX=yfac,FUN=mean,na.rm=T)
   datum$censuswaipft    = qapply(X=datum$censuswaipft   ,DIM=1,INDEX=yfac,FUN=mean,na.rm=T)
   datum$censustaipft    = qapply(X=datum$censustaipft   ,DIM=1,INDEX=yfac,FUN=mean,na.rm=T)
   datum$censusagbpft    = qapply(X=datum$censusagbpft   ,DIM=1,INDEX=yfac,FUN=mean,na.rm=T)
   datum$censusbapft     = qapply(X=datum$censusbapft    ,DIM=1,INDEX=yfac,FUN=mean,na.rm=T)
   #----- LU arrays.   The "+1" column contains the total. --------------------------------#
   cat ("      * Aggregating the annual mean of LU variables...","\n")
   datum$agblu           = qapply(X=datum$agblu          ,DIM=1,INDEX=yfac,FUN=mean,na.rm=T)
   datum$lailu           = qapply(X=datum$lailu          ,DIM=1,INDEX=yfac,FUN=mean,na.rm=T)
   datum$gpplu           = qapply(X=datum$gpplu          ,DIM=1,INDEX=yfac,FUN=mean,na.rm=T)
   datum$npplu           = qapply(X=datum$npplu          ,DIM=1,INDEX=yfac,FUN=mean,na.rm=T)
   datum$arealu          = qapply(X=datum$arealu         ,DIM=1,INDEX=yfac,FUN=mean,na.rm=T)
   datum$basarealu       = qapply(X=datum$basarealu      ,DIM=1,INDEX=yfac,FUN=mean,na.rm=T)
   #----- Miscellaneous arrays. -----------------------------------------------------------#
   cat ("      * Aggregating the annual mean of DIST variables...","\n")
   datum$dist            = qapply(X=datum$dist           ,DIM=1,INDEX=yfac,FUN=mean,na.rm=T)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Here we find the monthly means for month, then compute the standard deviation.   #
   #---------------------------------------------------------------------------------------#
   cat ("    - Finding the monthly mean...","\n")
   cat ("      * Aggregating the monthly mean...","\n")
   mont12mn               = list()
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
   cat ("      * Standard deviation of the monthly means...","\n")
   mont12sd               = list()
   mont12sd$gpp         = tapply(X=datum$gpp          ,INDEX=mfac      ,FUN=sd  ,na.rm=T)
   mont12sd$npp         = tapply(X=datum$npp          ,INDEX=mfac      ,FUN=sd  ,na.rm=T)
   mont12sd$nep         = tapply(X=datum$nep          ,INDEX=mfac      ,FUN=sd  ,na.rm=T)
   mont12sd$plresp      = tapply(X=datum$plresp       ,INDEX=mfac      ,FUN=sd  ,na.rm=T)
   mont12sd$leaf.resp   = tapply(X=datum$leaf.resp    ,INDEX=mfac      ,FUN=sd  ,na.rm=T)
   mont12sd$root.resp   = tapply(X=datum$root.resp    ,INDEX=mfac      ,FUN=sd  ,na.rm=T)
   mont12sd$growth.resp = tapply(X=datum$growth.resp  ,INDEX=mfac      ,FUN=sd  ,na.rm=T)
   mont12sd$hetresp     = tapply(X=datum$hetresp      ,INDEX=mfac      ,FUN=sd  ,na.rm=T)
   mont12sd$cwdresp     = tapply(X=datum$cwdresp      ,INDEX=mfac      ,FUN=sd  ,na.rm=T)
   mont12sd$cflxca      = tapply(X=datum$cflxca       ,INDEX=mfac      ,FUN=sd  ,na.rm=T)
   mont12sd$cflxst      = tapply(X=datum$cflxst       ,INDEX=mfac      ,FUN=sd  ,na.rm=T)
   mont12sd$nee         = tapply(X=datum$nee          ,INDEX=mfac      ,FUN=sd  ,na.rm=T)
   mont12sd$reco        = tapply(X=datum$reco         ,INDEX=mfac      ,FUN=sd  ,na.rm=T)
   mont12sd$hflxca      = tapply(X=datum$hflxca       ,INDEX=mfac      ,FUN=sd  ,na.rm=T)
   mont12sd$hflxlc      = tapply(X=datum$hflxlc       ,INDEX=mfac      ,FUN=sd  ,na.rm=T)
   mont12sd$hflxwc      = tapply(X=datum$hflxwc       ,INDEX=mfac      ,FUN=sd  ,na.rm=T)
   mont12sd$hflxgc      = tapply(X=datum$hflxgc       ,INDEX=mfac      ,FUN=sd  ,na.rm=T)
   mont12sd$wflxca      = tapply(X=datum$wflxca       ,INDEX=mfac      ,FUN=sd  ,na.rm=T)
   mont12sd$qwflxca     = tapply(X=datum$qwflxca      ,INDEX=mfac      ,FUN=sd  ,na.rm=T)
   mont12sd$wflxlc      = tapply(X=datum$wflxlc       ,INDEX=mfac      ,FUN=sd  ,na.rm=T)
   mont12sd$wflxwc      = tapply(X=datum$wflxwc       ,INDEX=mfac      ,FUN=sd  ,na.rm=T)
   mont12sd$wflxgc      = tapply(X=datum$wflxgc       ,INDEX=mfac      ,FUN=sd  ,na.rm=T)
   mont12sd$evap        = tapply(X=datum$evap         ,INDEX=mfac      ,FUN=sd  ,na.rm=T)
   mont12sd$transp      = tapply(X=datum$transp       ,INDEX=mfac      ,FUN=sd  ,na.rm=T)
   mont12sd$rain        = tapply(X=datum$rain         ,INDEX=mfac      ,FUN=sd  ,na.rm=T)
   mont12sd$atm.temp    = tapply(X=datum$atm.temp     ,INDEX=mfac      ,FUN=sd  ,na.rm=T)
   mont12sd$rshort      = tapply(X=datum$rshort       ,INDEX=mfac      ,FUN=sd  ,na.rm=T)
   mont12sd$rshortup    = tapply(X=datum$rshortup     ,INDEX=mfac      ,FUN=sd  ,na.rm=T)
   mont12sd$rlong       = tapply(X=datum$rlong        ,INDEX=mfac      ,FUN=sd  ,na.rm=T)
   mont12sd$rlongup     = tapply(X=datum$rlongup      ,INDEX=mfac      ,FUN=sd  ,na.rm=T)
   # mont12sd$par.tot     = tapply(X=datum$par.tot      ,INDEX=mfac      ,FUN=sd  ,na.rm=T)
   mont12sd$parup       = tapply(X=datum$parup        ,INDEX=mfac      ,FUN=sd  ,na.rm=T)
   mont12sd$rnet        = tapply(X=datum$rnet         ,INDEX=mfac      ,FUN=sd  ,na.rm=T)
   mont12sd$albedo      = tapply(X=datum$albedo       ,INDEX=mfac      ,FUN=sd  ,na.rm=T)
   mont12sd$atm.shv     = tapply(X=datum$atm.shv      ,INDEX=mfac      ,FUN=sd  ,na.rm=T)
   mont12sd$atm.co2     = tapply(X=datum$atm.co2      ,INDEX=mfac      ,FUN=sd  ,na.rm=T)
   mont12sd$atm.prss    = tapply(X=datum$atm.prss     ,INDEX=mfac      ,FUN=sd  ,na.rm=T)
   mont12sd$atm.vels    = tapply(X=datum$atm.vels     ,INDEX=mfac      ,FUN=sd  ,na.rm=T)
   mont12sd$ustar       = tapply(X=datum$ustar        ,INDEX=mfac      ,FUN=sd  ,na.rm=T)
   mont12sd$soil.temp   = qapply(X=datum$soil.temp    ,INDEX=mfac,DIM=1,FUN=sd  ,na.rm=T)
   mont12sd$soil.water  = qapply(X=datum$soil.water   ,INDEX=mfac,DIM=1,FUN=sd  ,na.rm=T)
   mont12sd$soil.mstpot = qapply(X=datum$soil.mstpot  ,INDEX=mfac,DIM=1,FUN=sd  ,na.rm=T)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Here we find the monthly means for month, then compute the standard deviation.   #
   #---------------------------------------------------------------------------------------#
   cat ("    - Finding the annual mean...","\n")
   cat ("      * Aggregating the annual mean...","\n")
   year12mn             = list()
   year12mn$gpp         = tapply(X=datum$gpp          ,INDEX=yfac      ,FUN=mean,na.rm=T)
   year12mn$npp         = tapply(X=datum$npp          ,INDEX=yfac      ,FUN=mean,na.rm=T)
   year12mn$nep         = tapply(X=datum$nep          ,INDEX=yfac      ,FUN=mean,na.rm=T)
   year12mn$plresp      = tapply(X=datum$plresp       ,INDEX=yfac      ,FUN=mean,na.rm=T)
   year12mn$leaf.resp   = tapply(X=datum$leaf.resp    ,INDEX=yfac      ,FUN=mean,na.rm=T)
   year12mn$root.resp   = tapply(X=datum$root.resp    ,INDEX=yfac      ,FUN=mean,na.rm=T)
   year12mn$growth.resp = tapply(X=datum$growth.resp  ,INDEX=yfac      ,FUN=mean,na.rm=T)
   year12mn$hetresp     = tapply(X=datum$hetresp      ,INDEX=yfac      ,FUN=mean,na.rm=T)
   year12mn$cwdresp     = tapply(X=datum$cwdresp      ,INDEX=yfac      ,FUN=mean,na.rm=T)
   year12mn$cflxca      = tapply(X=datum$cflxca       ,INDEX=yfac      ,FUN=mean,na.rm=T)
   year12mn$cflxst      = tapply(X=datum$cflxst       ,INDEX=yfac      ,FUN=mean,na.rm=T)
   year12mn$nee         = tapply(X=datum$nee          ,INDEX=yfac      ,FUN=mean,na.rm=T)
   year12mn$reco        = tapply(X=datum$reco         ,INDEX=yfac      ,FUN=mean,na.rm=T)
   year12mn$hflxca      = tapply(X=datum$hflxca       ,INDEX=yfac      ,FUN=mean,na.rm=T)
   year12mn$hflxlc      = tapply(X=datum$hflxlc       ,INDEX=yfac      ,FUN=mean,na.rm=T)
   year12mn$hflxwc      = tapply(X=datum$hflxwc       ,INDEX=yfac      ,FUN=mean,na.rm=T)
   year12mn$hflxgc      = tapply(X=datum$hflxgc       ,INDEX=yfac      ,FUN=mean,na.rm=T)
   year12mn$wflxca      = tapply(X=datum$wflxca       ,INDEX=yfac      ,FUN=mean,na.rm=T)
   year12mn$qwflxca     = tapply(X=datum$qwflxca      ,INDEX=yfac      ,FUN=mean,na.rm=T)
   year12mn$wflxlc      = tapply(X=datum$wflxlc       ,INDEX=yfac      ,FUN=mean,na.rm=T)
   year12mn$wflxwc      = tapply(X=datum$wflxwc       ,INDEX=yfac      ,FUN=mean,na.rm=T)
   year12mn$wflxgc      = tapply(X=datum$wflxgc       ,INDEX=yfac      ,FUN=mean,na.rm=T)
   year12mn$evap        = tapply(X=datum$evap         ,INDEX=yfac      ,FUN=mean,na.rm=T)
   year12mn$transp      = tapply(X=datum$transp       ,INDEX=yfac      ,FUN=mean,na.rm=T)
   year12mn$rain        = tapply(X=datum$rain         ,INDEX=yfac      ,FUN=sum ,na.rm=T)
   year12mn$atm.temp    = tapply(X=datum$atm.temp     ,INDEX=yfac      ,FUN=mean,na.rm=T)
   year12mn$rshort      = tapply(X=datum$rshort       ,INDEX=yfac      ,FUN=mean,na.rm=T)
   year12mn$rshortup    = tapply(X=datum$rshortup     ,INDEX=yfac      ,FUN=mean,na.rm=T)
   year12mn$rlong       = tapply(X=datum$rlong        ,INDEX=yfac      ,FUN=mean,na.rm=T)
   year12mn$rlongup     = tapply(X=datum$rlongup      ,INDEX=yfac      ,FUN=mean,na.rm=T)
   # year12mn$par.tot     = tapply(X=datum$par.tot      ,INDEX=yfac      ,FUN=mean,na.rm=T)
   year12mn$parup       = tapply(X=datum$parup        ,INDEX=yfac      ,FUN=mean,na.rm=T)
   year12mn$rnet        = tapply(X=datum$rnet         ,INDEX=yfac      ,FUN=mean,na.rm=T)
   year12mn$albedo      = tapply(X=datum$albedo       ,INDEX=yfac      ,FUN=mean,na.rm=T)
   year12mn$atm.shv     = tapply(X=datum$atm.shv      ,INDEX=yfac      ,FUN=mean,na.rm=T)
   year12mn$atm.co2     = tapply(X=datum$atm.co2      ,INDEX=yfac      ,FUN=mean,na.rm=T)
   year12mn$atm.prss    = tapply(X=datum$atm.prss     ,INDEX=yfac      ,FUN=mean,na.rm=T)
   year12mn$atm.vels    = tapply(X=datum$atm.vels     ,INDEX=yfac      ,FUN=mean,na.rm=T)
   year12mn$ustar       = tapply(X=datum$ustar        ,INDEX=yfac      ,FUN=mean,na.rm=T)
   year12mn$soil.temp   = qapply(X=datum$soil.temp    ,INDEX=yfac,DIM=1,FUN=mean,na.rm=T)
   year12mn$soil.water  = qapply(X=datum$soil.water   ,INDEX=yfac,DIM=1,FUN=mean,na.rm=T)
   year12mn$soil.mstpot = qapply(X=datum$soil.mstpot  ,INDEX=yfac,DIM=1,FUN=mean,na.rm=T)
   cat ("      * Standard deviation of the annual mean...","\n")
   year12sd             = list()
   year12sd$gpp         = tapply(X=datum$gpp          ,INDEX=yfac      ,FUN=sd  ,na.rm=T)
   year12sd$npp         = tapply(X=datum$npp          ,INDEX=yfac      ,FUN=sd  ,na.rm=T)
   year12sd$nep         = tapply(X=datum$nep          ,INDEX=yfac      ,FUN=sd  ,na.rm=T)
   year12sd$plresp      = tapply(X=datum$plresp       ,INDEX=yfac      ,FUN=sd  ,na.rm=T)
   year12sd$leaf.resp   = tapply(X=datum$leaf.resp    ,INDEX=yfac      ,FUN=sd  ,na.rm=T)
   year12sd$root.resp   = tapply(X=datum$root.resp    ,INDEX=yfac      ,FUN=sd  ,na.rm=T)
   year12sd$growth.resp = tapply(X=datum$growth.resp  ,INDEX=yfac      ,FUN=sd  ,na.rm=T)
   year12sd$hetresp     = tapply(X=datum$hetresp      ,INDEX=yfac      ,FUN=sd  ,na.rm=T)
   year12sd$cwdresp     = tapply(X=datum$cwdresp      ,INDEX=yfac      ,FUN=sd  ,na.rm=T)
   year12sd$cflxca      = tapply(X=datum$cflxca       ,INDEX=yfac      ,FUN=sd  ,na.rm=T)
   year12sd$cflxst      = tapply(X=datum$cflxst       ,INDEX=yfac      ,FUN=sd  ,na.rm=T)
   year12sd$nee         = tapply(X=datum$nee          ,INDEX=yfac      ,FUN=sd  ,na.rm=T)
   year12sd$reco        = tapply(X=datum$reco         ,INDEX=yfac      ,FUN=sd  ,na.rm=T)
   year12sd$hflxca      = tapply(X=datum$hflxca       ,INDEX=yfac      ,FUN=sd  ,na.rm=T)
   year12sd$hflxlc      = tapply(X=datum$hflxlc       ,INDEX=yfac      ,FUN=sd  ,na.rm=T)
   year12sd$hflxwc      = tapply(X=datum$hflxwc       ,INDEX=yfac      ,FUN=sd  ,na.rm=T)
   year12sd$hflxgc      = tapply(X=datum$hflxgc       ,INDEX=yfac      ,FUN=sd  ,na.rm=T)
   year12sd$wflxca      = tapply(X=datum$wflxca       ,INDEX=yfac      ,FUN=sd  ,na.rm=T)
   year12sd$qwflxca     = tapply(X=datum$qwflxca      ,INDEX=yfac      ,FUN=sd  ,na.rm=T)
   year12sd$wflxlc      = tapply(X=datum$wflxlc       ,INDEX=yfac      ,FUN=sd  ,na.rm=T)
   year12sd$wflxwc      = tapply(X=datum$wflxwc       ,INDEX=yfac      ,FUN=sd  ,na.rm=T)
   year12sd$wflxgc      = tapply(X=datum$wflxgc       ,INDEX=yfac      ,FUN=sd  ,na.rm=T)
   year12sd$evap        = tapply(X=datum$evap         ,INDEX=yfac      ,FUN=sd  ,na.rm=T)
   year12sd$transp      = tapply(X=datum$transp       ,INDEX=yfac      ,FUN=sd  ,na.rm=T)
   year12sd$rain        = tapply(X=datum$rain*12      ,INDEX=yfac      ,FUN=sd  ,na.rm=T)
   year12sd$atm.temp    = tapply(X=datum$atm.temp     ,INDEX=yfac      ,FUN=sd  ,na.rm=T)
   year12sd$rshort      = tapply(X=datum$rshort       ,INDEX=yfac      ,FUN=sd  ,na.rm=T)
   year12sd$rshortup    = tapply(X=datum$rshortup     ,INDEX=yfac      ,FUN=sd  ,na.rm=T)
   year12sd$rlong       = tapply(X=datum$rlong        ,INDEX=yfac      ,FUN=sd  ,na.rm=T)
   year12sd$rlongup     = tapply(X=datum$rlongup      ,INDEX=yfac      ,FUN=sd  ,na.rm=T)
   # year12sd$par.tot     = tapply(X=datum$par.tot      ,INDEX=yfac      ,FUN=sd  ,na.rm=T)
   year12sd$parup       = tapply(X=datum$parup        ,INDEX=yfac      ,FUN=sd  ,na.rm=T)
   year12sd$rnet        = tapply(X=datum$rnet         ,INDEX=yfac      ,FUN=sd  ,na.rm=T)
   year12sd$albedo      = tapply(X=datum$albedo       ,INDEX=yfac      ,FUN=sd  ,na.rm=T)
   year12sd$atm.shv     = tapply(X=datum$atm.shv      ,INDEX=yfac      ,FUN=sd  ,na.rm=T)
   year12sd$atm.co2     = tapply(X=datum$atm.co2      ,INDEX=yfac      ,FUN=sd  ,na.rm=T)
   year12sd$atm.prss    = tapply(X=datum$atm.prss     ,INDEX=yfac      ,FUN=sd  ,na.rm=T)
   year12sd$atm.vels    = tapply(X=datum$atm.vels     ,INDEX=yfac      ,FUN=sd  ,na.rm=T)
   year12sd$ustar       = tapply(X=datum$ustar        ,INDEX=yfac      ,FUN=sd  ,na.rm=T)
   year12sd$soil.temp   = qapply(X=datum$soil.temp    ,INDEX=yfac,DIM=1,FUN=sd  ,na.rm=T)
   year12sd$soil.water  = qapply(X=datum$soil.water   ,INDEX=yfac,DIM=1,FUN=sd  ,na.rm=T)
   year12sd$soil.mstpot = qapply(X=datum$soil.mstpot  ,INDEX=yfac,DIM=1,FUN=sd  ,na.rm=T)
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #      Here we find the Mean diurnal cycle for each month, then compute the standard    #
   # deviation.                                                                            #
   #---------------------------------------------------------------------------------------#
   cat ("    - Aggregating the annual mean of the diurnal cycle...","\n")
   dcyc12mn             =list()
   dcyc12mn$gpp         =qapply(X=dcycmean$gpp         ,INDEX=yfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$npp         =qapply(X=dcycmean$npp         ,INDEX=yfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$plresp      =qapply(X=dcycmean$plresp      ,INDEX=yfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$leaf.resp   =qapply(X=dcycmean$leaf.resp   ,INDEX=yfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$root.resp   =qapply(X=dcycmean$root.resp   ,INDEX=yfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$hetresp     =qapply(X=dcycmean$hetresp     ,INDEX=yfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$cwdresp     =qapply(X=dcycmean$cwdresp     ,INDEX=yfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$nep         =qapply(X=dcycmean$nep         ,INDEX=yfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$nee         =qapply(X=dcycmean$nee         ,INDEX=yfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$reco        =qapply(X=dcycmean$reco        ,INDEX=yfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$cflxca      =qapply(X=dcycmean$cflxca      ,INDEX=yfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$cflxst      =qapply(X=dcycmean$cflxst      ,INDEX=yfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$hflxca      =qapply(X=dcycmean$hflxca      ,INDEX=yfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$hflxlc      =qapply(X=dcycmean$hflxlc      ,INDEX=yfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$hflxwc      =qapply(X=dcycmean$hflxwc      ,INDEX=yfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$hflxgc      =qapply(X=dcycmean$hflxgc      ,INDEX=yfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$wflxca      =qapply(X=dcycmean$wflxca      ,INDEX=yfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$qwflxca     =qapply(X=dcycmean$qwflxca     ,INDEX=yfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$wflxlc      =qapply(X=dcycmean$wflxlc      ,INDEX=yfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$wflxwc      =qapply(X=dcycmean$wflxwc      ,INDEX=yfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$wflxgc      =qapply(X=dcycmean$wflxgc      ,INDEX=yfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$evap        =qapply(X=dcycmean$evap        ,INDEX=yfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$transp      =qapply(X=dcycmean$transp      ,INDEX=yfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$atm.temp    =qapply(X=dcycmean$atm.temp    ,INDEX=yfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$can.temp    =qapply(X=dcycmean$can.temp    ,INDEX=yfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$leaf.temp   =qapply(X=dcycmean$leaf.temp   ,INDEX=yfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$wood.temp   =qapply(X=dcycmean$wood.temp   ,INDEX=yfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$gnd.temp    =qapply(X=dcycmean$gnd.temp    ,INDEX=yfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$atm.shv     =qapply(X=dcycmean$atm.shv     ,INDEX=yfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$can.shv     =qapply(X=dcycmean$can.shv     ,INDEX=yfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$gnd.shv     =qapply(X=dcycmean$gnd.shv     ,INDEX=yfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$atm.vpd     =qapply(X=dcycmean$atm.vpd     ,INDEX=yfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$can.vpd     =qapply(X=dcycmean$can.vpd     ,INDEX=yfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$leaf.vpd    =qapply(X=dcycmean$leaf.vpd    ,INDEX=yfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$atm.co2     =qapply(X=dcycmean$atm.co2     ,INDEX=yfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$can.co2     =qapply(X=dcycmean$can.co2     ,INDEX=yfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$atm.prss    =qapply(X=dcycmean$atm.prss    ,INDEX=yfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$can.prss    =qapply(X=dcycmean$can.prss    ,INDEX=yfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$atm.vels    =qapply(X=dcycmean$atm.vels    ,INDEX=yfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$ustar       =qapply(X=dcycmean$ustar       ,INDEX=yfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$fs.open     =qapply(X=dcycmean$fs.open     ,INDEX=yfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$rain        =qapply(X=dcycmean$rain        ,INDEX=yfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$rshort      =qapply(X=dcycmean$rshort      ,INDEX=yfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$rshort.beam =qapply(X=dcycmean$rshort.beam ,INDEX=yfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$rshort.diff =qapply(X=dcycmean$rshort.diff ,INDEX=yfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$rshort.gnd  =qapply(X=dcycmean$rshort.gnd  ,INDEX=yfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$rshortup    =qapply(X=dcycmean$rshortup    ,INDEX=yfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$rlong       =qapply(X=dcycmean$rlong       ,INDEX=yfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$rlong.gnd   =qapply(X=dcycmean$rlong.gnd   ,INDEX=yfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$rlongup     =qapply(X=dcycmean$rlongup     ,INDEX=yfac,DIM=1,FUN=mean,na.rm=T)
   # dcyc12mn$par.tot     =qapply(X=dcycmean$par.tot     ,INDEX=yfac,DIM=1,FUN=mean,na.rm=T)
   # dcyc12mn$par.beam    =qapply(X=dcycmean$par.beam    ,INDEX=yfac,DIM=1,FUN=mean,na.rm=T)
   # dcyc12mn$par.diff    =qapply(X=dcycmean$par.diff    ,INDEX=yfac,DIM=1,FUN=mean,na.rm=T)
   # dcyc12mn$par.gnd     =qapply(X=dcycmean$par.gnd     ,INDEX=yfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$parup       =qapply(X=dcycmean$parup       ,INDEX=yfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$rnet        =qapply(X=dcycmean$rnet        ,INDEX=yfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$albedo      =qapply(X=dcycmean$albedo      ,INDEX=yfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$albedo.beam =qapply(X=dcycmean$albedo.beam ,INDEX=yfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$albedo.diff =qapply(X=dcycmean$albedo.diff ,INDEX=yfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$rlong.albedo=qapply(X=dcycmean$rlong.albedo,INDEX=yfac,DIM=1,FUN=mean,na.rm=T)
   #----- Find the mean sum of squares. ---------------------------------------------------#
   dcyc12sd             =list()
   dcyc12sd$gpp         =qapply(X=dcycmean$gpp         ,INDEX=yfac,DIM=1,FUN=sd  ,na.rm=T)
   dcyc12sd$npp         =qapply(X=dcycmean$npp         ,INDEX=yfac,DIM=1,FUN=sd  ,na.rm=T)
   dcyc12sd$plresp      =qapply(X=dcycmean$plresp      ,INDEX=yfac,DIM=1,FUN=sd  ,na.rm=T)
   dcyc12sd$leaf.resp   =qapply(X=dcycmean$leaf.resp   ,INDEX=yfac,DIM=1,FUN=sd  ,na.rm=T)
   dcyc12sd$root.resp   =qapply(X=dcycmean$root.resp   ,INDEX=yfac,DIM=1,FUN=sd  ,na.rm=T)
   dcyc12sd$hetresp     =qapply(X=dcycmean$hetresp     ,INDEX=yfac,DIM=1,FUN=sd  ,na.rm=T)
   dcyc12sd$cwdresp     =qapply(X=dcycmean$cwdresp     ,INDEX=yfac,DIM=1,FUN=sd  ,na.rm=T)
   dcyc12sd$nep         =qapply(X=dcycmean$nep         ,INDEX=yfac,DIM=1,FUN=sd  ,na.rm=T)
   dcyc12sd$nee         =qapply(X=dcycmean$nee         ,INDEX=yfac,DIM=1,FUN=sd  ,na.rm=T)
   dcyc12sd$reco        =qapply(X=dcycmean$reco        ,INDEX=yfac,DIM=1,FUN=sd  ,na.rm=T)
   dcyc12sd$cflxca      =qapply(X=dcycmean$cflxca      ,INDEX=yfac,DIM=1,FUN=sd  ,na.rm=T)
   dcyc12sd$cflxst      =qapply(X=dcycmean$cflxst      ,INDEX=yfac,DIM=1,FUN=sd  ,na.rm=T)
   dcyc12sd$hflxca      =qapply(X=dcycmean$hflxca      ,INDEX=yfac,DIM=1,FUN=sd  ,na.rm=T)
   dcyc12sd$hflxlc      =qapply(X=dcycmean$hflxlc      ,INDEX=yfac,DIM=1,FUN=sd  ,na.rm=T)
   dcyc12sd$hflxwc      =qapply(X=dcycmean$hflxwc      ,INDEX=yfac,DIM=1,FUN=sd  ,na.rm=T)
   dcyc12sd$hflxgc      =qapply(X=dcycmean$hflxgc      ,INDEX=yfac,DIM=1,FUN=sd  ,na.rm=T)
   dcyc12sd$wflxca      =qapply(X=dcycmean$wflxca      ,INDEX=yfac,DIM=1,FUN=sd  ,na.rm=T)
   dcyc12sd$qwflxca     =qapply(X=dcycmean$qwflxca     ,INDEX=yfac,DIM=1,FUN=sd  ,na.rm=T)
   dcyc12sd$wflxlc      =qapply(X=dcycmean$wflxlc      ,INDEX=yfac,DIM=1,FUN=sd  ,na.rm=T)
   dcyc12sd$wflxwc      =qapply(X=dcycmean$wflxwc      ,INDEX=yfac,DIM=1,FUN=sd  ,na.rm=T)
   dcyc12sd$wflxgc      =qapply(X=dcycmean$wflxgc      ,INDEX=yfac,DIM=1,FUN=sd  ,na.rm=T)
   dcyc12sd$evap        =qapply(X=dcycmean$evap        ,INDEX=yfac,DIM=1,FUN=sd  ,na.rm=T)
   dcyc12sd$transp      =qapply(X=dcycmean$transp      ,INDEX=yfac,DIM=1,FUN=sd  ,na.rm=T)
   dcyc12sd$atm.temp    =qapply(X=dcycmean$atm.temp    ,INDEX=yfac,DIM=1,FUN=sd  ,na.rm=T)
   dcyc12sd$can.temp    =qapply(X=dcycmean$can.temp    ,INDEX=yfac,DIM=1,FUN=sd  ,na.rm=T)
   dcyc12sd$leaf.temp   =qapply(X=dcycmean$leaf.temp   ,INDEX=yfac,DIM=1,FUN=sd  ,na.rm=T)
   dcyc12sd$wood.temp   =qapply(X=dcycmean$wood.temp   ,INDEX=yfac,DIM=1,FUN=sd  ,na.rm=T)
   dcyc12sd$gnd.temp    =qapply(X=dcycmean$gnd.temp    ,INDEX=yfac,DIM=1,FUN=sd  ,na.rm=T)
   dcyc12sd$atm.shv     =qapply(X=dcycmean$atm.shv     ,INDEX=yfac,DIM=1,FUN=sd  ,na.rm=T)
   dcyc12sd$can.shv     =qapply(X=dcycmean$can.shv     ,INDEX=yfac,DIM=1,FUN=sd  ,na.rm=T)
   dcyc12sd$gnd.shv     =qapply(X=dcycmean$gnd.shv     ,INDEX=yfac,DIM=1,FUN=sd  ,na.rm=T)
   dcyc12sd$atm.vpd     =qapply(X=dcycmean$atm.vpd     ,INDEX=yfac,DIM=1,FUN=sd  ,na.rm=T)
   dcyc12sd$can.vpd     =qapply(X=dcycmean$can.vpd     ,INDEX=yfac,DIM=1,FUN=sd  ,na.rm=T)
   dcyc12sd$leaf.vpd    =qapply(X=dcycmean$leaf.vpd    ,INDEX=yfac,DIM=1,FUN=sd  ,na.rm=T)
   dcyc12sd$atm.co2     =qapply(X=dcycmean$atm.co2     ,INDEX=yfac,DIM=1,FUN=sd  ,na.rm=T)
   dcyc12sd$can.co2     =qapply(X=dcycmean$can.co2     ,INDEX=yfac,DIM=1,FUN=sd  ,na.rm=T)
   dcyc12sd$atm.prss    =qapply(X=dcycmean$atm.prss    ,INDEX=yfac,DIM=1,FUN=sd  ,na.rm=T)
   dcyc12sd$can.prss    =qapply(X=dcycmean$can.prss    ,INDEX=yfac,DIM=1,FUN=sd  ,na.rm=T)
   dcyc12sd$atm.vels    =qapply(X=dcycmean$atm.vels    ,INDEX=yfac,DIM=1,FUN=sd  ,na.rm=T)
   dcyc12sd$ustar       =qapply(X=dcycmean$ustar       ,INDEX=yfac,DIM=1,FUN=sd  ,na.rm=T)
   dcyc12sd$fs.open     =qapply(X=dcycmean$fs.open     ,INDEX=yfac,DIM=1,FUN=sd  ,na.rm=T)
   dcyc12sd$rain        =qapply(X=dcycmean$rain        ,INDEX=yfac,DIM=1,FUN=sd  ,na.rm=T)
   dcyc12sd$rshort      =qapply(X=dcycmean$rshort      ,INDEX=yfac,DIM=1,FUN=sd  ,na.rm=T)
   dcyc12sd$rshort.beam =qapply(X=dcycmean$rshort.beam ,INDEX=yfac,DIM=1,FUN=sd  ,na.rm=T)
   dcyc12sd$rshort.diff =qapply(X=dcycmean$rshort.diff ,INDEX=yfac,DIM=1,FUN=sd  ,na.rm=T)
   dcyc12sd$rshort.gnd  =qapply(X=dcycmean$rshort.gnd  ,INDEX=yfac,DIM=1,FUN=sd  ,na.rm=T)
   dcyc12sd$rshortup    =qapply(X=dcycmean$rshortup    ,INDEX=yfac,DIM=1,FUN=sd  ,na.rm=T)
   dcyc12sd$rlong       =qapply(X=dcycmean$rlong       ,INDEX=yfac,DIM=1,FUN=sd  ,na.rm=T)
   dcyc12sd$rlong.gnd   =qapply(X=dcycmean$rlong.gnd   ,INDEX=yfac,DIM=1,FUN=sd  ,na.rm=T)
   dcyc12sd$rlongup     =qapply(X=dcycmean$rlongup     ,INDEX=yfac,DIM=1,FUN=sd  ,na.rm=T)
   # dcyc12sd$par.tot     =qapply(X=dcycmean$par.tot     ,INDEX=yfac,DIM=1,FUN=sd  ,na.rm=T)
   # dcyc12sd$par.beam    =qapply(X=dcycmean$par.beam    ,INDEX=yfac,DIM=1,FUN=sd  ,na.rm=T)
   # dcyc12sd$par.diff    =qapply(X=dcycmean$par.diff    ,INDEX=yfac,DIM=1,FUN=sd  ,na.rm=T)
   # dcyc12sd$par.gnd     =qapply(X=dcycmean$par.gnd     ,INDEX=yfac,DIM=1,FUN=sd  ,na.rm=T)
   dcyc12sd$parup       =qapply(X=dcycmean$parup       ,INDEX=yfac,DIM=1,FUN=sd  ,na.rm=T)
   dcyc12sd$rnet        =qapply(X=dcycmean$rnet        ,INDEX=yfac,DIM=1,FUN=sd  ,na.rm=T)
   dcyc12sd$albedo      =qapply(X=dcycmean$albedo      ,INDEX=yfac,DIM=1,FUN=sd  ,na.rm=T)
   dcyc12sd$albedo.beam =qapply(X=dcycmean$albedo.beam ,INDEX=yfac,DIM=1,FUN=sd  ,na.rm=T)
   dcyc12sd$albedo.diff =qapply(X=dcycmean$albedo.diff ,INDEX=yfac,DIM=1,FUN=sd  ,na.rm=T)
   dcyc12sd$rlong.albedo=qapply(X=dcycmean$rlong.albedo,INDEX=yfac,DIM=1,FUN=sd  ,na.rm=T)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Remove all elements of the DBH/PFT class that do not have a single valid cohort   #
   # at any given time.                                                                    #
   #---------------------------------------------------------------------------------------#
   empty = datum$nplantpftdbh == 0
   datum$agbpftdbh       [empty] = NA
   datum$basareapftdbh   [empty] = NA
   datum$laipftdbh       [empty] = NA
   datum$waipftdbh       [empty] = NA
   datum$taipftdbh       [empty] = NA
   datum$gpppftdbh       [empty] = NA
   datum$npppftdbh       [empty] = NA
   datum$mcopftdbh       [empty] = NA
   datum$cbapftdbh       [empty] = NA
   datum$cbalightpftdbh  [empty] = NA
   datum$cbamoistpftdbh  [empty] = NA
   datum$cbarelpftdbh    [empty] = NA
   datum$ldroppftdbh     [empty] = NA
   datum$fsopftdbh       [empty] = NA
   datum$demandpftdbh    [empty] = NA
   datum$supplypftdbh    [empty] = NA
   datum$mortpftdbh      [empty] = NA
   datum$agemortpftdbh   [empty] = NA
   datum$ncbmortpftdbh   [empty] = NA
   datum$tfallmortpftdbh [empty] = NA
   datum$coldmortpftdbh  [empty] = NA
   datum$distmortpftdbh  [empty] = NA
   datum$growthpftdbh    [empty] = NA
   datum$nplantpftdbh    [empty] = NA
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Convert mortality and recruitment so it is scaled between 0 and 100%.             #
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
         cat("      +",description,"time series for all PFTs...","\n")

         #----- Load variable -------------------------------------------------------------#
         if (vnam %in% names(datum)){
            thisvar = datum[[vnam]]
            if (plog){
               #----- Eliminate non-positive values in case it is a log plot. -------------#
               thisvar[thisvar <= 0] = NA
            }#end if
         }else{
            thisvar = matrix(NA,ncol=npft+1,nrow=nyears)
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
            ylimit = pretty.xylim(u = thisvar[,selpft],fracexp=scalleg,is.log=plog)
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



            letitre = paste(description,lieu,sep=" - ")
            cols    = pft$colour[selpft]
            legs    = pft$name  [selpft]
            plot(x=datum$toyear,y=thisvar[,1],type="n",main=letitre,ylim=ylimit
                ,xlab="Year",ylab=unit,cex.main=0.7,log=xylog)

            if (drought.mark){
               for (n in 1:ndrought){
                  rect(xleft  = drought[[n]][1],ybottom = ydrought[1]
                      ,xright = drought[[n]][2],ytop    = ydrought[2]
                      ,col    = "gray84",border=NA)
               }#end for
            }#end if
            
            if (plotgrid){ 
               abline(v=axTicks(side=1),h=axTicks(side=2),col="gray52",lty="solid")
            }#end if
            for (n in 1:(npft+1)){
               if (selpft[n]){
                  lines(datum$toyear,thisvar[,n],type="l",col=pft$colour[n],lwd=lwidth)
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
      vnam        = thistspftdbh$vnam
      description = thistspftdbh$desc
      unit        = thistspftdbh$unit
      plog        = thistspftdbh$plog
      plotit      = thistspftdbh$plt
      
      #----- Load variable ----------------------------------------------------------------#
      if (vnam %in% names(datum)){
         thisvar = datum[[vnam]]
         if (plog){
            xylog="y"
            thisvar[thisvar <= 0] = NA
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
         outdir = paste(outpref,"tsdbh",sep="/")
         if (! file.exists(outdir)) dir.create(outdir)
         outvar = paste(outdir,vnam,sep="/")
         if (! file.exists(outvar)) dir.create(outvar)
         #---------------------------------------------------------------------------------#

         cat("      +",description,"time series for DBH class...","\n")


         #---------------------------------------------------------------------------------#
         #     Find the limit, make some room for the legend, and in case the field is a   #
         # constant, nudge the limits so the plot command will not complain.               #
         #---------------------------------------------------------------------------------#
         ylimit = pretty.xylim(u=thisvar[,,pftuse],fracexp=scalleg,is.log=plog)
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
               plot(x=datum$toyear,y=thisvar[,1,p],type="n",main=letitre,ylim=ylimit
                   ,xlab="Time",ylab=unit,cex.main=0.7,log=xylog)
               if (drought.mark){
                  for (n in 1:ndrought){
                     rect(xleft  = drought[[n]][1],ybottom = ydrought[1]
                         ,xright = drought[[n]][2],ytop    = ydrought[2]
                         ,col    = "gray84",border=NA)
                  }#end for
               }#end if
               if (plotgrid){ 
                  abline(v=axTicks(side=1),h=axTicks(side=2),col="gray52",lty="solid")
               }#end if
               for (d in seq(from=1,to=ndbh+1,by=1)){
                  lines(datum$toyear,thisvar[,d,p],type="l",col=dbhcols[d],lwd=lwidth)
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
   cat("    + Year-by-year comparisons of monthly means...","\n")
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

      plotit       = plotit && vname %in% ls() && vname %in% names(mont12mn)

      if (plotit){
         #---------------------------------------------------------------------------------#
         #    Copy the observations to a scratch variable.                                 #
         #---------------------------------------------------------------------------------#
         thisvar     = datum[[vname]]
         thismean    = mont12mn[[vname]]
         thissdev    = mont12sd[[vname]]
         if (length(mont12sd[[vname]]) == 0){
            thissdev = 0. * thismean
         }else{
            thissdev = mont12sd[[vname]]
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
         outdir   = paste(outpref,"ycomp",sep="/")
         outvar   = paste(outdir,vname,sep="/")
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
         ylimit = pretty.xylim(u=ylimit,fracexp=scalleg,is.log=FALSE)
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Loop over all years, and make one plot per year.                            #
         #---------------------------------------------------------------------------------#
         for (y in 1:nyears){
            #----- Retrieve the year and the variable for this year. ----------------------#
            year.now = datum$toyear[y]
            cyear    = sprintf("%4.4i",year.now)
            var.year = thisvar[yfac == year.now]
            #------------------------------------------------------------------------------#


            #----- Loop over formats. -----------------------------------------------------#
            for (o in 1:nout){
               fichier = paste(outvar,"/",vname,"-",cyear,".",outform[o],sep="")
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
               letitre = paste(description," - ",lieu,"\n","Monthly mean - ",cyear,sep="")
               plot(x=montmont,y=var.year,type="n",main=letitre,xlab="Time"
                   ,ylim=ylimit,ylab=paste("[",unit,"]",sep=""),log=plog,xaxt="n"
                   ,cex.main=cex.main)
               axis(side=1,at=montplot$levels,labels=montplot$labels,padj=montplot$padj)
               if (plotgrid){ 
                  abline(v=montplot$levels,h=axTicks(side=2),col="gray52",lty="solid")
               }#end if
               if (plotsd){
                  polygon(x=mod.x.poly,y=mod.y.poly,col=errcolours[2],angle=angle[2]
                         ,density=dens[1],lty="solid",lwd=shwd[1])
               }#end if
               points(x=montmont,y=var.year,col=lcolours[1],lwd=llwd[1],type=ltype
                     ,pch=16,cex=1.0)
               points(x=montmont,y=thismean,col=lcolours[2],lwd=llwd[2],type=ltype
                     ,pch=16,cex=1.0)
               if (plotsd){
                  legend( x       = legpos
                        , inset   = 0.01
                        , legend  = c(cyear,paste("Mean: ",yeara,"-",yearz,sep=""))
                        , fill    = errcolours
                        , angle   = angle
                        , density = dens
                        , lwd     = llwd
                        , col     = lcolours
                        , bg      = "white"
                        , title   = "Shaded areas = 1 SD"
                        , cex     = 1.0
                        , pch     = 16
                        )#end legend
               }else{
                  legend( x      = legpos
                        , inset  = 0.05
                        , legend = c(cyear,paste("Mean: ",yeara,"-",yearz,sep=""))
                        , col    = lcolours
                        , lwd    = llwd
                        , cex    = 1.0
                        , pch    = 16
                        )#end legend
               }#end if
               if (outform[o] == "x11"){
                  locator(n=1)
                  dev.off()
               }else{
                  dev.off()
               }#end if
               #---------------------------------------------------------------------------#
            } #end for outform
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
         cat("      +",description,"time series for all LUs...","\n")



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
            ylimit = pretty.xylim(u=thisvar[,sellu],fracexp=scalleg,is.log=plog)
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

            letitre = paste(description,lieu,sep=" - ")
            cols    = lucols[sellu]
            legs    = lunames[sellu]
            plot(datum$toyear,thisvar[,1],type="n",main=letitre,ylim=ylimit
                ,xlab="Year",ylab=unit,cex.main=0.7)

            if (drought.mark){
               for (n in 1:ndrought){
                  rect(xleft  = drought[[n]][1],ybottom = ydrought[1]
                      ,xright = drought[[n]][2],ytop    = ydrought[2]
                      ,col    = "gray84",border=NA)
               }#end for
            }#end if
            if (plotgrid){ 
               abline(v=axTicks(side=1),h=axTicks(side=2),col="gray52",lty="solid")
            }#end if
            for (n in 1:(nlu+1)){
               if (sellu[n]){
                  lines(datum$toyear,thisvar[,n],type="l",col=lucols[n],lwd=lwidth)
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
                  ylimit = range(c(ylimit,datum$dist[,ilu,jlu]),na.rm=TRUE)
               }#end if
            }#end for
         }#end for
         ylimit   = pretty.xylim(u=ylimit,fracexp=scalleg,is.log=FALSE)
         ydrought = c( ylimit[1] - 0.5 * diff(ylimit), ylimit[2] + 0.5 * diff(ylimit) )
         #---------------------------------------------------------------------------------#

         letitre = paste("Disturbance rates",lieu,sep=" - ")
         cols    = NULL
         legs    = NULL
         plot(datum$toyear,datum$dist[,1,1],type="n",main=letitre,ylim=ylimit
             ,xlab="Year",ylab="[1/yr]",cex.main=0.7)
            if (drought.mark){
               for (n in 1:ndrought){
                  rect(xleft  = drought[[n]][1],ybottom = ydrought[1]
                      ,xright = drought[[n]][2],ytop    = ydrought[2]
                      ,col    = "gray84",border=NA)
               }#end for
            }#end if
            if (plotgrid){ 
               abline(v=axTicks(side=1),h=axTicks(side=2),col="gray52",lty="solid")
            }#end if
         n = 0
         for (jlu in 1:nlu){
            for (ilu in 1:nlu){
               n = n + 1
               if (seldist[ilu,jlu]){
                  cols = c(cols,distcols[n])
                  legs = c(legs,distnames[n])
                  lines(datum$toyear,datum$dist[,ilu,jlu],type="l"
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
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #   Plot the time series diagrams showing annual means.                                 #
   #---------------------------------------------------------------------------------------#
   cat("      * Plot some time series...","\n")
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
      plotit       = tsernow$plt & all(vnames %in% names(year12mn))

      if (plotit){

         #---------------------------------------------------------------------------------#
         #    Check whether the time series directory exists.  If not, create it.          #
         #---------------------------------------------------------------------------------#
         outdir = paste(outpref,"tseries",sep="/")
         if (! file.exists(outdir)) dir.create(outdir)
         cat("      +",theme,"time series for several variables...","\n")


         #----- Define the number of layers. ----------------------------------------------#
         nlayers   = length(vnames)
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Find the limit, make some room for the legend, and in case the field is a   #
         # constant, nudge the limits so the plot command will not complain.               #
         #---------------------------------------------------------------------------------#
         ylimit    = NULL
         for (l in 1:nlayers){
            thisvar = year12mn[[vnames[l]]]
            ylimit  = range(c(ylimit,thisvar),na.rm=TRUE)
         }#end for
         ylimit = pretty.xylim(u=ylimit,fracexp=scalleg,is.log=plog)
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
            thisvar = year12mn[[vnames[1]]]

            letitre = paste(theme," - ",lieu," \n"," Time series: ",theme,sep="")

            plot(x=datum$toyear,y=thisvar,type="n",main=letitre,xlab="Year"
                ,ylim=ylimit,ylab=paste("[",unit,"]",sep=""),log=xylog,cex.main=cex.main)
            if (drought.mark){
               for (n in 1:ndrought){
                  rect(xleft  = drought[[n]][1],ybottom = ydrought[1]
                      ,xright = drought[[n]][2],ytop    = ydrought[2]
                      ,col    = "gray84",border=NA)
               }#end for
            }#end if
            if (plotgrid){ 
               abline(v=axTicks(side=1),h=axTicks(side=2),col="gray52",lty="solid")
            }#end if
            for (l in 1:nlayers){
               thisvar = year12mn[[vnames[l]]]
               points(x=datum$toyear,y=thisvar,col=lcolours[l],lwd=llwd[l],type=ltype
                     ,pch=16,cex=0.8)
            }#end for
            legend(x=legpos,inset=inset,legend=description,col=lcolours,lwd=llwd,cex=0.8)
            if (outform[o] == "x11"){
               locator(n=1)
               dev.off()
            }else{
               dev.off()
            }#end if
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
         cat("      +",theme,"diurnal cycle for several variables...","\n")


         #----- Define the number of layers. ----------------------------------------------#
         nlayers   = length(vnames)
         ylimit    = NULL
         for (l in 1:nlayers){
            thisvar = dcyc12mn[[vnames[l]]]
            ylimit  = c(ylimit,thisvar)
         }#end for
         ylimit = pretty.xylim(u=ylimit,fracexp=scalleg,is.log=length(grep("y",plog)) > 0)
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #      Loop over all months.                                                      #
         #---------------------------------------------------------------------------------#
         for (pmon in 1:12){
            cmon    = sprintf("%2.2i",pmon)
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
               legend(x=legpos,inset=0.05,legend=description,col=lcolours,lwd=llwd)
               if (outform[o] == "x11"){
                  locator(n=1)
                  dev.off()
               }else{
                  dev.off()
               }#end if
               #---------------------------------------------------------------------------#
            }#end for outform
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
         monaxis  = sort(unique(datum$year))
         soilaxis = slz
         nmon     = length(monaxis)
         nsoil    = nzg

         #----- Convert the vector data into an array. ------------------------------------#
         vararr  = year12mn[[vnam]]

         #----- Copy the first and the last year to make the edges buffered. --------------#
         first    = vararr[1,]
         first    = c(first,first[nzg],first[nzg])

         last     = vararr[nyears,]
         last     = c(last[1],last[1],last)
         #----------------------------------------------------------------------------------#



         #----- Bind first and last year to the array, to make the edges buffered. ---------#
         varbuff  = cbind(vararr[,1],vararr,vararr[,nzg])
         varbuff  = rbind(last,varbuff,first)
         #----------------------------------------------------------------------------------#



         #----------------------------------------------------------------------------------#
         #   Expand the month and year axes.  Make the -------------------------------------------#
         yearaxis = c(yeara-1,datum$toyear,yearz+1)
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
            sombreado(x=yearaxis,y=soilaxis,z=varbuff,levels=vlevels,nlevels=vnlev
                     ,color.palette=get(vcscheme)
                     ,plot.title=title(main=letitre,xlab="Month",ylab="Soil depth [m]"
                                      ,cex.main=0.7)
                     ,key.title=title(main=unit,cex.main=0.8)
                     ,key.log=pnlog
                     ,plot.axes={axis(side=1)
                                 axis(side=2,at=zat,labels=znice)
                                 if (hovgrid){
                                    abline(h=zat,v=axTicks(1),col="gray52",lty="dotted")
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
   #      Bar plot by DBH class.                                                           #
   #---------------------------------------------------------------------------------------#
   cat("    + Bar plot by DBH classes...","\n")
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
         thisvnam                  = datum[[vnam]]
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
         for (y in 1:nyears){

            #----- Find which year we are plotting. ---------------------------------------#
            cyear     = sprintf("%4.4i",datum$toyear[y])
            yy        = as.numeric(cyear)
            #------------------------------------------------------------------------------#


            #----- Loop over output formats. ----------------------------------------------#
            for (o in 1:nout){
               #------ Open the plot. -----------------------------------------------------#
               fichier = paste(outdir,"/",vnam,"-",cyear,"-",suffix,".",outform[o],sep="")
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
               letitre = paste(lieu,"\n",description," - Year : ",cyear,sep="")
               lexlab  = "DBH Classes"
               leylab  = paste(description," [",unit,"]",sep="")
               #---------------------------------------------------------------------------#


               #----- Plot all monthly means together. ------------------------------------#
               barplot(height=t(thisvnam[y,,]),names.arg=dbhnames[1:ndbh],width=1.0
                      ,main=letitre,xlab=lexlab,ylab=leylab,ylim=ylimit,legend.text=FALSE
                      ,beside=(! stacked),col=pftcol.use,log=xylog
                      ,border="gray23",xpd=FALSE,cex.main=cex.main)
               if (plotgrid & (! stacked)){
                  xgrid=0.5+(1:ndbh)*(1+npftuse)
                  abline(v=xgrid,col="gray46",lty="solid")
               }#end if
               box()
               legend(x="topleft",inset=0.01,legend=pftname.use,fill=pftcol.use
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
}#end for places
#q("no")
