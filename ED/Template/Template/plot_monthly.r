#----- Here is the user-defined variable section. -----------------------------------------#
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
outform        = thisoutform           # Formats for output file.  Supported formats are:
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

sasfixlimits   = FALSE          # Should I use a fixed scale for size and age-structure
                                # plots? (FALSE will set a suitable scale for each year)

ncolsfc        = 200            # Target number of colours for filled contour plots.
fcgrid         = TRUE           # Should I include a grid on the filled contour plots?

ncolshov       = 200            # Target number of colours for Hovmoller diagrams.
hovgrid        = TRUE           # Should I include a grid on the Hovmoller plots?

legwhere       = "topleft"      # Where should I place the legend?
inset          = 0.01           # inset distance between legend and edge of plot region.
legbg          = "white"        # Legend background colour.
scalleg        = 0.20
cex.main       = 0.8             # Scale coefficient for the title
slz.min        = -5.0           # Find the deepest depth that trees access water.

theta           = 315.                    # Azimuth for perspective projection
phi             = 30.                     # Vertical angle for perspective projection
ltheta          = -210.                   # Azimuth angle for light
shade           = 0.125                   # Shade intensity
expz            = 0.5                     # Expansion factor for Z axis
gcol            = c("lightblue","white")  # Colours for the 50's style floor
cexmin          = 0.5                     # Minimum "head" size of the lollipop
cexmax          = 3.0                     # Maximum "head" size of the lollipop

ylnudge         = 0.05                    # Nudging factor for ylimit
ptype          = "l"                      # Type of plot
ptyped         = "p"                      # Type of plot
ptypeb         = "o"                      # Type of plot

tserdist        = TRUE                    # Time series of disturbance rates

drought.mark    = mydroughtmark
drought.yeara   = mydroughtyeara
drought.yearz   = mydroughtyearz
months.drought  = mymonthsdrought

#----- Loading some packages. -------------------------------------------------------------#
library(hdf5)
library(chron)
library(scatterplot3d)
library(lattice)
library(maps)
library(mapdata)
library(akima)
library(Hmisc)

#----- In case there is some graphic still opened. ----------------------------------------#
graphics.off()

#----- Set how many formats we must output. -----------------------------------------------#
outform = tolower(outform)
nout = length(outform)

#----- Avoid unecessary and extremely annoying beeps. -------------------------------------#
options(locatorBell=FALSE)


#----- Load some files with functions. ----------------------------------------------------#
source(paste(srcdir,"atlas.r"           ,sep="/"))
source(paste(srcdir,"charutils.r"       ,sep="/"))
source(paste(srcdir,"census.r"          ,sep="/"))
source(paste(srcdir,"cloudy.r"          ,sep="/"))
source(paste(srcdir,"error.bar.r"       ,sep="/"))
source(paste(srcdir,"globdims.r"        ,sep="/"))
source(paste(srcdir,"locations.r"       ,sep="/"))
source(paste(srcdir,"muitas.r"          ,sep="/"))
source(paste(srcdir,"plotsize.r"        ,sep="/"))
source(paste(srcdir,"pmonthly_varlist.r",sep="/"))
source(paste(srcdir,"pretty.log.r"      ,sep="/"))
source(paste(srcdir,"pretty.time.r"     ,sep="/"))
source(paste(srcdir,"qapply.r"          ,sep="/"))
source(paste(srcdir,"rconstants.r"      ,sep="/"))
source(paste(srcdir,"soilutils.r"       ,sep="/"))
source(paste(srcdir,"sombreado.r"       ,sep="/"))
source(paste(srcdir,"southammap.r"      ,sep="/"))
source(paste(srcdir,"thermlib.r"        ,sep="/"))
source(paste(srcdir,"timeutils.r"       ,sep="/"))
#----- These should be called after the others. --------------------------------------------#
source(paste(srcdir,"pft.coms.r"       ,sep="/"))


#----- Load observations. -----------------------------------------------------------------#
obsrfile = paste(srcdir,"LBA_MIP.v7.RData",sep="/")
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
   outpref = paste(outmain,"monthly",sep="/")
   lieu    = thispoi$lieu
   iata    = thispoi$iata
   suffix  = thispoi$iata
   yeara   = thispoi$yeara
   yearz   = thispoi$yearz
   meszz   = thispoi$monz

   #----- Decide how frequently the cohort-level variables should be saved. ---------------#
   if (yearz - yeara + 1 <= nyears.long){
      sasmonth = sasmonth.short
   }else{
      sasmonth = sasmonth.long
   }#end if



   #----- Create the directories in case they don't exist. --------------------------------#
   if (! file.exists(outmain)) dir.create(outmain)
   if (! file.exists(outpref)) dir.create(outpref)

   #----- Print a banner to entretain the user. -------------------------------------------#
   print(paste(" + Post-processing output from ",lieu,"...",sep=""))


   #---------------------------------------------------------------------------------------#
   #     Flush all variables that will hold the data.                                      #
   #---------------------------------------------------------------------------------------#
   totmon      = (yearz-yeara-1)*12+meszz+(12-monthbeg+1)
   #----- Size (DBH) and age arrays. ------------------------------------------------------#
   agbpftdbh       = array(data=0.,dim=c(totmon,ndbh+1,npft))
   laipftdbh       = array(data=0.,dim=c(totmon,ndbh+1,npft))
   waipftdbh       = array(data=0.,dim=c(totmon,ndbh+1,npft))
   taipftdbh       = array(data=0.,dim=c(totmon,ndbh+1,npft))
   gpppftdbh       = array(data=0.,dim=c(totmon,ndbh+1,npft))
   npppftdbh       = array(data=0.,dim=c(totmon,ndbh+1,npft))
   mcopftdbh       = array(data=0.,dim=c(totmon,ndbh+1,npft))
   cbapftdbh       = array(data=0.,dim=c(totmon,ndbh+1,npft))
   cbalightpftdbh  = array(data=0.,dim=c(totmon,ndbh+1,npft))
   cbamoistpftdbh  = array(data=0.,dim=c(totmon,ndbh+1,npft))
   cbarelpftdbh    = array(data=0.,dim=c(totmon,ndbh+1,npft))
   ldrpftdbh       = array(data=0.,dim=c(totmon,ndbh+1,npft))
   fsopftdbh       = array(data=0.,dim=c(totmon,ndbh+1,npft))
   demandpftdbh    = array(data=0.,dim=c(totmon,ndbh+1,npft))
   supplypftdbh    = array(data=0.,dim=c(totmon,ndbh+1,npft))
   nplantpftdbh    = array(data=0.,dim=c(totmon,ndbh+1,npft))
   mortpftdbh      = array(data=0.,dim=c(totmon,ndbh+1,npft))
   agemortpftdbh   = array(data=0.,dim=c(totmon,ndbh+1,npft))
   ncbmortpftdbh   = array(data=0.,dim=c(totmon,ndbh+1,npft))
   tfallmortpftdbh = array(data=0.,dim=c(totmon,ndbh+1,npft))
   coldmortpftdbh  = array(data=0.,dim=c(totmon,ndbh+1,npft))
   distmortpftdbh  = array(data=0.,dim=c(totmon,ndbh+1,npft))
   growthpftdbh    = array(data=0.,dim=c(totmon,ndbh+1,npft))
   #----- PFT arrays.   The "+1" column contains the total. -------------------------------#
   agbpft          = matrix(data=0,nrow=totmon,ncol=npft+1)
   bseedspft       = matrix(data=0,nrow=totmon,ncol=npft+1)
   nplantpft       = matrix(data=0,nrow=totmon,ncol=npft+1)
   laipft          = matrix(data=0,nrow=totmon,ncol=npft+1)
   waipft          = matrix(data=0,nrow=totmon,ncol=npft+1)
   taipft          = matrix(data=0,nrow=totmon,ncol=npft+1)
   gpppft          = matrix(data=0,nrow=totmon,ncol=npft+1)
   npppft          = matrix(data=0,nrow=totmon,ncol=npft+1)
   mcopft          = matrix(data=0,nrow=totmon,ncol=npft+1)
   cbapft          = matrix(data=0,nrow=totmon,ncol=npft+1)
   cbalightpft     = matrix(data=0,nrow=totmon,ncol=npft+1)
   cbamoistpft     = matrix(data=0,nrow=totmon,ncol=npft+1)
   cbarelpft       = matrix(data=0,nrow=totmon,ncol=npft+1)
   ldroppft        = matrix(data=0,nrow=totmon,ncol=npft+1)
   fsopft          = matrix(data=0,nrow=totmon,ncol=npft+1)
   balivepft       = matrix(data=0,nrow=totmon,ncol=npft+1)
   bdeadpft        = matrix(data=0,nrow=totmon,ncol=npft+1)
   bleafpft        = matrix(data=0,nrow=totmon,ncol=npft+1)
   brootpft        = matrix(data=0,nrow=totmon,ncol=npft+1)
   bswoodpft       = matrix(data=0,nrow=totmon,ncol=npft+1)
   bstorepft       = matrix(data=0,nrow=totmon,ncol=npft+1)
   basareapft      = matrix(data=0,nrow=totmon,ncol=npft+1)
   leafresppft     = matrix(data=0,nrow=totmon,ncol=npft+1)
   rootresppft     = matrix(data=0,nrow=totmon,ncol=npft+1)
   growthresppft   = matrix(data=0,nrow=totmon,ncol=npft+1)
   mortpft         = matrix(data=0,nrow=totmon,ncol=npft+1)
   agemortpft      = matrix(data=0,nrow=totmon,ncol=npft+1)
   ncbmortpft      = matrix(data=0,nrow=totmon,ncol=npft+1)
   tfallmortpft    = matrix(data=0,nrow=totmon,ncol=npft+1)
   coldmortpft     = matrix(data=0,nrow=totmon,ncol=npft+1)
   distmortpft     = matrix(data=0,nrow=totmon,ncol=npft+1)
   recrpft         = matrix(data=0,nrow=totmon,ncol=npft+1)
   growthpft       = matrix(data=0,nrow=totmon,ncol=npft+1)
   censusnplantpft = matrix(data=0,nrow=totmon,ncol=npft+1)
   censuslaipft    = matrix(data=0,nrow=totmon,ncol=npft+1)
   censuswaipft    = matrix(data=0,nrow=totmon,ncol=npft+1)
   censustaipft    = matrix(data=0,nrow=totmon,ncol=npft+1)
   censusagbpft    = matrix(data=0,nrow=totmon,ncol=npft+1)
   censusbapft     = matrix(data=0,nrow=totmon,ncol=npft+1)

   #----- LU arrays.   The "+1" column contains the total. --------------------------------#
   agblu           = matrix(data=0,nrow=totmon,ncol=nlu+1)
   lailu           = matrix(data=0,nrow=totmon,ncol=nlu+1)
   gpplu           = matrix(data=0,nrow=totmon,ncol=nlu+1)
   npplu           = matrix(data=0,nrow=totmon,ncol=nlu+1)
   arealu          = matrix(data=0,nrow=totmon,ncol=nlu+1)
   basarealu       = matrix(data=0,nrow=totmon,ncol=nlu+1)
   #----- Miscellaneous arrays. -----------------------------------------------------------#
   dist           = array(NA,dim=c(totmon,nlu,nlu))
   #----- Polygon level vectors. ----------------------------------------------------------#
   fast.soil.c     = NULL
   slow.soil.c     = NULL
   struct.soil.c   = NULL
   gpp             = NULL
   npp             = NULL
   plresp          = NULL
   leaf.resp       = NULL
   root.resp       = NULL
   growth.resp     = NULL
   hetresp         = NULL
   reco            = NULL
   mco             = NULL
   npp             = NULL
   cba             = NULL
   cbalight        = NULL
   cbamoist        = NULL
   cbarel          = NULL
   ldrop           = NULL
   nep             = NULL
   nee             = NULL
   cflxca          = NULL
   cflxst          = NULL
   evap            = NULL
   transp          = NULL
   ustar           = NULL
   atm.vels        = NULL
   atm.prss        = NULL
   atm.temp        = NULL
   can.prss        = NULL
   can.temp        = NULL
   atm.co2         = NULL
   can.co2         = NULL
   leaf.temp       = NULL
   wood.temp       = NULL
   atm.shv         = NULL
   can.shv         = NULL
   atm.vpd         = NULL
   can.vpd         = NULL
   leaf.vpd        = NULL
   can.co2         = NULL
   hflxca          = NULL
   qwflxca         = NULL
   wflxca          = NULL
   agb             = NULL
   nplant          = NULL
   lai             = NULL
   wai             = NULL
   tai             = NULL
   area            = NULL
   rain            = NULL
   gnd.temp        = NULL
   gnd.shv         = NULL
   workload        = NULL
   specwork        = NULL
   fs.open         = NULL
   hflxgc          = NULL
   hflxlc          = NULL
   hflxwc          = NULL
   wflxgc          = NULL
   wflxlc          = NULL
   wflxwc          = NULL
   et              = NULL
   rshort          = NULL
   rshort.beam     = NULL
   rshort.diff     = NULL
   rlong           = NULL
   rshort.gnd      = NULL
   rlong.gnd       = NULL
   rlongup         = NULL
   albedo          = NULL
   albedo.beam     = NULL
   albedo.diff     = NULL
   rlong.albedo    = NULL
   paw             = NULL
   smpot           = NULL
   npat.global     = NULL
   ncoh.global     = NULL
   mmsqu.gpp       = NULL
   mmsqu.plresp    = NULL
   mmsqu.leaf.resp = NULL
   mmsqu.root.resp = NULL
   mmsqu.plresp    = NULL
   mmsqu.hetresp   = NULL
   mmsqu.reco      = NULL
   mmsqu.cflxca    = NULL
   mmsqu.cflxst    = NULL
   mmsqu.hflxca    = NULL
   mmsqu.hflxlc    = NULL
   mmsqu.hflxwc    = NULL
   mmsqu.hflxgc    = NULL
   mmsqu.wflxca    = NULL
   mmsqu.qwflxca   = NULL
   mmsqu.wflxlc    = NULL
   mmsqu.wflxwc    = NULL
   mmsqu.wflxgc    = NULL
   mmsqu.evap      = NULL
   mmsqu.transp    = NULL

   #----- Cohort level lists. -------------------------------------------------------------#
   lightco      = list()
   lightbeamco  = list()
   lightdiffco  = list()
   parlco       = list()
   gppco        = list()
   gpplco       = list()
   respco       = list()
   nppco        = list()
   cbaco        = list()
   cbarelco     = list()
   cbalightco   = list()
   cbamoistco   = list()
   mcostco      = list()
   ldropco      = list()
   agbco        = list()
   fsoco        = list()
   nplantco     = list()
   pftco        = list()
   dbhco        = list()
   laico        = list()
   waico        = list()
   taico        = list()
   ageco        = list()
   areaco       = list()
   demandco     = list()
   supplyco     = list()
   heightco     = list()
   baco         = list()
   baliveco     = list()
   bdeadco      = list()
   bleafco      = list()
   brootco      = list()
   bswoodco     = list()
   bstoreco     = list()
   growthco     = list()
   mortco       = list()
   agemortco    = list()
   ncbmortco    = list()
   tfallmortco  = list()
   coldmortco   = list()
   distmortco   = list()
   recruitco    = list()

   n            = 0
   m            = 0
   thismonth    = NULL
   monnum       = NULL
   myear        = NULL

   first.time   = TRUE


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




   #----- Loop over years. ----------------------------------------------------------------#
   for (year in yeara:yearz){
       if (year == yeara){
          firstmonth = monthbeg
       }else{
          firstmonth = 1
       }#end if
       if (year == yearz){
          lastmonth = meszz
       }else{
          lastmonth = 12
       }#end if
       print (paste("    - Reading data from year ",year,"...",sep=""))

       #----- Loop over months. -----------------------------------------------------------#
       for (month in firstmonth:lastmonth){
          m = m + 1
          
          #----- Find out the previous month. ---------------------------------------------#
          lastmonth = 12 * (month == 1) + month - 1

          #----- Build the month and year vector. -----------------------------------------#
          monnum = c(monnum,month)
          myear  = c(myear,year)

          #----- Build the file name. -----------------------------------------------------#
          cmonth = substring(100+month,2,3)
          ddd    = daymax(month,year)
          myfile = paste(inpref,"-Q-",year,"-",cmonth,"-00-000000-g01.h5",sep="")
          #--------------------------------------------------------------------------------#



          #----- Read data and close connection immediately after. ------------------------#
          mymont = hdf5load(file=myfile,load=FALSE,verbosity=0,tidy=TRUE)
          #--------------------------------------------------------------------------------#


          #----- Build the time. ----------------------------------------------------------#
          thismonth = c(thismonth,chron(paste(month,1,year,sep="/"),times="0:0:0"))
          #--------------------------------------------------------------------------------#


          #----- Define the number of soil layers. ----------------------------------------#
          nzg      = mymont$NZG
          nzs      = mymont$NZS
          ndcycle  = mymont$NDCYCLE
          #--------------------------------------------------------------------------------#


          #--------------------------------------------------------------------------------#
          #      If this is the first time, allocate and initialise the mean diurnal cycle #
          # arrays.                                                                        #
          #--------------------------------------------------------------------------------#
          if (first.time){
             first.time          = FALSE

             #----- Find which soil are we solving, and save properties into soil.prop. ---#
             isoilflg   = mymont$ISOILFLG
             slz        = mymont$SLZ
             slxsand    = mymont$SLXSAND
             slxclay    = mymont$SLXCLAY
             ntext      = mymont$NTEXT.SOIL[nzg]

             soil.prop  = soil.params(ntext,isoilflg,slxsand,slxclay)
             dslz       = diff(c(slz,0))
             soil.depth = rev(cumsum(rev(dslz)))
             soil.dry   = rev(cumsum(rev(soil.prop$soilcp * wdns * dslz)))
             soil.poro  = rev(cumsum(rev(soil.prop$slmsts * wdns * dslz)))


             #----- Find the layers we care about. ----------------------------------------#
             sel        = slz < slz.min
             if (any(sel)){
                ka      = which.max(slz[sel])
             }else{
                ka      = 1
             }#end if
             kz         = nzg
             #-----------------------------------------------------------------------------#



             #----- Mean diurnal cycle. ---------------------------------------------------#
             dcycmean                = list()
             dcycmean$gpp            = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$npp            = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$plresp         = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$leaf.resp      = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$root.resp      = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$growth.resp    = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$hetresp        = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$nep            = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$nee            = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$reco           = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$cflxca         = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$cflxst         = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$hflxca         = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$hflxlc         = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$hflxwc         = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$hflxgc         = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$qwflxca        = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$wflxca         = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$wflxlc         = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$wflxwc         = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$wflxgc         = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$evap           = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$transp         = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$atm.temp       = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$can.temp       = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$leaf.temp      = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$wood.temp      = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$gnd.temp       = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$atm.shv        = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$can.shv        = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$gnd.shv        = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$atm.vpd        = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$can.vpd        = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$leaf.vpd       = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$atm.co2        = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$can.co2        = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$atm.prss       = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$can.prss       = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$atm.vels       = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$ustar          = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$fs.open        = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$rain           = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$rshort         = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$rshort.beam    = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$rshort.diff    = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$rlong          = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$rshort.gnd     = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$rlong.gnd      = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$rlongup        = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$albedo         = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$albedo.beam    = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$albedo.diff    = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$rlong.albedo   = matrix(data=0,nrow=totmon,ncol=ndcycle)

             dcycmsqu             = list()
             dcycmsqu$gpp         = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmsqu$npp         = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmsqu$plresp      = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmsqu$leaf.resp   = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmsqu$root.resp   = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmsqu$hetresp     = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmsqu$nep         = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmsqu$nee         = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmsqu$reco        = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmsqu$cflxca      = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmsqu$cflxst      = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmsqu$hflxca      = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmsqu$hflxlc      = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmsqu$hflxwc      = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmsqu$hflxgc      = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmsqu$qwflxca     = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmsqu$wflxca      = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmsqu$wflxlc      = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmsqu$wflxwc      = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmsqu$wflxgc      = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmsqu$transp      = matrix(data=0,nrow=totmon,ncol=ndcycle)

             soil.water           = matrix(data=0,nrow=totmon,ncol=nzg)
             soil.temp            = matrix(data=0,nrow=totmon,ncol=nzg)
             soil.mstpot          = matrix(data=0,nrow=totmon,ncol=nzg)

          }#end if
          #--------------------------------------------------------------------------------#

          #----- Load the total number of patches and cohorts. ----------------------------#
          npat.global = c(npat.global, mymont$NPATCHES.GLOBAL)
          ncoh.global = c(ncoh.global, mymont$NCOHORTS.GLOBAL)
          #--------------------------------------------------------------------------------#



          #----- Load the simple variables. -----------------------------------------------#
          gpp             = c(gpp              ,   mymont$MMEAN.GPP                      )
          plresp          = c(plresp           ,   mymont$MMEAN.PLRESP                   )
          leaf.resp       = c(leaf.resp        ,   mymont$MMEAN.LEAF.RESP                )
          root.resp       = c(root.resp        ,   mymont$MMEAN.ROOT.RESP                )
          growth.resp     = c(growth.resp      ,   mymont$MMEAN.GROWTH.RESP              )
          hetresp         = c(hetresp          ,   mymont$MMEAN.RH                       )
          nep             = c(nep              ,   mymont$MMEAN.NEP                      )
          nee             = c(nee              , - mymont$MMEAN.CARBON.AC                
                                                 + mymont$MMEAN.CARBON.ST                )
          reco            = c(reco             ,   mymont$MMEAN.PLRESP                   
                                                 + mymont$MMEAN.RH                       )
          cflxca          = c(cflxca           , - mymont$MMEAN.CARBON.AC                )
          cflxst          = c(cflxst           ,   mymont$MMEAN.CARBON.ST                )
          hflxca          = c(hflxca           , - mymont$MMEAN.SENSIBLE.AC              )
          hflxlc          = c(hflxlc           ,   mymont$MMEAN.SENSIBLE.LC              )
          hflxwc          = c(hflxwc           ,   mymont$MMEAN.SENSIBLE.WC              )
          hflxgc          = c(hflxgc           ,   mymont$MMEAN.SENSIBLE.GC              )
          qwflxca         = c(qwflxca          , - mymont$MMEAN.VAPOR.AC    
                                                 * alvli(mymont$MMEAN.CAN.TEMP)          )
          wflxca          = c(wflxca           , - mymont$MMEAN.VAPOR.AC    * day.sec    )
          wflxlc          = c(wflxlc           ,   mymont$MMEAN.VAPOR.LC    * day.sec    )
          wflxwc          = c(wflxwc           ,   mymont$MMEAN.VAPOR.WC    * day.sec    )
          wflxgc          = c(wflxgc           ,   mymont$MMEAN.VAPOR.GC    * day.sec    )
          evap            = c(evap             ,   mymont$MMEAN.EVAP        * day.sec    )
          transp          = c(transp           ,   mymont$MMEAN.TRANSP      * day.sec    )

          mmsqu.gpp       = c(mmsqu.gpp       , mymont$MMSQU.GPP                         )
          mmsqu.plresp    = c(mmsqu.plresp    , mymont$MMSQU.PLRESP                      )
          mmsqu.leaf.resp = c(mmsqu.leaf.resp , mymont$MMSQU.PLRESP                      )
          mmsqu.root.resp = c(mmsqu.root.resp , mymont$MMSQU.PLRESP                      )
          mmsqu.hetresp   = c(mmsqu.hetresp   , mymont$MMSQU.RH                          )
          mmsqu.cflxca    = c(mmsqu.cflxca    , mymont$MMSQU.CARBON.AC                   )
          mmsqu.cflxst    = c(mmsqu.cflxst    , mymont$MMSQU.CARBON.ST                   )
          mmsqu.hflxca    = c(mmsqu.hflxca    , mymont$MMSQU.SENSIBLE.AC                 )
          mmsqu.hflxlc    = c(mmsqu.hflxlc    , mymont$MMSQU.SENSIBLE.LC                 )
          mmsqu.hflxwc    = c(mmsqu.hflxwc    , mymont$MMSQU.SENSIBLE.WC                 )
          mmsqu.hflxgc    = c(mmsqu.hflxgc    , mymont$MMSQU.SENSIBLE.GC                 )
          mmsqu.wflxca    = c(mmsqu.wflxca    , mymont$MMSQU.VAPOR.AC  
                                              * day.sec * day.sec                        )
          mmsqu.qwflxca   = c(mmsqu.qwflxca   , mymont$MMSQU.VAPOR.AC  
                                              * alvli(mymont$MMEAN.CAN.TEMP)
                                              * alvli(mymont$MMEAN.CAN.TEMP)              )
          mmsqu.wflxlc    = c(mmsqu.wflxlc    , mymont$MMSQU.VAPOR.LC  
                                              * day.sec * day.sec                        )
          mmsqu.wflxwc    = c(mmsqu.wflxwc    , mymont$MMSQU.VAPOR.WC
                                              * day.sec * day.sec                        )
          mmsqu.wflxgc    = c(mmsqu.wflxgc    , mymont$MMSQU.VAPOR.GC
                                              * day.sec * day.sec                        )
          mmsqu.evap      = c(mmsqu.evap      , mymont$MMSQU.EVAP
                                              * day.sec * day.sec                        )
          mmsqu.transp    = c(mmsqu.transp    , mymont$MMSQU.TRANSP
                                              * day.sec * day.sec                        )

          ustar         = c(ustar        ,mymont$MMEAN.USTAR                             )

          atm.vels      = c(atm.vels     ,mymont$MMEAN.ATM.VELS                          )
          atm.prss      = c(atm.prss     ,mymont$MMEAN.ATM.PRSS  * 0.01                  )
          atm.temp      = c(atm.temp     ,mymont$MMEAN.ATM.TEMP  - t00                   )
          atm.shv       = c(atm.shv      ,mymont$MMEAN.ATM.SHV   * kg2g                  )
          atm.co2       = c(atm.co2      ,mymont$MMEAN.ATM.CO2                           )
          atm.vpd       = c(atm.vpd      ,mymont$MMEAN.ATM.VPDEF * 0.01                  )

          can.prss      = c(can.prss     ,mymont$MMEAN.CAN.PRSS  * 0.01                  )
          can.temp      = c(can.temp     ,mymont$MMEAN.CAN.TEMP  - t00                   )
          can.shv       = c(can.shv      ,mymont$MMEAN.CAN.SHV   * kg2g                  )
          can.co2       = c(can.co2      ,mymont$MMEAN.CAN.CO2                           )
          can.vpd       = c(can.vpd      ,mymont$MMEAN.CAN.VPDEF * 0.01                  )

          gnd.temp      = c(gnd.temp     ,mymont$MMEAN.GND.TEMP  - t00                   )
          gnd.shv       = c(gnd.shv      ,mymont$MMEAN.GND.SHV   * kg2g                  )

          leaf.temp     = c(leaf.temp    ,mymont$MMEAN.LEAF.TEMP  - t00                  )
          leaf.vpd      = c(leaf.vpd     ,mymont$MMEAN.LEAF.VPDEF * 0.01                 )
          wood.temp     = c(wood.temp    ,mymont$MMEAN.WOOD.TEMP  - t00                  )
          rain          = c(rain         ,mymont$MMEAN.PCPG*ddd  * day.sec               )

          fs.open       = c(fs.open      ,mymont$MMEAN.FS.OPEN                           )
          rshort        = c(rshort       ,mymont$MMEAN.RSHORT                            )
          rshort.beam   = c(rshort.beam  ,mymont$MMEAN.RSHORT - mymont$MMEAN.RSHORT.DIFF )
          rshort.diff   = c(rshort.diff  ,mymont$MMEAN.RSHORT.DIFF                       )
          rlong         = c(rlong        ,mymont$MMEAN.RLONG                             )
          rshort.gnd    = c(rshort.gnd   ,mymont$MMEAN.RSHORT.GND                        )
          rlong.gnd     = c(rlong.gnd    ,mymont$MMEAN.RLONG.GND                         )
          rlongup       = c(rlongup      ,mymont$MMEAN.RLONGUP                           )
          albedo        = c(albedo       ,mymont$MMEAN.ALBEDO                            )
          albedo.beam   = c(albedo.beam  ,mymont$MMEAN.ALBEDO.BEAM                       )
          albedo.diff   = c(albedo.diff  ,mymont$MMEAN.ALBEDO.DIFFUSE                    )
          rlong.albedo  = c(rlong.albedo ,mymont$MMEAN.RLONG.ALBEDO                      )
          #--------------------------------------------------------------------------------#


          #----- Find averaged soil properties. -------------------------------------------#
          swater.now  = rev(cumsum(rev(mymont$MMEAN.SOIL.WATER * wdns * dslz)))
          smoist.avg  = swater.now / (wdns * soil.depth)
          paw.now     = 100. * ( ( swater.now[ka] - soil.dry [ka] )
                               / ( soil.poro [ka] - soil.dry [ka] ) )
          smpot.now   = ( - smoist2mpot(smoist=smoist.avg[ka],mysoil=soil.prop)
                        * 0.001 * grav )
          paw         = c(   paw,   paw.now )
          smpot       = c( smpot, smpot.now )
          #--------------------------------------------------------------------------------#


          #------ Read in soil properties. ------------------------------------------------#
          soil.temp  [m,] =   mymont$MMEAN.SOIL.TEMP - t00
          soil.water [m,] =   mymont$MMEAN.SOIL.WATER
          soil.mstpot[m,] = - mymont$MMEAN.SOIL.MSTPOT
          #--------------------------------------------------------------------------------#



          #----- Read workload, and retrieve only the current month. ----------------------#
          workload   = c(workload  , mymont$WORKLOAD[month])
          specwork   = c(specwork  , mymont$WORKLOAD[month]/sum(mymont$SIPA.N,na.rm=TRUE))
          #--------------------------------------------------------------------------------#



          #----- Disturbance rates. -------------------------------------------------------#
          dist  [m,1:nlu,1:nlu] = mymont$DISTURBANCE.RATES[,1:nlu,1:nlu]
          #--------------------------------------------------------------------------------#


          #--------------------------------------------------------------------------------#
          #       Read the mean diurnal cycle and the mean sum of the squares.             #
          #--------------------------------------------------------------------------------#
          dcycmean$gpp         [m,] = mymont$QMEAN.GPP
          dcycmean$plresp      [m,] = mymont$QMEAN.PLRESP
          dcycmean$npp         [m,] = mymont$QMEAN.GPP - mymont$QMEAN.PLRESP
          dcycmean$leaf.resp   [m,] = mymont$QMEAN.LEAF.RESP
          dcycmean$root.resp   [m,] = mymont$QMEAN.ROOT.RESP
          dcycmean$hetresp     [m,] = mymont$QMEAN.RH
          dcycmean$nep         [m,] = mymont$QMEAN.NEP
          dcycmean$nee         [m,] = - mymont$QMEAN.CARBON.AC + mymont$QMEAN.CARBON.ST
          dcycmean$reco        [m,] = mymont$QMEAN.PLRESP + mymont$QMEAN.RH
          dcycmean$cflxca      [m,] = - mymont$QMEAN.CARBON.AC
          dcycmean$cflxst      [m,] = - mymont$QMEAN.CARBON.ST
          dcycmean$hflxca      [m,] = - mymont$QMEAN.SENSIBLE.AC
          dcycmean$hflxlc      [m,] = mymont$QMEAN.SENSIBLE.LC
          dcycmean$hflxwc      [m,] = mymont$QMEAN.SENSIBLE.WC
          dcycmean$hflxgc      [m,] = mymont$QMEAN.SENSIBLE.GC
          dcycmean$wflxca      [m,] = - mymont$QMEAN.VAPOR.AC         * day.sec
          dcycmean$qwflxca     [m,] = ( - mymont$QMEAN.VAPOR.AC
                                        * alvli(mymont$QMEAN.CAN.TEMP) )
          dcycmean$wflxlc      [m,] = mymont$QMEAN.VAPOR.LC           * day.sec
          dcycmean$wflxwc      [m,] = mymont$QMEAN.VAPOR.WC           * day.sec
          dcycmean$wflxgc      [m,] = mymont$QMEAN.VAPOR.GC           * day.sec
          dcycmean$evap        [m,] = ( mymont$QMEAN.VAPOR.GC
                                      + mymont$QMEAN.VAPOR.WC
                                      + mymont$QMEAN.VAPOR.LC )       * day.sec
          dcycmean$transp      [m,] = mymont$QMEAN.TRANSP             * day.sec
          dcycmean$atm.temp    [m,] = mymont$QMEAN.ATM.TEMP           - t00
          dcycmean$can.temp    [m,] = mymont$QMEAN.CAN.TEMP           - t00
          dcycmean$leaf.temp   [m,] = mymont$QMEAN.LEAF.TEMP          - t00
          dcycmean$wood.temp   [m,] = mymont$QMEAN.WOOD.TEMP          - t00
          dcycmean$gnd.temp    [m,] = mymont$QMEAN.GND.TEMP           - t00
          dcycmean$atm.shv     [m,] = mymont$QMEAN.ATM.SHV            * kg2g
          dcycmean$can.shv     [m,] = mymont$QMEAN.CAN.SHV            * kg2g
          dcycmean$gnd.shv     [m,] = mymont$QMEAN.GND.SHV            * kg2g
          dcycmean$atm.vpd     [m,] = mymont$QMEAN.ATM.VPDEF          * 0.01
          dcycmean$can.vpd     [m,] = mymont$QMEAN.CAN.VPDEF          * 0.01
          dcycmean$leaf.vpd    [m,] = mymont$QMEAN.LEAF.VPDEF         * 0.01
          dcycmean$atm.co2     [m,] = mymont$QMEAN.ATM.CO2
          dcycmean$can.co2     [m,] = mymont$QMEAN.CAN.CO2
          dcycmean$atm.vels    [m,] = mymont$QMEAN.ATM.VELS
          dcycmean$ustar       [m,] = mymont$QMEAN.USTAR
          dcycmean$atm.prss    [m,] = mymont$QMEAN.ATM.PRSS * 0.01
          dcycmean$can.prss    [m,] = mymont$QMEAN.CAN.PRSS * 0.01
          dcycmean$fs.open     [m,] = mymont$QMEAN.FS.OPEN
          dcycmean$rain        [m,] = mymont$QMEAN.PCPG               * day.sec
          dcycmean$rshort      [m,] = mymont$QMEAN.RSHORT
          dcycmean$rshort.beam [m,] = mymont$QMEAN.RSHORT - mymont$QMEAN.RSHORT.DIFF
          dcycmean$rshort.diff [m,] = mymont$QMEAN.RSHORT.DIFF
          dcycmean$rlong       [m,] = mymont$QMEAN.RLONG
          dcycmean$rshort.gnd  [m,] = mymont$QMEAN.RSHORT.GND
          dcycmean$rlong.gnd   [m,] = mymont$QMEAN.RLONG.GND
          dcycmean$rlongup     [m,] = mymont$QMEAN.RLONGUP
          dcycmean$albedo      [m,] = mymont$QMEAN.ALBEDO
          dcycmean$albedo.beam [m,] = mymont$QMEAN.ALBEDO.BEAM
          dcycmean$albedo.diff [m,] = mymont$QMEAN.ALBEDO.DIFFUSE
          dcycmean$rlong.albedo[m,] = mymont$QMEAN.RLONG.ALBEDO

          dcycmsqu$gpp         [m,] = mymont$QMSQU.GPP
          dcycmsqu$plresp      [m,] = mymont$QMSQU.PLRESP
          dcycmsqu$leaf.resp   [m,] = mymont$QMSQU.LEAF.RESP
          dcycmsqu$root.resp   [m,] = mymont$QMSQU.ROOT.RESP
          dcycmsqu$hetresp     [m,] = mymont$QMSQU.RH
          dcycmsqu$nep         [m,] = mymont$QMSQU.NEP
          dcycmsqu$cflxca      [m,] = mymont$QMSQU.CARBON.AC
          dcycmsqu$cflxst      [m,] = mymont$QMSQU.CARBON.ST
          dcycmsqu$hflxca      [m,] = mymont$QMSQU.SENSIBLE.AC
          dcycmsqu$hflxlc      [m,] = mymont$QMSQU.SENSIBLE.LC
          dcycmsqu$hflxwc      [m,] = mymont$QMSQU.SENSIBLE.WC
          dcycmsqu$hflxgc      [m,] = mymont$QMSQU.SENSIBLE.GC
          dcycmsqu$wflxca      [m,] = mymont$QMSQU.VAPOR.AC    * day.sec * day.sec
          dcycmsqu$qwflxca     [m,] = ( mymont$QMSQU.VAPOR.AC    
                                      * alvli(mymont$QMEAN.CAN.TEMP)
                                      * alvli(mymont$QMEAN.CAN.TEMP) )
          dcycmsqu$wflxlc      [m,] = mymont$QMSQU.VAPOR.WC    * day.sec * day.sec
          dcycmsqu$wflxwc      [m,] = mymont$QMSQU.VAPOR.LC    * day.sec * day.sec
          dcycmsqu$wflxgc      [m,] = mymont$QMSQU.VAPOR.GC    * day.sec * day.sec
          dcycmsqu$transp      [m,] = mymont$QMSQU.TRANSP      * day.sec * day.sec
          #--------------------------------------------------------------------------------#


          #---- Read in the site-level area. ----------------------------------------------#
          areasi     = mymont$AREA.SI
          npatches   = mymont$SIPA.N
          #--------------------------------------------------------------------------------#


          #----- Read a few patch-level variables. ----------------------------------------#
          areapa     = mymont$AREA * rep(areasi,times=npatches)
          lupa       = mymont$DIST.TYPE
          agepa      = mymont$AGE
          #--------------------------------------------------------------------------------#


          #--------------------------------------------------------------------------------#
          #     Get the soil carbon.                                                       #
          #--------------------------------------------------------------------------------#
          fast.soil.c   = c(fast.soil.c  ,sum(mymont$FAST.SOIL.C       * areapa))
          slow.soil.c   = c(slow.soil.c  ,sum(mymont$SLOW.SOIL.C       * areapa))
          struct.soil.c = c(struct.soil.c,sum(mymont$STRUCTURAL.SOIL.C * areapa))
          #--------------------------------------------------------------------------------#


          #--------------------------------------------------------------------------------#
          #    If this is a biomass initialisation, or a run with anthropogenic            #
          # disturbance, we must jitter the age so we can distinguish the patches.         #
          #--------------------------------------------------------------------------------#
          sameage        = duplicated(agepa)
          agepa[sameage] = jitter(x=agepa[sameage],amount=0.4)
          #--------------------------------------------------------------------------------#


          #--------------------------------------------------------------------------------#
          #     Read the cohort-level variables.  Because empty patchs do exist (deserts), #
          # we must check whether there is any cohort to be read.  If not, assign NA to    #
          # all variables.                                                                 #
          #--------------------------------------------------------------------------------#
          ncohorts   = mymont$PACO.N
          if (any (ncohorts > 0)){
             areaconow  = rep(areapa,times=ncohorts)

             #----- Define the land use classes. ------------------------------------------#
             luconow    = rep(lupa,times=ncohorts)

             #----- Define the DBH classes. -----------------------------------------------#
             dbhconow   = mymont$DBH
             dbhbks     = c(-Inf,seq(from=ddbh,to=(ndbh-1)*ddbh,by=ddbh),Inf)
             dbhcut     = cut(dbhconow,breaks=dbhbks)
             dbhlevs    = levels(dbhcut)
             dbhfac     = match(dbhcut,dbhlevs)

             #----- Define the age classes. -----------------------------------------------#
             ageconow   = rep(x=agepa,times=ncohorts)
             agebks     = c(-Inf,seq(from=dage,to=(nage-1)*dage,by=dage),Inf)
             agecut     = cut(ageconow,breaks=agebks)
             agelevs    = levels(agecut)
             agefac     = match(agecut,agelevs)

             agepacut   = cut(agepa,breaks=agebks)
             agepafac   = match(agepacut,agelevs)
             areapaage  = tapply(X=areapa,INDEX=agepafac,sum,na.rm=TRUE)
             areaage    = areapaage[as.character(agefac)]

             pftconow        = mymont$PFT
             nplantconow     = mymont$NPLANT
             heightconow     = mymont$HITE
             baconow         = mymont$BA.CO
             agbconow        = mymont$AGB.CO
             bseedsconow     = mymont$BSEEDS.CO
             laiconow        = mymont$LAI.CO
             waiconow        = mymont$WAI.CO
             taiconow        = laiconow + waiconow
             gppconow        = mymont$MMEAN.GPP.CO
             leafrespconow   = mymont$MMEAN.LEAF.RESP.CO
             rootrespconow   = mymont$MMEAN.ROOT.RESP.CO
             growthrespconow = mymont$MMEAN.GROWTH.RESP.CO
             respconow       = ( mymont$MMEAN.LEAF.RESP.CO   + mymont$MMEAN.ROOT.RESP.CO
                               + mymont$MMEAN.GROWTH.RESP.CO + mymont$MMEAN.STORAGE.RESP.CO
                               + mymont$MMEAN.VLEAF.RESP.CO  )
             nppconow        = gppconow-respconow
             cbaconow        = mymont$MMEAN.CB
             cbalightconow   = mymont$CB.LIGHTMAX[,lastmonth]
             cbamoistconow   = mymont$CB.MOISTMAX[,lastmonth]
             cbarelconow     = mymont$CBR.BAR
             mcostconow      = ( mymont$MMEAN.LEAF.MAINTENANCE
                               + mymont$MMEAN.ROOT.MAINTENANCE )
             ldropconow      = mymont$MMEAN.LEAF.DROP.CO
             fsoconow        = mymont$MMEAN.FS.OPEN.CO
             lightconow      = mymont$MMEAN.LIGHT.LEVEL
             lightbeamconow  = mymont$MMEAN.LIGHT.LEVEL.BEAM
             lightdiffconow  = mymont$MMEAN.LIGHT.LEVEL.DIFF
             parlconow       = mymont$MMEAN.PAR.L

             baliveconow     = mymont$BALIVE
             bdeadconow      = mymont$BDEAD
             bleafconow      = mymont$BLEAF
             brootconow      = mymont$BROOT
             bswoodconow     = mymont$BSAPWOODA+mymont$BSAPWOODB
             bstoreconow     = mymont$BSTORAGE


             sel               = laiconow > 1.e-10
             demandconow       = mymont$MMEAN.PSI.OPEN     * laiconow * day.sec
             supplyconow       = mymont$MMEAN.WATER.SUPPLY * day.sec
             gpplconow         = gppconow
             gpplconow  [sel]  = nplantconow[sel] * gppconow[sel] / laiconow[sel]
             gpplconow  [!sel] = 0.


             mortconow       = rowSums(mymont$MMEAN.MORT.RATE)
             agemortconow    = mymont$MMEAN.MORT.RATE[,1]
             ncbmortconow    = mymont$MMEAN.MORT.RATE[,2]
             tfallmortconow  = mymont$MMEAN.MORT.RATE[,3]
             coldmortconow   = mymont$MMEAN.MORT.RATE[,4]
             distmortconow   = mymont$MMEAN.MORT.RATE[,5]
             recruitconow    = mymont$RECRUIT.DBH
             growthconow     = 100. * mymont$DLNDBH.DT
             #-----------------------------------------------------------------------------#
          }else{
             #----- Make everything NA. ---------------------------------------------------#
             areaconow       = NA
             luconow         = NA
             dbhconow        = NA
             dbhbks          = NA
             dbhcut          = NA
             dbhlevs         = NA
             dbhfac          = NA
             ageconow        = NA
             agebks          = NA
             agecut          = NA
             agelevs         = NA
             agefac          = NA
             agepacut        = NA
             agepafac        = NA
             areapaage       = NA
             areaage         = NA
             pftconow        = NA
             nplantconow     = NA
             heightconow     = NA
             agbconow        = NA
             baconow         = NA
             bseedsconow     = NA
             laiconow        = NA
             waiconow        = NA
             taiconow        = NA
             gppconow        = NA
             gpplconow       = NA
             leafrespconow   = NA
             rootrespconow   = NA
             growthrespconow = NA
             respconow       = NA
             nppconow        = NA
             cbaconow        = NA
             cbalightconow   = NA
             cbamoistconow   = NA
             mcostconow      = NA
             ldropconow      = NA
             cbarelconow     = NA
             fsoconow        = NA
             lightconow      = NA
             lightbeamconow  = NA
             lightdiffconow  = NA
             parlconow       = NA
             demandconow     = NA
             supplyconow     = NA
             baliveconow     = NA
             bdeadconow      = NA
             bleafconow      = NA
             brootconow      = NA
             bswoodconow     = NA
             bstoreconow     = NA
             mortconow       = NA
             agemortconow    = NA
             ncbmortconow    = NA
             tfallmortconow  = NA
             coldmortconow   = NA
             distmortconow   = NA
             recruitconow    = NA
             growthconow     = NA
          }#end if
          #--------------------------------------------------------------------------------#


          #--------------------------------------------------------------------------------#
          #     The following two variables are used to scale "intensive" properties       #
          # (whatever/plant) to "extensive" (whatever/m2).  Sometimes they may be used to  #
          # build weighted averages.                                                       #
          #--------------------------------------------------------------------------------#
          w.nplant = nplantconow * areaconow
          w.lai    = laiconow    * areaconow
          #--------------------------------------------------------------------------------#


          #----- Define some classes that can be defined even when no cohorts exist. ------#
          agebks     = c(-Inf,seq(from=dage,to=(nage-1)*dage,by=dage),Inf)
          agepacut   = cut(agepa,breaks=agebks)
          agepafac   = match(agepacut,agelevs)
          areapaage  = tapply(X=areapa,INDEX=agepafac,sum,na.rm=TRUE)
          areaage    = areapaage[as.character(agefac)]
          #--------------------------------------------------------------------------------#



          #--------------------------------------------------------------------------------#
          #     Build the PFT arrays.                                                      #
          #--------------------------------------------------------------------------------#
          for (p in 1:npft){
             if (all(is.na(pftconow))){
                sel      = rep(FALSE,times=length(pftconow))
             }else{
                sel      = pftconow == p
             }#end if
             if (any(sel)){

                #----- "Extensive" variables, add them. -----------------------------------#
                nplantpft    [m,p] = sum(nplantconow[sel] * areaconow[sel])
                laipft       [m,p] = sum(laiconow   [sel] * areaconow[sel])
                waipft       [m,p] = sum(waiconow   [sel] * areaconow[sel])
                taipft       [m,p] = sum(taiconow   [sel] * areaconow[sel])
                #--------------------------------------------------------------------------#

                #----- "Intensive" variables, use nplant or LAI to make them extensive. ---#
                basareapft   [m,p] = sum( w.nplant[sel] * baconow        [sel] )
                agbpft       [m,p] = sum( w.nplant[sel] * agbconow       [sel] )
                bseedspft    [m,p] = sum( w.nplant[sel] * bseedsconow    [sel] )
                gpppft       [m,p] = sum( w.nplant[sel] * gppconow       [sel] )
                npppft       [m,p] = sum( w.nplant[sel] * nppconow       [sel] )
                mcopft       [m,p] = sum( w.nplant[sel] * mcostconow     [sel] )
                cbapft       [m,p] = sum( w.nplant[sel] * cbaconow       [sel] )
                cbalightpft  [m,p] = sum( w.nplant[sel] * cbalightconow  [sel] )
                cbamoistpft  [m,p] = sum( w.nplant[sel] * cbamoistconow  [sel] )
                ldroppft     [m,p] = sum( w.nplant[sel] * ldropconow     [sel] )
                balivepft    [m,p] = sum( w.nplant[sel] * baliveconow    [sel] )
                bdeadpft     [m,p] = sum( w.nplant[sel] * bdeadconow     [sel] )
                bleafpft     [m,p] = sum( w.nplant[sel] * bleafconow     [sel] )
                brootpft     [m,p] = sum( w.nplant[sel] * brootconow     [sel] )
                bswoodpft    [m,p] = sum( w.nplant[sel] * bswoodconow    [sel] )
                bstorepft    [m,p] = sum( w.nplant[sel] * bstoreconow    [sel] )
                leafresppft  [m,p] = sum( w.nplant[sel] * leafrespconow  [sel] )
                rootresppft  [m,p] = sum( w.nplant[sel] * rootrespconow  [sel] )
                growthresppft[m,p] = sum( w.nplant[sel] * growthrespconow[sel] )

                #----- Fso is added here, but it will be scaled back to intensive. --------#
                fsopft       [m,p] = sum( w.lai   [sel] * fsoconow       [sel] )
                #--------------------------------------------------------------------------#

                #----- cbarel is added here, but it will be scaled back to intensive. -----#
                cbarelpft    [m,p] = sum( w.nplant[sel] * cbarelconow    [sel] )
                #--------------------------------------------------------------------------#
             }#end if (any(sel))
          }#end for (p in 1:npft)
          #------ Find the total. ---------------------------------------------------------#
          nplantpft    [m,npft+1] = sum(nplantpft    [m,1:npft],na.rm=TRUE)
          laipft       [m,npft+1] = sum(laipft       [m,1:npft],na.rm=TRUE)
          waipft       [m,npft+1] = sum(waipft       [m,1:npft],na.rm=TRUE)
          taipft       [m,npft+1] = sum(taipft       [m,1:npft],na.rm=TRUE)
          agbpft       [m,npft+1] = sum(agbpft       [m,1:npft],na.rm=TRUE)
          bseedspft    [m,npft+1] = sum(bseedspft    [m,1:npft],na.rm=TRUE)
          gpppft       [m,npft+1] = sum(gpppft       [m,1:npft],na.rm=TRUE)
          npppft       [m,npft+1] = sum(npppft       [m,1:npft],na.rm=TRUE)
          mcopft       [m,npft+1] = sum(mcopft       [m,1:npft],na.rm=TRUE)
          cbapft       [m,npft+1] = sum(cbapft       [m,1:npft],na.rm=TRUE)
          ldroppft     [m,npft+1] = sum(ldroppft     [m,1:npft],na.rm=TRUE)
          balivepft    [m,npft+1] = sum(balivepft    [m,1:npft],na.rm=TRUE)
          bdeadpft     [m,npft+1] = sum(bdeadpft     [m,1:npft],na.rm=TRUE)
          bleafpft     [m,npft+1] = sum(bleafpft     [m,1:npft],na.rm=TRUE)
          brootpft     [m,npft+1] = sum(brootpft     [m,1:npft],na.rm=TRUE)
          bswoodpft    [m,npft+1] = sum(bswoodpft    [m,1:npft],na.rm=TRUE)
          bstorepft    [m,npft+1] = sum(bstorepft    [m,1:npft],na.rm=TRUE)
          basareapft   [m,npft+1] = sum(basareapft   [m,1:npft],na.rm=TRUE)
          leafresppft  [m,npft+1] = sum(leafresppft  [m,1:npft],na.rm=TRUE)
          rootresppft  [m,npft+1] = sum(rootresppft  [m,1:npft],na.rm=TRUE)
          growthresppft[m,npft+1] = sum(growthresppft[m,1:npft],na.rm=TRUE)
          #--------------------------------------------------------------------------------#


         #----- Find the average Fsopen for each PFT. -------------------------------------#
         for (p in 1:npft){
            if (laipft[m,p] != 0){
               fsopft[m,p] = fsopft[m,p] / laipft[m,p]
            }#end if
         }#end for
         if (laipft[m,npft+1] != 0){
            fsopft[m,npft+1] = sum(fsopft[m,1:npft] * laipft[m,1:npft]) / laipft[m,npft+1]
         }#end for
         #---------------------------------------------------------------------------------#


         #----- Find the average Fsopen for each PFT. -------------------------------------#
         for (p in 1:npft){
            if (nplantpft[m,p] != 0){
               cbarelpft[m,p] = cbarelpft[m,p] / nplantpft[m,p]
            }#end if
         }#end for
         if (nplantpft[m,npft+1] != 0){
            cbarelpft[m,npft+1] = ( sum(cbarelpft[m,1:npft] * nplantpft[m,1:npft]) 
                                  / nplantpft[m,npft+1] )
         }#end for
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Find the growth, mortality, and recruitment rates for each PFT, and the     #
         # global rates.  We only use the cohorts that were flagged as 1 or 2 (which means #
         # that their DBH is greater than 10 cm).                                          #
         #---------------------------------------------------------------------------------#
         for (p in 1:npft){
            if (all(is.na(pftconow))){
               sel.hgt      = rep(FALSE,times=length(pftconow))
               sel.dbh      = rep(FALSE,times=length(pftconow))
               sel.re2      = rep(FALSE,times=length(pftconow))
            }else{
               sel.hgt      = pftconow == p & heightconow  >  1.5
               sel.dbh      = pftconow == p & recruitconow >  0
               sel.re2      = pftconow == p & recruitconow == 2
            }#end if
            #------------------------------------------------------------------------------#


            #----- Census LAI, WAI, and TAI neglect plants shorter than 1.5m. -------------#
            if (any(sel.hgt)){
               censuslaipft    [m,p] = sum(laiconow[sel.hgt] * areaconow[sel.hgt])
               censuswaipft    [m,p] = sum(waiconow[sel.hgt] * areaconow[sel.hgt])
               censustaipft    [m,p] = sum(taiconow[sel.hgt] * areaconow[sel.hgt])
            }#end if
            #------------------------------------------------------------------------------#



            #----- Census AGB and Basal area ignore plants whose DBH is less than 10cm. ---#
            if (any(sel.dbh)){
               censusbapft [m,p] = sum( w.nplant[sel.dbh] * baconow  [sel.dbh] )
               censusagbpft[m,p] = sum( w.nplant[sel.dbh] * agbconow [sel.dbh] )
            }#end if
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #     Find the PFT-level mortality rates.                                      #
            #------------------------------------------------------------------------------#
            if (any(sel.re2)){
               #---- This is the number of survivors. -------------------------------------#
               survivor          = sum(w.nplant[sel.re2])
               previous          = sum(w.nplant[sel.re2] * exp(mortconow[sel.re2]) )
               mortpft[m,p]      = 100. * log( previous / survivor )

               survivor          = sum(w.nplant[sel.re2])
               previous          = sum(w.nplant[sel.re2] * exp(agemortconow[sel.re2]) )
               agemortpft[m,p]   = 100. * log( previous / survivor )

               survivor          = sum(w.nplant[sel.re2])
               previous          = sum(w.nplant[sel.re2] * exp(ncbmortconow[sel.re2]) )
               ncbmortpft[m,p]   = 100. * log( previous / survivor )

               survivor          = sum(w.nplant[sel.re2])
               previous          = sum(w.nplant[sel.re2] * exp(tfallmortconow[sel.re2]) )
               tfallmortpft[m,p] = 100. * log( previous / survivor )

               survivor          = sum(w.nplant[sel.re2])
               previous          = sum(w.nplant[sel.re2] * exp(coldmortconow[sel.re2]) )
               coldmortpft[m,p]  = 100. * log( previous / survivor )

               survivor          = sum(w.nplant[sel.re2])
               previous          = sum(w.nplant[sel.re2] * exp(distmortconow[sel.re2]) )
               distmortpft[m,p]  = 100. * log( previous / survivor )
            }#end if
            if (any(sel.dbh) & any(sel.re2)){
               #---- This is the number of survivors. -------------------------------------#
               population   = sum(w.nplant[sel.dbh])
               established  = sum(w.nplant[sel.re2])
               recrpft[m,p] = 100. * log( population / established) / 12.0
            }#end if
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Growth rates are found only for established cohorts.                     #
            #------------------------------------------------------------------------------#
            if (any(sel.re2)){
               censusnplantpft [m,p] = sum( w.nplant[sel.re2] )
               growthpft[m,p]        = sum( w.nplant[sel.re2] * growthconow [sel.re2] )
            }#end if
            #------------------------------------------------------------------------------#
         }#end for
         #---------------------------------------------------------------------------------#




         #---------------------------------------------------------------------------------#
         #      Find the total or the mean amongst all PFTs.                               #
         #---------------------------------------------------------------------------------#
         #----- "Extensive" properties, just add them up. ---------------------------------#
         censusnplantpft    [m,npft+1] = sum(censusnplantpft    [m,1:npft],na.rm=TRUE)
         censuslaipft       [m,npft+1] = sum(censuslaipft       [m,1:npft],na.rm=TRUE)
         censuswaipft       [m,npft+1] = sum(censuswaipft       [m,1:npft],na.rm=TRUE)
         censustaipft       [m,npft+1] = sum(censustaipft       [m,1:npft],na.rm=TRUE)
         censusagbpft       [m,npft+1] = sum(censusagbpft       [m,1:npft],na.rm=TRUE)
         censusbapft        [m,npft+1] = sum(censusbapft        [m,1:npft],na.rm=TRUE)
         #----- Growth rates are averaged, not added. -------------------------------------#
         if (censusnplantpft[m,npft+1] > 0.){
            growthpft[m,npft+1] = ( sum(growthpft[m,1:npft] * censusnplantpft[m,1:npft])
                                  / censusnplantpft[m,npft+1] )
         }#end if
         #---------------------------------------------------------------------------------#
         #      Mortality and recruitment rates are done by looking at the entire          #
         # population rather than averaged.                                                #
         #---------------------------------------------------------------------------------#
         if (all(is.na(pftconow))){
            sel.dbh = rep(FALSE,times=length(pftconow))
            sel.re2 = rep(FALSE,times=length(pftconow))
         }else{
            sel.dbh = recruitconow >  0
            sel.re2 = recruitconow == 2
         }#end if

         #----- w.nplant will convert intensive variables into extensive. -----------------#
         w.nplant          = nplantconow * areaconow
         #---------------------------------------------------------------------------------#

         if (any(sel.re2)){
            #---- This is the number of survivors. ----------------------------------------#
            survivor               = sum( w.nplant[sel.re2] )
            previous               = sum( w.nplant[sel.re2] * exp(mortconow[sel.re2]))
            mortpft[m,npft+1]      = 100. * log( previous / survivor )

            survivor               = sum( w.nplant[sel.re2] )
            previous               = sum( w.nplant[sel.re2] * exp(agemortconow[sel.re2]))
            agemortpft  [m,npft+1] = 100. * log( previous / survivor )

            survivor               = sum( w.nplant[sel.re2] )
            previous               = sum( w.nplant[sel.re2] * exp(ncbmortconow[sel.re2]))
            ncbmortpft  [m,npft+1] = 100. * log( previous / survivor )

            survivor               = sum( w.nplant[sel.re2] )
            previous               = sum( w.nplant[sel.re2] * exp(tfallmortconow[sel.re2]))
            tfallmortpft[m,npft+1] = 100. * log( previous / survivor )

            survivor               = sum( w.nplant[sel.re2] )
            previous               = sum( w.nplant[sel.re2] * exp(coldmortconow[sel.re2]))
            coldmortpft [m,npft+1] = 100. * log( previous / survivor )

            survivor               = sum( w.nplant[sel.re2] )
            previous               = sum( w.nplant[sel.re2] * exp(distmortconow[sel.re2]))
            distmortpft [m,npft+1] = 100. * log( previous / survivor )
         }#end if

         if (any(sel.dbh) & any(sel.re2)){
            #---- This is the number of survivors. ----------------------------------------#
            population        = sum( w.nplant[sel.dbh] )
            established       = sum( w.nplant[sel.re2] )
            recrpft[m,npft+1] = 100. * log( population / established ) / 12.0
         }#end if
         #---------------------------------------------------------------------------------#






         #---------------------------------------------------------------------------------#
         #     Build the LU arrays.                                                        #
         #---------------------------------------------------------------------------------#
         for (l in 1:nlu){
            selpa    = lupa    == l
            if (all(is.na(luconow))){
               sel      = rep(FALSE,times=length(luconow))
            }else{
               sel      = luconow == l
            }#end if
            
            if (any(sel)){
               lailu    [m,l] = sum( laiconow[sel] * areaconow[sel] )
               basarealu[m,l] = sum( w.nplant[sel] * baconow  [sel] )
               agblu [m,l]    = sum( w.nplant[sel] * agbconow [sel] )
               gpplu [m,l]    = sum( w.nplant[sel] * gppconow [sel] )
               npplu [m,l]    = sum( w.nplant[sel] * nppconow [sel] )
            }#end if
            arealu [m,l]   = arealu [m,l] + sum(areapa[selpa])
         }#end for
         #------ Find the total. ----------------------------------------------------------#
         lailu    [m,nlu+1] = sum(lailu    [m,1:nlu], na.rm=TRUE)
         basarealu[m,nlu+1] = sum(basarealu[m,1:nlu], na.rm=TRUE)
         agblu    [m,nlu+1] = sum(agblu    [m,1:nlu], na.rm=TRUE)
         gpplu    [m,nlu+1] = sum(gpplu    [m,1:nlu], na.rm=TRUE)
         npplu    [m,nlu+1] = sum(npplu    [m,1:nlu], na.rm=TRUE)
         arealu   [m,nlu+1] = sum(arealu   [m,1:nlu], na.rm=TRUE)
         #---------------------------------------------------------------------------------#




         #---------------------------------------------------------------------------------#
         #     Build the size (DBH) structure arrays.                                      #
         #---------------------------------------------------------------------------------#
         for (d in 1:ndbh){
            if (all(is.na(dbhfac))){
               seldbh  = rep(FALSE,times=length(dbhfac))
            }else{
               seldbh  = dbhfac == d
            }#end if
            for (p in 1:npft){
               selpft   = pftconow == p
               sel      = selpft & seldbh
               if (any(sel)){
                  #----- Extensive properties. --------------------------------------------#
                  laipftdbh      [m,d,p] = sum( laiconow    [sel] * areaconow     [sel] )
                  waipftdbh      [m,d,p] = sum( waiconow    [sel] * areaconow     [sel] )
                  taipftdbh      [m,d,p] = sum( taiconow    [sel] * areaconow     [sel] )
                  nplantpftdbh   [m,d,p] = sum( nplantconow [sel] * areaconow     [sel] )
                  demandpftdbh   [m,d,p] = sum( demandconow [sel] * areaconow     [sel] )
                  supplypftdbh   [m,d,p] = sum( supplyconow [sel] * areaconow     [sel] )
                  #----- Intensive properties, use nplant to make them extensive. ---------#
                  agbpftdbh      [m,d,p] = sum( w.nplant    [sel] * agbconow      [sel] )
                  gpppftdbh      [m,d,p] = sum( w.nplant    [sel] * gppconow      [sel] )
                  npppftdbh      [m,d,p] = sum( w.nplant    [sel] * nppconow      [sel] )
                  mcopftdbh      [m,d,p] = sum( w.nplant    [sel] * mcostconow    [sel] )
                  cbapftdbh      [m,d,p] = sum( w.nplant    [sel] * cbaconow      [sel] )
                  cbalightpftdbh [m,d,p] = sum( w.nplant    [sel] * cbalightconow [sel] )
                  cbamoistpftdbh [m,d,p] = sum( w.nplant    [sel] * cbamoistconow [sel] )
                  ldrpftdbh      [m,d,p] = sum( w.nplant    [sel] * ldropconow    [sel] )
                  #----- FSO, use LAI to scale them, we will normalise outside the loop. --#
                  fsopftdbh      [m,d,p] = sum( laiconow    [sel] * fsoconow      [sel] )
                  #----- CBAREL, like FSO, but we use nplant to scale them. ---------------#
                  cbarelpftdbh   [m,d,p] = sum( nplantconow [sel] * cbarelconow   [sel] )
                  #------------------------------------------------------------------------#
               }#end if

               #---------------------------------------------------------------------------#
               #    For mortality and growth, we keep deleting the tiny guys because they  #
               # skew the rates quite significantly.                                       #
               #---------------------------------------------------------------------------#
               sel = sel & recruitconow == 2
               if (any(sel)){
                  growthpftdbh [m,d,p] = sum( w.nplant [sel] * growthconow [sel] )

                  #---- This is the number of survivors and living before. ----------------#
                  survivor               = sum( w.nplant[sel] )
                  previous               = sum( w.nplant[sel] * exp(mortconow[sel]) )
                  mortpftdbh[m,d,p]      = 100. * log( previous / survivor )

                  survivor               = sum( w.nplant[sel] )
                  previous               = sum( w.nplant[sel] * exp(agemortconow[sel]) )
                  agemortpftdbh  [m,d,p] = 100. * log( previous / survivor )

                  survivor               = sum( w.nplant[sel] )
                  previous               = sum( w.nplant[sel] * exp(ncbmortconow[sel]) )
                  ncbmortpftdbh  [m,d,p] = 100. * log( previous / survivor )

                  survivor               = sum( w.nplant[sel] )
                  previous               = sum( w.nplant[sel] * exp(tfallmortconow[sel]) )
                  tfallmortpftdbh[m,d,p] = 100. * log( previous / survivor )

                  survivor               = sum( w.nplant[sel] )
                  previous               = sum( w.nplant[sel] * exp(coldmortconow[sel]) )
                  coldmortpftdbh [m,d,p] = 100. * log( previous / survivor )

                  survivor               = sum( w.nplant[sel] )
                  previous               = sum( w.nplant[sel] * exp(distmortconow[sel]) )
                  distmortpftdbh [m,d,p] = 100. * log( previous / survivor )
               }#end if
               #---------------------------------------------------------------------------#
            }#end for PFT
            #------------------------------------------------------------------------------#
         }#end for DBH
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Either sum or normalise variables.                                          #
         #---------------------------------------------------------------------------------#
         for (p in 1:npft){

            #---- Find the total for this PFT and store at class ndbh + 1... --------------#
            laipftdbh    [m,ndbh+1,p] = sum(laipftdbh    [m,1:ndbh,p])
            waipftdbh    [m,ndbh+1,p] = sum(waipftdbh    [m,1:ndbh,p])
            taipftdbh    [m,ndbh+1,p] = sum(taipftdbh    [m,1:ndbh,p])
            nplantpftdbh [m,ndbh+1,p] = sum(nplantpftdbh [m,1:ndbh,p])
            agbpftdbh    [m,ndbh+1,p] = sum(agbpftdbh    [m,1:ndbh,p])
            gpppftdbh    [m,ndbh+1,p] = sum(gpppftdbh    [m,1:ndbh,p])
            npppftdbh    [m,ndbh+1,p] = sum(npppftdbh    [m,1:ndbh,p])
            mcopftdbh    [m,ndbh+1,p] = sum(mcopftdbh    [m,1:ndbh,p])
            cbapftdbh    [m,ndbh+1,p] = sum(cbapftdbh    [m,1:ndbh,p])
            ldrpftdbh    [m,ndbh+1,p] = sum(ldrpftdbh    [m,1:ndbh,p])
            demandpftdbh [m,ndbh+1,p] = sum(demandpftdbh [m,1:ndbh,p])
            supplypftdbh [m,ndbh+1,p] = sum(supplypftdbh [m,1:ndbh,p])
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #     Find the average Fsopen for each DBH class and amongst all classes.      #
            #------------------------------------------------------------------------------#
            for (d in 1:ndbh){
               #---------------------------------------------------------------------------#
               #      FSopen is scaled by LAI.                                             #
               #---------------------------------------------------------------------------#
               if (laipftdbh[m,d,p] != 0){
                  fsopftdbh[m,d,p] = fsopftdbh[m,d,p] / laipftdbh[m,d,p]
               }#end if
            }#end for
            if (laipftdbh[m,ndbh+1,p] != 0){
               fsopftdbh[m,ndbh+1,p] = ( sum( fsopftdbh[m,1:ndbh,p]
                                            * laipftdbh[m,1:ndbh,p] )
                                       / laipftdbh[m,ndbh+1,p] )
            }#end for
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #     Find the average relative carbon balance for each DBH class and amongst  #
            # all classes.                                                                 #
            #------------------------------------------------------------------------------#
            for (d in 1:ndbh){
               #---------------------------------------------------------------------------#
               #      CBAREL is scaled by Nplant.                                          #
               #---------------------------------------------------------------------------#
               if (nplantpftdbh[m,d,p] != 0){
                  cbarelpftdbh[m,d,p] = cbarelpftdbh[m,d,p] / nplantpftdbh[m,d,p]
               }#end if
            }#end for
            if (nplantpftdbh[m,ndbh+1,p] != 0){
               cbarelpftdbh[m,ndbh+1,p] = ( sum( cbarelpftdbh[m,1:ndbh,p]
                                               * nplantpftdbh[m,1:ndbh,p] )
                                          / nplantpftdbh[m,ndbh+1,p] )
            }#end for
            #------------------------------------------------------------------------------#


            for (d in 2:ndbh){
               #---------------------------------------------------------------------------#
               #      Growth rate is scaled by NPlant.                                     #
               #---------------------------------------------------------------------------#
               nplant.use = 0.
               if (nplantpftdbh[m,d,p] != 0){
                  nplant.use          = nplant.use + nplantpftdbh[m,d,p]
                  growthpftdbh[m,d,p] = growthpftdbh[m,d,p] / nplantpftdbh[m,d,p]
               }#end if
            }#end for
            if (nplant.use != 0){
               growthpftdbh[m,ndbh+1,p] = ( sum( growthpftdbh[m,2:ndbh,p]
                                               * nplantpftdbh[m,2:ndbh,p] )
                                          / nplant.use )
            }#end for
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Find the mortality amongst all DBH classes for this PFT.  This will      #
            # include the small cohorts too.                                               #
            #------------------------------------------------------------------------------#
            if (all(is.na(pftconow))){
               sel = rep(FALSE,times=length(pftconow))
            }else{
               sel = pftconow == p & recruitconow == 2
            }#end if
            #---- This is the number of survivors and living before. ----------------------#
            survivor               = sum( w.nplant[sel] )
            previous               = sum( w.nplant[sel] * exp(mortconow[sel]) )
            mortpftdbh[m,ndbh+1,p] = 100. * log( previous / survivor)

            survivor                    = sum( w.nplant[sel] )
            previous                    = sum( w.nplant[sel] * exp(agemortconow[sel]) )
            agemortpftdbh[m,ndbh+1,p]   = 100. * log( previous / survivor)

            survivor                    = sum( w.nplant[sel] )
            previous                    = sum( w.nplant[sel] * exp(ncbmortconow[sel]) )
            ncbmortpftdbh[m,ndbh+1,p]   = 100. * log( previous / survivor)

            survivor                    = sum( w.nplant[sel] )
            previous                    = sum( w.nplant[sel] * exp(tfallmortconow[sel]) )
            tfallmortpftdbh[m,ndbh+1,p] = 100. * log( previous / survivor)

            survivor                    = sum( w.nplant[sel] )
            previous                    = sum( w.nplant[sel] * exp(coldmortconow[sel]) )
            coldmortpftdbh[m,ndbh+1,p]  = 100. * log( previous / survivor)

            survivor                    = sum( w.nplant[sel] )
            previous                    = sum( w.nplant[sel] * exp(distmortconow[sel]) )
            distmortpftdbh[m,ndbh+1,p]  = 100. * log( previous / survivor)
         }#end for
         #---------------------------------------------------------------------------------#




         #---------------------------------------------------------------------------------#
         #       Build the derived variables.                                              #
         #---------------------------------------------------------------------------------#
         npp    = c(npp   ,sum(npppft   [m,1:npft]) )
         mco    = c(mco   ,sum(mcopft   [m,1:npft]) )
         cba    = c(cba   ,sum(cbapft   [m,1:npft]) )
         nplant = c(nplant,sum(nplantpft[m,1:npft]) )
         lai    = c(lai   ,sum(laipft   [m,1:npft]) )
         wai    = c(wai   ,sum(waipft   [m,1:npft]) )
         tai    = c(tai   ,sum(taipft   [m,1:npft]) )
         agb    = c(agb   ,sum(agbpft   [m,1:npft]) )
         ldrop  = c(ldrop ,sum(ldroppft [m,1:npft]) )
         #---------------------------------------------------------------------------------#




         #---------------------------------------------------------------------------------#
         #      Build the cohort-level lists if this is the right month.                   #
         #---------------------------------------------------------------------------------#
         if (month %in% sasmonth){
            cyear  = substring(10000 + year,2,5)
            cmonth = substring(100+month,2,3)
            labwhen     = paste("y",cyear,"m",cmonth,sep="")
            #----- Binding the current cohorts. -------------------------------------------#
            lightco     [[labwhen]] = lightconow
            lightbeamco [[labwhen]] = lightbeamconow
            lightdiffco [[labwhen]] = lightdiffconow
            parlco      [[labwhen]] = parlconow
            gppco       [[labwhen]] = gppconow
            gpplco      [[labwhen]] = gpplconow
            respco      [[labwhen]] = respconow
            nppco       [[labwhen]] = nppconow
            cbaco       [[labwhen]] = cbaconow
            cbalightco  [[labwhen]] = cbalightconow
            cbamoistco  [[labwhen]] = cbamoistconow
            cbarelco    [[labwhen]] = cbarelconow
            mcostco     [[labwhen]] = mcostconow
            agbco       [[labwhen]] = agbconow
            fsoco       [[labwhen]] = fsoconow
            nplantco    [[labwhen]] = nplantconow * areaconow
            heightco    [[labwhen]] = heightconow
            baco        [[labwhen]] = nplantconow * baconow * areaconow
            pftco       [[labwhen]] = pftconow
            dbhco       [[labwhen]] = dbhconow
            laico       [[labwhen]] = laiconow
            waico       [[labwhen]] = waiconow
            taico       [[labwhen]] = taiconow
            ageco       [[labwhen]] = ageconow
            areaco      [[labwhen]] = areaconow
            demandco    [[labwhen]] = demandconow
            supplyco    [[labwhen]] = supplyconow
            baliveco    [[labwhen]] = baliveconow
            bdeadco     [[labwhen]] = bdeadconow
            bleafco     [[labwhen]] = bleafconow
            brootco     [[labwhen]] = brootconow
            bswoodco    [[labwhen]] = bswoodconow
            bstoreco    [[labwhen]] = bstoreconow
            growthco    [[labwhen]] = growthconow
            mortco      [[labwhen]] = mortconow
            agemortco   [[labwhen]] = agemortconow
            ncbmortco   [[labwhen]] = ncbmortconow
            tfallmortco [[labwhen]] = tfallmortconow
            coldmortco  [[labwhen]] = coldmortconow
            distmortco  [[labwhen]] = distmortconow
            recruitco   [[labwhen]] = recruitconow
         } #end if month=sasmonth
         #---------------------------------------------------------------------------------#
      }# end for, month
   }#end for, year



   #----- Find the date variables. --------------------------------------------------------#
   thismonth = chron(thismonth,out.format=c(dates="day-mon-yr",times=NULL))
   mmonth    = months(thismonth)
   mfac      = nummonths(thismonth)
   moncensus = census(pop=mfac,categ=seq(from=1,to=12,by=1))
   moncnt    = matrix(data=rep(moncensus,times=ndcycle),nrow=12)


   #---------------------------------------------------------------------------------------#
   #      Here we find the monthly means for month, then compute the standard deviation.   #
   #---------------------------------------------------------------------------------------#
   print ("    - Finding the monthly mean...")
   print ("      * Aggregating the monthly mean...")
   mont12mn             = list()
   mont12mn$gpp         = tapply(X=gpp          ,INDEX=mfac      ,FUN=mean,na.rm=TRUE)
   mont12mn$npp         = tapply(X=npp          ,INDEX=mfac      ,FUN=mean,na.rm=TRUE)
   mont12mn$nep         = tapply(X=nep          ,INDEX=mfac      ,FUN=mean,na.rm=TRUE)
   mont12mn$plresp      = tapply(X=plresp       ,INDEX=mfac      ,FUN=mean,na.rm=TRUE)
   mont12mn$leaf.resp   = tapply(X=leaf.resp    ,INDEX=mfac      ,FUN=mean,na.rm=TRUE)
   mont12mn$root.resp   = tapply(X=root.resp    ,INDEX=mfac      ,FUN=mean,na.rm=TRUE)
   mont12mn$growth.resp = tapply(X=growth.resp  ,INDEX=mfac      ,FUN=mean,na.rm=TRUE)
   mont12mn$hetresp     = tapply(X=hetresp      ,INDEX=mfac      ,FUN=mean,na.rm=TRUE)
   mont12mn$cflxca      = tapply(X=cflxca       ,INDEX=mfac      ,FUN=mean,na.rm=TRUE)
   mont12mn$cflxst      = tapply(X=cflxst       ,INDEX=mfac      ,FUN=mean,na.rm=TRUE)
   mont12mn$nee         = tapply(X=nee          ,INDEX=mfac      ,FUN=mean,na.rm=TRUE)
   mont12mn$reco        = tapply(X=reco         ,INDEX=mfac      ,FUN=mean,na.rm=TRUE)
   mont12mn$hflxca      = tapply(X=hflxca       ,INDEX=mfac      ,FUN=mean,na.rm=TRUE)
   mont12mn$hflxlc      = tapply(X=hflxlc       ,INDEX=mfac      ,FUN=mean,na.rm=TRUE)
   mont12mn$hflxwc      = tapply(X=hflxwc       ,INDEX=mfac      ,FUN=mean,na.rm=TRUE)
   mont12mn$hflxgc      = tapply(X=hflxgc       ,INDEX=mfac      ,FUN=mean,na.rm=TRUE)
   mont12mn$wflxca      = tapply(X=wflxca       ,INDEX=mfac      ,FUN=mean,na.rm=TRUE)
   mont12mn$qwflxca     = tapply(X=qwflxca      ,INDEX=mfac      ,FUN=mean,na.rm=TRUE)
   mont12mn$wflxlc      = tapply(X=wflxlc       ,INDEX=mfac      ,FUN=mean,na.rm=TRUE)
   mont12mn$wflxwc      = tapply(X=wflxwc       ,INDEX=mfac      ,FUN=mean,na.rm=TRUE)
   mont12mn$wflxgc      = tapply(X=wflxgc       ,INDEX=mfac      ,FUN=mean,na.rm=TRUE)
   mont12mn$evap        = tapply(X=evap         ,INDEX=mfac      ,FUN=mean,na.rm=TRUE)
   mont12mn$transp      = tapply(X=transp       ,INDEX=mfac      ,FUN=mean,na.rm=TRUE)
   mont12mn$rain        = tapply(X=rain         ,INDEX=mfac      ,FUN=mean,na.rm=TRUE)
   mont12mn$atm.temp    = tapply(X=atm.temp     ,INDEX=mfac      ,FUN=mean,na.rm=TRUE)
   mont12mn$rshort      = tapply(X=rshort       ,INDEX=mfac      ,FUN=mean,na.rm=TRUE)
   mont12mn$rlong       = tapply(X=rlong        ,INDEX=mfac      ,FUN=mean,na.rm=TRUE)
   mont12mn$atm.shv     = tapply(X=atm.shv      ,INDEX=mfac      ,FUN=mean,na.rm=TRUE)
   mont12mn$atm.co2     = tapply(X=atm.co2      ,INDEX=mfac      ,FUN=mean,na.rm=TRUE)
   mont12mn$atm.prss    = tapply(X=atm.prss     ,INDEX=mfac      ,FUN=mean,na.rm=TRUE)
   mont12mn$atm.vels    = tapply(X=atm.vels     ,INDEX=mfac      ,FUN=mean,na.rm=TRUE)
   mont12mn$ustar       = tapply(X=ustar        ,INDEX=mfac      ,FUN=mean,na.rm=TRUE)
   mont12mn$soil.temp   = qapply(X=soil.temp    ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   mont12mn$soil.water  = qapply(X=soil.water   ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   mont12mn$soil.mstpot = qapply(X=soil.mstpot  ,INDEX=mfac,DIM=1,FUN=mean,na.rm=TRUE)
   #----- Find the mean sum of squares. ---------------------------------------------------#
   print ("      * Aggregating the monthly mean sum of squares...")
   mont12sq           = list()
   mont12sq$gpp       = tapply(X=mmsqu.gpp       ,INDEX=mfac     ,FUN=mean,na.rm=TRUE)
   mont12sq$plresp    = tapply(X=mmsqu.plresp    ,INDEX=mfac     ,FUN=mean,na.rm=TRUE)
   mont12sq$leaf.resp = tapply(X=mmsqu.leaf.resp ,INDEX=mfac     ,FUN=mean,na.rm=TRUE)
   mont12sq$root.resp = tapply(X=mmsqu.root.resp ,INDEX=mfac     ,FUN=mean,na.rm=TRUE)
   mont12sq$hetresp   = tapply(X=mmsqu.hetresp   ,INDEX=mfac     ,FUN=mean,na.rm=TRUE)
   mont12sq$cflxca    = tapply(X=mmsqu.cflxca    ,INDEX=mfac     ,FUN=mean,na.rm=TRUE)
   mont12sq$cflxst    = tapply(X=mmsqu.cflxst    ,INDEX=mfac     ,FUN=mean,na.rm=TRUE)
   mont12sq$hflxca    = tapply(X=mmsqu.hflxca    ,INDEX=mfac     ,FUN=mean,na.rm=TRUE)
   mont12sq$hflxlc    = tapply(X=mmsqu.hflxlc    ,INDEX=mfac     ,FUN=mean,na.rm=TRUE)
   mont12sq$hflxwc    = tapply(X=mmsqu.hflxwc    ,INDEX=mfac     ,FUN=mean,na.rm=TRUE)
   mont12sq$hflxgc    = tapply(X=mmsqu.hflxgc    ,INDEX=mfac     ,FUN=mean,na.rm=TRUE)
   mont12sq$wflxca    = tapply(X=mmsqu.wflxca    ,INDEX=mfac     ,FUN=mean,na.rm=TRUE)
   mont12sq$qwflxca   = tapply(X=mmsqu.qwflxca   ,INDEX=mfac     ,FUN=mean,na.rm=TRUE)
   mont12sq$wflxlc    = tapply(X=mmsqu.wflxlc    ,INDEX=mfac     ,FUN=mean,na.rm=TRUE)
   mont12sq$wflxwc    = tapply(X=mmsqu.wflxwc    ,INDEX=mfac     ,FUN=mean,na.rm=TRUE)
   mont12sq$wflxgc    = tapply(X=mmsqu.wflxgc    ,INDEX=mfac     ,FUN=mean,na.rm=TRUE)
   mont12sq$evap      = tapply(X=mmsqu.evap      ,INDEX=mfac     ,FUN=mean,na.rm=TRUE)
   mont12sq$transp    = tapply(X=mmsqu.transp    ,INDEX=mfac     ,FUN=mean,na.rm=TRUE)
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
   print ("      * Finding the standard deviation...")
   srnorm1 = sqrt(1./(1. - 1. / moncensus))
   srnorm1[!is.finite(srnorm1)] = 0.
   mont12sd            = list()
   mont12sd$gpp        = sqrt(mont12sq$gpp        - mont12mn$gpp^2        ) * srnorm1
   mont12sd$plresp     = sqrt(mont12sq$plresp     - mont12mn$plresp^2     ) * srnorm1
   mont12sd$leaf.resp  = sqrt(mont12sq$leaf.resp  - mont12mn$leaf.resp^2  ) * srnorm1
   mont12sd$root.resp  = sqrt(mont12sq$root.resp  - mont12mn$root.resp^2  ) * srnorm1
   mont12sd$hetresp    = sqrt(mont12sq$hetresp    - mont12mn$hetresp^2    ) * srnorm1
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
   #---------------------------------------------------------------------------------------#
   #     Set standard deviations that became NaN to 0.  This usually happens when we run   #
   # the post-processing script when the simulation hasn't run for more than 2 years.  We  #
   # can't find the standard deviation because the number of degrees of freedom becomes 0. #
   #---------------------------------------------------------------------------------------#
   mont12sd$gpp        [!is.finite(mont12mn$gpp        )] = 0.
   mont12sd$plresp     [!is.finite(mont12mn$plresp     )] = 0.
   mont12sd$leaf.resp  [!is.finite(mont12mn$leaf.resp  )] = 0.
   mont12sd$root.resp  [!is.finite(mont12mn$root.resp  )] = 0.
   mont12sd$hetresp    [!is.finite(mont12mn$hetresp    )] = 0.
   mont12sd$cflxca     [!is.finite(mont12mn$cflxca     )] = 0.
   mont12sd$cflxst     [!is.finite(mont12mn$cflxst     )] = 0.
   mont12sd$hflxca     [!is.finite(mont12mn$hflxca     )] = 0.
   mont12sd$hflxlc     [!is.finite(mont12mn$hflxlc     )] = 0.
   mont12sd$hflxlc     [!is.finite(mont12mn$hflxwc     )] = 0.
   mont12sd$hflxgc     [!is.finite(mont12mn$hflxgc     )] = 0.
   mont12sd$wflxca     [!is.finite(mont12mn$wflxca     )] = 0.
   mont12sd$qwflxca    [!is.finite(mont12mn$qwflxca    )] = 0.
   mont12sd$wflxlc     [!is.finite(mont12mn$wflxlc     )] = 0.
   mont12sd$wflxwc     [!is.finite(mont12mn$wflxwc     )] = 0.
   mont12sd$wflxgc     [!is.finite(mont12mn$wflxgc     )] = 0.
   mont12sd$evap       [!is.finite(mont12mn$evap       )] = 0.
   mont12sd$transp     [!is.finite(mont12mn$transp     )] = 0.
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
   print ("    - Aggregating the monthly mean of the diurnal cycle...")
   dcyc12mn             =list()
   dcyc12mn$gpp         =qapply(X=dcycmean$gpp         ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$npp         =qapply(X=dcycmean$npp         ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$plresp      =qapply(X=dcycmean$plresp      ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$leaf.resp   =qapply(X=dcycmean$leaf.resp   ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$root.resp   =qapply(X=dcycmean$root.resp   ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$hetresp     =qapply(X=dcycmean$hetresp     ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
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
   dcyc12mn$rlong       =qapply(X=dcycmean$rlong       ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$rshort.gnd  =qapply(X=dcycmean$rshort.gnd  ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$rlong.gnd   =qapply(X=dcycmean$rlong.gnd   ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$rlongup     =qapply(X=dcycmean$rlongup     ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$albedo      =qapply(X=dcycmean$albedo      ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$albedo.beam =qapply(X=dcycmean$albedo.beam ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$albedo.diff =qapply(X=dcycmean$albedo.diff ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12mn$rlong.albedo=qapply(X=dcycmean$rlong.albedo,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)

   #----- Find the mean sum of squares. ---------------------------------------------------#
   print ("    - Aggregating the monthly mean sum of squares...")
   dcyc12sq            = list()
   dcyc12sq$gpp        = qapply(X=dcycmsqu$gpp       ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12sq$plresp     = qapply(X=dcycmsqu$plresp    ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12sq$leaf.resp  = qapply(X=dcycmsqu$leaf.resp ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12sq$root.resp  = qapply(X=dcycmsqu$root.resp ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
   dcyc12sq$hetresp    = qapply(X=dcycmsqu$hetresp   ,INDEX=mfac,DIM=1,FUN=mean,na.rm=T)
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
   print ("    - Finding the standard deviation...")
   srnorm1 = sqrt(1./(1. - 1. / moncnt))
   srnorm1[!is.finite(srnorm1)] = 0.
   dcyc12sd            = list()
   dcyc12sd$gpp        = sqrt(dcyc12sq$gpp      -dcyc12mn$gpp^2             )*srnorm1
   dcyc12sd$plresp     = sqrt(dcyc12sq$plresp   -dcyc12mn$plresp^2          )*srnorm1
   dcyc12sd$leaf.resp  = sqrt(dcyc12sq$leaf.resp-dcyc12mn$leaf.resp^2       )*srnorm1
   dcyc12sd$root.resp  = sqrt(dcyc12sq$root.resp-dcyc12mn$root.resp^2       )*srnorm1
   dcyc12sd$hetresp    = sqrt(dcyc12sq$hetresp  -dcyc12mn$hetresp^2         )*srnorm1
   dcyc12sd$nep        = sqrt(dcyc12sq$nep      -dcyc12mn$nep^2             )*srnorm1
   dcyc12sd$cflxca     = sqrt(dcyc12sq$cflxca   -dcyc12mn$cflxca^2          )*srnorm1
   dcyc12sd$cflxst     = sqrt(dcyc12sq$cflxst   -dcyc12mn$cflxst^2          )*srnorm1
   dcyc12sd$hflxca     = sqrt(dcyc12sq$hflxca   -dcyc12mn$hflxca^2          )*srnorm1
   dcyc12sd$hflxlc     = sqrt(dcyc12sq$hflxlc   -dcyc12mn$hflxlc^2          )*srnorm1
   dcyc12sd$hflxwc     = sqrt(dcyc12sq$hflxwc   -dcyc12mn$hflxwc^2          )*srnorm1
   dcyc12sd$hflxgc     = sqrt(dcyc12sq$hflxgc   -dcyc12mn$hflxgc^2          )*srnorm1
   dcyc12sd$wflxca     = sqrt(dcyc12sq$wflxca   -dcyc12mn$wflxca^2          )*srnorm1
   dcyc12sd$qwflxca    = sqrt(dcyc12sq$qwflxca  -dcyc12mn$qwflxca^2         )*srnorm1
   dcyc12sd$wflxlc     = sqrt(dcyc12sq$wflxlc   -dcyc12mn$wflxlc^2          )*srnorm1
   dcyc12sd$wflxwc     = sqrt(dcyc12sq$wflxwc   -dcyc12mn$wflxwc^2          )*srnorm1
   dcyc12sd$wflxgc     = sqrt(dcyc12sq$wflxgc   -dcyc12mn$wflxgc^2          )*srnorm1
   dcyc12sd$transp     = sqrt(dcyc12sq$transp   -dcyc12mn$transp^2          )*srnorm1
   #---------------------------------------------------------------------------------------#
   #     Set standard deviations that became NaN to 0.  This usually happens when we run   #
   # the post-processing script when the simulation hasn't run for more than 2 years.  We  #
   # can't find the standard deviation because the number of degrees of freedom becomes 0. #
   #---------------------------------------------------------------------------------------#
   dcyc12sd$gpp        [!is.finite(dcyc12sd$gpp       )] = 0.
   dcyc12sd$plresp     [!is.finite(dcyc12sd$plresp    )] = 0.
   dcyc12sd$leaf.resp  [!is.finite(dcyc12sd$leaf.resp )] = 0.
   dcyc12sd$root.resp  [!is.finite(dcyc12sd$root.resp )] = 0.
   dcyc12sd$hetresp    [!is.finite(dcyc12sd$hetresp   )] = 0.
   dcyc12sd$nep        [!is.finite(dcyc12sd$nep       )] = 0.
   dcyc12sd$cflxca     [!is.finite(dcyc12sd$cflxca    )] = 0.
   dcyc12sd$cflxst     [!is.finite(dcyc12sd$cflxst    )] = 0.
   dcyc12sd$hflxca     [!is.finite(dcyc12sd$hflxca    )] = 0.
   dcyc12sd$hflxlc     [!is.finite(dcyc12sd$hflxlc    )] = 0.
   dcyc12sd$hflxwc     [!is.finite(dcyc12sd$hflxwc    )] = 0.
   dcyc12sd$hflxgc     [!is.finite(dcyc12sd$hflxgc    )] = 0.
   dcyc12sd$wflxca     [!is.finite(dcyc12sd$wflxca    )] = 0.
   dcyc12sd$qwflxca    [!is.finite(dcyc12sd$qwflxca   )] = 0.
   dcyc12sd$wflxlc     [!is.finite(dcyc12sd$wflxlc    )] = 0.
   dcyc12sd$wflxwc     [!is.finite(dcyc12sd$wflxwc    )] = 0.
   dcyc12sd$wflxgc     [!is.finite(dcyc12sd$wflxgc    )] = 0.
   dcyc12sd$transp     [!is.finite(dcyc12sd$transp    )] = 0.
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
   empty = nplantpftdbh == 0
   agbpftdbh      [empty] = NA
   laipftdbh      [empty] = NA
   waipftdbh      [empty] = NA
   taipftdbh      [empty] = NA
   gpppftdbh      [empty] = NA
   npppftdbh      [empty] = NA
   mcopftdbh      [empty] = NA
   cbapftdbh      [empty] = NA
   ldrpftdbh      [empty] = NA
   fsopftdbh      [empty] = NA
   demandpftdbh   [empty] = NA
   supplypftdbh   [empty] = NA
   nplantpftdbh   [empty] = NA
   #---------------------------------------------------------------------------------------#



   #----- Find which PFTs, land uses and transitions we need to consider ------------------#
   pftave  = colMeans(agbpft,na.rm=TRUE)
   luave   = colMeans(agblu ,na.rm=TRUE)
   distave = matrix(NA,nrow=3,ncol=3)
   for (jlu in 1:nlu){
      for (ilu in 1:nlu){
          distave[ilu,jlu] = mean(dist[,ilu,jlu])
      }#end for
   }#end for
   selpft  = pftave  > 0.
   sellu   = luave   > 0.
   seldist = distave > 0.


   #----- Determine the last data available. ----------------------------------------------#
   tlast = length(thismonth)

   #---------------------------------------------------------------------------------------#
   #      Define a suitable scale for those time series that uses thismonth...             #
   #---------------------------------------------------------------------------------------#
   whenplot6 = pretty.time(thismonth,n=6)
   whenplot8 = pretty.time(thismonth,n=8)
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
   print ("    - Plotting figures...")
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
         print (paste("      +",description,"time series for all PFTs..."))

         #----- Load variable -------------------------------------------------------------#
         thisvar = get(vnam)
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
            if (plog){
               xylog   = "y"
               ylimit  = range(log(thisvar[,selpft]),na.rm=TRUE)
               if (any(! is.finite(ylimit)) || (ylimit[1] == ylimit[2] && ylimit[1] == 0)){
                  ylimit = c(-1,1)
               }else if (ylimit[1] == ylimit[2] ){
                  ylimit[1] = ylimit[1] * ( 1. - sign(ylimit[1]) * ylnudge)
                  ylimit[2] = ylimit[2] * ( 1. + sign(ylimit[2]) * ylnudge)
               }else{
                  ylimit[2] = ylimit[2] + scalleg * (ylimit[2] - ylimit[1])
               }#end if
               ydrought  = c(ylimit[1] - 0.5 * (ylimit[2]-ylimit[1])
                            ,ylimit[2] + 0.5 * (ylimit[2]-ylimit[1]))
               ylimit    = exp(ylimit)
            }else{
               ylimit  = range(thisvar[,selpft],na.rm=TRUE)
               xylog=""
               if (ylimit[1] == ylimit[2] && ylimit[1] == 0){
                  ylimit = c(-1,1)
               }else if(ylimit[1] == ylimit[2] ){
                  ylimit[1] = ylimit[1] * ( 1. - sign(ylimit[1]) * ylnudge)
                  ylimit[2] = ylimit[2] * ( 1. + sign(ylimit[2]) * ylnudge)
               }#end if
               ylimit[2] = ylimit[2] + scalleg * (ylimit[2] - ylimit[1])
               ydrought  = c(ylimit[1] - 0.5 * (ylimit[2]-ylimit[1])
                            ,ylimit[2] + 0.5 * (ylimit[2]-ylimit[1]))
            }#end if
            #------------------------------------------------------------------------------#



            letitre = paste(description,lieu,sep=" - ")
            cols    = pft$colour[selpft]
            legs    = pft$name  [selpft]
            plot(x=thismonth,y=thisvar[,1],type="n",main=letitre,ylim=ylimit
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
                  lines(thismonth,thisvar[,n],type="l",col=pft$colour[n],lwd=lwidth)
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
   pftuse  = which(apply(X=agbpftdbh,MARGIN=3,FUN=sum,na.rm=TRUE) > 0.)
   for (v in 1:ntspftdbh){
      thistspftdbh   = tspftdbh[[v]]
      vnam        = thistspftdbh$vnam
      description = thistspftdbh$desc
      unit        = thistspftdbh$unit
      plog        = thistspftdbh$plog
      plotit      = thistspftdbh$plt
      
      #----- Load variable ----------------------------------------------------------------#
      thisvar = get(vnam)
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

         print (paste("      +",description,"time series for DBH class..."))


         #---------------------------------------------------------------------------------#
         #     Find the limit, make some room for the legend, and in case the field is a   #
         # constant, nudge the limits so the plot command will not complain.               #
         #---------------------------------------------------------------------------------#
         if (plog){
            xylog   = "y"
            ylimit  = range(log(thisvar[,,pftuse]),na.rm=TRUE)
            if (any(! is.finite(ylimit)) || (ylimit[1] == ylimit[2] && ylimit[1] == 0)){
               ylimit = c(-1,1)
            }else if (ylimit[1] == ylimit[2] ){
               ylimit[1] = ylimit[1] * ( 1. - sign(ylimit[1]) * ylnudge)
               ylimit[2] = ylimit[2] * ( 1. + sign(ylimit[2]) * ylnudge)
            }else{
               ylimit[2] = ylimit[2] + scalleg * (ylimit[2] - ylimit[1])
            }#end if
            ydrought  = c(ylimit[1] - 0.5 * (ylimit[2]-ylimit[1])
                         ,ylimit[2] + 0.5 * (ylimit[2]-ylimit[1]))
            ylimit    = exp(ylimit)
         }else{
            ylimit  = range(thisvar[,,pftuse],na.rm=TRUE)
            xylog=""
            if (ylimit[1] == ylimit[2] && ylimit[1] == 0){
               ylimit = c(-1,1)
            }else if(ylimit[1] == ylimit[2] ){
               ylimit[1] = ylimit[1] * ( 1. - sign(ylimit[1]) * ylnudge)
               ylimit[2] = ylimit[2] * ( 1. + sign(ylimit[2]) * ylnudge)
            }#end if
            ylimit[2] = ylimit[2] + scalleg * (ylimit[2] - ylimit[1])
            ydrought  = c(ylimit[1] - 0.5 * (ylimit[2]-ylimit[1])
                         ,ylimit[2] + 0.5 * (ylimit[2]-ylimit[1]))
         }#end if
         #---------------------------------------------------------------------------------#

         for (p in pftuse){

            cpp = substring(100+p,2,3)
            pftlab = paste("pft-",cpp,sep="")

            print (paste("        -",pft$name[p]))


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
               plot(x=thismonth,y=thisvar[,1,p],type="n",main=letitre,ylim=ylimit
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
                  lines(thismonth,thisvar[,d,p],type="l",col=dbhcols[d],lwd=lwidth)
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
   print(paste("    + Comparisons of time series (model vs. observations)..."))
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

      plotit       = plotit && obsnow %in% ls()





      if (plotit){
         #---------------------------------------------------------------------------------#
         #    Copy the observations to a scratch variable.                                 #
         #---------------------------------------------------------------------------------#
         thisobs = get(obsnow)
         mnvar   = paste("emean",vname,sep=".")
         obsmean = thisobs[[mnvar]]
         obswhen = thisobs$tomonth
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #    Check whether the time series directory exists.  If not, create it.          #
         #---------------------------------------------------------------------------------#
         outdir   = paste(outpref,"compemean",sep="/")
         if (! file.exists(outdir)) dir.create(outdir)
         print (paste("      - ",description,"comparison..."))
         #---------------------------------------------------------------------------------#



         #----- Define the number of layers. ----------------------------------------------#
         sel       = thismonth >= min(obswhen) & thismonth <= max(obswhen)
         thiswhen  = thismonth
         thismean  = get(vname)
         #---------------------------------------------------------------------------------# 



         #----- Find the plot range. ------------------------------------------------------#
         ylimit    = range(c(thismean[sel],obsmean),na.rm=TRUE)
         #----- Expand the upper range in so the legend doesn't hide things. --------------#
         if (ylimit[1] == ylimit[2]  & ylimit[1] == 0){
            ylimit[1] = -1
            ylimit[2] =  1
         }else if (ylimit[1] == ylimit[2] & ylimit[1] > 0){
            ylimit[2] = (1.0+scalleg) * ylimit[1]
         }else if (ylimit[1] == ylimit[2] & ylimit[1] < 0){
            ylimit[2] = (1.0-scalleg) * ylimit[1]
         }else{
            ylimit[2] = ylimit[2] + scalleg * (ylimit[2] - ylimit[1])
         }#end if
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
            legend(x=legpos,inset=0.01,legend=c("Model","Observation")
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
   print(paste("    + Comparisons of monthly means (model vs. observations)..."))
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
         #    Check whether the time series directory exists.  If not, create it.          #
         #---------------------------------------------------------------------------------#
         outdir   = paste(outpref,"compmmean",sep="/")
         if (! file.exists(outdir)) dir.create(outdir)
         print (paste("      - ",description,"comparison..."))
         #---------------------------------------------------------------------------------#



         #----- Define the number of layers. ----------------------------------------------#
         thismean  = mont12mn[[vname]]
         thissdev  = mont12sd[[vname]]
         #---------------------------------------------------------------------------------# 



         #---------------------------------------------------------------------------------# 
         #    Some variables have no standard deviation in the model.  Make them 0 if this #
         # is the case.                                                                    #
         #---------------------------------------------------------------------------------# 
         if (length(thissdev) == 0){
            thissdev = 0. * thismean
         }#end if
         #---------------------------------------------------------------------------------# 



         #----- Find the plot range. ------------------------------------------------------#
         if (plotsd){
            ylimit    = range(c(thismean + thissdev ,thismean - thissdev
                               ,obsmean  + obssdev  ,obsmean  - obssdev    ),na.rm=TRUE)
         }else{
            ylimit    = range(c(thismean,obsmean),na.rm=TRUE)
         }#end if
         #----- Expand the upper range in so the legend doesn't hide things. --------------#
         if (ylimit[1] == ylimit[2]  & ylimit[1] == 0){
            ylimit[1] = -1
            ylimit[2] =  1
         }else if (ylimit[1] == ylimit[2] & ylimit[1] > 0){
            ylimit[2] = (1.0+scalleg) * ylimit[1]
         }else if (ylimit[1] == ylimit[2] & ylimit[1] < 0){
            ylimit[2] = (1.0-scalleg) * ylimit[1]
         }else{
            ylimit[2] = ylimit[2] + scalleg * (ylimit[2] - ylimit[1])
         }#end if
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
               err.x = c(montmont,rev(montmont),NA,montmont,rev(montmont))
               err.y = c(thismean + thissdev,rev(thismean) - rev(thissdev),NA
                        ,obsmean  + obssdev ,rev(obsmean ) - rev(obssdev )   )
               polygon(x=err.x,y=err.y,col=errcolours,angle=angle,density=dens
                      ,lty="solid",lwd=shwd)
            }#end if
            points(x=montmont,y=thismean,col=lcolours[1],lwd=llwd[1],type=ltype
                  ,pch=16,cex=1.0)
            points(x=montmont,y=obsmean ,col=lcolours[2],lwd=llwd[2],type=ltype
                  ,pch=16,cex=1.0)
            if (plotsd){
               legend(x=legpos,inset=0.05,legend=c("Model","Observation")
                     ,fill=errcolours,angle=angle,density=dens,lwd=llwd,col=lcolours
                     ,bg="white",title="Shaded areas = 1 SD",cex=1.0,pch=16)
            }else{
               legend(x=legpos,inset=0.05,legend=c("Model","Observation")
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
   print(paste("      * Comparisons of mean diurnal cycle (model vs. observations)..."))
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
         #    Check whether the time series directory exists.  If not, create it.          #
         #---------------------------------------------------------------------------------#
         outdir   = paste(outpref,"compdcyc",sep="/")
         if (! file.exists(outdir)) dir.create(outdir)
         outtheme = paste(outdir,vname,sep="/")
         if (! file.exists(outtheme)) dir.create(outtheme)
         print (paste("      +",description,"comparison..."))
         #---------------------------------------------------------------------------------# 



         #----- Define the number of layers. ----------------------------------------------#
         thismean  = dcyc12mn[[vname]]
         thissdev  = dcyc12sd[[vname]]
         #---------------------------------------------------------------------------------# 



         #---------------------------------------------------------------------------------# 
         #    Some variables have no standard deviation in the model.  Make them 0 if this #
         # is the case.                                                                    #
         #---------------------------------------------------------------------------------# 
         if (length(thissdev) == 0){
            thissdev = 0. * thismean
         }#end if
         #---------------------------------------------------------------------------------# 


         #----- Append the last hour before the first one. --------------------------------#
         thismean  = cbind(thismean[,ndcycle],thismean)
         thissdev  = cbind(thissdev[,ndcycle],thissdev)
         #---------------------------------------------------------------------------------# 


         #----- Find the plot range. ------------------------------------------------------#
         if (plotsd){
            ylimit    = range(c(thismean + thissdev ,thismean - thissdev
                               ,obsmean  + obssdev  ,obsmean  - obssdev    ),na.rm=TRUE)
         }else{
            ylimit    = range(c(thismean,obsmean),na.rm=TRUE)
         }#end if
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
            print (paste("        > ",description," time series - ",namemon,"...",sep=""))

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
                  err.x = c(thisday,rev(thisday),NA,thisday,rev(thisday))
                  err.y = c(thismean[pmon,] + thissdev[pmon,]
                           ,rev(thismean[pmon,]) - rev(thissdev[pmon,])
                           ,NA
                           ,obsmean[pmon,]      + obssdev[pmon,]
                           ,rev(obsmean[pmon,]) - rev(obssdev[pmon,]))
                  polygon(x=err.x,y=err.y,col=errcolours,angle=angle,density=dens
                         ,lty="solid",lwd=shwd)
               }#end if
               points(x=thisday,y=thismean[pmon,],col=lcolours[1]
                     ,lwd=llwd[1],type=ltype,pch=16,cex=1.0)
               points(x=thisday,y=obsmean[pmon,],col=lcolours[2]
                     ,lwd=llwd[2],type=ltype,pch=16,cex=1.0)
               if (plotsd){
                  legend(x=legpos,inset=0.05,legend=c("Model","Observation")
                        ,fill=errcolours,angle=angle,density=dens,lwd=llwd,col=lcolours
                        ,bg="white",title="Shaded areas = 1 SD",cex=1.0,pch=16)
               }else{
                  legend(x=legpos,inset=0.05,legend=c("Model","Observation")
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
         print (paste("      +",description,"time series for all LUs..."))



         #----- Load variable -------------------------------------------------------------#
         thisvar = get(vnam)
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
            if (plog){
               xylog   = "y"
               ylimit  = range(log(thisvar[,selpft]),na.rm=TRUE)
               if (any(! is.finite(ylimit)) || (ylimit[1] == ylimit[2] && ylimit[1] == 0)){
                  ylimit = c(-1,1)
               }else if (ylimit[1] == ylimit[2] ){
                  ylimit[1] = ylimit[1] * ( 1. - sign(ylimit[1]) * ylnudge)
                  ylimit[2] = ylimit[2] * ( 1. + sign(ylimit[2]) * ylnudge)
               }else{
                  ylimit[2] = ylimit[2] + scalleg * (ylimit[2] - ylimit[1])
               }#end if
               ydrought  = c(ylimit[1] - 0.5 * (ylimit[2]-ylimit[1])
                            ,ylimit[2] + 0.5 * (ylimit[2]-ylimit[1]))
               ylimit    = exp(ylimit)
            }else{
               ylimit  = range(thisvar[,selpft],na.rm=TRUE)
               xylog=""
               if (ylimit[1] == ylimit[2] && ylimit[1] == 0){
                  ylimit = c(-1,1)
               }else if(ylimit[1] == ylimit[2] ){
                  ylimit[1] = ylimit[1] * ( 1. - sign(ylimit[1]) * ylnudge)
                  ylimit[2] = ylimit[2] * ( 1. + sign(ylimit[2]) * ylnudge)
               }#end if
               ylimit[2] = ylimit[2] + scalleg * (ylimit[2] - ylimit[1])
               ydrought  = c(ylimit[1] - 0.5 * (ylimit[2]-ylimit[1])
                            ,ylimit[2] + 0.5 * (ylimit[2]-ylimit[1]))
            }#end if
            #------------------------------------------------------------------------------#

            letitre = paste(description,lieu,sep=" - ")
            cols    = lucols[sellu]
            legs    = lunames[sellu]
            plot(thismonth,thisvar[,1],type="n",main=letitre,ylim=ylimit
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
                  lines(thismonth,thisvar[,n],type="l",col=lucols[n],lwd=lwidth)
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
      print (paste("      + Disturbance rate time series for all disturbances..."))
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
               if (seldist[ilu,jlu]) ylimit = range(c(ylimit,dist[,ilu,jlu]),na.rm=TRUE)
            }#end for
         }#end for
         if (ylimit[1] == ylimit[2] && ylimit[1] == 0){
            ylimit = c(-1,1)
         }else if(ylimit[1] == ylimit[2] ){
            ylimit[1] = ylimit[1] * ( 1. - sign(ylimit[1]) * ylnudge)
            ylimit[2] = ylimit[2] * ( 1. + sign(ylimit[2]) * ylnudge)
         }#end if
         ylimit[2] = ylimit[2] + scalleg * (ylimit[2] - ylimit[1])
         ydrought  = c(ylimit[1] - 0.5 * (ylimit[2]-ylimit[1])
                      ,ylimit[2] + 0.5 * (ylimit[2]-ylimit[1]))
         #---------------------------------------------------------------------------------#

         letitre = paste("Disturbance rates",lieu,sep=" - ")
         cols    = NULL
         legs    = NULL
         plot(thismonth,dist[,1,1],type="n",main=letitre,ylim=ylimit
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
                  lines(thismonth,dist[,ilu,jlu],type="l",col=distcols[n],lwd=lwidth)
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
   print(paste("      * Plot some time series figures..."))
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
         print (paste("      +",theme,"time series for several variables..."))


         #----- Define the number of layers. ----------------------------------------------#
         nlayers   = length(vnames)
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Find the limit, make some room for the legend, and in case the field is a   #
         # constant, nudge the limits so the plot command will not complain.               #
         #---------------------------------------------------------------------------------#
         ylimit    = NULL
         for (l in 1:nlayers){
            thisvar = get(vnames[l])
            ylimit  = range(c(ylimit,thisvar),na.rm=TRUE)
         }#end for
         if (plog){
            xylog   = "y"
            ylimit  = range(ylimit,na.rm=TRUE)
            if (any(! is.finite(ylimit)) || (ylimit[1] == ylimit[2] && ylimit[1] == 0)){
               ylimit = c(-1,1)
            }else if (ylimit[1] == ylimit[2] ){
               ylimit[1] = ylimit[1] * ( 1. - sign(ylimit[1]) * ylnudge)
               ylimit[2] = ylimit[2] * ( 1. + sign(ylimit[2]) * ylnudge)
            }else{
               ylimit[2] = ylimit[2] + scalleg * (ylimit[2] - ylimit[1])
            }#end if
            ydrought  = c(ylimit[1] - 0.5 * (ylimit[2]-ylimit[1])
                         ,ylimit[2] + 0.5 * (ylimit[2]-ylimit[1]))
            ylimit    = exp(ylimit)
         }else{
            xylog=""
            if (ylimit[1] == ylimit[2] && ylimit[1] == 0){
               ylimit = c(-1,1)
            }else if(ylimit[1] == ylimit[2] ){
               ylimit[1] = ylimit[1] * ( 1. - sign(ylimit[1]) * ylnudge)
               ylimit[2] = ylimit[2] * ( 1. + sign(ylimit[2]) * ylnudge)
            }#end if
            ylimit[2] = ylimit[2] + scalleg * (ylimit[2] - ylimit[1])
            ydrought  = c(ylimit[1] - 0.5 * (ylimit[2]-ylimit[1])
                         ,ylimit[2] + 0.5 * (ylimit[2]-ylimit[1]))
         }#end if
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Check if the directory exists.  If not, create it.                          #
         #---------------------------------------------------------------------------------#
         print (paste("        > ",theme," time series ...",sep=""))

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
            thisvar = get(vnames[1])

            letitre = paste(theme," - ",lieu," \n"," Time series: ",theme,sep="")

            plot(x=thismonth,y=thisvar,type="n",main=letitre,xlab="Time"
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
               thisvar = get(vnames[l])
               points(x=thismonth,y=thisvar,col=lcolours[l]
                     ,lwd=llwd[l],type=ltype,pch=16,cex=0.8)
            }#end for
            legend(x=legpos,inset=0.05,legend=description,col=lcolours,lwd=llwd,cex=0.8)
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
   print(paste("      * Plot some climatology of diurnal cycle..."))
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
         print (paste("      +",description,"diurnal cycle for several variables..."))


         #----- Define the number of layers. ----------------------------------------------#
         nlayers   = length(vnames)
         ylimit    = NULL
         for (l in 1:nlayers){
            thisvar = dcyc12mn[[vnames[l]]]
            ylimit  = range(c(ylimit,thisvar),na.rm=TRUE)
         }#end for
         if (ylimit[1] == ylimit[2]  & ylimit[1] == 0){
            ylimit[1] = -1
            ylimit[2] =  1
         }else if (ylimit[1] == ylimit[2] & ylimit[1] > 0){
            ylimit[2] = (1.0+scalleg) * ylimit[1]
         }else if (ylimit[1] == ylimit[2] & ylimit[1] < 0){
            ylimit[2] = (1.0-scalleg) * ylimit[1]
         }else{
            ylimit[2] = ylimit[2] + scalleg * (ylimit[2] - ylimit[1])
         }#end if
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
            print (paste("        > ",theme," time series - ",namemon,"...",sep=""))

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
         print (paste("      + Climatology profile of ",description,"..."))

         #----- Find the number of rows and columns, and the axes. ------------------------#
         monaxis  = sort(unique(monnum))
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
         print (paste("      + Time series profile of ",description,"..."))

         #----- Find the number of rows and columns, and the axes. ------------------------#
         timeaxis  = thismonth
         soilaxis  = slz
         nmon      = length(timeaxis)
         nsoil     = nzg

         #----- Convert the vector data into an array. ------------------------------------#
         vararr  = get(vnam)

         #----- Copy Decembers ans Januaries to make the edges buffered. ------------------#
         first    = vararr[1,]
         first    = c(first,first[nzg],first[nzg])

         last     = vararr[totmon,]
         last     = c(last[1],last[1],last)

         #----- Bind first and last year to the array, to make the edges buffered. --------#
         varbuff  = cbind(vararr[,1],vararr,vararr[,nzg])
         varbuff  = rbind(first,varbuff,last)

         #---------------------------------------------------------------------------------#
         #      Expand the month and year axes.  Make the first and last time equal time   #
         # steps.                                                                          #
         #---------------------------------------------------------------------------------#
         dwhen    = as.numeric(thismonth[2]-thismonth[1])
         whenaxis = c(chron(as.numeric(thismonth[1]-dwhen))
                     ,timeaxis
                     ,chron(as.numeric(thismonth[totmon]+dwhen)))
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
      sel      = myear >= yearaa & myear <= yearzz
      twoyears = sum(sel) >= 24

      if (plotit && twoyears){

         #---------------------------------------------------------------------------------#
         #     Check if the directory exists.  If not, create it.                          #
         #---------------------------------------------------------------------------------#
         outdir  =  paste(outpref,"hovmoller",sep="/")
         if (! file.exists(outdir)) dir.create(outdir)
         print (paste("      +",description,"Hovmoller time series ..."))

         #----- Load this variable into "thisvar". ----------------------------------------#
         thisvar = get(vnam)

         #----- Find the number of rows and columns, and the axes. ------------------------#
         monaxis = sort(unique(monnum[sel]))
         yraxis  = sort(unique(myear[sel]))
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
         print (paste("      +",description,"time series of diurnal cycle..."))

         #----- Load this variable into "thisvar". ----------------------------------------#
         vararr   = dcycmean[[vnam]]

         #----- Copy Decembers ans Januaries to make the edges buffered. ------------------#
         firsthr  = vararr[,1]
         firsthr  = c(firsthr,firsthr[totmon],firsthr[totmon])

         lasthr   = vararr[,ndcycle]
         lasthr   = c(lasthr[1],lasthr[1],lasthr)

         #----- Bind first and last year to the array, to make the edges buffered. --------#
         varbuff  = rbind(vararr[1,],vararr,vararr[totmon,])
         varbuff  = cbind(lasthr,varbuff,firsthr)

         #----- Expand the month and year axes. -------------------------------------------#
         hraxis    = seq(from=0,to=ndcycle+1,by=1) * 24 / ndcycle
         dwhen     = thismonth[2]-thismonth[1]
         whenaxis  = c(thismonth[1]-dwhen,thismonth,thismonth[totmon]+dwhen)
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
      thisbplot   = bplot[[v]]
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
         print (paste("      +",description,"box plot..."))

         #----- Load this variable into "thisvar". ----------------------------------------#
         thisvar = get(vnam)

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
            ylimit  = range(thisvar, na.rm=TRUE)
            if (ylimit[1] == ylimit[2] && ylimit[1] == 0){
               ylimit = c(-1,1)
            }else if(ylimit[1] == ylimit[2] ){
               ylimit[1] = ylimit[1] * ( 1. - sign(ylimit[1]) * ylnudge)
               ylimit[2] = ylimit[2] * ( 1. + sign(ylimit[2]) * ylnudge)
            }#end if
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
   print(paste("    + Bar plot by DBH classes..."))
   monbplot    = which(nummonths(thismonth) %in% sasmonth)
   nmonbplot   = length(monbplot)
   pftuse      = which(apply(X=nplantpftdbh,MARGIN=3,FUN=sum,na.rm=TRUE) > 0.)
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
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #      Check whether to plot this 
      #------------------------------------------------------------------------------------#
      if (plotit){
         print (paste("      - ",description,"..."))


         #---------------------------------------------------------------------------------#
         #     Retrieve the variable, and keep only the part that is usable.               #
         #---------------------------------------------------------------------------------#
         thisvnam                  = get(vnam)[monbplot,,]
         thisvnam                  = thisvnam [,,pftuse]
         thisvnam                  = thisvnam [,-(ndbh+1),]
         
         thisvnam[is.na(thisvnam)] = 0.
         thiswhen                  = thismonth[monbplot]
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
         #----- Expand the upper range in so the legend doesn't hide things. --------------#
         if (ylimit[1] == ylimit[2]  & ylimit[1] == 0){
            ylimit[1] = -1
            ylimit[2] =  1
         }else if (ylimit[1] == ylimit[2] & ylimit[1] > 0){
            ylimit[2] = (1.0+scalleg) * ylimit[1]
         }else if (ylimit[1] == ylimit[2] & ylimit[1] < 0){
            ylimit[2] = (1.0-scalleg) * ylimit[1]
         }else{
            ylimit[2] = ylimit[2] + scalleg * (ylimit[2] - ylimit[1])
         }#end if
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
                      ,beside=(! stacked),col=pftcol.use
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

         print (paste("      +",description,"size and age structure plot..."))

         #---------------------------------------------------------------------------------#
         #     Check if the directory exists.  If not, create it.                          #
         #---------------------------------------------------------------------------------#
         sasdir = paste(outpref,"sas",sep="/")
         if (! file.exists(sasdir)) dir.create(sasdir)
         outdir = paste(sasdir,vnam,sep="/")
         if (! file.exists(outdir)) dir.create(outdir)

         #----- Load this list into "thislist". -------------------------------------------#
         varco =  get(vnam)


         for (ww in names(ageco)){

            #----- Find which year we are plotting. ---------------------------------------#
            cmonth   = substring(ww,7,8)
            thisyear = substring(ww,2,5)
            mm       = as.numeric(cmonth)
            yy       = as.numeric(thisyear)

            #----- Retrieve variable list, age, DBH, and PFT for this year. ---------------#
            ageww   = ageco[[ww]]
            dbhww   = dbhco[[ww]]
            pftww   = pftco[[ww]]
            varww   = varco[[ww]]
            popww   = nplantco[[ww]] * areaco[[ww]]

            #------------------------------------------------------------------------------#
            #     We only plot the SAS figures when the polygon is not an absolute desert. #
            #------------------------------------------------------------------------------#
            if (any (! is.na(varww))){
               #---------------------------------------------------------------------------#
               #      Find the range.  If the user wants the range to be fixed, then use   #
               # the global range, otherwise, simply use the range for this year.          #
               #---------------------------------------------------------------------------#
               if (sasfixlimits){
                  xlimit  = range(ageco            , na.rm=TRUE)
                  ylimit  = range(dbhco            , na.rm=TRUE)
                  zlimit  = range(varco            , na.rm=TRUE)
                  popmin  = min  (nplantco * areaco, na.rm=TRUE)
                  popmax  = max  (nplantco * areaco, na.rm=TRUE)
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
                  legend(x="bottomright",inset=0.01,legend=pftleg,fill=colleg
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
