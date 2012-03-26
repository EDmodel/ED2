#----- Here is the user-defined variable section. -----------------------------------------#
here           = "thispath"                                  # Current directory.
srcdir         = "/n/moorcroft_data/mlongo/util/Rsc"      # Source  directory.
outroot        = "thisoutroot" # Source  directory.
myplaces       = c("thispoly")
#------------------------------------------------------------------------------------------#
#     Initial and final times, they must be character vectors of size 2, the first one     #
# with m/d/y, and the second one with h:m:s".                                              #
#------------------------------------------------------------------------------------------#
whena          = c("thismontha/thisdatea/thisyeara","thishoura:thisminua:00")
whenz          = c("thismonthz/thisdatez/thisyearz","thishourz:thisminuz:00")
ptype          = "l"                  # Type of plot
ptyped         = "p"                  # Type of plot
ptypeb         = "o"                  # Type of plot

outform        = "png"           # Formats for output file.  Supported formats are:
                                 #   - "X11" - for printing on screen
                                 #   - "eps" - for postscript printing
                                 #   - "png" - for PNG printing

cex.main       = 0.8             # Scale coefficient for the title

byeold         = TRUE           # Remove old files of the given format?

depth          = 96             # PNG resolution, in pixels per inch
paper          = "letter"       # Paper size, to define the plot shape
ptsz           = 14             # Font size.
lwidth         = 2.5            # Line width
plotgrid       = TRUE           # Should I plot the grid in the background? 

legwhere       = "topleft"      # Where should I place the legend?
inset          = 0.05           # inset distance between legend and edge of plot region.
legbg          = "white"        # Legend background colour.

scalleg        = 0.32           # Increase in y scale to fit the legend.
ncolshov       = 200            # Target number of colours for Hovmoller diagrams.
hovgrid        = TRUE           # Should I include a grid on the Hovmoller plots?

#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#     List of possible plots. In case you don't want some of them, simply switch plt to F. #
#------------------------------------------------------------------------------------------#
#----- Time series plots. -----------------------------------------------------------------#
nphov = 28
phovdi01 = list(vnam   = c("gpp","plresp","hetresp","cflxac")
               ,desc   = c("GPP","Plant resp.","Het. resp.","Atm->Canopy")
               ,colour = c("forestgreen","chartreuse","sienna","deepskyblue")
               ,lwd    = c(2.0,2.0,2.0,2.0)
               ,type   = ptype
               ,plog   = ""
               ,prefix = "carbflux"
               ,theme  = "Carbon fluxes"
               ,unit   = "umol/m2/s"
               ,legpos = "topleft"
               ,plt    = TRUE)
phovdi02 = list(vnam   = c("met.rshort","met.rlong","gnd.rshort"
                          ,"gnd.rlong","qwflxca","hflxca")
               ,desc   = c("Met SW","Met LW","Ground SW","Ground LW","Latent","Sensible")
               ,colour = c("goldenrod","limegreen","deepskyblue","chartreuse"
                          ,"midnightblue","firebrick")
               ,lwd    = c(2.0,2.0,2.0,2.0,2.0,2.0)
               ,type   = ptype
               ,plog   = ""
               ,prefix = "eneflux"
               ,theme  = "Energy fluxes"
               ,unit   = "W/m2"
               ,legpos = "topleft"
               ,plt    = TRUE)
phovdi03 = list(vnam   = c("wflxgc","wflxca","wflxlc","wflxwc","transp")
               ,desc   = c("Ground->Canopy","Canopy->Air","Leaf->Canopy","Wood->Canopy"
                          ,"Transpiration")
               ,colour = c("firebrick","midnightblue","chartreuse","goldenrod"
                          ,"darkolivegreen")
               ,lwd    = c(2.0,2.0,2.0,2.0,2.0)
               ,type   = ptype
               ,plog   = ""
               ,prefix = "h2oflux"
               ,theme  = "Water fluxes"
               ,unit   = "kg/m2/day"
               ,legpos = "topleft"
               ,plt    = TRUE)
phovdi04 = list(vnam   = c("hflxgc","hflxca","hflxlc","hflxwc")
               ,desc   = c("Ground->Canopy","Canopy->Air","Leaf->Canopy","Wood->Canopy")
               ,colour = c("firebrick","midnightblue","chartreuse","goldenrod")
               ,lwd    = c(2.0,2.0,2.0,2.0)
               ,type   = ptype
               ,plog   = ""
               ,prefix = "sensflux"
               ,theme  = "Sensible heat fluxes"
               ,unit   = "W/m2"
               ,legpos = "topleft"
               ,plt    = TRUE)
phovdi05 = list(vnam   = c("atm.temp","can.temp","leaf.temp","wood.temp","sfc.temp")
               ,desc   = c("Atmosphere","Canopy air","Leaf","Wood","Surface")
               ,colour = c("deepskyblue","gray21","chartreuse","goldenrod","sienna")
               ,lwd    = c(2.0,2.0,2.0,2.0,2.0)
               ,type   = ptype
               ,plog   = ""
               ,prefix = "temperature"
               ,theme  = "Temperature"
               ,unit   = "degC"
               ,legpos = "topleft"
               ,plt    = TRUE)
phovdi06 = list(vnam   = c("atm.shv","can.shv","sfc.shv")
               ,desc   = c("Atmosphere","Canopy air","Surface")
               ,colour = c("deepskyblue","gray21","sienna")
               ,lwd    = c(2.0,2.0,2.0)
               ,type   = ptype
               ,plog   = ""
               ,prefix = "h2ovapour"
               ,theme  = "Water vapour mixing ratio"
               ,unit   = "g/kg"
               ,legpos = "topleft"
               ,plt    = TRUE)
phovdi07 = list(vnam   = c("atm.relhum","can.relhum")
               ,desc   = c("Atmosphere","Canopy air")
               ,colour = c("aquamarine","navy")
               ,lwd    = c(2.0,2.0)
               ,type   = ptype
               ,plog   = ""
               ,prefix = "relhum"
               ,theme  = "Relative humidity"
               ,unit   = "%"
               ,legpos = "topleft"
               ,plt    = TRUE)
phovdi08 = list(vnam   = c("atm.supsat","can.supsat")
               ,desc   = c("Atmosphere","Canopy air")
               ,colour = c("steelblue","navy")
               ,lwd    = c(2.0,2.0)
               ,type   = ptype
               ,plog   = ""
               ,prefix = "supsat"
               ,theme  = "Super-saturation"
               ,unit   = "%"
               ,legpos = "topleft"
               ,plt    = TRUE)
phovdi09 = list(vnam   = c("atm.co2","can.co2")
               ,desc   = c("Atmosphere","Canopy air")
               ,colour = c("chartreuse","darkolivegreen")
               ,lwd    = c(2.0,2.0)
               ,type   = ptype
               ,plog   = ""
               ,prefix = "co2con"
               ,theme  = "CO2 mixing ratio"
               ,unit   = "umol/mol"
               ,legpos = "topleft"
               ,plt    = TRUE)
phovdi10 = list(vnam   = c("atm.prss","can.prss")
               ,desc   = c("Atmosphere","Canopy air")
               ,colour = c("violetred3","purple3")
               ,lwd    = c(2.0,2.0)
               ,type   = ptype
               ,plog   = ""
               ,prefix = "press"
               ,theme  = "Pressure"
               ,unit   = "hPa"
               ,legpos = "topleft"
               ,plt    = TRUE)
phovdi11 = list(vnam   = c("atm.theiv","can.theiv")
               ,desc   = c("Atmosphere","Canopy air")
               ,colour = c("orange","firebrick")
               ,lwd    = c(2.0,2.0)
               ,type   = ptype
               ,plog   = ""
               ,prefix = "theiv"
               ,theme  = "Equivalent potential temperature"
               ,unit   = "K"
               ,legpos = "topleft"
               ,plt    = TRUE)
phovdi12 = list(vnam   = c("atm.thetav","can.thetav")
               ,desc   = c("Atmosphere","Canopy air")
               ,colour = c("lawngreen","forestgreen")
               ,lwd    = c(2.0,2.0)
               ,type   = ptype
               ,plog   = ""
               ,prefix = "thetav"
               ,theme  = "Virtual potential temperature"
               ,unit   = "K"
               ,legpos = "topleft"
               ,plt    = TRUE)
phovdi13 = list(vnam   = c("prec","intercept","throughfall","wshed")
               ,desc   = c("Precipitation","Interception","Throughfall","Shedding")
               ,colour = c("midnightblue","forestgreen","dodgerblue","aquamarine")
               ,lwd    = c(2.5,2.0,2.0,2.0)
               ,type   = ptype
               ,plog   = ""
               ,prefix = "prec"
               ,theme  = "Precipitation rate"
               ,unit   = "mm/hr"
               ,legpos = "topleft"
               ,plt    = TRUE)
phovdi14 = list(vnam   = c("atm.vels","ustar")
               ,desc   = c("Wind speed","Friction vel.")
               ,colour = c("midnightblue","goldenrod")
               ,lwd    = c(2.0,2.0)
               ,type   = ptype
               ,plog   = ""
               ,prefix = "wind"
               ,theme  = "Wind speed"
               ,unit   = "m/s"
               ,legpos = "topleft"
               ,plt    = TRUE)
phovdi15 = list(vnam   = c("soilcp","soilwp","soilfc","slmsts","soil.water")
               ,desc   = c("Dry soil","Wilting Point","Field capacity","Saturation"
                          ,"Soil moisture")
               ,colour = c("firebrick","goldenrod","steelblue","midnightblue","olivedrab")
               ,lwd    = c(2.0,2.0,2.0,2.0,2.0)
               ,type   = ptype
               ,plog   = ""
               ,prefix = "smoist"
               ,theme  = "Soil moisture"
               ,unit   = "m3/m3"
               ,legpos = "topleft"
               ,plt    = TRUE)
phovdi16 = list(vnam   = c("zeta","ri.bulk")
               ,desc   = c("Height scale","Richardson")
               ,colour = c("goldenrod","steelblue")
               ,lwd    = c(2.0,2.0)
               ,type   = ptype
               ,plog   = ""
               ,prefix = "turb"
               ,theme  = "Turbulence terms"
               ,unit   = "n/d"
               ,legpos = "topleft"
               ,plt    = TRUE)
phovdi17 = list(vnam   = c("ksn")
               ,desc   = c("Ponding layers")
               ,colour = c("steelblue")
               ,lwd    = c(2.0)
               ,type   = ptype
               ,plog   = ""
               ,prefix = "pond"
               ,theme  = "Ponding layers"
               ,unit   = "layers"
               ,legpos = "topleft"
               ,plt    = TRUE)
phovdi18 = list(vnam   = c("hdid")
               ,desc   = c("Time step")
               ,colour = c("midnightblue")
               ,lwd    = c(2.0)
               ,type   = ptypeb
               ,plog   = "y"
               ,prefix = "tstep"
               ,theme  = "Time step"
               ,unit   = "s"
               ,legpos = "topleft"
               ,plt    = TRUE)
phovdi19 = list(vnam   = c("flag.wflxgc")
               ,desc   = c("Flag")
               ,colour = c("purple")
               ,lwd    = c(2.0)
               ,type   = ptyped
               ,plog   = ""
               ,prefix = "flagwflx"
               ,theme  = "Flag of water flux"
               ,unit   = "--"
               ,legpos = "topleft"
               ,plt    = TRUE)
phovdi20 = list(vnam   = c("atm.rhos","can.rhos")
               ,desc   = c("Atmosphere","Canopy air")
               ,colour = c("sienna","goldenrod")
               ,lwd    = c(2.0,2.0)
               ,type   = ptype
               ,plog   = ""
               ,prefix = "dens"
               ,theme  = "Density"
               ,unit   = "kg/m3"
               ,legpos = "topleft"
               ,plt    = TRUE)
phovdi21 = list(vnam   = c("atm.idgas","can.idgas")
               ,desc   = c("Atmosphere","Canopy air")
               ,colour = c("midnightblue","deepskyblue")
               ,lwd    = c(2.0,2.0)
               ,type   = ptype
               ,plog   = ""
               ,prefix = "idgas"
               ,theme  = "Ideal gas departure"
               ,unit   = "%"
               ,legpos = "topleft"
               ,plt    = TRUE)
phovdi22 = list(vnam   = c("flag.sfcwater")
               ,desc   = c("Flag")
               ,colour = c("steelblue")
               ,lwd    = c(2.0)
               ,type   = ptyped
               ,plog   = ""
               ,prefix = "flagsfcw"
               ,theme  = "Flag of temporary surface water"
               ,unit   = "--"
               ,legpos = "topleft"
               ,plt    = TRUE)
phovdi23 = list(vnam   = c("cum.step")
               ,desc   = c("Time step")
               ,colour = c("darkolivegreen")
               ,lwd    = c(2.0)
               ,type   = ptype
               ,plog   = ""
               ,prefix = "cum.step"
               ,theme  = "Cumulative number of steps"
               ,unit   = "#"
               ,legpos = "topleft"
               ,plt    = TRUE)
phovdi24 = list(vnam   = c("ggbare","ggveg","ggnet","ggold")
               ,desc   = c("Bare","Vegetation","Net","0.2*u*")
               ,colour = c("darkorange3","lawngreen","goldenrod","steelblue")
               ,lwd    = c(2.0,2.0,2.0)
               ,type   = ptype
               ,plog   = ""
               ,prefix = "ggnd"
               ,theme  = "Ground to canopy conductance"
               ,unit   = "m/s"
               ,legpos = "topleft"
               ,plt    = TRUE)
phovdi25 = list(vnam   = c("rgbare","rgveg","rgnet","rgold")
               ,desc   = c("Bare","Vegetation","Net","0.2*u*")
               ,colour = c("darkorange3","lawngreen","goldenrod","steelblue")
               ,lwd    = c(2.0,2.0,2.0,2.0)
               ,type   = ptype
               ,plog   = ""
               ,prefix = "rgnd"
               ,theme  = "Ground to canopy resistance"
               ,unit   = "s/m"
               ,legpos = "topleft"
               ,plt    = TRUE)
phovdi26 = list(vnam   = c("par.top","par.beam.top","par.diff.top"
                          ,"par.bot","par.beam.bot","par.diff.bot")
               ,desc   = c("Top","Beam Top","Diffuse Top"
                          ,"Bottom","Beam Bottom","Diffuse Bottom")
               ,colour = c("firebrick","sienna","gold"
                          ,"midnightblue","steelblue","deepskyblue")
               ,lwd    = c(2.0,2.0,2.0,2.0,2.0,2.0)
               ,type   = ptype
               ,plog   = ""
               ,prefix = "parrad"
               ,theme  = "PAR radiation"
               ,unit   = "W/m2"
               ,legpos = "topleft"
               ,plt    = TRUE)
phovdi27 = list(vnam   = c("nir.top","nir.beam.top","nir.diff.top"
                          ,"nir.bot","nir.beam.bot","nir.diff.bot")
               ,desc   = c("Top","Beam Top","Diffuse Top"
                          ,"Bottom","Beam Bottom","Diffuse Bottom")
               ,colour = c("firebrick","sienna","gold"
                          ,"midnightblue","steelblue","deepskyblue")
               ,lwd    = c(2.0,2.0,2.0,2.0,2.0,2.0)
               ,type   = ptype
               ,plog   = ""
               ,prefix = "nirrad"
               ,theme  = "NIR radiation"
               ,unit   = "W/m2"
               ,legpos = "topleft"
               ,plt    = TRUE)
phovdi28 = list(vnam   = c("swv.top","swv.beam.top","swv.diff.top"
                          ,"swv.bot","swv.beam.bot","swv.diff.bot")
               ,desc   = c("Top","Beam Top","Diffuse Top"
                          ,"Bottom","Beam Bottom","Diffuse Bottom")
               ,colour = c("firebrick","sienna","gold"
                          ,"midnightblue","steelblue","deepskyblue")
               ,lwd    = c(2.0,2.0,2.0,2.0,2.0,2.0)
               ,type   = ptype
               ,plog   = ""
               ,prefix = "swvrad"
               ,theme  = "SWV radiation"
               ,unit   = "W/m2"
               ,legpos = "topleft"
               ,plt    = TRUE)
#------------------------------------------------------------------------------------------#





#----- Define the PFT names. --------------------------------------------------------------#
pftnames = c("C4 Grass","Early Tropical","Mid Tropical","Late Tropical","Temp. C3 Grass"
             ,"North Pine","South Pine","Late Conifer","Early Temperate","Mid Temperate"
             ,"Late Temperate","C3 Pasture","C3 Crop","C4 Pasture","C4 Crop"
             ,"C3 Grass","Araucaria")
#------------------------------------------------------------------------------------------#



#----- Loading some packages. -------------------------------------------------------------#
library(hdf5)
library(chron)
library(scatterplot3d)
library(lattice)
library(maps)
library(mapdata)
library(akima)
#------------------------------------------------------------------------------------------#



#----- In case there is some graphic still opened. ----------------------------------------#
graphics.off()
#------------------------------------------------------------------------------------------#



#----- Setting how many formats we must output. -------------------------------------------#
outform = tolower(outform)
nout = length(outform)
#------------------------------------------------------------------------------------------#



#----- Avoiding unecessary and extremely annoying beeps. ----------------------------------#
options(locatorBell=FALSE)
#------------------------------------------------------------------------------------------#



#----- Loading some files with functions. -------------------------------------------------#
source(paste(srcdir,"atlas.r"      ,sep="/"))
source(paste(srcdir,"globdims.r"   ,sep="/"))
source(paste(srcdir,"locations.r"  ,sep="/"))
source(paste(srcdir,"muitas.r"     ,sep="/"))
source(paste(srcdir,"pretty.log.r" ,sep="/"))
source(paste(srcdir,"pretty.time.r",sep="/"))
source(paste(srcdir,"plotsize.r"   ,sep="/"))
source(paste(srcdir,"qapply.r"     ,sep="/"))
source(paste(srcdir,"rconstants.r" ,sep="/"))
source(paste(srcdir,"sombreado.r"  ,sep="/"))
source(paste(srcdir,"southammap.r" ,sep="/"))
source(paste(srcdir,"thermlib.r"   ,sep="/"))
source(paste(srcdir,"timeutils.r"  ,sep="/"))
#------------------------------------------------------------------------------------------#



#----- Defining plot window size ----------------------------------------------------------#
size = plotsize(proje=FALSE,paper=paper)
#------------------------------------------------------------------------------------------#



#----- Defining the initial and final time ------------------------------------------------#
whena = chron(dates=whena[1],times=whena[2])
whenz = chron(dates=whenz[1],times=whenz[2])

#----- Time series for the patch ----------------------------------------------------------#
phovdi      = list()
namesphovdi = NULL
for (s in 1:nphov){
  sss         = substring(100+s,2,3)
  hdhd        = paste("phovdi",sss,sep="")
  namesphovdi = c(namesphovdi,hdhd)
  phovdi[[s]] = get(hdhd)
} #end for
names(phovdi) = namesphovdi
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Big place loop starts here...                                                        #
#------------------------------------------------------------------------------------------#
for (place in myplaces){

   #----- Retrieve default information about this place and set up some variables. --------#
   thispoi    = locations(where=place,here=here)
   inpref     = paste(here,place,sep="/")
   outpref    = thispoi$pathout
   lieu       = thispoi$lieu
   suffix     = thispoi$iata


   #----- Print the banner to entretain the user. -----------------------------------------#
   print (paste("  + ",thispoi$lieu,"...",sep=""))

   #----- Make the main output directory in case it doesn't exist. ------------------------#
   if (! file.exists(outroot)) dir.create(outroot)
   outmain = paste(outroot,place,sep="/")
   if (! file.exists(outmain)) dir.create(outmain)
   outdir = paste(outmain,"runge-kutta",sep="/")
   if (! file.exists(outdir)) dir.create(outdir)


   #----- Determine the number of patches. ------------------------------------------------#
   filelist  = dir(inpref)
   mypatches = length(grep("thermo_state_prk4_patch_",filelist))



   #---------------------------------------------------------------------------------------#
   #    Patch loop.                                                                        #
   #---------------------------------------------------------------------------------------#
   for (ipa in 1:mypatches){
      #----- Find the character version of the patch number. ------------------------------#
      cipa = substring(10000+ipa,2,5)

      print (paste("    - Patch # ",ipa,"...",sep=""))
      #----- Define the output directory. -------------------------------------------------#
      patchdir  = paste(outdir,paste("patch_",cipa,sep=""),sep="/")
      if (! file.exists(patchdir)) dir.create(patchdir)

      #----- Define the input file name. --------------------------------------------------#
      inputfile = paste(inpref,paste("thermo_state_prk4_patch_",cipa,".txt",sep=""),sep="/")
      print(paste("      * Open file:",inputfile))

      #----- Read the file, just to grab the header. --------------------------------------#
      vnames   = scan(file=inputfile,what="raw",nlines=1,quiet=TRUE)
      nvars    = length(vnames)
      for (v in 1:nvars){
         aux          = tolower(vnames[v])
         saux         = strsplit(aux,split="")[[1]]
         uscore       = which(saux == "_")
         saux[uscore] = "."
         vnames[v]    = paste(saux,collapse="")
      }#end for
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Read the input file, this time reading all the data and skipping the first     #
      # line.                                                                              #
      #------------------------------------------------------------------------------------#
      aux                 = as.numeric(scan(file=inputfile,what="numeric",skip=1
                                           ,quiet=TRUE))
      cpatch              = matrix(aux,ncol=nvars,byrow=TRUE)
      dimnames(cpatch)    = list(NULL,vnames)
      cpatch              = data.frame(cpatch)


      #----- Cumulative number of time steps. ---------------------------------------------#
      print(paste("      * Determine the number of steps..."))
      cpatch$cum.step = seq(from=1,to=length(cpatch$month),by=1)

      #----- Reduce the size of the file to be the period of interest only. ---------------#
      print(paste("      * Reduce data to the period of interest..."))
      when   = chron( chron(dates=paste(cpatch$month,cpatch$day,cpatch$year,sep="/"))
                    + cpatch$time/day.sec, out.format=c(dates="m/d/y",times="h:m:s"))
      sel    = when >= whena & when <= whenz
      cpatch = cpatch[sel,]
      when   = when[sel]

      #----- Re-scale or re-define some variables. ----------------------------------------#
      print(paste("      * Re-scale/re-define some variables..."))


      #----- Ideal gas departure. ---------------------------------------------------------#
      cpatch$atm.idgas = 100. * ( cpatch$atm.prss
                                / ( cpatch$atm.rhos * rdry * cpatch$atm.temp 
                                  * (1. + epim1 *cpatch$atm.shv))
                                - 1.0)
      cpatch$can.idgas = 100. * ( cpatch$can.prss
                                / ( cpatch$can.rhos * rdry * cpatch$can.temp 
                                  * (1. + epim1 *cpatch$can.shv))
                                - 1.0)


      #----- Total heterotrophic respiration. ---------------------------------------------#
      cpatch$hetresp         = cpatch$cwdrh        + cpatch$soilrh

      #----- Water flux in kg/m2/day. -----------------------------------------------------#
      cpatch$wflxlc         =   cpatch$wflxlc       * day.sec
      cpatch$wflxwc         =   cpatch$wflxwc       * day.sec
      cpatch$wflxgc         =   cpatch$wflxgc       * day.sec
      cpatch$wflxac         =   cpatch$wflxac       * day.sec
      cpatch$wflxca         = - cpatch$wflxac
      cpatch$transp         =   cpatch$transp       * day.sec

      #----- Canopy -> Atmosphere fluxes in W/m2. -----------------------------------------#
      cpatch$qwflxca         = - cpatch$wflxac      * alvli(cpatch$can.temp) / day.sec
      cpatch$hflxca          = - cpatch$hflxac
      cpatch$cflxca          = - cpatch$cflxac

      #----- Virtual potential temperature. -----------------------------------------------#
      cpatch$atm.thetav      = cpatch$atm.theta    * (1. + epim1 * cpatch$atm.shv)
      cpatch$can.thetav      = cpatch$can.theta    * (1. + epim1 * cpatch$can.shv)

      #----- Temperatures in Celsius. -----------------------------------------------------#
      cpatch$atm.temp        = cpatch$atm.temp     - t00
      cpatch$can.temp        = cpatch$can.temp     - t00
      cpatch$leaf.temp       = cpatch$leaf.temp    - t00
      cpatch$wood.temp       = cpatch$wood.temp    - t00
      cpatch$sfc.temp        = cpatch$sfc.temp     - t00

      #----- Specific humidity in g/kg. ---------------------------------------------------#
      cpatch$sfc.shv         = cpatch$sfc.shv      * 1000.
      cpatch$atm.shv         = cpatch$atm.shv      * 1000.
      cpatch$can.shv         = cpatch$can.shv      * 1000.

      #----- Relative humidity in %. ------------------------------------------------------#
      cpatch$atm.relhum      = cpatch$atm.relhum   * 100.
      cpatch$can.relhum      = cpatch$can.relhum   * 100.

      #----- Super-saturation in %, forced to be zero if it is not saturated. -------------#
      cpatch$atm.supsat      = cpatch$atm.relhum   - 100.
      cpatch$can.supsat      = cpatch$can.relhum   - 100.
      sel                    = cpatch$atm.supsat < 0.
      cpatch$atm.supsat[sel] = 0.
      sel                    = cpatch$can.supsat < 0.
      cpatch$can.supsat[sel] = 0.

      #----- Pressure in hPa. -------------------------------------------------------------#
      cpatch$atm.prss        = cpatch$atm.prss     * 0.01
      cpatch$can.prss        = cpatch$can.prss     * 0.01

      #----- Density in kg/m3. ------------------------------------------------------------#
      cpatch$atm.rhos        = cpatch$atm.rhos
      cpatch$can.rhos        = cpatch$can.rhos

      #----- Precipitation-related fluxes in mm/hr. ---------------------------------------#
      cpatch$prec            = cpatch$atm.prate   * 3600.
      cpatch$intercept       = cpatch$intercept   * 3600.
      cpatch$throughfall     = cpatch$throughfall * 3600.
      cpatch$wshed           = cpatch$wshed       * 3600.
      #------------------------------------------------------------------------------------#


      #----- Original ground -> canopy conductance. ---------------------------------------#
      cpatch$ggold           = 0.2 * cpatch$ustar
      #------------------------------------------------------------------------------------#


      #----- Ground -> canopy resistance. -------------------------------------------------#
      cpatch$rgbare = 1. / cpatch$ggbare
      cpatch$rgveg  = 1. / cpatch$ggveg
      cpatch$rgnet  = 1. / cpatch$ggnet
      cpatch$rgold  = 1. / cpatch$ggold
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Total shortwave radiation.                                                     #
      #------------------------------------------------------------------------------------#
      cpatch$par.top      = cpatch$par.beam.top + cpatch$par.diff.top
      cpatch$par.bot      = cpatch$par.beam.bot + cpatch$par.diff.bot

      cpatch$nir.top      = cpatch$nir.beam.top + cpatch$nir.diff.top
      cpatch$nir.bot      = cpatch$nir.beam.bot + cpatch$nir.diff.bot

      cpatch$swv.top      = cpatch$par.top      + cpatch$nir.top
      cpatch$swv.beam.top = cpatch$par.beam.top + cpatch$nir.beam.top
      cpatch$swv.diff.top = cpatch$par.diff.top + cpatch$nir.diff.top
      cpatch$swv.bot      = cpatch$par.bot      + cpatch$nir.bot
      cpatch$swv.beam.bot = cpatch$par.beam.bot + cpatch$nir.beam.bot
      cpatch$swv.diff.bot = cpatch$par.diff.bot + cpatch$nir.diff.bot
      #------------------------------------------------------------------------------------#

      #------------------------------------------------------------------------------------#
      #      Define a nice grid for time.                                                  #
      #------------------------------------------------------------------------------------#
      whenout = pretty.time(when,n=8)
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #   Plot the time series diagrams showing months and years.                          #
      #------------------------------------------------------------------------------------#
      print(paste("      * Plot some patch-level figures..."))
      for (hh in 1:nphov){

         #----- Retrieve variable information from the list. ------------------------------#
         phovdinow    = phovdi[[hh]]
         vnames       = phovdinow$vnam  
         description  = phovdinow$desc  
         lcolours     = phovdinow$colour
         llwd         = phovdinow$lwd
         ltype        = phovdinow$type
         plog         = phovdinow$plog
         prefix       = phovdinow$prefix
         theme        = phovdinow$theme 
         unit         = phovdinow$unit  
         legpos       = phovdinow$legpos
         plotit       = phovdinow$plt   
    
         if (plotit){


            #----- Define the number of layers. -------------------------------------------#
            nlayers   = length(vnames)
            ylimit    = range(cpatch[vnames],na.rm=TRUE)
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

            #------------------------------------------------------------------------------#
            #     Check if the directory exists.  If not, create it.                       #
            #------------------------------------------------------------------------------#
            print (paste("        > ",theme," time series ...",sep=""))

            #----- Loop over formats. -----------------------------------------------------#
            for (o in 1:nout){
               fichier = paste(patchdir,"/",prefix,"-patch-",cipa,"-",suffix
                              ,".",outform[o],sep="")
               if(outform[o] == "x11"){
                  X11(width=size$width,height=size$height,pointsize=ptsz)
               }else if(outform[o] == "png"){
                  png(filename=fichier,width=size$width*depth,height=size$height*depth
                     ,pointsize=ptsz,res=depth)
               }else if(outform[o] == "eps"){
                  postscript(file=fichier,width=size$width,height=size$height
                            ,pointsize=ptsz,paper=paper)
               }#end if

               letitre = paste(theme," - ",thispoi$lieu,"(Patch ",ipa,")",
                               " \n"," Time series: ",theme,sep="")

               plot(x=when,y=cpatch[[vnames[1]]],type="n",main=letitre,xlab="Time"
                   ,ylim=ylimit,ylab=paste("[",unit,"]",sep=""),log=plog,xaxt="n"
                   ,cex.main=cex.main)
               axis(side=1,at=whenout$levels,labels=whenout$labels,padj=whenout$padj)
               if (hovgrid){
                   abline(h=axTicks(side=2),v=whenout$levels,col="gray66",lty="dotted")
               }#end if
               for (l in 1:nlayers){
                  points(x=when,y=cpatch[[vnames[l]]],col=lcolours[l]
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
         }#end if plotit
      }#end for nphov
      #------------------------------------------------------------------------------------#

   }#end for (ipa in patches)
   #---------------------------------------------------------------------------------------#

}#end for (place in myplaces)
#------------------------------------------------------------------------------------------#
