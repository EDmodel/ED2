#----- Here is the user-defined variable section. -----------------------------------------#
here           = "thispath"                                  # Current directory.
srcdir         = "/n/Moorcroft_Lab/Users/mlongo/util/Rsc" # Source  directory.
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
nphov = 25
phovdi01 = list(vnam   = c("gpp","plresp","hetresp","cflxac")
               ,desc   = c("GPP","Plant resp.","Het. resp.","Atm->Canopy")
               ,colour = c("forestgreen","chartreuse","sienna","deepskyblue")
               ,lwd    = c(1.5,1.5,1.5,1.5)
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
               ,lwd    = c(1.5,1.5,1.5,1.5,1.5,1.5)
               ,type   = ptype
               ,plog   = ""
               ,prefix = "eneflux"
               ,theme  = "Energy fluxes"
               ,unit   = "W/m2"
               ,legpos = "topleft"
               ,plt    = TRUE)
phovdi03 = list(vnam   = c("wflxgc","wflxca","wflxlc","wflxwc","transp","dewgnd")
               ,desc   = c("Ground->Canopy","Canopy->Air","Leaf->Canopy","Wood->Canopy"
                          ,"Transpiration","Dew")
               ,colour = c("firebrick","midnightblue","chartreuse","goldenrod"
                          ,"darkolivegreen","deepskyblue")
               ,lwd    = c(1.5,1.5,1.5,1.5,1.5,1.5)
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
               ,lwd    = c(1.5,1.5,1.5,1.5)
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
               ,lwd    = c(1.5,1.5,1.5,1.5,1.5)
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
               ,lwd    = c(1.5,1.5,1.5)
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
               ,lwd    = c(1.5,1.5)
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
               ,lwd    = c(1.5,1.5)
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
               ,lwd    = c(1.5,1.5)
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
               ,lwd    = c(1.5,1.5)
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
               ,lwd    = c(1.5,1.5)
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
               ,lwd    = c(1.5,1.5)
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
               ,lwd    = c(2.5,1.5,1.5,1.5)
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
               ,lwd    = c(1.5,1.5)
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
               ,lwd    = c(1.5,1.5,1.5,1.5,1.5)
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
               ,lwd    = c(1.5,1.5)
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
               ,lwd    = c(1.5)
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
               ,lwd    = c(1.5)
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
               ,lwd    = c(1.5)
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
               ,lwd    = c(1.5,1.5)
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
               ,lwd    = c(1.5,1.5)
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
               ,lwd    = c(1.5)
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
               ,lwd    = c(1.5)
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
               ,lwd    = c(1.5,1.5,1.5)
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
               ,lwd    = c(1.5,1.5,1.5,1.5)
               ,type   = ptype
               ,plog   = ""
               ,prefix = "rgnd"
               ,theme  = "Ground to canopy resistance"
               ,unit   = "s/m"
               ,legpos = "topleft"
               ,plt    = TRUE)
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     List of possible cohort-level plots.  In case you don't want some of them, switch    #
# plt to F.                                                                                #
#------------------------------------------------------------------------------------------#
#----- Time series plots. -----------------------------------------------------------------#
nchov = 13
chovdi01 = list(vnam   = c("gpp","leaf.resp","root.resp","growth.resp","storage.resp"
                          ,"vleaf.resp")
               ,desc   = c("GPP","R (Leaf)","R (root)","R (growth)","R (storage)"
                          ,"R (vleaf)")
               ,colour = c("forestgreen","chartreuse","sienna","steelblue","goldenrod"
                          ,"purple3")
               ,cohlev = c(TRUE,TRUE,TRUE,TRUE,TRUE,TRUE)
               ,lwd    = c(1.5,1.5,1.5,1.5,1.5,1.5)
               ,type   = ptype
               ,prefix = "carbflux"
               ,theme  = "Carbon fluxes"
               ,unit   = "umol/m2/s"
               ,legpos = "topleft"
               ,plt    = TRUE)
chovdi02 = list(vnam   = c("atm.temp","can.temp","leaf.temp","wood.temp","sfc.temp")
               ,desc   = c("Atmosphere","Canopy air","Leaf","Wood","Surface")
               ,colour = c("deepskyblue","gray21","chartreuse","goldenrod","sienna")
               ,cohlev = c(FALSE,FALSE,TRUE,TRUE,FALSE)
               ,lwd    = c(1.5,1.5,1.5,1.5,1.5)
               ,type   = ptype
               ,prefix = "temperature"
               ,theme  = "Temperature"
               ,unit   = "degC"
               ,legpos = "topleft"
               ,plt    = TRUE)
chovdi03 = list(vnam   = c("atm.shv","can.shv","sfc.shv","lint.shv")
               ,desc   = c("Atmosphere","Canopy air","Ground","Intercellular")
               ,colour = c("deepskyblue","gray21","sienna","chartreuse","forestgreen")
               ,cohlev = c(FALSE,FALSE,FALSE,TRUE)
               ,lwd    = c(1.5,1.5,1.5,1.5)
               ,type   = ptype
               ,prefix = "h2ovap"
               ,theme  = "Water vapour mixing ratio"
               ,unit   = "g/kg"
               ,legpos = "topleft"
               ,plt    = TRUE)
chovdi04 = list(vnam   = c("wood.rbw","leaf.rbw","rsw","rsw.open","rsw.clos")
               ,desc   = c("Wood Bnd. layer","Leaf Bnd. layer","Stomatal"
                          ,"Stomatal (open)","Stomatal (closed)")
               ,colour = c("goldenrod","limegreen","midnightblue","steelblue","sienna")
               ,cohlev = c(TRUE,TRUE,TRUE,TRUE,TRUE)
               ,lwd    = c(1.5,1.5,1.5,1.5,1.5)
               ,type   = ptype
               ,prefix = "resist"
               ,theme  = "Resistance"
               ,unit   = "s/m"
               ,legpos = "topleft"
               ,plt    = TRUE)
chovdi05 = list(vnam   = c("fs.open")
               ,desc   = c("Net")
               ,colour = c("midnightblue")
               ,cohlev = c(TRUE)
               ,lwd    = c(1.5)
               ,type   = ptype
               ,prefix = "fsopen"
               ,theme  = "Fraction of stomata that are open"
               ,unit   = "n/d"
               ,legpos = "topleft"
               ,plt    = TRUE)
chovdi06 = list(vnam   = c("atm.vels","veg.wind","ustar")
               ,desc   = c("Free atmosphere","Vegetation","Friction vel.")
               ,colour = c("deepskyblue","forestgreen","sienna")
               ,cohlev = c(FALSE,TRUE,FALSE)
               ,lwd    = c(1.5,1.5,1.5)
               ,type   = ptype
               ,prefix = "wind"
               ,theme  = "Wind speed"
               ,unit   = "m/s"
               ,legpos = "topleft"
               ,plt    = TRUE)
chovdi07 = list(vnam   = c("wood.gbw","leaf.gbw","gsw","gsw.open","gsw.clos")
               ,desc   = c("Wood Bnd. layer","Leaf Bnd. layer","Stomatal"
                          ,"Stomatal (open)","Stomatal (closed)")
               ,colour = c("goldenrod","limegreen","midnightblue","steelblue","sienna")
               ,cohlev = c(TRUE,TRUE,TRUE,TRUE,TRUE)
               ,lwd    = c(1.5,1.5,1.5,1.5,1.5)
               ,type   = ptype
               ,prefix = "condct"
               ,theme  = "Conductance"
               ,unit   = "mol/m2_leaf/s"
               ,legpos = "topleft"
               ,plt    = TRUE)
chovdi08 = list(vnam   = c("wood.gbw.mmos","leaf.gbw.mmos","gsw.clos.mmos")
               ,desc   = c("Wood Bnd. layer","Leaf Bnd. layer","Stomatal (closed)")
               ,colour = c("goldenrod","limegreen","sienna")
               ,cohlev = c(TRUE,TRUE,TRUE)
               ,lwd    = c(1.5,1.5,1.5)
               ,type   = ptype
               ,prefix = "condctmmos"
               ,theme  = "Conductance"
               ,unit   = "mm/s"
               ,legpos = "topleft"
               ,plt    = TRUE)
chovdi09 = list(vnam   = c("leaf.liquid","leaf.frozen","wood.liquid","wood.frozen")
               ,desc   = c("Liquid (Leaf)","Ice (Leaf)","Liquid (Wood)","Ice (Wood)")
               ,colour = c("midnightblue","deepskyblue","sienna","goldenrod")
               ,cohlev = c(TRUE,TRUE,TRUE,TRUE)
               ,lwd    = c(1.5,1.5,1.5,1.5)
               ,type   = ptype
               ,prefix = "water"
               ,theme  = "Water on vegetation surface"
               ,unit   = "kg_h2o/m2"
               ,legpos = "topleft"
               ,plt    = TRUE)
chovdi10 = list(vnam   = c("leaf.energy","wood.energy")
               ,desc   = c("Leaf","Wood")
               ,colour = c("forestgreen","goldenrod")
               ,cohlev = c(TRUE,TRUE)
               ,lwd    = c(1.5,1.5)
               ,type   = ptype
               ,prefix = "energy"
               ,theme  = "Internal energy"
               ,unit   = "J/m2"
               ,legpos = "topleft"
               ,plt    = TRUE)
chovdi11 = list(vnam   = c("leaf.hcap","leaf.h2o.hcap","wood.hcap","wood.h2o.hcap")
               ,desc   = c("Leaf","Leaf Water","Wood","Wood Water")
               ,colour = c("lawngreen","deepskyblue","sienna","goldenrod")
               ,cohlev = c(TRUE,TRUE,TRUE,TRUE)
               ,lwd    = c(1.5,1.5,1.5,1.5)
               ,type   = ptype
               ,prefix = "hcap"
               ,theme  = "Heat capacity"
               ,unit   = "J/m2/K"
               ,legpos = "topleft"
               ,plt    = TRUE)
chovdi12 = list(vnam   = c("rshort.l","rlong.l","hflxlc","qwflxlc","qwshed","qtransp"
                          ,"qintercepted")
               ,desc   = c("Shortwave","Longwave","Sensible","Evaporation","Shedding"
                          ,"Transpiration","Intercepted")
               ,colour = c("yellow3","goldenrod","firebrick","midnightblue","royalblue"
                          ,"lawngreen","deepskyblue")
               ,cohlev = c(TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE)
               ,lwd    = c(1.5,1.5,1.5,1.5,1.5,1.5,1.5)
               ,type   = ptype
               ,prefix = "leafenergy"
               ,theme  = "Energy fluxes"
               ,unit   = "W/m2"
               ,legpos = "topleft"
               ,plt    = TRUE)
chovdi13 = list(vnam   = c("rshort.w","rlong.w","hflxwc","qwflxwc","qwshed","qtransp"
                          ,"qintercepted")
               ,desc   = c("Shortwave","Longwave","Sensible","Evaporation","Shedding"
                          ,"Transpiration","Intercepted")
               ,colour = c("yellow3","goldenrod","firebrick","midnightblue","royalblue"
                          ,"lawngreen","deepskyblue")
               ,cohlev = c(TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE)
               ,lwd    = c(1.5,1.5,1.5,1.5,1.5,1.5,1.5)
               ,type   = ptype
               ,prefix = "woodenergy"
               ,theme  = "Energy fluxes"
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
source(paste(srcdir,"atlas.r",sep="/"))
source(paste(srcdir,"globdims.r",sep="/"))
source(paste(srcdir,"locations.r",sep="/"))
source(paste(srcdir,"muitas.r",sep="/"))
source(paste(srcdir,"pretty.log.r",sep="/"))
source(paste(srcdir,"pretty.time.r",sep="/"))
source(paste(srcdir,"plotsize.r",sep="/"))
source(paste(srcdir,"qapply.r",sep="/"))
source(paste(srcdir,"rconstants.r",sep="/"))
source(paste(srcdir,"sombreado.r",sep="/"))
source(paste(srcdir,"southammap.r",sep="/"))
source(paste(srcdir,"timeutils.r",sep="/"))
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

#----- Time series for the cohort ---------------------------------------------------------#
chovdi      = list()
nameschovdi = NULL
for (s in 1:nchov){
  sss         = substring(100+s,2,3)
  hdhd        = paste("chovdi",sss,sep="")
  nameschovdi = c(nameschovdi,hdhd)
  chovdi[[s]] = get(hdhd)
} #end for
names(chovdi) = nameschovdi
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
      cpatch$wflxlc          =   cpatch$wflxlc       * day.sec
      cpatch$wflxwc          =   cpatch$wflxwc       * day.sec
      cpatch$wflxgc          =   cpatch$wflxgc       * day.sec
      cpatch$wflxac          =   cpatch$wflxac       * day.sec
      cpatch$wflxca          = - cpatch$wflxac
      cpatch$transp          =   cpatch$transp       * day.sec
      cpatch$dewgnd          =   cpatch$dewgnd       * day.sec

      #----- Canopy -> Atmosphere fluxes in W/m2. -----------------------------------------#
      cpatch$qwflxca         = -cpatch$wflxac      * alvl / day.sec
      cpatch$hflxca          = -cpatch$hflxac
      cpatch$cflxca          = -cpatch$cflxac

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

               plot(x=when,y=cpatch[[vnames[1]]],type="n",main=letitre,xlab="Time",
                    ylim=ylimit,ylab=paste("[",unit,"]",sep=""),log=plog,xaxt="n")
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


      #------------------------------------------------------------------------------------#
      #    Now we loop over the cohorts.                                                   #
      #------------------------------------------------------------------------------------#
      filelist  = dir(inpref)
      mycohorts = length(grep(paste("thermo_state_crk4_patch_",cipa,sep=""),filelist))

      #------------------------------------------------------------------------------------#
      #    Create an output table with some cohort properties.                             #
      #------------------------------------------------------------------------------------#
      cohtable        = data.frame(matrix(NA,nrow=mycohorts,ncol=9))
      names(cohtable) = c("ico","leaf.resolve","wood.resolve","pft","nplant","height"
                         ,"lai","wai","crown.area")

      #------------------------------------------------------------------------------------#
      #    Cohort loop.                                                                    #
      #------------------------------------------------------------------------------------#
      for (ico in 1:mycohorts){
         #----- Find the character version of the cohort number. --------------------------#
         cico = substring(10000+ico,2,5)

         print (paste("        > Cohort # ",ico,"...",sep=""))
         #----- Define the output directory. ----------------------------------------------#
         cohortdir  = paste(patchdir,paste("cohort_",cico,sep=""),sep="/")
         if (! file.exists(cohortdir)) dir.create(cohortdir)

         #----- Define the input file name. -----------------------------------------------#
         inputfile = paste(inpref,paste("thermo_state_crk4_patch_",cipa,"_cohort_",cico
                                       ,".txt",sep="")
                          ,sep="/")
         print(paste("        * Open file:",inputfile))

         #----- Read the file, just to grab the header. -----------------------------------#
         vnames   = scan(file=inputfile,what="raw",nlines=1,quiet=TRUE)
         nvars    = length(vnames)
         for (v in 1:nvars){
            aux          = tolower(vnames[v])
            saux         = strsplit(aux,split="")[[1]]
            uscore       = which(saux == "_")
            saux[uscore] = "."
            vnames[v]    = paste(saux,collapse="")
         }#end for
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #     Read the input file, this time reading all the data and skipping the first  #
         # line.                                                                           #
         #---------------------------------------------------------------------------------#
         aux                 = as.numeric(scan(file=inputfile,what="numeric",skip=1
                                              ,quiet=TRUE))
         ccohort              = matrix(aux,ncol=nvars,byrow=TRUE)
         dimnames(ccohort)    = list(NULL,vnames)
         ccohort              = data.frame(ccohort)

         #----- Fill in the table. --------------------------------------------------------#
         cohtable[ico,] = c(ico,ccohort$leaf.resolve[1],ccohort$wood.resolve[1]
                           ,ccohort$pft[1],ccohort$nplant[1],ccohort$height[1]
                           ,ccohort$lai[1],ccohort$wai[1],ccohort$crown.area[1])

         #----- Reduce the size of the file to be the period of interest only. ------------#
         print(paste("          * Reduce data to the period of interest..."))
         when    = chron( chron(dates=paste(ccohort$month,ccohort$day,ccohort$year,sep="/"))
                        + ccohort$time/day.sec, out.format=c(dates="m/d/y",times="h:m:s"))
         sel     = when >= whena & when <= whenz
         ccohort = ccohort[sel,]
         when    = when[sel]

         #----- Save some variables that don't change over time. --------------------------#
         pft    = ccohort$pft[1]
         height = signif(x=ccohort$height[1],digits=3)


         #----- Skip this cohort if its nplant and LAI are tiny. --------------------------#
         if (all(ccohort$leaf.resolve == 1 | ccohort$wood.resolve == 1)){


            #----- Net water conductance. -------------------------------------------------#
            ccohort$gsw             = ( ccohort$gsw.open * ccohort$fs.open 
                                      + ccohort$gsw.clos * (1. - ccohort$fs.open) )

            #----- Convert conductance to resistance [s/m] for plotting. ------------------#
            ccohort$leaf.rbw        = ( cpatch$can.rhos * 0.001 * cpatch$can.shv 
                                      / ccohort$leaf.gbw )
            ccohort$wood.rbw        = ( cpatch$can.rhos * 0.001 * cpatch$can.shv 
                                      / ccohort$wood.gbw )
            ccohort$rsw             = cpatch$can.rhos * 0.001 * cpatch$can.shv / ccohort$gsw
            ccohort$rsw.open        = cpatch$can.rhos * 0.001 * cpatch$can.shv / ccohort$gsw.open
            ccohort$rsw.clos        = cpatch$can.rhos * 0.001 * cpatch$can.shv / ccohort$gsw.clos
            ccohort$leaf.rbw[!is.finite(ccohort$leaf.rbw)] = 0.
            ccohort$wood.rbw[!is.finite(ccohort$wood.rbw)] = 0.
            ccohort$rsw[!is.finite(ccohort$rsw)] = 0.
            ccohort$rsw.open[!is.finite(ccohort$rsw.open)] = 0.
            ccohort$rsw.clos[!is.finite(ccohort$rsw.clos)] = 0.

            #----- Convert the conductances from kg_H2O/m2/s to mm/s. ---------------------#
            ccohort$leaf.gbw.mmos      = ( ccohort$leaf.gbw * 1000.
                                         / (cpatch$can.rhos * 0.001 * cpatch$can.shv) )
            ccohort$wood.gbw.mmos      = ( ccohort$wood.gbw * 1000.
                                         / (cpatch$can.rhos * 0.001 * cpatch$can.shv) )
            ccohort$gsw.mmos           = ( ccohort$gsw      * 1000.
                                         / (cpatch$can.rhos * 0.001 *  cpatch$can.shv) )
            ccohort$gsw.open.mmos      = ( ccohort$gsw.open * 1000.
                                         / (cpatch$can.rhos * 0.001 *  cpatch$can.shv) )
            ccohort$gsw.clos.mmos      = ( ccohort$gsw.clos * 1000.
                                         / (cpatch$can.rhos * 0.001 *  cpatch$can.shv) )

            #----- Convert the conductances from kg_H2O/m2/s to mol/m2/s. -----------------#
            ccohort$leaf.gbw        = ccohort$leaf.gbw * mmh2oi
            ccohort$wood.gbw        = ccohort$wood.gbw * mmh2oi
            ccohort$gsw             = ccohort$gsw      * mmh2oi
            ccohort$gsw.open        = ccohort$gsw.open * mmh2oi
            ccohort$gsw.clos        = ccohort$gsw.clos * mmh2oi

            #----- Temperatures in Celsius. -----------------------------------------------#
            ccohort$leaf.temp       = ccohort$leaf.temp     - t00
            ccohort$wood.temp       = ccohort$wood.temp     - t00

            #----- Specific humidity in g/kg. ---------------------------------------------#
            ccohort$lint.shv        = ccohort$lint.shv      * 1000.

            #----- Water sitting on leaf surface. -----------------------------------------#
            ccohort$leaf.liquid      = ccohort$leaf.water * ccohort$leaf.fliq
            ccohort$leaf.frozen      = ccohort$leaf.water * (1. - ccohort$leaf.fliq)
            ccohort$wood.liquid      = ccohort$wood.water * ccohort$wood.fliq
            ccohort$wood.frozen      = ccohort$wood.water * (1. - ccohort$wood.fliq)

            #----- Water heat capacity. ---------------------------------------------------#
            ccohort$leaf.h2o.hcap    = ( ccohort$leaf.liquid * cliq
                                       + ccohort$leaf.frozen * cice )
            ccohort$wood.h2o.hcap    = ( ccohort$wood.liquid * cliq
                                       + ccohort$wood.frozen * cice )


            #------------------------------------------------------------------------------#
            #      Define a nice grid for time.                                            #
            #------------------------------------------------------------------------------#
            whenout = pretty.time(when,n=8)
            #------------------------------------------------------------------------------#

            #------------------------------------------------------------------------------#
            #   Plot the time series diagrams showing months and years.                    #
            #------------------------------------------------------------------------------#
            print(paste("          * Plot some figures..."))
            for (hh in 1:nchov){

               #----- Retrieve variable information from the list. ------------------------#
               chovdinow   = chovdi[[hh]]
               vnames      = chovdinow$vnam  
               description = chovdinow$desc  
               lcolours    = chovdinow$colour
               lcohlev     = chovdinow$cohlev
               llwd        = chovdinow$lwd
               ltype       = chovdinow$type
               prefix      = chovdinow$prefix
               theme       = chovdinow$theme 
               unit        = chovdinow$unit  
               legpos      = chovdinow$legpos
               plotit      = chovdinow$plt   
       
               if (plotit){


                  #----- Define the number of layers. -------------------------------------#
                  nlayers   = length(vnames)
                  tmp       = NULL
                  for (ll in 1: nlayers){
                     if (lcohlev[ll]){
                        tmp = c(tmp,ccohort[vnames[ll]])
                     }else{
                        tmp = c(tmp,cpatch[vnames[ll]])
                     }#end if
                  }#end if
                  ylimit    = range(tmp,na.rm=TRUE)
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

                  #------------------------------------------------------------------------#
                  #     Check if the directory exists.  If not, create it.                 #
                  #------------------------------------------------------------------------#
                  print (paste("            # ",theme," time series ...",sep=""))

                  #----- Loop over formats. -----------------------------------------------#
                  for (o in 1:nout){
                     fichier = paste(cohortdir,"/",prefix,"-cohort-",cico,"-",suffix
                                    ,".",outform[o],sep="")
                     if(outform[o] == "x11"){
                        X11(width=size$width,height=size$height,pointsize=ptsz)
                     }else if(outform[o] == "png"){
                        png(filename=fichier,width=size$width*depth
                           ,height=size$height*depth,pointsize=ptsz,res=depth)
                     }else if(outform[o] == "eps"){
                        postscript(file=fichier,width=size$width,height=size$height
                                  ,pointsize=ptsz,paper=paper)
                     }#end if

                     letitre = paste(theme," - ",thispoi$lieu, 
                                    " \n"," Time series: ",theme,
                                    " \n"," Patch - ",ipa,"     Cohort - ",ico,
                                    " \n"," PFT: ",pftnames[pft], " - Height: ",height,"m",
                                    sep="")

                     if (lcohlev[1]){
                        thisvar = ccohort[[vnames[1]]]
                     }else{
                        thisvar = cpatch[[vnames[1]]]
                     }#end if

                     plot(x=when,y=thisvar,type="n",main=letitre,xlab="Time",
                          ylim=ylimit,ylab=paste("[",unit,"]",sep=""),xaxt="n",cex.main=0.8)
                     axis(side=1,at=whenout$levels,labels=whenout$labels,padj=whenout$padj)
                     if (hovgrid){
                         abline(h=axTicks(side=2),v=whenout$levels,col="gray66"
                               ,lty="dotted")
                     }#end if
                     for (l in 1:nlayers){

                        if (lcohlev[l]){
                           thisvar = ccohort[[vnames[l]]]
                        }else{
                           thisvar = cpatch[[vnames[l]]]
                        }#end if
                        points(x=when,y=thisvar,col=lcolours[l]
                              ,lwd=llwd[l],type=ltype,pch=16,cex=0.9)
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
            }#end for nhov
            #------------------------------------------------------------------------------#
         }else{
            #----- Re-scale or re-define some variables. ----------------------------------#
            print(paste("        * Skipping sparse cohort..."))
         }#end if (any(resolvable == 0))


      }#end for (ico in cohorts)
      #------------------------------------------------------------------------------------#


   }#end for (ipa in patches)
   #---------------------------------------------------------------------------------------#

}#end for (place in myplaces)
#------------------------------------------------------------------------------------------#
