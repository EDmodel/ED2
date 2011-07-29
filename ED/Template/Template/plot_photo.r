#----- Here is the user-defined variable section. -----------------------------------------#
here           = "thispath"                                  # Current directory.
srcdir         = "/n/Moorcroft_Lab/Users/mlongo/util/Rsc" # Source  directory.
outroot        = "thisoutroot" # Source  directory.
myplaces       = c("thispoly")
iphoto         = myphysiol
iallom         = myallom

#------------------------------------------------------------------------------------------#
#     Initial and final times, they must be character vectors of size 2, the first one     #
# with m/d/y, and the second one with h:m:s".                                              #
#------------------------------------------------------------------------------------------#
whena          = c("thismontha/thisdatea/thisyeara","thishoura:thisminua:00")
whenz          = c("thismonthz/thisdatez/thisyearz","thishourz:thisminuz:00")
ptype          = "l"                  # Type of plot
ptyped         = "p"                  # Type of plot

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

lai.min        = 1.e-5          # Minimum leaf area index
nplant.min     = 4.e-8          # Minimum nplant.

#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#     List of possible plots. In case you don't want some of them, simply switch plt to F. #
#------------------------------------------------------------------------------------------#
#----- Time series plots. -----------------------------------------------------------------#
nhov = 17
hovdi01 = list(vnam   = c("gpp","leaf.resp","a.open","a.clos")
              ,desc   = c("GPP","Leaf resp.","A (open)","A (closed)")
              ,colour = c("forestgreen","sienna","chartreuse","goldenrod")
              ,lwd    = c(1.5,1.5,1.5,1.5)
              ,type   = ptype
              ,prefix = "carbflux"
              ,theme  = "Carbon fluxes"
              ,unit   = "umol/m2_leaf/s"
              ,legpos = "topleft"
              ,plt    = TRUE)
hovdi02 = list(vnam   = c("parv","util.parv","parv.min")
              ,desc   = c("Cohort PAR","Used PAR","Daytime PAR")
              ,colour = c("darkorange","goldenrod","firebrick")
              ,lwd    = c(1.5,1.5)
              ,type   = ptype
              ,prefix = "parflux"
              ,theme  = "PAR fluxes"
              ,unit   = "umol/m2/s"
              ,legpos = "topleft"
              ,plt    = TRUE)
hovdi03 = list(vnam   = c("atm.temp","can.temp","leaf.temp","wood.temp","ground.temp")
              ,desc   = c("Atmosphere","Canopy air","Leaf","Wood","Surface")
              ,colour = c("deepskyblue","gray21","chartreuse","goldenrod","sienna")
              ,lwd    = c(1.5,1.5,1.5,1.5,1.5)
              ,type   = ptype
              ,prefix = "temperature"
              ,theme  = "Temperature"
              ,unit   = "degC"
              ,legpos = "topleft"
              ,plt    = TRUE)
hovdi04 = list(vnam   = c("atm.shv","can.shv","ground.shv","lsfc.shv.clos","lint.shv")
              ,desc   = c("Atmosphere","Canopy air","Ground","Leaf Surface","Intercellular")
              ,colour = c("deepskyblue","gray21","sienna","chartreuse","forestgreen")
              ,lwd    = c(1.5,1.5,1.5,1.5,1.5)
              ,type   = ptype
              ,prefix = "h2oclosed"
              ,theme  = "Water vapour mixing ratio (Closed stomata)"
              ,unit   = "g/kg"
              ,legpos = "topleft"
              ,plt    = TRUE)
hovdi05 = list(vnam   = c("atm.shv","can.shv","ground.shv","lsfc.shv.open","lint.shv")
              ,desc   = c("Atmosphere","Canopy air","Ground","Leaf Surface","Intercellular")
              ,colour = c("deepskyblue","gray21","sienna","chartreuse","forestgreen")
              ,lwd    = c(1.5,1.5,1.5,1.5,1.5)
              ,type   = ptype
              ,prefix = "h2oopen"
              ,theme  = "Water vapour mixing ratio (Open stomata)"
              ,unit   = "g/kg"
              ,legpos = "topleft"
              ,plt    = TRUE)
hovdi06 = list(vnam   = c("atm.co2","can.co2","lsfc.co2.clos","lint.co2.clos","compp")
              ,desc   = c("Atmosphere","Canopy air","Leaf surface","Intercellular"
                         ,"GPP comp. point")
              ,colour = c("deepskyblue","gray21","chartreuse","forestgreen","yellowgreen")
              ,lwd    = c(1.5,1.5,1.5,1.5,1.5)
              ,type   = ptype
              ,prefix = "co2closed"
              ,theme  = "CO2 mixing ratio (Closed stomata)"
              ,unit   = "umol/mol"
              ,legpos = "topleft"
              ,plt    = TRUE)
hovdi07 = list(vnam   = c("atm.co2","can.co2","lsfc.co2.open","lint.co2.open","compp")
              ,desc   = c("Atmosphere","Canopy air","Leaf surface","Intercellular"
                         ,"GPP comp. point")
              ,colour = c("deepskyblue","gray21","chartreuse","forestgreen","yellowgreen")
              ,lwd    = c(1.5,1.5,1.5,1.5)
              ,type   = ptype
              ,prefix = "co2open"
              ,theme  = "CO2 mixing ratio (open stomata)"
              ,unit   = "umol/mol"
              ,legpos = "topleft"
              ,plt    = TRUE)
hovdi08 = list(vnam   = c("prec")
              ,desc   = c("Precipitation")
              ,colour = c("midnightblue")
              ,lwd    = c(2.5)
              ,type   = ptype
              ,prefix = "prec"
              ,theme  = "Precipitation rate"
              ,unit   = "mm/hr"
              ,legpos = "topleft"
              ,plt    = TRUE)
hovdi09 = list(vnam   = c("wood.rbw","leaf.rbw","rsw","rsw.open","rsw.clos")
              ,desc   = c("Wood Bnd. layer","Leaf Bnd. layer","Stomatal"
                         ,"Stomatal (open)","Stomatal (closed)")
              ,colour = c("goldenrod","limegreen","midnightblue","steelblue","sienna")
              ,lwd    = c(1.5,1.5,1.5,1.5,1.5)
              ,type   = ptype
              ,prefix = "resist"
              ,theme  = "Resistance"
              ,unit   = "s/m"
              ,legpos = "topleft"
              ,plt    = TRUE)
hovdi10 = list(vnam   = c("fsw","fsn","fs.open")
              ,desc   = c("Water","Nitrogen","Net")
              ,colour = c("deepskyblue","goldenrod","midnightblue")
              ,lwd    = c(1.5,1.5,1.5)
              ,type   = ptype
              ,prefix = "fsopen"
              ,theme  = "Fraction of stomata that are open"
              ,unit   = "n/d"
              ,legpos = "topleft"
              ,plt    = TRUE)
hovdi11 = list(vnam   = c("psi.open","psi.clos","h2o.supply")
              ,desc   = c("Psi (open)","Psi (closed)","Supply")
              ,colour = c("midnightblue","steelblue","olivedrab")
              ,lwd    = c(1.5,1.5,1.5)
              ,type   = ptype
              ,prefix = "psi"
              ,theme  = "Water flux"
              ,unit   = "kg/m2_leaf/day"
              ,legpos = "topleft"
              ,plt    = TRUE)
hovdi12 = list(vnam   = c("atm.wind","veg.wind","ustar")
              ,desc   = c("Free atmosphere","Vegetation","Friction vel.")
              ,colour = c("deepskyblue","forestgreen","sienna")
              ,lwd    = c(1.5,1.5,1.5)
              ,type   = ptype
              ,prefix = "wind"
              ,theme  = "Wind speed"
              ,unit   = "m/s"
              ,legpos = "topleft"
              ,plt    = TRUE)
hovdi13 = list(vnam   = c("limit.flag")
              ,desc   = c("Photo limitation")
              ,colour = c("olivedrab")
              ,lwd    = c(1.5,1.5,1.5)
              ,type   = ptyped
              ,prefix = "flglim"
              ,theme  = "Limitation flag"
              ,unit   = "n/d"
              ,legpos = "topleft"
              ,plt    = TRUE)
hovdi14 = list(vnam   = c("light.lim","rubisco.lim","low.co2.lim")
              ,desc   = c("Light","Rubisco","Low CO2")
              ,colour = c("lawngreen","goldenrod","steelblue")
              ,lwd    = c(1.5,1.5,1.5)
              ,type   = ptype
              ,prefix = "limitation"
              ,theme  = "Photosynthesis limitations"
              ,unit   = "umol/m2/s"
              ,legpos = "topleft"
              ,plt    = TRUE)
hovdi15 = list(vnam   = c("vm")
              ,desc   = c("Maximum capacity")
              ,colour = c("olivedrab")
              ,lwd    = c(1.5)
              ,type   = ptype
              ,prefix = "vm"
              ,theme  = "Maximum Rubisco capacity"
              ,unit   = "umol/m2/s"
              ,legpos = "topleft"
              ,plt    = TRUE)
hovdi16 = list(vnam   = c("wood.gbw","leaf.gbw","gsw","gsw.open","gsw.clos")
              ,desc   = c("Wood Bnd. layer","Leaf Bnd. layer","Stomatal"
                         ,"Stomatal (open)","Stomatal (closed)")
              ,colour = c("goldenrod","limegreen","midnightblue","steelblue","sienna")
              ,lwd    = c(1.5,1.5,1.5,1.5,1.5)
              ,type   = ptype
              ,prefix = "condct"
              ,theme  = "Conductance"
              ,unit   = "mol/m2_leaf/s"
              ,legpos = "topleft"
              ,plt    = TRUE)
hovdi17 = list(vnam   = c("wood.gbw.mmos","leaf.gbw.mmos","gsw.clos.mmos")
              ,desc   = c("Wood Bnd. layer","Leaf Bnd. layer","Stomatal (closed)")
              ,colour = c("goldenrod","limegreen","sienna")
              ,lwd    = c(1.5,1.5,1.5)
              ,type   = ptype
              ,prefix = "condctmmos"
              ,theme  = "Conductance"
              ,unit   = "mm/s"
              ,legpos = "topleft"
              ,plt    = TRUE)


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
nout    = length(outform)
#------------------------------------------------------------------------------------------#



#----- Avoiding unecessary and extremely annoying beeps. ----------------------------------#
options(locatorBell=FALSE)
#------------------------------------------------------------------------------------------#



#----- Load some files with global constants. ---------------------------------------------#
source(paste(srcdir,"allometry.r",sep="/"))
source(paste(srcdir,"rconstants.r",sep="/"))
#----- Load some files with useful functions. ---------------------------------------------#
source(paste(srcdir,"arrhenius.r",sep="/"))
source(paste(srcdir,"atlas.r",sep="/"))
source(paste(srcdir,"collatz.r",sep="/"))
source(paste(srcdir,"globdims.r",sep="/"))
source(paste(srcdir,"locations.r",sep="/"))
source(paste(srcdir,"muitas.r",sep="/"))
source(paste(srcdir,"pft.coms.r",sep="/"))
source(paste(srcdir,"pretty.log.r",sep="/"))
source(paste(srcdir,"pretty.time.r",sep="/"))
source(paste(srcdir,"plotsize.r",sep="/"))
source(paste(srcdir,"qapply.r",sep="/"))
source(paste(srcdir,"sombreado.r",sep="/"))
source(paste(srcdir,"southammap.r",sep="/"))
source(paste(srcdir,"timeutils.r",sep="/"))
#------------------------------------------------------------------------------------------#



#----- Defining plot window size ----------------------------------------------------------#
size = plotsize(proje=FALSE,paper=paper)
#------------------------------------------------------------------------------------------#



#----- Define the PFT names. --------------------------------------------------------------#
pftnames = c("C4 Grass","Early Tropical","Mid Tropical","Late Tropical","Temp. C3 Grass"
             ,"North Pine","South Pine","Late Conifer","Early Temperate","Mid Temperate"
             ,"Late Temperate","C3 Pasture","C3 Crop","C4 Pasture","C4 Crop"
             ,"C3 Grass","Araucaria")
#------------------------------------------------------------------------------------------#


#----- Defining the initial and final time ------------------------------------------------#
whena = chron(dates=whena[1],times=whena[2])
whenz = chron(dates=whenz[1],times=whenz[2])

#----- Time series for the cohort. --------------------------------------------------------#
hovdi      = list()
nameshovdi = NULL
for (s in 1:nhov){
  sss        = substring(100+s,2,3)
  hdhd       = paste("hovdi",sss,sep="")
  nameshovdi = c(nameshovdi,hdhd)
  hovdi[[s]] = get(hdhd)
} #end for
names(hovdi) = nameshovdi
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

   if (! file.exists(outroot)) dir.create(outroot)
   outmain = paste(outroot,place,sep="/")
   if (! file.exists(outmain)) dir.create(outmain)
   outdir = paste(outmain,"photosynthesis",sep="/")
   if (! file.exists(outdir)) dir.create(outdir)

   filelist  = dir(inpref)
   allpatches = substring(filelist,19,22)
   allpatches = as.numeric(allpatches)
   allpatches[!is.finite(allpatches)] = NA
   mypatches = max(allpatches,na.rm=TRUE)
   #---------------------------------------------------------------------------------------#
   #    Patch loop.                                                                        #
   #---------------------------------------------------------------------------------------#
   for (ipa in 1:mypatches){
      cipa = substring(10000+ipa,2,5)

      print (paste("    - Patch # ",ipa,"...",sep=""))

      patchdir = paste(outdir,paste("patch_",cipa,sep=""),sep="/")
      if (! file.exists(patchdir)) dir.create(patchdir)

      filelist  = dir(inpref)
      mycohorts = length(grep(paste("photo_state_patch_",cipa,sep=""),filelist))

      #------------------------------------------------------------------------------------#
      #    Cohort loop.                                                                    #
      #------------------------------------------------------------------------------------#
      for (ico in 1:mycohorts){
         #----- Find the character version of the cohort number. --------------------------#
         cico = substring(10000+ico,2,5)

         print (paste("      > Cohort # ",ico,"...",sep=""))
         #----- Define the output directory. ----------------------------------------------#
         cohortdir  = paste(patchdir,paste("cohort_",cico,sep=""),sep="/")
         if (! file.exists(cohortdir)) dir.create(cohortdir)

         #----- Define the input file name. -----------------------------------------------#
         inputfile = paste(inpref,paste("photo_state_patch_",cipa,"_cohort_",cico,".txt"
                                       ,sep="")
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

         #----- Reduce the size of the file to be the period of interest only. ------------#
         print(paste("        * Reduce data to the period of interest..."))
         when   = chron( chron(dates=paste(ccohort$month,ccohort$day,ccohort$year,sep="/"))
                       + ccohort$time/day.sec, out.format=c(dates="m/d/y",times="h:m:s"))
         sel    = when >= whena & when <= whenz
         ccohort = ccohort[sel,]
         when   = when[sel]

         #----- Save some variables that don't change over time. --------------------------#
         ipft   = ccohort$pft[1]
         height = signif(x=ccohort$height[1],digits=3)


         #----- Skip this cohort if its nplant and LAI are tiny. --------------------------#
         if (ccohort$lai >= lai.min && ccohort$nplant >= nplant.min){

            #----- Re-scale or re-define some variables. ----------------------------------#
            print(paste("        * Re-scale/re-define some variables..."))

            #----- PAR flux in W/m2. ------------------------------------------------------#
            ccohort$parv.min       = (ccohort$parv * 0. + 0.5 * Watts.2.Ein)


            #------------------------------------------------------------------------------#
            #    Limitations due to different constrains.                                  #
            #------------------------------------------------------------------------------#
            if (ipft == 1){
               ccohort$light.lim   = ccohort$util.parv 
               ccohort$rubisco.lim = ccohort$vm
               ccohort$low.co2.lim = ( klowco2 * ccohort$vm * ccohort$lint.co2.open
                                     / mol.2.umol)
            }else{
               if (iphoto %in% 0:1){
                  kco2  = ( mol.2.umol 
                          * arrhenius(ccohort$leaf.temp,kco2.ref.ibis,kco2.hor.ibis) )
                  ko2   = arrhenius(ccohort$leaf.temp,ko2.ref.ibis,ko2.hor.ibis)
               }else if (iphoto %in% 2:3){
                  kco2  = ( mol.2.umol 
                          * collatz(ccohort$leaf.temp,kco2.ref.coll,kco2.base.coll) )
                  ko2   = arrhenius(ccohort$leaf.temp,ko2.ref.coll ,ko2.base.coll)
               }#end if

               ccohort$light.lim   = ( ccohort$util.parv 
                                     * ( (ccohort$lint.co2.open - ccohort$compp)
                                       / (ccohort$lint.co2.open + 2. * ccohort$compp)) )

               ccohort$rubisco.lim = ccohort$vm * ( (ccohort$lint.co2.open - ccohort$compp)
                                                  / ( ccohort$lint.co2.open
                                                    + kco2*(1.+0.209/ko2)))
               ccohort$low.co2.lim = ccohort$lint.co2.open + NA
            }#end if


            #----- Net water conductance. -------------------------------------------------#
            ccohort$gsw             = ( ccohort$gsw.open * ccohort$fs.open 
                                      + ccohort$gsw.clos * (1. - ccohort$fs.open) )

            #----- Convert conductance to resistance [s/m] for plotting. ------------------#
            ccohort$leaf.rbw        = ( ccohort$can.rhos * ccohort$can.shv 
                                      / ccohort$leaf.gbw )
            ccohort$wood.rbw        = ( ccohort$can.rhos * ccohort$can.shv 
                                      / ccohort$wood.gbw )
            ccohort$rsw             = ccohort$can.rhos * ccohort$can.shv / ccohort$gsw
            ccohort$rsw.open        = ccohort$can.rhos * ccohort$can.shv / ccohort$gsw.open
            ccohort$rsw.clos        = ccohort$can.rhos * ccohort$can.shv / ccohort$gsw.clos
            ccohort$leaf.rbw[!is.finite(ccohort$leaf.rbw)] = 0.
            ccohort$wood.rbw[!is.finite(ccohort$wood.rbw)] = 0.
            ccohort$rsw[!is.finite(ccohort$rsw)] = 0.
            ccohort$rsw.open[!is.finite(ccohort$rsw.open)] = 0.
            ccohort$rsw.clos[!is.finite(ccohort$rsw.clos)] = 0.

            #----- Convert the conductances from kg_H2O/m2/s to mol/m2/s. -----------------#
            ccohort$leaf.gbw.mmos    = ( ccohort$leaf.gbw  * 1000.
                                       / (ccohort$can.rhos * ccohort$can.shv) )
            ccohort$wood.gbw.mmos    = ( ccohort$wood.gbw  * 1000.
                                       / (ccohort$can.rhos * ccohort$can.shv) )
            ccohort$gsw.mmos         = ( ccohort$gsw      * 1000.
                                       / (ccohort$can.rhos * ccohort$can.shv) )
            ccohort$gsw.open.mmos    = ( ccohort$gsw.open * 1000.
                                       / (ccohort$can.rhos * ccohort$can.shv) )
            ccohort$gsw.clos.mmos    = ( ccohort$gsw.clos * 1000.
                                       / (ccohort$can.rhos * ccohort$can.shv) )

            #----- Convert the conductances from kg_H2O/m2/s to mol/m2/s. -----------------#
            ccohort$leaf.gbw        = ccohort$leaf.gbw * mmh2oi
            ccohort$wood.gbw        = ccohort$wood.gbw * mmh2oi
            ccohort$gsw             = ccohort$gsw      * mmh2oi
            ccohort$gsw.open        = ccohort$gsw.open * mmh2oi
            ccohort$gsw.clos        = ccohort$gsw.clos * mmh2oi

            #----- Temperatures in Celsius. -----------------------------------------------#
            ccohort$atm.temp        = ccohort$atm.temp     - t00
            ccohort$can.temp        = ccohort$can.temp     - t00
            ccohort$leaf.temp       = ccohort$leaf.temp    - t00
            ccohort$wood.temp       = ccohort$wood.temp    - t00
            ccohort$ground.temp     = ccohort$ground.temp  - t00

            #----- Specific humidity in g/kg. ---------------------------------------------#
            ccohort$ground.shv      = ccohort$ground.shv    * 1000.
            ccohort$atm.shv         = ccohort$atm.shv       * 1000.
            ccohort$can.shv         = ccohort$can.shv       * 1000.
            ccohort$lsfc.shv.open   = ccohort$lsfc.shv.open * 1000.
            ccohort$lsfc.shv.clos   = ccohort$lsfc.shv.clos * 1000.
            ccohort$lint.shv        = ccohort$lint.shv      * 1000.

            #----- Pressure in hPa. -------------------------------------------------------#
            ccohort$atm.prss        = ccohort$atm.prss     * 0.01
            ccohort$can.prss        = ccohort$can.prss     * 0.01

            #----- Precipitation-related fluxes in mm/hr. ---------------------------------#
            ccohort$prec            = ccohort$pcpg         * 3600.

            #----- GPP and leaf respiration in µmol/m2leaf/s. -----------------------------#
            sel = ccohort$lai > 0.
            ccohort$gpp[sel]       = ccohort$gpp[sel]       / ccohort$lai[sel]
            ccohort$leaf.resp[sel] = ccohort$leaf.resp[sel] / ccohort$lai[sel]


            #----- Water demand terms in kg/m2leaf/day. -----------------------------------#
            sel = ccohort$lai > 0.
            ccohort$psi.open  [sel] = ccohort$psi.open  [sel] * day.sec / ccohort$lai[sel]
            ccohort$psi.clos  [sel] = ccohort$psi.clos  [sel] * day.sec / ccohort$lai[sel]
            ccohort$h2o.supply[sel] = ccohort$h2o.supply[sel] * day.sec / ccohort$lai[sel]


            #------------------------------------------------------------------------------#
            #      Define a nice grid for time.                                            #
            #------------------------------------------------------------------------------#
            whenout = pretty.time(when,n=8)
            #------------------------------------------------------------------------------#




            #------------------------------------------------------------------------------#
            #   Plot the time series diagrams showing months and years.                    #
            #------------------------------------------------------------------------------#
            print(paste("        * Plot some figures..."))
            for (hh in 1:nhov){

               #----- Retrieve variable information from the list. ------------------------#
               hovdinow    = hovdi[[hh]]
               vnames      = hovdinow$vnam  
               description = hovdinow$desc  
               lcolours    = hovdinow$colour
               llwd        = hovdinow$lwd
               ltype       = hovdinow$type
               prefix      = hovdinow$prefix
               theme       = hovdinow$theme 
               unit        = hovdinow$unit  
               legpos      = hovdinow$legpos
               plotit      = hovdinow$plt   
       
               if (plotit){


                  #----- Define the number of layers. -------------------------------------#
                  nlayers   = length(vnames)
                  ylimit    = range(ccohort[vnames],na.rm=TRUE)
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
                  print (paste("          # ",theme," time series ...",sep=""))

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
                                    " \n"," PFT: ",pftnames[ipft], " - Height: ",height,"m",
                                    sep="")

                     plot(x=when,y=ccohort[[vnames[1]]],type="n",main=letitre,xlab="Time"
                         ,ylim=ylimit,ylab=paste("[",unit,"]",sep=""),xaxt="n"
                         ,cex.main=cex.main)
                     axis(side=1,at=whenout$levels,labels=whenout$labels,padj=whenout$padj)
                     if (hovgrid){
                         abline(h=axTicks(side=2),v=whenout$levels,col="gray66"
                               ,lty="dotted")
                     }#end if
                     for (l in 1:nlayers){
                        points(x=when,y=ccohort[[vnames[l]]],col=lcolours[l]
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
         }#end if (lai > lai.min && nplant > nplant.min)


      }#end for (ico in cohorts)
      #------------------------------------------------------------------------------------#
   }#end for (ipa in mypatches)
}#end for (place in myplaces)
#------------------------------------------------------------------------------------------#
