#==========================================================================================#
#==========================================================================================#
#     Leave these commands at the beginning.  They will refresh the session.               #
#------------------------------------------------------------------------------------------#
rm(list=ls())
graphics.off()
options(warn=0)
gc()
#==========================================================================================#
#==========================================================================================#



#==========================================================================================#
#==========================================================================================#
#      Here is the user defined variable section.                                          #
#------------------------------------------------------------------------------------------#

#----- Paths. -----------------------------------------------------------------------------#
here           = "thispath"    # Current directory.
srcdir         = "thisrscpath" # Source  directory.
outroot        = "thisoutroot" # Directory for figures
#------------------------------------------------------------------------------------------#



#----- Name of the simulations. -----------------------------------------------------------#
myplaces       = c("thispoly")
budget_midfix  = "budget_state_patch" # Text that defines a budget file.
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Initial and final times, they must be character vectors of size 2, the first one     #
# with m/d/y, and the second one with h:m:s".                                              #
#------------------------------------------------------------------------------------------#
whena          = c("thismontha/thisdatea/thisyeara","thishoura:thisminua:00")
whenz          = c("thismonthz/thisdatez/thisyearz","thishourz:thisminuz:00")
ptype          = "l"                  # Type of plot
#------------------------------------------------------------------------------------------#



#----- Plot options. ----------------------------------------------------------------------#
outform        = thisoutform   # Formats for output file.  Supported formats are:
                               #   - "X11"    - for printing on screen
                               #   - "quartz" - for printing on Mac OS screen
                               #   - "eps"    - for postscript printing
                               #   - "png"    - for PNG printing
                               #   - "tif"    - for TIFF printing
                               #   - "pdf"    - for PDF printing
depth          = 96            # PNG resolution, in pixels per inch
ptsz           = 22            # Font size.
paper          = "square"      # Paper
lwidth         = 2.5           # Line width
plotgrid       = FALSE         # Should I plot the grid in the background? 
cex.main       = 0.8           # Scale coefficient for the title
ibackground    = mybackground  # Background settings (check load_everything.r)
f.leg          = 1/6           # Fraction of plotting area for legend.
same.lty       = TRUE          # Use the same lty for all curves (FALSE uses budget)
lty.all        = "solid"       # Line type for all curves in case same.lty = TRUE
lwd.all        = 3.0           # Line width for all curves in case same.lty = TRUE
density.terms  = 30            # Polygon line density for all terms (except for residual)
density.resid  = -1            # Polygon line density for for residuals
alpha.terms    = 0.70          # Polygon opacity for all terms (except for residual)
alpha.resid    = 1.00          # Polygon opacity for residual
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#  nyears.min.usemo  -- Minimum number of years to aggregate data to monthly.              #
#                       When time span exceeds this threshold, elapsed time will be shown  #
#                       in years.                                                          #
#  nmonths.min.usedy -- Minimum number of months to aggregate data to daily.               #
#                       When time span exceeds this threshold, elapsed time will be shown  #
#                       in months.                                                         #
#------------------------------------------------------------------------------------------#
nyears.min.usemo  = 6
nmonths.min.usedy = 4
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#     List of possible plots. In case you don't want some of them, simply switch plt to F. #
# Important, always leave the residual as the last variable to plot.                       #
#------------------------------------------------------------------------------------------#


#----- Time series plots. -----------------------------------------------------------------#
n           = 0
budget      = list()
n           = n + 1
budget[[n]] = list( vnam    = c("co2.nee.rel","co2.eddy.flux.rel","co2.veg.dyn.rel"
                               ,"co2.dstorage.rel","co2.residual.rel")
                  , desc    = c("NEE","Eddy flux","Capacity change","Storage change"
                               ,"Residual")
                  , colour  = c("#811F9E","#1BA2F7","#0E6E81","#CB003D","grey30")
                  , lty     = c("longdash","dashed","twodash","dotdash","solid")
                  , lwd     = c(3.0,3.0,3.0,3.0,3.0)
                  , range   = c(TRUE,TRUE,TRUE,FALSE,TRUE)
                  , type    = ptype
                  , prefix  = "co2.rel"
                  , theme   = "Accumulated Carbon dioxide (relative)"
                  , unit    = "pc"
                  , mult    = "1.0"
                  , ylim    = NA # c(-100,100)
                  , average = FALSE
                  , poly    = TRUE
                  , plt     = TRUE
                  )#end list
n           = n + 1
budget[[n]] = list( vnam    = c("crb.eddy.flux.rel","crb.veg.dyn.rel","crb.netseed.rel"
                               ,"crb.dstorage.rel","crb.residual.rel")
                  , desc    = c("Eddy flux","Capacity change","Net seed flux"
                               ,"Storage change","Residual")
                  , colour  = c("#1BA2F7","#0E6E81","#F87856","#CB003D","grey30")
                  , lty     = c("dashed","twodash","dotted","dotdash","solid")
                  , lwd     = c(3.0,3.0,4.0,3.0,3.0)
                  , range   = c(TRUE,TRUE,TRUE,FALSE,TRUE)
                  , type    = ptype
                  , prefix  = "crb.rel"
                  , theme   = "Accumulated Carbon (relative)"
                  , unit    = "pc"
                  , mult    = "1.0"
                  , ylim    = NA # c(-0.01,0.01)
                  , average = FALSE
                  , poly    = TRUE
                  , plt     = TRUE
                  )#end list
n           = n + 1
budget[[n]] = list( vnam    = c("ent.precip.rel","ent.netrad.rel","ent.prss.eff.rel"
                               ,"ent.eddy.flux.rel","ent.runoff.rel","ent.veg.dyn.rel"
                               ,"ent.dstorage.rel","ent.residual.rel")
                  , desc    = c("Rainfall","Net Radiation","Pressure effect","Eddy flux"
                               ,"Total runoff","Capacity change","Storage change","Residual")
                  , colour  = c("#2BD2DB","#F87856","#CCCA3D","#1BA2F7"
                               ,"#811F9E","#0E6E81","#CB003D","grey30")
                  , lty     = c("dotdash","dashed","longdash","dashed"
                               ,"dotted","twodash","dotdash","solid")
                  , lwd     = c(3.0,3.0,3.0,3.0,4.0,3.0,3.0,3.0)
                  , range   = c(TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,FALSE,TRUE)
                  , type    = ptype
                  , prefix  = "ent.rel"
                  , theme   = "Accumulated enthalpy (relative)"
                  , unit    = "pc"
                  , mult    = "1.0"
                  , ylim    = NA # c(-10.,10.)
                  , average = FALSE
                  , poly    = TRUE
                  , plt     = TRUE
                  )#end list
n           = n + 1
budget[[n]] = list( vnam    = c("h2o.precip.rel","h2o.eddy.flux.rel","h2o.runoff.rel"
                               ,"h2o.veg.dyn.rel","h2o.dstorage.rel","h2o.residual.rel")
                  , desc    = c("Rainfall","Eddy flux","Total runoff"
                               ,"Capacity change","Storage change","Residual")
                  , colour  = c("#2BD2DB","#1BA2F7","#811F9E","#0E6E81","#CB003D","grey30")
                  , lty     = c("dotdash","dashed","dotted","twodash","dotdash","solid")
                  , lwd     = c(3.0,3.0,4.0,3.0,3.0,3.0)
                  , range   = c(TRUE,TRUE,TRUE,TRUE,FALSE,TRUE)
                  , type    = ptype
                  , prefix  = "h2o.rel"
                  , theme   = "Accumulated water (relative)"
                  , unit    = "pc"
                  , mult    = "1.0"
                  , ylim    = NA # c(-10,10)
                  , average = FALSE
                  , poly    = TRUE
                  , plt     = TRUE
                  )#end list
n           = n + 1
budget[[n]] = list( vnam    = c("co2.nee.cum","co2.eddy.flux.cum","co2.veg.dyn.cum"
                               ,"co2.dstorage.cum","co2.residual.cum")
                  , desc    = c("NEE","Eddy flux","Capacity change","Storage change"
                               ,"Residual")
                  , colour  = c("#811F9E","#1BA2F7","#0E6E81","#CB003D","grey30")
                  , lty     = c("longdash","dashed","twodash","dotdash","solid")
                  , lwd     = c(3.0,3.0,3.0,3.0,3.0)
                  , range   = c(TRUE,TRUE,TRUE,FALSE,TRUE)
                  , type    = ptype
                  , prefix  = "co2.cum"
                  , theme   = "Accumulated Carbon dioxide"
                  , unit    = "molcom2"
                  , mult    = "1.e-6"
                  , ylim    = NA # c(-100,100)
                  , average = FALSE
                  , poly    = TRUE
                  , plt     = TRUE
                  )#end list
n           = n + 1
budget[[n]] = list( vnam    = c("crb.eddy.flux.cum","crb.veg.dyn.cum","crb.netseed.cum"
                               ,"crb.dstorage.cum","crb.residual.cum")
                  , desc    = c("Eddy flux","Capacity change","Net seed flux"
                               ,"Storage change","Residual")
                  , colour  = c("#1BA2F7","#0E6E81","#F87856","#CB003D","grey30")
                  , lty     = c("dashed","twodash","dotted","dotdash","solid")
                  , lwd     = c(3.0,3.0,4.0,3.0,3.0)
                  , range   = c(TRUE,TRUE,TRUE,FALSE,TRUE)
                  , type    = ptype
                  , prefix  = "crb.cum"
                  , theme   = "Accumulated Carbon"
                  , unit    = "kgcom2"
                  , mult    = "1.0"
                  , ylim    = NA # c(-0.01,0.01)
                  , average = FALSE
                  , poly    = TRUE
                  , plt     = TRUE
                  )#end list
n           = n + 1
budget[[n]] = list( vnam    = c("ent.precip.cum","ent.netrad.cum","ent.prss.eff.cum"
                               ,"ent.eddy.flux.cum","ent.runoff.cum","ent.veg.dyn.cum"
                               ,"ent.dstorage.cum","ent.residual.cum")
                  , desc    = c("Rainfall","Net Radiation","Pressure effect","Eddy flux"
                               ,"Total runoff","Capacity change","Storage change"
                               ,"Residual")
                  , colour  = c("#2BD2DB","#F87856","#CCCA3D","#1BA2F7"
                               ,"#811F9E","#0E6E81","#CB003D","grey30")
                  , lty     = c("dotdash","dashed","longdash","dashed"
                               ,"dotted","twodash","dotdash","solid")
                  , lwd     = c(3.0,3.0,3.0,3.0,4.0,3.0,3.0,3.0)
                  , range   = c(TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,FALSE,TRUE)
                  , type    = ptype
                  , prefix  = "ent.cum"
                  , theme   = "Accumulated enthalpy"
                  , mult    = "1.e-9"
                  , unit    = "Gjom2"
                  , ylim    = NA
                  , average = FALSE
                  , poly    = TRUE
                  , plt     = TRUE
                  )#end list
n           = n + 1
budget[[n]] = list( vnam    = c("h2o.precip.cum","h2o.eddy.flux.cum","h2o.runoff.cum"
                               ,"h2o.veg.dyn.cum","h2o.dstorage.cum","h2o.residual.cum")
                  , desc    = c("Rainfall","Eddy flux","Total runoff"
                               ,"Capacity change","Storage change","Residual")
                  , colour  = c("#2BD2DB","#1BA2F7","#811F9E","#0E6E81","#CB003D","grey30")
                  , lty     = c("dotdash","dashed","dotted","twodash","dotdash","solid")
                  , lwd     = c(3.0,3.0,4.0,3.0,3.0,3.0)
                  , range   = c(TRUE,TRUE,TRUE,TRUE,FALSE,TRUE)
                  , type    = ptype
                  , prefix  = "h2o.cum"
                  , theme   = "Accumulated water"
                  , unit    = "Mgwom2"
                  , mult    = "1.e-3"
                  , ylim    = NA
                  , average = FALSE
                  , poly    = TRUE
                  , plt     = TRUE
                  )#end list
n           = n + 1
budget[[n]] = list( vnam    = c("co2.nee","co2.eddy.flux","co2.veg.dyn"
                               ,"co2.dstorage","co2.residual")
                  , desc    = c("NEE","Eddy flux","Capacity change","Storage change"
                               ,"Residual")
                  , colour  = c("#811F9E","#1BA2F7","#0E6E81","#CB003D","grey30")
                  , lty     = c("longdash","dashed","twodash","dotdash","solid")
                  , lwd     = c(3.0,3.0,3.0,3.0,3.0)
                  , range   = c(TRUE,TRUE,TRUE,FALSE,TRUE)
                  , type    = ptype
                  , prefix  = "co2.avg"
                  , theme   = "Carbon dioxide budget"
                  , unit    = "umolcom2os"
                  , mult    = "1.0"
                  , ylim    = NA # c(-100,100)
                  , average = FALSE
                  , poly    = TRUE
                  , plt     = TRUE
                  )#end list
n           = n + 1
budget[[n]] = list( vnam    = c("crb.eddy.flux","crb.veg.dyn","crb.netseed"
                               ,"crb.dstorage","crb.residual")
                  , desc    = c("Eddy flux","Capacity change","Net seed flux"
                               ,"Storage change","Residual")
                  , colour  = c("#1BA2F7","#0E6E81","#F87856","#CB003D","grey30")
                  , lty     = c("dashed","twodash","dotted","dotdash","solid")
                  , lwd     = c(3.0,3.0,4.0,3.0,3.0)
                  , range   = c(TRUE,TRUE,TRUE,FALSE,TRUE)
                  , type    = ptype
                  , prefix  = "crb.avg"
                  , theme   = "Carbon budget"
                  , unit    = "kgcom2oday"
                  , mult    = "day.sec"
                  , ylim    = NA # c(-0.01,0.01)
                  , average = FALSE
                  , poly    = TRUE
                  , plt     = TRUE
                  )#end list
n           = n + 1
budget[[n]] = list( vnam    = c("ent.precip","ent.netrad","ent.prss.eff"
                               ,"ent.eddy.flux","ent.runoff","ent.veg.dyn"
                               ,"ent.dstorage","ent.residual")
                  , desc    = c("Rainfall","Net Radiation","Pressure effect","Eddy flux"
                               ,"Total runoff","Capacity change","Storage change"
                               ,"Residual")
                  , colour  = c("#2BD2DB","#F87856","#CCCA3D","#1BA2F7"
                               ,"#811F9E","#0E6E81","#CB003D","grey30")
                  , lty     = c("dotdash","dashed","longdash","dashed"
                               ,"dotted","twodash","dotdash","solid")
                  , lwd     = c(3.0,3.0,3.0,3.0,4.0,3.0,3.0,3.0)
                  , range   = c(TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,FALSE,TRUE)
                  , type    = ptype
                  , prefix  = "ent.avg"
                  , theme   = "Enthalpy budget"
                  , mult    = "1.0"
                  , unit    = "wom2"
                  , ylim    = NA
                  , average = TRUE
                  , poly    = FALSE
                  , plt     = TRUE
                  )#end list
n           = n + 1
budget[[n]] = list( vnam    = c("h2o.precip","h2o.eddy.flux","h2o.runoff"
                               ,"h2o.veg.dyn","h2o.dstorage","h2o.residual")
                  , desc    = c("Rainfall","Eddy flux","Total runoff"
                               ,"Capacity change","Storage change","Residual")
                  , colour  = c("#2BD2DB","#1BA2F7","#811F9E","#0E6E81","#CB003D","grey30")
                  , lty     = c("dotdash","dashed","dotted","twodash","dotdash","solid")
                  , lwd     = c(3.0,3.0,4.0,3.0,3.0,3.0)
                  , range   = c(TRUE,TRUE,TRUE,TRUE,FALSE,TRUE)
                  , type    = ptype
                  , prefix  = "h2o.avg"
                  , theme   = "Water budget"
                  , unit    = "kgwom2oday"
                  , mult    = "day.sec"
                  , ylim    = NA
                  , average = TRUE
                  , poly    = FALSE
                  , plt     = TRUE
                  )#end list
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



#----- Setting how many formats we must output. -------------------------------------------#
outform = tolower(outform)
nout = length(outform)
#------------------------------------------------------------------------------------------#



#----- Define plot window size ------------------------------------------------------------#
f.ext   = f.leg / (1. - f.leg)
ey.size = plotsize(proje=FALSE,paper=paper,extendfc="lat",extfactor=f.ext)
#------------------------------------------------------------------------------------------#



#----- Define the initial and final time --------------------------------------------------#
whena = chron(dates=whena[1],times=whena[2])
whenz = chron(dates=whenz[1],times=whenz[2])
#------------------------------------------------------------------------------------------#



#----- Time series for the patch ----------------------------------------------------------#
nbudget      = length(budget)
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Big place loop starts here...                                                        #
#------------------------------------------------------------------------------------------#
for (place in myplaces){

   #----- Retrieve default information about this place and set up some variables. --------#
   thispoi    = locations(where=place,here=here)
   inpref     = file.path(here,place,"analy")
   outpref    = thispoi$pathout
   lieu       = thispoi$lieu
   suffix     = thispoi$iata
   #---------------------------------------------------------------------------------------#


   #----- Print the banner to entretain the user. -----------------------------------------#
   cat0("  + ",thispoi$lieu,".")
   #---------------------------------------------------------------------------------------#

   #----- Make the main output directory in case it doesn't exist. ------------------------#
   if (! file.exists(outroot)) dir.create(outroot)
   outmain = file.path(outroot,place)
   if (! file.exists(outmain)) dir.create(outmain)
   outdir  = file.path(outmain,"budget")
   if (! file.exists(outdir)) dir.create(outdir)
   #---------------------------------------------------------------------------------------#


   #----- Determine the number of patches. ------------------------------------------------#
   filelist  = dir(inpref)
   mypatches = length(grep(budget_midfix,filelist))
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #    Patch loop.                                                                        #
   #---------------------------------------------------------------------------------------#
   for (ipa in sequence(mypatches)){
      #----- Find the character version of the patch number. ------------------------------#
      cipa = sprintf("%4.4i",ipa)

      cat0("    - Patch # ",ipa,".")
      #----- Define the output directory. -------------------------------------------------#
      patchdir  = file.path(outdir,paste0("patch_",cipa))
      if (! file.exists(patchdir)) dir.create(patchdir)
      #------------------------------------------------------------------------------------#


      #----- Define the input file name. --------------------------------------------------#
      inputfile = file.path(inpref,paste0(place,"_",budget_midfix,"_",cipa,".txt"))
      cat0("      * Open file: ",basename(inputfile))
      #------------------------------------------------------------------------------------#


      #----- Read the file, just to grab the header. --------------------------------------#
      vnames   = scan(file=inputfile,what="raw",nlines=1,quiet=TRUE)
      nvars    = length(vnames)
      for (v in sequence(nvars)){
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
      #------------------------------------------------------------------------------------#


      #----- Reduce the size of the file to be the period of interest only. ---------------#
      cat0("      * Reduce data to the period of interest.")
      when   = chron( chron(dates=paste(cpatch$month,cpatch$day,cpatch$year,sep="/"))
                    + cpatch$time/day.sec, out.format=c(dates="m/d/y",times="h:m:s"))
      sel    = when >= whena & when <= whenz
      cpatch = cpatch[sel,,drop=FALSE]
      when   = when[sel]
      nwhen  = length(when)
      dtime  = mean((when[-1] - when[-nwhen])* day.sec)
      #------------------------------------------------------------------------------------#



      #----- Find additional time-related information. ------------------------------------#
      w.day     = numdays  (when)
      w.month   = nummonths(when)
      w.year    = numyears (when)
      w.today   = chron(paste(w.month,w.day,w.year,sep="/"))
      w.tomonth = chron(paste(w.month,    1,w.year,sep="/"))
      e.hours   = as.numeric(when-when[1]) * 24.
      e.days    = as.numeric(w.today-w.today[1])
      e.months  = 12 * (w.year - w.year[1]) + w.month - w.month[1]
      e.years   = floor(e.months/12)
      nmonths   = max(e.months)
      nyears    = max(e.years)
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #      Define time since beginning of the simulation.                                #
      #------------------------------------------------------------------------------------#
      if (nyears >= nyears.min.usemo){
         #----- Aggregate data by month. --------------------------------------------------#
         e.aggr   = e.months
         e.uniq   = unique(e.months)
         e.at     = pretty.elapsed(x=e.uniq,base=12)
         e.diff   = mean(diff(e.at))
         #---------------------------------------------------------------------------------#


         #----- Decide whether to show results in months or in years. ---------------------#
         if (e.diff >= 12){
            e.mult = 1. / 12
            e.unit = untab$yr
         }else{
            e.mult = 1.
            e.unit = untab$month
         }#end if (e.diff >= 12)
         #---------------------------------------------------------------------------------#
      }else if (nmonths >= nmonths.min.usedy){
         #------ Aggregate data by day. ---------------------------------------------------#
         e.aggr   = e.days
         e.uniq   = unique(e.days)
         e.at     = pretty.elapsed(x=e.uniq,base=30)
         e.diff   = mean(diff(e.at))
         #---------------------------------------------------------------------------------#


         #----- Decide whether to show results in days or in months. ----------------------#
         if (e.diff >= 30){
            e.mult = 1. / 30.
            e.unit = untab$month
         }else{
            e.mult = 1.
            e.unit = untab$day
         }#end if (e.diff >= 30)
         #---------------------------------------------------------------------------------#
      }else{
         #------ Aggregate data by day. ---------------------------------------------------#
         e.aggr   = e.hours
         e.uniq   = unique(e.hours)
         e.at     = pretty.elapsed(x=e.hours,base=24)
         e.diff   = mean(diff(e.at))
         #---------------------------------------------------------------------------------#


         #----- Decide whether to show results in days or in months. ----------------------#
         if (e.diff >= 24){
            e.mult = 1. / 24.
            e.unit = untab$day
         }else{
            e.mult = 1.
            e.unit = untab$hr
         }#end if (e.diff >= 24)
         #---------------------------------------------------------------------------------#
      }#end if
      #------------------------------------------------------------------------------------#



      #----- Additional plot annotation settings. -----------------------------------------#
      e.limit  = pretty.xylim(e.uniq)
      e.labels = sprintf("%g",e.at*e.mult)
      e.desc   = "Simulation time"
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Reduce the number of variables, and standardise signs:                         #
      # > 0 : entering the system;                                                         #
      # < 0 : leaving the system.                                                          #
      #------------------------------------------------------------------------------------#
      cat0("      * Reduce the number of terms.")
      #----- CO2. -------------------------------------------------------------------------#
      cpatch$co2.eddy.flux = - cpatch$co2.loss2atm
      cpatch$co2.nee       = - cpatch$co2.nep
      cpatch$co2.veg.dyn   =   cpatch$co2.zcan.eff
      cpatch$co2.residual  =   cpatch$co2.residual + cpatch$co2.dens.eff
      #----- Carbon. ----------------------------------------------------------------------#
      cpatch$crb.eddy.flux = - cpatch$crb.loss2atm
      cpatch$crb.veg.dyn   =   cpatch$crb.zcan.eff
      cpatch$crb.netseed   =   cpatch$crb.seedrain - cpatch$crb.loss2yield
      cpatch$crb.residual  =   cpatch$crb.residual + cpatch$crb.dens.eff
      #----- Enthalpy. --------------------------------------------------------------------#
      cpatch$ent.eddy.flux = - cpatch$ent.loss2atm
      cpatch$ent.runoff    = - cpatch$ent.runoff   - cpatch$ent.drainage
      cpatch$ent.veg.dyn   =   cpatch$ent.zcan.eff + cpatch$ent.hcap.eff
      cpatch$ent.residual  =   cpatch$ent.residual + cpatch$ent.dens.eff
      #----- Water. -----------------------------------------------------------------------#
      cpatch$h2o.eddy.flux = - cpatch$h2o.loss2atm
      cpatch$h2o.runoff    = - cpatch$h2o.runoff   - cpatch$h2o.drainage
      cpatch$h2o.veg.dyn   =   cpatch$h2o.zcan.eff
      cpatch$h2o.residual  =   cpatch$h2o.residual + cpatch$h2o.dens.eff
      #------------------------------------------------------------------------------------#


      #----- Re-scale or re-define some variables. ----------------------------------------#
      cat0("      * Find the average storage.")
      co2.scale = mean(cpatch$co2.storage)
      crb.scale = mean(cpatch$crb.storage)
      ent.scale = mean(cpatch$ent.storage)
      h2o.scale = mean(cpatch$h2o.storage)
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #      Re-scale or re-define some variables.                                         #
      #------------------------------------------------------------------------------------#
      cat0("      * Define the cumulative sum of all budget terms.")
      #----- Cumulative terms for all CO2 budget. -----------------------------------------#
      cpatch$co2.dstorage.cum  = cumsum( cpatch$co2.dstorage  * dtime )
      cpatch$co2.nee.cum       = cumsum( cpatch$co2.nee       * dtime )
      cpatch$co2.eddy.flux.cum = cumsum( cpatch$co2.eddy.flux * dtime )
      cpatch$co2.veg.dyn.cum   = cumsum( cpatch$co2.veg.dyn   * dtime )
      cpatch$co2.residual.cum  = cumsum( cpatch$co2.residual  * dtime )
      #----- Cumulative terms for all carbon budget. --------------------------------------#
      cpatch$crb.dstorage.cum  = cumsum( cpatch$crb.dstorage  * dtime )
      cpatch$crb.eddy.flux.cum = cumsum( cpatch$crb.eddy.flux * dtime )
      cpatch$crb.veg.dyn.cum   = cumsum( cpatch$crb.veg.dyn   * dtime )
      cpatch$crb.netseed.cum   = cumsum( cpatch$crb.netseed   * dtime )
      cpatch$crb.residual.cum  = cumsum( cpatch$crb.residual  * dtime )
      #----- Cumulative terms for all enthalpy budget. ------------------------------------#
      cpatch$ent.dstorage.cum  = cumsum( cpatch$ent.dstorage  * dtime )
      cpatch$ent.precip.cum    = cumsum( cpatch$ent.precip    * dtime )
      cpatch$ent.netrad.cum    = cumsum( cpatch$ent.netrad    * dtime )
      cpatch$ent.prss.eff.cum  = cumsum( cpatch$ent.prss.eff  * dtime )
      cpatch$ent.eddy.flux.cum = cumsum( cpatch$ent.eddy.flux * dtime )
      cpatch$ent.runoff.cum    = cumsum( cpatch$ent.runoff    * dtime )
      cpatch$ent.veg.dyn.cum   = cumsum( cpatch$ent.veg.dyn   * dtime )
      cpatch$ent.residual.cum  = cumsum( cpatch$ent.residual  * dtime )
      #----- Cumulative terms for all water budget. ---------------------------------------#
      cpatch$h2o.dstorage.cum  = cumsum( cpatch$h2o.dstorage  * dtime )
      cpatch$h2o.precip.cum    = cumsum( cpatch$h2o.precip    * dtime )
      cpatch$h2o.eddy.flux.cum = cumsum( cpatch$h2o.eddy.flux * dtime )
      cpatch$h2o.runoff.cum    = cumsum( cpatch$h2o.runoff    * dtime )
      cpatch$h2o.veg.dyn.cum   = cumsum( cpatch$h2o.veg.dyn   * dtime )
      cpatch$h2o.residual.cum  = cumsum( cpatch$h2o.residual  * dtime )
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #      Re-scale or re-define some variables.                                         #
      #------------------------------------------------------------------------------------#
      cat0("      * Define the relative cumulative sum of the terms.")
      #----- Cumulative terms for all CO2 budget. -----------------------------------------#
      cpatch$co2.dstorage.rel  = 100. * cpatch$co2.dstorage.cum  / co2.scale
      cpatch$co2.nee.rel       = 100. * cpatch$co2.nee.cum       / co2.scale
      cpatch$co2.eddy.flux.rel = 100. * cpatch$co2.eddy.flux.cum / co2.scale
      cpatch$co2.veg.dyn.rel   = 100. * cpatch$co2.veg.dyn.cum   / co2.scale
      cpatch$co2.residual.rel  = 100. * cpatch$co2.residual.cum  / co2.scale
      #----- Cumulative terms for all carbon budget. --------------------------------------#
      cpatch$crb.dstorage.rel  = 100. * cpatch$crb.dstorage.cum  / crb.scale
      cpatch$crb.eddy.flux.rel = 100. * cpatch$crb.eddy.flux.cum / crb.scale
      cpatch$crb.veg.dyn.rel   = 100. * cpatch$crb.veg.dyn.cum   / crb.scale
      cpatch$crb.netseed.rel   = 100. * cpatch$crb.netseed.cum   / crb.scale
      cpatch$crb.residual.rel  = 100. * cpatch$crb.residual.cum  / crb.scale
      #----- Cumulative terms for all enthalpy budget. ------------------------------------#
      cpatch$ent.dstorage.rel  = 100. * cpatch$ent.dstorage.cum  / ent.scale
      cpatch$ent.precip.rel    = 100. * cpatch$ent.precip.cum    / ent.scale
      cpatch$ent.netrad.rel    = 100. * cpatch$ent.netrad.cum    / ent.scale
      cpatch$ent.prss.eff.rel  = 100. * cpatch$ent.prss.eff.cum  / ent.scale
      cpatch$ent.eddy.flux.rel = 100. * cpatch$ent.eddy.flux.cum / ent.scale
      cpatch$ent.runoff.rel    = 100. * cpatch$ent.runoff.cum    / ent.scale
      cpatch$ent.veg.dyn.rel   = 100. * cpatch$ent.veg.dyn.cum   / ent.scale
      cpatch$ent.residual.rel  = 100. * cpatch$ent.residual.cum  / ent.scale
      #----- Cumulative terms for all water budget. ---------------------------------------#
      cpatch$h2o.dstorage.rel  = 100. * cpatch$h2o.dstorage.cum  / h2o.scale
      cpatch$h2o.precip.rel    = 100. * cpatch$h2o.precip.cum    / h2o.scale
      cpatch$h2o.eddy.flux.rel = 100. * cpatch$h2o.eddy.flux.cum / h2o.scale
      cpatch$h2o.runoff.rel    = 100. * cpatch$h2o.runoff.cum    / h2o.scale
      cpatch$h2o.veg.dyn.rel   = 100. * cpatch$h2o.veg.dyn.cum   / h2o.scale
      cpatch$h2o.residual.rel  = 100. * cpatch$h2o.residual.cum  / h2o.scale
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #      Find the summary variables for budgets.  They will be used for plotting and   #
      # to save a summary ascii file.                                                      #
      #------------------------------------------------------------------------------------#
      budget.summ = data.frame(elapsed = e.uniq)
      cat0("      * Summarise budget.")
      for (bb in sequence(nbudget)){

         #----- Retrieve variable information from the list. ------------------------------#
         budget.now   = budget[[bb]]
         vnames       = budget.now$vnam
         averageit    = budget.now$average
         multit       = budget.now$mult
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Append variables.                                                           #
         #---------------------------------------------------------------------------------#
         nlayers   = length(vnames)
         for (l in sequence(nlayers)){
            v.layer = vnames[l]
            #----- Aggregate data into daily groups, to reduce file size. -----------------#
            if (averageit){
               ynow = tapply(X=cpatch[[v.layer]],INDEX=e.aggr,FUN=mean,na.rm=TRUE)
            }else{
               ynow = tapply(X=cpatch[[v.layer]],INDEX=e.aggr,FUN=max ,na.rm=TRUE)
            }#end if (aggrit)
            budget.summ[[v.layer]] = eval(parse(text=paste0("ynow * ",multit)))
            #------------------------------------------------------------------------------#
         }#end for (l in sequence(nlayers))
         #---------------------------------------------------------------------------------#
      }#end for (bb in sequence(nbudget))
      #------------------------------------------------------------------------------------#


      #----- Save summary file. -----------------------------------------------------------#
      summfile = file.path(outroot,paste0(place,"_summ-",budget_midfix,"_",cipa,".csv"))
      dummy    = write.table( x         = budget.summ
                            , file      = summfile
                            , append    = FALSE
                            , quote     = FALSE
                            , sep       = ","
                            , row.names = FALSE
                            , col.names = TRUE
                            )#end write.table
      #------------------------------------------------------------------------------------#







      #------------------------------------------------------------------------------------#
      #   Plot the time series diagrams showing months and years.                          #
      #------------------------------------------------------------------------------------#
      cat0("      * Plot some patch-level figures.")
      for (bb in sequence(nbudget)){

         #----- Retrieve variable information from the list. ------------------------------#
         budget.now   = budget[[bb]]
         vnames       = budget.now$vnam
         description  = budget.now$desc
         lcolours     = budget.now$colour
         llty         = budget.now$lty
         llwd         = budget.now$lwd
         lrange       = budget.now$range
         ltype        = budget.now$type
         prefix       = budget.now$prefix
         theme        = budget.now$theme
         unit         = budget.now$unit
         y.limit      = budget.now$ylim
         legpos       = budget.now$legpos
         plotit       = budget.now$plt
         #---------------------------------------------------------------------------------#

         if (plotit){
            #----- Define the number of layers. -------------------------------------------#
            nlayers   = length(vnames)
            namerange = vnames[lrange]
            langle    = c((sequence(nlayers-1)-0.5)*180 / (nlayers-1),0)
            ldensity  = c(rep(density.terms,nlayers-1),density.resid)
            lalpha    = c(rep(alpha.terms,nlayers-1),alpha.resid)
            #------------------------------------------------------------------------------#



            #----- Set the limits for the y axis. -----------------------------------------#
            if (any(is.na(y.limit))){
               y.range = c(-1.,1.) * max(abs(budget.summ[,namerange]),na.rm=TRUE)
               y.limit = pretty.xylim(u=y.range,fracexp=0.0,is.log=FALSE)
            }#end if (any(is.na(y.limit)))
            y.at     = pretty(y.limit)
            y.labels = sprintf("%g",y.at)
            #------------------------------------------------------------------------------#


            #----- Make plot annotations. -------------------------------------------------#
            letitre = paste0(thispoi$lieu
                            ,"\n Patch: ",ipa,";  Time series: ",theme)
            lex     = desc.unit(desc=e.desc,unit=e.unit       )
            ley     = desc.unit(desc=theme,unit=untab[[unit]])
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Check if the directory exists.  If not, create it.                       #
            #------------------------------------------------------------------------------#
            cat0("        > ",theme," time series.")

            #----- Loop over formats. -----------------------------------------------------#
            for (o in sequence(nout)){
               #----- Open file. ----------------------------------------------------------#
               fichier = paste0(prefix,"-tseries-patch-",cipa,"-",suffix,".",outform[o])
               fichier = file.path(patchdir,fichier)
               dummy   = open.plot( fichier = fichier
                                  , outform = outform[o]
                                  , size    = ey.size
                                  , ptsz    = ptsz
                                  , depth   = depth
                                  )#end open.plot
               #---------------------------------------------------------------------------#


               #----- Split the window into two. ------------------------------------------#
               par(par.user)
               layout( mat     = rbind(2,1)
                     , heights = c(1.-f.leg,f.leg)
                     )#end layout
               #---------------------------------------------------------------------------#



               #---- Plot the legend. -----------------------------------------------------#
               par(mar=c(0.1,4.6,0.1,0.6))
               plot.new()
               plot.window(xlim=c(0,1),ylim=c(0,1))
               legend( x      = "center"
                     , inset  = 0.0
                     , legend = description
                     , col    = lcolours
                     , lty    = if (same.lty){lty.all}else{llty}
                     , lwd    = if (same.lty){lwd.all}else{llwd}
                     , ncol   = 2
                     , cex    = 0.8 * cex.ptsz
                     , xpd    = TRUE
                     , bty    = "n"
                     )#end legend
               #---------------------------------------------------------------------------#


               #---- Open the window and plot the axes and annotations. -------------------#
               par(mar=c(5.1,5.1,3.1,0.6))
               plot.new()
               plot.window(xlim=e.limit,ylim=y.limit,log="")
               axis(side=1,las=1,at=e.at,labels=e.labels)
               axis(side=2,las=1,at=y.at,labels=y.labels)
               if (plotgrid) abline(h=y.at,v=e.at,col=grid.colour,lty="longdash")
               box()
               title(main=letitre,cex.main=0.6)
               title(xlab=lex)
               title(ylab=ley,line=3.5)
               #---------------------------------------------------------------------------#



               #----- Plot the curves. ----------------------------------------------------#
               for (l in sequence(nlayers)){
                  v.layer = vnames[l]
                  lines( x    = budget.summ$elapsed
                       , y    = budget.summ[[v.layer]]
                       , col  = lcolours[l]
                       , lty  = if (same.lty){lty.all}else{llty[l]}
                       , lwd  = if (same.lty){lwd.all}else{llwd[l]}
                       , type = ltype
                       , pch  = 16
                       )#end lines
               }#end for
               #---------------------------------------------------------------------------#


               #----- Close the device. ---------------------------------------------------#
               dummy = close.plot(outform=outform[o])
               #---------------------------------------------------------------------------#
            }#end for outform
            #------------------------------------------------------------------------------#






            #------------------------------------------------------------------------------#
            #     Check if the directory exists.  If not, create it.                       #
            #------------------------------------------------------------------------------#
            cat0("        > ",theme," polygons.")

            #----- Loop over formats. -----------------------------------------------------#
            for (o in sequence(nout)){
               #----- Open file. ----------------------------------------------------------#
               fichier = paste0(prefix,"-polygon-patch-",cipa,"-",suffix,".",outform[o])
               fichier = file.path(patchdir,fichier)
               dummy   = open.plot( fichier = fichier
                                  , outform = outform[o]
                                  , size    = ey.size
                                  , ptsz    = ptsz
                                  , depth   = depth
                                  )#end open.plot
               #---------------------------------------------------------------------------#


               #----- Split the window into two. ------------------------------------------#
               par(par.user)
               layout( mat     = rbind(2,1)
                     , heights = c(1.-f.leg,f.leg)
                     )#end layout
               #---------------------------------------------------------------------------#



               #---- Plot the legend. -----------------------------------------------------#
               par(mar=c(0.1,4.6,0.1,0.6))
               plot.new()
               plot.window(xlim=c(0,1),ylim=c(0,1))
               legend( x       = "center"
                     , inset   = 0.0
                     , legend  = description
                     , fill    = alpha(lcolours,lalpha)
                     , border  = lcolours
                     , angle   = langle
                     , density = ldensity
                     , ncol    = 2
                     , cex     = 0.8 * cex.ptsz
                     , xpd     = TRUE
                     , bty     = "n"
                     )#end legend
               #---------------------------------------------------------------------------#


               #---- Open the window and plot the axes and annotations. -------------------#
               par(mar=c(5.1,5.1,3.1,0.6))
               plot.new()
               plot.window(xlim=e.limit,ylim=y.limit,log="")
               axis(side=1,las=1,at=e.at,labels=e.labels)
               axis(side=2,las=1,at=y.at,labels=y.labels)
               if (plotgrid) abline(h=y.at,v=e.at,col=grid.colour,lty="longdash")
               box()
               title(main=letitre,cex.main=0.6)
               title(xlab=lex)
               title(ylab=ley,line=3.5)
               #---------------------------------------------------------------------------#



               #----- Plot the curves. ----------------------------------------------------#
               for (l in sequence(nlayers)){
                  v.layer = vnames[l]
                  ep.x    = c(budget.summ$elapsed   ,rev(budget.summ$elapsed      ))
                  ep.y    = c(budget.summ[[v.layer]],rev(0.*budget.summ[[v.layer]]))
                  
                  epolygon( x       = ep.x
                          , y       = ep.y
                          , col     = alpha(lcolours[l],lalpha[l])
                          , border  = lcolours[l]
                          , density = ldensity[l]
                          , lwd     = 0.5
                          , angle   = langle  [l]
                          )#end epolygon
                  lines( x    = budget.summ$elapsed
                       , y    = budget.summ[[v.layer]]
                       , col  = lcolours[l]
                       , lty  = if (same.lty){lty.all}else{llty[l]}
                       , lwd  = if (same.lty){lwd.all}else{llwd[l]}
                       , type = ltype
                       , pch  = 16
                       )#end lines
               }#end for
               #---------------------------------------------------------------------------#


               #----- Close the device. ---------------------------------------------------#
               dummy = close.plot(outform=outform[o])
               #---------------------------------------------------------------------------#
            }#end for outform
            #------------------------------------------------------------------------------#
         }#end if plotit
         #---------------------------------------------------------------------------------#
      }#end for nphov
      #------------------------------------------------------------------------------------#
   }#end for (ipa in sequence(npatches))
   #---------------------------------------------------------------------------------------#
}#end for (place in myplaces)
#------------------------------------------------------------------------------------------#
